//! Format dispatch for input images.
//!
//! arcsec reads three container formats. They agree on nothing except that a
//! two-dimensional grid of numbers comes out the other end, so each has its own
//! reader and this module picks between them:
//!
//! | Format | Detected by | Reader |
//! |---|---|---|
//! | FITS | `SIMPLE  =` | [`crate::fits_io`], via CFITSIO |
//! | XISF | `XISF0100` | [`crate::xisf_io`], via the `xisf` crate |
//! | ASDF | `#ASDF ` | [`crate::asdf_io`], via the `asdf-rs` crate |
//!
//! Detection reads the file's first bytes rather than trusting the extension: a
//! `.fit`, `.fits`, `.fts` or no extension at all are all common for FITS, and an
//! image that has been renamed should still solve.

use std::io::Read;
use std::path::Path;

use arcsec_core::types::{ImageBuffer, WcsSolution};
use arcsec_core::wcs::TanWcs;

use crate::{asdf_io, fits_io, xisf_io};

/// A container format arcsec can read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImageFormat {
    Fits,
    Xisf,
    Asdf,
}

impl ImageFormat {
    pub const fn name(self) -> &'static str {
        match self {
            Self::Fits => "FITS",
            Self::Xisf => "XISF",
            Self::Asdf => "ASDF",
        }
    }
}

/// Identify a file's format from its leading bytes.
///
/// Unrecognised files are rejected here rather than handed to CFITSIO. Before
/// rsfitsio 0.470.3, `fits_open_image` panicked instead of setting its status when
/// it could not open a file (cruzzil/rsfitsio#136), so a fallback to FITS turned a
/// mistyped path into exit 101. It returns an error now, but naming the formats
/// arcsec reads is still a clearer message than CFITSIO's.
///
/// Two things reach CFITSIO deliberately:
///
/// - Compressed FITS. CFITSIO decompresses gzip, bzip2 and Unix `compress`
///   transparently, so those magic numbers are reported as [`ImageFormat::Fits`].
///   (gzip and `compress` need rsfitsio 0.470.3, which fixed their magic numbers.)
/// - Extended filename syntax (`image.fits[1]`, `image.fits[col>3]`). Those paths
///   do not name a file on disk, so when the open fails and the path carries a
///   `[`, it is passed through for CFITSIO to parse.
pub fn detect_format(path: &Path) -> Result<ImageFormat, String> {
    let mut head = [0u8; 16];
    let n = match std::fs::File::open(path).and_then(|mut f| f.read(&mut head)) {
        Ok(n) => n,
        Err(e) => {
            // CFITSIO extended filename syntax never names a real file.
            if path.to_string_lossy().contains('[') {
                return Ok(ImageFormat::Fits);
            }
            return Err(format!("{}: {e}", path.display()));
        }
    };
    let head = &head[..n];

    if head.starts_with(b"XISF0100") {
        Ok(ImageFormat::Xisf)
    } else if head.starts_with(b"#ASDF ") {
        Ok(ImageFormat::Asdf)
    } else if head.starts_with(b"SIMPLE  =")
        || head.starts_with(b"\x1f\x8b")  // gzip
        || head.starts_with(b"BZh")      // bzip2
        || head.starts_with(b"\x1f\x9d")
    // Unix compress
    {
        Ok(ImageFormat::Fits)
    } else {
        Err(format!(
            "{}: not a FITS, XISF or ASDF image (no recognised magic bytes)",
            path.display()
        ))
    }
}

// ── Limits on what a file may ask for ─────────────────────────────────────────

/// The most pixels (times channels) arcsec will allocate for one image: 2³⁰, a
/// 32768 × 32768 frame, 4 GiB as `f32`.
///
/// A header can declare any dimensions, and a compressed image (tile-compressed
/// FITS, a compressed XISF block, a zlib ASDF block) can legitimately be far
/// larger than its file, so the size on disk cannot bound it. This is several
/// times the largest single-sensor astronomical camera, so it refuses only files
/// that would exhaust memory rather than solve. Fuzzing builds use a much smaller
/// limit so that the rejection path is exercised instead of the allocator.
#[cfg(not(fuzzing))]
pub const MAX_IMAGE_PIXELS: usize = 1 << 30;
#[cfg(fuzzing)]
pub const MAX_IMAGE_PIXELS: usize = 1 << 22;

/// `width × height × channels`, refused if it overflows or exceeds
/// [`MAX_IMAGE_PIXELS`].
pub fn checked_pixel_count(width: usize, height: usize, channels: usize) -> Result<usize, String> {
    width
        .checked_mul(height)
        .and_then(|n| n.checked_mul(channels))
        .filter(|&n| n <= MAX_IMAGE_PIXELS)
        .ok_or_else(|| {
            format!(
                "image dimensions {width}×{height}×{channels} exceed the {MAX_IMAGE_PIXELS}-pixel \
                 limit"
            )
        })
}

/// A zeroed pixel buffer, or an error rather than an abort if memory runs out.
pub fn try_alloc_pixels(npix: usize) -> Result<Vec<f32>, String> {
    let mut data = Vec::new();
    data.try_reserve_exact(npix)
        .map_err(|_| format!("not enough memory for {npix} pixels"))?;
    data.resize(npix, 0.0);
    Ok(data)
}

// ── Reader boundary ───────────────────────────────────────────────────────────

std::thread_local! {
    /// Set while a [`guarded`] reader runs, so the panic hook stays quiet.
    static IN_READER: core::cell::Cell<bool> = const { core::cell::Cell::new(false) };
    /// Where the last quietened panic happened, for the error message.
    static PANIC_SITE: core::cell::RefCell<String> = const { core::cell::RefCell::new(String::new()) };
}

/// Run a third-party reader, turning a panic into an error.
///
/// CFITSIO (rsfitsio), `xisf` and `asdf-rs` parse untrusted bytes, and a malformed
/// file has been found to panic inside them. A panic would otherwise end the
/// process with Rust's exit code 101 and no `.ini`, where the ASTAP contract says
/// an unreadable image is exit 16, and a header keyword that panics would sink an
/// image that is otherwise fine. The panic's message and location go into the
/// error, which is what a bug report needs, rather than being printed by the hook
/// as an apparent crash.
fn guarded<T>(path: &Path, read: impl FnOnce() -> Result<T, String>) -> Result<T, String> {
    // Under the fuzzer a panic in a reader is a finding, to be fixed where it is.
    #[cfg(fuzzing)]
    if true {
        let _ = path;
        return read();
    }
    static HOOK: std::sync::Once = std::sync::Once::new();
    HOOK.call_once(|| {
        let previous = std::panic::take_hook();
        std::panic::set_hook(Box::new(move |info| {
            if IN_READER.with(core::cell::Cell::get) {
                let site = info
                    .location()
                    .map(|l| format!(" at {}:{}", l.file(), l.line()))
                    .unwrap_or_default();
                PANIC_SITE.with(|s| *s.borrow_mut() = site);
            } else {
                previous(info);
            }
        }));
    });
    let outer = IN_READER.with(|g| g.replace(true));
    let result = std::panic::catch_unwind(core::panic::AssertUnwindSafe(read));
    IN_READER.with(|g| g.set(outer));
    result.unwrap_or_else(|panic| {
        Err(format!(
            "{}: the reader failed on this file, which is probably corrupt ({}{})",
            path.display(),
            panic_message(panic.as_ref()),
            PANIC_SITE.with(core::cell::RefCell::take)
        ))
    })
}

/// The text of a panic payload, for an error message.
pub fn panic_message(payload: &(dyn core::any::Any + Send)) -> &str {
    payload
        .downcast_ref::<&str>()
        .copied()
        .or_else(|| payload.downcast_ref::<String>().map(String::as_str))
        .unwrap_or("panic")
}

/// Load an image's pixels, whatever container they arrived in.
pub fn read_image(path: &Path) -> Result<ImageBuffer, String> {
    let format = detect_format(path)?;
    let img = guarded(path, || match format {
        ImageFormat::Fits => fits_io::read_fits_image(path),
        ImageFormat::Xisf => xisf_io::read_xisf_image(path),
        ImageFormat::Asdf => asdf_io::read_asdf_image(path),
    })?;
    // Every reader promises this; the solver indexes `data[y * width + x]`.
    if img.data.len() != img.width * img.height {
        return Err(format!(
            "{}: the reader returned {} samples for a {}×{} image",
            path.display(),
            img.data.len(),
            img.width,
            img.height
        ));
    }
    Ok(img)
}

/// An approximate pointing from the file's metadata, in degrees.
pub fn read_ra_dec(path: &Path) -> Option<(f64, f64)> {
    let format = detect_format(path).ok()?;
    guarded(path, || {
        Ok(match format {
            ImageFormat::Fits => fits_io::read_fits_ra_dec(path),
            ImageFormat::Xisf => xisf_io::read_xisf_ra_dec(path),
            ImageFormat::Asdf => asdf_io::read_asdf_ra_dec(path),
        })
    })
    .inspect_err(|e| log::warn!("{e}"))
    .ok()
    .flatten()
}

/// Plate scale in arcsec/pixel from the file's metadata.
pub fn read_pixel_scale(path: &Path) -> Option<f64> {
    let format = detect_format(path).ok()?;
    guarded(path, || {
        Ok(match format {
            ImageFormat::Fits => fits_io::read_fits_pixel_scale(path),
            ImageFormat::Xisf => xisf_io::read_xisf_pixel_scale(path),
            ImageFormat::Asdf => asdf_io::read_asdf_pixel_scale(path),
        })
    })
    .inspect_err(|e| log::warn!("{e}"))
    .ok()
    .flatten()
}

/// Image dimensions without loading the pixels.
pub fn read_dimensions(path: &Path) -> Option<(u32, u32)> {
    let format = detect_format(path).ok()?;
    guarded(path, || {
        Ok(match format {
            ImageFormat::Fits => fits_io::read_fits_dimensions(path),
            ImageFormat::Xisf => xisf_io::read_xisf_dimensions(path),
            ImageFormat::Asdf => asdf_io::read_asdf_dimensions(path),
        })
    })
    .inspect_err(|e| log::warn!("{e}"))
    .ok()
    .flatten()
}

/// A TAN WCS already in the file's header (from an earlier solve), if any.
///
/// Used by `--extract` for the RA and Dec columns. Like ASTAP, only a CD matrix
/// counts: a header with CDELT but no CD gives none, and SIP keywords are ignored.
pub fn read_header_wcs(path: &Path) -> Option<TanWcs> {
    let format = detect_format(path).ok()?;
    guarded(path, || {
        Ok(match format {
            ImageFormat::Fits => fits_io::read_fits_header_wcs(path),
            ImageFormat::Xisf => xisf_io::read_xisf_header_wcs(path),
            ImageFormat::Asdf => None,
        })
    })
    .inspect_err(|e| log::warn!("{e}"))
    .ok()
    .flatten()
}

/// Number of colour channels the image had before it was reduced to one.
pub fn read_channels(path: &Path) -> usize {
    let Ok(format) = detect_format(path) else {
        return 1;
    };
    guarded(path, || {
        Ok(match format {
            ImageFormat::Fits => fits_io::read_fits_channels(path),
            ImageFormat::Xisf => xisf_io::read_xisf_channels(path),
            ImageFormat::Asdf => 1,
        })
    })
    .inspect_err(|e| log::warn!("{e}"))
    .unwrap_or(1)
}

/// Write the solution back into the image file (`--update`).
///
/// Only FITS supports this. XISF and ASDF would each need the whole file
/// rewritten — neither crate edits headers in place — and silently rewriting a
/// user's master light to add four keywords is not a reasonable default. The
/// `.wcs` and `.ini` sidecars are written for every format regardless, and
/// `main` reports this as a warning rather than a failure.
pub fn update_wcs(path: &Path, wcs: &WcsSolution) -> Result<(), String> {
    match detect_format(path)? {
        ImageFormat::Fits => guarded(path, || fits_io::update_fits_wcs(path, wcs)),
        other => Err(format!(
            "--update writes the solution into the image header, which is only \
             supported for FITS (this file is {}). The .wcs and .ini files were \
             still written.",
            other.name()
        )),
    }
}

// ── Shared metadata interpretation ──────────────────────────────────────────────

/// Plate scale from the optics keywords, shared by every format that carries them.
///
/// `plate_scale [arcsec/px] = pixel_size[µm] × binning / focal_length[mm] × 206.265`
///
/// Kept here rather than duplicated per reader so the formula and the
/// what-counts-as-missing rules stay in one place: XISF files written by
/// PixInsight carry the same FOCALLEN/XPIXSZ/XBINNING keywords a FITS file does.
pub fn pixel_scale_from(
    focallen_mm: Option<f64>,
    xpixsz_um: Option<f64>,
    xbinning: Option<f64>,
) -> Option<f64> {
    let usable = |v: &f64| v.is_finite() && *v > 0.0;
    let fl = focallen_mm.filter(usable)?;
    let ps = xpixsz_um.filter(usable)?;
    let bin = xbinning.filter(usable).unwrap_or(1.0);
    // The field size and the search step are derived from this, so a scale no
    // instrument has (FOCALLEN = 1E300, say, from a corrupt header) is no guess at
    // all: it would make the search step vanishingly small.
    Some(ps * bin / fl * 206.265).filter(|s| PLAUSIBLE_PIXEL_SCALE.contains(s))
}

/// Pixel scales, arcseconds per pixel, that some real instrument could have: from
/// a large telescope's adaptive optics to an all-sky fisheye, with a wide margin.
const PLAUSIBLE_PIXEL_SCALE: core::ops::RangeInclusive<f64> = 1e-3..=1e4;

/// A TAN WCS from header keywords, shared by the formats that carry them.
///
/// `key` looks a numeric keyword up by name. CRVAL1/2, CRPIX1/2 and a non-zero
/// `CD1_1` are required; missing off-diagonal terms are zero.
pub fn tan_wcs_from(mut key: impl FnMut(&str) -> Option<f64>) -> Option<TanWcs> {
    let cd1_1 = key("CD1_1").filter(|v| *v != 0.0)?;
    let wcs = TanWcs {
        ra0: key("CRVAL1")?.to_radians(),
        dec0: key("CRVAL2")?.to_radians(),
        crpix1: key("CRPIX1")?,
        crpix2: key("CRPIX2")?,
        cd: [
            [cd1_1, key("CD1_2").unwrap_or(0.0)],
            [key("CD2_1").unwrap_or(0.0), key("CD2_2").unwrap_or(0.0)],
        ],
        sip: None,
    };
    let finite = [wcs.ra0, wcs.dec0, wcs.crpix1, wcs.crpix2]
        .into_iter()
        .chain(wcs.cd.into_iter().flatten())
        .all(f64::is_finite);
    finite.then_some(wcs)
}

/// Pointing from the usual keyword pair, preferring the telescope's own report
/// over a reference pixel that may belong to a previous solve.
pub fn ra_dec_from(
    ra: Option<f64>,
    dec: Option<f64>,
    crval1: Option<f64>,
    crval2: Option<f64>,
) -> Option<(f64, f64)> {
    let finite = |v: &f64| v.is_finite();
    let (ra, dec) = (ra.filter(finite), dec.filter(finite));
    let (crval1, crval2) = (crval1.filter(finite), crval2.filter(finite));
    Some((ra.or(crval1)?, dec.or(crval2)?))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(name: &str, bytes: &[u8]) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(name);
        std::fs::write(&p, bytes).unwrap();
        p
    }

    #[test]
    fn detects_by_content_not_extension() {
        let f = write("arcsec_fmt_a.dat", b"SIMPLE  =                    T");
        assert_eq!(detect_format(&f), Ok(ImageFormat::Fits));
        let x = write("arcsec_fmt_b.fits", b"XISF0100\x00\x00\x00\x00");
        assert_eq!(
            detect_format(&x),
            Ok(ImageFormat::Xisf),
            "extension must not win"
        );
        let a = write("arcsec_fmt_c.fits", b"#ASDF 1.0.0\n#ASDF_STANDARD 1.6.0\n");
        assert_eq!(detect_format(&a), Ok(ImageFormat::Asdf));
        for p in [f, x, a] {
            std::fs::remove_file(p).ok();
        }
    }

    #[test]
    fn unknown_and_missing_files_are_rejected() {
        // Neither may reach CFITSIO: rsfitsio panics rather than reporting a
        // status when it cannot open a file, which costs us the exit-16 contract.
        let u = write("arcsec_fmt_d.bin", b"not any known magic");
        assert!(
            detect_format(&u).is_err(),
            "junk must not be taken for FITS"
        );
        std::fs::remove_file(&u).ok();
        assert!(
            detect_format(Path::new("/nonexistent/arcsec/none.fits")).is_err(),
            "a missing file must not be taken for FITS"
        );
    }

    #[test]
    fn compressed_fits_and_extended_syntax_still_reach_cfitsio() {
        // CFITSIO decompresses these itself, so the magic bytes are its, not ours.
        let g = write("arcsec_fmt_e.fits", b"\x1f\x8b\x08\x00rest-is-deflate");
        assert_eq!(detect_format(&g), Ok(ImageFormat::Fits), "gzip");
        std::fs::remove_file(&g).ok();
        let b = write("arcsec_fmt_f.fits", b"BZh91AY&SY");
        assert_eq!(detect_format(&b), Ok(ImageFormat::Fits), "bzip2");
        std::fs::remove_file(&b).ok();
        // Extended filename syntax names no file on disk; CFITSIO parses it.
        assert_eq!(
            detect_format(Path::new("/nonexistent/arcsec/none.fits[1]")),
            Ok(ImageFormat::Fits)
        );
    }

    #[test]
    fn pixel_scale_formula() {
        // 3.76 µm pixels on a 250 mm lens: 3.76 / 250 * 206.265 ≈ 3.10"/px
        let ps = pixel_scale_from(Some(250.0), Some(3.76), Some(1.0)).unwrap();
        assert!((ps - 3.1022).abs() < 1e-3, "got {ps}");
        // Binning multiplies it.
        let ps2 = pixel_scale_from(Some(250.0), Some(3.76), Some(2.0)).unwrap();
        assert!((ps2 - 2.0 * ps).abs() < 1e-9);
        // Missing or nonsense optics yield nothing rather than a wrong guess.
        assert!(pixel_scale_from(None, Some(3.76), None).is_none());
        assert!(pixel_scale_from(Some(250.0), None, None).is_none());
        assert!(pixel_scale_from(Some(0.0), Some(3.76), None).is_none());
        assert!(pixel_scale_from(Some(250.0), Some(-1.0), None).is_none());
        // Nor do values no instrument has, which a corrupt header can hold: they
        // would make the search step vanish (or the field cover the sky).
        assert!(pixel_scale_from(Some(1e300), Some(3.76), None).is_none());
        assert!(pixel_scale_from(Some(1e-300), Some(3.76), None).is_none());
        assert!(pixel_scale_from(Some(f64::INFINITY), Some(3.76), None).is_none());
        assert!(pixel_scale_from(Some(250.0), Some(f64::NAN), None).is_none());
        // An unusable binning is ignored, as an absent one is.
        assert_eq!(
            pixel_scale_from(Some(250.0), Some(3.76), Some(f64::INFINITY)),
            pixel_scale_from(Some(250.0), Some(3.76), None)
        );
        // Absent binning is 1, not zero.
        assert!(pixel_scale_from(Some(250.0), Some(3.76), None).is_some());
    }

    #[test]
    fn a_header_wcs_needs_a_cd_matrix() {
        let full = |k: &str| match k {
            "CRVAL1" => Some(150.0),
            "CRVAL2" => Some(40.0),
            "CRPIX1" => Some(100.5),
            "CRPIX2" => Some(80.5),
            "CD1_1" => Some(-2e-4),
            "CD2_2" => Some(2e-4),
            _ => None,
        };
        let w = tan_wcs_from(full).unwrap();
        assert_eq!(w.cd, [[-2e-4, 0.0], [0.0, 2e-4]]);
        let (ra, dec) = w.pixel_to_sky(100.5, 80.5);
        assert!((ra.to_degrees() - 150.0).abs() < 1e-9 && (dec.to_degrees() - 40.0).abs() < 1e-9);
        // CDELT alone, or a zero CD1_1, is not a solution.
        assert!(tan_wcs_from(|k| if k == "CD1_1" { None } else { full(k) }).is_none());
        assert!(tan_wcs_from(|k| if k == "CD1_1" { Some(0.0) } else { full(k) }).is_none());
        assert!(tan_wcs_from(|k| if k == "CRPIX2" { None } else { full(k) }).is_none());
    }

    #[test]
    fn ra_dec_prefers_the_telescope_over_crval() {
        assert_eq!(
            ra_dec_from(Some(10.0), Some(20.0), Some(99.0), Some(99.0)),
            Some((10.0, 20.0))
        );
        assert_eq!(
            ra_dec_from(None, None, Some(30.0), Some(40.0)),
            Some((30.0, 40.0))
        );
        // A half-present pair is not a pointing.
        assert_eq!(ra_dec_from(Some(10.0), None, None, None), None);
        // Nor is one that is not a number; CRVAL can stand in for it.
        assert_eq!(ra_dec_from(Some(f64::NAN), Some(20.0), None, None), None);
        assert_eq!(
            ra_dec_from(Some(f64::INFINITY), Some(20.0), Some(30.0), Some(40.0)),
            Some((30.0, 20.0))
        );
    }

    #[test]
    fn a_panicking_reader_is_an_error_naming_where() {
        let p = Path::new("x.fits");
        let err = guarded::<()>(p, || panic!("boom {}", 42)).unwrap_err();
        assert!(
            err.starts_with("x.fits: ") && err.contains("boom 42"),
            "{err}"
        );
        assert!(err.contains("image_io.rs"), "location missing: {err}");
        // The hook is restored for panics outside a reader: nested use too.
        assert_eq!(guarded(p, || guarded(p, || Ok(7))), Ok(7));
        assert!(!IN_READER.with(core::cell::Cell::get));
    }

    #[test]
    fn pixel_counts_are_checked_before_allocating() {
        assert_eq!(checked_pixel_count(4, 3, 2), Ok(24));
        assert!(checked_pixel_count(usize::MAX, 2, 1).is_err());
        assert!(checked_pixel_count(1 << 20, 1 << 20, 1).is_err());
        assert!(checked_pixel_count(MAX_IMAGE_PIXELS, 1, 1).is_ok());
        assert!(checked_pixel_count(MAX_IMAGE_PIXELS, 1, 2).is_err());
        assert_eq!(try_alloc_pixels(5).unwrap(), vec![0.0; 5]);
    }
}
