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

use crate::{asdf_io, fits_io, xisf_io};

/// A container format arcsec can read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImageFormat {
    Fits,
    Xisf,
    Asdf,
}

impl ImageFormat {
    pub fn name(self) -> &'static str {
        match self {
            ImageFormat::Fits => "FITS",
            ImageFormat::Xisf => "XISF",
            ImageFormat::Asdf => "ASDF",
        }
    }
}

/// Identify a file's format from its leading bytes.
///
/// Unrecognised files are rejected here rather than handed to CFITSIO. This used
/// to fall back to FITS on the theory that CFITSIO would report the error for us,
/// but `rsfitsio::fits_open_image` panics instead of setting its status when it
/// cannot open a file, so the fallback turned a mistyped path into a "Null Pointer"
/// panic and exit 101 instead of the documented exit 16. See cruzzil/rsfitsio#136.
///
/// Two things reach CFITSIO deliberately:
///
/// - Compressed FITS. CFITSIO decompresses gzip, bzip2 and Unix `compress`
///   transparently, so those magic numbers are reported as [`ImageFormat::Fits`].
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

/// Load an image's pixels, whatever container they arrived in.
pub fn read_image(path: &Path) -> Result<ImageBuffer, String> {
    match detect_format(path)? {
        ImageFormat::Fits => fits_io::read_fits_image(path),
        ImageFormat::Xisf => xisf_io::read_xisf_image(path),
        ImageFormat::Asdf => asdf_io::read_asdf_image(path),
    }
}

/// An approximate pointing from the file's metadata, in degrees.
pub fn read_ra_dec(path: &Path) -> Option<(f64, f64)> {
    match detect_format(path).ok()? {
        ImageFormat::Fits => fits_io::read_fits_ra_dec(path),
        ImageFormat::Xisf => xisf_io::read_xisf_ra_dec(path),
        ImageFormat::Asdf => asdf_io::read_asdf_ra_dec(path),
    }
}

/// Plate scale in arcsec/pixel from the file's metadata.
pub fn read_pixel_scale(path: &Path) -> Option<f64> {
    match detect_format(path).ok()? {
        ImageFormat::Fits => fits_io::read_fits_pixel_scale(path),
        ImageFormat::Xisf => xisf_io::read_xisf_pixel_scale(path),
        ImageFormat::Asdf => asdf_io::read_asdf_pixel_scale(path),
    }
}

/// Image dimensions without loading the pixels.
pub fn read_dimensions(path: &Path) -> Option<(u32, u32)> {
    match detect_format(path).ok()? {
        ImageFormat::Fits => fits_io::read_fits_dimensions(path),
        ImageFormat::Xisf => xisf_io::read_xisf_dimensions(path),
        ImageFormat::Asdf => asdf_io::read_asdf_dimensions(path),
    }
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
        ImageFormat::Fits => fits_io::update_fits_wcs(path, wcs),
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
    let fl = focallen_mm.filter(|v| *v > 0.0)?;
    let ps = xpixsz_um.filter(|v| *v > 0.0)?;
    let bin = xbinning.filter(|v| *v > 0.0).unwrap_or(1.0);
    Some(ps * bin / fl * 206.265)
}

/// Pointing from the usual keyword pair, preferring the telescope's own report
/// over a reference pixel that may belong to a previous solve.
pub fn ra_dec_from(
    ra: Option<f64>,
    dec: Option<f64>,
    crval1: Option<f64>,
    crval2: Option<f64>,
) -> Option<(f64, f64)> {
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
        // Absent binning is 1, not zero.
        assert!(pixel_scale_from(Some(250.0), Some(3.76), None).is_some());
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
    }
}
