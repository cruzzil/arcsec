use alloc::ffi::CString;
use std::path::Path;

use libc::{c_char, c_int, c_long};
use rsfitsio::aliases::rust_api::*;
use rsfitsio::fitsio::{LONGLONG, READONLY, READWRITE, fitsfile};

use arcsec_core::types::{ImageBuffer, WcsSolution};

use crate::image_io;

/// Reinterpret a byte slice as a `c_char` slice for the rsfitsio wrappers.
///
/// `libc::c_char` is `i8` on x86_64 and `u8` on aarch64, so `b"KEY\0"` literals
/// cannot be passed directly on every platform. Same size and alignment either way.
#[inline]
const fn cc(b: &[u8]) -> &[c_char] {
    // Safety: c_char is i8 or u8; identical layout, and we only read.
    unsafe { core::slice::from_raw_parts(b.as_ptr().cast::<c_char>(), b.len()) }
}

/// An open CFITSIO file, closed when dropped.
///
/// Every reader used to open, read and then close by hand on each exit path, and
/// a `?` in between leaked the handle. Closing in `Drop` makes that impossible.
struct FitsFile(Option<Box<fitsfile>>);

impl FitsFile {
    /// Open the primary image of `path` with `mode` (`READONLY`/`READWRITE`).
    ///
    /// Callers reach this only through `image_io`, which has already checked that
    /// the file exists and is FITS: `rsfitsio` panics rather than reporting a
    /// status for a file it cannot open (cruzzil/rsfitsio#136).
    fn open(path: &Path, mode: c_int) -> Result<Self, String> {
        let path_str = path.to_str().ok_or("non-UTF-8 path")?;
        let cpath = CString::new(path_str).map_err(|e| e.to_string())?;
        let mut fptr: Option<Box<fitsfile>> = None;
        let mut status: c_int = 0;
        fits_open_image(&mut fptr, cc(cpath.to_bytes_with_nul()), mode, &mut status);
        let file = Self(fptr);
        if status != 0 {
            return Err(format!("fits_open_image failed: status {status}"));
        }
        if file.0.is_none() {
            return Err("fits_open_image returned no file".to_string());
        }
        Ok(file)
    }

    fn fp(&mut self) -> &mut fitsfile {
        self.0
            .as_deref_mut()
            .unwrap_or_else(|| unreachable!("FitsFile is only built around an open file"))
    }

    /// A header keyword as a double, or `None` if absent or not numeric.
    fn key_f64(&mut self, key: &[u8]) -> Option<f64> {
        let mut val = 0.0f64;
        let mut st: c_int = 0;
        fits_read_key_dbl(self.fp(), cc(key), &mut val, None, &mut st);
        (st == 0).then_some(val)
    }

    /// A header keyword as an integer, or `None` if absent or not numeric.
    fn key_i64(&mut self, key: &[u8]) -> Option<i64> {
        // fits_read_key_lng writes a c_long, which is i32 on Windows and i64 on
        // unix, so the destination must be c_long and widen on the way out.
        let mut val: c_long = 0;
        let mut st: c_int = 0;
        fits_read_key_lng(self.fp(), cc(key), &mut val, None, &mut st);
        // A no-op where c_long is already i64, which clippy flags on those targets.
        #[allow(clippy::useless_conversion)]
        let val = i64::from(val);
        (st == 0).then_some(val)
    }

    /// Close now and report CFITSIO's status, which for a file opened for writing
    /// is where a failed flush of the updated header shows up.
    fn close(mut self) -> c_int {
        let mut status: c_int = 0;
        if let Some(b) = self.0.take() {
            fits_close_file(b, &mut status);
        }
        status
    }
}

impl Drop for FitsFile {
    fn drop(&mut self) {
        if let Some(b) = self.0.take() {
            let mut status: c_int = 0;
            fits_close_file(b, &mut status);
        }
    }
}

/// Read RA and DEC (in degrees) from a FITS header.
/// Tries `RA`/`DEC` first (telescope pointing, NINA style), then `CRVAL1`/`CRVAL2`.
pub fn read_fits_ra_dec(path: &Path) -> Option<(f64, f64)> {
    let mut f = FitsFile::open(path, READONLY).ok()?;
    image_io::ra_dec_from(
        f.key_f64(b"RA\0"),
        f.key_f64(b"DEC\0"),
        f.key_f64(b"CRVAL1\0"),
        f.key_f64(b"CRVAL2\0"),
    )
}

/// Read NAXIS1/NAXIS2 from a FITS header without loading the pixel data.
///
/// Used by `arcsec catalog recommend --like`, which only needs the field size.
pub fn read_fits_dimensions(path: &Path) -> Option<(u32, u32)> {
    let mut f = FitsFile::open(path, READONLY).ok()?;
    let w = u32::try_from(f.key_i64(b"NAXIS1\0")?).ok()?;
    let h = u32::try_from(f.key_i64(b"NAXIS2\0")?).ok()?;
    Some((w, h))
}

/// Read pixel scale (arcsec/pixel, accounting for XBINNING) from FITS header.
/// Returns `None` if FOCALLEN or XPIXSZ are missing or non-positive.
pub fn read_fits_pixel_scale(path: &Path) -> Option<f64> {
    let mut f = FitsFile::open(path, READONLY).ok()?;
    image_io::pixel_scale_from(
        f.key_f64(b"FOCALLEN\0"),
        f.key_f64(b"XPIXSZ\0"),
        f.key_f64(b"XBINNING\0"),
    )
}

/// Open a FITS image and return its pixel data as an `ImageBuffer`.
///
/// A data cube (NAXIS = 3, e.g. an RGB image) is read as its first plane.
pub fn read_fits_image(path: &Path) -> Result<ImageBuffer, String> {
    let mut f = FitsFile::open(path, READONLY)?;
    let mut status: c_int = 0;

    // Read image header: SIMPLE, BITPIX, NAXIS, NAXIS1/2, PCOUNT, GCOUNT, EXTEND
    let mut simple: c_int = 0;
    let mut bitpix: c_int = 0;
    let mut naxis: c_int = 0;
    let mut naxes: [c_long; 9] = [0; 9];
    let mut pcount: c_long = 0;
    let mut gcount: c_long = 0;
    let mut extend: c_int = 0;
    fits_read_imghdr(
        f.fp(),
        9,
        &mut simple,
        &mut bitpix,
        Some(&mut naxis),
        Some(&mut naxes[..]),
        &mut pcount,
        &mut gcount,
        &mut extend,
        &mut status,
    );
    if status != 0 {
        return Err(format!("fits_read_imghdr failed: status {status}"));
    }
    if naxis < 2 {
        return Err(format!("FITS image has only {naxis} axes (need ≥ 2)"));
    }

    let dim = |n: c_long| usize::try_from(n).ok().filter(|&d| d > 0);
    let (Some(width), Some(height)) = (dim(naxes[0]), dim(naxes[1])) else {
        return Err(format!(
            "invalid image dimensions {}×{}",
            naxes[0], naxes[1]
        ));
    };
    let npix = width
        .checked_mul(height)
        .ok_or_else(|| format!("FITS image dimensions {width}×{height} overflow"))?;
    let mut data: Vec<f32> = vec![0.0f32; npix];

    // fits_read_img_flt = ffgpve_safe: reads directly as f32, any BITPIX
    fits_read_img_flt(
        f.fp(),
        1, // group 1
        1, // firstelem 1-based
        npix as LONGLONG,
        0.0f32, // null pixel replacement
        &mut data,
        None, // anynul
        &mut status,
    );
    if status != 0 {
        return Err(format!("fits_read_img_flt failed: status {status}"));
    }

    Ok(ImageBuffer {
        data,
        width,
        height,
    })
}

// ── FITS card helpers ─────────────────────────────────────────────────────────

/// Format a float as a right-justified 20-char FITS scientific notation value.
/// e.g.  1.608758172485E+002  or -5.952377960677E+001
///
/// A non-finite value has no exponent to reformat and is written as Rust prints it
/// (`NaN`, `inf`), right-justified, rather than panicking.
fn fits_f64(v: f64) -> String {
    let raw = format!("{v:.12E}");
    let Some((mantissa, exp_str)) = raw.split_once('E') else {
        return format!("{raw:>20}");
    };
    let exp: i32 = exp_str.parse().unwrap_or(0);
    // {:>20} right-justifies; {:+04} gives e.g. "+002" or "-004"
    format!("{:>20}", format!("{mantissa}E{exp:+04}"))
}

/// Build a float FITS card: `KEYWORD = <20-char value> / <comment>` (80 chars).
fn dbl_card(keyword: &str, value: f64, comment: &str) -> String {
    let card = format!("{keyword:<8}= {} / {comment:<47}", fits_f64(value));
    format!("{:<80}", &card[..80.min(card.len())])
}

/// Build a string FITS card: `KEYWORD = '<value>'<pad>  / <comment>` (80 chars).
fn str_card(keyword: &str, value: &str, comment: &str) -> String {
    let quoted = format!("'{value}'");
    let card = format!("{keyword:<8}= {quoted:<20} / {comment:<47}");
    format!("{:<80}", &card[..80.min(card.len())])
}

/// Build a logical FITS card: `KEYWORD =                    T / <comment>` (80 chars).
fn log_card(keyword: &str, value: bool, comment: &str) -> String {
    let v = if value { "T" } else { "F" };
    let card = format!("{keyword:<8}= {v:>20} / {comment:<47}");
    format!("{:<80}", &card[..80.min(card.len())])
}

// ── Public output writers ─────────────────────────────────────────────────────

/// Write WCS solution to a `.wcs` ASCII FITS-header file (ASTAP-compatible format).
pub fn write_wcs_file(path: &Path, wcs: &WcsSolution) -> std::io::Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;

    let ra_deg = wcs.ra0.to_degrees();
    let dec_deg = wcs.dec0.to_degrees();

    let cards = [
        str_card(
            "CTYPE1",
            "RA---TAN",
            "first parameter RA,    projection TANgential",
        ),
        str_card(
            "CTYPE2",
            "DEC--TAN",
            "second parameter DEC,  projection TANgential",
        ),
        str_card("CUNIT1", "deg     ", "Unit of coordinates"),
        dbl_card("CRPIX1", wcs.crpix1, "X of reference pixel"),
        dbl_card("CRPIX2", wcs.crpix2, "Y of reference pixel"),
        dbl_card("CRVAL1", ra_deg, "RA of reference pixel (deg)"),
        dbl_card("CRVAL2", dec_deg, "DEC of reference pixel (deg)"),
        // Use abs CDELT (unsigned pixel scale, matching ASTAP convention)
        dbl_card("CDELT1", wcs.cdelt1.abs(), "X pixel size (deg)"),
        dbl_card("CDELT2", wcs.cdelt2.abs(), "Y pixel size (deg)"),
        dbl_card("CROTA1", wcs.crota2, "Image twist of X axis        (deg)"),
        dbl_card("CROTA2", wcs.crota2, "Image twist of Y axis        (deg)"),
        dbl_card(
            "CD1_1",
            wcs.cd1_1,
            "CD matrix to convert (x,y) to (Ra, Dec)",
        ),
        dbl_card(
            "CD1_2",
            wcs.cd1_2,
            "CD matrix to convert (x,y) to (Ra, Dec)",
        ),
        dbl_card(
            "CD2_1",
            wcs.cd2_1,
            "CD matrix to convert (x,y) to (Ra, Dec)",
        ),
        dbl_card(
            "CD2_2",
            wcs.cd2_2,
            "CD matrix to convert (x,y) to (Ra, Dec)",
        ),
        log_card("PLTSOLVD", true, SOLVED_BY),
    ];

    for card in &cards {
        writeln!(f, "{card}")?;
    }
    writeln!(f, "{:<80}", "END")?;

    Ok(())
}

/// Write an `.ini` summary file with key=value pairs (ASTAP-compatible layout).
///
/// The keys and their order follow ASTAP's `.ini`. N.I.N.A. reads `CRPIX1/2` and the
/// `CD` matrix from it unconditionally (and nothing from the `.wcs`), so a solve
/// without them is a failed solve there. `NSTARS`, `NQUADS` and `RMS` are arcsec's
/// own additions; readers look keys up by name.
pub fn write_ini_file(
    path: &Path,
    wcs: &WcsSolution,
    nstars: usize,
    cmdline: &str,
) -> std::io::Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "PLTSOLVD=T")?;
    writeln!(f, "CRPIX1={:.13E}", wcs.crpix1)?;
    writeln!(f, "CRPIX2={:.13E}", wcs.crpix2)?;
    writeln!(f, "CRVAL1={:.9}", wcs.ra0.to_degrees())?;
    writeln!(f, "CRVAL2={:.9}", wcs.dec0.to_degrees())?;
    writeln!(f, "CDELT1={:.13E}", wcs.cdelt1)?;
    writeln!(f, "CDELT2={:.13E}", wcs.cdelt2)?;
    // One rotation, as in the .wcs file, which also writes it as both.
    writeln!(f, "CROTA1={:.4}", wcs.crota2)?;
    writeln!(f, "CROTA2={:.4}", wcs.crota2)?;
    writeln!(f, "CD1_1={:.13E}", wcs.cd1_1)?;
    writeln!(f, "CD1_2={:.13E}", wcs.cd1_2)?;
    writeln!(f, "CD2_1={:.13E}", wcs.cd2_1)?;
    writeln!(f, "CD2_2={:.13E}", wcs.cd2_2)?;
    writeln!(f, "CMDLINE={cmdline}")?;
    writeln!(f, "NSTARS={nstars}")?;
    writeln!(f, "NQUADS={}", wcs.stars_matched)?;
    writeln!(f, "RMS={:.4}", wcs.residual_rms)?;
    Ok(())
}

/// Write the `.ini` for a failed solve: `PLTSOLVD=F` and the command line, as
/// ASTAP does.
pub fn write_unsolved_ini_file(path: &Path, cmdline: &str) -> std::io::Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "PLTSOLVD=F")?;
    writeln!(f, "CMDLINE={cmdline}")?;
    Ok(())
}

/// `PLTSOLVD` comment, naming the version that wrote the solution.
const SOLVED_BY: &str = concat!("Astrometric solved by arcsec ", env!("CARGO_PKG_VERSION"));

/// Write WCS keywords back into the FITS file header in-place (`--update` flag).
pub fn update_fits_wcs(path: &Path, wcs: &WcsSolution) -> Result<(), String> {
    // CFITSIO wants NUL-terminated strings; every literal here is ASCII without
    // an interior NUL, so the conversion cannot fail.
    fn c(s: &str) -> CString {
        CString::new(s).unwrap_or_default()
    }

    let mut f = FitsFile::open(path, READWRITE).map_err(|e| format!("{e} (opening for update)"))?;
    let fp = f.fp();
    let mut status: c_int = 0;
    let decim: c_int = 12;

    // CFITSIO routines do nothing once `status` is non-zero, so the first failure
    // is the one reported, and the file is still closed on drop.
    for (key, val, com) in [
        (
            "CTYPE1",
            "RA---TAN",
            "first parameter RA,    projection TANgential",
        ),
        (
            "CTYPE2",
            "DEC--TAN",
            "second parameter DEC,  projection TANgential",
        ),
        ("CUNIT1", "deg", "Unit of coordinates"),
        ("CUNIT2", "deg", "Unit of coordinates"),
    ] {
        let (key, val, com) = (c(key), c(val), c(com));
        fits_update_key_str(
            fp,
            cc(key.as_bytes_with_nul()),
            cc(val.as_bytes_with_nul()),
            Some(cc(com.as_bytes_with_nul())),
            &mut status,
        );
    }
    for (key, val, com) in [
        ("CRPIX1", wcs.crpix1, "X of reference pixel"),
        ("CRPIX2", wcs.crpix2, "Y of reference pixel"),
        (
            "CRVAL1",
            wcs.ra0.to_degrees(),
            "RA of reference pixel (deg)",
        ),
        (
            "CRVAL2",
            wcs.dec0.to_degrees(),
            "DEC of reference pixel (deg)",
        ),
        ("CDELT1", wcs.cdelt1.abs(), "X pixel size (deg)"),
        ("CDELT2", wcs.cdelt2.abs(), "Y pixel size (deg)"),
        ("CROTA1", wcs.crota2, "Image twist of X axis (deg)"),
        ("CROTA2", wcs.crota2, "Image twist of Y axis (deg)"),
        ("CD1_1", wcs.cd1_1, "CD matrix element"),
        ("CD1_2", wcs.cd1_2, "CD matrix element"),
        ("CD2_1", wcs.cd2_1, "CD matrix element"),
        ("CD2_2", wcs.cd2_2, "CD matrix element"),
    ] {
        let (key, com) = (c(key), c(com));
        fits_update_key_dbl(
            fp,
            cc(key.as_bytes_with_nul()),
            val,
            decim,
            Some(cc(com.as_bytes_with_nul())),
            &mut status,
        );
    }
    let (key, com) = (c("PLTSOLVD"), c(SOLVED_BY));
    fits_update_key_log(
        fp,
        cc(key.as_bytes_with_nul()),
        1,
        Some(cc(com.as_bytes_with_nul())),
        &mut status,
    );
    let close_status = f.close();

    if status == 0 && close_status != 0 {
        status = close_status;
    }
    if status != 0 {
        return Err(format!("update_fits_wcs failed: status {status}"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use arcsec_core::types::PlateConstants;

    fn solution() -> WcsSolution {
        WcsSolution {
            ra0: 15f64.to_radians(),
            dec0: 0.0,
            crpix1: 1450.5,
            crpix2: 1450.5,
            cd1_1: -3.448e-4,
            cd1_2: -1.8e-8,
            cd2_1: -1.1e-8,
            cd2_2: 3.448e-4,
            cdelt1: -3.448e-4,
            cdelt2: 3.448e-4,
            crota2: -0.0019,
            residual_rms: 0.7,
            stars_matched: 144,
            plate: PlateConstants {
                a: 1.0,
                b: 0.0,
                c: 0.0,
                d: 0.0,
                e: 1.0,
                f: 0.0,
            },
            mag_limit: 20.0,
            search_dist_deg: 0.0,
            step_distances: Vec::new(),
            raw_matches: 144,
        }
    }

    /// Parse an .ini the way N.I.N.A.'s ASTAP solver does: split each non-blank line
    /// on the first `=`, into a dictionary that rejects duplicate keys.
    fn read_like_nina(path: &Path) -> std::collections::HashMap<String, String> {
        let mut dict = std::collections::HashMap::new();
        for line in std::fs::read_to_string(path).unwrap().lines() {
            if line.trim().is_empty() {
                continue;
            }
            let (k, v) = line.split_once('=').expect("every line has a '='");
            assert!(
                dict.insert(k.to_string(), v.to_string()).is_none(),
                "duplicate {k}"
            );
        }
        dict
    }

    #[test]
    fn solved_ini_has_every_key_nina_reads() {
        let path = std::env::temp_dir().join(format!("arcsec_ini_{}.ini", std::process::id()));
        write_ini_file(&path, &solution(), 500, "arcsec -f x.fits -fov 1").unwrap();
        let dict = read_like_nina(&path);
        std::fs::remove_file(&path).ok();

        assert_eq!(dict["PLTSOLVD"], "T");
        // N.I.N.A. indexes these directly and parses them as invariant-culture doubles.
        for key in [
            "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "CD1_1", "CD1_2", "CD2_1", "CD2_2",
        ] {
            let v = dict.get(key).unwrap_or_else(|| panic!("{key} missing"));
            assert!(v.parse::<f64>().is_ok_and(f64::is_finite), "{key}={v}");
        }
        assert_eq!(dict["CRPIX1"].parse::<f64>().unwrap(), 1450.5);
        assert_eq!(dict["CMDLINE"], "arcsec -f x.fits -fov 1");
    }

    #[test]
    fn unsolved_ini_says_so() {
        let path = std::env::temp_dir().join(format!("arcsec_unini_{}.ini", std::process::id()));
        write_unsolved_ini_file(&path, "arcsec -f x.fits").unwrap();
        let dict = read_like_nina(&path);
        std::fs::remove_file(&path).ok();
        assert_eq!(dict["PLTSOLVD"], "F");
    }
}
