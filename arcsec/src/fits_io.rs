use alloc::ffi::CString;
use std::path::Path;

use libc::{c_char, c_int, c_long};
use rsfitsio::aliases::rust_api::*;
use rsfitsio::fitsio::{LONGLONG, READONLY, READWRITE, fitsfile};

use arcsec_core::types::{ImageBuffer, WcsSolution};

/// Reinterpret a byte slice as a `c_char` slice for the rsfitsio wrappers.
///
/// `libc::c_char` is `i8` on x86_64 and `u8` on aarch64, so `b"KEY\0"` literals
/// cannot be passed directly on every platform. Same size and alignment either way.
#[inline]
fn cc(b: &[u8]) -> &[c_char] {
    // Safety: c_char is i8 or u8; identical layout, and we only read.
    unsafe { core::slice::from_raw_parts(b.as_ptr() as *const c_char, b.len()) }
}

/// Read RA and DEC (in degrees) from a FITS header.
/// Tries `RA`/`DEC` first (telescope pointing, NINA style), then `CRVAL1`/`CRVAL2`.
pub fn read_fits_ra_dec(path: &Path) -> Option<(f64, f64)> {
    let path_str = path.to_str()?;
    let cpath = CString::new(path_str).ok()?;

    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    fits_open_image(
        &mut fptr,
        cc(cpath.to_bytes_with_nul()),
        READONLY,
        &mut status,
    );
    if status != 0 {
        return None;
    }
    let fp = fptr.as_deref_mut()?;

    let read_dbl = |fp: &mut fitsfile, key: &[u8]| -> Option<f64> {
        let mut val = 0.0f64;
        let mut st: c_int = 0;
        fits_read_key_dbl(fp, cc(key), &mut val, None, &mut st);
        if st == 0 { Some(val) } else { None }
    };

    // Read both before returning: a `?` here would leak the open handle, and a
    // header with no pointing is the normal case for the blind solver, so this
    // path is taken routinely rather than exceptionally.
    let ra = read_dbl(fp, b"RA\0").or_else(|| read_dbl(fp, b"CRVAL1\0"));
    let dec = read_dbl(fp, b"DEC\0").or_else(|| read_dbl(fp, b"CRVAL2\0"));

    let mut close_status: c_int = 0;
    if let Some(b) = fptr {
        fits_close_file(b, &mut close_status);
    }

    Some((ra?, dec?))
}

/// Read NAXIS1/NAXIS2 from a FITS header without loading the pixel data.
///
/// Used by `arcsec catalog recommend --like`, which only needs the field size.
pub fn read_fits_dimensions(path: &Path) -> Option<(u32, u32)> {
    let path_str = path.to_str()?;
    let cpath = CString::new(path_str).ok()?;

    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    fits_open_image(
        &mut fptr,
        cc(cpath.to_bytes_with_nul()),
        READONLY,
        &mut status,
    );
    if status != 0 {
        return None;
    }
    let fp = fptr.as_deref_mut()?;

    let read_lng = |fp: &mut fitsfile, key: &[u8]| -> Option<i64> {
        // fits_read_key_lng writes a c_long, which is i32 on Windows and i64 on
        // unix, so the destination must be c_long and widen on the way out.
        let mut val: c_long = 0;
        let mut st: c_int = 0;
        fits_read_key_lng(fp, cc(key), &mut val, None, &mut st);
        if st == 0 { Some(val as i64) } else { None }
    };
    let w = read_lng(fp, b"NAXIS1\0");
    let h = read_lng(fp, b"NAXIS2\0");

    let mut close_status: c_int = 0;
    if let Some(b) = fptr {
        fits_close_file(b, &mut close_status);
    }
    Some((w? as u32, h? as u32))
}

/// Read pixel scale (arcsec/pixel, accounting for XBINNING) from FITS header.
/// Returns `None` if FOCALLEN or XPIXSZ are missing or non-positive.
pub fn read_fits_pixel_scale(path: &Path) -> Option<f64> {
    let path_str = path.to_str()?;
    let cpath = CString::new(path_str).ok()?;

    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    fits_open_image(
        &mut fptr,
        cc(cpath.to_bytes_with_nul()),
        READONLY,
        &mut status,
    );
    if status != 0 {
        return None;
    }
    let fp = fptr.as_deref_mut()?;

    let read_dbl = |fp: &mut fitsfile, key: &[u8]| -> Option<f64> {
        let mut val = 0.0f64;
        let mut st: c_int = 0;
        fits_read_key_dbl(fp, cc(key), &mut val, None, &mut st);
        if st == 0 { Some(val) } else { None }
    };

    let focallen = read_dbl(fp, b"FOCALLEN\0");
    let xpixsz = read_dbl(fp, b"XPIXSZ\0");
    let xbinning = read_dbl(fp, b"XBINNING\0").unwrap_or(1.0);

    let mut close_status: c_int = 0;
    if let Some(b) = fptr {
        fits_close_file(b, &mut close_status);
    }

    let fl = focallen.filter(|&v| v > 0.0)?;
    let ps = xpixsz.filter(|&v| v > 0.0)?;
    // plate_scale [arcsec/px] = pixel_size_µm × binning / focal_length_mm × 206.265
    Some(ps * xbinning / fl * 206.265)
}

/// Open a FITS image and return its pixel data as an `ImageBuffer`.
pub fn read_fits_image(path: &Path) -> Result<ImageBuffer, String> {
    let path_str = path.to_str().ok_or("non-UTF-8 path")?;
    let cpath = CString::new(path_str).map_err(|e| e.to_string())?;

    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;

    // Open file — cc() adapts the &[u8] literal to the platform's c_char
    fits_open_image(
        &mut fptr,
        cc(cpath.to_bytes_with_nul()),
        READONLY,
        &mut status,
    );
    if status != 0 {
        return Err(format!("fits_open_image failed: status {status}"));
    }

    let fp = fptr.as_deref_mut().ok_or("null fptr after open")?;

    // Read image header: SIMPLE, BITPIX, NAXIS, NAXIS1/2, PCOUNT, GCOUNT, EXTEND
    let mut simple: c_int = 0;
    let mut bitpix: c_int = 0;
    let mut naxis: c_int = 0;
    let mut naxes: [c_long; 9] = [0; 9];
    let mut pcount: c_long = 0;
    let mut gcount: c_long = 0;
    let mut extend: c_int = 0;
    fits_read_imghdr(
        fp,
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
        let _ = fptr.map(|b| {
            let mut s = 0i32;
            fits_close_file(b, &mut s);
        });
        return Err(format!("fits_read_imghdr failed: status {status}"));
    }

    if naxis < 2 {
        let _ = fptr.map(|b| {
            let mut s = 0i32;
            fits_close_file(b, &mut s);
        });
        return Err(format!("FITS image has only {naxis} axes (need ≥ 2)"));
    }

    let width = naxes[0] as usize;
    let height = naxes[1] as usize;
    if width == 0 || height == 0 {
        let _ = fptr.map(|b| {
            let mut s = 0i32;
            fits_close_file(b, &mut s);
        });
        return Err(format!("invalid image dimensions {width}×{height}"));
    }

    // Re-borrow (fptr was partially moved in error paths above — here we're past those)
    let fp = fptr.as_deref_mut().ok_or("null fptr")?;

    let npix = width * height;
    let mut data: Vec<f32> = vec![0.0f32; npix];

    // fits_read_img_flt = ffgpve_safe: reads directly as f32, any BITPIX
    fits_read_img_flt(
        fp,
        1, // group 1
        1, // firstelem 1-based
        npix as LONGLONG,
        0.0f32, // null pixel replacement
        &mut data,
        None, // anynul
        &mut status,
    );

    // Close (consume the Box)
    let mut close_status: c_int = 0;
    if let Some(b) = fptr {
        fits_close_file(b, &mut close_status);
    }

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
fn fits_f64(v: f64) -> String {
    let raw = format!("{:.12E}", v);
    let (mantissa, exp_str) = raw.split_once('E').unwrap();
    let exp: i32 = exp_str.parse().unwrap_or(0);
    // {:>20} right-justifies; {:+04} gives e.g. "+002" or "-004"
    format!("{:>20}", format!("{}E{:+04}", mantissa, exp))
}

/// Build a float FITS card: `KEYWORD = <20-char value> / <comment>` (80 chars).
fn dbl_card(keyword: &str, value: f64, comment: &str) -> String {
    let card = format!("{:<8}= {} / {:<47}", keyword, fits_f64(value), comment);
    format!("{:<80}", &card[..80.min(card.len())])
}

/// Build a string FITS card: `KEYWORD = '<value>'<pad>  / <comment>` (80 chars).
fn str_card(keyword: &str, value: &str, comment: &str) -> String {
    let quoted = format!("'{value}'");
    let card = format!("{:<8}= {:<20} / {:<47}", keyword, quoted, comment);
    format!("{:<80}", &card[..80.min(card.len())])
}

/// Build a logical FITS card: `KEYWORD =                    T / <comment>` (80 chars).
fn log_card(keyword: &str, value: bool, comment: &str) -> String {
    let v = if value { "T" } else { "F" };
    let card = format!("{:<8}= {:>20} / {:<47}", keyword, v, comment);
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
        log_card("PLTSOLVD", true, "Astrometric solved by arcsec 0.1.0"),
    ];

    for card in &cards {
        writeln!(f, "{card}")?;
    }
    writeln!(f, "{:<80}", "END")?;

    Ok(())
}

/// Write an `.ini` summary file with key=value pairs (ASTAP-compatible layout).
pub fn write_ini_file(path: &Path, wcs: &WcsSolution, nstars: usize) -> std::io::Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "PLTSOLVD=T")?;
    writeln!(f, "CRVAL1={:.9}", wcs.ra0.to_degrees())?;
    writeln!(f, "CRVAL2={:.9}", wcs.dec0.to_degrees())?;
    writeln!(f, "CDELT1={:.13E}", wcs.cdelt1)?;
    writeln!(f, "CDELT2={:.13E}", wcs.cdelt2)?;
    writeln!(f, "CROTA2={:.4}", wcs.crota2)?;
    writeln!(f, "NSTARS={nstars}")?;
    writeln!(f, "NQUADS={}", wcs.stars_matched)?;
    writeln!(f, "RMS={:.4}", wcs.residual_rms)?;
    Ok(())
}

/// Write WCS keywords back into the FITS file header in-place (`--update` flag).
pub fn update_fits_wcs(path: &Path, wcs: &WcsSolution) -> Result<(), String> {
    let path_str = path.to_str().ok_or("non-UTF-8 path")?;
    let cpath = CString::new(path_str).map_err(|e| e.to_string())?;

    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;

    fits_open_image(
        &mut fptr,
        cc(cpath.to_bytes_with_nul()),
        READWRITE,
        &mut status,
    );
    if status != 0 {
        return Err(format!("fits_open_image (RW) failed: status {status}"));
    }

    {
        let fp = fptr.as_deref_mut().ok_or("null fptr")?;
        let decim: c_int = 12;

        // Helper closures capture fp + status
        let upd_s = |fp: &mut fitsfile, key: &[u8], val: &[u8], com: &[u8], st: &mut c_int| {
            fits_update_key_str(fp, cc(key), cc(val), Some(cc(com)), st);
        };
        let upd_d = |fp: &mut fitsfile, key: &[u8], val: f64, com: &[u8], st: &mut c_int| {
            fits_update_key_dbl(fp, cc(key), val, decim, Some(cc(com)), st);
        };
        let upd_l = |fp: &mut fitsfile, key: &[u8], val: bool, com: &[u8], st: &mut c_int| {
            fits_update_key_log(fp, cc(key), if val { 1 } else { 0 }, Some(cc(com)), st);
        };

        upd_s(
            fp,
            b"CTYPE1\0",
            b"RA---TAN\0",
            b"first parameter RA,    projection TANgential\0",
            &mut status,
        );
        upd_s(
            fp,
            b"CTYPE2\0",
            b"DEC--TAN\0",
            b"second parameter DEC,  projection TANgential\0",
            &mut status,
        );
        upd_s(
            fp,
            b"CUNIT1\0",
            b"deg\0",
            b"Unit of coordinates\0",
            &mut status,
        );
        upd_s(
            fp,
            b"CUNIT2\0",
            b"deg\0",
            b"Unit of coordinates\0",
            &mut status,
        );
        upd_d(
            fp,
            b"CRPIX1\0",
            wcs.crpix1,
            b"X of reference pixel\0",
            &mut status,
        );
        upd_d(
            fp,
            b"CRPIX2\0",
            wcs.crpix2,
            b"Y of reference pixel\0",
            &mut status,
        );
        upd_d(
            fp,
            b"CRVAL1\0",
            wcs.ra0.to_degrees(),
            b"RA of reference pixel (deg)\0",
            &mut status,
        );
        upd_d(
            fp,
            b"CRVAL2\0",
            wcs.dec0.to_degrees(),
            b"DEC of reference pixel (deg)\0",
            &mut status,
        );
        upd_d(
            fp,
            b"CDELT1\0",
            wcs.cdelt1.abs(),
            b"X pixel size (deg)\0",
            &mut status,
        );
        upd_d(
            fp,
            b"CDELT2\0",
            wcs.cdelt2.abs(),
            b"Y pixel size (deg)\0",
            &mut status,
        );
        upd_d(
            fp,
            b"CROTA1\0",
            wcs.crota2,
            b"Image twist of X axis (deg)\0",
            &mut status,
        );
        upd_d(
            fp,
            b"CROTA2\0",
            wcs.crota2,
            b"Image twist of Y axis (deg)\0",
            &mut status,
        );
        upd_d(
            fp,
            b"CD1_1\0",
            wcs.cd1_1,
            b"CD matrix element\0",
            &mut status,
        );
        upd_d(
            fp,
            b"CD1_2\0",
            wcs.cd1_2,
            b"CD matrix element\0",
            &mut status,
        );
        upd_d(
            fp,
            b"CD2_1\0",
            wcs.cd2_1,
            b"CD matrix element\0",
            &mut status,
        );
        upd_d(
            fp,
            b"CD2_2\0",
            wcs.cd2_2,
            b"CD matrix element\0",
            &mut status,
        );
        upd_l(
            fp,
            b"PLTSOLVD\0",
            true,
            b"Astrometric solved by arcsec 0.1.0\0",
            &mut status,
        );
    }

    let mut close_status: c_int = 0;
    if let Some(b) = fptr {
        fits_close_file(b, &mut close_status);
    }

    if status != 0 {
        return Err(format!("update_fits_wcs failed: status {status}"));
    }
    Ok(())
}
