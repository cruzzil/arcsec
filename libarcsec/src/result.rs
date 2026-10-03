//! `arcsec_result`: a solution, and what C can read out of it.

use alloc::ffi::CString;
use core::ffi::c_char;

use arcsec_core::auto::Solved;
use arcsec_core::wcs::TanWcs;
use arcsec_core::wcs::sip::SIP_TERMS;

use crate::error::{Failure, guard, guard_with};
use crate::util::{Versioned, put, write_versioned};

/// Highest SIP order an `arcsec_wcs` can hold. arcsec fits order 3 today; the
/// arrays are sized for more so that higher orders need no ABI change.
pub const ARCSEC_SIP_MAX_ORDER: usize = 9;
/// Side of the SIP coefficient arrays: `ARCSEC_SIP_MAX_ORDER + 1`.
pub const ARCSEC_SIP_SIZE: usize = ARCSEC_SIP_MAX_ORDER + 1;

/// A plate solution, owned by the library. Read it with the `arcsec_result_*`
/// functions and release it with `arcsec_result_free`. Immutable, so it may be
/// read from several threads at once.
pub struct arcsec_result {
    solved: Solved,
    database: CString,
    binning: usize,
    elapsed_seconds: f64,
}

impl arcsec_result {
    pub(crate) fn new(
        solved: Solved,
        database: &str,
        binning: usize,
        elapsed_seconds: f64,
    ) -> Self {
        Self {
            solved,
            database: CString::new(database).unwrap_or_default(),
            binning,
            elapsed_seconds,
        }
    }
}

/// The WCS of a solution, as FITS keywords. Pixel coordinates are 1-based FITS
/// pixels of the image as passed (unbinned), with row 1 the first row in FITS
/// order (see `ARCSEC_IMAGE_TOP_DOWN`).
///
/// Set `struct_size = sizeof(arcsec_wcs)` before calling `arcsec_result_wcs`.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct arcsec_wcs {
    /// `sizeof(arcsec_wcs)`, set by the caller.
    pub struct_size: usize,
    /// CRVAL1: right ascension of the reference pixel, degrees.
    pub crval1: f64,
    /// CRVAL2: declination of the reference pixel, degrees.
    pub crval2: f64,
    /// CRPIX1: reference pixel column (the image centre).
    pub crpix1: f64,
    /// CRPIX2: reference pixel row.
    pub crpix2: f64,
    /// CD1_1, degrees per pixel.
    pub cd1_1: f64,
    /// CD1_2, degrees per pixel.
    pub cd1_2: f64,
    /// CD2_1, degrees per pixel.
    pub cd2_1: f64,
    /// CD2_2, degrees per pixel.
    pub cd2_2: f64,
    /// CDELT1, degrees per pixel, carrying the parity: negative for an image with
    /// the sky's usual orientation, positive for a mirrored one. With CDELT2 and
    /// CROTA1/2 these are the old-style keywords as `astap_cli` (and the `arcsec`
    /// command line) write them; the CD matrix is the authoritative WCS.
    pub cdelt1: f64,
    /// CDELT2, degrees per pixel (positive).
    pub cdelt2: f64,
    /// CROTA1, degrees: the rotation of the image's +X axis. It differs from
    /// `crota2` only when the plate is slightly skewed.
    pub crota1: f64,
    /// CROTA2, degrees: the rotation of the +Y axis, with the FITS sign (a frame
    /// rotated 90° east of north reads +90).
    pub crota2: f64,
    /// Pixel scale, arcseconds per pixel.
    pub pixel_scale_arcsec: f64,
    /// RMS residual of the matched stars, arcseconds.
    pub rms_arcsec: f64,
    /// Stars matched to the catalogue in the final verification.
    pub matched_stars: u32,
    /// Nonzero if the image is mirrored (det CD > 0).
    pub mirrored: i32,
    /// Order of the SIP polynomials below; 0 if there are none.
    pub sip_order: i32,
    /// Reserved; 0.
    pub reserved: i32,
    /// SIP `A_p_q` at `sip_a[p][q]` (p + q ≤ `sip_order`, others 0): pixel → sky
    /// correction along x, in pixels, applied to offsets from CRPIX.
    pub sip_a: [[f64; ARCSEC_SIP_SIZE]; ARCSEC_SIP_SIZE],
    /// SIP `B_p_q`: pixel → sky correction along y.
    pub sip_b: [[f64; ARCSEC_SIP_SIZE]; ARCSEC_SIP_SIZE],
    /// SIP `AP_p_q`: sky → pixel correction along x.
    pub sip_ap: [[f64; ARCSEC_SIP_SIZE]; ARCSEC_SIP_SIZE],
    /// SIP `BP_p_q`: sky → pixel correction along y.
    pub sip_bp: [[f64; ARCSEC_SIP_SIZE]; ARCSEC_SIP_SIZE],
}

// SAFETY: repr(C), starts with struct_size, all fields numeric.
unsafe impl Versioned for arcsec_wcs {
    fn defaults() -> Self {
        Self {
            struct_size: core::mem::size_of::<Self>(),
            crval1: 0.0,
            crval2: 0.0,
            crpix1: 0.0,
            crpix2: 0.0,
            cd1_1: 0.0,
            cd1_2: 0.0,
            cd2_1: 0.0,
            cd2_2: 0.0,
            cdelt1: 0.0,
            cdelt2: 0.0,
            crota1: 0.0,
            crota2: 0.0,
            pixel_scale_arcsec: 0.0,
            rms_arcsec: 0.0,
            matched_stars: 0,
            mirrored: 0,
            sip_order: 0,
            reserved: 0,
            sip_a: [[0.0; ARCSEC_SIP_SIZE]; ARCSEC_SIP_SIZE],
            sip_b: [[0.0; ARCSEC_SIP_SIZE]; ARCSEC_SIP_SIZE],
            sip_ap: [[0.0; ARCSEC_SIP_SIZE]; ARCSEC_SIP_SIZE],
            sip_bp: [[0.0; ARCSEC_SIP_SIZE]; ARCSEC_SIP_SIZE],
        }
    }
}

/// How a solve went, beyond the WCS. Set `struct_size` before the call.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct arcsec_solve_info {
    /// `sizeof(arcsec_solve_info)`, set by the caller.
    pub struct_size: usize,
    /// Binning factor the image was solved at.
    pub binning: u32,
    /// Pattern matches found before outlier rejection.
    pub raw_matches: u32,
    /// Distance from the hint to the solved centre, degrees.
    pub search_distance_deg: f64,
    /// Faintest catalogue magnitude used.
    pub mag_limit: f64,
    /// Wall-clock time of the solve, seconds.
    pub elapsed_seconds: f64,
    /// Nonzero if an Astrometry.net index estimated the position first.
    pub has_index_estimate: i32,
    /// Reserved; 0.
    pub reserved: i32,
    /// That estimate's right ascension, degrees.
    pub index_ra_deg: f64,
    /// That estimate's declination, degrees.
    pub index_dec_deg: f64,
}

// SAFETY: repr(C), starts with struct_size, all fields numeric.
unsafe impl Versioned for arcsec_solve_info {
    fn defaults() -> Self {
        Self {
            struct_size: core::mem::size_of::<Self>(),
            binning: 0,
            raw_matches: 0,
            search_distance_deg: 0.0,
            mag_limit: 0.0,
            elapsed_seconds: 0.0,
            has_index_estimate: 0,
            reserved: 0,
            index_ra_deg: 0.0,
            index_dec_deg: 0.0,
        }
    }
}

/// A detected star and the catalogue star it was matched to.
#[repr(C)]
#[derive(Debug, Clone, Copy, Default)]
pub struct arcsec_matched_star {
    /// Column, 1-based FITS pixels of the unbinned image.
    pub x: f64,
    /// Row, 1-based FITS pixels.
    pub y: f64,
    /// Catalogue right ascension, degrees.
    pub ra_deg: f64,
    /// Catalogue declination, degrees.
    pub dec_deg: f64,
}

/// Borrow the result behind a C pointer.
///
/// # Safety
///
/// `r` must be NULL or a live pointer from a solve function.
unsafe fn borrow<'a>(r: *const arcsec_result) -> Result<&'a arcsec_result, Failure> {
    // SAFETY: per the contract, a non-NULL r is a live Box<arcsec_result>.
    unsafe { r.as_ref() }.ok_or_else(|| Failure::invalid("result is NULL"))
}

/// Copy the solution's WCS into `out`, whose `struct_size` you have set.
///
/// # Safety
///
/// `result` must be NULL or a live result; `out` must be NULL or point to
/// `out->struct_size` writable bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn arcsec_result_wcs(
    result: *const arcsec_result,
    out: *mut arcsec_wcs,
) -> crate::arcsec_status {
    guard(|| {
        // SAFETY: forwarded contract.
        let r = unsafe { borrow(result) }?;
        let w = &r.solved.wcs;
        let mut v = arcsec_wcs::defaults();
        v.crval1 = w.ra0.to_degrees();
        v.crval2 = w.dec0.to_degrees();
        v.crpix1 = w.crpix1;
        v.crpix2 = w.crpix2;
        v.cd1_1 = w.cd1_1;
        v.cd1_2 = w.cd1_2;
        v.cd2_1 = w.cd2_1;
        v.cd2_2 = w.cd2_2;
        v.cdelt1 = w.cdelt1;
        v.cdelt2 = w.cdelt2;
        v.crota1 = w.crota1();
        v.crota2 = w.crota2;
        v.pixel_scale_arcsec = (w.cd1_1 * w.cd2_2 - w.cd1_2 * w.cd2_1).abs().sqrt() * 3600.0;
        v.rms_arcsec = w.residual_rms;
        v.matched_stars = u32::try_from(w.stars_matched).unwrap_or(u32::MAX);
        v.mirrored = i32::from(w.cd1_1 * w.cd2_2 - w.cd1_2 * w.cd2_1 > 0.0);
        if let Some(sip) = &w.sip {
            v.sip_order = i32::try_from(arcsec_core::wcs::sip::SIP_ORDER).unwrap_or(0);
            for (k, &(p, q)) in SIP_TERMS.iter().enumerate() {
                let (p, q) = (p as usize, q as usize);
                v.sip_a[p][q] = sip.a[k];
                v.sip_b[p][q] = sip.b[k];
                v.sip_ap[p][q] = sip.ap[k];
                v.sip_bp[p][q] = sip.bp[k];
            }
        }
        // SAFETY: forwarded contract.
        unsafe { write_versioned(out, &v, "arcsec_wcs", core::mem::size_of::<arcsec_wcs>()) }
    })
}

/// Copy facts about the solve into `out`, whose `struct_size` you have set.
///
/// # Safety
///
/// As `arcsec_result_wcs`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn arcsec_result_info(
    result: *const arcsec_result,
    out: *mut arcsec_solve_info,
) -> crate::arcsec_status {
    guard(|| {
        // SAFETY: forwarded contract.
        let r = unsafe { borrow(result) }?;
        let w = &r.solved.wcs;
        let mut v = arcsec_solve_info::defaults();
        v.binning = u32::try_from(r.binning).unwrap_or(u32::MAX);
        v.raw_matches = u32::try_from(w.raw_matches).unwrap_or(u32::MAX);
        v.search_distance_deg = w.search_dist_deg;
        v.mag_limit = w.mag_limit;
        v.elapsed_seconds = r.elapsed_seconds;
        if let Some((ra, dec)) = r.solved.index_estimate {
            v.has_index_estimate = 1;
            v.index_ra_deg = ra.to_degrees();
            v.index_dec_deg = dec.to_degrees();
        }
        // SAFETY: forwarded contract.
        unsafe {
            write_versioned(
                out,
                &v,
                "arcsec_solve_info",
                core::mem::size_of::<arcsec_solve_info>(),
            )
        }
    })
}

/// The star database the solve used (e.g. `"d50"`), valid until the result is
/// freed. NULL if `result` is NULL.
///
/// # Safety
///
/// `result` must be NULL or a live result.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn arcsec_result_database(result: *const arcsec_result) -> *const c_char {
    guard_with(core::ptr::null(), || {
        // SAFETY: forwarded contract.
        Ok(unsafe { borrow(result) }?.database.as_ptr())
    })
    .unwrap_or_else(|null| null)
}

/// Copy up to `capacity` of the matched star pairs into `out` and return how many
/// there are in all (so a call with `capacity` 0 asks for the count). `out` may
/// be NULL when `capacity` is 0. Returns 0 if `result` is NULL.
///
/// # Safety
///
/// `result` must be NULL or a live result; `out` must be NULL or have room for
/// `capacity` elements.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn arcsec_result_matched_stars(
    result: *const arcsec_result,
    out: *mut arcsec_matched_star,
    capacity: usize,
) -> usize {
    guard_with(0, || {
        // SAFETY: forwarded contract.
        let r = unsafe { borrow(result) }?;
        let pairs = &r.solved.wcs.matched_stars;
        if !out.is_null() {
            for (i, m) in pairs.iter().take(capacity).enumerate() {
                let star = arcsec_matched_star {
                    x: m.x,
                    y: m.y,
                    ra_deg: m.ra.to_degrees(),
                    dec_deg: m.dec.to_degrees(),
                };
                // SAFETY: i < capacity elements fit in out.
                unsafe { out.add(i).write_unaligned(star) };
            }
        }
        Ok(pairs.len())
    })
    .unwrap_or_else(|zero| zero)
}

/// Map a 1-based FITS pixel to the sky (degrees, RA in [0, 360)), through the
/// solution including its SIP terms.
///
/// # Safety
///
/// `result` must be NULL or a live result; `ra_deg` and `dec_deg` must be NULL or
/// writable.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn arcsec_result_pixel_to_sky(
    result: *const arcsec_result,
    x: f64,
    y: f64,
    ra_deg: *mut f64,
    dec_deg: *mut f64,
) -> crate::arcsec_status {
    guard(|| {
        // SAFETY: forwarded contract.
        let r = unsafe { borrow(result) }?;
        if !(x.is_finite() && y.is_finite()) {
            return Err(Failure::invalid("pixel coordinates must be finite"));
        }
        let (ra, dec) = TanWcs::from(&r.solved.wcs).pixel_to_sky(x, y);
        // SAFETY: forwarded contract.
        unsafe {
            put(ra_deg, ra.to_degrees());
            put(dec_deg, dec.to_degrees());
        }
        Ok(())
    })
}

/// Map a sky position (degrees) to a 1-based FITS pixel through the solution.
/// `ARCSEC_INVALID_ARGUMENT` for a position on the far side of the sky, which has
/// no pixel.
///
/// # Safety
///
/// As `arcsec_result_pixel_to_sky`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn arcsec_result_sky_to_pixel(
    result: *const arcsec_result,
    ra_deg: f64,
    dec_deg: f64,
    x: *mut f64,
    y: *mut f64,
) -> crate::arcsec_status {
    guard(|| {
        // SAFETY: forwarded contract.
        let r = unsafe { borrow(result) }?;
        let (px, py) = TanWcs::from(&r.solved.wcs)
            .sky_to_pixel(ra_deg.to_radians(), dec_deg.to_radians())
            .ok_or_else(|| {
                Failure::invalid("that position is not on the image's side of the sky")
            })?;
        // SAFETY: forwarded contract.
        unsafe {
            put(x, px);
            put(y, py);
        }
        Ok(())
    })
}

/// Write the solution as FITS header cards into `buf`, as `snprintf` does (at
/// most `len - 1` bytes and a NUL; returns the full length, so a value `>= len`
/// means `buf` was too small — call with `len` 0 to size it). The cards are
/// 80-character records run together without newlines, as CFITSIO's
/// `fits_hdr2str` makes them, ending with `END`: CTYPE1/2 (`RA---TAN`, or
/// `RA---TAN-SIP` with SIP), CUNIT1, CRPIX1/2, CRVAL1/2, CDELT1/2, CROTA1/2,
/// CD1_1..CD2_2, the SIP keywords if fitted, and PLTSOLVD — the content of the
/// `.wcs` file `arcsec` and ASTAP write, ready for wcslib's `wcspih`. Returns 0
/// if `result` is NULL.
///
/// # Safety
///
/// `result` must be NULL or a live result; `buf` NULL or `len` writable bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn arcsec_result_fits_header(
    result: *const arcsec_result,
    buf: *mut c_char,
    len: usize,
) -> usize {
    guard_with(0, || {
        // SAFETY: forwarded contract.
        let r = unsafe { borrow(result) }?;
        let mut text = arcsec_io::fits_io::wcs_header_cards(&r.solved.wcs).concat();
        text.push_str(&format!("{:<80}", "END"));
        // SAFETY: forwarded contract.
        Ok(unsafe { crate::util::write_c_string(text.as_bytes(), buf, len) })
    })
    .unwrap_or_else(|zero| zero)
}

/// Free a result. NULL is ignored.
///
/// # Safety
///
/// `result` must be NULL or a live result, not used again afterwards.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn arcsec_result_free(result: *mut arcsec_result) {
    let _ = guard(|| {
        if !result.is_null() {
            // SAFETY: a live result is a Box made by into_raw, freed exactly once.
            drop(unsafe { Box::from_raw(result) });
        }
        Ok(())
    });
}
