//! `arcsec_solve_options`: what the caller knows and wants, read into an
//! `arcsec_core::auto::SolveRequest`.

use core::ffi::{c_char, c_int, c_void};

use arcsec_core::auto::{ScaleSearch, SolveRequest};
use arcsec_core::cancel::CancelToken;
use arcsec_core::pipeline::{SearchSpeed, SolveMethod};
use arcsec_io::header::HeaderCards;

use crate::error::{Failure, Outcome, guard};
use crate::util::{Versioned, path_arg, read_versioned, str_arg};

/// `arcsec_solve_options::method`: ASTAP-style four-star quads (the default).
pub const ARCSEC_METHOD_QUADS: u32 = 0;
/// `arcsec_solve_options::method`: three-star triangles.
pub const ARCSEC_METHOD_TETRA: u32 = 1;

/// Progress callback: `fraction` in [0, 1] through the current `stage`, or -1
/// when the stage's extent is unknown. Stages: "detecting stars", "searching"
/// (the fraction of search positions within the radius started so far; a solve
/// usually ends well before 1), "blind index". `stage` is a static string.
///
/// Called from the solver's worker threads, possibly several at once, and at most
/// a few hundred times per solve, so it must be thread-safe. It must not call back
/// into arcsec, throw or `longjmp`.
pub type arcsec_progress_fn =
    Option<unsafe extern "C" fn(user: *mut c_void, fraction: f64, stage: *const c_char)>;

/// Cancellation callback: return nonzero to stop the solve.
///
/// Called often (at every search position) and from the solver's worker threads,
/// possibly several at once, so it must be fast and thread-safe — reading a flag
/// is ideal. It must not call back into arcsec, throw or `longjmp`.
pub type arcsec_cancel_fn = Option<unsafe extern "C" fn(user: *mut c_void) -> c_int>;

/// How to solve. Initialise with `arcsec_solve_options_init`, which sets
/// `struct_size` and every default, then change what you need. Strings are
/// UTF-8, NUL-terminated, read during the call only, and NULL (or "") for unset.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct arcsec_solve_options {
    /// `sizeof(arcsec_solve_options)`; set by the init function.
    pub struct_size: usize,

    /// Nonzero: `ra_deg`/`dec_deg` are the approximate image centre. Zero: no
    /// hint (from these fields; the FITS header may still give one).
    pub has_hint: i32,
    /// Approximate right ascension of the image centre, degrees.
    pub ra_deg: f64,
    /// Approximate declination of the image centre, degrees.
    pub dec_deg: f64,
    /// How far from the hint to search, degrees. Default 180 (the whole sky). 0
    /// tries the hint position only.
    pub search_radius_deg: f64,

    /// Field of view along the image *height*, degrees (ASTAP's `-fov`); 0 if
    /// unknown. Takes precedence over `pixel_scale_arcsec`.
    pub fov_deg: f64,
    /// Pixel scale of the image as passed (unbinned), arcseconds per pixel; 0 if
    /// unknown. With neither this nor `fov_deg` nor header optics, 1″/px is
    /// assumed, which is usually wrong: pass the scale if you know it.
    pub pixel_scale_arcsec: f64,
    /// Optional FITS header text: 80-character cards run together (as CFITSIO's
    /// `fits_hdr2str` writes them) or one per line. Used only for what the fields
    /// above leave unset: the hint from RA/DEC (or CRVAL1/2), the pixel scale from
    /// FOCALLEN, XPIXSZ and XBINNING — the keywords the CLI reads from a file.
    pub fits_header: *const c_char,

    /// Star database directory. NULL: the directory `arcsec catalog install`
    /// uses (see `arcsec_default_catalog_dir`) if it holds a database,
    /// else the working directory.
    pub catalog_dir: *const c_char,
    /// Star database name (`"d50"`, `"g05"`, ...). NULL: the densest installed
    /// database whose range suits the field size.
    pub database: *const c_char,
    /// Blind index to use: an arcsec `.arcsecix` file, a directory holding one,
    /// or a directory of Astrometry.net `index-*.fits` files; NULL for none (but
    /// see `auto_index`). Without a hint it is tried first and solves blind. With
    /// one, it is a fallback: the solver first searches a few fields round the
    /// hint (five fields, at least 1°, within `search_radius_deg`) and consults
    /// the index only if that finds nothing — unless `index_first` is set.
    pub index_path: *const c_char,
    /// Nonzero: with a hint, try `index_path` before searching round the hint, as
    /// `arcsec --index` does. Default 0.
    pub index_first: i32,
    /// Nonzero (default): when the search radius is wide (at least 10° and more
    /// than five fields), consult an arcsec index installed in the catalogue
    /// directory, as the CLI does. Zero: never use an index not named in
    /// `index_path`.
    pub auto_index: i32,

    /// SIP distortion order wanted: 0 or 1 (default) for a linear solution, 2 or
    /// more to fit SIP polynomials. arcsec currently fits order 3 whatever the
    /// value, and returns them only when the distortion is significant; the
    /// result's `sip_order` says what was fitted.
    pub sip_order: i32,
    /// Most image stars to use. Default 500.
    pub max_stars: u32,
    /// Binning factor before solving; 0 (default) chooses one from the scale.
    pub downsample: u32,
    /// Worker threads for this solve; 0 (default) uses one per core.
    pub threads: u32,
    /// `ARCSEC_METHOD_QUADS` (default) or `ARCSEC_METHOD_TETRA`.
    pub method: u32,
    /// Nonzero: read a catalogue window twice the field at every search position
    /// (ASTAP's `-speed slow`). Slower; helps fields with many stars.
    pub slow: i32,
    /// Nonzero: even out a raw one-shot-colour (Bayer) mosaic first (ASTAP's
    /// `-check y`). Ignored for multi-channel images.
    pub check_pattern: i32,
    /// Pattern-matching tolerance. Default 0.007.
    pub quad_tolerance: f64,
    /// Smallest star to use, half-flux diameter in arcseconds. Default 1.5.
    pub hfd_min_arcsec: f64,

    /// Optional cancellation callback, polled during the solve.
    pub cancel: arcsec_cancel_fn,
    /// Passed to `cancel`.
    pub cancel_user: *mut c_void,
    /// Optional progress callback.
    pub progress: arcsec_progress_fn,
    /// Passed to `progress`.
    pub progress_user: *mut c_void,
}

/// Size of the first ABI's options struct.
pub(crate) const OPTIONS_V1_SIZE: usize = core::mem::size_of::<arcsec_solve_options>();

// SAFETY: repr(C), starts with struct_size, and every field (numbers, raw
// pointers, an optional function pointer) is valid for any bit pattern.
unsafe impl Versioned for arcsec_solve_options {
    fn defaults() -> Self {
        Self {
            struct_size: core::mem::size_of::<Self>(),
            has_hint: 0,
            ra_deg: 0.0,
            dec_deg: 0.0,
            search_radius_deg: 180.0,
            fov_deg: 0.0,
            pixel_scale_arcsec: 0.0,
            fits_header: core::ptr::null(),
            catalog_dir: core::ptr::null(),
            database: core::ptr::null(),
            index_path: core::ptr::null(),
            index_first: 0,
            auto_index: 1,
            sip_order: 0,
            max_stars: 500,
            downsample: 0,
            threads: 0,
            method: ARCSEC_METHOD_QUADS,
            slow: 0,
            check_pattern: 0,
            quad_tolerance: 0.007,
            hfd_min_arcsec: 1.5,
            cancel: None,
            cancel_user: core::ptr::null_mut(),
            progress: None,
            progress_user: core::ptr::null_mut(),
        }
    }
}

/// Fill `opts` with the defaults (those of the `arcsec` command line) and set its
/// `struct_size`. Returns `ARCSEC_INVALID_ARGUMENT` if `opts` is NULL.
///
/// # Safety
///
/// `opts` must be NULL or point to a writable `arcsec_solve_options`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn arcsec_solve_options_init(
    opts: *mut arcsec_solve_options,
) -> crate::arcsec_status {
    guard(|| {
        if opts.is_null() {
            return Err(Failure::invalid("options is NULL"));
        }
        // SAFETY: non-NULL and writable per the contract.
        unsafe { opts.write_unaligned(arcsec_solve_options::defaults()) };
        Ok(())
    })
}

/// A user pointer handed back to C on other threads. The caller promised, by
/// passing a callback documented as thread-safe, that this is fine.
#[derive(Clone, Copy)]
struct UserPtr(*mut c_void);
// SAFETY: see above; arcsec only passes the pointer back, never dereferences it.
unsafe impl Send for UserPtr {}
// SAFETY: as for Send.
unsafe impl Sync for UserPtr {}

/// Options read and checked, still owning nothing of the caller's.
pub(crate) struct Options {
    pub(crate) request: SolveRequest,
    pub(crate) check_pattern: bool,
}

/// Header-derived fallbacks for what the options leave unset.
#[derive(Default)]
pub(crate) struct HeaderHints {
    /// Pointing, degrees.
    pub(crate) ra_dec: Option<(f64, f64)>,
    /// Pixel scale, arcseconds per pixel.
    pub(crate) pixel_scale: Option<f64>,
}

impl HeaderHints {
    fn from_cards(h: &HeaderCards) -> Self {
        Self {
            ra_dec: h.ra_dec(),
            pixel_scale: h.pixel_scale(),
        }
    }

    /// `self`, with gaps filled from `other`.
    pub(crate) fn or(self, other: Self) -> Self {
        Self {
            ra_dec: self.ra_dec.or(other.ra_dec),
            pixel_scale: self.pixel_scale.or(other.pixel_scale),
        }
    }
}

/// Read and check the caller's options. `file_hints` are what an image file's own
/// header says (for `arcsec_solve_file`); the options' `fits_header`
/// comes before them.
///
/// `stop` is the solver's own cancel flag, polled with the caller's callback.
///
/// # Safety
///
/// `opts` must be NULL or point to an initialised `arcsec_solve_options` whose
/// string pointers are NULL or NUL-terminated.
pub(crate) unsafe fn read_options(
    opts: *const arcsec_solve_options,
    file_hints: HeaderHints,
    stop: alloc::sync::Arc<core::sync::atomic::AtomicBool>,
) -> Outcome<Options> {
    let o = if opts.is_null() {
        // NULL means "all defaults": convenient, and not ambiguous.
        arcsec_solve_options::defaults()
    } else {
        // SAFETY: forwarded contract.
        unsafe { read_versioned(opts, "arcsec_solve_options", OPTIONS_V1_SIZE) }?
    };
    let finite_nonneg = |v: f64| v.is_finite() && v >= 0.0;

    // SAFETY (for the string reads below): forwarded contract.
    let header = unsafe { str_arg(o.fits_header, "fits_header") }?
        .map(|t| HeaderHints::from_cards(&HeaderCards::parse(t)))
        .unwrap_or_default()
        .or(file_hints);

    let hint = if o.has_hint != 0 {
        if !(o.ra_deg.is_finite() && o.dec_deg.is_finite() && o.dec_deg.abs() <= 90.0) {
            return Err(Failure::invalid(format!(
                "hint RA {} Dec {} is not a sky position",
                o.ra_deg, o.dec_deg
            )));
        }
        Some((o.ra_deg, o.dec_deg))
    } else {
        header.ra_dec
    };
    if !finite_nonneg(o.search_radius_deg) {
        return Err(Failure::invalid(format!(
            "search_radius_deg {} must be zero or positive",
            o.search_radius_deg
        )));
    }
    if !finite_nonneg(o.fov_deg) || !finite_nonneg(o.pixel_scale_arcsec) {
        return Err(Failure::invalid(
            "fov_deg and pixel_scale_arcsec must be zero (unknown) or positive",
        ));
    }
    if !(o.quad_tolerance.is_finite() && o.quad_tolerance > 0.0 && o.quad_tolerance < 1.0) {
        return Err(Failure::invalid(format!(
            "quad_tolerance {} is out of range",
            o.quad_tolerance
        )));
    }
    if !finite_nonneg(o.hfd_min_arcsec) {
        return Err(Failure::invalid("hfd_min_arcsec must not be negative"));
    }
    if o.max_stars == 0 {
        return Err(Failure::invalid("max_stars must be at least 1"));
    }
    let method = match o.method {
        ARCSEC_METHOD_QUADS => SolveMethod::Quads,
        ARCSEC_METHOD_TETRA => SolveMethod::Tetra,
        m => return Err(Failure::invalid(format!("unknown method {m}"))),
    };
    let fov_height = (o.fov_deg > 0.0).then(|| o.fov_deg.to_radians());
    let pixel_scale = if o.pixel_scale_arcsec > 0.0 {
        Some(o.pixel_scale_arcsec)
    } else if fov_height.is_some() {
        None
    } else {
        header.pixel_scale
    };

    let user = UserPtr(o.cancel_user);
    let cancel = match o.cancel {
        Some(cb) => CancelToken::with_poll(move || {
            let user = user;
            // SAFETY: the caller supplied a thread-safe callback taking this
            // pointer (see arcsec_cancel_fn).
            stop.load(core::sync::atomic::Ordering::Relaxed) || unsafe { cb(user.0) } != 0
        }),
        None => CancelToken::with_poll(move || stop.load(core::sync::atomic::Ordering::Relaxed)),
    };
    let cancel = match o.progress {
        Some(cb) => {
            let user = UserPtr(o.progress_user);
            cancel.with_progress(move |stage, fraction| {
                let user = user;
                let text = alloc::ffi::CString::new(stage).unwrap_or_default();
                // SAFETY: the caller supplied a thread-safe callback taking this
                // pointer (see arcsec_progress_fn); the string outlives the call.
                unsafe { cb(user.0, fraction, text.as_ptr()) };
            })
        }
        None => cancel,
    };

    let request = SolveRequest {
        hint: hint.map(|(ra, dec)| (ra.to_radians(), dec.to_radians())),
        fov_height,
        pixel_scale,
        scale_search: ScaleSearch::IfUnknown,
        search_radius: o.search_radius_deg.to_radians(),
        downsample: Some(o.downsample as usize),
        // SAFETY: forwarded contract.
        db_path: unsafe { path_arg(o.catalog_dir, "catalog_dir") }?,
        // SAFETY: forwarded contract.
        db_name: unsafe { str_arg(o.database, "database") }?
            .filter(|s| !s.is_empty())
            .map(str::to_ascii_lowercase),
        // SAFETY: forwarded contract.
        index: unsafe { path_arg(o.index_path, "index_path") }?,
        index_first: o.index_first != 0,
        auto_index: o.auto_index != 0,
        hfd_min_arcsec: o.hfd_min_arcsec,
        quad_tolerance: o.quad_tolerance,
        max_stars: o.max_stars as usize,
        method,
        speed: if o.slow != 0 {
            SearchSpeed::Slow
        } else {
            SearchSpeed::Auto
        },
        threads: o.threads as usize,
        sip: o.sip_order >= 2,
        cancel: Some(cancel),
    };
    Ok(Options {
        request,
        check_pattern: o.check_pattern != 0,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::sync::Arc;
    use core::sync::atomic::AtomicBool;

    fn stop() -> Arc<AtomicBool> {
        Arc::new(AtomicBool::new(false))
    }

    fn init() -> arcsec_solve_options {
        let mut o = core::mem::MaybeUninit::<arcsec_solve_options>::uninit();
        assert_eq!(
            unsafe { arcsec_solve_options_init(o.as_mut_ptr()) },
            crate::arcsec_status::ARCSEC_OK
        );
        unsafe { o.assume_init() }
    }

    #[test]
    fn defaults_are_the_cli_defaults() {
        let o = init();
        assert_eq!(o.struct_size, core::mem::size_of::<arcsec_solve_options>());
        let r = unsafe { read_options(&raw const o, HeaderHints::default(), stop()) }
            .unwrap()
            .request;
        let cli = SolveRequest::default();
        assert_eq!(r.hint, cli.hint);
        assert_eq!(r.search_radius, cli.search_radius);
        assert_eq!(r.max_stars, cli.max_stars);
        assert_eq!(r.quad_tolerance, cli.quad_tolerance);
        assert_eq!(r.hfd_min_arcsec, cli.hfd_min_arcsec);
        assert!(r.auto_index && !r.sip);
        assert_eq!(r.downsample, Some(0));
        // NULL options are the defaults too.
        assert!(unsafe { read_options(core::ptr::null(), HeaderHints::default(), stop()) }.is_ok());
        assert_eq!(
            unsafe { arcsec_solve_options_init(core::ptr::null_mut()) },
            crate::arcsec_status::ARCSEC_INVALID_ARGUMENT
        );
    }

    #[test]
    fn explicit_values_beat_the_header() {
        let header = c"RA      = 10.0\nDEC     = 20.0\nFOCALLEN= 500\nXPIXSZ  = 5\n";
        let mut o = init();
        o.fits_header = header.as_ptr();
        let r = unsafe { read_options(&raw const o, HeaderHints::default(), stop()) }
            .unwrap()
            .request;
        let (ra, dec) = r.hint.unwrap();
        assert!((ra.to_degrees() - 10.0).abs() < 1e-12 && (dec.to_degrees() - 20.0).abs() < 1e-12);
        assert!((r.pixel_scale.unwrap() - 5.0 / 500.0 * 206.265).abs() < 1e-12);

        o.has_hint = 1;
        o.ra_deg = 30.0;
        o.dec_deg = -40.0;
        o.fov_deg = 1.5;
        let r = unsafe { read_options(&raw const o, HeaderHints::default(), stop()) }
            .unwrap()
            .request;
        assert!((r.hint.unwrap().0.to_degrees() - 30.0).abs() < 1e-12);
        assert!(
            r.pixel_scale.is_none(),
            "an explicit FOV wins over header optics"
        );
        assert!((r.fov_height.unwrap().to_degrees() - 1.5).abs() < 1e-12);

        // The file's own header is the last resort.
        let mut o = init();
        let file = HeaderHints {
            ra_dec: Some((1.0, 2.0)),
            pixel_scale: Some(3.0),
        };
        o.pixel_scale_arcsec = 4.0;
        let r = unsafe { read_options(&raw const o, file, stop()) }
            .unwrap()
            .request;
        assert!((r.hint.unwrap().1.to_degrees() - 2.0).abs() < 1e-12);
        assert_eq!(r.pixel_scale, Some(4.0));
    }

    #[test]
    fn nonsense_is_refused() {
        let base = init();
        let mut cases: Vec<(&str, arcsec_solve_options)> = Vec::new();
        let mut o = base;
        o.has_hint = 1;
        o.dec_deg = 91.0;
        cases.push(("dec", o));
        let mut o = base;
        o.has_hint = 1;
        o.ra_deg = f64::NAN;
        cases.push(("ra", o));
        let mut o = base;
        o.search_radius_deg = -1.0;
        cases.push(("radius", o));
        let mut o = base;
        o.fov_deg = f64::INFINITY;
        cases.push(("fov", o));
        let mut o = base;
        o.pixel_scale_arcsec = -2.0;
        cases.push(("scale", o));
        let mut o = base;
        o.quad_tolerance = 0.0;
        cases.push(("tolerance", o));
        let mut o = base;
        o.max_stars = 0;
        cases.push(("stars", o));
        let mut o = base;
        o.method = 7;
        cases.push(("method", o));
        let mut o = base;
        o.struct_size = 0;
        cases.push(("struct_size", o));
        let bad_utf8 = [0xffu8, 0];
        let mut o = base;
        o.database = bad_utf8.as_ptr().cast();
        cases.push(("utf8", o));
        for (name, o) in cases {
            let r = unsafe { read_options(&raw const o, HeaderHints::default(), stop()) };
            assert!(r.is_err(), "{name} accepted");
        }
    }

    #[test]
    fn the_cancel_callback_and_flag_both_stop() {
        unsafe extern "C" fn yes(user: *mut c_void) -> c_int {
            // SAFETY: the test passes a pointer to an AtomicBool.
            c_int::from(
                unsafe { &*user.cast::<AtomicBool>() }.load(core::sync::atomic::Ordering::Relaxed),
            )
        }
        let flag = AtomicBool::new(false);
        let mut o = init();
        o.cancel = Some(yes);
        o.cancel_user = core::ptr::from_ref(&flag).cast_mut().cast();
        let s = stop();
        let r = unsafe { read_options(&raw const o, HeaderHints::default(), Arc::clone(&s)) }
            .unwrap()
            .request;
        let token = r.cancel.unwrap();
        assert!(!token.is_cancelled());
        flag.store(true, core::sync::atomic::Ordering::Relaxed);
        assert!(token.is_cancelled());

        let r = unsafe { read_options(&init(), HeaderHints::default(), Arc::clone(&s)) }
            .unwrap()
            .request;
        let token = r.cancel.unwrap();
        assert!(!token.is_cancelled());
        s.store(true, core::sync::atomic::Ordering::Relaxed);
        assert!(token.is_cancelled());
    }

    #[test]
    fn an_older_smaller_struct_gets_defaults_for_the_rest() {
        // A caller that knows only the fields up to the hint.
        let mut o = init();
        o.struct_size = core::mem::offset_of!(arcsec_solve_options, search_radius_deg);
        o.has_hint = 1;
        o.ra_deg = 5.0;
        o.max_stars = 0; // beyond its struct_size: must not be read
        // The minimum is the first ABI's full size, so this is refused...
        assert!(unsafe { read_options(&raw const o, HeaderHints::default(), stop()) }.is_err());
        // ...while a larger (newer) struct is read up to what we know.
        let mut big = [0u8; core::mem::size_of::<arcsec_solve_options>() + 64];
        let mut o = init();
        o.struct_size = big.len();
        o.max_stars = 77;
        unsafe {
            core::ptr::copy_nonoverlapping(
                core::ptr::from_ref(&o).cast::<u8>(),
                big.as_mut_ptr(),
                core::mem::size_of::<arcsec_solve_options>(),
            );
        }
        let r = unsafe { read_options(big.as_ptr().cast(), HeaderHints::default(), stop()) }
            .unwrap()
            .request;
        assert_eq!(r.max_stars, 77);
    }
}
