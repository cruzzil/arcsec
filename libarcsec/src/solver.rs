//! `arcsec_solver`: the handle a solve runs on, and the solve entry points.

use alloc::sync::Arc;
use core::ffi::c_char;
use core::sync::atomic::{AtomicBool, Ordering};
use std::path::Path;
use std::sync::{Mutex, PoisonError};
use std::time::Instant;

use arcsec_core::ImageBuffer;
use arcsec_core::auto::{Plan, SolveRequest};
use arcsec_io::image_io;

use crate::error::{Failure, Outcome, busy, file_error, guard, guard_with};
use crate::image::{arcsec_image, read_image};
use crate::options::{HeaderHints, Options, arcsec_solve_options, read_options};
use crate::result::arcsec_result;
use crate::util::{path_arg, put};

/// A solver: the handle a solve runs on, and through which it can be cancelled.
///
/// It holds no catalogue or image data between solves; one per thread that
/// solves is the intended use. A handle runs one solve at a time: a second
/// `arcsec_solve` on it from another thread, while one is running, returns
/// `ARCSEC_BUSY` at once rather than waiting. `arcsec_solver_cancel` may be
/// called from any thread at any time.
pub struct arcsec_solver {
    /// Held for the duration of a solve.
    running: Mutex<()>,
    /// The cancel flag of the solve in progress, if any.
    stop: Mutex<Option<Arc<AtomicBool>>>,
}

/// Create a solver. Returns NULL only if memory is exhausted. Free it with
/// `arcsec_solver_free`.
#[unsafe(no_mangle)]
pub extern "C" fn arcsec_solver_new() -> *mut arcsec_solver {
    guard_with(core::ptr::null_mut(), || {
        Ok(Box::into_raw(Box::new(arcsec_solver {
            running: Mutex::new(()),
            stop: Mutex::new(None),
        })))
    })
    .unwrap_or_else(|null| null)
}

/// Free a solver. NULL is ignored. It must not be running a solve: free it after
/// the solve call has returned (cancel it first to make that quick).
///
/// # Safety
///
/// `solver` must be NULL or a live solver, not used again afterwards.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn arcsec_solver_free(solver: *mut arcsec_solver) {
    let _ = guard(|| {
        if !solver.is_null() {
            // SAFETY: a live solver is a Box made by into_raw, freed exactly once.
            drop(unsafe { Box::from_raw(solver) });
        }
        Ok(())
    });
}

/// Ask the solve running on `solver` to stop. It returns `ARCSEC_CANCELLED` soon
/// after (at its next checkpoint, typically within milliseconds). Thread-safe;
/// does nothing if no solve is running, and does not affect a solve started later.
///
/// # Safety
///
/// `solver` must be NULL or a live solver.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn arcsec_solver_cancel(solver: *const arcsec_solver) {
    let _ = guard(|| {
        // SAFETY: forwarded contract.
        if let Some(s) = unsafe { solver.as_ref() } {
            let stop = s.stop.lock().unwrap_or_else(PoisonError::into_inner);
            if let Some(flag) = stop.as_ref() {
                flag.store(true, Ordering::Relaxed);
            }
        }
        Ok(())
    });
}

/// Clears the solver's stop flag when a solve ends, however it ends.
struct Running<'a> {
    solver: &'a arcsec_solver,
}

impl Drop for Running<'_> {
    fn drop(&mut self) {
        *self
            .solver
            .stop
            .lock()
            .unwrap_or_else(PoisonError::into_inner) = None;
    }
}

/// Run one solve on `solver`: claim it, read the options, then `body`.
///
/// # Safety
///
/// `solver` NULL or live; `opts` as for [`read_options`]; `out` NULL or writable.
unsafe fn run_solve(
    solver: *const arcsec_solver,
    opts: *const arcsec_solve_options,
    out: *mut *mut arcsec_result,
    hints: impl FnOnce() -> Outcome<HeaderHints>,
    load: impl FnOnce() -> Outcome<(ImageBuffer, usize)>,
) -> Outcome<()> {
    // SAFETY: forwarded contract.
    unsafe { put(out, core::ptr::null_mut()) };
    if out.is_null() {
        return Err(Failure::invalid("result output pointer is NULL"));
    }
    // SAFETY: forwarded contract.
    let solver = unsafe { solver.as_ref() }.ok_or_else(|| Failure::invalid("solver is NULL"))?;
    let _claim = match solver.running.try_lock() {
        Ok(g) => g,
        Err(std::sync::TryLockError::Poisoned(p)) => p.into_inner(),
        Err(std::sync::TryLockError::WouldBlock) => return Err(busy()),
    };
    let stop = Arc::new(AtomicBool::new(false));
    *solver.stop.lock().unwrap_or_else(PoisonError::into_inner) = Some(Arc::clone(&stop));
    let _running = Running { solver };
    let t0 = Instant::now();

    // SAFETY: forwarded contract.
    let Options {
        request,
        check_pattern,
    } = unsafe { read_options(opts, hints()?, stop) }?;
    let (mut img, channels) = load()?;
    if request
        .cancel
        .as_ref()
        .is_some_and(arcsec_core::cancel::CancelToken::is_cancelled)
    {
        return Err(crate::error::cancelled());
    }
    if check_pattern {
        if channels > 1 {
            log::info!("Skipping check pattern filter. This filter works only for raw OSC images!");
        } else if img.check_pattern_filter() {
            log::info!("Applying check pattern filter.");
        }
    }
    let result = solve_buffer(&mut img, &request, t0)?;
    // SAFETY: out is non-NULL (checked) and writable per the contract.
    unsafe { put(out, Box::into_raw(Box::new(result))) };
    Ok(())
}

/// Normalise, plan and solve, as the CLI does after reading the file.
fn solve_buffer(
    img: &mut ImageBuffer,
    request: &SolveRequest,
    t0: Instant,
) -> Outcome<arcsec_result> {
    // Replace non-finite pixels and, for float data in physical units or
    // normalised to 0..1, rescale into the range the histogram background
    // estimator needs. 16-bit data is left untouched.
    arcsec_core::with_max_threads(request.threads, || img.normalize_for_detection());
    let plan = Plan::new(request, img.width, img.height)?;
    log::info!(
        "Solving a {}x{} image with star database {} for a {:.2}° field, binning {}",
        img.width,
        img.height,
        plan.params.db_name.to_uppercase(),
        plan.params.fov.to_degrees(),
        plan.binning
    );
    let solved = plan.solve(img)?;
    Ok(arcsec_result::new(
        solved,
        &plan.params.db_name,
        plan.binning,
        t0.elapsed().as_secs_f64(),
    ))
}

/// Solve an image held in memory.
///
/// On `ARCSEC_OK`, `*out` is a new result to free with `arcsec_result_free`;
/// otherwise `*out` is NULL and `arcsec_last_error` says why. `opts` may
/// be NULL for all defaults.
///
/// Blocks until the solve ends; run it on a worker thread and cancel it from
/// another with `arcsec_solver_cancel` or the options' cancel callback.
///
/// # Safety
///
/// `solver` must be NULL or a live solver; `image` must be NULL or describe
/// readable memory as documented on `arcsec_image`; `opts` must be NULL or an
/// initialised options struct; `out` must be NULL or writable.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn arcsec_solve(
    solver: *const arcsec_solver,
    image: *const arcsec_image,
    opts: *const arcsec_solve_options,
    out: *mut *mut arcsec_result,
) -> crate::arcsec_status {
    guard(|| {
        // SAFETY: forwarded contract.
        unsafe {
            run_solve(
                solver,
                opts,
                out,
                || Ok(HeaderHints::default()),
                || read_image(image),
            )
        }
    })
}

/// Solve an image file: FITS (including compressed), XISF or ASDF, told apart by
/// content. The file's header supplies the hint and pixel scale where the
/// options (and their `fits_header`) do not, exactly as the `arcsec` command
/// line reads them. Colour images are averaged to one channel, except FITS cubes,
/// which are read as their first plane (as the CLI does).
///
/// Results, errors and threading are as for `arcsec_solve`; an unreadable file
/// is `ARCSEC_FILE_ERROR`.
///
/// # Safety
///
/// `path` must be NULL or a NUL-terminated string (UTF-8; any bytes on Unix); the
/// rest as for `arcsec_solve`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn arcsec_solve_file(
    solver: *const arcsec_solver,
    path: *const c_char,
    opts: *const arcsec_solve_options,
    out: *mut *mut arcsec_result,
) -> crate::arcsec_status {
    guard(|| {
        // SAFETY: forwarded contract.
        let path = unsafe { path_arg(path, "path") }?
            .ok_or_else(|| Failure::invalid("path is NULL or empty"))?;
        let p: &Path = &path;
        // SAFETY: forwarded contract.
        unsafe {
            run_solve(
                solver,
                opts,
                out,
                || {
                    image_io::detect_format(p).map_err(file_error)?;
                    Ok(HeaderHints {
                        ra_dec: image_io::read_ra_dec(p),
                        pixel_scale: image_io::read_pixel_scale(p),
                    })
                },
                || {
                    let img = image_io::read_image(p).map_err(file_error)?;
                    Ok((img, image_io::read_channels(p)))
                },
            )
        }
    })
}
