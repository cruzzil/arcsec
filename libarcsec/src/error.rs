//! Status codes, the per-thread error message, and the panic barrier.

use alloc::ffi::CString;
use core::cell::RefCell;
use core::ffi::{c_char, c_int};
use core::panic::AssertUnwindSafe;
use std::sync::Once;

use arcsec_core::ArcsecError;

/// What a call did. The solve codes have the values of the `arcsec` (and ASTAP)
/// command-line exit codes, so a host that already handles those can map them
/// one to one.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum arcsec_status {
    /// Success; for a solve, a verified solution.
    ARCSEC_OK = 0,
    /// The search finished without a verified solution (CLI exit 1).
    ARCSEC_NO_SOLUTION = 1,
    /// Too few stars were detected to solve (CLI exit 2).
    ARCSEC_INSUFFICIENT_STARS = 2,
    /// The image file could not be read, or is not FITS, XISF or ASDF (CLI exit 16).
    ARCSEC_FILE_ERROR = 16,
    /// The star database or blind index was not found (CLI exit 32).
    ARCSEC_DATABASE_NOT_FOUND = 32,
    /// A star database or index file could not be read (CLI exit 33).
    ARCSEC_DATABASE_ERROR = 33,
    /// An argument was NULL, out of range or inconsistent, or a struct's
    /// `struct_size` was not set.
    ARCSEC_INVALID_ARGUMENT = 100,
    /// The solve was cancelled (`arcsec_solver_cancel` or the options'
    /// cancel callback).
    ARCSEC_CANCELLED = 101,
    /// The solver handle is already running a solve on another thread.
    ARCSEC_BUSY = 102,
    /// A bug in arcsec: an internal error (a Rust panic) was caught at the library
    /// boundary. The library is still usable; please report it with the message
    /// from `arcsec_last_error`.
    ARCSEC_INTERNAL_ERROR = 199,
}

use arcsec_status::{
    ARCSEC_BUSY, ARCSEC_CANCELLED, ARCSEC_DATABASE_ERROR, ARCSEC_DATABASE_NOT_FOUND,
    ARCSEC_FILE_ERROR, ARCSEC_INSUFFICIENT_STARS, ARCSEC_INTERNAL_ERROR, ARCSEC_INVALID_ARGUMENT,
    ARCSEC_NO_SOLUTION, ARCSEC_OK,
};

/// A failed call: the status to return and the message to leave for
/// `arcsec_last_error`.
#[derive(Debug)]
pub(crate) struct Failure {
    pub(crate) status: arcsec_status,
    pub(crate) message: String,
}

impl Failure {
    pub(crate) fn new(status: arcsec_status, message: impl Into<String>) -> Self {
        Self {
            status,
            message: message.into(),
        }
    }

    pub(crate) fn invalid(message: impl Into<String>) -> Self {
        Self::new(ARCSEC_INVALID_ARGUMENT, message)
    }
}

impl From<ArcsecError> for Failure {
    fn from(e: ArcsecError) -> Self {
        let status = match &e {
            ArcsecError::InsufficientStars { .. } => ARCSEC_INSUFFICIENT_STARS,
            ArcsecError::CatalogNotFound(_) | ArcsecError::IndexNotFound(_) => {
                ARCSEC_DATABASE_NOT_FOUND
            }
            ArcsecError::CatalogIo(_) => ARCSEC_DATABASE_ERROR,
            ArcsecError::InvalidParameter(_) => ARCSEC_INVALID_ARGUMENT,
            ArcsecError::Cancelled => ARCSEC_CANCELLED,
            // InsufficientQuads, Singular, BadSolution, OutsideSearchRadius, and
            // anything added later: the search ended without a solution.
            _ => ARCSEC_NO_SOLUTION,
        };
        let message = match &e {
            ArcsecError::InsufficientQuads { .. } => "no solution found".to_string(),
            ArcsecError::CatalogNotFound(p) => {
                format!("star database not found in {}", p.display())
            }
            _ => e.to_string(),
        };
        Self { status, message }
    }
}

pub(crate) type Outcome<T> = Result<T, Failure>;

std::thread_local! {
    /// The message of the last failed call on this thread.
    static LAST_ERROR: RefCell<CString> = RefCell::new(CString::default());
    /// Where and why the last panic on this thread happened, from the panic hook.
    static LAST_PANIC: RefCell<Option<String>> = const { RefCell::new(None) };
}

pub(crate) fn set_last_error(message: &str) {
    // An interior NUL cannot be represented; cut the message there.
    let bytes = message.as_bytes();
    let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
    let c = CString::new(&bytes[..end]).unwrap_or_default();
    LAST_ERROR.with(|e| *e.borrow_mut() = c);
}

/// Route panic reports to the log callback (if one is set) instead of stderr, and
/// keep the location for the error message. Installed once; a panic hook belongs
/// to this library's own copy of the Rust runtime, so it does not affect a host
/// that is itself written in Rust.
fn install_panic_hook() {
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        let previous = std::panic::take_hook();
        std::panic::set_hook(Box::new(move |info| {
            let text = format!("{info}");
            LAST_PANIC.with(|p| *p.borrow_mut() = Some(text.clone()));
            if !crate::logging::emit_error(&text) {
                previous(info);
            }
        }));
    });
}

/// The text of a panic payload.
fn payload_text(payload: &(dyn core::any::Any + Send)) -> String {
    if let Some(s) = payload.downcast_ref::<&str>() {
        (*s).to_string()
    } else if let Some(s) = payload.downcast_ref::<String>() {
        s.clone()
    } else {
        "unknown panic".to_string()
    }
}

/// Run an entry point's body: its result, or `fallback` after a panic. A failure
/// or panic leaves its message for `arcsec_last_error`.
pub(crate) fn guard_with<T>(fallback: T, body: impl FnOnce() -> Outcome<T>) -> Result<T, T> {
    install_panic_hook();
    match std::panic::catch_unwind(AssertUnwindSafe(body)) {
        Ok(Ok(v)) => Ok(v),
        Ok(Err(f)) => {
            set_last_error(&f.message);
            Err(fallback)
        }
        Err(payload) => {
            let detail = LAST_PANIC
                .with(|p| p.borrow_mut().take())
                .unwrap_or_else(|| payload_text(&*payload));
            set_last_error(&format!("internal error (please report): {detail}"));
            // Dropping a payload can itself panic; that must not escape either.
            let _ = std::panic::catch_unwind(AssertUnwindSafe(move || drop(payload)));
            Err(fallback)
        }
    }
}

/// [`guard_with`] for the common case: a body that returns only a status.
pub(crate) fn guard(body: impl FnOnce() -> Outcome<()>) -> arcsec_status {
    // The fallback is replaced below by the failure's own status; a panic keeps it.
    let mut status = ARCSEC_INTERNAL_ERROR;
    let r = guard_with((), || {
        body().inspect_err(|f| {
            status = f.status;
        })
    });
    match r {
        Ok(()) => ARCSEC_OK,
        Err(()) => status,
    }
}

/// The message describing the most recent call on this thread that failed (did
/// not return `ARCSEC_OK`, or returned NULL or an error value). An empty string
/// if none has. Successful calls leave it unchanged.
///
/// The string belongs to the library and stays valid until the next failing
/// call on the same thread; copy it if you need it longer. Never NULL.
#[unsafe(no_mangle)]
pub extern "C" fn arcsec_last_error() -> *const c_char {
    // The CString lives in this thread's slot until replaced, so the pointer
    // outlives the borrow.
    LAST_ERROR.with(|e| e.borrow().as_ptr())
}

/// A short static description of a status code (e.g. "no solution"). Unknown
/// values give "unknown status". Never NULL; do not free.
#[unsafe(no_mangle)]
pub extern "C" fn arcsec_status_string(status: c_int) -> *const c_char {
    let s: &'static str = match status {
        0 => "ok\0",
        1 => "no solution\0",
        2 => "insufficient stars\0",
        16 => "file error\0",
        32 => "database not found\0",
        33 => "database read error\0",
        100 => "invalid argument\0",
        101 => "cancelled\0",
        102 => "solver busy\0",
        199 => "internal error\0",
        _ => "unknown status\0",
    };
    s.as_ptr().cast()
}

/// The status for a solver that is already busy.
pub(crate) fn busy() -> Failure {
    Failure::new(
        ARCSEC_BUSY,
        "this solver is already running a solve on another thread",
    )
}

/// The status for an unreadable image file.
pub(crate) fn file_error(message: impl Into<String>) -> Failure {
    Failure::new(ARCSEC_FILE_ERROR, message)
}

/// The status for a cancelled call.
pub(crate) fn cancelled() -> Failure {
    Failure::new(ARCSEC_CANCELLED, "cancelled")
}
