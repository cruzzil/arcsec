//! The log callback: arcsec's progress messages, delivered to C.

use alloc::ffi::CString;
use core::ffi::{c_char, c_int, c_void};
use core::sync::atomic::{AtomicI32, Ordering};
use std::sync::{Once, PoisonError, RwLock};

use log::{Level, LevelFilter, Log, Metadata, Record};

use crate::error::guard;

/// Message severities, most severe first. `arcsec_set_log_callback` passes
/// messages at or above the level it is given.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum arcsec_log_level {
    /// No messages.
    ARCSEC_LOG_OFF = 0,
    /// Errors: internal errors caught at the library boundary.
    ARCSEC_LOG_ERROR = 1,
    /// Warnings: a blind index that could not be read, a fallback to the hint.
    ARCSEC_LOG_WARN = 2,
    /// Progress: stars found, databases chosen, search positions, the solution.
    /// What `arcsec --progress` prints.
    ARCSEC_LOG_INFO = 3,
    /// Diagnostics.
    ARCSEC_LOG_DEBUG = 4,
    /// Detailed diagnostics.
    ARCSEC_LOG_TRACE = 5,
}

/// Log callback: `level` is an `arcsec_log_level`, `message` a NUL-terminated
/// UTF-8 line without a trailing newline, valid only during the call.
///
/// It may be called from any thread, including the solver's worker threads and
/// several at once, so it must be thread-safe. It must not call
/// `arcsec_set_log_callback`, throw or `longjmp`.
pub type arcsec_log_fn =
    Option<unsafe extern "C" fn(user: *mut c_void, level: c_int, message: *const c_char)>;

#[derive(Clone, Copy)]
struct Sink {
    callback: unsafe extern "C" fn(*mut c_void, c_int, *const c_char),
    user: *mut c_void,
}

// SAFETY: the caller registered a callback documented to be callable from any
// thread with this user pointer; arcsec never dereferences the pointer.
unsafe impl Send for Sink {}
// SAFETY: as for Send.
unsafe impl Sync for Sink {}

static SINK: RwLock<Option<Sink>> = RwLock::new(None);
/// The most verbose level the caller asked for (an `arcsec_log_level`).
static MAX_LEVEL: AtomicI32 = AtomicI32::new(0);

struct Bridge;

impl Log for Bridge {
    fn enabled(&self, m: &Metadata<'_>) -> bool {
        m.level() <= log::max_level()
    }

    fn log(&self, r: &Record<'_>) {
        let level = match r.level() {
            Level::Error => 1,
            Level::Warn => 2,
            // The search's per-position lines: thousands in a wide search, so
            // debug here, though the CLI prints them with its other progress.
            Level::Info if r.target() == arcsec_core::pipeline::SEARCH_LOG_TARGET => 4,
            Level::Info => 3,
            Level::Debug => 4,
            Level::Trace => 5,
        };
        if level <= MAX_LEVEL.load(Ordering::Relaxed) {
            deliver(level, &r.args().to_string());
        }
    }

    fn flush(&self) {}
}

/// Hand a message to the callback, if there is one. Returns whether there was.
fn deliver(level: c_int, text: &str) -> bool {
    // Copy the sink out so the lock is not held while C runs.
    let sink = *SINK.read().unwrap_or_else(PoisonError::into_inner);
    let Some(sink) = sink else {
        return false;
    };
    let bytes: Vec<u8> = text.bytes().filter(|&b| b != 0).collect();
    let msg = CString::new(bytes).unwrap_or_default();
    // SAFETY: the registered callback, with the user pointer registered with it,
    // and a NUL-terminated string alive for the duration of the call.
    unsafe { (sink.callback)(sink.user, level, msg.as_ptr()) };
    true
}

/// Report an error message (a caught panic) through the callback, if one is set
/// and errors are enabled. Returns whether it was delivered.
pub(crate) fn emit_error(text: &str) -> bool {
    MAX_LEVEL.load(Ordering::Relaxed) >= 1 && deliver(1, text)
}

/// Send arcsec's log messages to `callback` (NULL to stop), for messages at
/// `max_level` (an `arcsec_log_level`) and more severe. Process-wide: one
/// callback serves every solver and thread, and messages from concurrent solves
/// interleave. Off until called.
///
/// Thread-safe, but a message already being delivered may still reach the
/// previous callback just after this returns; keep the old `user` data alive
/// until the solves that might log to it have finished.
///
/// # Safety
///
/// `callback` must be NULL or a function safe to call from any thread with
/// `user`, as described on `arcsec_log_fn`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn arcsec_set_log_callback(
    callback: arcsec_log_fn,
    user: *mut c_void,
    max_level: c_int,
) {
    let _ = guard(|| {
        static INSTALL: Once = Once::new();
        INSTALL.call_once(|| {
            // Fails only if a logger is already set, which in this library's own
            // copy of the log crate nothing else does.
            let _ = log::set_logger(&Bridge);
        });
        *SINK.write().unwrap_or_else(PoisonError::into_inner) =
            callback.map(|callback| Sink { callback, user });
        let filter = match (callback.is_some(), max_level) {
            (false, _) | (true, ..=0) => LevelFilter::Off,
            (true, 1) => LevelFilter::Error,
            (true, 2) => LevelFilter::Warn,
            (true, 3) => LevelFilter::Info,
            (true, 4) => LevelFilter::Debug,
            (true, _) => LevelFilter::Trace,
        };
        MAX_LEVEL.store(
            if callback.is_some() {
                max_level.clamp(0, 5)
            } else {
                0
            },
            Ordering::Relaxed,
        );
        log::set_max_level(filter);
        Ok(())
    });
}
