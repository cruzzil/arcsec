//! Cooperative cancellation of a running solve.
//!
//! A solve can take seconds (a wide spiral search, a blind index pass), and a host
//! application — an imaging suite with a Stop button — needs to end one early
//! without killing a thread. Cancellation here is cooperative: the solver polls a
//! [`CancelToken`] at its natural checkpoints (each spiral position, each blind
//! pass, each index hypothesis) and returns [`ArcsecError::Cancelled`] at the next
//! one after the token fires. Polling is one relaxed atomic load, so the cost is
//! nothing measurable.
//!
//! The token is *ambient* rather than a parameter: [`with_token`] installs it for
//! the duration of a closure on the calling thread, and the solvers read it with
//! [`current`] when they start, then hand it to any worker threads they spawn. That
//! keeps every existing signature (and every `SolveParams` literal) unchanged, and a
//! call made without a token behaves exactly as before.
//!
//! ```
//! use arcsec_core::cancel::{CancelToken, with_token};
//!
//! let token = CancelToken::new();
//! let stop = token.clone(); // hand this to the UI thread
//! stop.cancel();
//! let cancelled = with_token(&token, || arcsec_core::cancel::is_cancelled());
//! assert!(cancelled);
//! ```
//!
//! [`ArcsecError::Cancelled`]: crate::ArcsecError::Cancelled

use alloc::sync::Arc;
use core::cell::RefCell;
use core::fmt;
use core::sync::atomic::{AtomicBool, Ordering};

/// A poll function: returns `true` once the solve should stop.
type Poll = dyn Fn() -> bool + Send + Sync;

struct Inner {
    flag: AtomicBool,
    poll: Option<Box<Poll>>,
}

/// A shared, thread-safe "please stop" flag.
///
/// Clones share the flag: cancel any clone and every holder sees it. Optionally it
/// also polls a caller-supplied function (see [`CancelToken::with_poll`]), so a
/// host that already keeps its own stop flag need not mirror it into this one.
#[derive(Clone)]
pub struct CancelToken(Arc<Inner>);

impl fmt::Debug for CancelToken {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("CancelToken")
            .field("cancelled", &self.0.flag.load(Ordering::Relaxed))
            .field("poll", &self.0.poll.is_some())
            .finish()
    }
}

impl Default for CancelToken {
    fn default() -> Self {
        Self::new()
    }
}

impl CancelToken {
    /// A token that has not been cancelled.
    #[must_use]
    pub fn new() -> Self {
        Self(Arc::new(Inner {
            flag: AtomicBool::new(false),
            poll: None,
        }))
    }

    /// A token that is also cancelled once `poll` returns `true`.
    ///
    /// `poll` is called from whichever thread reaches a checkpoint, worker threads
    /// included, possibly from several at once, so it must be cheap and
    /// thread-safe. Once it has returned `true` the token stays cancelled and
    /// `poll` is not called again.
    #[must_use]
    pub fn with_poll(poll: impl Fn() -> bool + Send + Sync + 'static) -> Self {
        Self(Arc::new(Inner {
            flag: AtomicBool::new(false),
            poll: Some(Box::new(poll)),
        }))
    }

    /// Ask every solve using this token to stop at its next checkpoint.
    pub fn cancel(&self) {
        self.0.flag.store(true, Ordering::Relaxed);
    }

    /// Whether the token has been cancelled (or its poll function says so).
    #[must_use]
    pub fn is_cancelled(&self) -> bool {
        if self.0.flag.load(Ordering::Relaxed) {
            return true;
        }
        if let Some(poll) = &self.0.poll
            && poll()
        {
            self.0.flag.store(true, Ordering::Relaxed);
            return true;
        }
        false
    }
}

std::thread_local! {
    static CURRENT: RefCell<Option<CancelToken>> = const { RefCell::new(None) };
}

/// Restores the previous ambient token when dropped, panics included.
struct Restore(Option<CancelToken>);

impl Drop for Restore {
    fn drop(&mut self) {
        let prev = self.0.take();
        CURRENT.with(|c| *c.borrow_mut() = prev);
    }
}

/// Run `f` with `token` as this thread's ambient cancellation token.
///
/// Solves started inside `f` on this thread poll `token`, and pass it on to the
/// worker threads they spawn. Calls nest; the previous token is restored when `f`
/// returns or unwinds.
pub fn with_token<R>(token: &CancelToken, f: impl FnOnce() -> R) -> R {
    let prev = CURRENT.with(|c| c.borrow_mut().replace(token.clone()));
    let _restore = Restore(prev);
    f()
}

/// Like [`with_token`], but with no token at all when `token` is `None`.
pub fn with_optional<R>(token: Option<&CancelToken>, f: impl FnOnce() -> R) -> R {
    match token {
        Some(t) => with_token(t, f),
        None => f(),
    }
}

/// This thread's ambient token, if one is installed.
#[must_use]
pub fn current() -> Option<CancelToken> {
    CURRENT.with(|c| c.borrow().clone())
}

/// Whether this thread's ambient token (if any) has been cancelled.
#[must_use]
pub fn is_cancelled() -> bool {
    CURRENT.with(|c| c.borrow().as_ref().is_some_and(CancelToken::is_cancelled))
}

/// `Some(token)` cancelled, as a predicate the hot loops can capture by reference.
pub(crate) fn fired(token: Option<&CancelToken>) -> bool {
    token.is_some_and(CancelToken::is_cancelled)
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::sync::atomic::AtomicUsize;

    #[test]
    fn clones_share_the_flag() {
        let a = CancelToken::new();
        let b = a.clone();
        assert!(!a.is_cancelled());
        b.cancel();
        assert!(a.is_cancelled());
    }

    #[test]
    fn a_poll_function_cancels_and_sticks() {
        let calls = Arc::new(AtomicUsize::new(0));
        let c = Arc::clone(&calls);
        let t = CancelToken::with_poll(move || c.fetch_add(1, Ordering::Relaxed) >= 2);
        assert!(!t.is_cancelled());
        assert!(!t.is_cancelled());
        assert!(t.is_cancelled());
        assert!(t.is_cancelled());
        assert_eq!(calls.load(Ordering::Relaxed), 3, "not polled once fired");
    }

    #[test]
    fn the_ambient_token_is_scoped_and_nests() {
        assert!(current().is_none());
        let outer = CancelToken::new();
        let inner = CancelToken::new();
        inner.cancel();
        with_token(&outer, || {
            assert!(!is_cancelled());
            with_token(&inner, || assert!(is_cancelled()));
            assert!(!is_cancelled(), "outer restored");
        });
        assert!(current().is_none());
        assert!(!is_cancelled());
    }

    #[test]
    fn the_ambient_token_is_restored_after_a_panic() {
        let t = CancelToken::new();
        let r = std::panic::catch_unwind(core::panic::AssertUnwindSafe(|| {
            with_token(&t, || panic!("boom"));
        }));
        assert!(r.is_err());
        assert!(current().is_none());
    }

    #[test]
    fn other_threads_do_not_see_it() {
        let t = CancelToken::new();
        t.cancel();
        with_token(&t, || {
            let seen = std::thread::spawn(is_cancelled).join().unwrap();
            assert!(!seen);
        });
    }
}
