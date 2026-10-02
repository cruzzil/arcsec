//! Core astrometry for the [arcsec] plate solver.
//!
//! Given the pixels of an astronomical image and an approximate pointing, this crate
//! works out exactly where the image lies on the sky and returns a FITS-style WCS
//! solution. It implements the star-pattern matching approach introduced by ASTAP and
//! reads ASTAP's star databases (`.1476`, `.290` and `.001`). For hint-free ("blind")
//! solving it builds its own pattern index from those databases ([`mod@index`],
//! [`pipeline::index_solve()`]), and also reads Astrometry.net index files.
//!
//! The pipeline, driven by [`pipeline::solve_image`]:
//!
//! 1. **Detection** ([`detection`]) — background and noise estimation, then a
//!    multi-pass star finder measuring centroid, HFD and SNR.
//! 2. **Patterns** ([`quads`]) — 4-star quads described by five distance ratios.
//! 3. **Search** ([`pipeline`]) — a square spiral around the hint; at each position
//!    catalogue stars ([`catalog`]) are projected onto the tangent plane
//!    ([`math::coords`]), turned into quads and matched against the image.
//! 4. **Fit and verify** ([`math::lsq`], [`wcs`]) — a least-squares plate fit,
//!    checked star by star before it is accepted.
//!
//! Angles are radians throughout the API unless a name says otherwise.
//!
//! # Example
//!
//! ```no_run
//! use std::path::PathBuf;
//! use arcsec_core::ImageBuffer;
//! use arcsec_core::pipeline::{SearchSpeed, SolveMethod, SolveParams, solve_image};
//!
//! # fn load_pixels() -> ImageBuffer { ImageBuffer::new(4096, 4096) }
//! let img: ImageBuffer = load_pixels(); // row-major f32 pixels from your FITS reader
//! let params = SolveParams {
//!     ra_hint: 83.82_f64.to_radians(),
//!     dec_hint: (-5.39_f64).to_radians(),
//!     fov: 1.2_f64.to_radians(),
//!     search_radius: 10.0_f64.to_radians(),
//!     quad_tolerance: 0.007,
//!     hfd_min: 1.5,
//!     max_stars: 500,
//!     db_path: PathBuf::from("/usr/share/astap/data"),
//!     db_name: "d50".into(),
//!     binning: 1,
//!     method: SolveMethod::Quads,
//!     threads: 0,
//!     speed: SearchSpeed::Auto,
//! };
//! let wcs = solve_image(&img, &params)?;
//! println!(
//!     "centre RA {:.4}°, Dec {:.4}°, scale {:.2}\"/px",
//!     wcs.ra0.to_degrees(),
//!     wcs.dec0.to_degrees(),
//!     wcs.cdelt2 * 3600.0
//! );
//! # Ok::<(), arcsec_core::ArcsecError>(())
//! ```
//!
//! Progress is reported through the [`log`] crate at `info` level; install any
//! logger to see it.
//!
//! [arcsec]: https://github.com/cruzzil/arcsec
//! [`log`]: https://docs.rs/log

#![warn(missing_docs)]
// Library-only API hygiene, on top of the workspace lints (which the CLI shares).
#![warn(
    clippy::must_use_candidate,
    clippy::return_self_not_must_use,
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::doc_markdown
)]

extern crate alloc;

use core::sync::atomic::{AtomicUsize, Ordering};

/// Process-wide worker-thread limit. 0 = one per available core.
static MAX_THREADS: AtomicUsize = AtomicUsize::new(0);

/// Set the maximum number of worker threads any stage may use.
///
/// Detection bands, the background histogram, the pixel-range scan and the spiral
/// search each spawn their own workers, so a single knob has to reach all of them;
/// threading a parameter through every signature would be worse. `1` makes the whole
/// solve single-threaded, which is what you want when running many solves in
/// parallel yourself, or when profiling.
pub fn set_max_threads(n: usize) {
    MAX_THREADS.store(n, Ordering::Relaxed);
}

std::thread_local! {
    /// Per-thread override of [`MAX_THREADS`]; 0 = none. See [`with_max_threads`].
    static LOCAL_MAX_THREADS: core::cell::Cell<usize> = const { core::cell::Cell::new(0) };
}

/// Resolve the thread limit: this thread's override from [`with_max_threads`], else
/// the process-wide value from [`set_max_threads`], else one per core.
#[must_use]
pub fn max_threads() -> usize {
    match LOCAL_MAX_THREADS.with(core::cell::Cell::get) {
        0 => match MAX_THREADS.load(Ordering::Relaxed) {
            0 => std::thread::available_parallelism().map_or(1, core::num::NonZero::get),
            n => n,
        },
        n => n,
    }
}

/// Run `f` with the thread limit set to `n` on this thread only (0 = no override).
///
/// [`set_max_threads`] is process-wide, which is right for a command-line tool and
/// wrong for a library host running several solves at once with different
/// budgets. Every stage reads the limit on the thread that called the solver, and
/// the solvers that hand work to their own threads pass the override on, so a
/// solve inside `f` keeps to `n` threads whatever the rest of the process does.
/// The previous value is restored when `f` returns or unwinds.
pub fn with_max_threads<R>(n: usize, f: impl FnOnce() -> R) -> R {
    struct Restore(usize);
    impl Drop for Restore {
        fn drop(&mut self) {
            LOCAL_MAX_THREADS.with(|c| c.set(self.0));
        }
    }
    let _restore = Restore(LOCAL_MAX_THREADS.with(|c| c.replace(n)));
    f()
}

/// This thread's override from [`with_max_threads`], 0 if none: what a solver
/// passes on to a thread it spawns.
#[must_use]
pub fn local_max_threads() -> usize {
    LOCAL_MAX_THREADS.with(core::cell::Cell::get)
}

pub mod auto;
pub mod cancel;
pub mod catalog;
pub mod detection;
pub mod error;
pub mod index;
pub mod math;
pub mod pipeline;
pub mod quads;
pub mod types;
pub mod wcs;

#[cfg(test)]
mod test_support;

pub use catalog::{AnetIndex, AnetIndexEntry, AnetStar, load_anet_index, peek_anet_scale};
pub use error::{ArcsecError, Result};
pub use pipeline::{BlindSolveParams, blind_solve};
pub use types::{
    ImageBuffer, MatchedStar, PairedPositions, PlateConstants, Quad, QuadList, Star, StarList,
    WcsSolution,
};
