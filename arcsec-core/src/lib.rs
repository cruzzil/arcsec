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

/// Resolve the thread limit: the configured value, or one per core if unset.
pub fn max_threads() -> usize {
    match MAX_THREADS.load(Ordering::Relaxed) {
        0 => std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1),
        n => n,
    }
}

pub mod catalog;
pub mod detection;
pub mod error;
pub mod math;
pub mod pipeline;
pub mod quads;
pub mod types;
pub mod wcs;

pub use catalog::{AnetIndex, AnetIndexEntry, AnetStar, load_anet_index, peek_anet_scale};
pub use error::{ArcsecError, Result};
pub use pipeline::{BlindSolveParams, blind_solve};
pub use types::{
    ImageBuffer, PairedPositions, PlateConstants, Quad, QuadList, Star, StarList, WcsSolution,
};
