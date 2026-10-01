//! arcsec's own blind index: disc-anchored 4-star patterns, built from an installed
//! ASTAP star database and memory-mapped at solve time.
//!
//! * [`pattern`] — the descriptor, its hash key and the tangent-plane geometry,
//!   shared by both sides.
//! * [`build`] — the builder ([`build_index`]).
//! * [`mod@format`] — the `ARCSECIX` file ([`BlindIndex`]).
//!
//! The solver that uses it is [`mod@crate::pipeline::index_solve`]. Design and
//! measurements: `docs/offline-index.md`.

pub mod build;
pub mod format;
pub mod pattern;

pub use build::{
    BuildParams, BuildProgress, DEFAULT_TIERS, TierSpec, build_index, default_index_path,
    tier_fov_range,
};
pub use format::{BlindIndex, BuiltIndex, TierInfo, is_blind_index};
