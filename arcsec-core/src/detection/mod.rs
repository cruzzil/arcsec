//! Star detection: background/noise estimation, the multi-pass star finder, and
//! the solve-free image analysis behind `--analyse` and `--extract`.

pub mod analyse;
pub mod background;
pub mod stars;

pub use analyse::{Analysis, MeasuredStar, analyse_image};
pub use background::{Background, get_background};
pub use stars::{find_stars, find_stars_with_background};
