//! Star detection: background/noise estimation and the multi-pass star finder.

pub mod background;
pub mod stars;

pub use background::{Background, get_background};
pub use stars::{find_stars, find_stars_with_background};
