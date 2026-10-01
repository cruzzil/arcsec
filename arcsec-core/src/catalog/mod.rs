//! Star catalogue readers: ASTAP `.1476`/`.290`/`.001` databases and
//! Astrometry.net index files, plus the sky tilings used to find the right files.

pub mod anet;
pub mod areas;
pub mod areas_290;
pub mod format_001;
pub mod format_1476;

pub use anet::{AnetIndex, AnetIndexEntry, AnetStar, load_anet_index, peek_anet_scale};
pub use areas::{DEC_BOUNDARIES_1476, area_and_boundaries_1476, filename_1476, find_areas_1476};
pub use areas_290::{DEC_BOUNDARIES_290, area_nr_290, filename_290, find_areas_290};
pub use format_001::read_001_file;
pub use format_1476::{
    CatalogLayout, CatalogStar, catalog_present, detect_layout, for_each_star_in_dec_band,
    read_area_file, read_catalog_stars,
};
