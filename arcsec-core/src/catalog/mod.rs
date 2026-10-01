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
    CatalogLayout, CatalogStar, catalog_present, detect_layout, read_area_file, read_catalog_stars,
};

/// The star density, in stars per square degree, that the database `db_name` is
/// built to, or `None` if it is not known.
///
/// ASTAP's databases are density-limited: every tile holds the brightest Gaia stars
/// up to a fixed number per square degree, and ASTAP reads that number from the
/// name, two digits times 100 (`d80` 8000, `d50` 5000, `d20` 2000, `d05`, `v05` and
/// `g05` 500, `v50` 5000). The exceptions are the old 17/18 databases, which were
/// magnitude-limited (no figure), and `w08`, whose digits are its magnitude limit,
/// 8: its single all-sky file holds 41 265 stars, one per square degree.
///
/// A field can hold no more catalogue stars than this density times its area, so
/// detecting more image stars than that only adds patterns that cannot match.
#[must_use]
pub fn database_density(db_name: &str) -> Option<f64> {
    let name = db_name.to_ascii_lowercase();
    if name == "w08" {
        return Some(1.0);
    }
    let n: u32 = name.get(1..3)?.parse().ok()?;
    match n {
        0 | 17 | 18 => None,
        n => Some(f64::from(n) * 100.0),
    }
}

#[cfg(test)]
mod density_tests {
    use super::database_density;

    #[test]
    fn densities_follow_the_database_names() {
        for (name, want) in [
            ("d80", Some(8000.0)),
            ("D50", Some(5000.0)),
            ("d20", Some(2000.0)),
            ("d05", Some(500.0)),
            ("v50", Some(5000.0)),
            ("v05", Some(500.0)),
            ("g05", Some(500.0)),
            ("w08", Some(1.0)),
            ("t50", Some(5000.0)),
            ("v17", None),
            ("g18", None),
            ("x", None),
            ("", None),
            ("dxx", None),
        ] {
            assert_eq!(database_density(name), want, "{name}");
        }
    }
}
