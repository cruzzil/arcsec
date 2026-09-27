//! Choosing a star database: which directory (`-d`) and which database (`-D`)
//! when the user does not say.

use std::fs;
use std::path::{Path, PathBuf};

use crate::catalog_cmd;

/// The star database directory to use without `-d`.
///
/// The directory `arcsec catalog install` writes to, so a user who installed a
/// catalogue never has to say where it went — provided a star database is actually
/// there. Otherwise the working directory, which is ASTAP's behaviour. Blind
/// indexes do not count: a managed directory holding only `anet-4100` would
/// otherwise hide the ASTAP databases the user keeps in the working directory.
pub fn default_db_path() -> PathBuf {
    let managed = catalog_cmd::default_dir();
    if catalog_cmd::REGISTRY
        .iter()
        .filter(|e| matches!(e.files, catalog_cmd::Files::AstapDb { .. }))
        .any(|e| catalog_cmd::is_installed(&managed, e))
    {
        managed
    } else {
        PathBuf::from(".")
    }
}

/// Field-of-view range each ASTAP database is built for, and how much we prefer it
/// when several are eligible (higher = denser, so a better fit).
///
/// Ranges are as published on the ASTAP download page. The D-series stops at 6°; G05
/// and W08 exist precisely to cover wider fields, and before `.290`/`.001` support
/// they could not be read at all, which is why fields beyond ~6° never solved.
const DB_FOV_RANGES: &[(&str, f64, f64, u8)] = &[
    // prefix, min FOV (deg), max FOV (deg), preference
    ("d80", 0.15, 6.0, 8),
    ("v50", 0.20, 6.0, 6),
    ("d50", 0.20, 6.0, 5),
    ("d20", 0.30, 6.0, 4),
    ("v05", 0.60, 6.0, 3),
    ("d05", 0.60, 6.0, 2),
    ("g05", 3.00, 20.0, 7),
    ("w08", 20.0, 80.0, 7),
];

/// Database prefix of a star database file name, if it is one.
///
/// Every supported format is `<prefix>_<cell>.<ext>`; see
/// [`catalog_cmd::ASTAP_EXTS`].
fn db_prefix(file_name: &str) -> Option<&str> {
    if catalog_cmd::ASTAP_EXTS
        .iter()
        .any(|(ext, _)| file_name.ends_with(ext))
    {
        file_name.split('_').next()
    } else {
        None
    }
}

/// Every database prefix present in `db_path`, in any supported format.
fn available_dbs(db_path: &Path) -> Vec<String> {
    let mut out: Vec<String> = fs::read_dir(db_path)
        .into_iter()
        .flatten()
        .flatten()
        .filter_map(|e| db_prefix(&e.file_name().to_string_lossy()).map(str::to_string))
        .collect();
    out.sort();
    out.dedup();
    out
}

/// Pick the installed database best suited to a `fov_deg` field.
pub fn select_db_for_fov(db_path: &Path, fov_deg: f64) -> Option<String> {
    select_from(&available_dbs(db_path), fov_deg)
}

/// Pick from the databases in `present` for a `fov_deg` field.
///
/// Prefers the densest database whose published range contains the field, then any
/// database whose range merely comes closest — so an unusual field size still gets
/// the least-bad option rather than nothing.
fn select_from(present: &[String], fov_deg: f64) -> Option<String> {
    if present.is_empty() {
        return None;
    }
    let installed = DB_FOV_RANGES
        .iter()
        .filter(|(prefix, ..)| present.iter().any(|p| p == prefix));

    let mut best: Option<(u8, &str)> = None;
    for &(prefix, lo, hi, pref) in installed.clone() {
        if fov_deg >= lo && fov_deg <= hi && best.is_none_or(|(bp, _)| pref > bp) {
            best = Some((pref, prefix));
        }
    }
    if let Some((_, prefix)) = best {
        return Some(prefix.to_string());
    }

    // Nothing covers this field: take the database whose range is nearest.
    let mut fallback: Option<(f64, &str)> = None;
    for &(prefix, lo, hi, _) in installed {
        let dist = if fov_deg < lo {
            lo - fov_deg
        } else {
            fov_deg - hi
        };
        if fallback.is_none_or(|(bd, _)| dist < bd) {
            fallback = Some((dist, prefix));
        }
    }
    fallback
        .map(|(_, p)| p.to_string())
        .or_else(|| present.first().cloned())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn owned(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| (*s).to_string()).collect()
    }

    #[test]
    fn prefixes_come_from_every_database_format() {
        assert_eq!(db_prefix("d50_0101.1476"), Some("d50"));
        assert_eq!(db_prefix("g05_0101.290"), Some("g05"));
        assert_eq!(db_prefix("w08_0101.001"), Some("w08"));
        assert_eq!(db_prefix("index-4107.fits"), None);
        assert_eq!(db_prefix("readme.txt"), None);
    }

    #[test]
    fn the_densest_covering_database_wins() {
        let all = owned(&["d05", "d50", "d80", "g05", "w08"]);
        assert_eq!(select_from(&all, 1.0).as_deref(), Some("d80"));
        assert_eq!(select_from(&all, 10.0).as_deref(), Some("g05"));
        assert_eq!(select_from(&all, 40.0).as_deref(), Some("w08"));
        // 3°–6° is covered by both the D-series and G05; D80 is denser.
        assert_eq!(select_from(&all, 4.0).as_deref(), Some("d80"));
    }

    #[test]
    fn an_uncovered_field_gets_the_nearest_range() {
        let d = owned(&["d50", "g05"]);
        assert_eq!(select_from(&d, 0.05).as_deref(), Some("d50"));
        assert_eq!(select_from(&d, 30.0).as_deref(), Some("g05"));
    }

    #[test]
    fn unknown_prefixes_are_a_last_resort() {
        assert_eq!(select_from(&owned(&["h18"]), 1.0).as_deref(), Some("h18"));
        assert_eq!(select_from(&[], 1.0), None);
    }
}
