//! `-i/--index`: rank Astrometry.net index files and run the blind solver over them
//! to estimate the image position for the catalogue solve.

use alloc::sync::Arc;
use std::fs;
use std::path::{Path, PathBuf};

use arcsec_core::ArcsecError;
use arcsec_core::catalog::{load_anet_index, peek_anet_scale};
use arcsec_core::pipeline::{BlindSolveParams, blind_solve};
use arcsec_core::types::ImageBuffer;

/// Index files tried per solve. They run on separate threads, so total blind time is
/// `max(t_index0, t_index1)` rather than the sum.
const BLIND_MAX_INDEXES: usize = 2;

/// What the blind stage concluded.
pub enum BlindOutcome {
    /// Position estimate (RA, Dec) in radians, from the best-scoring index.
    Found(f64, f64),
    /// No index produced a position.
    NotFound,
    /// Nothing solved, and at least one index reported too few stars.
    InsufficientStars { found: usize, required: usize },
}

/// Is `name` an Astrometry.net index file name (`index-*.fits`)?
fn is_index_name(name: &str) -> bool {
    name.starts_with("index-") && name.ends_with(".fits")
}

/// Collect and rank astrometry.net index files for blind solving.
///
/// For a single file, returns `[path]`. For a directory, peeks every
/// `index-*.fits` header, filters out files whose scale range is incompatible
/// with `fov_deg`, and returns the remainder sorted by closeness of scale
/// midpoint to `fov_deg / 2` (best match first).
///
/// When `fov_deg <= 0` the scale filter is skipped and files are returned
/// in ascending filename order.
pub fn collect_index_files(path: &Path, fov_deg: f64) -> Vec<PathBuf> {
    if path.is_file() {
        return vec![path.to_path_buf()];
    }
    let Ok(rd) = fs::read_dir(path) else {
        return vec![];
    };
    let mut raw: Vec<PathBuf> = rd
        .filter_map(Result::ok)
        .map(|e| e.path())
        .filter(|p| {
            p.is_file()
                && p.file_name()
                    .and_then(|n| n.to_str())
                    .is_some_and(is_index_name)
        })
        .collect();
    raw.sort();

    if fov_deg <= 0.0 {
        return raw;
    }

    // Peek each file's scale range; filter + rank by match to FOV.
    // A "compatible" index has some overlap with [fov*0.2, fov*1.5].
    let fov_lo = fov_deg * 0.2;
    let fov_hi = fov_deg * 1.5;
    let ideal = fov_deg / 2.0;

    let mut ranked: Vec<(PathBuf, f64)> = raw
        .into_iter()
        .filter_map(|p| {
            // An unreadable file is skipped silently.
            let (_dq, lo_rad, hi_rad) = peek_anet_scale(&p).ok()?;
            let lo_deg = lo_rad.to_degrees();
            let hi_deg = hi_rad.to_degrees();
            if hi_deg < fov_lo || lo_deg > fov_hi {
                log::info!(
                    "Blind: skipping {} (scale {lo_deg:.2}°–{hi_deg:.2}°, outside [{fov_lo:.2}°–{fov_hi:.2}°])",
                    p.display(),
                );
                return None;
            }
            let mid = (lo_deg + hi_deg) / 2.0;
            Some((p, (mid - ideal).abs()))
        })
        .collect();

    // Ascending distance to the ideal scale: best first.
    ranked.sort_by(|a, b| a.1.total_cmp(&b.1));
    ranked.into_iter().map(|(p, _)| p).collect()
}

/// Run the blind solver over the best-ranked `index_files` in parallel and keep the
/// highest-scoring position.
///
/// The number of concurrent indexes is capped by the thread limit, so
/// `--threads 1` stays genuinely single-threaded.
pub fn estimate_position(
    img: &ImageBuffer,
    index_files: &[PathBuf],
    params: &BlindSolveParams,
) -> BlindOutcome {
    let max_indexes = BLIND_MAX_INDEXES.min(arcsec_core::max_threads().max(1));
    let img = Arc::new(img.clone());

    let handles: Vec<_> = index_files
        .iter()
        .take(max_indexes)
        .map(|idx_path| {
            let idx_path = idx_path.clone();
            let img = Arc::clone(&img);
            let params = params.clone();
            std::thread::spawn(move || -> Result<(f64, f64, usize), ArcsecError> {
                log::info!("Blind: trying index {}", idx_path.display());
                let anet_index = load_anet_index(&idx_path)?;
                let res = blind_solve(&img, &anet_index, &params);
                if let Ok((ra, dec, score)) = &res {
                    log::info!(
                        "Blind: {} → RA={:.3}° Dec={:.3}° score={score}",
                        idx_path
                            .file_name()
                            .map(|n| n.to_string_lossy().into_owned())
                            .unwrap_or_default(),
                        ra.to_degrees(),
                        dec.to_degrees(),
                    );
                }
                res
            })
        })
        .collect();

    let mut best: Option<(f64, f64, usize)> = None;
    let mut too_few: Option<(usize, usize)> = None;
    for handle in handles {
        match handle.join() {
            Ok(Ok((ra, dec, score))) => {
                if best.is_none_or(|(_, _, s)| score > s) {
                    best = Some((ra, dec, score));
                }
            }
            Ok(Err(ArcsecError::InsufficientStars { found, required })) => {
                too_few = Some((found, required));
            }
            Ok(Err(e)) => log::info!("Blind: did not solve: {e}"),
            Err(_) => log::info!("Blind: thread panicked"),
        }
    }

    // One index running short of stars must not discard another's solution.
    match (best, too_few) {
        (Some((ra, dec, _)), _) => BlindOutcome::Found(ra, dec),
        (None, Some((found, required))) => BlindOutcome::InsufficientStars { found, required },
        (None, None) => BlindOutcome::NotFound,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn index_names_are_recognised() {
        assert!(is_index_name("index-4107.fits"));
        assert!(is_index_name("index-5200-07.fits"));
        assert!(!is_index_name("index-4107.fits.part"));
        assert!(!is_index_name("d50_0101.1476"));
    }

    #[test]
    fn a_missing_path_yields_no_indexes() {
        assert!(collect_index_files(Path::new("/nonexistent/arcsec/indexes"), 1.0).is_empty());
    }
}
