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

// ── arcsec's own index ──────────────────────────────────────────────────────────

/// Relative uncertainty allowed on a pixel scale taken from `--fov` or the header.
/// FOCALLEN and XPIXSZ are often a few percent off (reducers, binning written
/// inconsistently); 20% covers that without letting the vote spread.
const SCALE_SLACK: f64 = 1.2;

/// Pixel scales searched when nothing gives one, arcseconds per pixel.
const SCALE_UNKNOWN: (f64, f64) = (0.3, 60.0);

/// The arcsec blind index `path` names: the file itself, or the first `*.arcsecix`
/// in a directory. `None` when there is none (the path may still hold
/// Astrometry.net files).
pub fn find_arcsec_index(path: &Path) -> Option<PathBuf> {
    if path.is_file() {
        return arcsec_core::index::is_blind_index(path).then(|| path.to_path_buf());
    }
    crate::catalog_cmd::index_cmd::index_files(path)
        .into_iter()
        .find(|p| arcsec_core::index::is_blind_index(p))
}

/// Search radius, in fields, that the spiral covers before an automatically found
/// index is consulted. Inside it the result is exactly the spiral's, so a usable
/// hint solves as it always did; beyond it the spiral's cost grows with the square
/// of the radius, while the index's does not grow at all.
pub const AUTO_SPIRAL_FIELDS: f64 = 5.0;

/// The smallest stage-one spiral radius, radians (1°).
const AUTO_SPIRAL_MIN: f64 = 1.0 * core::f64::consts::PI / 180.0;

/// An arcsec index to use, and whether the user named it.
pub struct OwnIndex {
    path: PathBuf,
    explicit: bool,
}

/// The arcsec index for this solve: the one `--index` names, if it names one;
/// otherwise, when the search radius reaches past [`AUTO_SPIRAL_FIELDS`] fields,
/// one installed in the catalogue directory or beside the star database.
pub fn arcsec_index_for(
    explicit: Option<&PathBuf>,
    template: &arcsec_core::pipeline::SolveParams,
) -> Option<OwnIndex> {
    if let Some(p) = explicit {
        return find_arcsec_index(p).map(|path| OwnIndex {
            path,
            explicit: true,
        });
    }
    if template.search_radius <= stage_one_radius(template) {
        return None;
    }
    find_arcsec_index(&crate::catalog_cmd::default_dir())
        .or_else(|| find_arcsec_index(&template.db_path))
        .map(|path| OwnIndex {
            path,
            explicit: false,
        })
}

fn stage_one_radius(template: &arcsec_core::pipeline::SolveParams) -> f64 {
    (template.fov * AUTO_SPIRAL_FIELDS).max(AUTO_SPIRAL_MIN)
}

/// Solve with an arcsec index; `None` if nothing verified, and the caller runs
/// the ordinary search.
///
/// Named with `--index`, the index is tried first. Found automatically, the
/// spiral first searches [`AUTO_SPIRAL_FIELDS`] fields round the hint (if there is
/// a hint), which returns exactly what the full search would for any field that
/// close; only then is the index consulted, limited to `-r` round the hint unless
/// the radius is the whole sky.
///
/// `scale` is arcseconds per pixel of the image as solved (after binning);
/// `scale_known` says whether it came from the user or the header rather than the
/// 1″/px fallback.
pub fn index_stage(
    img: &ImageBuffer,
    ix: &OwnIndex,
    template: &arcsec_core::pipeline::SolveParams,
    has_hint: bool,
    scale: f64,
    scale_known: bool,
) -> Option<arcsec_core::types::WcsSolution> {
    use core::f64::consts::PI;
    if !ix.explicit && has_hint {
        let r0 = stage_one_radius(template);
        log::info!(
            "Searching {:.1}° round the hint before the blind index.",
            r0.to_degrees()
        );
        let near = arcsec_core::pipeline::SolveParams {
            search_radius: r0,
            ..template.clone()
        };
        if let Ok(w) = arcsec_core::pipeline::solve_image(img, &near) {
            return Some(w);
        }
    }
    // Named with --index the solve is blind, as with Astrometry.net files; found
    // automatically it stands in for the rest of the spiral, so it keeps to -r.
    let within = (!ix.explicit && has_hint && template.search_radius < PI).then_some((
        template.ra_hint,
        template.dec_hint,
        template.search_radius + template.fov,
    ));
    let mut wcs = solve_with_arcsec_index(img, &ix.path, template, scale, scale_known, within)?;
    // The hinted solve started at the index's hypothesis; report the distance from
    // the user's start position, as the spiral would.
    let (s1, c1) = template.dec_hint.sin_cos();
    let (s2, c2) = wcs.dec0.sin_cos();
    let cos_d = (s1 * s2 + c1 * c2 * (wcs.ra0 - template.ra_hint).cos()).clamp(-1.0, 1.0);
    wcs.search_dist_deg = cos_d.acos().to_degrees();
    Some(wcs)
}

/// Solve with an arcsec blind index; `None` if it found nothing that verified.
fn solve_with_arcsec_index(
    img: &ImageBuffer,
    path: &Path,
    template: &arcsec_core::pipeline::SolveParams,
    scale: f64,
    scale_known: bool,
    within: Option<(f64, f64, f64)>,
) -> Option<arcsec_core::types::WcsSolution> {
    let t0 = std::time::Instant::now();
    let index = match arcsec_core::index::BlindIndex::open(path) {
        Ok(ix) => ix,
        Err(e) => {
            eprintln!("Blind index {}: {e}", path.display());
            return None;
        }
    };
    let (scale_lo, scale_hi) = if scale_known {
        (scale / SCALE_SLACK, scale * SCALE_SLACK)
    } else {
        SCALE_UNKNOWN
    };
    log::info!(
        "Blind index {} ({} patterns), pixel scale {scale_lo:.3}–{scale_hi:.3}\"/px",
        path.display(),
        index.n_patterns()
    );
    let params = arcsec_core::pipeline::IndexSolveParams {
        scale_lo,
        scale_hi,
        within,
    };
    match arcsec_core::pipeline::index_solve(img, &index, template, &params) {
        Ok((wcs, stats)) => {
            log::info!(
                "Blind index: solved in {:.2} s (hypothesis rank {:?}, score {}, {} hinted solves)",
                t0.elapsed().as_secs_f64(),
                stats.accepted_rank,
                stats.best_score,
                stats.verified
            );
            Some(wcs)
        }
        Err(e) => {
            log::info!(
                "Blind index: no solution after {:.2} s: {e}",
                t0.elapsed().as_secs_f64()
            );
            None
        }
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
