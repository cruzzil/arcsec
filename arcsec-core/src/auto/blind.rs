//! Blind solving for [`super::Plan`]: rank Astrometry.net index files and run the
//! blind solver over them to estimate the image position for the catalogue solve,
//! and decide when and how arcsec's own index is used.

use alloc::sync::Arc;
use std::fs;
use std::path::{Path, PathBuf};

use crate::ArcsecError;
use crate::catalog::{load_anet_index, peek_anet_scale};
use crate::pipeline::{BlindSolveParams, blind_solve};
use crate::types::ImageBuffer;

/// Index files tried per solve. They run on separate threads, so total blind time is
/// `max(t_index0, t_index1)` rather than the sum.
const BLIND_MAX_INDEXES: usize = 2;

/// What the blind stage concluded.
pub(crate) enum BlindOutcome {
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
#[must_use]
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
/// `--threads 1` stays genuinely single-threaded. Each index thread inherits the
/// caller's thread limit and cancellation token.
pub(crate) fn estimate_position(
    img: &ImageBuffer,
    index_files: &[PathBuf],
    params: &BlindSolveParams,
) -> BlindOutcome {
    let max_indexes = BLIND_MAX_INDEXES.min(crate::max_threads().max(1));
    let img = Arc::new(img.clone());
    let cancel = crate::cancel::current();
    let local_threads = crate::local_max_threads();

    let handles: Vec<_> = index_files
        .iter()
        .take(max_indexes)
        .map(|idx_path| {
            let idx_path = idx_path.clone();
            let img = Arc::clone(&img);
            let params = params.clone();
            let cancel = cancel.clone();
            std::thread::spawn(move || -> Result<(f64, f64, usize), ArcsecError> {
                log::info!("Blind: trying index {}", idx_path.display());
                let anet_index = load_anet_index(&idx_path)?;
                let res = crate::cancel::with_optional(cancel.as_ref(), || {
                    crate::with_max_threads(local_threads, || {
                        blind_solve(&img, &anet_index, &params)
                    })
                });
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
#[must_use]
pub fn find_arcsec_index(path: &Path) -> Option<PathBuf> {
    if path.is_file() {
        return crate::index::is_blind_index(path).then(|| path.to_path_buf());
    }
    super::index_files(path)
        .into_iter()
        .find(|p| crate::index::is_blind_index(p))
}

/// Search radius, in fields, that the spiral covers before an automatically found
/// index is consulted. Inside it the result is exactly the spiral's, so a usable
/// hint solves as it always did; beyond it the spiral's cost grows with the square
/// of the radius, while the index's does not grow at all.
pub(crate) const AUTO_SPIRAL_FIELDS: f64 = 5.0;

/// The smallest stage-one spiral radius, radians (1°).
const AUTO_SPIRAL_MIN: f64 = 1.0 * core::f64::consts::PI / 180.0;

/// Smallest `-r` at which an installed index is consulted automatically, radians
/// (10°). Below it a failed search is cheap anyway (about a second on the corpus),
/// and consulting the index roughly doubled it for nothing; above it the spiral's
/// cost dominates and the index adds a few percent to a failure.
const AUTO_MIN_RADIUS: f64 = 10.0 * core::f64::consts::PI / 180.0;

/// An arcsec index to use, and whether the user named it.
pub(crate) struct OwnIndex {
    path: PathBuf,
    explicit: bool,
}

/// The arcsec index for this solve: the one `--index` names, if it names one;
/// otherwise, when the search radius reaches past [`AUTO_SPIRAL_FIELDS`] fields and
/// is at least 10°, one installed in the catalogue directory or beside the star
/// database.
pub(crate) fn arcsec_index_for(
    explicit: Option<&PathBuf>,
    template: &crate::pipeline::SolveParams,
) -> Option<OwnIndex> {
    if let Some(p) = explicit {
        return find_arcsec_index(p).map(|path| OwnIndex {
            path,
            explicit: true,
        });
    }
    if template.search_radius <= stage_one_radius(template)
        || template.search_radius < AUTO_MIN_RADIUS
    {
        return None;
    }
    find_arcsec_index(&super::default_catalog_dir())
        .or_else(|| find_arcsec_index(&template.db_path))
        .map(|path| OwnIndex {
            path,
            explicit: false,
        })
}

fn stage_one_radius(template: &crate::pipeline::SolveParams) -> f64 {
    (template.fov * AUTO_SPIRAL_FIELDS).max(AUTO_SPIRAL_MIN)
}

/// What [`index_stage`] concluded.
pub(crate) enum IndexOutcome {
    /// A verified solution within the search radius (or anywhere, for a named
    /// index).
    Solved(Box<crate::types::WcsSolution>),
    /// The index verified the field outside the search radius, at this distance
    /// from the hint (degrees): the ordinary search cannot find it within `-r`,
    /// and there is no need to run it.
    Elsewhere(f64),
    /// Nothing verified; the caller runs the ordinary search.
    NotFound,
    /// The solve was cancelled.
    Cancelled,
}

/// Solve with an arcsec index.
///
/// Named with `--index`, the index is tried first. Found automatically, the
/// spiral first searches [`AUTO_SPIRAL_FIELDS`] fields round the hint (if there is
/// a hint), which returns exactly what the full search would for any field that
/// close; only then is the index consulted, limited to `-r` round the hint unless
/// the radius is the whole sky.
///
/// When that finds nothing, the index is asked again without the limit. A field it
/// verifies more than [`ELSEWHERE_FIELDS`] fields beyond `-r` is somewhere the rest
/// of the spiral cannot reach, and the spiral would spend all of its time (the
/// bulk of a failed search: 6000 positions for a 0.2° field at `-r 10`) finding
/// nothing; [`IndexOutcome::Elsewhere`] tells the caller to stop. The answer is
/// still "no solution", as `-r` requires: the field is not reported.
///
/// `scale` is arcseconds per pixel of the image as solved (after binning);
/// `scale_known` says whether it came from the user or the header rather than the
/// 1″/px fallback.
pub(crate) fn index_stage(
    img: &ImageBuffer,
    ix: &OwnIndex,
    template: &crate::pipeline::SolveParams,
    has_hint: bool,
    scale: f64,
    scale_known: bool,
) -> IndexOutcome {
    use core::f64::consts::PI;
    if !ix.explicit && has_hint {
        let r0 = stage_one_radius(template);
        log::info!(
            "Searching {:.1}° round the hint before the blind index.",
            r0.to_degrees()
        );
        let near = crate::pipeline::SolveParams {
            search_radius: r0,
            ..template.clone()
        };
        match crate::pipeline::solve_image(img, &near) {
            Ok(w) => return IndexOutcome::Solved(Box::new(w)),
            Err(ArcsecError::Cancelled) => return IndexOutcome::Cancelled,
            Err(_) => {}
        }
    }
    // Named with --index the solve is blind, as with Astrometry.net files; found
    // automatically it stands in for the rest of the spiral, so it keeps to -r.
    let within = (!ix.explicit && has_hint && template.search_radius < PI).then_some((
        template.ra_hint,
        template.dec_hint,
        template.search_radius + template.fov,
    ));
    let from_hint = |w: &crate::types::WcsSolution| {
        crate::math::coords::ang_sep(w.ra0, w.dec0, template.ra_hint, template.dec_hint)
    };
    if let Some(mut wcs) =
        solve_with_arcsec_index(img, &ix.path, template, scale, scale_known, within)
    {
        // The hinted solve started at the index's hypothesis; report the distance
        // from the user's start position, as the spiral would.
        wcs.search_dist_deg = from_hint(&wcs).to_degrees();
        return IndexOutcome::Solved(Box::new(wcs));
    }
    if within.is_some()
        && let Some(wcs) =
            solve_with_arcsec_index(img, &ix.path, template, scale, scale_known, None)
    {
        let sep = from_hint(&wcs);
        if beyond_reach(sep, template) {
            log::info!(
                "Blind index: the field is at RA={:.4}° Dec={:.4}°, {:.1}° from the hint, \
                 outside the search radius; not searching it.",
                wcs.ra0.to_degrees(),
                wcs.dec0.to_degrees(),
                sep.to_degrees()
            );
            return IndexOutcome::Elsewhere(sep.to_degrees());
        }
    }
    if crate::cancel::is_cancelled() {
        return IndexOutcome::Cancelled;
    }
    IndexOutcome::NotFound
}

/// How many fields past `-r` a field the index verifies must lie for the
/// ordinary search to be skipped. The spiral reads catalogue windows up to half a
/// step and up to a field beyond `-r`, so a field just outside the radius could
/// still be matched there; two fields is clear of that.
const ELSEWHERE_FIELDS: f64 = 2.0;

/// Whether a field centred `sep` radians from the hint is out of the ordinary
/// search's reach.
fn beyond_reach(sep: f64, template: &crate::pipeline::SolveParams) -> bool {
    sep > template.search_radius + ELSEWHERE_FIELDS * template.fov
}

/// Solve with an arcsec blind index; `None` if it found nothing that verified.
fn solve_with_arcsec_index(
    img: &ImageBuffer,
    path: &Path,
    template: &crate::pipeline::SolveParams,
    scale: f64,
    scale_known: bool,
    within: Option<(f64, f64, f64)>,
) -> Option<crate::types::WcsSolution> {
    let t0 = std::time::Instant::now();
    let index = match crate::index::BlindIndex::open(path) {
        Ok(ix) => ix,
        Err(e) => {
            log::warn!("Blind index {}: {e}", path.display());
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
    let params = crate::pipeline::IndexSolveParams {
        scale_lo,
        scale_hi,
        within,
    };
    match crate::pipeline::index_solve(img, &index, template, &params) {
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
    fn only_a_field_two_fields_past_the_radius_is_beyond_reach() {
        use crate::pipeline::{SearchSpeed, SolveMethod, SolveParams};
        let deg = f64::to_radians;
        let t = SolveParams {
            ra_hint: 0.0,
            dec_hint: 0.0,
            fov: deg(0.5),
            search_radius: deg(10.0),
            quad_tolerance: 0.007,
            hfd_min: 1.5,
            max_stars: 500,
            db_path: PathBuf::new(),
            db_name: "d80".into(),
            binning: 1,
            method: SolveMethod::Quads,
            speed: SearchSpeed::Auto,
            threads: 1,
        };
        assert!(!beyond_reach(deg(5.0), &t));
        assert!(!beyond_reach(deg(10.9), &t));
        assert!(beyond_reach(deg(11.1), &t));
        assert!(beyond_reach(deg(40.0), &t));
    }

    #[test]
    fn a_missing_path_yields_no_indexes() {
        assert!(collect_index_files(Path::new("/nonexistent/arcsec/indexes"), 1.0).is_empty());
    }
}
