// Full plate-solving pipeline.

use core::f64::consts::PI;
use std::io;
use std::path::PathBuf;

use crate::catalog::read_catalog_stars;
use crate::detection::get_background;
use crate::detection::stars::find_stars_with_background;
use crate::error::{ArcsecError, Result};
use crate::math::coords::{ang_sep, equatorial_standard};
use crate::math::lsq::solve_plate_constants;
use crate::quads::{
    TETRA_TOL_FACTOR, bijective_filter, build_quads, build_quads_presorted, build_triangles,
    extract_star_pairs, extract_triangle_pairs, filter_by_scale, filter_triangles_by_scale,
    find_matches_sorted, find_triangle_matches, vote_filter,
};
use crate::types::{PairedPositions, PlateConstants, Star, StarList, WcsSolution};
use crate::wcs::output::derive_wcs;

use super::spiral::SpiralSearch;

/// Which pattern-matching algorithm to use in the catalog spiral loop.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum SolveMethod {
    /// ASTAP-style 5-ratio quad matching with vote_filter (default).
    #[default]
    Quads,
    /// TETRA 2-ratio triangle matching with bijective filter.
    Tetra,
}

/// Parameters for `solve_image`.
pub struct SolveParams {
    /// Approximate RA of image centre (radians, hint only).
    pub ra_hint: f64,
    /// Approximate DEC of image centre (radians, hint only).
    pub dec_hint: f64,
    /// Image field of view (square side, radians). Used as the spiral step size.
    pub fov: f64,
    /// Maximum search radius from the hint position (radians).
    pub search_radius: f64,
    /// Quad ratio matching tolerance (ASTAP default ≈ 0.007).
    pub quad_tolerance: f64,
    /// Minimum HFD for valid stars (pixels).
    pub hfd_min: f64,
    /// Maximum number of image stars to detect.
    pub max_stars: usize,
    /// Path to the catalog directory.
    pub db_path: PathBuf,
    /// Catalog name prefix (e.g. `"d20"`, `"d80"`).
    pub db_name: String,
    /// Pixel binning factor applied before solving (1 = none, 2 = 2×2, ...).
    /// WCS output is scaled back to original image pixel coordinates.
    pub binning: usize,
    /// Pattern-matching algorithm for the catalog spiral.
    pub method: SolveMethod,
    /// Worker threads for the spiral search. 0 = one per available core.
    ///
    /// Spiral positions are independent, so they are evaluated a batch at a time
    /// across this many threads. Results are identical to the serial search: within
    /// a batch the lowest spiral index still wins, so the first position that
    /// verifies is the one returned, exactly as before.
    pub threads: usize,
}

/// Iterative sigma-clipping: fit plate constants, reject pairs with large residuals,
/// re-fit until stable or fewer than `min_count` pairs remain.
///
/// First pass uses a 10-pixel absolute threshold (in catalog arcsec) to cut the
/// large residuals of false-positive triangle matches. Subsequent passes apply
/// `sigma × rms` clipping until the set is stable.
fn sigma_clip_pairs(
    mut img_pos: Vec<(f64, f64)>,
    mut cat_pos: Vec<(f64, f64)>,
    sigma: f64,
    min_count: usize,
) -> PairedPositions {
    let mut first_pass = true;
    for _ in 0..10 {
        if img_pos.len() < min_count.max(3) {
            break;
        }
        let plate = match solve_plate_constants(&img_pos, &cat_pos) {
            Ok(p) => p,
            Err(_) => break,
        };
        let residuals: Vec<f64> = img_pos
            .iter()
            .zip(cat_pos.iter())
            .map(|(&(xi, yi), &(xc, yc))| {
                let xp = plate.a * xi + plate.b * yi + plate.c;
                let yp = plate.d * xi + plate.e * yi + plate.f;
                ((xp - xc).powi(2) + (yp - yc).powi(2)).sqrt()
            })
            .collect();
        let rms = (residuals.iter().map(|r| r * r).sum::<f64>() / residuals.len() as f64).sqrt();
        let threshold = if first_pass {
            first_pass = false;
            // 10 px in catalog-arcsec: generous cut for large FP residuals on first pass.
            let cdelt = (plate.a.powi(2) + plate.d.powi(2)).sqrt();
            (10.0 * cdelt).max(10.0)
        } else {
            sigma * rms
        };
        let before = img_pos.len();
        let mut new_img = Vec::with_capacity(before);
        let mut new_cat = Vec::with_capacity(before);
        for ((&ip, &cp), &r) in img_pos.iter().zip(cat_pos.iter()).zip(residuals.iter()) {
            if r <= threshold {
                new_img.push(ip);
                new_cat.push(cp);
            }
        }
        if new_img.len() == before {
            break; // stable — no more outliers
        }
        img_pos = new_img;
        cat_pos = new_cat;
    }
    (img_pos, cat_pos)
}

/// Minimum number of individually matched stars required to believe a solution.
///
/// Correct solves typically match 200-375 stars, so this is deliberately loose;
/// its job is to reject the handful-of-coincidences case. Together with
/// MIN_VERIFY_SPREAD it separates two otherwise identical-looking results: M31 at
/// 2 degrees (22 stars, spread 0.221, rms 0.65", rotation wrong by 1.56 degrees)
/// from the Dec -88 field (46 stars, spread 0.207, rms 0.66", correct to 2.3").
const MIN_VERIFIED_STARS: usize = 30;
/// Match radii (pixels) used by successive verification passes, coarse to fine.
const VERIFY_RADII: [f64; 3] = [6.0, 3.0, 2.0];
/// Minimum spread of the matched stars, as a fraction of the image half-diagonal.
///
/// A count threshold alone is not enough: matches clustered in one part of the
/// frame (the core of a bright galaxy, say) pin the position but leave rotation
/// and scale essentially free. M31 at 2 degrees passed with 22 matched stars and
/// a 1.56-degree rotation error, which is 154" at the field corners.
const MIN_VERIFY_SPREAD: f64 = 0.20;

/// One verification pass: the re-fitted plate, the number of matched stars, the
/// per-star RMS in arcsec, and the spatial spread of the matches.
type VerifyPass = (PlateConstants, usize, f64, f64);

/// Project the catalogue onto the image with a candidate plate solution, match
/// individual stars, and re-fit on those matches.
///
/// The quad matcher only ever produces quad *centroids*, so the plate fit is built
/// from a handful of averaged positions and nothing ever checks that the individual
/// stars agree. This does that check: invert the plate to map every catalogue star
/// into pixel space, pair each with the nearest detected star, re-fit on the pairs,
/// and repeat with a shrinking radius.
///
/// Returns `(refined_plate, n_matched_stars, rms_arcsec)`, or `None` if the plate is
/// degenerate or too few stars agree.
fn verify_and_refit(
    img_stars: &StarList,
    cat_stars: &StarList,
    plate: &PlateConstants,
    img_w: usize,
    img_h: usize,
) -> Option<(PlateConstants, usize, f64)> {
    if img_stars.is_empty() || cat_stars.is_empty() {
        return None;
    }

    // Uniform grid over the detected stars for nearest-neighbour lookup.
    let (mut min_x, mut min_y) = (f64::INFINITY, f64::INFINITY);
    let (mut max_x, mut max_y) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    for st in &img_stars.0 {
        min_x = min_x.min(st.x);
        max_x = max_x.max(st.x);
        min_y = min_y.min(st.y);
        max_y = max_y.max(st.y);
    }
    if !(min_x.is_finite() && min_y.is_finite() && max_x > min_x && max_y > min_y) {
        return None;
    }
    let cell = VERIFY_RADII[0].max(1.0);
    let nx = (((max_x - min_x) / cell).ceil() as usize + 1).max(1);
    let ny = (((max_y - min_y) / cell).ceil() as usize + 1).max(1);
    let mut grid: Vec<Vec<u32>> = vec![Vec::new(); nx * ny];
    for (i, st) in img_stars.0.iter().enumerate() {
        let gx = ((st.x - min_x) / cell) as usize;
        let gy = ((st.y - min_y) / cell) as usize;
        grid[gy.min(ny - 1) * nx + gx.min(nx - 1)].push(i as u32);
    }

    let mut current = plate.clone();
    let mut best: Option<VerifyPass> = None;

    for &radius in VERIFY_RADII.iter() {
        let det = current.a * current.e - current.b * current.d;
        if det.abs() < 1e-12 {
            return None;
        }
        let r2 = radius * radius;

        let mut img_pos: Vec<(f64, f64)> = Vec::new();
        let mut cat_pos: Vec<(f64, f64)> = Vec::new();
        let mut used = vec![false; img_stars.len()];

        for cs in &cat_stars.0 {
            // Invert  xi = a*x + b*y + c ;  eta = d*x + e*y + f
            let dx = cs.x - current.c;
            let dy = cs.y - current.f;
            let px = (current.e * dx - current.b * dy) / det;
            let py = (-current.d * dx + current.a * dy) / det;
            if px < min_x - radius
                || px > max_x + radius
                || py < min_y - radius
                || py > max_y + radius
            {
                continue;
            }

            let gx = (((px - min_x) / cell) as isize).clamp(0, nx as isize - 1);
            let gy = (((py - min_y) / cell) as isize).clamp(0, ny as isize - 1);
            let mut best_i: Option<usize> = None;
            let mut best_d2 = r2;
            for oy in -1isize..=1 {
                for ox in -1isize..=1 {
                    let cx = gx + ox;
                    let cy = gy + oy;
                    if cx < 0 || cy < 0 || cx >= nx as isize || cy >= ny as isize {
                        continue;
                    }
                    for &i in &grid[cy as usize * nx + cx as usize] {
                        let i = i as usize;
                        if used[i] {
                            continue;
                        }
                        let st = &img_stars.0[i];
                        let d2 = (st.x - px) * (st.x - px) + (st.y - py) * (st.y - py);
                        if d2 < best_d2 {
                            best_d2 = d2;
                            best_i = Some(i);
                        }
                    }
                }
            }
            if let Some(i) = best_i {
                used[i] = true; // one-to-one: a detected star backs at most one catalogue star
                img_pos.push((img_stars.0[i].x, img_stars.0[i].y));
                cat_pos.push((cs.x, cs.y));
            }
        }

        if img_pos.len() < 4 {
            break;
        }
        let refined = match solve_plate_constants(&img_pos, &cat_pos) {
            Ok(p) => p,
            Err(_) => break,
        };
        let mut sq = 0.0;
        for (&(xi, yi), &(xc, yc)) in img_pos.iter().zip(cat_pos.iter()) {
            let xp = refined.a * xi + refined.b * yi + refined.c;
            let yp = refined.d * xi + refined.e * yi + refined.f;
            sq += (xp - xc).powi(2) + (yp - yc).powi(2);
        }
        let rms = (sq / img_pos.len() as f64).sqrt();
        // Spread of the matched stars about their own centroid, as a fraction of the
        // image half-diagonal. Matches clustered in one corner leave rotation free.
        let n = img_pos.len() as f64;
        let mx = img_pos.iter().map(|p| p.0).sum::<f64>() / n;
        let my = img_pos.iter().map(|p| p.1).sum::<f64>() / n;
        let var = img_pos
            .iter()
            .map(|&(x, y)| (x - mx) * (x - mx) + (y - my) * (y - my))
            .sum::<f64>()
            / n;
        let half_diag = 0.5 * ((img_w * img_w + img_h * img_h) as f64).sqrt();
        let spread = var.sqrt() / half_diag;
        log::debug!(
            "verify: {} stars, spread {:.3}, rms {:.2}\"",
            img_pos.len(),
            spread,
            rms
        );

        best = Some((refined.clone(), img_pos.len(), rms, spread));
        current = refined;
    }

    best.filter(|&(_, n, _, spread)| n >= MIN_VERIFIED_STARS && spread >= MIN_VERIFY_SPREAD)
        .map(|(p, n, r, _)| (p, n, r))
}

/// Everything a spiral position needs that does not change between positions.
struct SpiralCtx<'a> {
    params: &'a SolveParams,
    img: &'a crate::types::ImageBuffer,
    stars: &'a StarList,
    img_quads: &'a crate::types::QuadList,
    img_tris: &'a crate::quads::TriangleList,
    nrstars_image: usize,
    nrstars_required: usize,
    oversize: f64,
    min_quads: usize,
    step_size: f64,
}

/// A spiral position that produced a verified solution.
struct PositionOutcome {
    idx: usize,
    ra_db: f64,
    dec_db: f64,
    sep_deg: f64,
    plate: PlateConstants,
    n_verified: usize,
    rms: f64,
    n_matched: usize,
    n_raw: usize,
    mag_limit: f64,
}

/// Result of trying one spiral position: the angular distance if the catalogue was
/// actually read there (for the ASTAP-style progress line), and the solution if one
/// verified.
struct PositionTry {
    sep_deg: Option<f64>,
    outcome: Option<PositionOutcome>,
}

impl PositionTry {
    const NONE: Self = PositionTry {
        sep_deg: None,
        outcome: None,
    };
}

/// Evaluate a single spiral position. Pure with respect to `ctx`, so positions can
/// be run concurrently.
fn try_position(ctx: &SpiralCtx<'_>, idx: usize, sx: i32, sy: i32) -> PositionTry {
    let params = ctx.params;
    let step_size = ctx.step_size;

    let dec_db_raw = params.dec_hint + step_size * sy as f64;
    let (dec_db, flip) = if dec_db_raw > PI / 2.0 {
        (PI - dec_db_raw, PI)
    } else if dec_db_raw < -PI / 2.0 {
        (-PI - dec_db_raw, PI)
    } else {
        (dec_db_raw, 0.0)
    };

    let extra = if dec_db > 0.0 {
        step_size * 0.5
    } else {
        -step_size * 0.5
    };
    let ra_offset = step_size * sx as f64 / (dec_db - extra).cos();
    if ra_offset > PI / 2.0 + step_size * 0.5 || ra_offset < -PI / 2.0 {
        return PositionTry::NONE;
    }

    let ra_db = (flip + params.ra_hint + ra_offset).rem_euclid(2.0 * PI);
    let sep = ang_sep(ra_db, dec_db, params.ra_hint, params.dec_hint);
    if sep > params.search_radius + step_size / 2.0 {
        return PositionTry::NONE;
    }

    let cat_raw = match read_catalog_stars(
        &params.db_path,
        &params.db_name,
        ra_db,
        dec_db,
        params.fov * ctx.oversize,
        ctx.nrstars_required,
    ) {
        Ok(v) if !v.is_empty() => v,
        Ok(_) => return PositionTry::NONE,
        Err(ArcsecError::CatalogIo(ref e)) if e.kind() == io::ErrorKind::NotFound => {
            return PositionTry::NONE;
        }
        Err(_) => return PositionTry::NONE,
    };

    let sep_deg = sep.to_degrees();
    let mag_limit = cat_raw
        .iter()
        .map(|s| s.mag)
        .fold(f64::NEG_INFINITY, f64::max);
    log::info!(
        "Search {}, [{},{}], position: {}  Down to magn {:.1}  {} database stars  {} database quads to compare.",
        idx,
        sx,
        sy,
        format_radec(ra_db, dec_db),
        mag_limit,
        cat_raw.len(),
        cat_raw.len(),
    );

    let mut cat_stars: Vec<Star> = cat_raw
        .iter()
        .map(|s| {
            let (x, y) = equatorial_standard(ra_db, dec_db, s.ra, s.dec, 1.0);
            Star {
                x,
                y,
                snr: 1.0,
                hfd: 2.0,
            }
        })
        .collect();
    cat_stars.sort_unstable_by(|a, b| a.x.total_cmp(&b.x));
    let cat_star_list = StarList(cat_stars);

    let failed = PositionTry {
        sep_deg: Some(sep_deg),
        outcome: None,
    };

    let (img_pos, cat_pos, n_matched, n_raw) = match params.method {
        SolveMethod::Quads => {
            let mut cat_quads = build_quads_presorted(&cat_star_list, ctx.nrstars_image);
            if cat_quads.is_empty() {
                return failed;
            }
            crate::quads::r#match::sort_catalog_quads(&mut cat_quads);
            let raw = find_matches_sorted(ctx.img_quads, &cat_quads, params.quad_tolerance);
            let n_raw = raw.len();
            log::info!("Found {} references", n_raw);
            let mut filtered = vote_filter(ctx.img_quads, &cat_quads, &raw, params.quad_tolerance);
            if filtered.len() < ctx.min_quads {
                let (by_scale, _) = filter_by_scale(&raw, params.quad_tolerance);
                if by_scale.len() > filtered.len() {
                    filtered = by_scale;
                }
            }
            if filtered.len() < ctx.min_quads {
                return failed;
            }
            let (ip, cp) = extract_star_pairs(ctx.img_quads, &cat_quads, &filtered);
            (ip, cp, filtered.len(), n_raw)
        }
        SolveMethod::Tetra => {
            let cat_tris = build_triangles(&cat_star_list);
            if cat_tris.is_empty() {
                return failed;
            }
            let tol = params.quad_tolerance * TETRA_TOL_FACTOR;
            let raw = find_triangle_matches(ctx.img_tris, &cat_tris, tol);
            let n_raw = raw.len();
            log::info!("Found {} triangle references", n_raw);
            let biject = bijective_filter(&raw, ctx.img_tris, &cat_tris);
            let (filtered, _) = filter_triangles_by_scale(&biject, params.quad_tolerance);
            if filtered.len() < ctx.min_quads {
                return failed;
            }
            let (ip, cp) = extract_triangle_pairs(ctx.img_tris, &cat_tris, &filtered);
            let (ip, cp) = sigma_clip_pairs(ip, cp, 3.0, ctx.min_quads);
            if ip.len() < ctx.min_quads {
                return failed;
            }
            let n_clean = ip.len();
            (ip, cp, n_clean, n_raw)
        }
    };

    let plate = match solve_plate_constants(&img_pos, &cat_pos) {
        Ok(p) => p,
        Err(_) => return failed,
    };

    let (plate, n_verified, rms) = match verify_and_refit(
        ctx.stars,
        &cat_star_list,
        &plate,
        ctx.img.width,
        ctx.img.height,
    ) {
        Some(v) => v,
        None => {
            log::info!("Verification failed at this position; continuing search.");
            return failed;
        }
    };
    log::info!(
        "Verified {} stars against the catalogue, residual {:.2}\"",
        n_verified,
        rms
    );

    PositionTry {
        sep_deg: Some(sep_deg),
        outcome: Some(PositionOutcome {
            idx,
            ra_db,
            dec_db,
            sep_deg,
            plate,
            n_verified,
            rms,
            n_matched,
            n_raw,
            mag_limit,
        }),
    }
}

/// Solve the WCS for an image using the .1476 catalog.
///
/// All progress is emitted via the `log` crate at INFO level — callers install
/// whichever logger backend they need (file, stderr, both, or none).
pub fn solve_image(img: &crate::types::ImageBuffer, params: &SolveParams) -> Result<WcsSolution> {
    // Check the database up front. Every spiral position swallows a missing-file
    // error as "nothing catalogued here", so without this a wrong -d/-D reads
    // nothing everywhere and surfaces as InsufficientQuads - exit 1, "no
    // solution" - when the image is fine and the database is the problem.
    if !crate::catalog::catalog_present(&params.db_path, &params.db_name) {
        return Err(ArcsecError::CatalogNotFound(params.db_path.clone()));
    }

    // --- Phase A: star detection ---
    let bg = get_background(img, params.max_stars);
    log::info!("Start finding stars");
    let (stars, stars_raw) = find_stars_with_background(
        img,
        &bg,
        params.hfd_min,
        params.max_stars,
        img.width,
        img.height,
    );
    log::info!(
        "{} stars found of the requested {}. Background value is {:.0}. \
         Detection level used {:.0} above background. Star level is {:.0} above background. \
         Noise level is {:.0}",
        stars_raw,
        params.max_stars,
        bg.mean,
        bg.star_level,
        bg.star_level,
        bg.noise,
    );
    if stars_raw > params.max_stars {
        log::info!("Selecting the {} brightest stars only.", params.max_stars);
    }

    // Detection is not trimmed. Stars beyond `-s` are faint enough to be absent
    // from the catalog, which once corrupted 3-NN quads badly enough to justify
    // dropping all but the brightest half; quad redundancy and star-level
    // verification absorb that now, and halving the list halved the quad count.
    // The catalog still reads the full requested depth.

    let nrstars_image = stars.len();
    if nrstars_image < 5 {
        return Err(ArcsecError::InsufficientStars {
            found: nrstars_image,
            required: 5,
        });
    }

    // --- Phase B: image pattern building ---
    let img_quads = build_quads(&stars, nrstars_image);
    let nr_quads = img_quads.len();

    let img_tris = if params.method == SolveMethod::Tetra {
        build_triangles(&stars)
    } else {
        crate::quads::TriangleList::default()
    };

    let patterns_empty = match params.method {
        SolveMethod::Quads => nr_quads == 0,
        SolveMethod::Tetra => img_tris.is_empty(),
    };
    if patterns_empty {
        return Err(ArcsecError::InsufficientQuads {
            found: 0,
            required: 3,
        });
    }

    let min_quads: usize = 3 + nrstars_image / 140;

    let oversize: f64 = if nrstars_image < 35 {
        2.0
    } else if nrstars_image > 140 {
        1.0
    } else {
        2.0 * (35.0 / nrstars_image as f64).sqrt()
    };

    // Use the full catalog depth regardless of how many image stars we trimmed.
    let nrstars_required = (params.max_stars as f64 * oversize * oversize).round() as usize;
    let step_size = params.fov;
    let fov_deg = step_size.to_degrees();
    let max_distance = (params.search_radius / step_size + 2.0) as i32;

    log::info!(
        "{} stars, {} quads selected in the image. {} database stars, {} database quads required \
         for the {:.2}d square search window. Step size {:.2}d. Oversize {:.2}",
        nrstars_image,
        nr_quads,
        nrstars_required,
        nrstars_required,
        fov_deg * oversize,
        fov_deg,
        oversize,
    );

    // --- Phase C: spiral search ---
    //
    // Spiral positions are independent, so they are evaluated a batch at a time
    // across a thread pool. Semantics are unchanged from the serial search: within a
    // batch the lowest spiral index wins, and batches are processed in order, so the
    // position returned is exactly the one the serial loop would have returned. The
    // only cost is evaluating the rest of a batch after its first success.
    let ctx = SpiralCtx {
        params,
        img,
        stars: &stars,
        img_quads: &img_quads,
        img_tris: &img_tris,
        nrstars_image,
        nrstars_required,
        oversize,
        min_quads,
        step_size,
    };

    let n_threads = if params.threads > 0 {
        params.threads
    } else {
        crate::max_threads()
    }
    .clamp(1, 64);

    let positions: Vec<(i32, i32)> = SpiralSearch::new(max_distance).collect();
    let mut step_distances: Vec<f64> = Vec::new();

    let mut winner: Option<PositionOutcome> = None;
    let mut start_idx = 0usize;
    while start_idx < positions.len() && winner.is_none() {
        // The first position is the hint itself and usually solves outright, so try it
        // on its own: spawning a pool for it would cost more than it saves.
        let batch_len = if start_idx == 0 {
            1
        } else {
            n_threads.min(positions.len() - start_idx)
        };
        let batch = &positions[start_idx..start_idx + batch_len];

        let tries: Vec<PositionTry> = if n_threads == 1 || batch.len() == 1 {
            batch
                .iter()
                .enumerate()
                .map(|(k, &(sx, sy))| try_position(&ctx, start_idx + k, sx, sy))
                .collect()
        } else {
            std::thread::scope(|scope| {
                let handles: Vec<_> = batch
                    .iter()
                    .enumerate()
                    .map(|(k, &(sx, sy))| {
                        let ctx = &ctx;
                        scope.spawn(move || try_position(ctx, start_idx + k, sx, sy))
                    })
                    .collect();
                handles
                    .into_iter()
                    // A dead worker must not read as "nothing matched here":
                    // the spiral would move on and the solve would fail for a
                    // reason with no trace anywhere.
                    .map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
                    .collect()
            })
        };

        for t in tries {
            if let Some(d) = t.sep_deg {
                step_distances.push(d);
            }
            if let Some(o) = t.outcome
                && winner.as_ref().is_none_or(|w| o.idx < w.idx)
            {
                winner = Some(o);
            }
        }

        start_idx += batch_len;
    }

    if let Some(o) = winner {
        log::info!(
            "{} of {} patterns selected matching within {:.3} tolerance.",
            o.n_matched,
            o.n_raw,
            params.quad_tolerance,
        );

        let mut wcs = derive_wcs(o.ra_db, o.dec_db, &o.plate, img.width, img.height);
        if params.binning > 1 {
            let b = params.binning as f64;
            wcs.crpix1 = (wcs.crpix1 - 0.5) * b + 0.5;
            wcs.crpix2 = (wcs.crpix2 - 0.5) * b + 0.5;
            wcs.cd1_1 /= b;
            wcs.cd1_2 /= b;
            wcs.cd2_1 /= b;
            wcs.cd2_2 /= b;
            wcs.cdelt1 /= b;
            wcs.cdelt2 /= b;
        }
        wcs.residual_rms = o.rms;
        wcs.stars_matched = o.n_verified;
        wcs.raw_matches = o.n_raw;
        wcs.plate = o.plate;
        wcs.mag_limit = o.mag_limit;
        wcs.search_dist_deg = o.sep_deg;
        wcs.step_distances = step_distances;
        return Ok(wcs);
    }

    Err(ArcsecError::InsufficientQuads {
        found: 0,
        required: min_quads,
    })
}

/// Format RA (radians) and Dec (radians) as ASTAP-style "HH: MM  SS.S ±DDd MM  SS".
pub fn format_radec(ra_rad: f64, dec_rad: f64) -> String {
    let ra_h = ra_rad.to_degrees() / 15.0;
    let ra_h = ra_h.rem_euclid(24.0);
    let h = ra_h as u32;
    let ra_min = (ra_h - h as f64) * 60.0;
    let m = ra_min as u32;
    let s = (ra_min - m as f64) * 60.0;

    let dec_deg = dec_rad.to_degrees();
    let sign = if dec_deg < 0.0 { '-' } else { '+' };
    let dec_abs = dec_deg.abs();
    let dd = dec_abs as u32;
    let dec_min = (dec_abs - dd as f64) * 60.0;
    let dm = dec_min as u32;
    let ds = (dec_min - dm as f64) * 60.0;

    format!("{h}: {m:02}  {s:.1} {sign}{dd}d {dm:02}  {ds:.0}")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::coords::{ang_sep, standard_equatorial};
    use crate::types::{ImageBuffer, PlateConstants};
    use crate::wcs::output::derive_wcs;
    use core::f64::consts::PI;

    fn deg(d: f64) -> f64 {
        d * PI / 180.0
    }

    fn make_test_scene(
        n_stars: usize,
        ra_center: f64,
        dec_center: f64,
        cdelt_arcsec: f64,
        width: usize,
        height: usize,
    ) -> (ImageBuffer, Vec<(f64, f64)>, PlateConstants) {
        let mut data = vec![100.0f32; width * height];
        let mut catalog_sky: Vec<(f64, f64)> = Vec::new();
        let stars_per_row = (n_stars as f64).sqrt().ceil() as usize;
        let spacing = 40.0;
        let cx = (width as f64 - 1.0) / 2.0;
        let cy = (height as f64 - 1.0) / 2.0;
        let a = cdelt_arcsec;
        let c = -a * cx;
        let e = cdelt_arcsec;
        let f_offset = -e * cy;
        let plate = PlateConstants {
            a,
            b: 0.0,
            c,
            d: 0.0,
            e,
            f: f_offset,
        };
        let mut count = 0;
        'outer: for row in 0..stars_per_row {
            for col in 0..stars_per_row {
                if count >= n_stars {
                    break 'outer;
                }
                let px = 20.0 + col as f64 * spacing;
                let py = 20.0 + row as f64 * spacing;
                if px >= width as f64 - 20.0 || py >= height as f64 - 20.0 {
                    continue;
                }
                let x_std = a * px + c;
                let y_std = e * py + f_offset;
                let (ra, dec) = standard_equatorial(ra_center, dec_center, x_std, y_std, 1.0);
                catalog_sky.push((ra, dec));
                let sigma = 2.0;
                let amp = 30000.0f32;
                for dy in -8i32..=8 {
                    for dx in -8i32..=8 {
                        let x = (px as i32 + dx) as usize;
                        let y = (py as i32 + dy) as usize;
                        if x < width && y < height {
                            let r2 = (dx * dx + dy * dy) as f64 / (2.0 * sigma * sigma);
                            data[y * width + x] += amp * (-r2).exp() as f32;
                        }
                    }
                }
                count += 1;
            }
        }
        let img = ImageBuffer {
            data,
            width,
            height,
        };
        (img, catalog_sky, plate)
    }

    #[test]
    fn derive_wcs_recovers_position() {
        let ra_center = deg(45.0);
        let dec_center = deg(30.0);
        let (img, _cat, plate) = make_test_scene(16, ra_center, dec_center, 2.0, 300, 300);
        let wcs = derive_wcs(ra_center, dec_center, &plate, img.width, img.height);
        let sep_arcsec = ang_sep(wcs.ra0, wcs.dec0, ra_center, dec_center) * (180.0 / PI * 3600.0);
        assert!(sep_arcsec < 0.5, "centre offset = {sep_arcsec} arcsec");
    }

    #[test]
    fn spiral_covers_origin_first() {
        assert_eq!(SpiralSearch::new(5).next(), Some((0, 0)));
    }

    #[test]
    fn oversize_formula_limits() {
        for n in [10, 35, 70, 140, 200] {
            let ov: f64 = if n < 35 {
                2.0
            } else if n > 140 {
                1.0
            } else {
                2.0 * (35.0 / n as f64).sqrt()
            };
            assert!((1.0..=2.0).contains(&ov), "oversize={ov} for n={n}");
        }
    }

    #[test]
    fn format_radec_roundtrip() {
        let s = format_radec(deg(160.875), deg(-59.524));
        assert!(s.contains("10:"), "RA hours: {s}");
        assert!(s.contains('-'), "dec sign: {s}");
    }
}
