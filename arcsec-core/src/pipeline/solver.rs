//! Full plate-solving pipeline.

use core::f64::consts::PI;
use std::path::PathBuf;

use crate::catalog::read_catalog_stars;
use crate::catalog::{CatalogLayout, CatalogStar};
use crate::detection::get_background;
use crate::detection::stars::find_stars_with_background;
use crate::error::{ArcsecError, Result};
use crate::math::coords::{ang_sep, equatorial_standard, standard_equatorial};
use crate::math::lsq::{fit_affine, solve_plate_constants};
use crate::quads::{
    TETRA_TOL_FACTOR, bijective_filter, build_quads, build_quads_presorted, build_triangles,
    extract_star_pairs, extract_triangle_pairs, filter_by_scale, filter_triangles_by_scale,
    find_matches_sorted, find_triangle_matches, vote_filter,
};
use crate::types::{MatchedStar, PairedPositions, PlateConstants, Star, StarList, WcsSolution};
use crate::wcs::output::derive_wcs;

use super::distortion::{Pair, Refined, StarGrid, best_linear, max_departure_px, refine};
use super::spiral::SpiralSearch;

/// Which pattern-matching algorithm to use in the catalog spiral loop.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SolveMethod {
    /// ASTAP-style 5-ratio quad matching with `vote_filter` (default).
    #[default]
    Quads,
    /// TETRA 2-ratio triangle matching with bijective filter.
    Tetra,
}

/// How much sky the spiral search reads around each position (ASTAP's `-speed`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SearchSpeed {
    /// Size the catalogue window by the image's star count: twice the field for an
    /// image with fewer than 35 stars, falling to the field itself above 140.
    #[default]
    Auto,
    /// Always read a window twice the field, so neighbouring spiral positions
    /// overlap and a field that straddles two of them is still seen whole. Four
    /// times the catalogue stars per position, so slower; it helps images with
    /// many stars that the auto window only just misses.
    Slow,
}

/// Parameters for [`solve_image`].
#[derive(Debug, Clone)]
pub struct SolveParams {
    /// Approximate RA of image centre (radians, hint only).
    pub ra_hint: f64,
    /// Approximate DEC of image centre (radians, hint only).
    pub dec_hint: f64,
    /// Image field of view along its longer side (radians). Used as the spiral step
    /// size, and, with the image's size, for the pixel scale a solution is expected
    /// to have.
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
    /// Catalogue window per spiral position.
    pub speed: SearchSpeed,
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
/// large residuals of false pattern matches. Subsequent passes apply
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
        // Unchecked: the first fit is made on the contaminated set, and gross
        // outliers can skew it past solve_plate_constants' similarity check even though
        // clipping them is exactly what would fix it.
        let Ok(plate) = fit_affine(&img_pos, &cat_pos) else {
            break;
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
            // Gross outliers drag a least-squares fit towards themselves and inflate
            // every residual, so a fixed cut can keep them. The median residual is
            // not moved by a minority of outliers: allow 3 sigma of it (1.4826 x the
            // median absolute residual estimates sigma) when that is larger.
            let mut sorted = residuals.clone();
            sorted.sort_unstable_by(f64::total_cmp);
            let median = sorted[sorted.len() / 2];
            (10.0 * cdelt).max(10.0).max(3.0 * 1.4826 * median)
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

/// Fit the plate to the matched pattern centroids, sigma-clipping them first.
///
/// Returns the plate and the number of pairs it was fitted to, or `None` if fewer
/// than `min_count` survive the clipping or the fit is not a similarity.
///
/// The clipping matters for both methods. A few wrong patterns that land in the
/// winning vote cell drag the unweighted fit far enough that the similarity check
/// refuses it, and the position is abandoned although most patterns agree. On the
/// Coalsack (`obj_coalsack`) the first position, the right one, kept 18 quads whose
/// plain fit was refused; clipped, the fit verified 91 stars. With no wrong pairs
/// the clipping removes little, and the star-level refit that follows makes the
/// difference immaterial.
fn fit_pattern_pairs(
    img_pos: Vec<(f64, f64)>,
    cat_pos: Vec<(f64, f64)>,
    min_count: usize,
) -> Option<(PlateConstants, usize)> {
    let (img_pos, cat_pos) = sigma_clip_pairs(img_pos, cat_pos, 3.0, min_count);
    if img_pos.len() < min_count {
        return None;
    }
    let plate = solve_plate_constants(&img_pos, &cat_pos).ok()?;
    Some((plate, img_pos.len()))
}

/// Minimum number of individually matched stars required to believe a solution.
///
/// Correct solves typically match 200-375 stars, so this is deliberately loose;
/// its job is to reject the handful-of-coincidences case. Together with
/// `MIN_VERIFY_SPREAD` it separates two otherwise identical-looking results: M31 at
/// 2 degrees (22 stars, spread 0.221, rms 0.65", rotation wrong by 1.56 degrees)
/// from the Dec -88 field (46 stars, spread 0.207, rms 0.66", correct to 2.3").
const MIN_VERIFIED_STARS: usize = 30;
/// Fewest verified stars accepted for an image with `nrstars_image` detections:
/// [`MIN_VERIFIED_STARS`], relaxed to 15% of the detections for a sparse image, but
/// never below 10. Unchanged from 200 detections up.
///
/// A sparse frame (a short exposure, a narrow band, a small field) cannot match 30
/// stars when it shows only 40. Below 30 the solution must also pass
/// [`Acceptance::accepts`]' scale and residual checks.
fn min_verified_stars(nrstars_image: usize) -> usize {
    nrstars_image
        .saturating_mul(15)
        .div_ceil(100)
        .clamp(10, MIN_VERIFIED_STARS)
}
/// Below [`MIN_VERIFIED_STARS`] matched stars, the fitted pixel scale must be
/// within this fraction of the one the hint implies.
const RELAXED_SCALE_TOL: f64 = 0.10;
/// Below [`MIN_VERIFIED_STARS`] matched stars, the largest star-level rms, in
/// pixels. The genuine relaxed solves on the corpus have 0.16-0.46 px; the false
/// positive the relaxed count alone lets through (`ls2_25`, 12 stars) has 2.9 px,
/// and a scale 1.36 times the hint's.
const RELAXED_MAX_RMS_PX: f64 = 0.5;

/// When a verified plate is believed.
struct Acceptance {
    /// Fewest matched stars ([`min_verified_stars`]).
    min_stars: usize,
    /// The pixel scale the hint implies, arcsec per (binned) pixel.
    expected_scale: f64,
}

impl Acceptance {
    fn new(nrstars_image: usize, params: &SolveParams, img: &crate::types::ImageBuffer) -> Self {
        Self {
            min_stars: min_verified_stars(nrstars_image),
            expected_scale: params.fov.to_degrees() * 3600.0
                / img.width.max(img.height).max(1) as f64,
        }
    }

    /// Enough stars, spread over the frame; and if fewer than
    /// [`MIN_VERIFIED_STARS`], the hint's pixel scale and a sub-half-pixel fit.
    ///
    /// A handful of chance coincidences can be fitted by some plate at some scale;
    /// they are not fitted at the scale the optics give, to a fraction of a pixel.
    fn accepts(&self, v: &Verified, spread: f64) -> bool {
        if v.n() < self.min_stars || spread < MIN_VERIFY_SPREAD {
            return false;
        }
        if v.n() >= MIN_VERIFIED_STARS {
            return true;
        }
        let p = &v.plate;
        let scale = (p.a * p.e - p.b * p.d).abs().sqrt();
        let ok = (scale / self.expected_scale - 1.0).abs() <= RELAXED_SCALE_TOL
            && v.rms <= RELAXED_MAX_RMS_PX * scale;
        log::info!(
            "{} stars verified, scale {:.4}\"/px against {:.4} expected, residual {:.2} px: {}",
            v.n(),
            scale,
            self.expected_scale,
            v.rms / scale,
            if ok { "accepted" } else { "refused" }
        );
        ok
    }
}

/// Match radii (pixels) used by successive verification passes, coarse to fine.
const VERIFY_RADII: [f64; 3] = [6.0, 3.0, 2.0];
/// Minimum spread of the matched stars, as a fraction of the image half-diagonal.
///
/// A count threshold alone is not enough: matches clustered in one part of the
/// frame (the core of a bright galaxy, say) pin the position but leave rotation
/// and scale essentially free. M31 at 2 degrees passed with 22 matched stars and
/// a 1.56-degree rotation error, which is 154" at the field corners.
const MIN_VERIFY_SPREAD: f64 = 0.20;

/// A verified plate: the re-fitted plate constants, the per-star RMS in arcsec,
/// and the star pairs the fit was made from.
struct Verified {
    plate: PlateConstants,
    rms: f64,
    /// Detected star positions, pixels of the solved (binned) image, 0-based.
    img_pos: Vec<(f64, f64)>,
    /// The catalogue star each was paired with, in standard coordinates (arcsec)
    /// about the plane the plate maps into.
    cat_pos: Vec<(f64, f64)>,
}

impl Verified {
    /// Number of individually matched stars.
    fn n(&self) -> usize {
        self.img_pos.len()
    }
}

/// Project the catalogue onto the image with a candidate plate solution, match
/// individual stars, and re-fit on those matches.
///
/// The quad matcher only ever produces quad *centroids*, so the plate fit is built
/// from a handful of averaged positions and nothing ever checks that the individual
/// stars agree. This does that check: invert the plate to map every catalogue star
/// into pixel space, pair each with the nearest detected star, re-fit on the pairs,
/// and repeat with a shrinking radius.
///
/// Returns the refined plate with its matched pairs, or `None` if the plate is
/// degenerate or `accept` refuses the result (too few stars agree, the matches
/// are too clustered, or a sparse match is at the wrong scale or fits loosely).
fn verify_and_refit(
    img_stars: &StarList,
    cat_stars: &StarList,
    plate: &PlateConstants,
    img_w: usize,
    img_h: usize,
    accept: &Acceptance,
) -> Option<Verified> {
    if img_stars.is_empty() || cat_stars.is_empty() {
        return None;
    }

    // Uniform grid over the detected stars for nearest-neighbour lookup.
    let grid = StarGrid::new(img_stars, VERIFY_RADII[0])?;

    let mut current = plate.clone();
    // The last pass that fitted, with the spread of its matches.
    let mut best: Option<(Verified, f64)> = None;

    for &radius in &VERIFY_RADII {
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
            if !grid.near(px, py, radius) {
                continue;
            }
            if let Some(i) = grid.nearest(px, py, r2, &used) {
                used[i] = true; // one-to-one: a detected star backs at most one catalogue star
                img_pos.push(grid.pos(i));
                cat_pos.push((cs.x, cs.y));
            }
        }

        if img_pos.len() < 4 {
            break;
        }
        let Ok(refined) = solve_plate_constants(&img_pos, &cat_pos) else {
            break;
        };
        let mut sq = 0.0;
        for (&(xi, yi), &(xc, yc)) in img_pos.iter().zip(cat_pos.iter()) {
            let xp = refined.a * xi + refined.b * yi + refined.c;
            let yp = refined.d * xi + refined.e * yi + refined.f;
            sq += (xp - xc).powi(2) + (yp - yc).powi(2);
        }
        let rms = (sq / img_pos.len() as f64).sqrt();
        // Matches clustered in one corner leave rotation free.
        let spread = spread_of(&img_pos, img_w, img_h);
        log::debug!(
            "verify: {} stars, spread {:.3}, rms {:.2}\"",
            img_pos.len(),
            spread,
            rms
        );

        current = refined.clone();
        best = Some((
            Verified {
                plate: refined,
                rms,
                img_pos,
                cat_pos,
            },
            spread,
        ));
    }

    best.filter(|(v, spread)| accept.accepts(v, *spread))
        .map(|(v, _)| v)
}

/// The most image stars the database can match in this field: its density
/// ([`crate::catalog::database_density`]) times the image's area, or `-s` if that
/// is smaller or the density is unknown. ASTAP's "database limit".
///
/// The area is `fov²` scaled by the aspect ratio, `fov` being the long side.
fn density_star_limit(params: &SolveParams, img: &crate::types::ImageBuffer) -> usize {
    let Some(density) = crate::catalog::database_density(&params.db_name) else {
        return params.max_stars;
    };
    let fov_deg = params.fov.to_degrees();
    let (w, h) = (img.width as f64, img.height as f64);
    let area = fov_deg * fov_deg * w.min(h) / w.max(h).max(1.0);
    let cap = (density * area).round();
    if cap < params.max_stars as f64 {
        cap as usize
    } else {
        params.max_stars
    }
}

/// Everything a spiral position needs that does not change between positions.
struct SpiralCtx<'a> {
    params: &'a SolveParams,
    img: &'a crate::types::ImageBuffer,
    stars: &'a StarList,
    img_quads: &'a crate::types::QuadList,
    img_tris: &'a crate::quads::TriangleList,
    nrstars_image: usize,
    /// The most image stars worth using: `-s`, or fewer if the database cannot hold
    /// that many in the field ([`density_star_limit`]).
    star_limit: usize,
    nrstars_required: usize,
    oversize: f64,
    min_quads: usize,
    step_size: f64,
    accept: Acceptance,
    /// Long side over short side of the image.
    aspect: f64,
}

/// A spiral position that produced a verified solution.
struct PositionOutcome {
    idx: usize,
    ra_db: f64,
    dec_db: f64,
    sep_deg: f64,
    verified: Verified,
    n_matched: usize,
    n_raw: usize,
    mag_limit: f64,
    /// The field is distorted beyond what the model can follow: the search
    /// stops here, without a solution.
    refused: bool,
}

/// Result of trying one spiral position: the angular distance if the catalogue was
/// actually read there (for the ASTAP-style progress line), and the solution if one
/// verified.
struct PositionTry {
    sep_deg: Option<f64>,
    outcome: Option<PositionOutcome>,
}

impl PositionTry {
    const NONE: Self = Self {
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

    // Any read failure (a missing tile included) counts as "nothing catalogued
    // here"; `solve_image` has already checked that the database exists at all.
    let cat_raw = match read_catalog_stars(
        &params.db_path,
        &params.db_name,
        ra_db,
        dec_db,
        params.fov * ctx.oversize,
        ctx.nrstars_required,
    ) {
        Ok(v) if !v.is_empty() => v,
        Ok(_) | Err(_) => return PositionTry::NONE,
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

    let (img_pos, cat_pos, n_raw) = match params.method {
        SolveMethod::Quads => {
            let mut cat_quads = build_quads_presorted(&cat_star_list, ctx.nrstars_image);
            if ctx.nrstars_image < ctx.star_limit {
                add_density_matched_quads(ctx, &cat_raw, ra_db, dec_db, &mut cat_quads);
            }
            if cat_quads.is_empty() {
                return failed;
            }
            crate::quads::r#match::sort_catalog_quads(&mut cat_quads);
            let raw = find_matches_sorted(ctx.img_quads, &cat_quads, params.quad_tolerance);
            let n_raw = raw.len();
            log::info!("Found {n_raw} references");
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
            (ip, cp, n_raw)
        }
        SolveMethod::Tetra => {
            let cat_tris = build_triangles(&cat_star_list);
            if cat_tris.is_empty() {
                return failed;
            }
            let tol = params.quad_tolerance * TETRA_TOL_FACTOR;
            let raw = find_triangle_matches(ctx.img_tris, &cat_tris, tol);
            let n_raw = raw.len();
            log::info!("Found {n_raw} triangle references");
            let biject = bijective_filter(&raw, ctx.img_tris, &cat_tris);
            let (filtered, _) = filter_triangles_by_scale(&biject, params.quad_tolerance);
            if filtered.len() < ctx.min_quads {
                return failed;
            }
            let (ip, cp) = extract_triangle_pairs(ctx.img_tris, &cat_tris, &filtered);
            (ip, cp, n_raw)
        }
    };

    // The matched patterns, kept as seeds for the distortion model: wherever in
    // the frame they fall, they are correspondences whatever the plate does there.
    let seeds = Seeds {
        img: img_pos.clone(),
        cat: cat_pos.clone(),
        ra: ra_db,
        dec: dec_db,
    };
    let Some((plate, n_matched)) = fit_pattern_pairs(img_pos, cat_pos, ctx.min_quads) else {
        return failed;
    };

    let found = |verified, ra_db, dec_db, refused| PositionTry {
        sep_deg: Some(sep_deg),
        outcome: Some(PositionOutcome {
            idx,
            ra_db,
            dec_db,
            sep_deg,
            verified,
            n_matched,
            n_raw,
            mag_limit,
            refused,
        }),
    };

    let Some(verified) = verify_and_refit(
        ctx.stars,
        &cat_star_list,
        &plate,
        ctx.img.width,
        ctx.img.height,
        &ctx.accept,
    ) else {
        log::info!("Verification failed at this position; continuing search.");
        if n_matched >= STRONG_VOTE
            && let Some((verified, ra_c, dec_c)) =
                second_chance(ctx, &cat_raw, &seeds, &plate, ra_db, dec_db)
        {
            return found(verified, ra_c, dec_c, false);
        }
        return failed;
    };
    log::info!(
        "Verified {} stars against the catalogue, residual {:.2}\"",
        verified.n(),
        verified.rms
    );

    let (verified, ra_db, dec_db) = recentre(ctx, &cat_raw, verified, ra_db, dec_db);
    match model_distortion(ctx, &cat_raw, &seeds, verified, ra_db, dec_db) {
        Modelled::Linear(v) => found(v, ra_db, dec_db, false),
        Modelled::Distorted(v, ra_c, dec_c) => found(v, ra_c, dec_c, false),
        Modelled::Refused(v) => found(v, ra_db, dec_db, true),
    }
}

/// Pattern pairs a position must keep, after clipping, for [`second_chance`].
///
/// A wrong position rarely keeps more than a handful; the two corpus images the
/// second chance solves kept 60 and more.
const STRONG_VOTE: usize = 50;

/// Catalogue stars in standard coordinates about `(ra, dec)`, brightest first.
fn project(cat_raw: &[CatalogStar], ra: f64, dec: f64) -> StarList {
    StarList(
        cat_raw
            .iter()
            .map(|s| {
                let (x, y) = equatorial_standard(ra, dec, s.ra, s.dec, 1.0);
                Star {
                    x,
                    y,
                    snr: 1.0,
                    hfd: 2.0,
                }
            })
            .collect(),
    )
}

/// The matched pattern centroids of a position, in the plane about `(ra, dec)`.
struct Seeds {
    img: Vec<(f64, f64)>,
    cat: Vec<(f64, f64)>,
    ra: f64,
    dec: f64,
}

impl Seeds {
    /// The pairs with their catalogue side moved to the plane about `(ra, dec)`.
    fn in_plane(&self, ra: f64, dec: f64) -> Vec<Pair> {
        self.img
            .iter()
            .zip(&self.cat)
            .map(|(&i, &(x, y))| {
                if ra == self.ra && dec == self.dec {
                    return (i, (x, y));
                }
                let (sra, sdec) = standard_equatorial(self.ra, self.dec, x, y, 1.0);
                (i, equatorial_standard(ra, dec, sra, sdec, 1.0))
            })
            .collect()
    }
}

/// Fit the distortion model ([`refine`]) from a linear plate in the plane about
/// `(ra, dec)`.
fn fit_distortion(
    ctx: &SpiralCtx<'_>,
    cat_raw: &[CatalogStar],
    seeds: &Seeds,
    plate: &PlateConstants,
    ra: f64,
    dec: f64,
) -> Option<Refined> {
    // The spiral position's catalogue window need not cover the image: with the
    // hint a third of a field off, a third of the frame has no catalogue stars,
    // and the model cannot reach it. So read the catalogue again about the image
    // centre, wide enough for its corners at any rotation, at the same density.
    // One read per solve (or per strong position that failed to verify).
    let (w, h) = (ctx.img.width as f64, ctx.img.height as f64);
    let (xs, ys) = (
        plate.a * (w - 1.0) * 0.5 + plate.b * (h - 1.0) * 0.5 + plate.c,
        plate.d * (w - 1.0) * 0.5 + plate.e * (h - 1.0) * 0.5 + plate.f,
    );
    let (ra_c, dec_c) = standard_equatorial(ra, dec, xs, ys, 1.0);
    let window = w.hypot(h) / w.max(h);
    let cat = match read_catalog_stars(
        &ctx.params.db_path,
        &ctx.params.db_name,
        ra_c,
        dec_c,
        ctx.params.fov * window,
        (ctx.params.max_stars as f64 * window * window).round() as usize,
    ) {
        Ok(v) if !v.is_empty() => project(&v, ra, dec),
        _ => project(cat_raw, ra, dec),
    };
    let grid = StarGrid::new(ctx.stars, VERIFY_RADII[0])?;
    let r = refine(
        &grid,
        &cat,
        &seeds.in_plane(ra, dec),
        plate,
        ctx.img.width,
        ctx.img.height,
        VERIFY_RADII[VERIFY_RADII.len() - 1],
    )?;
    log::info!(
        "Distortion model: {} terms, {} stars within {} px, rms {:.2} px, F {:.1} over linear, {} of 9 cells",
        r.model.n_terms,
        r.img_pos.len(),
        VERIFY_RADII[VERIFY_RADII.len() - 1],
        r.rms / r.model.scale(),
        r.f_linear,
        r.cells
    );
    Some(r)
}

/// The linear plate closest to `r`'s model over the frame ([`best_linear`]), in
/// the tangent plane at the image centre, with the model's star pairs carried into
/// that plane. Returns it as a verification record, with the new tangent point.
fn linear_from_model(
    ctx: &SpiralCtx<'_>,
    r: &Refined,
    ra: f64,
    dec: f64,
) -> Option<(Verified, f64, f64)> {
    let (w, h) = (ctx.img.width, ctx.img.height);
    let (xs, ys) = r
        .model
        .apply((w as f64 - 1.0) * 0.5, (h as f64 - 1.0) * 0.5);
    let (ra_c, dec_c) = standard_equatorial(ra, dec, xs, ys, 1.0);
    // The plate of the model's plane is not linear in the centre's; carry the
    // model across point by point.
    let moved = |(x, y): (f64, f64)| {
        let (sra, sdec) = standard_equatorial(ra, dec, x, y, 1.0);
        equatorial_standard(ra_c, dec_c, sra, sdec, 1.0)
    };
    let plate = best_linear(|x, y| moved(r.model.apply(x, y)), w, h)?;
    let cat_pos = r.cat_pos.iter().map(|&p| moved(p)).collect();
    Some((
        Verified {
            plate,
            rms: r.rms,
            img_pos: r.img_pos.clone(),
            cat_pos,
        },
        ra_c,
        dec_c,
    ))
}

/// What [`model_distortion`] made of a verified position.
enum Modelled {
    /// No distortion worth reporting: the verified plate, unchanged.
    Linear(Verified),
    /// The linear plate closest to the distortion model, with its tangent point.
    Distorted(Verified, f64, f64),
    /// Strong distortion the model cannot follow over the frame: the linear plate
    /// would be wrong at the edges, so the solve is refused.
    Refused(Verified),
}

/// F statistic of the distortion model over a linear plate, on the final pairs,
/// needed for it to change the reported plate.
///
/// On the corpus, survey images with no distortion reach F = 4–24 (catalogue and
/// centroid systematics a cubic can fit); every real distortion that changes a
/// result is above 30 except one 2.2° TESS crop (F = 7). See test-images.md §7.9.
const MIN_REPORT_F: f64 = 30.0;

/// A distortion model changes the reported plate only if the verified plate is
/// this far (pixels) from it somewhere in the frame. ZTF's real distortion of under
/// a pixel (F up to 47) stays below it, and those solves are unchanged.
const MIN_DEPARTURE_PX: f64 = 1.0;

/// When the model cannot be used (its pairs leave part of the frame empty), a
/// cubic fitted to the wide-radius pairs at least this significant ...
const REFUSE_F: f64 = 100.0;
/// ... and this far (pixels) from the verified plate where it has stars means the
/// linear plate is a fit to part of a strongly distorted field: refuse it. No solve
/// on the corpus comes near (largest F 21 at ≥ 2 px); fields with 4 px or more of
/// distortion and a third of the frame empty all reach both.
const REFUSE_DEPARTURE_PX: f64 = 3.0;

/// After a linear plate verified, fit the distortion model and decide what to
/// report.
///
/// The model replaces the verified plate when it is significant
/// ([`MIN_REPORT_F`]), its pairs cover the frame (all nine cells of a 3×3 grid for
/// a cubic, seven for a quadratic), it moves some part of the frame by
/// [`MIN_DEPARTURE_PX`] or more, and it pairs at least 90% as many stars as the
/// plate. Otherwise the verified plate stands, unless the field is visibly and
/// strongly distorted where it has stars ([`REFUSE_F`], [`REFUSE_DEPARTURE_PX`]).
fn model_distortion(
    ctx: &SpiralCtx<'_>,
    cat_raw: &[CatalogStar],
    seeds: &Seeds,
    verified: Verified,
    ra: f64,
    dec: f64,
) -> Modelled {
    let Some(r) = fit_distortion(ctx, cat_raw, seeds, &verified.plate, ra, dec) else {
        return Modelled::Linear(verified);
    };
    let (w, h) = (ctx.img.width, ctx.img.height);
    let departure = max_departure_px(&r.model, &verified.plate, w, h);
    log::info!(
        "Distortion: verified plate departs {departure:.2} px from the model; {} stars against {} verified",
        r.img_pos.len(),
        verified.n(),
    );
    let min_cells = if r.model.n_terms == 10 { 9 } else { 7 };
    let usable = r.model.n_terms > 3
        && r.cells >= min_cells
        && r.f_linear >= MIN_REPORT_F
        && r.img_pos.len() * 10 >= verified.n() * 9;
    if usable {
        if departure >= MIN_DEPARTURE_PX
            && let Some((v, ra_c, dec_c)) = linear_from_model(ctx, &r, ra, dec)
        {
            log::info!("Reporting the linear plate closest to the distortion model.");
            return Modelled::Distorted(v, ra_c, dec_c);
        }
        return Modelled::Linear(verified);
    }
    let (wide_f, wide_dep) = r.unmodelled();
    log::info!("Where the stars are: a cubic with F {wide_f:.1}, {wide_dep:.2} px from the plate.");
    if wide_f >= REFUSE_F && wide_dep >= REFUSE_DEPARTURE_PX {
        log::info!(
            "The field is distorted by {wide_dep:.1} px where it has stars, and the distortion \
             cannot be modelled over the whole frame: refusing a linear solution."
        );
        return Modelled::Refused(verified);
    }
    Modelled::Linear(verified)
}

/// Spread of matched stars about their centroid, as a fraction of the image
/// half-diagonal.
fn spread_of(img_pos: &[(f64, f64)], img_w: usize, img_h: usize) -> f64 {
    let n = img_pos.len() as f64;
    let mx = img_pos.iter().map(|p| p.0).sum::<f64>() / n;
    let my = img_pos.iter().map(|p| p.1).sum::<f64>() / n;
    let var = img_pos
        .iter()
        .map(|&(x, y)| (x - mx) * (x - mx) + (y - my) * (y - my))
        .sum::<f64>()
        / n;
    let half_diag = 0.5 * ((img_w * img_w + img_h * img_h) as f64).sqrt();
    var.sqrt() / half_diag
}

/// A position whose patterns agree strongly but whose linear plate did not verify:
/// fit the distortion model from the pattern plate and verify that instead, by
/// the same acceptance rules. A strongly distorted field can leave too few stars
/// within 2 pixels of any linear plate (`s_tess_b_pincush`, before the model,
/// verified 11 at the right position and then accepted a neighbour's wrong plate).
fn second_chance(
    ctx: &SpiralCtx<'_>,
    cat_raw: &[CatalogStar],
    seeds: &Seeds,
    plate: &PlateConstants,
    ra: f64,
    dec: f64,
) -> Option<(Verified, f64, f64)> {
    log::info!("Strong pattern match: retrying verification with a distortion model.");
    let r = fit_distortion(ctx, cat_raw, seeds, plate, ra, dec)?;
    // As for a reported model: it must reach the frame it will be extrapolated to.
    if r.cells < if r.model.n_terms == 10 { 9 } else { 7 } {
        log::info!("The distortion model's stars do not cover the frame.");
        return None;
    }
    let probe = Verified {
        plate: r.model.linear_part(),
        rms: r.rms,
        img_pos: r.img_pos.clone(),
        cat_pos: r.cat_pos.clone(),
    };
    let spread = spread_of(&r.img_pos, ctx.img.width, ctx.img.height);
    if !ctx.accept.accepts(&probe, spread) {
        log::info!("The distortion model did not verify either.");
        return None;
    }
    log::info!(
        "Verified {} stars with the distortion model.",
        r.img_pos.len()
    );
    linear_from_model(ctx, &r, ra, dec)
}

/// How many times denser than the image the catalogue read must be before
/// [`add_density_matched_quads`] adds anything.
///
/// Every image the added quads solve on the corpus had a catalogue at least 3.7
/// times denser than itself (the density-matched count at most 0.27 of the read);
/// below 2.5 the full-depth quads already share the image's neighbourhoods closely
/// enough, and the extra quads made a failing search on a 159-star image 50% slower.
const DENSITY_MATCH_MIN_RATIO: f64 = 2.5;

/// Add the quads of a catalogue star list as dense as the image's.
///
/// The catalogue is read to the depth `-s` asks for, so when the image yields fewer
/// stars than that the catalogue is denser than the image, and its quads join
/// neighbours the image never detected. The window's brightest
/// `n · oversize² · long/short` catalogue stars have the image's density (the
/// window is square, `oversize` fields of the long side across), so their quads are
/// built over the same neighbourhoods as the image's. They are added to the
/// full-depth quads, not substituted: reading only that depth loses more images
/// than it gains, those whose faint detections are real (test-images.md §7.8).
/// Verification still runs against the full-depth list.
///
/// Only when the catalogue read is at least [`DENSITY_MATCH_MIN_RATIO`] times
/// denser than the image: the extra quads are matched at every spiral position, so
/// on a search that fails everywhere they cost what they add to the quad count.
fn add_density_matched_quads(
    ctx: &SpiralCtx<'_>,
    cat_raw: &[CatalogStar],
    ra_db: f64,
    dec_db: f64,
    cat_quads: &mut crate::types::QuadList,
) {
    let k = (ctx.nrstars_image as f64 * ctx.oversize * ctx.oversize * ctx.aspect).round() as usize;
    if k < 5 || (k as f64) * DENSITY_MATCH_MIN_RATIO > cat_raw.len() as f64 {
        return;
    }
    // `cat_raw` is brightest first.
    let mut sub: Vec<Star> = cat_raw[..k]
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
    sub.sort_unstable_by(|a, b| a.x.total_cmp(&b.x));
    let extra = build_quads_presorted(&StarList(sub), ctx.nrstars_image);
    // A quad of the same four stars can come out of both lists: count it once.
    // Centroid and size to a milliarcsecond identify it.
    let key = |q: &crate::types::Quad| {
        (
            (q.center_x * 1000.0).round() as i64,
            (q.center_y * 1000.0).round() as i64,
            (q.d1 * 1000.0).round() as i64,
        )
    };
    let seen: std::collections::HashSet<_> = cat_quads.0.iter().map(key).collect();
    let before = cat_quads.len();
    cat_quads
        .0
        .extend(extra.0.into_iter().filter(|q| !seen.contains(&key(q))));
    log::info!(
        "{} more database quads from its {k} brightest stars, the image's density.",
        cat_quads.len() - before
    );
}

/// Refit a verified plate in the tangent plane at the image centre.
///
/// The plate constants are a linear map from pixels to the tangent plane at the
/// spiral position the catalogue was projected about. Pixels map linearly onto a
/// tangent plane only at the optical axis, so away from it the fit absorbs the
/// projection's curvature as a rotation and shear, which grow with the distance
/// from the field and with declination. Star-level RMS stays small, because the fit
/// is good *in that plane*, but the CD matrix derived from it is wrong at the image
/// centre: with the hint 0.3 fields off, most corpus solves were out by 5-1600" at
/// the corners.
///
/// So once a position verifies, move the tangent point to the image centre, pair
/// stars as the verified plate predicts them, fit those pairs in the new plane,
/// and verify again. Twice, since the centre moves slightly with
/// the new fit. If a pass fails to verify, the previous solution is kept: this can
/// only improve a solve, never lose one.
fn recentre(
    ctx: &SpiralCtx<'_>,
    cat_raw: &[CatalogStar],
    mut verified: Verified,
    mut ra_db: f64,
    mut dec_db: f64,
) -> (Verified, f64, f64) {
    let (w, h) = (ctx.img.width as f64, ctx.img.height as f64);
    let (cx, cy) = ((w - 1.0) * 0.5, (h - 1.0) * 0.5);
    let apply =
        |p: &PlateConstants, x: f64, y: f64| (p.a * x + p.b * y + p.c, p.d * x + p.e * y + p.f);

    for _ in 0..2 {
        let plate = &verified.plate;
        let (xs, ys) = apply(plate, cx, cy);
        // Already centred to well under a milliarcsecond: nothing to gain.
        if xs.hypot(ys) < 1e-3 {
            break;
        }
        let (ra0, dec0) = standard_equatorial(ra_db, dec_db, xs, ys, 1.0);

        // Starting plate for the new tangent plane. Mapping one tangent plane onto
        // another is far from linear over a wide field (a 10-degree field 3 degrees
        // off moves by ~65" under a straight-line fit), so rather than carry the
        // plate across, pair stars exactly as the verified plate predicts them and
        // fit those pairs against their positions in the new plane.
        let det = plate.a * plate.e - plate.b * plate.d;
        if det.abs() < 1e-12 {
            break;
        }
        let r2 = VERIFY_RADII[0] * VERIFY_RADII[0];
        let mut used = vec![false; ctx.stars.len()];
        let mut img_pos = Vec::new();
        let mut new_pos = Vec::new();
        let mut cat = Vec::with_capacity(cat_raw.len());
        for s in cat_raw {
            let (nx, ny) = equatorial_standard(ra0, dec0, s.ra, s.dec, 1.0);
            cat.push(Star {
                x: nx,
                y: ny,
                snr: 1.0,
                hfd: 2.0,
            });
            let (ox, oy) = equatorial_standard(ra_db, dec_db, s.ra, s.dec, 1.0);
            let (dx, dy) = (ox - plate.c, oy - plate.f);
            let px = (plate.e * dx - plate.b * dy) / det;
            let py = (-plate.d * dx + plate.a * dy) / det;
            let nearest = ctx
                .stars
                .0
                .iter()
                .enumerate()
                .filter(|&(i, _)| !used[i])
                .map(|(i, st)| (i, (st.x - px).powi(2) + (st.y - py).powi(2)))
                .filter(|&(_, d2)| d2 < r2)
                .min_by(|a, b| a.1.total_cmp(&b.1));
            if let Some((i, _)) = nearest {
                used[i] = true;
                img_pos.push((ctx.stars.0[i].x, ctx.stars.0[i].y));
                new_pos.push((nx, ny));
            }
        }
        let Ok(guess) = solve_plate_constants(&img_pos, &new_pos) else {
            break;
        };
        let cat = StarList(cat);
        let Some(v) = verify_and_refit(
            ctx.stars,
            &cat,
            &guess,
            ctx.img.width,
            ctx.img.height,
            &ctx.accept,
        ) else {
            log::info!("Re-centring on the image centre did not verify; keeping the fit.");
            break;
        };
        log::info!(
            "Re-centred on the image centre: verified {} stars, residual {:.2}\"",
            v.n(),
            v.rms
        );
        (verified, ra_db, dec_db) = (v, ra0, dec0);
    }
    (verified, ra_db, dec_db)
}

/// Solve the WCS for an image against an ASTAP star database.
///
/// Walks a square spiral out from the hint in steps of one field of view, and
/// returns the first position whose quad match survives star-by-star verification.
/// If `params.binning > 1`, `img` is taken to be the binned image and the returned
/// CRPIX/CD/CDELT are scaled back to the unbinned pixel grid.
///
/// All progress is emitted via the `log` crate at INFO level — callers install
/// whichever logger backend they need (file, stderr, both, or none).
///
/// # Errors
///
/// - [`ArcsecError::InvalidParameter`] if `fov` is not positive and finite, or
///   `search_radius` is negative or not finite.
/// - [`ArcsecError::CatalogNotFound`] if `db_path` holds no database called `db_name`.
/// - [`ArcsecError::InsufficientStars`] if fewer than 5 stars are detected.
/// - [`ArcsecError::InsufficientQuads`] if no spiral position yields a verified match.
pub fn solve_image(img: &crate::types::ImageBuffer, params: &SolveParams) -> Result<WcsSolution> {
    // The spiral steps by one FOV out to the search radius, so a zero, negative or
    // NaN FOV would make the step count infinite (and saturate to i32::MAX).
    if !(params.fov.is_finite() && params.fov > 0.0) {
        return Err(ArcsecError::InvalidParameter(format!(
            "field of view must be positive, got {} rad",
            params.fov
        )));
    }
    if !(params.search_radius.is_finite() && params.search_radius >= 0.0) {
        return Err(ArcsecError::InvalidParameter(format!(
            "search radius must be non-negative, got {} rad",
            params.search_radius
        )));
    }

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

    // Detection is not trimmed to a fraction. Stars beyond `-s` are faint enough to
    // be absent from the catalog, which once corrupted 3-NN quads badly enough to
    // justify dropping all but the brightest half; quad redundancy and star-level
    // verification absorb that now, and halving the list halved the quad count.
    // The catalog still reads the full requested depth.
    //
    // It is trimmed to what the database can hold, though (ASTAP's "database
    // limit"). A database is density-limited, so a small field holds only so many
    // catalogue stars, however many the image shows: 973 detections on a 0.21°
    // field (`rnd_080`) against ~100 catalogue stars build quads from stars the
    // catalogue has never heard of, and the field cannot match.
    let star_limit = density_star_limit(params, img);
    let mut stars = stars;
    if stars.len() > star_limit {
        stars.0.sort_by(|a, b| b.snr.total_cmp(&a.snr));
        stars.0.truncate(star_limit);
        log::info!(
            "Database limit for this field is {star_limit} stars; using the {star_limit} brightest."
        );
    }

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

    let oversize: f64 = match params.speed {
        SearchSpeed::Auto if nrstars_image < 35 => 2.0,
        SearchSpeed::Auto if nrstars_image > 140 => 1.0,
        SearchSpeed::Auto => 2.0 * (35.0 / nrstars_image as f64).sqrt(),
        // As ASTAP, never more than one database tile: a larger window could reach
        // past the neighbouring tile, which the tile lookup does not cover.
        SearchSpeed::Slow => {
            let max_fov_deg = match crate::catalog::detect_layout(&params.db_path, &params.db_name)
            {
                CatalogLayout::Areas1476 => 5.142_857_143_f64,
                CatalogLayout::Areas290 => 9.53,
                CatalogLayout::AllSky001 => 180.0,
            };
            2.0_f64.min(max_fov_deg.to_radians() / params.fov).max(1.0)
        }
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
        star_limit,
        nrstars_required,
        oversize,
        min_quads,
        step_size,
        accept: Acceptance::new(nrstars_image, params, img),
        aspect: img.width.max(img.height) as f64 / img.width.min(img.height).max(1) as f64,
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

    if let Some(o) = winner.as_ref().filter(|o| o.refused) {
        log::info!(
            "No solution: the field at search position {} is too distorted for a linear plate.",
            o.idx
        );
        return Err(ArcsecError::InsufficientQuads {
            found: 0,
            required: min_quads,
        });
    }
    if let Some(o) = winner {
        log::info!(
            "{} of {} patterns selected matching within {:.3} tolerance.",
            o.n_matched,
            o.n_raw,
            params.quad_tolerance,
        );

        let v = o.verified;
        let mut wcs = derive_wcs(o.ra_db, o.dec_db, &v.plate, img.width, img.height);
        // The verified pairs, on the original image's pixel grid: a binned pixel
        // centre at 0-based `x` is at `(x + 0.5) * b + 0.5` in unbinned FITS pixels.
        let b = params.binning.max(1) as f64;
        wcs.matched_stars = v
            .img_pos
            .iter()
            .zip(&v.cat_pos)
            .map(|(&(x, y), &(sx, sy))| {
                let (ra, dec) = standard_equatorial(o.ra_db, o.dec_db, sx, sy, 1.0);
                MatchedStar {
                    x: (x + 0.5) * b + 0.5,
                    y: (y + 0.5) * b + 0.5,
                    ra,
                    dec,
                }
            })
            .collect();
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
        wcs.residual_rms = v.rms;
        wcs.stars_matched = v.n();
        wcs.raw_matches = o.n_raw;
        wcs.plate = v.plate;
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

/// Format RA (radians) as `astap_cli` prints it: `"HH: MM  SS.S"`, each field at
/// least two digits (ASTAP's `prepare_ra(ra, ': ')`).
#[must_use]
pub fn format_ra(ra_rad: f64) -> String {
    // Round once, at the printed precision, and only then split into fields.
    // Splitting first and letting `{:.1}` round the seconds printed 59.96 s as
    // "60.0" without carrying into the minutes (and 23:59:59.96 as "23: 59  60.0").
    // ASTAP carries, but prints 23:59:59.96 as "24: 00  00.0"; this wraps to 00h.
    const TENTHS_PER_DAY: f64 = 24.0 * 36_000.0;
    let ra_tenths = ((ra_rad.to_degrees() / 15.0 * 36_000.0)
        .round()
        .rem_euclid(TENTHS_PER_DAY)) as u64;
    let h = ra_tenths / 36_000;
    let m = ra_tenths / 600 % 60;
    let s = ra_tenths % 600 / 10;
    let tenths = ra_tenths % 10;
    format!("{h:02}: {m:02}  {s:02}.{tenths}")
}

/// Format Dec (radians) as `astap_cli` prints it: `"±DDd MM  SS"`, each field at
/// least two digits (ASTAP's `prepare_dec(dec, 'd ')`).
#[must_use]
pub fn format_dec(dec_rad: f64) -> String {
    let dec_deg = dec_rad.to_degrees();
    let sign = if dec_deg < 0.0 { '-' } else { '+' };
    let dec_secs = (dec_deg.abs() * 3600.0).round() as u64;
    let dd = dec_secs / 3600;
    let dm = dec_secs / 60 % 60;
    let ds = dec_secs % 60;
    format!("{sign}{dd:02}d {dm:02}  {ds:02}")
}

/// Format RA and Dec (radians) as `astap_cli`'s `Solution found:` line does:
/// `"HH: MM  SS.S ±DDd MM  SS"`. (Its `Start position:` line puts a comma between
/// the two; see [`format_ra`] and [`format_dec`].)
#[must_use]
pub fn format_radec(ra_rad: f64, dec_rad: f64) -> String {
    format!("{} {}", format_ra(ra_rad), format_dec(dec_rad))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::coords::{ang_sep, standard_equatorial};
    use crate::test_support::{
        Rng, SkySpec, SkyStar, TempDir, TruthWcs, random_sky, render, write_001_db, write_290_db,
        write_1476_db,
    };
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
    fn format_radec_carries_rounded_seconds() {
        // 1h 59m 59.97s must round up to 2h 00m 00.0s, not print "60.0" seconds.
        let ra = deg((1.0 + 59.0 / 60.0 + 59.97 / 3600.0) * 15.0);
        // +10° 59' 59.7" rounds to +11° 00' 00".
        let dec = deg(10.0 + 59.0 / 60.0 + 59.7 / 3600.0);
        assert_eq!(format_radec(ra, dec), "02: 00  00.0 +11d 00  00");
        // RA just short of 24h wraps to 0h.
        let s = format_radec(deg(359.999_999_9), deg(-0.5));
        assert_eq!(s, "00: 00  00.0 -00d 30  00");
        // An ordinary value is unchanged by the rewrite.
        assert_eq!(
            format_radec(deg((5.0 + 35.0 / 60.0 + 17.3 / 3600.0) * 15.0), deg(-5.39)),
            "05: 35  17.3 -05d 23  24"
        );
    }

    /// Byte-for-byte what `astap_cli` prints (checked against 2026.07.30):
    /// `Start position: 04: 20  00.0, +35d 00  00` and
    /// `Solution found: 04: 20  00.0 +35d 00  00`. Every field is at least two
    /// digits wide, as ASTAP's `LeadingZero` makes it.
    #[test]
    fn ra_and_dec_are_formatted_as_astap_cli_prints_them() {
        let ra = deg(65.0); // 4h 20m
        let dec = deg(35.0);
        assert_eq!(format_ra(ra), "04: 20  00.0");
        assert_eq!(format_dec(dec), "+35d 00  00");
        assert_eq!(format_radec(ra, dec), "04: 20  00.0 +35d 00  00");
        assert_eq!(
            format_radec(
                deg((13.0 + 7.0 / 60.0 + 9.25 / 3600.0) * 15.0),
                -deg(89.0 + 1.0 / 60.0 + 2.0 / 3600.0)
            ),
            "13: 07  09.3 -89d 01  02"
        );
        assert_eq!(format_dec(deg(-0.0001)), "-00d 00  00");
    }

    #[test]
    fn solve_image_rejects_a_non_positive_fov() {
        let img = ImageBuffer::new(64, 64);
        let params = SolveParams {
            ra_hint: 0.0,
            dec_hint: 0.0,
            fov: 0.0,
            search_radius: 0.1,
            quad_tolerance: 0.007,
            hfd_min: 1.5,
            max_stars: 500,
            db_path: std::path::PathBuf::from("/nonexistent"),
            db_name: "d50".into(),
            binning: 1,
            method: SolveMethod::Quads,
            threads: 1,
            speed: SearchSpeed::Auto,
        };
        assert!(matches!(
            solve_image(&img, &params),
            Err(ArcsecError::InvalidParameter(_))
        ));
    }

    // ── Plate-fit helpers ─────────────────────────────────────────────────────

    /// A known similarity transform (pixels → catalogue arcsec), with a flip.
    fn known_plate() -> PlateConstants {
        let (s, r) = (3.2_f64, 0.61_f64);
        PlateConstants {
            a: -s * r.cos(),
            b: s * r.sin(),
            c: 640.0,
            d: s * r.sin(),
            e: s * r.cos(),
            f: -512.0,
        }
    }

    fn apply(p: &PlateConstants, (x, y): (f64, f64)) -> (f64, f64) {
        (p.a * x + p.b * y + p.c, p.d * x + p.e * y + p.f)
    }

    fn plate_close(p: &PlateConstants, q: &PlateConstants, tol: f64) -> bool {
        [
            (p.a, q.a),
            (p.b, q.b),
            (p.c, q.c),
            (p.d, q.d),
            (p.e, q.e),
            (p.f, q.f),
        ]
        .iter()
        .all(|(u, v)| (u - v).abs() <= tol)
    }

    /// The acceptance rule for a well-populated image, for the plate of
    /// [`known_plate`].
    const STRICT: Acceptance = Acceptance {
        min_stars: MIN_VERIFIED_STARS,
        expected_scale: 3.2,
    };

    fn star_at(x: f64, y: f64) -> Star {
        Star {
            x,
            y,
            snr: 50.0,
            hfd: 2.5,
        }
    }

    /// 40 exact pairs under `known_plate`, then five pairs whose catalogue side is
    /// displaced by `outlier(k)`.
    fn pairs_with_outliers(outlier: impl Fn(usize, (f64, f64)) -> (f64, f64)) -> PairedPositions {
        let plate = known_plate();
        let mut rng = Rng::new(7);
        let mut img = Vec::new();
        let mut cat = Vec::new();
        for _ in 0..40 {
            let p = (rng.range(0.0, 500.0), rng.range(0.0, 500.0));
            img.push(p);
            cat.push(apply(&plate, p));
        }
        for k in 0..5 {
            let p = (rng.range(0.0, 500.0), rng.range(0.0, 500.0));
            img.push(p);
            cat.push(outlier(k, apply(&plate, p)));
        }
        (img, cat)
    }

    #[test]
    fn sigma_clip_pairs_rejects_outliers_and_keeps_the_rest() {
        // Five wrong pairings, each ~100" (30 px) from where the plate puts them.
        let (img, cat) = pairs_with_outliers(|k, (x, y)| {
            let a = k as f64 * 1.3;
            (x + 100.0 * a.cos(), y + 100.0 * a.sin())
        });
        let (ci, cc) = sigma_clip_pairs(img, cat, 3.0, 3);
        assert_eq!(ci.len(), 40, "all and only the true pairs survive");
        let fit = solve_plate_constants(&ci, &cc).unwrap();
        assert!(plate_close(&fit, &known_plate(), 1e-6), "{fit:?}");
    }

    /// `sigma_clip_pairs` gives up as soon as a fit fails, and the first fit is made
    /// on the contaminated set. Five gross outliers in 45 pairs are enough to skew
    /// that fit past the similarity check in `solve_plate_constants`
    /// (`BadSolution`), so nothing is clipped and all 45 come
    /// back. `try_position` would then refit the same contaminated set, fail the
    /// same check, and abandon a position whose 40 good pairs would have solved it.
    /// So the clipper's first fit must be unchecked.
    #[test]
    fn sigma_clip_pairs_rejects_gross_outliers() {
        let (img, cat) =
            pairs_with_outliers(|k, _| (1000.0 + 150.0 * k as f64, -900.0 + 70.0 * k as f64));
        assert!(matches!(
            solve_plate_constants(&img, &cat),
            Err(ArcsecError::BadSolution { .. })
        ));
        let (ci, _) = sigma_clip_pairs(img, cat, 3.0, 3);
        assert_eq!(ci.len(), 40, "the five gross outliers should be clipped");
    }

    /// The pattern-pair fit of both methods: a set whose plain fit is refused
    /// (gross outliers skew it off a similarity) still yields the true plate, from
    /// the 40 good pairs; one with too few good pairs left yields nothing.
    #[test]
    fn fit_pattern_pairs_recovers_a_plate_the_plain_fit_refuses() {
        let (img, cat) =
            pairs_with_outliers(|k, _| (1000.0 + 150.0 * k as f64, -900.0 + 70.0 * k as f64));
        assert!(solve_plate_constants(&img, &cat).is_err());
        let (plate, n) = fit_pattern_pairs(img.clone(), cat.clone(), 3).expect("clipped fit");
        assert_eq!(n, 40);
        assert!(plate_close(&plate, &known_plate(), 1e-6), "{plate:?}");
        // Clean pairs: nothing to clip, the same plate as a plain fit.
        let (plate, n) = fit_pattern_pairs(img[..40].to_vec(), cat[..40].to_vec(), 3).unwrap();
        assert_eq!(n, 40);
        assert!(plate_close(&plate, &known_plate(), 1e-6));
        // More pairs demanded than survive the clipping.
        assert!(fit_pattern_pairs(img, cat, 41).is_none());
    }

    #[test]
    fn sigma_clip_pairs_leaves_too_few_pairs_alone() {
        let img = vec![(0.0, 0.0), (1.0, 0.0)];
        let cat = vec![(5.0, 5.0), (9.0, 9.0)];
        let (ci, cc) = sigma_clip_pairs(img.clone(), cat.clone(), 3.0, 3);
        assert_eq!((ci, cc), (img, cat));
    }

    #[test]
    fn verify_and_refit_recovers_the_plate_from_a_rough_guess() {
        let truth = known_plate();
        let mut rng = Rng::new(11);
        let mut img_stars = Vec::new();
        let mut cat_stars = Vec::new();
        for _ in 0..60 {
            let (x, y) = (rng.range(5.0, 395.0), rng.range(5.0, 295.0));
            img_stars.push(star_at(x, y));
            let (cx, cy) = apply(&truth, (x, y));
            cat_stars.push(star_at(cx, cy));
        }
        // Catalogue stars that fall outside the frame must be ignored, not paired.
        for k in 0..20 {
            let (cx, cy) = apply(&truth, (-300.0 - 10.0 * k as f64, 900.0));
            cat_stars.push(star_at(cx, cy));
        }
        // Start 2 px and a little rotation away from the truth.
        let mut rough = truth.clone();
        rough.c += 2.0 * truth.a;
        rough.f += 2.0 * truth.e;
        rough.b += 0.01;
        let v = verify_and_refit(
            &StarList(img_stars),
            &StarList(cat_stars),
            &rough,
            400,
            300,
            &STRICT,
        )
        .expect("a correct plate must verify");
        assert_eq!(v.n(), 60);
        assert_eq!(v.cat_pos.len(), 60);
        assert!(v.rms < 1e-6, "rms {}", v.rms);
        assert!(plate_close(&v.plate, &truth, 1e-6), "{:?}", v.plate);
        // Each pair is a star and its own catalogue entry.
        for (&(x, y), &(cx, cy)) in v.img_pos.iter().zip(&v.cat_pos) {
            let (px, py) = apply(&truth, (x, y));
            assert!((px - cx).hypot(py - cy) < 1e-6);
        }
    }

    #[test]
    fn verify_and_refit_rejects_too_few_or_clustered_matches() {
        let truth = known_plate();
        let mut rng = Rng::new(12);
        let build = |pts: &[(f64, f64)]| {
            let img = StarList(pts.iter().map(|&(x, y)| star_at(x, y)).collect());
            let cat = StarList(
                pts.iter()
                    .map(|&p| apply(&truth, p))
                    .map(|(x, y)| star_at(x, y))
                    .collect(),
            );
            (img, cat)
        };

        // 20 well-spread stars: fewer than MIN_VERIFIED_STARS.
        let few: Vec<_> = (0..20)
            .map(|_| (rng.range(0.0, 400.0), rng.range(0.0, 300.0)))
            .collect();
        let (img, cat) = build(&few);
        assert!(verify_and_refit(&img, &cat, &truth, 400, 300, &STRICT).is_none());

        // 80 stars, all in one 40-pixel corner: rotation is unconstrained.
        let clustered: Vec<_> = (0..80)
            .map(|_| (rng.range(0.0, 40.0), rng.range(0.0, 40.0)))
            .collect();
        let (img, cat) = build(&clustered);
        assert!(verify_and_refit(&img, &cat, &truth, 400, 300, &STRICT).is_none());

        // The same 80 spread over the frame pass.
        let spread: Vec<_> = (0..80)
            .map(|_| (rng.range(0.0, 400.0), rng.range(0.0, 300.0)))
            .collect();
        let (img, cat) = build(&spread);
        assert!(verify_and_refit(&img, &cat, &truth, 400, 300, &STRICT).is_some());

        // Degenerate inputs.
        let empty = StarList::default();
        assert!(verify_and_refit(&empty, &cat, &truth, 400, 300, &STRICT).is_none());
        let mut singular = truth.clone();
        singular.a = 0.0;
        singular.b = 0.0;
        assert!(verify_and_refit(&img, &cat, &singular, 400, 300, &STRICT).is_none());
    }

    // ── End-to-end solves against synthetic catalogues ────────────────────────

    #[derive(Clone, Copy)]
    enum Db {
        Areas1476,
        Areas290,
        AllSky001,
    }

    /// A rendered field and the database it was drawn from.
    struct Scene {
        dir: TempDir,
        img: ImageBuffer,
        truth: TruthWcs,
        /// Every star drawn, in and around the frame.
        sky: Vec<SkyStar>,
    }

    /// Render ~`n_in_frame` stars through `truth` and write the surrounding sky
    /// (six fields wide, so offset hints still find their stars) as a database.
    fn scene(truth: TruthWcs, db: Db, n_in_frame: usize, seed: u64) -> Scene {
        let mut rng = Rng::new(seed);
        let scale_deg = truth.cd[1].hypot(truth.cd[3]);
        let (w_deg, h_deg) = (
            truth.width as f64 * scale_deg,
            truth.height as f64 * scale_deg,
        );
        let side = 6.0 * w_deg.max(h_deg);
        let sky = random_sky(
            &mut rng,
            &SkySpec {
                ra0: truth.ra0,
                dec0: truth.dec0,
                side_deg: side,
                n: (n_in_frame as f64 * side * side / (w_deg * h_deg)) as usize,
                min_sep_deg: 12.0 * scale_deg,
                mag_lo: 10.0,
                mag_hi: 14.5,
            },
        );
        // A PSF of ~1.3 px on a 5"/px frame, scaled so binned frames stay sampled.
        let sigma = 1.3 * 5.0 / (scale_deg * 3600.0);
        let img = render(
            &truth,
            &sky,
            sigma.max(1.3),
            1000.0,
            8.0,
            30_000.0,
            &mut rng,
        );
        let dir = TempDir::new("solve");
        match db {
            Db::Areas1476 => write_1476_db(dir.path(), "t50", &sky),
            Db::Areas290 => write_290_db(dir.path(), "t50", &sky),
            Db::AllSky001 => write_001_db(dir.path(), "t50", &sky),
        }
        Scene {
            dir,
            img,
            truth,
            sky,
        }
    }

    /// Parameters that fail every check that needs a database.
    fn params_for_blank() -> SolveParams {
        SolveParams {
            ra_hint: 0.0,
            dec_hint: 0.0,
            fov: deg(1.0),
            search_radius: 0.0,
            quad_tolerance: 0.007,
            hfd_min: 1.5,
            max_stars: 500,
            db_path: std::path::PathBuf::from("/nonexistent"),
            db_name: "d50".into(),
            binning: 1,
            method: SolveMethod::Quads,
            threads: 1,
            speed: SearchSpeed::Auto,
        }
    }

    fn params_for(s: &Scene, ra_hint: f64, dec_hint: f64) -> SolveParams {
        SolveParams {
            ra_hint,
            dec_hint,
            fov: (s.truth.height as f64 * s.truth.cd[1].hypot(s.truth.cd[3])).to_radians(),
            search_radius: deg(2.0),
            quad_tolerance: 0.007,
            hfd_min: 1.5,
            max_stars: 500,
            db_path: s.dir.path().to_path_buf(),
            db_name: "t50".into(),
            binning: 1,
            method: SolveMethod::Quads,
            threads: 1,
            speed: SearchSpeed::Auto,
        }
    }

    fn assert_solved(s: &Scene, wcs: &WcsSolution, tol_arcsec: f64) {
        let err = s.truth.max_error_arcsec(wcs);
        assert!(
            err < tol_arcsec,
            "worst centre/corner error {err:.3}\" (matched {}, rms {:.3})",
            wcs.stars_matched,
            wcs.residual_rms
        );
        assert!(wcs.stars_matched >= 10);
        // Star-level residual under a third of a pixel.
        let scale_arcsec = s.truth.cd[1].hypot(s.truth.cd[3]) * 3600.0;
        assert!(
            wcs.residual_rms < 0.3 * scale_arcsec,
            "rms {}",
            wcs.residual_rms
        );
        assert!(wcs.raw_matches > 0);
        assert_matches_agree(wcs, 0.3, 1.0);
        assert!(wcs.mag_limit > 10.0 && wcs.mag_limit <= 14.5);
        assert!(
            wcs.cdelt1 < 0.0 && wcs.cdelt2 > 0.0,
            "CDELT sign convention"
        );
    }

    /// The verified pairs agree with the solution: RMS under `rms_px` pixels, and
    /// none further off than the final verification radius (`binning` pixels each).
    fn assert_matches_agree(wcs: &WcsSolution, rms_px: f64, binning: f64) {
        assert_eq!(wcs.matched_stars.len(), wcs.stars_matched);
        assert!(wcs.sip.is_none(), "solve_image never fits SIP");
        let tan = crate::wcs::TanWcs::from(wcs);
        let mut sq = 0.0;
        for m in &wcs.matched_stars {
            let (x, y) = tan.sky_to_pixel(m.ra, m.dec).unwrap();
            let d = (x - m.x).hypot(y - m.y);
            assert!(
                d < VERIFY_RADII[VERIFY_RADII.len() - 1] * binning,
                "pair at ({:.2},{:.2}) projects to ({x:.2},{y:.2})",
                m.x,
                m.y
            );
            sq += d * d;
        }
        let rms = (sq / wcs.matched_stars.len() as f64).sqrt();
        assert!(rms < rms_px, "pair rms {rms} px");
    }

    #[test]
    fn solves_a_1476_database_from_an_offset_hint() {
        let truth = TruthWcs::new(deg(84.3), deg(-5.2), 5.0, 23.0, false, 400, 320);
        let s = scene(truth, Db::Areas1476, 130, 1);
        // Hint roughly one field away in each axis: the spiral has to move.
        let mut p = params_for(&s, deg(84.3 + 0.6), deg(-5.2 - 0.45));
        p.threads = 4;
        let wcs = solve_image(&s.img, &p).expect("solve");
        assert_solved(&s, &wcs, 1.0);
        assert!(wcs.search_dist_deg > 0.1, "solved at the hint itself?");
        assert!(wcs.step_distances.len() > 1);
        // The pixel scale and rotation come back too.
        assert!((wcs.cdelt2 * 3600.0 - 5.0).abs() < 0.01, "{}", wcs.cdelt2);
        assert!((wcs.crota2 - 23.0).abs() < 0.05, "crota2 {}", wcs.crota2);
    }

    #[test]
    fn solves_a_mirrored_image_on_a_290_database() {
        let truth = TruthWcs::new(deg(201.0), deg(47.5), 6.0, 160.0, true, 360, 360);
        let s = scene(truth, Db::Areas290, 120, 2);
        let wcs = solve_image(&s.img, &params_for(&s, truth.ra0, truth.dec0)).expect("solve");
        assert_solved(&s, &wcs, 1.0);
        assert!(wcs.search_dist_deg < 1e-9, "should solve at the hint");
        // A mirrored image has det(CD) > 0.
        assert!(wcs.cd1_1 * wcs.cd2_2 - wcs.cd1_2 * wcs.cd2_1 > 0.0);
    }

    #[test]
    fn solves_across_ra_zero_with_an_all_sky_001_database() {
        // The field straddles RA 0h, so its catalogue stars sit either side of 2π.
        let truth = TruthWcs::new(deg(0.05), deg(21.0), 5.0, -70.0, false, 360, 300);
        let s = scene(truth, Db::AllSky001, 120, 3);
        let wcs = solve_image(&s.img, &params_for(&s, truth.ra0, truth.dec0)).expect("solve");
        assert_solved(&s, &wcs, 1.0);
    }

    #[test]
    fn solves_across_ra_zero_with_a_1476_database() {
        let truth = TruthWcs::new(deg(359.97), deg(-33.0), 5.0, 95.0, false, 360, 300);
        let s = scene(truth, Db::Areas1476, 120, 4);
        let wcs = solve_image(&s.img, &params_for(&s, truth.ra0, truth.dec0)).expect("solve");
        assert_solved(&s, &wcs, 1.0);
    }

    #[test]
    fn solves_a_field_near_the_celestial_pole() {
        let truth = TruthWcs::new(deg(40.0), deg(88.9), 5.0, 10.0, false, 360, 300);
        let s = scene(truth, Db::Areas1476, 120, 5);
        let wcs = solve_image(&s.img, &params_for(&s, truth.ra0, truth.dec0)).expect("solve");
        assert_solved(&s, &wcs, 1.0);
    }

    /// The plate constants are fitted in the tangent plane of the spiral position
    /// that matched (`ra_db`, `dec_db`), but `derive_wcs` then moves CRVAL to the
    /// image centre and keeps the CD matrix unchanged, as though the two tangent
    /// planes were the same. They are not, and the error grows linearly with the
    /// distance between the matched spiral position and the true field centre.
    ///
    /// Measured on this 1.5° × 1.25° field at 15"/px: worst-corner error 0.35" with
    /// the hint on the centre, 7.5" at 0.2° off, 14.5" at 0.4°, 21" at 0.6° (1.4 px),
    /// while the star-level RMS stays at 0.7-0.85" throughout — the verification
    /// cannot see it, because it runs in the same (offset) tangent plane. Spiral
    /// positions land up to half a step (half a field) from the truth, and a blind
    /// estimate can be a whole field off, so this is well inside normal use. 5" is
    /// the corner error `scripts/benchmark.py` counts as a false positive.
    #[test]
    fn accuracy_does_not_depend_on_the_hint_offset() {
        let truth = TruthWcs::new(deg(150.0), deg(30.0), 15.0, 20.0, false, 360, 300);
        let s = scene(truth, Db::Areas1476, 120, 21);
        let off = 0.4;
        let p = params_for(&s, deg(150.0 + off / deg(30.0).cos()), deg(30.0 + off));
        let wcs = solve_image(&s.img, &p).expect("solve");
        assert!(wcs.search_dist_deg < 1e-9, "solved at the hint");
        let err = s.truth.max_error_arcsec(&wcs);
        assert!(
            err < 5.0,
            "worst corner error {err:.2}\" with a {off}° hint offset"
        );
    }

    #[test]
    fn solves_with_the_tetra_method() {
        let truth = TruthWcs::new(deg(150.0), deg(2.0), 5.0, 45.0, false, 360, 300);
        let s = scene(truth, Db::Areas1476, 110, 6);
        let mut p = params_for(&s, truth.ra0, truth.dec0);
        p.method = SolveMethod::Tetra;
        let wcs = solve_image(&s.img, &p).expect("solve");
        assert_solved(&s, &wcs, 1.0);
    }

    #[test]
    fn slow_speed_solves_from_an_offset_hint() {
        let truth = TruthWcs::new(deg(84.3), deg(-5.2), 5.0, 23.0, false, 400, 320);
        let s = scene(truth, Db::Areas1476, 130, 1);
        let mut p = params_for(&s, deg(84.3 + 0.6), deg(-5.2 - 0.45));
        p.speed = SearchSpeed::Slow;
        let wcs = solve_image(&s.img, &p).expect("solve");
        assert_solved(&s, &wcs, 1.0);
    }

    #[test]
    fn binned_solve_is_reported_on_the_unbinned_pixel_grid() {
        // Render at full resolution, then solve the 2×2-binned frame.
        let truth = TruthWcs::new(deg(10.0), deg(40.0), 2.5, 30.0, false, 720, 600);
        let s = scene(truth, Db::Areas1476, 120, 7);
        let binned = s.img.bin_image(2);
        assert_eq!((binned.width, binned.height), (360, 300));
        let mut p = params_for(&s, truth.ra0, truth.dec0);
        p.binning = 2;
        let wcs = solve_image(&binned, &p).expect("solve");
        // crpix is the centre of the unbinned frame, and the scale is unbinned.
        assert!((wcs.crpix1 - 360.5).abs() < 1e-9, "crpix1 {}", wcs.crpix1);
        assert!((wcs.crpix2 - 300.5).abs() < 1e-9, "crpix2 {}", wcs.crpix2);
        assert!((wcs.cdelt2 * 3600.0 - 2.5).abs() < 0.01, "{}", wcs.cdelt2);
        let err = s.truth.max_error_arcsec(&wcs);
        assert!(err < 2.0, "worst corner error {err:.3}\"");
        // The pairs are on the unbinned grid too: one binned pixel is two of these.
        assert_matches_agree(&wcs, 0.6, 2.0);
    }

    #[test]
    fn the_star_limit_is_the_database_density_times_the_field_area() {
        let params = |fov_deg: f64, db: &str, max_stars: usize| SolveParams {
            fov: deg(fov_deg),
            max_stars,
            db_name: db.into(),
            ..params_for_blank()
        };
        let square = ImageBuffer::new(200, 200);
        let wide = ImageBuffer::new(400, 200);
        // d80 on a 0.2° square field: 8000 × 0.04 = 320 stars.
        assert_eq!(density_star_limit(&params(0.2, "d80", 500), &square), 320);
        // The same long side on a 2:1 frame covers half the area.
        assert_eq!(density_star_limit(&params(0.2, "d80", 500), &wide), 160);
        // Never more than -s.
        assert_eq!(density_star_limit(&params(1.0, "d80", 500), &square), 500);
        assert_eq!(density_star_limit(&params(0.2, "d80", 100), &square), 100);
        // g05 (500/deg²) binds only below ~1°, w08 (1/deg²) on all but the widest.
        assert_eq!(density_star_limit(&params(0.8, "g05", 500), &square), 320);
        assert_eq!(density_star_limit(&params(20.0, "w08", 500), &wide), 200);
        // Unknown density: -s.
        assert_eq!(density_star_limit(&params(0.1, "v17", 500), &square), 500);
    }

    /// A field showing far more stars than the database holds there: with every
    /// detection the image quads are built from stars the catalogue does not have,
    /// and it does not solve; capped at the database's density it does.
    #[test]
    fn a_frame_deeper_than_the_database_solves_at_the_database_limit() {
        // 0.5° x 0.42° at 3"/px: 0.21 deg², ~450 stars in the frame.
        let truth = TruthWcs::new(deg(250.0), deg(36.0), 3.0, 12.0, false, 600, 500);
        let s = scene(truth, Db::Areas1476, 450, 31);
        // The database holds only the brightest 200 per square degree (the
        // scene's sky is 3° on a side), ~42 of them in the frame.
        let mut sky = s.sky.clone();
        sky.sort_by(|a, b| a.mag.total_cmp(&b.mag));
        sky.truncate(200 * 9);
        write_1476_db(s.dir.path(), "t02", &sky);
        // The same stars under a name whose density is unknown, so no limit.
        write_1476_db(s.dir.path(), "t17", &sky);

        let mut p = params_for(&s, truth.ra0, truth.dec0);
        p.fov = (600.0 * 3.0 / 3600.0_f64).to_radians();
        p.search_radius = 0.0;
        p.db_name = "t17".into();
        assert!(
            solve_image(&s.img, &p).is_err(),
            "every detection: should not match"
        );
        p.db_name = "t02".into();
        let wcs = solve_image(&s.img, &p).expect("solve at the database limit");
        assert!(s.truth.max_error_arcsec(&wcs) < 2.0);
    }

    #[test]
    fn min_verified_stars_relaxes_only_for_sparse_images() {
        for (n, want) in [
            (0, 10),
            (5, 10),
            (66, 10),
            (67, 11),
            (100, 15),
            (193, 29),
            (194, 30),
            (200, 30),
            (500, 30),
            (usize::MAX, 30),
        ] {
            assert_eq!(min_verified_stars(n), want, "{n} detections");
        }
    }

    /// Below 30 matches a solution must be at the hint's scale, to half a pixel.
    #[test]
    fn a_sparse_match_must_have_the_expected_scale_and_a_tight_fit() {
        let truth = known_plate(); // 3.2"/px
        let verified = |n: usize, rms_px: f64, scale: f64| {
            let mut plate = truth.clone();
            for c in [&mut plate.a, &mut plate.b, &mut plate.d, &mut plate.e] {
                *c *= scale;
            }
            Verified {
                plate,
                rms: rms_px * 3.2 * scale,
                img_pos: vec![(0.0, 0.0); n],
                cat_pos: vec![(0.0, 0.0); n],
            }
        };
        let accept = Acceptance {
            min_stars: 12,
            expected_scale: 3.2,
        };
        // Enough stars: scale and residual are not looked at.
        assert!(accept.accepts(&verified(30, 3.9, 1.36), 0.5));
        // Sparse, right scale, tight fit.
        assert!(accept.accepts(&verified(12, 0.3, 1.0), 0.5));
        assert!(accept.accepts(&verified(20, 0.49, 1.09), 0.5));
        assert!(accept.accepts(&verified(20, 0.49, 0.91), 0.5));
        // Sparse and wrong: scale 10% off, a loose fit, too few, too clustered.
        assert!(!accept.accepts(&verified(20, 0.3, 1.11), 0.5));
        assert!(!accept.accepts(&verified(20, 0.3, 0.89), 0.5));
        assert!(!accept.accepts(&verified(29, 0.51, 1.0), 0.5));
        assert!(!accept.accepts(&verified(11, 0.1, 1.0), 0.5));
        assert!(!accept.accepts(&verified(20, 0.1, 1.0), 0.1));
        // The false positive the relaxed count alone lets through (ls2_25).
        assert!(!accept.accepts(&verified(12, 2.9, 1.36), 0.5));
    }

    /// A frame showing only its ~20 brightest stars against a deep catalogue: it
    /// solves with fewer than 30 matches when the hint's scale is right, and is
    /// refused when the scale the hint implies is 20% off.
    #[test]
    fn a_sparse_frame_solves_at_the_hint_scale_only() {
        let truth = TruthWcs::new(deg(30.0), deg(-12.0), 5.0, 40.0, false, 360, 300);
        let s = scene(truth, Db::Areas1476, 150, 41);
        let mut bright = s.sky.clone();
        bright.sort_by(|a, b| a.mag.total_cmp(&b.mag));
        let bright: Vec<SkyStar> = bright
            .into_iter()
            .filter(|st| {
                s.truth
                    .sky_to_pixel(st.ra, st.dec)
                    .is_some_and(|(x, y)| (5.0..355.0).contains(&x) && (5.0..295.0).contains(&y))
            })
            .take(22)
            .collect();
        let mut rng = Rng::new(42);
        let img = render(&s.truth, &bright, 1.3, 1000.0, 8.0, 30_000.0, &mut rng);
        let mut p = params_for(&s, truth.ra0, truth.dec0);
        p.fov = (360.0 * 5.0 / 3600.0_f64).to_radians(); // the long side
        p.search_radius = 0.0;
        let wcs = solve_image(&img, &p).expect("sparse solve");
        assert!(
            wcs.stars_matched < MIN_VERIFIED_STARS,
            "{}",
            wcs.stars_matched
        );
        assert!(s.truth.max_error_arcsec(&wcs) < 2.0);

        p.fov *= 1.2;
        assert!(matches!(
            solve_image(&img, &p),
            Err(ArcsecError::InsufficientQuads { .. })
        ));
    }

    /// A shallow frame, ~40 stars, against a catalogue twelve times deeper: the
    /// image's quads join neighbours the deep catalogue's do not, and only the
    /// density-matched catalogue quads match them (without them this does not
    /// solve).
    #[test]
    fn a_shallow_frame_matches_the_density_matched_catalogue_quads() {
        let truth = TruthWcs::new(deg(140.0), deg(55.0), 5.0, -25.0, true, 360, 300);
        let s = scene(truth, Db::Areas1476, 500, 51);
        let mut bright = s.sky.clone();
        bright.sort_by(|a, b| a.mag.total_cmp(&b.mag));
        let in_frame = |st: &SkyStar| {
            s.truth
                .sky_to_pixel(st.ra, st.dec)
                .is_some_and(|(x, y)| (0.0..360.0).contains(&x) && (0.0..300.0).contains(&y))
        };
        let n_frame = bright.iter().filter(|st| in_frame(st)).count();
        bright.truncate(bright.len() * 40 / n_frame.max(1));
        let mut rng = Rng::new(52);
        let img = render(&s.truth, &bright, 1.3, 1000.0, 8.0, 30_000.0, &mut rng);
        let mut p = params_for(&s, truth.ra0, truth.dec0);
        p.fov = (360.0 * 5.0 / 3600.0_f64).to_radians();
        p.search_radius = 0.0;
        let wcs = solve_image(&img, &p).expect("shallow solve");
        assert!(s.truth.max_error_arcsec(&wcs) < 2.0);
    }

    /// A 1024 × 768 field at 10"/px with `corner_px` of radial distortion at the
    /// corners (positive: pincushion).
    fn distorted_scene(corner_px: f64, seed: u64) -> Scene {
        let truth = TruthWcs::new(deg(84.3), deg(-5.2), 10.0, 23.0, false, 1024, 768)
            .with_corner_distortion(corner_px);
        scene(truth, Db::Areas1476, 300, seed)
    }

    /// Worst error (arcsec) of a solution with its SIP terms against the truth,
    /// over the centre, the corners and two edge midpoints.
    fn sip_error_arcsec(s: &Scene, wcs: &WcsSolution) -> f64 {
        let tan = crate::wcs::TanWcs::from(wcs);
        let (w, h) = (s.truth.width as f64 - 1.0, s.truth.height as f64 - 1.0);
        let mut worst: f64 = 0.0;
        for (fx, fy) in [
            (0.5, 0.5),
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 1.0),
            (1.0, 1.0),
            (0.5, 0.0),
            (0.0, 0.5),
        ] {
            let (x, y) = (w * fx, h * fy);
            let (ra_t, dec_t) = s.truth.pixel_to_sky(x, y);
            let (ra_s, dec_s) = tan.pixel_to_sky(x + 1.0, y + 1.0);
            let sep = crate::test_support::separation(ra_t, dec_t, ra_s, dec_s);
            worst = worst.max(sep.to_degrees() * 3600.0);
        }
        worst
    }

    #[test]
    fn a_distorted_field_reports_the_best_linear_plate_over_the_frame() {
        // 30 px of pincushion at the corners. The 2 px verification keeps only the
        // stars a linear plate fits, around the centre, and a plate fitted to them
        // alone is ~190" out at the corners; the best linear plate over the frame
        // is ~103" out, and that is what should be reported.
        for (hint_ra, hint_dec) in [(84.3, -5.2), (84.3 + 0.9, -5.2 - 0.7)] {
            let s = distorted_scene(30.0, 7);
            let wcs =
                solve_image(&s.img, &params_for(&s, deg(hint_ra), deg(hint_dec))).expect("solve");
            let floor = s.truth.linear_floor_arcsec();
            let err = s.truth.max_error_arcsec(&wcs);
            assert!(floor > 80.0, "floor {floor:.1}\"");
            assert!(
                err < floor + 5.0,
                "corner error {err:.1}\" against a linear floor of {floor:.1}\""
            );
            assert!(wcs.sip.is_none(), "solve_image never fits SIP");
            // The model's pairs reach the corners, so --sip can follow the distortion.
            assert!(wcs.stars_matched > 150, "{} stars", wcs.stars_matched);
            let mut with_sip = wcs.clone();
            with_sip.sip = crate::wcs::fit_sip(&wcs, 1024, 768);
            assert!(with_sip.sip.is_some(), "the distortion is significant");
            let sip_err = sip_error_arcsec(&s, &with_sip);
            assert!(sip_err < 3.0, "SIP error {sip_err:.2}\"");
        }
    }

    /// [`distorted_scene`] with the right third of the frame blanked to the
    /// background, as a nebula or a dark cloud would leave it.
    fn part_empty_scene(corner_px: f64) -> Scene {
        let mut s = distorted_scene(corner_px, 7);
        let w = s.img.width;
        for y in 0..s.img.height {
            for x in (2 * w / 3)..w {
                s.img.data[y * w + x] = 1000.0;
            }
        }
        s
    }

    #[test]
    fn strong_distortion_that_cannot_be_modelled_over_the_frame_is_refused() {
        // 30 px of barrel distortion, and stars in only two thirds of the frame: a
        // cubic cannot be trusted over the empty third, and the linear plate the
        // verified stars give is ~230" out at the corners against a floor of
        // ~124". Reporting it would be a false positive.
        let s = part_empty_scene(30.0);
        let r = solve_image(&s.img, &params_for(&s, deg(84.3), deg(-5.2)));
        assert!(
            matches!(r, Err(ArcsecError::InsufficientQuads { .. })),
            "{:?}",
            r.map(|w| s.truth.max_error_arcsec(&w))
        );
    }

    #[test]
    fn an_undistorted_field_with_an_empty_third_still_solves() {
        let s = part_empty_scene(0.0);
        let wcs = solve_image(&s.img, &params_for(&s, deg(84.3), deg(-5.2))).expect("solve");
        assert_solved(&s, &wcs, 1.0);
    }

    #[test]
    fn mild_distortion_is_modelled_too() {
        // 3 px at the corners: a linear plate fitted to the verified stars is
        // already close to the floor, but not as close as the best one.
        let s = distorted_scene(3.0, 11);
        let wcs = solve_image(&s.img, &params_for(&s, deg(84.3), deg(-5.2))).expect("solve");
        let floor = s.truth.linear_floor_arcsec();
        let err = s.truth.max_error_arcsec(&wcs);
        assert!(
            err < floor + 2.0,
            "corner error {err:.1}\" against a floor of {floor:.1}\""
        );
    }

    #[test]
    fn a_field_absent_from_the_catalogue_does_not_solve() {
        // The image shows one random sky, the database holds a different one at the
        // same place: nothing may verify, however many quads happen to match.
        let truth = TruthWcs::new(deg(120.0), deg(-40.0), 5.0, 0.0, false, 360, 300);
        let s = scene(truth, Db::Areas1476, 120, 8);
        let decoy = TempDir::new("decoy");
        let mut rng = Rng::new(99);
        let other = random_sky(
            &mut rng,
            &SkySpec {
                ra0: truth.ra0,
                dec0: truth.dec0,
                side_deg: 3.0,
                n: 4000,
                min_sep_deg: 0.015,
                mag_lo: 10.0,
                mag_hi: 14.5,
            },
        );
        write_1476_db(decoy.path(), "t50", &other);
        let mut p = params_for(&s, truth.ra0, truth.dec0);
        p.db_path = decoy.path().to_path_buf();
        p.search_radius = deg(0.5);
        match solve_image(&s.img, &p) {
            Err(ArcsecError::InsufficientQuads { found: 0, required }) => {
                assert!(required >= 3);
            }
            other => panic!("expected InsufficientQuads, got {other:?}"),
        }
    }

    #[test]
    fn a_corrupt_catalogue_tile_is_skipped_not_fatal() {
        let truth = TruthWcs::new(deg(120.0), deg(-40.0), 5.0, 0.0, false, 360, 300);
        let s = scene(truth, Db::Areas1476, 120, 9);
        // Declare an unsupported record size in every tile.
        for entry in std::fs::read_dir(s.dir.path()).unwrap() {
            let path = entry.unwrap().path();
            let mut bytes = std::fs::read(&path).unwrap();
            bytes[109] = 7;
            std::fs::write(&path, bytes).unwrap();
        }
        let mut p = params_for(&s, truth.ra0, truth.dec0);
        p.search_radius = 0.0;
        assert!(matches!(
            solve_image(&s.img, &p),
            Err(ArcsecError::InsufficientQuads { .. })
        ));
    }

    #[test]
    fn a_blank_frame_reports_insufficient_stars() {
        let dir = TempDir::new("blank");
        write_1476_db(dir.path(), "t50", &[]);
        let mut rng = Rng::new(3);
        let img = ImageBuffer {
            data: (0..200 * 200)
                .map(|_| (1000.0 + 5.0 * rng.gauss()) as f32)
                .collect(),
            width: 200,
            height: 200,
        };
        let p = SolveParams {
            ra_hint: 0.0,
            dec_hint: 0.0,
            fov: deg(0.3),
            search_radius: deg(1.0),
            quad_tolerance: 0.007,
            hfd_min: 1.5,
            max_stars: 500,
            db_path: dir.path().to_path_buf(),
            db_name: "t50".into(),
            binning: 1,
            method: SolveMethod::Quads,
            threads: 1,
            speed: SearchSpeed::Auto,
        };
        match solve_image(&img, &p) {
            Err(ArcsecError::InsufficientStars { found, required: 5 }) => assert!(found < 5),
            other => panic!("expected InsufficientStars, got {other:?}"),
        }
    }

    #[test]
    fn a_missing_database_is_reported_before_any_detection() {
        let dir = TempDir::new("nodb");
        let p = SolveParams {
            ra_hint: 0.0,
            dec_hint: 0.0,
            fov: deg(1.0),
            search_radius: deg(1.0),
            quad_tolerance: 0.007,
            hfd_min: 1.5,
            max_stars: 500,
            db_path: dir.path().to_path_buf(),
            db_name: "d50".into(),
            binning: 1,
            method: SolveMethod::Quads,
            threads: 1,
            speed: SearchSpeed::Auto,
        };
        match solve_image(&ImageBuffer::new(64, 64), &p) {
            Err(ArcsecError::CatalogNotFound(path)) => assert_eq!(path, dir.path()),
            other => panic!("expected CatalogNotFound, got {other:?}"),
        }
    }

    #[test]
    fn solve_image_rejects_a_bad_search_radius_or_fov() {
        let base = SolveParams {
            ra_hint: 0.0,
            dec_hint: 0.0,
            fov: deg(1.0),
            search_radius: 0.1,
            quad_tolerance: 0.007,
            hfd_min: 1.5,
            max_stars: 500,
            db_path: std::path::PathBuf::from("/nonexistent"),
            db_name: "d50".into(),
            binning: 1,
            method: SolveMethod::Quads,
            threads: 1,
            speed: SearchSpeed::Auto,
        };
        let img = ImageBuffer::new(64, 64);
        for (fov, radius) in [
            (f64::NAN, 0.1),
            (-1.0, 0.1),
            (f64::INFINITY, 0.1),
            (0.01, -0.1),
            (0.01, f64::NAN),
            (0.01, f64::INFINITY),
        ] {
            let p = SolveParams {
                fov,
                search_radius: radius,
                ..base.clone()
            };
            assert!(
                matches!(solve_image(&img, &p), Err(ArcsecError::InvalidParameter(_))),
                "fov {fov}, radius {radius}"
            );
        }
    }

    #[test]
    fn format_radec_roundtrip() {
        let s = format_radec(deg(160.875), deg(-59.524));
        assert!(s.contains("10:"), "RA hours: {s}");
        assert!(s.contains('-'), "dec sign: {s}");
    }
}
