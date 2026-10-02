//! Blind astrometric solver using an Astrometry.net index file.
//!
//! Algorithm:
//! 1. Detect stars in the image; build all image quads (triangles or 4-star quads,
//!    depending on index DIMQUADS).
//! 2. For each image quad, compute the code using the astrometry.net formula.
//!    Both image parities are tried (CDELT1<0 normal, CDELT1>0 flipped).
//!    Canonical form: `CX ≤ 0.5` for triangles; `CX + DX ≤ 1` and `CX ≤ DX` for
//!    quads (see `make_quad4`).
//! 3. Scale-filter: only try image quads whose A-B axis pixel length
//!    corresponds to the index's angular scale range given the image FOV.
//!    (Skipped when `fov_deg == 0` — truly blind solve.)
//! 4. For each (image quad, index quad) code match, compute a precise
//!    WCS hypothesis from the star affine transform and vote in 0.1° bins.
//! 5. Verify vote cells, most-voted first, by projecting index stars; return the
//!    cell with the highest verified star-match count (stopping early once a cell
//!    reaches `EARLY_STOP_SCORE`).

use core::f64::consts::PI;

use super::sky_votes::{SkyVotes, medoid};
use crate::catalog::anet::AnetIndex;
use crate::detection::get_background;
use crate::detection::stars::find_stars_with_background;
use crate::error::{ArcsecError, Result};
use crate::math::coords::equatorial_standard;
use crate::math::lsq::solve_plate_constants;
use crate::types::StarList;
use crate::wcs::output::derive_wcs;

/// Parameters for [`blind_solve`].
#[derive(Debug, Clone)]
pub struct BlindSolveParams {
    /// Code-space matching tolerance (Euclidean distance between codes).
    pub quad_tolerance: f64,
    /// Minimum HFD for valid stars (pixels).
    pub hfd_min: f64,
    /// Maximum number of image stars to detect.
    pub max_stars: usize,
    /// Binning already applied to the image. Informational only: the blind solver
    /// returns a sky position, which binning does not affect.
    pub binning: usize,
    /// Image height in degrees (used for scale filtering). 0 = auto (no filter).
    pub fov_deg: f64,
}

// ── Image quad types ────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct ImageEntry {
    /// [CX, CY] for triangles, [CX, CY, DX, DY] for quads; unused slots are 0.
    code: [f64; 4],
    /// Pixel positions of stars [A, B, C[, D]]; only first `n_stars` are valid.
    stars: [(f64, f64); 4],
    n_stars: usize,
    /// Pixel length of the A-B axis, for scale filtering.
    d_ab_px: f64,
}

/// Compute the (CX, CY) code of point P relative to the axis A→B, where
/// `(adx, ady)` is A→P, `(abx, aby)` is A→B and `scale` is `|AB|²`.
///
/// The frame is astrometry.net's 45°-rotated one, so that A maps to (0, 0) and B to
/// (1, 1): `cos_t = (aby+abx)/|AB|²`, `sin_t = (aby-abx)/|AB|²`.
///
/// `parity_flip = false` (normal):  east = −x (CDELT1 < 0, standard FITS).
///   `code_x = −adx·sin + ady·cos`   (east-dominant)
///   `code_y =  adx·cos + ady·sin`
///
/// `parity_flip = true`  (flipped): east = +x (CDELT1 > 0).
///   `code_x =  adx·cos + ady·sin`
///   `code_y = −adx·sin + ady·cos`
fn code_for_point(
    adx: f64,
    ady: f64,
    abx: f64,
    aby: f64,
    scale: f64,
    parity_flip: bool,
) -> (f64, f64) {
    let cos_theta = (aby + abx) / scale;
    let sin_theta = (aby - abx) / scale;
    if parity_flip {
        (
            adx * cos_theta + ady * sin_theta,
            -adx * sin_theta + ady * cos_theta,
        )
    } else {
        (
            -adx * sin_theta + ady * cos_theta,
            adx * cos_theta + ady * sin_theta,
        )
    }
}

/// Build a DIMQUADS=3 image triangle.
fn make_triangle(
    a: (f64, f64),
    b: (f64, f64),
    c: (f64, f64),
    parity_flip: bool,
) -> Option<ImageEntry> {
    let (xa, ya) = a;
    let (xb, yb) = b;
    let abx = xb - xa;
    let aby = yb - ya;
    let scale = abx * abx + aby * aby;
    if scale < 1.0 {
        return None;
    }
    let d_ab_px = scale.sqrt();

    let (mut cx, mut cy) = code_for_point(c.0 - xa, c.1 - ya, abx, aby, scale, parity_flip);
    let (mut pa, mut pb) = (a, b);

    // Canonical form: CX ≤ 0.5
    if cx > 0.5 {
        core::mem::swap(&mut pa, &mut pb);
        cx = 1.0 - cx;
        cy = 1.0 - cy;
    }

    Some(ImageEntry {
        code: [cx, cy, 0.0, 0.0],
        stars: [pa, pb, c, (0.0, 0.0)],
        n_stars: 3,
        d_ab_px,
    })
}

/// Build a DIMQUADS=4 image quad.
///
/// Stars are provided as [A, B, C, D] where A-B is the longest pair in the group.
fn make_quad4(
    a: (f64, f64),
    b: (f64, f64),
    c: (f64, f64),
    d: (f64, f64),
    parity_flip: bool,
) -> Option<ImageEntry> {
    let (xa, ya) = a;
    let (xb, yb) = b;
    let abx = xb - xa;
    let aby = yb - ya;
    let scale = abx * abx + aby * aby;
    if scale < 1.0 {
        return None;
    }
    let d_ab_px = scale.sqrt();

    let (mut cx, mut cy) = code_for_point(c.0 - xa, c.1 - ya, abx, aby, scale, parity_flip);
    let (mut dx, mut dy) = code_for_point(d.0 - xa, d.1 - ya, abx, aby, scale, parity_flip);
    let (mut pa, mut pb, mut pc, mut pd) = (a, b, c, d);

    // Canonical step 1: CX + DX ≤ 1 (swap A↔B, invert all codes).
    //
    // Not `cx > 0.5`, which is the *triangle* rule and is what this used to test.
    // For a 4-star quad astrometry.net canonicalises on the sum, and the two are
    // not the same set: cx ≤ 0.5 ∧ cx ≤ dx admits (0.4, 0.7), whose sum is 1.1,
    // while every quad in the index satisfies cx + dx ≤ 1. Measured against the
    // 4100 series, 100% of 7.9M index quads are inside cx + dx ≤ 1, and the old
    // rule put 21% of image quads outside it — permanently unmatchable.
    //
    // The swap sends (cx, dx) to (1-cx, 1-dx), so a sum above 1 becomes a sum
    // below it; step 2 then swaps C↔D, which leaves the sum alone. Both
    // invariants therefore hold on exit.
    if cx + dx > 1.0 {
        core::mem::swap(&mut pa, &mut pb);
        cx = 1.0 - cx;
        cy = 1.0 - cy;
        dx = 1.0 - dx;
        dy = 1.0 - dy;
    }

    // Canonical step 2: CX ≤ DX (if equal: CY ≤ DY) — swap C↔D
    #[allow(clippy::float_cmp)] // exact tie-break, as astrometry.net does it
    let swap_cd = dx < cx || (dx == cx && dy < cy);
    if swap_cd {
        core::mem::swap(&mut pc, &mut pd);
        core::mem::swap(&mut cx, &mut dx);
        core::mem::swap(&mut cy, &mut dy);
    }

    Some(ImageEntry {
        code: [cx, cy, dx, dy],
        stars: [pa, pb, pc, pd],
        n_stars: 4,
        d_ab_px,
    })
}

#[inline]
fn dist_sq(p: (f64, f64), q: (f64, f64)) -> f64 {
    let dx = p.0 - q.0;
    let dy = p.1 - q.1;
    dx * dx + dy * dy
}

/// Build all image triangles (DIMQUADS=3) from the N brightest detected stars.
fn build_triangles(stars: &StarList, n: usize, parity_flip: bool) -> Vec<ImageEntry> {
    let n = stars.len().min(n);
    let mut out = Vec::new();
    if n < 3 {
        return out;
    }
    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                let a = (stars.0[i].x, stars.0[i].y);
                let b = (stars.0[j].x, stars.0[j].y);
                let c = (stars.0[k].x, stars.0[k].y);
                // Each triple gives 3 axis choices (A-B, A-C, B-C)
                if let Some(e) = make_triangle(a, b, c, parity_flip) {
                    out.push(e);
                }
                if let Some(e) = make_triangle(a, c, b, parity_flip) {
                    out.push(e);
                }
                if let Some(e) = make_triangle(b, c, a, parity_flip) {
                    out.push(e);
                }
            }
        }
    }
    out
}

/// Build all image quads (DIMQUADS=4) from the N brightest detected stars.
///
/// For each combination of 4 stars, the A-B pair is the one with maximum
/// pixel separation (matching astrometry.net's convention).
fn build_quads4(stars: &StarList, n: usize, parity_flip: bool) -> Vec<ImageEntry> {
    const PAIR_INDICES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

    let n = stars.len().min(n);
    let mut out = Vec::new();
    if n < 4 {
        return out;
    }

    let pts: Vec<(f64, f64)> = stars.0[..n].iter().map(|s| (s.x, s.y)).collect();

    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                for l in (k + 1)..n {
                    let group = [pts[i], pts[j], pts[k], pts[l]];

                    // Find pair with maximum distance → becomes the A-B axis.
                    let (ai, bi) = PAIR_INDICES
                        .iter()
                        .copied()
                        .max_by(|&(a1, b1), &(a2, b2)| {
                            dist_sq(group[a1], group[b1]).total_cmp(&dist_sq(group[a2], group[b2]))
                        })
                        .unwrap();

                    // The other two, in ascending order.
                    let mut rest = (0..4).filter(|&x| x != ai && x != bi);
                    let (Some(ci), Some(di)) = (rest.next(), rest.next()) else {
                        unreachable!("four stars minus a pair leaves two");
                    };

                    if let Some(e) =
                        make_quad4(group[ai], group[bi], group[ci], group[di], parity_flip)
                    {
                        out.push(e);
                    }
                }
            }
        }
    }
    out
}

/// Build all image entries for a given parity and DIMQUADS value.
fn build_image_entries(
    stars: &StarList,
    n: usize,
    parity_flip: bool,
    dim_quads: usize,
) -> Vec<ImageEntry> {
    match dim_quads {
        3 => build_triangles(stars, n, parity_flip),
        4 => build_quads4(stars, n, parity_flip),
        _ => vec![],
    }
}

// ── WCS hypothesis ─────────────────────────────────────────────────────────────

struct HypEntry {
    est_ra: f64,
    est_dec: f64,
    ref_ra: f64,
    ref_dec: f64,
    plate: crate::types::PlateConstants,
    /// RA/Dec of the catalog stars used to derive this WCS (quad stars).
    /// These are excluded from `verify_score` so they don't inflate the count.
    quad_cat_ra: [f64; 4],
    quad_cat_dec: [f64; 4],
    n_quad: usize,
}

fn hyp_from_entry(
    img_entry: &ImageEntry,
    idx_entry: &crate::catalog::anet::AnetIndexEntry,
    img_w: usize,
    img_h: usize,
) -> Option<HypEntry> {
    let ref_ra = idx_entry.center_ra;
    let ref_dec = idx_entry.center_dec;
    let n = img_entry.n_stars.min(idx_entry.n_stars);

    let mut img_px = [(0.0f64, 0.0f64); 4];
    let mut cat_xy = [(0.0f64, 0.0f64); 4];
    for i in 0..n {
        img_px[i] = img_entry.stars[i];
        cat_xy[i] = equatorial_standard(
            ref_ra,
            ref_dec,
            idx_entry.star_ra[i],
            idx_entry.star_dec[i],
            1.0,
        );
    }

    let plate = solve_plate_constants(&img_px[..n], &cat_xy[..n]).ok()?;
    let wcs = derive_wcs(ref_ra, ref_dec, &plate, img_w, img_h);

    let mut quad_cat_ra = [0.0f64; 4];
    let mut quad_cat_dec = [0.0f64; 4];
    quad_cat_ra[..n].copy_from_slice(&idx_entry.star_ra[..n]);
    quad_cat_dec[..n].copy_from_slice(&idx_entry.star_dec[..n]);

    Some(HypEntry {
        est_ra: wcs.ra0,
        est_dec: wcs.dec0,
        ref_ra,
        ref_dec,
        plate,
        quad_cat_ra,
        quad_cat_dec,
        n_quad: n,
    })
}

fn sky_to_px(
    h: &HypEntry,
    star_ra: f64,
    star_dec: f64,
    det: f64,
    img_w: usize,
    img_h: usize,
) -> Option<(f64, f64)> {
    let (xs, ys) = equatorial_standard(h.ref_ra, h.ref_dec, star_ra, star_dec, 1.0);
    let xc = xs - h.plate.c;
    let yc = ys - h.plate.f;
    let px = (h.plate.e * xc - h.plate.b * yc) / det;
    let py = (h.plate.a * yc - h.plate.d * xc) / det;
    if px >= 0.0 && px < img_w as f64 && py >= 0.0 && py < img_h as f64 {
        Some((px, py))
    } else {
        None
    }
}

/// Count how many index stars project within `sqrt(match_px_sq)` of any detected star,
/// excluding the catalog stars that were used to derive this WCS hypothesis.
/// `stars_by_dec` must be sorted ascending by DEC for binary-search pre-filtering.
fn verify_score(
    h: &HypEntry,
    stars_by_dec: &[crate::catalog::anet::AnetStar],
    det_stars: &[(f64, f64)],
    fov_rad: f64,
    match_px_sq: f64,
    img_w: usize,
    img_h: usize,
) -> usize {
    let det = h.plate.a * h.plate.e - h.plate.b * h.plate.d;
    if det.abs() < 1e-15 {
        return 0;
    }

    // Pre-filter to DEC band — reduces work from O(n_total) to O(n_in_strip).
    let half_fov = fov_rad * 1.1;
    let dec_lo = h.est_dec - half_fov;
    let dec_hi = h.est_dec + half_fov;
    let lo = stars_by_dec.partition_point(|s| s.dec < dec_lo);
    let hi = stars_by_dec.partition_point(|s| s.dec <= dec_hi);
    let stars_in_band = &stars_by_dec[lo..hi];

    let min_cos_sep = half_fov.cos();
    let (sin_e, cos_e) = h.est_dec.sin_cos();

    let mut score = 0usize;
    'star: for star in stars_in_band {
        let cos_sep = sin_e * star.dec.sin() + cos_e * star.dec.cos() * (star.ra - h.est_ra).cos();
        if cos_sep < min_cos_sep {
            continue;
        }

        // Skip the catalog stars that were used to derive this WCS hypothesis —
        // they always project back to their matched pixel positions and would inflate
        // the score for false positives equally as much as for true positives.
        for i in 0..h.n_quad {
            let dra = (star.ra - h.quad_cat_ra[i]).abs();
            let ddec = (star.dec - h.quad_cat_dec[i]).abs();
            if dra < 1e-9 && ddec < 1e-9 {
                continue 'star;
            }
        }

        let Some((px, py)) = sky_to_px(h, star.ra, star.dec, det, img_w, img_h) else {
            continue;
        };

        let nearest_sq = det_stars.iter().fold(f64::INFINITY, |acc, &(dx, dy)| {
            let d = (dx - px) * (dx - px) + (dy - py) * (dy - py);
            acc.min(d)
        });
        if nearest_sq <= match_px_sq {
            score += 1;
        }
    }
    score
}

// ── Core single-parity pass ────────────────────────────────────────────────────

const N_ENTRY_STARS: usize = 30;
const VOTE_STEP: f64 = 0.1 * PI / 180.0; // 0.1° bins
/// Vote bin width in `ln(pixel scale)`: hypotheses whose scales differ by more than
/// about 5% are not agreeing on a field.
const VOTE_LOG_SCALE_STEP: f64 = 0.05;
const MATCH_PX: f64 = 5.0;
// Quad stars are excluded from verify_score; empirical false-positive max ≈ 17
// across 20 test fields.  MIN_VERIFY_SCORE=18 eliminates all observed false
// positives while keeping the weakest true positive (Piscis-Aus, score=18).
// EARLY_STOP: once any cell scores this high we are very confident — stop early.
const MIN_VERIFY_SCORE: usize = 18;
const EARLY_STOP_SCORE: usize = 20;

/// Run one full blind-solve pass for a single parity assumption.
///
/// `d_px` is the permitted A-B axis length in pixels: the index's angular scale band
/// mapped through the image scale.
///
/// Returns the best (ra, dec) estimate and its verification score.
fn run_blind_pass(
    img: &crate::types::ImageBuffer,
    index: &AnetIndex,
    stars: &StarList,
    parity_flip: bool,
    fov_rad: f64,
    d_px: core::ops::RangeInclusive<f64>,
    tol: f64,
) -> (f64, f64, usize) {
    let n_entry_stars = stars.len().min(N_ENTRY_STARS);
    let all_entries = build_image_entries(stars, n_entry_stars, parity_flip, index.dim_quads);

    let img_entries: Vec<ImageEntry> = if *d_px.start() <= 0.0 && *d_px.end() == f64::INFINITY {
        all_entries
    } else {
        all_entries
            .into_iter()
            .filter(|e| d_px.contains(&e.d_ab_px))
            .collect()
    };

    let parity_label = if parity_flip { "flipped" } else { "normal" };
    let n_total = if index.dim_quads == 4 {
        let n = n_entry_stars;
        n * (n - 1) * (n - 2) * (n - 3) / 24 // C(n,4)
    } else {
        n_entry_stars * (n_entry_stars - 1) * (n_entry_stars - 2) / 6 * 3 // C(n,3) * 3
    };
    log::info!(
        "Blind ({}): {}/{} entry-axis pairs in scale range [{:.0},{:.0}] px.",
        parity_label,
        img_entries.len(),
        n_total,
        d_px.start(),
        d_px.end(),
    );

    if img_entries.len() < 3 {
        return (0.0, 0.0, 0);
    }

    // ── Phase C: match and accumulate WCS hypotheses ──────────────────────────
    let mut votes: SkyVotes<HypEntry> = SkyVotes::new(VOTE_STEP, VOTE_LOG_SCALE_STEP);
    let mut n_matches = 0usize;
    let mut n_hyp = 0usize;

    let mut hits_scratch: Vec<usize> = Vec::with_capacity(32);
    for (k, ie) in img_entries.iter().enumerate() {
        // The caller reports the cancellation; an empty pass is all it needs.
        if k % 64 == 0 && crate::cancel::is_cancelled() {
            return (0.0, 0.0, 0);
        }
        index.find_code_matches_into(&ie.code, tol, &mut hits_scratch);
        n_matches += hits_scratch.len();

        for &idx in &hits_scratch {
            let idx_entry = &index.entries[idx];
            if let Some(h) = hyp_from_entry(ie, idx_entry, img.width, img.height) {
                let scale = (h.plate.a * h.plate.e - h.plate.b * h.plate.d).abs().sqrt();
                votes.add(h.est_ra, h.est_dec, scale, h);
                n_hyp += 1;
            }
        }
    }

    log::info!(
        "Blind ({}): {} entries, {} code matches → {} WCS hypotheses, {} vote cells.",
        parity_label,
        img_entries.len(),
        n_matches,
        n_hyp,
        votes.len(),
    );

    if n_hyp == 0 {
        return (0.0, 0.0, 0);
    }

    // Sort a copy of index stars by DEC once; verify_score uses binary search
    // to restrict to the DEC band around each hypothesis (~55× speedup).
    let mut stars_by_dec = index.stars.clone();
    stars_by_dec.sort_unstable_by(|a, b| a.dec.total_cmp(&b.dec));

    // ── Phase D: verify top-K cells ───────────────────────────────────────────
    let det_stars: Vec<(f64, f64)> = stars.0.iter().map(|s| (s.x, s.y)).collect();

    // Regions of agreement, strongest first (see `sky_votes`): each is represented
    // by the medoid of its strongest bucket, so a region is not scored on a stray
    // member. The order is deterministic, which the early stop below relies on.
    let regions = votes.regions(usize::MAX);

    let mut best_score = 0usize;
    let mut best_ra = 0.0f64;
    let mut best_dec = 0.0f64;

    'outer: for region in &regions {
        if crate::cancel::is_cancelled() {
            return (0.0, 0.0, 0);
        }
        let h = &region.members[medoid(region.members, |h| (h.est_ra, h.est_dec))];
        let sc = verify_score(
            h,
            &stars_by_dec,
            &det_stars,
            fov_rad,
            MATCH_PX * MATCH_PX,
            img.width,
            img.height,
        );
        log::debug!(
            "Blind verify ({}): RA={:.2}° Dec={:.2}° votes={} score={}",
            parity_label,
            h.est_ra.to_degrees(),
            h.est_dec.to_degrees(),
            region.votes,
            sc,
        );
        if sc > best_score {
            best_score = sc;
            best_ra = h.est_ra;
            best_dec = h.est_dec;
            if best_score >= EARLY_STOP_SCORE {
                break 'outer;
            }
        }
    }

    log::info!(
        "Blind ({}): best verified score = {} (threshold {}, regions={}).",
        parity_label,
        best_score,
        MIN_VERIFY_SCORE,
        regions.len(),
    );

    (best_ra, best_dec, best_score)
}

// ── Public API ─────────────────────────────────────────────────────────────────

/// Blind-solve one index file, returning `(ra_rad, dec_rad, verify_score)`.
///
/// Tries both parities; fails if best score < `MIN_VERIFY_SCORE`.
/// Call for each candidate index file and keep the result with the highest score.
///
/// # Errors
///
/// - [`ArcsecError::InsufficientStars`] if fewer than 5 stars are detected.
/// - [`ArcsecError::InsufficientQuads`] if no hypothesis verifies; `found` is the
///   best verification score reached.
pub fn blind_solve(
    img: &crate::types::ImageBuffer,
    index: &AnetIndex,
    params: &BlindSolveParams,
) -> Result<(f64, f64, usize)> {
    // ── Phase A: detect stars ─────────────────────────────────────────────────
    crate::cancel::progress(crate::cancel::stage::BLIND_INDEX, -1.0);
    let bg = get_background(img, params.max_stars);
    let (stars, stars_raw) = find_stars_with_background(
        img,
        &bg,
        params.hfd_min,
        params.max_stars,
        img.width,
        img.height,
    );
    let mut stars = stars;
    // Patterns are built from the first stars in the list, which must therefore be
    // the brightest. Detection only sorts by SNR when it has more than `max_stars`
    // to trim; below that the list is in scan order (cascade level, then rows), and
    // a sparse field would build every pattern from the top strip of the frame.
    stars.0.sort_by(|a, b| b.snr.total_cmp(&a.snr));
    if stars_raw > params.max_stars {
        stars.0.truncate((params.max_stars / 2).max(50));
    }
    let n = stars.len();
    if n < 5 {
        return Err(ArcsecError::InsufficientStars {
            found: n,
            required: 5,
        });
    }

    log::info!("Blind: {n} stars detected.");

    let (min_d_px, max_d_px) = if params.fov_deg > 0.0 {
        let fov = params.fov_deg.to_radians();
        let lo = index.scale_lo * img.height as f64 / fov * 0.8;
        let hi = index.scale_hi * img.height as f64 / fov * 1.2;
        (lo, hi)
    } else {
        (0.0, f64::INFINITY)
    };

    let verify_fov = if params.fov_deg > 0.0 {
        params.fov_deg.to_radians()
    } else {
        index.scale_hi * 2.0
    };

    let tol = params.quad_tolerance;

    let mut best_ra = 0.0f64;
    let mut best_dec = 0.0f64;
    let mut best_score = 0usize;

    for parity_flip in [false, true] {
        if crate::cancel::is_cancelled() {
            return Err(ArcsecError::Cancelled);
        }
        let (ra, dec, score) = run_blind_pass(
            img,
            index,
            &stars,
            parity_flip,
            verify_fov,
            min_d_px..=max_d_px,
            tol,
        );
        if score > best_score {
            best_score = score;
            best_ra = ra;
            best_dec = dec;
        }
        if best_score >= EARLY_STOP_SCORE {
            break;
        }
    }
    if crate::cancel::is_cancelled() {
        return Err(ArcsecError::Cancelled);
    }

    if best_score >= MIN_VERIFY_SCORE {
        log::info!(
            "Blind: estimated centre RA={:.3}°, Dec={:.3}° (score={best_score})",
            best_ra.to_degrees(),
            best_dec.to_degrees(),
        );
        Ok((best_ra, best_dec, best_score))
    } else {
        Err(ArcsecError::InsufficientQuads {
            found: best_score,
            required: MIN_VERIFY_SCORE,
        })
    }
}

#[cfg(test)]
mod canonical_form {
    use super::*;

    /// Deterministic xorshift, so this test cannot flake.
    fn rng() -> impl FnMut() -> f64 {
        let mut seed = 0x2545F4914F6CDD1Du64;
        move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    /// Image quads must land in the same code region the index occupies.
    ///
    /// Every quad in the astrometry.net 4100 series satisfies `cx + dx <= 1`
    /// (checked against all 7.9M of them), so an image quad outside that region
    /// cannot match anything at any tolerance. Canonicalising on `cx <= 0.5` —
    /// which is the rule for 3-star triangles — put 21% of image quads outside it.
    #[test]
    fn image_quads_lie_in_the_index_code_region() {
        let mut next = rng();
        let (mut total, mut outside) = (0usize, 0usize);
        for _ in 0..4000 {
            let p: Vec<(f64, f64)> = (0..4).map(|_| (next() * 1000.0, next() * 1000.0)).collect();
            for flip in [false, true] {
                if let Some(e) = make_quad4(p[0], p[1], p[2], p[3], flip) {
                    total += 1;
                    let (cx, dx) = (e.code[0], e.code[2]);
                    if cx + dx > 1.0 + 1e-9 {
                        outside += 1;
                    }
                    assert!(cx <= dx + 1e-9, "second invariant broken: cx={cx} dx={dx}");
                }
            }
        }
        assert!(
            total > 1000,
            "probe built too few quads to be meaningful: {total}"
        );
        assert_eq!(
            outside, 0,
            "{outside} of {total} image quads outside cx + dx <= 1"
        );
    }

    /// Triangles keep the `cx <= 0.5` rule: it is the correct one for 3-star codes.
    #[test]
    fn triangles_keep_the_half_plane_rule() {
        let mut next = rng();
        let mut total = 0usize;
        for _ in 0..4000 {
            let p: Vec<(f64, f64)> = (0..3).map(|_| (next() * 1000.0, next() * 1000.0)).collect();
            if let Some(e) = make_triangle(p[0], p[1], p[2], false) {
                total += 1;
                assert!(e.code[0] <= 0.5 + 1e-9, "cx={} above 0.5", e.code[0]);
            }
        }
        assert!(total > 1000, "probe built too few triangles: {total}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::catalog::anet::AnetIndexEntry;
    use crate::test_support::{
        RawIndex, Rng, SkySpec, SkyStar, TruthWcs, anet_code, random_sky, render, separation,
    };
    use crate::types::{ImageBuffer, Star};

    fn deg(d: f64) -> f64 {
        d.to_radians()
    }

    fn close(a: &[f64], b: &[f64], tol: f64) -> bool {
        a.iter().zip(b).all(|(x, y)| (x - y).abs() <= tol)
    }

    /// `n` random detector positions and where a WCS puts them on the sky.
    fn group_on_detector(wcs: &TruthWcs, rng: &mut Rng, n: usize) -> crate::types::PairedPositions {
        let px: Vec<(f64, f64)> = (0..n)
            .map(|_| (rng.range(20.0, 380.0), rng.range(20.0, 280.0)))
            .collect();
        let sky = px.iter().map(|&(x, y)| wcs.pixel_to_sky(x, y)).collect();
        (sky, px)
    }

    /// The image-side code must equal the code astrometry.net computes on the sky,
    /// for either parity, whatever the rotation. This pins the sign conventions in
    /// `code_for_point` against an independent implementation (`anet_code`).
    #[test]
    fn image_codes_agree_with_sky_codes_for_both_parities() {
        let mut rng = Rng::new(21);
        let (mut built, mut differs) = (0, 0);
        for trial in 0..300 {
            let mirrored = trial % 2 == 1;
            let wcs = TruthWcs::new(
                rng.range(0.0, 2.0 * PI),
                rng.range(-1.4, 1.4),
                rng.range(1.0, 30.0),
                rng.range(-180.0, 180.0),
                mirrored,
                400,
                300,
            );
            let (sky, px) = group_on_detector(&wcs, &mut rng, 4);
            let Some((order, sky_code)) = anet_code(&sky, &[0, 1, 2, 3]) else {
                continue;
            };
            built += 1;
            // Build the image quad the way build_quads4 does: A-B = longest pair.
            let pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
            let &(ai, bi) = pairs
                .iter()
                .max_by(|p, q| dist_sq(px[p.0], px[p.1]).total_cmp(&dist_sq(px[q.0], px[q.1])))
                .unwrap();
            let rest: Vec<usize> = (0..4).filter(|&k| k != ai && k != bi).collect();
            let e = make_quad4(px[ai], px[bi], px[rest[0]], px[rest[1]], mirrored).unwrap();
            assert!(
                close(&e.code, &sky_code, 2e-3),
                "trial {trial}: image {:?} vs sky {:?}",
                e.code,
                sky_code
            );
            // The same stars, in the same canonical order.
            for (k, &id) in order.iter().enumerate() {
                assert!(dist_sq(e.stars[k], px[id as usize]) < 1e-18, "star {k}");
            }
            // The wrong parity gives a different code (it is the mirror image).
            let wrong = make_quad4(px[ai], px[bi], px[rest[0]], px[rest[1]], !mirrored).unwrap();
            if !close(&wrong.code, &sky_code, 1e-2) {
                differs += 1;
            }
        }
        assert!(built > 100, "too few buildable quads: {built}");
        assert!(
            differs * 10 > built * 9,
            "parity barely matters? {differs}/{built}"
        );
    }

    #[test]
    fn triangle_codes_agree_with_sky_codes() {
        let mut rng = Rng::new(22);
        let mut built = 0;
        for trial in 0..200 {
            let mirrored = trial % 3 == 0;
            let wcs = TruthWcs::new(1.0, 0.4, 10.0, rng.range(-180.0, 180.0), mirrored, 400, 300);
            let (sky, px) = group_on_detector(&wcs, &mut rng, 3);
            let Some((order, sky_code)) = anet_code(&sky, &[0, 1, 2]) else {
                continue; // an acute triangle: the index would not hold it
            };
            built += 1;
            // Image side: the axis is the longest side, as the index has it.
            let (a, b, c) = (order[0] as usize, order[1] as usize, order[2] as usize);
            let e = make_triangle(px[a], px[b], px[c], mirrored).unwrap();
            assert!(
                close(&e.code[..2], &sky_code[..2], 2e-3),
                "trial {trial}: image {:?} vs sky {:?}",
                e.code,
                sky_code
            );
        }
        assert!(built > 50, "too few buildable triangles: {built}");
    }

    #[test]
    fn codes_are_invariant_under_similarity_transforms() {
        let mut rng = Rng::new(23);
        for _ in 0..500 {
            let p: Vec<(f64, f64)> = (0..4)
                .map(|_| (rng.range(0.0, 500.0), rng.range(0.0, 500.0)))
                .collect();
            let (s, r) = (rng.range(0.2, 5.0), rng.range(-PI, PI));
            let (tx, ty) = (rng.range(-1e3, 1e3), rng.range(-1e3, 1e3));
            let q: Vec<(f64, f64)> = p
                .iter()
                .map(|&(x, y)| {
                    (
                        s * (x * r.cos() - y * r.sin()) + tx,
                        s * (x * r.sin() + y * r.cos()) + ty,
                    )
                })
                .collect();
            let (Some(e1), Some(e2)) = (
                make_quad4(p[0], p[1], p[2], p[3], false),
                make_quad4(q[0], q[1], q[2], q[3], false),
            ) else {
                continue;
            };
            assert!(
                close(&e1.code, &e2.code, 1e-9),
                "{:?} vs {:?}",
                e1.code,
                e2.code
            );
            assert!((e2.d_ab_px - s * e1.d_ab_px).abs() < 1e-6);
            // Mirroring x is exactly a parity flip.
            let m: Vec<(f64, f64)> = p.iter().map(|&(x, y)| (-x, y)).collect();
            let e3 = make_quad4(m[0], m[1], m[2], m[3], true).unwrap();
            assert!(
                close(&e1.code, &e3.code, 1e-9),
                "{:?} vs {:?}",
                e1.code,
                e3.code
            );
        }
    }

    #[test]
    fn degenerate_axes_are_rejected() {
        assert!(make_triangle((5.0, 5.0), (5.5, 5.0), (9.0, 9.0), false).is_none());
        assert!(make_quad4((5.0, 5.0), (5.0, 5.5), (1.0, 1.0), (2.0, 2.0), false).is_none());
        assert!(build_image_entries(&StarList::default(), 30, false, 4).is_empty());
        let one = StarList(vec![Star {
            x: 1.0,
            y: 1.0,
            snr: 1.0,
            hfd: 1.0,
        }]);
        assert!(build_image_entries(&one, 30, false, 5).is_empty());
    }

    #[test]
    fn image_entry_counts() {
        let mut rng = Rng::new(24);
        let stars = StarList(
            (0..9)
                .map(|_| Star {
                    x: rng.range(0.0, 400.0),
                    y: rng.range(0.0, 400.0),
                    snr: 10.0,
                    hfd: 2.0,
                })
                .collect(),
        );
        // C(9,3) triangles, each on all three axes; C(9,4) quads. Only the first
        // `n` stars take part.
        assert_eq!(build_triangles(&stars, 9, false).len(), 84 * 3);
        assert_eq!(build_quads4(&stars, 9, false).len(), 126);
        assert_eq!(build_quads4(&stars, 6, true).len(), 15);
        assert!(build_quads4(&stars, 3, false).is_empty());
        assert!(build_triangles(&stars, 2, false).is_empty());
    }

    // ── Hypotheses and verification ─────────────────────────────────────────────

    /// The index entry for a group of sky stars, as the index builder writes it,
    /// and the image entry that detecting those stars at `px` would produce.
    fn matched_pair(sky: &[(f64, f64)], px: &[(f64, f64)]) -> (ImageEntry, AnetIndexEntry) {
        let raw = RawIndex::build(4, sky, &[vec![0, 1, 2, 3]]);
        let idx = raw.to_index().entries.remove(0);
        let at = |k: usize| px[raw.quads[0][k] as usize];
        let img = ImageEntry {
            code: idx.code,
            stars: [at(0), at(1), at(2), at(3)],
            n_stars: 4,
            d_ab_px: 0.0,
        };
        (img, idx)
    }

    #[test]
    fn a_matched_quad_gives_the_true_field_centre() {
        let wcs = TruthWcs::new(deg(250.0), deg(-60.0), 4.0, 33.0, false, 400, 300);
        let mut rng = Rng::new(25);
        let (sky, px) = group_on_detector(&wcs, &mut rng, 4);
        let (img, idx) = matched_pair(&sky, &px);
        let h = hyp_from_entry(&img, &idx, 400, 300).unwrap();
        let (ra_c, dec_c) = wcs.pixel_to_sky(199.5, 149.5);
        let err = separation(h.est_ra, h.est_dec, ra_c, dec_c).to_degrees() * 3600.0;
        assert!(err < 0.5, "centre error {err}\"");
    }

    #[test]
    fn verify_score_counts_projected_stars_but_not_the_quad_itself() {
        let wcs = TruthWcs::new(deg(10.0), deg(10.0), 5.0, 0.0, false, 400, 300);
        let mut rng = Rng::new(26);
        // 40 stars on a jittered 8 × 5 grid, so no two are within the match radius.
        let px: Vec<(f64, f64)> = (0..40)
            .map(|k| {
                (
                    30.0 + 45.0 * (k % 8) as f64 + rng.range(-5.0, 5.0),
                    30.0 + 55.0 * (k / 8) as f64 + rng.range(-5.0, 5.0),
                )
            })
            .collect();
        let sky: Vec<(f64, f64)> = px.iter().map(|&(x, y)| wcs.pixel_to_sky(x, y)).collect();
        // The quad: two opposite corners and two central stars, inside their circle.
        let (img, idx) = matched_pair(
            &[sky[0], sky[18], sky[21], sky[39]],
            &[px[0], px[18], px[21], px[39]],
        );
        let h = hyp_from_entry(&img, &idx, 400, 300).unwrap();
        let mut by_dec: Vec<crate::catalog::AnetStar> = sky
            .iter()
            .map(|&(ra, dec)| crate::catalog::AnetStar { ra, dec })
            .collect();
        // A star far outside the field and one just off the frame edge.
        by_dec.push(crate::catalog::AnetStar {
            ra: deg(100.0),
            dec: deg(10.0),
        });
        let (ra, dec) = wcs.pixel_to_sky(-30.0, 150.0);
        by_dec.push(crate::catalog::AnetStar { ra, dec });
        by_dec.sort_by(|a, b| a.dec.total_cmp(&b.dec));

        let fov = deg(300.0 * 5.0 / 3600.0);
        assert_eq!(verify_score(&h, &by_dec, &px, fov, 25.0, 400, 300), 36);
        // Only half the stars were detected: half the score.
        let detected: Vec<(f64, f64)> = px.iter().copied().step_by(2).collect();
        let s = verify_score(&h, &by_dec, &detected, fov, 25.0, 400, 300);
        assert_eq!(s, 18, "score {s}");
        // A singular plate scores nothing.
        let mut flat = h;
        flat.plate.a = 0.0;
        flat.plate.b = 0.0;
        assert_eq!(verify_score(&flat, &by_dec, &px, fov, 25.0, 400, 300), 0);
    }

    // ── End to end ───────────────────────────────────────────────────────────────

    const FIELD_RA: f64 = 310.0;
    const FIELD_DEC: f64 = 44.0;

    fn field_spec(seedless_n: usize) -> SkySpec {
        SkySpec {
            ra0: deg(FIELD_RA),
            dec0: deg(FIELD_DEC),
            side_deg: 2.0,
            n: seedless_n,
            min_sep_deg: 12.0 * 6.0 / 3600.0,
            mag_lo: 10.0,
            mag_hi: 14.5,
        }
    }

    /// A rendered field and an index holding quads of its brightest stars, plus a
    /// decoy region elsewhere on the sky with its own quads.
    fn blind_scene(truth: &TruthWcs, dim_quads: usize, seed: u64) -> (ImageBuffer, RawIndex) {
        let mut rng = Rng::new(seed);
        let field = random_sky(&mut rng, &field_spec(1200));
        let img = render(truth, &field, 1.3, 1000.0, 8.0, 30_000.0, &mut rng);
        let decoy = random_sky(
            &mut rng,
            &SkySpec {
                ra0: deg(80.0),
                dec0: deg(-20.0),
                ..field_spec(1200)
            },
        );
        (img, index_over(truth, &field, &decoy, dim_quads))
    }

    /// Quads from the brightest ten in-frame stars of `field` and of the decoy's
    /// central region; every star of both goes in the star table.
    fn index_over(
        truth: &TruthWcs,
        field: &[SkyStar],
        decoy: &[SkyStar],
        dim_quads: usize,
    ) -> RawIndex {
        let mut all: Vec<SkyStar> = field.to_vec();
        all.extend_from_slice(decoy);
        let sky: Vec<(f64, f64)> = all.iter().map(|s| (s.ra, s.dec)).collect();
        let brightest = |range: core::ops::Range<usize>, inside: &dyn Fn(&SkyStar) -> bool| {
            let mut ids: Vec<u32> = range
                .filter(|&i| inside(&all[i]))
                .map(|i| i as u32)
                .collect();
            ids.sort_by(|&a, &b| all[a as usize].mag.total_cmp(&all[b as usize].mag));
            ids.truncate(10);
            ids
        };
        let in_frame = |s: &SkyStar| {
            truth
                .sky_to_pixel(s.ra, s.dec)
                .is_some_and(|(x, y)| x > 20.0 && y > 20.0 && x < 380.0 && y < 300.0)
        };
        let decoy_centre = |s: &SkyStar| separation(s.ra, s.dec, deg(80.0), deg(-20.0)) < deg(0.25);
        let mut groups = Vec::new();
        for ids in [
            brightest(0..field.len(), &in_frame),
            brightest(field.len()..all.len(), &decoy_centre),
        ] {
            combos(&ids, dim_quads, &mut Vec::new(), 0, &mut groups);
        }
        RawIndex::build(dim_quads, &sky, &groups)
    }

    fn combos(ids: &[u32], k: usize, cur: &mut Vec<u32>, from: usize, out: &mut Vec<Vec<u32>>) {
        if cur.len() == k {
            out.push(cur.clone());
            return;
        }
        for i in from..ids.len() {
            cur.push(ids[i]);
            combos(ids, k, cur, i + 1, out);
            cur.pop();
        }
    }

    /// `max_stars` is kept below the ~95 stars each scene holds, so detection sorts
    /// its list by SNR; see `blind_solve_uses_the_brightest_stars_when_few_are_found`
    /// for what happens otherwise.
    fn params(fov_deg: f64) -> BlindSolveParams {
        BlindSolveParams {
            quad_tolerance: 0.007,
            hfd_min: 1.5,
            max_stars: 60,
            binning: 1,
            fov_deg,
        }
    }

    fn truth(rot: f64, mirrored: bool) -> TruthWcs {
        TruthWcs::new(deg(FIELD_RA), deg(FIELD_DEC), 6.0, rot, mirrored, 400, 320)
    }

    fn centre_error_arcsec(t: &TruthWcs, ra: f64, dec: f64) -> f64 {
        let (ra_c, dec_c) = t.pixel_to_sky(199.5, 159.5);
        separation(ra, dec, ra_c, dec_c).to_degrees() * 3600.0
    }

    #[test]
    fn blind_solve_finds_the_field_with_quads() {
        let t = truth(57.0, false);
        let (img, raw) = blind_scene(&t, 4, 31);
        let (ra, dec, score) =
            blind_solve(&img, &raw.to_index(), &params(320.0 * 6.0 / 3600.0)).expect("blind solve");
        let err = centre_error_arcsec(&t, ra, dec);
        assert!(err < 5.0, "estimate {err:.2}\" from the centre");
        assert!(score >= MIN_VERIFY_SCORE, "score {score}");
    }

    /// The CLI's `--index` flow, end to end: an index written as a FITS file and
    /// loaded back, a blind estimate from it, then the catalogue solver seeded with
    /// that estimate and a narrowed radius.
    #[test]
    fn index_file_to_blind_estimate_to_catalogue_solve() {
        use crate::pipeline::solver::{SolveMethod, SolveParams, solve_image};
        use crate::test_support::{TempDir, write_1476_db};

        let t = truth(-35.0, false);
        let mut rng = Rng::new(41);
        let field = random_sky(&mut rng, &field_spec(1200));
        let img = render(&t, &field, 1.3, 1000.0, 8.0, 30_000.0, &mut rng);
        let decoy = random_sky(
            &mut rng,
            &SkySpec {
                ra0: deg(80.0),
                dec0: deg(-20.0),
                ..field_spec(1200)
            },
        );
        let raw = index_over(&t, &field, &decoy, 4);

        let dir = TempDir::new("index-flow");
        let index_path = dir.path().join("index-9999.fits");
        std::fs::write(&index_path, raw.fits_bytes()).unwrap();
        let index = crate::catalog::load_anet_index(&index_path).expect("load index");
        let fov_deg = 320.0 * 6.0 / 3600.0;
        let (ra, dec, _) = blind_solve(&img, &index, &params(fov_deg)).expect("blind solve");

        write_1476_db(dir.path(), "t50", &field);
        let wcs = solve_image(
            &img,
            &SolveParams {
                ra_hint: ra,
                dec_hint: dec,
                fov: deg(fov_deg),
                search_radius: deg(fov_deg * 2.0),
                quad_tolerance: 0.007,
                hfd_min: 1.5,
                max_stars: 500,
                db_path: dir.path().to_path_buf(),
                db_name: "t50".into(),
                binning: 1,
                method: SolveMethod::Quads,
                threads: 1,
                speed: crate::pipeline::SearchSpeed::Auto,
            },
        )
        .expect("catalogue solve");
        let err = t.max_error_arcsec(&wcs);
        assert!(err < 1.0, "worst centre/corner error {err:.3}\"");
    }

    #[test]
    fn blind_solve_finds_the_field_with_triangles_and_no_scale_hint() {
        let t = truth(-15.0, false);
        let (img, raw) = blind_scene(&t, 3, 32);
        let (ra, dec, _) = blind_solve(&img, &raw.to_index(), &params(0.0)).expect("blind solve");
        // A three-star fit extrapolated to the centre is looser than a quad's; the
        // estimate only has to land the catalogue solver within a field.
        let err = centre_error_arcsec(&t, ra, dec);
        assert!(err < 20.0, "estimate {err:.2}\" from the centre");
    }

    /// The image entries are built from the first `N_ENTRY_STARS` (30) detected
    /// stars, documented as "the N brightest". But `find_stars_with_background` only
    /// sorts by SNR when it has more than `max_stars` to trim; otherwise the list is
    /// in detection order — cascade level, then raster order down the frame. With
    /// the default `-s 500` and a field holding fewer stars than that, the blind
    /// solver therefore builds its patterns from the stars nearest the top of the
    /// frame, not the brightest, and the index (built from bright stars) rarely
    /// shares a pattern with them.
    ///
    /// Here ~95 stars are detected. With `max_stars = 60` (sorted) this field solves
    /// with a score in the 90s; with `max_stars = 300` (unsorted) not one of the
    /// index's 210 in-field quads has all four stars among the 30 used, and the
    /// best score is 5.
    #[test]
    fn blind_solve_uses_the_brightest_stars_when_few_are_found() {
        let t = truth(0.0, false);
        let (img, raw) = blind_scene(&t, 4, 33);
        let index = raw.to_index();
        let sorted = blind_solve(&img, &index, &params(0.533));
        assert!(sorted.is_ok(), "control: solves when detection sorts");
        let unsorted = blind_solve(
            &img,
            &index,
            &BlindSolveParams {
                max_stars: 300,
                ..params(0.533)
            },
        );
        assert!(unsorted.is_ok(), "{unsorted:?}");
    }

    #[test]
    fn blind_solve_finds_a_mirrored_field() {
        let t = truth(-120.0, true);
        let (img, raw) = blind_scene(&t, 4, 33);
        let (ra, dec, _) = blind_solve(&img, &raw.to_index(), &params(0.533)).expect("blind solve");
        let err = centre_error_arcsec(&t, ra, dec);
        assert!(err < 5.0, "estimate {err:.2}\" from the centre");
    }

    #[test]
    fn blind_solve_rejects_a_field_the_index_does_not_hold() {
        // The index holds a different sky at the same place, plus the decoy.
        let t = truth(0.0, false);
        let (img, _) = blind_scene(&t, 4, 35);
        let mut rng = Rng::new(36);
        let elsewhere = random_sky(&mut rng, &field_spec(1200));
        let raw = index_over(&t, &elsewhere, &[], 4);
        match blind_solve(&img, &raw.to_index(), &params(0.533)) {
            Err(ArcsecError::InsufficientQuads { found, required }) => {
                assert_eq!(required, MIN_VERIFY_SCORE);
                assert!(found < MIN_VERIFY_SCORE);
            }
            other => panic!("expected InsufficientQuads, got {other:?}"),
        }
    }

    #[test]
    fn blind_solve_needs_stars() {
        let t = truth(0.0, false);
        let (_, raw) = blind_scene(&t, 4, 37);
        let img = ImageBuffer {
            data: vec![1000.0; 200 * 200],
            width: 200,
            height: 200,
        };
        assert!(matches!(
            blind_solve(&img, &raw.to_index(), &params(0.5)),
            Err(ArcsecError::InsufficientStars { required: 5, .. })
        ));
    }
}
