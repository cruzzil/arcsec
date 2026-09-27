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
    let mut vote_map: std::collections::HashMap<(i32, i32), Vec<HypEntry>> =
        std::collections::HashMap::new();
    let mut n_matches = 0usize;
    let mut n_hyp = 0usize;

    let mut hits_scratch: Vec<usize> = Vec::with_capacity(32);
    for ie in &img_entries {
        index.find_code_matches_into(&ie.code, tol, &mut hits_scratch);
        n_matches += hits_scratch.len();

        for &idx in &hits_scratch {
            let idx_entry = &index.entries[idx];
            if let Some(h) = hyp_from_entry(ie, idx_entry, img.width, img.height) {
                let ra_bin = (h.est_ra / VOTE_STEP) as i32;
                let dec_bin = ((h.est_dec + PI / 2.0) / VOTE_STEP) as i32;
                vote_map.entry((ra_bin, dec_bin)).or_default().push(h);
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
        vote_map.len(),
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

    let mut vote_cells: Vec<_> = vote_map.iter().collect();
    // Ties break on the (ra_bin, dec_bin) key. sort_by_key is stable, so without
    // this the order among equally-voted cells came from HashMap iteration, which
    // is randomly seeded per process; the verification loop below takes strictly
    // greater scores and stops early, so two runs on one image could disagree.
    vote_cells.sort_by(|a, b| b.1.len().cmp(&a.1.len()).then_with(|| a.0.cmp(b.0)));

    let mut best_score = 0usize;
    let mut best_ra = 0.0f64;
    let mut best_dec = 0.0f64;

    'outer: for (_, hyps) in &vote_cells {
        let h = &hyps[0];
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
            hyps.len(),
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
        "Blind ({}): best verified score = {} (threshold {}, cells={}).",
        parity_label,
        best_score,
        MIN_VERIFY_SCORE,
        vote_cells.len(),
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
