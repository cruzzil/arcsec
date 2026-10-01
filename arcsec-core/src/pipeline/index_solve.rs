//! Blind solving with arcsec's own index ([`crate::index`]).
//!
//! 1. **Detect** stars, brightest first.
//! 2. **Image patterns** — 4-star groups the index is likely to hold: all 4-subsets
//!    of the brightest stars, spread over the frame (at most two per cell of a 6×6
//!    grid), plus, round every bright star, the five brightest stars within a ladder
//!    of window radii — the image's version of the index's disc groups.
//! 3. **Lookup** — each pattern's descriptor keys are looked up in the tiers whose
//!    disc size suits the field. A candidate whose implied pixel scale is outside
//!    the allowed range is dropped; otherwise the canonical vertex order gives the
//!    star correspondence, a 4-point affine fit gives the field (rejected if it is
//!    not close to a similarity), and the field votes in (RA, Dec, ln scale).
//! 4. **Rank** — the vote regions (`sky_votes`) are scored by projecting
//!    the index's stars through each hypothesis and counting those that land on a
//!    detected star, refitting once on the matches.
//! 5. **Accept** — the best-scoring hypotheses, in order, go to the ordinary hinted
//!    solver ([`super::solve_image`]) with the hypothesis as the hint, its scale as
//!    the field size and no search radius. Its star-level verification is the only
//!    acceptance test, so a blind solve is held to exactly the hinted solver's
//!    standard. The index never decides on its own that a field is found.

use core::f64::consts::PI;
use std::collections::{HashMap, HashSet};

use super::sky_votes::{SkyVotes, medoid};
use super::solver::{SolveParams, solve_image};
use crate::detection::get_background;
use crate::detection::stars::find_stars_with_background;
use crate::error::{ArcsecError, Result};
use crate::index::format::{BlindIndex, TierInfo};
use crate::index::pattern::{
    Tangent, canonical, descriptor, fit_affine, invert_affine, orderings, probe_keys, shape_error,
    unit,
};
use crate::types::{ImageBuffer, WcsSolution};

/// Brightest stars (spread over the frame) whose every 4-subset is tried.
const N_GLOBAL: usize = 20;
/// Bright stars around which local windows are built.
const N_WINDOW_ANCHORS: usize = 150;
/// Window radii, as divisors of the image's short side.
const WINDOW_DIVISORS: [f64; 7] = [16.0, 11.0, 8.0, 5.6, 4.0, 2.8, 2.0];
/// Patterns smaller than this (longest edge, pixels) are too noisy to look up.
const MIN_PATTERN_PX: f64 = 12.0;
/// Largest departure from a similarity transform a quad match may show.
const MAX_SHAPE_ERROR: f64 = 0.06;
/// Vote regions scored by projection.
const MAX_REGIONS: usize = 3000;
/// Hypotheses handed to the hinted solver, at most.
const MAX_VERIFY: usize = 6;
/// Ranking score (projected stars landing on detections, plus [`VOTE_WEIGHT`] per
/// extra agreeing quad) a hypothesis needs before it is worth a hinted solve. The
/// hinted solver's own threshold (30 verified stars) is what accepts.
const MIN_COARSE: usize = 8;
/// Ranking weight of each vote beyond the first.
const VOTE_WEIGHT: usize = 2;
/// Leading hypotheses checked against the star database before any hinted solve.
const CHECK_TOP: usize = 24;
/// Significance against the star database a hypothesis needs to be worth a hinted
/// solve. A true field scores tens; chance, about zero.
const MIN_CHECK: usize = 5;
/// Catalogue stars read for that check.
const CHECK_STARS: usize = 300;

/// Parameters for [`index_solve`].
#[derive(Debug, Clone)]
pub struct IndexSolveParams {
    /// Smallest plausible pixel scale of the image as solved (after binning),
    /// arcseconds per pixel.
    pub scale_lo: f64,
    /// Largest plausible pixel scale, arcseconds per pixel.
    pub scale_hi: f64,
    /// Only consider fields centred within this distance (radians) of
    /// (`ra`, `dec`): `Some((ra, dec, radius))` for a search limited by `-r`
    /// around a hint, `None` for the whole sky.
    pub within: Option<(f64, f64, f64)>,
}

/// What an index solve did, for logs and benchmarks.
#[derive(Debug, Clone, Default)]
pub struct IndexSolveStats {
    /// Tiers searched.
    pub tiers: usize,
    /// Image patterns looked up.
    pub patterns: usize,
    /// Index entries returned by the lookups.
    pub candidates: usize,
    /// Field hypotheses that passed the scale and shape checks.
    pub hypotheses: usize,
    /// Vote regions scored.
    pub regions: usize,
    /// Hinted solves attempted.
    pub verified: usize,
    /// Rank (0-based) of the accepted hypothesis among the scored regions.
    pub accepted_rank: Option<usize>,
    /// Projection score of the accepted (or else the best) hypothesis.
    pub best_score: usize,
}

/// A field hypothesis: an affine map from image pixels to the standard coordinates
/// (radians) of a tangent plane, and the field centre and scale it implies.
#[derive(Debug, Clone, Copy)]
struct Hyp {
    ra: f64,
    dec: f64,
    /// Radians per pixel.
    scale: f64,
    tp_ra: f64,
    tp_dec: f64,
    aff: [f64; 6],
}

/// Detected stars on a grid, for nearest-neighbour queries.
struct DetGrid {
    cell: f64,
    cells: HashMap<(i32, i32), Vec<u32>>,
    pts: Vec<(f64, f64)>,
}

impl DetGrid {
    fn new(pts: Vec<(f64, f64)>, cell: f64) -> Self {
        let mut cells: HashMap<(i32, i32), Vec<u32>> = HashMap::new();
        for (i, &(x, y)) in pts.iter().enumerate() {
            cells
                .entry(((x / cell).floor() as i32, (y / cell).floor() as i32))
                .or_default()
                .push(i as u32);
        }
        Self { cell, cells, pts }
    }

    /// Nearest detection within `r` (≤ the cell size) of (x, y).
    fn nearest(&self, x: f64, y: f64, r: f64) -> Option<u32> {
        let (cx, cy) = (
            (x / self.cell).floor() as i32,
            (y / self.cell).floor() as i32,
        );
        let mut best: Option<(f64, u32)> = None;
        for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(v) = self.cells.get(&(cx + dx, cy + dy)) {
                    for &i in v {
                        let (px, py) = self.pts[i as usize];
                        let d = (px - x).powi(2) + (py - y).powi(2);
                        if d <= r * r && best.is_none_or(|(bd, _)| d < bd) {
                            best = Some((d, i));
                        }
                    }
                }
            }
        }
        best.map(|b| b.1)
    }
}

/// Detection cap passed to the star finder: high enough never to trim by SNR.
const DETECT_ALL: usize = 20_000;

/// Background-subtracted flux in a circular aperture of 1.5 HFD (2–25 px) round
/// (x, y), the background being the median of an annulus outside it.
fn aperture_flux(img: &ImageBuffer, x: f64, y: f64, hfd: f64) -> f64 {
    let r = (1.5 * hfd).clamp(2.0, 25.0);
    let (r_in, r_out) = (r + 3.0, r + 7.0);
    let (w, h) = (img.width as i64, img.height as i64);
    let (x0, y0) = (x.round() as i64, y.round() as i64);
    let ro = r_out.ceil() as i64;
    let mut ann: Vec<f32> = Vec::new();
    let mut inner: Vec<f32> = Vec::new();
    for py in (y0 - ro).max(0)..=(y0 + ro).min(h - 1) {
        for px in (x0 - ro).max(0)..=(x0 + ro).min(w - 1) {
            let d = ((px as f64 - x).powi(2) + (py as f64 - y).powi(2)).sqrt();
            let v = img.data[(py * w + px) as usize];
            if d <= r {
                inner.push(v);
            } else if d >= r_in && d <= r_out {
                ann.push(v);
            }
        }
    }
    if ann.is_empty() {
        return 0.0;
    }
    let mid = ann.len() / 2;
    let bg = f64::from(*ann.select_nth_unstable_by(mid, f32::total_cmp).1);
    inner.iter().map(|&v| f64::from(v) - bg).sum()
}

/// The 4-star groups of detections to look up, as sorted index sets.
fn image_quads(pts: &[(f64, f64)], w: f64, h: f64) -> Vec<[usize; 4]> {
    let mut seen: HashSet<[usize; 4]> = HashSet::new();
    let mut out = Vec::new();
    let mut push = |mut q: [usize; 4]| {
        q.sort_unstable();
        if seen.insert(q) {
            out.push(q);
        }
    };

    // Global: the brightest stars, at most two per cell of a 6×6 grid so a bright
    // cluster or a nebula's knots cannot take every slot.
    let mut per_cell: HashMap<(i32, i32), u32> = HashMap::new();
    let picks: Vec<usize> = (0..pts.len())
        .filter(|&i| {
            let (x, y) = pts[i];
            let c = (
                ((x / w * 6.0) as i32).clamp(0, 5),
                ((y / h * 6.0) as i32).clamp(0, 5),
            );
            let n = per_cell.entry(c).or_insert(0);
            *n += 1;
            *n <= 2
        })
        .take(N_GLOBAL)
        .collect();
    let n = picks.len();
    for a in 0..n {
        for b in a + 1..n {
            for c in b + 1..n {
                for d in c + 1..n {
                    push([picks[a], picks[b], picks[c], picks[d]]);
                }
            }
        }
    }

    // Local windows: the five brightest detections within a radius of each bright
    // star, over a ladder of radii.
    let short = w.min(h);
    for &(ax, ay) in pts.iter().take(N_WINDOW_ANCHORS) {
        for div in WINDOW_DIVISORS {
            let r2 = (short / div).powi(2);
            let win: Vec<usize> = (0..pts.len())
                .filter(|&j| (pts[j].0 - ax).powi(2) + (pts[j].1 - ay).powi(2) <= r2)
                .take(5)
                .collect();
            let m = win.len();
            if m < 4 {
                continue;
            }
            for a in 0..m {
                for b in a + 1..m {
                    for c in b + 1..m {
                        for d in c + 1..m {
                            push([win[a], win[b], win[c], win[d]]);
                        }
                    }
                }
            }
        }
    }
    out
}

/// Field hypotheses from one image pattern.
#[allow(clippy::too_many_arguments)]
fn hypotheses_for(
    index: &BlindIndex,
    tiers: &[TierInfo],
    img: &[(f64, f64); 4],
    centre: (f64, f64),
    scale_lo: f64,
    scale_hi: f64,
    keys: &mut Vec<u64>,
    out: &mut Vec<Hyp>,
) -> usize {
    let Some((order, totals, edge_px)) = canonical(img) else {
        return 0;
    };
    if edge_px < MIN_PATTERN_PX {
        return 0;
    }
    let ordered = order.map(|k| img[k]);
    probe_keys(&descriptor(&ordered), keys);
    let alts = orderings(order, totals);
    let mut n_cand = 0;
    for tier in tiers {
        // No pattern of the tier is wider than its disc's diameter.
        if edge_px * scale_lo > 2.0 * tier.radius * 1.02 {
            continue;
        }
        for &k in keys.iter() {
            for pi in index.lookup(tier, k) {
                n_cand += 1;
                let Some(stars) = index.quad(pi) else {
                    continue;
                };
                let us = stars.map(|s| unit(f64::from(s.ra), f64::from(s.dec)));
                let csum = [
                    us[0][0] + us[1][0] + us[2][0] + us[3][0],
                    us[0][1] + us[1][1] + us[2][1] + us[3][1],
                    us[0][2] + us[1][2] + us[2][2] + us[3][2],
                ];
                let Some(tp) = Tangent::at(csum) else {
                    continue;
                };
                let mut cat = [(0.0, 0.0); 4];
                let mut ok = true;
                for (c, u) in cat.iter_mut().zip(&us) {
                    match tp.project(u) {
                        Some(xy) => *c = xy,
                        None => ok = false,
                    }
                }
                if !ok {
                    continue;
                }
                let Some((_, _, edge_cat)) = canonical(&cat) else {
                    continue;
                };
                let scale = edge_cat / edge_px;
                if scale < scale_lo || scale > scale_hi {
                    continue;
                }
                // The best-fitting of the plausible correspondences.
                let mut best: Option<(f64, [f64; 6])> = None;
                for o in &alts {
                    let p = o.map(|k| img[k]);
                    let Some(m) = fit_affine(&p, &cat) else {
                        continue;
                    };
                    let err = shape_error(&m);
                    let s = (m[0] * m[4] - m[1] * m[3]).abs().sqrt();
                    if err <= MAX_SHAPE_ERROR
                        && (s / scale - 1.0).abs() <= MAX_SHAPE_ERROR
                        && best.is_none_or(|(e, _)| err < e)
                    {
                        best = Some((err, m));
                    }
                }
                let Some((_, m)) = best else { continue };
                let (xi, eta) = (
                    m[0] * centre.0 + m[1] * centre.1 + m[2],
                    m[3] * centre.0 + m[4] * centre.1 + m[5],
                );
                let (ra, dec) = tp.deproject(xi, eta);
                let (tp_ra, tp_dec) = tp.deproject(0.0, 0.0);
                out.push(Hyp {
                    ra,
                    dec,
                    scale: (m[0] * m[4] - m[1] * m[3]).abs().sqrt(),
                    tp_ra,
                    tp_dec,
                    aff: m,
                });
            }
        }
    }
    n_cand
}

/// A detection paired with a projected index star: pixel position, standard
/// coordinates.
type Pair = ((f64, f64), (f64, f64));

/// Score a hypothesis by projecting index stars onto the detections; refit once on
/// the matches and rescore. Returns the score and the refined hypothesis.
fn score(
    index: &BlindIndex,
    deepest_tier: u8,
    h: &Hyp,
    det: &DetGrid,
    w: f64,
    hgt: f64,
) -> (usize, Hyp) {
    let diag = w.hypot(hgt);
    let radius = 0.5 * diag * h.scale * 1.1;
    let Some(tp) = Tangent::at_radec(h.tp_ra, h.tp_dec) else {
        return (0, *h);
    };
    let mut cat: Vec<(f64, f64)> = Vec::new();
    let cos_r = radius.cos();
    let cu = unit(h.ra, h.dec);
    index.stars_near(h.ra, h.dec, radius, |s| {
        if s.tier > deepest_tier {
            return;
        }
        let u = unit(f64::from(s.ra), f64::from(s.dec));
        if u[0] * cu[0] + u[1] * cu[1] + u[2] * cu[2] < cos_r {
            return;
        }
        if let Some(xy) = tp.project(&u) {
            cat.push(xy);
        }
    });

    // Matches, as a significance over what chance gives: `n` projected stars each
    // land within `r` of one of the detections with probability `det_density·πr²`.
    // Without this a hypothesis at a far too coarse scale, which projects hundreds
    // of stars into the frame, outscores the true one on raw count.
    let density = det.pts.len() as f64 / (w * hgt);
    let significance = |m: usize, n_in: usize, r: f64| -> usize {
        let e = n_in as f64 * (density * PI * r * r).min(1.0);
        ((m as f64 - e) / (e + 1.0).sqrt()).max(0.0).round() as usize
    };
    let count = |aff: &[f64; 6], r: f64, pairs: &mut Vec<Pair>| -> (usize, usize) {
        pairs.clear();
        let Some(inv) = invert_affine(aff) else {
            return (0, 0);
        };
        let mut used: HashSet<u32> = HashSet::new();
        let mut n_in = 0;
        for &(xi, eta) in &cat {
            let x = inv[0] * xi + inv[1] * eta + inv[2];
            let y = inv[3] * xi + inv[4] * eta + inv[5];
            if x < 0.0 || y < 0.0 || x >= w || y >= hgt {
                continue;
            }
            n_in += 1;
            if let Some(i) = det.nearest(x, y, r)
                && used.insert(i)
            {
                pairs.push((det.pts[i as usize], (xi, eta)));
            }
        }
        (pairs.len(), n_in)
    };

    let mut pairs = Vec::new();
    let r1 = (0.012 * diag).clamp(4.0, det.cell);
    let (n1, in1) = count(&h.aff, r1, &mut pairs);
    let s1 = significance(n1, in1, r1);
    if n1 < 6 {
        return (s1, *h);
    }
    let (p, q): (Vec<_>, Vec<_>) = pairs.iter().copied().unzip();
    let Some(m) = fit_affine(&p, &q) else {
        return (s1, *h);
    };
    if shape_error(&m) > MAX_SHAPE_ERROR {
        return (s1, *h);
    }
    let r2 = (0.004 * diag).clamp(3.0, r1);
    let (n2, in2) = count(&m, r2, &mut pairs);
    let c = ((w - 1.0) / 2.0, (hgt - 1.0) / 2.0);
    let (ra, dec) = tp.deproject(
        m[0] * c.0 + m[1] * c.1 + m[2],
        m[3] * c.0 + m[4] * c.1 + m[5],
    );
    let refined = Hyp {
        ra,
        dec,
        scale: (m[0] * m[4] - m[1] * m[3]).abs().sqrt(),
        aff: m,
        ..*h
    };
    (significance(n2, in2, r2), refined)
}

/// Significance of a hypothesis against the star database: its brightest
/// [`CHECK_STARS`] stars in the field, projected through the hypothesis, counted on
/// the detections and compared with chance. `None` if the catalogue cannot be read
/// (the caller then does not filter on it).
fn catalog_check(template: &SolveParams, hy: &Hyp, det: &DetGrid, w: f64, h: f64) -> Option<usize> {
    let fov = hy.scale * w.max(h);
    let stars = crate::catalog::read_catalog_stars(
        &template.db_path,
        &template.db_name,
        hy.ra,
        hy.dec,
        fov,
        CHECK_STARS,
    )
    .ok()?;
    let tp = Tangent::at_radec(hy.tp_ra, hy.tp_dec)?;
    let inv = invert_affine(&hy.aff)?;
    let r = (0.012 * w.hypot(h)).clamp(4.0, det.cell);
    let mut used: HashSet<u32> = HashSet::new();
    let (mut n_in, mut m) = (0usize, 0usize);
    for s in &stars {
        let Some((xi, eta)) = tp.project(&unit(s.ra, s.dec)) else {
            continue;
        };
        let x = inv[0] * xi + inv[1] * eta + inv[2];
        let y = inv[3] * xi + inv[4] * eta + inv[5];
        if x < 0.0 || y < 0.0 || x >= w || y >= h {
            continue;
        }
        n_in += 1;
        if let Some(i) = det.nearest(x, y, r)
            && used.insert(i)
        {
            m += 1;
        }
    }
    let e = n_in as f64 * (det.pts.len() as f64 / (w * h) * PI * r * r).min(1.0);
    Some(((m as f64 - e) / (e + 1.0).sqrt()).max(0.0).round() as usize)
}

/// Development aid: with `ARCSEC_INDEX_DEBUG=ra,dec,crpix1,crpix2,cd11,cd12,cd21,cd22`
/// (degrees, FITS pixels) describing the true WCS, report how many of the in-play
/// tiers' patterns lie in the image, are detected, hash the same, and were
/// generated as image quads.
fn debug_patterns(
    index: &BlindIndex,
    tiers: &[TierInfo],
    pts: &[(f64, f64)],
    quads: &[[usize; 4]],
    w: f64,
    h: f64,
) {
    let Ok(t) = std::env::var("ARCSEC_INDEX_DEBUG") else {
        return;
    };
    let v: Vec<f64> = t.split(',').filter_map(|x| x.trim().parse().ok()).collect();
    if v.len() != 8 {
        return;
    }
    let Some(tp) = Tangent::at_radec(v[0].to_radians(), v[1].to_radians()) else {
        return;
    };
    let det = v[4] * v[7] - v[5] * v[6];
    let to_px = |ra: f64, dec: f64| -> Option<(f64, f64)> {
        let (xi, eta) = tp.project(&unit(ra, dec))?;
        // FITS xi is east-positive too, but CD maps to (RA-ish, Dec) with RA
        // increasing east: xi_fits = xi.
        let (x, y) = (xi.to_degrees(), eta.to_degrees());
        let dx = (v[7] * x - v[5] * y) / det;
        let dy = (-v[6] * x + v[4] * y) / det;
        Some((v[2] - 1.0 + dx, v[3] - 1.0 + dy))
    };
    let qset: HashSet<[usize; 4]> = quads.iter().copied().collect();
    let grid = DetGrid::new(pts.to_vec(), 8.0);
    {
        // Star-level: detected fraction of the index's stars in the frame, by mag.
        let mut bins = [(0usize, 0usize, 0usize); 20];
        let (ra0, dec0) = (v[0].to_radians(), v[1].to_radians());
        index.stars_near(
            ra0,
            dec0,
            (w.hypot(h) * 0.5 * v[7].abs()).to_radians(),
            |s| {
                let Some((x, y)) = to_px(f64::from(s.ra), f64::from(s.dec)) else {
                    return;
                };
                if x < 5.0 || y < 5.0 || x > w - 5.0 || y > h - 5.0 {
                    return;
                }
                let b = ((f64::from(s.mag) / 100.0).max(0.0) as usize).min(19);
                bins[b].0 += 1;
                if let Some(i) = grid.nearest(x, y, 4.0) {
                    bins[b].1 += 1;
                    bins[b].2 += i as usize;
                }
            },
        );
        for (m, b) in bins.iter().enumerate() {
            if b.0 > 0 {
                log::info!(
                    "Index debug: mag {m}: {} in frame, {} detected, mean rank {}",
                    b.0,
                    b.1,
                    b.2 / b.1.max(1)
                );
            }
        }
    }
    let mut keys = Vec::new();
    for t in tiers {
        let (mut inside, mut detected, mut key_ok, mut generated) = (0, 0, 0, 0);
        let mut ranks: Vec<usize> = Vec::new();
        for pi in t.patterns() {
            let Some(stars) = index.quad(pi) else {
                continue;
            };
            let px: Vec<Option<(f64, f64)>> = stars
                .iter()
                .map(|s| to_px(f64::from(s.ra), f64::from(s.dec)))
                .collect();
            if !px
                .iter()
                .all(|p| p.is_some_and(|(x, y)| x > 5.0 && y > 5.0 && x < w - 5.0 && y < h - 5.0))
            {
                continue;
            }
            inside += 1;
            let ids: Vec<Option<u32>> = px
                .iter()
                .map(|p| p.and_then(|(x, y)| grid.nearest(x, y, 4.0)))
                .collect();
            if ids.iter().any(Option::is_none) {
                continue;
            }
            detected += 1;
            let mut q = [0usize; 4];
            for (k, i) in ids.iter().enumerate() {
                q[k] = i.unwrap_or(0) as usize;
            }
            ranks.push(*q.iter().max().unwrap_or(&0));
            let img = q.map(|i| pts[i]);
            let Some((o, _, _)) = canonical(&img) else {
                continue;
            };
            probe_keys(&descriptor(&o.map(|k| img[k])), &mut keys);
            let mut sorted = q;
            sorted.sort_unstable();
            if keys.iter().any(|&k| index.lookup(t, k).contains(&pi)) {
                key_ok += 1;
                if qset.contains(&sorted) {
                    generated += 1;
                }
            }
        }
        ranks.sort_unstable();
        log::info!(
            "Index debug: tier {:.2}°: {inside} patterns in frame, {detected} detected, \
             {key_ok} same key, {generated} generated; brightness rank of faintest member \
             (median) {:?}",
            t.radius.to_degrees(),
            ranks.get(ranks.len() / 2)
        );
    }
}

/// Development aid: with `ARCSEC_INDEX_TRUTH=<ra°>,<dec°>` in the environment, log
/// where the hypothesis nearest the true centre ranked. Costs nothing otherwise.
fn debug_truth(scored: &[(usize, usize, usize, Hyp)], fov: f64) {
    let Ok(t) = std::env::var("ARCSEC_INDEX_TRUTH") else {
        return;
    };
    let v: Vec<f64> = t.split(',').filter_map(|x| x.trim().parse().ok()).collect();
    if v.len() != 2 {
        return;
    }
    let tu = unit(v[0].to_radians(), v[1].to_radians());
    let near = scored.iter().enumerate().find(|(_, s)| {
        let u = unit(s.3.ra, s.3.dec);
        (u[0] * tu[0] + u[1] * tu[1] + u[2] * tu[2])
            .clamp(-1.0, 1.0)
            .acos()
            < 0.25 * fov
    });
    match near {
        Some((rank, s)) => log::info!(
            "Index truth: rank {rank} of {}, score {}, votes {}",
            scored.len(),
            s.0,
            s.1
        ),
        None => log::info!("Index truth: no hypothesis within a quarter field of the truth"),
    }
}

/// Solve `img` with no position hint, using `index` for the position and the
/// hinted solver (configured by `template`, whose hint, field size and search
/// radius are replaced) to accept it.
///
/// # Errors
///
/// [`ArcsecError::InsufficientStars`] with fewer than 6 detected stars;
/// [`ArcsecError::InsufficientQuads`] when no hypothesis verifies (`found` is the
/// best projection score); [`ArcsecError::InvalidParameter`] for a bad scale range;
/// errors from the hinted solver's catalogue access.
pub fn index_solve(
    img: &ImageBuffer,
    index: &BlindIndex,
    template: &SolveParams,
    params: &IndexSolveParams,
) -> Result<(WcsSolution, IndexSolveStats)> {
    let as2rad = PI / 180.0 / 3600.0;
    if !(params.scale_lo > 0.0 && params.scale_hi >= params.scale_lo && params.scale_hi.is_finite())
    {
        return Err(ArcsecError::InvalidParameter(format!(
            "pixel scale range {}..{} arcsec/px",
            params.scale_lo, params.scale_hi
        )));
    }
    let (scale_lo, scale_hi) = (params.scale_lo * as2rad, params.scale_hi * as2rad);
    let mut stats = IndexSolveStats::default();

    // ── Detect ────────────────────────────────────────────────────────────────
    // Every detection, not the `-s` brightest by SNR: the index's patterns are
    // made of each region's brightest stars, and SNR is a poor brightness order at
    // the bright end (a saturated star's flat top and wide aperture give it a lower
    // SNR than a fainter, sharper one), so the list is re-ranked by flux.
    let bg = get_background(img, template.max_stars);
    let (stars, _) = find_stars_with_background(
        img,
        &bg,
        template.hfd_min,
        DETECT_ALL,
        img.width,
        img.height,
    );
    if stars.len() < 6 {
        return Err(ArcsecError::InsufficientStars {
            found: stars.len(),
            required: 6,
        });
    }
    let mut ranked: Vec<(f64, (f64, f64))> = stars
        .0
        .iter()
        .map(|s| (aperture_flux(img, s.x, s.y, s.hfd), (s.x, s.y)))
        .collect();
    ranked.sort_by(|a, b| b.0.total_cmp(&a.0));
    ranked.truncate(template.max_stars.max(N_WINDOW_ANCHORS));
    let pts: Vec<(f64, f64)> = ranked.iter().map(|r| r.1).collect();
    let (w, h) = (img.width as f64, img.height as f64);
    let short = w.min(h);

    // ── Tiers whose discs suit this field ─────────────────────────────────────
    let (f_lo, f_hi) = (short * scale_lo, short * scale_hi);
    let tiers: Vec<(usize, TierInfo)> = index
        .tiers()
        .iter()
        .copied()
        .enumerate()
        .filter(|(_, t)| f_hi >= 2.0 * t.radius && f_lo <= 16.0 * t.radius)
        .collect();
    stats.tiers = tiers.len();
    let Some(&(deepest, _)) = tiers.last() else {
        log::info!(
            "Index: no tier suits a {:.2}°–{:.2}° field.",
            f_lo.to_degrees(),
            f_hi.to_degrees()
        );
        return Err(ArcsecError::InsufficientQuads {
            found: 0,
            required: MIN_COARSE,
        });
    };
    let tier_infos: Vec<TierInfo> = tiers.iter().map(|t| t.1).collect();

    // ── Patterns → hypotheses (parallel over patterns) ────────────────────────
    let quads = image_quads(&pts, w, h);
    stats.patterns = quads.len();
    debug_patterns(index, &tier_infos, &pts, &quads, w, h);
    let centre = ((w - 1.0) / 2.0, (h - 1.0) / 2.0);
    let threads = crate::max_threads().max(1);
    let chunk = quads.len().div_ceil(threads).max(1);
    let parts: Vec<(usize, Vec<Hyp>)> = std::thread::scope(|scope| {
        let handles: Vec<_> = quads
            .chunks(chunk)
            .map(|part| {
                let (pts, tier_infos) = (&pts, &tier_infos);
                scope.spawn(move || {
                    let mut keys = Vec::new();
                    let mut out = Vec::new();
                    let mut n = 0;
                    for q in part {
                        let p = q.map(|i| pts[i]);
                        n += hypotheses_for(
                            index, tier_infos, &p, centre, scale_lo, scale_hi, &mut keys, &mut out,
                        );
                    }
                    (n, out)
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
            .collect()
    });

    let s_mid = (scale_lo * scale_hi).sqrt();
    let step = (0.05 * w.max(h) * s_mid).clamp(0.002f64.to_radians(), 1f64.to_radians());
    let mut votes: SkyVotes<Hyp> = SkyVotes::new(step, 0.05);
    for (n, hyps) in parts {
        stats.candidates += n;
        stats.hypotheses += hyps.len();
        for hy in hyps {
            if let Some((ra, dec, r)) = params.within {
                let (u, v) = (unit(ra, dec), unit(hy.ra, hy.dec));
                if u[0] * v[0] + u[1] * v[1] + u[2] * v[2] < r.min(PI).cos() {
                    continue;
                }
            }
            votes.add(hy.ra, hy.dec, hy.scale, hy);
        }
    }
    let regions = votes.regions(MAX_REGIONS);
    stats.regions = regions.len();
    log::info!(
        "Index: {} stars, {} tiers, {} patterns, {} candidates, {} hypotheses, {} regions.",
        pts.len(),
        stats.tiers,
        stats.patterns,
        stats.candidates,
        stats.hypotheses,
        stats.regions
    );

    // ── Rank by projection (parallel) ─────────────────────────────────────────
    let det = DetGrid::new(
        pts.iter().copied().take(300).collect(),
        (0.012 * w.hypot(h)).max(8.0),
    );
    let reps: Vec<(usize, Hyp)> = regions
        .iter()
        .map(|r| (r.votes, r.members[medoid(r.members, |h| (h.ra, h.dec))]))
        .collect();
    let chunk = reps.len().div_ceil(threads).max(1);
    let mut scored: Vec<(usize, usize, usize, Hyp)> = std::thread::scope(|scope| {
        let handles: Vec<_> = reps
            .chunks(chunk)
            .enumerate()
            .map(|(ci, part)| {
                let det = &det;
                scope.spawn(move || {
                    part.iter()
                        .enumerate()
                        .map(|(k, (votes, hy))| {
                            let (s, refined) = score(index, deepest as u8, hy, det, w, h);
                            // Independent quads agreeing on a field are strong
                            // evidence on their own; the projection score counts
                            // only the index's few bright stars in the frame.
                            (
                                s + VOTE_WEIGHT * (votes - 1),
                                *votes,
                                ci * chunk + k,
                                refined,
                            )
                        })
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
            .collect()
    });
    scored.sort_by(|a, b| b.0.cmp(&a.0).then(b.1.cmp(&a.1)).then(a.2.cmp(&b.2)));
    stats.best_score = scored.first().map_or(0, |s| s.0);
    debug_truth(&scored, w.max(h) * s_mid);

    // ── Check against the star database, then accept through the hinted solver ─
    // A hinted solve that fails costs a detection pass and a few spiral positions,
    // so the leading hypotheses are first checked against the full catalogue (the
    // index holds only a few bright stars per field) and tried in that order.
    let mut checked: Vec<(usize, usize, usize, usize, Hyp)> = scored
        .iter()
        .enumerate()
        .take(CHECK_TOP)
        .filter(|(_, s)| s.0 >= MIN_COARSE)
        .map(|(rank, &(sc, votes, _, hy))| {
            let c = catalog_check(template, &hy, &det, w, h).unwrap_or(MIN_CHECK);
            (c, rank, sc, votes, hy)
        })
        .collect();
    checked.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
    let mut last_err = None;
    for &(check, rank, sc, votes, hy) in checked.iter().take(MAX_VERIFY) {
        if check < MIN_CHECK {
            break;
        }
        stats.verified += 1;
        let fov = hy.scale * w.max(h);
        log::info!(
            "Index: verifying #{rank}: RA={:.4}° Dec={:.4}° scale={:.3}\"/px score={sc} votes={votes} check={check}",
            hy.ra.to_degrees(),
            hy.dec.to_degrees(),
            hy.scale / as2rad
        );
        let p = SolveParams {
            ra_hint: hy.ra,
            dec_hint: hy.dec,
            fov,
            search_radius: 0.0,
            ..template.clone()
        };
        match solve_image(img, &p) {
            Ok(wcs) => {
                stats.accepted_rank = Some(rank);
                stats.best_score = sc;
                return Ok((wcs, stats));
            }
            Err(e @ (ArcsecError::CatalogNotFound(_) | ArcsecError::CatalogIo(_))) => {
                return Err(e);
            }
            Err(e) => last_err = Some(e),
        }
    }
    log::info!(
        "Index: no hypothesis verified (best score {}, {} tried).",
        stats.best_score,
        stats.verified
    );
    let _ = last_err;
    Err(ArcsecError::InsufficientQuads {
        found: stats.best_score,
        required: MIN_COARSE,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::{BuildParams, TierSpec, build_index};
    use crate::pipeline::solver::{SearchSpeed, SolveMethod};
    use crate::test_support::{Rng, SkySpec, TempDir, TruthWcs, random_sky, render, write_1476_db};

    fn deg(d: f64) -> f64 {
        d.to_radians()
    }

    const W: usize = 600;
    const H: usize = 480;
    const SCALE: f64 = 6.0;

    /// A field and a decoy region in one database, its index, and an image of the
    /// field.
    fn scene(t: &TruthWcs, seed: u64, dir: &std::path::Path) -> (ImageBuffer, BlindIndex) {
        let mut rng = Rng::new(seed);
        let spec = SkySpec {
            ra0: t.ra0,
            dec0: t.dec0,
            side_deg: 2.5,
            n: 2500,
            min_sep_deg: 12.0 * SCALE / 3600.0,
            mag_lo: 8.0,
            mag_hi: 14.5,
        };
        let mut sky = random_sky(&mut rng, &spec);
        let img = render(t, &sky, 1.3, 1000.0, 8.0, 30_000.0, &mut rng);
        sky.extend(random_sky(
            &mut rng,
            &SkySpec {
                ra0: deg(80.0),
                dec0: deg(-20.0),
                ..spec
            },
        ));
        write_1476_db(dir, "t80", &sky);
        let built = build_index(
            &BuildParams {
                db_path: dir.to_path_buf(),
                db_name: "t80".into(),
                tiers: vec![TierSpec {
                    radius_deg: 0.15,
                    mag_cap: 14.5,
                    members: 6,
                }],
                threads: 2,
            },
            |_| {},
        )
        .unwrap();
        let path = dir.join("t80.arcsecix");
        built.write(&path).unwrap();
        (img, BlindIndex::open(&path).unwrap())
    }

    fn template(dir: &std::path::Path) -> SolveParams {
        SolveParams {
            ra_hint: 0.0,
            dec_hint: 0.0,
            fov: deg(W as f64 * SCALE / 3600.0),
            search_radius: 0.0,
            quad_tolerance: 0.007,
            hfd_min: 1.5,
            max_stars: 500,
            db_path: dir.to_path_buf(),
            db_name: "t80".into(),
            binning: 1,
            method: SolveMethod::Quads,
            speed: SearchSpeed::Auto,
            threads: 2,
        }
    }

    #[test]
    fn finds_a_field_blind_and_verifies_it() {
        for (rot, mirrored, seed) in [(30.0, false, 5), (-100.0, true, 6)] {
            let dir = TempDir::new("ixsolve");
            let t = TruthWcs::new(deg(310.0), deg(44.0), SCALE, rot, mirrored, W, H);
            let (img, ix) = scene(&t, seed, dir.path());
            let (wcs, stats) = index_solve(
                &img,
                &ix,
                &template(dir.path()),
                &IndexSolveParams {
                    scale_lo: SCALE / 1.2,
                    scale_hi: SCALE * 1.2,
                    within: None,
                },
            )
            .unwrap_or_else(|e| panic!("rot {rot} mirrored {mirrored}: {e}"));
            let err = t.max_error_arcsec(&wcs);
            assert!(err < 2.0, "worst error {err:.2}\" ({stats:?})");
        }
    }
}
