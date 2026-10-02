//! Catalogue-seeded quad search: a fallback matcher that does not rank image stars.
//!
//! The quad matcher pairs image and catalogue quads built from each side's own
//! brightest stars, so it needs the two brightness rankings to agree: when the
//! image's brightest stars are not the catalogue's (saturated discs, nebulosity, a
//! cluster core the catalogue resolves and the image does not, a passband far from
//! Gaia's), the two quad sets share too few neighbourhoods to match. This search
//! uses brightness on the catalogue side only, where it can be trusted.
//!
//! Quads are built from the catalogue's brightest stars, each with three of its four
//! nearest bright neighbours. On the image side every detected star takes part,
//! whatever its brightness: a table of all image star pairs sorted by length, and a
//! position hash. For each catalogue quad, the image pairs as long as the quad's
//! widest pair (at the expected pixel scale) are looked up by binary search; each
//! pair, in both orders and both parities, fixes a similarity transform, under
//! which the quad's other two stars must land on detected stars. A transform that
//! passes is scored by how many of the catalogue's brightest stars it puts on a
//! detection, and one that scores enough is handed to the caller's verification.
//!
//! The idea and its constants follow seiza's rank-robust fallback
//! (<https://github.com/theatrus/seiza>, `docs/design/rank-robust-matching.md`,
//! Apache-2.0); the code is arcsec's own.

use crate::types::PlateConstants;

/// Image star pairs shorter than this (pixels) are not indexed: too short to fix a
/// rotation.
const MIN_PAIR_PX: f64 = 8.0;

/// A transform's census must find at least this many times the hits expected by
/// chance (catalogue stars in the frame times the chance of a detection within the
/// census radius of a random point). On a 2.2° TESS crop, 470 detections on
/// 384 × 384 pixels, a random point has a detection within 3.75 px 14% of the time,
/// and wrong transforms passed a fixed floor of 10 by the thousand.
pub const CENSUS_SIGNIFICANCE: f64 = 2.0;

/// A uniform grid over image positions, for "is there a star within `tol` of this
/// point" lookups. Compressed-row layout: one allocation for the whole grid.
struct PosGrid {
    min_x: f64,
    min_y: f64,
    inv_cell: f64,
    nx: usize,
    ny: usize,
    /// `start[c]..start[c + 1]` indexes `items` for cell `c`.
    start: Vec<u32>,
    items: Vec<u32>,
}

impl PosGrid {
    fn new(pos: &[(f64, f64)], cell: f64) -> Self {
        let (mut min_x, mut min_y) = (f64::INFINITY, f64::INFINITY);
        let (mut max_x, mut max_y) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
        for &(x, y) in pos {
            min_x = min_x.min(x);
            max_x = max_x.max(x);
            min_y = min_y.min(y);
            max_y = max_y.max(y);
        }
        if pos.is_empty() {
            (min_x, min_y, max_x, max_y) = (0.0, 0.0, 0.0, 0.0);
        }
        let cell = cell.max(0.5);
        let inv_cell = 1.0 / cell;
        let nx = ((max_x - min_x) * inv_cell) as usize + 1;
        let ny = ((max_y - min_y) * inv_cell) as usize + 1;
        let cell_of = |x: f64, y: f64| {
            let gx = (((x - min_x) * inv_cell) as usize).min(nx - 1);
            let gy = (((y - min_y) * inv_cell) as usize).min(ny - 1);
            gy * nx + gx
        };
        let mut start = vec![0u32; nx * ny + 1];
        for &(x, y) in pos {
            start[cell_of(x, y) + 1] += 1;
        }
        for c in 0..nx * ny {
            start[c + 1] += start[c];
        }
        let mut fill = start.clone();
        let mut items = vec![0u32; pos.len()];
        for (i, &(x, y)) in pos.iter().enumerate() {
            let c = cell_of(x, y);
            items[fill[c] as usize] = i as u32;
            fill[c] += 1;
        }
        Self {
            min_x,
            min_y,
            inv_cell,
            nx,
            ny,
            start,
            items,
        }
    }

    /// The nearest star within `tol` of `(x, y)`.
    fn nearest(&self, pos: &[(f64, f64)], x: f64, y: f64, tol: f64) -> Option<usize> {
        let gx0 = ((x - tol - self.min_x) * self.inv_cell).floor();
        let gy0 = ((y - tol - self.min_y) * self.inv_cell).floor();
        let gx1 = ((x + tol - self.min_x) * self.inv_cell).floor();
        let gy1 = ((y + tol - self.min_y) * self.inv_cell).floor();
        if gx1 < 0.0 || gy1 < 0.0 || gx0 >= self.nx as f64 || gy0 >= self.ny as f64 {
            return None;
        }
        let (gx0, gy0) = (gx0.max(0.0) as usize, gy0.max(0.0) as usize);
        let gx1 = (gx1 as usize).min(self.nx - 1);
        let gy1 = (gy1 as usize).min(self.ny - 1);
        let mut best = None;
        let mut best_d2 = tol * tol;
        for gy in gy0..=gy1 {
            let row = gy * self.nx;
            for c in row + gx0..=row + gx1 {
                for &i in &self.items[self.start[c] as usize..self.start[c + 1] as usize] {
                    let (px, py) = pos[i as usize];
                    let d2 = (px - x) * (px - x) + (py - y) * (py - y);
                    if d2 <= best_d2 {
                        best_d2 = d2;
                        best = Some(i as usize);
                    }
                }
            }
        }
        best
    }
}

/// Cells of at least `probe_tol`, few enough (about a million) that the bit map
/// stays in cache: probes land all over the frame, and on a 4300-pixel frame a
/// map of 2.5-pixel cells (3 million bits) made every probe a cache miss.
const NEAR_MAX_CELLS_PER_SIDE: f64 = 1024.0;

/// One bit per cell, set for every cell within one cell of a star; the cell is at
/// least `probe_tol` wide, so a point whose cell is clear has no star within
/// `probe_tol`.
struct NearMap {
    min_x: f32,
    min_y: f32,
    inv_cell: f32,
    nx: usize,
    ny: usize,
    bits: Vec<u64>,
}

impl NearMap {
    fn new(pos: &[(f64, f64)], grid: &PosGrid, probe_tol: f64) -> Self {
        let (w, h) = (
            grid.nx as f64 / grid.inv_cell,
            grid.ny as f64 / grid.inv_cell,
        );
        // A little over `probe_tol`, so that rounding the f32 probes cannot carry
        // a point within `probe_tol` of a star two cells from it.
        let cell = (1.01 * probe_tol)
            .max(w.max(h) / NEAR_MAX_CELLS_PER_SIDE)
            .max(0.5);
        let inv_cell = 1.0 / cell;
        let nx = (w * inv_cell) as usize + 1;
        let ny = (h * inv_cell) as usize + 1;
        let mut bits = vec![0u64; (nx * ny).div_ceil(64)];
        for &(x, y) in pos {
            let gx = ((x - grid.min_x) * inv_cell) as usize;
            let gy = ((y - grid.min_y) * inv_cell) as usize;
            for cy in gy.saturating_sub(1)..=(gy + 1).min(ny - 1) {
                for cx in gx.saturating_sub(1)..=(gx + 1).min(nx - 1) {
                    let c = cy * nx + cx;
                    bits[c / 64] |= 1 << (c % 64);
                }
            }
        }
        Self {
            min_x: grid.min_x as f32,
            min_y: grid.min_y as f32,
            inv_cell: inv_cell as f32,
            nx,
            ny,
            bits,
        }
    }

    /// False if no star can be within `probe_tol` of `(x, y)`.
    #[inline]
    fn maybe(&self, x: f32, y: f32) -> bool {
        let fx = (x - self.min_x) * self.inv_cell;
        let fy = (y - self.min_y) * self.inv_cell;
        if !(fx >= 0.0 && fy >= 0.0) {
            // Within one cell outside the map can still be within tol of a star.
            return fx > -1.0 && fy > -1.0;
        }
        let (cx, cy) = (fx as usize, fy as usize);
        if cx >= self.nx || cy >= self.ny {
            return cx <= self.nx && cy <= self.ny;
        }
        let c = cy * self.nx + cx;
        self.bits[c / 64] & (1 << (c % 64)) != 0
    }
}

/// The image side of the search: every star's position, a position hash, and the
/// star pairs sorted by length. Built once per image.
pub struct ImageIndex {
    pos: Vec<(f64, f64)>,
    grid: PosGrid,
    /// A point whose cell here is clear has no star within `probe_tol`.
    near: NearMap,
    /// The pairs' lengths, sorted, for the binary search ...
    pair_len: Vec<f32>,
    /// ... each pair's first star and the vector to its second, for the probes ...
    pair_vec: Vec<[f32; 4]>,
    /// ... and the two stars.
    pair_idx: Vec<(u32, u32)>,
    probe_tol: f64,
}

impl ImageIndex {
    /// Index `pos` (pixels) for probes within `probe_tol` pixels, and the pairs
    /// among them no longer than `max_pair_px`.
    #[must_use]
    pub fn new(pos: Vec<(f64, f64)>, max_pair_px: f64, probe_tol: f64) -> Self {
        let grid = PosGrid::new(&pos, probe_tol);
        let near = NearMap::new(&pos, &grid, probe_tol);
        // Sorted by x so only pairs within max_pair_px in x are examined.
        let mut order: Vec<u32> = (0..pos.len() as u32).collect();
        order.sort_unstable_by(|&a, &b| pos[a as usize].0.total_cmp(&pos[b as usize].0));
        let max2 = max_pair_px * max_pair_px;
        let min2 = MIN_PAIR_PX * MIN_PAIR_PX;
        let mut pairs = Vec::new();
        for (k, &a) in order.iter().enumerate() {
            let (ax, ay) = pos[a as usize];
            for &b in &order[k + 1..] {
                let (bx, by) = pos[b as usize];
                if bx - ax > max_pair_px {
                    break;
                }
                let d2 = (bx - ax) * (bx - ax) + (by - ay) * (by - ay);
                if d2 >= min2 && d2 <= max2 {
                    pairs.push((d2.sqrt() as f32, a, b));
                }
            }
        }
        pairs.sort_unstable_by(|p, q| p.0.total_cmp(&q.0));
        let pair_len = pairs.iter().map(|p| p.0).collect();
        let pair_vec = pairs
            .iter()
            .map(|&(_, a, b)| {
                let ((ax, ay), (bx, by)) = (pos[a as usize], pos[b as usize]);
                [ax as f32, ay as f32, (bx - ax) as f32, (by - ay) as f32]
            })
            .collect();
        let pair_idx = pairs.iter().map(|&(_, a, b)| (a, b)).collect();
        Self {
            pos,
            grid,
            near,
            pair_len,
            pair_vec,
            pair_idx,
            probe_tol,
        }
    }

    /// Number of indexed stars.
    #[must_use]
    pub fn len(&self) -> usize {
        self.pos.len()
    }

    /// Whether no star is indexed.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.pos.is_empty()
    }

    /// Number of indexed pairs.
    #[must_use]
    pub fn n_pairs(&self) -> usize {
        self.pair_len.len()
    }

    fn hit(&self, x: f64, y: f64, tol: f64) -> Option<usize> {
        self.grid.nearest(&self.pos, x, y, tol)
    }
}

/// What the search is looking for.
#[derive(Debug, Clone)]
pub struct SeedParams {
    /// Expected pixel scale, catalogue units (arcsec) per pixel.
    pub scale: f64,
    /// Fractional tolerance on the scale.
    pub scale_tol: f64,
    /// Image width and height, pixels.
    pub width: f64,
    /// Image height, pixels.
    pub height: f64,
    /// Catalogue stars that seed quads (the brightest).
    pub seed_stars: usize,
    /// Most catalogue quads examined.
    pub max_quads: usize,
    /// Catalogue stars (the brightest) a transform is scored on.
    pub census_stars: usize,
    /// Fewest of them a transform must put on a detection, when at least
    /// `1.5 × min_census` of them land in the frame; fewer in frame lower it to
    /// two thirds of those, but never below 4. In a dense frame it is raised to
    /// [`CENSUS_SIGNIFICANCE`] times the hits expected by chance.
    pub min_census: usize,
    /// Work charged to the budget for each candidate handed to the verification.
    pub verify_cost: u64,
}

/// A candidate found by [`search`]: a plate (pixel → catalogue coordinates) and the
/// star pairs its census found.
#[derive(Debug, Clone)]
pub struct Candidate {
    /// Plate from the similarity transform, refitted to the census pairs.
    pub plate: PlateConstants,
    /// Image positions of the census pairs.
    pub img: Vec<(f64, f64)>,
    /// Their catalogue positions.
    pub cat: Vec<(f64, f64)>,
}

/// A similarity transform from catalogue coordinates to pixels:
/// `z_pix = s · w + t`, with `w = z_cat` or its mirror image `conj(z_cat)`.
#[derive(Debug, Clone, Copy)]
struct Similarity {
    sr: f64,
    si: f64,
    tr: f64,
    ti: f64,
    mirrored: bool,
}

impl Similarity {
    /// The transform taking catalogue `p1`, `p2` to pixels `a`, `b`. (The search
    /// computes the same inline, with the per-quad terms hoisted.)
    #[cfg(test)]
    fn from_pair(
        p1: (f64, f64),
        p2: (f64, f64),
        a: (f64, f64),
        b: (f64, f64),
        mirrored: bool,
    ) -> Option<Self> {
        let flip = |p: (f64, f64)| if mirrored { (p.0, -p.1) } else { p };
        let (w1, w2) = (flip(p1), flip(p2));
        let (dwr, dwi) = (w2.0 - w1.0, w2.1 - w1.1);
        let den = dwr * dwr + dwi * dwi;
        if den <= 0.0 {
            return None;
        }
        let (dzr, dzi) = (b.0 - a.0, b.1 - a.1);
        // s = dz / dw
        let sr = (dzr * dwr + dzi * dwi) / den;
        let si = (dzi * dwr - dzr * dwi) / den;
        let tr = a.0 - (sr * w1.0 - si * w1.1);
        let ti = a.1 - (sr * w1.1 + si * w1.0);
        Some(Self {
            sr,
            si,
            tr,
            ti,
            mirrored,
        })
    }

    #[inline]
    fn apply(&self, p: (f64, f64)) -> (f64, f64) {
        let wy = if self.mirrored { -p.1 } else { p.1 };
        (
            self.sr * p.0 - self.si * wy + self.tr,
            self.sr * wy + self.si * p.0 + self.ti,
        )
    }

    /// The inverse, as plate constants (pixel → catalogue).
    fn plate(&self) -> Option<PlateConstants> {
        // pixel = M · (x, wy) + t with M = [[sr, -si], [si, sr]]; cat y = ±wy.
        let det = self.sr * self.sr + self.si * self.si;
        if det <= 0.0 {
            return None;
        }
        let (ir, ii) = (self.sr / det, -self.si / det); // 1/s
        // w = (z - t) / s
        let (a, b, c) = (ir, -ii, -(ir * self.tr - ii * self.ti));
        let (d, e, f) = (ii, ir, -(ir * self.ti + ii * self.tr));
        Some(if self.mirrored {
            PlateConstants {
                a,
                b,
                c,
                d: -d,
                e: -e,
                f: -f,
            }
        } else {
            PlateConstants { a, b, c, d, e, f }
        })
    }
}

/// Catalogue quads: each of the brightest `seed_stars` with three of its four nearest
/// neighbours among them, as indexes into `cat`, the widest pair first.
fn seed_quads(cat: &[(f64, f64)], seed_stars: usize, max_quads: usize) -> Vec<[usize; 4]> {
    let n = cat.len().min(seed_stars);
    let mut quads = Vec::new();
    for a in 0..n {
        let mut near: Vec<(f64, usize)> = (0..n)
            .filter(|&b| b != a)
            .map(|b| {
                let d = (cat[b].0 - cat[a].0).hypot(cat[b].1 - cat[a].1);
                (d, b)
            })
            .collect();
        near.sort_unstable_by(|p, q| p.0.total_cmp(&q.0));
        near.truncate(4);
        if near.len() < 3 {
            continue;
        }
        for skip in 0..near.len() {
            let mut q = [a; 4];
            let mut k = 1;
            for (m, &(_, b)) in near.iter().enumerate() {
                if m != skip && k < 4 {
                    q[k] = b;
                    k += 1;
                }
            }
            if k < 4 {
                continue;
            }
            // Widest pair first: it fixes the transform best.
            let mut widest = (0, 1, -1.0);
            for i in 0..4 {
                for j in i + 1..4 {
                    let d = (cat[q[i]].0 - cat[q[j]].0).hypot(cat[q[i]].1 - cat[q[j]].1);
                    if d > widest.2 {
                        widest = (i, j, d);
                    }
                }
            }
            let rest: Vec<usize> = (0..4).filter(|&k| k != widest.0 && k != widest.1).collect();
            quads.push([q[widest.0], q[widest.1], q[rest[0]], q[rest[1]]]);
            if quads.len() >= max_quads {
                break;
            }
        }
        if quads.len() >= max_quads {
            break;
        }
    }
    quads
}

/// The longest pair, in pixels, any seed quad can need: for sizing the image's
/// pair table.
#[must_use]
pub fn max_backbone_px(cat: &[(f64, f64)], p: &SeedParams) -> f64 {
    let quads = seed_quads(cat, p.seed_stars, p.max_quads);
    let longest = quads
        .iter()
        .map(|q| (cat[q[0]].0 - cat[q[1]].0).hypot(cat[q[0]].1 - cat[q[1]].1))
        .fold(0.0, f64::max);
    longest / (p.scale * (1.0 - p.scale_tol)) + 1.0
}

/// Search for a transform that puts the catalogue's bright stars on detections.
///
/// `cat` holds catalogue positions in standard coordinates (arcsec), brightest first.
/// Every transform whose census passes is passed to `accept` (the caller's
/// verification); the first it accepts is returned. Each probe — one transform tried
/// on a quad's third star — and each census star scored costs one unit of `budget`;
/// the search stops when it runs out.
pub fn search(
    index: &ImageIndex,
    cat: &[(f64, f64)],
    p: &SeedParams,
    budget: &mut u64,
    mut accept: impl FnMut(&Candidate) -> bool,
) -> Option<Candidate> {
    if index.len() < 4 || cat.len() < 4 {
        return None;
    }
    let tol = index.probe_tol;
    let census: Vec<(f64, f64)> = cat.iter().take(p.census_stars).copied().collect();
    let margin = 1.5 * tol;
    let density = index.len() as f64 / (p.width * p.height).max(1.0);
    let p_chance = 1.0 - (-density * core::f64::consts::PI * margin * margin).exp();
    let in_frame = |q: (f64, f64)| {
        q.0 >= -margin && q.1 >= -margin && q.0 < p.width + margin && q.1 < p.height + margin
    };
    // Transforms already scored: the same image stars reached from another quad.
    let mut seen: std::collections::HashSet<(i64, i64, i64, bool)> = Default::default();

    for quad in seed_quads(cat, p.seed_stars, p.max_quads) {
        let (p1, p2, p3, p4) = (cat[quad[0]], cat[quad[1]], cat[quad[2]], cat[quad[3]]);
        let backbone = (p2.0 - p1.0).hypot(p2.1 - p1.1);
        let lo = (backbone / (p.scale * (1.0 + p.scale_tol)) - tol) as f32;
        let hi = (backbone / (p.scale * (1.0 - p.scale_tol)) + tol) as f32;
        let from = index.pair_len.partition_point(|&l| l < lo);
        let to = index.pair_len.partition_point(|&l| l <= hi);
        for mirrored in [false, true] {
            // Per quad and parity: with w = z_cat (or its mirror image), the
            // transform taking w1, w2 to pixels a, b has s = (b - a) / (w2 - w1),
            // and puts w3 at a + s (w3 - w1).
            let flip = |q: (f64, f64)| if mirrored { (q.0, -q.1) } else { q };
            let (w1, w2, w3, w4) = (flip(p1), flip(p2), flip(p3), flip(p4));
            let (dwr, dwi) = (w2.0 - w1.0, w2.1 - w1.1);
            let den = dwr * dwr + dwi * dwi;
            if den <= 0.0 {
                continue;
            }
            let (ir, ii) = (dwr / den, -dwi / den); // 1 / (w2 - w1)
            let (d3, d4) = ((w3.0 - w1.0, w3.1 - w1.1), (w4.0 - w1.0, w4.1 - w1.1));
            // The same in f32 for the first test, which rejects almost every probe.
            let (ir32, ii32) = (ir as f32, ii as f32);
            let (d3x, d3y, d4x, d4y) = (d3.0 as f32, d3.1 as f32, d4.0 as f32, d4.1 as f32);
            for (k, &[ax, ay, dx, dy]) in index.pair_vec[from..to].iter().enumerate() {
                // From a to b, s = d / (w2 - w1) puts w3 at a + s (w3 - w1); from b
                // to a, s is negated and w3 lands at b - s (w3 - w1).
                let (sr, si) = (dx * ir32 - dy * ii32, dx * ii32 + dy * ir32);
                let (e3x, e3y) = (sr * d3x - si * d3y, sr * d3y + si * d3x);
                let (e4x, e4y) = (sr * d4x - si * d4y, sr * d4y + si * d4x);
                let (bx, by) = (ax + dx, ay + dy);
                for (forward, zx, zy, sign) in [(true, ax, ay, 1.0f32), (false, bx, by, -1.0)] {
                    if *budget == 0 {
                        return None;
                    }
                    *budget -= 1;
                    if !(index.near.maybe(zx + sign * e3x, zy + sign * e3y)
                        && index.near.maybe(zx + sign * e4x, zy + sign * e4y))
                    {
                        continue;
                    }
                    // Exactly, in f64, from the stars themselves.
                    let (i, j) = index.pair_idx[from + k];
                    let (a, b) = (index.pos[i as usize], index.pos[j as usize]);
                    let (za, zb) = if forward { (a, b) } else { (b, a) };
                    let (dzr, dzi) = (zb.0 - za.0, zb.1 - za.1);
                    let (sr, si) = (dzr * ir - dzi * ii, dzr * ii + dzi * ir);
                    let q3 = (za.0 + sr * d3.0 - si * d3.1, za.1 + sr * d3.1 + si * d3.0);
                    let q4 = (za.0 + sr * d4.0 - si * d4.1, za.1 + sr * d4.1 + si * d4.0);
                    if index.hit(q3.0, q3.1, tol).is_none() || index.hit(q4.0, q4.1, tol).is_none()
                    {
                        continue;
                    }
                    let t = Similarity {
                        sr,
                        si,
                        tr: za.0 - (sr * w1.0 - si * w1.1),
                        ti: za.1 - (sr * w1.1 + si * w1.0),
                        mirrored,
                    };
                    let key = (
                        (t.tr / 2.0).round() as i64,
                        (t.ti / 2.0).round() as i64,
                        (t.si.atan2(t.sr) * 200.0).round() as i64,
                        mirrored,
                    );
                    if !seen.insert(key) {
                        continue;
                    }
                    // Census of the brightest catalogue stars in the frame.
                    *budget = budget.saturating_sub(census.len() as u64);
                    let mut n_in = 0usize;
                    let mut img = Vec::new();
                    let mut catp = Vec::new();
                    let mut used = std::collections::HashSet::new();
                    for &c in &census {
                        let q = t.apply(c);
                        if !in_frame(q) {
                            continue;
                        }
                        n_in += 1;
                        if let Some(m) = index.hit(q.0, q.1, margin)
                            && used.insert(m)
                        {
                            img.push(index.pos[m]);
                            catp.push(c);
                        }
                    }
                    // The quad's own four stars are hits by construction.
                    let need = (p.min_census.min((n_in * 2 / 3).max(4)) as f64)
                        .max(CENSUS_SIGNIFICANCE * p_chance * n_in.saturating_sub(4) as f64 + 4.0);
                    if (img.len() as f64) < need {
                        continue;
                    }
                    let plate = crate::math::lsq::fit_affine(&img, &catp)
                        .ok()
                        .or_else(|| t.plate());
                    let Some(plate) = plate else { continue };
                    let cand = Candidate {
                        plate,
                        img,
                        cat: catp,
                    };
                    *budget = budget.saturating_sub(p.verify_cost);
                    if accept(&cand) {
                        return Some(cand);
                    }
                }
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (*seed >> 11) as f64 / (1u64 << 53) as f64
    }

    #[test]
    fn similarity_maps_the_pair_and_inverts_to_a_plate() {
        for mirrored in [false, true] {
            let t = Similarity::from_pair(
                (10.0, 20.0),
                (110.0, -5.0),
                (300.0, 400.0),
                (350.0, 480.0),
                mirrored,
            )
            .unwrap();
            let a = t.apply((10.0, 20.0));
            let b = t.apply((110.0, -5.0));
            assert!((a.0 - 300.0).abs() < 1e-9 && (a.1 - 400.0).abs() < 1e-9);
            assert!((b.0 - 350.0).abs() < 1e-9 && (b.1 - 480.0).abs() < 1e-9);
            let pl = t.plate().unwrap();
            let q = (37.0, -12.0);
            let z = t.apply(q);
            let back = (
                pl.a * z.0 + pl.b * z.1 + pl.c,
                pl.d * z.0 + pl.e * z.1 + pl.f,
            );
            assert!((back.0 - q.0).abs() < 1e-9 && (back.1 - q.1).abs() < 1e-9);
            // A mirrored transform has a negative determinant.
            let det = pl.a * pl.e - pl.b * pl.d;
            assert_eq!(det < 0.0, mirrored);
        }
    }

    /// A field whose image list is shuffled in brightness and mostly unrelated to
    /// the catalogue still yields the true transform.
    #[test]
    fn finds_the_transform_whatever_the_image_ranking() {
        let mut seed = 7u64;
        let (w, h) = (1000.0, 800.0);
        let scale = 2.0; // arcsec per pixel
        let rot = 0.6f64;
        let (c, s) = (rot.cos(), rot.sin());
        // Catalogue: 200 stars over the field, brightest first.
        let cat: Vec<(f64, f64)> = (0..200)
            .map(|_| {
                (
                    (lcg(&mut seed) - 0.5) * w * scale,
                    (lcg(&mut seed) - 0.5) * h * scale,
                )
            })
            .collect();
        let to_pix = |p: (f64, f64)| {
            let (x, y) = (p.0 / scale, p.1 / scale);
            (c * x - s * y + w / 2.0, s * x + c * y + h / 2.0)
        };
        // Image: a third of the catalogue (every third star), with 1000 unrelated
        // detections, in no particular order.
        let mut pos: Vec<(f64, f64)> = cat.iter().step_by(3).map(|&p| to_pix(p)).collect();
        for _ in 0..1000 {
            pos.push((lcg(&mut seed) * w, lcg(&mut seed) * h));
        }
        for i in (1..pos.len()).rev() {
            let j = (lcg(&mut seed) * (i + 1) as f64) as usize;
            pos.swap(i, j);
        }
        let p = SeedParams {
            scale,
            scale_tol: 0.05,
            width: w,
            height: h,
            seed_stars: 100,
            max_quads: 400,
            census_stars: 100,
            min_census: 10,
            verify_cost: 0,
        };
        let index = ImageIndex::new(pos, max_backbone_px(&cat, &p), 2.5);
        let mut budget = 50_000_000u64;
        // Stand-in for the solver's verification: a third of the 100 brightest
        // catalogue stars are in the image, and the true transform finds them.
        let found = search(&index, &cat, &p, &mut budget, |c| c.img.len() >= 25).expect("found");
        let pl = &found.plate;
        // Check the plate against the truth at the frame corners.
        for &(x, y) in &[(0.0, 0.0), (w, 0.0), (0.0, h), (w, h)] {
            let got = (pl.a * x + pl.b * y + pl.c, pl.d * x + pl.e * y + pl.f);
            let (dx, dy) = (x - w / 2.0, y - h / 2.0);
            let want = ((c * dx + s * dy) * scale, (-s * dx + c * dy) * scale);
            assert!(
                (got.0 - want.0).hypot(got.1 - want.1) < 2.0 * scale,
                "corner ({x},{y}) {got:?} vs {want:?}"
            );
        }
    }

    #[test]
    fn the_budget_bounds_a_search_that_cannot_succeed() {
        let mut seed = 11u64;
        let cat: Vec<(f64, f64)> = (0..200)
            .map(|_| (lcg(&mut seed) * 2000.0, lcg(&mut seed) * 2000.0))
            .collect();
        let pos: Vec<(f64, f64)> = (0..800)
            .map(|_| (lcg(&mut seed) * 1000.0, lcg(&mut seed) * 1000.0))
            .collect();
        let p = SeedParams {
            scale: 2.0,
            scale_tol: 0.1,
            width: 1000.0,
            height: 1000.0,
            seed_stars: 100,
            max_quads: 400,
            census_stars: 100,
            min_census: 10,
            verify_cost: 0,
        };
        let index = ImageIndex::new(pos, max_backbone_px(&cat, &p), 2.5);
        let mut budget = 1_000u64;
        assert!(search(&index, &cat, &p, &mut budget, |_| false).is_none());
        assert_eq!(budget, 0, "a refused search runs until the budget is spent");
    }
}
