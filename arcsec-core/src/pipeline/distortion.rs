//! Distortion-aware matching: a polynomial plate, fitted by re-matching the whole
//! catalogue against the image as the model improves (astrometry.net's "tweak").
//!
//! The spiral search verifies a *linear* plate: catalogue stars are mapped into the
//! image and paired with detections within 6, 3 and finally 2 pixels. Where the
//! optics distort the field by more than that, only the stars in the region the
//! linear plate happens to fit survive, and the refit is made to that region alone.
//! On a camera lens or a 12° TESS frame that leaves the corners hundreds or
//! thousands of arcseconds out, although the field is the right one.
//!
//! [`refine`] fixes that. Starting from the verified plate it pairs every catalogue
//! star with a detection within a radius the image's star density allows, fits a
//! polynomial plate (linear, quadratic or cubic, the order chosen by F-tests and by
//! how much of the frame the pairs cover) with outliers clipped, and repeats with
//! the improved model, so that pairs spread outwards from the region the plate
//! first fitted. The final pairs are those within the tight verification radius of
//! the model.
//!
//! The model is internal. It decides which stars match and which plate is reported,
//! but the WCS written by default stays linear, as ASTAP's is: see
//! [`best_linear`].

use crate::types::{PairedPositions, PlateConstants, StarList};

/// An image position (pixels) and the standard coordinates (arcsec) it is paired with.
pub(super) type Pair = ((f64, f64), (f64, f64));

/// A uniform grid over the detected stars, for nearest-neighbour lookups.
pub(super) struct StarGrid<'a> {
    stars: &'a StarList,
    pub(super) min_x: f64,
    pub(super) min_y: f64,
    pub(super) max_x: f64,
    pub(super) max_y: f64,
    cell: f64,
    nx: usize,
    ny: usize,
    cells: Vec<Vec<u32>>,
}

impl<'a> StarGrid<'a> {
    /// Bucket `stars` into square cells of side `cell` pixels. `None` if the stars
    /// do not span a rectangle (fewer than two distinct x or y values).
    pub(super) fn new(stars: &'a StarList, cell: f64) -> Option<Self> {
        let (mut min_x, mut min_y) = (f64::INFINITY, f64::INFINITY);
        let (mut max_x, mut max_y) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
        for st in &stars.0 {
            min_x = min_x.min(st.x);
            max_x = max_x.max(st.x);
            min_y = min_y.min(st.y);
            max_y = max_y.max(st.y);
        }
        if !(min_x.is_finite() && min_y.is_finite() && max_x > min_x && max_y > min_y) {
            return None;
        }
        let cell = cell.max(1.0);
        let nx = (((max_x - min_x) / cell).ceil() as usize + 1).max(1);
        let ny = (((max_y - min_y) / cell).ceil() as usize + 1).max(1);
        let mut cells: Vec<Vec<u32>> = vec![Vec::new(); nx * ny];
        for (i, st) in stars.0.iter().enumerate() {
            let gx = ((st.x - min_x) / cell) as usize;
            let gy = ((st.y - min_y) / cell) as usize;
            cells[gy.min(ny - 1) * nx + gx.min(nx - 1)].push(i as u32);
        }
        Some(Self {
            stars,
            min_x,
            min_y,
            max_x,
            max_y,
            cell,
            nx,
            ny,
            cells,
        })
    }

    /// Number of stars in the grid.
    pub(super) fn len(&self) -> usize {
        self.stars.len()
    }

    /// Position of star `i`.
    pub(super) fn pos(&self, i: usize) -> (f64, f64) {
        let s = &self.stars.0[i];
        (s.x, s.y)
    }

    /// Whether `(px, py)` lies within `margin` of the stars' bounding box.
    pub(super) fn near(&self, px: f64, py: f64, margin: f64) -> bool {
        px >= self.min_x - margin
            && px <= self.max_x + margin
            && py >= self.min_y - margin
            && py <= self.max_y + margin
    }

    /// The nearest star not yet `used` within `√r2` of `(px, py)`, scanning cells in
    /// row order and keeping the first of equally near stars.
    pub(super) fn nearest(&self, px: f64, py: f64, r2: f64, used: &[bool]) -> Option<usize> {
        let span = (r2.sqrt() / self.cell).ceil().max(1.0) as isize;
        let gx = (((px - self.min_x) / self.cell) as isize).clamp(0, self.nx as isize - 1);
        let gy = (((py - self.min_y) / self.cell) as isize).clamp(0, self.ny as isize - 1);
        let mut best_i: Option<usize> = None;
        let mut best_d2 = r2;
        for oy in -span..=span {
            for ox in -span..=span {
                let cx = gx + ox;
                let cy = gy + oy;
                if cx < 0 || cy < 0 || cx >= self.nx as isize || cy >= self.ny as isize {
                    continue;
                }
                for &i in &self.cells[cy as usize * self.nx + cx as usize] {
                    let i = i as usize;
                    if used[i] {
                        continue;
                    }
                    let st = &self.stars.0[i];
                    let d2 = (st.x - px) * (st.x - px) + (st.y - py) * (st.y - py);
                    if d2 < best_d2 {
                        best_d2 = d2;
                        best_i = Some(i);
                    }
                }
            }
        }
        best_i
    }
}

/// A polynomial plate: pixels → standard coordinates (arcsec), each axis a
/// polynomial in `u = (x - cx) / s`, `v = (y - cy) / s` over the first `n_terms`
/// of [`SIP_TERMS`](crate::wcs::sip::SIP_TERMS) (3: linear, 6: quadratic, 10: cubic).
#[derive(Debug, Clone, PartialEq)]
pub(super) struct PolyPlate {
    cx: f64,
    cy: f64,
    s: f64,
    pub(super) n_terms: usize,
    kx: [f64; 10],
    ky: [f64; 10],
}

/// The monomials `uᵖ vᵠ` of [`SIP_TERMS`](crate::wcs::sip::SIP_TERMS) at `(u, v)`.
fn monomials(u: f64, v: f64) -> [f64; 10] {
    let (u2, v2) = (u * u, v * v);
    [1.0, u, v, u2, u * v, v2, u2 * u, u2 * v, u * v2, v2 * v]
}

impl PolyPlate {
    /// The linear plate `p` in this representation.
    pub(super) fn linear(p: &PlateConstants, cx: f64, cy: f64, s: f64) -> Self {
        let mut kx = [0.0; 10];
        let mut ky = [0.0; 10];
        kx[0] = p.a * cx + p.b * cy + p.c;
        kx[1] = p.a * s;
        kx[2] = p.b * s;
        ky[0] = p.d * cx + p.e * cy + p.f;
        ky[1] = p.d * s;
        ky[2] = p.e * s;
        Self {
            cx,
            cy,
            s,
            n_terms: 3,
            kx,
            ky,
        }
    }

    /// Pixel → standard coordinates (arcsec).
    pub(super) fn apply(&self, x: f64, y: f64) -> (f64, f64) {
        let m = monomials((x - self.cx) / self.s, (y - self.cy) / self.s);
        let mut xi = 0.0;
        let mut eta = 0.0;
        for ((kx, ky), m) in self.kx.iter().zip(&self.ky).zip(&m).take(self.n_terms) {
            xi += kx * m;
            eta += ky * m;
        }
        (xi, eta)
    }

    /// The linear part as plate constants (pixels → arcsec): the plate the
    /// polynomial reduces to at the frame centre, constant included.
    pub(super) fn linear_part(&self) -> PlateConstants {
        let (a, b) = (self.kx[1] / self.s, self.kx[2] / self.s);
        let (d, e) = (self.ky[1] / self.s, self.ky[2] / self.s);
        PlateConstants {
            a,
            b,
            c: self.kx[0] - a * self.cx - b * self.cy,
            d,
            e,
            f: self.ky[0] - d * self.cx - e * self.cy,
        }
    }

    /// Standard coordinates → pixel, by fixed-point iteration on the linear part.
    /// `None` if the linear part is singular or the iteration does not settle.
    pub(super) fn to_pixel(&self, xi: f64, eta: f64) -> Option<(f64, f64)> {
        let l = self.linear_part();
        let det = l.a * l.e - l.b * l.d;
        if det.abs() < 1e-12 || !det.is_finite() {
            return None;
        }
        let inv = |dx: f64, dy: f64| ((l.e * dx - l.b * dy) / det, (-l.d * dx + l.a * dy) / det);
        let (mut x, mut y) = inv(xi - l.c, eta - l.f);
        if self.n_terms <= 3 {
            return Some((x, y));
        }
        for _ in 0..12 {
            let (fx, fy) = self.apply(x, y);
            let (dx, dy) = inv(xi - fx, eta - fy);
            x += dx;
            y += dy;
            if dx.abs() + dy.abs() < 1e-4 {
                return Some((x, y));
            }
        }
        None
    }

    /// The pixel scale at the centre, arcsec per pixel.
    pub(super) fn scale(&self) -> f64 {
        let l = self.linear_part();
        (l.a * l.e - l.b * l.d).abs().sqrt()
    }
}

/// Least-squares polynomial plate over the first `n_terms` of [`SIP_TERMS`](crate::wcs::sip::SIP_TERMS).
fn fit_poly(
    img: &[(f64, f64)],
    cat: &[(f64, f64)],
    n_terms: usize,
    frame: (f64, f64, f64),
) -> Option<PolyPlate> {
    let (cx, cy, s) = frame;
    if img.len() < n_terms + 2 {
        return None;
    }
    // Normal equations, solved by Cholesky: with `u`, `v` scaled to about ±1 the
    // cubic's monomials are well conditioned, and this is several times cheaper
    // than the Givens solver over a few hundred pairs, fitted many times a solve.
    let mut ata = [[0.0f64; 10]; 10];
    let mut atb = [[0.0f64; 10]; 2];
    for (&(x, y), &(xi, eta)) in img.iter().zip(cat) {
        let m = monomials((x - cx) / s, (y - cy) / s);
        for i in 0..n_terms {
            atb[0][i] += m[i] * xi;
            atb[1][i] += m[i] * eta;
            for j in 0..=i {
                ata[i][j] += m[i] * m[j];
            }
        }
    }
    let [kx, ky] = cholesky_solve(&ata, &atb, n_terms)?;
    if !kx.iter().chain(&ky).all(|v| v.is_finite()) {
        return None;
    }
    Some(PolyPlate {
        cx,
        cy,
        s,
        n_terms,
        kx,
        ky,
    })
}

/// Solve `A k = b` for both right-hand sides, `A` symmetric positive definite and
/// given by its lower triangle (first `n` rows and columns). `None` if it is not
/// positive definite to working precision.
#[allow(clippy::needless_range_loop)] // index arithmetic is the algorithm
fn cholesky_solve(a: &[[f64; 10]; 10], b: &[[f64; 10]; 2], n: usize) -> Option<[[f64; 10]; 2]> {
    let mut l = [[0.0f64; 10]; 10];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i][j];
            for k in 0..j {
                sum -= l[i][k] * l[j][k];
            }
            if i == j {
                if sum.is_nan() || sum <= 1e-12 * a[i][i].abs().max(f64::MIN_POSITIVE) {
                    return None;
                }
                l[i][i] = sum.sqrt();
            } else {
                l[i][j] = sum / l[j][j];
            }
        }
    }
    let mut out = [[0.0f64; 10]; 2];
    for (rhs, k) in b.iter().zip(out.iter_mut()) {
        let mut y = [0.0f64; 10];
        for i in 0..n {
            let mut sum = rhs[i];
            for j in 0..i {
                sum -= l[i][j] * y[j];
            }
            y[i] = sum / l[i][i];
        }
        for i in (0..n).rev() {
            let mut sum = y[i];
            for j in (i + 1)..n {
                sum -= l[j][i] * k[j];
            }
            k[i] = sum / l[i][i];
        }
    }
    Some(out)
}

/// Residual sum of squares (arcsec²) of the pairs under `m`.
fn rss(m: &PolyPlate, img: &[(f64, f64)], cat: &[(f64, f64)]) -> f64 {
    img.iter()
        .zip(cat)
        .map(|(&(x, y), &(xi, eta))| {
            let (px, py) = m.apply(x, y);
            (px - xi) * (px - xi) + (py - eta) * (py - eta)
        })
        .sum()
}

/// Residuals (arcsec) of every pair under `m`.
fn residuals(m: &PolyPlate, img: &[(f64, f64)], cat: &[(f64, f64)]) -> Vec<f64> {
    img.iter()
        .zip(cat)
        .map(|(&(x, y), &(xi, eta))| {
            let (px, py) = m.apply(x, y);
            ((px - xi) * (px - xi) + (py - eta) * (py - eta)).sqrt()
        })
        .collect()
}

/// A polynomial fit with outliers clipped: the model, which pairs it kept, and its
/// residual sum of squares over them (arcsec²).
struct RobustFit {
    model: PolyPlate,
    keep: Vec<bool>,
    rss: f64,
    n: usize,
}

/// Fit `n_terms`, then drop pairs further than `max(floor, 3σ)` from the model
/// (σ from the median residual, so a minority of wrong pairs cannot inflate it)
/// and refit, until the kept set stops changing.
fn robust_fit(
    img: &[(f64, f64)],
    cat: &[(f64, f64)],
    n_terms: usize,
    frame: (f64, f64, f64),
    floor: f64,
) -> Option<RobustFit> {
    let mut keep = vec![true; img.len()];
    let mut model = None;
    for _ in 0..6 {
        let (ki, kc): (Vec<_>, Vec<_>) = img
            .iter()
            .zip(cat)
            .zip(&keep)
            .filter(|&(_, &k)| k)
            .map(|((&i, &c), _)| (i, c))
            .unzip();
        let m = fit_poly(&ki, &kc, n_terms, frame)?;
        let r = residuals(&m, img, cat);
        let mut kept: Vec<f64> = r
            .iter()
            .zip(&keep)
            .filter(|&(_, &k)| k)
            .map(|(&v, _)| v)
            .collect();
        let mid = kept.len() / 2;
        let sigma = 1.4826 * *kept.select_nth_unstable_by(mid, f64::total_cmp).1;
        let threshold = (3.0 * sigma).max(floor);
        let next: Vec<bool> = r.iter().map(|&v| v <= threshold).collect();
        model = Some(m);
        if next == keep {
            break;
        }
        keep = next;
    }
    let model = model?;
    let r = residuals(&model, img, cat);
    let (mut rss, mut n) = (0.0, 0);
    for (&v, &k) in r.iter().zip(&keep) {
        if k {
            rss += v * v;
            n += 1;
        }
    }
    (n >= n_terms + 2).then_some(RobustFit {
        model,
        keep,
        rss,
        n,
    })
}

/// Cells of a `GRID` × `GRID` division of the frame holding at least
/// `COVERAGE_MIN` pairs.
fn cells_covered(img: &[(f64, f64)], keep: &[bool], w: f64, h: f64) -> usize {
    const GRID: usize = 3;
    const COVERAGE_MIN: usize = 3;
    let mut cells = [0usize; GRID * GRID];
    let cell =
        |v: f64, size: f64| ((v / size * GRID as f64).floor().max(0.0) as usize).min(GRID - 1);
    for (&(x, y), &k) in img.iter().zip(keep) {
        if k {
            cells[cell(y, h) * GRID + cell(x, w)] += 1;
        }
    }
    cells.iter().filter(|&&n| n >= COVERAGE_MIN).count()
}

/// The F statistic of `hi` (more terms) over `lo` on the same pairs.
fn f_stat(lo_rss: f64, hi_rss: f64, lo_terms: usize, hi_terms: usize, n: usize) -> f64 {
    let extra = 2.0 * (hi_terms - lo_terms) as f64;
    let dof = (2 * n).saturating_sub(2 * hi_terms).max(1) as f64;
    if hi_rss <= 0.0 {
        return if lo_rss > 0.0 { f64::INFINITY } else { 0.0 };
    }
    ((lo_rss - hi_rss) / extra) / (hi_rss / dof)
}

/// The F statistic a higher-order term set must reach over the next lower one.
pub(super) const MIN_DISTORTION_F: f64 = 4.0;

/// Chance-match budget for the wide matching radius: on average at most this many
/// unrelated detections within it of a catalogue star.
const CHANCE_MATCHES: f64 = 0.1;
/// Bounds of the wide matching radius, pixels.
const WIDE_RADIUS_PX: (f64, f64) = (4.0, 24.0);
/// Re-match/refit rounds.
const TWEAK_ROUNDS: usize = 10;

/// The outcome of [`refine`].
pub(super) struct Refined {
    /// The polynomial plate.
    pub(super) model: PolyPlate,
    /// Pairs within the tight radius of the model: detections (pixels)…
    pub(super) img_pos: Vec<(f64, f64)>,
    /// …and their catalogue stars (standard coordinates, arcsec).
    pub(super) cat_pos: Vec<(f64, f64)>,
    /// RMS of the model over those pairs, arcsec.
    pub(super) rms: f64,
    /// F statistic of the model over a linear plate on the same pairs; 0 for a
    /// linear model.
    pub(super) f_linear: f64,
    /// Cells of the 3×3 frame grid the final pairs reach.
    pub(super) cells: usize,
    /// The last round's wide-radius pairs, for [`Refined::unmodelled`].
    wide: PairedPositions,
    /// The plate the model was started from.
    start: PlateConstants,
    /// Centre and half-size of the frame, as the model's `(cx, cy, s)`.
    frame: (f64, f64, f64),
    /// The tight radius, pixels.
    tight_px: f64,
}

impl Refined {
    /// Distortion where the stars are, whether or not the model could be fitted
    /// over the frame: a cubic fitted to the last wide-radius pairs whatever their
    /// coverage. Returns its F statistic over a linear fit to the same pairs, and
    /// the largest distance (pixels) between it and the starting plate at them.
    pub(super) fn unmodelled(&self) -> (f64, f64) {
        let scale = self.model.scale();
        let (img, cat) = (&self.wide.0, &self.wide.1);
        let Some(cubic) = robust_fit(img, cat, 10, self.frame, self.tight_px * scale) else {
            return (0.0, 0.0);
        };
        let (ki, kc): (Vec<_>, Vec<_>) = img
            .iter()
            .zip(cat)
            .zip(&cubic.keep)
            .filter(|&(_, &k)| k)
            .map(|((&i, &c), _)| (i, c))
            .unzip();
        let Some(lin) = fit_poly(&ki, &kc, 3, self.frame) else {
            return (0.0, 0.0);
        };
        let f = f_stat(rss(&lin, &ki, &kc), cubic.rss, 3, 10, cubic.n);
        let p = &self.start;
        let dep = ki
            .iter()
            .map(|&(x, y)| {
                let (mx, my) = cubic.model.apply(x, y);
                let (lx, ly) = (p.a * x + p.b * y + p.c, p.d * x + p.e * y + p.f);
                (mx - lx).hypot(my - ly) / scale
            })
            .fold(0.0, f64::max);
        (f, dep)
    }
}

/// Pair catalogue stars with detections under `model`, one-to-one, nearest first
/// in catalogue order.
fn match_stars(
    grid: &StarGrid<'_>,
    cat: &StarList,
    model: &PolyPlate,
    radius: f64,
) -> PairedPositions {
    let r2 = radius * radius;
    let mut used = vec![false; grid.len()];
    let mut img_pos = Vec::new();
    let mut cat_pos = Vec::new();
    for cs in &cat.0 {
        let Some((px, py)) = model.to_pixel(cs.x, cs.y) else {
            continue;
        };
        if !grid.near(px, py, radius) {
            continue;
        }
        if let Some(i) = grid.nearest(px, py, r2, &used) {
            used[i] = true;
            img_pos.push(grid.pos(i));
            cat_pos.push((cs.x, cs.y));
        }
    }
    (img_pos, cat_pos)
}

/// Fit a distortion model from the linear plate `start` (pixels of a `w` × `h`
/// image → standard coordinates of `cat`), re-matching the catalogue as it improves.
///
/// Returns `None` if too few pairs are found to fit anything.
pub(super) fn refine(
    grid: &StarGrid<'_>,
    cat: &StarList,
    seeds: &[Pair],
    start: &PlateConstants,
    w: usize,
    h: usize,
    tight_px: f64,
) -> Option<Refined> {
    let (wf, hf) = (w as f64, h as f64);
    let frame = ((wf - 1.0) * 0.5, (hf - 1.0) * 0.5, 0.5 * wf.max(hf));
    let mut model = PolyPlate::linear(start, frame.0, frame.1, frame.2);
    let scale = model.scale();
    if !(scale.is_finite() && scale > 0.0) {
        return None;
    }
    let density = grid.len() as f64 / (wf * hf).max(1.0);
    let wide = (CHANCE_MATCHES / (core::f64::consts::PI * density))
        .sqrt()
        .clamp(WIDE_RADIUS_PX.0, WIDE_RADIUS_PX.1);

    let mut last_wide = (Vec::new(), Vec::new());
    for round in 0..TWEAK_ROUNDS {
        let (mut ip, mut cp) = match_stars(grid, cat, &model, wide);
        last_wide = (ip.clone(), cp.clone());
        for &(i, c) in seeds {
            ip.push(i);
            cp.push(c);
        }
        let fit = choose_order(&ip, &cp, frame, (wf, hf), tight_px * scale)?;
        let moved = model_distance(&model, &fit.model, w, h) / scale;
        log::debug!(
            "tweak round {round}: {} pairs within {wide:.1} px, {} kept, {} terms, moved {moved:.2} px",
            ip.len(),
            fit.n,
            fit.model.n_terms,
        );
        model = fit.model;
        if moved < 0.05 {
            break;
        }
    }

    // Final pairs: within the tight radius of the model, the model refitted on
    // them at the same order if that pairs more stars.
    let (mut ip, mut cp) = match_stars(grid, cat, &model, tight_px);
    if let Some(fit) = robust_fit(&ip, &cp, model.n_terms, frame, tight_px * scale) {
        let (ip2, cp2) = match_stars(grid, cat, &fit.model, tight_px);
        if ip2.len() >= ip.len() {
            (model, ip, cp) = (fit.model, ip2, cp2);
        }
    }
    if ip.len() < model.n_terms + 2 {
        return None;
    }
    let model_rss = rss(&model, &ip, &cp);
    let rms = (model_rss / ip.len() as f64).sqrt();
    let f_linear = if model.n_terms > 3 {
        fit_poly(&ip, &cp, 3, frame).map_or(0.0, |lin| {
            f_stat(rss(&lin, &ip, &cp), model_rss, 3, model.n_terms, ip.len())
        })
    } else {
        0.0
    };
    let cells = cells_covered(&ip, &vec![true; ip.len()], wf, hf);
    Some(Refined {
        model,
        img_pos: ip,
        cat_pos: cp,
        rms,
        f_linear,
        cells,
        wide: last_wide,
        start: start.clone(),
        frame,
        tight_px,
    })
}

/// Largest distance (arcsec) between two models over a 5 × 5 grid of the frame.
fn model_distance(a: &PolyPlate, b: &PolyPlate, w: usize, h: usize) -> f64 {
    let (wf, hf) = ((w as f64 - 1.0).max(1.0), (h as f64 - 1.0).max(1.0));
    let mut worst: f64 = 0.0;
    for i in 0..5 {
        for j in 0..5 {
            let (x, y) = (wf * i as f64 / 4.0, hf * j as f64 / 4.0);
            let (ax, ay) = a.apply(x, y);
            let (bx, by) = b.apply(x, y);
            worst = worst.max((ax - bx).hypot(ay - by));
        }
    }
    worst
}

/// Fit the pairs, of a `size.0` × `size.1` image, at the highest order they
/// support, clipping at no less than `floor` arcsec: a quadratic or a cubic only if
/// the kept pairs reach seven of the nine cells of a 3×3 grid over the frame, and
/// each only if it beats the best lower order by [`MIN_DISTORTION_F`].
///
/// Seven cells, not nine, for the cubic: this is the bootstrap, and a cubic from
/// most of the frame is what carries the pairs into the rest. Reporting the model
/// needs all nine (`model_distortion` in the solver).
fn choose_order(
    img: &[(f64, f64)],
    cat: &[(f64, f64)],
    frame: (f64, f64, f64),
    size: (f64, f64),
    floor: f64,
) -> Option<RobustFit> {
    let (w, h) = size;
    let mut best = robust_fit(img, cat, 3, frame, floor)?;
    // A quadratic need not beat the linear plate for a cubic to: radial distortion
    // about the centre is odd, all cubic. So each order is tried against the best
    // so far.
    for n_terms in [6, 10] {
        let Some(hi) = robust_fit(img, cat, n_terms, frame, floor) else {
            continue;
        };
        if cells_covered(img, &hi.keep, w, h) < 7 {
            continue;
        }
        // Compare on the pairs the higher order kept.
        let (ki, kc): (Vec<_>, Vec<_>) = img
            .iter()
            .zip(cat)
            .zip(&hi.keep)
            .filter(|&(_, &k)| k)
            .map(|((&i, &c), _)| (i, c))
            .unzip();
        let Some(lo) = fit_poly(&ki, &kc, best.model.n_terms, frame) else {
            continue;
        };
        let lo_rss = rss(&lo, &ki, &kc);
        let f = f_stat(lo_rss, hi.rss, best.model.n_terms, n_terms, hi.n);
        if f >= MIN_DISTORTION_F {
            best = hi;
        }
    }
    Some(best)
}

/// Grid the frame is sampled on for [`best_linear`] and [`max_departure_px`]:
/// `N` × `N` points, corners included, as the benchmark's linear floor is computed.
const FRAME_GRID: usize = 9;

/// The points of the [`FRAME_GRID`] over a `w` × `h` image, 0-based pixels.
fn frame_grid(w: usize, h: usize) -> impl Iterator<Item = (f64, f64)> {
    let (wf, hf) = ((w as f64 - 1.0).max(1.0), (h as f64 - 1.0).max(1.0));
    let step = (FRAME_GRID - 1) as f64;
    (0..FRAME_GRID).flat_map(move |i| {
        (0..FRAME_GRID).map(move |j| (wf * i as f64 / step, hf * j as f64 / step))
    })
}

/// The linear plate closest over the whole frame to the mapping `sky` (pixels →
/// standard coordinates): a least-squares fit to it on the [`FRAME_GRID`].
///
/// This is what a solver that reports a linear WCS should report for a distorted
/// field. It spreads the distortion over the frame, where a plate fitted to the
/// stars a tight verification radius keeps fits the region they come from and
/// lets the corners go: up to twice as far out as this one. It is also the
/// reference the benchmark scores against (`fitslite.best_linear`).
pub(super) fn best_linear(
    sky: impl Fn(f64, f64) -> (f64, f64),
    w: usize,
    h: usize,
) -> Option<PlateConstants> {
    let (img, cat): (Vec<_>, Vec<_>) = frame_grid(w, h).map(|(x, y)| ((x, y), sky(x, y))).unzip();
    crate::math::lsq::fit_affine(&img, &cat).ok()
}

/// Largest distance (pixels) over the [`FRAME_GRID`] between where `model` and the
/// linear plate `p` put the same pixel on the sky.
pub(super) fn max_departure_px(model: &PolyPlate, p: &PlateConstants, w: usize, h: usize) -> f64 {
    let scale = (p.a * p.e - p.b * p.d).abs().sqrt();
    frame_grid(w, h)
        .map(|(x, y)| {
            let (mx, my) = model.apply(x, y);
            let (lx, ly) = (p.a * x + p.b * y + p.c, p.d * x + p.e * y + p.f);
            (mx - lx).hypot(my - ly) / scale
        })
        .fold(0.0, f64::max)
}
