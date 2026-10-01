//! The geometry shared by the index builder and the index solver: the 4-star
//! pattern descriptor, its quantised hash key, the canonical vertex order that
//! gives a star correspondence, and the tangent-plane projection.
//!
//! The descriptor is the six pairwise distances of a quad, sorted, divided by the
//! largest — the five ratios ASTAP's quads use too (`quads/`). It is invariant to
//! translation, rotation, scale and reflection, so one index entry serves both image
//! parities. Unlike astrometry.net's code space it does not say which star is which,
//! so the vertices are put in a canonical order (ascending total distance to the
//! other three), which both sides compute the same way; near-ties in that order are
//! the one ambiguity, and [`orderings`] enumerates the alternatives on the image
//! side. A wrong correspondence fails the affine shape check, so trying an extra
//! ordering costs a few multiplications, never a false match.
//!
//! **Any change to [`BINS`], [`descriptor`], [`key`] or [`canonical`] changes which
//! entries an index holds or how they hash, and must bump the format version in
//! `format.rs`.**

/// Quantisation bins per descriptor dimension. A bin is 1/128 = 0.0078 wide, about
/// the hinted solver's default quad tolerance (0.007).
pub const BINS: f64 = 128.0;

/// Probe the neighbouring bin when a ratio lies within this distance of a bin edge.
/// Measurement noise can only move a ratio across an edge it is already close to,
/// so only those dimensions need a second probe: typically one to four keys per
/// quad instead of 3⁵.
pub const PROBE_EPS: f64 = 0.0025;

/// Relative difference in total distance below which two vertices count as tied in
/// the canonical order, so both orders are tried.
const TIE_EPS: f64 = 0.015;

/// Pairs of a quad's vertices, in the order the six distances are taken.
const PAIRS: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

/// The canonical vertex order of a quad and its longest edge.
///
/// Vertices are sorted by their total distance to the other three, ascending (ties
/// by original position, so the result is deterministic). Returns the order and the
/// sorted totals, or `None` for a degenerate quad (all points coincident).
#[must_use]
pub fn canonical(p: &[(f64, f64); 4]) -> Option<([usize; 4], [f64; 4], f64)> {
    let mut totals = [0.0f64; 4];
    let mut max_edge = 0.0f64;
    for &(i, j) in &PAIRS {
        let d = (p[i].0 - p[j].0).hypot(p[i].1 - p[j].1);
        totals[i] += d;
        totals[j] += d;
        max_edge = max_edge.max(d);
    }
    if max_edge <= 0.0 || !max_edge.is_finite() {
        return None;
    }
    let mut order = [0usize, 1, 2, 3];
    order.sort_by(|&a, &b| totals[a].total_cmp(&totals[b]).then(a.cmp(&b)));
    let sorted = order.map(|i| totals[i]);
    Some((order, sorted, max_edge))
}

/// The canonical order plus the alternatives that a near-tie makes possible: each
/// adjacent pair whose totals differ by less than `TIE_EPS` (1.5 %) may be swapped. At
/// most eight orders (three independent adjacent swaps).
#[must_use]
pub fn orderings(order: [usize; 4], totals: [f64; 4]) -> Vec<[usize; 4]> {
    let mut out = vec![order];
    for k in 0..3 {
        let tied = (totals[k + 1] - totals[k]) <= TIE_EPS * totals[k + 1];
        if tied {
            let n = out.len();
            for i in 0..n {
                let mut o = out[i];
                o.swap(k, k + 1);
                if !out.contains(&o) {
                    out.push(o);
                }
            }
        }
    }
    out
}

/// The five sorted distance ratios of a quad, each in `(0, 1]`.
#[must_use]
pub fn descriptor(p: &[(f64, f64); 4]) -> [f64; 5] {
    let mut e = [0.0f64; 6];
    for (k, &(i, j)) in PAIRS.iter().enumerate() {
        e[k] = (p[i].0 - p[j].0).hypot(p[i].1 - p[j].1);
    }
    e.sort_by(f64::total_cmp);
    let m = e[5].max(1e-300);
    [e[0] / m, e[1] / m, e[2] / m, e[3] / m, e[4] / m]
}

#[inline]
fn bin(v: f64) -> i64 {
    ((v * BINS) as i64).clamp(0, BINS as i64 - 1)
}

/// The descriptor's home bucket: each ratio floored into one of [`BINS`] bins,
/// packed eight bits per dimension, first ratio most significant.
#[must_use]
pub fn key(d: &[f64; 5]) -> u64 {
    d.iter().fold(0u64, |k, &v| (k << 8) | bin(v) as u64)
}

/// Every key a measured descriptor could have had in the index: the home bucket,
/// plus the neighbouring bin in any dimension within [`PROBE_EPS`] of an edge.
/// Written into `out`, which is cleared first.
pub fn probe_keys(d: &[f64; 5], out: &mut Vec<u64>) {
    out.clear();
    let mut opts = [[0i64; 2]; 5];
    let mut n_opts = [1usize; 5];
    for (dim, &v) in d.iter().enumerate() {
        let b = bin(v);
        opts[dim][0] = b;
        let frac = v * BINS - b as f64;
        if frac < PROBE_EPS * BINS && b > 0 {
            opts[dim][1] = b - 1;
            n_opts[dim] = 2;
        } else if frac > 1.0 - PROBE_EPS * BINS && b < BINS as i64 - 1 {
            opts[dim][1] = b + 1;
            n_opts[dim] = 2;
        }
    }
    let total: usize = n_opts.iter().product();
    for combo in 0..total {
        let mut rest = combo;
        let mut k = 0u64;
        for dim in 0..5 {
            let pick = rest % n_opts[dim];
            rest /= n_opts[dim];
            k = (k << 8) | opts[dim][pick] as u64;
        }
        out.push(k);
    }
}

// ── Sphere ↔ tangent plane ─────────────────────────────────────────────────────

/// Unit vector of (`ra`, `dec`), radians.
#[inline]
#[must_use]
pub fn unit(ra: f64, dec: f64) -> [f64; 3] {
    let (sr, cr) = ra.sin_cos();
    let (sd, cd) = dec.sin_cos();
    [cd * cr, cd * sr, sd]
}

#[inline]
fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// A gnomonic tangent plane: its centre and local east and north unit vectors.
/// Standard coordinates are radians, east and north positive.
#[derive(Debug, Clone, Copy)]
pub struct Tangent {
    /// Unit vector of the tangent point.
    pub c: [f64; 3],
    e: [f64; 3],
    n: [f64; 3],
}

impl Tangent {
    /// The plane touching the sphere at unit vector `c` (need not be normalised).
    /// `None` at the poles' exact axis, where east is undefined.
    #[must_use]
    pub fn at(c: [f64; 3]) -> Option<Self> {
        let norm = dot(&c, &c).sqrt();
        if norm <= 0.0 {
            return None;
        }
        let c = [c[0] / norm, c[1] / norm, c[2] / norm];
        let h = c[0].hypot(c[1]);
        if h < 1e-12 {
            return None;
        }
        let e = [-c[1] / h, c[0] / h, 0.0];
        let n = [-c[2] * e[1], c[2] * e[0], h];
        Some(Self { c, e, n })
    }

    /// The plane touching (`ra`, `dec`).
    #[must_use]
    pub fn at_radec(ra: f64, dec: f64) -> Option<Self> {
        Self::at(unit(ra, dec))
    }

    /// Standard coordinates of unit vector `u`; `None` behind the plane.
    #[inline]
    #[must_use]
    pub fn project(&self, u: &[f64; 3]) -> Option<(f64, f64)> {
        let depth = dot(u, &self.c);
        if depth <= 1e-9 {
            return None;
        }
        Some((dot(u, &self.e) / depth, dot(u, &self.n) / depth))
    }

    /// (RA, Dec) radians of standard coordinates (`xi`, `eta`).
    #[must_use]
    pub fn deproject(&self, xi: f64, eta: f64) -> (f64, f64) {
        let v = [
            self.c[0] + xi * self.e[0] + eta * self.n[0],
            self.c[1] + xi * self.e[1] + eta * self.n[1],
            self.c[2] + xi * self.e[2] + eta * self.n[2],
        ];
        let norm = dot(&v, &v).sqrt();
        let ra = v[1].atan2(v[0]).rem_euclid(2.0 * core::f64::consts::PI);
        (ra, (v[2] / norm).clamp(-1.0, 1.0).asin())
    }
}

/// Least-squares affine map `q = A·p + t` from points `p` to points `q`, returned as
/// `[a, b, c, d, e, f]` with `qx = a·px + b·py + c`, `qy = d·px + e·py + f`. `None`
/// when the points are collinear or fewer than three.
#[must_use]
pub fn fit_affine(p: &[(f64, f64)], q: &[(f64, f64)]) -> Option<[f64; 6]> {
    let n = p.len().min(q.len());
    if n < 3 {
        return None;
    }
    // Centre both sets for conditioning, then solve the 2×2 normal equations.
    let (mut mx, mut my, mut nx, mut ny) = (0.0, 0.0, 0.0, 0.0);
    for i in 0..n {
        mx += p[i].0;
        my += p[i].1;
        nx += q[i].0;
        ny += q[i].1;
    }
    let k = n as f64;
    let (mx, my, nx, ny) = (mx / k, my / k, nx / k, ny / k);
    let (mut sxx, mut sxy, mut syy) = (0.0, 0.0, 0.0);
    let (mut ux, mut uy, mut vx, mut vy) = (0.0, 0.0, 0.0, 0.0);
    for i in 0..n {
        let (x, y) = (p[i].0 - mx, p[i].1 - my);
        let (u, v) = (q[i].0 - nx, q[i].1 - ny);
        sxx += x * x;
        sxy += x * y;
        syy += y * y;
        ux += u * x;
        uy += u * y;
        vx += v * x;
        vy += v * y;
    }
    let det = sxx * syy - sxy * sxy;
    if det.abs() <= 1e-12 * (sxx * syy).max(1e-300) {
        return None;
    }
    let a = (ux * syy - uy * sxy) / det;
    let b = (uy * sxx - ux * sxy) / det;
    let d = (vx * syy - vy * sxy) / det;
    let e = (vy * sxx - vx * sxy) / det;
    Some([a, b, nx - a * mx - b * my, d, e, ny - d * mx - e * my])
}

/// How far an affine map is from a similarity (rotation + uniform scale, either
/// parity): the larger of the column-length mismatch and the columns' cosine.
/// Zero for a perfect similarity.
#[must_use]
pub fn shape_error(m: &[f64; 6]) -> f64 {
    let c1 = m[0].hypot(m[3]);
    let c2 = m[1].hypot(m[4]);
    if c1 <= 0.0 || c2 <= 0.0 {
        return f64::INFINITY;
    }
    let ortho = (m[0] * m[1] + m[3] * m[4]).abs() / (c1 * c2);
    ortho.max((c1 / c2 - 1.0).abs())
}

/// Invert an affine map; `None` if singular.
#[must_use]
pub fn invert_affine(m: &[f64; 6]) -> Option<[f64; 6]> {
    let det = m[0] * m[4] - m[1] * m[3];
    if det.abs() < 1e-300 {
        return None;
    }
    let (a, b, d, e) = (m[4] / det, -m[1] / det, -m[3] / det, m[0] / det);
    Some([a, b, -(a * m[2] + b * m[5]), d, e, -(d * m[2] + e * m[5])])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quad() -> [(f64, f64); 4] {
        [(0.0, 0.0), (10.0, 1.0), (3.0, 7.0), (8.0, 9.5)]
    }

    fn transform(p: &[(f64, f64); 4], s: f64, th: f64, flip: bool) -> [(f64, f64); 4] {
        let (st, ct) = th.sin_cos();
        p.map(|(x, y)| {
            let x = if flip { -x } else { x };
            (s * (ct * x - st * y) + 100.0, s * (st * x + ct * y) - 40.0)
        })
    }

    #[test]
    fn descriptor_and_key_are_similarity_invariant() {
        let q = quad();
        let k = key(&descriptor(&q));
        for (s, th, flip) in [(2.0, 0.3, false), (0.01, 2.0, true), (37.0, -1.0, true)] {
            let t = transform(&q, s, th, flip);
            let d = descriptor(&t);
            let mut keys = Vec::new();
            probe_keys(&d, &mut keys);
            assert!(keys.contains(&k), "{s} {th} {flip}");
        }
    }

    #[test]
    fn the_canonical_order_gives_the_correspondence() {
        let q = quad();
        let t = transform(&q, 3.0, 1.1, true);
        let (o1, _, _) = canonical(&q).unwrap();
        let (o2, _, _) = canonical(&t).unwrap();
        // The same original vertex sits in each canonical slot.
        assert_eq!(o1, o2);
        let p: Vec<_> = o1.iter().map(|&i| q[i]).collect();
        let r: Vec<_> = o2.iter().map(|&i| t[i]).collect();
        let m = fit_affine(&p, &r).unwrap();
        assert!(shape_error(&m) < 1e-9);
    }

    #[test]
    fn near_ties_offer_both_orders() {
        // A square: every total ties.
        let sq = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        let (o, t, _) = canonical(&sq).unwrap();
        assert!(orderings(o, t).len() >= 4);
        let (o, t, _) = canonical(&quad()).unwrap();
        assert!(!orderings(o, t).is_empty());
    }

    #[test]
    fn probing_covers_a_value_just_across_an_edge() {
        // Bin centres, far from any edge, except the one under test.
        let mut d = [25.5 / BINS, 51.5 / BINS, 0.0, 76.5 / BINS, 115.5 / BINS];
        d[2] = 64.0 / BINS + 0.001; // just above an edge
        let mut keys = Vec::new();
        probe_keys(&d, &mut keys);
        let mut below = d;
        below[2] = 64.0 / BINS - 0.0005;
        assert!(keys.contains(&key(&below)));
        assert_eq!(keys.len(), 2);
    }

    #[test]
    fn tangent_round_trip() {
        let t = Tangent::at_radec(1.0, 0.5).unwrap();
        let (ra, dec) = (1.01, 0.49);
        let (x, y) = t.project(&unit(ra, dec)).unwrap();
        let (r2, d2) = t.deproject(x, y);
        assert!((r2 - ra).abs() < 1e-12 && (d2 - dec).abs() < 1e-12);
        // East is +xi, north is +eta.
        let (x, _) = t.project(&unit(1.001, 0.5)).unwrap();
        assert!(x > 0.0);
        let (_, y) = t.project(&unit(1.0, 0.501)).unwrap();
        assert!(y > 0.0);
    }

    #[test]
    fn affine_inverse_round_trips() {
        let m = [2.0, 0.5, 3.0, -0.4, 1.5, -7.0];
        let i = invert_affine(&m).unwrap();
        let (x, y) = (3.3, -1.2);
        let (u, v) = (m[0] * x + m[1] * y + m[2], m[3] * x + m[4] * y + m[5]);
        let (x2, y2) = (i[0] * u + i[1] * v + i[2], i[3] * u + i[4] * v + i[5]);
        assert!((x - x2).abs() < 1e-12 && (y - y2).abs() < 1e-12);
    }
}
