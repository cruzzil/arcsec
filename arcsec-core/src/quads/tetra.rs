// TETRA triangle-based star pattern matching.
//
// Uses 3-star triangles characterised by two normalized side-length ratios
// (s2/s1, s3/s1 where s1 ≥ s2 ≥ s3).  Fewer patterns than 4-star quads and
// a 2-D instead of 5-D feature space, so faster to build and less
// discriminating.
//
// FALSE-POSITIVE RATE: with N features and per-feature tolerance t, the brute-
// force false-positive rate scales as (2t)^N.  For quads (N=5) the default
// tolerance 0.007 gives (0.014)^5 ≈ 5e-11 per pair — excellent.  For
// triangles (N=2) the same tolerance gives (0.014)^2 ≈ 2e-4 — 4 million times
// worse.  Use TETRA_TOL_FACTOR to scale the tolerance down accordingly.
//
// Reference: Mortari et al. "TETRA: A Celestial Body Detection Algorithm",
// Journal of Guidance Control and Dynamics, 2004.

/// Multiply quad_tolerance by this factor when calling find_triangle_matches
/// to achieve a comparable false-positive rate to the 5-ratio quad matcher.
/// (0.007^(1/2.5) ≈ 0.007^0.4 ≈ 0.3  — equal false-positive rate heuristic)
pub const TETRA_TOL_FACTOR: f64 = 0.3;

use crate::types::{PairedPositions, StarList};

// ── Public types ──────────────────────────────────────────────────────────────

/// A 3-star triangle pattern with two normalized side-length ratios.
#[derive(Debug, Clone)]
pub struct Triangle {
    /// Normalized side ratios: [s2/s1, s3/s1] where s1 ≥ s2 ≥ s3.
    pub ratios: [f64; 2],
    pub center_x: f64,
    pub center_y: f64,
    /// Largest side (absolute, pixels or standard coords) — used for scale.
    pub d_max: f64,
}

/// List of triangles.
#[derive(Debug, Clone, Default)]
pub struct TriangleList(pub Vec<Triangle>);

impl TriangleList {
    pub fn len(&self) -> usize {
        self.0.len()
    }
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

/// A pair of matching triangles: image index, catalog index, scale ratio.
#[derive(Debug, Clone, Copy)]
pub struct TriMatch {
    pub img_idx: usize,
    pub cat_idx: usize,
    /// d_max_img / d_max_cat — pixel/catalog scale ratio for this match.
    pub scale_ratio: f64,
}

// ── Builder ───────────────────────────────────────────────────────────────────

/// Build a triangle from three (x,y) points, or `None` if degenerate.
fn make_triangle(p1: (f64, f64), p2: (f64, f64), p3: (f64, f64)) -> Option<Triangle> {
    let dist = |a: (f64, f64), b: (f64, f64)| ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt();

    let mut sides = [dist(p1, p2), dist(p1, p3), dist(p2, p3)];
    // Sort descending
    if sides[1] > sides[0] {
        sides.swap(0, 1);
    }
    if sides[2] > sides[0] {
        sides.swap(0, 2);
    }
    if sides[2] > sides[1] {
        sides.swap(1, 2);
    }

    let s1 = sides[0];
    if s1 < 1e-10 {
        return None;
    }

    Some(Triangle {
        ratios: [sides[1] / s1, sides[2] / s1],
        center_x: (p1.0 + p2.0 + p3.0) / 3.0,
        center_y: (p1.1 + p2.1 + p3.1) / 3.0,
        d_max: s1,
    })
}

/// Build triangles from a star list.
///
/// For ≤ 45 stars: all C(N,3) combinations.
/// For > 45 stars: each star paired with every pair of its K=7 nearest neighbours
/// → C(7,2)=21 triangles each (the star itself is always the third vertex).
pub fn build_triangles(stars: &StarList) -> TriangleList {
    let n = stars.len();
    if n < 3 {
        return TriangleList::default();
    }

    if n <= 45 {
        build_triangles_all_combos(stars)
    } else {
        build_triangles_nn(stars)
    }
}

fn build_triangles_all_combos(stars: &StarList) -> TriangleList {
    let n = stars.len();
    let mut tris = Vec::with_capacity(n * (n - 1) * (n - 2) / 6);
    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                let p1 = (stars.0[i].x, stars.0[i].y);
                let p2 = (stars.0[j].x, stars.0[j].y);
                let p3 = (stars.0[k].x, stars.0[k].y);
                if let Some(t) = make_triangle(p1, p2, p3) {
                    tris.push(t);
                }
            }
        }
    }
    TriangleList(tris)
}

fn build_triangles_nn(stars: &StarList) -> TriangleList {
    const K: usize = 7;
    /// Triangles emitted per star: every pair drawn from its K neighbours, with
    /// the star itself as the third vertex. C(7,2) = 21, not C(7,3) = 35.
    const TRIS_PER_STAR: usize = K * (K - 1) / 2;
    let n = stars.len();
    let mut tris = Vec::with_capacity(n * TRIS_PER_STAR);

    for i in 0..n {
        let xi = stars.0[i].x;
        let yi = stars.0[i].y;

        // Find K nearest neighbours (insertion sort)
        let mut nn_idx = [0usize; K];
        let mut nn_d2 = [f64::MAX; K];

        for j in 0..n {
            if j == i {
                continue;
            }
            let dx = stars.0[j].x - xi;
            let dy = stars.0[j].y - yi;
            let d2 = dx * dx + dy * dy;
            if d2 <= 1.0 {
                continue;
            }
            if d2 < nn_d2[K - 1] {
                let mut pos = K - 1;
                while pos > 0 && d2 < nn_d2[pos - 1] {
                    pos -= 1;
                }
                for m in (pos..K - 1).rev() {
                    nn_d2[m + 1] = nn_d2[m];
                    nn_idx[m + 1] = nn_idx[m];
                }
                nn_d2[pos] = d2;
                nn_idx[pos] = j;
            }
        }

        // Count valid neighbours
        let valid = nn_d2.iter().take_while(|&&d| d < f64::MAX).count();
        if valid < 2 {
            continue;
        }

        // All C(valid, 2) pairs with star i
        for a in 0..valid {
            for b in (a + 1)..valid {
                let p1 = (xi, yi);
                let p2 = (stars.0[nn_idx[a]].x, stars.0[nn_idx[a]].y);
                let p3 = (stars.0[nn_idx[b]].x, stars.0[nn_idx[b]].y);
                if let Some(t) = make_triangle(p1, p2, p3) {
                    tris.push(t);
                }
            }
        }
    }
    TriangleList(tris)
}

// ── Matching ──────────────────────────────────────────────────────────────────

/// Find all triangle pairs whose 2 ratios agree within `tolerance`.
pub fn find_triangle_matches(
    img: &TriangleList,
    cat: &TriangleList,
    tolerance: f64,
) -> Vec<TriMatch> {
    let mut matches = Vec::new();
    for (i, it) in img.0.iter().enumerate() {
        for (j, ct) in cat.0.iter().enumerate() {
            if ct.d_max < 1e-10 {
                continue;
            }
            if (it.ratios[0] - ct.ratios[0]).abs() <= tolerance
                && (it.ratios[1] - ct.ratios[1]).abs() <= tolerance
            {
                matches.push(TriMatch {
                    img_idx: i,
                    cat_idx: j,
                    scale_ratio: it.d_max / ct.d_max,
                });
            }
        }
    }
    matches
}

/// Bijective (1-to-1) filter: for each image triangle keep only its closest catalog
/// triangle (by L1 ratio distance), and for each catalog triangle keep only its closest
/// image triangle.  Eliminates the many-to-one false matches that plague brute-force
/// triangle comparison when 2 ratio features produce multiple near-equal hits.
pub fn bijective_filter(
    matches: &[TriMatch],
    img: &TriangleList,
    cat: &TriangleList,
) -> Vec<TriMatch> {
    use std::collections::HashMap;

    // Step 1: best catalog match for each image triangle
    let mut best_img: HashMap<usize, (usize, f64)> = HashMap::new();
    for m in matches {
        let it = &img.0[m.img_idx];
        let ct = &cat.0[m.cat_idx];
        let d = (it.ratios[0] - ct.ratios[0]).abs() + (it.ratios[1] - ct.ratios[1]).abs();
        let e = best_img.entry(m.img_idx).or_insert((m.cat_idx, f64::MAX));
        if d < e.1 {
            *e = (m.cat_idx, d);
        }
    }

    // Step 2: best image triangle for each catalog triangle (among step-1 winners)
    let mut best_cat: HashMap<usize, (usize, f64)> = HashMap::new();
    for (&img_idx, &(cat_idx, d)) in &best_img {
        let e = best_cat.entry(cat_idx).or_insert((img_idx, f64::MAX));
        if d < e.1 {
            *e = (img_idx, d);
        }
    }

    // Step 3: collect mutual best pairs
    best_cat
        .into_iter()
        .map(|(cat_idx, (img_idx, _))| TriMatch {
            img_idx,
            cat_idx,
            scale_ratio: img.0[img_idx].d_max / cat.0[cat_idx].d_max,
        })
        .collect()
}

/// Remove scale outliers from triangle matches (same logic as quad filter_by_scale).
pub fn filter_triangles_by_scale(
    matches: &[TriMatch],
    quad_tolerance: f64,
) -> (Vec<TriMatch>, f64) {
    if matches.is_empty() {
        return (vec![], 0.0);
    }
    let mut ratios: Vec<f64> = matches.iter().map(|m| m.scale_ratio).collect();
    ratios.sort_by(|a, b| a.total_cmp(b));
    let mid = ratios.len() / 2;
    let med = if ratios.len().is_multiple_of(2) {
        (ratios[mid - 1] + ratios[mid]) * 0.5
    } else {
        ratios[mid]
    };
    if med < 1e-12 {
        return (vec![], 0.0);
    }
    let tol = quad_tolerance * med;
    let filtered = matches
        .iter()
        .filter(|m| (m.scale_ratio - med).abs() <= tol)
        .copied()
        .collect();
    (filtered, med)
}

/// Extract (image_xy, catalog_xy) centre pairs from matched triangles.
pub fn extract_triangle_pairs(
    img: &TriangleList,
    cat: &TriangleList,
    matches: &[TriMatch],
) -> PairedPositions {
    let mut img_pos = Vec::with_capacity(matches.len());
    let mut cat_pos = Vec::with_capacity(matches.len());
    for m in matches {
        let it = &img.0[m.img_idx];
        let ct = &cat.0[m.cat_idx];
        img_pos.push((it.center_x, it.center_y));
        cat_pos.push((ct.center_x, ct.center_y));
    }
    (img_pos, cat_pos)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Star;

    fn make_stars(coords: &[(f64, f64)]) -> StarList {
        StarList(
            coords
                .iter()
                .map(|&(x, y)| Star {
                    x,
                    y,
                    snr: 100.0,
                    hfd: 2.0,
                })
                .collect(),
        )
    }

    #[test]
    fn make_triangle_equilateral() {
        // Equilateral triangle: all sides equal → both ratios = 1.0
        let t = make_triangle((0.0, 0.0), (2.0, 0.0), (1.0, 1.732)).unwrap();
        assert!((t.ratios[0] - 1.0).abs() < 0.01, "r0={}", t.ratios[0]);
        assert!((t.ratios[1] - 1.0).abs() < 0.01, "r1={}", t.ratios[1]);
        assert!(t.d_max > 1.9 && t.d_max < 2.1);
    }

    #[test]
    fn make_triangle_right_angle() {
        // 3-4-5 right triangle: sides 5, 4, 3 → ratios [4/5, 3/5]
        let t = make_triangle((0.0, 0.0), (4.0, 0.0), (0.0, 3.0)).unwrap();
        assert!((t.d_max - 5.0).abs() < 1e-10);
        assert!((t.ratios[0] - 0.8).abs() < 1e-10, "r0={}", t.ratios[0]);
        assert!((t.ratios[1] - 0.6).abs() < 1e-10, "r1={}", t.ratios[1]);
    }

    #[test]
    fn build_triangles_all_combos_count() {
        // 6 stars → C(6,3) = 20 triangles
        let coords: Vec<(f64, f64)> = (0..6).map(|i| (i as f64 * 50.0, i as f64 * 30.0)).collect();
        let stars = make_stars(&coords);
        let tl = build_triangles(&stars);
        assert_eq!(tl.len(), 20, "expected 20 triangles, got {}", tl.len());
    }

    #[test]
    fn identical_triangles_match() {
        let coords = [(0.0, 0.0), (100.0, 0.0), (50.0, 86.6), (200.0, 0.0)];
        let img_stars = make_stars(&coords[..3]);
        let cat_stars = make_stars(&coords[..3]);
        let img_tris = build_triangles(&img_stars);
        let cat_tris = build_triangles(&cat_stars);
        let m = find_triangle_matches(&img_tris, &cat_tris, 0.01);
        assert!(!m.is_empty(), "should find at least one match");
        for mm in &m {
            assert!(
                (mm.scale_ratio - 1.0).abs() < 1e-6,
                "scale={}",
                mm.scale_ratio
            );
        }
    }

    #[test]
    fn scaled_triangles_match() {
        // Same geometry, catalog is 2× smaller
        let img_coords = vec![(0.0, 0.0), (200.0, 0.0), (100.0, 173.2)];
        let cat_coords = vec![(0.0, 0.0), (100.0, 0.0), (50.0, 86.6)];
        let img_tris = build_triangles(&make_stars(&img_coords));
        let cat_tris = build_triangles(&make_stars(&cat_coords));
        let m = find_triangle_matches(&img_tris, &cat_tris, 0.01);
        assert!(!m.is_empty());
        for mm in &m {
            assert!(
                (mm.scale_ratio - 2.0).abs() < 1e-6,
                "scale={}",
                mm.scale_ratio
            );
        }
    }

    #[test]
    fn filter_removes_outlier() {
        let matches = vec![
            TriMatch {
                img_idx: 0,
                cat_idx: 0,
                scale_ratio: 2.0,
            },
            TriMatch {
                img_idx: 1,
                cat_idx: 1,
                scale_ratio: 2.02,
            },
            TriMatch {
                img_idx: 2,
                cat_idx: 2,
                scale_ratio: 5.0,
            }, // outlier
        ];
        let (filt, med) = filter_triangles_by_scale(&matches, 0.1);
        assert!(med > 1.9 && med < 2.1, "med={med}");
        assert_eq!(filt.len(), 2);
    }
}
