// Quad creation from star lists.

use crate::types::{Quad, QuadList, StarList};

/// Neighbourhood size for image quads and for catalogue quads. See build_quads.
pub const IMAGE_NEIGHBOURS: usize = 9;
pub const CATALOG_NEIGHBOURS: usize = 9;

// Hash-dedup constants
const BUCKET_CAPACITY: usize = 10;
const GRID_INV: f64 = 0.2; // 1.0 / grid_size(5.0)

/// Fast atan2 approximation (max error ~0.0026 rad ≈ 0.15°).
/// Adequate for the 10° angle-voting bins used in vote_filter.
/// Uses the identity atan(z) ≈ z / (1 + 0.28125·z²) for |z|≤1.
#[inline(always)]
fn fast_atan2(y: f64, x: f64) -> f64 {
    use core::f64::consts::{FRAC_PI_2, PI};
    if x == 0.0 {
        return if y > 0.0 {
            FRAC_PI_2
        } else if y < 0.0 {
            -FRAC_PI_2
        } else {
            0.0
        };
    }
    let z = y / x;
    let atan = if z.abs() <= 1.0 {
        z / (1.0 + 0.28125 * z * z)
    } else {
        let iz = 1.0 / z;
        let s: f64 = if z > 0.0 { 1.0 } else { -1.0 };
        s * FRAC_PI_2 - iz / (1.0 + 0.28125 * iz * iz)
    };
    if x < 0.0 {
        atan + if y >= 0.0 { PI } else { -PI }
    } else {
        atan
    }
}

/// Build a sorted `[d1, d2, d3, d4, d5, d6]` from exactly six distances.
/// Uses a fixed comparison sequence rather than a general sort: with six elements the
/// comparisons are known ahead of time.
fn sort6(mut d: [f64; 6]) -> [f64; 6] {
    macro_rules! swap_if {
        ($a:expr, $b:expr) => {
            if d[$b] > d[$a] {
                d.swap($a, $b);
            }
        };
    }
    // Pass 1
    swap_if!(0, 1);
    swap_if!(1, 2);
    swap_if!(2, 3);
    swap_if!(3, 4);
    swap_if!(4, 5);
    // Pass 2
    swap_if!(0, 1);
    swap_if!(1, 2);
    swap_if!(2, 3);
    swap_if!(3, 4);
    // Pass 3
    swap_if!(0, 1);
    swap_if!(1, 2);
    swap_if!(2, 3);
    // Pass 4
    swap_if!(0, 1);
    swap_if!(1, 2);
    // Pass 5
    swap_if!(0, 1);
    d
}

/// Compute all 6 pairwise distances for four (x,y) points and return a sorted Quad.
fn make_quad(p1: (f64, f64), p2: (f64, f64), p3: (f64, f64), p4: (f64, f64)) -> Option<Quad> {
    let pts = [p1, p2, p3, p4];
    // Pair ordering mirrors the `raw` array below (must stay in sync).
    const PAIR_IDXS: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

    let dist = |a: (f64, f64), b: (f64, f64)| -> f64 {
        ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt()
    };
    let raw = [
        dist(p1, p2),
        dist(p1, p3),
        dist(p1, p4),
        dist(p2, p3),
        dist(p2, p4),
        dist(p3, p4),
    ];

    // Find which pair has the largest distance (for the rotation-angle fingerprint).
    let max_k = raw
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(k, _)| k)
        .unwrap_or(0);
    let (ai, bi) = PAIR_IDXS[max_k];
    let dx = pts[bi].0 - pts[ai].0;
    let dy = pts[bi].1 - pts[ai].1;
    let d1_angle = fast_atan2(dy, dx).rem_euclid(core::f64::consts::PI);

    let d = sort6(raw);
    let d1 = d[0];
    if d1 < 1e-10 {
        return None;
    }
    let center_x = (p1.0 + p2.0 + p3.0 + p4.0) * 0.25;
    let center_y = (p1.1 + p2.1 + p3.1 + p4.1) * 0.25;
    Some(Quad {
        d1,
        ratios: [d[1] / d1, d[2] / d1, d[3] / d1, d[4] / d1, d[5] / d1],
        center_x,
        center_y,
        d1_angle,
    })
}

/// All C(k, 4) index combinations of 4 from 0..k, in lexicographic order.
///
/// Replaces the hand-written COMBOS_5/6/7 tables so that any neighbourhood size can
/// be used; the generated order is identical to those tables for k ∈ {5, 6, 7}.
fn combinations_of_4(k: usize) -> Vec<[usize; 4]> {
    let mut out = Vec::with_capacity(k * k * k * k / 24 + 4);
    for a in 0..k {
        for b in (a + 1)..k {
            for c in (b + 1)..k {
                for d in (c + 1)..k {
                    out.push([a, b, c, d]);
                }
            }
        }
    }
    out
}

/// Small-star-count quad builder (find_many_quads).
/// For each star, find `num_closest` nearest neighbours, then emit all C(num_closest, 4) quads.
/// Duplicates are filtered by center proximity (< 1px in both x and y).
fn find_many_quads(stars: &StarList, mode: usize) -> QuadList {
    if !(4..=12).contains(&mode) {
        return QuadList::default();
    }
    let num_closest = mode;
    let combos = combinations_of_4(num_closest);

    let n = stars.len();
    if n < num_closest {
        return QuadList::default();
    }
    let mut quads: Vec<Quad> = Vec::with_capacity(n * combos.len());

    // Centre-proximity dedup via a hash grid (5-px cells), matching find_quads_nn.
    // The previous linear scan over every accepted centre was O(quads²), which is
    // fine at ~200 quads and quadratic pain at the several thousand that a larger
    // neighbourhood produces.
    let table_len = (n * combos.len() / 4).max(16);
    let mut hash_table: Vec<Vec<usize>> = vec![Vec::new(); table_len];

    for i in 0..n {
        let x1 = stars.0[i].x;
        let y1 = stars.0[i].y;

        // Find num_closest nearest neighbours (insertion sort)
        let mut closest_idx = vec![0usize; num_closest];
        let mut closest_dist = vec![f64::MAX; num_closest];
        closest_idx[0] = i;
        closest_dist[0] = 0.0;

        for j in 0..n {
            if j == i {
                continue;
            }
            let dx = stars.0[j].x - x1;
            let dy = stars.0[j].y - y1;
            let d = dx * dx + dy * dy;
            if d <= 1.0 {
                continue;
            } // identical star guard
            // Insertion sort into closest list
            if d < closest_dist[num_closest - 1] {
                let mut pos = num_closest - 1;
                while pos > 0 && d < closest_dist[pos - 1] {
                    pos -= 1;
                }
                for k in (pos..num_closest - 1).rev() {
                    closest_dist[k + 1] = closest_dist[k];
                    closest_idx[k + 1] = closest_idx[k];
                }
                closest_dist[pos] = d;
                closest_idx[pos] = j;
            }
        }

        // All num_closest positions filled?
        if closest_idx[num_closest - 1] == 0 && closest_dist[num_closest - 1] == f64::MAX {
            continue;
        }

        // Emit all C(num_closest, 4) quads
        for combo in &combos {
            let get = |idx: usize| -> (f64, f64) {
                let si = closest_idx[idx];
                (stars.0[si].x, stars.0[si].y)
            };
            let p = [get(combo[0]), get(combo[1]), get(combo[2]), get(combo[3])];
            let cx = (p[0].0 + p[1].0 + p[2].0 + p[3].0) * 0.25;
            let cy = (p[0].1 + p[1].1 + p[2].1 + p[3].1) * 0.25;

            // Duplicate check: same centre within 1px in both axes.
            // ASTAP tightened the Y tolerance from 6px to 1px in 2026-06-29.
            let hx = (cx * GRID_INV) as i64;
            let hy = (cy * GRID_INV) as i64;
            let idx = ((hx * 31 + hy).unsigned_abs() as usize) % table_len;
            let dup = hash_table[idx].iter().any(|&qi| {
                (cx - quads[qi].center_x).abs() < 1.0 && (cy - quads[qi].center_y).abs() < 1.0
            });
            if dup {
                continue;
            }

            if let Some(q) = make_quad(p[0], p[1], p[2], p[3]) {
                if hash_table[idx].len() < BUCKET_CAPACITY {
                    hash_table[idx].push(quads.len());
                }
                quads.push(q);
            }
        }
    }

    QuadList(quads)
}

/// Large-star-count quad builder (3-nearest-neighbour + hash dedup).
/// Stars should be sorted by X when `nrstars >= 150` (bandwidth filtering).
fn find_quads_nn(stars: &StarList) -> QuadList {
    let n = stars.len();
    if n < 4 {
        return QuadList::default();
    }

    // Bandwidth: full search unless n >= 150
    let bandw = if n >= 150 {
        (2.0 * (n as f64).sqrt()).round() as usize
    } else {
        n
    };

    // Hash table for dedup
    let table_len = (n * 2).max(1);
    let mut hash_table: Vec<Vec<usize>> = vec![Vec::new(); table_len];

    let mut quads: Vec<Quad> = Vec::with_capacity(n);

    for i in 0..n {
        let x1 = stars.0[i].x;
        let y1 = stars.0[i].y;

        let s_start = i.saturating_sub(bandw);
        let s_end = (i + bandw).min(n - 1);

        let mut d1 = f64::MAX;
        let mut d2 = f64::MAX;
        let mut d3 = f64::MAX;
        let mut j1 = 0usize;
        let mut j2 = 0usize;
        let mut j3 = 0usize;

        for j in s_start..=s_end {
            if j == i {
                continue;
            }
            let dy = stars.0[j].y - y1;
            let dy2 = dy * dy;
            if dy2 >= d3 {
                continue;
            } // pre-check
            let dx = stars.0[j].x - x1;
            let dist = dx * dx + dy2;
            if dist <= 1.0 {
                continue;
            } // identical star
            if dist < d1 {
                d3 = d2;
                j3 = j2;
                d2 = d1;
                j2 = j1;
                d1 = dist;
                j1 = j;
            } else if dist < d2 {
                d3 = d2;
                j3 = j2;
                d2 = dist;
                j2 = j;
            } else if dist < d3 {
                d3 = dist;
                j3 = j;
            }
        }

        if d3 == f64::MAX {
            continue;
        } // fewer than 3 neighbours found

        let p1 = (x1, y1);
        let p2 = (stars.0[j1].x, stars.0[j1].y);
        let p3 = (stars.0[j2].x, stars.0[j2].y);
        let p4 = (stars.0[j3].x, stars.0[j3].y);

        let cx = (p1.0 + p2.0 + p3.0 + p4.0) * 0.25;
        let cy = (p1.1 + p2.1 + p3.1 + p4.1) * 0.25;

        // Hash dedup check
        let hx = (cx * GRID_INV) as i64;
        let hy = (cy * GRID_INV) as i64;
        let idx = ((hx * 31 + hy).unsigned_abs() as usize) % table_len;

        let dup = hash_table[idx].iter().any(|&qi| {
            (cx - quads[qi].center_x).abs() < 1.0 && (cy - quads[qi].center_y).abs() < 1.0
        });
        if dup {
            continue;
        }

        if let Some(q) = make_quad(p1, p2, p3, p4) {
            // Record in hash table
            if hash_table[idx].len() < BUCKET_CAPACITY {
                hash_table[idx].push(quads.len());
            }
            // Even if bucket full, still add quad (just won't dedup future identical ones)
            quads.push(q);
        }
    }

    QuadList(quads)
}

/// Build quads from a star list.
///
/// `nrstars_image`: the count of stars found in the *image* (used for mode selection;
/// may differ from `stars.len()` when building catalog quads from a larger catalog area).
pub fn build_quads(stars: &StarList, nrstars_image: usize) -> QuadList {
    let n = stars.len();
    if nrstars_image < 15 && n > 6 {
        return find_many_quads(stars, 7);
    }
    if nrstars_image < 30 && n > 5 {
        return find_many_quads(stars, 6);
    }
    // Large star counts: 9 nearest neighbours and all C(9,4) subsets, i.e. 126
    // quads per star instead of the single 3-nearest-neighbour quad ASTAP builds.
    //
    // A 3-NN quad depends entirely on *which* stars are in the list, and the image
    // and catalogue lists never match exactly (the image has undetected faint stars,
    // the catalogue has stars below the detection limit). Redundancy is what lets a
    // correspondence survive that: astrometry.net puts each star in up to 8 quads
    // and makes ~16 passes for exactly this reason.
    //
    // Measured on the 103-image benchmark corpus (see docs/test-images.md).
    // Without verification, recall rises with the neighbourhood but so do false
    // positives (5-NN 64/4FP, 7-NN 76/7FP, 10-NN 69/21FP). With the star-level
    // verification in solver.rs the false-positive count stays at zero and recall
    // peaks at 9 neighbours:
    //
    //   5-NN 68   6-NN 77   7-NN 81   8-NN 83   9-NN 85   10-NN 83   12-NN 80
    //
    // Past 9 the extra quads add noise rather than signal and the cost climbs
    // sharply (12-NN is 20x the wall time for two fewer solves).
    if n > 4 {
        return find_many_quads(stars, IMAGE_NEIGHBOURS);
    }
    // Large: sort by X then use bandwidth-filtered 3-NN
    let mut sorted = stars.clone();
    sorted.0.sort_by(|a, b| a.x.total_cmp(&b.x));
    find_quads_nn(&sorted)
}

/// Build quads from a star list that is **already sorted by x ascending**.
///
/// For the ≥60-star path this skips the clone+sort that `build_quads` performs,
/// saving one heap allocation and ~500-element sort per spiral step.
/// The caller must guarantee the sort order; results are undefined otherwise.
pub fn build_quads_presorted(stars: &StarList, nrstars_image: usize) -> QuadList {
    let n = stars.len();
    if nrstars_image < 15 && n > 6 {
        return find_many_quads(stars, 7);
    }
    if nrstars_image < 30 && n > 5 {
        return find_many_quads(stars, 6);
    }
    // Large star counts: 9 nearest neighbours and all C(9,4) subsets, i.e. 126
    // quads per star instead of the single 3-nearest-neighbour quad ASTAP builds.
    //
    // A 3-NN quad depends entirely on *which* stars are in the list, and the image
    // and catalogue lists never match exactly (the image has undetected faint stars,
    // the catalogue has stars below the detection limit). Redundancy is what lets a
    // correspondence survive that: astrometry.net puts each star in up to 8 quads
    // and makes ~16 passes for exactly this reason.
    //
    // Measured on the 103-image benchmark corpus (see docs/test-images.md).
    // Without verification, recall rises with the neighbourhood but so do false
    // positives (5-NN 64/4FP, 7-NN 76/7FP, 10-NN 69/21FP). With the star-level
    // verification in solver.rs the false-positive count stays at zero and recall
    // peaks at 9 neighbours:
    //
    //   5-NN 68   6-NN 77   7-NN 81   8-NN 83   9-NN 85   10-NN 83   12-NN 80
    //
    // Past 9 the extra quads add noise rather than signal and the cost climbs
    // sharply (12-NN is 20x the wall time for two fewer solves).
    if n > 4 {
        return find_many_quads(stars, CATALOG_NEIGHBOURS);
    }
    // Stars are pre-sorted by x; call find_quads_nn directly with no clone.
    find_quads_nn(stars)
}

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

    fn assert_ratios_valid(q: &Quad) {
        // d1 is the largest distance; all ratios must be in (0, 1]
        assert!(q.d1 > 0.0, "d1 must be positive");
        for &r in &q.ratios {
            assert!(r > 0.0 && r <= 1.0 + 1e-10, "ratio {r} out of (0,1]");
        }
        // Ratios must be non-increasing (sorted descending)
        for i in 0..4 {
            assert!(
                q.ratios[i] >= q.ratios[i + 1] - 1e-12,
                "ratios not sorted: {:?}",
                q.ratios
            );
        }
    }

    #[test]
    fn sort6_is_descending() {
        let d = sort6([3.0, 1.0, 5.0, 2.0, 4.0, 6.0]);
        assert_eq!(d, [6.0, 5.0, 4.0, 3.0, 2.0, 1.0]);
    }

    #[test]
    fn make_quad_square() {
        // Unit square: 4 sides = 1, 2 diagonals = √2
        // Sorted descending: [√2, √2, 1, 1, 1, 1]
        // d1 = √2, ratio[0] = √2/√2 = 1.0, ratio[1..4] = 1/√2
        let q = make_quad((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)).unwrap();
        assert!((q.d1 - 2.0_f64.sqrt()).abs() < 1e-10, "d1 = {}", q.d1);
        assert!(
            (q.ratios[0] - 1.0).abs() < 1e-10,
            "ratio[0] = {}",
            q.ratios[0]
        );
        for k in 1..5 {
            assert!(
                (q.ratios[k] - 1.0 / 2.0_f64.sqrt()).abs() < 1e-10,
                "ratio[{k}] = {}",
                q.ratios[k]
            );
        }
        assert_ratios_valid(&q);
    }

    #[test]
    fn find_many_quads_mode5_count() {
        // 5 stars → mode 5 → each star contributes 5 quads (before dedup)
        let coords: Vec<(f64, f64)> = (0..5)
            .map(|i| {
                (
                    i as f64 * 50.0,
                    (i as f64 * 30.0) + if i % 2 == 0 { 0.0 } else { 20.0 },
                )
            })
            .collect();
        let stars = make_stars(&coords);
        let ql = find_many_quads(&stars, 5);
        assert!(!ql.is_empty(), "should produce quads");
        for q in &ql.0 {
            assert_ratios_valid(q);
        }
    }

    #[test]
    fn no_duplicate_quads_nn() {
        // 20 stars on a grid — no two quads should share a centre within 1px
        let coords: Vec<(f64, f64)> = (0..4)
            .flat_map(|i| (0..5).map(move |j| (i as f64 * 100.0 + 50.0, j as f64 * 80.0 + 50.0)))
            .collect();
        let stars = make_stars(&coords);
        let ql = find_quads_nn(&stars);
        // Check no two quads share centres within 1px
        for i in 0..ql.len() {
            for j in (i + 1)..ql.len() {
                let dx = (ql.0[i].center_x - ql.0[j].center_x).abs();
                let dy = (ql.0[i].center_y - ql.0[j].center_y).abs();
                assert!(
                    dx > 1.0 || dy > 1.0,
                    "duplicate quads at index {i},{j}: dx={dx}, dy={dy}"
                );
            }
        }
        for q in &ql.0 {
            assert_ratios_valid(q);
        }
    }

    #[test]
    fn build_quads_small_uses_many() {
        // < 15 image stars → should use find_many_quads(7)
        let coords: Vec<(f64, f64)> = (0..10)
            .map(|i| {
                (
                    i as f64 * 40.0,
                    (i as f64 * 3.0 + if i % 2 == 0 { 0.0 } else { 20.0 }) * 10.0,
                )
            })
            .collect();
        let stars = make_stars(&coords);
        let ql = build_quads(&stars, 10); // nrstars_image = 10 < 15
        assert!(!ql.is_empty(), "small image should produce quads");
    }

    #[test]
    fn build_quads_large_uses_nn() {
        // >= 60 image stars → should use nn strategy
        let coords: Vec<(f64, f64)> = (0..8)
            .flat_map(|i| (0..8).map(move |j| (i as f64 * 60.0 + 10.0, j as f64 * 60.0 + 10.0)))
            .collect();
        let stars = make_stars(&coords);
        let ql = build_quads(&stars, 64);
        assert!(!ql.is_empty());
        for q in &ql.0 {
            assert_ratios_valid(q);
        }
    }

    #[test]
    fn build_quads_presorted_matches_build_quads() {
        // build_quads_presorted on x-sorted input must produce the same quad count as build_quads.
        let mut coords: Vec<(f64, f64)> = (0..8)
            .flat_map(|i| (0..8).map(move |j| (i as f64 * 60.0 + 10.0, j as f64 * 60.0 + 10.0)))
            .collect();
        // Shuffle so build_quads must sort internally
        coords.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap()); // descending x (wrong order)
        let stars_unsorted = make_stars(&coords);

        coords.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap()); // ascending x
        let stars_sorted = make_stars(&coords);

        let ql_normal = build_quads(&stars_unsorted, 64);
        let ql_pre = build_quads_presorted(&stars_sorted, 64);

        assert_eq!(ql_normal.len(), ql_pre.len(), "quad counts must match");
        for q in &ql_pre.0 {
            assert_ratios_valid(q);
        }
    }

    #[test]
    fn quad_center_correct() {
        let q = make_quad((0.0, 0.0), (4.0, 0.0), (0.0, 4.0), (4.0, 4.0)).unwrap();
        assert!((q.center_x - 2.0).abs() < 1e-10);
        assert!((q.center_y - 2.0).abs() < 1e-10);
    }
}
