//! Quad creation from star lists.

use crate::types::{Quad, QuadList, StarList};

/// Neighbourhood size for image quads. See `build_quads`.
pub const IMAGE_NEIGHBOURS: usize = 9;
/// Neighbourhood size for catalogue quads. See `build_quads_presorted`.
pub const CATALOG_NEIGHBOURS: usize = 9;

// Hash-dedup constants
const BUCKET_CAPACITY: usize = 10;
const GRID_INV: f64 = 0.2; // 1.0 / grid_size(5.0)

/// Fast atan2 approximation (max error ~0.0026 rad ≈ 0.15°).
/// Adequate for the 10° angle-voting bins used in `vote_filter`.
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

/// `a.rem_euclid(PI)`, bit for bit, without its `fmod` call for |a| < π, where
/// `a % PI` is `a` exactly: `fast_atan2` returns -π..=π, so that is nearly always.
#[inline(always)]
fn mod_pi(a: f64) -> f64 {
    use core::f64::consts::PI;
    let r = if a.abs() < PI { a } else { a % PI };
    if r < 0.0 { r + PI } else { r }
}

/// Build a sorted `[d1, d2, d3, d4, d5, d6]` (largest first) from exactly six
/// distances.
///
/// A sorting network of `max`/`min` pairs, which compile to branch-free
/// instructions: the comparisons in a quad's distances go either way at random, so
/// a compare-and-swap that branches mispredicts half the time, at every one of the
/// tens of thousands of quads built per spiral position. Any correct sort gives the
/// same array.
fn sort6(d: [f64; 6]) -> [f64; 6] {
    let mut d = d;
    // The optimal 12-comparator network for six inputs.
    for (a, b) in [
        (0, 5),
        (1, 3),
        (2, 4),
        (1, 2),
        (3, 4),
        (0, 3),
        (2, 5),
        (0, 1),
        (2, 3),
        (4, 5),
        (1, 2),
        (3, 4),
    ] {
        let (hi, lo) = (d[a].max(d[b]), d[a].min(d[b]));
        d[a] = hi;
        d[b] = lo;
    }
    d
}

/// Compute all 6 pairwise distances for four (x,y) points and return a sorted Quad.
fn make_quad(p1: (f64, f64), p2: (f64, f64), p3: (f64, f64), p4: (f64, f64)) -> Option<Quad> {
    // Pair ordering mirrors the `raw` array below (must stay in sync).
    const PAIR_IDXS: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
    let pts = [p1, p2, p3, p4];

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
        .map_or(0, |(k, _)| k);
    let (ai, bi) = PAIR_IDXS[max_k];
    let dx = pts[bi].0 - pts[ai].0;
    let dy = pts[bi].1 - pts[ai].1;
    let d1_angle = mod_pi(fast_atan2(dy, dx));

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

/// The centres of the quads accepted so far, hashed into buckets of at most
/// [`BUCKET_CAPACITY`] for the duplicate check.
///
/// Each bucket's centres sit together in one array rather than in a `Vec` per
/// bucket, looked up through the quad list: the catalogue quads are rebuilt at
/// every spiral position, which meant thousands of small allocations each time
/// and a cache miss per comparison. Which centres a bucket holds, and so which
/// quads are kept, is unchanged.
struct CentreTable {
    /// Entries in each bucket.
    len: Vec<u8>,
    /// Each bucket's centres.
    xy: Vec<[(f64, f64); BUCKET_CAPACITY]>,
}

impl CentreTable {
    fn new(buckets: usize) -> Self {
        Self {
            len: vec![0; buckets],
            xy: vec![[(0.0, 0.0); BUCKET_CAPACITY]; buckets],
        }
    }

    /// Whether bucket `b` holds a centre within 1 unit of `(cx, cy)` in both axes.
    fn near(&self, b: usize, cx: f64, cy: f64) -> bool {
        self.xy[b][..usize::from(self.len[b])]
            .iter()
            .any(|&(x, y)| (cx - x).abs() < 1.0 && (cy - y).abs() < 1.0)
    }

    /// Record a centre in bucket `b`, unless the bucket is full.
    fn insert(&mut self, b: usize, cx: f64, cy: f64) {
        let n = usize::from(self.len[b]);
        if n < BUCKET_CAPACITY {
            self.xy[b][n] = (cx, cy);
            self.len[b] += 1;
        }
    }
}

/// All C(k, 4) index combinations of 4 from 0..k, in lexicographic order.
///
/// Replaces the hand-written `COMBOS_5/6/7` tables so that any neighbourhood size can
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

/// A star list in x order, for nearest-neighbour searches.
struct XOrder {
    /// Star indices in increasing x.
    by_x: Vec<usize>,
    /// Each star's place in `by_x`.
    rank: Vec<usize>,
}

impl XOrder {
    fn new(stars: &StarList) -> Self {
        let mut by_x: Vec<usize> = (0..stars.len()).collect();
        by_x.sort_by(|&a, &b| stars.0[a].x.total_cmp(&stars.0[b].x));
        let mut rank = vec![0; by_x.len()];
        for (r, &i) in by_x.iter().enumerate() {
            rank[i] = r;
        }
        Self { by_x, rank }
    }

    /// Fill `closest_dist` / `closest_idx` (squared distance and index, nearest
    /// first) with star `i` itself and its `len - 1` nearest neighbours, skipping
    /// stars within 1 unit (the same star twice). Slots left unfilled keep
    /// `f64::MAX` and index 0.
    ///
    /// The list holds the nearest by `(distance, index)`, which is exactly what a
    /// scan of every star in index order keeps, but only the stars near `i` in x
    /// are looked at: walking out from `i` in x order, once the list is full a star
    /// further away in x than the furthest kept cannot get in. Catalogue quads are
    /// built at every spiral position, and the full scan was a third of it.
    fn nearest(
        &self,
        stars: &StarList,
        i: usize,
        closest_dist: &mut [f64],
        closest_idx: &mut [usize],
    ) {
        closest_dist.fill(f64::MAX);
        closest_idx.fill(0);
        closest_dist[0] = 0.0;
        closest_idx[0] = i;
        let (x1, y1) = (stars.0[i].x, stars.0[i].y);
        let last = closest_dist.len() - 1;
        let n = self.by_x.len();
        let (mut up, mut down) = (self.rank[i] + 1, self.rank[i]);
        loop {
            let reach = closest_dist[last];
            // A NaN x ends the walk, which loses nothing: NaN sorts to the ends of
            // `by_x`, and a NaN distance never enters the list.
            let up_ok = up < n && (stars.0[self.by_x[up]].x - x1).powi(2) <= reach;
            let down_ok = down > 0 && (x1 - stars.0[self.by_x[down - 1]].x).powi(2) <= reach;
            if !up_ok && !down_ok {
                break;
            }
            for (go, j) in [(up_ok, up), (down_ok, down.wrapping_sub(1))] {
                if go {
                    let sj = &stars.0[self.by_x[j]];
                    let (dx, dy) = (sj.x - x1, sj.y - y1);
                    let d = dx * dx + dy * dy;
                    if d > 1.0 {
                        insert_neighbour(closest_dist, closest_idx, d, self.by_x[j]);
                    }
                }
            }
            up += usize::from(up_ok);
            down -= usize::from(down_ok);
        }
    }
}

/// Insert neighbour `j` at squared distance `d` into a nearest-first list, which
/// keeps the smallest by `(distance, index)`; slot 0 holds the star itself.
fn insert_neighbour(closest_dist: &mut [f64], closest_idx: &mut [usize], d: f64, j: usize) {
    let last = closest_dist.len() - 1;
    let before =
        |pos: usize| d < closest_dist[pos] || (d == closest_dist[pos] && j < closest_idx[pos]);
    if !before(last) {
        return;
    }
    let mut pos = last;
    while pos > 0 && before(pos - 1) {
        pos -= 1;
    }
    for k in (pos..last).rev() {
        closest_dist[k + 1] = closest_dist[k];
        closest_idx[k + 1] = closest_idx[k];
    }
    closest_dist[pos] = d;
    closest_idx[pos] = j;
}

/// Small-star-count quad builder (`find_many_quads`).
/// For each star, find `num_closest` nearest neighbours, then emit all `C(num_closest, 4)` quads.
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
    let mut centres = CentreTable::new(table_len);

    let x_order = XOrder::new(stars);

    for i in 0..n {
        // Find num_closest nearest neighbours.
        let mut closest_idx = vec![0usize; num_closest];
        let mut closest_dist = vec![f64::MAX; num_closest];
        x_order.nearest(stars, i, &mut closest_dist, &mut closest_idx);

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
            if centres.near(idx, cx, cy) {
                continue;
            }

            if let Some(q) = make_quad(p[0], p[1], p[2], p[3]) {
                // The quad's own centre, which is (cx, cy) up to rounding.
                centres.insert(idx, q.center_x, q.center_y);
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

        for (j, sj) in (s_start..=s_end).zip(&stars.0[s_start..=s_end]) {
            if j == i {
                continue;
            }
            let dy = sj.y - y1;
            let dy2 = dy * dy;
            if dy2 >= d3 {
                continue;
            } // pre-check
            let dx = sj.x - x1;
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
#[must_use]
pub fn build_quads(stars: &StarList, nrstars_image: usize) -> QuadList {
    if let Some(quads) = build_quads_neighbourhood(stars, nrstars_image, IMAGE_NEIGHBOURS) {
        return quads;
    }
    // Four stars or fewer: bandwidth-filtered 3-NN, which wants the list sorted by x.
    let mut sorted = stars.clone();
    sorted.0.sort_by(|a, b| a.x.total_cmp(&b.x));
    find_quads_nn(&sorted)
}

/// Build quads from a star list that is **already sorted by x ascending**.
///
/// Identical to [`build_quads`] (with [`CATALOG_NEIGHBOURS`]) except that the
/// small-list 3-NN path skips its clone and sort. The caller must guarantee the sort
/// order; results are unspecified otherwise.
#[must_use]
pub fn build_quads_presorted(stars: &StarList, nrstars_image: usize) -> QuadList {
    build_quads_neighbourhood(stars, nrstars_image, CATALOG_NEIGHBOURS)
        .unwrap_or_else(|| find_quads_nn(stars))
}

/// The all-subsets-of-k-neighbours builders, or `None` when the list is too short for
/// any of them and the caller should fall back to 3-NN.
fn build_quads_neighbourhood(
    stars: &StarList,
    nrstars_image: usize,
    neighbours: usize,
) -> Option<QuadList> {
    let n = stars.len();
    if nrstars_image < 15 && n > 6 {
        return Some(find_many_quads(stars, 7));
    }
    if nrstars_image < 30 && n > 5 {
        return Some(find_many_quads(stars, 6));
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
    //
    // Note that `find_many_quads` needs at least `neighbours` stars, so a list of
    // 5..neighbours stars that reaches this point yields no quads at all.
    if n > 4 {
        return Some(find_many_quads(stars, neighbours));
    }
    None
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
    fn sort6_sorts_every_arrangement_with_and_without_ties() {
        // Every 0/1 input sorted means the network sorts everything (the 0-1
        // principle); random inputs with ties check it against the library sort.
        for bits in 0u32..64 {
            let d: [f64; 6] = core::array::from_fn(|k| f64::from((bits >> k) & 1));
            let mut want = d;
            want.sort_by(|a, b| b.total_cmp(a));
            assert_eq!(sort6(d), want, "{bits:06b}");
        }
        let mut rng = crate::test_support::Rng::new(4);
        for _ in 0..1000 {
            let d: [f64; 6] = core::array::from_fn(|_| (rng.uniform() * 4.0).floor());
            let mut want = d;
            want.sort_by(|a, b| b.total_cmp(a));
            assert_eq!(sort6(d), want);
        }
    }

    #[test]
    fn mod_pi_is_rem_euclid() {
        use core::f64::consts::PI;
        let mut rng = crate::test_support::Rng::new(8);
        let edges = [
            0.0,
            -0.0,
            PI,
            -PI,
            PI / 2.0,
            -PI / 2.0,
            1e-300,
            -1e-300,
            3.0 * PI,
        ];
        let random = (0..5000).map(|_| rng.range(-PI, PI));
        for a in edges.into_iter().chain(random) {
            assert_eq!(mod_pi(a).to_bits(), a.rem_euclid(PI).to_bits(), "{a}");
        }
    }

    /// The neighbour list a scan of every star in index order keeps (the search
    /// `XOrder::nearest` replaced).
    fn nearest_by_full_scan(stars: &StarList, i: usize, k: usize) -> (Vec<f64>, Vec<usize>) {
        let mut idx = vec![0usize; k];
        let mut dist = vec![f64::MAX; k];
        idx[0] = i;
        dist[0] = 0.0;
        for (j, sj) in stars.0.iter().enumerate() {
            if j == i {
                continue;
            }
            let d = (sj.x - stars.0[i].x).powi(2) + (sj.y - stars.0[i].y).powi(2);
            if d <= 1.0 || d >= dist[k - 1] {
                continue;
            }
            let mut pos = k - 1;
            while pos > 0 && d < dist[pos - 1] {
                pos -= 1;
            }
            dist.insert(pos, d);
            idx.insert(pos, j);
            dist.truncate(k);
            idx.truncate(k);
        }
        (dist, idx)
    }

    #[test]
    fn the_x_ordered_neighbour_search_keeps_what_a_full_scan_keeps() {
        let mut rng = crate::test_support::Rng::new(9);
        for case in 0..30 {
            let n = 1 + (rng.next_u64() % 150) as usize;
            // Integer coordinates on a small grid give plenty of exact ties in
            // distance and in x, and some stars on top of each other.
            let side = if case % 2 == 0 { 12.0 } else { 1000.0 };
            let stars = StarList(
                (0..n)
                    .map(|_| Star {
                        x: (rng.uniform() * side).floor(),
                        y: (rng.uniform() * side).floor(),
                        snr: 1.0,
                        hfd: 2.0,
                    })
                    .collect(),
            );
            let order = XOrder::new(&stars);
            for k in [4, 7, 9] {
                for i in 0..n {
                    let mut dist = vec![0.0; k];
                    let mut idx = vec![0; k];
                    order.nearest(&stars, i, &mut dist, &mut idx);
                    let (want_d, want_i) = nearest_by_full_scan(&stars, i, k);
                    assert_eq!(
                        (dist, idx),
                        (want_d, want_i),
                        "case {case}, k {k}, star {i}"
                    );
                }
            }
        }
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
        for (k, r) in q.ratios.iter().enumerate().skip(1) {
            assert!((r - 1.0 / 2.0_f64.sqrt()).abs() < 1e-10, "ratio[{k}] = {r}");
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

    use crate::quads::r#match::{
        extract_star_pairs, filter_by_scale, find_matches_sorted, sort_catalog_quads,
    };
    use crate::quads::vote::vote_filter;
    use crate::test_support::Rng;
    use core::f64::consts::PI;

    /// `(x, y) → s·R(r)·(±x, y) + t`: a similarity, mirrored when `flip`.
    fn similarity(p: &[(f64, f64)], s: f64, r: f64, flip: bool) -> Vec<(f64, f64)> {
        let m = if flip { -1.0 } else { 1.0 };
        p.iter()
            .map(|&(x, y)| {
                let x = m * x;
                (
                    s * (x * r.cos() - y * r.sin()) + 250.0,
                    s * (x * r.sin() + y * r.cos()) - 90.0,
                )
            })
            .collect()
    }

    fn random_points(rng: &mut Rng, n: usize, side: f64) -> Vec<(f64, f64)> {
        (0..n)
            .map(|_| (rng.range(0.0, side), rng.range(0.0, side)))
            .collect()
    }

    /// The doc comment claims a maximum error of ~0.0026 rad; the approximation
    /// `z / (1 + 0.28125 z²)` actually peaks near 0.0049 rad (0.28°) around
    /// |z| ≈ 0.6. Either is negligible against the 10° vote bins, so this pins the
    /// real bound rather than the documented one.
    #[test]
    fn fast_atan2_is_accurate_in_every_quadrant() {
        let mut worst = 0.0f64;
        for k in 0..20_000 {
            let t = -PI + 2.0 * PI * k as f64 / 20_000.0;
            for rad in [0.01, 1.0, 1e4] {
                let (y, x) = (rad * t.sin(), rad * t.cos());
                worst = worst.max((fast_atan2(y, x) - y.atan2(x)).abs());
            }
        }
        assert!(worst < 5e-3, "max error {worst}");
        assert!(worst > 4e-3, "better than expected: {worst}");
        assert_eq!(fast_atan2(0.0, 0.0), 0.0);
        assert_eq!(fast_atan2(1.0, 0.0), core::f64::consts::FRAC_PI_2);
        assert_eq!(fast_atan2(-1.0, 0.0), -core::f64::consts::FRAC_PI_2);
    }

    /// Ratios are invariant under rotation, scale, translation and reflection, and
    /// under any ordering of the four stars; `d1` scales and `d1_angle` rotates.
    #[test]
    fn quad_fingerprint_is_similarity_invariant() {
        let mut rng = Rng::new(2);
        for trial in 0..300 {
            let p = random_points(&mut rng, 4, 100.0);
            let (s, r, flip) = (rng.range(0.05, 20.0), rng.range(-PI, PI), trial % 2 == 1);
            let q = similarity(&p, s, r, flip);
            let a = make_quad(p[0], p[1], p[2], p[3]).unwrap();
            let b = make_quad(q[3], q[1], q[0], q[2]).unwrap();
            for k in 0..5 {
                assert!((a.ratios[k] - b.ratios[k]).abs() < 1e-9, "ratio {k}");
            }
            assert!((b.d1 - s * a.d1).abs() < 1e-6 * b.d1);
            let (cx, cy) = similarity(&[(a.center_x, a.center_y)], s, r, flip)[0];
            assert!((b.center_x - cx).abs() < 1e-6 && (b.center_y - cy).abs() < 1e-6);
            // Angle of the longest side, mod π, turns with the image, to within two
            // fast_atan2 errors. Under a flip it is reflected instead.
            let want = if flip {
                (r + PI - a.d1_angle).rem_euclid(PI)
            } else {
                (a.d1_angle + r).rem_euclid(PI)
            };
            let d = (b.d1_angle - want).rem_euclid(PI);
            assert!(d.min(PI - d) < 1e-2, "angle {} vs {want}", b.d1_angle);
        }
    }

    /// The pattern chain the solver runs — build quads on both sides, match on
    /// ratios, vote on scale and rotation, pair the centres — recovers a known
    /// transform from two star lists that share most, but not all, of their stars.
    #[test]
    fn quad_pipeline_recovers_a_similarity_for_either_parity() {
        let mut rng = Rng::new(3);
        for flip in [false, true] {
            let shared = random_points(&mut rng, 80, 1000.0);
            let (s, r) = (2.7, 0.9);
            let mut img = shared.clone();
            img.extend(random_points(&mut rng, 8, 1000.0)); // image-only stars
            let mut cat = similarity(&shared, s, r, flip);
            cat.extend(similarity(&random_points(&mut rng, 8, 1000.0), s, r, flip));
            cat.sort_by(|a, b| a.0.total_cmp(&b.0));

            let img_q = build_quads(&make_stars(&img), img.len());
            let mut cat_q = build_quads_presorted(&make_stars(&cat), img.len());
            sort_catalog_quads(&mut cat_q);
            let raw = find_matches_sorted(&img_q, &cat_q, 0.005);
            let voted = vote_filter(&img_q, &cat_q, &raw, 0.005);
            assert!(
                voted.len() > 50,
                "{} votes of {} raw",
                voted.len(),
                raw.len()
            );
            assert!(voted.iter().all(|m| (m.scale_ratio * s - 1.0).abs() < 0.01));
            let (by_scale, med) = filter_by_scale(&raw, 0.005);
            assert!((med * s - 1.0).abs() < 1e-3 && !by_scale.is_empty());

            let (ip, cp) = extract_star_pairs(&img_q, &cat_q, &voted);
            let plate = crate::math::lsq::solve_plate_constants(&ip, &cp).unwrap();
            let m = if flip { -1.0 } else { 1.0 };
            let want = [m * s * r.cos(), -s * r.sin(), m * s * r.sin(), s * r.cos()];
            let got = [plate.a, plate.b, plate.d, plate.e];
            for k in 0..4 {
                assert!(
                    (got[k] - want[k]).abs() < 1e-3,
                    "flip {flip}: {got:?} vs {want:?}"
                );
            }
            assert!((plate.c - 250.0).abs() < 0.5 && (plate.f + 90.0).abs() < 0.5);
        }
    }

    #[test]
    fn builder_mode_follows_the_star_count() {
        let mut rng = Rng::new(4);
        let p = random_points(&mut rng, 40, 500.0);
        let stars = make_stars(&p);
        // Under 15 image stars: all C(7,4) quads per star (deduplicated).
        let q7 = build_quads(&make_stars(&p[..12]), 12);
        assert!(!q7.is_empty() && q7.len() <= 12 * 35);
        // 15..30: C(6,4) per star.
        let q6 = build_quads(&make_stars(&p[..20]), 20);
        assert!(!q6.is_empty() && q6.len() <= 20 * 15);
        // Otherwise the 9-neighbour builder.
        let q9 = build_quads(&stars, 40);
        assert!(q9.len() > q6.len());
        // Four stars or fewer fall back to 3-NN, which needs at least four.
        assert_eq!(build_quads(&make_stars(&p[..4]), 4).len(), 1);
        assert!(build_quads(&make_stars(&p[..3]), 3).is_empty());
        // A 5..8 star list at the 9-neighbour stage has too few stars for it.
        assert!(build_quads(&make_stars(&p[..6]), 40).is_empty());
        for q in q7.0.iter().chain(&q6.0).chain(&q9.0) {
            assert_ratios_valid(q);
        }
    }
}
