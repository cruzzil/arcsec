//! Quad pattern matching.

use crate::types::{PairedPositions, Quad, QuadList};

/// A pair of matching quads: image index and catalog index.
#[derive(Debug, Clone, Copy)]
pub struct QuadMatch {
    /// Index into the image `QuadList`.
    pub img_idx: usize,
    /// Index into the catalogue `QuadList`.
    pub cat_idx: usize,
    /// `d1_img` / `d1_cat` — the pixel/catalog scale ratio for this match
    pub scale_ratio: f64,
}

/// Compare two quads: all 5 ratios must agree within `tolerance`.
fn ratios_match(img: &Quad, cat: &Quad, tolerance: f64) -> bool {
    for k in 0..5 {
        if (img.ratios[k] - cat.ratios[k]).abs() > tolerance {
            return false;
        }
    }
    true
}

/// Find all quad pairs where the 5 normalized ratios agree within `quad_tolerance`.
///
/// Returns a list of `QuadMatch` with raw scale ratios.
#[must_use]
pub fn find_matches(
    img_quads: &QuadList,
    cat_quads: &QuadList,
    quad_tolerance: f64,
) -> Vec<QuadMatch> {
    let mut matches = Vec::new();

    for (i, iq) in img_quads.0.iter().enumerate() {
        for (j, cq) in cat_quads.0.iter().enumerate() {
            if cq.d1 < 1e-10 {
                continue;
            }
            if ratios_match(iq, cq, quad_tolerance) {
                matches.push(QuadMatch {
                    img_idx: i,
                    cat_idx: j,
                    scale_ratio: iq.d1 / cq.d1,
                });
            }
        }
    }

    matches
}

/// Which of the five ratios the catalogue is sorted and binary-searched on.
///
/// The ratios are `d2/d1 .. d6/d1` with `d` sorted descending, so `ratios[0]` is the
/// one closest to 1 and the most tightly clustered - the worst possible index key,
/// because the +/- tolerance window around it catches a large slice of the catalogue.
/// `ratios[4] = d6/d1` is the smallest and most widely spread, so its window is far
/// narrower and fewer candidates need the full five-ratio check.
pub const INDEX_RATIO: usize = 4;

/// Sort catalogue quads into the order [`find_matches_sorted`] requires.
///
/// Call this rather than spelling the sort out: which ratio indexes the search is
/// an implementation detail, and getting it wrong drops matches without any error.
pub fn sort_catalog_quads(cat_quads: &mut QuadList) {
    cat_quads
        .0
        .sort_unstable_by(|a, b| a.ratios[INDEX_RATIO].total_cmp(&b.ratios[INDEX_RATIO]));
}

/// Find all quad pairs using binary search on catalog quads pre-sorted by `ratios[INDEX_RATIO]`.
///
/// `cat_quads` must be sorted ascending by `ratios[INDEX_RATIO]` before calling —
/// use [`sort_catalog_quads`], which is the only supported way to establish it.
/// For each image quad, binary-searches to find only the catalog quads in the
/// `ratios[INDEX_RATIO]` window, then checks the remaining 4 ratios — reducing
/// O(n·m) to O(n·(log m + hits)).
///
/// Sorting by the wrong ratio fails silently: `partition_point` on unsorted data
/// returns an arbitrary split and matches are simply dropped.
#[must_use]
pub fn find_matches_sorted(
    img_quads: &QuadList,
    cat_quads: &QuadList,
    quad_tolerance: f64,
) -> Vec<QuadMatch> {
    let codes = CatalogCodes::build(cat_quads);
    find_matches_indexed(img_quads, cat_quads, &codes, quad_tolerance)
}

/// Compact copy of the catalogue quads' ratios for cache-efficient scanning.
///
/// `Quad` is 72 bytes (nine f64s), so scanning the tolerance window touches far more
/// memory than the comparison needs. Storing the five ratios as f32 in a parallel
/// array is 20 bytes per quad, so several times more candidates fit in cache - the
/// same trick `catalog::anet` already uses for its 9 MB code array. f32 has ~7
/// significant digits, which is an order of magnitude finer than the 0.007 default
/// tolerance, so the narrowing is harmless; the surviving candidates are re-checked
/// against the full-precision f64 ratios anyway.
pub struct CatalogCodes {
    ratios: Vec<[f32; 5]>,
}

impl CatalogCodes {
    /// `cat_quads` must already be sorted ascending by `ratios[INDEX_RATIO]`.
    #[must_use]
    pub fn build(cat_quads: &QuadList) -> Self {
        let ratios = cat_quads
            .0
            .iter()
            .map(|q| {
                [
                    q.ratios[0] as f32,
                    q.ratios[1] as f32,
                    q.ratios[2] as f32,
                    q.ratios[3] as f32,
                    q.ratios[4] as f32,
                ]
            })
            .collect();
        Self { ratios }
    }
}

/// As `find_matches_sorted`, but reusing a prebuilt [`CatalogCodes`].
#[must_use]
pub fn find_matches_indexed(
    img_quads: &QuadList,
    cat_quads: &QuadList,
    codes: &CatalogCodes,
    quad_tolerance: f64,
) -> Vec<QuadMatch> {
    let cat = &cat_quads.0;
    let ratios = &codes.ratios;
    debug_assert_eq!(cat.len(), ratios.len());
    let tol = quad_tolerance as f32;
    // f32 rounding can move a value by up to ~1e-7 relative; widen the window by a
    // hair so a borderline true match is never dropped before the f64 re-check.
    let tol_pad = tol * 1.000_01 + f32::EPSILON;
    let mut matches = Vec::new();

    for (i, iq) in img_quads.0.iter().enumerate() {
        let key = iq.ratios[INDEX_RATIO] as f32;
        let lo = key - tol_pad;
        let hi = key + tol_pad;

        let start = ratios.partition_point(|r| r[INDEX_RATIO] < lo);
        let ir = [
            iq.ratios[0] as f32,
            iq.ratios[1] as f32,
            iq.ratios[2] as f32,
            iq.ratios[3] as f32,
            iq.ratios[4] as f32,
        ];

        for (offset, cr) in ratios[start..].iter().enumerate() {
            if cr[INDEX_RATIO] > hi {
                break;
            }
            // Cheap f32 pass over the compact array; only survivors touch `Quad`.
            let mut ok = true;
            for k in 0..5 {
                if (ir[k] - cr[k]).abs() > tol_pad {
                    ok = false;
                    break;
                }
            }
            if !ok {
                continue;
            }

            let j = start + offset;
            let cq = &cat[j];
            if cq.d1 < 1e-10 {
                continue;
            }
            // Re-check at full precision.
            let mut exact = true;
            for k in 0..5 {
                if (iq.ratios[k] - cq.ratios[k]).abs() > quad_tolerance {
                    exact = false;
                    break;
                }
            }
            if exact {
                matches.push(QuadMatch {
                    img_idx: i,
                    cat_idx: j,
                    scale_ratio: iq.d1 / cq.d1,
                });
            }
        }
    }

    matches
}

/// Compute the median of a slice (sorts a copy).
pub(crate) fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut v = values.to_vec();
    v.sort_by(f64::total_cmp);
    let mid = v.len() / 2;
    if v.len().is_multiple_of(2) {
        (v[mid - 1] + v[mid]) * 0.5
    } else {
        v[mid]
    }
}

/// Filter `matches` by removing scale outliers.
///
/// Computes the median `d1_img / d1_cat` ratio, then retains only matches
/// where `|ratio - median| <= quad_tolerance * median`.
///
/// Returns `(filtered_matches, median_ratio)`.
#[must_use]
pub fn filter_by_scale(matches: &[QuadMatch], quad_tolerance: f64) -> (Vec<QuadMatch>, f64) {
    if matches.is_empty() {
        return (vec![], 0.0);
    }

    let ratios: Vec<f64> = matches.iter().map(|m| m.scale_ratio).collect();
    let med = median(&ratios);

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

/// Extract (`image_xy`, `catalog_xy`) position pairs from matching quads.
///
/// The image coordinates are quad centre positions (pixels).
/// The catalog coordinates are quad centre positions in whatever space the catalog
/// quads were built in (standard projection coordinates when calling from the pipeline).
///
/// Returns `(img_positions, cat_positions)` suitable for `solve_plate_constants`.
#[must_use]
pub fn extract_star_pairs(
    img_quads: &QuadList,
    cat_quads: &QuadList,
    matches: &[QuadMatch],
) -> PairedPositions {
    let mut img_pos = Vec::with_capacity(matches.len());
    let mut cat_pos = Vec::with_capacity(matches.len());
    for m in matches {
        let iq = &img_quads.0[m.img_idx];
        let cq = &cat_quads.0[m.cat_idx];
        img_pos.push((iq.center_x, iq.center_y));
        cat_pos.push((cq.center_x, cq.center_y));
    }
    (img_pos, cat_pos)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Quad;

    fn make_quad(d1: f64, ratios: [f64; 5], cx: f64, cy: f64) -> Quad {
        Quad {
            d1,
            ratios,
            center_x: cx,
            center_y: cy,
            d1_angle: 0.0,
        }
    }

    #[test]
    fn identical_quads_match() {
        let r = [0.9, 0.8, 0.7, 0.6, 0.5];
        let img = QuadList(vec![make_quad(100.0, r, 50.0, 50.0)]);
        let cat = QuadList(vec![make_quad(50.0, r, 0.01, 0.02)]);
        let matches = find_matches(&img, &cat, 0.01);
        assert_eq!(matches.len(), 1);
        assert!((matches[0].scale_ratio - 2.0).abs() < 1e-10);
    }

    #[test]
    fn sorted_matches_agree_with_linear() {
        // Build several image and catalog quads; verify sorted variant finds same matches.
        let img_rs = [
            [0.95, 0.85, 0.75, 0.65, 0.55],
            [0.80, 0.70, 0.60, 0.50, 0.40],
            [0.60, 0.50, 0.40, 0.30, 0.20],
        ];
        let cat_rs = [
            [0.951, 0.849, 0.751, 0.651, 0.549], // matches img[0]
            [0.70, 0.60, 0.50, 0.40, 0.30],      // no match
            [0.801, 0.699, 0.601, 0.499, 0.401], // matches img[1]
            [0.60, 0.50, 0.40, 0.30, 0.21],      // close to img[2] but last ratio off
        ];
        let img = QuadList(
            img_rs
                .iter()
                .map(|&r| make_quad(100.0, r, 0.0, 0.0))
                .collect(),
        );
        let cat_quads: Vec<Quad> = cat_rs
            .iter()
            .map(|&r| make_quad(50.0, r, 0.0, 0.0))
            .collect();
        let tol = 0.005;

        // Linear reference
        let cat_unsorted = QuadList(cat_quads.clone());
        let mut linear = find_matches(&img, &cat_unsorted, tol);
        linear.sort_by_key(|m| (m.img_idx, m.cat_idx));

        // Sorted variant — sort by ratios[0] ascending first
        // Must be INDEX_RATIO, not ratios[0]: this fixture happens to be
        // co-monotonic in both, so sorting by the wrong one still passed.
        let mut cat_quads = QuadList(cat_quads);
        sort_catalog_quads(&mut cat_quads);
        let cat_quads = cat_quads.0;
        let cat_sorted = QuadList(cat_quads);
        let mut sorted = find_matches_sorted(&img, &cat_sorted, tol);
        sorted.sort_by_key(|m| (m.img_idx, m.scale_ratio.to_bits()));

        assert_eq!(linear.len(), sorted.len(), "match counts must agree");
    }

    #[test]
    fn tolerant_match() {
        let r_img = [0.9, 0.8, 0.7, 0.6, 0.5];
        let r_cat = [0.901, 0.799, 0.701, 0.601, 0.499];
        let img = QuadList(vec![make_quad(100.0, r_img, 50.0, 50.0)]);
        let cat = QuadList(vec![make_quad(50.0, r_cat, 0.01, 0.02)]);
        // Tolerance 0.002: all diffs are 0.001, should match
        let matches = find_matches(&img, &cat, 0.002);
        assert_eq!(matches.len(), 1);
        // Tolerance 0.0005: all diffs are 0.001, should NOT match
        let no_matches = find_matches(&img, &cat, 0.0005);
        assert_eq!(no_matches.len(), 0);
    }

    #[test]
    fn no_match_on_ratio_mismatch() {
        let r1 = [0.9, 0.8, 0.7, 0.6, 0.5];
        let r2 = [0.9, 0.8, 0.7, 0.6, 0.3]; // last ratio differs by 0.2
        let img = QuadList(vec![make_quad(100.0, r1, 50.0, 50.0)]);
        let cat = QuadList(vec![make_quad(50.0, r2, 0.01, 0.02)]);
        let matches = find_matches(&img, &cat, 0.01);
        assert_eq!(matches.len(), 0);
    }

    #[test]
    fn filter_removes_scale_outliers() {
        // Two matches with scale 2.0, one outlier with scale 4.0
        let matches = vec![
            QuadMatch {
                img_idx: 0,
                cat_idx: 0,
                scale_ratio: 2.0,
            },
            QuadMatch {
                img_idx: 1,
                cat_idx: 1,
                scale_ratio: 2.05,
            },
            QuadMatch {
                img_idx: 2,
                cat_idx: 2,
                scale_ratio: 4.0,
            },
        ];
        let (filtered, med) = filter_by_scale(&matches, 0.1);
        // median should be ~2.0 or 2.05
        assert!(med > 1.9 && med < 2.1, "median = {med}");
        // Outlier at 4.0 is > 10% from median of ~2.0 → filtered out
        assert_eq!(filtered.len(), 2, "expected 2 matches after filtering");
    }

    #[test]
    fn filter_empty_input() {
        let (filtered, med) = filter_by_scale(&[], 0.1);
        assert!(filtered.is_empty());
        assert_eq!(med, 0.0);
    }

    #[test]
    fn extract_pairs_positions() {
        let r = [0.9, 0.8, 0.7, 0.6, 0.5];
        let img = QuadList(vec![make_quad(100.0, r, 50.0, 60.0)]);
        let cat = QuadList(vec![make_quad(50.0, r, 0.1, 0.2)]);
        let matches = vec![QuadMatch {
            img_idx: 0,
            cat_idx: 0,
            scale_ratio: 2.0,
        }];
        let (ip, cp) = extract_star_pairs(&img, &cat, &matches);
        assert_eq!(ip, vec![(50.0, 60.0)]);
        assert_eq!(cp, vec![(0.1, 0.2)]);
    }

    #[test]
    fn median_even_count() {
        let v = [1.0, 3.0, 5.0, 7.0];
        assert_eq!(median(&v), 4.0);
    }

    #[test]
    fn median_odd_count() {
        let v = [1.0, 3.0, 5.0];
        assert_eq!(median(&v), 3.0);
    }
}
