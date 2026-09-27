//! Vote-accumulator quad filter.
//!
//! After the standard O(M×N) ratio test finds raw quad matches, each match
//! implies a plate scale and a rotation angle.  True matches all vote for the
//! same (scale, rotation) cell; false matches scatter across the grid.
//!
//! Implementation: 2-D `HashMap` accumulator keyed by
//! `(round(scale / SCALE_STEP), round(angle_diff / ANGLE_STEP))`
//! so no explicit range assumption is needed.  The peak cell's matches are
//! returned, then optionally narrowed by scale with the existing median filter.

use core::f64::consts::PI;
use std::collections::HashMap;

use super::r#match::{QuadMatch, filter_by_scale};
use crate::types::QuadList;

/// Scale bin width: 5 % of plate scale.
const SCALE_STEP: f64 = 0.05;
/// Angle bin width: 10° (in radians).
const ANGLE_STEP: f64 = PI / 18.0;

/// Replace `filter_by_scale` with a 2-D (scale × rotation) vote accumulator.
///
/// Each raw match votes for `(scale_bin, angle_bin)`.  The bin with the most
/// votes is the winning transformation; its matches are returned.  A subsequent
/// median-scale pass removes any residual outliers.
///
/// Two angle-key formulas are tried in parallel:
///   • **Direct** (det > 0, pure rotation):    key = `(φ_img − φ_cat) mod π`
///   • **Reflected** (det < 0, FITS convention): key = `(φ_img + φ_cat) mod π`
///
/// For a direct transform, all true matches share the same diff-key.
/// For a reflected transform (which is the common FITS case — CDELT1 is
/// negative so RA increases opposite to pixel x), the diff-key scatters
/// while the sum-key is constant.  By searching both grids we handle both
/// cases without knowing the handedness in advance.
#[must_use]
pub fn vote_filter(
    img_quads: &QuadList,
    cat_quads: &QuadList,
    raw_matches: &[QuadMatch],
    quad_tolerance: f64,
) -> Vec<QuadMatch> {
    if raw_matches.is_empty() {
        return vec![];
    }

    let n_angle_bins = (PI / ANGLE_STEP).round() as i32;

    // Two grids: one for direct (rotation-only) transforms, one for reflected.
    let mut grid_direct: HashMap<(i32, i32), Vec<usize>> = HashMap::new();
    let mut grid_reflected: HashMap<(i32, i32), Vec<usize>> = HashMap::new();

    for (idx, m) in raw_matches.iter().enumerate() {
        let iq = &img_quads.0[m.img_idx];
        let cq = &cat_quads.0[m.cat_idx];

        let scale_key = (m.scale_ratio / SCALE_STEP).round() as i32;

        // Direct: consistent when the pixel→catalog transform is a pure rotation.
        let diff = (iq.d1_angle - cq.d1_angle).rem_euclid(PI);
        let diff_key = (diff / ANGLE_STEP).round() as i32 % n_angle_bins;

        // Reflected: consistent when the transform includes a handedness flip
        // (det < 0).  φ_cat = (C − φ_img) mod π  ⟹  φ_img + φ_cat = C (const).
        let sum = (iq.d1_angle + cq.d1_angle).rem_euclid(PI);
        let sum_key = (sum / ANGLE_STEP).round() as i32 % n_angle_bins;

        grid_direct
            .entry((scale_key, diff_key))
            .or_default()
            .push(idx);
        grid_reflected
            .entry((scale_key, sum_key))
            .or_default()
            .push(idx);
    }

    // Pick the larger peak across both grids.
    //
    // Ties are broken on the bin key, not on iteration order. HashMap iteration
    // is seeded randomly per process, so `max_by_key` on the count alone returned
    // whichever tied bin happened to come first — and ties are commonest exactly
    // where there are few raw matches, i.e. in the marginal solves. The same
    // image could then solve on one run and not the next.
    let peak_direct = grid_direct
        .iter()
        .max_by_key(|(k, v)| (v.len(), **k))
        .map(|(_, v)| v.clone())
        .unwrap_or_default();
    let peak_reflected = grid_reflected
        .iter()
        .max_by_key(|(k, v)| (v.len(), **k))
        .map(|(_, v)| v.clone())
        .unwrap_or_default();

    let peak_indices = if peak_reflected.len() >= peak_direct.len() {
        peak_reflected
    } else {
        peak_direct
    };

    // Sorted so the downstream least-squares accumulation sees a fixed order:
    // the GIVENS sweep is float-order-dependent, so even a won bin could produce
    // a slightly different fit run to run.
    let mut peak_indices = peak_indices;
    peak_indices.sort_unstable();
    let peak_matches: Vec<QuadMatch> = peak_indices.iter().map(|&i| raw_matches[i]).collect();

    // Secondary scale pass to remove edge-of-bin outliers
    let (refined, _) = filter_by_scale(&peak_matches, quad_tolerance);
    if refined.is_empty() {
        peak_matches
    } else {
        refined
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{Quad, QuadList};

    fn quad(d1: f64, ratios: [f64; 5], cx: f64, cy: f64, angle: f64) -> Quad {
        Quad {
            d1,
            ratios,
            center_x: cx,
            center_y: cy,
            d1_angle: angle,
        }
    }

    #[test]
    fn peak_bin_selected() {
        // 3 matches at scale≈2, angle≈0.5 rad; 1 false match at scale≈5, angle≈2.0 rad
        let r = [0.9, 0.8, 0.7, 0.6, 0.5];
        let img = QuadList(vec![
            quad(200.0, r, 50.0, 50.0, 0.52),
            quad(202.0, r, 60.0, 60.0, 0.50),
            quad(198.0, r, 70.0, 70.0, 0.54),
            quad(500.0, r, 80.0, 80.0, 2.00), // false
        ]);
        let cat = QuadList(vec![
            quad(100.0, r, 0.1, 0.1, 0.52),
            quad(101.0, r, 0.2, 0.2, 0.50),
            quad(99.0, r, 0.3, 0.3, 0.54),
            quad(100.0, r, 0.4, 0.4, 2.00),
        ]);
        let raw = vec![
            QuadMatch {
                img_idx: 0,
                cat_idx: 0,
                scale_ratio: 2.00,
            },
            QuadMatch {
                img_idx: 1,
                cat_idx: 1,
                scale_ratio: 2.00,
            },
            QuadMatch {
                img_idx: 2,
                cat_idx: 2,
                scale_ratio: 2.00,
            },
            QuadMatch {
                img_idx: 3,
                cat_idx: 3,
                scale_ratio: 5.00,
            }, // outlier
        ];
        let filtered = vote_filter(&img, &cat, &raw, 0.1);
        // Should keep the 3 true matches, discard the false one
        assert_eq!(
            filtered.len(),
            3,
            "expected 3 matches, got {}",
            filtered.len()
        );
        for m in &filtered {
            assert!(
                (m.scale_ratio - 2.0).abs() < 0.1,
                "unexpected scale {}",
                m.scale_ratio
            );
        }
    }

    /// 30 true matches (one scale, one rotation) hidden among 300 random ones, for
    /// a direct and for a reflected transform: the vote returns the true ones and
    /// at most a stray coincidence.
    #[test]
    fn vote_finds_the_true_transform_among_many_false_matches() {
        let mut rng = crate::test_support::Rng::new(8);
        let r = [0.9, 0.8, 0.7, 0.6, 0.5];
        for reflected in [false, true] {
            let (scale, rot) = (3.3, 1.2);
            let mut img = Vec::new();
            let mut cat = Vec::new();
            let mut raw = Vec::new();
            for k in 0..330 {
                let cat_angle = rng.range(0.0, PI);
                let is_true = k % 11 == 0;
                let (s, img_angle) = if is_true {
                    let a = if reflected {
                        (rot - cat_angle).rem_euclid(PI)
                    } else {
                        (cat_angle + rot).rem_euclid(PI)
                    };
                    (scale * rng.range(0.999, 1.001), a)
                } else {
                    (rng.range(0.5, 10.0), rng.range(0.0, PI))
                };
                cat.push(quad(10.0, r, 0.0, 0.0, cat_angle));
                img.push(quad(10.0 * s, r, 0.0, 0.0, img_angle));
                raw.push(QuadMatch {
                    img_idx: k,
                    cat_idx: k,
                    scale_ratio: s,
                });
            }
            let out = vote_filter(&QuadList(img), &QuadList(cat), &raw, 0.007);
            let n_true = out.iter().filter(|m| m.img_idx % 11 == 0).count();
            assert_eq!(n_true, 30, "reflected {reflected}");
            assert!(out.len() <= 31, "reflected {reflected}: {} kept", out.len());
            // Output is in input order, so the downstream fit is deterministic.
            assert!(out.windows(2).all(|w| w[0].img_idx < w[1].img_idx));
        }
    }

    #[test]
    fn empty_input_returns_empty() {
        let img = QuadList(vec![]);
        let cat = QuadList(vec![]);
        let result = vote_filter(&img, &cat, &[], 0.007);
        assert!(result.is_empty());
    }

    #[test]
    fn single_match_returned() {
        let r = [0.9, 0.8, 0.7, 0.6, 0.5];
        let img = QuadList(vec![quad(100.0, r, 10.0, 10.0, 1.0)]);
        let cat = QuadList(vec![quad(50.0, r, 0.1, 0.1, 1.0)]);
        let raw = vec![QuadMatch {
            img_idx: 0,
            cat_idx: 0,
            scale_ratio: 2.0,
        }];
        let result = vote_filter(&img, &cat, &raw, 0.1);
        assert_eq!(result.len(), 1);
    }
}
