//! Searching pixel scales: what the catalogue search does when it does not know
//! the image's scale, or (if asked) when the scale it was given does not solve.
//!
//! The spiral steps one field at a time and reads a catalogue window one field
//! wide, both derived from the pixel scale, so a scale wrong by about 2× reads the
//! wrong stars at every position: too few to match an image that is really wider,
//! too many, and too bright, for one that is really narrower. The quads themselves
//! are scale-free; what the scale buys is a catalogue that looks like the image.
//!
//! So a scale search is a ladder of hypotheses, √2 apart, each a short catalogue
//! search round the hint ([`LADDER_FIELDS`]); a true scale lies within 2^¼ (19 %) of
//! one of them, and the measured tolerance of a single solve is wider than that
//! (docs/plate-solving.md §10.3e).

use core::f64::consts::SQRT_2;

/// Ratio between neighbouring hypotheses: astrometry.net's index scales step by
/// the same factor.
pub const SCALE_STEP: f64 = SQRT_2;

/// Steps of [`SCALE_STEP`] below and above 1″/px searched when nothing gives the
/// scale: 0.25″/px to 64″/px, from a long focal length binned to a camera lens. The
/// blind index searches 0.3–60″/px for the same case.
pub const UNKNOWN_STEPS: (i32, i32) = (-4, 12);

/// Steps either side of a scale that was given (or read from the header) and did
/// not solve, with [`ScaleSearch::AlsoIfWrong`]: a quarter to four times it.
pub const WRONG_STEPS: i32 = 4;

/// How far round the hint each hypothesis is searched, in fields of that
/// hypothesis: the hint's own position and the ring of eight round it. The search
/// at the assumed scale then continues out to `-r`, as it always has.
pub const LADDER_FIELDS: f64 = 1.0;

/// When the catalogue search tries other pixel scales.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum ScaleSearch {
    /// Never: an unknown scale is taken to be 1″/px, as before 0.6.
    Never,
    /// When neither a field of view nor header optics give the scale (ASTAP's
    /// `-fov 0` with no FOCALLEN/XPIXSZ): a ladder of scales round 1″/px, then the
    /// usual search at 1″/px. The default.
    #[default]
    IfUnknown,
    /// As [`Self::IfUnknown`], and also when a solve at a known scale finds
    /// nothing: a ladder from a quarter to four times that scale (the CLI's
    /// `--fov-search`). Off by default, so a correctly configured solve that fails
    /// costs no more than it did.
    AlsoIfWrong,
}

/// One scale hypothesis: `centre × SCALE_STEP^step`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hypothesis {
    /// Steps of [`SCALE_STEP`] from the centre of the ladder.
    pub step: i32,
    /// Pixel scale, arcseconds per (unbinned) pixel.
    pub scale: f64,
}

/// The hypotheses `centre × √2^k` for `k` in `lo..=hi`, most likely first.
///
/// "Most likely" is nearest the centre (the scale assumed or given), in steps;
/// of two equally near, the larger scale first, because a wider field is the
/// cheaper search (fewer catalogue stars per degree, fewer positions) and the
/// common way to have no scale at all is a camera lens. `skip_centre` leaves the
/// centre out, for when it has already been searched.
#[must_use]
pub fn ladder(centre: f64, lo: i32, hi: i32, skip_centre: bool) -> Vec<Hypothesis> {
    let mut steps: Vec<i32> = (lo..=hi).filter(|&k| !(skip_centre && k == 0)).collect();
    steps.sort_by_key(|&k| (k.unsigned_abs(), k < 0));
    steps
        .into_iter()
        .map(|step| Hypothesis {
            step,
            scale: centre * SCALE_STEP.powi(step),
        })
        .collect()
}

/// The fraction by which a solved scale may differ from the one the solve started
/// from before the solution says so ([`inaccurate_scale_warning`]), as `astap_cli`.
pub const INACCURATE_SCALE: f64 = 0.05;

/// `astap_cli`'s warning when the solved scale is not the one it started from:
/// `Warning scale was inaccurate! Set FOV=1.00d, scale=1.2"`, printed after the
/// solution and written to the `.ini` as `WARNING`.
///
/// `start_scale` is the scale the solve started from and `solved_scale` the
/// solution's, both arcseconds per unbinned pixel; `height` is the unbinned image
/// height in pixels. The FOV is the solution's image height, as `-fov` takes it.
#[must_use]
pub fn inaccurate_scale_warning(
    start_scale: f64,
    solved_scale: f64,
    height: usize,
) -> Option<String> {
    let ratio = start_scale / solved_scale;
    if !ratio.is_finite() || (ratio - 1.0).abs() <= INACCURATE_SCALE {
        return None;
    }
    Some(format!(
        "Warning scale was inaccurate! Set FOV={:.2}d, scale={:.1}\"",
        solved_scale * height as f64 / 3600.0,
        solved_scale
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn steps(v: &[Hypothesis]) -> Vec<i32> {
        v.iter().map(|h| h.step).collect()
    }

    #[test]
    fn the_ladder_starts_at_the_centre_and_alternates_outwards() {
        let l = ladder(1.0, -2, 3, false);
        assert_eq!(steps(&l), [0, 1, -1, 2, -2, 3]);
        assert!((l[1].scale - SQRT_2).abs() < 1e-12);
        assert!((l[2].scale - 1.0 / SQRT_2).abs() < 1e-12);
        assert!((l[5].scale - 2.0 * SQRT_2).abs() < 1e-12);
    }

    #[test]
    fn the_unknown_ladder_covers_a_quarter_to_64_arcsec_per_pixel() {
        let (lo, hi) = UNKNOWN_STEPS;
        let l = ladder(1.0, lo, hi, false);
        assert_eq!(l.len(), 17);
        let min = l.iter().map(|h| h.scale).fold(f64::INFINITY, f64::min);
        let max = l.iter().map(|h| h.scale).fold(0.0, f64::max);
        assert!((min - 0.25).abs() < 1e-12 && (max - 64.0).abs() < 1e-9);
        // Past the smaller end, only larger scales remain, in order.
        assert_eq!(steps(&l)[8..], [-4, 5, 6, 7, 8, 9, 10, 11, 12]);
        // Every scale in the range is within 2^(1/4) of a hypothesis.
        let mut s = 0.25_f64;
        while s <= 64.0 {
            let nearest = l
                .iter()
                .map(|h| (h.scale / s).ln().abs())
                .fold(f64::INFINITY, f64::min);
            assert!(nearest <= SCALE_STEP.ln() / 2.0 + 1e-12, "{s}");
            s *= 1.01;
        }
    }

    #[test]
    fn a_given_scale_that_failed_is_not_tried_again() {
        let l = ladder(2.0, -WRONG_STEPS, WRONG_STEPS, true);
        assert_eq!(steps(&l), [1, -1, 2, -2, 3, -3, 4, -4]);
        assert!((l[6].scale - 8.0).abs() < 1e-9 && (l[7].scale - 0.5).abs() < 1e-12);
    }

    #[test]
    fn the_warning_follows_astap() {
        assert_eq!(inaccurate_scale_warning(1.25, 1.241, 2900), None);
        assert_eq!(
            inaccurate_scale_warning(1.0, 1.2857, 1400).as_deref(),
            Some("Warning scale was inaccurate! Set FOV=0.50d, scale=1.3\"")
        );
        assert!(inaccurate_scale_warning(1.31, 1.241, 2900).is_some());
        assert!(inaccurate_scale_warning(1.17, 1.241, 2900).is_some());
        assert_eq!(inaccurate_scale_warning(1.0, 0.0, 100), None);
    }
}
