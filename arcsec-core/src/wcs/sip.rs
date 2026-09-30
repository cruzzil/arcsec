//! SIP (Simple Imaging Polynomial) distortion: fitting, and pixel/sky transforms
//! that honour it.
//!
//! SIP (Shupe et al. 2005, ADASS XIV) adds a polynomial to the linear TAN WCS. With
//! `u = x - CRPIX1` and `v = y - CRPIX2`,
//!
//! ```text
//! (ξ, η) = CD · (u + A(u, v),  v + B(u, v))          pixel → sky
//! (u, v) = (U + AP(U, V),  V + BP(U, V)),  (U, V) = CD⁻¹ · (ξ, η)      sky → pixel
//! ```
//!
//! where `(ξ, η)` are the TAN standard coordinates in degrees and each polynomial is
//! `Σ c_pq · uᵖ vᵠ` over `p + q ≤ 3`.
//!
//! The fit follows ASTAP's `add_sip`, so the keywords arcsec writes carry the same
//! meaning as `astap_cli -sip`'s: a full cubic including the constant and linear
//! terms (they absorb what is left of the linear fit, since CD is not refitted), with
//! the inverse fitted separately from the same star pairs rather than by inverting
//! the forward polynomial. Where ASTAP fits quad centroids, arcsec fits the
//! individual verified stars ([`WcsSolution::matched_stars`]), of which there are
//! many more.

use crate::math::lsq::lsq_fit;
use crate::types::WcsSolution;

/// Polynomial order of the fitted distortion (`A_ORDER`, `B_ORDER`, `AP_ORDER`,
/// `BP_ORDER`).
pub const SIP_ORDER: u32 = 3;

/// The exponents `(p, q)` of the terms `uᵖ vᵠ`, in the order the coefficient arrays of
/// [`Sip`] hold them. It is also the order ASTAP writes the `A_p_q` keywords in.
pub const SIP_TERMS: [(i32, i32); 10] = [
    (0, 0),
    (1, 0),
    (0, 1),
    (2, 0),
    (1, 1),
    (0, 2),
    (3, 0),
    (2, 1),
    (1, 2),
    (0, 3),
];

/// Fewest star pairs a SIP fit is attempted from, as in ASTAP. Each axis has ten
/// coefficients.
pub const MIN_SIP_STARS: usize = 20;

/// The F statistic a cubic must reach over a linear fit to be kept.
///
/// With 14 extra parameters and a few hundred stars, chance alone exceeds 2.6
/// once in a thousand fits; 4 is stricter still, since a false distortion costs
/// accuracy while a real one (optics) is detected with F in the hundreds.
pub const MIN_SIP_F: f64 = 4.0;

/// The frame is divided into `COVERAGE_GRID` × `COVERAGE_GRID` cells for the
/// coverage test, and each must hold at least `COVERAGE_MIN` matched stars.
const COVERAGE_GRID: usize = 3;
/// See [`COVERAGE_GRID`].
const COVERAGE_MIN: usize = 2;

/// The largest correction a fit may make at an image corner, as a fraction of the
/// image half-diagonal. A cubic fitted to stars that do not reach the corners can
/// extrapolate wildly there; real optical distortion stays far below this.
const MAX_CORNER_CORRECTION: f64 = 0.05;

/// A star pair for fitting: its measured pixel offset from CRPIX, and the offset
/// the linear WCS puts its catalogue position at.
type Pair = ((f64, f64), (f64, f64));

/// Third-order SIP polynomials, coefficients indexed as [`SIP_TERMS`].
#[derive(Debug, Clone, PartialEq)]
pub struct Sip {
    /// `A_p_q`: pixel → sky correction along the first axis, pixels.
    pub a: [f64; 10],
    /// `B_p_q`: pixel → sky correction along the second axis, pixels.
    pub b: [f64; 10],
    /// `AP_p_q`: sky → pixel correction along the first axis, pixels.
    pub ap: [f64; 10],
    /// `BP_p_q`: sky → pixel correction along the second axis, pixels.
    pub bp: [f64; 10],
}

/// `Σ c_pq · uᵖ vᵠ`.
fn poly(c: &[f64; 10], u: f64, v: f64) -> f64 {
    SIP_TERMS
        .iter()
        .zip(c)
        .map(|(&(p, q), &k)| k * u.powi(p) * v.powi(q))
        .sum()
}

impl Sip {
    /// Pixel offsets from CRPIX → the linear (undistorted) pixel offsets that CD
    /// maps onto the sky.
    #[must_use]
    pub fn pixel_to_linear(&self, u: f64, v: f64) -> (f64, f64) {
        (u + poly(&self.a, u, v), v + poly(&self.b, u, v))
    }

    /// Linear pixel offsets (`CD⁻¹ · (ξ, η)`) → pixel offsets from CRPIX.
    #[must_use]
    pub fn linear_to_pixel(&self, u: f64, v: f64) -> (f64, f64) {
        (u + poly(&self.ap, u, v), v + poly(&self.bp, u, v))
    }
}

/// A TAN WCS with optional SIP distortion, for mapping between pixels and the sky.
#[derive(Debug, Clone, PartialEq)]
pub struct TanWcs {
    /// CRVAL1, radians.
    pub ra0: f64,
    /// CRVAL2, radians.
    pub dec0: f64,
    /// CRPIX1, 1-based FITS pixels.
    pub crpix1: f64,
    /// CRPIX2, 1-based FITS pixels.
    pub crpix2: f64,
    /// `[[CD1_1, CD1_2], [CD2_1, CD2_2]]`, degrees per pixel.
    pub cd: [[f64; 2]; 2],
    /// Distortion, if any.
    pub sip: Option<Sip>,
}

impl From<&WcsSolution> for TanWcs {
    fn from(w: &WcsSolution) -> Self {
        Self {
            ra0: w.ra0,
            dec0: w.dec0,
            crpix1: w.crpix1,
            crpix2: w.crpix2,
            cd: [[w.cd1_1, w.cd1_2], [w.cd2_1, w.cd2_2]],
            sip: w.sip.clone(),
        }
    }
}

impl TanWcs {
    /// 1-based FITS pixel → (RA, Dec) in radians, RA in `[0, 2π)`.
    #[must_use]
    pub fn pixel_to_sky(&self, x: f64, y: f64) -> (f64, f64) {
        let (u, v) = (x - self.crpix1, y - self.crpix2);
        let (u, v) = self
            .sip
            .as_ref()
            .map_or((u, v), |s| s.pixel_to_linear(u, v));
        let xi = (self.cd[0][0] * u + self.cd[0][1] * v).to_radians();
        let eta = (self.cd[1][0] * u + self.cd[1][1] * v).to_radians();
        let (sd0, cd0) = self.dec0.sin_cos();
        let denom = cd0 - eta * sd0;
        let ra = (self.ra0 + xi.atan2(denom)).rem_euclid(core::f64::consts::TAU);
        let dec = (sd0 + eta * cd0).atan2(denom.hypot(xi));
        (ra, dec)
    }

    /// (RA, Dec) in radians → 1-based FITS pixel, or `None` for a point on the far
    /// side of the tangent plane or a singular CD matrix.
    #[must_use]
    pub fn sky_to_pixel(&self, ra: f64, dec: f64) -> Option<(f64, f64)> {
        let (u, v) = self.sky_to_linear(ra, dec)?;
        let (u, v) = self
            .sip
            .as_ref()
            .map_or((u, v), |s| s.linear_to_pixel(u, v));
        Some((u + self.crpix1, v + self.crpix2))
    }

    /// (RA, Dec) → `CD⁻¹ · (ξ, η)`: the pixel offsets from CRPIX the linear WCS alone
    /// puts it at.
    fn sky_to_linear(&self, ra: f64, dec: f64) -> Option<(f64, f64)> {
        let (sd0, cd0) = self.dec0.sin_cos();
        let (sd, cd) = dec.sin_cos();
        let (sdr, cdr) = (ra - self.ra0).sin_cos();
        let h = sd * sd0 + cd * cd0 * cdr;
        if h <= 0.0 {
            return None;
        }
        let xi = (cd * sdr / h).to_degrees();
        let eta = ((sd * cd0 - cd * sd0 * cdr) / h).to_degrees();
        let [[a, b], [c, d]] = self.cd;
        let det = a * d - b * c;
        if det == 0.0 || !det.is_finite() {
            return None;
        }
        Some(((d * xi - b * eta) / det, (a * eta - c * xi) / det))
    }
}

/// Fit third-order SIP polynomials to the solution's verified star pairs.
///
/// `width` and `height` are the image's (unbinned) dimensions. Returns `None`, and
/// the solution should stay linear, when:
///
/// - there are fewer than [`MIN_SIP_STARS`] pairs, or they leave part of the frame
///   empty (see `covers_the_frame`);
/// - the fit is singular, or would move an image corner by more than 5% of the
///   half-diagonal (the sign of a cubic extrapolating beyond its stars);
/// - the distortion is not significant: the cubic's 14 extra terms must explain
///   the residuals of a linear fit better than chance would, by an F-test at
///   [`MIN_SIP_F`]. On an undistorted field a cubic only fits the centroid noise,
///   and that noise is largest where it matters most, in the corners: over the
///   benchmark corpus (survey images, reprojected, so distortion-free) fitting
///   one regardless made the worst-corner error 20-60% worse, as it does for
///   `astap_cli -sip`.
///
/// One round of 3σ clipping drops pairs the fit cannot explain before the final fit.
#[must_use]
pub fn fit_sip(wcs: &WcsSolution, width: usize, height: usize) -> Option<Sip> {
    let tan = TanWcs {
        sip: None,
        ..TanWcs::from(wcs)
    };
    // (measured offset from CRPIX, where the linear WCS puts the catalogue star)
    let mut pairs: Vec<Pair> = wcs
        .matched_stars
        .iter()
        .filter_map(|m| {
            let linear = tan.sky_to_linear(m.ra, m.dec)?;
            Some(((m.x - wcs.crpix1, m.y - wcs.crpix2), linear))
        })
        .collect();
    if pairs.len() < MIN_SIP_STARS {
        log::info!(
            "Not enough stars for calculating SIP: {} of {MIN_SIP_STARS}.",
            pairs.len()
        );
        return None;
    }
    if !covers_the_frame(&wcs.matched_stars, width, height) {
        log::info!("Not calculating SIP: the matched stars do not cover the whole image.");
        return None;
    }

    // Fit in coordinates scaled to about ±1, so the cubic terms are not 10¹⁰ times
    // the constant one; the coefficients are scaled back afterwards.
    let scale = 0.5 * (width.max(height).max(2) as f64);
    let (mut a, mut b) = fit_pair(&pairs, scale, SIP_TERMS.len())?;

    let residual = |&((u, v), (lu, lv)): &Pair, a: &[f64; 10], b: &[f64; 10]| {
        (u + poly(a, u, v) - lu).hypot(v + poly(b, u, v) - lv)
    };
    let rms = (pairs
        .iter()
        .map(|p| residual(p, &a, &b).powi(2))
        .sum::<f64>()
        / pairs.len() as f64)
        .sqrt();
    let before = pairs.len();
    pairs.retain(|p| residual(p, &a, &b) <= 3.0 * rms);
    if pairs.len() < before && pairs.len() >= MIN_SIP_STARS {
        (a, b) = fit_pair(&pairs, scale, SIP_TERMS.len())?;
    }

    // Is the distortion real? Compare with the best constant-plus-linear fit.
    let rss =
        |a: &[f64; 10], b: &[f64; 10]| pairs.iter().map(|p| residual(p, a, b).powi(2)).sum::<f64>();
    let rss_cubic = rss(&a, &b);
    let (la, lb) = fit_pair(&pairs, scale, 3)?;
    let rss_linear = rss(&la, &lb);
    let extra = 2.0 * (SIP_TERMS.len() - 3) as f64;
    let dof = 2.0 * pairs.len().saturating_sub(SIP_TERMS.len()).max(1) as f64;
    let f = ((rss_linear - rss_cubic) / extra) / (rss_cubic / dof);
    if !(f >= MIN_SIP_F || rss_cubic <= 0.0 && rss_linear > 0.0) {
        log::info!(
            "No significant distortion (F = {f:.1} over a linear fit, need {MIN_SIP_F}); \
             leaving the solution linear."
        );
        return None;
    }

    // The inverse, sky → pixel, from the same pairs with the roles swapped.
    let swapped: Vec<_> = pairs.iter().map(|&(m, l)| (l, m)).collect();
    let (ap, bp) = fit_pair(&swapped, scale, SIP_TERMS.len())?;
    let sip = Sip { a, b, ap, bp };

    let (hw, hh) = (0.5 * width as f64, 0.5 * height as f64);
    let limit = MAX_CORNER_CORRECTION * hw.hypot(hh);
    for (cu, cv) in [(-hw, -hh), (hw, -hh), (-hw, hh), (hw, hh)] {
        let (lu, lv) = sip.pixel_to_linear(cu, cv);
        let (pu, pv) = sip.linear_to_pixel(cu, cv);
        let worst = (lu - cu).hypot(lv - cv).max((pu - cu).hypot(pv - cv));
        if worst.is_nan() || worst > limit {
            log::info!("SIP fit rejected: it moves an image corner by {worst:.1} px.");
            return None;
        }
    }
    log::info!(
        "SIP distortion fitted to {} stars: residual {:.2} px, {:.2} px without it (F = {f:.1}).",
        pairs.len(),
        (rss_cubic / pairs.len() as f64).sqrt(),
        (rss_linear / pairs.len() as f64).sqrt(),
    );
    Some(sip)
}

/// Whether the matched stars reach every part of the frame: every cell of a 3×3
/// grid holds at least two of them.
///
/// A cubic is pinned down only where it has stars; beyond them it extrapolates,
/// and a fit to stars covering half the frame can move the empty half's corners by
/// tens of pixels while fitting its own stars perfectly. ASTAP needs no such test
/// because its catalogue read always covers the frame; arcsec's verified pairs can
/// fall short of it (a field that straddles a catalogue tile boundary is matched
/// only on one side), and a linear solution is better than a wrong cubic.
fn covers_the_frame(stars: &[crate::types::MatchedStar], width: usize, height: usize) -> bool {
    let mut cells = [0usize; COVERAGE_GRID * COVERAGE_GRID];
    let cell = |v: f64, size: usize| {
        let f = (v - 0.5) / size.max(1) as f64;
        ((f * COVERAGE_GRID as f64).floor().max(0.0) as usize).min(COVERAGE_GRID - 1)
    };
    for m in stars {
        cells[cell(m.y, height) * COVERAGE_GRID + cell(m.x, width)] += 1;
    }
    cells.iter().all(|&n| n >= COVERAGE_MIN)
}

/// Least-squares polynomial `target - source = poly(source)` along each axis, over
/// the first `terms` of [`SIP_TERMS`] (3: constant and linear; 10: the full cubic),
/// with the coordinates divided by `scale` during the fit.
fn fit_pair(pairs: &[Pair], scale: f64, terms: usize) -> Option<([f64; 10], [f64; 10])> {
    let columns: Vec<Vec<f64>> = SIP_TERMS[..terms]
        .iter()
        .map(|&(p, q)| {
            pairs
                .iter()
                .map(|&((u, v), _)| (u / scale).powi(p) * (v / scale).powi(q))
                .collect()
        })
        .collect();
    let du: Vec<f64> = pairs.iter().map(|&((u, _), (tu, _))| tu - u).collect();
    let dv: Vec<f64> = pairs.iter().map(|&((_, v), (_, tv))| tv - v).collect();
    let cu = lsq_fit(&columns, &du).ok()?;
    let cv = lsq_fit(&columns, &dv).ok()?;

    let mut a = [0.0; 10];
    let mut b = [0.0; 10];
    for (k, &(p, q)) in SIP_TERMS[..terms].iter().enumerate() {
        let s = scale.powi(p + q);
        a[k] = cu[k] / s;
        b[k] = cv[k] / s;
    }
    (a.iter().chain(&b).all(|c| c.is_finite())).then_some((a, b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::Rng;
    use crate::types::{MatchedStar, PlateConstants};

    const W: usize = 3000;
    const H: usize = 2000;

    /// A 1.2"/px solution at (150°, +40°), rotated slightly, centred on the frame.
    fn linear() -> WcsSolution {
        let s = 1.2 / 3600.0;
        let rot = 0.3_f64.to_radians();
        let (sr, cr) = rot.sin_cos();
        let mut w = crate::wcs::derive_wcs(
            150f64.to_radians(),
            40f64.to_radians(),
            &PlateConstants {
                a: 1.0,
                b: 0.0,
                c: 0.0,
                d: 0.0,
                e: 1.0,
                f: 0.0,
            },
            W,
            H,
        );
        w.cd1_1 = -s * cr;
        w.cd1_2 = s * sr;
        w.cd2_1 = s * sr;
        w.cd2_2 = s * cr;
        w.ra0 = 150f64.to_radians();
        w.dec0 = 40f64.to_radians();
        w
    }

    /// A barrel-like distortion: a radial cubic plus a small quadratic tilt, of
    /// about 6 px at the corners.
    fn truth_sip() -> Sip {
        let k = 6.0 / 1800f64.powi(3);
        let mut a = [0.0; 10];
        let mut b = [0.0; 10];
        // u (u² + v²) and v (u² + v²)
        a[6] = k; // u³
        a[8] = k; // u v²
        b[7] = k; // u² v
        b[9] = k; // v³
        a[3] = 2e-7; // u²
        b[4] = -1e-7; // u v
        Sip {
            a,
            b,
            ap: [0.0; 10],
            bp: [0.0; 10],
        }
    }

    /// Stars scattered over the frame, their catalogue positions taken through the
    /// distorted WCS, with `noise` px of measuring error.
    fn matches(truth: &TanWcs, n: usize, noise: f64, seed: u64) -> Vec<MatchedStar> {
        let mut rng = Rng::new(seed);
        (0..n)
            .map(|_| {
                let (x, y) = (rng.range(1.0, W as f64), rng.range(1.0, H as f64));
                let (ra, dec) = truth.pixel_to_sky(x, y);
                MatchedStar {
                    x: x + noise * rng.gauss(),
                    y: y + noise * rng.gauss(),
                    ra,
                    dec,
                }
            })
            .collect()
    }

    #[test]
    fn recovers_a_known_distortion() {
        let mut wcs = linear();
        let truth = TanWcs {
            sip: Some(truth_sip()),
            ..TanWcs::from(&wcs)
        };
        wcs.matched_stars = matches(&truth, 300, 0.05, 1);
        let sip = fit_sip(&wcs, W, H).expect("the fit must succeed");
        wcs.sip = Some(sip);
        let fitted = TanWcs::from(&wcs);
        let linear_only = TanWcs {
            sip: None,
            ..fitted.clone()
        };

        // Across the frame, corners included, the fitted WCS agrees with the truth
        // to a small fraction of a pixel, where the linear one is out by pixels.
        let arcsec_px = 1.2;
        let mut worst_sip: f64 = 0.0;
        let mut worst_lin: f64 = 0.0;
        for &(x, y) in &[
            (1.0, 1.0),
            (W as f64, 1.0),
            (1.0, H as f64),
            (W as f64, H as f64),
            (1500.0, 1000.0),
            (700.0, 1600.0),
        ] {
            let (tr, td) = truth.pixel_to_sky(x, y);
            let sep = |(r, d): (f64, f64)| {
                crate::math::coords::ang_sep(r, d, tr, td).to_degrees() * 3600.0 / arcsec_px
            };
            worst_sip = worst_sip.max(sep(fitted.pixel_to_sky(x, y)));
            worst_lin = worst_lin.max(sep(linear_only.pixel_to_sky(x, y)));
        }
        assert!(worst_sip < 0.05, "SIP error {worst_sip} px");
        assert!(worst_lin > 3.0, "linear error {worst_lin} px");

        // The inverse polynomial takes the sky back to the pixel it came from.
        for &(x, y) in &[(1.0, 1.0), (W as f64, H as f64), (400.0, 1500.0)] {
            let (ra, dec) = fitted.pixel_to_sky(x, y);
            let (bx, by) = fitted.sky_to_pixel(ra, dec).unwrap();
            assert!((bx - x).hypot(by - y) < 0.05, "({x},{y}) -> ({bx},{by})");
        }
    }

    #[test]
    fn an_undistorted_field_stays_linear() {
        let mut wcs = linear();
        let truth = TanWcs::from(&wcs);
        for seed in 0..20 {
            wcs.matched_stars = matches(&truth, 200, 0.3, 100 + seed);
            assert!(fit_sip(&wcs, W, H).is_none(), "seed {seed}");
        }
    }

    #[test]
    fn a_small_distortion_under_noise_is_still_found() {
        // 2 px at the corners, under 0.3 px of centroid noise.
        let mut small = truth_sip();
        for c in small.a.iter_mut().chain(small.b.iter_mut()) {
            *c /= 3.0;
        }
        let mut wcs = linear();
        let truth = TanWcs {
            sip: Some(small),
            ..TanWcs::from(&wcs)
        };
        wcs.matched_stars = matches(&truth, 200, 0.3, 7);
        assert!(fit_sip(&wcs, W, H).is_some());
    }

    #[test]
    fn outliers_are_clipped() {
        let mut wcs = linear();
        let truth = TanWcs {
            sip: Some(truth_sip()),
            ..TanWcs::from(&wcs)
        };
        wcs.matched_stars = matches(&truth, 300, 0.05, 3);
        // Five mismatches, 4 px off.
        for m in wcs.matched_stars.iter_mut().take(5) {
            m.x += 4.0;
        }
        wcs.sip = fit_sip(&wcs, W, H);
        let fitted = TanWcs::from(&wcs);
        let (tr, td) = truth.pixel_to_sky(W as f64, H as f64);
        let (r, d) = fitted.pixel_to_sky(W as f64, H as f64);
        let err = crate::math::coords::ang_sep(r, d, tr, td).to_degrees() * 3600.0 / 1.2;
        assert!(err < 0.1, "corner error {err} px");
    }

    #[test]
    fn too_few_or_clustered_stars_give_no_fit() {
        let mut wcs = linear();
        let truth = TanWcs::from(&wcs);
        wcs.matched_stars = matches(&truth, MIN_SIP_STARS - 1, 0.1, 4);
        assert!(fit_sip(&wcs, W, H).is_none());

        // Plenty of stars, but only on the left half of the frame.
        wcs.matched_stars = matches(&truth, 300, 0.05, 6);
        wcs.matched_stars.retain(|m| m.x < 0.5 * W as f64);
        assert!(wcs.matched_stars.len() > 100);
        assert!(fit_sip(&wcs, W, H).is_none());

        // Stars only in a small central patch, with a lot of noise: the cubic is
        // unconstrained at the corners and the fit must be refused.
        let mut rng = Rng::new(5);
        wcs.matched_stars = (0..40)
            .map(|_| {
                let (x, y) = (rng.range(1450.0, 1550.0), rng.range(950.0, 1050.0));
                let (ra, dec) = truth.pixel_to_sky(x, y);
                MatchedStar {
                    x: x + 2.0 * rng.gauss(),
                    y: y + 2.0 * rng.gauss(),
                    ra,
                    dec,
                }
            })
            .collect();
        assert!(fit_sip(&wcs, W, H).is_none());
    }

    #[test]
    fn tan_round_trips_without_sip() {
        let w = TanWcs::from(&linear());
        let (ra, dec) = w.pixel_to_sky(w.crpix1, w.crpix2);
        assert!((ra - 150f64.to_radians()).abs() < 1e-12);
        assert!((dec - 40f64.to_radians()).abs() < 1e-12);
        let (x, y) = w.sky_to_pixel(ra + 0.001, dec - 0.002).unwrap();
        let (r2, d2) = w.pixel_to_sky(x, y);
        assert!((r2 - ra - 0.001).abs() < 1e-12 && (d2 - dec + 0.002).abs() < 1e-12);
        // The antipode is behind the tangent plane.
        assert!(w.sky_to_pixel(ra + core::f64::consts::PI, -dec).is_none());
    }
}
