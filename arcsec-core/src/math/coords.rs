//! Coordinate transforms between equatorial (RA/Dec) and tangent-plane (standard) coordinates.
//!
//! Gnomonic (TAN) projection, per Calabretta & Greisen (2002), "Representations of
//! celestial coordinates in FITS", A&A 395, 1077, section 5.1.

use core::f64::consts::PI;

/// Factor: cdelt [arcsec/px] → radians/px.
/// `cdelt / (3600 * 180/π)` = `cdelt * π / (3600 * 180)`
#[inline]
fn arcsec_to_rad(cdelt: f64) -> f64 {
    cdelt * PI / (3600.0 * 180.0)
}

/// Convert equatorial (ra, dec) to tangent-plane (x, y) in pixels.
///
/// - `ra0`, `dec0`: field center (radians)
/// - `ra`, `dec`: star position (radians)
/// - `cdelt`: pixel scale in arcsec/pixel (use 1.0 for unitless standard coords)
///
/// Returns `(x, y)` in pixels from center (or standard coords when `cdelt = 1.0`).
///
#[must_use]
pub fn equatorial_standard(ra0: f64, dec0: f64, ra: f64, dec: f64, cdelt: f64) -> (f64, f64) {
    let (sin_dec0, cos_dec0) = dec0.sin_cos();
    let (sin_dec, cos_dec) = dec.sin_cos();
    let (sin_dra, cos_dra) = (ra - ra0).sin_cos();

    let scale = arcsec_to_rad(cdelt);
    // dv = projection_factor * cdelt/(3600*180/π) = projection * scale
    let dv = (cos_dec0 * cos_dec * cos_dra + sin_dec0 * sin_dec) * scale;

    let xx = -cos_dec * sin_dra / dv;
    let yy = -(sin_dec0 * cos_dec * cos_dra - cos_dec0 * sin_dec) / dv;
    (xx, yy)
}

/// Convert tangent-plane (x, y) in pixels to equatorial (ra, dec) in radians.
///
/// - `ra0`, `dec0`: field center (radians)
/// - `x`, `y`: pixel offset from center
/// - `cdelt`: pixel scale in arcsec/pixel
///
/// Returns `(ra, dec)` in radians, with RA normalised to [0, 2π).
///
#[must_use]
pub fn standard_equatorial(ra0: f64, dec0: f64, x: f64, y: f64, cdelt: f64) -> (f64, f64) {
    let (sin_dec0, cos_dec0) = dec0.sin_cos();
    let scale = arcsec_to_rad(cdelt);
    let xs = x * scale;
    let ys = y * scale;

    let mut ra = ra0 + (-xs).atan2(cos_dec0 - ys * sin_dec0);
    if ra >= 2.0 * PI {
        ra -= 2.0 * PI;
    }
    if ra < 0.0 {
        ra += 2.0 * PI;
    }

    let dec = ((sin_dec0 + ys * cos_dec0) / (1.0 + xs * xs + ys * ys).sqrt()).asin();
    (ra, dec)
}

/// Angular separation between two positions (radians in, radians out).
#[must_use]
pub fn ang_sep(ra1: f64, dec1: f64, ra2: f64, dec2: f64) -> f64 {
    let (sin_d1, cos_d1) = dec1.sin_cos();
    let (sin_d2, cos_d2) = dec2.sin_cos();
    let cos_sep = (sin_d1 * sin_d2 + cos_d1 * cos_d2 * (ra1 - ra2).cos()).clamp(-1.0, 1.0);
    cos_sep.acos()
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;

    const DEG: f64 = PI / 180.0;
    const ARCSEC: f64 = DEG / 3600.0;
    // 1e-9 radians ≈ 0.2 milliarcseconds: tight but achievable for small offsets
    const TOL: f64 = 1e-9;

    /// Round-trip: `equatorial_standard` → `standard_equatorial` recovers original position.
    /// Offsets kept small (< 20 arcmin at 2.5 arcsec/px = 480 px) for numerical stability.
    #[test]
    fn round_trip() {
        let ra0 = 101.29 * DEG;
        let dec0 = -16.72 * DEG;
        let cdelt = 2.5; // arcsec/pixel

        let test_cases = [
            (101.30 * DEG, -16.71 * DEG), // ~80 arcsec offset
            (101.00 * DEG, -17.00 * DEG), // ~23 arcmin offset
            (101.29 * DEG, -16.72 * DEG), // center itself
        ];

        for (ra, dec) in test_cases {
            let (x, y) = equatorial_standard(ra0, dec0, ra, dec, cdelt);
            let (ra2, dec2) = standard_equatorial(ra0, dec0, x, y, cdelt);
            // For RA, account for possible 2π wrap
            let dra = (ra - ra2).rem_euclid(2.0 * PI);
            let dra = dra.min(2.0 * PI - dra);
            assert!(
                dra < TOL,
                "RA round-trip failed: {ra:.9} vs {ra2:.9} (diff {dra:.2e})"
            );
            assert!(
                (dec - dec2).abs() < TOL,
                "Dec round-trip failed: {dec:.9} vs {dec2:.9}"
            );
        }
    }

    /// RA wrap-around: star just below 360° relative to field center near 0°.
    /// Angular separation is ~0.03° = 108 arcsec across the 0° boundary.
    #[test]
    fn ra_wraparound() {
        let ra0 = 0.02 * DEG;
        let dec0 = 10.0 * DEG;
        let cdelt = 1.0; // 1 arcsec/px

        let ra = 359.99 * DEG; // 0.03° short of 360°, so 0.01° from ra0
        let dec = 10.0 * DEG;
        let (x, y) = equatorial_standard(ra0, dec0, ra, dec, cdelt);
        let (ra2, dec2) = standard_equatorial(ra0, dec0, x, y, cdelt);

        // Near-0° boundary: trig accumulates ~1e-8 rad error, still fine for plate solving
        let sep = ang_sep(ra, dec, ra2, dec2);
        assert!(sep < 1e-7, "RA wraparound sep = {sep:.2e} rad");
        assert!(
            (dec - dec2).abs() < TOL,
            "dec mismatch {}",
            (dec - dec2).abs()
        );
    }

    /// Near-pole round-trip: small offset (10 arcsec) near dec = 89°.
    #[test]
    fn near_pole() {
        let ra0 = 90.0 * DEG;
        let dec0 = 89.0 * DEG;
        let cdelt = 1.0;

        // Small 10 arcsec offset in dec to stay numerically safe
        let ra = 90.0 * DEG;
        let dec = 89.0 * DEG + 10.0 * ARCSEC;
        let (x, y) = equatorial_standard(ra0, dec0, ra, dec, cdelt);
        let (ra2, dec2) = standard_equatorial(ra0, dec0, x, y, cdelt);
        let sep = ang_sep(ra, dec, ra2, dec2);
        assert!(sep < TOL, "near-pole sep = {sep:.2e} rad");
        let _ = (x, y, ra2, dec2);
    }

    /// `ang_sep` agrees with direct formula for known pairs.
    #[test]
    fn ang_sep_known() {
        // Same point → 0
        assert!(ang_sep(1.0, 0.5, 1.0, 0.5) < 1e-15);

        // Antipodal points → π
        let sep = ang_sep(0.0, PI / 2.0, PI, -PI / 2.0);
        assert!((sep - PI).abs() < 1e-12, "antipodal sep = {sep}");

        // NCP to a point 10° away in dec
        let sep = ang_sep(0.0, PI / 2.0, 0.0, PI / 2.0 - 10.0 * DEG);
        assert!((sep - 10.0 * DEG).abs() < 1e-12, "10deg from pole: {sep}");
    }

    /// center → exactly zero offset
    #[test]
    fn center_is_zero() {
        let ra0 = 45.0 * DEG;
        let dec0 = 30.0 * DEG;
        let (x, y) = equatorial_standard(ra0, dec0, ra0, dec0, 1.0);
        assert!(x.abs() < 1e-12, "x at center = {x}");
        assert!(y.abs() < 1e-12, "y at center = {y}");
    }

    /// 10 arcsec offset in dec → ~10 pixel offset at cdelt=1 arcsec/px
    #[test]
    fn small_dec_offset() {
        let ra0 = 100.0 * DEG;
        let dec0 = 30.0 * DEG;
        let cdelt = 1.0;

        let ra = ra0;
        let dec = dec0 + 10.0 * ARCSEC;
        let (x, y) = equatorial_standard(ra0, dec0, ra, dec, cdelt);
        // Tangent-plane handedness: +Dec gives +y, and +RA (east) gives -x.
        // derive_wcs depends on this when it negates cd1_1/cd1_2 and leaves
        // cd2_1/cd2_2 positive.
        assert!(x.abs() < 0.001, "x should be ~0, got {x}");
        assert!((y - 10.0).abs() < 0.001, "y should be ~+10, got {y}");

        // The other axis of the same convention, previously untested.
        let (x, y) = equatorial_standard(ra0, dec0, ra0 + 10.0 * ARCSEC, dec0, cdelt);
        assert!(x < 0.0, "+RA (east) should give -x, got {x}");
        assert!(y.abs() < 0.001, "y should be ~0, got {y}");
    }

    /// Against the textbook gnomonic projection in `test_support` (xi east, eta
    /// north): `equatorial_standard` is the same projection with x = -xi, for
    /// tangent points everywhere including both poles and across RA 0, and the
    /// round trip closes to micro-arcseconds.
    #[test]
    fn matches_an_independent_projection_everywhere() {
        use crate::test_support::{Rng, gnomonic, separation};
        let mut rng = Rng::new(17);
        let rad_per_arcsec = PI / (180.0 * 3600.0);
        let mut n = 0;
        for k in 0..3000 {
            // Tangent points: uniform, plus both poles exactly and RA ≈ 0/2π.
            let (ra0, dec0) = match k {
                0 => (1.0, PI / 2.0),
                1 => (4.0, -PI / 2.0),
                2 => (2.0 * PI - 1e-12, 0.3),
                _ => (rng.range(0.0, 2.0 * PI), rng.range(-1.0, 1.0).asin()),
            };
            // A star up to 10° away in a random direction.
            let (xi, eta) = (rng.range(-0.17, 0.17), rng.range(-0.17, 0.17));
            let (ra, dec) = crate::test_support::inverse_gnomonic(ra0, dec0, xi, eta);

            let (x, y) = equatorial_standard(ra0, dec0, ra, dec, 1.0);
            let (xi2, eta2) = gnomonic(ra0, dec0, ra, dec).unwrap();
            assert!(
                (x * rad_per_arcsec + xi2).abs() < 1e-12
                    && (y * rad_per_arcsec - eta2).abs() < 1e-12,
                "tangent ({ra0}, {dec0}): ({x}, {y}) vs ({xi2}, {eta2})"
            );
            let (ra3, dec3) = standard_equatorial(ra0, dec0, x, y, 1.0);
            assert!((0.0..2.0 * PI).contains(&ra3));
            assert!(
                separation(ra, dec, ra3, dec3) < 1e-11,
                "round trip at ({ra0}, {dec0})"
            );
            n += 1;
        }
        assert_eq!(n, 3000);
    }

    #[test]
    fn ang_sep_is_a_metric_on_the_sphere() {
        use crate::test_support::{Rng, separation};
        let mut rng = Rng::new(18);
        for _ in 0..1000 {
            let p: Vec<(f64, f64)> = (0..3)
                .map(|_| (rng.range(0.0, 2.0 * PI), rng.range(-1.0, 1.0).asin()))
                .collect();
            let d = |i: usize, j: usize| ang_sep(p[i].0, p[i].1, p[j].0, p[j].1);
            assert!((d(0, 1) - d(1, 0)).abs() < 1e-12);
            assert!(d(0, 2) <= d(0, 1) + d(1, 2) + 1e-12);
            assert!((0.0..=PI).contains(&d(0, 1)));
            // Agrees with the haversine form away from the tiny-angle regime.
            assert!((d(0, 1) - separation(p[0].0, p[0].1, p[1].0, p[1].1)).abs() < 1e-7);
        }
        assert!((ang_sep(0.0, PI / 2.0, 3.0, -PI / 2.0) - PI).abs() < 1e-12);
    }
}
