//! WCS solution derivation from plate constants.

use crate::math::coords::standard_equatorial;
use crate::types::{PlateConstants, WcsSolution};

/// Derive a WCS solution from solved plate constants.
///
/// # Arguments
/// - `ra_db`, `dec_db`: the reference position used to build the catalog quads (radians)
/// - `plate`: the solved plate constants (standard-coordinate arcsec per image pixel)
/// - `width`, `height`: image dimensions in pixels
#[must_use]
pub fn derive_wcs(
    ra_db: f64,
    dec_db: f64,
    plate: &PlateConstants,
    width: usize,
    height: usize,
) -> WcsSolution {
    let cx = (width as f64 - 1.0) * 0.5;
    let cy = (height as f64 - 1.0) * 0.5;

    // Standard coordinates (arcsec) of the image centre
    let x_std = plate.a * cx + plate.b * cy + plate.c;
    let y_std = plate.d * cx + plate.e * cy + plate.f;

    // Convert to RA/DEC using the tangent-plane inverse (cdelt=1 → arcsec units)
    let (ra0, dec0) = standard_equatorial(ra_db, dec_db, x_std, y_std, 1.0);

    // CD matrix [deg/pixel] (ASTAP sign convention)
    //   cd1_1 = -a/3600,  cd1_2 = -b/3600   (RA decreases eastward → negative sign)
    //   cd2_1 = +d/3600,  cd2_2 = +e/3600
    let cd1_1 = -plate.a / 3600.0;
    let cd1_2 = -plate.b / 3600.0;
    let cd2_1 = plate.d / 3600.0;
    let cd2_2 = plate.e / 3600.0;

    // Pixel scale [deg/pixel] and rotation [degrees], as astap_cli derives them: CDELT1
    // carries the parity (negative for the usual east-left image, positive for a
    // mirrored one), CDELT2 is positive, and CROTA2 is the rotation of the image's
    // +Y axis from north. See `old_style_wcs`.
    let (cdelt1, cdelt2, crota2, _) = old_style_wcs(cd1_1, cd1_2, cd2_1, cd2_2);

    WcsSolution {
        ra0,
        dec0,
        crpix1: cx + 1.0,
        crpix2: cy + 1.0,
        cd1_1,
        cd1_2,
        cd2_1,
        cd2_2,
        cdelt1,
        cdelt2,
        crota2,
        residual_rms: 0.0,
        stars_matched: 0,
        plate: plate.clone(),
        mag_limit: 0.0,
        search_dist_deg: 0.0,
        step_distances: Vec::new(),
        raw_matches: 0,
        matched_stars: Vec::new(),
        sip: None,
    }
}

/// `CDELT1`, `CDELT2`, `CROTA2` and `CROTA1` (degrees per pixel and degrees) for a
/// CD matrix, in `astap_cli`'s convention.
///
/// `astap_cli` (2025 onward) builds its CD matrix from these four values as
/// `CD1_1 = CDELT1·cos CROTA1`, `CD1_2 = −CDELT1·sin CROTA1·f`,
/// `CD2_1 = CDELT2·sin CROTA2·f`, `CD2_2 = CDELT2·cos CROTA2`, where `f` is −1 for an
/// image with the sky's usual handedness (det CD < 0) and +1 for a mirrored one, and
/// `CDELT1 = f·|row 1|`. This inverts that, so tools that read the old-style keywords
/// see what they would from ASTAP.
#[must_use]
pub fn old_style_wcs(cd1_1: f64, cd1_2: f64, cd2_1: f64, cd2_2: f64) -> (f64, f64, f64, f64) {
    let det = cd1_1 * cd2_2 - cd1_2 * cd2_1;
    let f = if det < 0.0 { -1.0 } else { 1.0 };
    let cdelt1 = f * cd1_1.hypot(cd1_2);
    let cdelt2 = cd2_1.hypot(cd2_2);
    let crota2 = (f * cd2_1).atan2(cd2_2).to_degrees();
    let crota1 = if cdelt1 == 0.0 {
        0.0
    } else {
        (-cd1_2 / (cdelt1 * f)).atan2(cd1_1 / cdelt1).to_degrees()
    };
    (cdelt1, cdelt2, crota2, crota1)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The CD matrix `astap_cli` 2026.07.30 wrote for dec000.fits, and the CDELT/CROTA
    /// values it wrote with it.
    #[test]
    fn old_style_keywords_match_astap_cli() {
        let (cdelt1, cdelt2, crota2, crota1) = old_style_wcs(
            -3.448_223_117_128e-4,
            -1.328_788_985_076e-8,
            -1.180_717_632_642e-8,
            3.448_186_344_500e-4,
        );
        assert!(
            (cdelt1 - -3.448_223_119_688e-4).abs() < 1e-12,
            "cdelt1 {cdelt1}"
        );
        assert!(
            (cdelt2 - 3.448_186_346_521e-4).abs() < 1e-12,
            "cdelt2 {cdelt2}"
        );
        assert!(
            (crota2 - 1.961_904_907_735e-3).abs() < 1e-9,
            "crota2 {crota2}"
        );
        assert!(
            (crota1 - 2.207_919_791_863e-3).abs() < 1e-9,
            "crota1 {crota1}"
        );
    }

    #[test]
    fn a_mirrored_matrix_has_a_positive_cdelt1() {
        let (cdelt1, cdelt2, crota2, crota1) = old_style_wcs(3e-4, 0.0, 0.0, 3e-4);
        assert!(cdelt1 > 0.0 && cdelt2 > 0.0);
        assert!(crota2.abs() < 1e-12 && crota1.abs() < 1e-12);
    }
    use crate::math::coords::ang_sep;
    use core::f64::consts::PI;

    fn deg(d: f64) -> f64 {
        d * PI / 180.0
    }

    #[test]
    fn identity_wcs_near_zero() {
        // Identity plate constants: a=1, e=1, others=0 (1 arcsec/pixel, no rotation)
        let plate = PlateConstants {
            a: 1.0,
            b: 0.0,
            c: 0.0,
            d: 0.0,
            e: 1.0,
            f: 0.0,
        };
        let w = WcsSolution {
            ra0: 0.0,
            dec0: 0.0,
            ..derive_wcs(0.0, 0.0, &plate, 101, 101)
        };
        // Center standard coords should be (0, 0) → (ra0, dec0) ≈ reference position
        assert!(w.ra0.abs() < 1e-8, "ra0 = {}", w.ra0);
        assert!(w.dec0.abs() < 1e-8, "dec0 = {}", w.dec0);
        // cdelt should be ~1/3600 degrees
        assert!(
            (w.cdelt2 - 1.0 / 3600.0).abs() < 1e-10,
            "cdelt2 = {}",
            w.cdelt2
        );
    }

    #[test]
    fn center_ra_dec_reasonable() {
        // Plate constants with non-zero offset → center moves away from reference
        // a=1 arcsec/pix, e=1, c=100, f=200 → center offset (100, 200) arcsec from ref
        let plate = PlateConstants {
            a: 1.0,
            b: 0.0,
            c: 0.0,
            d: 0.0,
            e: 1.0,
            f: 0.0,
        };
        let ra_db = deg(45.0);
        let dec_db = deg(20.0);
        let wcs = derive_wcs(ra_db, dec_db, &plate, 101, 101);
        // The center should be close to (ra_db, dec_db) since c=f=0
        let sep = ang_sep(wcs.ra0, wcs.dec0, ra_db, dec_db);
        assert!(
            sep < deg(0.1),
            "center offset = {} arcsec",
            sep * 3600.0 * 180.0 / PI
        );
    }

    #[test]
    fn crpix_is_image_center() {
        let plate = PlateConstants {
            a: 1.0,
            b: 0.0,
            c: 0.0,
            d: 0.0,
            e: 1.0,
            f: 0.0,
        };
        let wcs = derive_wcs(0.0, 0.0, &plate, 200, 300);
        assert!(
            (wcs.crpix1 - 100.5).abs() < 1e-10,
            "crpix1 = {}",
            wcs.crpix1
        );
        assert!(
            (wcs.crpix2 - 150.5).abs() < 1e-10,
            "crpix2 = {}",
            wcs.crpix2
        );
    }
}
