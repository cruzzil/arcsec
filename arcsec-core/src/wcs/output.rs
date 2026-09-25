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

    // Pixel scale [deg/pixel] and rotation [degrees]
    let cdelt1 = -(cd1_1 * cd1_1 + cd1_2 * cd1_2).sqrt(); // FITS: cdelt1 < 0 for east-to-right
    let cdelt2 = (cd2_1 * cd2_1 + cd2_2 * cd2_2).sqrt();
    let crota2 = cd2_1.atan2(cd2_2).to_degrees();

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
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
