//! Sky area lookup for the .290 catalog format (HNSKY/ASTAP 290-tile grid).
//!
//! The .290 databases (G05 for 3°–20° fields, W08 for 20°–80°) carry exactly the
//! same 110-byte header and packed star records as .1476 — only the sky tiling
//! differs, so `format_1476::read_area_file` reads both unchanged.
//!
//! Unlike the 1476 grid, whose rings are equal-width in declination, the 290 grid is
//! **equal-area**: every cell subtends the same solid angle. With ring RA counts
//!
//! ```text
//! 1, 4, 8, 12, 16, 20, 24, 28, 32, 32, 28, 24, 20, 16, 12, 8, 4, 1   (= 290)
//! ```
//!
//! and the two polar caps counting half a cell each (total weight 289), the ring
//! boundaries satisfy
//!
//! ```text
//! sin(dec_k) = -1 + 2 * (cumulative weight up to k) / 289
//! ```
//!
//! That formula was derived from the shipped G05 files and reproduces the observed
//! per-ring declination ranges to better than 0.001° (the residual is just the
//! granularity of the brightest star in each ring).

use core::f64::consts::PI;

/// RA cell count for each of the 18 declination rings, south to north.
pub const RING_N_290: [usize; 18] = [
    1, 4, 8, 12, 16, 20, 24, 28, 32, 32, 28, 24, 20, 16, 12, 8, 4, 1,
];

/// 1-based area number of the first cell in each ring.
pub const RING_BASE_290: [usize; 18] = [
    1, 2, 6, 14, 26, 42, 62, 86, 114, 146, 178, 206, 230, 250, 266, 278, 286, 290,
];

/// DEC boundaries for the 290 grid (radians, 19 values from -π/2 to +π/2).
/// Index 0 = -90°, index 18 = +90°. Adjacent pair defines one DEC ring.
pub const DEC_BOUNDARIES_290: [f64; 19] = [
    -90.0_f64 * PI / 180.0,
    -85.232244043 * PI / 180.0,
    -75.663487557 * PI / 180.0,
    -65.992866371 * PI / 180.0,
    -56.144973872 * PI / 180.0,
    -46.031630674 * PI / 180.0,
    -35.543077453 * PI / 180.0,
    -24.533481154 * PI / 180.0,
    -12.794405888 * PI / 180.0,
    0.0,
    12.794405888 * PI / 180.0,
    24.533481154 * PI / 180.0,
    35.543077453 * PI / 180.0,
    46.031630674 * PI / 180.0,
    56.144973872 * PI / 180.0,
    65.992866371 * PI / 180.0,
    75.663487557 * PI / 180.0,
    85.232244043 * PI / 180.0,
    90.0_f64 * PI / 180.0,
];

/// Ring index (0..17) containing `dec`.
fn ring_of_dec(dec: f64) -> usize {
    for ring in (0usize..18).rev() {
        if dec >= DEC_BOUNDARIES_290[ring] {
            return ring;
        }
    }
    0
}

/// 1-based area number of the cell containing `(ra, dec)`.
#[must_use]
pub fn area_nr_290(ra: f64, dec: f64) -> usize {
    let ring = ring_of_dec(dec);
    let n_ra = RING_N_290[ring];
    let rot = ra.rem_euclid(2.0 * PI) * n_ra as f64 / (2.0 * PI);
    RING_BASE_290[ring] + (rot.floor() as usize).min(n_ra - 1)
}

/// Filename segment for a 290 area number (e.g. area 1 → `"0101.290"`).
///
/// Format is `{ring:02}{cell:02}.290`, both 1-based — the same shape the 1476 grid
/// uses, which is why the south-pole cell is `0101` in both and can be probed to
/// tell the two layouts apart.
#[must_use]
pub fn filename_290(area_nr: usize) -> String {
    let area = area_nr.clamp(1, 290);
    let mut ring = 17usize;
    for r in 0..18 {
        if area >= RING_BASE_290[r] && area < RING_BASE_290[r] + RING_N_290[r] {
            ring = r;
            break;
        }
    }
    let cell = area - RING_BASE_290[ring] + 1;
    format!("{:02}{:02}.290", ring + 1, cell)
}

/// Enumerate every 290 area overlapping a square field of `fov` radians on a side.
///
/// The 1476 reader samples the field's four corners, which is correct only when the
/// field is no taller than one ring. The .290 databases exist for 3°–80° fields,
/// which span many rings and many RA cells, so this walks the rings the field
/// touches and, within each, the RA cells it touches. A field reaching a pole takes
/// that ring whole.
///
/// Areas are returned in south-to-north, increasing-RA order and are unique.
#[must_use]
pub fn find_areas_290(ra: f64, dec: f64, fov: f64) -> Vec<usize> {
    let half = (fov * 0.5).clamp(0.0, PI);
    let dec_lo = (dec - half).max(-PI / 2.0);
    let dec_hi = (dec + half).min(PI / 2.0);

    let ring_lo = ring_of_dec(dec_lo);
    let ring_hi = ring_of_dec(dec_hi);

    // A field that touches a pole wraps in RA, so no RA restriction applies at all.
    let touches_pole = dec_hi >= DEC_BOUNDARIES_290[18] - 1e-12
        || dec_lo <= DEC_BOUNDARIES_290[0] + 1e-12
        || dec + half >= PI / 2.0
        || dec - half <= -PI / 2.0;

    let mut out = Vec::new();
    for ring in ring_lo..=ring_hi {
        let n_ra = RING_N_290[ring];
        let base = RING_BASE_290[ring];

        if n_ra == 1 || touches_pole {
            for c in 0..n_ra {
                out.push(base + c);
            }
            continue;
        }

        // Half-width in RA at the declination of this ring nearest the field centre,
        // where the field is widest in RA terms.
        let d_near = dec.clamp(DEC_BOUNDARIES_290[ring], DEC_BOUNDARIES_290[ring + 1]);
        let cos_d = d_near.cos();
        if cos_d < 1e-6 {
            for c in 0..n_ra {
                out.push(base + c);
            }
            continue;
        }
        let ra_half = half / cos_d;
        if ra_half >= PI {
            for c in 0..n_ra {
                out.push(base + c);
            }
            continue;
        }

        let step = 2.0 * PI / n_ra as f64;
        let c_lo = ((ra - ra_half).rem_euclid(2.0 * PI) / step).floor() as i64;
        let c_hi = ((ra + ra_half).rem_euclid(2.0 * PI) / step).floor() as i64;
        // Walk forward from c_lo to c_hi the short way round, so RA wrap is handled
        // without special cases.
        let span = (c_hi - c_lo).rem_euclid(n_ra as i64);
        for k in 0..=span {
            let c = (c_lo + k).rem_euclid(n_ra as i64) as usize;
            out.push(base + c);
        }
    }

    out.sort_unstable();
    out.dedup();
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn deg(d: f64) -> f64 {
        d * PI / 180.0
    }

    #[test]
    fn ring_table_totals_290() {
        assert_eq!(RING_N_290.iter().sum::<usize>(), 290);
        for r in 0..18 {
            assert_eq!(
                RING_BASE_290[r] + RING_N_290[r],
                if r == 17 { 291 } else { RING_BASE_290[r + 1] },
                "ring {r} base/count are inconsistent"
            );
        }
    }

    #[test]
    fn dec_boundaries_are_equal_area() {
        // Every cell must subtend the same solid angle, with half-weight poles.
        let weights: Vec<f64> = (0..18)
            .map(|r| {
                if r == 0 || r == 17 {
                    0.5
                } else {
                    RING_N_290[r] as f64
                }
            })
            .collect();
        let total: f64 = weights.iter().sum();
        let mut cum = 0.0;
        for r in 0..18 {
            cum += weights[r];
            let want = -1.0 + 2.0 * cum / total;
            let got = DEC_BOUNDARIES_290[r + 1].sin();
            assert!(
                (want - got).abs() < 1e-9,
                "ring {r} boundary: sin want {want}, got {got}"
            );
        }
    }

    #[test]
    fn poles_and_equator_land_in_the_right_areas() {
        assert_eq!(area_nr_290(0.0, deg(-89.9)), 1);
        assert_eq!(area_nr_290(0.0, deg(89.9)), 290);
        // Ring 9 (index 8) is the last southern ring, 32 cells starting at 114.
        assert_eq!(area_nr_290(0.0, deg(-0.001)), 114);
        // Ring 10 (index 9) starts at the equator, 32 cells from 146.
        assert_eq!(area_nr_290(0.0, 0.0), 146);
    }

    #[test]
    fn filenames_round_trip() {
        for area in 1..=290usize {
            let f = filename_290(area);
            assert!(f.ends_with(".290"));
            let ring: usize = f[0..2].parse().unwrap();
            let cell: usize = f[2..4].parse().unwrap();
            assert!((1..=18).contains(&ring), "area {area} -> {f}");
            assert!(
                cell >= 1 && cell <= RING_N_290[ring - 1],
                "area {area} -> {f}"
            );
            assert_eq!(RING_BASE_290[ring - 1] + cell - 1, area);
        }
        assert_eq!(filename_290(1), "0101.290");
        assert_eq!(filename_290(290), "1801.290");
        assert_eq!(filename_290(2), "0201.290");
        assert_eq!(filename_290(5), "0204.290");
    }

    #[test]
    fn small_field_takes_few_areas_and_contains_its_own_centre() {
        for &(ra, dec) in &[
            (0.0, 0.0),
            (deg(180.0), deg(30.0)),
            (deg(275.0), deg(-45.0)),
            (deg(10.0), deg(70.0)),
        ] {
            let areas = find_areas_290(ra, dec, deg(4.0));
            assert!(
                areas.contains(&area_nr_290(ra, dec)),
                "field at ({ra}, {dec}) does not include its own centre cell"
            );
            assert!(areas.len() <= 8, "4° field took {} areas", areas.len());
        }
    }

    #[test]
    fn wide_field_spans_many_rings() {
        // A 60° field must reach well beyond one ring — this is the case the 1476
        // four-corner sampling cannot express.
        let areas = find_areas_290(deg(180.0), 0.0, deg(60.0));
        assert!(
            areas.len() > 20,
            "60° field took only {} areas",
            areas.len()
        );
        assert!(areas.contains(&area_nr_290(deg(180.0), 0.0)));
    }

    #[test]
    fn polar_field_takes_whole_rings() {
        let areas = find_areas_290(0.0, deg(88.0), deg(20.0));
        assert!(areas.contains(&290), "must include the north polar cap");
        // The ring below the cap has 4 cells; all of them are in RA range near a pole.
        for c in 0..RING_N_290[16] {
            assert!(
                areas.contains(&(RING_BASE_290[16] + c)),
                "polar field missed cell {c} of ring 17"
            );
        }
    }

    #[test]
    fn ra_wrap_is_handled() {
        // Field straddling RA = 0 must take cells on both sides of the seam.
        let areas = find_areas_290(deg(0.5), 0.0, deg(8.0));
        let n_ra = RING_N_290[9];
        let base = RING_BASE_290[9];
        assert!(areas.contains(&base), "missing the first cell of the ring");
        assert!(
            areas.contains(&(base + n_ra - 1)),
            "missing the last cell of the ring (RA wrap)"
        );
    }

    #[test]
    fn all_sky_field_takes_every_area() {
        let areas = find_areas_290(0.0, 0.0, deg(180.0));
        assert_eq!(areas.len(), 290);
    }
}
