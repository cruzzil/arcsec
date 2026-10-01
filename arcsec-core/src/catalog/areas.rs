//! Sky area lookup for the .1476 catalog format.
//!
//! The sky is divided into 1476 tiles of ~5° × 5°.
//! Rows of constant DEC are indexed 0 (south pole) to 35 (north pole).
//! Each row is subdivided into `n_ra_cells` RA cells (see `RING_TABLE`).

use core::f64::consts::PI;

/// DEC boundaries for the 1476 grid (radians, 37 values from -π/2 to +π/2).
/// Index 0 = -90°, index 36 = +90°.  Adjacent pair defines one DEC ring.
pub const DEC_BOUNDARIES_1476: [f64; 37] = [
    -90.0_f64 * PI / 180.0,
    -87.42857143 * PI / 180.0,
    -82.28571429 * PI / 180.0,
    -77.14285714 * PI / 180.0,
    -72.0 * PI / 180.0,
    -66.85714286 * PI / 180.0,
    -61.71428571 * PI / 180.0,
    -56.57142857 * PI / 180.0,
    -51.42857143 * PI / 180.0,
    -46.28571429 * PI / 180.0,
    -41.14285714 * PI / 180.0,
    -36.0 * PI / 180.0,
    -30.85714286 * PI / 180.0,
    -25.71428571 * PI / 180.0,
    -20.57142857 * PI / 180.0,
    -15.42857143 * PI / 180.0,
    -10.28571429 * PI / 180.0,
    -5.142857143 * PI / 180.0,
    0.0,
    5.142857143 * PI / 180.0,
    10.28571429 * PI / 180.0,
    15.42857143 * PI / 180.0,
    20.57142857 * PI / 180.0,
    25.71428571 * PI / 180.0,
    30.85714286 * PI / 180.0,
    36.0 * PI / 180.0,
    41.14285714 * PI / 180.0,
    46.28571429 * PI / 180.0,
    51.42857143 * PI / 180.0,
    56.57142857 * PI / 180.0,
    61.71428571 * PI / 180.0,
    66.85714286 * PI / 180.0,
    72.0 * PI / 180.0,
    77.14285714 * PI / 180.0,
    82.28571429 * PI / 180.0,
    87.42857143 * PI / 180.0,
    90.0_f64 * PI / 180.0,
];

/// Per-ring table: (`n_ra_cells`, `base_area_1indexed`).
/// Index 0 = south pole (1 cell), index 35 = north pole (1 cell).
/// All other rings subdivide the RA circle into `n_ra_cells` equal cells.
const RING_TABLE: [(usize, usize); 36] = [
    (1, 1),     // ring  0: south pole
    (3, 2),     // ring  1
    (9, 5),     // ring  2
    (15, 14),   // ring  3
    (21, 29),   // ring  4
    (27, 50),   // ring  5
    (33, 77),   // ring  6
    (38, 110),  // ring  7
    (43, 148),  // ring  8
    (48, 191),  // ring  9
    (52, 239),  // ring 10
    (56, 291),  // ring 11
    (60, 347),  // ring 12
    (63, 407),  // ring 13
    (65, 470),  // ring 14
    (67, 535),  // ring 15
    (68, 602),  // ring 16
    (69, 670),  // ring 17
    (69, 739),  // ring 18
    (68, 808),  // ring 19
    (67, 876),  // ring 20
    (65, 943),  // ring 21
    (63, 1008), // ring 22
    (60, 1071), // ring 23
    (56, 1131), // ring 24
    (52, 1187), // ring 25
    (48, 1239), // ring 26
    (43, 1287), // ring 27
    (38, 1330), // ring 28
    (33, 1368), // ring 29
    (27, 1401), // ring 30
    (21, 1428), // ring 31
    (15, 1449), // ring 32
    (9, 1464),  // ring 33
    (3, 1473),  // ring 34
    (1, 1476),  // ring 35: north pole
];

/// Boundary distances from a point (ra, dec) to the 4 edges of its area cell.
#[derive(Debug, Clone, Copy)]
pub struct AreaBounds {
    /// 1-indexed area number in the 1476 grid.
    pub area_nr: usize,
    /// Distance to the cell's eastern edge (radians on the sky).
    pub space_east: f64,
    /// Distance to the cell's western edge (radians on the sky).
    pub space_west: f64,
    /// Distance to the cell's northern edge (radians).
    pub space_north: f64,
    /// Distance to the cell's southern edge (radians).
    pub space_south: f64,
}

/// Find the 1476 area number and bounding distances for a given (ra, dec).
#[must_use]
pub fn area_and_boundaries_1476(ra: f64, dec: f64) -> AreaBounds {
    let cos_dec = dec.cos();

    // North pole
    if dec > DEC_BOUNDARIES_1476[35] {
        return AreaBounds {
            area_nr: 1476,
            space_east: PI * 2.0,
            space_west: PI * 2.0,
            space_north: DEC_BOUNDARIES_1476[36] - DEC_BOUNDARIES_1476[35],
            space_south: dec - DEC_BOUNDARIES_1476[35],
        };
    }

    // Rings 34 down to 1, then fallthrough to south pole
    // Iterate north to south: the first ring whose boundary lies below `dec` is the match.
    for ring in (1usize..=34).rev() {
        if dec > DEC_BOUNDARIES_1476[ring] {
            let (n_ra, base) = RING_TABLE[ring];
            let rot = ra * n_ra as f64 / (2.0 * PI);
            let area_nr = base + rot.floor() as usize;
            let frac = rot.fract();
            let ra_step = 2.0 * PI / n_ra as f64;

            let north_boundary = DEC_BOUNDARIES_1476[ring + 1];

            return AreaBounds {
                area_nr,
                space_east: ra_step * (1.0 - frac) * cos_dec,
                space_west: ra_step * frac * cos_dec,
                space_north: north_boundary - dec,
                space_south: dec - DEC_BOUNDARIES_1476[ring],
            };
        }
    }

    // South pole fallback
    AreaBounds {
        area_nr: 1,
        space_east: PI * 2.0,
        space_west: PI * 2.0,
        space_north: DEC_BOUNDARIES_1476[1] - dec,
        space_south: DEC_BOUNDARIES_1476[1] - DEC_BOUNDARIES_1476[0],
    }
}

/// Every 1476 area whose declination ring overlaps `[dec_lo, dec_hi]` (radians):
/// whole rings, in ascending area order. For reading a declination band of the
/// whole sky, as the blind-index builder does.
#[must_use]
pub fn areas_in_dec_band_1476(dec_lo: f64, dec_hi: f64) -> Vec<usize> {
    let mut out = Vec::new();
    for (ring, &(n_ra, base)) in RING_TABLE.iter().enumerate() {
        if DEC_BOUNDARIES_1476[ring + 1] >= dec_lo && DEC_BOUNDARIES_1476[ring] <= dec_hi {
            out.extend(base..base + n_ra);
        }
    }
    out
}

/// Return the filename segment for a 1476 area number (e.g. area 1 → "0101.1476").
/// Filename format: `{ring:02}{cell:02}.1476` where ring and cell are 1-based.
#[must_use]
pub fn filename_1476(area_nr: usize) -> String {
    let area = area_nr.clamp(1, 1476);

    // Find the ring by searching the RING_TABLE
    let mut ring_idx = 0usize;
    for (r, &(n_ra, base)) in RING_TABLE.iter().enumerate() {
        if area >= base && area < base + n_ra {
            ring_idx = r;
            break;
        }
        if r == 35 {
            ring_idx = 35; // north pole
        }
    }

    let (n_ra, base) = RING_TABLE[ring_idx];
    let cell_1indexed = if n_ra == 1 { 1 } else { area - base + 1 };
    format!("{:02}{:02}.1476", ring_idx + 1, cell_1indexed)
}

/// Find up to 4 distinct 1476 areas that overlap a square FOV.
///
/// `fov` is the side length of the image in radians (must be ≤ 5.14°).
///
/// Returns `(area_nr, fraction)` pairs, where fraction is the approximate fraction
/// of the FOV covered by that area.  Areas with fraction < 0.01 are omitted.
/// Duplicate areas are deduplicated.
#[must_use]
pub fn find_areas_1476(ra: f64, dec: f64, fov: f64) -> Vec<(usize, f64)> {
    let fov = fov.min(5.142857_f64.to_radians()); // one ring height: four-corner sampling cannot express a taller field
    let fov_half = fov * 0.5;

    let dec_n = dec + fov_half;
    let dec_s = dec - fov_half;

    let cos_n = dec_n.cos().max(1e-6);
    let cos_s = dec_s.cos().max(1e-6);

    let mut ra_wn = ra - fov_half / cos_n;
    if ra_wn < 0.0 {
        ra_wn += 2.0 * PI;
    }
    let mut ra_en = ra + fov_half / cos_n;
    if ra_en >= 2.0 * PI {
        ra_en -= 2.0 * PI;
    }
    let mut ra_ws = ra - fov_half / cos_s;
    if ra_ws < 0.0 {
        ra_ws += 2.0 * PI;
    }
    let mut ra_es = ra + fov_half / cos_s;
    if ra_es >= 2.0 * PI {
        ra_es -= 2.0 * PI;
    }

    // 4 corners: NE, NW, SE, SW
    let corners = [
        (ra_en, dec_n),
        (ra_wn, dec_n),
        (ra_es, dec_s),
        (ra_ws, dec_s),
    ];

    let fov2 = fov * fov;
    let corner_fracs = |_b: AreaBounds, space_h_toward: f64, space_v_toward: f64| -> f64 {
        let h = space_h_toward.min(fov);
        let v = space_v_toward.min(fov);
        h * v / fov2
    };

    let b: Vec<AreaBounds> = corners
        .iter()
        .map(|&(r, d)| area_and_boundaries_1476(r, d))
        .collect();

    // Fractions: corner coverage fractions towards the image centre
    // NE corner: coverage toward W and S
    let f0 = corner_fracs(b[0], b[0].space_west, b[0].space_south);
    // NW corner: coverage toward E and S
    let f1 = corner_fracs(b[1], b[1].space_east, b[1].space_south);
    // SE corner: coverage toward W and N
    let f2 = corner_fracs(b[2], b[2].space_west, b[2].space_north);
    // SW corner: coverage toward E and N
    let f3 = corner_fracs(b[3], b[3].space_east, b[3].space_north);

    let raw = [
        (b[0].area_nr, f0),
        (b[1].area_nr, f1),
        (b[2].area_nr, f2),
        (b[3].area_nr, f3),
    ];

    // Dedup: same area → keep only first entry (merge fractions implicitly via the sum)
    let mut seen = [0usize; 4];
    let mut result = Vec::with_capacity(4);
    for (area, frac) in raw {
        if frac < 0.01 {
            continue;
        }
        if seen.contains(&area) {
            continue;
        }
        seen[result.len()] = area;
        result.push((area, frac));
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn deg(d: f64) -> f64 {
        d * PI / 180.0
    }

    #[test]
    fn south_pole_area_is_1() {
        let b = area_and_boundaries_1476(0.0, deg(-89.0));
        assert_eq!(b.area_nr, 1);
    }

    #[test]
    fn north_pole_area_is_1476() {
        let b = area_and_boundaries_1476(0.0, deg(89.0));
        assert_eq!(b.area_nr, 1476);
    }

    #[test]
    fn equator_area_in_valid_range() {
        let b = area_and_boundaries_1476(deg(180.0), 0.0);
        assert!(b.area_nr >= 670 && b.area_nr <= 807, "area={}", b.area_nr);
    }

    #[test]
    fn filename_south_pole() {
        assert_eq!(filename_1476(1), "0101.1476");
    }

    #[test]
    fn filename_ring1_cell2() {
        // Ring 1 starts at area 2; cell 2 = area 3
        assert_eq!(filename_1476(3), "0202.1476");
    }

    #[test]
    fn filename_north_pole() {
        assert_eq!(filename_1476(1476), "3601.1476");
    }

    #[test]
    fn find_areas_small_fov_returns_at_least_one() {
        let areas = find_areas_1476(deg(45.0), deg(20.0), deg(2.0));
        assert!(!areas.is_empty(), "should find at least one area");
        for (a, f) in &areas {
            assert!(*a >= 1 && *a <= 1476, "area out of range: {a}");
            assert!(*f >= 0.01, "fraction {f} below threshold");
        }
    }

    #[test]
    fn find_areas_no_duplicates() {
        let areas = find_areas_1476(deg(0.0), deg(0.0), deg(4.0));
        let mut seen = std::collections::HashSet::new();
        for (a, _) in areas {
            assert!(seen.insert(a), "duplicate area {a}");
        }
    }

    #[test]
    fn area_nr_agrees_with_ring_offsets() {
        // Area 2 = ring 1, cell 1 → '0201.1476'
        assert_eq!(filename_1476(2), "0201.1476");
        // Area 5 = ring 2, cell 1 → '0301.1476'
        assert_eq!(filename_1476(5), "0301.1476");
        // Area 14 = ring 3, cell 1 → '0401.1476'
        assert_eq!(filename_1476(14), "0401.1476");
    }

    #[test]
    fn boundary_distances_positive() {
        // For a point well inside a ring the boundary distances should all be positive
        let b = area_and_boundaries_1476(deg(90.0), deg(30.0));
        assert!(b.space_north > 0.0, "space_north = {}", b.space_north);
        assert!(b.space_south > 0.0, "space_south = {}", b.space_south);
        assert!(b.space_east > 0.0, "space_east  = {}", b.space_east);
        assert!(b.space_west > 0.0, "space_west  = {}", b.space_west);
    }
}
