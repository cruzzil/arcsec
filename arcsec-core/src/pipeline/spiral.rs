//! Square spiral iterator.
//! Sequence starting at (0,0): (1,0), (1,1), (0,1), (-1,1), (-1,0), (-1,-1),
//! (0,-1), (1,-1), (2,-1), (2,0), ...
//!
//! The iterator stops when `spiral_x > max_distance`.

/// Iterates (`spiral_x`, `spiral_y`) over a square spiral centred at the origin.
/// Stops when `spiral_x` exceeds `max_distance`.
pub struct SpiralSearch {
    x: i32,
    y: i32,
    dx: i32,
    dy: i32,
    count: usize,
    max_distance: i32,
    done: bool,
}

impl SpiralSearch {
    /// A spiral that stops once `spiral_x` would exceed `max_distance`.
    #[must_use]
    pub fn new(max_distance: i32) -> Self {
        Self {
            x: 0,
            y: 0,
            dx: 0,
            dy: -1,
            count: 0,
            max_distance,
            done: false,
        }
    }
}

impl Iterator for SpiralSearch {
    type Item = (i32, i32);

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        if self.count != 0 {
            // Turning-point check (ASTAP condition exactly)
            if self.x == self.y
                || (self.x < 0 && self.x == -self.y)
                || (self.x > 0 && self.x == 1 - self.y)
            {
                let t = self.dx;
                self.dx = -self.dy;
                self.dy = t;
            }
            self.x += self.dx;
            self.y += self.dy;
        }

        self.count += 1;

        if self.x > self.max_distance {
            self.done = true;
            return None;
        }

        Some((self.x, self.y))
    }
}

/// The number of positions [`SpiralSearch::new`]`(max_distance)` yields: every ring
/// out to `max_distance`, `(2m + 1)²`, or none for a negative distance.
///
/// Computed rather than counted so the solver need not hold the spiral in memory:
/// a search radius of many fields (or a tiny field from a corrupt header, which
/// saturates `max_distance` at `i32::MAX`) would otherwise allocate for billions
/// of positions before trying the first.
#[must_use]
pub fn spiral_len(max_distance: i32) -> u64 {
    match u64::try_from(max_distance) {
        Ok(m) => (2 * m + 1) * (2 * m + 1),
        Err(_) => 0,
    }
}

/// The position [`SpiralSearch`] yields at `index` (0-based), without walking
/// the spiral to it.
///
/// Ring `k ≥ 1` holds indices `(2k-1)² .. (2k+1)²` and is walked as four sides of
/// `2k` positions each: up the east side from `(k, 1-k)`, west along the north
/// side, down the west side and east along the south side, ending at `(k, -k)`.
#[must_use]
pub fn spiral_position(index: u64) -> (i32, i32) {
    if index == 0 {
        return (0, 0);
    }
    // In u64 until the ring's start is subtracted: (2k - 1)² overflows i64 for the
    // outermost rings an index can name.
    let ku = index.isqrt().div_ceil(2);
    let offset = index - (2 * ku - 1) * (2 * ku - 1);
    let (side, t) = (offset / (2 * ku), (offset % (2 * ku)) as i64);
    let k = ku as i64;
    let (x, y) = match side {
        0 => (k, 1 - k + t),
        1 => (k - 1 - t, k),
        2 => (-k, k - 1 - t),
        _ => (1 - k + t, -k),
    };
    // Every ring a `u64` index can reach (k < 2³²) fits; the solver's spiral
    // never goes past `i32::MAX` rings.
    (x as i32, y as i32)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn closed_form_matches_the_walk() {
        for max in -2..=40 {
            let walked: Vec<(i32, i32)> = SpiralSearch::new(max).collect();
            assert_eq!(walked.len() as u64, spiral_len(max), "max {max}");
            for (i, &p) in walked.iter().enumerate() {
                assert_eq!(spiral_position(i as u64), p, "index {i}");
            }
        }
    }

    #[test]
    fn a_saturated_spiral_is_counted_not_built() {
        // What a field of 1e-300 rad gives: (radius / fov + 2) as i32.
        let m = (1.0f64 / 1e-300 + 2.0) as i32;
        assert_eq!(m, i32::MAX);
        let n = spiral_len(m);
        assert_eq!(n, 4_294_967_295u64 * 4_294_967_295);
        // The last position is the south-east corner of the outermost ring.
        assert_eq!(spiral_position(n - 1), (i32::MAX, -i32::MAX));
    }

    #[test]
    fn first_few_positions() {
        let v: Vec<_> = SpiralSearch::new(2).collect();
        // First item is (0,0), then (1,0), (1,1), (0,1), (-1,1), (-1,0), (-1,-1), (0,-1),
        // (1,-1), (2,-1), (2,0), ...
        assert_eq!(v[0], (0, 0));
        assert_eq!(v[1], (1, 0));
        assert_eq!(v[2], (1, 1));
        assert_eq!(v[3], (0, 1));
        assert_eq!(v[4], (-1, 1));
        assert_eq!(v[5], (-1, 0));
        assert_eq!(v[6], (-1, -1));
        assert_eq!(v[7], (0, -1));
        assert_eq!(v[8], (1, -1));
    }

    #[test]
    fn stops_after_max_distance() {
        let max = 3;
        let positions: Vec<_> = SpiralSearch::new(max).collect();
        for &(x, _) in &positions {
            assert!(x <= max, "x={x} exceeded max_distance={max}");
        }
    }

    #[test]
    fn max_distance_zero_gives_only_origin() {
        let v: Vec<_> = SpiralSearch::new(0).collect();
        assert_eq!(v, vec![(0, 0)]);
    }

    #[test]
    fn covers_expected_area() {
        // A spiral with max_distance=2 should cover all points with x in [-2, 2]
        let positions: Vec<_> = SpiralSearch::new(2).collect();
        // All positions must have x in [-2, 2]
        for &(x, y) in &positions {
            assert!((-2..=2).contains(&x), "x={x}");
            let _ = y; // y range not strictly limited by ASTAP
        }
        // Should start at origin
        assert_eq!(positions[0], (0, 0));
    }
}
