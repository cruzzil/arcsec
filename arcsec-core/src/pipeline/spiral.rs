// Square spiral iterator.
// Sequence starting at (0,0): (1,0), (1,1), (0,1), (-1,1), (-1,0), (-1,-1),
// (0,-1), (1,-1), (2,-1), (2,0), ...
//
// The iterator stops when `spiral_x > max_distance`.

/// Iterates (spiral_x, spiral_y) over a square spiral centred at the origin.
/// Stops when spiral_x exceeds `max_distance`.
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

#[cfg(test)]
mod tests {
    use super::*;

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
