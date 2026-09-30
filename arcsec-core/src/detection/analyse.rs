//! Image analysis without solving: star count, median HFD and a per-star table.
//!
//! This is ASTAP's `analyse_image`, behind its `-analyse`, `-extract` and `-extract2`
//! options. It differs from the solver's detection ([`crate::detection::stars`]) in
//! ways that matter to anyone comparing the two:
//!
//! - it works on the full-resolution image, never a binned one;
//! - the minimum SNR is the caller's (`snr_min`), not the solver's fixed 10, and there
//!   is no minimum HFD beyond the 0.8 px hot-pixel floor;
//! - each detection pass starts afresh. The solver accumulates stars over its passes,
//!   whereas here a pass that finds too few stars is discarded and the next, lower,
//!   threshold scans the whole image again;
//! - the last pass uses a single threshold of `max(snr_min, 7)` × noise over the
//!   whole frame, not the solver's grid of local backgrounds;
//! - nothing is trimmed: every star the final pass finds is reported.
//!
//! The per-star measurement is shared with the solver ([`measure_star_with_flux`]).

use crate::detection::background::{Background, get_background};
use crate::detection::stars::{measure_star_with_flux, median_f64};
use crate::types::ImageBuffer;

/// A star found by [`analyse_image`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeasuredStar {
    /// Sub-pixel column of the centroid (0-based).
    pub x: f64,
    /// Sub-pixel row of the centroid (0-based).
    pub y: f64,
    /// Half-flux diameter, pixels.
    pub hfd: f64,
    /// Signal-to-noise ratio of the aperture flux.
    pub snr: f64,
    /// Background-subtracted flux in the measuring aperture, in pixel units (ADU).
    pub flux: f64,
}

/// The result of [`analyse_image`].
#[derive(Debug, Clone)]
pub struct Analysis {
    /// Every star the final detection pass found, in scan order (row by row).
    pub stars: Vec<MeasuredStar>,
    /// The background and noise the detection thresholds were derived from.
    pub background: Background,
}

impl Analysis {
    /// Median HFD of the stars, in pixels, or `None` if there are none.
    #[must_use]
    pub fn hfd_median(&self) -> Option<f64> {
        if self.stars.is_empty() {
            return None;
        }
        let mut hfds: Vec<f64> = self.stars.iter().map(|s| s.hfd).collect();
        Some(median_f64(&mut hfds))
    }
}

/// Find and measure the stars of `img` without solving it.
///
/// Detection runs in up to four passes at falling thresholds, exactly as ASTAP's
/// `analyse_image` does:
///
/// 1. the bright-star level from the histogram, if it is above 30 × noise;
/// 2. the fainter histogram level, on the same condition;
/// 3. 30 × noise, skipped when `snr_min` is 30 or more;
/// 4. `max(snr_min, 7)` × noise.
///
/// It stops at the first pass that finds `max_stars` stars or more (or after the
/// last), and returns that pass's stars. A star counts if its SNR exceeds `snr_min`
/// and its HFD lies in (0.8, 30] pixels.
///
/// Single-threaded; images smaller than 3×3 yield no stars.
#[must_use]
pub fn analyse_image(img: &ImageBuffer, snr_min: f64, max_stars: usize) -> Analysis {
    let background = get_background(img, max_stars);
    let (w, h) = (img.width, img.height);
    if w < 3 || h < 3 || img.data.len() < w * h {
        return Analysis {
            stars: Vec::new(),
            background,
        };
    }

    let noise = background.noise;
    let mut retries = 4u8;
    let stars = loop {
        let mut level = background.star_level;
        if retries == 4 && background.star_level <= 30.0 * noise {
            retries = 3;
        }
        if retries == 3 {
            if background.star_level2 > 30.0 * noise {
                level = background.star_level2;
            } else {
                retries = 2;
            }
        }
        if retries == 2 {
            level = 30.0 * noise;
            if snr_min >= 30.0 {
                retries = 1;
            }
        }
        if retries == 1 {
            level = snr_min.max(7.0) * noise;
        }

        let found = scan(img, &background, level, snr_min);
        retries -= 1;
        if found.len() >= max_stars || retries == 0 {
            break found;
        }
    };

    Analysis { stars, background }
}

/// One detection pass over the whole image (less a one-pixel border) at
/// `level` above the background.
fn scan(img: &ImageBuffer, bg: &Background, level: f64, snr_min: f64) -> Vec<MeasuredStar> {
    let (w, h) = (img.width, img.height);
    let detect_abs = bg.mean + level;
    let hot_abs = bg.mean + 4.0 * bg.noise;
    // Pixels already claimed by a star found in this pass.
    let mut taken = vec![false; w * h];
    let mut out = Vec::new();

    for fy in 1..h - 1 {
        let row = fy * w;
        for fx in 1..w - 1 {
            if (img.data[row + fx] as f64) <= detect_abs || taken[row + fx] {
                continue;
            }
            // A hot pixel stands alone: a star lights at least two of the four
            // pixels around it.
            let lit = [row + fx - 1, row + fx + 1, row - w + fx, row + w + fx]
                .iter()
                .filter(|&&i| img.data[i] as f64 > hot_abs)
                .count();
            if lit < 2 {
                continue;
            }

            let Some((star, flux)) = measure_star_with_flux(img, fx as i32, fy as i32) else {
                continue;
            };
            if !(star.hfd <= 30.0 && star.snr > snr_min && star.hfd > 0.8) {
                continue;
            }
            let xci = star.x.round() as i64;
            let yci = star.y.round() as i64;
            if xci >= 0
                && yci >= 0
                && (xci as usize) < w
                && (yci as usize) < h
                && taken[yci as usize * w + xci as usize]
            {
                continue; // the same star, reached again from another of its pixels
            }

            // Claim the star's disc so its other pixels do not seed it again.
            let radius = (3.0 * star.hfd).round() as i64;
            for n in -radius..=radius {
                for m in -radius..=radius {
                    let (xi, yi) = (xci + m, yci + n);
                    if xi >= 0
                        && yi >= 0
                        && (xi as usize) < w
                        && (yi as usize) < h
                        && m * m + n * n <= radius * radius
                    {
                        taken[yi as usize * w + xi as usize] = true;
                    }
                }
            }
            out.push(MeasuredStar {
                x: star.x,
                y: star.y,
                hfd: star.hfd,
                snr: star.snr,
                flux,
            });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Gaussian star added into `img`.
    fn add_star(img: &mut ImageBuffer, cx: f64, cy: f64, sigma: f64, peak: f64) {
        let rs = (5.0 * sigma).ceil() as i64;
        for dy in -rs..=rs {
            for dx in -rs..=rs {
                let (x, y) = (cx.round() as i64 + dx, cy.round() as i64 + dy);
                if x < 0 || y < 0 || x as usize >= img.width || y as usize >= img.height {
                    continue;
                }
                let r2 = (x as f64 - cx).powi(2) + (y as f64 - cy).powi(2);
                img.data[y as usize * img.width + x as usize] +=
                    (peak * (-r2 / (2.0 * sigma * sigma)).exp()) as f32;
            }
        }
    }

    fn noisy(w: usize, h: usize, bg: f64, sigma: f64) -> ImageBuffer {
        let mut rng = crate::test_support::Rng::new(3);
        ImageBuffer {
            data: (0..w * h)
                .map(|_| (bg + sigma * rng.gauss()) as f32)
                .collect(),
            width: w,
            height: h,
        }
    }

    #[test]
    fn finds_every_star_with_its_hfd_and_position() {
        let mut img = noisy(400, 300, 1000.0, 10.0);
        let planted = [
            (50.3, 60.7),
            (200.0, 150.0),
            (330.6, 40.2),
            (120.0, 250.5),
            (300.0, 220.0),
        ];
        for &(x, y) in &planted {
            add_star(&mut img, x, y, 1.5, 5000.0);
        }
        let a = analyse_image(&img, 30.0, 500);
        assert_eq!(a.stars.len(), planted.len(), "{:?}", a.stars);
        for s in &a.stars {
            assert!(
                planted
                    .iter()
                    .any(|&(x, y)| (s.x - x).abs() < 0.2 && (s.y - y).abs() < 0.2),
                "{s:?} is not a planted star"
            );
            // A Gaussian's HFD is about 2.35 sigma.
            assert!((s.hfd - 2.35 * 1.5).abs() < 0.6, "hfd {}", s.hfd);
            assert!(s.snr > 30.0 && s.flux > 0.0);
        }
        // Scan order: row by row.
        assert!(a.stars.windows(2).all(|p| p[0].y.round() <= p[1].y.round()));
        let median = a.hfd_median().unwrap();
        assert!((median - 2.35 * 1.5).abs() < 0.6, "median {median}");
    }

    #[test]
    fn snr_min_sets_the_faintest_star_reported() {
        let mut img = noisy(300, 300, 1000.0, 10.0);
        add_star(&mut img, 80.0, 80.0, 1.5, 5000.0);
        // Peak 12 sigma: above the 7 sigma last-pass threshold, but an aperture SNR
        // of about 25.
        add_star(&mut img, 200.0, 200.0, 2.0, 120.0);
        let strict = analyse_image(&img, 30.0, 500);
        assert_eq!(strict.stars.len(), 1);
        let loose = analyse_image(&img, 5.0, 500);
        assert_eq!(loose.stars.len(), 2, "{:?}", loose.stars);
    }

    #[test]
    fn a_blank_frame_has_no_stars_and_no_median() {
        let img = noisy(200, 200, 1000.0, 10.0);
        let a = analyse_image(&img, 30.0, 500);
        assert!(a.stars.is_empty());
        assert_eq!(a.hfd_median(), None);
        // Too small to scan at all.
        assert!(
            analyse_image(&ImageBuffer::new(2, 2), 30.0, 500)
                .stars
                .is_empty()
        );
    }

    #[test]
    fn a_hot_pixel_is_not_a_star() {
        let mut img = noisy(100, 100, 1000.0, 10.0);
        img.data[50 * 100 + 50] = 60000.0;
        assert!(analyse_image(&img, 10.0, 500).stars.is_empty());
    }
}
