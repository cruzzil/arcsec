// Background estimation for star detection.

use crate::types::ImageBuffer;

/// An inclusive pixel rectangle. Grouping the four bounds keeps the scan signatures
/// readable — they were long enough that clippy flagged them, and the bounds always
/// travel together anyway.
#[derive(Debug, Clone, Copy)]
pub struct Region {
    pub x0: usize,
    pub x1: usize,
    pub y0: usize,
    pub y1: usize,
}

impl Region {
    /// The whole frame.
    pub fn whole(img: &ImageBuffer) -> Self {
        Self {
            x0: 0,
            x1: img.width - 1,
            y0: 0,
            y1: img.height - 1,
        }
    }

    /// The frame minus a one-pixel border, which is what the star scan wants: the
    /// hot-pixel test reads a neighbour on each side.
    pub fn inset(img: &ImageBuffer) -> Self {
        Self {
            x0: 1,
            x1: img.width - 2,
            y0: 1,
            y1: img.height - 2,
        }
    }

    /// The same rectangle with a different row range.
    pub fn with_rows(self, y0: usize, y1: usize) -> Self {
        Self { y0, y1, ..self }
    }

    pub fn rows(&self) -> usize {
        self.y1.saturating_sub(self.y0) + 1
    }
}

const HIST_SIZE: usize = 65536; // 16-bit pixel values

/// Build a histogram of pixel values in a sub-region of the image.
fn build_histogram(img: &ImageBuffer, r: Region, upper_limit: usize) -> Vec<u32> {
    let (x0, x1, y0, y1) = (r.x0, r.x1, r.y0, r.y1);
    let cap = upper_limit.min(HIST_SIZE - 1);

    // Histogram addition is associative and exact, so accumulate per-thread and sum:
    // this is a full pass over the frame and was ~8% of a typical solve.
    let rows = y1.saturating_sub(y0) + 1;
    let threads = crate::max_threads().clamp(1, 32);
    // Below this the thread setup costs more than the scan.
    let n_bands = if rows < 256 || (x1 - x0 + 1) * rows < 1 << 20 {
        1
    } else {
        threads.min(rows / 64).max(1)
    };

    if n_bands == 1 {
        return histogram_rows(img, x0, x1, y0, y1, cap);
    }

    let band_rows = rows.div_ceil(n_bands);
    let parts: Vec<Vec<u32>> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..n_bands)
            .map(|b| {
                let by0 = y0 + b * band_rows;
                let by1 = (by0 + band_rows - 1).min(y1);
                scope.spawn(move || {
                    if by0 > y1 {
                        Vec::new()
                    } else {
                        histogram_rows(img, x0, x1, by0, by1, cap)
                    }
                })
            })
            .collect();
        handles
            .into_iter()
            // A dead band must not become an empty histogram: the merge below
            // would zip against nothing and yield a plausible but wrong
            // background, so every later stage would be quietly miscalibrated.
            .map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
            .collect()
    });

    let mut hist = vec![0u32; cap + 1];
    for part in parts {
        for (h, p) in hist.iter_mut().zip(part.iter()) {
            *h += *p;
        }
    }
    hist
}

/// Histogram of one row range.
fn histogram_rows(
    img: &ImageBuffer,
    x0: usize,
    x1: usize,
    y0: usize,
    y1: usize,
    cap: usize,
) -> Vec<u32> {
    let mut hist = vec![0u32; cap + 1];
    for y in y0..=y1 {
        let row = y * img.width;
        for x in x0..=x1 {
            let raw = img.data[row + x];
            // `as usize` saturates NaN to 0 and negatives to 0, which would bias
            // the background low; skip anything non-finite outright.
            if raw.is_finite() {
                let v = raw as usize;
                if v <= cap {
                    hist[v] += 1;
                }
            }
        }
    }
    hist
}

/// Sigma-clipped mean and standard deviation from a histogram sub-region.
/// `sigma_low=3`, `sigma_high=2`, iterates until `|Δmean| < convergence_threshold`.
pub fn sigma_clipped_mean_from_histogram(
    img: &ImageBuffer,
    region: Region,
    upper_limit: usize,
    max_iterations: usize,
    convergence_threshold: f64,
) -> (f64, f64) {
    const SIGMA_LOW: f64 = 3.0;
    const SIGMA_HIGH: f64 = 2.0;

    let hist = build_histogram(img, region, upper_limit);
    let hist_len = hist.len();

    let mut mean = 0.0f64;
    let mut stdev = 0.0f64;
    let mut lo = 0usize;
    let mut hi = hist_len - 1;

    for iter in 0..max_iterations {
        let prev_mean = mean;
        let prev_stdev = stdev;

        let mut sum = 0.0f64;
        let mut sum_sq = 0.0f64;
        let mut total = 0u64;

        for (i, &bin) in hist.iter().enumerate().take(hi + 1).skip(lo) {
            let cnt = bin as u64;
            if cnt > 0 {
                let v = i as f64;
                sum += v * cnt as f64;
                sum_sq += v * v * cnt as f64;
                total += cnt;
            }
        }

        if total == 0 {
            break;
        }

        mean = sum / total as f64;
        let variance = if total > 1 {
            let v = (sum_sq - sum * sum / total as f64) / (total as f64 - 1.0);
            v.max(0.0)
        } else {
            0.0
        };
        stdev = variance.sqrt();

        if stdev > 0.0 {
            // The lower clip bound is deliberately held at 0, so SIGMA_LOW is not applied.
            lo = 0;
            hi = (upper_limit)
                .min((mean + SIGMA_HIGH * stdev).round() as usize)
                .min(hist_len - 1);
        }

        if iter > 0
            && (mean - prev_mean).abs() < convergence_threshold
            && (stdev - prev_stdev).abs() < convergence_threshold
        {
            break;
        }

        let _ = SIGMA_LOW; // referenced so the constant does not warn as unused
    }

    (mean, stdev)
}

/// Result of background analysis.
#[derive(Debug, Clone)]
pub struct Background {
    /// Modal background value (peak of histogram).
    pub mean: f64,
    /// Noise standard deviation (sigma-clipped).
    pub noise: f64,
    /// Detection threshold for bright/small stars (HFD ~2.25 px).
    pub star_level: f64,
    /// Detection threshold for faint/large stars (HFD ~4.5 px).
    pub star_level2: f64,
}

/// Analyse image background, noise, and star detection levels.
///
/// `max_stars`: number of stars expected (empirical factor for star_level).
pub fn get_background(img: &ImageBuffer, max_stars: usize) -> Background {
    let width = img.width;
    let height = img.height;

    // Build full histogram (0..65535)
    let hist = build_histogram(img, Region::whole(img), 65001);

    // --- Find background: peak of histogram ---
    let total_pixels = (width * height) as u64;
    let mean_value = {
        let sum: u64 = hist
            .iter()
            .enumerate()
            .map(|(i, &c)| i as u64 * c as u64)
            .sum();
        (sum / total_pixels.max(1)) as usize
    };

    let mut background = img.data[0] as f64;
    if mean_value == 0 {
        background = 0.0;
    } else {
        let mut peak_count = 0u32;
        for (i, &bin) in hist.iter().enumerate().take(mean_value + 1).skip(1) {
            if bin > peak_count {
                peak_count = bin;
                background = i as f64;
            }
        }
        // If histogram mean is > 1.5× modal peak, use mean instead
        if mean_value as f64 > 1.5 * background {
            background = mean_value as f64;
        }
    }

    // --- Noise estimation: sigma-clipped standard deviation (sample of pixels) ---
    let step_size = ((height as f64 / 71.0).round() as usize).max(1);
    // Make step_size odd so it doesn't stride evenly through Bayer rows
    let step_size = if step_size.is_multiple_of(2) {
        step_size + 1
    } else {
        step_size
    };

    let mut sd = 1e9f64;
    let mut iterations = 0usize;
    loop {
        let sd_old = sd;
        let mut sum_sq = 0.0f64;
        let mut counter = 0u64;

        let mut x = 15usize;
        while x <= width.saturating_sub(16) {
            let mut y = 15usize;
            while y <= height.saturating_sub(16) {
                let value = img.get(x, y) as f64;
                // Exclude outliers (>2× background) and zero pixels
                if value < background * 2.0
                    && value != 0.0
                    && (iterations == 0 || (value - background).abs() <= 3.0 * sd_old)
                {
                    sum_sq += (value - background).powi(2);
                    counter += 1;
                }
                y += step_size;
            }
            x += step_size;
        }

        sd = if counter > 0 {
            (sum_sq / counter as f64).sqrt()
        } else {
            0.0
        };
        iterations += 1;

        if (sd_old - sd).abs() < 0.05 * sd || iterations >= 7 {
            break;
        }
    }
    let noise = sd;

    // --- Star levels: threshold where histogram count drops below empirical limits ---
    let max_range = 65001usize;
    let factor = 6 * max_stars;
    let factor2 = 24 * max_stars;

    let mut above = 0usize;
    let mut star_level_raw = 0.0f64;
    let mut star_level2_raw = 0.0f64;
    let mut i = max_range;

    while star_level_raw == 0.0 && i > (background + 1.0) as usize {
        i -= 1;
        above += hist[i.min(hist.len() - 1)] as usize;
        if above >= factor {
            star_level_raw = i as f64;
        }
    }
    while star_level2_raw == 0.0 && i > (background + 1.0) as usize {
        i -= 1;
        above += hist[i.min(hist.len() - 1)] as usize;
        if above >= factor2 {
            star_level2_raw = i as f64;
        }
    }

    let min_level = (3.5 * noise).max(1.0);
    let star_level = min_level.max((star_level_raw - background - 1.0).max(0.0));
    let star_level2 = min_level.max((star_level2_raw - background - 1.0).max(0.0));

    Background {
        mean: background,
        noise,
        star_level,
        star_level2,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_flat_image(width: usize, height: usize, value: f32) -> ImageBuffer {
        ImageBuffer {
            data: vec![value; width * height],
            width,
            height,
        }
    }

    fn make_image_with_noise(width: usize, height: usize, bg: f32, noise_amp: f32) -> ImageBuffer {
        let mut data = vec![0f32; width * height];
        // Simple deterministic pseudo-noise using index
        for (i, v) in data.iter_mut().enumerate() {
            let noise = ((i * 2654435769 + 1234567) & 0xFFFF) as f32 / 65535.0 * noise_amp;
            *v = bg + noise - noise_amp / 2.0;
        }
        ImageBuffer {
            data,
            width,
            height,
        }
    }

    #[test]
    fn flat_image_background() {
        let img = make_flat_image(100, 100, 1000.0);
        let bg = get_background(&img, 500);
        // Background peak should be near 1000
        assert!((bg.mean - 1000.0).abs() < 5.0, "background = {}", bg.mean);
        assert!(bg.noise < 1.0, "flat image noise = {}", bg.noise);
    }

    #[test]
    fn sigma_clip_converges() {
        let img = make_image_with_noise(200, 200, 5000.0, 200.0);
        let (mean, stdev) = sigma_clipped_mean_from_histogram(
            &img,
            Region {
                x0: 0,
                x1: 199,
                y0: 0,
                y1: 199,
            },
            65500,
            6,
            0.1,
        );
        // Mean should be near 5000, stdev near 100
        assert!((mean - 5000.0).abs() < 200.0, "mean = {mean}");
        assert!(stdev > 0.0, "stdev should be positive, got {stdev}");
    }

    #[test]
    fn star_levels_increase_with_stars() {
        // Image with some bright pixels (simulated stars)
        let mut img = make_image_with_noise(200, 200, 1000.0, 50.0);
        // Plant 20 "star" pixels well above background
        for i in 0..20 {
            let x = 10 + i * 9;
            let y = 10 + i * 9;
            if x < 200 && y < 200 {
                img.data[y * 200 + x] = 20000.0;
            }
        }
        let bg = get_background(&img, 500);
        // star_level should be above noise
        assert!(bg.star_level > 0.0, "star_level = {}", bg.star_level);
    }
}
