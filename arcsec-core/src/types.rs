/// Row-major pixel buffer. Access: `data[y * width + x]`
#[derive(Debug, Clone)]
pub struct ImageBuffer {
    pub data: Vec<f32>,
    pub width: usize,
    pub height: usize,
}

impl ImageBuffer {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            data: vec![0.0; width * height],
            width,
            height,
        }
    }

    #[inline]
    pub fn get(&self, x: usize, y: usize) -> f32 {
        self.data[y * self.width + x]
    }

    #[inline]
    pub fn get_checked(&self, x: i32, y: i32) -> Option<f32> {
        if x >= 0 && y >= 0 && (x as usize) < self.width && (y as usize) < self.height {
            Some(self.data[y as usize * self.width + x as usize])
        } else {
            None
        }
    }

    /// Replace non-finite pixels and rescale the data into a range the
    /// histogram-based background estimator can actually resolve.
    ///
    /// `detection::get_background` bins pixel values into a 65536-entry integer
    /// histogram (`value as usize`). That silently assumes the data is already
    /// 16-bit-ADU-like. Float FITS in physical units breaks the assumption:
    /// SDSS frames are nanomaggies spanning roughly [-0.15, 4.5] and DESI Legacy
    /// cutouts [-0.008, 27], so every pixel lands in bins 0..4, the background and
    /// noise estimates collapse to zero, and no stars are found at all.
    ///
    /// Returns the (scale, offset) actually applied, or `None` when the data was
    /// left untouched.
    ///
    /// Rules:
    /// * Non-finite pixels (NaN / ±inf — routine in drizzled, reprojected and
    ///   edge-of-survey data) are replaced with the finite minimum, i.e. treated
    ///   as background rather than as sources.
    /// * Data already spanning ≥ `NATIVE_RANGE` counts is left exactly as-is, so
    ///   camera output and 16-bit survey images keep their current behaviour.
    /// * Otherwise the 99.9th percentile is mapped to ~20000 counts, which leaves
    ///   the background and noise several hundred counts wide. Pixels above that
    ///   percentile simply extend past 20000; nothing is clipped.
    pub fn normalize_for_detection(&mut self) -> Option<(f32, f32)> {
        /// Span (in counts) above which the data is assumed to be ADU-like already.
        const NATIVE_RANGE: f32 = 4096.0;
        /// Where the 99.9th percentile lands after rescaling.
        const TARGET_P999: f32 = 20000.0;
        /// Floor value for the rescaled minimum, keeping everything positive.
        const FLOOR: f32 = 100.0;

        // min/max/non-finite count is an exact reduction, so split it across threads:
        // this is a full pass over the frame before anything else happens.
        let n_threads = crate::max_threads().clamp(1, 32);
        let chunk = self.data.len().div_ceil(n_threads.max(1)).max(1 << 18);
        let (mut lo, mut hi, mut n_bad) = (f32::INFINITY, f32::NEG_INFINITY, 0usize);
        if self.data.len() <= chunk {
            for &v in &self.data {
                if v.is_finite() {
                    if v < lo {
                        lo = v;
                    }
                    if v > hi {
                        hi = v;
                    }
                } else {
                    n_bad += 1;
                }
            }
        } else {
            let parts: Vec<(f32, f32, usize)> = std::thread::scope(|scope| {
                let handles: Vec<_> = self
                    .data
                    .chunks(chunk)
                    .map(|c| {
                        scope.spawn(move || {
                            let (mut l, mut h, mut b) = (f32::INFINITY, f32::NEG_INFINITY, 0usize);
                            for &v in c {
                                if v.is_finite() {
                                    if v < l {
                                        l = v;
                                    }
                                    if v > h {
                                        h = v;
                                    }
                                } else {
                                    b += 1;
                                }
                            }
                            (l, h, b)
                        })
                    })
                    .collect();
                handles
                    .into_iter()
                    // A dead chunk must not report n_bad = 0: if it held the
                    // only non-finite pixels, the caller would skip the
                    // replacement pass and feed NaNs straight into detection.
                    .map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
                    .collect()
            });
            for (l, h, b) in parts {
                lo = lo.min(l);
                hi = hi.max(h);
                n_bad += b;
            }
        }
        if !lo.is_finite() || !hi.is_finite() {
            // Nothing finite at all — zero the buffer so downstream code sees a
            // flat frame instead of NaN.
            self.data.iter_mut().for_each(|v| *v = 0.0);
            return None;
        }
        if n_bad > 0 {
            self.data
                .iter_mut()
                .filter(|v| !v.is_finite())
                .for_each(|v| *v = lo);
        }

        if hi - lo >= NATIVE_RANGE {
            return None; // already ADU-like; leave the values alone
        }

        // 99.9th percentile from a strided sample — the maximum is usually a
        // saturated star and would waste most of the dynamic range on it.
        let stride = (self.data.len() / 200_000).max(1);
        let mut sample: Vec<f32> = self.data.iter().step_by(stride).copied().collect();
        if sample.is_empty() {
            return None;
        }
        sample.sort_unstable_by(|a, b| a.total_cmp(b));
        let p999 = sample[((sample.len() as f64 * 0.999) as usize).min(sample.len() - 1)];

        let span = p999 - lo;
        if span <= 0.0 {
            return None;
        }
        let scale = TARGET_P999 / span;
        let offset = FLOOR - lo * scale;
        for v in self.data.iter_mut() {
            *v = *v * scale + offset;
        }
        Some((scale, offset))
    }

    /// Downsample by averaging N×N pixel blocks. Returns a clone if factor ≤ 1.
    pub fn bin_image(&self, factor: usize) -> ImageBuffer {
        if factor <= 1 {
            return self.clone();
        }
        let new_w = self.width / factor;
        let new_h = self.height / factor;
        let fac_sq = (factor * factor) as f32;
        let mut data = vec![0.0f32; new_w * new_h];
        for ny in 0..new_h {
            for nx in 0..new_w {
                let mut sum = 0.0f32;
                for dy in 0..factor {
                    for dx in 0..factor {
                        let ox = nx * factor + dx;
                        let oy = ny * factor + dy;
                        sum += self.data[oy * self.width + ox];
                    }
                }
                data[ny * new_w + nx] = sum / fac_sq;
            }
        }
        ImageBuffer {
            data,
            width: new_w,
            height: new_h,
        }
    }
}

/// A detected star in image coordinates.
#[derive(Debug, Clone)]
pub struct Star {
    pub x: f64,
    pub y: f64,
    pub snr: f64,
    pub hfd: f64,
}

/// Ordered list of detected stars.
#[derive(Debug, Clone, Default)]
pub struct StarList(pub Vec<Star>);

impl StarList {
    pub fn len(&self) -> usize {
        self.0.len()
    }
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

/// A 4-star quad pattern with distance ratios.
#[derive(Debug, Clone)]
pub struct Quad {
    /// Largest inter-star distance (absolute, in pixels or standard coords).
    pub d1: f64,
    /// d2/d1 through d6/d1: the 5 normalized distance ratios.
    pub ratios: [f64; 5],
    pub center_x: f64,
    pub center_y: f64,
    /// Direction angle (radians, 0..π) of the longest pair, used by vote_filter.
    pub d1_angle: f64,
}

/// List of quads.
#[derive(Debug, Clone, Default)]
pub struct QuadList(pub Vec<Quad>);

impl QuadList {
    pub fn len(&self) -> usize {
        self.0.len()
    }
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

/// Matched image and catalogue positions, ready for `solve_plate_constants`.
///
/// The two vectors are parallel: element `i` of the first is the image position whose
/// counterpart is element `i` of the second. Which space the catalogue side is in
/// depends on the caller — the pipeline passes tangent-plane standard coordinates.
pub type PairedPositions = (Vec<(f64, f64)>, Vec<(f64, f64)>);

/// The 6 plate constants for the linear WCS transformation.
/// X_ref = a·x + b·y + c
/// Y_ref = d·x + e·y + f
/// Units depend on context: when called with equatorial_standard coords (cdelt=1),
/// a/b/d/e are in standard coordinate units per pixel and c/f are offsets.
#[derive(Debug, Clone)]
pub struct PlateConstants {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub d: f64,
    pub e: f64,
    pub f: f64,
}

/// The final WCS solution for an image.
#[derive(Debug, Clone)]
pub struct WcsSolution {
    pub ra0: f64,    // CRVAL1 in radians
    pub dec0: f64,   // CRVAL2 in radians
    pub crpix1: f64, // reference pixel X (1-based FITS convention)
    pub crpix2: f64, // reference pixel Y (1-based FITS convention)
    pub cd1_1: f64,  // CD matrix elements (degrees/pixel)
    pub cd1_2: f64,
    pub cd2_1: f64,
    pub cd2_2: f64,
    pub cdelt1: f64, // degrees/pixel
    pub cdelt2: f64,
    pub crota2: f64,       // rotation in degrees
    pub residual_rms: f64, // arcsec
    pub stars_matched: usize,
    /// The 6 plate constants in arcsec/pixel units (as returned by solve_plate_constants).
    pub plate: PlateConstants,
    /// Faintest catalog star magnitude used in the matching step.
    pub mag_limit: f64,
    /// Angular distance (degrees) from the search hint to the winning spiral position.
    pub search_dist_deg: f64,
    /// Angular distance to each tried spiral position that yielded a catalog read (degrees).
    /// Used to generate the "Xd,Yd,..." progress line in ASTAP-compatible output.
    pub step_distances: Vec<f64>,
    /// Total quad matches found before scale-outlier filtering (the M in "N of M quads").
    pub raw_matches: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn buf(data: Vec<f32>, w: usize, h: usize) -> ImageBuffer {
        ImageBuffer {
            data,
            width: w,
            height: h,
        }
    }

    #[test]
    fn normalize_leaves_adu_like_data_untouched() {
        // A 16-bit-ish frame: background 1000, a star at 30000.
        let mut d = vec![1000.0f32; 64];
        d[10] = 30000.0;
        let mut img = buf(d.clone(), 8, 8);
        assert!(img.normalize_for_detection().is_none());
        assert_eq!(img.data, d, "ADU-like data must not be rescaled");
    }

    #[test]
    fn normalize_expands_small_range_float_data() {
        // SDSS-like nanomaggies: everything would land in histogram bins 0..4.
        let mut d = vec![0.0f32; 1000];
        for (i, v) in d.iter_mut().enumerate() {
            *v = -0.15 + (i % 7) as f32 * 0.01; // noise around a small background
        }
        d[500] = 4.5; // a star
        let mut img = buf(d, 100, 10);
        let (scale, _offset) = img.normalize_for_detection().expect("should rescale");
        assert!(scale > 1000.0, "scale={scale} should open the range up");
        let lo = img.data.iter().cloned().fold(f32::INFINITY, f32::min);
        let hi = img.data.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        assert!(lo >= 0.0, "minimum {lo} should be non-negative");
        assert!(hi - lo > 4096.0, "range {} should be resolvable", hi - lo);
    }

    #[test]
    fn normalize_replaces_non_finite_with_the_finite_minimum() {
        let mut img = buf(
            vec![
                f32::NAN,
                5000.0,
                1000.0,
                f32::INFINITY,
                20000.0,
                f32::NEG_INFINITY,
            ],
            3,
            2,
        );
        img.normalize_for_detection();
        assert!(
            img.data.iter().all(|v| v.is_finite()),
            "no NaN/inf may survive"
        );
        assert_eq!(img.data[0], 1000.0, "NaN becomes the finite minimum");
        assert_eq!(img.data[3], 1000.0);
        assert_eq!(img.data[5], 1000.0);
    }

    #[test]
    fn normalize_handles_an_all_nan_frame() {
        let mut img = buf(vec![f32::NAN; 16], 4, 4);
        assert!(img.normalize_for_detection().is_none());
        assert!(img.data.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn normalize_handles_a_flat_frame() {
        let mut img = buf(vec![7.0f32; 16], 4, 4);
        img.normalize_for_detection();
        assert!(img.data.iter().all(|v| v.is_finite()));
    }
}
