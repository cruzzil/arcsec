//! Plain data types shared across the pipeline.

/// Row-major greyscale pixel buffer. Access: `data[y * width + x]`.
///
/// `data.len()` must equal `width * height`; the detection code indexes on that
/// assumption.
#[derive(Debug, Clone)]
pub struct ImageBuffer {
    /// Pixel values, row by row from the top-left.
    pub data: Vec<f32>,
    /// Columns.
    pub width: usize,
    /// Rows.
    pub height: usize,
}

impl ImageBuffer {
    /// A zero-filled image of the given size.
    #[must_use]
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            data: vec![0.0; width * height],
            width,
            height,
        }
    }

    /// Pixel at `(x, y)`.
    ///
    /// # Panics
    ///
    /// Panics if the coordinates are outside the buffer.
    #[inline]
    #[must_use]
    pub fn get(&self, x: usize, y: usize) -> f32 {
        self.data[y * self.width + x]
    }

    /// Pixel at `(x, y)`, or `None` if the coordinates are outside the image.
    #[inline]
    #[must_use]
    pub fn get_checked(&self, x: i32, y: i32) -> Option<f32> {
        if x >= 0 && y >= 0 && (x as usize) < self.width && (y as usize) < self.height {
            Some(self.data[y as usize * self.width + x as usize])
        } else {
            None
        }
    }

    /// Even out the four pixel phases of a Bayer matrix (ASTAP's check-pattern
    /// filter, `-check y`).
    ///
    /// A raw one-shot-colour frame has a colour filter over each pixel in a
    /// repeating 2×2 pattern, so its pixels alternate in brightness, and the
    /// detector sees the checkerboard rather than the stars. This scales each of the
    /// four phases (even/odd column × even/odd row) so its mean over the central
    /// quarter of the frame matches the brightest phase's, then rounds to whole
    /// numbers, as ASTAP does. Only meaningful on an unbinned raw mosaic.
    ///
    /// Returns `false`, leaving the image untouched, if a phase has no positive
    /// mean to scale by (an image under 2×2 pixels, or a blank one).
    pub fn check_pattern_filter(&mut self) -> bool {
        let (w, h) = (self.width, self.height);
        let mut sum = [0.0f64; 4];
        let mut count = [0u64; 4];
        let phase = |x: usize, y: usize| (x & 1) + 2 * (y & 1);
        for y in h / 4..=(h * 3 / 4).min(h.saturating_sub(1)) {
            for x in w / 4..=(w * 3 / 4).min(w.saturating_sub(1)) {
                sum[phase(x, y)] += f64::from(self.data[y * w + x]);
                count[phase(x, y)] += 1;
            }
        }
        let mut mean = [0.0f64; 4];
        for k in 0..4 {
            if count[k] == 0 {
                return false;
            }
            mean[k] = sum[k] / count[k] as f64;
            if !(mean[k] > 0.0 && mean[k].is_finite()) {
                return false;
            }
        }
        let max = mean.iter().copied().fold(f64::MIN, f64::max);
        let factor = mean.map(|m| max / m);
        for y in 0..h {
            for x in 0..w {
                let f = factor[phase(x, y)];
                if f != 1.0 {
                    let v = &mut self.data[y * w + x];
                    *v = (f64::from(*v) * f).round_ties_even() as f32;
                }
            }
        }
        true
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
        sample.sort_unstable_by(f32::total_cmp);
        let p999 = sample[((sample.len() as f64 * 0.999) as usize).min(sample.len() - 1)];

        let span = p999 - lo;
        if span <= 0.0 {
            return None;
        }
        let scale = TARGET_P999 / span;
        let offset = FLOOR - lo * scale;
        for v in &mut self.data {
            *v = *v * scale + offset;
        }
        Some((scale, offset))
    }

    /// Downsample by averaging N×N pixel blocks. Returns a clone if factor ≤ 1.
    ///
    /// Partial blocks at the right and bottom edges are dropped.
    #[must_use]
    pub fn bin_image(&self, factor: usize) -> Self {
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
        Self {
            data,
            width: new_w,
            height: new_h,
        }
    }
}

/// A detected star in image coordinates.
#[derive(Debug, Clone)]
pub struct Star {
    /// Sub-pixel column of the centroid (0-based).
    pub x: f64,
    /// Sub-pixel row of the centroid (0-based).
    pub y: f64,
    /// Signal-to-noise ratio of the star's aperture flux.
    pub snr: f64,
    /// Half-flux diameter in pixels.
    pub hfd: f64,
}

/// Ordered list of detected stars.
#[derive(Debug, Clone, Default)]
pub struct StarList(pub Vec<Star>);

impl StarList {
    /// Number of stars.
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }
    /// Whether the list is empty.
    #[must_use]
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
    /// Mean x of the four stars.
    pub center_x: f64,
    /// Mean y of the four stars.
    pub center_y: f64,
    /// Direction angle (radians, 0..π) of the longest pair, used by `vote_filter`.
    pub d1_angle: f64,
}

/// List of quads.
#[derive(Debug, Clone, Default)]
pub struct QuadList(pub Vec<Quad>);

impl QuadList {
    /// Number of quads.
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }
    /// Whether the list is empty.
    #[must_use]
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
///
/// `X_ref = a·x + b·y + c` and `Y_ref = d·x + e·y + f`.
///
/// Units depend on context: when called with `equatorial_standard` coords (cdelt=1),
/// a/b/d/e are in standard coordinate units (arcsec) per pixel and c/f are offsets.
#[derive(Debug, Clone)]
pub struct PlateConstants {
    /// `∂X_ref/∂x`.
    pub a: f64,
    /// `∂X_ref/∂y`.
    pub b: f64,
    /// `X_ref` at pixel (0, 0).
    pub c: f64,
    /// `∂Y_ref/∂x`.
    pub d: f64,
    /// `∂Y_ref/∂y`.
    pub e: f64,
    /// `Y_ref` at pixel (0, 0).
    pub f: f64,
}

/// The final WCS solution for an image.
#[derive(Debug, Clone)]
pub struct WcsSolution {
    /// CRVAL1: RA of the reference pixel, radians.
    pub ra0: f64,
    /// CRVAL2: Dec of the reference pixel, radians.
    pub dec0: f64,
    /// CRPIX1: reference pixel X (1-based FITS convention; the image centre).
    pub crpix1: f64,
    /// CRPIX2: reference pixel Y (1-based FITS convention; the image centre).
    pub crpix2: f64,
    /// `CD1_1`, degrees/pixel.
    pub cd1_1: f64,
    /// `CD1_2`, degrees/pixel.
    pub cd1_2: f64,
    /// `CD2_1`, degrees/pixel.
    pub cd2_1: f64,
    /// `CD2_2`, degrees/pixel.
    pub cd2_2: f64,
    /// CDELT1, degrees/pixel, signed as `astap_cli` writes it: negative for an image with
    /// the sky's usual handedness (east to the left with north up), positive for a
    /// mirrored one. See [`crate::wcs::output::old_style_wcs`].
    pub cdelt1: f64,
    /// CDELT2, degrees/pixel. Always positive.
    pub cdelt2: f64,
    /// CROTA2, degrees: the rotation of the image's +Y axis from north, as `astap_cli`
    /// reports it. [`WcsSolution::crota1`] gives the +X axis's.
    pub crota2: f64,
    /// RMS residual of the verified star matches, arcsec.
    pub residual_rms: f64,
    /// Number of individual stars matched during verification.
    pub stars_matched: usize,
    /// The 6 plate constants in arcsec/pixel units (as returned by `solve_plate_constants`).
    ///
    /// These are in the pixel units of the image that was solved, i.e. *after*
    /// binning; the CRPIX/CD/CDELT fields have been scaled back to the original image.
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
    /// The star pairs the final fit was verified on: each detected star with the
    /// catalogue star it was identified as. What [`crate::wcs::sip::fit_sip`] fits.
    ///
    /// For a field the solver found distorted, these are the pairs its distortion
    /// model verified, over the whole frame: they agree with that model to the
    /// verification radius, not with the linear CD matrix, which is the closest
    /// linear approximation to it.
    pub matched_stars: Vec<MatchedStar>,
    /// SIP distortion polynomials on top of the linear solution, if fitted.
    ///
    /// [`crate::pipeline::solve_image`] leaves this `None`; add it with
    /// [`crate::wcs::sip::fit_sip`].
    pub sip: Option<crate::wcs::sip::Sip>,
}

impl WcsSolution {
    /// CROTA1, degrees: the rotation of the image's +X axis, as `astap_cli` reports it.
    /// It differs from [`WcsSolution::crota2`] only when the plate is slightly skewed.
    #[must_use]
    pub fn crota1(&self) -> f64 {
        crate::wcs::output::old_style_wcs(self.cd1_1, self.cd1_2, self.cd2_1, self.cd2_2).3
    }
}

/// A detected star paired with the catalogue star it was identified as.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MatchedStar {
    /// Column of the detected star: 1-based FITS pixels on the unbinned image.
    pub x: f64,
    /// Row of the detected star: 1-based FITS pixels on the unbinned image.
    pub y: f64,
    /// Catalogue right ascension, radians.
    pub ra: f64,
    /// Catalogue declination, radians.
    pub dec: f64,
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
    fn check_pattern_filter_evens_out_a_bayer_mosaic() {
        // An RGGB-like mosaic: R 1000, G 2000, B 500, on a 40×30 frame.
        let (w, h) = (40, 30);
        let level = |x: usize, y: usize| match (x & 1, y & 1) {
            (0, 0) => 1000.0,
            (1, 1) => 500.0,
            _ => 2000.0,
        };
        let mut img = buf((0..w * h).map(|i| level(i % w, i / w)).collect(), w, h);
        // A "star" on a red pixel, outside the central quarter the means come from,
        // is scaled with its phase.
        img.data[2 * w + 2] = 1500.0;
        assert!(img.check_pattern_filter());
        for (i, &v) in img.data.iter().enumerate() {
            if i == 2 * w + 2 {
                assert_eq!(v, 3000.0);
            } else {
                assert_eq!(v, 2000.0, "pixel {i}");
            }
        }
        // Nothing to scale by: left alone.
        let mut blank = buf(vec![0.0; 16], 4, 4);
        assert!(!blank.check_pattern_filter());
        assert!(!buf(vec![5.0], 1, 1).check_pattern_filter());
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
        let lo = img.data.iter().copied().fold(f32::INFINITY, f32::min);
        let hi = img.data.iter().copied().fold(f32::NEG_INFINITY, f32::max);
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

    /// A frame big enough for the threaded min/max scan (over 2¹⁸ pixels per
    /// chunk) must give exactly what the serial scan gives: the global minimum
    /// found in one chunk replaces non-finite pixels found in another.
    #[test]
    fn normalize_large_frames_match_the_serial_rules() {
        let (w, h) = (1024usize, 1024usize);
        let mut data: Vec<f32> = (0..w * h).map(|i| 0.5 + (i % 997) as f32 * 1e-3).collect();
        data[5] = -0.25; // the minimum, in the first chunk
        data[w * h - 3] = f32::NAN; // non-finite, in the last chunk
        data[w * h / 2 + 7] = f32::NEG_INFINITY;
        let mut img = buf(data, w, h);
        let (scale, offset) = img.normalize_for_detection().expect("rescaled");
        assert!(img.data.iter().all(|v| v.is_finite()));
        // Non-finite pixels became the minimum, which maps to the floor of 100.
        let floor = -0.25 * scale + offset;
        assert!((floor - 100.0).abs() < 1e-2, "floor {floor}");
        assert_eq!(img.data[w * h - 3], img.data[5]);
        assert_eq!(img.data[w * h / 2 + 7], img.data[5]);
        // And a large ADU-like frame is left alone.
        let mut adu = buf((0..w * h).map(|i| (i % 60_000) as f32).collect(), w, h);
        assert!(adu.normalize_for_detection().is_none());
        assert_eq!(adu.data[59_999], 59_999.0);
    }

    #[test]
    fn get_checked_rejects_everything_outside_the_frame() {
        let img = buf((0..12).map(|v| v as f32).collect(), 4, 3);
        assert_eq!(img.get_checked(0, 0), Some(0.0));
        assert_eq!(img.get_checked(3, 2), Some(11.0));
        assert_eq!(img.get(1, 2), 9.0);
        for (x, y) in [(-1, 0), (0, -1), (4, 0), (0, 3), (i32::MIN, i32::MAX)] {
            assert_eq!(img.get_checked(x, y), None, "({x}, {y})");
        }
    }

    #[test]
    fn bin_image_averages_blocks_and_drops_partial_edges() {
        // 5×3, values x + 10y.
        let img = buf(
            (0..15)
                .map(|i| (i % 5) as f32 + 10.0 * (i / 5) as f32)
                .collect(),
            5,
            3,
        );
        let b = img.bin_image(2);
        assert_eq!((b.width, b.height), (2, 1));
        // Block (0,0): 0, 1, 10, 11 → 5.5; block (1,0): 2, 3, 12, 13 → 7.5.
        assert_eq!(b.data, vec![5.5, 7.5]);
        // Factor 1 (and 0) is a copy.
        assert_eq!(img.bin_image(1).data, img.data);
        assert_eq!(img.bin_image(0).data, img.data);
        // Binning conserves the mean over whole blocks.
        let big = buf((0..64 * 48).map(|i| (i * 7 % 101) as f32).collect(), 64, 48);
        let b4 = big.bin_image(4);
        let mean = |v: &[f32]| v.iter().map(|&x| f64::from(x)).sum::<f64>() / v.len() as f64;
        assert!((mean(&b4.data) - mean(&big.data)).abs() < 1e-3);
    }

    #[test]
    fn list_helpers() {
        assert!(StarList::default().is_empty());
        assert!(QuadList::default().is_empty());
        let s = StarList(vec![Star {
            x: 0.0,
            y: 0.0,
            snr: 1.0,
            hfd: 1.0,
        }]);
        assert_eq!(s.len(), 1);
        assert!(!s.is_empty());
    }
}
