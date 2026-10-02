//! Star detection: HFD measurement and multi-pass star finder.

use std::collections::HashMap;

use crate::detection::background::{
    Background, Region, get_background, sigma_clipped_mean_from_histogram,
};
use crate::types::{ImageBuffer, Star, StarList};

/// The four numbers that decide what counts as a star in one scan pass.
///
/// They are set together and passed together through every layer of the detection
/// cascade, so carrying them as one value keeps those signatures short.
#[derive(Debug, Clone, Copy)]
struct Thresholds {
    /// Sky level subtracted before the threshold test.
    background: f64,
    /// Background noise σ, used by the hot-pixel test.
    noise: f64,
    /// Height above background a pixel must exceed to seed a candidate.
    detection_level: f64,
    /// Candidates with an HFD at or below this are rejected (pixels).
    hfd_min: f64,
}

const ANNULUS_RS: i32 = 14; // search radius passed to HFD (annulus outer = rs+1)
const RASTER_STEPS: usize = 12; // sub-grid divisions for the fallback retry

/// `round(sqrt(n))` for a non-negative integer, without floating point.
///
/// This sits in the innermost loop of `measure_star`, once per pixel of the
/// aperture box, and the libm `round` call alone was 6% of total runtime.
///
/// `round(sqrt(n)) == k` exactly when `(k-0.5)^2 <= n < (k+0.5)^2`, i.e. when
/// `k^2 - k + 1 <= n <= k^2 + k` over the integers.
///
/// The aperture loop calls this once per pixel with `n = i*i + j*j` bounded by
/// `2 * (ANNULUS_RS + 2)^2 = 512`, so that range is a compile-time table and the
/// integer sqrt (which showed at 3.6% of a profile) never runs there.
#[inline(always)]
fn round_sqrt(n: i32) -> usize {
    if (n as usize) < ROUND_SQRT_TABLE.len() {
        return ROUND_SQRT_TABLE[n as usize] as usize;
    }
    let r = (n as u32).isqrt();
    if n as u32 > r * r + r {
        (r + 1) as usize
    } else {
        r as usize
    }
}

// Covers `measure_large`'s box too: 2 * (LARGE_RS + 1)^2 = 2178.
const ROUND_SQRT_N: usize = 2304;
static ROUND_SQRT_TABLE: [u8; ROUND_SQRT_N] = build_round_sqrt_table();

const fn build_round_sqrt_table() -> [u8; ROUND_SQRT_N] {
    let mut t = [0u8; ROUND_SQRT_N];
    let mut n = 0usize;
    while n < ROUND_SQRT_N {
        let mut r = 0usize;
        while (r + 1) * (r + 1) <= n {
            r += 1;
        }
        t[n] = if n > r * r + r {
            (r + 1) as u8
        } else {
            r as u8
        };
        n += 1;
    }
    t
}

/// Compute HFD, SNR and sub-pixel centroid for a candidate star.
///
/// Returns `None` if the pixel region is too close to the border, or the candidate
/// fails star quality checks (not boxed, single hot pixel, too large).
#[must_use]
pub fn measure_star(img: &ImageBuffer, x1: i32, y1: i32) -> Option<Star> {
    measure::<false>(img, x1, y1).ok().map(|(star, _)| star)
}

/// As [`measure_star`], but with `astap_cli`'s current test for a star disc, and
/// also returning the star's flux: the background-subtracted sum over the
/// measuring aperture, in the image's pixel units.
///
/// The one difference in which stars pass is the disc test (a star whose
/// illuminated pixels fill too little of its aperture is taken for a blend):
/// `astap_cli` compares against 35% of `(2r - 2)²`, where the solver's detection
/// uses 35% of `(2r)²`. The solver keeps its stricter test, which the benchmark
/// corpus is calibrated on; the analysis behind `--analyse` and `--extract`
/// reports stars, and should report the ones ASTAP does.
#[must_use]
#[inline]
pub fn measure_star_with_flux(img: &ImageBuffer, x1: i32, y1: i32) -> Option<(Star, f64)> {
    measure::<true>(img, x1, y1).ok()
}

/// [`measure`] as the detection scan calls it, kept out of line: inlined, it
/// swells the scan's per-pixel loop, which every pixel of the frame runs through,
/// and on a frame that solves in 0.1 s the scan alone took 10 ms longer.
#[inline(never)]
fn measure_for_scan(img: &ImageBuffer, x1: i32, y1: i32) -> Result<(Star, f64), Reject> {
    measure::<false>(img, x1, y1)
}

/// Why [`measure`] turned a candidate down.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Reject {
    /// Too close to the frame edge, too faint, a single hot pixel, or not finite.
    NotAStar,
    /// The star's isophote reaches the edge of the measuring box.
    TooLarge,
    /// The illuminated pixels fill too little of the aperture.
    NotADisc,
}

/// The measurement behind [`measure_star`] and [`measure_star_with_flux`].
/// `CLI_DISC` selects `astap_cli`'s disc test; see [`measure_star_with_flux`].
#[inline(always)]
fn measure<const CLI_DISC: bool>(
    img: &ImageBuffer,
    x1: i32,
    y1: i32,
) -> Result<(Star, f64), Reject> {
    /// Annulus buffer size. The annulus is a fixed size (rs = `ANNULUS_RS`), about
    /// 91 pixels, so it lives on the stack: this runs once per candidate and a
    /// crowded field has tens of thousands of them.
    const BG_CAP: usize = 160;

    let width = img.width as i32;
    let height = img.height as i32;
    let mut rs = ANNULUS_RS;

    let r2 = rs + 1;
    if x1 - r2 <= 0 || x1 + r2 >= width - 1 || y1 - r2 <= 0 || y1 + r2 >= height - 1 {
        return Err(Reject::NotAStar);
    }

    // --- Annulus background ---
    let r1_sq = rs * rs;
    let r2_sq = r2 * r2;
    let mut bg_buf = [0.0f64; BG_CAP];
    let mut bg_len = 0usize;
    // Walk the annulus by row-runs instead of testing every pixel of the enclosing
    // box: the box is 31x31 = 961 pixels and only ~91 are in the annulus, so the
    // naive form spent 90% of its work rejecting. For a row `j` the condition
    // `r1_sq < i*i + j*j <= r2_sq` is just `min_abs <= |i| <= hi`, i.e. two
    // contiguous runs. Same pixels in the same order, so `bg_buf` is bit-identical.
    for j in -r2..=r2 {
        let jj = j * j;
        let hi2 = r2_sq - jj;
        if hi2 < 0 {
            continue;
        }
        let hi = (hi2 as u32).isqrt() as i32;
        let lo2 = r1_sq - jj;
        // Smallest |i| with i*i > lo2.
        let min_abs = if lo2 < 0 {
            0
        } else {
            (lo2 as u32).isqrt() as i32 + 1
        };
        if min_abs > hi {
            continue;
        }
        let row = (y1 + j) as usize * img.width;
        let push = |i: i32, bg_buf: &mut [f64; BG_CAP], bg_len: &mut usize| {
            if *bg_len < BG_CAP {
                bg_buf[*bg_len] = img.data[row + (x1 + i) as usize] as f64;
                *bg_len += 1;
            }
        };
        if min_abs == 0 {
            for i in -hi..=hi {
                push(i, &mut bg_buf, &mut bg_len);
            }
        } else {
            for i in -hi..=-min_abs {
                push(i, &mut bg_buf, &mut bg_len);
            }
            for i in min_abs..=hi {
                push(i, &mut bg_buf, &mut bg_len);
            }
        }
    }
    let bg_pixels = &mut bg_buf[..bg_len];
    let star_bg = median_f64(bg_pixels);
    for v in bg_pixels.iter_mut() {
        *v = (*v - star_bg).abs();
    }
    let mad = median_f64(bg_pixels);
    let sd_bg = (mad * 1.4826).max(1.0);

    // --- Iterative centroid with shrinking aperture ---
    let mut xc = x1 as f64;
    let mut yc = y1 as f64;

    // Shrink the aperture until the star is boxed or rs <= 1.
    //
    // The original re-scanned the whole (2rs+1)^2 box on every shrink step. The
    // per-pixel test does not depend on rs, so the box of half-width r is exactly
    // the union of Chebyshev rings 0..r: accumulate each ring once from a single
    // scan of the widest box, and the cascade becomes prefix sums over at most 15
    // rings. A star that shrinks 14 -> 6 used to read 841+625+441+289+169 = 2365
    // pixels; it now reads 841.
    let cut = 3.0 * sd_bg;
    let mut ring_val = [0.0f64; ANNULUS_RS as usize + 1];
    let mut ring_vx = [0.0f64; ANNULUS_RS as usize + 1];
    let mut ring_vy = [0.0f64; ANNULUS_RS as usize + 1];
    let mut ring_n = [0usize; ANNULUS_RS as usize + 1];

    for j in -rs..=rs {
        let row = (y1 + j) as usize * img.width;
        let aj = j.unsigned_abs() as usize;
        for i in -rs..=rs {
            let val = img.data[row + (x1 + i) as usize] as f64 - star_bg;
            if val > cut {
                let k = aj.max(i.unsigned_abs() as usize);
                ring_val[k] += val;
                ring_vx[k] += val * i as f64;
                ring_vy[k] += val * j as f64;
                ring_n[k] += 1;
            }
        }
    }

    let centroid_ok = loop {
        let r = rs as usize;
        let mut sum_val = 0.0f64;
        let mut sum_vx = 0.0f64;
        let mut sum_vy = 0.0f64;
        let mut sig_count = 0usize;
        for k in 0..=r {
            sum_val += ring_val[k];
            sum_vx += ring_vx[k];
            sum_vy += ring_vy[k];
            sig_count += ring_n[k];
        }

        if sum_val <= 12.0 * sd_bg {
            break false; // too noisy
        }

        let cx = x1 as f64 + sum_vx / sum_val;
        let cy = y1 as f64 + sum_vy / sum_val;

        let rs_f = rs as f64;
        if cx - rs_f < 0.0
            || cx + rs_f > width as f64 - 1.0
            || cy - rs_f < 0.0
            || cy + rs_f > height as f64 - 1.0
        {
            break false;
        }

        xc = cx;
        yc = cy;

        let side = (2 * rs as usize + 1).pow(2);
        let boxed = sig_count >= 2 * side / 9;
        if boxed {
            break true;
        }
        if sig_count <= 1 {
            break false; // single hot pixel
        }
        if rs > 4 {
            rs -= 2;
        } else if rs > 1 {
            rs -= 1;
        } else {
            break true;
        }
    };

    if !centroid_ok {
        return Err(Reject::NotAStar);
    }

    rs += 2; // extra margin around the star

    // --- Distance histogram for r_aperture detection ---
    //
    // The two loops below sample at (xc + i, yc + j) for integer i, j, so the
    // bilinear weights are the *same* for every sample - only the integer base
    // moves. Hoisting them makes each sample four loads and four multiply-adds with
    // constant weights, instead of recomputing the truncation, the fractional parts
    // and the bounds test per pixel. The valid i/j range is likewise computed once
    // rather than being rejected per pixel inside value_subpixel.
    let iw = img.width;
    let xt0 = xc as i32;
    let yt0 = yc as i32;
    let xf = xc - xt0 as f64;
    let yf = yc - yt0 as f64;
    let (w00, w01) = ((1.0 - xf) * (1.0 - yf), xf * (1.0 - yf));
    let (w10, w11) = ((1.0 - xf) * yf, xf * yf);
    // value_subpixel accepted xt in (0, width-2) and yt in (0, height-2).
    let i_lo_bound = 1 - xt0;
    let i_hi_bound = width - 3 - xt0;
    let j_lo_bound = 1 - yt0;
    let j_hi_bound = height - 3 - yt0;
    let sample = |i: i32, j: i32| -> f64 {
        let base = (yt0 + j) as usize * iw + (xt0 + i) as usize;
        let d = &img.data;
        w00 * d[base] as f64
            + w01 * d[base + 1] as f64
            + w10 * d[base + iw] as f64
            + w11 * d[base + iw + 1] as f64
    };

    let rs_clamped = rs.min(50) as usize;
    let mut dist_hist_buf = [0i32; 51];
    let dist_hist = &mut dist_hist_buf[..=rs_clamped];

    for j in (-rs).max(j_lo_bound)..=rs.min(j_hi_bound) {
        let jj = j * j;
        for i in (-rs).max(i_lo_bound)..=rs.min(i_hi_bound) {
            let d = round_sqrt(i * i + jj);
            if d <= rs_clamped {
                let val = sample(i, j) - star_bg;
                if val > cut {
                    dist_hist[d] += 1;
                }
            }
        }
    }

    // Walk outward to find r_aperture (where histogram drops to <10% of peak)
    let mut r_aperture = 0usize;
    let mut dist_top = 0i32;
    let mut hist_started = false;
    let mut illuminated = 0i32;
    loop {
        illuminated += dist_hist[r_aperture];
        if dist_hist[r_aperture] > 0 {
            hist_started = true;
        }
        if dist_hist[r_aperture] > dist_top {
            dist_top = dist_hist[r_aperture];
        }
        if r_aperture >= rs_clamped
            || (hist_started && dist_hist[r_aperture] <= (0.1 * dist_top as f64) as i32)
        {
            break;
        }
        r_aperture += 1;
    }

    if r_aperture >= rs_clamped {
        return Err(Reject::TooLarge); // star is larger than detection box
    }
    let disc_side = if CLI_DISC {
        2 * r_aperture - 2
    } else {
        2 * r_aperture
    };
    if r_aperture > 2 && (illuminated as f64) < 0.35 * disc_side.pow(2) as f64 {
        return Err(Reject::NotADisc); // not a disk — likely overlapping stars
    }

    // --- HFD and SNR calculation ---
    let mut sum_val = 0.0f64;
    let mut sum_val_r = 0.0f64;

    let ra = r_aperture as i32;
    for j in (-ra).max(j_lo_bound)..=ra.min(j_hi_bound) {
        let jj = j * j;
        for i in (-ra).max(i_lo_bound)..=ra.min(i_hi_bound) {
            let val = sample(i, j) - star_bg;
            let r = ((i * i + jj) as f64).sqrt();
            sum_val += val;
            sum_val_r += val * r;
        }
    }

    let flux = sum_val.max(1e-5);
    let hfd = (2.0 * sum_val_r / flux).max(0.7);
    let snr =
        flux / (flux + (r_aperture as f64).powi(2) * core::f64::consts::PI * sd_bg.powi(2)).sqrt();

    // A NaN anywhere in the pixel neighbourhood (chip gaps, masked columns, the
    // edge of a reprojected survey cutout) propagates into the centroid, HFD and
    // SNR. Such a candidate is not a star, and letting it through used to make the
    // brightness sort panic on a comparator that is not a total order.
    if !(xc.is_finite() && yc.is_finite() && snr.is_finite() && hfd.is_finite()) {
        return Err(Reject::NotAStar);
    }

    Ok((
        Star {
            x: xc,
            y: yc,
            snr,
            hfd,
        },
        flux,
    ))
}

/// A detection and the radius (pixels) of the area marked out around it, which
/// later candidates inside are not measured again.
struct Found {
    star: Star,
    mark: i32,
    /// Measured by [`measure_large`].
    large: bool,
}

/// Half-width of the box in which a star too large for the ordinary measurement
/// is measured again ([`measure_large`]).
const LARGE_RS: i32 = 32;
/// [`measure_large`]'s isophote: this fraction of the star's peak above the local
/// background, or 3σ if that is higher.
const LARGE_PEAK_FRACTION: f64 = 0.05;
/// Largest ratio of the eigenvalues of the second moments (the squared axis
/// ratio) [`measure_large`] accepts.
const LARGE_MAX_ELONGATION: f64 = 2.5;

/// Measure a bright star that [`measure`] turned down as too large or not a disc.
///
/// The ordinary measurement takes the star's extent at 3σ of the local background
/// in a box 14 pixels from the seed. A bright star on a photographic plate is a
/// saturated disc 10-40 pixels across, and a bright star in an undersampled
/// wide-field camera (TESS, 21″/px) has faint wings that reach past the box at 3σ
/// although its core is 3 pixels wide; both are refused, and those are exactly the
/// stars the catalogue's brightest are. Here the star is re-centred on its
/// brightest pixel, measured in a box of [`LARGE_RS`], and its extent taken at
/// [`LARGE_PEAK_FRACTION`] of its peak: the core, not the wings.
///
/// A source still too large, not a disc, or elongated (a trail, a galaxy, a
/// nebula knot) is refused.
///
/// Every pixel of a large source above the threshold is a seed, so the area each
/// measurement covered is remembered in `visited` (cells of [`VISITED_CELL`]
/// pixels) and later seeds there are passed over: a source is measured once.
fn measure_large(
    img: &ImageBuffer,
    x1: i32,
    y1: i32,
    detect_abs: f64,
    visited: &mut Visited,
) -> Option<(Star, i32)> {
    const SEARCH: i32 = LARGE_RS / 2;
    if visited.get(x1, y1) {
        return None;
    }
    let (w, h) = (img.width as i32, img.height as i32);
    // The seed is the first pixel above the threshold in raster order, usually
    // the top edge of the star: re-centre on the brightest pixel near it.
    if x1 - SEARCH < 0 || x1 + SEARCH >= w || y1 - SEARCH < 0 || y1 + SEARCH >= h {
        return None;
    }
    let (mut px, mut py, mut peak) = (x1, y1, f32::NEG_INFINITY);
    for y in y1 - SEARCH..=y1 + SEARCH {
        let row = &img.data[y as usize * img.width..][..img.width];
        for x in x1 - SEARCH..=x1 + SEARCH {
            if row[x as usize] > peak {
                (px, py, peak) = (x, y, row[x as usize]);
            }
        }
    }
    if visited.get(px, py) {
        visited.set(x1, y1);
        return None;
    }
    let (found, radius) = measure_large_at(img, px, py, f64::from(peak), detect_abs);
    // Mark what this measurement covered: the star's own marked area, or, if it
    // was refused, its aperture (the whole box if it was too large for it).
    visited.set_square(px, py, radius.max(VISITED_CELL));
    visited.set(x1, y1);
    found
}

/// Side, in pixels, of the cells [`measure_large`] remembers as measured.
const VISITED_CELL: i32 = 4;

/// One bit per [`VISITED_CELL`]-pixel cell of the image: where [`measure_large`]
/// has already been run.
struct Visited {
    bits: Vec<u64>,
    nx: i32,
    ny: i32,
}

impl Visited {
    fn new(img: &ImageBuffer) -> Self {
        let nx = (img.width as i32).div_euclid(VISITED_CELL) + 1;
        let ny = (img.height as i32).div_euclid(VISITED_CELL) + 1;
        Self {
            bits: vec![0; (nx as usize * ny as usize).div_ceil(64)],
            nx,
            ny,
        }
    }

    /// The bit of the cell holding pixel `(x, y)`, if it is in the image.
    fn index(&self, x: i32, y: i32) -> Option<usize> {
        let (cx, cy) = (x.div_euclid(VISITED_CELL), y.div_euclid(VISITED_CELL));
        (cx >= 0 && cy >= 0 && cx < self.nx && cy < self.ny).then(|| (cy * self.nx + cx) as usize)
    }

    fn get(&self, x: i32, y: i32) -> bool {
        self.index(x, y)
            .is_some_and(|i| self.bits[i / 64] & (1 << (i % 64)) != 0)
    }

    fn set(&mut self, x: i32, y: i32) {
        if let Some(i) = self.index(x, y) {
            self.bits[i / 64] |= 1 << (i % 64);
        }
    }

    /// Every cell within `r` pixels (in x and y) of `(x, y)`.
    fn set_square(&mut self, x: i32, y: i32, r: i32) {
        for yy in ((y - r).max(0)..=y + r).step_by(VISITED_CELL as usize) {
            for xx in ((x - r).max(0)..=x + r).step_by(VISITED_CELL as usize) {
                self.set(xx, yy);
            }
            self.set(x + r, yy);
        }
        for xx in ((x - r).max(0)..=x + r).step_by(VISITED_CELL as usize) {
            self.set(xx, y + r);
        }
    }
}

/// [`measure_large`] about the brightest pixel `(px, py)`, value `peak`: the star
/// and its marked radius, if accepted, and the radius the measurement covered.
fn measure_large_at(
    img: &ImageBuffer,
    px: i32,
    py: i32,
    peak: f64,
    detect_abs: f64,
) -> (Option<(Star, i32)>, i32) {
    const R_OUT: i32 = LARGE_RS + 1;
    let (w, h) = (img.width as i32, img.height as i32);
    let at = |x: i32, y: i32| f64::from(img.data[y as usize * img.width + x as usize]);
    let refused = |r: i32| (None, r);
    if px - R_OUT <= 0 || px + R_OUT >= w - 1 || py - R_OUT <= 0 || py + R_OUT >= h - 1 {
        return refused(LARGE_RS / 2);
    }

    // Local background and noise from the annulus LARGE_RS < r <= LARGE_RS + 1.
    let mut ann: Vec<f64> = Vec::with_capacity(256);
    for j in -R_OUT..=R_OUT {
        for i in -R_OUT..=R_OUT {
            let d2 = i * i + j * j;
            if d2 > LARGE_RS * LARGE_RS && d2 <= R_OUT * R_OUT {
                ann.push(at(px + i, py + j));
            }
        }
    }
    let bg = median_f64(&mut ann);
    for v in &mut ann {
        *v = (*v - bg).abs();
    }
    let sd = (median_f64(&mut ann) * 1.4826).max(1.0);
    let height = peak - bg;
    if height.is_nan() || height <= 10.0 * sd {
        return refused(LARGE_RS / 2);
    }
    let cut = (3.0 * sd).max(LARGE_PEAK_FRACTION * height);

    // Extent: walk out in rings until the count above the cut falls to a tenth
    // of the fullest ring, as `measure` does.
    // Also, per ring, the pixels above the detection threshold, and all pixels:
    // the star's area is marked out to where its wings fall below the threshold,
    // so that the scan does not take them for stars of their own.
    let mut hist = [0u32; LARGE_RS as usize + 1];
    let mut above = [0u32; LARGE_RS as usize + 1];
    let mut ring = [0u32; LARGE_RS as usize + 1];
    for j in -LARGE_RS..=LARGE_RS {
        for i in -LARGE_RS..=LARGE_RS {
            let d = round_sqrt(i * i + j * j);
            if d <= LARGE_RS as usize {
                let v = at(px + i, py + j);
                ring[d] += 1;
                hist[d] += u32::from(v - bg > cut);
                above[d] += u32::from(v > detect_abs);
            }
        }
    }
    let (mut r_ap, mut top, mut illuminated) = (0usize, 0u32, 0u32);
    loop {
        illuminated += hist[r_ap];
        top = top.max(hist[r_ap]);
        if r_ap >= LARGE_RS as usize || hist[r_ap] * 10 <= top {
            break;
        }
        r_ap += 1;
    }
    let ra = r_ap as i32;
    if r_ap >= LARGE_RS as usize {
        return refused(LARGE_RS);
    }
    if r_ap > 2 && f64::from(illuminated) < 0.35 * ((2 * r_ap).pow(2)) as f64 {
        return refused(ra);
    }

    // Centroid and second moments of the pixels above the cut within the aperture,
    // flux, HFD and SNR over it. The brightest pixel of a saturated disc can be
    // anywhere on its plateau, and an aperture about it would cut the disc
    // unevenly: re-centre the aperture on the centroid until it settles.
    let (mut cx, mut cy) = (px, py);
    let mut moments = [0.0f64; 8];
    for _ in 0..4 {
        if cx - ra < 1 || cx + ra >= w - 1 || cy - ra < 1 || cy + ra >= h - 1 {
            return refused(ra);
        }
        let [
            mut sw,
            mut sx,
            mut sy,
            mut sxx,
            mut syy,
            mut sxy,
            mut flux,
            mut flux_r,
        ] = [0.0; 8];
        for j in -ra..=ra {
            for i in -ra..=ra {
                let v = at(cx + i, cy + j) - bg;
                flux += v;
                flux_r += v * f64::from(i * i + j * j).sqrt();
                if v > cut && i * i + j * j <= ra * ra {
                    let (fi, fj) = (f64::from(i), f64::from(j));
                    sw += v;
                    sx += v * fi;
                    sy += v * fj;
                    sxx += v * fi * fi;
                    syy += v * fj * fj;
                    sxy += v * fi * fj;
                }
            }
        }
        moments = [sw, sx, sy, sxx, syy, sxy, flux, flux_r];
        if sw.is_nan() || sw <= 0.0 {
            return refused(ra);
        }
        let (nx, ny) = (cx + (sx / sw).round() as i32, cy + (sy / sw).round() as i32);
        if (nx, ny) == (cx, cy) {
            break;
        }
        (cx, cy) = (nx, ny);
    }
    let [sw, sx, sy, sxx, syy, sxy, flux, flux_r] = moments;
    if !(sw > 0.0 && flux > 0.0) {
        return refused(ra);
    }
    let (mx, my) = (sx / sw, sy / sw);
    let (cxx, cyy, cxy) = (sxx / sw - mx * mx, syy / sw - my * my, sxy / sw - mx * my);
    let tr = cxx + cyy;
    let disc = ((cxx - cyy).powi(2) + 4.0 * cxy * cxy).sqrt();
    let (l1, l2) = (0.5 * (tr + disc), 0.5 * (tr - disc));
    if l2 <= 0.0 || l1 > LARGE_MAX_ELONGATION * l2 {
        return refused(ra);
    }
    let hfd = (2.0 * flux_r / flux).max(0.7);
    let snr = flux / (flux + (r_ap as f64).powi(2) * core::f64::consts::PI * sd * sd).sqrt();
    let star = Star {
        x: f64::from(cx) + mx,
        y: f64::from(cy) + my,
        snr,
        hfd,
    };
    // Marked out to the first ring less than half above the detection threshold,
    // and never further than the ordinary measurement would mark.
    let mark = (r_ap..=LARGE_RS as usize)
        .find(|&d| above[d] * 2 < ring[d])
        .unwrap_or(LARGE_RS as usize) as i32
        + 1;
    let mark = mark.min((3.0 * hfd).round() as i32);
    let ok = star.x.is_finite() && star.y.is_finite() && snr.is_finite() && hfd.is_finite();
    (ok.then_some((star, mark)), mark.max(ra))
}

/// Detect stars in an image using the ASTAP 4-retry strategy.
///
/// Retry passes (matching ASTAP exactly):
/// - Pass 4: `detection_level = star_level` (if > 30 × noise)
/// - Pass 3: `detection_level = star_level2` (if > 30 × noise)
/// - Pass 2: `detection_level = 30 × noise`
/// - Pass 1: 12×12 grid, local background recalc, `detection_level = 7 × noise`
///
/// Stars are marked in a mask after detection to avoid double-counting.
/// After collection, trims to the `max_stars` brightest by SNR.
#[must_use]
pub fn find_stars(img: &ImageBuffer, hfd_min: f64, max_stars: usize) -> StarList {
    let w = img.width;
    let h = img.height;
    let bg = get_background(img, max_stars);

    find_stars_with_background(img, &bg, hfd_min, max_stars, w, h).0
}

/// As [`find_stars`], with a pre-computed background.
///
/// `w` and `h` must equal `img.width` and `img.height`. Images smaller than 3×3
/// yield no stars.
///
/// Returns `(stars, raw_count)` where `raw_count` is the total found before
/// trimming to `max_stars` — used to emit ASTAP-style progress messages.
#[must_use]
pub fn find_stars_with_background(
    img: &ImageBuffer,
    bg: &Background,
    hfd_min: f64,
    max_stars: usize,
    w: usize,
    h: usize,
) -> (StarList, usize) {
    debug_assert_eq!((w, h), (img.width, img.height));
    let mut stars = detect_all(img, bg, hfd_min, max_stars);
    let raw_count = stars.len();

    // Trim to max_stars brightest by SNR
    if stars.len() > max_stars {
        stars.sort_by(|a, b| b.snr.total_cmp(&a.snr));
        stars.truncate(max_stars);
    }

    (StarList(stars), raw_count)
}

/// As [`find_stars_with_background`], and also the brightest `deep` of every star
/// the cascade found (by SNR, however many that is beyond `max_stars`).
///
/// The cascade stops at the first level that brings the count to `max_stars`, and
/// that level is scanned whole, so it usually finds more than `max_stars`; the
/// solver's catalogue-seeded fallback uses them.
#[must_use]
pub fn find_stars_and_deep(
    img: &ImageBuffer,
    bg: &Background,
    hfd_min: f64,
    max_stars: usize,
    deep: usize,
) -> (StarList, usize, StarList) {
    let mut stars = detect_all(img, bg, hfd_min, max_stars);
    let raw_count = stars.len();
    let mut more = stars.clone();
    more.sort_by(|a, b| b.snr.total_cmp(&a.snr));
    more.truncate(deep);
    if stars.len() > max_stars {
        stars.sort_by(|a, b| b.snr.total_cmp(&a.snr));
        stars.truncate(max_stars);
    }
    (StarList(stars), raw_count, StarList(more))
}

/// Every star the detection cascade finds, in the order found.
fn detect_all(img: &ImageBuffer, bg: &Background, hfd_min: f64, max_stars: usize) -> Vec<Star> {
    let (w, h) = (img.width, img.height);
    // The scan works on the frame minus a one-pixel border; anything smaller has
    // nothing to scan, and the inset region's bounds would underflow.
    if w < 3 || h < 3 || img.data.len() < w * h {
        return Vec::new();
    }
    let mut stars: Vec<Star> = Vec::with_capacity(max_stars + 1000);
    // img_sa: persistent star-area map. 1 = already detected, 0 = free.
    // Using a single marker prevents double-detection across retry passes.
    let mut img_sa = vec![0u8; w * h];

    let noise = bg.noise;
    let background = bg.mean;

    // Cascade through detection levels.
    // Each level runs if the previous level found too few stars.
    // Stars only `measure_large` accepts are bright stars the cascade has always
    // missed; they do not count towards the stars that end it, so the same levels
    // run, and find the same stars, as without them.
    let mut n_large = 0usize;
    let mut level = 4u8;
    while stars.len() - n_large < max_stars && level >= 1 {
        let mut pass_stars: Vec<Found> = Vec::new();

        match level {
            4 => {
                if bg.star_level > 30.0 * noise {
                    detect_pass(
                        img,
                        &mut img_sa,
                        Thresholds {
                            background,
                            noise,
                            detection_level: bg.star_level,
                            hfd_min,
                        },
                        Region::inset(img),
                        &mut pass_stars,
                    );
                }
                // If star_level too low or nothing found, fall through to level 3
            }
            3 => {
                if bg.star_level2 > 30.0 * noise {
                    detect_pass(
                        img,
                        &mut img_sa,
                        Thresholds {
                            background,
                            noise,
                            detection_level: bg.star_level2,
                            hfd_min,
                        },
                        Region::inset(img),
                        &mut pass_stars,
                    );
                }
            }
            2 => {
                detect_pass(
                    img,
                    &mut img_sa,
                    Thresholds {
                        background,
                        noise,
                        detection_level: 30.0 * noise,
                        hfd_min,
                    },
                    Region::inset(img),
                    &mut pass_stars,
                );
            }
            1 => {
                // Grid-based fallback with local background per section
                let (steps_x, steps_y) = if h < w {
                    (
                        RASTER_STEPS,
                        ((RASTER_STEPS as f64 * h as f64 / w as f64).round() as usize).max(1),
                    )
                } else {
                    (
                        ((RASTER_STEPS as f64 * w as f64 / h as f64).round() as usize).max(1),
                        RASTER_STEPS,
                    )
                };
                for yy in 0..=steps_y {
                    for xx in 0..=steps_x {
                        let sx =
                            1 + (w as f64 * xx as f64 / (steps_x as f64 + 1.0)).round() as usize;
                        let ex = (w - 2).min(
                            (w as f64 * (xx + 1) as f64 / (steps_x as f64 + 1.0)).round() as usize,
                        );
                        let sy =
                            1 + (h as f64 * yy as f64 / (steps_y as f64 + 1.0)).round() as usize;
                        let ey = (h - 2).min(
                            (h as f64 * (yy + 1) as f64 / (steps_y as f64 + 1.0)).round() as usize,
                        );
                        if ex <= sx || ey <= sy {
                            continue;
                        }
                        // `max`, not `min`, is deliberate: the 2x-background term
                        // can only raise the limit, so each cell histograms the full
                        // 16-bit range. Swapping it to `min` reads like the obvious
                        // fix and was measured over the 103-image corpus: same 90/103
                        // and 0 false positives, but total runtime rose from 56.6s to
                        // 62.7s and accuracy was a wash (15 images marginally worse,
                        // 10 better, all sub-arcsec). Do not re-try it.
                        let upper = (65500usize).max((background as usize).saturating_mul(2));
                        let (local_bg, local_noise) = sigma_clipped_mean_from_histogram(
                            img,
                            Region {
                                x0: sx,
                                x1: ex,
                                y0: sy,
                                y1: ey,
                            },
                            upper,
                            6,
                            0.1,
                        );
                        detect_pass(
                            img,
                            &mut img_sa,
                            Thresholds {
                                background: local_bg,
                                noise: local_noise,
                                detection_level: 7.0 * local_noise,
                                hfd_min,
                            },
                            Region {
                                x0: sx,
                                x1: ex,
                                y0: sy,
                                y1: ey,
                            },
                            &mut pass_stars,
                        );
                    }
                }
            }
            _ => {}
        }

        n_large += pass_stars.iter().filter(|f| f.large).count();
        stars.extend(pass_stars.into_iter().map(|f| f.star));
        level -= 1;
    }

    stars
}

/// Scan the whole image with the shared marker map.
///
/// `img_sa[i] == 1` means the pixel is already claimed by a detected star.
fn detect_pass_serial(
    img: &ImageBuffer,
    img_sa: &mut [u8],
    thr: Thresholds,
    region: Region,
    out: &mut Vec<Found>,
) {
    let w = img.width;
    let mut m = FullMarkers { data: img_sa, w };
    detect_pass_scan(img, &mut m, thr, region, out);
}

/// Scan one band with a band-local marker map.
fn detect_pass_banded(
    img: &ImageBuffer,
    markers: &mut BandMarkers<'_>,
    thr: Thresholds,
    region: Region,
    out: &mut Vec<Found>,
) {
    detect_pass_scan(img, markers, thr, region, out);
}

/// Rows of overlap between detection bands.
///
/// A star's exclusion disc has radius `3 * hfd` and `hfd` is capped at 30, so 90
/// rows is the widest a marking can reach. Bands overlap by that much so a star
/// near a boundary is seen whole by at least one band.
const BAND_OVERLAP: usize = 90;

/// Run `detect_pass_serial` over horizontal bands in parallel.
///
/// Detection is the largest single cost in a solve and the frame splits naturally,
/// but the `img_sa` "already detected" map makes the serial pass order-dependent.
/// Each band therefore gets its own marker buffer covering its rows plus
/// `BAND_OVERLAP`, and the results are merged with a positional dedup - a star that
/// straddles a boundary is found by both neighbours and kept once.
///
/// The shared `img_sa` is still updated afterwards so later cascade levels skip what
/// earlier ones found, exactly as before.
///
/// One consequence worth knowing: because the exclusion marking is order-dependent,
/// the star list can differ very slightly near band boundaries depending on how many
/// bands were used. Across the 103-image corpus this changes no solve outcome and no
/// false-positive count, and shifts the fitted position by at most 0.06" (median 0),
/// but it does mean a solve is not bit-reproducible across different `--threads`
/// values. Use `--threads 1` when you need an exactly repeatable result.
fn detect_pass(
    img: &ImageBuffer,
    img_sa: &mut [u8],
    thr: Thresholds,
    region: Region,
    out: &mut Vec<Found>,
) {
    /// Dedup hash cell size (pixels) for merging band results.
    const CELL: f64 = 2.0;

    let threads = crate::max_threads().clamp(1, 32);
    let (y0, y1) = (region.y0, region.y1);
    let rows = region.rows();
    // Bands must be comfortably taller than the overlap or the duplicated work
    // swamps the parallelism.
    let min_band = BAND_OVERLAP * 4;
    let n_bands = threads.min(rows / min_band.max(1)).max(1);

    if n_bands <= 1 {
        detect_pass_serial(img, img_sa, thr, region, out);
        return;
    }

    let band_rows = rows.div_ceil(n_bands);
    let w = img.width;

    let results: Vec<Vec<Found>> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..n_bands)
            .map(|b| {
                let by0 = y0 + b * band_rows;
                let by1 = (by0 + band_rows - 1).min(y1);
                // Grow downward only: a star centred above the band was already
                // handled by the band above, and the dedup catches the overlap.
                let ey1 = (by1 + BAND_OVERLAP).min(y1);
                let sa_y0 = by0.saturating_sub(BAND_OVERLAP);
                let sa_y1 = ey1;
                scope.spawn(move || {
                    if by0 > y1 {
                        return Vec::new();
                    }
                    // Local marker buffer covering only this band's span.
                    let mut local = vec![0u8; w * (sa_y1 - sa_y0 + 1)];
                    let mut view = BandMarkers {
                        data: &mut local,
                        y0: sa_y0,
                        w,
                    };
                    let mut band_out = Vec::new();
                    detect_pass_banded(
                        img,
                        &mut view,
                        thr,
                        region.with_rows(by0, ey1),
                        &mut band_out,
                    );
                    band_out
                })
            })
            .collect();
        handles
            .into_iter()
            // A dead band must not silently contribute zero stars.
            .map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
            .collect()
    });

    // Merge with a positional dedup. A linear scan is O(n^2) and a crowded field
    // yields tens of thousands of stars, while a dense grid over a 4300px frame
    // would allocate millions of buckets - so key a hash map by a 2-pixel cell. A
    // duplicate is within 1 px, so it can only be in the same or an adjacent cell.
    let mut cells: HashMap<(i32, i32), Vec<u32>> = HashMap::new();
    let mut merged: Vec<Found> = Vec::new();

    for band in results {
        'star: for found in band {
            let st = &found.star;
            let xci = st.x.round() as usize;
            let yci = st.y.round() as usize;
            if xci < img.width && yci < img.height && img_sa[yci * img.width + xci] == 1 {
                continue; // found by an earlier cascade level
            }
            let gx = (st.x / CELL) as i32;
            let gy = (st.y / CELL) as i32;
            for oy in -1i32..=1 {
                for ox in -1i32..=1 {
                    if let Some(bucket) = cells.get(&(gx + ox, gy + oy)) {
                        for &i in bucket {
                            let o = &merged[i as usize].star;
                            if (o.x - st.x).abs() < 1.0 && (o.y - st.y).abs() < 1.0 {
                                continue 'star;
                            }
                        }
                    }
                }
            }
            cells.entry((gx, gy)).or_default().push(merged.len() as u32);
            merged.push(found);
        }
    }

    // Mark the accepted stars in the shared map so later cascade levels skip them.
    for found in &merged {
        let st = &found.star;
        let xci = st.x.round() as usize;
        let yci = st.y.round() as usize;
        let radius = found.mark;
        let sqr_r = radius * radius;
        for n in -radius..=radius {
            for m in -radius..=radius {
                if m * m + n * n <= sqr_r {
                    let xi = xci as i32 + m;
                    let yi = yci as i32 + n;
                    if xi >= 0 && yi >= 0 && (xi as usize) < img.width && (yi as usize) < img.height
                    {
                        img_sa[yi as usize * img.width + xi as usize] = 1;
                    }
                }
            }
        }
    }

    out.extend(merged);
}

/// Read/write access to the "already detected" map, in whole-image coordinates.
trait Markers {
    fn get(&self, x: usize, y: usize) -> u8;
    fn set(&mut self, x: usize, y: usize);
}

/// The whole-image marker map.
struct FullMarkers<'a> {
    data: &'a mut [u8],
    w: usize,
}

impl Markers for FullMarkers<'_> {
    #[inline]
    fn get(&self, x: usize, y: usize) -> u8 {
        self.data[y * self.w + x]
    }
    #[inline]
    fn set(&mut self, x: usize, y: usize) {
        self.data[y * self.w + x] = 1;
    }
}

impl Markers for BandMarkers<'_> {
    #[inline]
    fn get(&self, x: usize, y: usize) -> u8 {
        BandMarkers::get(self, x, y)
    }
    #[inline]
    fn set(&mut self, x: usize, y: usize) {
        BandMarkers::set(self, x, y);
    }
}

/// A band-local view of the star-area map, addressed in whole-image coordinates.
struct BandMarkers<'a> {
    data: &'a mut [u8],
    y0: usize,
    w: usize,
}

impl BandMarkers<'_> {
    #[inline]
    fn get(&self, x: usize, y: usize) -> u8 {
        if y < self.y0 {
            return 0;
        }
        let i = (y - self.y0) * self.w + x;
        if i < self.data.len() { self.data[i] } else { 0 }
    }
    #[inline]
    fn set(&mut self, x: usize, y: usize) {
        if y < self.y0 {
            return;
        }
        let i = (y - self.y0) * self.w + x;
        if i < self.data.len() {
            self.data[i] = 1;
        }
    }
}

fn detect_pass_scan<M: Markers>(
    img: &ImageBuffer,
    img_sa: &mut M,
    thr: Thresholds,
    region: Region,
    out: &mut Vec<Found>,
) {
    let Thresholds {
        background,
        noise,
        detection_level,
        hfd_min,
    } = thr;
    let (x0, x1, y0, y1) = (region.x0, region.x1, region.y0, region.y1);
    let w = img.width;
    let h = img.height;

    // Absolute thresholds, so the per-pixel test is one compare against the raw
    // value rather than a subtract-then-compare. This loop runs over every pixel of
    // the frame - ~18 million of them on a 4300px image.
    let detect_abs = background + detection_level;
    // Where `measure_large` has already been run. Allocated on first use: most
    // frames never need it.
    let mut large_visited: Option<Visited> = None;
    let hot_abs = background + 4.0 * noise;

    for fy in y0..=y1 {
        let row = fy * w;
        for fx in x0..=x1 {
            // Threshold first: it rejects almost every pixel and needs no lookup,
            // whereas the marker check is a bounds-checked (and, for a band, offset)
            // load. `||` short-circuits and both operands are pure, so the order is
            // ours to choose.
            if (img.data[row + fx] as f64) <= detect_abs || img_sa.get(fx, fy) == 1 {
                continue;
            }
            // Hot-pixel check: at least 2 of 4 neighbours above 4×noise
            let mut star_pixels = 0u8;
            if fx > 0 && img.data[row + fx - 1] as f64 > hot_abs {
                star_pixels += 1;
            }
            if fx + 1 < w && img.data[row + fx + 1] as f64 > hot_abs {
                star_pixels += 1;
            }
            if fy > 0 && img.data[row - w + fx] as f64 > hot_abs {
                star_pixels += 1;
            }
            if fy + 1 < h && img.data[row + w + fx] as f64 > hot_abs {
                star_pixels += 1;
            }
            if star_pixels < 2 {
                continue;
            }

            let measured = match measure_for_scan(img, fx as i32, fy as i32) {
                Ok((star, _)) => {
                    let mark = (3.0 * star.hfd).round() as i32;
                    Some((star, mark, false))
                }
                Err(Reject::TooLarge | Reject::NotADisc) => {
                    // A bright star: measure it again at its own size.
                    let visited = large_visited.get_or_insert_with(|| Visited::new(img));
                    measure_large(img, fx as i32, fy as i32, detect_abs, visited)
                        .map(|(star, mark)| (star, mark, true))
                }
                Err(Reject::NotAStar) => None,
            };
            if let Some((star, radius, large)) = measured
                && star.snr > 10.0
                && star.hfd > hfd_min
                && star.hfd <= 30.0
            {
                let xci = star.x.round() as usize;
                let yci = star.y.round() as usize;

                // Skip if the computed centroid is already in a marked area
                if xci < w && yci < h && img_sa.get(xci, yci) == 1 {
                    continue;
                }

                // Mark circular star area to prevent re-detection
                let sqr_r = radius * radius;
                for n in -radius..=radius {
                    for m in -radius..=radius {
                        if m * m + n * n <= sqr_r {
                            let xi = xci as i32 + m;
                            let yi = yci as i32 + n;
                            if xi >= 0 && yi >= 0 && (xi as usize) < w && (yi as usize) < h {
                                img_sa.set(xi as usize, yi as usize);
                            }
                        }
                    }
                }
                out.push(Found {
                    star,
                    mark: radius,
                    large,
                });
            }
        }
    }
}

/// In-place median (quickselect).
///
/// `measure_star` calls this twice per candidate star - once for the annulus
/// background and once for its MAD - and a crowded field produces tens of
/// thousands of candidates. A full sort is O(n log n) for a single order
/// statistic; `select_nth_unstable_by` is O(n), and the two sorts were 18% of
/// total runtime in a flamegraph of a 1-degree field.
pub(crate) fn median_f64(v: &mut [f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    let len = v.len();
    let mid = len / 2;
    let (lo, nth, _) = v.select_nth_unstable_by(mid, f64::total_cmp);
    if len % 2 == 1 {
        *nth
    } else {
        // Even length: average the two central values. Everything <= the nth is in
        // `lo`, so the value below the median is that partition's maximum.
        let lower = lo.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        (lower + *nth) / 2.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[allow(unused_imports)]
    use core::f64::consts::PI;

    /// Generate a synthetic Gaussian star PSF centred at (cx, cy) with given sigma.
    fn add_star(img: &mut ImageBuffer, cx: f64, cy: f64, sigma: f64, peak: f32) {
        let rs = (4.0 * sigma).ceil() as i32;
        for dy in -rs..=rs {
            for dx in -rs..=rs {
                let x = cx + dx as f64;
                let y = cy + dy as f64;
                if x >= 0.0 && y >= 0.0 && (x as usize) < img.width && (y as usize) < img.height {
                    let r2 =
                        (dx as f64 * dx as f64 + dy as f64 * dy as f64) / (2.0 * sigma * sigma);
                    let flux = peak as f64 * (-r2).exp();
                    let xi = x as usize;
                    let yi = y as usize;
                    img.data[yi * img.width + xi] += flux as f32;
                }
            }
        }
    }

    fn make_background_image(width: usize, height: usize, bg: f32, noise: f32) -> ImageBuffer {
        let mut data = vec![0f32; width * height];
        for (i, v) in data.iter_mut().enumerate() {
            // Deterministic "noise" using a simple hash
            let hash = ((i as u64)
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407)) as u32;
            *v = bg + (hash as f32 / u32::MAX as f32 - 0.5) * 2.0 * noise;
        }
        ImageBuffer {
            data,
            width,
            height,
        }
    }

    #[test]
    fn detect_five_stars() {
        let mut img = make_background_image(256, 256, 1000.0, 50.0);
        let positions = [
            (60.0, 60.0),
            (120.0, 80.0),
            (80.0, 160.0),
            (180.0, 120.0),
            (200.0, 200.0),
        ];
        for &(x, y) in &positions {
            add_star(&mut img, x, y, 2.0, 8000.0);
        }
        let stars = find_stars(&img, 0.5, 500);
        assert_eq!(stars.len(), 5, "expected 5 stars, found {}", stars.len());
        // Every detected star should be within 1.5 px of a planted star
        for s in &stars.0 {
            let near = positions
                .iter()
                .any(|&(px, py)| ((s.x - px).powi(2) + (s.y - py).powi(2)).sqrt() < 1.5);
            assert!(
                near,
                "star at ({:.1},{:.1}) not near any planted star",
                s.x, s.y
            );
        }
    }

    #[test]
    fn hot_pixel_rejected() {
        let mut img = make_background_image(100, 100, 1000.0, 30.0);
        // Single hot pixel — should be rejected
        img.data[50 * 100 + 50] = 60000.0;
        let stars = find_stars(&img, 0.5, 500);
        assert_eq!(stars.len(), 0, "hot pixel should not be detected as a star");
    }

    #[test]
    fn snr_cutoff_respected() {
        let mut img = make_background_image(128, 128, 1000.0, 100.0);
        // Very faint star that won't clear SNR>10
        add_star(&mut img, 64.0, 64.0, 2.0, 200.0);
        let stars = find_stars(&img, 0.5, 500);
        assert_eq!(stars.len(), 0, "faint star should be below SNR threshold");
    }

    #[test]
    fn hfd_cutoff_respected() {
        // Tiny sharp star below hfd_min should be rejected
        let mut img = make_background_image(128, 128, 1000.0, 50.0);
        add_star(&mut img, 64.0, 64.0, 0.3, 10000.0); // very sharp → HFD ~0.7 px
        let stars = find_stars(&img, 1.5, 500); // hfd_min = 1.5
        // May or may not detect; if detected HFD must be > 1.5
        for s in &stars.0 {
            assert!(s.hfd > 1.5, "hfd = {} should be > hfd_min 1.5", s.hfd);
        }
    }

    #[test]
    fn max_stars_limit() {
        let mut img = make_background_image(512, 512, 1000.0, 30.0);
        // Plant 20 bright stars
        for i in 0..20 {
            let x = 50.0 + (i as f64 % 5.0) * 90.0;
            let y = 50.0 + (i as f64 / 5.0).floor() * 90.0;
            add_star(&mut img, x, y, 2.0, 12000.0);
        }
        let stars = find_stars(&img, 0.5, 10);
        assert!(
            stars.len() <= 10,
            "should not exceed max_stars=10, got {}",
            stars.len()
        );
    }

    #[test]
    fn measure_star_rejects_nan_neighbourhood() {
        // A NaN patch used to propagate into HFD/SNR and then panic the
        // brightness sort ("comparison function does not implement a total order").
        let w = 64;
        let mut img = ImageBuffer::new(w, w);
        img.data.fill(100.0);
        // A bright blob with NaNs punched through it.
        for dy in -3i32..=3 {
            for dx in -3i32..=3 {
                let x = (32 + dx) as usize;
                let y = (32 + dy) as usize;
                img.data[y * w + x] = 30000.0;
            }
        }
        img.data[32 * w + 32] = f32::NAN;
        img.data[32 * w + 33] = f32::NAN;
        let s = measure_star(&img, 32, 32);
        assert!(
            s.is_none_or(|st| st.x.is_finite()
                && st.y.is_finite()
                && st.snr.is_finite()
                && st.hfd.is_finite()),
            "measure_star must never return a non-finite star"
        );
    }

    #[test]
    fn find_stars_survives_nan_pixels() {
        let w = 128;
        let mut img = ImageBuffer::new(w, w);
        img.data.fill(100.0);
        for (i, &(cx, cy)) in [(30usize, 30usize), (80, 40), (50, 90), (100, 100)]
            .iter()
            .enumerate()
        {
            let amp = 20000.0 + i as f32 * 1000.0;
            for dy in -2i32..=2 {
                for dx in -2i32..=2 {
                    img.data[(cy as i32 + dy) as usize * w + (cx as i32 + dx) as usize] = amp;
                }
            }
        }
        // Scatter NaNs across the frame, including inside a star.
        for k in 0..200 {
            img.data[(k * 37) % (w * w)] = f32::NAN;
        }
        img.data[30 * w + 30] = f32::NAN;
        // Must not panic.
        let bg = get_background(&img, 100);
        let _ = find_stars_with_background(&img, &bg, 1.0, 100, w, w);
    }

    #[test]
    fn median_f64_ignores_nothing_but_never_panics_on_nan() {
        let mut v = vec![3.0, f64::NAN, 1.0, 2.0];
        let _ = median_f64(&mut v); // total_cmp: must not panic
    }

    #[test]
    fn round_sqrt_matches_float() {
        for n in 0..6000i32 {
            let want = ((n as f64).sqrt().round()) as usize;
            assert_eq!(round_sqrt(n), want, "n={n}");
        }
    }

    #[test]
    fn median_f64_basic() {
        let mut v = vec![3.0, 1.0, 2.0];
        assert_eq!(median_f64(&mut v), 2.0);
        let mut v2 = vec![4.0, 1.0, 3.0, 2.0];
        assert_eq!(median_f64(&mut v2), 2.5);
    }

    /// Degenerate image sizes must yield no stars rather than panic.
    #[test]
    fn tiny_images_yield_no_stars() {
        for (w, h) in [(0, 0), (1, 1), (2, 2), (0, 5), (5, 0), (2, 10)] {
            let img = ImageBuffer::new(w, h);
            let stars = find_stars(&img, 1.0, 100);
            assert!(stars.is_empty(), "{w}x{h}");
        }
    }

    /// A frame tall enough to be split into detection bands (the split needs 4 ×
    /// `BAND_OVERLAP` rows per band), with a star every 25 rows so several sit on
    /// or near every band boundary. Each must be found exactly once, where it is,
    /// and the banded pass must agree with the serial one.
    #[test]
    fn banded_detection_finds_each_star_once() {
        let (w, h) = (120usize, 2000usize);
        let mut img = make_background_image(w, h, 1000.0, 10.0);
        let mut truth = Vec::new();
        for k in 0..78 {
            // Whole-pixel centres: the `add_star` helper samples on integer offsets.
            let (x, y) = (30.0 + (k % 3) as f64 * 30.0, 20.0 + k as f64 * 25.0);
            add_star(&mut img, x, y, 1.5, 8000.0);
            truth.push((x, y));
        }
        let bg = get_background(&img, 500);
        let thr = Thresholds {
            background: bg.mean,
            noise: bg.noise,
            detection_level: 30.0 * bg.noise,
            hfd_min: 0.8,
        };

        let mut banded = Vec::new();
        detect_pass(
            &img,
            &mut vec![0u8; w * h],
            thr,
            Region::inset(&img),
            &mut banded,
        );
        let mut serial = Vec::new();
        detect_pass_serial(
            &img,
            &mut vec![0u8; w * h],
            thr,
            Region::inset(&img),
            &mut serial,
        );

        for (name, found) in [("banded", &banded), ("serial", &serial)] {
            assert_eq!(found.len(), truth.len(), "{name}");
            for &(x, y) in &truth {
                let n = found
                    .iter()
                    .filter(|f| (f.star.x - x).abs() < 0.2 && (f.star.y - y).abs() < 0.2)
                    .count();
                assert_eq!(n, 1, "{name}: star at ({x}, {y}) found {n} times");
            }
        }

        // The public entry point marks what it found, so a second cascade level
        // does not find the same stars again.
        let (stars, raw) = find_stars_with_background(&img, &bg, 0.8, 500, w, h);
        assert_eq!(raw, truth.len());
        assert_eq!(stars.len(), truth.len());
    }

    /// With more candidates than `max_stars`, the list is trimmed to the highest
    /// SNR and sorted by it; with fewer, it is returned in detection order.
    #[test]
    fn trimming_keeps_the_highest_snr() {
        let mut img = make_background_image(300, 300, 1000.0, 10.0);
        for k in 0..30 {
            let (x, y) = (30.0 + (k % 6) as f64 * 45.0, 30.0 + (k / 6) as f64 * 55.0);
            add_star(&mut img, x, y, 1.5, 500.0 + 300.0 * k as f32);
        }
        let bg = get_background(&img, 10);
        let (top, raw) = find_stars_with_background(&img, &bg, 0.8, 10, 300, 300);
        // The cascade stops at the first level that finds enough, so `raw` counts
        // that level's stars, not every star in the frame.
        assert!(raw > 10 && raw <= 30, "raw {raw}");
        assert_eq!(top.len(), 10);
        assert!(top.0.windows(2).all(|w| w[0].snr >= w[1].snr));
        // The ten kept are the ten brightest planted, k = 20..30.
        for s in &top.0 {
            let k = ((s.y - 30.0) / 55.0).round() * 6.0 + ((s.x - 30.0) / 45.0).round();
            assert!(k >= 20.0, "kept k = {k}: {s:?}");
        }
    }

    /// A saturated disc as a photographic plate records a bright star: flat at
    /// `level` out to `radius`, then a Gaussian edge of width `edge`.
    fn add_saturated_disc(img: &mut ImageBuffer, cx: f64, cy: f64, radius: f64, level: f32) {
        let edge = 2.0;
        let rs = (radius + 5.0 * edge).ceil() as i32;
        for dy in -rs..=rs {
            for dx in -rs..=rs {
                let (x, y) = (cx.round() as i32 + dx, cy.round() as i32 + dy);
                if x < 0 || y < 0 || x as usize >= img.width || y as usize >= img.height {
                    continue;
                }
                let r = (f64::from(x) - cx).hypot(f64::from(y) - cy);
                let v = if r <= radius {
                    1.0
                } else {
                    (-(r - radius).powi(2) / (2.0 * edge * edge)).exp()
                };
                let i = y as usize * img.width + x as usize;
                img.data[i] = (img.data[i] + level * v as f32).min(30_000.0);
            }
        }
    }

    /// A bright star too large for the 14-pixel measuring box (a saturated disc 24
    /// pixels across) used to be refused; it is measured again in a larger box,
    /// at its true centre, and the ordinary stars around it are found as before.
    #[test]
    fn a_saturated_disc_is_measured_in_a_larger_box() {
        let mut img = make_background_image(300, 300, 1000.0, 20.0);
        add_saturated_disc(&mut img, 150.3, 140.6, 12.0, 25_000.0);
        let small = [(50.0, 50.0), (250.0, 60.0), (60.0, 250.0), (240.0, 240.0)];
        for &(x, y) in &small {
            add_star(&mut img, x, y, 1.5, 3000.0);
        }
        let stars = find_stars(&img, 0.8, 500);
        let big = stars
            .0
            .iter()
            .find(|s| (s.x - 150.3).hypot(s.y - 140.6) < 1.0)
            .expect("the saturated disc is found, centred");
        assert!(
            stars.0.iter().all(|s| s.snr <= big.snr),
            "and is the brightest"
        );
        for &(x, y) in &small {
            assert!(stars.0.iter().any(|s| (s.x - x).hypot(s.y - y) < 1.0));
        }
        assert_eq!(stars.len(), small.len() + 1, "{stars:?}");
    }

    /// A trail is not measured as a large star, however bright: the second-moment
    /// test refuses it.
    #[test]
    fn a_bright_trail_is_not_a_large_star() {
        let mut img = make_background_image(200, 200, 1000.0, 20.0);
        for k in 0..60 {
            add_star(
                &mut img,
                70.0 + k as f64,
                100.0 + 0.2 * k as f64,
                1.5,
                20_000.0,
            );
        }
        let stars = find_stars(&img, 0.5, 500);
        assert!(stars.is_empty(), "{stars:?}");
    }

    /// Stars only the large measurement finds do not end the detection cascade:
    /// the same levels run and find the same ordinary stars as without them.
    #[test]
    fn large_stars_do_not_cut_the_cascade_short() {
        let mut img = make_background_image(400, 400, 1000.0, 10.0);
        // Ten bright saturated discs, found at the first level.
        for k in 0..10 {
            add_saturated_disc(&mut img, 40.0 + 35.0 * k as f64, 40.0, 8.0, 25_000.0);
        }
        // Thirty faint stars, below the first levels' thresholds.
        let mut faint = Vec::new();
        for k in 0..30 {
            let (x, y) = (40.0 + (k % 6) as f64 * 60.0, 120.0 + (k / 6) as f64 * 55.0);
            add_star(&mut img, x, y, 1.5, 400.0);
            faint.push((x, y));
        }
        let bg = get_background(&img, 10);
        let (stars, _) = find_stars_with_background(&img, &bg, 0.8, 10, 400, 400);
        let all = detect_all(&img, &bg, 0.8, 10);
        // With max_stars = 10 the ten discs alone would have ended it.
        let n_faint = all
            .iter()
            .filter(|s| faint.iter().any(|&(x, y)| (s.x - x).hypot(s.y - y) < 1.0))
            .count();
        assert_eq!(n_faint, faint.len());
        assert_eq!(stars.len(), 10);
    }

    /// The deep list holds every star found, brightest first; the ordinary list is
    /// what `find_stars_with_background` returns.
    #[test]
    fn the_deep_list_extends_the_ordinary_one() {
        let mut img = make_background_image(300, 300, 1000.0, 10.0);
        for k in 0..30 {
            let (x, y) = (30.0 + (k % 6) as f64 * 45.0, 30.0 + (k / 6) as f64 * 55.0);
            add_star(&mut img, x, y, 1.5, 500.0 + 300.0 * k as f32);
        }
        let bg = get_background(&img, 10);
        let (top, raw) = find_stars_with_background(&img, &bg, 0.8, 10, 300, 300);
        let (top2, raw2, deep) = find_stars_and_deep(&img, &bg, 0.8, 10, 1000);
        assert_eq!(raw, raw2);
        assert_eq!(deep.len(), raw);
        assert!(deep.0.windows(2).all(|w| w[0].snr >= w[1].snr));
        for (a, b) in top.0.iter().zip(&top2.0) {
            assert_eq!((a.x, a.y), (b.x, b.y));
        }
        for (a, b) in top.0.iter().zip(&deep.0) {
            assert_eq!((a.x, a.y), (b.x, b.y));
        }
    }
}
