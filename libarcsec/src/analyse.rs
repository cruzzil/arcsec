//! `arcsec_analyse`: star count and HFD without solving (the CLI's `--analyse`).

use arcsec_core::detection::analyse_image;

use crate::error::{Failure, guard};
use crate::image::{arcsec_image, read_image};
use crate::util::{Versioned, write_versioned};

/// A star measured by `arcsec_analyse`.
#[repr(C)]
#[derive(Debug, Clone, Copy, Default)]
pub struct arcsec_star {
    /// Centroid column, 1-based FITS pixels.
    pub x: f64,
    /// Centroid row, 1-based FITS pixels.
    pub y: f64,
    /// Half-flux diameter, pixels.
    pub hfd: f64,
    /// Signal-to-noise ratio.
    pub snr: f64,
    /// Background-subtracted flux, in the units of the (normalised) pixels.
    pub flux: f64,
}

/// What `arcsec_analyse` found. Set `struct_size` before the call.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct arcsec_analysis {
    /// `sizeof(arcsec_analysis)`, set by the caller.
    pub struct_size: usize,
    /// Stars found.
    pub star_count: u64,
    /// Median HFD of those stars, pixels; 0 if there are none.
    pub hfd_median: f64,
    /// Background level, in the units of the (normalised) pixels.
    pub background: f64,
    /// Background noise (standard deviation), same units.
    pub noise: f64,
}

// SAFETY: repr(C), starts with struct_size, all fields numeric.
unsafe impl Versioned for arcsec_analysis {
    fn defaults() -> Self {
        Self {
            struct_size: core::mem::size_of::<Self>(),
            star_count: 0,
            hfd_median: 0.0,
            background: 0.0,
            noise: 0.0,
        }
    }
}

/// Detect and measure the stars of an image without solving it, as `arcsec
/// --analyse` does: their number and median HFD (ASTAP's `HFD_MEDIAN` and
/// `STARS`), and optionally the stars themselves.
///
/// `snr_min` is the smallest signal-to-noise ratio counted (0 for ASTAP's 30);
/// `max_stars` the number of stars detection aims for (0 for 500). Up to
/// `stars_capacity` stars are copied to `stars` (which may be NULL when the
/// capacity is 0), in detection order; `out->star_count` is the full number.
///
/// # Safety
///
/// `image` as for `arcsec_solve`; `out` must be NULL or point to
/// `out->struct_size` writable bytes; `stars` must be NULL or have room for
/// `stars_capacity` elements.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn arcsec_analyse(
    image: *const arcsec_image,
    snr_min: f64,
    max_stars: u32,
    out: *mut arcsec_analysis,
    stars: *mut arcsec_star,
    stars_capacity: usize,
) -> crate::arcsec_status {
    guard(|| {
        if out.is_null() {
            return Err(Failure::invalid("analysis output is NULL"));
        }
        if !(snr_min.is_finite() && snr_min >= 0.0) {
            return Err(Failure::invalid("snr_min must be zero or positive"));
        }
        // SAFETY: forwarded contract.
        let (mut img, _) = unsafe { read_image(image) }?;
        img.normalize_for_detection();
        let snr_min = if snr_min == 0.0 { 30.0 } else { snr_min };
        let max_stars = if max_stars == 0 {
            500
        } else {
            max_stars as usize
        };
        let a = analyse_image(&img, snr_min, max_stars);
        let v = arcsec_analysis {
            star_count: a.stars.len() as u64,
            hfd_median: a.hfd_median().unwrap_or(0.0),
            background: a.background.mean,
            noise: a.background.noise,
            ..arcsec_analysis::defaults()
        };
        if !stars.is_null() {
            for (i, s) in a.stars.iter().take(stars_capacity).enumerate() {
                let star = arcsec_star {
                    x: s.x + 1.0,
                    y: s.y + 1.0,
                    hfd: s.hfd,
                    snr: s.snr,
                    flux: s.flux,
                };
                // SAFETY: i < stars_capacity elements fit.
                unsafe { stars.add(i).write_unaligned(star) };
            }
        }
        // SAFETY: forwarded contract.
        unsafe {
            write_versioned(
                out,
                &v,
                "arcsec_analysis",
                core::mem::size_of::<arcsec_analysis>(),
            )
        }
    })
}
