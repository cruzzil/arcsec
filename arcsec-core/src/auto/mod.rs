//! Everything the `arcsec` command line decides for the user, as a library.
//!
//! [`crate::pipeline::solve_image`] needs to be told everything: the field size in
//! radians, the binning, the database directory and name, the minimum star size in
//! (binned) pixels. A user — of the CLI, or of an application embedding arcsec —
//! knows less than that and expects the rest to be worked out the way the CLI does
//! it. This module is that working-out, shared by the `arcsec` binary and the C
//! library so the two cannot drift:
//!
//! - where catalogues live ([`default_catalog_dir`], [`default_db_path`]);
//! - which star database suits a field ([`select_db_for_fov`]);
//! - how far to bin ([`choose_binning`]);
//! - which pixel scales to try when the scale is not known ([`ScaleSearch`],
//!   [`ladder`]);
//! - when to use a blind index, and how ([`find_arcsec_index`],
//!   [`collect_index_files`]);
//! - and the whole solve built from those decisions: a [`SolveRequest`] becomes a
//!   [`Plan`], and [`Plan::solve`] runs it.
//!
//! ```no_run
//! use arcsec_core::auto::{Plan, SolveRequest};
//! # fn load() -> arcsec_core::ImageBuffer { arcsec_core::ImageBuffer::new(4096, 3072) }
//! let mut img = load();
//! img.normalize_for_detection();
//! let request = SolveRequest {
//!     hint: Some((83.82_f64.to_radians(), (-5.39_f64).to_radians())),
//!     pixel_scale: Some(1.1),
//!     search_radius: 10.0_f64.to_radians(),
//!     ..SolveRequest::default()
//! };
//! let plan = Plan::new(&request, img.width, img.height)?;
//! let solved = plan.solve(&img)?;
//! println!("RA {:.4}°", solved.wcs.ra0.to_degrees());
//! # Ok::<(), arcsec_core::ArcsecError>(())
//! ```

mod blind;
mod db;
mod scale;

use alloc::borrow::Cow;
use core::f64::consts::PI;
use std::path::{Path, PathBuf};

pub use blind::{
    SOURCES, collect_index_files, depth_rank, find_arcsec_index, preferred_index,
    wants_installed_index,
};
pub use db::{
    ASTAP_EXTS, DB_FOV_RANGES, available_dbs, default_db_path, has_star_database, select_db_for_fov,
};
pub use scale::{
    Hypothesis, INACCURATE_SCALE, LADDER_FIELDS, SCALE_STEP, ScaleSearch, UNKNOWN_STEPS,
    WRONG_STEPS, inaccurate_scale_warning, ladder,
};

use crate::cancel::CancelToken;
use crate::error::{ArcsecError, Result};
use crate::pipeline::solver::{
    Detected, ScaleTrust, detect, search_in_order, solve_detected, solve_image_with,
};
use crate::pipeline::{BlindSolveParams, SearchSpeed, SolveMethod, SolveParams, solve_image};
use crate::types::{ImageBuffer, WcsSolution};

// ── Where things live ───────────────────────────────────────────────────────────

/// Where catalogues are kept, in priority order:
///
/// 1. `$ARCSEC_CATALOG_DIR`, if set — for people who keep them on another disk.
/// 2. `$XDG_DATA_HOME/arcsec/catalogs` on Linux, or the platform equivalent:
///    `~/Library/Application Support/arcsec/catalogs` on macOS,
///    `%LOCALAPPDATA%\arcsec\catalogs` on Windows.
/// 3. `~/.arcsec/catalogs` if the home directory cannot be resolved any other way.
///
/// The point is that a user who runs `arcsec catalog install d50` never has to know
/// this path, and the solver looks here without being told.
#[must_use]
pub fn default_catalog_dir() -> PathBuf {
    default_catalog_dir_from(|k| std::env::var(k).ok())
}

/// [`default_catalog_dir`] with the environment supplied by `var`, so it can be
/// tested without mutating the real process environment.
fn default_catalog_dir_from(var: impl Fn(&str) -> Option<String>) -> PathBuf {
    let get = |k: &str| var(k).filter(|v| !v.is_empty()).map(PathBuf::from);

    if let Some(p) = get("ARCSEC_CATALOG_DIR") {
        return p;
    }

    let platform = if cfg!(target_os = "windows") {
        get("LOCALAPPDATA").map(|p| p.join("arcsec").join("catalogs"))
    } else if cfg!(target_os = "macos") {
        get("HOME").map(|h| {
            h.join("Library")
                .join("Application Support")
                .join("arcsec")
                .join("catalogs")
        })
    } else {
        get("XDG_DATA_HOME")
            .map(|p| p.join("arcsec").join("catalogs"))
            .or_else(|| {
                get("HOME").map(|h| {
                    h.join(".local")
                        .join("share")
                        .join("arcsec")
                        .join("catalogs")
                })
            })
    };

    platform.unwrap_or_else(|| {
        get("HOME").or_else(|| get("USERPROFILE")).map_or_else(
            || PathBuf::from("catalogs"),
            |h| h.join(".arcsec").join("catalogs"),
        )
    })
}

/// Every arcsec blind index (`*.arcsecix`) directly in `dir`, sorted by name.
/// Empty when there are none or `dir` cannot be read.
#[must_use]
pub fn index_files(dir: &Path) -> Vec<PathBuf> {
    let mut v: Vec<PathBuf> = std::fs::read_dir(dir)
        .map(|rd| {
            rd.filter_map(core::result::Result::ok)
                .map(|e| e.path())
                .filter(|p| {
                    p.extension()
                        .is_some_and(|e| e == crate::index::format::EXTENSION)
                })
                .collect()
        })
        .unwrap_or_default();
    v.sort();
    v
}

// ── Binning ─────────────────────────────────────────────────────────────────────

/// Smallest image side, in (binned) pixels, that reaches the solver. Detection
/// cannot run on a one-pixel-wide image, and would otherwise panic on it.
pub const MIN_SOLVE_DIM: usize = 2;

/// The binning factor: `requested` if given and non-zero, else automatic.
///
/// Automatic binning brings a sampling finer than 1"/px back to about 1"/px, up to
/// 16×. Either way the factor is capped so the binned image keeps at least
/// [`MIN_SOLVE_DIM`] pixels a side: binning past the image size leaves nothing to
/// detect in, and used to crash.
#[must_use]
pub fn choose_binning(
    requested: Option<usize>,
    arcsec_per_px: f64,
    width: usize,
    height: usize,
) -> usize {
    let binning = match requested {
        Some(0) | None => {
            if arcsec_per_px < 1.0 {
                (1.0 / arcsec_per_px).round().clamp(1.0, 16.0) as usize
            } else {
                1
            }
        }
        Some(z) => z,
    };
    binning.min((width.min(height) / MIN_SOLVE_DIM).max(1))
}

// ── The request and its plan ────────────────────────────────────────────────────

/// What the caller knows about an image and how it wants it solved, in the
/// caller's terms. Every field has the CLI's default ([`SolveRequest::default`]).
#[derive(Debug, Clone)]
pub struct SolveRequest {
    /// Approximate centre (RA, Dec), radians. `None` for no hint: the search
    /// starts at (0, 0), and an installed blind index is used as for any wide
    /// search.
    pub hint: Option<(f64, f64)>,
    /// Field of view along the image *height*, radians (ASTAP's `-fov`). Takes
    /// precedence over [`Self::pixel_scale`].
    pub fov_height: Option<f64>,
    /// Pixel scale of the unbinned image, arcseconds per pixel (for instance from
    /// the FOCALLEN and XPIXSZ header keywords). Without it or a field of view, 1″/px
    /// is assumed, and [`Self::scale_search`] says whether other scales are tried.
    pub pixel_scale: Option<f64>,
    /// When the catalogue search tries other pixel scales: by default only when
    /// neither [`Self::fov_height`] nor [`Self::pixel_scale`] gives one.
    pub scale_search: ScaleSearch,
    /// Search radius around the hint, radians. Negative or NaN means 0.
    pub search_radius: f64,
    /// Binning factor; `None` or `Some(0)` chooses one ([`choose_binning`]).
    pub downsample: Option<usize>,
    /// Star database directory; `None` for [`default_db_path`].
    pub db_path: Option<PathBuf>,
    /// Star database name (`"d50"`); `None` to choose by field size
    /// ([`select_db_for_fov`]).
    pub db_name: Option<String>,
    /// Blind index to use (the CLI's `--index`): an arcsec `.arcsecix` file, a
    /// directory holding one, or Astrometry.net `index-*.fits` files. `None` still
    /// uses an arcsec index installed in the catalogue directory for a wide search,
    /// unless [`Self::auto_index`] is off.
    pub index: Option<PathBuf>,
    /// With a hint, whether the index named in [`Self::index`] is tried before
    /// the search round the hint (the CLI's `--index`: true, the default), or only
    /// after the spiral has searched a few fields round it (false), as an
    /// installed index is. Without a hint the index always comes first.
    pub index_first: bool,
    /// Consult an arcsec index installed in the catalogue directory (or beside the
    /// star database) when the search is wider than a few fields and at least 10°.
    /// On by default, as in the CLI.
    pub auto_index: bool,
    /// Minimum star size (HFD), arcseconds.
    pub hfd_min_arcsec: f64,
    /// Pattern-matching tolerance.
    pub quad_tolerance: f64,
    /// Maximum number of image stars to use.
    pub max_stars: usize,
    /// Pattern-matching algorithm.
    pub method: SolveMethod,
    /// Catalogue window per spiral position.
    pub speed: SearchSpeed,
    /// Worker threads for this solve; 0 for the process-wide limit
    /// ([`crate::max_threads`]).
    pub threads: usize,
    /// Fit SIP distortion polynomials to the solution ([`crate::wcs::fit_sip`]).
    pub sip: bool,
    /// Stop early when this token is cancelled ([`ArcsecError::Cancelled`]).
    pub cancel: Option<CancelToken>,
}

impl Default for SolveRequest {
    /// The CLI's defaults: no hint, a 180° radius (the whole sky), 500 stars,
    /// tolerance 0.007, minimum HFD 1.5″, automatic binning and database.
    fn default() -> Self {
        Self {
            hint: None,
            fov_height: None,
            pixel_scale: None,
            scale_search: ScaleSearch::default(),
            search_radius: PI,
            downsample: None,
            db_path: None,
            db_name: None,
            index: None,
            index_first: true,
            auto_index: true,
            hfd_min_arcsec: 1.5,
            quad_tolerance: 0.007,
            max_stars: 500,
            method: SolveMethod::Quads,
            speed: SearchSpeed::Auto,
            threads: 0,
            sip: false,
            cancel: None,
        }
    }
}

/// A [`SolveRequest`] with every decision made, for an image of a given size.
///
/// Built by [`Plan::new`] without touching the pixels, so a caller can report what
/// will happen (the CLI prints it, as ASTAP does) before [`Plan::solve`] runs.
#[derive(Debug, Clone)]
pub struct Plan {
    /// Start of the search (RA, Dec), radians: the hint, or (0, 0).
    pub start: (f64, f64),
    /// Whether the request had a hint.
    pub has_hint: bool,
    /// Pixel scale of the unbinned image, arcseconds per pixel.
    pub arcsec_per_px: f64,
    /// Whether the scale came from the request rather than the 1″/px fallback.
    pub scale_known: bool,
    /// Field of view along the image height, radians.
    pub fov_height: f64,
    /// Binning factor the image is solved at.
    pub binning: usize,
    /// Unbinned image size, pixels.
    pub image_size: (usize, usize),
    /// Minimum star size, arcseconds (as requested).
    pub hfd_min_arcsec: f64,
    /// The catalogue solve's parameters: the start, the field along the longer
    /// side, the database, and the minimum HFD in binned pixels.
    pub params: SolveParams,
    index: Option<PathBuf>,
    index_first: bool,
    auto_index: bool,
    sip: bool,
    cancel: Option<CancelToken>,
    /// The request, for the plans of other scales ([`ScaleSearch`]).
    request: SolveRequest,
}

/// Something [`Plan::solve_with`] reports while it runs.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum Event {
    /// The Astrometry.net blind solver estimated the position (RA, Dec), radians;
    /// the catalogue search starts there.
    IndexEstimate(f64, f64),
}

/// A solution from [`Plan::solve`].
#[derive(Debug, Clone)]
pub struct Solved {
    /// The WCS, on the unbinned image's pixel grid; SIP included if requested
    /// and the fit was worthwhile.
    pub wcs: WcsSolution,
    /// The Astrometry.net blind solver's position estimate, if one was used.
    pub index_estimate: Option<(f64, f64)>,
}

impl Plan {
    /// Make every decision for a `width` × `height` image (unbinned).
    ///
    /// # Errors
    ///
    /// [`ArcsecError::InvalidParameter`] for an empty image or a field of view or
    /// pixel scale that is not positive and finite.
    pub fn new(req: &SolveRequest, width: usize, height: usize) -> Result<Self> {
        if width == 0 || height == 0 {
            return Err(ArcsecError::InvalidParameter(format!(
                "image is {width}×{height} pixels"
            )));
        }
        let bad = |v: f64| !(v.is_finite() && v > 0.0);
        if req.fov_height.is_some_and(bad) || req.pixel_scale.is_some_and(bad) {
            return Err(ArcsecError::InvalidParameter(
                "field of view and pixel scale must be positive".into(),
            ));
        }
        // Priority: an explicit field of view > header optics > 1"/px fallback.
        // `fov` (the field along the longer side) is what database selection and
        // the search window use.
        let naxis = width.max(height) as f64;
        let h = height as f64;
        let (arcsec_per_px, scale_known) = match (req.fov_height, req.pixel_scale) {
            (Some(fov_h), _) => (fov_h.to_degrees() * 3600.0 / h, true),
            (None, Some(ps)) => (ps, true),
            (None, None) => (1.0, false),
        };
        let fov = match req.fov_height {
            Some(fov_h) => fov_h * (naxis / h),
            None => (naxis * arcsec_per_px / 3600.0).to_radians(),
        };
        let binning = choose_binning(req.downsample, arcsec_per_px, width, height);

        let db_path = req.db_path.clone().unwrap_or_else(default_db_path);
        // Resolved from the field size: the D-series covers 0.15°–6°, G05 3°–20°
        // and W08 20°–80°.
        let db_name = req.db_name.clone().unwrap_or_else(|| {
            select_db_for_fov(&db_path, fov.to_degrees()).unwrap_or_else(|| "d80".to_string())
        });
        let hfd_min = (req.hfd_min_arcsec / (binning as f64 * arcsec_per_px)).max(0.8);
        let start = req.hint.unwrap_or((0.0, 0.0));
        Ok(Self {
            start,
            has_hint: req.hint.is_some(),
            arcsec_per_px,
            scale_known,
            fov_height: fov * (h / naxis),
            binning,
            image_size: (width, height),
            hfd_min_arcsec: req.hfd_min_arcsec,
            params: SolveParams {
                ra_hint: start.0,
                dec_hint: start.1,
                fov,
                // ASTAP treats a negative (or NaN) radius as zero and solves at the
                // start position; `max` maps NaN to 0.
                search_radius: req.search_radius.max(0.0),
                quad_tolerance: req.quad_tolerance,
                hfd_min,
                max_stars: req.max_stars,
                db_path,
                db_name,
                binning,
                method: req.method,
                speed: req.speed,
                threads: req.threads,
            },
            index: req.index.clone(),
            index_first: req.index_first,
            auto_index: req.auto_index,
            sip: req.sip,
            cancel: req.cancel.clone(),
            request: req.clone(),
        })
    }

    /// `astap_cli`'s warning for a solution whose scale is not the one the solve
    /// started from ([`inaccurate_scale_warning`]): the scale given, read from the
    /// header, or assumed.
    #[must_use]
    pub fn scale_warning(&self, wcs: &WcsSolution) -> Option<String> {
        let solved = (wcs.cd1_1 * wcs.cd2_2 - wcs.cd1_2 * wcs.cd2_1).abs().sqrt() * 3600.0;
        inaccurate_scale_warning(self.arcsec_per_px, solved, self.image_size.1)
    }

    /// Whether the catalogue search will try other scales should this one not
    /// solve at once: always when the scale is unknown (unless
    /// [`ScaleSearch::Never`]), and after a failure with
    /// [`ScaleSearch::AlsoIfWrong`].
    #[must_use]
    pub fn searches_scales(&self) -> bool {
        match self.request.scale_search {
            ScaleSearch::Never => false,
            ScaleSearch::IfUnknown => !self.scale_known,
            ScaleSearch::AlsoIfWrong => true,
        }
    }

    /// Size of the image as solved, after binning.
    #[must_use]
    pub fn binned_size(&self) -> (usize, usize) {
        let b = self.binning.max(1);
        (self.image_size.0 / b, self.image_size.1 / b)
    }

    /// Solve `img`, the unbinned image this plan was made for, with its pixels
    /// already normalised ([`ImageBuffer::normalize_for_detection`]).
    ///
    /// # Errors
    ///
    /// Everything [`solve_image`] returns, and:
    /// - [`ArcsecError::InvalidParameter`] if `img` is not the size planned for;
    /// - [`ArcsecError::InsufficientStars`] if the binned image is too small;
    /// - [`ArcsecError::IndexNotFound`] if [`SolveRequest::index`] names nothing;
    /// - [`ArcsecError::OutsideSearchRadius`] if an index found the field beyond
    ///   the search radius;
    /// - [`ArcsecError::Cancelled`] if the request's token fired.
    pub fn solve(&self, img: &ImageBuffer) -> Result<Solved> {
        self.solve_with(img, |_| {})
    }

    /// [`Plan::solve`], reporting [`Event`]s to `on_event` as they happen.
    ///
    /// # Errors
    ///
    /// As [`Plan::solve`].
    pub fn solve_with(&self, img: &ImageBuffer, on_event: impl FnMut(Event)) -> Result<Solved> {
        if (img.width, img.height) != self.image_size || img.data.len() != img.width * img.height {
            return Err(ArcsecError::InvalidParameter(format!(
                "image is {}×{}, planned for {}×{}",
                img.width, img.height, self.image_size.0, self.image_size.1
            )));
        }
        if !self.scale_known {
            log::warn!(
                "No pixel scale given (a field of view, or FOCALLEN and XPIXSZ in the header): {}",
                if self.searches_scales() {
                    "searching scales from 0.25 to 64\"/px round the hint, then 1\"/px"
                } else {
                    "assuming 1\"/px, which will not solve unless it is roughly right"
                }
            );
        }
        crate::cancel::with_optional(self.cancel.as_ref(), || {
            crate::with_max_threads(self.params.threads, || self.run(img, on_event))
        })
    }

    fn run(&self, unbinned: &ImageBuffer, mut on_event: impl FnMut(Event)) -> Result<Solved> {
        let img = unbinned;
        let (bw, bh) = self.binned_size();
        if bw < MIN_SOLVE_DIM || bh < MIN_SOLVE_DIM {
            return Err(ArcsecError::InsufficientStars {
                found: 0,
                required: 5,
            });
        }
        let img: Cow<'_, ImageBuffer> = if self.binning > 1 {
            log::info!(
                "Creating grayscale x {0} binning image for solving/star alignment.",
                self.binning
            );
            Cow::Owned(img.bin_image(self.binning))
        } else {
            Cow::Borrowed(img)
        };
        let template = &self.params;
        let (ra_hint, dec_hint) = self.start;

        // arcsec's own blind index, named by `index` or, for a search wider than a
        // few fields, found in the catalogue directory: it finds the field and the
        // hinted solver accepts it (see blind::index_stage). With Astrometry.net
        // files, the blind solver estimates the position first, and that estimate
        // becomes the hint for the catalogue spiral solver.
        // A hint, and an index named only as a fallback: the index waits until the
        // spiral has searched round the hint, exactly as an installed one does.
        let fallback = self.has_hint && !self.index_first;
        let own_index = blind::arcsec_index_for(self.index.as_ref(), template, self.auto_index)
            .map(|mut ix| {
                ix.explicit &= !fallback;
                ix
            });
        let index_wcs = match own_index.as_ref().map(|ix| {
            blind::index_stage(
                &img,
                ix,
                template,
                self.has_hint,
                self.arcsec_per_px * self.binning as f64,
                self.scale_known,
            )
        }) {
            Some(blind::IndexOutcome::Solved(w)) => Some(*w),
            Some(blind::IndexOutcome::Elsewhere(separation_deg)) => {
                return Err(ArcsecError::OutsideSearchRadius { separation_deg });
            }
            Some(blind::IndexOutcome::Cancelled) => return Err(ArcsecError::Cancelled),
            Some(blind::IndexOutcome::NotFound) | None => None,
        };

        let mut index_estimate = None;
        let (ra, dec, search_radius) = match &self.index {
            None => (ra_hint, dec_hint, template.search_radius),
            _ if index_wcs.is_some() => (ra_hint, dec_hint, template.search_radius),
            Some(_) if own_index.is_some() => {
                log::warn!("Blind index found no verified position. Falling back to hint.");
                (ra_hint, dec_hint, template.search_radius)
            }
            Some(idx_root) => {
                let index_files = collect_index_files(idx_root, template.fov.to_degrees());
                if index_files.is_empty() {
                    return Err(ArcsecError::IndexNotFound(idx_root.clone()));
                }
                // As a fallback, the index waits for a search round the hint.
                if fallback {
                    let near = SolveParams {
                        search_radius: blind::stage_one_radius(template)
                            .min(template.search_radius),
                        ..template.clone()
                    };
                    log::info!(
                        "Searching {:.1}° round the hint before the blind index.",
                        near.search_radius.to_degrees()
                    );
                    match solve_image(&img, &near) {
                        Ok(mut wcs) => {
                            if self.sip {
                                wcs.sip =
                                    crate::wcs::fit_sip(&wcs, self.image_size.0, self.image_size.1);
                            }
                            return Ok(Solved {
                                wcs,
                                index_estimate: None,
                            });
                        }
                        Err(
                            e @ (ArcsecError::Cancelled
                            | ArcsecError::CatalogNotFound(_)
                            | ArcsecError::CatalogIo(_)
                            | ArcsecError::InsufficientStars { .. }),
                        ) => return Err(e),
                        Err(_) => {}
                    }
                }
                let params = BlindSolveParams {
                    quad_tolerance: template.quad_tolerance,
                    hfd_min: template.hfd_min,
                    max_stars: template.max_stars,
                    binning: self.binning,
                    // The blind scale filter maps this through the image height.
                    fov_deg: self.fov_height.to_degrees(),
                };
                match blind::estimate_position(&img, &index_files, &params) {
                    blind::BlindOutcome::Found(ra, dec) => {
                        index_estimate = Some((ra, dec));
                        on_event(Event::IndexEstimate(ra, dec));
                        // Narrow the catalog search so the spiral checks step 0 (the
                        // blind position) and at most a few neighbours: the blind
                        // position is off by at most one image width, so 2× fov is a
                        // generous ceiling.
                        (ra, dec, (template.fov * 2.0).max(5.0_f64.to_radians()))
                    }
                    blind::BlindOutcome::InsufficientStars { found, required } => {
                        return Err(ArcsecError::InsufficientStars { found, required });
                    }
                    blind::BlindOutcome::NotFound => {
                        if crate::cancel::is_cancelled() {
                            return Err(ArcsecError::Cancelled);
                        }
                        log::warn!(
                            "Blind position estimate failed for all index files. Falling back to hint."
                        );
                        (ra_hint, dec_hint, template.search_radius)
                    }
                }
            }
        };

        let mut wcs = match index_wcs {
            Some(w) => w,
            None => self.catalogue_solve(unbinned, &img, (ra, dec), search_radius)?,
        };
        if self.sip {
            wcs.sip = crate::wcs::fit_sip(&wcs, self.image_size.0, self.image_size.1);
        }
        Ok(Solved {
            wcs,
            index_estimate,
        })
    }
}

impl Plan {
    /// The catalogue search from `start` out to `radius`
    /// ([`Self::catalogue_search`]), then, if the solution's scale is more than
    /// [`INACCURATE_SCALE`] from the one the search used, the same field again at
    /// the solved scale ([`Self::refine`]).
    fn catalogue_solve(
        &self,
        unbinned: &ImageBuffer,
        img: &ImageBuffer,
        start: (f64, f64),
        radius: f64,
    ) -> Result<WcsSolution> {
        let wcs = self.catalogue_search(unbinned, img, start, radius)?;
        Ok(self.refine(unbinned, &wcs).unwrap_or(wcs))
    }

    /// `wcs` solved again at its own pixel scale, at its own centre, when that
    /// scale is more than [`INACCURATE_SCALE`] from the scale this plan searched
    /// with; `None` when it is not, or when the second solve does not verify, or
    /// lands elsewhere (more than a tenth of a field away).
    ///
    /// The scale sets the catalogue window, its depth and the density the star
    /// list is trimmed to, and the verification and distortion model work within
    /// that window: a solution found at half the true scale has been checked
    /// against the catalogue of the middle quarter of the frame. Solved again at
    /// the right scale it is the solution a correct scale would have given (on the
    /// benchmark, identical to it in 93 of 95 images, and never more than 0.002″
    /// apart at the corners; without this, up to 1.7″ worse, and three near-misses
    /// past the 5″ limit). The cost is one detection and one catalogue position,
    /// and only when the scale was off.
    fn refine(&self, unbinned: &ImageBuffer, wcs: &WcsSolution) -> Option<WcsSolution> {
        let solved = (wcs.cd1_1 * wcs.cd2_2 - wcs.cd1_2 * wcs.cd2_1).abs().sqrt() * 3600.0;
        let off = (solved / self.arcsec_per_px - 1.0).abs();
        if off.is_nan() || off <= INACCURATE_SCALE {
            return None;
        }
        let (w, h) = self.image_size;
        let req = SolveRequest {
            hint: Some((wcs.ra0, wcs.dec0)),
            fov_height: None,
            pixel_scale: Some(solved),
            scale_search: ScaleSearch::Never,
            search_radius: 0.0,
            index: None,
            auto_index: false,
            sip: false,
            cancel: None,
            ..self.request.clone()
        };
        let p = Plan::new(&req, w, h).ok()?;
        let img: Cow<'_, ImageBuffer> = if p.binning > 1 {
            Cow::Owned(unbinned.bin_image(p.binning))
        } else {
            Cow::Borrowed(unbinned)
        };
        let mut r = solve_image_with(&img, &p.params, ScaleTrust::Hypothesis).ok()?;
        let sep = crate::math::coords::ang_sep(r.ra0, r.dec0, wcs.ra0, wcs.dec0);
        log::info!(
            "Solved again at {solved:.3}\"/px: {} stars verified (first {}), {:.1}\" from the first solution",
            r.stars_matched,
            wcs.stars_matched,
            sep.to_degrees() * 3600.0
        );
        if sep > 0.1 * p.params.fov {
            return None;
        }
        // The search that found the field is the one to report.
        r.search_dist_deg = wcs.search_dist_deg;
        r.step_distances.clone_from(&wcs.step_distances);
        Some(r)
    }

    /// The catalogue search from `start` out to `radius`, trying other pixel
    /// scales as [`SolveRequest::scale_search`] asks.
    ///
    /// With the scale unknown, the ladder of [`UNKNOWN_STEPS`] round the assumed
    /// 1″/px is searched near the hint first ([`Self::scale_ladder`]), and then the
    /// search at 1″/px runs as it always did, out to `radius`. With a known scale
    /// and [`ScaleSearch::AlsoIfWrong`], a search that finds nothing is followed by
    /// the ladder of [`WRONG_STEPS`] round that scale.
    fn catalogue_search(
        &self,
        unbinned: &ImageBuffer,
        img: &ImageBuffer,
        start: (f64, f64),
        radius: f64,
    ) -> Result<WcsSolution> {
        let params = SolveParams {
            ra_hint: start.0,
            dec_hint: start.1,
            search_radius: radius,
            ..self.params.clone()
        };
        if !self.searches_scales() {
            return solve_image(img, &params);
        }
        if !self.scale_known {
            let (lo, hi) = UNKNOWN_STEPS;
            let hyps = ladder(self.arcsec_per_px, lo, hi, false);
            if let Some(w) = self.scale_ladder(unbinned, img, start, radius, &hyps)? {
                return Ok(w);
            }
            log::info!(
                "No scale solved near the hint; searching {:.1}° at {:.2}\"/px.",
                radius.to_degrees(),
                self.arcsec_per_px
            );
            return solve_image(img, &params);
        }
        match solve_image(img, &params) {
            Err(e @ ArcsecError::InsufficientQuads { .. }) => {
                let hyps = ladder(self.arcsec_per_px, -WRONG_STEPS, WRONG_STEPS, true);
                log::info!(
                    "No solution at {:.3}\"/px; trying other scales round the hint.",
                    self.arcsec_per_px
                );
                self.scale_ladder(unbinned, img, start, radius, &hyps)?
                    .ok_or(e)
            }
            r => r,
        }
    }

    /// Search each scale hypothesis in `hyps` (most likely first) within
    /// [`LADDER_FIELDS`] fields of `start`, at most `radius`, and return the
    /// solution of the first that verifies.
    ///
    /// Hypotheses whose field no star database covers (the named one, or any
    /// installed) to within a factor 2 are left out, and so is one
    /// whose binned image would be too small to detect stars in. They are run on
    /// the solve's threads, in order: the first alone with all of them, the rest
    /// one thread each, and none starts once one has solved; the earliest that
    /// verifies wins, so the answer does not depend on the thread count (one
    /// after it that is still running is cancelled). Each is
    /// a [`ScaleTrust::Hypothesis`] search: at least 30 verified stars, and no
    /// catalogue-seeded fallback.
    ///
    /// The cost of a search that finds nothing is bounded by the ladder: at most
    /// one star detection per binning and minimum star size, and nine catalogue
    /// positions per hypothesis.
    fn scale_ladder(
        &self,
        unbinned: &ImageBuffer,
        img: &ImageBuffer,
        start: (f64, f64),
        radius: f64,
        hyps: &[Hypothesis],
    ) -> Result<Option<WcsSolution>> {
        let (w, h) = self.image_size;
        let naxis = w.max(h) as f64;
        let installed = available_dbs(&self.params.db_path);
        let covers = |fov_deg: f64| {
            let fits = |name: &str| {
                DB_FOV_RANGES.iter().find(|r| r.0 == name).map_or(
                    // A database arcsec has no range for: let it try.
                    self.request.db_name.is_some(),
                    |&(_, lo, hi, _)| fov_deg >= lo / 2.0 && fov_deg <= hi * 2.0,
                )
            };
            match &self.request.db_name {
                Some(name) => fits(name),
                None => installed.iter().any(|d| fits(d)),
            }
        };
        let plans: Vec<(Hypothesis, Plan)> = hyps
            .iter()
            .filter_map(|&hy| {
                let fov = (naxis * hy.scale / 3600.0).to_radians();
                if !covers(fov.to_degrees()) {
                    return None;
                }
                let req = SolveRequest {
                    hint: Some(start),
                    fov_height: None,
                    pixel_scale: Some(hy.scale),
                    scale_search: ScaleSearch::Never,
                    search_radius: radius.min(LADDER_FIELDS * fov),
                    index: None,
                    auto_index: false,
                    sip: false,
                    cancel: None,
                    ..self.request.clone()
                };
                let p = Plan::new(&req, w, h).ok()?;
                let (bw, bh) = p.binned_size();
                (bw >= MIN_SOLVE_DIM && bh >= MIN_SOLVE_DIM).then_some((hy, p))
            })
            .collect();
        if plans.is_empty() {
            return Ok(None);
        }
        log::info!(
            "Trying {} pixel scales, {:.3}–{:.3}\"/px, {} field round the hint each.",
            plans.len(),
            plans
                .iter()
                .map(|(h, _)| h.scale)
                .fold(f64::INFINITY, f64::min),
            plans.iter().map(|(h, _)| h.scale).fold(0.0, f64::max),
            LADDER_FIELDS
        );

        // Each binning the ladder needs, made once.
        let mut binned: Vec<(usize, Cow<'_, ImageBuffer>)> =
            vec![(self.binning, Cow::Borrowed(img))];
        for (_, p) in &plans {
            if binned.iter().all(|(b, _)| *b != p.binning) {
                let b = if p.binning > 1 {
                    Cow::Owned(unbinned.bin_image(p.binning))
                } else {
                    Cow::Borrowed(unbinned)
                };
                binned.push((p.binning, b));
            }
        }
        let image_for = |b: usize| {
            binned.iter().find(|(bb, _)| *bb == b).map_or_else(
                || unreachable!("binning {b} was made above"),
                |(_, i)| i.as_ref(),
            )
        };

        let n_threads = if self.params.threads > 0 {
            self.params.threads
        } else {
            crate::max_threads()
        }
        .clamp(1, 64);
        // Stars are detected once per binning and minimum star size, and shared by
        // the hypotheses with both: the minimum size bottoms out at 0.8 px, so all
        // but the finest few scales at a binning share one detection. The first
        // hypothesis detects with every thread; the rest share them out.
        let key = |p: &Plan| (p.binning, p.params.hfd_min.to_bits());
        let mut keys: Vec<(usize, u64)> = plans.iter().map(|(_, p)| key(p)).collect();
        keys.sort_unstable();
        keys.dedup();
        let detections: Vec<std::sync::OnceLock<Detected>> =
            keys.iter().map(|_| std::sync::OnceLock::new()).collect();
        let detect_threads = (n_threads / keys.len()).max(1);

        let cancel = crate::cancel::current();
        // The lowest hypothesis that has solved. A later one still running can no
        // longer win, so it is stopped at its next checkpoint rather than finished.
        let found = alloc::sync::Arc::new(core::sync::atomic::AtomicUsize::new(usize::MAX));
        let (_, winner) = search_in_order(plans.len(), n_threads, |i| {
            use core::sync::atomic::Ordering::Relaxed;
            let (hy, p) = &plans[i];
            if crate::cancel::fired(cancel.as_ref()) || found.load(Relaxed) < i {
                return (None, None);
            }
            let token = {
                let (found, outer) = (alloc::sync::Arc::clone(&found), cancel.clone());
                CancelToken::with_poll(move || {
                    found.load(Relaxed) < i || crate::cancel::fired(outer.as_ref())
                })
            };
            let threads = if i == 0 { n_threads } else { 1 };
            let params = SolveParams {
                threads,
                ..p.params.clone()
            };
            log::info!(
                "Scale hypothesis {:.3}\"/px ({:.2}° field, binning {}, star database {})",
                hy.scale,
                p.params.fov.to_degrees(),
                p.binning,
                p.params.db_name.to_uppercase()
            );
            let img = image_for(p.binning);
            let k = keys
                .binary_search(&key(p))
                .unwrap_or_else(|_| unreachable!());
            let solved = crate::cancel::with_token(&token, || {
                let stars = detections[k].get_or_init(|| {
                    let t = if i == 0 { n_threads } else { detect_threads };
                    crate::with_max_threads(t, || detect(img, &params))
                });
                crate::with_max_threads(threads, || {
                    solve_detected(img, &params, ScaleTrust::Hypothesis, stars)
                })
            });
            if solved.is_ok() {
                found.fetch_min(i, Relaxed);
            }
            (None, solved.ok().map(|w| (hy.scale, w)))
        });
        if crate::cancel::fired(cancel.as_ref()) {
            return Err(ArcsecError::Cancelled);
        }
        Ok(winner.map(|(_, (scale, wcs))| {
            log::info!("Solved at the scale hypothesis {scale:.3}\"/px.");
            wcs
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn catalog_dir_is_namespaced() {
        // Deliberately not asserting the exact path: it is platform dependent.
        let d = default_catalog_dir();
        assert!(
            d.to_string_lossy().contains("arcsec"),
            "catalogue dir should be namespaced: {}",
            d.display()
        );
    }

    #[test]
    fn env_override_wins() {
        let env = |k: &str| match k {
            "ARCSEC_CATALOG_DIR" => Some("/data/catalogs".to_string()),
            "HOME" => Some("/home/u".to_string()),
            _ => None,
        };
        assert_eq!(
            default_catalog_dir_from(env),
            PathBuf::from("/data/catalogs")
        );
    }

    #[test]
    fn empty_variables_count_as_unset() {
        let env = |k: &str| match k {
            "ARCSEC_CATALOG_DIR" | "XDG_DATA_HOME" | "LOCALAPPDATA" => Some(String::new()),
            "HOME" => Some("/home/u".to_string()),
            _ => None,
        };
        let d = default_catalog_dir_from(env);
        assert!(d.starts_with("/home/u"), "got {}", d.display());
        assert!(d.ends_with("catalogs"));
    }

    #[test]
    fn no_home_at_all_still_yields_a_path() {
        assert_eq!(
            default_catalog_dir_from(|_| None),
            PathBuf::from("catalogs")
        );
    }

    #[test]
    fn binning_is_automatic_below_one_arcsec_per_pixel() {
        assert_eq!(choose_binning(None, 2.0, 4000, 3000), 1);
        assert_eq!(choose_binning(Some(0), 0.5, 4000, 3000), 2);
        assert_eq!(choose_binning(None, 0.01, 4000, 3000), 16);
        assert_eq!(choose_binning(Some(3), 2.0, 4000, 3000), 3);
    }

    #[test]
    fn binning_never_exceeds_the_image() {
        assert_eq!(choose_binning(Some(100), 1.0, 4, 4), 2);
        assert_eq!(choose_binning(Some(usize::MAX), 1.0, 4, 4), 2);
        assert_eq!(choose_binning(None, 0.01, 10, 50), 5);
        assert_eq!(choose_binning(Some(4), 1.0, 1, 1), 1);
    }

    #[test]
    fn a_plan_follows_the_cli_rules() {
        let dir = crate::test_support::TempDir::new("auto_plan");
        let req = SolveRequest {
            hint: Some((1.0, 0.5)),
            fov_height: Some(1.0_f64.to_radians()),
            search_radius: -3.0,
            db_path: Some(dir.path().to_path_buf()),
            ..SolveRequest::default()
        };
        let p = Plan::new(&req, 4000, 2000).unwrap();
        assert_eq!(p.start, (1.0, 0.5));
        assert!(p.has_hint && p.scale_known);
        // Field along the longer side; scale from the height.
        assert!((p.params.fov.to_degrees() - 2.0).abs() < 1e-12);
        assert!((p.arcsec_per_px - 1.8).abs() < 1e-12);
        assert!((p.fov_height.to_degrees() - 1.0).abs() < 1e-12);
        assert_eq!(p.binning, 1);
        assert_eq!(p.params.search_radius, 0.0, "a negative radius is zero");
        assert_eq!(p.params.db_name, "d80", "nothing installed: the fallback");

        // No scale at all: 1"/px, unknown; fine sampling bins.
        let p = Plan::new(
            &SolveRequest {
                pixel_scale: Some(0.25),
                db_path: Some(dir.path().to_path_buf()),
                ..SolveRequest::default()
            },
            1000,
            1000,
        )
        .unwrap();
        assert_eq!(p.binning, 4);
        assert_eq!(p.binned_size(), (250, 250));
        assert_eq!(p.start, (0.0, 0.0));
        assert!(!p.has_hint);

        assert!(Plan::new(&SolveRequest::default(), 0, 10).is_err());
        let bad = SolveRequest {
            pixel_scale: Some(f64::NAN),
            ..SolveRequest::default()
        };
        assert!(Plan::new(&bad, 10, 10).is_err());
    }

    #[test]
    fn a_plan_refuses_an_image_of_another_size() {
        let dir = crate::test_support::TempDir::new("auto_size");
        let req = SolveRequest {
            db_path: Some(dir.path().to_path_buf()),
            ..SolveRequest::default()
        };
        let p = Plan::new(&req, 100, 80).unwrap();
        let err = p.solve(&ImageBuffer::new(80, 100)).unwrap_err();
        assert!(matches!(err, ArcsecError::InvalidParameter(_)), "{err}");
    }

    #[test]
    fn index_files_lists_only_arcsec_indexes() {
        let dir = crate::test_support::TempDir::new("auto_ix");
        for name in [
            "b.arcsecix",
            "a.arcsecix",
            "index-4107.fits",
            "d50_0101.1476",
        ] {
            std::fs::write(dir.path().join(name), b"x").unwrap();
        }
        let names: Vec<_> = index_files(dir.path())
            .iter()
            .map(|p| p.file_name().unwrap().to_string_lossy().into_owned())
            .collect();
        assert_eq!(names, ["a.arcsecix", "b.arcsecix"]);
        assert!(index_files(Path::new("/nonexistent/arcsec")).is_empty());
        assert!(has_star_database(dir.path()));
        assert_eq!(available_dbs(dir.path()), ["d50"]);
    }
}
