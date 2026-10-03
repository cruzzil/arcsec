//! Shared support for the fuzz targets.
//!
//! The image readers come from the `arcsec-io` crate. The archive extraction lives
//! in the `arcsec` binary crate, which has no library target to depend on, so its
//! source file is compiled into this crate directly.

// As in arcsec's main.rs: `alloc` must be declared before `alloc::` paths resolve.
extern crate alloc;

pub use arcsec_io::{asdf_io, fits_io, image_io, xisf_io};

/// The CLI's `catalog_cmd/fetch.rs`, which refers to its parent as `super`; this
/// crate's root stands in for that parent.
#[allow(dead_code)]
#[path = "../../arcsec/src/catalog_cmd/fetch.rs"]
pub mod fetch;

/// Stand-in for the CLI's byte-count formatter, which only progress output uses.
#[must_use]
pub fn human(bytes: u64) -> String {
    format!("{bytes} B")
}

use std::path::{Path, PathBuf};
use std::sync::OnceLock;

/// A directory private to this process, on a RAM disk where there is one.
///
/// Every reader under test takes a path, so each input is written to a file first.
/// The process id keeps parallel fuzzing jobs (`-jobs`, `-fork`) apart.
pub fn scratch_dir() -> &'static Path {
    static DIR: OnceLock<PathBuf> = OnceLock::new();
    DIR.get_or_init(|| {
        let shm = Path::new("/dev/shm");
        let base = if shm.is_dir() {
            shm.to_path_buf()
        } else {
            std::env::temp_dir()
        };
        let dir = base.join(format!("arcsec-fuzz-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("create the fuzz scratch directory");
        dir
    })
}

/// Write `data` to `name` in the scratch directory and return its path.
pub fn write_input(name: &str, data: &[u8]) -> PathBuf {
    let path = scratch_dir().join(name);
    std::fs::write(&path, data).expect("write the fuzz input");
    path
}

/// Overwrite the start of `data` with `magic`, so a target spends its time inside
/// one reader rather than in format detection's rejection path.
#[must_use]
pub fn with_magic(data: &[u8], magic: &[u8]) -> Vec<u8> {
    let mut v = data.to_vec();
    if v.len() < magic.len() {
        v.resize(magic.len(), 0);
    }
    v[..magic.len()].copy_from_slice(magic);
    v
}

/// Everything the CLI does with an input image before solving it: format
/// detection, the pixel read, every metadata read, and the detection normalisation.
/// With `update`, also the `--update` header write that follows a solve.
pub fn exercise_image(data: &[u8], name: &str, update: bool) {
    let path = write_input(name, data);
    let _ = image_io::detect_format(&path);
    let _ = image_io::read_ra_dec(&path);
    let _ = image_io::read_pixel_scale(&path);
    let _ = image_io::read_dimensions(&path);
    let _ = image_io::read_header_wcs(&path);
    let _ = image_io::read_channels(&path);
    if let Ok(mut img) = image_io::read_image(&path) {
        assert_eq!(img.data.len(), img.width * img.height, "pixel count");
        img.normalize_for_detection();
        if update {
            let _ = image_io::update_wcs(&path, &solution());
        }
    }
}

/// A small synthetic star field (256 × 256, 40 stars, sky 1000 ADU with noise),
/// for the targets that drive a solver with a fuzzed index. Built once.
pub fn synthetic_field() -> &'static arcsec_core::types::ImageBuffer {
    static IMG: OnceLock<arcsec_core::types::ImageBuffer> = OnceLock::new();
    IMG.get_or_init(|| {
        let (w, h) = (256usize, 256usize);
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut data: Vec<f32> = (0..w * h).map(|_| 1000.0 + 20.0 * next() as f32).collect();
        for _ in 0..40 {
            let (cx, cy) = (8.0 + next() * 240.0, 8.0 + next() * 240.0);
            let peak = 2000.0 + next() * 30_000.0;
            for y in (cy as usize - 6)..(cy as usize + 6) {
                for x in (cx as usize - 6)..(cx as usize + 6) {
                    let r2 = (x as f64 - cx).powi(2) + (y as f64 - cy).powi(2);
                    data[y * w + x] += (peak * (-r2 / 4.5).exp()) as f32;
                }
            }
        }
        arcsec_core::types::ImageBuffer {
            data,
            width: w,
            height: h,
        }
    })
}

/// A plausible solution, for the `--update` path.
fn solution() -> arcsec_core::types::WcsSolution {
    arcsec_core::types::WcsSolution {
        ra0: 1.0,
        dec0: 0.5,
        crpix1: 100.5,
        crpix2: 80.5,
        cd1_1: -3.4e-4,
        cd1_2: 1e-8,
        cd2_1: -1e-8,
        cd2_2: 3.4e-4,
        cdelt1: -3.4e-4,
        cdelt2: 3.4e-4,
        crota2: 0.1,
        residual_rms: 0.5,
        stars_matched: 50,
        plate: arcsec_core::types::PlateConstants {
            a: 1.0,
            b: 0.0,
            c: 0.0,
            d: 0.0,
            e: 1.0,
            f: 0.0,
        },
        mag_limit: 15.0,
        search_dist_deg: 0.0,
        step_distances: Vec::new(),
        raw_matches: 50,
        matched_stars: Vec::new(),
        sip: None,
    }
}
