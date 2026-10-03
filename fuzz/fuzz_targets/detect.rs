//! Pixel data whatever its values: what reaches the solver once a reader has
//! accepted a file. NaN, infinities, constant frames, extreme ranges and tiny or
//! oddly shaped images go through the same steps as in `main`: normalisation, the
//! check-pattern filter, binning, `--analyse`, and a hinted solve (against an empty
//! catalogue, so it runs detection and pattern building and then finds nothing).
//!
//! Input: width and height (one byte each, plus one), a sample format byte, then
//! samples, repeated to fill the frame.
#![no_main]

use std::path::PathBuf;

use arcsec_core::detection::analyse_image;
use arcsec_core::pipeline::{SearchSpeed, SolveMethod, SolveParams, solve_image};
use arcsec_core::types::ImageBuffer;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let [w, h, fmt, rest @ ..] = data else {
        return;
    };
    let (width, height) = (usize::from(*w) + 1, usize::from(*h) + 1);
    let samples: Vec<f32> = match fmt % 3 {
        0 => rest
            .as_chunks::<4>()
            .0
            .iter()
            .map(|b| f32::from_le_bytes(*b))
            .collect(),
        1 => rest
            .as_chunks::<2>()
            .0
            .iter()
            .map(|b| f32::from(u16::from_le_bytes(*b)))
            .collect(),
        _ => rest.iter().map(|&b| f32::from(b)).collect(),
    };
    if samples.is_empty() {
        return;
    }
    let data: Vec<f32> = samples
        .iter()
        .copied()
        .cycle()
        .take(width * height)
        .collect();
    let mut img = ImageBuffer {
        data,
        width,
        height,
    };
    arcsec_core::set_max_threads(1);
    img.normalize_for_detection();
    let mut filtered = img.clone();
    filtered.check_pattern_filter();
    let binned = img.bin_image(usize::from(fmt / 3 % 4) + 1);
    let _ = analyse_image(&img, 10.0, 200);

    // An empty catalogue: present, so the solve proceeds, but every tile is short.
    let db = arcsec_fuzz::scratch_dir().join("emptydb");
    std::fs::create_dir_all(&db).unwrap();
    let tile: PathBuf = db.join("d50_0101.1476");
    if !tile.exists() {
        std::fs::write(&tile, b"").unwrap();
    }
    let params = SolveParams {
        ra_hint: 1.0,
        dec_hint: 0.3,
        fov: 0.01,
        search_radius: 0.0,
        quad_tolerance: 0.007,
        hfd_min: 0.8,
        max_stars: 300,
        db_path: db,
        db_name: "d50".to_string(),
        binning: 1,
        method: if fmt & 0x80 == 0 {
            SolveMethod::Quads
        } else {
            SolveMethod::Tetra
        },
        speed: SearchSpeed::Auto,
        threads: 1,
    };
    // main refuses anything smaller before solving (MIN_SOLVE_DIM).
    if binned.width >= 2 && binned.height >= 2 {
        let _ = solve_image(&binned, &params);
    }
});
