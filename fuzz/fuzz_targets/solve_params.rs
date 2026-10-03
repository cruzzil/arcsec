//! The numbers a file's header (RA, DEC, FOCALLEN, XPIXSZ ...) or the command line
//! hands the solver: any `f64`, NaN and infinities included. A small synthetic
//! field is solved against the star database in `$ARCSEC_FUZZ_DB` (an ASTAP
//! directory; default `~/star_database`, database `d80`), or against an empty one if
//! that is absent, with the hint, field size, radius, tolerance and limits taken
//! from the input. The search radius is limited to a few fields so that one input
//! cannot run for hours; that the spiral itself is bounded is a unit test.
#![no_main]

use std::path::PathBuf;

use arcsec_core::pipeline::{
    IndexSolveParams, SearchSpeed, SolveMethod, SolveParams, index_solve, solve_image,
};
use libfuzzer_sys::fuzz_target;

fn db() -> (PathBuf, String) {
    let dir = std::env::var_os("ARCSEC_FUZZ_DB").map_or_else(
        || PathBuf::from(std::env::var_os("HOME").unwrap_or_default()).join("star_database"),
        PathBuf::from,
    );
    if arcsec_core::catalog::catalog_present(&dir, "d80") {
        return (dir, "d80".into());
    }
    let empty = arcsec_fuzz::scratch_dir().join("emptydb");
    std::fs::create_dir_all(&empty).unwrap();
    let _ = std::fs::write(empty.join("d80_0101.1476"), b"");
    (empty, "d80".into())
}

fuzz_target!(|data: &[u8]| {
    let mut words = data
        .as_chunks::<8>()
        .0
        .iter()
        .map(|b| f64::from_le_bytes(*b));
    let mut next = || words.next().unwrap_or(0.0);
    let (ra, dec, fov, radius) = (next(), next(), next(), next());
    let (tol, hfd, scale_lo, scale_hi) = (next(), next(), next(), next());
    let flags = data.last().copied().unwrap_or(0);
    arcsec_core::set_max_threads(1);

    let (db_path, db_name) = db();
    let params = SolveParams {
        ra_hint: ra,
        dec_hint: dec,
        fov,
        // A few fields at most: NaN stays NaN, which the solver must refuse.
        search_radius: if radius.is_nan() {
            radius
        } else {
            radius.min(fov.abs() * 3.0)
        },
        quad_tolerance: tol,
        hfd_min: hfd,
        max_stars: usize::from(flags) * 4,
        db_path,
        db_name,
        binning: usize::from(flags >> 6) + 1,
        method: if flags & 1 == 0 {
            SolveMethod::Quads
        } else {
            SolveMethod::Tetra
        },
        speed: if flags & 2 == 0 {
            SearchSpeed::Auto
        } else {
            SearchSpeed::Slow
        },
        threads: 1,
    };
    let img = arcsec_fuzz::synthetic_field();
    let _ = solve_image(img, &params);

    if flags & 4 != 0 {
        // An index solve against an empty index: the parameter checks and the
        // hand-off to the hinted solver, without a real index.
        let path = arcsec_fuzz::scratch_dir().join("empty.arcsecix");
        if !path.exists() {
            arcsec_core::index::BuiltIndex {
                star_dir: vec![0; arcsec_core::index::format::STAR_BANDS as usize + 1],
                ..Default::default()
            }
            .write(&path)
            .unwrap();
        }
        if let Ok(ix) = arcsec_core::index::BlindIndex::open(&path) {
            let within = (flags & 8 != 0).then_some((ra, dec, radius));
            let p = IndexSolveParams {
                scale_lo,
                scale_hi,
                within,
            };
            let _ = index_solve(img, &ix, &params, &p);
        }
    }
});
