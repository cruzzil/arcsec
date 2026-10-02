//! Astrometry.net index files (FITS binary tables read through CFITSIO's memory
//! driver): the header peek used to rank them, the full load, code lookups, and a
//! blind solve of a small synthetic field against whatever loaded.
#![no_main]

use arcsec_core::{BlindSolveParams, blind_solve, load_anet_index, peek_anet_scale};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    arcsec_core::set_max_threads(1);
    let path = arcsec_fuzz::write_input("index.fits", data);
    let _ = peek_anet_scale(&path);
    let Ok(index) = load_anet_index(&path) else {
        return;
    };
    assert_eq!(index.entries.len(), index.codes.len());
    for code in [[0.1, 0.2, 0.3, 0.4], [0.5; 4], [f64::NAN; 4]] {
        for i in index.find_code_matches(&code, 0.01) {
            assert!(i < index.entries.len());
        }
    }
    let params = BlindSolveParams {
        quad_tolerance: 0.008,
        hfd_min: 0.8,
        max_stars: 100,
        binning: 1,
        fov_deg: 0.5,
    };
    let _ = blind_solve(arcsec_fuzz::synthetic_field(), &index, &params);
});
