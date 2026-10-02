//! arcsec's own blind index (`ARCSECIX`, memory-mapped): open, the full checksum
//! validation of `arcsec catalog verify`, every accessor, and an index solve of a
//! small synthetic field. The header CRC is recomputed over the mutated bytes, so
//! the fuzzer reaches the section checks and lookups rather than stopping at it.
#![no_main]

use std::path::PathBuf;

use arcsec_core::index::BlindIndex;
use arcsec_core::pipeline::{IndexSolveParams, SearchSpeed, SolveMethod, SolveParams, index_solve};
use libfuzzer_sys::fuzz_target;

fn crc32(bytes: &[u8]) -> u32 {
    let mut c = 0xFFFF_FFFFu32;
    for &b in bytes {
        c ^= u32::from(b);
        for _ in 0..8 {
            c = if c & 1 != 0 {
                0xEDB8_8320 ^ (c >> 1)
            } else {
                c >> 1
            };
        }
    }
    !c
}

fuzz_target!(|data: &[u8]| {
    arcsec_core::set_max_threads(1);
    let mut data = arcsec_fuzz::with_magic(data, b"ARCSECIX");
    // Most of the time, make the header checksum right so the body is reached.
    if data.len() >= 256 && data[255] & 1 == 0 {
        let crc = crc32(&data[..252]);
        data[252..256].copy_from_slice(&crc.to_le_bytes());
    }
    let path = arcsec_fuzz::write_input("index.arcsecix", &data);
    let Ok(ix) = BlindIndex::open(&path) else {
        return;
    };
    let _ = ix.validate();
    let _ = (
        ix.n_stars(),
        ix.n_patterns(),
        ix.source(),
        ix.built_unix(),
        ix.file_size(),
    );
    for t in ix.tiers() {
        let _ = ix.tier_bytes(t);
        for key in [0, 1, 0x1234_5678, u64::MAX] {
            for p in ix.lookup(t, key).take(64) {
                let _ = ix.quad(p);
            }
        }
    }
    for i in [
        0,
        1,
        ix.n_stars().saturating_sub(1),
        ix.n_stars(),
        usize::MAX,
    ] {
        let _ = ix.star(i);
    }
    for i in [
        0,
        ix.n_patterns().saturating_sub(1),
        ix.n_patterns(),
        usize::MAX / 16,
    ] {
        let _ = ix.quad(i);
    }
    for (ra, dec, r) in [
        (0.0, 0.0, 0.01),
        (6.0, -1.5, 0.2),
        (3.0, 1.57, 3.2),
        (1.0, 0.3, f64::NAN),
    ] {
        ix.stars_near(ra, dec, r, |_| {});
    }

    // The solve re-detects the synthetic field every time, which costs a hundred
    // times what the index work does; run it for one input in eight.
    if data.len() < 256 || data[254] & 7 != 0 {
        return;
    }
    let template = SolveParams {
        ra_hint: 0.0,
        dec_hint: 0.0,
        fov: 0.01,
        search_radius: 0.0,
        quad_tolerance: 0.007,
        hfd_min: 0.8,
        max_stars: 100,
        db_path: PathBuf::from("/nonexistent/arcsec-fuzz"),
        db_name: "d50".to_string(),
        binning: 1,
        method: SolveMethod::Quads,
        speed: SearchSpeed::Auto,
        threads: 1,
    };
    let params = IndexSolveParams {
        scale_lo: 0.5,
        scale_hi: 20.0,
        within: None,
    };
    let _ = index_solve(arcsec_fuzz::synthetic_field(), &ix, &template, &params);
});
