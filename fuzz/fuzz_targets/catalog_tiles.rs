//! ASTAP star-database tiles: `.1476` and `.290` packed records and the `.001`
//! all-sky file. The first byte picks the layout; the rest is the tile. The tile is
//! written as the database's first area (`_0101`), which covers the south pole in
//! both tilings, and read the ways the solver and the index builder read it.
#![no_main]

use arcsec_core::catalog::{for_each_star_in_dec_band, read_area_file, read_catalog_stars};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Some((&which, tile)) = data.split_first() else {
        return;
    };
    let (db, ext) = match which % 3 {
        0 => ("d50", "1476"),
        1 => ("g05", "290"),
        _ => ("w08", "001"),
    };
    // One directory per layout: the layout is detected from which files exist.
    let dir = arcsec_fuzz::scratch_dir().join(ext);
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join(format!("{db}_0101.{ext}"));
    std::fs::write(&path, tile).unwrap();

    let deg = f64::to_radians;
    if ext == "001" {
        let _ = arcsec_core::catalog::read_001_file(&path, 1.0, deg(-60.0), deg(30.0), 0.5, 500);
    } else {
        // The whole tile, then a small window.
        let _ = read_area_file(&path, 0.0, 0.0, 4.0 * core::f64::consts::PI, 1.0, 10_000);
        let _ = read_area_file(&path, 1.0, deg(-89.0), deg(1.0), deg(-89.0).cos(), 50);
    }
    for (ra, dec, fov) in [(0.0, -89.0, 2.0), (3.0, -85.0, 5.0), (1.0, -80.0, 15.0)] {
        if let Ok(stars) = read_catalog_stars(&dir, db, ra, deg(dec), deg(fov), 400) {
            assert!(stars.len() <= 400);
            assert!(stars.windows(2).all(|w| w[0].mag <= w[1].mag) || ext == "001");
        }
    }
    let _ = for_each_star_in_dec_band(&dir, db, deg(-90.0), deg(-75.0), 12.0, |_| {});
});
