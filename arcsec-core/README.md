# arcsec-core

[![crates.io](https://img.shields.io/crates/v/arcsec-core.svg)](https://crates.io/crates/arcsec-core)
[![docs.rs](https://docs.rs/arcsec-core/badge.svg)](https://docs.rs/arcsec-core)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/cruzzil/arcsec/blob/main/LICENSE)

The solving library behind [arcsec](https://github.com/cruzzil/arcsec), an astrometric
plate solver: given the pixels of an astronomical image and a rough idea of where it
points, it finds the exact position, scale and rotation as a WCS solution.

**Most people want the command-line tool instead**, which reads FITS, XISF and ASDF
files, installs star catalogues, writes ASTAP-compatible output files and is a drop-in
replacement for `astap_cli`:

```bash
cargo install arcsec
```

This crate is for Rust programs that want to solve in-process. It does not read image
files: you supply the pixels.

## What is in it

- `pipeline::solve_image` — the catalogue solve. Detects stars, builds ASTAP-style
  four-star quads, and searches outwards from the hint position in a spiral, matching
  against an ASTAP star database (`.1476`, `.290` or `.001` files). Returns a
  `WcsSolution`.
- `blind_solve` with `load_anet_index` — a position estimate with no hint, from an
  Astrometry.net index file. Returns an approximate RA/Dec to pass to `solve_image` as
  its hint.
- `catalog` — readers for the ASTAP star databases (`read_catalog_stars`,
  `detect_layout`) and the Astrometry.net index format (`load_anet_index`,
  `peek_anet_scale`).
- `types` — `ImageBuffer` (row-major `f32` greyscale pixels), `WcsSolution`, and the
  star and quad types used between stages.
- `set_max_threads` — one process-wide limit on worker threads for every stage;
  `1` makes a solve single-threaded.

The lower-level modules (`detection`, `quads`, `math`, `wcs`) are public too, but are
shaped by the pipeline's needs rather than designed as a general-purpose API.

**Angles are radians** in the parameters and in the solution's centre (`ra0`,
`dec0`); only the FITS-style terms of `WcsSolution` (CD matrix, `CDELT`, `CROTA2`) are
in degrees, as they are written to a header. Star databases come from ASTAP; the arcsec CLI can install
them (`arcsec catalog install d50`), or use an existing ASTAP directory.

## Example

```rust
use arcsec_core::ImageBuffer;
use arcsec_core::pipeline::{SolveMethod, SolveParams, solve_image};

fn main() -> arcsec_core::Result<()> {
    // Greyscale pixels, row by row from the top-left: data[y * width + x].
    let (data, width, height) = load_pixels();
    let mut img = ImageBuffer { data, width, height };
    // Rescales float data in physical units into the range detection expects;
    // leaves 16-bit camera data alone.
    img.normalize_for_detection();

    let params = SolveParams {
        ra_hint: 83.8_f64.to_radians(),
        dec_hint: (-5.4_f64).to_radians(),
        fov: 1.5_f64.to_radians(),          // image height
        search_radius: 10.0_f64.to_radians(),
        quad_tolerance: 0.007,
        hfd_min: 1.5,                       // minimum star size, in pixels
        max_stars: 500,
        db_path: "/path/to/star_databases".into(),
        db_name: "d50".into(),
        binning: 1,
        method: SolveMethod::Quads,
        threads: 0,                         // one per core
    };

    let wcs = solve_image(&img, &params)?;
    println!(
        "centre RA {:.5}°, Dec {:.5}°, {:.3}\"/px",
        wcs.ra0.to_degrees(),
        wcs.dec0.to_degrees(),
        wcs.cdelt2.abs() * 3600.0,
    );
    Ok(())
}

fn load_pixels() -> (Vec<f32>, usize, usize) {
    todo!("read the image with a FITS or other image library")
}
```

The `WcsSolution` carries the reference point (`ra0`, `dec0`, `crpix1`, `crpix2`), the
CD matrix, `cdelt1`/`cdelt2` and `crota2`, in the original image's pixels even when the
image was binned, plus the fit statistics. Errors are `ArcsecError`; the CLI maps them to ASTAP's exit codes.

## Status

The API follows the needs of the arcsec CLI and may change between 0.x releases; see
the [changelog](https://github.com/cruzzil/arcsec/blob/main/CHANGELOG.md). How the
solver works is described in
[docs/plate-solving.md](https://github.com/cruzzil/arcsec/blob/main/docs/plate-solving.md).

The approach is ASTAP's, by Han Kleijn; see the
[arcsec README](https://github.com/cruzzil/arcsec#credit) for credit.

## License

MIT.
