# arcsec

[![CI](https://github.com/cruzzil/arcsec/actions/workflows/ci.yml/badge.svg)](https://github.com/cruzzil/arcsec/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/arcsec.svg)](https://crates.io/crates/arcsec)
[![docs.rs](https://docs.rs/arcsec-core/badge.svg)](https://docs.rs/arcsec-core)
[![codecov](https://codecov.io/gh/cruzzil/arcsec/graph/badge.svg?token=dBUWxzz6KM)](https://codecov.io/gh/cruzzil/arcsec)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

An astrometric plate solver, written in Rust. Give it an astronomical image and it
works out where the telescope was pointing.

It works as a **drop-in replacement for ASTAP's command-line solver** (`astap_cli`):
same flags, same output files, same exit codes, same star databases. Imaging software
that drives ASTAP, such as [N.I.N.A.](https://nighttime-imaging.eu/), can use it
unchanged.

**[cruzzil.github.io/arcsec](https://cruzzil.github.io/arcsec/)** has downloads,
set-up guides and the full command-line reference.

## Features

- **Solves fields from 0.15° to 80°**, picking the right star database for the field.
- **Blind solving** with no position and no pixel scale, from an index built locally
  from the star database you install. Astrometry.net index files work too.
- **Catalogue management built in:** `arcsec catalog install` downloads a star database
  to a place the solver already looks.
- **Reads FITS** (including compressed), **XISF** (PixInsight) and **ASDF**
  (Roman/astropy).
- **Writes ASTAP's `.wcs` and `.ini` files**, and can update the FITS header in place,
  with optional SIP distortion terms.
- **Tested against hundreds of images with known answers** from public sky surveys,
  every solution checked at the centre and all four corners
  ([results](docs/test-images.md)).
- **Usable as a library** from Rust (`arcsec-core`) and from C or C++ (`libarcsec`).

## Install

Download a prebuilt binary for Linux (x86-64, arm64), macOS (Apple silicon) or Windows
(x86-64) from the [latest release](https://github.com/cruzzil/arcsec/releases/latest),
unpack it, and put `arcsec` on your `PATH`.

Or build it with Rust 1.96 or newer (a C compiler is also needed, for the `ring`
dependency used by HTTPS downloads):

```bash
cargo install arcsec
```

## Quick start

```bash
arcsec catalog recommend --fov 2.5   # which star database suits a 2.5° field?
arcsec catalog install d50           # download it, and build its blind index
arcsec -f image.fits                 # solve
```

arcsec takes the rough position and pixel scale from the image header (`RA`/`DEC`,
`FOCALLEN`, `XPIXSZ`). You can also give them yourself:

```bash
arcsec -f image.fits --ra 5.58 --spd 82.0 --fov 1.5 -r 10
```

`--ra` is in hours, `--spd` is 90 + Dec in degrees, `--fov` is the image height in
degrees and `-r` is the search radius in degrees. ASTAP's single-dash spellings (`-fov`,
`-ra`, `-spd`) work too. With no position at all, arcsec solves blind using the index
that `catalog install` built.

A successful solve writes `image.wcs` (the WCS as a FITS header) and `image.ini`
(ASTAP's summary), and `--update` also writes the solution into a FITS image's own
header. Exit codes
match ASTAP's: `0` solved, `1` no solution, `2` too few stars, `16` image error, `32`
database not found, `33` database read error.

See the [command-line reference](https://cruzzil.github.io/arcsec/reference/cli/) for
every option, star measurement without solving (`--analyse`, `--extract`), and SIP
distortion fitting.

## Using it with other software

- **N.I.N.A.:** choose ASTAP as the plate solver and type the full path to `arcsec.exe`
  as the ASTAP location. N.I.N.A.'s file browser only lists files named `astap.exe`, so
  type the path rather than browsing. [Set-up guide](https://cruzzil.github.io/arcsec/nina/)
- **Siril:** solve with `arcsec -f image.fit --update`, then open the file in Siril,
  which treats it as plate solved. [Guide](https://cruzzil.github.io/arcsec/siril/)
- **Anything that runs `astap_cli`:** point it at `arcsec` instead.

## Star databases

arcsec reads ASTAP's star databases, so an existing ASTAP install can share its
databases with arcsec: set `ARCSEC_CATALOG_DIR` (or pass `-d`) to ASTAP's folder. The
[catalogue guide](https://cruzzil.github.io/arcsec/which-catalogue/) explains which
database suits which field size.

## Libraries

| Crate | For |
|---|---|
| [`arcsec-core`](https://crates.io/crates/arcsec-core) | The solver, for Rust programs ([README](arcsec-core/README.md)) |
| [`arcsec-io`](https://crates.io/crates/arcsec-io) | Reading FITS, XISF and ASDF images; writing `.wcs`/`.ini` |
| [`arcsec-catalogue`](https://crates.io/crates/arcsec-catalogue) | Installing and managing star databases and the blind index ([README](arcsec-catalogue/README.md)) |
| [`libarcsec`](libarcsec/README.md) | A C library and header, also shipped prebuilt with each release |

## Documentation

- [docs/catalogues.md](docs/catalogues.md): catalogues, field sizes and install locations
- [docs/plate-solving.md](docs/plate-solving.md): how plate solving works, and exactly
  what arcsec does
- [docs/test-images.md](docs/test-images.md): the benchmark corpus and results
- [docs/offline-index.md](docs/offline-index.md): the blind index
- [CHANGELOG.md](CHANGELOG.md) and [CONTRIBUTING.md](CONTRIBUTING.md)

## Credit

arcsec is built on the star-pattern approach of **[ASTAP](https://www.hnsky.org/astap.htm)**
by Han Kleijn, who
[documents the algorithm openly](https://www.hnsky.org/astap_astrometric_solving.htm),
and it reads ASTAP's star database formats. arcsec is an independent Rust
implementation and contains no ASTAP source; ASTAP itself is licensed under the MPL 2.0.

## License

MIT, see [LICENSE](LICENSE).
