# arcsec

[![CI](https://github.com/cruzzil/arcsec/actions/workflows/ci.yml/badge.svg)](https://github.com/cruzzil/arcsec/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/cruzzil/arcsec/graph/badge.svg?token=MjXKzC5keQ)](https://codecov.io/gh/cruzzil/arcsec)
[![Dependency status](https://deps.rs/repo/github/cruzzil/arcsec/status.svg)](https://deps.rs/repo/github/cruzzil/arcsec)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

An astrometric plate solver written in Rust. Give it an astronomical image and it
works out where the telescope was pointing, writing a WCS solution.

```bash
arcsec -f image.fits
```

Reads **FITS**, **XISF** (PixInsight) and **ASDF** (Roman/astropy) — the format is
detected from the file's magic bytes, not its extension.

## Credit

arcsec would not exist without **[ASTAP](https://www.hnsky.org/astap.htm)** by
Han Kleijn. Two debts in particular:

- **The algorithm.** ASTAP's star-pattern approach — describing a quad of four stars
  by five normalised distance ratios, which are invariant under rotation, scaling and
  flipping — is the idea arcsec is built on. Kleijn documents it openly at
  [*ASTAP star pattern recognition algorithm and astrometric (plate) solving*](https://www.hnsky.org/astap_astrometric_solving.htm).
- **The catalogue file formats.** arcsec reads ASTAP's `.1476`, `.290` and `.001` star
  database files, so the same catalogues serve both programs. The `.290` and `.001`
  layouts are not documented upstream and were reverse engineered from the shipped
  files.

arcsec also mirrors ASTAP's command-line flags, stdout format, output files and exit
codes, so it can be dropped into an existing workflow in place of `astap_cli`.

ASTAP itself is licensed under the Mozilla Public License 2.0. arcsec is an
independent implementation in Rust and contains no ASTAP source.

## Installing catalogues

The solver needs a star database. `arcsec` can fetch one for you:

```bash
arcsec catalog recommend --fov 2.5    # which catalogue suits a 2.5° field?
arcsec catalog install d50            # install it
arcsec catalog list
```

Catalogues land in a per-platform directory that the solver searches by default, so
after an install you need neither `-d` nor `-D`. Override it with `--dir` or the
`ARCSEC_CATALOG_DIR` environment variable. See [docs/catalogues.md](docs/catalogues.md)
for which catalogue suits which field size.

## Solving

With a rough idea of where the image points (from the FITS header, or given explicitly):

```bash
arcsec -f image.fits --ra 5.58 --spd 82.0 --fov 1.5
```

Blind, with no hint at all, using Astrometry.net index files:

```bash
arcsec catalog install anet-4100
arcsec -f image.fits -i ~/.local/share/arcsec/catalogs/anet-4100
```

Useful flags: `--progress` logs to stderr, `--log` writes `<base>.log`, `--update`
writes the solution back into a FITS header, `--wcs` writes an Astrometry.net-style
`.wcs` file, and `--threads 1` gives a genuinely single-threaded run.

Exit codes follow ASTAP: `0` solved, `1` no solution, `2` too few stars, `16` file
error, `32` database not found, `33` database read error.

## Documentation

- [docs/catalogues.md](docs/catalogues.md) — which catalogue for which field, and where
  they install.
- [docs/plate-solving.md](docs/plate-solving.md) — how plate solving works in general,
  and precisely what arcsec does, with flowcharts and the constants table.
- [docs/test-images.md](docs/test-images.md) — the benchmark corpus and measured results.
- [docs/offline-index.md](docs/offline-index.md) — design notes for a pre-computed quad
  index.

## Building

```bash
cargo build --release
cargo test
```

No external dependencies beyond the crate graph; the release profile keeps debug symbols
for profiling.

## License

MIT — see [LICENSE](LICENSE).
