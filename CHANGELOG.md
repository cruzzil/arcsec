# Changelog

All notable changes to arcsec are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html). The `arcsec` and
`arcsec-core` crates share one version number.

The ASTAP-compatible command line - flags, stdout format, output files and exit codes -
is part of the public interface: a change to any of them is listed here, and a breaking
one is called out as such.

## [Unreleased]

## [0.1.0] - 2026-09-26

First public release.

### Added

- **Catalogue spiral solve**, the default mode. Detects stars, describes them as
  ASTAP-style four-star quads (five normalised distance ratios), and matches them
  against an ASTAP star database while spiralling outwards from an approximate
  position taken from the command line (`--ra`, `--spd`) or the image header. Matches
  are filtered by a scale-and-rotation vote and checked star by star before the
  plate constants are fitted by least squares.
- **Blind solve** with `-i/--index`: Astrometry.net index files (a single file or a
  directory) supply a position estimate with no hint at all, which the catalogue solver
  then refines. The index files best matched to the image's field of view are chosen
  automatically, and up to two are searched in parallel.
- `--method tetra` (experimental): an alternative matcher using three-star triangles
  instead of quads. It currently solves far fewer images than the default (20 of the 98
  solvable benchmark images, against 90) and is not recommended yet.
- **Input formats**: FITS, XISF (PixInsight) and ASDF (Roman/astropy), detected from
  the file's contents rather than its extension. Pointing and pixel scale are taken
  from the image's metadata where present and interpreted the same way for every
  format.
- **ASTAP compatibility**: the same flags (`-f`, `-r`, `--fov`, `--ra`, `--spd`, `-s`,
  `-t`, `-m`, `-z`, `-d`, `-D`, `-o`, `--wcs`, `--log`, `--update`, `--progress`), the
  same stdout report, the same `.wcs` and `.ini` output files, and the same exit codes:
  0 solved, 1 no solution, 2 too few stars, 16 file error, 32 database not found,
  33 database read error; a command-line usage error exits 1 rather than clap's usual 2,
  which ASTAP uses for "too few stars". As with ASTAP, a failed solve still writes an `.ini` holding
  `PLTSOLVD=F`, for tools that poll that file. ASTAP's single-dash spellings of the
  long options (`-fov`, `-ra`, `-spd`, ...) are accepted too. ASTAP options that are not implemented (`--sip`, `--check`,
  `--analyse`, `--extract`, `--extract2`, `--speed slow`) are refused with an error
  rather than silently ignored.
- **Star databases**: ASTAP's `.1476` files (D80, D50, D20, D05, V50), the
  `.290` files (G05, and V05) and the all-sky `.001` file (W08), covering fields from
  about 0.15° to 80°.
- **Automatic database selection**: without `-D`, the densest installed database whose
  field-of-view range contains the image is used.
- **`arcsec catalog`** subcommand with `list`, `recommend`, `install`, `remove`,
  `verify` and `path`. Installs ASTAP databases and Astrometry.net index sets
  (`anet-4100`, `anet-5200`) into a per-platform directory that the solver searches by
  default, overridable with `--dir` or `ARCSEC_CATALOG_DIR`.
- **Downsampling** with `-z`, including automatic selection with `-z 0`; the solution is
  reported in the original image's pixel coordinates.
- `--threads` to limit worker threads across detection, the search and blind solving;
  `--threads 1` runs single-threaded.
- `--update` writes the solution into the FITS header in place (FITS input only).
- Compressed FITS in bzip2 and Unix `compress` form is read directly. gzip-compressed
  FITS is not yet supported and is refused with exit 16.
- Builds with stable Rust 1.96 or later (the minimum supported Rust version).
- Pre-built binaries for Linux (x86-64, arm64), macOS (arm64) and Windows (x86-64).

### Performance

- On the 103-image benchmark corpus described in `docs/test-images.md`, arcsec solves
  90 images correctly with no false positives, against 47 for ASTAP CLI-2026.07.30 given
  the same hints.

[Unreleased]: https://github.com/cruzzil/arcsec/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/cruzzil/arcsec/releases/tag/v0.1.0
