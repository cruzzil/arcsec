# arcsec

[![CI](https://github.com/cruzzil/arcsec/actions/workflows/ci.yml/badge.svg)](https://github.com/cruzzil/arcsec/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/arcsec.svg)](https://crates.io/crates/arcsec)
[![docs.rs](https://docs.rs/arcsec-core/badge.svg)](https://docs.rs/arcsec-core)
[![codecov](https://codecov.io/gh/cruzzil/arcsec/graph/badge.svg?token=dBUWxzz6KM)](https://codecov.io/gh/cruzzil/arcsec)
[![Dependency status](https://deps.rs/repo/github/cruzzil/arcsec/status.svg)](https://deps.rs/repo/github/cruzzil/arcsec)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

An astrometric plate solver written in Rust. Give it an astronomical image and it
works out where the telescope was pointing, writing a WCS solution.

**Website: [cruzzil.github.io/arcsec](https://cruzzil.github.io/arcsec/)** — downloads,
setting up N.I.N.A. and Siril, and [which catalogue you need](https://cruzzil.github.io/arcsec/which-catalogue/).

```bash
arcsec -f image.fits
```

- **A drop-in replacement for `astap_cli`**, so imaging software that drives ASTAP -
  [N.I.N.A.](https://nighttime-imaging.eu/) in particular - can use arcsec instead,
  unchanged. The same solving flags, stdout report, `.wcs` and `.ini` files and exit
  codes, and it reads ASTAP's star databases. ASTAP's single-dash spellings (`-fov`,
  `-ra`, `-spd`) work as well as `--fov`, `--ra`, `--spd`. See
  [Using arcsec from N.I.N.A.](#using-arcsec-from-nina)
- **Tested against hundreds of images with known answers,** from ten public sky
  surveys and telescope archives plus simulated camera faults, with every solution
  checked at the centre and all four corners. Methods and results are in
  [docs/test-images.md](docs/test-images.md).
- **Reads FITS, XISF (PixInsight) and ASDF (Roman/astropy).** The format is detected
  from the file's contents, not its extension.
- **Fields from 0.15° to 80°**, choosing the right database for the field size
  automatically.
- **Blind solving** with no position hint — and no pixel scale either — using an index
  that `arcsec catalog index build` makes from the star database you already have, in
  minutes, with nothing to download. Astrometry.net index files work too.
- **Catalogue management built in**: `arcsec catalog install d50` downloads and
  unpacks a star database into a directory the solver already knows about.

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

## Installing

**Prebuilt binaries** for Linux (x86-64, arm64), macOS (Apple silicon) and Windows
(x86-64) are attached to each [GitHub release](https://github.com/cruzzil/arcsec/releases),
with SHA-256 checksums. Unpack the archive and put `arcsec` (or `arcsec.exe`) on your
`PATH`. The Linux binaries are built on Ubuntu 22.04 and need glibc 2.35 or newer.

**From crates.io**, with Rust 1.96 or newer:

```bash
cargo install arcsec
```

**From source:**

```bash
git clone https://github.com/cruzzil/arcsec
cd arcsec
cargo build --release          # binary at target/release/arcsec
```

Building needs a C compiler (`cc` on Linux and macOS, the MSVC build tools on Windows)
for one dependency: [`ring`](https://crates.io/crates/ring), the cryptography behind the
HTTPS that `arcsec catalog install` downloads over, which compiles some C and assembly.
Nothing else is required: FITS support comes from
[rsfitsio](https://crates.io/crates/rsfitsio), a Rust port of CFITSIO, so no system
CFITSIO is needed. The `arcsec-core` library on its own is pure Rust.

## Quick start

The solver needs a star database. `arcsec` can fetch one for you:

```bash
arcsec catalog recommend --fov 2.5    # which catalogue suits a 2.5° field?
arcsec catalog install d05            # install what it suggested
arcsec -f image.fits                  # solve
```

Catalogues land in a per-platform directory that the solver searches by default, so
after an install you need neither `-d` nor `-D`. Override it with `--dir` or the
`ARCSEC_CATALOG_DIR` environment variable. If you already have ASTAP databases, point
`ARCSEC_CATALOG_DIR` (or `-d`) at them; both programs can share one directory. See
[docs/catalogues.md](docs/catalogues.md) for which catalogue suits which field size.

## Using arcsec from N.I.N.A.

N.I.N.A. runs ASTAP as a command-line program and reads the `.ini` file it writes, so
arcsec takes its place without any change on N.I.N.A.'s side:

1. Install a star database for arcsec, e.g. `arcsec catalog install d50` (run
   `arcsec catalog recommend --fov <your field height in degrees>` to choose). N.I.N.A.
   does not pass `-d`, so arcsec uses its own catalogue directory. To reuse the
   databases an existing ASTAP install already has instead, set `ARCSEC_CATALOG_DIR`
   to ASTAP's folder (by default `C:\Program Files\astap`).
2. In N.I.N.A., under **Options > Plate Solving**, choose **ASTAP** as the plate
   solver (and as the blind solver, if you like), and set **ASTAP location** to
   `arcsec.exe`. Type or paste the full path into the field, e.g.
   `C:\Users\<you>\.cargo\bin\arcsec.exe` after `cargo install`: N.I.N.A.'s file
   browser for this setting only shows files named `astap.exe`, so it cannot select
   `arcsec.exe`.

N.I.N.A.'s own settings - search radius, downsampling, maximum stars - are passed
through as the corresponding ASTAP options. Its `-fov` is the image height, which is
what arcsec takes it to be too.

## Using arcsec from Siril

[Siril](https://siril.org/) (tested with 1.4.4) cannot run ASTAP's command-line solver:
its plate solver is its own, or a local Astrometry.net `solve-field`, which is given a
star list rather than the image. Instead, solve the FITS file with arcsec and write the
solution into its header, then open it in Siril:

```bash
arcsec -f M101_stacked.fit --update
```

Siril reads the solution when it loads the file and treats the image as plate solved
(annotations, PCC/SPCC; `platesolve` reports it is already solved). To use Siril's
astrometric registration, solve each frame this way before building the sequence, and
let Siril's sequence plate solve skip the frames that are already solved. Siril needs
the `CD` or `PC` matrix, which arcsec writes; SIP terms from `--sip` are applied too.
See [Use with Siril](https://cruzzil.github.io/arcsec/siril/) for details.

## Solving

arcsec needs a rough position and the field size. It takes the position from the
`RA`/`DEC` header keywords in degrees (falling back to `CRVAL1`/`CRVAL2`), and the pixel
scale from `FOCALLEN`, `XPIXSZ` and `XBINNING`, as most capture software writes them.
XISF files carry the same keywords; ASDF files are read from their metadata tree.
Either can be given explicitly instead:

```bash
arcsec -f image.fits --ra 5.58 --spd 82.0 --fov 1.5 -r 10
```

`--ra` is in hours, `--spd` is south pole distance (90 + Dec) in degrees, `--fov` is the
image height in degrees, and `-r` is the search radius around that position in degrees
(default 180, the whole sky). Without `FOCALLEN`/`XPIXSZ` in the header, give `--fov`:
otherwise arcsec assumes 1″ per pixel.

Blind, with no position hint, using arcsec's own index:

```bash
arcsec catalog index build                       # once: ~2 min, ~290 MB from D80
arcsec -f image.fits -i "$(arcsec catalog path)"
```

The index is built from the installed star database (fields 0.3°–30° by default;
`--min-fov 0.15` for D80's narrowest), so nothing is downloaded. Every position it finds
is verified by the ordinary solver before it is reported. With an index installed,
searches of `-r` 10° or more that reach past five fields round the hint — N.I.N.A.'s
blind mode sends `-r 180` — use it automatically once the first five fields have failed. See
[docs/offline-index.md](docs/offline-index.md).

`-i` also takes Astrometry.net index files (`arcsec catalog install anet-4100`, fields of
about 0.7° and wider), or a directory of `index-*.fits`, and picks the ones whose scale
suits the field. Those only estimate the position, which is then refined against the
star database.

Other useful flags:

| Flag | Effect |
|---|---|
| `-o <base>` | Base path for the output files (default: the image path without its extension) |
| `-D <name>` | Force a database (`d80`, `d50`, `g05`, `w08`, ...) instead of choosing by field size |
| `-d <dir>` | Star database directory, for this run only |
| `-z <n>` | Bin the image n×n before solving; `0` or absent chooses automatically |
| `-s <n>` | Maximum number of stars to use (default 500) |
| `--update` | Write the solution into the FITS header in place (FITS only) |
| `--progress` | Log each step to stderr |
| `--log` | Write the same log to `<base>.log` |
| `--threads <n>` | Limit worker threads; `--threads 1` is genuinely single-threaded |
| `--sip` | Add SIP distortion terms to the solution, when the field shows distortion (below) |
| `--speed slow` | Read twice the field at every search position, for more overlap between them |
| `--check y` | Even out a raw one-shot-colour (Bayer) frame before solving; unbinned raw OSC only |

`arcsec --help` lists everything. ASTAP's single-dash spellings (`-sip`, `-speed slow`,
`-check y`, `-analyse 30`, ...) work as they do there.

**SIP distortion.** With `--sip`, arcsec fits third-order SIP polynomials (`A_p_q`,
`B_p_q` and the inverse `AP_p_q`, `BP_p_q`, `CTYPE` `RA---TAN-SIP`) to the individually
matched stars and writes them to the `.wcs` file and, with `--update`, the FITS header;
the `.ini` stays linear, as ASTAP's does. Unlike `astap_cli`, which adds a cubic
whenever it has 20 stars, arcsec keeps one only if the distortion is statistically
real and the matched stars cover the whole frame. On an undistorted field a cubic fits
only the centroid noise, and on most images of the (distortion-free) benchmark corpus it
made the worst corner worse, for `astap_cli -sip` as for an unconditional fit in arcsec.
On a frame with lens distortion it removes it: a DSS field warped by 8 px of barrel
distortion at the corners goes from 6.1″ worst-corner error to 1.2″.

### Measuring stars without solving

As `astap_cli`, and with the same output, for focusing and quality scripts:

```bash
arcsec -f image.fits --analyse 30     # HFD_MEDIAN=4.5 and STARS=733 on stdout
arcsec -f image.fits --extract 30     # the same, and every star to image.csv
arcsec -f image.fits --extract2 30    # solve, then image.csv with RA and Dec too
```

The value is the minimum SNR of a star (`0` means 30). `--analyse` and `--extract` do no
solve and need no catalogue; they write no `.ini` or `.wcs`. The CSV goes next to the
image whatever `-o` says, as ASTAP's does, with the header
`x,y,hfd,snr,flux,ra[0..360],dec[0..360]`: positions in 1-based FITS pixels, HFD in
pixels, flux in ADU above the background, RA and Dec in degrees (with `--extract`, only
if the header already holds a CD-matrix WCS). `--extract2` writes it whether or not the
solve succeeds, and fits SIP as `--sip` does. `-s` bounds the detection passes as it does
for solving. On Windows `--analyse` reports in its exit code too, as ASTAP does:
`round(HFD × 100) × 1 000 000 + stars`.

### Output files

On a successful solve arcsec writes, next to the image (or at `-o <base>`):

- `<base>.wcs` — the solution as a FITS header (`CRVAL`, `CRPIX`, `CD`, `CDELT`,
  `CROTA`, and SIP terms with `--sip`), as ASTAP and Astrometry.net write it. Always
  written; `--wcs` is accepted for compatibility.
- `<base>.ini` — ASTAP's summary: `PLTSOLVD`, `CRVAL1/2`, `CDELT1/2`, `CROTA2` and the
  fit statistics.

When a solve fails, `<base>.ini` is still written, holding `PLTSOLVD=F` and the command
line, as ASTAP does; tools that poll the `.ini` rely on it.

With `--update` the same keywords are also written into the FITS image's own header,
after removing any `PC` matrix and SIP terms left there by an earlier solve.
The report on stdout follows ASTAP's layout.

### Exit codes

As ASTAP's:

| Code | Meaning |
|---|---|
| `0` | Solved |
| `1` | No solution (also: a command-line usage error) |
| `2` | Not enough stars detected |
| `16` | Image file error (missing, unreadable or unrecognised) |
| `32` | Star database or index files not found |
| `33` | Star database read error |

## Documentation

- [docs/catalogues.md](docs/catalogues.md) — which catalogue for which field, and where
  they install.
- [docs/plate-solving.md](docs/plate-solving.md) — how plate solving works in general,
  and precisely what arcsec does, with flowcharts and the constants table.
- [docs/test-images.md](docs/test-images.md) — the benchmark corpus and measured results.
- [docs/offline-index.md](docs/offline-index.md) — arcsec's blind index: design, file
  format and measured results.
- [CHANGELOG.md](CHANGELOG.md) — what changed in each release.
- [CONTRIBUTING.md](CONTRIBUTING.md) — building, testing, benchmarking and releasing.

The solving library is published separately as
[`arcsec-core`](https://crates.io/crates/arcsec-core), for use from other Rust programs;
see [its README](arcsec-core/README.md). C and C++ programs can use **libarcsec**, a
shared or static library with a C header, released beside the program; see
[libarcsec/README.md](libarcsec/README.md).

## License

MIT — see [LICENSE](LICENSE).
