# Changelog

All notable changes to arcsec are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html). The `arcsec` and
`arcsec-core` crates share one version number.

The ASTAP-compatible command line - flags, stdout format, output files and exit codes -
is part of the public interface: a change to any of them is listed here, and a breaking
one is called out as such.

## [Unreleased]

### Added

- **Fuzz targets** for every file arcsec reads (FITS, XISF and ASDF images, the ASTAP
  star databases, Astrometry.net and arcsec index files, catalogue archives) and for the
  numbers that reach the solver, in a `fuzz/` crate run with cargo-fuzz, with a weekly
  advisory CI job. See CONTRIBUTING.md, "Fuzzing". The fixes below came from it.

### Changed

- `-t` (quad tolerance) above 0.1 is refused with exit 1. The ratios it compares lie in
  0..1, so a large tolerance matched every pattern to every other, never solved, and with
  `--method tetra` could ask for gigabytes.

### Fixed

- The Astrometry.net index loader read each table at a fixed row width, so a table
  declaring narrower rows could make it allocate up to 16 times the file's size before
  the read failed. Row counts are now checked against the file (with overflow checks)
  before allocating, and a negative row count is refused.
- **A malformed image is an unreadable file (exit 16), never a crash.** Header
  dimensions that the file is too short to hold, or beyond 2³⁰ pixels, are refused before
  anything is allocated for them; a 5 KB FITS file could ask for 4 GB, and an XISF or
  ASDF geometry could overflow. Panics inside the FITS reader (rsfitsio 0.470.3 panics on a
  header value that is not text and on several kinds of corrupt tile-compressed image) and
  the XISF and ASDF readers are caught at the reader and reported, with where they
  happened; a keyword that cannot be read counts as missing rather than failing an
  otherwise good image.
- A corrupt header can no longer stall or exhaust the solver. A pixel scale from
  FOCALLEN/XPIXSZ outside 0.001–10 000″/px, or an RA/Dec that is not a number, is ignored
  as missing. The search spiral is no longer built in memory before it is walked: a
  narrow field with the default `-r 180` allocated hundreds of megabytes up front, and a
  vanishing field size aborted; a radius of more than a million fields is refused.
- A very large `-s` no longer reserves memory in proportion to it.
- `arcsec catalog install`: a corrupt xz payload in a `.deb` is an error rather than a
  panic, and each file is written beside its name and renamed into place. Rewriting a
  star-database tile in place could crash a solve that had it memory-mapped (SIGBUS),
  an interrupted install left a short tile that then read as corrupt, and a symbolic
  link at a tile's name was written through.
- Any panic that remains is reported as an internal error with exit 1 and a
  `PLTSOLVD=F` `.ini`, as for any failed solve, rather than Rust's exit 101 and no `.ini`.

## [0.4.0] - 2026-10-02

### Added

- **Catalogue-seeded fallback search** for fields whose brightness ranking disagrees
  with the catalogue's (crowded and nebulous fields, saturated bright stars, cluster
  cores, infrared passbands). When the spiral finds nothing, quads are built from the
  catalogue's brightest stars about the hint and looked for among the brightest 2000
  detections without ranking them — a length-sorted table of image star pairs and a
  position hash, as in [seiza](https://github.com/theatrus/seiza)'s rank-robust
  fallback — and any plate found is verified exactly as a spiral position's. It runs
  once, after the spiral, under a fixed work budget (about a third of a second of one
  core), so it changes nothing the spiral solves.

### Changed

- **Failed searches are several times faster, with identical results.** A search that
  finds nothing visits every spiral position out to `-r`, and each position rebuilt and
  matched tens of thousands of catalogue quads; that work is now 4–8× cheaper on one
  thread. The image quads are bucketed once per solve on two of their ratios, so a
  catalogue quad is compared only with the few image quads it could match (80 % of a
  failed search went on the comparisons); the neighbour search, duplicate check and
  distance sort that build each position's quads no longer scan every star, chase
  pointers or branch at random, and the catalogue quads are no longer sorted at every
  position. Every result on the 635-image corpus, true-centre and offset hint, with and
  without an index installed, is the same as before, down to the bytes of the `.wcs`
  files. On one thread the 98-image benchmark subset takes 115 s instead of 581 s, and
  the slowest failure in it 45 s instead of over 300 s.
- Spiral positions are handed to the worker threads one at a time instead of in
  batches, so no thread waits for the slowest position of a batch. The answer is the
  serial search's, as before. The line of step distances printed without `--progress`
  now lists the positions up to the solution whatever the thread count (one thread
  printed exactly that before; more threads printed to the end of the solution's batch).
- With a blind index installed and `-r` of 10° or more, a hint far from the field no
  longer costs the whole spiral: when the index finds nothing within `-r`, it is asked
  once more without the limit, and if it verifies the field more than two fields beyond
  `-r`, the search stops there. The result is unchanged, "No solution found." and exit
  code 1, as `-r` requires; stderr says how far away the index placed the field.
- **Bright stars too large for the measuring box are measured.** Detection refused any
  star whose 3σ isophote reached 14 pixels from its seed — the saturated discs of a
  photographic plate, and the wide wings of an undersampled TESS star — and those are
  the catalogue's brightest stars. They are now measured again in a 32-pixel box, at
  5% of their peak, and refused only if still too large, not a disc, or elongated.
  They do not count towards the `-s` stars that end the detection cascade.
- **Verification needs significance as well as a count.** A plate must verify at least
  four times as many stars as chance would pair at the frame's star density, with an
  rms within the last match radius. In a dense frame (a 2.2° TESS crop of 500 stars on
  384 × 384 pixels) a wrong plate can otherwise reach the 30-star minimum by chance;
  every correct solve on the benchmark corpus verifies at least 7.9 times the chance
  count.
- When the distortion model's full-frame match pairs at least 1.5 times as many stars
  as the verified plate, the linear plate is refitted to those pairs: with the hint a
  third of a field off, the verified plate fitted only the part of the frame its
  search position's catalogue covered.
- Benchmark corpus (635 images), against 0.3.0: with the true-centre hint 558 → 589
  correct, false positives 1 → 1, no image lost; with the hint 0.3 fields off 537 → 580,
  1 → 1. Tier-D controls still all refused. Of the 22 fields seiza solved and arcsec did
  not, 17 now solve. On one thread the median solved image takes 2% longer (wide
  fields full of saturated stars up to twice as long), a failed full spiral 0.5%
  longer, and a quick failure 0.1–0.3 s longer for the fallback.

## [0.3.0] - 2026-10-02

### Added

- Using arcsec with Siril: a README section and a website page,
  [Use with Siril](https://cruzzil.github.io/arcsec/siril/). Siril 1.4 has no setting
  for an external ASTAP solver, so the page covers solving with `--update` and letting
  Siril read the solution from the header, for single images and for sequences; tested
  with Siril 1.4.4 on Linux.
- Distortion handling in the catalogue solve. Once a position verifies, the solver
  re-reads the catalogue about the image centre and fits a polynomial plate (up to
  cubic, order chosen by F-tests and frame coverage) by re-matching every catalogue
  star as the model improves, astrometry.net `tweak` style. When the field is
  measurably distorted (F ≥ 30 over a linear fit, a pixel or more of difference, the
  pairs covering the frame) the reported plate is the linear plate closest to the model
  over the whole frame instead of the one the centre's stars give, and the model's star
  pairs, which reach the corners, become `WcsSolution::matched_stars`. The written WCS
  stays linear and ASTAP-compatible; `--sip` now fits its SIP terms to those pairs. On
  the 635-image corpus: TESS 9 → 23 of 42 correct (the 12° frames' corners from
  2000–2800″ off to the ~1000″ linear floor; with `--sip` 15–41″), WISE with the offset
  hint 1 → 20 of 25, the synthetic lens set 3 → 10 of 10. Undistorted fields are
  unchanged.
- A position whose quads agree strongly (≥ 50 pairs) but whose linear plate fails
  verification is retried with the distortion model before the search moves on.
- **Blind index built from the installed star database.** `arcsec catalog index build`
  writes `<db>.arcsecix` (fields 0.3°–30° by default, `--min-fov`/`--max-fov` to
  choose; about 2 minutes and 290 MB from D80, nothing downloaded), and
  `arcsec catalog index info` describes it; `catalog list` and `catalog verify` include
  it. `-i` accepts it (a file, or a directory holding one) and solves blind: on the
  benchmark corpus it finds 473 of the 558 fields the hinted solver finds with the true
  centre (399 with the default 0.3° index, which drops the narrowest tier), typically
  in half a second; on fields of 0.6° and wider the Astrometry.net 4100 series finds
  70 of 254 to the index's 216, at 26 s against 0.8 s. Every position it reports passes the ordinary solver's star-level
  verification. Without `--fov` or FOCALLEN/XPIXSZ it searches pixel scales of
  0.3–60″/px instead of assuming 1″/px.
- With an index installed, a search of `-r` 10° or more reaching more than five fields
  from the hint (such as N.I.N.A.'s blind mode, `-r 180`) tries the index once the first
  five fields have failed, and falls back to the full search if it finds nothing; below
  10°, and within five fields, the result is unchanged. No new flags; the ASTAP-compatible command line is unchanged.
- Library: `arcsec_core::index` (format, builder, `BlindIndex`) and
  `pipeline::index_solve`; `catalog::for_each_star_in_dec_band`.
- Benchmark tooling: `scripts/benchmark.py --blind` (hint at the antipode, so only a
  blind index can find the field), `--blind-index`, `--no-fov`, and
  [seiza](https://github.com/theatrus/seiza) as a third solver (`--seiza`) alongside
  ASTAP, with timing options (`--order`, `--taskset`). The comparison is in
  `docs/test-images.md` §9.

### Changed

- A strongly distorted field that the model cannot follow over the whole frame (a
  significant cubic 3 px or more from the verified plate where there are stars, but
  too little of the frame covered to fit it there) is now refused (exit 1, no
  solution) rather than reported with a linear plate fitted to part of it. No corpus
  image is affected.

- The catalogue solve uses at most as many image stars as the database can hold in the
  field, its density times the field's area (ASTAP's "database limit"): the brightest
  `min(-s, density × area)` detections. Small, crowded fields with d80 (below ~0.25°)
  no longer build their quads from stars the catalogue does not have. Library:
  `catalog::database_density`.
- Sparse images solve. When the image yields fewer stars than it may use, the catalogue
  read is denser than the image; when it is at least 2.5 times denser, the catalogue
  spiral now also builds quads from the catalogue's brightest stars at the image's
  density and adds them to the full-depth ones. And an image with fewer than 194 detections needs fewer than 30 matched stars,
  15 % of its detections but at least 10; below 30 the solution must also have the
  pixel scale the hint implies (within 10 %) and a star-level rms of at most 0.5 px.
  Narrow SkyMapper frames and LCO frames with few stars gain most.
  `SolveParams::fov` is documented as the long side, which is what the CLI passes; the
  scale check relies on it.

- The Astrometry.net blind path ranks its vote cells by (RA, Dec, ln scale) with the RA
  bin widened by 1/cos δ, smooths each over its neighbours and verifies the medoid of
  the strongest bucket rather than the first hypothesis of each cell.

### Fixed

- The plate fit's similarity check compared the lengths of the matrix *rows*, which a
  sheared plate passes: a tier-D control solved to a plate stretching the image three
  times more one way than the other. `solve_plate_constants` now requires the ratio of
  the plate's two singular values to be at most 1.08 (`math::lsq::plate_anisotropy`,
  `MAX_PLATE_ANISOTROPY`); the largest on any correct corpus solve is 1.027.
  `ArcsecError::BadSolution::ratio` now carries that singular-value ratio rather than
  the squared row-norm ratio, and its message changes accordingly.
- A few wrong quads in the winning vote could drag the plate fit off a similarity, and
  the search abandoned the right position. The quad path now sigma-clips the matched
  quad centroids before fitting, as the triangle path already did. Nebulous and crowded fields (Coalsack, B68, M16), coarse DSS and SHASSA fields
  and TESS frames gain most.
- `--update` now removes the `PC` matrix and SIP terms of an earlier solution before
  writing its own. Left in place, a `PC` matrix takes precedence over the new `CD` matrix
  in wcslib, astropy and most other readers, and combined with arcsec's `CDELT` it
  described a mirrored field: re-solving an image Siril had already solved (Siril writes
  `PC` + `CDELT` with SIP) put the corners of a 1° frame about a degree out for every
  reader except Siril. `astap_cli -update` leaves these keywords behind too.

## [0.2.0] - 2026-10-01

### Added

Every ASTAP command-line option is now implemented; none is refused any more. Output
follows `astap_cli`'s, checked against it.

- `--analyse <snr_min>`: measure without solving. Prints `HFD_MEDIAN=` and `STARS=`,
  writes no `.ini` or `.wcs`, needs no catalogue. On Windows the exit code also carries
  the result, `round(HFD × 100) × 1 000 000 + stars`, as ASTAP's does.
- `--extract <snr_min>`: as `--analyse`, and writes every star to `<image>.csv`
  (`x,y,hfd,snr,flux,ra[0..360],dec[0..360]`; RA/Dec when the header already holds a WCS).
  As in ASTAP the CSV goes next to the image whatever `-o` says.
- `--extract2 <snr_min>`: solve, then write the same CSV with every star's RA and Dec,
  whether or not the solve succeeded.
- `--sip`: third-order SIP distortion terms (`A_p_q`, `B_p_q`, `AP_p_q`, `BP_p_q`,
  `CTYPE RA---TAN-SIP`) in the `.wcs` file and with `--update`, in ASTAP's layout. Unlike
  ASTAP they are added only when the distortion is statistically significant and the
  matched stars cover the frame, since on an undistorted field a cubic only adds noise
  at the corners; `--sip n` turns it off, as in ASTAP. `--extract2` implies it.
- `--speed slow`: read a catalogue window twice the field at every search position.
- `--check y` (or a bare `--check`): even out the Bayer pattern of a raw one-shot-colour
  frame before solving.
- Library: `detection::analyse_image`, `wcs::sip` (`fit_sip`, `Sip`, `TanWcs`),
  `ImageBuffer::check_pattern_filter`, `pipeline::SearchSpeed`, and the verified star
  pairs of a solve in `WcsSolution::matched_stars`.
- Benchmark tooling only (no change to the solver or its command line): an expanded
  635-image benchmark corpus (`scripts/corpus.tsv`, fetched by `scripts/fetch-corpus.py`)
  drawn from ten public archives plus simulated camera artefacts and more negative
  controls, and `scripts/benchmark.py --corpus` with per-tier, per-source and per-FOV
  breakdowns and SIP/TPV/SIN-aware truth. The original 103 images remain the `v1`
  subset. See `docs/test-images.md`.
- A project website, [cruzzil.github.io/arcsec](https://cruzzil.github.io/arcsec/):
  downloads, setting up N.I.N.A., a "which catalogue do I need?" picker, the catalogue
  guide, the command-line reference and an FAQ. It is now the crates' homepage.

### Changed

- **Breaking for `arcsec-core` users:** `SolveParams` has a new `speed` field and `WcsSolution` new
  `matched_stars` and `sip` fields, so code that builds them with struct literals needs
  to set them (`SearchSpeed::Auto`, `Vec::new()`, `None`).
- `--sip`, `--check` and `--speed` take ASTAP's optional values (`-sip n`, `-check y`,
  `-speed slow`); anything else is a usage error.

### Fixed

- The catalogue read now returns the field's brightest stars from every database tile
  the field overlaps. The `.1476` reader (D-series) filled its whole star budget from the
  first tile, so a field straddling a tile boundary had catalogue stars on one side only
  and its fit rested on half the frame; the `.290` reader (G05) gave each tile a fixed
  share, under-sampling a tile that covered most of the field. On the 635-image
  benchmark corpus this solves 29 more images with the true centre as hint (475 → 504
  of 596) and 30 more with the hint 0.3 fields off, with fewer false positives in both,
  and it is no slower. It changes which catalogue stars a multi-tile field is matched
  against, so a few marginal solves change either way (docs/test-images.md §7.7).
- ASTAP compatibility: the `Start position:` and `Solution found:` lines on stdout now
  match `astap_cli` byte for byte. Every field is two digits wide and the start position
  has ASTAP's comma: `Start position: 04: 20  00.0, +35d 00  00` where arcsec printed
  `4: 20  0.0 +35d 00  0`. N.I.N.A. reads only the `.ini` file, so it never saw the
  difference, but a script that parses stdout could. One deliberate difference remains:
  an RA that rounds up to 24h prints as `00:` where ASTAP prints `24:`.

## [0.1.2] - 2026-09-28

### Fixed

- gzip-compressed FITS (`.fits.gz`) is read directly instead of being refused. This
  needs rsfitsio 0.470.3, which fixed its compression magic numbers; the same fix covers
  Unix `compress` (`.Z`) files, which 0.1.0 listed as supported but could not open.
- A corrupt or empty Astrometry.net index file in the `-i` directory is skipped instead
  of crashing the blind solve (rsfitsio 0.470.3 returns an error where it panicked).

## [0.1.1] - 2026-09-28

### Fixed

- **Rotation error when the hint is off-centre.** The plate was fitted in the tangent
  plane of the search position that matched rather than at the image centre, so the
  reported rotation (and so the corners) drifted with the distance from the hint and
  with declination: 0.1° at Dec −20 to 1.7° at Dec −80 with the hint 0.3 fields off,
  and up to 1600″ at the corners of a 10° field. The centre position was unaffected.
  Solves are now refitted at the image centre. With the benchmark's hint 0.3 fields off,
  false positives fall from 47 to 0 and correct tier-A solves rise from 8 to 55.
- `--method tetra`: outlier clipping gave up when gross outliers skewed its first fit,
  and so abandoned positions it could solve. Tetra now solves 29 of the 98 solvable
  benchmark images, up from 20.
- Blind solving builds its patterns from the brightest detected stars on sparse fields
  too; they were taken in scan order when fewer than `-s` stars were detected.

### Added

- Much wider test coverage of `arcsec-core` (library lines 55% to 96%), including
  end-to-end solves against synthetic catalogues in all three database formats and
  synthetic Astrometry.net indexes.

## [0.1.0] - 2026-09-27

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
- **N.I.N.A. support**: arcsec can be set as N.I.N.A.'s "ASTAP" solver. `--fov` is the
  image height, as N.I.N.A. sends it; the `.ini` carries the keys N.I.N.A. reads
  (`CRPIX1/2` and the `CD` matrix, as well as ASTAP's others); and the Windows
  executable has a version resource, without which N.I.N.A. refuses automatic
  downsampling.
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

[Unreleased]: https://github.com/cruzzil/arcsec/compare/v0.4.0...HEAD
[0.4.0]: https://github.com/cruzzil/arcsec/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/cruzzil/arcsec/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/cruzzil/arcsec/compare/v0.1.2...v0.2.0
[0.1.2]: https://github.com/cruzzil/arcsec/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/cruzzil/arcsec/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/cruzzil/arcsec/releases/tag/v0.1.0
