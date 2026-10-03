# Contributing to arcsec

Bug reports, fixes and improvements are welcome. This file covers how to build and
test the project, the rules the code has to keep, and how releases are made.

## Setting up

Use current stable Rust, installed with [rustup](https://rustup.rs/). The repository
does not pin a toolchain. CI runs on the latest stable release, so if CI reports a
clippy lint you do not see locally, run `rustup update` first.

Two other toolchains matter:

- **The MSRV** (minimum supported Rust version) is `rust-version` in the root
  `Cargo.toml`, currently 1.96. The code must build and pass its tests on it; CI checks
  this on every push. It is set by the dependencies (rsfitsio needs 1.96), so it rises
  only when a dependency update requires it, and that is a deliberate change recorded
  in the changelog.
- **Nightly** is optional, for clippy lints that have not reached stable yet:
  `cargo +nightly clippy --workspace --all-targets --all-features -- -D warnings`. CI
  runs this as an advisory job that does not block merging. The code must never
  require nightly.

There are no system dependencies beyond a C compiler, which `ring` (the TLS
cryptography behind the catalogue downloader) uses to build bundled C and assembly.

## Building and testing

These are the commands CI runs; if they pass locally, CI should pass too.

```bash
cargo fmt --all --check
cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --workspace --locked
cargo test --workspace --locked
cargo +1.96 test --workspace --locked      # MSRV; rustup toolchain install 1.96 first
```

Unit tests live next to the code in `#[cfg(test)] mod tests` blocks; there is no
`tests/` directory. To run a subset:

```bash
cargo test -p arcsec-core lsq                              # one crate, name filter
cargo test -- --nocapture spiral_covers_origin_first       # one test, with output
```

CI also runs the CLI against inputs that need no catalogue, and checks its exit codes:
see the "Check the CLI surface" step in `.github/workflows/ci.yml`.

## Benchmarking

A plate solve cannot be tested without star catalogues and real images, neither of
which is in the repository, so the unit tests do not show whether a change still
solves. Before merging anything that touches detection, matching, the search or the
fit, or any hot path, run the benchmark corpus and compare against `main`:

```bash
scripts/fetch-test-images.sh                 # the 103-image corpus, about 2.2 GB
arcsec catalog install d80                   # or point --db at an existing ASTAP directory
cargo build --release
scripts/benchmark.py --db ~/.local/share/arcsec/catalogs --auto-db
scripts/benchmark.py --db ~/.local/share/arcsec/catalogs --auto-db --offset-hint 0.3
```

`catalog install` also builds the blind index into the catalogue directory, and the
solver consults an installed index automatically at `-r` 10° and above. The benchmark's
default `-r 5` never does; to measure the spiral alone at wider radii, install with
`--no-index`, or make sure neither the catalogue directory nor the `--db` directory
holds a `*.arcsecix` (the solver looks in both).

For a change with wider reach — detection, the fit, anything that might behave
differently on real cameras, wide fields or unusual formats — run the expanded corpus too
(635 images from ten archives plus simulated camera artefacts, about 6.5 GB beyond v1;
the fetch is resumable and takes about an hour and a quarter):

```bash
scripts/fetch-corpus.py                      # hard-links the v1 images already fetched
scripts/benchmark.py --corpus --db ~/.local/share/arcsec/catalogs --auto-db
scripts/benchmark.py --corpus --db ~/.local/share/arcsec/catalogs --auto-db --offset-hint 0.3
```

It prints breakdowns by tier, source and field size; `--set v1` restricts it to the
original 103 so numbers stay comparable with earlier results.

The `--offset-hint 0.3` runs start each solve 0.3 fields away from the true centre, as a
mount's reported position would. A perfect hint once hid a rotation error that only showed up
off-centre (docs/test-images.md §6.7).

`benchmark.py` scores every solve against the true WCS at the centre and all four
corners, and counts any solve more than 5″ (or one pixel, if larger) out at a corner as a
false positive; for truths with SIP/TPV distortion it allows for what a linear plate
cannot reach and reports those cases as INEXACT. The
current standing, and the results of changes that did not work and are not worth
repeating, are in [docs/test-images.md](docs/test-images.md). A change that adds a
false positive is a regression however many new solves it brings.

`benchmark.py` can solve every image with other solvers too and print a head-to-head:
`--astap ~/astap_cli` for ASTAP (with `--auto-db` it lets ASTAP pick its own database by
field size; `--astap-db d80` pins one) and `--seiza <binary> --seiza-data <dir>` for
[seiza](https://github.com/theatrus/seiza). Each gets the same position, field size and
radius and is scored by the same rules, false positives included. `--blind-index <index>`
compares blind solving (arcsec with its own index or Astrometry.net index files, seiza
with its own index); `--blind` alone moves the hint to the antipode without adding `-i`. For
timings use `--jobs 1`, alternate `--order` between rounds, and compare like with like:
`--threads 1 --seiza-threads 1 --taskset <cpu>` pins all three to one core. The method and
the latest numbers are in [docs/test-images.md §9](docs/test-images.md#9-results-arcsec-vs-astap-vs-seiza).

`scripts/bench_all.sh` compares speed and positions against ASTAP's `astap_cli` over
`resources/*.fits`. The hot paths depend on presorting, binary search, memory-mapped
catalogues and avoiding per-step allocation, and speed is easy to lose: include
timings in the pull request for any change to one.

Read [docs/plate-solving.md](docs/plate-solving.md) before changing detection,
matching or the WCS fit.

## Rules the code keeps

- **ASTAP compatibility is a contract.** arcsec is a drop-in replacement for
  `astap_cli`, so its command-line flags, stdout format, output files (`.wcs`, `.ini`)
  and exit codes (0, 1, 2, 16, 32, 33) must not change. New options are fine; changing
  what an existing one means is a breaking change. ASTAP flags that are not implemented
  are refused with an error rather than silently ignored.
- **Angles are radians inside `arcsec-core`.** Degrees, hours and arcseconds appear
  only at the CLI boundary and in output files. `--ra` is in hours and `--spd` is
  south-pole distance (90 + Dec) in degrees, both ASTAP conventions.
- **CDELT and CROTA follow astap_cli.** `CDELT1` carries the image's parity (negative
  for the sky's usual handedness, positive when mirrored), `CDELT2` is positive, and
  `CROTA1`/`CROTA2` are the rotations of the +X and +Y axes in the FITS (Calabretta &
  Greisen) sense, exactly as `astap_cli` derives them from the CD matrix
  (`wcs::output::old_style_wcs`). The CD matrix is what readers should use; these
  old-style keywords are there for ASTAP compatibility.
- **`core` and `alloc` before `std`.** Clippy denies `std_instead_of_core` and
  `std_instead_of_alloc` across the workspace, which keeps operating-system
  dependencies visible and out of the pure maths.
- **Threads go through `arcsec_core::max_threads()`**, so `--threads 1` stays genuinely
  single-threaded.
- Image data is row-major `f32`, indexed `data[y * width + x]`.

## Website

The project website lives in `site/`: an
[Astro](https://astro.build/) and [Starlight](https://starlight.astro.build/) project
with its own `package.json`, needing Node.js 22.12 or newer.

```bash
cd site
npm ci
npm run dev        # http://localhost:4321/arcsec/
npm run build && npm run check-links
```

The Website workflow (`.github/workflows/site.yml`) builds every pull request that
touches `site/`, and deploys to GitHub Pages on a push to `main` and after each release.
The site describes the CLI, so a change to options, exit codes, output files or the
catalogue list needs a matching change there: [site/README.md](site/README.md) lists
what to keep in step.

## Pull requests

- Keep each pull request to one change, and explain why as well as what. For a change
  to solving, include benchmark results before and after.
- Add tests for new behaviour where it can be tested without external data.
- Add a line to the `[Unreleased]` section of [CHANGELOG.md](CHANGELOG.md) for anything
  a user would notice: new options, changed behaviour, fixed bugs, a raised MSRV.
  Internal refactoring does not need an entry.
- Commit messages: a short summary line in the imperative ("Fix exit code for a
  missing index"), then a body saying why if it is not obvious.

## Releasing

The `arcsec` and `arcsec-core` crates share one version, set in `[workspace.package]`
in the root `Cargo.toml`.

1. Set the new version in `[workspace.package]` and in the `arcsec-core` entry of
   `[workspace.dependencies]`, then run `cargo check` so `Cargo.lock` follows.
2. In `CHANGELOG.md`, rename `[Unreleased]` to the new version with today's date, add
   a fresh empty `[Unreleased]` above it, and update the comparison links at the bottom.
3. Commit, and let CI pass on `main`.
4. Tag and push: `git tag -a v0.2.0 -m "arcsec 0.2.0" && git push origin v0.2.0`.
   The Release workflow checks that the tag matches the crate version, builds the
   binaries for Linux (x86-64, arm64), macOS (arm64) and Windows (x86-64), and
   publishes a GitHub Release with checksums and the changelog section as its notes.
5. Publish to crates.io, `arcsec-core` first since the CLI depends on it. Either run the
   "Publish to crates.io" workflow from the Actions tab on the tag (dry run first), or
   locally from a clean checkout of the tag:

   ```bash
   cargo publish --workspace --locked --dry-run
   cargo publish --workspace --locked
   ```

   `--workspace` publishes in dependency order and waits for `arcsec-core` to reach the
   index before publishing `arcsec`. A version on crates.io cannot be replaced, only
   yanked, so check the dry run.
