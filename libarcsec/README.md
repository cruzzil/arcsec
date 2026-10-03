# libarcsec — the arcsec plate solver as a C library

libarcsec lets a C or C++ program solve images in-process with
[arcsec](https://github.com/cruzzil/arcsec): hand it pixels (or a file) and a rough
idea of the field, get back a FITS WCS. It makes the same decisions the `arcsec`
command line does — which star database suits the field, how far to bin, when a blind
index helps — because it calls the same code (`arcsec_core::auto`), and it reads the
same catalogues from the same place (`arcsec catalog install`).

The API is declared in [`include/arcsec.h`](include/arcsec.h); every function and field
is documented there. This page is the overview: building, linking, the rules for
memory, threads and errors, and how the ABI is versioned.

## Contents

- [Building](#building)
- [Linking](#linking)
- [A first program](#a-first-program)
- [The API at a glance](#the-api-at-a-glance)
- [Memory ownership](#memory-ownership)
- [Threads](#threads)
- [Errors](#errors)
- [Versioning](#versioning)
- [Platform notes](#platform-notes)
- [How the library is made safe](#how-the-library-is-made-safe)
- [Testing](#testing)

## Building

You need Rust (1.96 or newer, from [rustup](https://rustup.rs/)) and a C compiler.

```bash
git clone https://github.com/cruzzil/arcsec && cd arcsec
libarcsec/dist.sh /tmp/arcsec        # build, and lay out an install tree there
sudo cp -a /tmp/arcsec/. /usr/local/ # optional: install it
sudo ldconfig                        # Linux, after installing
```

`dist.sh` builds the `dist` profile (optimised, no debug info) and produces:

| Path | Linux | macOS | Windows (MSVC) |
|---|---|---|---|
| `include/` | `arcsec.h` | `arcsec.h` | `arcsec.h` |
| `lib/` (shared) | `libarcsec.so.1`, `libarcsec.so` → it | `libarcsec.dylib` | `arcsec.lib` (import library) |
| `bin/` | | | `arcsec.dll` |
| `lib/` (static) | `libarcsec.a` | `libarcsec.a` | `arcsec_static.lib` |
| `lib/pkgconfig/` | `arcsec.pc` | `arcsec.pc` | `arcsec.pc` |
| `lib/cmake/arcsec/` | `arcsecConfig.cmake`, `arcsecConfigVersion.cmake` | same | same |

The pkg-config and CMake files locate everything relative to themselves, so the tree
works wherever it is copied. `dist.sh --target <triple>` cross-builds (the Rust target
must be installed); `--profile release` keeps debug info.

Release archives (`arcsec-lib-vX.Y.Z-<target>.tar.gz`/`.zip` on the GitHub release)
hold the same tree, built for Linux x86-64 and arm64 (glibc 2.35 or newer), macOS
arm64 and Windows x86-64.

With plain cargo: `cargo build --release -p libarcsec` leaves `libarcsec.so`,
`libarcsec.dylib` or `arcsec.dll`, and the static library, in `target/release/`.

## Linking

pkg-config (Meson's `dependency('arcsec')`, autotools, Make):

```bash
cc myapp.c $(pkg-config --cflags --libs arcsec)
cc myapp.c $(pkg-config --cflags arcsec) /path/to/libarcsec.a $(pkg-config --static --libs-only-l arcsec)
```

CMake:

```cmake
find_package(arcsec 0.4 REQUIRED)            # CMAKE_PREFIX_PATH=<install tree>
target_link_libraries(myapp PRIVATE arcsec::arcsec)          # shared
target_link_libraries(myapp PRIVATE arcsec::arcsec_static)   # static
```

Use the **shared library** if your program also links CFITSIO or bzip2. arcsec reads
FITS with a Rust port of CFITSIO that defines CFITSIO's C functions under their real
names (`ffopen`, `fits_open_file`, ...), as its bzip2 does (`BZ2_*`). The shared
library exports only the `arcsec_*` API, so the two never meet; the static library
cannot hide them, and linking it next to the real CFITSIO collides.

## A first program

```c
#include <arcsec.h>
#include <stdio.h>

int solve(const float *pixels, uint32_t width, uint32_t height) {
    if (arcsec_abi_version() != ARCSEC_ABI_VERSION) return -1; /* wrong library */

    arcsec_image image = {0};
    image.struct_size = sizeof image;
    image.data = pixels;              /* row 0 = FITS row 1 */
    image.pixel_type = ARCSEC_PIXEL_F32;
    image.width = width;
    image.height = height;

    arcsec_solve_options opts;
    arcsec_solve_options_init(&opts); /* sets struct_size and the CLI's defaults */
    opts.has_hint = 1;
    opts.ra_deg = 83.82;
    opts.dec_deg = -5.39;
    opts.search_radius_deg = 10.0;
    opts.pixel_scale_arcsec = 1.2;    /* or fov_deg = image height in degrees */

    arcsec_solver *solver = arcsec_solver_new();
    arcsec_result *result = NULL;
    arcsec_status st = arcsec_solve(solver, &image, &opts, &result);
    if (st == ARCSEC_OK) {
        arcsec_wcs wcs = {0};
        wcs.struct_size = sizeof wcs;
        arcsec_result_wcs(result, &wcs);
        printf("CRVAL %.5f %.5f, %.3f\"/px, %u stars\n", wcs.crval1, wcs.crval2,
               wcs.pixel_scale_arcsec, wcs.matched_stars);
    } else {
        fprintf(stderr, "%s: %s\n", arcsec_status_string(st), arcsec_last_error());
    }
    arcsec_result_free(result);       /* NULL is fine */
    arcsec_solver_free(solver);
    return st;
}
```

[`tests/c/smoke.c`](tests/c/smoke.c) exercises every function and is a fuller example.

## The API at a glance

| Area | Functions and types |
|---|---|
| Version | `arcsec_version()`, `arcsec_abi_version()`, `ARCSEC_ABI_VERSION` |
| Solving | `arcsec_solver_new/free`, `arcsec_solve` (pixels), `arcsec_solve_file` (FITS, XISF, ASDF), `arcsec_solver_cancel` |
| Input | `arcsec_image` (u8/u16/i16/u32/i32/float/double; planes or interleaved-by-plane; strides; `ARCSEC_IMAGE_TOP_DOWN`), `arcsec_solve_options` + `arcsec_solve_options_init` |
| Result | `arcsec_result_wcs` (CRVAL, CRPIX, CD; CDELT1/2 and CROTA1/2 as `astap_cli` writes them; SIP A/B/AP/BP up to order 9, RMS, matched stars), `arcsec_result_info`, `arcsec_result_matched_stars`, `arcsec_result_fits_header` (80-character cards), `arcsec_result_pixel_to_sky` / `_sky_to_pixel`, `arcsec_result_database`, `arcsec_result_free` |
| Analysis | `arcsec_analyse` (star count, median HFD, the stars), like `arcsec --analyse` |
| Catalogues | `arcsec_default_catalog_dir`, `arcsec_default_database_dir`, `arcsec_has_star_database`, `arcsec_select_database`, `arcsec_has_blind_index` |
| Diagnostics | `arcsec_last_error`, `arcsec_status_string`, `arcsec_set_log_callback`; per solve, `progress` and `cancel` callbacks in the options |

What the options mean, in brief (the header has the detail):

- **Where.** `has_hint`, `ra_deg`, `dec_deg` (degrees) and `search_radius_deg`
  (default 180, the whole sky). Without a hint, the search starts at RA 0, Dec 0 — so
  give one, or an index.
- **Scale.** `pixel_scale_arcsec` or `fov_deg` (the image *height*, as ASTAP's `-fov`).
  Without either, `fits_header` keywords FOCALLEN/XPIXSZ/XBINNING are used, and failing
  those 1″/px is assumed (with a warning in the log). Pass the scale.
- **Catalogues.** `catalog_dir` (NULL: where `arcsec catalog install` puts them),
  `database` (NULL: chosen by field size, D80 … W08), `index_path` (a blind index: tried
  first without a hint, as a fallback with one unless `index_first`), `auto_index`
  (consult an installed arcsec index for searches ≥ 10° wide; default on).
- **Solver.** `max_stars` (500), `downsample` (0 = automatic), `threads` (0 = every
  core, per call), `sip_order` (≥ 2 fits SIP, currently order 3), `slow`,
  `check_pattern` (Bayer mosaics), `quad_tolerance`, `hfd_min_arcsec`, `method`.
- **Callbacks.** `cancel(user)` returns nonzero to stop; `progress(user, fraction,
  stage)` reports "detecting stars", "searching" (fraction of the search started) and
  "blind index".

Coordinates: pixel coordinates in and out are **1-based FITS pixels of the image as
passed**, unbinned, with row 1 the first row of the buffer — FITS order, as CFITSIO
reads an image and Siril keeps it. A buffer whose first row is the *top* of the picture
sets `ARCSEC_IMAGE_TOP_DOWN`. Angles are degrees; RA in [0, 360).

## Memory ownership

- **Nothing is freed across the boundary except through the library's own `*_free`
  functions.** `arcsec_solver_new` and the solve functions allocate an
  `arcsec_solver` / `arcsec_result`; free them with `arcsec_solver_free` /
  `arcsec_result_free` (both accept NULL). Never `free()` them.
- **Everything else is caller-owned.** Images, options and strings you pass are read
  during the call and not kept. Structs and arrays the library fills (`arcsec_wcs`,
  `arcsec_solve_info`, `arcsec_analysis`, star arrays, string buffers) are yours; you
  allocate them and say how big they are.
- **Library strings are borrowed.** `arcsec_version()` and `arcsec_status_string()`
  return static strings; `arcsec_result_database()` lives as long as its result;
  `arcsec_last_error()` until the next failing call on the same thread.
- Strings into buffers follow `snprintf`: the return value is the full length, so call
  with `len` 0 to size the buffer, and a return `>= len` means it was cut.

## Threads

- **Separate solver handles may solve at the same time on separate threads.** Each
  call has its own options, thread budget (`threads`, which applies to that call only,
  not the process) and cancellation. Siril-style parallel sequences — several threads,
  each with its own handle and `threads = 1` — are the intended use.
- **One handle runs one solve at a time.** A second `arcsec_solve` on a busy handle
  returns `ARCSEC_BUSY` immediately rather than waiting.
- **`arcsec_solver_cancel` may be called from any thread**, at any time, including
  while the solve is running. It affects only the solve in progress. Don't free a
  handle while a solve on it is running: cancel, wait for the call to return, then free.
- **Results are immutable**: any number of threads may read one at once.
- **Callbacks run on arcsec's worker threads**, possibly several at once: the cancel,
  progress and log callbacks must be thread-safe, quick, and must not call back into
  arcsec. They must not throw C++ exceptions or `longjmp` out (Rust aborts the process
  rather than unwinding through them).
- **The log callback is process-wide** (`arcsec_set_log_callback`): one for all
  solvers, so messages from concurrent solves interleave. INFO is a handful of lines per
  solve; the per-position search messages arrive at DEBUG.
- `arcsec_last_error()` is per thread: read it on the thread that made the failing call.

## Errors

Every call that can fail returns an `arcsec_status` (or NULL, 0 or -1, as its
documentation says) and leaves a message for `arcsec_last_error()`. The solve codes
have the values of arcsec's (and ASTAP's) exit codes:

| Status | Value | Meaning |
|---|---|---|
| `ARCSEC_OK` | 0 | Solved (or the call succeeded) |
| `ARCSEC_NO_SOLUTION` | 1 | The search ended without a verified solution |
| `ARCSEC_INSUFFICIENT_STARS` | 2 | Too few stars detected |
| `ARCSEC_FILE_ERROR` | 16 | The image file could not be read |
| `ARCSEC_DATABASE_NOT_FOUND` | 32 | No star database (or index) where you said |
| `ARCSEC_DATABASE_ERROR` | 33 | A catalogue or index file is unreadable |
| `ARCSEC_INVALID_ARGUMENT` | 100 | NULL, out-of-range, or `struct_size` not set |
| `ARCSEC_CANCELLED` | 101 | Cancelled through the handle or the callback |
| `ARCSEC_BUSY` | 102 | The handle is already solving on another thread |
| `ARCSEC_INTERNAL_ERROR` | 199 | A bug in arcsec, caught at the boundary; please report it |

New codes may be added; treat unknown values as failures.

## Versioning

- **`ARCSEC_ABI_VERSION`** (1) changes only for a change that would break a program
  compiled against an older header: a function removed or re-typed, a field moved. It
  is the number in the Linux soname (`libarcsec.so.1`), so an incompatible library is
  never loaded in place of a compatible one. Check `arcsec_abi_version() ==
  ARCSEC_ABI_VERSION` at start-up if you load the library dynamically.
- **Within an ABI version, additions are compatible.** New functions, new status
  codes and enum values, and new fields *at the end* of a struct. Every struct you pass
  in or get filled starts with `struct_size`: the library reads and writes only as much
  as both sides know, and gives new fields their defaults when an older program does
  not have them. Always set it to `sizeof` (or use `arcsec_solve_options_init`).
- **`arcsec_version()`** is the arcsec release (`"0.4.0"`) the library was built from;
  the C library is released with the CLI, from the same tag.

libarcsec is new: until arcsec 1.0, an ABI change is possible in a minor release. It
will bump `ARCSEC_ABI_VERSION` and the soname, and be listed in the changelog.

## Platform notes

- **Linux.** The soname is `libarcsec.so.1`, and the library exports only `arcsec_*`
  symbols. Install with `ldconfig`, or set `LD_LIBRARY_PATH`/an rpath for a private copy.
- **macOS.** The install name is `@rpath/libarcsec.dylib`; `arcsec.pc` adds
  `-Wl,-rpath,<libdir>`. In an app bundle, copy the dylib into `Frameworks/` and add
  `@executable_path/../Frameworks` to the rpath. The dylib still exports the Rust
  CFITSIO's symbols (the Linux build hides them with a linker option that has no
  drop-in macOS equivalent alongside rustc's own export list): with two-level
  namespaces nothing is interposed at run time, but put the real CFITSIO *before*
  libarcsec on the link line so your own calls bind to it.
- **Windows.** `arcsec.dll` is built with MSVC and a static C runtime, so it needs no
  redistributable and works from MinGW-built programs too (Siril's MSYS2 build): link
  `arcsec.lib`, or the DLL directly with MinGW. As on macOS, list CFITSIO before
  arcsec on the link line. Paths are UTF-8.

## How the library is made safe

All of arcsec's `unsafe` FFI code is in this crate; `arcsec-core` and `arcsec-io` hold
the solver and stay free of it. Every entry point:

- runs under `catch_unwind`: a Rust panic becomes `ARCSEC_INTERNAL_ERROR` with the
  panic's message, never an unwind into C (panic reports go to the log callback when
  one is set, not to stderr);
- checks pointers for NULL and sizes for overflow before use, bounds image dimensions
  and strides, refuses an image larger than memory allows rather than aborting, and
  rejects non-finite or out-of-range numbers;
- reads caller structs through their `struct_size` with unaligned copies, so a struct
  from an older or newer header is read correctly and never overrun;
- hands out only `Box`ed handles that its own `*_free` functions release.

Callbacks are `extern "C"` function pointers; the `user` pointers are passed back and
never dereferenced.

## Testing

```bash
cargo test -p libarcsec                 # the FFI tests, the header check, and the C test
libarcsec/tests/c/run-checks.sh           # shared, static, ASan+UBSan+LSan, valgrind (Linux)
ARCSEC_BLESS=1 cargo test -p libarcsec header_is_current   # regenerate include/arcsec.h
```

The header is generated by [cbindgen](https://github.com/mozilla/cbindgen) from the
Rust source (`cbindgen.toml`) and committed; a test fails if it is out of date. The C
test (`tests/c_api.rs`) builds a synthetic star field and its database, compiles
`tests/c/smoke.c` with the platform's C compiler against the library, and runs it.
