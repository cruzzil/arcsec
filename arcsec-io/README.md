# arcsec-io

Image readers for the [arcsec](https://github.com/cruzzil/arcsec) plate solver: FITS
(via CFITSIO), XISF (PixInsight) and ASDF (Roman, astropy) pixels, plus the pointing
and pixel scale their headers carry, and ASTAP's `.wcs`/`.ini` writers.

It is shared by the `arcsec` command line and the arcsec C library. Most people want
the command-line tool (`cargo install arcsec`); Rust programs that solve in-process
want [`arcsec-core`](https://crates.io/crates/arcsec-core), and this crate only if
they also want arcsec to read the image file. Its API follows those two users and
may change between minor versions.

Licensed under the MIT licence.
