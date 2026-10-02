//! Image input and output for the [arcsec] plate solver.
//!
//! [`arcsec_core`] solves pixels; this crate gets them out of files. It reads three
//! container formats, telling them apart by their first bytes rather than the
//! extension ([`image_io::detect_format`]):
//!
//! | Format | Reader |
//! |---|---|
//! | FITS (also gzip/bzip2/compress-compressed) | [`fits_io`], via CFITSIO (`rsfitsio`) |
//! | XISF (PixInsight) | [`xisf_io`], via the `xisf` crate |
//! | ASDF (Roman, astropy) | [`asdf_io`], via the `asdf-rs` crate |
//!
//! Each supplies the pixels as an [`arcsec_core::ImageBuffer`] and what the header
//! says about pointing and pixel scale; [`image_io`] dispatches between them and
//! holds the keyword interpretation they share, and [`header`] applies the same
//! rules to FITS header text an application already has in memory.
//!
//! It also writes ASTAP's `.wcs` and `.ini` sidecar files and updates a FITS header
//! in place ([`fits_io`]).
//!
//! The `arcsec` command line and the arcsec C library (`libarcsec`) both read
//! images through this crate. Its API follows their needs and may change between
//! minor versions.
//!
//! [arcsec]: https://github.com/cruzzil/arcsec

// `alloc` is not in the extern prelude for a crate that links std, so it has to
// be declared before alloc:: paths can be written.
extern crate alloc;

pub mod asdf_io;
pub mod fits_io;
pub mod header;
pub mod image_io;
pub mod xisf_io;
