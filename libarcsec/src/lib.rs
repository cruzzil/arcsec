//! libarcsec: the arcsec plate solver as a C library.
//!
//! Everything here is the boundary between C and Rust, and nothing else: the
//! solving is `arcsec_core`'s, the decisions the command line makes for a user
//! (database by field size, binning, blind index use) are
//! `arcsec_core::auto`'s, and image files are read by `arcsec_io`. This crate
//! is where all of arcsec's `unsafe` FFI code lives, so that those crates stay
//! free of it.
//!
//! The C API is declared in `include/arcsec.h`, which is generated from this
//! source by cbindgen and checked by a test (`ARCSEC_BLESS=1 cargo test -p
//! libarcsec header` regenerates it). `README.md` describes the API for C
//! programmers: ownership, threading, errors and versioning.
//!
//! # Rules every entry point follows
//!
//! - It never unwinds into C: its body runs under [`std::panic::catch_unwind`]
//!   ([`error::guard`]), and a panic becomes `ARCSEC_INTERNAL_ERROR`.
//! - It checks every pointer for NULL and every size for overflow before use, and
//!   reads caller structs through their `struct_size` ([`util::read_versioned`]),
//!   so a struct from an older or newer header is handled, not misread.
//! - Memory crosses the boundary only as opaque handles freed by their own
//!   `*_free` function, or as caller-owned buffers that this side fills.
//! - Callbacks into C are `extern "C"`: a C++ exception or `longjmp` out of one
//!   is not supported (Rust aborts rather than unwinding through such a frame).

// C names: arcsec_image, arcsec_status, ARCSEC_OK.
#![allow(non_camel_case_types)]
// Every extern "C" fn takes raw pointers by design; each documents its contract
// under "# Safety", which clippy::missing_safety_doc (on by default) enforces.
#![allow(clippy::not_unsafe_ptr_arg_deref)]

extern crate alloc;

mod analyse;
mod catalog;
mod error;
mod image;
mod logging;
mod options;
mod result;
mod solver;
mod util;

#[cfg(test)]
mod tests;

use core::ffi::c_char;

pub use analyse::{arcsec_analyse, arcsec_analysis, arcsec_star};
pub use catalog::{
    arcsec_default_catalog_dir, arcsec_default_database_dir, arcsec_has_blind_index,
    arcsec_has_star_database, arcsec_select_database,
};
pub use error::{arcsec_last_error, arcsec_status, arcsec_status_string};
pub use image::{ARCSEC_IMAGE_TOP_DOWN, arcsec_image, arcsec_pixel_type};
pub use logging::{arcsec_log_fn, arcsec_log_level, arcsec_set_log_callback};
pub use options::{
    ARCSEC_METHOD_QUADS, ARCSEC_METHOD_TETRA, arcsec_cancel_fn, arcsec_progress_fn,
    arcsec_solve_options, arcsec_solve_options_init,
};
pub use result::{
    ARCSEC_SIP_MAX_ORDER, ARCSEC_SIP_SIZE, arcsec_matched_star, arcsec_result,
    arcsec_result_database, arcsec_result_fits_header, arcsec_result_free, arcsec_result_info,
    arcsec_result_matched_stars, arcsec_result_pixel_to_sky, arcsec_result_sky_to_pixel,
    arcsec_result_wcs, arcsec_solve_info, arcsec_wcs,
};
pub use solver::{
    arcsec_solve, arcsec_solve_file, arcsec_solver, arcsec_solver_cancel, arcsec_solver_free,
    arcsec_solver_new,
};

/// Version of the C ABI this library implements.
///
/// It changes only when a change would break a program compiled against an older
/// header: a function removed or its signature changed, a struct field moved or
/// retyped. Adding functions, adding fields at the end of a struct (callers set
/// `struct_size`), and adding status codes or enum values do not change it.
pub const ARCSEC_ABI_VERSION: u32 = 1;

/// The library version, e.g. `"0.4.0"`, as a NUL-terminated static string.
static VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), "\0");

/// The library's version string, e.g. `"0.4.0"` (the arcsec release it was built
/// from). Static; never NULL; do not free.
#[unsafe(no_mangle)]
pub extern "C" fn arcsec_version() -> *const c_char {
    VERSION.as_ptr().cast()
}

/// The ABI version the library implements (`ARCSEC_ABI_VERSION` when it was
/// built). Compare it with the `ARCSEC_ABI_VERSION` your program was compiled
/// with: a different value means the header and the library do not match.
#[unsafe(no_mangle)]
pub extern "C" fn arcsec_abi_version() -> u32 {
    ARCSEC_ABI_VERSION
}
