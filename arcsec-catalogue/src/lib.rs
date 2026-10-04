//! Catalogue management for the [arcsec] plate solver.
//!
//! arcsec solves against star databases (ASTAP's D05 … W08) and, for blind solving,
//! a blind index built from one of them or Astrometry.net index files. This crate is
//! everything about those files that is not solving:
//!
//! - what can be installed, and how to recognise it on disk ([`registry`]);
//! - where catalogues live ([`default_dir`]) and which database suits a field
//!   ([`select_db_for_fov`], re-exported from [`arcsec_core::auto`], which the
//!   solver itself uses);
//! - downloading and unpacking them (`install`, `fetch`; the `download`
//!   feature);
//! - removing and checking them ([`remove`], [`verify`]);
//! - planning, costing and building arcsec's own blind index ([`index`]).
//!
//! Nothing here prints or exits: results and errors are typed. Long operations
//! (a download, an index build) report through a callback, and stop early when the
//! thread's [`arcsec_core::cancel`] token is cancelled, as a solve does.
//!
//! ```no_run
//! use arcsec_catalogue::{default_dir, registry, select_db_for_fov};
//!
//! let dir = default_dir();
//! for e in registry::REGISTRY {
//!     let state = if registry::is_installed(&dir, e) { "installed" } else { "-" };
//!     println!("{:<10} {:<10} {state}", e.id, e.purpose.label());
//! }
//! println!("for a 1.5° field: {:?}", select_db_for_fov(&dir, 1.5));
//! ```
//!
//! # Features
//!
//! - `download` (default): `fetch` and `install`, over HTTPS with `ureq` and
//!   `rustls`, and zip, tar and xz extraction. Without it the crate pulls in no
//!   network or TLS code; everything else still works.
//!
//! The `arcsec` command line is built on this crate; its API follows that user and
//! may change between minor versions.
//!
//! [arcsec]: https://github.com/cruzzil/arcsec

// `alloc` is not in the extern prelude for a crate that links std, so it has to
// be declared before alloc:: paths can be written.
extern crate alloc;

mod error;
#[cfg(feature = "download")]
pub mod fetch;
mod format;
pub mod index;
#[cfg(feature = "download")]
mod install;
mod manage;
pub mod registry;
pub mod sys;

/// A temporary directory and a tiny star database, for tests: compiled for this
/// crate's own tests, and with the `test-support` feature for the command line's.
/// Not a stable API.
#[cfg(any(test, feature = "test-support"))]
#[doc(hidden)]
pub mod test_support;

pub use error::{Error, Result};
pub use format::{duration, gigabytes, human_bytes};
#[cfg(feature = "download")]
pub use install::{InstallEvent, InstallOptions, Installed, install};
pub use manage::{
    CatalogueCheck, IndexCheck, IndexHealth, Problem, Removed, Verification, download_disk,
    recommend, remove, verify,
};

/// Where catalogues are kept and the solver looks by default: `$ARCSEC_CATALOG_DIR`,
/// else the platform's data directory (`~/.local/share/arcsec/catalogs` on Linux).
/// See [`arcsec_core::auto::default_catalog_dir`] for the full rules. The directory
/// may not exist yet.
#[must_use]
pub fn default_dir() -> std::path::PathBuf {
    arcsec_core::auto::default_catalog_dir()
}

/// Choosing a star database. These live in [`arcsec_core::auto`] because a solve
/// makes the same choice when it is not told; they are re-exported here so that a
/// program managing catalogues finds them in one place.
pub use arcsec_core::auto::{
    ASTAP_EXTS, DB_FOV_RANGES, available_dbs, default_db_path, has_star_database, select_db_for_fov,
};
