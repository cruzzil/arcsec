//! Where the catalogues are, and what is installed there.

use core::ffi::{c_char, c_int};

use arcsec_core::auto;

use crate::error::{Failure, guard_with};
use crate::util::{path_arg, path_bytes, write_c_string};

/// Write the directory `arcsec catalog install` puts catalogues in
/// (`$ARCSEC_CATALOG_DIR`, else e.g. `~/.local/share/arcsec/catalogs`) into
/// `buf`, as `snprintf` does: at most `len - 1` bytes and a NUL. Returns the
/// path's full length, so a value `>= len` means `buf` was too small (call with
/// `len` 0 to size it). The directory may not exist.
///
/// # Safety
///
/// `buf` must be NULL or point to `len` writable bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn arcsec_default_catalog_dir(buf: *mut c_char, len: usize) -> usize {
    guard_with(0, || {
        let p = auto::default_catalog_dir();
        // SAFETY: forwarded contract.
        Ok(unsafe { write_c_string(&path_bytes(&p), buf, len) })
    })
    .unwrap_or_else(|zero| zero)
}

/// Like `arcsec_default_catalog_dir`, but the directory a solve uses when
/// `catalog_dir` is NULL: the catalogue directory if it holds a star database,
/// else the working directory (`"."`).
///
/// # Safety
///
/// As `arcsec_default_catalog_dir`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn arcsec_default_database_dir(buf: *mut c_char, len: usize) -> usize {
    guard_with(0, || {
        let p = auto::default_db_path();
        // SAFETY: forwarded contract.
        Ok(unsafe { write_c_string(&path_bytes(&p), buf, len) })
    })
    .unwrap_or_else(|zero| zero)
}

/// Whether `dir` (NULL: the default database directory) holds a star database
/// arcsec can solve with: 1 if so, 0 if not, -1 for a bad argument.
///
/// # Safety
///
/// `dir` must be NULL or a NUL-terminated string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn arcsec_has_star_database(dir: *const c_char) -> c_int {
    guard_with(-1, || {
        // SAFETY: forwarded contract.
        let dir = unsafe { path_arg(dir, "dir") }?.unwrap_or_else(auto::default_db_path);
        Ok(c_int::from(auto::has_star_database(&dir)))
    })
    .unwrap_or_else(|e| e)
}

/// Choose the star database for a field `fov_deg` degrees across (the longer
/// side) from those in `dir` (NULL: the default database directory), as a solve
/// does when `database` is NULL, and write its name (e.g. `"d50"`) into `buf` as
/// `arcsec_default_catalog_dir` does. Returns the name's length, 0 if `dir`
/// holds no database (or on a bad argument; see `arcsec_last_error`).
///
/// # Safety
///
/// `dir` must be NULL or a NUL-terminated string; `buf` NULL or `len` writable
/// bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn arcsec_select_database(
    dir: *const c_char,
    fov_deg: f64,
    buf: *mut c_char,
    len: usize,
) -> usize {
    guard_with(0, || {
        if !(fov_deg.is_finite() && fov_deg > 0.0) {
            return Err(Failure::invalid("fov_deg must be positive"));
        }
        // SAFETY: forwarded contract.
        let dir = unsafe { path_arg(dir, "dir") }?.unwrap_or_else(auto::default_db_path);
        Ok(auto::select_db_for_fov(&dir, fov_deg).map_or(0, |name| {
            // SAFETY: forwarded contract.
            unsafe { write_c_string(name.as_bytes(), buf, len) }
        }))
    })
    .unwrap_or_else(|zero| zero)
}

/// Whether a blind index is available at `path` (NULL: the default catalogue
/// directory): 1 for an arcsec index (`*.arcsecix`, the file or one in the
/// directory), 2 for Astrometry.net `index-*.fits` files, 0 for none, -1 for a
/// bad argument. With an index, a solve needs no position hint.
///
/// # Safety
///
/// `path` must be NULL or a NUL-terminated string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn arcsec_has_blind_index(path: *const c_char) -> c_int {
    guard_with(-1, || {
        // SAFETY: forwarded contract.
        let path = unsafe { path_arg(path, "path") }?.unwrap_or_else(auto::default_catalog_dir);
        if auto::find_arcsec_index(&path).is_some() {
            return Ok(1);
        }
        let anet = if path.is_file() {
            path.file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.starts_with("index-") && n.ends_with(".fits"))
        } else {
            !auto::collect_index_files(&path, 0.0).is_empty()
        };
        Ok(if anet { 2 } else { 0 })
    })
    .unwrap_or_else(|e| e)
}
