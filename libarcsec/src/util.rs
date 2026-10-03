//! Reading caller memory: versioned structs, strings, paths; and writing strings
//! back into caller buffers.

use core::ffi::{CStr, c_char};
use std::path::PathBuf;

use crate::error::{Failure, Outcome};

/// A `#[repr(C)]` struct that starts with a `size_t struct_size` and whose every
/// field is valid for any bit pattern (integers, floats, raw and function
/// pointers; no `bool`, no Rust enums), so a struct partly copied from the caller
/// and partly from its defaults is always a valid value.
///
/// # Safety
///
/// Implementors must be `#[repr(C)]`, begin with a `usize` `struct_size` field,
/// and contain only fields for which every bit pattern is valid.
pub(crate) unsafe trait Versioned: Copy {
    /// The value to start from: the defaults for fields an older caller's struct
    /// does not have.
    fn defaults() -> Self;
}

/// Copy a caller's versioned struct.
///
/// Reads `struct_size` first, then copies the smaller of it and our size over
/// [`Versioned::defaults`]: a caller built against an older header (smaller
/// struct) gets defaults for the fields it lacks, and a newer one (larger struct)
/// has the fields we do not know about ignored. A `struct_size` smaller than
/// `min_size` (the first ABI's size) means the caller did not initialise it.
///
/// # Safety
///
/// `ptr` must be NULL or point to at least `struct_size` readable bytes (it may
/// be unaligned).
pub(crate) unsafe fn read_versioned<T: Versioned>(
    ptr: *const T,
    what: &str,
    min_size: usize,
) -> Outcome<T> {
    if ptr.is_null() {
        return Err(Failure::invalid(format!("{what} is NULL")));
    }
    // SAFETY: the caller guarantees at least the leading size_t is readable.
    let size = unsafe { ptr.cast::<usize>().read_unaligned() };
    if size < min_size {
        return Err(Failure::invalid(format!(
            "{what}.struct_size is {size}; set it to sizeof({what}) ({}), or call the init function",
            core::mem::size_of::<T>()
        )));
    }
    let mut out = T::defaults();
    let n = size.min(core::mem::size_of::<T>());
    // SAFETY: n bytes are readable at ptr (n <= struct_size) and writable at out
    // (n <= size_of::<T>()); the regions cannot overlap (out is a local); and any
    // bit pattern is a valid T (the Versioned contract).
    unsafe {
        core::ptr::copy_nonoverlapping(
            ptr.cast::<u8>(),
            core::ptr::from_mut(&mut out).cast::<u8>(),
            n,
        );
    }
    Ok(out)
}

/// Copy `value` into a caller's versioned output struct, whose `struct_size` the
/// caller has set: as much of it as both sides know about. `struct_size` itself
/// is left as the caller set it.
///
/// # Safety
///
/// `out` must be NULL or point to at least `struct_size` writable bytes.
pub(crate) unsafe fn write_versioned<T: Versioned>(
    out: *mut T,
    value: &T,
    what: &str,
    min_size: usize,
) -> Outcome<()> {
    if out.is_null() {
        return Err(Failure::invalid(format!("{what} output is NULL")));
    }
    // SAFETY: the caller guarantees the leading size_t is readable.
    let size = unsafe { out.cast::<usize>().read_unaligned() };
    if size < min_size {
        return Err(Failure::invalid(format!(
            "{what}.struct_size is {size}; set it to sizeof({what}) ({}) before the call",
            core::mem::size_of::<T>()
        )));
    }
    let n = size.min(core::mem::size_of::<T>());
    let skip = core::mem::size_of::<usize>();
    // SAFETY: bytes skip..n lie within both the caller's struct_size and our T;
    // value is a separate local, so the regions do not overlap.
    unsafe {
        core::ptr::copy_nonoverlapping(
            core::ptr::from_ref(value).cast::<u8>().add(skip),
            out.cast::<u8>().add(skip),
            n - skip,
        );
    }
    Ok(())
}

/// A NUL-terminated UTF-8 string argument; `None` for NULL.
///
/// # Safety
///
/// `p` must be NULL or point to a NUL-terminated string.
pub(crate) unsafe fn str_arg<'a>(p: *const c_char, what: &str) -> Outcome<Option<&'a str>> {
    if p.is_null() {
        return Ok(None);
    }
    // SAFETY: non-NULL and NUL-terminated per the contract.
    let c = unsafe { CStr::from_ptr(p) };
    c.to_str()
        .map(Some)
        .map_err(|_| Failure::invalid(format!("{what} is not valid UTF-8")))
}

/// A path argument; `None` for NULL or an empty string.
///
/// On Unix the bytes are taken as they are, so a path that is not UTF-8 still
/// works; elsewhere it must be UTF-8.
///
/// # Safety
///
/// `p` must be NULL or point to a NUL-terminated string.
pub(crate) unsafe fn path_arg(p: *const c_char, what: &str) -> Outcome<Option<PathBuf>> {
    if p.is_null() {
        return Ok(None);
    }
    // SAFETY: non-NULL and NUL-terminated per the contract.
    let c = unsafe { CStr::from_ptr(p) };
    if c.is_empty() {
        return Ok(None);
    }
    #[cfg(unix)]
    {
        use std::os::unix::ffi::OsStrExt as _;
        let _ = what;
        Ok(Some(PathBuf::from(std::ffi::OsStr::from_bytes(
            c.to_bytes(),
        ))))
    }
    #[cfg(not(unix))]
    {
        c.to_str()
            .map(|s| Some(PathBuf::from(s)))
            .map_err(|_| Failure::invalid(format!("{what} is not valid UTF-8")))
    }
}

/// The bytes of a path for C: exact on Unix, UTF-8 elsewhere.
pub(crate) fn path_bytes(p: &std::path::Path) -> Vec<u8> {
    #[cfg(unix)]
    {
        use std::os::unix::ffi::OsStrExt as _;
        p.as_os_str().as_bytes().to_vec()
    }
    #[cfg(not(unix))]
    {
        p.to_string_lossy().into_owned().into_bytes()
    }
}

/// Write `bytes` into a caller buffer as `snprintf` would: at most `len - 1`
/// bytes and a NUL (nothing at all when `buf` is NULL or `len` is 0). Returns the
/// full length without the NUL, so a return value `>= len` means the buffer was
/// too small and the output was cut.
///
/// # Safety
///
/// `buf` must be NULL or point to `len` writable bytes.
pub(crate) unsafe fn write_c_string(bytes: &[u8], buf: *mut c_char, len: usize) -> usize {
    if !buf.is_null() && len > 0 {
        let n = bytes.len().min(len - 1);
        // SAFETY: n + 1 <= len bytes are writable at buf; bytes is Rust memory
        // that cannot overlap the caller's buffer.
        unsafe {
            core::ptr::copy_nonoverlapping(bytes.as_ptr(), buf.cast::<u8>(), n);
            buf.add(n).write(0);
        }
    }
    bytes.len()
}

/// Write `value` through an optional out-pointer.
///
/// # Safety
///
/// `out` must be NULL or valid for a write of `T`.
pub(crate) unsafe fn put<T>(out: *mut T, value: T) {
    if !out.is_null() {
        // SAFETY: non-NULL and writable per the contract; unaligned is allowed.
        unsafe { out.write_unaligned(value) };
    }
}
