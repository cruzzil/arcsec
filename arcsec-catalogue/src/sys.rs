//! What the machine has to spare: free disk space and memory, for the checks made
//! before building a blind index or downloading a catalogue.
//!
//! Each probe returns `None` when the platform cannot say, and every caller treats
//! that as "unknown, do not block": a missing number must never stop an install
//! that would have worked.

use std::path::Path;

/// Bytes free to an unprivileged user on the file system holding `path`. A path
/// that does not exist yet (a catalogue directory about to be created) is measured
/// at its nearest existing ancestor.
pub fn free_space(path: &Path) -> Option<u64> {
    let mut p = path;
    loop {
        if p.exists() {
            return free_space_at(p);
        }
        p = p.parent()?;
        if p.as_os_str().is_empty() {
            return free_space_at(Path::new("."));
        }
    }
}

/// Bytes of memory available to a new process without swapping, as best the
/// platform can tell: `MemAvailable` on Linux, free physical memory on Windows,
/// and on macOS (which does not report an equivalent cheaply) total memory.
pub fn available_memory() -> Option<u64> {
    available_memory_impl()
}

#[cfg(unix)]
fn free_space_at(p: &Path) -> Option<u64> {
    use std::os::unix::ffi::OsStrExt as _;
    let c = alloc::ffi::CString::new(p.as_os_str().as_bytes()).ok()?;
    let mut st = core::mem::MaybeUninit::<libc::statvfs>::uninit();
    // Safety: `c` is a valid NUL-terminated path and `st` is a properly sized
    // out-parameter, read only after statvfs reports success.
    let rc = unsafe { libc::statvfs(c.as_ptr(), st.as_mut_ptr()) };
    if rc != 0 {
        return None;
    }
    let st = unsafe { st.assume_init() };
    #[allow(clippy::useless_conversion)] // the field widths differ between platforms
    let bytes = u64::from(st.f_bavail).saturating_mul(u64::from(st.f_frsize));
    Some(bytes)
}

#[cfg(windows)]
fn free_space_at(p: &Path) -> Option<u64> {
    use std::os::windows::ffi::OsStrExt as _;
    use windows_sys::Win32::Storage::FileSystem::GetDiskFreeSpaceExW;
    let wide: Vec<u16> = p.as_os_str().encode_wide().chain([0]).collect();
    let mut avail = 0u64;
    // Safety: `wide` is NUL-terminated; the out-pointers are valid or null.
    let ok = unsafe {
        GetDiskFreeSpaceExW(
            wide.as_ptr(),
            &raw mut avail,
            core::ptr::null_mut(),
            core::ptr::null_mut(),
        )
    };
    (ok != 0).then_some(avail)
}

#[cfg(not(any(unix, windows)))]
fn free_space_at(_: &Path) -> Option<u64> {
    None
}

#[cfg(target_os = "linux")]
fn available_memory_impl() -> Option<u64> {
    let s = std::fs::read_to_string("/proc/meminfo").ok()?;
    parse_meminfo(&s)
}

#[cfg(target_os = "macos")]
fn available_memory_impl() -> Option<u64> {
    let mut v: u64 = 0;
    let mut len = core::mem::size_of::<u64>();
    // Safety: the name is NUL-terminated and `v`/`len` describe a u64 buffer.
    let rc = unsafe {
        libc::sysctlbyname(
            c"hw.memsize".as_ptr(),
            (&raw mut v).cast(),
            &raw mut len,
            core::ptr::null_mut(),
            0,
        )
    };
    (rc == 0 && v > 0).then_some(v)
}

#[cfg(windows)]
fn available_memory_impl() -> Option<u64> {
    use windows_sys::Win32::System::SystemInformation::{GlobalMemoryStatusEx, MEMORYSTATUSEX};
    // Safety: MEMORYSTATUSEX is plain data; dwLength must be set before the call.
    let mut m: MEMORYSTATUSEX = unsafe { core::mem::zeroed() };
    m.dwLength = core::mem::size_of::<MEMORYSTATUSEX>() as u32;
    let ok = unsafe { GlobalMemoryStatusEx(&raw mut m) };
    (ok != 0).then_some(m.ullAvailPhys)
}

#[cfg(not(any(target_os = "linux", target_os = "macos", windows)))]
fn available_memory_impl() -> Option<u64> {
    None
}

/// `MemAvailable` from the text of `/proc/meminfo`, in bytes.
#[cfg_attr(not(target_os = "linux"), allow(dead_code))]
fn parse_meminfo(s: &str) -> Option<u64> {
    let line = s.lines().find(|l| l.starts_with("MemAvailable:"))?;
    let kb: u64 = line.split_whitespace().nth(1)?.parse().ok()?;
    Some(kb * 1024)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn meminfo_parses() {
        let s = "MemTotal:       32768000 kB\nMemFree:          100000 kB\nMemAvailable:   25000000 kB\n";
        assert_eq!(parse_meminfo(s), Some(25_000_000 * 1024));
        assert_eq!(parse_meminfo("MemTotal: 1 kB\n"), None);
    }

    #[test]
    fn free_space_walks_up_to_an_existing_directory() {
        let tmp = std::env::temp_dir();
        let here = free_space(&tmp);
        if cfg!(any(unix, windows)) {
            assert!(here.is_some_and(|b| b > 0));
        }
        let deep = tmp.join("arcsec-no-such-dir").join("a").join("b");
        assert_eq!(free_space(&deep).is_some(), here.is_some());
    }
}
