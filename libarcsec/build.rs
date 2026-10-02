//! Build settings for the C library:
//!
//! - On ELF platforms, gives the shared library the soname `libarcsec.so.<ABI>`, so
//!   a program linked against ABI 1 is never loaded with an incompatible library.
//!   `ABI` must equal `ARCSEC_ABI_VERSION` in src/lib.rs (a test checks).
//! - On ELF platforms, keeps every symbol but the `arcsec_*` API out of the shared
//!   library's dynamic symbol table. rustc exports the `#[no_mangle]` functions of
//!   every dependency from a cdylib, and rsfitsio (arcsec's FITS reader, a Rust
//!   port of CFITSIO) and libbz2-rs-sys define CFITSIO's and bzip2's C API under
//!   their real names. A host that links the real CFITSIO — Siril does — would
//!   otherwise have the two interpose each other at load time. The dependencies
//!   reach the linker as archives (rlibs), so `--exclude-libs` hides them; with
//!   LTO they would not, which is why no profile here enables it.
//! - Passes the target and host triples to the crate's tests, which compile a C
//!   program against the library with the `cc` crate.

/// The ABI version; see `ARCSEC_ABI_VERSION`.
const ABI: u32 = 1;

fn main() {
    let os = std::env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    let elf = matches!(
        os.as_str(),
        "linux"
            | "android"
            | "freebsd"
            | "netbsd"
            | "openbsd"
            | "dragonfly"
            | "illumos"
            | "solaris"
    );
    if elf {
        println!("cargo:rustc-cdylib-link-arg=-Wl,-soname,libarcsec.so.{ABI}");
        println!("cargo:rustc-cdylib-link-arg=-Wl,--exclude-libs=ALL");
    }
    println!("cargo:rustc-env=ARCSEC_BUILD_ABI={ABI}");
    println!(
        "cargo:rustc-env=ARCSEC_BUILD_SONAME={}",
        if elf {
            format!("libarcsec.so.{ABI}")
        } else {
            String::new()
        }
    );
    for var in ["TARGET", "HOST"] {
        if let Ok(v) = std::env::var(var) {
            println!("cargo:rustc-env=ARCSEC_BUILD_{var}={v}");
        }
    }
    println!("cargo:rerun-if-changed=build.rs");
}
