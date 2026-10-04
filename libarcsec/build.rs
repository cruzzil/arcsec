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
//! - On macOS, sets the dylib's install name to `@rpath/libarcsec.dylib`. rustc
//!   leaves it as the path the library was built at, which no installed copy has;
//!   with `@rpath` a program finds it through its own rpath (`arcsec.pc` adds one).
//! - On MSVC, writes the DLL's debug info to its own PDB in `OUT_DIR`. The
//!   library target is named `arcsec` so that the files are `arcsec.dll` and
//!   `arcsec.lib`, but the CLI binary is `arcsec` too, and both would otherwise
//!   link `target/<profile>/deps/arcsec.pdb`; when a workspace build links them
//!   at the same moment, one fails with LNK1201. Cargo still warns about the
//!   predicted collision (rust-lang/cargo#6313); the files no longer collide.
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
    if matches!(os.as_str(), "macos" | "ios") {
        println!("cargo:rustc-cdylib-link-arg=-Wl,-install_name,@rpath/libarcsec.dylib");
    }
    if std::env::var("CARGO_CFG_TARGET_ENV").as_deref() == Ok("msvc")
        && let Ok(out) = std::env::var("OUT_DIR")
    {
        let pdb = std::path::Path::new(&out).join("arcsec.pdb");
        println!("cargo:rustc-cdylib-link-arg=/PDB:{}", pdb.display());
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
