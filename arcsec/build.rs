//! Embeds a Windows version resource (FileVersion, ProductVersion, ...) in
//! `arcsec.exe`.
//!
//! N.I.N.A. checks the ASTAP executable it is pointed at with
//! `FileVersionInfo.GetVersionInfo(...).FileVersion`, and takes a missing version as
//! an ASTAP older than 0.9.1.0, which it then refuses to run with its default
//! automatic downsampling (`-z 0`). A Rust binary has no version resource unless one
//! is added here, so without this arcsec fails N.I.N.A.'s validation before solving.

fn main() {
    // `cfg(windows)` is the host; the resource compiler it needs (rc.exe from the
    // Windows SDK) only exists there. The target check skips it for anything else.
    #[cfg(windows)]
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("windows") {
        let mut res = winresource::WindowsResource::new();
        res.set("FileDescription", "arcsec astrometric plate solver");
        res.set("OriginalFilename", "arcsec.exe");
        if let Err(e) = res.compile() {
            // A build without the resource still works everywhere except the
            // N.I.N.A. check, so warn rather than fail.
            println!("cargo:warning=could not embed the Windows version resource: {e}");
        }
    }
}
