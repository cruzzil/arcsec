//! The C test: compile `tests/c/smoke.c` against the built library, the way a C
//! program would use it, and run it on a synthetic field.
//!
//! This is an integration test (the repository otherwise keeps its tests in-file)
//! because it needs the library built as a C library, which cargo does for an
//! integration test's dependencies and not for unit tests.
//!
//! Environment, for the sanitizer and valgrind runs (`tests/c/run-checks.sh`):
//!
//! - `ARCSEC_C_CFLAGS`: extra compiler flags, e.g. `-fsanitize=address`.
//! - `ARCSEC_C_RUNNER`: a command to run the program under, e.g.
//!   `valgrind --error-exitcode=1 --leak-check=full`.
//! - `ARCSEC_C_LINK=static`: link `libarcsec.a` instead of the shared library.
//! - `ARCSEC_C_KEEP=<dir>`: leave the program and its fixture in `<dir>`.

use std::path::{Path, PathBuf};
use std::process::Command;

use arcsec_core::test_support::{TempDir, TruthWcs, fits_f32_bytes, synthetic_field};

/// The directory cargo put the library in: the test executable is
/// `<target>/<profile>/deps/c_api-<hash>`.
fn profile_dir() -> PathBuf {
    let exe = std::env::current_exe().unwrap();
    exe.parent().and_then(Path::parent).unwrap().to_path_buf()
}

/// The first of `names` found in the profile directory's `deps`, or in the
/// profile directory itself. `deps` first: that is where the build for this test
/// run puts the library, while the copy beside it is only refreshed by `cargo
/// build` and may be stale.
fn find(names: &[&str]) -> Option<PathBuf> {
    let dir = profile_dir();
    [dir.join("deps"), dir]
        .iter()
        .flat_map(|d| names.iter().map(move |n| d.join(n)))
        .find(|p| p.is_file())
}

/// On ELF platforms the shared library exports the C API and nothing else: in
/// particular not the CFITSIO functions of the Rust FITS reader inside it, which a
/// host linking the real CFITSIO would otherwise pick up (see build.rs).
#[test]
fn only_the_api_is_exported() {
    if env!("ARCSEC_BUILD_SONAME").is_empty() {
        return;
    }
    let lib = find(&["libarcsec.so"]).expect("the library was not built");
    let Ok(out) = Command::new("nm")
        .args(["-D", "--defined-only"])
        .arg(&lib)
        .output()
    else {
        eprintln!("nm is not available; not checked");
        return;
    };
    let text = String::from_utf8_lossy(&out.stdout);
    let names: Vec<&str> = text
        .lines()
        .filter_map(|l| l.split_whitespace().nth(2))
        .collect();
    assert!(names.contains(&"arcsec_solve"), "{text}");
    let foreign: Vec<&&str> = names.iter().filter(|n| !n.starts_with("arcsec_")).collect();
    assert!(foreign.is_empty(), "exported beyond the API: {foreign:?}");
}

#[test]
fn the_c_program_passes() {
    let keep = std::env::var_os("ARCSEC_C_KEEP").map(PathBuf::from);
    let tmp = TempDir::new("capi-c");
    let work = keep.clone().unwrap_or_else(|| tmp.path().to_path_buf());
    std::fs::create_dir_all(&work).unwrap();

    // ── The fixture ───────────────────────────────────────────────────────────
    let (ra, dec, scale) = (84.3_f64, -5.2_f64, 5.0_f64);
    let truth = TruthWcs::new(
        ra.to_radians(),
        dec.to_radians(),
        scale,
        23.0,
        false,
        400,
        320,
    );
    let catalog = work.join("catalog");
    std::fs::create_dir_all(&catalog).unwrap();
    let img = synthetic_field(&catalog, "d50", &truth, 130, 1);
    let raw: Vec<u8> = img.data.iter().flat_map(|v| v.to_ne_bytes()).collect();
    std::fs::write(work.join("field.f32"), raw).unwrap();
    std::fs::write(
        work.join("field.fits"),
        fits_f32_bytes(
            &img,
            &[
                ("RA", "84.45"),
                ("DEC", "-5.25"),
                ("FOCALLEN", "206.265"),
                ("XPIXSZ", "5.0"),
            ],
        ),
    )
    .unwrap();

    // ── Compile ───────────────────────────────────────────────────────────────
    let manifest = Path::new(env!("CARGO_MANIFEST_DIR"));
    let src = manifest.join("tests").join("c").join("smoke.c");
    let include = manifest.join("include");
    let target = env!("ARCSEC_BUILD_TARGET");
    let tool = cc::Build::new()
        .target(target)
        .host(env!("ARCSEC_BUILD_HOST"))
        .opt_level(1)
        .debug(true)
        .cargo_metadata(false)
        .cargo_warnings(false)
        .try_get_compiler()
        .expect("no C compiler");
    let exe = work.join(if cfg!(windows) { "smoke.exe" } else { "smoke" });
    let static_link = std::env::var("ARCSEC_C_LINK").is_ok_and(|v| v == "static");
    let mut cc = tool.to_command();
    cc.current_dir(&work);
    let lib_dir;
    if tool.is_like_msvc() {
        let lib = if static_link {
            find(&["arcsec.lib"])
        } else {
            find(&["arcsec.dll.lib"])
        }
        .expect("the library was not built");
        lib_dir = lib.parent().unwrap().to_path_buf();
        cc.arg("/nologo")
            .arg("/W3")
            .arg(format!("/I{}", include.display()))
            .arg(&src)
            .arg(format!("/Fe{}", exe.display()))
            .arg(&lib);
        if static_link {
            // What a Rust static library needs from the system on Windows.
            cc.args([
                "ws2_32.lib",
                "userenv.lib",
                "ntdll.lib",
                "bcrypt.lib",
                "advapi32.lib",
            ]);
        }
    } else {
        let lib = if static_link {
            find(&["libarcsec.a"])
        } else if cfg!(target_os = "macos") {
            find(&["libarcsec.dylib"])
        } else if cfg!(windows) {
            find(&["libarcsec.dll.a", "arcsec.dll.lib"])
        } else {
            find(&["libarcsec.so"])
        }
        .expect("the library was not built");
        // On ELF platforms the library's soname is libarcsec.so.<ABI>, which is
        // the name the program will look for: install it under that name, as a
        // package would.
        let soname = env!("ARCSEC_BUILD_SONAME");
        let lib = if static_link || soname.is_empty() {
            lib
        } else {
            let installed = work.join(soname);
            std::fs::copy(&lib, &installed).unwrap();
            installed
        };
        lib_dir = lib.parent().unwrap().to_path_buf();
        cc.args(["-std=c99", "-Wall", "-Wextra", "-Werror"])
            .arg(format!("-I{}", include.display()))
            .arg(&src)
            .arg("-o")
            .arg(&exe)
            .arg(&lib);
        if !static_link && !cfg!(windows) {
            cc.arg(format!("-Wl,-rpath,{}", lib_dir.display()));
        }
        if static_link {
            if cfg!(target_os = "macos") {
                cc.args(["-framework", "CoreFoundation", "-framework", "Security"]);
            } else if !cfg!(windows) {
                cc.args(["-lpthread", "-ldl"]);
            }
        }
        cc.arg("-lm");
    }
    if let Ok(flags) = std::env::var("ARCSEC_C_CFLAGS") {
        cc.args(flags.split_whitespace());
    }
    let out = cc.output().expect("could not run the C compiler");
    assert!(
        out.status.success(),
        "compiling smoke.c failed:\n{:?}\n{}{}",
        cc,
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );

    // ── Run ───────────────────────────────────────────────────────────────────
    let runner: Vec<String> = std::env::var("ARCSEC_C_RUNNER")
        .map(|r| r.split_whitespace().map(String::from).collect())
        .unwrap_or_default();
    let mut run = match runner.split_first() {
        Some((first, rest)) => {
            let mut c = Command::new(first);
            c.args(rest).arg(&exe);
            c
        }
        None => Command::new(&exe),
    };
    run.args([
        catalog.to_str().unwrap(),
        work.join("field.f32").to_str().unwrap(),
        work.join("field.fits").to_str().unwrap(),
        "400",
        "320",
        &ra.to_string(),
        &dec.to_string(),
        &scale.to_string(),
    ]);
    if cfg!(windows) {
        // The DLL is found on PATH.
        let path = std::env::var_os("PATH").unwrap_or_default();
        let mut dirs = vec![lib_dir.clone()];
        dirs.extend(std::env::split_paths(&path));
        run.env("PATH", std::env::join_paths(dirs).unwrap());
    }
    let out = run.output().expect("could not run the C program");
    let stdout = String::from_utf8_lossy(&out.stdout);
    let stderr = String::from_utf8_lossy(&out.stderr);
    println!("{stdout}");
    eprintln!("{stderr}");
    assert!(
        out.status.success(),
        "the C program failed ({}):\n{stdout}\n{stderr}",
        out.status
    );
    assert!(stdout.contains("all C API checks passed"));
}
