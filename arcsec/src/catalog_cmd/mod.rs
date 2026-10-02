//! `arcsec catalog` — install and manage star catalogues.
//!
//! arcsec needs a star catalogue to solve anything, and a *photometric* catalogue to
//! do colour calibration. Both live behind download pages that are easy to get wrong,
//! so this puts them one command away and in one place:
//!
//! ```text
//! arcsec catalog list                  # what exists, what is installed
//! arcsec catalog recommend --fov 1.5   # what this rig needs
//! arcsec catalog install d50 v05       # fetch and unpack
//! arcsec catalog path                  # where they went
//! ```
//!
//! Everything lands in a per-platform data directory (see [`default_dir`]) which the
//! solver reads by default, so `-d` is only needed to override it.

mod fetch;
pub mod index_cmd;
mod registry;

use std::ffi::OsString;
use std::fs;
use std::io::Write as _;
use std::path::{Path, PathBuf};

use registry::{
    Archive, Entry, Purpose, astap_file_count, expected_file_count, files_of, find, installed_size,
    is_complete, loose_files,
};
use registry::{REGISTRY, is_installed};

use arcsec_io::image_io;

use crate::cli;

/// Parse and run `arcsec catalog <subcommand>`; returns the process exit code.
///
/// `args` starts at `catalog` itself, i.e. the command line without `arcsec`.
pub fn run(args: impl IntoIterator<Item = OsString>) -> i32 {
    let m = cli::catalog_command().get_matches_from(args);

    let dir = m
        .get_one::<PathBuf>("dir")
        .cloned()
        .unwrap_or_else(default_dir);

    let names = |sm: &clap::ArgMatches| -> Vec<String> {
        sm.get_many::<String>("names")
            .map_or_else(Vec::new, |v| v.cloned().collect())
    };

    let result = match m.subcommand() {
        Some(("list", _)) => {
            cmd_list(&dir);
            Ok(())
        }
        Some(("path", _)) => {
            println!("{}", dir.display());
            Ok(())
        }
        Some(("recommend", sm)) => {
            let fov = match (sm.get_one::<f64>("fov"), sm.get_one::<PathBuf>("like")) {
                (Some(f), _) => Some(*f),
                (None, Some(img)) => image_io::read_pixel_scale(img).and_then(|ps| {
                    image_io::read_dimensions(img).map(|(w, h)| ps * f64::from(w.max(h)) / 3600.0)
                }),
                _ => None,
            };
            match fov {
                Some(f) if f > 0.0 => {
                    cmd_recommend(&dir, f, sm.get_flag("photometry"));
                    Ok(())
                }
                _ => Err(
                    "give --fov <degrees>, or --like <image> with FOCALLEN and XPIXSZ in its header"
                        .to_string(),
                ),
            }
        }
        Some(("install", sm)) => {
            cmd_install(&dir, &names(sm), sm.get_flag("yes"), sm.get_flag("keep"))
        }
        Some(("remove", sm)) => cmd_remove(&dir, &names(sm), sm.get_flag("yes")),
        Some(("verify", _)) => cmd_verify(&dir),
        Some(("index", sm)) => match sm.subcommand() {
            Some(("build", b)) => index_cmd::cmd_build(
                &dir,
                b.get_one::<PathBuf>("db"),
                b.get_one::<String>("name"),
                *b.get_one::<f64>("min-fov").unwrap_or(&0.3),
                *b.get_one::<f64>("max-fov").unwrap_or(&30.0),
                b.get_one::<PathBuf>("out"),
                *b.get_one::<usize>("threads").unwrap_or(&0),
            ),
            Some(("info", i)) => index_cmd::cmd_info(&dir, i.get_one::<PathBuf>("file")),
            _ => Err("unknown index subcommand".to_string()),
        },
        _ => Err("unknown subcommand".to_string()),
    };

    match result {
        Ok(()) => 0,
        Err(e) => {
            eprintln!("Error: {e}");
            1
        }
    }
}

// ── Install location ────────────────────────────────────────────────────────────

/// Where catalogues are kept: `$ARCSEC_CATALOG_DIR`, else the platform's data
/// directory. The rules live in [`arcsec_core::auto::default_catalog_dir`], which
/// the solver (and the C library) use too.
pub fn default_dir() -> PathBuf {
    arcsec_core::auto::default_catalog_dir()
}

/// Bytes, in decimal units, to one decimal place: `901.3 MB`.
pub fn human(bytes: u64) -> String {
    const U: [&str; 5] = ["B", "kB", "MB", "GB", "TB"];
    let mut v = bytes as f64;
    let mut i = 0;
    while v >= 1000.0 && i < U.len() - 1 {
        v /= 1000.0;
        i += 1;
    }
    if i == 0 {
        format!("{bytes} {}", U[i])
    } else {
        format!("{v:.1} {}", U[i])
    }
}

/// `"  [installed]"` if `e` is in `dir`, else nothing.
fn installed_tag(dir: &Path, e: &Entry) -> &'static str {
    if is_installed(dir, e) {
        "  [installed]"
    } else {
        ""
    }
}

// ── Commands ────────────────────────────────────────────────────────────────────

fn cmd_list(dir: &Path) {
    println!("Catalogue directory: {}", dir.display());
    println!("  (override with --dir, or the ARCSEC_CATALOG_DIR environment variable)\n");
    println!(
        "{:<11} {:<11} {:>9}  {:<13} STATUS",
        "NAME", "PURPOSE", "DOWNLOAD", "FIELDS"
    );
    for e in REGISTRY {
        let fov = e
            .fov
            .map_or_else(|| "—".to_string(), |(lo, hi)| format!("{lo}°–{hi}°"));
        let status = if is_installed(dir, e) {
            let sz = installed_size(dir, e);
            if sz > 0 {
                format!("installed ({})", human(sz))
            } else {
                "installed".to_string()
            }
        } else {
            "not installed".to_string()
        };
        println!(
            "{:<11} {:<11} {:>9}  {fov:<13} {status}",
            e.id,
            e.purpose.label(),
            human(e.bytes),
        );
    }
    println!("\nDescriptions:");
    for e in REGISTRY {
        println!("  {:<11} {}", e.id, e.desc);
    }
    let indexes = index_cmd::index_files(dir);
    println!("\nBlind indexes (built locally with `arcsec catalog index build`):");
    if indexes.is_empty() {
        println!("  none");
    }
    for p in indexes {
        if let Err(e) = index_cmd::describe(&p) {
            println!("  {}: {e}", p.display());
        }
    }
}

/// Suggest catalogues for a field size.
fn cmd_recommend(dir: &Path, fov_deg: f64, want_photometry: bool) {
    println!("For a {fov_deg:.2}° field:\n");
    let pick = |purpose: Purpose| -> Option<&'static Entry> {
        // Prefer the smallest download whose range covers the field, so the advice
        // does not push a gigabyte on someone who does not need it.
        REGISTRY
            .iter()
            .filter(|e| e.purpose == purpose)
            .filter(|e| matches!(e.fov, Some((lo, hi)) if fov_deg >= lo && fov_deg <= hi))
            .min_by_key(|e| e.bytes)
    };

    match pick(Purpose::Solving) {
        Some(e) => println!(
            "  solving      {:<10} {:>9}  {}{}",
            e.id,
            human(e.bytes),
            e.desc,
            installed_tag(dir, e)
        ),
        None => println!("  solving      no catalogue covers this field size"),
    }
    if want_photometry {
        match pick(Purpose::Photometry) {
            Some(e) => println!(
                "  photometry   {:<10} {:>9}  {}{}",
                e.id,
                human(e.bytes),
                e.desc,
                installed_tag(dir, e)
            ),
            None => println!("  photometry   no catalogue covers this field size"),
        }
    }
    match pick(Purpose::BlindIndex) {
        Some(e) => println!(
            "  blind        {:<10} {:>9}  optional: solve with no position hint{}",
            e.id,
            human(e.bytes),
            installed_tag(dir, e)
        ),
        None => println!("  blind        no index set covers this field size"),
    }

    let ids: Vec<&str> = [
        pick(Purpose::Solving),
        want_photometry.then(|| pick(Purpose::Photometry)).flatten(),
    ]
    .into_iter()
    .flatten()
    .filter(|e| !is_installed(dir, e))
    .map(|e| e.id)
    .collect();
    if !ids.is_empty() {
        println!("\n  arcsec catalog install {}", ids.join(" "));
    }
}

/// Ask on stderr and read a yes/no answer from stdin; anything but `y`/`yes`,
/// including end of input, is a no.
fn confirm() -> bool {
    eprint!("Continue? [y/N] ");
    let _ = std::io::stderr().flush();
    let mut line = String::new();
    std::io::stdin().read_line(&mut line).is_ok()
        && matches!(line.trim().to_ascii_lowercase().as_str(), "y" | "yes")
}

fn cmd_remove(dir: &Path, ids: &[String], assume_yes: bool) -> Result<(), String> {
    let mut wanted = Vec::new();
    for id in ids {
        let e = find(id)
            .ok_or_else(|| format!("unknown catalogue '{id}' (try: arcsec catalog list)"))?;
        if !wanted.iter().any(|w: &&Entry| w.id == e.id) {
            wanted.push(e);
        }
    }

    // The catalogue directory may be shared with ASTAP (docs/catalogues.md §6), in
    // which case these are ASTAP's files too, so say exactly what will go first.
    println!("Removing from {}\n", dir.display());
    for e in &wanted {
        let files = files_of(dir, e);
        let bytes: u64 = files
            .iter()
            .filter_map(|p| fs::metadata(p).ok())
            .map(|m| m.len())
            .sum();
        println!(
            "  {:<11} {:>5} files  {:>9}",
            e.id,
            files.len(),
            human(bytes)
        );
    }
    println!();
    if !assume_yes && !confirm() {
        println!("Cancelled.");
        return Ok(());
    }

    for e in wanted {
        let files = files_of(dir, e);
        let mut freed = 0u64;
        for p in &files {
            freed += fs::metadata(p).map_or(0, |m| m.len());
            fs::remove_file(p).map_err(|err| format!("{}: {err}", p.display()))?;
        }
        println!(
            "removed {}: {} files, {} freed",
            e.id,
            files.len(),
            human(freed)
        );
    }
    Ok(())
}

/// Check that every installed catalogue looks structurally sound.
fn cmd_verify(dir: &Path) -> Result<(), String> {
    let mut problems = 0;
    let mut checked = 0;
    for e in REGISTRY {
        if !is_installed(dir, e) {
            continue;
        }
        checked += 1;
        let files = files_of(dir, e);
        let mut bad = Vec::new();

        for p in &files {
            match fs::metadata(p) {
                // Every format has at least a 110-byte header (or a FITS block).
                Ok(m) if m.len() < 120 => bad.push(format!("{} is truncated", p.display())),
                Err(err) => bad.push(format!("{}: {err}", p.display())),
                _ => {}
            }
        }

        // The file count is fixed: by the grid for an ASTAP database, by the set
        // for downloaded indexes.
        if let Some(want) = astap_file_count(dir, e).or_else(|| expected_file_count(e))
            && files.len() != want
        {
            bad.push(format!("expected {want} files, found {}", files.len()));
        }

        if bad.is_empty() {
            println!(
                "  {:<11} ok  ({} files, {})",
                e.id,
                files.len(),
                human(installed_size(dir, e))
            );
        } else {
            problems += bad.len();
            println!("  {:<11} PROBLEMS:", e.id);
            for b in bad {
                println!("      {b}");
            }
        }
    }
    for p in index_cmd::index_files(dir) {
        checked += 1;
        let name = p
            .file_name()
            .map_or_else(String::new, |n| n.to_string_lossy().into_owned());
        let res = arcsec_core::index::BlindIndex::open(&p).and_then(|ix| {
            ix.validate()?;
            Ok(ix)
        });
        match res {
            Ok(ix) => println!(
                "  {name:<11} ok  ({} patterns, {})",
                ix.n_patterns(),
                human(ix.file_size() as u64)
            ),
            Err(e) => {
                problems += 1;
                println!(
                    "  {name:<11} PROBLEM: {e} - rebuild it with `arcsec catalog index build`"
                );
            }
        }
    }
    if checked == 0 {
        println!("No catalogues installed in {}", dir.display());
        return Ok(());
    }
    if problems > 0 {
        return Err(format!(
            "{problems} problem(s) found; re-run `arcsec catalog install <name>`"
        ));
    }
    Ok(())
}

fn cmd_install(dir: &Path, ids: &[String], assume_yes: bool, keep: bool) -> Result<(), String> {
    let mut wanted: Vec<&'static Entry> = Vec::new();
    for id in ids {
        let e = find(id).ok_or_else(|| {
            format!("unknown catalogue '{id}'. Run `arcsec catalog list` to see the options.")
        })?;
        if is_complete(dir, e) {
            println!("{}: already installed, skipping", e.id);
            continue;
        }
        if !wanted.iter().any(|w| w.id == e.id) {
            wanted.push(e);
        }
    }
    if wanted.is_empty() {
        return Ok(());
    }

    let total: u64 = wanted.iter().map(|e| e.bytes).sum();
    println!("Installing into {}\n", dir.display());
    for e in &wanted {
        println!("  {:<11} {:>9}  {}", e.id, human(e.bytes), e.desc);
    }
    println!("\nTotal download: {}", human(total));
    if !assume_yes && !confirm() {
        println!("Cancelled.");
        return Ok(());
    }

    fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;

    for e in wanted {
        println!("\n{} — {}", e.id, e.desc);
        match e.archive {
            Archive::Loose => install_loose(dir, e)?,
            Archive::Zip | Archive::Deb => install_archive(dir, e, keep)?,
        }
        if is_installed(dir, e) {
            println!("  {} installed ({})", e.id, human(installed_size(dir, e)));
        } else {
            return Err(format!(
                "{}: extraction finished but no catalogue files appeared in {}",
                e.id,
                dir.display()
            ));
        }
    }
    println!("\nDone. The solver uses {} by default.", dir.display());
    Ok(())
}

/// Download each file of a `Loose` set that is not already present.
///
/// A file that fails is reported and the rest continue; the set as a whole then
/// fails, so a script sees it, and a re-run fetches only what is missing.
fn install_loose(dir: &Path, e: &Entry) -> Result<(), String> {
    let files = loose_files(e);
    let mut done = 0;
    for (i, (url, name)) in files.iter().enumerate() {
        let dest = dir.join(name);
        if dest.is_file() {
            done += 1;
            continue;
        }
        let label = format!("{} [{}/{}] {name}", e.id, i + 1, files.len());
        match fetch::download(url, &dest, &label) {
            Ok(()) => done += 1,
            Err(err) => eprintln!("  warning: {err}"),
        }
    }
    println!("  {done}/{} index files present", files.len());
    if done < files.len() {
        return Err(format!(
            "{}: {} of {} files could not be downloaded; run the install again to fetch them",
            e.id,
            files.len() - done,
            files.len()
        ));
    }
    Ok(())
}

/// Download a `.zip` or `.deb` and unpack its catalogue files into `dir`.
///
/// Files are extracted into a staging directory first and moved into place only
/// once the whole archive has unpacked. An install interrupted part way therefore
/// leaves nothing that looks installed — `is_installed` probes a single file, so a
/// half-extracted database would otherwise be reported as present, and skipped by
/// the next install.
fn install_archive(dir: &Path, e: &Entry, keep: bool) -> Result<(), String> {
    let ext = if e.archive == Archive::Zip {
        "zip"
    } else {
        "deb"
    };
    let tmp = dir.join(format!(".{}-download.{ext}", e.id));
    if !tmp.is_file() {
        fetch::download(e.url, &tmp, e.id)?;
    }

    let staging = dir.join(format!(".{}-staging", e.id));
    let _ = fs::remove_dir_all(&staging);
    fs::create_dir_all(&staging).map_err(|err| format!("{}: {err}", staging.display()))?;

    eprintln!("  extracting ...");
    let wanted = |name: &str| e.files.owns(name);
    let extracted = if e.archive == Archive::Zip {
        fetch::extract_zip(&tmp, &staging, &wanted)
    } else {
        fetch::extract_deb(&tmp, &staging, &wanted)
    };
    let moved = extracted.and_then(|n| {
        for entry in
            fs::read_dir(&staging).map_err(|err| format!("{}: {err}", staging.display()))?
        {
            let from = entry.map_err(|err| err.to_string())?.path();
            let Some(name) = from.file_name() else {
                continue;
            };
            let to = dir.join(name);
            fs::rename(&from, &to).map_err(|err| format!("{}: {err}", to.display()))?;
        }
        Ok(n)
    });
    let _ = fs::remove_dir_all(&staging);
    let n = moved.map_err(|err| {
        format!(
            "{err}\n  (the archive is kept at {}; delete it to download afresh)",
            tmp.display()
        )
    })?;

    if keep {
        println!("  kept archive at {}", tmp.display());
    } else {
        let _ = fs::remove_file(&tmp);
    }
    println!("  extracted {n} files");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn human_sizes_read_sensibly() {
        assert_eq!(human(500), "500 B");
        assert_eq!(human(1_500_000), "1.5 MB");
        assert_eq!(human(901_300_000), "901.3 MB");
        assert_eq!(human(1_213_400_000), "1.2 GB");
    }
}
