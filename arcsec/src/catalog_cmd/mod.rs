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
pub mod plan;
mod prompt;
mod registry;
mod sys;
#[cfg(test)]
mod testutil;

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};

pub use registry::{ASTAP_EXTS, Files, REGISTRY, is_installed};
use registry::{
    Archive, Entry, Purpose, astap_file_count, expected_file_count, files_of, find, installed_size,
    is_complete, loose_files,
};

use plan::{
    Estimate, Existing, Machine, Plan, SOURCES, check_disk, concerns, depth_rank,
    keep_existing_range, rebuild_reason,
};
use prompt::{Prompter, Terminal, install_questions};

use crate::{cli, image_io};

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
        Some(("install", sm)) => cmd_install(
            &dir,
            &names(sm),
            sm.get_flag("yes"),
            sm.get_flag("keep"),
            &IndexOpts {
                skip: sm.get_flag("no-index"),
                min_fov: sm.get_one::<f64>("index-min-fov").copied(),
                max_fov: sm.get_one::<f64>("index-max-fov").copied(),
            },
            &mut Terminal,
        ),
        Some(("remove", sm)) => cmd_remove(
            &dir,
            &names(sm),
            sm.get_flag("yes"),
            sm.get_flag("keep-index"),
            &mut Terminal,
        ),
        Some(("verify", _)) => cmd_verify(&dir),
        Some(("index", sm)) => match sm.subcommand() {
            Some(("build", b)) => index_cmd::cmd_build(
                &dir,
                b.get_one::<PathBuf>("db"),
                b.get_one::<String>("name"),
                b.get_one::<f64>("min-fov").copied(),
                b.get_one::<f64>("max-fov").copied(),
                b.get_one::<PathBuf>("out"),
                *b.get_one::<usize>("threads").unwrap_or(&0),
                b.get_flag("yes"),
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

/// Where catalogues are kept, in priority order:
///
/// 1. `$ARCSEC_CATALOG_DIR`, if set — for people who keep them on another disk.
/// 2. `$XDG_DATA_HOME/arcsec/catalogs` on Linux, or the platform equivalent:
///    `~/Library/Application Support/arcsec/catalogs` on macOS,
///    `%LOCALAPPDATA%\arcsec\catalogs` on Windows.
/// 3. `~/.arcsec/catalogs` if the home directory cannot be resolved any other way.
///
/// The point is that a user who runs `arcsec catalog install d50` never has to know
/// this path, and the solver looks here without being told.
pub fn default_dir() -> PathBuf {
    default_dir_from(|k| std::env::var(k).ok())
}

/// [`default_dir`] with the environment supplied by `var`, so it can be tested
/// without mutating the real process environment.
fn default_dir_from(var: impl Fn(&str) -> Option<String>) -> PathBuf {
    let get = |k: &str| var(k).filter(|v| !v.is_empty()).map(PathBuf::from);

    if let Some(p) = get("ARCSEC_CATALOG_DIR") {
        return p;
    }

    let platform = if cfg!(target_os = "windows") {
        get("LOCALAPPDATA").map(|p| p.join("arcsec").join("catalogs"))
    } else if cfg!(target_os = "macos") {
        get("HOME").map(|h| {
            h.join("Library")
                .join("Application Support")
                .join("arcsec")
                .join("catalogs")
        })
    } else {
        get("XDG_DATA_HOME")
            .map(|p| p.join("arcsec").join("catalogs"))
            .or_else(|| {
                get("HOME").map(|h| {
                    h.join(".local")
                        .join("share")
                        .join("arcsec")
                        .join("catalogs")
                })
            })
    };

    platform.unwrap_or_else(|| {
        get("HOME").or_else(|| get("USERPROFILE")).map_or_else(
            || PathBuf::from("catalogs"),
            |h| h.join(".arcsec").join("catalogs"),
        )
    })
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
    println!("\nBlind index (arcsec's own, built from a star database; nothing to download):");
    print_index_status(dir);
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

/// The blind-index lines of `catalog list`: each index with its source and any
/// staleness, or what building one would take.
fn print_index_status(dir: &Path) {
    let files = index_cmd::index_files(dir);
    let used = Existing::preferred(dir).map(|e| e.path);
    for p in &files {
        let name = p
            .file_name()
            .map_or_else(String::new, |n| n.to_string_lossy().into_owned());
        let Some(ex) = Existing::open(p) else {
            let why = arcsec_core::index::BlindIndex::open(p)
                .err()
                .map_or_else(String::new, |e| e.to_string());
            println!("  {name}: unusable ({why}); rebuild it with `arcsec catalog index build`");
            continue;
        };
        let size = fs::metadata(p).map_or(0, |m| m.len());
        let tag = if files.len() > 1 && used.as_ref() == Some(p) {
            "  [used by the solver]"
        } else {
            ""
        };
        println!(
            "  {name:<15} from {}, fields {:.2}°–{:.0}°, {}{tag}",
            ex.source.to_uppercase(),
            ex.coverage.0,
            ex.coverage.1,
            human(size)
        );
        if let Some(note) = index_cmd::freshness_note(&ex, dir) {
            println!("      {note}");
        }
    }
    let sources = index_cmd::installed_sources(dir);
    if files.is_empty() {
        match index_cmd::build_suggestion(dir, "") {
            Some(s) => println!("  none: {s}"),
            None => println!("  none (installing a solving database builds one)"),
        }
    } else if let Some(plan) = Plan::for_databases(&sources, None, None)
        && let Some(why) = rebuild_reason(Existing::preferred(dir).as_ref(), &plan, dir)
    {
        println!("  rebuild suggested: {why}; run `arcsec catalog index build`");
    }
    for p in index_cmd::stale_parts(dir) {
        println!(
            "  {} is left from an interrupted build; delete it, or the next build will",
            p.display()
        );
    }
}

/// Index files in `dir` built from database `db` (by their header).
fn indexes_built_from(dir: &Path, db: &str) -> Vec<PathBuf> {
    index_cmd::index_files(dir)
        .into_iter()
        .filter(|p| Existing::open(p).is_some_and(|e| e.source.eq_ignore_ascii_case(db)))
        .collect()
}

fn cmd_remove(
    dir: &Path,
    ids: &[String],
    assume_yes: bool,
    keep_index: bool,
    prompter: &mut dyn Prompter,
) -> Result<(), String> {
    let mut wanted = Vec::new();
    for id in ids {
        let e = find(id)
            .ok_or_else(|| format!("unknown catalogue '{id}' (try: arcsec catalog list)"))?;
        if !wanted.iter().any(|w: &&Entry| w.id == e.id) {
            wanted.push(e);
        }
    }

    // A blind index built from a database goes with it unless --keep-index: it
    // still solves on its own, but nothing would update or verify it any more, and
    // it is usually the bigger file.
    let indexes: Vec<PathBuf> = if keep_index {
        Vec::new()
    } else {
        wanted
            .iter()
            .filter(|e| e.purpose == Purpose::Solving)
            .flat_map(|e| indexes_built_from(dir, e.id))
            .collect()
    };

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
    for p in &indexes {
        println!(
            "  {:<23} {:>9}  blind index built from it (keep it with --keep-index)",
            p.file_name()
                .map_or_else(String::new, |n| n.to_string_lossy().into_owned()),
            human(fs::metadata(p).map_or(0, |m| m.len()))
        );
    }
    println!();
    if !assume_yes && !prompter.ask("Continue?") {
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
    for p in &indexes {
        let freed = fs::metadata(p).map_or(0, |m| m.len());
        fs::remove_file(p).map_err(|err| format!("{}: {err}", p.display()))?;
        println!("removed {}: {} freed", p.display(), human(freed));
    }
    if Existing::preferred(dir).is_none()
        && let Some(s) = index_cmd::build_suggestion(dir, "")
    {
        println!("\nNo blind index is left; {s}.");
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
    let index_files = index_cmd::index_files(dir);
    for p in &index_files {
        checked += 1;
        let name = p
            .file_name()
            .map_or_else(String::new, |n| n.to_string_lossy().into_owned());
        let res = arcsec_core::index::BlindIndex::open(p).and_then(|ix| {
            ix.validate()?;
            Ok(ix)
        });
        match res {
            Ok(ix) => {
                let ex = Existing::open(p);
                match ex.as_ref().map(|ex| plan::freshness(ex, dir)) {
                    Some(plan::Freshness::Changed) => {
                        problems += 1;
                        println!(
                            "  {name:<11} STALE: {} has changed since it was built - rebuild it with `arcsec catalog index build`",
                            ix.source().to_uppercase()
                        );
                    }
                    _ => {
                        println!(
                            "  {name:<11} ok  ({} patterns, {})",
                            ix.n_patterns(),
                            human(ix.file_size() as u64)
                        );
                        if let Some(note) = ex
                            .as_ref()
                            .and_then(|ex| index_cmd::freshness_note(ex, dir))
                        {
                            println!("      {note}");
                        }
                    }
                }
            }
            Err(e) => {
                problems += 1;
                println!(
                    "  {name:<11} PROBLEM: {e} - rebuild it with `arcsec catalog index build`"
                );
            }
        }
    }
    let sources = index_cmd::installed_sources(dir);
    if let Some(plan) = Plan::for_databases(&sources, None, None) {
        let existing = Existing::preferred(dir);
        match (existing.as_ref(), index_cmd::build_suggestion(dir, "")) {
            (None, Some(s)) => println!("  blind index none: {s}"),
            (Some(_), _) => {
                if let Some(why) = rebuild_reason(existing.as_ref(), &plan, dir)
                    && !why.contains("has changed")
                {
                    println!("  blind index: {why}; `arcsec catalog index build` rebuilds it");
                }
            }
            (None, None) => {}
        }
    }
    for p in index_cmd::stale_parts(dir) {
        println!(
            "  note: {} is left from an interrupted build; delete it, or the next build will",
            p.display()
        );
    }
    if checked == 0 {
        println!("No catalogues installed in {}", dir.display());
        return Ok(());
    }
    if problems > 0 {
        return Err(format!(
            "{problems} problem(s) found; re-run `arcsec catalog install <name>`, or `arcsec catalog index build` for an index"
        ));
    }
    Ok(())
}

/// `catalog install`'s blind-index options.
#[derive(Debug, Clone, Copy, Default)]
pub struct IndexOpts {
    /// `--no-index`: download only.
    pub skip: bool,
    /// `--index-min-fov`.
    pub min_fov: Option<f64>,
    /// `--index-max-fov`.
    pub max_fov: Option<f64>,
}

/// A blind index `install` will build once the downloads are in.
#[derive(Debug, Clone)]
struct IndexAction {
    plan: Plan,
    est: Estimate,
    /// Why: no index yet, a deeper database, a stale or too-narrow index.
    reason: String,
    /// Indexes in the directory the new one supersedes (removed after it is built).
    replaces: Vec<PathBuf>,
}

/// The index to build after installing so that `dbs` (every solving database that
/// will then be in `dir`) are served, or `None` if the one there already does.
fn index_action(dir: &Path, dbs: &[&str], opts: &IndexOpts, threads: usize) -> Option<IndexAction> {
    if opts.skip {
        return None;
    }
    let base = Plan::for_databases(dbs, opts.min_fov, opts.max_fov)?;
    let existing = Existing::preferred(dir);
    let reason = rebuild_reason(existing.as_ref(), &base, dir)?;
    let plan = keep_existing_range(
        base,
        existing.as_ref(),
        opts.min_fov.is_some(),
        opts.max_fov.is_some(),
    );
    // Our own indexes, by their default names: the new one replaces any built from
    // another database, since the solver uses only one.
    let replaces = SOURCES
        .iter()
        .map(|db| arcsec_core::index::default_index_path(dir, db))
        .filter(|p| p.is_file())
        .collect();
    let est = Estimate::of(&plan, threads);
    Some(IndexAction {
        plan,
        est,
        reason,
        replaces,
    })
}

/// Disk an install of `e` needs at its peak: an archive and its extracted files
/// side by side (assumed no smaller than the archive), or the loose files.
fn download_disk(e: &Entry) -> u64 {
    match e.archive {
        Archive::Loose => e.bytes,
        Archive::Zip | Archive::Deb => e.bytes.saturating_mul(2),
    }
}

fn cmd_install(
    dir: &Path,
    ids: &[String],
    assume_yes: bool,
    keep: bool,
    index: &IndexOpts,
    prompter: &mut dyn Prompter,
) -> Result<(), String> {
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

    // The solving databases there will be once this is done, deepest first.
    let mut dbs = index_cmd::installed_sources(dir);
    for e in &wanted {
        if let Some(s) = SOURCES.iter().find(|s| **s == e.id)
            && !dbs.contains(s)
        {
            dbs.push(s);
        }
    }
    dbs.sort_by_key(|d| depth_rank(d));
    let machine = Machine::probe(dir, 0);
    let action = index_action(dir, &dbs, index, machine.threads);
    if wanted.is_empty() && action.is_none() {
        return Ok(());
    }

    println!("Installing into {}\n", dir.display());
    let total: u64 = wanted.iter().map(|e| e.bytes).sum();
    if !wanted.is_empty() {
        for e in &wanted {
            println!("  {:<11} {:>9}  {}", e.id, human(e.bytes), e.desc);
        }
        println!("\nTotal download: {}", human(total));
    }
    if let Some(a) = &action {
        println!(
            "{}uild a blind index from {} ({}):",
            if wanted.is_empty() { "B" } else { "Then b" },
            a.plan.label(),
            a.reason
        );
        println!("  {}", a.est.summary());
        for r in &a.replaces {
            println!(
                "  replacing {}",
                r.file_name()
                    .map_or_else(String::new, |n| n.to_string_lossy().into_owned())
            );
        }
        let narrow = if a.plan.source == "d80" && a.plan.min_fov > 0.15 {
            "; --index-min-fov 0.15 for D80's narrowest fields"
        } else {
            ""
        };
        println!("  (--no-index to skip it{narrow})");
    }

    let need = wanted.iter().map(|e| download_disk(e)).sum::<u64>()
        + action.as_ref().map_or(0, |a| a.est.bytes);
    check_disk(dir, need, machine.free_disk)?;
    let worries = action
        .as_ref()
        .map_or_else(Vec::new, |a| concerns(&a.est, &machine));
    index_cmd::print_concerns(&worries);
    println!();

    let answer = install_questions(
        !wanted.is_empty(),
        action.is_some(),
        !worries.is_empty(),
        assume_yes,
        prompter,
    );
    if !answer.download && !answer.build_index {
        println!("Cancelled.");
        return Ok(());
    }

    if answer.download {
        fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
        for e in &wanted {
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
    }

    if let Some(a) = &action {
        if answer.build_index {
            build_after_install(dir, a, machine.threads)?;
        } else {
            println!(
                "\nSkipped the blind index. Build it later with `arcsec catalog index build`."
            );
        }
    }
    println!("\nDone. The solver uses {} by default.", dir.display());
    Ok(())
}

/// Build the index `install` planned, then remove the ones it supersedes.
fn build_after_install(dir: &Path, a: &IndexAction, threads: usize) -> Result<(), String> {
    let out = arcsec_core::index::default_index_path(dir, &a.plan.source);
    println!("\nBlind index from {}", a.plan.label());
    index_cmd::run_build(&a.plan, dir, &out, threads, &a.est).map_err(|e| {
        format!(
            "{e}\n  The catalogues are installed; build the index later with `arcsec catalog index build`."
        )
    })?;
    for r in a.replaces.iter().filter(|r| **r != out) {
        match fs::remove_file(r) {
            Ok(()) => println!("  removed {} (superseded)", r.display()),
            Err(e) => eprintln!("  warning: could not remove {}: {e}", r.display()),
        }
    }
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
    fn default_dir_is_namespaced() {
        // Deliberately not asserting the exact path: it is platform dependent.
        let d = default_dir();
        assert!(
            d.to_string_lossy().contains("arcsec"),
            "catalogue dir should be namespaced: {}",
            d.display()
        );
    }

    #[test]
    fn env_override_wins() {
        let env = |k: &str| match k {
            "ARCSEC_CATALOG_DIR" => Some("/data/catalogs".to_string()),
            "HOME" => Some("/home/u".to_string()),
            _ => None,
        };
        assert_eq!(default_dir_from(env), PathBuf::from("/data/catalogs"));
    }

    #[test]
    fn empty_variables_count_as_unset() {
        let env = |k: &str| match k {
            "ARCSEC_CATALOG_DIR" | "XDG_DATA_HOME" | "LOCALAPPDATA" => Some(String::new()),
            "HOME" => Some("/home/u".to_string()),
            _ => None,
        };
        let d = default_dir_from(env);
        assert!(d.starts_with("/home/u"), "got {}", d.display());
        assert!(d.ends_with("catalogs"));
    }

    #[test]
    fn no_home_at_all_still_yields_a_path() {
        assert_eq!(default_dir_from(|_| None), PathBuf::from("catalogs"));
    }

    #[test]
    fn human_sizes_read_sensibly() {
        assert_eq!(human(500), "500 B");
        assert_eq!(human(1_500_000), "1.5 MB");
        assert_eq!(human(901_300_000), "901.3 MB");
        assert_eq!(human(1_213_400_000), "1.2 GB");
    }

    use prompt::tests::Scripted;
    use testutil::{TempDir, write_001_db};

    fn ids(v: &[&str]) -> Vec<String> {
        v.iter().map(ToString::to_string).collect()
    }

    /// `install` of a database that is already there (so nothing is downloaded)
    /// builds the missing index: asked once, honouring --yes, --no-index and end
    /// of input; then a second install has nothing to do.
    #[test]
    fn install_builds_a_missing_index_for_an_existing_database() {
        let dir = TempDir::new("inst_ix");
        let d = dir.path();
        write_001_db(d, "w08", 2000, 3);
        let ix = d.join("w08.arcsecix");
        let none = IndexOpts::default();

        // No terminal, no --yes: cancelled, as installs always were.
        let mut p = Scripted::new(&[], false);
        cmd_install(d, &ids(&["w08"]), false, false, &none, &mut p).unwrap();
        assert_eq!(p.asked, ["Build the blind index now?"]);
        assert!(!ix.exists());

        // --no-index: nothing to do, nothing asked.
        let mut p = Scripted::new(&[], true);
        let skip = IndexOpts {
            skip: true,
            ..IndexOpts::default()
        };
        cmd_install(d, &ids(&["w08"]), false, false, &skip, &mut p).unwrap();
        assert!(p.asked.is_empty() && !ix.exists());

        // Answered yes: built.
        let mut p = Scripted::new(&[true], true);
        cmd_install(d, &ids(&["w08"]), false, false, &none, &mut p).unwrap();
        assert!(ix.is_file());
        let built = arcsec_core::index::BlindIndex::open(&ix).unwrap();
        assert_eq!(built.source(), "w08");
        assert!(built.source_stamp().is_recorded());

        // Up to date: no question at all.
        let mut p = Scripted::new(&[], true);
        cmd_install(d, &ids(&["w08"]), false, false, &none, &mut p).unwrap();
        assert!(p.asked.is_empty());
        cmd_verify(d).unwrap();

        // The database changes: verify reports the index stale, and the next
        // install rebuilds it without asking under --yes.
        write_001_db(d, "w08", 2001, 4);
        assert!(cmd_verify(d).unwrap_err().contains("problem"));
        let mut p = Scripted::new(&[], false);
        cmd_install(d, &ids(&["w08"]), true, false, &none, &mut p).unwrap();
        assert!(p.asked.is_empty());
        cmd_verify(d).unwrap();
    }

    #[test]
    fn a_deeper_database_supersedes_the_old_index() {
        let dir = TempDir::new("inst_deeper");
        let d = dir.path();
        write_001_db(d, "w08", 2000, 5);
        let mut p = Scripted::new(&[], false);
        cmd_install(
            d,
            &ids(&["w08"]),
            true,
            false,
            &IndexOpts::default(),
            &mut p,
        )
        .unwrap();
        assert!(d.join("w08.arcsecix").is_file());

        // G05 appears (here as a tiny all-sky file: the layout does not matter).
        write_001_db(d, "g05", 2500, 6);
        let opts = IndexOpts::default();
        let a = index_action(d, &["g05", "w08"], &opts, 2).unwrap();
        assert_eq!(a.plan.source, "g05");
        assert!(a.reason.contains("deeper"), "{}", a.reason);
        // Covers both: G05's 3° up to W08's 80°.
        assert_eq!((a.plan.min_fov, a.plan.max_fov), (3.0, 80.0));
        cmd_install(d, &ids(&["g05"]), true, false, &opts, &mut p).unwrap();
        assert!(d.join("g05.arcsecix").is_file());
        assert!(!d.join("w08.arcsecix").exists(), "superseded index removed");
        assert_eq!(Existing::preferred(d).unwrap().source, "g05");
    }

    #[test]
    fn remove_takes_the_index_built_from_the_database_unless_told_not_to() {
        let dir = TempDir::new("rm_ix");
        let d = dir.path();
        write_001_db(d, "w08", 2000, 9);
        let mut p = Scripted::new(&[], false);
        cmd_install(
            d,
            &ids(&["w08"]),
            true,
            false,
            &IndexOpts::default(),
            &mut p,
        )
        .unwrap();
        let ix = d.join("w08.arcsecix");
        assert!(ix.is_file());

        // End of input: nothing removed.
        let mut p = Scripted::new(&[], false);
        cmd_remove(d, &ids(&["w08"]), false, false, &mut p).unwrap();
        assert!(ix.is_file() && d.join("w08_0101.001").is_file());

        // --keep-index keeps it.
        cmd_remove(d, &ids(&["w08"]), true, true, &mut p).unwrap();
        assert!(ix.is_file() && !d.join("w08_0101.001").exists());

        // Otherwise it goes with the database.
        write_001_db(d, "w08", 2000, 9);
        cmd_remove(d, &ids(&["w08"]), true, false, &mut p).unwrap();
        assert!(!ix.exists());
    }

    #[test]
    fn install_index_options_shape_the_plan() {
        let dir = TempDir::new("inst_opts");
        let d = dir.path();
        let opts = IndexOpts {
            skip: false,
            min_fov: Some(0.15),
            max_fov: None,
        };
        let a = index_action(d, &["d80"], &opts, 8).unwrap();
        assert_eq!(a.plan.tiers.last().unwrap().radius_deg, 0.06);
        assert!(a.est.bytes > 600_000_000, "{:?}", a.est);
        assert!(
            index_action(d, &[], &opts, 8).is_none(),
            "no solving database"
        );
        assert!(index_action(d, &["d80"], &IndexOpts { skip: true, ..opts }, 8).is_none());
    }

    #[test]
    fn disk_needed_for_downloads_counts_archive_and_contents() {
        let d80 = find("d80").unwrap();
        assert_eq!(download_disk(d80), 2 * d80.bytes);
        let anet = find("anet-4100").unwrap();
        assert_eq!(download_disk(anet), anet.bytes);
    }
}
