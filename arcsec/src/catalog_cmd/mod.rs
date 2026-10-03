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
//!
//! The work itself is the `arcsec-catalogue` library's; this module parses the
//! command line, asks the questions, and prints.

mod download;
pub mod index_cmd;
mod prompt;

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};

use arcsec_catalogue::index::{
    self as index, Existing, IndexAction, IndexOptions, Machine, Plan, Rebuild, check_disk,
    concerns, rebuild_reason,
};
use arcsec_catalogue::registry::{
    Entry, Purpose, REGISTRY, files_of, find, installed_size, is_complete, is_installed,
};
use arcsec_catalogue::{
    CatalogueCheck, Error, IndexHealth, InstallEvent, InstallOptions, download_disk, duration,
    human_bytes as human, recommend,
};
use prompt::{Prompter, Terminal, install_questions};

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
        Some(("install", sm)) => cmd_install(
            &dir,
            &names(sm),
            sm.get_flag("yes"),
            sm.get_flag("keep"),
            &IndexOptions {
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

/// Where catalogues are kept: `$ARCSEC_CATALOG_DIR`, else the platform's data
/// directory. The rules live in [`arcsec_core::auto::default_catalog_dir`], which
/// the solver (and the C library) use too.
pub fn default_dir() -> PathBuf {
    arcsec_catalogue::default_dir()
}

/// The message for a library error, with what the command line adds: the commands
/// that fix it.
fn message(e: &Error) -> String {
    match e {
        Error::Unpack { archive, source } => format!(
            "{source}\n  (the archive is kept at {}; delete it to download afresh)",
            archive.display()
        ),
        Error::NoDatabase { dir } => format!(
            "no star database in {}; install one (`arcsec catalog install d50`) or pass --db and --name",
            dir.display()
        ),
        Error::NotWritable { path, source } => format!(
            "cannot write the index to {}: {source}\n  If this is a shared or system folder (ASTAP's under Program Files, say), \
             run this once as administrator, or build elsewhere with -o and pass that file to the solver with -i.",
            path.display()
        ),
        e => e.to_string(),
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
    let pick = |purpose: Purpose| recommend(purpose, fov_deg);

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
    // The blind index is built from the solving database on install, so there is
    // nothing to choose; say what it will cost. The default index stops at 0.3°
    // fields, so a narrower field needs the 0.06° tier, asked for at install.
    let narrow = fov_deg < NARROW_INDEX_FOV;
    let index_flag = if narrow { " --index-min-fov 0.15" } else { "" };
    if let Some(e) = pick(Purpose::Solving)
        && let Some(plan) = Plan::for_databases(&[e.id], narrow.then_some(0.15), None)
    {
        let est = index::Estimate::of(&plan, arcsec_core::max_threads());
        println!(
            "  blind index  built from {} on install{index_flag}, no download: ~{}, {}",
            e.id,
            human(est.bytes),
            duration(est.secs)
        );
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
        println!("\n  arcsec catalog install {}{index_flag}", ids.join(" "));
    }
}

/// Fields narrower than this (degrees) need the blind index's deepest tier, which
/// the default index leaves out; `recommend` adds `--index-min-fov 0.15` for them.
const NARROW_INDEX_FOV: f64 = 0.3;

/// The blind-index lines of `catalog list`: each index with its source and any
/// staleness, or what building one would take.
fn print_index_status(dir: &Path) {
    let files = index::index_files(dir);
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
    let sources = index::installed_sources(dir);
    if files.is_empty() {
        match index_cmd::build_suggestion(dir, "") {
            Some(s) => println!("  none: {s}"),
            None => println!("  none (installing a solving database builds one)"),
        }
    } else if let Some(plan) = Plan::for_databases(&sources, None, None)
        && let Some(why) = rebuild_reason(Existing::preferred(dir).as_ref(), &plan, dir)
        // A changed database is already flagged STALE on the index's own line.
        && !matches!(why, Rebuild::Changed(_))
    {
        println!("  rebuild suggested: {why}; run `arcsec catalog index build`");
    }
    for p in index::stale_parts(dir) {
        println!(
            "  {} is left from an interrupted build; delete it, or the next build will",
            p.display()
        );
    }
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
        index::indexes_removed_with(dir, &wanted)
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
        let removed = arcsec_catalogue::remove(dir, e).map_err(|err| message(&err))?;
        println!(
            "removed {}: {} files, {} freed",
            e.id,
            removed.files,
            human(removed.bytes)
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
    let v = arcsec_catalogue::verify(dir);
    for c in &v.catalogues {
        print_catalogue_check(c);
    }
    for ix in &v.indexes {
        let name = ix
            .path
            .file_name()
            .map_or_else(String::new, |n| n.to_string_lossy().into_owned());
        match &ix.health {
            IndexHealth::Stale { source } => {
                println!(
                    "  {name:<12} STALE: {} has changed since it was built - rebuild it with `arcsec catalog index build`",
                    source.to_uppercase()
                );
            }
            IndexHealth::Sound {
                patterns,
                bytes,
                existing,
                ..
            } => {
                println!("  {name:<12} ok  ({patterns} patterns, {})", human(*bytes));
                if let Some(note) = existing
                    .as_ref()
                    .and_then(|ex| index_cmd::freshness_note(ex, dir))
                {
                    println!("      {note}");
                }
            }
            IndexHealth::Broken(e) => {
                println!(
                    "  {name:<12} PROBLEM: {e} - rebuild it with `arcsec catalog index build`"
                );
            }
            _ => {}
        }
    }
    let sources = index::installed_sources(dir);
    if let Some(plan) = Plan::for_databases(&sources, None, None) {
        let existing = Existing::preferred(dir);
        match (existing.as_ref(), index_cmd::build_suggestion(dir, "")) {
            (None, Some(s)) => println!("  blind index none: {s}"),
            (Some(_), _) => {
                if let Some(why) = rebuild_reason(existing.as_ref(), &plan, dir)
                    && !matches!(why, Rebuild::Changed(_))
                {
                    println!("  blind index: {why}; `arcsec catalog index build` rebuilds it");
                }
            }
            (None, None) => {}
        }
    }
    for p in index::stale_parts(dir) {
        println!(
            "  note: {} is left from an interrupted build; delete it, or the next build will",
            p.display()
        );
    }
    if v.is_empty() {
        println!("No catalogues installed in {}", dir.display());
        return Ok(());
    }
    let problems = v.problems();
    if problems > 0 {
        return Err(format!(
            "{problems} problem(s) found; re-run `arcsec catalog install <name>`, or `arcsec catalog index build` for an index"
        ));
    }
    Ok(())
}

/// One catalogue's line (and problems) in `catalog verify`.
fn print_catalogue_check(c: &CatalogueCheck) {
    if c.problems.is_empty() {
        println!(
            "  {:<12} ok  ({} files, {})",
            c.entry.id,
            c.files,
            human(c.bytes)
        );
        return;
    }
    println!("  {:<12} PROBLEMS:", c.entry.id);
    for b in &c.problems {
        println!("      {b}");
    }
}

fn cmd_install(
    dir: &Path,
    ids: &[String],
    assume_yes: bool,
    keep: bool,
    index: &IndexOptions,
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
    let dbs = index::sources_after_install(dir, &wanted);
    let machine = Machine::probe(dir, 0);
    let action = index::index_action(dir, &dbs, index, machine.threads);
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
    check_disk(dir, need, machine.free_disk).map_err(|e| message(&e))?;
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
            install_one(dir, e, keep)?;
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

/// Download and unpack one catalogue, reporting as it goes.
fn install_one(dir: &Path, e: &'static Entry, keep: bool) -> Result<(), String> {
    let mut progress = download::Progress::default();
    let result = arcsec_catalogue::install(
        dir,
        e,
        &InstallOptions { keep_archive: keep },
        &mut |event| match event {
            InstallEvent::Download { label, event } => progress.event(label, event),
            InstallEvent::FileFailed(err) => eprintln!("  warning: {}", message(err)),
            InstallEvent::FilesPresent { present, total } => {
                println!("  {present}/{total} index files present");
            }
            InstallEvent::Extracting => eprintln!("  extracting ..."),
            InstallEvent::Decompressing(name) => eprintln!("  decompressing {name} ..."),
            InstallEvent::ArchiveKept(p) => println!("  kept archive at {}", p.display()),
            InstallEvent::Extracted(n) => println!("  extracted {n} files"),
            _ => {}
        },
    );
    match result {
        Ok(got) => {
            println!("  {} installed ({})", e.id, human(got.bytes));
            Ok(())
        }
        Err(err) => Err(message(&err)),
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    use arcsec_catalogue::test_support::{TempDir, write_001_db};
    use prompt::tests::Scripted;

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
        let none = IndexOptions::default();

        // No terminal, no --yes: cancelled, as installs always were.
        let mut p = Scripted::new(&[], false);
        cmd_install(d, &ids(&["w08"]), false, false, &none, &mut p).unwrap();
        assert_eq!(p.asked, ["Build the blind index now?"]);
        assert!(!ix.exists());

        // --no-index: nothing to do, nothing asked.
        let mut p = Scripted::new(&[], true);
        let skip = IndexOptions {
            skip: true,
            ..IndexOptions::default()
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
            &IndexOptions::default(),
            &mut p,
        )
        .unwrap();
        assert!(d.join("w08.arcsecix").is_file());

        // G05 appears (here as a tiny all-sky file: the layout does not matter).
        write_001_db(d, "g05", 2500, 6);
        let opts = IndexOptions::default();
        let a = index::index_action(d, &["g05", "w08"], &opts, 2).unwrap();
        assert_eq!(a.plan.source, "g05");
        assert!(matches!(a.reason, Rebuild::Deeper { .. }), "{}", a.reason);
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
            &IndexOptions::default(),
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
}
