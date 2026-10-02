//! `arcsec catalog index build|info` — arcsec's own blind index, built from a star
//! database the user already has. See `docs/offline-index.md`.
//!
//! `catalog install` builds one automatically after downloading a solving database
//! (see [`super::cmd_install`]); this module holds the build itself, its progress
//! report, and how an installed index is described.

use std::io::IsTerminal as _;
use std::path::{Path, PathBuf};
use std::time::Instant;

use arcsec_core::catalog::catalog_present;
use arcsec_core::index::{
    BlindIndex, BuildParams, BuildProgress, DEFAULT_TIERS, TierSpec, build_index,
    default_index_path, tier_fov_range,
};

use super::human;
use super::plan::{
    Estimate, Existing, Freshness, Machine, Plan, SOURCES, check_disk, concerns, duration,
    freshness,
};
use super::prompt::{Terminal, build_question};

/// Solving databases installed in `dir`, deepest first.
pub fn installed_sources(dir: &Path) -> Vec<&'static str> {
    SOURCES
        .into_iter()
        .filter(|n| catalog_present(dir, n))
        .collect()
}

/// The tiers needed for fields (short side) from `min_fov` to `max_fov` degrees:
/// every tier whose field range overlaps that span, except that of the tiers
/// already reaching below `min_fov` only the widest is kept, and of those reaching
/// above `max_fov` only the narrowest — so both ends are covered and no deeper or
/// wider tier than necessary is built. A span that single tiers cover whole gets
/// just those. Widest first.
#[must_use]
pub fn tiers_for(min_fov: f64, max_fov: f64) -> Vec<TierSpec> {
    let range = |t: &TierSpec| tier_fov_range(t.radius_deg);
    let overlapping: Vec<TierSpec> = DEFAULT_TIERS
        .into_iter()
        .filter(|t| range(t).0 <= max_fov && range(t).1 >= min_fov)
        .collect();
    // A span narrow enough for single tiers to cover it whole: those tiers.
    let spanning: Vec<TierSpec> = overlapping
        .iter()
        .copied()
        .filter(|t| range(t).0 <= min_fov && range(t).1 >= max_fov)
        .collect();
    if !spanning.is_empty() {
        return spanning;
    }
    let widest_below = overlapping
        .iter()
        .filter(|t| range(t).0 <= min_fov)
        .map(|t| t.radius_deg)
        .fold(f64::NEG_INFINITY, f64::max);
    let narrowest_above = overlapping
        .iter()
        .filter(|t| range(t).1 >= max_fov)
        .map(|t| t.radius_deg)
        .fold(f64::INFINITY, f64::min);
    overlapping
        .into_iter()
        .filter(|t| {
            let (lo, hi) = range(t);
            #[allow(clippy::float_cmp)] // comparing a value with itself, copied
            let keep_lo = lo > min_fov || t.radius_deg == widest_below;
            #[allow(clippy::float_cmp)]
            let keep_hi = hi < max_fov || t.radius_deg == narrowest_above;
            keep_lo && keep_hi
        })
        .collect()
}

/// The temporary file [`arcsec_core::index::BuiltIndex::write`] writes `out` through.
pub fn part_path(out: &Path) -> PathBuf {
    out.with_extension(format!("{}.part", arcsec_core::index::format::EXTENSION))
}

/// The plan for `catalog index build`: `--name` (or the deepest installed database)
/// as the source, and `--min-fov`/`--max-fov` or else the fields every installed
/// database covers.
///
/// # Errors
///
/// A message if there is no database, the named one is missing, or the range is
/// empty or served by no tier.
pub fn plan_for_build(
    db_path: &Path,
    name: Option<&str>,
    min_fov: Option<f64>,
    max_fov: Option<f64>,
) -> Result<Plan, String> {
    let installed = installed_sources(db_path);
    let plan = match name.map(str::to_ascii_lowercase) {
        Some(n) if installed.first().is_none_or(|f| *f != n) => {
            if !catalog_present(db_path, &n) {
                return Err(format!(
                    "database {n} not found in {}",
                    db_path.display()
                ));
            }
            let (lo, hi) = super::plan::default_fields(&n);
            Plan::new(&n, min_fov.unwrap_or(lo), max_fov.unwrap_or(hi))
        }
        _ => Plan::for_databases(&installed, min_fov, max_fov).ok_or_else(|| {
            format!(
                "no star database in {}; install one (`arcsec catalog install d50`) or pass --db and --name",
                db_path.display()
            )
        })?,
    };
    if !(plan.min_fov > 0.0 && plan.max_fov >= plan.min_fov) {
        return Err(format!(
            "bad field range {}°–{}°",
            plan.min_fov, plan.max_fov
        ));
    }
    if plan.tiers.is_empty() {
        return Err(format!(
            "no tier fits fields {}°–{}° from {}",
            plan.min_fov,
            plan.max_fov,
            plan.source.to_uppercase()
        ));
    }
    Ok(plan)
}

/// Print a plan's tiers, indented.
pub fn print_tiers(plan: &Plan) {
    for t in &plan.tiers {
        let (lo, hi) = tier_fov_range(t.radius_deg);
        println!(
            "    disc {:>5}°  mag ≤ {:>4}  groups of {}  (fields {lo:.2}°–{hi:.1}°)",
            t.radius_deg, t.mag_cap, t.members
        );
    }
}

/// Print the notice for a build with concerns.
pub fn print_concerns(concerns: &[String]) {
    if concerns.is_empty() {
        return;
    }
    println!("\nNote: this blind index is a large job:");
    for c in concerns {
        println!("  - {c}");
    }
}

/// `arcsec catalog index build`.
#[allow(clippy::too_many_arguments)]
pub fn cmd_build(
    cat_dir: &Path,
    db: Option<&PathBuf>,
    name: Option<&String>,
    min_fov: Option<f64>,
    max_fov: Option<f64>,
    out: Option<&PathBuf>,
    threads: usize,
    assume_yes: bool,
) -> Result<(), String> {
    let db_path = db.cloned().unwrap_or_else(|| cat_dir.to_path_buf());
    let plan = plan_for_build(&db_path, name.map(String::as_str), min_fov, max_fov)?;
    let out = out
        .cloned()
        .unwrap_or_else(|| default_index_path(cat_dir, &plan.source));
    let out_dir = out
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(cat_dir);
    let machine = Machine::probe(out_dir, threads);
    let est = Estimate::of(&plan, machine.threads);

    println!(
        "Building a blind index from {} in {}",
        plan.source.to_uppercase(),
        db_path.display()
    );
    println!(
        "  fields {}°–{}° (short side), {} tiers:",
        plan.min_fov,
        plan.max_fov,
        plan.tiers.len()
    );
    print_tiers(&plan);
    println!("  output {}", out.display());
    println!("  estimate: {}", est.summary());

    check_disk(out_dir, est.bytes, machine.free_disk)?;
    let concerns = concerns(&est, &machine);
    print_concerns(&concerns);
    if !build_question(!concerns.is_empty(), assume_yes, &mut Terminal) {
        println!("Cancelled.");
        return Ok(());
    }
    run_build(&plan, &db_path, &out, machine.threads, &est).map(|_| ())
}

/// Build `plan` from the database in `db_path` and write it to `out`, reporting
/// progress on stderr. Returns the file size.
///
/// Nothing is written until the whole index is built in memory, and then it goes to
/// a `.part` file renamed into place, so an interrupted build never leaves a file
/// that looks like an index. A `.part` left by an earlier interrupted write is
/// removed first, and on Unix an interrupt during the write removes its own.
///
/// # Errors
///
/// A message if the database cannot be read or the file cannot be written.
pub fn run_build(
    plan: &Plan,
    db_path: &Path,
    out: &Path,
    threads: usize,
    est: &Estimate,
) -> Result<u64, String> {
    let part = part_path(out);
    if part.exists() {
        let _ = std::fs::remove_file(&part);
    }
    if let Some(parent) = out.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
    }
    // Find out now, not after the build, if the file cannot be written: an ASTAP
    // folder under Program Files or /opt is read-only to a normal user.
    std::fs::File::create(&part)
        .and_then(|_| std::fs::remove_file(&part))
        .map_err(|e| {
            format!(
                "cannot write the index to {}: {e}\n  If this is a shared or system folder (ASTAP's under Program Files, say), \
                 run this once as administrator, or build elsewhere with -o and pass that file to the solver with -i.",
                out.display()
            )
        })?;

    let t0 = Instant::now();
    let mut progress = Progress::new(est, t0);
    let built = build_index(
        &BuildParams {
            db_path: db_path.to_path_buf(),
            db_name: plan.source.clone(),
            tiers: plan.tiers.clone(),
            threads,
        },
        |p| progress.event(p),
    )
    .map_err(|e| format!("build failed: {e}"))?;

    {
        let _guard = interrupt::RemoveOnInterrupt::new(&part);
        built
            .write(out)
            .map_err(|e| format!("writing {}: {e}", out.display()))?;
    }
    let size = std::fs::metadata(out).map_or(0, |m| m.len());
    println!(
        "Wrote {} ({}, {} patterns, {} stars) in {:.1} s (estimated {}, {}).",
        out.display(),
        human(size),
        built.keys.len(),
        built.stars.len(),
        t0.elapsed().as_secs_f64(),
        human(est.bytes),
        duration(est.secs)
    );
    Ok(size)
}

/// Build progress on stderr: a line per tier, with the time left estimated from
/// the cost model's share of each tier, corrected by how fast the build is
/// actually going. On a terminal the current tier's line is redrawn in place.
struct Progress {
    shares: Vec<f64>,
    est_secs: f64,
    t0: Instant,
    t_tier: Instant,
    tier: usize,
    of: usize,
    label: String,
    done_share: f64,
    tty: bool,
    last_draw: Option<Instant>,
}

impl Progress {
    fn new(est: &Estimate, t0: Instant) -> Self {
        Self {
            shares: est.tier_share.clone(),
            est_secs: est.secs,
            t0,
            t_tier: t0,
            tier: 0,
            of: 0,
            label: String::new(),
            done_share: 0.0,
            tty: std::io::stderr().is_terminal(),
            last_draw: None,
        }
    }

    /// `about 40 s left`, given the fraction of the work done.
    fn left(&self, frac: f64) -> String {
        let elapsed = self.t0.elapsed().as_secs_f64();
        let secs = if frac >= 0.1 {
            elapsed * (1.0 - frac) / frac
        } else {
            (self.est_secs - elapsed).max(self.est_secs * (1.0 - frac) * 0.5)
        };
        if secs < 10.0 {
            "a few seconds left".to_string()
        } else {
            format!("about {} left", duration(secs).trim_start_matches('~'))
        }
    }

    fn event(&mut self, p: &BuildProgress) {
        match p {
            BuildProgress::Tier { index, of, spec } => {
                self.tier = *index;
                self.of = *of;
                self.t_tier = Instant::now();
                self.label = format!("  tier {}/{of} (disc {}°)", index + 1, spec.radius_deg);
                eprint!("{}", self.label);
                if !self.tty {
                    eprint!(" ");
                }
            }
            BuildProgress::Strip { done, of, .. } => {
                if self.tty {
                    // At most four redraws a second.
                    if self
                        .last_draw
                        .is_some_and(|t| t.elapsed().as_millis() < 250)
                    {
                        return;
                    }
                    self.last_draw = Some(Instant::now());
                    let share = self.shares.get(self.tier).copied().unwrap_or(0.0);
                    let frac =
                        (self.done_share + share * (*done as f64 / (*of).max(1) as f64)).min(1.0);
                    eprint!(
                        "\r\x1b[K{}  {done}/{of} strips, {:.0} % done, {}",
                        self.label,
                        frac * 100.0,
                        self.left(frac)
                    );
                } else if *of > 1 && done % 6 == 0 {
                    eprint!(".");
                }
            }
            BuildProgress::TierDone { info, stars_read } => {
                self.done_share += self.shares.get(self.tier).copied().unwrap_or(0.0);
                let left = if self.tier + 1 < self.of {
                    format!("; {}", self.left(self.done_share.min(1.0)))
                } else {
                    String::new()
                };
                if self.tty {
                    eprint!("\r\x1b[K{}", self.label);
                }
                eprintln!(
                    " {} anchors, {} patterns, {} ({stars_read} stars read, {:.1} s){left}",
                    info.n_anchors,
                    info.n_patterns,
                    human(info.n_patterns * 24),
                    self.t_tier.elapsed().as_secs_f64()
                );
            }
        }
    }
}

/// Removing a half-written index when the process is interrupted.
mod interrupt {
    use std::path::Path;

    /// While alive, SIGINT/SIGTERM/SIGHUP remove the file before the process dies
    /// of the signal. A signal the process was started ignoring (`nohup`, a
    /// background job) stays ignored. A no-op off Unix, where the next build
    /// removes a stale `.part` instead.
    pub struct RemoveOnInterrupt {
        /// The dispositions replaced, restored on drop; empty if not armed.
        #[cfg(unix)]
        previous: Vec<(libc::c_int, libc::sighandler_t)>,
    }

    #[cfg(unix)]
    mod imp {
        use core::sync::atomic::{AtomicPtr, Ordering};

        /// The path to remove, as a C string from `CString::into_raw`; null when
        /// disarmed.
        pub static PATH: AtomicPtr<libc::c_char> = AtomicPtr::new(core::ptr::null_mut());
        pub const SIGNALS: [libc::c_int; 3] = [libc::SIGINT, libc::SIGTERM, libc::SIGHUP];

        /// Makes only async-signal-safe calls: unlink, signal, raise.
        pub extern "C" fn on_signal(sig: libc::c_int) {
            let p = PATH.load(Ordering::SeqCst);
            // Safety: `p` is null or a C string that stays allocated while armed.
            unsafe {
                if !p.is_null() {
                    libc::unlink(p);
                }
                libc::signal(sig, libc::SIG_DFL);
                libc::raise(sig);
            }
        }
    }

    impl RemoveOnInterrupt {
        #[cfg(unix)]
        pub fn new(path: &Path) -> Self {
            use core::sync::atomic::Ordering;
            use std::os::unix::ffi::OsStrExt as _;
            let Ok(c) = alloc::ffi::CString::new(path.as_os_str().as_bytes()) else {
                return Self {
                    previous: Vec::new(),
                };
            };
            imp::PATH.store(c.into_raw(), Ordering::SeqCst);
            let handler: extern "C" fn(libc::c_int) = imp::on_signal;
            let mut previous = Vec::new();
            for s in imp::SIGNALS {
                // Safety: the handler only makes async-signal-safe calls; an
                // ignored signal is put straight back to ignored.
                unsafe {
                    let prev = libc::signal(s, handler as libc::sighandler_t);
                    if prev == libc::SIG_IGN {
                        libc::signal(s, libc::SIG_IGN);
                    } else {
                        previous.push((s, prev));
                    }
                }
            }
            Self { previous }
        }

        #[cfg(not(unix))]
        pub fn new(_: &Path) -> Self {
            Self {}
        }
    }

    impl Drop for RemoveOnInterrupt {
        fn drop(&mut self) {
            #[cfg(unix)]
            {
                use core::sync::atomic::Ordering;
                for &(s, prev) in &self.previous {
                    // Safety: restores the disposition `new` replaced.
                    unsafe {
                        libc::signal(s, prev);
                    }
                }
                let p = imp::PATH.swap(core::ptr::null_mut(), Ordering::SeqCst);
                if !p.is_null() {
                    // Safety: `p` came from `CString::into_raw` in `new`.
                    drop(unsafe { alloc::ffi::CString::from_raw(p) });
                }
            }
        }
    }
}

/// Index files in `dir`: every `*.arcsecix`, in name order.
pub fn index_files(dir: &Path) -> Vec<PathBuf> {
    let mut v: Vec<PathBuf> = std::fs::read_dir(dir)
        .map(|rd| {
            rd.filter_map(Result::ok)
                .map(|e| e.path())
                .filter(|p| {
                    p.extension()
                        .is_some_and(|e| e == arcsec_core::index::format::EXTENSION)
                })
                .collect()
        })
        .unwrap_or_default();
    v.sort();
    v
}

/// Leftovers of interrupted index writes in `dir` (`*.arcsecix.part`).
pub fn stale_parts(dir: &Path) -> Vec<PathBuf> {
    let suffix = format!(".{}.part", arcsec_core::index::format::EXTENSION);
    let mut v: Vec<PathBuf> = std::fs::read_dir(dir)
        .map(|rd| {
            rd.filter_map(Result::ok)
                .map(|e| e.path())
                .filter(|p| {
                    p.file_name()
                        .is_some_and(|n| n.to_string_lossy().ends_with(&suffix))
                })
                .collect()
        })
        .unwrap_or_default();
    v.sort();
    v
}

/// A phrase for an index's freshness against the database in `db_dir`, or `None`
/// when there is nothing to say (current).
pub fn freshness_note(ex: &Existing, db_dir: &Path) -> Option<String> {
    match freshness(ex, db_dir) {
        Freshness::Current => None,
        Freshness::Changed => Some(format!(
            "STALE: {} has changed since this index was built; rebuild it with `arcsec catalog index build`",
            ex.source.to_uppercase()
        )),
        Freshness::SourceMissing => Some(format!(
            "built from {}, which is not installed here (the index still works on its own)",
            ex.source.to_uppercase()
        )),
        Freshness::Unrecorded => Some(
            "built before arcsec recorded its source database's version, so it cannot be checked for staleness"
                .to_string(),
        ),
    }
}

/// Print one index's description.
pub fn describe(path: &Path) -> Result<(), String> {
    let ix = BlindIndex::open(path).map_err(|e| e.to_string())?;
    let (lo, hi) = super::plan::coverage_of(ix.tiers().iter().map(|t| t.radius.to_degrees()));
    println!(
        "{}: from {}, fields {lo:.2}°–{hi:.0}°, {} patterns, {} stars, {}",
        path.display(),
        ix.source().to_uppercase(),
        ix.n_patterns(),
        ix.n_stars(),
        human(ix.file_size() as u64)
    );
    for t in ix.tiers() {
        let r = t.radius.to_degrees();
        let (lo, hi) = tier_fov_range(r);
        println!(
            "    disc {r:>5.2}°  mag ≤ {:>4.1}  {:>9} anchors {:>10} patterns {:>9}  fields {lo:.2}°–{hi:.1}°",
            t.mag_cap,
            t.n_anchors,
            t.n_patterns,
            human(ix.tier_bytes(t)),
        );
    }
    if let Some(ex) = Existing::open(path)
        && let Some(note) = freshness_note(&ex, path.parent().unwrap_or(Path::new(".")))
    {
        println!("    {note}");
    }
    Ok(())
}

/// The one-line suggestion for a directory with databases but no index: the
/// command and what it will cost. `None` if no solving database is installed.
///
/// `flags` is appended to the command (`--db <dir>` when the database is not in the
/// catalogue directory).
pub fn build_suggestion(db_dir: &Path, flags: &str) -> Option<String> {
    let plan = Plan::for_databases(&installed_sources(db_dir), None, None)?;
    let est = Estimate::of(&plan, arcsec_core::max_threads());
    Some(format!(
        "build one with `arcsec catalog index build{flags}` (from {}: {})",
        plan.label(),
        est.summary()
    ))
}

/// `arcsec catalog index info`.
pub fn cmd_info(dir: &Path, file: Option<&PathBuf>) -> Result<(), String> {
    let files = file.map_or_else(|| index_files(dir), |f| vec![f.clone()]);
    if files.is_empty() {
        match build_suggestion(dir, "") {
            Some(s) => println!("No blind index in {}; {s}.", dir.display()),
            None => println!(
                "No blind index in {}. Install a star database first (`arcsec catalog install d50`); it builds one.",
                dir.display()
            ),
        }
        return Ok(());
    }
    for f in files {
        describe(&f)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::catalog_cmd::testutil::{TempDir, write_001_db};

    #[test]
    fn tier_choice_covers_the_range_without_extra_tiers() {
        let r = |a, b| -> Vec<f64> { tiers_for(a, b).iter().map(|t| t.radius_deg).collect() };
        // 0.3°–30°: 3° reaches 36°, 0.1° reaches down to 0.25°.
        assert_eq!(r(0.3, 30.0), vec![3.0, 1.5, 0.75, 0.4, 0.2, 0.1]);
        // A narrow band: the tiers that each cover all of it.
        assert_eq!(r(1.0, 2.0), vec![0.4, 0.2]);
        // Down to the D-series floor.
        assert_eq!(*r(0.15, 1.0).last().unwrap(), 0.06);
    }

    #[test]
    fn the_part_file_is_the_one_the_writer_uses() {
        assert_eq!(
            part_path(Path::new("/c/d80.arcsecix")),
            PathBuf::from("/c/d80.arcsecix.part")
        );
    }

    #[test]
    fn build_plans_default_to_what_is_installed() {
        let dir = TempDir::new("ix_plan");
        let d = dir.path();
        assert!(plan_for_build(d, None, None, None).is_err(), "no database");
        write_001_db(d, "w08", 50, 1);
        let p = plan_for_build(d, None, None, None).unwrap();
        assert_eq!(
            (p.source.as_str(), p.min_fov, p.max_fov),
            ("w08", 10.0, 80.0)
        );
        // A named database that is not the deepest: its own defaults.
        write_001_db(d, "g05", 50, 2);
        let p = plan_for_build(d, None, None, None).unwrap();
        assert_eq!(
            (p.source.as_str(), p.min_fov, p.max_fov),
            ("g05", 3.0, 80.0)
        );
        let p = plan_for_build(d, Some("W08"), None, None).unwrap();
        assert_eq!((p.source.as_str(), p.max_fov), ("w08", 80.0));
        let p = plan_for_build(d, None, Some(1.0), Some(5.0)).unwrap();
        assert_eq!((p.min_fov, p.max_fov), (1.0, 5.0));
        assert!(plan_for_build(d, Some("d80"), None, None).is_err());
        assert!(plan_for_build(d, None, Some(5.0), Some(1.0)).is_err());
    }

    /// The whole path on a tiny synthetic all-sky database: plan, estimate, build,
    /// a verified file with its source stamp, and staleness once the database
    /// changes. A stale `.part` from an interrupted write is cleared first.
    #[test]
    fn a_tiny_database_builds_a_stamped_index_that_goes_stale() {
        let dir = TempDir::new("ix_e2e");
        let d = dir.path();
        write_001_db(d, "w08", 3000, 7);
        let plan = plan_for_build(d, None, None, None).unwrap();
        let est = Estimate::of(&plan, 2);
        let out = default_index_path(d, &plan.source);
        std::fs::write(part_path(&out), b"half an index").unwrap();
        assert_eq!(stale_parts(d).len(), 1);

        let size = run_build(&plan, d, &out, 2, &est).unwrap();
        assert!(size > 256 && size == std::fs::metadata(&out).unwrap().len());
        assert!(stale_parts(d).is_empty(), "the stale .part is gone");
        let ix = BlindIndex::open(&out).unwrap();
        ix.validate().unwrap();
        assert!(ix.n_patterns() > 0);
        assert_eq!(ix.tiers().len(), 3);

        let ex = Existing::open(&out).unwrap();
        assert!(ex.stamp.is_recorded());
        assert_eq!(freshness(&ex, d), Freshness::Current);
        assert_eq!(freshness_note(&ex, d), None);
        assert_eq!(Existing::preferred(d).unwrap().path, out);
        assert_eq!(
            super::super::plan::rebuild_reason(Some(&ex), &plan, d),
            None
        );

        // A new copy of the database: the index is stale.
        write_001_db(d, "w08", 3001, 8);
        assert_eq!(freshness(&ex, d), Freshness::Changed);
        assert!(freshness_note(&ex, d).unwrap().starts_with("STALE"));
    }
}
