//! `arcsec catalog index build|info` — arcsec's own blind index, built from a star
//! database the user already has. See `docs/offline-index.md`.
//!
//! `catalog install` builds one automatically after downloading a solving database
//! (see `cmd_install`); this module runs the build with its progress report, and
//! describes an installed index. What to build and the build itself are
//! `arcsec_catalogue::index`'s.

use std::io::IsTerminal as _;
use std::path::{Path, PathBuf};
use std::time::Instant;

use arcsec_catalogue::index::{
    self as index, Concern, Estimate, Existing, Freshness, Machine, Plan, check_disk, concerns,
    coverage_of, freshness, index_files, plan_for_build,
};
use arcsec_catalogue::{duration, human_bytes as human};
use arcsec_core::index::{BlindIndex, BuildProgress, default_index_path, tier_fov_range};

use super::message;
use super::prompt::{Terminal, build_question};

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
pub fn print_concerns(concerns: &[Concern]) {
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
    let plan = plan_for_build(&db_path, name.map(String::as_str), min_fov, max_fov)
        .map_err(|e| message(&e))?;
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

    check_disk(out_dir, est.bytes, machine.free_disk).map_err(|e| message(&e))?;
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
    index::prepare_output(out).map_err(|e| message(&e))?;

    let t0 = Instant::now();
    let mut progress = Progress::new(est, t0);
    let built =
        index::build(plan, db_path, threads, |p| progress.event(p)).map_err(|e| message(&e))?;

    let size = {
        let _guard = interrupt::RemoveOnInterrupt::new(&index::part_path(out));
        index::write(&built, out).map_err(|e| message(&e))?
    };
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
    let (lo, hi) = coverage_of(ix.tiers().iter().map(|t| t.radius.to_degrees()));
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
    let (plan, est) = index::suggestion(db_dir)?;
    Some(format!(
        "build one with `arcsec catalog index build{flags}` (from {}: {})",
        plan.label(),
        est.summary()
    ))
}

/// A one-line stderr hint when a search is wide enough that an installed index
/// would be used, but there is none: say how to build one and what it costs. The
/// solve itself is unchanged; nothing is built during a solve.
pub fn missing_index_hint(
    explicit: Option<&PathBuf>,
    template: &arcsec_core::pipeline::SolveParams,
) -> Option<String> {
    use arcsec_core::auto::{find_arcsec_index, wants_installed_index};
    if explicit.is_some() || !wants_installed_index(template) {
        return None;
    }
    let cat_dir = super::default_dir();
    if find_arcsec_index(&cat_dir).is_some() || find_arcsec_index(&template.db_path).is_some() {
        return None;
    }
    let same = |a: &Path, b: &Path| match (a.canonicalize(), b.canonicalize()) {
        (Ok(x), Ok(y)) => x == y,
        _ => a == b,
    };
    let db_flag = if same(&template.db_path, &cat_dir) {
        String::new()
    } else {
        format!(" --db {}", template.db_path.display())
    };
    let s = build_suggestion(&template.db_path, &db_flag)?;
    Some(format!(
        "Hint: no blind index is installed, so a search this wide can take minutes; {s}."
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
