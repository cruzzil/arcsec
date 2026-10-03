//! arcsec's own blind index: which one to build for the databases installed, what
//! it will cost, whether the one there is current, and building it. See
//! `docs/offline-index.md` in the arcsec repository.
//!
//! The index format and the builder itself are in [`arcsec_core::index`]; this
//! module decides what to build and handles the file around the build. A typical
//! install:
//!
//! ```no_run
//! use arcsec_catalogue::{default_dir, index};
//!
//! let dir = default_dir();
//! let dbs = index::installed_sources(&dir);
//! let plan = index::Plan::for_databases(&dbs, None, None).expect("a solving database");
//! let machine = index::Machine::probe(&dir, 0);
//! let est = index::Estimate::of(&plan, machine.threads);
//! index::check_disk(&dir, est.bytes, machine.free_disk)?;
//! let out = arcsec_core::index::default_index_path(&dir, &plan.source);
//! let built = index::build_file(&plan, &dir, &out, machine.threads, |_| {})?;
//! println!("{} bytes, {} patterns", built.bytes, built.patterns);
//! # Ok::<(), arcsec_catalogue::Error>(())
//! ```

mod plan;

use std::path::{Path, PathBuf};

use arcsec_core::ArcsecError;
use arcsec_core::catalog::catalog_present;
use arcsec_core::index::{
    BuildParams, BuildProgress, BuiltIndex, DEFAULT_TIERS, TierSpec, build_index,
    default_index_path, tier_fov_range,
};

pub use arcsec_core::auto::{SOURCES, depth_rank, index_files, preferred_index};
pub use plan::{
    Concern, DISK_MARGIN, Estimate, Existing, Freshness, Machine, NOTICE_BYTES,
    NOTICE_RAM_FRACTION, NOTICE_SECS, Plan, Rebuild, check_disk, concerns, coverage_of,
    default_fields, disk_needed, freshness, keep_existing_range, rebuild_reason,
};

use crate::error::{Error, Result};
use crate::registry::{Entry, Purpose};

/// Solving databases installed in `dir`, deepest first.
#[must_use]
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

/// The temporary file [`BuiltIndex::write`] writes `out` through.
#[must_use]
pub fn part_path(out: &Path) -> PathBuf {
    out.with_extension(format!("{}.part", arcsec_core::index::format::EXTENSION))
}

/// The plan for building an index from the databases in `db_path`: `name` (or the
/// deepest installed database) as the source, and `min_fov`/`max_fov` or else the
/// fields every installed database covers.
///
/// # Errors
///
/// [`Error::NoDatabase`] if there is no solving database, [`Error::DatabaseMissing`]
/// if `name` is not installed, [`Error::BadFieldRange`] or [`Error::NoTier`] if the
/// range is empty or served by no tier.
pub fn plan_for_build(
    db_path: &Path,
    name: Option<&str>,
    min_fov: Option<f64>,
    max_fov: Option<f64>,
) -> Result<Plan> {
    let installed = installed_sources(db_path);
    let plan = match name.map(str::to_ascii_lowercase) {
        Some(n) if installed.first().is_none_or(|f| *f != n) => {
            if !catalog_present(db_path, &n) {
                return Err(Error::DatabaseMissing {
                    name: n,
                    dir: db_path.to_path_buf(),
                });
            }
            let (lo, hi) = default_fields(&n);
            Plan::new(&n, min_fov.unwrap_or(lo), max_fov.unwrap_or(hi))
        }
        _ => {
            Plan::for_databases(&installed, min_fov, max_fov).ok_or_else(|| Error::NoDatabase {
                dir: db_path.to_path_buf(),
            })?
        }
    };
    if !(plan.min_fov > 0.0 && plan.max_fov >= plan.min_fov) {
        return Err(Error::BadFieldRange {
            min: plan.min_fov,
            max: plan.max_fov,
        });
    }
    if plan.tiers.is_empty() {
        return Err(Error::NoTier {
            min: plan.min_fov,
            max: plan.max_fov,
            source: plan.source,
        });
    }
    Ok(plan)
}

/// The plan for the databases in `db_dir` and its estimated cost on
/// [`arcsec_core::max_threads`] threads: what building an index there with no
/// options would do. `None` if no solving database is installed.
#[must_use]
pub fn suggestion(db_dir: &Path) -> Option<(Plan, Estimate)> {
    let plan = Plan::for_databases(&installed_sources(db_dir), None, None)?;
    let est = Estimate::of(&plan, arcsec_core::max_threads());
    Some((plan, est))
}

// ── Building ───────────────────────────────────────────────────────────────────

/// Get ready to write an index to `out`: remove a `.part` file an interrupted write
/// left ([`part_path`]), create the directory, and check that a file can be
/// created there — before a build that may take minutes, not after it.
///
/// # Errors
///
/// [`Error::Io`] if the directory cannot be created, [`Error::NotWritable`] if no
/// file can be created in it (an ASTAP folder under Program Files or `/opt` is
/// read-only to a normal user).
pub fn prepare_output(out: &Path) -> Result<()> {
    let part = part_path(out);
    if part.exists() {
        let _ = std::fs::remove_file(&part);
    }
    if let Some(parent) = out.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent).map_err(Error::io(parent))?;
    }
    std::fs::File::create(&part)
        .and_then(|_| std::fs::remove_file(&part))
        .map_err(|source| Error::NotWritable {
            path: out.to_path_buf(),
            source,
        })
}

/// Build `plan` from the database in `db_path` on `threads` worker threads (0 =
/// [`arcsec_core::max_threads`]), reporting each tier and strip to `on_progress`.
/// The index is built in memory; nothing is written.
///
/// Stops with [`Error::Cancelled`] once the thread's [`arcsec_core::cancel`] token
/// is cancelled (checked between declination strips).
///
/// # Errors
///
/// [`Error::Build`] if the database cannot be read, or [`Error::Cancelled`].
pub fn build(
    plan: &Plan,
    db_path: &Path,
    threads: usize,
    on_progress: impl FnMut(&BuildProgress),
) -> Result<BuiltIndex> {
    build_index(
        &BuildParams {
            db_path: db_path.to_path_buf(),
            db_name: plan.source.clone(),
            tiers: plan.tiers.clone(),
            threads,
        },
        on_progress,
    )
    .map_err(|e| match e {
        ArcsecError::Cancelled => Error::Cancelled,
        e => Error::Build(e),
    })
}

/// Write a built index to `out`, through a `.part` file renamed into place so an
/// interrupted write never leaves a file that looks like an index. Returns the
/// file's size.
///
/// # Errors
///
/// [`Error::Write`] if the file cannot be written.
pub fn write(built: &BuiltIndex, out: &Path) -> Result<u64> {
    built.write(out).map_err(|source| Error::Write {
        path: out.to_path_buf(),
        source,
    })?;
    Ok(std::fs::metadata(out).map_or(0, |m| m.len()))
}

/// What [`build_file`] wrote.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Built {
    /// Size of the file, bytes.
    pub bytes: u64,
    /// Patterns in it.
    pub patterns: usize,
    /// Stars in it.
    pub stars: usize,
}

/// [`prepare_output`], [`build`] and [`write()`]: build `plan` from the database in
/// `db_path` and write it to `out`.
///
/// # Errors
///
/// As those three.
pub fn build_file(
    plan: &Plan,
    db_path: &Path,
    out: &Path,
    threads: usize,
    on_progress: impl FnMut(&BuildProgress),
) -> Result<Built> {
    prepare_output(out)?;
    let built = build(plan, db_path, threads, on_progress)?;
    let bytes = write(&built, out)?;
    Ok(Built {
        bytes,
        patterns: built.keys.len(),
        stars: built.stars.len(),
    })
}

// ── The index alongside installs and removals ──────────────────────────────────

/// Leftovers of interrupted index writes in `dir` (`*.arcsecix.part`), by name.
#[must_use]
pub fn stale_parts(dir: &Path) -> Vec<PathBuf> {
    let suffix = format!(".{}.part", arcsec_core::index::format::EXTENSION);
    let mut v: Vec<PathBuf> = std::fs::read_dir(dir)
        .map(|rd| {
            rd.filter_map(core::result::Result::ok)
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

/// Index files in `dir` built from database `db` (by their header).
#[must_use]
pub fn indexes_built_from(dir: &Path, db: &str) -> Vec<PathBuf> {
    index_files(dir)
        .into_iter()
        .filter(|p| Existing::open(p).is_some_and(|e| e.source.eq_ignore_ascii_case(db)))
        .collect()
}

/// Options for the index an install builds.
#[derive(Debug, Clone, Copy, Default)]
pub struct IndexOptions {
    /// Build no index.
    pub skip: bool,
    /// Smallest field to serve (short side, degrees), instead of the databases'.
    pub min_fov: Option<f64>,
    /// Largest field to serve, likewise.
    pub max_fov: Option<f64>,
}

/// A blind index an install should build once the downloads are in.
#[derive(Debug, Clone)]
pub struct IndexAction {
    /// What to build.
    pub plan: Plan,
    /// What it costs.
    pub est: Estimate,
    /// Why: no index yet, a deeper database, a stale or too-narrow index.
    pub reason: Rebuild,
    /// Indexes in the directory the new one supersedes, to remove once it is built.
    pub replaces: Vec<PathBuf>,
}

/// The solving databases `dir` will hold once `adding` are installed, deepest
/// first.
#[must_use]
pub fn sources_after_install(dir: &Path, adding: &[&Entry]) -> Vec<&'static str> {
    let mut dbs = installed_sources(dir);
    for e in adding {
        if let Some(s) = SOURCES.iter().find(|s| **s == e.id)
            && !dbs.contains(s)
        {
            dbs.push(s);
        }
    }
    dbs.sort_by_key(|d| depth_rank(d));
    dbs
}

/// The index to build after an install so that `dbs` (every solving database that
/// will then be in `dir`, as [`sources_after_install`] gives) are served, costed on
/// `threads` threads; `None` if the one there already serves them, or
/// `opts.skip`.
#[must_use]
pub fn index_action(
    dir: &Path,
    dbs: &[&str],
    opts: &IndexOptions,
    threads: usize,
) -> Option<IndexAction> {
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
        .map(|db| default_index_path(dir, db))
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

/// The indexes in `dir` that go when `removing` are removed: those built from a
/// solving database among them. They would still solve, but nothing would update
/// or verify them any more, and they are usually the bigger files.
#[must_use]
pub fn indexes_removed_with(dir: &Path, removing: &[&Entry]) -> Vec<PathBuf> {
    removing
        .iter()
        .filter(|e| e.purpose == Purpose::Solving)
        .flat_map(|e| indexes_built_from(dir, e.id))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::registry::find;
    use crate::test_support::{TempDir, write_001_db};

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
        assert!(matches!(
            plan_for_build(d, None, None, None),
            Err(Error::NoDatabase { .. })
        ));
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
        assert!(matches!(
            plan_for_build(d, Some("d80"), None, None),
            Err(Error::DatabaseMissing { .. })
        ));
        assert!(matches!(
            plan_for_build(d, None, Some(5.0), Some(1.0)),
            Err(Error::BadFieldRange { .. })
        ));
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
        let out = default_index_path(d, &plan.source);
        std::fs::write(part_path(&out), b"half an index").unwrap();
        assert_eq!(stale_parts(d).len(), 1);

        let mut tiers = 0;
        let built = build_file(&plan, d, &out, 2, |p| {
            if matches!(p, BuildProgress::TierDone { .. }) {
                tiers += 1;
            }
        })
        .unwrap();
        assert_eq!(tiers, 3);
        assert!(built.bytes > 256 && built.bytes == std::fs::metadata(&out).unwrap().len());
        assert!(stale_parts(d).is_empty(), "the stale .part is gone");
        let ix = arcsec_core::index::BlindIndex::open(&out).unwrap();
        ix.validate().unwrap();
        assert!(ix.n_patterns() > 0 && ix.n_patterns() == built.patterns);
        assert_eq!(ix.tiers().len(), 3);

        let ex = Existing::open(&out).unwrap();
        assert!(ex.stamp.is_recorded());
        assert_eq!(freshness(&ex, d), Freshness::Current);
        assert_eq!(Existing::preferred(d).unwrap().path, out);
        assert_eq!(rebuild_reason(Some(&ex), &plan, d), None);
        assert_eq!(indexes_built_from(d, "W08"), vec![out.clone()]);
        assert_eq!(
            indexes_removed_with(d, &[find("w08").unwrap()]),
            vec![out.clone()]
        );
        assert!(indexes_removed_with(d, &[find("v05").unwrap()]).is_empty());

        // A new copy of the database: the index is stale.
        write_001_db(d, "w08", 3001, 8);
        assert_eq!(freshness(&ex, d), Freshness::Changed);
    }

    #[test]
    fn a_cancelled_build_writes_nothing() {
        let dir = TempDir::new("ix_cancel");
        let d = dir.path();
        write_001_db(d, "w08", 500, 3);
        let plan = plan_for_build(d, None, None, None).unwrap();
        let out = default_index_path(d, &plan.source);
        let token = arcsec_core::cancel::CancelToken::new();
        token.cancel();
        let r = arcsec_core::cancel::with_token(&token, || build_file(&plan, d, &out, 1, |_| {}));
        assert!(matches!(r, Err(Error::Cancelled)), "{r:?}");
        assert!(!out.exists() && stale_parts(d).is_empty());
    }

    #[test]
    fn an_unwritable_destination_is_found_before_the_build() {
        let dir = TempDir::new("ix_ro");
        let blocker = dir.path().join("file");
        std::fs::write(&blocker, b"x").unwrap();
        // A directory that cannot exist: its parent is a file.
        let r = prepare_output(&blocker.join("sub").join("d80.arcsecix"));
        assert!(matches!(r, Err(Error::Io { .. })), "{r:?}");
    }

    #[test]
    fn install_index_options_shape_the_plan() {
        let dir = TempDir::new("inst_opts");
        let d = dir.path();
        let opts = IndexOptions {
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
        assert!(index_action(d, &["d80"], &IndexOptions { skip: true, ..opts }, 8).is_none());
    }

    #[test]
    fn the_databases_after_an_install_are_deepest_first() {
        let dir = TempDir::new("after");
        let d = dir.path();
        write_001_db(d, "w08", 50, 1);
        let adding = [
            find("g05").unwrap(),
            find("v05").unwrap(),
            find("d50").unwrap(),
        ];
        assert_eq!(sources_after_install(d, &adding), ["d50", "g05", "w08"]);
        assert_eq!(sources_after_install(d, &[]), ["w08"]);
    }
}
