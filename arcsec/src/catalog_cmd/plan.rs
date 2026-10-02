//! Which blind index to build for the databases a user has, and what it will cost.
//!
//! Every solving database can seed arcsec's blind index (docs/offline-index.md), and
//! the tiers down to a 0.1° disc come out practically identical whichever one is
//! used: they draw on stars to magnitude 14.2, which even the 500 stars/deg² D05 and
//! G05 hold almost everywhere. So the choice is mostly *which fields to cover*, and
//! that follows the databases installed: a D-series database solves from 0.15–0.6°
//! up, G05 from 3°, W08 from 20°. The index is built from the deepest database
//! installed, for the union of their field ranges.
//!
//! The cost estimate is a model fitted to real builds (2026-10-02, 24-thread x86-64
//! desktop, `docs/offline-index.md` §2.5): per tier, the stars the builder reads, the
//! patterns it keeps and the stars it stores, measured on D80, D05, G05 and W08.
//! D20 and D50 use D80's numbers, which over-estimate them slightly.

use std::path::{Path, PathBuf};

use arcsec_core::index::{BlindIndex, SourceStamp, TierSpec, tier_fov_range};

use super::human;
use super::index_cmd::tiers_for;

/// Databases an index can be built from, deepest first: the first installed one is
/// the default source, and an index built from an earlier one is preferred.
pub const SOURCES: [&str; 6] = ["d80", "d50", "d20", "d05", "g05", "w08"];

/// Position of `db` in [`SOURCES`] (deeper is smaller); unknown names sort last.
pub fn depth_rank(db: &str) -> usize {
    SOURCES
        .iter()
        .position(|s| s.eq_ignore_ascii_case(db))
        .unwrap_or(SOURCES.len())
}

/// Fields (short side, degrees) an index built for database `db` covers by default.
///
/// The low end is where the database itself stops being useful for verification
/// (a blind position still has to be confirmed by the hinted solver against the
/// same database), except that D80 and D50 stop at 0.3°: the 0.06° tier that
/// D80's 0.15–0.3° fields need costs 410 MB of the 698 MB index and finds
/// half of those fields blind (docs/offline-index.md §7.1), so it is offered, not
/// built unasked. The high end is 30° for everything below W08: the 3° and 1.5°
/// tiers cost under 3 MB together.
pub fn default_fields(db: &str) -> (f64, f64) {
    match db.to_ascii_lowercase().as_str() {
        "d05" => (0.6, 30.0),
        "g05" => (3.0, 30.0),
        "w08" => (10.0, 80.0),
        _ => (0.3, 30.0),
    }
}

/// The widest tier disc (degrees) a database is too shallow to fill. W08 stops at
/// magnitude 8, so a 1.5° disc (magnitude cap 9.2) would be built from an
/// incomplete star list whose groups an image cannot rebuild.
fn shallowest_useless_radius(db: &str) -> Option<f64> {
    db.eq_ignore_ascii_case("w08").then_some(1.5)
}

/// What to build: a source database and a field range, and the tiers for it.
#[derive(Debug, Clone, PartialEq)]
pub struct Plan {
    /// Database the stars come from.
    pub source: String,
    /// Smallest field (short side, degrees).
    pub min_fov: f64,
    /// Largest field (short side, degrees).
    pub max_fov: f64,
    /// Tiers, widest first.
    pub tiers: Vec<TierSpec>,
}

impl Plan {
    /// A plan for `source` covering `min_fov`–`max_fov`, minus any tier the database
    /// is too shallow for.
    pub fn new(source: &str, min_fov: f64, max_fov: f64) -> Self {
        let mut tiers = tiers_for(min_fov, max_fov);
        if let Some(r) = shallowest_useless_radius(source) {
            tiers.retain(|t| t.radius_deg > r);
        }
        Self {
            source: source.to_ascii_lowercase(),
            min_fov,
            max_fov,
            tiers,
        }
    }

    /// The plan for a set of installed solving databases: built from the deepest,
    /// covering every one's default field range. `min`/`max` override the range.
    pub fn for_databases(dbs: &[&str], min: Option<f64>, max: Option<f64>) -> Option<Self> {
        let source = SOURCES.into_iter().find(|s| dbs.contains(s))?;
        let lo = dbs
            .iter()
            .map(|d| default_fields(d).0)
            .fold(f64::INFINITY, f64::min);
        let hi = dbs
            .iter()
            .map(|d| default_fields(d).1)
            .fold(f64::NEG_INFINITY, f64::max);
        Some(Self::new(source, min.unwrap_or(lo), max.unwrap_or(hi)))
    }

    /// The fields the plan's tiers actually serve, which is a little wider than
    /// asked for: (smallest, largest) short side in degrees.
    pub fn coverage(&self) -> (f64, f64) {
        coverage_of(self.tiers.iter().map(|t| t.radius_deg))
    }

    /// One line: `D80, fields 0.3°–30°`.
    pub fn label(&self) -> String {
        format!(
            "{}, fields {}°–{}°",
            self.source.to_uppercase(),
            trim(self.min_fov),
            trim(self.max_fov)
        )
    }
}

/// Fields served by tiers of these disc radii (degrees).
pub fn coverage_of(radii: impl Iterator<Item = f64>) -> (f64, f64) {
    radii
        .map(tier_fov_range)
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), (a, b)| {
            (lo.min(a), hi.max(b))
        })
}

/// `0.3`, `30`, `0.15`: a degree value without trailing zeros.
fn trim(v: f64) -> String {
    let s = format!("{v:.3}");
    s.trim_end_matches('0').trim_end_matches('.').to_string()
}

// ── Cost model ─────────────────────────────────────────────────────────────────

/// What one tier costs to build from a given database.
#[derive(Debug, Clone, Copy)]
struct TierCost {
    /// Patterns kept.
    patterns: f64,
    /// Stars the tier stores (before merging with other tiers).
    stars: f64,
    /// Stars read from the database, margins included.
    read: f64,
}

/// Measured on D80 (Gaia to 8000 stars/deg²): disc radius → cost.
const D80_COSTS: [(f64, TierCost); 9] = [
    (12.0, tc(1_166.0, 440.0, 974.0)),
    (6.0, tc(5_212.0, 1_864.0, 5_165.0)),
    (3.0, tc(21_404.0, 7_595.0, 26_921.0)),
    (1.5, tc(87_162.0, 30_471.0, 143_476.0)),
    (0.75, tc(349_669.0, 122_537.0, 625_508.0)),
    (0.4, tc(1_185_086.0, 421_308.0, 2_006_011.0)),
    (0.2, tc(3_639_264.0, 1_427_214.0, 4_121_724.0)),
    (0.1, tc(4_428_404.0, 4_285_114.0, 13_777_880.0)),
    (0.06, tc(12_900_826.0, 12_303_302.0, 51_962_040.0)),
];

const fn tc(patterns: f64, stars: f64, read: f64) -> TierCost {
    TierCost {
        patterns,
        stars,
        read,
    }
}

/// A tier's cost when built from `db`.
fn tier_cost(db: &str, radius: f64) -> TierCost {
    let base = D80_COSTS
        .iter()
        .min_by(|a, b| (a.0 - radius).abs().total_cmp(&(b.0 - radius).abs()))
        .map_or(tc(0.0, 0.0, 0.0), |c| c.1);
    let near = |r: f64| (radius - r).abs() < 1e-6;
    match db.to_ascii_lowercase().as_str() {
        // 500 stars/deg²: the same tiers down to 0.1°, read from fewer stars; the
        // 0.06° disc holds too few stars for every anchor to keep a group. Measured.
        "d05" | "g05" if near(0.1) => TierCost {
            read: 10_960_000.0,
            ..base
        },
        "d05" | "g05" if near(0.06) => tc(10_110_000.0, 10_540_000.0, 18_160_000.0),
        // Magnitude 8: only the widest tiers are complete. Measured.
        "w08" if radius < 2.0 => tc(base.patterns.min(55_000.0), 15_000.0, 45_000.0),
        _ => base,
    }
}

/// Fraction of the shallower tiers' stars that are not already in the deepest
/// tier (fitted: the default D80 build stores 4.52 M stars against 4.29 M in its
/// deepest tier and 6.29 M summed over all of them).
const SHARED_STAR_FRACTION: f64 = 0.116;

/// Build time per star read, seconds: a serial part and a part spread over the
/// worker threads. Fitted to fourteen builds (1–24 threads, 4–107 s); within ±30 %
/// but for one 2-thread run 30 % slower than modelled.
const SECS_PER_STAR_SERIAL: f64 = 1.06e-6;
const SECS_PER_STAR_PARALLEL: f64 = 2.5e-6;

/// Peak resident memory: a base plus bytes per pattern and per stored star. Fitted
/// to the same builds (peak RSS 90 MB–1.28 GB); within ±5% above 200 MB.
const RAM_BASE: f64 = 60e6;
const RAM_PER_PATTERN: f64 = 37.0;
const RAM_PER_STAR: f64 = 29.0;

/// Index file layout constants (see `arcsec_core::index::format`).
const FILE_FIXED: f64 = 256.0 + 721.0 * 4.0;
const FILE_PER_TIER: f64 = 40.0;
const FILE_PER_PATTERN: f64 = 24.0;
const FILE_PER_STAR: f64 = 12.0;

/// What building an index is expected to cost.
#[derive(Debug, Clone, PartialEq)]
pub struct Estimate {
    /// Size of the index file, bytes.
    pub bytes: u64,
    /// Build time, seconds.
    pub secs: f64,
    /// Peak memory, bytes.
    pub ram: u64,
    /// Patterns.
    pub patterns: u64,
    /// Fraction of the build time each tier takes, in the plan's order.
    pub tier_share: Vec<f64>,
}

impl Estimate {
    /// Estimate the cost of `plan` on `threads` worker threads.
    pub fn of(plan: &Plan, threads: usize) -> Self {
        let threads = threads.max(1) as f64;
        let costs: Vec<TierCost> = plan
            .tiers
            .iter()
            .map(|t| tier_cost(&plan.source, t.radius_deg))
            .collect();
        let patterns: f64 = costs.iter().map(|c| c.patterns).sum();
        let deepest = costs.iter().map(|c| c.stars).fold(0.0, f64::max);
        let all: f64 = costs.iter().map(|c| c.stars).sum();
        let stars = deepest + SHARED_STAR_FRACTION * (all - deepest);
        let per_star = SECS_PER_STAR_SERIAL + SECS_PER_STAR_PARALLEL / threads;
        let tier_secs: Vec<f64> = costs.iter().map(|c| c.read * per_star).collect();
        let secs: f64 = tier_secs.iter().sum();
        let tier_share = tier_secs
            .iter()
            .map(|t| if secs > 0.0 { t / secs } else { 0.0 })
            .collect();
        Self {
            bytes: (FILE_FIXED
                + FILE_PER_TIER * plan.tiers.len() as f64
                + FILE_PER_PATTERN * patterns
                + FILE_PER_STAR * stars) as u64,
            secs,
            ram: (RAM_BASE + RAM_PER_PATTERN * patterns + RAM_PER_STAR * stars) as u64,
            patterns: patterns as u64,
            tier_share,
        }
    }

    /// `~287 MB on disk, ~30 s, ~0.6 GB memory`.
    pub fn summary(&self) -> String {
        format!(
            "~{} on disk, {}, ~{} memory",
            human(self.bytes),
            duration(self.secs),
            gigabytes(self.ram)
        )
    }
}

/// A duration for a prompt, rounded the way a person would say it.
pub fn duration(secs: f64) -> String {
    if secs < 10.0 {
        "under 10 s".to_string()
    } else if secs < 55.0 {
        format!("~{} s", ((secs / 5.0).round() * 5.0) as u64)
    } else if secs < 3600.0 {
        format!("~{} min", ((secs / 60.0).round() as u64).max(1))
    } else {
        let m = (secs / 60.0).round() as u64;
        format!("~{} h {} min", m / 60, m % 60)
    }
}

/// `0.6 GB`, or `90 MB` below 0.1 GB.
pub fn gigabytes(bytes: u64) -> String {
    if bytes < 100_000_000 {
        human(bytes)
    } else {
        format!("{:.1} GB", bytes as f64 / 1e9)
    }
}

// ── When to say more than the usual prompt ─────────────────────────────────────

/// An index bigger than this on disk gets its own notice and question.
pub const NOTICE_BYTES: u64 = 1_000_000_000;
/// A build expected to take longer than this, seconds, likewise.
pub const NOTICE_SECS: f64 = 300.0;
/// Likewise a build whose peak memory exceeds this fraction of available memory.
pub const NOTICE_RAM_FRACTION: f64 = 0.5;

/// Headroom kept free on the disk beyond the index itself: 10 % plus this.
pub const DISK_MARGIN: u64 = 100_000_000;

/// What the machine has to spare (see `sys`), injectable for tests.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Machine {
    /// Free bytes on the catalogue's disk, if known.
    pub free_disk: Option<u64>,
    /// Available memory, if known.
    pub avail_ram: Option<u64>,
    /// Worker threads the build would use.
    pub threads: usize,
}

impl Machine {
    /// Probe the real machine for building into `dir` with `threads` (0 = all).
    pub fn probe(dir: &Path, threads: usize) -> Self {
        Self {
            free_disk: super::sys::free_space(dir),
            avail_ram: super::sys::available_memory(),
            threads: if threads == 0 {
                arcsec_core::max_threads()
            } else {
                threads
            },
        }
    }
}

/// Bytes that must be free to write `bytes` of new files.
pub fn disk_needed(bytes: u64) -> u64 {
    bytes + bytes / 10 + DISK_MARGIN
}

/// `Err` with a message if `needed` bytes will not fit on a disk with `free` bytes
/// (unknown free space never blocks).
pub fn check_disk(dir: &Path, needed: u64, free: Option<u64>) -> Result<(), String> {
    match free {
        Some(f) if f < disk_needed(needed) => Err(format!(
            "not enough free space in {}: this needs about {} and {} is free",
            dir.display(),
            human(disk_needed(needed)),
            human(f)
        )),
        _ => Ok(()),
    }
}

/// Why a build deserves its own notice, if it does: one phrase per reason.
pub fn concerns(est: &Estimate, m: &Machine) -> Vec<String> {
    let mut v = Vec::new();
    if est.bytes > NOTICE_BYTES {
        v.push(format!("it needs {} of disk", human(est.bytes)));
    }
    if est.secs > NOTICE_SECS {
        v.push(format!(
            "it takes {} on {} threads",
            duration(est.secs),
            m.threads
        ));
    }
    if let Some(avail) = m.avail_ram {
        if est.ram as f64 > avail as f64 {
            v.push(format!(
                "it needs ~{} of memory and only {} is available; the build may fail or swap heavily",
                gigabytes(est.ram),
                gigabytes(avail)
            ));
        } else if est.ram as f64 > NOTICE_RAM_FRACTION * avail as f64 {
            v.push(format!(
                "it needs ~{} of memory, more than half of the {} available",
                gigabytes(est.ram),
                gigabytes(avail)
            ));
        }
    }
    v
}

// ── Existing indexes ───────────────────────────────────────────────────────────

/// An index already in the catalogue directory.
#[derive(Debug, Clone)]
pub struct Existing {
    /// The file.
    pub path: PathBuf,
    /// Database it was built from.
    pub source: String,
    /// Fingerprint of that database at build time.
    pub stamp: SourceStamp,
    /// Fields its tiers serve (short side, degrees).
    pub coverage: (f64, f64),
}

impl Existing {
    /// Read an index's header; `None` if it is not a usable index.
    pub fn open(path: &Path) -> Option<Self> {
        let ix = BlindIndex::open(path).ok()?;
        Some(Self {
            path: path.to_path_buf(),
            source: ix.source().to_ascii_lowercase(),
            stamp: ix.source_stamp(),
            coverage: coverage_of(ix.tiers().iter().map(|t| t.radius.to_degrees())),
        })
    }

    /// The index the solver uses in `dir`: the one built from the deepest database
    /// (by [`SOURCES`]), then by file name.
    pub fn preferred(dir: &Path) -> Option<Self> {
        let mut v: Vec<Self> = super::index_cmd::index_files(dir)
            .iter()
            .filter_map(|p| Self::open(p))
            .collect();
        v.sort_by(|a, b| {
            depth_rank(&a.source)
                .cmp(&depth_rank(&b.source))
                .then_with(|| a.path.cmp(&b.path))
        });
        v.into_iter().next()
    }

    /// Whether its tiers serve every field of `plan` (with 1 % slack).
    pub fn covers(&self, plan: &Plan) -> bool {
        let (lo, hi) = plan.coverage();
        self.coverage.0 <= lo * 1.01 && self.coverage.1 >= hi * 0.99
    }
}

/// Whether an index still matches the database it was built from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Freshness {
    /// The database is installed and unchanged since the build.
    Current,
    /// Its files differ from those the index was built from.
    Changed,
    /// The database is not in the directory (the index still works on its own).
    SourceMissing,
    /// The index predates source fingerprints (arcsec 0.4 and earlier).
    Unrecorded,
}

/// Compare an index's recorded source with the database in `db_dir`.
pub fn freshness(ex: &Existing, db_dir: &Path) -> Freshness {
    if !arcsec_core::catalog::catalog_present(db_dir, &ex.source) {
        return Freshness::SourceMissing;
    }
    if !ex.stamp.is_recorded() {
        return Freshness::Unrecorded;
    }
    match SourceStamp::of_database(db_dir, &ex.source) {
        Ok(now) if now == ex.stamp => Freshness::Current,
        _ => Freshness::Changed,
    }
}

/// Why the index in a directory should be (re)built.
#[derive(Debug, Clone, PartialEq)]
pub enum Rebuild {
    /// There is no index.
    Missing,
    /// A deeper database than the index's source is installed.
    Deeper {
        /// The deeper database.
        new: String,
        /// The index's source.
        old: String,
    },
    /// The source database has changed since the build.
    Changed(String),
    /// The index does not serve every installed database's fields.
    Narrow {
        /// Fields the index serves.
        have: (f64, f64),
        /// Fields wanted.
        want: (f64, f64),
    },
}

impl core::fmt::Display for Rebuild {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Missing => write!(f, "no blind index is installed yet"),
            Self::Deeper { new, old } => write!(
                f,
                "{} is deeper than {}, which the current index was built from",
                new.to_uppercase(),
                old.to_uppercase()
            ),
            Self::Changed(db) => write!(
                f,
                "{} has changed since the index was built",
                db.to_uppercase()
            ),
            Self::Narrow { have, want } => write!(
                f,
                "the current index serves fields {}°–{}°, not {}°–{}°",
                trim(have.0),
                trim(have.1),
                trim(want.0),
                trim(want.1)
            ),
        }
    }
}

/// Why the index in a directory should be (re)built for `plan`, or `None` if the one
/// there already serves it.
pub fn rebuild_reason(existing: Option<&Existing>, plan: &Plan, db_dir: &Path) -> Option<Rebuild> {
    let Some(ex) = existing else {
        return Some(Rebuild::Missing);
    };
    if depth_rank(&plan.source) < depth_rank(&ex.source) {
        return Some(Rebuild::Deeper {
            new: plan.source.clone(),
            old: ex.source.clone(),
        });
    }
    if freshness(ex, db_dir) == Freshness::Changed {
        return Some(Rebuild::Changed(ex.source.clone()));
    }
    if !ex.covers(plan) {
        return Some(Rebuild::Narrow {
            have: ex.coverage,
            want: plan.coverage(),
        });
    }
    None
}

/// Widen `plan` to keep whatever the existing index already served, so a rebuild
/// prompted by a new database never drops tiers the user chose (`--min-fov 0.15`).
/// Bounds the user set explicitly are kept as given.
pub fn keep_existing_range(
    plan: Plan,
    existing: Option<&Existing>,
    min_given: bool,
    max_given: bool,
) -> Plan {
    let Some(ex) = existing else {
        return plan;
    };
    let lo = if min_given {
        plan.min_fov
    } else {
        plan.min_fov.min(ex.coverage.0)
    };
    let hi = if max_given {
        plan.max_fov
    } else {
        plan.max_fov.max(ex.coverage.1 / 1.2)
    };
    if (lo - plan.min_fov).abs() < 1e-9 && (hi - plan.max_fov).abs() < 1e-9 {
        return plan;
    }
    let widened = Plan::new(&plan.source, lo, hi);
    // Only widen if it changes the tiers; otherwise keep the user-facing range.
    if widened.tiers == plan.tiers {
        plan
    } else {
        widened
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::catalog_cmd::testutil::TempDir;
    use arcsec_core::index::{BuiltIndex, TierInfo};

    fn radii(p: &Plan) -> Vec<f64> {
        p.tiers.iter().map(|t| t.radius_deg).collect()
    }

    #[test]
    fn plans_follow_the_installed_databases() {
        let p = Plan::for_databases(&["d80"], None, None).unwrap();
        assert_eq!(p.source, "d80");
        assert_eq!(radii(&p), vec![3.0, 1.5, 0.75, 0.4, 0.2, 0.1]);

        // Deepest source, union of the ranges.
        let p = Plan::for_databases(&["g05", "w08", "d50"], None, None).unwrap();
        assert_eq!(p.source, "d50");
        assert_eq!((p.min_fov, p.max_fov), (0.3, 80.0));
        assert_eq!(radii(&p), vec![12.0, 6.0, 3.0, 1.5, 0.75, 0.4, 0.2, 0.1]);

        assert_eq!(
            radii(&Plan::for_databases(&["d05"], None, None).unwrap()),
            vec![3.0, 1.5, 0.75, 0.4, 0.2]
        );
        assert_eq!(
            radii(&Plan::for_databases(&["g05"], None, None).unwrap()),
            vec![3.0, 1.5, 0.75]
        );
        // W08 is too shallow for anything narrower than the 3° disc.
        assert_eq!(
            radii(&Plan::for_databases(&["w08"], None, None).unwrap()),
            vec![12.0, 6.0, 3.0]
        );
        assert_eq!(
            radii(&Plan::for_databases(&["w08"], Some(0.3), None).unwrap()),
            vec![12.0, 6.0, 3.0]
        );
        // Overrides.
        let p = Plan::for_databases(&["d80"], Some(0.15), None).unwrap();
        assert_eq!(*radii(&p).last().unwrap(), 0.06);
        // Photometric or blind-only sets give no plan.
        assert!(Plan::for_databases(&["v05", "anet-4100"], None, None).is_none());
        assert!(Plan::for_databases(&[], None, None).is_none());
    }

    /// The estimator against the builds it was fitted to (24 threads unless noted):
    /// sizes are exact to a few percent, memory to ±10 %, time to ±35 %.
    #[test]
    fn the_estimate_matches_measured_builds() {
        let within = |what: &str, est: f64, actual: f64, tol: f64| {
            let r = est / actual;
            assert!(
                (1.0 - tol..=1.0 + tol).contains(&r),
                "{what}: estimate {est:.3e}, measured {actual:.3e} (ratio {r:.2})"
            );
        };
        // (db, min, max, threads, bytes, secs, peak RSS)
        let builds = [
            ("d80", 0.3, 30.0, 24, 287.3e6, 22.0, 548e6),
            ("d80", 0.3, 30.0, 4, 287.3e6, 36.9, 488e6),
            ("d80", 0.15, 30.0, 24, 697.5e6, 107.3, 1276e6),
            ("d05", 0.6, 30.0, 24, 144.5e6, 6.55, 303e6),
            ("d05", 0.3, 30.0, 24, 287.3e6, 20.5, 544e6),
            ("g05", 0.15, 144.0, 24, 610.0e6, 46.0, 1097e6),
        ];
        for (db, lo, hi, threads, bytes, secs, ram) in builds {
            let e = Estimate::of(&Plan::new(db, lo, hi), threads);
            let what = format!("{db} {lo}–{hi} @{threads}");
            within(&format!("{what} size"), e.bytes as f64, bytes, 0.03);
            within(&format!("{what} time"), e.secs, secs, 0.35);
            within(&format!("{what} memory"), e.ram as f64, ram, 0.15);
        }
        // Tiny builds stay tiny.
        let w = Estimate::of(&Plan::for_databases(&["w08"], None, None).unwrap(), 4);
        assert!(w.bytes < 2_000_000 && w.secs < 1.0, "{w:?}");
        let g = Estimate::of(&Plan::for_databases(&["g05"], None, None).unwrap(), 4);
        within("g05 default size", g.bytes as f64, 12.5e6, 0.05);
        // Tier shares add up.
        let d = Estimate::of(&Plan::new("d80", 0.15, 30.0), 8);
        assert!((d.tier_share.iter().sum::<f64>() - 1.0).abs() < 1e-9);
        assert!(
            *d.tier_share.last().unwrap() > 0.5,
            "the deepest tier dominates"
        );
    }

    #[test]
    fn durations_read_naturally() {
        assert_eq!(duration(3.0), "under 10 s");
        assert_eq!(duration(22.0), "~20 s");
        assert_eq!(duration(107.0), "~2 min");
        assert_eq!(duration(782.0), "~13 min");
        assert_eq!(duration(5400.0), "~1 h 30 min");
        assert_eq!(gigabytes(1_276_000_000), "1.3 GB");
        assert_eq!(gigabytes(90_000_000), "90.0 MB");
    }

    fn est(bytes: u64, secs: f64, ram: u64) -> Estimate {
        Estimate {
            bytes,
            secs,
            ram,
            patterns: 0,
            tier_share: vec![],
        }
    }

    #[test]
    fn notices_trigger_on_size_time_and_memory_only() {
        let roomy = Machine {
            free_disk: Some(500_000_000_000),
            avail_ram: Some(16_000_000_000),
            threads: 8,
        };
        assert!(concerns(&est(287_000_000, 30.0, 600_000_000), &roomy).is_empty());
        assert_eq!(
            concerns(&est(1_500_000_000, 30.0, 6e8 as u64), &roomy).len(),
            1
        );
        assert_eq!(
            concerns(&est(2e8 as u64, 900.0, 6e8 as u64), &roomy).len(),
            1
        );
        // A Raspberry Pi: 1.2 GB free, so a 0.6 GB build is over half.
        let pi = Machine {
            avail_ram: Some(1_200_000_000),
            ..roomy
        };
        let c = concerns(&est(287_000_000, 30.0, 650_000_000), &pi);
        assert!(c.len() == 1 && c[0].contains("more than half"), "{c:?}");
        let c = concerns(&est(287_000_000, 30.0, 1_300_000_000), &pi);
        assert!(c[0].contains("may fail"), "{c:?}");
        // Unknown memory never raises a concern.
        let unknown = Machine {
            avail_ram: None,
            ..roomy
        };
        assert!(concerns(&est(287_000_000, 30.0, u64::MAX / 2), &unknown).is_empty());
    }

    #[test]
    fn the_disk_check_refuses_only_when_it_knows() {
        let d = Path::new("/cat");
        assert!(check_disk(d, 287_000_000, Some(10_000_000_000)).is_ok());
        let e = check_disk(d, 287_000_000, Some(300_000_000)).unwrap_err();
        assert!(e.contains("not enough free space"), "{e}");
        assert!(check_disk(d, 287_000_000, None).is_ok());
        assert_eq!(disk_needed(1_000_000_000), 1_200_000_000);
    }

    /// Write a small index with the given tiers and source, as the builder would.
    fn write_index(path: &Path, source: &str, radii: &[f64], stamp: SourceStamp) {
        let ix = BuiltIndex {
            tiers: radii
                .iter()
                .map(|r| TierInfo {
                    radius: r.to_radians(),
                    mag_cap: 10.0,
                    members: 6,
                    first_pattern: 0,
                    n_patterns: 0,
                    n_anchors: 0,
                })
                .collect(),
            star_dir: vec![0; arcsec_core::index::format::STAR_BANDS as usize + 1],
            source: source.into(),
            source_stamp: stamp,
            ..BuiltIndex::default()
        };
        ix.write(path).unwrap();
    }

    #[test]
    fn rebuilds_are_asked_for_only_when_something_changed() {
        let dir = TempDir::new("plan_rebuild");
        let d = dir.path();
        std::fs::write(d.join("d80_0101.1476"), vec![7u8; 200]).unwrap();
        let stamp = SourceStamp::of_database(d, "d80").unwrap();
        let plan = Plan::for_databases(&["d80"], None, None).unwrap();

        assert_eq!(rebuild_reason(None, &plan, d), Some(Rebuild::Missing));

        let p = d.join("d80.arcsecix");
        write_index(&p, "d80", &[3.0, 1.5, 0.75, 0.4, 0.2, 0.1], stamp);
        let ex = Existing::open(&p).unwrap();
        assert_eq!(freshness(&ex, d), Freshness::Current);
        assert_eq!(rebuild_reason(Some(&ex), &plan, d), None);

        // W08 added: the index stops at 36°, the plan wants 80°.
        let wide = Plan::for_databases(&["d80", "w08"], None, None).unwrap();
        let r = rebuild_reason(Some(&ex), &wide, d).unwrap();
        assert!(matches!(r, Rebuild::Narrow { .. }), "{r:?}");
        assert!(
            r.to_string().contains("serves fields 0.25°–36°, not"),
            "{r}"
        );

        // The database changed underneath it.
        std::fs::write(d.join("d80_0101.1476"), vec![8u8; 200]).unwrap();
        assert_eq!(freshness(&ex, d), Freshness::Changed);
        assert_eq!(
            rebuild_reason(Some(&ex), &plan, d),
            Some(Rebuild::Changed("d80".into()))
        );

        // An index from before fingerprints: not stale, just unrecorded.
        write_index(
            &p,
            "d80",
            &[3.0, 1.5, 0.75, 0.4, 0.2, 0.1],
            SourceStamp::default(),
        );
        let old = Existing::open(&p).unwrap();
        assert_eq!(freshness(&old, d), Freshness::Unrecorded);
        assert_eq!(rebuild_reason(Some(&old), &plan, d), None);

        // A G05 index when D80 arrives: rebuild from the deeper database.
        let g = d.join("g05.arcsecix");
        write_index(&g, "g05", &[3.0, 1.5, 0.75], SourceStamp::default());
        let gx = Existing::open(&g).unwrap();
        assert_eq!(freshness(&gx, d), Freshness::SourceMissing);
        assert!(matches!(
            rebuild_reason(Some(&gx), &plan, d),
            Some(Rebuild::Deeper { .. })
        ));

        // The solver prefers the deepest source, whatever the names.
        assert_eq!(Existing::preferred(d).unwrap().source, "d80");
        std::fs::remove_file(&p).unwrap();
        assert_eq!(Existing::preferred(d).unwrap().source, "g05");
    }

    #[test]
    fn a_rebuild_keeps_tiers_the_user_added() {
        let dir = TempDir::new("plan_keep");
        let p = dir.path().join("d80.arcsecix");
        // Built with --min-fov 0.15.
        write_index(
            &p,
            "d80",
            &[3.0, 1.5, 0.75, 0.4, 0.2, 0.1, 0.06],
            SourceStamp::default(),
        );
        let ex = Existing::open(&p).unwrap();
        let plan = Plan::for_databases(&["d80", "w08"], None, None).unwrap();
        let kept = keep_existing_range(plan.clone(), Some(&ex), false, false);
        assert_eq!(*radii(&kept).last().unwrap(), 0.06);
        assert_eq!(radii(&kept)[0], 12.0);
        // An explicit --index-min-fov wins.
        let given = keep_existing_range(plan, Some(&ex), true, false);
        assert_eq!(*radii(&given).last().unwrap(), 0.1);
        // Nothing to keep: unchanged.
        let p2 = Plan::for_databases(&["d80"], None, None).unwrap();
        assert_eq!(keep_existing_range(p2.clone(), None, false, false), p2);
    }
}
