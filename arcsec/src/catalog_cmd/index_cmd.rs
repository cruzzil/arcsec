//! `arcsec catalog index build|info` — arcsec's own blind index, built from a star
//! database the user already has. See `docs/offline-index.md`.

use std::path::{Path, PathBuf};
use std::time::Instant;

use arcsec_core::catalog::catalog_present;
use arcsec_core::index::{
    BlindIndex, BuildParams, BuildProgress, DEFAULT_TIERS, TierSpec, build_index,
    default_index_path, tier_fov_range,
};

use super::human;

/// Databases an index can be built from, deepest first: the first installed one is
/// the default source.
const SOURCES: [&str; 6] = ["d80", "d50", "d20", "d05", "g05", "w08"];

/// The deepest installed database in `dir`, if any.
pub fn default_source(dir: &Path) -> Option<&'static str> {
    SOURCES.into_iter().find(|n| catalog_present(dir, n))
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

/// `arcsec catalog index build`.
pub fn cmd_build(
    cat_dir: &Path,
    db: Option<&PathBuf>,
    name: Option<&String>,
    min_fov: f64,
    max_fov: f64,
    out: Option<&PathBuf>,
    threads: usize,
) -> Result<(), String> {
    if !(min_fov > 0.0 && max_fov >= min_fov) {
        return Err(format!("bad field range {min_fov}°–{max_fov}°"));
    }
    let db_path = db.cloned().unwrap_or_else(|| cat_dir.to_path_buf());
    let db_name = match name {
        Some(n) => n.to_lowercase(),
        None => default_source(&db_path)
            .ok_or_else(|| {
                format!(
                    "no star database in {}; install one (`arcsec catalog install d50`) or pass --db and --name",
                    db_path.display()
                )
            })?
            .to_string(),
    };
    if !catalog_present(&db_path, &db_name) {
        return Err(format!(
            "database {db_name} not found in {}",
            db_path.display()
        ));
    }
    let out = out
        .cloned()
        .unwrap_or_else(|| default_index_path(cat_dir, &db_name));
    let tiers = tiers_for(min_fov, max_fov);
    if tiers.is_empty() {
        return Err("no tier fits that field range".into());
    }

    println!(
        "Building a blind index from {} in {}",
        db_name.to_uppercase(),
        db_path.display()
    );
    println!(
        "  fields {min_fov}°–{max_fov}° (short side), {} tiers:",
        tiers.len()
    );
    for t in &tiers {
        let (lo, hi) = tier_fov_range(t.radius_deg);
        println!(
            "    disc {:>5}°  mag ≤ {:>4}  groups of {}  (fields {lo:.2}°–{hi:.1}°)",
            t.radius_deg, t.mag_cap, t.members
        );
    }
    println!("  output {}", out.display());

    let t0 = Instant::now();
    let mut t_tier = Instant::now();
    let built = build_index(
        &BuildParams {
            db_path,
            db_name: db_name.clone(),
            tiers,
            threads,
        },
        |p| match p {
            BuildProgress::Tier { index, of, spec } => {
                t_tier = Instant::now();
                eprint!("  tier {}/{of} (disc {}°) ", index + 1, spec.radius_deg);
            }
            BuildProgress::Strip { done, of, .. } => {
                if *of > 1 && done % 6 == 0 {
                    eprint!(".");
                }
            }
            BuildProgress::TierDone { info, stars_read } => {
                eprintln!(
                    " {} anchors, {} patterns, {} ({stars_read} stars read, {:.1} s)",
                    info.n_anchors,
                    info.n_patterns,
                    human(info.n_patterns * 24),
                    t_tier.elapsed().as_secs_f64()
                );
            }
        },
    )
    .map_err(|e| format!("build failed: {e}"))?;
    if let Some(parent) = out.parent() {
        std::fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
    }
    built
        .write(&out)
        .map_err(|e| format!("writing {}: {e}", out.display()))?;
    let size = std::fs::metadata(&out).map_or(0, |m| m.len());
    println!(
        "Wrote {} ({}, {} patterns, {} stars) in {:.1} s.",
        out.display(),
        human(size),
        built.keys.len(),
        built.stars.len(),
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

/// Index files to describe: `file`, or every `*.arcsecix` in `dir`.
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

/// Print one index's description.
pub fn describe(path: &Path) -> Result<(), String> {
    let ix = BlindIndex::open(path).map_err(|e| e.to_string())?;
    println!(
        "{}: from {}, {} patterns, {} stars, {}",
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
    Ok(())
}

/// `arcsec catalog index info`.
pub fn cmd_info(dir: &Path, file: Option<&PathBuf>) -> Result<(), String> {
    let files = file.map_or_else(|| index_files(dir), |f| vec![f.clone()]);
    if files.is_empty() {
        println!(
            "No blind index in {}. Build one with `arcsec catalog index build`.",
            dir.display()
        );
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
}
