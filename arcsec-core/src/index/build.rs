//! Building a blind index from an installed ASTAP star database.
//!
//! **Disc-anchored patterns.** For each tier — a disc radius `r` and a magnitude
//! cap — a star *anchors* patterns only if it is the brightest star within `r` of
//! itself (among stars no fainter than the cap). Its group is itself plus the next
//! brightest stars of its disc, and every 4-subset of the group is a pattern. An
//! image wide enough to contain the disc sees the same locally-brightest stars, so
//! it can rebuild the same group without knowing where it is; the 4-subsets give
//! redundancy against one member going undetected or ranking differently in the
//! image. The idea comes from seiza's blind index (Apache-2.0); this is an
//! independent implementation over ASTAP's databases.
//!
//! No two anchors share a pattern: a pattern contains its anchor, and an anchor
//! cannot be inside another anchor's disc (one of the two would be brighter). So
//! nothing needs de-duplicating within a tier.
//!
//! **Memory.** The sky is processed a declination strip at a time, reading the
//! strip plus a margin of `r` from the database, so even the deepest tier never
//! holds the whole catalogue: peak memory is the strip's stars plus the tier's
//! patterns (24 bytes each).

use core::f64::consts::{FRAC_PI_2, PI};
use std::collections::HashMap;
use std::path::Path;

use super::format::{BuiltIndex, IndexStar, STAR_BANDS, TierInfo, star_band};
use super::pattern::{Tangent, canonical, descriptor, key, unit};
use crate::catalog::format_1476::for_each_star_in_dec_band;
use crate::error::Result;

/// One tier to build.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TierSpec {
    /// Disc radius, degrees.
    pub radius_deg: f64,
    /// Faintest magnitude used.
    pub mag_cap: f64,
    /// Group size, anchor included (5 gives 5 patterns per anchor, 6 gives 15).
    pub members: usize,
}

/// The tier ladder. Each disc radius is paired with the magnitude at which a disc
/// of that size holds roughly 15–20 stars on average, so the anchor's group is
/// drawn from stars an image of that field will show. Radii step by about 2×, so
/// any field from ~0.15° up sees at least one tier at a comfortable size (a tier
/// is used for fields between about 2.5 and 12 disc radii across).
///
/// Measured on D80 (Gaia BP): mag ≤ 6.1 → 0.1 stars/deg², ≤ 9.2 → 7, ≤ 12.7 →
/// 150, ≤ 14.2 → 600, ≤ 16 → 1900. The wider tiers use 6-star groups (15
/// patterns per anchor); the two deepest use 5 (5 patterns), which keeps the
/// index size in check where anchors are most numerous.
pub const DEFAULT_TIERS: [TierSpec; 9] = [
    TierSpec {
        radius_deg: 12.0,
        mag_cap: 4.6,
        members: 6,
    },
    TierSpec {
        radius_deg: 6.0,
        mag_cap: 6.1,
        members: 6,
    },
    TierSpec {
        radius_deg: 3.0,
        mag_cap: 7.6,
        members: 6,
    },
    TierSpec {
        radius_deg: 1.5,
        mag_cap: 9.2,
        members: 6,
    },
    TierSpec {
        radius_deg: 0.75,
        mag_cap: 10.7,
        members: 6,
    },
    TierSpec {
        radius_deg: 0.4,
        mag_cap: 11.8,
        members: 6,
    },
    TierSpec {
        radius_deg: 0.2,
        mag_cap: 12.7,
        members: 6,
    },
    TierSpec {
        radius_deg: 0.1,
        mag_cap: 14.2,
        members: 5,
    },
    TierSpec {
        radius_deg: 0.06,
        mag_cap: 16.0,
        members: 5,
    },
];

/// The smallest image field (degrees, short side) a tier serves, and the largest.
#[must_use]
pub fn tier_fov_range(radius_deg: f64) -> (f64, f64) {
    (2.5 * radius_deg, 12.0 * radius_deg)
}

/// What to build.
#[derive(Debug, Clone)]
pub struct BuildParams {
    /// Directory holding the database.
    pub db_path: std::path::PathBuf,
    /// Database name (`d80`, `g05`, ...).
    pub db_name: String,
    /// Tiers, any order (they are stored widest first).
    pub tiers: Vec<TierSpec>,
    /// Worker threads; 0 = [`crate::max_threads`].
    pub threads: usize,
}

/// Progress reported by [`build_index`].
#[derive(Debug, Clone)]
pub enum BuildProgress {
    /// A tier is starting.
    Tier {
        /// Its position (0-based) and the number of tiers.
        index: usize,
        /// Number of tiers.
        of: usize,
        /// The tier.
        spec: TierSpec,
    },
    /// A declination strip of the current tier is done.
    Strip {
        /// Strips done.
        done: usize,
        /// Strips in the tier.
        of: usize,
        /// Patterns so far in the tier.
        patterns: usize,
    },
    /// A tier is finished.
    TierDone {
        /// The finished tier.
        info: TierInfo,
        /// Stars read from the database for it.
        stars_read: usize,
    },
}

/// A star of the strip being processed.
#[derive(Clone, Copy)]
struct Src {
    u: [f64; 3],
    ra: f64,
    dec: f64,
    mag: f64,
}

/// Grid of a strip's stars in cells about `r` on a side, for disc queries.
struct Grid {
    r: f64,
    cells: HashMap<(i32, i32), Vec<u32>>,
}

impl Grid {
    fn band(&self, dec: f64) -> i32 {
        ((dec + FRAC_PI_2) / self.r).floor() as i32
    }
    fn ra_cells(&self, band: i32) -> i32 {
        let lo = f64::from(band) * self.r - FRAC_PI_2;
        let hi = lo + self.r;
        // Narrowest cell at the band's edge nearest the equator, so a cell is never
        // narrower than r on the sky.
        let cos_max = if lo <= 0.0 && hi >= 0.0 {
            1.0
        } else {
            lo.cos().max(hi.cos())
        };
        ((2.0 * PI * cos_max / self.r).floor() as i32).max(1)
    }
    fn cell(&self, ra: f64, dec: f64) -> (i32, i32) {
        let b = self.band(dec);
        let n = self.ra_cells(b);
        let c =
            ((ra.rem_euclid(2.0 * PI) / (2.0 * PI) * f64::from(n)).floor() as i32).rem_euclid(n);
        (b, c)
    }
    fn new(r: f64, stars: &[Src]) -> Self {
        let mut g = Self {
            r,
            cells: HashMap::new(),
        };
        for (i, s) in stars.iter().enumerate() {
            let k = g.cell(s.ra, s.dec);
            g.cells.entry(k).or_default().push(i as u32);
        }
        g
    }
    /// Every star index that could lie within `r` of (`ra`, `dec`).
    fn near(&self, ra: f64, dec: f64, mut f: impl FnMut(u32)) {
        let b0 = self.band(dec);
        for b in b0 - 1..=b0 + 1 {
            let n = self.ra_cells(b);
            let lo = f64::from(b) * self.r - FRAC_PI_2;
            let hi = lo + self.r;
            if hi < -FRAC_PI_2 - self.r || lo > FRAC_PI_2 + self.r {
                continue;
            }
            // Worst-case RA half-width of a disc of radius r anywhere in this band.
            let d_extreme = dec.abs().max(lo.abs()).max(hi.abs()).min(FRAC_PI_2);
            let cos_d = d_extreme.cos();
            let all = cos_d * PI <= self.r * 1.01 || n <= 3;
            let (c_lo, c_hi) = if all {
                (0, n - 1)
            } else {
                let half = (self.r / cos_d).min(PI);
                let w = 2.0 * PI / f64::from(n);
                let c = ra.rem_euclid(2.0 * PI) / w;
                ((c - half / w).floor() as i32, (c + half / w).floor() as i32)
            };
            let span = (c_hi - c_lo + 1).min(n);
            for k in 0..span {
                let c = (c_lo + k).rem_euclid(n);
                if let Some(v) = self.cells.get(&(b, c)) {
                    for &i in v {
                        f(i);
                    }
                }
            }
        }
    }
}

/// Patterns found for one anchor: its group (strip-local indices) and the keys of
/// its 4-subsets with their members in canonical order.
type AnchorOut = (Vec<u32>, Vec<(u64, [u32; 4])>);

/// Build the group and patterns of star `i`, if it is an anchor.
fn anchor_patterns(
    i: usize,
    stars: &[Src],
    grid: &Grid,
    cos_r: f64,
    members: usize,
) -> Option<AnchorOut> {
    let s = &stars[i];
    let mut group: Vec<u32> = Vec::new();
    let mut brighter = false;
    grid.near(s.ra, s.dec, |j| {
        if brighter || j as usize == i {
            return;
        }
        let o = &stars[j as usize];
        if s.u[0] * o.u[0] + s.u[1] * o.u[1] + s.u[2] * o.u[2] >= cos_r {
            if (j as usize) < i {
                brighter = true;
            } else {
                group.push(j);
            }
        }
    });
    if brighter || group.len() < 3 {
        return None;
    }
    group.sort_unstable();
    group.truncate(members - 1);
    group.insert(0, i as u32);

    let m = group.len();
    let mut pats = Vec::new();
    for a in 0..m {
        for b in a + 1..m {
            for c in b + 1..m {
                for d in c + 1..m {
                    let ids = [group[a], group[b], group[c], group[d]];
                    let mut csum = [0.0f64; 3];
                    for &id in &ids {
                        let u = stars[id as usize].u;
                        csum[0] += u[0];
                        csum[1] += u[1];
                        csum[2] += u[2];
                    }
                    let Some(tp) = Tangent::at(csum) else {
                        continue;
                    };
                    let mut p = [(0.0, 0.0); 4];
                    let mut ok = true;
                    for (k, &id) in ids.iter().enumerate() {
                        match tp.project(&stars[id as usize].u) {
                            Some(xy) => p[k] = xy,
                            None => ok = false,
                        }
                    }
                    if !ok {
                        continue;
                    }
                    let Some((order, _, _)) = canonical(&p) else {
                        continue;
                    };
                    let ordered = order.map(|k| p[k]);
                    pats.push((key(&descriptor(&ordered)), order.map(|k| ids[k])));
                }
            }
        }
    }
    Some((group, pats))
}

/// Read one strip's stars (core plus margin), sorted brightest first with a total
/// order, so "brighter" means the same thing in every strip.
fn read_strip(params: &BuildParams, lo: f64, hi: f64, cap: f64) -> Result<Vec<Src>> {
    let mut v = Vec::new();
    for_each_star_in_dec_band(&params.db_path, &params.db_name, lo, hi, cap, |s| {
        v.push(Src {
            u: unit(s.ra, s.dec),
            ra: s.ra,
            dec: s.dec,
            mag: s.mag,
        });
    })?;
    v.sort_by(|a, b| {
        a.mag
            .total_cmp(&b.mag)
            .then(a.ra.total_cmp(&b.ra))
            .then(a.dec.total_cmp(&b.dec))
    });
    // A star present twice (the band reader is called once per strip, but a
    // database may hold an exact duplicate) would block its own anchoring.
    v.dedup_by(|a, b| a.ra == b.ra && a.dec == b.dec);
    Ok(v)
}

/// Build an index. Deterministic: the same database and tiers give the same file
/// (apart from the build time), whatever the thread count.
///
/// # Errors
///
/// [`crate::ArcsecError::CatalogIo`] if a database file cannot be read.
pub fn build_index(
    params: &BuildParams,
    mut progress: impl FnMut(&BuildProgress),
) -> Result<BuiltIndex> {
    let threads = if params.threads > 0 {
        params.threads
    } else {
        crate::max_threads()
    }
    .max(1);

    let mut tiers = params.tiers.clone();
    tiers.sort_by(|a, b| b.radius_deg.total_cmp(&a.radius_deg));

    let mut out = BuiltIndex {
        source: params.db_name.clone(),
        ..BuiltIndex::default()
    };
    // Stars as appended (duplicates across strips are merged at the end).
    let mut stars: Vec<IndexStar> = Vec::new();

    for (ti, spec) in tiers.iter().enumerate() {
        progress(&BuildProgress::Tier {
            index: ti,
            of: tiers.len(),
            spec: *spec,
        });
        let r = spec.radius_deg.to_radians();
        let cos_r = r.cos();
        // Wide tiers have few stars: one strip. Deep ones go 5° at a time.
        let strip = if spec.mag_cap <= 11.0 {
            PI
        } else {
            5f64.to_radians()
        };
        let n_strips = (PI / strip).ceil() as usize;
        let first_pattern = out.keys.len();
        let mut tier_pats: Vec<(u64, [u32; 4])> = Vec::new();
        let mut n_anchors = 0u64;
        let mut stars_read = 0usize;

        for si in 0..n_strips {
            let lo = -FRAC_PI_2 + si as f64 * strip;
            let hi = (lo + strip).min(FRAC_PI_2);
            let src = read_strip(params, lo - r, hi + r, spec.mag_cap)?;
            stars_read += src.len();
            let grid = Grid::new(r, &src);
            // Anchors are the strip's own stars; the margin only supplies
            // neighbours. The top strip owns +90° itself.
            let core: Vec<usize> = (0..src.len())
                .filter(|&i| src[i].dec >= lo && (src[i].dec < hi || si + 1 == n_strips))
                .collect();
            let chunk = core.len().div_ceil(threads).max(1);
            let results: Vec<Vec<AnchorOut>> = std::thread::scope(|scope| {
                let handles: Vec<_> = core
                    .chunks(chunk)
                    .map(|part| {
                        let (src, grid) = (&src, &grid);
                        scope.spawn(move || {
                            part.iter()
                                .filter_map(|&i| anchor_patterns(i, src, grid, cos_r, spec.members))
                                .collect::<Vec<_>>()
                        })
                    })
                    .collect();
                handles
                    .into_iter()
                    .map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
                    .collect()
            });

            // Strip-local star index → global (pre-merge) index.
            let mut local: HashMap<u32, u32> = HashMap::new();
            for (group, pats) in results.into_iter().flatten() {
                n_anchors += 1;
                for &l in &group {
                    local.entry(l).or_insert_with(|| {
                        let s = &src[l as usize];
                        stars.push(IndexStar {
                            ra: s.ra as f32,
                            dec: s.dec as f32,
                            mag: (s.mag * 100.0).round().clamp(-32768.0, 32767.0) as i16,
                            tier: ti as u8,
                        });
                        (stars.len() - 1) as u32
                    });
                }
                for (k, ids) in pats {
                    tier_pats.push((k, ids.map(|l| local[&l])));
                }
            }
            progress(&BuildProgress::Strip {
                done: si + 1,
                of: n_strips,
                patterns: tier_pats.len(),
            });
        }

        tier_pats.sort_unstable();
        let info = TierInfo {
            radius: r,
            mag_cap: spec.mag_cap as f32,
            members: spec.members as u32,
            first_pattern: first_pattern as u64,
            n_patterns: tier_pats.len() as u64,
            n_anchors,
        };
        out.keys.extend(tier_pats.iter().map(|p| p.0));
        out.quads.extend(tier_pats.iter().map(|p| p.1));
        drop(tier_pats);
        out.tiers.push(info);
        progress(&BuildProgress::TierDone { info, stars_read });
    }

    // Merge stars: sort by (band, RA, Dec), collapse exact duplicates (the same
    // database record reached from two strips or two tiers, keeping the widest
    // tier), and renumber the quads.
    let mut order: Vec<u32> = (0..stars.len() as u32).collect();
    let band_of = |s: &IndexStar| star_band(f64::from(s.dec));
    order.sort_unstable_by(|&a, &b| {
        let (x, y) = (&stars[a as usize], &stars[b as usize]);
        band_of(x)
            .cmp(&band_of(y))
            .then(x.ra.total_cmp(&y.ra))
            .then(x.dec.total_cmp(&y.dec))
            .then(x.tier.cmp(&y.tier))
    });
    let mut remap = vec![0u32; stars.len()];
    let mut merged: Vec<IndexStar> = Vec::with_capacity(stars.len());
    for &i in &order {
        let s = stars[i as usize];
        match merged.last() {
            Some(m) if m.ra.to_bits() == s.ra.to_bits() && m.dec.to_bits() == s.dec.to_bits() => {}
            _ => merged.push(s),
        }
        remap[i as usize] = (merged.len() - 1) as u32;
    }
    drop(stars);
    for q in &mut out.quads {
        *q = q.map(|i| remap[i as usize]);
    }
    let mut dir = vec![0u32; STAR_BANDS as usize + 1];
    for s in &merged {
        dir[star_band(f64::from(s.dec)) as usize + 1] += 1;
    }
    for b in 1..dir.len() {
        dir[b] += dir[b - 1];
    }
    out.star_dir = dir;
    out.stars = merged;
    Ok(out)
}

/// Default index file for database `db_name` in directory `dir`.
#[must_use]
pub fn default_index_path(dir: &Path, db_name: &str) -> std::path::PathBuf {
    dir.join(format!("{db_name}.{}", super::format::EXTENSION))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::{Rng, SkyStar, TempDir, write_1476_db};

    fn deg(d: f64) -> f64 {
        d.to_radians()
    }

    /// Six stars inside a 0.1° disc round (ra, dec), brightest at the centre.
    fn disc(ra: f64, dec: f64) -> Vec<SkyStar> {
        let offs = [
            (0.0, 0.0),
            (0.05, 0.01),
            (-0.03, 0.04),
            (0.02, -0.06),
            (-0.06, -0.02),
            (0.01, 0.07),
        ];
        offs.iter()
            .enumerate()
            .map(|(i, &(dx, dy))| SkyStar {
                ra: deg(ra + dx / deg(dec).cos()).rem_euclid(2.0 * PI),
                dec: deg(dec + dy),
                mag: 9.0 + 0.3 * i as f64,
            })
            .collect()
    }

    fn build(dir: &Path, stars: &[SkyStar], threads: usize) -> BuiltIndex {
        write_1476_db(dir, "t", stars);
        build_index(
            &BuildParams {
                db_path: dir.to_path_buf(),
                db_name: "t".into(),
                tiers: vec![TierSpec {
                    radius_deg: 0.1,
                    mag_cap: 14.0, // deep: processed in 5° strips
                    members: 6,
                }],
                threads,
            },
            |_| {},
        )
        .unwrap()
    }

    #[test]
    fn a_lone_disc_gives_all_its_four_subsets_across_strip_and_ra_seams() {
        // At a strip boundary (Dec 0°), across RA 0°, and in open sky.
        for (ra, dec) in [(0.0, 0.0), (120.0, 30.0), (240.0, -5.0)] {
            let dir = TempDir::new("ixbuild_disc");
            let ix = build(dir.path(), &disc(ra, dec), 2);
            assert_eq!(ix.tiers[0].n_anchors, 1, "({ra}, {dec})");
            assert_eq!(ix.keys.len(), 15, "({ra}, {dec})");
            assert_eq!(ix.stars.len(), 6);
            for q in &ix.quads {
                let mut s = q.to_vec();
                s.sort_unstable();
                s.dedup();
                assert_eq!(s.len(), 4, "four distinct stars");
            }
        }
    }

    #[test]
    fn only_the_brightest_star_of_a_disc_anchors() {
        let dir = TempDir::new("ixbuild_two");
        let mut stars = disc(50.0, 20.0);
        stars.extend(disc(50.0, 20.4));
        let ix = build(dir.path(), &stars, 1);
        assert_eq!(ix.tiers[0].n_anchors, 2);
        // A brighter star beside the second disc's centre takes over its anchoring.
        stars.push(SkyStar {
            ra: deg(50.0),
            dec: deg(20.405),
            mag: 5.0,
        });
        let dir2 = TempDir::new("ixbuild_two_b");
        let ix2 = build(dir2.path(), &stars, 1);
        assert_eq!(ix2.tiers[0].n_anchors, 2);
        assert!(ix2.stars.iter().any(|s| s.mag == 500));
    }

    #[test]
    fn the_result_does_not_depend_on_the_thread_count() {
        let mut rng = Rng::new(7);
        let stars: Vec<SkyStar> = (0..3000)
            .map(|_| SkyStar {
                ra: deg(rng.range(10.0, 14.0)),
                dec: deg(rng.range(-6.0, 3.0)),
                mag: rng.range(6.0, 14.0),
            })
            .collect();
        let (d1, d2) = (TempDir::new("ixbuild_t1"), TempDir::new("ixbuild_t2"));
        let a = build(d1.path(), &stars, 1);
        let b = build(d2.path(), &stars, 5);
        assert!(a.keys.len() > 100);
        assert_eq!(a.keys, b.keys);
        assert_eq!(a.quads, b.quads);
        assert_eq!(a.stars, b.stars);
        assert_eq!(a.star_dir, b.star_dir);
        assert!(a.keys.windows(2).all(|w| w[0] <= w[1]));
        assert_eq!(*a.star_dir.last().unwrap() as usize, a.stars.len());
    }
}
