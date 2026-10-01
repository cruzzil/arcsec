//! A vote accumulator over field hypotheses: (RA, Dec, log pixel scale).
//!
//! The blind solvers turn every pattern match into a full field hypothesis — a
//! centre on the sky and a pixel scale — and the true field is the one that many
//! independent matches agree on. This collects those hypotheses into buckets and
//! ranks regions of agreement:
//!
//! * **Buckets are square on the sky.** The RA bin width is divided by `cos(Dec)`
//!   of its declination band, so a bucket near a pole spans as much sky as one on
//!   the equator. Binning raw RA split one field's votes over many buckets near
//!   the poles.
//! * **Scale is a third axis.** Two matches that imply the same centre at very
//!   different pixel scales are not agreeing on a field.
//! * **Neighbour smoothing.** Correct matches imply centres that scatter across
//!   a bucket edge, so each bucket is scored by the votes of its 3×3×3
//!   neighbourhood. Uniform junk stays flat; a split peak recombines.
//! * **The representative comes from the strongest bucket.** A smoothed region's
//!   hypothesis is taken from whichever bucket in its neighbourhood holds the most
//!   direct votes, not from the bucket the smoothing happened to centre on. Taking
//!   it from the centre bucket — or, as `blind.rs` once did, the first hypothesis
//!   ever deposited in a cell — can score a region on a stray member.
//! * **Non-maximum suppression.** A strong region floods its neighbours with the
//!   same smoothed sum, so once a region is ranked its surroundings are dropped and
//!   each region is verified once.

use core::f64::consts::PI;
use std::collections::{HashMap, HashSet};

/// Bucket key: (Dec band, RA bin within the band, log-scale bin).
type Key = (i32, i32, i32);

/// A ranked region of agreement.
#[derive(Debug)]
pub(crate) struct Region<'a, T> {
    /// Votes in the region's 3×3×3 neighbourhood.
    pub votes: usize,
    /// Hypotheses of the strongest bucket in the neighbourhood, in insertion order.
    pub members: &'a [T],
}

/// Hypotheses binned by (RA, Dec, ln scale).
pub(crate) struct SkyVotes<T> {
    /// Bucket side on the sky, radians.
    step: f64,
    /// Bucket width in `ln(scale)`.
    log_step: f64,
    buckets: HashMap<Key, Vec<T>>,
}

impl<T> SkyVotes<T> {
    /// An empty accumulator with `step` (radians) buckets on the sky and
    /// `log_step` buckets in natural-log pixel scale.
    pub(crate) fn new(step: f64, log_step: f64) -> Self {
        Self {
            step,
            log_step,
            buckets: HashMap::new(),
        }
    }

    /// Number of non-empty buckets.
    pub(crate) fn len(&self) -> usize {
        self.buckets.len()
    }

    /// RA bins in the Dec band `band`: at least one, and the band's bins exactly
    /// tile 2π so the last bin wraps onto the first.
    fn ra_bins(&self, band: i32) -> i32 {
        let dec_mid = (f64::from(band) + 0.5) * self.step - PI / 2.0;
        let width = self.step / dec_mid.cos().max(1e-6);
        ((2.0 * PI / width).floor() as i32).max(1)
    }

    fn key(&self, ra: f64, dec: f64, scale: f64) -> Key {
        let n_bands = (PI / self.step).ceil() as i32;
        let band = (((dec + PI / 2.0) / self.step).floor() as i32).clamp(0, n_bands - 1);
        let n_ra = self.ra_bins(band);
        let ra = ra.rem_euclid(2.0 * PI);
        let ra_bin = ((ra / (2.0 * PI) * f64::from(n_ra)).floor() as i32).rem_euclid(n_ra);
        let s_bin = (scale.max(1e-12).ln() / self.log_step).floor() as i32;
        (band, ra_bin, s_bin)
    }

    /// Centre of a bucket on the sky.
    fn centre(&self, key: Key) -> (f64, f64) {
        let dec = (f64::from(key.0) + 0.5) * self.step - PI / 2.0;
        let n_ra = self.ra_bins(key.0);
        let ra = (f64::from(key.1) + 0.5) * 2.0 * PI / f64::from(n_ra);
        (ra, dec)
    }

    /// The buckets within `reach` bins of `key` on the sky and `s_reach` in scale,
    /// `key` included. Neighbouring bands are found by stepping the bucket centre,
    /// since their RA bins are a different width.
    fn neighbours(&self, key: Key, reach: i32, s_reach: i32) -> Vec<Key> {
        let (ra, dec) = self.centre(key);
        let mut out = Vec::with_capacity(((2 * reach + 1).pow(2) * (2 * s_reach + 1)) as usize);
        for dd in -reach..=reach {
            let d = dec + f64::from(dd) * self.step;
            if !(-PI / 2.0..=PI / 2.0).contains(&d) {
                continue;
            }
            let cos_d = d.cos().max(1e-6);
            for dr in -reach..=reach {
                let r = ra + f64::from(dr) * self.step / cos_d;
                let (band, ra_bin, _) = self.key(r, d, 1.0);
                for ds in -s_reach..=s_reach {
                    let k = (band, ra_bin, key.2 + ds);
                    if !out.contains(&k) {
                        out.push(k);
                    }
                }
            }
        }
        out
    }

    /// Record a hypothesis centred at (`ra`, `dec`) radians with pixel `scale`
    /// (any positive unit, used only through its logarithm).
    pub(crate) fn add(&mut self, ra: f64, dec: f64, scale: f64, item: T) {
        let key = self.key(ra, dec, scale);
        self.buckets.entry(key).or_default().push(item);
    }

    /// Regions of agreement, strongest first, at most `limit` of them.
    ///
    /// Each bucket is scored by the votes in its 3×3×3 neighbourhood and
    /// represented by the neighbourhood's strongest bucket; after a region is
    /// taken, every bucket within two bins of it on the sky and one in scale is
    /// suppressed. Ties break on the bucket key, so the order is deterministic
    /// whatever the hash map's iteration order.
    pub(crate) fn regions(&self, limit: usize) -> Vec<Region<'_, T>> {
        let mut smoothed: Vec<(usize, Key, Key)> = self
            .buckets
            .keys()
            .map(|&key| {
                let mut sum = 0usize;
                let mut best: Option<(usize, Key)> = None;
                for n in self.neighbours(key, 1, 1) {
                    if let Some(v) = self.buckets.get(&n) {
                        sum += v.len();
                        let better =
                            best.is_none_or(|(bn, bk)| v.len() > bn || (v.len() == bn && n < bk));
                        if better {
                            best = Some((v.len(), n));
                        }
                    }
                }
                let (_, rep) = best.unwrap_or((0, key));
                (sum, key, rep)
            })
            .collect();
        smoothed.sort_by(|a, b| b.0.cmp(&a.0).then_with(|| a.1.cmp(&b.1)));

        let mut taken: HashSet<Key> = HashSet::new();
        let mut out = Vec::new();
        for (votes, key, rep) in smoothed {
            if out.len() >= limit {
                break;
            }
            if taken.contains(&key) || taken.contains(&rep) {
                continue;
            }
            taken.extend(self.neighbours(key, 2, 1));
            taken.extend(self.neighbours(rep, 1, 1));
            if let Some(members) = self.buckets.get(&rep) {
                out.push(Region {
                    votes,
                    members: members.as_slice(),
                });
            }
        }
        out
    }
}

/// Index of the member closest, in total, to all the others: a robust
/// representative of a bucket that holds a few stray hypotheses among agreeing
/// ones. Considers at most the first 64 members, so the cost stays bounded.
pub(crate) fn medoid<T>(members: &[T], pos: impl Fn(&T) -> (f64, f64)) -> usize {
    let m = &members[..members.len().min(64)];
    if m.len() <= 2 {
        return 0;
    }
    let pts: Vec<[f64; 3]> = m
        .iter()
        .map(|t| {
            let (ra, dec) = pos(t);
            [dec.cos() * ra.cos(), dec.cos() * ra.sin(), dec.sin()]
        })
        .collect();
    let mut best = (f64::INFINITY, 0usize);
    for (i, p) in pts.iter().enumerate() {
        let total: f64 = pts
            .iter()
            .map(|q| {
                let d = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
            })
            .sum();
        if total < best.0 {
            best = (total, i);
        }
    }
    best.1
}

#[cfg(test)]
mod tests {
    use super::*;

    fn deg(d: f64) -> f64 {
        d.to_radians()
    }

    #[test]
    fn a_split_peak_recombines_and_outranks_a_single_dense_bucket() {
        let mut v = SkyVotes::new(deg(0.1), 0.05);
        // Six votes straddling a bucket edge at RA 10.0°, against four in one bucket.
        for i in 0..6 {
            let ra = 10.0 + if i % 2 == 0 { -0.01 } else { 0.01 };
            v.add(deg(ra), deg(20.05), 2.0, i);
        }
        for i in 0..4 {
            v.add(deg(200.0), deg(-30.05), 2.0, 100 + i);
        }
        let r = v.regions(10);
        assert_eq!(r[0].votes, 6);
        assert!(r[0].members.iter().all(|&m| m < 100));
        assert_eq!(r[1].votes, 4);
        // Suppression leaves one region per peak.
        assert_eq!(r.len(), 2);
    }

    #[test]
    fn ra_bins_are_square_near_the_pole() {
        let mut v = SkyVotes::new(deg(0.1), 0.05);
        // At Dec 89.5° one degree of RA is 0.0087° of sky: a raw-RA bin would put
        // these in ten different buckets; a cos(Dec)-scaled one keeps them together.
        for i in 0..10 {
            v.add(deg(100.0 + f64::from(i) * 0.5), deg(89.55), 1.0, i);
        }
        assert!(v.len() <= 2, "{} buckets", v.len());
        assert_eq!(v.regions(5)[0].votes, 10);
    }

    #[test]
    fn ra_wraps_at_zero() {
        let mut v = SkyVotes::new(deg(0.1), 0.05);
        v.add(deg(359.99), deg(0.0), 1.0, 0);
        v.add(deg(0.01), deg(0.0), 1.0, 1);
        v.add(deg(0.03), deg(0.0), 1.0, 2);
        assert_eq!(v.regions(5)[0].votes, 3);
    }

    #[test]
    fn scale_separates_hypotheses_at_one_position() {
        let mut v = SkyVotes::new(deg(0.1), 0.05);
        for i in 0..3 {
            v.add(deg(50.0), deg(10.0), 1.0, i);
        }
        for i in 0..3 {
            v.add(deg(50.0), deg(10.0), 2.0, 10 + i);
        }
        let r = v.regions(5);
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].votes, 3);
    }

    #[test]
    fn the_medoid_ignores_a_stray_member() {
        let pts = [(0.0, 0.0), (deg(5.0), 0.0), (0.0001, 0.0), (0.0002, 0.0)];
        let i = medoid(&pts, |&p| p);
        assert!(i == 2 || i == 3 || i == 0);
        assert_ne!(i, 1);
    }
}
