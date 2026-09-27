//! Fixtures shared by the in-file unit tests (compiled only under `cfg(test)`).
//!
//! A plate solve cannot be exercised without star catalogues and images, and none
//! ship with the repository, so this module makes both: a synthetic sky, a truth
//! WCS that projects it onto a detector, a renderer that draws Gaussian stars with
//! noise, and writers for ASTAP's `.1476`, `.290` and `.001` database formats.
//!
//! The gnomonic projection here is written out from the textbook formulas rather
//! than borrowed from [`crate::math::coords`], so a solve checked against it is
//! checked against an independent implementation.

use core::f64::consts::PI;
use core::sync::atomic::{AtomicUsize, Ordering};
use std::path::{Path, PathBuf};

use crate::catalog::areas::{area_and_boundaries_1476, filename_1476};
use crate::catalog::areas_290::{area_nr_290, filename_290};
use crate::types::{ImageBuffer, WcsSolution};

// ── Random numbers ─────────────────────────────────────────────────────────────

/// Deterministic xorshift64* generator, so no test can flake on its inputs.
pub(crate) struct Rng(u64);

impl Rng {
    /// A generator seeded with `seed` (zero is remapped; xorshift cannot leave it).
    pub(crate) fn new(seed: u64) -> Self {
        Self(if seed == 0 {
            0x9E37_79B9_7F4A_7C15
        } else {
            seed
        })
    }

    /// Next raw 64-bit value.
    pub(crate) fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// Uniform in `[0, 1)`.
    pub(crate) fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform in `[lo, hi)`.
    pub(crate) fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.uniform()
    }

    /// Standard normal deviate (Box-Muller).
    pub(crate) fn gauss(&mut self) -> f64 {
        let u1 = self.uniform().max(1e-300);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
    }
}

// ── Temporary directories ──────────────────────────────────────────────────────

/// A fresh directory under the system temp dir, removed again on drop.
///
/// Named by process id, a caller tag and a counter, so parallel tests (and
/// parallel test binaries) never share one.
pub(crate) struct TempDir(PathBuf);

impl TempDir {
    /// Create a new empty directory.
    pub(crate) fn new(tag: &str) -> Self {
        static COUNTER: AtomicUsize = AtomicUsize::new(0);
        let n = COUNTER.fetch_add(1, Ordering::Relaxed);
        let path =
            std::env::temp_dir().join(format!("arcsec-core-test-{}-{tag}-{n}", std::process::id()));
        drop(std::fs::remove_dir_all(&path));
        std::fs::create_dir_all(&path).expect("create temp dir");
        Self(path)
    }

    /// The directory.
    pub(crate) fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        drop(std::fs::remove_dir_all(&self.0));
    }
}

// ── Independent gnomonic projection ────────────────────────────────────────────

/// Sky → standard coordinates (radians; xi east, eta north). `None` behind the
/// tangent point.
pub(crate) fn gnomonic(ra0: f64, dec0: f64, ra: f64, dec: f64) -> Option<(f64, f64)> {
    let dra = ra - ra0;
    let cos_c = dec0.sin() * dec.sin() + dec0.cos() * dec.cos() * dra.cos();
    if cos_c <= 1e-9 {
        return None;
    }
    let xi = dec.cos() * dra.sin() / cos_c;
    let eta = (dec0.cos() * dec.sin() - dec0.sin() * dec.cos() * dra.cos()) / cos_c;
    Some((xi, eta))
}

/// Standard coordinates (radians; xi east, eta north) → sky, RA in `[0, 2π)`.
pub(crate) fn inverse_gnomonic(ra0: f64, dec0: f64, xi: f64, eta: f64) -> (f64, f64) {
    let rho = xi.hypot(eta);
    if rho < 1e-300 {
        return (ra0.rem_euclid(2.0 * PI), dec0);
    }
    let c = rho.atan();
    let (sin_c, cos_c) = c.sin_cos();
    let dec = (cos_c * dec0.sin() + eta * sin_c * dec0.cos() / rho).asin();
    let ra = ra0 + (xi * sin_c).atan2(rho * dec0.cos() * cos_c - eta * dec0.sin() * sin_c);
    (ra.rem_euclid(2.0 * PI), dec)
}

/// Great-circle separation (radians), by the haversine formula (well conditioned
/// for the tiny separations the accuracy checks deal in).
pub(crate) fn separation(ra1: f64, dec1: f64, ra2: f64, dec2: f64) -> f64 {
    let s_dec = ((dec2 - dec1) * 0.5).sin();
    let s_ra = ((ra2 - ra1) * 0.5).sin();
    let h = s_dec * s_dec + dec1.cos() * dec2.cos() * s_ra * s_ra;
    2.0 * h.sqrt().min(1.0).asin()
}

// ── Truth WCS ──────────────────────────────────────────────────────────────────

/// A TAN WCS with its reference pixel at the image centre, in the solver's pixel
/// convention: 0-based `(x, y)` with `data[y * width + x]`, FITS pixel = index + 1.
#[derive(Debug, Clone, Copy)]
pub(crate) struct TruthWcs {
    /// Reference RA (radians).
    pub(crate) ra0: f64,
    /// Reference Dec (radians).
    pub(crate) dec0: f64,
    /// CD matrix, degrees per pixel: `[cd1_1, cd1_2, cd2_1, cd2_2]`.
    pub(crate) cd: [f64; 4],
    /// Columns.
    pub(crate) width: usize,
    /// Rows.
    pub(crate) height: usize,
}

impl TruthWcs {
    /// A WCS with `scale` arcsec/pixel rotated by `rot_deg`. `mirrored = false`
    /// gives the usual sky orientation (east left of north, det(CD) < 0);
    /// `true` flips the x axis, as a diagonal or mirror in the light path does.
    pub(crate) fn new(
        ra0: f64,
        dec0: f64,
        scale_arcsec: f64,
        rot_deg: f64,
        mirrored: bool,
        width: usize,
        height: usize,
    ) -> Self {
        let s = scale_arcsec / 3600.0;
        let (sin_r, cos_r) = rot_deg.to_radians().sin_cos();
        let px = if mirrored { 1.0 } else { -1.0 };
        Self {
            ra0,
            dec0,
            cd: [px * s * cos_r, s * sin_r, -px * s * sin_r, s * cos_r],
            width,
            height,
        }
    }

    fn centre(&self) -> (f64, f64) {
        (
            (self.width as f64 - 1.0) * 0.5,
            (self.height as f64 - 1.0) * 0.5,
        )
    }

    /// Sky position of 0-based pixel `(x, y)`.
    pub(crate) fn pixel_to_sky(&self, x: f64, y: f64) -> (f64, f64) {
        let (cx, cy) = self.centre();
        let (dx, dy) = (x - cx, y - cy);
        let xi = (self.cd[0] * dx + self.cd[1] * dy).to_radians();
        let eta = (self.cd[2] * dx + self.cd[3] * dy).to_radians();
        inverse_gnomonic(self.ra0, self.dec0, xi, eta)
    }

    /// 0-based pixel of a sky position (may lie outside the frame).
    pub(crate) fn sky_to_pixel(&self, ra: f64, dec: f64) -> Option<(f64, f64)> {
        let (xi, eta) = gnomonic(self.ra0, self.dec0, ra, dec)?;
        let (xi, eta) = (xi.to_degrees(), eta.to_degrees());
        let det = self.cd[0] * self.cd[3] - self.cd[1] * self.cd[2];
        let dx = (self.cd[3] * xi - self.cd[1] * eta) / det;
        let dy = (-self.cd[2] * xi + self.cd[0] * eta) / det;
        let (cx, cy) = self.centre();
        Some((cx + dx, cy + dy))
    }

    /// Worst disagreement (arcsec) between this WCS and a solved one, over the
    /// centre and the four corners — a centre-only check hides scale and rotation
    /// error, which is what `scripts/benchmark.py` checks too.
    pub(crate) fn max_error_arcsec(&self, sol: &WcsSolution) -> f64 {
        let (w, h) = (self.width as f64 - 1.0, self.height as f64 - 1.0);
        [(w * 0.5, h * 0.5), (0.0, 0.0), (w, 0.0), (0.0, h), (w, h)]
            .iter()
            .map(|&(x, y)| {
                let (ra_t, dec_t) = self.pixel_to_sky(x, y);
                let (ra_s, dec_s) = solution_pixel_to_sky(sol, x, y);
                separation(ra_t, dec_t, ra_s, dec_s).to_degrees() * 3600.0
            })
            .fold(0.0, f64::max)
    }
}

/// Evaluate a solver [`WcsSolution`] at 0-based pixel `(x, y)` using the standard
/// FITS TAN definition of its CRVAL/CRPIX/CD keywords.
pub(crate) fn solution_pixel_to_sky(sol: &WcsSolution, x: f64, y: f64) -> (f64, f64) {
    let dx = x + 1.0 - sol.crpix1;
    let dy = y + 1.0 - sol.crpix2;
    let xi = (sol.cd1_1 * dx + sol.cd1_2 * dy).to_radians();
    let eta = (sol.cd2_1 * dx + sol.cd2_2 * dy).to_radians();
    inverse_gnomonic(sol.ra0, sol.dec0, xi, eta)
}

// ── Synthetic sky and image ────────────────────────────────────────────────────

/// A catalogue star: RA/Dec in radians, magnitude.
#[derive(Debug, Clone, Copy)]
pub(crate) struct SkyStar {
    pub(crate) ra: f64,
    pub(crate) dec: f64,
    pub(crate) mag: f64,
}

/// The shape of a synthetic star field.
#[derive(Debug, Clone, Copy)]
pub(crate) struct SkySpec {
    /// Centre of the field (radians).
    pub(crate) ra0: f64,
    /// Centre of the field (radians).
    pub(crate) dec0: f64,
    /// Side of the square (tangent-plane degrees) the stars are spread over.
    pub(crate) side_deg: f64,
    /// Number of stars wanted.
    pub(crate) n: usize,
    /// No two stars closer than this (degrees). Blended pairs centroid badly, which
    /// would make the accuracy checks measure the fixture instead of the solver.
    pub(crate) min_sep_deg: f64,
    /// Brightest magnitude.
    pub(crate) mag_lo: f64,
    /// Faintest magnitude.
    pub(crate) mag_hi: f64,
}

/// Stars spread uniformly (with a minimum separation) over a square around the
/// centre, magnitudes weighted towards the faint end as real star counts are.
/// Returns fewer than `n` only if the square cannot hold that many.
pub(crate) fn random_sky(rng: &mut Rng, spec: &SkySpec) -> Vec<SkyStar> {
    let half = (spec.side_deg * 0.5).to_radians();
    let sep = spec.min_sep_deg.to_radians();
    // Dart throwing, bucketed on a grid of `sep`-sized cells.
    let cell = |v: f64| (v / sep.max(1e-12)).floor() as i64;
    let mut grid: std::collections::HashMap<(i64, i64), Vec<(f64, f64)>> =
        std::collections::HashMap::new();
    let mut out = Vec::with_capacity(spec.n);
    let mut attempts = 0usize;
    while out.len() < spec.n && attempts < 50 * spec.n.max(1) {
        attempts += 1;
        let xi = rng.range(-half, half);
        let eta = rng.range(-half, half);
        let (gx, gy) = (cell(xi), cell(eta));
        let crowded = (-1..=1).any(|ox| {
            (-1..=1).any(|oy| {
                grid.get(&(gx + ox, gy + oy))
                    .is_some_and(|b| b.iter().any(|&(x, y)| (x - xi).hypot(y - eta) < sep))
            })
        });
        if crowded {
            continue;
        }
        grid.entry((gx, gy)).or_default().push((xi, eta));
        let (ra, dec) = inverse_gnomonic(spec.ra0, spec.dec0, xi, eta);
        let mag = spec.mag_lo + (spec.mag_hi - spec.mag_lo) * rng.uniform().powf(0.4);
        out.push(SkyStar { ra, dec, mag });
    }
    out
}

/// Draw `stars` through `wcs`: Gaussian PSFs of width `sigma_px` on a flat
/// `background` with Gaussian read noise `noise`. A magnitude-10 star peaks at
/// `peak10` counts above the background.
pub(crate) fn render(
    wcs: &TruthWcs,
    stars: &[SkyStar],
    sigma_px: f64,
    background: f64,
    noise: f64,
    peak10: f64,
    rng: &mut Rng,
) -> ImageBuffer {
    let (w, h) = (wcs.width, wcs.height);
    let mut buf = vec![0.0f64; w * h];
    let r = (5.0 * sigma_px).ceil() as i64;
    for s in stars {
        let Some((x, y)) = wcs.sky_to_pixel(s.ra, s.dec) else {
            continue;
        };
        if x < -10.0 || y < -10.0 || x > w as f64 + 10.0 || y > h as f64 + 10.0 {
            continue;
        }
        let peak = peak10 * 10f64.powf(-0.4 * (s.mag - 10.0));
        let (xi, yi) = (x.round() as i64, y.round() as i64);
        for py in (yi - r).max(0)..=(yi + r).min(h as i64 - 1) {
            for px in (xi - r).max(0)..=(xi + r).min(w as i64 - 1) {
                let d2 = (px as f64 - x).powi(2) + (py as f64 - y).powi(2);
                buf[py as usize * w + px as usize] +=
                    peak * (-d2 / (2.0 * sigma_px * sigma_px)).exp();
            }
        }
    }
    let data = buf
        .into_iter()
        .map(|v| (v + background + noise * rng.gauss()).clamp(0.0, 65_000.0) as f32)
        .collect();
    ImageBuffer {
        data,
        width: w,
        height: h,
    }
}

// ── ASTAP database writers ─────────────────────────────────────────────────────

/// Encode one star as a packed 5- or 6-byte area-file record, returning the record
/// and the header-record key (`dec9`, magnitude byte) it needs in front of it.
fn encode_record(s: &SkyStar, record_size: usize) -> ((u8, u8), Vec<u8>) {
    let mut ra_raw = (s.ra.rem_euclid(2.0 * PI) / (2.0 * PI) * 16_777_215.0).round() as u32;
    if ra_raw >= 0xFF_FF_FF {
        ra_raw = 0; // 2π is 0, and 0xFFFFFF is the header-record sentinel
    }
    let dec_raw = (s.dec / (PI * 0.5) * 8_388_607.0).round() as i32;
    let dec9 = dec_raw >> 16; // arithmetic shift: -128..=127
    let key = (
        (dec9 + 128) as u8,
        (s.mag * 10.0 + 16.0).round().clamp(0.0, 255.0) as u8,
    );
    let mut rec = vec![
        (ra_raw & 0xFF) as u8,
        ((ra_raw >> 8) & 0xFF) as u8,
        ((ra_raw >> 16) & 0xFF) as u8,
        (dec_raw & 0xFF) as u8,
        ((dec_raw >> 8) & 0xFF) as u8,
    ];
    rec.resize(record_size, 0);
    (key, rec)
}

/// The bytes of one `.1476`/`.290` area file: the 110-byte header, then the stars
/// brightest first, with a header record wherever the magnitude or the high
/// declination byte changes — the layout ASTAP writes.
pub(crate) fn area_file_bytes(stars: &[SkyStar], record_size: usize) -> Vec<u8> {
    let mut sorted = stars.to_vec();
    sorted.sort_by(|a, b| a.mag.total_cmp(&b.mag));
    let mut out = vec![b' '; 110];
    out[109] = record_size as u8;
    let mut current: Option<(u8, u8)> = None;
    for s in &sorted {
        let (key, rec) = encode_record(s, record_size);
        if current != Some(key) {
            let mut hdr = vec![0xFF, 0xFF, 0xFF, key.0, key.1];
            hdr.resize(record_size, 0);
            out.extend_from_slice(&hdr);
            current = Some(key);
        }
        out.extend_from_slice(&rec);
    }
    out
}

/// Write `stars` as a `.1476` database called `name` in `dir`: one file per
/// occupied area, plus the south-pole `0101` tile that layout detection probes.
pub(crate) fn write_1476_db(dir: &Path, name: &str, stars: &[SkyStar]) {
    let mut by_area: alloc::collections::BTreeMap<usize, Vec<SkyStar>> =
        alloc::collections::BTreeMap::new();
    by_area.entry(1).or_default();
    for s in stars {
        let area = area_and_boundaries_1476(s.ra, s.dec).area_nr;
        by_area.entry(area).or_default().push(*s);
    }
    for (area, list) in by_area {
        let path = dir.join(format!("{name}_{}", filename_1476(area)));
        std::fs::write(path, area_file_bytes(&list, 5)).expect("write area file");
    }
}

/// Write `stars` as a `.290` database called `name` in `dir`.
pub(crate) fn write_290_db(dir: &Path, name: &str, stars: &[SkyStar]) {
    let mut by_area: alloc::collections::BTreeMap<usize, Vec<SkyStar>> =
        alloc::collections::BTreeMap::new();
    by_area.entry(1).or_default();
    for s in stars {
        by_area
            .entry(area_nr_290(s.ra, s.dec))
            .or_default()
            .push(*s);
    }
    for (area, list) in by_area {
        let path = dir.join(format!("{name}_{}", filename_290(area)));
        std::fs::write(path, area_file_bytes(&list, 6)).expect("write area file");
    }
}

/// The bytes of a `.001` all-sky file: a u32 count, then f32 triples brightest first.
pub(crate) fn file_001_bytes(stars: &[SkyStar]) -> Vec<u8> {
    let mut sorted = stars.to_vec();
    sorted.sort_by(|a, b| a.mag.total_cmp(&b.mag));
    let mut out = (sorted.len() as u32).to_le_bytes().to_vec();
    for s in &sorted {
        out.extend_from_slice(&((s.mag * 10.0) as f32).to_le_bytes());
        out.extend_from_slice(&(s.ra as f32).to_le_bytes());
        out.extend_from_slice(&(s.dec as f32).to_le_bytes());
    }
    out
}

/// Write `stars` as a `.001` database called `name` in `dir`.
pub(crate) fn write_001_db(dir: &Path, name: &str, stars: &[SkyStar]) {
    std::fs::write(dir.join(format!("{name}_0101.001")), file_001_bytes(stars))
        .expect("write .001 file");
}
