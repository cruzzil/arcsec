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

// Test scaffolding, public only under the `test-support` feature: the library's
// API-hygiene lints are for the real API.
#![allow(
    missing_docs,
    clippy::must_use_candidate,
    clippy::missing_panics_doc,
    clippy::missing_errors_doc,
    clippy::doc_markdown,
    clippy::return_self_not_must_use
)]

use core::f64::consts::PI;
use core::sync::atomic::{AtomicUsize, Ordering};
use std::path::{Path, PathBuf};

use crate::catalog::areas::{area_and_boundaries_1476, filename_1476};
use crate::catalog::areas_290::{area_nr_290, filename_290};
use crate::types::{ImageBuffer, WcsSolution};

// ── Random numbers ─────────────────────────────────────────────────────────────

/// Deterministic xorshift64* generator, so no test can flake on its inputs.
pub struct Rng(u64);

impl Rng {
    /// A generator seeded with `seed` (zero is remapped; xorshift cannot leave it).
    pub fn new(seed: u64) -> Self {
        Self(if seed == 0 {
            0x9E37_79B9_7F4A_7C15
        } else {
            seed
        })
    }

    /// Next raw 64-bit value.
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// Uniform in `[0, 1)`.
    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform in `[lo, hi)`.
    pub fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.uniform()
    }

    /// Standard normal deviate (Box-Muller).
    pub fn gauss(&mut self) -> f64 {
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
pub struct TempDir(PathBuf);

impl TempDir {
    /// Create a new empty directory.
    pub fn new(tag: &str) -> Self {
        static COUNTER: AtomicUsize = AtomicUsize::new(0);
        let n = COUNTER.fetch_add(1, Ordering::Relaxed);
        let path =
            std::env::temp_dir().join(format!("arcsec-core-test-{}-{tag}-{n}", std::process::id()));
        drop(std::fs::remove_dir_all(&path));
        std::fs::create_dir_all(&path).expect("create temp dir");
        Self(path)
    }

    /// The directory.
    pub fn path(&self) -> &Path {
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
pub fn gnomonic(ra0: f64, dec0: f64, ra: f64, dec: f64) -> Option<(f64, f64)> {
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
pub fn inverse_gnomonic(ra0: f64, dec0: f64, xi: f64, eta: f64) -> (f64, f64) {
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
pub fn separation(ra1: f64, dec1: f64, ra2: f64, dec2: f64) -> f64 {
    let s_dec = ((dec2 - dec1) * 0.5).sin();
    let s_ra = ((ra2 - ra1) * 0.5).sin();
    let h = s_dec * s_dec + dec1.cos() * dec2.cos() * s_ra * s_ra;
    2.0 * h.sqrt().min(1.0).asin()
}

// ── Truth WCS ──────────────────────────────────────────────────────────────────

/// A TAN WCS with its reference pixel at the image centre, in the solver's pixel
/// convention: 0-based `(x, y)` with `data[y * width + x]`, FITS pixel = index + 1.
#[derive(Debug, Clone, Copy)]
pub struct TruthWcs {
    /// Reference RA (radians).
    pub ra0: f64,
    /// Reference Dec (radians).
    pub dec0: f64,
    /// CD matrix, degrees per pixel: `[cd1_1, cd1_2, cd2_1, cd2_2]`.
    pub cd: [f64; 4],
    /// Columns.
    pub width: usize,
    /// Rows.
    pub height: usize,
    /// Radial distortion about the frame centre, in the SIP sense (pixel → sky):
    /// the pixel `ρ` from the centre sees the sky the linear WCS puts at
    /// `ρ (1 + radial ρ²)`. Positive: the field's edges are squeezed (barrel);
    /// negative: stretched (pincushion); 0 is a pure TAN. A cubic in pixels, so
    /// SIP and arcsec's distortion model can represent it exactly.
    pub radial: f64,
}

impl TruthWcs {
    /// A WCS with `scale` arcsec/pixel rotated by `rot_deg`. `mirrored = false`
    /// gives the usual sky orientation (east left of north, det(CD) < 0);
    /// `true` flips the x axis, as a diagonal or mirror in the light path does.
    pub fn new(
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
            radial: 0.0,
        }
    }

    /// The same WCS with radial distortion of `corner_px` pixels at the frame's
    /// corners (barrel if positive).
    pub fn with_corner_distortion(mut self, corner_px: f64) -> Self {
        let (cx, cy) = self.centre();
        let r = cx.hypot(cy);
        self.radial = corner_px / (r * r * r);
        self
    }

    fn centre(&self) -> (f64, f64) {
        (
            (self.width as f64 - 1.0) * 0.5,
            (self.height as f64 - 1.0) * 0.5,
        )
    }

    /// Sky position of 0-based pixel `(x, y)`.
    pub fn pixel_to_sky(&self, x: f64, y: f64) -> (f64, f64) {
        let (cx, cy) = self.centre();
        let (dx, dy) = (x - cx, y - cy);
        let f = 1.0 + self.radial * (dx * dx + dy * dy);
        let (dx, dy) = (dx * f, dy * f);
        let xi = (self.cd[0] * dx + self.cd[1] * dy).to_radians();
        let eta = (self.cd[2] * dx + self.cd[3] * dy).to_radians();
        inverse_gnomonic(self.ra0, self.dec0, xi, eta)
    }

    /// 0-based pixel of a sky position (may lie outside the frame).
    pub fn sky_to_pixel(&self, ra: f64, dec: f64) -> Option<(f64, f64)> {
        let (xi, eta) = gnomonic(self.ra0, self.dec0, ra, dec)?;
        let (xi, eta) = (xi.to_degrees(), eta.to_degrees());
        let det = self.cd[0] * self.cd[3] - self.cd[1] * self.cd[2];
        let dx = (self.cd[3] * xi - self.cd[1] * eta) / det;
        let dy = (-self.cd[2] * xi + self.cd[0] * eta) / det;
        let (dx, dy) = self.distort(dx, dy)?;
        let (cx, cy) = self.centre();
        Some((cx + dx, cy + dy))
    }

    /// Linear pixel offset from the centre → the pixel offset that sees it: the
    /// root of `ρ + k ρ³ = ρ_lin` along the same direction, by Newton's method.
    /// `None` past the fold of a pincushion (`k < 0`), where `ρ + k ρ³` stops
    /// growing: nothing out there is imaged.
    fn distort(&self, lx: f64, ly: f64) -> Option<(f64, f64)> {
        let rl = lx.hypot(ly);
        let k = self.radial;
        if rl == 0.0 || k == 0.0 {
            return Some((lx, ly));
        }
        if k < 0.0 {
            let fold = (-1.0 / (3.0 * k)).sqrt();
            if rl >= 0.95 * (fold + k * fold * fold * fold) {
                return None;
            }
        }
        let mut r = rl;
        for _ in 0..60 {
            r -= (r + k * r * r * r - rl) / (1.0 + 3.0 * k * r * r);
        }
        Some((lx * r / rl, ly * r / rl))
    }

    /// Worst corner error (arcsec) of the best linear WCS over the frame: the
    /// least-squares linear fit to the truth on a 9 × 9 grid, as
    /// `scripts/fitslite.py`'s `best_linear` computes the benchmark's floor.
    pub fn linear_floor_arcsec(&self) -> f64 {
        const N: usize = 9;
        let (w, h) = (self.width as f64 - 1.0, self.height as f64 - 1.0);
        let (cx, cy) = self.centre();
        let mut img = Vec::new();
        let mut sky = Vec::new();
        for i in 0..N {
            for j in 0..N {
                let (x, y) = (w * i as f64 / (N - 1) as f64, h * j as f64 / (N - 1) as f64);
                let (ra, dec) = self.pixel_to_sky(x, y);
                let (xi, eta) = gnomonic(self.ra0, self.dec0, ra, dec).unwrap();
                img.push((x - cx, y - cy));
                sky.push((xi, eta));
            }
        }
        let lin = crate::math::lsq::fit_affine(&img, &sky).unwrap();
        [(0.0, 0.0), (w, 0.0), (0.0, h), (w, h)]
            .iter()
            .map(|&(x, y)| {
                let (u, v) = (x - cx, y - cy);
                let (xi, eta) = (lin.a * u + lin.b * v + lin.c, lin.d * u + lin.e * v + lin.f);
                let (ra_l, dec_l) = inverse_gnomonic(self.ra0, self.dec0, xi, eta);
                let (ra_t, dec_t) = self.pixel_to_sky(x, y);
                separation(ra_t, dec_t, ra_l, dec_l).to_degrees() * 3600.0
            })
            .fold(0.0, f64::max)
    }

    /// Worst disagreement (arcsec) between this WCS and a solved one, over the
    /// centre and the four corners — a centre-only check hides scale and rotation
    /// error, which is what `scripts/benchmark.py` checks too.
    pub fn max_error_arcsec(&self, sol: &WcsSolution) -> f64 {
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
pub fn solution_pixel_to_sky(sol: &WcsSolution, x: f64, y: f64) -> (f64, f64) {
    let dx = x + 1.0 - sol.crpix1;
    let dy = y + 1.0 - sol.crpix2;
    let xi = (sol.cd1_1 * dx + sol.cd1_2 * dy).to_radians();
    let eta = (sol.cd2_1 * dx + sol.cd2_2 * dy).to_radians();
    inverse_gnomonic(sol.ra0, sol.dec0, xi, eta)
}

// ── Synthetic sky and image ────────────────────────────────────────────────────

/// A catalogue star: RA/Dec in radians, magnitude.
#[derive(Debug, Clone, Copy)]
pub struct SkyStar {
    pub ra: f64,
    pub dec: f64,
    pub mag: f64,
}

/// The shape of a synthetic star field.
#[derive(Debug, Clone, Copy)]
pub struct SkySpec {
    /// Centre of the field (radians).
    pub ra0: f64,
    /// Centre of the field (radians).
    pub dec0: f64,
    /// Side of the square (tangent-plane degrees) the stars are spread over.
    pub side_deg: f64,
    /// Number of stars wanted.
    pub n: usize,
    /// No two stars closer than this (degrees). Blended pairs centroid badly, which
    /// would make the accuracy checks measure the fixture instead of the solver.
    pub min_sep_deg: f64,
    /// Brightest magnitude.
    pub mag_lo: f64,
    /// Faintest magnitude.
    pub mag_hi: f64,
}

/// Stars spread uniformly (with a minimum separation) over a square around the
/// centre, magnitudes weighted towards the faint end as real star counts are.
/// Returns fewer than `n` only if the square cannot hold that many.
pub fn random_sky(rng: &mut Rng, spec: &SkySpec) -> Vec<SkyStar> {
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
pub fn render(
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
pub fn area_file_bytes(stars: &[SkyStar], record_size: usize) -> Vec<u8> {
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
pub fn write_1476_db(dir: &Path, name: &str, stars: &[SkyStar]) {
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
pub fn write_290_db(dir: &Path, name: &str, stars: &[SkyStar]) {
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
pub fn file_001_bytes(stars: &[SkyStar]) -> Vec<u8> {
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
pub fn write_001_db(dir: &Path, name: &str, stars: &[SkyStar]) {
    std::fs::write(dir.join(format!("{name}_0101.001")), file_001_bytes(stars))
        .expect("write .001 file");
}

// ── Astrometry.net index fixtures ──────────────────────────────────────────────

/// Lower end of the astrometry.net code range, `0.5 - √2/2`.
pub const ANET_CODE_LO: f64 = 0.5 - core::f64::consts::FRAC_1_SQRT_2;
/// Upper end of the astrometry.net code range, `0.5 + √2/2`.
pub const ANET_CODE_HI: f64 = 0.5 + core::f64::consts::FRAC_1_SQRT_2;

/// An index as astrometry.net's builder would lay it out, before FITS encoding.
#[derive(Debug, Clone)]
pub struct RawIndex {
    /// Stars per quad: 3 or 4.
    pub dim_quads: usize,
    /// Every index star, `(ra, dec)` radians.
    pub stars: Vec<(f64, f64)>,
    /// Per quad, `dim_quads` star indices in canonical order (A, B, C[, D]).
    pub quads: Vec<Vec<u32>>,
    /// Per quad, its code (`2 * (dim_quads - 2)` used slots).
    pub codes: Vec<[f64; 4]>,
    /// Smallest A-B separation (radians).
    pub scale_lo: f64,
    /// Largest A-B separation (radians).
    pub scale_hi: f64,
}

fn unit(ra: f64, dec: f64) -> [f64; 3] {
    [dec.cos() * ra.cos(), dec.cos() * ra.sin(), dec.sin()]
}

/// The canonical astrometry.net code of a star group, computed on the sky.
///
/// Written from astrometry.net's own recipe, independently of the image-side code
/// in `pipeline::blind`: A and B are the most separated pair; every star is
/// projected onto the tangent plane at the A-B midpoint (x east, y north); the
/// frame is rotated so A lands on (0, 0) and B on (1, 1); then the invariants are
/// enforced — mean C/D x-code ≤ 0.5 (else swap A and B), and C/D sorted by x-code.
///
/// Returns the star indices in canonical order and the code, or `None` for a
/// degenerate group or one whose other stars fall outside the circle on A-B as
/// diameter, which astrometry.net does not build.
pub fn anet_code(sky: &[(f64, f64)], group: &[u32]) -> Option<(Vec<u32>, [f64; 4])> {
    let n = group.len();
    let pos = |i: u32| sky[i as usize];
    // Most separated pair becomes A-B.
    let mut best = (0usize, 1usize, -1.0f64);
    for i in 0..n {
        for j in (i + 1)..n {
            let (a, b) = (pos(group[i]), pos(group[j]));
            let d = separation(a.0, a.1, b.0, b.1);
            if d > best.2 {
                best = (i, j, d);
            }
        }
    }
    let (ia, ib, _) = best;
    let mut order = vec![group[ia], group[ib]];
    order.extend((0..n).filter(|&k| k != ia && k != ib).map(|k| group[k]));

    let (ua, ub) = (
        unit(pos(order[0]).0, pos(order[0]).1),
        unit(pos(order[1]).0, pos(order[1]).1),
    );
    let m = [ua[0] + ub[0], ua[1] + ub[1], ua[2] + ub[2]];
    let norm = (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt();
    if norm < 1e-12 {
        return None;
    }
    let ra_m = m[1].atan2(m[0]);
    let dec_m = (m[2] / norm).asin();
    let xy: Vec<(f64, f64)> = order
        .iter()
        .map(|&i| gnomonic(ra_m, dec_m, pos(i).0, pos(i).1))
        .collect::<Option<_>>()?;

    let (abx, aby) = (xy[1].0 - xy[0].0, xy[1].1 - xy[0].1);
    let scale = abx * abx + aby * aby;
    if scale < 1e-24 {
        return None;
    }
    let (cos_t, sin_t) = ((aby + abx) / scale, (aby - abx) / scale);
    let mut codes: Vec<(f64, f64, u32)> = xy[2..]
        .iter()
        .zip(&order[2..])
        .map(|(&(x, y), &id)| {
            let (dx, dy) = (x - xy[0].0, y - xy[0].1);
            (dx * cos_t + dy * sin_t, -dx * sin_t + dy * cos_t, id)
        })
        .collect();

    // The builder only accepts groups whose other stars lie inside the circle on
    // A-B as diameter, which is what bounds codes to [0.5 - √2/2, 0.5 + √2/2].
    if codes
        .iter()
        .any(|c| (c.0 - 0.5).powi(2) + (c.1 - 0.5).powi(2) > 0.5)
    {
        return None;
    }
    let mean_x = codes.iter().map(|c| c.0).sum::<f64>() / codes.len() as f64;
    if mean_x > 0.5 {
        order.swap(0, 1);
        for c in &mut codes {
            c.0 = 1.0 - c.0;
            c.1 = 1.0 - c.1;
        }
    }
    codes.sort_by(|p, q| p.0.total_cmp(&q.0).then(p.1.total_cmp(&q.1)));

    let mut code = [0.0; 4];
    for (k, c) in codes.iter().enumerate() {
        code[2 * k] = c.0;
        code[2 * k + 1] = c.1;
        order[2 + k] = c.2;
    }
    Some((order, code))
}

impl RawIndex {
    /// An index over `sky`, with one entry per star group in `groups`.
    pub fn build(dim_quads: usize, sky: &[(f64, f64)], groups: &[Vec<u32>]) -> Self {
        let mut quads = Vec::new();
        let mut codes = Vec::new();
        let (mut lo, mut hi) = (f64::INFINITY, 0.0f64);
        for g in groups {
            assert_eq!(g.len(), dim_quads);
            let Some((order, code)) = anet_code(sky, g) else {
                continue;
            };
            let (a, b) = (sky[order[0] as usize], sky[order[1] as usize]);
            let d = separation(a.0, a.1, b.0, b.1);
            lo = lo.min(d);
            hi = hi.max(d);
            quads.push(order);
            codes.push(code);
        }
        Self {
            dim_quads,
            stars: sky.to_vec(),
            quads,
            codes,
            scale_lo: lo,
            scale_hi: hi,
        }
    }

    /// The in-memory [`crate::catalog::AnetIndex`] that loading this file would give,
    /// built without going through FITS.
    pub fn to_index(&self) -> crate::catalog::AnetIndex {
        use crate::catalog::anet::{AnetIndex, AnetIndexEntry, AnetStar};
        let stars: Vec<AnetStar> = self
            .stars
            .iter()
            .map(|&(ra, dec)| AnetStar { ra, dec })
            .collect();
        let mut entries: Vec<AnetIndexEntry> = self
            .quads
            .iter()
            .zip(&self.codes)
            .map(|(q, &code)| {
                let mut star_ra = [0.0; 4];
                let mut star_dec = [0.0; 4];
                let mut c = [0.0f64; 3];
                for (k, &i) in q.iter().enumerate() {
                    let (ra, dec) = self.stars[i as usize];
                    star_ra[k] = ra;
                    star_dec[k] = dec;
                    let u = unit(ra, dec);
                    c = [c[0] + u[0], c[1] + u[1], c[2] + u[2]];
                }
                AnetIndexEntry {
                    code,
                    n_stars: q.len(),
                    star_ra,
                    star_dec,
                    center_ra: c[1].atan2(c[0]).rem_euclid(2.0 * PI),
                    center_dec: c[2].atan2(c[0].hypot(c[1])),
                }
            })
            .collect();
        entries.sort_by(|a, b| a.code[0].total_cmp(&b.code[0]));
        let codes = entries.iter().map(|e| e.code.map(|v| v as f32)).collect();
        AnetIndex {
            entries,
            codes,
            stars,
            scale_lo: self.scale_lo,
            scale_hi: self.scale_hi,
            dim_quads: self.dim_quads,
        }
    }

    /// Encode as an astrometry.net index FITS file.
    pub fn fits_bytes(&self) -> Vec<u8> {
        let n_dims = 2 * (self.dim_quads - 2);
        let code_scale = 65535.0 / (ANET_CODE_HI - ANET_CODE_LO);
        let mut fits = FitsWriter::default();
        fits.primary(&[
            ("DIMQUADS", self.dim_quads.to_string()),
            ("NQUADS", self.quads.len().to_string()),
            ("NSTARS", self.stars.len().to_string()),
            ("SCALE_U", format!("{:.15E}", self.scale_hi)),
            ("SCALE_L", format!("{:.15E}", self.scale_lo)),
        ]);

        let quad_bytes: Vec<u8> = self
            .quads
            .iter()
            .flatten()
            .flat_map(|i| i.to_le_bytes())
            .collect();
        fits.table("quads", 4 * self.dim_quads, &quad_bytes);

        // Range table: n lo values, n hi values, then the scale.
        let mut range = Vec::new();
        for _ in 0..n_dims {
            range.extend_from_slice(&ANET_CODE_LO.to_le_bytes());
        }
        for _ in 0..n_dims {
            range.extend_from_slice(&ANET_CODE_HI.to_le_bytes());
        }
        range.extend_from_slice(&code_scale.to_le_bytes());
        fits.table("kdtree_range_codes", 8, &range);

        let code_bytes: Vec<u8> = self
            .codes
            .iter()
            .flat_map(|c| {
                c[..n_dims].iter().flat_map(|&v| {
                    (((v - ANET_CODE_LO) * code_scale)
                        .round()
                        .clamp(0.0, 65535.0) as u16)
                        .to_le_bytes()
                })
            })
            .collect();
        fits.table("kdtree_data_codes", 2 * n_dims, &code_bytes);

        let star_bytes: Vec<u8> = self
            .stars
            .iter()
            .flat_map(|&(ra, dec)| {
                unit(ra, dec).map(|v| {
                    (((v + 1.0) * (f64::from(u32::MAX) / 2.0)).round() as u32).to_le_bytes()
                })
            })
            .flatten()
            .collect();
        fits.table("kdtree_data_stars", 12, &star_bytes);
        // Real indexes carry more tables after the stars (sweep, magnitudes); the
        // loader cannot reach an index's final HDU (see anet.rs tests), so end with
        // one it does not need, as they do.
        fits.table("sweep", 1, &vec![0u8; self.stars.len()]);
        fits.bytes
    }
}

/// Minimal FITS writer: a primary header and single-column `nA` binary tables,
/// which is all an astrometry.net index uses.
#[derive(Default)]
pub struct FitsWriter {
    /// The file so far.
    pub bytes: Vec<u8>,
}

impl FitsWriter {
    fn card(&mut self, text: &str) {
        let mut c = text.as_bytes().to_vec();
        c.resize(80, b' ');
        self.bytes.extend_from_slice(&c);
    }

    fn value(&mut self, key: &str, value: &str) {
        self.card(&format!("{key:<8}= {value:>20}"));
    }

    fn end_header(&mut self) {
        self.card("END");
        let pad = self.bytes.len().next_multiple_of(2880);
        self.bytes.resize(pad, b' ');
    }

    /// Primary HDU with no data and the given integer/float keywords.
    pub fn primary(&mut self, keys: &[(&str, String)]) {
        self.value("SIMPLE", "T");
        self.value("BITPIX", "8");
        self.value("NAXIS", "0");
        self.value("EXTEND", "T");
        for (k, v) in keys {
            self.value(k, v);
        }
        self.end_header();
    }

    /// A binary table with one `row_bytes`-wide raw-byte column named `name`.
    pub fn table(&mut self, name: &str, row_bytes: usize, data: &[u8]) {
        assert_eq!(data.len() % row_bytes, 0);
        self.card("XTENSION= 'BINTABLE'");
        self.value("BITPIX", "8");
        self.value("NAXIS", "2");
        self.value("NAXIS1", &row_bytes.to_string());
        self.value("NAXIS2", &(data.len() / row_bytes).to_string());
        self.value("PCOUNT", "0");
        self.value("GCOUNT", "1");
        self.value("TFIELDS", "1");
        self.card(&format!("TTYPE1  = '{name}'"));
        self.card(&format!("TFORM1  = '{row_bytes}A'"));
        self.end_header();
        self.bytes.extend_from_slice(data);
        let pad = self.bytes.len().next_multiple_of(2880);
        self.bytes.resize(pad, 0);
    }
}

// ── A ready-made field, for other crates' end-to-end tests ─────────────────────

/// A synthetic field: an image of a random sky through `truth`, with that sky
/// (six fields wide, so offset hints still find their stars) written to `dir` as
/// a `.1476` database called `db_name`. About `n_in_frame` stars fall in the frame.
pub fn synthetic_field(
    dir: &Path,
    db_name: &str,
    truth: &TruthWcs,
    n_in_frame: usize,
    seed: u64,
) -> ImageBuffer {
    let mut rng = Rng::new(seed);
    let scale_deg = truth.cd[1].hypot(truth.cd[3]);
    let (w_deg, h_deg) = (
        truth.width as f64 * scale_deg,
        truth.height as f64 * scale_deg,
    );
    let side = 6.0 * w_deg.max(h_deg);
    let sky = random_sky(
        &mut rng,
        &SkySpec {
            ra0: truth.ra0,
            dec0: truth.dec0,
            side_deg: side,
            n: (n_in_frame as f64 * side * side / (w_deg * h_deg)) as usize,
            min_sep_deg: 12.0 * scale_deg,
            mag_lo: 10.0,
            mag_hi: 14.5,
        },
    );
    let sigma = (1.3 * 5.0 / (scale_deg * 3600.0)).max(1.3);
    let img = render(truth, &sky, sigma, 1000.0, 8.0, 30_000.0, &mut rng);
    write_1476_db(dir, db_name, &sky);
    img
}

/// The bytes of a FITS file holding `img` as 32-bit floats (BITPIX -32), with
/// extra header `cards` (`("FOCALLEN", "206.265")`, ...). Row 0 of `img` is the
/// first row of the file.
pub fn fits_f32_bytes(img: &ImageBuffer, cards: &[(&str, &str)]) -> Vec<u8> {
    let mut hdr = String::new();
    let mut card = |s: String| hdr.push_str(&format!("{s:<80}"));
    card(format!("{:<8}= {:>20}", "SIMPLE", "T"));
    card(format!("{:<8}= {:>20}", "BITPIX", "-32"));
    card(format!("{:<8}= {:>20}", "NAXIS", "2"));
    card(format!("{:<8}= {:>20}", "NAXIS1", img.width));
    card(format!("{:<8}= {:>20}", "NAXIS2", img.height));
    for (k, v) in cards {
        card(format!("{k:<8}= {v:>20}"));
    }
    card("END".to_string());
    let mut out = hdr.into_bytes();
    out.resize(out.len().next_multiple_of(2880), b' ');
    for v in &img.data {
        out.extend_from_slice(&v.to_be_bytes());
    }
    out.resize(out.len().next_multiple_of(2880), 0);
    out
}
