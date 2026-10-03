//! Fixtures for the catalogue manager's tests: a temporary directory and a tiny
//! star database. The all-sky `.001` layout is simple enough to write here, and is
//! a real database format, so the registry and the index builder treat it exactly
//! as an installed one.

use core::sync::atomic::{AtomicUsize, Ordering};
use std::path::{Path, PathBuf};

/// A fresh directory under the system temp dir, removed on drop.
pub struct TempDir(PathBuf);

impl TempDir {
    /// Create a new empty directory; `tag` goes in its name.
    #[must_use]
    pub fn new(tag: &str) -> Self {
        static N: AtomicUsize = AtomicUsize::new(0);
        let n = N.fetch_add(1, Ordering::Relaxed);
        let p =
            std::env::temp_dir().join(format!("arcsec-cat-test-{}-{tag}-{n}", std::process::id()));
        drop(std::fs::remove_dir_all(&p));
        std::fs::create_dir_all(&p).expect("create temp dir");
        Self(p)
    }

    /// The directory.
    #[must_use]
    pub fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        drop(std::fs::remove_dir_all(&self.0));
    }
}

/// Write an all-sky `.001` database called `name` in `dir`: `n` stars spread
/// pseudo-randomly over the sky, magnitudes 1–7.5 (deterministic for a seed).
pub fn write_001_db(dir: &Path, name: &str, n: usize, seed: u64) {
    let mut x = seed.max(1);
    let mut next = move || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        (x >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut stars: Vec<(f32, f32, f32)> = (0..n)
        .map(|_| {
            let ra = next() * core::f64::consts::TAU;
            let dec = (2.0 * next() - 1.0).asin();
            let mag = 1.0 + 6.5 * next().sqrt();
            ((mag * 10.0) as f32, ra as f32, dec as f32)
        })
        .collect();
    stars.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut out = (stars.len() as u32).to_le_bytes().to_vec();
    for (m, r, d) in stars {
        out.extend_from_slice(&m.to_le_bytes());
        out.extend_from_slice(&r.to_le_bytes());
        out.extend_from_slice(&d.to_le_bytes());
    }
    std::fs::write(dir.join(format!("{name}_0101.001")), out).expect("write .001");
}
