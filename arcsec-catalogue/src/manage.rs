//! Choosing, removing and checking installed catalogues.

use core::fmt;
use std::fs;
use std::path::{Path, PathBuf};

use arcsec_core::ArcsecError;

use crate::error::{Error, Result};
use crate::index::{Existing, Freshness, freshness, index_files};
use crate::registry::{
    Archive, Entry, Purpose, REGISTRY, astap_file_count, expected_file_count, files_of,
    installed_size, is_installed,
};

/// The catalogue to suggest for `purpose` and a field `fov_deg` degrees across: the
/// smallest download whose published range covers the field, so the advice does
/// not push a gigabyte on someone who does not need it. `None` if none covers it.
#[must_use]
pub fn recommend(purpose: Purpose, fov_deg: f64) -> Option<&'static Entry> {
    REGISTRY
        .iter()
        .filter(|e| e.purpose == purpose)
        .filter(|e| matches!(e.fov, Some((lo, hi)) if fov_deg >= lo && fov_deg <= hi))
        .min_by_key(|e| e.bytes)
}

/// Disk an install of `e` needs at its peak: an archive and its extracted files
/// side by side (assumed no smaller than the archive), or the loose files.
#[must_use]
pub fn download_disk(e: &Entry) -> u64 {
    match e.archive {
        Archive::Loose => e.bytes,
        Archive::Zip | Archive::Deb => e.bytes.saturating_mul(2),
    }
}

/// What [`remove`] removed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Removed {
    /// Files deleted.
    pub files: usize,
    /// Bytes freed.
    pub bytes: u64,
}

/// Delete every file of catalogue `e` from `dir`. Blind indexes built from it are
/// not touched; see [`crate::index::indexes_removed_with`].
///
/// The directory may be shared with ASTAP, in which case these are ASTAP's files
/// too: [`crate::registry::files_of`] lists exactly what would go.
///
/// # Errors
///
/// [`Error::Io`] for the first file that cannot be deleted; the ones before it are
/// gone.
pub fn remove(dir: &Path, e: &Entry) -> Result<Removed> {
    let files = files_of(dir, e);
    let mut freed = 0u64;
    for p in &files {
        freed += fs::metadata(p).map_or(0, |m| m.len());
        fs::remove_file(p).map_err(Error::io(p))?;
    }
    Ok(Removed {
        files: files.len(),
        bytes: freed,
    })
}

/// Something wrong with an installed catalogue.
#[derive(Debug)]
#[non_exhaustive]
pub enum Problem {
    /// A file is too short to be one of the catalogue's (every format has at least
    /// a 110-byte header or a FITS block).
    Truncated(PathBuf),
    /// A file could not be examined.
    Unreadable {
        /// The file.
        path: PathBuf,
        /// What the operating system said.
        error: std::io::Error,
    },
    /// The catalogue has the wrong number of files: its grid (an ASTAP database) or
    /// its set (downloaded indexes) fixes how many.
    FileCount {
        /// Files it should have.
        expected: usize,
        /// Files found.
        found: usize,
    },
}

impl fmt::Display for Problem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Truncated(p) => write!(f, "{} is truncated", p.display()),
            Self::Unreadable { path, error } => write!(f, "{}: {error}", path.display()),
            Self::FileCount { expected, found } => {
                write!(f, "expected {expected} files, found {found}")
            }
        }
    }
}

/// [`verify`]'s finding for one installed catalogue.
#[derive(Debug)]
pub struct CatalogueCheck {
    /// The catalogue.
    pub entry: &'static Entry,
    /// Its files present.
    pub files: usize,
    /// Bytes they occupy.
    pub bytes: u64,
    /// What is wrong; empty if it looks sound.
    pub problems: Vec<Problem>,
}

/// [`verify`]'s finding for one blind index file.
#[derive(Debug)]
pub struct IndexCheck {
    /// The file.
    pub path: PathBuf,
    /// What was found.
    pub health: IndexHealth,
}

/// The state of a blind index file.
#[derive(Debug)]
#[non_exhaustive]
pub enum IndexHealth {
    /// It opens and passes its integrity checks.
    Sound {
        /// Patterns in it.
        patterns: usize,
        /// File size, bytes.
        bytes: u64,
        /// Its header, if readable as an [`Existing`] index.
        existing: Option<Existing>,
        /// How it stands against its source database: never
        /// [`Freshness::Changed`], which is [`IndexHealth::Stale`] instead.
        freshness: Option<Freshness>,
    },
    /// Sound, but its source database has changed since it was built: it should
    /// be rebuilt.
    Stale {
        /// The source database.
        source: String,
    },
    /// It cannot be opened or fails its checks.
    Broken(ArcsecError),
}

impl IndexHealth {
    /// Whether this is a problem to report: stale or broken.
    #[must_use]
    pub const fn is_problem(&self) -> bool {
        !matches!(self, Self::Sound { .. })
    }
}

/// What [`verify`] found in a directory.
#[derive(Debug, Default)]
pub struct Verification {
    /// Every installed catalogue of the registry, in registry order.
    pub catalogues: Vec<CatalogueCheck>,
    /// Every blind index file (`*.arcsecix`), by name.
    pub indexes: Vec<IndexCheck>,
}

impl Verification {
    /// The number of problems found: per catalogue problem, and per stale or
    /// broken index.
    #[must_use]
    pub fn problems(&self) -> usize {
        self.catalogues
            .iter()
            .map(|c| c.problems.len())
            .sum::<usize>()
            + self
                .indexes
                .iter()
                .filter(|i| i.health.is_problem())
                .count()
    }

    /// Whether nothing at all is installed.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.catalogues.is_empty() && self.indexes.is_empty()
    }
}

/// Check that every catalogue installed in `dir` looks structurally sound: no
/// truncated files and the right number of them, and every blind index readable,
/// intact and built from the database now installed.
#[must_use]
pub fn verify(dir: &Path) -> Verification {
    let mut v = Verification::default();
    for e in REGISTRY {
        if !is_installed(dir, e) {
            continue;
        }
        let files = files_of(dir, e);
        let mut problems = Vec::new();
        for p in &files {
            match fs::metadata(p) {
                Ok(m) if m.len() < 120 => problems.push(Problem::Truncated(p.clone())),
                Err(error) => problems.push(Problem::Unreadable {
                    path: p.clone(),
                    error,
                }),
                _ => {}
            }
        }
        if let Some(want) = astap_file_count(dir, e).or_else(|| expected_file_count(e))
            && files.len() != want
        {
            problems.push(Problem::FileCount {
                expected: want,
                found: files.len(),
            });
        }
        v.catalogues.push(CatalogueCheck {
            entry: e,
            files: files.len(),
            bytes: installed_size(dir, e),
            problems,
        });
    }
    for path in index_files(dir) {
        let opened = arcsec_core::index::BlindIndex::open(&path).and_then(|ix| {
            ix.validate()?;
            Ok(ix)
        });
        let health = match opened {
            Ok(ix) => {
                let existing = Existing::open(&path);
                let fresh = existing.as_ref().map(|ex| freshness(ex, dir));
                if fresh == Some(Freshness::Changed) {
                    IndexHealth::Stale {
                        source: ix.source().to_string(),
                    }
                } else {
                    IndexHealth::Sound {
                        patterns: ix.n_patterns(),
                        bytes: ix.file_size() as u64,
                        existing,
                        freshness: fresh,
                    }
                }
            }
            Err(e) => IndexHealth::Broken(e),
        };
        v.indexes.push(IndexCheck { path, health });
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::registry::find;
    use crate::test_support::{TempDir, write_001_db};

    #[test]
    fn disk_needed_for_downloads_counts_archive_and_contents() {
        let d80 = find("d80").unwrap();
        assert_eq!(download_disk(d80), 2 * d80.bytes);
        let anet = find("anet-4100").unwrap();
        assert_eq!(download_disk(anet), anet.bytes);
    }

    #[test]
    fn recommendations_are_the_smallest_covering_download() {
        let id = |p, f| recommend(p, f).map(|e| e.id);
        assert_eq!(id(Purpose::Solving, 1.5), Some("d05"));
        assert_eq!(id(Purpose::Solving, 0.25), Some("d50"));
        assert_eq!(id(Purpose::Solving, 10.0), Some("g05"));
        assert_eq!(id(Purpose::Solving, 40.0), Some("w08"));
        assert_eq!(id(Purpose::Solving, 0.1), None);
        assert_eq!(id(Purpose::Photometry, 1.5), Some("v05"));
    }

    #[test]
    fn verify_and_remove_a_tiny_database() {
        let dir = TempDir::new("verify");
        let d = dir.path();
        assert!(verify(d).is_empty());
        write_001_db(d, "w08", 50, 1);
        let v = verify(d);
        assert_eq!(v.catalogues.len(), 1);
        assert_eq!(v.catalogues[0].entry.id, "w08");
        assert_eq!(v.problems(), 0, "{v:?}");

        // A truncated file, and an index that is not one.
        std::fs::write(d.join("w08_0101.001"), b"short").unwrap();
        std::fs::write(d.join("w08.arcsecix"), b"not an index").unwrap();
        let v = verify(d);
        assert!(matches!(
            v.catalogues[0].problems[..],
            [Problem::Truncated(_)]
        ));
        assert!(matches!(v.indexes[0].health, IndexHealth::Broken(_)));
        assert_eq!(v.problems(), 2);

        let r = remove(d, find("w08").unwrap()).unwrap();
        assert_eq!(r, Removed { files: 1, bytes: 5 });
        assert!(!is_installed(d, find("w08").unwrap()));
    }
}
