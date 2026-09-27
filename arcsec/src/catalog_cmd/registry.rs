//! Every catalogue `arcsec catalog` knows how to install, and how to recognise its
//! files on disk.

use std::fs;
use std::path::{Path, PathBuf};

/// What a catalogue is for. Users pick by task, not by file format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Purpose {
    /// Star database for the catalogue (spiral) solver.
    Solving,
    /// Carries photometric magnitudes and colour, for photometric colour calibration.
    Photometry,
    /// Astrometry.net index files, for blind solving with no position hint.
    BlindIndex,
}

impl Purpose {
    pub const fn label(self) -> &'static str {
        match self {
            Self::Solving => "solving",
            Self::Photometry => "photometry",
            Self::BlindIndex => "blind",
        }
    }
}

/// How a catalogue's bytes arrive.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Archive {
    /// A `.zip`, files at the archive root.
    Zip,
    /// A Debian package: an `ar` archive whose `data.tar.xz` holds the files.
    Deb,
    /// Plain files, downloaded individually (see [`loose_files`]).
    Loose,
}

pub struct Entry {
    /// Short name the user types.
    pub id: &'static str,
    pub purpose: Purpose,
    /// Download URL. For `Loose` sets this is a descriptive template only; the
    /// individual URLs come from [`loose_files`].
    pub url: &'static str,
    pub archive: Archive,
    /// Approximate download size in bytes, for the confirmation prompt.
    pub bytes: u64,
    /// Field-of-view range this catalogue is built for, in degrees.
    pub fov: Option<(f64, f64)>,
    /// How to recognise this catalogue's files on disk.
    pub files: Files,
    pub desc: &'static str,
}

/// File naming of an installed catalogue.
///
/// ASTAP databases are `<prefix>_RRCC.<ext>` where the extension says which sky grid
/// they use — and which extension a given database ships is not something to guess.
/// V05, for instance, is `.290` despite covering the same field range as the `.1476`
/// D-series, and being wrong about it made a successful install report as a failure.
/// So probe for the first cell under any known extension instead.
///
/// Astrometry.net indexes are recognised by their exact file names - the names
/// [`loose_files`] downloads - not by a prefix. The 4100 and 5200 sets share a
/// directory with each other and possibly with index files arcsec did not install
/// (the full 5200 series, 5203-5206, ...), and `remove` deletes whatever a set
/// claims, so a prefix match would delete files that are not arcsec's.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Files {
    AstapDb { prefix: &'static str },
    AnetIndex { set: AnetSet },
}

/// The Astrometry.net index sets arcsec can install.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnetSet {
    /// `index-4107.fits` ... `index-4119.fits`.
    Tycho4100,
    /// `index-5200-00.fits` ... `index-5202-47.fits` (the LITE files).
    GaiaLite5200,
}

impl AnetSet {
    /// Is `name` exactly one of this set's file names?
    fn owns(self, name: &str) -> bool {
        let Some(stem) = name
            .strip_prefix("index-")
            .and_then(|n| n.strip_suffix(".fits"))
        else {
            return false;
        };
        let two_digits = |s: &str| -> Option<u32> {
            (s.len() == 2 && s.bytes().all(|b| b.is_ascii_digit()))
                .then(|| s.parse().ok())
                .flatten()
        };
        match self {
            Self::Tycho4100 => stem
                .strip_prefix("41")
                .and_then(two_digits)
                .is_some_and(|scale| (7..=19).contains(&scale)),
            Self::GaiaLite5200 => stem.split_once('-').is_some_and(|(series, hp)| {
                matches!(series, "5200" | "5201" | "5202")
                    && two_digits(hp).is_some_and(|hp| hp < 48)
            }),
        }
    }
}

/// Extensions an ASTAP star database can use, and the file count each implies.
pub const ASTAP_EXTS: &[(&str, usize)] = &[(".1476", 1476), (".290", 290), (".001", 1)];

impl Files {
    /// Does `name` belong to this catalogue?
    pub fn owns(&self, name: &str) -> bool {
        match *self {
            Self::AstapDb { prefix } => {
                name.strip_prefix(prefix)
                    .is_some_and(|rest| rest.starts_with('_'))
                    && ASTAP_EXTS.iter().any(|(e, _)| name.ends_with(e))
            }
            Self::AnetIndex { set } => set.owns(name),
        }
    }

    /// A file whose presence means the catalogue is installed.
    fn probe(&self, dir: &Path) -> Option<PathBuf> {
        match *self {
            Self::AstapDb { prefix } => ASTAP_EXTS
                .iter()
                .map(|(e, _)| dir.join(format!("{prefix}_0101{e}")))
                .find(|p| p.exists()),
            Self::AnetIndex { .. } => installed_files(dir, self).into_iter().next(),
        }
    }
}

/// Every catalogue arcsec knows how to install.
///
/// Sizes and URLs are from the ASTAP and astrometry.net download pages, verified
/// 2026-09-04. Note that D80, V50 and V05 are published only as `.deb`/`.exe`/`.pkg`
/// — there is no `.zip` — which is why `Archive::Deb` support exists at all.
pub const REGISTRY: &[Entry] = &[
    // ── Star databases for solving ──────────────────────────────────────────────
    Entry {
        id: "d05",
        purpose: Purpose::Solving,
        url: "https://sourceforge.net/projects/astap-program/files/star_databases/d05_star_database.zip/download",
        archive: Archive::Zip,
        bytes: 102_200_000,
        fov: Some((0.6, 6.0)),
        files: Files::AstapDb { prefix: "d05" },
        desc: "Gaia DR3 to 500 stars/deg². Smallest useful solving database.",
    },
    Entry {
        id: "d20",
        purpose: Purpose::Solving,
        url: "https://sourceforge.net/projects/astap-program/files/star_databases/d20_star_database.zip/download",
        archive: Archive::Zip,
        bytes: 399_600_000,
        fov: Some((0.3, 6.0)),
        files: Files::AstapDb { prefix: "d20" },
        desc: "Gaia DR3 to 2000 stars/deg².",
    },
    Entry {
        id: "d50",
        purpose: Purpose::Solving,
        url: "https://sourceforge.net/projects/astap-program/files/star_databases/d50_star_database.zip/download",
        archive: Archive::Zip,
        bytes: 901_300_000,
        fov: Some((0.2, 6.0)),
        files: Files::AstapDb { prefix: "d50" },
        desc: "Gaia DR3 to 5000 stars/deg². The usual choice for solving.",
    },
    Entry {
        id: "d80",
        purpose: Purpose::Solving,
        url: "https://sourceforge.net/projects/astap-program/files/star_databases/d80_star_database.deb/download",
        archive: Archive::Deb,
        bytes: 1_213_400_000,
        fov: Some((0.15, 6.0)),
        files: Files::AstapDb { prefix: "d80" },
        desc: "Gaia DR3 to 8000 stars/deg². Densest; needed below ~0.2° fields.",
    },
    Entry {
        id: "g05",
        purpose: Purpose::Solving,
        url: "https://sourceforge.net/projects/astap-program/files/star_databases/g05_star_database.zip/download",
        archive: Archive::Zip,
        bytes: 101_600_000,
        fov: Some((3.0, 20.0)),
        files: Files::AstapDb { prefix: "g05" },
        desc: "Wide fields, 3°–20°. The D-series stops at 6°.",
    },
    Entry {
        id: "w08",
        purpose: Purpose::Solving,
        url: "https://sourceforge.net/projects/astap-program/files/star_databases/w08_star_database_mag08_astap.zip/download",
        archive: Archive::Zip,
        bytes: 330_000,
        fov: Some((20.0, 80.0)),
        files: Files::AstapDb { prefix: "w08" },
        desc: "Very wide fields, 20°–80°, to magnitude 8. Tiny.",
    },
    // ── Photometric catalogues (colour calibration) ─────────────────────────────
    Entry {
        id: "v05",
        purpose: Purpose::Photometry,
        url: "https://sourceforge.net/projects/astap-program/files/star_databases/v05_star_database.deb/download",
        archive: Archive::Deb,
        bytes: 116_900_000,
        fov: Some((0.6, 6.0)),
        files: Files::AstapDb { prefix: "v05" },
        desc: "Johnson-V magnitudes plus Gaia BP-RP colour, 500 stars/deg².",
    },
    Entry {
        id: "v50",
        purpose: Purpose::Photometry,
        url: "https://sourceforge.net/projects/astap-program/files/star_databases/v50_star_database.deb/download",
        archive: Archive::Deb,
        bytes: 1_011_000_000,
        fov: Some((0.2, 6.0)),
        files: Files::AstapDb { prefix: "v50" },
        desc: "Johnson-V plus BP-RP colour, 5000 stars/deg². Deeper photometry.",
    },
    // ── Astrometry.net indexes for blind solving ────────────────────────────────
    Entry {
        id: "anet-4100",
        purpose: Purpose::BlindIndex,
        url: "https://data.astrometry.net/4100/index-41{:02}.fits",
        archive: Archive::Loose,
        bytes: 355_500_000,
        fov: Some((0.7, 180.0)),
        files: Files::AnetIndex {
            set: AnetSet::Tycho4100,
        },
        desc: "Tycho-2 blind indexes, scales 07–19 (fields ~0.7° and wider).",
    },
    Entry {
        id: "anet-5200",
        purpose: Purpose::BlindIndex,
        url: "https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/index-52{:02}-{:02}.fits",
        archive: Archive::Loose,
        bytes: 8_800_000_000,
        fov: Some((0.1, 2.0)),
        files: Files::AnetIndex {
            set: AnetSet::GaiaLite5200,
        },
        desc: "Gaia LITE blind indexes 5200/5201/5202, 48 HEALPix each. Large.",
    },
];

/// Look up a catalogue by the name the user typed, ignoring case.
pub fn find(id: &str) -> Option<&'static Entry> {
    REGISTRY.iter().find(|e| e.id.eq_ignore_ascii_case(id))
}

/// `(url, file name)` of every file in a `Loose` entry; empty for any other.
pub fn loose_files(e: &Entry) -> Vec<(String, String)> {
    match e.id {
        "anet-4100" => (7..=19)
            .map(|scale| {
                (
                    format!("https://data.astrometry.net/4100/index-41{scale:02}.fits"),
                    format!("index-41{scale:02}.fits"),
                )
            })
            .collect(),
        "anet-5200" => [5200u32, 5201, 5202]
            .into_iter()
            .flat_map(|series| (0..48).map(move |hp| (series, hp)))
            .map(|(series, hp)| {
                (
                    format!(
                        "https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/index-{series}-{hp:02}.fits"
                    ),
                    format!("index-{series}-{hp:02}.fits"),
                )
            })
            .collect(),
        _ => Vec::new(),
    }
}

/// Is `e` present in `dir`? (At least partly, for a `Loose` set.)
pub fn is_installed(dir: &Path, e: &Entry) -> bool {
    e.files.probe(dir).is_some()
}

/// Is every file of `e` present in `dir`?
///
/// Differs from [`is_installed`] only for `Loose` sets, where an interrupted
/// install leaves some of the files; `install` must then fetch the rest rather
/// than report the set as already installed.
pub fn is_complete(dir: &Path, e: &Entry) -> bool {
    match e.archive {
        Archive::Loose => loose_files(e)
            .iter()
            .all(|(_, name)| dir.join(name).is_file()),
        Archive::Zip | Archive::Deb => is_installed(dir, e),
    }
}

/// Files belonging to `files` that are present in `dir`, in name order.
fn installed_files(dir: &Path, files: &Files) -> Vec<PathBuf> {
    let mut out: Vec<PathBuf> = fs::read_dir(dir)
        .into_iter()
        .flatten()
        .flatten()
        .map(|d| d.path())
        .filter(|p| {
            p.file_name()
                .is_some_and(|n| files.owns(&n.to_string_lossy()))
        })
        .collect();
    out.sort();
    out
}

/// Files belonging to `e` that are present in `dir`.
pub fn files_of(dir: &Path, e: &Entry) -> Vec<PathBuf> {
    installed_files(dir, &e.files)
}

/// Bytes an installed catalogue occupies on disk.
pub fn installed_size(dir: &Path, e: &Entry) -> u64 {
    files_of(dir, e)
        .iter()
        .filter_map(|p| fs::metadata(p).ok())
        .map(|m| m.len())
        .sum()
}

/// How many files a complete install of `e` has, when that is known before looking
/// at the disk. For an ASTAP database it depends on the grid, so see `verify`.
pub fn expected_file_count(e: &Entry) -> Option<usize> {
    match e.archive {
        Archive::Loose => Some(loose_files(e).len()),
        Archive::Zip | Archive::Deb => None,
    }
}

/// The number of files an installed ASTAP database should have, from the extension
/// of the file that [`is_installed`] found.
pub fn astap_file_count(dir: &Path, e: &Entry) -> Option<usize> {
    let Files::AstapDb { .. } = e.files else {
        return None;
    };
    let probe = e.files.probe(dir)?;
    let name = probe.file_name()?.to_string_lossy().into_owned();
    ASTAP_EXTS
        .iter()
        .find(|(ext, _)| name.ends_with(ext))
        .map(|(_, n)| *n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn registry_is_well_formed() {
        for e in REGISTRY {
            assert!(!e.id.is_empty());
            assert!(!e.desc.is_empty(), "{} has no description", e.id);
            assert!(e.bytes > 0, "{} has no size", e.id);
            assert!(e.url.starts_with("https://"), "{} is not https", e.id);
            if let Some((lo, hi)) = e.fov {
                assert!(lo < hi, "{} has an inverted FOV range", e.id);
            }
            if let Files::AstapDb { prefix } = e.files {
                assert_eq!(prefix, e.id, "{}: prefix should match the id", e.id);
            }
            // Every file a Loose set downloads must be recognised as its own, or
            // install would report success on files that list/remove cannot see.
            for (url, name) in loose_files(e) {
                assert!(url.ends_with(&name), "{url} does not fetch {name}");
                assert!(e.files.owns(&name), "{}: {name} is not recognised", e.id);
            }
        }
        let mut ids: Vec<&str> = REGISTRY.iter().map(|e| e.id).collect();
        ids.sort_unstable();
        let n = ids.len();
        ids.dedup();
        assert_eq!(ids.len(), n, "duplicate catalogue id");
    }

    #[test]
    fn every_solving_field_size_has_a_catalogue() {
        for fov in [0.2, 0.5, 1.0, 3.0, 5.0, 10.0, 30.0, 60.0] {
            let any = REGISTRY.iter().any(|e| {
                e.purpose == Purpose::Solving
                    && matches!(e.fov, Some((lo, hi)) if fov >= lo && fov <= hi)
            });
            assert!(any, "no solving catalogue covers a {fov}° field");
        }
    }

    #[test]
    fn loose_url_sets_are_complete() {
        let a = loose_files(find("anet-4100").unwrap());
        assert_eq!(a.len(), 13, "4100 series should be scales 07..19");
        assert!(a[0].0.ends_with("index-4107.fits"));
        let b = loose_files(find("anet-5200").unwrap());
        assert_eq!(b.len(), 3 * 48, "5200 series is 3 scales x 48 healpix");
        assert!(b.iter().all(|(u, _)| u.starts_with("https://")));
    }

    #[test]
    fn index_sets_do_not_claim_each_others_files() {
        let a = find("anet-4100").unwrap().files;
        let b = find("anet-5200").unwrap().files;
        assert!(a.owns("index-4107.fits") && !b.owns("index-4107.fits"));
        assert!(b.owns("index-5200-07.fits") && !a.owns("index-5200-07.fits"));
        assert!(!a.owns("index-4107.fits.part"));
    }

    #[test]
    fn index_sets_claim_only_the_files_they_install() {
        // Other series share the directory and must survive `remove`.
        for e in ["anet-4100", "anet-5200"].map(|id| find(id).unwrap()) {
            for (_, name) in loose_files(e) {
                assert!(e.files.owns(&name), "{}: {name} not recognised", e.id);
            }
        }
        let a = find("anet-4100").unwrap().files;
        let b = find("anet-5200").unwrap().files;
        for foreign in [
            "index-4106.fits",
            "index-4120.fits",
            "index-4107-00.fits",
            "index-5203-00.fits",
            "index-5206-47.fits",
            "index-5200-48.fits",
            "index-5200-7.fits",
            "index-5200.fits",
        ] {
            assert!(!a.owns(foreign) && !b.owns(foreign), "{foreign} claimed");
        }
    }

    #[test]
    fn astap_prefixes_need_the_separator() {
        let d5 = Files::AstapDb { prefix: "d5" };
        assert!(!d5.owns("d50_0101.1476"), "d5 must not claim d50's files");
        let d50 = Files::AstapDb { prefix: "d50" };
        assert!(d50.owns("d50_0101.1476"));
        assert!(!d50.owns("d50_0101.1476.part"));
    }
}
