//! Installing a catalogue from the registry: download, unpack, move into place.

use std::fs;
use std::path::{Path, PathBuf};

use crate::error::{Error, Result};
use crate::fetch::{self, DownloadEvent};
use crate::registry::{Archive, Entry, installed_size, is_installed, loose_files};

/// Options for [`install`].
#[derive(Debug, Clone, Copy, Default)]
pub struct InstallOptions {
    /// Keep a downloaded `.zip` or `.deb` in the catalogue directory after
    /// unpacking it (as `.<id>-download.<ext>`), instead of deleting it.
    pub keep_archive: bool,
}

/// What [`install`] reports while it runs, in order.
#[derive(Debug)]
#[non_exhaustive]
pub enum InstallEvent<'a> {
    /// Progress of a download. `label` names it: the catalogue's id for an
    /// archive, `anet-4100 [3/13] index-4109.fits` for a file of a set.
    Download {
        /// What is downloading.
        label: &'a str,
        /// How far it has got.
        event: DownloadEvent,
    },
    /// A file of a set downloaded file by file failed. The others continue, and the
    /// install then fails with [`Error::IncompleteSet`].
    FileFailed(&'a Error),
    /// Every file of a set has been tried: this many of the set are present.
    FilesPresent {
        /// Files present.
        present: usize,
        /// Files in the set.
        total: usize,
    },
    /// The archive has downloaded and is being unpacked.
    Extracting,
    /// A `.deb`'s compressed payload (named) is being decompressed.
    Decompressing(&'a str),
    /// The archive was kept, at this path ([`InstallOptions::keep_archive`]).
    ArchiveKept(&'a Path),
    /// The archive has been unpacked: this many files are in place.
    Extracted(usize),
}

/// What [`install`] installed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Installed {
    /// Files the catalogue has in the directory now.
    pub files: usize,
    /// Bytes they occupy.
    pub bytes: u64,
}

/// Download catalogue `e` into `dir` (created if need be) and unpack it there,
/// reporting progress to `on_event`.
///
/// A `.zip` or `.deb` is downloaded beside the catalogue as `.<id>-download.<ext>`
/// (kept on failure, so a retry reuses it), and its files are extracted into a
/// staging directory and moved into place only once the whole archive has
/// unpacked: an install interrupted part way leaves nothing that looks installed.
/// A set of loose files (the Astrometry.net indexes) is fetched file by file,
/// skipping those already present.
///
/// Building the blind index afterwards is separate: see [`crate::index`].
///
/// Stops with [`Error::Cancelled`] once the thread's [`arcsec_core::cancel`] token
/// is cancelled, keeping partial downloads to resume.
///
/// ```no_run
/// use arcsec_catalogue::{InstallOptions, default_dir, install, registry};
/// use arcsec_core::cancel::{CancelToken, with_token};
///
/// let d50 = registry::find("d50").expect("a known catalogue");
/// let token = CancelToken::new(); // cancel a clone from another thread to stop
/// let installed = with_token(&token, || {
///     install(&default_dir(), d50, &InstallOptions::default(), &mut |event| {
///         println!("{event:?}");
///     })
/// })?;
/// println!("{} files, {} bytes", installed.files, installed.bytes);
/// # Ok::<(), arcsec_catalogue::Error>(())
/// ```
///
/// # Errors
///
/// A download or unpacking failure, [`Error::IncompleteSet`] if some files of a set
/// could not be fetched, [`Error::NothingInstalled`] if the archive held none of
/// the catalogue's files, or [`Error::Cancelled`].
pub fn install(
    dir: &Path,
    e: &'static Entry,
    opts: &InstallOptions,
    on_event: &mut dyn FnMut(InstallEvent<'_>),
) -> Result<Installed> {
    fs::create_dir_all(dir).map_err(Error::io(dir))?;
    match e.archive {
        Archive::Loose => install_loose(dir, e, on_event)?,
        Archive::Zip | Archive::Deb => install_archive(dir, e, opts.keep_archive, on_event)?,
    }
    if !is_installed(dir, e) {
        return Err(Error::NothingInstalled {
            id: e.id,
            dir: dir.to_path_buf(),
        });
    }
    Ok(Installed {
        files: crate::registry::files_of(dir, e).len(),
        bytes: installed_size(dir, e),
    })
}

/// Download each file of a `Loose` set that is not already present.
///
/// A file that fails is reported and the rest continue; the set as a whole then
/// fails, so a script sees it, and a re-run fetches only what is missing.
fn install_loose(
    dir: &Path,
    e: &'static Entry,
    on_event: &mut dyn FnMut(InstallEvent<'_>),
) -> Result<()> {
    let files = loose_files(e);
    let mut done = 0;
    for (i, (url, name)) in files.iter().enumerate() {
        let dest = dir.join(name);
        if dest.is_file() {
            done += 1;
            continue;
        }
        let label = format!("{} [{}/{}] {name}", e.id, i + 1, files.len());
        let result = fetch::download(url, &dest, &label, &mut |event| {
            on_event(InstallEvent::Download {
                label: &label,
                event,
            });
        });
        match result {
            Ok(()) => done += 1,
            Err(Error::Cancelled) => return Err(Error::Cancelled),
            Err(err) => on_event(InstallEvent::FileFailed(&err)),
        }
    }
    on_event(InstallEvent::FilesPresent {
        present: done,
        total: files.len(),
    });
    if done < files.len() {
        return Err(Error::IncompleteSet {
            id: e.id,
            present: done,
            total: files.len(),
        });
    }
    Ok(())
}

/// Where `install` downloads the archive of `e` in `dir`.
fn archive_path(dir: &Path, e: &Entry) -> PathBuf {
    let ext = if e.archive == Archive::Zip {
        "zip"
    } else {
        "deb"
    };
    dir.join(format!(".{}-download.{ext}", e.id))
}

/// Download a `.zip` or `.deb` and unpack its catalogue files into `dir`.
///
/// Files are extracted into a staging directory first and moved into place only
/// once the whole archive has unpacked. An install interrupted part way therefore
/// leaves nothing that looks installed — `is_installed` probes a single file, so a
/// half-extracted database would otherwise be reported as present, and skipped by
/// the next install.
fn install_archive(
    dir: &Path,
    e: &'static Entry,
    keep: bool,
    on_event: &mut dyn FnMut(InstallEvent<'_>),
) -> Result<()> {
    let tmp = archive_path(dir, e);
    if !tmp.is_file() {
        fetch::download(e.url, &tmp, e.id, &mut |event| {
            on_event(InstallEvent::Download { label: e.id, event });
        })?;
    }

    let staging = dir.join(format!(".{}-staging", e.id));
    let _ = fs::remove_dir_all(&staging);
    fs::create_dir_all(&staging).map_err(Error::io(&staging))?;

    on_event(InstallEvent::Extracting);
    let wanted = |name: &str| e.files.owns(name);
    let extracted = if e.archive == Archive::Zip {
        fetch::extract_zip(&tmp, &staging, &wanted)
    } else {
        fetch::extract_deb(&tmp, &staging, &wanted, &mut |name| {
            on_event(InstallEvent::Decompressing(name));
        })
    };
    let moved = extracted.and_then(|n| {
        for entry in fs::read_dir(&staging).map_err(Error::io(&staging))? {
            let from = entry
                .map_err(|err| Error::Malformed(err.to_string()))?
                .path();
            let Some(name) = from.file_name() else {
                continue;
            };
            let to = dir.join(name);
            fs::rename(&from, &to).map_err(Error::io(to))?;
        }
        Ok(n)
    });
    let _ = fs::remove_dir_all(&staging);
    let n = moved.map_err(|err| match err {
        Error::Cancelled => Error::Cancelled,
        err => Error::Unpack {
            archive: tmp.clone(),
            source: Box::new(err),
        },
    })?;

    if keep {
        on_event(InstallEvent::ArchiveKept(&tmp));
    } else {
        let _ = fs::remove_file(&tmp);
    }
    on_event(InstallEvent::Extracted(n));
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::registry::find;
    use crate::test_support::TempDir;
    use std::io::Write as _;

    /// An archive already downloaded (as an interrupted install leaves it) is
    /// unpacked without a download, and the install reports each step.
    #[test]
    fn a_downloaded_archive_is_unpacked_into_place() {
        let dir = TempDir::new("inst_zip");
        let d = dir.path();
        let w08 = find("w08").unwrap();
        {
            let f = fs::File::create(archive_path(d, w08)).unwrap();
            let mut w = zip::ZipWriter::new(f);
            let o = zip::write::SimpleFileOptions::default();
            for name in ["w08/w08_0101.001", "readme.txt"] {
                w.start_file(name, o).unwrap();
                w.write_all(&[0u8; 200]).unwrap();
            }
            w.finish().unwrap();
        }
        let mut seen = Vec::new();
        let got = install(d, w08, &InstallOptions { keep_archive: true }, &mut |ev| {
            seen.push(format!("{ev:?}"));
        })
        .unwrap();
        assert_eq!(
            got,
            Installed {
                files: 1,
                bytes: 200
            }
        );
        assert!(d.join("w08_0101.001").is_file() && !d.join("readme.txt").exists());
        assert_eq!(seen.len(), 3, "{seen:?}");
        assert!(
            seen[0] == "Extracting" && seen[2] == "Extracted(1)",
            "{seen:?}"
        );
        assert!(archive_path(d, w08).is_file(), "kept");
        assert!(!d.join(".w08-staging").exists());
    }

    /// An archive with none of the catalogue's files is a failure, not a success.
    #[test]
    fn an_archive_without_the_catalogue_fails() {
        let dir = TempDir::new("inst_empty");
        let d = dir.path();
        let w08 = find("w08").unwrap();
        {
            let mut w = zip::ZipWriter::new(fs::File::create(archive_path(d, w08)).unwrap());
            w.start_file("readme.txt", zip::write::SimpleFileOptions::default())
                .unwrap();
            w.finish().unwrap();
        }
        let r = install(d, w08, &InstallOptions::default(), &mut |_| {});
        assert!(matches!(r, Err(Error::NothingInstalled { .. })), "{r:?}");
        assert!(!archive_path(d, w08).exists(), "not kept");
    }

    /// A corrupt archive is kept for the next attempt and named in the error.
    #[test]
    fn a_corrupt_archive_is_kept() {
        let dir = TempDir::new("inst_bad");
        let d = dir.path();
        let w08 = find("w08").unwrap();
        fs::write(archive_path(d, w08), b"not a zip").unwrap();
        let r = install(d, w08, &InstallOptions::default(), &mut |_| {});
        let Err(Error::Unpack { archive, .. }) = r else {
            panic!("{r:?}");
        };
        assert!(archive.is_file());
    }
}
