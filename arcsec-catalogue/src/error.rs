//! What can go wrong managing catalogues.

use core::fmt;
use std::io;
use std::path::PathBuf;

use arcsec_core::ArcsecError;

use crate::format::human_bytes;

/// A catalogue operation that failed. The `Display` text is a complete sentence
/// fragment naming the file or download concerned, fit to show a user.
#[derive(Debug)]
#[non_exhaustive]
pub enum Error {
    /// A file or directory could not be read, written, renamed or removed.
    Io {
        /// The file or directory.
        path: PathBuf,
        /// What the operating system said.
        source: io::Error,
    },
    /// A download could not be made or was cut off: a network, TLS or read error.
    Network {
        /// What was being downloaded (a catalogue id, or a file of a set).
        label: String,
        /// The error.
        message: String,
    },
    /// The server answered a download with an HTTP status other than success.
    Http {
        /// What was being downloaded.
        label: String,
        /// The status code.
        status: u16,
    },
    /// A resumed download came back starting at the wrong byte. The partial file
    /// has been discarded; downloading again starts afresh.
    BadResume {
        /// What was being downloaded.
        label: String,
        /// The `Content-Range` the server sent, or `no Content-Range`.
        range: String,
    },
    /// The connection closed before the length the server announced. The partial
    /// file is kept, and downloading again resumes it.
    Truncated {
        /// What was being downloaded.
        label: String,
        /// Bytes received.
        written: u64,
        /// Bytes announced.
        total: u64,
    },
    /// An archive is not one this crate can unpack, or is corrupt.
    Archive {
        /// The archive.
        path: PathBuf,
        /// What is wrong with it.
        message: String,
    },
    /// A member of an archive could not be read.
    Malformed(String),
    /// Unpacking a downloaded archive failed. The archive is left in place, so the
    /// next install reuses it rather than downloading it again; delete it to
    /// download afresh.
    Unpack {
        /// The downloaded archive, kept.
        archive: PathBuf,
        /// Why unpacking failed.
        source: Box<Self>,
    },
    /// Some files of a set downloaded file by file failed. Those that arrived are
    /// kept, and installing again fetches only the rest.
    IncompleteSet {
        /// The catalogue.
        id: &'static str,
        /// Files present.
        present: usize,
        /// Files in the set.
        total: usize,
    },
    /// An archive unpacked, but none of the catalogue's files were in it.
    NothingInstalled {
        /// The catalogue.
        id: &'static str,
        /// The catalogue directory.
        dir: PathBuf,
    },
    /// The disk has too little space for the job, with the usual margin.
    NoSpace {
        /// The directory being written.
        dir: PathBuf,
        /// Bytes needed, margin included.
        needed: u64,
        /// Bytes free.
        free: u64,
    },
    /// No solving database is installed to build a blind index from.
    NoDatabase {
        /// Where it was looked for.
        dir: PathBuf,
    },
    /// The database named to build a blind index from is not installed.
    DatabaseMissing {
        /// The database (`d50`, ...).
        name: String,
        /// Where it was looked for.
        dir: PathBuf,
    },
    /// A blind index field range that is empty or not positive.
    BadFieldRange {
        /// Smallest field, degrees.
        min: f64,
        /// Largest field, degrees.
        max: f64,
    },
    /// No tier of the blind index serves the field range asked for.
    NoTier {
        /// Smallest field, degrees.
        min: f64,
        /// Largest field, degrees.
        max: f64,
        /// The source database.
        source: String,
    },
    /// The blind index file cannot be created where it was to go; found before the
    /// build starts, not after.
    NotWritable {
        /// The index file.
        path: PathBuf,
        /// What the operating system said.
        source: io::Error,
    },
    /// Building a blind index failed (a database file could not be read).
    Build(ArcsecError),
    /// A built blind index could not be written.
    Write {
        /// The index file.
        path: PathBuf,
        /// The error.
        source: ArcsecError,
    },
    /// The thread's [`arcsec_core::cancel`] token was cancelled. A download keeps
    /// its partial file to resume; nothing else is left half written.
    Cancelled,
}

/// A catalogue operation's result.
pub type Result<T> = core::result::Result<T, Error>;

impl Error {
    /// [`Error::Io`] for `path`.
    pub(crate) fn io(path: impl Into<PathBuf>) -> impl FnOnce(io::Error) -> Self {
        let path = path.into();
        move |source| Self::Io { path, source }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io { path, source } => write!(f, "{}: {source}", path.display()),
            Self::Network { label, message } => write!(f, "{label}: {message}"),
            Self::Http { label, status } => write!(f, "{label}: HTTP {status}"),
            Self::BadResume { label, range } => write!(
                f,
                "{label}: the server resumed at the wrong offset ({range}); the partial \
                 download was discarded, run the install again"
            ),
            Self::Truncated {
                label,
                written,
                total,
            } => write!(
                f,
                "{label}: download ended at {} of {}; run the install again to resume",
                human_bytes(*written),
                human_bytes(*total)
            ),
            Self::Archive { path, message } => write!(f, "{}: {message}", path.display()),
            Self::Malformed(message) => f.write_str(message),
            Self::Unpack { source, .. } => write!(f, "{source}"),
            Self::IncompleteSet { id, present, total } => write!(
                f,
                "{id}: {} of {total} files could not be downloaded; run the install again to fetch them",
                total - present
            ),
            Self::NothingInstalled { id, dir } => write!(
                f,
                "{id}: extraction finished but no catalogue files appeared in {}",
                dir.display()
            ),
            Self::NoSpace { dir, needed, free } => write!(
                f,
                "not enough free space in {}: this needs about {} and {} is free",
                dir.display(),
                human_bytes(*needed),
                human_bytes(*free)
            ),
            Self::NoDatabase { dir } => write!(f, "no star database in {}", dir.display()),
            Self::DatabaseMissing { name, dir } => {
                write!(f, "database {name} not found in {}", dir.display())
            }
            Self::BadFieldRange { min, max } => write!(f, "bad field range {min}°–{max}°"),
            Self::NoTier { min, max, source } => write!(
                f,
                "no tier fits fields {min}°–{max}° from {}",
                source.to_uppercase()
            ),
            Self::NotWritable { path, source } => {
                write!(f, "cannot write the index to {}: {source}", path.display())
            }
            Self::Build(e) => write!(f, "build failed: {e}"),
            Self::Write { path, source } => write!(f, "writing {}: {source}", path.display()),
            Self::Cancelled => f.write_str("cancelled"),
        }
    }
}

impl core::error::Error for Error {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Io { source, .. } | Self::NotWritable { source, .. } => Some(source),
            Self::Unpack { source, .. } => Some(source.as_ref()),
            Self::Build(e) | Self::Write { source: e, .. } => Some(e),
            _ => None,
        }
    }
}
