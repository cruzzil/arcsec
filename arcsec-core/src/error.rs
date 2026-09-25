//! The crate's error type.

use core::fmt;

/// Everything that can stop a solve.
///
/// The variants map onto the ASTAP-compatible exit codes the `arcsec` CLI returns
/// (2 for too few stars, 1 for no solution, 32/33 for database problems), which is
/// why "no match" is reported through [`ArcsecError::InsufficientQuads`] rather than
/// a dedicated variant. New variants may be added, so match with a wildcard arm.
#[derive(Debug)]
#[non_exhaustive]
pub enum ArcsecError {
    /// The least-squares system was singular (degenerate or too few points).
    Singular,
    /// Fewer stars were detected than the solver needs.
    InsufficientStars {
        /// Stars actually detected.
        found: usize,
        /// Minimum the solver requires.
        required: usize,
    },
    /// Pattern matching did not produce a verified solution.
    ///
    /// Returned both when too few patterns could be built and when the search
    /// finished without any position verifying.
    InsufficientQuads {
        /// Patterns (or, for the blind solver, verification score) achieved.
        found: usize,
        /// Minimum required.
        required: usize,
    },
    /// The plate fit's X and Y scales disagree by more than 10%.
    BadSolution {
        /// Ratio of the squared X scale to the squared Y scale.
        ratio: f64,
    },
    /// The star database directory does not contain the named database.
    CatalogNotFound(std::path::PathBuf),
    /// A catalogue or index file could not be read or parsed.
    CatalogIo(std::io::Error),
    /// A caller-supplied parameter is out of range (e.g. a non-positive field of view).
    InvalidParameter(String),
}

impl fmt::Display for ArcsecError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Singular => write!(f, "singular matrix in LSQ solver"),
            Self::InsufficientStars { found, required } => {
                write!(f, "insufficient stars: found {found}, required {required}")
            }
            Self::InsufficientQuads { found, required } => {
                write!(f, "insufficient quads: found {found}, required {required}")
            }
            Self::BadSolution { ratio } => {
                write!(
                    f,
                    "bad solution: xy scale ratio {ratio:.4} not in [0.9, 1.1]"
                )
            }
            Self::CatalogNotFound(p) => write!(f, "catalog not found: {}", p.display()),
            Self::CatalogIo(e) => write!(f, "catalog I/O error: {e}"),
            Self::InvalidParameter(msg) => write!(f, "invalid parameter: {msg}"),
        }
    }
}

impl core::error::Error for ArcsecError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::CatalogIo(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for ArcsecError {
    fn from(e: std::io::Error) -> Self {
        Self::CatalogIo(e)
    }
}

/// `Result` specialised to [`ArcsecError`].
pub type Result<T> = core::result::Result<T, ArcsecError>;
