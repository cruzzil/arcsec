use core::fmt;

#[derive(Debug)]
pub enum ArcsecError {
    Singular,
    InsufficientStars { found: usize, required: usize },
    InsufficientQuads { found: usize, required: usize },
    BadSolution { ratio: f64 },
    CatalogNotFound(std::path::PathBuf),
    CatalogIo(std::io::Error),
}

impl fmt::Display for ArcsecError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ArcsecError::Singular => write!(f, "singular matrix in LSQ solver"),
            ArcsecError::InsufficientStars { found, required } => {
                write!(f, "insufficient stars: found {found}, required {required}")
            }
            ArcsecError::InsufficientQuads { found, required } => {
                write!(f, "insufficient quads: found {found}, required {required}")
            }
            ArcsecError::BadSolution { ratio } => {
                write!(
                    f,
                    "bad solution: xy scale ratio {ratio:.4} not in [0.9, 1.1]"
                )
            }
            ArcsecError::CatalogNotFound(p) => write!(f, "catalog not found: {}", p.display()),
            ArcsecError::CatalogIo(e) => write!(f, "catalog I/O error: {e}"),
        }
    }
}

impl core::error::Error for ArcsecError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            ArcsecError::CatalogIo(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for ArcsecError {
    fn from(e: std::io::Error) -> Self {
        ArcsecError::CatalogIo(e)
    }
}

pub type Result<T> = core::result::Result<T, ArcsecError>;
