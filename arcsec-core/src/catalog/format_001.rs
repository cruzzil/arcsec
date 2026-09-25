//! Reader for ASTAP's .001 all-sky star database format (W08).
//!
//! Unlike .1476 and .290, which tile the sky into per-area files of packed 5-byte
//! records, .001 is a single whole-sky file of plain little-endian f32 triples. W08
//! is the only database shipped in this form: it targets 20°–80° fields, where any
//! tiling would be read in its entirety anyway, and at magnitude 8 the whole sky is
//! only ~41k stars.
//!
//! Layout, derived from the shipped `w08_star_database_mag08_astap` file and
//! confirmed against known stars — the first two records decode to Sirius
//! (mag -1.5, RA 101.28°, Dec -16.72°) and Canopus (mag -0.6, RA 95.99°, Dec -52.70°):
//!
//! ```text
//! u32  star_count                     (little-endian)
//! star_count × {
//!     f32  magnitude × 10
//!     f32  right ascension, radians
//!     f32  declination, radians
//! }
//! ```
//!
//! Records are ordered brightest first, so a truncated read is still the brightest
//! subset. `star_count * 12 + 4` accounts for the file exactly.

use core::f64::consts::PI;
use std::io;
use std::path::Path;

use memmap2::Mmap;

use super::format_1476::CatalogStar;
use crate::error::{ArcsecError, Result};

/// Bytes per star record: three f32 fields.
const RECORD_LEN: usize = 12;
/// Bytes of header before the first record: the u32 star count.
const HEADER_LEN: usize = 4;

/// Read stars from a `.001` all-sky file that fall within a square field.
///
/// Arguments mirror `format_1476::read_area_file` so the two are interchangeable
/// from the caller's point of view.
///
/// # Errors
///
/// [`ArcsecError::CatalogIo`] if the file cannot be opened or mapped, or its
/// declared star count does not match its length.
pub fn read_001_file(
    file_path: &Path,
    telescope_ra: f64,
    telescope_dec: f64,
    field_diameter: f64,
    cos_telescope_dec: f64,
    max_stars: usize,
) -> Result<Vec<CatalogStar>> {
    let file = std::fs::File::open(file_path).map_err(ArcsecError::CatalogIo)?;
    // Safety: read-only mapping; the catalog files are never written concurrently.
    let mmap = unsafe { Mmap::map(&file).map_err(ArcsecError::CatalogIo)? };

    if mmap.len() < HEADER_LEN {
        return Ok(vec![]);
    }
    let declared = u32::from_le_bytes([mmap[0], mmap[1], mmap[2], mmap[3]]) as usize;
    let available = (mmap.len() - HEADER_LEN) / RECORD_LEN;
    if declared != available {
        return Err(ArcsecError::CatalogIo(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "{}: header declares {declared} stars but the file holds {available}",
                file_path.display()
            ),
        )));
    }

    let half_diam = field_diameter * 0.5;
    let body = &mmap[HEADER_LEN..];
    let mut stars = Vec::new();

    for i in 0..available {
        let p = i * RECORD_LEN;
        let mag = f32::from_le_bytes([body[p], body[p + 1], body[p + 2], body[p + 3]]);
        let ra = f32::from_le_bytes([body[p + 4], body[p + 5], body[p + 6], body[p + 7]]) as f64;
        let dec = f32::from_le_bytes([body[p + 8], body[p + 9], body[p + 10], body[p + 11]]) as f64;

        // Square FOV filter, matching format_1476::read_area_file exactly.
        let mut delta_ra = (ra - telescope_ra).abs();
        if delta_ra > PI {
            delta_ra = 2.0 * PI - delta_ra;
        }
        if delta_ra * cos_telescope_dec < half_diam && (dec - telescope_dec).abs() < half_diam {
            stars.push(CatalogStar {
                ra,
                dec,
                mag: mag as f64 / 10.0,
            });
            if stars.len() >= max_stars {
                break;
            }
        }
    }

    Ok(stars)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a synthetic .001 file in memory.
    fn synth(stars: &[(f64, f64, f64)]) -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(&(stars.len() as u32).to_le_bytes());
        for &(mag, ra, dec) in stars {
            buf.extend_from_slice(&((mag * 10.0) as f32).to_le_bytes());
            buf.extend_from_slice(&(ra as f32).to_le_bytes());
            buf.extend_from_slice(&(dec as f32).to_le_bytes());
        }
        buf
    }

    fn write_temp(bytes: &[u8], name: &str) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(name);
        std::fs::write(&p, bytes).unwrap();
        p
    }

    #[test]
    fn reads_stars_inside_the_field_and_rejects_the_rest() {
        let bytes = synth(&[
            (-1.5, 1.0, 0.0),  // inside
            (2.0, 1.001, 0.0), // inside
            (3.0, 2.5, 0.0),   // far in RA
            (4.0, 1.0, 1.0),   // far in Dec
        ]);
        let p = write_temp(&bytes, "arcsec_test_001_basic.001");
        let got = read_001_file(&p, 1.0, 0.0, 0.05, 1.0, 100).unwrap();
        assert_eq!(got.len(), 2);
        assert!(
            (got[0].mag - -1.5).abs() < 1e-6,
            "mag decoded as {}",
            got[0].mag
        );
        assert!((got[0].ra - 1.0).abs() < 1e-6);
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn respects_max_stars() {
        let bytes = synth(&[(1.0, 1.0, 0.0), (2.0, 1.0, 0.0), (3.0, 1.0, 0.0)]);
        let p = write_temp(&bytes, "arcsec_test_001_cap.001");
        let got = read_001_file(&p, 1.0, 0.0, 0.05, 1.0, 2).unwrap();
        assert_eq!(got.len(), 2);
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn ra_wrap_is_handled() {
        // A field at RA ~0 must pick up a star just below 2π.
        let bytes = synth(&[(1.0, 2.0 * PI - 0.001, 0.0)]);
        let p = write_temp(&bytes, "arcsec_test_001_wrap.001");
        let got = read_001_file(&p, 0.0, 0.0, 0.05, 1.0, 10).unwrap();
        assert_eq!(got.len(), 1, "star across the RA=0 seam was rejected");
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn truncated_file_is_an_error_not_garbage() {
        let mut bytes = synth(&[(1.0, 1.0, 0.0), (2.0, 1.0, 0.0)]);
        bytes.truncate(bytes.len() - 5);
        let p = write_temp(&bytes, "arcsec_test_001_trunc.001");
        assert!(read_001_file(&p, 1.0, 0.0, 0.05, 1.0, 10).is_err());
        std::fs::remove_file(&p).ok();
    }
}
