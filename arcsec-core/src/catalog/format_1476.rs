// Binary reader for ASTAP's packed Gaia catalog formats, .1476 and .290.
//
// Both carry an identical 110-byte ASCII header (byte 109 = record size) followed by
// the same packed records, so this reader handles both; only the sky tiling differs
// (`areas.rs` for 1476, `areas_290.rs` for 290).
// Record format (5 bytes):
//   ra7, ra8, ra9  — 24-bit LE unsigned integer for RA
//   dec7, dec8     — low 2 bytes of a 24-bit two's-complement DEC integer
//
// Header records have ra_raw == 0xFF_FF_FF and carry the high DEC byte (dec9_storage)
// and magnitude for the following group of star records.
//
// RA  = ra_raw / (2^24 - 1) * 2π
// DEC = signed_24bit_int(dec9_storage:dec8:dec7) / (2^23 - 1) * (π/2)

use core::f64::consts::PI;
use std::io;
use std::path::Path;

use memmap2::Mmap;

use super::areas::filename_1476;
use crate::error::{ArcsecError, Result};

const RA_SCALE: f64 = 2.0 * PI / 16_777_215.0; // 2π / (2^24 - 1)
const DEC_SCALE: f64 = PI * 0.5 / 8_388_607.0; // π/2 / (2^23 - 1)
const HEADER_SENTINEL: u32 = 0xFF_FF_FF;

/// A catalog star — just RA/DEC in radians and magnitude.
#[derive(Debug, Clone)]
pub struct CatalogStar {
    pub ra: f64,
    pub dec: f64,
    pub mag: f64,
}

/// Read all stars from a single .1476 area file that fall within a square FOV.
///
/// - `file_path`: full path to the area file (e.g. `/db/d20_0101.1476`)
/// - `telescope_ra/dec`: centre of the field (radians)
/// - `field_diameter`: diameter of the square area to collect (radians)
/// - `cos_telescope_dec`: pre-computed cos(telescope_dec) for fast RA delta check
/// - `max_stars`: maximum stars to return (stops reading after this many)
pub fn read_area_file(
    file_path: &Path,
    telescope_ra: f64,
    telescope_dec: f64,
    field_diameter: f64,
    cos_telescope_dec: f64,
    max_stars: usize,
) -> Result<Vec<CatalogStar>> {
    let file = std::fs::File::open(file_path).map_err(ArcsecError::CatalogIo)?;

    // Memory-map the file: zero-copy, avoids the read_to_end heap allocation
    // and the kernel→user memcpy that was 18% of flamegraph samples.
    // Safety: we only read; the catalog files are never written concurrently.
    let mmap = unsafe { Mmap::map(&file).map_err(ArcsecError::CatalogIo)? };

    if mmap.len() < 110 {
        return Ok(vec![]);
    }

    // 110-byte header; last byte encodes record size
    let record_size = if mmap[109] == b' ' {
        11
    } else {
        mmap[109] as usize
    };

    // We only handle record_size 5 (RA+DEC) and 6 (RA+DEC+Gaia colour)
    if record_size != 5 && record_size != 6 {
        return Err(ArcsecError::CatalogIo(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("unsupported record_size {record_size}"),
        )));
    }

    // Parse directly from the memory-mapped slice — no copy, no extra allocation.
    let data = &mmap[110..];

    let half_diam = field_diameter * 0.5;
    let mut dec9_storage: i32 = 0; // high DEC byte, set by header records
    let mut current_mag = 0.0f64;
    let mut stars = Vec::new();
    let mut pos = 0;

    while pos + record_size <= data.len() {
        let ra7 = data[pos];
        let ra8 = data[pos + 1];
        let ra9 = data[pos + 2];
        let dec7 = data[pos + 3];
        let dec8 = data[pos + 4];
        pos += record_size;

        let ra_raw = (ra7 as u32) | ((ra8 as u32) << 8) | ((ra9 as u32) << 16);

        if ra_raw == HEADER_SENTINEL {
            // Magnitude header record
            current_mag = (dec8 as f64 - 16.0) / 10.0;
            dec9_storage = dec7 as i32 - 128; // signed byte
            continue;
        }

        // Normal star record
        let ra2 = ra_raw as f64 * RA_SCALE;
        let dec_raw = (dec9_storage << 16) | ((dec8 as i32) << 8) | (dec7 as i32);
        let dec2 = dec_raw as f64 * DEC_SCALE;

        // FOV filter (square window, not circular)
        let mut delta_ra = (ra2 - telescope_ra).abs();
        if delta_ra > PI {
            delta_ra = 2.0 * PI - delta_ra;
        }

        if delta_ra * cos_telescope_dec < half_diam && (dec2 - telescope_dec).abs() < half_diam {
            stars.push(CatalogStar {
                ra: ra2,
                dec: dec2,
                mag: current_mag,
            });
            if stars.len() >= max_stars {
                break;
            }
        }
    }

    Ok(stars)
}

/// Read catalog stars for a field of view by looking up the relevant areas.
///
/// - `db_path`: directory containing the catalog files
/// - `db_name`: catalog name prefix (e.g. `"d20"` for `d20_0101.1476`)
/// - `telescope_ra/dec`: pointing centre (radians)
/// - `fov`: square FOV side length (radians); capped to 5.14°
/// - `max_stars`: max stars to return in total
pub fn read_catalog_stars(
    db_path: &Path,
    db_name: &str,
    telescope_ra: f64,
    telescope_dec: f64,
    fov: f64,
    max_stars: usize,
) -> Result<Vec<CatalogStar>> {
    match detect_layout(db_path, db_name) {
        CatalogLayout::Areas290 => read_catalog_stars_290(
            db_path,
            db_name,
            telescope_ra,
            telescope_dec,
            fov,
            max_stars,
        ),
        CatalogLayout::AllSky001 => super::format_001::read_001_file(
            &db_path.join(format!("{db_name}_0101.001")),
            telescope_ra,
            telescope_dec,
            fov * 1.05,
            telescope_dec.cos().max(1e-6),
            max_stars,
        ),
        CatalogLayout::Areas1476 => read_catalog_stars_1476(
            db_path,
            db_name,
            telescope_ra,
            telescope_dec,
            fov,
            max_stars,
        ),
    }
}

/// Which sky tiling a database directory uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CatalogLayout {
    /// 1476 tiles across 36 equal-declination rings (D80, D50, D20, D05, V50, V05).
    Areas1476,
    /// 290 equal-area tiles across 18 rings (G05, for 3°–20° fields).
    Areas290,
    /// A single all-sky file of f32 triples (W08, for 20°–80° fields).
    AllSky001,
}

/// Detect the layout of `db_name` in `db_path` by probing for its first area file,
/// which is named `0101` in every format. Defaults to 1476 when none is present, so
/// a missing database still reports through the usual not-found path.
/// Whether `db_path` holds a database called `db_name` in any supported layout.
///
/// [`detect_layout`] cannot answer this: it falls back to [`CatalogLayout::Areas1476`]
/// when it recognises nothing, so a missing database is indistinguishable from a
/// 1476 one whose tiles are all absent. Callers check this first so a wrong `-d`
/// or `-D` reports "database not found" rather than reading nothing from every
/// spiral position and concluding the image is unsolvable.
pub fn catalog_present(db_path: &Path, db_name: &str) -> bool {
    ["1476", "290", "001"]
        .iter()
        .any(|ext| db_path.join(format!("{db_name}_0101.{ext}")).exists())
}

pub fn detect_layout(db_path: &Path, db_name: &str) -> CatalogLayout {
    if db_path.join(format!("{db_name}_0101.1476")).exists() {
        CatalogLayout::Areas1476
    } else if db_path.join(format!("{db_name}_0101.290")).exists() {
        CatalogLayout::Areas290
    } else if db_path.join(format!("{db_name}_0101.001")).exists() {
        CatalogLayout::AllSky001
    } else {
        CatalogLayout::Areas1476
    }
}

/// Read stars from a .290 database.
///
/// Two differences from the 1476 path, both forced by the field sizes these
/// databases exist for (G05 3°–20°, W08 20°–80°):
///
/// * every overlapping area is enumerated rather than sampling four corners, which
///   only works when the field fits inside one declination ring;
/// * `max_stars` is shared out across the areas and the combined list is then cut by
///   magnitude. Filling the budget area by area would take every star from the
///   southern edge of a 60° field and none from the north, because areas are visited
///   in declination order. Records within an area are already magnitude-ordered, so
///   a per-area slice is that area's brightest.
fn read_catalog_stars_290(
    db_path: &Path,
    db_name: &str,
    telescope_ra: f64,
    telescope_dec: f64,
    fov: f64,
    max_stars: usize,
) -> Result<Vec<CatalogStar>> {
    let areas = super::areas_290::find_areas_290(telescope_ra, telescope_dec, fov);
    if areas.is_empty() {
        return Ok(vec![]);
    }

    let field_diameter = fov * 1.05;
    let cos_dec = telescope_dec.cos().max(1e-6);
    // Read a little over the fair share so a sparse area does not starve the total.
    let per_area = (max_stars / areas.len()).max(16) * 2;

    let mut all_stars: Vec<CatalogStar> = Vec::new();
    for area_nr in areas {
        let fname = super::areas_290::filename_290(area_nr);
        let file_path = db_path.join(format!("{db_name}_{fname}"));
        match read_area_file(
            &file_path,
            telescope_ra,
            telescope_dec,
            field_diameter,
            cos_dec,
            per_area,
        ) {
            Ok(mut stars) => all_stars.append(&mut stars),
            Err(ArcsecError::CatalogIo(ref e)) if e.kind() == io::ErrorKind::NotFound => continue,
            Err(e) => return Err(e),
        }
    }

    if all_stars.len() > max_stars {
        all_stars.sort_unstable_by(|a, b| a.mag.total_cmp(&b.mag));
        all_stars.truncate(max_stars);
    }
    Ok(all_stars)
}

/// Read stars from a .1476 database. Unchanged ASTAP behaviour: four-corner area
/// sampling, FOV capped at one ring height, budget filled in area order.
fn read_catalog_stars_1476(
    db_path: &Path,
    db_name: &str,
    telescope_ra: f64,
    telescope_dec: f64,
    fov: f64,
    max_stars: usize,
) -> Result<Vec<CatalogStar>> {
    let areas = super::areas::find_areas_1476(telescope_ra, telescope_dec, fov);
    if areas.is_empty() {
        return Ok(vec![]);
    }

    // We read a slightly oversized area to ensure we cover the full FOV plus margins
    let field_diameter = fov * 1.05; // small margin
    let cos_dec = telescope_dec.cos().max(1e-6);

    let mut all_stars: Vec<CatalogStar> = Vec::new();

    for (area_nr, _frac) in areas {
        let fname = filename_1476(area_nr);
        let file_path = db_path.join(format!("{}_{}", db_name, fname));

        match read_area_file(
            &file_path,
            telescope_ra,
            telescope_dec,
            field_diameter,
            cos_dec,
            max_stars,
        ) {
            Ok(mut stars) => all_stars.append(&mut stars),
            Err(ArcsecError::CatalogIo(ref e)) if e.kind() == io::ErrorKind::NotFound => {
                // Area file not present — non-fatal (sparse databases)
                continue;
            }
            Err(e) => return Err(e),
        }

        if all_stars.len() >= max_stars {
            break;
        }
    }

    all_stars.truncate(max_stars);
    Ok(all_stars)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn deg(d: f64) -> f64 {
        d * PI / 180.0
    }

    /// Build a minimal synthetic .1476 file in memory.
    /// Emits a header record before each star so dec9_storage is always correct.
    /// Returns raw bytes: 110-byte file header + (header_record + star_record) pairs.
    fn make_synthetic_file(record_size: u8, stars: &[(f64, f64, f64)]) -> Vec<u8> {
        // File header: 110 bytes, last byte = record_size
        let mut buf = vec![0u8; 110];
        buf[109] = record_size;

        for &(ra, dec, mag) in stars {
            // Encode RA
            let ra_raw = (ra / (2.0 * PI) * 16_777_215.0).round() as u32;

            // Encode DEC (24-bit two's complement)
            let dec_raw = (dec / (PI * 0.5) * 8_388_607.0).round() as i32;
            let dec7 = (dec_raw & 0xFF) as u8;
            let dec8 = ((dec_raw >> 8) & 0xFF) as u8;
            let dec9: i8 = ((dec_raw >> 16) & 0xFF) as i8;

            // Header record (magnitude + dec9 high byte)
            let mag_enc = (mag * 10.0 + 16.0).round() as u8;
            buf.extend_from_slice(&[0xFF, 0xFF, 0xFF]);
            buf.push(dec9 as u8 + 128); // dec7 in header = dec9_storage + 128
            buf.push(mag_enc); // dec8 in header = mag + 16
            if record_size == 6 {
                buf.push(0);
            }

            // Star record
            buf.push((ra_raw & 0xFF) as u8);
            buf.push(((ra_raw >> 8) & 0xFF) as u8);
            buf.push(((ra_raw >> 16) & 0xFF) as u8);
            buf.push(dec7);
            buf.push(dec8);
            if record_size == 6 {
                buf.push(0);
            }
        }

        buf
    }

    #[test]
    fn decode_equatorial_star() {
        // Encode a star at RA=0, DEC=0
        let bytes = make_synthetic_file(5, &[(0.0, 0.0, 1.0)]);
        let path = std::env::temp_dir().join("test_decode.1476");
        std::fs::write(&path, &bytes).unwrap();

        let stars = read_area_file(&path, 0.0, 0.0, deg(10.0), 1.0, 100).unwrap();
        assert!(!stars.is_empty(), "should decode at least one star");
        let s = &stars[0];
        assert!(s.ra.abs() < 1e-4, "ra={}", s.ra);
        assert!(s.dec.abs() < 1e-4, "dec={}", s.dec);
        assert!((s.mag - 1.0).abs() < 0.1, "mag={}", s.mag);
    }

    #[test]
    fn header_record_not_returned_as_star() {
        // File with only a header record should return 0 stars
        let mut bytes = vec![0u8; 110];
        bytes[109] = 5; // record_size = 5
        // Header record
        bytes.extend_from_slice(&[0xFF, 0xFF, 0xFF, 128, 17]);
        let path = std::env::temp_dir().join("test_header_only.1476");
        std::fs::write(&path, &bytes).unwrap();

        let stars = read_area_file(&path, 0.0, 0.0, deg(10.0), 1.0, 100).unwrap();
        assert_eq!(
            stars.len(),
            0,
            "header record must not be returned as a star"
        );
    }

    #[test]
    fn fov_filter_excludes_distant_stars() {
        // Plant one star at (0,0) and one far away at (1 rad, 1 rad)
        let bytes = make_synthetic_file(5, &[(0.0, 0.0, 1.0), (1.0, 0.0, 2.0)]);
        let path = std::env::temp_dir().join("test_fov_filter.1476");
        std::fs::write(&path, &bytes).unwrap();

        // Only look in a narrow 1° FOV around (0,0)
        let stars = read_area_file(&path, 0.0, 0.0, deg(1.0), 1.0, 100).unwrap();
        assert_eq!(stars.len(), 1, "far star should be excluded");
    }

    #[test]
    fn ra_roundtrip_accuracy() {
        // Test RA encoding/decoding round-trip accuracy
        let test_ra = deg(123.456);
        let bytes = make_synthetic_file(5, &[(test_ra, 0.0, 1.0)]);
        let path = std::env::temp_dir().join("test_ra_roundtrip.1476");
        std::fs::write(&path, &bytes).unwrap();

        let stars = read_area_file(&path, test_ra, 0.0, deg(5.0), 1.0, 100).unwrap();
        assert!(!stars.is_empty(), "star not found");
        // RA resolution is ~0.077 arcsec = 3.7e-7 rad
        assert!(
            (stars[0].ra - test_ra).abs() < 1e-5,
            "ra error = {}",
            stars[0].ra - test_ra
        );
    }

    #[test]
    fn record_size_6_works() {
        // 6-byte record format (with Gaia colour byte)
        let bytes = make_synthetic_file(6, &[(0.5, 0.1, 2.0)]);
        let path = std::env::temp_dir().join("test_record6.1476");
        std::fs::write(&path, &bytes).unwrap();

        let stars = read_area_file(&path, 0.5, 0.1, deg(5.0), 0.995_f64.cos(), 100).unwrap();
        assert!(!stars.is_empty(), "should decode record_size=6 star");
    }

    #[test]
    fn missing_file_returns_empty_catalog() {
        let db_path = std::env::temp_dir();
        let result = read_catalog_stars(&db_path, "nonexistent_db", 0.0, 0.0, deg(2.0), 100);
        // With every area file missing, the read must report no stars rather than
        // inventing any; a propagated I/O error is equally acceptable.
        if let Ok(stars) = result {
            assert!(
                stars.is_empty(),
                "missing database yielded {} stars",
                stars.len()
            );
        }
    }
}
