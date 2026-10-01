//! Binary reader for ASTAP's packed Gaia catalog formats, .1476 and .290.
//!
//! Both carry an identical 110-byte ASCII header (byte 109 = record size) followed by
//! the same packed records, so this reader handles both; only the sky tiling differs
//! (`areas.rs` for 1476, `areas_290.rs` for 290).
//!
//! Record format (5 bytes; a 6th, Gaia colour, is present in some databases and
//! ignored here):
//!
//! ```text
//! ra7, ra8, ra9  — 24-bit LE unsigned integer for RA
//! dec7, dec8     — low 2 bytes of a 24-bit two's-complement DEC integer
//! ```
//!
//! Header records have `ra_raw == 0xFF_FF_FF` and carry the high DEC byte
//! (`dec9_storage`) and magnitude for the following group of star records.
//!
//! ```text
//! RA  = ra_raw / (2^24 - 1) * 2π
//! DEC = signed_24bit_int(dec9_storage:dec8:dec7) / (2^23 - 1) * (π/2)
//! ```

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
    /// Right ascension, radians.
    pub ra: f64,
    /// Declination, radians.
    pub dec: f64,
    /// Magnitude (Gaia BP for the ASTAP databases).
    pub mag: f64,
}

/// Read all stars from a single .1476 area file that fall within a square FOV.
///
/// - `file_path`: full path to the area file (e.g. `/db/d20_0101.1476`)
/// - `telescope_ra/dec`: centre of the field (radians)
/// - `field_diameter`: diameter of the square area to collect (radians)
/// - `cos_telescope_dec`: pre-computed `cos(telescope_dec)` for fast RA delta check
/// - `max_stars`: maximum stars to return (stops reading after this many)
///
/// # Errors
///
/// [`ArcsecError::CatalogIo`] if the file cannot be opened or mapped (a missing
/// file is `io::ErrorKind::NotFound`), or its header declares an unsupported
/// record size.
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
    let record_size = record_size(&mmap)?;

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
///
/// Missing area files are skipped silently (sparse databases are normal).
///
/// # Errors
///
/// [`ArcsecError::CatalogIo`] if an area file exists but cannot be read or is
/// malformed.
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

/// Call `f` for every star of `db_name` brighter than or equal to `mag_limit` whose
/// declination lies in `[dec_lo, dec_hi]` (radians), all the way round in RA.
///
/// This is the whole-sky read the blind-index builder needs, a declination band at
/// a time so that a deep magnitude limit never holds the whole catalogue in memory.
/// Every area file is sorted brightest first, so each is read only down to
/// `mag_limit`. Stars arrive in no particular order. Works on all three layouts.
///
/// # Errors
///
/// [`ArcsecError::CatalogIo`] if an area file exists but cannot be read or is
/// malformed. Missing area files are skipped, as everywhere else.
pub fn for_each_star_in_dec_band(
    db_path: &Path,
    db_name: &str,
    dec_lo: f64,
    dec_hi: f64,
    mag_limit: f64,
    mut f: impl FnMut(&CatalogStar),
) -> Result<()> {
    let paths: Vec<std::path::PathBuf> = match detect_layout(db_path, db_name) {
        CatalogLayout::AllSky001 => {
            // Brightest first and small (W08 is ~41k stars): read whole.
            let all = super::format_001::read_001_file(
                &db_path.join(format!("{db_name}_0101.001")),
                0.0,
                0.0,
                4.0 * PI,
                1.0,
                usize::MAX,
            )?;
            for s in all.iter().filter(|s| s.dec >= dec_lo && s.dec <= dec_hi) {
                if s.mag > mag_limit {
                    break;
                }
                f(s);
            }
            return Ok(());
        }
        CatalogLayout::Areas290 => super::areas_290::areas_in_dec_band_290(dec_lo, dec_hi)
            .into_iter()
            .map(|a| db_path.join(format!("{db_name}_{}", super::areas_290::filename_290(a))))
            .collect(),
        CatalogLayout::Areas1476 => super::areas::areas_in_dec_band_1476(dec_lo, dec_hi)
            .into_iter()
            .map(|a| db_path.join(format!("{db_name}_{}", filename_1476(a))))
            .collect(),
    };
    for path in paths {
        let cursor = match AreaCursor::open(&path) {
            Ok(Some(c)) => c,
            Ok(None) => continue,
            Err(ArcsecError::CatalogIo(ref e)) if e.kind() == io::ErrorKind::NotFound => continue,
            Err(e) => return Err(e),
        };
        let data = &cursor.mmap[..];
        let rs = cursor.record_size;
        let mut pos = cursor.pos;
        let mut dec9 = 0i32;
        let mut mag = 0.0f64;
        while pos + rs <= data.len() {
            let r = &data[pos..pos + 5];
            pos += rs;
            let ra_raw = (r[0] as u32) | ((r[1] as u32) << 8) | ((r[2] as u32) << 16);
            if ra_raw == HEADER_SENTINEL {
                mag = header_mag(r[4]);
                if mag > mag_limit {
                    break;
                }
                dec9 = r[3] as i32 - 128;
                continue;
            }
            let dec_raw = (dec9 << 16) | ((r[4] as i32) << 8) | (r[3] as i32);
            let dec = dec_raw as f64 * DEC_SCALE;
            if dec >= dec_lo && dec <= dec_hi {
                f(&CatalogStar {
                    ra: ra_raw as f64 * RA_SCALE,
                    dec,
                    mag,
                });
            }
        }
    }
    Ok(())
}

/// Which sky tiling a database directory uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CatalogLayout {
    /// 1476 tiles across 36 equal-declination rings (D80, D50, D20, D05, V50).
    Areas1476,
    /// 290 equal-area tiles across 18 rings (G05 for 3°–20° fields; also V05).
    Areas290,
    /// A single all-sky file of f32 triples (W08, for 20°–80° fields).
    AllSky001,
}

/// Whether `db_path` holds a database called `db_name` in any supported layout.
///
/// [`detect_layout`] cannot answer this: it falls back to [`CatalogLayout::Areas1476`]
/// when it recognises nothing, so a missing database is indistinguishable from a
/// 1476 one whose tiles are all absent. Callers check this first so a wrong `-d`
/// or `-D` reports "database not found" rather than reading nothing from every
/// spiral position and concluding the image is unsolvable.
#[must_use]
pub fn catalog_present(db_path: &Path, db_name: &str) -> bool {
    ["1476", "290", "001"]
        .iter()
        .any(|ext| db_path.join(format!("{db_name}_0101.{ext}")).exists())
}

/// Detect the layout of `db_name` in `db_path`.
///
/// Probes for its first area file, which is named `0101` in every format. Defaults
/// to 1476 when none is present, so a missing database still reports through the
/// usual not-found path.
#[must_use]
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
/// Every overlapping area is enumerated rather than sampling four corners, which
/// only works when the field fits inside one declination ring; the databases exist
/// for 3°–80° fields, which span many rings. The areas are then read together, a
/// magnitude step at a time, by [`read_brightest`].
fn read_catalog_stars_290(
    db_path: &Path,
    db_name: &str,
    telescope_ra: f64,
    telescope_dec: f64,
    fov: f64,
    max_stars: usize,
) -> Result<Vec<CatalogStar>> {
    let areas = super::areas_290::find_areas_290(telescope_ra, telescope_dec, fov);
    read_brightest(
        areas.into_iter().map(|area_nr| {
            let fname = super::areas_290::filename_290(area_nr);
            db_path.join(format!("{db_name}_{fname}"))
        }),
        telescope_ra,
        telescope_dec,
        fov,
        max_stars,
    )
}

/// Read stars from a .1476 database: ASTAP's four-corner area sampling, FOV capped
/// at one ring height, and the (up to four) areas read together by
/// [`read_brightest`].
fn read_catalog_stars_1476(
    db_path: &Path,
    db_name: &str,
    telescope_ra: f64,
    telescope_dec: f64,
    fov: f64,
    max_stars: usize,
) -> Result<Vec<CatalogStar>> {
    let areas = super::areas::find_areas_1476(telescope_ra, telescope_dec, fov);
    read_brightest(
        areas.into_iter().map(|(area_nr, _frac)| {
            let fname = filename_1476(area_nr);
            db_path.join(format!("{db_name}_{fname}"))
        }),
        telescope_ra,
        telescope_dec,
        fov,
        max_stars,
    )
}

/// The brightest `max_stars` stars in the square window of side `fov * 1.05` round
/// the pointing, drawn from every area file in `paths`, whichever of them holds
/// them, and returned brightest first.
///
/// Each area file is brightest first, in groups of one magnitude step (0.1), so
/// the files are read in step: the brightest unread group of every file, then the
/// next, until the window holds `max_stars` stars. Each file is read only as deep
/// as the field's own magnitude limit, so a file that the field barely clips costs
/// no more than one the field sits inside.
///
/// This replaces filling the budget file by file, which took every star from the
/// first file and none from the rest, so a field straddling a tile boundary had
/// catalogue stars on one side only (plate-solving.md §11.11). ASTAP shares the
/// budget in proportion to each area's share of the field instead, which gives the
/// brightest stars only where the sky is uniformly dense; reading by magnitude
/// gives them everywhere.
///
/// The cut usually falls part-way through the last magnitude group read. The stars
/// kept from that group are the first read: from the files in the order given, and
/// within a file south to north (records in a group are ordered by their high Dec
/// byte). With a single file that is exactly the old, and ASTAP's, choice. Choosing
/// them by a hash of position instead, so the cut is spatially even, made no
/// measurable difference to the benchmark and was not kept (test-images.md §7.7).
///
/// Missing area files are skipped (sparse databases are normal); one that exists
/// but cannot be read is an error.
fn read_brightest(
    paths: impl Iterator<Item = std::path::PathBuf>,
    telescope_ra: f64,
    telescope_dec: f64,
    fov: f64,
    max_stars: usize,
) -> Result<Vec<CatalogStar>> {
    let window = Window {
        ra: telescope_ra,
        dec: telescope_dec,
        // A slightly oversized area, to cover the full FOV plus margins.
        half: fov * 1.05 * 0.5,
        cos_dec: telescope_dec.cos().max(1e-6),
    };

    let mut cursors = Vec::new();
    for path in paths {
        match AreaCursor::open(&path) {
            Ok(Some(c)) => cursors.push(c),
            Ok(None) => {}
            Err(ArcsecError::CatalogIo(ref e)) if e.kind() == io::ErrorKind::NotFound => {}
            Err(e) => return Err(e),
        }
    }

    let mut stars: Vec<CatalogStar> = Vec::new();
    while stars.len() < max_stars {
        let Some(group) = cursors
            .iter()
            .filter_map(AreaCursor::next_mag)
            .min_by(f64::total_cmp)
        else {
            break; // every file exhausted
        };
        for c in &mut cursors {
            c.read_through(group, &window, &mut stars);
        }
    }

    // Stable, so stars of the magnitude the cut falls in keep their read order.
    stars.sort_by(|a, b| a.mag.total_cmp(&b.mag));
    stars.truncate(max_stars);
    Ok(stars)
}

/// The square field window the readers collect: centre, half side and `cos(dec)`,
/// all radians.
struct Window {
    ra: f64,
    dec: f64,
    half: f64,
    cos_dec: f64,
}

impl Window {
    fn contains(&self, ra: f64, dec: f64) -> bool {
        let mut delta_ra = (ra - self.ra).abs();
        if delta_ra > PI {
            delta_ra = 2.0 * PI - delta_ra;
        }
        delta_ra * self.cos_dec < self.half && (dec - self.dec).abs() < self.half
    }
}

/// A read position in one memory-mapped area file.
struct AreaCursor {
    mmap: Mmap,
    record_size: usize,
    /// Byte offset of the next unread record.
    pos: usize,
    /// High Dec byte and magnitude from the last header record.
    dec9: i32,
    mag: f64,
}

impl AreaCursor {
    /// Map `path` and check its header. `None` for a file too short to hold one.
    fn open(path: &Path) -> Result<Option<Self>> {
        let file = std::fs::File::open(path).map_err(ArcsecError::CatalogIo)?;
        // Safety: we only read; the catalog files are never written concurrently.
        let mmap = unsafe { Mmap::map(&file).map_err(ArcsecError::CatalogIo)? };
        if mmap.len() < 110 {
            return Ok(None);
        }
        let record_size = record_size(&mmap)?;
        Ok(Some(Self {
            mmap,
            record_size,
            pos: 110,
            dec9: 0,
            mag: 0.0,
        }))
    }

    /// The magnitude of the next unread record: its own if it is a header,
    /// otherwise that of the group it continues. `None` at the end of the file.
    fn next_mag(&self) -> Option<f64> {
        let r = self.mmap.get(self.pos..self.pos + self.record_size)?;
        Some(if r[..3] == [0xFF; 3] {
            header_mag(r[4])
        } else {
            self.mag
        })
    }

    /// Read every record up to and including magnitude `limit`, appending the stars
    /// inside `window` to `out`. Stops before the first header fainter than `limit`.
    fn read_through(&mut self, limit: f64, window: &Window, out: &mut Vec<CatalogStar>) {
        let data = &self.mmap[..];
        let rs = self.record_size;
        while self.pos + rs <= data.len() {
            let r = &data[self.pos..self.pos + 5];
            let ra_raw = (r[0] as u32) | ((r[1] as u32) << 8) | ((r[2] as u32) << 16);
            if ra_raw == HEADER_SENTINEL {
                let mag = header_mag(r[4]);
                if mag > limit {
                    return;
                }
                self.mag = mag;
                self.dec9 = r[3] as i32 - 128; // signed byte
            } else {
                let ra = ra_raw as f64 * RA_SCALE;
                let dec_raw = (self.dec9 << 16) | ((r[4] as i32) << 8) | (r[3] as i32);
                let dec = dec_raw as f64 * DEC_SCALE;
                if window.contains(ra, dec) {
                    out.push(CatalogStar {
                        ra,
                        dec,
                        mag: self.mag,
                    });
                }
            }
            self.pos += rs;
        }
    }
}

/// Magnitude carried by a header record's last byte.
fn header_mag(byte: u8) -> f64 {
    (byte as f64 - 16.0) / 10.0
}

/// Record size from a mapped area file's 110-byte header (byte 109): 5 (RA + Dec)
/// or 6 (plus Gaia colour); anything else is unsupported.
fn record_size(mmap: &[u8]) -> Result<usize> {
    let record_size = if mmap[109] == b' ' {
        11
    } else {
        mmap[109] as usize
    };
    if record_size != 5 && record_size != 6 {
        return Err(ArcsecError::CatalogIo(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("unsupported record_size {record_size}"),
        )));
    }
    Ok(record_size)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn deg(d: f64) -> f64 {
        d * PI / 180.0
    }

    /// Build a minimal synthetic .1476 file in memory.
    /// Emits a header record before each star so `dec9_storage` is always correct.
    /// Returns raw bytes: 110-byte file header + (`header_record` + `star_record`) pairs.
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

    // ── Record decoding, using the ASTAP-layout writer ─────────────────────────

    use crate::test_support::{
        Rng, SkySpec, SkyStar, TempDir, area_file_bytes, random_sky, separation, write_001_db,
        write_290_db, write_1476_db,
    };

    fn write_area(dir: &TempDir, name: &str, bytes: &[u8]) -> std::path::PathBuf {
        let path = dir.path().join(name);
        std::fs::write(&path, bytes).unwrap();
        path
    }

    /// Stars across the whole declination range, both RA edges and a spread of
    /// magnitudes, so every high-DEC-byte value and header-record transition is hit.
    fn all_sky_sample() -> Vec<SkyStar> {
        let mut out = Vec::new();
        for (k, dec_deg) in [-89.99, -60.0, -30.5, -0.001, 0.0, 0.001, 29.0, 61.3, 89.99]
            .iter()
            .enumerate()
        {
            for (j, ra_deg) in [0.0, 0.0001, 123.456, 359.9999].iter().enumerate() {
                // A distinct magnitude per star, so each can be found again by it.
                out.push(SkyStar {
                    ra: deg(*ra_deg),
                    dec: deg(*dec_deg),
                    mag: -1.5 + 0.3 * (4 * k + j) as f64,
                });
            }
        }
        out
    }

    #[test]
    fn every_declination_and_magnitude_round_trips() {
        let dir = TempDir::new("f1476-rt");
        for record_size in [5usize, 6] {
            let stars = all_sky_sample();
            let path = write_area(&dir, "all.1476", &area_file_bytes(&stars, record_size));
            // A field as wide as the sky, so the filter passes everything.
            let got = read_area_file(&path, PI, 0.0, 4.0 * PI, 1.0, usize::MAX).unwrap();
            assert_eq!(got.len(), stars.len(), "record size {record_size}");
            for want in &stars {
                let found = got
                    .iter()
                    .find(|g| (g.mag - want.mag).abs() < 0.051)
                    .unwrap_or_else(|| panic!("lost {want:?}"));
                // Packing resolution: 2π/2²⁴ in RA, (π/2)/2²³ in Dec, 0.1 in magnitude.
                assert!((found.dec - want.dec).abs() < 2e-7, "{found:?} vs {want:?}");
                let dra = (found.ra - want.ra).rem_euclid(2.0 * PI);
                assert!(dra.min(2.0 * PI - dra) < 4e-7, "{found:?} vs {want:?}");
                assert!((0.0..2.0 * PI).contains(&found.ra));
            }
        }
    }

    #[test]
    fn stars_come_back_brightest_first_and_max_stars_truncates() {
        let dir = TempDir::new("f1476-max");
        let mut rng = Rng::new(5);
        let stars: Vec<SkyStar> = (0..50)
            .map(|_| SkyStar {
                ra: rng.range(1.0, 1.01),
                dec: rng.range(0.2, 0.21),
                mag: rng.range(5.0, 15.0),
            })
            .collect();
        let path = write_area(&dir, "a.1476", &area_file_bytes(&stars, 5));
        let got = read_area_file(&path, 1.005, 0.205, deg(5.0), 0.205f64.cos(), 10).unwrap();
        assert_eq!(got.len(), 10);
        assert!(got.windows(2).all(|w| w[0].mag <= w[1].mag));
        let mut mags: Vec<f64> = stars.iter().map(|s| s.mag).collect();
        mags.sort_by(f64::total_cmp);
        assert!((got[9].mag - mags[9]).abs() < 0.051);
    }

    #[test]
    fn the_square_window_wraps_at_ra_zero_and_scales_with_declination() {
        let dir = TempDir::new("f1476-win");
        let stars = [
            SkyStar {
                ra: deg(359.8),
                dec: deg(60.0),
                mag: 5.0,
            },
            SkyStar {
                ra: deg(0.3),
                dec: deg(60.0),
                mag: 5.0,
            },
            SkyStar {
                ra: deg(2.5),
                dec: deg(60.0),
                mag: 5.0,
            }, // 1.25° on the sky
            SkyStar {
                ra: deg(0.0),
                dec: deg(61.1),
                mag: 5.0,
            }, // outside in Dec
        ];
        let path = write_area(&dir, "w.1476", &area_file_bytes(&stars, 5));
        // A 2°-wide field at Dec 60: ±1° in Dec, ±2° in RA.
        let got = read_area_file(&path, 0.0, deg(60.0), deg(2.0), deg(60.0).cos(), 100).unwrap();
        assert_eq!(got.len(), 2, "{got:?}");
    }

    #[test]
    fn short_and_malformed_area_files() {
        let dir = TempDir::new("f1476-bad");
        // Missing: NotFound, which the database readers treat as "no tile here".
        match read_area_file(&dir.path().join("x.1476"), 0.0, 0.0, 1.0, 1.0, 10) {
            Err(ArcsecError::CatalogIo(e)) => assert_eq!(e.kind(), io::ErrorKind::NotFound),
            other => panic!("{other:?}"),
        }
        // Shorter than the header: no stars, not an error.
        let short = write_area(&dir, "short.1476", &[b' '; 60]);
        assert!(
            read_area_file(&short, 0.0, 0.0, 1.0, 1.0, 10)
                .unwrap()
                .is_empty()
        );
        // Unsupported record sizes, including the legacy ' ' (11-byte) marker.
        for marker in [b' ', 7, 0] {
            let mut bytes = vec![b' '; 110];
            bytes[109] = marker;
            bytes.extend_from_slice(&[0; 22]);
            let path = write_area(&dir, "odd.1476", &bytes);
            match read_area_file(&path, 0.0, 0.0, 1.0, 1.0, 10) {
                Err(ArcsecError::CatalogIo(e)) => {
                    assert_eq!(e.kind(), io::ErrorKind::InvalidData);
                }
                other => panic!("marker {marker}: {other:?}"),
            }
        }
        // A truncated final record is ignored, not read past the end.
        let mut bytes = area_file_bytes(
            &[SkyStar {
                ra: 0.1,
                dec: 0.1,
                mag: 3.0,
            }],
            5,
        );
        bytes.extend_from_slice(&[1, 2, 3]);
        let path = write_area(&dir, "trunc.1476", &bytes);
        assert_eq!(
            read_area_file(&path, 0.1, 0.1, 0.1, 1.0, 10).unwrap().len(),
            1
        );
    }

    // ── Database-level reads ────────────────────────────────────────────────────

    fn field(seed: u64, ra: f64, dec: f64, side_deg: f64, n: usize) -> Vec<SkyStar> {
        random_sky(
            &mut Rng::new(seed),
            &SkySpec {
                ra0: deg(ra),
                dec0: deg(dec),
                side_deg,
                n,
                min_sep_deg: 0.0,
                mag_lo: 6.0,
                mag_hi: 16.0,
            },
        )
    }

    /// Whether a star lies in the readers' square window of side `side`.
    fn in_window(ra0: f64, dec0: f64, side: f64, ra: f64, dec: f64) -> bool {
        let half = side * 0.5;
        let mut dra = (ra - ra0).abs();
        if dra > PI {
            dra = 2.0 * PI - dra;
        }
        dra * dec0.cos() < half && (dec - dec0).abs() < half
    }

    /// Every star of the field proper comes back, once, and nothing from outside
    /// the reader's 5% margin. (Stars *inside* the margin are returned only when
    /// their tile overlaps the unmargined field, because both tile finders are
    /// given `fov` rather than `fov * 1.05`; so the margin is best-effort.)
    fn assert_field_read(all: &[SkyStar], got: &[CatalogStar], ra: f64, dec: f64, fov: f64) {
        for s in all.iter().filter(|s| in_window(ra, dec, fov, s.ra, s.dec)) {
            let n = got
                .iter()
                .filter(|g| separation(g.ra, g.dec, s.ra, s.dec) < 1e-6)
                .count();
            assert_eq!(n, 1, "{s:?} returned {n} times");
        }
        for g in got {
            assert!(in_window(ra, dec, fov * 1.05, g.ra, g.dec), "{g:?} outside");
        }
    }

    #[test]
    fn layout_detection_and_presence() {
        let dir = TempDir::new("layout");
        assert!(!catalog_present(dir.path(), "d50"));
        assert_eq!(detect_layout(dir.path(), "d50"), CatalogLayout::Areas1476);
        write_290_db(dir.path(), "g05", &[]);
        write_001_db(dir.path(), "w08", &[]);
        write_1476_db(dir.path(), "d50", &[]);
        assert_eq!(detect_layout(dir.path(), "g05"), CatalogLayout::Areas290);
        assert_eq!(detect_layout(dir.path(), "w08"), CatalogLayout::AllSky001);
        assert_eq!(detect_layout(dir.path(), "d50"), CatalogLayout::Areas1476);
        for name in ["g05", "w08", "d50"] {
            assert!(catalog_present(dir.path(), name), "{name}");
        }
        assert!(!catalog_present(dir.path(), "d80"));
    }

    #[test]
    fn a_1476_database_returns_the_whole_field_across_tile_boundaries() {
        // RA 0 and a declination ring boundary (5.14°) both run through this field.
        let (ra, dec, fov) = (0.0, deg(5.2), deg(2.0));
        let all = field(1, 0.0, 5.2, 4.0, 3000);
        let dir = TempDir::new("db1476");
        write_1476_db(dir.path(), "d50", &all);
        let got = read_catalog_stars(dir.path(), "d50", ra, dec, fov, usize::MAX).unwrap();
        assert_field_read(&all, &got, ra, dec, fov);
        // With a budget, it stops there.
        let capped = read_catalog_stars(dir.path(), "d50", ra, dec, fov, 25).unwrap();
        assert_eq!(capped.len(), 25);
    }

    /// A 12° G05-style field at Dec -30 with a 290 database under it. It touches
    /// six tiles, but tile 75 holds 76 of the field's 100 brightest stars.
    fn g05_field() -> (TempDir, Vec<CatalogStar>, (f64, f64, f64)) {
        let (ra, dec, fov) = (deg(200.0), deg(-30.0), deg(12.0));
        let all = field(2, 200.0, -30.0, 20.0, 4000);
        let dir = TempDir::new("db290");
        write_290_db(dir.path(), "g05", &all);
        let full = read_catalog_stars(dir.path(), "g05", ra, dec, fov, usize::MAX).unwrap();
        assert_field_read(&all, &full, ra, dec, fov);
        (dir, full, (ra, dec, fov))
    }

    #[test]
    fn a_290_database_reads_every_tile_and_keeps_to_its_budget() {
        let (dir, _, (ra, dec, fov)) = g05_field();
        let got = read_catalog_stars(dir.path(), "g05", ra, dec, fov, 100).unwrap();
        assert_eq!(got.len(), 100);
        assert!(got.windows(2).all(|w| w[0].mag <= w[1].mag), "sorted");
        // Not filled tile by tile from the south: both ends of the field appear.
        assert!(got.iter().any(|s| s.dec < dec - deg(3.0)));
        assert!(got.iter().any(|s| s.dec > dec + deg(3.0)));
    }

    /// Tile 75 holds 76 of the field's 100 brightest stars. The reader used to give
    /// each tile a fixed share of the budget, `(max_stars / n_tiles).max(16) * 2`,
    /// so tile 75 contributed only 32 and the tiles clipping the field's edges filled
    /// the rest with fainter stars (a cut at mag 10.8 where the field's 100th
    /// brightest is mag 9.0): a catalogue dense at the edges and sparse in the middle.
    #[test]
    fn a_290_database_returns_the_brightest_stars_of_the_field() {
        let (dir, full, (ra, dec, fov)) = g05_field();
        let got = read_catalog_stars(dir.path(), "g05", ra, dec, fov, 100).unwrap();
        let mut mags: Vec<f64> = full.iter().map(|s| s.mag).collect();
        mags.sort_by(f64::total_cmp);
        let faintest = got.iter().map(|s| s.mag).fold(f64::MIN, f64::max);
        assert!(
            faintest <= mags[99] + 0.05,
            "cut at mag {faintest}, but the field's 100th brightest is {}",
            mags[99]
        );
        assert_brightest(&full, &got, 100);
    }

    /// `got` is `n` stars of `full` (a whole-field read), brightest first, and holds
    /// every star of `full` strictly brighter than the faintest it returns: the
    /// field's `n` brightest, whatever tiles they came from.
    fn assert_brightest(full: &[CatalogStar], got: &[CatalogStar], n: usize) {
        assert_eq!(got.len(), n);
        assert!(
            got.windows(2).all(|w| w[0].mag <= w[1].mag),
            "brightest first"
        );
        let cut = got[n - 1].mag;
        for s in full.iter().filter(|s| s.mag < cut) {
            assert!(
                got.iter()
                    .any(|g| separation(g.ra, g.dec, s.ra, s.dec) < 1e-9),
                "{s:?} is brighter than the cut at {cut} but was not returned"
            );
        }
        for g in got {
            assert!(
                full.iter()
                    .any(|s| separation(g.ra, g.dec, s.ra, s.dec) < 1e-9),
                "{g:?} is not in the field"
            );
        }
    }

    /// A field whose four corners fall in four different 1476 tiles (RA 0 and the
    /// 5.14° ring boundary cross in it). The budget used to be filled from the first
    /// tile alone, leaving the other three quarters of the field without catalogue
    /// stars (plate-solving.md §11.11); now it is the field's brightest, wherever
    /// they are.
    #[test]
    fn a_1476_budget_is_shared_across_every_tile_of_the_field() {
        let (ra, dec, fov) = (deg(0.3), deg(5.0), deg(2.0));
        assert_eq!(super::super::areas::find_areas_1476(ra, dec, fov).len(), 4);
        let all = field(5, 0.3, 5.0, 4.0, 4000);
        let dir = TempDir::new("db1476-budget");
        write_1476_db(dir.path(), "d50", &all);
        let full = read_catalog_stars(dir.path(), "d50", ra, dec, fov, usize::MAX).unwrap();
        assert_field_read(&all, &full, ra, dec, fov);
        for n in [25, 200, 600] {
            let got = read_catalog_stars(dir.path(), "d50", ra, dec, fov, n).unwrap();
            assert_brightest(&full, &got, n);
        }
        // Every quarter of the field is represented in proportion to its area.
        let got = read_catalog_stars(dir.path(), "d50", ra, dec, fov, 400).unwrap();
        let ring = deg(5.142_857);
        for (east, north) in [(false, false), (false, true), (true, false), (true, true)] {
            let n = got
                .iter()
                .filter(|s| (s.ra < PI) == east && (s.dec > ring) == north)
                .count();
            assert!(n > 20, "quarter east={east} north={north} has {n} of 400");
        }
    }

    #[test]
    fn an_all_sky_001_database_is_read_through_the_same_entry_point() {
        let (ra, dec, fov) = (deg(90.0), deg(-50.0), deg(30.0));
        let all = field(3, 90.0, -50.0, 60.0, 2000);
        let dir = TempDir::new("db001");
        write_001_db(dir.path(), "w08", &all);
        let got = read_catalog_stars(dir.path(), "w08", ra, dec, fov, usize::MAX).unwrap();
        assert_field_read(&all, &got, ra, dec, fov);
    }

    #[test]
    fn a_corrupt_tile_is_an_error_but_a_missing_one_is_not() {
        let (ra, dec, fov) = (deg(40.0), deg(20.0), deg(1.0));
        let all = field(4, 40.0, 20.0, 2.0, 500);
        for layout in ["1476", "290"] {
            let dir = TempDir::new("corrupt");
            if layout == "1476" {
                write_1476_db(dir.path(), "db", &all);
            } else {
                write_290_db(dir.path(), "db", &all);
            }
            // A clean read first; then corrupt every tile, then delete them.
            assert!(
                !read_catalog_stars(dir.path(), "db", ra, dec, fov, 1000)
                    .unwrap()
                    .is_empty()
            );
            let tiles: Vec<_> = std::fs::read_dir(dir.path())
                .unwrap()
                .map(|e| e.unwrap().path())
                .filter(|p| !p.to_string_lossy().contains("_0101."))
                .collect();
            for t in &tiles {
                let mut b = std::fs::read(t).unwrap();
                b[109] = 9;
                std::fs::write(t, b).unwrap();
            }
            assert!(
                matches!(
                    read_catalog_stars(dir.path(), "db", ra, dec, fov, 1000),
                    Err(ArcsecError::CatalogIo(_))
                ),
                "{layout}"
            );
            for t in &tiles {
                std::fs::remove_file(t).unwrap();
            }
            assert!(
                read_catalog_stars(dir.path(), "db", ra, dec, fov, 1000)
                    .unwrap()
                    .is_empty()
            );
        }
    }
}
