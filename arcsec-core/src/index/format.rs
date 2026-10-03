//! The on-disk blind index: `ARCSECIX` version 1.
//!
//! One file, memory-mapped, little-endian throughout, every section 8-byte aligned.
//! Opening it checks the header (magic, version, byte-order marker, its own CRC) and
//! that every section lies inside the file with the size its counts imply — a few
//! hundred bytes of work, so a 1 GB index opens instantly. The section CRCs are
//! checked only by [`BlindIndex::validate`] (`arcsec catalog verify`); lookups
//! bounds-check every star reference they follow, so a corrupt body can make a solve
//! fail but cannot make it read out of bounds.
//!
//! ```text
//! header, 256 bytes
//!   0  magic        [u8; 8]  "ARCSECIX"
//!   8  version      u32      1
//!  12  byte order   u32      0x0A0B0C0D as written (reads 0x0D0C0B0A if swapped)
//!  16  header_len   u32      256
//!  20  bins         u32      descriptor bins per dimension (128)
//!  24  n_tiers      u32
//!  28  star_bands   u32      declination bands in the star directory
//!  32  n_stars      u64
//!  40  n_patterns   u64
//!  48  built_unix   i64
//!  56  source       [u8; 16] database the stars came from, NUL-padded ("d80")
//!  72  sections     5 × { offset u64, length u64, crc32 u32, reserved u32 }
//!                   tiers, stars, star directory, keys, quads
//! 192  source hash  u64      fingerprint of the source database's files (see
//!                            [`SourceStamp`]); 0 if not recorded
//! 200  source bytes u64      total size of those files
//! 208  source files u32      how many there were; 0 = no stamp recorded
//! 212  reserved     zero
//! 252  header crc32 over bytes 0..252
//!
//! tiers        n_tiers × 40 bytes, widest first:
//!                radius f64 (rad), mag_cap f32, members u32,
//!                first_pattern u64, n_patterns u64, n_anchors u64
//! stars        n_stars × 12 bytes, sorted by (declination band, RA):
//!                ra f32 (rad), dec f32 (rad), mag i16 (×100), tier u8, 0 u8
//! star dir     (star_bands + 1) × u32: first star of each declination band
//! keys         n_patterns × u64, sorted ascending within each tier
//! quads        n_patterns × 4 × u32: star indices in canonical vertex order
//! ```
//!
//! A tier's patterns are the contiguous range `first_pattern .. +n_patterns` of
//! both `keys` and `quads`.
//!
//! The source stamp (bytes 192–211) was added after the first release of version 1,
//! in space that was reserved and written as zero, so files without it are still
//! version 1 and read as "not recorded"; older readers ignore it.

use core::ops::Range;
use std::fs::File;
use std::io::{self, BufWriter, Seek as _, Write};
use std::path::Path;

use memmap2::Mmap;

use crate::error::{ArcsecError, Result};

/// File magic.
pub const MAGIC: &[u8; 8] = b"ARCSECIX";
/// Format version this code reads and writes.
pub const VERSION: u32 = 1;
/// Conventional file name in the catalogue directory: `<db>.arcsecix`.
pub const EXTENSION: &str = "arcsecix";

const HEADER_LEN: usize = 256;
const BYTE_ORDER: u32 = 0x0A0B_0C0D;
const N_SECTIONS: usize = 5;
const SECTIONS_AT: usize = 72;
const HEADER_CRC_AT: usize = 252;
const STAMP_AT: usize = 192;
const TIER_LEN: usize = 40;
const STAR_LEN: usize = 12;
const QUAD_LEN: usize = 16;

/// One scale tier of the index.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TierInfo {
    /// Disc radius, radians: patterns are drawn from stars within this of an anchor,
    /// so no pattern is wider than twice it.
    pub radius: f64,
    /// Faintest magnitude considered for this tier.
    pub mag_cap: f32,
    /// Stars per anchor group (anchor included).
    pub members: u32,
    /// First pattern of the tier in the key and quad arrays.
    pub first_pattern: u64,
    /// Patterns in the tier.
    pub n_patterns: u64,
    /// Anchors (groups) that produced patterns.
    pub n_anchors: u64,
}

impl TierInfo {
    /// The tier's pattern range.
    #[must_use]
    pub fn patterns(&self) -> Range<usize> {
        self.first_pattern as usize..(self.first_pattern + self.n_patterns) as usize
    }
}

/// A star of the index.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IndexStar {
    /// Right ascension, radians.
    pub ra: f32,
    /// Declination, radians.
    pub dec: f32,
    /// Magnitude × 100.
    pub mag: i16,
    /// The widest tier any pattern using this star belongs to.
    pub tier: u8,
}

/// Which copy of a star database an index was built from, so a later check can tell
/// whether the database has changed since (a new ASTAP release, a re-download, a
/// different directory's files).
///
/// The fingerprint covers each database file's name, size and first
/// [`SourceStamp::HEAD_BYTES`] bytes, not its modification time: copying a database
/// to another disk changes every mtime without changing a star. All zero when not
/// recorded (indexes built before the stamp existed).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SourceStamp {
    /// Number of database files.
    pub files: u32,
    /// Their total size in bytes.
    pub bytes: u64,
    /// FNV-1a over every file's name, size and head, in name order.
    pub hash: u64,
}

impl SourceStamp {
    /// Bytes of each file's head folded into the hash.
    pub const HEAD_BYTES: usize = 4096;

    /// Whether the index recorded its source at all.
    #[must_use]
    pub fn is_recorded(&self) -> bool {
        self.files > 0
    }

    /// Fingerprint the files of database `db_name` in `db_path`: every
    /// `<db_name>_*.{1476,290,001}`. Reads a few kilobytes per file.
    ///
    /// # Errors
    ///
    /// [`ArcsecError::CatalogIo`] if the directory or a file cannot be read.
    pub fn of_database(db_path: &Path, db_name: &str) -> Result<Self> {
        use std::io::Read as _;
        let prefix = format!("{db_name}_");
        let mut names: Vec<String> = std::fs::read_dir(db_path)
            .map_err(ArcsecError::CatalogIo)?
            .filter_map(core::result::Result::ok)
            .filter_map(|e| e.file_name().into_string().ok())
            .filter(|n| {
                n.starts_with(&prefix) && [".1476", ".290", ".001"].iter().any(|x| n.ends_with(x))
            })
            .collect();
        names.sort();
        let mut st = Self::default();
        let mut h = Fnv::new();
        let mut buf = vec![0u8; Self::HEAD_BYTES];
        for n in &names {
            let mut f = File::open(db_path.join(n)).map_err(ArcsecError::CatalogIo)?;
            let len = f.metadata().map_err(ArcsecError::CatalogIo)?.len();
            let mut got = 0;
            while got < buf.len() {
                match f.read(&mut buf[got..]) {
                    Ok(0) => break,
                    Ok(k) => got += k,
                    Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                    Err(e) => return Err(ArcsecError::CatalogIo(e)),
                }
            }
            h.update(n.as_bytes());
            h.update(&[0]);
            h.update(&len.to_le_bytes());
            h.update(&buf[..got]);
            st.files += 1;
            st.bytes += len;
        }
        st.hash = h.0;
        Ok(st)
    }
}

/// 64-bit FNV-1a: small, dependency-free, and stable across platforms and releases,
/// which is all a fingerprint stored in a file needs.
struct Fnv(u64);

impl Fnv {
    fn new() -> Self {
        Self(0xCBF2_9CE4_8422_2325)
    }
    fn update(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 ^= u64::from(b);
            self.0 = self.0.wrapping_mul(0x0100_0000_01B3);
        }
    }
}

/// An index assembled in memory, ready to write (the builder's output).
#[derive(Debug, Default)]
pub struct BuiltIndex {
    /// Tiers, widest first.
    pub tiers: Vec<TierInfo>,
    /// Stars, sorted by declination band then RA (see [`star_band`]).
    pub stars: Vec<IndexStar>,
    /// First star of each declination band, `star_bands + 1` entries.
    pub star_dir: Vec<u32>,
    /// Pattern keys, sorted within each tier.
    pub keys: Vec<u64>,
    /// Pattern stars, parallel to `keys`.
    pub quads: Vec<[u32; 4]>,
    /// Source database name.
    pub source: String,
    /// Fingerprint of the source database's files.
    pub source_stamp: SourceStamp,
}

/// Declination bands in the star directory: quarter-degree strips.
pub const STAR_BANDS: u32 = 720;

/// The star-directory band holding declination `dec` (radians).
#[must_use]
pub fn star_band(dec: f64) -> u32 {
    let b = ((dec + core::f64::consts::FRAC_PI_2) / core::f64::consts::PI * f64::from(STAR_BANDS))
        .floor();
    (b.max(0.0) as u32).min(STAR_BANDS - 1)
}

// ── CRC-32 (IEEE 802.3, as zip and PNG use) ────────────────────────────────────

const fn crc_table() -> [u32; 256] {
    let mut t = [0u32; 256];
    let mut i = 0;
    while i < 256 {
        let mut c = i as u32;
        let mut k = 0;
        while k < 8 {
            c = if c & 1 != 0 {
                0xEDB8_8320 ^ (c >> 1)
            } else {
                c >> 1
            };
            k += 1;
        }
        t[i] = c;
        i += 1;
    }
    t
}

static CRC_TABLE: [u32; 256] = crc_table();

/// Incremental CRC-32.
#[derive(Clone, Copy)]
struct Crc(u32);

impl Crc {
    fn new() -> Self {
        Self(0xFFFF_FFFF)
    }
    fn update(&mut self, bytes: &[u8]) {
        let mut c = self.0;
        for &b in bytes {
            c = CRC_TABLE[((c ^ u32::from(b)) & 0xFF) as usize] ^ (c >> 8);
        }
        self.0 = c;
    }
    fn finish(self) -> u32 {
        !self.0
    }
}

fn crc32(bytes: &[u8]) -> u32 {
    let mut c = Crc::new();
    c.update(bytes);
    c.finish()
}

// ── Writing ────────────────────────────────────────────────────────────────────

/// A writer that tracks its offset and the CRC of the current section.
struct SectionWriter<W: Write> {
    w: W,
    pos: u64,
    crc: Crc,
}

impl<W: Write> SectionWriter<W> {
    fn put(&mut self, b: &[u8]) -> io::Result<()> {
        self.w.write_all(b)?;
        self.crc.update(b);
        self.pos += b.len() as u64;
        Ok(())
    }
    /// Start a section: pad to 8 bytes (padding is outside every section's CRC).
    fn begin(&mut self) -> io::Result<u64> {
        let pad = (8 - self.pos % 8) % 8;
        self.w.write_all(&[0u8; 8][..pad as usize])?;
        self.pos += pad;
        self.crc = Crc::new();
        Ok(self.pos)
    }
}

impl BuiltIndex {
    /// Total bytes the file will occupy (to within section padding).
    #[must_use]
    pub fn file_size(&self) -> u64 {
        (HEADER_LEN
            + self.tiers.len() * TIER_LEN
            + self.stars.len() * STAR_LEN
            + self.star_dir.len() * 4
            + self.keys.len() * 8
            + self.quads.len() * QUAD_LEN
            + 32) as u64
    }

    /// Write the index to `path`, via a temporary file renamed into place, so a
    /// reader never sees a half-written index.
    ///
    /// # Errors
    ///
    /// [`ArcsecError::CatalogIo`] if the file cannot be written.
    pub fn write(&self, path: &Path) -> Result<()> {
        let tmp = path.with_extension(format!("{EXTENSION}.part"));
        self.write_to(&tmp).map_err(ArcsecError::CatalogIo)?;
        std::fs::rename(&tmp, path).map_err(ArcsecError::CatalogIo)
    }

    fn write_to(&self, path: &Path) -> io::Result<()> {
        let f = File::create(path)?;
        let mut sw = SectionWriter {
            w: BufWriter::with_capacity(1 << 20, f),
            pos: 0,
            crc: Crc::new(),
        };
        // Placeholder header; rewritten at the end with the section table.
        sw.put(&[0u8; HEADER_LEN])?;
        let mut sections = [(0u64, 0u64, 0u32); N_SECTIONS];

        let off = sw.begin()?;
        for t in &self.tiers {
            sw.put(&t.radius.to_le_bytes())?;
            sw.put(&t.mag_cap.to_le_bytes())?;
            sw.put(&t.members.to_le_bytes())?;
            sw.put(&t.first_pattern.to_le_bytes())?;
            sw.put(&t.n_patterns.to_le_bytes())?;
            sw.put(&t.n_anchors.to_le_bytes())?;
        }
        sections[0] = (off, sw.pos - off, sw.crc.finish());

        let off = sw.begin()?;
        for s in &self.stars {
            let mut r = [0u8; STAR_LEN];
            r[0..4].copy_from_slice(&s.ra.to_le_bytes());
            r[4..8].copy_from_slice(&s.dec.to_le_bytes());
            r[8..10].copy_from_slice(&s.mag.to_le_bytes());
            r[10] = s.tier;
            sw.put(&r)?;
        }
        sections[1] = (off, sw.pos - off, sw.crc.finish());

        let off = sw.begin()?;
        for &d in &self.star_dir {
            sw.put(&d.to_le_bytes())?;
        }
        sections[2] = (off, sw.pos - off, sw.crc.finish());

        let off = sw.begin()?;
        for &k in &self.keys {
            sw.put(&k.to_le_bytes())?;
        }
        sections[3] = (off, sw.pos - off, sw.crc.finish());

        let off = sw.begin()?;
        for q in &self.quads {
            let mut r = [0u8; QUAD_LEN];
            for (i, &s) in q.iter().enumerate() {
                r[i * 4..i * 4 + 4].copy_from_slice(&s.to_le_bytes());
            }
            sw.put(&r)?;
        }
        sections[4] = (off, sw.pos - off, sw.crc.finish());

        let mut h = [0u8; HEADER_LEN];
        h[0..8].copy_from_slice(MAGIC);
        h[8..12].copy_from_slice(&VERSION.to_le_bytes());
        h[12..16].copy_from_slice(&BYTE_ORDER.to_le_bytes());
        h[16..20].copy_from_slice(&(HEADER_LEN as u32).to_le_bytes());
        h[20..24].copy_from_slice(&(super::pattern::BINS as u32).to_le_bytes());
        h[24..28].copy_from_slice(&(self.tiers.len() as u32).to_le_bytes());
        h[28..32].copy_from_slice(&STAR_BANDS.to_le_bytes());
        h[32..40].copy_from_slice(&(self.stars.len() as u64).to_le_bytes());
        h[40..48].copy_from_slice(&(self.keys.len() as u64).to_le_bytes());
        let now = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_secs() as i64);
        h[48..56].copy_from_slice(&now.to_le_bytes());
        let src = self.source.as_bytes();
        let n = src.len().min(16);
        h[56..56 + n].copy_from_slice(&src[..n]);
        h[STAMP_AT..STAMP_AT + 8].copy_from_slice(&self.source_stamp.hash.to_le_bytes());
        h[STAMP_AT + 8..STAMP_AT + 16].copy_from_slice(&self.source_stamp.bytes.to_le_bytes());
        h[STAMP_AT + 16..STAMP_AT + 20].copy_from_slice(&self.source_stamp.files.to_le_bytes());
        for (i, (o, l, c)) in sections.iter().enumerate() {
            let at = SECTIONS_AT + i * 24;
            h[at..at + 8].copy_from_slice(&o.to_le_bytes());
            h[at + 8..at + 16].copy_from_slice(&l.to_le_bytes());
            h[at + 16..at + 20].copy_from_slice(&c.to_le_bytes());
        }
        let crc = crc32(&h[..HEADER_CRC_AT]);
        h[HEADER_CRC_AT..].copy_from_slice(&crc.to_le_bytes());

        let mut f = sw.w.into_inner().map_err(io::IntoInnerError::into_error)?;

        f.seek(io::SeekFrom::Start(0))?;
        f.write_all(&h)?;
        f.sync_all()
    }
}

// ── Reading ────────────────────────────────────────────────────────────────────

#[inline]
fn u32_at(b: &[u8], at: usize) -> u32 {
    u32::from_le_bytes([b[at], b[at + 1], b[at + 2], b[at + 3]])
}

#[inline]
fn u64_at(b: &[u8], at: usize) -> u64 {
    let mut x = [0u8; 8];
    x.copy_from_slice(&b[at..at + 8]);
    u64::from_le_bytes(x)
}

#[inline]
fn f32_at(b: &[u8], at: usize) -> f32 {
    f32::from_bits(u32_at(b, at))
}

fn invalid(path: &Path, why: impl core::fmt::Display) -> ArcsecError {
    ArcsecError::CatalogIo(io::Error::new(
        io::ErrorKind::InvalidData,
        format!("{}: not a usable arcsec blind index: {why}", path.display()),
    ))
}

/// A memory-mapped blind index.
pub struct BlindIndex {
    map: Mmap,
    tiers: Vec<TierInfo>,
    sections: [(usize, usize, u32); N_SECTIONS],
    n_stars: usize,
    n_patterns: usize,
    star_bands: u32,
    source: String,
    source_stamp: SourceStamp,
    built_unix: i64,
}

impl core::fmt::Debug for BlindIndex {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("BlindIndex")
            .field("source", &self.source)
            .field("n_stars", &self.n_stars)
            .field("n_patterns", &self.n_patterns)
            .field("tiers", &self.tiers)
            .finish_non_exhaustive()
    }
}

/// Whether `path` starts with the index magic. Reads eight bytes.
#[must_use]
pub fn is_blind_index(path: &Path) -> bool {
    use std::io::Read as _;
    let mut m = [0u8; 8];
    File::open(path)
        .and_then(|mut f| f.read_exact(&mut m))
        .is_ok()
        && &m == MAGIC
}

impl BlindIndex {
    /// Map and check an index file. Cheap whatever the file's size: see the module
    /// documentation for what is and is not checked here.
    ///
    /// # Errors
    ///
    /// [`ArcsecError::CatalogIo`] if the file cannot be opened, is not an index, is
    /// another version or byte order, or its header is inconsistent.
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path).map_err(ArcsecError::CatalogIo)?;
        // Safety: read-only mapping. As with the star databases, the file is not
        // expected to change while a solve is reading it; the writer replaces it by
        // rename, which leaves an existing mapping intact.
        let map = unsafe { Mmap::map(&file) }.map_err(ArcsecError::CatalogIo)?;
        if map.len() < HEADER_LEN || &map[..8] != MAGIC {
            return Err(invalid(path, "bad magic"));
        }
        let version = u32_at(&map, 8);
        if version != VERSION {
            return Err(invalid(
                path,
                format!("format version {version}, this build reads {VERSION} - rebuild it"),
            ));
        }
        if u32_at(&map, 12) != BYTE_ORDER {
            return Err(invalid(path, "written with the other byte order"));
        }
        if crc32(&map[..HEADER_CRC_AT]) != u32_at(&map, HEADER_CRC_AT) {
            return Err(invalid(path, "header checksum mismatch"));
        }
        if u32_at(&map, 16) as usize != HEADER_LEN
            || f64::from(u32_at(&map, 20)) != super::pattern::BINS
        {
            return Err(invalid(
                path,
                "unsupported header length or descriptor bins",
            ));
        }
        let n_tiers = u32_at(&map, 24) as usize;
        let star_bands = u32_at(&map, 28);
        let n_stars = usize::try_from(u64_at(&map, 32)).map_err(|_| invalid(path, "too large"))?;
        let n_patterns =
            usize::try_from(u64_at(&map, 40)).map_err(|_| invalid(path, "too large"))?;
        let built_unix = u64_at(&map, 48) as i64;
        let source = String::from_utf8_lossy(&map[56..72])
            .trim_end_matches('\0')
            .to_string();
        let source_stamp = SourceStamp {
            hash: u64_at(&map, STAMP_AT),
            bytes: u64_at(&map, STAMP_AT + 8),
            files: u32_at(&map, STAMP_AT + 16),
        };

        let mut sections = [(0usize, 0usize, 0u32); N_SECTIONS];
        let expected = [
            n_tiers.checked_mul(TIER_LEN),
            n_stars.checked_mul(STAR_LEN),
            (star_bands as usize + 1).checked_mul(4),
            n_patterns.checked_mul(8),
            n_patterns.checked_mul(QUAD_LEN),
        ];
        for (i, s) in sections.iter_mut().enumerate() {
            let at = SECTIONS_AT + i * 24;
            let off = usize::try_from(u64_at(&map, at)).map_err(|_| invalid(path, "too large"))?;
            let len =
                usize::try_from(u64_at(&map, at + 8)).map_err(|_| invalid(path, "too large"))?;
            if Some(len) != expected[i] || off.checked_add(len).is_none_or(|e| e > map.len()) {
                return Err(invalid(
                    path,
                    format!("section {i} is truncated or mis-sized"),
                ));
            }
            *s = (off, len, u32_at(&map, at + 16));
        }

        let tb = &map[sections[0].0..sections[0].0 + sections[0].1];
        let mut tiers = Vec::with_capacity(n_tiers);
        for t in 0..n_tiers {
            let r = &tb[t * TIER_LEN..(t + 1) * TIER_LEN];
            let tier = TierInfo {
                radius: f64::from_bits(u64_at(r, 0)),
                mag_cap: f32_at(r, 8),
                members: u32_at(r, 12),
                first_pattern: u64_at(r, 16),
                n_patterns: u64_at(r, 24),
                n_anchors: u64_at(r, 32),
            };
            let end = tier.first_pattern.checked_add(tier.n_patterns);
            if !(tier.radius.is_finite() && tier.radius > 0.0)
                || end.is_none_or(|e| e > n_patterns as u64)
            {
                return Err(invalid(path, format!("tier {t} is inconsistent")));
            }
            tiers.push(tier);
        }

        Ok(Self {
            map,
            tiers,
            sections,
            n_stars,
            n_patterns,
            star_bands,
            source,
            source_stamp,
            built_unix,
        })
    }

    /// Recompute every section checksum: a full read of the file.
    ///
    /// # Errors
    ///
    /// [`ArcsecError::CatalogIo`] naming the first section whose checksum differs.
    pub fn validate(&self) -> Result<()> {
        const NAMES: [&str; N_SECTIONS] = ["tiers", "stars", "star directory", "keys", "quads"];
        for (i, &(off, len, crc)) in self.sections.iter().enumerate() {
            if crc32(&self.map[off..off + len]) != crc {
                return Err(ArcsecError::CatalogIo(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("blind index: {} section checksum mismatch", NAMES[i]),
                )));
            }
        }
        Ok(())
    }

    /// Tiers, widest first.
    #[must_use]
    pub fn tiers(&self) -> &[TierInfo] {
        &self.tiers
    }

    /// Number of stars.
    #[must_use]
    pub fn n_stars(&self) -> usize {
        self.n_stars
    }

    /// Number of patterns, all tiers.
    #[must_use]
    pub fn n_patterns(&self) -> usize {
        self.n_patterns
    }

    /// Database the index was built from.
    #[must_use]
    pub fn source(&self) -> &str {
        &self.source
    }

    /// Fingerprint of the database the index was built from; not recorded
    /// ([`SourceStamp::is_recorded`] false) in indexes built before it existed.
    #[must_use]
    pub fn source_stamp(&self) -> SourceStamp {
        self.source_stamp
    }

    /// Build time, seconds since the Unix epoch.
    #[must_use]
    pub fn built_unix(&self) -> i64 {
        self.built_unix
    }

    /// File size in bytes.
    #[must_use]
    pub fn file_size(&self) -> usize {
        self.map.len()
    }

    /// Bytes of the file belonging to one tier's patterns (keys and quads).
    #[must_use]
    pub fn tier_bytes(&self, t: &TierInfo) -> u64 {
        t.n_patterns * (8 + QUAD_LEN as u64)
    }

    #[inline]
    fn key(&self, i: usize) -> u64 {
        u64_at(&self.map, self.sections[3].0 + i * 8)
    }

    /// The patterns of `tier` whose key is exactly `key`.
    #[must_use]
    pub fn lookup(&self, tier: &TierInfo, key: u64) -> Range<usize> {
        let r = tier.patterns();
        let (mut lo, mut hi) = (r.start, r.end);
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if self.key(mid) < key {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        let start = lo;
        let mut hi = r.end;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if self.key(mid) <= key {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        start..lo
    }

    /// The four stars of pattern `i` in canonical order; `None` if the pattern
    /// index or any star index is out of range (a corrupt file).
    #[must_use]
    pub fn quad(&self, i: usize) -> Option<[IndexStar; 4]> {
        if i >= self.n_patterns {
            return None;
        }
        let at = self.sections[4].0 + i * QUAD_LEN;
        let mut out = [IndexStar {
            ra: 0.0,
            dec: 0.0,
            mag: 0,
            tier: 0,
        }; 4];
        for (k, s) in out.iter_mut().enumerate() {
            *s = self.star(u32_at(&self.map, at + k * 4) as usize)?;
        }
        Some(out)
    }

    /// Star `i`; `None` if out of range.
    #[must_use]
    pub fn star(&self, i: usize) -> Option<IndexStar> {
        if i >= self.n_stars {
            return None;
        }
        let at = self.sections[1].0 + i * STAR_LEN;
        let b = &self.map[at..at + STAR_LEN];
        Some(IndexStar {
            ra: f32_at(b, 0),
            dec: f32_at(b, 4),
            mag: i16::from_le_bytes([b[8], b[9]]),
            tier: b[10],
        })
    }

    /// Call `f` for every star within `radius` (radians) of (`ra`, `dec`), and a few
    /// just outside it: candidates come from whole directory bands and a RA window
    /// widened for the band's declination, and the caller does the exact test.
    pub fn stars_near(&self, ra: f64, dec: f64, radius: f64, mut f: impl FnMut(&IndexStar)) {
        use core::f64::consts::{FRAC_PI_2, PI};
        if self.star_bands != STAR_BANDS {
            return;
        }
        let lo_b = star_band((dec - radius).max(-FRAC_PI_2));
        let hi_b = star_band((dec + radius).min(FRAC_PI_2));
        let dir_at = self.sections[2].0;
        for band in lo_b..=hi_b {
            let s0 = u32_at(&self.map, dir_at + band as usize * 4) as usize;
            let s1 =
                (u32_at(&self.map, dir_at + (band as usize + 1) * 4) as usize).min(self.n_stars);
            if s0 >= s1 {
                continue;
            }
            // The band's worst-case cos(dec), for the RA half-width.
            let b_lo = f64::from(band) / f64::from(STAR_BANDS) * PI - FRAC_PI_2;
            let b_hi = b_lo + PI / f64::from(STAR_BANDS);
            let cos_min = b_lo.cos().min(b_hi.cos()).max(0.0);
            let half = if cos_min * PI <= radius || dec.abs() + radius >= FRAC_PI_2 {
                PI
            } else {
                (radius / cos_min).min(PI)
            };
            if half >= PI {
                for i in s0..s1 {
                    if let Some(s) = self.star(i) {
                        f(&s);
                    }
                }
                continue;
            }
            let ra0 = (ra - half).rem_euclid(2.0 * PI);
            let ra1 = (ra + half).rem_euclid(2.0 * PI);
            let windows: &[(f64, f64)] = if ra0 <= ra1 {
                &[(ra0, ra1)]
            } else {
                &[(ra0, 2.0 * PI), (0.0, ra1)]
            };
            for &(w0, w1) in windows {
                // Binary search for the first star with RA >= w0.
                let (mut lo, mut hi) = (s0, s1);
                while lo < hi {
                    let mid = lo + (hi - lo) / 2;
                    let r = self.star(mid).map_or(f32::INFINITY, |s| s.ra);
                    if f64::from(r) < w0 {
                        lo = mid + 1;
                    } else {
                        hi = mid;
                    }
                }
                for i in lo..s1 {
                    let Some(s) = self.star(i) else { break };
                    if f64::from(s.ra) > w1 {
                        break;
                    }
                    f(&s);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::TempDir;

    fn sample() -> BuiltIndex {
        let stars: Vec<IndexStar> = (0..10)
            .map(|i| IndexStar {
                ra: 0.1 * i as f32,
                dec: 0.2,
                mag: 1000 + i,
                tier: 0,
            })
            .collect();
        let mut dir = vec![0u32; STAR_BANDS as usize + 1];
        let b = star_band(0.2) as usize;
        for (i, d) in dir.iter_mut().enumerate() {
            *d = if i <= b { 0 } else { 10 };
        }
        BuiltIndex {
            tiers: vec![TierInfo {
                radius: 0.01,
                mag_cap: 12.0,
                members: 5,
                first_pattern: 0,
                n_patterns: 3,
                n_anchors: 1,
            }],
            stars,
            star_dir: dir,
            keys: vec![5, 7, 7],
            quads: vec![[0, 1, 2, 3], [1, 2, 3, 4], [5, 6, 7, 8]],
            source: "d80".into(),
            source_stamp: SourceStamp::default(),
        }
    }

    #[test]
    fn round_trips_and_looks_up() {
        let dir = TempDir::new("arcsecix_rt");
        let p = dir.path().join("t.arcsecix");
        sample().write(&p).unwrap();
        assert!(is_blind_index(&p));
        let ix = BlindIndex::open(&p).unwrap();
        ix.validate().unwrap();
        assert_eq!(ix.source(), "d80");
        assert_eq!(ix.n_stars(), 10);
        let t = ix.tiers()[0];
        assert_eq!(ix.lookup(&t, 7), 1..3);
        assert_eq!(ix.lookup(&t, 5), 0..1);
        assert!(ix.lookup(&t, 6).is_empty());
        assert_eq!(ix.quad(2).unwrap()[3].mag, 1008);
        let mut n = 0;
        ix.stars_near(0.3, 0.2, 0.11, |_| n += 1);
        assert!((3..=5).contains(&n), "{n}");
    }

    #[test]
    fn rejects_bad_magic_truncation_and_corruption() {
        let dir = TempDir::new("arcsecix_bad");
        let p = dir.path().join("t.arcsecix");
        sample().write(&p).unwrap();
        let good = std::fs::read(&p).unwrap();

        let mut b = good.clone();
        b[0] = b'X';
        std::fs::write(&p, &b).unwrap();
        assert!(BlindIndex::open(&p).is_err());

        std::fs::write(&p, &good[..good.len() - 20]).unwrap();
        assert!(BlindIndex::open(&p).is_err(), "truncated");

        let mut b = good.clone();
        b[30] ^= 1; // inside the header: header CRC
        std::fs::write(&p, &b).unwrap();
        assert!(BlindIndex::open(&p).is_err());

        let mut b = good.clone();
        let n = b.len();
        b[n - 3] ^= 0x40; // inside the quads: opens, fails validation
        std::fs::write(&p, &b).unwrap();
        let ix = BlindIndex::open(&p).unwrap();
        assert!(ix.validate().is_err());
    }

    #[test]
    fn a_corrupt_star_reference_is_refused_not_followed() {
        let dir = TempDir::new("arcsecix_ref");
        let p = dir.path().join("t.arcsecix");
        let mut s = sample();
        s.quads[0] = [0, 1, 2, 999];
        s.write(&p).unwrap();
        let ix = BlindIndex::open(&p).unwrap();
        assert!(ix.quad(0).is_none());
        assert!(ix.quad(99).is_none());
    }

    #[test]
    fn the_source_stamp_round_trips_and_an_unstamped_file_reads_as_unrecorded() {
        let dir = TempDir::new("arcsecix_stamp");
        let p = dir.path().join("t.arcsecix");
        sample().write(&p).unwrap();
        let ix = BlindIndex::open(&p).unwrap();
        assert!(!ix.source_stamp().is_recorded(), "zeros mean not recorded");

        let mut s = sample();
        s.source_stamp = SourceStamp {
            files: 1476,
            bytes: 1_300_000_000,
            hash: 0x0123_4567_89AB_CDEF,
        };
        s.write(&p).unwrap();
        let ix = BlindIndex::open(&p).unwrap();
        ix.validate().unwrap();
        assert_eq!(ix.source_stamp(), s.source_stamp);
    }

    #[test]
    fn the_database_stamp_sees_size_and_content_but_not_other_files() {
        let dir = TempDir::new("arcsecix_dbstamp");
        let d = dir.path();
        std::fs::write(d.join("t_0101.1476"), vec![1u8; 5000]).unwrap();
        std::fs::write(d.join("t_0201.1476"), vec![2u8; 300]).unwrap();
        std::fs::write(d.join("u_0101.1476"), b"another database").unwrap();
        std::fs::write(d.join("t.arcsecix"), b"not a database file").unwrap();
        let a = SourceStamp::of_database(d, "t").unwrap();
        assert_eq!((a.files, a.bytes), (2, 5300));
        assert!(a.is_recorded());
        assert_eq!(
            a,
            SourceStamp::of_database(d, "t").unwrap(),
            "deterministic"
        );

        // Unrelated files do not count.
        std::fs::write(d.join("u_0201.1476"), b"more").unwrap();
        assert_eq!(a, SourceStamp::of_database(d, "t").unwrap());

        // Same size, different head: a different database.
        let mut b = vec![1u8; 5000];
        b[10] = 9;
        std::fs::write(d.join("t_0101.1476"), &b).unwrap();
        let c = SourceStamp::of_database(d, "t").unwrap();
        assert_eq!(c.bytes, a.bytes);
        assert_ne!(c.hash, a.hash);

        // A missing database: an empty, unrecorded stamp.
        assert!(!SourceStamp::of_database(d, "zz").unwrap().is_recorded());
    }

    #[test]
    fn crc_matches_the_standard_check_value() {
        assert_eq!(crc32(b"123456789"), 0xCBF4_3926);
    }
}
