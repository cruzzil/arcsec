//! Astrometry.net index file reader.
//!
//! Reads the raw `quads`, `kdtree_data_stars`, and `kdtree_data_codes` binary-table
//! extensions from a FITS index file (e.g. `index-4112.fits`, DIMQUADS=3).
//!
//! All three use TFORM='nA' (raw bytes). astrometry.net writes them in the builder's
//! native byte order; every distributed index is little-endian, and that is what
//! this reader assumes:
//!
//! ```text
//! Stars:  3 × u32 LE (x,y,z) on unit sphere, fixed-point mapped [-1,+1]→[0,2³²-1]
//! Quads:  dim_quads × u32 LE — star indices A, B, C[, D] into the star table
//! Codes:  n_code_dims × u16 LE — decoded via range table
//! ```
//!
//! In astrometry.net's convention A = `star[0]` and B = `star[1]` are the two
//! most-separated stars, and define the AB frame.
//!
//! DIMQUADS=3 (triangles): 2 code dims (CX, CY) — position of C in the AB frame.
//! DIMQUADS=4 (quads):     4 code dims (CX, CY, DX, DY) — C and D in the AB frame.
//!
//! Canonical form (swap A↔B and invert all codes if violated): CX ≤ 0.5 for
//! triangles, CX + DX ≤ 1 for quads. For DIMQUADS=4 additionally CX ≤ DX (swap C↔D
//! if violated). See `pipeline::blind::make_quad4` for why the two rules differ.

use core::f64::consts::PI;
use std::io;
use std::path::Path;

use libc::{c_char, c_int, c_long, c_void};
use rsfitsio::aliases::rust_api::{
    fits_close_file, fits_get_colnum, fits_get_num_rowsll, fits_movabs_hdu, fits_movnam_hdu,
    fits_movrel_hdu, fits_open_image, fits_open_memfile, fits_read_col_byt, fits_read_key_dbl,
    fits_read_key_lng, fits_read_key_str,
};
use rsfitsio::fitsio::{ANY_HDU, FLEN_VALUE, LONGLONG, READONLY, fitsfile};

use crate::error::ArcsecError;
use crate::math::coords::ang_sep;

// ── Public types ───────────────────────────────────────────────────────────────

/// A star in the index, in sky coordinates (radians).
#[derive(Debug, Clone, Copy)]
pub struct AnetStar {
    /// Right ascension, radians in `[0, 2π)`.
    pub ra: f64,
    /// Declination, radians.
    pub dec: f64,
}

/// An index entry (triangle or quad) from the file.
///
/// `n_stars` is 3 for DIMQUADS=3 entries, 4 for DIMQUADS=4.
/// Only the first `n_stars` elements of `star_ra/dec` are valid.
/// `code` holds `n_code_dims = 2*(dim_quads-2)` code values in `[0]`, `[1]`, ...;
/// unused slots are 0.0.
#[derive(Debug, Clone)]
pub struct AnetIndexEntry {
    /// Code values: [CX, CY] for triangles; [CX, CY, DX, DY] for quads.
    pub code: [f64; 4],
    /// Number of valid stars (3 or 4).
    pub n_stars: usize,
    /// RA of constituent stars in file order.
    pub star_ra: [f64; 4],
    /// Dec of constituent stars in file order.
    pub star_dec: [f64; 4],
    /// Approximate centroid RA (radians).
    pub center_ra: f64,
    /// Approximate centroid Dec (radians).
    pub center_dec: f64,
}

/// In-memory Astrometry.net index, ready for code lookups.
pub struct AnetIndex {
    /// Entries sorted by `code[0]` ascending for binary-search lookup.
    pub entries: Vec<AnetIndexEntry>,
    /// Compact parallel code array (f32, 16 bytes/entry) for cache-efficient search.
    /// Sorted identically to `entries`; `find_code_matches` scans this 9 MB array
    /// (fits in L3) instead of the 70 MB `entries` array (all DRAM cache misses).
    pub codes: Vec<[f32; 4]>,
    /// All catalog stars in the index (for WCS verification).
    pub stars: Vec<AnetStar>,
    /// Smallest quad scale (A-B separation) in this index file, radians.
    pub scale_lo: f64,
    /// Largest quad scale (A-B separation) in this index file, radians.
    pub scale_hi: f64,
    /// Stars per quad entry (DIMQUADS in the header; 3 or 4).
    pub dim_quads: usize,
}

// ── FITS helpers ───────────────────────────────────────────────────────────────

/// Reinterpret a byte slice as a `c_char` slice for the rsfitsio wrappers.
///
/// `libc::c_char` is `i8` on `x86_64` and `u8` on aarch64, so `b"KEY\0"` literals
/// cannot be passed directly on every platform. The two types always have the same
/// size and alignment, so the cast is a no-op at runtime.
#[inline]
fn cc(b: &[u8]) -> &[c_char] {
    // Safety: c_char is i8 or u8; identical layout, and we only read.
    unsafe { core::slice::from_raw_parts(b.as_ptr().cast::<c_char>(), b.len()) }
}

/// Decode a NUL-terminated `c_char` buffer filled in by CFITSIO into a String.
fn cstr_to_string(buf: &[c_char]) -> String {
    let end = buf.iter().position(|&c| c == 0).unwrap_or(buf.len());
    let bytes: Vec<u8> = buf[..end].iter().map(|&c| c as u8).collect();
    String::from_utf8_lossy(&bytes).trim().to_string()
}

/// Move to a named HDU by matching TTYPE1 or EXTNAME keyword.
fn move_to_hdu(fp: &mut fitsfile, target: &[u8]) -> Result<(), String> {
    {
        let mut status: c_int = 0;
        fits_movnam_hdu(fp, ANY_HDU, cc(target), 0, &mut status);
        if status == 0 {
            return Ok(());
        }
    }

    let target_str = String::from_utf8_lossy(target)
        .trim_end_matches('\0')
        .to_lowercase();

    {
        let mut status: c_int = 0;
        fits_movabs_hdu(fp, 1, None, &mut status);
    }

    loop {
        let mut status: c_int = 0;
        fits_movrel_hdu(fp, 1, None, &mut status);
        if status != 0 {
            return Err(format!("HDU with TTYPE1='{target_str}' not found"));
        }

        let mut val = vec![0 as c_char; FLEN_VALUE];
        let mut st: c_int = 0;
        fits_read_key_str(fp, cc(b"TTYPE1\0"), &mut val, None, &mut st);
        if st == 0 {
            let name = cstr_to_string(&val).to_lowercase();
            if name == target_str {
                return Ok(());
            }
        }
    }
}

fn get_num_rows(fp: &mut fitsfile) -> Result<usize, String> {
    let mut nrows: LONGLONG = 0;
    let mut st: c_int = 0;
    fits_get_num_rowsll(fp, &mut nrows, &mut st);
    if st == 0 {
        Ok(nrows as usize)
    } else {
        Err(format!("fits_get_num_rows: {st}"))
    }
}

fn get_colnum(fp: &mut fitsfile, name: &[u8]) -> Result<c_int, String> {
    let mut col: c_int = 0;
    let mut st: c_int = 0;
    fits_get_colnum(fp, 0, cc(name), &mut col, &mut st);
    if st == 0 {
        Ok(col)
    } else {
        Err(format!("fits_get_colnum: {st}"))
    }
}

fn read_raw_bytes(fp: &mut fitsfile, col: c_int, n: usize) -> Result<Vec<u8>, String> {
    let mut bytes = vec![0u8; n];
    let mut st: c_int = 0;
    // nelem is LONGLONG, not c_long — the two differ on Windows, where c_long is i32.
    fits_read_col_byt(fp, col, 1, 1, n as LONGLONG, 0, &mut bytes, None, &mut st);
    if st == 0 {
        Ok(bytes)
    } else {
        Err(format!("fits_read_col_byt: {st}"))
    }
}

// ── Star parsing ───────────────────────────────────────────────────────────────

/// Convert a u32 in [0, 2³²-1] to a f64 in [-1, +1].
fn u32_to_unit(v: u32) -> f64 {
    v as f64 / (u32::MAX as f64 / 2.0) - 1.0
}

fn parse_stars(bytes: &[u8]) -> Vec<AnetStar> {
    bytes
        .as_chunks::<12>()
        .0
        .iter()
        .map(|chunk| {
            let xf = u32_to_unit(u32::from_le_bytes(chunk[0..4].try_into().unwrap()));
            let yf = u32_to_unit(u32::from_le_bytes(chunk[4..8].try_into().unwrap()));
            let zf = u32_to_unit(u32::from_le_bytes(chunk[8..12].try_into().unwrap()));
            let ra = yf.atan2(xf).rem_euclid(2.0 * PI);
            let dec = zf.atan2((xf * xf + yf * yf).sqrt());
            AnetStar { ra, dec }
        })
        .collect()
}

// ── Code and quad-index parsing ────────────────────────────────────────────────

/// Parse `bytes` as little-endian u32s. Rows of `dim_quads` are contiguous, so the
/// caller walks the result with `chunks_exact(dim_quads)`; one flat vector avoids an
/// allocation per quad (millions of them in a large index).
fn parse_quad_indices(bytes: &[u8]) -> Vec<u32> {
    bytes
        .as_chunks::<4>()
        .0
        .iter()
        .map(|&b| u32::from_le_bytes(b))
        .collect()
}

/// Parse code bytes: `n_code_dims` × u16 LE per row, decoded with `lo + u16 / scale`.
/// Returns codes as `[f64; 4]` with unused slots set to 0.0.
fn parse_codes(bytes: &[u8], n_code_dims: usize, code_lo: f64, code_scale: f64) -> Vec<[f64; 4]> {
    let row_bytes = n_code_dims * 2;
    bytes
        .chunks_exact(row_bytes)
        .map(|chunk| {
            let mut code = [0.0f64; 4];
            for (i, slot) in code.iter_mut().take(n_code_dims).enumerate() {
                let off = i * 2;
                let u = u16::from_le_bytes(chunk[off..off + 2].try_into().unwrap());
                *slot = code_lo + u as f64 / code_scale;
            }
            code
        })
        .collect()
}

// ── Entry builder ──────────────────────────────────────────────────────────────

/// Build an `AnetIndexEntry` from star indices and pre-computed code values.
fn entry_from_sky(
    stars: &[AnetStar],
    indices: &[u32],
    code: [f64; 4],
    n_stars: usize,
) -> Option<AnetIndexEntry> {
    if !(3..=4).contains(&n_stars) {
        return None;
    }

    let mut star_ra = [0.0f64; 4];
    let mut star_dec = [0.0f64; 4];
    for i in 0..n_stars {
        let s = stars.get(indices[i] as usize)?;
        star_ra[i] = s.ra;
        star_dec[i] = s.dec;
    }

    let a = stars.get(indices[0] as usize)?;
    let b = stars.get(indices[1] as usize)?;
    if ang_sep(a.ra, a.dec, b.ra, b.dec) < 1e-20 {
        return None;
    }

    // Centroid: vector average of all n_stars on the unit sphere.
    let (mut sx, mut sy, mut sz) = (0.0, 0.0, 0.0);
    for i in 0..n_stars {
        let cos_d = star_dec[i].cos();
        sx += cos_d * star_ra[i].cos();
        sy += cos_d * star_ra[i].sin();
        sz += star_dec[i].sin();
    }
    let center_ra = sy.atan2(sx).rem_euclid(2.0 * PI);
    let center_dec = sz.atan2((sx * sx + sy * sy).sqrt());

    Some(AnetIndexEntry {
        code,
        n_stars,
        star_ra,
        star_dec,
        center_ra,
        center_dec,
    })
}

// ── Public API ─────────────────────────────────────────────────────────────────

impl AnetIndex {
    /// Number of code dimensions: 2 for DIMQUADS=3, 4 for DIMQUADS=4.
    #[must_use]
    pub fn n_code_dims(&self) -> usize {
        2 * self.dim_quads.saturating_sub(2)
    }

    /// Find index entries whose code is within `tol` (Euclidean in `n_code_dims` space)
    /// of `code[0..n_code_dims()]`.
    ///
    /// Binary-searches the compact f32 `codes` array (9 MB, L3-resident) rather
    /// than the full `entries` array (70 MB, DRAM). Reduces cache-miss rate ~7.5×
    /// and avoids loading star RA/Dec data for non-matching entries.
    #[must_use]
    pub fn find_code_matches(&self, code: &[f64; 4], tol: f64) -> Vec<usize> {
        let mut out = Vec::new();
        self.find_code_matches_into(code, tol, &mut out);
        out
    }

    /// Like `find_code_matches` but reuses a caller-provided buffer (cleared first).
    /// Avoids repeated malloc/free in hot loops; uses the compact codes[] array for L3 locality.
    pub fn find_code_matches_into(&self, code: &[f64; 4], tol: f64, out: &mut Vec<usize>) {
        out.clear();
        let n = self.n_code_dims().min(4);
        let lo = (code[0] - tol) as f32;
        let hi = (code[0] + tol) as f32;
        let start = self.codes.partition_point(|c| c[0] < lo);
        let end = self.codes.partition_point(|c| c[0] <= hi);
        let tol_sq = (tol * tol) as f32;
        let code_f32 = [
            code[0] as f32,
            code[1] as f32,
            code[2] as f32,
            code[3] as f32,
        ];
        for (i, c) in (start..end).zip(&self.codes[start..end]) {
            let d1 = c[1] - code_f32[1];
            if d1 * d1 > tol_sq {
                continue;
            }
            let d0 = c[0] - code_f32[0];
            let mut dist_sq = d0 * d0 + d1 * d1;
            if n > 2 {
                let d2 = c[2] - code_f32[2];
                if d2 * d2 > tol_sq {
                    continue;
                }
                dist_sq += d2 * d2;
                if n > 3 {
                    let d3 = c[3] - code_f32[3];
                    if d3 * d3 > tol_sq {
                        continue;
                    }
                    dist_sq += d3 * d3;
                }
            }
            if dist_sq <= tol_sq {
                out.push(i);
            }
        }
    }
}

// ── FITS loader ────────────────────────────────────────────────────────────────

fn read_key_int(fp: &mut fitsfile, key: &[u8]) -> i64 {
    // fits_read_key_lng writes a c_long, which is i32 on Windows and i64 on unix,
    // so the destination must be c_long and widen afterwards rather than be i64.
    let mut val: c_long = 0;
    let mut st: c_int = 0;
    fits_read_key_lng(fp, cc(key), &mut val, None, &mut st);
    if st == 0 { val as i64 } else { 0 }
}

fn read_key_dbl(fp: &mut fitsfile, key: &[u8]) -> Option<f64> {
    let mut val = 0.0f64;
    let mut st: c_int = 0;
    fits_read_key_dbl(fp, cc(key), &mut val, None, &mut st);
    if st == 0 { Some(val) } else { None }
}

/// Read the code scale parameters from `kdtree_range_codes`.
///
/// For `n_code_dims` dimensions the range table has `2*n_code_dims + 1` rows of 8 bytes:
///   rows 0..(n-1)         → lo[i] for each dim  (all equal in astrometry.net)
///   rows n..(2n-1)        → hi[i] for each dim
///   row  2n               → scale  (same for all dims)
/// Returns (lo[0], scale).
fn read_code_range(fp: &mut fitsfile, n_code_dims: usize) -> (f64, f64) {
    let n_rows = 2 * n_code_dims + 1;
    let total_bytes = n_rows * 8;

    let mut st: c_int = 0;
    fits_movabs_hdu(fp, 1, None, &mut st);

    loop {
        let mut status: c_int = 0;
        fits_movrel_hdu(fp, 1, None, &mut status);
        if status != 0 {
            break;
        }

        let mut val = vec![0 as c_char; FLEN_VALUE];
        let mut st2: c_int = 0;
        fits_read_key_str(fp, cc(b"TTYPE1\0"), &mut val, None, &mut st2);
        if st2 == 0 {
            let name = cstr_to_string(&val).to_lowercase();
            if name == "kdtree_range_codes" {
                let mut col: c_int = 0;
                let mut cst: c_int = 0;
                fits_get_colnum(fp, 0, cc(b"kdtree_range_codes\0"), &mut col, &mut cst);
                if cst != 0 {
                    break;
                }

                let mut bytes = vec![0u8; total_bytes];
                let mut rstat: c_int = 0;
                fits_read_col_byt(
                    fp,
                    col,
                    1,
                    1,
                    total_bytes as LONGLONG,
                    0,
                    &mut bytes,
                    None,
                    &mut rstat,
                );
                if rstat != 0 {
                    break;
                }

                let lo = f64::from_le_bytes(bytes[0..8].try_into().unwrap());
                let scale_off = 2 * n_code_dims * 8;
                let scale = f64::from_le_bytes(bytes[scale_off..scale_off + 8].try_into().unwrap());
                return (lo, scale);
            }
        }
    }

    // Fallback: use the values we know from the comment in the file.
    (-0.207_107, 46_340.2)
}

/// What `load_anet_index` pulls out of the file before the CFITSIO handle is closed.
struct AnetParts {
    stars: Vec<AnetStar>,
    /// Flat star indices, `dim_quads` per quad.
    quad_indices: Vec<u32>,
    codes: Vec<[f64; 4]>,
    n_quad_rows: usize,
    scale_lo: f64,
    scale_hi: f64,
    dim_quads: usize,
}

/// Load an Astrometry.net FITS index file and return an `AnetIndex`.
///
/// The entire file is pre-read into memory with a single `read()` syscall, then
/// handed to CFITSIO's in-memory driver (`fits_open_memfile`).  CFITSIO parses
/// the FITS structure as usual but accesses data through a memory pointer with
/// zero `read()` syscalls — eliminating the ~30 s of WSL2 per-syscall overhead
/// seen with file-backed I/O (CFITSIO's default 2880-byte block reads).
///
/// # Errors
///
/// [`ArcsecError::CatalogIo`] if the file cannot be read, is not a FITS file, has a
/// DIMQUADS other than 3 or 4, or lacks the expected binary-table extensions.
pub fn load_anet_index(path: &Path) -> Result<AnetIndex, ArcsecError> {
    // Stub realloc: never called for a read-only fixed-size buffer (deltasize=0).
    unsafe extern "C" fn no_realloc(_p: *mut c_void, _n: usize) -> *mut c_void {
        core::ptr::null_mut()
    }

    // One bulk read — Vec stays alive for the entire CFITSIO session below.
    let mut file_bytes = std::fs::read(path).map_err(ArcsecError::CatalogIo)?;
    let mut buf_size = file_bytes.len();
    // Taken as *mut because that is what fits_open_memfile's signature wants; the
    // file is opened READONLY with deltasize = 0 and a no-op realloc, so CFITSIO
    // never writes through it or tries to grow it.
    let mut buf_ptr: *mut c_void = file_bytes.as_mut_ptr().cast::<c_void>();

    let name = cc(b"arcsec_index\0");

    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    fits_open_memfile(
        &mut fptr,
        name,
        READONLY,
        // buffptr is *mut *mut c_void: address of the data pointer.
        &raw mut buf_ptr,
        &mut buf_size,
        0, // deltasize = 0: never grow the read-only buffer
        no_realloc,
        &mut status,
    );
    if status != 0 {
        return Err(ArcsecError::CatalogIo(io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "fits_open_memfile failed for {}: status {status}",
                path.display()
            ),
        )));
    }

    let fp = fptr
        .as_deref_mut()
        .ok_or_else(|| ArcsecError::CatalogIo(std::io::Error::other("null fptr")))?;

    // Everything that can fail runs inside this closure so the handle is closed
    // exactly once, on success and on every error path. Previously a dozen `?`
    // sat between the open and the close, and `main` tries several index files
    // in turn, so a malformed one leaked a CFITSIO handle per attempt.
    let parsed = (|| -> Result<AnetParts, ArcsecError> {
        // ── Primary header metadata ───────────────────────────────────────────────
        let dim_quads = read_key_int(fp, b"DIMQUADS\0").max(3) as usize;
        let n_quads = read_key_int(fp, b"NQUADS\0") as usize;
        let n_stars = read_key_int(fp, b"NSTARS\0") as usize;
        let scale_lo = read_key_dbl(fp, b"SCALE_L\0").unwrap_or(0.0);
        let scale_hi = read_key_dbl(fp, b"SCALE_U\0").unwrap_or(PI);

        log::info!(
            "Anet index: DIMQUADS={}, {} quads, {} stars, scale {:.2}°–{:.2}°",
            dim_quads,
            n_quads,
            n_stars,
            scale_lo.to_degrees(),
            scale_hi.to_degrees()
        );

        if !(3..=4).contains(&dim_quads) {
            return Err(ArcsecError::CatalogIo(io::Error::new(
                io::ErrorKind::Unsupported,
                format!("DIMQUADS={dim_quads} not supported (only 3 or 4)"),
            )));
        }

        let n_code_dims = 2 * (dim_quads - 2); // 2 for triangles, 4 for quads

        let io_err =
            |msg: String| ArcsecError::CatalogIo(io::Error::new(io::ErrorKind::InvalidData, msg));

        // ── Code range parameters ─────────────────────────────────────────────────
        let (code_lo, code_scale) = read_code_range(fp, n_code_dims);
        log::debug!("Code range: lo={code_lo:.6}, scale={code_scale:.3e}");

        // ── Read star positions from "kdtree_data_stars" ──────────────────────────
        move_to_hdu(fp, b"kdtree_data_stars\0").map_err(&io_err)?;
        let n_star_rows = get_num_rows(fp).map_err(&io_err)?;
        let star_col = get_colnum(fp, b"kdtree_data_stars\0").map_err(&io_err)?;
        let star_bytes = read_raw_bytes(fp, star_col, n_star_rows * 12).map_err(&io_err)?;
        let stars = parse_stars(&star_bytes);

        // ── Read quad star indices from "quads" ───────────────────────────────────
        move_to_hdu(fp, b"quads\0").map_err(&io_err)?;
        let n_quad_rows = get_num_rows(fp).map_err(&io_err)?;
        let quad_col = get_colnum(fp, b"quads\0").map_err(&io_err)?;
        let row_bytes = dim_quads * 4;
        let quad_bytes = read_raw_bytes(fp, quad_col, n_quad_rows * row_bytes).map_err(&io_err)?;
        let quad_indices = parse_quad_indices(&quad_bytes);

        // ── Read codes from "kdtree_data_codes" ──────────────────────────────────
        move_to_hdu(fp, b"kdtree_data_codes\0").map_err(&io_err)?;
        let n_code_rows = get_num_rows(fp).map_err(&io_err)?;
        let code_col = get_colnum(fp, b"kdtree_data_codes\0").map_err(&io_err)?;
        let code_bytes =
            read_raw_bytes(fp, code_col, n_code_rows * n_code_dims * 2).map_err(&io_err)?;
        let codes = parse_codes(&code_bytes, n_code_dims, code_lo, code_scale);

        Ok(AnetParts {
            stars,
            quad_indices,
            codes,
            n_quad_rows,
            scale_lo,
            scale_hi,
            dim_quads,
        })
    })();

    // Close the in-memory FITS handle before dropping file_bytes, whatever happened.
    let mut cst: c_int = 0;
    if let Some(b) = fptr {
        fits_close_file(b, &mut cst);
    }
    let AnetParts {
        stars,
        quad_indices,
        codes,
        n_quad_rows,
        scale_lo,
        scale_hi,
        dim_quads,
    } = parsed?;

    // ── Build entry descriptors ───────────────────────────────────────────────
    let mut entries: Vec<AnetIndexEntry> = Vec::with_capacity(n_quad_rows);
    for (i, indices) in quad_indices.chunks_exact(dim_quads).enumerate() {
        let code = codes.get(i).copied().unwrap_or([0.0; 4]);
        if let Some(entry) = entry_from_sky(&stars, indices, code, dim_quads) {
            entries.push(entry);
        }
    }

    // Sort ascending on code[0] for binary search in find_code_matches.
    entries.sort_by(|a, b| a.code[0].total_cmp(&b.code[0]));

    log::info!("Anet index: {} entries loaded and sorted.", entries.len());

    // Build compact f32 code array for cache-efficient searching.
    let codes: Vec<[f32; 4]> = entries
        .iter()
        .map(|e| {
            [
                e.code[0] as f32,
                e.code[1] as f32,
                e.code[2] as f32,
                e.code[3] as f32,
            ]
        })
        .collect();

    Ok(AnetIndex {
        entries,
        codes,
        stars,
        scale_lo,
        scale_hi,
        dim_quads,
    })
}

/// Read only the primary header of an index file.
///
/// Returns `(dim_quads, scale_lo_rad, scale_hi_rad)` without loading any data tables.
/// Very fast — use this before `load_anet_index` to filter candidates.
///
/// # Errors
///
/// [`ArcsecError::CatalogIo`] if the path is not UTF-8 or CFITSIO cannot open it.
pub fn peek_anet_scale(path: &Path) -> Result<(usize, f64, f64), ArcsecError> {
    let path_str = path.to_str().ok_or_else(|| {
        ArcsecError::CatalogIo(io::Error::new(
            io::ErrorKind::InvalidInput,
            "non-UTF-8 path",
        ))
    })?;
    let cpath = format!("{path_str}\0");

    let mut fptr: Option<Box<fitsfile>> = None;
    let mut status: c_int = 0;
    fits_open_image(&mut fptr, cc(cpath.as_bytes()), READONLY, &mut status);
    if status != 0 {
        return Err(ArcsecError::CatalogIo(io::Error::new(
            io::ErrorKind::NotFound,
            format!("cannot open {}", path.display()),
        )));
    }

    let fp = fptr
        .as_deref_mut()
        .ok_or_else(|| ArcsecError::CatalogIo(std::io::Error::other("null fptr")))?;

    let dim_quads = read_key_int(fp, b"DIMQUADS\0").max(3) as usize;
    let scale_lo = read_key_dbl(fp, b"SCALE_L\0").unwrap_or(0.0);
    let scale_hi = read_key_dbl(fp, b"SCALE_U\0").unwrap_or(PI);

    let mut cst: c_int = 0;
    if let Some(b) = fptr {
        fits_close_file(b, &mut cst);
    }

    Ok((dim_quads, scale_lo, scale_hi))
}

// ── Unit tests ─────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_index_entry(cx: f64, cy: f64) -> AnetIndexEntry {
        AnetIndexEntry {
            code: [cx, cy, 0.0, 0.0],
            n_stars: 3,
            star_ra: [0.0; 4],
            star_dec: [0.0; 4],
            center_ra: 0.0,
            center_dec: 0.0,
        }
    }

    #[test]
    fn u32_to_unit_endpoints() {
        assert!((u32_to_unit(0) - (-1.0)).abs() < 1e-9);
        assert!((u32_to_unit(u32::MAX) - 1.0).abs() < 1e-9);
        assert!((u32_to_unit(u32::MAX / 2) - 0.0).abs() < 0.01);
    }

    #[test]
    fn find_code_matches_basic() {
        let entries = vec![make_index_entry(0.5, 0.3), make_index_entry(0.8, 0.4)];
        let codes: Vec<[f32; 4]> = entries
            .iter()
            .map(|e| {
                [
                    e.code[0] as f32,
                    e.code[1] as f32,
                    e.code[2] as f32,
                    e.code[3] as f32,
                ]
            })
            .collect();
        let index = AnetIndex {
            entries,
            codes,
            stars: vec![],
            scale_lo: 0.0,
            scale_hi: PI,
            dim_quads: 3,
        };
        let hits = index.find_code_matches(&[0.5, 0.3, 0.0, 0.0], 0.02);
        assert_eq!(hits.len(), 1, "expected 1 match, got {}", hits.len());
        let no_hits = index.find_code_matches(&[0.5, 0.6, 0.0, 0.0], 0.02);
        assert_eq!(no_hits.len(), 0);
    }

    // ── Loading real FITS files ────────────────────────────────────────────────

    use crate::test_support::{
        ANET_CODE_HI, ANET_CODE_LO, FitsWriter, RawIndex, Rng, SkySpec, TempDir, random_sky,
        separation,
    };

    /// A small index: 300 stars around (40°, -25°), `n_quads` random groups.
    fn raw_index(dim_quads: usize, n_quads: usize, seed: u64) -> RawIndex {
        let mut rng = Rng::new(seed);
        let sky: Vec<(f64, f64)> = random_sky(
            &mut rng,
            &SkySpec {
                ra0: 40f64.to_radians(),
                dec0: (-25f64).to_radians(),
                side_deg: 4.0,
                n: 300,
                min_sep_deg: 0.01,
                mag_lo: 8.0,
                mag_hi: 12.0,
            },
        )
        .iter()
        .map(|s| (s.ra, s.dec))
        .collect();
        let groups: Vec<Vec<u32>> = (0..n_quads)
            .map(|_| {
                let mut g: Vec<u32> = Vec::new();
                while g.len() < dim_quads {
                    let i = (rng.next_u64() % sky.len() as u64) as u32;
                    if !g.contains(&i) {
                        g.push(i);
                    }
                }
                g
            })
            .collect();
        RawIndex::build(dim_quads, &sky, &groups)
    }

    fn write(dir: &TempDir, name: &str, bytes: &[u8]) -> std::path::PathBuf {
        let path = dir.path().join(name);
        std::fs::write(&path, bytes).unwrap();
        path
    }

    fn assert_loads_faithfully(raw: &RawIndex) {
        let dir = TempDir::new("anet");
        let path = write(&dir, "index.fits", &raw.fits_bytes());
        let idx = load_anet_index(&path).expect("load");

        assert_eq!(idx.dim_quads, raw.dim_quads);
        assert_eq!(idx.n_code_dims(), 2 * (raw.dim_quads - 2));
        assert!((idx.scale_lo - raw.scale_lo).abs() < 1e-12);
        assert!((idx.scale_hi - raw.scale_hi).abs() < 1e-12);

        // Stars: 32-bit fixed-point unit vectors, good to a few milliarcseconds.
        assert_eq!(idx.stars.len(), raw.stars.len());
        for (s, &(ra, dec)) in idx.stars.iter().zip(&raw.stars) {
            assert!(separation(s.ra, s.dec, ra, dec) < 1e-8, "star moved");
            assert!((0.0..2.0 * PI).contains(&s.ra));
        }

        // Entries come back sorted on code[0], with the parallel f32 array.
        assert_eq!(idx.entries.len(), raw.quads.len());
        assert_eq!(idx.codes.len(), idx.entries.len());
        assert!(idx.entries.windows(2).all(|w| w[0].code[0] <= w[1].code[0]));
        for (e, c) in idx.entries.iter().zip(&idx.codes) {
            assert!((0..4).all(|k| (e.code[k] as f32 - c[k]).abs() < 1e-6));
        }

        // Every quad is present with its stars in canonical order and its code
        // decoded to within the u16 quantum.
        let quantum = (ANET_CODE_HI - ANET_CODE_LO) / 65535.0;
        let expected = raw.to_index();
        for want in &expected.entries {
            let got = idx
                .entries
                .iter()
                .find(|e| {
                    (0..want.n_stars).all(|k| {
                        separation(
                            e.star_ra[k],
                            e.star_dec[k],
                            want.star_ra[k],
                            want.star_dec[k],
                        ) < 1e-8
                    })
                })
                .expect("quad lost in loading");
            assert_eq!(got.n_stars, raw.dim_quads);
            for k in 0..idx.n_code_dims() {
                assert!(
                    (got.code[k] - want.code[k]).abs() <= quantum,
                    "code[{k}] {} vs {}",
                    got.code[k],
                    want.code[k]
                );
            }
            assert!(
                separation(
                    got.center_ra,
                    got.center_dec,
                    want.center_ra,
                    want.center_dec
                ) < 1e-8
            );
        }

        // A lookup by one of its own codes finds the entry.
        let probe = &idx.entries[idx.entries.len() / 2];
        let hits = idx.find_code_matches(&probe.code, 1e-6);
        assert!(hits.iter().any(|&i| idx.entries[i].code == probe.code));

        // The header peek agrees with the full load.
        let (dq, lo, hi) = peek_anet_scale(&path).expect("peek");
        assert_eq!(dq, raw.dim_quads);
        assert!((lo - raw.scale_lo).abs() < 1e-12 && (hi - raw.scale_hi).abs() < 1e-12);
    }

    #[test]
    fn loads_a_quad_index() {
        assert_loads_faithfully(&raw_index(4, 400, 1));
    }

    #[test]
    fn loads_a_triangle_index() {
        assert_loads_faithfully(&raw_index(3, 400, 2));
    }

    #[test]
    fn find_code_matches_uses_every_code_dimension() {
        let idx = raw_index(4, 400, 3).to_index();
        let target = idx.entries[123].code;
        // Brute force over the entries must agree with the binary-searched scan.
        for tol in [0.001, 0.01, 0.05] {
            for probe in [target, [0.3, 0.4, 0.5, 0.6], [0.0, 0.0, 1.0, 1.0]] {
                let mut want: Vec<usize> = idx
                    .entries
                    .iter()
                    .enumerate()
                    .filter(|(_, e)| {
                        let d2: f64 = (0..4).map(|k| (e.code[k] - probe[k]).powi(2)).sum();
                        d2.sqrt() <= tol * 0.999
                    })
                    .map(|(i, _)| i)
                    .collect();
                let mut got = idx.find_code_matches(&probe, tol);
                want.sort_unstable();
                got.sort_unstable();
                // f32 rounding may admit an entry sitting exactly on the boundary.
                assert!(want.iter().all(|i| got.contains(i)), "tol {tol}");
                assert!(got.len() <= want.len() + 1, "tol {tol}");
            }
        }
    }

    #[test]
    fn a_missing_range_table_falls_back_to_the_standard_code_range() {
        // Rebuild the file without kdtree_range_codes: the loader must use the
        // standard astrometry.net range, which is what the writer used anyway.
        let raw = raw_index(4, 50, 4);
        let bytes = raw.fits_bytes();
        let dir = TempDir::new("anet-norange");
        let mut w = FitsWriter::default();
        copy_hdus_except(&bytes, "kdtree_range_codes", &mut w);
        let path = write(&dir, "index.fits", &w.bytes);
        let idx = load_anet_index(&path).expect("load");
        let want = raw.to_index();
        for (a, b) in idx.entries.iter().zip(&want.entries) {
            assert!((a.code[0] - b.code[0]).abs() < 1e-4);
        }
    }

    /// Copy every HDU of a FITS byte stream except the table called `skip`.
    fn copy_hdus_except(bytes: &[u8], skip: &str, out: &mut FitsWriter) {
        let mut pos = 0;
        while pos < bytes.len() {
            let start = pos;
            let mut naxis1 = 0usize;
            let mut naxis2 = 0usize;
            let mut name = String::new();
            let mut is_ext = false;
            loop {
                let card = core::str::from_utf8(&bytes[pos..pos + 80]).unwrap();
                pos += 80;
                let val = card
                    .get(10..)
                    .unwrap_or("")
                    .trim()
                    .trim_matches('\'')
                    .trim();
                match card.get(..8).unwrap_or("").trim() {
                    "XTENSION" => is_ext = true,
                    "NAXIS1" => naxis1 = val.parse().unwrap(),
                    "NAXIS2" => naxis2 = val.parse().unwrap(),
                    "TTYPE1" => name = val.to_string(),
                    "END" => break,
                    _ => {}
                }
            }
            pos = pos.next_multiple_of(2880);
            if is_ext {
                pos += (naxis1 * naxis2).next_multiple_of(2880);
            }
            if name != skip {
                out.bytes.extend_from_slice(&bytes[start..pos]);
            }
        }
    }

    fn err_kind(r: Result<AnetIndex, ArcsecError>) -> io::ErrorKind {
        match r {
            Err(ArcsecError::CatalogIo(e)) => e.kind(),
            Err(e) => panic!("expected CatalogIo, got {e:?}"),
            Ok(_) => panic!("expected an error"),
        }
    }

    /// `peek_anet_scale` is documented to return `CatalogIo` when CFITSIO cannot
    /// open the file, but rsfitsio's `fits_open_image` (`ffiopn_safer`,
    /// cfileio.rs:887) unwraps the file handle after a failed open and panics with
    /// "Null Pointer" instead of returning the status. The CLI's
    /// `collect_index_files` peeks every `index-*.fits` in the directory and means
    /// to skip bad ones with `.ok()?`, so one empty or corrupt index file (an
    /// interrupted download, say) crashes the whole run.
    #[test]
    #[ignore = "bug: peek_anet_scale panics inside rsfitsio on a missing, empty or non-FITS file"]
    fn peek_anet_scale_reports_unreadable_files() {
        let dir = TempDir::new("anet-peek");
        let missing = dir.path().join("absent.fits");
        let junk = write(&dir, "junk.fits", &[0x5Au8; 4000]);
        let empty = write(&dir, "empty.fits", &[]);
        for path in [missing, junk, empty] {
            let r = std::panic::catch_unwind(|| peek_anet_scale(&path));
            assert!(
                matches!(r, Ok(Err(ArcsecError::CatalogIo(_)))),
                "{} did not give CatalogIo",
                path.display()
            );
        }
    }

    /// An index whose final HDU is one the loader needs cannot be loaded: moving to
    /// the last HDU of the in-memory file fails with CFITSIO status 113
    /// (`READ_ERROR`), where the same move on the same bytes opened from disk
    /// succeeds. Distributed astrometry.net indexes end with tables the loader never
    /// reads (sweep, magnitudes), so this does not bite today, but an index trimmed
    /// to the tables arcsec uses fails with a misleading "HDU … not found".
    #[test]
    #[ignore = "bug: load_anet_index cannot reach the last HDU of an in-memory FITS file"]
    fn load_anet_index_reads_a_required_table_in_the_last_hdu() {
        let raw = raw_index(4, 50, 8);
        let mut w = FitsWriter::default();
        copy_hdus_except(&raw.fits_bytes(), "sweep", &mut w);
        let dir = TempDir::new("anet-last");
        let path = write(&dir, "index.fits", &w.bytes);
        let idx = load_anet_index(&path).expect("kdtree_data_stars is the last HDU");
        assert_eq!(idx.stars.len(), raw.stars.len());
    }

    #[test]
    fn unreadable_indexes_are_errors_not_panics() {
        let dir = TempDir::new("anet-bad");

        // No such file.
        let missing = dir.path().join("absent.fits");
        assert_eq!(err_kind(load_anet_index(&missing)), io::ErrorKind::NotFound);

        // Not FITS at all.
        let junk = write(&dir, "junk.fits", &[0x5Au8; 4000]);
        assert!(matches!(
            load_anet_index(&junk),
            Err(ArcsecError::CatalogIo(_))
        ));

        // An empty file.
        let empty = write(&dir, "empty.fits", &[]);
        assert!(matches!(
            load_anet_index(&empty),
            Err(ArcsecError::CatalogIo(_))
        ));

        // DIMQUADS outside 3..=4. (Each file ends with an unused table: the loader
        // cannot read a file's last HDU, see the ignored test above.)
        let mut w = FitsWriter::default();
        w.primary(&[("DIMQUADS", "5".into())]);
        w.table("sweep", 1, &[0]);
        let five = write(&dir, "dim5.fits", &w.bytes);
        assert_eq!(err_kind(load_anet_index(&five)), io::ErrorKind::Unsupported);

        // A valid primary header but none of the tables.
        let mut w = FitsWriter::default();
        w.primary(&[("DIMQUADS", "4".into())]);
        w.table("sweep", 1, &[0]);
        let bare = write(&dir, "bare.fits", &w.bytes);
        assert_eq!(err_kind(load_anet_index(&bare)), io::ErrorKind::InvalidData);

        // Each required table missing in turn.
        let full = raw_index(4, 20, 5).fits_bytes();
        for table in ["kdtree_data_stars", "quads", "kdtree_data_codes"] {
            let mut w = FitsWriter::default();
            copy_hdus_except(&full, table, &mut w);
            let path = write(&dir, &format!("no-{table}.fits"), &w.bytes);
            assert_eq!(
                err_kind(load_anet_index(&path)),
                io::ErrorKind::InvalidData,
                "without {table}"
            );
        }
    }

    #[test]
    fn quads_pointing_at_missing_or_coincident_stars_are_dropped() {
        let stars = vec![
            AnetStar { ra: 1.0, dec: 0.5 },
            AnetStar { ra: 1.01, dec: 0.5 },
            AnetStar { ra: 1.0, dec: 0.51 },
            AnetStar { ra: 1.0, dec: 0.5 }, // same place as star 0
        ];
        let code = [0.1, 0.2, 0.3, 0.4];
        assert!(entry_from_sky(&stars, &[0, 1, 2], code, 3).is_some());
        assert!(entry_from_sky(&stars, &[0, 1, 9], code, 3).is_none());
        assert!(entry_from_sky(&stars, &[0, 3, 1], code, 3).is_none());
        assert!(entry_from_sky(&stars, &[0, 1], code, 2).is_none());
        assert!(entry_from_sky(&stars, &[0, 1, 2, 3, 0], code, 5).is_none());

        // Loading a file with such a quad keeps the good ones.
        let mut raw = raw_index(4, 30, 6);
        raw.quads[3] = vec![0, 1, 2, 99_999];
        let dir = TempDir::new("anet-dangling");
        let path = write(&dir, "index.fits", &raw.fits_bytes());
        let idx = load_anet_index(&path).unwrap();
        assert_eq!(idx.entries.len(), raw.quads.len() - 1);
    }

    #[test]
    fn star_decoding_covers_the_whole_sphere() {
        // The poles, the equator either side of RA 0, and a southern star.
        let pts = [
            (0.0, PI / 2.0),
            (3.0, -PI / 2.0),
            (0.0, 0.0),
            (2.0 * PI - 1e-6, 0.0),
            (4.0, -1.2),
        ];
        let mut bytes = Vec::new();
        for &(ra, dec) in &pts {
            for v in [dec.cos() * ra.cos(), dec.cos() * ra.sin(), dec.sin()] {
                let u = ((v + 1.0) * (f64::from(u32::MAX) / 2.0)).round() as u32;
                bytes.extend_from_slice(&u.to_le_bytes());
            }
        }
        let stars = parse_stars(&bytes);
        assert_eq!(stars.len(), pts.len());
        for (s, &(ra, dec)) in stars.iter().zip(&pts) {
            assert!((s.dec - dec).abs() < 1e-8, "dec {} vs {dec}", s.dec);
            if dec.abs() < 1.5 {
                assert!(separation(s.ra, s.dec, ra, dec) < 1e-8);
            }
        }
        // A trailing partial record is ignored rather than misread.
        bytes.extend_from_slice(&[1, 2, 3]);
        assert_eq!(parse_stars(&bytes).len(), pts.len());
    }

    #[test]
    fn code_parsing_decodes_every_dimension() {
        let lo = -0.25;
        let scale = 1000.0;
        let mut bytes = Vec::new();
        for u in [0u16, 250, 500, 65_535] {
            bytes.extend_from_slice(&u.to_le_bytes());
        }
        let four = parse_codes(&bytes, 4, lo, scale);
        assert_eq!(four.len(), 1);
        assert_eq!(four[0], [-0.25, 0.0, 0.25, 65.285]);
        let two = parse_codes(&bytes, 2, lo, scale);
        assert_eq!(two, vec![[-0.25, 0.0, 0.0, 0.0], [0.25, 65.285, 0.0, 0.0]]);
        assert_eq!(
            parse_quad_indices(&[1, 0, 0, 0, 2, 1, 0, 0, 9]),
            vec![1, 258]
        );
    }
}
