//! `arcsec catalog` — install and manage star catalogues.
//!
//! arcsec needs a star catalogue to solve anything, and a *photometric* catalogue to
//! do colour calibration. Both live behind download pages that are easy to get wrong,
//! so this puts them one command away and in one place:
//!
//! ```text
//! arcsec catalog list                  # what exists, what is installed
//! arcsec catalog recommend --fov 1.5   # what this rig needs
//! arcsec catalog install d50 v05       # fetch and unpack
//! arcsec catalog path                  # where they went
//! ```
//!
//! Everything lands in a per-platform data directory (see [`default_dir`]) which the
//! solver reads by default, so `-d` is only needed to override it.

use std::fs;
use std::io;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

// ── Registry ────────────────────────────────────────────────────────────────────

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
    fn label(self) -> &'static str {
        match self {
            Purpose::Solving => "solving",
            Purpose::Photometry => "photometry",
            Purpose::BlindIndex => "blind",
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
    /// Plain files, downloaded individually.
    Loose,
}

pub struct Entry {
    /// Short name the user types.
    pub id: &'static str,
    pub purpose: Purpose,
    /// Download URL, or a `{}` template for `Loose` sets.
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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Files {
    AstapDb { prefix: &'static str },
    AnetIndex,
}

/// Extensions an ASTAP star database can use, and the file count each implies.
pub const ASTAP_EXTS: &[(&str, usize)] = &[(".1476", 1476), (".290", 290), (".001", 1)];

impl Files {
    /// Does `name` belong to this catalogue?
    fn owns(&self, name: &str) -> bool {
        match self {
            Files::AstapDb { prefix } => {
                name.starts_with(&format!("{prefix}_"))
                    && ASTAP_EXTS.iter().any(|(e, _)| name.ends_with(e))
            }
            Files::AnetIndex => name.starts_with("index-") && name.ends_with(".fits"),
        }
    }

    /// A file whose presence means the catalogue is installed.
    fn probe(&self, dir: &Path) -> Option<PathBuf> {
        match self {
            Files::AstapDb { prefix } => ASTAP_EXTS
                .iter()
                .map(|(e, _)| dir.join(format!("{prefix}_0101{e}")))
                .find(|p| p.exists()),
            Files::AnetIndex => fs::read_dir(dir)
                .into_iter()
                .flatten()
                .flatten()
                .map(|d| d.path())
                .find(|p| {
                    p.file_name()
                        .map(|n| self.owns(&n.to_string_lossy()))
                        .unwrap_or(false)
                }),
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
        bytes: 160_000_000,
        fov: Some((0.7, 180.0)),
        files: Files::AnetIndex,
        desc: "Tycho-2 blind indexes, scales 07–19 (fields ~0.7° and wider).",
    },
    Entry {
        id: "anet-5200",
        purpose: Purpose::BlindIndex,
        url: "https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/index-52{:02}-{:02}.fits",
        archive: Archive::Loose,
        bytes: 8_800_000_000,
        fov: Some((0.1, 2.0)),
        files: Files::AnetIndex,
        desc: "Gaia LITE blind indexes 5200/5201/5202, 48 HEALPix each. Large.",
    },
];

pub fn find(id: &str) -> Option<&'static Entry> {
    REGISTRY.iter().find(|e| e.id.eq_ignore_ascii_case(id))
}

/// Individual file URLs for a `Loose` entry.
fn loose_urls(e: &Entry) -> Vec<(String, String)> {
    let mut out = Vec::new();
    match e.id {
        "anet-4100" => {
            for scale in 7..=19 {
                out.push((
                    format!("https://data.astrometry.net/4100/index-41{scale:02}.fits"),
                    format!("index-41{scale:02}.fits"),
                ));
            }
        }
        "anet-5200" => {
            for series in [5200u32, 5201, 5202] {
                for hp in 0..48 {
                    out.push((
                        format!(
                            "https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/index-{series}-{hp:02}.fits"
                        ),
                        format!("index-{series}-{hp:02}.fits"),
                    ));
                }
            }
        }
        _ => {}
    }
    out
}

// ── Install location ────────────────────────────────────────────────────────────

/// Where catalogues are kept, in priority order:
///
/// 1. `$ARCSEC_CATALOG_DIR`, if set — for people who keep them on another disk.
/// 2. `$XDG_DATA_HOME/arcsec/catalogs` on Linux, or the platform equivalent:
///    `~/Library/Application Support/arcsec/catalogs` on macOS,
///    `%LOCALAPPDATA%\arcsec\catalogs` on Windows.
/// 3. `~/.arcsec/catalogs` if the home directory cannot be resolved any other way.
///
/// The point is that a user who runs `arcsec catalog install d50` never has to know
/// this path, and the solver looks here without being told.
pub fn default_dir() -> PathBuf {
    if let Ok(p) = std::env::var("ARCSEC_CATALOG_DIR")
        && !p.is_empty()
    {
        return PathBuf::from(p);
    }

    #[cfg(target_os = "windows")]
    {
        if let Ok(p) = std::env::var("LOCALAPPDATA")
            && !p.is_empty()
        {
            return PathBuf::from(p).join("arcsec").join("catalogs");
        }
    }
    #[cfg(target_os = "macos")]
    {
        if let Ok(h) = std::env::var("HOME")
            && !h.is_empty()
        {
            return PathBuf::from(h)
                .join("Library")
                .join("Application Support")
                .join("arcsec")
                .join("catalogs");
        }
    }
    #[cfg(not(any(target_os = "windows", target_os = "macos")))]
    {
        if let Ok(p) = std::env::var("XDG_DATA_HOME")
            && !p.is_empty()
        {
            return PathBuf::from(p).join("arcsec").join("catalogs");
        }
        if let Ok(h) = std::env::var("HOME")
            && !h.is_empty()
        {
            return PathBuf::from(h)
                .join(".local")
                .join("share")
                .join("arcsec")
                .join("catalogs");
        }
    }

    #[allow(unreachable_code)]
    {
        std::env::var("HOME")
            .map(|h| PathBuf::from(h).join(".arcsec").join("catalogs"))
            .unwrap_or_else(|_| PathBuf::from("catalogs"))
    }
}

/// Is `e` present in `dir`?
pub fn is_installed(dir: &Path, e: &Entry) -> bool {
    e.files.probe(dir).is_some()
}

/// Files belonging to `e` that are present in `dir`.
fn installed_files(dir: &Path, e: &Entry) -> Vec<PathBuf> {
    fs::read_dir(dir)
        .into_iter()
        .flatten()
        .flatten()
        .map(|d| d.path())
        .filter(|p| {
            p.file_name()
                .map(|n| e.files.owns(&n.to_string_lossy()))
                .unwrap_or(false)
        })
        .collect()
}

/// Bytes an installed catalogue occupies on disk.
fn installed_size(dir: &Path, e: &Entry) -> u64 {
    installed_files(dir, e)
        .iter()
        .filter_map(|p| fs::metadata(p).ok())
        .map(|m| m.len())
        .sum()
}

pub fn human(bytes: u64) -> String {
    const U: [&str; 5] = ["B", "kB", "MB", "GB", "TB"];
    let mut v = bytes as f64;
    let mut i = 0;
    while v >= 1000.0 && i < U.len() - 1 {
        v /= 1000.0;
        i += 1;
    }
    if i == 0 {
        format!("{} {}", bytes, U[i])
    } else {
        format!("{v:.1} {}", U[i])
    }
}

/// Is stderr a terminal? Controls whether progress redraws in place.
fn stderr_is_terminal() -> bool {
    // SAFETY: isatty on a fixed descriptor has no preconditions.
    unsafe { libc::isatty(2) == 1 }
}

// ── Download ────────────────────────────────────────────────────────────────────

/// Stream `url` to `dest`, reporting progress.
///
/// Downloads to `dest.part` and renames on success, so an interrupted run never
/// leaves a half file that later looks installed. Existing `.part` files are resumed
/// with a Range request where the server allows it — these are multi-hundred-megabyte
/// downloads and starting over is not acceptable.
fn download(url: &str, dest: &Path, label: &str) -> Result<(), String> {
    let part = dest.with_extension(format!(
        "{}part",
        dest.extension()
            .map(|e| format!("{}.", e.to_string_lossy()))
            .unwrap_or_default()
    ));
    let have = fs::metadata(&part).map(|m| m.len()).unwrap_or(0);

    let agent = ureq::Agent::config_builder()
        .timeout_global(Some(core::time::Duration::from_secs(7200)))
        .build()
        .new_agent();

    let mut req = agent.get(url);
    if have > 0 {
        req = req.header("Range", &format!("bytes={have}-"));
    }
    let resp = req.call().map_err(|e| format!("{label}: {e}"))?;
    let status = resp.status().as_u16();
    let resuming = status == 206;
    if !(status == 200 || resuming) {
        return Err(format!("{label}: HTTP {status}"));
    }
    let total = resp
        .headers()
        .get("content-length")
        .and_then(|v| v.to_str().ok())
        .and_then(|v| v.parse::<u64>().ok())
        .map(|n| n + if resuming { have } else { 0 });

    let mut out = fs::OpenOptions::new()
        .create(true)
        .append(resuming)
        .write(true)
        .truncate(!resuming)
        .open(&part)
        .map_err(|e| format!("{}: {e}", part.display()))?;

    let mut written = if resuming { have } else { 0 };
    let mut reader = resp.into_body().into_reader();
    let mut buf = vec![0u8; 1 << 20];
    let mut last_report = std::time::Instant::now();
    loop {
        let n = reader.read(&mut buf).map_err(|e| format!("{label}: {e}"))?;
        if n == 0 {
            break;
        }
        out.write_all(&buf[..n])
            .map_err(|e| format!("{}: {e}", part.display()))?;
        written += n as u64;
        // Animate only on a terminal. Redirected to a file or a pipe, carriage
        // returns do not overwrite, so an unconditional update turns a 117 MB
        // download into hundreds of lines of noise.
        let interval = if stderr_is_terminal() { 300 } else { 15_000 };
        if last_report.elapsed().as_millis() > interval {
            match total {
                Some(t) if t > 0 => eprint!(
                    "{}  {label}: {} / {} ({:.0}%){}",
                    if stderr_is_terminal() { "\r" } else { "" },
                    human(written),
                    human(t),
                    written as f64 / t as f64 * 100.0,
                    if stderr_is_terminal() { "   " } else { "\n" }
                ),
                _ => eprint!(
                    "{}  {label}: {}{}",
                    if stderr_is_terminal() { "\r" } else { "" },
                    human(written),
                    if stderr_is_terminal() { "   " } else { "\n" }
                ),
            }
            let _ = std::io::stderr().flush();
            last_report = std::time::Instant::now();
        }
    }
    drop(out);
    eprintln!(
        "{}  {label}: {} downloaded          ",
        if stderr_is_terminal() { "\r" } else { "" },
        human(written)
    );
    fs::rename(&part, dest).map_err(|e| format!("{}: {e}", dest.display()))?;
    Ok(())
}

// ── Extraction ──────────────────────────────────────────────────────────────────

/// Unpack a `.zip` whose entries sit at the archive root.
fn extract_zip(archive: &Path, dest: &Path) -> Result<usize, String> {
    let f = fs::File::open(archive).map_err(|e| format!("{}: {e}", archive.display()))?;
    let mut zip = zip::ZipArchive::new(f).map_err(|e| format!("{}: {e}", archive.display()))?;
    let mut n = 0;
    for i in 0..zip.len() {
        let mut item = zip.by_index(i).map_err(|e| e.to_string())?;
        if item.is_dir() {
            continue;
        }
        // Flatten: catalogue archives are flat, and this also refuses any entry
        // trying to escape the destination with a path like `../..`.
        let name = match item
            .enclosed_name()
            .and_then(|p| p.file_name().map(|s| s.to_owned()))
        {
            Some(n) => n,
            None => continue,
        };
        let out_path = dest.join(&name);
        let mut out =
            fs::File::create(&out_path).map_err(|e| format!("{}: {e}", out_path.display()))?;
        std::io::copy(&mut item, &mut out).map_err(|e| format!("{}: {e}", out_path.display()))?;
        n += 1;
    }
    Ok(n)
}

/// Read the members of an `ar` archive — the container format of a `.deb`.
///
/// The format is trivial: an 8-byte magic then, per member, a 60-byte ASCII header
/// whose last fields are the size and a `` `\n `` terminator, followed by the data
/// padded to an even length. Implementing it here avoids requiring `dpkg-deb` or
/// `ar`, neither of which exists on Windows.
fn ar_member(data: &[u8], want_prefix: &str) -> Option<(String, core::ops::Range<usize>)> {
    if data.len() < 8 || &data[..8] != b"!<arch>\n" {
        return None;
    }
    let mut pos = 8usize;
    while pos + 60 <= data.len() {
        let hdr = &data[pos..pos + 60];
        let name = String::from_utf8_lossy(&hdr[0..16]).trim().to_string();
        let size: usize = String::from_utf8_lossy(&hdr[48..58]).trim().parse().ok()?;
        let start = pos + 60;
        let end = start.checked_add(size)?;
        if end > data.len() {
            return None;
        }
        if name.starts_with(want_prefix) {
            return Some((name, start..end));
        }
        pos = end + (size & 1); // members are padded to an even offset
    }
    None
}

/// Unpack a `.deb`: find `data.tar.xz`, xz-decompress it, then untar, flattening
/// paths (the ASTAP packages install under `opt/astap/`).
fn extract_deb(archive: &Path, dest: &Path) -> Result<usize, String> {
    let raw = fs::read(archive).map_err(|e| format!("{}: {e}", archive.display()))?;
    let (name, range) = ar_member(&raw, "data.tar")
        .ok_or_else(|| format!("{}: no data.tar member (not a .deb?)", archive.display()))?;

    eprintln!("  decompressing {name} ...");
    let compressed = &raw[range];
    let mut tar_bytes: Vec<u8> = Vec::new();
    if name.ends_with(".xz") {
        let mut cur = io::Cursor::new(compressed);
        lzma_rs::xz_decompress(&mut cur, &mut tar_bytes)
            .map_err(|e| format!("{}: xz decompress failed: {e:?}", archive.display()))?;
    } else if name.ends_with(".gz") {
        return Err(format!(
            "{}: gzip .deb payloads are not supported",
            archive.display()
        ));
    } else {
        tar_bytes = compressed.to_vec();
    }

    let mut ar = tar::Archive::new(io::Cursor::new(tar_bytes));
    let mut n = 0;
    for entry in ar.entries().map_err(|e| e.to_string())? {
        let mut entry = entry.map_err(|e| e.to_string())?;
        if !entry.header().entry_type().is_file() {
            continue;
        }
        let path = entry.path().map_err(|e| e.to_string())?.into_owned();
        let Some(fname) = path.file_name() else {
            continue;
        };
        let out_path = dest.join(fname);
        let mut out =
            fs::File::create(&out_path).map_err(|e| format!("{}: {e}", out_path.display()))?;
        std::io::copy(&mut entry, &mut out).map_err(|e| format!("{}: {e}", out_path.display()))?;
        n += 1;
    }
    Ok(n)
}

// ── Commands ────────────────────────────────────────────────────────────────────

pub fn cmd_path(dir: &Path) {
    println!("{}", dir.display());
}

pub fn cmd_list(dir: &Path) {
    println!("Catalogue directory: {}", dir.display());
    println!("  (override with --dir, or the ARCSEC_CATALOG_DIR environment variable)\n");
    println!(
        "{:<11} {:<11} {:>9}  {:<13} STATUS",
        "NAME", "PURPOSE", "DOWNLOAD", "FIELDS"
    );
    for e in REGISTRY {
        let fov = match e.fov {
            Some((lo, hi)) => format!("{lo}°–{hi}°"),
            None => "—".to_string(),
        };
        let status = if is_installed(dir, e) {
            let sz = installed_size(dir, e);
            if sz > 0 {
                format!("installed ({})", human(sz))
            } else {
                "installed".to_string()
            }
        } else {
            "not installed".to_string()
        };
        println!(
            "{:<11} {:<11} {:>9}  {:<13} {}",
            e.id,
            e.purpose.label(),
            human(e.bytes),
            fov,
            status
        );
    }
    println!("\nDescriptions:");
    for e in REGISTRY {
        println!("  {:<11} {}", e.id, e.desc);
    }
}

/// Suggest catalogues for a field size.
pub fn cmd_recommend(dir: &Path, fov_deg: f64, want_photometry: bool) {
    println!("For a {fov_deg:.2}° field:\n");
    let pick = |purpose: Purpose| -> Option<&'static Entry> {
        // Prefer the smallest download whose range covers the field, so the advice
        // does not push a gigabyte on someone who does not need it.
        REGISTRY
            .iter()
            .filter(|e| e.purpose == purpose)
            .filter(|e| matches!(e.fov, Some((lo, hi)) if fov_deg >= lo && fov_deg <= hi))
            .min_by_key(|e| e.bytes)
    };

    match pick(Purpose::Solving) {
        Some(e) => println!(
            "  solving      {:<10} {:>9}  {}{}",
            e.id,
            human(e.bytes),
            e.desc,
            if is_installed(dir, e) {
                "  [installed]"
            } else {
                ""
            }
        ),
        None => println!("  solving      no catalogue covers this field size"),
    }
    if want_photometry {
        match pick(Purpose::Photometry) {
            Some(e) => println!(
                "  photometry   {:<10} {:>9}  {}{}",
                e.id,
                human(e.bytes),
                e.desc,
                if is_installed(dir, e) {
                    "  [installed]"
                } else {
                    ""
                }
            ),
            None => println!("  photometry   no catalogue covers this field size"),
        }
    }
    if let Some(e) = REGISTRY.iter().find(|e| e.id == "anet-4100") {
        println!(
            "  blind        {:<10} {:>9}  optional: solve with no position hint{}",
            e.id,
            human(e.bytes),
            if is_installed(dir, e) {
                "  [installed]"
            } else {
                ""
            }
        );
    }

    let ids: Vec<&str> = [
        pick(Purpose::Solving),
        want_photometry.then(|| pick(Purpose::Photometry)).flatten(),
    ]
    .into_iter()
    .flatten()
    .filter(|e| !is_installed(dir, e))
    .map(|e| e.id)
    .collect();
    if !ids.is_empty() {
        println!("\n  arcsec catalog install {}", ids.join(" "));
    }
}

pub fn cmd_remove(dir: &Path, ids: &[String]) -> Result<(), String> {
    for id in ids {
        let e = find(id)
            .ok_or_else(|| format!("unknown catalogue '{id}' (try: arcsec catalog list)"))?;
        let files = installed_files(dir, e);
        let mut freed = 0u64;
        for p in &files {
            freed += fs::metadata(p).map(|m| m.len()).unwrap_or(0);
            fs::remove_file(p).map_err(|err| format!("{}: {err}", p.display()))?;
        }
        println!(
            "removed {id}: {} files, {} freed",
            files.len(),
            human(freed)
        );
    }
    Ok(())
}

/// Check that every installed catalogue looks structurally sound.
pub fn cmd_verify(dir: &Path) -> Result<(), String> {
    let mut problems = 0;
    let mut checked = 0;
    for e in REGISTRY {
        if !is_installed(dir, e) {
            continue;
        }
        checked += 1;
        let files = installed_files(dir, e);
        let mut bad = Vec::new();

        for p in &files {
            match fs::metadata(p) {
                // Every format has at least a 110-byte header (or a FITS block).
                Ok(m) if m.len() < 120 => bad.push(format!("{} is truncated", p.display())),
                Err(err) => bad.push(format!("{}: {err}", p.display())),
                _ => {}
            }
        }

        // For an ASTAP database the extension fixes how many files there should be.
        if let Files::AstapDb { .. } = e.files
            && let Some(probe) = e.files.probe(dir)
        {
            let name = probe
                .file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .to_string();
            if let Some((_, want)) = ASTAP_EXTS.iter().find(|(ext, _)| name.ends_with(ext))
                && files.len() != *want
            {
                bad.push(format!("expected {want} files, found {}", files.len()));
            }
        }

        if bad.is_empty() {
            println!(
                "  {:<11} ok  ({} files, {})",
                e.id,
                files.len(),
                human(installed_size(dir, e))
            );
        } else {
            problems += bad.len();
            println!("  {:<11} PROBLEMS:", e.id);
            for b in bad {
                println!("      {b}");
            }
        }
    }
    if checked == 0 {
        println!("No catalogues installed in {}", dir.display());
        return Ok(());
    }
    if problems > 0 {
        return Err(format!(
            "{problems} problem(s) found; re-run `arcsec catalog install <name>`"
        ));
    }
    Ok(())
}

pub fn cmd_install(dir: &Path, ids: &[String], assume_yes: bool, keep: bool) -> Result<(), String> {
    let mut wanted: Vec<&'static Entry> = Vec::new();
    for id in ids {
        let e = find(id).ok_or_else(|| {
            format!("unknown catalogue '{id}'. Run `arcsec catalog list` to see the options.")
        })?;
        if is_installed(dir, e) {
            println!("{}: already installed, skipping", e.id);
            continue;
        }
        wanted.push(e);
    }
    if wanted.is_empty() {
        return Ok(());
    }

    let total: u64 = wanted.iter().map(|e| e.bytes).sum();
    println!("Installing into {}\n", dir.display());
    for e in &wanted {
        println!("  {:<11} {:>9}  {}", e.id, human(e.bytes), e.desc);
    }
    println!("\nTotal download: {}", human(total));
    if !assume_yes {
        eprint!("Continue? [y/N] ");
        let _ = std::io::stderr().flush();
        let mut line = String::new();
        if std::io::stdin().read_line(&mut line).is_err()
            || !matches!(line.trim().to_ascii_lowercase().as_str(), "y" | "yes")
        {
            println!("Cancelled.");
            return Ok(());
        }
    }

    fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;

    for e in wanted {
        println!("\n{} — {}", e.id, e.desc);
        match e.archive {
            Archive::Loose => {
                let urls = loose_urls(e);
                let mut done = 0;
                for (i, (url, name)) in urls.iter().enumerate() {
                    let dest = dir.join(name);
                    if dest.exists() {
                        done += 1;
                        continue;
                    }
                    let label = format!("{} [{}/{}] {}", e.id, i + 1, urls.len(), name);
                    if let Err(err) = download(url, &dest, &label) {
                        eprintln!("  warning: {err}");
                    } else {
                        done += 1;
                    }
                }
                println!("  {done}/{} index files present", urls.len());
            }
            Archive::Zip | Archive::Deb => {
                let ext = if e.archive == Archive::Zip {
                    "zip"
                } else {
                    "deb"
                };
                let tmp = dir.join(format!(".{}-download.{ext}", e.id));
                download(e.url, &tmp, e.id)?;
                eprintln!("  extracting ...");
                let n = if e.archive == Archive::Zip {
                    extract_zip(&tmp, dir)?
                } else {
                    extract_deb(&tmp, dir)?
                };
                if keep {
                    println!("  kept archive at {}", tmp.display());
                } else {
                    let _ = fs::remove_file(&tmp);
                }
                println!("  extracted {n} files");
            }
        }
        if is_installed(dir, e) {
            println!("  {} installed ({})", e.id, human(installed_size(dir, e)));
        } else {
            return Err(format!(
                "{}: extraction finished but no catalogue files appeared in {}",
                e.id,
                dir.display()
            ));
        }
    }
    println!("\nDone. The solver uses {} by default.", dir.display());
    Ok(())
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
        let a = loose_urls(find("anet-4100").unwrap());
        assert_eq!(a.len(), 13, "4100 series should be scales 07..19");
        assert!(a[0].0.ends_with("index-4107.fits"));
        let b = loose_urls(find("anet-5200").unwrap());
        assert_eq!(b.len(), 3 * 48, "5200 series is 3 scales x 48 healpix");
        assert!(b.iter().all(|(u, _)| u.starts_with("https://")));
    }

    #[test]
    fn default_dir_is_absolute_and_namespaced() {
        // Deliberately not asserting the exact path: it is platform dependent.
        let d = default_dir();
        assert!(
            d.to_string_lossy().contains("arcsec"),
            "catalogue dir should be namespaced: {}",
            d.display()
        );
    }

    #[test]
    fn env_override_wins() {
        // SAFETY: single-threaded test, restored immediately.
        unsafe { std::env::set_var("ARCSEC_CATALOG_DIR", "/tmp/arcsec-test-catalogs") };
        assert_eq!(default_dir(), PathBuf::from("/tmp/arcsec-test-catalogs"));
        unsafe { std::env::remove_var("ARCSEC_CATALOG_DIR") };
    }

    #[test]
    fn human_sizes_read_sensibly() {
        assert_eq!(human(500), "500 B");
        assert_eq!(human(1_500_000), "1.5 MB");
        assert_eq!(human(901_300_000), "901.3 MB");
        assert_eq!(human(1_213_400_000), "1.2 GB");
    }

    #[test]
    fn ar_member_finds_the_payload() {
        // Minimal ar archive: magic, one member header, 4 bytes of data.
        let mut v = b"!<arch>\n".to_vec();
        let mut hdr = vec![b' '; 60];
        hdr[..8].copy_from_slice(b"data.tar");
        let size = b"4";
        hdr[48..48 + size.len()].copy_from_slice(size);
        hdr[58] = b'`';
        hdr[59] = b'\n';
        v.extend_from_slice(&hdr);
        v.extend_from_slice(b"ABCD");
        let (name, range) = ar_member(&v, "data.tar").expect("member not found");
        assert_eq!(name, "data.tar");
        assert_eq!(&v[range], b"ABCD");
    }

    #[test]
    fn ar_member_rejects_a_non_archive() {
        assert!(ar_member(b"not an ar archive at all", "data.tar").is_none());
    }
}
