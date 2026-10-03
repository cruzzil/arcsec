//! Downloading and unpacking catalogues.

use std::fs;
use std::io::{self, BufReader, BufWriter, IsTerminal as _, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use super::human;

/// Is stderr a terminal? Controls whether progress redraws in place.
fn stderr_is_terminal() -> bool {
    io::stderr().is_terminal()
}

/// `name` with `suffix` appended: `a.zip` + `.part` = `a.zip.part`.
fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}

/// Where the resumable partial download of `dest` is kept.
fn part_path(dest: &Path) -> PathBuf {
    with_suffix(dest, ".part")
}

/// The complete length a 416 response reports (`Content-Range: bytes */N`).
fn unsatisfied_range_length(content_range: Option<&str>) -> Option<u64> {
    content_range?
        .trim()
        .strip_prefix("bytes */")?
        .trim()
        .parse()
        .ok()
}

/// Does a 206 response's `Content-Range` start where our partial file ends?
fn range_starts_at(content_range: Option<&str>, have: u64) -> bool {
    content_range
        .and_then(|v| v.trim().strip_prefix("bytes "))
        .and_then(|v| v.split('-').next())
        .and_then(|v| v.trim().parse::<u64>().ok())
        == Some(have)
}

/// Stream `url` to `dest`, reporting progress.
///
/// Downloads to `dest.part` and renames on success, so an interrupted run never
/// leaves a half file that later looks installed. Existing `.part` files are resumed
/// with a Range request where the server allows it — these are multi-hundred-megabyte
/// downloads and starting over is not acceptable.
pub fn download(url: &str, dest: &Path, label: &str) -> Result<(), String> {
    let part = part_path(dest);
    let have = fs::metadata(&part).map_or(0, |m| m.len());

    // Statuses are checked here rather than turned into errors by ureq, so that a
    // 416 on a resume can be told apart from a real failure.
    let agent = ureq::Agent::config_builder()
        .timeout_global(Some(core::time::Duration::from_secs(7200)))
        .http_status_as_error(false)
        .build()
        .new_agent();

    let mut req = agent.get(url);
    if have > 0 {
        req = req.header("Range", &format!("bytes={have}-"));
    }
    let resp = req.call().map_err(|e| format!("{label}: {e}"))?;
    let status = resp.status().as_u16();

    // 416: the partial file is not a prefix the server can continue. If the server
    // says the file is exactly as long as what we have, the download finished but
    // was never renamed (killed between the sync and the rename), so keep it.
    // Otherwise the file changed upstream or the part is corrupt: start over.
    if status == 416 && have > 0 {
        let range = resp
            .headers()
            .get("content-range")
            .and_then(|v| v.to_str().ok());
        if unsatisfied_range_length(range) == Some(have) {
            fs::rename(&part, dest).map_err(|e| format!("{}: {e}", dest.display()))?;
            eprintln!("  {label}: already downloaded ({})", human(have));
            return Ok(());
        }
        fs::remove_file(&part).map_err(|e| format!("{}: {e}", part.display()))?;
        return download(url, dest, label);
    }

    let resuming = status == 206;
    if !(status == 200 || resuming) {
        return Err(format!("{label}: HTTP {status}"));
    }
    if resuming {
        let range = resp
            .headers()
            .get("content-range")
            .and_then(|v| v.to_str().ok());
        if !range_starts_at(range, have) {
            let _ = fs::remove_file(&part);
            return Err(format!(
                "{label}: the server resumed at the wrong offset ({}); the partial \
                 download was discarded, run the install again",
                range.unwrap_or("no Content-Range")
            ));
        }
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

    let tty = stderr_is_terminal();
    let (cr, tail) = if tty { ("\r", "   ") } else { ("", "\n") };
    // Animate only on a terminal. Redirected to a file or a pipe, carriage returns
    // do not overwrite, so an unconditional update turns a 117 MB download into
    // hundreds of lines of noise.
    let interval_ms = if tty { 300 } else { 15_000 };

    let mut written = if resuming { have } else { 0 };
    let mut reader = resp.into_body().into_reader();
    let mut buf = vec![0u8; 1 << 20];
    let mut last_report = Instant::now();
    loop {
        let n = reader.read(&mut buf).map_err(|e| format!("{label}: {e}"))?;
        if n == 0 {
            break;
        }
        out.write_all(&buf[..n])
            .map_err(|e| format!("{}: {e}", part.display()))?;
        written += n as u64;
        if last_report.elapsed().as_millis() > interval_ms {
            match total {
                Some(t) if t > 0 => eprint!(
                    "{cr}  {label}: {} / {} ({:.0}%){tail}",
                    human(written),
                    human(t),
                    written as f64 / t as f64 * 100.0,
                ),
                _ => eprint!("{cr}  {label}: {}{tail}", human(written)),
            }
            let _ = io::stderr().flush();
            last_report = Instant::now();
        }
    }
    out.sync_all()
        .map_err(|e| format!("{}: {e}", part.display()))?;
    drop(out);
    eprintln!("{cr}  {label}: {} downloaded          ", human(written));

    // A connection that closes early can look like a clean end of stream. Keep the
    // partial file so the next run resumes it, but do not promote it.
    if let Some(t) = total
        && written != t
    {
        return Err(format!(
            "{label}: download ended at {} of {}; run the install again to resume",
            human(written),
            human(t)
        ));
    }
    fs::rename(&part, dest).map_err(|e| format!("{}: {e}", dest.display()))?;
    Ok(())
}

// ── Extraction ──────────────────────────────────────────────────────────────────

/// Copy one archive member to `dest/<name>`, flattening its path.
///
/// Only the final component of `path` is used, so a member such as `../../x` or
/// `/etc/x` still lands inside `dest` ("zip slip"); a path with no normal final
/// component (`..`, `/`) is skipped. Members `wanted` rejects are skipped too, so
/// packaging files (a `.deb`'s `copyright`, a zip's readme) never reach the
/// catalogue directory the solver reads.
///
/// The member is written beside its final name and renamed into place. The solver
/// memory-maps star-database files, and a file truncated under a mapping faults the
/// process reading it (SIGBUS), so a reinstall must never rewrite a tile in place
/// while a solve may be reading it; nor may an interrupted one leave a short tile
/// that reads as corrupt. Renaming also replaces a symbolic link at the name rather
/// than writing through it.
fn extract_member(
    mut reader: impl Read,
    path: &Path,
    dest: &Path,
    wanted: &dyn Fn(&str) -> bool,
) -> Result<bool, String> {
    let Some(name) = path.file_name() else {
        return Ok(false);
    };
    if !wanted(&name.to_string_lossy()) {
        return Ok(false);
    }
    let out_path = dest.join(name);
    let part = part_path(&out_path);
    let written = (|| {
        let mut out = BufWriter::new(fs::File::create(&part)?);
        io::copy(&mut reader, &mut out)?;
        out.flush()?;
        drop(out);
        fs::rename(&part, &out_path)
    })();
    if let Err(e) = written {
        let _ = fs::remove_file(&part);
        return Err(format!("{}: {e}", out_path.display()));
    }
    Ok(true)
}

/// Unpack the wanted files of a `.zip` into `dest`, flattening paths.
pub fn extract_zip(
    archive: &Path,
    dest: &Path,
    wanted: &dyn Fn(&str) -> bool,
) -> Result<usize, String> {
    let f = fs::File::open(archive).map_err(|e| format!("{}: {e}", archive.display()))?;
    let mut zip = zip::ZipArchive::new(BufReader::new(f))
        .map_err(|e| format!("{}: {e}", archive.display()))?;
    let mut n = 0;
    for i in 0..zip.len() {
        let item = zip.by_index(i).map_err(|e| e.to_string())?;
        if item.is_dir() {
            continue;
        }
        // `enclosed_name` refuses absolute paths and `..` escapes outright.
        let Some(path) = item.enclosed_name() else {
            continue;
        };
        if extract_member(item, &path, dest, wanted)? {
            n += 1;
        }
    }
    Ok(n)
}

/// Locate a member of an `ar` archive — the container format of a `.deb`.
///
/// Returns the member's name and its byte range. The format is trivial: an 8-byte
/// magic then, per member, a 60-byte ASCII header whose last fields are the size
/// and a `` `\n `` terminator, followed by the data padded to an even length.
/// Implementing it here avoids requiring `dpkg-deb` or `ar`, neither of which
/// exists on Windows.
fn ar_member<R: Read + Seek>(
    r: &mut R,
    want_prefix: &str,
) -> io::Result<Option<(String, u64, u64)>> {
    let len = r.seek(SeekFrom::End(0))?;
    r.seek(SeekFrom::Start(0))?;
    let mut magic = [0u8; 8];
    if len < 8 {
        return Ok(None);
    }
    r.read_exact(&mut magic)?;
    if &magic != b"!<arch>\n" {
        return Ok(None);
    }
    let mut pos = 8u64;
    while pos + 60 <= len {
        r.seek(SeekFrom::Start(pos))?;
        let mut hdr = [0u8; 60];
        r.read_exact(&mut hdr)?;
        // GNU ar terminates names with '/'; dpkg pads them with spaces.
        let name = String::from_utf8_lossy(&hdr[0..16])
            .trim()
            .trim_end_matches('/')
            .to_string();
        let Ok(size) = String::from_utf8_lossy(&hdr[48..58]).trim().parse::<u64>() else {
            return Ok(None);
        };
        let start = pos + 60;
        let Some(end) = start.checked_add(size).filter(|&e| e <= len) else {
            return Ok(None);
        };
        if name.starts_with(want_prefix) {
            return Ok(Some((name, start, size)));
        }
        pos = end + (size & 1); // members are padded to an even offset
    }
    Ok(None)
}

/// Unpack a `.deb`: find `data.tar.xz`, xz-decompress it, then untar the wanted
/// files, flattening paths (the ASTAP packages install under `opt/astap/`).
///
/// Streams throughout: the payload is decompressed to a temporary file beside the
/// archive rather than into memory, because D80's is over a gigabyte and a
/// Raspberry Pi is a common place to run this.
pub fn extract_deb(
    archive: &Path,
    dest: &Path,
    wanted: &dyn Fn(&str) -> bool,
) -> Result<usize, String> {
    let io_err = |e: io::Error| format!("{}: {e}", archive.display());
    let mut f = fs::File::open(archive).map_err(io_err)?;
    let (name, start, size) = ar_member(&mut f, "data.tar")
        .map_err(io_err)?
        .ok_or_else(|| format!("{}: no data.tar member (not a .deb?)", archive.display()))?;
    f.seek(SeekFrom::Start(start)).map_err(io_err)?;
    let payload = f.take(size);

    if name == "data.tar" {
        return untar(payload, dest, wanted);
    }
    if !name.ends_with(".xz") {
        return Err(format!(
            "{}: {name}: only xz-compressed .deb payloads are supported",
            archive.display()
        ));
    }

    eprintln!("  decompressing {name} ...");
    let tar_path = with_suffix(archive, ".tar");
    let result = (|| {
        let mut tar_out = BufWriter::new(
            fs::File::create(&tar_path).map_err(|e| format!("{}: {e}", tar_path.display()))?,
        );
        // lzma-rs overflows an addition on a corrupt stream footer (a backward
        // size of u32::MAX), which panics in a build with overflow checks; a
        // corrupt download is an error, not a crash.
        let mut input = BufReader::new(payload);
        std::panic::catch_unwind(core::panic::AssertUnwindSafe(|| {
            lzma_rs::xz_decompress(&mut input, &mut tar_out)
        }))
        .map_err(|_| {
            format!(
                "{}: xz decompress failed: corrupt stream",
                archive.display()
            )
        })?
        .map_err(|e| format!("{}: xz decompress failed: {e:?}", archive.display()))?;
        tar_out
            .flush()
            .map_err(|e| format!("{}: {e}", tar_path.display()))?;
        drop(tar_out);
        let tar_in =
            fs::File::open(&tar_path).map_err(|e| format!("{}: {e}", tar_path.display()))?;
        untar(BufReader::new(tar_in), dest, wanted)
    })();
    let _ = fs::remove_file(&tar_path);
    result
}

/// Extract the wanted regular files of a tar stream into `dest`, flattening paths.
fn untar(reader: impl Read, dest: &Path, wanted: &dyn Fn(&str) -> bool) -> Result<usize, String> {
    let mut ar = tar::Archive::new(reader);
    let mut n = 0;
    for entry in ar.entries().map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        // Regular files only: a symlink or hard link member could otherwise point
        // the next write outside `dest`.
        if !entry.header().entry_type().is_file() {
            continue;
        }
        let path = entry.path().map_err(|e| e.to_string())?.into_owned();
        if extract_member(entry, &path, dest, wanted)? {
            n += 1;
        }
    }
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn scratch(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("arcsec-fetch-{name}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        d
    }

    fn ar_with(members: &[(&str, &[u8])]) -> Vec<u8> {
        let mut v = b"!<arch>\n".to_vec();
        for (name, data) in members {
            let mut hdr = vec![b' '; 60];
            hdr[..name.len()].copy_from_slice(name.as_bytes());
            let size = data.len().to_string();
            hdr[48..48 + size.len()].copy_from_slice(size.as_bytes());
            hdr[58] = b'`';
            hdr[59] = b'\n';
            v.extend_from_slice(&hdr);
            v.extend_from_slice(data);
            if data.len() % 2 == 1 {
                v.push(b'\n');
            }
        }
        v
    }

    #[test]
    fn ar_member_finds_the_payload() {
        let v = ar_with(&[
            ("debian-binary", b"2.0\n"),
            ("control.tar.xz", b"abc"),
            ("data.tar", b"ABCD"),
        ]);
        let (name, start, size) = ar_member(&mut Cursor::new(&v), "data.tar")
            .unwrap()
            .expect("member not found");
        assert_eq!(name, "data.tar");
        let (s, n) = (
            usize::try_from(start).unwrap(),
            usize::try_from(size).unwrap(),
        );
        assert_eq!(&v[s..s + n], b"ABCD");
    }

    #[test]
    fn ar_member_strips_gnu_name_terminators() {
        let v = ar_with(&[("data.tar.xz/", b"XZ")]);
        let (name, _, _) = ar_member(&mut Cursor::new(&v), "data.tar")
            .unwrap()
            .unwrap();
        assert_eq!(name, "data.tar.xz");
    }

    #[test]
    fn ar_member_rejects_a_non_archive() {
        let mut c = Cursor::new(b"not an ar archive at all".to_vec());
        assert!(ar_member(&mut c, "data.tar").unwrap().is_none());
        assert!(
            ar_member(&mut Cursor::new(Vec::new()), "data.tar")
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn ar_member_rejects_a_truncated_member() {
        let mut v = ar_with(&[("data.tar", b"ABCD")]);
        v.truncate(v.len() - 2);
        assert!(
            ar_member(&mut Cursor::new(&v), "data.tar")
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn zip_extraction_flattens_and_cannot_escape() {
        let dir = scratch("zip");
        let archive = dir.join("a.zip");
        {
            let mut w = zip::ZipWriter::new(fs::File::create(&archive).unwrap());
            let o = zip::write::SimpleFileOptions::default();
            for name in ["d50/d50_0101.1476", "../../d50_0201.1476", "readme.txt"] {
                w.start_file(name, o).unwrap();
                w.write_all(b"data").unwrap();
            }
            w.finish().unwrap();
        }
        let out = dir.join("out");
        fs::create_dir_all(&out).unwrap();
        let n = extract_zip(&archive, &out, &|n| n.ends_with(".1476")).unwrap();
        assert!(
            out.join("d50_0101.1476").is_file(),
            "nested member flattened"
        );
        assert!(!out.join("readme.txt").exists(), "unwanted member skipped");
        assert!(
            !dir.join("d50_0201.1476").exists(),
            "escaped the destination"
        );
        assert!(n <= 2);
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn tar_extraction_flattens_and_cannot_escape() {
        let dir = scratch("tar");
        let mut b = tar::Builder::new(Vec::new());
        for name in ["opt/astap/v05_0101.290", "usr/share/doc/copyright"] {
            let mut h = tar::Header::new_gnu();
            h.set_size(4);
            h.set_mode(0o644);
            h.set_cksum();
            b.append_data(&mut h, name, &b"data"[..]).unwrap();
        }
        // An escaping name, written raw: the builder itself refuses `..`.
        let mut h = tar::Header::new_gnu();
        h.as_old_mut().name[..20].copy_from_slice(b"../../v05_0201.290\0\0");
        h.set_size(4);
        h.set_mode(0o644);
        h.set_cksum();
        b.append(&h, &b"data"[..]).unwrap();
        let bytes = b.into_inner().unwrap();

        let out = dir.join("out");
        fs::create_dir_all(&out).unwrap();
        let n = untar(Cursor::new(bytes), &out, &|n| n.ends_with(".290")).unwrap();
        assert!(out.join("v05_0101.290").is_file());
        assert!(!out.join("copyright").exists());
        assert!(
            !dir.join("v05_0201.290").exists(),
            "escaped the destination"
        );
        assert!(n <= 2);
        fs::remove_dir_all(&dir).ok();
    }

    /// Reinstalling replaces each file rather than rewriting it: a solve may have
    /// the old one memory-mapped, and truncating a mapped file faults the reader.
    /// A link at the name is replaced too, not written through.
    #[test]
    fn extraction_replaces_files_instead_of_rewriting_them() {
        let dir = scratch("replace");
        let out = dir.join("out");
        fs::create_dir_all(&out).unwrap();
        let tile = out.join("d50_0101.1476");
        fs::write(&tile, b"old tile").unwrap();
        let old = fs::File::open(&tile).unwrap();
        #[cfg(unix)]
        {
            fs::write(dir.join("victim"), b"untouched").unwrap();
            std::os::unix::fs::symlink(dir.join("victim"), out.join("d50_0201.1476")).unwrap();
        }
        let mut b = tar::Builder::new(Vec::new());
        for name in ["d50_0101.1476", "d50_0201.1476"] {
            let mut h = tar::Header::new_gnu();
            h.set_size(8);
            h.set_mode(0o644);
            h.set_cksum();
            b.append_data(&mut h, name, &b"new tile"[..]).unwrap();
        }
        let n = untar(Cursor::new(b.into_inner().unwrap()), &out, &|_| true).unwrap();
        assert_eq!(n, 2);
        assert_eq!(fs::read(&tile).unwrap(), b"new tile");
        // The file a reader already had open still holds what it held.
        let mut was = String::new();
        let mut old = old;
        old.read_to_string(&mut was).unwrap();
        assert_eq!(was, "old tile");
        #[cfg(unix)]
        {
            assert_eq!(fs::read(dir.join("victim")).unwrap(), b"untouched");
            assert!(
                fs::symlink_metadata(out.join("d50_0201.1476"))
                    .unwrap()
                    .is_file()
            );
        }
        let left: Vec<_> = fs::read_dir(&out)
            .unwrap()
            .map(|e| e.unwrap().file_name())
            .collect();
        assert_eq!(left.len(), 2, "temporary files left behind: {left:?}");
        fs::remove_dir_all(&dir).ok();
    }

    /// A `.deb` whose xz payload has a footer with a backward size of `u32::MAX`
    /// (found by fuzzing): lzma-rs computes `(backward_size + 1) << 2`, which
    /// overflows. The install must fail with an error, not a panic.
    #[test]
    fn a_corrupt_xz_payload_is_an_error() {
        let dir = scratch("xz");
        let xz: &[u8] = b"\xfd7zXZ\x00\x00\x04\xe6\xd6\xb4F\x00\x00\x00\x00\x1c\xdfD!x\x00\x00\
            \xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\x00\x00";
        let deb = dir.join("bad.deb");
        fs::write(&deb, ar_with(&[("data.tar.xz", xz)])).unwrap();
        let out = dir.join("out");
        fs::create_dir_all(&out).unwrap();
        let err = extract_deb(&deb, &out, &|_| true).unwrap_err();
        assert!(err.contains("xz decompress failed"), "{err}");
        assert!(
            !dir.join("bad.deb.tar").exists(),
            "temporary tar left behind"
        );
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn resumes_only_at_the_right_offset() {
        assert!(range_starts_at(Some("bytes 100-199/200"), 100));
        assert!(!range_starts_at(Some("bytes 0-199/200"), 100));
        assert!(!range_starts_at(None, 100));
        assert_eq!(unsatisfied_range_length(Some("bytes */200")), Some(200));
        assert_eq!(unsatisfied_range_length(Some("bytes 0-1/200")), None);
        assert_eq!(unsatisfied_range_length(None), None);
    }

    #[test]
    fn part_files_keep_the_whole_name() {
        assert_eq!(
            part_path(Path::new("d/.d50-download.zip")),
            Path::new("d/.d50-download.zip.part")
        );
        assert_eq!(
            part_path(Path::new("d/index-4107.fits")),
            Path::new("d/index-4107.fits.part")
        );
    }
}
