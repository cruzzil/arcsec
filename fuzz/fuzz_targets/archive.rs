//! Catalogue archive extraction (`arcsec catalog install`): zip, and `.deb` (an ar
//! archive holding `data.tar` or `data.tar.xz`). The first byte picks the form;
//! for the `.deb` forms other than raw, the rest is wrapped in an ar archive as
//! its `data.tar` member so the fuzzer reaches tar and xz directly.
//!
//! Besides not crashing, extraction must never write outside the destination
//! ("zip slip"), never create anything but regular files there, and never keep a
//! member the filter rejected.
#![no_main]

use std::fs;
use std::path::Path;

use arcsec_fuzz::fetch::{extract_deb, extract_zip};
use libfuzzer_sys::fuzz_target;

fn ar_with(name: &str, payload: &[u8]) -> Vec<u8> {
    let mut v = b"!<arch>\n".to_vec();
    let mut hdr = vec![b' '; 60];
    hdr[..name.len()].copy_from_slice(name.as_bytes());
    let size = payload.len().to_string();
    hdr[48..48 + size.len()].copy_from_slice(size.as_bytes());
    hdr[58] = b'`';
    hdr[59] = b'\n';
    v.extend_from_slice(&hdr);
    v.extend_from_slice(payload);
    if payload.len() % 2 == 1 {
        v.push(b'\n');
    }
    v
}

fn wanted(name: &str) -> bool {
    !name.starts_with("skip")
}

fn names(dir: &Path) -> Vec<String> {
    let mut v: Vec<String> = fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    v.sort();
    v
}

fuzz_target!(|data: &[u8]| {
    let Some((&which, rest)) = data.split_first() else {
        return;
    };
    let root = arcsec_fuzz::scratch_dir().join("archive");
    let _ = fs::remove_dir_all(&root);
    let out = root.join("out");
    fs::create_dir_all(&out).unwrap();
    let archive = root.join("a.bin");
    let result = match which % 4 {
        0 => {
            fs::write(&archive, rest).unwrap();
            extract_zip(&archive, &out, &wanted)
        }
        1 => {
            fs::write(&archive, rest).unwrap();
            extract_deb(&archive, &out, &wanted)
        }
        2 => {
            fs::write(&archive, ar_with("data.tar", rest)).unwrap();
            extract_deb(&archive, &out, &wanted)
        }
        _ => {
            fs::write(&archive, ar_with("data.tar.xz", rest)).unwrap();
            extract_deb(&archive, &out, &wanted)
        }
    };
    // Nothing may land beside the destination, and the decompressed payload is
    // cleaned up whether or not extraction succeeded.
    assert_eq!(
        names(&root),
        ["a.bin", "out"],
        "wrote outside the destination"
    );
    let got = names(&out);
    for name in &got {
        let meta = fs::symlink_metadata(out.join(name)).unwrap();
        assert!(meta.is_file(), "{name} is not a regular file");
        assert!(wanted(name), "{name} was filtered out but extracted");
    }
    if let Ok(n) = result {
        assert!(n >= got.len().min(1) || got.is_empty());
    }
});
