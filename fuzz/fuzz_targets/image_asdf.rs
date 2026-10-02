//! ASDF images through the `asdf-rs` crate, as `arcsec -f image.asdf` reads them.
#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let data = arcsec_fuzz::with_magic(data, b"#ASDF ");
    arcsec_fuzz::exercise_image(&data, "img.asdf", false);
});
