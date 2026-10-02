//! XISF images through the `xisf` crate, as `arcsec -f image.xisf` reads them.
#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let data = arcsec_fuzz::with_magic(data, b"XISF0100");
    arcsec_fuzz::exercise_image(&data, "img.xisf", false);
});
