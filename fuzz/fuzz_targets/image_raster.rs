//! JPEG, PNG and TIFF images through `zune-jpeg`, `png` and `tiff`, as
//! `arcsec -f image.jpg` reads them. The magic bytes are left to the input: the
//! seeds carry all three, and format detection picks the decoder.
#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    arcsec_fuzz::exercise_image(data, "img.raster", false);
});
