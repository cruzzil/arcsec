//! FITS images through CFITSIO (rsfitsio), as `arcsec -f image.fits` reads them,
//! including gzip/bzip2/compress-wrapped and tile-compressed files, and the
//! `--update` header write.
#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    arcsec_fuzz::exercise_image(data, "img.fits", true);
});
