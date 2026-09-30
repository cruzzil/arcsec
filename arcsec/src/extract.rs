//! `--analyse`, `--extract` and `--extract2`: star statistics and star lists.
//!
//! The output mirrors `astap_cli` exactly, because scripts parse it:
//!
//! - stdout gets two lines, `HFD_MEDIAN=4.5` (one decimal) and `STARS=733`;
//! - the star list goes to the *image's* path with a `.csv` extension. `-o` does not
//!   move it, in ASTAP or here;
//! - the CSV header is always `x,y,hfd,snr,flux,ra[0..360],dec[0..360]`, even when
//!   the rows have no RA and Dec; rows are `x,y,hfd,snr,flux[,ra,dec]` with x and y
//!   as 1-based FITS pixels to four decimals, HFD to four, SNR and flux rounded to
//!   integers, RA and Dec in degrees to eight; lines end with the platform's line
//!   ending and the file with an empty line;
//! - with `--analyse` on Windows the exit code carries the result:
//!   `round(HFD × 100) × 1 000 000 + stars`. Elsewhere it is 0, since exit codes
//!   there only reach 255.

use std::path::{Path, PathBuf};

use arcsec_core::detection::{Analysis, MeasuredStar, analyse_image};
use arcsec_core::types::ImageBuffer;
use arcsec_core::wcs::TanWcs;

/// What ASTAP reports as the median HFD when there is none (no stars): the
/// largest value whose `--analyse` exit code still fits in a signed 32-bit int.
pub const NO_HFD: f64 = 21.47;

/// The minimum SNR an `--analyse`/`--extract`/`--extract2` value asks for; 0 means
/// ASTAP's default of 30.
pub fn snr_min(value: f64) -> f64 {
    if value == 0.0 { 30.0 } else { value }
}

/// Where the star list goes: the image path with its extension replaced by `.csv`
/// (`M31.fits` → `M31.csv`), whatever `-o` says, as ASTAP does.
pub fn csv_path(image: &Path) -> PathBuf {
    let mut s = image.with_extension("").into_os_string();
    s.push(".csv");
    PathBuf::from(s)
}

/// Analyse `img` and print ASTAP's two report lines. Returns the analysis and the
/// median HFD printed.
pub fn analyse_and_report(img: &ImageBuffer, snr_min: f64, max_stars: usize) -> (Analysis, f64) {
    let analysis = analyse_image(img, snr_min, max_stars);
    let hfd = analysis.hfd_median().unwrap_or(NO_HFD);
    println!("HFD_MEDIAN={hfd:.1}");
    println!("STARS={}", analysis.stars.len());
    (analysis, hfd)
}

/// The exit code `astap_cli -analyse` reports its result in on Windows.
///
/// ASTAP computes it in 64 bits and hands it to a 32-bit exit status, so an HFD
/// above 21.47 wraps; this wraps the same way.
pub fn analyse_exit_code(hfd_median: f64, stars: usize) -> i32 {
    let hfd = (hfd_median * 100.0).round_ties_even() as i64;
    (hfd * 1_000_000 + stars as i64) as i32
}

#[cfg(windows)]
const LINE_END: &str = "\r\n";
#[cfg(not(windows))]
const LINE_END: &str = "\n";

/// The CSV text for `stars`, with RA and Dec columns if `wcs` is given.
pub fn csv_text(stars: &[MeasuredStar], wcs: Option<&TanWcs>) -> String {
    use core::fmt::Write as _;
    let mut out = String::with_capacity(64 * (stars.len() + 1));
    out.push_str("x,y,hfd,snr,flux,ra[0..360],dec[0..360]");
    out.push_str(LINE_END);
    for s in stars {
        let (x, y) = (s.x + 1.0, s.y + 1.0);
        let _ = write!(
            out,
            "{x:.4},{y:.4},{:.4},{},{}",
            s.hfd,
            s.snr.round_ties_even() as i64,
            s.flux.round_ties_even() as i64
        );
        if let Some(w) = wcs {
            let (ra, dec) = w.pixel_to_sky(x, y);
            let _ = write!(out, ",{:.8},{:.8}", ra.to_degrees(), dec.to_degrees());
        }
        out.push_str(LINE_END);
    }
    // ASTAP writes the text with writeln, which ends the file with one more.
    out.push_str(LINE_END);
    out
}

/// Write the star list for `--extract`/`--extract2`.
pub fn write_csv(path: &Path, stars: &[MeasuredStar], wcs: Option<&TanWcs>) -> std::io::Result<()> {
    std::fs::write(path, csv_text(stars, wcs))
}

/// `--extract2`: analyse the image after the solve and write the star list, with
/// RA and Dec from the solution (or, if the solve failed, from a WCS already in
/// the header). The solve's own exit code stands, as in ASTAP.
pub struct Extract2 {
    /// Minimum SNR, after the 0 → 30 default.
    pub snr_min: f64,
    /// `-s`, which also bounds the analysis passes.
    pub max_stars: usize,
    /// Where the CSV goes.
    pub csv: PathBuf,
    /// The full-resolution image as solved (before any binning).
    pub img: ImageBuffer,
    /// A WCS already in the image header, for when the solve fails.
    pub header_wcs: Option<TanWcs>,
}

impl Extract2 {
    /// Analyse and write the CSV. Problems are warnings: the solve result, and
    /// its exit code, stand either way.
    pub fn run(&self, solution: Option<&TanWcs>) {
        let analysis = analyse_image(&self.img, self.snr_min, self.max_stars);
        let wcs = solution.or(self.header_wcs.as_ref());
        if let Err(e) = write_csv(&self.csv, &analysis.stars, wcs) {
            eprintln!("Warning: could not write {}: {e}", self.csv.display());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn star(x: f64, y: f64) -> MeasuredStar {
        MeasuredStar {
            x,
            y,
            hfd: 4.66784,
            snr: 274.5,
            flux: 377_475.6,
        }
    }

    #[test]
    fn csv_matches_astap_layout() {
        let text = csv_text(&[star(1588.67292, 14.80514)], None);
        let expected = format!(
            "x,y,hfd,snr,flux,ra[0..360],dec[0..360]{LINE_END}\
             1589.6729,15.8051,4.6678,274,377476{LINE_END}{LINE_END}"
        );
        assert_eq!(
            text, expected,
            "SNR 274.5 rounds half to even, as Pascal does"
        );
    }

    #[test]
    fn csv_adds_ra_and_dec_with_a_wcs() {
        let wcs = TanWcs {
            ra0: 210.8f64.to_radians(),
            dec0: 54.35f64.to_radians(),
            crpix1: 1589.6729,
            crpix2: 15.8051,
            cd: [[-3.45e-4, 0.0], [0.0, 3.45e-4]],
            sip: None,
        };
        let text = csv_text(&[star(1588.6729, 14.8051)], Some(&wcs));
        let row = text.lines().nth(1).unwrap();
        assert_eq!(
            row,
            "1589.6729,15.8051,4.6678,274,377476,210.80000000,54.35000000"
        );
        // An empty list is the header and the closing blank line.
        assert_eq!(
            csv_text(&[], Some(&wcs)),
            format!("x,y,hfd,snr,flux,ra[0..360],dec[0..360]{LINE_END}{LINE_END}")
        );
    }

    #[test]
    fn csv_goes_next_to_the_image() {
        assert_eq!(
            csv_path(Path::new("dir/M31.fits")),
            Path::new("dir/M31.csv")
        );
        assert_eq!(
            csv_path(Path::new("dir/30.00s_0018.fit")),
            Path::new("dir/30.00s_0018.csv")
        );
        assert_eq!(csv_path(Path::new("image")), Path::new("image.csv"));
    }

    #[test]
    fn windows_exit_code_packs_hfd_and_count() {
        assert_eq!(analyse_exit_code(4.5, 733), 450_000_733);
        assert_eq!(analyse_exit_code(NO_HFD, 0), 2_147_000_000);
        // Past 21.47 it wraps, as ASTAP's does.
        assert_eq!(analyse_exit_code(25.0, 7), (2_500_000_007_i64 as i32));
    }

    #[test]
    fn zero_snr_means_thirty() {
        assert_eq!(snr_min(0.0), 30.0);
        assert_eq!(snr_min(12.0), 12.0);
    }
}
