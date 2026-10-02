// `alloc` is not in the extern prelude for a crate that links std, so it has to
// be declared before alloc:: paths can be written.
extern crate alloc;

mod asdf_io;
mod blind;
mod catalog_cmd;
mod cli;
mod db_select;
mod extract;
mod fits_io;
mod image_io;
mod logger;
mod xisf_io;

use core::f64::consts::PI;
use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::process;
use std::time::Instant;

use arcsec_core::ArcsecError;
use arcsec_core::pipeline::{
    BlindSolveParams, SearchSpeed, SolveMethod, SolveParams, format_dec, format_ra, format_radec,
    solve_image,
};
use arcsec_core::types::{ImageBuffer, WcsSolution};
use arcsec_core::wcs::{TanWcs, fit_sip};
use clap::ArgMatches;

use crate::blind::BlindOutcome;
use crate::cli::VERSION;
use crate::extract::Extract2;

/// Smallest image side, in (binned) pixels, that reaches the solver. Detection
/// cannot run on a one-pixel-wide image, and would otherwise panic on it.
const MIN_SOLVE_DIM: usize = 2;

/// The unsolved `.ini` and command line, once known, for [`main`]'s panic handler.
static UNSOLVED_INI: std::sync::Mutex<Option<(PathBuf, String)>> = std::sync::Mutex::new(None);

/// Run arcsec, mapping a panic to a clean failure.
///
/// A panic is a bug, but whatever caused it, a program driving arcsec (N.I.N.A.,
/// Ekos, a script) should still see an ASTAP exit code and a `PLTSOLVD=F` `.ini`
/// rather than Rust's exit 101 and nothing. The panic hook has already printed the
/// message and location by the time it is caught here. Exit 1, "no solution", is
/// the generic failure; the image readers catch their own panics first and report
/// exit 16, an unreadable file, instead. Worker threads re-raise their panics on
/// the main thread (see `search_in_order`), so those are caught here too.
fn main() {
    if let Err(panic) = std::panic::catch_unwind(run) {
        eprintln!(
            "Error: internal error ({}). This is a bug in arcsec; please report it.",
            image_io::panic_message(panic.as_ref())
        );
        if let Some((ini, cmdline)) = UNSOLVED_INI.lock().ok().and_then(|mut u| u.take()) {
            let _ = fits_io::write_unsolved_ini_file(&ini, &cmdline);
        }
        process::exit(1);
    }
}

fn run() {
    // `arcsec catalog ...` is handled by its own parser. Dispatching on argv[1]
    // before clap sees it keeps the ASTAP-compatible flag form (`arcsec -f x.fits`)
    // completely untouched — no subcommand can shadow a flag, and no flag parsing
    // changes shape because a subcommand exists.
    if std::env::args_os().nth(1).is_some_and(|a| a == "catalog") {
        process::exit(catalog_cmd::run(std::env::args_os().skip(1)));
    }

    let argv: Vec<OsString> = std::env::args_os().collect();
    // clap exits 2 on a usage error, which is ASTAP's "insufficient stars"; a tool
    // reading the exit code would misreport it, so use 1 (no solution) instead.
    // --help and --version are not errors and keep clap's exit 0.
    let matches = cli::solver_command()
        .try_get_matches_from(cli::normalize_astap_args(argv.clone()))
        .unwrap_or_else(|e| {
            if e.use_stderr() {
                let _ = e.print();
                process::exit(1);
            }
            e.exit()
        });

    // ── File path ────────────────────────────────────────────────────────────
    let Some(file) = matches.get_one::<PathBuf>("file") else {
        eprintln!("Error: -f <filename> is required");
        process::exit(16);
    };

    // Apply the worker-thread limit before anything touches pixel data: the
    // normalise pass runs during the FITS read, well before the solve. Reaches
    // detection bands, the background histogram, the pixel-range scan, the spiral,
    // and the blind index passes.
    let threads = arg::<usize>(&matches, "threads");
    arcsec_core::set_max_threads(threads);

    let do_progress = matches.get_flag("progress");
    let out_base = output_base(file, matches.get_one::<PathBuf>("output"));
    let mut unsolved = Unsolved {
        ini_path: with_extension(&out_base, "ini"),
        cmdline: argv
            .iter()
            .map(|a| a.to_string_lossy())
            .collect::<Vec<_>>()
            .join(" "),
        extract2: None,
    };
    if let Ok(mut u) = UNSOLVED_INI.lock() {
        *u = Some((unsolved.ini_path.clone(), unsolved.cmdline.clone()));
    }

    let log_path = matches
        .get_flag("log")
        .then(|| with_extension(&out_base, "log"));
    logger::install(do_progress, log_path.as_deref());
    if do_progress || log_path.is_some() {
        let line: Vec<String> = argv
            .iter()
            .map(|a| a.to_string_lossy().into_owned())
            .collect();
        log::info!("{}", line.join(" "));
    }

    // ── Database ─────────────────────────────────────────────────────────────
    let db_path: PathBuf = matches
        .get_one::<PathBuf>("database")
        .cloned()
        .unwrap_or_else(db_select::default_db_path);
    let db_abbrev: Option<String> = matches.get_one::<String>("db-abbrev").cloned();
    log::info!("Creating grayscale image for solving");

    // ── Read image (FITS, XISF or ASDF) ──────────────────────────────────────
    let mut img = match image_io::read_image(file) {
        Ok(i) => i,
        Err(e) => {
            eprintln!("Error reading image: {e}");
            process::exit(16);
        }
    };

    // Replace non-finite pixels and, for float data in physical units (survey
    // cutouts in nanomaggies, Jy/beam, ...), rescale into the ADU-like range the
    // histogram background estimator needs. 16-bit camera data is left untouched.
    if let Some((scale, offset)) = img.normalize_for_detection() {
        log::info!(
            "Rescaled pixel data for detection: value * {scale:.6} + {offset:.3}              (float input with a range too small for the 16-bit histogram)"
        );
    }

    let max_stars = arg::<usize>(&matches, "stars");

    // ── Analyse only (--analyse, --extract) ──────────────────────────────────
    // No solve, no .ini or .wcs: just the report, and for --extract the star list.
    let analyse = matches.get_one::<f64>("analyse").copied();
    let extract = matches.get_one::<f64>("extract").copied();
    if analyse.is_some() || extract.is_some() {
        run_analysis(file, &img, analyse, extract, max_stars);
    }

    // ── Check-pattern filter (--check) ───────────────────────────────────────
    if matches.get_one::<String>("check").is_some_and(|v| v == "y") {
        if image_io::read_channels(file) > 1 {
            log::info!("Skipping check pattern filter. This filter works only for raw OSC images!");
        } else if img.check_pattern_filter() {
            log::info!("Applying check pattern filter.");
        }
    }

    // --extract2 analyses the full-resolution image after the solve, whether or
    // not it succeeds, so keep a copy before binning.
    unsolved.extract2 = matches.get_one::<f64>("extract2").map(|&v| Extract2 {
        snr_min: extract::snr_min(v),
        max_stars,
        csv: extract::csv_path(file),
        img: img.clone(),
        header_wcs: image_io::read_header_wcs(file),
    });
    // SIP is on with --sip (but not --sip n), and always for --extract2, whose
    // RA/Dec columns ASTAP computes through it.
    let want_sip =
        matches.get_one::<String>("sip").is_some_and(|v| v != "n") || unsolved.extract2.is_some();
    let speed = match matches.get_one::<String>("speed").map(String::as_str) {
        Some("slow") => SearchSpeed::Slow,
        _ => SearchSpeed::Auto,
    };

    let (ra_hint_rad, dec_hint_rad) = pointing_hint(&matches, file);

    // ── Pixel scale and FOV ──────────────────────────────────────────────────
    // Priority: explicit --fov flag > header FOCALLEN/XPIXSZ > 1"/px fallback.
    //
    // `--fov` is the image *height*, as ASTAP defines it and as N.I.N.A. sends it
    // (`FoVH`). `fov_rad` below is the field along the longer side, which is what
    // database selection and the search window use; for a square image the two are
    // the same number.
    let fov_hint = matches.get_one::<f64>("fov").copied().unwrap_or(0.0);
    let naxis = img.width.max(img.height) as f64;
    let height = img.height as f64;
    let mut scale_known = true;
    let (arcsec_per_px, fov_rad) = if fov_hint > 0.0 {
        let fov_height = fov_hint * PI / 180.0;
        let ps = fov_height.to_degrees() * 3600.0 / height;
        (ps, fov_height * (naxis / height))
    } else {
        let ps = image_io::read_pixel_scale(file).unwrap_or_else(|| {
            scale_known = false;
            1.0
        });
        let fov = naxis * ps / 3600.0 * PI / 180.0;
        (ps, fov)
    };

    let (image_w, image_h) = (img.width, img.height);
    let binning = choose_binning(
        matches.get_one::<u32>("downsample").copied(),
        arcsec_per_px,
        img.width,
        img.height,
    );
    let img = if binning > 1 {
        log::info!("Creating grayscale x {binning} binning image for solving/star alignment.");
        img.bin_image(binning)
    } else {
        img
    };

    // ── Star database ────────────────────────────────────────────────────────
    // Resolved here rather than earlier because the right database depends on the
    // field size: the D-series covers 0.15°–6°, G05 3°–20° and W08 20°–80°.
    let db_name: String = db_abbrev.unwrap_or_else(|| {
        db_select::select_db_for_fov(&db_path, fov_rad.to_degrees())
            .unwrap_or_else(|| "d80".to_string())
    });
    log::info!(
        "Using star database {} for a {:.2}° field",
        db_name.to_uppercase(),
        fov_rad.to_degrees()
    );

    let hfd_min_arcsec = arg::<f64>(&matches, "hfd-min");
    let hfd_min = (hfd_min_arcsec / (binning as f64 * arcsec_per_px)).max(0.8);
    // ASTAP treats a negative (or NaN) radius as zero and solves at the start
    // position; the library rejects it, so clamp here. `max` maps NaN to 0.
    let search_radius_rad = arg::<f64>(&matches, "radius").max(0.0) * PI / 180.0;
    let quad_tol = arg::<f64>(&matches, "tolerance");

    // ── Solve header (always printed to stdout, like ASTAP) ──────────────────
    println!("arcsec astrometric solver version {VERSION}");
    println!(
        "Search radius: {:.0} degrees, ",
        search_radius_rad.to_degrees()
    );
    // ASTAP separates RA and Dec with a comma here, but not on "Solution found".
    println!(
        "Start position: {}, {}",
        format_ra(ra_hint_rad),
        format_dec(dec_hint_rad)
    );
    println!(
        "Image height: {:.2} degrees",
        (fov_rad * (height / naxis)).to_degrees()
    );
    println!("Binning: {binning}x{binning}");
    println!(
        "Image dimensions: {}x{}",
        img.width * binning,
        img.height * binning
    );
    println!("Quad tolerance: {quad_tol:.3}");
    println!("Minimum star size: {hfd_min_arcsec:.1}\"");
    println!(
        "Speed: {}",
        if speed == SearchSpeed::Slow {
            "slow"
        } else {
            "normal"
        }
    );

    if img.width < MIN_SOLVE_DIM || img.height < MIN_SOLVE_DIM {
        eprintln!(
            "Insufficient stars: the image is too small to solve ({}x{} pixels)",
            img.width, img.height
        );
        unsolved.exit(2);
    }

    // ── Solve ────────────────────────────────────────────────────────────────
    let t0 = Instant::now();

    let method = match matches.get_one::<String>("method").map(String::as_str) {
        Some("tetra") => SolveMethod::Tetra,
        _ => SolveMethod::Quads,
    };

    let template = SolveParams {
        ra_hint: ra_hint_rad,
        dec_hint: dec_hint_rad,
        fov: fov_rad,
        search_radius: search_radius_rad,
        quad_tolerance: quad_tol,
        hfd_min,
        max_stars,
        db_path,
        db_name,
        binning,
        method,
        speed,
        threads,
    };

    // arcsec's own blind index, named by --index or, for a search wider than a
    // few fields, found in the catalogue directory: it finds the field and the
    // hinted solver accepts it (see blind::index_stage). With Astrometry.net
    // files, the blind solver estimates the position first, and that estimate
    // becomes the hint for the catalogue spiral solver.
    let has_hint = matches.get_one::<f64>("ra").is_some()
        || matches.get_one::<f64>("spd").is_some()
        || image_io::read_ra_dec(file).is_some();
    let own_index = blind::arcsec_index_for(matches.get_one::<PathBuf>("index"), &template);
    let index_wcs = match own_index.as_ref().map(|ix| {
        blind::index_stage(
            &img,
            ix,
            &template,
            has_hint,
            arcsec_per_px * binning as f64,
            scale_known,
        )
    }) {
        Some(blind::IndexOutcome::Solved(w)) => Some(*w),
        Some(blind::IndexOutcome::Elsewhere(sep_deg)) => {
            eprintln!(
                "The blind index places this field {sep_deg:.1}° from the start position, \
                 outside the search radius."
            );
            println!("No solution found.");
            unsolved.exit(1);
        }
        Some(blind::IndexOutcome::NotFound) | None => None,
    };

    let (ra, dec, search_radius) = match matches.get_one::<PathBuf>("index") {
        None => (ra_hint_rad, dec_hint_rad, search_radius_rad),
        _ if index_wcs.is_some() => (ra_hint_rad, dec_hint_rad, search_radius_rad),
        Some(_) if own_index.is_some() => {
            if index_wcs.is_none() {
                eprintln!("Blind index found no verified position. Falling back to hint.");
            }
            (ra_hint_rad, dec_hint_rad, search_radius_rad)
        }
        Some(idx_root) => {
            let index_files = blind::collect_index_files(idx_root, fov_rad.to_degrees());
            if index_files.is_empty() {
                eprintln!("No index files found at {}", idx_root.display());
                unsolved.exit(32);
            }
            let params = BlindSolveParams {
                quad_tolerance: quad_tol,
                hfd_min,
                max_stars,
                binning,
                // The blind scale filter maps this through the image height.
                fov_deg: (fov_rad * (height / naxis)).to_degrees(),
            };
            match blind::estimate_position(&img, &index_files, &params) {
                BlindOutcome::Found(ra, dec) => {
                    println!(
                        "Index position estimate: RA={:.3}°, Dec={:.3}°",
                        ra.to_degrees(),
                        dec.to_degrees()
                    );
                    // Narrow the catalog search so the spiral checks step 0 (the
                    // blind position) and at most a few neighbours: the blind
                    // position is off by at most one image width, so 2× fov is a
                    // generous ceiling.
                    (ra, dec, (fov_rad * 2.0).max(5.0_f64.to_radians()))
                }
                BlindOutcome::InsufficientStars { found, required } => {
                    eprintln!("Insufficient stars: found {found}, required {required}");
                    unsolved.exit(2);
                }
                BlindOutcome::NotFound => {
                    eprintln!(
                        "Blind position estimate failed for all index files. Falling back to hint."
                    );
                    (ra_hint_rad, dec_hint_rad, search_radius_rad)
                }
            }
        }
    };

    let mut wcs = match index_wcs {
        Some(w) => w,
        None => run_catalog_solve(
            &unsolved,
            &img,
            &SolveParams {
                ra_hint: ra,
                dec_hint: dec,
                search_radius,
                ..template
            },
        ),
    };
    if want_sip {
        wcs.sip = fit_sip(&wcs, image_w, image_h);
    }
    let elapsed_s = t0.elapsed().as_secs_f64();

    print_solution(
        &wcs,
        quad_tol,
        elapsed_s,
        (ra_hint_rad, dec_hint_rad),
        do_progress,
    );

    // ── Write output files ───────────────────────────────────────────────────
    // Solved: a later panic (in --update or --extract2) must not mark it unsolved.
    if let Ok(mut u) = UNSOLVED_INI.lock() {
        u.take();
    }
    let wcs_path = with_extension(&out_base, "wcs");
    let ini_path = &unsolved.ini_path;
    if let Err(e) = fits_io::write_wcs_file(&wcs_path, &wcs) {
        eprintln!("Warning: could not write {}: {e}", wcs_path.display());
    }
    if let Err(e) = fits_io::write_ini_file(ini_path, &wcs, max_stars, &unsolved.cmdline) {
        eprintln!("Warning: could not write {}: {e}", ini_path.display());
    }

    // ── Update input FITS header in-place (--update) ─────────────────────────
    if matches.get_flag("update") {
        if let Err(e) = image_io::update_wcs(file, &wcs) {
            eprintln!("Warning: --update failed: {e}");
        } else {
            log::info!("WCS written to FITS header");
        }
    }

    if let Some(job) = &unsolved.extract2 {
        job.run(Some(&TanWcs::from(&wcs)));
    }

    process::exit(0);
}

/// `--analyse` / `--extract`: report the median HFD and star count, write the star
/// list for `--extract`, and exit without solving.
///
/// Given both, `--extract`'s minimum SNR is the one used, as in ASTAP. The exit code
/// is 0, except that on Windows `--analyse` reports its result in it (see
/// [`extract::analyse_exit_code`]).
fn run_analysis(
    file: &Path,
    img: &ImageBuffer,
    analyse: Option<f64>,
    extract: Option<f64>,
    max_stars: usize,
) -> ! {
    let snr_min = extract::snr_min(extract.or(analyse).unwrap_or(0.0));
    let (analysis, hfd) = extract::analyse_and_report(img, snr_min, max_stars);
    if extract.is_some() {
        let csv = extract::csv_path(file);
        let header_wcs = image_io::read_header_wcs(file);
        if let Err(e) = extract::write_csv(&csv, &analysis.stars, header_wcs.as_ref()) {
            eprintln!("Error: could not write {}: {e}", csv.display());
            process::exit(16);
        }
    }
    if cfg!(windows) && analyse.is_some() {
        process::exit(extract::analyse_exit_code(hfd, analysis.stars.len()));
    }
    process::exit(0);
}

/// A flag with a default value, which is therefore always present.
fn arg<T: Clone + Send + Sync + 'static>(matches: &ArgMatches, id: &str) -> T {
    matches
        .get_one::<T>(id)
        .cloned()
        .unwrap_or_else(|| unreachable!("--{id} has a default value"))
}

/// Base path for the output files: `-o` if given, else the image path without its
/// last extension (`30.00s_0018.fits` → `30.00s_0018`).
fn output_base(file: &Path, output: Option<&PathBuf>) -> PathBuf {
    output.map_or_else(|| file.with_extension(""), Clone::clone)
}

/// `base` + `.ext`, appended rather than substituted.
///
/// `Path::with_extension` would treat a dot inside the stem as the start of an
/// extension, turning `30.00s_0018` into `30.wcs`.
fn with_extension(base: &Path, ext: &str) -> PathBuf {
    let mut s = base.as_os_str().to_owned();
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

/// The search centre (RA, Dec) in radians.
///
/// `--ra` is in hours (0–24) and `--spd` is south pole distance in degrees
/// (0–180, SPD = 90 + Dec), both ASTAP conventions. If neither is given, the
/// pointing comes from the image header; failing that, (0, 0).
fn pointing_hint(matches: &ArgMatches, file: &Path) -> (f64, f64) {
    let cli_ra = matches.get_one::<f64>("ra").copied();
    let cli_spd = matches.get_one::<f64>("spd").copied();
    if cli_ra.is_some() || cli_spd.is_some() {
        let ra = cli_ra.map_or(0.0, |h| h * PI / 12.0);
        let dec = cli_spd.map_or(0.0, |spd| (spd - 90.0).clamp(-90.0, 90.0) * PI / 180.0);
        (ra, dec)
    } else if let Some((ra_deg, dec_deg)) = image_io::read_ra_dec(file) {
        (ra_deg * PI / 180.0, dec_deg * PI / 180.0)
    } else {
        (0.0, 0.0)
    }
}

/// The binning factor: `-z` if given and non-zero, else automatic.
///
/// Automatic binning brings a sampling finer than 1"/px back to about 1"/px, up to
/// 16×. Either way the factor is capped so the binned image keeps at least
/// [`MIN_SOLVE_DIM`] pixels a side: binning past the image size leaves nothing to
/// detect in, and used to crash.
fn choose_binning(
    requested: Option<u32>,
    arcsec_per_px: f64,
    width: usize,
    height: usize,
) -> usize {
    let binning = match requested {
        Some(0) | None => {
            if arcsec_per_px < 1.0 {
                (1.0 / arcsec_per_px).round().clamp(1.0, 16.0) as usize
            } else {
                1
            }
        }
        Some(z) => usize::try_from(z).unwrap_or(usize::MAX),
    };
    binning.min((width.min(height) / MIN_SOLVE_DIM).max(1))
}

/// Run the catalog-based (spiral search) solver, exiting with the ASTAP exit code
/// if it fails.
fn run_catalog_solve(unsolved: &Unsolved, img: &ImageBuffer, params: &SolveParams) -> WcsSolution {
    match solve_image(img, params) {
        Ok(w) => w,
        Err(ArcsecError::InsufficientStars { found, required }) => {
            eprintln!("Insufficient stars: found {found}, required {required}");
            unsolved.exit(2);
        }
        Err(ArcsecError::InsufficientQuads { .. }) => {
            println!("No solution found.");
            unsolved.exit(1);
        }
        Err(ArcsecError::CatalogNotFound(p)) => {
            eprintln!("Star database not found: {}", p.display());
            unsolved.exit(32);
        }
        Err(ArcsecError::CatalogIo(e)) => {
            eprintln!("Star database read error: {e}");
            unsolved.exit(33);
        }
        Err(e) => {
            eprintln!("Solver error: {e}");
            unsolved.exit(1);
        }
    }
}

/// How to report a solve that ended without a solution.
///
/// ASTAP writes an `.ini` holding `PLTSOLVD=F` whenever a solve fails, and tools
/// that drive it (N.I.N.A., Ekos, ...) poll that file rather than the exit code, so
/// arcsec does the same before exiting.
///
/// `--extract2` writes its star list however the solve ends, so that job lives
/// here too.
struct Unsolved {
    ini_path: PathBuf,
    cmdline: String,
    extract2: Option<Extract2>,
}

impl Unsolved {
    fn exit(&self, code: i32) -> ! {
        if let Err(e) = fits_io::write_unsolved_ini_file(&self.ini_path, &self.cmdline) {
            eprintln!("Warning: could not write {}: {e}", self.ini_path.display());
        }
        if let Some(job) = &self.extract2 {
            job.run(None);
        }
        process::exit(code);
    }
}

/// Print the solution report to stdout, in ASTAP's format.
fn print_solution(
    wcs: &WcsSolution,
    quad_tol: f64,
    elapsed_s: f64,
    (ra_hint_rad, dec_hint_rad): (f64, f64),
    progress: bool,
) {
    // Step-distance progress line (only without --progress; ASTAP style).
    if !progress && !wcs.step_distances.is_empty() {
        let dots: Vec<String> = wcs
            .step_distances
            .iter()
            .map(|d| format!("{d:.0}d"))
            .collect();
        println!("{},", dots.join(","));
    }

    let p = &wcs.plate;
    let n_matched = wcs.stars_matched;
    let n_raw = wcs.raw_matches;
    println!("{n_matched} of {n_raw} quads selected matching within {quad_tol:.3} tolerance.");
    println!(
        "Solution[\"] x:={:.6}*x+ {:.6}*y+ {:.6},  y:={:.6}*x+ {:.6}*y+ {:.6}",
        p.a, p.b, p.c, p.d, p.e, p.f
    );

    let sol_str = format_radec(wcs.ra0, wcs.dec0);
    println!("Solution found: {sol_str}");
    log::info!("Solution found: {sol_str}");

    let delta_str = if wcs.search_dist_deg >= 1.0 {
        format!("{:.1}d", wcs.search_dist_deg)
    } else {
        format!("{:.1}\"", wcs.search_dist_deg * 3600.0)
    };
    let dra_arcsec = (wcs.ra0 - ra_hint_rad) * wcs.dec0.cos() * (180.0 / PI) * 3600.0;
    let ddec_arcsec = (wcs.dec0 - dec_hint_rad) * (180.0 / PI) * 3600.0;
    println!(
        "Solved in {elapsed_s:.1} sec. Δ was {delta_str}.  Mount Δα={dra_arcsec:.1}\",  Δδ={ddec_arcsec:.1}\".  Used stars down to magnitude: {:.1}",
        wcs.mag_limit,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_base_strips_only_the_last_extension() {
        let b = |f: &str| output_base(Path::new(f), None);
        assert_eq!(b("dir/30.00s_0018.fits"), Path::new("dir/30.00s_0018"));
        assert_eq!(b("light.fit"), Path::new("light"));
        // A dot in a directory name is not an extension.
        assert_eq!(
            b("/data/2026.09.25/image"),
            Path::new("/data/2026.09.25/image")
        );
        let o = PathBuf::from("out/solved");
        assert_eq!(output_base(Path::new("x.fits"), Some(&o)), o);
    }

    #[test]
    fn output_files_append_their_extension() {
        assert_eq!(
            with_extension(Path::new("dir/30.00s_0018"), "wcs"),
            Path::new("dir/30.00s_0018.wcs")
        );
    }

    #[cfg(unix)]
    #[test]
    fn non_utf8_image_paths_survive() {
        use std::os::unix::ffi::OsStrExt as _;
        let f = Path::new(std::ffi::OsStr::from_bytes(b"dir/caf\xe9.fits"));
        let base = output_base(f, None);
        assert_eq!(base.as_os_str().as_bytes(), b"dir/caf\xe9");
        assert_eq!(
            with_extension(&base, "ini").as_os_str().as_bytes(),
            b"dir/caf\xe9.ini"
        );
    }

    #[test]
    fn binning_is_automatic_below_one_arcsec_per_pixel() {
        assert_eq!(choose_binning(None, 2.0, 4000, 3000), 1);
        assert_eq!(choose_binning(Some(0), 0.5, 4000, 3000), 2);
        assert_eq!(choose_binning(None, 0.01, 4000, 3000), 16);
        assert_eq!(choose_binning(Some(3), 2.0, 4000, 3000), 3);
    }

    #[test]
    fn binning_never_exceeds_the_image() {
        assert_eq!(choose_binning(Some(100), 1.0, 4, 4), 2);
        assert_eq!(choose_binning(Some(u32::MAX), 1.0, 4, 4), 2);
        assert_eq!(choose_binning(None, 0.01, 10, 50), 5);
        assert_eq!(choose_binning(Some(4), 1.0, 1, 1), 1);
    }
}
