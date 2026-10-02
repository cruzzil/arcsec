// `alloc` is not in the extern prelude for a crate that links std, so it has to
// be declared before alloc:: paths can be written.
extern crate alloc;

mod catalog_cmd;
mod cli;
mod extract;
mod logger;

use core::f64::consts::PI;
use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::process;
use std::time::Instant;

use arcsec_core::ArcsecError;
use arcsec_core::auto::{Event, MIN_SOLVE_DIM, Plan, SolveRequest};
use arcsec_core::pipeline::{SearchSpeed, SolveMethod, format_dec, format_ra, format_radec};
use arcsec_core::types::{ImageBuffer, WcsSolution};
use arcsec_core::wcs::TanWcs;
use arcsec_io::{fits_io, image_io};
use clap::ArgMatches;

use crate::cli::VERSION;
use crate::extract::Extract2;

fn main() {
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
        .unwrap_or_else(arcsec_core::auto::default_db_path);
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

    let hint = pointing_hint(&matches, file);
    let (ra_hint_rad, dec_hint_rad) = hint.unwrap_or((0.0, 0.0));

    // ── Plan: pixel scale, FOV, binning, database ────────────────────────────
    // Priority: explicit --fov flag > header FOCALLEN/XPIXSZ > 1"/px fallback.
    //
    // `--fov` is the image *height*, as ASTAP defines it and as N.I.N.A. sends it
    // (`FoVH`). The plan's `params.fov` is the field along the longer side, which
    // is what database selection and the search window use; for a square image the
    // two are the same number. The database is resolved from that field size: the
    // D-series covers 0.15°–6°, G05 3°–20° and W08 20°–80°.
    let fov_hint = matches.get_one::<f64>("fov").copied().filter(|v| *v > 0.0);
    let quad_tol = arg::<f64>(&matches, "tolerance");
    let request = SolveRequest {
        hint,
        fov_height: fov_hint.map(f64::to_radians),
        pixel_scale: if fov_hint.is_some() {
            None
        } else {
            image_io::read_pixel_scale(file)
        },
        search_radius: arg::<f64>(&matches, "radius").to_radians(),
        downsample: matches
            .get_one::<u32>("downsample")
            .map(|&z| usize::try_from(z).unwrap_or(usize::MAX)),
        db_path: Some(db_path),
        db_name: db_abbrev,
        index: matches.get_one::<PathBuf>("index").cloned(),
        index_first: true,
        auto_index: true,
        hfd_min_arcsec: arg::<f64>(&matches, "hfd-min"),
        quad_tolerance: quad_tol,
        max_stars,
        method: match matches.get_one::<String>("method").map(String::as_str) {
            Some("tetra") => SolveMethod::Tetra,
            _ => SolveMethod::Quads,
        },
        speed,
        threads,
        sip: want_sip,
        cancel: None,
    };
    let plan = match Plan::new(&request, img.width, img.height) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("Solver error: {e}");
            unsolved.exit(1);
        }
    };
    let binning = plan.binning;
    let (binned_w, binned_h) = plan.binned_size();
    log::info!(
        "Using star database {} for a {:.2}° field",
        plan.params.db_name.to_uppercase(),
        plan.params.fov.to_degrees()
    );

    // ── Solve header (always printed to stdout, like ASTAP) ──────────────────
    println!("arcsec astrometric solver version {VERSION}");
    println!(
        "Search radius: {:.0} degrees, ",
        plan.params.search_radius.to_degrees()
    );
    // ASTAP separates RA and Dec with a comma here, but not on "Solution found".
    println!(
        "Start position: {}, {}",
        format_ra(ra_hint_rad),
        format_dec(dec_hint_rad)
    );
    println!("Image height: {:.2} degrees", plan.fov_height.to_degrees());
    println!("Binning: {binning}x{binning}");
    println!(
        "Image dimensions: {}x{}",
        binned_w * binning,
        binned_h * binning
    );
    println!("Quad tolerance: {quad_tol:.3}");
    println!("Minimum star size: {:.1}\"", plan.hfd_min_arcsec);
    println!(
        "Speed: {}",
        if speed == SearchSpeed::Slow {
            "slow"
        } else {
            "normal"
        }
    );

    if binned_w < MIN_SOLVE_DIM || binned_h < MIN_SOLVE_DIM {
        eprintln!(
            "Insufficient stars: the image is too small to solve ({binned_w}x{binned_h} pixels)"
        );
        unsolved.exit(2);
    }

    // ── Solve ────────────────────────────────────────────────────────────────
    // The plan runs arcsec's own blind index when --index names one or the search
    // is wide enough to want one, the Astrometry.net blind solver when --index
    // names those files, and the catalogue spiral search; see arcsec_core::auto.
    let t0 = Instant::now();
    let solved = plan.solve_with(&img, |event| {
        if let Event::IndexEstimate(ra, dec) = event {
            println!(
                "Index position estimate: RA={:.3}°, Dec={:.3}°",
                ra.to_degrees(),
                dec.to_degrees()
            );
        }
    });
    let wcs = match solved {
        Ok(s) => s.wcs,
        Err(e) => report_failure(&unsolved, e),
    };
    let elapsed_s = t0.elapsed().as_secs_f64();

    print_solution(
        &wcs,
        quad_tol,
        elapsed_s,
        (ra_hint_rad, dec_hint_rad),
        do_progress,
    );

    // ── Write output files ───────────────────────────────────────────────────
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
/// pointing comes from the image header; failing that, none (the search starts
/// at (0, 0)).
fn pointing_hint(matches: &ArgMatches, file: &Path) -> Option<(f64, f64)> {
    let cli_ra = matches.get_one::<f64>("ra").copied();
    let cli_spd = matches.get_one::<f64>("spd").copied();
    if cli_ra.is_some() || cli_spd.is_some() {
        let ra = cli_ra.map_or(0.0, |h| h * PI / 12.0);
        let dec = cli_spd.map_or(0.0, |spd| (spd - 90.0).clamp(-90.0, 90.0) * PI / 180.0);
        Some((ra, dec))
    } else {
        image_io::read_ra_dec(file)
            .map(|(ra_deg, dec_deg)| (ra_deg * PI / 180.0, dec_deg * PI / 180.0))
    }
}

/// Report a failed solve as ASTAP would and exit with its code.
fn report_failure(unsolved: &Unsolved, err: ArcsecError) -> ! {
    match err {
        ArcsecError::InsufficientStars { found, required } => {
            eprintln!("Insufficient stars: found {found}, required {required}");
            unsolved.exit(2);
        }
        ArcsecError::InsufficientQuads { .. } => {
            println!("No solution found.");
            unsolved.exit(1);
        }
        ArcsecError::OutsideSearchRadius { separation_deg } => {
            eprintln!(
                "The blind index places this field {separation_deg:.1}° from the start position, \
                 outside the search radius."
            );
            println!("No solution found.");
            unsolved.exit(1);
        }
        ArcsecError::IndexNotFound(p) => {
            eprintln!("No index files found at {}", p.display());
            unsolved.exit(32);
        }
        ArcsecError::CatalogNotFound(p) => {
            eprintln!("Star database not found: {}", p.display());
            unsolved.exit(32);
        }
        ArcsecError::CatalogIo(e) => {
            eprintln!("Star database read error: {e}");
            unsolved.exit(33);
        }
        e => {
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
}
