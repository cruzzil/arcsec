// `alloc` is not in the extern prelude for a crate that links std, so it has to
// be declared before alloc:: paths can be written.
extern crate alloc;

mod asdf_io;
mod catalog_cmd;
mod fits_io;
mod image_io;
mod xisf_io;

use alloc::sync::Arc;
use std::fs;
use std::io::Write as IoWrite;
use std::path::{Path, PathBuf};
use std::process;
use std::sync::Mutex;
use std::time::Instant;

use arcsec_core::catalog::{load_anet_index, peek_anet_scale};
use arcsec_core::pipeline::{
    BlindSolveParams, SolveMethod, SolveParams, blind_solve, format_radec, solve_image,
};
use clap::{Arg, Command, value_parser};
use env_logger::Logger as EnvLogger;
use log::{LevelFilter, Log, Metadata, Record};

const VERSION: &str = "0.1.0";

// ── Logger: env_logger for stderr, thin wrapper adds optional file sink ──────

struct ArcsecLogger {
    /// env_logger backend — handles format + stderr (Some when -progress).
    stderr: Option<EnvLogger>,
    /// Optional file sink for -log.
    file: Option<Arc<Mutex<fs::File>>>,
}

impl Log for ArcsecLogger {
    fn enabled(&self, m: &Metadata<'_>) -> bool {
        self.stderr
            .as_ref()
            .map_or(self.file.is_some(), |l| l.enabled(m))
    }

    fn log(&self, r: &Record<'_>) {
        if let Some(l) = &self.stderr {
            l.log(r);
        }
        if let Some(f) = &self.file {
            let line = format!("{}  {}", hms_now(), r.args());
            if let Ok(mut g) = f.lock() {
                let _ = writeln!(g, "{line}");
            }
        }
    }

    fn flush(&self) {
        if let Some(l) = &self.stderr {
            l.flush();
        }
    }
}

fn main() {
    // `arcsec catalog ...` is handled by its own parser. Dispatching on argv[1]
    // before clap sees it keeps the ASTAP-compatible flag form (`arcsec -f x.fits`)
    // completely untouched — no subcommand can shadow a flag, and no flag parsing
    // changes shape because a subcommand exists.
    if std::env::args().nth(1).as_deref() == Some("catalog") {
        std::process::exit(run_catalog_cli());
    }

    let matches = Command::new("arcsec")
        .version(VERSION)
        .about("Astrometric solver — ASTAP-compatible CLI")
        .arg(
            Arg::new("file")
                .short('f')
                .value_parser(value_parser!(PathBuf))
                .help("filename {fits, xisf, asdf files}"),
        )
        .arg(
            Arg::new("radius")
                .short('r')
                .value_parser(value_parser!(f64))
                .allow_hyphen_values(true)
                .default_value("180")
                .help("radius_area_to_search[degrees]"),
        )
        .arg(
            Arg::new("fov")
                .long("fov")
                .value_parser(value_parser!(f64))
                .allow_hyphen_values(true)
                .help("diameter_field[degrees] {enter zero for auto}"),
        )
        .arg(
            Arg::new("ra")
                .long("ra")
                .value_parser(value_parser!(f64))
                .allow_hyphen_values(true)
                .help("right_ascension[hours]"),
        )
        .arg(
            Arg::new("spd")
                .long("spd")
                .value_parser(value_parser!(f64))
                .allow_hyphen_values(true)
                .help("south_pole_distance[degrees]"),
        )
        .arg(
            Arg::new("stars")
                .short('s')
                .value_parser(value_parser!(usize))
                .default_value("500")
                .help("max_number_of_stars  {default 500}"),
        )
        .arg(
            Arg::new("tolerance")
                .short('t')
                .value_parser(value_parser!(f64))
                .default_value("0.007")
                .help("quad_tolerance  {default 0.007}"),
        )
        .arg(
            Arg::new("hfd-min")
                .short('m')
                .value_parser(value_parser!(f64))
                .default_value("1.5")
                .help("minimum_star_size[\"]  {default 1.5}"),
        )
        .arg(
            Arg::new("downsample")
                .short('z')
                .value_parser(value_parser!(u32))
                .help("downsample_factor[0,1,2,3,4,..] {0 for auto}"),
        )
        .arg(
            Arg::new("check")
                .long("check")
                .action(clap::ArgAction::SetTrue)
                .help("[not implemented] Apply check pattern filter prior to solving"),
        )
        .arg(
            Arg::new("database")
                .short('d')
                .value_parser(value_parser!(PathBuf))
                .help("path {star database directory; default: the arcsec catalogue directory}"),
        )
        .arg(
            Arg::new("db-abbrev")
                .short('D')
                .help("abbreviation[d80,d50,g05,w08,...] {Specify a star database}"),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .value_parser(value_parser!(PathBuf))
                .help("file {Name the output files with this base path & file name}"),
        )
        .arg(
            Arg::new("sip")
                .long("sip")
                .action(clap::ArgAction::SetTrue)
                .help("[not implemented] Add SIP (Simple Image Polynomial) coefficients"),
        )
        .arg(Arg::new("speed").long("speed").help("mode[auto/slow] {only auto is implemented}"))
        .arg(
            Arg::new("index")
                .long("index")
                .short('i')
                .value_parser(value_parser!(PathBuf))
                .help("astrometry.net index file or directory of index-*.fits files for blind solving"),
        )
        .arg(
            Arg::new("threads")
                .long("threads")
                .value_parser(value_parser!(usize))
                .default_value("0")
                .help("max worker threads for detection and the search {0 = one per core}"),
        )
        .arg(
            Arg::new("method")
                .long("method")
                .default_value("quads")
                .help("catalog matching method: quads (default) or tetra"),
        )
        .arg(
            Arg::new("wcs")
                .long("wcs")
                .action(clap::ArgAction::SetTrue)
                .help("Write a .wcs file in similar format as Astrometry.net {always written}"),
        )
        .arg(
            Arg::new("log")
                .long("log")
                .action(clap::ArgAction::SetTrue)
                .help("Write the solver log to a .log text file"),
        )
        .arg(
            Arg::new("update")
                .long("update")
                .action(clap::ArgAction::SetTrue)
                .help("Add the solution to the input fits file header"),
        )
        .arg(
            Arg::new("progress")
                .long("progress")
                .action(clap::ArgAction::SetTrue)
                .help("Log all progress steps and messages"),
        )
        .arg(
            Arg::new("analyse")
                .long("analyse")
                .value_parser(value_parser!(f64))
                .help("[not implemented] snr_min {Analyse only and report median HFD and number of stars}"),
        )
        .arg(
            Arg::new("extract")
                .long("extract")
                .value_parser(value_parser!(f64))
                .help("[not implemented] snr_min {As -analyse but export star info to .csv}"),
        )
        .arg(
            Arg::new("extract2")
                .long("extract2")
                .value_parser(value_parser!(f64))
                .help("[not implemented] snr_min {Solve and export star info including ra, dec to .csv}"),
        )
        .get_matches();

    // ── Unimplemented flags ──────────────────────────────────────────────────
    // These parse, for ASTAP command-line compatibility, but nothing reads them.
    // Accepting an option and then ignoring what it asked for is worse than
    // refusing it: --sip silently returned a solution with no SIP coefficients,
    // and --analyse ran a full solve and wrote output files.
    //
    // Not listed here, because they ask for what already happens: --wcs (the .wcs
    // file is always written) and --speed auto (the only mode there is).
    for (flag, what) in [
        ("check", "--check (check-pattern filter)"),
        ("sip", "--sip (SIP distortion coefficients)"),
        ("analyse", "--analyse (analyse-only mode)"),
        ("extract", "--extract (star list export)"),
        ("extract2", "--extract2 (solved star list export)"),
    ] {
        // analyse/extract/extract2 parse as f64, check/sip are bare flags.
        let given = match flag {
            "check" | "sip" => matches.get_flag(flag),
            _ => matches.get_one::<f64>(flag).is_some(),
        };
        if given {
            eprintln!("Error: {what} is not implemented");
            process::exit(1);
        }
    }
    if let Some(speed) = matches.get_one::<String>("speed")
        && speed != "auto"
    {
        eprintln!("Error: --speed {speed} is not implemented (only 'auto')");
        process::exit(1);
    }

    // ── File path ────────────────────────────────────────────────────────────
    let file: &PathBuf = matches.get_one("file").unwrap_or_else(|| {
        eprintln!("Error: -f <filename> is required");
        process::exit(16);
    });

    // Apply the worker-thread limit before anything touches pixel data: the
    // normalise pass runs during the FITS read, well before the solve. Reaches
    // detection bands, the background histogram, the pixel-range scan, the spiral,
    // and the blind index passes.
    let threads = *matches.get_one::<usize>("threads").unwrap();
    arcsec_core::set_max_threads(threads);

    let do_log = matches.get_flag("log");
    let do_progress = matches.get_flag("progress");

    // ── Output base path ─────────────────────────────────────────────────────
    // Strip the last extension from the image path as a string to avoid
    // PathBuf::with_extension treating dots in the stem (e.g. "30.00s_0018")
    // as an extension separator.
    let out_base: PathBuf = if let Some(o) = matches.get_one::<PathBuf>("output") {
        o.clone()
    } else {
        let p = file.to_string_lossy();
        let base = p.rfind('.').map_or(&*p, |i| &p[..i]);
        PathBuf::from(base)
    };

    // ── Open log file early ──────────────────────────────────────────────────
    let log_file: Option<Arc<Mutex<fs::File>>> = if do_log {
        let lp = PathBuf::from(format!("{}.log", out_base.display()));
        match fs::File::create(&lp) {
            Ok(f) => Some(Arc::new(Mutex::new(f))),
            Err(e) => {
                eprintln!("Warning: cannot create {}: {e}", lp.display());
                None
            }
        }
    } else {
        None
    };

    // ── Install logger ───────────────────────────────────────────────────────
    let stderr_logger: Option<EnvLogger> = if do_progress {
        let mut b = env_logger::Builder::new();
        b.filter_level(LevelFilter::Info);
        b.format(|buf, r| writeln!(buf, "{}  {}", hms_now(), r.args()));
        Some(b.build())
    } else {
        None
    };
    let logger = Box::new(ArcsecLogger {
        stderr: stderr_logger,
        file: log_file.clone(),
    });
    log::set_boxed_logger(logger).ok();
    log::set_max_level(if do_log || do_progress {
        LevelFilter::Info
    } else {
        LevelFilter::Off
    });

    // Log command invocation
    if do_log || do_progress {
        let argv: Vec<String> = std::env::args().collect();
        log::info!("{}", argv.join(" "));
    }

    // ── Database ─────────────────────────────────────────────────────────────
    // Without -d, use the directory `arcsec catalog install` writes to, so a user
    // who installed a catalogue never has to say where it went. Fall back to the
    // working directory (ASTAP's behaviour) if nothing is installed there.
    let db_path: PathBuf = matches
        .get_one::<PathBuf>("database")
        .cloned()
        .unwrap_or_else(|| {
            let managed = catalog_cmd::default_dir();
            if catalog_cmd::REGISTRY
                .iter()
                .any(|e| catalog_cmd::is_installed(&managed, e))
            {
                managed
            } else {
                PathBuf::from(".")
            }
        });
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
    let img = img;

    // ── RA / SPD / FOV ───────────────────────────────────────────────────────
    // -ra in hours (0-24), -spd in south-pole-distance degrees (0-180, SPD = 90 + dec)
    // If neither is given, fall back to RA/DEC keywords in the FITS header.
    let cli_ra = matches.get_one::<f64>("ra").copied();
    let cli_spd = matches.get_one::<f64>("spd").copied();
    let (ra_hint_rad, dec_hint_rad) = if cli_ra.is_some() || cli_spd.is_some() {
        let ra = cli_ra
            .map(|h| h * core::f64::consts::PI / 12.0)
            .unwrap_or(0.0);
        let dec = cli_spd
            .map(|spd| (spd - 90.0).clamp(-90.0, 90.0) * core::f64::consts::PI / 180.0)
            .unwrap_or(0.0);
        (ra, dec)
    } else if let Some((ra_deg, dec_deg)) = image_io::read_ra_dec(file) {
        (
            ra_deg * core::f64::consts::PI / 180.0,
            dec_deg * core::f64::consts::PI / 180.0,
        )
    } else {
        (0.0, 0.0)
    };

    // ── Pixel scale and FOV ──────────────────────────────────────────────────
    // Priority: explicit --fov flag > FITS header FOCALLEN/XPIXSZ > naxis/3600 fallback.
    let fov_hint = matches.get_one::<f64>("fov").copied().unwrap_or(0.0);
    let naxis = img.width.max(img.height) as f64;
    let (arcsec_per_px, fov_rad) = if fov_hint > 0.0 {
        let fov = fov_hint * core::f64::consts::PI / 180.0;
        let ps = fov.to_degrees() * 3600.0 / naxis;
        (ps, fov)
    } else if let Some(ps) = image_io::read_pixel_scale(file) {
        let fov = naxis * ps / 3600.0 * core::f64::consts::PI / 180.0;
        (ps, fov)
    } else {
        let ps = 1.0;
        let fov = naxis * ps / 3600.0 * core::f64::consts::PI / 180.0;
        (ps, fov)
    };

    // ── Auto-binning ─────────────────────────────────────────────────────────
    // Bin when height > 2500 px OR pixel scale < 1 arcsec/px.
    let downsample_arg = matches.get_one::<u32>("downsample").copied();
    let binning: usize = match downsample_arg {
        Some(0) | None => {
            if arcsec_per_px < 1.0 {
                (1.0 / arcsec_per_px).round().clamp(1.0, 16.0) as usize
            } else {
                1
            }
        }
        Some(z) => z as usize,
    };

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
        select_db_for_fov(&db_path, fov_rad.to_degrees()).unwrap_or_else(|| "d80".to_string())
    });
    log::info!(
        "Using star database {} for a {:.2}° field",
        db_name.to_uppercase(),
        fov_rad.to_degrees()
    );

    let hfd_min_arcsec = *matches.get_one::<f64>("hfd-min").unwrap();
    let hfd_min = (hfd_min_arcsec / (binning as f64 * arcsec_per_px)).max(0.8);

    let search_radius_rad =
        *matches.get_one::<f64>("radius").unwrap() * core::f64::consts::PI / 180.0;
    let quad_tol = *matches.get_one::<f64>("tolerance").unwrap();
    let max_stars = *matches.get_one::<usize>("stars").unwrap();

    // ── Solve header (always printed to stdout, like ASTAP) ──────────────────
    println!("arcsec astrometric solver version {VERSION}");
    println!(
        "Search radius: {:.0} degrees, ",
        search_radius_rad.to_degrees()
    );
    println!(
        "Start position: {}",
        format_radec(ra_hint_rad, dec_hint_rad)
    );
    println!("Image height: {:.2} degrees", fov_rad.to_degrees());
    println!("Binning: {binning}x{binning}");
    println!(
        "Image dimensions: {}x{}",
        img.width * binning,
        img.height * binning
    );
    println!("Quad tolerance: {quad_tol:.3}");
    println!("Minimum star size: {hfd_min_arcsec:.1}\"");
    println!("Speed: normal");

    // ── Solve ────────────────────────────────────────────────────────────────
    let t0 = Instant::now();

    let solve_method = match matches.get_one::<String>("method").map(|s| s.as_str()) {
        Some("tetra") => SolveMethod::Tetra,
        _ => SolveMethod::Quads,
    };

    // If --index is given, run the blind solver to estimate image position, then
    // use that estimate as the RA/Dec hint for the catalog spiral solver.
    let index_path = matches.get_one::<PathBuf>("index").cloned();

    // When blind succeeds, narrow the catalog search radius to fov_deg so the
    // spiral solver checks step 0 (blind position) and at most a few neighbours.
    let mut search_radius_rad = search_radius_rad;

    let (effective_ra, effective_dec) = if let Some(ref idx_root) = index_path {
        let index_files = collect_index_files(idx_root, fov_rad.to_degrees());
        if index_files.is_empty() {
            eprintln!("No index files found at {}", idx_root.display());
            process::exit(32);
        }

        let blind_params = BlindSolveParams {
            quad_tolerance: quad_tol,
            hfd_min,
            max_stars,
            binning,
            fov_deg: fov_rad.to_degrees(),
        };

        // Run up to BLIND_MAX_INDEXES index files in parallel so total blind time
        // equals max(t_index0, t_index1, ...) instead of the sum. Capped by the
        // thread limit so --threads 1 stays genuinely single-threaded.
        const BLIND_MAX_INDEXES: usize = 2;
        let blind_max_indexes = BLIND_MAX_INDEXES.min(arcsec_core::max_threads().max(1));
        let img_arc = Arc::new(img.clone());

        let handles: Vec<_> = index_files
            .iter()
            .take(blind_max_indexes)
            .map(|idx_path| {
                let idx_path = idx_path.clone();
                let img_ref = Arc::clone(&img_arc);
                let bp = blind_params.clone();
                std::thread::spawn(
                    move || -> Result<(f64, f64, usize), arcsec_core::ArcsecError> {
                        log::info!("Blind: trying index {}", idx_path.display());
                        let anet_index = load_anet_index(&idx_path)?;
                        let res = blind_solve(&img_ref, &anet_index, &bp);
                        if let Ok((ra, dec, score)) = &res {
                            log::info!(
                                "Blind: {} → RA={:.3}° Dec={:.3}° score={}",
                                idx_path
                                    .file_name()
                                    .map(|n| n.to_string_lossy().into_owned())
                                    .unwrap_or_default(),
                                ra.to_degrees(),
                                dec.to_degrees(),
                                score
                            );
                        }
                        res
                    },
                )
            })
            .collect();

        let mut best_blind: Option<(f64, f64, usize)> = None;
        let mut stars_ok = true;

        for handle in handles {
            match handle.join() {
                Ok(Ok((ra, dec, score))) => {
                    if best_blind.as_ref().is_none_or(|&(_, _, s)| score > s) {
                        best_blind = Some((ra, dec, score));
                    }
                }
                Ok(Err(arcsec_core::ArcsecError::InsufficientStars { found, required })) => {
                    eprintln!("Insufficient stars: found {found}, required {required}");
                    stars_ok = false;
                }
                Ok(Err(e)) => {
                    log::info!("Blind: did not solve: {e}");
                }
                Err(_) => {
                    log::info!("Blind: thread panicked");
                }
            }
        }

        if !stars_ok {
            process::exit(2);
        }

        let blind_result = best_blind.map(|(ra, dec, _score)| (ra, dec));

        match blind_result {
            Some((ra, dec)) => {
                println!(
                    "Index position estimate: RA={:.3}°, Dec={:.3}°",
                    ra.to_degrees(),
                    dec.to_degrees()
                );
                // Narrow the catalog search: blind position uncertainty is at most
                // one image-width, so 2× fov is a generous ceiling.
                search_radius_rad = (fov_rad * 2.0).max(5.0_f64.to_radians());
                (ra, dec)
            }
            None => {
                eprintln!(
                    "Blind position estimate failed for all index files. Falling back to hint."
                );
                (ra_hint_rad, dec_hint_rad)
            }
        }
    } else {
        (ra_hint_rad, dec_hint_rad)
    };

    let wcs = run_catalog_solve(
        &img,
        effective_ra,
        effective_dec,
        fov_rad,
        search_radius_rad,
        quad_tol,
        hfd_min,
        max_stars,
        db_path,
        db_name,
        binning,
        solve_method,
        threads,
    );

    let elapsed_s = t0.elapsed().as_secs_f64();

    // ── Step-distance progress line (only without -progress; ASTAP style) ────
    if !do_progress && !wcs.step_distances.is_empty() {
        let dots: Vec<String> = wcs
            .step_distances
            .iter()
            .map(|d| format!("{:.0}d", d))
            .collect();
        println!("{},", dots.join(","));
    }

    // ── Solution output ───────────────────────────────────────────────────────
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
    let dra_arcsec =
        (wcs.ra0 - ra_hint_rad) * wcs.dec0.cos() * (180.0 / core::f64::consts::PI) * 3600.0;
    let ddec_arcsec = (wcs.dec0 - dec_hint_rad) * (180.0 / core::f64::consts::PI) * 3600.0;
    println!(
        "Solved in {:.1} sec. Δ was {}.  Mount Δα={:.1}\",  Δδ={:.1}\".  Used stars down to magnitude: {:.1}",
        elapsed_s, delta_str, dra_arcsec, ddec_arcsec, wcs.mag_limit,
    );

    // ── Write output files ───────────────────────────────────────────────────
    let wcs_path = PathBuf::from(format!("{}.wcs", out_base.display()));
    let ini_path = PathBuf::from(format!("{}.ini", out_base.display()));
    if let Err(e) = fits_io::write_wcs_file(&wcs_path, &wcs) {
        eprintln!("Warning: could not write {}: {e}", wcs_path.display());
    }
    if let Err(e) = fits_io::write_ini_file(&ini_path, &wcs, max_stars) {
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

    process::exit(0);
}

/// Parse and run `arcsec catalog <subcommand>`.
fn run_catalog_cli() -> i32 {
    let m = Command::new("arcsec catalog")
        .about("Install and manage star catalogues")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .arg(
            Arg::new("dir")
                .long("dir")
                .global(true)
                .value_parser(value_parser!(PathBuf))
                .help("catalogue directory {default: the per-platform data directory}"),
        )
        .subcommand(Command::new("list").about("Show every catalogue and whether it is installed"))
        .subcommand(Command::new("path").about("Print the catalogue directory"))
        .subcommand(
            Command::new("recommend")
                .about("Suggest catalogues for a field size")
                .arg(
                    Arg::new("fov")
                        .long("fov")
                        .value_parser(value_parser!(f64))
                        .help("field of view in degrees"),
                )
                .arg(
                    Arg::new("like")
                        .long("like")
                        .value_parser(value_parser!(PathBuf))
                        .help("take the field size from this FITS file instead"),
                )
                .arg(
                    Arg::new("photometry")
                        .long("photometry")
                        .action(clap::ArgAction::SetTrue)
                        .help("also suggest a photometric catalogue (colour calibration)"),
                ),
        )
        .subcommand(
            Command::new("install")
                .about("Download and unpack catalogues")
                .arg(
                    Arg::new("names")
                        .required(true)
                        .num_args(1..)
                        .help("catalogue names, e.g. d50 v05 g05"),
                )
                .arg(
                    Arg::new("yes")
                        .long("yes")
                        .short('y')
                        .action(clap::ArgAction::SetTrue)
                        .help("do not ask for confirmation"),
                )
                .arg(
                    Arg::new("keep")
                        .long("keep-archive")
                        .action(clap::ArgAction::SetTrue)
                        .help("keep the downloaded archive after extracting"),
                ),
        )
        .subcommand(
            Command::new("remove")
                .about("Delete installed catalogue files")
                .arg(Arg::new("names").required(true).num_args(1..)),
        )
        .subcommand(
            Command::new("verify")
                .about("Check installed catalogues for missing or truncated files"),
        )
        .get_matches_from(std::env::args().skip(1));

    let dir = m
        .get_one::<PathBuf>("dir")
        .cloned()
        .unwrap_or_else(catalog_cmd::default_dir);

    let result = match m.subcommand() {
        Some(("list", _)) => {
            catalog_cmd::cmd_list(&dir);
            Ok(())
        }
        Some(("path", _)) => {
            catalog_cmd::cmd_path(&dir);
            Ok(())
        }
        Some(("recommend", sm)) => {
            let fov = match (sm.get_one::<f64>("fov"), sm.get_one::<PathBuf>("like")) {
                (Some(f), _) => Some(*f),
                (None, Some(img)) => image_io::read_pixel_scale(img).and_then(|ps| {
                    image_io::read_dimensions(img).map(|(w, h)| ps * w.max(h) as f64 / 3600.0)
                }),
                _ => None,
            };
            match fov {
                Some(f) if f > 0.0 => {
                    catalog_cmd::cmd_recommend(&dir, f, sm.get_flag("photometry"));
                    Ok(())
                }
                _ => Err(
                    "give --fov <degrees>, or --like <image.fits> with FOCALLEN and XPIXSZ in its header"
                        .to_string(),
                ),
            }
        }
        Some(("install", sm)) => {
            let names: Vec<String> = sm
                .get_many::<String>("names")
                .map(|v| v.cloned().collect())
                .unwrap_or_default();
            catalog_cmd::cmd_install(&dir, &names, sm.get_flag("yes"), sm.get_flag("keep"))
        }
        Some(("remove", sm)) => {
            let names: Vec<String> = sm
                .get_many::<String>("names")
                .map(|v| v.cloned().collect())
                .unwrap_or_default();
            catalog_cmd::cmd_remove(&dir, &names)
        }
        Some(("verify", _)) => catalog_cmd::cmd_verify(&dir),
        _ => Err("unknown subcommand".to_string()),
    };

    match result {
        Ok(()) => 0,
        Err(e) => {
            eprintln!("Error: {e}");
            1
        }
    }
}

/// Run the catalog-based (spiral search) solver and return WCS or exit.
#[allow(clippy::too_many_arguments)]
fn run_catalog_solve(
    img: &arcsec_core::types::ImageBuffer,
    ra_hint_rad: f64,
    dec_hint_rad: f64,
    fov_rad: f64,
    search_radius_rad: f64,
    quad_tol: f64,
    hfd_min: f64,
    max_stars: usize,
    db_path: std::path::PathBuf,
    db_name: String,
    binning: usize,
    method: SolveMethod,
    threads: usize,
) -> arcsec_core::types::WcsSolution {
    let params = SolveParams {
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
        threads,
    };
    match solve_image(img, &params) {
        Ok(w) => w,
        Err(arcsec_core::ArcsecError::InsufficientStars { found, required }) => {
            eprintln!("Insufficient stars: found {found}, required {required}");
            process::exit(2);
        }
        Err(arcsec_core::ArcsecError::InsufficientQuads { .. }) => {
            println!("No solution found.");
            process::exit(1);
        }
        Err(arcsec_core::ArcsecError::CatalogNotFound(p)) => {
            eprintln!("Star database not found: {}", p.display());
            process::exit(32);
        }
        Err(arcsec_core::ArcsecError::CatalogIo(e)) => {
            eprintln!("Star database read error: {e}");
            process::exit(33);
        }
        Err(e) => {
            eprintln!("Solver error: {e}");
            process::exit(1);
        }
    }
}

/// Collect and rank astrometry.net index files for blind solving.
///
/// For a single file, returns `[path]`. For a directory, peeks every
/// `index-*.fits` header, filters out files whose scale range is incompatible
/// with `fov_deg`, and returns the remainder sorted by closeness of scale
/// midpoint to `fov_deg / 2` (best match first).
///
/// When `fov_deg == 0` the full scale filter is skipped and files are returned
/// in ascending filename order.
fn collect_index_files(path: &Path, fov_deg: f64) -> Vec<PathBuf> {
    if path.is_file() {
        return vec![path.to_path_buf()];
    }
    if !path.is_dir() {
        return vec![];
    }

    let raw: Vec<PathBuf> = match fs::read_dir(path) {
        Err(_) => return vec![],
        Ok(rd) => {
            let mut v: Vec<PathBuf> = rd
                .filter_map(|e| e.ok())
                .map(|e| e.path())
                .filter(|p| {
                    p.is_file()
                        && p.file_name()
                            .and_then(|n| n.to_str())
                            .is_some_and(|s| s.starts_with("index-") && s.ends_with(".fits"))
                })
                .collect();
            v.sort();
            v
        }
    };

    if fov_deg <= 0.0 {
        return raw;
    }

    // Peek each file's scale range; filter + rank by match to FOV.
    // A "compatible" index has some overlap with [fov*0.2, fov*1.5].
    let fov_lo = fov_deg * 0.2;
    let fov_hi = fov_deg * 1.5;
    let ideal = fov_deg / 2.0;

    let mut ranked: Vec<(PathBuf, f64)> = raw
        .into_iter()
        .filter_map(|p| {
            match peek_anet_scale(&p) {
                Ok((_dq, lo_rad, hi_rad)) => {
                    let lo_deg = lo_rad.to_degrees();
                    let hi_deg = hi_rad.to_degrees();
                    // Skip indexes with no overlap in [fov_lo, fov_hi].
                    if hi_deg < fov_lo || lo_deg > fov_hi {
                        log::info!(
                            "Blind: skipping {} (scale {:.2}°–{:.2}°, outside [{:.2}°–{:.2}°])",
                            p.display(),
                            lo_deg,
                            hi_deg,
                            fov_lo,
                            fov_hi,
                        );
                        None
                    } else {
                        let mid = (lo_deg + hi_deg) / 2.0;
                        let dist = (mid - ideal).abs();
                        Some((p, dist))
                    }
                }
                Err(_) => None, // unreadable file — skip silently
            }
        })
        .collect();

    // Sort ascending by distance to ideal scale (best first).
    ranked.sort_by(|a, b| a.1.total_cmp(&b.1));
    ranked.into_iter().map(|(p, _)| p).collect()
}

/// Field-of-view range each ASTAP database is built for, and how much we prefer it
/// when several are eligible (higher = denser, so a better fit).
///
/// Ranges are as published on the ASTAP download page. The D-series stops at 6°; G05
/// and W08 exist precisely to cover wider fields, and before `.290`/`.001` support
/// they could not be read at all, which is why fields beyond ~6° never solved.
const DB_FOV_RANGES: &[(&str, f64, f64, u8)] = &[
    // prefix, min FOV (deg), max FOV (deg), preference
    ("d80", 0.15, 6.0, 8),
    ("v50", 0.20, 6.0, 6),
    ("d50", 0.20, 6.0, 5),
    ("d20", 0.30, 6.0, 4),
    ("v05", 0.60, 6.0, 3),
    ("d05", 0.60, 6.0, 2),
    ("g05", 3.00, 20.0, 7),
    ("w08", 20.0, 80.0, 7),
];

/// Every database prefix present in `db_path`, in any supported format.
fn available_dbs(db_path: &Path) -> Vec<String> {
    let mut out: Vec<String> = fs::read_dir(db_path)
        .into_iter()
        .flatten()
        .flatten()
        .filter_map(|e| {
            let n = e.file_name();
            let s = n.to_string_lossy();
            if s.ends_with(".1476") || s.ends_with(".290") || s.ends_with(".001") {
                s.split('_').next().map(|p| p.to_string())
            } else {
                None
            }
        })
        .collect();
    out.sort();
    out.dedup();
    out
}

/// Pick the installed database best suited to a `fov_deg` field.
///
/// Prefers the densest database whose published range contains the field, then any
/// database whose range merely comes closest — so an unusual field size still gets
/// the least-bad option rather than nothing.
fn select_db_for_fov(db_path: &Path, fov_deg: f64) -> Option<String> {
    let present = available_dbs(db_path);
    if present.is_empty() {
        return None;
    }

    let mut best: Option<(u8, &str)> = None;
    for &(prefix, lo, hi, pref) in DB_FOV_RANGES {
        if !present.iter().any(|p| p == prefix) {
            continue;
        }
        if fov_deg >= lo && fov_deg <= hi && best.is_none_or(|(bp, _)| pref > bp) {
            best = Some((pref, prefix));
        }
    }
    if let Some((_, prefix)) = best {
        return Some(prefix.to_string());
    }

    // Nothing covers this field: take the database whose range is nearest.
    let mut fallback: Option<(f64, &str)> = None;
    for &(prefix, lo, hi, _) in DB_FOV_RANGES {
        if !present.iter().any(|p| p == prefix) {
            continue;
        }
        let dist = if fov_deg < lo {
            lo - fov_deg
        } else {
            fov_deg - hi
        };
        if fallback.is_none_or(|(bd, _)| dist < bd) {
            fallback = Some((dist, prefix));
        }
    }
    fallback
        .map(|(_, p)| p.to_string())
        .or_else(|| present.first().cloned())
}

/// Current wall-clock time as "HH:MM:SS" (UTC seconds mod 86400).
fn hms_now() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let s = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0)
        % 86400;
    format!("{:02}:{:02}:{:02}", s / 3600, (s % 3600) / 60, s % 60)
}
