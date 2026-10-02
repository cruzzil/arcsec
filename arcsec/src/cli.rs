//! Command-line definitions: the ASTAP-compatible solver flags and the
//! `arcsec catalog` subcommands.
//!
//! The flag names, the stdout format, the output files and the exit codes are a
//! compatibility contract with ASTAP; change them only deliberately.

use std::ffi::OsString;
use std::path::PathBuf;

use clap::{Arg, ArgAction, Command, value_parser};

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// ASTAP's multi-letter options, which it spells with a single dash (`-fov 1.2`).
///
/// clap would read `-fov` as `-f ov`, so [`normalize_astap_args`] rewrites these to
/// the double-dash form before parsing. Every one of them was a parse error before
/// the rewrite, so accepting them changes nothing for the double-dash spelling.
const ASTAP_SINGLE_DASH: &[&str] = &[
    "fov", "ra", "spd", "check", "sip", "speed", "wcs", "log", "update", "progress", "analyse",
    "extract", "extract2",
];

/// Rewrite ASTAP's single-dash long options (`-fov`, `-ra`, `-spd`, `-progress`, ...)
/// to the `--fov` form clap expects, so arcsec accepts an `astap_cli` command line
/// verbatim.
///
/// Only a token that is exactly one of those names is rewritten. The first element
/// (the program name) is left alone, and so is everything after a `--`.
pub fn normalize_astap_args<I>(args: I) -> Vec<OsString>
where
    I: IntoIterator<Item = OsString>,
{
    let mut out = Vec::new();
    let mut rest_verbatim = false;
    for (i, arg) in args.into_iter().enumerate() {
        if i == 0 || rest_verbatim {
            out.push(arg);
            continue;
        }
        if arg == "--" {
            rest_verbatim = true;
            out.push(arg);
            continue;
        }
        let rewritten = arg
            .to_str()
            .and_then(|s| s.strip_prefix('-'))
            .filter(|name| ASTAP_SINGLE_DASH.contains(name))
            .map(|name| OsString::from(format!("--{name}")));
        out.push(rewritten.unwrap_or(arg));
    }
    out
}

/// The solver's command line.
pub fn solver_command() -> Command {
    Command::new("arcsec")
        .version(VERSION)
        .about("Astrometric plate solver with an ASTAP-compatible command line")
        .after_help(
            "ASTAP's single-dash spellings (-fov, -ra, -spd, -progress, ...) are accepted too.\n\
             Catalogues: run `arcsec catalog --help`.\n\
             Exit codes: 0 solved, 1 no solution or usage error, 2 too few stars, 16 file error,\n\
             32 star database or index not found, 33 star database read error.",
        )
        .arg(
            Arg::new("file")
                .short('f')
                .value_name("FILE")
                .value_parser(value_parser!(PathBuf))
                .help("Image to solve: FITS, XISF or ASDF"),
        )
        .arg(
            Arg::new("radius")
                .short('r')
                .value_name("DEG")
                .value_parser(value_parser!(f64))
                .allow_hyphen_values(true)
                .default_value("180")
                .help("Radius of the area to search, in degrees"),
        )
        .arg(
            Arg::new("fov")
                .long("fov")
                .value_name("DEG")
                .value_parser(value_parser!(f64))
                .allow_hyphen_values(true)
                .help("Field height in degrees; 0 or absent to take it from the header"),
        )
        .arg(
            Arg::new("ra")
                .long("ra")
                .value_name("HOURS")
                .value_parser(value_parser!(f64))
                .allow_hyphen_values(true)
                .help("Right ascension of the search centre, in hours"),
        )
        .arg(
            Arg::new("spd")
                .long("spd")
                .value_name("DEG")
                .value_parser(value_parser!(f64))
                .allow_hyphen_values(true)
                .help("South pole distance of the search centre (90 + Dec), in degrees"),
        )
        .arg(
            Arg::new("stars")
                .short('s')
                .value_name("N")
                .value_parser(value_parser!(usize))
                .default_value("500")
                .help("Maximum number of stars to use"),
        )
        .arg(
            Arg::new("tolerance")
                .short('t')
                .value_name("TOL")
                .value_parser(value_parser!(f64))
                .default_value("0.007")
                .help("Quad matching tolerance"),
        )
        .arg(
            Arg::new("hfd-min")
                .short('m')
                .value_name("ARCSEC")
                .value_parser(value_parser!(f64))
                .default_value("1.5")
                .help("Minimum star size (HFD), in arcseconds"),
        )
        .arg(
            Arg::new("downsample")
                .short('z')
                .value_name("FACTOR")
                .value_parser(value_parser!(u32))
                .help("Downsample (bin) factor; 0 or absent for automatic"),
        )
        .arg(
            Arg::new("check")
                .long("check")
                .value_name("y|n")
                .num_args(0..=1)
                .default_missing_value("y")
                .value_parser(["y", "n"])
                .help("Even out a raw one-shot-colour (Bayer) image before solving; for unbinned raw OSC frames only"),
        )
        .arg(
            Arg::new("database")
                .short('d')
                .value_name("DIR")
                .value_parser(value_parser!(PathBuf))
                .help("Star database directory [default: the `arcsec catalog` directory if it holds one, else the current directory]"),
        )
        .arg(
            Arg::new("db-abbrev")
                .short('D')
                .value_name("NAME")
                .help("Star database to use: d80, d50, g05, w08, ... [default: chosen by field size]"),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .value_name("BASE")
                .value_parser(value_parser!(PathBuf))
                .help("Base path for the output files [default: the image path without its extension]"),
        )
        .arg(
            Arg::new("sip")
                .long("sip")
                .value_name("y|n")
                .num_args(0..=1)
                .default_missing_value("y")
                .value_parser(["y", "n"])
                .help("Add third-order SIP (Simple Imaging Polynomial) distortion terms to the .wcs file and --update, when the field shows significant distortion"),
        )
        .arg(
            Arg::new("speed")
                .long("speed")
                .value_name("MODE")
                .value_parser(["auto", "slow"])
                .default_value("auto")
                .help("Search mode: slow reads twice the field at every search position, for more overlap"),
        )
        .arg(
            Arg::new("index")
                .long("index")
                .short('i')
                .value_name("PATH")
                .value_parser(value_parser!(PathBuf))
                .help("Blind index for solving with no position: an arcsec index (.arcsecix, see `arcsec catalog index build`), an Astrometry.net index file, or a directory holding either"),
        )
        .arg(
            Arg::new("threads")
                .long("threads")
                .value_name("N")
                .value_parser(value_parser!(usize))
                .default_value("0")
                .help("Maximum worker threads; 0 for one per core"),
        )
        .arg(
            Arg::new("method")
                .long("method")
                .value_name("METHOD")
                .value_parser(["quads", "tetra"])
                .default_value("quads")
                .help("Catalogue matching method; tetra (three-star triangles) is experimental"),
        )
        .arg(
            Arg::new("wcs")
                .long("wcs")
                .action(ArgAction::SetTrue)
                .help("Write an Astrometry.net-style .wcs file (always written; accepted for compatibility)"),
        )
        .arg(
            Arg::new("log")
                .long("log")
                .action(ArgAction::SetTrue)
                .help("Write the solver log to <BASE>.log"),
        )
        .arg(
            Arg::new("update")
                .long("update")
                .action(ArgAction::SetTrue)
                .help("Write the solution into the input FITS header"),
        )
        .arg(
            Arg::new("progress")
                .long("progress")
                .action(ArgAction::SetTrue)
                .help("Log every progress step to stderr"),
        )
        .arg(
            Arg::new("analyse")
                .long("analyse")
                .value_name("SNR_MIN")
                .value_parser(value_parser!(f64))
                .help("Analyse only, without solving: print the median HFD and number of stars (SNR_MIN 0 means 30)"),
        )
        .arg(
            Arg::new("extract")
                .long("extract")
                .value_name("SNR_MIN")
                .value_parser(value_parser!(f64))
                .help("Analyse only, and write every star found to <IMAGE>.csv (SNR_MIN 0 means 30)"),
        )
        .arg(
            Arg::new("extract2")
                .long("extract2")
                .value_name("SNR_MIN")
                .value_parser(value_parser!(f64))
                .help("Solve (with SIP), then write every star found, with its RA and Dec, to <IMAGE>.csv"),
        )
}

/// `arcsec catalog ...`.
pub fn catalog_command() -> Command {
    Command::new("catalog")
        .bin_name("arcsec catalog")
        .about("Install and manage star catalogues")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .arg(
            Arg::new("dir")
                .long("dir")
                .global(true)
                .value_name("DIR")
                .value_parser(value_parser!(PathBuf))
                .help("Catalogue directory [default: $ARCSEC_CATALOG_DIR, else the per-platform data directory]"),
        )
        .subcommand(Command::new("list").about("Show every catalogue and whether it is installed"))
        .subcommand(Command::new("path").about("Print the catalogue directory"))
        .subcommand(
            Command::new("recommend")
                .about("Suggest catalogues for a field size")
                .arg(
                    Arg::new("fov")
                        .long("fov")
                        .value_name("DEG")
                        .value_parser(value_parser!(f64))
                        .help("Field of view in degrees"),
                )
                .arg(
                    Arg::new("like")
                        .long("like")
                        .value_name("IMAGE")
                        .value_parser(value_parser!(PathBuf))
                        .help("Take the field size from this image's header instead"),
                )
                .arg(
                    Arg::new("photometry")
                        .long("photometry")
                        .action(ArgAction::SetTrue)
                        .help("Also suggest a photometric catalogue (colour calibration)"),
                ),
        )
        .subcommand(
            Command::new("install")
                .about("Download and unpack catalogues")
                .arg(
                    Arg::new("names")
                        .value_name("NAME")
                        .required(true)
                        .num_args(1..)
                        .help("Catalogue names, e.g. d50 v05 g05"),
                )
                .arg(
                    Arg::new("yes")
                        .long("yes")
                        .short('y')
                        .action(ArgAction::SetTrue)
                        .help("Do not ask for confirmation"),
                )
                .arg(
                    Arg::new("keep")
                        .long("keep-archive")
                        .action(ArgAction::SetTrue)
                        .help("Keep the downloaded archive after extracting"),
                )
                .arg(
                    Arg::new("no-index")
                        .long("no-index")
                        .action(ArgAction::SetTrue)
                        .conflicts_with_all(["index-min-fov", "index-max-fov"])
                        .help("Do not build the blind index after installing a solving database"),
                )
                .arg(
                    Arg::new("index-min-fov")
                        .long("index-min-fov")
                        .value_name("DEG")
                        .value_parser(value_parser!(f64))
                        .help("Smallest field (short side, degrees) the blind index serves [default: from the databases, e.g. 0.3 for D80]"),
                )
                .arg(
                    Arg::new("index-max-fov")
                        .long("index-max-fov")
                        .value_name("DEG")
                        .value_parser(value_parser!(f64))
                        .help("Largest field (short side, degrees) the blind index serves [default: from the databases, 30 or 80 with W08]"),
                )
                .after_help(
                    "Installing a solving database (d05, d20, d50, d80, g05, w08) also builds arcsec's \
                     blind index from it, which lets the solver find a field with no position hint. \
                     The confirmation shows its size, build time and memory first; a large build \
                     (over 1 GB, over 5 minutes, or more than half the free memory) is asked about \
                     separately.",
                ),
        )
        .subcommand(
            Command::new("remove")
                .about("Delete installed catalogue files")
                .arg(
                    Arg::new("names")
                        .value_name("NAME")
                        .required(true)
                        .num_args(1..)
                        .help("Catalogue names, as shown by `arcsec catalog list`"),
                )
                .arg(
                    Arg::new("yes")
                        .long("yes")
                        .short('y')
                        .action(ArgAction::SetTrue)
                        .help("Do not ask for confirmation"),
                )
                .arg(
                    Arg::new("keep-index")
                        .long("keep-index")
                        .action(ArgAction::SetTrue)
                        .help("Keep the blind index built from a removed database"),
                ),
        )
        .subcommand(
            Command::new("verify")
                .about("Check installed catalogues for missing or truncated files"),
        )
        .subcommand(
            Command::new("index")
                .about("Build or inspect arcsec's blind index (no download: built from an installed star database)")
                .subcommand_required(true)
                .arg_required_else_help(true)
                .subcommand(
                    Command::new("build")
                        .about("Build a blind index from an installed star database")
                        .arg(
                            Arg::new("db")
                                .long("db")
                                .value_name("DIR")
                                .value_parser(value_parser!(PathBuf))
                                .help("Directory holding the star database [default: the catalogue directory]"),
                        )
                        .arg(
                            Arg::new("name")
                                .long("name")
                                .short('D')
                                .value_name("DB")
                                .help("Database to build from, e.g. d80 [default: the deepest installed]"),
                        )
                        .arg(
                            Arg::new("min-fov")
                                .long("min-fov")
                                .value_name("DEG")
                                .value_parser(value_parser!(f64))
                                .help("Smallest field (short side, degrees) to support; smaller fields need much larger indexes [default: from the installed databases: 0.3 for D20-D80, 0.6 D05, 3 G05, 10 W08]"),
                        )
                        .arg(
                            Arg::new("max-fov")
                                .long("max-fov")
                                .value_name("DEG")
                                .value_parser(value_parser!(f64))
                                .help("Largest field (short side, degrees) to support [default: 30, or 80 with W08]"),
                        )
                        .arg(
                            Arg::new("yes")
                                .long("yes")
                                .short('y')
                                .action(ArgAction::SetTrue)
                                .help("Do not ask before a large build (over 1 GB, 5 minutes, or half the free memory)"),
                        )
                        .arg(
                            Arg::new("out")
                                .long("out")
                                .short('o')
                                .value_name("FILE")
                                .value_parser(value_parser!(PathBuf))
                                .help("Output file [default: <catalogue dir>/<db>.arcsecix]"),
                        )
                        .arg(
                            Arg::new("threads")
                                .long("threads")
                                .value_name("N")
                                .value_parser(value_parser!(usize))
                                .default_value("0")
                                .help("Worker threads; 0 for one per core"),
                        ),
                )
                .subcommand(
                    Command::new("info")
                        .about("Describe a blind index: tiers, sizes, source")
                        .arg(
                            Arg::new("file")
                                .value_name("FILE")
                                .value_parser(value_parser!(PathBuf))
                                .help("Index file [default: every *.arcsecix in the catalogue directory]"),
                        ),
                ),
        )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn norm(args: &[&str]) -> Vec<String> {
        normalize_astap_args(args.iter().map(OsString::from))
            .into_iter()
            .map(|s| s.into_string().unwrap())
            .collect()
    }

    #[test]
    fn command_definitions_are_valid() {
        solver_command().debug_assert();
        catalog_command().debug_assert();
    }

    #[test]
    fn astap_single_dash_options_are_rewritten() {
        assert_eq!(
            norm(&[
                "arcsec",
                "-f",
                "x.fits",
                "-fov",
                "1.2",
                "-ra",
                "5.5",
                "-spd",
                "97",
                "-progress"
            ]),
            [
                "arcsec",
                "-f",
                "x.fits",
                "--fov",
                "1.2",
                "--ra",
                "5.5",
                "--spd",
                "97",
                "--progress"
            ]
        );
    }

    #[test]
    fn everything_else_is_left_alone() {
        let args = [
            "arcsec",
            "-f",
            "-fov.fits",
            "-r",
            "-30",
            "--fov",
            "0",
            "-z",
            "2",
            "-D",
            "d50",
        ];
        assert_eq!(norm(&args), args);
        // After `--`, nothing is an option.
        assert_eq!(norm(&["arcsec", "--", "-fov"]), ["arcsec", "--", "-fov"]);
        // The program name is never rewritten.
        assert_eq!(norm(&["-log"]), ["-log"]);
    }

    #[test]
    fn astap_command_line_parses() {
        let m = solver_command()
            .try_get_matches_from(normalize_astap_args(
                [
                    "astap_cli",
                    "-f",
                    "a.fits",
                    "-fov",
                    "0",
                    "-r",
                    "30",
                    "-ra",
                    "5.5",
                    "-spd",
                    "97",
                    "-z",
                    "0",
                    "-s",
                    "500",
                    "-wcs",
                    "-log",
                    "-update",
                    "-speed",
                    "auto",
                ]
                .map(OsString::from),
            ))
            .expect("an astap_cli command line must parse");
        assert_eq!(m.get_one::<f64>("fov"), Some(&0.0));
        assert_eq!(m.get_one::<f64>("spd"), Some(&97.0));
        assert!(m.get_flag("wcs") && m.get_flag("log") && m.get_flag("update"));
    }

    fn parse(args: &[&str]) -> Result<clap::ArgMatches, clap::Error> {
        solver_command().try_get_matches_from(normalize_astap_args(
            args.iter().copied().map(OsString::from),
        ))
    }

    #[test]
    fn sip_and_check_take_an_optional_y_or_n() {
        let m = parse(&["astap_cli", "-f", "a.fits", "-sip", "-check", "y", "-wcs"]).unwrap();
        assert_eq!(m.get_one::<String>("sip").map(String::as_str), Some("y"));
        assert_eq!(m.get_one::<String>("check").map(String::as_str), Some("y"));
        assert!(m.get_flag("wcs"), "a following flag is not the value");
        let m = parse(&["astap_cli", "-sip", "n", "-check"]).unwrap();
        assert_eq!(m.get_one::<String>("sip").map(String::as_str), Some("n"));
        assert_eq!(m.get_one::<String>("check").map(String::as_str), Some("y"));
        let m = parse(&["astap_cli", "-f", "a.fits"]).unwrap();
        assert!(m.get_one::<String>("sip").is_none() && m.get_one::<String>("check").is_none());
        assert!(parse(&["astap_cli", "-sip", "maybe"]).is_err());
    }

    #[test]
    fn speed_is_auto_or_slow() {
        let m = parse(&["astap_cli", "-speed", "slow"]).unwrap();
        assert_eq!(
            m.get_one::<String>("speed").map(String::as_str),
            Some("slow")
        );
        let m = parse(&["astap_cli"]).unwrap();
        assert_eq!(
            m.get_one::<String>("speed").map(String::as_str),
            Some("auto")
        );
        assert!(parse(&["astap_cli", "-speed", "fast"]).is_err());
    }

    #[test]
    fn analyse_options_take_a_minimum_snr() {
        let m = parse(&[
            "astap_cli",
            "-f",
            "a.fits",
            "-analyse",
            "30",
            "-extract",
            "0",
        ])
        .unwrap();
        assert_eq!(m.get_one::<f64>("analyse"), Some(&30.0));
        assert_eq!(m.get_one::<f64>("extract"), Some(&0.0));
        let m = parse(&["astap_cli", "-extract2", "20"]).unwrap();
        assert_eq!(m.get_one::<f64>("extract2"), Some(&20.0));
    }

    #[test]
    fn an_unknown_method_is_rejected() {
        let r = solver_command().try_get_matches_from(["arcsec", "--method", "triangles"]);
        assert!(r.is_err());
    }
}
