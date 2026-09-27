//! The solver's logger: `env_logger` for `--progress` on stderr, plus an optional
//! file sink for `--log`.
//!
//! With neither flag the `log` crate is switched off entirely (`LevelFilter::Off`),
//! so the `log::info!` calls throughout `arcsec-core` cost nothing by default.

use alloc::sync::Arc;
use std::fs;
use std::io::Write as _;
use std::path::Path;
use std::sync::Mutex;

use env_logger::Logger as EnvLogger;
use log::{LevelFilter, Log, Metadata, Record};

struct ArcsecLogger {
    /// `env_logger` backend: handles format and stderr (`Some` with `--progress`).
    stderr: Option<EnvLogger>,
    /// Optional file sink for `--log`.
    file: Option<Arc<Mutex<fs::File>>>,
}

impl Log for ArcsecLogger {
    fn enabled(&self, m: &Metadata<'_>) -> bool {
        self.stderr
            .as_ref()
            .map_or_else(|| self.file.is_some(), |l| l.enabled(m))
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

/// Install the global logger.
///
/// `log_path` is the `--log` destination, if requested; a file that cannot be
/// created is reported as a warning and logging continues without it.
pub fn install(progress: bool, log_path: Option<&Path>) {
    let file = log_path.and_then(|lp| match fs::File::create(lp) {
        Ok(f) => Some(Arc::new(Mutex::new(f))),
        Err(e) => {
            eprintln!("Warning: cannot create {}: {e}", lp.display());
            None
        }
    });

    let stderr = progress.then(|| {
        let mut b = env_logger::Builder::new();
        b.filter_level(LevelFilter::Info);
        b.format(|buf, r| writeln!(buf, "{}  {}", hms_now(), r.args()));
        b.build()
    });

    log::set_boxed_logger(Box::new(ArcsecLogger { stderr, file })).ok();
    log::set_max_level(if progress || log_path.is_some() {
        LevelFilter::Info
    } else {
        LevelFilter::Off
    });
}

/// Current wall-clock time as "HH:MM:SS" (UTC seconds mod 86400).
fn hms_now() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let s = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_secs())
        % 86400;
    format!("{:02}:{:02}:{:02}", s / 3600, (s % 3600) / 60, s % 60)
}
