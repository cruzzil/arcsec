//! Download progress on stderr.

use std::io::{self, IsTerminal as _, Write as _};
use std::time::Instant;

use arcsec_catalogue::fetch::DownloadEvent;
use arcsec_catalogue::human_bytes as human;

/// Draws one download's progress. On a terminal the line is redrawn in place a
/// few times a second; redirected to a file or a pipe, where carriage returns do
/// not overwrite, it is a line every 15 s, so a 117 MB download is not hundreds of
/// lines of noise.
#[derive(Default)]
pub struct Progress {
    /// Set when a transfer starts: whether stderr is a terminal, and when progress
    /// was last drawn.
    state: Option<(bool, Instant)>,
}

impl Progress {
    /// Report `event` of the download called `label`.
    pub fn event(&mut self, label: &str, event: DownloadEvent) {
        match event {
            DownloadEvent::Started { .. } => {
                self.state = Some((io::stderr().is_terminal(), Instant::now()));
            }
            DownloadEvent::Progress { written, total } => {
                let Some((tty, last_report)) = &mut self.state else {
                    return;
                };
                let (cr, tail) = ends(*tty);
                let interval_ms = if *tty { 300 } else { 15_000 };
                if last_report.elapsed().as_millis() > interval_ms {
                    match total {
                        Some(t) if t > 0 => eprint!(
                            "{cr}  {label}: {} / {} ({:.0}%){tail}",
                            human(written),
                            human(t),
                            written as f64 / t as f64 * 100.0,
                        ),
                        _ => eprint!("{cr}  {label}: {}{tail}", human(written)),
                    }
                    let _ = io::stderr().flush();
                    *last_report = Instant::now();
                }
            }
            DownloadEvent::Finished { written } => {
                let tty = self.state.is_some_and(|(tty, _)| tty);
                let (cr, _) = ends(tty);
                eprintln!("{cr}  {label}: {} downloaded          ", human(written));
            }
            DownloadEvent::AlreadyComplete { bytes } => {
                eprintln!("  {label}: already downloaded ({})", human(bytes));
            }
        }
    }
}

/// What starts and ends a progress line: a carriage return to redraw in place on
/// a terminal, a newline per report elsewhere.
const fn ends(tty: bool) -> (&'static str, &'static str) {
    if tty { ("\r", "   ") } else { ("", "\n") }
}
