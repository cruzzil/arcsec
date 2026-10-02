//! Yes/no questions, and the rules for when to ask them.
//!
//! The rules are pure functions over a [`Prompter`], so the non-interactive paths
//! (`--yes`, no terminal, end of input) are tested without a terminal.

use std::io::{IsTerminal as _, Write as _};

/// Something that can answer a yes/no question.
pub trait Prompter {
    /// Ask `question`; anything but an explicit yes is a no.
    fn ask(&mut self, question: &str) -> bool;
    /// Whether a person is there to answer (stdin is a terminal).
    fn interactive(&self) -> bool;
}

/// The real terminal: asks on stderr, reads stdin. End of input (a pipe, a
/// scheduler, CI) reads as no, as it always has for `install` and `remove`.
pub struct Terminal;

impl Prompter for Terminal {
    fn ask(&mut self, question: &str) -> bool {
        eprint!("{question} [y/N] ");
        let _ = std::io::stderr().flush();
        let mut line = String::new();
        std::io::stdin().read_line(&mut line).is_ok()
            && matches!(line.trim().to_ascii_lowercase().as_str(), "y" | "yes")
    }

    fn interactive(&self) -> bool {
        std::io::stdin().is_terminal()
    }
}

/// What `catalog install` will do after asking.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InstallAnswer {
    /// Download the catalogues.
    pub download: bool,
    /// Build the blind index afterwards.
    pub build_index: bool,
}

/// The questions `catalog install` asks, all before anything is downloaded so that
/// nobody comes back from a long download to find a question waiting.
///
/// * With `--yes`, everything planned goes ahead (the notices are still printed).
/// * Downloads are confirmed by one question, which also covers an ordinary index
///   build listed in the same summary.
/// * An index build with `concerns` (big, slow, or short of memory) gets a question
///   of its own, so the download can be accepted and the index declined.
/// * With nothing to download, the index build is the only question.
pub fn install_questions(
    downloads: bool,
    index_planned: bool,
    has_concerns: bool,
    assume_yes: bool,
    p: &mut dyn Prompter,
) -> InstallAnswer {
    if assume_yes {
        return InstallAnswer {
            download: downloads,
            build_index: index_planned,
        };
    }
    if downloads && !p.ask("Continue?") {
        return InstallAnswer {
            download: false,
            build_index: false,
        };
    }
    let build_index = index_planned
        && if downloads && !has_concerns {
            true
        } else if downloads {
            p.ask("Build the blind index as well?")
        } else {
            p.ask("Build the blind index now?")
        };
    InstallAnswer {
        download: downloads,
        build_index,
    }
}

/// Whether `catalog index build` goes ahead. An ordinary build just runs, as it
/// always has. One with concerns asks first when someone is there to answer; with
/// no terminal it goes ahead after printing the notice, so scripts that already
/// run `index build` keep working.
pub fn build_question(has_concerns: bool, assume_yes: bool, p: &mut dyn Prompter) -> bool {
    if !has_concerns || assume_yes || !p.interactive() {
        return true;
    }
    p.ask("Build it?")
}

#[cfg(test)]
pub mod tests {
    use super::*;

    /// Answers from a script; `None` once it runs out, like end of input.
    pub struct Scripted {
        pub answers: Vec<bool>,
        pub asked: Vec<String>,
        pub tty: bool,
    }

    impl Scripted {
        pub fn new(answers: &[bool], tty: bool) -> Self {
            Self {
                answers: answers.iter().rev().copied().collect(),
                asked: Vec::new(),
                tty,
            }
        }
    }

    impl Prompter for Scripted {
        fn ask(&mut self, q: &str) -> bool {
            self.asked.push(q.to_string());
            self.answers.pop().unwrap_or(false)
        }
        fn interactive(&self) -> bool {
            self.tty
        }
    }

    fn ans(download: bool, build_index: bool) -> InstallAnswer {
        InstallAnswer {
            download,
            build_index,
        }
    }

    #[test]
    fn yes_skips_every_question() {
        let mut p = Scripted::new(&[], false);
        assert_eq!(
            install_questions(true, true, true, true, &mut p),
            ans(true, true)
        );
        assert_eq!(
            install_questions(false, true, false, true, &mut p),
            ans(false, true)
        );
        assert!(build_question(true, true, &mut p));
        assert!(p.asked.is_empty());
    }

    #[test]
    fn an_ordinary_index_rides_on_the_download_question() {
        let mut p = Scripted::new(&[true], true);
        assert_eq!(
            install_questions(true, true, false, false, &mut p),
            ans(true, true)
        );
        assert_eq!(p.asked, ["Continue?"]);
    }

    #[test]
    fn a_big_index_gets_its_own_question_and_can_be_declined() {
        let mut p = Scripted::new(&[true, false], true);
        assert_eq!(
            install_questions(true, true, true, false, &mut p),
            ans(true, false)
        );
        assert_eq!(p.asked.len(), 2);
        let mut p = Scripted::new(&[true, true], true);
        assert_eq!(
            install_questions(true, true, true, false, &mut p),
            ans(true, true)
        );
    }

    #[test]
    fn declining_the_download_declines_everything() {
        let mut p = Scripted::new(&[false], true);
        assert_eq!(
            install_questions(true, true, true, false, &mut p),
            ans(false, false)
        );
        assert_eq!(p.asked.len(), 1, "no second question after a no");
    }

    #[test]
    fn end_of_input_without_yes_cancels_as_before() {
        // Piped or scheduled, no --yes: every question reads as no.
        let mut p = Scripted::new(&[], false);
        assert_eq!(
            install_questions(true, true, false, false, &mut p),
            ans(false, false)
        );
        let mut p = Scripted::new(&[], false);
        assert_eq!(
            install_questions(false, true, false, false, &mut p),
            ans(false, false)
        );
    }

    #[test]
    fn index_only_installs_ask_once() {
        let mut p = Scripted::new(&[true], true);
        assert_eq!(
            install_questions(false, true, false, false, &mut p),
            ans(false, true)
        );
        assert_eq!(p.asked, ["Build the blind index now?"]);
    }

    #[test]
    fn index_build_asks_only_about_big_builds_and_only_at_a_terminal() {
        let mut p = Scripted::new(&[], true);
        assert!(
            build_question(false, false, &mut p),
            "ordinary: no question"
        );
        assert!(p.asked.is_empty());
        let mut p = Scripted::new(&[false], true);
        assert!(!build_question(true, false, &mut p));
        let mut p = Scripted::new(&[], false);
        assert!(build_question(true, false, &mut p), "scripts keep working");
        assert!(p.asked.is_empty());
    }
}
