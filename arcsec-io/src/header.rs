//! FITS header text held in memory.
//!
//! An application embedding the solver usually has the image's header already —
//! Siril keeps it as one string of 80-character cards, others as lines — and wants
//! the solver to read pointing and pixel scale from it exactly as it would from the
//! file. [`HeaderCards`] parses that text and answers with the same rules the file
//! readers use ([`crate::image_io::ra_dec_from`], [`crate::image_io::pixel_scale_from`],
//! [`crate::image_io::tan_wcs_from`]).

use arcsec_core::wcs::TanWcs;

use crate::image_io;

/// Length of one FITS header card.
const CARD: usize = 80;

/// The keyword cards of a FITS header, parsed from text.
#[derive(Debug, Clone, Default)]
pub struct HeaderCards {
    /// `(KEYWORD, raw value)` in header order; the value is the text between `= `
    /// and the comment, trimmed.
    cards: Vec<(String, String)>,
}

impl HeaderCards {
    /// Parse header text: 80-character cards run together (as CFITSIO's
    /// `fits_hdr2str` writes them), or one card per line. Anything that is not a
    /// `KEYWORD = value` card (COMMENT, HISTORY, blank, END) is skipped, as is
    /// everything after `END`.
    #[must_use]
    pub fn parse(text: &str) -> Self {
        let lines: Vec<&str> = if text.contains('\n') {
            text.lines().collect()
        } else {
            let mut v = Vec::new();
            let mut rest = text;
            while !rest.is_empty() {
                let mut cut = rest.len().min(CARD);
                while !rest.is_char_boundary(cut) {
                    cut -= 1;
                }
                if cut == 0 {
                    break;
                }
                let (card, tail) = rest.split_at(cut);
                v.push(card);
                rest = tail;
            }
            v
        };
        let mut cards = Vec::new();
        for line in lines {
            let line = line.trim_end_matches('\r');
            let key = line.get(..8.min(line.len())).unwrap_or("").trim();
            if key == "END" {
                break;
            }
            // A value card has "= " in columns 9-10.
            let Some(rest) = line.get(8..).and_then(|r| r.strip_prefix('=')) else {
                continue;
            };
            if key.is_empty() {
                continue;
            }
            cards.push((key.to_ascii_uppercase(), value_text(rest).to_string()));
        }
        Self { cards }
    }

    /// Number of keyword cards.
    #[must_use]
    pub fn len(&self) -> usize {
        self.cards.len()
    }

    /// Whether there are no keyword cards.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.cards.is_empty()
    }

    /// The raw value of `key` (case-insensitive), as written: quotes kept.
    #[must_use]
    pub fn raw(&self, key: &str) -> Option<&str> {
        self.cards
            .iter()
            .find(|(k, _)| k.eq_ignore_ascii_case(key))
            .map(|(_, v)| v.as_str())
    }

    /// A numeric keyword. A quoted number counts (CFITSIO converts those too), and
    /// so does a Fortran `D` exponent.
    #[must_use]
    pub fn number(&self, key: &str) -> Option<f64> {
        let v = self.raw(key)?.trim().trim_matches('\'').trim();
        v.parse::<f64>()
            .ok()
            .or_else(|| v.replace(['D', 'd'], "E").parse::<f64>().ok())
            .filter(|x| x.is_finite())
    }

    /// Pointing in degrees: RA/DEC, else CRVAL1/CRVAL2, as the file readers read it.
    #[must_use]
    pub fn ra_dec(&self) -> Option<(f64, f64)> {
        image_io::ra_dec_from(
            self.number("RA"),
            self.number("DEC"),
            self.number("CRVAL1"),
            self.number("CRVAL2"),
        )
    }

    /// Pixel scale in arcseconds per pixel from FOCALLEN, XPIXSZ and XBINNING.
    #[must_use]
    pub fn pixel_scale(&self) -> Option<f64> {
        image_io::pixel_scale_from(
            self.number("FOCALLEN"),
            self.number("XPIXSZ"),
            self.number("XBINNING"),
        )
    }

    /// A TAN WCS already in the header, if it has a CD matrix.
    #[must_use]
    pub fn tan_wcs(&self) -> Option<TanWcs> {
        image_io::tan_wcs_from(|k| self.number(k))
    }
}

/// The value part of a card after `=`: up to the `/` that starts the comment,
/// outside any quoted string.
fn value_text(rest: &str) -> &str {
    let mut in_str = false;
    for (i, c) in rest.char_indices() {
        match c {
            '\'' => in_str = !in_str,
            '/' if !in_str => return rest[..i].trim(),
            _ => {}
        }
    }
    rest.trim()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn card(s: &str) -> String {
        format!("{s:<80}")
    }

    #[test]
    fn parses_run_together_cards() {
        let text = [
            "SIMPLE  =                    T / conforms",
            "RA      =     83.8221 / [deg] pointing",
            "DEC     = -5.3911 / [deg]",
            "OBJECT  = 'M42 / Orion' / a slash inside a string",
            "FOCALLEN=                530.0",
            "XPIXSZ  =                 3.76",
            "XBINNING=                    2",
            "COMMENT  RA = 1",
            "END",
            "DEC     = 99",
        ]
        .map(card)
        .concat();
        let h = HeaderCards::parse(&text);
        assert_eq!(h.ra_dec(), Some((83.8221, -5.3911)));
        assert_eq!(h.raw("object"), Some("'M42 / Orion'"));
        let ps = h.pixel_scale().unwrap();
        assert!((ps - 3.76 * 2.0 / 530.0 * 206.265).abs() < 1e-12);
        assert!(h.tan_wcs().is_none());
        assert_eq!(h.len(), 7, "COMMENT and anything after END are skipped");
    }

    #[test]
    fn parses_lines_and_odd_numbers() {
        let text = "CRVAL1  = '150.5'\nCRVAL2  = 2.25D+01 / Fortran exponent\r\nCRPIX1  = 100\nCRPIX2  = 80\nCD1_1   = -1E-4\nCD2_2   = 1E-4\n";
        let h = HeaderCards::parse(text);
        assert_eq!(h.ra_dec(), Some((150.5, 22.5)));
        let w = h.tan_wcs().unwrap();
        assert_eq!(w.cd, [[-1e-4, 0.0], [0.0, 1e-4]]);
        assert_eq!(h.number("CRPIX1"), Some(100.0));
        assert!(h.pixel_scale().is_none());
    }

    #[test]
    fn junk_is_harmless() {
        for text in [
            "",
            "=",
            "é",
            &"x".repeat(1000),
            "RA      = 'abc'\nDEC     = nan",
        ] {
            let h = HeaderCards::parse(text);
            assert!(h.ra_dec().is_none());
        }
        // Multi-byte text that does not fall on card boundaries.
        let h = HeaderCards::parse(&"é".repeat(100));
        assert!(h.is_empty());
    }
}
