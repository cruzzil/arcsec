//! Conversion of plate constants into a FITS WCS solution, and SIP distortion.

pub mod output;
pub mod sip;

pub use output::derive_wcs;
pub use sip::{Sip, TanWcs, fit_sip};
