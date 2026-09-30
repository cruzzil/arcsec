//! XISF input, via the `xisf` crate.
//!
//! XISF is PixInsight's native format: an XML header describing one or more
//! images, followed by their sample data. Two things make it easy to slot in
//! beside FITS:
//!
//! * files written by PixInsight carry the original `<FITSKeyword>` elements, so
//!   the pointing and optics metadata is read with the same rules as FITS;
//! * the crate hands back typed samples, so the only work here is converting
//!   whatever sample format the file uses into the `f32` the solver wants.

use std::path::Path;

use arcsec_core::types::ImageBuffer;
use xisf::{ImageRef, SampleFormat, XisfFile};

use crate::image_io;

/// Open a file and hand back its first image, which is the one to solve.
///
/// A monolithic XISF may hold several images (PixInsight writes masks and
/// previews alongside the light), and the first `<Image>` element is the primary
/// one by the spec's ordering.
fn open_primary(path: &Path) -> Result<(XisfFile, usize), String> {
    let file = XisfFile::open(path).map_err(|e| format!("XISF open failed: {e}"))?;
    if file.images().is_empty() {
        return Err(format!("{}: XISF file contains no images", path.display()));
    }
    Ok((file, 0))
}

/// Width and height of an image element.
///
/// `geometry()` is fastest-varying first and excludes the channel count, so for a
/// 2-D image it reads `[width, height]`.
fn dimensions_of(image: &ImageRef<'_>) -> Result<(usize, usize), String> {
    let g = image.geometry();
    if g.len() < 2 {
        return Err(format!(
            "XISF image has {} dimension(s), need at least 2",
            g.len()
        ));
    }
    let w = g[0] as usize;
    let h = g[1] as usize;
    if w == 0 || h == 0 {
        return Err(format!("invalid XISF image dimensions {w}×{h}"));
    }
    Ok((w, h))
}

/// Read the samples as `f32`, whatever they are stored as.
///
/// `ImageRef::read` deliberately refuses to reinterpret one sample format as
/// another, so the format is matched first and converted here. Integer samples are
/// left on their own scale rather than normalised: the detector's background and
/// noise estimate work from the raw distribution, and
/// `ImageBuffer::normalize_for_detection` rescales whatever arrives.
fn samples_as_f32(image: &ImageRef<'_>) -> Result<Vec<f32>, String> {
    let e = |err: xisf::Error| format!("XISF read failed: {err}");
    Ok(match image.sample_format() {
        SampleFormat::UInt8 => image
            .read_planar::<u8>()
            .map_err(e)?
            .into_iter()
            .map(f32::from)
            .collect(),
        SampleFormat::UInt16 => image
            .read_planar::<u16>()
            .map_err(e)?
            .into_iter()
            .map(f32::from)
            .collect(),
        SampleFormat::UInt32 => image
            .read_planar::<u32>()
            .map_err(e)?
            .into_iter()
            .map(|v| v as f32)
            .collect(),
        SampleFormat::UInt64 => image
            .read_planar::<u64>()
            .map_err(e)?
            .into_iter()
            .map(|v| v as f32)
            .collect(),
        SampleFormat::Float32 => image.read_planar::<f32>().map_err(e)?,
        SampleFormat::Float64 => image
            .read_planar::<f64>()
            .map_err(e)?
            .into_iter()
            .map(|v| v as f32)
            .collect(),
        other => {
            return Err(format!(
                "XISF sample format {other:?} is not supported for solving \
                 (complex images have no single brightness per pixel)"
            ));
        }
    })
}

/// Load an XISF image as a greyscale buffer.
///
/// Colour images are averaged across channels rather than reduced to one, which
/// costs nothing and gives the detector a better signal-to-noise ratio than any
/// single channel would. `read_planar` returns channels whole and in order
/// regardless of how the file interleaves them, so the average is a strided sum.
pub fn read_xisf_image(path: &Path) -> Result<ImageBuffer, String> {
    let (file, idx) = open_primary(path)?;
    let images = file.images();
    let image = &images[idx];

    let (width, height) = dimensions_of(image)?;
    let channels = image.channels().max(1) as usize;
    let npix = width
        .checked_mul(height)
        .ok_or_else(|| format!("XISF image dimensions {width}×{height} overflow"))?;

    let mut samples = samples_as_f32(image)?;

    // `offset` is a pedestal added to every sample, which the spec says must be
    // subtracted to get zero-based values. It is uniform, so it does not change
    // which pixels are stars, but leaving it in would misreport the background.
    if let Some(pedestal) = image.attributes().offset
        && pedestal != 0.0
    {
        let p = pedestal as f32;
        for v in &mut samples {
            *v -= p;
        }
    }
    if samples.len() < npix * channels {
        return Err(format!(
            "XISF image declares {width}×{height}×{channels} samples but only {} were read",
            samples.len()
        ));
    }

    let data = if channels == 1 {
        samples
    } else {
        let inv = 1.0 / channels as f32;
        (0..npix)
            .map(|i| (0..channels).map(|c| samples[c * npix + i]).sum::<f32>() * inv)
            .collect()
    };

    Ok(ImageBuffer {
        data,
        width,
        height,
    })
}

/// Dimensions without decoding the pixel data.
pub fn read_xisf_dimensions(path: &Path) -> Option<(u32, u32)> {
    let (file, idx) = open_primary(path).ok()?;
    let images = file.images();
    let (w, h) = dimensions_of(&images[idx]).ok()?;
    Some((w as u32, h as u32))
}

/// A numeric FITS keyword of the image, matched case-insensitively.
///
/// PixInsight preserves the keywords of whatever it opened, so a light frame
/// converted from FITS still carries RA, DEC, FOCALLEN and the rest.
fn keyword(image: &ImageRef<'_>, name: &str) -> Option<f64> {
    image
        .fits_keywords()
        .into_iter()
        .find(|(k, _, _)| k.trim().eq_ignore_ascii_case(name))
        .and_then(|(_, v, _)| parse_keyword_value(&v))
}

/// Parse a FITS keyword value written as text.
///
/// XISF stores keyword values as they appeared in the FITS card, so a string
/// value keeps its quotes and a numeric one may carry trailing comment padding.
fn parse_keyword_value(raw: &str) -> Option<f64> {
    let v = raw.trim().trim_matches('\'').trim();
    v.parse::<f64>().ok()
}

/// Pointing from the file's FITS keywords, in degrees.
pub fn read_xisf_ra_dec(path: &Path) -> Option<(f64, f64)> {
    let (file, idx) = open_primary(path).ok()?;
    let images = file.images();
    let img = &images[idx];
    image_io::ra_dec_from(
        keyword(img, "RA"),
        keyword(img, "DEC"),
        keyword(img, "CRVAL1"),
        keyword(img, "CRVAL2"),
    )
}

/// A TAN WCS from the file's FITS keywords, if it carries one.
pub fn read_xisf_header_wcs(path: &Path) -> Option<arcsec_core::wcs::TanWcs> {
    let (file, idx) = open_primary(path).ok()?;
    let images = file.images();
    let img = &images[idx];
    image_io::tan_wcs_from(|name| keyword(img, name))
}

/// Colour channels of the primary image.
pub fn read_xisf_channels(path: &Path) -> usize {
    open_primary(path).map_or(1, |(file, idx)| {
        file.images()[idx].channels().max(1) as usize
    })
}

/// Plate scale in arcsec/pixel from the file's FITS keywords.
pub fn read_xisf_pixel_scale(path: &Path) -> Option<f64> {
    let (file, idx) = open_primary(path).ok()?;
    let images = file.images();
    let img = &images[idx];
    image_io::pixel_scale_from(
        keyword(img, "FOCALLEN"),
        keyword(img, "XPIXSZ"),
        keyword(img, "XBINNING"),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use xisf::{
        BlockOptions, Bounds, ColorSpace, FitsKeyword, Image, PendingImage, PixelStorage,
        SampleFormat, Writer,
    };

    /// Build an XISF file in memory and write it to a temporary path.
    fn write_xisf(
        name: &str,
        w: u64,
        h: u64,
        channels: u64,
        format: SampleFormat,
        bytes: Vec<u8>,
        keywords: Vec<FitsKeyword>,
    ) -> std::path::PathBuf {
        let image = Image {
            dimensions: vec![w, h],
            channels,
            sample_format: format,
            color_space: if channels >= 3 {
                ColorSpace::Rgb
            } else {
                ColorSpace::Gray
            },
            pixel_storage: PixelStorage::Planar,
            // xisf 0.5.1 requires the real formats to declare a representable
            // range on write: unlike the integer formats there is no default to
            // fall back on. The reader never applies bounds to pixel values -
            // it only parses the attribute - so this affects the fixture, not
            // what arcsec does with a real file.
            bounds: if format.requires_bounds() {
                Some(Bounds {
                    low: 0.0,
                    high: 1.0,
                })
            } else {
                None
            },
            id: None,
            uuid: None,
            image_type: None,
            offset: None,
            orientation: None,
        };
        let mut pending = PendingImage::new(image, bytes, BlockOptions::default());
        pending.fits_keywords = keywords;
        let mut writer = Writer::new();
        writer.add_image(pending).expect("add_image");
        let p = std::env::temp_dir().join(name);
        std::fs::write(&p, writer.to_bytes().expect("to_bytes")).unwrap();
        p
    }

    fn kw(name: &str, value: &str) -> FitsKeyword {
        FitsKeyword {
            name: name.to_string(),
            value: value.to_string(),
            comment: String::new(),
        }
    }

    #[test]
    fn reads_a_uint16_greyscale_image() {
        // 4x3, value == index, little-endian u16.
        let mut bytes = Vec::new();
        for i in 0u16..12 {
            bytes.extend_from_slice(&i.to_le_bytes());
        }
        let p = write_xisf(
            "arcsec_x_u16.xisf",
            4,
            3,
            1,
            SampleFormat::UInt16,
            bytes,
            vec![],
        );
        let img = read_xisf_image(&p).expect("read");
        assert_eq!((img.width, img.height), (4, 3));
        assert_eq!(img.data.len(), 12);
        assert_eq!(img.data[0], 0.0);
        assert_eq!(img.data[11], 11.0);
        assert_eq!(read_xisf_dimensions(&p), Some((4, 3)));
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn reads_a_float32_image() {
        let mut bytes = Vec::new();
        for i in 0..6 {
            bytes.extend_from_slice(&(i as f32 * 0.5).to_le_bytes());
        }
        let p = write_xisf(
            "arcsec_x_f32.xisf",
            3,
            2,
            1,
            SampleFormat::Float32,
            bytes,
            vec![],
        );
        let img = read_xisf_image(&p).expect("read");
        assert_eq!(img.data[3], 1.5);
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn colour_channels_are_averaged() {
        // 2x1 RGB, planar: R = [0, 0], G = [3, 30], B = [6, 60].
        // Averages: [3, 30].
        let vals: [u16; 6] = [0, 0, 3, 30, 6, 60];
        let mut bytes = Vec::new();
        for v in vals {
            bytes.extend_from_slice(&v.to_le_bytes());
        }
        let p = write_xisf(
            "arcsec_x_rgb.xisf",
            2,
            1,
            3,
            SampleFormat::UInt16,
            bytes,
            vec![],
        );
        let img = read_xisf_image(&p).expect("read");
        assert_eq!((img.width, img.height), (2, 1));
        assert_eq!(img.data.len(), 2, "colour must collapse to one plane");
        assert!((img.data[0] - 3.0).abs() < 1e-6, "got {}", img.data[0]);
        assert!((img.data[1] - 30.0).abs() < 1e-6, "got {}", img.data[1]);
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn reads_pointing_and_optics_from_fits_keywords() {
        let p = write_xisf(
            "arcsec_x_kw.xisf",
            2,
            2,
            1,
            SampleFormat::UInt16,
            vec![0; 8],
            vec![
                kw("RA", "83.822"),
                kw("DEC", "-5.391"),
                kw("FOCALLEN", "250.0"),
                kw("XPIXSZ", "3.76"),
                kw("XBINNING", "1"),
            ],
        );
        let (ra, dec) = read_xisf_ra_dec(&p).expect("pointing");
        assert!((ra - 83.822).abs() < 1e-6);
        assert!((dec + 5.391).abs() < 1e-6);
        let ps = read_xisf_pixel_scale(&p).expect("scale");
        assert!((ps - 3.1022).abs() < 1e-3, "got {ps}");
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn falls_back_to_crval_when_there_is_no_telescope_pointing() {
        let p = write_xisf(
            "arcsec_x_crval.xisf",
            2,
            2,
            1,
            SampleFormat::UInt16,
            vec![0; 8],
            vec![kw("CRVAL1", "10.5"), kw("CRVAL2", "41.2")],
        );
        assert_eq!(read_xisf_ra_dec(&p), Some((10.5, 41.2)));
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn missing_metadata_is_absent_rather_than_wrong() {
        let p = write_xisf(
            "arcsec_x_bare.xisf",
            2,
            2,
            1,
            SampleFormat::UInt16,
            vec![0; 8],
            vec![],
        );
        assert_eq!(read_xisf_ra_dec(&p), None);
        assert_eq!(read_xisf_pixel_scale(&p), None);
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn quoted_keyword_values_parse() {
        assert_eq!(parse_keyword_value("'83.822'"), Some(83.822));
        assert_eq!(parse_keyword_value("  -5.391 "), Some(-5.391));
        assert_eq!(parse_keyword_value("'M42'"), None);
    }

    #[test]
    fn a_file_with_no_images_is_an_error() {
        let writer = Writer::new();
        let p = std::env::temp_dir().join("arcsec_x_empty.xisf");
        std::fs::write(&p, writer.to_bytes().expect("to_bytes")).unwrap();
        let err = read_xisf_image(&p).expect_err("should refuse");
        assert!(err.contains("no images"), "unhelpful error: {err}");
        std::fs::remove_file(&p).ok();
    }
}
