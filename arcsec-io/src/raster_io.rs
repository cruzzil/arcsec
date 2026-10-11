//! JPEG, PNG and TIFF input, via `zune-jpeg`, `png` and `tiff`.
//!
//! These are the formats of DSLR and phone exports, all-sky and meteor cameras,
//! and stacking tools' "save as" menus. They carry pixels and nothing else arcsec
//! can use: no pointing and no optics keywords, so a solve needs `--ra`/`--spd`
//! on the command line or a blind index, and the pixel scale from `--fov` or the
//! scale search.
//!
//! # Orientation
//!
//! A raster's first row is the top of the picture, where a FITS image's first row
//! is the bottom (FITS viewers draw row 1 at the bottom). ASTAP flips these formats
//! as it loads them, so that the image appears the same way up as it does in any
//! picture viewer and pixel 1,1 of the solution is the picture's bottom-left
//! corner. arcsec does the same, which keeps a `.wcs` for a JPEG interchangeable
//! with ASTAP's. The EXIF orientation tag is ignored, as ASTAP ignores it: the
//! solution describes the pixels as stored.
//!
//! # Values
//!
//! Colour is averaged to greyscale as for XISF, an alpha channel is dropped, and
//! integer samples keep their own scale (an 8-bit image stays 0–255);
//! `ImageBuffer::normalize_for_detection` rescales a range too narrow for the
//! detector's histogram, as it does for float FITS.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use arcsec_core::types::ImageBuffer;

use crate::image_io::{self, ImageFormat, MAX_IMAGE_PIXELS};

/// Load a JPEG, PNG or TIFF as a greyscale buffer, bottom row first.
pub fn read_raster_image(path: &Path, format: ImageFormat) -> Result<ImageBuffer, String> {
    match format {
        ImageFormat::Jpeg => read_jpeg(path),
        ImageFormat::Png => read_png(path),
        ImageFormat::Tiff => read_tiff(path),
        other => Err(format!("{} is not a raster format", other.name())),
    }
}

/// Width and height from the file's header, without decoding the pixels.
pub fn read_raster_dimensions(path: &Path, format: ImageFormat) -> Option<(u32, u32)> {
    let (w, h, _) = raster_header(path, format).ok()?;
    Some((w, h))
}

/// Colour channels before the reduction to grey, alpha not counted.
pub fn read_raster_channels(path: &Path, format: ImageFormat) -> usize {
    raster_header(path, format).map_or(1, |(_, _, c)| c)
}

/// Width, height and colour channels from the header alone.
fn raster_header(path: &Path, format: ImageFormat) -> Result<(u32, u32, usize), String> {
    match format {
        ImageFormat::Jpeg => {
            let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
            let mut dec =
                zune_jpeg::JpegDecoder::new(zune_jpeg::zune_core::bytestream::ZCursor::new(&bytes));
            dec.decode_headers()
                .map_err(|e| format!("JPEG read failed: {e:?}"))?;
            let info = dec.info().ok_or("JPEG read failed: no header")?;
            let colour = dec
                .input_colorspace()
                .map_or(1, |c| if c.num_components() > 1 { 3 } else { 1 });
            Ok((u32::from(info.width), u32::from(info.height), colour))
        }
        ImageFormat::Png => {
            let reader = png_decoder(path)?
                .read_info()
                .map_err(|e| format!("PNG read failed: {e}"))?;
            let info = reader.info();
            Ok((info.width, info.height, png_colour(info.color_type).0))
        }
        ImageFormat::Tiff => {
            let mut dec = tiff_decoder(path)?;
            let (w, h) = dec.dimensions().map_err(tiff_err)?;
            let ct = dec.colortype().map_err(tiff_err)?;
            Ok((w, h, tiff_channels(ct).map_or(1, |(_, colour)| colour)))
        }
        other => Err(format!("{} is not a raster format", other.name())),
    }
}

// ── Shared: interleaved or planar samples to a flipped greyscale image ────────

/// Where sample `c` of pixel `i` lives in a decoded buffer.
#[derive(Clone, Copy)]
enum Layout {
    /// `stride` samples per pixel, side by side.
    Chunky { stride: usize },
    /// Each channel whole, `plane` samples apart.
    Planar { plane: usize },
}

impl Layout {
    #[inline]
    const fn index(self, i: usize, c: usize) -> usize {
        match self {
            Self::Chunky { stride } => i * stride + c,
            Self::Planar { plane } => c * plane + i,
        }
    }
}

/// Average the first `colour` channels of every pixel and turn the image upside
/// down, so that the last row decoded (the bottom of the picture) comes first.
///
/// `sample` reads one value as `f32`; it is generic so that a 16-bit TIFF is
/// converted pixel by pixel rather than copied whole into floats first.
fn to_grey(
    width: usize,
    height: usize,
    colour: usize,
    layout: Layout,
    sample: impl Fn(usize) -> f32,
) -> Result<ImageBuffer, String> {
    let mut data = image_io::try_alloc_pixels(width * height)?;
    let inv = 1.0 / colour as f32;
    for (out_row, row) in data.chunks_exact_mut(width).enumerate() {
        let src_row = height - 1 - out_row;
        for (x, v) in row.iter_mut().enumerate() {
            let i = src_row * width + x;
            *v = if colour == 1 {
                sample(layout.index(i, 0))
            } else {
                (0..colour).map(|c| sample(layout.index(i, c))).sum::<f32>() * inv
            };
        }
    }
    Ok(ImageBuffer {
        data,
        width,
        height,
    })
}

/// Refuse dimensions that are empty or too large before any decoder allocates.
fn check_dims(width: u32, height: u32, samples: usize) -> Result<(usize, usize), String> {
    let (w, h) = (width as usize, height as usize);
    if w == 0 || h == 0 {
        return Err(format!("invalid image dimensions {w}×{h}"));
    }
    image_io::checked_pixel_count(w, h, samples)?;
    Ok((w, h))
}

/// The most bytes a decoder may allocate for one image: [`MAX_IMAGE_PIXELS`]
/// samples of the widest type a raster can hold (`f64`).
const MAX_DECODE_BYTES: usize = MAX_IMAGE_PIXELS.saturating_mul(8);

// ── JPEG ──────────────────────────────────────────────────────────────────────

/// A JPEG is decoded straight to luma: its colour is stored as YCbCr, so the Y
/// plane is the brightness without converting to RGB and back.
fn read_jpeg(path: &Path) -> Result<ImageBuffer, String> {
    use zune_jpeg::zune_core::bytestream::ZCursor;
    use zune_jpeg::zune_core::colorspace::ColorSpace;
    use zune_jpeg::zune_core::options::DecoderOptions;

    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    // The decoder's default refuses images over 16384 px on a side; the pixel
    // count is what bounds memory here.
    let options = DecoderOptions::default()
        .set_max_width(1 << 16)
        .set_max_height(1 << 16)
        .jpeg_set_out_colorspace(ColorSpace::Luma);
    let mut dec = zune_jpeg::JpegDecoder::new_with_options(ZCursor::new(&bytes), options);
    let e = |err| format!("JPEG read failed: {err:?}");
    dec.decode_headers().map_err(e)?;
    let info = dec.info().ok_or("JPEG read failed: no header")?;
    let (w, h) = check_dims(u32::from(info.width), u32::from(info.height), 1)?;
    let pixels = dec.decode().map_err(e)?;
    if pixels.len() < w * h {
        return Err(format!(
            "JPEG decoded {} samples for a {w}×{h} image",
            pixels.len()
        ));
    }
    to_grey(w, h, 1, Layout::Chunky { stride: 1 }, |i| {
        f32::from(pixels[i])
    })
}

// ── PNG ───────────────────────────────────────────────────────────────────────

fn png_decoder(path: &Path) -> Result<png::Decoder<BufReader<File>>, String> {
    let f = File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut dec = png::Decoder::new_with_limits(
        BufReader::new(f),
        png::Limits {
            bytes: MAX_DECODE_BYTES,
        },
    );
    // Palette and sub-byte greyscale come out as whole 8-bit samples; 16-bit
    // samples are left 16-bit.
    dec.set_transformations(png::Transformations::EXPAND);
    dec.set_ignore_text_chunk(true);
    Ok(dec)
}

/// Colour channels to average and samples per pixel.
const fn png_colour(ct: png::ColorType) -> (usize, usize) {
    match ct {
        png::ColorType::Grayscale => (1, 1),
        png::ColorType::GrayscaleAlpha => (1, 2),
        // Indexed is expanded to RGB, or RGBA with a transparency chunk; the
        // reader's output colour type says which.
        png::ColorType::Rgb | png::ColorType::Indexed => (3, 3),
        png::ColorType::Rgba => (3, 4),
    }
}

fn read_png(path: &Path) -> Result<ImageBuffer, String> {
    let e = |err: png::DecodingError| format!("PNG read failed: {err}");
    let mut reader = png_decoder(path)?.read_info().map_err(e)?;
    let (ct, depth) = reader.output_color_type();
    let (colour, stride) = png_colour(ct);
    let info = reader.info();
    let (w, h) = check_dims(info.width, info.height, stride)?;
    let size = reader
        .output_buffer_size()
        .ok_or("PNG read failed: image too large")?;
    let mut buf = Vec::new();
    buf.try_reserve_exact(size)
        .map_err(|_| format!("not enough memory for a {w}×{h} PNG"))?;
    buf.resize(size, 0);
    reader.next_frame(&mut buf).map_err(e)?;
    let layout = Layout::Chunky { stride };
    match depth {
        png::BitDepth::Sixteen => to_grey(w, h, colour, layout, |i| {
            f32::from(u16::from_be_bytes([buf[2 * i], buf[2 * i + 1]]))
        }),
        png::BitDepth::Eight => to_grey(w, h, colour, layout, |i| f32::from(buf[i])),
        other => Err(format!(
            "PNG read failed: {other:?}-bit samples after expansion"
        )),
    }
}

// ── TIFF ──────────────────────────────────────────────────────────────────────

fn tiff_decoder(path: &Path) -> Result<tiff::decoder::Decoder<BufReader<File>>, String> {
    let f = File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut limits = tiff::decoder::Limits::default();
    limits.decoding_buffer_size = MAX_DECODE_BYTES;
    Ok(tiff::decoder::Decoder::new(BufReader::new(f))
        .map_err(tiff_err)?
        .with_limits(limits))
}

fn tiff_err(e: tiff::TiffError) -> String {
    format!("TIFF read failed: {e}")
}

/// Samples per pixel and how many of them are colour (the rest is alpha), or
/// `None` for a colour model with no brightness to average.
const fn tiff_channels(ct: tiff::ColorType) -> Option<(usize, usize)> {
    use tiff::ColorType as C;
    match ct {
        C::Gray(_) => Some((1, 1)),
        C::GrayA(_) => Some((2, 1)),
        C::RGB(_) => Some((3, 3)),
        C::RGBA(_) => Some((4, 3)),
        // Extra bands (a camera's raw channels, a multispectral stack) are all
        // brightness; average them.
        C::Multiband { num_samples, .. } => Some((num_samples as usize, num_samples as usize)),
        _ => None,
    }
}

/// The first image in the file (a TIFF may also hold thumbnails and
/// reduced-resolution copies, which follow it).
fn read_tiff(path: &Path) -> Result<ImageBuffer, String> {
    use tiff::decoder::DecodingResult as D;

    let mut dec = tiff_decoder(path)?;
    let (w, h) = dec.dimensions().map_err(tiff_err)?;
    let ct = dec.colortype().map_err(tiff_err)?;
    let (stride, colour) = tiff_channels(ct).ok_or_else(|| {
        format!("TIFF colour model {ct:?} is not supported for solving (use greyscale or RGB)")
    })?;
    if stride == 0 {
        return Err("TIFF image has no samples per pixel".into());
    }
    let (w, h) = check_dims(w, h, stride)?;

    let mut result = D::U8(Vec::new());
    let layout = dec.read_image_to_buffer(&mut result).map_err(tiff_err)?;
    let npix = w * h;
    let layout = if layout.planes > 1 {
        Layout::Planar { plane: npix }
    } else {
        Layout::Chunky { stride }
    };
    let need = npix * stride;
    macro_rules! grey {
        ($v:expr, $conv:expr) => {{
            let v = $v;
            if v.len() < need {
                return Err(format!(
                    "TIFF decoded {} samples for a {w}×{h}×{stride} image (planar \
                     images larger than the decoder's buffer limit are not supported)",
                    v.len()
                ));
            }
            to_grey(w, h, colour, layout, |i| $conv(v[i]))
        }};
    }
    match &result {
        D::U8(v) => grey!(v, f32::from),
        D::U16(v) => grey!(v, f32::from),
        D::U32(v) => grey!(v, |x: u32| x as f32),
        D::U64(v) => grey!(v, |x: u64| x as f32),
        D::I8(v) => grey!(v, f32::from),
        D::I16(v) => grey!(v, f32::from),
        D::I32(v) => grey!(v, |x: i32| x as f32),
        D::I64(v) => grey!(v, |x: i64| x as f32),
        D::F16(v) => grey!(v, f32::from),
        D::F32(v) => grey!(v, |x: f32| x),
        D::F64(v) => grey!(v, |x: f64| x as f32),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("arcsec_raster_{}_{name}", std::process::id()))
    }

    /// A 3×2 test pattern: value = 10·row + col, row 0 the top of the picture.
    fn pattern(w: usize, h: usize) -> Vec<u16> {
        (0..h)
            .flat_map(|y| (0..w).map(move |x| (10 * y + x) as u16))
            .collect()
    }

    /// Row `y` of the flipped output is row `h-1-y` of the picture.
    fn assert_flipped(img: &ImageBuffer, w: usize, h: usize, scale: f32) {
        assert_eq!((img.width, img.height), (w, h));
        for y in 0..h {
            for x in 0..w {
                let want = (10 * (h - 1 - y) + x) as f32 * scale;
                assert_eq!(img.data[y * w + x], want, "pixel {x},{y}");
            }
        }
    }

    fn write_png(
        path: &Path,
        w: u32,
        h: u32,
        ct: png::ColorType,
        depth: png::BitDepth,
        bytes: &[u8],
    ) {
        let f = File::create(path).unwrap();
        let mut enc = png::Encoder::new(std::io::BufWriter::new(f), w, h);
        enc.set_color(ct);
        enc.set_depth(depth);
        enc.write_header().unwrap().write_image_data(bytes).unwrap();
    }

    #[test]
    fn png_grey16_is_read_bottom_row_first() {
        let p = tmp("g16.png");
        let bytes: Vec<u8> = pattern(3, 2).iter().flat_map(|v| v.to_be_bytes()).collect();
        write_png(
            &p,
            3,
            2,
            png::ColorType::Grayscale,
            png::BitDepth::Sixteen,
            &bytes,
        );
        assert_eq!(image_io::detect_format(&p), Ok(ImageFormat::Png));
        let img = image_io::read_image(&p).unwrap();
        assert_flipped(&img, 3, 2, 1.0);
        assert_eq!(image_io::read_dimensions(&p), Some((3, 2)));
        assert_eq!(image_io::read_channels(&p), 1);
        assert_eq!(image_io::read_ra_dec(&p), None);
        assert_eq!(image_io::read_pixel_scale(&p), None);
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn png_rgba_averages_colour_and_drops_alpha() {
        let p = tmp("rgba.png");
        // Each pixel (v, 2v, 3v, 255): the mean of the colour is 2v.
        let bytes: Vec<u8> = pattern(3, 2)
            .iter()
            .flat_map(|&v| [v as u8, 2 * v as u8, 3 * v as u8, 255])
            .collect();
        write_png(&p, 3, 2, png::ColorType::Rgba, png::BitDepth::Eight, &bytes);
        let img = image_io::read_image(&p).unwrap();
        assert_flipped(&img, 3, 2, 2.0);
        assert_eq!(image_io::read_channels(&p), 3);
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn tiff_grey16_and_rgb_float() {
        use tiff::encoder::{TiffEncoder, colortype};

        let p = tmp("g16.tif");
        let mut enc = TiffEncoder::new(File::create(&p).unwrap()).unwrap();
        enc.write_image::<colortype::Gray16>(3, 2, &pattern(3, 2))
            .unwrap();
        drop(enc);
        assert_eq!(image_io::detect_format(&p), Ok(ImageFormat::Tiff));
        assert_flipped(&image_io::read_image(&p).unwrap(), 3, 2, 1.0);
        assert_eq!(image_io::read_dimensions(&p), Some((3, 2)));
        std::fs::remove_file(&p).ok();

        // 32-bit float RGB, as Siril and PixInsight export it.
        let p = tmp("rgbf.tif");
        let rgb: Vec<f32> = pattern(3, 2)
            .iter()
            .flat_map(|&v| {
                let v = f32::from(v);
                [v, 2.0 * v, 3.0 * v]
            })
            .collect();
        let mut enc = TiffEncoder::new(File::create(&p).unwrap()).unwrap();
        enc.write_image::<colortype::RGB32Float>(3, 2, &rgb)
            .unwrap();
        drop(enc);
        assert_flipped(&image_io::read_image(&p).unwrap(), 3, 2, 2.0);
        assert_eq!(image_io::read_channels(&p), 3);
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn jpeg_is_read_as_luma_bottom_row_first() {
        // Top half dark, bottom half bright: flat 8×8 blocks survive the
        // compression to within a level or two.
        let (w, h) = (16u16, 16u16);
        let grey: Vec<u8> = (0..h)
            .flat_map(|y| (0..w).map(move |_| if y < 8 { 20 } else { 230 }))
            .collect();
        let p = tmp("grey.jpg");
        jpeg_encoder::Encoder::new_file(&p, 95)
            .unwrap()
            .encode(&grey, w, h, jpeg_encoder::ColorType::Luma)
            .unwrap();
        assert_eq!(image_io::detect_format(&p), Ok(ImageFormat::Jpeg));
        let img = image_io::read_image(&p).unwrap();
        assert_eq!((img.width, img.height), (16, 16));
        assert_eq!(image_io::read_dimensions(&p), Some((16, 16)));
        // The bright bottom of the picture is the first row after the flip.
        assert!(
            (img.data[0] - 230.0).abs() < 4.0,
            "first row: {}",
            img.data[0]
        );
        assert!(
            (img.data[15 * 16] - 20.0).abs() < 4.0,
            "last row: {}",
            img.data[15 * 16]
        );
        std::fs::remove_file(&p).ok();

        // Colour is reduced to its luma: a grey RGB picture reads as that grey.
        let rgb: Vec<u8> = grey.iter().flat_map(|&v| [v, v, v]).collect();
        jpeg_encoder::Encoder::new_file(&p, 95)
            .unwrap()
            .encode(&rgb, w, h, jpeg_encoder::ColorType::Rgb)
            .unwrap();
        let img = image_io::read_image(&p).unwrap();
        assert!(
            (img.data[0] - 230.0).abs() < 4.0,
            "first row: {}",
            img.data[0]
        );
        assert_eq!(image_io::read_channels(&p), 3);
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn corrupt_rasters_are_errors() {
        for (name, magic) in [
            ("bad.png", &b"\x89PNG\r\n\x1a\n\0\0\0\rIHDR"[..]),
            ("bad.tif", &b"II*\0\x08\0\0\0\xff\xff"[..]),
            ("bad.jpg", &b"\xff\xd8\xff\xe0\0\x10JFIF"[..]),
        ] {
            let p = tmp(name);
            std::fs::write(&p, magic).unwrap();
            assert!(image_io::read_image(&p).is_err(), "{name}");
            assert_eq!(image_io::read_dimensions(&p), None, "{name}");
            std::fs::remove_file(p).ok();
        }
    }
}
