//! `arcsec_image`: the caller's pixels, read into an [`ImageBuffer`].

use core::ffi::c_void;

use arcsec_core::ImageBuffer;

use crate::error::{Failure, Outcome};
use crate::util::{Versioned, read_versioned};

/// Sample types an `arcsec_image` can hold, in the machine's native byte order.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum arcsec_pixel_type {
    /// `uint8_t`.
    ARCSEC_PIXEL_U8 = 1,
    /// `uint16_t` (Siril's `WORD`, most camera data).
    ARCSEC_PIXEL_U16 = 2,
    /// `int16_t` (FITS BITPIX 16 without BZERO applied).
    ARCSEC_PIXEL_I16 = 3,
    /// `uint32_t`.
    ARCSEC_PIXEL_U32 = 4,
    /// `int32_t`.
    ARCSEC_PIXEL_I32 = 5,
    /// `float`; any range (data normalised to 0..1 is fine).
    ARCSEC_PIXEL_F32 = 6,
    /// `double`.
    ARCSEC_PIXEL_F64 = 7,
}

/// `arcsec_image::flags`: row 0 of the buffer is the *top* of the picture (the
/// last row of a FITS image). Without it, row 0 is FITS row 1, the bottom, as
/// CFITSIO reads it and Siril keeps it.
pub const ARCSEC_IMAGE_TOP_DOWN: u32 = 1;

/// An image in the caller's memory. Not copied by the caller: arcsec reads it
/// during the call and keeps no pointer to it afterwards.
///
/// Zero it, set `struct_size = sizeof(arcsec_image)`, then fill in the fields.
/// Multi-channel images are averaged to one channel before solving.
///
/// The pixel at column `x`, row `y` of channel `c` is read from
/// `plane(c) + y * row_stride + x * sizeof(sample)`, where `plane(c)` is
/// `planes[c]` if `planes` is set, else `data + c * plane_stride`. Strides are in
/// bytes; 0 means tightly packed. Pointers need no particular alignment.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct arcsec_image {
    /// `sizeof(arcsec_image)`.
    pub struct_size: usize,
    /// The first sample of channel 0, when `planes` is NULL.
    pub data: *const c_void,
    /// Optional: `channels` pointers, one per channel plane (Siril's `pdata`).
    /// Takes precedence over `data` and `plane_stride`.
    pub planes: *const *const c_void,
    /// Sample type: one of `arcsec_pixel_type`.
    pub pixel_type: u32,
    /// Columns.
    pub width: u32,
    /// Rows.
    pub height: u32,
    /// Channels (planes); 0 is taken as 1.
    pub channels: u32,
    /// Bytes from one row to the next; 0 = `width * sizeof(sample)`.
    pub row_stride: usize,
    /// Bytes from one plane to the next (ignored with `planes`); 0 =
    /// `height * row_stride`.
    pub plane_stride: usize,
    /// `ARCSEC_IMAGE_*` flags.
    pub flags: u32,
    /// Reserved; set to 0.
    pub reserved: u32,
}

/// Size of the first ABI's `arcsec_image`.
pub(crate) const IMAGE_V1_SIZE: usize = core::mem::size_of::<arcsec_image>();

// SAFETY: repr(C), starts with struct_size, and every field (integers and raw
// pointers) is valid for any bit pattern.
unsafe impl Versioned for arcsec_image {
    fn defaults() -> Self {
        Self {
            struct_size: core::mem::size_of::<Self>(),
            data: core::ptr::null(),
            planes: core::ptr::null(),
            pixel_type: 0,
            width: 0,
            height: 0,
            channels: 1,
            row_stride: 0,
            plane_stride: 0,
            flags: 0,
            reserved: 0,
        }
    }
}

/// Largest side accepted, pixels. Far beyond any sensor; it exists so that a
/// garbage width cannot make the stride arithmetic, or the allocation, absurd.
const MAX_SIDE: u32 = 1 << 20;
/// Largest pixel count accepted (4 Gpx): the f32 copy is 16 GB at this size.
const MAX_PIXELS: usize = 1 << 32;
/// Most channels accepted.
const MAX_CHANNELS: u32 = 64;

/// A sample type: how to read one from memory and widen it.
trait Sample: Copy {
    fn to_f32(self) -> f32;
}

macro_rules! sample {
    ($($t:ty),*) => {$(
        impl Sample for $t {
            #[inline]
            #[allow(clippy::cast_lossless, clippy::cast_precision_loss)]
            fn to_f32(self) -> f32 {
                self as f32
            }
        }
    )*};
}
sample!(u8, u16, i16, u32, i32, f32, f64);

/// Where the image's samples are, once checked.
struct Layout {
    width: usize,
    height: usize,
    row_stride: usize,
    top_down: bool,
}

/// The caller's image, converted to the solver's row-major greyscale `f32`
/// (row 0 = FITS row 1), with its channel count.
///
/// # Safety
///
/// `img` must be NULL or point to an `arcsec_image` (at least `struct_size`
/// readable bytes) whose pointers address every sample the layout describes.
pub(crate) unsafe fn read_image(img: *const arcsec_image) -> Outcome<(ImageBuffer, usize)> {
    // SAFETY: forwarded contract.
    let img = unsafe { read_versioned(img, "arcsec_image", IMAGE_V1_SIZE) }?;
    let elem = match img.pixel_type {
        1 => 1,
        2 | 3 => 2,
        4..=6 => 4,
        7 => 8,
        t => return Err(Failure::invalid(format!("unknown pixel_type {t}"))),
    };
    if img.width == 0 || img.height == 0 || img.width > MAX_SIDE || img.height > MAX_SIDE {
        return Err(Failure::invalid(format!(
            "image size {}x{} is out of range (1..={MAX_SIDE} a side)",
            img.width, img.height
        )));
    }
    let channels = img.channels.max(1);
    if channels > MAX_CHANNELS {
        return Err(Failure::invalid(format!("{channels} channels is too many")));
    }
    let (width, height, channels) = (img.width as usize, img.height as usize, channels as usize);
    if width * height > MAX_PIXELS {
        return Err(Failure::invalid("image has too many pixels"));
    }
    let packed_row = width * elem; // < 2^20 * 8: cannot overflow
    let row_stride = if img.row_stride == 0 {
        packed_row
    } else {
        img.row_stride
    };
    if row_stride < packed_row {
        return Err(Failure::invalid(format!(
            "row_stride {row_stride} is less than a row ({packed_row} bytes)"
        )));
    }
    let plane_bytes = row_stride
        .checked_mul(height)
        .filter(|&b| isize::try_from(b).is_ok())
        .ok_or_else(|| Failure::invalid("row_stride * height overflows"))?;
    let plane_stride = if img.plane_stride == 0 {
        plane_bytes
    } else {
        img.plane_stride
    };
    if img.planes.is_null() {
        if img.data.is_null() {
            return Err(Failure::invalid("arcsec_image has neither data nor planes"));
        }
        if channels > 1 && plane_stride < plane_bytes {
            return Err(Failure::invalid(format!(
                "plane_stride {plane_stride} is less than a plane ({plane_bytes} bytes)"
            )));
        }
        plane_stride
            .checked_mul(channels - 1)
            .and_then(|b| b.checked_add(plane_bytes))
            .filter(|&b| isize::try_from(b).is_ok())
            .ok_or_else(|| Failure::invalid("plane_stride * channels overflows"))?;
    }
    let layout = Layout {
        width,
        height,
        row_stride,
        top_down: img.flags & ARCSEC_IMAGE_TOP_DOWN != 0,
    };

    let mut data: Vec<f32> = Vec::new();
    data.try_reserve_exact(width * height)
        .map_err(|_| Failure::invalid("not enough memory for the image"))?;
    data.resize(width * height, 0.0);
    for c in 0..channels {
        let base: *const u8 = if img.planes.is_null() {
            // SAFETY: offset within the caller's buffer, checked above not to
            // overflow isize; the caller guarantees the planes exist.
            unsafe { img.data.cast::<u8>().add(c * plane_stride) }
        } else {
            // SAFETY: the caller guarantees `channels` plane pointers.
            let p = unsafe { img.planes.add(c).read_unaligned() };
            if p.is_null() {
                return Err(Failure::invalid(format!("planes[{c}] is NULL")));
            }
            p.cast::<u8>()
        };
        // SAFETY: base addresses a plane of `height` rows of `row_stride` bytes,
        // each holding `width` samples of the declared type (caller's contract).
        unsafe {
            match img.pixel_type {
                1 => accumulate::<u8>(&mut data, base, &layout),
                2 => accumulate::<u16>(&mut data, base, &layout),
                3 => accumulate::<i16>(&mut data, base, &layout),
                4 => accumulate::<u32>(&mut data, base, &layout),
                5 => accumulate::<i32>(&mut data, base, &layout),
                6 => accumulate::<f32>(&mut data, base, &layout),
                _ => accumulate::<f64>(&mut data, base, &layout),
            }
        }
    }
    if channels > 1 {
        let inv = 1.0 / channels as f32;
        data.iter_mut().for_each(|v| *v *= inv);
    }
    Ok((
        ImageBuffer {
            data,
            width,
            height,
        },
        channels,
    ))
}

/// Add one plane into `out`, flipping rows if the caller's are top-down.
///
/// # Safety
///
/// `base` must address `layout.height` rows, `layout.row_stride` bytes apart,
/// each holding `layout.width` readable samples of type `T` (any alignment).
unsafe fn accumulate<T: Sample>(out: &mut [f32], base: *const u8, layout: &Layout) {
    for y in 0..layout.height {
        let src_row = if layout.top_down {
            layout.height - 1 - y
        } else {
            y
        };
        // SAFETY: within the plane per the contract.
        let row = unsafe { base.add(src_row * layout.row_stride) }.cast::<T>();
        let dst = &mut out[y * layout.width..(y + 1) * layout.width];
        for (x, d) in dst.iter_mut().enumerate() {
            // SAFETY: sample x of the row is readable; read_unaligned allows any
            // alignment.
            let v = unsafe { row.add(x).read_unaligned() };
            *d += v.to_f32();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn image(data: *const c_void, t: arcsec_pixel_type, w: u32, h: u32) -> arcsec_image {
        arcsec_image {
            data,
            pixel_type: t as u32,
            width: w,
            height: h,
            ..arcsec_image::defaults()
        }
    }

    #[test]
    fn every_sample_type_reads() {
        let u8s = [1u8, 2, 3, 4, 5, 6];
        let (b, _) = unsafe {
            read_image(&image(
                u8s.as_ptr().cast(),
                arcsec_pixel_type::ARCSEC_PIXEL_U8,
                3,
                2,
            ))
        }
        .unwrap();
        assert_eq!(b.data, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert_eq!((b.width, b.height), (3, 2));
        let i16s = [-1i16, 2];
        let (b, _) = unsafe {
            read_image(&image(
                i16s.as_ptr().cast(),
                arcsec_pixel_type::ARCSEC_PIXEL_I16,
                2,
                1,
            ))
        }
        .unwrap();
        assert_eq!(b.data, [-1.0, 2.0]);
        let f64s = [0.25f64, 0.5];
        let (b, _) = unsafe {
            read_image(&image(
                f64s.as_ptr().cast(),
                arcsec_pixel_type::ARCSEC_PIXEL_F64,
                1,
                2,
            ))
        }
        .unwrap();
        assert_eq!(b.data, [0.25, 0.5]);
        let u32s = [7u32, 8];
        let (b, _) = unsafe {
            read_image(&image(
                u32s.as_ptr().cast(),
                arcsec_pixel_type::ARCSEC_PIXEL_U32,
                2,
                1,
            ))
        }
        .unwrap();
        assert_eq!(b.data, [7.0, 8.0]);
    }

    #[test]
    fn strides_planes_flips_and_averaging() {
        // Two 2x2 u16 planes in one buffer, each row padded to 3 samples, planes
        // padded by one row.
        #[rustfmt::skip]
        let buf: [u16; 18] = [
            10, 20, 0,
            30, 40, 0,
            0, 0, 0,
            30, 40, 0,
            50, 60, 0,
            0, 0, 0,
        ];
        let mut img = image(
            buf.as_ptr().cast(),
            arcsec_pixel_type::ARCSEC_PIXEL_U16,
            2,
            2,
        );
        img.channels = 2;
        img.row_stride = 6;
        img.plane_stride = 18;
        let (b, ch) = unsafe { read_image(&raw const img) }.unwrap();
        assert_eq!(ch, 2);
        assert_eq!(b.data, [20.0, 30.0, 40.0, 50.0]);

        img.flags = ARCSEC_IMAGE_TOP_DOWN;
        let (b, _) = unsafe { read_image(&raw const img) }.unwrap();
        assert_eq!(b.data, [40.0, 50.0, 20.0, 30.0]);

        // The same through plane pointers.
        let planes = [buf.as_ptr().cast::<c_void>(), buf[9..].as_ptr().cast()];
        let mut img2 = image(core::ptr::null(), arcsec_pixel_type::ARCSEC_PIXEL_U16, 2, 2);
        img2.planes = planes.as_ptr();
        img2.channels = 2;
        img2.row_stride = 6;
        let (b, _) = unsafe { read_image(&raw const img2) }.unwrap();
        assert_eq!(b.data, [20.0, 30.0, 40.0, 50.0]);
    }

    #[test]
    fn bad_layouts_are_refused() {
        let px = [0f32; 4];
        let ok = image(
            px.as_ptr().cast(),
            arcsec_pixel_type::ARCSEC_PIXEL_F32,
            2,
            2,
        );
        let cases: [(&str, arcsec_image); 8] = [
            (
                "null",
                image(core::ptr::null(), arcsec_pixel_type::ARCSEC_PIXEL_F32, 2, 2),
            ),
            (
                "type",
                arcsec_image {
                    pixel_type: 99,
                    ..ok
                },
            ),
            ("zero width", arcsec_image { width: 0, ..ok }),
            (
                "huge",
                arcsec_image {
                    width: u32::MAX,
                    ..ok
                },
            ),
            (
                "short stride",
                arcsec_image {
                    row_stride: 4,
                    ..ok
                },
            ),
            (
                "size",
                arcsec_image {
                    struct_size: 8,
                    ..ok
                },
            ),
            (
                "plane stride",
                arcsec_image {
                    channels: 2,
                    plane_stride: 4,
                    ..ok
                },
            ),
            (
                "overflow",
                arcsec_image {
                    row_stride: usize::MAX / 2,
                    ..ok
                },
            ),
        ];
        for (name, img) in cases {
            assert!(unsafe { read_image(&raw const img) }.is_err(), "{name}");
        }
        assert!(unsafe { read_image(core::ptr::null()) }.is_err());
        let null_plane = [core::ptr::null::<c_void>()];
        let img = arcsec_image {
            data: core::ptr::null(),
            planes: null_plane.as_ptr(),
            ..ok
        };
        assert!(unsafe { read_image(&raw const img) }.is_err());
    }
}
