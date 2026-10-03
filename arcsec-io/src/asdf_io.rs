//! ASDF input, via the `asdf-rs` crate.
//!
//! ASDF is a YAML tree describing the data followed by binary blocks holding it.
//! Unlike FITS and XISF it has no fixed place for "the image" — where the pixels
//! live is a matter of the writing convention. Roman puts them at `roman/data`,
//! astropy's own writer at `data`, and a hand-built file anywhere at all.
//!
//! So the array is located by search rather than by assumption: try the paths the
//! established conventions use, then fall back to walking the tree for the largest
//! two-dimensional array. That is what makes an arbitrary `.asdf` solvable instead
//! of only the ones written by a tool we happened to special-case.

use std::path::Path;

use arcsec_core::types::ImageBuffer;
use asdf::{AsdfFile, Ndarray, Tree, Value};

use crate::image_io;

/// Tree paths that conventionally hold the science image, tried in order.
const IMAGE_PATHS: &[&str] = &[
    "roman/data", // Nancy Grace Roman datamodels
    "data",       // astropy / plain writers
    "sci/data",
    "science/data",
    "primary/data",
];

/// Depth limit for the fallback tree walk.
///
/// Deep enough for the datamodels in the wild, shallow enough that a pathological
/// or cyclic-looking tree cannot spin.
const MAX_WALK_DEPTH: usize = 8;

/// The image array's shape as (width, height).
///
/// ASDF shapes are C-ordered — slowest-varying first — so a 2-D array is
/// `[rows, columns]` and the image is `width = columns`, `height = rows`. That is
/// the opposite of XISF's ordering, and getting it backwards transposes the image
/// silently, so it is stated once here and used everywhere.
fn shape_2d(array: &Ndarray) -> Option<(usize, usize)> {
    let dims: Vec<u64> = array.shape.iter().copied().collect::<Option<Vec<u64>>>()?;
    let [.., rows, cols] = dims.as_slice() else {
        return None;
    };
    let (rows, cols) = (*rows, *cols);
    if rows == 0 || cols == 0 {
        return None;
    }
    Some((cols as usize, rows as usize))
}

/// Number of elements a shape covers, or `None` if it is not at least 2-D.
fn element_count(array: &Ndarray) -> Option<u64> {
    let dims: Vec<u64> = array.shape.iter().copied().collect::<Option<Vec<u64>>>()?;
    if dims.len() < 2 {
        return None;
    }
    dims.iter().try_fold(1u64, |a, d| a.checked_mul(*d))
}

/// Walk the tree for the largest array that could be an image.
///
/// "Largest" rather than "first" because a datamodel carries several arrays — data
/// quality maps, error planes, variance — and the science image is the one with
/// the most pixels among those sharing the largest footprint.
fn find_largest_array(value: Value<'_>, depth: usize, best: &mut Option<(u64, Ndarray)>) {
    if depth > MAX_WALK_DEPTH {
        return;
    }
    if let Some(array) = value.as_ndarray()
        && let Some(n) = element_count(&array)
        && best.as_ref().is_none_or(|(bn, _)| n > *bn)
    {
        *best = Some((n, array));
        return;
    }
    if value.is_mapping() {
        for (_, child) in value.entries() {
            find_largest_array(child, depth + 1, best);
        }
    } else if value.is_sequence() {
        for child in value.items() {
            find_largest_array(child, depth + 1, best);
        }
    }
}

/// Locate the science array in a tree.
fn locate_image(tree: &Tree) -> Result<Ndarray, String> {
    for path in IMAGE_PATHS {
        if let Some(array) = tree.get(path).and_then(|v| v.as_ndarray())
            && shape_2d(&array).is_some()
        {
            return Ok(array);
        }
    }
    let mut best = None;
    if let Some(root) = tree.root() {
        find_largest_array(root, 0, &mut best);
    }
    match best {
        Some((_, array)) if shape_2d(&array).is_some() => Ok(array),
        _ => Err(
            "ASDF file contains no two-dimensional array to solve. Tried the \
             conventional paths (roman/data, data, sci/data) and searched the tree."
                .to_string(),
        ),
    }
}

/// Open a file and its tree together, since every read needs both.
fn open_tree(path: &Path) -> Result<(AsdfFile, Tree), String> {
    let file = AsdfFile::open(path).map_err(|e| format!("ASDF open failed: {e}"))?;
    let tree = file
        .tree()
        .map_err(|e| format!("ASDF tree could not be read: {e}"))?
        .ok_or_else(|| format!("{}: ASDF file has no tree", path.display()))?;
    Ok((file, tree))
}

/// Load an ASDF image as a greyscale buffer.
///
/// Arrays of more than two dimensions — a Roman ramp, say — are read as their
/// first plane, matching what the FITS reader does with a data cube.
pub fn read_asdf_image(path: &Path) -> Result<ImageBuffer, String> {
    let (file, tree) = open_tree(path)?;
    let array = locate_image(&tree)?;
    let (width, height) =
        shape_2d(&array).ok_or_else(|| "ASDF array is not two-dimensional".to_string())?;
    // The whole array is decoded, every plane of a cube included, so it is the
    // array's size that must be bounded before the crate allocates for it.
    let elements = element_count(&array)
        .and_then(|n| usize::try_from(n).ok())
        .unwrap_or(usize::MAX);
    image_io::checked_pixel_count(elements, 1, 1)?;

    let values = file
        .read_array_f64(&array)
        .map_err(|e| format!("ASDF array could not be read as numbers: {e}"))?;

    let npix = width
        .checked_mul(height)
        .ok_or_else(|| format!("ASDF image dimensions {width}×{height} overflow"))?;
    if values.len() < npix {
        return Err(format!(
            "ASDF array declares {width}×{height} = {npix} elements but only {} were read",
            values.len()
        ));
    }

    let data: Vec<f32> = values[..npix].iter().map(|v| *v as f32).collect();
    Ok(ImageBuffer {
        data,
        width,
        height,
    })
}

/// Dimensions without reading the blocks.
pub fn read_asdf_dimensions(path: &Path) -> Option<(u32, u32)> {
    let (_, tree) = open_tree(path).ok()?;
    let array = locate_image(&tree).ok()?;
    let (w, h) = shape_2d(&array)?;
    Some((u32::try_from(w).ok()?, u32::try_from(h).ok()?))
}

/// Read a number from the first tree path that has one.
fn number_at(tree: &Tree, paths: &[&str]) -> Option<f64> {
    paths
        .iter()
        .find_map(|p| tree.get(p).and_then(|v| v.as_f64()))
}

/// Pointing from the tree, in degrees.
///
/// The paths are the ones the Roman and JWST datamodels use. A file without any of
/// them yields nothing, and the caller falls back to `--ra`/`--spd`: guessing a
/// pointing is worse than admitting there isn't one.
pub fn read_asdf_ra_dec(path: &Path) -> Option<(f64, f64)> {
    let (_, tree) = open_tree(path).ok()?;
    let ra = number_at(
        &tree,
        &[
            "roman/meta/wcsinfo/ra_ref",
            "roman/meta/pointing/ra_v1",
            "meta/wcsinfo/ra_ref",
            "meta/pointing/ra_v1",
            "meta/ra",
            "ra",
        ],
    );
    let dec = number_at(
        &tree,
        &[
            "roman/meta/wcsinfo/dec_ref",
            "roman/meta/pointing/dec_v1",
            "meta/wcsinfo/dec_ref",
            "meta/pointing/dec_v1",
            "meta/dec",
            "dec",
        ],
    );
    let crval1 = number_at(&tree, &["meta/wcsinfo/crval1", "wcsinfo/crval1"]);
    let crval2 = number_at(&tree, &["meta/wcsinfo/crval2", "wcsinfo/crval2"]);
    image_io::ra_dec_from(ra, dec, crval1, crval2)
}

/// Plate scale in arcsec/pixel from the tree.
///
/// Prefers a scale the file states outright; otherwise derives it from the optics
/// keywords if the writer carried them across from FITS.
pub fn read_asdf_pixel_scale(path: &Path) -> Option<f64> {
    let (_, tree) = open_tree(path).ok()?;
    if let Some(scale) = number_at(
        &tree,
        &[
            "roman/meta/wcsinfo/pixel_scale",
            "meta/wcsinfo/pixel_scale",
            "meta/pixel_scale",
            "pixel_scale",
        ],
    )
    .filter(|v| v.is_finite() && *v > 0.0)
    {
        return Some(scale);
    }
    image_io::pixel_scale_from(
        number_at(
            &tree,
            &["meta/instrument/focal_length", "meta/focallen", "focallen"],
        ),
        number_at(
            &tree,
            &["meta/instrument/pixel_size", "meta/xpixsz", "xpixsz"],
        ),
        number_at(
            &tree,
            &["meta/instrument/binning", "meta/xbinning", "xbinning"],
        ),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use asdf::AsdfBuilder;

    /// Write an ASDF file with one shaped array at `path_in_tree`.
    fn write_asdf(
        name: &str,
        path_in_tree: &str,
        shape: &[u64],
        values: &[f64],
        scalars: &[(&str, f64)],
    ) -> std::path::PathBuf {
        let mut b = AsdfBuilder::new();
        b.set_array_shaped(path_in_tree, values, shape)
            .expect("set_array_shaped");
        for (k, v) in scalars {
            b.set_f64(k, *v).expect("set_f64");
        }
        let p = std::env::temp_dir().join(name);
        b.write_to_path(&p).expect("write_to_path");
        p
    }

    #[test]
    fn reads_a_two_dimensional_array_at_the_conventional_path() {
        // 3 rows x 4 columns => 4 wide, 3 high.
        let values: Vec<f64> = (0..12).map(|i| i as f64).collect();
        let p = write_asdf("arcsec_a_data.asdf", "data", &[3, 4], &values, &[]);
        let img = read_asdf_image(&p).expect("read");
        assert_eq!(
            (img.width, img.height),
            (4, 3),
            "ASDF shape is [rows, cols]; width must be the last axis"
        );
        assert_eq!(img.data[0], 0.0);
        assert_eq!(img.data[11], 11.0);
        assert_eq!(read_asdf_dimensions(&p), Some((4, 3)));
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn finds_the_roman_convention_path() {
        let values: Vec<f64> = vec![1.0; 6];
        let p = write_asdf("arcsec_a_roman.asdf", "roman/data", &[2, 3], &values, &[]);
        let img = read_asdf_image(&p).expect("read");
        assert_eq!((img.width, img.height), (3, 2));
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn falls_back_to_searching_the_tree() {
        let values: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let p = write_asdf(
            "arcsec_a_odd.asdf",
            "some/unconventional/place/pixels",
            &[4, 5],
            &values,
            &[],
        );
        let img = read_asdf_image(&p).expect("read");
        assert_eq!((img.width, img.height), (5, 4));
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn the_largest_array_wins_the_search() {
        // A datamodel carries several planes; the science image is the biggest.
        let mut b = AsdfBuilder::new();
        b.set_array_shaped("odd/dq", &[0.0f64; 6], &[2, 3]).unwrap();
        b.set_array_shaped("odd/science", &vec![7.0f64; 30], &[5, 6])
            .unwrap();
        b.set_array_shaped("odd/err", &[0.0f64; 12], &[3, 4])
            .unwrap();
        let p = std::env::temp_dir().join("arcsec_a_multi.asdf");
        b.write_to_path(&p).unwrap();
        let img = read_asdf_image(&p).expect("read");
        assert_eq!((img.width, img.height), (6, 5));
        assert_eq!(img.data[0], 7.0, "picked the wrong array");
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn integer_arrays_convert() {
        let mut b = AsdfBuilder::new();
        let values: Vec<u16> = (0..12).collect();
        b.set_array_shaped("data", &values, &[3, 4]).unwrap();
        let p = std::env::temp_dir().join("arcsec_a_u16.asdf");
        b.write_to_path(&p).unwrap();
        let img = read_asdf_image(&p).expect("read");
        assert_eq!(img.data[5], 5.0);
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn reads_pointing_from_the_tree() {
        let values = vec![0.0f64; 4];
        let p = write_asdf(
            "arcsec_a_meta.asdf",
            "data",
            &[2, 2],
            &values,
            &[
                ("meta/wcsinfo/ra_ref", 83.822),
                ("meta/wcsinfo/dec_ref", -5.391),
            ],
        );
        let (ra, dec) = read_asdf_ra_dec(&p).expect("pointing");
        assert!((ra - 83.822).abs() < 1e-9);
        assert!((dec + 5.391).abs() < 1e-9);
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn reads_a_stated_pixel_scale_and_derives_one_otherwise() {
        let v = vec![0.0f64; 4];
        let p = write_asdf(
            "arcsec_a_ps.asdf",
            "data",
            &[2, 2],
            &v,
            &[("meta/wcsinfo/pixel_scale", 0.11)],
        );
        assert_eq!(read_asdf_pixel_scale(&p), Some(0.11));
        std::fs::remove_file(&p).ok();

        let p2 = write_asdf(
            "arcsec_a_optics.asdf",
            "data",
            &[2, 2],
            &v,
            &[
                ("meta/instrument/focal_length", 250.0),
                ("meta/instrument/pixel_size", 3.76),
            ],
        );
        let ps = read_asdf_pixel_scale(&p2).expect("derived scale");
        assert!((ps - 3.1022).abs() < 1e-3, "got {ps}");
        std::fs::remove_file(&p2).ok();
    }

    #[test]
    fn a_file_with_no_array_is_a_clear_error() {
        let mut b = AsdfBuilder::new();
        b.set_str("name", "no pixels here").unwrap();
        let p = std::env::temp_dir().join("arcsec_a_none.asdf");
        b.write_to_path(&p).unwrap();
        let err = read_asdf_image(&p).expect_err("should refuse");
        assert!(err.contains("two-dimensional"), "unhelpful error: {err}");
        assert_eq!(read_asdf_dimensions(&p), None);
        std::fs::remove_file(&p).ok();
    }

    /// A shape no block could fill, declared over a 16-byte block: refused before
    /// the array is decoded, whichever plane would have been used.
    #[test]
    fn an_absurd_shape_is_refused_before_reading() {
        let p = write_asdf(
            "arcsec_a_huge.asdf",
            "data",
            &[2, 2],
            &[1.0, 2.0, 3.0, 4.0],
            &[],
        );
        let good = std::fs::read(&p).unwrap();
        let text = String::from_utf8_lossy(&good).into_owned();
        let shape = text.find("shape: [2, 2]").expect("fixture shape");
        for huge in ["shape: [100000, 100000]", "shape: [65536, 65536, 65536]"] {
            let mut bytes = good[..shape].to_vec();
            bytes.extend_from_slice(huge.as_bytes());
            bytes.extend_from_slice(&good[shape + "shape: [2, 2]".len()..]);
            std::fs::write(&p, &bytes).unwrap();
            let err = read_asdf_image(&p).expect_err(huge);
            assert!(err.contains("limit"), "{huge}: {err}");
        }
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn a_one_dimensional_array_is_not_an_image() {
        let mut b = AsdfBuilder::new();
        b.set_array("data", &[1.0f64; 10]).unwrap();
        let p = std::env::temp_dir().join("arcsec_a_1d.asdf");
        b.write_to_path(&p).unwrap();
        assert!(read_asdf_image(&p).is_err());
        std::fs::remove_file(&p).ok();
    }
}
