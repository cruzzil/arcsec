//! GIVENS rotation-based least-squares solver.
//! Reference: Montenbruck & Pfleger, "Astronomy on the Personal Computer"

use crate::error::{ArcsecError, Result};
use crate::types::PlateConstants;

/// Solves the overdetermined system A·x = b with 3 unknowns via GIVENS rotations.
///
/// `a_matrix`: column-major, indexed `[col][row]`.
///   - `a_matrix[0]` = image x-coordinates
///   - `a_matrix[1]` = image y-coordinates
///   - `a_matrix[2]` = constant 1.0 for each equation
///
/// `b_matrix`: right-hand side (reference positions for one axis).
///
/// Returns the solution vector `[coeff_x, coeff_y, coeff_const]`.
///
/// # Errors
///
/// [`ArcsecError::Singular`] if the system is degenerate, has fewer equations than
/// unknowns, or the columns and `b_matrix` differ in length.
pub fn lsq_fit(a_matrix: &[Vec<f64>], b_matrix: &[f64]) -> Result<Vec<f64>> {
    const TINY: f64 = 1e-10;

    let nr_columns = a_matrix.len();
    let nr_equations = b_matrix.len();
    // Fewer equations than unknowns cannot be solved, and the elimination below
    // indexes row j of column j, so it would also run off the end of the columns.
    if nr_columns == 0
        || nr_equations < nr_columns
        || a_matrix.iter().any(|col| col.len() != nr_equations)
    {
        return Err(ArcsecError::Singular);
    }

    // Duplicate matrices so the caller's originals are not modified
    let mut temp: Vec<Vec<f64>> = a_matrix.to_vec();
    let mut b: Vec<f64> = b_matrix.to_vec();

    // Forward elimination via GIVENS rotations
    for j in 0..nr_columns {
        for i in (j + 1)..nr_equations {
            if temp[j][i] == 0.0 {
                continue;
            }
            let (p, q) = if temp[j][j].abs() < TINY * temp[j][i].abs() {
                // Near-zero pivot: swap rows j and i with sign change
                let old_ji = temp[j][i];
                temp[j][j] = -old_ji;
                temp[j][i] = 0.0;
                (0.0_f64, 1.0_f64)
            } else {
                let mut h = (temp[j][j] * temp[j][j] + temp[j][i] * temp[j][i]).sqrt();
                if temp[j][j] < 0.0 {
                    h = -h;
                }
                let p = temp[j][j] / h;
                let q = -temp[j][i] / h;
                temp[j][j] = h;
                temp[j][i] = 0.0;
                (p, q)
            };

            for col in temp.iter_mut().skip(j + 1) {
                let h = p * col[j] - q * col[i];
                col[i] = q * col[j] + p * col[i];
                col[j] = h;
            }
            let h = p * b[j] - q * b[i];
            b[i] = q * b[j] + p * b[i];
            b[j] = h;
        }
    }

    // Back substitution
    let mut x = vec![0.0f64; nr_columns];
    for i in (0..nr_columns).rev() {
        let mut h = b[i];
        for k in (i + 1)..nr_columns {
            h -= temp[k][i] * x[k];
        }
        if temp[i][i].abs() <= 1e-30 {
            return Err(ArcsecError::Singular);
        }
        x[i] = h / temp[i][i];
    }
    Ok(x)
}

/// Solves for all 6 plate constants by calling `lsq_fit` for each axis.
/// Validates that the X and Y pixel scales are within 10% of each other.
///
/// `img_xy[i]` and `ref_xy[i]` must be the same star; at least three pairs are needed.
///
/// # Errors
///
/// - [`ArcsecError::Singular`] if the slices differ in length, hold fewer than three
///   pairs, or the points are degenerate (e.g. collinear).
/// - [`ArcsecError::BadSolution`] if the X and Y scales disagree by more than 10%.
pub fn solve_plate_constants(
    img_xy: &[(f64, f64)],
    ref_xy: &[(f64, f64)],
) -> Result<PlateConstants> {
    let n = img_xy.len();
    if ref_xy.len() != n {
        return Err(ArcsecError::Singular);
    }

    // Build column-major A matrix: [x_pixels, y_pixels, 1.0]
    let col_x: Vec<f64> = img_xy.iter().map(|&(x, _)| x).collect();
    let col_y: Vec<f64> = img_xy.iter().map(|&(_, y)| y).collect();
    let col_ones: Vec<f64> = vec![1.0; n];
    let a_matrix = vec![col_x, col_y, col_ones];

    let b_x: Vec<f64> = ref_xy.iter().map(|&(x, _)| x).collect();
    let b_y: Vec<f64> = ref_xy.iter().map(|&(_, y)| y).collect();

    let sol_x = lsq_fit(&a_matrix, &b_x)?;
    let sol_y = lsq_fit(&a_matrix, &b_y)?;

    // Check that X and Y pixel scales agree within 10% (a larger disagreement means the fit is not a similarity transform)
    let xy_sqr_ratio =
        (sol_x[0].powi(2) + sol_x[1].powi(2)) / (1e-8 + sol_y[0].powi(2) + sol_y[1].powi(2));
    if !(0.9..=1.1).contains(&xy_sqr_ratio) {
        return Err(ArcsecError::BadSolution {
            ratio: xy_sqr_ratio,
        });
    }

    Ok(PlateConstants {
        a: sol_x[0],
        b: sol_x[1],
        c: sol_x[2],
        d: sol_y[0],
        e: sol_y[1],
        f: sol_y[2],
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(a: f64, b: f64, tol: f64) {
        assert!(
            (a - b).abs() < tol,
            "expected {b:.6} got {a:.6} (diff {:.2e})",
            (a - b).abs()
        );
    }

    /// 2D star grid: 4×4 = 16 stars, well-distributed (not collinear)
    fn grid_stars() -> Vec<(f64, f64)> {
        (0..4)
            .flat_map(|i| (0..4).map(move |j| (i as f64 * 100.0 + 50.0, j as f64 * 100.0 + 50.0)))
            .collect()
    }

    /// Grid of stars: ref == image → identity transform [1, 0, 0, 0, 1, 0]
    #[test]
    fn identity_transform() {
        let stars = grid_stars();
        let pc = solve_plate_constants(&stars, &stars).unwrap();
        assert_close(pc.a, 1.0, 1e-10);
        assert_close(pc.b, 0.0, 1e-10);
        assert_close(pc.c, 0.0, 1e-10);
        assert_close(pc.d, 0.0, 1e-10);
        assert_close(pc.e, 1.0, 1e-10);
        assert_close(pc.f, 0.0, 1e-10);
    }

    /// Pure translation: ref = image + (100, 200)
    #[test]
    fn pure_translation() {
        let img = grid_stars();
        let rf: Vec<(f64, f64)> = img.iter().map(|&(x, y)| (x + 100.0, y + 200.0)).collect();
        let pc = solve_plate_constants(&img, &rf).unwrap();
        assert_close(pc.a, 1.0, 1e-10);
        assert_close(pc.b, 0.0, 1e-10);
        assert_close(pc.c, 100.0, 1e-8);
        assert_close(pc.d, 0.0, 1e-10);
        assert_close(pc.e, 1.0, 1e-10);
        assert_close(pc.f, 200.0, 1e-8);
    }

    /// Pure scale 2×: ref = 2 * image
    #[test]
    fn pure_scale() {
        let img = grid_stars();
        let rf: Vec<(f64, f64)> = img.iter().map(|&(x, y)| (2.0 * x, 2.0 * y)).collect();
        let pc = solve_plate_constants(&img, &rf).unwrap();
        assert_close(pc.a, 2.0, 1e-8);
        assert_close(pc.b, 0.0, 1e-8);
        assert_close(pc.d, 0.0, 1e-8);
        assert_close(pc.e, 2.0, 1e-8);
    }

    /// 90° rotation: `ref_x = -img_y`, `ref_y = img_x`
    #[test]
    fn rotation_90() {
        let img = grid_stars();
        let rf: Vec<(f64, f64)> = img.iter().map(|&(x, y)| (-y, x)).collect();
        let pc = solve_plate_constants(&img, &rf).unwrap();
        assert_close(pc.a, 0.0, 1e-8);
        assert_close(pc.b, -1.0, 1e-8);
        assert_close(pc.d, 1.0, 1e-8);
        assert_close(pc.e, 0.0, 1e-8);
    }

    /// Degenerate input (all stars at same point) → Singular error
    #[test]
    fn singular_input() {
        let img: Vec<(f64, f64)> = vec![(1.0, 1.0); 5];
        let rf: Vec<(f64, f64)> = vec![(2.0, 2.0); 5];
        assert!(matches!(
            solve_plate_constants(&img, &rf),
            Err(ArcsecError::Singular)
        ));
    }

    /// Bad scale ratio → `BadSolution` error
    #[test]
    fn bad_solution_ratio() {
        // X scale = 1, Y scale = 10 → ratio = 0.01 → out of [0.9, 1.1]
        let img: Vec<(f64, f64)> = (0..10)
            .map(|i| (i as f64 * 50.0 + 1.0, i as f64 * 3.0 + 1.0))
            .collect();
        let rf: Vec<(f64, f64)> = img.iter().map(|&(x, y)| (x, y * 10.0)).collect();
        assert!(matches!(
            solve_plate_constants(&img, &rf),
            Err(ArcsecError::BadSolution { .. })
        ));
    }

    /// `lsq_fit` directly: simple 1-unknown system ax = b → x = b/a
    #[test]
    fn lsq_fit_1d() {
        // Overdetermined: 5 equations, 1 unknown (simplify: constant 1, no x/y)
        let a = vec![vec![2.0f64, 2.0, 2.0, 2.0, 2.0]];
        let b = vec![6.0, 6.0, 6.0, 6.0, 6.0];
        let x = lsq_fit(&a, &b).unwrap();
        assert_close(x[0], 3.0, 1e-10);
    }

    /// Fewer than three points, or mismatched slices, must be an error, not a panic.
    #[test]
    fn too_few_or_mismatched_points_are_errors() {
        for n in 0..3 {
            let img: Vec<(f64, f64)> = (0..n).map(|i| (i as f64, 2.0 * i as f64)).collect();
            assert!(
                matches!(
                    solve_plate_constants(&img, &img),
                    Err(ArcsecError::Singular)
                ),
                "n = {n}"
            );
        }
        let img = grid_stars();
        assert!(matches!(
            solve_plate_constants(&img, &img[1..]),
            Err(ArcsecError::Singular)
        ));
        assert!(matches!(lsq_fit(&[], &[]), Err(ArcsecError::Singular)));
    }

    /// A known similarity (with a flip) plus Gaussian noise: the fit recovers the
    /// transform to within what the noise allows, and its residuals match the
    /// noise level.
    #[test]
    fn recovers_a_noisy_similarity_transform() {
        let mut rng = crate::test_support::Rng::new(4);
        let (s, r) = (2.37_f64, -0.83_f64);
        let truth = [
            -s * r.cos(),
            s * r.sin(),
            1234.5,
            s * r.sin(),
            s * r.cos(),
            -987.6,
        ];
        let sigma = 0.5;
        let img: Vec<(f64, f64)> = (0..400)
            .map(|_| (rng.range(0.0, 4000.0), rng.range(0.0, 3000.0)))
            .collect();
        let cat: Vec<(f64, f64)> = img
            .iter()
            .map(|&(x, y)| {
                (
                    truth[0] * x + truth[1] * y + truth[2] + sigma * rng.gauss(),
                    truth[3] * x + truth[4] * y + truth[5] + sigma * rng.gauss(),
                )
            })
            .collect();
        let p = solve_plate_constants(&img, &cat).unwrap();
        let got = [p.a, p.b, p.c, p.d, p.e, p.f];
        // Slope errors ~ σ / (√n · spread) ≈ 1e-5; offsets ~ σ·few/√n ≈ 0.1.
        for k in [0, 1, 3, 4] {
            assert_close(got[k], truth[k], 1e-4);
        }
        assert_close(p.c, truth[2], 0.3);
        assert_close(p.f, truth[5], 0.3);
        let rms = (img
            .iter()
            .zip(&cat)
            .map(|(&(x, y), &(u, v))| {
                (p.a * x + p.b * y + p.c - u).powi(2) + (p.d * x + p.e * y + p.f - v).powi(2)
            })
            .sum::<f64>()
            / img.len() as f64
            / 2.0)
            .sqrt();
        assert!((rms - sigma).abs() < 0.05, "per-axis rms {rms}");
    }

    /// Exactly three non-collinear points determine the transform exactly.
    #[test]
    fn three_points_are_enough() {
        let img = [(0.0, 0.0), (10.0, 0.0), (0.0, 10.0)];
        let cat = [(5.0, 5.0), (5.0, 25.0), (-15.0, 5.0)]; // 90° rotation, scale 2
        let p = solve_plate_constants(&img, &cat).unwrap();
        assert_close(p.a, 0.0, 1e-12);
        assert_close(p.b, -2.0, 1e-12);
        assert_close(p.d, 2.0, 1e-12);
        assert_close(p.e, 0.0, 1e-12);
        assert_close(p.c, 5.0, 1e-12);
        assert_close(p.f, 5.0, 1e-12);
        // Collinear points cannot.
        let line = [(0.0, 0.0), (1.0, 1.0), (2.0, 2.0), (3.0, 3.0)];
        assert!(solve_plate_constants(&line, &line).is_err());
    }
}
