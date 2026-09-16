//! Integer matrix utilities for transition matrix computations.
//!
//! These are extracted from combinatoric-core's transition module and
//! provide the linear algebra foundation for basis conversions in Sym, QSym, etc.

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, ToPrimitive, Zero};

/// Identity matrix of size n.
pub fn identity_matrix(n: usize) -> Vec<Vec<i64>> {
    let mut m = vec![vec![0i64; n]; n];
    for i in 0..n {
        m[i][i] = 1;
    }
    m
}

/// Matrix multiplication: `C = A * B`.
///
/// Both inputs must be rectangular and have compatible dimensions. Products
/// are accumulated exactly and this bounded API panics if an output entry does
/// not fit in `i64`.
pub fn mat_mul(a: &[Vec<i64>], b: &[Vec<i64>]) -> Vec<Vec<i64>> {
    let (rows, inner) = matrix_dimensions(a, "left matrix");
    let (b_rows, cols) = matrix_dimensions(b, "right matrix");
    if rows == 0 {
        return vec![];
    }
    assert_eq!(
        inner, b_rows,
        "matrix dimensions do not match for multiplication"
    );
    let mut exact_result = vec![vec![BigInt::zero(); cols]; rows];
    for i in 0..rows {
        for k in 0..inner {
            if a[i][k] == 0 {
                continue;
            }
            for j in 0..cols {
                exact_result[i][j] += BigInt::from(a[i][k]) * BigInt::from(b[k][j]);
            }
        }
    }
    exact_result
        .into_iter()
        .enumerate()
        .map(|(i, row)| {
            row.into_iter()
                .enumerate()
                .map(|(j, entry)| {
                    entry.to_i64().unwrap_or_else(|| {
                        panic!("matrix product entry ({i},{j}) does not fit in i64")
                    })
                })
                .collect()
        })
        .collect()
}

/// Transpose a rectangular matrix.
pub fn transpose(m: &[Vec<i64>]) -> Vec<Vec<i64>> {
    let (rows, cols) = matrix_dimensions(m, "matrix");
    if rows == 0 {
        return vec![];
    }
    let mut t = vec![vec![0i64; rows]; cols];
    for i in 0..rows {
        for j in 0..cols {
            t[j][i] = m[i][j];
        }
    }
    t
}

/// Invert a square integer matrix using exact Gaussian elimination over
/// arbitrary-precision rationals.
///
/// This is a bounded-output API: it panics for a singular matrix, a
/// non-integral inverse, or an integral inverse entry outside `i64`.
pub fn invert_integer_matrix(m: &[Vec<i64>]) -> Vec<Vec<i64>> {
    let n = m.len();
    if n == 0 {
        return vec![];
    }

    let (_, cols) = matrix_dimensions(m, "matrix");
    assert_eq!(n, cols, "matrix inversion requires a square matrix");

    type Q = Ratio<BigInt>;

    let mut aug: Vec<Vec<Q>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut row = Vec::with_capacity(2 * n);
        for j in 0..n {
            row.push(Q::from_integer(BigInt::from(m[i][j])));
        }
        for j in 0..n {
            row.push(if i == j {
                Q::from_integer(BigInt::one())
            } else {
                Q::from_integer(BigInt::zero())
            });
        }
        aug.push(row);
    }

    // Forward elimination
    for col in 0..n {
        let mut pivot = None;
        for row in col..n {
            if !aug[row][col].is_zero() {
                pivot = Some(row);
                break;
            }
        }
        let pivot = pivot.expect("matrix is singular");
        aug.swap(col, pivot);

        let diag = aug[col][col].clone();
        for j in 0..2 * n {
            aug[col][j] = aug[col][j].clone() / diag.clone();
        }

        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = aug[row][col].clone();
            for j in 0..2 * n {
                let val = aug[col][j].clone() * factor.clone();
                aug[row][j] = aug[row][j].clone() - val;
            }
        }
    }

    // Extract inverse (should be integer for our use cases)
    let mut inv = vec![vec![0i64; n]; n];
    for i in 0..n {
        for j in 0..n {
            let val = &aug[i][n + j];
            assert!(
                val.denom().is_one(),
                "inverse matrix entry ({},{}) is not integer: {}",
                i,
                j,
                val
            );
            inv[i][j] = val
                .numer()
                .to_i64()
                .unwrap_or_else(|| panic!("inverse matrix entry ({i},{j}) does not fit in i64"));
        }
    }

    inv
}

fn matrix_dimensions(matrix: &[Vec<i64>], name: &str) -> (usize, usize) {
    let Some(first_row) = matrix.first() else {
        return (0, 0);
    };
    let cols = first_row.len();
    assert!(
        matrix.iter().all(|row| row.len() == cols),
        "{name} rows have inconsistent lengths"
    );
    (matrix.len(), cols)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identity() {
        let id = identity_matrix(3);
        assert_eq!(id[0], vec![1, 0, 0]);
        assert_eq!(id[1], vec![0, 1, 0]);
        assert_eq!(id[2], vec![0, 0, 1]);
    }

    #[test]
    fn test_mat_mul_identity() {
        let id = identity_matrix(2);
        let a = vec![vec![1, 2], vec![3, 4]];
        assert_eq!(mat_mul(&a, &id), a);
        assert_eq!(mat_mul(&id, &a), a);
    }

    #[test]
    fn test_transpose() {
        let m = vec![vec![1, 2, 3], vec![4, 5, 6]];
        let t = transpose(&m);
        assert_eq!(t, vec![vec![1, 4], vec![2, 5], vec![3, 6]]);
    }

    #[test]
    fn test_invert() {
        // [[2, 1], [1, 1]] has inverse [[1, -1], [-1, 2]]
        let m = vec![vec![2, 1], vec![1, 1]];
        let inv = invert_integer_matrix(&m);
        assert_eq!(inv, vec![vec![1, -1], vec![-1, 2]]);

        // Verify: M * M^{-1} = I
        let product = mat_mul(&m, &inv);
        assert_eq!(product, identity_matrix(2));
    }

    #[test]
    fn test_bigint_intermediates_preserve_representable_results() {
        let a = i64::MAX - 1;
        let m = vec![vec![a, a - 1], vec![a + 1, a]];
        let inv = invert_integer_matrix(&m);

        assert_eq!(inv, vec![vec![a, -(a - 1)], vec![-(a + 1), a]]);
        assert_eq!(mat_mul(&m, &inv), identity_matrix(2));
    }

    #[test]
    #[should_panic(expected = "matrix dimensions do not match for multiplication")]
    fn test_mat_mul_rejects_incompatible_dimensions() {
        let _ = mat_mul(&[vec![1, 2]], &[vec![1, 2]]);
    }

    #[test]
    #[should_panic(expected = "matrix rows have inconsistent lengths")]
    fn test_transpose_rejects_ragged_matrix() {
        let _ = transpose(&[vec![1], vec![2, 3]]);
    }

    #[test]
    #[should_panic(expected = "matrix inversion requires a square matrix")]
    fn test_invert_rejects_nonsquare_matrix() {
        let _ = invert_integer_matrix(&[vec![1, 0, 0], vec![0, 1, 0]]);
    }
}
