//! Positive (semi)definiteness of symmetric `BigInt` matrices.
//!
//! The default [`is_positive_definite`] applies
//! Sylvester's criterion with BigInt Bareiss minors for small matrices and tries
//! the modular CRT minors first from
//! [`MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD`]
//! on.  Semidefiniteness uses exact rational elimination.

use super::determinants::{bareiss_leading_principal_minors_bigint, is_square_bigint_matrix};
use super::modular::modular_leading_principal_minors_bigint;
use super::Q;
use num_bigint::BigInt;
use num_traits::Zero;

/// Matrix dimension at which the default positive-definiteness check switches
/// from BigInt Bareiss to modular CRT reconstruction.
pub const MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD: usize = 30;

/// Check if a symmetric BigInt matrix is positive definite.
///
/// Uses Sylvester's criterion.  Small matrices use fraction-free BigInt
/// Bareiss elimination; matrices of size
/// [`MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD`] or larger try the modular
/// CRT path first and fall back to Bareiss if the CRT bounds do not certify in
/// time. The explicit [`is_positive_definite_bareiss`] and
/// [`is_positive_definite_modular`] entry points are kept for benchmarking.
pub fn is_positive_definite(mat: &[Vec<BigInt>]) -> bool {
    if !is_symmetric_bigint_matrix(mat) {
        return false;
    }
    if mat.len() >= MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD {
        is_positive_definite_modular(mat).unwrap_or_else(|| is_positive_definite_bareiss(mat))
    } else {
        is_positive_definite_bareiss(mat)
    }
}

/// Check positive definiteness via BigInt Bareiss leading principal minors.
pub fn is_positive_definite_bareiss(mat: &[Vec<BigInt>]) -> bool {
    if !is_symmetric_bigint_matrix(mat) {
        return false;
    }
    bareiss_leading_principal_minors_bigint(mat)
        .map(|minors| minors.into_iter().all(|minor| minor > BigInt::zero()))
        .unwrap_or(false)
}

/// Try to check positive definiteness via modular CRT.
///
/// This is a second implementation of Sylvester's criterion.  It reconstructs
/// the leading principal minors from determinants computed over prime fields,
/// using Hadamard bounds to certify when CRT reconstruction is exact.
///
/// Returns `None` when the matrix is not square or when the accumulated CRT
/// modulus does not certify all leading minors within the prime budget.
pub fn is_positive_definite_modular(mat: &[Vec<BigInt>]) -> Option<bool> {
    if !is_square_bigint_matrix(mat) {
        return None;
    }
    if !is_symmetric_bigint_matrix(mat) {
        return Some(false);
    }
    modular_leading_principal_minors_bigint(mat)
        .map(|minors| minors.into_iter().all(|minor| minor > BigInt::zero()))
}

/// Check if a symmetric BigInt matrix is positive semi-definite.
///
/// A symmetric matrix is positive semi-definite iff all diagonal pivots
/// are non-negative, and zero pivots have zero entries in their column below.
pub fn is_positive_semidefinite(mat: &[Vec<BigInt>]) -> bool {
    let n = mat.len();
    if !is_symmetric_bigint_matrix(mat) {
        return false;
    }
    let qmat = bigint_to_q(mat);
    let mut a = qmat;

    for k in 0..n {
        if a[k][k] < Q::zero() {
            return false;
        }
        if a[k][k].is_zero() {
            // All entries below must also be zero
            if a[(k + 1)..n].iter().any(|row| !row[k].is_zero()) {
                return false;
            }
            continue;
        }
        let pivot = a[k][k].clone();
        let pivot_row = a[k].clone();
        for row in &mut a[(k + 1)..n] {
            if row[k].is_zero() {
                continue;
            }
            let factor = row[k].clone() / pivot.clone();
            for (target, source) in row[(k + 1)..n].iter_mut().zip(&pivot_row[(k + 1)..n]) {
                let sub = factor.clone() * source.clone();
                *target -= sub;
            }
            row[k] = Q::zero();
        }
    }
    true
}

fn is_symmetric_bigint_matrix(mat: &[Vec<BigInt>]) -> bool {
    let n = mat.len();
    if !is_square_bigint_matrix(mat) {
        return false;
    }
    (0..n).all(|i| ((i + 1)..n).all(|j| mat[i][j] == mat[j][i]))
}

fn bigint_to_q(mat: &[Vec<BigInt>]) -> Vec<Vec<Q>> {
    mat.iter()
        .map(|row| row.iter().map(|v| Q::from_integer(v.clone())).collect())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::linalg::test_support::bi_mat;
    use num_traits::One;

    fn shifted_hilbert_like_matrix(size: usize) -> Vec<Vec<BigInt>> {
        let mut mat = vec![vec![BigInt::zero(); size]; size];
        for (i, row) in mat.iter_mut().enumerate() {
            for (j, entry) in row.iter_mut().enumerate() {
                *entry = BigInt::from((i + j + 1) as i64);
            }
            row[i] += BigInt::from((size * size) as i64);
        }
        mat
    }

    fn large_positive_diagonal_matrix(size: usize) -> Vec<Vec<BigInt>> {
        let diagonal = BigInt::one() << 700usize;
        let mut mat = vec![vec![BigInt::zero(); size]; size];
        for (i, row) in mat.iter_mut().enumerate() {
            row[i] = diagonal.clone();
        }
        mat
    }

    // -----------------------------------------------------------------------
    // Positive definiteness
    // -----------------------------------------------------------------------

    #[test]
    fn test_pd_identity() {
        let m = bi_mat(&[&[1, 0, 0], &[0, 1, 0], &[0, 0, 1]]);
        assert!(is_positive_definite(&m));
    }

    #[test]
    fn test_pd_2x2() {
        // [[2, 1], [1, 2]] -> pivots: 2, 2 - 1/2 = 3/2 > 0
        let m = bi_mat(&[&[2, 1], &[1, 2]]);
        assert!(is_positive_definite(&m));
        assert!(is_positive_definite_bareiss(&m));
        assert_eq!(is_positive_definite_modular(&m), Some(true));
    }

    #[test]
    fn test_pd_rejects_nonsymmetric_matrix() {
        let m = bi_mat(&[&[1, 2], &[0, 1]]);

        assert!(!is_positive_definite(&m));
        assert!(!is_positive_definite_bareiss(&m));
        assert_eq!(is_positive_definite_modular(&m), Some(false));
    }

    #[test]
    fn test_pd_default_agrees_with_paths_across_threshold() {
        for size in [
            MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD - 1,
            MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD,
        ] {
            let m = shifted_hilbert_like_matrix(size);
            let bareiss = is_positive_definite_bareiss(&m);
            let modular = is_positive_definite_modular(&m).unwrap();
            assert_eq!(bareiss, modular);
            assert_eq!(is_positive_definite(&m), modular);
        }
    }

    #[test]
    fn test_pd_default_falls_back_when_modular_bounds_do_not_certify() {
        let m = large_positive_diagonal_matrix(MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD);

        assert!(modular_leading_principal_minors_bigint(&m).is_none());
        assert_eq!(is_positive_definite_modular(&m), None);
        assert!(is_positive_definite_bareiss(&m));
        assert!(is_positive_definite(&m));
    }

    #[test]
    fn test_pd_negative_definite() {
        // [[-1, 0], [0, -1]] -> first pivot = -1 < 0
        let m = bi_mat(&[&[-1, 0], &[0, -1]]);
        assert!(!is_positive_definite(&m));
        assert_eq!(is_positive_definite_modular(&m), Some(false));
    }

    #[test]
    fn test_pd_zero_matrix() {
        let m = bi_mat(&[&[0, 0], &[0, 0]]);
        assert!(!is_positive_definite(&m)); // PD requires strict positivity
    }

    #[test]
    fn test_pd_indefinite() {
        // [[1, 0], [0, -1]]
        let m = bi_mat(&[&[1, 0], &[0, -1]]);
        assert!(!is_positive_definite(&m));
        assert_eq!(is_positive_definite_modular(&m), Some(false));
    }

    // -----------------------------------------------------------------------
    // Positive semi-definiteness
    // -----------------------------------------------------------------------

    #[test]
    fn test_psd_identity() {
        let m = bi_mat(&[&[1, 0], &[0, 1]]);
        assert!(is_positive_semidefinite(&m));
    }

    #[test]
    fn test_psd_zero_matrix() {
        let m = bi_mat(&[&[0, 0], &[0, 0]]);
        assert!(is_positive_semidefinite(&m));
    }

    #[test]
    fn test_psd_rank_deficient() {
        // [[1, 1], [1, 1]] -> eigenvalues 0, 2 -> PSD
        let m = bi_mat(&[&[1, 1], &[1, 1]]);
        assert!(is_positive_semidefinite(&m));
    }

    #[test]
    fn test_psd_indefinite() {
        // [[1, 2], [2, 1]] -> eigenvalues 3, -1 -> not PSD
        let m = bi_mat(&[&[1, 2], &[2, 1]]);
        assert!(!is_positive_semidefinite(&m));
    }

    #[test]
    fn test_psd_rejects_nonsquare_matrix() {
        let m = bi_mat(&[&[1, 0, 0], &[0, 1, 0]]);

        assert!(!is_positive_semidefinite(&m));
    }

    #[test]
    fn test_psd_rejects_nonsymmetric_matrix() {
        let m = bi_mat(&[&[1, 0], &[2, 1]]);

        assert!(!is_positive_semidefinite(&m));
    }
}
