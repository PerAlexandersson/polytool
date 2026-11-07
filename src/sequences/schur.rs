//! Rectangular Schur polynomials with one elementary coordinate varied.
//!
//! The canonical generator evaluates the dual Jacobi–Trudi determinant directly,
//! before specializing the moving elementary coordinate. Coefficients are exact
//! integers in ascending degree order; the empty determinant is one.

use crate::linalg::bareiss_leading_principal_minors_polynomial_bigint;
use crate::Polynomial;
use num_bigint::BigInt;
use num_traits::{One, Zero};

/// Generate the rectangular elementary-coordinate deformation through width `max_n`.
///
/// Row `n` is the determinant of the `n` by `n` matrix with entry
/// `e_{height-i+j}(alphabet) + (-1)^(moving-height) t [height-i+j = moving]`,
/// where `i,j` are zero based and elementary functions outside `0..=alphabet.len()`
/// are zero. In particular the substitution occurs before finite-alphabet
/// specialization, so `moving > alphabet.len()` is allowed.
///
/// Returns `None` unless the alphabet is strictly positive, `1 <= height <= d`,
/// and `moving >= 1`, or if exact Bareiss division fails. For valid input the
/// constant terms are positive rectangular Schur specializations, so all pivots
/// are nonzero. This computation is not itself an all-width real-rootedness proof.
/// Complexity is that of one polynomial Bareiss elimination of size `max_n`.
pub fn rectangular_elementary_substitution_polynomials_bigint(
    alphabet: &[BigInt],
    height: usize,
    moving: usize,
    max_n: usize,
) -> Option<Vec<Vec<BigInt>>> {
    if height == 0
        || height > alphabet.len()
        || moving == 0
        || alphabet.iter().any(|a| a <= &BigInt::zero())
    {
        return None;
    }
    let mut elementary = vec![BigInt::zero(); alphabet.len() + 1];
    elementary[0] = BigInt::one();
    for (i, a) in alphabet.iter().enumerate() {
        for k in (1..=i + 1).rev() {
            let term = a * &elementary[k - 1];
            elementary[k] += term;
        }
    }
    let sign = if moving % 2 == height % 2 {
        BigInt::one()
    } else {
        -BigInt::one()
    };
    let matrix: Vec<Vec<_>> = (0..max_n)
        .map(|i| {
            (0..max_n)
                .map(|j| {
                    let index = height.checked_add(j)?.checked_sub(i);
                    let constant = index
                        .and_then(|k| elementary.get(k))
                        .cloned()
                        .unwrap_or_default();
                    Some(Polynomial::new(vec![
                        constant,
                        if index == Some(moving) {
                            sign.clone()
                        } else {
                            BigInt::zero()
                        },
                    ]))
                })
                .collect::<Option<Vec<_>>>()
        })
        .collect::<Option<_>>()?;
    let mut rows = vec![vec![BigInt::one()]];
    rows.extend(
        bareiss_leading_principal_minors_polynomial_bigint(&matrix)?
            .iter()
            .map(|p| p.coeffs().to_vec()),
    );
    Some(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rows(d: usize, r: usize, j: usize, n: usize) -> Vec<Vec<BigInt>> {
        rectangular_elementary_substitution_polynomials_bigint(&vec![1.into(); d], r, j, n).unwrap()
    }

    #[test]
    fn rectangular_examples_and_empty_determinant() {
        let expected = [
            vec![1],
            vec![3, 1],
            vec![6, 6, 1],
            vec![10, 21, 9, 1],
            vec![15, 56, 45, 12, 1],
        ];
        assert_eq!(
            rows(3, 2, 2, 4),
            expected
                .iter()
                .map(|p| p.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>())
                .collect::<Vec<_>>()
        );
        assert_eq!(rows(4, 3, 2, 4)[4], vec![35.into(), 36.into(), 1.into()]);
    }

    #[test]
    fn substitution_above_alphabet_and_invalid_inputs() {
        assert_eq!(rows(2, 1, 3, 3)[3], vec![4.into(), 1.into()]);
        assert!(rectangular_elementary_substitution_polynomials_bigint(&[], 1, 1, 0).is_none());
        assert!(
            rectangular_elementary_substitution_polynomials_bigint(&[0.into()], 1, 1, 0).is_none()
        );
    }

    #[test]
    fn constants_agree_with_rectangular_hook_content() {
        // s_(n,n)(1,1,1) = binomial(n+2,2), independently of moving coordinate.
        for j in 1..=5 {
            for (n, p) in rows(3, 2, j, 8).iter().enumerate() {
                assert_eq!(p[0], BigInt::from((n + 1) * (n + 2) / 2));
            }
        }
    }
}
