//! Fraction-free inversion, rational fallback, shape and singularity contracts.

use exact_algebra::linear_systems::{invert_rational_matrix, MatrixInverseError};
use num_rational::BigRational;

#[test]
fn rational_inverse_and_typed_failures() {
    let matrix = vec![vec![q(1, 2), q(1, 3)], vec![q(2, 3), q(1, 2)]];
    let expected = vec![vec![q(18, 1), q(-12, 1)], vec![q(-24, 1), q(18, 1)]];
    assert_eq!(invert_rational_matrix(&matrix), Ok(expected));
    assert_eq!(invert_rational_matrix(&[]), Ok(vec![]));
    assert_eq!(
        invert_rational_matrix(&[vec![q(1, 1), q(2, 1)]]),
        Err(MatrixInverseError::InvalidDimensions)
    );
    assert_eq!(
        invert_rational_matrix(&[vec![q(1, 1), q(2, 1)], vec![q(2, 1), q(4, 1)]]),
        Err(MatrixInverseError::Singular)
    );
}

#[cfg(feature = "test-hooks")]
#[test]
fn rational_fallback_gives_the_same_inverse() {
    let matrix = vec![vec![q(0, 1), q(2, 3)], vec![q(3, 5), q(1, 7)]];
    let expected = invert_rational_matrix(&matrix).unwrap();
    let _guard = exact_algebra::integer_elimination::ForcedFallback::engage();
    assert_eq!(invert_rational_matrix(&matrix), Ok(expected));
}

fn q(numerator: i64, denominator: i64) -> BigRational {
    BigRational::new(numerator.into(), denominator.into())
}
