//! Tests for the enclosing module.

use super::*;

fn q(v: i64) -> Q {
    Q::from_integer(BigInt::from(v))
}

#[test]
fn polynomial_roots_are_exact() {
    // (1 + x)^3 and (x^2 + 3x + 1)^2.
    let cube = vec_poly::pow(&[q(1), q(1)], 3);
    assert_eq!(monic_polynomial_root(&cube, 3), Some(vec![q(1), q(1)]));
    let square = vec_poly::pow(&[q(1), q(3), q(1)], 2);
    assert_eq!(
        monic_polynomial_root(&square, 2),
        Some(vec![q(1), q(3), q(1)])
    );
    assert_eq!(
        monic_polynomial_root(&[q(0), q(0), q(1)], 2),
        Some(vec![q(0), q(1)])
    );
    assert_eq!(monic_polynomial_root(&[q(1), q(0), q(1)], 2), None);
    assert_eq!(monic_polynomial_root(&[q(1), q(1)], 2), None);
    assert_eq!(monic_polynomial_root(&[q(5)], 1), None);
    assert_eq!(monic_polynomial_root(&[q(1)], 3), Some(vec![q(1)]));
}

#[test]
fn root_sides_follow_the_coefficient_signs() {
    assert_eq!(fixed_root_side(&[q(1), q(2), q(1)]), RootSide::Negative);
    assert_eq!(fixed_root_side(&[q(-3)]), RootSide::Negative);
    assert_eq!(fixed_root_side(&[q(0), q(1), q(1)]), RootSide::NonPositive);
    assert_eq!(fixed_root_side(&[q(-1), q(0), q(1)]), RootSide::Other);
}
