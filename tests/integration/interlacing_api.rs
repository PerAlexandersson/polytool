//! Public interlacing API tests.

use num_bigint::BigInt;
use polytool::{check_interlacing, check_interlacing_bigint_coeffs, check_weak_interlacing};

#[test]
fn strict_interlacing_public_api_is_sign_invariant() {
    let f = [-2, 1]; // t - 2
    let g = [3, -4, 1]; // (t - 1)(t - 3)

    assert_eq!(check_interlacing(&f, &g), Some(true));
    assert_eq!(check_interlacing(&[-f[0], -f[1]], &g), Some(true));
    assert_eq!(check_interlacing(&f, &[-g[0], -g[1], -g[2]]), Some(true));
    assert_eq!(
        check_interlacing(&[-f[0], -f[1]], &[-g[0], -g[1], -g[2]]),
        Some(true)
    );
}

#[test]
fn directed_weak_interlacing_distinguishes_same_degree_order() {
    // Roots -1 and 1 respectively, so 1 + t is to the left of -1 + t.
    assert_eq!(check_weak_interlacing(&[1, 1], &[-1, 1]), Some(true));
    assert_eq!(check_weak_interlacing(&[-1, 1], &[1, 1]), Some(false));
}

#[test]
fn weak_interlacing_rejects_nonreal_common_factor() {
    // p = t^2 + 1 and q = t(t^2 + 1) share a non-real common factor.
    assert_eq!(
        check_weak_interlacing(&[1, 0, 1], &[0, 1, 0, 1]),
        Some(false)
    );
}

#[test]
fn invalid_strict_interlacing_degree_relation_returns_none() {
    assert_eq!(check_interlacing(&[1, 0, 1], &[1, 1]), None);
}

#[test]
fn bigint_interlacing_handles_coefficients_beyond_i64() {
    let center = BigInt::from(10).pow(20);
    let f = vec![-&center, BigInt::from(1)];
    let g = vec![
        (&center - 1u32) * (&center + 1u32),
        -BigInt::from(2) * &center,
        BigInt::from(1),
    ];

    assert_eq!(check_interlacing_bigint_coeffs(&f, &g), Some(true));
}
