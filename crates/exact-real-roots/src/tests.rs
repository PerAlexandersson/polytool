//! Regression tests and independent exact root-counting oracles.

use crate::arithmetic::{
    degree, poly_derivative, poly_mul, primitive_positive, pseudo_remainder_positive_multiplier,
    pseudo_remainder_positive_multiplier_into, strip_initial_zeros_bigint, trim, trim_slice,
};
use crate::coefficient_tests::low_degree_is_real_rooted;
use crate::prs::{
    count_positive_roots_from_sturm_chain, derivative_sturm_chain, derivative_sturm_summary,
    signed_remainder_sequence, signed_remainder_sequence_subresultant, sturm_chain_from_squarefree,
    variations_at_neg_infinity, variations_at_pos_infinity, variations_at_zero_plus,
};
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

use super::*;
use num_rational::BigRational;

#[test]
fn newton_inequalities_hold_for_real_rooted_polynomials_of_any_sign() {
    let big = |v: &[i64]| v.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
    // (t - 1)(t - 2), -(1 + t)^2 and (t + 1)(t - 3)(t - 4) are real-rooted.
    for p in [vec![2, -3, 1], vec![-1, -2, -1], vec![12, 5, -6, 1]] {
        assert!(
            satisfies_newton_inequalities_bigint_coeffs(&big(&p)),
            "{p:?}"
        );
    }
    // 1 - t + t^2 has complex roots and violates the inequality at k = 1.
    assert!(!satisfies_newton_inequalities_bigint_coeffs(&big(&[
        1, -1, 1
    ])));
}

fn b(v: &[i64]) -> Vec<BigInt> {
    v.iter().map(|&x| BigInt::from(x)).collect()
}

/// Independently compute `|lc(b)|^(deg(a)-deg(b)+1) * rem(a, b)` over Q.
///
/// This intentionally does not use pseudo-division, so it can detect a
/// regression in the shared integer pseudo-remainder kernel.
fn rational_scaled_remainder(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let a = trim_slice(a);
    let b = trim_slice(b);
    let Some(a_degree) = degree(&a) else {
        return Vec::new();
    };
    let b_degree = degree(&b).expect("rational remainder divisor must be nonzero");
    if a_degree < b_degree {
        return a;
    }

    let divisor_sign = if b[b_degree].is_negative() {
        -BigInt::one()
    } else {
        BigInt::one()
    };
    let divisor: Vec<BigRational> = b
        .iter()
        .map(|coefficient| BigRational::from_integer(coefficient * &divisor_sign))
        .collect();
    let mut remainder: Vec<BigRational> =
        a.iter().cloned().map(BigRational::from_integer).collect();

    while let Some(remainder_degree) = remainder.iter().rposition(|c| !c.is_zero()) {
        if remainder_degree < b_degree {
            break;
        }
        let shift = remainder_degree - b_degree;
        let quotient_term = remainder[remainder_degree].clone() / &divisor[b_degree];
        for (j, coefficient) in divisor.iter().enumerate() {
            remainder[j + shift] -= &quotient_term * coefficient;
        }
        while remainder.last().is_some_and(BigRational::is_zero) {
            remainder.pop();
        }
    }

    let scale = b[b_degree].abs().pow((a_degree - b_degree + 1) as u32);
    trim(
        remainder
            .into_iter()
            .map(|coefficient| {
                let scaled = coefficient * BigRational::from_integer(scale.clone());
                assert!(scaled.denom().is_one());
                scaled.numer().clone()
            })
            .collect(),
    )
}

fn mul_i64(a: &[i64], b: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; a.len() + b.len() - 1];
    for (i, &ca) in a.iter().enumerate() {
        for (j, &cb) in b.iter().enumerate() {
            r[i + j] += ca * cb;
        }
    }
    while r.last() == Some(&0) {
        r.pop();
    }
    r
}

fn binomial(n: usize, k: usize) -> BigInt {
    let k = k.min(n - k);
    (1..=k).fold(BigInt::one(), |value, j| {
        value * BigInt::from(n - k + j) / BigInt::from(j)
    })
}

fn generated_eulerian(n: usize) -> Vec<BigInt> {
    let mut row = vec![BigInt::one()];
    for current_n in 2..=n {
        let mut next = vec![BigInt::zero(); current_n];
        for k in 0..current_n {
            if k < row.len() {
                next[k] += &row[k] * BigInt::from(k + 1);
            }
            if k > 0 {
                next[k] += &row[k - 1] * BigInt::from(current_n - k);
            }
        }
        row = next;
    }
    row
}

fn generated_type_b_eulerian(n: usize) -> Vec<BigInt> {
    let mut row = vec![BigInt::one()];
    for current_n in 1..=n {
        let mut next = vec![BigInt::zero(); current_n + 1];
        for k in 0..=current_n {
            if k < row.len() {
                next[k] += &row[k] * BigInt::from(2 * k + 1);
            }
            if k > 0 {
                next[k] += &row[k - 1] * BigInt::from(2 * current_n - 2 * k + 1);
            }
        }
        row = next;
    }
    row
}

fn generated_narayana(n: usize) -> Vec<BigInt> {
    (0..n)
        .map(|k| binomial(n, k + 1) * binomial(n, k) / BigInt::from(n))
        .collect()
}

#[test]
fn test_pseudo_remainder_basic() {
    // (x^3 - 1) mod (x - 1) = 0.
    let r = pseudo_remainder_positive_multiplier(&b(&[-1, 0, 0, 1]), &b(&[-1, 1]));
    assert!(r.is_empty());

    // (x^2 + 1) mod (x - 2) = 5.  The divisor is monic, so prem = rem.
    let r = pseudo_remainder_positive_multiplier(&b(&[1, 0, 1]), &b(&[-2, 1]));
    assert_eq!(r, b(&[5]));
}

#[test]
fn test_pseudo_remainder_matches_independent_rational_oracle() {
    let mut cases = vec![
        // Negative-leading divisor.
        (b(&[7, -3, 5, 0, 2]), b(&[2, 0, -3])),
        // The exact remainder is constant, forcing a sparse degree drop.
        (b(&[-3, 0, -2, 15, 0, 3]), b(&[-2, 0, 0, 3])),
        (b(&[-3, 0, -2, 15, 0, 3]), b(&[2, 0, 0, -3])),
        // Dividend degree below divisor degree, with trailing zeros.
        (b(&[5, -2, 0, 0]), b(&[1, 0, 4, 0])),
        // A constant divisor must leave the zero remainder.
        (b(&[9, -4, 7, 1]), b(&[-3])),
    ];

    let mut state = 0xa54f_f53a_5f1d_36f1u64;
    for case_index in 0..512usize {
        let a_degree = case_index % 8;
        let b_degree = 1 + (case_index * 5 % 6);
        let mut a = Vec::with_capacity(a_degree + 1);
        let mut divisor = Vec::with_capacity(b_degree + 1);
        for _ in 0..=a_degree {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            a.push(BigInt::from(((state >> 32) % 17) as i64 - 8));
        }
        for _ in 0..=b_degree {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            divisor.push(BigInt::from(((state >> 32) % 17) as i64 - 8));
        }
        if a[a_degree].is_zero() {
            a[a_degree] = BigInt::one();
        }
        if divisor[b_degree].is_zero() {
            divisor[b_degree] = BigInt::one();
        }
        if case_index % 2 == 1 && divisor[b_degree].is_positive() {
            for coefficient in &mut divisor {
                *coefficient = -std::mem::take(coefficient);
            }
        }
        cases.push((a, divisor));
    }

    for (a, divisor) in cases {
        let expected = rational_scaled_remainder(&a, &divisor);
        assert_eq!(
            pseudo_remainder_positive_multiplier(&a, &divisor),
            expected,
            "a={a:?}, divisor={divisor:?}"
        );

        let mut scratch = vec![(BigInt::one() << 256usize) - 1; a.len() + 11];
        pseudo_remainder_positive_multiplier_into(&mut scratch, &a, &divisor);
        assert_eq!(scratch, expected, "a={a:?}, divisor={divisor:?}");
    }
}

#[test]
fn test_squarefree_part_and_degree() {
    // (1+x)^3 has square-free part 1+x.
    assert_eq!(squarefree_degree_bigint_coeffs(&b(&[1, 3, 3, 1])), 1);

    // (x-1)^2 (x-2) = -2 + 5x - 4x^2 + x^3, square-free degree 2.
    assert_eq!(squarefree_degree_bigint_coeffs(&b(&[-2, 5, -4, 1])), 2);
}

#[test]
fn test_count_real_roots() {
    assert_eq!(count_real_roots_prs_bigint_coeffs(&b(&[1, 0, 1])), 0);
    assert_eq!(count_real_roots_prs_bigint_coeffs(&b(&[-6, 11, -6, 1])), 3);
    assert_eq!(count_real_roots_prs_bigint_coeffs(&b(&[1, 3, 3, 1])), 1);
    assert_eq!(count_real_roots_prs_bigint_coeffs(&b(&[0, -1, 0, 1])), 3);
}

#[test]
fn test_a046741_degree_eight_cross_checks_exact_paths() {
    let polynomial = b(&[1, 22, 191, 838, 1982, 2504, 1577, 424, 34]);
    let prs_count = count_real_roots_prs_bigint_coeffs(&polynomial);
    let uspensky_count = count_real_roots_uspensky_bigint_coeffs(&polynomial);
    let derivative = poly_derivative(&polynomial);
    let primitive_chain = signed_remainder_sequence(&polynomial, &derivative);
    let primitive_count = variations_at_neg_infinity(&primitive_chain)
        .saturating_sub(variations_at_pos_infinity(&primitive_chain));

    assert_eq!(prs_count, 8);
    assert_eq!(uspensky_count, 8);
    assert_eq!(primitive_count, 8);
    assert!(satisfies_newton_inequalities_bigint_coeffs(&polynomial));
    assert_eq!(
        is_real_rooted_one_signed_bigint_coeffs(&polynomial),
        Some(true)
    );
    assert!(is_real_rooted_fast_bigint_coeffs(&polynomial));
}

#[test]
fn test_count_positive_roots() {
    // (x-1)(x-2)(x+3) = 6 - 7x + x^3 has two positive roots.
    assert_eq!(
        count_positive_roots_prs_bigint_coeffs(&b(&[6, -7, 0, 1])),
        2
    );
    assert_eq!(count_positive_roots_prs_bigint_coeffs(&b(&[1, 2, 1])), 0);
}

#[test]
fn test_tarski_query_constant_is_real_root_count() {
    let p = b(&[-6, 11, -6, 1]); // (x-1)(x-2)(x-3)
    assert_eq!(tarski_query_prs_bigint_coeffs(&p, &b(&[1])), 3);
    assert_eq!(tarski_query_prs_bigint_coeffs(&p, &b(&[-1])), -3);

    let q = b(&[1, 0, 1]); // x^2+1 has no real roots.
    assert_eq!(tarski_query_prs_bigint_coeffs(&q, &b(&[1])), 0);
}

#[test]
fn test_root_sign_counts_basic() {
    let p = b(&[-6, 11, -6, 1]); // roots 1,2,3

    assert_eq!(
        count_root_signs_prs_bigint_coeffs(&p, &b(&[-2, 1])), // x-2
        RootSignCounts {
            positive: 1,
            negative: 1,
            zero: 1,
            total: 3,
        }
    );
    assert_eq!(
        count_root_signs_prs_bigint_coeffs(&p, &b(&[-4, 1])), // x-4
        RootSignCounts {
            positive: 0,
            negative: 3,
            zero: 0,
            total: 3,
        }
    );
    assert_eq!(
        count_root_signs_prs_bigint_coeffs(&p, &b(&[2, -3, 1])), // (x-1)(x-2)
        RootSignCounts {
            positive: 1,
            negative: 0,
            zero: 2,
            total: 3,
        }
    );
}

#[test]
fn test_root_sign_counts_squarefree_root_part() {
    let repeated = b(&[4, -12, 13, -6, 1]); // (x-1)^2 (x-2)^2
    assert_eq!(
        count_root_signs_prs_bigint_coeffs(&repeated, &b(&[-1, 1])), // x-1
        RootSignCounts {
            positive: 1,
            negative: 0,
            zero: 1,
            total: 2,
        }
    );
}

#[test]
fn test_one_signed_real_rootedness() {
    assert_eq!(
        is_real_rooted_one_signed_bigint_coeffs(&b(&[1, 4, 6, 4, 1])),
        Some(true)
    );
    assert_eq!(
        is_real_rooted_one_signed_bigint_coeffs(&b(&[1, 43, 196, 168, 23, 1])),
        Some(false)
    );
    assert_eq!(
        is_real_rooted_one_signed_bigint_coeffs(&b(&[1, -2, 1])),
        None
    );
}

#[test]
fn test_newton_filter() {
    assert!(satisfies_newton_inequalities_bigint_coeffs(&b(&[
        1, 4, 6, 4, 1
    ])));
    assert!(!satisfies_newton_inequalities_bigint_coeffs(&b(&[
        1, 1, 10, 1
    ])));
}

#[test]
fn test_kurtz_condition_filter() {
    // (t + 1)(t + 10)(t + 100) = 1000 + 1110t + 111t^2 + t^3.
    let widely_spaced = b(&[1000, 1110, 111, 1]);
    assert!(satisfies_kurtz_condition_bigint_coeffs(&widely_spaced));
    assert_eq!(
        is_real_rooted_one_signed_bigint_coeffs(&widely_spaced),
        Some(true)
    );

    // (1+t)^4 is real-rooted, but the Kurtz condition is intentionally
    // much stronger than real-rootedness.
    assert!(!satisfies_kurtz_condition_bigint_coeffs(&b(&[
        1, 4, 6, 4, 1
    ])));
    assert!(!satisfies_kurtz_condition_bigint_coeffs(&b(&[1, 0, 1])));
    assert!(!satisfies_kurtz_condition_bigint_coeffs(&b(&[
        -1000, -1110, -111, -1
    ])));
    assert!(!satisfies_kurtz_condition_bigint_coeffs(&b(&[1, 10])));
}

#[test]
fn test_prs_real_rootedness_general() {
    assert!(is_real_rooted_prs_bigint_coeffs(&b(&[-6, 11, -6, 1])));
    assert!(!is_real_rooted_prs_bigint_coeffs(&b(&[1, 0, 1])));
    assert!(is_real_rooted_prs_bigint_coeffs(&b(&[0, 0, 1, 2, 1]))); // x^2(1+x)^2
}

#[test]
fn test_low_degree_closed_form_edge_cases() {
    let cases = [
        (vec![], true),
        (vec![3], true),
        (vec![1, 0, 1], false),
        (vec![1, -2, 1], true),
        (vec![1, 0, 2, 0, 1], false), // (x^2 + 1)^2
        (vec![1, 0, -2, 0, 1], true), // (x^2 - 1)^2
        (mul_i64(&[1, -2, 1], &[1, 0, 1]), false),
        (mul_i64(&mul_i64(&[1, -2, 1], &[-2, 1]), &[-3, 1]), true),
        (vec![4, 0, -5, 0, 1], true), // roots -2, -1, 1, 2
        (vec![1, 0, 0, 0, 1], false), // no real roots
    ];

    for (coeffs, expected) in cases {
        assert_eq!(low_degree_is_real_rooted(&b(&coeffs)), Some(expected));
    }
}

#[test]
fn test_low_degree_closed_forms_agree_with_prs_exhaustively() {
    for nominal_degree in 0usize..=4 {
        let total = 7usize.pow((nominal_degree + 1) as u32);
        for mut mask in 0..total {
            let mut coeffs = Vec::with_capacity(nominal_degree + 1);
            for _ in 0..=nominal_degree {
                coeffs.push(BigInt::from((mask % 7) as i64 - 3));
                mask /= 7;
            }

            assert_eq!(
                low_degree_is_real_rooted(&coeffs),
                Some(is_real_rooted_prs_bigint_coeffs(&coeffs)),
                "coeffs={coeffs:?}"
            );
        }
    }
}

#[test]
fn test_low_degree_closed_forms_agree_with_prs_on_large_coefficients() {
    let mut state = 0x6a09_e667_f3bc_c909u64;
    for _ in 0..4096 {
        let mut coeffs = Vec::with_capacity(5);
        for _ in 0..5 {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            coeffs.push(BigInt::from(((state >> 32) % 2_000_001) as i64 - 1_000_000));
        }
        if coeffs[4].is_zero() {
            coeffs[4] = BigInt::one();
        }

        assert_eq!(
            low_degree_is_real_rooted(&coeffs),
            Some(is_real_rooted_prs_bigint_coeffs(&coeffs)),
            "coeffs={coeffs:?}"
        );
    }
}

#[test]
fn test_derivative_prs_reuses_gcd_chain_on_small_grid() {
    for degree in 0usize..=4 {
        let total = 5usize.pow((degree + 1) as u32);
        for mut mask in 0..total {
            let mut coeffs = Vec::with_capacity(degree + 1);
            for _ in 0..=degree {
                coeffs.push(BigInt::from((mask % 5) as i64 - 2));
                mask /= 5;
            }

            let (direct, direct_squarefree_degree) = derivative_sturm_chain(&coeffs);
            let summary = derivative_sturm_summary(&coeffs, true, true);
            let squarefree = squarefree_part_primitive_bigint_coeffs(&coeffs);
            let reference_squarefree_degree = crate::arithmetic::degree(&squarefree).unwrap_or(0);
            let reference = sturm_chain_from_squarefree(squarefree);

            assert_eq!(direct_squarefree_degree, reference_squarefree_degree);
            assert_eq!(summary.squarefree_degree, reference_squarefree_degree);
            assert_eq!(
                variations_at_neg_infinity(&direct)
                    .saturating_sub(variations_at_pos_infinity(&direct)),
                variations_at_neg_infinity(&reference)
                    .saturating_sub(variations_at_pos_infinity(&reference)),
                "coeffs={coeffs:?}"
            );
            assert_eq!(
                summary.real_roots,
                variations_at_neg_infinity(&reference)
                    .saturating_sub(variations_at_pos_infinity(&reference)),
                "coeffs={coeffs:?}"
            );
            assert_eq!(
                variations_at_zero_plus(&direct)
                    .saturating_sub(variations_at_pos_infinity(&direct)),
                variations_at_zero_plus(&reference)
                    .saturating_sub(variations_at_pos_infinity(&reference)),
                "coeffs={coeffs:?}"
            );
            assert_eq!(
                summary.positive_roots,
                variations_at_zero_plus(&reference)
                    .saturating_sub(variations_at_pos_infinity(&reference)),
                "coeffs={coeffs:?}"
            );
        }
    }
}

#[test]
fn test_adaptive_small_degree_prs_agrees_with_primitive_prs() {
    let mut state = 0xbb67_ae85_84ca_a73bu64;
    for d in 5usize..=7 {
        let mut cases = Vec::new();
        for _ in 0..512 {
            let mut coefficients = Vec::with_capacity(d + 1);
            for _ in 0..=d {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                coefficients.push(BigInt::from(((state >> 32) % 101) as i64 - 50));
            }
            if coefficients[d].is_zero() {
                coefficients[d] = BigInt::one();
            }
            cases.push(coefficients);
        }

        // Include all-real inputs with repeated roots, since the final
        // remainder degree also supplies the square-free degree.
        for repeated_root in 1..=8usize {
            let mut polynomial = vec![BigInt::one()];
            for root_index in 0..d {
                let root = if root_index < 2 {
                    repeated_root
                } else {
                    repeated_root + root_index
                };
                polynomial = poly_mul(&polynomial, &[BigInt::from(root), BigInt::one()]);
            }
            cases.push(polynomial);
        }

        for coefficients in cases {
            let p = primitive_positive(coefficients);
            let derivative = poly_derivative(&p);
            let primitive_chain = signed_remainder_sequence(&p, &derivative);
            let primitive_gcd_degree = primitive_chain
                .last()
                .and_then(|last| degree(last))
                .unwrap_or(0);
            let primitive_squarefree_degree = d - primitive_gcd_degree;
            let (adaptive_chain, adaptive_squarefree_degree) = derivative_sturm_chain(&p);
            let summary = derivative_sturm_summary(&p, true, true);

            assert_eq!(adaptive_squarefree_degree, primitive_squarefree_degree);
            assert_eq!(summary.squarefree_degree, primitive_squarefree_degree);
            assert_eq!(
                variations_at_neg_infinity(&adaptive_chain)
                    .saturating_sub(variations_at_pos_infinity(&adaptive_chain)),
                variations_at_neg_infinity(&primitive_chain)
                    .saturating_sub(variations_at_pos_infinity(&primitive_chain)),
                "degree={d}, coeffs={p:?}"
            );
            assert_eq!(
                summary.real_roots,
                variations_at_neg_infinity(&primitive_chain)
                    .saturating_sub(variations_at_pos_infinity(&primitive_chain)),
                "degree={d}, coeffs={p:?}"
            );
            assert_eq!(
                summary.positive_roots,
                variations_at_zero_plus(&primitive_chain)
                    .saturating_sub(variations_at_pos_infinity(&primitive_chain)),
                "degree={d}, coeffs={p:?}"
            );
        }
    }
}

#[test]
fn test_subresultant_prs_agrees_with_primitive_prs() {
    let mut state = 0x3c6e_f372_fe94_f82bu64;
    for d in 8usize..=14 {
        let mut cases = Vec::new();
        for _ in 0..256 {
            let mut coefficients = Vec::with_capacity(d + 1);
            for _ in 0..=d {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                coefficients.push(BigInt::from(((state >> 32) % 101) as i64 - 50));
            }
            if coefficients[d].is_zero() {
                coefficients[d] = BigInt::one();
            }
            cases.push(coefficients);
        }

        let mut repeated = vec![BigInt::one()];
        for root_index in 0..d {
            let root = if root_index < 3 { 2 } else { root_index + 1 };
            repeated = poly_mul(&repeated, &[BigInt::from(root), BigInt::one()]);
        }
        cases.push(repeated);

        for coefficients in cases {
            let p = primitive_positive(coefficients);
            let derivative = poly_derivative(&p);
            let primitive_chain = signed_remainder_sequence(&p, &derivative);
            let subresultant_chain = signed_remainder_sequence_subresultant(&p, &derivative);
            let summary = derivative_sturm_summary(&p, true, true);
            let primitive_gcd_degree = primitive_chain
                .last()
                .and_then(|last| degree(last))
                .unwrap_or(0);
            let primitive_squarefree_degree = d - primitive_gcd_degree;

            assert_eq!(
                subresultant_chain.last().and_then(|last| degree(last)),
                primitive_chain.last().and_then(|last| degree(last)),
                "degree={d}, coeffs={p:?}"
            );
            assert_eq!(
                variations_at_neg_infinity(&subresultant_chain)
                    .saturating_sub(variations_at_pos_infinity(&subresultant_chain)),
                variations_at_neg_infinity(&primitive_chain)
                    .saturating_sub(variations_at_pos_infinity(&primitive_chain)),
                "degree={d}, coeffs={p:?}"
            );
            assert_eq!(
                count_positive_roots_from_sturm_chain(&subresultant_chain),
                count_positive_roots_from_sturm_chain(&primitive_chain),
                "degree={d}, coeffs={p:?}"
            );
            assert_eq!(
                summary.squarefree_degree, primitive_squarefree_degree,
                "degree={d}, coeffs={p:?}"
            );
            assert_eq!(
                summary.real_roots,
                variations_at_neg_infinity(&primitive_chain)
                    .saturating_sub(variations_at_pos_infinity(&primitive_chain)),
                "degree={d}, coeffs={p:?}"
            );
            assert_eq!(
                summary.positive_roots,
                count_positive_roots_from_sturm_chain(&primitive_chain),
                "degree={d}, coeffs={p:?}"
            );
        }
    }
}

#[test]
fn test_uspensky_root_counts() {
    assert_eq!(
        count_positive_roots_uspensky_bigint_coeffs(&b(&[6, -7, 0, 1])),
        2
    ); // (x-1)(x-2)(x+3)
    assert_eq!(
        count_real_roots_uspensky_bigint_coeffs(&b(&[-6, 11, -6, 1])),
        3
    );
    assert_eq!(count_real_roots_uspensky_bigint_coeffs(&b(&[1, 0, 1])), 0);
    assert_eq!(
        count_positive_roots_uspensky_bigint_coeffs(&b(&[-2, 7, -7, 2])),
        3
    ); // (x-1)(x-2)(2x-1), including dyadic boundary roots
    assert_eq!(
        count_real_roots_uspensky_bigint_coeffs(&b(&[0, 0, 1, -2, 1])),
        2
    ); // x^2(x-1)^2 has distinct roots 0 and 1
}

#[test]
fn test_uspensky_handles_large_close_roots() {
    let center = BigInt::from(10u32).pow(20);
    let p = vec![
        &center * (&center + BigInt::one()),
        -(BigInt::from(2) * &center + BigInt::one()),
        BigInt::one(),
    ];
    assert_eq!(count_positive_roots_uspensky_bigint_coeffs(&p), 2);
    assert!(is_real_rooted_uspensky_bigint_coeffs(&p));
}

#[test]
fn test_uspensky_agrees_with_prs_on_small_grid() {
    for degree in 0usize..=4 {
        let total = 5usize.pow((degree + 1) as u32);
        for mut mask in 0..total {
            let mut coeffs = Vec::with_capacity(degree + 1);
            for _ in 0..=degree {
                coeffs.push((mask % 5) as i64 - 2);
                mask /= 5;
            }

            let bigint = b(&coeffs);
            assert_eq!(
                count_real_roots_uspensky_bigint_coeffs(&bigint),
                count_real_roots_prs_bigint_coeffs(&bigint),
                "coeffs={coeffs:?}"
            );
        }
    }
}

#[test]
fn test_uspensky_i64_wrapper_agrees_with_bigint_and_prs_on_grid() {
    // Every polynomial of degree <= 4 with coefficients in -2..=2: the `i64` wrapper
    // must equal the `BigInt` Uspensky test and the independent PRS test.
    for degree in 0usize..=4 {
        let total = 5usize.pow((degree + 1) as u32);
        for mut mask in 0..total {
            let mut coeffs = Vec::with_capacity(degree + 1);
            for _ in 0..=degree {
                coeffs.push((mask % 5) as i64 - 2);
                mask /= 5;
            }
            let bigint = b(&coeffs);
            let expected = is_real_rooted_prs_bigint_coeffs(&bigint);
            assert_eq!(is_real_rooted_uspensky_i64(&coeffs), expected, "{coeffs:?}");
            assert_eq!(is_real_rooted_uspensky_bigint_coeffs(&bigint), expected);
        }
    }
    // Known cases: Eulerian A_4 (A008292) is real-rooted, t^2 + 1 and the
    // classical non-real-rooted 1 + 43t + 196t^2 + 168t^3 + 23t^4 + t^5 are not.
    assert!(is_real_rooted_uspensky_i64(&[1, 11, 11, 1]));
    assert!(!is_real_rooted_uspensky_i64(&[1, 0, 1]));
    assert!(!is_real_rooted_uspensky_i64(&[1, 43, 196, 168, 23, 1]));
    assert!(is_real_rooted_uspensky_i64(&[]));
    assert!(is_real_rooted_uspensky_i64(&[7]));
    // Repeated roots count once: (t - 1)^2 (t + 2) and t^3 are real-rooted.
    assert!(is_real_rooted_uspensky_i64(&[2, -3, 0, 1]));
    assert!(is_real_rooted_uspensky_i64(&[0, 0, 0, 1]));
}

#[test]
fn test_primitive_sturm_max_coefficient_bits_known_chains() {
    // Hand/independent computation of the primitive Sturm chains:
    //   t^2 - 2:            [-2, 0, 1], [0, 1], [1]                         -> 2 bits
    //   (t-1)(t-2)(t-3):    [-6, 11, -6, 1], [11, -12, 3], [-2, 1], [1]     -> 4 bits
    //   Eulerian A_4:       [1, 11, 11, 1], [11, 22, 3], [7, 11], [1]       -> 5 bits
    assert_eq!(primitive_sturm_max_coefficient_bits(&b(&[-2, 0, 1])), 2);
    assert_eq!(
        primitive_sturm_max_coefficient_bits(&b(&[-6, 11, -6, 1])),
        4
    );
    assert_eq!(primitive_sturm_max_coefficient_bits(&b(&[1, 11, 11, 1])), 5);
    // Zero polynomial and constants.
    assert_eq!(primitive_sturm_max_coefficient_bits(&[]), 0);
    assert_eq!(primitive_sturm_max_coefficient_bits(&b(&[0, 0])), 0);
    assert_eq!(primitive_sturm_max_coefficient_bits(&b(&[1000])), 1);
    // The chain is built from the primitive squarefree part: scaling by content
    // and repeating factors do not change it ((t-1)^3 has squarefree part t - 1).
    let scale = BigInt::from(10).pow(30);
    let scaled: Vec<BigInt> = b(&[-6, 11, -6, 1]).iter().map(|c| c * &scale).collect();
    assert_eq!(primitive_sturm_max_coefficient_bits(&scaled), 4);
    assert_eq!(primitive_sturm_max_coefficient_bits(&b(&[-1, 3, -3, 1])), 1);
    // The chain always contains the squarefree polynomial itself, so the value is
    // at least the bit size of its largest coefficient.
    let mut product = vec![BigInt::one()];
    for root in 1..=12 {
        product = poly_mul(&product, &[BigInt::from(-root), BigInt::one()]);
    }
    let largest = product.iter().map(|c| c.bits()).max().unwrap();
    assert!(primitive_sturm_max_coefficient_bits(&product) >= largest);
}

#[test]
fn test_uspensky_agrees_with_prs_on_deterministic_higher_degree_cases() {
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    for degree in 5usize..=8 {
        for _ in 0..64 {
            let mut coeffs = Vec::with_capacity(degree + 1);
            for _ in 0..=degree {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                coeffs.push(((state >> 32) % 21) as i64 - 10);
            }
            if coeffs[degree] == 0 {
                coeffs[degree] = 1;
            }

            let bigint = b(&coeffs);
            assert_eq!(
                count_real_roots_uspensky_bigint_coeffs(&bigint),
                count_real_roots_prs_bigint_coeffs(&bigint),
                "coeffs={coeffs:?}"
            );
        }
    }
}

#[test]
fn test_agrees_with_existing_examples() {
    let cases = [
        (vec![1, 2, 1], true),
        (vec![1, 11, 11, 1], true),
        (vec![1, 0, 1], false),
        (vec![1, 43, 196, 168, 23, 1], false),
        (vec![-15, 23, -9, 1], true),
    ];
    for (coeffs, expected) in cases {
        assert_eq!(
            is_real_rooted_fast_bigint_coeffs(&b(&coeffs)),
            expected,
            "coeffs={coeffs:?}"
        );
    }
}

#[test]
fn test_products_with_known_real_roots() {
    let mut p = vec![1i64];
    for a in 1..=8 {
        p = mul_i64(&p, &[a, 1]); // (x+a)
        assert!(is_real_rooted_fast_bigint_coeffs(&b(&p)));
    }
    let complex_factor = mul_i64(&p, &[1, 0, 1]);
    assert!(!is_real_rooted_fast_bigint_coeffs(&b(&complex_factor)));
}

#[test]
fn test_one_signed_path_agrees_with_both_backends_on_benchmarked_families() {
    let eulerian = generated_eulerian(36);
    let eulerian_below_cutoff = generated_eulerian(35);
    let narayana = generated_narayana(36);
    let type_b_eulerian = generated_type_b_eulerian(35);

    let mut touchard = vec![BigInt::one()];
    for _ in 0..36 {
        let mut next = vec![BigInt::zero(); touchard.len() + 1];
        for (k, coefficient) in touchard.iter().enumerate() {
            next[k + 1] += coefficient;
            if k > 0 {
                next[k] += coefficient * BigInt::from(k);
            }
        }
        touchard = trim(next);
    }
    let touchard = trim_slice(strip_initial_zeros_bigint(&touchard));

    let mut product = vec![BigInt::one()];
    for root in 1..=35 {
        product = poly_mul(&product, &[BigInt::from(root), BigInt::one()]);
    }

    let mut dented = eulerian.clone();
    let middle = dented.len() / 2;
    dented[middle] = &dented[middle] * 3 / 4;

    for (polynomial, expected) in [
        (&eulerian, true),
        (&eulerian_below_cutoff, true),
        (&type_b_eulerian, true),
        (&touchard, true),
        (&narayana, true),
        (&product, true),
        (&dented, false),
    ] {
        assert_eq!(is_real_rooted_fast_bigint_coeffs(polynomial), expected);
        assert_eq!(is_real_rooted_prs_bigint_coeffs(polynomial), expected);
        assert_eq!(is_real_rooted_uspensky_bigint_coeffs(polynomial), expected);
    }
}

fn dyadic_to_rational(value: &DyadicRational) -> BigRational {
    BigRational::new(value.numerator().clone(), value.denominator())
}

fn rational_value(coeffs: &[BigInt], point: &BigRational) -> BigRational {
    coeffs.iter().rev().fold(BigRational::zero(), |acc, c| {
        acc * point + BigRational::from_integer(c.clone())
    })
}

/// Check the documented isolation contract against independent evaluation.
fn assert_valid_isolation(coeffs: &[BigInt]) -> Vec<IsolatedRealRoot> {
    let roots = isolate_real_roots_bigint_coeffs(coeffs);
    assert_eq!(roots.len(), count_real_roots_prs_bigint_coeffs(coeffs));
    let squarefree = squarefree_part_primitive_bigint_coeffs(coeffs);
    for root in &roots {
        match root {
            IsolatedRealRoot::Exact(value) => {
                assert!(rational_value(coeffs, &dyadic_to_rational(value)).is_zero());
            }
            IsolatedRealRoot::Interval { lower, upper } => {
                assert!(lower < upper);
                let low = rational_value(&squarefree, &dyadic_to_rational(lower));
                let high = rational_value(&squarefree, &dyadic_to_rational(upper));
                assert!(!low.is_zero() && !high.is_zero());
                // A simple root of the square-free part changes its sign.
                assert!(low.is_negative() != high.is_negative());
            }
        }
    }
    for pair in roots.windows(2) {
        let strict = matches!(pair[0], IsolatedRealRoot::Exact(_))
            || matches!(pair[1], IsolatedRealRoot::Exact(_));
        if strict {
            assert!(pair[0].upper() < pair[1].lower());
        } else {
            assert!(pair[0].upper() <= pair[1].lower());
        }
    }
    roots
}

#[test]
fn dyadic_rational_is_reduced_and_ordered() {
    let half = DyadicRational::new(BigInt::from(4), 3);
    assert_eq!(half.numerator(), &BigInt::from(1));
    assert_eq!(half.exponent(), 1);
    assert_eq!(DyadicRational::new(BigInt::zero(), 7).exponent(), 0);
    let one = DyadicRational::from_integer(BigInt::one());
    assert!(half < one);
    assert_eq!(half.midpoint(&one), DyadicRational::new(BigInt::from(3), 2));
    assert_eq!(half.to_string(), "1/2^1");
}

#[test]
fn isolates_constant_and_irrational_roots() {
    assert!(isolate_real_roots_bigint_coeffs(&[]).is_empty());
    assert!(isolate_real_roots_bigint_coeffs(&b(&[5])).is_empty());
    assert!(isolate_real_roots_bigint_coeffs(&b(&[1, 0, 1])).is_empty());

    // x^2 - 2 has the two irrational roots +-sqrt(2).
    let roots = assert_valid_isolation(&b(&[-2, 0, 1]));
    assert_eq!(roots.len(), 2);
    assert!(roots
        .iter()
        .all(|root| matches!(root, IsolatedRealRoot::Interval { .. })));
    assert!(roots[0].upper() <= &DyadicRational::from_integer(BigInt::zero()));
}

#[test]
fn isolates_exact_dyadic_roots_hit_by_bisection() {
    // x (x - 1) (2x - 1) (x + 3): the bisection midpoints hit 0 and 1/2.
    let mut p = b(&[1]);
    for factor in [[0, 1], [-1, 1], [-1, 2], [3, 1]] {
        p = poly_mul(&p, &b(&factor));
    }
    let roots = assert_valid_isolation(&p);
    assert_eq!(roots.len(), 4);
    assert!(roots
        .iter()
        .any(|root| *root == IsolatedRealRoot::Exact(DyadicRational::from_integer(0.into()))));
}

#[test]
fn isolates_repeated_and_close_roots() {
    // (3x - 1)^2 (3x - 2) (x^2 + 1) (1000x - 333)
    let mut p = b(&[1]);
    for factor in [
        vec![-1, 3],
        vec![-1, 3],
        vec![-2, 3],
        vec![1, 0, 1],
        vec![-333, 1000],
    ] {
        p = poly_mul(&p, &b(&factor));
    }
    let roots = assert_valid_isolation(&p);
    assert_eq!(roots.len(), 3);

    // A product of many linear factors, including negative roots.
    let mut wilkinson = b(&[1]);
    for root in -6..=7 {
        wilkinson = poly_mul(&wilkinson, &b(&[-root, 1]));
    }
    assert_eq!(assert_valid_isolation(&wilkinson).len(), 14);
}
