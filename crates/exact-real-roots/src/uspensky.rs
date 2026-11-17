//! Independent Uspensky/Descartes root counting over the integers.
//! Coefficients are ascending; root isolation uses the same strict root bound.

use crate::arithmetic::{
    degree, exact_div, primitive_keep_sign, sign_i8, squarefree_part_primitive_bigint_coeffs,
    strip_initial_zeros_bigint, trim, trim_slice,
};
use crate::coefficient_tests::alternating_neg_argument;
use crate::prs::sign_variations;
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

fn ceil_div_positive(numerator: &BigInt, denominator: &BigInt) -> BigInt {
    debug_assert!(!numerator.is_negative());
    debug_assert!(denominator.is_positive());
    if numerator.is_zero() {
        BigInt::zero()
    } else {
        (numerator + denominator - BigInt::one()) / denominator
    }
}

/// A strict integer version of Fujiwara's root bound.
///
/// Every complex root of `p` has modulus strictly smaller than the returned
/// integer.  Using the exponent `1 / (degree - i)` is important here: the
/// simpler Cauchy bound can be exponentially too large for products of linear
/// factors, which in turn creates needless Uspensky subdivision levels.
pub(super) fn strict_fujiwara_root_bound(coeffs: &[BigInt]) -> BigInt {
    let d = degree(coeffs).expect("root bound requires a nonzero polynomial");
    if d == 0 {
        return BigInt::one();
    }
    let lc = coeffs[d].abs();
    let mut maximum_root = BigInt::zero();

    for (i, coefficient) in coeffs.iter().enumerate().take(d) {
        if coefficient.is_zero() {
            continue;
        }
        let quotient_ceiling = ceil_div_positive(&coefficient.abs(), &lc);
        let exponent = u32::try_from(d - i).expect("polynomial degree should fit in u32");
        let mut root_ceiling = quotient_ceiling.nth_root(exponent);
        if root_ceiling.pow(exponent) < quotient_ceiling {
            root_ceiling += BigInt::one();
        }
        maximum_root = maximum_root.max(root_ceiling);
    }

    BigInt::from(2) * maximum_root.max(BigInt::one()) + BigInt::one()
}

fn scale_argument(coeffs: &[BigInt], scale: &BigInt) -> Vec<BigInt> {
    let mut power = BigInt::one();
    let mut result = Vec::with_capacity(coeffs.len());
    for coefficient in coeffs {
        result.push(coefficient * &power);
        power *= scale;
    }
    trim(result)
}

fn shift_argument(coeffs: &[BigInt], shift: &BigInt) -> Vec<BigInt> {
    let Some(d) = degree(coeffs) else {
        return Vec::new();
    };
    let mut result: Vec<BigInt> = Vec::new();

    for coefficient in coeffs.iter().take(d + 1).rev() {
        let mut next = vec![BigInt::zero(); result.len() + 1];
        for (i, value) in result.iter().enumerate() {
            next[i] += value * shift;
            next[i + 1] += value;
        }
        next[0] += coefficient;
        result = next;
    }
    trim(result)
}

/// Return `(1+x)^d p(x/(1+x))`, where `d = degree(p)`.
///
/// Positive roots of the result correspond to roots of `p` in `(0, 1)`.
fn transform_unit_interval(coeffs: &[BigInt]) -> Vec<BigInt> {
    let Some(d) = degree(coeffs) else {
        return Vec::new();
    };
    let mut result = vec![BigInt::zero(); d + 1];

    for (k, coefficient) in coeffs.iter().enumerate().take(d + 1) {
        if coefficient.is_zero() {
            continue;
        }
        let remaining = d - k;
        let mut binomial = BigInt::one();
        for offset in 0..=remaining {
            result[k + offset] += coefficient * &binomial;
            if offset < remaining {
                binomial = binomial * BigInt::from(remaining - offset) / BigInt::from(offset + 1);
            }
        }
    }
    trim(result)
}

/// Return `p(x+1)`.
///
/// Positive roots of the result correspond to roots of `p` in `(1, +infinity)`.
fn shift_argument_by_one(coeffs: &[BigInt]) -> Vec<BigInt> {
    shift_argument(coeffs, &BigInt::one())
}

fn value_at_one(coeffs: &[BigInt]) -> BigInt {
    coeffs.iter().sum()
}

fn value_at_integer(coeffs: &[BigInt], value: &BigInt) -> BigInt {
    coeffs
        .iter()
        .rev()
        .fold(BigInt::zero(), |acc, coefficient| acc * value + coefficient)
}

fn descartes_sign_variations(coeffs: &[BigInt]) -> usize {
    sign_variations(coeffs.iter().map(sign_i8))
}

/// Count roots represented by one open Uspensky interval transform.
fn count_roots_in_uspensky_node(initial: Vec<BigInt>) -> usize {
    let linear_midpoint_factor = vec![-BigInt::one(), BigInt::one()];
    let mut stack = vec![initial];
    let mut roots = 0usize;

    while let Some(mut node) = stack.pop() {
        let variations = descartes_sign_variations(&node);
        if variations <= 1 {
            roots += variations;
            continue;
        }

        if value_at_one(&node).is_zero() {
            roots += 1;
            node = primitive_keep_sign(exact_div(&node, &linear_midpoint_factor));
            let remaining_variations = descartes_sign_variations(&node);
            if remaining_variations <= 1 {
                roots += remaining_variations;
                continue;
            }
        }

        let left = primitive_keep_sign(transform_unit_interval(&node));
        let right = primitive_keep_sign(shift_argument_by_one(&node));
        stack.push(right);
        stack.push(left);
    }
    roots
}

/// Transform the open interval `(left, right)` to `(0, +infinity)`.
fn transform_open_integer_interval(
    coeffs: &[BigInt],
    left: &BigInt,
    right: &BigInt,
) -> Vec<BigInt> {
    debug_assert!(left < right);
    let shifted = shift_argument(coeffs, left);
    let width = right - left;
    let scaled = scale_argument(&shifted, &width);
    primitive_keep_sign(transform_unit_interval(&scaled))
}

/// Count roots in `(1, +infinity)` using dyadic magnitude bands.
fn count_roots_above_one_uspensky_squarefree(coeffs: &[BigInt]) -> usize {
    let mut p = trim_slice(coeffs);
    if degree(&p).unwrap_or(0) == 0 {
        return 0;
    }
    debug_assert!(!value_at_one(&p).is_zero());

    let bound = strict_fujiwara_root_bound(&p);
    let mut left = BigInt::one();
    let mut roots = 0usize;

    while left < bound {
        let doubled = BigInt::from(2) * &left;
        let right = doubled.min(bound.clone());

        if right < bound && value_at_integer(&p, &right).is_zero() {
            roots += 1;
            p = primitive_keep_sign(exact_div(&p, &[-right.clone(), BigInt::one()]));
        }

        let transformed = transform_open_integer_interval(&p, &left, &right);
        roots += count_roots_in_uspensky_node(transformed);
        left = right;
    }
    roots
}

/// Count positive roots of a square-free polynomial with nonzero constant term.
///
/// Roots greater than one are grouped into dyadic magnitude bands.  Roots in
/// `(0, 1)` are treated as reciprocals of roots greater than one of the reversed
/// polynomial.  This avoids placing roots of very different magnitudes into a
/// single bounded interval, which is particularly important for reciprocal
/// families such as Eulerian polynomials.
fn count_positive_roots_uspensky_squarefree(coeffs: &[BigInt]) -> usize {
    let mut p = trim_slice(coeffs);
    let Some(d) = degree(&p) else {
        return 0;
    };
    if d == 0 {
        return 0;
    }
    debug_assert!(!p[0].is_zero());

    let initial_variations = descartes_sign_variations(&p);
    if initial_variations <= 1 {
        return initial_variations;
    }

    let mut roots = 0usize;
    if value_at_one(&p).is_zero() {
        roots += 1;
        p = primitive_keep_sign(exact_div(&p, &[-BigInt::one(), BigInt::one()]));
    }

    roots += count_roots_above_one_uspensky_squarefree(&p);
    let reciprocal: Vec<BigInt> = p.into_iter().rev().collect();
    roots + count_roots_above_one_uspensky_squarefree(&reciprocal)
}

/// Count distinct positive roots by exact Uspensky/Descartes subdivision.
///
/// The polynomial is replaced by its square-free part first.  The method uses
/// only integer arithmetic: no finite fields, floating point, or approximate
/// roots are involved.
pub fn count_positive_roots_uspensky_bigint_coeffs(coeffs: &[BigInt]) -> usize {
    let squarefree = squarefree_part_primitive_bigint_coeffs(coeffs);
    let without_zero = strip_initial_zeros_bigint(&squarefree);
    count_positive_roots_uspensky_squarefree(without_zero)
}

/// Count distinct real roots by exact Uspensky/Descartes subdivision.
pub fn count_real_roots_uspensky_bigint_coeffs(coeffs: &[BigInt]) -> usize {
    let squarefree = squarefree_part_primitive_bigint_coeffs(coeffs);
    let zero_root = usize::from(squarefree.first().is_some_and(BigInt::is_zero));
    let without_zero = strip_initial_zeros_bigint(&squarefree);
    if degree(without_zero).unwrap_or(0) == 0 {
        return zero_root;
    }

    let positive = count_positive_roots_uspensky_squarefree(without_zero);
    let reflected = alternating_neg_argument(without_zero);
    let negative = count_positive_roots_uspensky_squarefree(&reflected);
    zero_root + positive + negative
}

/// Exact real-rootedness via Uspensky's method and Descartes' rule of signs.
pub fn is_real_rooted_uspensky_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let squarefree = squarefree_part_primitive_bigint_coeffs(coeffs);
    let squarefree_degree = degree(&squarefree).unwrap_or(0);
    let zero_root = usize::from(squarefree.first().is_some_and(BigInt::is_zero));
    let without_zero = strip_initial_zeros_bigint(&squarefree);
    let positive = count_positive_roots_uspensky_squarefree(without_zero);
    let reflected = alternating_neg_argument(without_zero);
    let negative = count_positive_roots_uspensky_squarefree(&reflected);
    zero_root + positive + negative == squarefree_degree
}

/// Convenience Uspensky wrapper for `i64` coefficient vectors.
pub fn is_real_rooted_uspensky_i64(coeffs: &[i64]) -> bool {
    let coeffs: Vec<BigInt> = coeffs.iter().map(|&c| BigInt::from(c)).collect();
    is_real_rooted_uspensky_bigint_coeffs(&coeffs)
}
