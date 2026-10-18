//! Exact real-root counting over the integers.
//!
//! This module provides the default exact backend for boolean
//! real-rootedness. It uses a primitive pseudo-remainder Sturm sequence for
//! general inputs and adaptively selects an exact Uspensky/Descartes path for
//! a conservative class of large one-signed inputs. Primitive PRS keeps
//! intermediate coefficient growth much smaller than a naive rational PRS,
//! while both paths avoid the large exact PSD matrices used by the Bézout
//! criterion.
//!
//! The positive-coefficient path is intended for combinatorial polynomials.  If
//! all coefficients have one sign, roots can only be non-positive; after
//! removing powers of `t`, real-rootedness of `f(t)` is equivalent to all roots
//! of `f(-t)` being positive.
//!
//! Public functions in this module count distinct roots of the square-free
//! part.  Real-rootedness tests compare that count with the square-free degree,
//! so repeated real roots are handled without needing multiplicities from the
//! Sturm sequence.
//!
//! An independent Uspensky/Descartes implementation is also available for
//! exact comparison.  It bounds the roots, applies homographic subdivision,
//! and uses Descartes' rule of signs to discard or certify open intervals.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

/// Exact sign distribution of one polynomial evaluated at the real roots of another.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RootSignCounts {
    /// Number of distinct real roots `r` of the root polynomial with `q(r)>0`.
    pub positive: usize,
    /// Number of distinct real roots `r` of the root polynomial with `q(r)<0`.
    pub negative: usize,
    /// Number of distinct real roots `r` of the root polynomial with `q(r)=0`.
    pub zero: usize,
    /// Total number of distinct real roots of the square-free part of the root polynomial.
    pub total: usize,
}

fn trim(mut p: Vec<BigInt>) -> Vec<BigInt> {
    while p.last().is_some_and(|c| c.is_zero()) {
        p.pop();
    }
    p
}

fn trim_slice(p: &[BigInt]) -> Vec<BigInt> {
    trim(p.to_vec())
}

fn degree(p: &[BigInt]) -> Option<usize> {
    p.iter().rposition(|c| !c.is_zero())
}

fn leading(p: &[BigInt]) -> BigInt {
    degree(p).map(|d| p[d].clone()).unwrap_or_else(BigInt::zero)
}

fn sign_i8(x: &BigInt) -> i8 {
    if x.is_positive() {
        1
    } else if x.is_negative() {
        -1
    } else {
        0
    }
}

fn content(p: &[BigInt]) -> BigInt {
    let mut g = BigInt::zero();
    for c in p {
        if !c.is_zero() {
            let a = c.abs();
            g = if g.is_zero() { a } else { g.gcd(&a) };
        }
    }
    g
}

fn primitive_keep_sign(p: Vec<BigInt>) -> Vec<BigInt> {
    let p = trim(p);
    let g = content(&p);
    if g.is_zero() || g.is_one() {
        p
    } else {
        trim(p.into_iter().map(|c| c / &g).collect())
    }
}

fn primitive_positive(p: Vec<BigInt>) -> Vec<BigInt> {
    let mut p = primitive_keep_sign(p);
    if leading(&p).is_negative() {
        for c in &mut p {
            *c = -c.clone();
        }
    }
    p
}

fn poly_neg(p: &[BigInt]) -> Vec<BigInt> {
    p.iter().map(|c| -c).collect()
}

fn poly_derivative(p: &[BigInt]) -> Vec<BigInt> {
    if p.len() <= 1 {
        return vec![];
    }
    trim(
        p.iter()
            .enumerate()
            .skip(1)
            .map(|(i, c)| c * BigInt::from(i))
            .collect(),
    )
}

fn poly_scale(p: &[BigInt], a: &BigInt) -> Vec<BigInt> {
    if a.is_zero() || p.is_empty() {
        return vec![];
    }
    trim(p.iter().map(|c| c * a).collect())
}

fn poly_mul(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r = vec![BigInt::zero(); a.len() + b.len() - 1];
    for (i, ca) in a.iter().enumerate() {
        if ca.is_zero() {
            continue;
        }
        for (j, cb) in b.iter().enumerate() {
            if !cb.is_zero() {
                r[i + j] += ca * cb;
            }
        }
    }
    trim(r)
}

fn poly_sub_shifted_scaled(
    a: &[BigInt],
    b: &[BigInt],
    shift: usize,
    scale: &BigInt,
) -> Vec<BigInt> {
    let n = a.len().max(b.len() + shift);
    let mut r = vec![BigInt::zero(); n];
    for (i, c) in a.iter().enumerate() {
        r[i] += c;
    }
    for (i, c) in b.iter().enumerate() {
        r[i + shift] -= c * scale;
    }
    trim(r)
}

/// Pseudo-remainder with a positive multiplier.
///
/// If `b_pos` is `b` multiplied by a sign so that its leading coefficient is
/// positive, this returns `lc(b_pos)^k * rem(a, b_pos)` for a nonnegative `k`.
/// The scalar is positive, so it does not change Sturm sign variations.
fn pseudo_remainder_positive_multiplier(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let da = match degree(a) {
        Some(d) => d,
        None => return vec![],
    };
    let db = degree(b).expect("pseudo_remainder: division by zero polynomial");
    if da < db {
        return trim_slice(a);
    }

    let mut b_pos = trim_slice(b);
    if leading(&b_pos).is_negative() {
        for c in &mut b_pos {
            *c = -c.clone();
        }
    }
    let lc = leading(&b_pos);
    debug_assert!(lc.is_positive());

    let mut r = trim_slice(a);
    let mut e = da - db + 1;
    while let Some(dr) = degree(&r) {
        if dr < db {
            break;
        }
        let c = leading(&r);
        let shift = dr - db;
        let scaled_r = poly_scale(&r, &lc);
        r = poly_sub_shifted_scaled(&scaled_r, &b_pos, shift, &c);
        e -= 1;
    }
    for _ in 0..e {
        r = poly_scale(&r, &lc);
    }
    trim(r)
}

fn primitive_gcd(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let mut r0 = primitive_positive(trim_slice(a));
    let mut r1 = primitive_positive(trim_slice(b));
    if r0.is_empty() {
        return r1;
    }
    if r1.is_empty() {
        return r0;
    }

    while !r1.is_empty() {
        let rem = primitive_positive(pseudo_remainder_positive_multiplier(&r0, &r1));
        r0 = r1;
        r1 = rem;
    }
    primitive_positive(r0)
}

fn exact_div(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let da = match degree(a) {
        Some(d) => d,
        None => return vec![],
    };
    let db = degree(b).expect("exact_div: division by zero polynomial");
    if da < db {
        return vec![];
    }

    let mut rem = trim_slice(a);
    let mut q = vec![BigInt::zero(); da - db + 1];
    let lc_b = leading(b);
    for k in (0..=da - db).rev() {
        let coeff = rem[db + k].clone();
        if coeff.is_zero() {
            continue;
        }
        assert!(
            (&coeff % &lc_b).is_zero(),
            "exact_div: non-exact leading quotient"
        );
        let qk = coeff / &lc_b;
        q[k] = qk.clone();
        for (j, bj) in b.iter().enumerate().take(db + 1) {
            rem[j + k] -= &qk * bj;
        }
    }
    debug_assert!(trim(rem).is_empty(), "exact_div: nonzero remainder");
    trim(q)
}

fn strip_initial_zeros_bigint(coeffs: &[BigInt]) -> &[BigInt] {
    match coeffs.iter().position(|c| !c.is_zero()) {
        Some(i) => &coeffs[i..],
        None => &[],
    }
}

/// Return the primitive, positive-leading square-free part.
///
/// This crate-private entry point is shared by the public square-free API and
/// the root-counting backends so they use the same fraction-free integer GCD.
pub(crate) fn squarefree_part_primitive_bigint_coeffs(coeffs: &[BigInt]) -> Vec<BigInt> {
    let p = primitive_positive(trim_slice(coeffs));
    if degree(&p).unwrap_or(0) == 0 {
        return p;
    }
    let dp = poly_derivative(&p);
    if dp.is_empty() {
        return p;
    }
    let g = primitive_gcd(&p, &dp);
    if degree(&g).unwrap_or(0) == 0 {
        p
    } else {
        primitive_positive(exact_div(&p, &g))
    }
}

fn sturm_chain_from_squarefree(p0: Vec<BigInt>) -> Vec<Vec<BigInt>> {
    if p0.is_empty() || degree(&p0).unwrap_or(0) == 0 {
        return vec![p0];
    }
    let p1 = primitive_keep_sign(poly_derivative(&p0));
    let mut chain = vec![p0, p1];

    loop {
        let n = chain.len();
        let prem = pseudo_remainder_positive_multiplier(&chain[n - 2], &chain[n - 1]);
        if prem.is_empty() {
            break;
        }
        let next = primitive_keep_sign(poly_neg(&primitive_keep_sign(prem)));
        if next.is_empty() {
            break;
        }
        chain.push(next);
    }
    chain
}

fn sturm_chain_squarefree(coeffs: &[BigInt]) -> Vec<Vec<BigInt>> {
    sturm_chain_from_squarefree(squarefree_part_primitive_bigint_coeffs(coeffs))
}

fn signed_remainder_sequence(a: &[BigInt], b: &[BigInt]) -> Vec<Vec<BigInt>> {
    let p0 = primitive_positive(trim_slice(a));
    if p0.is_empty() || degree(&p0).unwrap_or(0) == 0 {
        return vec![p0];
    }
    let p1 = primitive_keep_sign(trim_slice(b));
    if p1.is_empty() {
        return vec![p0];
    }

    let mut chain = vec![p0, p1];
    loop {
        let n = chain.len();
        let prem = pseudo_remainder_positive_multiplier(&chain[n - 2], &chain[n - 1]);
        if prem.is_empty() {
            break;
        }
        let next = primitive_keep_sign(poly_neg(&primitive_keep_sign(prem)));
        if next.is_empty() {
            break;
        }
        chain.push(next);
    }
    chain
}

/// Build the signed PRS of `p` and `p'` and recover the square-free degree
/// from its final nonzero member, which is the gcd up to a nonzero scalar.
///
/// Sturm variations are unchanged when the common gcd is divided out of every
/// member.  Consequently this one sequence both counts the distinct roots of
/// `p` and supplies `deg(p / gcd(p, p'))`, without first constructing the
/// square-free quotient and starting a second Euclidean pass.
fn derivative_sturm_chain(coeffs: &[BigInt]) -> (Vec<Vec<BigInt>>, usize) {
    let p = primitive_positive(trim_slice(coeffs));
    let Some(p_degree) = degree(&p) else {
        return (Vec::new(), 0);
    };
    if p_degree == 0 {
        return (vec![p], 0);
    }

    let derivative = poly_derivative(&p);
    let chain = signed_remainder_sequence(&p, &derivative);
    let gcd_degree = chain.last().and_then(|last| degree(last)).unwrap_or(0);
    (chain, p_degree.saturating_sub(gcd_degree))
}

fn sign_at_pos_infinity(p: &[BigInt]) -> i8 {
    sign_i8(&leading(p))
}

fn sign_at_neg_infinity(p: &[BigInt]) -> i8 {
    let d = match degree(p) {
        Some(d) => d,
        None => return 0,
    };
    let s = sign_i8(&leading(p));
    if d % 2 == 0 {
        s
    } else {
        -s
    }
}

fn sign_at_zero_plus(p: &[BigInt]) -> i8 {
    for c in p {
        let s = sign_i8(c);
        if s != 0 {
            return s;
        }
    }
    0
}

fn sign_variations<I>(signs: I) -> usize
where
    I: IntoIterator<Item = i8>,
{
    let mut prev = 0i8;
    let mut changes = 0usize;
    for s in signs {
        if s == 0 {
            continue;
        }
        if prev != 0 && prev != s {
            changes += 1;
        }
        prev = s;
    }
    changes
}

fn variations_at_neg_infinity(chain: &[Vec<BigInt>]) -> usize {
    sign_variations(chain.iter().map(|p| sign_at_neg_infinity(p)))
}

fn variations_at_pos_infinity(chain: &[Vec<BigInt>]) -> usize {
    sign_variations(chain.iter().map(|p| sign_at_pos_infinity(p)))
}

fn variations_at_zero_plus(chain: &[Vec<BigInt>]) -> usize {
    sign_variations(chain.iter().map(|p| sign_at_zero_plus(p)))
}

fn alternating_neg_argument(coeffs: &[BigInt]) -> Vec<BigInt> {
    trim(
        coeffs
            .iter()
            .enumerate()
            .map(|(i, c)| if i % 2 == 0 { c.clone() } else { -c })
            .collect(),
    )
}

fn all_one_sign(coeffs: &[BigInt]) -> bool {
    let mut seen_pos = false;
    let mut seen_neg = false;
    for c in coeffs {
        if c.is_positive() {
            seen_pos = true;
        } else if c.is_negative() {
            seen_neg = true;
        }
    }
    !(seen_pos && seen_neg)
}

fn make_nonnegative(mut coeffs: Vec<BigInt>) -> Vec<BigInt> {
    if coeffs.iter().any(|c| c.is_negative()) {
        for c in &mut coeffs {
            *c = -c.clone();
        }
    }
    coeffs
}

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
fn strict_fujiwara_root_bound(coeffs: &[BigInt]) -> BigInt {
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

/// Degree of the square-free part of a `BigInt` polynomial.
///
/// The zero polynomial and nonzero constants both return `0`.
pub fn squarefree_degree_bigint_coeffs(coeffs: &[BigInt]) -> usize {
    derivative_sturm_chain(coeffs).1
}

/// Count distinct real roots using a primitive pseudo-remainder Sturm sequence.
///
/// The input may have repeated roots. The signed PRS of `p` and `p'` ends in
/// their gcd; its endpoint variations count the distinct roots directly.
pub fn count_real_roots_prs_bigint_coeffs(coeffs: &[BigInt]) -> usize {
    let (chain, _) = derivative_sturm_chain(coeffs);
    if chain.is_empty() {
        return 0;
    }
    variations_at_neg_infinity(&chain).saturating_sub(variations_at_pos_infinity(&chain))
}

/// Count distinct positive roots using a primitive pseudo-remainder Sturm sequence.
///
/// The count is over the open interval `(0, +infinity)`.
pub fn count_positive_roots_prs_bigint_coeffs(coeffs: &[BigInt]) -> usize {
    let (chain, _) = derivative_sturm_chain(coeffs);
    count_positive_roots_from_sturm_chain(&chain)
}

fn count_positive_roots_from_sturm_chain(chain: &[Vec<BigInt>]) -> usize {
    if chain.is_empty() {
        return 0;
    }
    variations_at_zero_plus(chain).saturating_sub(variations_at_pos_infinity(chain))
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

/// Compute the Sturm--Tarski query of `query` at the real roots of `root_poly`.
///
/// This returns `# {r : query(r)>0} - # {r : query(r)<0}`, where `r` runs over
/// the distinct real roots of the square-free part of `root_poly`.  The
/// implementation uses the signed pseudo-remainder sequence of
/// `P` and `P' * query`, so it avoids isolating roots.
pub fn tarski_query_prs_bigint_coeffs(root_poly: &[BigInt], query: &[BigInt]) -> isize {
    let p = squarefree_part_primitive_bigint_coeffs(root_poly);
    if p.is_empty() || degree(&p).unwrap_or(0) == 0 || trim_slice(query).is_empty() {
        return 0;
    }
    let p_prime_query = poly_mul(&poly_derivative(&p), query);
    if p_prime_query.is_empty() {
        return 0;
    }
    let chain = signed_remainder_sequence(&p, &p_prime_query);
    let variations_neg = variations_at_neg_infinity(&chain);
    let variations_pos = variations_at_pos_infinity(&chain);
    variations_neg as isize - variations_pos as isize
}

/// Count how `test_poly` signs distribute over the distinct real roots of `root_poly`.
///
/// The returned counts are exact.  Roots are not isolated: the function uses
/// Sturm--Tarski queries for `test_poly` and `test_poly^2`.
pub fn count_root_signs_prs_bigint_coeffs(
    root_poly: &[BigInt],
    test_poly: &[BigInt],
) -> RootSignCounts {
    let total = count_real_roots_prs_bigint_coeffs(root_poly);
    if total == 0 {
        return RootSignCounts {
            positive: 0,
            negative: 0,
            zero: 0,
            total: 0,
        };
    }
    let query = trim_slice(test_poly);
    if query.is_empty() {
        return RootSignCounts {
            positive: 0,
            negative: 0,
            zero: total,
            total,
        };
    }

    let signed = tarski_query_prs_bigint_coeffs(root_poly, &query);
    let nonzero = tarski_query_prs_bigint_coeffs(root_poly, &poly_mul(&query, &query));
    debug_assert!(nonzero >= 0);
    debug_assert!(nonzero >= signed.abs());
    debug_assert_eq!((nonzero + signed) % 2, 0);
    debug_assert!(usize::try_from(nonzero).is_ok_and(|n| n <= total));

    let positive = usize::try_from((nonzero + signed) / 2).unwrap_or(0);
    let negative = usize::try_from((nonzero - signed) / 2).unwrap_or(0);
    let zero = total.saturating_sub(positive + negative);
    RootSignCounts {
        positive,
        negative,
        zero,
        total,
    }
}

/// Check Newton's inequalities for a nonnegative coefficient sequence.
///
/// This is a cheap necessary condition for real-rootedness.  It is exact and
/// works with arbitrary-size coefficients.
pub fn satisfies_newton_inequalities_bigint(coeffs: &[BigInt]) -> bool {
    let p = trim_slice(strip_initial_zeros_bigint(coeffs));
    let d = match degree(&p) {
        Some(d) if d >= 2 => d,
        _ => return true,
    };
    if p.iter().any(|c| c.is_negative()) {
        return false;
    }

    let mut binom = vec![BigInt::one(); d + 1];
    for k in 1..=d {
        binom[k] = &binom[k - 1] * BigInt::from(d - k + 1) / BigInt::from(k);
    }

    for k in 1..d {
        let lhs = p[k].pow(2) * &binom[k - 1] * &binom[k + 1];
        let rhs = &p[k - 1] * &p[k + 1] * &binom[k] * &binom[k];
        if lhs < rhs {
            return false;
        }
    }
    true
}

/// Check Kurtz's sufficient condition for distinct real roots.
///
/// If `f(t)=sum a_i t^i` has positive coefficients and
/// `a_i^2 > 4 a_{i-1} a_{i+1}` for every `1 <= i < d`, then all roots of `f`
/// are distinct and real.  Since all coefficients are positive, those roots
/// are automatically negative.  The condition is sufficient, not necessary.
///
/// This returns `false` for degree `0` and `1`, where the criterion has no
/// inequalities to check.
pub fn satisfies_kurtz_condition_bigint(coeffs: &[BigInt]) -> bool {
    let p = trim_slice(strip_initial_zeros_bigint(coeffs));
    let d = match degree(&p) {
        Some(d) if d >= 2 => d,
        _ => return false,
    };
    if p.iter().take(d + 1).any(|c| !c.is_positive()) {
        return false;
    }

    for i in 1..d {
        if p[i].pow(2) <= BigInt::from(4) * &p[i - 1] * &p[i + 1] {
            return false;
        }
    }
    true
}

/// Exact real-rootedness via primitive integer Sturm/PRS root counting.
///
/// This is a general method: it works for mixed-sign coefficients as well as
/// one-signed combinatorial polynomials.  It is often useful as a fallback when
/// Bézout matrices become large.
pub fn is_real_rooted_prs_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let (chain, squarefree_degree) = derivative_sturm_chain(coeffs);
    if squarefree_degree <= 1 {
        return true;
    }
    variations_at_neg_infinity(&chain).saturating_sub(variations_at_pos_infinity(&chain))
        == squarefree_degree
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

/// Whether the benchmarked Uspensky path is likely to beat primitive PRS.
///
/// Degree alone is not a reliable discriminator: PRS is substantially faster
/// on Narayana, Chebyshev, Hermite, and evenly spaced linear-factor products.
/// The families where Uspensky wins consistently in our benchmark (ordinary
/// and type-B Eulerian, and Touchard after removing its zero root) have equal
/// endpoint coefficients and interior coefficients at least `4^degree` times
/// as large. This exact comparison is invariant under scalar multiplication.
fn prefers_uspensky_for_one_signed(coeffs: &[BigInt]) -> bool {
    const MIN_USPENSKY_DEGREE: usize = 35;

    let d = match degree(coeffs) {
        Some(d) if d >= MIN_USPENSKY_DEGREE => d,
        _ => return false,
    };
    if coeffs[0] != coeffs[d] || !coeffs[0].is_positive() {
        return false;
    }

    let Some(growth_shift) = d.checked_mul(2) else {
        return false;
    };
    let large_interior_threshold = &coeffs[0] << growth_shift;
    coeffs[1..d]
        .iter()
        .any(|coefficient| coefficient >= &large_interior_threshold)
}

/// Exact real-rootedness optimized for one-signed coefficient polynomials.
///
/// Returns `None` if the nonzero coefficients do not all have the same sign.
/// Powers of `t` are removed first, since zero is already a real root.
pub fn is_real_rooted_one_signed_bigint_coeffs(coeffs: &[BigInt]) -> Option<bool> {
    let p = trim_slice(strip_initial_zeros_bigint(coeffs));
    let d = match degree(&p) {
        Some(d) => d,
        None => return Some(true),
    };
    if d <= 1 {
        return Some(true);
    }
    if !all_one_sign(&p) {
        return None;
    }

    let p = make_nonnegative(p);
    if satisfies_kurtz_condition_bigint(&p) {
        return Some(true);
    }
    if !satisfies_newton_inequalities_bigint(&p) {
        return Some(false);
    }

    let use_uspensky = prefers_uspensky_for_one_signed(&p);
    let (positive_roots, squarefree_degree) = if use_uspensky {
        let squarefree = squarefree_part_primitive_bigint_coeffs(&p);
        let squarefree_degree = degree(&squarefree).unwrap_or(0);
        let transformed = alternating_neg_argument(&squarefree);
        (
            count_positive_roots_uspensky_squarefree(&transformed),
            squarefree_degree,
        )
    } else {
        let transformed = alternating_neg_argument(&p);
        let (chain, squarefree_degree) = derivative_sturm_chain(&transformed);
        (
            count_positive_roots_from_sturm_chain(&chain),
            squarefree_degree,
        )
    };
    Some(positive_roots == squarefree_degree)
}

/// Exact real-rootedness with an adaptive one-signed path and PRS fallback.
///
/// The one-signed path first applies coefficient filters, then selects between
/// primitive PRS and Uspensky using a conservative benchmark-derived
/// heuristic. Mixed-sign inputs use PRS. This avoids constructing Bézout/PSD
/// matrices entirely and is the backend used by the public
/// `real_rootedness::is_real_rooted_bigint_coeffs` default.
pub fn is_real_rooted_fast_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    if let Some(rr) = is_real_rooted_one_signed_bigint_coeffs(coeffs) {
        rr
    } else {
        is_real_rooted_prs_bigint_coeffs(coeffs)
    }
}

/// Convenience wrapper for `i64` coefficient vectors.
pub fn is_real_rooted_fast_i64(coeffs: &[i64]) -> bool {
    let coeffs: Vec<BigInt> = coeffs.iter().map(|&c| BigInt::from(c)).collect();
    is_real_rooted_fast_bigint_coeffs(&coeffs)
}

/// The largest coefficient bit-size seen in a primitive Sturm/PRS chain.
///
/// This is useful for benchmarking and diagnosing coefficient swell.
pub fn primitive_sturm_max_coefficient_bits(coeffs: &[BigInt]) -> u64 {
    sturm_chain_squarefree(coeffs)
        .iter()
        .flat_map(|p| p.iter())
        .map(|c| c.bits())
        .max()
        .unwrap_or(0)
}

/// Try to reduce a BigInt vector to i64 coefficients.
pub fn bigint_coeffs_to_i64(coeffs: &[BigInt]) -> Option<Vec<i64>> {
    coeffs.iter().map(ToPrimitive::to_i64).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn b(v: &[i64]) -> Vec<BigInt> {
        v.iter().map(|&x| BigInt::from(x)).collect()
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
        assert!(satisfies_newton_inequalities_bigint(&b(&[1, 4, 6, 4, 1])));
        assert!(!satisfies_newton_inequalities_bigint(&b(&[1, 1, 10, 1])));
    }

    #[test]
    fn test_kurtz_condition_filter() {
        // (t + 1)(t + 10)(t + 100) = 1000 + 1110t + 111t^2 + t^3.
        let widely_spaced = b(&[1000, 1110, 111, 1]);
        assert!(satisfies_kurtz_condition_bigint(&widely_spaced));
        assert_eq!(
            is_real_rooted_one_signed_bigint_coeffs(&widely_spaced),
            Some(true)
        );

        // (1+t)^4 is real-rooted, but the Kurtz condition is intentionally
        // much stronger than real-rootedness.
        assert!(!satisfies_kurtz_condition_bigint(&b(&[1, 4, 6, 4, 1])));
        assert!(!satisfies_kurtz_condition_bigint(&b(&[1, 0, 1])));
        assert!(!satisfies_kurtz_condition_bigint(&b(&[
            -1000, -1110, -111, -1
        ])));
        assert!(!satisfies_kurtz_condition_bigint(&b(&[1, 10])));
    }

    #[test]
    fn test_prs_real_rootedness_general() {
        assert!(is_real_rooted_prs_bigint_coeffs(&b(&[-6, 11, -6, 1])));
        assert!(!is_real_rooted_prs_bigint_coeffs(&b(&[1, 0, 1])));
        assert!(is_real_rooted_prs_bigint_coeffs(&b(&[0, 0, 1, 2, 1]))); // x^2(1+x)^2
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
                let squarefree = squarefree_part_primitive_bigint_coeffs(&coeffs);
                let reference_squarefree_degree = super::degree(&squarefree).unwrap_or(0);
                let reference = sturm_chain_from_squarefree(squarefree);

                assert_eq!(direct_squarefree_degree, reference_squarefree_degree);
                assert_eq!(
                    variations_at_neg_infinity(&direct)
                        .saturating_sub(variations_at_pos_infinity(&direct)),
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
    fn test_adaptive_uspensky_selector_on_benchmarked_families() {
        let eulerian = crate::sequences::eulerian_polynomials_bigint(36)
            .pop()
            .unwrap();
        let eulerian_below_cutoff = crate::sequences::eulerian_polynomials_bigint(35)
            .pop()
            .unwrap();
        let narayana = crate::sequences::narayana_polynomials_bigint(36)
            .pop()
            .unwrap();
        let type_b_eulerian = crate::sequences::type_b_eulerian_polynomials_bigint(35)
            .pop()
            .unwrap();

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

        assert!(prefers_uspensky_for_one_signed(&eulerian));
        assert!(!prefers_uspensky_for_one_signed(&eulerian_below_cutoff));
        assert!(prefers_uspensky_for_one_signed(&type_b_eulerian));
        assert!(prefers_uspensky_for_one_signed(&touchard));
        assert!(!prefers_uspensky_for_one_signed(&narayana));
        assert!(!prefers_uspensky_for_one_signed(&product));

        let scaled_eulerian = poly_scale(&eulerian, &BigInt::from(7));
        assert!(prefers_uspensky_for_one_signed(&scaled_eulerian));

        assert!(is_real_rooted_fast_bigint_coeffs(&eulerian));
        assert!(is_real_rooted_fast_bigint_coeffs(&type_b_eulerian));
        assert!(is_real_rooted_fast_bigint_coeffs(&touchard));
    }

    #[test]
    fn test_prs_agrees_with_rational_sturm_on_small_grid() {
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
                    is_real_rooted_prs_bigint_coeffs(&bigint),
                    crate::real_rootedness::is_real_rooted_sturm_bigint_coeffs(&bigint),
                    "coeffs={coeffs:?}"
                );
            }
        }
    }
}
