//! Signed pseudo-remainder sequences and exact Sturm--Tarski root counts.
//! Uses the integer arithmetic helpers and shares Sturm chains with isolation.

#[cfg(test)]
use crate::arithmetic::trim;
use crate::arithmetic::{
    degree, leading, poly_derivative, poly_mul, poly_neg, primitive_keep_sign, primitive_positive,
    pseudo_remainder_positive_multiplier, pseudo_remainder_positive_multiplier_into, sign_i8,
    squarefree_part_primitive_bigint_coeffs, trim_slice,
};
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

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

pub(super) fn sturm_chain_from_squarefree(p0: Vec<BigInt>) -> Vec<Vec<BigInt>> {
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

pub(super) fn signed_remainder_sequence(a: &[BigInt], b: &[BigInt]) -> Vec<Vec<BigInt>> {
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

/// Build the same signed pseudo-remainder sequence without removing content.
///
/// This avoids coefficient GCDs at small fixed degrees, where the chain is
/// short enough that coefficient swell costs less than repeated normalization.
#[cfg(test)]
fn signed_remainder_sequence_unreduced(a: &[BigInt], b: &[BigInt]) -> Vec<Vec<BigInt>> {
    // The derivative driver already supplies a primitive, positive-leading
    // first member. Avoid repeating its coefficient-content GCD here.
    let p0 = trim_slice(a);
    debug_assert!(p0.is_empty() || leading(&p0).is_positive());
    if p0.is_empty() || degree(&p0).unwrap_or(0) == 0 {
        return vec![p0];
    }
    let p1 = trim_slice(b);
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
        chain.push(poly_neg(&prem));
    }
    chain
}

/// Build a signed subresultant PRS using the Collins--Brown scaling rule.
///
/// Each pseudo-remainder has a predictable factor derived from earlier
/// leading coefficients and degree drops. Dividing by that factor is much
/// cheaper than computing the full coefficient content at every step. We use
/// its positive absolute value and negate the remainder, preserving the
/// Sturm sign convention.
#[cfg(test)]
pub(super) fn signed_remainder_sequence_subresultant(
    a: &[BigInt],
    b: &[BigInt],
) -> Vec<Vec<BigInt>> {
    // The derivative driver already supplies a primitive, positive-leading
    // first member. Avoid repeating its coefficient-content GCD here.
    let p0 = trim_slice(a);
    debug_assert!(p0.is_empty() || leading(&p0).is_positive());
    if p0.is_empty() || degree(&p0).unwrap_or(0) == 0 {
        return vec![p0];
    }
    let p1 = primitive_keep_sign(trim_slice(b));
    if p1.is_empty() {
        return vec![p0];
    }

    let mut chain = vec![p0, p1];
    let mut g = BigInt::one();
    let mut h = BigInt::one();

    loop {
        let n = chain.len();
        let previous_degree = degree(&chain[n - 2]).expect("nonzero PRS member");
        let current_degree = degree(&chain[n - 1]).expect("nonzero PRS member");
        let delta = previous_degree - current_degree;
        let prem = pseudo_remainder_positive_multiplier(&chain[n - 2], &chain[n - 1]);
        if prem.is_empty() {
            break;
        }

        let divisor = &g * h.pow(delta as u32);
        debug_assert!(divisor.is_positive());
        let mut next = Vec::with_capacity(prem.len());
        for coefficient in prem {
            debug_assert!((&coefficient % &divisor).is_zero());
            next.push(-(coefficient / &divisor));
        }
        let next = trim(next);
        if next.is_empty() {
            break;
        }

        g = leading(&chain[n - 1]).abs();
        h = if delta == 1 {
            g.clone()
        } else {
            let numerator = g.pow(delta as u32);
            let denominator = h.pow((delta - 1) as u32);
            debug_assert!((&numerator % &denominator).is_zero());
            numerator / denominator
        };
        debug_assert!(h.is_positive());
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
#[cfg(test)]
pub(super) fn derivative_sturm_chain(coeffs: &[BigInt]) -> (Vec<Vec<BigInt>>, usize) {
    let p = primitive_positive(trim_slice(coeffs));
    let Some(p_degree) = degree(&p) else {
        return (Vec::new(), 0);
    };
    if p_degree == 0 {
        return (vec![p], 0);
    }

    let derivative = poly_derivative(&p);
    // Very short chains do not repay the scaling bookkeeping. From degree
    // six onward, predictable subresultant divisions avoid both coefficient
    // swell and a coefficient-content GCD at every step.
    let chain = if (2..=5).contains(&p_degree) {
        signed_remainder_sequence_unreduced(&p, &derivative)
    } else {
        signed_remainder_sequence_subresultant(&p, &derivative)
    };
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

pub(super) fn sign_variations<I>(signs: I) -> usize
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

pub(super) fn variations_at_neg_infinity(chain: &[Vec<BigInt>]) -> usize {
    sign_variations(chain.iter().map(|p| sign_at_neg_infinity(p)))
}

pub(super) fn variations_at_pos_infinity(chain: &[Vec<BigInt>]) -> usize {
    sign_variations(chain.iter().map(|p| sign_at_pos_infinity(p)))
}

#[cfg(test)]
pub(super) fn variations_at_zero_plus(chain: &[Vec<BigInt>]) -> usize {
    sign_variations(chain.iter().map(|p| sign_at_zero_plus(p)))
}

#[derive(Clone, Copy, Debug, Default)]
struct RunningVariations {
    previous: i8,
    changes: usize,
}

impl RunningVariations {
    fn push(&mut self, sign: i8) {
        if sign == 0 {
            return;
        }
        if self.previous != 0 && self.previous != sign {
            self.changes += 1;
        }
        self.previous = sign;
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(super) struct DerivativeSturmSummary {
    pub(super) real_roots: usize,
    pub(super) positive_roots: usize,
    pub(super) squarefree_degree: usize,
}

fn observe_sturm_member(
    polynomial: &[BigInt],
    count_real_roots: bool,
    count_positive_roots: bool,
    negative_infinity: &mut RunningVariations,
    positive_infinity: &mut RunningVariations,
    zero_plus: &mut RunningVariations,
) {
    if count_real_roots {
        negative_infinity.push(sign_at_neg_infinity(polynomial));
    }
    if count_real_roots || count_positive_roots {
        positive_infinity.push(sign_at_pos_infinity(polynomial));
    }
    if count_positive_roots {
        zero_plus.push(sign_at_zero_plus(polynomial));
    }
}

fn negate_in_place(polynomial: &mut [BigInt]) {
    for coefficient in polynomial {
        *coefficient = -std::mem::take(coefficient);
    }
}

/// Stream the signed PRS of `p` and `p'` into the endpoint data Sturm needs.
///
/// Only two consecutive sequence members and one reusable remainder buffer are
/// retained. The buffer becomes the polynomial from two steps earlier after
/// each rotation, so both its coefficient vector and its `BigInt` limb storage
/// can be reused by the next pseudo-division.
pub(super) fn derivative_sturm_summary(
    coeffs: &[BigInt],
    count_real_roots: bool,
    count_positive_roots: bool,
) -> DerivativeSturmSummary {
    let mut previous = primitive_positive(trim_slice(coeffs));
    let Some(polynomial_degree) = degree(&previous) else {
        return DerivativeSturmSummary::default();
    };
    if polynomial_degree == 0 {
        return DerivativeSturmSummary::default();
    }

    let use_subresultant_scaling = polynomial_degree >= 6;
    let derivative = poly_derivative(&previous);
    let mut current = if use_subresultant_scaling {
        primitive_keep_sign(derivative)
    } else {
        derivative
    };
    debug_assert!(!current.is_empty());

    let mut negative_infinity = RunningVariations::default();
    let mut positive_infinity = RunningVariations::default();
    let mut zero_plus = RunningVariations::default();
    observe_sturm_member(
        &previous,
        count_real_roots,
        count_positive_roots,
        &mut negative_infinity,
        &mut positive_infinity,
        &mut zero_plus,
    );
    observe_sturm_member(
        &current,
        count_real_roots,
        count_positive_roots,
        &mut negative_infinity,
        &mut positive_infinity,
        &mut zero_plus,
    );

    let mut scratch = Vec::new();
    let mut g = BigInt::one();
    let mut h = BigInt::one();

    loop {
        let previous_degree = degree(&previous).expect("nonzero PRS member");
        let current_degree = degree(&current).expect("nonzero PRS member");
        let delta = previous_degree - current_degree;
        pseudo_remainder_positive_multiplier_into(&mut scratch, &previous, &current);
        if scratch.is_empty() {
            break;
        }

        if use_subresultant_scaling {
            let divisor = &g * h.pow(delta as u32);
            debug_assert!(divisor.is_positive());
            for coefficient in &mut scratch {
                debug_assert!((&*coefficient % &divisor).is_zero());
                *coefficient /= &divisor;
            }

            g = leading(&current).abs();
            h = if delta == 1 {
                g.clone()
            } else {
                let numerator = g.pow(delta as u32);
                let denominator = h.pow((delta - 1) as u32);
                debug_assert!((&numerator % &denominator).is_zero());
                numerator / denominator
            };
            debug_assert!(h.is_positive());
        }
        negate_in_place(&mut scratch);

        observe_sturm_member(
            &scratch,
            count_real_roots,
            count_positive_roots,
            &mut negative_infinity,
            &mut positive_infinity,
            &mut zero_plus,
        );
        std::mem::swap(&mut previous, &mut current);
        std::mem::swap(&mut current, &mut scratch);
    }

    let gcd_degree = degree(&current).unwrap_or(0);
    DerivativeSturmSummary {
        real_roots: negative_infinity
            .changes
            .saturating_sub(positive_infinity.changes),
        positive_roots: zero_plus.changes.saturating_sub(positive_infinity.changes),
        squarefree_degree: polynomial_degree.saturating_sub(gcd_degree),
    }
}

/// Degree of the square-free part of a `BigInt` polynomial.
///
/// The zero polynomial and nonzero constants both return `0`.
pub fn squarefree_degree_bigint_coeffs(coeffs: &[BigInt]) -> usize {
    derivative_sturm_summary(coeffs, false, false).squarefree_degree
}

/// Count distinct real roots using a signed pseudo-remainder Sturm sequence.
///
/// The input may have repeated roots. The signed PRS of `p` and `p'` ends in
/// their gcd; its endpoint variations count the distinct roots directly.
pub fn count_real_roots_prs_bigint_coeffs(coeffs: &[BigInt]) -> usize {
    derivative_sturm_summary(coeffs, true, false).real_roots
}

/// Count distinct positive roots using a signed pseudo-remainder Sturm sequence.
///
/// The count is over the open interval `(0, +infinity)`.
pub fn count_positive_roots_prs_bigint_coeffs(coeffs: &[BigInt]) -> usize {
    derivative_sturm_summary(coeffs, false, true).positive_roots
}

#[cfg(test)]
pub(super) fn count_positive_roots_from_sturm_chain(chain: &[Vec<BigInt>]) -> usize {
    if chain.is_empty() {
        return 0;
    }
    variations_at_zero_plus(chain).saturating_sub(variations_at_pos_infinity(chain))
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

/// Exact real-rootedness via integer Sturm/PRS root counting.
///
/// This is a general method: it works for mixed-sign coefficients as well as
/// one-signed combinatorial polynomials.  It is often useful as a fallback when
/// Bézout matrices become large.
pub fn is_real_rooted_prs_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let summary = derivative_sturm_summary(coeffs, true, false);
    if summary.squarefree_degree <= 1 {
        return true;
    }
    summary.real_roots == summary.squarefree_degree
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
