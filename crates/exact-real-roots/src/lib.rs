//! Exact real-root counting over the integers.
//!
//! This module provides the default exact backend for boolean
//! real-rootedness. Degrees at most four use exact discriminant and quartic
//! invariant formulas. Degree five uses an unreduced signed pseudo-remainder
//! sequence; from degree six, a Collins--Brown subresultant PRS removes known
//! exact factors without computing coefficient-content GCDs. One-signed
//! inputs first try cheap coefficient filters. These paths avoid the large
//! exact PSD matrices used by the Bézout criterion.
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
//!
//! # Entry points
//!
//! - Decide real-rootedness: [`is_real_rooted_fast_bigint_coeffs`] (adaptive
//!   default, also [`is_real_rooted_fast_i64`]);
//!   [`is_real_rooted_prs_bigint_coeffs`] and
//!   [`is_real_rooted_uspensky_bigint_coeffs`] force one backend.
//! - Count distinct real roots: [`count_real_roots_prs_bigint_coeffs`],
//!   [`count_positive_roots_prs_bigint_coeffs`] and the Uspensky counterparts;
//!   [`squarefree_degree_bigint_coeffs`] gives the square-free degree to compare
//!   against.
//! - Isolate roots: [`isolate_real_roots_bigint_coeffs`] returns sorted
//!   [`IsolatedRealRoot`] values with exact dyadic endpoints.
//! - Signs at roots: [`count_root_signs_prs_bigint_coeffs`] and
//!   [`tarski_query_prs_bigint_coeffs`].
//! - Coefficient-only filters: [`satisfies_newton_inequalities_bigint_coeffs`]
//!   (necessary) and [`satisfies_kurtz_condition_bigint_coeffs`] (sufficient).
//!
//! # Conventions
//!
//! Polynomials are `&[BigInt]` (or `&[i64]` in the convenience wrappers) in
//! ascending degree order: `coeffs[i]` is the coefficient of `t^i`.  Trailing
//! zeros are ignored, and the zero polynomial and nonzero constants have no
//! roots and count as real-rooted.  Repeated roots are counted once, so counts
//! are numbers of *distinct* roots.  All results are exact.

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

/// Pseudo-remainder with a positive multiplier.
///
/// If `b_pos` is `b` multiplied by a sign so that its leading coefficient is
/// positive, this returns `lc(b_pos)^k * rem(a, b_pos)` for a nonnegative `k`.
/// The scalar is positive, so it does not change Sturm sign variations.
fn pseudo_remainder_positive_multiplier(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let mut output = Vec::new();
    pseudo_remainder_positive_multiplier_into(&mut output, a, b);
    output
}

/// In-place-buffer variant of [`pseudo_remainder_positive_multiplier`].
///
/// Existing coefficient slots are updated with `clone_from`, allowing their
/// limb allocations to be reused by a streaming remainder sequence.
fn pseudo_remainder_positive_multiplier_into(output: &mut Vec<BigInt>, a: &[BigInt], b: &[BigInt]) {
    let da = match degree(a) {
        Some(d) => d,
        None => {
            output.clear();
            return;
        }
    };
    let db = degree(b).expect("pseudo_remainder: division by zero polynomial");
    if da < db {
        output.resize_with(a.len(), BigInt::zero);
        for (target, source) in output.iter_mut().zip(a) {
            target.clone_from(source);
        }
        output.truncate(da + 1);
        return;
    }

    let b = &b[..=db];
    let divisor_is_negative = b[db].is_negative();
    let lc = b[db].abs();
    debug_assert!(lc.is_positive());

    output.resize_with(da + 1, BigInt::zero);
    for (target, source) in output.iter_mut().zip(&a[..=da]) {
        target.clone_from(source);
    }
    output.truncate(da + 1);

    let mut e = da - db + 1;
    while let Some(dr) = degree(output) {
        if dr < db {
            break;
        }
        // The leading term cancels exactly in
        // `lc(b) * output - lc(output) * b`. Move its coefficient out rather
        // than cloning it, and leave that slot zero instead of performing the
        // two large-integer operations whose results would cancel.
        let c = std::mem::take(&mut output[dr]);
        let shift = dr - db;
        for coefficient in &mut output[..dr] {
            *coefficient *= &lc;
        }
        for (j, coefficient) in b[..db].iter().enumerate() {
            if !coefficient.is_zero() {
                if divisor_is_negative {
                    output[j + shift] += &c * coefficient;
                } else {
                    output[j + shift] -= &c * coefficient;
                }
            }
        }
        while output.last().is_some_and(BigInt::is_zero) {
            output.pop();
        }
        e -= 1;
    }
    if e > 0 && !output.is_empty() {
        let remaining_scale = lc.pow(e as u32);
        for coefficient in output.iter_mut() {
            *coefficient *= &remaining_scale;
        }
    }
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
/// This lower-level entry point is exposed so polynomial libraries can reuse
/// the same fraction-free integer GCD as the root-counting backends.
#[doc(hidden)]
pub fn squarefree_part_primitive_bigint_coeffs(coeffs: &[BigInt]) -> Vec<BigInt> {
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
fn signed_remainder_sequence_subresultant(a: &[BigInt], b: &[BigInt]) -> Vec<Vec<BigInt>> {
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
fn derivative_sturm_chain(coeffs: &[BigInt]) -> (Vec<Vec<BigInt>>, usize) {
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

#[cfg(test)]
fn variations_at_zero_plus(chain: &[Vec<BigInt>]) -> usize {
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
struct DerivativeSturmSummary {
    real_roots: usize,
    positive_roots: usize,
    squarefree_degree: usize,
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
fn derivative_sturm_summary(
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

/// Return whether a polynomial of degree at most four has only real roots.
///
/// The discriminant formulas below are evaluated only up to a positive
/// factor, since only their sign is needed.  The quartic invariants also cover
/// zero-discriminant cases, so repeated roots do not require a GCD fallback.
fn low_degree_is_real_rooted(coeffs: &[BigInt]) -> Option<bool> {
    let p = trim_slice(coeffs);
    let Some(d) = degree(&p) else {
        return Some(true);
    };

    match d {
        0 | 1 => Some(true),
        2 => {
            // b^2 - 4ac >= 0.
            let discriminant = &p[1] * &p[1] - BigInt::from(4) * &p[2] * &p[0];
            Some(!discriminant.is_negative())
        }
        3 => {
            // For ax^3 + bx^2 + cx + d, the usual discriminant has the
            // same sign as 4 Delta_0^3 - Delta_1^2, since their ratio is
            // the positive quantity 27a^2.
            let delta_0 = &p[2] * &p[2] - BigInt::from(3) * &p[3] * &p[1];
            let delta_1 = BigInt::from(2) * &p[2] * &p[2] * &p[2]
                - BigInt::from(9) * &p[3] * &p[2] * &p[1]
                + BigInt::from(27) * &p[3] * &p[3] * &p[0];
            let discriminant_numerator =
                BigInt::from(4) * &delta_0 * &delta_0 * &delta_0 - &delta_1 * &delta_1;
            Some(!discriminant_numerator.is_negative())
        }
        4 => {
            // Binary-quartic invariants.  The usual discriminant has the
            // same sign as 4 I^3 - J^2 (the ratio is the positive constant
            // 27).  When it is positive, the quartic has four real roots
            // exactly when P < 0 and D < 0.
            let i =
                BigInt::from(12) * &p[4] * &p[0] - BigInt::from(3) * &p[3] * &p[1] + &p[2] * &p[2];
            let j = BigInt::from(72) * &p[4] * &p[2] * &p[0]
                + BigInt::from(9) * &p[3] * &p[2] * &p[1]
                - BigInt::from(27) * &p[4] * &p[1] * &p[1]
                - BigInt::from(27) * &p[3] * &p[3] * &p[0]
                - BigInt::from(2) * &p[2] * &p[2] * &p[2];
            let discriminant_numerator = BigInt::from(4) * &i * &i * &i - &j * &j;

            if discriminant_numerator.is_negative() {
                return Some(false);
            }

            let quartic_p = BigInt::from(8) * &p[4] * &p[2] - BigInt::from(3) * &p[3] * &p[3];
            if quartic_p.is_positive()
                || (discriminant_numerator.is_positive() && quartic_p.is_zero())
            {
                return Some(false);
            }

            let quartic_d = BigInt::from(64) * &p[4] * &p[4] * &p[4] * &p[0]
                - BigInt::from(16) * &p[4] * &p[4] * &p[2] * &p[2]
                + BigInt::from(16) * &p[4] * &p[3] * &p[3] * &p[2]
                - BigInt::from(16) * &p[4] * &p[4] * &p[3] * &p[1]
                - BigInt::from(3) * &p[3] * &p[3] * &p[3] * &p[3];

            if discriminant_numerator.is_positive() {
                Some(quartic_d.is_negative())
            } else {
                // Strict inequalities give one double and two simple real
                // roots. Equality includes two double roots and the triple-
                // or quadruple-root cases.
                Some(!quartic_p.is_positive() && !quartic_d.is_positive())
            }
        }
        _ => None,
    }
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
pub fn satisfies_newton_inequalities_bigint_coeffs(coeffs: &[BigInt]) -> bool {
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

/// Deprecated name of [`satisfies_newton_inequalities_bigint_coeffs`].
#[deprecated(
    since = "0.2.0",
    note = "renamed to `satisfies_newton_inequalities_bigint_coeffs`"
)]
pub fn satisfies_newton_inequalities_bigint(coeffs: &[BigInt]) -> bool {
    satisfies_newton_inequalities_bigint_coeffs(coeffs)
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
pub fn satisfies_kurtz_condition_bigint_coeffs(coeffs: &[BigInt]) -> bool {
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

/// Deprecated name of [`satisfies_kurtz_condition_bigint_coeffs`].
#[deprecated(
    since = "0.2.0",
    note = "renamed to `satisfies_kurtz_condition_bigint_coeffs`"
)]
pub fn satisfies_kurtz_condition_bigint(coeffs: &[BigInt]) -> bool {
    satisfies_kurtz_condition_bigint_coeffs(coeffs)
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

/// Exact real-rootedness optimized for one-signed coefficient polynomials.
///
/// Returns `None` if the nonzero coefficients do not all have the same sign.
/// Powers of `t` are removed first, since zero is already a real root.  After
/// that, polynomials of degree at most 1 always give `Some(true)`, whatever the
/// signs.
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
    if let Some(real_rooted) = low_degree_is_real_rooted(&p) {
        return Some(real_rooted);
    }

    let p = make_nonnegative(p);
    if satisfies_kurtz_condition_bigint_coeffs(&p) {
        return Some(true);
    }
    if !satisfies_newton_inequalities_bigint_coeffs(&p) {
        return Some(false);
    }

    // After stripping zero roots, a one-signed polynomial has only negative
    // real roots, so counting all distinct real roots (signs at +-infinity
    // only) is equivalent to, and cheaper than, counting positive roots of
    // p(-t). The subresultant PRS beats the Uspensky/Descartes path on every
    // benchmarked one-signed family (Eulerian, type B Eulerian, Touchard,
    // Narayana, linear-factor products; see `polytool bench real-rootedness`).
    let summary = derivative_sturm_summary(&p, true, false);
    Some(summary.real_roots == summary.squarefree_degree)
}

/// Exact real-rootedness with an adaptive one-signed path and PRS fallback.
///
/// The one-signed path first applies the Kurtz (sufficient) and Newton
/// (necessary) coefficient filters, then counts distinct real roots with the
/// subresultant PRS; mixed-sign inputs use PRS directly. (An earlier
/// heuristic switched large palindromic inputs to Uspensky; after the
/// Collins--Brown PRS it was about 10--26 times slower on those families.) This avoids constructing Bézout/PSD
/// matrices entirely and is the backend used by the public
/// `real_rootedness::is_real_rooted_bigint_coeffs` default.
pub fn is_real_rooted_fast_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    if let Some(real_rooted) = low_degree_is_real_rooted(coeffs) {
        return real_rooted;
    }
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

/// An exact dyadic rational `numerator / 2^exponent`.
///
/// Values are kept in lowest terms: either `exponent == 0` or the numerator
/// is odd. Root isolation uses dyadic endpoints because bisection stays in
/// this set and signs can be evaluated with integer arithmetic only.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct DyadicRational {
    numerator: BigInt,
    exponent: u32,
}

impl DyadicRational {
    /// Construct `numerator / 2^exponent`, reduced to lowest terms.
    pub fn new(numerator: BigInt, exponent: u32) -> Self {
        let mut numerator = numerator;
        let mut exponent = exponent;
        if numerator.is_zero() {
            exponent = 0;
        }
        while exponent > 0 && numerator.is_even() {
            numerator /= 2;
            exponent -= 1;
        }
        DyadicRational {
            numerator,
            exponent,
        }
    }

    /// The integer `value` as a dyadic rational.
    pub fn from_integer(value: BigInt) -> Self {
        DyadicRational {
            numerator: value,
            exponent: 0,
        }
    }

    /// The reduced numerator.
    pub fn numerator(&self) -> &BigInt {
        &self.numerator
    }

    /// The reduced power of two in the denominator.
    pub fn exponent(&self) -> u32 {
        self.exponent
    }

    /// The denominator `2^exponent`.
    pub fn denominator(&self) -> BigInt {
        BigInt::one() << self.exponent
    }

    /// Numerator after rescaling to the denominator `2^exponent`, which must
    /// be at least `self.exponent()`.
    fn numerator_at_exponent(&self, exponent: u32) -> BigInt {
        debug_assert!(exponent >= self.exponent);
        &self.numerator << (exponent - self.exponent)
    }

    /// The exact midpoint `(self + other) / 2`.
    pub fn midpoint(&self, other: &Self) -> Self {
        let exponent = self.exponent.max(other.exponent);
        let sum = self.numerator_at_exponent(exponent) + other.numerator_at_exponent(exponent);
        DyadicRational::new(sum, exponent + 1)
    }

    /// `self + numerator / 2^exponent`.
    fn add_dyadic(&self, numerator: &BigInt, exponent: u32) -> Self {
        let common = self.exponent.max(exponent);
        let other = numerator << (common - exponent);
        DyadicRational::new(self.numerator_at_exponent(common) + other, common)
    }

    /// Sign of `polynomial(self)` for ascending integer coefficients.
    ///
    /// The value is multiplied by the positive factor `2^(exponent * degree)`,
    /// so only integer arithmetic is needed.
    fn sign_of(&self, polynomial: &[BigInt]) -> i8 {
        let Some(d) = degree(polynomial) else {
            return 0;
        };
        let denominator = self.denominator();
        let mut accumulator = polynomial[d].clone();
        let mut denominator_power = BigInt::one();
        for coefficient in polynomial[..d].iter().rev() {
            denominator_power *= &denominator;
            accumulator = accumulator * &self.numerator + coefficient * &denominator_power;
        }
        sign_i8(&accumulator)
    }
}

impl PartialOrd for DyadicRational {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for DyadicRational {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let exponent = self.exponent.max(other.exponent);
        self.numerator_at_exponent(exponent)
            .cmp(&other.numerator_at_exponent(exponent))
    }
}

impl std::fmt::Display for DyadicRational {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.exponent == 0 {
            write!(f, "{}", self.numerator)
        } else {
            write!(f, "{}/2^{}", self.numerator, self.exponent)
        }
    }
}

/// One isolated distinct real root of an integer polynomial.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IsolatedRealRoot {
    /// The root is exactly this dyadic rational.
    Exact(DyadicRational),
    /// Exactly one distinct root lies in the open interval `(lower, upper)`;
    /// neither endpoint is a root.
    Interval {
        /// Lower endpoint (not a root).
        lower: DyadicRational,
        /// Upper endpoint (not a root).
        upper: DyadicRational,
    },
}

impl IsolatedRealRoot {
    /// Smallest point of the closed isolating set.
    pub fn lower(&self) -> &DyadicRational {
        match self {
            IsolatedRealRoot::Exact(root) => root,
            IsolatedRealRoot::Interval { lower, .. } => lower,
        }
    }

    /// Largest point of the closed isolating set.
    pub fn upper(&self) -> &DyadicRational {
        match self {
            IsolatedRealRoot::Exact(root) => root,
            IsolatedRealRoot::Interval { upper, .. } => upper,
        }
    }
}

fn sturm_variations_at(chain: &[Vec<BigInt>], point: &DyadicRational) -> usize {
    sign_variations(chain.iter().map(|p| point.sign_of(p)))
}

enum IsolationWork {
    Interval {
        lower: DyadicRational,
        upper: DyadicRational,
        lower_variations: usize,
        upper_variations: usize,
    },
    Emit(IsolatedRealRoot),
}

/// Isolate the distinct real roots of an integer polynomial.
///
/// Coefficients are in ascending degree order. The result has one entry per
/// distinct real root, sorted increasingly. Each entry is either an exact
/// dyadic root or an open interval with dyadic endpoints containing exactly
/// one root; interval endpoints are never roots. Consecutive entries satisfy
/// `previous.upper() <= next.lower()`, with strict inequality whenever one
/// of them is [`IsolatedRealRoot::Exact`], so neighbouring intervals may share
/// a (non-root) endpoint but never overlap.
///
/// The method is Sturm bisection on the primitive square-free part, starting
/// from a strict Fujiwara bound, with exact integer sign evaluation. Constant
/// polynomials, including the zero polynomial, return an empty vector: the
/// caller must treat the zero polynomial separately.
pub fn isolate_real_roots_bigint_coeffs(coeffs: &[BigInt]) -> Vec<IsolatedRealRoot> {
    let p = squarefree_part_primitive_bigint_coeffs(coeffs);
    if degree(&p).unwrap_or(0) == 0 {
        return Vec::new();
    }
    let bound = strict_fujiwara_root_bound(&p);
    let chain = sturm_chain_from_squarefree(p);

    let lower = DyadicRational::from_integer(-bound.clone());
    let upper = DyadicRational::from_integer(bound);
    let mut stack = vec![IsolationWork::Interval {
        lower_variations: sturm_variations_at(&chain, &lower),
        upper_variations: sturm_variations_at(&chain, &upper),
        lower,
        upper,
    }];
    let mut roots = Vec::new();

    while let Some(work) = stack.pop() {
        let (lower, upper, lower_variations, upper_variations) = match work {
            IsolationWork::Emit(root) => {
                roots.push(root);
                continue;
            }
            IsolationWork::Interval {
                lower,
                upper,
                lower_variations,
                upper_variations,
            } => (lower, upper, lower_variations, upper_variations),
        };
        let count = lower_variations.saturating_sub(upper_variations);
        if count == 0 {
            continue;
        }
        if count == 1 {
            roots.push(IsolatedRealRoot::Interval { lower, upper });
            continue;
        }

        let middle = lower.midpoint(&upper);
        if middle.sign_of(&chain[0]) != 0 {
            let middle_variations = sturm_variations_at(&chain, &middle);
            stack.push(IsolationWork::Interval {
                lower: middle.clone(),
                upper,
                lower_variations: middle_variations,
                upper_variations,
            });
            stack.push(IsolationWork::Interval {
                lower,
                upper: middle,
                lower_variations,
                upper_variations: middle_variations,
            });
            continue;
        }

        // The midpoint is a root: shrink a symmetric window around it until it
        // isolates that root, then continue on both sides of the window.
        let width_exponent = lower.exponent().max(upper.exponent());
        let width = upper.numerator_at_exponent(width_exponent)
            - lower.numerator_at_exponent(width_exponent);
        let mut half_window_exponent = width_exponent + 2;
        let (left, right, left_variations, right_variations) = loop {
            let left = middle.add_dyadic(&-&width, half_window_exponent);
            let right = middle.add_dyadic(&width, half_window_exponent);
            if left.sign_of(&chain[0]) != 0 && right.sign_of(&chain[0]) != 0 {
                let left_variations = sturm_variations_at(&chain, &left);
                let right_variations = sturm_variations_at(&chain, &right);
                if left_variations.saturating_sub(right_variations) == 1 {
                    break (left, right, left_variations, right_variations);
                }
            }
            half_window_exponent += 1;
        };
        stack.push(IsolationWork::Interval {
            lower: right,
            upper,
            lower_variations: right_variations,
            upper_variations,
        });
        stack.push(IsolationWork::Emit(IsolatedRealRoot::Exact(middle)));
        stack.push(IsolationWork::Interval {
            lower,
            upper: left,
            lower_variations,
            upper_variations: left_variations,
        });
    }
    roots
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::BigRational;

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
                let reference_squarefree_degree = super::degree(&squarefree).unwrap_or(0);
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
}
