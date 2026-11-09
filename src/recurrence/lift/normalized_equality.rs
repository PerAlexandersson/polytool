//! All-index equality of two recurrence-defined sequences up to a
//! normalization: an index shift, dropped leading rows, a power of `x`
//! with an affine exponent and a hypergeometric per-row scalar.
//!
//! [`certify_normalized_recurrence_equality`] proves, for every index
//! `n >= n1`, that the rows of a *target* record are
//!
//! ```text
//! T_n = c(n) x^(alpha n + beta) S_(n + s),     n1 = (target first index) + d,
//! ```
//!
//! with `S` the rows of a *source* record, or reports why it could not.
//! The normalization ([`RecurrenceNormalization`]) is applied in a fixed
//! order: the index shift `s`, the `d` dropped leading target rows, the
//! power `x^(alpha n + beta)` (integers `alpha`, `beta`; a negative
//! exponent is an exact division) and the scalar `c(n)` with `c(n1)` a
//! nonzero rational and `c(n) / c(n - 1) = p(n) / q(n)` for `n > n1`
//! (polynomials `p`, `q` in `n`).  The rows before `n1` are recorded as the
//! dropped target rows (zero rows, or exceptional rows that a transfer
//! checks exactly).
//!
//! # The argument
//!
//! 0. **Left-hand factors.**  A rule whose left-hand factor is `a(n) h(x)`
//!    with `deg h >= 1` and `h` dividing every right-hand coefficient
//!    exactly is replaced by the rule with left-hand factor `a(n)` (the same
//!    rows, since `h` cancels in `Q[x]`); forward images often carry such a
//!    factor and [`certify_recurrence_equality`] needs a leader free of `x`
//!    up to a power of `x`.  This applies to the source and to the target.
//! 1. **The normalized rule.**  From the source rule
//!    `f(n) S_n = sum_j c_j(n, x) D^(d_j) S_(n - r_j)` the rule of
//!    `U_n = c(n) x^(e(n)) S_(n + s)`, `e(n) = alpha n + beta`, is computed
//!    exactly: the index shift substitutes `n + s` for `n`; then
//!    `D^d (x^m V) = sum_i C(d, i) (m)_i x^(m - i) D^(d - i) V` with the
//!    falling factorial `(m)_i` of `m = -e(n - r)`, and multiplying by
//!    `x^(e(n))` leaves the index-free powers `x^(alpha r - i)` (the whole
//!    rule is multiplied by `x^M` when one is negative); finally
//!    `c(n) / c(n - r) = p(n) ... p(n - r + 1) / (q(n) ... q(n - r + 1))`
//!    and the rule is multiplied by `q(n) ... q(n - R + 1)` (`R` the
//!    largest lag).  Multiplying an identity keeps it, so `U` satisfies the
//!    normalized rule wherever the source rule holds at `n + s`.  Its seeds
//!    are `U_(n1), ..., U_(n1 + K - 1)` computed exactly from the source
//!    rows (`K` covers the largest lag and the source's seeds), and its
//!    left-hand factor `x^M f(n + s) q(n) ... q(n - R + 1)` is nonzero at
//!    every row it computes (checked by the equality certificate), so the
//!    record generates exactly `U` (over `Q[x, 1/x]`, an integral domain).
//! 2. **The scalar is nonzero.**  `c(n1) != 0` and `p(m) / q(m) != 0` (with
//!    `q(m) != 0`) for every integer `m > n1`, decided exactly by
//!    [`crate::integer_range_sign`], so `c(n) != 0` for every `n >= n1`.
//! 3. **Equality.**  [`certify_recurrence_equality`] proves that the
//!    normalized record and the target's tail (the target rule from `n1`
//!    with the target rows as seeds) define the same rows for every index.
//!
//! Then `T_n = U_n` for every `n >= n1`.  The certificate is replayed by
//! recomputation ([`NormalizedRecurrenceEqualityCertificate::verify`]); the
//! search options only choose the normalization.
//!
//! # The search
//!
//! The shifts `0, -1, 1, -2, 2, ...` (`|s|` at most
//! [`NormalizedEqualityOptions::max_index_shift`]) and the dropped rows
//! `d = 0, 1, ...` are tried in that order.  For each, the first rows are
//! compared: zero rows must coincide, the exponent `e(n)` is the difference
//! of the valuations (it must be affine), the stripped rows must be
//! proportional, and the ratios `c(n) / c(n - 1)` must fit a rational
//! function of degree at most
//! [`NormalizedEqualityOptions::max_scalar_ratio_degree`] (an exact kernel,
//! with one more equation than unknowns).  A candidate that matches every
//! compared row is certified as above; the first certified one is kept.
//!
//! # Transfer of real-rootedness
//!
//! [`transfer_real_rootedness_along_normalization`] carries a proved claim
//! "every row from `f` on is zero or real-rooted with every root in `R`"
//! along a certified normalization, in either direction.  A nonzero scalar
//! changes no root and an index shift only relabels rows; a power of `x`
//! adds or removes the root `0`, so the region becomes the convex hull of
//! `R` and `0` unless `alpha = beta = 0`.  Rows of the conclusion that the
//! identity or the premise does not cover (the dropped rows, or rows before
//! the premise's first index) are checked exactly against the region; the
//! first index of the conclusion is the first row from which all are
//! covered or checked.
//!
//! Related: [`certify_recurrence_equality`] (the inner certificate),
//! [`crate::integer_range_sign`] (the nonvanishing of the scalar).

use super::*;
use crate::integer_range_sign::{
    decide_index_rational_function_sign, IndexSignCondition, IndexSignDecision, IntegerRange,
};
use crate::linalg::{rational_kernel_with_certificate, RationalKernelLimits};
use crate::vec_poly;
use num_bigint::BigInt;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;

type Q = BigRational;

/// JSON schema tag of [`NormalizedRecurrenceEqualityCertificate`].
pub const NORMALIZED_RECURRENCE_EQUALITY_SCHEMA: &str =
    "polytool.recurrence-normalized-equality.v1";

/// JSON schema tag of [`NormalizedRealRootednessTransfer`].
pub const NORMALIZED_REAL_ROOTEDNESS_TRANSFER_SCHEMA: &str =
    "polytool.normalized-real-rootedness-transfer.v1";

/// Largest number of rows a transfer checks exactly.
const MAX_TRANSFER_CHECKED_ROWS: usize = 64;

/// Largest index shift, dropped-row count and exponent accepted when
/// certifying a recorded normalization (a guard against hostile records).
const MAX_RECORDED_NORMALIZATION_SIZE: i64 = 64;

/// Number of rows beyond the seeds on which the normalized record is
/// compared with the directly normalized source rows (a self-check of the
/// rule transformation, not part of the proof).
const SELF_CHECK_ROWS: usize = 4;

/// Options of [`certify_normalized_recurrence_equality`] (the search only;
/// a certificate does not depend on them).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NormalizedEqualityOptions {
    /// Largest `|s|` of the index shift.  Default 4.
    pub max_index_shift: usize,
    /// Largest number of dropped leading target rows.  Default 3.
    pub max_dropped_rows: usize,
    /// Largest degree of `p` and `q` in the scalar ratio
    /// `c(n) / c(n - 1) = p(n) / q(n)`.  Default 2.
    pub max_scalar_ratio_degree: usize,
    /// Number of target rows compared in the search.  Default 16.
    pub compared_rows: usize,
    /// Least number of aligned rows a candidate must match.  Default 6.
    pub minimum_aligned_rows: usize,
    /// Largest number of candidates sent to certification.  Default 4.
    pub max_certification_attempts: usize,
}

impl Default for NormalizedEqualityOptions {
    fn default() -> Self {
        Self {
            max_index_shift: 4,
            max_dropped_rows: 3,
            max_scalar_ratio_degree: 2,
            compared_rows: 16,
            minimum_aligned_rows: 6,
            max_certification_attempts: 4,
        }
    }
}

/// A normalization `T_n = c(n) x^(alpha n + beta) S_(n + s)` for
/// `n >= n1 = (target first index) + dropped_leading_rows` (module
/// documentation).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceNormalization {
    /// The index shift `s`: the target row `n` uses the source row `n + s`.
    pub index_shift: i64,
    /// The number `d` of leading target rows outside the identity.
    pub dropped_leading_rows: usize,
    /// The slope `alpha` of the exponent of `x`.
    pub x_power_slope: i64,
    /// The offset `beta` of the exponent of `x`.
    pub x_power_offset: i64,
    /// `c(n1)`, a nonzero rational (exact string).
    pub scalar_first_value: String,
    /// `p`, ascending coefficients in `n` (exact strings), of the ratio
    /// `c(n) / c(n - 1) = p(n) / q(n)`.
    pub scalar_ratio_numerator: Vec<String>,
    /// `q`, ascending coefficients in `n` (exact strings).
    pub scalar_ratio_denominator: Vec<String>,
}

impl RecurrenceNormalization {
    /// The identity normalization (plain equality).
    pub fn identity() -> Self {
        Self {
            index_shift: 0,
            dropped_leading_rows: 0,
            x_power_slope: 0,
            x_power_offset: 0,
            scalar_first_value: "1".into(),
            scalar_ratio_numerator: vec!["1".into()],
            scalar_ratio_denominator: vec!["1".into()],
        }
    }

    /// Whether the power of `x` is `x^0` at every index.
    pub fn has_trivial_x_power(&self) -> bool {
        self.x_power_slope == 0 && self.x_power_offset == 0
    }

    /// The normalization in words.
    pub fn describe(&self) -> String {
        let scalar = |v: &[String]| {
            parse_coefficients(v)
                .map(|p| format_index_polynomial_text(&p))
                .unwrap_or_else(|_| "?".into())
        };
        format!(
            "T_n = c(n) x^({} n + {}) S_(n + {}), dropping {} leading row(s); c(n1) = {}, \
             c(n) / c(n - 1) = ({}) / ({})",
            self.x_power_slope,
            self.x_power_offset,
            self.index_shift,
            self.dropped_leading_rows,
            self.scalar_first_value,
            scalar(&self.scalar_ratio_numerator),
            scalar(&self.scalar_ratio_denominator)
        )
    }
}

/// A certificate that two records define the same rows up to a
/// [`RecurrenceNormalization`] for every index from `n1` on (module
/// documentation), schema `polytool.recurrence-normalized-equality.v1`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NormalizedRecurrenceEqualityCertificate {
    /// [`NORMALIZED_RECURRENCE_EQUALITY_SCHEMA`].
    pub schema: String,
    /// The source record `S`.
    pub source: RecurrenceJson,
    /// The target record `T`.
    pub target: RecurrenceJson,
    /// `certified` or `not_certified` (no proof found; not a disproof).
    pub status: String,
    /// The normalization (omitted when the search found no candidate).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub normalization: Option<RecurrenceNormalization>,
    /// `n1`, the first index of the identity.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub identity_first_index: Option<usize>,
    /// The target rows before `n1` (exact strings, ascending coefficients).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub dropped_target_rows: Vec<Vec<String>>,
    /// The normalized source record (rows `c(n) x^(e(n)) S_(n + s)` from
    /// `n1`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub normalized_source: Option<RecurrenceJson>,
    /// The target's tail: the target rule from `n1` with target rows as
    /// seeds.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_tail: Option<RecurrenceJson>,
    /// The exact decision of `p(m) / q(m) != 0` for every `m > n1`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub scalar_nonzero: Option<IndexSignDecision>,
    /// The equality of the target's tail (defining) and the normalized
    /// source (alternative).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub equality: Option<RecurrenceEqualityCertificate>,
    /// What the certificate establishes.
    pub statement: String,
    /// Why candidates failed.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub reasons: Vec<String>,
}

impl NormalizedRecurrenceEqualityCertificate {
    /// Whether the normalized equality is proved.
    pub fn is_certified(&self) -> bool {
        self.status == "certified"
    }

    /// `n1` and the normalization of a certified certificate.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::VerificationFailed`] when not certified.
    pub fn certified_normalization(
        &self,
    ) -> Result<(usize, &RecurrenceNormalization), RecurrenceError> {
        match (
            self.is_certified(),
            self.identity_first_index,
            &self.normalization,
        ) {
            (true, Some(first), Some(normalization)) => Ok((first, normalization)),
            _ => Err(RecurrenceError::VerificationFailed(
                "the normalized equality is not certified".into(),
            )),
        }
    }

    /// Replay: recompute the certification of the recorded normalization
    /// (or, without one, check that nothing is claimed) and compare.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::VerificationFailed`] when the certificate differs
    /// from its recomputation, has another schema or claims without a
    /// normalization; errors of the records.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        if self.schema != NORMALIZED_RECURRENCE_EQUALITY_SCHEMA {
            return Err(RecurrenceError::VerificationFailed(format!(
                "unsupported schema `{}`",
                self.schema
            )));
        }
        let Some(normalization) = &self.normalization else {
            if self.is_certified()
                || self.identity_first_index.is_some()
                || self.normalized_source.is_some()
                || self.equality.is_some()
            {
                return Err(RecurrenceError::VerificationFailed(
                    "a normalized equality without a normalization claims nothing".into(),
                ));
            }
            return Ok(());
        };
        let mut recomputed =
            certify_recurrence_normalization(&self.source, &self.target, normalization)?;
        // The reasons of earlier search candidates are not recomputed.
        let earlier = self.reasons.len().saturating_sub(recomputed.reasons.len());
        let mut reasons = self.reasons[..earlier].to_vec();
        reasons.append(&mut recomputed.reasons);
        recomputed.reasons = reasons;
        if &recomputed != self {
            return Err(RecurrenceError::VerificationFailed(
                "the normalized equality certificate differs from its recomputation".into(),
            ));
        }
        Ok(())
    }
}

/// `p` without trailing zeros (the zero polynomial is empty).
fn strip(p: &[Q]) -> Vec<Q> {
    let mut out = p.to_vec();
    vec_poly::strip_trailing_zeros(&mut out);
    out
}

fn q_int(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

fn parse_coefficients(values: &[String]) -> Result<Vec<Q>, RecurrenceError> {
    values
        .iter()
        .map(|v| parse_rational_coeff(v).map_err(RecurrenceError::from))
        .collect()
}

fn format_coefficients(values: &[Q]) -> Vec<String> {
    let trimmed = strip(values);
    if trimmed.is_empty() {
        return vec!["0".into()];
    }
    trimmed.iter().map(format_rational_coeff).collect()
}

/// `p(n)` as text.
fn format_index_polynomial_text(p: &[Q]) -> String {
    let p = strip(p);
    if p.is_empty() {
        return "0".into();
    }
    let mut parts = Vec::new();
    for (i, c) in p.iter().enumerate().rev() {
        if c.is_zero() {
            continue;
        }
        let c = format_rational_coeff(c);
        parts.push(match i {
            0 => c,
            1 => format!("{c} n"),
            _ => format!("{c} n^{i}"),
        });
    }
    parts.join(" + ")
}

/// `p(value)` for an integer value.
fn evaluate_index_polynomial(p: &[Q], value: i64) -> Q {
    vec_poly::evaluate(p, &q_int(value))
}

/// `p(n + delta)`.
fn shift_index_polynomial(p: &[Q], delta: i64) -> Vec<Q> {
    let column =
        BivarPoly::from_coefficients_normalized(p.iter().map(|c| vec![c.clone()]).collect());
    let shifted = column.shift_index(delta);
    strip(
        &shifted
            .coeffs
            .iter()
            .map(|row| row.first().cloned().unwrap_or_else(Q::zero))
            .collect::<Vec<_>>(),
    )
}

/// `p(n) b(n, x)`.
fn bivar_times_index_polynomial(b: &BivarPoly, p: &[Q]) -> BivarPoly {
    let width = b.coeffs.iter().map(Vec::len).max().unwrap_or(0).max(1);
    let mut out = vec![vec![Q::zero(); width]; b.coeffs.len() + p.len()];
    for (i, row) in b.coeffs.iter().enumerate() {
        for (j, c) in row.iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            for (k, a) in p.iter().enumerate() {
                out[i + k][j] += c * a;
            }
        }
    }
    BivarPoly::from_coefficients_normalized(out)
}

/// `x^k b(n, x)`.
fn bivar_times_x_power(b: &BivarPoly, k: usize) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(
        b.coeffs
            .iter()
            .map(|row| {
                let mut shifted = vec![Q::zero(); k];
                shifted.extend(row.iter().cloned());
                shifted
            })
            .collect(),
    )
}

fn bivar_add(a: &BivarPoly, b: &BivarPoly) -> BivarPoly {
    let rows = a.coeffs.len().max(b.coeffs.len());
    let width = a
        .coeffs
        .iter()
        .chain(&b.coeffs)
        .map(Vec::len)
        .max()
        .unwrap_or(0)
        .max(1);
    let mut out = vec![vec![Q::zero(); width]; rows.max(1)];
    for source in [a, b] {
        for (i, row) in source.coeffs.iter().enumerate() {
            for (j, c) in row.iter().enumerate() {
                out[i][j] += c;
            }
        }
    }
    BivarPoly::from_coefficients_normalized(out)
}

fn binomial(n: usize, k: usize) -> Q {
    let mut value = Q::one();
    for i in 0..k {
        value = value * q_int((n - i) as i64) / q_int((i + 1) as i64);
    }
    value
}

/// The falling factorial `(m)_i = m (m - 1) ... (m - i + 1)` of the affine
/// polynomial `m(n) = a n + b`, as a polynomial in `n`.
fn falling_factorial(a: i64, b: i64, i: usize) -> Vec<Q> {
    let mut out = vec![Q::one()];
    for t in 0..i {
        out = vec_poly::mul(&out, &[q_int(b - t as i64), q_int(a)]);
    }
    out
}

/// `p / h` when `h` divides `p` exactly in `Q[x]`.
fn exact_polynomial_division(p: &[Q], h: &[Q]) -> Option<Vec<Q>> {
    let h = strip(h);
    let mut rest = strip(p);
    let lead = h.last()?.clone();
    if rest.is_empty() {
        return Some(Vec::new());
    }
    if rest.len() < h.len() {
        return None;
    }
    let mut quotient = vec![Q::zero(); rest.len() - h.len() + 1];
    while rest.len() >= h.len() && !rest.is_empty() {
        let shift = rest.len() - h.len();
        let factor = rest.last().expect("nonempty") / &lead;
        for (k, c) in h.iter().enumerate() {
            rest[shift + k] -= c * &factor;
        }
        quotient[shift] = factor;
        rest = strip(&rest);
    }
    rest.is_empty().then_some(quotient)
}

/// The rule with a left-hand factor `a(n) h(x)` (`deg h >= 1`) replaced by
/// `a(n)` when `h` divides every right-hand coefficient exactly (the same
/// rows: `h` cancels in the integral domain `Q[x]`); otherwise the rule.
/// [`certify_recurrence_equality`] needs a leader free of `x` (up to a
/// power of `x`), and forward images often carry such a factor.
fn rule_with_separable_left_hand_factor_cancelled(rule: &Recurrence) -> Recurrence {
    let Some(f) = &rule.denominator else {
        return rule.clone();
    };
    let rows: Vec<Vec<Q>> = f.coeffs.iter().map(|r| strip(r)).collect();
    let Some(h) = rows.iter().find(|r| !r.is_empty()) else {
        return rule.clone();
    };
    if h.len() < 2 {
        return rule.clone();
    }
    let lead = h.last().expect("nonempty").clone();
    let h: Vec<Q> = h.iter().map(|c| c / &lead).collect();
    let mut a = Vec::with_capacity(rows.len());
    for row in &rows {
        match exact_polynomial_division(row, &h) {
            Some(quotient) if quotient.len() <= 1 => {
                a.push(quotient.first().cloned().unwrap_or_else(Q::zero))
            }
            _ => return rule.clone(),
        }
    }
    let mut terms = Vec::with_capacity(rule.terms.len());
    for term in &rule.terms {
        let mut coeffs = Vec::with_capacity(term.coeff.coeffs.len());
        for row in &term.coeff.coeffs {
            match exact_polynomial_division(row, &h) {
                Some(quotient) => coeffs.push(quotient),
                None => return rule.clone(),
            }
        }
        let mut term = term.clone();
        term.coeff = BivarPoly::from_coefficients_normalized(coeffs);
        terms.push(term);
    }
    let a = strip(&a);
    let denominator =
        BivarPoly::from_coefficients_normalized(a.iter().map(|c| vec![c.clone()]).collect());
    Recurrence {
        terms,
        denominator: (!denominator.is_one()).then_some(denominator),
        inhomogeneous: rule.inhomogeneous.clone(),
    }
}

/// Add the term `(offset, derivative, sign, coefficient)`, merging equal
/// keys.
fn push_term(terms: &mut Vec<RecurrenceTerm>, term: RecurrenceTerm) {
    if term.coeff.is_zero() {
        return;
    }
    if let Some(existing) = terms.iter_mut().find(|t| {
        t.offset == term.offset && t.deriv_order == term.deriv_order && t.sign == term.sign
    }) {
        existing.coeff = bivar_add(&existing.coeff, &term.coeff);
    } else {
        terms.push(term);
    }
}

/// The parsed scalar of a normalization: `(c(n1), p, q)`.
struct Scalar {
    first_value: Q,
    numerator: Vec<Q>,
    denominator: Vec<Q>,
}

impl Scalar {
    fn parse(normalization: &RecurrenceNormalization) -> Result<Self, RecurrenceError> {
        Ok(Self {
            first_value: parse_rational_coeff(&normalization.scalar_first_value)?,
            numerator: strip(&parse_coefficients(&normalization.scalar_ratio_numerator)?),
            denominator: strip(&parse_coefficients(
                &normalization.scalar_ratio_denominator,
            )?),
        })
    }

    fn is_constant_one_ratio(&self) -> bool {
        self.numerator.len() == 1
            && self.denominator.len() == 1
            && self.numerator == self.denominator
    }

    /// `c(n1), ..., c(n1 + count - 1)`; `None` when `q` vanishes.
    fn values(&self, n1: usize, count: usize) -> Option<Vec<Q>> {
        let mut out = Vec::with_capacity(count);
        let mut value = self.first_value.clone();
        for k in 0..count {
            if k > 0 {
                let n = (n1 + k) as i64;
                let q = evaluate_index_polynomial(&self.denominator, n);
                if q.is_zero() {
                    return None;
                }
                value = value * evaluate_index_polynomial(&self.numerator, n) / q;
            }
            out.push(value.clone());
        }
        Some(out)
    }
}

/// The rule of `U_n = c(n) x^(alpha n + beta) S_(n + s)` (module
/// documentation, step 1).
fn normalized_rule(
    rule: &Recurrence,
    normalization: &RecurrenceNormalization,
    scalar: &Scalar,
) -> Result<Recurrence, String> {
    if rule.inhomogeneous.as_ref().is_some_and(|g| !g.is_zero()) {
        return Err("the source rule has an inhomogeneous term".into());
    }
    let shifted = rule.shift_index(normalization.index_shift);
    let alpha = normalization.x_power_slope;
    let beta = normalization.x_power_offset;
    // x-power: (offset, derivative, sign, exponent, coefficient).
    let mut powered: Vec<(RecurrenceTerm, i64)> = Vec::new();
    for term in &shifted.terms {
        let r = term.offset as i64;
        for i in 0..=term.deriv_order {
            // m(n) = -e(n - r) = -alpha n + (alpha r - beta).
            let factor = vec_poly::mul(
                &[binomial(term.deriv_order, i)],
                &falling_factorial(-alpha, alpha * r - beta, i),
            );
            if vec_poly::is_zero(&factor) {
                continue;
            }
            powered.push((
                RecurrenceTerm {
                    offset: term.offset,
                    deriv_order: term.deriv_order - i,
                    sign: term.sign,
                    coeff: bivar_times_index_polynomial(&term.coeff, &factor),
                },
                alpha * r - i as i64,
            ));
        }
    }
    let lift = powered
        .iter()
        .map(|(_, exponent)| (-exponent).max(0))
        .max()
        .unwrap_or(0) as usize;
    let one = BivarPoly::from_coefficients_normalized(vec![vec![Q::one()]]);
    let mut denominator = bivar_times_x_power(shifted.denominator.as_ref().unwrap_or(&one), lift);
    // Scalar: multiply by q(n) ... q(n - R + 1).
    let largest_lag = powered.iter().map(|(t, _)| t.offset).max().unwrap_or(0);
    let scalar_trivial = scalar.is_constant_one_ratio();
    let product = |from: usize, to: usize, p: &[Q]| -> Vec<Q> {
        let mut out = vec![Q::one()];
        for k in from..to {
            out = vec_poly::mul(&out, &shift_index_polynomial(p, -(k as i64)));
        }
        out
    };
    if !scalar_trivial {
        denominator = bivar_times_index_polynomial(
            &denominator,
            &product(0, largest_lag, &scalar.denominator),
        );
    }
    let mut terms = Vec::new();
    for (mut term, exponent) in powered {
        term.coeff = bivar_times_x_power(&term.coeff, (exponent + lift as i64) as usize);
        if !scalar_trivial {
            let factor = vec_poly::mul(
                &product(0, term.offset, &scalar.numerator),
                &product(term.offset, largest_lag, &scalar.denominator),
            );
            term.coeff = bivar_times_index_polynomial(&term.coeff, &factor);
        }
        push_term(&mut terms, term);
    }
    if terms.is_empty() {
        return Err("the normalized rule has no term".into());
    }
    Ok(Recurrence {
        terms,
        denominator: (!denominator.is_one()).then_some(denominator),
        inhomogeneous: None,
    })
}

/// The exact rows of a record from its first index, `count` of them.
fn record_rows(record: &RecurrenceJson, count: usize) -> Result<Vec<Vec<Q>>, String> {
    let (rule, first_index, seeds) = record
        .to_recurrence_parts()
        .map_err(|error| error.to_string())?;
    if rule.inhomogeneous.as_ref().is_some_and(|g| !g.is_zero()) {
        return Err("the rule has an inhomogeneous term".into());
    }
    rule.generate_rows_rational(&seeds, first_index, count.max(seeds.len()))
        .map(|rows| {
            rows.into_iter()
                .take(count)
                .map(|row| strip(&row))
                .collect()
        })
        .map_err(|error| format!("the rows cannot be generated: {error}"))
}

/// The lowest power of `x` of a nonzero row.
fn valuation(row: &[Q]) -> usize {
    row.iter().position(|c| !c.is_zero()).unwrap_or(0)
}

/// `c x^e row` (exact; `None` when a negative power does not divide).
fn scale_row(row: &[Q], c: &Q, e: i64) -> Option<Vec<Q>> {
    let row = strip(row);
    if row.is_empty() {
        return Some(row);
    }
    let scaled: Vec<Q> = row.iter().map(|a| a * c).collect();
    if e >= 0 {
        let mut out = vec![Q::zero(); e as usize];
        out.extend(scaled);
        Some(strip(&out))
    } else {
        let drop = (-e) as usize;
        if valuation(&row) < drop {
            return None;
        }
        Some(strip(&scaled[drop..]))
    }
}

fn row_strings(row: &[Q]) -> Vec<String> {
    strip(row).iter().map(format_rational_coeff).collect()
}

/// Certify one recorded normalization (module documentation): the result
/// is `certified` or `not_certified` with the reasons.  Used by the search
/// and by [`NormalizedRecurrenceEqualityCertificate::verify`].
///
/// # Errors
///
/// Malformed records or normalizations ([`RecurrenceError::InvalidInput`]).
pub fn certify_recurrence_normalization(
    source: &RecurrenceJson,
    target: &RecurrenceJson,
    normalization: &RecurrenceNormalization,
) -> Result<NormalizedRecurrenceEqualityCertificate, RecurrenceError> {
    let (source_rule, source_first, source_seeds) = source.to_recurrence_parts()?;
    let (target_rule, target_first, target_seeds) = target.to_recurrence_parts()?;
    let source_rule = rule_with_separable_left_hand_factor_cancelled(&source_rule);
    let target_rule = rule_with_separable_left_hand_factor_cancelled(&target_rule);
    let scalar = Scalar::parse(normalization)?;
    for (value, name) in [
        (normalization.index_shift, "index shift"),
        (
            normalization.dropped_leading_rows as i64,
            "dropped row count",
        ),
        (normalization.x_power_slope, "exponent slope"),
        (normalization.x_power_offset, "exponent offset"),
    ] {
        if value.abs() > MAX_RECORDED_NORMALIZATION_SIZE {
            return Err(RecurrenceError::InvalidInput(format!(
                "the {name} {value} exceeds {MAX_RECORDED_NORMALIZATION_SIZE} in absolute value"
            )));
        }
    }
    let n1 = target_first + normalization.dropped_leading_rows;
    let mut certificate = NormalizedRecurrenceEqualityCertificate {
        schema: NORMALIZED_RECURRENCE_EQUALITY_SCHEMA.into(),
        source: source.clone(),
        target: target.clone(),
        status: "not_certified".into(),
        normalization: Some(normalization.clone()),
        identity_first_index: Some(n1),
        dropped_target_rows: Vec::new(),
        normalized_source: None,
        target_tail: None,
        scalar_nonzero: None,
        equality: None,
        statement: String::new(),
        reasons: Vec::new(),
    };
    let fail = |mut certificate: NormalizedRecurrenceEqualityCertificate, reason: String| {
        certificate.statement = format!(
            "not certified for the normalization {} (not a disproof)",
            normalization.describe()
        );
        certificate.reasons.push(reason);
        Ok(certificate)
    };
    let source_index = n1 as i64 + normalization.index_shift;
    if source_index < source_first as i64 {
        return fail(
            certificate,
            format!(
                "the identity would use the source row {source_index}, before its first index \
                 {source_first}"
            ),
        );
    }
    // Step 2: the scalar is nonzero for every n >= n1.
    if scalar.first_value.is_zero() || scalar.numerator.is_empty() || scalar.denominator.is_empty()
    {
        return fail(certificate, "the scalar c(n1) or its ratio is zero".into());
    }
    let decision = decide_index_rational_function_sign(
        &scalar.numerator,
        &scalar.denominator,
        &IntegerRange::starting_at(BigInt::from(n1 + 1)),
        IndexSignCondition::Nonzero,
    );
    let holds = decision.holds();
    let description = decision.describe("p(n)/q(n)");
    certificate.scalar_nonzero = Some(decision);
    if !holds {
        return fail(
            certificate,
            format!("the scalar ratio is not proved nonzero for n > {n1}: {description}"),
        );
    }
    // Step 1: the normalized rule and its seeds.
    let rule = match normalized_rule(&source_rule, normalization, &scalar) {
        Ok(rule) => rule,
        Err(reason) => return fail(certificate, reason),
    };
    let lag = rule.max_offset().max(1);
    let source_rule_start = (source_first + source_seeds.len()) as i64;
    let seed_count = lag.max((source_rule_start - source_index).max(0) as usize);
    let checked = seed_count + SELF_CHECK_ROWS;
    let last_source = (source_index - source_first as i64) as usize + checked;
    let source_rows = match record_rows(source, last_source) {
        Ok(rows) => rows,
        Err(reason) => return fail(certificate, format!("the source: {reason}")),
    };
    let Some(scalars) = scalar.values(n1, checked) else {
        return fail(certificate, "the scalar ratio has a pole".into());
    };
    let mut direct = Vec::with_capacity(checked);
    for (k, c) in scalars.iter().enumerate() {
        let n = (n1 + k) as i64;
        let row = &source_rows[(n + normalization.index_shift) as usize - source_first];
        let exponent = normalization.x_power_slope * n + normalization.x_power_offset;
        match scale_row(row, c, exponent) {
            Some(row) => direct.push(row),
            None => {
                return fail(
                    certificate,
                    format!(
                        "x^{} does not divide the source row {} (exact)",
                        -exponent,
                        n + normalization.index_shift
                    ),
                )
            }
        }
    }
    let normalized =
        RecurrenceJson::from_recurrence_rational(&rule, n1, &direct[..seed_count], None);
    match rule.generate_rows_rational(&direct[..seed_count], n1, checked) {
        Ok(rows) => {
            let rows: Vec<Vec<Q>> = rows.iter().map(|r| strip(r)).collect();
            if rows != direct {
                return fail(
                    certificate,
                    "the normalized rule does not reproduce the normalized source rows (internal \
                     check)"
                        .into(),
                );
            }
        }
        Err(error) => {
            return fail(
                certificate,
                format!("the normalized rule cannot be replayed: {error}"),
            )
        }
    }
    certificate.normalized_source = Some(normalized.clone());
    // The target's tail and its dropped rows.
    let tail_count = normalization.dropped_leading_rows + target_seeds.len().max(1);
    let target_rows = match record_rows(target, tail_count) {
        Ok(rows) => rows,
        Err(reason) => return fail(certificate, format!("the target: {reason}")),
    };
    certificate.dropped_target_rows = target_rows[..normalization.dropped_leading_rows]
        .iter()
        .map(|r| row_strings(r))
        .collect();
    let tail_seeds = &target_rows[normalization.dropped_leading_rows..];
    let tail = RecurrenceJson::from_recurrence_rational(
        &target_rule,
        n1,
        &tail_seeds[..target_seeds.len().max(1).min(tail_seeds.len())],
        None,
    );
    certificate.target_tail = Some(tail.clone());
    // Step 3: equality.
    let equality = certify_recurrence_equality(&tail, &normalized)?;
    let certified = equality.is_certified();
    if !certified {
        let reasons = equality.reasons.join("; ");
        certificate.equality = Some(equality);
        return fail(
            certificate,
            format!("the equality of the target's tail and the normalized source: {reasons}"),
        );
    }
    certificate.equality = Some(equality);
    certificate.status = "certified".into();
    certificate.statement = format!(
        "for every n >= {n1}: T_n = c(n) x^({alpha} n + {beta}) S_(n + {s}) with c({n1}) = {c} \
         and c(n) / c(n - 1) = ({p}) / ({q}) nonzero for every n > {n1} (exact); the \
         normalized source rule is computed exactly from the source rule and certified equal \
         to the target rule from n = {n1} on (all-index identity, no finite fit){dropped}",
        alpha = normalization.x_power_slope,
        beta = normalization.x_power_offset,
        s = normalization.index_shift,
        c = normalization.scalar_first_value,
        p = format_index_polynomial_text(&scalar.numerator),
        q = format_index_polynomial_text(&scalar.denominator),
        dropped = if normalization.dropped_leading_rows == 0 {
            String::new()
        } else {
            format!(
                "; the target rows {target_first}, ..., {} are outside the identity and \
                 recorded exactly",
                n1 - 1
            )
        }
    );
    Ok(certificate)
}

/// The search order of the index shifts: `0, -1, 1, -2, 2, ...`.
fn shift_order(bound: usize) -> Vec<i64> {
    let mut out = vec![0];
    for k in 1..=bound as i64 {
        out.push(-k);
        out.push(k);
    }
    out
}

/// A rational function `p / q` with `deg p, deg q <= degree` through the
/// points `(n, value)` (exact kernel; one more equation than needed).
fn fit_ratio(points: &[(i64, Q)], max_degree: usize) -> Option<(Vec<Q>, Vec<Q>)> {
    if points.iter().all(|(_, v)| v.is_one()) {
        return Some((vec![Q::one()], vec![Q::one()]));
    }
    for degree in 0..=max_degree {
        let unknowns = 2 * (degree + 1);
        if points.len() < unknowns {
            break;
        }
        let matrix: Vec<Vec<Q>> = points
            .iter()
            .map(|(n, value)| {
                let mut row = Vec::with_capacity(unknowns);
                let mut power = Q::one();
                let mut powers = Vec::with_capacity(degree + 1);
                for _ in 0..=degree {
                    powers.push(power.clone());
                    power *= q_int(*n);
                }
                row.extend(powers.iter().cloned());
                row.extend(powers.iter().map(|p| -(p * value)));
                row
            })
            .collect();
        let limits = RationalKernelLimits {
            max_entries: 4_096,
            max_operations: 1_000_000,
            max_coefficient_bits: 4_096,
        };
        let Ok(kernel) = rational_kernel_with_certificate(&matrix, unknowns, limits) else {
            continue;
        };
        let Some(vector) = kernel.basis.first() else {
            continue;
        };
        let p = strip(&vector[..=degree]);
        let mut q = strip(&vector[degree + 1..]);
        if q.is_empty() || p.is_empty() {
            continue;
        }
        if points
            .iter()
            .any(|(n, _)| evaluate_index_polynomial(&q, *n).is_zero())
        {
            continue;
        }
        let lead = q.last().cloned().expect("nonzero");
        q = q.iter().map(|c| c / &lead).collect();
        let p: Vec<Q> = p.iter().map(|c| c / &lead).collect();
        return Some((p, q));
    }
    None
}

/// The normalization matching the compared rows for a shift and a number
/// of dropped rows, if any (module documentation, "The search").
fn match_rows(
    source_rows: &[Vec<Q>],
    source_first: usize,
    target_rows: &[Vec<Q>],
    target_first: usize,
    shift: i64,
    dropped: usize,
    options: &NormalizedEqualityOptions,
) -> Option<RecurrenceNormalization> {
    let n1 = target_first + dropped;
    // (n, Some((exponent, scalar))) for nonzero rows, None for zero rows.
    let mut aligned: Vec<(i64, Option<(i64, Q)>)> = Vec::new();
    for (offset, target) in target_rows.iter().enumerate().skip(dropped) {
        let n = (target_first + offset) as i64;
        let m = n + shift;
        if m < source_first as i64 {
            return None;
        }
        let Some(source) = source_rows.get(m as usize - source_first) else {
            break;
        };
        match (target.is_empty(), source.is_empty()) {
            (true, true) => aligned.push((n, None)),
            (false, false) => {
                let (vt, vs) = (valuation(target), valuation(source));
                let (t, s) = (&target[vt..], &source[vs..]);
                if t.len() != s.len() {
                    return None;
                }
                let c = &t[0] / &s[0];
                if t.iter().zip(s).any(|(a, b)| a != &(b * &c)) {
                    return None;
                }
                aligned.push((n, Some((vt as i64 - vs as i64, c))));
            }
            _ => return None,
        }
    }
    if aligned.len() < options.minimum_aligned_rows {
        return None;
    }
    let nonzero: Vec<(i64, i64, Q)> = aligned
        .iter()
        .filter_map(|(n, v)| v.as_ref().map(|(e, c)| (*n, *e, c.clone())))
        .collect();
    if nonzero.len() < 3 {
        return None;
    }
    // Affine exponent.
    let (n_a, e_a, _) = &nonzero[0];
    let (n_b, e_b, _) = &nonzero[1];
    let (rise, run) = (e_b - e_a, n_b - n_a);
    if run == 0 || rise % run != 0 {
        return None;
    }
    let slope = rise / run;
    let offset = e_a - slope * n_a;
    if nonzero.iter().any(|(n, e, _)| slope * n + offset != *e) {
        return None;
    }
    // Hypergeometric scalar from consecutive nonzero rows.
    let ratios: Vec<(i64, Q)> = nonzero
        .windows(2)
        .filter(|w| w[1].0 == w[0].0 + 1)
        .map(|w| (w[1].0, &w[1].2 / &w[0].2))
        .collect();
    let (p, q) = fit_ratio(&ratios, options.max_scalar_ratio_degree)?;
    // c(n1) from the first nonzero row.
    let (n_first, _, c_first) = &nonzero[0];
    let mut first_value = c_first.clone();
    for m in (n1 as i64 + 1..=*n_first).rev() {
        let numerator = evaluate_index_polynomial(&p, m);
        if numerator.is_zero() {
            return None;
        }
        first_value = first_value * evaluate_index_polynomial(&q, m) / numerator;
    }
    let normalization = RecurrenceNormalization {
        index_shift: shift,
        dropped_leading_rows: dropped,
        x_power_slope: slope,
        x_power_offset: offset,
        scalar_first_value: format_rational_coeff(&first_value),
        scalar_ratio_numerator: format_coefficients(&p),
        scalar_ratio_denominator: format_coefficients(&q),
    };
    // Every compared row must follow the scalar.
    let scalar = Scalar::parse(&normalization).ok()?;
    let values = scalar.values(n1, aligned.len())?;
    for ((n, v), c) in aligned.iter().zip(&values) {
        debug_assert!(*n >= n1 as i64);
        if let Some((_, observed)) = v {
            if observed != c {
                return None;
            }
        }
    }
    Some(normalization)
}

/// Search the normalization menu and certify the first candidate that
/// matches the compared rows (module documentation).  The result is
/// `certified` or `not_certified` with the reasons.
///
/// # Errors
///
/// Malformed records ([`RecurrenceJson::to_recurrence_parts`]).
pub fn certify_normalized_recurrence_equality(
    source: &RecurrenceJson,
    target: &RecurrenceJson,
    options: &NormalizedEqualityOptions,
) -> Result<NormalizedRecurrenceEqualityCertificate, RecurrenceError> {
    source.to_recurrence_parts()?;
    target.to_recurrence_parts()?;
    let mut certificate = NormalizedRecurrenceEqualityCertificate {
        schema: NORMALIZED_RECURRENCE_EQUALITY_SCHEMA.into(),
        source: source.clone(),
        target: target.clone(),
        status: "not_certified".into(),
        normalization: None,
        identity_first_index: None,
        dropped_target_rows: Vec::new(),
        normalized_source: None,
        target_tail: None,
        scalar_nonzero: None,
        equality: None,
        statement: String::new(),
        reasons: Vec::new(),
    };
    let target_count = options.compared_rows + options.max_dropped_rows;
    let source_count = target_count + 2 * options.max_index_shift + 1;
    let rows = record_rows(source, source_count)
        .map_err(|r| format!("the source: {r}"))
        .and_then(|s| {
            record_rows(target, target_count)
                .map(|t| (s, t))
                .map_err(|r| format!("the target: {r}"))
        });
    let (source_rows, target_rows) = match rows {
        Ok(rows) => rows,
        Err(reason) => {
            certificate.statement = "not certified: the rows cannot be generated".into();
            certificate.reasons.push(reason);
            return Ok(certificate);
        }
    };
    let mut reasons = Vec::new();
    let mut attempts = 0;
    for shift in shift_order(options.max_index_shift) {
        for dropped in 0..=options.max_dropped_rows {
            if attempts >= options.max_certification_attempts {
                break;
            }
            let Some(normalization) = match_rows(
                &source_rows,
                source.first_index,
                &target_rows,
                target.first_index,
                shift,
                dropped,
                options,
            ) else {
                continue;
            };
            attempts += 1;
            let mut attempt = certify_recurrence_normalization(source, target, &normalization)?;
            if attempt.is_certified() {
                reasons.append(&mut attempt.reasons);
                attempt.reasons = reasons;
                return Ok(attempt);
            }
            reasons.extend(
                attempt
                    .reasons
                    .iter()
                    .map(|r| format!("{}: {r}", normalization.describe())),
            );
        }
    }
    certificate.statement = if attempts == 0 {
        "not certified: no normalization of the menu matches the compared rows".into()
    } else {
        "not certified: no matching normalization was certified (not a disproof)".into()
    };
    certificate.reasons = reasons;
    Ok(certificate)
}

/// The direction of a [`NormalizedRealRootednessTransfer`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NormalizationTransferDirection {
    /// From a claim about the source to one about the target.
    SourceToTarget,
    /// From a claim about the target to one about the source.
    TargetToSource,
}

/// A real-rootedness claim carried along a certified normalized equality
/// (module documentation, "Transfer of real-rootedness"), schema
/// `polytool.normalized-real-rootedness-transfer.v1`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NormalizedRealRootednessTransfer {
    /// [`NORMALIZED_REAL_ROOTEDNESS_TRANSFER_SCHEMA`].
    pub schema: String,
    /// The certified normalized equality.
    pub equality: NormalizedRecurrenceEqualityCertificate,
    /// The direction.
    pub direction: NormalizationTransferDirection,
    /// The premise: every row of the premise sequence from this index on
    /// is zero or real-rooted with every root in `premise_region`.
    pub premise_from_index: usize,
    /// The premise's region (index-independent).
    pub premise_region: RootInterval,
    /// The conclusion's first index.
    pub conclusion_from_index: usize,
    /// The conclusion's region.
    pub conclusion_region: RootInterval,
    /// Rows of the conclusion sequence checked exactly (not covered by the
    /// identity and the premise), ascending.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub checked_rows: Vec<usize>,
    /// What the transfer establishes.
    pub statement: String,
}

impl NormalizedRealRootednessTransfer {
    /// The record whose rows the conclusion is about.
    pub fn conclusion_record(&self) -> &RecurrenceJson {
        match self.direction {
            NormalizationTransferDirection::SourceToTarget => &self.equality.target,
            NormalizationTransferDirection::TargetToSource => &self.equality.source,
        }
    }

    /// The record whose rows the premise is about.
    pub fn premise_record(&self) -> &RecurrenceJson {
        match self.direction {
            NormalizationTransferDirection::SourceToTarget => &self.equality.source,
            NormalizationTransferDirection::TargetToSource => &self.equality.target,
        }
    }

    /// Replay: verify the equality and recompute the transfer.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::VerificationFailed`] when anything differs.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        if self.schema != NORMALIZED_REAL_ROOTEDNESS_TRANSFER_SCHEMA {
            return Err(RecurrenceError::VerificationFailed(format!(
                "unsupported schema `{}`",
                self.schema
            )));
        }
        self.equality.verify()?;
        let recomputed = transfer_real_rootedness_along_normalization(
            &self.equality,
            self.direction,
            self.premise_from_index,
            &self.premise_region,
        )?;
        if &recomputed != self {
            return Err(RecurrenceError::VerificationFailed(
                "the real-rootedness transfer differs from its recomputation".into(),
            ));
        }
        Ok(())
    }
}

/// The exact interval of an index-independent region.
pub(super) fn algebraic_region(
    region: &RootInterval,
) -> Result<AlgebraicInterval, RecurrenceError> {
    if !region.is_index_independent() {
        return Err(RecurrenceError::Unsupported(format!(
            "the region {region} depends on the index"
        )));
    }
    region.check_exact_endpoints()?;
    let end = |end: &Option<RootEndpoint>| -> Result<Option<AlgebraicEndpoint>, RecurrenceError> {
        end.as_ref()
            .map(|e| {
                let value = match &e.exact {
                    Some(exact) => exact.to_quadratic_number()?,
                    None => QuadraticNumber::rational(parse_rational_coeff(&e.value)?),
                };
                Ok(AlgebraicEndpoint {
                    value,
                    closed: e.closed,
                })
            })
            .transpose()
    };
    Ok(AlgebraicInterval {
        lower: end(&region.lower)?,
        upper: end(&region.upper)?,
    })
}

/// The smallest interval containing `region` and `0` (an endpoint that
/// moves becomes the closed rational `0`).
fn hull_with_zero(region: &RootInterval) -> Result<RootInterval, RecurrenceError> {
    let exact = algebraic_region(region)?;
    let zero = QuadraticNumber::rational(Q::zero());
    if exact.contains(&zero) {
        return Ok(region.clone());
    }
    let mut out = region.clone();
    let closed_zero = Some(RootEndpoint::rational(&Q::zero(), true));
    if exact
        .upper
        .as_ref()
        .is_some_and(|u| u.value.cmp_rational(&Q::zero()) != Ordering::Greater)
    {
        out.upper = closed_zero;
    } else {
        out.lower = closed_zero;
    }
    Ok(out)
}

/// Carry a real-rootedness claim along a certified normalized equality
/// (module documentation, "Transfer of real-rootedness").  The premise
/// "every premise row `P_m`, `m >= premise_from_index`, is zero or
/// real-rooted with every root in `premise_region`" gives the conclusion
/// recorded in the result.
///
/// # Errors
///
/// [`RecurrenceError::VerificationFailed`] when the equality is not
/// certified, [`RecurrenceError::Unsupported`] for an index-dependent
/// region, and errors of the records.
pub fn transfer_real_rootedness_along_normalization(
    equality: &NormalizedRecurrenceEqualityCertificate,
    direction: NormalizationTransferDirection,
    premise_from_index: usize,
    premise_region: &RootInterval,
) -> Result<NormalizedRealRootednessTransfer, RecurrenceError> {
    let (n1, normalization) = equality.certified_normalization()?;
    let region = if normalization.has_trivial_x_power() {
        algebraic_region(premise_region)?;
        premise_region.clone()
    } else {
        hull_with_zero(premise_region)?
    };
    let exact = algebraic_region(&region)?;
    let s = normalization.index_shift;
    // The first conclusion row covered by the identity and the premise, in
    // the conclusion's own indexing, and the conclusion record.
    let (covered_from, record) = match direction {
        NormalizationTransferDirection::SourceToTarget => (
            (n1 as i64).max(premise_from_index as i64 - s),
            &equality.target,
        ),
        NormalizationTransferDirection::TargetToSource => {
            (n1.max(premise_from_index) as i64 + s, &equality.source)
        }
    };
    let first = record.first_index as i64;
    let covered_from = covered_from.max(first);
    let uncovered = (covered_from - first) as usize;
    let mut checked_rows = Vec::new();
    let mut conclusion_from = covered_from as usize;
    if uncovered > 0 {
        let start = uncovered.saturating_sub(MAX_TRANSFER_CHECKED_ROWS);
        let rows = record_rows(record, uncovered).map_err(RecurrenceError::InvalidInput)?;
        // Extend downwards while the rows pass.
        for offset in (start..uncovered).rev() {
            let row = &rows[offset];
            if !row.is_empty() && !exact.contains_all_roots_of(row) {
                break;
            }
            conclusion_from = first as usize + offset;
            checked_rows.push(conclusion_from);
        }
        checked_rows.reverse();
    }
    let (premise_name, conclusion_name) = match direction {
        NormalizationTransferDirection::SourceToTarget => ("S", "T"),
        NormalizationTransferDirection::TargetToSource => ("T", "S"),
    };
    let statement = format!(
        "{premise_name}_m is zero or real-rooted with every root in {premise} for m >= \
         {premise_from_index}; by the certified normalization (nonzero scalar, {power}) every \
         row {conclusion_name}_n, n >= {conclusion_from}, is zero or real-rooted with every root \
         in {conclusion}{checked}",
        premise = format_root_interval_exact(premise_region),
        power = if normalization.has_trivial_x_power() {
            "no power of x".to_string()
        } else {
            "a power of x adds or removes only the root 0".to_string()
        },
        conclusion = format_root_interval_exact(&region),
        checked = if checked_rows.is_empty() {
            String::new()
        } else {
            format!(" (rows {checked_rows:?} checked exactly)")
        }
    );
    Ok(NormalizedRealRootednessTransfer {
        schema: NORMALIZED_REAL_ROOTEDNESS_TRANSFER_SCHEMA.into(),
        equality: equality.clone(),
        direction,
        premise_from_index,
        premise_region: premise_region.clone(),
        conclusion_from_index: conclusion_from,
        conclusion_region: region,
        checked_rows,
        statement,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn falling_factorials_and_shifts_are_exact() {
        // (m)_2 with m = 2n + 1: (2n + 1)(2n) = 4n^2 + 2n.
        assert_eq!(
            falling_factorial(2, 1, 2),
            vec![q_int(0), q_int(2), q_int(4)]
        );
        // (n + 1)^2 at n - 1 is n^2.
        let p = vec![q_int(1), q_int(2), q_int(1)];
        assert_eq!(
            shift_index_polynomial(&p, -1),
            vec![q_int(0), q_int(0), q_int(1)]
        );
    }

    #[test]
    fn ratios_are_fitted_with_minimal_degree() {
        let points: Vec<(i64, Q)> = (2..9)
            .map(|n| (n, Q::new(n.into(), (n - 1).into())))
            .collect();
        let (p, q) = fit_ratio(&points, 2).unwrap();
        assert_eq!(p, vec![q_int(0), q_int(1)]);
        assert_eq!(q, vec![q_int(-1), q_int(1)]);
        let constant: Vec<(i64, Q)> = (1..5).map(|n| (n, q_int(-2))).collect();
        assert_eq!(
            fit_ratio(&constant, 2).unwrap(),
            (vec![q_int(-2)], vec![q_int(1)])
        );
    }

    #[test]
    fn the_hull_with_zero_moves_one_end() {
        let region = RootInterval::from_real(&crate::real_rootedness::RealInterval {
            lower: Some(q_int(-3)),
            lower_closed: true,
            upper: Some(q_int(-1)),
            upper_closed: false,
        });
        let hull = hull_with_zero(&region).unwrap();
        assert_eq!(hull.upper, Some(RootEndpoint::rational(&Q::zero(), true)));
        assert_eq!(hull.lower, region.lower);
        let ray =
            RootInterval::from_real(&crate::real_rootedness::RealInterval::at_most(Q::zero()));
        assert_eq!(hull_with_zero(&ray).unwrap(), ray);
    }
}
