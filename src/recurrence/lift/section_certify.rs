//! Pseudo-division of an input relation by a transported parent system, the
//! gauge and index shift of a candidate parent, and the forward certificate
//! [`VeroneseForwardDivisionCertificate`].
//!
//! # Pseudo-division
//!
//! Let `Phi = phi(L)` be the transported parent system (an `r x r` operator
//! matrix, `section_operators`) of a parent relation `L` whose leader is
//! `g(n) S^k` with `g` free of `x` and of derivatives, so that the shift-`k`
//! part of `Phi` is `g(n) S^k I`.  A row `rho` of `r` operators is divided as
//! follows: while the largest shift `m` of the remainder is at least `k`,
//! take the coefficient row `c` of `S^m` (shift-free operators), multiply the
//! remainder and all recorded quotients on the left by `g(n + m - k)`,
//! subtract `c o S^(m-k) o Phi` and record `c` as the quotient of
//! `S^(m-k)`.  The shift-`m` terms cancel because `c` commutes with
//! `g(n + m - k)` (both are free of `S`, and `g` is free of `y` and `D_y`).
//! The result is the identity
//!
//! ```text
//! w(n) rho = sum_nu H_nu o S^nu o Phi + R
//! ```
//!
//! with `w` the product of the multipliers (a polynomial in `n` only), shift
//! -free rows `H_nu` and a remainder `R` of shift `< k`.
//!
//! # The forward certificate
//!
//! For the forward certificate `rho = e_j^T M` with `M` the input relation,
//! so `rho F = M F_j` for a section vector `F`.  When `R = 0`,
//! [`VeroneseForwardDivisionCertificate::verify`] recomputes the identity
//! with exact operator arithmetic; what it proves is stated in
//! `certify_veronese_section_lift`.
//!
//! # Gauge and shift
//!
//! A candidate `Q_n = lambda(n) V_(r, j)(P_(n+s))` uses the gauged parent
//! `P~_n = lambda(n) P_(n+s)`.  With `rho(n) = lambda(n)/lambda(n-1) =
//! a(n)/b(n)`, a parent relation `sum_(i,d) c_(i,d)(n, x) D^d P_(n+i) = 0`
//! with leader shift `k` gives, after the substitution `n -> n + s`,
//! multiplication by `lambda(n+k) prod_(t=1..k) b(n+t)` and
//! `lambda(n+k)/lambda(n+i) = prod_(t=i+1..k) a(n+t)/b(n+t)`,
//!
//! ```text
//! sum_(i,d) c_(i,d)(n+s, x) prod_(t=i+1..k) a(n+t) prod_(t=1..i) b(n+t) D^d P~_(n+i) = 0,
//! ```
//!
//! valid wherever the original relation holds at `n + s` and `lambda` is
//! finite and nonzero at `n, ..., n+k`.  Its leader
//! `g(n+s) prod_(t=1..k) b(n+t)` is still free of `x`.

use super::super::implicit::algebra::{q, Biv, DOp};
use super::super::implicit::{biv_from_json, biv_json};
use super::section_operators::{check_term_budget, transport_to_sections, OperatorMatrix};
use super::*;
use crate::integer_range_sign::{
    decide_index_polynomial_sign, IndexSignCondition, IndexSignDecision, IntegerRange,
};
use crate::vec_poly;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

type Q = BigRational;

/// The normalisation `lambda(n)` of a candidate `Q_n = lambda(n) V_(r, j)(P_(n+s))`.
///
/// Polynomials in `n` are ascending exact rational strings (integers or
/// `num/den`).  Every form requires `lambda` finite and nonzero at every
/// index `n >= n0` (the input's first index); this is decided exactly by
/// [`crate::decide_index_polynomial_sign`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum VeroneseSectionGauge {
    /// `lambda(n) = 1`.
    None,
    /// `lambda(n0) = initial_value` and `lambda(n)/lambda(n-1) =
    /// numerator(n)/denominator(n)` for `n > n0`; both polynomials must be
    /// nonzero at every integer `n > n0`.  For example `lambda(n) = 1/n!`
    /// from `n0 = 0` is `numerator = [1]`, `denominator = [0, 1]`,
    /// `initial_value = 1`, and `lambda(n) = 2^n` is `numerator = [2]`.
    Ratio {
        /// `a(n)`, ascending coefficients.
        numerator: Vec<String>,
        /// `b(n)`, ascending coefficients.
        denominator: Vec<String>,
        /// `lambda(n0)`, nonzero.
        initial_value: String,
    },
    /// `lambda(n) = numerator(n)/denominator(n)`; both polynomials must be
    /// nonzero at every integer `n >= n0`.
    Closed {
        /// Ascending coefficients.
        numerator: Vec<String>,
        /// Ascending coefficients.
        denominator: Vec<String>,
    },
}

/// A gauge in ratio form: `lambda(n0) = initial` and
/// `lambda(n)/lambda(n-1) = numerator(n)/denominator(n)` for `n > n0`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct RatioGauge {
    pub(super) numerator: Vec<Q>,
    pub(super) denominator: Vec<Q>,
    pub(super) initial: Q,
}

fn parse_polynomial(values: &[String], name: &str) -> Result<Vec<Q>, RecurrenceError> {
    let parsed = values
        .iter()
        .map(|value| {
            parse_rational_coeff(value)
                .map_err(|error| RecurrenceError::from(error).with_context(format!("gauge {name}")))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let trimmed = vec_poly::trim(&parsed);
    if trimmed.iter().all(Zero::is_zero) {
        return Err(RecurrenceError::InvalidInput(format!(
            "the gauge {name} is the zero polynomial"
        )));
    }
    Ok(trimmed)
}

/// `p(n) != 0` for every integer `n >= from`, decided exactly
/// ([`crate::decide_index_polynomial_sign`]).
fn require_no_roots_from(p: &[Q], from: usize, what: &str) -> Result<(), RecurrenceError> {
    match decide_index_polynomial_sign(
        p,
        &IntegerRange::starting_at(from),
        IndexSignCondition::Nonzero,
    ) {
        IndexSignDecision::Holds { .. } => Ok(()),
        IndexSignDecision::FailsAt { index, .. } => Err(RecurrenceError::InvalidInput(format!(
            "the gauge {what} vanishes at n = {index}; it must be nonzero for every n >= {from}"
        ))),
        IndexSignDecision::Undecided { reason } => Err(RecurrenceError::BudgetExhausted(format!(
            "the gauge {what} could not be checked to be nonzero for every n >= {from}: {reason}"
        ))),
    }
}

impl VeroneseSectionGauge {
    /// The ratio form for first index `first_index`, after checking that
    /// `lambda` is finite and nonzero at every `n >= first_index`.
    pub(super) fn ratio_form(&self, first_index: usize) -> Result<RatioGauge, RecurrenceError> {
        match self {
            Self::None => Ok(RatioGauge {
                numerator: vec![Q::one()],
                denominator: vec![Q::one()],
                initial: Q::one(),
            }),
            Self::Ratio {
                numerator,
                denominator,
                initial_value,
            } => {
                let numerator = parse_polynomial(numerator, "numerator")?;
                let denominator = parse_polynomial(denominator, "denominator")?;
                let initial = parse_rational_coeff(initial_value)
                    .map_err(|error| RecurrenceError::from(error).with_context("gauge initial"))?;
                if initial.is_zero() {
                    return Err(RecurrenceError::InvalidInput(
                        "the gauge initial value lambda(n0) must be nonzero".into(),
                    ));
                }
                require_no_roots_from(&numerator, first_index + 1, "ratio numerator")?;
                require_no_roots_from(&denominator, first_index + 1, "ratio denominator")?;
                Ok(RatioGauge {
                    numerator,
                    denominator,
                    initial,
                })
            }
            Self::Closed {
                numerator,
                denominator,
            } => {
                let numerator = parse_polynomial(numerator, "numerator")?;
                let denominator = parse_polynomial(denominator, "denominator")?;
                require_no_roots_from(&numerator, first_index, "numerator")?;
                require_no_roots_from(&denominator, first_index, "denominator")?;
                let at = |p: &[Q], n: i64| vec_poly::evaluate(p, &q(n));
                let n0 = i64::try_from(first_index)
                    .map_err(|_| RecurrenceError::LimitExceeded("first index overflow".into()))?;
                let initial = at(&numerator, n0) / at(&denominator, n0);
                // lambda(n)/lambda(n-1) = a(n) b(n-1) / (b(n) a(n-1)).
                let previous = |p: &[Q]| vec_poly::compose(p, &[q(-1), Q::one()]);
                Ok(RatioGauge {
                    numerator: vec_poly::mul(&numerator, &previous(&denominator)),
                    denominator: vec_poly::mul(&denominator, &previous(&numerator)),
                    initial,
                })
            }
        }
    }
}

impl RatioGauge {
    /// `lambda(first_index), ..., lambda(first_index + count - 1)`.
    pub(super) fn values(&self, first_index: usize, count: usize) -> Vec<Q> {
        let mut values = Vec::with_capacity(count);
        let mut current = self.initial.clone();
        for offset in 0..count {
            if offset > 0 {
                let n = Q::from_integer((first_index + offset).into());
                current = current * vec_poly::evaluate(&self.numerator, &n)
                    / vec_poly::evaluate(&self.denominator, &n);
            }
            values.push(current.clone());
        }
        values
    }
}

/// A polynomial in `n` as a [`Biv`] without `y`.
pub(super) fn n_polynomial(p: &[Q]) -> Biv {
    Biv::new(p.iter().map(|c| vec![c.clone()]).collect())
}

/// The coefficients of a [`Biv`] free of `y` as a polynomial in `n`.
pub(super) fn n_coefficients(p: &Biv) -> Vec<Q> {
    vec_poly::trim(
        &p.rows()
            .iter()
            .map(|row| row.first().cloned().unwrap_or_else(Q::zero))
            .collect::<Vec<_>>(),
    )
}

/// The leader `g(n)` of a relation with leader shift `k`, required to be a
/// single derivative-free term free of `x`.
pub(super) fn x_free_leader(relation: &DOp) -> Result<(usize, Biv), RecurrenceError> {
    let k = relation.max_shift();
    if k == 0 {
        return Err(RecurrenceError::Unsupported(
            "the parent relation has no positive shift".into(),
        ));
    }
    let column = relation.column(k);
    match column.0.get(&(0, 0)) {
        Some(leader) if column.0.len() == 1 && leader.is_y_free() => Ok((k, leader.clone())),
        _ => Err(RecurrenceError::Unsupported(
            "the parent leader must be a derivative-free coefficient g(n) free of x".into(),
        )),
    }
}

/// The gauged and shifted parent relation of the module documentation.
pub(super) fn gauge_parent_relation(
    parent: &DOp,
    index_shift: i64,
    gauge: &RatioGauge,
) -> Result<DOp, RecurrenceError> {
    let (k, _) = x_free_leader(parent)?;
    let numerator = n_polynomial(&gauge.numerator);
    let denominator = n_polynomial(&gauge.denominator);
    let mut factors = Vec::with_capacity(k + 1);
    for i in 0..=k {
        let mut factor = Biv::one();
        for t in i + 1..=k {
            factor = factor.mul(&numerator.shift_n(t as i64));
        }
        for t in 1..=i {
            factor = factor.mul(&denominator.shift_n(t as i64));
        }
        factors.push(factor);
    }
    let mut out = DOp::default();
    for ((i, d), c) in &parent.0 {
        out.add_term((*i, *d), c.shift_n(index_shift).mul(&factors[*i]));
    }
    Ok(out)
}

/// `w rho = sum_nu H_nu o S^nu o Phi + R`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct ForwardDivision {
    /// `w(n)`, free of `y`.
    pub(super) multiplier: Biv,
    /// `nu -> H_nu`, rows of shift-free operators.
    pub(super) quotients: BTreeMap<usize, Vec<DOp>>,
    /// `R`, every entry of shift `< k`.
    pub(super) remainder: Vec<DOp>,
}

impl ForwardDivision {
    pub(super) fn remainder_is_zero(&self) -> bool {
        self.remainder.iter().all(DOp::is_zero)
    }
}

fn left_multiply(factor: &Biv, op: &DOp) -> DOp {
    let mut out = DOp::default();
    for (key, c) in &op.0 {
        out.add_term(*key, factor.mul(c));
    }
    out
}

/// `sum_b row[b] o S^nu o system[b][a]` for every column `a`.
fn row_times_shifted_system(row: &[DOp], nu: usize, system: &OperatorMatrix) -> Vec<DOp> {
    let size = system.size();
    (0..size)
        .map(|a| {
            let mut sum = DOp::default();
            for (b, c) in row.iter().enumerate() {
                if c.is_zero() || system.entries[b][a].is_zero() {
                    continue;
                }
                sum = sum.add(&c.compose(&system.entries[b][a].shifted(nu)));
            }
            sum
        })
        .collect()
}

fn row_term_count(row: &[DOp]) -> usize {
    row.iter().map(|op| op.0.len()).sum()
}

/// Check that the shift-`k` part of `system` is `g(n) S^k I` and nothing has
/// a larger shift.
fn check_system_leader(
    system: &OperatorMatrix,
    k: usize,
    leader: &Biv,
) -> Result<(), RecurrenceError> {
    for (b, row) in system.entries.iter().enumerate() {
        for (a, entry) in row.iter().enumerate() {
            let expected = if a == b {
                DOp::single(0, 0, leader.clone())
            } else {
                DOp::default()
            };
            if entry.max_shift() > k || entry.column(k) != expected {
                return Err(RecurrenceError::VerificationFailed(
                    "the transported parent system does not have the leader g(n) S^k I".into(),
                ));
            }
        }
    }
    Ok(())
}

/// Pseudo-divide `row` by `system` (module documentation).
pub(super) fn pseudo_divide_row(
    row: &[DOp],
    system: &OperatorMatrix,
    k: usize,
    leader: &Biv,
    max_terms: usize,
) -> Result<ForwardDivision, RecurrenceError> {
    check_system_leader(system, k, leader)?;
    let mut remainder = row.to_vec();
    let mut multiplier = Biv::one();
    let mut quotients: BTreeMap<usize, Vec<DOp>> = BTreeMap::new();
    loop {
        if remainder.iter().all(DOp::is_zero) {
            break;
        }
        let top = remainder.iter().map(DOp::max_shift).max().unwrap_or(0);
        if top < k {
            break;
        }
        let nu = top - k;
        let factor = leader.shift_n(nu as i64);
        let coefficient_row: Vec<DOp> = remainder.iter().map(|op| op.column(top)).collect();
        let subtracted = row_times_shifted_system(&coefficient_row, nu, system);
        remainder = remainder
            .iter()
            .zip(&subtracted)
            .map(|(r, s)| left_multiply(&factor, r).sub(s))
            .collect();
        for quotient in quotients.values_mut() {
            *quotient = quotient.iter().map(|h| left_multiply(&factor, h)).collect();
        }
        let entry = quotients
            .entry(nu)
            .or_insert_with(|| vec![DOp::default(); row.len()]);
        for (h, c) in entry.iter_mut().zip(&coefficient_row) {
            *h = h.add(c);
        }
        multiplier = multiplier.mul(&factor);
        if remainder
            .iter()
            .any(|op| op.max_shift() >= top && !op.column(top).is_zero())
        {
            return Err(RecurrenceError::VerificationFailed(
                "pseudo-division did not cancel the leading shift".into(),
            ));
        }
        let terms = row_term_count(&remainder)
            + quotients.values().map(|q| row_term_count(q)).sum::<usize>();
        check_term_budget(terms, max_terms)?;
    }
    Ok(ForwardDivision {
        multiplier,
        quotients,
        remainder,
    })
}

/// Whether `w rho - sum_nu H_nu o S^nu o Phi - R` vanishes exactly.
pub(super) fn division_identity_holds(
    row: &[DOp],
    system: &OperatorMatrix,
    division: &ForwardDivision,
) -> bool {
    let mut difference: Vec<DOp> = row
        .iter()
        .zip(&division.remainder)
        .map(|(r, remainder)| left_multiply(&division.multiplier, r).sub(remainder))
        .collect();
    for (nu, quotient) in &division.quotients {
        let product = row_times_shifted_system(quotient, *nu, system);
        difference = difference
            .iter()
            .zip(&product)
            .map(|(d, p)| d.sub(p))
            .collect();
    }
    difference.iter().all(DOp::is_zero)
}

/// `e_j^T M`: the row with `M` in column `residue`.
pub(super) fn residue_row(relation: &DOp, modulus: usize, residue: usize) -> Vec<DOp> {
    let mut row = vec![DOp::default(); modulus];
    row[residue] = relation.clone();
    row
}

/// One entry of a quotient row `H_nu` of a [`VeroneseForwardDivisionCertificate`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VeroneseDivisionQuotient {
    /// The shift `nu` of `H_nu o S^nu o Phi`.
    pub shift: usize,
    /// The column `b` of the row `H_nu` (the parent-system row it multiplies).
    pub component: usize,
    /// `H_nu[b] = sum_d h_d(n, y) D_y^d`.
    pub differentiations: Vec<OperatorTermJson>,
}

/// The exact operator identity `w(n) e_j^T M = sum_nu H_nu o S^nu o phi(L~)`
/// (remainder zero) of a certified Veronese-section lift.
///
/// `M` is the input relation (variable `y`), `L~` the gauged and shifted
/// parent relation (variable `x`) and `phi` the section transport at the
/// modulus.  [`verify`](Self::verify) recomputes `L~` from the parent
/// relation, the shift and the gauge, transports it and checks the identity
/// with exact Ore-algebra arithmetic.  The identity alone shows that
/// `w(n) (M V_j(P~))_n = 0` wherever the parent relation holds at the bases
/// `n + nu`; the index domains, the exceptional bases (integer roots of `w`
/// and bases below the parent's domain) and the identifying rows recorded
/// here complete the all-index argument of `certify_veronese_section_lift`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VeroneseForwardDivisionCertificate {
    /// The modulus `r`.
    pub modulus: usize,
    /// The residue `j`.
    pub residue: usize,
    /// `M`: the input rule as a relation at base index `n` (leader shift
    /// `input_leader_shift`), coefficients in `(n, y)`.
    pub input_relation: ImplicitRelationJson,
    /// The largest offset `R` of the input rule.
    pub input_leader_shift: usize,
    /// The parent rule as a relation at base index `n`, coefficients in `(n, x)`.
    pub parent_relation: ImplicitRelationJson,
    /// The index shift `s` of `P~_n = lambda(n) P_(n+s)`.
    pub index_shift: i64,
    /// The gauge `lambda`.
    pub gauge: VeroneseSectionGauge,
    /// `L~`: the gauged and shifted parent relation.
    pub gauged_parent_relation: ImplicitRelationJson,
    /// Its leader shift `k`; its leader `g~(n) S^k` is free of `x`.
    pub parent_leader_shift: usize,
    /// `w(n)`, ascending exact strings in `n`.
    pub multiplier: Vec<String>,
    /// The nonzero entries of the quotient rows `H_nu`.
    pub quotients: Vec<VeroneseDivisionQuotient>,
    /// First base index `n` at which the input rule holds for the input
    /// sequence (`n0 + K - R`, with `K` identifying rows).
    pub input_valid_from_base: usize,
    /// First base index at which `L~` holds for the gauged parent.
    pub parent_valid_from_base: usize,
    /// Base indices where the identity gives nothing (roots of `w`, bases
    /// below `parent_valid_from_base`); `M V_j(P~)` was evaluated exactly
    /// and found zero at each of them.
    pub exceptional_bases: Vec<usize>,
    /// Indices `n0, ..., n0 + K - 1` where `Q_n = V_j(P~_n)` was checked
    /// exactly.
    pub identifying_indices: Vec<usize>,
    /// What the certificate proves.
    pub statement: String,
}

fn shift_free_to_json(op: &DOp) -> Vec<OperatorTermJson> {
    op.0.iter()
        .map(|((_, d), c)| OperatorTermJson {
            deriv_order: *d,
            coeff: biv_json(c),
        })
        .collect()
}

fn shift_free_from_json(terms: &[OperatorTermJson]) -> Result<DOp, RecurrenceError> {
    let mut op = DOp::default();
    for term in terms {
        op.add_term((0, term.deriv_order), biv_from_json(&term.coeff)?);
    }
    Ok(op)
}

impl VeroneseForwardDivisionCertificate {
    pub(super) fn quotient_entries(division: &ForwardDivision) -> Vec<VeroneseDivisionQuotient> {
        let mut out = Vec::new();
        for (nu, row) in &division.quotients {
            for (component, op) in row.iter().enumerate() {
                if !op.is_zero() {
                    out.push(VeroneseDivisionQuotient {
                        shift: *nu,
                        component,
                        differentiations: shift_free_to_json(op),
                    });
                }
            }
        }
        out
    }

    /// Recompute the gauged parent relation and the transported system, and
    /// check `w e_j^T M = sum_nu H_nu o S^nu o phi(L~)` exactly.
    ///
    /// This replays the operator identity only; the index-domain argument
    /// (seeds, exceptional bases, nonvanishing leaders) is documented in
    /// `certify_veronese_section_lift` and recorded in the other fields.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        let failed = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
        if !(2..=MAX_VERONESE_CERTIFY_MODULUS).contains(&self.modulus)
            || self.residue >= self.modulus
        {
            return failed("the certificate has an invalid modulus or residue");
        }
        if self.index_shift.unsigned_abs() > MAX_VERONESE_CERTIFY_INDEX_SHIFT as u64 {
            return failed("the certificate index shift is out of range");
        }
        let input = ImplicitRelation::from_json(&self.input_relation)?.op();
        let parent = ImplicitRelation::from_json(&self.parent_relation)?.op();
        let gauged = ImplicitRelation::from_json(&self.gauged_parent_relation)?.op();
        if input.max_shift() != self.input_leader_shift
            || input.column(self.input_leader_shift).max_derivative() != 0
        {
            return failed("the input relation does not have the recorded derivative-free leader");
        }
        if input.max_derivative() > MAX_VERONESE_CERTIFY_DERIVATIVE_ORDER
            || parent.max_derivative() > MAX_VERONESE_CERTIFY_DERIVATIVE_ORDER
        {
            return failed("the certificate exceeds the derivative-order bound");
        }
        let Some(&first_index) = self.identifying_indices.first() else {
            return failed("the certificate lists no identifying indices");
        };
        let ratio = self.gauge.ratio_form(first_index)?;
        let recomputed = gauge_parent_relation(&parent, self.index_shift, &ratio)?;
        if recomputed != gauged {
            return failed("the gauged parent relation does not match its recomputation");
        }
        let (k, leader) = x_free_leader(&gauged)?;
        if k != self.parent_leader_shift {
            return failed("the parent leader shift does not match");
        }
        let system = transport_to_sections(&gauged, self.modulus, MAX_VERONESE_OPERATOR_TERMS)?;
        let row = residue_row(&input, self.modulus, self.residue);
        let multiplier = n_polynomial(&decode_poly(&self.multiplier)?);
        let mut quotients: BTreeMap<usize, Vec<DOp>> = BTreeMap::new();
        for entry in &self.quotients {
            if entry.component >= self.modulus {
                return failed("a quotient component is out of range");
            }
            let row = quotients
                .entry(entry.shift)
                .or_insert_with(|| vec![DOp::default(); self.modulus]);
            row[entry.component] =
                row[entry.component].add(&shift_free_from_json(&entry.differentiations)?);
        }
        let division = ForwardDivision {
            multiplier,
            quotients,
            remainder: vec![DOp::default(); self.modulus],
        };
        if division.multiplier.is_zero() {
            return failed("the multiplier w is zero");
        }
        check_system_leader(&system, k, &leader)?;
        if division_identity_holds(&row, &system, &division) {
            Ok(())
        } else {
            failed("the forward-division identity w e_j^T M = sum H_nu S^nu phi(L~) fails")
        }
    }
}
