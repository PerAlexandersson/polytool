//! Shared input handling of the Hadamard-factor and row-factor lifts:
//! validation of a scalar rule with its identifying rows, exact replay, the
//! leader nonvanishing check that makes the rows unique, and the JSON form of
//! shift-free operators.
//!
//! A rule with a left-hand factor `x^e f(n)` (`e >= 1`, `f` free of `x`)
//! whose `x^e` divides every right-hand coefficient is first replaced by the
//! exactly equivalent rule `f(n) P_n = sum (c / x^e) D^d P_(n-r)` (same rows;
//! [`crate::recurrence::LeftHandFactorNormalization`]), so the prepared rule
//! and its relation have a leader free of `x`; the Veronese lift's
//! preparation does the same.  A factor `x^e f(n)` that does not divide out
//! is rejected with the reason.  A replay re-runs the preparation on the
//! recorded rule and so re-derives the normalization.

use crate::recurrence::implicit::algebra::DOp;
use crate::recurrence::implicit::{biv_from_json, biv_json};
use crate::recurrence::lift::parent_theorems::rule_with_x_free_left_hand_factor;
use crate::recurrence::root_domain::bivariate_leader_nonvanishing;
use crate::recurrence::LeftHandFactorNormalization;
use crate::recurrence::{
    BivarPoly, ImplicitRelation, OperatorTermJson, Recurrence, RecurrenceError, RecurrenceTerm,
};
use crate::vec_poly;
use crate::BigRational;

type Q = BigRational;

/// Largest derivative order of a rule accepted by the two lifts.
pub(in crate::recurrence::lift) const MAX_LIFT_DERIVATIVE_ORDER: usize = 6;

/// Largest number of identifying rows of a rule.
pub(in crate::recurrence::lift) const MAX_LIFT_IDENTIFYING_ROWS: usize = 64;

/// A scalar rule validated for the lifts ([`prepare_lift_rule`]), with its
/// identifying seeds: the rows `first_index, ..., first_index + K - 1` that
/// determine every later row, because the leader of the rule is nonzero at
/// every later target index (decided exactly).
///
/// `rule` is the normalized rule (a left-hand factor `x^e f(n)` divided out
/// when `x^e` divides the right-hand side); it has the same rows as the rule
/// that was prepared.
#[derive(Clone, Debug)]
pub struct PreparedLiftRule {
    /// The validated, normalized rule.
    pub rule: Recurrence,
    /// The index of the first identifying seed.
    pub first_index: usize,
    /// The `K` identifying seeds (`Recurrence::generation_initial_count`),
    /// trimmed; `seeds[0]` is row `first_index`.
    pub seeds: Vec<Vec<Q>>,
    /// The rule at base index `n` (`ImplicitRelation::from_recurrence`).
    pub(crate) relation: DOp,
    /// The leader shift `R` (the order of the rule).
    pub leader_shift: usize,
}

impl PreparedLiftRule {
    /// Rows `first_index, ..., first_index + count - 1`, trimmed.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::InvalidInput`] when the rule cannot be replayed
    /// from its seeds.
    pub fn rows(&self, count: usize) -> Result<Vec<Vec<Q>>, RecurrenceError> {
        self.rule
            .generate_rows_rational(&self.seeds, self.first_index, count.max(self.seeds.len()))
            .map(|rows| {
                rows.iter()
                    .take(count)
                    .map(|row| vec_poly::trim(row))
                    .collect()
            })
            .map_err(|error| {
                RecurrenceError::InvalidInput(format!("the recurrence cannot be replayed: {error}"))
            })
    }

    /// The first base index at which the relation holds for the generated
    /// rows: `first_index + K - R`.
    pub(in crate::recurrence::lift) fn valid_from_base(&self) -> usize {
        self.first_index + self.seeds.len() - self.leader_shift
    }
}

/// `rule` with a left-hand factor `x^e f(n)` normalized to `f(n)` when `x^e`
/// divides every right-hand coefficient (module documentation), with the
/// normalization record; `rule` itself otherwise.
///
/// # Errors
///
/// [`RecurrenceError::Unsupported`] naming the `name` recurrence when the
/// left-hand factor is `x^e f(n)` but `x^e` does not divide the right-hand
/// side.
pub(in crate::recurrence::lift) fn x_free_left_hand_factor_rule(
    rule: &Recurrence,
    name: &str,
) -> Result<(Recurrence, Option<LeftHandFactorNormalization>), RecurrenceError> {
    rule_with_x_free_left_hand_factor(rule).map_err(|error| match error {
        RecurrenceError::Unsupported(reason) => {
            RecurrenceError::Unsupported(format!("the {name} recurrence: {reason}"))
        }
        other => other,
    })
}

/// Validate `rule` with `rows` (`rows[0]` is row `first_index`), after
/// normalizing a left-hand factor `x^e f(n)` whose `x^e` divides every
/// right-hand coefficient to `f(n)` (the prepared rule is the normalized
/// one; [`crate::recurrence::LeftHandFactorNormalization`]): positive
/// offsets, no `(-1)^n` factor, homogeneous, derivative order at most 6, at
/// most 64 identifying rows, supplied rows reproduced by the rule, and
/// the leader `f(m, .)` a nonzero polynomial at every target
/// `m >= first_index + K` (decided exactly), so that the identifying rows
/// determine every row.
///
/// `rows` must contain at least the identifying rows; `name` names the rule
/// in error messages (for example `"parent"`).
///
/// # Errors
///
/// [`RecurrenceError::Unsupported`], [`RecurrenceError::InvalidInput`] or
/// [`RecurrenceError::LimitExceeded`] naming the first condition above that
/// fails.
pub fn prepare_lift_rule(
    rule: &Recurrence,
    rows: &[Vec<Q>],
    first_index: usize,
    name: &str,
) -> Result<PreparedLiftRule, RecurrenceError> {
    let (rule, _) = x_free_left_hand_factor_rule(rule, name)?;
    let rule = &rule;
    if rule.terms.is_empty() {
        return Err(RecurrenceError::Unsupported(format!(
            "the {name} recurrence has no terms"
        )));
    }
    if let Some(position) = rule.terms.iter().position(|t| t.offset == 0) {
        return Err(RecurrenceError::InvalidInput(format!(
            "{name} recurrence term {position} has offset 0; offsets must be positive"
        )));
    }
    if rule.terms.iter().any(RecurrenceTerm::is_alternating) {
        return Err(RecurrenceError::Unsupported(format!(
            "the {name} recurrence has a (-1)^n factor"
        )));
    }
    if rule.inhomogeneous.as_ref().is_some_and(|g| !g.is_zero()) {
        return Err(RecurrenceError::Unsupported(format!(
            "the {name} recurrence has an inhomogeneous term"
        )));
    }
    if let Some(term) = rule
        .terms
        .iter()
        .find(|t| t.deriv_order > MAX_LIFT_DERIVATIVE_ORDER)
    {
        return Err(RecurrenceError::LimitExceeded(format!(
            "the {name} recurrence has derivative order {} > {MAX_LIFT_DERIVATIVE_ORDER}",
            term.deriv_order
        )));
    }
    if rule.denominator.as_ref().is_some_and(BivarPoly::is_zero) {
        return Err(RecurrenceError::InvalidInput(format!(
            "the {name} recurrence denominator is zero"
        )));
    }
    let order = rule.max_offset();
    if rows.len() < order {
        return Err(RecurrenceError::InvalidInput(format!(
            "the {name} recurrence has order {order} but only {} initial rows were supplied",
            rows.len()
        )));
    }
    let count = rule.generation_initial_count(first_index, rows.len());
    if count > MAX_LIFT_IDENTIFYING_ROWS {
        return Err(RecurrenceError::LimitExceeded(format!(
            "the {name} recurrence needs {count} identifying rows, more than \
             {MAX_LIFT_IDENTIFYING_ROWS}"
        )));
    }
    let normalized: Vec<Vec<Q>> = rows.iter().map(|row| vec_poly::trim(row)).collect();
    let seeds = normalized[..count].to_vec();
    let generated = rule
        .generate_rows_rational(&seeds, first_index, normalized.len())
        .map_err(|error| {
            RecurrenceError::InvalidInput(format!(
                "the {name} recurrence cannot be replayed: {error}"
            ))
        })?;
    for (offset, (supplied, generated)) in normalized.iter().zip(&generated).enumerate() {
        if supplied != &vec_poly::trim(generated) {
            return Err(RecurrenceError::InvalidInput(format!(
                "supplied {name} row {} disagrees with the row generated by its recurrence",
                first_index + offset
            )));
        }
    }
    let leader = rule.denominator.clone().unwrap_or(BivarPoly {
        coeffs: vec![vec![Q::from_integer(1.into())]],
    });
    if !bivariate_leader_nonvanishing(&leader, first_index + count) {
        return Err(RecurrenceError::Unsupported(format!(
            "the {name} leader is not proved to be a nonzero polynomial at every target index \
             >= {}, so the identifying rows need not determine the sequence",
            first_index + count
        )));
    }
    let (relation, leader_shift) = ImplicitRelation::from_recurrence(rule)?;
    Ok(PreparedLiftRule {
        rule: rule.clone(),
        first_index,
        seeds,
        relation: relation.op(),
        leader_shift,
    })
}

/// JSON form of a shift-free operator.
pub(in crate::recurrence::lift) fn shift_free_json(op: &DOp) -> Vec<OperatorTermJson> {
    op.0.iter()
        .map(|((_, d), c)| OperatorTermJson {
            deriv_order: *d,
            coeff: biv_json(c),
        })
        .collect()
}

/// The shift-free operator of its JSON form.
pub(in crate::recurrence::lift) fn shift_free_from_json(
    terms: &[OperatorTermJson],
) -> Result<DOp, RecurrenceError> {
    let mut op = DOp::default();
    for term in terms {
        if term.deriv_order > MAX_LIFT_DERIVATIVE_ORDER * 4 {
            return Err(RecurrenceError::LimitExceeded(
                "an operator term exceeds the derivative-order bound".into(),
            ));
        }
        op.add_term((0, term.deriv_order), biv_from_json(&term.coeff)?);
    }
    Ok(op)
}
