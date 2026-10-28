//! Shared input handling of the Hadamard-factor and row-factor lifts:
//! validation of a scalar rule with its identifying rows, exact replay, the
//! leader nonvanishing check that makes the rows unique, and the JSON form of
//! shift-free operators.

use crate::recurrence::implicit::algebra::DOp;
use crate::recurrence::implicit::{biv_from_json, biv_json};
use crate::recurrence::root_domain::bivariate_leader_nonvanishing;
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

/// A validated scalar rule with its identifying seeds and relation.
#[derive(Clone, Debug)]
pub(in crate::recurrence::lift) struct PreparedRule {
    pub(in crate::recurrence::lift) rule: Recurrence,
    pub(in crate::recurrence::lift) first_index: usize,
    /// The `K` identifying seeds (`Recurrence::generation_initial_count`).
    pub(in crate::recurrence::lift) seeds: Vec<Vec<Q>>,
    /// The rule at base index `n` (`ImplicitRelation::from_recurrence`).
    pub(in crate::recurrence::lift) relation: DOp,
    /// The leader shift `R` (the order).
    pub(in crate::recurrence::lift) leader_shift: usize,
}

impl PreparedRule {
    /// Rows `first_index, ..., first_index + count - 1`, trimmed.
    pub(in crate::recurrence::lift) fn rows(
        &self,
        count: usize,
    ) -> Result<Vec<Vec<Q>>, RecurrenceError> {
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

/// Validate `rule` with `rows` (`rows[0]` is row `first_index`): positive
/// offsets, no `(-1)^n` factor, homogeneous, derivative order at most
/// [`MAX_LIFT_DERIVATIVE_ORDER`], supplied rows reproduced by the rule, and
/// the leader `f(m, .)` a nonzero polynomial at every target
/// `m >= first_index + K` (decided exactly), so that the identifying rows
/// determine every row.
pub(in crate::recurrence::lift) fn prepare_rule(
    rule: &Recurrence,
    rows: &[Vec<Q>],
    first_index: usize,
    name: &str,
) -> Result<PreparedRule, RecurrenceError> {
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
    Ok(PreparedRule {
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

/// `factor o op` for a coefficient `factor` free of `S` and `D`.
pub(in crate::recurrence::lift) fn left_multiply(
    factor: &crate::recurrence::implicit::algebra::Biv,
    op: &DOp,
) -> DOp {
    let mut out = DOp::default();
    for (key, c) in &op.0 {
        out.add_term(*key, factor.mul(c));
    }
    out
}
