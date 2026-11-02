//! Transport of a rule for a coefficient family `R_n(u)` through a basis
//! map: the operator system `Phi = sum_i e_i(n, C) S^i` satisfied by the
//! state vectors `V(R_n) = (T(R_n), T(u R_n), ..., T(u^(r-1) R_n))`
//! (module documentation of `basis_map`, "Transport").

use super::super::section_operators::{check_term_budget, OperatorMatrix};
use super::map::BasisMap;
use crate::recurrence::implicit::algebra::{Biv, DOp};
use crate::recurrence::{ImplicitRelation, Recurrence, RecurrenceError};

/// `e(n)`: the coefficient of `u^j` of `e(n, u)`, as a polynomial in `n`
/// (a [`Biv`] free of the second variable).
fn u_coefficient(e: &Biv, j: usize) -> Biv {
    Biv::new(
        e.rows()
            .iter()
            .map(|row| vec![row.get(j).cloned().unwrap_or_else(num_traits::Zero::zero)])
            .collect(),
    )
}

/// `Phi = sum_i e_i(n, C) S^i` for the relation `sum_i e_i(n, u) S^i` of
/// the coefficient family (derivative free) and the companion matrix `C`
/// of `map`.
///
/// # Errors
///
/// [`RecurrenceError::Unsupported`] when the relation has a derivative in
/// `u` or the map has no companion form, and
/// [`RecurrenceError::LimitExceeded`] over the term budget.
pub(super) fn transported_system(
    relation: &DOp,
    map: &BasisMap,
    max_terms: usize,
) -> Result<OperatorMatrix, RecurrenceError> {
    if relation.max_derivative() > 0 {
        return Err(RecurrenceError::Unsupported(
            "the coefficient family's rule has derivatives in u; only multiplication by u is \
             intertwined with a fixed operator"
                .into(),
        ));
    }
    let companion = map.companion_matrix()?;
    let size = companion.size();
    let max_power = relation
        .0
        .values()
        .filter_map(Biv::y_degree)
        .max()
        .unwrap_or(0);
    let mut powers = vec![OperatorMatrix::identity(size)];
    for _ in 0..max_power {
        let next = powers.last().expect("nonempty").compose(&companion);
        check_term_budget(next.term_count(), max_terms)?;
        powers.push(next);
    }
    let mut system = OperatorMatrix::zero(size);
    for ((shift, _), e) in &relation.0 {
        for (j, power) in powers.iter().enumerate() {
            let coefficient = u_coefficient(e, j);
            if coefficient.is_zero() {
                continue;
            }
            for b in 0..size {
                for a in 0..size {
                    for ((_, d), c) in &power.entries[b][a].0 {
                        system.entries[b][a].add_term((*shift, *d), coefficient.mul(c));
                    }
                }
            }
        }
        check_term_budget(system.term_count(), max_terms)?;
    }
    Ok(system)
}

/// For a map of order one (`C = [A]`), the scalar transported rule
/// `f(m) H_m = ...` (target index `m`) of `Phi`, with its order.
///
/// # Errors
///
/// [`RecurrenceError::Unsupported`] for a map of higher order or a system
/// without a derivative-free leader.
pub(super) fn transported_scalar_rule(
    system: &OperatorMatrix,
) -> Result<(Recurrence, usize), RecurrenceError> {
    if system.size() != 1 {
        return Err(RecurrenceError::Unsupported(
            "only a map of order one gives a scalar transported rule directly".into(),
        ));
    }
    ImplicitRelation::from_op(&system.entries[0][0])
        .to_successor_rule()
        .ok_or_else(|| {
            RecurrenceError::Unsupported(
                "the transported relation has no derivative-free leader".into(),
            )
        })
}
