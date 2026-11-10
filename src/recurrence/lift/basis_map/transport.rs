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
        .filter_map(Biv::variable_degree)
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
                        system.entries[b][a].add_term((*shift, *d), coefficient.multiply(c));
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

/// How a [`Transport`] was obtained.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum TransportKind {
    /// The companion transport of an index-free map (`Phi = sum_i e_i(n, C) S^i`).
    Companion,
    /// The index-side transport through a Riordan coefficient family
    /// (`super::riordan`): a scalar rule of the image.
    Riordan,
}

/// The transported system with what the identity certificate needs: its
/// leader `g(n) S^k I`, the first base from which `Phi W = 0` holds, and
/// the number of image rows from the first index that determine `H` by
/// `Phi` (used by the reverse division).
#[derive(Debug, Clone)]
pub(super) struct Transport {
    pub(super) system: OperatorMatrix,
    pub(super) order: usize,
    pub(super) leader: Biv,
    pub(super) valid_from_base: usize,
    pub(super) identifying_rows: usize,
    pub(super) kind: TransportKind,
}

/// The transported system of the coefficient family through `map`: the
/// companion transport when the map has a companion form (unchanged), else
/// the index-side transport when the coefficient family is a Riordan
/// family.
///
/// # Errors
///
/// [`RecurrenceError::Unsupported`] when neither applies (with both
/// reasons), [`RecurrenceError::LimitExceeded`] over the term budget, and
/// the errors of the constructions.
pub(super) fn transport_coefficient_family(
    coefficient: &crate::recurrence::lift::hadamard_factor::prepare::PreparedLiftRule,
    map: &BasisMap,
    max_terms: usize,
) -> Result<Transport, RecurrenceError> {
    use crate::recurrence::lift::section_certify::x_free_leader;
    let companion_reason = match map.companion_operators() {
        Ok(_) => {
            let system = transported_system(&coefficient.relation, map, max_terms)?;
            let (order, leader) = x_free_leader(&coefficient.relation)?;
            return Ok(Transport {
                system,
                order,
                leader,
                valid_from_base: coefficient.valid_from_base(),
                identifying_rows: coefficient.seeds.len(),
                kind: TransportKind::Companion,
            });
        }
        Err(RecurrenceError::Unsupported(reason)) => reason,
        Err(error) => return Err(error),
    };
    if let Some(transport) =
        super::even_derivative::exponential_riordan_transport(coefficient, map)?
    {
        check_term_budget(transport.system.term_count(), max_terms)?;
        return Ok(transport);
    }
    let transport = match super::riordan::riordan_transport(coefficient, map) {
        Ok(transport) => transport,
        Err(RecurrenceError::Unsupported(reason)) => {
            return Err(RecurrenceError::Unsupported(format!(
                "{companion_reason}; the index-side transport does not apply either: {reason}"
            )))
        }
        Err(error) => return Err(error),
    };
    check_term_budget(transport.relation.term_count(), max_terms)?;
    let (order, leader) = x_free_leader(&transport.relation)?;
    // H is determined by Phi from the rows before valid_from_base + order and
    // after every integer root of the leader at a base >= valid_from_base.
    let leader_values: Vec<_> = leader
        .rows()
        .iter()
        .map(|row| row.first().cloned().unwrap_or_else(num_traits::Zero::zero))
        .collect();
    let last_root = crate::recurrence::implicit::algebra::roots_in_domain(
        &leader_values,
        transport.valid_from_base,
        None,
    )
    .into_iter()
    .max();
    let mut identifying = transport.valid_from_base + order - coefficient.first_index;
    if let Some(root) = last_root {
        identifying = identifying.max(root + order + 1 - coefficient.first_index);
    }
    Ok(Transport {
        system: OperatorMatrix {
            entries: vec![vec![transport.relation]],
        },
        order,
        leader,
        valid_from_base: transport.valid_from_base,
        identifying_rows: identifying,
        kind: TransportKind::Riordan,
    })
}
