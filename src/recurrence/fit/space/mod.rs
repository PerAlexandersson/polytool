//! Complete bounded recurrence spaces and replayable exact certificates.
//!
//! Coefficients are evaluated in the original index, never silently shifted.
//! A supplied monic constant-index source gives all-index identities via the
//! `h*k` residual bound. Without a source, positives remain finite fits, but
//! dual witnesses can still prove absence in the requested support.
//!
//! # Finite certificate bound
//!
//! Write `chi = S^k - sum_i a_i(X) S^i` for the monic source operator.
//! Since its coefficients are independent of n, `chi` commutes with shifts,
//! multiplication by X, and its own coefficient derivatives. Differentiating
//! `chi P = 0` gives `chi D(P) = -D(chi) P`; induction gives
//! `chi^(s+1) D^s(P) = 0`. The commutator with multiplication by n is
//! `S (d chi / d S)`, which commutes with chi. Thus multiplying by n raises
//! the required power by at most one. Every column `X^v n^u S^j D^s(P)` is
//! annihilated by `chi^(u+s+1)`; their sum is annihilated by `chi^h`.
//! This is a monic recurrence of order h*k. Its first h*k zero values on
//! the declared domain therefore imply all subsequent values are zero.
//! No factorization, distinct-root assumption, or invertible oldest
//! coefficient is required. A source inferred from rows remains an assumption.
//!
//! With the opt-in `coefficient_rescaling` field the same machinery runs on
//! rescaled rows; the `fit::rescaling` module documents its residual bound.

use super::{
    parse_rational_coeff, BigRational as Q, BivarPoly, Recurrence, RecurrenceCoefficientRescaling,
    RecurrenceJson, RecurrenceRescalingBound, RecurrenceSign, RecurrenceTerm,
};
use crate::linalg::{rational_kernel_with_certificate, RationalKernelError, RationalKernelLimits};
use crate::recurrence::RecurrenceError;
use crate::{sturm::SturmChain, Polynomial};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

mod propagation;
pub use propagation::{IndexedLagSource, RecurrencePropagationCertificate};

mod types;
pub use types::*;

type Poly = Polynomial<Q>;
fn is_false(value: &bool) -> bool {
    !*value
}

fn encode_polys(polys: &[Poly]) -> ExactRows {
    polys
        .iter()
        .map(|p| p.coeffs().iter().map(ToString::to_string).collect())
        .collect()
}
pub(in crate::recurrence) fn encode_matrix(matrix: &[Vec<Q>]) -> ExactRows {
    matrix
        .iter()
        .map(|r| r.iter().map(ToString::to_string).collect())
        .collect()
}

/// Return every identity in each requested finite coefficient support.
/// Positive certificates are conditional on the source if supplied; otherwise
/// they describe finite fits. Negative dual certificates need only the rows.
pub fn find_recurrence_spaces(
    request: &RecurrenceSpaceRequest,
) -> Result<RecurrenceSpaces, RecurrenceError> {
    let input = prepare(request)?;
    let symbol = row_symbol(input.rescaling.is_some());
    // A contradictory source is invalid input, not an incomplete search.
    source_rows(&input, input.rows.len(), &request.limits)?;
    let mut results = Vec::new();
    for shape in &request.shapes {
        let data = match prepare_shape(shape, &input, request, None, false) {
            Ok(data) => data,
            Err(error @ RecurrenceError::BudgetExhausted(_)) => {
                let mut result = empty_result(shape_target(shape), Vec::new(), 0);
                result.reason = Some(error.to_string());
                results.push(result);
                continue;
            }
            Err(error) => return Err(error),
        };
        let ShapeData {
            target,
            cols,
            h,
            count,
            rows,
            required,
            automatic,
            propagation,
            rescaling_bound,
        } = data;
        let mut result = empty_result(target, cols.clone(), h);
        result.automatic_degree_bound = automatic;
        result.propagation = propagation;
        result.rescaling_bound = rescaling_bound;
        result.residual_count = count;
        result.last_required_source_row = required
            .checked_sub(1)
            .and_then(|n| input.origin.checked_add(n));
        if count == 0 && result.propagation.is_none() {
            result.status = "insufficient_rows".into();
            result.reason = Some("no complete residual rows in the declared domain".into());
            results.push(result);
            continue;
        }
        result.all_checked_rows_zero = rows.iter().all(Poly::is_zero);
        let matrix = match residual_matrix(
            &rows,
            input.origin,
            request.domain_start,
            count,
            target,
            &cols,
            &request.limits,
        ) {
            Ok(m) => m,
            Err(e) => {
                result.reason = Some(e.to_string());
                results.push(result);
                continue;
            }
        };
        result.equations = matrix.len();
        let kernel = match rational_kernel_with_certificate(
            &matrix,
            cols.len(),
            kernel_limits(&request.limits),
        ) {
            Ok(k) => k,
            Err(RationalKernelError::BudgetExceeded(e)) => {
                result.reason = Some(e.into());
                results.push(result);
                continue;
            }
            Err(_) => {
                return Err(RecurrenceError::InvalidInput(
                    "invalid residual matrix".into(),
                ))
            }
        };
        let leading: Vec<_> = cols
            .iter()
            .enumerate()
            .filter_map(|(i, c)| (c.lag == target).then_some(i))
            .collect();
        let projection: Vec<_> = leading
            .iter()
            .map(|&i| kernel.basis.iter().map(|v| v[i].clone()).collect())
            .collect();
        let projection_kernel = match rational_kernel_with_certificate(
            &projection,
            kernel.basis.len(),
            kernel_limits(&request.limits),
        ) {
            Ok(k) => k,
            Err(_) => {
                result.reason = Some("projection budget".into());
                results.push(result);
                continue;
            }
        };
        result.successor_projection_rank = projection_kernel.pivot_columns.len();
        let absent = result.successor_projection_rank == 0;
        let mut duals = Vec::new();
        if absent {
            for &i in &leading {
                let p = kernel
                    .pivot_columns
                    .iter()
                    .position(|&c| c == i)
                    .ok_or(RecurrenceError::InvalidInput("missing dual pivot".into()))?;
                duals.push(RecurrenceSpaceDual {
                    coordinate: i,
                    weights: kernel.inverse_minor[p]
                        .iter()
                        .map(ToString::to_string)
                        .collect(),
                });
            }
        }
        // An independent set of projected basis columns supplies several natural
        // successor choices; the full kernel remains available for combinations.
        for &i in &projection_kernel.pivot_columns {
            result
                .representatives
                .push(representative(&cols, &kernel.basis[i], target, i, symbol));
        }
        result.status = if absent {
            "absent_in_shape"
        } else {
            "successor_exists"
        }
        .into();
        result.evidence = if absent {
            "rows_only_obstruction"
        } else if input.a.is_some() || input.indexed_a.is_some() {
            "conditional_on_source"
        } else {
            "exact_fit"
        }
        .into();
        result.certificate = Some(RecurrenceSpaceCertificate {
            basis: encode_matrix(&kernel.basis),
            free_columns: kernel.free_columns,
            pivot_columns: kernel.pivot_columns,
            independent_rows: kernel.independent_rows,
            inverse_minor: encode_matrix(&kernel.inverse_minor),
            absence_duals: duals,
        });
        if request.include_leader_coverage {
            result.leader_coverage = leader_coverage(&result, &request.limits).unwrap_or(None);
        }
        results.push(result);
    }
    let frontier = frontier(request, &results);
    Ok(RecurrenceSpaces {
        schema: "polytool.recurrence-spaces.v1".into(),
        request: request.clone(),
        results,
        frontier,
    })
}

fn frontier(request: &RecurrenceSpaceRequest, results: &[RecurrenceSpaceResult]) -> Vec<usize> {
    let key = |i: usize| match request.shapes[i] {
        RecurrenceSpaceShape::Box {
            length,
            derivative_order,
            index_degree,
            variable_degree,
            leading_index_degree,
            leading_variable_degree,
        } => Some((
            length,
            derivative_order,
            (
                index_degree,
                variable_degree,
                leading_index_degree,
                leading_variable_degree,
            ),
        )),
        RecurrenceSpaceShape::AutomaticVariableDegree {
            length,
            derivative_order,
            index_degree,
            leading_index_degree,
            leading_variable_degree,
        } => Some((
            length,
            derivative_order,
            (
                index_degree,
                usize::MAX,
                leading_index_degree,
                leading_variable_degree,
            ),
        )),
        _ => None,
    };
    (0..results.len())
        .filter(|&i| {
            results[i].status == "successor_exists"
                && !(0..results.len()).any(|j| {
                    if results[j].status != "successor_exists" {
                        return false;
                    }
                    match (key(i), key(j)) {
                        (Some((r, s, b)), Some((r2, s2, b2))) => {
                            b == b2 && r2 <= r && s2 <= s && (r2 < r || s2 < s)
                        }
                        _ => false,
                    }
                })
        })
        .collect()
}

/// Replay a report without rerunning nullspace search. The caller supplies
/// the expected request so changed assumptions cannot pass unnoticed.
/// This is exact certificate checking, not a Lean proof or authentication.
pub fn verify_recurrence_spaces(
    request: &RecurrenceSpaceRequest,
    report: &RecurrenceSpaces,
) -> Result<(), RecurrenceError> {
    if report.schema != "polytool.recurrence-spaces.v1"
        || &report.request != request
        || report.results.len() != request.shapes.len()
    {
        return Err(RecurrenceError::VerificationFailed(
            "report/request mismatch".into(),
        ));
    }
    let input = prepare(request)?;
    let symbol = row_symbol(input.rescaling.is_some());
    source_rows(&input, input.rows.len(), &request.limits)?;
    for (shape, result) in request.shapes.iter().zip(&report.results) {
        if result.certificate.is_none() {
            if !["incomplete_budget", "insufficient_rows"].contains(&result.status.as_str())
                || result.evidence != "none"
                || !result.representatives.is_empty()
                || result.leader_coverage.is_some()
                || (result.rescaling_bound.is_some() && request.coefficient_rescaling.is_none())
            {
                return Err(RecurrenceError::VerificationFailed(
                    "uncertified claim".into(),
                ));
            }
            continue;
        }
        let ShapeData {
            target,
            cols,
            h,
            count,
            rows,
            required,
            automatic,
            propagation,
            rescaling_bound,
        } = prepare_shape(shape, &input, request, result.propagation.as_ref(), true)?;
        if target != result.target_lag
            || cols != result.monomials
            || h != result.h
            || automatic != result.automatic_degree_bound
            || propagation != result.propagation
            || rescaling_bound != result.rescaling_bound
        {
            return Err(RecurrenceError::VerificationFailed(
                "support or bound mismatch".into(),
            ));
        }
        if (count == 0 && propagation.is_none())
            || count != result.residual_count
            || result.last_required_source_row
                != required
                    .checked_sub(1)
                    .and_then(|n| input.origin.checked_add(n))
        {
            return Err(RecurrenceError::VerificationFailed(
                "residual range mismatch".into(),
            ));
        }
        if result.all_checked_rows_zero != rows.iter().all(Poly::is_zero) {
            return Err(RecurrenceError::VerificationFailed(
                "degeneracy mismatch".into(),
            ));
        }
        let matrix = residual_matrix(
            &rows,
            input.origin,
            request.domain_start,
            count,
            target,
            &cols,
            &request.limits,
        )?;
        if matrix.len() != result.equations {
            return Err(RecurrenceError::VerificationFailed(
                "equation count mismatch".into(),
            ));
        }
        let cert = result.certificate.as_ref().unwrap();
        let basis = check_kernel_certificate(&matrix, cols.len(), cert, &request.limits)?;
        let rank = cert.pivot_columns.len();
        let dimension = basis.len();
        let leading: Vec<_> = cols
            .iter()
            .enumerate()
            .filter_map(|(i, c)| (c.lag == target).then_some(i))
            .collect();
        let projection: Vec<_> = leading
            .iter()
            .map(|&i| basis.iter().map(|v| v[i].clone()).collect())
            .collect();
        // Only the small target projection is reduced; the full kernel is
        // checked by multiplication and the supplied rank minor.
        let projected = rational_kernel_with_certificate(
            &projection,
            dimension,
            kernel_limits(&request.limits),
        )
        .map_err(|error| match error {
            RationalKernelError::BudgetExceeded(_) => {
                RecurrenceError::BudgetExhausted("projection verification budget".into())
            }
            RationalKernelError::InvalidDimensions => {
                RecurrenceError::VerificationFailed("projection verification dimensions".into())
            }
        })?;
        let projection_rank = projected.pivot_columns.len();
        if projection_rank != result.successor_projection_rank {
            return Err(RecurrenceError::VerificationFailed(
                "successor projection mismatch".into(),
            ));
        }
        let absent = projection_rank == 0;
        let evidence = if absent {
            "rows_only_obstruction"
        } else if input.a.is_some() || input.indexed_a.is_some() {
            "conditional_on_source"
        } else {
            "exact_fit"
        };
        if result.status
            != if absent {
                "absent_in_shape"
            } else {
                "successor_exists"
            }
            || result.evidence != evidence
        {
            return Err(RecurrenceError::VerificationFailed(
                "claim/evidence mismatch".into(),
            ));
        }
        if absent {
            if cert.absence_duals.len() != leading.len() {
                return Err(RecurrenceError::VerificationFailed(
                    "missing absence duals".into(),
                ));
            }
            for (&coordinate, dual) in leading.iter().zip(&cert.absence_duals) {
                let weights: Vec<_> = dual
                    .weights
                    .iter()
                    .map(|s| checked_q(s, &request.limits))
                    .collect::<Result<_, _>>()?;
                if dual.coordinate != coordinate || weights.len() != rank {
                    return Err(RecurrenceError::VerificationFailed(
                        "invalid dual dimensions".into(),
                    ));
                }
                for (c, _) in cols.iter().enumerate() {
                    let sum: Q = (0..rank)
                        .map(|i| &weights[i] * &matrix[cert.independent_rows[i]][c])
                        .sum();
                    if sum != if c == coordinate { Q::one() } else { Q::zero() } {
                        return Err(RecurrenceError::VerificationFailed(
                            "invalid absence dual".into(),
                        ));
                    }
                }
            }
        } else if !cert.absence_duals.is_empty() {
            return Err(RecurrenceError::VerificationFailed(
                "spurious absence duals".into(),
            ));
        }
        let expected: Vec<_> = projected
            .pivot_columns
            .iter()
            .map(|&i| representative(&cols, &basis[i], target, i, symbol))
            .collect();
        if result.representatives != expected {
            return Err(RecurrenceError::VerificationFailed(
                "representative or singular-rank mismatch".into(),
            ));
        }
        if result.leader_coverage.is_some()
            && (!request.include_leader_coverage
                || result.leader_coverage != leader_coverage(result, &request.limits)?)
        {
            return Err(RecurrenceError::VerificationFailed(
                "leader coverage mismatch".into(),
            ));
        }
    }
    if report.frontier != frontier(request, &report.results) {
        return Err(RecurrenceError::VerificationFailed(
            "frontier mismatch".into(),
        ));
    }
    Ok(())
}

mod input;
pub(in crate::recurrence) use input::residual_matrix;
use input::*;
mod analysis;
use analysis::*;
pub(in crate::recurrence) use analysis::{check_kernel_certificate, kernel_limits};
#[cfg(test)]
#[path = "tests.rs"]
mod coverage_tests;
