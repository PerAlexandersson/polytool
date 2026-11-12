//! Exact kernels, representatives and coefficient conditions.

use super::*;

// ---------------------------------------------------------------------------
// Kernels
// ---------------------------------------------------------------------------

fn certificate_of(kernel: &crate::linalg::RationalKernel) -> RecurrenceSpaceCertificate {
    RecurrenceSpaceCertificate {
        basis: encode_matrix(&kernel.basis),
        free_columns: kernel.free_columns.clone(),
        pivot_columns: kernel.pivot_columns.clone(),
        independent_rows: kernel.independent_rows.clone(),
        inverse_minor: encode_matrix(&kernel.inverse_minor),
        absence_duals: Vec::new(),
    }
}

/// The homogeneous residual matrix (every column with the same sign).
pub(super) fn homogeneous_matrix(
    rows: &[Poly],
    origin: usize,
    start: usize,
    count: usize,
    cols: &[RecurrenceMonomial],
    limits: &RecurrenceSpaceLimits,
) -> Result<Vec<Vec<Q>>, RecurrenceError> {
    // No column has lag usize::MAX, so every column enters negatively and the
    // kernel is that of the homogeneous ansatz.
    residual_matrix(rows, origin, start, count, usize::MAX, cols, limits)
}

/// The condition matrix: one row per base, one column per monomial.
pub(super) fn condition_matrix(
    rows: &[Vec<Q>],
    origin: usize,
    start: usize,
    count: usize,
    cols: &[ConditionMonomial],
    limits: &RecurrenceSpaceLimits,
) -> Result<Vec<Vec<Q>>, RecurrenceError> {
    if count
        .checked_mul(cols.len())
        .is_none_or(|e| e > limits.max_matrix_entries)
    {
        return Err(RecurrenceError::BudgetExhausted(
            "condition matrix entry budget".into(),
        ));
    }
    (start..start + count)
        .map(|n| {
            cols.iter()
                .map(|c| {
                    let row =
                        rows.get(n + c.shift - origin)
                            .ok_or(RecurrenceError::InvalidInput(
                                "missing row for a condition".into(),
                            ))?;
                    let value = row
                        .get(c.coefficient_index)
                        .cloned()
                        .unwrap_or_else(Q::zero);
                    let power = u32::try_from(c.index_degree).map_err(|_| {
                        RecurrenceError::LimitExceeded("index degree overflow".into())
                    })?;
                    Ok(value * Q::from_integer(BigInt::from(n).pow(power)))
                })
                .collect()
        })
        .collect()
}

pub(super) fn relation_from_vector(cols: &[RecurrenceMonomial], vector: &[Q]) -> ImplicitRelation {
    let mut op = DOp::default();
    for (c, value) in cols.iter().zip(vector) {
        if !value.is_zero() {
            op.add_term(
                (c.lag, c.derivative_order),
                Biv::monomial(c.index_degree, c.variable_degree, value.clone()),
            );
        }
    }
    ImplicitRelation::from_op(&op)
}

pub(super) fn condition_from_vector(
    cols: &[ConditionMonomial],
    vector: &[Q],
) -> CoefficientCondition {
    let mut terms = Vec::new();
    for (c, value) in cols.iter().zip(vector) {
        if value.is_zero() {
            continue;
        }
        let mut coeff = vec![Q::zero(); c.index_degree + 1];
        coeff[c.index_degree] = value.clone();
        terms.push(ConditionTerm {
            shift: c.shift,
            coefficient_index: c.coefficient_index,
            coeff,
        });
    }
    CoefficientCondition { terms }.canonical()
}

pub(super) fn leader_columns(cols: &[RecurrenceMonomial], r: usize) -> Vec<usize> {
    cols.iter()
        .enumerate()
        .filter_map(|(i, c)| (r == 0 || c.lag == r).then_some(i))
        .collect()
}

/// Rank of the kernel basis projected on `leading`, and the basis indices of
/// an independent set of projections.
pub(super) fn projection(
    basis: &[Vec<Q>],
    leading: &[usize],
    limits: &RecurrenceSpaceLimits,
) -> Result<Vec<usize>, RecurrenceError> {
    let matrix = leading
        .iter()
        .map(|&i| basis.iter().map(|v| v[i].clone()).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let kernel = rational_kernel_with_certificate(&matrix, basis.len(), kernel_limits(limits))
        .map_err(|_| RecurrenceError::BudgetExhausted("projection budget".to_string()))?;
    Ok(kernel.pivot_columns)
}

// ---------------------------------------------------------------------------
// The search
// ---------------------------------------------------------------------------

fn empty_result(support: &ImplicitSearchSupport) -> ImplicitSupportResult {
    ImplicitSupportResult {
        support: support.clone(),
        status: "incomplete_budget".into(),
        evidence: "none".into(),
        reason: None,
        monomials: Vec::new(),
        annihilator_power: None,
        residual_count: 0,
        last_required_row: None,
        equations: 0,
        kernel_dimension: 0,
        leader_projection_rank: 0,
        certificate: None,
        completeness: None,
        representatives: Vec::new(),
    }
}

pub(super) fn annihilator_power(theorem: Theorem, cols: &[RecurrenceMonomial]) -> Option<usize> {
    let us = cols
        .iter()
        .map(|c| c.index_degree.checked_add(c.derivative_order))
        .collect::<Option<Vec<_>>>()?
        .into_iter()
        .max()
        .unwrap_or(0);
    match theorem {
        Theorem::Identity | Theorem::Transported { .. } => us.checked_add(1),
        Theorem::Derivative(m) => us.checked_add(m)?.checked_add(1),
        Theorem::Factorial(m) => {
            RecurrenceCoefficientRescaling::FactorialPower { power: m }.annihilator_power(cols)
        }
        Theorem::Weight { m: 0, total, .. } => us.checked_add(total)?.checked_add(1),
        Theorem::Weight { m, total, .. } => {
            RecurrenceCoefficientRescaling::FactorialPower { power: m }
                .annihilator_power(cols)?
                .checked_add(total)
        }
    }
}

/// `(K', e)` of the condition theorem: the truncation index `K'` (the
/// conditions on `[x^j] Q`, `j <= k`, depend on `P_n mod x^(K'+1)`) and the
/// extra index degree `e` contributed by the map; `None` when no condition
/// theorem applies.
pub(super) fn condition_theorem(theorem: Theorem, k: usize) -> Option<(usize, usize)> {
    match theorem {
        Theorem::Identity | Theorem::Factorial(_) => Some((k, 0)),
        Theorem::Derivative(m) => Some((k.checked_add(m)?, 0)),
        Theorem::Transported { regular_at_zero } => regular_at_zero.then_some((k, 0)),
        Theorem::Weight { index, .. } => Some((k, index)),
    }
}

pub(super) fn run_support(
    context: &mut Context<'_>,
    support: &ImplicitSearchSupport,
) -> Result<ImplicitSupportResult, RecurrenceError> {
    let limits = &context.request.limits;
    let mut result = empty_result(support);
    let cols = match support.monomials(limits) {
        Ok(cols) => cols,
        Err(error @ RecurrenceError::BudgetExhausted(_)) => {
            result.reason = Some(error.to_string());
            return Ok(result);
        }
        Err(error) => return Err(error),
    };
    result.monomials = cols.clone();
    let r = support.max_shift;
    let offset = context.start - context.origin;
    let (count, complete) = match (context.theorem, context.order) {
        (Some(theorem), Some(k)) => {
            let Some(h) = annihilator_power(theorem, &cols) else {
                result.reason = Some("annihilator power overflow".into());
                return Ok(result);
            };
            result.annihilator_power = Some(h);
            match h.checked_mul(k) {
                Some(count) => (count, true),
                None => {
                    result.reason = Some("residual bound overflow".into());
                    return Ok(result);
                }
            }
        }
        _ => (context.rows.len().saturating_sub(offset + r), false),
    };
    let needed = offset + count + r;
    if complete {
        if let Err(error) = context.ensure(needed) {
            if matches!(error, RecurrenceError::BudgetExhausted(_)) {
                result.reason = Some(error.to_string());
                return Ok(result);
            }
            return Err(error);
        }
    }
    result.residual_count = count;
    if count == 0 || context.rows.len() < needed {
        result.status = "insufficient_rows".into();
        result.reason = Some("no complete residual in the declared domain".into());
        return Ok(result);
    }
    result.last_required_row = Some(context.origin + needed - 1);
    let polys = context.polys();
    let matrix =
        match homogeneous_matrix(&polys, context.origin, context.start, count, &cols, limits) {
            Ok(m) => m,
            Err(error) => {
                result.reason = Some(error.to_string());
                return Ok(result);
            }
        };
    result.equations = matrix.len();
    let kernel = match rational_kernel_with_certificate(&matrix, cols.len(), kernel_limits(limits))
    {
        Ok(k) => k,
        Err(RationalKernelError::BudgetExceeded(e)) => {
            result.reason = Some(e.into());
            return Ok(result);
        }
        Err(_) => {
            return Err(RecurrenceError::InvalidInput(
                "invalid residual matrix".into(),
            ))
        }
    };
    let leading = leader_columns(&cols, r);
    let pivots = match projection(&kernel.basis, &leading, limits) {
        Ok(p) => p,
        Err(error) => {
            result.reason = Some(error.to_string());
            return Ok(result);
        }
    };
    result.kernel_dimension = kernel.basis.len();
    result.leader_projection_rank = pivots.len();
    result.certificate = Some(certificate_of(&kernel));
    if pivots.is_empty() {
        result.status = "absent_in_support".into();
        result.evidence = "rows_only_obstruction".into();
        return Ok(result);
    }
    result.status = "relations_found".into();
    result.evidence = if complete {
        CERTIFIED_COMPLETE_SPACE_EVIDENCE
    } else {
        FINITE_FIT_EVIDENCE
    }
    .into();
    if complete {
        result.completeness = Some(format!(
            "{}. Every relation in this support that vanishes at the {count} bases n = {}..{} \
             vanishes for all n >= {}; the kernel is the complete space of all-index relations \
             in the support.",
            context_basis(context),
            context.start,
            context.start + count - 1,
            context.start
        ));
    }
    for &index in pivots.iter().take(MAX_IMPLICIT_SEARCH_REPRESENTATIVES) {
        let representative =
            representative(context, &cols, &kernel.basis[index], index, r, complete)?;
        result.representatives.push(representative);
    }
    Ok(result)
}

fn context_basis(context: &Context<'_>) -> String {
    match (context.theorem, context.order) {
        (Some(Theorem::Identity), Some(k)) => {
            format!("P satisfies a constant-index monic pure-lag rule of order {k}")
        }
        (Some(Theorem::Factorial(m)), Some(k)) => format!(
            "Q = T(P), T(x^j) = x^j/(j!)^{m}, with P satisfying a constant-index monic pure-lag \
             rule of order {k}"
        ),
        (Some(Theorem::Derivative(m)), Some(k)) => format!(
            "Q = D^{m} P with P satisfying a constant-index monic pure-lag rule of order {k}"
        ),
        (Some(Theorem::Transported { .. }), Some(k)) => format!(
            "Q = T(P) for a Möbius-class map T, so Q satisfies the transported monic rule of \
             order {k} with index-free coefficients in Q(x)"
        ),
        (Some(Theorem::Weight { m, total, .. }), Some(k)) => format!(
            "Q_n = lambda(n, theta) P_n with lambda = (k!)^(-{m}) times a polynomial of total \
             degree {total}, P satisfying a constant-index monic pure-lag rule of order {k}"
        ),
        _ => String::new(),
    }
}

fn representative(
    context: &mut Context<'_>,
    cols: &[RecurrenceMonomial],
    vector: &[Q],
    basis_index: usize,
    r: usize,
    complete: bool,
) -> Result<ImplicitSearchRepresentative, RecurrenceError> {
    let relation = relation_from_vector(cols, vector).canonical();
    let symbol = context.symbol;
    let mut first_base = context.start;
    let mut evidence = if complete {
        CERTIFIED_COMPLETE_SPACE_EVIDENCE
    } else {
        FINITE_FIT_EVIDENCE
    };
    let mut division = None;
    if !complete {
        if let (Some((l, l_first)), Some(t)) = (&context.defining, &context.forward) {
            if let Some(certificate) = find_reconstruction_certificate(&relation, t, l) {
                certificate.verify_identity()?;
                evidence = CERTIFIED_BY_DIVISION_EVIDENCE;
                first_base = first_base.max(*l_first);
                division = Some(certificate);
            }
        }
    }
    // Sanity: the relation holds on every available row of its domain.
    let end = (context.origin + context.rows.len()).saturating_sub(relation.max_shift());
    for base in first_base..end {
        if relation.holds_at(&context.rows, context.origin, base) == Some(false) {
            if complete || division.is_some() {
                return Err(RecurrenceError::InvalidInput(format!(
                    "internal: a certified relation fails at base {base}"
                )));
            }
            return Err(RecurrenceError::InvalidInput(format!(
                "internal: a fitted relation fails at base {base}"
            )));
        }
    }
    let note = match evidence {
        CERTIFIED_COMPLETE_SPACE_EVIDENCE => "member of a certified complete space",
        CERTIFIED_BY_DIVISION_EVIDENCE => {
            "fitted, then proved by right division by the defining rule"
        }
        _ => "fitted to the rows (finite evidence)",
    };
    let described = DescribedRelation::new(&relation, first_base, evidence, symbol).with_note(note);
    let seed_end = (first_base + r).max(context.start + 1);
    if context.rows.len() < seed_end - context.origin {
        return Err(RecurrenceError::InvalidInput(
            "internal: missing seed rows".into(),
        ));
    }
    let seeds = context.rows[context.start - context.origin..seed_end - context.origin]
        .iter()
        .map(|row| encode_poly(row))
        .collect::<Vec<_>>();
    let mut description = GeneratingDescription {
        first_index: context.start,
        seeds,
        relations: vec![described.clone()],
        exceptional_index_conditions: Vec::new(),
        boundary_conditions: Vec::new(),
    };
    let mut assessment = assess_generating_description(&description)?;
    let mut searches = Vec::new();
    if context.request.conditions && r > 0 && assessment.uniqueness.status != "certified" {
        if let Some(indices) = valuation_exceptional_indices(&relation, first_base) {
            if !indices.is_empty()
                && indices
                    .iter()
                    .all(|&k| k <= MAX_CONDITION_COEFFICIENT_INDEX)
            {
                for k in indices {
                    let search = condition_search(context, k, r, first_base)?;
                    if let Some(condition) = &search.condition {
                        description.boundary_conditions.push(condition.clone());
                    }
                    searches.push(search);
                }
                assessment = assess_generating_description(&description)?;
            }
        }
    }
    Ok(ImplicitSearchRepresentative {
        basis_index,
        relation: described,
        division_certificate: division,
        condition_searches: searches,
        description,
        assessment,
    })
}

fn condition_columns(r: usize, k: usize, degree: usize) -> Vec<ConditionMonomial> {
    let mut cols = Vec::new();
    for shift in 0..=r {
        for coefficient_index in 0..=k {
            for index_degree in 0..=degree {
                cols.push(ConditionMonomial {
                    shift,
                    coefficient_index,
                    index_degree,
                });
            }
        }
    }
    cols.sort();
    cols
}

fn condition_search(
    context: &mut Context<'_>,
    k: usize,
    r: usize,
    first_base: usize,
) -> Result<ImplicitConditionSearch, RecurrenceError> {
    let degree = context.request.condition_index_degree;
    let cols = condition_columns(r, k, degree);
    let mut search = ImplicitConditionSearch {
        coefficient_index: k,
        status: "incomplete_budget".into(),
        evidence: "none".into(),
        reason: None,
        monomials: cols.clone(),
        residual_count: 0,
        certificate: None,
        chosen_basis_index: None,
        condition: None,
    };
    let limits = context.request.limits.clone();
    if cols.len() > limits.max_columns || degree > limits.max_degree {
        search.reason = Some("condition column budget".into());
        return Ok(search);
    }
    let offset = context.start - context.origin;
    let condition = context.theorem.and_then(|t| condition_theorem(t, k));
    let (count, complete) = match (condition, context.order) {
        (Some((truncation, extra)), Some(order)) => match truncation
            .checked_add(1)
            .and_then(|t| order.checked_mul(t))
            .and_then(|c| c.checked_mul(degree.checked_add(extra)?.checked_add(1)?))
        {
            Some(count) => (count, true),
            None => {
                search.reason = Some("condition residual bound overflow".into());
                return Ok(search);
            }
        },
        _ => (context.rows.len().saturating_sub(offset + r), false),
    };
    let needed = offset + count + r;
    if complete {
        if let Err(error) = context.ensure(needed) {
            if matches!(error, RecurrenceError::BudgetExhausted(_)) {
                search.reason = Some(error.to_string());
                return Ok(search);
            }
            return Err(error);
        }
    }
    search.residual_count = count;
    if count == 0 || context.rows.len() < needed {
        search.status = "insufficient_rows".into();
        search.reason = Some("no complete condition residual in the declared domain".into());
        return Ok(search);
    }
    let matrix = match condition_matrix(
        &context.rows,
        context.origin,
        context.start,
        count,
        &cols,
        &limits,
    ) {
        Ok(m) => m,
        Err(error) => {
            search.reason = Some(error.to_string());
            return Ok(search);
        }
    };
    let kernel = match rational_kernel_with_certificate(&matrix, cols.len(), kernel_limits(&limits))
    {
        Ok(kernel) => kernel,
        Err(RationalKernelError::BudgetExceeded(e)) => {
            search.reason = Some(e.into());
            return Ok(search);
        }
        Err(_) => {
            return Err(RecurrenceError::InvalidInput(
                "invalid condition matrix".into(),
            ))
        }
    };
    search.certificate = Some(certificate_of(&kernel));
    if kernel.basis.is_empty() {
        search.status = "absent_in_support".into();
        search.evidence = "rows_only_obstruction".into();
        return Ok(search);
    }
    let evidence = if complete {
        CERTIFIED_COMPLETE_SPACE_EVIDENCE
    } else {
        FINITE_FIT_EVIDENCE
    };
    // A basis vector whose free column is [x^k] of the row n + r has its top
    // shift terms at indices <= k; prefer the smallest power of n.
    for (i, &free) in kernel.free_columns.iter().enumerate() {
        let column = cols[free];
        if column.shift != r || column.coefficient_index != k {
            continue;
        }
        let condition = condition_from_vector(&cols, &kernel.basis[i]);
        let pins = condition
            .terms
            .iter()
            .find(|t| t.shift == r && t.coefficient_index == k)
            .is_some_and(|t| {
                super::super::algebra::roots_in_domain(&t.coeff, first_base, None).is_empty()
            });
        if !pins {
            continue;
        }
        let note = if complete {
            "member of a certified complete condition space (truncated source recurrence)"
        } else {
            "fitted to the rows (finite evidence)"
        };
        search.status = "condition_found".into();
        search.evidence = evidence.into();
        search.chosen_basis_index = Some(i);
        search.condition = Some(
            DescribedCondition::new(&condition, context.start, evidence, Some(k), context.symbol)
                .with_note(note),
        );
        return Ok(search);
    }
    search.status = "no_covering_condition".into();
    search.evidence = evidence.into();
    search.reason = Some(format!(
        "no condition in the support pins [x^{k}] of the row n + {r} with a coefficient that \
         is nonzero on the domain"
    ));
    Ok(search)
}
