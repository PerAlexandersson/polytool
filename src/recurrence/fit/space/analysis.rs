//! Recurrence-space analysis, representatives and kernel certificates.

use super::*;

fn dot(a: &[Q], b: &[Q]) -> Q {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

pub(super) fn singular_ranks(
    cols: &[RecurrenceMonomial],
    vector: &[Q],
    target: usize,
) -> (Option<Vec<String>>, bool) {
    let mut by_x: BTreeMap<usize, Vec<Q>> = BTreeMap::new();
    for (c, v) in cols
        .iter()
        .zip(vector)
        .filter(|(c, v)| c.lag == target && !v.is_zero())
    {
        let a = by_x.entry(c.variable_degree).or_default();
        a.resize(a.len().max(c.index_degree + 1), Q::zero());
        a[c.index_degree] += v;
    }
    let requires = by_x.keys().any(|&x| x > 0);
    (
        integer_zeros_of_gcd(
            &by_x
                .values()
                .map(|p| Poly::new(p.clone()))
                .collect::<Vec<_>>(),
        ),
        requires,
    )
}

pub(super) fn integer_zeros_of_gcd(polynomials: &[Poly]) -> Option<Vec<String>> {
    if polynomials.iter().any(|p| p.degree() == Some(0)) {
        return Some(vec![]);
    }
    if polynomials.iter().any(|p| {
        p.degree().unwrap_or(0) > 8
            || p.coeffs()
                .iter()
                .any(|c| c.numer().bits().max(c.denom().bits()) > 128)
    }) {
        return None;
    }
    let mut gcd = Poly::zero();
    for p in polynomials {
        gcd = gcd.gcd(p);
    }
    if gcd.is_zero()
        || gcd.degree().unwrap_or(0) > 8
        || gcd
            .coeffs()
            .iter()
            .any(|c| c.numer().bits().max(c.denom().bits()) > 128)
    {
        return None;
    }
    let integer = crate::vec_poly::clear_denominators(gcd.coeffs());
    let chain = SturmChain::from_bigint_coeffs(&integer);
    let mut roots = BTreeSet::new();
    for (lo, hi) in chain.isolate_roots(&Q::new(BigInt::one(), BigInt::from(2))) {
        let first = lo.numer().div_ceil(lo.denom());
        let last = hi.numer().div_floor(hi.denom());
        if first <= last
            && !first.is_negative()
            && gcd.evaluate(&Q::from_integer(first.clone())).is_zero()
        {
            roots.insert(first);
        }
    }
    Some(roots.into_iter().map(|n| n.to_string()).collect())
}

pub(super) fn grouped_coefficients(
    cols: &[RecurrenceMonomial],
    vector: &[Q],
) -> BTreeMap<(usize, usize), BivarPoly> {
    let mut groups: BTreeMap<(usize, usize), Vec<Vec<Q>>> = BTreeMap::new();
    for (c, v) in cols.iter().zip(vector).filter(|(_, v)| !v.is_zero()) {
        let p = groups.entry((c.lag, c.derivative_order)).or_default();
        p.resize(p.len().max(c.index_degree + 1), Vec::new());
        let width = p[c.index_degree].len().max(c.variable_degree + 1);
        p[c.index_degree].resize(width, Q::zero());
        p[c.index_degree][c.variable_degree] += v;
    }
    groups
        .into_iter()
        .map(|(k, p)| (k, BivarPoly::from_coefficients_normalized(p)))
        .collect()
}

pub(super) fn representative(
    cols: &[RecurrenceMonomial],
    vector: &[Q],
    target: usize,
    index: usize,
    symbol: &str,
) -> RecurrenceSpaceRepresentative {
    let (singular_ranks, requires_x_division) = singular_ranks(cols, vector, target);
    let groups = grouped_coefficients(cols, vector);
    let side = |leading: bool| {
        groups
            .iter()
            .filter(|((lag, _), _)| (*lag == target) == leading)
            .map(|((lag, d), p)| {
                let row = if *d == 0 {
                    format!("{symbol}[n+{lag}]")
                } else {
                    format!("D^{d}({symbol}[n+{lag}])")
                };
                format!("({}) {row}", p.to_string().replace('t', "X"))
            })
            .collect::<Vec<_>>()
            .join(" + ")
    };
    let rhs = side(false);
    RecurrenceSpaceRepresentative {
        basis_index: index,
        coefficients: vector.iter().map(ToString::to_string).collect(),
        formula: format!(
            "{} = {}",
            side(true),
            if rhs.is_empty() { "0" } else { &rhs }
        ),
        singular_ranks,
        requires_x_division,
    }
}

pub(super) fn empty_result(
    target: usize,
    cols: Vec<RecurrenceMonomial>,
    h: usize,
) -> RecurrenceSpaceResult {
    RecurrenceSpaceResult {
        status: "incomplete_budget".into(),
        evidence: "none".into(),
        reason: None,
        target_lag: target,
        monomials: cols,
        h,
        residual_count: 0,
        last_required_source_row: None,
        equations: 0,
        successor_projection_rank: 0,
        certificate: None,
        representatives: Vec::new(),
        all_checked_rows_zero: false,
        automatic_degree_bound: None,
        propagation: None,
        leader_coverage: None,
        rescaling_bound: None,
    }
}

/// Row name used in formulas: `Q` for rescaled rows, `P` for source rows.
pub(super) fn row_symbol(rescaled: bool) -> &'static str {
    if rescaled {
        "Q"
    } else {
        "P"
    }
}

pub(super) fn shape_target(shape: &RecurrenceSpaceShape) -> usize {
    match shape {
        RecurrenceSpaceShape::Box { length, .. }
        | RecurrenceSpaceShape::AutomaticVariableDegree { length, .. } => *length,
        RecurrenceSpaceShape::Explicit { target_lag, .. } => *target_lag,
    }
}
pub(super) struct ShapeData {
    pub(super) cols: Vec<RecurrenceMonomial>,
    pub(super) target: usize,
    pub(super) h: usize,
    pub(super) count: usize,
    pub(super) rows: Vec<Poly>,
    pub(super) required: usize,
    pub(super) automatic: Option<RecurrenceAutomaticDegreeBound>,
    pub(super) propagation: Option<RecurrencePropagationCertificate>,
    pub(super) rescaling_bound: Option<RecurrenceRescalingBound>,
}
pub(super) fn prepare_shape(
    shape: &RecurrenceSpaceShape,
    input: &Input,
    request: &RecurrenceSpaceRequest,
    supplied: Option<&RecurrencePropagationCertificate>,
    verifying: bool,
) -> Result<ShapeData, RecurrenceError> {
    let lim = &request.limits;
    let mut resolved = match shape {
        RecurrenceSpaceShape::AutomaticVariableDegree {
            length,
            derivative_order,
            index_degree,
            leading_index_degree,
            leading_variable_degree,
        } => {
            if input.a.is_none() && input.indexed_a.is_none() {
                return Err(RecurrenceError::InvalidInput(
                    "automatic degree requires a defining source".into(),
                ));
            }
            if input.rescaling.is_some() {
                return Err(RecurrenceError::Unsupported(
                    "automatic_variable_degree is not supported with coefficient_rescaling".into(),
                ));
            }
            RecurrenceSpaceShape::Box {
                length: *length,
                derivative_order: *derivative_order,
                index_degree: *index_degree,
                variable_degree: 0,
                leading_index_degree: *leading_index_degree,
                leading_variable_degree: *leading_variable_degree,
            }
        }
        other => other.clone(),
    };
    let (target, mut cols) = shape_columns(&resolved, lim)?;
    let h = if let Some(rescaling) = &input.rescaling {
        rescaling
            .annihilator_power(&cols)
            .ok_or(RecurrenceError::LimitExceeded("rescaled h overflow".into()))?
    } else {
        cols.iter()
            .map(|c| c.index_degree.checked_add(c.derivative_order))
            .collect::<Option<Vec<_>>>()
            .ok_or(RecurrenceError::LimitExceeded(
                "support degree overflow".into(),
            ))?
            .into_iter()
            .max()
            .unwrap_or(0)
            .checked_add(1)
            .ok_or(RecurrenceError::LimitExceeded("h overflow".into()))?
    };
    let propagation = if let Some(a) = &input.indexed_a {
        if verifying {
            let cert = supplied.ok_or(RecurrenceError::InvalidInput(
                "missing indexed-source propagation certificate".into(),
            ))?;
            propagation::verify(a, &cols, request.domain_start, cert, lim)?;
            Some(cert.clone())
        } else {
            Some(propagation::build(a, &cols, request.domain_start, lim)?)
        }
    } else {
        if supplied.is_some() {
            return Err(RecurrenceError::InvalidInput(
                "unexpected propagation certificate".into(),
            ));
        }
        None
    };
    let offset = request.domain_start - input.origin;
    let count = if let Some(cert) = &propagation {
        cert.regular_from
            .checked_sub(request.domain_start)
            .and_then(|n| n.checked_add(cert.coefficients.len() - 1))
            .ok_or(RecurrenceError::LimitExceeded(
                "propagation bound overflow".into(),
            ))?
    } else if let Some(a) = &input.a {
        h.checked_mul(a.len())
            .ok_or(RecurrenceError::LimitExceeded(
                "residual bound overflow".into(),
            ))?
    } else {
        input
            .rows
            .len()
            .saturating_sub(offset.saturating_add(target))
    };
    let required = offset
        .checked_add(count)
        .and_then(|n| n.checked_add(target))
        .ok_or(RecurrenceError::LimitExceeded("row bound overflow".into()))?;
    let mut rows = source_rows(input, required, lim)?;
    let mut rescaling_bound = None;
    if let Some(rescaling) = input.rescaling {
        // Rescaling reaches this point only with a constant-index source.
        let order = input.a.as_ref().map_or(0, Vec::len);
        rows = rows
            .iter()
            .map(|p| {
                rescaling
                    .rescale_coefficients(p.coeffs(), lim)
                    .map(Poly::new)
            })
            .collect::<Result<_, _>>()?;
        let initial = rows
            .get(
                offset
                    ..offset
                        .checked_add(target)
                        .ok_or(RecurrenceError::LimitExceeded("row bound overflow".into()))?,
            )
            .ok_or(RecurrenceError::InvalidInput(
                "missing rescaled initial rows".into(),
            ))?;
        rescaling_bound = Some(RecurrenceRescalingBound::new(
            rescaling,
            order,
            &cols,
            encode_polys(initial),
            request.domain_start,
        )?);
    }
    let mut automatic = None;
    if matches!(shape, RecurrenceSpaceShape::AutomaticVariableDegree { .. }) {
        let mut rhs_degree = 0;
        let mut leader_degree = 0;
        let mut rhs_nonzero = false;
        let mut work = lim.max_operations;
        for n in request.domain_start
            ..request
                .domain_start
                .checked_add(count)
                .ok_or(RecurrenceError::LimitExceeded("rank overflow".into()))?
        {
            for c in &cols {
                let row = rows
                    .get(n - input.origin + c.lag)
                    .ok_or(RecurrenceError::InvalidInput("missing bound row".into()))?;
                let cost = row
                    .coeffs()
                    .len()
                    .checked_mul(c.derivative_order.saturating_add(2))
                    .ok_or(RecurrenceError::BudgetExhausted(
                        "bound arithmetic budget".into(),
                    ))?;
                work = work
                    .checked_sub(cost)
                    .ok_or(RecurrenceError::BudgetExhausted(
                        "bound arithmetic budget".into(),
                    ))?;
                let p = column_at(c, &rows, n, input.origin, lim)?;
                if c.lag == target {
                    leader_degree = leader_degree.max(p.degree().unwrap_or(0));
                } else {
                    rhs_degree = rhs_degree.max(p.degree().unwrap_or(0));
                    rhs_nonzero |= !p.is_zero();
                }
            }
        }
        let rank_upper_bound = if rhs_nonzero {
            count.min(cols.iter().filter(|c| c.lag < target).count())
        } else {
            0
        };
        let variable_degree = if rank_upper_bound == 0 {
            0
        } else {
            let first = leader_degree
                .checked_add((rank_upper_bound - 1).checked_mul(rhs_degree).ok_or(
                    RecurrenceError::BudgetExhausted("automatic degree budget".into()),
                )?)
                .ok_or(RecurrenceError::BudgetExhausted(
                    "automatic degree budget".into(),
                ))?;
            let second = rank_upper_bound
                .checked_mul(rhs_degree)
                .ok_or(RecurrenceError::BudgetExhausted(
                    "automatic degree budget".into(),
                ))?
                .saturating_sub(1);
            first.max(second)
        };
        automatic = Some(RecurrenceAutomaticDegreeBound {
            rhs_degree,
            leader_degree,
            rank_upper_bound,
            variable_degree,
        });
        if let RecurrenceSpaceShape::Box {
            variable_degree: v, ..
        } = &mut resolved
        {
            *v = variable_degree;
        }
        cols = shape_columns(&resolved, lim)
            .map_err(|error| {
                error.with_context(format!("automatic RHS X-degree {variable_degree}"))
            })?
            .1;
    }
    Ok(ShapeData {
        cols,
        target,
        h,
        count,
        rows,
        required,
        automatic,
        propagation,
        rescaling_bound,
    })
}

pub(super) fn leader_coverage(
    result: &RecurrenceSpaceResult,
    limits: &RecurrenceSpaceLimits,
) -> Result<Option<RecurrenceLeaderCoverage>, RecurrenceError> {
    if result.representatives.is_empty() {
        return Ok(None);
    }
    let mut combined: BTreeMap<RecurrenceMonomial, Q> = BTreeMap::new();
    let mut multipliers = Vec::new();
    let mut leaders = Vec::new();
    let mut remaining = limits.max_operations;
    for rep in &result.representatives {
        let v = rep
            .coefficients
            .iter()
            .map(|s| checked_q(s, limits))
            .collect::<Result<Vec<_>, _>>()?;
        let group = grouped_coefficients(&result.monomials, &v);
        let leader = &group[&(result.target_lag, 0)];
        multipliers.push(encode_matrix(&leader.coeffs));
        let width = leader.coeffs.iter().map(Vec::len).max().unwrap_or(0);
        for x in 0..width {
            leaders.push(Poly::new(
                leader
                    .coeffs
                    .iter()
                    .map(|r| r.get(x).cloned().unwrap_or_else(Q::zero))
                    .collect(),
            ));
        }
        for (a, ac) in result
            .monomials
            .iter()
            .zip(&v)
            .filter(|(a, c)| a.lag == result.target_lag && !c.is_zero())
        {
            for (b, bc) in result
                .monomials
                .iter()
                .zip(&v)
                .filter(|(_, c)| !c.is_zero())
            {
                remaining = remaining
                    .checked_sub(1)
                    .ok_or(RecurrenceError::BudgetExhausted(
                        "leader combination arithmetic budget".into(),
                    ))?;
                let n = a.index_degree.checked_add(b.index_degree).ok_or(
                    RecurrenceError::BudgetExhausted("leader combination degree budget".into()),
                )?;
                let x = a.variable_degree.checked_add(b.variable_degree).ok_or(
                    RecurrenceError::BudgetExhausted("leader combination degree budget".into()),
                )?;
                if n > limits.max_degree || x > limits.max_degree {
                    return Err(RecurrenceError::BudgetExhausted(
                        "leader combination degree budget".into(),
                    ));
                }
                let key = RecurrenceMonomial {
                    index_degree: n,
                    variable_degree: x,
                    ..*b
                };
                let c = combined.entry(key).or_insert_with(Q::zero);
                *c += ac * bc;
                if c.numer().bits().max(c.denom().bits()) > limits.max_coefficient_bits {
                    return Err(RecurrenceError::BudgetExhausted(
                        "leader combination coefficient budget".into(),
                    ));
                }
                if combined.len() > limits.max_matrix_entries {
                    return Err(RecurrenceError::BudgetExhausted(
                        "leader combination entry budget".into(),
                    ));
                }
            }
        }
    }
    combined.retain(|_, v| !v.is_zero());
    let (cols, vector): (Vec<_>, Vec<_>) = combined.into_iter().unzip();
    let rep = representative(
        &cols,
        &vector,
        result.target_lag,
        0,
        row_symbol(result.rescaling_bound.is_some()),
    );
    Ok(Some(RecurrenceLeaderCoverage {
        basis_indices: result
            .representatives
            .iter()
            .map(|p| p.basis_index)
            .collect(),
        multipliers,
        monomials: cols,
        coefficients: rep.coefficients,
        formula: rep.formula,
        common_singular_ranks: integer_zeros_of_gcd(&leaders),
        requires_x_division: rep.requires_x_division,
    }))
}

pub(in crate::recurrence) fn kernel_limits(l: &RecurrenceSpaceLimits) -> RationalKernelLimits {
    RationalKernelLimits {
        max_entries: l.max_matrix_entries,
        max_operations: l.max_operations,
        max_coefficient_bits: l.max_coefficient_bits,
    }
}

/// Check an exact kernel certificate against `matrix` (`width` columns):
/// dimensions, every basis vector is a null vector with the unit
/// free-column witness, and `inverse_minor` inverts the rank minor. Returns
/// the parsed basis. Shared by the complete-space and implicit-relation
/// verifiers.
pub(in crate::recurrence) fn check_kernel_certificate(
    matrix: &[Vec<Q>],
    width: usize,
    cert: &RecurrenceSpaceCertificate,
    limits: &RecurrenceSpaceLimits,
) -> Result<Vec<Vec<Q>>, RecurrenceError> {
    let basis = parse_matrix(&cert.basis, limits)?;
    let inverse = parse_matrix(&cert.inverse_minor, limits)?;
    let rank = cert.pivot_columns.len();
    let dimension = basis.len();
    let verify_work = matrix
        .len()
        .checked_mul(width)
        .and_then(|x| x.checked_mul(dimension + rank + 1))
        .and_then(|x| rank.checked_pow(3).and_then(|y| x.checked_add(y)))
        .ok_or(RecurrenceError::LimitExceeded(
            "verification work overflow".into(),
        ))?;
    if verify_work > limits.max_operations {
        return Err(RecurrenceError::BudgetExhausted(
            "certificate verification arithmetic budget".into(),
        ));
    }
    let partition: BTreeSet<_> = cert
        .pivot_columns
        .iter()
        .chain(&cert.free_columns)
        .copied()
        .collect();
    if rank + dimension != width
        || cert.free_columns.len() != dimension
        || partition != (0..width).collect()
        || cert.independent_rows.len() != rank
        || inverse.len() != rank
        || inverse.iter().any(|r| r.len() != rank)
        || cert.independent_rows.iter().any(|&i| i >= matrix.len())
        || basis.iter().any(|v| v.len() != width)
    {
        return Err(RecurrenceError::VerificationFailed(
            "invalid kernel dimensions".into(),
        ));
    }
    for (i, v) in basis.iter().enumerate() {
        if matrix.iter().any(|r| !dot(r, v).is_zero()) {
            return Err(RecurrenceError::VerificationFailed(
                "invalid null vector".into(),
            ));
        }
        for (j, &c) in cert.free_columns.iter().enumerate() {
            if v[c] != if i == j { Q::one() } else { Q::zero() } {
                return Err(RecurrenceError::VerificationFailed(
                    "invalid basis independence witness".into(),
                ));
            }
        }
    }
    for (i, inverse_row) in inverse.iter().enumerate() {
        for j in 0..rank {
            let sum: Q = (0..rank)
                .map(|t| &inverse_row[t] * &matrix[cert.independent_rows[t]][cert.pivot_columns[j]])
                .sum();
            if sum != if i == j { Q::one() } else { Q::zero() } {
                return Err(RecurrenceError::VerificationFailed(
                    "invalid rank minor inverse".into(),
                ));
            }
        }
    }
    Ok(basis)
}
