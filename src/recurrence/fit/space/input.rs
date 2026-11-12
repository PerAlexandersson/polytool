//! Recurrence-space input preparation and residual matrices.

use super::*;

pub(super) fn checked_q(s: &str, limits: &RecurrenceSpaceLimits) -> Result<Q, RecurrenceError> {
    if s.len() as u64 > limits.max_coefficient_bits.saturating_add(32) {
        return Err(RecurrenceError::BudgetExhausted(
            "coefficient text exceeds bit budget".into(),
        ));
    }
    let q = parse_rational_coeff(s)?;
    if q.numer().bits().max(q.denom().bits()) > limits.max_coefficient_bits {
        return Err(RecurrenceError::BudgetExhausted(
            "coefficient exceeds bit budget".into(),
        ));
    }
    Ok(q)
}
pub(super) fn parse_matrix(
    rows: &ExactRows,
    limits: &RecurrenceSpaceLimits,
) -> Result<Vec<Vec<Q>>, RecurrenceError> {
    if rows
        .len()
        .checked_mul(rows.iter().map(Vec::len).max().unwrap_or(0))
        .is_none_or(|n| n > limits.max_matrix_entries)
    {
        return Err(RecurrenceError::BudgetExhausted(
            "matrix exceeds entry budget".into(),
        ));
    }
    rows.iter()
        .map(|r| r.iter().map(|x| checked_q(x, limits)).collect())
        .collect()
}
fn parse_polys(
    rows: &ExactRows,
    limits: &RecurrenceSpaceLimits,
) -> Result<Vec<Poly>, RecurrenceError> {
    if rows.len() > limits.max_rows
        || rows
            .iter()
            .any(|r| r.len().saturating_sub(1) > limits.max_degree)
    {
        return Err(RecurrenceError::BudgetExhausted(
            "input rows exceed row/degree budget".into(),
        ));
    }
    Ok(parse_matrix(rows, limits)?
        .into_iter()
        .map(Poly::new)
        .collect())
}
pub(super) fn poly_budget(p: &Poly, limits: &RecurrenceSpaceLimits) -> Result<(), RecurrenceError> {
    if p.degree().unwrap_or(0) > limits.max_degree {
        return Err(RecurrenceError::BudgetExhausted(
            "polynomial degree budget".into(),
        ));
    }
    if p.coeffs()
        .iter()
        .any(|x| x.numer().bits().max(x.denom().bits()) > limits.max_coefficient_bits)
    {
        return Err(RecurrenceError::BudgetExhausted(
            "coefficient bit budget".into(),
        ));
    }
    Ok(())
}

pub(super) struct Input {
    pub(super) a: Option<Vec<Poly>>,
    pub(super) indexed_a: Option<Vec<BivarPoly>>,
    pub(super) rows: Vec<Poly>,
    pub(super) origin: usize,
    consistency: Vec<Poly>,
    pub(super) rescaling: Option<RecurrenceCoefficientRescaling>,
}

pub(super) fn prepare(request: &RecurrenceSpaceRequest) -> Result<Input, RecurrenceError> {
    let lim = &request.limits;
    if request.shapes.len() > lim.max_shapes {
        return Err(RecurrenceError::LimitExceeded(
            "too many requested shapes".into(),
        ));
    }
    if let Some(rescaling) = &request.coefficient_rescaling {
        rescaling.validate(lim)?;
        if request.source.is_none() && request.source_recurrence.is_none() {
            return Err(RecurrenceError::InvalidInput(
                "coefficient_rescaling requires a constant-index monic source \
                        (source or source_recurrence)"
                    .into(),
            ));
        }
    }
    if usize::from(request.source.is_some())
        + usize::from(request.source_recurrence.is_some())
        + usize::from(request.indexed_source.is_some())
        > 1
    {
        return Err(RecurrenceError::InvalidInput(
            "provide exactly one source form".into(),
        ));
    }
    if let Some(src) = &request.indexed_source {
        if src.coefficients.is_empty()
            || src.coefficients.len() > lim.max_rows
            || src.seeds.len() != src.coefficients.len()
        {
            return Err(RecurrenceError::InvalidInput(
                "indexed source needs exactly k seeds for positive order k".into(),
            ));
        }
        if request.domain_start < src.origin {
            return Err(RecurrenceError::InvalidInput(
                "domain starts before source origin".into(),
            ));
        }
        if request.rows.is_some() && request.row_origin != src.origin {
            return Err(RecurrenceError::InvalidInput(
                "consistency row_origin must equal source origin".into(),
            ));
        }
        return Ok(Input {
            a: None,
            indexed_a: Some(
                src.coefficients
                    .iter()
                    .map(|p| propagation::parse(p, lim))
                    .collect::<Result<_, _>>()?,
            ),
            rows: parse_polys(&src.seeds, lim)?,
            origin: src.origin,
            consistency: request
                .rows
                .as_ref()
                .map(|p| parse_polys(p, lim))
                .transpose()?
                .unwrap_or_default(),
            rescaling: None,
        });
    }
    if let Some(json) = &request.source_recurrence {
        parse_polys(&json.initial_polynomials, lim)?;
        if json.recurrence.terms.len() > lim.max_matrix_entries
            || json
                .recurrence
                .terms
                .iter()
                .any(|t| t.offset > lim.max_rows)
        {
            return Err(RecurrenceError::LimitExceeded(
                "source recurrence exceeds limits".into(),
            ));
        }
        for p in json
            .recurrence
            .terms
            .iter()
            .map(|t| &t.coeff)
            .chain(json.recurrence.denominator.iter())
            .chain(json.recurrence.inhomogeneous.iter())
        {
            if p.coeffs.len().saturating_sub(1) > lim.max_degree
                || p.coeffs
                    .iter()
                    .any(|r| r.len().saturating_sub(1) > lim.max_degree)
            {
                return Err(RecurrenceError::BudgetExhausted(
                    "source coefficient degree budget".into(),
                ));
            }
            parse_matrix(&p.coeffs, lim)?;
        }
    }
    let adapted = request
        .source_recurrence
        .as_ref()
        .map(ConstantLagSource::from_recurrence_json)
        .transpose()?;
    let source = request
        .source
        .as_ref()
        .or_else(|| adapted.as_ref().map(|x| &x.0));
    let supplied = request
        .rows
        .as_ref()
        .or_else(|| adapted.as_ref().map(|x| &x.1));
    if let Some(src) = source {
        let a = parse_polys(&src.coefficients, lim)?;
        let rows = parse_polys(&src.seeds, lim)?;
        if a.is_empty() || rows.len() != a.len() {
            return Err(RecurrenceError::InvalidInput(
                "source requires exactly k seeds for positive order k".into(),
            ));
        }
        if request.domain_start < src.origin {
            return Err(RecurrenceError::InvalidInput(
                "domain starts before the source origin".into(),
            ));
        }
        if request.rows.is_some() && request.row_origin != src.origin {
            return Err(RecurrenceError::InvalidInput(
                "consistency row_origin must equal source origin".into(),
            ));
        }
        let mut consistency = supplied
            .map(|r| parse_polys(r, lim))
            .transpose()?
            .unwrap_or_default();
        if let Some((_, initial)) = &adapted {
            let original = parse_polys(initial, lim)?;
            if !original.iter().zip(&consistency).all(|(a, b)| a == b) {
                return Err(RecurrenceError::InvalidInput(
                    "supplied rows disagree with source_recurrence initial rows".into(),
                ));
            }
            if original.len() > consistency.len() {
                consistency = original;
            }
        }
        Ok(Input {
            a: Some(a),
            indexed_a: None,
            rows,
            origin: src.origin,
            consistency,
            rescaling: request.coefficient_rescaling,
        })
    } else {
        let rows = supplied.ok_or(RecurrenceError::InvalidInput(
            "provide a source recurrence or finite rows".into(),
        ))?;
        if request.domain_start < request.row_origin {
            return Err(RecurrenceError::InvalidInput(
                "domain starts before row_origin".into(),
            ));
        }
        Ok(Input {
            a: None,
            indexed_a: None,
            rows: parse_polys(rows, lim)?,
            origin: request.row_origin,
            consistency: Vec::new(),
            rescaling: None,
        })
    }
}

pub(super) fn source_rows(
    input: &Input,
    count: usize,
    limits: &RecurrenceSpaceLimits,
) -> Result<Vec<Poly>, RecurrenceError> {
    let needed = count.max(input.consistency.len());
    if needed > limits.max_rows {
        return Err(RecurrenceError::BudgetExhausted("source row budget".into()));
    }
    let mut rows = input.rows.clone();
    let mut operations = limits.max_operations;
    if input.a.is_some() || input.indexed_a.is_some() {
        while rows.len() < needed {
            let order = input
                .a
                .as_ref()
                .map(Vec::len)
                .unwrap_or_else(|| input.indexed_a.as_ref().unwrap().len());
            let start = rows.len() - order;
            let evaluated;
            let a = if let Some(a) = &input.a {
                a
            } else {
                let rank =
                    input
                        .origin
                        .checked_add(start)
                        .ok_or(RecurrenceError::LimitExceeded(
                            "source rank overflow".into(),
                        ))?;
                let source = input.indexed_a.as_ref().unwrap();
                let work = source
                    .iter()
                    .flat_map(|p| &p.coeffs)
                    .map(Vec::len)
                    .sum::<usize>();
                operations =
                    operations
                        .checked_sub(work)
                        .ok_or(RecurrenceError::BudgetExhausted(
                            "source arithmetic budget".into(),
                        ))?;
                evaluated = source
                    .iter()
                    .map(|p| propagation::evaluate(p, rank, limits))
                    .collect::<Result<Vec<_>, _>>()?;
                &evaluated
            };
            let mut p = Poly::zero();
            for (i, c) in a.iter().enumerate() {
                let work = c
                    .coeffs()
                    .len()
                    .checked_mul(rows[start + i].coeffs().len())
                    .ok_or(RecurrenceError::LimitExceeded(
                        "source work overflow".into(),
                    ))?;
                operations =
                    operations
                        .checked_sub(work)
                        .ok_or(RecurrenceError::BudgetExhausted(
                            "source arithmetic budget".into(),
                        ))?;
                if c.degree()
                    .unwrap_or(0)
                    .saturating_add(rows[start + i].degree().unwrap_or(0))
                    > limits.max_degree
                {
                    return Err(RecurrenceError::BudgetExhausted(
                        "source degree budget".into(),
                    ));
                }
                p = p + c.clone() * rows[start + i].clone();
                poly_budget(&p, limits)?;
            }
            rows.push(p);
        }
        if !input.consistency.iter().zip(&rows).all(|(a, b)| a == b) {
            return Err(RecurrenceError::InvalidInput(
                "supplied rows disagree with the source recurrence and seeds".into(),
            ));
        }
    }
    Ok(rows)
}

pub(super) fn shape_columns(
    shape: &RecurrenceSpaceShape,
    limits: &RecurrenceSpaceLimits,
) -> Result<(usize, Vec<RecurrenceMonomial>), RecurrenceError> {
    let (target, mut cols) = match shape {
        RecurrenceSpaceShape::AutomaticVariableDegree { .. } => {
            return Err(RecurrenceError::InvalidInput(
                "unresolved automatic shape".into(),
            ))
        }
        RecurrenceSpaceShape::Explicit {
            target_lag,
            monomials,
        } => (*target_lag, monomials.clone()),
        RecurrenceSpaceShape::Box {
            length,
            derivative_order,
            index_degree,
            variable_degree,
            leading_index_degree,
            leading_variable_degree,
        } => {
            let dimensions = [
                *length,
                derivative_order.saturating_add(1),
                index_degree.saturating_add(1),
                variable_degree.saturating_add(1),
            ];
            let rhs = dimensions.iter().try_fold(1usize, |a, b| a.checked_mul(*b));
            let lhs = leading_index_degree.checked_add(1).and_then(|n| {
                leading_variable_degree
                    .checked_add(1)
                    .and_then(|x| n.checked_mul(x))
            });
            if rhs
                .and_then(|r| lhs.and_then(|l| r.checked_add(l)))
                .is_none_or(|n| n > limits.max_columns)
            {
                return Err(RecurrenceError::BudgetExhausted("column budget".into()));
            }
            let mut cols = Vec::new();
            for n in 0..=*leading_index_degree {
                for x in 0..=*leading_variable_degree {
                    cols.push(RecurrenceMonomial {
                        lag: *length,
                        derivative_order: 0,
                        index_degree: n,
                        variable_degree: x,
                    });
                }
            }
            for lag in 0..*length {
                for d in 0..=*derivative_order {
                    for n in 0..=*index_degree {
                        for x in 0..=*variable_degree {
                            cols.push(RecurrenceMonomial {
                                lag,
                                derivative_order: d,
                                index_degree: n,
                                variable_degree: x,
                            });
                        }
                    }
                }
            }
            (*length, cols)
        }
    };
    if target == 0
        || cols
            .iter()
            .any(|c| c.lag > target || (c.lag == target && c.derivative_order != 0))
    {
        return Err(RecurrenceError::InvalidInput(
            "invalid successor support".into(),
        ));
    }
    if cols.len() > limits.max_columns
        || cols.iter().any(|c| {
            [c.index_degree, c.variable_degree, c.derivative_order]
                .iter()
                .any(|d| *d > limits.max_degree)
        })
    {
        return Err(RecurrenceError::BudgetExhausted(
            "support degree/column budget".into(),
        ));
    }
    cols.sort();
    if cols.windows(2).any(|w| w[0] == w[1]) {
        return Err(RecurrenceError::InvalidInput(
            "duplicate support monomial".into(),
        ));
    }
    Ok((target, cols))
}

// Evaluate directly on rows; no generating-function numerator construction.
pub(super) fn column_at(
    c: &RecurrenceMonomial,
    rows: &[Poly],
    rank: usize,
    origin: usize,
    limits: &RecurrenceSpaceLimits,
) -> Result<Poly, RecurrenceError> {
    let row = rank
        .checked_sub(origin)
        .and_then(|n| n.checked_add(c.lag))
        .and_then(|i| rows.get(i))
        .ok_or(RecurrenceError::InvalidInput("missing row".into()))?;
    let mut p = if c.derivative_order > row.degree().unwrap_or(0) {
        Poly::zero()
    } else {
        row.nth_derivative(c.derivative_order)
    };
    if p.degree().unwrap_or(0).saturating_add(c.variable_degree) > limits.max_degree {
        return Err(RecurrenceError::BudgetExhausted(
            "residual degree budget".into(),
        ));
    }
    let power = u32::try_from(c.index_degree)
        .map_err(|_| RecurrenceError::LimitExceeded("index degree overflow".into()))?;
    let factor = Q::from_integer(BigInt::from(rank).pow(power));
    p = p.scale(&factor) * Poly::monomial(Q::one(), c.variable_degree);
    poly_budget(&p, limits)?;
    Ok(p)
}

/// The residual matrix of the support `cols` at the `count` base indices
/// from `start`: one row per `X`-coefficient of each residual, one column per
/// monomial. Columns at lag `target` enter positively, the others negatively.
pub(in crate::recurrence) fn residual_matrix(
    rows: &[Poly],
    origin: usize,
    start: usize,
    count: usize,
    target: usize,
    cols: &[RecurrenceMonomial],
    lim: &RecurrenceSpaceLimits,
) -> Result<Vec<Vec<Q>>, RecurrenceError> {
    let mut matrix = Vec::new();
    let mut operations = lim.max_operations;
    for offset in 0..count {
        let n = start
            .checked_add(offset)
            .ok_or(RecurrenceError::LimitExceeded("rank overflow".into()))?;
        for c in cols {
            let p = rows
                .get(
                    n.checked_sub(origin)
                        .and_then(|i| i.checked_add(c.lag))
                        .ok_or(RecurrenceError::LimitExceeded("rank overflow".into()))?,
                )
                .ok_or(RecurrenceError::InvalidInput("missing row".into()))?;
            let work = p
                .coeffs()
                .len()
                .checked_mul(c.derivative_order.saturating_add(2))
                .ok_or(RecurrenceError::LimitExceeded(
                    "residual work overflow".into(),
                ))?;
            operations = operations
                .checked_sub(work)
                .ok_or(RecurrenceError::BudgetExhausted(
                    "residual arithmetic budget".into(),
                ))?;
        }
        let values: Vec<_> = cols
            .iter()
            .map(|c| {
                column_at(c, rows, n, origin, lim).map(|p| {
                    if c.lag == target {
                        p
                    } else {
                        p.scale(&(-Q::one()))
                    }
                })
            })
            .collect::<Result<_, _>>()?;
        let degree = values.iter().filter_map(Poly::degree).max();
        if let Some(degree) = degree {
            if matrix
                .len()
                .checked_add(degree + 1)
                .and_then(|r| r.checked_mul(cols.len()))
                .is_none_or(|x| x > lim.max_matrix_entries)
            {
                return Err(RecurrenceError::BudgetExhausted(
                    "matrix entry budget".into(),
                ));
            }
            for x in 0..=degree {
                matrix.push(values.iter().map(|p| p.coeff(x)).collect());
            }
        }
    }
    Ok(matrix)
}
