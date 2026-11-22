//! Exact Weyl coordinate identities and bounded coefficient fitting.
use super::*;

pub(super) fn decode_biv(p: &BivarPolyJson) -> Result<Biv, RecurrenceError> {
    if p.coeffs.len() > 9 || p.coeffs.iter().any(|r| r.len() > 9) {
        return Err(invalid("bivariate input degrees must be at most eight"));
    }
    Ok(Biv::new(
        p.coeffs
            .iter()
            .map(|r| numbers(r))
            .collect::<Result<Vec<_>, _>>()?,
    ))
}
pub(super) fn encode_biv(p: &Biv) -> BivarPolyJson {
    BivarPolyJson {
        coeffs: p.rows().iter().map(|r| strings(r)).collect(),
    }
}
pub(super) fn decode_coordinates(
    matrix: &RefinementCoordinates,
    rows: usize,
    cols: usize,
) -> Result<OperatorMatrix, RecurrenceError> {
    if matrix.len() != rows || matrix.iter().any(|r| r.len() != cols) {
        return Err(invalid("coordinate matrix dimensions differ"));
    }
    matrix
        .iter()
        .map(|r| {
            r.iter()
                .map(|coeffs| {
                    if coeffs.len() > 3 {
                        return Err(invalid("coordinate derivative order exceeds two"));
                    }
                    let mut op = Op::default();
                    for (d, p) in coeffs.iter().enumerate() {
                        op.add_term((0, d), decode_biv(p)?);
                    }
                    Ok(op)
                })
                .collect()
        })
        .collect()
}
pub(super) fn encode_coordinates(matrix: &OperatorMatrix) -> RefinementCoordinates {
    matrix
        .iter()
        .map(|r| {
            r.iter()
                .map(|op| {
                    let degree = op.0.keys().map(|(_, d)| *d).max().unwrap_or(0);
                    (0..=degree)
                        .map(|d| encode_biv(&op.0.get(&(0, d)).cloned().unwrap_or_default()))
                        .collect()
                })
                .collect()
        })
        .collect()
}
pub(super) fn identity(n: usize) -> OperatorMatrix {
    (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    if i == j {
                        Op::identity()
                    } else {
                        Op::default()
                    }
                })
                .collect()
        })
        .collect()
}
pub(super) fn multiply(a: &OperatorMatrix, b: &OperatorMatrix) -> OperatorMatrix {
    (0..a.len())
        .map(|i| {
            (0..b[0].len())
                .map(|j| (0..b.len()).fold(Op::default(), |s, k| s.add(&a[i][k].compose(&b[k][j]))))
                .collect()
        })
        .collect()
}
pub(super) fn apply(matrix: &OperatorMatrix, state: &[Vec<Q>], index: usize) -> Vec<Vec<Q>> {
    let n = Q::from_integer(index.into());
    matrix
        .iter()
        .map(|row| {
            row.iter().zip(state).fold(vec![Q::zero()], |sum, (op, p)| {
                op.0.iter().fold(sum, |sum, ((_, d), c)| {
                    let mut derivative = p.clone();
                    for _ in 0..*d {
                        derivative = vp::derivative(&derivative);
                    }
                    vp::add(&sum, &vp::mul(&c.evaluate_index(&n), &derivative))
                })
            })
        })
        .collect()
}
pub(super) fn decode_system(
    json: &TransformVectorRecurrenceJson,
) -> Result<(OperatorMatrix, Biv, Vec<Vec<Q>>), RecurrenceError> {
    let n = json.matrix.len();
    if n == 0
        || n > 5
        || json.initial_state.len() != n
        || json.matrix.iter().any(|r| r.len() != n)
        || json.forcing.is_some()
    {
        return Err(invalid(
            "system must be square, homogeneous, dimension 1..=5, with a matching seed",
        ));
    }
    let matrix = decode_coordinates(&json.matrix, n, n)?;
    let denominator = decode_biv(&json.denominator)?;
    if denominator.is_zero() || denominator.variable_degree().is_some_and(|d| d > 0) {
        return Err(invalid(
            "system denominator must be nonzero and independent of x",
        ));
    }
    let seed = json
        .initial_state
        .iter()
        .map(|r| {
            if r.len() > 65 {
                return Err(invalid("seed degree exceeds 64"));
            }
            numbers(r)
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok((matrix, denominator, seed))
}
pub(super) fn encode_system(
    matrix: &OperatorMatrix,
    denominator: &Biv,
    seed: &[Vec<Q>],
    first: usize,
) -> TransformVectorRecurrenceJson {
    TransformVectorRecurrenceJson {
        first_index: first,
        initial_state: seed.iter().map(|p| strings(p)).collect(),
        denominator: encode_biv(denominator),
        matrix: encode_coordinates(matrix),
        forcing: None,
    }
}
pub(super) fn companion(
    record: &RecurrenceJson,
    maximum: usize,
) -> Result<TransformVectorRecurrenceJson, RecurrenceError> {
    let (rule, first, seeds) = record.to_recurrence_parts()?;
    if rule.inhomogeneous.is_some()
        || rule
            .terms
            .iter()
            .any(|t| t.is_alternating() || t.offset == 0)
    {
        return Err(invalid(
            "companion requires a homogeneous nonalternating rule with positive lags",
        ));
    }
    let lag = rule.terms.iter().map(|t| t.offset).max().unwrap_or(0);
    if lag == 0 || lag > maximum || seeds.len() < lag || seeds.len() > 64 {
        return Err(invalid("companion lag or seed count exceeds limits"));
    }
    let start = first
        .checked_add(seeds.len() - 1)
        .ok_or_else(|| invalid("index overflow"))?;
    let denominator = rule
        .denominator
        .as_ref()
        .map(Biv::from_bivar)
        .unwrap_or_else(Biv::one)
        .shift_index(1);
    if denominator.variable_degree().is_some_and(|d| d > 0) {
        return Err(invalid("companion requires an index-only denominator"));
    }
    let mut matrix = vec![vec![Op::default(); lag]; lag];
    for (i, row) in matrix.iter_mut().enumerate().take(lag - 1) {
        row[i + 1] = Op::single(0, 0, denominator.clone());
    }
    for term in &rule.terms {
        matrix[lag - 1][lag - term.offset].add_term(
            (0, term.deriv_order),
            Biv::from_bivar(&term.coeff).shift_index(1),
        );
    }
    Ok(encode_system(
        &matrix,
        &denominator,
        &seeds[seeds.len() - lag..],
        start,
    ))
}
pub(super) fn denominator_positive(p: &Biv, first: usize) -> bool {
    crate::decide_index_rational_function_sign(
        &p.variable_coefficient(0),
        &[Q::one()],
        &crate::IntegerRange::starting_at(first),
        crate::IndexSignCondition::Positive,
    )
    .holds()
}

fn monomial(n: usize, x: usize, d: usize) -> Op {
    let mut rows = vec![vec![Q::zero(); x + 1]; n + 1];
    rows[n][x] = Q::one();
    Op::single(0, d, Biv::new(rows))
}
fn flatten(ops: &[Op]) -> std::collections::BTreeMap<(usize, usize, usize, usize), Q> {
    let mut out = std::collections::BTreeMap::new();
    for (j, op) in ops.iter().enumerate() {
        for ((_, d), p) in &op.0 {
            for (n, row) in p.rows().iter().enumerate() {
                for (x, c) in row.iter().enumerate() {
                    if !c.is_zero() {
                        out.insert((j, *d, n, x), c.clone());
                    }
                }
            }
        }
    }
    out
}
/// Solve T(n+1) C(n)=M(n) T(n) as a polynomial Ore identity.
pub(super) fn fit_symbolic(
    t: &OperatorMatrix,
    c: &OperatorMatrix,
    limits: &RefinementSearchLimits,
) -> Result<Option<OperatorMatrix>, RecurrenceError> {
    let shifted = t
        .iter()
        .map(|r| r.iter().map(|op| op.shift_index(1)).collect())
        .collect();
    let lhs = multiply(&shifted, c);
    let dim = t.len();
    for derivative in 0..=limits.derivative_order {
        let support: Vec<_> = (0..dim)
            .flat_map(|k| {
                (0..=derivative).flat_map(move |d| {
                    (0..=limits.index_degree)
                        .flat_map(move |n| (0..=limits.variable_degree).map(move |x| (k, d, n, x)))
                })
            })
            .collect();
        let columns: Vec<_> = support
            .iter()
            .map(|&(k, d, n, x)| {
                flatten(
                    &t[k]
                        .iter()
                        .map(|op| monomial(n, x, d).compose(op))
                        .collect::<Vec<_>>(),
                )
            })
            .collect();
        let mut matrix = Vec::new();
        let mut failed = false;
        for target in &lhs {
            let rhs = flatten(target);
            let mut keys: std::collections::BTreeSet<_> = rhs.keys().copied().collect();
            for col in &columns {
                keys.extend(col.keys().copied());
            }
            if keys.len().saturating_mul(support.len()) > limits.maximum_fit_cells {
                return Err(invalid("symbolic fit exceeds maximum_fit_cells"));
            }
            let equations: Vec<Vec<_>> = keys
                .iter()
                .map(|key| {
                    columns
                        .iter()
                        .map(|col| col.get(key).cloned().unwrap_or_default())
                        .collect()
                })
                .collect();
            let values: Vec<_> = keys
                .iter()
                .map(|key| rhs.get(key).cloned().unwrap_or_default())
                .collect();
            let solution = if equations.is_empty() {
                Some(vec![Q::zero(); support.len()])
            } else {
                nonnegative_fit(&equations, &values, limits.maximum_supports)
            };
            let Some(solution) = solution else {
                failed = true;
                break;
            };
            let mut row = vec![Op::default(); dim];
            for (&(k, d, n, x), value) in support.iter().zip(solution) {
                row[k] = row[k].add(&monomial(n, x, d).scale(&value));
            }
            matrix.push(row);
        }
        if !failed {
            return Ok(Some(matrix));
        }
    }
    Ok(None)
}

/// Bounded active-support search. Each proposal is an exact solution; no
/// completeness claim is made when the support cap is reached. Keep the
/// original particular solution for diagnostics if no nonnegative one is found.
fn nonnegative_fit(equations: &[Vec<Q>], values: &[Q], maximum: usize) -> Option<Vec<Q>> {
    let width = equations.first()?.len();
    let original = crate::linalg::solve_linear_system_with_rank(equations, values)?.particular;
    if original.iter().all(|c| !c.is_negative()) {
        return Some(original);
    }
    let mut queue = std::collections::VecDeque::from([Vec::<usize>::new()]);
    let mut seen = std::collections::BTreeSet::new();
    let mut steps = 0;
    while let Some(excluded) = queue.pop_front() {
        if !seen.insert(excluded.clone()) {
            continue;
        }
        if steps >= maximum {
            break;
        }
        steps += 1;
        let support: Vec<_> = (0..width).filter(|j| !excluded.contains(j)).collect();
        if support.is_empty() {
            continue;
        }
        let reduced: Vec<Vec<_>> = equations
            .iter()
            .map(|r| support.iter().map(|&j| r[j].clone()).collect())
            .collect();
        let Some(solution) = crate::linalg::solve_linear_system_with_rank(&reduced, values) else {
            continue;
        };
        if solution.particular.iter().all(|c| !c.is_negative()) {
            let mut full = vec![Q::zero(); width];
            for (&j, c) in support.iter().zip(solution.particular) {
                full[j] = c;
            }
            return Some(full);
        }
        for (&j, c) in support.iter().zip(&solution.particular) {
            if c.is_negative() {
                let mut next = excluded.clone();
                next.push(j);
                next.sort_unstable();
                queue.push_back(next);
            }
        }
    }
    Some(original)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn active_support_finds_a_nonnegative_solution_beyond_the_particular_one() {
        let q = |n: i64| Q::from_integer(n.into());
        assert_eq!(
            nonnegative_fit(&[vec![q(1), q(-1)]], &[q(-1)], 4),
            Some(vec![q(0), q(1)])
        );
        assert_eq!(
            nonnegative_fit(&[vec![q(1), q(-1)]], &[q(-1)], 1),
            Some(vec![q(-1), q(0)])
        );
    }
}
