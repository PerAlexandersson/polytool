//! Exact hypotheses supporting the enclosing theorem.

use super::*;

/// `a_0 + a_1 n` with integer coefficients through every point, from the
/// first two; `None` when there are fewer than two points or no such line.
pub(super) fn fit_affine(points: &[(usize, Q)]) -> Option<Affine> {
    let (n0, y0) = points.first()?;
    let (n1, y1) = points.get(1)?;
    let slope = (y1 - y0) / Q::from_integer(BigInt::from(n1 - n0));
    let intercept = y0 - &slope * Q::from_integer(BigInt::from(*n0));
    let line = (intercept, slope);
    (affine_is_integral(&line)
        && points
            .iter()
            .all(|(n, y)| *y == &line.0 + &line.1 * Q::from_integer(BigInt::from(*n))))
    .then_some(line)
}

/// The valuation bound `val P_n >= v(n)` (so `[x^k] P_n = 0` for
/// `k < v(n)`): by induction on the rule when every term `c x^e D^d P_(n-r)`
/// has `e - d - v_1 r >= 0` and the rows before the rule start are checked
/// exactly; otherwise, for a constant valuation `v`, by the relation with
/// `p(n, k) != 0` for `0 <= k < v`.
pub(super) fn lower_support_hypothesis(
    shape: &RecurrenceShape,
    p: &LinearFactors,
    valuation: &Affine,
    earlier: &[(usize, &Vec<Q>)],
    first: usize,
) -> LiftHypothesis {
    let name = "lower_support";
    let statement = "[x^k] P_n = 0 for k < v(n), the valuation";
    let terms_hold = shape.terms.iter().all(|(&(lag, order), coefficient)| {
        (0..x_width(coefficient))
            .filter(|power| !x_column(coefficient, *power).is_empty())
            .all(|power| {
                Q::from_integer(BigInt::from(power))
                    - Q::from_integer(BigInt::from(order))
                    - &valuation.1 * Q::from_integer(BigInt::from(lag))
                    >= Q::zero()
            })
    });
    let rows_hold = earlier.iter().all(|(n, row)| {
        let bound = &valuation.0 + &valuation.1 * Q::from_integer(BigInt::from(*n));
        (0..row.len()).all(|k| row[k].is_zero() || Q::from_integer(BigInt::from(k)) >= bound)
    });
    if terms_hold && rows_hold {
        return exact_hypothesis(
            name,
            statement,
            Ok(format!(
                "induction on the rule: every term x^e D^d P_(n-r) has e - d - ({}) r >= 0, and \
                 the rows before the rule start are checked exactly",
                affine_text(&(valuation.1.clone(), Q::zero()))
            )),
        );
    }
    let Some(v) = (valuation.1.is_zero())
        .then(|| valuation.0.to_integer().to_usize())
        .flatten()
    else {
        return unresolved_hypothesis(
            name,
            statement,
            "neither the rule nor the relation bounds the valuation".into(),
        );
    };
    let mut product = IndexFunction::constant(Q::one());
    for k in 0..v {
        product = times(
            &product,
            &IndexFunction::polynomial(p.at(&affine_constant(k as i64))),
        );
    }
    sign_hypothesis(
        name,
        &format!("{statement}: p(n, k) != 0 for 0 <= k < {v}"),
        "prod_(k < v) p(n, k)",
        &product,
        first,
        IndexSignCondition::Nonzero,
    )
}

/// The degree bound `deg P_n <= d(n)` (so `[x^k] P_n = 0` for `k > d(n)`):
/// by induction on the rule when every term `c x^e D^d P_(n-r)` has
/// `e - d - d_1 r <= 0` and the rows before the rule start are checked
/// exactly; otherwise by the relation, with the factor `k - d(n) - 1` of `q`,
/// `c_p(n) != 0` and `u(n) + d(n) + 1 > 0` for every factor `k + u` of `p`.
pub(super) fn upper_support_hypotheses(
    shape: &RecurrenceShape,
    p: &LinearFactors,
    q: &LinearFactors,
    degree: &Affine,
    earlier: &[(usize, &Vec<Q>)],
    first: usize,
) -> Vec<LiftHypothesis> {
    let name = "upper_support";
    let statement = "[x^k] P_n = 0 for k > d(n), the degree";
    let terms_hold = shape.terms.iter().all(|(&(lag, order), coefficient)| {
        x_degree(coefficient).is_none_or(|power| {
            Q::from_integer(BigInt::from(power))
                - Q::from_integer(BigInt::from(order))
                - &degree.1 * Q::from_integer(BigInt::from(lag))
                <= Q::zero()
        })
    });
    let rows_hold = earlier.iter().all(|(n, row)| {
        Q::from_integer(BigInt::from(row.len().saturating_sub(1)))
            <= &degree.0 + &degree.1 * Q::from_integer(BigInt::from(*n))
    });
    if terms_hold && rows_hold {
        return vec![exact_hypothesis(
            name,
            statement,
            Ok(format!(
                "induction on the rule: every term x^e D^d P_(n-r) has e - d - ({}) r <= 0, and \
                 the rows before the rule start are checked exactly",
                affine_text(&(degree.1.clone(), Q::zero()))
            )),
        )];
    }
    let top = affine_neg(&affine_add(degree, &affine_constant(1)));
    let mut out = vec![exact_hypothesis(
        &format!("{name}_factor"),
        "q(n, k) has the factor k - d(n) - 1",
        if q.roots.contains(&top) {
            Ok(format!("d(n) = {}", affine_text(degree)))
        } else {
            Err((
                None,
                format!("no factor k - d(n) - 1, d(n) = {}", affine_text(degree)),
            ))
        },
    )];
    out.push(sign_hypothesis(
        &format!("{name}_scalar"),
        "c_p(n) != 0, the scalar factor of p",
        "c_p(n)",
        &IndexFunction::polynomial(p.scalar.clone()),
        first,
        IndexSignCondition::Nonzero,
    ));
    // A crude degree bound D(n) = D_0 + beta n from the rule: beta =
    // max (e - d) / r over the terms, D_0 covering the rows checked exactly.
    let beta = shape
        .terms
        .iter()
        .filter_map(|(&(lag, order), coefficient)| {
            x_degree(coefficient).map(|power| {
                (Q::from_integer(BigInt::from(power)) - Q::from_integer(BigInt::from(order)))
                    / Q::from_integer(BigInt::from(lag))
            })
        })
        .fold(
            degree.1.clone(),
            |acc, value| if value > acc { value } else { acc },
        );
    let crude = earlier
        .iter()
        .map(|(n, row)| {
            Q::from_integer(BigInt::from(row.len().saturating_sub(1)))
                - &beta * Q::from_integer(BigInt::from(*n))
        })
        .reduce(|a, b| if a > b { a } else { b })
        .map(|intercept| (intercept, beta.clone()));
    for (index, u) in p.roots.iter().enumerate() {
        let below = sign_hypothesis(
            &format!("{name}_root_{index}"),
            "u(n) + d(n) + 1 > 0, so p(n, k) != 0 for every integer k > d(n) and [x^k] P_n = 0 \
             there",
            "u(n) + d(n) + 1",
            &affine_function(&affine_add(&affine_add(u, degree), &affine_constant(1))),
            first,
            IndexSignCondition::Positive,
        );
        if below.status == LiftHypothesisStatus::Verified {
            out.push(below);
            continue;
        }
        let above = crude.as_ref().map(|crude| {
            sign_hypothesis(
                &format!("{name}_root_{index}"),
                &format!(
                    "u(n) + D(n) < 0 for the degree bound D(n) = {} of the rule (deg P_n <= D(n) \
                     by induction: D(n) >= e - d + D(n - r) for every term x^e D^d P_(n-r), the \
                     rows before checked exactly), so the root -u(n) of p lies above every \
                     degree and p(n, k) != 0 for d(n) < k <= D(n)",
                    IndexFunction::polynomial(vec![crude.0.clone(), crude.1.clone()]).text()
                ),
                "u(n) + D(n)",
                &affine_function(&affine_add(u, crude)),
                first,
                IndexSignCondition::Negative,
            )
        });
        match above {
            Some(above) if above.status == LiftHypothesisStatus::Verified => out.push(above),
            _ => out.push(below),
        }
    }
    out
}

/// The hypotheses of a matching (module documentation), decided for
/// `n >= first`.
pub(super) fn matching_hypotheses(
    ratio: &Ratio,
    matching: &Matching,
    first: usize,
    upper: Option<&Affine>,
) -> Vec<LiftHypothesis> {
    let mut out = Vec::new();
    let one = affine_constant(1);
    let size = complement_size(&ratio.numerators[matching.base_numerator])
        .unwrap_or_else(|| affine_constant(0));
    out.push(sign_hypothesis(
        "base_size",
        &format!(
            "the base factor {} and the factor j give binom(M, j) with M(n) >= 0",
            factor_text(&ratio.numerators[matching.base_numerator])
        ),
        "M(n)",
        &affine_function(&size),
        first,
        IndexSignCondition::Nonnegative,
    ));
    if let Some(upper) = upper {
        out.push(sign_hypothesis(
            "base_within_degree",
            "M(n) <= d(n) - v (reversal at the degree)",
            "d(n) - v - M(n)",
            &affine_function(&affine_add(upper, &affine_neg(&size))),
            first,
            IndexSignCondition::Nonnegative,
        ));
    }
    for index in &matching.complements {
        let other =
            complement_size(&ratio.numerators[*index]).unwrap_or_else(|| affine_constant(0));
        out.push(sign_hypothesis(
            &format!("complement_{index}"),
            &format!(
                "the factor {} gives k -> M'!/(M' - k)! with M'(n) >= M(n) (Polya--Schur)",
                factor_text(&ratio.numerators[*index])
            ),
            "M'(n) - M(n)",
            &affine_function(&affine_add(&other, &affine_neg(&size))),
            first,
            IndexSignCondition::Nonnegative,
        ));
    }
    for (numerator, denominator) in &matching.pairs {
        let (a, b) = (
            &ratio.numerators[*numerator],
            &ratio.denominators[*denominator],
        );
        out.push(sign_hypothesis(
            &format!("pair_{numerator}_{denominator}_length"),
            &format!(
                "({}) / ({}) pairs to the polynomial (k + b + 1) ... (k + a): a - b is an \
                 integer >= 0 (Laguerre)",
                factor_text(a),
                factor_text(b)
            ),
            "a(n) - b(n)",
            &affine_function(&affine_add(a, &affine_neg(b))),
            first,
            IndexSignCondition::Nonnegative,
        ));
        out.push(sign_hypothesis(
            &format!("pair_{numerator}_{denominator}_roots"),
            "b(n) + 1 > 0: every root of the paired polynomial is negative",
            "b(n) + 1",
            &affine_function(&affine_add(b, &one)),
            first,
            IndexSignCondition::Positive,
        ));
    }
    for index in &matching.reciprocals {
        let b = &ratio.denominators[*index];
        out.push(sign_hypothesis(
            &format!("reciprocal_{index}"),
            &format!(
                "the factor {} gives k -> 1/(b + 1)_k with b(n) + 1 > 0",
                factor_text(b)
            ),
            "b(n) + 1",
            &affine_function(&affine_add(b, &one)),
            first,
            IndexSignCondition::Positive,
        ));
    }
    let sign = if (matching.complements.len() + 1) % 2 == 0 {
        Q::one()
    } else {
        -Q::one()
    };
    out.push(sign_hypothesis(
        "scale_positive",
        "C_eff(n) = C(n) (-1)^(complements) > 0 (the scaling x -> C_eff x)",
        "C_eff(n)",
        &ratio.scale.scale(&sign),
        first,
        IndexSignCondition::Positive,
    ));
    out
}

/// `H(n) P_n = 0` for every `n >= from`: an invariant annihilator system
/// of the rule whose module contains `H` (module documentation).
pub(super) fn relation_hypothesis(
    operator: &DOp,
    prepared: &crate::recurrence::lift::hadamard_factor::prepare::PreparedLiftRule,
    from: usize,
    rows: &[Vec<Q>],
    first_index: usize,
) -> LiftHypothesis {
    let name = "row_relation";
    let statement = "H(n) P_n = 0 for every n >= the first base, H = p(n, theta) + x q(n, theta + \
                     1): an invariant row-annihilator system of the rule contains H";
    // A cheap rejection on the rows first.
    if let Some((n, _)) = rows
        .iter()
        .enumerate()
        .map(|(i, row)| (first_index + i, row))
        .filter(|(n, _)| *n >= from)
        .find(|(n, row)| !operator.apply_row(row, *n).iter().all(Zero::is_zero))
    {
        return exact_hypothesis(
            name,
            statement,
            Err((Some(n), format!("H does not kill P_{n}"))),
        );
    }
    let decomposition = match find_invariant_annihilator_system(
        operator,
        prepared,
        from,
        MAX_HADAMARD_OPERATOR_TERMS,
    ) {
        Ok(Some(decomposition)) => decomposition,
        Ok(None) => {
            return unresolved_hypothesis(
                name,
                statement,
                "no invariant row-annihilator system of the rule containing H was found".into(),
            )
        }
        Err(error) => {
            return unresolved_hypothesis(
                name,
                statement,
                format!("the annihilator search stopped: {error}"),
            )
        }
    };
    if decomposition.multiplier.is_zero() {
        return unresolved_hypothesis(name, statement, "the multiplier mu is zero".into());
    }
    let restart = decomposition.system.restart_bases(from);
    if let Err(error) =
        check_system_rows(&decomposition.system.generators, prepared, from, &restart)
    {
        return unresolved_hypothesis(name, statement, error.to_string());
    }
    let exceptional = vanishing_bases(&decomposition.multiplier, from);
    let last = exceptional.iter().copied().max().unwrap_or(from);
    if last >= first_index + rows.len() {
        return unresolved_hypothesis(
            name,
            statement,
            format!("the exceptional base {last} lies beyond the generated rows"),
        );
    }
    if let Some(n) = exceptional.iter().find(|n| {
        !operator
            .apply_row(&rows[*n - first_index], **n)
            .iter()
            .all(Zero::is_zero)
    }) {
        return exact_hypothesis(
            name,
            statement,
            Err((Some(*n), format!("H does not kill P_{n} at a zero of mu"))),
        );
    }
    exact_hypothesis(
        name,
        statement,
        Ok(format!(
            "invariant system with {} generator(s) from base {from} (restart bases {:?}), \
             mu H = sum_l Z_l A^(l), H P_m = 0 checked exactly at the zeros {:?} of mu",
            decomposition.system.generators.len(),
            restart,
            exceptional
        )),
    )
}

/// `U_0(n) != 0` with one sign for every `n >= start - max_lag` from the
/// shared leading-coefficient module: the scalar recurrence
/// `f(n) U_0(n) = sum_r sigma_r(n) U_0(n - r)` (valid for `n >= start`),
/// cleared of the denominators `L(n)` of the `sigma_r`, is a recurrence of
/// constant rows `f(n) L(n) U_0(n) = sum_r (L sigma_r)(n) U_0(n - r)` with
/// the seeds `U_0(start - max_lag), ..., U_0(start - 1)`, and a proved
/// certificate with degree zero and one sign from the first seed on is the
/// claim (the module also checks `f(n) L(n) != 0`).  `None` when it is not
/// proved.
pub(super) fn shared_first_coefficient_hypothesis(
    lag_coefficients: &BTreeMap<usize, IndexFunction>,
    left: &[Q],
    values: &[(usize, Q)],
    start: usize,
    max_lag: usize,
) -> Option<LiftHypothesis> {
    let first = start.checked_sub(max_lag)?;
    let seeds: Vec<Vec<Q>> = values
        .iter()
        .filter(|(n, _)| *n < start)
        .map(|(_, value)| vec![value.clone()])
        .collect();
    if seeds.len() != max_lag || lag_coefficients.is_empty() {
        return None;
    }
    let common = lag_coefficients
        .values()
        .fold(vec![Q::one()], |acc, c| vec_poly::mul(&acc, &c.denominator));
    let constant_in_x = |p: Vec<Q>| {
        BivarPoly::from_coefficients_normalized(p.into_iter().map(|c| vec![c]).collect())
    };
    let terms = lag_coefficients
        .iter()
        .map(|(&lag, c)| {
            let others = lag_coefficients
                .iter()
                .filter(|(other, _)| **other != lag)
                .fold(vec![Q::one()], |acc, (_, o)| {
                    vec_poly::mul(&acc, &o.denominator)
                });
            RecurrenceTerm::new(lag, 0, constant_in_x(vec_poly::mul(&c.numerator, &others)))
        })
        .collect();
    let scalar = Recurrence {
        terms,
        denominator: Some(constant_in_x(vec_poly::mul(left, &common))),
        inhomogeneous: None,
    };
    shared_degree_law_hypothesis(
        "first_coefficient",
        &format!(
            "the first coefficient U_0(n) is nonzero with one sign for every n >= {first} (its \
             scalar recurrence, as a family of constant rows, by the shared leading-coefficient \
             module)"
        ),
        (&scalar, first, &seeds),
        (first, Some(0)),
        |step| step == 0,
        |_| false,
    )
}

/// `U_0(n) != 0` with one sign for every `n >= start - max_lag`: a rational
/// `rho = B / A` (degrees at most [`MAX_FIRST_COEFFICIENT_DEGREE`]) with
/// `A(m) U_0(m) = B(m) U_0(m - 1)` on every computed pair (an exact fit),
/// the identity `sum_r sigma_r(n) prod_(i < r) A(n - i) / B(n - i) = f(n)`
/// of rational functions (so the hypergeometric sequence `H(m) = rho(m)
/// H(m - 1)` started at `U_0(start - max_lag)` solves the scalar recurrence
/// of `U_0`, whose leader `f(n)` does not vanish, and agrees with `U_0` on
/// the seeds, hence everywhere), `rho(m) > 0` for `m > start - max_lag`
/// and `U_0(start - max_lag) != 0`.
pub(super) fn first_coefficient_hypothesis(
    lag_coefficients: &BTreeMap<usize, IndexFunction>,
    left: &IndexFunction,
    values: &[(usize, Q)],
    start: usize,
    max_lag: usize,
) -> LiftHypothesis {
    let name = "first_coefficient";
    let statement = "the first coefficient U_0(n) is nonzero with one sign for every n from the \
                     induction seeds on (a hypergeometric solution U_0(n) = rho(n) U_0(n - 1) of \
                     its scalar recurrence with rho(n) > 0)";
    let fail = |evidence: String| exact_hypothesis(name, statement, Err((None, evidence)));
    if lag_coefficients.is_empty() {
        return fail("every lag coefficient of the scalar recurrence of U_0 vanishes".into());
    }
    let Some((first_seed, first_value)) = values.first() else {
        return fail("no computed first coefficients".into());
    };
    if first_value.is_zero() {
        return fail(format!("U_0({first_seed}) = 0"));
    }
    let pairs: Vec<(Q, &Q, &Q)> = values
        .windows(2)
        .map(|w| (Q::from_integer(BigInt::from(w[1].0)), &w[1].1, &w[0].1))
        .collect();
    let mut rho = None;
    for degree in 0..=MAX_FIRST_COEFFICIENT_DEGREE {
        let width = 2 * (degree + 1);
        if pairs.len() < width + 3 {
            break;
        }
        // A(m) U_0(m) - B(m) U_0(m - 1) = 0, unknowns (A_0..A_d, B_0..B_d).
        let equations: Vec<Vec<Q>> = pairs
            .iter()
            .map(|(m, current, previous)| {
                let mut equation = Vec::with_capacity(width);
                for i in 0..=degree {
                    equation.push(*current * m.pow(i as i32));
                }
                for i in 0..=degree {
                    equation.push(-(*previous * m.pow(i as i32)));
                }
                equation
            })
            .collect();
        let integral: Vec<Vec<BigInt>> = equations
            .iter()
            .map(|equation| vec_poly::clear_denominators(equation))
            .collect();
        let Some(vector) = exact_kernel_vectors_by_group(&integral, &vec![0; width], 1)
            .into_iter()
            .next()
            .flatten()
        else {
            continue;
        };
        let solves = equations.iter().all(|equation| {
            equation
                .iter()
                .zip(&vector)
                .fold(Q::zero(), |acc, (a, b)| acc + a * b)
                .is_zero()
        });
        let a = vec_poly::trim(&vector[..=degree]);
        let b = vec_poly::trim(&vector[degree + 1..]);
        if solves && !vec_poly::is_zero(&a) && !vec_poly::is_zero(&b) {
            rho = Some(IndexFunction::new(b, a));
            break;
        }
    }
    // Without a hypergeometric solution: positive lag coefficients and
    // seeds of one sign keep U_0 of that sign.
    let positivity = || -> Option<String> {
        let seeds: Vec<&Q> = values
            .iter()
            .filter(|(n, _)| *n < start)
            .map(|(_, value)| value)
            .collect();
        let one_sign = seeds.len() == max_lag
            && (seeds.iter().all(|v| v.is_positive()) || seeds.iter().all(|v| v.is_negative()));
        let positive = lag_coefficients.values().all(|coefficient| {
            matches!(
                coefficient
                    .div(left)
                    .decide(start, IndexSignCondition::Positive),
                IndexSignDecision::Holds { .. }
            )
        });
        (one_sign && positive).then(|| {
            format!(
                "every sigma_r(n) / f(n) > 0 for n >= {start} (exact sign checker: {}) and the \
                 seeds U_0({}), ..., U_0({}) are nonzero of one sign",
                lag_coefficients
                    .iter()
                    .map(|(lag, c)| format!("sigma_{lag}/f = {}", c.div(left).text()))
                    .collect::<Vec<_>>()
                    .join(", "),
                start - max_lag,
                start - 1
            )
        })
    };
    let Some(rho) = rho else {
        if let Some(evidence) = positivity() {
            return exact_hypothesis(name, statement, Ok(evidence));
        }
        return unresolved_hypothesis(
            name,
            statement,
            format!(
                "no ratio U_0(n) / U_0(n - 1) of degree at most {MAX_FIRST_COEFFICIENT_DEGREE} \
                 fits the computed first coefficients, and the lag coefficients are not all \
                 positive"
            ),
        );
    };
    // sum_r sigma_r(n) prod_(i < r) 1 / rho(n - i) - f(n) = 0.
    let mut sum = left.negated();
    for (lag, coefficient) in lag_coefficients {
        let mut term = coefficient.clone();
        for i in 0..*lag {
            term = term.div(&rho.shifted(-(i as i64)));
        }
        sum = sum.add(&term);
    }
    if !vec_poly::is_zero(&sum.numerator) {
        if let Some(evidence) = positivity() {
            return exact_hypothesis(name, statement, Ok(evidence));
        }
        return unresolved_hypothesis(
            name,
            statement,
            format!(
                "the fitted ratio rho(n) = {} does not solve the scalar recurrence of U_0",
                rho.text()
            ),
        );
    }
    let from = start - max_lag + 1;
    let mut hypothesis = sign_hypothesis(
        name,
        statement,
        "rho(n)",
        &rho,
        from,
        IndexSignCondition::Positive,
    );
    hypothesis.evidence = format!(
        "U_0(n) = rho(n) U_0(n - 1) from U_0({first_seed}) = {first_value}, the identity \
         sum_r sigma_r(n) prod_(i < r) 1/rho(n - i) = f(n) holds exactly; {}",
        hypothesis.evidence
    );
    hypothesis
}
