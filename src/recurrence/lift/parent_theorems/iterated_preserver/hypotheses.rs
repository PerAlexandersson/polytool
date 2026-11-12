//! Exact hypotheses supporting the enclosing theorem.

use super::*;

/// `Lambda(n, deg P_(n-1)) != 0` for every `n >= N` (module documentation,
/// hypothesis 3).
pub(super) fn degree_step_hypothesis(
    terms: &StepCoefficients,
    n_first: usize,
    start: usize,
    start_degree: usize,
) -> LiftHypothesis {
    let name = "degree_step";
    let statement = format!(
        "deg P_n = deg P_(n-1) + delta for every n >= {n_first}: the coefficient Lambda(n, e) of \
         x^(e + delta) in f(n) T_n(x^e) is nonzero at e = deg P_(n-1)"
    );
    let Some(delta) = degree_step(terms) else {
        return unresolved_hypothesis(name, &statement, "the step operator is zero".into());
    };
    if delta < 0 {
        return unresolved_hypothesis(
            name,
            &statement,
            format!("the step lowers the degree (delta = {delta}), so rows eventually vanish"),
        );
    }
    // e(n) = deg P_(n-1) = start_degree + delta (n - 1 - start).
    let previous_degree = previous_degree_polynomial(start, start_degree, delta);
    let lambda = leading_lambda(terms, delta, &[Q::zero(), Q::one()], &previous_degree);
    let mut hypothesis = sign_hypothesis(
        name,
        &statement,
        "Lambda(n, deg P_(n-1))",
        &IndexFunction::polynomial(lambda),
        n_first,
        IndexSignCondition::Nonzero,
    );
    hypothesis.evidence = format!(
        "delta = {delta}, deg P_(n-1) = {}: {}",
        format_index_polynomial(&previous_degree, "n"),
        hypothesis.evidence
    );
    hypothesis
}

/// `Lambda(n(t), e(t))` as a polynomial in `t`, with `n(t)` and `e(t)`
/// polynomials in `t`: the coefficient of `x^(e + delta)` in
/// `f(n) T_n(x^e)`.
pub(super) fn leading_lambda(
    terms: &StepCoefficients,
    delta: i64,
    index: &[Q],
    degree: &[Q],
) -> Vec<Q> {
    let mut lambda = vec![Q::zero()];
    for (d, c) in terms {
        let Some(x_power) = x_degree(c) else {
            continue;
        };
        if x_power as i64 - *d as i64 != delta {
            continue;
        }
        let mut product = vec_poly::compose(&x_column(c, x_power), index);
        for t in 0..*d {
            let factor = vec_poly::sub(degree, &[Q::from_integer(BigInt::from(t))]);
            product = vec_poly::mul(&product, &factor);
        }
        lambda = vec_poly::add(&lambda, &product);
    }
    lambda
}

/// The degree step of the residue class of the start row `class_start`
/// (lag `L >= 2`): `Lambda(s + L t, deg P_s + delta (t - 1)) != 0` for every
/// `t >= 1`.
pub(super) fn class_degree_step_hypothesis(
    terms: &StepCoefficients,
    lag: usize,
    class_start: usize,
    class_degree: usize,
    class: usize,
) -> LiftHypothesis {
    let name = format!("degree_step_class_{class}");
    let statement = format!(
        "deg P_n = deg P_(n-{lag}) + delta for every n = {class_start} + {lag} t, t >= 1: the \
         coefficient Lambda(n, e) of x^(e + delta) in f(n) T_n(x^e) is nonzero at e = deg \
         P_(n-{lag}) = deg P_{class_start} + delta (t - 1)"
    );
    let Some(delta) = degree_step(terms) else {
        return unresolved_hypothesis(&name, &statement, "the step operator is zero".into());
    };
    if delta < 0 {
        return unresolved_hypothesis(
            &name,
            &statement,
            format!("the step lowers the degree (delta = {delta}), so rows eventually vanish"),
        );
    }
    let to_q = |value: i64| Q::from_integer(BigInt::from(value));
    let index = [to_q(class_start as i64), to_q(lag as i64)];
    let degree = vec_poly::trim(&[to_q(class_degree as i64 - delta), to_q(delta)]);
    let lambda = leading_lambda(terms, delta, &index, &degree);
    let mut hypothesis = sign_hypothesis(
        &name,
        &statement,
        "Lambda(n(t), e(t)) in t",
        &IndexFunction::polynomial(lambda),
        1,
        IndexSignCondition::Nonzero,
    );
    // The sign checker runs in t: report a failure at the row index
    // n = s + L t.
    if let LiftHypothesisStatus::Failed {
        witness: Some(witness),
    } = &mut hypothesis.status
    {
        *witness = BigInt::from(class_start) + BigInt::from(lag) * &*witness;
    }
    hypothesis.evidence = format!(
        "delta = {delta}, n = {class_start} + {lag} t, deg P_(n-{lag}) = {}: {}",
        format_index_polynomial(&degree, "t"),
        hypothesis.evidence.replace(" n = ", " t = ")
    );
    hypothesis
}

/// `e(n) = deg P_(n-1) = start_degree + delta (n - 1 - start)`.
fn previous_degree_polynomial(start: usize, start_degree: usize, delta: i64) -> Vec<Q> {
    let delta_q = Q::from_integer(BigInt::from(delta));
    vec_poly::trim(&[
        Q::from_integer(BigInt::from(start_degree))
            - &delta_q * Q::from_integer(BigInt::from(start + 1)),
        delta_q,
    ])
}

/// An upper bound of `deg P_(n-L)` valid in every residue class:
/// `max_c (D_c - delta (s_c + L)/L) + delta n / L` (class `c` starting at
/// `s_c = N - L + c` with degree `D_c`); for `L = 1` this is the exact
/// [`previous_degree_polynomial`].
pub(super) fn previous_degree_bound(
    n_first: usize,
    lag: usize,
    class_degrees: &[usize],
    delta: i64,
) -> Vec<Q> {
    if lag == 1 {
        return previous_degree_polynomial(n_first - 1, class_degrees[0], delta);
    }
    let lag_q = Q::from_integer(BigInt::from(lag));
    let slope = Q::from_integer(BigInt::from(delta)) / &lag_q;
    let constant = class_degrees
        .iter()
        .enumerate()
        .map(|(c, &degree)| {
            Q::from_integer(BigInt::from(degree))
                - &slope * (Q::from_integer(BigInt::from(n_first - lag + c)) + &lag_q)
        })
        .max()
        .unwrap_or_else(Q::zero);
    vec_poly::trim(&[constant, slope])
}

/// The row hypotheses: the class starts in the working region `region`
/// (in the coordinates) and the earlier seeds in the certified region
/// `claimed` (in `P` coordinates).
pub(super) fn row_hypotheses(
    family: &Family<'_>,
    region: ClosedRegion,
    claimed: ClosedRegion,
) -> Vec<LiftHypothesis> {
    let mut out = Vec::new();
    let start = family.start;
    let letter = family.row_letter();
    if family.lag == 1 {
        let seed = &family.start_rows[0];
        out.push(exact_hypothesis(
            "start_row_in_region",
            &format!(
                "{letter}_{start} is nonzero and real-rooted with every root in {}",
                region.text()
            ),
            if region.contains_roots_of(seed) {
                Ok(format!(
                    "{letter}_{start} has degree {}; checked exactly",
                    family.start_degree
                ))
            } else {
                Err((
                    Some(start),
                    format!(
                        "{letter}_{start} is zero, not real-rooted, or has a root outside {} \
                         (exact)",
                        region.text()
                    ),
                ))
            },
        ));
    } else {
        let first_start = family.n_first - family.lag;
        let failing = family
            .start_rows
            .iter()
            .position(|row| !region.contains_roots_of(row));
        out.push(exact_hypothesis(
            "start_rows_in_region",
            &format!(
                "the class starts {letter}_{first_start}, ..., {letter}_{start} are nonzero and \
                 real-rooted with every root in {}",
                region.text()
            ),
            match failing {
                None => Ok(format!(
                    "checked exactly for {letter}_{first_start}, ..., {letter}_{start}"
                )),
                Some(offset) => Err((
                    Some(first_start + offset),
                    format!(
                        "{letter}_{} is zero, not real-rooted, or has a root outside {} (exact)",
                        first_start + offset,
                        region.text()
                    ),
                )),
            },
        ));
    }
    let earlier = &family.seeds[..family.seeds.len() - family.lag];
    if !earlier.is_empty() {
        let first = family.first_index;
        let last = first + earlier.len() - 1;
        let failing = earlier
            .iter()
            .position(|row| !is_zero_row(row) && !claimed.contains_roots_of(row));
        let zero_rows: Vec<String> = earlier
            .iter()
            .enumerate()
            .filter(|(_, row)| is_zero_row(row))
            .map(|(offset, _)| format!("P_{}", first + offset))
            .collect();
        let zero_note = if zero_rows.is_empty() {
            String::new()
        } else {
            format!(
                "; {} zero (the conclusion starts after the last zero seed row)",
                zero_rows.join(", ")
            )
        };
        out.push(exact_hypothesis(
            "earlier_seeds_in_region",
            &format!(
                "every nonzero seed row among P_{first}, ..., P_{last} is real-rooted with every \
                 root in {}",
                claimed.text()
            ),
            match failing {
                None => Ok(format!(
                    "checked exactly for P_{first}, ..., P_{last}{zero_note}"
                )),
                Some(offset) => Err((
                    Some(first + offset),
                    format!(
                        "P_{} is not real-rooted, or has a root outside {} (exact)",
                        first + offset,
                        claimed.text()
                    ),
                )),
            },
        ));
    }
    out
}

/// The multiplier hypotheses: `a_n = g h` real-rooted with every root in
/// the region for every `n >= N`.
pub(super) fn multiplier_hypotheses(
    decomposition: &Decomposition,
    region: ClosedRegion,
    n_first: usize,
) -> Vec<LiftHypothesis> {
    let mut out = Vec::new();
    if !decomposition.factors.is_empty() {
        let denominators: BTreeSet<_> = decomposition.index_factor.iter().map(|(_, q)| q).collect();
        for (i, q) in denominators.into_iter().enumerate() {
            out.push(sign_hypothesis(
                &format!("multiplier_coefficients_defined_{i}"),
                "the common multiplier coefficients are defined for every rule index",
                "multiplier coefficient denominator",
                &IndexFunction::polynomial(q.clone()),
                n_first,
                IndexSignCondition::Nonzero,
            ));
        }
    }
    let g = &decomposition.fixed_factor;
    if vec_poly::degree(g).is_some_and(|d| d >= 1) {
        let text = format_index_polynomial(g, "x");
        out.push(exact_hypothesis(
            "multiplier_fixed_factor_in_region",
            &format!(
                "the multiplier factor g(x) = {text} is real-rooted with every root in {}",
                region.text()
            ),
            if region.contains_roots_of(g) {
                Ok("exact check of the fixed polynomial".into())
            } else {
                Err((
                    None,
                    format!(
                        "g(x) = {text} is not real-rooted with every root in {} (exact)",
                        region.text()
                    ),
                ))
            },
        ));
    }
    if decomposition.index_factor.len() >= 2 {
        let (prefix, reflect, nonpositive) = match region {
            ClosedRegion::NonpositiveRay => ("multiplier_index_factor", false, true),
            ClosedRegion::NonnegativeRay => ("multiplier_index_factor_reflected", true, true),
            ClosedRegion::RealLine => ("multiplier_index_factor", false, false),
        };
        let coefficients = decomposition.index_factor_functions(reflect);
        let (hypotheses, _) =
            polynomial_hypotheses(prefix, "x", &coefficients, n_first, nonpositive);
        out.extend(hypotheses);
    }
    out
}

/// The hypotheses of `S_n` for the region, with the theorem text and the
/// Lean declarations.
pub(super) fn preserver_hypotheses(
    operator: &PreserverOperator,
    region: ClosedRegion,
    n_first: usize,
    previous_degree: Option<&[Q]>,
) -> (Vec<LiftHypothesis>, String, Vec<String>) {
    if let Some((parts, denominator)) = weyl_step_of(operator) {
        return weyl_step_preserver_hypotheses(&parts, &denominator, region, n_first);
    }
    let transfer = operator_transfer(operator, n_first);
    let mut out = transfer.hypotheses.clone();
    let mut theorem = transfer.theorem.clone();
    let mut lean_declarations = transfer.lean_declarations.clone();
    if let (PreserverOperator::EulerStep { raising, .. }, Some(previous_degree)) =
        (operator, previous_degree)
    {
        // The degree-local form: the step acts on P_(n-1) of degree
        // e(n) = deg P_(n-1), so u + v k > 0 is needed for k <= e(n) only
        // (linear in k: at k = 0, which is u > 0, and at k = e(n)); the
        // menu's v >= 0 is dropped.
        if let Some(hypothesis) = degree_local_raising_hypothesis(raising, previous_degree, n_first)
        {
            out.retain(|h| h.name != "euler_step_raising_slope_nonnegative");
            out.push(hypothesis);
            theorem = EULER_STEP_DEGREE_LOCAL.into();
            lean_declarations.push("RealRooted.EulerBidiagonal.generalRowsDep_spec".into());
        }
    }
    if let Some((_, denominator)) = operator.weyl_parts() {
        if vec_poly::degree(&denominator).is_some_and(|d| d >= 1) {
            out.push(sign_hypothesis(
                "preserver_coefficients_defined",
                &format!(
                    "the coefficients of S_n are finite at every n >= {n_first} (their common \
                     denominator does not vanish)"
                ),
                "delta_S(n)",
                &IndexFunction::polynomial(denominator),
                n_first,
                IndexSignCondition::Nonzero,
            ));
        }
    }
    let ray_image_in_ray = transfer
        .from_nonpositive
        .as_ref()
        .is_some_and(|image| region_contains(&nonpositive_ray(), image));
    let reflection_invariant = matches!(
        operator,
        PreserverOperator::EulerPolynomial { .. }
            | PreserverOperator::PowerTimesDerivative { .. }
            | PreserverOperator::DerivativeOfPowerMultiple { .. }
    );
    let (holds, how) = match region {
        ClosedRegion::NonpositiveRay => (
            ray_image_in_ray,
            "the menu image of (-inf, 0] lies in (-inf, 0]".to_string(),
        ),
        ClosedRegion::NonnegativeRay => (
            ray_image_in_ray && reflection_invariant,
            "S_n commutes up to sign with x -> -x (a polynomial in theta, x^m D^j or D^j x^m) \
             and the menu image of (-inf, 0] lies in (-inf, 0], so [0, inf) maps to [0, inf)"
                .to_string(),
        ),
        ClosedRegion::RealLine => (
            transfer.from_real_line,
            "the menu theorem holds for real-rooted polynomials anywhere on the line".to_string(),
        ),
    };
    out.push(exact_hypothesis(
        "preserver_keeps_region",
        &format!(
            "S_n maps every polynomial that is zero or real-rooted with every root in {} to zero or \
             such a polynomial, for every n >= {n_first}",
            region.text()
        ),
        if holds {
            Ok(how)
        } else {
            Err((None, format!("not given by the menu theorem: {how} fails")))
        },
    ));
    (out, theorem, lean_declarations)
}

/// The coefficients `(A, B, C)` and the denominator of a Weyl-step operator
/// `S_n = (A + B D + C D^2) / delta_S(n)`.
pub(super) fn weyl_step_of(operator: &PreserverOperator) -> Option<([Bivariate; 3], Vec<Q>)> {
    if !matches!(operator, PreserverOperator::WeylOperator { .. }) {
        return None;
    }
    let (numerator, denominator) = operator.weyl_parts()?;
    Some((weyl_step_parts(&numerator)?, denominator))
}

/// The hypotheses of Brändén's Weyl step for the region (module
/// documentation, "Brändén's Weyl step").
fn weyl_step_preserver_hypotheses(
    parts: &[Bivariate; 3],
    denominator: &[Q],
    region: ClosedRegion,
    n_first: usize,
) -> (Vec<LiftHypothesis>, String, Vec<String>) {
    let mut out = weyl_step_hypotheses(parts, n_first);
    if vec_poly::degree(denominator).is_some_and(|d| d >= 1) {
        out.push(sign_hypothesis(
            "preserver_coefficients_defined",
            &format!(
                "the coefficients of S_n are finite at every n >= {n_first} (their common \
                 denominator does not vanish)"
            ),
            "delta_S(n)",
            &IndexFunction::polynomial(denominator.to_vec()),
            n_first,
            IndexSignCondition::Nonzero,
        ));
    }
    let how = match region {
        ClosedRegion::RealLine => {
            "Branden's Corollary 3.5 holds for real-rooted polynomials anywhere on the line"
        }
        ClosedRegion::NonpositiveRay => {
            out.extend(weyl_step_region_hypotheses(parts, false, n_first));
            "with coefficients of one sign, S_n maps a polynomial with coefficients of one sign \
             to one with coefficients of one sign, real-rooted by Corollary 3.5, so with every \
             root in (-inf, 0]"
        }
        ClosedRegion::NonnegativeRay => {
            out.extend(weyl_step_region_hypotheses(parts, true, n_first));
            "the reflected operator (x -> -x) has coefficients of one sign, so S_n keeps [0, inf) \
             by Corollary 3.5 and the reflection"
        }
    };
    out.push(exact_hypothesis(
        "preserver_keeps_region",
        &format!(
            "S_n maps every polynomial that is zero or real-rooted with every root in {} to zero or \
             such a polynomial, for every n >= {n_first}",
            region.text()
        ),
        Ok(how.to_string()),
    ));
    (out, BRANDEN_WEYL_STEP.to_string(), Vec::new())
}

const EULER_STEP_DEGREE_LOCAL: &str = "The degree-local Euler step: T = kappa (theta + a)(theta + \
     b) + x (u + v theta), theta = x D, with kappa, a, b > 0, u + v k > 0 for every k <= d and a \
     negative comparison defect kappa (a b - (a + b + 1)^2/4) + (u - (a + b + 1) v/2) s - v^2 \
     s^2/(4 kappa) < 0 for every s <= 0, maps every polynomial of degree <= d with nonnegative \
     coefficients and only real roots to such a polynomial (RealRooted \
     IsPFPolynomial.generalStep_of_natDegree_le, a half-line theorem: the symbol of T has \
     discriminant -4 kappa s^2 times the comparison defect at s, so its real-rootedness is \
     required for s <= 0 only, and the hypothesis (i) of P. Branden, Linear transformations \
     preserving the Polya frequency property, Trans. Amer. Math. Soc. 358 (2006), Corollary \
     3.5, for every real s, is not needed); here d = deg P_(n-1), so u(n) + v(n) k > 0 is \
     required for 0 <= k <= deg P_(n-1) only (both ends, the condition being linear in k); \
     by linearity a row with nonpositive coefficients is handled through its negative";

/// `u(n) + v(n) e(n) > 0` for every `n >= n_first`, with `e(n) = deg
/// P_(n-1)` (with `u > 0` from the menu, `u + v k > 0` for `0 <= k <= e(n)`).
fn degree_local_raising_hypothesis(
    raising: &[IndexRationalFunction],
    previous_degree: &[Q],
    n_first: usize,
) -> Option<LiftHypothesis> {
    let [u, v] = raising else {
        return None;
    };
    let ((pu, qu), (pv, qv)) = (u.parse().ok()?, v.parse().ok()?);
    let numerator = vec_poly::add(
        &vec_poly::mul(&pu, &qv),
        &vec_poly::mul(&vec_poly::mul(&pv, previous_degree), &qu),
    );
    let function = IndexFunction::new(numerator, vec_poly::mul(&qu, &qv));
    Some(sign_hypothesis(
        "euler_step_raising_positive_up_to_degree",
        &format!(
            "u(n) + v(n) deg P_(n-1) > 0 for every n >= {n_first} (with u(n) > 0: u + v k > 0 \
             for 0 <= k <= deg P_(n-1), the degree-local hypothesis)"
        ),
        "u(n) + v(n) deg P_(n-1)",
        &function,
        n_first,
        IndexSignCondition::Positive,
    ))
}
