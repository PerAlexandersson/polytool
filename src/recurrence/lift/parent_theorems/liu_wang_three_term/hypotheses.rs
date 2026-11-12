//! Exact hypotheses supporting the enclosing theorem.

use super::*;

/// The degree-law hypotheses (module documentation, "The degree law").
pub(super) fn degree_law_hypotheses(
    alpha: usize,
    beta: usize,
    seed_degrees: (usize, usize),
    (c1, cd, c2, f): (&Bivariate, Option<&Bivariate>, &Bivariate, &[Q]),
    n_first: usize,
    rows: Option<&[Vec<Q>]>,
) -> Vec<LiftHypothesis> {
    let delta = seed_degrees.1 - seed_degrees.0;
    // (next state, the power k of x with u_n = [x^k] A_n + d [x^(k+1)] V_n
    // (d = deg P_(n-1)), the power of x of v_n = [x^j] B_n)
    let table =
        |from: usize| -> Result<(usize, Option<usize>, Option<usize>), String> {
            Ok(match (alpha, beta, from) {
                (1, 2, 0) | (0, 2, 0) => return Err(
                    "deg_x c_2 = 2 and the previous rows have equal degrees: the degree may grow \
                     by two"
                        .into(),
                ),
                (1, 2, _) => (1, Some(1), Some(2)),
                (1, 1, 0) => (1, Some(1), Some(1)),
                (1, _, _) => (1, Some(1), None),
                (0, 0, 0) => (0, Some(0), Some(0)),
                (0, 0, _) => (0, Some(0), None),
                (0, 1, 0) => (1, None, Some(1)),
                (0, 1, _) => (0, Some(0), Some(1)),
                (0, 2, _) => (1, None, Some(2)),
                _ => return Err("unsupported coefficient degrees".into()),
            })
        };
    let mut transitions: Vec<(usize, usize, Option<usize>, Option<usize>)> = Vec::new();
    let mut state = delta;
    while !transitions.iter().any(|t| t.0 == state) {
        match table(state) {
            Ok((to, k, j)) => {
                transitions.push((state, to, k, j));
                state = to;
            }
            Err(reason) => {
                return vec![unresolved_hypothesis(
                    "degree_law",
                    &format!("deg P_n - deg P_(n-1) is 0 or 1 for every n >= {n_first}"),
                    reason,
                )];
            }
        }
    }
    let single_state = transitions.len() == 1 && transitions[0].0 == transitions[0].1;
    let mut extra = Vec::new();
    let steps: Vec<DegreeStep> = transitions
        .iter()
        .map(|&(from, to, k, j)| {
            let (u, u_text) = match k {
                None => (IndexFunction::constant(Q::zero()), "0".to_string()),
                Some(k) => {
                    let derivative_column = cd.map(|c| x_column(c, k + 1)).unwrap_or_default();
                    if derivative_column.is_empty() {
                        (column(c1, k, f), format!("[x^{k}] A_n"))
                    } else {
                        // deg P_(n-1) = d_0 + to (n - N) on one state, else at
                        // least d_0 (with [x^(k+1)] V_n >= 0).
                        let d0 = Q::from_integer(BigInt::from(seed_degrees.1));
                        let degree: Vec<Q> = if single_state && to == 1 {
                            vec![d0 - Q::from_integer(BigInt::from(n_first)), Q::one()]
                        } else {
                            vec![d0]
                        };
                        if !single_state {
                            extra.push(sign_hypothesis(
                                &format!("derivative_degree_coefficient_nonnegative_{k}"),
                                &format!(
                                    "[x^{}] V_n >= 0 for every n >= {n_first} (so that \
                                     [x^{k}] A_n + deg P_(n-1) [x^{}] V_n is at least its \
                                     value at deg P_(n-1) = {})",
                                    k + 1,
                                    k + 1,
                                    seed_degrees.1
                                ),
                                &format!("[x^{}] V_n", k + 1),
                                &column(cd.expect("a derivative term"), k + 1, f),
                                n_first,
                                IndexSignCondition::Nonnegative,
                            ));
                        }
                        let numerator = vec_poly::add(
                            &vec_poly::trim(&x_column(c1, k)),
                            &vec_poly::mul(&degree, &derivative_column),
                        );
                        (
                            IndexFunction::new(numerator, f.to_vec()),
                            format!(
                                "[x^{k}] A_n + {} [x^{}] V_n",
                                if single_state && to == 1 {
                                    "deg P_(n-1)"
                                } else {
                                    "(a lower bound for deg P_(n-1))"
                                },
                                k + 1
                            ),
                        )
                    }
                }
            };
            let (v, v_text) = match j {
                None => (IndexFunction::constant(Q::zero()), "0".to_string()),
                Some(j) => (column(c2, j, f), format!("[x^{j}] B_n")),
            };
            DegreeStep {
                from,
                to,
                u,
                v,
                u_text,
                v_text,
            }
        })
        .collect();
    let statement = format!(
        "deg P_n - deg P_(n-1) is 0 or 1 and P_n, P_(n-2) have positive leading coefficients \
         (after the sign normalization) for every n >= {n_first}"
    );
    let schedule = match (steps.len(), single_state) {
        (1, true) if steps[0].to == 1 => "full growth (deg P_n = deg P_(n-1) + 1)".to_string(),
        (1, true) => "constant degree".to_string(),
        _ => format!(
            "degree steps (deg P_(n-1) - deg P_(n-2) -> deg P_n - deg P_(n-1)) {}",
            steps
                .iter()
                .map(|s| format!("{} -> {}", s.from, s.to))
                .collect::<Vec<_>>()
                .join(", ")
        ),
    };
    // The leading coefficients l_k (the coefficient of x^(d_k) for the
    // degrees d_k of the schedule) of the computed rows, k >= N - 2.
    let leads: Vec<Q> = match rows {
        Some(rows) => {
            let mut degrees = vec![seed_degrees.0, seed_degrees.1];
            let mut state = delta;
            while degrees.len() < rows.len() {
                state = steps
                    .iter()
                    .find(|s| s.from == state)
                    .map_or(state, |s| s.to);
                degrees.push(degrees[degrees.len() - 1] + state);
            }
            rows.iter()
                .zip(&degrees)
                .map(|(row, &degree)| vec_poly::coeff(row, degree))
                .collect()
        }
        None => Vec::new(),
    };
    steps
        .iter()
        .map(|step| {
            let name = if steps.len() == 1 {
                "degree_law".to_string()
            } else {
                format!("degree_law_from_step_{}", step.from)
            };
            let label = format!(
                "l_n = ({}) l_(n-1) + ({}) l_(n-2) when deg P_(n-1) - deg P_(n-2) = {} (then deg \
                 P_n - deg P_(n-1) = {})",
                step.u_text, step.v_text, step.from, step.to
            );
            let mut hypothesis = positive_sequence_hypothesis(
                &name,
                &format!("{statement}: {schedule}; {label}"),
                "l_n",
                &step.u,
                &step.v,
                n_first,
                &leads,
                single_state,
                false,
            );
            hypothesis.evidence = format!("{label}: {}", hypothesis.evidence);
            hypothesis
        })
        .chain(extra)
        .collect()
}

/// `u_k > 0` for every `k >= first - 2` for the sequence
/// `u_n = s(n) u_(n-1) + t(n) u_(n-2)` (`n >= first`), with positive
/// `u_(first-2)`, `u_(first-1)` and `values = u_(first-2), u_(first-1), ...`
/// computed exactly (possibly empty): by `s, t >= 0`, `s + t > 0`, or (when
/// `allow_ratio`) a positive-ratio induction (module documentation, "The
/// degree law").  A computed value `<= 0` refutes the statement when
/// `refutes` (the values are the sequence itself), and otherwise only ends
/// the attempt.
#[allow(clippy::too_many_arguments)]
pub(super) fn positive_sequence_hypothesis(
    name: &str,
    statement: &str,
    sequence: &str,
    s: &IndexFunction,
    t: &IndexFunction,
    first: usize,
    values: &[Q],
    allow_ratio: bool,
    refutes: bool,
) -> LiftHypothesis {
    if let Some(offset) = values.iter().position(|v| !v.is_positive()) {
        let index = first - 2 + offset;
        let evidence = format!(
            "{sequence} = {} at n = {index} (exact)",
            format_rational_coeff(&values[offset])
        );
        return if refutes {
            exact_hypothesis(name, statement, Err((Some(index), evidence)))
        } else {
            unresolved_hypothesis(
                name,
                statement,
                format!("{evidence}, so the sufficient test fails (not a disproof)"),
            )
        };
    }
    let s_decision = s.decide(first, IndexSignCondition::Nonnegative);
    let t_decision = t.decide(first, IndexSignCondition::Nonnegative);
    let sum_decision = s.add(t).decide(first, IndexSignCondition::Positive);
    if s_decision.holds() && t_decision.holds() && sum_decision.holds() {
        return exact_hypothesis(
            name,
            statement,
            Ok(format!(
                "positive seeds and nonnegative step coefficients: {}; {}; {}",
                s_decision.describe("s_n"),
                t_decision.describe("t_n"),
                sum_decision.describe("s_n + t_n")
            )),
        );
    }
    let direct = format!(
        "{}; {}; {}",
        s_decision.describe("s_n"),
        t_decision.describe("t_n"),
        sum_decision.describe("s_n + t_n")
    );
    if !allow_ratio {
        return unresolved_hypothesis(
            name,
            statement,
            format!(
                "the step coefficients are not all nonnegative (not a disproof): {direct} (s_n = \
                 {}, t_n = {})",
                s.text(),
                t.text()
            ),
        );
    }
    let ratios: Vec<(usize, Q)> = values
        .windows(2)
        .enumerate()
        .skip(1)
        .map(|(offset, pair)| (first - 1 + offset, &pair[1] / &pair[0]))
        .collect();
    let Some(first_ratio) = ratios.first().map(|(_, r)| r.clone()) else {
        return unresolved_hypothesis(
            name,
            statement,
            format!("no computed value of {sequence} at n = {first} (not a disproof): {direct}"),
        );
    };
    let t_nonpositive = t.decide(first + 1, IndexSignCondition::Nonpositive);
    if !t_nonpositive.holds() {
        return unresolved_hypothesis(
            name,
            statement,
            format!(
                "neither nonnegative step coefficients ({direct}) nor t_n <= 0 for a \
                 positive-ratio \
                 induction ({}) (not a disproof)",
                t_nonpositive.describe("t_n")
            ),
        );
    }
    let mut candidates = Vec::new();
    if let Some(fitted) = fit_index_rational_function(&ratios) {
        candidates.push(("the exact ratio", fitted));
    }
    candidates.push((
        "s_n / 2",
        s.scale(&Q::new(BigInt::from(1), BigInt::from(2))),
    ));
    let mut failures = Vec::new();
    for (label, g) in candidates {
        let positive = g.decide(first, IndexSignCondition::Positive);
        if !positive.holds() {
            failures.push(format!(
                "g(n) = {} ({label}): {}",
                g.text(),
                positive.describe("g")
            ));
            continue;
        }
        let Some(g_first) = g.at(first) else {
            failures.push(format!(
                "g(n) = {} ({label}): undefined at n = {first}",
                g.text()
            ));
            continue;
        };
        if first_ratio < g_first {
            failures.push(format!(
                "g(n) = {} ({label}): u_{first}/u_{} = {} < g({first}) = {}",
                g.text(),
                first - 1,
                format_rational_coeff(&first_ratio),
                format_rational_coeff(&g_first)
            ));
            continue;
        }
        let condition = s.add(&t.div(&g.shifted(-1))).sub(&g);
        let decision = condition.decide(first + 1, IndexSignCondition::Nonnegative);
        if !decision.holds() {
            failures.push(format!(
                "g(n) = {} ({label}): {}",
                g.text(),
                decision.describe("s_n + t_n / g(n-1) - g(n)")
            ));
            continue;
        }
        let mut hypothesis = LiftHypothesis::from_sign_decision(
            name,
            statement,
            decision,
            "s_n + t_n / g(n-1) - g(n)",
        );
        hypothesis.evidence = format!(
            "positive-ratio induction with g(n) = {} ({label}): {}; u_{first}/u_{} = {} >= \
             g({first}) = {} exactly; {}; and {}",
            g.text(),
            positive.describe("g(n)"),
            first - 1,
            format_rational_coeff(&first_ratio),
            format_rational_coeff(&g_first),
            t_nonpositive.describe("t_n"),
            hypothesis.evidence
        );
        return hypothesis;
    }
    unresolved_hypothesis(
        name,
        statement,
        format!(
            "no candidate closes the positive-ratio induction (not a disproof): {}",
            failures.join("; ")
        ),
    )
}

/// "`B_n(x) <= 0` for every real `x` and every `n >= first`" (module
/// documentation, "The sign condition").
fn everywhere_nonpositive_hypothesis(
    label: &str,
    c2: &Bivariate,
    f: &[Q],
    first: usize,
) -> LiftHypothesis {
    let name = "lag_coefficient_nonpositive";
    let statement = format!("{label}(x) <= 0 for every real x and every n >= {first}");
    match x_degree(c2) {
        Some(0) => sign_hypothesis(
            name,
            &statement,
            label,
            &column(c2, 0, f),
            first,
            IndexSignCondition::Nonpositive,
        ),
        Some(2) => {
            let (p0, p2) = (column(c2, 0, f), column(c2, 2, f));
            let discriminant = IndexFunction::new(
                vec_poly::sub(
                    &vec_poly::mul(&x_column(c2, 1), &x_column(c2, 1)),
                    &vec_poly::scale(
                        &vec_poly::mul(&x_column(c2, 2), &x_column(c2, 0)),
                        &Q::from_integer(BigInt::from(4)),
                    ),
                ),
                vec_poly::mul(f, f),
            );
            let decisions = [
                (
                    format!("p_2 = [x^2] {label}"),
                    p2.decide(first, IndexSignCondition::Nonpositive),
                ),
                (
                    format!("p_0 = [x^0] {label}"),
                    p0.decide(first, IndexSignCondition::Nonpositive),
                ),
                (
                    "p_1^2 - 4 p_2 p_0".to_string(),
                    discriminant.decide(first, IndexSignCondition::Nonpositive),
                ),
            ];
            let evidence: Vec<String> = decisions
                .iter()
                .map(|(label, decision)| decision.describe(label))
                .collect();
            if decisions.iter().all(|(_, decision)| decision.holds()) {
                exact_hypothesis(
                    name,
                    &statement,
                    Ok(format!(
                        "{label} = p_2 x^2 + p_1 x + p_0 with p_2 <= 0, p_0 <= 0 and p_1^2 <= \
                         4 p_2 \
                         p_0 (so {label} <= 0 on the real line): {}",
                        evidence.join("; ")
                    )),
                )
            } else {
                unresolved_hypothesis(
                    name,
                    &statement,
                    format!(
                        "the quadratic test p_2 <= 0, p_0 <= 0, p_1^2 <= 4 p_2 p_0 does not hold \
                         for every n (not a disproof): {}",
                        evidence.join("; ")
                    ),
                )
            }
        }
        _ => unresolved_hypothesis(
            name,
            &statement,
            format!("{label} has degree 1 in x, so it is not <= 0 on the whole line"),
        ),
    }
}

/// The sign hypothesis "`c / f <= 0` on the real line, or on `x <= 0`" for
/// the multiplier `label` (`B_n` or `V_n`), and the case that holds.
pub(super) fn multiplier_sign_hypothesis(
    name: &str,
    label: &str,
    c: &Bivariate,
    f: &[Q],
    n_first: usize,
) -> (LiftHypothesis, Option<LagSign>) {
    let mut everywhere = everywhere_nonpositive_hypothesis(label, c, f, n_first);
    if everywhere.status == LiftHypothesisStatus::Verified {
        everywhere.name = name.into();
        everywhere.statement = format!("{label}(x) <= 0 for every real x and every n >= {n_first}");
        return (everywhere, Some(LagSign::Everywhere));
    }
    let mut ray = negative_axis_sign_hypothesis(
        name,
        &format!(
            "{label}(x) <= 0 for every real x, or for every x <= 0 (with a root location), for \
             every n >= {n_first}"
        ),
        label,
        c,
        f,
        n_first,
        false,
    );
    if ray.status == LiftHypothesisStatus::Verified {
        ray.evidence = format!("on x <= 0 (x < 0 and continuity): {}", ray.evidence);
        (ray, Some(LagSign::NonpositiveAxis))
    } else {
        ray.evidence = format!(
            "neither everywhere ({}) nor on x <= 0 ({})",
            everywhere.evidence, ray.evidence
        );
        (ray, None)
    }
}

/// `factor c(n, x) / x` for a `c` divisible by `x`.
pub(super) fn scale_divided_by_x(c: &Bivariate, factor: usize) -> Bivariate {
    let factor = Q::from_integer(BigInt::from(factor));
    c.iter()
        .map(|row| row.iter().skip(1).map(|value| value * &factor).collect())
        .collect()
}

/// Whether `c_2 = -kappa(n) (x + v)^2`-type: degree two with identically
/// vanishing discriminant.
pub(super) fn is_negative_square(c2: &Bivariate) -> bool {
    x_degree(c2) == Some(2) && {
        let discriminant = vec_poly::sub(
            &vec_poly::mul(&x_column(c2, 1), &x_column(c2, 1)),
            &vec_poly::scale(
                &vec_poly::mul(&x_column(c2, 2), &x_column(c2, 0)),
                &Q::from_integer(BigInt::from(4)),
            ),
        );
        vec_poly::is_zero(&vec_poly::trim(&discriminant))
    }
}

/// The root-location hypotheses (module documentation, "Root location"):
/// the first method whose hypotheses all hold, else the hypotheses of the
/// root-sign induction (for the evidence).
#[allow(clippy::too_many_arguments)]
pub(super) fn root_location_hypotheses(
    seed_previous: &[Q],
    seed_current: &[Q],
    c1: &Bivariate,
    cd: Option<&Bivariate>,
    c2: &Bivariate,
    f: &[Q],
    n_first: usize,
    generated: Option<&[Vec<Q>]>,
    row: &dyn Fn(usize) -> String,
) -> (Vec<LiftHypothesis>, Option<RootLocation>) {
    let (previous_index, current_index) = (n_first - 2, n_first - 1);
    // Root-sign induction.
    let negative = has_real_roots_below(seed_previous, &Q::zero())
        && has_real_roots_below(seed_current, &Q::zero());
    let seed_hypothesis = exact_hypothesis(
        "seed_roots_negative",
        &format!(
            "{} and {} have every root < 0",
            row(previous_index),
            row(current_index)
        ),
        if negative {
            Ok("checked exactly".into())
        } else {
            Err((
                Some(current_index),
                "a seed of the induction has a root >= 0 (exact)".into(),
            ))
        },
    );
    let constants: Vec<Q> = match generated {
        Some(rows) => rows
            .iter()
            .map(|r| r.first().cloned().unwrap_or_else(Q::zero))
            .collect(),
        None => vec![
            seed_previous.first().cloned().unwrap_or_else(Q::zero),
            seed_current.first().cloned().unwrap_or_else(Q::zero),
        ],
    };
    let constant_hypothesis = positive_sequence_hypothesis(
        "constant_terms_positive",
        &format!(
            "P_n(0) > 0 for every n >= {n_first} (normalized rows), with P_n(0) = A_n(0) \
             P_(n-1)(0) + B_n(0) P_(n-2)(0)"
        ),
        "P_n(0)",
        &column(c1, 0, f),
        &column(c2, 0, f),
        n_first,
        &constants,
        true,
        true,
    );
    let mut root_sign = vec![seed_hypothesis, constant_hypothesis];
    // With a derivative term, P_n(0) = A_n(0) P_(n-1)(0) + B_n(0) P_(n-2)(0)
    // needs V_n(0) = 0.
    if let Some(c) = cd {
        root_sign.push(exact_hypothesis(
            "derivative_coefficient_vanishes_at_zero",
            &format!("V_n(0) = 0 for every n >= {n_first}"),
            if x_column(c, 0).is_empty() {
                Ok("[x^0] c_d is the zero polynomial".into())
            } else {
                Err((None, "[x^0] c_d is not the zero polynomial".into()))
            },
        ));
    }
    if root_sign
        .iter()
        .all(|h| h.status == LiftHypothesisStatus::Verified)
    {
        return (root_sign, Some(RootLocation::Negative));
    }
    // Nonnegative coefficients.
    let seeds_nonnegative = seed_previous
        .iter()
        .chain(seed_current)
        .all(|c| !c.is_negative());
    let mut evidence = Vec::new();
    let mut multipliers_nonnegative = true;
    let mut multipliers = vec![("A_n", c1), ("B_n", c2)];
    if let Some(c) = cd {
        multipliers.push(("V_n", c));
    }
    for (label, c) in multipliers {
        for j in 0..x_width(c) {
            if x_column(c, j).is_empty() {
                continue;
            }
            let decision = column(c, j, f).decide(n_first, IndexSignCondition::Nonnegative);
            evidence.push(decision.describe(&format!("[x^{j}] {label}")));
            multipliers_nonnegative &= decision.holds();
        }
    }
    if seeds_nonnegative && multipliers_nonnegative {
        let hypothesis = exact_hypothesis(
            "nonnegative_coefficients",
            &format!(
                "{} and {} and every coefficient [x^k] A_n, [x^k] B_n (and [x^k] V_n) for n >= \
                 {n_first} are nonnegative, so every row has nonnegative coefficients and (being \
                 real-rooted) every root <= 0",
                row(previous_index),
                row(current_index)
            ),
            Ok(format!("seeds checked exactly; {}", evidence.join("; "))),
        );
        return (
            vec![hypothesis],
            Some(RootLocation::NonnegativeCoefficients),
        );
    }
    let mut hypotheses = root_sign;
    hypotheses.push(unresolved_hypothesis(
        "nonnegative_coefficients",
        "the seeds and the multipliers A_n, B_n have nonnegative coefficients",
        format!(
            "{}{}",
            if seeds_nonnegative {
                ""
            } else {
                "a seed has a negative coefficient; "
            },
            evidence.join("; ")
        ),
    ));
    (hypotheses, None)
}

/// The requested region must contain the certified one.
pub(super) fn requested_region_hypothesis(
    requested: &RealInterval,
    certified: Option<&RealInterval>,
    zero_is_root: bool,
    current_index: usize,
) -> LiftHypothesis {
    let name = "requested_root_region";
    let Some(certified) = certified else {
        return unresolved_hypothesis(
            name,
            &format!(
                "the requested region {} contains the certified region",
                open_root_interval(requested)
            ),
            "no root region is certified (B_n <= 0 on the whole line; the root location does not \
             hold or is not shown)"
                .into(),
        );
    };
    let certified_text = RootInterval::from_real(certified).to_string();
    let statement = format!(
        "the requested region {} contains {certified_text}, the region this theorem certifies",
        open_root_interval(requested)
    );
    // The requested region is open at its finite ends.
    let upper_ok = requested
        .upper
        .as_ref()
        .is_none_or(|requested_upper| match &certified.upper {
            None => false,
            Some(upper) => {
                upper < requested_upper || (upper == requested_upper && !certified.upper_closed)
            }
        });
    let lower_ok = requested
        .lower
        .as_ref()
        .is_none_or(|requested_lower| match &certified.lower {
            None => false,
            Some(lower) => {
                lower > requested_lower || (lower == requested_lower && !certified.lower_closed)
            }
        });
    if lower_ok && upper_ok {
        exact_hypothesis(name, &statement, Ok("exact".into()))
    } else if zero_is_root && requested.upper.as_ref().is_some_and(Zero::is_zero) {
        exact_hypothesis(
            name,
            &statement,
            Err((
                Some(current_index),
                format!("0 is a root of P_{current_index}, outside the open region (exact)"),
            )),
        )
    } else {
        unresolved_hypothesis(
            name,
            &statement,
            format!("this theorem certifies the region {certified_text} only"),
        )
    }
}

/// The RealRooted route and declarations mirrored by the check.
pub(super) fn lean_correspondence(
    full_growth: bool,
    lag_sign: Option<LagSign>,
    location: Option<RootLocation>,
    window: Option<&RootWindow>,
) -> (Option<&'static str>, Vec<String>) {
    let (route, names): (Option<&str>, Vec<&str>) = if !full_growth {
        (
            None,
            vec!["RealRooted.MaWang.strictInterl_of_interlaces_evalCoeff_nonpos"],
        )
    } else if let Some(window) = window {
        match (&window.lower, &window.upper) {
            (Some(_), Some(_)) => (
                Some(LEAN_ROUTE_WINDOW),
                vec!["RealRooted.threeTerm_interlaces_of_roots_mem_Icc"],
            ),
            (None, Some(_)) => (
                Some(LEAN_ROUTE_WINDOW),
                vec!["RealRooted.threeTerm_interlaces_of_roots_le"],
            ),
            _ => (None, vec!["RealRooted.threeTerm_interlaces_of_window"]),
        }
    } else {
        match (lag_sign, location) {
            (Some(LagSign::Everywhere), _) => (
                Some(LEAN_ROUTE_SIGN),
                vec![
                    "RealRooted.threeTerm_interlaces_of_eval_nonpos",
                    "RealRooted.threeTerm_interlaces",
                    "RealRooted.eval_nonpos_of_natDegree_le_two",
                ],
            ),
            (_, Some(RootLocation::NonnegativeCoefficients)) => (
                Some(LEAN_ROUTE_SIGN),
                vec![
                    "RealRooted.threeTerm_interlaces_of_nonnegCoeffs",
                    "RealRooted.threeTerm_interlaces",
                    "RealRooted.threeTerm_hasNonnegCoeffs",
                ],
            ),
            _ => (
                Some(LEAN_ROUTE_SIGN),
                vec![
                    "RealRooted.threeTerm_interlaces_of_nonnegCoeffs",
                    "RealRooted.threeTerm_interlaces",
                ],
            ),
        }
    };
    (route, names.into_iter().map(str::to_string).collect())
}
