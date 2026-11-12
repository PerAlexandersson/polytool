//! Root-window candidates and exact window hypotheses.

use super::*;

// ---------------------------------------------------------------------------
// Root windows
// ---------------------------------------------------------------------------

/// A closed root window `W` with at least one finite end.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct RootWindow {
    pub(super) lower: Option<Q>,
    pub(super) upper: Option<Q>,
}

impl RootWindow {
    pub(super) fn interval(&self) -> RealInterval {
        RealInterval {
            lower: self.lower.clone(),
            lower_closed: self.lower.is_some(),
            upper: self.upper.clone(),
            upper_closed: self.upper.is_some(),
        }
    }

    pub(super) fn text(&self) -> String {
        RootInterval::from_real(&self.interval()).to_string()
    }
}

/// The candidate windows: the rays and gaps of the rational roots of `w` for
/// a separable `c_2 = gamma(n) w(x)`.
fn window_candidates(c2: &Bivariate) -> Vec<RootWindow> {
    let Some((_, w)) = separate(c2) else {
        return Vec::new();
    };
    let Ok(factorization) = crate::factorization::irreducible_factorization_rational_coeffs(&w)
    else {
        return Vec::new();
    };
    let mut ends: Vec<Q> = factorization
        .factors
        .iter()
        .filter(|(factor, _)| factor.len() == 2)
        .map(|(factor, _)| -Q::new(factor[0].clone(), factor[1].clone()))
        .collect();
    ends.sort();
    ends.dedup();
    let mut windows = Vec::new();
    for end in &ends {
        windows.push(RootWindow {
            lower: None,
            upper: Some(end.clone()),
        });
        windows.push(RootWindow {
            lower: Some(end.clone()),
            upper: None,
        });
    }
    for pair in ends.windows(2) {
        windows.push(RootWindow {
            lower: Some(pair[0].clone()),
            upper: Some(pair[1].clone()),
        });
    }
    windows
}

/// `c(n, a + b t)`, as a bivariate polynomial in `n` and `t`.
fn substitute_bivariate(c: &Bivariate, a: &Q, b: &Q) -> Bivariate {
    c.iter()
        .map(|row| {
            let row = vec_poly::trim(row);
            if vec_poly::is_zero(&row) {
                Vec::new()
            } else {
                let mut substituted = vec_poly::compose(&row, &[a.clone(), b.clone()]);
                vec_poly::strip_trailing_zeros(&mut substituted);
                substituted
            }
        })
        .collect()
}

/// `c(n, value) / f(n)` as a function of `n`.
fn value_at(c: &Bivariate, value: &Q, f: &[Q]) -> IndexFunction {
    let column: Vec<Q> = c
        .iter()
        .map(|row| vec_poly::evaluate(&vec_poly::trim(row), value))
        .collect();
    IndexFunction::new(column, f.to_vec())
}

/// Whether every `t`-coefficient of `c / f` satisfies `condition` for every
/// `n >= first`, with the decisions as evidence (stopping at the first
/// failure).
fn all_columns_decided(
    c: &Bivariate,
    f: &[Q],
    first: usize,
    condition: IndexSignCondition,
    label: &str,
    evidence: &mut Vec<String>,
) -> bool {
    for j in 0..x_width(c) {
        let column = x_column(c, j);
        if column.is_empty() {
            continue;
        }
        let decision = IndexFunction::new(column, f.to_vec()).decide(first, condition);
        evidence.push(decision.describe(&format!("[t^{j}] {label}")));
        if !decision.holds() {
            return false;
        }
    }
    true
}

/// One sign condition with its evidence.
fn decided(
    function: &IndexFunction,
    first: usize,
    condition: IndexSignCondition,
    label: &str,
    evidence: &mut Vec<String>,
) -> bool {
    let decision = function.decide(first, condition);
    evidence.push(decision.describe(label));
    decision.holds()
}

/// Whether the nonzero real-rooted `p` has every root in the closed window.
fn has_roots_in_window(p: &[Q], window: &RootWindow) -> bool {
    has_real_roots_in_region(
        p,
        &RealInterval {
            lower: window.lower.clone(),
            lower_closed: true,
            upper: window.upper.clone(),
            upper_closed: true,
        },
    )
}

/// The hypotheses of the window version for `window` (module documentation,
/// "Root windows").
#[allow(clippy::too_many_arguments)]
pub(super) fn window_hypotheses(
    window: &RootWindow,
    c1: &Bivariate,
    c2: &Bivariate,
    f: &[Q],
    n_first: usize,
    seeds: (&[Q], &[Q]),
    full_growth: bool,
    row: &dyn Fn(usize) -> String,
) -> Vec<LiftHypothesis> {
    let w = window.text();
    let mut hypotheses = Vec::new();
    // B_n <= 0 on W.
    let statement = format!("B_n(x) <= 0 for every x in {w} and every n >= {n_first}");
    let mut evidence = Vec::new();
    let one = Q::one();
    let below = window.upper.as_ref().map(|u| {
        let mut e = Vec::new();
        let holds = all_columns_decided(
            &substitute_bivariate(c2, u, &-one.clone()),
            f,
            n_first,
            IndexSignCondition::Nonpositive,
            &format!("B_n({} - t)", format_rational_coeff(u)),
            &mut e,
        );
        (holds, e)
    });
    let above = window.lower.as_ref().map(|l| {
        let mut e = Vec::new();
        let holds = all_columns_decided(
            &substitute_bivariate(c2, l, &one),
            f,
            n_first,
            IndexSignCondition::Nonpositive,
            &format!("B_n({} + t)", format_rational_coeff(l)),
            &mut e,
        );
        (holds, e)
    });
    let mut sign_ok = false;
    for (holds, e) in [below, above].into_iter().flatten() {
        evidence.extend(e);
        if holds {
            sign_ok = true;
            break;
        }
    }
    if !sign_ok {
        if let (Some(l), Some(u)) = (&window.lower, &window.upper) {
            // Convex B_n with B_n(L), B_n(U) <= 0.
            let mut e = Vec::new();
            sign_ok = x_degree(c2).is_some_and(|d| d <= 2)
                && decided(
                    &column(c2, 2, f),
                    n_first,
                    IndexSignCondition::Nonnegative,
                    "[x^2] B_n",
                    &mut e,
                )
                && decided(
                    &value_at(c2, l, f),
                    n_first,
                    IndexSignCondition::Nonpositive,
                    &format!("B_n({})", format_rational_coeff(l)),
                    &mut e,
                )
                && decided(
                    &value_at(c2, u, f),
                    n_first,
                    IndexSignCondition::Nonpositive,
                    &format!("B_n({})", format_rational_coeff(u)),
                    &mut e,
                );
            if sign_ok {
                evidence = vec![format!(
                    "B_n is convex in x with B_n <= 0 at both ends: {}",
                    e.join("; ")
                )];
            } else {
                evidence.extend(e);
            }
        }
    }
    hypotheses.push(if sign_ok {
        exact_hypothesis(
            "lag_coefficient_nonpositive",
            &statement,
            Ok(evidence.join("; ")),
        )
    } else {
        unresolved_hypothesis(
            "lag_coefficient_nonpositive",
            &statement,
            format!(
                "the sufficient tests fail (not a disproof): {}",
                evidence.join("; ")
            ),
        )
    });
    // Above U: A_n > 0 and B_n >= 0.
    if let Some(u) = &window.upper {
        let mut e = Vec::new();
        let u_text = format_rational_coeff(u);
        let holds = all_columns_decided(
            &substitute_bivariate(c1, u, &one),
            f,
            n_first,
            IndexSignCondition::Nonnegative,
            &format!("A_n({u_text} + t)"),
            &mut e,
        ) && decided(
            &value_at(c1, &(u + &one), f),
            n_first,
            IndexSignCondition::Positive,
            &format!("A_n({})", format_rational_coeff(&(u + &one))),
            &mut e,
        ) && all_columns_decided(
            &substitute_bivariate(c2, u, &one),
            f,
            n_first,
            IndexSignCondition::Nonnegative,
            &format!("B_n({u_text} + t)"),
            &mut e,
        );
        let statement = format!(
            "A_n(x) > 0 and B_n(x) >= 0 for every x > {u_text} and every n >= {n_first} (no row \
             gets a root > {u_text})"
        );
        hypotheses.push(if holds {
            exact_hypothesis("window_upper_end", &statement, Ok(e.join("; ")))
        } else {
            unresolved_hypothesis(
                "window_upper_end",
                &statement,
                format!("the monomial tests fail (not a disproof): {}", e.join("; ")),
            )
        });
    }
    // Below L: A_n < 0 and B_n >= 0, with full degree growth.
    if let Some(l) = &window.lower {
        let l_text = format_rational_coeff(l);
        let statement = format!(
            "deg P_n = deg P_(n-1) + 1, A_n(x) < 0 and B_n(x) >= 0 for every x < {l_text} and \
             every n >= {n_first} (no row gets a root < {l_text})"
        );
        let mut e = Vec::new();
        let holds = full_growth
            && all_columns_decided(
                &substitute_bivariate(c1, l, &-one.clone()),
                f,
                n_first,
                IndexSignCondition::Nonpositive,
                &format!("A_n({l_text} - t)"),
                &mut e,
            )
            && decided(
                &value_at(c1, &(l - &one), f),
                n_first,
                IndexSignCondition::Negative,
                &format!("A_n({})", format_rational_coeff(&(l - &one))),
                &mut e,
            )
            && all_columns_decided(
                &substitute_bivariate(c2, l, &-one.clone()),
                f,
                n_first,
                IndexSignCondition::Nonnegative,
                &format!("B_n({l_text} - t)"),
                &mut e,
            );
        if !full_growth {
            e.insert(0, "the degree law is not full growth from the seeds".into());
        }
        hypotheses.push(if holds {
            exact_hypothesis("window_lower_end", &statement, Ok(e.join("; ")))
        } else {
            unresolved_hypothesis(
                "window_lower_end",
                &statement,
                format!("the tests fail (not a disproof): {}", e.join("; ")),
            )
        });
    }
    let (previous_index, current_index) = (n_first - 2, n_first - 1);
    let inside = has_roots_in_window(seeds.0, window) && has_roots_in_window(seeds.1, window);
    hypotheses.push(exact_hypothesis(
        "seed_roots_in_window",
        &format!(
            "{} and {} have every root in {w}",
            row(previous_index),
            row(current_index)
        ),
        if inside {
            Ok("checked exactly".into())
        } else {
            Err((
                Some(current_index),
                format!("a seed of the induction has a root outside {w} (exact)"),
            ))
        },
    ));
    hypotheses
}

/// The first candidate window whose hypotheses all hold.
#[allow(clippy::too_many_arguments)]
pub(super) fn proved_window(
    c1: &Bivariate,
    c2: &Bivariate,
    f: &[Q],
    n_first: usize,
    seeds: (&[Q], &[Q]),
    full_growth: bool,
    row: &dyn Fn(usize) -> String,
) -> Option<(RootWindow, Vec<LiftHypothesis>)> {
    window_candidates(c2).into_iter().find_map(|window| {
        let hypotheses = window_hypotheses(&window, c1, c2, f, n_first, seeds, full_growth, row);
        hypotheses
            .iter()
            .all(|h| h.status == LiftHypothesisStatus::Verified)
            .then_some((window, hypotheses))
    })
}
