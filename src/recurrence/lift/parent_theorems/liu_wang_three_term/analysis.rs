//! Liu–Wang three-term analysis in the chosen coordinates.

use super::*;

/// The theorem in the given coordinates.
pub(super) fn analyze(
    mut report: RealRootednessTheoremReport,
    coordinates: &Coordinates,
    c1: &Bivariate,
    c2: &Bivariate,
    input: &AnalysisInput<'_>,
) -> RealRootednessTheoremReport {
    let AnalysisInput {
        rule,
        first_index,
        seeds,
        options,
        f,
        rows,
        derivative,
    } = *input;
    let n_first = first_index + seeds.len();
    let previous_index = n_first - 2;
    let current_index = n_first - 1;
    let coordinate_text = coordinates.text(previous_index);
    let Some((c1, c2)) = coordinates.transform_rule(c1, c2) else {
        return report.not_applicable("the coordinate change does not apply to the rule");
    };
    let (c1, c2) = (&c1, &c2);
    let Some(beta) = x_degree(c2) else {
        return report.not_applicable("the lag term of the rule is zero");
    };
    let gamma = derivative.and_then(x_degree);
    // The x-degree of the multiplier of P_(n-1), the derivative term
    // counting one less.
    let alpha = x_degree(c1)
        .unwrap_or(0)
        .max(gamma.map_or(0, |g| g.saturating_sub(1)));
    if alpha > 1 || beta > 2 || gamma.is_some_and(|g| g > 2) {
        return report.not_applicable(format!(
            "deg_x c_1 > 1, deg_x c_d > 2 or deg_x c_2 > 2{}: the degree may grow by more than \
             one",
            coordinate_text
                .as_ref()
                .map_or(String::new(), |text| format!(" (after {text})"))
        ));
    }
    let variable = if coordinate_text.is_some() { "Q" } else { "P" };
    let mut hypotheses: Vec<LiftHypothesis> = Vec::new();
    hypotheses.extend(denominator_hypothesis(f, n_first));
    let transformed = (
        coordinates.transform_row(&seeds[seeds.len() - 2], previous_index, previous_index),
        coordinates.transform_row(&seeds[seeds.len() - 1], current_index, previous_index),
    );
    let (previous, current) = match transformed {
        (Some(previous), Some(current))
            if !vec_poly::is_zero(&previous) && !vec_poly::is_zero(&current) =>
        {
            (previous, current)
        }
        (previous, _) => {
            let index = if previous.is_some_and(|p| !vec_poly::is_zero(&p)) {
                current_index
            } else {
                previous_index
            };
            hypotheses.push(exact_hypothesis(
                "seed_rows",
                &format!(
                    "{variable}_{previous_index} and {variable}_{current_index} are nonzero \
                     polynomials"
                ),
                Err((
                    Some(index),
                    format!("{variable}_{index} is zero or not a polynomial"),
                )),
            ));
            report.normal_form = format!(
                "P_n = A_n(x) P_(n-1) + B_n(x) P_(n-2) for n >= {n_first}{}",
                coordinate_text
                    .as_ref()
                    .map_or(String::new(), |text| format!("; {text}"))
            );
            report.hypotheses = hypotheses;
            return report.finish(placeholder_conclusion(first_index));
        }
    };

    // Normalizations: the power of x and the signs.  With a derivative term
    // the power of x is divided out only when x divides c_d; then
    // Q_n = P_n / x^m satisfies the rule with c_1 + m c_d / x.
    let derivative_vanishes_at_zero = derivative.is_none_or(|c| x_column(c, 0).is_empty());
    let power = if derivative_vanishes_at_zero {
        order_at_zero(&previous)
            .unwrap_or(0)
            .min(order_at_zero(&current).unwrap_or(0))
    } else {
        0
    };
    let c1 = &match derivative {
        Some(c) if power > 0 => add_bivariate(c1, &scale_divided_by_x(c, power)),
        _ => c1.clone(),
    };
    let previous_row = previous[power..].to_vec();
    let current_row = current[power..].to_vec();
    let previous_lead = previous_row.last().cloned().unwrap_or_else(Q::zero);
    let current_lead = current_row.last().cloned().unwrap_or_else(Q::zero);
    let negate_all = current_lead.is_negative();
    let alternate = previous_lead.is_negative() != current_lead.is_negative();
    // s_index (true: multiply by -1), with s_(N-1) the sign of the current
    // leading coefficient and s_n = -s_(n-1) when alternating.
    let row_sign = |index: usize| -> bool {
        let parity_flip = alternate && (index - previous_index) % 2 == 0;
        negate_all != parity_flip
    };
    let signed = |row: &[Q], index: usize| -> Vec<Q> {
        if row_sign(index) {
            vec_poly::neg(row)
        } else {
            row.to_vec()
        }
    };
    let seed_previous = signed(&previous_row, previous_index);
    let seed_current = signed(&current_row, current_index);
    let c1n = if alternate {
        negate_bivariate(c1)
    } else {
        c1.clone()
    };
    let cdn = derivative.map(|c| {
        if alternate {
            negate_bivariate(c)
        } else {
            c.clone()
        }
    });
    let cdn = cdn.as_ref();
    let row = |index: usize| {
        if power == 0 {
            format!("{variable}_{index}")
        } else {
            format!("{variable}_{index} / x^{power}")
        }
    };
    let multipliers = format!(
        "A_n(x) = {}, {}B_n(x) = {}",
        x_polynomial_text(&c1n, f),
        cdn.map_or(String::new(), |c| format!(
            "V_n(x) = {}, ",
            x_polynomial_text(c, f)
        )),
        x_polynomial_text(c2, f)
    );
    let mut notes = Vec::new();
    if let Some(text) = &coordinate_text {
        notes.push(format!("{text}, which satisfies the rule below"));
    }
    if power > 0 {
        notes.push(format!(
            "{variable}_n = x^{power} R_n for n >= {previous_index} (x^{power} divides \
             {variable}_{previous_index} and {variable}_{current_index}; R_n satisfies the same \
             rule{})",
            if derivative.is_some() {
                format!(" with A_n increased by {power} V_n / x")
            } else {
                String::new()
            }
        ));
    }
    if negate_all || alternate {
        notes.push(format!(
            "the rows are multiplied by {} so that the seeds have positive leading coefficients{}",
            match (negate_all, alternate) {
                (true, false) => "-1".to_string(),
                (_, true) => format!(
                    "{}(-1)^(n - {current_index})",
                    if negate_all { "-" } else { "" }
                ),
                _ => "1".to_string(),
            },
            if alternate {
                " (A_n and V_n change sign)"
            } else {
                ""
            }
        ));
    }
    report.normal_form = format!(
        "{variable}_n = A_n(x) {variable}_(n-1) + {}B_n(x) {variable}_(n-2) for n >= {n_first}, \
         with {multipliers}{}",
        if cdn.is_some() {
            format!("V_n(x) {variable}'_(n-1) + ")
        } else {
            String::new()
        },
        if notes.is_empty() {
            String::new()
        } else {
            format!("; {}", notes.join("; "))
        }
    );

    // The seeds.
    let previous_degree = vec_poly::degree(&seed_previous).unwrap_or(0);
    let current_degree = vec_poly::degree(&seed_current).unwrap_or(0);
    let degree_difference = current_degree.checked_sub(previous_degree);
    let seed_outcome = if !is_real_rooted_rational(&seed_previous) {
        Err((
            Some(previous_index),
            format!("{} is not real-rooted (exact)", row(previous_index)),
        ))
    } else if !is_real_rooted_rational(&seed_current) {
        Err((
            Some(current_index),
            format!("{} is not real-rooted (exact)", row(current_index)),
        ))
    } else if !matches!(degree_difference, Some(0 | 1)) {
        Err((
            Some(current_index),
            format!(
                "deg {} = {current_degree} and deg {} = {previous_degree} do not differ by 0 or 1",
                row(current_index),
                row(previous_index)
            ),
        ))
    } else if crate::check_weak_interlacing_bigint_coeffs(
        &integer_coefficients(&seed_previous),
        &integer_coefficients(&seed_current),
    ) != Some(true)
    {
        Err((
            Some(current_index),
            format!(
                "{} does not interlace {} (exact weak interlacing test)",
                row(previous_index),
                row(current_index)
            ),
        ))
    } else {
        Ok(format!(
            "degrees {previous_degree} and {current_degree}; real-rooted and {} << {} checked \
             exactly",
            row(previous_index),
            row(current_index)
        ))
    };
    let seeds_ok = seed_outcome.is_ok();
    hypotheses.push(exact_hypothesis(
        "seed_rows",
        &format!(
            "{} and {} are nonzero and real-rooted, of degrees e and e or e + 1, and {} << {} \
             (weak interlacing)",
            row(previous_index),
            row(current_index),
            row(previous_index),
            row(current_index)
        ),
        seed_outcome,
    ));

    // The normalized rows from the induction start on.
    let generated: Option<Vec<Vec<Q>>> = rows.and_then(|rows| {
        rows.iter()
            .enumerate()
            .map(|(offset, p)| {
                let index = previous_index + offset;
                let q = coordinates.transform_row(p, index, previous_index)?;
                if q.len() > power && q[..power].iter().any(|c| !c.is_zero()) {
                    return None;
                }
                let quotient = if q.len() > power {
                    q[power..].to_vec()
                } else {
                    vec![Q::zero()]
                };
                Some(signed(&quotient, index))
            })
            .collect()
    });

    // The degree law.
    let local_degree_law = match degree_difference.filter(|_| seeds_ok) {
        None => vec![unresolved_hypothesis(
            "degree_law",
            &format!(
                "deg {variable}_n - deg {variable}_(n-1) is 0 or 1 and {variable}_n, \
                 {variable}_(n-2) have leading coefficients of one sign for every n >= {n_first}"
            ),
            "not checked: the seeds do not satisfy the seed hypothesis".into(),
        )],
        Some(_) => degree_law_hypotheses(
            alpha,
            beta,
            (previous_degree, current_degree),
            (&c1n, cdn, c2, f),
            n_first,
            generated.as_deref(),
        ),
    };
    let expects_full_growth = degree_difference == Some(1) && (alpha == 1 || beta == 2);
    // In the identity coordinates, the shared leading-coefficient module
    // where the degree-step state machine is too narrow: degree steps 0 or 1
    // and positive leading coefficients after the sign normalization.
    let shared_degree_law = (seeds_ok
        && matches!(coordinates, Coordinates::Identity)
        && local_degree_law
            .iter()
            .any(|h| h.status != LiftHypothesisStatus::Verified))
    .then(|| {
        shared_degree_law_hypothesis(
            "degree_law",
            &format!(
                "deg P_n - deg P_(n-1) is 0 or 1 and P_n, P_(n-2) have positive leading \
                 coefficients (after the sign normalization) for every n >= {n_first}"
            ),
            (rule, first_index, seeds),
            (previous_index, None),
            // Full growth from the seeds (used by the root windows) needs
            // every step to be one.
            |step| step == 1 || (step == 0 && !expects_full_growth),
            row_sign,
        )
    })
    .flatten();
    match shared_degree_law {
        Some(hypothesis) => hypotheses.push(hypothesis),
        None => hypotheses.extend(local_degree_law),
    }

    // The sign condition.
    let (sign_hypothesis, lag_sign) =
        multiplier_sign_hypothesis("lag_coefficient_nonpositive", "B_n", c2, f, n_first);
    let (derivative_sign_hypothesis, lag_sign) = match cdn {
        None => (None, lag_sign),
        Some(c) => {
            let (hypothesis, sign) = multiplier_sign_hypothesis(
                "derivative_coefficient_nonpositive",
                "V_n",
                c,
                f,
                n_first,
            );
            let combined = match (lag_sign, sign) {
                (Some(LagSign::Everywhere), Some(LagSign::Everywhere)) => Some(LagSign::Everywhere),
                (Some(_), Some(_)) => Some(LagSign::NonpositiveAxis),
                _ => None,
            };
            (Some(hypothesis), combined)
        }
    };
    if cdn.is_some() {
        report.reference = format!("{REFERENCE}{THEOREM_2_3}");
    }
    if lag_sign == Some(LagSign::Everywhere) && is_negative_square(c2) {
        report.reference = format!("{REFERENCE}{PROPOSITION_3_3}");
    }
    // Root location, or a root window when neither sign case closes.
    let (location_hypotheses, location) = root_location_hypotheses(
        &seed_previous,
        &seed_current,
        &c1n,
        cdn,
        c2,
        f,
        n_first,
        generated.as_deref(),
        &row,
    );
    let full_growth = degree_difference == Some(1) && (alpha == 1 || beta == 2);
    let core = lag_sign == Some(LagSign::Everywhere) || (lag_sign.is_some() && location.is_some());
    let window = if core || cdn.is_some() {
        None
    } else {
        proved_window(
            &c1n,
            c2,
            f,
            n_first,
            (&seed_previous, &seed_current),
            full_growth,
            &row,
        )
    };
    let stripped_region = match &window {
        Some((window, window_hypotheses)) => {
            hypotheses.extend(window_hypotheses.iter().cloned());
            Some(window.interval())
        }
        None => {
            hypotheses.push(sign_hypothesis);
            hypotheses.extend(derivative_sign_hypothesis);
            match (lag_sign, location) {
                (_, Some(_)) => hypotheses.extend(location_hypotheses),
                (Some(LagSign::Everywhere), None) => {}
                (_, None) => hypotheses.push(unresolved_hypothesis(
                    "root_location",
                    &format!(
                        "every root of every row {variable}_n, n >= {previous_index}, is <= 0 \
                         (needed for B_n <= 0 on x <= 0 only), or the rows stay in a root \
                         window"
                    ),
                    location_hypotheses
                        .iter()
                        .filter(|h| h.status != LiftHypothesisStatus::Verified)
                        .map(|h| format!("{}: {}", h.name, h.evidence))
                        .collect::<Vec<_>>()
                        .join("; "),
                )),
            }
            location.map(|location| RealInterval {
                lower: None,
                lower_closed: false,
                upper: Some(Q::zero()),
                upper_closed: location == RootLocation::NonnegativeCoefficients,
            })
        }
    };
    // The root 0 of the divided power of x, then the coordinate change.
    let transformed_region = stripped_region.map(|region| {
        if power > 0 {
            hull_with_point(region, &Q::zero())
        } else {
            region
        }
    });
    let region = coordinates.region_of_rows(transformed_region);

    // Earlier rows (in the coordinates of P).
    let earlier_rows = &seeds[..seeds.len() - 2];
    let zero_rows = match &region {
        None => {
            hypotheses.extend(earlier_seeds_hypothesis(earlier_rows, first_index, None));
            Vec::new()
        }
        Some(region) => {
            let (hypothesis, zero_rows) = earlier_rows_in_region(earlier_rows, first_index, region);
            hypotheses.extend(hypothesis);
            zero_rows
        }
    };
    let identity = *coordinates == Coordinates::Identity;
    if let Some(requested) = &options.root_region {
        hypotheses.push(requested_region_hypothesis(
            requested,
            region.as_ref(),
            power > 0 && identity,
            current_index,
        ));
    }
    report.hypotheses = hypotheses;
    // No RealRooted route mirrors the derivative-lag form with weak
    // interlacing; none is recorded then.
    let (route, declarations) = if cdn.is_some() {
        (None, Vec::new())
    } else {
        lean_correspondence(
            identity && full_growth,
            lag_sign,
            location,
            window.as_ref().map(|(window, _)| window),
        )
    };
    report.lean_route = route.map(str::to_string);
    report.lean_declarations = declarations;

    let root_region = region.as_ref().map(RootInterval::from_real);
    let conclusion_first = zero_rows.last().map_or(first_index, |last| last + 1);
    let zero_text = if zero_rows.is_empty() {
        String::new()
    } else {
        let list: Vec<String> = zero_rows.iter().map(|j| format!("P_{j}")).collect();
        format!(
            "{} {} zero; ",
            list.join(", "),
            if zero_rows.len() == 1 { "is" } else { "are" }
        )
    };
    let interlacing_text = if identity {
        format!(
            "; for n >= {current_index}, P_(n-1) interlaces P_n (weakly, in the sense of Liu and \
             Wang), and deg P_n - deg P_(n-1) is 0 or 1 for n >= {n_first}"
        )
    } else {
        format!(
            "; for n >= {current_index}, Q_(n-1) interlaces Q_n (weakly, in the sense of Liu and \
             Wang; no interlacing of the P_n is claimed)"
        )
    };
    let statement = format!(
        "{zero_text}every {}P_n, n >= {first_index}, is nonzero and \
         real-rooted{}{interlacing_text}{}",
        if zero_rows.is_empty() { "" } else { "other " },
        root_region.as_ref().map_or(String::new(), |r| format!(
            " with every root in {r}{}",
            if region.as_ref().is_some_and(within_nonpositive_ray) {
                " (so its coefficients have one sign)"
            } else {
                ""
            }
        )),
        if power > 0 && identity {
            format!("; P_n = x^{power} R_n for n >= {previous_index} with R_n(0) != 0")
        } else {
            String::new()
        }
    );
    report.finish(RealRootednessTheoremConclusion {
        first_index: conclusion_first,
        zero_residue_classes: None,
        interlacing_from: identity.then_some(current_index),
        interlacing: identity.then_some(RowInterlacing::Weak),
        simple_roots: false,
        root_region,
        statement,
    })
}
