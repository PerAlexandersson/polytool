//! Seeded half-step lag windows: a reachable-orbit theorem, rather than a
//! universal companion-matrix preserver. All coefficient signs are checked
//! on the whole integer rule range; extra seeds have exact positive lag steps.

use super::*;

const REFERENCE: &str = "Seeded half-step lag window theorem: A_m=1 for m<k, \
    A_(q+k)=(1+xD)A_(q+k-1)+c_q x A_q, k>=1, c_q>0. \
    Every row is PF, degree floor(m/k), with positive leading coefficient, \
    and rows at distance at most k interlace (HalfStepRows theorem).";

/// Exact seed and lag data for the seeded half-step window theorem.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HalfStepLagCertificate {
    /// Positive lag k; degree of row at relative index m is floor(m/k).
    pub lag: usize,
    /// Positive rational lag weights for supplied seeds at relative indices k onward.
    pub seed_weights: Vec<String>,
}

/// Apply the seeded half-step lag window theorem. Coefficients of the
/// normalized predecessor are exactly 1 and xD, and the delayed coefficient
/// is c(n)x with c(n)>0. The first k seeds must be 1. No finite fit is used
/// as an all-index identity; the supplied recurrence is the defining rule.
///
/// # Errors
/// Invalid recurrence records or seed normalization.
pub fn check_half_step_lag_theorem(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    check_real_rootedness_theorem(
        RealRootednessTheorem::HalfStepLag,
        rule,
        first_index,
        seeds,
        options,
    )
}

pub(super) const ROUTE: RoutePolicy = RoutePolicy {
    reference: REFERENCE,
    phase: DiscoveryPhase::Direct,
    after_direct_eligible: None,
    replay: ParentReplayPolicy::RecordedRoute,
    discovery_replay: DiscoveryReplayPolicy::Direct,
    check_at_supplied_seeds: half_step_lag_report,
    start_policy: StartPolicy::OwnStarts,
    is_preferred: any_proof,
    certificate_name: Some("seeded half-step lag window"),
    discovery: DiscoveryPolicy::WhenApplicable,
    in_residue_class_menu: true,
    in_index_section_menu: true,
    in_degree_law_reversal_menu: true,
    parent: ParentPolicy {
        stage: ParentStage::WhileUnproved,
        region: ParentRegion::OpenNonpositiveRay,
        record_only_if_applicable: true,
    },
};

fn recognize(rule: &Recurrence) -> Result<(usize, IndexFunction), String> {
    let mut shape = recurrence_shape(rule)?;
    let form = "not a seeded half-step lag operator (1+xD) at lag one plus c(n)x at lag k>=1";
    let f = shape.denominator.clone();
    // Subtract the Euler operator before inspecting the remaining lag. This
    // also handles k=1, where the two order-zero coefficients are combined.
    for (key, power) in [((1, 0), 0), ((1, 1), 1)] {
        let mut expected = vec![vec![Q::zero(); power + 1]; f.len()];
        for (i, value) in f.iter().enumerate() {
            expected[i][power] = value.clone();
        }
        let actual = shape.terms.remove(&key).ok_or(form)?;
        let remainder = add_bivariate(
            &actual,
            &expected.iter().map(|r| vec_poly::neg(r)).collect(),
        );
        if x_degree(&remainder).is_some() {
            shape.terms.insert(key, remainder);
        }
    }
    if shape.terms.len() != 1 {
        return Err(form.into());
    }
    let (&(lag, derivative), coefficient) = shape.terms.iter().next().ok_or(form)?;
    if lag == 0
        || derivative != 0
        || x_degree(coefficient) != Some(1)
        || !x_column(coefficient, 0).is_empty()
    {
        return Err(form.into());
    }
    Ok((lag, IndexFunction::new(x_column(coefficient, 1), f)))
}

fn theta_plus_one(p: &[Q]) -> Vec<Q> {
    p.iter()
        .enumerate()
        .map(|(i, c)| c * Q::from_integer(BigInt::from(i) + 1))
        .collect()
}

fn seed_weights(seeds: &[Vec<Q>], lag: usize) -> Option<Vec<String>> {
    if seeds.len() < lag || seeds[..lag].iter().any(|p| vec_poly::trim(p) != [Q::one()]) {
        return None;
    }
    let mut weights = Vec::new();
    for m in lag..seeds.len() {
        let rest = vec_poly::sub(&seeds[m], &theta_plus_one(&seeds[m - 1]));
        let weight = rest.get(1)?.clone(); // Every delayed row has constant 1.
        if !weight.is_positive()
            || vec_poly::trim(&rest)
                != vec_poly::trim(&vec_poly::scale(&vec_poly::shift(&seeds[m - lag]), &weight))
        {
            return None;
        }
        weights.push(weight.to_string());
    }
    Some(weights)
}

pub(super) fn verify_certificate_in_report(
    certificate: &HalfStepLagCertificate,
    report: &RealRootednessTheoremReport,
) -> Result<(), RecurrenceError> {
    let (rule, _, seeds) = report.recurrence.to_recurrence_parts()?;
    let (rule, _) = rule_with_x_free_left_hand_factor(&rule)?;
    let (lag, _) = recognize(&rule).map_err(RecurrenceError::VerificationFailed)?;
    if certificate.lag != lag
        || seed_weights(&seeds, lag).as_ref() != Some(&certificate.seed_weights)
    {
        return Err(RecurrenceError::VerificationFailed(
            "half-step lag seed certificate differs".into(),
        ));
    }
    Ok(())
}

fn half_step_lag_report(
    rule: &Recurrence,
    first: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    let mut report = RealRootednessTheoremReport::new(
        RealRootednessTheorem::HalfStepLag,
        rule,
        first,
        seeds,
        options,
    );
    let (lag, weight) = match recognize(rule) {
        Ok(x) => x,
        Err(reason) => return Ok(report.not_applicable(reason)),
    };
    if seeds.len() < lag {
        return Ok(report.not_applicable("fewer seeds than the positive lag"));
    }
    let start = first + seeds.len();
    report.normal_form = format!("A_m=P_({first}+m); A_m=1 for m<{lag}; A_(q+{lag})=(1+xD)A_(q+{lag}-1)+c_q x A_q; c_q={} evaluated at n={first}+q+{lag}", weight.text());
    report.lean_route = Some("halfStepLag".into());
    report.lean_declarations = [
        "RealRooted.HalfStepRows.strictInterl_add_of_le",
        "RealRooted.HalfStepRows.isPFPolynomial",
        "RealRooted.HalfStepRows.natDegree_eq_div",
    ]
    .map(str::to_string)
    .to_vec();
    report
        .hypotheses
        .extend(denominator_hypothesis(&weight.denominator, start));
    report.hypotheses.push(sign_hypothesis(
        "positive_lag_weight",
        "c(n)>0 throughout the integer rule range",
        "c(n)",
        &weight,
        start,
        IndexSignCondition::Positive,
    ));
    let weights = seed_weights(seeds, lag);
    report.hypotheses.push(exact_hypothesis(
        "seeded_window",
        "first k rows equal 1; every later seed is a positive lag step",
        weights
            .as_ref()
            .map(|_| "exact coefficient identities".into())
            .ok_or((Some(first), "seeded window hypotheses fail".into())),
    ));
    if let Some(region) = &options.root_region {
        let contains =
            region.lower.is_none() && region.upper.as_ref().is_none_or(|u| !u.is_negative());
        report.hypotheses.push(if contains {
            exact_hypothesis(
                "requested_root_region",
                "requested region contains (-inf,0)",
                Ok("exact".into()),
            )
        } else {
            unresolved_hypothesis(
                "requested_root_region",
                "requested region contains (-inf,0)",
                "the theorem certifies (-inf,0)".into(),
            )
        });
    }
    report.certificate = weights.map(|seed_weights| {
        RouteCertificate::HalfStepLag(HalfStepLagCertificate { lag, seed_weights })
    });
    Ok(report.finish(RealRootednessTheoremConclusion {
        first_index: first, zero_residue_classes: None, interlacing_from: Some(first+1), interlacing: Some(RowInterlacing::Weak), simple_roots: false, root_region: Some(open_root_interval(&RealInterval { lower: None, lower_closed: false, upper: Some(Q::zero()), upper_closed: false })),
        statement: format!("Every P_n is PF with constant term 1, positive leading coefficient, degree floor((n-{first})/{lag}) and negative real roots; P_i interlaces P_j whenever {first}<=i<=j<=i+{lag}. Seeded half-step window theorem, no universal companion-preserver claim."),
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn associated_bell_and_lag_two_families_prove_and_replay() {
        for id in ["A059022", "A059023", "A008299", "A124324", "A136011"] {
            let entry = crate::oeis::by_id(id).unwrap();
            let (rule, seeds) = entry.recurrence_parts().unwrap();
            let first = entry.recurrence_first_index;
            let report = check_half_step_lag_theorem(
                &rule,
                first,
                &seeds,
                &RealRootednessTheoremOptions::default(),
            )
            .unwrap();
            assert!(report.is_proved(), "{id}: {report:?}");
            report.verify().unwrap();
            let mut tampered = report.clone();
            if let Some(RouteCertificate::HalfStepLag(c)) = &mut tampered.certificate {
                c.lag += 1;
            }
            assert!(tampered.verify().is_err());
        }
    }
    #[test]
    fn rejects_nonpositive_weights_and_bad_seeds() {
        let entry = crate::oeis::by_id("A059022").unwrap();
        let (mut rule, mut seeds) = entry.recurrence_parts().unwrap();
        let first = entry.recurrence_first_index;
        rule.terms.last_mut().unwrap().coeff =
            BivarPoly::from_coefficients_normalized(vec![vec![Q::zero(), -Q::one()]]);
        assert!(!check_half_step_lag_theorem(
            &rule,
            first,
            &seeds,
            &RealRootednessTheoremOptions::default()
        )
        .unwrap()
        .is_proved());
        seeds[0] = vec![Q::from_integer(2.into())];
        assert!(!check_half_step_lag_theorem(
            &rule,
            first,
            &seeds,
            &RealRootednessTheoremOptions::default()
        )
        .unwrap()
        .is_proved());
    }
}

#[cfg(test)]
mod boundary_tests {
    use super::*;
    fn q(n: i64) -> Q {
        Q::from_integer(n.into())
    }
    #[test]
    fn lag_one_and_arbitrary_positive_rational_index_weights() {
        let rule = Recurrence {
            terms: vec![
                RecurrenceTerm::new(
                    1,
                    0,
                    BivarPoly::from_coefficients_normalized(vec![
                        vec![q(1), q(1)],
                        vec![q(0), q(1)],
                    ]),
                ),
                RecurrenceTerm::new(
                    1,
                    1,
                    BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(1)]]),
                ),
            ],
            denominator: None,
            inhomogeneous: None,
        };
        let report = check_half_step_lag_theorem(
            &rule,
            7,
            &[vec![q(1)]],
            &RealRootednessTheoremOptions::default(),
        )
        .unwrap();
        assert!(report.is_proved());
        report.verify().unwrap();
        // (n+1) P_n = (n+1)(1+xD)P_(n-1)+x P_(n-2).
        let rule = Recurrence {
            terms: vec![
                RecurrenceTerm::new(
                    1,
                    0,
                    BivarPoly::from_coefficients_normalized(vec![vec![q(1)], vec![q(1)]]),
                ),
                RecurrenceTerm::new(
                    1,
                    1,
                    BivarPoly::from_coefficients_normalized(vec![
                        vec![q(0), q(1)],
                        vec![q(0), q(1)],
                    ]),
                ),
                RecurrenceTerm::new(
                    2,
                    0,
                    BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(1)]]),
                ),
            ],
            denominator: Some(BivarPoly::from_coefficients_normalized(vec![
                vec![q(1)],
                vec![q(1)],
            ])),
            inhomogeneous: None,
        };
        let seeds = rule
            .generate_rows_rational(&[vec![q(1)], vec![q(1)]], 0, 5)
            .unwrap();
        let report =
            check_half_step_lag_theorem(&rule, 0, &seeds, &RealRootednessTheoremOptions::default())
                .unwrap();
        assert!(report.is_proved(), "{report:?}");
        assert_eq!(
            report
                .half_step_lag_certificate()
                .unwrap()
                .seed_weights
                .len(),
            3
        );
        report.verify().unwrap();
        let mut bad = seeds;
        bad[3][1] += q(1);
        assert!(!check_half_step_lag_theorem(
            &rule,
            0,
            &bad,
            &RealRootednessTheoremOptions::default()
        )
        .unwrap()
        .is_proved());
    }
}
