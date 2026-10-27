//! Tests of the root-domain contract, theorem adapters and reports.

use super::*;
use crate::real_rootedness::isolate_real_roots_rational_coeffs;

fn r(n: i64, d: i64) -> Q {
    Q::new(BigInt::from(n), BigInt::from(d))
}

fn strings(values: &[&str]) -> Vec<String> {
    values.iter().map(|s| s.to_string()).collect()
}

fn f_schedule(first_index: usize, first_degree: usize, step: i64) -> RecurrenceTransformation {
    RecurrenceTransformation::FPolynomial {
        schedule: DegreeSchedule::new(first_index, first_degree, step),
    }
}

fn eulerian_request(root_domain: RootDomainRequest) -> RecurrenceTransformRequest {
    RecurrenceTransformRequest {
        transformation: f_schedule(1, 0, 1),
        source_recurrence: None,
        oeis_id: Some("A008292".into()),
        rows: None,
        row_origin: None,
        row_count: Some(8),
        search: None,
        limits: TransformLimits::default(),
        root_domain: Some(root_domain),
    }
}

fn rows_request(
    transformation: RecurrenceTransformation,
    rows: &[&[&str]],
    origin: usize,
    root_domain: RootDomainRequest,
) -> RecurrenceTransformRequest {
    RecurrenceTransformRequest {
        transformation,
        source_recurrence: None,
        oeis_id: None,
        rows: Some(rows.iter().map(|row| strings(row)).collect()),
        row_origin: Some(origin),
        row_count: None,
        search: None,
        limits: TransformLimits::default(),
        root_domain: Some(root_domain),
    }
}

fn closed(lower: &str, upper: &str) -> RootInterval {
    RootInterval {
        lower: Some(RootEndpoint {
            value: lower.into(),
            per_index: None,
            closed: true,
        }),
        upper: Some(RootEndpoint {
            value: upper.into(),
            per_index: None,
            closed: true,
        }),
    }
}

fn claims_with(report: &RootDomainReport, evidence: RootEvidence) -> Vec<&RootDomainClaim> {
    report
        .claims
        .iter()
        .filter(|c| c.evidence == evidence)
        .collect()
}

/// Acceptance: canonical Eulerian A008292 under the f-polynomial map becomes
/// `Q_(n+1) = (x(1+x) Q_n)'`, `Q_1 = 1`; Rolle closes the induction on
/// `[-1, 0]` and the Möbius equivalence transfers it to every Eulerian row.
#[test]
fn eulerian_rolle_induction_on_default_window() {
    let report = transform_recurrence(&eulerian_request(RootDomainRequest::default())).unwrap();
    let root = report.root_domain.as_ref().unwrap();
    let theorem = root.recurrence_theorem.as_ref().unwrap();
    assert!(theorem.matched, "{theorem:?}");
    assert_eq!(
        theorem.factor.as_deref(),
        Some(&strings(&["0", "1", "1"])[..])
    );
    assert_eq!(theorem.derivative_order, Some(1));
    assert_eq!(theorem.index_scalar.as_deref(), Some(&strings(&["1"])[..]));
    assert_eq!(
        theorem.degree_formula.as_deref(),
        Some("deg Q_n = 0 + 1 (n - 1) for n >= 1")
    );
    assert_eq!(theorem.attempts.len(), 1);
    assert_eq!(theorem.attempts[0].origin, "default_window");
    assert!(theorem.attempts[0].success);
    assert_eq!(theorem.attempts[0].interval.to_string(), "[-1, 0]");

    let proposed = claims_with(root, RootEvidence::ProposedInterval);
    assert_eq!(proposed.len(), 1);
    assert_eq!(proposed[0].property, RootProperty::RealRootsIn);
    assert!(!proposed[0].conditional_on.is_empty());

    let induction = claims_with(root, RootEvidence::CompletedRecurrenceInduction);
    assert_eq!(induction.len(), 1);
    let claim = induction[0];
    assert_eq!(claim.subject, RootSubject::Transformed);
    assert_eq!(claim.property, RootProperty::RealRootedIn);
    assert_eq!(claim.interval.as_ref().unwrap().to_string(), "[-1, 0]");
    assert_eq!((claim.first_index, claim.last_index), (1, None));
    assert!(claim.conditional_on.is_empty());

    let family = claims_with(root, RootEvidence::CertifiedFamilyInvariant);
    assert_eq!(family.len(), 1);
    assert_eq!(family[0].subject, RootSubject::Source);
    assert_eq!(family[0].property, RootProperty::RealRootedIn);
    // deg Q_n = d_n, so no root is added at x = 1.
    assert_eq!(
        family[0].interval.as_ref().unwrap().to_string(),
        "(-inf, 0]"
    );
    assert!(report.conditions[0].status == CERTIFIED_ALL_INDICES);

    // The map implication is an equivalence with certified hypotheses.
    let implication = &root.implications[0];
    assert_eq!(implication.theorem, "real_mobius_equivalence");
    assert_eq!(implication.direction, RootImplicationDirection::Equivalent);
    assert_eq!(
        implication.hypotheses[0].status,
        HypothesisStatus::CertifiedAllIndices
    );
    assert_eq!(
        root.implications[1].theorem,
        "nonnegative_coefficient_window"
    );

    // Finite observations are labelled as such and agree with the claims.
    let observed = claims_with(root, RootEvidence::FiniteExactObservation);
    assert!(observed.iter().all(|c| c.last_index == Some(8)));
    assert!(root
        .observations
        .iter()
        .all(|o| o.real_rooted && o.index >= 1));

    // Replay and JSON round trip keep the root-domain section.
    verify_recurrence_transform(&report).unwrap();
    let json = serde_json::to_string(&report).unwrap();
    let back: RecurrenceTransformReport = serde_json::from_str(&json).unwrap();
    assert_eq!(back, report);
    verify_recurrence_transform(&back).unwrap();
    let mut tampered = report.clone();
    tampered.root_domain.as_mut().unwrap().claims.pop();
    assert!(verify_recurrence_transform(&tampered).is_err());
}

#[test]
fn eulerian_root_hull_without_default_window() {
    let request = RootDomainRequest {
        default_window: false,
        ..RootDomainRequest::default()
    };
    let report = transform_recurrence(&eulerian_request(request)).unwrap();
    let root = report.root_domain.unwrap();
    let theorem = root.recurrence_theorem.clone().unwrap();
    let attempt = theorem.attempts.last().unwrap();
    assert_eq!(attempt.origin, "root_hull");
    assert!(attempt.success);
    // The hull of g = x(1+x) and Q_1 = 1 contains [-1, 0].
    let hull = attempt.interval.at(0).unwrap();
    assert!(RealInterval::closed(r(-1, 1), r(0, 1)).is_subset_of(&hull));
    assert!(claims_with(&root, RootEvidence::ProposedInterval).is_empty());
    assert_eq!(
        claims_with(&root, RootEvidence::CertifiedFamilyInvariant).len(),
        1
    );
}

#[test]
fn bounded_failed_interval_search() {
    let request = RootDomainRequest {
        candidate_intervals: vec![closed("0", "1"), closed("-1/2", "0"), closed("-2", "2")],
        default_window: false,
        max_candidate_intervals: 2,
        ..RootDomainRequest::default()
    };
    let report = transform_recurrence(&eulerian_request(request)).unwrap();
    let root = report.root_domain.unwrap();
    let theorem = root.recurrence_theorem.clone().unwrap();
    assert_eq!(theorem.attempts.len(), 2);
    assert!(theorem.attempts.iter().all(|a| !a.success));
    assert!(theorem.attempts[0].detail.contains("g = "));
    assert!(root.notes.iter().any(|n| n.contains("were not tried")));
    assert!(root
        .notes
        .iter()
        .any(|n| n.contains("no candidate interval")));
    assert!(claims_with(&root, RootEvidence::CompletedRecurrenceInduction).is_empty());
    assert!(claims_with(&root, RootEvidence::CertifiedFamilyInvariant).is_empty());
}

/// `1 + z + z^2`: the fixed window holds, but the transformation does not
/// prove splitting.
#[test]
fn window_does_not_prove_splitting() {
    let request = rows_request(
        f_schedule(0, 2, 0),
        &[&["1", "1", "1"]],
        0,
        RootDomainRequest::default(),
    );
    let report = transform_recurrence(&request).unwrap();
    assert_eq!(report.transformed_states[0][0], strings(&["1", "3", "3"]));
    let root = report.root_domain.unwrap();
    assert!(root.recurrence_theorem.is_none());
    let proposed = claims_with(&root, RootEvidence::ProposedInterval);
    assert_eq!(proposed.len(), 1);
    assert_eq!(
        proposed[0].interval.as_ref().unwrap().to_string(),
        "[-1, 0]"
    );
    assert!(root.observations.iter().all(|o| !o.real_rooted));
    assert!(root.claims.iter().all(
        |c| c.property != RootProperty::RealRootedIn && c.property != RootProperty::RealRooted
    ));
    assert!(root.implications[1]
        .converse
        .as_ref()
        .unwrap()
        .contains("1 + z + z^2"));
}

/// `P_n = x + n`: the roots escape to minus infinity while the transformed
/// roots stay in `(-1, 0]`; index-dependent intervals and rays.
#[test]
fn escaping_roots_index_dependent_interval_and_rays() {
    let source = RecurrenceJson::from_recurrence_rational(
        &Recurrence {
            terms: vec![RecurrenceTerm::new(
                1,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![Q::one()]]),
            )],
            denominator: None,
            inhomogeneous: Some(BivarPoly::from_coefficients_normalized(vec![
                vec![Q::one()],
            ])),
        },
        0,
        &[vec![Q::zero(), Q::one()]],
        None,
    );
    let minus_n = RootEndpoint {
        value: "0".into(),
        per_index: Some("-1".into()),
        closed: true,
    };
    let at_minus_n = RootInterval {
        lower: Some(minus_n.clone()),
        upper: Some(minus_n),
    };
    let ray = RootInterval {
        lower: None,
        upper: Some(RootEndpoint {
            value: "0".into(),
            per_index: None,
            closed: true,
        }),
    };
    let mut request = RecurrenceTransformRequest::from_source(f_schedule(0, 1, 0), source, 12);
    request.root_domain = Some(RootDomainRequest {
        candidate_intervals: vec![at_minus_n.clone()],
        assumptions: vec![
            RootAssumption {
                subject: RootSubject::Source,
                component: 0,
                property: RootProperty::RealRootedIn,
                interval: Some(at_minus_n),
                first_index: None,
                last_index: None,
                note: Some("P_n = x + n".into()),
            },
            RootAssumption {
                subject: RootSubject::Source,
                component: 0,
                property: RootProperty::RealRootedIn,
                interval: Some(ray),
                first_index: None,
                last_index: None,
                note: None,
            },
        ],
        ..RootDomainRequest::default()
    });
    let report = transform_recurrence(&request).unwrap();
    let root = report.root_domain.as_ref().unwrap();
    // Q_n = n + (n + 1) x.
    assert_eq!(report.transformed_states[11][0], strings(&["11", "12"]));
    let source_obs = claims_with(root, RootEvidence::FiniteExactObservation)
        .into_iter()
        .find(|c| c.subject == RootSubject::Source)
        .unwrap();
    // The hull of the enclosures of the roots 0, -1, ..., -11.
    let hull = source_obs.interval.as_ref().unwrap().at(0).unwrap();
    assert!(RealInterval::closed(r(-11, 1), r(0, 1)).is_subset_of(&hull));
    assert!(hull.is_subset_of(&RealInterval::closed(r(-12, 1), r(0, 1))));
    let target_obs = claims_with(root, RootEvidence::FiniteExactObservation)
        .into_iter()
        .find(|c| c.subject == RootSubject::Transformed)
        .unwrap();
    let hull = target_obs.interval.as_ref().unwrap().at(0).unwrap();
    assert!(hull.is_subset_of(&RealInterval::closed(r(-1, 1), r(0, 1))));
    assert!(root
        .assumption_checks
        .iter()
        .all(|c| c.status == "consistent_with_rows" && c.checked_indices.len() == 12));
    // The Möbius image of [-n, -n] has endpoints rational in n, so its hull
    // over n >= 0, (-inf, 0], is mapped (to [-1, 0]); the ray maps exactly.
    let transferred: Vec<_> = claims_with(root, RootEvidence::SuppliedAssumption)
        .into_iter()
        .filter(|c| c.subject == RootSubject::Transformed)
        .collect();
    assert_eq!(transferred.len(), 2);
    assert_eq!(transferred[0].property, RootProperty::RealRootedIn);
    assert_eq!(
        transferred[0].interval.as_ref().unwrap().to_string(),
        "[-1, 0]"
    );
    assert!(transferred[0]
        .derivation
        .contains("hull over the index domain"));
    assert!(!transferred[1].derivation.contains("hull"));
    assert_eq!(
        transferred[1].interval.as_ref().unwrap().to_string(),
        "[-1, 0]"
    );
    assert!(transferred
        .iter()
        .all(|c| c.conditional_on[0].starts_with("assumption")));
    // The additive term is outside the Rolle form.
    let theorem = root.recurrence_theorem.as_ref().unwrap();
    assert!(!theorem.matched);
    assert_eq!(
        theorem.reason.as_deref(),
        Some("the rule has an additive term")
    );
    verify_recurrence_transform(&report).unwrap();
}

/// Gamma: `P` real-rooted iff every root of `G` is real and at most `1/4`.
#[test]
fn gamma_quarter_bound_and_insufficient_real_rootedness() {
    let quarter = RootInterval {
        lower: None,
        upper: Some(RootEndpoint {
            value: "1/4".into(),
            per_index: None,
            closed: true,
        }),
    };
    let assume = |property, interval| RootAssumption {
        subject: RootSubject::Transformed,
        component: 0,
        property,
        interval,
        first_index: None,
        last_index: None,
        note: None,
    };
    let request = rows_request(
        RecurrenceTransformation::Gamma {
            schedule: DegreeSchedule::new(0, 2, 0),
        },
        &[&["1", "1", "1"]],
        0,
        RootDomainRequest {
            assumptions: vec![
                assume(RootProperty::RealRootedIn, Some(quarter)),
                assume(RootProperty::RealRooted, None),
            ],
            ..RootDomainRequest::default()
        },
    );
    let report = transform_recurrence(&request).unwrap();
    // G = 1 - z.
    assert_eq!(report.transformed_states[0][0], strings(&["1", "-1"]));
    let root = report.root_domain.unwrap();
    let implication = &root.implications[0];
    assert_eq!(implication.theorem, "gamma_quarter_bound");
    assert!(implication.converse.as_ref().unwrap().contains("1 - z"));
    assert_eq!(root.implications.len(), 1, "no window for gamma");
    let source = root
        .observations
        .iter()
        .find(|o| o.subject == RootSubject::Source)
        .unwrap();
    let gamma = root
        .observations
        .iter()
        .find(|o| o.subject == RootSubject::Transformed)
        .unwrap();
    assert!(!source.real_rooted && gamma.real_rooted);
    assert_eq!(root.assumption_checks[0].status, "refuted_on_rows");
    assert_eq!(root.assumption_checks[1].status, "consistent_with_rows");
    // Real-rootedness of G alone transfers nothing to P.
    assert!(root.claims.iter().all(
        |c| c.subject != RootSubject::Source || c.evidence != RootEvidence::SuppliedAssumption
    ));
}

#[test]
fn symmetric_pair_needs_interlacing_and_is_one_directional() {
    let decomposition = |d| RecurrenceTransformation::SymmetricDecomposition {
        schedule: DegreeSchedule::new(0, d, 0),
    };
    // a = (1 + x)^3, b = x: both real-rooted, P = 1 + 3x + 4x^2 + x^3 is not.
    let report = transform_recurrence(&rows_request(
        decomposition(3),
        &[&["1", "3", "4", "1"]],
        0,
        RootDomainRequest::default(),
    ))
    .unwrap();
    assert_eq!(
        report.transformed_states[0],
        vec![strings(&["1", "3", "3", "1"]), strings(&["0", "1"])]
    );
    let root = report.root_domain.unwrap();
    let implication = &root.implications[0];
    assert_eq!(
        implication.direction,
        RootImplicationDirection::TargetToSource
    );
    let interlacing = implication.hypotheses.last().unwrap();
    assert_eq!(interlacing.status, HypothesisStatus::FailedOnRows);
    let observed = |subject, component| {
        root.observations
            .iter()
            .find(|o| o.subject == subject && o.component == component)
            .unwrap()
            .real_rooted
    };
    assert!(!observed(RootSubject::Source, 0));
    assert!(observed(RootSubject::Transformed, 0) && observed(RootSubject::Transformed, 1));

    // P = (1 + 2x)^2 at d = 2 is real-rooted, a = 1 + x + x^2 is not.
    let report = transform_recurrence(&rows_request(
        decomposition(2),
        &[&["1", "4", "4"]],
        0,
        RootDomainRequest::default(),
    ))
    .unwrap();
    assert_eq!(report.transformed_states[0][0], strings(&["1", "1", "1"]));
    let root = report.root_domain.unwrap();
    assert_eq!(
        root.implications[0].hypotheses.last().unwrap().status,
        HypothesisStatus::FailedOnRows
    );
    assert!(root.implications[0]
        .converse
        .as_ref()
        .unwrap()
        .contains("(1 + 2x)^2"));

    // Eulerian-like h-vector 1 + 4x + x^2 at d = 2: a = 1 + 4x + x^2, b = 0.
    let report = transform_recurrence(&rows_request(
        decomposition(2),
        &[&["1", "4", "1"]],
        0,
        RootDomainRequest::default(),
    ))
    .unwrap();
    assert_eq!(
        report.root_domain.unwrap().implications[0]
            .hypotheses
            .last()
            .unwrap()
            .status,
        HypothesisStatus::CheckedSuppliedRows
    );
}

/// `P = (xQ)'`: real-rooted `Q` gives real-rooted `P`, not conversely.
#[test]
fn derivative_lift_is_one_directional() {
    let q = vec![Q::one(), Q::one(), r(1, 3)];
    let p = vec_poly::derivative(&vec_poly::shift(&q));
    assert_eq!(p, vec![Q::one(), r(2, 1), Q::one()]);
    assert!(isolate_real_roots_rational_coeffs(&p).is_real_rooted());
    assert!(!isolate_real_roots_rational_coeffs(&q).is_real_rooted());
    // As a diagonal map lambda_k = k + 1 (P = (xQ)'), and its inverse.
    let up: Vec<Q> = (1..=9).map(|k| r(k, 1)).collect();
    let down: Vec<Q> = (1..=9).map(|k| r(1, k)).collect();
    let forward = jensen_fixed_cap_check(&up, 8).unwrap();
    assert_eq!(forward.status, JensenCapStatus::PreservesUpToCap);
    assert_eq!(forward.evidence, FIXED_CAP_EVIDENCE);
    assert!(forward.statement.contains("not a proof for degree 9"));
    let backward = jensen_fixed_cap_check(&down, 2).unwrap();
    assert_eq!(backward.status, JensenCapStatus::FailsAtCap);
    assert_eq!(backward.jensen_polynomial, strings(&["1", "1", "1/3"]));
    // Mixed signs and degree deficiency draw no conclusion.
    let mixed = [Q::one(), Q::zero(), r(-1, 1)];
    assert_eq!(
        jensen_fixed_cap_check(&mixed, 2).unwrap().status,
        JensenCapStatus::Inconclusive
    );
    let deficient = [Q::one(), Q::one(), Q::zero()];
    assert_eq!(
        jensen_fixed_cap_check(&deficient, 2).unwrap().status,
        JensenCapStatus::Inconclusive
    );
    assert!(jensen_fixed_cap_check(&up, 9).is_err());
}

#[test]
fn rolle_form_recognition() {
    let bivar = |rows: &[&[i64]]| {
        BivarPoly::from_coefficients_normalized(
            rows.iter()
                .map(|row| row.iter().map(|&c| r(c, 1)).collect())
                .collect(),
        )
    };
    // Q_m = (m + 1) D^2 (x^2 (1 + x) Q_(m-1)) = (m+1) [g'' Q + 2 g' Q' + g Q''].
    let g = [0, 0, 1, 1];
    let scaled = |p: &[i64], k: i64| -> Vec<Vec<i64>> {
        vec![
            p.iter().map(|c| c * k).collect(),
            p.iter().map(|c| c * k).collect(),
        ]
    };
    let g1 = [0, 2, 3];
    let g2 = [2, 6];
    let as_rows = |v: Vec<Vec<i64>>| bivar(&v.iter().map(Vec::as_slice).collect::<Vec<_>>());
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, as_rows(scaled(&g2, 1))),
            RecurrenceTerm::new(1, 1, as_rows(scaled(&g1, 2))),
            RecurrenceTerm::new(1, 2, as_rows(scaled(&g, 1))),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let form = match_rolle_form(&rule).unwrap();
    assert_eq!(form.order, 2);
    assert_eq!(form.factor, vec![r(0, 1), r(0, 1), r(1, 1), r(1, 1)]);
    assert_eq!(form.scalar, vec![r(1, 1), r(1, 1)]);
    assert!(scalar_nonvanishing(&form.scalar, 0));
    assert!(!scalar_nonvanishing(&[r(-1, 1), r(1, 1)], 2));
    // A wrong lower coefficient is rejected.
    let mut broken = rule.clone();
    broken.terms[0] = RecurrenceTerm::new(1, 0, as_rows(scaled(&g2, 2)));
    assert!(match_rolle_form(&broken).is_err());
    // Two lags are rejected.
    let mut lagged = rule;
    lagged.terms[0].offset = 2;
    assert!(match_rolle_form(&lagged).is_err());
}

#[test]
fn mobius_interval_images() {
    let zero = Q::zero();
    let one = Q::one();
    let ray = RealInterval::at_most(zero.clone());
    // Reversal y -> 1/y maps (-inf, 0] to (-inf, 0) (the pole 0 is dropped).
    let image = mobius_interval_image(&ray, &zero, &one, &one, &zero).unwrap();
    assert_eq!(image.to_string(), "(-inf, 0)");
    // f-polynomial inverse y -> y/(1 - y) maps (-inf, 0] to (-1, 0].
    let image = mobius_interval_image(&ray, &one, &zero, &-one.clone(), &one).unwrap();
    assert_eq!(image.to_string(), "(-1, 0]");
    // A pole in the interior has no interval image.
    let around = RealInterval::closed(r(-1, 1), r(1, 1));
    assert!(mobius_interval_image(&around, &zero, &one, &one, &zero).is_none());
    // A decreasing affine map swaps the ends.
    let image = mobius_interval_image(&around, &r(-2, 1), &one, &zero, &one).unwrap();
    assert_eq!(image.to_string(), "[-1, 3]");
}

#[test]
fn request_validation_and_backward_compatible_json() {
    let bad = RootDomainRequest {
        assumptions: vec![RootAssumption {
            subject: RootSubject::Source,
            component: 0,
            property: RootProperty::RealRooted,
            interval: Some(closed("0", "1")),
            first_index: None,
            last_index: None,
            note: None,
        }],
        ..RootDomainRequest::default()
    };
    assert!(bad.validate().is_err());
    let too_many = RootDomainRequest {
        max_observed_rows: ROOT_DOMAIN_MAX_ITEMS + 1,
        ..RootDomainRequest::default()
    };
    assert!(too_many.validate().is_err());
    let mut request = eulerian_request(RootDomainRequest::default());
    request.root_domain = None;
    let value = serde_json::to_value(&request).unwrap();
    assert!(value.get("root_domain").is_none());
    let report = transform_recurrence(&request).unwrap();
    assert!(serde_json::to_value(&report)
        .unwrap()
        .get("root_domain")
        .is_none());
    // `{}` selects the defaults.
    let parsed: RootDomainRequest = serde_json::from_str("{}").unwrap();
    assert_eq!(parsed, RootDomainRequest::default());
    let endpoint: RootEndpoint = serde_json::from_str(r#"{"value": "-1"}"#).unwrap();
    assert!(endpoint.closed);
}

#[test]
fn root_enclosure_requests() {
    let request = RootEnclosureRequest {
        // (x + 1)^3 x^2, x^2 - 2, 0
        polynomials: vec![
            strings(&["0", "0", "1", "3", "3", "1"]),
            strings(&["-2", "0", "1"]),
            Vec::new(),
        ],
        max_width: Some("1/1000".into()),
        interval: Some(closed("-1", "0")),
        first_index: 0,
    };
    let report = root_enclosures(&request).unwrap();
    assert_eq!(report.schema, ROOT_ENCLOSURES_SCHEMA);
    let first = &report.items[0];
    assert_eq!(first.valuation, 2);
    assert_eq!(first.roots.len(), 2);
    assert_eq!(first.roots[0].multiplicity, 3);
    assert_eq!(first.roots[1].lower, "0");
    assert_eq!(first.real_rooted_in_interval, Some(true));
    let second = &report.items[1];
    assert!(second.real_rooted);
    assert_eq!(second.real_roots_in_interval, Some(false));
    let zero = &report.items[2];
    assert_eq!(zero.degree, None);
    assert_eq!(zero.real_rooted_in_interval, Some(true));
    let bad = RootEnclosureRequest {
        max_width: Some("0".into()),
        ..request
    };
    assert!(root_enclosures(&bad).is_err());
}
