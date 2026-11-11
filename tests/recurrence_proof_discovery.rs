//! General theorem rules: parameter families, late failures, root domains and replay.

use polytool::real_rootedness::RealInterval;
use polytool::recurrence::*;

mod common;
use common::q;

fn b(rows: &[&[i64]]) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(
        rows.iter()
            .map(|r| r.iter().map(|&n| q(n)).collect())
            .collect(),
    )
}
fn product(coefficient: &[&[i64]], denominator: Option<&[&[i64]]>, lag: usize) -> Recurrence {
    Recurrence {
        terms: vec![RecurrenceTerm::new(lag, 0, b(coefficient))],
        denominator: denominator.map(b),
        inhomogeneous: None,
    }
}
fn options() -> RealRootednessTheoremOptions {
    RealRootednessTheoremOptions {
        maximum_start_delay: 0,
        ..Default::default()
    }
}
fn record(rule: &Recurrence, seeds: &[&[i64]], first: usize) -> RecurrenceJson {
    let seeds: Vec<_> = seeds
        .iter()
        .map(|r| r.iter().map(|&n| q(n)).collect())
        .collect();
    RecurrenceJson::from_recurrence_rational(rule, first, &seeds, None)
}
fn checked(record: &RecurrenceJson) -> RecurrenceProofDiscoveryReport {
    let report = discover_recurrence_proofs(record, &options()).unwrap();
    report.verify().unwrap();
    // The first eight theorems are always attempted, in order. The rational
    // PF Riordan route can also apply to the index-free rules in this file.
    let theorems: Vec<RealRootednessTheorem> = report.attempts.iter().map(|a| a.theorem).collect();
    assert_eq!(&theorems[..8], &RealRootednessTheorem::ALL[..8]);
    assert!(matches!(
        &theorems[8..],
        [] | [RealRootednessTheorem::RationalPfRiordan]
            | [RealRootednessTheorem::DegreeLawReversal]
            | [
                RealRootednessTheorem::RationalPfRiordan,
                RealRootednessTheorem::DegreeLawReversal
            ]
    ));
    assert_eq!(report.lean_replay, "not_attempted");
    report
}
fn product_check(rule: &Recurrence, seeds: &[&[i64]]) -> RealRootednessTheoremReport {
    let (_, index, seeds) = record(rule, seeds, 0).to_recurrence_parts().unwrap();
    let result = check_product_recurrence_theorem(rule, index, &seeds, &options()).unwrap();
    result.verify().unwrap();
    result
}

#[test]
fn parameterized_linear_products_allow_both_root_signs_and_scalars() {
    for slope in [-3, -1, 1, 4] {
        for constant in [-5, 0, 7] {
            // (n^2+1) P_n = (slope*x+n+constant) P_(n-1).
            let rule = product(&[&[constant, slope], &[1]], Some(&[&[1], &[0], &[1]]), 1);
            let report = checked(&record(&rule, &[&[1]], 0));
            assert_eq!(report.status, "proved");
            assert!(report.attempts[0].is_proved());
        }
    }
}

#[test]
fn quadratic_degree_drops_repeated_roots_and_irrational_real_roots() {
    let failed = product_check(&product(&[&[-1, -2, -2], &[0, 1, 1]], None, 1), &[&[1]]);
    assert!(!failed.is_proved());
    // (n-2)^2*x^2 - 1: a nonzero constant at n=2, real roots elsewhere.
    assert!(product_check(
        &product(&[&[-1, 0, 4], &[0, 0, -4], &[0, 0, 1]], None, 1),
        &[&[1]]
    )
    .is_proved());
    // x^2 - n is irreducible over Q(n), with real roots for all n>=1.
    assert!(product_check(&product(&[&[0, 0, 1], &[-1]], None, 1), &[&[1]]).is_proved());
    let repeated = product_check(&product(&[&[0, 0, 1], &[0, 2], &[1]], None, 1), &[&[1]]);
    assert!(repeated.is_proved());
    assert!(!repeated.conclusion.unwrap().simple_roots);
}

#[test]
fn late_failures_are_not_hidden_by_finite_prefixes() {
    // x^2+n-10000 is real-rooted for ten thousand steps only.
    let report = product_check(&product(&[&[-10000, 0, 1], &[1]], None, 1), &[&[1]]);
    assert!(report
        .hypotheses
        .iter()
        .any(|h| h.name.ends_with("discriminant")
            && h.status
                == LiftHypothesisStatus::Failed {
                    witness: Some(10001.into())
                }));
    let singular = product_check(
        &product(&[&[1, 1], &[1]], Some(&[&[-10000], &[1]]), 1),
        &[&[1]],
    );
    assert!(singular
        .hypotheses
        .iter()
        .any(|h| h.name == "left_factor_nonzero"
            && h.status
                == LiftHypothesisStatus::Failed {
                    witness: Some(10000.into())
                }));
    assert!(!product_check(&product(&[&[-10000, -10000], &[1, 1]], None, 1), &[&[1]]).is_proved());
}

#[test]
fn higher_degree_factorization_and_separable_multipliers() {
    // (x+n)(x+2n)(x+3n), expanded: no sequence-name information.
    let cubic = product(&[&[0, 0, 0, 1], &[0, 0, 6], &[0, 11], &[6]], None, 1);
    assert!(product_check(&cubic, &[&[1]]).is_proved());
    assert!(product_check(
        &product(&[&[0, -4, 0, 1], &[0, -4, 0, 1]], None, 1),
        &[&[1]]
    )
    .is_proved());
    // (x+n)(x^2+1) has a nonreal factor.
    assert!(!product_check(&product(&[&[0, 1, 0, 1], &[1, 0, 1]], None, 1), &[&[1]]).is_proved());
    let unsupported = product_check(&product(&[&[1, 0, 0, 1], &[0, 0, 0, 1]], None, 1), &[&[1]]);
    assert_eq!(unsupported.status, RealRootednessTheoremStatus::Unresolved);
}

#[test]
fn delayed_products_check_every_seed_and_cancel_left_factors() {
    let rule = product(&[&[1, 1], &[1]], None, 3);
    assert!(product_check(&rule, &[&[1], &[1, -1], &[2, 3, 1]]).is_proved());
    assert!(!product_check(&rule, &[&[1], &[1, 0, 1], &[1]]).is_proved());
    assert!(!product_check(&rule, &[&[1], &[0], &[1]]).is_proved());
    let canceled = product_check(
        &product(&[&[0, 1, 1], &[0, 1]], Some(&[&[0, 1]]), 1),
        &[&[1]],
    );
    assert!(canceled.is_proved());
    assert!(canceled.left_hand_factor_normalization.is_some());
}

#[test]
fn root_regions_and_transfer_seed_reuse() {
    let rule = product(&[&[1, 1], &[1]], None, 1);
    let seeds = vec![vec![q(1)]];
    let negative = RealRootednessTheoremOptions {
        root_region: Some(RealInterval {
            lower: None,
            lower_closed: false,
            upper: Some(q(0)),
            upper_closed: false,
        }),
        ..options()
    };
    let proof = check_product_recurrence_theorem(&rule, 0, &seeds, &negative).unwrap();
    assert!(proof.is_proved());
    proof.verify().unwrap();
    let parent = prove_parent_real_rootedness(
        &ParentRecurrenceRecord {
            rule: &rule,
            seeds: &seeds,
            first_index: 0,
        },
        0,
        &[],
        false,
    )
    .unwrap();
    assert_eq!(parent.route, Some(VeroneseParentProofRoute::LiuWang));
    assert!(parent.interlacing.is_some());
    assert_eq!(parent.status, VeroneseParentProofStatus::Proved);
    // New higher-degree product seeds do not acquire an interlacing claim.
    let cubic = product(&[&[0, 0, 0, 1], &[0, 0, 6], &[0, 11], &[6]], None, 1);
    let parent = prove_parent_real_rootedness(
        &ParentRecurrenceRecord {
            rule: &cubic,
            seeds: &seeds,
            first_index: 0,
        },
        0,
        &[],
        false,
    )
    .unwrap();
    assert_eq!(
        parent.route,
        Some(VeroneseParentProofRoute::ProductRecurrence)
    );
    assert!(parent.interlacing.is_none());
    let positive = product(&[&[-1, 1], &[-1]], None, 1);
    assert!(product_check(&positive, &[&[1]]).is_proved());
    assert!(
        !check_product_recurrence_theorem(&positive, 0, &seeds, &negative)
            .unwrap()
            .is_proved()
    );
    let quad = product(&[&[-1, 0, 4], &[0, 8], &[4]], None, 1);
    assert!(
        check_product_recurrence_theorem(&quad, 0, &seeds, &negative)
            .unwrap()
            .is_proved()
    );
}

#[test]
fn moving_quadratic_roots_obey_both_open_endpoints() {
    // (n+1)^2*x^2 - 1 has roots +/-1/(n+1); also test a negative leader.
    for scalar in [-1, 1] {
        let rule = product(
            &[&[-scalar, 0, scalar], &[0, 0, 2 * scalar], &[0, 0, scalar]],
            None,
            1,
        );
        let opts = RealRootednessTheoremOptions {
            root_region: Some(RealInterval {
                lower: Some(q(-1)),
                upper: Some(q(1)),
                lower_closed: false,
                upper_closed: false,
            }),
            ..options()
        };
        let report = check_product_recurrence_theorem(&rule, 0, &[vec![q(1)]], &opts).unwrap();
        assert!(report.is_proved());
        report.verify().unwrap();
        let mut boundary = opts;
        boundary.root_region.as_mut().unwrap().lower = Some(-q(1) / q(2));
        let failed = check_product_recurrence_theorem(&rule, 0, &[vec![q(1)]], &boundary).unwrap();
        assert!(!failed.is_proved());
        failed.verify().unwrap();
    }
}

#[test]
fn unsupported_product_attestations_still_require_exact_record_and_index() {
    // (1+n*x)^3 is beyond the nonconstant-leading cubic factorizer.
    let rule = product(&[&[1], &[0, 3], &[0, 0, 3], &[0, 0, 0, 1]], None, 1);
    let seeds = vec![vec![q(1)]];
    let parent = ParentRecurrenceRecord {
        rule: &rule,
        seeds: &seeds,
        first_index: 0,
    };
    let input = record(&rule, &[&[1]], 0);
    let attestation = |recurrence, first_index| VeroneseParentAttestation {
        label: "external hypothesis".into(),
        recurrence,
        first_index,
        statement: "all rows real-rooted with nonnegative coefficients".into(),
        source: "test".into(),
    };
    let check = |evidence: Vec<VeroneseParentAttestation>| {
        prove_parent_real_rootedness(&parent, 0, &evidence, false).unwrap()
    };
    assert_eq!(check(vec![]).status, VeroneseParentProofStatus::NotProved);
    assert_eq!(
        check(vec![attestation(
            record(&product(&[&[1, 1]], None, 1), &[&[1]], 0),
            0
        )])
        .status,
        VeroneseParentProofStatus::NotProved
    );
    assert_eq!(
        check(vec![attestation(input.clone(), 1)]).status,
        VeroneseParentProofStatus::NotProved
    );
    assert_eq!(
        check(vec![attestation(input, 0)]).status,
        VeroneseParentProofStatus::ConditionalOnAttestation
    );
    // Direct general discovery has no attestation route at all.
    assert_eq!(checked(&record(&rule, &[&[1]], 0)).status, "not_proved");
}

#[test]
fn general_differential_and_three_term_routes_are_discoverable() {
    for a in [1, 3, 7] {
        let mut differential = product(&[&[a, 1]], None, 1);
        differential
            .terms
            .push(RecurrenceTerm::new(1, 1, b(&[&[0, 1]])));
        let report = checked(&record(&differential, &[&[1]], 0));
        assert!(report
            .attempts
            .iter()
            .any(|r| r.theorem == RealRootednessTheorem::LiuWang && r.is_proved()));
        let mut three = product(&[&[0, 1]], None, 1);
        three.terms.push(RecurrenceTerm::new(2, 0, b(&[&[-a]])));
        let report = checked(&record(&three, &[&[1], &[0, 1]], 0));
        assert!(report
            .attempts
            .iter()
            .any(|r| r.theorem == RealRootednessTheorem::ThreeTermRecurrence && r.is_proved()));
    }
}

#[test]
fn discovery_replay_rejects_tampering_and_does_not_infer_a_proof_from_failure() {
    let input = record(&product(&[&[1, 1], &[1]], None, 1), &[&[1]], 0);
    let report = checked(&input);
    for change in 0..4 {
        let mut tampered = report.clone();
        match change {
            0 => tampered.attempts.clear(),
            1 => tampered.lean_replay = "verified".into(),
            2 => {
                tampered.attempts[0]
                    .conclusion
                    .as_mut()
                    .unwrap()
                    .simple_roots = true
            }
            _ => {
                tampered.recurrence.initial_polynomials[0] =
                    vec!["1".into(), "0".into(), "1".into()]
            }
        }
        assert!(tampered.verify().is_err());
    }
    let failure = checked(&record(&product(&[&[1, 0, 1]], None, 1), &[&[1]], 0));
    assert_eq!(failure.status, "not_proved");
    assert!(discover_recurrence_proofs(
        &record(&product(&[&[1]], None, 2), &[&[1]], 0),
        &options()
    )
    .is_err());
    assert!(discover_recurrence_proofs(
        &record(&product(&[&[1]], None, 1), &[&[1]], usize::MAX),
        &options()
    )
    .is_err());
}

#[test]
fn trivially_real_rooted_classifies_the_product_seed_class() {
    // (n^2+1) P_n = (x+n) P_(n-1): product class, proved.
    let rule = product(&[&[0, 1], &[1]], Some(&[&[1], &[0], &[1]]), 1);
    let rec = record(&rule, &[&[1]], 0);
    let (_, first, seeds) = rec.to_recurrence_parts().unwrap();
    assert!(is_trivially_real_rooted(&rule, first, &seeds).unwrap());
    let report = checked(&rec);
    assert!(report.trivially_real_rooted);
    assert_eq!(report.trivially_real_rooted_if_rule_holds, Some(true));
    assert!(report.attempts[0].trivially_real_rooted);
    assert!(report.attempts[1..]
        .iter()
        .all(|a| !a.trivially_real_rooted));

    // A fitted record is never classified from the fit alone.
    let mut fitted = rec.clone();
    fitted.search = Some(RecurrenceJsonSearch {
        recurrence_text: rule.to_string(),
        source_rows: 1,
        skip_prefix: 0,
        unknowns: 1,
        weighted_unknowns: 1,
        equations: 1,
        fit_polynomials: 1,
        verification_polynomials: 0,
        verification_equations: 0,
        candidates_tried: 1,
        options: RecurrenceOptionsJson::from(&RecurrenceOptions::default()),
        work_units: None,
    });
    let fitted_report = checked(&fitted);
    assert!(!fitted_report.trivially_real_rooted);
    assert_eq!(
        fitted_report.trivially_real_rooted_if_rule_holds,
        Some(true)
    );

    // A product whose multiplier has a non-real root (x^2 + 1) is not proved.
    let bad = product(&[&[1, 0, 1]], None, 1);
    let (_, first, seeds) = record(&bad, &[&[1]], 0).to_recurrence_parts().unwrap();
    assert!(!is_trivially_real_rooted(&bad, first, &seeds).unwrap());
    // Three-term families are proved, but not trivially.
    let mut three = product(&[&[0, 1]], None, 1);
    three.terms.push(RecurrenceTerm::new(2, 0, b(&[&[-1]])));
    let report = checked(&record(&three, &[&[1], &[0, 1]], 0));
    assert_eq!(report.status, "proved");
    assert!(!report.trivially_real_rooted);
    assert_eq!(report.trivially_real_rooted_if_rule_holds, Some(false));
}

/// Theorem reports serialized before `trivially_real_rooted` existed keep
/// their exact JSON (stored certificates are identified by a hash of it).
#[test]
fn theorem_report_serialization_omits_a_false_trivially_real_rooted() {
    // A product whose multiplier has a non-real root (x^2 + 1): not proved.
    let report = product_check(&product(&[&[1, 0, 1]], None, 1), &[&[1]]);
    assert!(!report.trivially_real_rooted);
    let json = serde_json::to_value(&report).unwrap();
    assert!(json.get("trivially_real_rooted").is_none());
    let reparsed: RealRootednessTheoremReport = serde_json::from_value(json.clone()).unwrap();
    assert_eq!(serde_json::to_value(&reparsed).unwrap(), json);
}
