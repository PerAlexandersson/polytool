//! Closed-endpoint and forward gamma proof regressions.

use crate::common::q;
use polytool::recurrence::*;

fn b(rows: &[&[i64]]) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(
        rows.iter()
            .map(|r| r.iter().map(|&n| q(n)).collect())
            .collect(),
    )
}

fn endpoint(a: i64, slope: i64, gamma: i64, upper: i64) -> Recurrence {
    Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, b(&[&[a - upper * slope, slope], &[-1]])),
            RecurrenceTerm::new(1, 1, b(&[&[-gamma * upper, gamma]])),
        ],
        denominator: None,
        inhomogeneous: None,
    }
}

#[test]
fn growing_endpoint_multiplicity_and_positive_coefficients() {
    for (d, slope, seed) in [
        (3, 1, 1),
        (4, 1, 1),
        (6, 1, 1),
        (5, 1, 1),
        (3, 3, 3),
        (4, 4, 4),
        (6, 6, 6),
        (7, 7, 7),
    ] {
        let rule = endpoint(d, slope, d, 0);
        let seeds = vec![vec![q(seed)]];
        let report = check_ma_wang_interval_theorem(&rule, 0, &seeds, &Default::default()).unwrap();
        assert!(report.is_proved(), "{report:?}");
        assert_eq!(
            report
                .conclusion
                .as_ref()
                .unwrap()
                .root_region
                .as_ref()
                .unwrap()
                .to_string(),
            "(-inf, 0]"
        );
        assert!(report.ma_wang_endpoint_certificate().is_some());
        report.verify().unwrap();
        for (n, row) in rule
            .generate_rows_rational(&seeds, 0, 32)
            .unwrap()
            .iter()
            .enumerate()
        {
            assert_eq!(row.iter().position(|c| c != &q(0)), Some(n / d as usize));
            assert!(row.iter().all(|c| c >= &q(0)));
            let integers: Vec<_> = row.iter().map(|c| c.to_integer()).collect();
            assert!(polytool::is_real_rooted_bigint_coeffs(&integers));
        }
    }
}

#[test]
fn endpoint_is_general_and_same_degree_is_checked() {
    // The lattice phase need not have a=d, and the endpoint need not be 0.
    let rule = endpoint(2, 2, 3, 5);
    let seeds = vec![vec![q(7)]];
    let report = check_ma_wang_interval_theorem(&rule, 0, &seeds, &Default::default()).unwrap();
    assert!(report.is_proved(), "{report:?}");
    assert_eq!(report.ma_wang_endpoint_certificate().unwrap().endpoint, "5");
    report.verify().unwrap();

    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, b(&[&[1]])),
            RecurrenceTerm::new(1, 1, b(&[&[0, 1]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let report =
        check_ma_wang_interval_theorem(&rule, 0, &[vec![q(1), q(1)]], &Default::default()).unwrap();
    assert!(report.is_proved(), "{report:?}");
    assert_eq!(
        report.ma_wang_endpoint_certificate().unwrap().degree_step,
        0
    );
    report.verify().unwrap();
}

#[test]
fn endpoint_rejects_invalid_hypotheses_and_tampering() {
    let seeds = vec![vec![q(1)]];
    let mut bad = endpoint(3, 1, 3, 0);
    bad.terms[0].coeff.coeffs[0][0] = q(5) / q(2);
    assert!(
        !check_ma_wang_interval_theorem(&bad, 0, &seeds, &Default::default())
            .unwrap()
            .is_proved()
    );
    let mut bad = endpoint(3, 1, -3, 0);
    assert!(
        !check_ma_wang_interval_theorem(&bad, 0, &seeds, &Default::default())
            .unwrap()
            .is_proved()
    );
    bad = endpoint(3, -1, 3, 0);
    assert!(
        !check_ma_wang_interval_theorem(&bad, 0, &seeds, &Default::default())
            .unwrap()
            .is_proved()
    );
    let rule = endpoint(3, 1, 3, 0);
    assert!(
        !check_ma_wang_interval_theorem(&rule, 0, &[vec![q(0)]], &Default::default())
            .unwrap()
            .is_proved()
    );
    let report = check_ma_wang_interval_theorem(&rule, 0, &seeds, &Default::default()).unwrap();
    for field in [
        "endpoint",
        "initial_multiplicity",
        "maximum_endpoint_hits",
        "maximum_hit_scan",
    ] {
        let mut json = serde_json::to_value(&report).unwrap();
        let slot = &mut json["certificate"]["ma_wang_interval"][field];
        *slot = if field == "endpoint" {
            serde_json::json!("1")
        } else {
            serde_json::json!(999999)
        };
        let bad: RealRootednessTheoremReport = serde_json::from_value(json).unwrap();
        assert!(bad.verify().is_err(), "{field}");
    }
    let mut json = serde_json::to_value(report).unwrap();
    json.as_object_mut().unwrap().remove("certificate");
    assert!(serde_json::from_value::<RealRootednessTheoremReport>(json)
        .unwrap()
        .verify()
        .is_err());
}

fn gamma_input() -> RecurrenceJson {
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, b(&[&[1, 4, 1], &[0, -1]])),
            RecurrenceTerm::new(2, 0, b(&[&[0, -1, -2, -1], &[0, 1, 2, 1]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    RecurrenceJson::from_recurrence_rational(&rule, 0, &[vec![q(1)], vec![q(1), q(3), q(1)]], None)
}

#[test]
fn forward_gamma_proves_and_replays_the_transported_family() {
    let record = gamma_input();
    let report = discover_recurrence_proofs(&record, &Default::default()).unwrap();
    assert_eq!(report.status, "proved");
    let closure = report.closure_lifts.as_ref().unwrap();
    assert_eq!(closure.attempts[0].stage, "forward_gamma");
    report.verify().unwrap();
    for pointer in [
        "/closure_lifts/attempts/0/certificate/certification/gamma_transport/conditions/0/status",
        "/closure_lifts/attempts/0/parent_proof/proof/statement",
    ] {
        let mut json = serde_json::to_value(&report).unwrap();
        *json.pointer_mut(pointer).unwrap() = serde_json::json!("changed");
        let bad: RecurrenceProofDiscoveryReport = serde_json::from_value(json).unwrap();
        assert!(bad.verify().is_err(), "{pointer}");
    }
}

#[test]
fn forward_gamma_checks_all_index_symmetry_and_budgets() {
    let mut record = gamma_input();
    let disabled = find_discovery_closure_lifts(
        &record,
        &Default::default(),
        &CatalogLiftOptions::disabled(),
    )
    .unwrap();
    assert!(!disabled.is_proved());
    // Both seeds remain palindromic but the changed rule breaks symmetry.
    record.recurrence.terms[0].coeff.coeffs[0][0] = "2".into();
    let search =
        find_discovery_closure_lifts(&record, &Default::default(), &Default::default()).unwrap();
    assert!(!search.attempts.iter().any(|a| a.stage == "forward_gamma"));
}

#[test]
fn forward_gamma_removes_a_common_power_and_rejects_negative_gamma() {
    let record = gamma_input();
    let transformed = transform_recurrence(&RecurrenceTransformRequest::from_source(
        RecurrenceTransformation::MultiplyPower { power: 2 },
        record,
        12,
    ))
    .unwrap();
    let shifted = transformed.transported.unwrap().scalar.unwrap();
    let search =
        find_discovery_closure_lifts(&shifted, &Default::default(), &Default::default()).unwrap();
    assert!(search.is_proved());
    let DiscoveryClosureCertificate::ForwardGamma { certification } =
        &search.attempts[0].certificate
    else {
        panic!("expected forward gamma");
    };
    assert!(certification.power_removal.is_some());
    search.verify(&Default::default()).unwrap();

    let mut negative = gamma_input();
    negative.initial_polynomials[1] = vec!["1".into(), "1".into(), "1".into()];
    let search =
        find_discovery_closure_lifts(&negative, &Default::default(), &Default::default()).unwrap();
    assert!(!search.attempts.iter().any(|a| a.stage == "forward_gamma"));
}
