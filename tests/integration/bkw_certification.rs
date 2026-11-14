//! Exact BKW pencils, theorem proofs, degeneracies and finite counterexamples.

use polytool::bkw::certified::{
    certify_real_rooted_pencil, certify_recurrence, CertificationOptions, LimitSetStatus,
};
use polytool::recurrence::RecurrenceJson;
use serde_json::json;

fn record(terms: &[(usize, &[i64])], seeds: &[&[i64]]) -> RecurrenceJson {
    serde_json::from_value(json!({"schema":"polytool.recurrence.v1","first_index":0,
        "initial_polynomials":seeds.iter().map(|r|r.iter().map(ToString::to_string).collect::<Vec<_>>()).collect::<Vec<_>>(),
        "recurrence":{"terms":terms.iter().map(|(offset,c)|json!({"offset":offset,"deriv_order":0,"sign":"none","coeff":{"coeffs":[c.iter().map(ToString::to_string).collect::<Vec<_>>()]}})).collect::<Vec<_>>()}
    })).unwrap()
}
fn options() -> CertificationOptions {
    CertificationOptions {
        maximum_row_index: 48,
        ..Default::default()
    }
}
fn strings(p: &[i64]) -> Vec<String> {
    p.iter().map(ToString::to_string).collect()
}

#[test]
fn quadratic_product_proves_and_positive_pencil_refutes() {
    let good = record(&[(1, &[0, 1]), (2, &[-1])], &[&[1], &[0, 1]]);
    let report = certify_recurrence(&good, &options()).unwrap();
    assert_eq!(report.result, "proved");
    assert_eq!(report.limit_set_status, LimitSetStatus::RealLimitSet);
    report.verify().unwrap();
    let bad = record(&[(1, &[0, 1]), (2, &[1])], &[&[1], &[0, 1]]);
    let report = certify_recurrence(&bad, &options()).unwrap();
    assert_eq!(report.limit_set_status, LimitSetStatus::NonrealLimitSet);
    assert_eq!(report.first_failing_row.as_ref().unwrap().index, 2);
    assert!(report.pencil.as_ref().unwrap().witness.is_some());
    report.verify().unwrap();
    let mut tampered = report;
    tampered
        .pencil
        .as_mut()
        .unwrap()
        .witness
        .as_mut()
        .unwrap()
        .parameter = "0".into();
    assert!(tampered.verify().is_err());
}

#[test]
fn pencils_handle_degree_drops_fixed_multiple_roots_and_irrational_boundaries() {
    // x^2( x^2-s ): the unreduced discriminant vanishes identically.
    let p =
        certify_real_rooted_pencil(&strings(&[0, 0, 0, 0, 1]), &strings(&[0, 0, -1]), "4").unwrap();
    assert!(p.all_real_rooted);
    // (1-s)x^2-1: roots escape to infinity at s=1.
    let p = certify_real_rooted_pencil(&strings(&[-1, 0, 1]), &strings(&[0, 0, -1]), "4").unwrap();
    assert!(!p.all_real_rooted);
    assert!(p.witness.is_some());
    // x^3-3x+s has an interior critical parameter s=2.
    let p = certify_real_rooted_pencil(&strings(&[0, -3, 0, 1]), &strings(&[1]), "4").unwrap();
    assert!(!p.all_real_rooted);
    assert!(!p.critical_intervals.is_empty());
    // A genuinely irrational critical parameter: x^3-6x+s, s=4sqrt(2).
    let p = certify_real_rooted_pencil(&strings(&[0, -6, 0, 1]), &strings(&[1]), "8").unwrap();
    assert!(!p.all_real_rooted);
    assert_eq!(p.critical_intervals.len(), 1);
}

#[test]
fn circle_theorem_includes_boundary_and_public_proofs_take_priority() {
    let cubic = record(
        &[(1, &[-3, 1]), (2, &[-3]), (3, &[-1])],
        &[&[1], &[-3, 1], &[6, -6, 1]],
    );
    let report = certify_recurrence(&cubic, &options()).unwrap();
    assert_eq!(report.result, "proved");
    assert!(report.conditional_on.is_empty());
    assert_eq!(report.limit_set_status, LimitSetStatus::RealLimitSet);
    assert_eq!(
        report.toeplitz_region.as_ref().unwrap().limit_real,
        Some(true)
    );
    assert_eq!(
        report.proof_route.as_deref(),
        Some("hessenberg_toeplitz_circle")
    );
    let circle = report.circle.as_ref().unwrap();
    assert_eq!(circle.radius, "1");
    assert_eq!(circle.polynomial, strings(&[2, 2]));
    assert_eq!(
        circle.provenance,
        "circle criterion for Hessenberg Toeplitz symbols (unpublished manuscript)"
    );
    assert_eq!(
        circle.conclusion,
        "every row has simple real zeros and consecutive rows strictly interlace"
    );
    assert_eq!(
        circle.lean_declarations,
        [
            "RealRooted.Toeplitz.Circle.strictInterl_of_circleNonneg",
            "RealRooted.Toeplitz.Circle.roots_card_of_circleNonneg",
        ]
    );
    let wire = serde_json::to_value(&report).unwrap();
    assert_eq!(
        wire["circle"]["lean_declarations"],
        json!(circle.lean_declarations)
    );
    let round_trip: polytool::bkw::certified::CertificationReport =
        serde_json::from_value(wire.clone()).unwrap();
    assert_eq!(round_trip, report);
    round_trip.verify().unwrap();
    let mut legacy = wire;
    legacy["circle"]
        .as_object_mut()
        .unwrap()
        .remove("lean_declarations");
    let legacy: polytool::bkw::certified::CertificationReport =
        serde_json::from_value(legacy).unwrap();
    assert!(legacy.circle.as_ref().unwrap().lean_declarations.is_empty());
    assert!(serde_json::to_value(&legacy).unwrap()["circle"]
        .get("lean_declarations")
        .is_none());
    assert!(legacy.verify().is_err());
    report.verify().unwrap();
    let mut altered = report.clone();
    altered.circle.as_mut().unwrap().lean_declarations[0] = "RealRooted.forged".into();
    assert!(altered.verify().is_err());
    let mut altered = report.clone();
    altered.result = "conditional".into();
    assert!(altered.verify().is_err());
    let mut altered = report;
    altered.circle.as_mut().unwrap().provenance = "published theorem".into();
    assert!(altered.verify().is_err());
    let mut shifted = cubic;
    shifted.initial_polynomials[2][0] = "7".into();
    let report = certify_recurrence(&shifted, &options()).unwrap();
    assert!(
        report.circle.is_none(),
        "every extra seed must match the standard expansion"
    );
    let chebyshev = record(&[(1, &[1, -2]), (2, &[-1])], &[&[1], &[1, -2]]);
    let report = certify_recurrence(&chebyshev, &options()).unwrap();
    assert_eq!(report.result, "proved");
    assert_eq!(
        report.proof_route.as_deref(),
        Some("quadratic_standard_seed_product")
    );
    assert!(report.circle.is_none());
    assert!(report.conditional_on.is_empty());
}

#[test]
fn degeneracies_precede_asymptotic_claims() {
    let cases = [
        // Repeated pole (1-xt)^2.
        record(&[(1, &[0, 2]), (2, &[0, 0, -1])], &[&[1], &[0, 2]]),
        // Constant phase, scalar times a monomial.
        record(&[(1, &[0, 1]), (2, &[0, 0, -1])], &[&[1], &[0, 1]]),
        // The three even/odd fixtures, with every exceptional seed retained.
        record(
            &[(2, &[2, 0, -1]), (4, &[-1])],
            &[&[1], &[-1, 1], &[1, -1, -1], &[-1, 2, 1, -1]],
        ),
        record(
            &[(2, &[-3, 0, 1]), (4, &[-1])],
            &[
                &[-1],
                &[-1, -1],
                &[-2, 0, 1],
                &[3, 3, -1, -1],
                &[5, 0, -5, 0, 1],
            ],
        ),
        record(
            &[(2, &[-3, 0, 1]), (4, &[-1])],
            &[&[1], &[-1, -1], &[0, 2, 1], &[1, 3, -1, -1]],
        ),
        record(&[(1, &[0, 1])], &[&[0]]),
    ];
    for r in cases {
        let report = certify_recurrence(&r, &options()).unwrap();
        assert_eq!(
            report.limit_set_status,
            LimitSetStatus::Degenerate,
            "{}",
            report.reason
        );
        report.verify().unwrap();
    }
    // Cancellation removes the inactive branch; the existing product route
    // proves the reduced x^n sequence and retains its exact proof report.
    let r = record(&[(1, &[1, 1]), (2, &[0, -1])], &[&[1], &[0, 1]]);
    let report = certify_recurrence(&r, &options()).unwrap();
    assert_eq!(report.result, "proved");
    assert!(report.product_proof.is_some());
    assert_eq!(
        report
            .generating_function
            .as_ref()
            .unwrap()
            .denominator
            .len(),
        2
    );
    report.verify().unwrap();
}

#[test]
fn exact_first_failures_at_eleven_and_seventeen() {
    let quartic = record(
        &[(3, &[0, 0, 0, 1]), (4, &[0, 0, 0, -1])],
        &[&[1], &[0, 1], &[0, 0, 1], &[0, 0, 0, 1]],
    );
    let report = certify_recurrence(&quartic, &options()).unwrap();
    assert_eq!(report.limit_set_status, LimitSetStatus::Inconclusive);
    let failure = report.first_failing_row.unwrap();
    assert_eq!(failure.index, 11);
    assert_eq!(
        failure.coefficients,
        strings(&[0, 0, 0, 0, 0, 0, 0, 0, 0, 3, -3, 1])
    );
    let quartic = record(
        &[
            (1, &[4, 1]),
            (2, &[-6, -1]),
            (3, &[4, -1, 1]),
            (4, &[-1, 1]),
        ],
        &[&[0], &[0, 1], &[0, 1, 1], &[0, 1, 4, 1], &[0, 1, 9, 8, 1]],
    );
    let report = certify_recurrence(&quartic, &options()).unwrap();
    assert_eq!(report.result, "refuted");
    let failure = report.first_failing_row.as_ref().unwrap();
    assert_eq!(failure.index, 17);
    assert_eq!(failure.root_counts.distinct_real_roots, 15);
    assert_eq!(failure.root_counts.squarefree_degree, 17);
    report.verify().unwrap();
}

#[test]
fn sparse_cubic_uses_the_tran_interval_without_claiming_a_proof() {
    let good = record(
        &[(1, &[0, 1]), (3, &[0, 0, -1])],
        &[&[1], &[0, 1], &[0, 0, 1]],
    );
    let report = certify_recurrence(&good, &options()).unwrap();
    assert_eq!(report.pencil.as_ref().unwrap().interval, ["0", "27/4"]);
    assert_eq!(report.limit_set_status, LimitSetStatus::RealLimitSet);
    assert_ne!(report.result, "proved");
    let bad = record(&[(1, &[1, 1]), (3, &[1])], &[&[1], &[1, 1], &[1, 2, 1]]);
    assert_eq!(
        certify_recurrence(&bad, &options())
            .unwrap()
            .limit_set_status,
        LimitSetStatus::NonrealLimitSet
    );
    // Even sparse order uses the same +sB sign, with endpoint 4^4/3^3.
    let fourth = record(
        &[(1, &[0, 1]), (4, &[-1])],
        &[&[1], &[0, 1], &[0, 0, 1], &[0, 0, 0, 1]],
    );
    let report = certify_recurrence(&fourth, &options()).unwrap();
    assert_eq!(report.pencil.as_ref().unwrap().interval, ["0", "256/27"]);
    assert_eq!(report.limit_set_status, LimitSetStatus::NonrealLimitSet);
    assert_eq!(report.first_failing_row.as_ref().unwrap().index, 4);
}

#[test]
fn cli_round_trip_and_invalid_radius() {
    let r = record(&[(1, &[0, 1]), (2, &[-1])], &[&[1], &[0, 1]]);
    let input = serde_json::to_string(&r).unwrap();
    let output = crate::common::run_polytool(&["bkw", "certify", "--json"], &input);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = String::from_utf8(output.stdout).unwrap();
    assert_eq!(
        serde_json::from_str::<serde_json::Value>(&text).unwrap()["result"],
        "proved"
    );
    assert!(
        crate::common::run_polytool(&["bkw", "certify", "--verify", "--json"], &text)
            .status
            .success()
    );
    assert_eq!(
        crate::common::run_polytool(&["bkw", "certify", "--radius", "0"], &input)
            .status
            .code(),
        Some(2)
    );
}

#[test]
fn small_shape_samples_agree_with_the_exact_regions() {
    // S1: A=x, M=e and H=d^2-(1+4c)e.
    for (c, d, e, good) in [
        (0, 0, -1, true),
        (-1, 0, -1, false),
        (1, 0, 1, false),
        (1, 3, -1, true),
        (0, 1, 0, true),
    ] {
        let r = record(&[(1, &[0, 1]), (2, &[e, d, c])], &[&[1], &[0, 1]]);
        let report = certify_recurrence(&r, &options()).unwrap();
        assert_eq!(report.pencil.as_ref().unwrap().all_real_rooted, good);
        assert_eq!(report.result == "proved", good);
    }
    // S2: A=x+1, Q=cx+d; a(ad-bc)=d-c.
    for (c, d, good) in [(1, 0, true), (1, 2, false), (-1, -2, true)] {
        let r = record(&[(1, &[1, 1]), (2, &[d, c])], &[&[1], &[1, 1]]);
        assert_eq!(
            certify_recurrence(&r, &options()).unwrap().result == "proved",
            good
        );
    }
}

#[test]
fn discovery_refutation_is_late_replayable_and_never_replaces_a_proof() {
    use polytool::recurrence::{
        discover_recurrence_proofs_with_stage_options, RealRootednessTheoremOptions,
        RecurrenceProofStageOptions,
    };
    let bad = record(
        &[(3, &[0, 0, 0, 1]), (4, &[0, 0, 0, -1])],
        &[&[1], &[0, 1], &[0, 0, 1], &[0, 0, 0, 1]],
    );
    let mut stages = RecurrenceProofStageOptions::default();
    stages.alternative_representations.max_candidates = 0;
    stages.catalog_lifts.max_candidates = 0;
    stages.affine_coordinates.max_coordinates = 0;
    let mut discovery = discover_recurrence_proofs_with_stage_options(
        &bad,
        &RealRootednessTheoremOptions {
            maximum_start_delay: 0,
            ..Default::default()
        },
        &stages,
    )
    .unwrap();
    assert_ne!(discovery.status, "proved");
    discovery.refute_if_possible().unwrap();
    assert_eq!(discovery.status, "refuted");
    discovery.verify().unwrap();
    let mut changed = discovery.clone();
    changed
        .bkw_refutation
        .as_mut()
        .unwrap()
        .first_failing_row
        .as_mut()
        .unwrap()
        .index = 12;
    assert!(changed.verify().is_err());
    let good = record(&[(1, &[0, 1]), (2, &[-1])], &[&[1], &[0, 1]]);
    let mut discovery =
        discover_recurrence_proofs_with_stage_options(&good, &Default::default(), &stages).unwrap();
    assert_eq!(discovery.status, "proved");
    discovery.refute_if_possible().unwrap();
    assert!(discovery.bkw_refutation.is_none());
    discovery.verify().unwrap();
    let input = serde_json::to_string(&bad).unwrap();
    let output = crate::common::run_polytool(
        &[
            "recurrence-prove",
            "--max-start-delay",
            "0",
            "--max-representation-candidates",
            "0",
            "--max-lift-candidates",
            "0",
            "--max-affine-coordinates",
            "0",
            "--json",
        ],
        &input,
    );
    assert_eq!(output.status.code(), Some(1));
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["status"], "refuted");
    assert_eq!(report["bkw_refutation"]["first_failing_row"]["index"], 11);
    let output = crate::common::run_polytool(
        &["recurrence-prove", "--verify", "--json"],
        &String::from_utf8(output.stdout).unwrap(),
    );
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}
