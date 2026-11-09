//! Factor gauge with a bounded residual (the RealRooted row routes
//! `powerForm`, `closedForm` and `factorPower`): catalog families, every
//! residual shape, negative cases, tampered certificates, discovery and the
//! parent provers.

use polytool::recurrence::*;

fn q(n: i64) -> BigRational {
    BigRational::from_integer(n.into())
}

fn b(rows: &[&[i64]]) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(
        rows.iter()
            .map(|r| r.iter().map(|&n| q(n)).collect())
            .collect(),
    )
}

fn rows(seeds: &[&[i64]]) -> Vec<Vec<BigRational>> {
    seeds
        .iter()
        .map(|r| r.iter().map(|&n| q(n)).collect())
        .collect()
}

fn options() -> RealRootednessTheoremOptions {
    RealRootednessTheoremOptions::default()
}

/// `P_n = sum_j c_j(x) P_(n-j)` with constant (in `n`) multipliers, given
/// as ascending coefficient lists in `x`, lag `1` first.
fn constant_rule(multipliers: &[&[i64]]) -> Recurrence {
    Recurrence {
        terms: multipliers
            .iter()
            .enumerate()
            .map(|(lag, c)| RecurrenceTerm::new(lag + 1, 0, b(&[c])))
            .collect(),
        denominator: None,
        inhomogeneous: None,
    }
}

fn catalog(id: &str) -> (Recurrence, Vec<Vec<BigRational>>, usize) {
    let entry = polytool::oeis::by_id(id).unwrap();
    let (rule, rows) = entry.recurrence_parts().unwrap();
    (rule, rows, entry.recurrence_first_index)
}

fn check(id: &str) -> RealRootednessTheoremReport {
    let (rule, seeds, first) = catalog(id);
    check_factor_gauge_theorem(&rule, first, &seeds, &options()).unwrap()
}

fn hypothesis<'a>(report: &'a RealRootednessTheoremReport, name: &str) -> &'a LiftHypothesis {
    report
        .hypotheses
        .iter()
        .find(|h| h.name == name)
        .unwrap_or_else(|| panic!("no hypothesis {name}"))
}

/// An edit of a report's JSON.
type Edit = Box<dyn Fn(&mut serde_json::Value)>;

fn tampered(report: &RealRootednessTheoremReport, edit: impl Fn(&mut serde_json::Value)) -> bool {
    let mut value = serde_json::to_value(report).unwrap();
    edit(&mut value);
    let report: RealRootednessTheoremReport = serde_json::from_value(value).unwrap();
    report.verify().is_err()
}

// ---------------------------------------------------------------------------
// Proved families
// ---------------------------------------------------------------------------

#[test]
fn the_catalog_targets_are_proved() {
    // (id, residual shape, gauge factor, certified region closed at 0)
    let expected = [
        (
            "A124376",
            FactorGaugeResidualShape::PalindromicQuartic,
            vec!["1", "1"],
            Some(false),
        ),
        (
            "A157169",
            FactorGaugeResidualShape::PalindromicQuartic,
            vec!["1", "1"],
            Some(false),
        ),
        (
            "A157171",
            FactorGaugeResidualShape::PalindromicQuartic,
            vec!["1", "1"],
            Some(false),
        ),
        (
            "A110200",
            FactorGaugeResidualShape::Linear,
            vec!["1", "1"],
            None,
        ),
        (
            "A133109",
            FactorGaugeResidualShape::Constant,
            vec!["0", "1"],
            Some(true),
        ),
    ];
    for (id, shape, gauge, closed) in expected {
        let report = check(id);
        assert!(report.is_proved(), "{id}: {:?}", report.open_hypotheses());
        let certificate = report.factor_gauge().unwrap();
        assert_eq!(certificate.residual_shape, shape, "{id}");
        assert_eq!(certificate.gauge_factor, gauge, "{id}");
        let conclusion = report.conclusion.as_ref().unwrap();
        assert_eq!(
            conclusion.first_index, report.recurrence.first_index,
            "{id}"
        );
        let region = conclusion.root_region.as_ref().map(|region| {
            let region = region.at(0).unwrap();
            assert!(region.lower.is_none());
            assert_eq!(region.upper, Some(q(0)));
            region.upper_closed
        });
        assert_eq!(region, closed, "{id}");
        report.verify().unwrap();
    }
}

#[test]
fn palindromic_quartics_use_the_u_substitution() {
    let report = check("A124376");
    for name in [
        "residual_outer_coefficient_sign",
        "residual_u_vertex_at_most_minus_two",
        "residual_u_value_at_minus_two_nonnegative",
        "residual_u_discriminant_nonnegative",
    ] {
        assert_eq!(
            hypothesis(&report, name).status,
            LiftHypothesisStatus::Verified,
            "{name}"
        );
    }
    assert_eq!(report.lean_route.as_deref(), Some("closedForm"));
    assert_eq!(
        report.lean_declarations,
        vec!["RealRooted.splits_palindromic_quartic"]
    );
    let certificate = report.factor_gauge().unwrap();
    // P_n = (1 + x)^(n - 4) S_n, S_n = a x^4 + b x^3 + c x^2 + b x + a with
    // a = 1, b = 2n + 2, c = (n^2 + 3n + 10)/2.
    assert_eq!(certificate.gauge_index, 4);
    assert_eq!(certificate.common_factor, vec!["1"]);
    assert_eq!(
        certificate.residual_coefficients.as_ref().unwrap(),
        &vec![
            vec!["1".to_string()],
            vec!["2".into(), "2".into()],
            vec!["5".into(), "3/2".into(), "1/2".into()],
            vec!["2".into(), "2".into()],
            vec!["1".into()],
        ]
    );
}

#[test]
fn exponential_coordinates_prove_nonzero_rows_by_a_dominant_root() {
    // A110200 (catalog first index 0): P_n = (1+x)^(n-2) (1+x) S_n with
    // S_n = (4^(n+1) - 1)/3 + (2^(n+1) - 1)^2 x.
    let report = check("A110200");
    let certificate = report.factor_gauge().unwrap();
    assert_eq!(certificate.common_factor, vec!["1", "1"]);
    assert!(certificate.residual_coefficients.is_none());
    let dominant = certificate.dominant_root.as_ref().unwrap();
    assert_eq!(dominant.coefficient, 1);
    assert_eq!(dominant.roots, vec!["4", "2", "1"]);
    assert_eq!(dominant.weights, vec!["4", "-4", "1"]);
    // No region is claimed for exponential coordinates.
    assert!(report.conclusion.as_ref().unwrap().root_region.is_none());
    let (rule, seeds, first) = catalog("A110200");
    certificate.verify(&rule, first, &seeds).unwrap();
    // A forged weight or threshold does not replay.
    let mut forged = certificate.clone();
    forged.dominant_root.as_mut().unwrap().weights[1] = "-3".into();
    assert!(forged.verify(&rule, first, &seeds).is_err());
    let mut forged = certificate.clone();
    forged.dominant_root.as_mut().unwrap().roots.swap(0, 1);
    assert!(forged.verify(&rule, first, &seeds).is_err());
}

#[test]
fn constant_residuals_with_a_linear_gauge_interlace() {
    // A133109: (n - 2) P_n = x P_(n-1) + x^2 P_(n-2) + (n - 1) x^3 P_(n-3),
    // monomial rows c_n x^(n-1).
    let report = check("A133109");
    let conclusion = report.conclusion.as_ref().unwrap();
    assert_eq!(conclusion.interlacing, Some(RowInterlacing::Weak));
    assert_eq!(conclusion.interlacing_from, Some(3));
    assert_eq!(report.lean_route.as_deref(), Some("powerForm"));
    assert!(report
        .lean_declarations
        .contains(&"RealRooted.interlaces_of_forall_eq_C_mul_pow".to_string()));
    assert_eq!(
        hypothesis(&report, "left_factor_nonzero").status,
        LiftHypothesisStatus::Verified
    );
    // The scalar coordinates are not polynomial: the rows are nonzero by the
    // shared leading-coefficient certificate.
    assert!(hypothesis(&report, "rows_nonzero")
        .evidence
        .contains("leading-coefficient"));
}

#[test]
fn the_lean_named_rows_are_proved_by_the_gauge() {
    for id in [
        "A106246", "A124932", "A125092", "A128540", "A141611", "A127647",
    ] {
        let report = check(id);
        assert!(report.is_proved(), "{id}: {:?}", report.open_hypotheses());
        report.verify().unwrap();
    }
    // A332496 starts with the zero row P_0 = 0: the record fails exactly at
    // that seed, and its tail from P_2 on is proved.
    let report = check("A332496");
    assert_eq!(report.status, RealRootednessTheoremStatus::HypothesisFailed);
    assert_eq!(
        hypothesis(&report, "earlier_seeds").status,
        LiftHypothesisStatus::Failed {
            witness: Some(0.into())
        }
    );
    let (rule, seeds, first) = catalog("A332496");
    let tail = check_factor_gauge_theorem(&rule, first + 1, &seeds[1..], &options()).unwrap();
    assert!(tail.is_proved(), "{:?}", tail.open_hypotheses());
}

#[test]
fn quadratic_and_cubic_pencils_are_decided_by_discriminants() {
    // P_n = 2 P_(n-1) - P_(n-2): P_n = A + n B (characteristic root 1,
    // twice), q = 1.
    let rule = constant_rule(&[&[2], &[-1]]);
    // A = (x+1)(x+3)(x+5), B = (x+2)(x+4)(x+6): interlacing, so A + n B is
    // real-rooted for n >= 0; P_1 = A + B.
    let a = [15, 23, 9, 1];
    let b_ = [48, 44, 12, 1];
    let p1: Vec<i64> = a.iter().zip(&b_).map(|(x, y)| x + y).collect();
    let seeds = rows(&[&a, &p1]);
    let report = check_factor_gauge_theorem(&rule, 0, &seeds, &options()).unwrap();
    assert!(report.is_proved(), "{:?}", report.open_hypotheses());
    let certificate = report.factor_gauge().unwrap();
    assert_eq!(certificate.residual_shape, FactorGaugeResidualShape::Cubic);
    assert_eq!(certificate.gauge_factor, vec!["1"]);
    assert_eq!(
        hypothesis(&report, "residual_cubic_discriminant_nonnegative").status,
        LiftHypothesisStatus::Verified
    );
    assert!(report.lean_route.is_none());
    report.verify().unwrap();

    // P_0 = x^2 - 1, P_1 = x^2: P_2 = x^2 + 1 is not real-rooted, and the
    // discriminant -4(n - 1) fails exactly at n = 2.
    let seeds = rows(&[&[-1, 0, 1], &[0, 0, 1]]);
    let report = check_factor_gauge_theorem(
        &rule,
        0,
        &seeds,
        &RealRootednessTheoremOptions {
            maximum_start_delay: 0,
            ..options()
        },
    )
    .unwrap();
    assert_eq!(report.status, RealRootednessTheoremStatus::HypothesisFailed);
    assert_eq!(
        hypothesis(&report, "residual_discriminant_nonnegative").status,
        LiftHypothesisStatus::Failed {
            witness: Some(2.into())
        }
    );
    report.verify().unwrap();
}

#[test]
fn a_requested_region_is_certified_or_left_open() {
    let (rule, seeds, first) = catalog("A124376");
    let negative_axis = RealRootednessTheoremOptions {
        root_region: Some(polytool::real_rootedness::RealInterval {
            lower: None,
            lower_closed: false,
            upper: Some(q(0)),
            upper_closed: false,
        }),
        ..options()
    };
    let report = check_factor_gauge_theorem(&rule, first, &seeds, &negative_axis).unwrap();
    assert!(report.is_proved());
    assert_eq!(
        hypothesis(&report, "requested_root_region").status,
        LiftHypothesisStatus::Verified
    );
    report.verify().unwrap();
    // A133109 has roots at 0: the open request is not certified.
    let (rule, seeds, first) = catalog("A133109");
    let report = check_factor_gauge_theorem(&rule, first, &seeds, &negative_axis).unwrap();
    assert_eq!(report.status, RealRootednessTheoremStatus::Unresolved);
    assert_eq!(
        hypothesis(&report, "requested_root_region").status,
        LiftHypothesisStatus::Unresolved
    );
}

// ---------------------------------------------------------------------------
// Negative cases
// ---------------------------------------------------------------------------

#[test]
fn rules_without_a_common_gauge_are_not_applicable() {
    let not_applicable = |rule: &Recurrence, seeds: &[Vec<BigRational>]| {
        let report = check_factor_gauge_theorem(rule, 0, seeds, &options()).unwrap();
        assert_eq!(
            report.status,
            RealRootednessTheoremStatus::NotApplicable,
            "{}",
            report.normal_form
        );
        assert!(report.factor_gauge().is_none());
        report.normal_form
    };
    // Chebyshev: 2x P_(n-1) - P_(n-2), q = x from lag 1 but -1 != x^2.
    let chebyshev = constant_rule(&[&[0, 2], &[-1]]);
    let reason = not_applicable(&chebyshev, &rows(&[&[1], &[0, 1]]));
    assert!(reason.contains("for one q"), "{reason}");
    // A first-order rule is the product route's case.
    let product = constant_rule(&[&[1, 1]]);
    assert!(not_applicable(&product, &rows(&[&[1]])).contains("first-order"));
    // A derivative term.
    let derivative = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 1, b(&[&[0, 1]])),
            RecurrenceTerm::new(2, 0, b(&[&[0, 0, 1]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    assert!(not_applicable(&derivative, &rows(&[&[1], &[1, 1]])).contains("derivative"));
    // The lag-2 multiplier x is not a square.
    let not_square = constant_rule(&[&[0], &[0, 1]]);
    let reason = not_applicable(&not_square, &rows(&[&[1], &[1, 1]]));
    assert!(reason.contains("power"), "{reason}");
    // A quartic residual that is not palindromic.
    let pencil = constant_rule(&[&[2], &[-1]]);
    let reason = not_applicable(&pencil, &rows(&[&[1, 2, 3, 4, 5], &[1, 0, 0, 0, 1]]));
    assert!(reason.contains("palindromic"), "{reason}");
    // A quadratic residual with exponential coordinates (roots 1 and 2).
    let exponential = constant_rule(&[&[3], &[-2]]);
    let reason = not_applicable(&exponential, &rows(&[&[1, 3, 1], &[2, 5, 1]]));
    assert!(reason.contains("not polynomials"), "{reason}");
}

#[test]
fn a_gauge_factor_that_is_not_real_rooted_fails() {
    // P_n = 2 (1 + x^2) P_(n-1) - (1 + x^2)^2 P_(n-2).
    let rule = constant_rule(&[&[2, 0, 2], &[-1, 0, -2, 0, -1]]);
    let report = check_factor_gauge_theorem(&rule, 0, &rows(&[&[1], &[1, 1]]), &options()).unwrap();
    assert_eq!(report.status, RealRootednessTheoremStatus::HypothesisFailed);
    assert!(matches!(
        hypothesis(&report, "gauge_factor_real_rooted").status,
        LiftHypothesisStatus::Failed { .. }
    ));
    report.verify().unwrap();
}

#[test]
fn zero_rows_fail_exactly() {
    // P_n = x^2 P_(n-2), P_0 = 1, P_1 = 0: every odd row is zero.
    let rule = constant_rule(&[&[0], &[0, 0, 1]]);
    let report = check_factor_gauge_theorem(&rule, 0, &rows(&[&[1], &[0]]), &options()).unwrap();
    assert_eq!(report.status, RealRootednessTheoremStatus::HypothesisFailed);
    assert_eq!(
        hypothesis(&report, "rows_nonzero").status,
        LiftHypothesisStatus::Failed {
            witness: Some(3.into())
        }
    );
}

// ---------------------------------------------------------------------------
// Tampering
// ---------------------------------------------------------------------------

#[test]
fn tampered_reports_do_not_replay() {
    let report = check("A124376");
    report.verify().unwrap();
    let edits: Vec<Edit> = vec![
        Box::new(|v| v["factor_gauge"]["gauge_factor"] = serde_json::json!(["2", "1"])),
        Box::new(|v| v["factor_gauge"]["common_factor"] = serde_json::json!(["1", "1"])),
        Box::new(|v| v["factor_gauge"]["residual_seeds"][0][1] = serde_json::json!("7")),
        Box::new(|v| v["factor_gauge"]["residual_coefficients"][2][0] = serde_json::json!("6")),
        Box::new(|v| v["factor_gauge"]["residual_shape"] = serde_json::json!("quadratic")),
        Box::new(|v| v["factor_gauge"]["scalar_multipliers"][0] = serde_json::json!(["4"])),
        Box::new(|v| v["factor_gauge"]["gauge_index"] = serde_json::json!(3)),
        Box::new(|v| {
            v["factor_gauge"]
                .as_object_mut()
                .unwrap()
                .remove("residual_coefficients");
        }),
        Box::new(|v| {
            v.as_object_mut().unwrap().remove("factor_gauge");
        }),
        Box::new(|v| v["hypotheses"][1]["evidence"] = serde_json::json!("trust me")),
        Box::new(|v| {
            v["conclusion"]["interlacing"] = serde_json::json!("weak");
            v["conclusion"]["interlacing_from"] = serde_json::json!(5);
        }),
        Box::new(|v| v["theorem"] = serde_json::json!("index_section_lift")),
    ];
    for (position, edit) in edits.iter().enumerate() {
        assert!(tampered(&report, edit), "edit {position} replays");
    }
    // The certificate alone: other seeds do not replay.
    let (rule, mut seeds, first) = catalog("A124376");
    let certificate = report.factor_gauge().unwrap();
    certificate.verify(&rule, first, &seeds).unwrap();
    seeds[3][1] = q(8);
    assert!(certificate.verify(&rule, first, &seeds).is_err());
}

// ---------------------------------------------------------------------------
// Discovery and the parent provers
// ---------------------------------------------------------------------------

#[test]
fn discovery_records_the_gauge_only_where_it_applies() {
    let (rule, seeds, first) = catalog("A124376");
    let record = RecurrenceJson::from_recurrence_rational(&rule, first, &seeds, None);
    let report = discover_recurrence_proofs_with_options(
        &record,
        &options(),
        &AlternativeRepresentationOptions {
            max_candidates: 0,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(report.status, "proved");
    assert_eq!(report.real_rooted_from_index, Some(0));
    let attempt = report.attempts.last().unwrap();
    assert_eq!(attempt.theorem, RealRootednessTheorem::FactorGauge);
    assert!(attempt.is_proved());
    report.verify().unwrap();
    // A report saved without the attempt still replays.
    let mut older = report.clone();
    older.attempts.pop();
    older.status = "not_proved".into();
    older.real_rooted_from_index = None;
    older.statement = report
        .statement
        .split(" Proved through")
        .next()
        .unwrap()
        .into();
    older.verify().unwrap();

    // Chebyshev rows: no gauge, so no attempt is recorded.
    let rule = constant_rule(&[&[0, 2], &[-1]]);
    let record = RecurrenceJson::from_recurrence_rational(&rule, 0, &rows(&[&[1], &[0, 1]]), None);
    let report = discover_recurrence_proofs(&record, &options()).unwrap();
    assert!(report
        .attempts
        .iter()
        .all(|a| a.theorem != RealRootednessTheorem::FactorGauge));
}

#[test]
fn the_parent_provers_use_the_gauge_for_an_open_parent() {
    let (rule, seeds, first) = catalog("A124376");
    let record = ParentRecurrenceRecord {
        rule: &rule,
        seeds: &seeds,
        first_index: first,
    };
    let proof = prove_parent_real_rootedness(&record, first, &[], false).unwrap();
    assert_eq!(proof.status, VeroneseParentProofStatus::Proved);
    assert_eq!(proof.route, Some(VeroneseParentProofRoute::FactorGauge));
    assert!(proof.interlacing.is_none());
    check_recorded_parent_proof(&record, first, &[], false, &proof).unwrap();
}
