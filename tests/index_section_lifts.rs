//! Index-section lifts (the RealRooted row route `pathPower`): catalog
//! families that are (reversed) every-`p`-th rows of a constant-coefficient
//! parent, families whose parent stays open, negative cases, tampered
//! certificates, discovery and the parent provers.

use polytool::recurrence::*;
use serde_json::json;

mod common;
use common::{q, rows};

fn b(rows: &[&[i64]]) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(
        rows.iter()
            .map(|r| r.iter().map(|&n| q(n)).collect())
            .collect(),
    )
}

fn options() -> RealRootednessTheoremOptions {
    RealRootednessTheoremOptions::default()
}

/// The catalog family: rule, all initial rows of the entry (the rule applies
/// after them, as in its `recurrence_data`) and first index.
fn catalog(id: &str) -> (Recurrence, Vec<Vec<BigRational>>, usize) {
    let entry = polytool::oeis::by_id(id).unwrap();
    let (rule, rows) = entry.recurrence_parts().unwrap();
    (rule, rows, entry.recurrence_first_index)
}

fn check(id: &str) -> RealRootednessTheoremReport {
    let (rule, seeds, first) = catalog(id);
    check_index_section_lift_theorem(&rule, first, &seeds, &options()).unwrap()
}

fn hypothesis<'a>(report: &'a RealRootednessTheoremReport, name: &str) -> &'a LiftHypothesis {
    report
        .hypotheses
        .iter()
        .find(|h| h.name == name)
        .unwrap_or_else(|| panic!("no hypothesis {name}"))
}

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
fn the_path_power_catalog_rows_are_proved() {
    // (id, period, reversed, exceptional rows)
    let expected = [
        ("A109955", 3, true, 0),
        ("A109960", 4, true, 0),
        ("A116088", 2, true, 0),
        ("A116089", 3, true, 0),
        ("A127893", 3, true, 0),
        ("A193636", 3, false, 0),
        ("A206294", 3, true, 1),
    ];
    for (id, period, reversed, exceptional) in expected {
        let report = check(id);
        assert!(report.is_proved(), "{id}: {:?}", report.open_hypotheses());
        let certificate = report.index_section().unwrap();
        assert_eq!(certificate.period, period, "{id}");
        assert_eq!(
            matches!(
                certificate.transform,
                IndexSectionTransform::Reversal { .. }
            ),
            reversed,
            "{id}"
        );
        assert_eq!(certificate.exceptional_rows, exceptional, "{id}");
        let selected = &certificate.parent_attempts[certificate.selected.unwrap()];
        assert_eq!(
            selected.theorem,
            RealRootednessTheorem::VariableDelay,
            "{id}"
        );
        let conclusion = report.conclusion.as_ref().unwrap();
        assert_eq!(conclusion.first_index, 0, "{id}");
        assert!(conclusion.interlacing.is_none(), "{id}");
        assert_eq!(report.lean_route.as_deref(), Some("pathPower"));
        assert_eq!(
            report.lean_declarations[0],
            "LinearRecurrence.ofMatrix_pow_isSolution_subsequence"
        );
        report.verify().unwrap();
        let (rule, seeds, first) = catalog(id);
        certificate.verify(&rule, first, &seeds).unwrap();
    }
}

#[test]
fn a_reversed_section_has_a_region_closed_at_zero() {
    // A116088: every row P_n, n >= 1, has the factor x, so the transferred
    // region keeps the end 0.
    let report = check("A116088");
    let region = report
        .conclusion
        .as_ref()
        .unwrap()
        .root_region
        .as_ref()
        .unwrap();
    let region = region.at(0).unwrap();
    assert_eq!(region.lower, None);
    assert_eq!(region.upper, Some(q(0)));
    assert!(region.upper_closed);
    // A direct section keeps the parent region (-inf, 0).
    let report = check("A193636");
    let region = report
        .conclusion
        .as_ref()
        .unwrap()
        .root_region
        .as_ref()
        .unwrap();
    let region = region.at(0).unwrap();
    assert_eq!(region.upper, Some(q(0)));
    assert!(!region.upper_closed);
}

#[test]
fn the_parent_of_a109955_is_the_fixed_delay_window() {
    let report = check("A109955");
    let certificate = report.index_section().unwrap();
    let (parent, first, seeds) = certificate.parent.to_recurrence_parts().unwrap();
    assert_eq!(first, 0);
    assert_eq!(seeds, rows(&[&[1], &[1], &[1]]));
    // Q_n = Q_(n-1) + x Q_(n-3): the reversed rows are every third row.
    let parent_rows = parent.generate_rows_rational(&seeds, 0, 13).unwrap();
    let (rule, seeds, first) = catalog("A109955");
    let input = rule.generate_rows_rational(&seeds, first, 5).unwrap();
    for (m, row) in input.iter().enumerate() {
        let mut reversed = parent_rows[3 * m].clone();
        reversed.resize(m + 1, q(0));
        reversed.reverse();
        assert_eq!(&reversed, row, "m = {m}");
    }
}

#[test]
fn a_requested_region_is_checked_against_the_transferred_one() {
    let (rule, seeds, first) = catalog("A193636");
    let negative = RealRootednessTheoremOptions {
        root_region: Some(polytool::real_rootedness::RealInterval {
            lower: None,
            lower_closed: false,
            upper: Some(q(0)),
            upper_closed: false,
        }),
        ..options()
    };
    let report = check_index_section_lift_theorem(&rule, first, &seeds, &negative).unwrap();
    assert!(report.is_proved(), "{:?}", report.open_hypotheses());
    report.verify().unwrap();
    // The reversed section of A116088 has roots at 0: the open request stays
    // open.
    let (rule, seeds, first) = catalog("A116088");
    let report = check_index_section_lift_theorem(&rule, first, &seeds, &negative).unwrap();
    assert_eq!(report.status, RealRootednessTheoremStatus::Unresolved);
    assert_eq!(
        hypothesis(&report, "requested_root_region").status,
        LiftHypothesisStatus::Unresolved
    );
    report.verify().unwrap();
}

// ---------------------------------------------------------------------------
// Open parents
// ---------------------------------------------------------------------------

#[test]
fn the_cubic_parent_rows_have_a_section_identity_but_an_open_parent() {
    // The parent 1/(1 - t - x t^2 - x t^3) needs the Tran--Zumba theorem.
    for id in ["A114123", "A184883", "A216182"] {
        let report = check(id);
        assert!(!report.is_proved(), "{id}");
        assert_ne!(report.status, RealRootednessTheoremStatus::NotApplicable);
        assert_eq!(
            hypothesis(&report, "section_identity").status,
            LiftHypothesisStatus::Verified
        );
        let parent = hypothesis(&report, "parent_real_rooted");
        assert_ne!(parent.status, LiftHypothesisStatus::Verified);
        assert!(parent.evidence.starts_with("parent open"), "{id}");
        let certificate = report.index_section().unwrap();
        assert_eq!(certificate.period, 2);
        assert!(certificate.selected.is_none());
        report.verify().unwrap();
    }
}

// ---------------------------------------------------------------------------
// Negative cases
// ---------------------------------------------------------------------------

#[test]
fn a_parent_with_zero_residue_classes_proves_nothing_here() {
    // A132918: P_n = x^2 P_(n-2), P_1 = 1, P_2 = 0, so every other row is
    // zero; a parent proved only with zero residue classes must not give
    // "every row nonzero".
    let report = check("A132918");
    assert!(!report.is_proved());
    let certificate = report.index_section().unwrap();
    assert!(certificate.selected.is_none());
    report.verify().unwrap();
}

#[test]
fn index_dependent_and_derivative_rules_do_not_apply() {
    // P_n = (n + x) P_(n-1) - P_(n-2).
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, b(&[&[0, 1], &[1]])),
            RecurrenceTerm::new(2, 0, b(&[&[-1]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let report =
        check_index_section_lift_theorem(&rule, 0, &rows(&[&[1], &[0, 1]]), &options()).unwrap();
    assert_eq!(report.status, RealRootednessTheoremStatus::NotApplicable);
    assert!(
        report.normal_form.contains("depends on n"),
        "{}",
        report.normal_form
    );
    assert!(report.index_section().is_none());
    report.verify().unwrap();
    // P_n = (1 + x) P_(n-1) + x P'_(n-1).
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, b(&[&[1, 1]])),
            RecurrenceTerm::new(1, 1, b(&[&[0, 1]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let report = check_index_section_lift_theorem(&rule, 0, &rows(&[&[1]]), &options()).unwrap();
    assert_eq!(report.status, RealRootednessTheoremStatus::NotApplicable);
    assert!(
        report.normal_form.contains("derivative"),
        "{}",
        report.normal_form
    );
}

#[test]
fn a_rule_without_a_section_parent_does_not_apply() {
    // Chebyshev: P_n = 2x P_(n-1) - P_(n-2); chi(w^p) has no divisor of
    // degree 2 with the resultant identity.
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, b(&[&[0, 2]])),
            RecurrenceTerm::new(2, 0, b(&[&[-1]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let report =
        check_index_section_lift_theorem(&rule, 0, &rows(&[&[1], &[0, 1]]), &options()).unwrap();
    assert_eq!(report.status, RealRootednessTheoremStatus::NotApplicable);
    assert!(report.normal_form.contains("no index-section parent"));
}

// ---------------------------------------------------------------------------
// Tampering
// ---------------------------------------------------------------------------

#[test]
fn tampered_index_section_certificates_are_rejected() {
    let report = check("A109955");
    report.verify().unwrap();
    let (rule, seeds, first) = catalog("A109955");
    let certificate = report.index_section().cloned().unwrap();
    // Another period.
    let mut wrong = certificate.clone();
    wrong.period = 2;
    assert!(wrong.verify(&rule, first, &seeds).is_err());
    // Another parent rule with the same seeds (a wrong section claim):
    // Q_n = Q_(n-1) + 2x Q_(n-3) fails the resultant identity.
    let mut wrong = certificate.clone();
    let parent = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, b(&[&[1]])),
            RecurrenceTerm::new(3, 0, b(&[&[0, 2]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    wrong.parent =
        RecurrenceJson::from_recurrence_rational(&parent, 0, &rows(&[&[1], &[1], &[1]]), None);
    assert!(wrong.verify(&rule, first, &seeds).is_err());
    // Other parent seeds.
    let mut wrong = certificate.clone();
    let (parent, _, _) = certificate.parent.to_recurrence_parts().unwrap();
    wrong.parent =
        RecurrenceJson::from_recurrence_rational(&parent, 0, &rows(&[&[1], &[1], &[1, 1]]), None);
    assert!(wrong.verify(&rule, first, &seeds).is_err());
    // Another transform.
    let mut wrong = certificate.clone();
    wrong.transform = IndexSectionTransform::Identity;
    assert!(wrong.verify(&rule, first, &seeds).is_err());
    let mut wrong = certificate.clone();
    wrong.transform = IndexSectionTransform::Reversal {
        degree_slope: 1,
        degree_offset: 1,
    };
    assert!(wrong.verify(&rule, first, &seeds).is_err());
    // Too many exceptional rows.
    let mut wrong = certificate.clone();
    wrong.exceptional_rows = 2;
    assert!(wrong.verify(&rule, first, &seeds).is_err());
    // The certificate does not fit another family.
    let (other_rule, other_seeds, other_first) = catalog("A193636");
    assert!(certificate
        .verify(&other_rule, other_first, &other_seeds)
        .is_err());

    // The same edits inside the serialized report.
    assert!(tampered(&report, |v| v["index_section"]["period"] = json!(2)));
    assert!(tampered(&report, |v| {
        v["index_section"]["parent"]["initial_polynomials"][2] = json!(["1", "1"])
    }));
    assert!(tampered(&report, |v| v["index_section"]
        ["exceptional_rows"] = json!(1)));
    assert!(tampered(&report, |v| {
        v["index_section"]["transform"] = json!({"kind": "identity"})
    }));
    assert!(tampered(&report, |v| v["index_section"]
        ["divisors_examined"] = json!(0)));
    assert!(tampered(&report, |v| v["index_section"]
        ["parent_maximum_start_delay"] =
        json!(99)));
    assert!(tampered(&report, |v| {
        v["index_section"]["selected"] = serde_json::Value::Null
    }));
    assert!(tampered(&report, |v| v["index_section"] = serde_json::Value::Null));
    assert!(tampered(&report, |v| v["conclusion"]["first_index"] = json!(1)));
    assert!(tampered(&report, |v| v["theorem"] =
        json!("residue_class_subsequences")));
}

// ---------------------------------------------------------------------------
// Discovery and the parent provers
// ---------------------------------------------------------------------------

#[test]
fn discovery_records_the_lift_only_where_a_section_exists() {
    let (rule, seeds, first) = catalog("A109955");
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
    assert_eq!(attempt.theorem, RealRootednessTheorem::IndexSectionLift);
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

    // Chebyshev rows: no section parent, so no attempt is recorded.
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, b(&[&[0, 2]])),
            RecurrenceTerm::new(2, 0, b(&[&[-1]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let record = RecurrenceJson::from_recurrence_rational(&rule, 0, &rows(&[&[1], &[0, 1]]), None);
    let report = discover_recurrence_proofs(&record, &options()).unwrap();
    assert!(report
        .attempts
        .iter()
        .all(|a| a.theorem != RealRootednessTheorem::IndexSectionLift));
}

#[test]
fn the_parent_provers_use_the_lift_for_an_open_parent() {
    let (rule, seeds, first) = catalog("A109955");
    let record = ParentRecurrenceRecord {
        rule: &rule,
        seeds: &seeds,
        first_index: first,
    };
    let proof = prove_parent_real_rootedness(&record, first, &[], false).unwrap();
    assert_eq!(proof.status, VeroneseParentProofStatus::Proved);
    assert_eq!(
        proof.route,
        Some(VeroneseParentProofRoute::IndexSectionLift)
    );
    assert!(proof.interlacing.is_none());
    check_recorded_parent_proof(&record, first, &[], false, &proof).unwrap();
}
