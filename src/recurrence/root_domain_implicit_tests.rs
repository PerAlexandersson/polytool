//! Tests of root-domain information on implicit transports.

use super::*;
use serde_json::json;

fn run(value: serde_json::Value) -> ImplicitTransformReport {
    let (request, root_domain) = parse_implicit_transform_input(value).unwrap();
    let report =
        implicit_transform_recurrence_with_root_domain(&request, root_domain.as_ref()).unwrap();
    verify_recurrence_implicit(&report).unwrap();
    report
}

fn claims(report: &ImplicitRootDomainReport, evidence: RootEvidence) -> Vec<&RootDomainClaim> {
    report
        .claims
        .iter()
        .filter(|c| c.evidence == evidence)
        .collect()
}

/// Binomial rows divided by `k!`: the defining rule `P_m = (1 + x) P_(m-1)`
/// closes a Rolle induction, and `1/k!` is a certified multiplier sequence,
/// so every `Q_n = sum_k binom(n,k) x^k/k!` is real-rooted.
#[test]
fn factorial_multiplier_transfers_the_source_induction() {
    let report = run(json!({
        "oeis_id": "A007318", "row_count": 8, "root_domain": {},
        "transformation": {"kind": "diagonal", "weight": {"kind": "factorial_power", "power": 1}}
    }));
    let root = report.root_domain.as_ref().unwrap();
    let implication = &root.implications[0];
    assert_eq!(implication.theorem, "polya_schur_multiplier_sequence");
    assert_eq!(
        implication.direction,
        RootImplicationDirection::SourceToTarget
    );
    assert!(implication
        .hypotheses
        .iter()
        .all(|h| h.status != HypothesisStatus::Violated));
    let theorem = root.source_rule_theorem.as_ref().unwrap();
    assert!(theorem.matched, "{theorem:?}");
    let induction = claims(root, RootEvidence::CompletedRecurrenceInduction);
    assert_eq!(induction.len(), 1);
    assert_eq!(induction[0].subject, RootSubject::Source);
    assert_eq!(
        induction[0].interval.as_ref().unwrap().to_string(),
        "[-1, -1]"
    );
    let family = claims(root, RootEvidence::CertifiedFamilyInvariant);
    assert_eq!(family.len(), 1);
    assert_eq!(family[0].subject, RootSubject::Transformed);
    assert_eq!(family[0].property, RootProperty::RealRooted);
    // Fixed-cap evidence at the largest observed degree 7.
    let jensen = root.jensen_fixed_cap.as_ref().unwrap();
    assert_eq!(jensen.cap, 7);
    assert_eq!(jensen.status, JensenCapStatus::PreservesUpToCap);
    assert_eq!(jensen.evidence, FIXED_CAP_EVIDENCE);
    assert!(root
        .observations
        .iter()
        .filter(|o| o.subject == RootSubject::Transformed)
        .all(|o| o.real_rooted));
}

/// A065826 divided by `k + 1` gives the Eulerian rows: `1/(k + 1)` is not a
/// multiplier sequence (the Jensen test fails at the cap), but its inverse
/// `k + 1` is, so only `Q_n => P_n` transfers.
#[test]
fn inverse_binomial_transfers_only_from_the_new_rows() {
    let report = run(json!({
        "oeis_id": "A065826", "row_count": 7, "search": {},
        "root_domain": {"assumptions": [
            {"subject": "transformed", "property": "real_rooted", "note": "Eulerian rows"},
            {"subject": "source", "property": "real_rooted"}
        ]},
        "transformation": {"kind": "diagonal", "weight": {"kind": "inverse_binomial", "m": 1}}
    }));
    let root = report.root_domain.as_ref().unwrap();
    assert_eq!(
        root.implications[0].direction,
        RootImplicationDirection::TargetToSource
    );
    assert_eq!(
        root.jensen_fixed_cap.as_ref().unwrap().status,
        JensenCapStatus::FailsAtCap
    );
    let supplied = claims(root, RootEvidence::SuppliedAssumption);
    // Two supplied claims and one transfer (Q => P); P => Q is not drawn.
    assert_eq!(supplied.len(), 3);
    let transferred = supplied
        .iter()
        .find(|c| c.derivation.starts_with("transferred"))
        .unwrap();
    assert_eq!(transferred.subject, RootSubject::Source);
    assert!(transferred.conditional_on[0].starts_with("assumption 0"));
    assert!(root
        .notes
        .iter()
        .any(|n| n.starts_with("assumption 1") && n.contains("not transferred")));
    assert!(root
        .assumption_checks
        .iter()
        .all(|c| c.status == "consistent_with_rows"));
    // The certified Eulerian fit is not of Rolle form.
    let fitted = root.fitted_rule_theorem.as_ref().unwrap();
    assert!(!fitted.matched);
    assert!(claims(root, RootEvidence::CertifiedFamilyInvariant).is_empty());
}

/// `P_n = (x Q_n)'` with `Q_n = (1 + x)^n`: the successor rule of `Q_n` is
/// certified by reconstruction, Rolle closes on `[-1, -1]`, and the
/// antecedent adapter puts the roots of `P_n` in the hull with the root of
/// `g = x`.
#[test]
fn derivative_antecedent_transfers_a_certified_fit() {
    // P_n = (1 + x)^n + n x (1 + x)^(n-1): P_m = 2(1 + x) P_(m-1) - (1 + x)^2 P_(m-2).
    let q = |n: i64| Q::from_integer(BigInt::from(n));
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(
                1,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![q(2), q(2)]]),
            ),
            RecurrenceTerm::new(
                2,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![q(-1), q(-2), q(-1)]]),
            ),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let source = serde_json::to_value(RecurrenceJson::from_recurrence_rational(
        &rule,
        0,
        &[vec![q(1)], vec![q(1), q(2)]],
        None,
    ))
    .unwrap();
    let report = run(json!({
        "source_recurrence": source, "row_count": 8, "search": {}, "root_domain": {},
        "transformation": {"kind": "derivative_antecedent", "g": ["0", "1"]}
    }));
    assert_eq!(report.transformed_rows[3], vec!["1", "3", "3", "1"]);
    let search = report.independent_search.as_ref().unwrap();
    assert_eq!(
        search.evidence, CERTIFIED_BY_RECONSTRUCTION_EVIDENCE,
        "{search:?}"
    );
    let root = report.root_domain.as_ref().unwrap();
    assert_eq!(root.implications[0].theorem, "rolle_derivative_antecedent");
    assert_eq!(
        root.implications[0].direction,
        RootImplicationDirection::TargetToSource
    );
    let fitted = root.fitted_rule_theorem.as_ref().unwrap();
    assert!(fitted.matched, "{fitted:?}");
    let family = claims(root, RootEvidence::CertifiedFamilyInvariant);
    assert_eq!(family.len(), 1);
    assert_eq!(family[0].subject, RootSubject::Source);
    assert_eq!(family[0].property, RootProperty::RealRootedIn);
    assert_eq!(family[0].interval.as_ref().unwrap().to_string(), "[-1, 0]");
    // The source rule has two lags, so it is not matched.
    assert!(!root.source_rule_theorem.as_ref().unwrap().matched);
}

/// Negative controls: an index-dependent weight gets no adapter, and a
/// non-real-rooted antecedent factor violates the Rolle hypothesis.
#[test]
fn implicit_adapters_with_violated_hypotheses() {
    let report = run(json!({
        "oeis_id": "A007318", "row_count": 6,
        "root_domain": {"assumptions": [{"subject": "source", "property": "real_rooted"}]},
        "transformation": {"kind": "diagonal", "weight": {"kind": "product", "scale": "1",
            "factors": [{"coeffs": {"coeffs": [["1", "1"], ["1"]]}, "exponent": 1}]}}
    }));
    let root = report.root_domain.as_ref().unwrap();
    assert_eq!(root.implications[0].target_property, "none");
    assert!(root.jensen_fixed_cap.is_none());
    assert!(root.notes.iter().any(|n| n.contains("depends on n")));
    assert!(claims(root, RootEvidence::CertifiedFamilyInvariant).is_empty());
    assert_eq!(claims(root, RootEvidence::SuppliedAssumption).len(), 1);

    let report = run(json!({
        "oeis_id": "A007318", "row_count": 6, "root_domain": {},
        "transformation": {"kind": "derivative_antecedent", "g": ["1", "0", "1"]}
    }));
    let root = report.root_domain.as_ref().unwrap();
    let hypothesis = &root.implications[0].hypotheses[0];
    assert_eq!(hypothesis.status, HypothesisStatus::Violated);
    assert!(claims(root, RootEvidence::CertifiedFamilyInvariant).is_empty());
}

#[test]
fn implicit_root_domain_is_optional_and_replayed() {
    let base = json!({
        "oeis_id": "A007318", "row_count": 5,
        "transformation": {"kind": "diagonal", "weight": {"kind": "factorial_power", "power": 1}}
    });
    let (request, root_domain) = parse_implicit_transform_input(base).unwrap();
    assert!(root_domain.is_none());
    let plain = implicit_transform_recurrence(&request).unwrap();
    assert!(plain.root_domain.is_none());
    assert!(serde_json::to_value(&plain)
        .unwrap()
        .get("root_domain")
        .is_none());
    let mut report = implicit_transform_recurrence_with_root_domain(
        &request,
        Some(&RootDomainRequest::default()),
    )
    .unwrap();
    let text = serde_json::to_string(&report).unwrap();
    let back: ImplicitTransformReport = serde_json::from_str(&text).unwrap();
    assert_eq!(back, report);
    verify_recurrence_implicit(&back).unwrap();
    report.root_domain.as_mut().unwrap().claims.pop();
    assert!(verify_recurrence_implicit(&report).is_err());
    assert!(parse_implicit_transform_input(json!({"root_domain": {"unknown": 1}})).is_err());
}
