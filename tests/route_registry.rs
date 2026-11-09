//! The registry of real-rooted theorem routes: names, the frozen route
//! order, the parent routes that run a theorem, and the JSON form of the
//! route certificates of a theorem report.

use polytool::recurrence::{
    check_real_rootedness_theorem, check_residue_class_subsequence_theorem,
    check_variable_delay_theorem, BigRational, RealRootednessTheorem, RealRootednessTheoremOptions,
    RealRootednessTheoremReport, Recurrence, RouteCertificate, VeroneseParentProofRoute,
    DISCOVERY_THEOREMS, ROW_ROUTE_THEOREMS,
};
use serde_json::Value;

/// The serialized name of a value of a unit-variant enum.
fn serialized<T: serde::Serialize>(value: T) -> String {
    serde_json::to_value(value)
        .unwrap()
        .as_str()
        .unwrap()
        .to_string()
}

/// Every parent route (the enum has no `ALL`).
const PARENT_ROUTES: [VeroneseParentProofRoute; 16] = [
    VeroneseParentProofRoute::PeriodicProduct,
    VeroneseParentProofRoute::ProductRecurrence,
    VeroneseParentProofRoute::PaperTheorem,
    VeroneseParentProofRoute::RootDomain,
    VeroneseParentProofRoute::ThreeTermRecurrence,
    VeroneseParentProofRoute::LiuWang,
    VeroneseParentProofRoute::LiuWangDerivativeLag,
    VeroneseParentProofRoute::MaWangInterval,
    VeroneseParentProofRoute::LiuWangThreeTerm,
    VeroneseParentProofRoute::IteratedPreserver,
    VeroneseParentProofRoute::HypergeometricRow,
    VeroneseParentProofRoute::VariableDelay,
    VeroneseParentProofRoute::ResidueClassSubsequences,
    VeroneseParentProofRoute::IndexSectionLift,
    VeroneseParentProofRoute::FactorGauge,
    VeroneseParentProofRoute::Attestation,
];

/// Fails to compile when a parent route is added without extending
/// [`PARENT_ROUTES`].
#[allow(dead_code)]
fn parent_routes_are_listed(route: VeroneseParentProofRoute) {
    use VeroneseParentProofRoute as R;
    match route {
        R::PeriodicProduct
        | R::ProductRecurrence
        | R::PaperTheorem
        | R::RootDomain
        | R::ThreeTermRecurrence
        | R::LiuWang
        | R::LiuWangDerivativeLag
        | R::MaWangInterval
        | R::LiuWangThreeTerm
        | R::IteratedPreserver
        | R::HypergeometricRow
        | R::VariableDelay
        | R::ResidueClassSubsequences
        | R::IndexSectionLift
        | R::FactorGauge
        | R::Attestation => {}
    }
}

#[test]
fn theorem_names_are_their_serialized_names() {
    for theorem in RealRootednessTheorem::ALL {
        assert_eq!(theorem.as_str(), serialized(theorem));
        let parsed: RealRootednessTheorem =
            serde_json::from_value(serde_json::Value::from(theorem.as_str())).unwrap();
        assert_eq!(parsed, theorem);
    }
    for route in PARENT_ROUTES {
        assert_eq!(route.as_str(), serialized(route));
    }
}

/// The route order is policy: it is the order of the parent cascade, of the
/// discovery attempts and of the residue-class and index-section menus.
#[test]
fn the_theorem_order_is_frozen() {
    let names: Vec<&str> = RealRootednessTheorem::ALL
        .iter()
        .map(|theorem| theorem.as_str())
        .collect();
    assert_eq!(
        names,
        [
            "product_recurrence",
            "three_term_recurrence",
            "liu_wang",
            "liu_wang_derivative_lag",
            "ma_wang_interval",
            "liu_wang_three_term",
            "iterated_preserver",
            "hypergeometric_row",
            "variable_delay",
            "residue_class_subsequences",
            "index_section_lift",
            "factor_gauge",
        ]
    );
}

#[test]
fn every_theorem_has_exactly_one_parent_route_of_the_same_name() {
    for theorem in RealRootednessTheorem::ALL {
        let routes: Vec<VeroneseParentProofRoute> = PARENT_ROUTES
            .into_iter()
            .filter(|route| route.theorem() == Some(theorem))
            .collect();
        assert_eq!(routes.len(), 1, "{theorem:?}");
        assert_eq!(routes[0].as_str(), theorem.as_str());
    }
    let without_theorem: Vec<&str> = PARENT_ROUTES
        .into_iter()
        .filter(|route| route.theorem().is_none())
        .map(VeroneseParentProofRoute::as_str)
        .collect();
    assert_eq!(
        without_theorem,
        [
            "periodic_product",
            "paper_theorem",
            "root_domain",
            "attestation"
        ]
    );
}

#[test]
fn the_discovery_lists_are_prefixes_of_the_theorem_order() {
    let listed: Vec<RealRootednessTheorem> = DISCOVERY_THEOREMS
        .iter()
        .chain(ROW_ROUTE_THEOREMS.iter())
        .copied()
        .collect();
    assert_eq!(listed, RealRootednessTheorem::ALL[..listed.len()]);
}

// ---------------------------------------------------------------------------
// Route certificates in JSON
// ---------------------------------------------------------------------------

/// The catalog family: rule, identifying seeds and first index.
fn catalog(id: &str) -> (Recurrence, Vec<Vec<BigRational>>, usize) {
    let entry = polytool::oeis::by_id(id).unwrap();
    let (rule, rows) = entry.recurrence_parts().unwrap();
    let count = rule.generation_initial_count(entry.recurrence_first_index, rows.len());
    (rule, rows[..count].to_vec(), entry.recurrence_first_index)
}

/// A proved fixed-delay report (certificate key `variable_delay`).
fn fixed_delay_report() -> RealRootednessTheoremReport {
    let (rule, seeds, first) = catalog("A201013");
    let report = check_variable_delay_theorem(&rule, first, &seeds, &Default::default()).unwrap();
    assert!(report.is_proved());
    report
}

/// A proved residue-class report (certificate key `residue_classes`).
fn residue_class_report() -> RealRootednessTheoremReport {
    let (rule, seeds, first) = catalog("A099605");
    let report =
        check_residue_class_subsequence_theorem(&rule, first, &seeds, &Default::default()).unwrap();
    assert!(report.is_proved());
    report
}

fn parse(value: Value) -> Result<RealRootednessTheoremReport, serde_json::Error> {
    serde_json::from_value(value)
}

#[test]
fn a_certificate_keeps_its_key_and_round_trips_byte_for_byte() {
    for (report, key, theorem) in [
        (
            fixed_delay_report(),
            "variable_delay",
            RealRootednessTheorem::VariableDelay,
        ),
        (
            residue_class_report(),
            "residue_classes",
            RealRootednessTheorem::ResidueClassSubsequences,
        ),
    ] {
        let certificate = report.certificate.as_ref().unwrap();
        assert_eq!(certificate.theorem(), theorem);
        let text = serde_json::to_string(&report).unwrap();
        // The certificate is the last field, after the Lean replay status.
        let lean_replay = text.find("\"lean_replay\":").unwrap();
        let certificate_key = text.find(&format!("\"{key}\":")).unwrap();
        assert!(lean_replay < certificate_key, "{key}");
        let parsed: RealRootednessTheoremReport = serde_json::from_str(&text).unwrap();
        assert_eq!(parsed, report);
        assert_eq!(serde_json::to_string(&parsed).unwrap(), text);
        parsed.verify().unwrap();
    }
    // The accessors see only the certificate of their own route.
    let report = fixed_delay_report();
    assert!(report.variable_delay().is_some());
    assert!(report.residue_classes().is_none());
    assert!(report.iterated_preserver().is_none());
    assert!(report.index_section().is_none());
    assert!(report.factor_gauge().is_none());
}

#[test]
fn a_report_without_a_certificate_has_no_certificate_key() {
    let (rule, seeds, first) = catalog("A201013");
    let report = check_real_rootedness_theorem(
        RealRootednessTheorem::LiuWang,
        &rule,
        first,
        &seeds,
        &RealRootednessTheoremOptions::default(),
    )
    .unwrap();
    assert!(report.certificate.is_none());
    let value = serde_json::to_value(&report).unwrap();
    for key in [
        "iterated_preserver",
        "variable_delay",
        "residue_classes",
        "index_section",
        "factor_gauge",
    ] {
        assert!(value.get(key).is_none(), "{key}");
    }
    assert_eq!(parse(value).unwrap(), report);
}

#[test]
fn two_certificate_keys_are_a_parse_error() {
    let mut value = serde_json::to_value(fixed_delay_report()).unwrap();
    let other = serde_json::to_value(residue_class_report()).unwrap();
    value["residue_classes"] = other["residue_classes"].clone();
    let error = parse(value).unwrap_err().to_string();
    assert!(error.contains("more than one route certificate"), "{error}");
}

#[test]
fn an_unknown_key_is_a_parse_error() {
    let mut value = serde_json::to_value(fixed_delay_report()).unwrap();
    value["certificate"] = value["variable_delay"].clone();
    let error = parse(value).unwrap_err().to_string();
    assert!(error.contains("unknown field `certificate`"), "{error}");
}

#[test]
fn a_malformed_certificate_is_a_parse_error() {
    let mut value = serde_json::to_value(fixed_delay_report()).unwrap();
    value["variable_delay"]["delay"] = Value::from("three");
    assert!(parse(value).is_err());
    let mut value = serde_json::to_value(fixed_delay_report()).unwrap();
    value["variable_delay"]["unexpected"] = Value::from(1);
    assert!(parse(value).is_err());
    // A certificate under the key of another route does not parse either.
    let mut value = serde_json::to_value(fixed_delay_report()).unwrap();
    let certificate = value["variable_delay"].take();
    value.as_object_mut().unwrap().remove("variable_delay");
    value["factor_gauge"] = certificate;
    assert!(parse(value).is_err());
}

#[test]
fn a_certificate_of_another_theorem_does_not_replay() {
    let mut report = fixed_delay_report();
    let other = residue_class_report();
    report.certificate = other.certificate.clone();
    assert!(matches!(
        report.certificate,
        Some(RouteCertificate::ResidueClasses(_))
    ));
    assert!(report.verify().is_err());
    report.certificate = None;
    assert!(report.verify().is_err());
}
