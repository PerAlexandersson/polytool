//! The affine-coordinate stage of proof discovery: A118424 and A118884
//! (`P_n = 2 P_(n-1) + (x - 1) P_(n-L)`, the fixed-delay window in
//! `u = x - 1`), hand-made affine images of the Eulerian polynomials (a
//! translation with the root region carried along, a negative slope with
//! the mirrored interlacing orientation), a family no coordinate proves,
//! the option that disables the stage, replay and tampering.

use polytool::oeis;
use polytool::real_rootedness::RealInterval;
use polytool::recurrence::{
    discover_recurrence_proofs, discover_recurrence_proofs_with_stage_options,
    AffineCoordinateOptions, AffineInterlacingOrientation, AlternativeRepresentationOptions,
    CatalogLiftOptions, RealRootednessTheorem, RealRootednessTheoremOptions, RecurrenceJson,
    RecurrenceProofDiscoveryReport, RecurrenceProofStageOptions, RootInterval,
};
use polytool::BigRational as Q;

fn catalog_record(id: &str) -> RecurrenceJson {
    let entry = oeis::by_id(id).expect("catalog entry");
    let (rule, seeds) = entry.recurrence_parts().expect("decodable");
    RecurrenceJson::from_recurrence_rational(&rule, entry.recurrence_first_index, &seeds, None)
}

fn record(json: serde_json::Value) -> RecurrenceJson {
    serde_json::from_value(json).expect("a recurrence record")
}

/// Discovery with only the theorems, the tail and the affine stage (the
/// two expensive stages in between are switched off).
fn discover_quickly(record: &RecurrenceJson) -> RecurrenceProofDiscoveryReport {
    discover_recurrence_proofs_with_stage_options(
        record,
        &RealRootednessTheoremOptions::default(),
        &RecurrenceProofStageOptions {
            alternative_representations: AlternativeRepresentationOptions {
                max_candidates: 0,
                ..Default::default()
            },
            catalog_lifts: CatalogLiftOptions::disabled(),
            affine_coordinates: AffineCoordinateOptions::default(),
        },
    )
    .expect("discovery runs")
}

fn rational(text: &str) -> Q {
    polytool::recurrence::parse_rational_coeff(text).unwrap()
}

/// Every one of the first `count` rows of the record is nonzero,
/// real-rooted and has every root in `region` (checked exactly).
fn rows_in_region(record: &RecurrenceJson, region: &RootInterval, count: usize) {
    let (rule, first, seeds) = record.to_recurrence_parts().unwrap();
    let rows = rule.generate_rows_rational(&seeds, first, count).unwrap();
    for (offset, row) in rows.iter().enumerate() {
        let interval: RealInterval = region.at(first + offset).unwrap();
        let isolation = polytool::real_rootedness::isolate_real_roots_rational_coeffs(row);
        assert!(
            !isolation.is_zero_polynomial() && isolation.real_rooted_in(&interval),
            "row {} = {row:?} has a root outside {interval:?}",
            first + offset
        );
    }
}

#[test]
fn a118424_and_a118884_are_fixed_delay_windows_in_u_equal_x_minus_1() {
    for id in ["A118424", "A118884"] {
        let original = catalog_record(id);
        let report = discover_quickly(&original);
        assert_eq!(report.status, "proved", "{id}");
        assert_eq!(report.real_rooted_from_index, Some(0), "{id}");
        let search = report.affine_coordinates.as_ref().expect("the stage ran");
        let conclusion = search.conclusion.as_ref().expect("a coordinate proves");
        assert_eq!(conclusion.coordinate.slope, "1");
        assert_eq!(conclusion.coordinate.intercept, "1");
        assert_eq!(conclusion.theorem, RealRootednessTheorem::VariableDelay);
        assert_eq!(
            conclusion.root_region,
            Some(RootInterval {
                lower: None,
                upper: Some(polytool::recurrence::RootEndpoint::rational(
                    &rational("1"),
                    false
                )),
            })
        );
        assert_eq!(
            conclusion.interlacing_orientation,
            Some(AffineInterlacingOrientation::Preserved)
        );
        assert!(
            report.statement.contains("Affine coordinate"),
            "{}",
            report.statement
        );
        rows_in_region(&original, conclusion.root_region.as_ref().unwrap(), 14);
        // The transformed record is u = x - 1: P_n(u + 1) = 2 P_(n-1)(u + 1)
        // + u P_(n-L)(u + 1).
        let proved = search.coordinates.last().unwrap();
        let (rule, _, seeds) = proved.recurrence.to_recurrence_parts().unwrap();
        assert_eq!(
            rule.terms[1].coeff.coeffs,
            vec![vec![rational("0"), rational("1")]]
        );
        // Constant seeds are unchanged.
        assert_eq!(seeds, original.to_recurrence_parts().unwrap().2);
        report.verify().unwrap();
        let json = serde_json::to_string(&report).unwrap();
        let back: RecurrenceProofDiscoveryReport = serde_json::from_str(&json).unwrap();
        back.verify().unwrap();
    }
}

#[test]
fn a118424_is_proved_with_the_default_stages() {
    let report = discover_recurrence_proofs(
        &catalog_record("A118424"),
        &RealRootednessTheoremOptions::default(),
    )
    .unwrap();
    assert_eq!(report.status, "proved");
    assert!(report.catalog_lifts.is_some(), "the lifts ran first");
    assert!(report
        .affine_coordinates
        .as_ref()
        .is_some_and(|search| search.is_proved()));
    report.verify().unwrap();
}

/// The Eulerian polynomials `A_n(x) = (1 + (n-1) x) A_(n-1) + x(1-x)
/// A'_(n-1)`, `A_1 = 1`, written in the variable `y` with `x = y + 2`.
fn shifted_eulerian() -> RecurrenceJson {
    record(serde_json::json!({
        "schema": "polytool.recurrence.v1",
        "first_index": 1,
        "initial_polynomials": [["1"]],
        "recurrence": {"terms": [
            {"offset": 1, "deriv_order": 0, "sign": "none",
             "coeff": {"coeffs": [["-1", "-1"], ["2", "1"]]}},
            {"offset": 1, "deriv_order": 1, "sign": "none",
             "coeff": {"coeffs": [["-2", "-3", "-1"]]}}
        ]}
    }))
}

#[test]
fn shifted_eulerian_rows_are_proved_with_the_region_mapped() {
    let original = shifted_eulerian();
    let report = discover_quickly(&original);
    assert_eq!(report.status, "proved");
    assert!(report.attempts.iter().all(|a| !a.is_proved()));
    let search = report.affine_coordinates.as_ref().unwrap();
    let conclusion = search.conclusion.as_ref().unwrap();
    assert_eq!(conclusion.coordinate.slope, "1");
    // The roots of A_n are negative, so those of the shifted rows lie below
    // -2; the region of the inner proof is carried to the variable y.
    let region = conclusion.root_region.as_ref().expect("a region");
    assert!(region.lower.is_none());
    let upper = rational(&region.upper.as_ref().unwrap().value);
    assert!(upper <= rational("-1"), "{region:?}");
    rows_in_region(&original, region, 12);
    report.verify().unwrap();

    // Tampering: a region moved inwards, a coordinate replaced, a seed of the
    // transformed record, a dropped coordinate, a flipped proved flag.
    let mut tampered = report.clone();
    let stage = tampered.affine_coordinates.as_mut().unwrap();
    stage
        .conclusion
        .as_mut()
        .unwrap()
        .root_region
        .as_mut()
        .unwrap()
        .upper
        .as_mut()
        .unwrap()
        .value = "-3".into();
    assert!(tampered.verify().is_err());

    let mut tampered = report.clone();
    let stage = tampered.affine_coordinates.as_mut().unwrap();
    let last = stage.coordinates.last_mut().unwrap();
    last.coordinate.intercept = "-5".into();
    assert!(tampered.verify().is_err());

    let mut tampered = report.clone();
    let stage = tampered.affine_coordinates.as_mut().unwrap();
    let last = stage.coordinates.last_mut().unwrap();
    last.recurrence.initial_polynomials[0] = vec!["2".into()];
    assert!(tampered.verify().is_err());

    let mut tampered = report.clone();
    let stage = tampered.affine_coordinates.as_mut().unwrap();
    if stage.coordinates.len() > 1 {
        stage.coordinates.remove(0);
        assert!(tampered.verify().is_err());
    }

    let mut tampered = report.clone();
    let stage = tampered.affine_coordinates.as_mut().unwrap();
    stage.coordinates.last_mut().unwrap().proved = false;
    stage.conclusion = None;
    assert!(tampered.verify().is_err());

    // A dropped stage leaves a summary its attempts do not support.
    let mut tampered = report.clone();
    tampered.affine_coordinates = None;
    assert!(tampered.verify().is_err());
}

/// A118424 reflected, `P_n(x) = T_n(-x)`: `P_n = 2 P_(n-1) + (-x - 1)
/// P_(n-3)`.  In `u = x + 1` the delayed weight is `-u` (no route); in
/// `x = -u - 1` it is `u`, the fixed-delay window.
fn reflected_a118424() -> RecurrenceJson {
    record(serde_json::json!({
        "schema": "polytool.recurrence.v1",
        "first_index": 0,
        "initial_polynomials": [["1"], ["2"], ["4"]],
        "recurrence": {"terms": [
            {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["2"]]}},
            {"offset": 3, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["-1", "-1"]]}}
        ]}
    }))
}

#[test]
fn a_negative_slope_mirrors_the_interlacing_and_reverses_the_region() {
    let original = reflected_a118424();
    let report = discover_quickly(&original);
    assert_eq!(report.status, "proved");
    let search = report.affine_coordinates.as_ref().unwrap();
    assert!(!search.coordinates[0].proved, "u = x + 1 gives no route");
    let conclusion = search.conclusion.as_ref().unwrap();
    assert_eq!(conclusion.coordinate.slope, "-1");
    assert_eq!(conclusion.coordinate.intercept, "-1");
    assert_eq!(conclusion.theorem, RealRootednessTheorem::VariableDelay);
    // (-inf, 0) for u is (-1, +inf) for x = -u - 1.
    let region = conclusion.root_region.as_ref().expect("a region");
    assert!(region.upper.is_none(), "{region:?}");
    let lower = region.lower.as_ref().unwrap();
    assert_eq!((lower.value.as_str(), lower.closed), ("-1", false));
    rows_in_region(&original, region, 14);
    assert_eq!(
        conclusion.interlacing_orientation,
        Some(AffineInterlacingOrientation::Mirrored)
    );
    assert!(conclusion.statement.contains("mirrored"));
    report.verify().unwrap();
    // The orientation is recomputed, not trusted.
    let mut tampered = report.clone();
    let stage = tampered.affine_coordinates.as_mut().unwrap();
    stage.conclusion.as_mut().unwrap().interlacing_orientation =
        Some(AffineInterlacingOrientation::Preserved);
    assert!(tampered.verify().is_err());
}

#[test]
fn a_family_no_coordinate_proves_and_the_disabled_stage() {
    // P_n = (x - 1) P_(n-1) + 4 P_(n-2): in u = x - 1 the rows are
    // u Q_(n-1) + 4 Q_(n-2), whose roots are not real (Q_2 = u^2 + 4).
    let negative = record(serde_json::json!({
        "schema": "polytool.recurrence.v1",
        "first_index": 0,
        "initial_polynomials": [["1"], ["-1", "1"]],
        "recurrence": {"terms": [
            {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["-1", "1"]]}},
            {"offset": 2, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["4"]]}}
        ]}
    }));
    let report = discover_quickly(&negative);
    assert_ne!(report.status, "proved");
    let search = report.affine_coordinates.as_ref().expect("the stage ran");
    assert!(search.conclusion.is_none());
    assert_eq!(search.coordinates.len(), search.candidates_available.min(8));
    assert!(search.coordinates.iter().all(|c| !c.proved));
    report.verify().unwrap();
    // A coordinate claimed proved without a proof is rejected.
    let mut tampered = report.clone();
    tampered.affine_coordinates.as_mut().unwrap().coordinates[0].proved = true;
    assert!(tampered.verify().is_err());

    // Switched off, the report has no stage and A118424 stays unproved.
    let off = discover_recurrence_proofs_with_stage_options(
        &catalog_record("A118424"),
        &RealRootednessTheoremOptions::default(),
        &RecurrenceProofStageOptions {
            alternative_representations: AlternativeRepresentationOptions {
                max_candidates: 0,
                ..Default::default()
            },
            catalog_lifts: CatalogLiftOptions::disabled(),
            affine_coordinates: AffineCoordinateOptions::disabled(),
        },
    )
    .unwrap();
    assert_ne!(off.status, "proved");
    assert!(off.affine_coordinates.is_none());
    off.verify().unwrap();
}

#[test]
fn a_requested_region_is_pulled_back_to_the_new_coordinate() {
    let original = catalog_record("A118424");
    let options = RealRootednessTheoremOptions {
        root_region: Some(RealInterval {
            lower: None,
            lower_closed: false,
            upper: Some(rational("3/2")),
            upper_closed: false,
        }),
        ..Default::default()
    };
    let report = discover_recurrence_proofs_with_stage_options(
        &original,
        &options,
        &RecurrenceProofStageOptions {
            alternative_representations: AlternativeRepresentationOptions {
                max_candidates: 0,
                ..Default::default()
            },
            catalog_lifts: CatalogLiftOptions::disabled(),
            affine_coordinates: AffineCoordinateOptions::default(),
        },
    )
    .unwrap();
    assert_eq!(report.status, "proved");
    let search = report.affine_coordinates.as_ref().unwrap();
    let proved = search.coordinates.last().unwrap();
    // The request (-inf, 3/2) for x is (-inf, 1/2) for u = x - 1.
    let requested = proved.attempts[0].requested_root_region.as_ref().unwrap();
    assert_eq!(requested.upper.as_ref().unwrap().value, "1/2");
    report.verify().unwrap();
}
