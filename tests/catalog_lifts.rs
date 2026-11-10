//! The catalog-lift stage of proof discovery: a normalized identification
//! with a catalog family (A123125, the Eulerian rows A008292 times `x`), a
//! backward basis-map lift (A067613 through the Hermite map of A060821),
//! depth one (the parent report has no lift stage), the option that
//! disables the stage, replay, and tampering (a changed conclusion, a parent
//! relabelled as attested or attested as proved, a dropped stage, a parent
//! with a lift stage of its own).

use polytool::oeis;
use polytool::recurrence::{
    discover_recurrence_proofs, discover_recurrence_proofs_with_stage_options,
    discover_recurrence_proofs_with_stages, AffineCoordinateOptions,
    AlternativeRepresentationOptions, CatalogLiftOptions, CatalogLiftRoute, CatalogParentProofKind,
    RealRootednessTheoremOptions, RecurrenceJson, RecurrenceProofDiscoveryReport,
    RecurrenceProofStageOptions,
};

/// The displayed record of a catalog entry (prefix rows and initial rows
/// as seeds), as `polytool oeis info --json` exports it.
fn catalog_record(id: &str) -> RecurrenceJson {
    let entry = oeis::by_id(id).expect("catalog entry");
    let (rule, seeds, first) = entry.recurrence_export_parts().expect("decodable");
    RecurrenceJson::from_recurrence_rational(&rule, usize::try_from(first).unwrap(), &seeds, None)
}

fn discover(id: &str) -> RecurrenceProofDiscoveryReport {
    discover_recurrence_proofs(
        &catalog_record(id),
        &RealRootednessTheoremOptions::default(),
    )
    .expect("discovery runs")
}

#[test]
fn a123125_is_the_eulerian_family_times_x() {
    let report = discover("A123125");
    assert_eq!(report.status, "proved");
    assert_eq!(report.real_rooted_from_index, Some(0));
    let search = report.catalog_lifts.as_ref().expect("the stage ran");
    let conclusion = search.conclusion.as_ref().expect("a lift concludes");
    assert_eq!(conclusion.route, CatalogLiftRoute::NormalizedIdentification);
    assert_eq!(conclusion.catalog_family.as_deref(), Some("A008292"));
    assert_eq!(conclusion.parent_proof, CatalogParentProofKind::Polytool);
    assert!(conclusion.conditional_on.is_empty());
    assert!(report.statement.contains("A008292"), "{}", report.statement);
    // The parent is proved by the general stages, never by another lift.
    let identification = search.identifications.last().unwrap();
    let parent = identification.parent.as_ref().expect("parent report");
    assert_eq!(parent.status, "proved");
    assert!(parent.catalog_lifts.is_none());
    report.verify().unwrap();
    // The serialized report replays.
    let json = serde_json::to_string(&report).unwrap();
    let back: RecurrenceProofDiscoveryReport = serde_json::from_str(&json).unwrap();
    back.verify().unwrap();

    // A changed conclusion (with a matching summary) is rejected.
    let mut tampered = report.clone();
    let lifts = tampered.catalog_lifts.as_mut().unwrap();
    lifts.conclusion.as_mut().unwrap().from_index = 2;
    tampered.real_rooted_from_index = Some(2);
    assert!(tampered.verify().is_err());
    // An identification relabelled as attested (A008292 is flagged) is not
    // its catalog flag premise.
    let mut tampered = report.clone();
    let identification = tampered
        .catalog_lifts
        .as_mut()
        .unwrap()
        .identifications
        .last_mut()
        .unwrap();
    identification.parent_proof = Some(CatalogParentProofKind::CatalogAttestation);
    identification.attestation = Some("oeis_catalog:A008292:lean_proved".into());
    assert!(tampered.verify().is_err());
    // Dropping the stage leaves a summary its attempts do not support.
    let mut tampered = report.clone();
    tampered.catalog_lifts = None;
    assert!(tampered.verify().is_err());
    // A parent report with a lift stage of its own (depth two) is rejected.
    let mut tampered = report.clone();
    let identification = tampered
        .catalog_lifts
        .as_mut()
        .unwrap()
        .identifications
        .last_mut()
        .unwrap();
    let parent = identification.parent.as_mut().unwrap();
    parent.catalog_lifts = report.catalog_lifts.clone();
    assert!(tampered.verify().is_err());
    // The identification is about another record.
    let mut tampered = report.clone();
    tampered.catalog_lifts.as_mut().unwrap().identifications[0].catalog_id = "A173018".into();
    assert!(tampered.verify().is_err());
}

#[test]
fn a067613_is_a_hermite_basis_map_image() {
    let report = discover("A067613");
    assert_eq!(report.status, "proved");
    let search = report.catalog_lifts.as_ref().expect("the stage ran");
    let conclusion = search.conclusion.as_ref().expect("a lift concludes");
    assert_eq!(conclusion.route, CatalogLiftRoute::BasisMap);
    assert_eq!(conclusion.basis_map.as_deref(), Some("A060821"));
    assert_eq!(conclusion.parent_proof, CatalogParentProofKind::Polytool);
    // Nonzero rows come from the record's leading-coefficient certificate.
    assert!(conclusion
        .leading_coefficient
        .as_ref()
        .is_some_and(|c| c.is_proved()));
    report.verify().unwrap();
    // A conclusion claimed as conditional is not what the attempts give.
    let mut tampered = report.clone();
    let conclusion = tampered
        .catalog_lifts
        .as_mut()
        .unwrap()
        .conclusion
        .as_mut()
        .unwrap();
    conclusion.parent_proof = CatalogParentProofKind::CatalogAttestation;
    assert!(tampered.verify().is_err());
}

#[test]
fn the_stage_can_be_disabled() {
    // The affine-coordinate stage is switched off too: it proves A123125 on
    // its own (iterated preserver in u = x - 1).
    let report = discover_recurrence_proofs_with_stage_options(
        &catalog_record("A123125"),
        &RealRootednessTheoremOptions::default(),
        &RecurrenceProofStageOptions {
            alternative_representations: AlternativeRepresentationOptions::default(),
            catalog_lifts: CatalogLiftOptions::disabled(),
            affine_coordinates: AffineCoordinateOptions::disabled(),
        },
    )
    .unwrap();
    assert_eq!(report.status, "not_proved");
    assert!(report.catalog_lifts.is_none());
    report.verify().unwrap();
    // With only the lifts off, the affine-coordinate stage proves it.
    let affine = discover_recurrence_proofs_with_stages(
        &catalog_record("A123125"),
        &RealRootednessTheoremOptions::default(),
        &AlternativeRepresentationOptions::default(),
        &CatalogLiftOptions::disabled(),
    )
    .unwrap();
    assert_eq!(affine.status, "proved");
    assert!(affine
        .affine_coordinates
        .as_ref()
        .is_some_and(|search| search.is_proved()));
    affine.verify().unwrap();
}

#[test]
fn a_proved_record_runs_no_lift_stage() {
    let report = discover("A008292");
    assert_eq!(report.status, "proved");
    assert!(report.catalog_lifts.is_none());
    // A lift stage recorded although the record is proved is rejected.
    let mut tampered = report.clone();
    tampered.catalog_lifts = discover("A123125").catalog_lifts;
    assert!(tampered.verify().is_err());
}

#[test]
fn a_lift_resting_on_a_catalog_flag_is_only_conditional() {
    // A008303 is normalized-equal to A263789 (peak polynomials), which
    // polytool does not prove but the catalog flags lean_proved.
    let report = discover("A008303");
    assert_eq!(report.status, "conditional_on_attestation");
    assert_eq!(report.real_rooted_from_index, None);
    let conclusion = report
        .catalog_lifts
        .as_ref()
        .and_then(|search| search.conclusion.as_ref())
        .expect("a conditional conclusion");
    assert_eq!(
        conclusion.parent_proof,
        CatalogParentProofKind::CatalogAttestation
    );
    assert_eq!(conclusion.catalog_family.as_deref(), Some("A263789"));
    assert_eq!(
        conclusion.conditional_on,
        ["oeis_catalog:A263789:lean_proved".to_string()]
    );
    assert!(report.statement.contains("Not proved by polytool"));
    report.verify().unwrap();
    // Claiming a polytool proof is rejected, in the summary or in the lift.
    let mut tampered = report.clone();
    tampered.status = "proved".into();
    tampered.real_rooted_from_index = Some(conclusion.from_index);
    assert!(tampered.verify().is_err());
    let mut tampered = report.clone();
    let lifts = tampered.catalog_lifts.as_mut().unwrap();
    for identification in &mut lifts.identifications {
        if identification.attestation.is_some() {
            identification.parent_proof = Some(CatalogParentProofKind::Polytool);
            identification.attestation = None;
        }
    }
    assert!(tampered.verify().is_err());
}
