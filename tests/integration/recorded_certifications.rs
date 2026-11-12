//! Replaying certifications recorded by an earlier, weaker version of the
//! parent provers: Veronese-section certifications whose parent was open or
//! attested then, and is proved outright by the iterated-preserver route now,
//! keep verifying with their weaker verdict, while a record claiming more
//! than the recomputation is refused.

use polytool::recurrence::{
    certify_veronese_section_lift, check_recorded_parent_proof, ParentRecurrenceRecord,
    RecurrenceJson, VeroneseParentProofRoute, VeroneseParentProofStatus,
    VeroneseSectionCertification, VeroneseSectionCertifyOptions, VeroneseSectionVerdict,
};
use serde::Deserialize;

#[derive(Deserialize)]
struct Recorded {
    input: RecurrenceJson,
    certification: VeroneseSectionCertification,
}

fn fixture(name: &str) -> Recorded {
    let path = format!(
        "{}/tests/fixtures/recorded-certifications/{name}",
        env!("CARGO_MANIFEST_DIR")
    );
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

/// The certification recomputed as it was recorded (16 input rows, no
/// root-domain route).
fn recompute(recorded: &Recorded) -> VeroneseSectionCertification {
    let (rule, first_index, seeds) = recorded.input.to_recurrence_parts().unwrap();
    let rows = rule
        .generate_rows_rational(&seeds, first_index, 16)
        .unwrap()
        .into_iter()
        .map(|mut row| {
            while row.last().is_some_and(num_traits::Zero::is_zero) {
                row.pop();
            }
            row
        })
        .collect::<Vec<_>>();
    certify_veronese_section_lift(
        &rule,
        &rows,
        first_index,
        &recorded.certification.candidate,
        &VeroneseSectionCertifyOptions {
            try_root_domain: false,
            ..VeroneseSectionCertifyOptions::default()
        },
    )
    .unwrap()
}

fn parent_record_check(recorded: &Recorded) -> Result<(), polytool::recurrence::RecurrenceError> {
    let parent = &recorded.certification.candidate.parent;
    let (rule, first_index, seeds) = parent.to_recurrence_parts().unwrap();
    check_recorded_parent_proof(
        &ParentRecurrenceRecord {
            rule: &rule,
            seeds: &seeds,
            first_index,
        },
        recorded.certification.input_first_index,
        &[],
        false,
        &recorded.certification.parent_proof,
    )
}

#[test]
fn weaker_recorded_parent_verdicts_keep_replaying() {
    for (name, status, verdict) in [
        (
            "veronese-parent-open.json",
            VeroneseParentProofStatus::NotProved,
            VeroneseSectionVerdict::IdentityCertifiedParentOpen,
        ),
        (
            "veronese-parent-attested.json",
            VeroneseParentProofStatus::ConditionalOnAttestation,
            VeroneseSectionVerdict::ConditionalOnAttestation,
        ),
    ] {
        let recorded = fixture(name);
        let certification = &recorded.certification;
        assert_eq!(certification.parent_proof.status, status, "{name}");
        assert_eq!(certification.verdict, verdict, "{name}");
        certification.verify().unwrap();
        let recomputed = recompute(&recorded);
        // The parent is now proved outright by the iterated-preserver route.
        assert_ne!(&recomputed, certification, "{name}");
        assert_eq!(recomputed.verdict, VeroneseSectionVerdict::Proved, "{name}");
        assert_eq!(
            recomputed.parent_proof.route,
            Some(VeroneseParentProofRoute::IteratedPreserver),
            "{name}"
        );
        certification
            .check_against_recomputation(&recomputed)
            .unwrap_or_else(|error| panic!("{name}: {error}"));
        recomputed.check_against_recomputation(&recomputed).unwrap();
        assert!(certification
            .parent_proof
            .is_supported_by(&recomputed.parent_proof));
        parent_record_check(&recorded).unwrap();

        // Claiming more than the recomputation: the old record with a proved
        // verdict but its old parent proof is inconsistent, and a parent
        // proof claiming another route is not supported.
        let mut tampered = certification.clone();
        tampered.verdict = VeroneseSectionVerdict::Proved;
        assert!(tampered.check_against_recomputation(&recomputed).is_err());
        let mut tampered = recomputed.clone();
        tampered.parent_proof.route = Some(VeroneseParentProofRoute::LiuWang);
        assert!(!tampered
            .parent_proof
            .is_supported_by(&recomputed.parent_proof));
        assert!(tampered.check_against_recomputation(&recomputed).is_err());
        // A changed statement of the weaker record is refused.
        let mut tampered = certification.clone();
        tampered.statement.push_str(" (edited)");
        assert!(tampered.check_against_recomputation(&recomputed).is_err());
        // An open record claiming an interlacing is refused.
        let mut tampered = recorded.certification.parent_proof.clone();
        tampered.interlacing = recomputed.parent_proof.interlacing.clone().or_else(|| {
            Some(polytool::recurrence::LiftInterlacing {
                kind: polytool::recurrence::RowInterlacing::Weak,
                orientation: polytool::recurrence::InterlacingOrientation::PreviousInterlacesNext,
                from_index: 1,
                justification: polytool::recurrence::LiftInterlacingJustification {
                    theorem: "none".into(),
                    premises: Vec::new(),
                },
            })
        });
        assert!(!tampered.is_supported_by(&recomputed.parent_proof));
    }
}

#[test]
fn a_recorded_proof_stronger_than_the_recomputation_is_refused() {
    // The recomputation with the old (weaker) record as the "current"
    // result: a proved record is not supported by an open or attested one
    // as a comparison; it is accepted only because its own route, replayed
    // alone, proves the parent, and refused when that route does not.
    let recorded = fixture("veronese-parent-attested.json");
    let recomputed = recompute(&recorded);
    assert!(!recomputed
        .parent_proof
        .is_supported_by(&recorded.certification.parent_proof));
    recomputed
        .check_against_recomputation(&recorded.certification)
        .unwrap();
    let mut wrong_route = recomputed.clone();
    wrong_route.parent_proof.route = Some(VeroneseParentProofRoute::PeriodicProduct);
    assert!(wrong_route
        .check_against_recomputation(&recorded.certification)
        .is_err());
    let open = fixture("veronese-parent-open.json");
    let recomputed = recompute(&open);
    assert!(!recomputed
        .parent_proof
        .is_supported_by(&open.certification.parent_proof));
}

/// A basis-map lift recorded when its coefficient family was only attested
/// (catalog flag), with the parent root region of the attestation: the
/// parent is now proved, in a smaller region, and the record still replays;
/// a record claiming the proved verdict with the old parent does not.
#[test]
fn a_basis_map_lift_with_a_weaker_recorded_parent_keeps_replaying() {
    let path = format!(
        "{}/tests/fixtures/recorded-certifications/basis-map-parent-attested.json",
        env!("CARGO_MANIFEST_DIR")
    );
    let recorded: polytool::recurrence::BasisMapLiftCertification =
        serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    assert_eq!(
        recorded.lift_summary.verdict,
        polytool::recurrence::LiftVerdict::ConditionalOnAttestation
    );
    recorded.verify().unwrap();
    let mut tampered = recorded.clone();
    tampered.lift_summary.verdict = polytool::recurrence::LiftVerdict::Proved;
    assert!(tampered.verify().is_err());
    let mut tampered = recorded.clone();
    tampered.lift_summary.root_transfer.input_root_region.upper = None;
    assert!(tampered.verify().is_err());
}

fn basis_map_fixture(name: &str) -> polytool::recurrence::BasisMapLiftCertification {
    let path = format!(
        "{}/tests/fixtures/recorded-certifications/{name}",
        env!("CARGO_MANIFEST_DIR")
    );
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

/// Basis-map lifts recorded with their coefficient family proved by the
/// root-domain and the paper-theorem routes.  The current provers settle
/// those families first by an earlier route (product induction), so the
/// default recomputation differs; the record replays through its own route
/// alone.  A record naming a route whose replay fails is refused.
#[test]
fn proofs_by_routes_the_cascade_no_longer_reaches_replay_by_their_route() {
    use polytool::recurrence::{certify_basis_map_lift, BasisMapLiftOptions, LiftParentProof};
    for (name, route) in [
        ("basis-map-parent-root-domain.json", "root_domain"),
        ("basis-map-parent-paper-theorem.json", "paper_theorem"),
    ] {
        let recorded = basis_map_fixture(name);
        let LiftParentProof::Proved {
            route: recorded_route,
            ..
        } = &recorded.lift_summary.parent_proof
        else {
            panic!("{name}: the recorded parent is not proved");
        };
        assert_eq!(
            recorded_route,
            &format!("basis_map_coefficient_roots({route})"),
            "{name}"
        );
        // The default recomputation now proves the family by another route.
        let (rule, n0, seeds) = recorded.input.to_recurrence_parts().unwrap();
        let rows = rule
            .generate_rows_rational(&seeds, n0, seeds.len() + 2)
            .unwrap();
        let recomputed = certify_basis_map_lift(
            &rule,
            &rows,
            n0,
            &recorded.candidate,
            &BasisMapLiftOptions::default(),
        )
        .unwrap();
        assert_ne!(
            recomputed.lift_summary.parent_proof, recorded.lift_summary.parent_proof,
            "{name}"
        );
        recorded
            .verify()
            .unwrap_or_else(|error| panic!("{name}: {error}"));
        // A recorded route that does not prove the family is refused.
        // (A route that also proves the family would be accepted: the claim
        // is then true; an attestation is never replayed.)
        for other in ["attestation", "no_such_route"] {
            let mut tampered = recorded.clone();
            if let LiftParentProof::Proved { route, .. } = &mut tampered.lift_summary.parent_proof {
                *route = format!("basis_map_coefficient_roots({other})");
            }
            assert!(tampered.verify().is_err(), "{name} as {other}");
        }
        // A changed statement is refused.
        let mut tampered = recorded.clone();
        if let LiftParentProof::Proved { statement, .. } = &mut tampered.lift_summary.parent_proof {
            statement.push_str(" (edited)");
        }
        assert!(tampered.verify().is_err(), "{name}");
    }
}

/// A parent-prover record proved by the paper-theorem route: the default
/// cascade now proves it by product induction first, and the record replays
/// through the paper-theorem route alone; claiming another route, or an
/// interlacing the route does not give, is refused.
#[test]
fn a_parent_proof_by_the_paper_theorem_replays_by_its_route() {
    #[derive(Deserialize)]
    struct Recorded {
        record: RecurrenceJson,
        try_root_domain: bool,
        proof: polytool::recurrence::VeroneseParentProof,
    }
    let path = format!(
        "{}/tests/fixtures/recorded-certifications/parent-proof-paper-theorem.json",
        env!("CARGO_MANIFEST_DIR")
    );
    let recorded: Recorded = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    assert_eq!(
        recorded.proof.route,
        Some(VeroneseParentProofRoute::PaperTheorem)
    );
    let (rule, first_index, seeds) = recorded.record.to_recurrence_parts().unwrap();
    let record = ParentRecurrenceRecord {
        rule: &rule,
        seeds: &seeds,
        first_index,
    };
    let current = polytool::recurrence::prove_parent_real_rootedness(
        &record,
        first_index,
        &[],
        recorded.try_root_domain,
    )
    .unwrap();
    assert_ne!(current.route, recorded.proof.route);
    assert!(!recorded.proof.is_supported_by(&current));
    check_recorded_parent_proof(
        &record,
        first_index,
        &[],
        recorded.try_root_domain,
        &recorded.proof,
    )
    .unwrap();
    let check = |proof: &polytool::recurrence::VeroneseParentProof| {
        check_recorded_parent_proof(&record, first_index, &[], recorded.try_root_domain, proof)
    };
    let mut tampered = recorded.proof.clone();
    tampered.route = Some(VeroneseParentProofRoute::RootDomain);
    assert!(check(&tampered).is_err());
    let mut tampered = recorded.proof.clone();
    tampered.route = Some(VeroneseParentProofRoute::ThreeTermRecurrence);
    assert!(check(&tampered).is_err());
    let mut tampered = recorded.proof.clone();
    tampered.statement.push_str(" (edited)");
    assert!(check(&tampered).is_err());
    let mut tampered = recorded.proof.clone();
    tampered.interlacing = match tampered.interlacing {
        Some(_) => None,
        None => current
            .interlacing
            .clone()
            .or(Some(polytool::recurrence::LiftInterlacing {
                kind: polytool::recurrence::RowInterlacing::Strict,
                orientation: polytool::recurrence::InterlacingOrientation::PreviousInterlacesNext,
                from_index: first_index + 1,
                justification: polytool::recurrence::LiftInterlacingJustification {
                    theorem: "none".into(),
                    premises: Vec::new(),
                },
            })),
    };
    assert!(check(&tampered).is_err());
}
