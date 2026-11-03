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
    // result: a proved record is not supported by an open or attested one.
    let recorded = fixture("veronese-parent-attested.json");
    let recomputed = recompute(&recorded);
    assert!(!recomputed
        .parent_proof
        .is_supported_by(&recorded.certification.parent_proof));
    assert!(recomputed
        .check_against_recomputation(&recorded.certification)
        .is_err());
    let open = fixture("veronese-parent-open.json");
    let recomputed = recompute(&open);
    assert!(!recomputed
        .parent_proof
        .is_supported_by(&open.certification.parent_proof));
}
