//! The [`LiftSummary`] of a preserver-image certification.

use super::certify::{
    PreserverImageCandidate, PreserverImageCertificate, PreserverImageIdentityStatus,
};
use super::interval::{interval_image, interval_route_label, region_contains};
use super::operator::OperatorForm;
use super::transfer::{denominator_hypothesis, operator_transfer};
use crate::recurrence::lift::lift_summary::{
    gauge_hypothesis, nonpositive_ray, parent_proof_summary, real_line,
};
use crate::recurrence::lift::veronese_algebra::classify_row;
use crate::recurrence::{
    LiftIdentityStatus, LiftOperation, LiftParent, LiftParentProof, LiftRootTransfer, LiftSummary,
    RealRootednessTheoremReport, RootInterval, VeroneseParentProof,
};

/// The parts of a certification its summary is computed from.
pub(super) struct PreserverImageSummaryParts<'a> {
    pub(super) candidate: &'a PreserverImageCandidate,
    pub(super) input_first_index: usize,
    pub(super) identity: &'a PreserverImageIdentityStatus,
    pub(super) parent_proof: &'a VeroneseParentProof,
    pub(super) whole_line_parent_proof: Option<&'a RealRootednessTheoremReport>,
    pub(super) interval_parent_proof: Option<&'a RealRootednessTheoremReport>,
}

/// Whether the parent row at `index` is real-rooted (so a failure of the
/// one-sign parent statement there is a sign failure only).
fn parent_row_is_real_rooted(candidate: &PreserverImageCandidate, index: usize) -> bool {
    let Ok((rule, first_index, rows)) = candidate.parent.to_recurrence_parts() else {
        return false;
    };
    if index < first_index || index - first_index > 256 {
        return false;
    }
    let count = rule.generation_initial_count(first_index, rows.len());
    rule.generate_rows_rational(
        &rows[..count.min(rows.len())],
        first_index,
        index - first_index + 1,
    )
    .ok()
    .and_then(|generated| generated.get(index - first_index).cloned())
    .is_some_and(|row| {
        let shape = classify_row(&row);
        shape.is_zero || shape.is_real_rooted
    })
}

/// The [`LiftSummary`] of a preserver-image lift: the identity maps one to
/// one; the hypotheses are the side conditions of the operator's theorem,
/// `lambda` finite and nonzero and (for operators with index-dependent
/// coefficients) the operator defined at every `n >= n0`, all decided
/// exactly; the parent statement is that of the shared parent provers
/// (roots in `(-inf, 0]`), or real-rootedness on the whole line by the
/// three-term theorem when the operator's theorem needs no more, or a
/// bounded interval `[L, U]` from an interval parent report: its image under
/// the operator (the `interval` module) when the ray statement is not
/// proved or the image is strictly smaller than the ray region, else the
/// line when the operator's theorem needs only real-rootedness.
pub(super) fn preserver_image_summary(parts: PreserverImageSummaryParts<'_>) -> LiftSummary {
    let candidate = parts.candidate;
    let n0 = parts.input_first_index;
    let s = candidate.index_shift;
    let transfer = operator_transfer(&candidate.operator, n0);
    let identity = match parts.identity {
        PreserverImageIdentityStatus::Certified { certificate } => LiftIdentityStatus::Certified {
            method: match certificate.as_ref() {
                PreserverImageCertificate::OperatorDivision(_) => "forward_division".into(),
                PreserverImageCertificate::HadamardConjugation(_) => "hadamard_conjugation".into(),
            },
        },
        PreserverImageIdentityStatus::SupportedByRows { rows, .. } => {
            LiftIdentityStatus::SupportedByRows { rows: *rows }
        }
        PreserverImageIdentityStatus::Refuted { index, reason } => LiftIdentityStatus::Refuted {
            index: Some(*index),
            reason: reason.clone(),
        },
    };
    let mut hypotheses = transfer.hypotheses.clone();
    hypotheses.push(gauge_hypothesis(&candidate.gauge, n0));
    if let Ok(OperatorForm::Weyl { denominator, .. }) = candidate.operator.form() {
        hypotheses.push(denominator_hypothesis(&denominator, n0));
    }
    let nonpositive_input = transfer.from_nonpositive.clone().unwrap_or_else(real_line);
    let base = parent_proof_summary(parts.parent_proof);
    let ray_hypothesis = "every parent row used is zero or real-rooted with all roots in \
                          (-inf, 0] (coefficients of one sign)"
        .to_string();
    // A proved interval parent covering every parent row used, its region
    // and the operator's image of it.
    let interval = parts
        .interval_parent_proof
        .filter(|report| {
            report.is_proved()
                && report
                    .conclusion
                    .as_ref()
                    .is_some_and(|c| c.first_index as i64 <= parts.input_first_index as i64 + s)
        })
        .and_then(|report| {
            let conclusion = report.conclusion.as_ref()?;
            let region = conclusion.root_region.clone()?;
            let image = region.at(0).ok().and_then(|concrete| {
                interval_image(
                    &candidate.operator,
                    &concrete,
                    transfer.from_nonpositive.as_ref(),
                )
            });
            Some((report, conclusion.statement.clone(), region, image))
        });
    let base_proved = matches!(base, LiftParentProof::Proved { .. });
    let interval_proof =
        |report: &RealRootednessTheoremReport, statement: &str| LiftParentProof::Proved {
            route: interval_route_label(report).into(),
            statement: statement.to_string(),
        };
    let interval_hypothesis = |region: &RootInterval| {
        format!("every parent row used is real-rooted with every root in {region}")
    };
    let preferred_interval = interval
        .as_ref()
        .and_then(|(report, statement, region, image)| {
            let image = image.clone()?;
            let preferable = !base_proved
                || (region_contains(&nonpositive_input, &image) && image != nonpositive_input);
            preferable.then(|| {
                (
                    interval_proof(report, statement),
                    region.clone(),
                    image,
                    interval_hypothesis(region),
                )
            })
        });
    let (parent_proof, parent_region, input_region, parent_hypothesis) =
        match (preferred_interval, base) {
            (Some(chosen), _) => chosen,
            (None, base @ (LiftParentProof::Proved { .. } | LiftParentProof::Attested { .. })) => {
                (base, nonpositive_ray(), nonpositive_input, ray_hypothesis)
            }
            (None, _)
                if transfer.from_real_line
                    && parts.whole_line_parent_proof.is_some_and(|report| {
                        report.conclusion.as_ref().is_some_and(|c| {
                            c.first_index as i64 <= parts.input_first_index as i64 + s
                        })
                    }) =>
            {
                let report = parts.whole_line_parent_proof.expect("checked");
                let statement = report
                    .conclusion
                    .as_ref()
                    .map(|c| c.statement.clone())
                    .unwrap_or_default();
                (
                    LiftParentProof::Proved {
                        route: "three_term_recurrence_whole_line".into(),
                        statement,
                    },
                    report
                        .conclusion
                        .as_ref()
                        .and_then(|c| c.root_region.clone())
                        .unwrap_or_else(real_line),
                    real_line(),
                    "every parent row used is real-rooted (anywhere on the line)".to_string(),
                )
            }
            (None, _) if transfer.from_real_line && interval.is_some() => {
                let (report, statement, region, _) = interval.as_ref().expect("checked");
                (
                    interval_proof(report, statement),
                    region.clone(),
                    real_line(),
                    interval_hypothesis(region),
                )
            }
            (
                None,
                LiftParentProof::Failed {
                    index: Some(index),
                    statement,
                },
            ) if parent_row_is_real_rooted(candidate, index) => (
                LiftParentProof::Open {
                    statement: format!(
                        "{statement} (the row is real-rooted with coefficients of both signs; \
                         real-rootedness of every parent row is not proved)"
                    ),
                },
                nonpositive_ray(),
                nonpositive_input,
                ray_hypothesis,
            ),
            (None, other) => (other, nonpositive_ray(), nonpositive_input, ray_hypothesis),
        };
    let root_transfer = LiftRootTransfer {
        theorem: transfer.theorem.clone(),
        parent_hypothesis,
        parent_root_region: Some(parent_region),
        input_conclusion: format!(
            "every Q_n, n >= {n0}, is zero or real-rooted with roots in {input_region}"
        ),
        input_root_region: input_region,
    };
    LiftSummary::new(
        LiftOperation::PreserverImage {
            operator: candidate.operator.clone(),
            description: candidate.operator.text(),
        },
        n0,
        LiftParent {
            recurrence: candidate.parent.clone(),
            relation: format!(
                "Q_n = lambda(n) T_n(P_(n{s:+})) with T = {}, for every n >= {n0}",
                candidate.operator.text()
            ),
            gauge: Some(candidate.gauge.clone()),
            index_shift: s,
        },
        identity,
        root_transfer,
        hypotheses,
        parent_proof,
    )
}
