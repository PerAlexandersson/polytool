//! Exact all-index evidence shared by consecutive-interlacing transfers.
//!
//! A sampled degree pattern is only a proposal: the existing leading-coefficient
//! induction must certify it. Enclosing lift verifiers recompute these records.

use super::lift_summary::{fixed_preserver_interlacing, proved_claim};
use crate::recurrence::*;
use crate::{vec_poly, BigRational};

fn degree_certificate(record: &RecurrenceJson) -> Option<LeadingCoefficientCertificate> {
    let (rule, first, seeds) = record.to_recurrence_parts().ok()?;
    let certificate = certify_leading_coefficient(
        &rule,
        first,
        &seeds,
        &LeadingCoefficientOptions {
            sample_rows: 12,
            ..Default::default()
        },
    )
    .ok()?;
    certificate.is_proved().then_some(certificate)
}

/// Nonzero rows and successor degrees from the row before `from` onwards.
pub(super) fn successor_degrees(
    record: &RecurrenceJson,
    from: usize,
) -> Option<LeadingCoefficientCertificate> {
    let certificate = degree_certificate(record)?;
    let law = certificate.degree_law.as_ref()?;
    let first = from.checked_sub(1)?;
    let last = first.max(law.start).checked_add(law.period)?;
    if last.checked_sub(first)? > 512 || !certificate.degree_is_exact {
        return None;
    }
    // Past this finite prefix the differences repeat by the certified law.
    for n in from..=last {
        if certificate.degree_at(n - 1)?.checked_add(1)? != certificate.degree_at(n)? {
            return None;
        }
    }
    Some(certificate)
}

/// Prove nonzero constant coefficients by an exact scalar recurrence.
/// Derivative terms must vanish at zero, so no coefficient is silently omitted.
pub(super) fn nonzero_constants(
    record: &RecurrenceJson,
    from: usize,
) -> Option<LeadingCoefficientCertificate> {
    let (rule, first, seeds) = record.to_recurrence_parts().ok()?;
    let (mut scalar, _) = super::rule_with_x_free_left_hand_factor(&rule).ok()?;
    let at_zero = |p: &BivarPoly| {
        BivarPoly::from_coefficients_normalized(
            p.coeffs
                .iter()
                .map(|row| vec![vec_poly::coeff(row, 0)])
                .collect(),
        )
    };
    if scalar
        .terms
        .iter()
        .any(|term| term.deriv_order > 0 && !at_zero(&term.coeff).is_zero())
    {
        return None;
    }
    scalar.terms.retain(|term| term.deriv_order == 0);
    for term in &mut scalar.terms {
        term.coeff = at_zero(&term.coeff);
    }
    scalar.denominator = scalar.denominator.as_ref().map(at_zero);
    scalar.inhomogeneous = scalar.inhomogeneous.as_ref().map(at_zero);
    let seeds: Vec<Vec<BigRational>> = seeds
        .iter()
        .map(|row| vec![vec_poly::coeff(row, 0)])
        .collect();
    let record = RecurrenceJson::from_recurrence_rational(&scalar, first, &seeds, None);
    let certificate = degree_certificate(&record)?;
    (certificate.holds_from? <= from && certificate.sign_at(from).is_some()).then_some(certificate)
}

pub(super) fn metadata(
    declarations: &[&str],
    degrees: Vec<LeadingCoefficientCertificate>,
    constants: Option<LeadingCoefficientCertificate>,
    hypotheses: Vec<String>,
) -> Box<LiftInterlacingTransfer> {
    Box::new(LiftInterlacingTransfer {
        lean_declarations: declarations.iter().map(|s| (*s).into()).collect(),
        degree_certificates: degrees,
        constant_coefficient_certificate: constants.map(Box::new),
        hypotheses,
    })
}

/// Keep the old HKO conclusion unless successor image degrees fix its order.
pub(super) fn fixed_preserver(
    summary: &LiftSummary,
    parent: Result<&LiftInterlacing, String>,
    preserver: String,
    input: &RecurrenceJson,
) -> LiftInterlacingClaim {
    let mut claim = fixed_preserver_interlacing(summary, parent, preserver);
    if let LiftInterlacingClaim::Proved { interlacing } = &mut claim {
        if let Some(degrees) = successor_degrees(input, interlacing.from_index) {
            interlacing.orientation = InterlacingOrientation::PreviousInterlacesNext;
            interlacing.justification.premises.push(
                "successor output degrees fix the orientation: Q_(n-1) << Q_n at every covered index".into(),
            );
            interlacing.justification.transfer = Some(metadata(
                &["RealRooted.PreservesRealRootedOrZero.interlaces",
                  "RealRooted.interlaces_map_of_pencil"],
                vec![degrees], None,
                vec!["the image rows are nonzero and have successor degrees at every consecutive index; the same global linear preserver acts on both parent rows".into()],
            ));
        }
    }
    claim
}

/// Build a weak transfer only after the existing identity and side conditions pass.
pub(super) fn transferred_claim(
    summary: &LiftSummary,
    parent: &LiftInterlacing,
    orientation: InterlacingOrientation,
    evidence: Box<LiftInterlacingTransfer>,
) -> LiftInterlacingClaim {
    if !matches!(summary.identity, LiftIdentityStatus::Certified { .. })
        || !matches!(summary.parent_proof, LiftParentProof::Proved { .. })
        || summary
            .hypotheses
            .iter()
            .any(|h| h.status != LiftHypothesisStatus::Verified)
    {
        return LiftInterlacingClaim::not_proved(
            "the all-index identity or transfer hypotheses are not proved",
        );
    }
    let theorem = evidence.lean_declarations.join("; ");
    let mut claim = proved_claim(
        summary,
        parent,
        RowInterlacing::Weak,
        orientation,
        &theorem,
        evidence.hypotheses.clone(),
    );
    if let LiftInterlacingClaim::Proved { interlacing } = &mut claim {
        interlacing.justification.transfer = Some(evidence);
    }
    claim
}

/// Recover a proved pair property without changing the chosen discovery route.
pub(super) fn discovery_interlacing(
    report: &RecurrenceProofDiscoveryReport,
) -> Option<LiftInterlacing> {
    if report.status != "proved" {
        return None;
    }
    let mut attempts: Vec<_> = report.attempts.iter().collect();
    if let Some(tail) = &report.tail {
        attempts.extend(&tail.attempts);
    }
    if let Some(search) = &report.alternative_representations {
        for alternative in search
            .alternatives
            .iter()
            .filter(|a| a.proved && a.equality.is_certified())
        {
            attempts.extend(&alternative.attempts);
        }
    }
    attempts
        .into_iter()
        .find_map(|attempt| {
            super::parent_proof::theorem_report_interlacing(attempt, attempt.theorem.as_str())
        })
        .or_else(|| {
            report
                .catalog_lifts
                .as_ref()?
                .conclusion
                .as_ref()?
                .interlacing
                .clone()
        })
        .or_else(|| {
            report
                .closure_lifts
                .as_ref()?
                .conclusion
                .as_ref()?
                .interlacing
                .clone()
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::implicit::algebra::q;

    fn record(coefficients: Vec<Vec<BigRational>>, seed: Vec<BigRational>) -> RecurrenceJson {
        let rule = Recurrence {
            terms: vec![RecurrenceTerm::new(
                1,
                0,
                BivarPoly::from_coefficients_normalized(coefficients),
            )],
            denominator: None,
            inhomogeneous: None,
        };
        RecurrenceJson::from_recurrence_rational(&rule, 0, &[seed], None)
    }

    #[test]
    fn degree_evidence_rejects_zero_rows_and_degree_jumps() {
        let linear = record(vec![vec![q(1), q(1)]], vec![q(1)]);
        assert!(successor_degrees(&linear, 1).is_some());
        let jump = record(vec![vec![q(1), q(0), q(1)]], vec![q(1)]);
        assert!(successor_degrees(&jump, 1).is_none());
        let zero = record(vec![vec![q(1), q(1)]], vec![q(0)]);
        assert!(successor_degrees(&zero, 1).is_none());
    }

    #[test]
    fn a_late_zero_of_the_constant_multiplier_is_not_finite_evidence() {
        let nonzero = record(vec![vec![q(1), q(1)], vec![q(1)]], vec![q(1)]);
        assert!(nonzero_constants(&nonzero, 0).is_some());
        let late_zero = record(vec![vec![q(-50), q(1)], vec![q(1)]], vec![q(1)]);
        assert!(nonzero_constants(&late_zero, 0).is_none());
        let (mut rule, first, seeds) = nonzero.to_recurrence_parts().unwrap();
        rule.terms.push(RecurrenceTerm::new(
            1,
            1,
            BivarPoly::from_coefficients_normalized(vec![vec![q(1)]]),
        ));
        let derivative = RecurrenceJson::from_recurrence_rational(&rule, first, &seeds, None);
        assert!(nonzero_constants(&derivative, 0).is_none());
    }
}
