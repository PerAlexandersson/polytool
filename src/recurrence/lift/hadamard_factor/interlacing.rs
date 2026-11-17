//! Positive diagonal and reversal-conjugated diagonal interlacing transfers.

use super::factor::{HadamardFactor, HadamardFactorTerm};
use crate::recurrence::lift::interlacing_transfer::{
    metadata, nonzero_constants, successor_degrees, transferred_claim,
};
use crate::recurrence::*;
use num_traits::Signed;

/// Try a stronger theorem, retaining the caller's existing claim otherwise.
pub(in crate::recurrence::lift) fn diagonal_interlacing(
    summary: &LiftSummary,
    factor: &HadamardFactor,
    parent: Result<&LiftInterlacing, String>,
    fallback: LiftInterlacingClaim,
) -> LiftInterlacingClaim {
    let Ok(parent) = parent else { return fallback };
    let Some(from) = super::super::lift_summary::input_from_index(
        parent.from_index,
        summary.parent.index_shift,
        summary.input_first_index,
    ) else {
        return fallback;
    };
    let Some(parent_from) = i128::try_from(from)
        .ok()
        .and_then(|n| n.checked_add(i128::from(summary.parent.index_shift)))
        .and_then(|n| usize::try_from(n).ok())
    else {
        return fallback;
    };
    let Some(degrees) = successor_degrees(&summary.parent.recurrence, parent_from) else {
        return fallback;
    };
    let positive_fixed = |term: &HadamardFactorTerm| match term {
        HadamardFactorTerm::InverseFactorial => true,
        HadamardFactorTerm::InverseRisingFactorial { parameter } => {
            parse_rational_coeff(parameter).is_ok_and(|p| p.is_positive())
        }
        _ => false,
    };
    if !factor.terms.is_empty() && factor.terms.iter().all(positive_fixed) {
        let claim = transferred_claim(summary, parent,
            InterlacingOrientation::PreviousInterlacesNext,
            metadata(&["RealRooted.Interlaces.diagonalOperator_of_pos"],
                vec![degrees], None,
                vec!["one fixed strictly positive multiplier sequence acts on both rows; positive diagonal weights preserve every degree and nonzero row; independent nonzero row gauges move no zero".into()]),
        );
        return if claim.is_proved() { claim } else { fallback };
    }
    let Some(degree) = degrees.degree_at(parent_from) else {
        return fallback;
    };
    let offset = degree as i128 - from as i128;
    let mut binomial = false;
    let mut complement = false;
    for term in &factor.terms {
        match term {
            t if positive_fixed(t) => {}
            HadamardFactorTerm::InverseComplementFactorial if offset <= 0 => {
                complement = true;
            }
            HadamardFactorTerm::Binomial {
                index_offset,
                coefficient_offset,
            } if *coefficient_offset >= 0
                && i128::from(*index_offset) - offset - i128::from(*coefficient_offset) >= 0 =>
            {
                binomial = true;
            }
            _ => return fallback,
        }
    }
    if !binomial && !complement {
        return fallback;
    }
    if summary.root_transfer.parent_root_region
        != Some(super::super::lift_summary::nonpositive_ray())
    {
        return fallback;
    }
    // The degree law is affine of slope one by successor_degrees, hence
    // d_(n+s)=n+offset for all n in this domain, not just sampled rows.
    let Some(constants) = nonzero_constants(&summary.parent.recurrence, parent_from - 1) else {
        return fallback;
    };
    let declarations = if binomial {
        vec![
            "RealRooted.Interlaces.diagonalOperator_choose_shift",
            "RealRooted.Interlaces.diagonalOperator_twoSided_of_pos",
        ]
    } else {
        vec![
            "RealRooted.Interlaces.reverse_diagonalOperator_reverse_of_pos",
            "RealRooted.Interlaces.diagonalOperator_twoSided_of_pos",
        ]
    };
    let claim = transferred_claim(summary, parent,
        InterlacingOrientation::PreviousInterlacesNext,
        metadata(&declarations, vec![degrees], Some(constants), vec![
            format!("the certified parent degrees satisfy d_(n+s)=n+({offset}) throughout the domain; all shifted reciprocal-factorial arguments are nonnegative on 0<=k<=d"),
            "the parent proof gives real-rooted rows with coefficients of one sign; the scalar constant-coefficient induction excludes zero constants, so each row can be normalized to PF; nonzero row normalizations move no zero".into(),
            "each binomial weight factors as a nonzero row factorial times reciprocal factorials in k+a and d-k+b with a,b>=0; products and positive reciprocal rising-factorial factors are fixed positive multiplier sequences on each side".into(),
        ]),
    );
    if claim.is_proved() {
        claim
    } else {
        fallback
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::lift::lift_summary::{
        nonpositive_ray, parent_proof_summary, proved_parent_interlacing,
    };

    fn summary(id: &str, shift: i64) -> (LiftSummary, VeroneseParentProof) {
        let entry = crate::oeis::by_id(id).unwrap();
        let (rule, seeds) = entry.recurrence_parts().unwrap();
        let record = RecurrenceJson::from_recurrence_rational(
            &rule,
            entry.recurrence_first_index,
            &seeds,
            None,
        );
        let proof = prove_parent_real_rootedness(
            &ParentRecurrenceRecord {
                rule: &rule,
                seeds: &seeds,
                first_index: entry.recurrence_first_index,
            },
            entry.recurrence_first_index,
            &[],
            false,
        )
        .unwrap();
        let summary = LiftSummary::new(
            LiftOperation::HadamardFactor {
                factor: "test".into(),
            },
            entry.recurrence_first_index,
            LiftParent {
                recurrence: record,
                relation: "certified diagonal identity".into(),
                gauge: None,
                index_shift: shift,
            },
            LiftIdentityStatus::Certified {
                method: "test".into(),
            },
            LiftRootTransfer {
                theorem: "multiplier".into(),
                parent_hypothesis: "PF up to row sign".into(),
                parent_root_region: Some(nonpositive_ray()),
                input_root_region: nonpositive_ray(),
                input_conclusion: "real-rooted".into(),
            },
            vec![],
            parent_proof_summary(&proof),
        );
        (summary, proof)
    }

    fn claim(
        summary: &LiftSummary,
        proof: &VeroneseParentProof,
        terms: Vec<HadamardFactorTerm>,
    ) -> LiftInterlacingClaim {
        diagonal_interlacing(
            summary,
            &HadamardFactor { terms },
            proved_parent_interlacing(proof),
            LiftInterlacingClaim::not_proved("unchanged"),
        )
    }

    #[test]
    fn positive_fixed_and_cross_rank_transfers_have_replayable_evidence() {
        let (summary, proof) = summary("A008292", 0);
        for (terms, expected) in [
            (
                vec![HadamardFactorTerm::InverseFactorial],
                "diagonalOperator_of_pos",
            ),
            (
                vec![HadamardFactorTerm::InverseComplementFactorial],
                "reverse_diagonalOperator_reverse_of_pos",
            ),
            (
                vec![
                    HadamardFactorTerm::InverseFactorial,
                    HadamardFactorTerm::InverseComplementFactorial,
                ],
                "diagonalOperator_twoSided_of_pos",
            ),
            (
                vec![HadamardFactorTerm::Binomial {
                    index_offset: 0,
                    coefficient_offset: 1,
                }],
                "diagonalOperator_choose_shift",
            ),
            (
                vec![
                    HadamardFactorTerm::Binomial {
                        index_offset: 1,
                        coefficient_offset: 1,
                    },
                    HadamardFactorTerm::InverseFactorial,
                ],
                "diagonalOperator_choose_shift",
            ),
        ] {
            let result = claim(&summary, &proof, terms);
            let interlacing = result.proved().expect("transfer");
            assert_eq!(interlacing.kind, RowInterlacing::Weak);
            assert_eq!(
                interlacing.orientation,
                InterlacingOrientation::PreviousInterlacesNext
            );
            let transfer = interlacing.justification.transfer.as_ref().unwrap();
            assert!(transfer
                .lean_declarations
                .iter()
                .any(|name| name.ends_with(expected)));
            for certificate in &transfer.degree_certificates {
                certificate.verify().unwrap();
            }
            if let Some(certificate) = &transfer.constant_coefficient_certificate {
                certificate.verify().unwrap();
            }
        }
    }

    #[test]
    fn failed_cross_rank_hypotheses_retain_the_old_claim() {
        let (mut summary, mut proof) = summary("A008292", 0);
        let binomial = || {
            vec![HadamardFactorTerm::Binomial {
                index_offset: 0,
                coefficient_offset: 2,
            }]
        };
        assert!(!claim(&summary, &proof, binomial()).is_proved()); // top weight zero
        assert!(!claim(&summary, &proof, vec![HadamardFactorTerm::Factorial]).is_proved());
        assert!(!claim(
            &summary,
            &proof,
            vec![HadamardFactorTerm::InverseRisingFactorial {
                parameter: "-1/2".into()
            }]
        )
        .is_proved());
        summary.identity = LiftIdentityStatus::SupportedByRows { rows: 100 };
        assert!(!claim(&summary, &proof, vec![HadamardFactorTerm::InverseFactorial]).is_proved());
        summary.identity = LiftIdentityStatus::Certified {
            method: "test".into(),
        };
        summary.root_transfer.parent_root_region = None;
        assert!(!claim(
            &summary,
            &proof,
            vec![HadamardFactorTerm::InverseComplementFactorial]
        )
        .is_proved());
        proof.interlacing = None;
        assert!(!claim(&summary, &proof, vec![HadamardFactorTerm::InverseFactorial]).is_proved());
        let (summary, proof) = self::summary("A048993", 0);
        // PF Stirling rows have zero constants after the first row.
        assert!(!claim(
            &summary,
            &proof,
            vec![HadamardFactorTerm::InverseComplementFactorial]
        )
        .is_proved());
    }
}
