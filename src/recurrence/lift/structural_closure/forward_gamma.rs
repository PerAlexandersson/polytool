//! Forward gamma discovery: prove the transported gamma rows, then expand.
//!
//! For G(z)=g product(z-zeta), zeta<=0, and degree r<=d/2,
//! P(x)=g(1+x)^(d-2r) product(x-zeta(1+x)^2).
//! For zeta<0 each quadratic has discriminant 1-4zeta>0 and two
//! negative roots; zeta=0 contributes x. Removed powers of x add only
//! zero roots. No interlacing or simplicity transfers between rows.
//! The transport certifies symmetry for all indices, not just sampled rows.

use super::*;

/// Exact forward gamma identity and the initial registry proof of its rows.
/// The enclosing closure attempt records a registry or bounded structural
/// proof of the gamma family and replays that proof as part of its chain.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ForwardGammaCertification {
    /// Optional exact division by a common power of x.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub power_removal: Option<Box<RecurrenceTransformReport>>,
    /// All-index symmetry and the exact transported gamma recurrence.
    pub gamma_transport: Box<RecurrenceTransformReport>,
    /// Registry proof of the gamma recurrence; an unresolved proof may be
    /// strengthened by the enclosing bounded structural closure proof.
    pub parent_proof: VeroneseParentProof,
    /// Bounded proof discovery on the transported gamma rule, with lifts
    /// and affine-coordinate searches disabled to prevent forward recursion.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gamma_discovery: Option<Box<RecurrenceProofDiscoveryReport>>,
    /// Uniform identity and root-transfer summary.
    pub lift_summary: LiftSummary,
}

fn valuation(row: &[Q]) -> Option<usize> {
    row.iter().position(|c| !c.is_zero())
}

fn palindromic_after_power(row: &[Q]) -> bool {
    let Some(v) = valuation(row) else {
        return false;
    };
    let d = vec_poly::degree(row).unwrap();
    row[v..=d].iter().eq(row[v..=d].iter().rev())
}

/// Cheap screen before any row generation, algebra, or parent search.
pub(super) fn eligible(record: &RecurrenceJson) -> bool {
    let Ok((rule, _, seeds)) = record.to_recurrence_parts() else {
        return false;
    };
    !seeds.is_empty()
        && seeds.len() <= ROWS
        && rule.terms.len() <= 32
        && rule.terms.iter().all(|t| {
            t.offset <= ROWS
                && t.deriv_order <= 2
                && t.coeff.coeffs.len() <= 4
                && t.coeff.coeffs.iter().all(|r| r.len() <= MAX_DEGREE + 1)
        })
        && seeds
            .iter()
            .all(|r| r.len() <= MAX_DEGREE + 1 && palindromic_after_power(r))
}

pub(super) fn certify(
    record: &RecurrenceJson,
    rows: &[Vec<Q>],
) -> Option<ForwardGammaCertification> {
    if !rows.iter().all(|r| palindromic_after_power(r)) {
        return None;
    }
    let power = rows
        .iter()
        .map(|r| valuation(r))
        .collect::<Option<Vec<_>>>()?
        .into_iter()
        .min()?;
    let power_removal = if power == 0 {
        None
    } else {
        let report = transform_recurrence(&RecurrenceTransformRequest::from_source(
            RecurrenceTransformation::DividePower { power },
            record.clone(),
            rows.len(),
        ))
        .ok()?;
        if !report.reconstruction_verified
            || report
                .conditions
                .iter()
                .any(|c| c.status != "certified_all_indices")
        {
            return None;
        }
        Some(Box::new(report))
    };
    let source = match &power_removal {
        Some(r) => r.transported.as_ref()?.scalar.as_ref()?.clone(),
        None => record.clone(),
    };
    let centers = rows
        .iter()
        .map(|r| Some(valuation(r)? + vec_poly::degree(r)? - 2 * power))
        .collect::<Option<Vec<_>>>()?;
    let start_degree = *centers.first()?;
    let step = i64::try_from(*centers.get(1)?).ok()? - i64::try_from(start_degree).ok()?;
    if step < 0
        || !centers
            .iter()
            .enumerate()
            .all(|(i, &d)| i64::try_from(d).ok() == Some(start_degree as i64 + step * i as i64))
    {
        return None;
    }
    let schedule = DegreeSchedule::new(record.first_index, start_degree, step);
    let gamma_transport = transform_recurrence(&RecurrenceTransformRequest::from_source(
        RecurrenceTransformation::Gamma {
            schedule: schedule.clone(),
        },
        source,
        rows.len(),
    ))
    .ok()?;
    if !gamma_transport.reconstruction_verified
        || gamma_transport
            .conditions
            .iter()
            .any(|c| c.status != "certified_all_indices")
    {
        return None;
    }
    let transported = gamma_transport.transported.as_ref()?;
    if !transported.holds_on_transformed_rows {
        return None;
    }
    let gamma = transported.scalar.as_ref()?;
    let (rule, first, seeds) = gamma.to_recurrence_parts().ok()?;
    // Negative gamma coefficients are a rejection screen, never a proof.
    if gamma_transport.transformed_states.iter().any(|s| {
        s[0].iter()
            .any(|c| parse_rational_coeff(c).is_ok_and(|c| c.is_negative()))
    }) {
        return None;
    }
    let parent_proof = prove_parent_real_rootedness(
        &ParentRecurrenceRecord {
            rule: &rule,
            first_index: first,
            seeds: &seeds,
        },
        first,
        &[],
        true,
    )
    .ok()?;
    let root_region = super::super::lift_summary::nonpositive_ray();
    let gamma_discovery = if parent_proof.status == VeroneseParentProofStatus::NotProved {
        let stages = RecurrenceProofStageOptions {
            catalog_lifts: CatalogLiftOptions::disabled(),
            affine_coordinates: AffineCoordinateOptions::disabled(),
            ..Default::default()
        };
        discover_recurrence_proofs_with_stage_options(gamma, &Default::default(), &stages)
            .ok()
            .filter(|p| accepted_discovery(p, first))
            .map(Box::new)
    } else {
        None
    };
    let lift_summary = LiftSummary::new(
        LiftOperation::GammaExpansion { symmetry_schedule: schedule, removed_power: power },
        record.first_index,
        LiftParent { recurrence: gamma.clone(), index_shift: 0, gauge: None,
            relation: format!("P_n(x) = x^{power} (1+x)^d_n G_n(x/(1+x)^2), for every n >= {}", record.first_index) },
        LiftIdentityStatus::Certified { method: "gamma_recurrence_transport".into() },
        LiftRootTransfer {
            theorem: "Gamma expansion: G(z)=g product_i(z-zeta_i), zeta_i<=0, gives P(x)=x^power g(1+x)^(d-2r) product_i(x-zeta_i(1+x)^2); each factor has only nonpositive real roots (quadratic discriminant 1-4zeta_i>0); gamma criterion, RealRooted IsGammaExpansion.isRealRooted_and_hasRootsNonpos_iff".into(),
            parent_hypothesis: "every gamma row is nonzero with nonnegative coefficients and only real roots <= 0".into(),
            parent_root_region: Some(root_region.clone()), input_root_region: root_region,
            input_conclusion: "every input row is nonzero and real-rooted with roots <= 0".into(),
        },
        vec![], super::super::lift_summary::parent_proof_summary(&parent_proof),
    );
    Some(ForwardGammaCertification {
        power_removal,
        gamma_transport: Box::new(gamma_transport),
        parent_proof,
        gamma_discovery,
        lift_summary,
    })
}

pub(super) fn accepted_discovery(proof: &RecurrenceProofDiscoveryReport, first: usize) -> bool {
    let in_ray = |r: &RealRootednessTheoremReport| {
        r.is_proved()
            && r.conclusion.as_ref().is_some_and(|c| {
                c.first_index <= first
                    && c.root_region.as_ref().is_some_and(|r| {
                        super::super::preserver_image::region_contains(
                            &super::super::lift_summary::nonpositive_ray(),
                            r,
                        )
                    })
            })
    };
    proof.status == "proved"
        && proof.real_rooted_from_index.is_some_and(|n| n <= first)
        && (proof.attempts.iter().any(in_ray)
            || proof.alternative_representations.as_ref().is_some_and(|s| {
                s.alternatives
                    .iter()
                    .any(|a| a.proved && a.attempts.iter().any(in_ray))
            }))
}
