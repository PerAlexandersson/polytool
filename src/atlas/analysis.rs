//! Certificate analysis: replay a certificate and read off the sequence it
//! is about, its claims and the premises it needs ([`analyze_certificate`]);
//! canonical definitions and the all-`n` agreement of two definitions with
//! the same rule ([`same_sequence_from`]).
//!
//! # What each certificate gives
//!
//! - **Theorem check** (a [`RealRootednessTheoremReport`], replayed by its
//!   `verify()` and required proved): every row from the conclusion's first
//!   index is real-rooted in the certified root region (the whole line when
//!   none is certified), and, when the conclusion states it, consecutive
//!   rows interlace from `interlacing_from` on (`P_(m-1) << P_m`, the
//!   orientation of [`crate::recurrence::RowInterlacing`]).  No premises.
//! - **Parent provers** (replayed by recomputing
//!   `prove_parent(record, first_index, no attestations, try_root_domain)`
//!   and comparing with the recorded proof, which must be proved by a
//!   polytool route, never by an attestation): every row from the record's
//!   first index is zero or real-rooted with coefficients of one sign,
//!   equivalently every root in `(-inf, 0]` (Descartes' rule of signs), plus
//!   the interlacing the route concludes.  No premises.
//! - **Lifts** (preserver images, basis maps, gamma polynomials; replayed
//!   by their `verify()`): accepted only when the identity is certified for
//!   every index, every hypothesis of the transfer theorem is verified and
//!   the summary does not report a failed parent.  The claim is the
//!   summary's input root region from the input's first index (and the
//!   input interlacing when proved); the premise is the parent record with
//!   the requirement "every row from `n0 + s` on is zero or real-rooted with
//!   every root in the transfer theorem's parent region".  The lift's own
//!   parent proof is not used: the premise must be met by an atlas node.
//! - **Basis-map theorem** (a [`BasisMapTheoremReport`], replayed and
//!   required proved): one [`AtlasClaim::BasisMapPreserver`] per conclusion,
//!   about the node whose rows are the map's basis.  No premises.
//! - **Definition equality** (a preserver-image certification of the
//!   identity substitution `x -> x` with index shift `0` and no gauge,
//!   whose identity is certified): the input definition equals the parent
//!   definition for every index `n >= n0` (the certificate's input first
//!   index).  No claims.
//! - **Obreschkoff combination**: see [`crate::atlas::obreschkoff`].

use super::error::AtlasError;
use super::rows::{nonpositive_ray, whole_line};
use super::types::{AtlasCertificate, AtlasClaim, AtlasPremiseRole};
use crate::recurrence::{
    prepare_rule, prove_parent, LiftHypothesisStatus, LiftIdentityStatus, LiftInterlacingClaim,
    LiftParentProof, LiftSummary, ParentRecord, PreserverImageCertificate,
    PreserverImageCertification, PreserverImageIdentityStatus, PreserverOperator, RecurrenceJson,
    RootInterval, VeroneseParentProofRoute, VeroneseParentProofStatus, VeroneseSectionGauge,
};
use crate::recurrence::{BasisMapIdentityStatus, BasisMapTheoremStatus, GammaLiftIdentityStatus};
use crate::recurrence::{InterlacingOrientation, RealRootednessTheoremReport};
use crate::BigRational;
use serde::{Deserialize, Serialize};

/// A premise a certificate needs: a node defining the same sequence as
/// `definition` (from an index at most the claim's first index) with a
/// proved claim at least as strong as `uses`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PremiseRequirement {
    /// The premise's definition as recorded in the certificate.
    pub definition: RecurrenceJson,
    /// The role.
    pub role: AtlasPremiseRole,
    /// The required property.
    pub uses: AtlasClaim,
}

/// What a replayed certificate proves.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CertificateAnalysis {
    /// The definition of the sequence the claims are about.
    pub input: RecurrenceJson,
    /// The claims.
    pub claims: Vec<AtlasClaim>,
    /// The premises.
    pub premises: Vec<PremiseRequirement>,
    /// For a definition equality: the definition the input is proved equal
    /// to (for every index from the input's first index on).
    pub equal_to: Option<RecurrenceJson>,
    /// The cited theorem with its reference.
    pub theorem: String,
    /// The RealRooted row route, when recorded.
    pub lean_route: Option<String>,
    /// RealRooted declarations, when recorded.
    pub lean_declarations: Vec<String>,
}

/// The canonical form of a definition: the rule after the left-hand factor
/// normalization of the lifts, with its first index and exactly the
/// identifying seeds (`Recurrence::generation_initial_count`), validated as
/// runnable (homogeneous, positive offsets, no `(-1)^n` factors, a leader
/// free of `x` and of derivatives that is nonzero at every target).
///
/// # Errors
///
/// [`AtlasError::Recurrence`] when the record does not parse or the rule is
/// not runnable in that sense.
pub fn canonical_definition(record: &RecurrenceJson) -> Result<RecurrenceJson, AtlasError> {
    let (rule, first_index, rows) = record.to_recurrence_parts()?;
    let prepared = prepare_rule(&rule, &rows, first_index, "sequence")?;
    Ok(RecurrenceJson::from_recurrence_rational(
        &prepared.rule,
        prepared.first_index,
        &prepared.seeds,
        None,
    ))
}

/// The first `count` rows of a definition (from its first index), trimmed.
///
/// # Errors
///
/// [`AtlasError::Recurrence`] when the record does not replay.
pub fn definition_rows(
    record: &RecurrenceJson,
    count: usize,
) -> Result<Vec<Vec<BigRational>>, AtlasError> {
    let (rule, first_index, seeds) = record.to_recurrence_parts()?;
    let rows = rule
        .generate_rows_rational(&seeds, first_index, count.max(seeds.len()))
        .map_err(|error| {
            AtlasError::Recurrence(crate::recurrence::RecurrenceError::InvalidInput(format!(
                "the definition cannot be replayed: {error}"
            )))
        })?;
    Ok(rows
        .into_iter()
        .take(count)
        .map(|mut row| {
            while row.last().is_some_and(num_traits::Zero::is_zero) {
                row.pop();
            }
            row
        })
        .collect())
}

/// The first index from which two definitions are proved to give the same
/// rows, or `None`.  Both must have the same rule (as recorded, after
/// [`canonical_definition`]); then from `start = max(first indices)` the
/// rows `P_start, ..., P_(start + K - 1)` (`K` the larger number of seeds)
/// are compared exactly, and every later row is computed by the same rule
/// from the same previous rows in both, so the sequences agree from
/// `start` on.
///
/// # Errors
///
/// [`AtlasError`] when a definition is not runnable.
pub fn same_sequence_from(
    a: &RecurrenceJson,
    b: &RecurrenceJson,
) -> Result<Option<usize>, AtlasError> {
    let a = canonical_definition(a)?;
    let b = canonical_definition(b)?;
    if a.recurrence != b.recurrence {
        return Ok(None);
    }
    let start = a.first_index.max(b.first_index);
    let seeds = a.initial_polynomials.len().max(b.initial_polynomials.len());
    let count = |record: &RecurrenceJson| start - record.first_index + seeds;
    let rows_a = definition_rows(&a, count(&a))?;
    let rows_b = definition_rows(&b, count(&b))?;
    let tail_a = &rows_a[start - a.first_index..];
    let tail_b = &rows_b[start - b.first_index..];
    Ok((tail_a == tail_b).then_some(start))
}

fn index_independent(region: Option<&RootInterval>) -> RootInterval {
    match region {
        Some(region) if region.is_index_independent() => region.clone(),
        _ => whole_line(),
    }
}

fn theorem_analysis(
    report: &RealRootednessTheoremReport,
) -> Result<CertificateAnalysis, AtlasError> {
    report.verify()?;
    let conclusion = match (&report.conclusion, report.is_proved()) {
        (Some(conclusion), true) => conclusion,
        _ => {
            return Err(AtlasError::VerificationFailed(
                "the theorem report is not proved".into(),
            ))
        }
    };
    let mut claims = vec![AtlasClaim::RealRooted {
        from_index: conclusion.first_index,
        region: index_independent(conclusion.root_region.as_ref()),
    }];
    if let (Some(from_index), Some(kind)) = (conclusion.interlacing_from, conclusion.interlacing) {
        claims.push(AtlasClaim::Interlacing {
            from_index,
            kind,
            orientation: InterlacingOrientation::PreviousInterlacesNext,
        });
    }
    Ok(CertificateAnalysis {
        input: report.recurrence.clone(),
        claims,
        premises: Vec::new(),
        equal_to: None,
        theorem: report.reference.clone(),
        lean_route: report.lean_route.clone(),
        lean_declarations: report.lean_declarations.clone(),
    })
}

fn route_name(route: VeroneseParentProofRoute) -> &'static str {
    match route {
        VeroneseParentProofRoute::PeriodicProduct => "periodic_product",
        VeroneseParentProofRoute::PaperTheorem => "paper_theorem",
        VeroneseParentProofRoute::RootDomain => "root_domain",
        VeroneseParentProofRoute::ThreeTermRecurrence => "three_term_recurrence",
        VeroneseParentProofRoute::LiuWang => "liu_wang",
        VeroneseParentProofRoute::LiuWangDerivativeLag => "liu_wang_derivative_lag",
        VeroneseParentProofRoute::MaWangInterval => "ma_wang_interval",
        VeroneseParentProofRoute::Attestation => "attestation",
    }
}

/// The name of a parent-prover route (`periodic_product`, `paper_theorem`,
/// `root_domain`, ...).
pub(crate) fn parent_route_name(route: VeroneseParentProofRoute) -> &'static str {
    route_name(route)
}

fn parent_provers_analysis(
    record: &RecurrenceJson,
    try_root_domain: bool,
    proof: &crate::recurrence::VeroneseParentProof,
) -> Result<CertificateAnalysis, AtlasError> {
    let (rule, first_index, seeds) = record.to_recurrence_parts()?;
    let recomputed = prove_parent(
        &ParentRecord {
            rule: &rule,
            seeds: &seeds,
            first_index,
        },
        first_index,
        &[],
        try_root_domain,
    )?;
    if &recomputed != proof {
        return Err(AtlasError::VerificationFailed(
            "the parent proof differs from its recomputation".into(),
        ));
    }
    let route = match (proof.status, proof.route) {
        (VeroneseParentProofStatus::Proved, Some(route))
            if route != VeroneseParentProofRoute::Attestation =>
        {
            route
        }
        _ => {
            return Err(AtlasError::VerificationFailed(
                "the parent provers did not prove the statement".into(),
            ))
        }
    };
    let mut claims = vec![AtlasClaim::RealRooted {
        from_index: first_index,
        region: nonpositive_ray(),
    }];
    if let Some(interlacing) = &proof.interlacing {
        claims.push(AtlasClaim::Interlacing {
            from_index: interlacing.from_index,
            kind: interlacing.kind,
            orientation: interlacing.orientation,
        });
    }
    Ok(CertificateAnalysis {
        input: record.clone(),
        claims,
        premises: Vec::new(),
        equal_to: None,
        theorem: format!(
            "polytool parent provers, route {}: {}",
            route_name(route),
            proof.statement
        ),
        lean_route: None,
        lean_declarations: Vec::new(),
    })
}

/// The analysis of a lift from its verified summary and the input record.
fn lift_analysis(
    summary: &LiftSummary,
    input: &RecurrenceJson,
    lean_route: Option<String>,
    lean_declarations: Vec<String>,
) -> Result<CertificateAnalysis, AtlasError> {
    if !matches!(summary.identity, LiftIdentityStatus::Certified { .. }) {
        return Err(AtlasError::VerificationFailed(
            "the identity of the lift is not certified".into(),
        ));
    }
    if summary
        .hypotheses
        .iter()
        .any(|h| h.status != LiftHypothesisStatus::Verified)
    {
        return Err(AtlasError::VerificationFailed(
            "a hypothesis of the transfer theorem is not verified".into(),
        ));
    }
    if matches!(summary.parent_proof, LiftParentProof::Failed { .. }) {
        return Err(AtlasError::VerificationFailed(
            "the lift reports a failing parent row".into(),
        ));
    }
    let parent_region = summary
        .root_transfer
        .parent_root_region
        .as_ref()
        .filter(|r| r.is_index_independent())
        .ok_or_else(|| {
            AtlasError::Unsupported("the transfer theorem names no parent root region".into())
        })?;
    if !summary
        .root_transfer
        .input_root_region
        .is_index_independent()
    {
        return Err(AtlasError::Unsupported(
            "the input root region depends on the index".into(),
        ));
    }
    let needed = summary.input_first_index as i64 + summary.parent.index_shift;
    let needed = usize::try_from(needed).map_err(|_| {
        AtlasError::InvalidRecord("the lift uses a parent row before index 0".into())
    })?;
    if input.first_index != summary.input_first_index {
        return Err(AtlasError::VerificationFailed(
            "the certificate's input does not start at the summary's first index".into(),
        ));
    }
    let mut claims = vec![AtlasClaim::RealRooted {
        from_index: summary.input_first_index,
        region: summary.root_transfer.input_root_region.clone(),
    }];
    if let LiftInterlacingClaim::Proved { interlacing } = &summary.input_interlacing {
        claims.push(AtlasClaim::Interlacing {
            from_index: interlacing.from_index,
            kind: interlacing.kind,
            orientation: interlacing.orientation,
        });
    }
    Ok(CertificateAnalysis {
        input: input.clone(),
        claims,
        premises: vec![PremiseRequirement {
            definition: summary.parent.recurrence.clone(),
            role: AtlasPremiseRole::Sequence,
            uses: AtlasClaim::RealRooted {
                from_index: needed,
                region: parent_region.clone(),
            },
        }],
        equal_to: None,
        theorem: summary.root_transfer.theorem.clone(),
        lean_route,
        lean_declarations,
    })
}

/// The certified input record of a preserver-image certification.
pub(crate) fn preserver_input(
    certification: &PreserverImageCertification,
) -> Result<&RecurrenceJson, AtlasError> {
    match &certification.identity {
        PreserverImageIdentityStatus::Certified { certificate } => Ok(match certificate.as_ref() {
            PreserverImageCertificate::OperatorDivision(c) => &c.input,
            PreserverImageCertificate::HadamardConjugation(c) => &c.input,
        }),
        _ => Err(AtlasError::VerificationFailed(
            "the identity of the preserver image is not certified".into(),
        )),
    }
}

/// Whether a preserver candidate is the identity `Q_n = P_n`.
pub(crate) fn is_identity_candidate(certification: &PreserverImageCertification) -> bool {
    let candidate = &certification.candidate;
    let identity_map = match &candidate.operator {
        PreserverOperator::Affine { slope, intercept } => {
            crate::recurrence::parse_rational_coeff(slope).is_ok_and(|s| s == super::rows::one())
                && crate::recurrence::parse_rational_coeff(intercept)
                    .is_ok_and(|b| b == BigRational::from_integer(0.into()))
        }
        _ => false,
    };
    identity_map && candidate.index_shift == 0 && candidate.gauge == VeroneseSectionGauge::None
}

/// Replay a certificate and read off what it proves (module documentation).
///
/// # Errors
///
/// [`AtlasError::VerificationFailed`] when the certificate does not replay
/// or proves nothing usable, [`AtlasError::Unsupported`] for a lift whose
/// transfer theorem is not stated through a parent root region, and the
/// errors of the replays.
pub fn analyze_certificate(
    certificate: &AtlasCertificate,
) -> Result<CertificateAnalysis, AtlasError> {
    match certificate {
        AtlasCertificate::TheoremCheck { report } => theorem_analysis(report),
        AtlasCertificate::ParentProvers {
            record,
            try_root_domain,
            proof,
        } => parent_provers_analysis(record, *try_root_domain, proof),
        AtlasCertificate::PreserverImage { certification } => {
            certification.verify()?;
            if is_identity_candidate(certification) {
                return Err(AtlasError::Unsupported(
                    "an identity map is a definition equality, not an image".into(),
                ));
            }
            let input = preserver_input(certification)?.clone();
            lift_analysis(
                &certification.lift_summary,
                &input,
                certification.lean_route.clone(),
                certification.lean_declarations.clone(),
            )
        }
        AtlasCertificate::BasisMap { certification } => {
            certification.verify()?;
            if !matches!(
                certification.identity,
                BasisMapIdentityStatus::Certified { .. }
            ) {
                return Err(AtlasError::VerificationFailed(
                    "the basis-map identity is not certified".into(),
                ));
            }
            lift_analysis(
                &certification.lift_summary,
                &certification.input,
                None,
                Vec::new(),
            )
        }
        AtlasCertificate::GammaPolynomial { certification } => {
            certification.verify()?;
            let input = match &certification.identity {
                GammaLiftIdentityStatus::Certified { certificate } => certificate.input.clone(),
                _ => {
                    return Err(AtlasError::VerificationFailed(
                        "the gamma identity is not certified".into(),
                    ))
                }
            };
            lift_analysis(&certification.lift_summary, &input, None, Vec::new())
        }
        AtlasCertificate::VeroneseSection {
            input,
            certification,
        } => {
            let candidate = &certification.candidate;
            if candidate.index_shift != 0 || candidate.gauge != VeroneseSectionGauge::None {
                return Err(AtlasError::Unsupported(
                    "only Veronese sections without index shift or gauge are recorded".into(),
                ));
            }
            let recomputed = super::grow::veronese_certification(
                input,
                &candidate.parent,
                candidate.modulus,
                candidate.residue,
            )?;
            if recomputed != **certification {
                return Err(AtlasError::VerificationFailed(
                    "the Veronese certification differs from its recomputation".into(),
                ));
            }
            certification.verify()?;
            if !matches!(
                certification.identity,
                crate::recurrence::VeroneseSectionIdentityStatus::Certified { .. }
            ) {
                return Err(AtlasError::VerificationFailed(
                    "the Veronese identity is not certified".into(),
                ));
            }
            let summary = certification.lift_summary.as_ref().ok_or_else(|| {
                AtlasError::VerificationFailed("the certification has no lift summary".into())
            })?;
            lift_analysis(summary, input, None, Vec::new())
        }
        AtlasCertificate::ObreschkoffCombination { certificate } => certificate.analyze(),
        AtlasCertificate::BasisMapTheorem { report } => {
            report.verify()?;
            if report.status != BasisMapTheoremStatus::Proved {
                return Err(AtlasError::VerificationFailed(
                    "the basis-map theorem is not proved".into(),
                ));
            }
            let claims = report
                .conclusions
                .iter()
                .map(|c| AtlasClaim::BasisMapPreserver {
                    input_interval: c.input_interval.to_string(),
                    output_interval: c.output_interval.to_string(),
                    statement: c.statement.clone(),
                })
                .collect();
            Ok(CertificateAnalysis {
                input: report.basis.clone(),
                claims,
                premises: Vec::new(),
                equal_to: None,
                theorem: report.theorem.clone(),
                lean_route: None,
                lean_declarations: Vec::new(),
            })
        }
        AtlasCertificate::DefinitionEquality { certification } => {
            certification.verify()?;
            if !is_identity_candidate(certification) {
                return Err(AtlasError::VerificationFailed(
                    "a definition equality must use the identity map".into(),
                ));
            }
            let input = preserver_input(certification)?.clone();
            if input.first_index != certification.input_first_index {
                return Err(AtlasError::VerificationFailed(
                    "the equality certificate's input does not start at its first index".into(),
                ));
            }
            Ok(CertificateAnalysis {
                input,
                claims: Vec::new(),
                premises: Vec::new(),
                equal_to: Some(certification.candidate.parent.clone()),
                theorem: "all-index identity Q_n = P_n certified by exact Ore-algebra division \
                          and identifying rows (polytool preserver-image certificate of the \
                          identity map)"
                    .into(),
                lean_route: None,
                lean_declarations: Vec::new(),
            })
        }
    }
}
