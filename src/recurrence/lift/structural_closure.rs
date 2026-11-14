//! Structural closures. Rows propose preimages; exact Ore identities
//! and registry proofs certify them. No sequence library is consulted.
//!
//! Weight and Euler multipliers use the Pólya--Schur multiplier-sequence
//! theorem (Pólya and Schur, J. Reine Angew. Math. 144 (1914), 89--113;
//! also Brändén, arXiv:1410.6601, Theorem 7.1). Constant differential
//! operators use the Hermite--Poulain theorem. Row-factor transfer uses the
//! union of roots of a product. Gamma transfer uses the reciprocal-pair
//! factorization explained in the `gamma_lift` module documentation.
//! The existing certifiers decide every side condition for all indices.
//! Forward gamma expansion instead transports a palindromic input to its
//! gamma recurrence and proves those rows before expanding back. It may
//! succeed before the expensive discovery lift searches, after direct checks.

use super::*;
use crate::vec_poly;
use crate::BigRational as Q;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};

const ROWS: usize = 12;
const MAX_DEPTH: usize = 2;
const MAX_DEGREE: usize = 24;

mod forward_gamma;
pub use forward_gamma::ForwardGammaCertification;

/// A named link in a structural closure proof; paths refer to this attempt.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ClosureProofStep {
    /// Mathematical step name.
    pub name: String,
    /// Identity, parent claim or cited closure theorem.
    pub statement: String,
    /// Relative path to the independently replayed certificate.
    pub certificate: String,
}

/// An existing exact identity certificate for an input-derived parent.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum DiscoveryClosureCertificate {
    /// Prove the transported gamma family and expand back to the input.
    ForwardGamma {
        /// All-index symmetry, transport and root-transfer certificate.
        certification: Box<ForwardGammaCertification>,
    },
    /// Coefficientwise removal of a fixed-menu weight.
    HadamardFactor {
        /// All-index conjugation certificate and side conditions.
        certification: Box<HadamardFactorLift>,
    },
    /// Ordinary product with an input-derived fixed factor to an affine power.
    RowFactor {
        /// All-index substitution/division certificate.
        certification: Box<RowFactorLift>,
    },
    /// Inverse gamma expansion with an inferred affine symmetry schedule.
    Gamma {
        /// All-index gamma transport and symmetry certificate.
        certification: Box<GammaLiftCertification>,
    },
    /// Inversion of a fixed differential or Euler preserver.
    PreserverImage {
        /// All-index Ore division certificate and operator hypotheses.
        certification: Box<PreserverImageCertification>,
    },
}
impl DiscoveryClosureCertificate {
    fn summary(&self) -> &LiftSummary {
        match self {
            Self::ForwardGamma { certification } => &certification.lift_summary,
            Self::HadamardFactor { certification } => &certification.lift_summary,
            Self::RowFactor { certification } => &certification.lift_summary,
            Self::Gamma { certification } => &certification.lift_summary,
            Self::PreserverImage { certification } => &certification.lift_summary,
        }
    }
    fn parent_proof(&self) -> &VeroneseParentProof {
        match self {
            Self::ForwardGamma { certification } => &certification.parent_proof,
            Self::HadamardFactor { certification } => &certification.parent_proof,
            Self::RowFactor { certification } => &certification.parent_proof,
            Self::Gamma { certification } => &certification.parent_proof,
            Self::PreserverImage { certification } => &certification.parent_proof,
        }
    }
    fn stage(&self) -> &'static str {
        match self {
            Self::ForwardGamma { .. } => "forward_gamma",
            Self::HadamardFactor { .. } => "hadamard_factor",
            Self::RowFactor { .. } => "row_factor",
            Self::Gamma { .. } => "gamma",
            Self::PreserverImage { .. } => "preserver_image",
        }
    }
}

/// Proof of the computed parent, with no attestation or stored-family lookup.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum StructuralClosureParentProof {
    /// Proof discovery on transported gamma rows, including certified
    /// alternative representations. Lift recursion is disabled inside it.
    Discovery {
        /// Replayed all-index gamma proof and its exact identities.
        proof: Box<RecurrenceProofDiscoveryReport>,
    },
    /// General theorem registry (including product induction and root domains).
    Registry {
        /// Replayed theorem attempts and conclusion.
        proof: Box<VeroneseParentProof>,
    },
    /// Another inverse closure within the fixed depth bound.
    Closure {
        /// Computed parent and its own chain certificate.
        search: Box<DiscoveryClosureSearch>,
    },
}

/// Nonzero, all-index conclusion of a structural closure chain.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StructuralClosureConclusion {
    /// Every row from this index is covered.
    pub from_index: usize,
    /// Certified root region.
    pub root_region: RootInterval,
    /// Human-readable conclusion.
    pub statement: String,
}

/// Exact nonvanishing argument for a closure image.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum StructuralClosureNonzero {
    /// Gamma expansion is injective and takes a nonzero gamma row to a
    /// nonzero row. Its positivity also follows from the positive gamma row.
    GammaExpansion {
        /// Gamma leading-sign certificate when the parent proof is a
        /// registry proof. A nested closure carries its own nonzero proof.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        gamma_leading: Option<Box<LeadingCoefficientCertificate>>,
        /// A direct theorem that also certifies nonzero, nonnegative gamma
        /// rows when the separate leading-coefficient menu cannot close.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        gamma_theorem: Option<Box<RealRootednessTheoremReport>>,
        /// Exact injectivity and positivity argument, replayed with the chain.
        reason: String,
    },
    /// The supplied recurrence has an all-index nonzero leading coefficient.
    InputLeading {
        /// Replayed leading-coefficient induction.
        certificate: Box<LeadingCoefficientCertificate>,
    },
    /// An injective fixed preserver applied to nonzero parent rows.
    InjectiveImage {
        /// Parent leading-coefficient induction.
        parent_leading: Box<LeadingCoefficientCertificate>,
        /// Exact polynomial injectivity argument.
        reason: String,
    },
}

/// One exact identity, its computed parent's proof, and a closure theorem.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DiscoveryClosureAttempt {
    /// General operation name.
    pub stage: String,
    /// Descriptive roles, never stored-family names.
    pub parents: Vec<String>,
    /// All-index identity and closure side conditions.
    pub certificate: DiscoveryClosureCertificate,
    /// Independent proof of the input-derived parent.
    pub parent_proof: StructuralClosureParentProof,
    /// Exact nonzero-row certificate, including injectivity when needed.
    pub nonzero: StructuralClosureNonzero,
    /// Named chain, derived again on replay.
    pub steps: Vec<ClosureProofStep>,
    /// Nonzero-row conclusion.
    pub conclusion: StructuralClosureConclusion,
}

/// Bounded structural search, reached after the previous discovery stages fail.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DiscoveryClosureSearch {
    /// Supplied recurrence, or its certified tail.
    pub recurrence: RecurrenceJson,
    /// Bound on identity candidates retained from each inverse operation.
    pub max_candidates: usize,
    /// Remaining closure depth, at most two.
    pub depth: usize,
    /// Identity candidates examined by the parent prover.
    pub candidates_tried: usize,
    /// Successful chains; the search stops at its first proof.
    pub attempts: Vec<DiscoveryClosureAttempt>,
    /// The first successful chain's conclusion.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub conclusion: Option<StructuralClosureConclusion>,
}
impl DiscoveryClosureSearch {
    /// Whether an unconditional chain was proved.
    pub fn is_proved(&self) -> bool {
        self.conclusion.is_some()
    }

    /// Replay inverse proposals, exact certificates, parent proofs and text.
    /// No saved verdict or finite row agreement is trusted.
    ///
    /// # Errors
    /// Returns a verification error for altered certificates, budgets or claims.
    pub fn verify(&self, options: &RealRootednessTheoremOptions) -> Result<(), RecurrenceError> {
        if self.depth > MAX_DEPTH
            || self.max_candidates > CatalogLiftOptions::default().max_candidates
        {
            return Err(RecurrenceError::VerificationFailed(
                "invalid structural search bounds".into(),
            ));
        }
        let expected = search(
            &self.recurrence,
            options,
            self.max_candidates,
            self.depth,
            None,
        )?;
        if *self != expected {
            return Err(RecurrenceError::VerificationFailed(
                "structural closure chain changed".into(),
            ));
        }
        Ok(())
    }
}

/// Compute inverse closures and prove their parents through the registry.
/// The menu consists of coefficient weights, input-derived row factors,
/// forward gamma transport, inverse gamma expansions, and inverses of `1 +/- D` and `theta + a`
/// (`a=1,2`). Recursion has depth two; no catalog or named generator is read.
/// The existing lift candidate bound also bounds each inverse stage; zero
/// disables it. A requested root region skips this whole-line discovery stage.
///
/// # Errors
/// Returns a typed recurrence error if the supplied rows cannot be generated.
pub fn find_discovery_closure_lifts(
    record: &RecurrenceJson,
    options: &RealRootednessTheoremOptions,
    bounds: &CatalogLiftOptions,
) -> Result<DiscoveryClosureSearch, RecurrenceError> {
    search(
        record,
        options,
        bounds
            .max_candidates
            .min(CatalogLiftOptions::default().max_candidates),
        MAX_DEPTH,
        None,
    )
}

/// A successful cheap forward gamma attempt may precede expensive lift
/// searches. Failed attempts leave no serialized search or changed budget.
pub(super) fn find_forward_gamma_closure(
    record: &RecurrenceJson,
    options: &RealRootednessTheoremOptions,
    bounds: &CatalogLiftOptions,
) -> Option<DiscoveryClosureSearch> {
    if !forward_gamma::eligible(record) || options.root_region.is_some() || !bounds.is_enabled() {
        return None;
    }
    let found = search(
        record,
        options,
        bounds
            .max_candidates
            .min(CatalogLiftOptions::default().max_candidates),
        MAX_DEPTH,
        Some("forward_gamma"),
    )
    .ok()?;
    found.is_proved().then_some(found)
}

fn positive_leading(c: &LeadingCoefficientCertificate) -> bool {
    c.is_proved()
        && c.initial_rows
            .iter()
            .all(|r| r.sign == Some(LeadingCoefficientSign::Positive))
        && c.degree_law.as_ref().is_some_and(|law| {
            law.leading_coefficient_signs
                .iter()
                .all(|s| *s == LeadingCoefficientSign::Positive)
        })
}

fn positive_closure(search: &DiscoveryClosureSearch) -> bool {
    search.attempts.last().is_some_and(|a| match &a.nonzero {
        StructuralClosureNonzero::InputLeading { certificate } => positive_leading(certificate),
        // The fixed operators used here have positive leading diagonal:
        // 1 +/- D leaves the leading coefficient, theta+a multiplies by k+a>0.
        StructuralClosureNonzero::InjectiveImage { parent_leading, .. } => {
            positive_leading(parent_leading)
        }
        StructuralClosureNonzero::GammaExpansion { .. } => true,
    })
}

// These two routes carry a stronger sign conclusion than the general parent
// interface: PF Riordan gives nonnegative rows, and Liu--Wang three-term
// proves positive leading signs after the sign normalization determined by
// its two induction seeds. Positive seeds make that normalization trivial.
// Its auxiliary coordinate maps are monic gauges or reversal; on rows with
// nonpositive roots they preserve the coefficient sign. Do not infer this
// extra conclusion from an arbitrary registry route's proved status.
fn positive_gamma_theorem(
    record: &RecurrenceJson,
    known: &VeroneseParentProof,
) -> Option<RealRootednessTheoremReport> {
    let theorem = known.route?.theorem()?;
    if !matches!(
        theorem,
        RealRootednessTheorem::LiuWangThreeTerm | RealRootednessTheorem::RationalPfRiordan
    ) {
        return None;
    }
    let (rule, first, seeds) = record.to_recurrence_parts().ok()?;
    let report =
        check_real_rootedness_theorem(theorem, &rule, first, &seeds, &Default::default()).ok()?;
    let conclusion = report.conclusion.as_ref()?;
    if !report.is_proved()
        || conclusion.first_index > first
        || conclusion.zero_residue_classes.is_some()
        || !super::preserver_image::region_contains(
            &super::lift_summary::nonpositive_ray(),
            conclusion.root_region.as_ref()?,
        )
    {
        return None;
    }
    let prefix = rule
        .generate_rows_rational(&seeds, first, seeds.len() + report.start_delay)
        .ok()?;
    if prefix
        .iter()
        .any(|r| vec_poly::is_zero(r) || r.iter().any(Signed::is_negative))
    {
        return None;
    }
    Some(report)
}

// Necessary finite PF tests reject proposals only. They never establish a
// parent theorem or an identity. Newton's inequalities cheaply reject many
// coefficient quotients before the bounded recurrence fit.
pub(super) fn plausibly_pf(rows: &[Vec<Q>]) -> bool {
    rows.iter().rev().take(2).all(|row| {
        let Some(sign) = row.iter().find(|c| !c.is_zero()).map(Signed::is_positive) else {
            return true;
        };
        if row.iter().any(|c| !c.is_zero() && c.is_positive() != sign) {
            return false;
        }
        let d = row.len().saturating_sub(1);
        (1..d).all(|k| {
            let lhs = &row[k] * &row[k] * Q::from_integer((k * (d - k)).into());
            let rhs = &row[k - 1] * &row[k + 1] * Q::from_integer(((k + 1) * (d - k + 1)).into());
            lhs >= rhs
        })
    })
}

fn fit(rows: &[Vec<Q>], first: usize) -> Option<RecurrenceJson> {
    if !plausibly_pf(rows) {
        return None;
    }
    let options = AdaptiveSearchOptions {
        max_rec_len: 2,
        max_var_deg: 2,
        max_idx_deg: 2,
        max_diff_deg: 1,
        max_denom_var_deg: 0,
        max_denom_idx_deg: 1,
        ..Default::default()
    };
    match find_recurrence_adaptive_rational_with_budget(
        rows,
        &options,
        AdaptiveSearchBudget::limited(32),
    ) {
        AdaptiveSearchOutcome::Found(found) => {
            let rule = found.recurrence.shift_index(1 - first as i64);
            let prepared = prepare_lift_rule(&rule, rows, first, "computed preimage").ok()?;
            Some(RecurrenceJson::from_recurrence_rational(
                &prepared.rule,
                first,
                &prepared.seeds,
                None,
            ))
        }
        _ => None,
    }
}

fn parent_proof(
    record: &RecurrenceJson,
    known: &VeroneseParentProof,
    options: &RealRootednessTheoremOptions,
    budget: usize,
    depth: usize,
) -> Option<StructuralClosureParentProof> {
    // The identity certifier already ran exactly this structural registry
    // policy. Reuse that proof rather than running the same cascade twice.
    if !known.conditional_on.is_empty()
        || known
            .attempts
            .iter()
            .any(|a| a.route == VeroneseParentProofRoute::Attestation)
    {
        return None;
    }
    if known.status == VeroneseParentProofStatus::Proved {
        return Some(StructuralClosureParentProof::Registry {
            proof: Box::new(known.clone()),
        });
    }
    if known.status == VeroneseParentProofStatus::Refuted {
        return None;
    }
    if depth > 1 {
        let nested = search(record, options, budget, depth - 1, None).ok()?;
        if nested.conclusion.as_ref().is_some_and(|c| {
            c.from_index <= record.first_index
                && super::preserver_image::region_contains(
                    &super::lift_summary::nonpositive_ray(),
                    &c.root_region,
                )
        }) {
            return Some(StructuralClosureParentProof::Closure {
                search: Box::new(nested),
            });
        }
    }
    None
}

fn add(
    search: &mut DiscoveryClosureSearch,
    certificate: DiscoveryClosureCertificate,
    options: &RealRootednessTheoremOptions,
) -> bool {
    let summary = certificate.summary();
    if !matches!(summary.identity, LiftIdentityStatus::Certified { .. })
        || summary
            .hypotheses
            .iter()
            .any(|h| h.status != LiftHypothesisStatus::Verified)
    {
        return false;
    }
    search.candidates_tried += 1;
    let parent_proof = match &certificate {
        DiscoveryClosureCertificate::ForwardGamma { certification }
            if certification.gamma_discovery.as_ref().is_some_and(|p| {
                forward_gamma::accepted_discovery(p, summary.parent.recurrence.first_index)
            }) =>
        {
            StructuralClosureParentProof::Discovery {
                proof: certification.gamma_discovery.as_ref().unwrap().clone(),
            }
        }
        _ => {
            let Some(proof) = parent_proof(
                &summary.parent.recurrence,
                certificate.parent_proof(),
                options,
                search.max_candidates,
                search.depth,
            ) else {
                return false;
            };
            proof
        }
    };
    let Ok((rule, first, seeds)) = search.recurrence.to_recurrence_parts() else {
        return false;
    };
    let leading = certify_leading_coefficient(&rule, first, &seeds, &Default::default()).ok();
    let forward = matches!(
        certificate,
        DiscoveryClosureCertificate::ForwardGamma { .. }
    );
    let (from, nonzero) = if let Some(leading) = leading
        .filter(|c| c.is_proved() && c.holds_from.is_some() && (!forward || positive_leading(c)))
    {
        (
            first.max(leading.holds_from.unwrap()),
            StructuralClosureNonzero::InputLeading {
                certificate: Box::new(leading),
            },
        )
    } else if forward {
        let (from, gamma_leading, gamma_theorem) = match &parent_proof {
            StructuralClosureParentProof::Closure { search: nested }
                if positive_closure(nested) =>
            {
                (nested.conclusion.as_ref().unwrap().from_index, None, None)
            }
            StructuralClosureParentProof::Registry { .. } => {
                let Ok((r, p0, s)) = summary.parent.recurrence.to_recurrence_parts() else {
                    return false;
                };
                let Ok(c) = certify_leading_coefficient(&r, p0, &s, &Default::default()) else {
                    return false;
                };
                if positive_leading(&c) {
                    let Some(from) = c.holds_from else {
                        return false;
                    };
                    (from, Some(Box::new(c)), None)
                } else if let Some(proof) =
                    positive_gamma_theorem(&summary.parent.recurrence, certificate.parent_proof())
                {
                    (
                        proof.conclusion.as_ref().unwrap().first_index,
                        None,
                        Some(Box::new(proof)),
                    )
                } else {
                    return false;
                }
            }
            _ => return false,
        };
        (first.max(from), StructuralClosureNonzero::GammaExpansion {
            gamma_leading,
                gamma_theorem,
            reason: "the gamma proof and its leading-sign/injectivity certificate give nonzero gamma rows with positive leading coefficients; roots <=0 imply nonnegative coefficients; gamma expansion is injective and its basis polynomials have nonnegative coefficients".into(),
        })
    } else if let DiscoveryClosureCertificate::PreserverImage { certification } = &certificate {
        // Only the fixed invertible operators constructed by this search reach
        // this branch: 1 +/- D has unchanged leading coefficient, and theta+a
        // has strictly positive diagonal k+a for every k>=0, a=1,2.
        let parent = &certification.candidate.parent;
        let Ok((rule, p0, seeds)) = parent.to_recurrence_parts() else {
            return false;
        };
        let Ok(leading) = certify_leading_coefficient(&rule, p0, &seeds, &Default::default())
        else {
            return false;
        };
        let Some(from) = leading.holds_from.filter(|_| leading.is_proved()) else {
            return false;
        };
        let reason = "the fixed differential operator has constant term one, or the Euler diagonal is k+a>0 (k>=0, a=1,2); the gauge is one".to_string();
        (
            first.max(from),
            StructuralClosureNonzero::InjectiveImage {
                parent_leading: Box::new(leading),
                reason,
            },
        )
    } else {
        return false;
    };
    let statement = format!("every P_n, n >= {from}, is nonzero and real-rooted by {} of its computed parent; nonvanishing follows from the recorded leading-coefficient/injectivity argument", certificate.stage());
    let conclusion = StructuralClosureConclusion {
        from_index: from,
        root_region: summary.root_transfer.input_root_region.clone(),
        statement,
    };
    let parent_statement = match &parent_proof {
        StructuralClosureParentProof::Discovery { proof } => proof.statement.clone(),
        StructuralClosureParentProof::Registry { proof } => proof.statement.clone(),
        StructuralClosureParentProof::Closure { search } => {
            search.conclusion.as_ref().unwrap().statement.clone()
        }
    };
    let mut steps = vec![
        ClosureProofStep {
            name: "All-index identity".into(),
            statement: summary.parent.relation.clone(),
            certificate: "certificate.certification.identity".into(),
        },
        ClosureProofStep {
            name: "Computed parent proof".into(),
            statement: parent_statement,
            certificate: "parent_proof".into(),
        },
        ClosureProofStep {
            name: "Closure theorem".into(),
            statement: summary
                .root_transfer
                .theorem
                .replace("Theorem 40 of P. Branden", "P. Branden"),
            certificate: "certificate.certification.lift_summary.root_transfer".into(),
        },
        ClosureProofStep {
            name: "Nonzero rows".into(),
            statement: format!("leading coefficient is nonzero for every n >= {from}"),
            certificate: "nonzero".into(),
        },
    ];
    if forward {
        steps[0].certificate = "certificate.certification.gamma_transport".into();
        steps[3].statement = format!("rows are nonzero for every n >= {from} by the recorded leading-coefficient or injective gamma expansion argument");
        steps.push(ClosureProofStep {
            name: "Gamma coefficients nonnegative".into(),
            statement: "the gamma proof puts all roots <=0; positive gamma leading sign follows from its nonzero chain, or from the input's positive leading sign and all-index symmetry (the lowest nonzero gamma coefficient equals the lowest input coefficient); hence every gamma coefficient is nonnegative".into(),
            certificate: "parent_proof and nonzero".into(),
        });
    }
    search.attempts.push(DiscoveryClosureAttempt {
        stage: certificate.stage().into(),
        parents: vec!["computed preimage".into()],
        certificate,
        parent_proof,
        nonzero,
        steps,
        conclusion: conclusion.clone(),
    });
    search.conclusion = Some(conclusion);
    true
}

fn inverse_derivative(row: &[Q], c: i64) -> Vec<Q> {
    let mut out = row.to_vec();
    if row.len() < 2 {
        return out;
    }
    for k in (0..row.len() - 1).rev() {
        let next = &out[k + 1] * Q::from_integer(((k + 1) as i64 * c).into());
        out[k] -= next;
    }
    vec_poly::trim(&out)
}

fn search(
    record: &RecurrenceJson,
    options: &RealRootednessTheoremOptions,
    budget: usize,
    depth: usize,
    only: Option<&str>,
) -> Result<DiscoveryClosureSearch, RecurrenceError> {
    let mut out = DiscoveryClosureSearch {
        recurrence: record.clone(),
        max_candidates: budget,
        depth,
        candidates_tried: 0,
        attempts: vec![],
        conclusion: None,
    };
    if budget == 0 || depth == 0 || options.root_region.is_some() {
        return Ok(out);
    }
    let (rule, first, seeds) = record.to_recurrence_parts()?;
    if seeds.len() > ROWS {
        return Ok(out);
    }
    let rows = rule
        .generate_rows_rational(&seeds, first, ROWS.max(seeds.len()))
        .map_err(|e| RecurrenceError::Unsupported(e.to_string()))?;
    if rows.iter().any(|r| r.len() > MAX_DEGREE + 1) {
        return Ok(out);
    }
    let enabled = |s| only.is_none_or(|v| v == s);
    // Forward expansion consumes one depth level. Inner parent searches do
    // not try another forward expansion, so it cannot recursively re-enter.
    if depth == MAX_DEPTH && enabled("forward_gamma") && forward_gamma::eligible(record) {
        if let Some(c) = forward_gamma::certify(record, &rows) {
            let mut trial = out.clone();
            if add(
                &mut trial,
                DiscoveryClosureCertificate::ForwardGamma {
                    certification: Box::new(c),
                },
                options,
            ) {
                return Ok(trial);
            }
        }
    }
    if enabled("hadamard_factor") {
        let settings = HadamardFactorLiftOptions {
            allow_attestations: false,
            try_root_domain: false,
            detection_rows: ROWS,
            parent_search_candidates: 32,
            max_factors: 33,
            ..Default::default()
        };
        if let Ok(report) = super::hadamard_factor::search::find_structural_hadamard_factor_lifts(
            &rule, &seeds, first, &settings,
        ) {
            for c in report.lifts.into_iter().take(budget) {
                if add(
                    &mut out,
                    DiscoveryClosureCertificate::HadamardFactor {
                        certification: Box::new(c),
                    },
                    options,
                ) {
                    return Ok(out);
                }
            }
        }
    }
    if enabled("row_factor") {
        let settings = RowFactorLiftOptions {
            allow_attestations: false,
            try_root_domain: false,
            ..Default::default()
        };
        if let Ok(report) = find_row_factor_lifts(&rule, &seeds, first, &settings) {
            for c in report.lifts.into_iter().take(budget) {
                if add(
                    &mut out,
                    DiscoveryClosureCertificate::RowFactor {
                        certification: Box::new(c),
                    },
                    options,
                ) {
                    return Ok(out);
                }
            }
        }
    }
    if enabled("gamma") {
        // Minimal affine majorants of twice the observed degree. Only the
        // subsequent all-index symmetry/transport certificate proves the schedule.
        let growth = rows
            .last()
            .map_or(0, |r| r.len().saturating_sub(1))
            .saturating_sub(rows.first().map_or(0, |r| r.len().saturating_sub(1)));
        let step = (2 * growth).div_ceil(rows.len() - 1).max(1);
        let degree = rows
            .iter()
            .enumerate()
            .map(|(i, r)| (2 * r.len().saturating_sub(2)).saturating_sub(step * i))
            .max()
            .unwrap_or(0);
        let schedule = DegreeSchedule::new(first, degree, step as i64);
        for scale in [Q::one(), Q::new(1.into(), 4.into())] {
            let transform = RecurrenceTransformation::Gamma {
                schedule: schedule.clone(),
            };
            let parent_rows: Option<Vec<_>> = rows
                .iter()
                .enumerate()
                .map(|(i, r)| {
                    let scaled: Vec<_> = r
                        .iter()
                        .enumerate()
                        .map(|(k, c)| c / scale.pow(k as i32))
                        .collect();
                    transform.inverse_row(first + i, &[scaled]).ok()
                })
                .collect();
            let Some(parent) = parent_rows.and_then(|r| fit(&r, first)) else {
                continue;
            };
            let candidate = GammaLiftCandidate {
                parent: GammaLiftParent {
                    name: "inverse gamma expansion".into(),
                    recurrence: parent,
                    symmetry_schedule: schedule.clone(),
                },
                index_shift: 0,
                scale: scale.to_string(),
                gauge: VeroneseSectionGauge::None,
            };
            if let Ok(c) = certify_gamma_lift(
                &rule,
                &seeds,
                first,
                &candidate,
                &GammaLiftOptions {
                    allow_attestations: false,
                    try_root_domain: false,
                    ..Default::default()
                },
            ) {
                if add(
                    &mut out,
                    DiscoveryClosureCertificate::Gamma {
                        certification: Box::new(c),
                    },
                    options,
                ) {
                    return Ok(out);
                }
            }
        }
    }
    if enabled("preserver_image") {
        for (derivative, a) in [(true, 1i64), (true, -1), (false, 1), (false, 2)] {
            let parent_rows: Vec<_> = rows
                .iter()
                .map(|r| {
                    if derivative {
                        inverse_derivative(r, a)
                    } else {
                        r.iter()
                            .enumerate()
                            .map(|(k, c)| c / Q::from_integer((k as i64 + a).into()))
                            .collect()
                    }
                })
                .collect();
            let Some(parent) = fit(&parent_rows, first) else {
                continue;
            };
            let constants = |a: i64| IndexRationalFunction::constant(&Q::from_integer(a.into()));
            let operator = if derivative {
                PreserverOperator::DerivativePolynomial {
                    coefficients: vec![constants(1i64), constants(a)],
                }
            } else {
                PreserverOperator::EulerPolynomial {
                    coefficients: vec![constants(a), constants(1)],
                }
            };
            let candidate = PreserverImageCandidate {
                parent_name: "inverse fixed preserver".into(),
                parent,
                index_shift: 0,
                gauge: VeroneseSectionGauge::None,
                operator,
            };
            if let Ok(c) = certify_preserver_image_lift(
                &rule,
                &seeds,
                first,
                &candidate,
                &PreserverImageOptions {
                    allow_attestations: false,
                    try_root_domain: false,
                    ..Default::default()
                },
            ) {
                if add(
                    &mut out,
                    DiscoveryClosureCertificate::PreserverImage {
                        certification: Box::new(c),
                    },
                    options,
                ) {
                    return Ok(out);
                }
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests;
