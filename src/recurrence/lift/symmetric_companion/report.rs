//! Report types of [`find_symmetric_companion_lifts`] and their
//! `polytool.recurrence-symmetric-companion.v1` JSON form.
//!
//! Every exact rational is an integer or `num/den` string; recurrences are
//! `polytool.recurrence.v1` records with their seeds.

use super::super::lift_summary::symmetric_companion_root_transfer;
use super::certificate::{
    decide_on_domain, sign_hypotheses, SymmetricCompanionCertificate,
    SymmetricCompanionCheckStatus, SymmetricCompanionIdentityStatus,
    SymmetricCompanionRealRootednessStatus, SymmetricCompanionVerification,
};
use super::operators::CompanionParameters;
use super::source::PreparedSource;
use crate::integer_range_sign::{IndexSignCondition, IndexSignDecision};
use crate::recurrence::implicit::algebra::poly_text;
use crate::recurrence::*;
use serde::{Deserialize, Serialize};

/// JSON schema tag of a serialized [`SymmetricCompanionLiftReport`].
pub const SYMMETRIC_COMPANION_LIFT_SCHEMA: &str = "polytool.recurrence-symmetric-companion.v1";

/// Where the symmetry schedule came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SymmetryScheduleOrigin {
    /// Supplied in the options.
    Supplied,
    /// Inferred from the searched rows (`d = valuation + degree` of the first
    /// nonzero row); a candidate that each certificate checks on the seeds.
    InferredFromRows,
}

/// The symmetry schedule `d_n = d_(n0) + (n - n0)`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymmetryScheduleRecord {
    /// `n0`.
    pub first_index: usize,
    /// `d_(n0)`; `None` when no schedule could be inferred.
    pub symmetry_degree_at_first_index: Option<usize>,
    /// Supplied or inferred; `None` when no schedule could be inferred.
    pub origin: Option<SymmetryScheduleOrigin>,
    /// The rule as text.
    pub rule: String,
}

impl SymmetryScheduleRecord {
    pub(super) fn new(
        first_index: usize,
        symmetry_degree_at_first_index: Option<usize>,
        origin: Option<SymmetryScheduleOrigin>,
    ) -> Self {
        let rule = match symmetry_degree_at_first_index {
            Some(d0) => format!("d_n = {d0} + (n - {first_index})"),
            None => "no schedule".into(),
        };
        Self {
            first_index,
            symmetry_degree_at_first_index,
            origin,
            rule,
        }
    }
}

/// Outcome of one ansatz degree.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SymmetricCompanionAttemptOutcome {
    /// One or more rational parameter points were found (candidates).
    CandidateFound,
    /// No parameters with `epsilon != 0` (at `n0`) in this ansatz fit the
    /// searched rows (zero linearized kernel, or no rational point).
    NoSolution,
    /// The constant terms rule out `R_n(0) = 0` for every `a_n` of this degree.
    ConstantTermsIncompatible,
    /// The rows do not determine the parameters (too few nonzero constant
    /// terms, or an exact system with infinitely many solutions).
    Underdetermined,
    /// A budget (`max_linearized_monomials` or the kernel limits) was reached.
    BudgetExceeded,
}

/// One ansatz degree of the search.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymmetricCompanionAnsatzAttempt {
    /// Degree in `n` of `a_n`, `c_n`, `epsilon_n`.
    pub parameter_degree_in_n: usize,
    /// Whether `a_n` was fitted from the constant terms (otherwise it was an
    /// unknown of the linearized system).
    pub a_from_constant_terms: bool,
    /// Number of unknown parameter coefficients.
    pub unknowns: usize,
    /// Number of monomials (columns) of the linearized system.
    pub monomials: usize,
    /// Number of distinct nonzero equations.
    pub equations: usize,
    /// Dimension of the exact kernel, when computed.
    pub kernel_dimension: Option<usize>,
    /// The outcome.
    pub outcome: SymmetricCompanionAttemptOutcome,
    /// Details of the outcome.
    pub detail: String,
}

/// Summary of the bounded search.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymmetricCompanionSearchSummary {
    /// The method, as text.
    pub method: String,
    /// Rows used by the search.
    pub fit_rows: usize,
    /// The attempts in increasing degree.
    pub attempts: Vec<SymmetricCompanionAnsatzAttempt>,
    /// No certified lift was found and some attempt was underdetermined or
    /// over budget, so a larger budget might succeed.
    pub exhausted_budget: bool,
}

/// One candidate lift with its certificate and the verifier's statuses.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymmetricCompanionLift {
    /// The replayable certificate.
    pub certificate: SymmetricCompanionCertificate,
    /// The result of [`verify_symmetric_companion_certificate`] on it.
    pub verification: SymmetricCompanionVerification,
    /// The companion as a formula in `A_n`, `A_(n+1)`.
    pub companion_formula: String,
    /// The parent rule as text.
    pub parent_rule: String,
    /// Comparison of the parent rule with the source rule.
    pub comparison: String,
    /// The uniform [`LiftSummary`] of this lift
    /// ([`SymmetricCompanionLift::recompute_lift_summary`]); always present in
    /// reports of this version, optional in JSON for records written before
    /// the field existed.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lift_summary: Option<LiftSummary>,
}

impl SymmetricCompanionLift {
    pub(super) fn new(
        certificate: SymmetricCompanionCertificate,
        verification: SymmetricCompanionVerification,
        source: &PreparedSource,
        parameters: &CompanionParameters,
    ) -> Self {
        let n = |p: &[BigRational]| format!("({})", poly_text(p, "n"));
        let d = parameters.schedule();
        let companion_formula = format!(
            "B_n = [A_(n+1) - U_n A_n] / (epsilon_n x), U_n = [a_n + (a_n + c_n d_n) x] + c_n x \
             (1 - x) D, with a_n = {}, c_n = {}, d_n = {}, epsilon_n = {}",
            n(&parameters.a),
            n(&parameters.c),
            n(&d),
            n(&parameters.epsilon)
        );
        let parent_rule = format!(
            "P_(n+1) = ({} + {} x) P_n + {} x (1 - x) D P_n, P_{} = A_{} + x B_{}",
            n(&parameters.a),
            n(&parameters.b()),
            n(&parameters.c),
            source.first_index,
            source.first_index,
            source.first_index
        );
        let source_derivatives = source
            .recurrence
            .terms
            .iter()
            .map(|t| t.deriv_order)
            .max()
            .unwrap_or(0);
        let comparison = format!(
            "parent: order 1, derivative order 1 (a first-order rule); source: order {}, \
             derivative order {source_derivatives}",
            source.order
        );
        let lift_summary = Some(lift_summary(&certificate, &verification));
        Self {
            certificate,
            verification,
            companion_formula,
            parent_rule,
            comparison,
            lift_summary,
        }
    }

    /// The [`LiftSummary`] of this lift, computed from the certificate and
    /// the verification (see the free function behind it).
    pub fn recompute_lift_summary(&self) -> LiftSummary {
        lift_summary(&self.certificate, &self.verification)
    }

    fn rank(&self) -> usize {
        match (
            self.verification.identity_status,
            self.verification.real_rootedness_status,
        ) {
            (_, SymmetricCompanionRealRootednessStatus::TheoremApplies) => 0,
            (SymmetricCompanionIdentityStatus::Certified, _) => 1,
            (SymmetricCompanionIdentityStatus::Candidate, _) => 2,
            (SymmetricCompanionIdentityStatus::Unsupported, _) => 3,
            (SymmetricCompanionIdentityStatus::Rejected, _) => 4,
        }
    }
}

/// The report of [`find_symmetric_companion_lifts`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymmetricCompanionLiftReport {
    /// [`SYMMETRIC_COMPANION_LIFT_SCHEMA`].
    pub schema: String,
    /// The source specification: the scalar recurrence and its seeds.
    pub source: RecurrenceJson,
    /// The symmetry schedule.
    pub symmetry_schedule: SymmetryScheduleRecord,
    /// The search summary.
    pub search: SymmetricCompanionSearchSummary,
    /// The candidates, best first: theorem applies, certified, candidate,
    /// unsupported, rejected.
    pub lifts: Vec<SymmetricCompanionLift>,
    /// Reasons why no (further) candidate exists.
    pub obstructions: Vec<String>,
}

impl SymmetricCompanionLiftReport {
    pub(super) fn new(source: RecurrenceJson, fit_rows: usize) -> Self {
        Self {
            schema: SYMMETRIC_COMPANION_LIFT_SCHEMA.into(),
            source,
            symmetry_schedule: SymmetryScheduleRecord::new(0, None, None),
            search: SymmetricCompanionSearchSummary {
                method: "fit a_n from the constant terms (R_n(0) = 0); linearize the eliminated \
                         second coupled equation G_n = 0 (every monomial carries epsilon) in the \
                         monomials of the remaining coefficients; recover a point from a \
                         one-dimensional exact kernel; certify each candidate"
                    .into(),
                fit_rows,
                attempts: Vec::new(),
                exhausted_budget: false,
            },
            lifts: Vec::new(),
            obstructions: Vec::new(),
        }
    }

    pub(super) fn finish(&mut self) {
        self.lifts.sort_by_key(SymmetricCompanionLift::rank);
        let any_certified = self.lifts.iter().any(|l| l.verification.is_certified());
        self.search.exhausted_budget = !any_certified
            && self.search.attempts.iter().any(|a| {
                matches!(
                    a.outcome,
                    SymmetricCompanionAttemptOutcome::Underdetermined
                        | SymmetricCompanionAttemptOutcome::BudgetExceeded
                )
            });
        if self.lifts.is_empty() {
            for attempt in &self.search.attempts {
                self.obstructions.push(format!(
                    "degree {}: {}",
                    attempt.parameter_degree_in_n, attempt.detail
                ));
            }
        }
    }

    /// The lifts whose identity is certified.
    pub fn certified_lifts(&self) -> impl Iterator<Item = &SymmetricCompanionLift> {
        self.lifts.iter().filter(|l| l.verification.is_certified())
    }

    /// The first lift for which the preservation theorem applies.
    pub fn proving_lift(&self) -> Option<&SymmetricCompanionLift> {
        self.lifts.iter().find(|l| l.verification.theorem_applies())
    }

    /// Serialize as pretty JSON.
    pub fn to_json_string(&self) -> Result<String, RecurrenceError> {
        serde_json::to_string_pretty(self)
            .map_err(|e| RecurrenceError::InvalidInput(format!("cannot serialize the report: {e}")))
    }
}

/// The [`LiftSummary`] of a certificate and its verification.
///
/// `certified` maps to a certified identity (method
/// `coupled_system_elimination`), `candidate` (no obligation fails and the
/// exact rows agree) to `supported_by_rows`, `unsupported` and `rejected` to
/// the matching statuses with the first such obligation as the reason.  The
/// parent statement is proved (route `interlacing_preservation_induction`)
/// exactly when the preservation theorem applies, and open otherwise.  The
/// hypotheses are the theorem's; the three sign hypotheses carry the exact
/// decision of [`crate::decide_index_polynomial_sign`].  The summary is
/// proved exactly when [`SymmetricCompanionVerification::theorem_applies`].
pub(super) fn lift_summary(
    certificate: &SymmetricCompanionCertificate,
    verification: &SymmetricCompanionVerification,
) -> LiftSummary {
    use SymmetricCompanionCheckStatus as Status;
    let n0 = certificate.source.first_index;
    let d0 = certificate.symmetry_degree_at_first_index;
    let first_with = |status: Status| {
        verification
            .identities
            .iter()
            .find(|check| check.status == status)
            .map(|check| format!("{}: {}", check.name, check.evidence))
            .unwrap_or_default()
    };
    let identity = match verification.identity_status {
        SymmetricCompanionIdentityStatus::Certified => LiftIdentityStatus::Certified {
            method: "coupled_system_elimination".into(),
        },
        SymmetricCompanionIdentityStatus::Candidate => LiftIdentityStatus::SupportedByRows {
            rows: verification.self_checked_rows,
        },
        SymmetricCompanionIdentityStatus::Unsupported => LiftIdentityStatus::Unsupported {
            reason: first_with(Status::Unsupported),
        },
        SymmetricCompanionIdentityStatus::Rejected => LiftIdentityStatus::Refuted {
            index: None,
            reason: first_with(Status::Failed),
        },
    };
    // The sign decisions behind the verifier's sign hypotheses.
    let decisions: Vec<(&str, IndexSignDecision, String)> = match certificate.parameters.values() {
        Ok((a, _, c, _, epsilon)) => {
            let parameters = CompanionParameters::new(a, c, epsilon, n0, d0);
            sign_hypotheses(&parameters)
                .into_iter()
                .map(|(name, _, p)| {
                    let decision = decide_on_domain(&p, n0, IndexSignCondition::Nonnegative);
                    (name, decision, format!("p(n) = {}", poly_text(&p, "n")))
                })
                .collect()
        }
        Err(_) => Vec::new(),
    };
    let hypotheses = verification
        .hypotheses
        .iter()
        .map(|check| {
            let status = match check.status {
                Status::Verified => LiftHypothesisStatus::Verified,
                Status::Failed => LiftHypothesisStatus::Failed { witness: None },
                Status::Unresolved | Status::Unsupported => LiftHypothesisStatus::Unresolved,
            };
            let decided = decisions
                .iter()
                .find(|(name, _, _)| *name == check.name)
                .map(|(_, decision, function)| {
                    LiftHypothesis::from_sign_decision(
                        &check.name,
                        &check.statement,
                        decision.clone(),
                        function,
                    )
                });
            match decided {
                // The verifier made the same decision; it adds the witness and
                // the replayable certificate.
                Some(decided)
                    if std::mem::discriminant(&decided.status)
                        == std::mem::discriminant(&status) =>
                {
                    decided
                }
                _ => LiftHypothesis {
                    name: check.name.clone(),
                    statement: check.statement.clone(),
                    status,
                    evidence: check.evidence.clone(),
                    sign_decision: None,
                },
            }
        })
        .collect();
    let parent_proof = if verification.theorem_applies() {
        LiftParentProof::Proved {
            route: "interlacing_preservation_induction".into(),
            statement: format!(
                "the preservation theorem carries the pair B_n <= A_n with nonnegative \
                 coefficients from n0 = {n0} on, so every P_n = A_n + x B_n (n >= {n0}) is zero \
                 or real-rooted"
            ),
        }
    } else {
        LiftParentProof::Open {
            statement: "the preservation induction needs the certified identity and every \
                        hypothesis; it is not established"
                .into(),
        }
    };
    LiftSummary::new(
        LiftOperation::SymmetricCompanion {
            symmetry_degree_at_first_index: d0,
        },
        n0,
        LiftParent {
            recurrence: certificate.parent.clone(),
            relation: format!(
                "A_n is the symmetric part at d_n = {d0} + (n - {n0}) of P_n = A_n + x B_n for \
                 every n >= {n0}"
            ),
            gauge: None,
            index_shift: 0,
        },
        identity,
        symmetric_companion_root_transfer(n0),
        hypotheses,
        parent_proof,
    )
}
