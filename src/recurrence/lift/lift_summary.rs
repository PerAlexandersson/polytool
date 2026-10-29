//! [`LiftSummary`]: one uniform view of every lift operation, embedded as the
//! additive field `lift_summary` in the operations' own reports.
//!
//! A lift writes the input as `Q_n = T(P_n)` for an operation `T` that
//! preserves real-rootedness under a root-location (or interlacing)
//! hypothesis.  The summary keeps the statuses apart:
//!
//! - [`LiftOperation`]: which `T`, with its parameters;
//! - [`LiftParent`]: the parent rule and seeds, gauge and index shift;
//! - [`LiftIdentityStatus`]: the all-index identity `Q_n = T(P_n)`
//!   (certified, supported by rows only, unsupported or refuted);
//! - [`LiftRootTransfer`]: the transfer theorem, the root region the parent
//!   needs and the region it gives for the input, as exact
//!   [`RootInterval`]s;
//! - [`LiftHypothesis`] entries: the remaining hypotheses of the transfer
//!   theorem, sign hypotheses with the exact decision of
//!   [`crate::integer_range_sign`];
//! - [`LiftParentProof`]: whether the parent statement is proved (and by
//!   which route), open, failed or attested;
//! - [`LiftVerdict`], derived by [`lift_verdict`] alone: proved exactly when
//!   the identity is certified, the parent is proved and every hypothesis is
//!   verified;
//! - `lean_replay`: always `"not_attempted"`.
//!
//! The summary is computed from each operation's detailed report, which
//! stays authoritative; a verdict that is not `proved` never says anything
//! against real-rootedness of the input.
//!
//! The mappings of the Veronese lifts are here
//! ([`VeroneseSectionLift::lift_summary`],
//! [`VeroneseSectionCertification::recompute_lift_summary`]); the
//! symmetric-companion mapping is in that module
//! ([`SymmetricCompanionLift::recompute_lift_summary`]).

use super::*;
use crate::integer_range_sign::exact_serde;
use crate::integer_range_sign::{
    decide_index_rational_function_sign, IndexSignCondition, IndexSignDecision, IntegerRange,
};
use num_bigint::BigInt;
use num_traits::Zero;
use serde::{Deserialize, Serialize};

/// The lift operation `T` of `Q_n = T(P_n)` and its parameters.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
#[non_exhaustive]
pub enum LiftOperation {
    /// The coefficient section `V_(r, j)(P)(y) = sum_k [x^(r k + j)] P(x) y^k`.
    VeroneseSection {
        /// The modulus `r >= 2`.
        modulus: usize,
        /// The residue `0 <= j < r`.
        residue: usize,
    },
    /// The symmetric part `A_n` of `P_n = A_n + x B_n` with the schedule
    /// `d_n = d_(n0) + (n - n0)`.
    SymmetricCompanion {
        /// `d_(n0)`.
        symmetry_degree_at_first_index: usize,
    },
    /// The Hadamard product `T(P)_n = lambda(n) sum_k u(n, k) [x^k] P_n x^k`
    /// with a known factor `u`.
    HadamardFactor {
        /// The factor `u(n, k)` as text, for example `binom(n+1, k+1)`.
        factor: String,
    },
    /// The row factorization `T(S)_n = F(x)^(e(n)) S_n` with
    /// `e(n) = exponent_slope n + exponent_offset`.
    RowFactor {
        /// `F`, ascending exact coefficients.
        factor: Vec<String>,
        /// The slope of `e(n)`.
        exponent_slope: usize,
        /// The offset of `e(n)`.
        exponent_offset: i64,
    },
    /// The rescaled gamma polynomial `T(P)_n = lambda(n) Gamma_n(c z)`, where
    /// `P_n = sum_j Gamma_(n,j) x^j (1 + x)^(d_n - 2j)` is symmetric at `d_n`.
    GammaPolynomial {
        /// The symmetry indices `d_n` of the parent rows.
        symmetry_schedule: DegreeSchedule,
        /// The scale `c != 0`, an exact rational string.
        scale: String,
    },
}

/// The parent of a lift.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LiftParent {
    /// The parent rule with its first index and seeds
    /// (`polytool.recurrence.v1`).
    pub recurrence: RecurrenceJson,
    /// The relation between input and parent, as text.
    pub relation: String,
    /// The normalisation `lambda(n)` of `Q_n = lambda(n) T(P_(n+s))`, for an
    /// operation that has one.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gauge: Option<VeroneseSectionGauge>,
    /// The index shift `s` (default 0).
    #[serde(default)]
    pub index_shift: i64,
}

/// Status of the all-index identity `Q_n = T(P_n)`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case", deny_unknown_fields)]
pub enum LiftIdentityStatus {
    /// Proved for every `n >= n0`; the operation's report holds the
    /// certificate.
    Certified {
        /// The certificate kind (`identifying_rows`, `forward_division`,
        /// `coupled_system_elimination`).
        method: String,
    },
    /// The compared rows agree but the identity is not proved: finite
    /// evidence only.
    SupportedByRows {
        /// Number of rows compared from `n0`.
        rows: usize,
    },
    /// The candidate leaves the class the certificate supports (for example
    /// a vanishing coupling); nothing is claimed either way.
    Unsupported {
        /// Why.
        reason: String,
    },
    /// The identity (or one of its obligations) fails exactly.
    Refuted {
        /// The absolute index of the first differing row, when the failure is
        /// a row.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        index: Option<usize>,
        /// What fails.
        reason: String,
    },
}

/// The transfer theorem of a lift and the root regions it relates.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LiftRootTransfer {
    /// The theorem, with its reference.
    pub theorem: String,
    /// What the parent must satisfy, in words.
    pub parent_hypothesis: String,
    /// The region that must contain every root of every parent row used
    /// (nonzero rows real-rooted); `None` when the hypothesis is not a parent
    /// root region (an interlacing hypothesis, for example).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub parent_root_region: Option<RootInterval>,
    /// The region that then contains every root of every nonzero input row,
    /// which is real-rooted.
    pub input_root_region: RootInterval,
    /// The conclusion for the input rows, in words.
    pub input_conclusion: String,
}

/// Status of one hypothesis of the transfer theorem.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum LiftHypothesisStatus {
    /// Established for every index of the domain.
    Verified,
    /// Refuted exactly.
    Failed {
        /// The smallest index where an index statement fails (exact integer
        /// string), when there is one.
        #[serde(
            default,
            skip_serializing_if = "Option::is_none",
            with = "exact_serde::optional_integer"
        )]
        witness: Option<BigInt>,
    },
    /// Neither established nor refuted (for example over a budget).
    Unresolved,
}

/// One hypothesis of the transfer theorem with its status and evidence.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LiftHypothesis {
    /// Stable snake-case name.
    pub name: String,
    /// What is claimed.
    pub statement: String,
    /// The status.
    pub status: LiftHypothesisStatus,
    /// How the status was obtained.
    pub evidence: String,
    /// For a sign hypothesis, the decision of
    /// [`crate::decide_index_rational_function_sign`] (with its replayable
    /// certificate when it holds).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sign_decision: Option<IndexSignDecision>,
}

impl LiftHypothesis {
    /// A sign hypothesis whose status is the checker's decision.
    pub(super) fn from_sign_decision(
        name: &str,
        statement: &str,
        decision: IndexSignDecision,
        function: &str,
    ) -> Self {
        let status = match &decision {
            IndexSignDecision::Holds { .. } => LiftHypothesisStatus::Verified,
            IndexSignDecision::FailsAt { index, .. } => LiftHypothesisStatus::Failed {
                witness: Some(index.clone()),
            },
            IndexSignDecision::Undecided { .. } => LiftHypothesisStatus::Unresolved,
        };
        Self {
            name: name.into(),
            statement: statement.into(),
            status,
            evidence: decision.describe(function),
            sign_decision: Some(decision),
        }
    }
}

/// The parent statement (every parent row used is zero or real-rooted in
/// the required region) and how it was settled.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case", deny_unknown_fields)]
pub enum LiftParentProof {
    /// Proved for every parent row used.
    Proved {
        /// The route (`periodic_product`, `paper_theorem`, `root_domain`,
        /// `interlacing_preservation_induction`).
        route: String,
        /// The statement as proved.
        statement: String,
    },
    /// No route settled it.
    Open {
        /// What is open.
        statement: String,
    },
    /// A parent row fails the statement exactly.
    Failed {
        /// The absolute index of the first failing parent row, when known.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        index: Option<usize>,
        /// What fails.
        statement: String,
    },
    /// Asserted by caller-supplied attestations only.
    Attested {
        /// The attestation labels.
        labels: Vec<String>,
        /// The attested statement.
        statement: String,
    },
}

impl LiftParentProof {
    /// The statement text of any variant.
    pub fn statement(&self) -> &str {
        match self {
            Self::Proved { statement, .. }
            | Self::Open { statement }
            | Self::Failed { statement, .. }
            | Self::Attested { statement, .. } => statement,
        }
    }
}

/// The verdict of a lift, derived by [`lift_verdict`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LiftVerdict {
    /// Identity certified, parent proved, every hypothesis verified: every
    /// input row `Q_n`, `n >= n0`, is zero or real-rooted in
    /// [`LiftRootTransfer::input_root_region`].
    Proved,
    /// As [`LiftVerdict::Proved`], but the parent statement is attested by
    /// the caller: the conclusion holds if the attestations do.
    ConditionalOnAttestation,
    /// Nothing is refuted, but something is not established.
    Open,
    /// The identity is refuted, the parent fails, or a hypothesis fails:
    /// this route gives no proof.  It says nothing against real-rootedness
    /// of the input.
    RouteFailed,
}

/// The value of [`LiftSummary::lean_replay`] in this version.
pub const LIFT_LEAN_REPLAY_NOT_ATTEMPTED: &str = "not_attempted";

/// The uniform summary of one lift (module documentation).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LiftSummary {
    /// The operation and its parameters.
    pub operation: LiftOperation,
    /// The input's first index `n0`.
    pub input_first_index: usize,
    /// The parent.
    pub parent: LiftParent,
    /// The identity status.
    pub identity: LiftIdentityStatus,
    /// The transfer theorem and the root regions.
    pub root_transfer: LiftRootTransfer,
    /// The hypotheses of the transfer theorem besides the parent statement.
    pub hypotheses: Vec<LiftHypothesis>,
    /// The parent statement.
    pub parent_proof: LiftParentProof,
    /// [`lift_verdict`] of the fields above.
    pub verdict: LiftVerdict,
    /// Always [`LIFT_LEAN_REPLAY_NOT_ATTEMPTED`] in this version.
    pub lean_replay: String,
}

/// The verdict of a lift from its statuses alone: [`LiftVerdict::Proved`]
/// exactly when the identity is certified, the parent is proved and every
/// hypothesis is verified; [`LiftVerdict::ConditionalOnAttestation`] when
/// the parent is attested instead; [`LiftVerdict::RouteFailed`] when the
/// identity is refuted, the parent fails or a hypothesis fails; otherwise
/// [`LiftVerdict::Open`].
pub fn lift_verdict(
    identity: &LiftIdentityStatus,
    parent_proof: &LiftParentProof,
    hypotheses: &[LiftHypothesis],
) -> LiftVerdict {
    let hypothesis_failed = hypotheses
        .iter()
        .any(|h| matches!(h.status, LiftHypothesisStatus::Failed { .. }));
    if matches!(identity, LiftIdentityStatus::Refuted { .. })
        || matches!(parent_proof, LiftParentProof::Failed { .. })
        || hypothesis_failed
    {
        return LiftVerdict::RouteFailed;
    }
    let certified = matches!(identity, LiftIdentityStatus::Certified { .. });
    let verified = hypotheses
        .iter()
        .all(|h| h.status == LiftHypothesisStatus::Verified);
    match parent_proof {
        LiftParentProof::Proved { .. } if certified && verified => LiftVerdict::Proved,
        LiftParentProof::Attested { .. } if certified && verified => {
            LiftVerdict::ConditionalOnAttestation
        }
        _ => LiftVerdict::Open,
    }
}

impl LiftSummary {
    /// Assemble a summary; the verdict is [`lift_verdict`] of the parts.
    pub fn new(
        operation: LiftOperation,
        input_first_index: usize,
        parent: LiftParent,
        identity: LiftIdentityStatus,
        root_transfer: LiftRootTransfer,
        hypotheses: Vec<LiftHypothesis>,
        parent_proof: LiftParentProof,
    ) -> Self {
        let verdict = lift_verdict(&identity, &parent_proof, &hypotheses);
        Self {
            operation,
            input_first_index,
            parent,
            identity,
            root_transfer,
            hypotheses,
            parent_proof,
            verdict,
            lean_replay: LIFT_LEAN_REPLAY_NOT_ATTEMPTED.into(),
        }
    }

    /// Whether the summary proves that every input row is zero or
    /// real-rooted in the input root region.
    pub fn is_proved(&self) -> bool {
        self.verdict == LiftVerdict::Proved
    }

    /// Whether the recorded verdict is [`lift_verdict`] of the recorded
    /// statuses and the Lean replay field has its only value.
    pub fn is_consistent(&self) -> bool {
        self.verdict == lift_verdict(&self.identity, &self.parent_proof, &self.hypotheses)
            && self.lean_replay == LIFT_LEAN_REPLAY_NOT_ATTEMPTED
    }
}

// ---------------------------------------------------------------------------
// Root regions and transfer theorems
// ---------------------------------------------------------------------------

/// The closed ray `(-infinity, 0]`.
pub(super) fn nonpositive_ray() -> RootInterval {
    RootInterval {
        lower: None,
        upper: Some(RootEndpoint {
            value: "0".into(),
            per_index: None,
            closed: true,
        }),
    }
}

/// The Veronese transfer theorem ([`find_veronese_section_lifts`]).
pub(super) fn veronese_root_transfer(conclusion: String) -> LiftRootTransfer {
    LiftRootTransfer {
        theorem: "Veronese sections of Polya frequency sequences (Aissen--Schoenberg--Whitney--\
                  Edrei): if P is real-rooted with coefficients of one sign, equivalently nonzero \
                  with all roots in (-inf, 0], every section V_(r,j)(P) is zero or real-rooted \
                  with roots in (-inf, 0]"
            .into(),
        parent_hypothesis: "every parent row used is zero or real-rooted with all roots in \
                            (-inf, 0] (coefficients of one sign)"
            .into(),
        parent_root_region: Some(nonpositive_ray()),
        input_root_region: nonpositive_ray(),
        input_conclusion: conclusion,
    }
}

/// The interlacing-preservation route of the symmetric-companion lifts.
pub(super) fn symmetric_companion_root_transfer(first_index: usize) -> LiftRootTransfer {
    LiftRootTransfer {
        theorem: "interlacing preservation of the symmetric-companion step with the \
                  interlacing-decomposition theorem of Branden--Solus (arXiv:1808.04141, \
                  Theorem 2.6): if a, b have nonnegative coefficients and b <= a (directed weak \
                  interlacing), then a + x b is real-rooted"
            .into(),
        parent_hypothesis: format!(
            "an interlacing hypothesis, not a parent root region: A_({first_index}), \
             B_({first_index}) have nonnegative coefficients with B <= A, d_n >= 1, a_n >= 0 and \
             0 <= delta_n <= c_n for every n >= {first_index}, so the step preserves the \
             interlacing pair (A_n, B_n) of P_n = A_n + x B_n"
        ),
        parent_root_region: None,
        input_root_region: nonpositive_ray(),
        input_conclusion: format!(
            "every A_n, n >= {first_index}, is zero or real-rooted (as are B_n and P_n); its \
             roots lie in (-inf, 0] because its coefficients are nonnegative"
        ),
    }
}

// ---------------------------------------------------------------------------
// Veronese v1
// ---------------------------------------------------------------------------

impl VeroneseSectionLift {
    /// The [`LiftSummary`] of this lift.  The identity is certified by the
    /// identifying-rows argument of [`find_veronese_section_lifts`]; the
    /// parent statement is proved exactly for
    /// [`VeroneseSectionLiftStatus::ProvedProductParent`], failed for a
    /// failing parent row and open otherwise; there are no further
    /// hypotheses.
    pub fn lift_summary(&self) -> LiftSummary {
        let n0 = self.parent_first_index;
        let parent_proof = match self.status {
            VeroneseSectionLiftStatus::ProvedProductParent => LiftParentProof::Proved {
                route: "periodic_product".into(),
                statement: self.status.statement(),
            },
            VeroneseSectionLiftStatus::ParentRealRootednessOpen => LiftParentProof::Open {
                statement: self.status.statement(),
            },
            VeroneseSectionLiftStatus::ParentRowZero { index }
            | VeroneseSectionLiftStatus::ParentRowMixedSigns { index }
            | VeroneseSectionLiftStatus::ParentRowNotRealRooted { index } => {
                LiftParentProof::Failed {
                    index: Some(index),
                    statement: self.status.statement(),
                }
            }
        };
        LiftSummary::new(
            LiftOperation::VeroneseSection {
                modulus: self.modulus,
                residue: self.residue,
            },
            n0,
            LiftParent {
                recurrence: self.parent_recurrence_json(),
                relation: format!(
                    "Q_n = V_({},{})(P_n) for every n >= {n0}",
                    self.modulus, self.residue
                ),
                gauge: None,
                index_shift: 0,
            },
            LiftIdentityStatus::Certified {
                method: "identifying_rows".into(),
            },
            veronese_root_transfer(format!(
                "every Q_n, n >= {n0}, is zero or real-rooted with roots in (-inf, 0]"
            )),
            Vec::new(),
            parent_proof,
        )
    }
}

// ---------------------------------------------------------------------------
// Veronese v2
// ---------------------------------------------------------------------------

/// The hypothesis `lambda(n)` finite and nonzero for every `n >= n0`, with
/// the exact sign decision of the gauge polynomials.
fn gauge_hypothesis(gauge: &VeroneseSectionGauge, n0: usize) -> LiftHypothesis {
    let name = "gauge_finite_nonzero";
    let statement = "lambda(n) is finite and nonzero for every n >= n0";
    let parse = |values: &[String]| -> Option<Vec<BigRational>> {
        values
            .iter()
            .map(|value| parse_rational_coeff(value).ok())
            .collect()
    };
    let unresolved = |evidence: &str| LiftHypothesis {
        name: name.into(),
        statement: statement.into(),
        status: LiftHypothesisStatus::Unresolved,
        evidence: evidence.into(),
        sign_decision: None,
    };
    match gauge {
        VeroneseSectionGauge::None => LiftHypothesis {
            name: name.into(),
            statement: statement.into(),
            status: LiftHypothesisStatus::Verified,
            evidence: "lambda(n) = 1".into(),
            sign_decision: None,
        },
        VeroneseSectionGauge::Closed {
            numerator,
            denominator,
        } => match (parse(numerator), parse(denominator)) {
            (Some(numerator), Some(denominator)) => LiftHypothesis::from_sign_decision(
                name,
                statement,
                decide_index_rational_function_sign(
                    &numerator,
                    &denominator,
                    &IntegerRange::starting_at(n0),
                    IndexSignCondition::Nonzero,
                ),
                "lambda(n) = numerator(n) / denominator(n)",
            ),
            _ => unresolved("the gauge polynomials do not parse"),
        },
        VeroneseSectionGauge::Ratio {
            numerator,
            denominator,
            initial_value,
        } => {
            let initial = parse_rational_coeff(initial_value).ok();
            match (parse(numerator), parse(denominator), initial) {
                (_, _, Some(initial)) if initial.is_zero() => LiftHypothesis {
                    name: name.into(),
                    statement: statement.into(),
                    status: LiftHypothesisStatus::Failed {
                        witness: Some(BigInt::from(n0)),
                    },
                    evidence: "lambda(n0) = 0".into(),
                    sign_decision: None,
                },
                (Some(numerator), Some(denominator), Some(_)) => {
                    let mut hypothesis = LiftHypothesis::from_sign_decision(
                        name,
                        statement,
                        decide_index_rational_function_sign(
                            &numerator,
                            &denominator,
                            &IntegerRange::starting_at(n0 + 1),
                            IndexSignCondition::Nonzero,
                        ),
                        "lambda(n) / lambda(n - 1) = numerator(n) / denominator(n)",
                    );
                    hypothesis.evidence = format!(
                        "lambda(n0) = {initial_value} != 0 and {}",
                        hypothesis.evidence
                    );
                    hypothesis
                }
                _ => unresolved("the gauge does not parse"),
            }
        }
    }
}

impl VeroneseSectionCertification {
    /// The [`LiftSummary`] of this certification, computed from the record:
    /// the identity status maps one to one, the parent proof by its status
    /// (a proved status counts only with a non-attestation route, exactly as
    /// in [`VeroneseSectionCertification::verdict`]), and the gauge
    /// hypothesis is decided exactly.  The summary is proved exactly when
    /// [`VeroneseSectionVerdict::Proved`] is the verdict.
    pub fn recompute_lift_summary(&self) -> LiftSummary {
        let candidate = &self.candidate;
        let n0 = self.input_first_index;
        let (r, j, s) = (candidate.modulus, candidate.residue, candidate.index_shift);
        let identity = match &self.identity {
            VeroneseSectionIdentityStatus::Certified { .. } => LiftIdentityStatus::Certified {
                method: "forward_division".into(),
            },
            VeroneseSectionIdentityStatus::SupportedByRows { rows } => {
                LiftIdentityStatus::SupportedByRows { rows: *rows }
            }
            VeroneseSectionIdentityStatus::Refuted { index } => LiftIdentityStatus::Refuted {
                index: Some(*index),
                reason: format!("Q_{index} differs from the candidate section at n = {index}"),
            },
        };
        let proof = &self.parent_proof;
        let statement = proof.statement.clone();
        let parent_proof = match proof.status {
            VeroneseParentProofStatus::Proved => match proof.route {
                Some(route) if route != VeroneseParentProofRoute::Attestation => {
                    LiftParentProof::Proved {
                        route: parent_route_label(route).into(),
                        statement,
                    }
                }
                _ => LiftParentProof::Open { statement },
            },
            VeroneseParentProofStatus::ConditionalOnAttestation => LiftParentProof::Attested {
                labels: proof.conditional_on.clone(),
                statement,
            },
            VeroneseParentProofStatus::NotProved => LiftParentProof::Open { statement },
            VeroneseParentProofStatus::Refuted => LiftParentProof::Failed {
                index: proof.failing_row_index,
                statement,
            },
        };
        LiftSummary::new(
            LiftOperation::VeroneseSection {
                modulus: r,
                residue: j,
            },
            n0,
            LiftParent {
                recurrence: candidate.parent.clone(),
                relation: format!("Q_n = lambda(n) V_({r},{j})(P_(n{s:+})) for every n >= {n0}"),
                gauge: Some(candidate.gauge.clone()),
                index_shift: s,
            },
            identity,
            veronese_root_transfer(format!(
                "every Q_n = lambda(n) V_({r},{j})(P_(n{s:+})), n >= {n0}, is zero or \
                 real-rooted with roots in (-inf, 0]"
            )),
            vec![gauge_hypothesis(&candidate.gauge, n0)],
            parent_proof,
        )
    }
}

/// The snake-case label of a parent route.
fn parent_route_label(route: VeroneseParentProofRoute) -> &'static str {
    match route {
        VeroneseParentProofRoute::PeriodicProduct => "periodic_product",
        VeroneseParentProofRoute::PaperTheorem => "paper_theorem",
        VeroneseParentProofRoute::RootDomain => "root_domain",
        VeroneseParentProofRoute::ThreeTermRecurrence => "three_term_recurrence",
        VeroneseParentProofRoute::LiuWang => "liu_wang",
        VeroneseParentProofRoute::Attestation => "attestation",
    }
}

// ---------------------------------------------------------------------------
// Hadamard-factor and row-factor lifts
// ---------------------------------------------------------------------------

/// The [`LiftParentProof`] of a parent-prover record: a proved status counts
/// only with a non-attestation route (as for the Veronese certification).
pub(super) fn parent_proof_summary(proof: &VeroneseParentProof) -> LiftParentProof {
    let statement = proof.statement.clone();
    match proof.status {
        VeroneseParentProofStatus::Proved => match proof.route {
            Some(route) if route != VeroneseParentProofRoute::Attestation => {
                LiftParentProof::Proved {
                    route: parent_route_label(route).into(),
                    statement,
                }
            }
            _ => LiftParentProof::Open { statement },
        },
        VeroneseParentProofStatus::ConditionalOnAttestation => LiftParentProof::Attested {
            labels: proof.conditional_on.clone(),
            statement,
        },
        VeroneseParentProofStatus::NotProved => LiftParentProof::Open { statement },
        VeroneseParentProofStatus::Refuted => LiftParentProof::Failed {
            index: proof.failing_row_index,
            statement,
        },
    }
}

/// The whole real line, the root region when nothing more is known.
pub(super) fn real_line() -> RootInterval {
    RootInterval {
        lower: None,
        upper: None,
    }
}

// ---------------------------------------------------------------------------
// Gamma lifts
// ---------------------------------------------------------------------------

/// The parts of a [`GammaLiftCertification`] its summary is computed from.
pub(super) struct GammaLiftSummaryParts<'a> {
    pub(super) candidate: &'a GammaLiftCandidate,
    pub(super) input_first_index: usize,
    pub(super) parent_symmetry: &'a GammaParentSymmetryStatus,
    pub(super) identity: &'a GammaLiftIdentityStatus,
    pub(super) parent_proof: &'a VeroneseParentProof,
}

/// The gamma criterion as a transfer theorem (`gamma_lift` module
/// documentation), with the input region `(-inf, 0]` for a positive scale
/// and `[0, inf)` for a negative one.
fn gamma_root_transfer(scale: Option<&BigRational>, first_index: usize) -> LiftRootTransfer {
    let (input_root_region, side) = match scale {
        Some(c) if c > &BigRational::zero() => (nonpositive_ray(), "(-inf, 0]"),
        Some(c) if c < &BigRational::zero() => (
            RootInterval {
                lower: Some(RootEndpoint {
                    value: "0".into(),
                    per_index: None,
                    closed: true,
                }),
                upper: None,
            },
            "[0, inf)",
        ),
        _ => (real_line(), "(-inf, inf)"),
    };
    LiftRootTransfer {
        theorem:
            "gamma criterion (RealRooted IsGammaExpansion.isRealRooted_and_hasRootsNonpos_iff): \
                  if P = sum_j gamma_j x^j (1 + x)^(d - 2j) with deg Gamma <= d/2, then P is \
                  nonzero and real-rooted with roots in (-inf, 0] if and only if \
                  Gamma(z) = sum_j gamma_j z^j is; Gamma(c z) with c != 0 has the roots of \
                  Gamma divided by c"
                .into(),
        parent_hypothesis: "every parent row used is zero or real-rooted with all roots in \
                            (-inf, 0] (coefficients of one sign), and symmetric at its index d_m"
            .into(),
        parent_root_region: Some(nonpositive_ray()),
        input_root_region,
        input_conclusion: format!(
            "every Q_n, n >= {first_index}, is zero or real-rooted with roots in {side}"
        ),
    }
}

/// The [`LiftSummary`] of a gamma lift: identity and parent proof map one to
/// one; the hypotheses are the parent symmetry (decided by the reversal
/// induction) and `lambda` finite and nonzero (decided exactly).
pub(super) fn gamma_lift_summary(parts: GammaLiftSummaryParts<'_>) -> LiftSummary {
    let candidate = parts.candidate;
    let n0 = parts.input_first_index;
    let s = candidate.index_shift;
    let identity = match parts.identity {
        GammaLiftIdentityStatus::Certified { .. } => LiftIdentityStatus::Certified {
            method: "forward_division".into(),
        },
        GammaLiftIdentityStatus::SupportedByRows { rows, .. } => {
            LiftIdentityStatus::SupportedByRows { rows: *rows }
        }
        GammaLiftIdentityStatus::Refuted { index, reason } => LiftIdentityStatus::Refuted {
            index: Some(*index),
            reason: reason.clone(),
        },
    };
    let symmetry_name = "parent_rows_symmetric";
    let symmetry_statement = "every parent row P_m, m >= p0, is symmetric at d_m";
    let symmetry = match parts.parent_symmetry {
        GammaParentSymmetryStatus::Certified { detail } => LiftHypothesis {
            name: symmetry_name.into(),
            statement: symmetry_statement.into(),
            status: LiftHypothesisStatus::Verified,
            evidence: detail.clone(),
            sign_decision: None,
        },
        GammaParentSymmetryStatus::SupportedByRows { rows, reason } => LiftHypothesis {
            name: symmetry_name.into(),
            statement: symmetry_statement.into(),
            status: LiftHypothesisStatus::Unresolved,
            evidence: format!("the compared rows (from {rows} input rows) are symmetric; {reason}"),
            sign_decision: None,
        },
        GammaParentSymmetryStatus::Refuted { index } => LiftHypothesis {
            name: symmetry_name.into(),
            statement: symmetry_statement.into(),
            status: LiftHypothesisStatus::Failed {
                witness: Some(BigInt::from(*index)),
            },
            evidence: format!("P_{index} is not symmetric at d_{index}"),
            sign_decision: None,
        },
    };
    let scale = parse_rational_coeff(&candidate.scale)
        .ok()
        .filter(|c| !c.is_zero());
    LiftSummary::new(
        LiftOperation::GammaPolynomial {
            symmetry_schedule: candidate.parent.symmetry_schedule.clone(),
            scale: candidate.scale.clone(),
        },
        n0,
        LiftParent {
            recurrence: candidate.parent.recurrence.clone(),
            relation: format!(
                "Q_n = lambda(n) Gamma_(n{s:+})(c z) with c = {} and \
                 P_m = sum_j Gamma_(m,j) x^j (1 + x)^(d_m - 2j), for every n >= {n0}",
                candidate.scale
            ),
            gauge: Some(candidate.gauge.clone()),
            index_shift: s,
        },
        identity,
        gamma_root_transfer(scale.as_ref(), n0),
        vec![symmetry, gauge_hypothesis(&candidate.gauge, n0)],
        parent_proof_summary(parts.parent_proof),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hypothesis(status: LiftHypothesisStatus) -> LiftHypothesis {
        LiftHypothesis {
            name: "h".into(),
            statement: "h".into(),
            status,
            evidence: String::new(),
            sign_decision: None,
        }
    }

    #[test]
    fn the_verdict_needs_every_part() {
        let certified = LiftIdentityStatus::Certified { method: "m".into() };
        let proved = LiftParentProof::Proved {
            route: "r".into(),
            statement: String::new(),
        };
        let verified = [hypothesis(LiftHypothesisStatus::Verified)];
        assert_eq!(
            lift_verdict(&certified, &proved, &verified),
            LiftVerdict::Proved
        );
        assert_eq!(lift_verdict(&certified, &proved, &[]), LiftVerdict::Proved);
        let unresolved = [hypothesis(LiftHypothesisStatus::Unresolved)];
        assert_eq!(
            lift_verdict(&certified, &proved, &unresolved),
            LiftVerdict::Open
        );
        let failed = [hypothesis(LiftHypothesisStatus::Failed { witness: None })];
        assert_eq!(
            lift_verdict(&certified, &proved, &failed),
            LiftVerdict::RouteFailed
        );
        let rows = LiftIdentityStatus::SupportedByRows { rows: 9 };
        assert_eq!(lift_verdict(&rows, &proved, &verified), LiftVerdict::Open);
        let attested = LiftParentProof::Attested {
            labels: vec!["a".into()],
            statement: String::new(),
        };
        assert_eq!(
            lift_verdict(&certified, &attested, &verified),
            LiftVerdict::ConditionalOnAttestation
        );
        assert_eq!(lift_verdict(&rows, &attested, &verified), LiftVerdict::Open);
        let failing = LiftParentProof::Failed {
            index: Some(3),
            statement: String::new(),
        };
        assert_eq!(
            lift_verdict(&certified, &failing, &verified),
            LiftVerdict::RouteFailed
        );
    }

    #[test]
    fn gauge_hypotheses_are_decided_exactly() {
        let strings = |values: &[&str]| values.iter().map(|s| s.to_string()).collect::<Vec<_>>();
        // lambda(n) = 1 / (n (n - 3)) from n0 = 4: verified.
        let closed = VeroneseSectionGauge::Closed {
            numerator: strings(&["1"]),
            denominator: strings(&["0", "-3", "1"]),
        };
        let verified = gauge_hypothesis(&closed, 4);
        assert_eq!(verified.status, LiftHypothesisStatus::Verified);
        assert!(verified.sign_decision.as_ref().is_some_and(|d| d.holds()));
        // From n0 = 2 the denominator vanishes at n = 3.
        assert_eq!(
            gauge_hypothesis(&closed, 2).status,
            LiftHypothesisStatus::Failed {
                witness: Some(BigInt::from(3))
            }
        );
        // lambda(n) / lambda(n - 1) = (2 n - 5) / 1 has no integer root.
        let ratio = VeroneseSectionGauge::Ratio {
            numerator: strings(&["-5", "2"]),
            denominator: strings(&["1"]),
            initial_value: "1".into(),
        };
        assert_eq!(
            gauge_hypothesis(&ratio, 0).status,
            LiftHypothesisStatus::Verified
        );
        let json = serde_json::to_value(gauge_hypothesis(&closed, 2)).unwrap();
        assert_eq!(json["status"]["kind"], "failed");
        assert_eq!(json["status"]["witness"], "3");
        assert_eq!(json["sign_decision"]["decision"], "fails_at");
    }
}
