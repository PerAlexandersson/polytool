//! The report types of [`find_hadamard_factor_lifts`] and their
//! `polytool.recurrence-hadamard-factor-lifts.v1` JSON form, with the
//! mapping to [`LiftSummary`].

use super::super::lift_summary::{nonpositive_ray, parent_proof_summary};
use super::certify::HadamardFactorCertificate;
use super::factor::{HadamardFactor, POLYA_SCHUR_THEOREM};
#[cfg(doc)]
use crate::recurrence::find_hadamard_factor_lifts;
use crate::recurrence::{
    LiftHypothesis, LiftHypothesisStatus, LiftIdentityStatus, LiftOperation, LiftParent,
    LiftRootTransfer, LiftSummary, RecurrenceError, RecurrenceJson, VeroneseParentProof,
    VeroneseSectionGauge,
};
use num_bigint::BigInt;
use serde::{Deserialize, Serialize};

/// JSON schema tag of a serialized [`HadamardFactorLiftReport`].
pub const HADAMARD_FACTOR_LIFTS_SCHEMA: &str = "polytool.recurrence-hadamard-factor-lifts.v1";

/// How the quotient coefficients `p(n, k)` were chosen where `u(n, k) = 0`
/// (there `q(n, k) = 0` and `p` is not determined by the input).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum HadamardHiddenCoefficients {
    /// Set to zero.
    Zero,
    /// Continued from the visible coefficients by a ratio
    /// `p(n, k + 1) / p(n, k) = (a0 + a1 n + a2 k) / (b0 + b1 n + b2 k)`
    /// fitted exactly to every row (detection only; the certificate does not
    /// depend on these values, which are multiplied by zero).
    Continuation,
}

/// Status of the identity `Q_n = lambda(n) sum_k u(n, k) p(n, k) x^k` for
/// every `n >= n0`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case", deny_unknown_fields)]
pub enum HadamardFactorIdentityStatus {
    /// Proved for every `n >= n0` (module documentation of
    /// `hadamard_factor`).
    Certified {
        /// The replayable certificate.
        certificate: Box<HadamardFactorCertificate>,
    },
    /// The compared rows agree, but the certificate does not apply: finite
    /// evidence only.
    SupportedByRows {
        /// Number of rows compared from `n0`.
        rows: usize,
        /// Why the certificate does not apply.
        reason: String,
    },
    /// `Q_index` differs from the candidate row.
    Refuted {
        /// The absolute index of the first differing row.
        index: usize,
    },
}

/// One Hadamard-factor lift `Q_n = lambda(n) sum_k u(n, k) [x^k] P_n x^k`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HadamardFactorLift {
    /// The factor `u`.
    pub factor: HadamardFactor,
    /// `u` as text.
    pub factor_text: String,
    /// The row normalisation `lambda(n)` (finite and nonzero for `n >= n0`).
    pub gauge: VeroneseSectionGauge,
    /// How the undetermined quotient coefficients were chosen.
    pub hidden_coefficients: HadamardHiddenCoefficients,
    /// The parent `P` (`polytool.recurrence.v1`, first index `n0`).
    pub parent: RecurrenceJson,
    /// Number of rows `Q_(n0), ...` compared exactly with the candidate.
    pub rows_compared: usize,
    /// The identity status.
    pub identity: HadamardFactorIdentityStatus,
    /// Why `k -> u(n, k)` is a nonnegative multiplier sequence for every
    /// `n >= n0`, or `None` when no transfer theorem is recorded for the
    /// factor (the lift is then open whatever the parent).
    pub multiplier_sequence: Option<String>,
    /// The parent statement (every `P_n` zero or real-rooted with
    /// coefficients of one sign) and the provers tried.
    pub parent_proof: VeroneseParentProof,
    /// The uniform [`LiftSummary`] (always equal to
    /// [`HadamardFactorLift::recompute_lift_summary`]).
    pub lift_summary: LiftSummary,
}

/// Why a factor of the menu gave no lift.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum HadamardFactorRejectionKind {
    /// The factor is not defined on the index domain.
    InvalidForDomain,
    /// `u(n, k)` is zero or undefined at a nonzero input coefficient.
    VanishesOnSupport,
    /// No parent recurrence smaller than the input was found for the
    /// quotient rows in the searched space.
    NoSmallerParent,
    /// The parent search ran out of its candidate budget.
    ParentSearchBudgetExhausted,
}

/// A factor that gave no lift (near misses are recorded here).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HadamardFactorRejection {
    /// The factor.
    pub factor: HadamardFactor,
    /// The factor as text.
    pub factor_text: String,
    /// The kind of rejection.
    pub kind: HadamardFactorRejectionKind,
    /// Details (first offending index, sizes compared).
    pub detail: String,
}

/// `polytool.recurrence-hadamard-factor-lifts.v1`: the result of
/// [`find_hadamard_factor_lifts`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HadamardFactorLiftReport {
    /// Always [`HADAMARD_FACTOR_LIFTS_SCHEMA`].
    pub schema: String,
    /// The input's first index `n0`.
    pub input_first_index: usize,
    /// Number of input rows used for detection.
    pub rows_examined: usize,
    /// Number of menu factors tried.
    pub factors_tried: usize,
    /// The lifts found.
    pub lifts: Vec<HadamardFactorLift>,
    /// The factors that gave no lift.
    pub rejected: Vec<HadamardFactorRejection>,
    /// Whether the factor budget stopped the menu early.
    pub exhausted_budget: bool,
    /// Further remarks.
    pub notes: Vec<String>,
}

impl HadamardFactorLift {
    /// The [`LiftSummary`] of this lift, computed from the record.
    pub fn recompute_lift_summary(&self) -> LiftSummary {
        hadamard_lift_summary(
            &self.factor_text,
            &self.gauge,
            &self.parent,
            &self.identity,
            self.multiplier_sequence.as_deref(),
            &self.parent_proof,
        )
    }

    /// Replay the certificate (when certified), check that it belongs to
    /// this lift, and that the recorded summary equals the recomputed one.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        if let HadamardFactorIdentityStatus::Certified { certificate } = &self.identity {
            certificate.verify()?;
            if certificate.factor != self.factor
                || certificate.gauge != self.gauge
                || certificate.parent != self.parent
            {
                return Err(RecurrenceError::VerificationFailed(
                    "the certificate does not belong to the lift".into(),
                ));
            }
        }
        if self.factor_text != self.factor.text()
            || self.lift_summary != self.recompute_lift_summary()
            || !self.lift_summary.is_consistent()
        {
            return Err(RecurrenceError::VerificationFailed(
                "the lift summary does not match the one recomputed from the record".into(),
            ));
        }
        Ok(())
    }
}

impl HadamardFactorLiftReport {
    /// Check the schema and replay every lift ([`HadamardFactorLift::verify`]).
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        if self.schema != HADAMARD_FACTOR_LIFTS_SCHEMA {
            return Err(RecurrenceError::InvalidInput(format!(
                "unsupported schema `{}`",
                self.schema
            )));
        }
        self.lifts.iter().try_for_each(HadamardFactorLift::verify)
    }

    /// The lifts whose summary is proved.
    pub fn proved_lifts(&self) -> impl Iterator<Item = &HadamardFactorLift> {
        self.lifts
            .iter()
            .filter(|lift| lift.lift_summary.is_proved())
    }
}

/// The [`LiftSummary`] of a Hadamard-factor lift from its parts.
pub(super) fn hadamard_lift_summary(
    u: &str,
    gauge: &VeroneseSectionGauge,
    parent: &RecurrenceJson,
    identity: &HadamardFactorIdentityStatus,
    multiplier_sequence: Option<&str>,
    parent_proof: &VeroneseParentProof,
) -> LiftSummary {
    let n0 = parent.first_index;
    let identity = match identity {
        HadamardFactorIdentityStatus::Certified { .. } => LiftIdentityStatus::Certified {
            method: "conjugation_division".into(),
        },
        HadamardFactorIdentityStatus::SupportedByRows { rows, .. } => {
            LiftIdentityStatus::SupportedByRows { rows: *rows }
        }
        HadamardFactorIdentityStatus::Refuted { index } => LiftIdentityStatus::Refuted {
            index: Some(*index),
            reason: format!("Q_{index} differs from lambda(n) u ⊙ P_n at n = {index}"),
        },
    };
    let mut hypotheses = vec![match &multiplier_sequence {
        Some(justification) => LiftHypothesis {
            name: "factor_multiplier_sequence".into(),
            statement: format!("k -> {u} is a nonnegative multiplier sequence for every n >= {n0}"),
            status: LiftHypothesisStatus::Verified,
            evidence: justification.to_string(),
            sign_decision: None,
        },
        None => LiftHypothesis {
            name: "factor_multiplier_sequence".into(),
            statement: format!(
                "k -> {u} preserves real-rootedness with roots in (-inf, 0] for every n >= \
                     {n0}"
            ),
            status: LiftHypothesisStatus::Unresolved,
            evidence: "no transfer theorem is recorded for this factor kind".into(),
            sign_decision: None,
        },
    }];
    if *gauge != VeroneseSectionGauge::None {
        let (status, evidence) = match gauge.ratio_form(n0) {
            Ok(_) => (
                LiftHypothesisStatus::Verified,
                "both ratio polynomials are nonzero at every integer n > n0 (exact sign \
                     checker) and lambda(n0) != 0"
                    .to_string(),
            ),
            Err(RecurrenceError::BudgetExhausted(reason)) => {
                (LiftHypothesisStatus::Unresolved, reason)
            }
            Err(error) => (
                LiftHypothesisStatus::Failed {
                    witness: Some(BigInt::from(n0)),
                },
                error.to_string(),
            ),
        };
        hypotheses.push(LiftHypothesis {
            name: "gauge_finite_nonzero".into(),
            statement: "lambda(n) is finite and nonzero for every n >= n0".into(),
            status,
            evidence,
            sign_decision: None,
        });
    }
    let theorem = if multiplier_sequence.is_some() {
        POLYA_SCHUR_THEOREM.to_string()
    } else {
        format!("no transfer theorem is recorded for the factor {u}")
    };
    LiftSummary::new(
        LiftOperation::HadamardFactor {
            factor: u.to_string(),
        },
        n0,
        LiftParent {
            recurrence: parent.clone(),
            relation: format!("Q_n = lambda(n) sum_k {u} [x^k]P_n x^k for every n >= {n0}"),
            gauge: (*gauge != VeroneseSectionGauge::None).then(|| gauge.clone()),
            index_shift: 0,
        },
        identity,
        LiftRootTransfer {
            theorem,
            parent_hypothesis: "every parent row used is zero or real-rooted with all roots \
                                    in (-inf, 0] (coefficients of one sign)"
                .into(),
            parent_root_region: Some(nonpositive_ray()),
            input_root_region: nonpositive_ray(),
            input_conclusion: format!(
                "every Q_n, n >= {n0}, is zero or real-rooted with roots in (-inf, 0]: a \
                     nonnegative multiplier sequence keeps real-rootedness and one coefficient \
                     sign"
            ),
        },
        hypotheses,
        parent_proof_summary(parent_proof),
    )
}
