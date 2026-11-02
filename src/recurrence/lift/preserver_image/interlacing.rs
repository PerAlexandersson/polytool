//! The interlacing claim of a preserver-image lift: whether consecutive
//! input rows `Q_(n-1)`, `Q_n` interlace, transferred from the parent's
//! proved interlacing by the Hermite--Kakeya--Obreschkoff argument.
//!
//! # The rule
//!
//! Interlacing is claimed only when all of the following hold; otherwise the
//! claim is `none` with the reason (never from finitely many rows):
//!
//! 1. the identity `Q_n = lambda(n) T(P_(n+s))` is certified for every
//!    `n >= n0`;
//! 2. the operator `T` does not depend on `n` and maps *every* real-rooted
//!    polynomial (of any degree) to a real-rooted polynomial or zero: a fixed
//!    `f(theta)` with real roots `<= 0` (Laguerre), a fixed real-rooted
//!    `phi(D)`, `D^j o x^m`, `x^m D^j`, an affine substitution, a reversal or
//!    Möbius map with a constant schedule `d_n = d`, or a hypergeometric
//!    multiplier free of `n` that is a multiplier sequence (`1/k!`,
//!    `1/(alpha)_k`); its side conditions are verified;
//! 3. a parent theorem report (the three-term or Liu--Wang checks, replayed)
//!    proves every parent row from its first index on real-rooted and
//!    `P_(m-1)` interlacing `P_m` for every `m >= m0` (its `interlacing_from`).
//!
//! **Argument.**  For `m >= m0` the rows `P_(m-1)`, `P_m` are real-rooted and
//! interlace, so by the Hermite--Kakeya--Obreschkoff theorem (Obreschkoff,
//! *Verteilung und Berechnung der Nullstellen reeller Polynome*, 1963; see
//! Fisk, *Polynomials, Roots, and Interlacing*) every combination
//! `alpha P_(m-1) + beta P_m`, `alpha, beta` real, is real-rooted or zero.
//! `T` is linear and preserves real-rootedness of every polynomial, so
//! `alpha T(P_(m-1)) + beta T(P_m) = T(alpha P_(m-1) + beta P_m)` is
//! real-rooted or zero for all `alpha, beta`; by the converse direction of
//! the same theorem `T(P_(m-1))` and `T(P_m)` have weakly interlacing zeros
//! (their degrees then differ by at most one; a zero row interlaces every
//! polynomial by convention).  The nonzero factors `lambda(n)` do not move
//! zeros.  With `m = n + s` this holds for every `n` with
//! `n - 1 >= n0` and `n + s >= m0`.
//!
//! **Kind and orientation.**  The converse gives weak interlacing and no
//! direction: the claim is `weak` with the orientation `unordered` (the
//! zeros of `Q_(n-1)` and `Q_n` alternate weakly; which one carries the
//! largest zero is not determined).  An index-dependent `T_n` gets no
//! claim: that would need a proved relation such as
//! `T_(n+1) = U_n o T_n` with `U_n` preserving interlacing, which is not
//! implemented.

use super::certify::{PreserverImageCandidate, PreserverImageIdentityStatus};
use super::operator::PreserverOperator;
use super::transfer::operator_transfer;
use crate::recurrence::{
    HadamardFactorTerm, LiftHypothesisStatus, RealRootednessTheoremReport, RowInterlacing,
};
use serde::{Deserialize, Serialize};

/// Whether consecutive input rows are claimed to interlace (module
/// documentation of the interlacing rule).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case", deny_unknown_fields)]
pub enum PreserverInterlacingClaim {
    /// No interlacing is claimed.
    None {
        /// Why the rule does not apply.
        reason: String,
    },
    /// `Q_(n-1)` and `Q_n` interlace for every `n >= from_index`.
    Claimed {
        /// The kind (always [`RowInterlacing::Weak`] from the
        /// Hermite--Kakeya--Obreschkoff argument).
        kind: RowInterlacing,
        /// The orientation: `unordered` (zeros alternate weakly, the
        /// direction is not determined).
        orientation: String,
        /// The first input index `n` with `Q_(n-1)`, `Q_n` interlacing; it
        /// holds for every later index.
        from_index: usize,
        /// The argument, with the parent's interlacing it starts from.
        justification: String,
    },
}

impl PreserverInterlacingClaim {
    fn none(reason: impl Into<String>) -> Self {
        Self::None {
            reason: reason.into(),
        }
    }

    /// Whether interlacing is claimed.
    pub fn is_claimed(&self) -> bool {
        matches!(self, Self::Claimed { .. })
    }
}

/// Why `operator` is not an `n`-independent preserver of every real-rooted
/// polynomial, or `None` when it is (its side conditions verified).
fn independent_preserver(operator: &PreserverOperator, n0: usize) -> Option<String> {
    let constant = |coefficients: &[super::operator::IndexRationalFunction]| {
        coefficients.iter().all(|c| {
            c.parse().is_ok_and(|(p, q)| {
                crate::vec_poly::degree(&p).is_none_or(|d| d == 0)
                    && crate::vec_poly::degree(&q).is_none_or(|d| d == 0)
            })
        })
    };
    let shape_reason = match operator {
        PreserverOperator::EulerPolynomial { coefficients }
        | PreserverOperator::DerivativePolynomial { coefficients } => {
            (!constant(coefficients)).then(|| "the operator coefficients depend on n".to_string())
        }
        PreserverOperator::DerivativeOfPowerMultiple { .. }
        | PreserverOperator::PowerTimesDerivative { .. }
        | PreserverOperator::Affine { .. } => None,
        PreserverOperator::Reversal { schedule } | PreserverOperator::Mobius { schedule, .. } => {
            (schedule.step != 0).then(|| "the schedule d_n depends on n".to_string())
        }
        PreserverOperator::HadamardMultiplier { factor } => factor
            .terms
            .iter()
            .any(|term| {
                !matches!(
                    term,
                    HadamardFactorTerm::InverseFactorial
                        | HadamardFactorTerm::InverseRisingFactorial { .. }
                )
            })
            .then(|| "the factor depends on n".to_string()),
        PreserverOperator::WeylOperator { .. } => {
            Some("a general Weyl operator has no preservation theorem".to_string())
        }
    };
    if shape_reason.is_some() {
        return shape_reason;
    }
    let transfer = operator_transfer(operator, n0);
    if !transfer.from_real_line
        || transfer
            .hypotheses
            .iter()
            .any(|h| h.status != LiftHypothesisStatus::Verified)
    {
        return Some(
            "the operator is not proved to preserve real-rootedness of every polynomial".into(),
        );
    }
    None
}

/// The interlacing claim of a certification (module documentation).
pub(super) fn interlacing_claim(
    candidate: &PreserverImageCandidate,
    input_first_index: usize,
    identity: &PreserverImageIdentityStatus,
    parent_report: Option<&RealRootednessTheoremReport>,
) -> PreserverInterlacingClaim {
    if !matches!(identity, PreserverImageIdentityStatus::Certified { .. }) {
        return PreserverInterlacingClaim::none("the identity is not certified");
    }
    if let Some(reason) = independent_preserver(&candidate.operator, input_first_index) {
        return PreserverInterlacingClaim::none(reason);
    }
    let Some(report) = parent_report.filter(|report| report.is_proved()) else {
        return PreserverInterlacingClaim::none(
            "no parent theorem report with an interlacing conclusion",
        );
    };
    let Some(conclusion) = &report.conclusion else {
        return PreserverInterlacingClaim::none("the parent report has no conclusion");
    };
    let (Some(kind), Some(parent_from)) = (conclusion.interlacing, conclusion.interlacing_from)
    else {
        return PreserverInterlacingClaim::none("the parent proof claims no interlacing");
    };
    let s = candidate.index_shift;
    let n0 = input_first_index as i64;
    if (conclusion.first_index as i64) > n0 + s {
        return PreserverInterlacingClaim::none(
            "the parent theorem does not cover every parent row used",
        );
    }
    let from = (parent_from as i64 - s).max(n0 + 1);
    let Ok(from_index) = usize::try_from(from) else {
        return PreserverInterlacingClaim::none("index overflow");
    };
    let parent_kind = match kind {
        RowInterlacing::Strict => "strictly",
        RowInterlacing::Weak => "weakly",
    };
    PreserverInterlacingClaim::Claimed {
        kind: RowInterlacing::Weak,
        orientation: "unordered".into(),
        from_index,
        justification: format!(
            "the parent theorem ({:?}) proves P_(m-1) {parent_kind} interlacing P_m for every \
             m >= {parent_from}; T = {} is independent of n and maps every real-rooted \
             polynomial to a real-rooted polynomial or zero, so by the Hermite--Kakeya--\
             Obreschkoff theorem (both directions) T(P_(m-1)) and T(P_m) have weakly \
             interlacing zeros; the nonzero gauge does not move zeros, so Q_(n-1) and Q_n \
             interlace weakly (orientation not determined) for every n >= {from_index}",
            report.theorem,
            candidate.operator.text()
        ),
    }
}
