//! Degree laws and signs of leading coefficients of recurrence rows, decided
//! for every index.
//!
//! For a family given by a rule
//!
//! ```text
//! f(n) P_n = sum_t s_t(n) c_t(n, x) D^(d_t) P_(n - j_t),      n >= N,
//! ```
//!
//! (`f` free of `x`, after the normalization of a left-hand factor
//! `x^e f(n)` described in the `parent_theorems` documentation; `s_t(n)` is
//! `1` or `(-1)^n`; no inhomogeneous term), a first index and seed rows
//! (`N = first_index + seeds`), [`certify_leading_coefficient`] decides a
//! **degree law** `deg P_n = D(n)` and the **sign of the leading coefficient**
//! `lead P_n` for every `n` from some index on, exactly.  The result is a
//! [`LeadingCoefficientCertificate`]: the law, the route that proves it, the
//! hypotheses (each an exact row computation or a sign statement about a
//! rational function of the index decided by
//! [`crate::decide_index_rational_function_sign`], with its replayable
//! certificate), and the RealRooted declarations of the route.
//! [`LeadingCoefficientCertificate::verify`] replays the recorded route.
//!
//! # Degree laws
//!
//! A [`RecurrenceDegreeLaw`] is periodic by residue class: from its start
//! `s` on, with period `p`, base degrees `b_0, ..., b_(p-1)` and growth `g`
//! per period,
//!
//! ```text
//! D(s + q p + r) = b_r + g q,          sign lead P_(s + q p + r) = sigma_r.
//! ```
//!
//! Affine laws are the case `p = 1`; half growth `D0 + floor((n + e) / 2)`
//! has `p = 2`, `g = 1`; Veronese sections and two-step products give other
//! periods.  The law is read off exactly computed rows (the smallest period,
//! then the earliest start, that fits every sampled row); the rows before the
//! induction start are recorded exactly
//! ([`LeadingCoefficientCertificate::initial_rows`]), and the law is then
//! proved for every later index by one of the routes below.
//!
//! # The top-coefficient recurrence
//!
//! Let `T = max(N, s + K)` be the induction start (`K` the largest lag) and
//! split `n >= T` into the classes `n = c + p m`, `m >= 0`,
//! `T <= c < T + p`.  On a class, `e_j = D(n) - D(n - j)` is constant and
//! `D(n - j) = D(c - j) + g m` is affine in `m`.  Assume every summand obeys
//! the degree bound `deg_x c_t(n, x) <= e_(j_t) + d_t` (the coefficients of
//! higher powers of `x` are the zero polynomial in `n`) and `e_j >= 0`.  If
//! the rows `P_(n-j)` have degree at most `D(n - j)`, the right-hand side has
//! degree at most `D(n)` and its coefficient of `x^(D(n))` is
//!
//! ```text
//! f(n) lambda(n) = sum_j M_j(n) lambda(n - j),
//! M_j(n) = sum_(t: j_t = j) s_t(n) [x^(e_j + d_t)] c_t(n, x) * (D(n - j))_(d_t),
//! ```
//!
//! with `lambda(n) = [x^(D(n))] P_n` and `(D)_d = D (D - 1) ... (D - d + 1)`
//! (RealRooted `LinRec.lagMultD`).  With `mu(n) = sigma(n) lambda(n)` and the
//! signed multipliers `m_j(n) = sigma(n) sigma(n - j) M_j(n) / f(n)` this is
//! `mu(n) = sum_j m_j(n) mu(n - j)`, and `mu > 0` gives both the degree law
//! (`lambda(n) != 0`) and the signs.  Every `m_j` is a rational function of
//! `m` on each class.
//!
//! # Routes
//!
//! The route names are those of the RealRooted tactic `rr_row_leadingCoeff_pos`
//! (<https://github.com/PerAlexandersson/RealRooted>,
//! `RealRooted/Tactic/ROUTES.md`):
//!
//! - `cone`: `m_j(n) >= 0` for every lag and `sum_j m_j(n) > 0`; then
//!   `mu(n) > 0` by induction from positive initial values
//!   (`RealRooted.LinRec.natDegree_eq_and_leadingCoeff_pos` for an affine
//!   law, `..._of_degreeLaw` otherwise);
//! - `recurrence`: multipliers of either sign with a growth-ratio invariant
//!   `mu(n) >= rho mu(n - 1)` for a rational `rho > 0`
//!   (`RealRooted.LinRec.natDegree_eq_and_leadingCoeff_pos_of_ratio`,
//!   `RealRooted.ScalarRec.pos_of_ratio`).  Writing the window
//!   `x_i = mu(n - K + i)`, `y_0 = x_0`, `y_(i+1) = x_(i+1) - rho x_i`, the
//!   step `mu(n) - rho x_(K-1) >= 0` holds for every window with
//!   `y_0 > 0`, `y_i >= 0` exactly when every coefficient
//!   `C_i(n) = sum_(t >= i) m_(K-t)(n) rho^(t-i) - rho^(K-i)`,
//!   `0 <= i < K`, is `>= 0`; these are the decided conditions, and the
//!   initial window must satisfy `mu(T - K) > 0`,
//!   `rho mu(t) <= mu(t + 1)`;
//! - `subseq`: the `cone` conditions for a periodic law
//!   whose lags are all multiples of the period, so that the residue classes
//!   evolve separately (`RealRooted.LinRec.degreeLaw_periodic_decomp`);
//! - `twoStep`: `P_n = q(x) P_(n-2)` with a fixed `q` of positive leading
//!   coefficient (`RealRooted.twoStepProduct_natDegree`,
//!   `RealRooted.twoStepProduct_leadingCoeff_pos`);
//! - `cancel`: first-order rules `P_n = A_n(x) P'_(n-1) + B_n(x) P_(n-1)`
//!   with `deg_x A <= m + 1`, `deg_x B <= m` (`m >= 1`) whose top terms either
//!   grow (`A_(m+1) D + B_m > 0`, step `m`) or cancel exactly
//!   (`A_(m+1) D + B_m = 0`, `A_(m+1) <= 0`, `A_m D + B_(m-1) > 0`, step
//!   `m - 1`), with `D = D(n - 1)`, while the coefficients stay nonnegative:
//!   `A_0 >= 0` and `B_i + r A_(i+1) >= 0` for `0 <= r <= D(n - 1)` (checked
//!   at both ends, the condition being affine in `r`)
//!   (`RealRooted.derivativeRecurrence_natDegree_eq_and_leadingCoeff_pos_of_cancel`).
//!   On cancellation the coefficient of `x^(D + m - 1)` is
//!   `-A_(m+1) [x^(D-1)] P_(n-1) + (A_m D + B_(m-1)) lead P_(n-1) > 0`.
//!   When the next level alone is not positive, the route carries the
//!   coefficient-ratio invariant `[x^(D-1)] P_n >= rho lead P_n` for a
//!   rational `rho > 0` (RealRooted
//!   `LinRec.natDegree_eq_and_leadingCoeff_pos_of_cancel`, as in the tactic
//!   `rr_row_cancel`): with `E_i(r) = B_i + r A_(i+1)` (`B_(-1) = 0`), the
//!   new top two coefficients are linear forms `alpha p_D + beta p_(D-1)`
//!   in the top two of `P_(n-1)` (a third term `-2 A_(m+1) p_(D-2) >= 0` on
//!   cancellation is dropped), and such a form is `>= 0` on the cone
//!   `p_D > 0`, `p_(D-1) >= rho p_D` exactly when `beta >= 0` and
//!   `alpha + rho beta >= 0`.  The decided conditions are
//!   `A_m D + B_(m-1) - rho A_(m+1) > 0` (positive top coefficient on
//!   cancellation) and these two for `[x^(D'-1)] P_n - rho lead P_n`, with
//!   `D(n - 1) >= 1` and the invariant on the row before the induction
//!   start;
//! - `coefficients`: finite-lag rules with at most first derivatives,
//!   `P_n = sum_j (B_j P_(n-j) + A_j P'_(n-j))`, under an affine degree
//!   envelope `D(n)` (`deg_x B_j <= e_j`, `deg_x A_j <= e_j + 1`):
//!   `[x^0] A_j >= 0` and `B_(j,i) + r A_(j,i+1) >= 0` for
//!   `0 <= r <= D(n - j)` keep every coefficient nonnegative
//!   (`RealRooted.LinRec.natDegree_le_and_hasNonnegCoeffs`); with
//!   `sum_j B_(j,0) > 0` and positive constant terms of the initial window,
//!   `P_n(0) >= sum_j B_(j,0) P_(n-j)(0) > 0`, so every row is nonzero and
//!   its leading coefficient is positive
//!   (`RealRooted.HasNonnegCoeffs.pos_leadingCoeff`).  The degree is only
//!   bounded by the envelope ([`LeadingCoefficientCertificate::degree_is_exact`]
//!   is `false`).
//!
//! Leading coefficients of negative sign are the positive case for
//! `sigma(n) P_n` (route `scale`, `RealRooted.hasPosLeadingCoeff_C_mul`).
//! The route `algebra` (sums, products, powers and derivatives of rows, a
//! statement about given polynomial expressions rather than a recurrence)
//! and the `cancel` route for rules with lags above one (the general case
//! of RealRooted `LinRec.natDegree_eq_and_leadingCoeff_pos_of_cancel`) are
//! not implemented.
//!
//! A certificate that is not proved claims nothing: a failed hypothesis of
//! a route is a failed sufficient condition, not a disproof, except that a
//! zero row among the computed rows is recorded exactly and moves
//! [`LeadingCoefficientCertificate::holds_from`] past it.
//!
//! Related: the theorem checks of `parent_theorems` use these certificates
//! where their own degree steps are too narrow; `polytool
//! recurrence-degree` is the command-line query.

use super::lift::rule_with_x_free_left_hand_factor;
use super::{
    format_rational_coeff, LiftHypothesis, LiftHypothesisStatus, Recurrence, RecurrenceError,
    RecurrenceJson, RecurrenceSign,
};
use crate::integer_range_sign::{
    decide_index_rational_function_sign, IndexSignCondition, IndexSignDecision, IntegerRange,
};
use crate::vec_poly;
use crate::BigRational;
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[cfg(test)]
mod tests;

type Q = BigRational;
/// Coefficient matrix `[i][j]` = coefficient of `n^i x^j`.
type Bivariate = Vec<Vec<Q>>;

/// JSON schema tag of [`LeadingCoefficientCertificate`].
pub const LEADING_COEFFICIENT_SCHEMA: &str = "polytool.recurrence-leading-coefficient.v1";

/// Default of [`LeadingCoefficientOptions::sample_rows`].
pub const DEFAULT_LEADING_COEFFICIENT_SAMPLE_ROWS: usize = 36;

/// Default of [`LeadingCoefficientOptions::maximum_period`].
pub const DEFAULT_LEADING_COEFFICIENT_MAXIMUM_PERIOD: usize = 6;

/// Largest degree or growth of a recorded law that
/// [`LeadingCoefficientCertificate::verify`] replays.
const MAXIMUM_RECORDED_DEGREE: usize = 1 << 20;

/// Largest number of law starts tried per period after the earliest one.
const EXTRA_LAW_STARTS: usize = 2;

/// Options of [`certify_leading_coefficient`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LeadingCoefficientOptions {
    /// Number of rows (from the first index, seeds included) computed to
    /// read off the degree law; more are computed when the seeds and the
    /// largest period need them.  Default
    /// [`DEFAULT_LEADING_COEFFICIENT_SAMPLE_ROWS`].
    pub sample_rows: usize,
    /// Largest period of a degree law tried.  Default
    /// [`DEFAULT_LEADING_COEFFICIENT_MAXIMUM_PERIOD`].
    pub maximum_period: usize,
}

impl Default for LeadingCoefficientOptions {
    fn default() -> Self {
        Self {
            sample_rows: DEFAULT_LEADING_COEFFICIENT_SAMPLE_ROWS,
            maximum_period: DEFAULT_LEADING_COEFFICIENT_MAXIMUM_PERIOD,
        }
    }
}

/// The sign of a nonzero leading coefficient.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LeadingCoefficientSign {
    /// `> 0`.
    Positive,
    /// `< 0`.
    Negative,
}

impl LeadingCoefficientSign {
    fn of(value: &Q) -> Option<Self> {
        if value.is_positive() {
            Some(Self::Positive)
        } else if value.is_negative() {
            Some(Self::Negative)
        } else {
            None
        }
    }

    fn factor(self) -> Q {
        match self {
            Self::Positive => Q::one(),
            Self::Negative => -Q::one(),
        }
    }

    fn text(self) -> &'static str {
        match self {
            Self::Positive => "positive",
            Self::Negative => "negative",
        }
    }
}

/// A degree law periodic by residue class (module documentation): for
/// `n = start + q period + r`, `0 <= r < period`,
/// `deg P_n = base_degrees[r] + growth_per_period q`, and the leading
/// coefficient has the sign `leading_coefficient_signs[r]`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceDegreeLaw {
    /// The first index the law describes.
    pub start: usize,
    /// The period `p >= 1` of the residue classes.
    pub period: usize,
    /// `D(start), ..., D(start + period - 1)`.
    pub base_degrees: Vec<usize>,
    /// `D(n + period) - D(n)`.
    pub growth_per_period: usize,
    /// The sign of `lead P_(start + r)` for each residue `r`.
    pub leading_coefficient_signs: Vec<LeadingCoefficientSign>,
}

impl RecurrenceDegreeLaw {
    fn is_well_formed(&self) -> bool {
        self.period >= 1
            && self.base_degrees.len() == self.period
            && self.leading_coefficient_signs.len() == self.period
    }

    /// `D(n)`, or `None` before the start.
    pub fn degree_at(&self, n: usize) -> Option<usize> {
        let offset = n.checked_sub(self.start)?;
        let (q, r) = (offset / self.period, offset % self.period);
        Some(self.base_degrees[r] + self.growth_per_period * q)
    }

    /// The sign of `lead P_n`, or `None` before the start.
    pub fn sign_at(&self, n: usize) -> Option<LeadingCoefficientSign> {
        let offset = n.checked_sub(self.start)?;
        Some(self.leading_coefficient_signs[offset % self.period])
    }

    /// Whether `D(n + 1) - D(n)` is the same for every `n >= start` (an
    /// affine law), and that step.
    pub fn affine_step(&self) -> Option<usize> {
        let steps: BTreeSet<i64> = (self.start + 1..=self.start + self.period)
            .map(|n| self.degree_at(n).unwrap() as i64 - self.degree_at(n - 1).unwrap() as i64)
            .collect();
        match steps.into_iter().collect::<Vec<_>>().as_slice() {
            [step] if *step >= 0 => Some(*step as usize),
            _ => None,
        }
    }

    /// Whether every leading coefficient from the start on has `sign`.
    pub fn has_constant_sign(&self, sign: LeadingCoefficientSign) -> bool {
        self.leading_coefficient_signs.iter().all(|s| *s == sign)
    }

    /// The law in words.
    pub fn text(&self) -> String {
        let signs: Vec<&str> = self
            .leading_coefficient_signs
            .iter()
            .map(|s| s.text())
            .collect();
        if self.period == 1 {
            return format!(
                "deg P_n = {} + {} (n - {}) with {} leading coefficient for n >= {}",
                self.base_degrees[0], self.growth_per_period, self.start, signs[0], self.start
            );
        }
        format!(
            "deg P_(s + {p} q + r) = b_r + {g} q with s = {s}, (b_0, ..., b_{last}) = ({b}) and \
             leading-coefficient signs ({signs}) by residue r",
            p = self.period,
            g = self.growth_per_period,
            s = self.start,
            last = self.period - 1,
            b = self
                .base_degrees
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(", "),
            signs = signs.join(", "),
        )
    }
}

/// The route of a [`LeadingCoefficientCertificate`] (module documentation),
/// named as in the RealRooted tactic `rr_row_leadingCoeff_pos`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum LeadingCoefficientRoute {
    /// Nonnegative signed multipliers with a positive sum.
    #[serde(rename = "cone")]
    Cone,
    /// A growth-ratio invariant `mu(n) >= rho mu(n - 1)`.
    #[serde(rename = "recurrence")]
    Recurrence,
    /// `cone` on residue classes that evolve separately.
    #[serde(rename = "subseq")]
    Subseq,
    /// `P_n = q(x) P_(n-2)` with a fixed `q`.
    #[serde(rename = "twoStep")]
    TwoStep,
    /// First-order rules whose top terms grow or cancel exactly.
    #[serde(rename = "cancel")]
    Cancel,
    /// Nonnegative coefficients within an affine degree envelope.
    #[serde(rename = "coefficients")]
    Coefficients,
}

impl LeadingCoefficientRoute {
    /// The RealRooted route name, for example `"twoStep"`.
    pub fn lean_route_name(self) -> &'static str {
        match self {
            Self::Cone => "cone",
            Self::Recurrence => "recurrence",
            Self::Subseq => "subseq",
            Self::TwoStep => "twoStep",
            Self::Cancel => "cancel",
            Self::Coefficients => "coefficients",
        }
    }
}

/// Outcome of [`certify_leading_coefficient`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LeadingCoefficientStatus {
    /// Every hypothesis of the route is verified: the law holds.
    Proved,
    /// No route closes (a failed sufficient condition, not a disproof).
    NotProved,
    /// The rule is outside the supported form, or its rows cannot be
    /// generated.
    NotApplicable,
}

/// An exactly computed row before the induction start.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LeadingCoefficientRow {
    /// The index `n`.
    pub index: usize,
    /// `deg P_n`, `None` for the zero row.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub degree: Option<usize>,
    /// The sign of `lead P_n`, `None` for the zero row.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sign: Option<LeadingCoefficientSign>,
}

/// The result of [`certify_leading_coefficient`] (schema
/// `polytool.recurrence-leading-coefficient.v1`).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LeadingCoefficientCertificate {
    /// [`LEADING_COEFFICIENT_SCHEMA`].
    pub schema: String,
    /// The family: rule, first index and seed rows (as supplied).
    pub recurrence: RecurrenceJson,
    /// The options of the search.
    pub options: LeadingCoefficientOptions,
    /// The outcome.
    pub status: LeadingCoefficientStatus,
    /// The rule in the normal form used, or why it is not applicable.
    pub normal_form: String,
    /// The route, present for a proved certificate and for the recorded
    /// failed attempt of a not-proved one.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub route: Option<LeadingCoefficientRoute>,
    /// The degree law of the route.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub degree_law: Option<RecurrenceDegreeLaw>,
    /// `true`: `deg P_n = D(n)`; `false` (route `coefficients`):
    /// `deg P_n <= D(n)` only.
    pub degree_is_exact: bool,
    /// The ratio `rho` of the `recurrence` route (`mu(n) >= rho mu(n - 1)`)
    /// or of the `cancel` route with the coefficient-ratio invariant
    /// (`[x^(D-1)] P_n >= rho lead P_n`); exact string in JSON.
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        with = "crate::integer_range_sign::exact_serde::optional_rational"
    )]
    pub ratio: Option<BigRational>,
    /// The induction start `T = max(N, law start + largest lag)`; the rows
    /// before it are [`Self::initial_rows`].
    pub induction_start: usize,
    /// The rows `P_first, ..., P_(T-1)`, computed exactly.
    pub initial_rows: Vec<LeadingCoefficientRow>,
    /// The hypotheses of the route in the order checked.
    pub hypotheses: Vec<LiftHypothesis>,
    /// For a proved certificate: every `P_n`, `n >= holds_from`, is nonzero
    /// with the degree and sign of the law (from the law start on) or of
    /// [`Self::initial_rows`] (before it).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub holds_from: Option<usize>,
    /// The RealRooted route name of [`Self::route`].
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lean_route: Option<String>,
    /// Fully qualified RealRooted declarations of the route.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub lean_declarations: Vec<String>,
    /// The conclusion in words, present exactly when proved.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub conclusion: Option<String>,
    /// For a proved exact law with one sign: the constant `c` with
    /// `lead P_n = c` for every `n >= T - K` (for example `1` for monic
    /// rows), when the initial window has leading coefficient `c` and the
    /// top-coefficient recurrence holds and fixes it (the degree bound on
    /// every class and `sum_j M_j(n) = f(n)`, the hypotheses
    /// [`Self::exact_leading_coefficient_hypotheses`]; never for a class
    /// on which the top terms of the `cancel` route cancel; RealRooted
    /// `LinRec.natDegree_eq_and_leadingCoeff_eq_of_degreeLaw`).  Exact
    /// string in JSON.
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        with = "crate::integer_range_sign::exact_serde::optional_rational"
    )]
    pub exact_leading_coefficient: Option<BigRational>,
    /// The hypotheses of [`Self::exact_leading_coefficient`] (empty when it
    /// is absent).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub exact_leading_coefficient_hypotheses: Vec<LiftHypothesis>,
    /// For a certificate that is not proved: one line per route attempt.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub attempts: Vec<String>,
}

impl LeadingCoefficientCertificate {
    /// Whether the law is proved.
    pub fn is_proved(&self) -> bool {
        self.status == LeadingCoefficientStatus::Proved
    }

    /// For a proved certificate with an exact law: `deg P_n`, for every
    /// `n >= holds_from`.
    pub fn degree_at(&self, n: usize) -> Option<usize> {
        if !self.is_proved() || n < self.holds_from? {
            return None;
        }
        let law = self.degree_law.as_ref()?;
        if !self.degree_is_exact {
            return None;
        }
        match law.degree_at(n) {
            Some(degree) => Some(degree),
            None => self.initial_row(n)?.degree,
        }
    }

    /// For a proved certificate: the sign of `lead P_n`, `n >= holds_from`.
    pub fn sign_at(&self, n: usize) -> Option<LeadingCoefficientSign> {
        if !self.is_proved() || n < self.holds_from? {
            return None;
        }
        match self.degree_law.as_ref()?.sign_at(n) {
            Some(sign) => Some(sign),
            None => self.initial_row(n)?.sign,
        }
    }

    fn initial_row(&self, n: usize) -> Option<&LeadingCoefficientRow> {
        self.initial_rows.iter().find(|row| row.index == n)
    }

    /// For a proved certificate with an exact law: the set of degree steps
    /// `deg P_n - deg P_(n-1)` over all `n > from` (`from >= holds_from`),
    /// and whether every leading coefficient `lead P_n`, `n >= from`, has
    /// the sign of `lead P_from`.
    #[deprecated(since = "0.6.0", note = "unused; it will be removed")]
    pub fn degree_steps_from(&self, from: usize) -> Option<(BTreeSet<i64>, bool)> {
        let law = self.degree_law.as_ref()?;
        let sign = self.sign_at(from)?;
        let mut steps = BTreeSet::new();
        let mut one_sign = true;
        // Past `max(from, law.start) + period` the steps and signs repeat.
        let last = from.max(law.start) + law.period + 1;
        for n in from + 1..=last {
            steps.insert(self.degree_at(n)? as i64 - self.degree_at(n - 1)? as i64);
            one_sign &= self.sign_at(n)? == sign;
        }
        Some((steps, one_sign))
    }

    /// Replay the certificate: recompute the recorded route (or, without
    /// one, the search) from the recorded family and options, compare, and
    /// replay every sign certificate independently.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::VerificationFailed`] when the record differs from
    /// its recomputation or a sign certificate does not replay; errors of
    /// the record itself are passed through.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        if self.schema != LEADING_COEFFICIENT_SCHEMA {
            return Err(RecurrenceError::VerificationFailed(format!(
                "unsupported schema `{}`",
                self.schema
            )));
        }
        let (rule, first_index, seeds) = self.recurrence.to_recurrence_parts()?;
        let recomputed = match (self.route, &self.degree_law) {
            (Some(route), Some(law)) => {
                if !law.is_well_formed() {
                    return Err(RecurrenceError::VerificationFailed(
                        "the recorded degree law is malformed".into(),
                    ));
                }
                // A searched law starts within the sampled rows and has a
                // period of at most the option (the envelope law has period
                // one); bounding the record keeps the replay finite and free
                // of overflow.
                let max_lag = rule.terms.iter().map(|t| t.offset).max().unwrap_or(1);
                let sampled = first_index
                    + seeds.len()
                    + self.options.sample_rows
                    + 3 * self.options.maximum_period
                    + 2 * max_lag
                    + 4;
                if law.start > sampled
                    || law.period > self.options.maximum_period.max(1)
                    || law.growth_per_period > MAXIMUM_RECORDED_DEGREE
                    || law
                        .base_degrees
                        .iter()
                        .any(|&d| d > MAXIMUM_RECORDED_DEGREE)
                {
                    return Err(RecurrenceError::VerificationFailed(
                        "the recorded degree law lies outside the searched range".into(),
                    ));
                }
                let mut certificate = certify_leading_coefficient_route(
                    &rule,
                    first_index,
                    &seeds,
                    &self.options,
                    law,
                    route,
                    self.ratio.as_ref(),
                )?;
                if !self.is_proved() {
                    certificate.status = self.status;
                    certificate.conclusion = None;
                    certificate.holds_from = None;
                    certificate.attempts = self.attempts.clone();
                }
                certificate
            }
            _ => certify_leading_coefficient(&rule, first_index, &seeds, &self.options)?,
        };
        if &recomputed != self {
            return Err(RecurrenceError::VerificationFailed(
                "the certificate differs from its recomputation".into(),
            ));
        }
        for hypothesis in self
            .hypotheses
            .iter()
            .chain(&self.exact_leading_coefficient_hypotheses)
        {
            if let Some(IndexSignDecision::Holds { certificate }) = &hypothesis.sign_decision {
                certificate.verify().map_err(|error| {
                    RecurrenceError::VerificationFailed(format!(
                        "hypothesis {}: {error}",
                        hypothesis.name
                    ))
                })?;
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Entry points
// ---------------------------------------------------------------------------

/// Decide a degree law and the signs of the leading coefficients of every
/// row of the family `rule`, `first_index`, `seeds` (module documentation):
/// the first proved route over the candidate laws, else a not-proved
/// certificate recording the first attempt and a line per attempt.
///
/// # Errors
///
/// [`RecurrenceError::InvalidInput`] when the family cannot be recorded
/// (as [`RecurrenceJson::to_recurrence_parts`]).  An unsupported rule, or
/// rows that cannot be generated, give a
/// [`LeadingCoefficientStatus::NotApplicable`] certificate.
///
/// # Example
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::recurrence::{
///     certify_leading_coefficient, BigRational, BivarPoly, LeadingCoefficientOptions,
///     LeadingCoefficientRoute, Recurrence, RecurrenceTerm,
/// };
///
/// let q = |v: i64| BigRational::from_integer(BigInt::from(v));
/// // P_n = P_(n-1) + x P_(n-2), P_0 = P_1 = 1: deg P_n = floor(n / 2).
/// let rule = Recurrence {
///     terms: vec![
///         RecurrenceTerm::new(1, 0, BivarPoly::from_coefficients_normalized(vec![vec![q(1)]])),
///         RecurrenceTerm::new(
///             2,
///             0,
///             BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(1)]]),
///         ),
///     ],
///     denominator: None,
///     inhomogeneous: None,
/// };
/// let certificate = certify_leading_coefficient(
///     &rule,
///     0,
///     &[vec![q(1)], vec![q(1)]],
///     &LeadingCoefficientOptions::default(),
/// )
/// .unwrap();
/// assert!(certificate.is_proved());
/// assert_eq!(certificate.route, Some(LeadingCoefficientRoute::Cone));
/// assert_eq!(certificate.degree_at(9), Some(4));
/// certificate.verify().unwrap();
/// ```
pub fn certify_leading_coefficient(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<BigRational>],
    options: &LeadingCoefficientOptions,
) -> Result<LeadingCoefficientCertificate, RecurrenceError> {
    let context = match Context::new(rule, first_index, seeds, options)? {
        Ok(context) => context,
        Err(not_applicable) => return Ok(*not_applicable),
    };
    let sample = context.sample_count();
    let rows = match context.rows(sample) {
        Ok(rows) => rows,
        Err(reason) => return Ok(context.not_applicable(reason)),
    };
    let laws = infer_degree_laws(
        &rows,
        context.first_index,
        options.maximum_period.max(1),
        context.shape.max_lag,
        context.shape.alternating,
    );
    let mut attempts = Vec::new();
    let mut first_attempt: Option<LeadingCoefficientCertificate> = None;
    let mut record =
        |certificate: LeadingCoefficientCertificate,
         first_attempt: &mut Option<LeadingCoefficientCertificate>| {
            attempts.push(attempt_line(&certificate));
            if first_attempt.is_none() {
                *first_attempt = Some(certificate);
            }
        };
    for law in &laws {
        for (route, ratio) in context.route_candidates(law, &rows) {
            let certificate = context.certify(law, route, ratio.as_ref(), &rows)?;
            if certificate.is_proved() {
                return Ok(certificate);
            }
            record(certificate, &mut first_attempt);
        }
    }
    // Nonnegative coefficients within an affine envelope.
    if let Some(law) = context.envelope_law(&rows) {
        let certificate =
            context.certify(&law, LeadingCoefficientRoute::Coefficients, None, &rows)?;
        if certificate.is_proved() {
            return Ok(certificate);
        }
        record(certificate, &mut first_attempt);
    }
    let Some(mut certificate) = first_attempt else {
        let reason = if laws.is_empty() {
            format!(
                "no degree law of period at most {} fits the {} computed rows",
                options.maximum_period,
                rows.len()
            )
        } else {
            "no route applies".to_string()
        };
        let mut certificate = context.not_applicable(reason);
        certificate.status = LeadingCoefficientStatus::NotProved;
        return Ok(certificate);
    };
    certificate.status = LeadingCoefficientStatus::NotProved;
    certificate.attempts = attempts;
    Ok(certificate)
}

/// The certificate of one route for a given degree law (and ratio for the
/// `recurrence` route, optionally for the `cancel` route), without a search; this is what
/// [`LeadingCoefficientCertificate::verify`] replays.
///
/// # Errors
///
/// As [`certify_leading_coefficient`]; [`RecurrenceError::InvalidInput`]
/// for a malformed law or a missing ratio of the `recurrence` route.
pub fn certify_leading_coefficient_route(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<BigRational>],
    options: &LeadingCoefficientOptions,
    law: &RecurrenceDegreeLaw,
    route: LeadingCoefficientRoute,
    ratio: Option<&BigRational>,
) -> Result<LeadingCoefficientCertificate, RecurrenceError> {
    if !law.is_well_formed() {
        return Err(RecurrenceError::InvalidInput(
            "the degree law needs a positive period and one degree and sign per residue".into(),
        ));
    }
    let context = match Context::new(rule, first_index, seeds, options)? {
        Ok(context) => context,
        Err(not_applicable) => return Ok(*not_applicable),
    };
    let needed = context.induction_start(law) - context.first_index;
    let rows = match context.rows(needed.max(context.seeds.len())) {
        Ok(rows) => rows,
        Err(reason) => return Ok(context.not_applicable(reason)),
    };
    context.certify(law, route, ratio, &rows)
}

fn attempt_line(certificate: &LeadingCoefficientCertificate) -> String {
    let route = certificate.route.map_or("none", |r| r.lean_route_name());
    let law = certificate
        .degree_law
        .as_ref()
        .map_or_else(String::new, |law| format!(" with {}", law.text()));
    let ratio = certificate.ratio.as_ref().map_or_else(String::new, |rho| {
        format!(", rho = {}", format_rational_coeff(rho))
    });
    let open: Vec<String> = certificate
        .hypotheses
        .iter()
        .filter(|h| h.status != LiftHypothesisStatus::Verified)
        .map(|h| format!("{}: {}", h.name, h.evidence))
        .take(2)
        .collect();
    format!("{route}{law}{ratio}: {}", open.join("; "))
}

// ---------------------------------------------------------------------------
// The rule
// ---------------------------------------------------------------------------

struct ShapeTerm {
    lag: usize,
    derivative_order: usize,
    alternating: bool,
    coefficients: Bivariate,
}

struct Shape {
    terms: Vec<ShapeTerm>,
    /// `f(n)`.
    left_factor: Vec<Q>,
    max_lag: usize,
    alternating: bool,
}

fn rule_shape(rule: &Recurrence) -> Result<Shape, String> {
    if rule.inhomogeneous.is_some() {
        return Err("the rule has an inhomogeneous term".into());
    }
    let left_factor = match &rule.denominator {
        None => vec![Q::one()],
        Some(denominator) => {
            let coefficients = denominator.normalized_coefficients();
            if (1..x_width(&coefficients)).any(|j| !column(&coefficients, j).is_empty()) {
                return Err("the left-hand factor f(n, x) depends on x".into());
            }
            let f = column(&coefficients, 0);
            if f.is_empty() {
                return Err("the left-hand factor is zero".into());
            }
            f
        }
    };
    let mut terms = Vec::new();
    for term in &rule.terms {
        let coefficients = term.coeff.normalized_coefficients();
        if x_degree(&coefficients).is_none() {
            continue;
        }
        if term.offset == 0 {
            return Err("a term refers to P_n itself (lag 0)".into());
        }
        terms.push(ShapeTerm {
            lag: term.offset,
            derivative_order: term.deriv_order,
            alternating: term.sign == RecurrenceSign::AlternatingN,
            coefficients,
        });
    }
    if terms.is_empty() {
        return Err("the rule has no nonzero term".into());
    }
    let max_lag = terms.iter().map(|t| t.lag).max().unwrap_or(1);
    let alternating = terms.iter().any(|t| t.alternating);
    Ok(Shape {
        terms,
        left_factor,
        max_lag,
        alternating,
    })
}

fn x_width(b: &Bivariate) -> usize {
    b.iter().map(Vec::len).max().unwrap_or(0)
}

/// `[x^j] b(n, x)` as a polynomial in `n`, without trailing zeros.
fn column(b: &Bivariate, j: usize) -> Vec<Q> {
    let mut column: Vec<Q> = b
        .iter()
        .map(|row| row.get(j).cloned().unwrap_or_else(Q::zero))
        .collect();
    vec_poly::strip_trailing_zeros(&mut column);
    column
}

fn x_degree(b: &Bivariate) -> Option<usize> {
    (0..x_width(b)).rev().find(|&j| !column(b, j).is_empty())
}

/// `p(c + period m)` as a polynomial in `m`.
fn on_class(p: &[Q], c: usize, period: usize) -> Vec<Q> {
    if p.is_empty() {
        return Vec::new();
    }
    vec_poly::compose(p, &[integer(c as i64), integer(period as i64)]).stripped()
}

/// Polynomials without trailing zeros (the zero polynomial is empty).
trait StripTrailingZeros {
    fn stripped(self) -> Self;
}

impl StripTrailingZeros for Vec<Q> {
    fn stripped(mut self) -> Self {
        vec_poly::strip_trailing_zeros(&mut self);
        self
    }
}

fn integer(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

/// `(base + growth m)_d = prod_(i < d) (base - i + growth m)` in `m`.
fn falling_factorial_on_class(base: usize, growth: usize, order: usize) -> Vec<Q> {
    let mut product = vec![Q::one()];
    for i in 0..order {
        product = vec_poly::mul(
            &product,
            &[integer(base as i64 - i as i64), integer(growth as i64)],
        );
    }
    product
}

fn leading_coefficient(row: &[Q]) -> Option<&Q> {
    row.iter().rev().find(|c| !c.is_zero())
}

fn row_degree(row: &[Q]) -> Option<usize> {
    row.iter().rposition(|c| !c.is_zero())
}

// ---------------------------------------------------------------------------
// Degree-law inference
// ---------------------------------------------------------------------------

/// Candidate laws read off the rows `rows[i] = P_(first + i)`: by increasing
/// period, the earliest start that fits every later sampled row and the next
/// few starts, with at least two periods and the largest lag of evidence.
fn infer_degree_laws(
    rows: &[Vec<Q>],
    first: usize,
    maximum_period: usize,
    max_lag: usize,
    alternating: bool,
) -> Vec<RecurrenceDegreeLaw> {
    let degrees: Vec<Option<usize>> = rows.iter().map(|r| row_degree(r)).collect();
    let signs: Vec<Option<LeadingCoefficientSign>> = rows
        .iter()
        .map(|r| leading_coefficient(r).and_then(LeadingCoefficientSign::of))
        .collect();
    let earliest = degrees
        .iter()
        .rposition(Option::is_none)
        .map_or(0, |zero| zero + 1);
    let mut laws = Vec::new();
    for period in 1..=maximum_period {
        if alternating && period % 2 == 1 {
            continue;
        }
        let evidence = 2 * period + max_lag + 2;
        if rows.len() < earliest + evidence {
            continue;
        }
        let last = rows.len() - 1;
        let growth = degrees[last].unwrap() as i64 - degrees[last - period].unwrap() as i64;
        if growth < 0 {
            continue;
        }
        // The earliest start of the suffix on which the law fits.
        let mut start = last - period;
        while start > earliest {
            let i = start - 1;
            let fits = degrees[i + period].unwrap() as i64 - degrees[i].unwrap() as i64 == growth
                && signs[i + period] == signs[i];
            if !fits {
                break;
            }
            start = i;
        }
        if start + evidence > rows.len() {
            continue;
        }
        for start in start..=(start + EXTRA_LAW_STARTS).min(rows.len() - evidence) {
            laws.push(RecurrenceDegreeLaw {
                start: first + start,
                period,
                base_degrees: (0..period).map(|r| degrees[start + r].unwrap()).collect(),
                growth_per_period: growth as usize,
                leading_coefficient_signs: (0..period).map(|r| signs[start + r].unwrap()).collect(),
            });
        }
    }
    laws
}

// ---------------------------------------------------------------------------
// The checks
// ---------------------------------------------------------------------------

struct Context {
    record: RecurrenceJson,
    options: LeadingCoefficientOptions,
    rule: Recurrence,
    first_index: usize,
    seeds: Vec<Vec<Q>>,
    shape: Shape,
    normal_form: String,
}

/// One residue class `n = c + period m`, `m >= 0`.
struct Class {
    c: usize,
    /// `e_j = D(c) - D(c - j)` for every lag of the rule.
    gaps: BTreeMap<usize, i64>,
}

impl Context {
    /// The context, or a not-applicable certificate.
    #[allow(clippy::type_complexity)]
    fn new(
        rule: &Recurrence,
        first_index: usize,
        seeds: &[Vec<Q>],
        options: &LeadingCoefficientOptions,
    ) -> Result<Result<Self, Box<LeadingCoefficientCertificate>>, RecurrenceError> {
        // The record in its normal form, so that a replay reproduces it.
        let (supplied, first_index, seeds) =
            RecurrenceJson::from_recurrence_rational(rule, first_index, seeds, None)
                .to_recurrence_parts()?;
        let record = RecurrenceJson::from_recurrence_rational(&supplied, first_index, &seeds, None);
        let empty_certificate = |reason: String| {
            Box::new(LeadingCoefficientCertificate {
                schema: LEADING_COEFFICIENT_SCHEMA.into(),
                recurrence: record.clone(),
                options: *options,
                status: LeadingCoefficientStatus::NotApplicable,
                normal_form: reason,
                route: None,
                degree_law: None,
                degree_is_exact: true,
                ratio: None,
                induction_start: first_index + seeds.len(),
                initial_rows: Vec::new(),
                hypotheses: Vec::new(),
                holds_from: None,
                lean_route: None,
                lean_declarations: Vec::new(),
                conclusion: None,
                exact_leading_coefficient: None,
                exact_leading_coefficient_hypotheses: Vec::new(),
                attempts: Vec::new(),
            })
        };
        let normalized = match rule_with_x_free_left_hand_factor(&supplied) {
            Ok((normalized, _)) => normalized,
            Err(error) => {
                let reason = match error {
                    RecurrenceError::Unsupported(reason) => reason,
                    other => other.to_string(),
                };
                return Ok(Err(empty_certificate(reason)));
            }
        };
        let shape = match rule_shape(&normalized) {
            Ok(shape) => shape,
            Err(reason) => return Ok(Err(empty_certificate(reason))),
        };
        if seeds.len() < shape.max_lag {
            return Ok(Err(empty_certificate(format!(
                "the rule has lag {} and needs as many seed rows",
                shape.max_lag
            ))));
        }
        let normal_form = normal_form_text(&shape, first_index + seeds.len());
        Ok(Ok(Self {
            record,
            options: *options,
            rule: normalized,
            first_index,
            seeds,
            shape,
            normal_form,
        }))
    }

    fn not_applicable(&self, reason: String) -> LeadingCoefficientCertificate {
        LeadingCoefficientCertificate {
            schema: LEADING_COEFFICIENT_SCHEMA.into(),
            recurrence: self.record.clone(),
            options: self.options,
            status: LeadingCoefficientStatus::NotApplicable,
            normal_form: reason,
            route: None,
            degree_law: None,
            degree_is_exact: true,
            ratio: None,
            induction_start: self.rule_start(),
            initial_rows: Vec::new(),
            hypotheses: Vec::new(),
            holds_from: None,
            lean_route: None,
            lean_declarations: Vec::new(),
            conclusion: None,
            exact_leading_coefficient: None,
            exact_leading_coefficient_hypotheses: Vec::new(),
            attempts: Vec::new(),
        }
    }

    /// `N`, the first index generated by the rule.
    fn rule_start(&self) -> usize {
        self.first_index + self.seeds.len()
    }

    fn sample_count(&self) -> usize {
        self.options.sample_rows.max(
            self.seeds.len() + 3 * self.options.maximum_period.max(1) + 2 * self.shape.max_lag + 4,
        )
    }

    fn rows(&self, count: usize) -> Result<Vec<Vec<Q>>, String> {
        self.rule
            .generate_rows_rational(&self.seeds, self.first_index, count)
            .map(|rows| {
                rows.into_iter()
                    .map(|mut row| {
                        vec_poly::strip_trailing_zeros(&mut row);
                        row
                    })
                    .collect()
            })
            .map_err(|error| format!("the rows cannot be generated: {error}"))
    }

    fn induction_start(&self, law: &RecurrenceDegreeLaw) -> usize {
        self.rule_start().max(law.start + self.shape.max_lag)
    }

    /// The routes worth trying for `law`, in order.
    fn route_candidates(
        &self,
        law: &RecurrenceDegreeLaw,
        rows: &[Vec<Q>],
    ) -> Vec<(LeadingCoefficientRoute, Option<Q>)> {
        let mut routes = vec![(LeadingCoefficientRoute::Cone, None)];
        for rho in self.ratio_candidates(law, rows) {
            routes.push((LeadingCoefficientRoute::Recurrence, Some(rho)));
        }
        if self.shape.max_lag == 1
            && law.has_constant_sign(LeadingCoefficientSign::Positive)
            && self
                .shape
                .terms
                .iter()
                .all(|t| t.derivative_order <= 1 && !t.alternating)
        {
            routes.push((LeadingCoefficientRoute::Cancel, None));
            for rho in self.cancel_ratio_candidates(law, rows) {
                routes.push((LeadingCoefficientRoute::Cancel, Some(rho)));
            }
        }
        routes
    }

    /// Exact candidates for `rho`: the smallest observed ratio
    /// `mu(n) / mu(n - 1)` past the induction start, its dyadic roundings
    /// down, `1` and `1/2`.
    fn ratio_candidates(&self, law: &RecurrenceDegreeLaw, rows: &[Vec<Q>]) -> Vec<Q> {
        let start = self.induction_start(law).saturating_sub(self.shape.max_lag);
        let mut values = Vec::new();
        for n in start.max(law.start)..self.first_index + rows.len() {
            let row = &rows[n - self.first_index];
            let (Some(lead), Some(sign)) = (leading_coefficient(row), law.sign_at(n)) else {
                return Vec::new();
            };
            values.push(lead * sign.factor());
        }
        let ratios: Vec<Q> = values
            .windows(2)
            .filter(|w| w[0].is_positive())
            .map(|w| &w[1] / &w[0])
            .collect();
        let Some(minimum) = ratios.iter().min().cloned() else {
            return Vec::new();
        };
        if !minimum.is_positive() {
            return Vec::new();
        }
        let mut candidates = rounded_down(&minimum);
        candidates.push(Q::one());
        candidates.push(Q::new(BigInt::one(), BigInt::from(2)));
        let mut seen = BTreeSet::new();
        candidates.retain(|c| seen.insert(c.clone()));
        candidates
    }

    /// Exact candidates for the `rho` of the `cancel` route: the smallest
    /// observed `[x^(D-1)] P_n / lead P_n` from the row before the induction
    /// start on, and its dyadic roundings down.
    fn cancel_ratio_candidates(&self, law: &RecurrenceDegreeLaw, rows: &[Vec<Q>]) -> Vec<Q> {
        let start = self.induction_start(law) - 1;
        let mut minimum: Option<Q> = None;
        for n in start..self.first_index + rows.len() {
            let row = &rows[n - self.first_index];
            let Some(degree) = row_degree(row).filter(|&d| d >= 1) else {
                return Vec::new();
            };
            let value = &row[degree - 1] / &row[degree];
            if minimum.as_ref().is_none_or(|current| value < *current) {
                minimum = Some(value);
            }
        }
        match minimum {
            Some(minimum) if minimum.is_positive() => rounded_down(&minimum),
            _ => Vec::new(),
        }
    }

    /// The affine envelope of the `coefficients` route: growth
    /// `g = max ceil((deg_x c_t - d_t) / j_t)`, start `N - K`, and the
    /// smallest `D0` covering the seed window.
    fn envelope_law(&self, rows: &[Vec<Q>]) -> Option<RecurrenceDegreeLaw> {
        if self.shape.alternating || self.shape.terms.iter().any(|t| t.derivative_order > 1) {
            return None;
        }
        let mut growth = 0i64;
        for term in &self.shape.terms {
            let excess = x_degree(&term.coefficients)? as i64 - term.derivative_order as i64;
            let lag = term.lag as i64;
            growth = growth.max((excess + lag - 1).div_euclid(lag));
        }
        let start = self.rule_start() - self.shape.max_lag;
        let mut base = 0i64;
        for (offset, row) in rows[start - self.first_index..self.rule_start() - self.first_index]
            .iter()
            .enumerate()
        {
            base = base.max(row_degree(row)? as i64 - growth * offset as i64);
        }
        Some(RecurrenceDegreeLaw {
            start,
            period: 1,
            base_degrees: vec![base.max(0) as usize],
            growth_per_period: growth as usize,
            leading_coefficient_signs: vec![LeadingCoefficientSign::Positive],
        })
    }

    fn classes(&self, law: &RecurrenceDegreeLaw) -> Vec<Class> {
        let start = self.induction_start(law);
        let lags: BTreeSet<usize> = self.shape.terms.iter().map(|t| t.lag).collect();
        (start..start + law.period)
            .map(|c| Class {
                c,
                gaps: lags
                    .iter()
                    .map(|&j| {
                        (
                            j,
                            law.degree_at(c).unwrap() as i64 - law.degree_at(c - j).unwrap() as i64,
                        )
                    })
                    .collect(),
            })
            .collect()
    }

    /// The summand coefficients on a class, with the factor `(-1)^n`
    /// resolved (the period is even when a term alternates).
    fn class_terms(&self, class: &Class) -> Vec<(usize, usize, Bivariate)> {
        self.shape
            .terms
            .iter()
            .map(|t| {
                let mut coefficients = t.coefficients.clone();
                if t.alternating && class.c % 2 == 1 {
                    for row in &mut coefficients {
                        for value in row {
                            *value = -value.clone();
                        }
                    }
                }
                (t.lag, t.derivative_order, coefficients)
            })
            .collect()
    }

    /// `f(c + period m)` in `m`.
    fn left_factor_on(&self, class: &Class, law: &RecurrenceDegreeLaw) -> Vec<Q> {
        on_class(&self.shape.left_factor, class.c, law.period)
    }

    /// The certificate of `route` for `law`.
    fn certify(
        &self,
        law: &RecurrenceDegreeLaw,
        route: LeadingCoefficientRoute,
        ratio: Option<&Q>,
        rows: &[Vec<Q>],
    ) -> Result<LeadingCoefficientCertificate, RecurrenceError> {
        let start = self.induction_start(law);
        let mut certificate = self.not_applicable(self.normal_form.clone());
        certificate.status = LeadingCoefficientStatus::NotProved;
        certificate.route = Some(route);
        certificate.degree_law = Some(law.clone());
        certificate.degree_is_exact = route != LeadingCoefficientRoute::Coefficients;
        certificate.ratio = ratio.cloned();
        certificate.induction_start = start;
        if law.start < self.first_index
            || (self.shape.alternating && law.period % 2 == 1)
            || rows.len() < start - self.first_index
        {
            return Err(RecurrenceError::InvalidInput(
                "the degree law starts before the first index, has an odd period for a rule \
                 with (-1)^n factors, or the rows are missing"
                    .into(),
            ));
        }
        let initial = &rows[..start - self.first_index];
        certificate.initial_rows = initial
            .iter()
            .enumerate()
            .map(|(i, row)| LeadingCoefficientRow {
                index: self.first_index + i,
                degree: row_degree(row),
                sign: leading_coefficient(row).and_then(LeadingCoefficientSign::of),
            })
            .collect();
        let mut hypotheses = Vec::new();
        hypotheses.push(self.initial_rows_hypothesis(law, route, ratio, initial));
        if self.shape.left_factor.len() > 1 {
            hypotheses.push(sign_hypothesis(
                "left_factor_nonzero",
                &format!("f(n) != 0 for every n >= {start}"),
                "f(n)",
                on_class(&self.shape.left_factor, start, 1),
                vec![Q::one()],
                start,
                1,
                IndexSignCondition::Nonzero,
            ));
        }
        let classes = self.classes(law);
        match route {
            LeadingCoefficientRoute::Cone
            | LeadingCoefficientRoute::Subseq
            | LeadingCoefficientRoute::TwoStep => {
                for class in &classes {
                    hypotheses.push(self.degree_bound_hypothesis(law, class, 0));
                    hypotheses.extend(self.cone_hypotheses(law, class));
                }
            }
            LeadingCoefficientRoute::Recurrence => {
                let Some(rho) = ratio.filter(|rho| rho.is_positive()) else {
                    return Err(RecurrenceError::InvalidInput(
                        "the recurrence route needs a positive ratio".into(),
                    ));
                };
                for class in &classes {
                    hypotheses.push(self.degree_bound_hypothesis(law, class, 0));
                    hypotheses.extend(self.ratio_hypotheses(law, class, rho));
                }
            }
            LeadingCoefficientRoute::Cancel => {
                if ratio.is_some_and(|rho| !rho.is_positive()) {
                    return Err(RecurrenceError::InvalidInput(
                        "the ratio of the cancel route must be positive".into(),
                    ));
                }
                hypotheses.extend(self.cancel_hypotheses(law, &classes, ratio));
            }
            LeadingCoefficientRoute::Coefficients => {
                for class in &classes {
                    hypotheses.push(self.degree_bound_hypothesis(law, class, 0));
                    hypotheses.extend(self.coefficient_hypotheses(law, class, true));
                }
            }
        }
        // The route label must fit the rule.
        let label_fits = match route {
            LeadingCoefficientRoute::Subseq => {
                law.period > 1 && self.shape.terms.iter().all(|t| t.lag % law.period == 0)
            }
            LeadingCoefficientRoute::TwoStep => {
                self.is_two_step() && law.has_constant_sign(LeadingCoefficientSign::Positive)
            }
            // Nonnegative coefficients: at most first derivatives, no
            // `(-1)^n` factor, positive leading coefficients (and first
            // order for `cancel`).
            LeadingCoefficientRoute::Cancel | LeadingCoefficientRoute::Coefficients => {
                !self.shape.alternating
                    && self.shape.terms.iter().all(|t| t.derivative_order <= 1)
                    && law.has_constant_sign(LeadingCoefficientSign::Positive)
                    && (route == LeadingCoefficientRoute::Coefficients || self.shape.max_lag == 1)
            }
            _ => true,
        };
        if !label_fits {
            hypotheses.push(LiftHypothesis {
                name: "route_shape".into(),
                statement: format!(
                    "the rule has the shape of the route {}",
                    route.lean_route_name()
                ),
                status: LiftHypothesisStatus::Failed { witness: None },
                evidence: "the rule does not have this shape".into(),
                sign_decision: None,
            });
        }
        certificate.hypotheses = hypotheses;
        if certificate
            .hypotheses
            .iter()
            .all(|h| h.status == LiftHypothesisStatus::Verified)
        {
            certificate.status = LeadingCoefficientStatus::Proved;
            let route = self.refined_route(law, route);
            certificate.route = Some(route);
            certificate.lean_route = Some(route.lean_route_name().into());
            certificate.lean_declarations = lean_declarations(route, law, ratio.is_some());
            let holds_from = certificate
                .initial_rows
                .iter()
                .rposition(|row| row.degree.is_none())
                .map_or(self.first_index, |i| self.first_index + i + 1);
            certificate.holds_from = Some(holds_from);
            certificate.conclusion = Some(self.conclusion(law, route, holds_from, &certificate));
            if certificate.degree_is_exact {
                if let Some((value, hypotheses)) = self.exact_leading_coefficient(law, rows) {
                    if let Some(conclusion) = certificate.conclusion.as_mut() {
                        conclusion.push_str(&format!(
                            "; lead P_n = {} for every n >= {}",
                            format_rational_coeff(&value),
                            start - self.shape.max_lag
                        ));
                    }
                    certificate.exact_leading_coefficient = Some(value);
                    certificate.exact_leading_coefficient_hypotheses = hypotheses;
                }
            }
        }
        Ok(certificate)
    }

    /// The constant leading coefficient `c` (see
    /// [`LeadingCoefficientCertificate::exact_leading_coefficient`]) with
    /// its verified hypotheses, when every one holds.
    fn exact_leading_coefficient(
        &self,
        law: &RecurrenceDegreeLaw,
        rows: &[Vec<Q>],
    ) -> Option<(Q, Vec<LiftHypothesis>)> {
        if !(law.has_constant_sign(LeadingCoefficientSign::Positive)
            || law.has_constant_sign(LeadingCoefficientSign::Negative))
        {
            return None;
        }
        let start = self.induction_start(law);
        let window = start - self.shape.max_lag;
        let value = leading_coefficient(&rows[window - self.first_index])?.clone();
        if (window..start).any(|n| leading_coefficient(&rows[n - self.first_index]) != Some(&value))
        {
            return None;
        }
        let mut hypotheses = Vec::new();
        for class in &self.classes(law) {
            // The top-coefficient recurrence needs the degree bound on every
            // class (it fails where the `cancel` route's top terms cancel:
            // there the coefficient also involves `[x^(D-1)] P_(n-1)`).
            let bound = self.degree_bound_hypothesis(law, class, 0);
            if bound.status != LiftHypothesisStatus::Verified {
                return None;
            }
            hypotheses.push(bound);
            let denominator = self.left_factor_on(class, law);
            let total = self
                .signed_multipliers(law, class)
                .iter()
                .fold(Vec::new(), |sum, m| vec_poly::add(&sum, m));
            let hypothesis = sign_hypothesis(
                &format!("top_multiplier_sum_is_one_class_{}", class.c),
                &format!(
                    "sum_j M_j(n) = f(n) for n = {} + {} m, so lead P_n = lead P_(n-j) for \
                     every lag",
                    class.c, law.period
                ),
                "sum_j M_j - f",
                vec_poly::sub(&total, &denominator).stripped(),
                vec![Q::one()],
                class.c,
                law.period,
                IndexSignCondition::Zero,
            );
            if hypothesis.status != LiftHypothesisStatus::Verified {
                return None;
            }
            hypotheses.push(hypothesis);
        }
        Some((value, hypotheses))
    }

    fn is_two_step(&self) -> bool {
        matches!(self.shape.terms.as_slice(), [term] if term.lag == 2
            && term.derivative_order == 0
            && !term.alternating
            && term.coefficients.iter().skip(1).all(|row| row.iter().all(Zero::is_zero)))
            && self.shape.left_factor.len() == 1
    }

    /// `cone` becomes `twoStep` or `subseq` when the rule has that shape.
    fn refined_route(
        &self,
        law: &RecurrenceDegreeLaw,
        route: LeadingCoefficientRoute,
    ) -> LeadingCoefficientRoute {
        match route {
            LeadingCoefficientRoute::Cone
                if self.is_two_step()
                    && law.has_constant_sign(LeadingCoefficientSign::Positive) =>
            {
                LeadingCoefficientRoute::TwoStep
            }
            LeadingCoefficientRoute::Cone
                if law.period > 1 && self.shape.terms.iter().all(|t| t.lag % law.period == 0) =>
            {
                LeadingCoefficientRoute::Subseq
            }
            other => other,
        }
    }

    fn conclusion(
        &self,
        law: &RecurrenceDegreeLaw,
        route: LeadingCoefficientRoute,
        holds_from: usize,
        certificate: &LeadingCoefficientCertificate,
    ) -> String {
        let relation = if certificate.degree_is_exact {
            "the degree law"
        } else {
            "the degree envelope (deg P_n <= D(n)), nonnegative coefficients and P_n(0) > 0, with"
        };
        let mut text = format!(
            "{relation} {} (route {}, induction from n = {})",
            law.text(),
            route.lean_route_name(),
            certificate.induction_start
        );
        if law.start > self.first_index {
            let earlier: Vec<String> = certificate
                .initial_rows
                .iter()
                .filter(|row| row.index < law.start)
                .map(|row| match (row.degree, row.sign) {
                    (Some(degree), Some(sign)) => {
                        format!("P_{}: degree {degree}, {}", row.index, sign.text())
                    }
                    _ => format!("P_{} = 0", row.index),
                })
                .collect();
            text.push_str(&format!("; earlier rows (exact): {}", earlier.join(", ")));
        }
        if holds_from > self.first_index {
            text.push_str(&format!("; every row from P_{holds_from} on is nonzero"));
        } else {
            text.push_str("; every row is nonzero");
        }
        text
    }

    // -- hypotheses -------------------------------------------------------

    fn initial_rows_hypothesis(
        &self,
        law: &RecurrenceDegreeLaw,
        route: LeadingCoefficientRoute,
        ratio: Option<&Q>,
        initial: &[Vec<Q>],
    ) -> LiftHypothesis {
        let start = self.induction_start(law);
        let window = start - self.shape.max_lag;
        let exact = route != LeadingCoefficientRoute::Coefficients;
        let needs_nonnegative = matches!(
            route,
            LeadingCoefficientRoute::Coefficients | LeadingCoefficientRoute::Cancel
        );
        let statement = format!(
            "the rows P_{}, ..., P_{} computed exactly: every row P_n with n >= {} has {} D(n) and \
             the sign of the law{}{}",
            self.first_index,
            start - 1,
            law.start,
            if exact { "degree" } else { "degree at most" },
            if needs_nonnegative {
                format!("; the rows from P_{window} on have nonnegative coefficients")
            } else {
                String::new()
            },
            if route == LeadingCoefficientRoute::Coefficients {
                " and a positive constant term"
            } else {
                ""
            },
        );
        let mut failure: Option<(usize, String)> = None;
        for (i, row) in initial.iter().enumerate() {
            let n = self.first_index + i;
            let Some(expected) = law.degree_at(n) else {
                continue;
            };
            let degree = row_degree(row);
            let sign = leading_coefficient(row).and_then(LeadingCoefficientSign::of);
            let degree_fits = match degree {
                None => false,
                Some(d) if exact => d == expected,
                Some(d) => d <= expected,
            };
            if !degree_fits || sign != law.sign_at(n) {
                failure = Some((n, format!("P_{n} has degree {degree:?} and sign {sign:?}")));
                break;
            }
            if needs_nonnegative && n >= window {
                if row.iter().any(Signed::is_negative) {
                    failure = Some((n, format!("P_{n} has a negative coefficient")));
                    break;
                }
                if route == LeadingCoefficientRoute::Coefficients
                    && !row.first().is_some_and(Signed::is_positive)
                {
                    failure = Some((n, format!("P_{n}(0) is not positive")));
                    break;
                }
            }
        }
        if failure.is_none() {
            if let (LeadingCoefficientRoute::Cancel, Some(rho)) = (route, ratio) {
                let n = start - 1;
                let row = &initial[n - self.first_index];
                let holds = row_degree(row)
                    .filter(|&d| d >= 1)
                    .is_some_and(|d| rho * &row[d] <= row[d - 1]);
                if !holds {
                    failure = Some((
                        n,
                        format!("the ratio base fails: [x^(D-1)] P_{n} < rho lead P_{n}"),
                    ));
                }
            }
            if let (LeadingCoefficientRoute::Recurrence, Some(rho)) = (route, ratio) {
                let values: Vec<Q> = (window..start)
                    .map(|n| {
                        let row = &initial[n - self.first_index];
                        leading_coefficient(row).cloned().unwrap_or_else(Q::zero)
                            * law.sign_at(n).map_or_else(Q::zero, |s| s.factor())
                    })
                    .collect();
                if let Some(t) = values.windows(2).position(|w| rho * &w[0] > w[1]) {
                    failure = Some((
                        window + t + 1,
                        format!(
                            "the ratio base fails: rho mu({}) > mu({}) for mu(n) = sign * lead P_n",
                            window + t,
                            window + t + 1
                        ),
                    ));
                }
            }
        }
        let (status, evidence) = match failure {
            None => (
                LiftHypothesisStatus::Verified,
                format!(
                    "checked exactly on {} rows{}",
                    initial.len(),
                    match (route, ratio) {
                        (LeadingCoefficientRoute::Recurrence, Some(rho)) => format!(
                            "; rho mu(t) <= mu(t + 1) for t = {window}, ..., {} with rho = {}",
                            start.saturating_sub(2),
                            format_rational_coeff(rho)
                        ),
                        (LeadingCoefficientRoute::Cancel, Some(rho)) => format!(
                            "; [x^(D-1)] P_{0} >= rho lead P_{0} with rho = {1}",
                            start - 1,
                            format_rational_coeff(rho)
                        ),
                        _ => String::new(),
                    }
                ),
            ),
            Some((n, reason)) => (
                LiftHypothesisStatus::Failed {
                    witness: Some(BigInt::from(n)),
                },
                format!("{reason} (exact)"),
            ),
        };
        LiftHypothesis {
            name: "initial_rows".into(),
            statement,
            status,
            evidence,
            sign_decision: None,
        }
    }

    /// The degree bound `deg_x c_t <= e_j + d_t + extra` and `e_j >= 0` on
    /// a class.
    fn degree_bound_hypothesis(
        &self,
        law: &RecurrenceDegreeLaw,
        class: &Class,
        extra: i64,
    ) -> LiftHypothesis {
        let c = class.c;
        let name = format!("degree_bound_class_{c}");
        let statement = format!(
            "for n = {c} + {} m: D(n - j) <= D(n) for every lag j, and deg_x c_t(n, x) <= \
             D(n) - D(n - j_t) + d_t{} for every summand",
            law.period,
            if extra > 0 { " + 1" } else { "" }
        );
        let mut failure = None;
        for (lag, order, coefficients) in self.class_terms(class) {
            let gap = class.gaps[&lag];
            if gap < 0 {
                failure = Some((None, format!("D({c}) < D({}) (lag {lag})", c - lag)));
                break;
            }
            let allowed = gap + order as i64 + extra;
            let top = x_degree(&coefficients).unwrap_or(0) as i64;
            if top > allowed {
                // The coefficient of x^top is a nonzero polynomial in n: it
                // is nonzero at some n of the class.
                let excess = column(&coefficients, top as usize);
                let witness = (0..=excess.len()).find_map(|m| {
                    let n = c + law.period * m;
                    (!vec_poly::evaluate(&excess, &integer(n as i64)).is_zero()).then_some(n)
                });
                failure = Some((
                    witness,
                    format!(
                        "the summand of lag {lag} and derivative order {order} has degree \
                         {top} > {allowed} in x"
                    ),
                ));
                break;
            }
        }
        let (status, evidence) = match failure {
            None => (
                LiftHypothesisStatus::Verified,
                format!(
                    "gaps {} (exact)",
                    class
                        .gaps
                        .iter()
                        .map(|(j, e)| format!("e_{j} = {e}"))
                        .collect::<Vec<_>>()
                        .join(", ")
                ),
            ),
            Some((witness, reason)) => (
                LiftHypothesisStatus::Failed {
                    witness: witness.map(BigInt::from),
                },
                format!("{reason} (exact)"),
            ),
        };
        LiftHypothesis {
            name,
            statement,
            status,
            evidence,
            sign_decision: None,
        }
    }

    /// The numerators (over `f(c + period m)`) of the signed multipliers
    /// `m_j` on a class, by lag `1..=K`.
    fn signed_multipliers(&self, law: &RecurrenceDegreeLaw, class: &Class) -> Vec<Vec<Q>> {
        let c = class.c;
        let mut multipliers = vec![Vec::new(); self.shape.max_lag + 1];
        for (lag, order, coefficients) in self.class_terms(class) {
            let gap = class.gaps[&lag];
            if gap < 0 {
                continue;
            }
            let top = column(&coefficients, gap as usize + order);
            let factor = falling_factorial_on_class(
                law.degree_at(c - lag).unwrap(),
                law.growth_per_period,
                order,
            );
            let sign = law.sign_at(c).unwrap().factor() * law.sign_at(c - lag).unwrap().factor();
            let summand = vec_poly::scale(
                &vec_poly::mul(&on_class(&top, c, law.period), &factor),
                &sign,
            );
            multipliers[lag] = vec_poly::add(&multipliers[lag], &summand).stripped();
        }
        multipliers
    }

    fn cone_hypotheses(&self, law: &RecurrenceDegreeLaw, class: &Class) -> Vec<LiftHypothesis> {
        let c = class.c;
        let denominator = self.left_factor_on(class, law);
        let multipliers = self.signed_multipliers(law, class);
        let mut hypotheses = Vec::new();
        let mut total = Vec::new();
        for (lag, numerator) in multipliers.iter().enumerate().skip(1) {
            if !class.gaps.contains_key(&lag) {
                continue;
            }
            total = vec_poly::add(&total, numerator);
            hypotheses.push(sign_hypothesis(
                &format!("top_multiplier_lag_{lag}_nonnegative_class_{c}"),
                &format!(
                    "the signed top-coefficient multiplier m_{lag}(n) >= 0 for n = {c} + {} m",
                    law.period
                ),
                &format!("m_{lag}"),
                numerator.clone(),
                denominator.clone(),
                c,
                law.period,
                IndexSignCondition::Nonnegative,
            ));
        }
        hypotheses.push(sign_hypothesis(
            &format!("top_multiplier_sum_positive_class_{c}"),
            &format!(
                "sum_j m_j(n) > 0 for n = {c} + {} m (the top coefficients stay positive)",
                law.period
            ),
            "sum_j m_j",
            total.stripped(),
            denominator,
            c,
            law.period,
            IndexSignCondition::Positive,
        ));
        hypotheses
    }

    fn ratio_hypotheses(
        &self,
        law: &RecurrenceDegreeLaw,
        class: &Class,
        rho: &Q,
    ) -> Vec<LiftHypothesis> {
        let c = class.c;
        let k = self.shape.max_lag;
        let denominator = self.left_factor_on(class, law);
        let multipliers = self.signed_multipliers(law, class);
        let power = |e: usize| (0..e).fold(Q::one(), |acc, _| acc * rho);
        (0..k)
            .map(|i| {
                // C_i = sum_(t = i)^(K-1) m_(K-t) rho^(t-i) - rho^(K-i), over f.
                let mut numerator = vec_poly::scale(&denominator, &-power(k - i));
                for t in i..k {
                    numerator = vec_poly::add(
                        &numerator,
                        &vec_poly::scale(&multipliers[k - t], &power(t - i)),
                    );
                }
                sign_hypothesis(
                    &format!("ratio_step_{i}_class_{c}"),
                    &format!(
                        "C_{i}(n) = sum_(t = {i})^{} m_(K-t)(n) rho^(t-{i}) - rho^{} >= 0 for \
                         n = {c} + {} m, rho = {} (the ratio invariant mu(n) >= rho mu(n - 1) \
                         propagates)",
                        k - 1,
                        k - i,
                        law.period,
                        format_rational_coeff(rho)
                    ),
                    &format!("C_{i}"),
                    numerator.stripped(),
                    denominator.clone(),
                    c,
                    law.period,
                    IndexSignCondition::Nonnegative,
                )
            })
            .collect()
    }

    /// The pairs `(B_j, A_j)` of non-derivative and first-derivative
    /// coefficients by lag, on a class.
    fn paired_terms(&self, class: &Class) -> BTreeMap<usize, (Bivariate, Bivariate)> {
        let mut pairs: BTreeMap<usize, (Bivariate, Bivariate)> = BTreeMap::new();
        for (lag, order, coefficients) in self.class_terms(class) {
            let entry = pairs.entry(lag).or_default();
            let target = if order == 0 {
                &mut entry.0
            } else {
                &mut entry.1
            };
            *target = add_bivariate(target, &coefficients);
        }
        pairs
    }

    /// Nonnegative coefficients: `[x^0] A_j >= 0` and `B_(j,i) + r A_(j,i+1)
    /// >= 0` at `r = 0` and `r = D(n - j)`; with `witness`, also
    /// `sum_j [x^0] B_j > 0`.
    fn coefficient_hypotheses(
        &self,
        law: &RecurrenceDegreeLaw,
        class: &Class,
        witness: bool,
    ) -> Vec<LiftHypothesis> {
        let c = class.c;
        let p = law.period;
        let denominator = self.left_factor_on(class, law);
        let mut hypotheses = Vec::new();
        let mut constant_terms = Vec::new();
        for (lag, (b, a)) in self.paired_terms(class) {
            let on = |poly: Vec<Q>| on_class(&poly, c, p);
            let previous = vec![
                integer(law.degree_at(c - lag).unwrap() as i64),
                integer(law.growth_per_period as i64),
            ];
            let a0 = on(column(&a, 0));
            if !a0.is_empty() {
                hypotheses.push(sign_hypothesis(
                    &format!("derivative_constant_nonnegative_lag_{lag}_class_{c}"),
                    &format!("[x^0] A_{lag}(n) >= 0 for n = {c} + {p} m"),
                    &format!("[x^0] A_{lag}"),
                    a0,
                    denominator.clone(),
                    c,
                    p,
                    IndexSignCondition::Nonnegative,
                ));
            }
            let top = x_degree(&b).max(x_degree(&a).map(|d| d.saturating_sub(1)));
            for i in 0..=top.unwrap_or(0) {
                let bi = on(column(&b, i));
                let ai = on(column(&a, i + 1));
                if bi.is_empty() && ai.is_empty() {
                    continue;
                }
                if i == 0 {
                    constant_terms = vec_poly::add(&constant_terms, &bi);
                }
                if !bi.is_empty() {
                    hypotheses.push(sign_hypothesis(
                        &format!("coefficient_multiplier_{i}_low_lag_{lag}_class_{c}"),
                        &format!("[x^{i}] B_{lag}(n) >= 0 for n = {c} + {p} m (r = 0)"),
                        &format!("[x^{i}] B_{lag}"),
                        bi.clone(),
                        denominator.clone(),
                        c,
                        p,
                        IndexSignCondition::Nonnegative,
                    ));
                }
                if !ai.is_empty() {
                    let upper = vec_poly::add(&bi, &vec_poly::mul(&ai, &previous));
                    hypotheses.push(sign_hypothesis(
                        &format!("coefficient_multiplier_{i}_high_lag_{lag}_class_{c}"),
                        &format!(
                            "[x^{i}] B_{lag}(n) + D(n - {lag}) [x^{}] A_{lag}(n) >= 0 for n = {c} \
                             + {p} m (r = D(n - {lag}))",
                            i + 1
                        ),
                        &format!("[x^{i}] B_{lag} + D(n - {lag}) [x^{}] A_{lag}", i + 1),
                        upper.stripped(),
                        denominator.clone(),
                        c,
                        p,
                        IndexSignCondition::Nonnegative,
                    ));
                }
            }
        }
        if witness {
            hypotheses.push(sign_hypothesis(
                &format!("constant_term_multiplier_positive_class_{c}"),
                &format!(
                    "sum_j [x^0] B_j(n) > 0 for n = {c} + {p} m (so P_n(0) >= sum_j [x^0] B_j(n) \
                     P_(n-j)(0) > 0)"
                ),
                "sum_j [x^0] B_j",
                constant_terms.stripped(),
                denominator,
                c,
                p,
                IndexSignCondition::Positive,
            ));
        }
        hypotheses
    }

    /// The `cancel` route on every class (first-order rules), with the
    /// coefficient-ratio invariant `[x^(D-1)] P_n >= rho lead P_n` when
    /// `ratio` is given (module documentation).
    fn cancel_hypotheses(
        &self,
        law: &RecurrenceDegreeLaw,
        classes: &[Class],
        ratio: Option<&Q>,
    ) -> Vec<LiftHypothesis> {
        let mut hypotheses = Vec::new();
        let Some(first_class) = classes.first() else {
            return hypotheses;
        };
        let pairs = self.paired_terms(first_class);
        let (b, a) = pairs.get(&1).cloned().unwrap_or_default();
        let m = x_degree(&b)
            .unwrap_or(0)
            .max(x_degree(&a).map_or(0, |d| d.saturating_sub(1)));
        let shape_ok = m >= 1;
        hypotheses.push(LiftHypothesis {
            name: "cancel_shape".into(),
            statement:
                "P_n = A_n(x) P'_(n-1) + B_n(x) P_(n-1) with deg_x A <= m + 1, deg_x B <= m \
                        for some m >= 1"
                    .into(),
            status: if shape_ok {
                LiftHypothesisStatus::Verified
            } else {
                LiftHypothesisStatus::Failed { witness: None }
            },
            evidence: format!("m = {m} (exact)"),
            sign_decision: None,
        });
        if !shape_ok {
            return hypotheses;
        }
        if ratio.is_some() {
            // The invariant needs a coefficient below the leading one.
            let lowest = classes
                .iter()
                .map(|class| (class.c - 1, law.degree_at(class.c - 1).unwrap()))
                .min_by_key(|(_, degree)| *degree);
            let positive = lowest.is_some_and(|(_, degree)| degree >= 1);
            hypotheses.push(LiftHypothesis {
                name: "cancel_positive_degree".into(),
                statement: "D(n - 1) >= 1 for every n from the induction start on (the ratio \
                            invariant compares [x^(D-1)] P_(n-1) with lead P_(n-1))"
                    .into(),
                status: if positive {
                    LiftHypothesisStatus::Verified
                } else {
                    LiftHypothesisStatus::Failed {
                        witness: lowest.map(|(n, _)| BigInt::from(n)),
                    }
                },
                evidence: format!(
                    "the law is nondecreasing; smallest D(n - 1) on the first classes: {:?} \
                     (exact)",
                    lowest.map(|(_, degree)| degree)
                ),
                sign_decision: None,
            });
        }
        for class in classes {
            let c = class.c;
            let p = law.period;
            let denominator = self.left_factor_on(class, law);
            let on = |poly: Vec<Q>| on_class(&poly, c, p);
            let previous = vec![
                integer(law.degree_at(c - 1).unwrap() as i64),
                integer(law.growth_per_period as i64),
            ];
            // E_i(D(n-1) - shift) = B_i + (D(n-1) - shift) A_(i+1), with
            // B_(-1) = 0: the coefficient of P_(n-1)[x^(D - shift)] in
            // [x^(D - shift + i)] P_n (with the factor 1 / f(n)).
            let level = |i: i64, shift: i64| -> Vec<Q> {
                let r = vec_poly::add(&previous, &[integer(-shift)]).stripped();
                let b_part = if i >= 0 {
                    on(column(&b, i as usize))
                } else {
                    Vec::new()
                };
                vec_poly::add(
                    &b_part,
                    &vec_poly::mul(&on(column(&a, (i + 1) as usize)), &r),
                )
                .stripped()
            };
            let a_top = on(column(&a, m + 1));
            let m_signed = m as i64;
            let top = level(m_signed, 0);
            let next = level(m_signed - 1, 0);
            let step = class.gaps[&1];
            let sign =
                |name: &str, statement: String, function: &str, numerator: Vec<Q>, condition| {
                    sign_hypothesis(
                        &format!("{name}_class_{c}"),
                        &format!("{statement} for n = {c} + {p} m'"),
                        function,
                        numerator,
                        denominator.clone(),
                        c,
                        p,
                        condition,
                    )
                };
            // `alpha p_D + beta p_(D-1) >= 0` on the cone `p_D > 0`,
            // `p_(D-1) >= rho p_D`: `beta >= 0` and `alpha + rho beta >= 0`.
            let invariant = |alpha: Vec<Q>, beta: Vec<Q>, rho: &Q, text: &str| {
                vec![
                    sign(
                        "ratio_invariant_slope",
                        format!(
                            "the coefficient of [x^(D-1)] P_(n-1) in rho^-1 [x^(D'-1)] P_n - \
                             lead P_n is >= 0 ({text})"
                        ),
                        "beta",
                        beta.clone(),
                        IndexSignCondition::Nonnegative,
                    ),
                    sign(
                        "ratio_invariant",
                        format!(
                            "alpha + rho beta >= 0 with rho = {}, so [x^(D'-1)] P_n >= rho lead \
                             P_n from [x^(D-1)] P_(n-1) >= rho lead P_(n-1) ({text})",
                            format_rational_coeff(rho)
                        ),
                        "alpha + rho beta",
                        vec_poly::add(&alpha, &vec_poly::scale(&beta, rho)).stripped(),
                        IndexSignCondition::Nonnegative,
                    ),
                ]
            };
            if step == m_signed {
                hypotheses.push(sign(
                    "top_grows",
                    format!("A_(m+1)(n) D(n-1) + B_m(n) > 0 (step m = {m})"),
                    "A_(m+1) D + B_m",
                    top.clone(),
                    IndexSignCondition::Positive,
                ));
                if let Some(rho) = ratio {
                    // lead P_n = E_m(D) p_D, [x^(D'-1)] P_n = E_(m-1)(D) p_D +
                    // E_m(D-1) p_(D-1).
                    let alpha = vec_poly::sub(&next, &vec_poly::scale(&top, rho));
                    hypotheses.extend(invariant(
                        alpha,
                        level(m_signed, 1),
                        rho,
                        "the top terms grow",
                    ));
                }
            } else if step + 1 == m_signed {
                hypotheses.push(sign(
                    "top_cancels",
                    "A_(m+1)(n) D(n-1) + B_m(n) = 0 (step m - 1)".into(),
                    "A_(m+1) D + B_m",
                    top,
                    IndexSignCondition::Zero,
                ));
                hypotheses.push(sign(
                    "top_derivative_nonpositive",
                    "A_(m+1)(n) <= 0".into(),
                    "A_(m+1)",
                    a_top.clone(),
                    IndexSignCondition::Nonpositive,
                ));
                match ratio {
                    None => hypotheses.push(sign(
                        "next_level_positive",
                        "A_m(n) D(n-1) + B_(m-1)(n) > 0".into(),
                        "A_m D + B_(m-1)",
                        next,
                        IndexSignCondition::Positive,
                    )),
                    Some(rho) => {
                        // lead P_n = E_(m-1)(D) p_D - A_(m+1) p_(D-1) > 0 on
                        // the cone; [x^(D'-1)] P_n = E_(m-2)(D) p_D +
                        // E_(m-1)(D-1) p_(D-1) - 2 A_(m+1) p_(D-2), the last
                        // term >= 0.
                        hypotheses.push(sign(
                            "next_level_with_ratio_positive",
                            format!(
                                "A_m(n) D(n-1) + B_(m-1)(n) - rho A_(m+1)(n) > 0 with rho = {}",
                                format_rational_coeff(rho)
                            ),
                            "A_m D + B_(m-1) - rho A_(m+1)",
                            vec_poly::sub(&next, &vec_poly::scale(&a_top, rho)).stripped(),
                            IndexSignCondition::Positive,
                        ));
                        let alpha =
                            vec_poly::sub(&level(m_signed - 2, 0), &vec_poly::scale(&next, rho));
                        let beta =
                            vec_poly::add(&level(m_signed - 1, 1), &vec_poly::scale(&a_top, rho));
                        hypotheses.extend(invariant(alpha, beta, rho, "the top terms cancel"));
                    }
                }
            } else {
                hypotheses.push(LiftHypothesis {
                    name: format!("cancel_step_class_{c}"),
                    statement: format!("the degree step on the class of {c} is m or m - 1"),
                    status: LiftHypothesisStatus::Failed {
                        witness: Some(BigInt::from(c)),
                    },
                    evidence: format!("the step is {step}, m = {m} (exact)"),
                    sign_decision: None,
                });
            }
            hypotheses.extend(self.coefficient_hypotheses(law, class, false));
        }
        hypotheses
    }
}

/// `value` and its roundings down to multiples of `1, 1/2, 1/4, 1/8` that
/// stay positive, without repetitions.
fn rounded_down(value: &Q) -> Vec<Q> {
    let mut candidates = vec![value.clone()];
    for k in 0..4u32 {
        let scale = integer(1 << k);
        let rounded = (value * &scale).floor() / scale;
        if rounded.is_positive() && !candidates.contains(&rounded) {
            candidates.push(rounded);
        }
    }
    candidates
}

fn add_bivariate(a: &Bivariate, b: &Bivariate) -> Bivariate {
    let rows = a.len().max(b.len());
    let width = x_width(a).max(x_width(b));
    (0..rows)
        .map(|i| {
            (0..width)
                .map(|j| {
                    let get = |m: &Bivariate| {
                        m.get(i)
                            .and_then(|row| row.get(j))
                            .cloned()
                            .unwrap_or_else(Q::zero)
                    };
                    get(a) + get(b)
                })
                .collect()
        })
        .collect()
}

/// A sign hypothesis on the class `n = c + period m`, decided for `m >= 0`;
/// a failure witness is reported as the index `n`.
#[allow(clippy::too_many_arguments)]
fn sign_hypothesis(
    name: &str,
    statement: &str,
    function: &str,
    numerator: Vec<Q>,
    denominator: Vec<Q>,
    c: usize,
    period: usize,
    condition: IndexSignCondition,
) -> LiftHypothesis {
    let decision = decide_index_rational_function_sign(
        &numerator,
        &denominator,
        &IntegerRange::starting_at(0),
        condition,
    );
    let status = match &decision {
        IndexSignDecision::Holds { .. } => LiftHypothesisStatus::Verified,
        IndexSignDecision::FailsAt { index, .. } => LiftHypothesisStatus::Failed {
            witness: Some(BigInt::from(c) + BigInt::from(period) * index),
        },
        IndexSignDecision::Undecided { .. } => LiftHypothesisStatus::Unresolved,
    };
    // The checker counts in `m`; a failure is reported at its index `n`.
    let evidence = match &decision {
        IndexSignDecision::FailsAt { index, value } => format!(
            "fails at m = {index}, that is n = {}: {} (exact)",
            BigInt::from(c) + BigInt::from(period) * index,
            match value {
                Some(value) => format!("{function} = {}", format_rational_coeff(value)),
                None => "the denominator vanishes".to_string(),
            }
        ),
        _ => decision.describe(&format!("{function}(m)")),
    };
    LiftHypothesis {
        name: name.into(),
        statement: statement.into(),
        status,
        evidence: format!("in m with n = {c} + {period} m: {evidence}"),
        sign_decision: Some(decision),
    }
}

fn lean_declarations(
    route: LeadingCoefficientRoute,
    law: &RecurrenceDegreeLaw,
    ratio: bool,
) -> Vec<String> {
    let affine = law.affine_step().is_some();
    let mut names: Vec<&str> = match route {
        LeadingCoefficientRoute::Cone if affine => {
            vec!["RealRooted.LinRec.natDegree_eq_and_leadingCoeff_pos"]
        }
        LeadingCoefficientRoute::Cone => {
            vec!["RealRooted.LinRec.natDegree_eq_and_leadingCoeff_pos_of_degreeLaw"]
        }
        LeadingCoefficientRoute::Recurrence => vec![
            "RealRooted.LinRec.natDegree_eq_and_leadingCoeff_pos_of_ratio",
            "RealRooted.ScalarRec.pos_of_ratio",
        ],
        LeadingCoefficientRoute::Subseq => vec![
            "RealRooted.LinRec.natDegree_eq_and_leadingCoeff_pos_of_degreeLaw",
            "RealRooted.LinRec.degreeLaw_periodic_decomp",
        ],
        LeadingCoefficientRoute::TwoStep => vec![
            "RealRooted.twoStepProduct_natDegree",
            "RealRooted.twoStepProduct_leadingCoeff_pos",
        ],
        LeadingCoefficientRoute::Cancel if ratio => {
            vec!["RealRooted.LinRec.natDegree_eq_and_leadingCoeff_pos_of_cancel"]
        }
        LeadingCoefficientRoute::Cancel => {
            vec!["RealRooted.derivativeRecurrence_natDegree_eq_and_leadingCoeff_pos_of_cancel"]
        }
        LeadingCoefficientRoute::Coefficients => vec![
            "RealRooted.LinRec.natDegree_le_and_hasNonnegCoeffs",
            "RealRooted.HasNonnegCoeffs.pos_leadingCoeff",
        ],
    };
    if !law.has_constant_sign(LeadingCoefficientSign::Positive) {
        names.push("RealRooted.hasPosLeadingCoeff_C_mul");
    }
    names.into_iter().map(String::from).collect()
}

fn normal_form_text(shape: &Shape, start: usize) -> String {
    let terms: Vec<String> = shape
        .terms
        .iter()
        .map(|t| {
            let derivative = match t.derivative_order {
                0 => String::new(),
                1 => "'".into(),
                d => format!("^({d})"),
            };
            format!(
                "{}c(n, x) P{derivative}_(n-{}) [deg_x c = {}]",
                if t.alternating { "(-1)^n " } else { "" },
                t.lag,
                x_degree(&t.coefficients).unwrap_or(0)
            )
        })
        .collect();
    let factor = if shape.left_factor.len() > 1 {
        "f(n) "
    } else {
        ""
    };
    format!("{factor}P_n = {} for n >= {start}", terms.join(" + "))
}
