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

mod shape;
use shape::*;
mod certification;
mod context;
mod hypotheses;
