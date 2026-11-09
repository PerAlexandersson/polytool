//! Index-section lifts: the family is a (possibly reversed) every-`p`-th
//! term section of a parent family with a constant-coefficient rule.  This
//! is the general, exact form of the RealRooted row route `pathPower` (and
//! of the `CubicParent` rows), and the inverse of the residue-class route
//! `subseq`.
//!
//! # The theorem
//!
//! Let the rule be `P_n = sum_(k=1)^r a_k(x) P_(n-k)` for `n >= N = n_0 + K`
//! (first index `n_0`, `K >= r` seed rows), with multipliers that depend on
//! `x` only (no index `n`, no derivative, no `(-1)^n`, a nonzero constant
//! left-hand factor divided out) and order `r` (`a_r != 0`).  The rows
//! `P_n`, `n >= n_0 + s_0` with `s_0 = K - r`, follow the pure order-`r`
//! rule; the rows before are exceptional.
//!
//! **Transformed family.**  For `s` in `{s_0, s_0 + 1}` (the number of
//! exceptional prefix rows) put
//!
//! - identity: `S_m = P_(n_0 + s + m)`;
//! - reversal: `S_m(x) = x^(d m + e) P_(n_0 + s + m)(1/x)`, with `d` the
//!   smallest integer such that `deg a_k <= d k` for every `k` (the degree
//!   slope of the rule) and `e >= 0` the smallest integer such that
//!   `deg P_(n_0 + s + m) <= d m + e` for `m < r`.
//!
//! Then `S_m = sum_k b_k(x) S_(m-k)` for every `m >= r`, with `b_k = a_k`
//! (identity) or `b_k(x) = x^(d k) a_k(1/x)` (reversal; a polynomial by the
//! choice of `d`): multiply the rule at `n = n_0 + s + m >= N` by
//! `x^(d m + e)` and substitute `1/x`.  The seeds `S_0, ..., S_(r-1)` are
//! polynomials (choice of `e`), hence so is every `S_m`.  The characteristic
//! polynomial of the transformed rule is `chi_T(z) = z^r - sum_k b_k(x)
//! z^(r-k)`, monic in `z` over `Q[x]`.
//!
//! **Parent.**  A parent of period `p >= 2` is a monic
//! `psi(w) = w^r - sum_k c_k(x) w^(r-k)` in `Q[x][w]` with
//!
//! ```text
//! Res_w(psi(w), z - w^p) = chi_T(z)
//! ```
//!
//! and polynomial seeds `Q_0, ..., Q_(r-1)` of the parent family
//! `Q_n = sum_k c_k(x) Q_(n-k)`, `n >= r`, such that `Q_(p m) = S_m` for
//! `m < r`.
//!
//! **The identity for every index** (an equal-order recurrence argument,
//! nothing fitted).  Let `C` be the companion matrix of `psi` acting on the
//! state `X_n = (Q_n, ..., Q_(n+r-1))`, so that `X_(n+1) = C X_n` for every
//! `n >= 0` and `X_(p m) = M^m X_0` with `M = C^p`.  The eigenvalues of `M`
//! are the `p`-th powers `alpha_i^p` of the roots of `psi`, so the
//! characteristic polynomial of `M` is `prod_i (z - alpha_i^p) =
//! Res_w(psi(w), z - w^p) = chi_T(z)`.  By Cayley--Hamilton,
//! `M^m = sum_k b_k M^(m-k)` for `m >= r`, and the first coordinate gives
//! `Q_(p m) = sum_k b_k Q_(p (m-k))`: the section `m -> Q_(p m)` satisfies
//! the order-`r` rule of `S` for every `m >= r`.  It agrees with `S` on the
//! `r` seeds `m < r` (checked exactly), so `Q_(p m) = S_m` for every
//! `m >= 0` by induction on `m`.  RealRooted proves the same subsequence
//! recurrence as `LinearRecurrence.ofMatrix_pow_isSolution_subsequence`.
//! The resultant is computed exactly from the power sums of the roots of
//! `psi` (Newton's identities over `Q[x]`): the power sums of the
//! `alpha_i^p` are the power sums of the `alpha_i` at multiples of `p`.
//!
//! **Real-rootedness.**  When a theorem of the parent menu proves that every
//! `Q_n`, `n >= F`, is nonzero and real-rooted, so is every `S_m` with
//! `p m >= F` (a subsequence), and so is `P_(n_0 + s + m)`: for the
//! reversal `P_(n_0 + s + m)(x) = x^(d m + e) S_m(1/x)`, whose nonzero roots
//! are the reciprocals of the nonzero roots of `S_m` and which has a root at
//! `0` exactly when `deg S_m < d m + e`.  The exceptional rows
//! `P_(n_0), ..., P_(n_0 + s - 1)` and the section rows with `p m < F` are
//! checked exactly.  The conclusion: every `P_n`, `n >= n_0`, is nonzero and
//! real-rooted.
//!
//! **Root region.**  For the identity, the parent's index-independent region
//! is the family's.  For the reversal, a parent region `(a, b)` inside
//! `(-inf, 0]` becomes `(1/b, 0]` (lower end `-inf` when `b = 0`), and a
//! region inside `[0, inf)` becomes `[0, 1/a)`; the end at `0` is closed,
//! since a degree drop of `S_m` gives `P_n(0) = 0` (as for OEIS A116088,
//! whose rows have the factor `x`); a parent region with `0` in its interior
//! gives no region.  A region is claimed only when every exactly checked row
//! lies in it.
//!
//! **Interlacing** is not concluded: consecutive rows of `P` are `p` parent
//! rows apart, and the parent theorems state the interlacing of consecutive
//! parent rows only.
//!
//! # The search (bounded, deterministic)
//!
//! For `p = 2, ..., MAX_INDEX_SECTION_PERIOD` (5) and `T` = identity, then
//! reversal, skipping `r p > 24` and the pairs for which `(-1)^r chi_T(0)` is
//! not a `p`-th power in `Q[x]` (it equals `((-1)^r psi(0))^p` for every
//! parent, so no parent exists then): factor `chi_T(w^p)` over `Q[x]`
//! ([`crate::factorization::irreducible_factorization_monic_bivariate_with_budget`],
//! recombination budget 4096) and enumerate its monic divisors `psi` of
//! `w`-degree `r` (products of the irreducible factors up to their
//! multiplicities, by increasing `x`-degree, then exponent vector).  At most
//! 64 divisors are examined in the whole search (each one resultant check;
//! the count is recorded as
//! [`IndexSectionCertificate::divisors_examined`]).  For each divisor with
//! the resultant identity and each `s`, the `r x r` system
//! `e_0^T M^m v = S_m`, `m < r`, is solved over `Q(x)` by fraction-free
//! elimination and must have a polynomial solution `v` (the parent seeds).
//! The first 6 section identities go to the parent menu; the first one whose
//! parent is proved is selected, and otherwise the first section identity is
//! reported with its open parent (status `unresolved` or
//! `hypothesis_failed`).  A rule with no section identity in this range is
//! `not_applicable`.
//!
//! The parent menu is product induction, three-term, Liu--Wang, Liu--Wang
//! derivative lag, Ma--Wang, Liu--Wang three-term, iterated preserver,
//! hypergeometric row, the fixed-delay window and residue classes (a
//! residue-class proof counts only without zero classes, whose rows stay
//! zero), each
//! with the induction-start retries of the options (recorded in
//! [`IndexSectionCertificate::parent_maximum_start_delay`]); for the
//! identity the requested root region is passed to the parent checks.
//!
//! # RealRooted
//!
//! The route `pathPower` (`RealRooted/Tactic/Row/PathPower.lean`) recognizes
//! direct, reflected (`RealRooted.Tactic.reflectedRecurrence_isSolution`),
//! shifted and padded (`RealRooted.Tactic.prependedRecurrence_isSolution`)
//! subsequences of one fixed parent, the path-power clique polynomials
//! (`RealRooted.splits_pathPowerCliquePolynomial`,
//! `RealRooted.splits_reciprocalShift_pathPowerCliquePolynomial`), and the
//! `CubicParent` rows (`RealRooted.CubicParent.splits_of_even_recurrence`,
//! `splits_of_even_reflected_recurrence`) use a Tran--Zumba parent.  Here
//! the parent is any constant-coefficient family found by factoring
//! `chi_T(w^p)` and proved by any theorem of the parent menu, and the
//! offset of the section is zero without loss of generality (a section
//! `Q_(p m + c)` is the section `Q'_(p m)` of the shifted parent
//! `Q'_n = Q_(n+c)`, with the same rule and polynomial seeds).

use super::super::veronese_algebra::{
    bivariate_to_coefficients, degree_in_second_variable, degree_in_shift, divide_monic_in_shift,
    format_shift_polynomial, BivariatePolynomial, UnivariatePolynomial,
};
use super::super::veronese_search::{
    collect_exponent_vectors, parent_recurrence_from_characteristic, polynomial_row,
    solve_section_system, SeedSolve,
};
use super::*;
use crate::factorization::irreducible_factorization_monic_bivariate_with_budget;
use crate::polynomial::Polynomial;

const REFERENCE: &str = "Index-section lift: if Q_n = sum c_k(x) Q_(n-k) has the \
     characteristic polynomial psi(w) with Res_w(psi(w), z - w^p) = chi(z), then m -> Q_(p m) \
     satisfies the rule of chi (Cayley-Hamilton for the p-th power of the companion matrix); a \
     subsequence of a real-rooted family is real-rooted and the reversal x^D S(1/x) preserves \
     real-rootedness; RealRooted route pathPower (RealRooted/Tactic/Row/PathPower.lean, \
     LinearRecurrence.ofMatrix_pow_isSolution_subsequence)";

/// Largest period `p` tried by [`check_index_section_lift_theorem`].
pub const MAX_INDEX_SECTION_PERIOD: usize = 5;

/// Largest `w`-degree `r p` of a factored `chi_T(w^p)`.
const MAXIMUM_SECTION_DEGREE: usize = 24;

/// Divisors `psi` examined (resultant checks) in the whole search.
const MAXIMUM_DIVISORS: usize = 64;

/// Section identities whose parent is sent to the parent menu.
const MAXIMUM_PARENT_PROOFS: usize = 6;

/// Recombination budget of the bivariate factorization.
const RECOMBINATION_BUDGET: usize = 4096;

/// Rows past the seeds compared as a defensive self-check of a section
/// identity (not part of the proof).
const SELF_CHECK_ROWS: usize = 3;

/// Largest induction-start delay a recorded certificate may ask the parent
/// checks to try.
pub(super) const MAXIMUM_PARENT_START_DELAY: usize = 16;

/// The theorems tried on the parent, in order.
const PARENT_MENU: [RealRootednessTheorem; 10] = [
    RealRootednessTheorem::ProductRecurrence,
    RealRootednessTheorem::ThreeTermRecurrence,
    RealRootednessTheorem::LiuWang,
    RealRootednessTheorem::LiuWangDerivativeLag,
    RealRootednessTheorem::MaWangInterval,
    RealRootednessTheorem::LiuWangThreeTerm,
    RealRootednessTheorem::IteratedPreserver,
    RealRootednessTheorem::HypergeometricRow,
    RealRootednessTheorem::VariableDelay,
    RealRootednessTheorem::ResidueClassSubsequences,
];

/// The transformation `T` of an [`IndexSectionCertificate`] (module
/// documentation, "Transformed family").
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum IndexSectionTransform {
    /// `S_m = P_(n_0 + s + m)`.
    Identity,
    /// `S_m(x) = x^(d m + e) P_(n_0 + s + m)(1/x)`.
    Reversal {
        /// `d`: every multiplier satisfies `deg a_k <= d k`.
        degree_slope: usize,
        /// `e`: `deg P_(n_0 + s + m) <= d m + e` for the seeds `m < r`.
        degree_offset: usize,
    },
}

/// The exact data of a [`check_index_section_lift_theorem`] report: the
/// section identity `P_(n_0 + s + m) = T(Q_(p m))` and the parent's
/// theorem reports.  [`IndexSectionCertificate::verify`] recomputes it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IndexSectionCertificate {
    /// The period `p >= 2`.
    pub period: usize,
    /// The transformation `T`.
    pub transform: IndexSectionTransform,
    /// `s`: the rows `P_(n_0), ..., P_(n_0 + s - 1)` precede the section
    /// family (checked exactly).
    pub exceptional_rows: usize,
    /// The parent family `Q`: the rule `Q_n = sum_k c_k(x) Q_(n-k)` of
    /// `psi(w) = w^r - sum_k c_k(x) w^(r-k)`, first index `0`, the seeds
    /// `Q_0, ..., Q_(r-1)`.
    pub parent: RecurrenceJson,
    /// Divisors `psi` examined by the search (deterministic work count).
    pub divisors_examined: usize,
    /// The largest induction-start delay the parent checks tried; replayed
    /// by [`RealRootednessTheoremReport::verify`].
    pub parent_maximum_start_delay: usize,
    /// The parent theorem checks that apply (status other than
    /// `not_applicable`), in menu order, up to the selected one.
    pub parent_attempts: Vec<RealRootednessTheoremReport>,
    /// Position in `parent_attempts` of the parent's proof, when one proves
    /// it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub selected: Option<usize>,
}

impl IndexSectionCertificate {
    /// Recompute the section identity for the family `rule`, `first_index`,
    /// `seeds` from the certificate alone: the rule has `x`-only multipliers
    /// of order `r`, the transform gives polynomial multipliers and seeds,
    /// the parent rule is the constant-coefficient rule of a monic `psi` of
    /// degree `r` with `Res_w(psi(w), z - w^p) = chi_T(z)` and
    /// `psi | chi_T(w^p)`, the parent seeds give `Q_(p m) = S_m` for `m < r`,
    /// and the selected parent proof is a proof of the recorded parent that
    /// replays ([`RealRootednessTheoremReport::verify`]).
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::VerificationFailed`] when a step fails; errors of
    /// the records themselves.
    pub fn verify(
        &self,
        rule: &Recurrence,
        first_index: usize,
        seeds: &[Vec<Q>],
    ) -> Result<(), RecurrenceError> {
        let failed = |what: String| Err(RecurrenceError::VerificationFailed(what));
        let multipliers = match constant_coefficient_multipliers(rule) {
            Ok(multipliers) => multipliers,
            Err(reason) => return failed(format!("index-section certificate: {reason}")),
        };
        let order = multipliers.len();
        if !(2..=MAX_INDEX_SECTION_PERIOD).contains(&self.period) {
            return failed(format!(
                "index-section period {} outside 2..={MAX_INDEX_SECTION_PERIOD}",
                self.period
            ));
        }
        let Some(base) = seeds.len().checked_sub(order) else {
            return failed("fewer seed rows than the order of the rule".into());
        };
        if self.exceptional_rows < base || self.exceptional_rows > base + 1 {
            return failed(format!(
                "{} exceptional rows, expected {base} or {}",
                self.exceptional_rows,
                base + 1
            ));
        }
        let transformed = match self.transform {
            IndexSectionTransform::Identity => Some(multipliers.clone()),
            IndexSectionTransform::Reversal { degree_slope, .. } => {
                reversed_multipliers(&multipliers, degree_slope)
            }
        };
        let Some(transformed) = transformed else {
            return failed("a multiplier has degree above the recorded degree slope".into());
        };
        let chi = characteristic_polynomial(&transformed);
        let (parent_rule, parent_first, parent_seeds) = self.parent.to_recurrence_parts()?;
        if parent_first != 0 || parent_seeds.len() != order {
            return failed(format!(
                "the parent record must start at 0 with {order} seed rows"
            ));
        }
        let parent_multipliers = match constant_coefficient_multipliers(&parent_rule) {
            Ok(m) if m.len() == order => m,
            _ => {
                return failed(format!(
                    "the parent rule is not a constant-coefficient rule of order {order}"
                ))
            }
        };
        let psi = characteristic_polynomial(&parent_multipliers);
        if power_characteristic_polynomial(&psi, self.period) != chi {
            return failed("Res_w(psi(w), z - w^p) differs from chi_T(z)".into());
        }
        let divides = divide_monic_in_shift(&substitute_shift_power(&chi, self.period), &psi)
            .is_some_and(|(_, remainder)| remainder.is_zero());
        if !divides {
            return failed("psi(w) does not divide chi_T(w^p)".into());
        }
        let rows = rule
            .generate_rows_rational(seeds, first_index, self.exceptional_rows + order)
            .map_err(|error| RecurrenceError::VerificationFailed(error.to_string()))?;
        let offset = match self.transform {
            IndexSectionTransform::Identity => None,
            IndexSectionTransform::Reversal {
                degree_slope,
                degree_offset,
            } => Some((degree_slope, degree_offset)),
        };
        let Some(targets) = transformed_rows(&rows[self.exceptional_rows..], offset) else {
            return failed("a transformed seed row is not a polynomial".into());
        };
        let parent_rows = parent_rule
            .generate_rows_rational(&parent_seeds, 0, self.period * (order - 1) + 1)
            .map_err(|error| RecurrenceError::VerificationFailed(error.to_string()))?;
        for (m, target) in targets.iter().enumerate() {
            if Polynomial::new(parent_rows[self.period * m].clone()) != *target {
                return failed(format!(
                    "Q_{} differs from the transformed row P_{}",
                    self.period * m,
                    first_index + self.exceptional_rows + m
                ));
            }
        }
        if let Some(position) = self.selected {
            let Some(proof) = self.parent_attempts.get(position) else {
                return failed("the certificate selects a missing parent attempt".into());
            };
            if !proves_every_row(proof)
                || proof.recurrence != self.parent
                || !PARENT_MENU.contains(&proof.theorem)
            {
                return failed("the selected parent attempt is not a proof of the parent".into());
            }
            proof.verify()?;
        }
        Ok(())
    }
}

/// Check the index-section lift (module documentation) for the family given
/// by `rule`, `first_index` and the seed rows `seeds` (the rule applies from
/// `N = first_index + seeds.len()` on).
///
/// It applies when the rule is `P_n = sum_k a_k(x) P_(n-k)` with multipliers
/// free of `n` (no derivative, no `(-1)^n`, a constant left-hand factor) and
/// the bounded search finds a section identity
/// `P_(first_index + s + m) = T(Q_(p m))` for every `m >= 0`, with `T` the
/// identity or a reversal and `Q` a polynomial parent family with a
/// constant-coefficient rule; it is proved when a theorem of the parent menu
/// proves the parent and the exceptional rows check exactly.  The
/// certificate is kept in [`RealRootednessTheoremReport::index_section`] and
/// replayed by [`RealRootednessTheoremReport::verify`].  No interlacing is
/// concluded.
///
/// # Errors
///
/// Errors of the record itself (seed normalization) and internal
/// consistency failures of the search.
///
/// # Example
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::recurrence::{
///     check_index_section_lift_theorem, BigRational, BivarPoly, Recurrence,
///     RecurrenceTerm, RealRootednessTheoremOptions,
/// };
///
/// let q = |v: i64| BigRational::from_integer(BigInt::from(v));
/// let constant = |row: Vec<i64>| {
///     BivarPoly::from_coefficients_normalized(vec![row.into_iter().map(q).collect()])
/// };
/// // P_n = (3 + x) P_(n-1) - 3 P_(n-2) + P_(n-3), P_0 = 1, P_1 = 1 + x,
/// // P_2 = 1 + 4x + x^2 (OEIS A109955): the reversed rows are every third row
/// // of Q_n = Q_(n-1) + x Q_(n-3).
/// let rule = Recurrence {
///     terms: vec![
///         RecurrenceTerm::new(1, 0, constant(vec![3, 1])),
///         RecurrenceTerm::new(2, 0, constant(vec![-3])),
///         RecurrenceTerm::new(3, 0, constant(vec![1])),
///     ],
///     denominator: None,
///     inhomogeneous: None,
/// };
/// let seeds = vec![vec![q(1)], vec![q(1), q(1)], vec![q(1), q(4), q(1)]];
/// let report = check_index_section_lift_theorem(
///     &rule,
///     0,
///     &seeds,
///     &RealRootednessTheoremOptions::default(),
/// )
/// .unwrap();
/// assert!(report.is_proved());
/// assert_eq!(report.index_section.as_ref().unwrap().period, 3);
/// report.verify().unwrap();
/// ```
pub fn check_index_section_lift_theorem(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    super::check_with_first_proving_start(
        RealRootednessTheorem::IndexSectionLift,
        rule,
        first_index,
        seeds,
        options,
    )
}

// ---------------------------------------------------------------------------
// Exact algebra
// ---------------------------------------------------------------------------

/// The multipliers `a_1, ..., a_r` of a rule `P_n = sum_k a_k(x) P_(n-k)`
/// (lags without a term give `0`), or why the rule has another form.
fn constant_coefficient_multipliers(
    rule: &Recurrence,
) -> Result<Vec<UnivariatePolynomial>, String> {
    let shape = recurrence_shape(rule)?;
    if shape.denominator.len() > 1 {
        return Err("the left-hand factor depends on n".into());
    }
    let leading = shape.denominator[0].clone();
    let order = shape.terms.keys().map(|&(lag, _)| lag).max().unwrap_or(0);
    if order == 0 {
        return Err("the rule has no term of positive lag".into());
    }
    let mut multipliers = vec![UnivariatePolynomial::zero(); order];
    for (&(lag, derivative), coefficient) in &shape.terms {
        if lag == 0 {
            return Err("the rule has a term of lag 0".into());
        }
        if derivative != 0 {
            return Err(format!(
                "the rule has a derivative term (lag {lag}, derivative order {derivative})"
            ));
        }
        if coefficient
            .iter()
            .skip(1)
            .any(|row| row.iter().any(|c| !c.is_zero()))
        {
            return Err(format!("the multiplier of lag {lag} depends on n"));
        }
        let row = coefficient.first().cloned().unwrap_or_default();
        multipliers[lag - 1] = Polynomial::new(vec_poly::scale(&row, &leading.recip()));
    }
    Ok(multipliers)
}

/// `z^r - sum_k b_k z^(r-k)` for `b = (b_1, ..., b_r)`.
fn characteristic_polynomial(multipliers: &[UnivariatePolynomial]) -> BivariatePolynomial {
    let order = multipliers.len();
    let mut coefficients = vec![UnivariatePolynomial::zero(); order + 1];
    coefficients[order] = UnivariatePolynomial::one();
    for (lag, multiplier) in multipliers.iter().enumerate() {
        coefficients[order - lag - 1] = -multiplier.clone();
    }
    Polynomial::new(coefficients)
}

/// `chi(w^period)`.
fn substitute_shift_power(chi: &BivariatePolynomial, period: usize) -> BivariatePolynomial {
    let degree = degree_in_shift(chi);
    let mut coefficients = vec![UnivariatePolynomial::zero(); degree * period + 1];
    for (power, coefficient) in chi.coeffs().iter().enumerate() {
        coefficients[power * period] = coefficient.clone();
    }
    Polynomial::new(coefficients)
}

/// `Res_w(psi(w), z - w^period) = prod_i (z - alpha_i^period)` for a monic
/// `psi` with roots `alpha_i`, from the power sums of the roots (Newton's
/// identities over `Q[x]`).
fn power_characteristic_polynomial(
    psi: &BivariatePolynomial,
    period: usize,
) -> BivariatePolynomial {
    let order = degree_in_shift(psi);
    let integer = |k: usize| Q::from_integer(BigInt::from(k));
    // psi = w^r + c_1 w^(r-1) + ... + c_r.
    let c: Vec<UnivariatePolynomial> = (0..=order).map(|j| psi.coeff(order - j)).collect();
    let mut power_sums = vec![UnivariatePolynomial::zero(); order * period + 1];
    for k in 1..=order * period {
        let mut value = UnivariatePolynomial::zero();
        for j in 1..=order.min(k - 1) {
            value = value + c[j].clone() * power_sums[k - j].clone();
        }
        if k <= order {
            value = value + c[k].scale(&integer(k));
        }
        power_sums[k] = -value;
    }
    // chi = z^r + d_1 z^(r-1) + ... + d_r with the power sums t_k = s_(k p).
    let t = |k: usize| power_sums[k * period].clone();
    let mut d = vec![UnivariatePolynomial::one()];
    for k in 1..=order {
        let mut value = t(k);
        for (j, dj) in d.iter().enumerate().skip(1) {
            value = value + dj.clone() * t(k - j);
        }
        d.push(value.scale(&(-integer(k).recip())));
    }
    Polynomial::new((0..=order).map(|power| d[order - power].clone()).collect())
}

/// The number of negative coefficients of the parent multipliers
/// `c_k(x)` of `psi(w) = w^r - sum_k c_k(x) w^(r-k)`.
fn negative_parent_coefficients(psi: &BivariatePolynomial) -> usize {
    let order = degree_in_shift(psi);
    (0..order)
        .flat_map(|power| psi.coeff(power).coeffs().to_vec())
        .filter(|c| c.is_positive())
        .count()
}

/// Whether `f` is `g^power` for some `g` in `Q[x]` (exact: the candidate
/// root is built coefficient by coefficient and its power compared).
fn is_perfect_power(f: &UnivariatePolynomial, power: usize) -> bool {
    let Some(degree) = f.degree() else {
        return true;
    };
    let valuation = f.coeffs().iter().position(|c| !c.is_zero()).unwrap_or(0);
    if degree % power != 0 || valuation % power != 0 {
        return false;
    }
    let shifted: Vec<Q> = f.coeffs()[valuation..].to_vec();
    let Some(lowest) = rational_root(&shifted[0], power) else {
        return false;
    };
    let length = (degree - valuation) / power + 1;
    let mut root = vec![lowest];
    let scale = Q::from_integer(BigInt::from(power)) * root[0].pow(power as i32 - 1);
    for k in 1..length {
        let mut truncated = root.clone();
        truncated.push(Q::zero());
        let mut product = vec![Q::one()];
        for _ in 0..power {
            product = vec_poly::mul(&product, &truncated);
            product.truncate(k + 1);
        }
        let current = product.get(k).cloned().unwrap_or_else(Q::zero);
        let target = shifted.get(k).cloned().unwrap_or_else(Q::zero);
        root.push((target - current) / &scale);
    }
    let mut product = vec![Q::one()];
    for _ in 0..power {
        product = vec_poly::mul(&product, &root);
    }
    vec_poly::trim(&product) == vec_poly::trim(&shifted)
}

/// The rational `power`-th root of `value`, when it is one.
fn rational_root(value: &Q, power: usize) -> Option<Q> {
    if value.is_negative() && power % 2 == 0 {
        return None;
    }
    let root_of = |n: &BigInt| {
        let root = n.abs().nth_root(power as u32);
        (root.pow(power as u32) == n.abs()).then_some(root)
    };
    let numerator = root_of(value.numer())?;
    let denominator = root_of(value.denom())?;
    let root = Q::new(numerator, denominator);
    Some(if value.is_negative() { -root } else { root })
}

/// The smallest `d` with `deg a_k <= d k` for every `k`.
fn degree_slope(multipliers: &[UnivariatePolynomial]) -> usize {
    multipliers
        .iter()
        .enumerate()
        .filter_map(|(lag, a)| a.degree().map(|degree| degree.div_ceil(lag + 1)))
        .max()
        .unwrap_or(0)
}

/// `b_k(x) = x^(d k) a_k(1/x)`, or `None` when some `deg a_k > d k`.
fn reversed_multipliers(
    multipliers: &[UnivariatePolynomial],
    slope: usize,
) -> Option<Vec<UnivariatePolynomial>> {
    multipliers
        .iter()
        .enumerate()
        .map(|(lag, a)| {
            if a.is_zero() {
                Some(UnivariatePolynomial::zero())
            } else {
                a.reverse_with_degree(slope * (lag + 1))
            }
        })
        .collect()
}

/// The transformed rows `S_m` of `rows[m] = P_(n_0 + s + m)`: unchanged,
/// or `x^(d m + e) P(1/x)` for `Some((d, e))`; `None` when a reversed row is
/// not a polynomial.
fn transformed_rows(
    rows: &[Vec<Q>],
    reversal: Option<(usize, usize)>,
) -> Option<Vec<UnivariatePolynomial>> {
    rows.iter()
        .enumerate()
        .map(|(m, row)| {
            let row = Polynomial::new(row.clone());
            match reversal {
                None => Some(row),
                Some(_) if row.is_zero() => Some(row),
                Some((slope, offset)) => row.reverse_with_degree(slope * m + offset),
            }
        })
        .collect()
}

/// The smallest `e >= 0` with `deg rows[m] <= d m + e`.
fn degree_offset(rows: &[Vec<Q>], slope: usize) -> usize {
    rows.iter()
        .enumerate()
        .filter_map(|(m, row)| {
            Polynomial::new(row.clone())
                .degree()
                .map(|degree| degree.saturating_sub(slope * m))
        })
        .max()
        .unwrap_or(0)
}

/// The parent seeds `v` with `Q_(p m) = targets[m]`, `m < r`, when the
/// system has a unique polynomial solution.
fn solve_parent_seeds(
    parent_rule: &Recurrence,
    period: usize,
    targets: &[UnivariatePolynomial],
) -> Result<Option<Vec<UnivariatePolynomial>>, RecurrenceError> {
    let order = targets.len();
    let count = period * (order - 1) + 1;
    let mut basis = Vec::with_capacity(order);
    for unit in 0..order {
        let seeds: Vec<Vec<Q>> = (0..order)
            .map(|j| vec![if j == unit { Q::one() } else { Q::zero() }])
            .collect();
        basis.push(
            parent_rule
                .generate_rows_rational(&seeds, 0, count)
                .map_err(|error| {
                    RecurrenceError::VerificationFailed(format!(
                        "the parent rule cannot be replayed: {error}"
                    ))
                })?,
        );
    }
    let matrix: Vec<Vec<UnivariatePolynomial>> = (0..order)
        .map(|m| {
            basis
                .iter()
                .map(|rows| Polynomial::new(rows[period * m].clone()))
                .collect()
        })
        .collect();
    Ok(match solve_section_system(&matrix, targets, order)? {
        SeedSolve::Solved(solution) => Some(solution),
        SeedSolve::NotDetermined | SeedSolve::NotPolynomial => None,
    })
}

// ---------------------------------------------------------------------------
// Search
// ---------------------------------------------------------------------------

/// The search order of a divisor `psi`: negative parent coefficients,
/// `x`-degree, exponent vector.
type DivisorOrder = (usize, usize, Vec<usize>);

/// A section identity found by the search, with its parent attempts.
struct SectionLift {
    period: usize,
    transform: IndexSectionTransform,
    exceptional_rows: usize,
    psi: BivariatePolynomial,
    parent_rule: Recurrence,
    parent_seeds: Vec<Vec<Q>>,
    attempts: Vec<RealRootednessTheoremReport>,
    selected: Option<usize>,
}

impl SectionLift {
    fn is_proved(&self) -> bool {
        self.selected.is_some()
    }
}

/// The parent menu on the parent family: the applying attempts up to the
/// selected one (a proof with an index-independent root region preferred).
fn prove_parent(
    rule: &Recurrence,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> (Vec<RealRootednessTheoremReport>, Option<usize>) {
    let mut attempts = Vec::new();
    let mut fallback = None;
    for theorem in PARENT_MENU {
        let Ok(attempt) = check_with_first_proving_start(theorem, rule, 0, seeds, options) else {
            continue;
        };
        if attempt.status == RealRootednessTheoremStatus::NotApplicable {
            continue;
        }
        let proved = proves_every_row(&attempt);
        let with_region = proved_region(&attempt).is_some();
        attempts.push(attempt);
        if proved && (options.root_region.is_some() || with_region) {
            let position = attempts.len() - 1;
            return (attempts, Some(position));
        }
        if proved && fallback.is_none() {
            fallback = Some(attempts.len() - 1);
        }
    }
    (attempts, fallback)
}

/// The bounded search (module documentation, "The search"): the selected
/// section lift (the first with a proved parent, else the first found) and
/// the number of divisors examined.
fn search_section_lifts(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
    multipliers: &[UnivariatePolynomial],
    options: &RealRootednessTheoremOptions,
) -> Result<(Option<SectionLift>, usize), RecurrenceError> {
    let order = multipliers.len();
    let base = seeds.len() - order;
    let Ok(rows) =
        rule.generate_rows_rational(seeds, first_index, base + 1 + order + SELF_CHECK_ROWS)
    else {
        return Ok((None, 0));
    };
    let slope = degree_slope(multipliers);
    let reversed = reversed_multipliers(multipliers, slope).expect("the slope bounds every degree");
    let mut examined = 0usize;
    let mut parent_proofs = 0usize;
    let mut first_found: Option<SectionLift> = None;
    'search: for period in 2..=MAX_INDEX_SECTION_PERIOD {
        if order * period > MAXIMUM_SECTION_DEGREE {
            break;
        }
        for reversal in [false, true] {
            let transformed = if reversal { &reversed } else { multipliers };
            let chi = characteristic_polynomial(transformed);
            // Necessary: (-1)^r chi_T(0) = ((-1)^r psi(0))^p is a p-th power
            // in Q[x]; skipping the factorization otherwise is exact.
            let constant = chi.coeff(0);
            let constant = if order % 2 == 1 { -constant } else { constant };
            if !is_perfect_power(&constant, period) {
                continue;
            }
            let substituted = substitute_shift_power(&chi, period);
            let Ok(factors) = irreducible_factorization_monic_bivariate_with_budget(
                &bivariate_to_coefficients(&substituted),
                RECOMBINATION_BUDGET,
            ) else {
                continue;
            };
            let factor_polynomials: Vec<BivariatePolynomial> = factors
                .iter()
                .map(|factor| {
                    Polynomial::new(
                        factor
                            .coefficients
                            .iter()
                            .map(|row| Polynomial::new(row.clone()))
                            .collect(),
                    )
                })
                .collect();
            let shift_degrees: Vec<usize> =
                factor_polynomials.iter().map(degree_in_shift).collect();
            let x_degrees: Vec<usize> = factor_polynomials
                .iter()
                .map(degree_in_second_variable)
                .collect();
            let multiplicities: Vec<usize> = factors.iter().map(|f| f.multiplicity).collect();
            let mut exponent_vectors = Vec::new();
            collect_exponent_vectors(
                &shift_degrees,
                &multiplicities,
                order,
                &mut vec![0; factors.len()],
                0,
                &mut exponent_vectors,
            );
            // Divisors whose parent rule has fewer negative coefficients come
            // first (a parent with nonnegative multipliers is the likely
            // provable one), then by x-degree and exponent vector.
            let mut divisors: Vec<(DivisorOrder, BivariatePolynomial)> = exponent_vectors
                .into_iter()
                .map(|exponents| {
                    let mut psi = BivariatePolynomial::one();
                    for (exponent, factor) in exponents.iter().zip(&factor_polynomials) {
                        psi = psi * factor.pow(*exponent);
                    }
                    let x_degree: usize =
                        exponents.iter().zip(&x_degrees).map(|(e, d)| e * d).sum();
                    (
                        (negative_parent_coefficients(&psi), x_degree, exponents),
                        psi,
                    )
                })
                .collect();
            divisors.sort_by(|a, b| a.0.cmp(&b.0));
            for (_, psi) in divisors {
                if examined >= MAXIMUM_DIVISORS {
                    break 'search;
                }
                examined += 1;
                if power_characteristic_polynomial(&psi, period) != chi {
                    continue;
                }
                let parent_rule = parent_recurrence_from_characteristic(&psi);
                // The shifted start `s_0 + 1` is tried only when `s_0` gives no
                // section identity (otherwise it is the section of the
                // shifted parent `Q_(n+p)`).
                for exceptional_rows in base..=base + 1 {
                    let section = &rows[exceptional_rows..exceptional_rows + order];
                    let transform = if reversal {
                        IndexSectionTransform::Reversal {
                            degree_slope: slope,
                            degree_offset: degree_offset(section, slope),
                        }
                    } else {
                        IndexSectionTransform::Identity
                    };
                    let offsets = match transform {
                        IndexSectionTransform::Identity => None,
                        IndexSectionTransform::Reversal {
                            degree_slope,
                            degree_offset,
                        } => Some((degree_slope, degree_offset)),
                    };
                    let checked =
                        &rows[exceptional_rows..exceptional_rows + order + SELF_CHECK_ROWS];
                    let Some(targets) = transformed_rows(checked, offsets) else {
                        return Err(RecurrenceError::VerificationFailed(
                            "a transformed row past the seeds is not a polynomial".into(),
                        ));
                    };
                    if targets.iter().all(Polynomial::is_zero) {
                        continue;
                    }
                    let Some(solution) =
                        solve_parent_seeds(&parent_rule, period, &targets[..order])?
                    else {
                        continue;
                    };
                    let parent_seeds: Vec<Vec<Q>> = solution.iter().map(polynomial_row).collect();
                    let parent_seeds = normalized_seed_rows(&parent_rule, 0, &parent_seeds)?;
                    let parent_rows = parent_rule
                        .generate_rows_rational(
                            &parent_seeds,
                            0,
                            period * (order + SELF_CHECK_ROWS - 1) + 1,
                        )
                        .map_err(|error| RecurrenceError::VerificationFailed(error.to_string()))?;
                    let agrees =
                        |m: usize| Polynomial::new(parent_rows[period * m].clone()) == targets[m];
                    if !(0..order).all(agrees) {
                        continue;
                    }
                    if !(order..order + SELF_CHECK_ROWS).all(agrees) {
                        return Err(RecurrenceError::VerificationFailed(format!(
                            "an index-section identity of period {period} holds on the seeds but \
                             not on the next rows"
                        )));
                    }
                    let mut lift = SectionLift {
                        period,
                        transform,
                        exceptional_rows,
                        psi: psi.clone(),
                        parent_rule: parent_rule.clone(),
                        parent_seeds,
                        attempts: Vec::new(),
                        selected: None,
                    };
                    if parent_proofs < MAXIMUM_PARENT_PROOFS {
                        parent_proofs += 1;
                        let parent_options = RealRootednessTheoremOptions {
                            root_region: if reversal {
                                None
                            } else {
                                options.root_region.clone()
                            },
                            maximum_start_delay: options.maximum_start_delay,
                        };
                        let (attempts, selected) =
                            prove_parent(&lift.parent_rule, &lift.parent_seeds, &parent_options);
                        lift.attempts = attempts;
                        lift.selected = selected;
                    }
                    if lift.is_proved() {
                        return Ok((Some(lift), examined));
                    }
                    if first_found.is_none() {
                        first_found = Some(lift);
                    }
                    break;
                }
            }
        }
    }
    Ok((first_found, examined))
}

// ---------------------------------------------------------------------------
// Root regions
// ---------------------------------------------------------------------------

/// Whether a parent report is proved with every row from its conclusion's
/// first index on nonzero: a residue-class proof with a zero class (whose
/// rows stay zero) does not count.
fn proves_every_row(report: &RealRootednessTheoremReport) -> bool {
    report.is_proved()
        && report
            .residue_classes
            .as_ref()
            .is_none_or(|certificate| certificate.classes.iter().all(|class| !class.zero))
        && report
            .conclusion
            .as_ref()
            .is_some_and(|conclusion| !conclusion.has_zero_residue_classes())
}

/// The region of the rows `x^D S(1/x)` from a region of the rows `S`
/// (module documentation, "Root region").
fn reversed_region(region: &RealInterval) -> Option<RealInterval> {
    if let Some(upper) = region.upper.as_ref().filter(|u| !u.is_positive()) {
        let (lower, lower_closed) = if upper.is_zero() {
            (None, false)
        } else {
            (Some(upper.recip()), region.upper_closed)
        };
        return Some(RealInterval {
            lower,
            lower_closed,
            upper: Some(Q::zero()),
            upper_closed: true,
        });
    }
    if let Some(lower) = region.lower.as_ref().filter(|l| !l.is_negative()) {
        let (upper, upper_closed) = if lower.is_zero() {
            (None, false)
        } else {
            (Some(lower.recip()), region.lower_closed)
        };
        return Some(RealInterval {
            lower: Some(Q::zero()),
            lower_closed: true,
            upper,
            upper_closed,
        });
    }
    None
}

// ---------------------------------------------------------------------------
// Report
// ---------------------------------------------------------------------------

fn transform_text(transform: IndexSectionTransform, start: usize) -> String {
    match transform {
        IndexSectionTransform::Identity => format!("S_m = P_({start} + m)"),
        IndexSectionTransform::Reversal {
            degree_slope,
            degree_offset,
        } => {
            let exponent = match (degree_slope, degree_offset) {
                (0, e) => e.to_string(),
                (1, 0) => "m".to_string(),
                (d, 0) => format!("{d} m"),
                (1, e) => format!("(m + {e})"),
                (d, e) => format!("({d} m + {e})"),
            };
            format!("S_m(x) = x^{exponent} P_({start} + m)(1/x)")
        }
    }
}

fn polynomial_text(p: &UnivariatePolynomial) -> String {
    format_shift_polynomial(&[p.coeffs().to_vec()], "x")
}

pub(super) fn index_section_report(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    let mut report = RealRootednessTheoremReport::new(
        RealRootednessTheorem::IndexSectionLift,
        REFERENCE,
        rule,
        first_index,
        seeds,
        options,
    );
    let multipliers = match constant_coefficient_multipliers(rule) {
        Ok(multipliers) => multipliers,
        Err(reason) => return Ok(report.not_applicable(reason)),
    };
    let order = multipliers.len();
    if seeds.len() < order {
        return Ok(report.not_applicable(format!(
            "the rule has lag {order} but only {} seed row(s)",
            seeds.len()
        )));
    }
    let (lift, examined) = search_section_lifts(rule, first_index, seeds, &multipliers, options)?;
    let Some(lift) = lift else {
        return Ok(report.not_applicable(format!(
            "no index-section parent: for p = 2, ..., {MAX_INDEX_SECTION_PERIOD} and T = identity \
             or reversal, no monic divisor psi(w) of chi_T(w^p) of degree {order} has \
             Res_w(psi(w), z - w^p) = chi_T(z) and polynomial parent seeds ({examined} \
             divisor(s) examined)"
        )));
    };
    let period = lift.period;
    let section_first = first_index + lift.exceptional_rows;
    let rule_start = first_index + seeds.len();
    let psi_text =
        format_shift_polynomial(&bivariate_to_coefficients(&lift.psi), "x").replace('L', "w");
    let parent_seed_text = lift
        .parent_seeds
        .iter()
        .map(|row| polynomial_text(&Polynomial::new(row.clone())))
        .collect::<Vec<_>>()
        .join(", ");
    let multiplier_text = multipliers
        .iter()
        .enumerate()
        .filter(|(_, a)| !a.is_zero())
        .map(|(lag, a)| format!("({}) P_(n-{})", polynomial_text(a), lag + 1))
        .collect::<Vec<_>>()
        .join(" + ");
    report.normal_form = format!(
        "P_n = {multiplier_text} for n >= {rule_start} (multipliers free of n, order {order}); \
         {}; S_m = Q_({period} m) for every m >= 0, where Q_n is the parent with characteristic \
         polynomial psi(w) = {psi_text} and seeds Q_0, ..., Q_{} = {parent_seed_text}",
        transform_text(lift.transform, section_first),
        order - 1
    );
    report.lean_route = Some("pathPower".into());
    let mut declarations = vec!["LinearRecurrence.ofMatrix_pow_isSolution_subsequence".to_string()];
    if matches!(lift.transform, IndexSectionTransform::Reversal { .. }) {
        declarations.push("RealRooted.Tactic.reflectedRecurrence_isSolution".into());
    }
    if lift.exceptional_rows > 0 {
        declarations.push("RealRooted.Tactic.prependedRecurrence_isSolution".into());
    }
    report.lean_declarations = declarations;

    let mut hypotheses = vec![exact_hypothesis(
        "section_identity",
        &format!(
            "P_({section_first} + m) = T(Q_({period} m)) for every m >= 0 (equal-order \
             recurrence argument)"
        ),
        Ok(format!(
            "Res_w(psi(w), z - w^{period}) equals chi_T(z) exactly (power sums of the roots of \
             psi), so m -> Q_({period} m) satisfies the transformed order-{order} rule for every \
             m >= {order} (Cayley-Hamilton for the {period}-th power of the companion matrix); \
             the seed system was solved over Q(x) with a polynomial solution and Q_({period} m) = \
             S_m for m < {order} (exact); {examined} divisor(s) of chi_T(w^p) examined"
        )),
    )];
    let parent_statement = "the parent family Q_n, n >= 0, is real-rooted by a theorem of the \
                            parent menu";
    let selected_proof = lift.selected.map(|position| &lift.attempts[position]);
    let unresolved = selected_proof.is_none()
        && lift
            .attempts
            .iter()
            .any(|a| a.status == RealRootednessTheoremStatus::Unresolved);
    let parent_outcome = match selected_proof {
        Some(proof) => Ok(format!(
            "{:?} on the parent: {}",
            proof.theorem,
            proof
                .conclusion
                .as_ref()
                .map_or(String::new(), |c| c.statement.clone())
        )),
        None => Err((
            None,
            if lift.attempts.is_empty() {
                "parent open: no theorem of the parent menu applies to the parent rule".to_string()
            } else {
                format!(
                    "parent open: no theorem of the parent menu proves the parent; open: {}",
                    lift.attempts
                        .iter()
                        .map(|a| format!("{:?}: {}", a.theorem, a.open_hypotheses().join("; ")))
                        .collect::<Vec<_>>()
                        .join(" | ")
                )
            },
        )),
    };
    hypotheses.push(match parent_outcome {
        Err((_, evidence)) if unresolved => {
            unresolved_hypothesis("parent_real_rooted", parent_statement, evidence)
        }
        outcome => exact_hypothesis("parent_real_rooted", parent_statement, outcome),
    });

    let mut region = None;
    if let Some(proof) = selected_proof {
        let conclusion = proof
            .conclusion
            .as_ref()
            .expect("a proved report concludes");
        let section_rows_checked = conclusion.first_index.div_ceil(period);
        let checked_count = lift.exceptional_rows + section_rows_checked;
        region = proved_region(proof).and_then(|parent_region| match lift.transform {
            IndexSectionTransform::Identity => Some(parent_region),
            IndexSectionTransform::Reversal { .. } => reversed_region(&parent_region),
        });
        if checked_count > 0 {
            let rows = rule
                .generate_rows_rational(seeds, first_index, checked_count)
                .map_err(|error| RecurrenceError::VerificationFailed(error.to_string()))?;
            let failing = rows.iter().position(|row| {
                let row = vec_poly::trim(row);
                vec_poly::is_zero(&row) || !is_real_rooted_rational(&row)
            });
            if region
                .as_ref()
                .is_some_and(|r| rows.iter().any(|row| !has_real_roots_in_closed_region(row, r)))
            {
                region = None;
            }
            let last = first_index + checked_count - 1;
            hypotheses.push(exact_hypothesis(
                "early_rows",
                &format!(
                    "the rows P_{first_index}, ..., P_{last} (exceptional rows and section rows \
                     before the parent's first proved index {}) are nonzero and real-rooted",
                    conclusion.first_index
                ),
                match failing {
                    None => Ok(format!(
                        "checked exactly for P_{first_index}, ..., P_{last}"
                    )),
                    Some(offset) => Err((
                        Some(first_index + offset),
                        format!(
                            "P_{} is zero or not real-rooted (exact)",
                            first_index + offset
                        ),
                    )),
                },
            ));
        }
    }
    if let Some(requested) = &options.root_region {
        let statement = format!(
            "the requested region {} contains the transferred root region",
            open_root_interval(requested)
        );
        hypotheses.push(match &region {
            Some(certified) if requested_region_contains(requested, certified) => {
                exact_hypothesis("requested_root_region", &statement, Ok("exact".into()))
            }
            Some(certified) => unresolved_hypothesis(
                "requested_root_region",
                &statement,
                format!(
                    "the lift certifies the region {} only",
                    RootInterval::from_real(certified)
                ),
            ),
            None => unresolved_hypothesis(
                "requested_root_region",
                &statement,
                "the lift certifies no index-independent root region".into(),
            ),
        });
    }
    report.hypotheses = hypotheses;
    let parent_record =
        RecurrenceJson::from_recurrence_rational(&lift.parent_rule, 0, &lift.parent_seeds, None);
    report.index_section = Some(IndexSectionCertificate {
        period,
        transform: lift.transform,
        exceptional_rows: lift.exceptional_rows,
        parent: parent_record,
        divisors_examined: examined,
        parent_maximum_start_delay: options.maximum_start_delay,
        parent_attempts: lift.attempts.clone(),
        selected: lift.selected,
    });
    let mut statement = format!(
        "every P_n, n >= {first_index}, is nonzero and real-rooted: P_({section_first} + m) = \
         T(Q_({period} m)) for every m >= 0 with {} and the parent psi(w) = {psi_text}",
        transform_text(lift.transform, section_first)
    );
    if let Some(proof) = selected_proof {
        statement.push_str(&format!(
            ", whose rows Q_n, n >= {}, are nonzero and real-rooted by {:?}",
            proof.conclusion.as_ref().map_or(0, |c| c.first_index),
            proof.theorem
        ));
    }
    if let Some(r) = &region {
        statement.push_str(&format!(
            " with every root in {}",
            RootInterval::from_real(r)
        ));
    }
    statement.push_str("; no interlacing of consecutive rows is concluded");
    Ok(report.finish(RealRootednessTheoremConclusion {
        first_index,
        interlacing_from: None,
        interlacing: None,
        simple_roots: false,
        root_region: region.as_ref().map(RootInterval::from_real),
        zero_residue_classes: None,
        statement,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(n: i64) -> Q {
        Q::from_integer(n.into())
    }

    fn x_polynomial(coefficients: &[i64]) -> UnivariatePolynomial {
        Polynomial::new(coefficients.iter().map(|&c| q(c)).collect())
    }

    #[test]
    fn power_characteristic_polynomial_is_the_resultant() {
        // psi = w^2 - w - x has roots a, b with a + b = 1, a b = -x;
        // a^2 + b^2 = 1 + 2x and a^2 b^2 = x^2, so the squares have
        // z^2 - (1 + 2x) z + x^2.
        let psi: BivariatePolynomial = Polynomial::new(vec![
            x_polynomial(&[0, -1]),
            x_polynomial(&[-1]),
            x_polynomial(&[1]),
        ]);
        let expected: BivariatePolynomial = Polynomial::new(vec![
            x_polynomial(&[0, 0, 1]),
            x_polynomial(&[-1, -2]),
            x_polynomial(&[1]),
        ]);
        assert_eq!(power_characteristic_polynomial(&psi, 2), expected);
        // Cubes: a^3 + b^3 = 1 + 3x, a^3 b^3 = -x^3.
        let expected: BivariatePolynomial = Polynomial::new(vec![
            x_polynomial(&[0, 0, 0, -1]),
            x_polynomial(&[-1, -3]),
            x_polynomial(&[1]),
        ]);
        assert_eq!(power_characteristic_polynomial(&psi, 3), expected);
        // psi divides chi(w^p) for its own power characteristic polynomial.
        let chi = power_characteristic_polynomial(&psi, 3);
        let (_, remainder) = divide_monic_in_shift(&substitute_shift_power(&chi, 3), &psi).unwrap();
        assert!(remainder.is_zero());
    }

    #[test]
    fn reversal_multipliers_follow_the_degree_slope() {
        // a_1 = 3 + x, a_2 = -3, a_3 = 1: slope 1, b_k = x^k a_k(1/x).
        let a = vec![
            x_polynomial(&[3, 1]),
            x_polynomial(&[-3]),
            x_polynomial(&[1]),
        ];
        assert_eq!(degree_slope(&a), 1);
        let b = reversed_multipliers(&a, 1).unwrap();
        assert_eq!(
            b,
            vec![
                x_polynomial(&[1, 3]),
                x_polynomial(&[0, 0, -3]),
                x_polynomial(&[0, 0, 0, 1])
            ]
        );
        assert!(reversed_multipliers(&[x_polynomial(&[0, 0, 1])], 1).is_none());
    }

    #[test]
    fn perfect_powers_are_recognized_exactly() {
        // 4x^4 = (2x^2)^2, not a cube; (1 + x)^3 / 8 = ((1 + x) / 2)^3.
        assert!(is_perfect_power(&x_polynomial(&[0, 0, 0, 0, 4]), 2));
        assert!(!is_perfect_power(&x_polynomial(&[0, 0, 0, 0, 4]), 3));
        let cube = Polynomial::new(vec![q(1), q(3), q(3), q(1)]).scale(&Q::new(1.into(), 8.into()));
        assert!(is_perfect_power(&cube, 3));
        assert!(!is_perfect_power(&x_polynomial(&[1, 3, 3, 2]), 3));
        assert!(is_perfect_power(&x_polynomial(&[-1]), 3));
        assert!(!is_perfect_power(&x_polynomial(&[-1]), 2));
        assert!(!is_perfect_power(&x_polynomial(&[2]), 2));
    }
}
