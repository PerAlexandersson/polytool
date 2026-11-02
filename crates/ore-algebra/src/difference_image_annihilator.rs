//! Runnable annihilators of images of sequences under maps with the
//! variable shift `E`, found as common left multiples, with identity
//! certificates.
//!
//! Let `F` be annihilated by `L` (`(L F)_n = 0` from some base on), a
//! [`DifferentialDifferenceOperator`] whose leader is `g(n) S^k` (free of
//! `x`, of derivatives and of `E`), and let
//! `G = step_last(... step_first(F))` for a list of
//! [`DifferenceRowMapStep`]s.  An annihilator `M` of `G` is **runnable**
//! when its top shift column is a single term `m_R(n, x) S^R` without
//! derivatives and without `E`: then `m_R(n, x) G_(n+R) = -(the other
//! terms)` determines `G_(n+R)` from the earlier rows wherever
//! `m_R(n, .)` is not the zero polynomial.
//!
//! # The certificate
//!
//! [`DifferenceImageAnnihilatorCertificate`] records `M`, `L`, the steps,
//! the [`DifferencePullbackBounds`] and the identity
//!
//! ```text
//! w(n) P = H o L,   (P, b) = pull_back_difference(M, steps, bounds),
//! ```
//!
//! with `w` a nonzero polynomial in `n` only, so the identity says that `M`
//! composed with the maps is a common left multiple of the maps and `L`.
//! Applied to `F` at a base `n` from which `L` holds,
//! `w(n) (P F)_n = sum h(n, x) (D^d E^e (L F)_(n+nu)) = 0`, hence
//! `(P F)_n = 0` where `w(n) != 0`, hence `(M G)_n = 0` where moreover
//! `b(n) != 0` (the contract of [`crate::pull_back_difference`]).  The
//! finitely many other bases are the caller's to check on rows.
//! [`DifferenceImageAnnihilatorCertificate::verify`] recomputes the pullback
//! and both sides of the identity exactly.
//!
//! # How `M` is found
//!
//! When every step has a push-forward ([`crate::push_forward_difference`])
//! and the pushed-forward operator is runnable, it is normalised and
//! certified.  Otherwise `M` is searched with an ansatz
//! `M = m_R S^R + sum_(i < R, d <= J, 0 <= e <= K) m_(i,d,e) D^d E^e S^i`
//! with coefficients of bounded degrees in `n` and `x`
//! ([`DifferenceAnsatzShape`]).  For fixed bounds the pullback and the
//! right pseudo-division by `L` with fixed levels
//! ([`crate::difference_right_pseudo_division`]) are linear in `M` over
//! polynomials in `n`, and the remainder vanishes exactly when the pullback
//! lies in the left ideal of `L` over `Q(n, x)`; so the condition is a
//! homogeneous linear system, solved exactly
//! ([`crate::exact_kernel_vectors_by_group`]) with the top-shift unknowns
//! last.  Shapes are tried by increasing order plus derivative order plus
//! `E`-order, and for each by increasing coefficient degrees
//! (`0, 1, 2, 4, ...`), skipping shapes dominated by one that already
//! failed, within a [`DifferenceImageAnnihilatorBudget`].  When nothing is
//! found the search reports exactly that; every returned annihilator
//! carries a certificate that has been verified.
//!
//! A found `M` is normalised as in [`crate::normalize_annihilator`]: divided
//! by the largest power of `x` and the largest polynomial in `n` dividing
//! all its coefficients (left factors), and scaled to coprime integer
//! coefficients with a positive leading coefficient; the unnormalised
//! operator is kept when the normalised one does not certify.

use crate::bivariate::BivariatePolynomial;
use crate::difference_division::{difference_index_only_leader, difference_right_pseudo_division};
use crate::difference_row_map::{
    pull_back_difference, push_forward_difference, DifferencePullbackBounds, DifferenceRowMapError,
    DifferenceRowMapStep,
};
use crate::differential_difference::DifferentialDifferenceOperator;
use crate::kernel::exact_kernel_vectors_by_group;
use crate::univariate as poly;
use crate::Rational;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::collections::{BTreeMap, HashMap};
use std::fmt;

/// The bounds of one ansatz
/// `M = m_R S^R + sum_(i < R, d <= J, 0 <= e <= K) m_(i,d,e) D^d E^e S^i`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct DifferenceAnsatzShape {
    /// The order `R` (top shift).
    pub order: usize,
    /// The largest derivative order `J` below the top shift.
    pub derivative_order: usize,
    /// The largest power `K` of `E` below the top shift.
    pub variable_shift_order: usize,
    /// The largest degree in `n` of a coefficient.
    pub index_degree: usize,
    /// The largest degree in `x` of a coefficient.
    pub variable_degree: usize,
}

impl DifferenceAnsatzShape {
    /// The number of unknown rational coefficients.
    pub fn unknowns(&self) -> usize {
        let monomials = (self.index_degree + 1) * (self.variable_degree + 1);
        let lower = (self.derivative_order + 1) * (self.variable_shift_order + 1);
        monomials * (self.order * lower + 1)
    }

    fn dominated_by(&self, other: &Self) -> bool {
        self.order <= other.order
            && self.derivative_order <= other.derivative_order
            && self.variable_shift_order <= other.variable_shift_order
            && self.index_degree <= other.index_degree
            && self.variable_degree <= other.variable_degree
    }

    /// The pullback bounds of the ansatz.
    pub fn bounds(&self) -> DifferencePullbackBounds {
        DifferencePullbackBounds {
            top_shift: self.order,
            derivative_order: self.derivative_order,
            lowest_variable_shift: 0,
            highest_variable_shift: self.variable_shift_order as i64,
            variable_degree: self.variable_degree,
        }
    }
}

/// The budget of [`find_difference_image_annihilator`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DifferenceImageAnnihilatorBudget {
    /// Largest order above the parent's order `k`.  Default 2.
    pub max_order_increase: usize,
    /// Largest derivative order of the ansatz.  Default 2.
    pub max_derivative_order: usize,
    /// Largest power of `E` of the ansatz.  Default 2.
    pub max_variable_shift_order: usize,
    /// Largest coefficient degree in `n`.  Default 4.
    pub max_index_degree: usize,
    /// Largest coefficient degree in `x`.  Default 4.
    pub max_variable_degree: usize,
    /// Largest number of unknowns of one ansatz.  Default 400.
    pub max_unknowns: usize,
    /// Largest number of ansatz shapes solved.  Default 60.
    pub max_shapes: usize,
}

impl Default for DifferenceImageAnnihilatorBudget {
    fn default() -> Self {
        Self {
            max_order_increase: 2,
            max_derivative_order: 2,
            max_variable_shift_order: 2,
            max_index_degree: 4,
            max_variable_degree: 4,
            max_unknowns: 400,
            max_shapes: 60,
        }
    }
}

/// Why [`find_difference_image_annihilator`] or a certificate check failed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DifferenceImageAnnihilatorError {
    /// The parent annihilator's top shift column is not a single term
    /// `g(n) S^k` free of `x`, of derivatives and of `E`.
    ParentLeaderNotIndexOnly,
    /// A step is invalid or does not apply.
    Step(DifferenceRowMapError),
    /// The recorded annihilator is not runnable.
    NotRunnable,
    /// The recorded bounds are below the annihilator's shape.
    BoundsTooSmall,
    /// The recorded multiplier is zero or involves `x`.
    InvalidMultiplier,
    /// The identity `w P = H o L` does not hold.
    IdentityFails,
}

impl fmt::Display for DifferenceImageAnnihilatorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ParentLeaderNotIndexOnly => f.write_str(
                "the parent annihilator's leader must be a single term g(n) S^k free of x, of \
                 derivatives and of E",
            ),
            Self::Step(error) => write!(f, "invalid step: {error}"),
            Self::NotRunnable => f.write_str(
                "the annihilator is not runnable (its top shift column is not a single term \
                 without derivatives and without E)",
            ),
            Self::BoundsTooSmall => f.write_str("the pullback bounds are below the shape"),
            Self::InvalidMultiplier => {
                f.write_str("the multiplier must be a nonzero polynomial in n only")
            }
            Self::IdentityFails => f.write_str("the identity w P = H o L does not hold"),
        }
    }
}

impl std::error::Error for DifferenceImageAnnihilatorError {}

impl From<DifferenceRowMapError> for DifferenceImageAnnihilatorError {
    fn from(error: DifferenceRowMapError) -> Self {
        Self::Step(error)
    }
}

/// Whether `op` is runnable: nonzero, with a top shift column that is a
/// single term without derivatives and without `E`.
pub fn is_difference_runnable(op: &DifferentialDifferenceOperator) -> bool {
    if op.is_zero() {
        return false;
    }
    let column = op.column(op.max_shift());
    column.term_count() == 1 && column.0.contains_key(&(0, 0, 0))
}

/// The certificate of a runnable image annihilator (module documentation).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DifferenceImageAnnihilatorCertificate {
    /// `M`, runnable, acting on the image `G`.
    pub annihilator: DifferentialDifferenceOperator,
    /// `L`, acting on the source `F`.
    pub parent_annihilator: DifferentialDifferenceOperator,
    /// The steps `F -> G`, applied first to last.
    pub steps: Vec<DifferenceRowMapStep>,
    /// The bounds of the pullback of `M`.
    pub bounds: DifferencePullbackBounds,
    /// `w(n)`, a nonzero polynomial in `n` only.
    pub multiplier: BivariatePolynomial,
    /// `H`.
    pub quotient: DifferentialDifferenceOperator,
}

impl DifferenceImageAnnihilatorCertificate {
    /// Recompute the pullback and check `w P = H o L` exactly, with `M`
    /// runnable, `w` a nonzero polynomial in `n` only and the parent leader
    /// free of `x`, of derivatives and of `E`.  Returns the base factor
    /// `b(n)`.
    pub fn verify(&self) -> Result<Vec<Rational>, DifferenceImageAnnihilatorError> {
        if !is_difference_runnable(&self.annihilator) {
            return Err(DifferenceImageAnnihilatorError::NotRunnable);
        }
        difference_index_only_leader(&self.parent_annihilator)
            .ok_or(DifferenceImageAnnihilatorError::ParentLeaderNotIndexOnly)?;
        if self.multiplier.is_zero() || !self.multiplier.is_variable_free() {
            return Err(DifferenceImageAnnihilatorError::InvalidMultiplier);
        }
        if !self.bounds.contains(&self.annihilator) {
            return Err(DifferenceImageAnnihilatorError::BoundsTooSmall);
        }
        let (pulled, base_factor) =
            pull_back_difference(&self.annihilator, &self.steps, self.bounds)?;
        let left = pulled.left_multiply(&self.multiplier);
        let right = self.quotient.compose(&self.parent_annihilator);
        if left != right {
            return Err(DifferenceImageAnnihilatorError::IdentityFails);
        }
        Ok(base_factor)
    }
}

/// How an annihilator was found.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DifferenceImageAnnihilatorMethod {
    /// Pushed forward through the steps.
    PushForward,
    /// A common left multiple from the ansatz of this shape.
    LeftMultiple(DifferenceAnsatzShape),
}

/// The outcome of [`find_difference_image_annihilator`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DifferenceImageAnnihilatorSearch {
    /// A runnable annihilator with a verified certificate.
    Found {
        /// The certificate.
        certificate: Box<DifferenceImageAnnihilatorCertificate>,
        /// How it was found.
        method: DifferenceImageAnnihilatorMethod,
    },
    /// No runnable annihilator within the budget (no statement about larger
    /// shapes).
    NotFoundWithinBudget {
        /// The number of ansatz shapes solved.
        shapes_tried: usize,
        /// The last shape solved.
        largest_shape: Option<DifferenceAnsatzShape>,
    },
}

/// Divide `op` by the largest left factors `x^e` and `c(n)` of all its
/// coefficients and scale it to coprime integer coefficients with a
/// positive leading coefficient.
pub fn normalize_difference_annihilator(
    op: &DifferentialDifferenceOperator,
) -> DifferentialDifferenceOperator {
    if op.is_zero() {
        return op.clone();
    }
    let valuation =
        op.0.values()
            .filter_map(BivariatePolynomial::variable_valuation)
            .min()
            .unwrap_or(0);
    let mut stripped = DifferentialDifferenceOperator::default();
    for (key, c) in &op.0 {
        stripped.add_term(*key, c.divide_by_variable_power(valuation));
    }
    let mut content = vec![Rational::zero()];
    for c in stripped.0.values() {
        for j in 0..=c.variable_degree().unwrap_or(0) {
            content = poly::gcd(&content, &c.variable_coefficient(j));
        }
    }
    let mut out = DifferentialDifferenceOperator::default();
    for (key, c) in &stripped.0 {
        let width = c.variable_degree().map_or(0, |d| d + 1);
        let columns: Vec<Vec<Rational>> = (0..width)
            .map(|j| poly::exact_quotient(&c.variable_coefficient(j), &content))
            .collect();
        let height = columns.iter().map(Vec::len).max().unwrap_or(0);
        let rows = (0..height)
            .map(|i| {
                columns
                    .iter()
                    .map(|column| column.get(i).cloned().unwrap_or_else(Rational::zero))
                    .collect()
            })
            .collect();
        out.add_term(*key, BivariatePolynomial::new(rows));
    }
    let values: Vec<&Rational> = out
        .0
        .values()
        .flat_map(|c| c.rows().iter().flatten())
        .collect();
    let denominators = values
        .iter()
        .fold(BigInt::one(), |acc, value| acc.lcm(value.denom()));
    let numerators = values
        .iter()
        .filter(|value| !value.is_zero())
        .fold(BigInt::zero(), |acc, value| {
            acc.gcd(&(value.numer() * (&denominators / value.denom())))
        });
    let mut scale = Rational::new(denominators, numerators);
    if out.leading().is_some_and(|lead| lead.is_negative()) {
        scale = -scale;
    }
    out.scale(&scale)
}

/// Certify a runnable candidate `M` (module documentation): pull back with
/// `bounds`, pseudo-divide by `L` and require a zero remainder.  `Ok(None)`
/// when `M` is not runnable or the remainder is not zero.
pub fn certify_difference_image_annihilator(
    candidate: &DifferentialDifferenceOperator,
    parent: &DifferentialDifferenceOperator,
    steps: &[DifferenceRowMapStep],
    bounds: DifferencePullbackBounds,
) -> Result<Option<DifferenceImageAnnihilatorCertificate>, DifferenceImageAnnihilatorError> {
    if !is_difference_runnable(candidate) {
        return Ok(None);
    }
    difference_index_only_leader(parent)
        .ok_or(DifferenceImageAnnihilatorError::ParentLeaderNotIndexOnly)?;
    if !bounds.contains(candidate) {
        return Err(DifferenceImageAnnihilatorError::BoundsTooSmall);
    }
    let (pulled, _) = pull_back_difference(candidate, steps, bounds)?;
    let Some(division) = difference_right_pseudo_division(&pulled, parent, bounds.top_shift) else {
        return Ok(None);
    };
    if !division.remainder.is_zero() {
        return Ok(None);
    }
    let certificate = DifferenceImageAnnihilatorCertificate {
        annihilator: candidate.clone(),
        parent_annihilator: parent.clone(),
        steps: steps.to_vec(),
        bounds,
        multiplier: division.multiplier,
        quotient: division.quotient,
    };
    // Defensive: the returned certificate always replays.
    Ok(certificate.verify().is_ok().then_some(certificate))
}

/// Normalise and certify, falling back to the unnormalised operator.
fn certify_normalized(
    candidate: &DifferentialDifferenceOperator,
    parent: &DifferentialDifferenceOperator,
    steps: &[DifferenceRowMapStep],
    fallback_bounds: DifferencePullbackBounds,
) -> Result<Option<DifferenceImageAnnihilatorCertificate>, DifferenceImageAnnihilatorError> {
    let normalized = normalize_difference_annihilator(candidate);
    for (operator, bounds) in [
        (&normalized, DifferencePullbackBounds::of(&normalized)),
        (&normalized, fallback_bounds),
        (candidate, DifferencePullbackBounds::of(candidate)),
        (candidate, fallback_bounds),
    ] {
        if !bounds.contains(operator) {
            continue;
        }
        // A step may not apply to every shape (a basis change needs an
        // operator free of D or of E); such a pair is skipped.
        match certify_difference_image_annihilator(operator, parent, steps, bounds) {
            Ok(Some(certificate)) => return Ok(Some(certificate)),
            Ok(None) | Err(DifferenceImageAnnihilatorError::Step(_)) => {}
            Err(error) => return Err(error),
        }
    }
    Ok(None)
}

/// The shapes within the budget, by increasing order plus derivative order
/// plus `E`-order (then order), and for each the coefficient degrees
/// `0, 1, 2, 4, ...`.
fn shapes(order: usize, budget: &DifferenceImageAnnihilatorBudget) -> Vec<DifferenceAnsatzShape> {
    let max_degree = budget.max_index_degree.max(budget.max_variable_degree);
    let mut degrees = vec![0, 1];
    while *degrees.last().expect("nonempty") < max_degree {
        let next = (degrees.last().expect("nonempty") * 2).min(max_degree);
        degrees.push(next);
    }
    let mut triples = Vec::new();
    for r in order.max(1)..=order + budget.max_order_increase {
        for d in 0..=budget.max_derivative_order {
            for k in 0..=budget.max_variable_shift_order {
                triples.push((r, d, k));
            }
        }
    }
    triples.sort_by_key(|(r, d, k)| (r - order + d + k, *r, *k));
    let mut out = Vec::new();
    for (r, d, k) in triples {
        for g in &degrees {
            let shape = DifferenceAnsatzShape {
                order: r,
                derivative_order: d,
                variable_shift_order: k,
                index_degree: (*g).min(budget.max_index_degree),
                variable_degree: (*g).min(budget.max_variable_degree),
            };
            if shape.unknowns() <= budget.max_unknowns && !out.contains(&shape) {
                out.push(shape);
            }
        }
    }
    out
}

/// One unknown: the monomial `n^a x^b D^d E^e S^i`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct Unknown {
    shift: usize,
    derivative: usize,
    variable_shift: i64,
    index_power: usize,
    variable_power: usize,
}

/// The outcome of one ansatz.
enum ShapeOutcome {
    /// A runnable candidate with zero remainder.
    Candidate(DifferentialDifferenceOperator),
    /// No kernel vector with a nonzero top-shift entry.
    NoCandidate,
    /// A step does not apply to the operators of this shape (a basis
    /// change needs an operator free of `D` or of `E`).
    Inapplicable,
}

/// Solve one ansatz.
fn solve_shape(
    shape: &DifferenceAnsatzShape,
    parent: &DifferentialDifferenceOperator,
    steps: &[DifferenceRowMapStep],
) -> Result<ShapeOutcome, DifferenceImageAnnihilatorError> {
    let mut unknowns = Vec::new();
    let mut groups = Vec::new();
    for i in 0..shape.order {
        for d in 0..=shape.derivative_order {
            for e in 0..=shape.variable_shift_order as i64 {
                for a in 0..=shape.index_degree {
                    for b in 0..=shape.variable_degree {
                        unknowns.push(Unknown {
                            shift: i,
                            derivative: d,
                            variable_shift: e,
                            index_power: a,
                            variable_power: b,
                        });
                        groups.push(0);
                    }
                }
            }
        }
    }
    for a in 0..=shape.index_degree {
        for b in 0..=shape.variable_degree {
            unknowns.push(Unknown {
                shift: shape.order,
                derivative: 0,
                variable_shift: 0,
                index_power: a,
                variable_power: b,
            });
            groups.push(1);
        }
    }
    let bounds = shape.bounds();
    // Remainders of the index-free monomials; that of n^a times a monomial
    // is n^a times it (left multiplication by n^a commutes with every
    // pullback step and with the division).
    type MonomialKey = (usize, usize, i64, usize);
    let mut remainders: HashMap<MonomialKey, DifferentialDifferenceOperator> = HashMap::new();
    for unknown in &unknowns {
        let key = (
            unknown.shift,
            unknown.derivative,
            unknown.variable_shift,
            unknown.variable_power,
        );
        if remainders.contains_key(&key) {
            continue;
        }
        let monomial = DifferentialDifferenceOperator::single(
            unknown.shift,
            unknown.derivative,
            unknown.variable_shift,
            BivariatePolynomial::monomial(0, unknown.variable_power, Rational::one()),
        );
        // The steps are validated, so an error means that a step does not
        // apply to this shape.
        let Ok((pulled, _)) = pull_back_difference(&monomial, steps, bounds) else {
            return Ok(ShapeOutcome::Inapplicable);
        };
        let division = difference_right_pseudo_division(&pulled, parent, shape.order)
            .ok_or(DifferenceImageAnnihilatorError::ParentLeaderNotIndexOnly)?;
        remainders.insert(key, division.remainder);
    }
    type EquationKey = (usize, usize, i64, usize, usize);
    let mut equations: BTreeMap<EquationKey, BTreeMap<usize, Rational>> = BTreeMap::new();
    for (column, unknown) in unknowns.iter().enumerate() {
        let remainder = &remainders[&(
            unknown.shift,
            unknown.derivative,
            unknown.variable_shift,
            unknown.variable_power,
        )];
        let factor = BivariatePolynomial::monomial(unknown.index_power, 0, Rational::one());
        for ((i, d, e), c) in &remainder.left_multiply(&factor).0 {
            for (p, row) in c.rows().iter().enumerate() {
                for (q, value) in row.iter().enumerate() {
                    if !value.is_zero() {
                        equations
                            .entry((*i, *d, *e, p, q))
                            .or_default()
                            .insert(column, value.clone());
                    }
                }
            }
        }
    }
    let rows: Vec<Vec<BigInt>> = equations
        .values()
        .map(|entries| {
            let denominator = entries
                .values()
                .fold(BigInt::one(), |acc, value| acc.lcm(value.denom()));
            let mut row = vec![BigInt::zero(); unknowns.len()];
            for (column, value) in entries {
                row[*column] = value.numer() * (&denominator / value.denom());
            }
            row
        })
        .collect();
    let vectors = exact_kernel_vectors_by_group(&rows, &groups, 2);
    let Some(vector) = vectors[1].as_ref() else {
        return Ok(ShapeOutcome::NoCandidate);
    };
    let mut candidate = DifferentialDifferenceOperator::default();
    for (unknown, value) in unknowns.iter().zip(vector) {
        if !value.is_zero() {
            candidate.add_term(
                (unknown.shift, unknown.derivative, unknown.variable_shift),
                BivariatePolynomial::monomial(
                    unknown.index_power,
                    unknown.variable_power,
                    value.clone(),
                ),
            );
        }
    }
    Ok(ShapeOutcome::Candidate(candidate))
}

/// A runnable annihilator of the image of a sequence annihilated by
/// `parent` under `steps`, with a verified certificate, or the report that
/// none was found within `budget` (module documentation).
///
/// # Errors
///
/// [`DifferenceImageAnnihilatorError::ParentLeaderNotIndexOnly`] for a
/// parent without a leader `g(n) S^k`, and
/// [`DifferenceImageAnnihilatorError::Step`] for an invalid step.
pub fn find_difference_image_annihilator(
    parent: &DifferentialDifferenceOperator,
    steps: &[DifferenceRowMapStep],
    budget: &DifferenceImageAnnihilatorBudget,
) -> Result<DifferenceImageAnnihilatorSearch, DifferenceImageAnnihilatorError> {
    let (order, _) = difference_index_only_leader(parent)
        .ok_or(DifferenceImageAnnihilatorError::ParentLeaderNotIndexOnly)?;
    for step in steps {
        step.validate()?;
    }
    let pushed = steps.iter().try_fold(parent.clone(), |op, step| {
        push_forward_difference(&op, step)
    });
    if let Some(pushed) = pushed.filter(is_difference_runnable) {
        let bounds = DifferencePullbackBounds::of(&pushed);
        if let Some(certificate) = certify_normalized(&pushed, parent, steps, bounds)? {
            return Ok(DifferenceImageAnnihilatorSearch::Found {
                certificate: Box::new(certificate),
                method: DifferenceImageAnnihilatorMethod::PushForward,
            });
        }
    }
    let mut failed: Vec<DifferenceAnsatzShape> = Vec::new();
    let mut tried = 0;
    let mut largest = None;
    for shape in shapes(order, budget) {
        if tried >= budget.max_shapes {
            break;
        }
        if failed.iter().any(|f| shape.dominated_by(f)) {
            continue;
        }
        let outcome = solve_shape(&shape, parent, steps)?;
        if matches!(outcome, ShapeOutcome::Inapplicable) {
            continue;
        }
        tried += 1;
        largest = Some(shape);
        match outcome {
            ShapeOutcome::Candidate(candidate) => {
                if let Some(certificate) =
                    certify_normalized(&candidate, parent, steps, shape.bounds())?
                {
                    return Ok(DifferenceImageAnnihilatorSearch::Found {
                        certificate: Box::new(certificate),
                        method: DifferenceImageAnnihilatorMethod::LeftMultiple(shape),
                    });
                }
                failed.push(shape);
            }
            ShapeOutcome::NoCandidate | ShapeOutcome::Inapplicable => failed.push(shape),
        }
    }
    Ok(DifferenceImageAnnihilatorSearch::NotFoundWithinBudget {
        shapes_tried: tried,
        largest_shape: largest,
    })
}
