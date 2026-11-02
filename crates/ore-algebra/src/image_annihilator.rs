//! Runnable annihilators of images of sequences, with identity
//! certificates.
//!
//! Let `F` be annihilated by `L` (`(L F)_n = 0` from some base on) whose
//! leader is `g(n) S^k` (free of `x` and of derivatives), and let
//! `G = step_last(... step_first(F))` for a list of [`RowMapStep`]s.  An
//! annihilator `M` of `G` is **runnable** when its top shift column is a
//! single term `m_R(n, x) S^R` without derivatives: then
//! `m_R(n, x) G_(n+R) = -(the other terms)` determines `G_(n+R)` from the
//! earlier rows wherever `m_R(n, .)` is not the zero polynomial.
//!
//! # The certificate
//!
//! [`ImageAnnihilatorCertificate`] records `M`, `L`, the steps, the
//! [`PullbackBounds`], optionally a row annihilator `A` of `F` (a
//! shift-free operator with `(A F)_n = 0` from some index on, proved by an
//! [`crate::AnnihilatorSystemCertificate`]) and the identity
//!
//! ```text
//! w(n, x) E = H o L + Z o A,   (E, b) = pull_back(M, steps, bounds),
//! ```
//!
//! with `w` nonzero (`Z = 0` and `w` a polynomial in `n` only without
//! `A`).  So `w(n, x) (E F)_n = sum h(n, x) D^d (L F)_(n+nu) +
//! sum z(n, x) D^d (A F)_(n+i) = 0` at every base `n` from which `L` and `A`
//! hold, hence `(E F)_n = 0` where `w(n, .)` is not the zero polynomial
//! (`Q[x]` has no zero divisors), hence `(M G)_n = 0` where moreover
//! `b(n) != 0` (the contract of [`crate::pull_back`]).  The finitely many
//! other bases are the caller's to check on rows.
//! [`ImageAnnihilatorCertificate::verify`] recomputes the pullback and both
//! sides of the identity exactly.
//!
//! The identity comes from the right pseudo-division `w_L E = H o L + R`
//! ([`crate::right_pseudo_division`]) and, when `R != 0` and `A` is given,
//! the reduction of `R` by `A` ([`crate::reduce_by_row_annihilator`]).
//!
//! # How `M` is found
//!
//! Without a Weyl-operator step, `L` is pushed forward through the steps
//! ([`crate::push_forward`]); when the result is runnable it is
//! normalised and certified.  Otherwise (and always for a Weyl operator)
//! `M` is searched as a common left multiple: an ansatz
//! `M = m_R S^R + sum_(i < R, d <= J) m_(i,d) D^d S^i` with coefficients
//! of bounded degree in `n` and `x` ([`AnsatzShape`]).  For fixed bounds the
//! pullback and the reductions with fixed levels are linear in `M`, so
//! the condition that the remainder vanishes (which says that the pullback
//! lies in the left ideal of `L` over `Q(n, x)`, or of `L` and `A`) is a
//! homogeneous linear system.
//! Its kernel is computed exactly ([`crate::exact_kernel_vectors_by_group`])
//! with the top-shift unknowns last, and a kernel vector with a nonzero
//! top-shift entry is taken.  Shapes are tried by increasing order plus
//! derivative order and, for each, increasing coefficient degrees
//! (`0, 1, 2, 4, ...`), skipping shapes dominated by one that already failed (for a
//! composition `M o T`, `S^r M` solves every larger shape when `M` solves a
//! smaller one; elsewhere the skip only narrows the search), within an
//! [`ImageAnnihilatorBudget`].  When nothing is found the search reports
//! exactly that ([`ImageAnnihilatorSearch::NotFoundWithinBudget`]); nothing
//! is guessed, and every returned annihilator carries a certificate that
//! has been verified.
//!
//! # Normalisation
//!
//! A found `M` is divided by the largest power of `x` and the largest
//! polynomial in `n` dividing all its coefficients (left factors, so the
//! quotient still annihilates `G` wherever the factor does not vanish, and
//! it is certified on its own), and scaled to integer coefficients with no
//! common factor and a positive leading coefficient.  When the normalised
//! operator does not certify, the unnormalised one is kept.

use crate::bivariate::BivariatePolynomial;
use crate::division::{index_only_leader, right_pseudo_division};
use crate::kernel::exact_kernel_vectors_by_group;
use crate::operator::OreOperator;
use crate::row_annihilator::reduce_by_row_annihilator;
use crate::row_map::{pull_back, push_forward, PullbackBounds, RowMapError, RowMapStep};
use crate::univariate as poly;
use crate::Rational;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::collections::{BTreeMap, HashMap};
use std::fmt;

/// The bounds of one ansatz `M = m_R S^R + sum_(i < R, d <= J) m_(i,d) D^d S^i`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct AnsatzShape {
    /// The order `R` (top shift).
    pub order: usize,
    /// The largest derivative order `J` below the top shift.
    pub derivative_order: usize,
    /// The largest degree in `n` of a coefficient.
    pub index_degree: usize,
    /// The largest degree in `x` of a coefficient.
    pub variable_degree: usize,
}

impl AnsatzShape {
    /// The number of unknown rational coefficients.
    pub fn unknowns(&self) -> usize {
        let monomials = (self.index_degree + 1) * (self.variable_degree + 1);
        monomials * (self.order * (self.derivative_order + 1) + 1)
    }

    pub(crate) fn dominated_by(&self, other: &Self) -> bool {
        self.order <= other.order
            && self.derivative_order <= other.derivative_order
            && self.index_degree <= other.index_degree
            && self.variable_degree <= other.variable_degree
    }

    fn bounds(&self) -> PullbackBounds {
        PullbackBounds {
            top_shift: self.order,
            derivative_order: self.derivative_order,
            variable_degree: self.variable_degree,
        }
    }
}

/// The budget of [`find_image_annihilator`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ImageAnnihilatorBudget {
    /// Largest order above the parent's order `k` (the order runs over
    /// `k, ..., k + max_order_increase`).  Default 2.
    pub max_order_increase: usize,
    /// Largest derivative order of the ansatz.  Default 4.
    pub max_derivative_order: usize,
    /// Largest coefficient degree in `n`.  Default 6.
    pub max_index_degree: usize,
    /// Largest coefficient degree in `x`.  Default 6.
    pub max_variable_degree: usize,
    /// Largest number of unknowns of one ansatz.  Default 400.
    pub max_unknowns: usize,
    /// Largest number of ansatz shapes solved.  Default 60.
    pub max_shapes: usize,
}

impl Default for ImageAnnihilatorBudget {
    fn default() -> Self {
        Self {
            max_order_increase: 2,
            max_derivative_order: 4,
            max_index_degree: 6,
            max_variable_degree: 6,
            max_unknowns: 400,
            max_shapes: 60,
        }
    }
}

/// Why [`find_image_annihilator`] or a certificate check failed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ImageAnnihilatorError {
    /// The parent annihilator's top shift column is not a single term
    /// `g(n) S^k` free of `x` and of derivatives.
    ParentLeaderNotIndexOnly,
    /// A step is invalid.
    Step(RowMapError),
    /// The recorded annihilator is not runnable.
    NotRunnable,
    /// The recorded bounds are below the annihilator's shape.
    BoundsTooSmall,
    /// The recorded multiplier is zero.
    InvalidMultiplier,
    /// The recorded row annihilator is zero or has a shift, or a row
    /// quotient is recorded without one.
    InvalidRowAnnihilator,
    /// The identity `w E = H o L + Z o A` does not hold.
    IdentityFails,
}

impl fmt::Display for ImageAnnihilatorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ParentLeaderNotIndexOnly => f.write_str(
                "the parent annihilator's leader must be a single term g(n) S^k free of x and \
                 of derivatives",
            ),
            Self::Step(error) => write!(f, "invalid step: {error}"),
            Self::NotRunnable => f.write_str(
                "the annihilator is not runnable (its top shift column is not a single term \
                 without derivatives)",
            ),
            Self::BoundsTooSmall => f.write_str("the pullback bounds are below the shape"),
            Self::InvalidMultiplier => f.write_str("the multiplier must be nonzero"),
            Self::InvalidRowAnnihilator => f.write_str(
                "the row annihilator must be nonzero and shift-free, and a row quotient needs one",
            ),
            Self::IdentityFails => f.write_str("the identity w E = H o L + Z o A does not hold"),
        }
    }
}

impl std::error::Error for ImageAnnihilatorError {}

impl From<RowMapError> for ImageAnnihilatorError {
    fn from(error: RowMapError) -> Self {
        Self::Step(error)
    }
}

/// Whether `op` is runnable: nonzero, with a top shift column that is a
/// single term without derivatives.
pub fn is_runnable(op: &OreOperator) -> bool {
    if op.is_zero() {
        return false;
    }
    let column = op.column(op.max_shift());
    column.term_count() == 1 && column.0.contains_key(&(0, 0))
}

/// The certificate of a runnable image annihilator (module documentation).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ImageAnnihilatorCertificate {
    /// `M`, runnable, acting on the image `G`.
    pub annihilator: OreOperator,
    /// `L`, acting on the source `F`.
    pub parent_annihilator: OreOperator,
    /// The steps `F -> G`, applied first to last.
    pub steps: Vec<RowMapStep>,
    /// The bounds of the pullback of `M`.
    pub bounds: PullbackBounds,
    /// A row annihilator `A` of the source used in the identity, if any
    /// (its invariance is a separate [`crate::AnnihilatorSystemCertificate`]).
    pub row_annihilator: Option<OreOperator>,
    /// `w(n, x)`, nonzero (a polynomial in `n` only without a row
    /// annihilator).
    pub multiplier: BivariatePolynomial,
    /// `H`.
    pub quotient: OreOperator,
    /// `Z` (zero without a row annihilator).
    pub row_quotient: OreOperator,
}

impl ImageAnnihilatorCertificate {
    /// The pullback `(E, b)` of the annihilator (recomputed).
    pub fn pulled_back(&self) -> Result<(OreOperator, Vec<Rational>), ImageAnnihilatorError> {
        let actual = PullbackBounds::of(&self.annihilator);
        if actual.top_shift > self.bounds.top_shift
            || actual.derivative_order > self.bounds.derivative_order
            || actual.variable_degree > self.bounds.variable_degree
        {
            return Err(ImageAnnihilatorError::BoundsTooSmall);
        }
        Ok(pull_back(&self.annihilator, &self.steps, self.bounds)?)
    }

    /// Recompute the pullback and check `w E = H o L + Z o A` exactly, with
    /// `M` runnable, `w` nonzero, `A` (if any) nonzero and shift-free, and
    /// the parent leader free of `x` and of derivatives.  Returns the base
    /// factor `b(n)`.
    pub fn verify(&self) -> Result<Vec<Rational>, ImageAnnihilatorError> {
        if !is_runnable(&self.annihilator) {
            return Err(ImageAnnihilatorError::NotRunnable);
        }
        index_only_leader(&self.parent_annihilator)
            .ok_or(ImageAnnihilatorError::ParentLeaderNotIndexOnly)?;
        if self.multiplier.is_zero() {
            return Err(ImageAnnihilatorError::InvalidMultiplier);
        }
        let row_part = match &self.row_annihilator {
            Some(a) if !a.is_zero() && a.max_shift() == 0 => self.row_quotient.compose(a),
            None if self.row_quotient.is_zero() => OreOperator::default(),
            _ => return Err(ImageAnnihilatorError::InvalidRowAnnihilator),
        };
        let (pulled, base_factor) = self.pulled_back()?;
        let left = pulled.left_multiply(&self.multiplier);
        let right = self
            .quotient
            .compose(&self.parent_annihilator)
            .add(&row_part);
        if left != right {
            return Err(ImageAnnihilatorError::IdentityFails);
        }
        Ok(base_factor)
    }
}

/// How an annihilator was found.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ImageAnnihilatorMethod {
    /// Pushed forward through the steps.
    PushForward,
    /// A common left multiple from the ansatz of this shape.
    LeftMultiple(AnsatzShape),
}

/// The outcome of [`find_image_annihilator`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ImageAnnihilatorSearch {
    /// A runnable annihilator with a verified certificate.
    Found {
        /// The certificate.
        certificate: Box<ImageAnnihilatorCertificate>,
        /// How it was found.
        method: ImageAnnihilatorMethod,
    },
    /// No runnable annihilator within the budget (no statement about larger
    /// shapes).
    NotFoundWithinBudget {
        /// The number of ansatz shapes solved.
        shapes_tried: usize,
        /// The largest shape solved, by number of unknowns.
        largest_shape: Option<AnsatzShape>,
    },
}

/// Divide `op` by the largest left factors `x^e` and `c(n)` of all its
/// coefficients and scale it to coprime integer coefficients with a positive
/// leading coefficient.
pub fn normalize_annihilator(op: &OreOperator) -> OreOperator {
    if op.is_zero() {
        return op.clone();
    }
    let stripped = op.divide_by_variable_power(op.variable_valuation());
    // Content in n: the gcd of every x-column of every coefficient.
    let mut content = vec![Rational::zero()];
    for c in stripped.0.values() {
        for j in 0..=c.variable_degree().unwrap_or(0) {
            content = poly::gcd(&content, &c.variable_coefficient(j));
        }
    }
    let mut out = OreOperator::default();
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
    // Coprime integer coefficients, positive leading coefficient.
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
/// `bounds`, pseudo-divide by `L`, reduce a nonzero remainder by the row
/// annihilator `A` when one is given, and require a zero remainder.
pub fn certify_image_annihilator(
    candidate: &OreOperator,
    parent: &OreOperator,
    row_annihilator: Option<&OreOperator>,
    steps: &[RowMapStep],
    bounds: PullbackBounds,
) -> Result<Option<ImageAnnihilatorCertificate>, ImageAnnihilatorError> {
    if !is_runnable(candidate) {
        return Ok(None);
    }
    index_only_leader(parent).ok_or(ImageAnnihilatorError::ParentLeaderNotIndexOnly)?;
    let actual = PullbackBounds::of(candidate);
    if actual.top_shift > bounds.top_shift
        || actual.derivative_order > bounds.derivative_order
        || actual.variable_degree > bounds.variable_degree
    {
        return Err(ImageAnnihilatorError::BoundsTooSmall);
    }
    let (pulled, _) = pull_back(candidate, steps, bounds)?;
    let Some(division) = right_pseudo_division(&pulled, parent, bounds.top_shift) else {
        return Ok(None);
    };
    let mut certificate = ImageAnnihilatorCertificate {
        annihilator: candidate.clone(),
        parent_annihilator: parent.clone(),
        steps: steps.to_vec(),
        bounds,
        row_annihilator: None,
        multiplier: division.multiplier,
        quotient: division.quotient,
        row_quotient: OreOperator::default(),
    };
    if !division.remainder.is_zero() {
        let Some(a) = row_annihilator else {
            return Ok(None);
        };
        let (order, _) =
            index_only_leader(parent).ok_or(ImageAnnihilatorError::ParentLeaderNotIndexOnly)?;
        let top = division.remainder.max_derivative();
        let Some(reduction) = reduce_by_row_annihilator(&division.remainder, a, order, top) else {
            return Ok(None);
        };
        if !reduction.remainder.is_zero() {
            return Ok(None);
        }
        certificate.row_annihilator = Some(a.clone());
        certificate.multiplier = certificate.multiplier.multiply(&reduction.multiplier);
        certificate.quotient = certificate.quotient.left_multiply(&reduction.multiplier);
        certificate.row_quotient = reduction.quotient;
    }
    // Defensive: the returned certificate always replays.
    Ok(certificate.verify().is_ok().then_some(certificate))
}

/// Normalise and certify, falling back to the unnormalised operator.
fn certify_normalized(
    candidate: &OreOperator,
    parent: &OreOperator,
    row_annihilator: Option<&OreOperator>,
    steps: &[RowMapStep],
    fallback_bounds: PullbackBounds,
) -> Result<Option<ImageAnnihilatorCertificate>, ImageAnnihilatorError> {
    let normalized = normalize_annihilator(candidate);
    for (operator, bounds) in [
        (&normalized, PullbackBounds::of(&normalized)),
        (&normalized, fallback_bounds),
        (candidate, PullbackBounds::of(candidate)),
        (candidate, fallback_bounds),
    ] {
        let actual = PullbackBounds::of(operator);
        if actual.top_shift > bounds.top_shift
            || actual.derivative_order > bounds.derivative_order
            || actual.variable_degree > bounds.variable_degree
        {
            continue;
        }
        if let Some(certificate) =
            certify_image_annihilator(operator, parent, row_annihilator, steps, bounds)?
        {
            return Ok(Some(certificate));
        }
    }
    Ok(None)
}

/// The shapes within the budget: the pairs (order, derivative order) by
/// increasing sum (then order), and for each the coefficient degrees
/// `0, 1, 2, 4, ...` up to the maximum (equal in `n` and `x`, capped by the
/// separate maxima).
pub(crate) fn shapes(order: usize, budget: &ImageAnnihilatorBudget) -> Vec<AnsatzShape> {
    let max_degree = budget.max_index_degree.max(budget.max_variable_degree);
    let mut degrees = vec![0, 1];
    while *degrees.last().expect("nonempty") < max_degree {
        let next = (degrees.last().expect("nonempty") * 2).min(max_degree);
        degrees.push(next);
    }
    let mut pairs = Vec::new();
    for r in order.max(1)..=order + budget.max_order_increase {
        for d in 0..=budget.max_derivative_order {
            pairs.push((r, d));
        }
    }
    pairs.sort_by_key(|(r, d)| (r - order + d, *r));
    let mut out = Vec::new();
    for (r, d) in pairs {
        for g in &degrees {
            let shape = AnsatzShape {
                order: r,
                derivative_order: d,
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

/// One unknown: the monomial `n^a x^b D^d S^i`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) struct Unknown {
    pub(crate) shift: usize,
    pub(crate) derivative: usize,
    pub(crate) index_power: usize,
    pub(crate) variable_power: usize,
}

impl Unknown {
    /// The monomial without its power of `n`.
    pub(crate) fn index_free_operator(&self) -> OreOperator {
        OreOperator::single(
            self.shift,
            self.derivative,
            BivariatePolynomial::monomial(0, self.variable_power, Rational::one()),
        )
    }
}

/// The unknowns of an ansatz, the top-shift ones last (group 1, the others
/// group 0).
pub(crate) fn ansatz_unknowns(shape: &AnsatzShape) -> (Vec<Unknown>, Vec<usize>) {
    let mut unknowns = Vec::new();
    let mut groups = Vec::new();
    for i in 0..shape.order {
        for d in 0..=shape.derivative_order {
            for a in 0..=shape.index_degree {
                for b in 0..=shape.variable_degree {
                    unknowns.push(Unknown {
                        shift: i,
                        derivative: d,
                        index_power: a,
                        variable_power: b,
                    });
                    groups.push(0);
                }
            }
        }
    }
    for a in 0..=shape.index_degree {
        for b in 0..=shape.variable_degree {
            unknowns.push(Unknown {
                shift: shape.order,
                derivative: 0,
                index_power: a,
                variable_power: b,
            });
            groups.push(1);
        }
    }
    (unknowns, groups)
}

/// The monomial of a remainder row an equation collects: entry, shift,
/// derivative order, power of `n`, power of `x`.
type EquationKey = (usize, usize, usize, usize, usize);

/// A kernel vector, with a nonzero top-shift entry, of the linear map
/// sending the coefficients of `M = sum_u c_u u` to the remainder row
/// `sum_u c_u n^(a_u) R(u)`, where `remainder_of(u)` is the remainder row
/// of the monomial `u` without its power of `n` (the reductions commute
/// with left multiplication by `n^a`); `None` when there is none.
pub(crate) fn kernel_candidate(
    unknowns: &[Unknown],
    groups: &[usize],
    remainder_of: impl Fn(&Unknown) -> Vec<OreOperator>,
) -> Option<OreOperator> {
    // Equations: (entry, shift, derivative, n power, x power) -> column -> value.
    let mut equations: BTreeMap<EquationKey, BTreeMap<usize, Rational>> = BTreeMap::new();
    for (column, unknown) in unknowns.iter().enumerate() {
        let factor = BivariatePolynomial::monomial(unknown.index_power, 0, Rational::one());
        for (entry, remainder) in remainder_of(unknown).iter().enumerate() {
            for ((i, d), c) in &remainder.left_multiply(&factor).0 {
                for (p, row) in c.rows().iter().enumerate() {
                    for (q, value) in row.iter().enumerate() {
                        if !value.is_zero() {
                            equations
                                .entry((entry, *i, *d, p, q))
                                .or_default()
                                .insert(column, value.clone());
                        }
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
    let vectors = exact_kernel_vectors_by_group(&rows, groups, 2);
    let vector = vectors[1].as_ref()?;
    let mut candidate = OreOperator::default();
    for (unknown, value) in unknowns.iter().zip(vector) {
        if !value.is_zero() {
            candidate.add_term(
                (unknown.shift, unknown.derivative),
                BivariatePolynomial::monomial(
                    unknown.index_power,
                    unknown.variable_power,
                    value.clone(),
                ),
            );
        }
    }
    Some(candidate)
}

/// Solve one ansatz: a runnable `M` with zero remainder, or `None`.
fn solve_shape(
    shape: &AnsatzShape,
    parent: &OreOperator,
    row_annihilator: Option<&OreOperator>,
    steps: &[RowMapStep],
) -> Result<Option<OreOperator>, ImageAnnihilatorError> {
    let (unknowns, groups) = ansatz_unknowns(shape);
    let bounds = shape.bounds();
    let (order, _) =
        index_only_leader(parent).ok_or(ImageAnnihilatorError::ParentLeaderNotIndexOnly)?;
    // Remainders of the index-free monomials x^b D^d S^i; that of n^a times
    // a monomial is n^a times it (left multiplication by n^a commutes with
    // every pullback and reduction step).
    let mut remainders: HashMap<(usize, usize, usize), OreOperator> = HashMap::new();
    for unknown in &unknowns {
        let key = (unknown.shift, unknown.derivative, unknown.variable_power);
        if remainders.contains_key(&key) {
            continue;
        }
        let (pulled, _) = pull_back(&unknown.index_free_operator(), steps, bounds)?;
        let division = right_pseudo_division(&pulled, parent, shape.order)
            .ok_or(ImageAnnihilatorError::ParentLeaderNotIndexOnly)?;
        remainders.insert(key, division.remainder);
    }
    if let Some(a) = row_annihilator {
        // One fixed derivative top for every monomial keeps it linear.
        let top = remainders
            .values()
            .map(OreOperator::max_derivative)
            .max()
            .unwrap_or(0);
        for remainder in remainders.values_mut() {
            let reduction = reduce_by_row_annihilator(remainder, a, order, top)
                .ok_or(ImageAnnihilatorError::InvalidRowAnnihilator)?;
            *remainder = reduction.remainder;
        }
    }
    Ok(kernel_candidate(&unknowns, &groups, |unknown| {
        vec![remainders[&(unknown.shift, unknown.derivative, unknown.variable_power)].clone()]
    }))
}

/// A runnable annihilator of the image of a sequence annihilated by
/// `parent` under `steps`, with a verified certificate, or the report that
/// none was found within `budget` (module documentation).
pub fn find_image_annihilator(
    parent: &OreOperator,
    row_annihilator: Option<&OreOperator>,
    steps: &[RowMapStep],
    budget: &ImageAnnihilatorBudget,
) -> Result<ImageAnnihilatorSearch, ImageAnnihilatorError> {
    let (order, _) =
        index_only_leader(parent).ok_or(ImageAnnihilatorError::ParentLeaderNotIndexOnly)?;
    for step in steps {
        step.validate()?;
    }
    if steps
        .iter()
        .all(|step| !matches!(step, RowMapStep::WeylOperator(_)))
    {
        let pushed = steps
            .iter()
            .try_fold(parent.clone(), |op, step| push_forward(&op, step));
        if let Some(pushed) = pushed.filter(is_runnable) {
            let bounds = PullbackBounds::of(&pushed);
            if let Some(certificate) =
                certify_normalized(&pushed, parent, row_annihilator, steps, bounds)?
            {
                return Ok(ImageAnnihilatorSearch::Found {
                    certificate: Box::new(certificate),
                    method: ImageAnnihilatorMethod::PushForward,
                });
            }
        }
    }
    let mut failed: Vec<AnsatzShape> = Vec::new();
    let mut tried = 0;
    let mut largest = None;
    for shape in shapes(order, budget) {
        if tried >= budget.max_shapes {
            break;
        }
        if failed.iter().any(|f| shape.dominated_by(f)) {
            continue;
        }
        tried += 1;
        largest = Some(shape);
        match solve_shape(&shape, parent, row_annihilator, steps)? {
            Some(candidate) => {
                if let Some(certificate) =
                    certify_normalized(&candidate, parent, row_annihilator, steps, shape.bounds())?
                {
                    return Ok(ImageAnnihilatorSearch::Found {
                        certificate: Box::new(certificate),
                        method: ImageAnnihilatorMethod::LeftMultiple(shape),
                    });
                }
                failed.push(shape);
            }
            None => failed.push(shape),
        }
    }
    Ok(ImageAnnihilatorSearch::NotFoundWithinBudget {
        shapes_tried: tried,
        largest_shape: largest,
    })
}
