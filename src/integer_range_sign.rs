//! Exact sign decisions for rational functions of an integer index.
//!
//! [`decide_index_rational_function_sign`] decides statements of the form
//! "`f(n) R 0` for every integer `n` in `[first, last]`" (or `[first, +inf)`),
//! where `f = p / q` with `p`, `q` polynomials in `n` with exact rational
//! coefficients (ascending order: entry `i` is the coefficient of `n^i`) and
//! `R` one of `> 0`, `>= 0`, `< 0`, `<= 0`, `!= 0`, `= 0`
//! ([`IndexSignCondition`]).  [`decide_index_polynomial_sign`] is the case
//! `q = 1`.  No floating point is used.
//!
//! # Conventions
//!
//! - A vanishing denominator at an index of the range is a failure for every
//!   condition ([`IndexSignDecision::FailsAt`] with `value: None`); nothing is
//!   cancelled, so `p / q` is not reduced first.
//! - The zero numerator gives `f = 0` wherever `q != 0`: the conditions
//!   `>= 0`, `<= 0` and `= 0` hold, the others fail at the first index.
//! - An empty range (`last < first`) satisfies every condition.
//! - A failure reports the smallest failing index of the range.
//!
//! # The algorithm
//!
//! 1. **Taylor shift (fast path).**  When every coefficient of `p(first + u)`
//!    has one sign, and likewise for `q` with a nonzero constant term, the sign
//!    of `f(first + u)` is known for `u = 0` and for every real `u > 0`; if
//!    the condition holds for both, this is the certificate
//!    ([`IndexSignMethod::TaylorShift`]).
//! 2. **Sign samples.**  Otherwise each nonzero polynomial is analysed by
//!    exact Sturm root counting ([`crate::sturm::SturmChain`]) on integer
//!    windows.  With `U` the ceiling of the Cauchy bound (every real root lies
//!    in `(-U, U)`), the part `[max(first, -U), min(last, U)]` of the range is
//!    bisected at integer midpoints: a window `(l, h)` without roots is not
//!    split, and a window of width one with roots is a root cell
//!    ([`IndexRootCell`]).  The sample indices are `first`, `last` (when
//!    finite), the window ends, and, in a root-free window both of whose ends
//!    are roots, the index `l + 1`.  Between two consecutive samples `s < t`
//!    with `t - s >= 2` the polynomial then has no root in `(s, t)` and does
//!    not vanish at `s` or at `t`, so its sign at every integer strictly
//!    between them is known; beyond the last sample of an unbounded range it
//!    has no root and keeps the (nonzero) sign of that sample.  The samples of
//!    `p` and `q` are merged (a refinement keeps both properties), `f` is
//!    evaluated exactly at every sample, and the sign at every other index is
//!    read off ([`IndexSignMethod::SignSamples`]).
//!
//! Both routes are complete: every input within the budget is decided.  The
//! number of Sturm counts is `O(deg * log(U))` per polynomial (each root
//! needs one bisection path), independent of the length of the range, so a
//! root like `n = 10^30` is located without enumerating the integers before
//! it.  [`IndexSignDecision::Undecided`] is returned only when the budget
//! ([`IndexSignBudget`]: a degree limit and a limit on the bisection steps)
//! is exceeded.
//!
//! # Certificates
//!
//! A [`IndexSignCertificate`] records the input, the condition and the method
//! data (the shifted coefficients, or the samples with the signs of `p` and
//! `q` and the root cells).  [`IndexSignCertificate::verify`] replays it from
//! the record alone: it recomputes the signs, recounts the roots between
//! consecutive samples and beyond the last one, and rechecks the condition at
//! every index class.
//!
//! # JSON
//!
//! Every type serializes with serde.  Integers (indices) and rationals are
//! exact strings (`"12"`, `"-3/4"`), as elsewhere in the crate; signs,
//! conditions and methods are snake-case tags.
//!
//! Related: [`crate::sturm`] (root counting), [`crate::recurrence`] (the lift
//! reports that carry these decisions as hypothesis evidence).

use crate::recurrence::{format_rational_coeff, parse_rational_coeff};
use crate::sturm::SturmChain;
use crate::vec_poly;
use crate::BigRational;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
use std::fmt;

type Q = BigRational;

/// Default largest degree of the numerator and the denominator accepted by
/// [`decide_index_rational_function_sign`].
pub const DEFAULT_INDEX_SIGN_MAXIMUM_DEGREE: usize = 64;

/// Default largest number of bisection steps (Sturm root counts on a window)
/// per decision.
pub const DEFAULT_INDEX_SIGN_MAXIMUM_BISECTION_STEPS: usize = 100_000;

/// The integers `n` with `first <= n <= last`; `last = None` is unbounded.
/// The range is empty when `last < first`.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IntegerRange {
    /// The first index (exact integer string in JSON).
    #[serde(with = "exact_serde::integer")]
    pub first: BigInt,
    /// The last index, or `None` (JSON `null`) for `+infinity`.
    #[serde(default, with = "exact_serde::optional_integer")]
    pub last: Option<BigInt>,
}

impl IntegerRange {
    /// The unbounded range `[first, +infinity)`.
    pub fn starting_at(first: impl Into<BigInt>) -> Self {
        Self {
            first: first.into(),
            last: None,
        }
    }

    /// The bounded range `[first, last]` (empty when `last < first`).
    pub fn closed(first: impl Into<BigInt>, last: impl Into<BigInt>) -> Self {
        Self {
            first: first.into(),
            last: Some(last.into()),
        }
    }

    /// Whether the range has no integer.
    pub fn is_empty(&self) -> bool {
        self.last.as_ref().is_some_and(|last| last < &self.first)
    }

    /// Whether `index` belongs to the range.
    pub fn contains(&self, index: &BigInt) -> bool {
        index >= &self.first && self.last.as_ref().is_none_or(|last| index <= last)
    }

    /// Whether the range has exactly one integer.
    fn is_single_index(&self) -> bool {
        self.last.as_ref() == Some(&self.first)
    }
}

impl fmt::Display for IntegerRange {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.last {
            Some(last) => write!(f, "[{}, {last}]", self.first),
            None => write!(f, "[{}, +inf)", self.first),
        }
    }
}

/// The sign of an exact value.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum IndexSign {
    /// `< 0`.
    Negative,
    /// `= 0`.
    Zero,
    /// `> 0`.
    Positive,
}

impl IndexSign {
    /// The sign of `value`.
    pub fn of(value: &BigRational) -> Self {
        if value.is_positive() {
            Self::Positive
        } else if value.is_negative() {
            Self::Negative
        } else {
            Self::Zero
        }
    }

    /// The sign of a product.
    pub fn times(self, other: Self) -> Self {
        match (self, other) {
            (Self::Zero, _) | (_, Self::Zero) => Self::Zero,
            (a, b) if a == b => Self::Positive,
            _ => Self::Negative,
        }
    }
}

/// The relation `f(n) R 0` to decide.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum IndexSignCondition {
    /// `f(n) > 0`.
    Positive,
    /// `f(n) >= 0`.
    Nonnegative,
    /// `f(n) < 0`.
    Negative,
    /// `f(n) <= 0`.
    Nonpositive,
    /// `f(n) != 0`.
    Nonzero,
    /// `f(n) = 0`.
    Zero,
}

impl IndexSignCondition {
    /// Whether a value of sign `sign` satisfies the condition.
    pub fn is_satisfied_by(self, sign: IndexSign) -> bool {
        match self {
            Self::Positive => sign == IndexSign::Positive,
            Self::Nonnegative => sign != IndexSign::Negative,
            Self::Negative => sign == IndexSign::Negative,
            Self::Nonpositive => sign != IndexSign::Positive,
            Self::Nonzero => sign != IndexSign::Zero,
            Self::Zero => sign == IndexSign::Zero,
        }
    }

    /// The relation as text: `"> 0"`, `">= 0"`, `"< 0"`, `"<= 0"`, `"!= 0"`
    /// or `"= 0"`.
    pub fn relation_text(self) -> &'static str {
        match self {
            Self::Positive => "> 0",
            Self::Nonnegative => ">= 0",
            Self::Negative => "< 0",
            Self::Nonpositive => "<= 0",
            Self::Nonzero => "!= 0",
            Self::Zero => "= 0",
        }
    }
}

/// Resource limits of [`decide_index_rational_function_sign_with_budget`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IndexSignBudget {
    /// Largest degree of the numerator and of the denominator (after
    /// trimming).  Default [`DEFAULT_INDEX_SIGN_MAXIMUM_DEGREE`].
    pub maximum_degree: usize,
    /// Largest number of bisection steps, each one exact Sturm root count on
    /// an integer window, over both polynomials.  Default
    /// [`DEFAULT_INDEX_SIGN_MAXIMUM_BISECTION_STEPS`].
    pub maximum_bisection_steps: usize,
}

impl Default for IndexSignBudget {
    fn default() -> Self {
        Self {
            maximum_degree: DEFAULT_INDEX_SIGN_MAXIMUM_DEGREE,
            maximum_bisection_steps: DEFAULT_INDEX_SIGN_MAXIMUM_BISECTION_STEPS,
        }
    }
}

/// The value signs at one sample index.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IndexSignSample {
    /// The index `n` (exact integer string in JSON).
    #[serde(with = "exact_serde::integer")]
    pub index: BigInt,
    /// The sign of `p(n)`.
    pub numerator_sign: IndexSign,
    /// The sign of `q(n)`.
    pub denominator_sign: IndexSign,
}

/// An open unit window `(lower, lower + 1)` of the range that contains
/// `root_count` distinct real roots of a polynomial.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IndexRootCell {
    /// The integer `lower` (exact integer string in JSON).
    #[serde(with = "exact_serde::integer")]
    pub lower: BigInt,
    /// The number of distinct real roots in `(lower, lower + 1)` (Sturm).
    pub root_count: usize,
}

/// How a [`IndexSignCertificate`] establishes the condition.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum IndexSignMethod {
    /// The range is empty, so the statement holds vacuously.
    EmptyRange,
    /// Every coefficient of `p(first + u)` has one sign, and likewise for
    /// `q(first + u)`, whose constant term is nonzero; this fixes the sign of
    /// `f(first + u)` at `u = 0` and for every real `u > 0`.
    TaylorShift {
        /// The coefficients of `p(first + u)`, ascending in `u`.
        #[serde(with = "exact_serde::rational_vector")]
        shifted_numerator: Vec<BigRational>,
        /// The coefficients of `q(first + u)`, ascending in `u`.
        #[serde(with = "exact_serde::rational_vector")]
        shifted_denominator: Vec<BigRational>,
    },
    /// Exact signs at sample indices that determine the sign at every index
    /// of the range (module documentation).
    SignSamples {
        /// The samples, by increasing index; the first is `first`, the last
        /// is `last` for a bounded range.
        samples: Vec<IndexSignSample>,
        /// The windows of the range containing roots of `p`, by increasing
        /// `lower`; integer roots appear as samples of sign zero instead.
        numerator_root_cells: Vec<IndexRootCell>,
        /// The same for `q`.
        denominator_root_cells: Vec<IndexRootCell>,
    },
}

/// A replayable proof that `p(n) / q(n) R 0` for every integer `n` of the
/// range.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IndexSignCertificate {
    /// `p`, ascending coefficients in `n`, trimmed (exact strings in JSON).
    #[serde(with = "exact_serde::rational_vector")]
    pub numerator: Vec<BigRational>,
    /// `q`, ascending coefficients in `n`, trimmed.
    #[serde(with = "exact_serde::rational_vector")]
    pub denominator: Vec<BigRational>,
    /// The range.
    pub range: IntegerRange,
    /// The relation.
    pub condition: IndexSignCondition,
    /// The method and its data.
    pub method: IndexSignMethod,
}

/// The decision of [`decide_index_rational_function_sign`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "decision", rename_all = "snake_case", deny_unknown_fields)]
pub enum IndexSignDecision {
    /// The condition holds at every index of the range.
    Holds {
        /// The replayable certificate.
        certificate: Box<IndexSignCertificate>,
    },
    /// The condition fails at `index`, the smallest failing index of the
    /// range.
    FailsAt {
        /// The failing index (exact integer string in JSON).
        #[serde(with = "exact_serde::integer")]
        index: BigInt,
        /// `f(index)`, or `None` (JSON `null`) when the denominator vanishes
        /// there.
        #[serde(default, with = "exact_serde::optional_rational")]
        value: Option<BigRational>,
    },
    /// The budget was exceeded; nothing is claimed.
    Undecided {
        /// Which limit was hit.
        reason: String,
    },
}

impl IndexSignDecision {
    /// Whether the condition was proved.
    pub fn holds(&self) -> bool {
        matches!(self, Self::Holds { .. })
    }

    /// The failing index, for [`IndexSignDecision::FailsAt`].
    pub fn failing_index(&self) -> Option<&BigInt> {
        match self {
            Self::FailsAt { index, .. } => Some(index),
            _ => None,
        }
    }

    /// A one-sentence account of the decision, naming the function `name`
    /// (for example `"a_n"`).
    pub fn describe(&self, name: &str) -> String {
        match self {
            Self::Holds { certificate } => {
                let condition = certificate.condition.relation_text();
                let range = &certificate.range;
                let method = match &certificate.method {
                    IndexSignMethod::EmptyRange => "the range is empty".to_string(),
                    IndexSignMethod::TaylorShift { .. } => format!(
                        "Taylor shift: the coefficients of {name} at n = {} + u have one sign",
                        range.first
                    ),
                    IndexSignMethod::SignSamples {
                        samples,
                        numerator_root_cells,
                        denominator_root_cells,
                    } => format!(
                        "exact Sturm root isolation: {} sample indices, {} root windows",
                        samples.len(),
                        numerator_root_cells.len() + denominator_root_cells.len()
                    ),
                };
                format!("{name} {condition} for every integer n in {range} (exact; {method})")
            }
            Self::FailsAt { index, value } => match value {
                Some(value) => format!(
                    "fails at n = {index}: {name} = {} there (exact)",
                    format_rational_coeff(value)
                ),
                None => format!("fails at n = {index}: the denominator of {name} vanishes there"),
            },
            Self::Undecided { reason } => format!("undecided (budget): {reason}"),
        }
    }
}

/// Why [`IndexSignCertificate::verify`] rejected a certificate.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum IndexSignCertificateError {
    /// The recorded method does not fit the input (an empty-range record for
    /// a nonempty range, a zero denominator, recorded shifted coefficients
    /// that differ from the recomputed ones, or coefficients of both signs).
    MethodDoesNotApply(String),
    /// The sample indices are not increasing, leave the range, or do not
    /// start at `first` (and end at `last` for a bounded range).
    SampleIndices(String),
    /// A recorded sign differs from the recomputed one.
    SampleSign {
        /// The sample index.
        index: BigInt,
    },
    /// A polynomial has a root strictly between two consecutive samples, or
    /// vanishes at both of them.
    RootBetweenSamples {
        /// The lower sample.
        lower: BigInt,
        /// The upper sample.
        upper: BigInt,
    },
    /// For an unbounded range, a polynomial vanishes at the last sample or
    /// has a root beyond it.
    RootAfterLastSample {
        /// The last sample.
        index: BigInt,
    },
    /// A root cell has a wrong or zero root count, or leaves the range.
    RootCell {
        /// The cell's `lower`.
        lower: BigInt,
    },
    /// The condition fails at this index.
    ConditionFails {
        /// The index.
        index: BigInt,
    },
}

impl fmt::Display for IndexSignCertificateError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MethodDoesNotApply(reason) => {
                write!(f, "the certificate method does not apply: {reason}")
            }
            Self::SampleIndices(reason) => write!(f, "invalid sample indices: {reason}"),
            Self::SampleSign { index } => write!(f, "wrong recorded sign at n = {index}"),
            Self::RootBetweenSamples { lower, upper } => write!(
                f,
                "a polynomial has a root in ({lower}, {upper}) or vanishes at both ends"
            ),
            Self::RootAfterLastSample { index } => write!(
                f,
                "a polynomial vanishes at the last sample {index} or has a root beyond it"
            ),
            Self::RootCell { lower } => {
                write!(f, "the root window ({lower}, {lower} + 1) is wrong")
            }
            Self::ConditionFails { index } => write!(f, "the condition fails at n = {index}"),
        }
    }
}

impl std::error::Error for IndexSignCertificateError {}

// ---------------------------------------------------------------------------
// Polynomial helpers
// ---------------------------------------------------------------------------

fn trimmed(p: &[Q]) -> Vec<Q> {
    let mut p = p.to_vec();
    vec_poly::strip_trailing_zeros(&mut p);
    p
}

fn degree(p: &[Q]) -> usize {
    p.len().saturating_sub(1)
}

fn evaluate_at(p: &[Q], index: &BigInt) -> Q {
    if p.is_empty() {
        return Q::zero();
    }
    vec_poly::evaluate(p, &Q::from_integer(index.clone()))
}

fn sign_at(p: &[Q], index: &BigInt) -> IndexSign {
    IndexSign::of(&evaluate_at(p, index))
}

/// `p(first + u)`, trimmed.
fn taylor_shift(p: &[Q], first: &BigInt) -> Vec<Q> {
    if p.is_empty() {
        return Vec::new();
    }
    trimmed(&vec_poly::compose(
        p,
        &[Q::from_integer(first.clone()), Q::one()],
    ))
}

/// The positive integer multiple of `p` with coprime integer coefficients
/// cleared of denominators (same roots, same signs).
fn integer_multiple(p: &[Q]) -> Vec<BigInt> {
    let common = p.iter().fold(BigInt::one(), |acc, c| acc.lcm(c.denom()));
    p.iter()
        .map(|c| (c * Q::from_integer(common.clone())).to_integer())
        .collect()
}

/// An integer `U >= 1` with every real root of the nonconstant polynomial
/// `p` in the open interval `(-U, U)`: the ceiling of the Cauchy bound
/// `1 + max_i |p_i / p_d|`.
fn cauchy_bound(p: &[Q]) -> BigInt {
    let leading = p.last().expect("a nonconstant polynomial").abs();
    let largest = p[..p.len() - 1]
        .iter()
        .map(|c| c.abs() / &leading)
        .max()
        .unwrap_or_else(Q::zero);
    (largest + Q::one()).ceil().to_integer()
}

/// A nonzero polynomial with its Sturm chain (for degree at least one).
struct CountedPolynomial<'a> {
    coefficients: &'a [Q],
    chain: Option<SturmChain>,
}

impl<'a> CountedPolynomial<'a> {
    fn new(coefficients: &'a [Q]) -> Self {
        let chain = (coefficients.len() >= 2)
            .then(|| SturmChain::from_bigint_coeffs(&integer_multiple(coefficients)));
        Self {
            coefficients,
            chain,
        }
    }

    /// Distinct real roots in the open interval `(lower, upper)`, `lower <= upper`.
    fn roots_between(&self, lower: &BigInt, upper: &BigInt) -> usize {
        match &self.chain {
            None => 0,
            Some(chain) => chain.count_roots_in(
                &Q::from_integer(lower.clone()),
                &Q::from_integer(upper.clone()),
            ),
        }
    }

    fn vanishes_at(&self, index: &BigInt) -> bool {
        evaluate_at(self.coefficients, index).is_zero()
    }
}

/// The sample indices and root cells of one nonzero polynomial on a nonempty
/// range (module documentation, step 2).
fn sample_polynomial(
    polynomial: &CountedPolynomial<'_>,
    range: &IntegerRange,
    steps: &mut usize,
    budget: &IndexSignBudget,
) -> Result<(BTreeSet<BigInt>, Vec<IndexRootCell>), String> {
    let mut samples = BTreeSet::new();
    samples.insert(range.first.clone());
    if let Some(last) = &range.last {
        samples.insert(last.clone());
    }
    let mut cells = Vec::new();
    if polynomial.chain.is_none() {
        return Ok((samples, cells));
    }
    let bound = cauchy_bound(polynomial.coefficients);
    let lower = std::cmp::max(range.first.clone(), -bound.clone());
    let upper = match &range.last {
        Some(last) => std::cmp::min(last.clone(), bound),
        None => bound,
    };
    if lower > upper {
        return Ok((samples, cells));
    }
    samples.insert(lower.clone());
    samples.insert(upper.clone());
    let mut windows = vec![(lower, upper)];
    while let Some((low, high)) = windows.pop() {
        if low == high {
            continue;
        }
        if *steps >= budget.maximum_bisection_steps {
            return Err(format!(
                "more than {} bisection steps",
                budget.maximum_bisection_steps
            ));
        }
        *steps += 1;
        let count = polynomial.roots_between(&low, &high);
        let width = &high - &low;
        if count == 0 {
            if width >= BigInt::from(2)
                && polynomial.vanishes_at(&low)
                && polynomial.vanishes_at(&high)
            {
                samples.insert(&low + 1);
            }
        } else if width.is_one() {
            cells.push(IndexRootCell {
                lower: low,
                root_count: count,
            });
        } else {
            let middle = (&low + &high).div_floor(&BigInt::from(2));
            samples.insert(middle.clone());
            windows.push((middle.clone(), high));
            windows.push((low, middle));
        }
    }
    cells.sort_by(|a, b| a.lower.cmp(&b.lower));
    Ok((samples, cells))
}

/// The sign of a polynomial at every index strictly between two consecutive
/// samples (one of the two ends is nonzero), or beyond the last sample.
fn interior_sign(left: IndexSign, right: Option<IndexSign>) -> IndexSign {
    match (left, right) {
        (IndexSign::Zero, Some(right)) => right,
        (left, _) => left,
    }
}

/// The smallest index of the range where the condition fails, given samples
/// that satisfy the gap and tail properties of the module documentation;
/// `(index, value)` with `value = None` for a vanishing denominator.
fn first_failure(
    numerator: &[Q],
    denominator: &[Q],
    samples: &[IndexSignSample],
    range: &IntegerRange,
    condition: IndexSignCondition,
) -> Option<(BigInt, Option<Q>)> {
    let value_at = |index: &BigInt| evaluate_at(numerator, index) / evaluate_at(denominator, index);
    for (position, sample) in samples.iter().enumerate() {
        if sample.denominator_sign == IndexSign::Zero {
            return Some((sample.index.clone(), None));
        }
        if !condition.is_satisfied_by(sample.numerator_sign.times(sample.denominator_sign)) {
            return Some((sample.index.clone(), Some(value_at(&sample.index))));
        }
        let next = samples.get(position + 1);
        let has_interior = match next {
            Some(next) => &next.index - &sample.index >= BigInt::from(2),
            None => range.last.is_none(),
        };
        if has_interior {
            let numerator_sign =
                interior_sign(sample.numerator_sign, next.map(|n| n.numerator_sign));
            let denominator_sign =
                interior_sign(sample.denominator_sign, next.map(|n| n.denominator_sign));
            if !condition.is_satisfied_by(numerator_sign.times(denominator_sign)) {
                let index = &sample.index + 1;
                let value = value_at(&index);
                return Some((index, Some(value)));
            }
        }
    }
    None
}

/// The sign classes `(sign at u = 0, sign for every real u > 0)` of a
/// shifted polynomial with coefficients of one sign; `None` otherwise.
fn taylor_classes(shifted: &[Q]) -> Option<(IndexSign, IndexSign)> {
    let Some(constant) = shifted.first() else {
        return Some((IndexSign::Zero, IndexSign::Zero));
    };
    let at_zero = IndexSign::of(constant);
    if shifted.iter().all(|c| !c.is_negative()) {
        Some((at_zero, IndexSign::Positive))
    } else if shifted.iter().all(|c| !c.is_positive()) {
        Some((at_zero, IndexSign::Negative))
    } else {
        None
    }
}

/// Whether the Taylor-shift data prove the condition on the nonempty range.
fn taylor_shift_proves(
    shifted_numerator: &[Q],
    shifted_denominator: &[Q],
    range: &IntegerRange,
    condition: IndexSignCondition,
) -> bool {
    let (Some((p0, p1)), Some((q0, q1))) = (
        taylor_classes(shifted_numerator),
        taylor_classes(shifted_denominator),
    ) else {
        return false;
    };
    if q0 == IndexSign::Zero || q1 == IndexSign::Zero {
        return false;
    }
    condition.is_satisfied_by(p0.times(q0))
        && (range.is_single_index() || condition.is_satisfied_by(p1.times(q1)))
}

// ---------------------------------------------------------------------------
// Entry points
// ---------------------------------------------------------------------------

/// Decide `p(n) R 0` for every integer `n` of `range`; `p` has ascending
/// exact rational coefficients.  The default [`IndexSignBudget`] applies.
///
/// # Example
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::{
///     decide_index_polynomial_sign, BigRational, IndexSignCondition, IndexSignDecision,
///     IntegerRange,
/// };
///
/// let q = |n: i64| BigRational::from_integer(BigInt::from(n));
/// // (n - 2)(n - 3) = n^2 - 5 n + 6 >= 0 at every integer, < 0 nowhere ...
/// let p = [q(6), q(-5), q(1)];
/// let range = IntegerRange::starting_at(0);
/// assert!(decide_index_polynomial_sign(&p, &range, IndexSignCondition::Nonnegative).holds());
/// // ... but it vanishes at n = 2.
/// let decision = decide_index_polynomial_sign(&p, &range, IndexSignCondition::Positive);
/// assert_eq!(decision.failing_index(), Some(&BigInt::from(2)));
/// if let IndexSignDecision::Holds { certificate } =
///     decide_index_polynomial_sign(&p, &range, IndexSignCondition::Nonnegative)
/// {
///     certificate.verify().unwrap();
/// }
/// ```
pub fn decide_index_polynomial_sign(
    polynomial: &[BigRational],
    range: &IntegerRange,
    condition: IndexSignCondition,
) -> IndexSignDecision {
    decide_index_rational_function_sign(polynomial, &[Q::one()], range, condition)
}

/// Decide `p(n) / q(n) R 0` for every integer `n` of `range` (module
/// documentation); the default [`IndexSignBudget`] applies.
pub fn decide_index_rational_function_sign(
    numerator: &[BigRational],
    denominator: &[BigRational],
    range: &IntegerRange,
    condition: IndexSignCondition,
) -> IndexSignDecision {
    decide_index_rational_function_sign_with_budget(
        numerator,
        denominator,
        range,
        condition,
        &IndexSignBudget::default(),
    )
}

/// [`decide_index_rational_function_sign`] with an explicit budget; the
/// result is [`IndexSignDecision::Undecided`] only when a limit is hit.
pub fn decide_index_rational_function_sign_with_budget(
    numerator: &[BigRational],
    denominator: &[BigRational],
    range: &IntegerRange,
    condition: IndexSignCondition,
    budget: &IndexSignBudget,
) -> IndexSignDecision {
    let numerator = trimmed(numerator);
    let denominator = trimmed(denominator);
    let certificate = |method| IndexSignDecision::Holds {
        certificate: Box::new(IndexSignCertificate {
            numerator: numerator.clone(),
            denominator: denominator.clone(),
            range: range.clone(),
            condition,
            method,
        }),
    };
    if range.is_empty() {
        return certificate(IndexSignMethod::EmptyRange);
    }
    if denominator.is_empty() {
        return IndexSignDecision::FailsAt {
            index: range.first.clone(),
            value: None,
        };
    }
    let largest = degree(&numerator).max(degree(&denominator));
    if largest > budget.maximum_degree {
        return IndexSignDecision::Undecided {
            reason: format!(
                "degree {largest} exceeds the limit {}",
                budget.maximum_degree
            ),
        };
    }

    let shifted_numerator = taylor_shift(&numerator, &range.first);
    let shifted_denominator = taylor_shift(&denominator, &range.first);
    if taylor_shift_proves(&shifted_numerator, &shifted_denominator, range, condition) {
        return certificate(IndexSignMethod::TaylorShift {
            shifted_numerator,
            shifted_denominator,
        });
    }

    let mut steps = 0;
    let mut samples = BTreeSet::new();
    let mut numerator_root_cells = Vec::new();
    let counted_numerator = CountedPolynomial::new(&numerator);
    let counted_denominator = CountedPolynomial::new(&denominator);
    if !numerator.is_empty() {
        match sample_polynomial(&counted_numerator, range, &mut steps, budget) {
            Ok((indices, cells)) => {
                samples.extend(indices);
                numerator_root_cells = cells;
            }
            Err(reason) => return IndexSignDecision::Undecided { reason },
        }
    }
    let denominator_root_cells =
        match sample_polynomial(&counted_denominator, range, &mut steps, budget) {
            Ok((indices, cells)) => {
                samples.extend(indices);
                cells
            }
            Err(reason) => return IndexSignDecision::Undecided { reason },
        };
    let samples: Vec<IndexSignSample> = samples
        .into_iter()
        .map(|index| IndexSignSample {
            numerator_sign: sign_at(&numerator, &index),
            denominator_sign: sign_at(&denominator, &index),
            index,
        })
        .collect();
    if let Some((index, value)) =
        first_failure(&numerator, &denominator, &samples, range, condition)
    {
        return IndexSignDecision::FailsAt { index, value };
    }
    certificate(IndexSignMethod::SignSamples {
        samples,
        numerator_root_cells,
        denominator_root_cells,
    })
}

// ---------------------------------------------------------------------------
// Replay
// ---------------------------------------------------------------------------

impl IndexSignCertificate {
    /// Replay the certificate from its record alone (module documentation):
    /// `Ok(())` proves `numerator(n) / denominator(n) R 0` for every integer
    /// `n` of the range.
    pub fn verify(&self) -> Result<(), IndexSignCertificateError> {
        let numerator = trimmed(&self.numerator);
        let denominator = trimmed(&self.denominator);
        let range = &self.range;
        match &self.method {
            IndexSignMethod::EmptyRange => {
                if range.is_empty() {
                    Ok(())
                } else {
                    Err(IndexSignCertificateError::MethodDoesNotApply(
                        "the range is not empty".into(),
                    ))
                }
            }
            _ if range.is_empty() => Ok(()),
            _ if denominator.is_empty() => Err(IndexSignCertificateError::MethodDoesNotApply(
                "the denominator is the zero polynomial".into(),
            )),
            IndexSignMethod::TaylorShift {
                shifted_numerator,
                shifted_denominator,
            } => {
                if &taylor_shift(&numerator, &range.first) != shifted_numerator
                    || &taylor_shift(&denominator, &range.first) != shifted_denominator
                {
                    return Err(IndexSignCertificateError::MethodDoesNotApply(
                        "the recorded shifted coefficients differ from the recomputed ones".into(),
                    ));
                }
                if taylor_shift_proves(
                    shifted_numerator,
                    shifted_denominator,
                    range,
                    self.condition,
                ) {
                    Ok(())
                } else {
                    Err(IndexSignCertificateError::MethodDoesNotApply(
                        "the shifted coefficients do not fix a sign satisfying the condition"
                            .into(),
                    ))
                }
            }
            IndexSignMethod::SignSamples {
                samples,
                numerator_root_cells,
                denominator_root_cells,
            } => {
                self.verify_samples(&numerator, &denominator, samples)?;
                for (polynomial, cells) in [
                    (&numerator, numerator_root_cells),
                    (&denominator, denominator_root_cells),
                ] {
                    if polynomial.is_empty() {
                        if let Some(cell) = cells.first() {
                            return Err(IndexSignCertificateError::RootCell {
                                lower: cell.lower.clone(),
                            });
                        }
                        continue;
                    }
                    let counted = CountedPolynomial::new(polynomial);
                    verify_root_cells(&counted, cells, range)?;
                    verify_gaps(&counted, samples, range)?;
                }
                match first_failure(&numerator, &denominator, samples, range, self.condition) {
                    Some((index, _)) => Err(IndexSignCertificateError::ConditionFails { index }),
                    None => Ok(()),
                }
            }
        }
    }

    fn verify_samples(
        &self,
        numerator: &[Q],
        denominator: &[Q],
        samples: &[IndexSignSample],
    ) -> Result<(), IndexSignCertificateError> {
        let range = &self.range;
        let indices_error =
            |reason: &str| Err(IndexSignCertificateError::SampleIndices(reason.into()));
        match (samples.first(), samples.last()) {
            (Some(first), Some(last)) => {
                if first.index != range.first {
                    return indices_error("the first sample is not the first index");
                }
                if range.last.as_ref().is_some_and(|end| &last.index != end) {
                    return indices_error("the last sample is not the last index");
                }
            }
            _ => return indices_error("no samples"),
        }
        if samples
            .windows(2)
            .any(|pair| pair[0].index >= pair[1].index)
        {
            return indices_error("the samples are not strictly increasing");
        }
        for sample in samples {
            if sample.numerator_sign != sign_at(numerator, &sample.index)
                || sample.denominator_sign != sign_at(denominator, &sample.index)
            {
                return Err(IndexSignCertificateError::SampleSign {
                    index: sample.index.clone(),
                });
            }
        }
        Ok(())
    }
}

/// Every recorded root cell lies in the range and has its exact root count.
fn verify_root_cells(
    polynomial: &CountedPolynomial<'_>,
    cells: &[IndexRootCell],
    range: &IntegerRange,
) -> Result<(), IndexSignCertificateError> {
    for cell in cells {
        let upper = &cell.lower + 1;
        if cell.root_count == 0
            || !range.contains(&cell.lower)
            || !range.contains(&upper)
            || polynomial.roots_between(&cell.lower, &upper) != cell.root_count
        {
            return Err(IndexSignCertificateError::RootCell {
                lower: cell.lower.clone(),
            });
        }
    }
    Ok(())
}

/// The gap and tail properties of the module documentation for one nonzero
/// polynomial and the merged samples.
fn verify_gaps(
    polynomial: &CountedPolynomial<'_>,
    samples: &[IndexSignSample],
    range: &IntegerRange,
) -> Result<(), IndexSignCertificateError> {
    for pair in samples.windows(2) {
        let (lower, upper) = (&pair[0].index, &pair[1].index);
        if upper - lower < BigInt::from(2) {
            continue;
        }
        if polynomial.roots_between(lower, upper) != 0
            || (polynomial.vanishes_at(lower) && polynomial.vanishes_at(upper))
        {
            return Err(IndexSignCertificateError::RootBetweenSamples {
                lower: lower.clone(),
                upper: upper.clone(),
            });
        }
    }
    if range.last.is_none() {
        let last = &samples.last().expect("samples were checked nonempty").index;
        let beyond = polynomial.chain.is_some() && {
            let bound = cauchy_bound(polynomial.coefficients);
            last < &bound && polynomial.roots_between(last, &bound) != 0
        };
        if beyond || polynomial.vanishes_at(last) {
            return Err(IndexSignCertificateError::RootAfterLastSample {
                index: last.clone(),
            });
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Exact JSON strings
// ---------------------------------------------------------------------------

/// Serde helpers writing integers and rationals as exact strings (shared
/// with the lift reports).
pub(crate) mod exact_serde {
    use super::*;
    use serde::de::Error;
    use serde::{Deserializer, Serializer};

    fn parse_integer<E: Error>(text: &str) -> Result<BigInt, E> {
        text.trim()
            .parse::<BigInt>()
            .map_err(|e| E::custom(format!("invalid integer `{text}`: {e}")))
    }

    fn parse_rational<E: Error>(text: &str) -> Result<BigRational, E> {
        parse_rational_coeff(text).map_err(|e| E::custom(format!("invalid rational `{text}`: {e}")))
    }

    pub(crate) mod integer {
        use super::*;

        pub(crate) fn serialize<S: Serializer>(
            value: &BigInt,
            serializer: S,
        ) -> Result<S::Ok, S::Error> {
            serializer.serialize_str(&value.to_string())
        }

        pub(crate) fn deserialize<'de, D: Deserializer<'de>>(
            deserializer: D,
        ) -> Result<BigInt, D::Error> {
            parse_integer(&String::deserialize(deserializer)?)
        }
    }

    pub(crate) mod optional_integer {
        use super::*;

        pub(crate) fn serialize<S: Serializer>(
            value: &Option<BigInt>,
            serializer: S,
        ) -> Result<S::Ok, S::Error> {
            value
                .as_ref()
                .map(ToString::to_string)
                .serialize(serializer)
        }

        pub(crate) fn deserialize<'de, D: Deserializer<'de>>(
            deserializer: D,
        ) -> Result<Option<BigInt>, D::Error> {
            Option::<String>::deserialize(deserializer)?
                .map(|text| parse_integer(&text))
                .transpose()
        }
    }

    pub(crate) mod optional_rational {
        use super::*;

        pub(crate) fn serialize<S: Serializer>(
            value: &Option<BigRational>,
            serializer: S,
        ) -> Result<S::Ok, S::Error> {
            value
                .as_ref()
                .map(format_rational_coeff)
                .serialize(serializer)
        }

        pub(crate) fn deserialize<'de, D: Deserializer<'de>>(
            deserializer: D,
        ) -> Result<Option<BigRational>, D::Error> {
            Option::<String>::deserialize(deserializer)?
                .map(|text| parse_rational(&text))
                .transpose()
        }
    }

    pub(crate) mod rational_vector {
        use super::*;

        pub(crate) fn serialize<S: Serializer>(
            values: &[BigRational],
            serializer: S,
        ) -> Result<S::Ok, S::Error> {
            values
                .iter()
                .map(format_rational_coeff)
                .collect::<Vec<_>>()
                .serialize(serializer)
        }

        pub(crate) fn deserialize<'de, D: Deserializer<'de>>(
            deserializer: D,
        ) -> Result<Vec<BigRational>, D::Error> {
            Vec::<String>::deserialize(deserializer)?
                .iter()
                .map(|text| parse_rational(text))
                .collect()
        }
    }
}

#[cfg(test)]
mod tests;
