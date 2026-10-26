//! Exact enumeration of real-rooted integer polynomials with prescribed
//! leading coefficients.
//!
//! This module implements the `AllRealRooted` derivative-tower search of
//! Greaves and Syatriadi (*Real-rooted integer polynomial enumeration
//! algorithms and interlacing polynomials via linear programming*,
//! arXiv:2504.09241v1, Lemma 2.5, Theorem 2.6, Proposition 2.8 and
//! Algorithms 1--3), which goes back to McKee and Smyth.  The paper's
//! floating-point endpoint step is replaced by certified root isolation and
//! exact endpoint tests, so every decision is made with exact integer and
//! rational arithmetic.
//!
//! # Entry points
//!
//! - [`enumerate_real_rooted_integer_polynomials_bigint`]: all real-rooted
//!   `p ∈ Z[x]` of degree `n` whose `t >= 3` leading coefficients are given.
//! - [`real_rooted_antiderivative_integer_constants_bigint`] and
//!   [`real_rooted_antiderivative_constant_multiples_bigint`]: the reusable
//!   kernel, the integers `C` (or the multiples `C = k m`) for which
//!   `P(x) + C` is real-rooted, where `P(x) = ∫_0^x q(y) dy`.
//! - [`antiderivative_exact_bigint`] and [`antiderivative_rational_bigint`]:
//!   zero-constant antiderivatives of ascending integer coefficient vectors.
//!
//! # Conventions
//!
//! The prescribed prefix is given in **descending** order, as in the paper:
//! for `p(x) = a_0 x^n + a_1 x^(n-1) + ... + a_n` the input is
//! `[a_0, a_1, ..., a_(t-1)]`.  Every returned polynomial uses Polytool's
//! usual **ascending** order: entry `k` is the coefficient of `x^k`, so the
//! vector is `[a_n, a_(n-1), ..., a_0]` and has length exactly `n + 1`.
//!
//! The leading coefficient may be negative: `p` is real-rooted exactly when
//! `-p` is, and the kernel normalizes the sign of `q` internally (if
//! `P + C` is real-rooted then so is `-P - C`).  No result is negated or
//! made monic; nonmonic prefixes are searched as given.
//!
//! Results are emitted in increasing lexicographic order of the free
//! coefficients `(a_t, a_(t+1), ..., a_n)`, which is the order of the
//! depth-first search.  Reducible polynomials and polynomials with repeated
//! roots are included.  "Real-rooted" follows the rest of Polytool: every
//! complex root is real, counted with multiplicity.
//!
//! # Mathematics
//!
//! For a candidate `p` define the scaled derivative tower
//! `p_i = d^(n-i) p / dx^(n-i)` for `t-1 <= i <= n`.  Then `deg p_i = i`,
//! `p_n = p`, the coefficient of `x^k` in `p_i` is
//! `a_(i-k) (n-i+k)! / k!`, and
//!
//! ```text
//! p_i(x) = ∫_0^x p_(i-1)(y) dy + (n-i)! a_i.
//! ```
//!
//! The first `t` coefficients determine `p_(t-1)`.  By Rolle's theorem every
//! `p_i` of a real-rooted `p` is real-rooted, so the search starts from
//! `p_(t-1)` (rejecting the prefix at once if it is not real-rooted) and, at
//! each level `i = t, ..., n`, adds every admissible constant `(n-i)! a_i`.
//! The antiderivative of `p_(i-1)` always has integer coefficients; the code
//! checks every division exactly instead of assuming it.
//!
//! **Lemma 2.5.** Let `q` be real-rooted and non-constant and
//! `P(x) = ∫_0^x q`.  Then `P + C` is real-rooted if and only if, for every
//! root `µ` of `q`,
//! `P(µ) + C >= 0` when `q'(µ) < 0`, `P(µ) + C <= 0` when `q'(µ) > 0`, and
//! `P(µ) + C = 0` when `q'(µ) = 0` (that is, when `µ` is a multiple root).
//! So the admissible constants form the closed interval
//! `[max -P(µ) over q'(µ) < 0, min -P(µ) over q'(µ) > 0]`, intersected with
//! the single point `-P(µ)` for each multiple root `µ`.  For `deg q >= 2` the
//! set is bounded (a real-rooted `q` of degree at least two either has a
//! multiple root or simple roots of both slopes).  Theorem 2.6 of the paper
//! states the same interval as the two middle elements of the multiset of
//! critical values `{-P(µ)}`, and Proposition 2.7 identifies that multiset
//! with the roots of `Disc_x(P(x) + C)`; using Lemma 2.5 directly avoids
//! sorting algebraic numbers and makes coincident critical values harmless.
//!
//! With `t >= 3` the start `p_(t-1)` has degree at least two, so every level
//! has a bounded interval and the search is finite (Proposition 2.8: the
//! roots satisfy `Σ λ_j^2 = (a_1^2 - 2 a_0 a_2) / a_0^2`, which bounds every
//! coefficient).  With only two prescribed coefficients the first step would
//! ask for the constants making a quadratic real-rooted, a half-line, and
//! the set can be infinite: `x^(n-2) (x^2 - k)` is real-rooted for every
//! integer `k >= 0`.  This is why the API requires `t >= 3`.
//!
//! The paper suggests pruning with the Newton power-sum bound of
//! Proposition 2.8.  Here every node's admissible range is computed exactly,
//! so no node outside it is ever generated; the necessary condition
//! `a_1^2 >= 2 a_0 a_2` is implied by the exact real-rootedness test of the
//! start polynomial `p_(t-1)` (whose own derivative `p_2` must be
//! real-rooted).  The bound is used instead by the brute-force oracle of the
//! test suite.
//!
//! # Certification of the endpoints
//!
//! The kernel never uses floating point.  For the derivative `q` (sign
//! normalized so its leading coefficient is positive):
//!
//! 1. The distinct real roots of `q` are isolated exactly with
//!    [`crate::root_count::isolate_real_roots_bigint_coeffs`] (Sturm
//!    bisection with dyadic endpoints).  If there are fewer distinct real
//!    roots than the degree of the square-free part, `q` is not real-rooted
//!    and no shift of `P` is real-rooted.
//! 2. Each root is classified exactly.  Let `s` be the square-free part and
//!    `g = q / s` (an exact division over `Z` by Gauss's lemma); the multiple
//!    roots of `q` are the roots of `g`.  A root isolated in an open
//!    interval `(l, u)` is multiple exactly when the square-free part of `g`
//!    changes sign on `(l, u)`, and for a simple root the sign of `q'(µ)` is
//!    the sign of `q(u)`.  Roots that the isolation returns exactly (dyadic
//!    rationals) are classified by exact evaluation.
//! 3. Each critical value `-P(µ)` gets a dyadic enclosure
//!    `-P(c) ± B w^2 / 8`, where `c` and `w` are the midpoint and width of
//!    the isolating interval and `B = Σ |q'_k| M^k` bounds `|q'|` on an
//!    interval `[-M, M]` containing every root.  The enclosure is valid
//!    because `q(µ) = 0` gives `|q(x)| <= B |x - µ|` and hence
//!    `|P(µ) - P(c)| <= B (µ - c)^2 / 2`.  The interval is bisected (exact
//!    integer sign evaluation of `s` at dyadic points) until the enclosure is
//!    *settled*: shorter than the required modulus `m` or free of multiples
//!    of `m`.  Only critical points that can still attain the extreme value
//!    are refined.
//! 4. The enclosure of `L` is `[max lower ends, max upper ends]` over the
//!    falling roots, and it lies inside the (settled) enclosure attaining the
//!    maximum upper end; likewise for `U` with minima over the rising roots.
//!    A multiple of `m` between the upper end of `L`'s enclosure and the
//!    lower end of `U`'s enclosure is certainly admissible, one outside
//!    `[lower end of L, upper end of U]` is certainly not, and at most one
//!    multiple remains undecided at each end.  It is decided by an exact
//!    real-rootedness test of the integer polynomial `D P + D C` (with `D`
//!    the common denominator of `P`; `D = 1` along the tower), which also
//!    certifies endpoints where `P + C` acquires a repeated root.  With a
//!    multiple root of `q`, the only candidate is the at most one multiple of
//!    `m` in the settled enclosure of its critical value, decided the same
//!    way.
//!
//! Every emitted polynomial is finally re-checked with the exact
//! real-rootedness test, and a failure is reported as
//! [`RealRootedEnumerationError::CertificationFailure`] rather than being
//! silently dropped.
//!
//! # Resource limits
//!
//! [`RealRootedEnumerationOptions`] bounds the number of examined tower nodes,
//! the number of returned polynomials and the bit length of tower
//! coefficients.  Hitting a limit returns the results found so far with
//! [`RealRootedEnumerationStatus::ResourceLimit`]; an empty
//! [`RealRootedEnumerationStatus::Complete`] result means the set is
//! certified empty.
//!
//! # Uses and limits of interpretation
//!
//! Typical uses are exhaustive searches for characteristic or generating
//! polynomials constrained by degree, trace (`a_1`), the quadratic
//! coefficient `a_2` (for example `tr A^2` of a symmetric matrix), symmetry,
//! or divisibility conditions that a caller filters afterwards.  A complete
//! enumeration is a certified finite computation for the given prefix; it is
//! not by itself a proof of any statement about an infinite family.
//!
//! The paper's linear-programming enumeration of interlacers is not part of
//! this module: its published form relies on numerical approximations.
//!
//! Related modules: [`crate::root_count`] (exact root counting and
//! isolation), [`crate::real_rootedness`] (discriminants, square-free parts),
//! [`crate::numeric`] (factorials).
//!
//! # Example
//!
//! ```
//! use num_bigint::BigInt;
//! use polytool::real_rooted_enumeration::{
//!     enumerate_real_rooted_integer_polynomials_bigint, RealRootedEnumerationOptions,
//!     RealRootedEnumerationStatus,
//! };
//!
//! // All real-rooted x^3 + 0 x^2 - 3 x + a_3 with integer a_3.
//! let prefix = [1, 0, -3].map(BigInt::from);
//! let result = enumerate_real_rooted_integer_polynomials_bigint(
//!     3,
//!     &prefix,
//!     &RealRootedEnumerationOptions::default(),
//! )
//! .unwrap();
//! assert_eq!(result.status, RealRootedEnumerationStatus::Complete);
//! // x^3 - 3x + c is real-rooted exactly for -2 <= c <= 2 (ascending output).
//! let constants: Vec<i64> = result
//!     .polynomials
//!     .iter()
//!     .map(|p| i64::try_from(&p[0]).unwrap())
//!     .collect();
//! assert_eq!(constants, vec![-2, -1, 0, 1, 2]);
//! assert_eq!(result.polynomials[0], [-2, -3, 0, 1].map(BigInt::from));
//! ```

use crate::numeric::factorial_bigint;
use crate::real_rootedness::squarefree_part_bigint_coeffs;
use crate::root_count::{
    is_real_rooted_fast_bigint_coeffs, isolate_real_roots_bigint_coeffs, DyadicRational,
    IsolatedRealRoot,
};
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};
use std::cmp::Ordering;
use std::fmt;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// A nonempty inclusive range `lower ..= upper` of integers.
///
/// Ranges can be astronomically large, so they are not materialized; use
/// [`InclusiveIntegerRange::iter`] for lazy iteration.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct InclusiveIntegerRange {
    /// Smallest element.
    pub lower: BigInt,
    /// Largest element (`lower <= upper`).
    pub upper: BigInt,
}

impl InclusiveIntegerRange {
    /// Number of integers in the range, `upper - lower + 1`.
    pub fn count(&self) -> BigInt {
        &self.upper - &self.lower + 1
    }

    /// Whether `value` lies in the range.
    pub fn contains(&self, value: &BigInt) -> bool {
        &self.lower <= value && value <= &self.upper
    }

    /// Lazily iterate over the range in increasing order.
    pub fn iter(&self) -> impl Iterator<Item = BigInt> + '_ {
        std::iter::successors(Some(self.lower.clone()), move |value| {
            (value < &self.upper).then(|| value + 1)
        })
    }
}

/// Error from [`antiderivative_exact_bigint`]: a coefficient is not divisible
/// by its new exponent, so the zero-constant antiderivative is not integral.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AntiderivativeDivisionError {
    /// The exponent `k + 1` that does not divide the coefficient of `x^k`.
    pub exponent: usize,
    /// The coefficient of `x^k` in the input (ascending index `k = exponent - 1`).
    pub coefficient: BigInt,
}

impl fmt::Display for AntiderivativeDivisionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "antiderivative is not integral: coefficient {} of x^{} is not divisible by {}",
            self.coefficient,
            self.exponent - 1,
            self.exponent
        )
    }
}

impl std::error::Error for AntiderivativeDivisionError {}

/// Unsupported input for the antiderivative-shift kernel.
///
/// These are never confused with an empty answer, which is `Ok(None)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AntiderivativeIntervalError {
    /// The derivative is the zero polynomial, so `P + C = C` has no degree.
    ZeroPolynomial,
    /// The derivative has degree 0 or 1, so the admissible set is all of `R`
    /// or a half-line and cannot be returned as a bounded range.
    UnboundedDegree {
        /// Degree of the supplied derivative.
        degree: usize,
    },
    /// The modulus of [`real_rooted_antiderivative_constant_multiples_bigint`]
    /// must be positive.
    NonPositiveModulus,
}

impl fmt::Display for AntiderivativeIntervalError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroPolynomial => write!(f, "the derivative is the zero polynomial"),
            Self::UnboundedDegree { degree } => write!(
                f,
                "derivative of degree {degree} gives an unbounded set of real-rooted shifts \
                 (degree at least 2 is required)"
            ),
            Self::NonPositiveModulus => write!(f, "the modulus must be positive"),
        }
    }
}

impl std::error::Error for AntiderivativeIntervalError {}

/// Resource budgets for [`enumerate_real_rooted_integer_polynomials_bigint`].
///
/// `None` means unlimited.  The default has no limits; adapters that expose
/// the search to untrusted input should set conservative limits.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RealRootedEnumerationOptions {
    /// Maximum number of tower nodes whose admissible-constant range is
    /// computed (see [`RealRootedEnumerationStats::nodes_examined`]).
    pub max_nodes: Option<u64>,
    /// Maximum number of polynomials returned.  Finding one more stops the
    /// search with [`RealRootedEnumerationLimit::Results`]; a search with
    /// exactly this many results is still complete.
    pub max_results: Option<usize>,
    /// Maximum bit length (of the absolute value) of any coefficient of any
    /// derivative-tower polynomial `p_i`, including the factorial scaling
    /// `(n-i)! a_i` of its constant term and the start polynomial
    /// `p_(t-1)`.
    pub max_coefficient_bits: Option<u64>,
}

/// Which budget stopped an incomplete enumeration.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RealRootedEnumerationLimit {
    /// [`RealRootedEnumerationOptions::max_nodes`].
    Nodes,
    /// [`RealRootedEnumerationOptions::max_results`].
    Results,
    /// [`RealRootedEnumerationOptions::max_coefficient_bits`].
    CoefficientBits,
}

impl RealRootedEnumerationLimit {
    /// Stable machine name: `max_nodes`, `max_results` or
    /// `max_coefficient_bits`.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Nodes => "max_nodes",
            Self::Results => "max_results",
            Self::CoefficientBits => "max_coefficient_bits",
        }
    }
}

/// Whether the enumeration is complete.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RealRootedEnumerationStatus {
    /// The returned list is the complete set (possibly empty).
    Complete,
    /// A budget was exhausted; the returned list is the prefix (in output
    /// order) of the complete set found so far.
    ResourceLimit(RealRootedEnumerationLimit),
}

impl RealRootedEnumerationStatus {
    /// Stable machine name: `complete` or `resource_limit`.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Complete => "complete",
            Self::ResourceLimit(_) => "resource_limit",
        }
    }
}

/// Deterministic search counters (no wall-clock times).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RealRootedEnumerationStats {
    /// Number of levels at which a coefficient is chosen, `n - t + 1`.
    pub derivative_levels: usize,
    /// Tower nodes `p_i` with `t-1 <= i < n` whose admissible-constant range
    /// was computed by the kernel.
    pub nodes_examined: u64,
    /// Examined nodes without any admissible constant in the required
    /// factorial divisibility class (dead ends).
    pub nodes_pruned: u64,
    /// Exact real-root isolations performed (one per examined node).
    pub exact_root_isolations: u64,
    /// Bisection steps used to refine isolating intervals until the
    /// critical-value enclosures were narrow enough.
    pub refinement_steps: u64,
    /// Exact real-rootedness tests of undecided endpoint candidates, plus the
    /// test of the start polynomial `p_(t-1)`.
    pub endpoint_tests: u64,
    /// Exact real-rootedness re-checks of emitted polynomials.
    pub final_verifications: u64,
    /// Number of polynomials returned.
    pub results_emitted: usize,
    /// `polynomials_by_level[j]` counts the real-rooted tower polynomials of
    /// degree `t - 1 + j` that were generated (entry 0 is the start
    /// polynomial; the last entry counts emitted results).
    pub polynomials_by_level: Vec<u64>,
}

/// Output of [`enumerate_real_rooted_integer_polynomials_bigint`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RealRootedEnumerationResult {
    /// The real-rooted polynomials, each as **ascending** coefficients
    /// `[a_n, ..., a_1, a_0]` of length `n + 1`, in increasing lexicographic
    /// order of `(a_t, ..., a_n)`.
    pub polynomials: Vec<Vec<BigInt>>,
    /// Whether the list is complete.
    pub status: RealRootedEnumerationStatus,
    /// Search counters.
    pub stats: RealRootedEnumerationStats,
}

/// Invalid input or an internal certification failure.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RealRootedEnumerationError {
    /// The degree must be at least 2.
    DegreeTooSmall {
        /// The requested degree.
        degree: usize,
    },
    /// The prefix length `t` must satisfy `3 <= t <= degree + 1`.
    PrefixLength {
        /// Number of supplied leading coefficients.
        prefix_len: usize,
        /// The requested degree.
        degree: usize,
    },
    /// The leading coefficient `a_0` is zero.
    ZeroLeadingCoefficient,
    /// An internal invariant failed (a tower antiderivative was not integral,
    /// a tower polynomial was not real-rooted, an emitted polynomial failed
    /// the final exact check, or the output order was violated).  This
    /// indicates a bug, never a property of the input.
    CertificationFailure {
        /// Description of the failed invariant.
        message: String,
    },
}

impl fmt::Display for RealRootedEnumerationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DegreeTooSmall { degree } => {
                write!(f, "degree must be at least 2, got {degree}")
            }
            Self::PrefixLength { prefix_len, degree } => write!(
                f,
                "expected between 3 and {} leading coefficients for degree {degree}, got {prefix_len}",
                degree + 1
            ),
            Self::ZeroLeadingCoefficient => write!(f, "the leading coefficient a_0 must be nonzero"),
            Self::CertificationFailure { message } => {
                write!(f, "internal certification failure: {message}")
            }
        }
    }
}

impl std::error::Error for RealRootedEnumerationError {}

// ---------------------------------------------------------------------------
// Antiderivatives
// ---------------------------------------------------------------------------

fn degree_of(p: &[BigInt]) -> Option<usize> {
    p.iter().rposition(|c| !c.is_zero())
}

/// Zero-constant antiderivative `∫_0^x q(y) dy` with integer coefficients.
///
/// Coefficients are ascending.  Every division by an exponent is checked:
/// if the coefficient of `x^k` is not divisible by `k + 1` the error names
/// it instead of truncating.  Trailing zeros are removed and the zero
/// polynomial is returned as `[0]`, matching [`crate::vec_poly`].
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::real_rooted_enumeration::antiderivative_exact_bigint;
///
/// let q = [0, 6, 3].map(BigInt::from); // 6x + 3x^2
/// assert_eq!(antiderivative_exact_bigint(&q).unwrap(), [0, 0, 3, 1].map(BigInt::from));
/// assert!(antiderivative_exact_bigint(&[BigInt::from(1), BigInt::from(1)]).is_err());
/// ```
pub fn antiderivative_exact_bigint(
    coeffs: &[BigInt],
) -> Result<Vec<BigInt>, AntiderivativeDivisionError> {
    let Some(d) = degree_of(coeffs) else {
        return Ok(vec![BigInt::zero()]);
    };
    let mut out = Vec::with_capacity(d + 2);
    out.push(BigInt::zero());
    for (k, c) in coeffs[..=d].iter().enumerate() {
        let exponent = BigInt::from(k + 1);
        let (quotient, remainder) = c.div_rem(&exponent);
        if !remainder.is_zero() {
            return Err(AntiderivativeDivisionError {
                exponent: k + 1,
                coefficient: c.clone(),
            });
        }
        out.push(quotient);
    }
    Ok(out)
}

/// Zero-constant antiderivative `∫_0^x q(y) dy` with exact rational
/// coefficients (ascending; the zero polynomial gives `[0]`).
///
/// ```
/// use num_bigint::BigInt;
/// use num_rational::BigRational;
/// use polytool::real_rooted_enumeration::antiderivative_rational_bigint;
///
/// let p = antiderivative_rational_bigint(&[BigInt::from(1), BigInt::from(1)]);
/// assert_eq!(p[2], BigRational::new(1.into(), 2.into()));
/// ```
pub fn antiderivative_rational_bigint(coeffs: &[BigInt]) -> Vec<BigRational> {
    let Some(d) = degree_of(coeffs) else {
        return vec![BigRational::zero()];
    };
    let mut out = Vec::with_capacity(d + 2);
    out.push(BigRational::zero());
    for (k, c) in coeffs[..=d].iter().enumerate() {
        out.push(BigRational::new(c.clone(), BigInt::from(k + 1)));
    }
    out
}

// ---------------------------------------------------------------------------
// Exact helpers
// ---------------------------------------------------------------------------

fn derivative(p: &[BigInt]) -> Vec<BigInt> {
    p.iter()
        .enumerate()
        .skip(1)
        .map(|(k, c)| c * BigInt::from(k))
        .collect()
}

/// Exact division `a / b` over `Z[x]`, or `None` if it is not exact.
fn exact_div(a: &[BigInt], b: &[BigInt]) -> Option<Vec<BigInt>> {
    let db = degree_of(b)?;
    let Some(da) = degree_of(a) else {
        return Some(vec![BigInt::zero()]);
    };
    if da < db {
        return None;
    }
    let mut rem: Vec<BigInt> = a[..=da].to_vec();
    let mut quotient = vec![BigInt::zero(); da - db + 1];
    for k in (0..=da - db).rev() {
        let (qk, r) = rem[k + db].div_rem(&b[db]);
        if !r.is_zero() {
            return None;
        }
        for (j, bj) in b[..=db].iter().enumerate() {
            rem[j + k] -= &qk * bj;
        }
        quotient[k] = qk;
    }
    rem.iter().all(Zero::is_zero).then_some(quotient)
}

/// The exact dyadic number `num / 2^exp` (not necessarily in lowest terms).
///
/// All kernel arithmetic stays in this set, so it needs only integer
/// multiplications and shifts.
#[derive(Clone, Debug)]
struct Dyadic {
    num: BigInt,
    exp: usize,
}

impl Dyadic {
    fn integer(value: BigInt) -> Self {
        Dyadic { num: value, exp: 0 }
    }

    fn from_rational(value: &DyadicRational) -> Self {
        Dyadic {
            num: value.numerator().clone(),
            exp: value.exponent() as usize,
        }
    }

    fn num_at(&self, exp: usize) -> BigInt {
        &self.num << (exp - self.exp)
    }

    fn add(&self, other: &Self) -> Self {
        let exp = self.exp.max(other.exp);
        Dyadic {
            num: self.num_at(exp) + other.num_at(exp),
            exp,
        }
    }

    fn sub(&self, other: &Self) -> Self {
        let exp = self.exp.max(other.exp);
        Dyadic {
            num: self.num_at(exp) - other.num_at(exp),
            exp,
        }
    }

    fn neg(&self) -> Self {
        Dyadic {
            num: -&self.num,
            exp: self.exp,
        }
    }

    fn sign(&self) -> Ordering {
        self.num.cmp(&BigInt::zero())
    }

    /// `floor(self / divisor)` for a positive integer divisor.
    fn floor_div(&self, divisor: &BigInt) -> BigInt {
        self.num.div_floor(&(divisor << self.exp))
    }

    /// `ceil(self / divisor)` for a positive integer divisor.
    fn ceil_div(&self, divisor: &BigInt) -> BigInt {
        -(-&self.num).div_floor(&(divisor << self.exp))
    }

    /// `self < value` for an integer `value`.
    fn lt_integer(&self, value: &BigInt) -> bool {
        self.num < (value << self.exp)
    }
}

impl PartialEq for Dyadic {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Dyadic {}

impl PartialOrd for Dyadic {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Dyadic {
    fn cmp(&self, other: &Self) -> Ordering {
        let exp = self.exp.max(other.exp);
        self.num_at(exp).cmp(&other.num_at(exp))
    }
}

/// `f(num / 2^exp)` exactly, for ascending integer coefficients `f`.
///
/// Horner's rule on the homogenized polynomial: the result has exponent
/// `exp * deg f`.
fn eval_dyadic(f: &[BigInt], num: &BigInt, exp: usize) -> Dyadic {
    let Some(d) = degree_of(f) else {
        return Dyadic::integer(BigInt::zero());
    };
    let mut acc = f[d].clone();
    for (k, c) in f[..d].iter().enumerate().rev() {
        acc = acc * num + (c << (exp * (d - k)));
    }
    Dyadic {
        num: acc,
        exp: exp * d,
    }
}

// ---------------------------------------------------------------------------
// The antiderivative-shift kernel
// ---------------------------------------------------------------------------

/// Counters collected by the kernel.
#[derive(Default)]
struct KernelCounters {
    isolations: u64,
    refinements: u64,
    endpoint_tests: u64,
}

enum KernelOutcome {
    /// The derivative itself is not real-rooted.
    NotRealRooted,
    /// No admissible multiple.
    Empty,
    /// Admissible multipliers `lower..=upper`.
    Range(BigInt, BigInt),
}

/// Slope class of a critical point of `P`, from Lemma 2.5.
#[derive(Clone, Copy, PartialEq, Eq)]
enum CriticalKind {
    /// `q'(µ) < 0`: requires `C >= -P(µ)`.
    Lower,
    /// `q'(µ) > 0`: requires `C <= -P(µ)`.
    Upper,
    /// Multiple root: requires `C = -P(µ)`.
    Forced,
}

enum RootLocation {
    /// The root is `num / 2^exp`.
    Exact { num: BigInt, exp: usize },
    /// Exactly one root of the square-free part lies in
    /// `(lower / 2^exp, upper / 2^exp)`; neither endpoint is a root.
    Interval {
        lower: BigInt,
        upper: BigInt,
        exp: usize,
        /// Sign of the square-free part at the lower endpoint (never zero).
        lower_sign: Ordering,
    },
}

/// A critical point `µ` of `P` with an exact enclosure of `-D P(µ)`.
struct CriticalPoint {
    kind: CriticalKind,
    location: RootLocation,
    lower: Dyadic,
    upper: Dyadic,
    /// The enclosure is shorter than `D m` or contains no multiple of it.
    settled: bool,
}

struct ShiftProblem {
    /// Square-free part of the derivative `q`.
    squarefree: Vec<BigInt>,
    /// `D P`, with `P = ∫_0^x q` and `D` its common denominator.
    scaled_antiderivative: Vec<BigInt>,
    /// `D` times an integer bound for `|q'|` on an interval containing every
    /// real root of `q`.
    slope_bound: BigInt,
    /// `D m`: enclosures of `-D P(µ)` are compared with multiples of this.
    scaled_modulus: BigInt,
}

impl ShiftProblem {
    /// Enclosure `[lower, upper]` of `-D P(µ)`.
    ///
    /// For an isolating interval `[l, u]` with midpoint `c` and width `w`,
    /// `q(µ) = 0` gives `|q(x)| <= |x - µ| max |q'|`, hence
    /// `|P(µ) - P(c)| <= max |q'| · w^2 / 8`.
    fn enclosure(&self, location: &RootLocation) -> (Dyadic, Dyadic) {
        match location {
            RootLocation::Exact { num, exp } => {
                let value = eval_dyadic(&self.scaled_antiderivative, num, *exp).neg();
                (value.clone(), value)
            }
            RootLocation::Interval {
                lower, upper, exp, ..
            } => {
                let center =
                    eval_dyadic(&self.scaled_antiderivative, &(lower + upper), exp + 1).neg();
                let width = upper - lower;
                let radius = Dyadic {
                    num: &self.slope_bound * &width * &width,
                    exp: 2 * exp + 3,
                };
                (center.sub(&radius), center.add(&radius))
            }
        }
    }

    /// Whether an enclosure is narrower than `D m` or contains no multiple
    /// of it.
    fn is_settled(&self, lower: &Dyadic, upper: &Dyadic) -> bool {
        upper.sub(lower).lt_integer(&self.scaled_modulus)
            || upper.floor_div(&self.scaled_modulus) < lower.ceil_div(&self.scaled_modulus)
    }

    fn update(&self, point: &mut CriticalPoint) {
        let (lower, upper) = self.enclosure(&point.location);
        point.settled = self.is_settled(&lower, &upper);
        point.lower = lower;
        point.upper = upper;
    }

    /// Bisect the isolating interval once and update the enclosure.
    fn refine(&self, point: &mut CriticalPoint, counters: &mut KernelCounters) {
        counters.refinements += 1;
        if let RootLocation::Interval {
            lower,
            upper,
            exp,
            lower_sign,
        } = &mut point.location
        {
            let mid = &*lower + &*upper;
            let mid_exp = *exp + 1;
            match eval_dyadic(&self.squarefree, &mid, mid_exp).sign() {
                Ordering::Equal => {
                    point.location = RootLocation::Exact {
                        num: mid,
                        exp: mid_exp,
                    };
                }
                sign if sign == *lower_sign => {
                    *lower = mid;
                    *upper <<= 1;
                    *exp = mid_exp;
                }
                _ => {
                    *lower <<= 1;
                    *upper = mid;
                    *exp = mid_exp;
                }
            }
        }
        self.update(point);
    }

    /// Exact test: is `P + multiple · m` real-rooted?
    fn is_admissible(&self, multiple: &BigInt, counters: &mut KernelCounters) -> bool {
        counters.endpoint_tests += 1;
        let mut shifted = self.scaled_antiderivative.clone();
        shifted[0] += &self.scaled_modulus * multiple;
        is_real_rooted_fast_bigint_coeffs(&shifted)
    }
}

/// Kernel for a derivative `q` with positive leading coefficient and degree
/// at least 2, and a positive modulus.
fn admissible_multiples_positive(
    q: &[BigInt],
    modulus: &BigInt,
    counters: &mut KernelCounters,
) -> KernelOutcome {
    let squarefree = squarefree_part_bigint_coeffs(q);
    let squarefree_degree = degree_of(&squarefree).unwrap_or(0);
    counters.isolations += 1;
    let roots = isolate_real_roots_bigint_coeffs(q);
    if roots.len() != squarefree_degree {
        return KernelOutcome::NotRealRooted;
    }

    // D P with integer coefficients; D = 1 along the derivative tower.
    let rational = antiderivative_rational_bigint(q);
    let denominator = rational
        .iter()
        .fold(BigInt::one(), |acc, c| acc.lcm(c.denom()));
    let scaled_antiderivative: Vec<BigInt> = rational
        .iter()
        .map(|c| c.numer() * (&denominator / c.denom()))
        .collect();

    // Every root lies in [-M, M] for an integer M bounding the endpoints.
    let root_bound = roots
        .iter()
        .map(|root| {
            let lower = root.lower();
            let upper = root.upper();
            let a = (lower.numerator().abs() >> lower.exponent()) + 1u32;
            let b = (upper.numerator().abs() >> upper.exponent()) + 1u32;
            a.max(b)
        })
        .max()
        .unwrap_or_else(BigInt::one);
    let q_prime = derivative(q);
    let abs_q_prime: Vec<BigInt> = q_prime.iter().map(Signed::abs).collect();
    let slope_bound = eval_dyadic(&abs_q_prime, &root_bound, 0).num * &denominator;

    let problem = ShiftProblem {
        squarefree,
        scaled_antiderivative,
        scaled_modulus: &denominator * modulus,
        slope_bound,
    };

    // Multiple roots of q are exactly the roots of g = q / squarefree.
    let multiple_part = exact_div(q, &problem.squarefree)
        .expect("the square-free part divides the polynomial over Z");
    let multiple_squarefree = if degree_of(&multiple_part).unwrap_or(0) > 0 {
        Some(squarefree_part_bigint_coeffs(&multiple_part))
    } else {
        None
    };

    let mut points: Vec<CriticalPoint> = Vec::with_capacity(roots.len());
    for root in &roots {
        let (kind, location) = match root {
            IsolatedRealRoot::Exact(r) => {
                let r = Dyadic::from_rational(r);
                let kind = match eval_dyadic(&q_prime, &r.num, r.exp).sign() {
                    Ordering::Less => CriticalKind::Lower,
                    Ordering::Greater => CriticalKind::Upper,
                    Ordering::Equal => CriticalKind::Forced,
                };
                (
                    kind,
                    RootLocation::Exact {
                        num: r.num,
                        exp: r.exp,
                    },
                )
            }
            IsolatedRealRoot::Interval { lower, upper } => {
                let exp = lower.exponent().max(upper.exponent()) as usize;
                let lower = Dyadic::from_rational(lower).num_at(exp);
                let upper = Dyadic::from_rational(upper).num_at(exp);
                let multiple = multiple_squarefree.as_ref().is_some_and(|g| {
                    eval_dyadic(g, &lower, exp).sign() != eval_dyadic(g, &upper, exp).sign()
                });
                let kind = if multiple {
                    CriticalKind::Forced
                } else if eval_dyadic(q, &upper, exp).sign() == Ordering::Less {
                    CriticalKind::Lower
                } else {
                    CriticalKind::Upper
                };
                let lower_sign = eval_dyadic(&problem.squarefree, &lower, exp).sign();
                (
                    kind,
                    RootLocation::Interval {
                        lower,
                        upper,
                        exp,
                        lower_sign,
                    },
                )
            }
        };
        let mut point = CriticalPoint {
            kind,
            location,
            lower: Dyadic::integer(BigInt::zero()),
            upper: Dyadic::integer(BigInt::zero()),
            settled: false,
        };
        problem.update(&mut point);
        points.push(point);
    }

    let dm = &problem.scaled_modulus;

    // With a multiple root the admissible set has at most one point.
    if let Some(point) = points
        .iter_mut()
        .find(|point| point.kind == CriticalKind::Forced)
    {
        while !point.settled {
            problem.refine(point, counters);
        }
        let mut k = point.lower.ceil_div(dm);
        let last = point.upper.floor_div(dm);
        while k <= last {
            if problem.is_admissible(&k, counters) {
                return KernelOutcome::Range(k.clone(), k);
            }
            k += 1;
        }
        return KernelOutcome::Empty;
    }

    // Simple roots only: [L, U] with L = max over Lower, U = min over Upper.
    // Refine only the points that can still attain the maximum (minimum).
    loop {
        let lower_max = points
            .iter()
            .filter(|p| p.kind == CriticalKind::Lower)
            .map(|p| p.lower.clone())
            .max()
            .expect("a real-rooted q of degree >= 2 with simple roots has a falling root");
        let upper_min = points
            .iter()
            .filter(|p| p.kind == CriticalKind::Upper)
            .map(|p| p.upper.clone())
            .min()
            .expect("a real-rooted q of degree >= 2 with simple roots has a rising root");
        let mut refined = false;
        for point in &mut points {
            let relevant = match point.kind {
                CriticalKind::Lower => point.upper >= lower_max,
                CriticalKind::Upper => point.lower <= upper_min,
                CriticalKind::Forced => false,
            };
            if relevant && !point.settled {
                problem.refine(point, counters);
                refined = true;
            }
        }
        if !refined {
            break;
        }
    }
    let lower_group = || points.iter().filter(|p| p.kind == CriticalKind::Lower);
    let upper_group = || points.iter().filter(|p| p.kind == CriticalKind::Upper);
    let l_lo = lower_group().map(|p| &p.lower).max().expect("nonempty");
    let l_hi = lower_group().map(|p| &p.upper).max().expect("nonempty");
    let u_lo = upper_group().map(|p| &p.lower).min().expect("nonempty");
    let u_hi = upper_group().map(|p| &p.upper).min().expect("nonempty");

    // Multiples k (of D m) with l_hi <= k D m <= u_lo are certainly
    // admissible; at most one undecided multiple remains at each end, since
    // the enclosure attaining the extreme is settled.
    let mut first = l_lo.ceil_div(dm);
    let mut last = u_hi.floor_div(dm);
    if first > last {
        return KernelOutcome::Empty;
    }
    let mut tested: Option<(BigInt, bool)> = None;
    if Dyadic::integer(&first * dm) < *l_hi {
        let ok = problem.is_admissible(&first, counters);
        tested = Some((first.clone(), ok));
        if !ok {
            first += 1;
        }
    }
    if first > last {
        return KernelOutcome::Empty;
    }
    if Dyadic::integer(&last * dm) > *u_lo {
        let ok = match &tested {
            Some((k, ok)) if k == &last => *ok,
            _ => problem.is_admissible(&last, counters),
        };
        if !ok {
            last -= 1;
        }
    }
    if first > last {
        KernelOutcome::Empty
    } else {
        KernelOutcome::Range(first, last)
    }
}

/// Validated kernel shared by the public entry points.
fn admissible_multiples(
    derivative_coefficients: &[BigInt],
    modulus: &BigInt,
    counters: &mut KernelCounters,
) -> Result<KernelOutcome, AntiderivativeIntervalError> {
    if !modulus.is_positive() {
        return Err(AntiderivativeIntervalError::NonPositiveModulus);
    }
    let Some(d) = degree_of(derivative_coefficients) else {
        return Err(AntiderivativeIntervalError::ZeroPolynomial);
    };
    if d < 2 {
        return Err(AntiderivativeIntervalError::UnboundedDegree { degree: d });
    }
    let q = &derivative_coefficients[..=d];
    if q[d].is_positive() {
        return Ok(admissible_multiples_positive(q, modulus, counters));
    }
    // P + C is real-rooted iff -P - C is; -P is the antiderivative of -q.
    let negated: Vec<BigInt> = q.iter().map(|c| -c).collect();
    Ok(
        match admissible_multiples_positive(&negated, modulus, counters) {
            KernelOutcome::Range(lo, hi) => KernelOutcome::Range(-hi, -lo),
            other => other,
        },
    )
}

/// Integers `C` for which `P(x) + C` is real-rooted, `P(x) = ∫_0^x q(y) dy`.
///
/// `derivative_coefficients` is `q` in ascending order (trailing zeros are
/// ignored); it may have any sign and need not have an integral
/// antiderivative.  The answer is `Ok(Some(range))` for a nonempty set,
/// `Ok(None)` when no integer works (including when `q` itself is not
/// real-rooted), and an error for unsupported input: the zero polynomial or
/// `deg q < 2`, where the set is unbounded.
///
/// The computation is exact; see the [module documentation](self) for how
/// the endpoints are certified.  Repeated roots of `q` force a single
/// candidate (Lemma 2.5), and coincident critical values need no special
/// treatment.
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::real_rooted_enumeration::real_rooted_antiderivative_integer_constants_bigint;
///
/// // q = 3x^2 - 3, P = x^3 - 3x: P + C is real-rooted iff -2 <= C <= 2.
/// let q = [-3, 0, 3].map(BigInt::from);
/// let range = real_rooted_antiderivative_integer_constants_bigint(&q).unwrap().unwrap();
/// assert_eq!((range.lower, range.upper), (BigInt::from(-2), BigInt::from(2)));
///
/// // q = x^2 + 1 is not real-rooted, so no shift is.
/// let q = [1, 0, 1].map(BigInt::from);
/// assert_eq!(real_rooted_antiderivative_integer_constants_bigint(&q), Ok(None));
/// ```
pub fn real_rooted_antiderivative_integer_constants_bigint(
    derivative_coefficients: &[BigInt],
) -> Result<Option<InclusiveIntegerRange>, AntiderivativeIntervalError> {
    real_rooted_antiderivative_constant_multiples_bigint(derivative_coefficients, &BigInt::one())
}

/// Multipliers `k` for which `P(x) + k·modulus` is real-rooted, where
/// `P(x) = ∫_0^x q(y) dy`.
///
/// This is [`real_rooted_antiderivative_integer_constants_bigint`] restricted
/// to a divisibility class; the enumerator uses `modulus = (n-i)!`.  The
/// returned range contains the multipliers `k`, not the constants `k·modulus`.
/// Refining only to the width of `modulus` makes large moduli cheaper than
/// filtering the full integer range.
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::real_rooted_enumeration::real_rooted_antiderivative_constant_multiples_bigint;
///
/// // P = x^3 - 3x, constants in [-2, 2] divisible by 2: k = -1, 0, 1.
/// let q = [-3, 0, 3].map(BigInt::from);
/// let range = real_rooted_antiderivative_constant_multiples_bigint(&q, &BigInt::from(2))
///     .unwrap()
///     .unwrap();
/// assert_eq!((range.lower, range.upper), (BigInt::from(-1), BigInt::from(1)));
/// ```
pub fn real_rooted_antiderivative_constant_multiples_bigint(
    derivative_coefficients: &[BigInt],
    modulus: &BigInt,
) -> Result<Option<InclusiveIntegerRange>, AntiderivativeIntervalError> {
    let mut counters = KernelCounters::default();
    Ok(
        match admissible_multiples(derivative_coefficients, modulus, &mut counters)? {
            KernelOutcome::Range(lower, upper) => Some(InclusiveIntegerRange { lower, upper }),
            KernelOutcome::Empty | KernelOutcome::NotRealRooted => None,
        },
    )
}

// ---------------------------------------------------------------------------
// The derivative-tower enumerator
// ---------------------------------------------------------------------------

/// Ascending coefficients of `p_i = d^(n-i) p / dx^(n-i)` from the descending
/// prefix `a_0, ..., a_i`: the coefficient of `x^k` is `a_(i-k) (n-i+k)!/k!`.
fn tower_start(degree: usize, prefix: &[BigInt]) -> Vec<BigInt> {
    let i = prefix.len() - 1;
    (0..=i)
        .map(|k| {
            let falling = ((k + 1)..=(degree - i + k))
                .fold(BigInt::one(), |acc, factor| acc * BigInt::from(factor));
            &prefix[i - k] * falling
        })
        .collect()
}

struct Frame {
    /// Degree `i` of the children `p_i`.
    level: usize,
    /// `∫_0^x p_(i-1)` with zero constant term.
    base: Vec<BigInt>,
    /// `(n-i)!`.
    scale: BigInt,
    next: BigInt,
    last: BigInt,
}

fn exceeds_bits(value: &BigInt, limit: Option<u64>) -> bool {
    limit.is_some_and(|bits| value.bits() > bits)
}

fn finish(
    polynomials: Vec<Vec<BigInt>>,
    status: RealRootedEnumerationStatus,
    mut stats: RealRootedEnumerationStats,
) -> Result<RealRootedEnumerationResult, RealRootedEnumerationError> {
    stats.results_emitted = polynomials.len();
    Ok(RealRootedEnumerationResult {
        polynomials,
        status,
        stats,
    })
}

fn certification_failure(message: String) -> RealRootedEnumerationError {
    RealRootedEnumerationError::CertificationFailure { message }
}

/// Enumerate all real-rooted integer polynomials of degree `degree` with the
/// prescribed leading coefficients.
///
/// `leading_coefficients_descending` is `[a_0, a_1, ..., a_(t-1)]` for
/// `p(x) = a_0 x^n + a_1 x^(n-1) + ... + a_n`, with `a_0 != 0` (of either
/// sign) and `3 <= t <= n + 1`.  The result lists every real-rooted
/// completion, reducible and repeated-root ones included, as **ascending**
/// coefficient vectors in increasing lexicographic order of
/// `(a_t, ..., a_n)`.  With `t = n + 1` the result is the input polynomial
/// itself if it is real-rooted and empty otherwise.
///
/// All decisions are exact (see the [module documentation](self)).  When a
/// budget in `options` is exhausted the status is
/// [`RealRootedEnumerationStatus::ResourceLimit`] and the list holds the
/// results found so far, which are an initial segment of the complete output.
///
/// # Errors
///
/// Invalid degree, prefix length or a zero leading coefficient; and
/// [`RealRootedEnumerationError::CertificationFailure`] if an internal
/// invariant fails (a bug, not a property of the input).
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::real_rooted_enumeration::{
///     enumerate_real_rooted_integer_polynomials_bigint, RealRootedEnumerationOptions,
/// };
///
/// // Real-rooted x^4 - 4x^2 + a_3 x + a_4: includes (x^2 - 2)^2 and x^2 (x^2 - 4).
/// let prefix = [1, 0, -4].map(BigInt::from);
/// let result = enumerate_real_rooted_integer_polynomials_bigint(
///     4,
///     &prefix,
///     &RealRootedEnumerationOptions::default(),
/// )
/// .unwrap();
/// assert!(result.polynomials.contains(&[4, 0, -4, 0, 1].map(BigInt::from).to_vec()));
/// assert!(result.polynomials.contains(&[0, 0, -4, 0, 1].map(BigInt::from).to_vec()));
/// ```
pub fn enumerate_real_rooted_integer_polynomials_bigint(
    degree: usize,
    leading_coefficients_descending: &[BigInt],
    options: &RealRootedEnumerationOptions,
) -> Result<RealRootedEnumerationResult, RealRootedEnumerationError> {
    let n = degree;
    let prefix = leading_coefficients_descending;
    if n < 2 {
        return Err(RealRootedEnumerationError::DegreeTooSmall { degree: n });
    }
    let t = prefix.len();
    if t < 3 || t > n + 1 {
        return Err(RealRootedEnumerationError::PrefixLength {
            prefix_len: t,
            degree: n,
        });
    }
    if prefix[0].is_zero() {
        return Err(RealRootedEnumerationError::ZeroLeadingCoefficient);
    }

    let mut stats = RealRootedEnumerationStats {
        derivative_levels: n + 1 - t,
        polynomials_by_level: vec![0; n + 2 - t],
        ..RealRootedEnumerationStats::default()
    };
    let mut polynomials: Vec<Vec<BigInt>> = Vec::new();

    let start = tower_start(n, prefix);
    if start
        .iter()
        .any(|c| exceeds_bits(c, options.max_coefficient_bits))
    {
        return finish(
            polynomials,
            RealRootedEnumerationStatus::ResourceLimit(RealRootedEnumerationLimit::CoefficientBits),
            stats,
        );
    }
    stats.endpoint_tests += 1;
    if !is_real_rooted_fast_bigint_coeffs(&start) {
        return finish(polynomials, RealRootedEnumerationStatus::Complete, stats);
    }
    stats.polynomials_by_level[0] = 1;

    let mut counters = KernelCounters::default();
    let mut stack: Vec<Frame> = Vec::new();
    let mut pending: Option<Vec<BigInt>> = Some(start);
    let mut status = RealRootedEnumerationStatus::Complete;

    'search: loop {
        // Expand a pending node p_(i-1), or emit it when it is p_n.
        if let Some(node) = pending.take() {
            let node_degree = node.len() - 1;
            if node_degree == n {
                if options.max_results == Some(polynomials.len()) {
                    status = RealRootedEnumerationStatus::ResourceLimit(
                        RealRootedEnumerationLimit::Results,
                    );
                    break 'search;
                }
                stats.final_verifications += 1;
                if !is_real_rooted_fast_bigint_coeffs(&node) {
                    return Err(certification_failure(format!(
                        "emitted polynomial {node:?} is not real-rooted"
                    )));
                }
                if let Some(previous) = polynomials.last() {
                    let increasing =
                        node.iter().rev().cmp(previous.iter().rev()) == Ordering::Greater;
                    if !increasing {
                        return Err(certification_failure(format!(
                            "output order violated at {node:?}"
                        )));
                    }
                }
                polynomials.push(node);
            } else {
                if options.max_nodes == Some(stats.nodes_examined) {
                    status = RealRootedEnumerationStatus::ResourceLimit(
                        RealRootedEnumerationLimit::Nodes,
                    );
                    break 'search;
                }
                stats.nodes_examined += 1;
                let level = node_degree + 1;
                let scale = factorial_bigint(n - level);
                let outcome =
                    admissible_multiples(&node, &scale, &mut counters).map_err(|error| {
                        certification_failure(format!("kernel rejected {node:?}: {error}"))
                    })?;
                match outcome {
                    KernelOutcome::NotRealRooted => {
                        return Err(certification_failure(format!(
                            "tower polynomial {node:?} is not real-rooted"
                        )));
                    }
                    KernelOutcome::Empty => stats.nodes_pruned += 1,
                    KernelOutcome::Range(first, last) => {
                        if exceeds_bits(&(&first * &scale), options.max_coefficient_bits)
                            || exceeds_bits(&(&last * &scale), options.max_coefficient_bits)
                        {
                            status = RealRootedEnumerationStatus::ResourceLimit(
                                RealRootedEnumerationLimit::CoefficientBits,
                            );
                            break 'search;
                        }
                        let base = antiderivative_exact_bigint(&node).map_err(|error| {
                            certification_failure(format!(
                                "tower antiderivative of {node:?}: {error}"
                            ))
                        })?;
                        stack.push(Frame {
                            level,
                            base,
                            scale,
                            next: first,
                            last,
                        });
                    }
                }
            }
        }

        // Produce the next child of the deepest open frame.
        loop {
            let Some(frame) = stack.last_mut() else {
                break 'search;
            };
            if frame.next > frame.last {
                stack.pop();
                continue;
            }
            let mut child = frame.base.clone();
            child[0] = &frame.next * &frame.scale;
            frame.next += 1;
            stats.polynomials_by_level[frame.level + 1 - t] += 1;
            pending = Some(child);
            break;
        }
    }

    // A result-limit stop counted a child that was not emitted.
    if status == RealRootedEnumerationStatus::ResourceLimit(RealRootedEnumerationLimit::Results) {
        stats.polynomials_by_level[n + 1 - t] -= 1;
    }
    stats.exact_root_isolations = counters.isolations;
    stats.refinement_steps = counters.refinements;
    stats.endpoint_tests += counters.endpoint_tests;
    finish(polynomials, status, stats)
}

#[cfg(test)]
mod tests;
