//! Exact real-root enclosures with multiplicities, and exact interval tests.
//!
//! [`isolate_real_roots_exact`] isolates the distinct real roots of a
//! polynomial with rational coefficients by rational (dyadic) intervals, using
//! the exact Sturm bisection of [`exact_real_roots::isolate_real_roots_bigint_coeffs`],
//! and attaches to every root its exact multiplicity from Yun's square-free
//! factorization over `Q`.  The result also records the degree, the valuation
//! (multiplicity of the root `0`) and the number of real roots counted with
//! multiplicity, so real-rootedness is decided exactly.
//!
//! [`RealInterval`] is an interval of the real line with rational or infinite
//! endpoints, each open or closed (rays are intervals with one infinite end).
//! [`RealRootIsolation::all_real_roots_in`] and
//! [`RealRootIsolation::real_rooted_in`] decide exactly whether every real root
//! (respectively every root, all of them real) lies in such an interval; roots
//! on an endpoint are located exactly, never by a floating-point comparison.
//!
//! Conventions: coefficients are in ascending degree order.  The zero
//! polynomial has no degree, no roots, and is treated as real-rooted (as in
//! [`crate::is_real_rooted_bigint_coeffs`]); interval predicates hold for it
//! vacuously.  A nonzero constant is real-rooted with no roots.
//!
//! An enclosure with `lower == upper` is an exact rational root.  Otherwise
//! exactly one distinct root lies in the open interval `(lower, upper)` and
//! neither endpoint is a root; neighbouring enclosures never overlap.  The
//! rational midpoint of an enclosure is only a representative: it is neither
//! a root nor a certified bound (contrast [`super::real_roots`], which returns
//! such midpoints).

use super::Q;
use crate::polynomial::Polynomial;
use exact_real_roots::{isolate_real_roots_bigint_coeffs, DyadicRational, IsolatedRealRoot};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::cmp::Ordering;
use std::fmt;

/// One distinct real root, enclosed exactly, with its multiplicity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RealRootEnclosure {
    /// Lower endpoint.  Equal to `upper` exactly when the root is this rational.
    pub lower: Q,
    /// Upper endpoint.  When `lower < upper` neither endpoint is a root and the
    /// root lies in the open interval `(lower, upper)`.
    pub upper: Q,
    /// Multiplicity of the root (at least 1).
    pub multiplicity: usize,
}

impl RealRootEnclosure {
    /// Whether the root is the rational `lower == upper`.
    pub fn is_exact(&self) -> bool {
        self.lower == self.upper
    }
}

/// An interval of the real line with rational or infinite endpoints.
///
/// `lower = None` is `-infinity` and `upper = None` is `+infinity`; infinite
/// ends are always open, whatever the `*_closed` flag says.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RealInterval {
    /// Lower endpoint; `None` means `-infinity`.
    pub lower: Option<Q>,
    /// Whether a finite lower endpoint belongs to the interval.
    pub lower_closed: bool,
    /// Upper endpoint; `None` means `+infinity`.
    pub upper: Option<Q>,
    /// Whether a finite upper endpoint belongs to the interval.
    pub upper_closed: bool,
}

impl RealInterval {
    /// The closed interval `[lower, upper]`.
    pub fn closed(lower: Q, upper: Q) -> Self {
        Self {
            lower: Some(lower),
            lower_closed: true,
            upper: Some(upper),
            upper_closed: true,
        }
    }

    /// The closed ray `(-infinity, upper]`.
    pub fn at_most(upper: Q) -> Self {
        Self {
            lower: None,
            lower_closed: false,
            upper: Some(upper),
            upper_closed: true,
        }
    }

    /// The closed ray `[lower, +infinity)`.
    pub fn at_least(lower: Q) -> Self {
        Self {
            lower: Some(lower),
            lower_closed: true,
            upper: None,
            upper_closed: false,
        }
    }

    /// The whole real line.
    pub fn whole_line() -> Self {
        Self {
            lower: None,
            lower_closed: false,
            upper: None,
            upper_closed: false,
        }
    }

    /// The single point `[point, point]`.
    pub fn point(point: Q) -> Self {
        Self::closed(point.clone(), point)
    }

    /// Whether the interval contains no point.
    pub fn is_empty(&self) -> bool {
        match (&self.lower, &self.upper) {
            (Some(a), Some(b)) => match a.cmp(b) {
                Ordering::Greater => true,
                Ordering::Equal => !(self.lower_closed && self.upper_closed),
                Ordering::Less => false,
            },
            _ => false,
        }
    }

    /// Whether `x` lies in the interval.
    pub fn contains(&self, x: &Q) -> bool {
        let above = match &self.lower {
            None => true,
            Some(a) => x > a || (self.lower_closed && x == a),
        };
        let below = match &self.upper {
            None => true,
            Some(b) => x < b || (self.upper_closed && x == b),
        };
        above && below
    }

    /// Whether the interval is contained in `other` (an empty interval is
    /// contained in everything).
    pub fn is_subset_of(&self, other: &Self) -> bool {
        if self.is_empty() {
            return true;
        }
        let lower_ok = match (&other.lower, &self.lower) {
            (None, _) => true,
            (Some(_), None) => false,
            (Some(b), Some(a)) => a > b || (a == b && (other.lower_closed || !self.lower_closed)),
        };
        let upper_ok = match (&other.upper, &self.upper) {
            (None, _) => true,
            (Some(_), None) => false,
            (Some(b), Some(a)) => a < b || (a == b && (other.upper_closed || !self.upper_closed)),
        };
        lower_ok && upper_ok
    }

    /// The smallest interval containing both intervals.
    pub fn hull(&self, other: &Self) -> Self {
        if self.is_empty() {
            return other.clone();
        }
        if other.is_empty() {
            return self.clone();
        }
        let (lower, lower_closed) = match (&self.lower, &other.lower) {
            (None, _) | (_, None) => (None, false),
            (Some(a), Some(b)) => match a.cmp(b) {
                Ordering::Less => (Some(a.clone()), self.lower_closed),
                Ordering::Greater => (Some(b.clone()), other.lower_closed),
                Ordering::Equal => (Some(a.clone()), self.lower_closed || other.lower_closed),
            },
        };
        let (upper, upper_closed) = match (&self.upper, &other.upper) {
            (None, _) | (_, None) => (None, false),
            (Some(a), Some(b)) => match a.cmp(b) {
                Ordering::Greater => (Some(a.clone()), self.upper_closed),
                Ordering::Less => (Some(b.clone()), other.upper_closed),
                Ordering::Equal => (Some(a.clone()), self.upper_closed || other.upper_closed),
            },
        };
        Self {
            lower_closed: lower_closed && lower.is_some(),
            lower,
            upper_closed: upper_closed && upper.is_some(),
            upper,
        }
    }
}

impl fmt::Display for RealInterval {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.lower {
            None => write!(f, "(-inf")?,
            Some(a) => write!(f, "{}{a}", if self.lower_closed { '[' } else { '(' })?,
        }
        match &self.upper {
            None => write!(f, ", +inf)"),
            Some(b) => write!(f, ", {b}{}", if self.upper_closed { ']' } else { ')' }),
        }
    }
}

/// Exact isolation of the real roots of one polynomial.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RealRootIsolation {
    /// Degree; `None` for the zero polynomial.
    pub degree: Option<usize>,
    /// Multiplicity of the root `0` (the valuation); `0` for the zero polynomial.
    pub valuation: usize,
    /// Number of real roots counted with multiplicity.
    pub real_root_count: usize,
    /// The distinct real roots in increasing order.
    pub roots: Vec<RealRootEnclosure>,
    /// Monic square-free part (`[1]` for constants, empty for zero); used for
    /// exact comparisons of roots with rational points.
    squarefree: Vec<Q>,
}

impl RealRootIsolation {
    /// Whether the polynomial is zero.
    pub fn is_zero_polynomial(&self) -> bool {
        self.degree.is_none()
    }

    /// Whether every root is real (vacuously true for zero and constants).
    pub fn is_real_rooted(&self) -> bool {
        self.degree.is_none_or(|d| d == self.real_root_count)
    }

    /// Exact position of a distinct root relative to the rational `point`.
    pub fn compare_root(&self, root: &RealRootEnclosure, point: &Q) -> Ordering {
        if root.is_exact() {
            return root.lower.cmp(point);
        }
        if point <= &root.lower {
            return Ordering::Greater;
        }
        if point >= &root.upper {
            return Ordering::Less;
        }
        let at_point = sign(&evaluate(&self.squarefree, point));
        if at_point == 0 {
            return Ordering::Equal;
        }
        if sign(&evaluate(&self.squarefree, &root.lower)) != at_point {
            Ordering::Less
        } else {
            Ordering::Greater
        }
    }

    /// Whether the distinct root lies in `interval` (decided exactly).
    pub fn root_in(&self, root: &RealRootEnclosure, interval: &RealInterval) -> bool {
        let above = match &interval.lower {
            None => true,
            Some(a) => match self.compare_root(root, a) {
                Ordering::Greater => true,
                Ordering::Equal => interval.lower_closed,
                Ordering::Less => false,
            },
        };
        let below = match &interval.upper {
            None => true,
            Some(b) => match self.compare_root(root, b) {
                Ordering::Less => true,
                Ordering::Equal => interval.upper_closed,
                Ordering::Greater => false,
            },
        };
        above && below
    }

    /// Whether every real root lies in `interval` (no claim about non-real
    /// roots; vacuous for the zero polynomial and constants).
    pub fn all_real_roots_in(&self, interval: &RealInterval) -> bool {
        self.roots.iter().all(|root| self.root_in(root, interval))
    }

    /// Whether the polynomial is real-rooted with every root in `interval`
    /// (vacuous for the zero polynomial and nonzero constants).
    pub fn real_rooted_in(&self, interval: &RealInterval) -> bool {
        self.is_real_rooted() && self.all_real_roots_in(interval)
    }

    /// The smallest closed interval with rational endpoints containing every
    /// enclosure, or `None` without real roots.
    pub fn enclosure_hull(&self) -> Option<RealInterval> {
        let first = self.roots.first()?;
        let last = self.roots.last()?;
        Some(RealInterval::closed(
            first.lower.clone(),
            last.upper.clone(),
        ))
    }

    /// Shrink every non-exact enclosure by exact bisection until its width is
    /// at most `max_width`; a midpoint that is a root becomes an exact root.
    ///
    /// # Panics
    ///
    /// Panics if `max_width` is not positive.
    pub fn refine(&mut self, max_width: &Q) {
        assert!(max_width.is_positive(), "max_width must be positive");
        let two = Q::from_integer(BigInt::from(2));
        for root in &mut self.roots {
            let mut lower_sign = sign(&evaluate(&self.squarefree, &root.lower));
            while !root.is_exact() && &(&root.upper - &root.lower) > max_width {
                let middle = (&root.lower + &root.upper) / &two;
                let middle_sign = sign(&evaluate(&self.squarefree, &middle));
                if middle_sign == 0 {
                    root.lower = middle.clone();
                    root.upper = middle;
                } else if middle_sign != lower_sign {
                    root.upper = middle;
                } else {
                    root.lower = middle;
                    lower_sign = middle_sign;
                }
            }
        }
    }
}

fn sign(value: &Q) -> i8 {
    if value.is_positive() {
        1
    } else if value.is_negative() {
        -1
    } else {
        0
    }
}

fn evaluate(p: &[Q], x: &Q) -> Q {
    p.iter().rev().fold(Q::zero(), |acc, c| acc * x + c)
}

fn dyadic(value: &DyadicRational) -> Q {
    Q::new(value.numerator().clone(), value.denominator())
}

/// Primitive integer multiple of a rational polynomial (nonzero input).
fn integer_multiple(p: &[Q]) -> Vec<BigInt> {
    let scale = p.iter().fold(BigInt::one(), |acc, c| acc.lcm(c.denom()));
    p.iter()
        .map(|c| (c * Q::from_integer(scale.clone())).to_integer())
        .collect()
}

/// Yun's square-free factorization: `factors[i]` is the monic product of the
/// distinct irreducible factors of multiplicity exactly `i + 1`.
fn squarefree_factors(p: &Polynomial<Q>) -> Vec<Polynomial<Q>> {
    let derivative = p.derivative();
    let common = p.gcd(&derivative);
    let mut c = p.exact_div(&common);
    let mut d = derivative.exact_div(&common) - c.derivative();
    let mut factors = Vec::new();
    while c.degree().is_some_and(|degree| degree > 0) {
        let a = c.gcd(&d);
        c = c.exact_div(&a);
        d = d.exact_div(&a) - c.derivative();
        factors.push(a);
    }
    factors
}

/// Enclosures are refined to width at most `1/DEFAULT_WIDTH_DENOMINATOR` by default.
const DEFAULT_WIDTH_DENOMINATOR: i64 = 1024;

/// Isolate the real roots of a polynomial with rational coefficients
/// (ascending order) exactly, with multiplicities.
///
/// Every enclosure is refined by exact bisection to width at most `1/1024`
/// (use [`RealRootIsolation::refine`] for narrower ones).  The root `0`, and
/// a root equal to the simplest rational of its enclosure (for example an
/// integer or `-1/2`), are reported exactly; other rational roots may keep an
/// open enclosure.
///
/// ```
/// use num_bigint::BigInt;
/// use num_rational::Ratio;
/// use polytool::real_rootedness::{isolate_real_roots_exact, RealInterval};
///
/// let q = |n: i64| Ratio::from_integer(BigInt::from(n));
/// // x (x + 1)^2 (x^2 - 2): roots -sqrt 2, -1 (double), 0, sqrt 2.
/// let p = [q(0), q(-2), q(-4), q(-1), q(2), q(1)];
/// let iso = isolate_real_roots_exact(&p);
/// assert!(iso.is_real_rooted());
/// assert_eq!(iso.valuation, 1);
/// assert_eq!(iso.roots.len(), 4);
/// assert_eq!(iso.roots[1].multiplicity, 2);
/// assert!(iso.real_rooted_in(&RealInterval::closed(q(-2), q(2))));
/// assert!(!iso.real_rooted_in(&RealInterval::closed(q(-1), q(2))));
/// ```
pub fn isolate_real_roots_exact(coeffs: &[Q]) -> RealRootIsolation {
    let polynomial = Polynomial::new(coeffs.to_vec());
    let Some(degree) = polynomial.degree() else {
        return RealRootIsolation {
            degree: None,
            valuation: 0,
            real_root_count: 0,
            roots: Vec::new(),
            squarefree: Vec::new(),
        };
    };
    let valuation = polynomial
        .coeffs()
        .iter()
        .position(|c| !c.is_zero())
        .expect("nonzero polynomial");
    if degree == 0 {
        return RealRootIsolation {
            degree: Some(0),
            valuation,
            real_root_count: 0,
            roots: Vec::new(),
            squarefree: vec![Q::one()],
        };
    }
    let factors = squarefree_factors(&polynomial);
    let squarefree = factors
        .iter()
        .fold(Polynomial::one(), |acc, f| acc * f.clone());
    let isolated = isolate_real_roots_bigint_coeffs(&integer_multiple(polynomial.coeffs()));
    let mut roots = Vec::with_capacity(isolated.len());
    for root in isolated {
        let (mut lower, mut upper) = match &root {
            IsolatedRealRoot::Exact(r) => (dyadic(r), dyadic(r)),
            IsolatedRealRoot::Interval { lower, upper } => (dyadic(lower), dyadic(upper)),
        };
        // The root 0 (valuation > 0) is reported exactly.
        if valuation > 0 && lower.is_negative() && upper.is_positive() {
            lower = Q::zero();
            upper = Q::zero();
        }
        let multiplicity = factors
            .iter()
            .position(|factor| {
                if lower == upper {
                    factor.evaluate(&lower).is_zero()
                } else {
                    sign(&factor.evaluate(&lower)) != sign(&factor.evaluate(&upper))
                }
            })
            .map(|i| i + 1)
            .expect("every distinct root belongs to one square-free factor");
        roots.push(RealRootEnclosure {
            lower,
            upper,
            multiplicity,
        });
    }
    let mut isolation = RealRootIsolation {
        degree: Some(degree),
        valuation,
        real_root_count: roots.iter().map(|r| r.multiplicity).sum(),
        roots,
        squarefree: squarefree.coeffs().to_vec(),
    };
    isolation.refine(&Q::new(
        BigInt::one(),
        BigInt::from(DEFAULT_WIDTH_DENOMINATOR),
    ));
    // A simple rational root (such as -1 or -1/2) is reported exactly.
    for root in &mut isolation.roots {
        if !root.is_exact() {
            let candidate = simplest_rational_between(&root.lower, &root.upper);
            if evaluate(&isolation.squarefree, &candidate).is_zero() {
                root.lower = candidate.clone();
                root.upper = candidate;
            }
        }
    }
    isolation
}

/// The rational with the smallest denominator in `[lower, upper]` (`lower <= upper`),
/// by continued fractions.
fn simplest_rational_between(lower: &Q, upper: &Q) -> Q {
    let floor = lower.floor();
    if &floor == lower || floor.clone() + Q::one() <= *upper {
        return if &floor == lower {
            floor
        } else {
            floor + Q::one()
        };
    }
    // Same integer part, lower not an integer: recurse on the reciprocals.
    let inner = simplest_rational_between(&(upper - &floor).recip(), &(lower - &floor).recip());
    floor + inner.recip()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(n: i64) -> Q {
        Q::from_integer(BigInt::from(n))
    }

    fn r(n: i64, d: i64) -> Q {
        Q::new(BigInt::from(n), BigInt::from(d))
    }

    fn poly(coeffs: &[i64]) -> Vec<Q> {
        coeffs.iter().map(|&c| q(c)).collect()
    }

    #[test]
    fn zero_constant_and_valuation_conventions() {
        let zero = isolate_real_roots_exact(&[]);
        assert!(zero.is_zero_polynomial() && zero.is_real_rooted());
        assert!(zero.real_rooted_in(&RealInterval::point(q(5))));
        let monomial = isolate_real_roots_exact(&poly(&[0, 0, 3]));
        assert_eq!(monomial.degree, Some(2));
        assert_eq!(monomial.valuation, 2);
        assert_eq!(monomial.roots.len(), 1);
        assert_eq!(monomial.roots[0].multiplicity, 2);
        assert!(monomial.roots[0].is_exact());
        let c = isolate_real_roots_exact(&poly(&[7]));
        assert!(c.is_real_rooted() && c.roots.is_empty() && c.valuation == 0);
    }

    #[test]
    fn repeated_and_endpoint_roots() {
        // (x + 1)^3 x^2 (2x + 1): roots -1 (3), -1/2, 0 (2).
        let p = vec_mul(&[&[1, 3, 3, 1], &[0, 0, 1], &[1, 2]]);
        let iso = isolate_real_roots_exact(&poly(&p));
        assert!(iso.is_real_rooted());
        assert_eq!(iso.real_root_count, 6);
        let multiplicities: Vec<_> = iso.roots.iter().map(|r| r.multiplicity).collect();
        assert_eq!(multiplicities, vec![3, 1, 2]);
        assert!(iso.roots.iter().all(RealRootEnclosure::is_exact));
        assert_eq!(iso.roots[1].lower, r(-1, 2));
        assert!(iso.real_rooted_in(&RealInterval::closed(q(-1), q(0))));
        let half_open = RealInterval {
            lower: Some(q(-1)),
            lower_closed: false,
            upper: Some(q(0)),
            upper_closed: true,
        };
        assert!(!iso.real_rooted_in(&half_open));
        assert!(iso.real_rooted_in(&RealInterval::at_most(q(0))));
        assert!(!iso.real_rooted_in(&RealInterval::at_most(r(-1, 3))));
        let mut refined = iso.clone();
        refined.refine(&r(1, 1000));
        assert!(refined
            .roots
            .iter()
            .all(|e| &e.upper - &e.lower <= r(1, 1000)));
    }

    #[test]
    fn irrational_roots_compare_exactly_with_points() {
        // x^2 - 2: compare with 1.414 and 1.415.
        let iso = isolate_real_roots_exact(&poly(&[-2, 0, 1]));
        let root = &iso.roots[1];
        assert_eq!(iso.compare_root(root, &r(1414, 1000)), Ordering::Greater);
        assert_eq!(iso.compare_root(root, &r(1415, 1000)), Ordering::Less);
        assert!(!iso.all_real_roots_in(&RealInterval::closed(r(-1414, 1000), q(2))));
        // 1 + x + x^2: no real roots; real roots vacuously in any interval but
        // not real-rooted.
        let complex = isolate_real_roots_exact(&poly(&[1, 1, 1]));
        assert!(complex.all_real_roots_in(&RealInterval::point(q(0))));
        assert!(!complex.is_real_rooted());
    }

    #[test]
    fn simplest_rationals() {
        assert_eq!(
            simplest_rational_between(&r(-2049, 2048), &r(-4095, 4096)),
            q(-1)
        );
        assert_eq!(simplest_rational_between(&r(3, 10), &r(4, 10)), r(1, 3));
        assert_eq!(simplest_rational_between(&r(-4, 10), &r(-3, 10)), r(-1, 3));
        assert_eq!(simplest_rational_between(&r(5, 2), &r(5, 2)), r(5, 2));
    }

    #[test]
    fn interval_algebra() {
        let a = RealInterval::closed(q(-1), q(0));
        let ray = RealInterval::at_most(q(0));
        assert!(a.is_subset_of(&ray) && !ray.is_subset_of(&a));
        assert_eq!(a.hull(&RealInterval::point(q(2))).to_string(), "[-1, 2]");
        assert_eq!(ray.to_string(), "(-inf, 0]");
        assert!(RealInterval {
            lower: Some(q(0)),
            lower_closed: false,
            upper: Some(q(0)),
            upper_closed: true
        }
        .is_empty());
    }

    fn vec_mul(factors: &[&[i64]]) -> Vec<i64> {
        factors
            .iter()
            .fold(vec![1], |acc, f| crate::vec_poly::mul(&acc, f))
    }
}
