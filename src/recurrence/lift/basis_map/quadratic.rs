//! Exact real quadratic numbers `a + b sqrt(d)` and intervals with such
//! endpoints, for the source and target intervals of the basis-map theorems
//! (for example the sharp radius `1 + sqrt(2)` of the A035607 map).
//!
//! Every comparison is exact: the sign of `a + b sqrt(d)` is decided by
//! comparing `a^2` with `b^2 d`, and a polynomial with rational coefficients
//! is evaluated in `Q(sqrt(d))`.  Rational enclosures are computed from the
//! integer square root and are only used for display and for rational outer
//! (or inner) approximations of an interval; they never decide anything.

use crate::polynomial::Polynomial;
use crate::real_rootedness::{isolate_real_roots_rational_coeffs, RealInterval};
use crate::recurrence::{format_rational_coeff, parse_rational_coeff, RecurrenceError};
use crate::recurrence::{RootEndpoint, RootInterval};
use crate::vec_poly;
use crate::BigRational;
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::fmt;

type Q = BigRational;

/// Number of bits of the rational enclosure recorded in the JSON form.
const ENCLOSURE_BITS: u32 = 32;

/// An exact real number `a + b sqrt(d)` with rationals `a`, `b` and `d > 0`.
///
/// The representation is normalized: either `b = 0` (a rational, stored with
/// `d = 1`) or `d` is a positive rational that is not the square of a
/// rational, so that the representation is unique up to the choice of `d`
/// (two numbers with different radicands are still compared exactly).  The
/// JSON form records the three exact parts as decimal strings, the readable
/// `text` and a rational enclosure `[lower_bound, upper_bound]` of width at
/// most `2^-32`; on input only the exact parts are read.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(try_from = "QuadraticNumberJson", into = "QuadraticNumberJson")]
pub struct QuadraticNumber {
    rational_part: Q,
    sqrt_coefficient: Q,
    radicand: Q,
}

/// The JSON form of a [`QuadraticNumber`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct QuadraticNumberJson {
    rational_part: String,
    sqrt_coefficient: String,
    radicand: String,
    #[serde(default)]
    text: String,
    #[serde(default)]
    lower_bound: String,
    #[serde(default)]
    upper_bound: String,
}

impl TryFrom<QuadraticNumberJson> for QuadraticNumber {
    type Error = RecurrenceError;

    fn try_from(json: QuadraticNumberJson) -> Result<Self, Self::Error> {
        Self::new(
            parse_rational_coeff(&json.rational_part)?,
            parse_rational_coeff(&json.sqrt_coefficient)?,
            parse_rational_coeff(&json.radicand)?,
        )
    }
}

impl From<QuadraticNumber> for QuadraticNumberJson {
    fn from(value: QuadraticNumber) -> Self {
        let (lower, upper) = value.enclosure(ENCLOSURE_BITS);
        Self {
            rational_part: format_rational_coeff(&value.rational_part),
            sqrt_coefficient: format_rational_coeff(&value.sqrt_coefficient),
            radicand: format_rational_coeff(&value.radicand),
            text: value.to_string(),
            lower_bound: format_rational_coeff(&lower),
            upper_bound: format_rational_coeff(&upper),
        }
    }
}

/// The rational square root of `value >= 0`, when it has one.
fn rational_square_root(value: &Q) -> Option<Q> {
    if value.is_negative() {
        return None;
    }
    let root = |n: &BigInt| {
        let r = n.sqrt();
        (&r * &r == *n).then_some(r)
    };
    Some(Q::new(root(value.numer())?, root(value.denom())?))
}

/// `floor(sqrt(value) * 2^bits) / 2^bits` for `value >= 0`.
fn square_root_floor(value: &Q, bits: u32) -> Q {
    // sqrt(p/q) = sqrt(p q) / q.
    let scale = BigInt::one() << bits;
    let radicand = value.numer() * value.denom() * &scale * &scale;
    Q::new(radicand.sqrt(), value.denom() * scale)
}

impl QuadraticNumber {
    /// The number `a + b sqrt(d)`, normalized.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::InvalidInput`] when `d < 0` and `b != 0`.
    pub fn new(a: Q, b: Q, d: Q) -> Result<Self, RecurrenceError> {
        if b.is_zero() {
            return Ok(Self::rational(a));
        }
        if d.is_negative() {
            return Err(RecurrenceError::InvalidInput(format!(
                "the radicand {} of a quadratic number is negative",
                format_rational_coeff(&d)
            )));
        }
        if let Some(root) = rational_square_root(&d) {
            return Ok(Self::rational(a + b * root));
        }
        Ok(Self {
            rational_part: a,
            sqrt_coefficient: b,
            radicand: d,
        })
    }

    /// The rational number `a`.
    pub fn rational(a: Q) -> Self {
        Self {
            rational_part: a,
            sqrt_coefficient: Q::zero(),
            radicand: Q::one(),
        }
    }

    /// `sqrt(d)` for a rational `d >= 0`.
    ///
    /// # Errors
    ///
    /// As [`QuadraticNumber::new`].
    pub fn square_root(d: Q) -> Result<Self, RecurrenceError> {
        Self::new(Q::zero(), Q::one(), d)
    }

    /// The rational part `a`.
    pub fn rational_part(&self) -> &Q {
        &self.rational_part
    }

    /// The coefficient `b` of `sqrt(d)` (zero for a rational).
    pub fn sqrt_coefficient(&self) -> &Q {
        &self.sqrt_coefficient
    }

    /// The radicand `d` (one for a rational).
    pub fn radicand(&self) -> &Q {
        &self.radicand
    }

    /// The value, when it is rational.
    pub fn as_rational(&self) -> Option<&Q> {
        self.sqrt_coefficient
            .is_zero()
            .then_some(&self.rational_part)
    }

    /// Whether the number is rational.
    pub fn is_rational(&self) -> bool {
        self.sqrt_coefficient.is_zero()
    }

    /// The exact sign.
    pub fn signum(&self) -> Ordering {
        sign_of(&self.rational_part, &self.sqrt_coefficient, &self.radicand)
    }

    /// `-self`.
    pub fn negated(&self) -> Self {
        Self {
            rational_part: -self.rational_part.clone(),
            sqrt_coefficient: -self.sqrt_coefficient.clone(),
            radicand: self.radicand.clone(),
        }
    }

    /// `self + c` for a rational `c`.
    pub fn add_rational(&self, c: &Q) -> Self {
        Self {
            rational_part: &self.rational_part + c,
            ..self.clone()
        }
    }

    /// `c self` for a rational `c`.
    pub fn scale(&self, c: &Q) -> Self {
        if c.is_zero() {
            return Self::rational(Q::zero());
        }
        Self {
            rational_part: &self.rational_part * c,
            sqrt_coefficient: &self.sqrt_coefficient * c,
            radicand: self.radicand.clone(),
        }
    }

    /// `self + other` when both lie in the same field `Q(sqrt(d))` (one of
    /// them rational, or the same radicand); `None` otherwise.
    pub fn checked_add(&self, other: &Self) -> Option<Self> {
        if other.is_rational() {
            return Some(self.add_rational(&other.rational_part));
        }
        if self.is_rational() {
            return Some(other.add_rational(&self.rational_part));
        }
        (self.radicand == other.radicand).then(|| {
            Self::new(
                &self.rational_part + &other.rational_part,
                &self.sqrt_coefficient + &other.sqrt_coefficient,
                self.radicand.clone(),
            )
            .expect("same nonnegative radicand")
        })
    }

    /// Exact comparison with a rational.
    pub fn cmp_rational(&self, value: &Q) -> Ordering {
        sign_of(
            &(&self.rational_part - value),
            &self.sqrt_coefficient,
            &self.radicand,
        )
    }

    /// Exact comparison with another quadratic number (any radicands).
    pub fn cmp_exact(&self, other: &Self) -> Ordering {
        if let Some(sum) = self.checked_add(&other.negated()) {
            return sum.signum();
        }
        // a1 + b1 sqrt(d1) - (a2 + b2 sqrt(d2)) = A - B with
        // A = (a1 - a2) + b1 sqrt(d1), B = b2 sqrt(d2).
        let left = Self {
            rational_part: &self.rational_part - &other.rational_part,
            sqrt_coefficient: self.sqrt_coefficient.clone(),
            radicand: self.radicand.clone(),
        };
        let right_sign = other.sqrt_coefficient.cmp(&Q::zero());
        let left_sign = left.signum();
        if left_sign == Ordering::Equal {
            return right_sign.reverse();
        }
        if left_sign != right_sign {
            return left_sign;
        }
        // Same sign: compare |A| with |B| through A^2 - B^2.
        let left_square = left.square();
        let right_square = &other.sqrt_coefficient * &other.sqrt_coefficient * &other.radicand;
        let difference = left_square.add_rational(&-right_square).signum();
        if left_sign == Ordering::Greater {
            difference
        } else {
            difference.reverse()
        }
    }

    /// `self^2`.
    fn square(&self) -> Self {
        let a = &self.rational_part;
        let b = &self.sqrt_coefficient;
        Self::new(
            a * a + b * b * &self.radicand,
            Q::from_integer(BigInt::from(2)) * a * b,
            self.radicand.clone(),
        )
        .expect("nonnegative radicand")
    }

    /// The value of the polynomial `p` (ascending rational coefficients) at
    /// this number, exactly.
    pub fn evaluate_polynomial(&self, p: &[Q]) -> Self {
        let mut a = Q::zero();
        let mut b = Q::zero();
        for coefficient in p.iter().rev() {
            // (a + b r)(x + y r) with r^2 = d, then + coefficient.
            let x = &self.rational_part;
            let y = &self.sqrt_coefficient;
            let next_a = &a * x + &b * y * &self.radicand + coefficient;
            let next_b = &a * y + &b * x;
            a = next_a;
            b = next_b;
        }
        Self::new(a, b, self.radicand.clone()).expect("nonnegative radicand")
    }

    /// A rational enclosure `lower <= self <= upper` of width at most
    /// `2^(1 - bits) |b|` (exact for a rational).
    pub fn enclosure(&self, bits: u32) -> (Q, Q) {
        if self.is_rational() {
            return (self.rational_part.clone(), self.rational_part.clone());
        }
        let floor = square_root_floor(&self.radicand, bits);
        let ceil = &floor + Q::new(BigInt::one(), BigInt::one() << bits);
        let low = &self.rational_part + &self.sqrt_coefficient * &floor;
        let high = &self.rational_part + &self.sqrt_coefficient * &ceil;
        if self.sqrt_coefficient.is_positive() {
            (low, high)
        } else {
            (high, low)
        }
    }

    /// The monic minimal polynomial over the rationals (ascending): `x - a`
    /// for a rational, else `(x - a)^2 - b^2 d`.
    pub fn minimal_polynomial(&self) -> Vec<Q> {
        if self.is_rational() {
            return vec![-self.rational_part.clone(), Q::one()];
        }
        let a = &self.rational_part;
        let b = &self.sqrt_coefficient;
        vec![
            a * a - b * b * &self.radicand,
            Q::from_integer(BigInt::from(-2)) * a,
            Q::one(),
        ]
    }
}

/// The exact sign of `a + b sqrt(d)`, `d >= 0`.
fn sign_of(a: &Q, b: &Q, d: &Q) -> Ordering {
    let zero = Q::zero();
    let sa = a.cmp(&zero);
    let sb = if d.is_zero() {
        Ordering::Equal
    } else {
        b.cmp(&zero)
    };
    if sb == Ordering::Equal {
        return sa;
    }
    if sa == Ordering::Equal || sa == sb {
        return sb;
    }
    // Opposite signs: compare a^2 with b^2 d.
    match (a * a).cmp(&(b * b * d)) {
        Ordering::Greater => sa,
        Ordering::Less => sb,
        Ordering::Equal => Ordering::Equal,
    }
}

impl fmt::Display for QuadraticNumber {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_rational() {
            return write!(f, "{}", format_rational_coeff(&self.rational_part));
        }
        let root = format!("sqrt({})", format_rational_coeff(&self.radicand));
        let b = &self.sqrt_coefficient;
        let magnitude = b.abs();
        let term = if magnitude.is_one() {
            root
        } else {
            format!("{} {root}", format_rational_coeff(&magnitude))
        };
        if self.rational_part.is_zero() {
            if b.is_negative() {
                write!(f, "-{term}")
            } else {
                write!(f, "{term}")
            }
        } else {
            let sign = if b.is_negative() { '-' } else { '+' };
            write!(
                f,
                "{} {sign} {term}",
                format_rational_coeff(&self.rational_part)
            )
        }
    }
}

/// One end of an [`AlgebraicInterval`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AlgebraicEndpoint {
    /// The exact value.
    pub value: QuadraticNumber,
    /// Whether the endpoint belongs to the interval.
    pub closed: bool,
}

/// An interval of the real line whose finite endpoints are
/// [`QuadraticNumber`]s; an omitted end is infinite.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AlgebraicInterval {
    /// Lower end; `None` means `-infinity`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lower: Option<AlgebraicEndpoint>,
    /// Upper end; `None` means `+infinity`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub upper: Option<AlgebraicEndpoint>,
}

fn closed(value: QuadraticNumber) -> Option<AlgebraicEndpoint> {
    Some(AlgebraicEndpoint {
        value,
        closed: true,
    })
}

impl AlgebraicInterval {
    /// The closed interval `[lower, upper]`.
    pub fn closed(lower: QuadraticNumber, upper: QuadraticNumber) -> Self {
        Self {
            lower: closed(lower),
            upper: closed(upper),
        }
    }

    /// The closed ray `(-infinity, upper]`.
    pub fn at_most(upper: QuadraticNumber) -> Self {
        Self {
            lower: None,
            upper: closed(upper),
        }
    }

    /// The closed ray `[lower, +infinity)`.
    pub fn at_least(lower: QuadraticNumber) -> Self {
        Self {
            lower: closed(lower),
            upper: None,
        }
    }

    /// The whole real line.
    pub fn whole_line() -> Self {
        Self {
            lower: None,
            upper: None,
        }
    }

    /// The closed ray `(-infinity, 0]`.
    pub fn nonpositive_ray() -> Self {
        Self::at_most(QuadraticNumber::rational(Q::zero()))
    }

    /// The interval with the rational endpoints of `interval`.
    pub fn from_real(interval: &RealInterval) -> Self {
        let end = |value: &Option<Q>, closed: bool| {
            value.as_ref().map(|v| AlgebraicEndpoint {
                value: QuadraticNumber::rational(v.clone()),
                closed,
            })
        };
        Self {
            lower: end(&interval.lower, interval.lower_closed),
            upper: end(&interval.upper, interval.upper_closed),
        }
    }

    /// Whether the interval has at most one point (decided exactly).
    pub fn is_degenerate(&self) -> bool {
        match (&self.lower, &self.upper) {
            (Some(a), Some(b)) => a.value.cmp_exact(&b.value) != Ordering::Less,
            _ => false,
        }
    }

    /// Whether the exact number `x` lies in the interval.
    pub fn contains(&self, x: &QuadraticNumber) -> bool {
        let above = self
            .lower
            .as_ref()
            .is_none_or(|a| match x.cmp_exact(&a.value) {
                Ordering::Greater => true,
                Ordering::Equal => a.closed,
                Ordering::Less => false,
            });
        let below = self
            .upper
            .as_ref()
            .is_none_or(|b| match x.cmp_exact(&b.value) {
                Ordering::Less => true,
                Ordering::Equal => b.closed,
                Ordering::Greater => false,
            });
        above && below
    }

    /// Whether every point of `self` lies in `other` (decided exactly).
    pub fn is_subset_of(&self, other: &Self) -> bool {
        let lower_ok = match (&self.lower, &other.lower) {
            (_, None) => true,
            (None, Some(_)) => false,
            (Some(a), Some(b)) => match a.value.cmp_exact(&b.value) {
                Ordering::Greater => true,
                Ordering::Equal => b.closed || !a.closed,
                Ordering::Less => false,
            },
        };
        let upper_ok = match (&self.upper, &other.upper) {
            (_, None) => true,
            (None, Some(_)) => false,
            (Some(a), Some(b)) => match a.value.cmp_exact(&b.value) {
                Ordering::Less => true,
                Ordering::Equal => b.closed || !a.closed,
                Ordering::Greater => false,
            },
        };
        lower_ok && upper_ok
    }

    /// The intersection with `other` (exact; possibly empty).
    pub fn intersection(&self, other: &Self) -> Self {
        let lower = match (&self.lower, &other.lower) {
            (None, b) => b.clone(),
            (a, None) => a.clone(),
            (Some(a), Some(b)) => match a.value.cmp_exact(&b.value) {
                Ordering::Greater => Some(a.clone()),
                Ordering::Less => Some(b.clone()),
                Ordering::Equal => Some(AlgebraicEndpoint {
                    value: a.value.clone(),
                    closed: a.closed && b.closed,
                }),
            },
        };
        let upper = match (&self.upper, &other.upper) {
            (None, b) => b.clone(),
            (a, None) => a.clone(),
            (Some(a), Some(b)) => match a.value.cmp_exact(&b.value) {
                Ordering::Less => Some(a.clone()),
                Ordering::Greater => Some(b.clone()),
                Ordering::Equal => Some(AlgebraicEndpoint {
                    value: a.value.clone(),
                    closed: a.closed && b.closed,
                }),
            },
        };
        Self { lower, upper }
    }

    /// Whether `self` lies in `(-infinity, 0]`.
    pub fn is_nonpositive(&self) -> bool {
        self.is_subset_of(&Self::nonpositive_ray())
    }

    /// A closed interval with rational endpoints containing `self`
    /// (irrational ends rounded outward to `bits` binary digits).
    pub fn rational_outer(&self, bits: u32) -> RealInterval {
        let lower = self.lower.as_ref().map(|e| e.value.enclosure(bits).0);
        let upper = self.upper.as_ref().map(|e| e.value.enclosure(bits).1);
        RealInterval {
            lower_closed: lower.is_some(),
            lower,
            upper_closed: upper.is_some(),
            upper,
        }
    }

    /// A closed interval with rational endpoints contained in `self`
    /// (irrational ends rounded inward to `bits` binary digits; rational
    /// ends kept, also when open, so the result is closed and may touch an
    /// open end of `self`: callers use it only where that end is closed).
    pub fn rational_inner(&self, bits: u32) -> RealInterval {
        let lower = self.lower.as_ref().map(|e| e.value.enclosure(bits).1);
        let upper = self.upper.as_ref().map(|e| e.value.enclosure(bits).0);
        RealInterval {
            lower_closed: lower.is_some(),
            lower,
            upper_closed: upper.is_some(),
            upper,
        }
    }

    /// The [`RootInterval`] of a rational outer approximation
    /// ([`AlgebraicInterval::rational_outer`] with 32 bits).
    pub fn to_root_interval(&self) -> RootInterval {
        let outer = self.rational_outer(ENCLOSURE_BITS);
        let end = |value: &Option<Q>, closed: bool| {
            value.as_ref().map(|v| RootEndpoint {
                value: format_rational_coeff(v),
                per_index: None,
                closed,
            })
        };
        let lower_closed = self
            .lower
            .as_ref()
            .is_none_or(|e| e.closed || !e.value.is_rational());
        let upper_closed = self
            .upper
            .as_ref()
            .is_none_or(|e| e.closed || !e.value.is_rational());
        RootInterval {
            lower: end(&outer.lower, lower_closed),
            upper: end(&outer.upper, upper_closed),
        }
    }

    /// Whether the nonzero polynomial `p` is real-rooted with every root in
    /// the interval (decided exactly; vacuous for nonzero constants, false
    /// for the zero polynomial).
    pub fn contains_all_roots_of(&self, p: &[Q]) -> bool {
        let p = vec_poly::trim(p);
        if vec_poly::is_zero(&p) {
            return false;
        }
        let isolation = isolate_real_roots_rational_coeffs(&p);
        if !isolation.is_real_rooted() {
            return false;
        }
        let squarefree = Polynomial::new(p.clone()).squarefree_part();
        let squarefree = squarefree.coeffs().to_vec();
        isolation.roots.iter().all(|root| {
            let compare = |end: &QuadraticNumber| -> Ordering {
                // Position of the root relative to `end`.
                if root.lower == root.upper {
                    return end.cmp_rational(&root.lower).reverse();
                }
                if end.cmp_rational(&root.lower) != Ordering::Greater {
                    return Ordering::Greater;
                }
                if end.cmp_rational(&root.upper) != Ordering::Less {
                    return Ordering::Less;
                }
                let at_end = end.evaluate_polynomial(&squarefree).signum();
                if at_end == Ordering::Equal {
                    return Ordering::Equal;
                }
                let at_lower = vec_poly::evaluate(&squarefree, &root.lower).cmp(&Q::zero());
                if at_lower != at_end {
                    Ordering::Less
                } else {
                    Ordering::Greater
                }
            };
            let above = self.lower.as_ref().is_none_or(|a| match compare(&a.value) {
                Ordering::Greater => true,
                Ordering::Equal => a.closed,
                Ordering::Less => false,
            });
            let below = self.upper.as_ref().is_none_or(|b| match compare(&b.value) {
                Ordering::Less => true,
                Ordering::Equal => b.closed,
                Ordering::Greater => false,
            });
            above && below
        })
    }
}

impl fmt::Display for AlgebraicInterval {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.lower {
            None => write!(f, "(-inf")?,
            Some(a) => write!(f, "{}{}", if a.closed { '[' } else { '(' }, a.value)?,
        }
        match &self.upper {
            None => write!(f, ", +inf)"),
            Some(b) => write!(f, ", {}{}", b.value, if b.closed { ']' } else { ')' }),
        }
    }
}

/// The integer `value` as a rational.
pub(super) fn integer(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(a: i64, b: i64) -> Q {
        Q::new(BigInt::from(a), BigInt::from(b))
    }

    #[test]
    fn signs_and_comparisons_are_exact() {
        let r = QuadraticNumber::new(q(1, 1), q(1, 1), q(2, 1)).unwrap();
        assert_eq!(r.to_string(), "1 + sqrt(2)");
        assert_eq!(r.cmp_rational(&q(12, 5)), Ordering::Greater);
        assert_eq!(r.cmp_rational(&q(5, 2)), Ordering::Less);
        // 3 - 2 sqrt(2) > 0, 1 - sqrt(2) < 0.
        assert_eq!(
            QuadraticNumber::new(q(3, 1), q(-2, 1), q(2, 1))
                .unwrap()
                .signum(),
            Ordering::Greater
        );
        assert_eq!(
            QuadraticNumber::new(q(1, 1), q(-1, 1), q(2, 1))
                .unwrap()
                .signum(),
            Ordering::Less
        );
        // A rational square folds: 1 + sqrt(9/4) = 5/2.
        let folded = QuadraticNumber::new(q(1, 1), q(1, 1), q(9, 4)).unwrap();
        assert_eq!(folded.as_rational(), Some(&q(5, 2)));
        // sqrt(2) + sqrt(3) versus 3: different radicands.
        let two = QuadraticNumber::square_root(q(2, 1)).unwrap();
        let three = QuadraticNumber::square_root(q(3, 1))
            .unwrap()
            .add_rational(&q(-3, 1));
        // sqrt(2) vs 3 - sqrt(3) = 1.268: greater.
        assert_eq!(two.cmp_exact(&three.negated()), Ordering::Greater);
        let (lo, hi) = r.enclosure(20);
        assert!(lo < hi && r.cmp_rational(&lo) == Ordering::Greater);
        assert_eq!(r.cmp_rational(&hi), Ordering::Less);
        // Minimal polynomial of 1 + sqrt(2): x^2 - 2x - 1.
        assert_eq!(r.minimal_polynomial(), vec![q(-1, 1), q(-2, 1), q(1, 1)]);
        assert!(r.evaluate_polynomial(&r.minimal_polynomial()).signum() == Ordering::Equal);
    }

    #[test]
    fn roots_are_located_against_irrational_endpoints() {
        let r = QuadraticNumber::new(q(-1, 1), q(-1, 1), q(2, 1)).unwrap();
        let interval = AlgebraicInterval::closed(r.clone(), QuadraticNumber::rational(q(0, 1)));
        // (u + 2)^2 u: roots -2, 0 inside [-1 - sqrt(2), 0].
        assert!(interval.contains_all_roots_of(&[q(0, 1), q(4, 1), q(4, 1), q(1, 1)]));
        // u + 5/2 has its root -5/2 < -2.414.
        assert!(!interval.contains_all_roots_of(&[q(5, 2), q(1, 1)]));
        // u^2 + 2u - 1 has the roots -1 +- sqrt(2): the endpoint itself and a positive root.
        assert!(!interval.contains_all_roots_of(&[q(-1, 1), q(2, 1), q(1, 1)]));
        let wider = AlgebraicInterval::closed(r, QuadraticNumber::rational(q(1, 1)));
        assert!(wider.contains_all_roots_of(&[q(-1, 1), q(2, 1), q(1, 1)]));
        let json = serde_json::to_string(&wider).unwrap();
        let back: AlgebraicInterval = serde_json::from_str(&json).unwrap();
        assert_eq!(back, wider);
    }
}
