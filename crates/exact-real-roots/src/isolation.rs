//! Exact dyadic endpoints and Sturm bisection for distinct real roots.
//! Uses the primitive square-free part and signed Sturm chains.

use crate::arithmetic::{degree, sign_i8, squarefree_part_primitive_bigint_coeffs};
use crate::prs::{sign_variations, sturm_chain_from_squarefree};
use crate::uspensky::strict_fujiwara_root_bound;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};

/// An exact dyadic rational `numerator / 2^exponent`.
///
/// Values are kept in lowest terms: either `exponent == 0` or the numerator
/// is odd. Root isolation uses dyadic endpoints because bisection stays in
/// this set and signs can be evaluated with integer arithmetic only.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct DyadicRational {
    numerator: BigInt,
    exponent: u32,
}

impl DyadicRational {
    /// Construct `numerator / 2^exponent`, reduced to lowest terms.
    pub fn new(numerator: BigInt, exponent: u32) -> Self {
        let mut numerator = numerator;
        let mut exponent = exponent;
        if numerator.is_zero() {
            exponent = 0;
        }
        while exponent > 0 && numerator.is_even() {
            numerator /= 2;
            exponent -= 1;
        }
        DyadicRational {
            numerator,
            exponent,
        }
    }

    /// The integer `value` as a dyadic rational.
    pub fn from_integer(value: BigInt) -> Self {
        DyadicRational {
            numerator: value,
            exponent: 0,
        }
    }

    /// The reduced numerator.
    pub fn numerator(&self) -> &BigInt {
        &self.numerator
    }

    /// The reduced power of two in the denominator.
    pub fn exponent(&self) -> u32 {
        self.exponent
    }

    /// The denominator `2^exponent`.
    pub fn denominator(&self) -> BigInt {
        BigInt::one() << self.exponent
    }

    /// Numerator after rescaling to the denominator `2^exponent`, which must
    /// be at least `self.exponent()`.
    fn numerator_at_exponent(&self, exponent: u32) -> BigInt {
        debug_assert!(exponent >= self.exponent);
        &self.numerator << (exponent - self.exponent)
    }

    /// The exact midpoint `(self + other) / 2`.
    pub fn midpoint(&self, other: &Self) -> Self {
        let exponent = self.exponent.max(other.exponent);
        let sum = self.numerator_at_exponent(exponent) + other.numerator_at_exponent(exponent);
        DyadicRational::new(sum, exponent + 1)
    }

    /// `self + numerator / 2^exponent`.
    fn add_dyadic(&self, numerator: &BigInt, exponent: u32) -> Self {
        let common = self.exponent.max(exponent);
        let other = numerator << (common - exponent);
        DyadicRational::new(self.numerator_at_exponent(common) + other, common)
    }

    /// Sign of `polynomial(self)` for ascending integer coefficients.
    ///
    /// The value is multiplied by the positive factor `2^(exponent * degree)`,
    /// so only integer arithmetic is needed.
    fn sign_of(&self, polynomial: &[BigInt]) -> i8 {
        let Some(d) = degree(polynomial) else {
            return 0;
        };
        let denominator = self.denominator();
        let mut accumulator = polynomial[d].clone();
        let mut denominator_power = BigInt::one();
        for coefficient in polynomial[..d].iter().rev() {
            denominator_power *= &denominator;
            accumulator = accumulator * &self.numerator + coefficient * &denominator_power;
        }
        sign_i8(&accumulator)
    }
}

impl PartialOrd for DyadicRational {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for DyadicRational {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let exponent = self.exponent.max(other.exponent);
        self.numerator_at_exponent(exponent)
            .cmp(&other.numerator_at_exponent(exponent))
    }
}

impl std::fmt::Display for DyadicRational {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.exponent == 0 {
            write!(f, "{}", self.numerator)
        } else {
            write!(f, "{}/2^{}", self.numerator, self.exponent)
        }
    }
}

/// One isolated distinct real root of an integer polynomial.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IsolatedRealRoot {
    /// The root is exactly this dyadic rational.
    Exact(DyadicRational),
    /// Exactly one distinct root lies in the open interval `(lower, upper)`;
    /// neither endpoint is a root.
    Interval {
        /// Lower endpoint (not a root).
        lower: DyadicRational,
        /// Upper endpoint (not a root).
        upper: DyadicRational,
    },
}

impl IsolatedRealRoot {
    /// Smallest point of the closed isolating set.
    pub fn lower(&self) -> &DyadicRational {
        match self {
            IsolatedRealRoot::Exact(root) => root,
            IsolatedRealRoot::Interval { lower, .. } => lower,
        }
    }

    /// Largest point of the closed isolating set.
    pub fn upper(&self) -> &DyadicRational {
        match self {
            IsolatedRealRoot::Exact(root) => root,
            IsolatedRealRoot::Interval { upper, .. } => upper,
        }
    }
}

fn sturm_variations_at(chain: &[Vec<BigInt>], point: &DyadicRational) -> usize {
    sign_variations(chain.iter().map(|p| point.sign_of(p)))
}

enum IsolationWork {
    Interval {
        lower: DyadicRational,
        upper: DyadicRational,
        lower_variations: usize,
        upper_variations: usize,
    },
    Emit(IsolatedRealRoot),
}

/// Isolate the distinct real roots of an integer polynomial.
///
/// Coefficients are in ascending degree order. The result has one entry per
/// distinct real root, sorted increasingly. Each entry is either an exact
/// dyadic root or an open interval with dyadic endpoints containing exactly
/// one root; interval endpoints are never roots. Consecutive entries satisfy
/// `previous.upper() <= next.lower()`, with strict inequality whenever one
/// of them is [`IsolatedRealRoot::Exact`], so neighbouring intervals may share
/// a (non-root) endpoint but never overlap.
///
/// The method is Sturm bisection on the primitive square-free part, starting
/// from a strict Fujiwara bound, with exact integer sign evaluation. Constant
/// polynomials, including the zero polynomial, return an empty vector: the
/// caller must treat the zero polynomial separately.
pub fn isolate_real_roots_bigint_coeffs(coeffs: &[BigInt]) -> Vec<IsolatedRealRoot> {
    let p = squarefree_part_primitive_bigint_coeffs(coeffs);
    if degree(&p).unwrap_or(0) == 0 {
        return Vec::new();
    }
    let bound = strict_fujiwara_root_bound(&p);
    let chain = sturm_chain_from_squarefree(p);

    let lower = DyadicRational::from_integer(-bound.clone());
    let upper = DyadicRational::from_integer(bound);
    let mut stack = vec![IsolationWork::Interval {
        lower_variations: sturm_variations_at(&chain, &lower),
        upper_variations: sturm_variations_at(&chain, &upper),
        lower,
        upper,
    }];
    let mut roots = Vec::new();

    while let Some(work) = stack.pop() {
        let (lower, upper, lower_variations, upper_variations) = match work {
            IsolationWork::Emit(root) => {
                roots.push(root);
                continue;
            }
            IsolationWork::Interval {
                lower,
                upper,
                lower_variations,
                upper_variations,
            } => (lower, upper, lower_variations, upper_variations),
        };
        let count = lower_variations.saturating_sub(upper_variations);
        if count == 0 {
            continue;
        }
        if count == 1 {
            roots.push(IsolatedRealRoot::Interval { lower, upper });
            continue;
        }

        let middle = lower.midpoint(&upper);
        if middle.sign_of(&chain[0]) != 0 {
            let middle_variations = sturm_variations_at(&chain, &middle);
            stack.push(IsolationWork::Interval {
                lower: middle.clone(),
                upper,
                lower_variations: middle_variations,
                upper_variations,
            });
            stack.push(IsolationWork::Interval {
                lower,
                upper: middle,
                lower_variations,
                upper_variations: middle_variations,
            });
            continue;
        }

        // The midpoint is a root: shrink a symmetric window around it until it
        // isolates that root, then continue on both sides of the window.
        let width_exponent = lower.exponent().max(upper.exponent());
        let width = upper.numerator_at_exponent(width_exponent)
            - lower.numerator_at_exponent(width_exponent);
        let mut half_window_exponent = width_exponent + 2;
        let (left, right, left_variations, right_variations) = loop {
            let left = middle.add_dyadic(&-&width, half_window_exponent);
            let right = middle.add_dyadic(&width, half_window_exponent);
            if left.sign_of(&chain[0]) != 0 && right.sign_of(&chain[0]) != 0 {
                let left_variations = sturm_variations_at(&chain, &left);
                let right_variations = sturm_variations_at(&chain, &right);
                if left_variations.saturating_sub(right_variations) == 1 {
                    break (left, right, left_variations, right_variations);
                }
            }
            half_window_exponent += 1;
        };
        stack.push(IsolationWork::Interval {
            lower: right,
            upper,
            lower_variations: right_variations,
            upper_variations,
        });
        stack.push(IsolationWork::Emit(IsolatedRealRoot::Exact(middle)));
        stack.push(IsolationWork::Interval {
            lower,
            upper: left,
            lower_variations,
            upper_variations: left_variations,
        });
    }
    roots
}
