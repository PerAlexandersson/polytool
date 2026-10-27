//! Sturm chains for real root isolation of polynomials with rational coefficients.
//!
//! [`SturmChain`] is built from the square-free part of an integer polynomial
//! (`i64` or `BigInt` coefficients in ascending order) with exact rational
//! arithmetic, and supports exact root counting on intervals
//! ([`SturmChain::count_roots_in`], [`SturmChain::count_real_roots`]) and
//! isolation of the real roots by rational intervals
//! ([`SturmChain::isolate_roots`]).  Repeated roots are counted once.  This is
//! the slower route for real-rootedness; the default exact tests live in
//! [`crate::real_rootedness`] and [`crate::root_count`].

use crate::BigRational;
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

type Q = BigRational;

/// A polynomial over Q stored as coefficient vector (ascending degree).
#[derive(Clone, Debug)]
struct QPoly {
    coeffs: Vec<Q>,
}

impl QPoly {
    fn new(mut coeffs: Vec<Q>) -> Self {
        while coeffs.last().is_some_and(|c| c.is_zero()) {
            coeffs.pop();
        }
        QPoly { coeffs }
    }

    fn zero() -> Self {
        QPoly { coeffs: vec![] }
    }

    fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    fn degree(&self) -> Option<usize> {
        if self.coeffs.is_empty() {
            None
        } else {
            Some(self.coeffs.len() - 1)
        }
    }

    fn leading(&self) -> Q {
        self.coeffs.last().cloned().unwrap_or_else(Q::zero)
    }

    fn evaluate(&self, x: &Q) -> Q {
        if self.is_zero() {
            return Q::zero();
        }
        let mut result = self.coeffs.last().unwrap().clone();
        for c in self.coeffs.iter().rev().skip(1) {
            result = result * x + c;
        }
        result
    }

    fn derivative(&self) -> Self {
        if self.coeffs.len() <= 1 {
            return Self::zero();
        }
        let coeffs: Vec<Q> = self.coeffs[1..]
            .iter()
            .enumerate()
            .map(|(i, c)| c * Q::from_integer(BigInt::from(i + 1)))
            .collect();
        Self::new(coeffs)
    }

    fn scale(&self, c: &Q) -> Self {
        Self::new(self.coeffs.iter().map(|a| a * c).collect())
    }

    fn neg(&self) -> Self {
        Self::new(self.coeffs.iter().map(|a| -a).collect())
    }

    /// Polynomial remainder: self mod other.
    fn rem(&self, other: &QPoly) -> QPoly {
        if other.is_zero() {
            panic!("division by zero polynomial");
        }
        let mut rem = self.coeffs.clone();
        let d_other = other.degree().unwrap();
        let lc_other = other.leading();

        while rem.len() > d_other {
            let lc_rem = rem.last().unwrap().clone();
            if lc_rem.is_zero() {
                rem.pop();
                continue;
            }
            let factor = &lc_rem / &lc_other;
            let shift = rem.len() - 1 - d_other;
            for (i, c) in other.coeffs.iter().enumerate() {
                rem[shift + i] = &rem[shift + i] - &factor * c;
            }
            rem.pop(); // remove leading zero
        }
        // Strip trailing zeros
        while rem.last().is_some_and(|c| c.is_zero()) {
            rem.pop();
        }
        QPoly::new(rem)
    }

    /// Make the polynomial square-free: p / gcd(p, p').
    fn square_free(&self) -> QPoly {
        if self.is_zero() {
            return self.clone();
        }
        let dp = self.derivative();
        if dp.is_zero() {
            return self.clone();
        }
        let g = poly_gcd(self, &dp);
        if g.degree() == Some(0) || g.is_zero() {
            return self.clone();
        }
        poly_exact_div(self, &g)
    }
}

/// GCD of two polynomials over Q via Euclidean algorithm.
fn poly_gcd(a: &QPoly, b: &QPoly) -> QPoly {
    let mut r0 = a.clone();
    let mut r1 = b.clone();
    while !r1.is_zero() {
        let rem = r0.rem(&r1);
        r0 = r1;
        r1 = rem;
    }
    // Make monic
    if !r0.is_zero() {
        let lc = r0.leading();
        r0 = r0.scale(&(Q::one() / lc));
    }
    r0
}

/// Exact polynomial division (assumes b divides a).
fn poly_exact_div(a: &QPoly, b: &QPoly) -> QPoly {
    if b.is_zero() {
        panic!("division by zero");
    }
    if a.is_zero() {
        return QPoly::zero();
    }
    let da = a.degree().unwrap();
    let db = b.degree().unwrap();
    if da < db {
        return QPoly::zero();
    }
    let mut rem = a.coeffs.clone();
    let lc_b = b.leading();
    let dq = da - db;
    let mut quot = vec![Q::zero(); dq + 1];

    for i in (0..=dq).rev() {
        let lc_rem = rem[i + db].clone();
        if lc_rem.is_zero() {
            continue;
        }
        let q = &lc_rem / &lc_b;
        quot[i] = q.clone();
        for (j, c) in b.coeffs.iter().enumerate() {
            rem[i + j] = &rem[i + j] - &q * c;
        }
    }
    QPoly::new(quot)
}

// ---------------------------------------------------------------------------
// Sturm chain
// ---------------------------------------------------------------------------

/// A Sturm chain for counting and isolating real roots of a square-free polynomial.
///
/// Built from the square-free part of the input, so multiple roots of the input
/// are counted once.  All computations are exact over `Q`.
pub struct SturmChain {
    chain: Vec<QPoly>,
}

impl SturmChain {
    /// Build a Sturm chain from a polynomial (given as i64 coefficients, ascending order).
    ///
    /// The polynomial is automatically made square-free.
    pub fn from_i64_coeffs(coeffs: &[i64]) -> Self {
        let qcoeffs: Vec<Q> = coeffs
            .iter()
            .map(|&c| Q::from_integer(BigInt::from(c)))
            .collect();
        let p = QPoly::new(qcoeffs).square_free();
        Self::build(p)
    }

    /// Build a Sturm chain from a polynomial with `BigInt` coefficients
    /// in ascending order.
    ///
    /// The polynomial is automatically made square-free.
    pub fn from_bigint_coeffs(coeffs: &[BigInt]) -> Self {
        let qcoeffs: Vec<Q> = coeffs.iter().map(|c| Q::from_integer(c.clone())).collect();
        let p = QPoly::new(qcoeffs).square_free();
        Self::build(p)
    }

    fn build(p: QPoly) -> Self {
        if p.is_zero() || p.degree() <= Some(0) {
            return SturmChain { chain: vec![p] };
        }
        let mut chain = Vec::new();
        let dp = p.derivative();
        chain.push(p);
        chain.push(dp);

        loop {
            let len = chain.len();
            let rem = chain[len - 2].rem(&chain[len - 1]);
            if rem.is_zero() {
                break;
            }
            chain.push(rem.neg());
        }
        SturmChain { chain }
    }

    /// Count sign changes in the Sturm chain evaluated at x.
    fn sign_changes_at(&self, x: &Q) -> usize {
        let signs: Vec<i8> = self
            .chain
            .iter()
            .map(|p| {
                let v = p.evaluate(x);
                if v.is_positive() {
                    1
                } else if v.is_negative() {
                    -1
                } else {
                    0
                }
            })
            .collect();

        let mut changes = 0;
        let mut prev = 0i8;
        for &s in &signs {
            if s == 0 {
                continue;
            }
            if prev != 0 && prev != s {
                changes += 1;
            }
            prev = s;
        }
        changes
    }

    /// Count distinct real roots in the open interval (a, b).
    ///
    /// Roots at the endpoints are not counted.  Expects `a <= b`; the result is not
    /// meaningful otherwise.
    pub fn count_roots_in(&self, a: &BigRational, b: &BigRational) -> usize {
        let sa = self.sign_changes_at(a);
        let sb = self.sign_changes_at(b);
        let mut count = sa.saturating_sub(sb);
        // Sturm's theorem counts roots in (a, b]. Subtract root at b if present.
        if count > 0 && self.chain[0].evaluate(b).is_zero() {
            count -= 1;
        }
        count
    }

    /// Total number of distinct real roots (the zero polynomial has none).
    pub fn count_real_roots(&self) -> usize {
        let bound = self.root_bound();
        let neg = -bound.clone();
        self.count_roots_in(&neg, &bound)
    }

    /// A Cauchy bound on the absolute value of any root.
    fn root_bound(&self) -> Q {
        let p = &self.chain[0];
        if p.is_zero() || p.degree() == Some(0) {
            return Q::one();
        }
        let lc = p.leading();
        let lc_abs = if lc.is_negative() {
            -lc.clone()
        } else {
            lc.clone()
        };
        let mut bound = Q::one();
        for c in &p.coeffs {
            let c_abs = if c.is_negative() {
                -c.clone()
            } else {
                c.clone()
            };
            let ratio = &c_abs / &lc_abs + Q::one();
            if ratio > bound {
                bound = ratio;
            }
        }
        bound
    }

    /// Degree of the square-free polynomial in this chain.
    pub fn square_free_degree(&self) -> usize {
        self.chain[0].degree().unwrap_or(0)
    }

    /// Find all real roots as isolating intervals, refined to width ≤ epsilon.
    /// Returns a vector of (lo, hi) pairs in increasing order.
    ///
    /// Each pair contains exactly one distinct root: in the open interval
    /// `(lo, hi)`, or `lo == hi` when the root is found exactly (a bisection
    /// midpoint).
    ///
    /// # Panics
    ///
    /// Panics if `epsilon` is not positive.
    pub fn isolate_roots(&self, epsilon: &BigRational) -> Vec<(BigRational, BigRational)> {
        assert!(
            epsilon.is_positive(),
            "root isolation epsilon must be positive"
        );
        let bound = self.root_bound();
        let neg_bound = -bound.clone();
        self.isolate_in(&neg_bound, &bound, epsilon)
    }

    fn isolate_in(&self, a: &Q, b: &Q, epsilon: &Q) -> Vec<(Q, Q)> {
        let count = self.count_roots_in(a, b); // open interval (a, b)
        if count == 0 {
            return vec![];
        }
        if count == 1 {
            // Refine this interval by bisection
            let mut lo = a.clone();
            let mut hi = b.clone();
            while &hi - &lo > *epsilon {
                let mid = (&lo + &hi) / Q::from_integer(BigInt::from(2));
                // If mid is an exact root, return it
                if self.chain[0].evaluate(&mid).is_zero() {
                    return vec![(mid.clone(), mid.clone())];
                }
                let left = self.count_roots_in(&lo, &mid);
                if left == 1 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            return vec![(lo, hi)];
        }
        // Multiple roots: bisect and recurse.
        let mid = (a + b) / Q::from_integer(BigInt::from(2));
        let mut roots = self.isolate_in(a, &mid, epsilon);
        // Check if mid itself is a root (it's excluded from both open sub-intervals).
        if self.chain[0].evaluate(&mid).is_zero() {
            roots.push((mid.clone(), mid.clone()));
        }
        roots.extend(self.isolate_in(&mid, b, epsilon));
        roots
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_count_roots_real_rooted() {
        // (1 + t)^3 = 1 + 3t + 3t^2 + t^3: square-free = 1+t, 1 real root
        let sc = SturmChain::from_i64_coeffs(&[1, 3, 3, 1]);
        assert_eq!(sc.count_real_roots(), 1);
    }

    #[test]
    fn test_count_roots_distinct() {
        // (t-1)(t-2)(t-3) = -6 + 11t - 6t^2 + t^3
        let sc = SturmChain::from_i64_coeffs(&[-6, 11, -6, 1]);
        assert_eq!(sc.count_real_roots(), 3);
    }

    #[test]
    fn test_count_roots_complex() {
        // t^2 + 1: no real roots
        let sc = SturmChain::from_i64_coeffs(&[1, 0, 1]);
        assert_eq!(sc.count_real_roots(), 0);
    }

    #[test]
    fn test_isolate_roots() {
        // (t-1)(t+2) = -2 + t + t^2... no: (t-1)(t+2) = t^2 + t - 2
        let sc = SturmChain::from_i64_coeffs(&[-2, 1, 1]);
        let eps = Q::new(BigInt::from(1), BigInt::from(1000));
        let roots = sc.isolate_roots(&eps);
        assert_eq!(roots.len(), 2);
        // Roots should be near -2 and 1
        let (lo0, _hi0) = &roots[0];
        assert!(lo0 < &Q::from_integer(BigInt::from(-1)));
        let (_lo1, hi1) = &roots[1];
        assert!(hi1 > &Q::from_integer(BigInt::from(0)));
    }

    #[test]
    #[should_panic(expected = "root isolation epsilon must be positive")]
    fn test_isolate_roots_rejects_zero_epsilon() {
        let sc = SturmChain::from_i64_coeffs(&[-2, 1, 1]);
        sc.isolate_roots(&Q::zero());
    }

    #[test]
    #[should_panic(expected = "root isolation epsilon must be positive")]
    fn test_isolate_roots_rejects_negative_epsilon() {
        let sc = SturmChain::from_i64_coeffs(&[-2, 1, 1]);
        sc.isolate_roots(&Q::from_integer(BigInt::from(-1)));
    }

    #[test]
    fn test_square_free() {
        // (1+t)^2 = 1 + 2t + t^2, square-free = 1+t
        let p = QPoly::new(vec![
            Q::from_integer(BigInt::from(1)),
            Q::from_integer(BigInt::from(2)),
            Q::from_integer(BigInt::from(1)),
        ]);
        let sf = p.square_free();
        assert_eq!(sf.degree(), Some(1));
    }
}
