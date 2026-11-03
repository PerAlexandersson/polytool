//! Dense univariate polynomials over the rationals as ascending coefficient
//! vectors (entry `i` is the coefficient of `t^i`).  The arithmetic returns
//! the canonical form `[0]` for the zero polynomial; [`strip_trailing_zeros`]
//! gives the empty vector instead.
//!
//! Most coefficients met in practice are integers.  The rational type
//! reduces every sum and product by a gcd with the denominator, which for
//! a denominator of one costs a binary gcd over every bit of the numerator;
//! [`rational_product`], [`add_rational_into`] and
//! [`subtract_rational_from`] skip that reduction when both operands are
//! integers (an integer over one is already reduced, so the result is the
//! same rational).

use crate::Rational;
use num_bigint::BigInt;
use num_traits::{One, Zero};

/// `a b`, without the gcd reduction when both are integers (the same
/// value as `a * b`).
#[inline]
pub(crate) fn rational_product(a: &Rational, b: &Rational) -> Rational {
    if a.denom().is_one() && b.denom().is_one() {
        Rational::new_raw(a.numer() * b.numer(), BigInt::one())
    } else {
        a * b
    }
}

/// `target += value`, without the gcd reduction when both are integers.
#[inline]
pub(crate) fn add_rational_into(target: &mut Rational, value: &Rational) {
    if target.denom().is_one() && value.denom().is_one() {
        let (mut numerator, denominator) = std::mem::take(target).into_raw();
        numerator += value.numer();
        *target = Rational::new_raw(numerator, denominator);
    } else {
        *target += value;
    }
}

/// `target -= value`, without the gcd reduction when both are integers.
#[inline]
pub(crate) fn subtract_rational_from(target: &mut Rational, value: &Rational) {
    if target.denom().is_one() && value.denom().is_one() {
        let (mut numerator, denominator) = std::mem::take(target).into_raw();
        numerator -= value.numer();
        *target = Rational::new_raw(numerator, denominator);
    } else {
        *target -= value;
    }
}

/// `target -= a b`, without gcd reductions when all three are integers.
#[inline]
pub(crate) fn subtract_product_from(target: &mut Rational, a: &Rational, b: &Rational) {
    if target.denom().is_one() && a.denom().is_one() && b.denom().is_one() {
        let (mut numerator, denominator) = std::mem::take(target).into_raw();
        numerator -= a.numer() * b.numer();
        *target = Rational::new_raw(numerator, denominator);
    } else {
        *target -= a * b;
    }
}

/// `target += a b`, without gcd reductions when all three are integers.
#[inline]
pub(crate) fn add_product_into(target: &mut Rational, a: &Rational, b: &Rational) {
    if target.denom().is_one() && a.denom().is_one() && b.denom().is_one() {
        let (mut numerator, denominator) = std::mem::take(target).into_raw();
        numerator += a.numer() * b.numer();
        *target = Rational::new_raw(numerator, denominator);
    } else {
        *target += a * b;
    }
}

/// Degree, or `None` for the zero polynomial.
pub(crate) fn degree(p: &[Rational]) -> Option<usize> {
    p.iter().rposition(|c| !c.is_zero())
}

/// Whether every coefficient is zero.
pub(crate) fn is_zero(p: &[Rational]) -> bool {
    p.iter().all(Zero::is_zero)
}

/// The polynomial without trailing zeros; zero becomes `[0]`.
pub(crate) fn trim(p: &[Rational]) -> Vec<Rational> {
    match degree(p) {
        Some(d) => p[..=d].to_vec(),
        None => vec![Rational::zero()],
    }
}

/// Remove trailing zeros in place; zero becomes the empty vector.
pub(crate) fn strip_trailing_zeros(p: &mut Vec<Rational>) {
    let len = degree(p).map_or(0, |d| d + 1);
    p.truncate(len);
}

fn canonical(mut p: Vec<Rational>) -> Vec<Rational> {
    strip_trailing_zeros(&mut p);
    if p.is_empty() {
        p.push(Rational::zero());
    }
    p
}

pub(crate) fn add(a: &[Rational], b: &[Rational]) -> Vec<Rational> {
    let (long, short) = if a.len() >= b.len() { (a, b) } else { (b, a) };
    let mut r = long.to_vec();
    for (x, y) in r.iter_mut().zip(short) {
        add_rational_into(x, y);
    }
    canonical(r)
}

pub(crate) fn mul(a: &[Rational], b: &[Rational]) -> Vec<Rational> {
    let (Some(da), Some(db)) = (degree(a), degree(b)) else {
        return vec![Rational::zero()];
    };
    let mut r = vec![Rational::zero(); da + db + 1];
    for (i, x) in a[..=da].iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        for (j, y) in b[..=db].iter().enumerate() {
            if !y.is_zero() {
                add_product_into(&mut r[i + j], x, y);
            }
        }
    }
    canonical(r)
}

pub(crate) fn scale(p: &[Rational], c: &Rational) -> Vec<Rational> {
    if c.is_zero() {
        return vec![Rational::zero()];
    }
    canonical(p.iter().map(|x| rational_product(x, c)).collect())
}

/// Multiply by `t^k`.
pub(crate) fn shift_by(p: &[Rational], k: usize) -> Vec<Rational> {
    let Some(d) = degree(p) else {
        return vec![Rational::zero()];
    };
    let mut r = vec![Rational::zero(); k];
    r.extend_from_slice(&p[..=d]);
    r
}

pub(crate) fn evaluate(p: &[Rational], x: &Rational) -> Rational {
    let Some(d) = degree(p) else {
        return Rational::zero();
    };
    let mut result = p[d].clone();
    for c in p[..d].iter().rev() {
        result = rational_product(&result, x);
        add_rational_into(&mut result, c);
    }
    result
}

/// `p(q(t))`.
pub(crate) fn compose(p: &[Rational], q: &[Rational]) -> Vec<Rational> {
    let Some(d) = degree(p) else {
        return vec![Rational::zero()];
    };
    let mut result = vec![p[d].clone()];
    for c in p[..d].iter().rev() {
        result = mul(&result, q);
        add_rational_into(&mut result[0], c);
        result = canonical(result);
    }
    result
}

/// The translate `p(t + v)` (Horner's scheme with the linear factor
/// `t + v`).
pub(crate) fn translate(p: &[Rational], v: &Rational) -> Vec<Rational> {
    let Some(d) = degree(p) else {
        return vec![Rational::zero()];
    };
    if v.is_zero() {
        return p[..=d].to_vec();
    }
    let mut result = vec![Rational::zero(); d + 1];
    result[0] = p[d].clone();
    for (len, c) in (1..).zip(p[..d].iter().rev()) {
        // result <- result (t + v) + c, with result of degree len - 1.
        for j in (0..=len).rev() {
            let mut value = if j < len && !result[j].is_zero() {
                rational_product(&result[j], v)
            } else {
                Rational::zero()
            };
            if j > 0 {
                add_rational_into(&mut value, &result[j - 1]);
            }
            result[j] = value;
        }
        add_rational_into(&mut result[0], c);
    }
    canonical(result)
}

pub(crate) fn derivative(p: &[Rational]) -> Vec<Rational> {
    let Some(d) = degree(p) else {
        return vec![Rational::zero()];
    };
    if d == 0 {
        return vec![Rational::zero()];
    }
    canonical(
        (1..=d)
            .map(|k| rational_product(&p[k], &Rational::from_integer(k.into())))
            .collect(),
    )
}

/// The remainder of `a` modulo a nonzero `b`.
fn remainder(a: &[Rational], b: &[Rational]) -> Vec<Rational> {
    let db = degree(b).expect("nonzero divisor");
    let lead = b[db].clone();
    let mut r = trim(a);
    while let Some(dr) = degree(&r) {
        if dr < db {
            break;
        }
        let factor = &r[dr] / &lead;
        for (j, c) in b[..=db].iter().enumerate() {
            subtract_product_from(&mut r[dr - db + j], &factor, c);
        }
        r = trim(&r);
    }
    r
}

/// The monic greatest common divisor (`[0]` when both are zero).
pub(crate) fn gcd(a: &[Rational], b: &[Rational]) -> Vec<Rational> {
    let (mut a, mut b) = (trim(a), trim(b));
    while !is_zero(&b) {
        let r = remainder(&a, &b);
        a = b;
        b = r;
    }
    match degree(&a) {
        Some(d) => {
            let lead = a[d].clone();
            scale(&a, &lead.recip())
        }
        None => vec![Rational::zero()],
    }
}

/// The quotient of an exact division by a nonzero `b`.
pub(crate) fn exact_quotient(a: &[Rational], b: &[Rational]) -> Vec<Rational> {
    let db = degree(b).expect("nonzero divisor");
    let lead = b[db].clone();
    let mut r = trim(a);
    let Some(da) = degree(&r) else {
        return vec![Rational::zero()];
    };
    if da < db {
        return vec![Rational::zero()];
    }
    let mut q = vec![Rational::zero(); da - db + 1];
    while let Some(dr) = degree(&r) {
        if dr < db {
            break;
        }
        let factor = &r[dr] / &lead;
        for (j, c) in b[..=db].iter().enumerate() {
            subtract_product_from(&mut r[dr - db + j], &factor, c);
        }
        q[dr - db] = factor;
        r = trim(&r);
    }
    trim(&q)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(numerator: i64, denominator: i64) -> Rational {
        Rational::new(numerator.into(), denominator.into())
    }

    /// The integer fast paths give the same rationals as the generic
    /// operators, for integer and non-integer operands of either sign.
    #[test]
    fn fast_paths_agree_with_the_rational_operators() {
        let big = Rational::from_integer(BigInt::from(10).pow(40u32) + BigInt::from(3));
        let values = [
            q(0, 1),
            q(1, 1),
            q(-7, 1),
            q(3, 4),
            q(-5, 6),
            big.clone(),
            -big,
        ];
        for a in &values {
            for b in &values {
                assert_eq!(rational_product(a, b), a * b);
                let mut sum = a.clone();
                add_rational_into(&mut sum, b);
                assert_eq!(sum, a + b);
                let mut difference = a.clone();
                subtract_rational_from(&mut difference, b);
                assert_eq!(difference, a - b);
                for c in &values {
                    let mut target = c.clone();
                    add_product_into(&mut target, a, b);
                    assert_eq!(target, c + a * b);
                    let mut target = c.clone();
                    subtract_product_from(&mut target, a, b);
                    assert_eq!(target, c - a * b);
                }
            }
        }
    }

    #[test]
    fn translate_agrees_with_composition() {
        let p = [q(1, 1), q(-2, 3), q(0, 1), q(5, 1), q(7, 2)];
        for v in [q(0, 1), q(1, 1), q(-3, 1), q(2, 5)] {
            assert_eq!(translate(&p, &v), compose(&p, &[v.clone(), q(1, 1)]));
        }
        assert_eq!(translate(&[q(0, 1)], &q(2, 1)), vec![q(0, 1)]);
        assert_eq!(translate(&[q(4, 1)], &q(2, 1)), vec![q(4, 1)]);
    }
}
