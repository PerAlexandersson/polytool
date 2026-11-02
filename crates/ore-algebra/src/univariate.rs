//! Dense univariate polynomials over the rationals as ascending coefficient
//! vectors (entry `i` is the coefficient of `t^i`).  The arithmetic returns
//! the canonical form `[0]` for the zero polynomial; [`strip_trailing_zeros`]
//! gives the empty vector instead.

use crate::Rational;
use num_traits::Zero;

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
        *x += y;
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
                r[i + j] += x * y;
            }
        }
    }
    canonical(r)
}

pub(crate) fn scale(p: &[Rational], c: &Rational) -> Vec<Rational> {
    if c.is_zero() {
        return vec![Rational::zero()];
    }
    canonical(p.iter().map(|x| x * c).collect())
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
        result = result * x + c;
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
        result[0] += c;
        result = canonical(result);
    }
    result
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
            .map(|k| &p[k] * Rational::from_integer(k.into()))
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
            r[dr - db + j] -= &factor * c;
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
            r[dr - db + j] -= &factor * c;
        }
        q[dr - db] = factor;
        r = trim(&r);
    }
    trim(&q)
}
