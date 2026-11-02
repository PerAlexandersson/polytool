//! The Euler form `sum x^s a_(i,s)(n, theta) S^i` of an operator, with the
//! Euler operator `theta = x D`, and the conversions to and from the
//! normal-ordered form.

use crate::bivariate::BivariatePolynomial;
use crate::operator::{rational, OreOperator};
use crate::univariate as poly;
use crate::Rational;
use num_traits::{One, Zero};
use std::collections::BTreeMap;

/// `sum x^s a_(i,s)(n, theta) S^i`, keyed by `(shift i, x power s)`, with
/// coefficients in `(n, theta)` (the variable `y` of the
/// [`BivariatePolynomial`] is `theta`).
pub type EulerOperator = BTreeMap<(usize, usize), BivariatePolynomial>;

/// Add `c` to the coefficient of `key` in an Euler-form operator, removing
/// zero coefficients.
pub fn add_euler_term(op: &mut EulerOperator, key: (usize, usize), c: BivariatePolynomial) {
    if c.is_zero() {
        return;
    }
    let entry = op.entry(key).or_default();
    *entry = entry.add(&c);
    if entry.is_zero() {
        op.remove(&key);
    }
}

/// The falling factorial `theta (theta - 1) ... (theta - d + 1)`
/// (ascending coefficients in `theta`).
pub fn falling_factorial(d: usize) -> Vec<Rational> {
    (0..d).fold(vec![Rational::one()], |acc, t| {
        poly::mul(&acc, &[-rational(t), Rational::one()])
    })
}

/// Stirling numbers of the second kind `S(m, d)` for `m <= max`.
fn stirling_second_kind(max: usize) -> Vec<Vec<Rational>> {
    let mut table = vec![vec![Rational::zero(); max + 1]; max + 1];
    table[0][0] = Rational::one();
    for m in 1..=max {
        for d in 1..=m {
            table[m][d] = &table[m - 1][d - 1] + rational(d) * &table[m - 1][d];
        }
    }
    table
}

/// `(J, x^J L)` in Euler form with `J` the largest derivative order of `L`,
/// using `x^d D^d = theta (theta - 1) ... (theta - d + 1)`.
pub fn to_euler_form(op: &OreOperator) -> (usize, EulerOperator) {
    let j_max = op.max_derivative();
    let mut out = EulerOperator::new();
    for ((i, d), c) in &op.0 {
        let falling = falling_factorial(*d);
        for (r, row) in c.rows().iter().enumerate() {
            for (j, cj) in row.iter().enumerate() {
                if cj.is_zero() {
                    continue;
                }
                let s = j + j_max - d;
                let mut rows = vec![Vec::new(); r + 1];
                rows[r] = poly::scale(&falling, cj);
                add_euler_term(&mut out, (*i, s), BivariatePolynomial::new(rows));
            }
        }
    }
    (j_max, out)
}

/// The normal-ordered form of an Euler-form operator, using
/// `theta^m = sum_d S(m, d) x^d D^d`.
pub fn from_euler_form(op: &EulerOperator) -> OreOperator {
    let max_m = op
        .values()
        .filter_map(BivariatePolynomial::variable_degree)
        .max()
        .unwrap_or(0);
    let table = stirling_second_kind(max_m);
    let mut out = OreOperator::default();
    for ((i, s), a) in op {
        for (r, row) in a.rows().iter().enumerate() {
            for (m, alpha) in row.iter().enumerate() {
                if alpha.is_zero() {
                    continue;
                }
                for (d, stirling) in table[m].iter().enumerate().take(m + 1) {
                    if stirling.is_zero() {
                        continue;
                    }
                    out.add_term(
                        (*i, d),
                        BivariatePolynomial::monomial(r, s + d, alpha * stirling),
                    );
                }
            }
        }
    }
    out
}
