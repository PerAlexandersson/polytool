//! Private rational-polynomial helpers (gcd, remainder, exact division and
//! conversion to primitive integer coefficients) over `Ratio<BigInt>`.
//!
//! Used by the weak-interlacing check to divide out common factors.

use super::Q;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{Signed, Zero};

pub(super) fn q_degree(p: &[Q]) -> usize {
    let zero = Q::from_integer(BigInt::from(0));
    p.iter().rposition(|c| *c != zero).unwrap_or(0)
}

pub(super) fn poly_gcd_q(a: &[Q], b: &[Q]) -> Vec<Q> {
    let zero = Q::from_integer(BigInt::from(0));
    let mut r0: Vec<Q> = a.to_vec();
    let mut r1: Vec<Q> = b.to_vec();

    // Trim trailing zeros
    while r0.last().is_some_and(|c| *c == zero) {
        r0.pop();
    }
    while r1.last().is_some_and(|c| *c == zero) {
        r1.pop();
    }

    while !r1.is_empty() && r1.iter().any(|c| *c != zero) {
        let rem = poly_rem_q(&r0, &r1);
        r0 = r1;
        r1 = rem;
    }
    // Make monic
    if !r0.is_empty() {
        let lc = r0.last().unwrap().clone();
        if lc != zero {
            let inv = Q::from_integer(BigInt::from(1)) / lc;
            for c in r0.iter_mut() {
                *c = c.clone() * inv.clone();
            }
        }
    }
    // Trim
    while r0.last().is_some_and(|c| *c == zero) {
        r0.pop();
    }
    if r0.is_empty() {
        r0.push(zero);
    }
    r0
}

fn poly_rem_q(a: &[Q], b: &[Q]) -> Vec<Q> {
    let zero = Q::from_integer(BigInt::from(0));
    if b.is_empty() || b.iter().all(|c| *c == zero) {
        panic!("division by zero polynomial");
    }
    let mut rem = a.to_vec();
    let db = q_degree(b);
    let lc_b = b[db].clone();

    while rem.len() > db + 1 || (rem.len() == db + 1 && q_degree(&rem) >= db) {
        let dr = q_degree(&rem);
        if dr < db {
            break;
        }
        let lc_r = rem[dr].clone();
        if lc_r == zero {
            rem.pop();
            continue;
        }
        let factor = lc_r / lc_b.clone();
        let shift = dr - db;
        for (i, c) in b.iter().enumerate() {
            rem[shift + i] = rem[shift + i].clone() - factor.clone() * c.clone();
        }
        // The leading term should now be zero; pop it
        while rem.last().is_some_and(|c| *c == zero) {
            rem.pop();
        }
    }
    if rem.is_empty() {
        rem.push(zero);
    }
    rem
}

pub(super) fn poly_exact_div_q(a: &[Q], b: &[Q]) -> Vec<Q> {
    let zero = Q::from_integer(BigInt::from(0));
    let da = q_degree(a);
    let db = q_degree(b);
    if da < db {
        return vec![zero];
    }
    let mut rem = a.to_vec();
    let lc_b = b[db].clone();
    let dq = da - db;
    let mut quot = vec![zero.clone(); dq + 1];

    for i in (0..=dq).rev() {
        let lc_rem = if i + db < rem.len() {
            rem[i + db].clone()
        } else {
            zero.clone()
        };
        if lc_rem == zero {
            continue;
        }
        let q_coeff = lc_rem / lc_b.clone();
        quot[i] = q_coeff.clone();
        for (j, c) in b.iter().enumerate() {
            if i + j < rem.len() {
                rem[i + j] = rem[i + j].clone() - q_coeff.clone() * c.clone();
            }
        }
    }
    while quot.last().is_some_and(|c| *c == zero) {
        quot.pop();
    }
    if quot.is_empty() {
        quot.push(zero);
    }
    quot
}

/// Convert a Q polynomial to primitive BigInt coefficients by clearing denominators.
pub(super) fn q_poly_to_primitive_bigint(p: &[Q]) -> Vec<BigInt> {
    let zero = Q::from_integer(BigInt::from(0));
    // LCM of all denominators
    let lcm = p
        .iter()
        .filter(|c| **c != zero)
        .fold(BigInt::from(1), |acc, c| acc.lcm(c.denom()));

    let lcm_q = Q::from_integer(lcm);
    let mut result: Vec<BigInt> = p
        .iter()
        .map(|c| {
            let scaled = c * &lcm_q;
            debug_assert_eq!(scaled.denom(), &BigInt::from(1));
            scaled.to_integer()
        })
        .collect();
    while result.last().is_some_and(BigInt::is_zero) {
        result.pop();
    }
    if result.is_empty() {
        return vec![BigInt::zero()];
    }

    let mut content = BigInt::zero();
    for coeff in &result {
        if coeff.is_zero() {
            continue;
        }
        let abs_coeff = coeff.abs();
        content = if content.is_zero() {
            abs_coeff
        } else {
            content.gcd(&abs_coeff)
        };
    }
    if content > BigInt::from(1) {
        for coeff in &mut result {
            *coeff /= &content;
        }
    }
    result
}
