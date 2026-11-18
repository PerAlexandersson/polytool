//! Minimal integer PRS glue for the real/imaginary gcd and multiplicities.
//!
//! The root-counting and squarefree kernels remain in exact-real-roots. Its
//! internal two-polynomial gcd is not public, so the split uses a primitive
//! pseudo-remainder gcd here, without depending on the consumer polynomial type.

use super::RealRootLayer;
use crate::BigRational as Q;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

fn trim(p: &mut Vec<BigInt>) {
    while p.last().is_some_and(BigInt::is_zero) {
        p.pop();
    }
}

fn primitive(mut p: Vec<BigInt>) -> Vec<BigInt> {
    trim(&mut p);
    let content = p.iter().fold(BigInt::zero(), |a, b| a.gcd(b));
    if !content.is_zero() {
        let divisor = if p.last().expect("nonzero polynomial").is_negative() {
            -content
        } else {
            content
        };
        for x in &mut p {
            *x /= &divisor;
        }
    }
    p
}

pub(super) fn gcd(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let mut a = primitive(a.to_vec());
    let mut b = primitive(b.to_vec());
    while !b.is_empty() {
        while a.len() >= b.len() {
            let shift = a.len() - b.len();
            let leading = a.last().expect("nonzero remainder").clone();
            let multiplier = b.last().expect("nonzero divisor");
            for x in &mut a {
                *x *= multiplier;
            }
            for (i, x) in b.iter().enumerate() {
                a[shift + i] -= &leading * x;
            }
            a = primitive(a);
        }
        (a, b) = (b, a);
    }
    a
}

pub(super) fn clear_pair(a: &[Q], b: &[Q]) -> (Vec<BigInt>, Vec<BigInt>) {
    let denominator = a
        .iter()
        .chain(b)
        .fold(BigInt::one(), |d, x| d.lcm(x.denom()));
    let clear = |p: &[Q]| {
        let mut result: Vec<BigInt> = p
            .iter()
            .map(|x| x.numer() * (&denominator / x.denom()))
            .collect();
        trim(&mut result);
        result
    };
    (clear(a), clear(b))
}

pub(super) fn divide_exact(a: &[BigInt], b: &[BigInt]) -> Vec<Q> {
    assert!(!b.is_empty(), "nonzero divisor");
    if a.len() < b.len() {
        assert!(a.is_empty(), "exact polynomial divisibility");
        return Vec::new();
    }
    let mut remainder: Vec<Q> = a.iter().cloned().map(Q::from_integer).collect();
    let divisor: Vec<Q> = b.iter().cloned().map(Q::from_integer).collect();
    let mut quotient = vec![Q::zero(); a.len() - b.len() + 1];
    for shift in (0..quotient.len()).rev() {
        let q = &remainder[shift + b.len() - 1] / divisor.last().expect("nonzero divisor");
        for (i, x) in divisor.iter().enumerate() {
            remainder[shift + i] -= &q * x;
        }
        quotient[shift] = q;
    }
    assert!(
        remainder.iter().all(Q::is_zero),
        "exact polynomial divisibility"
    );
    quotient
}

pub(super) fn count_real_roots_with_multiplicity(p: &[BigInt]) -> (usize, Vec<RealRootLayer>) {
    let mut p = primitive(p.to_vec());
    let mut count = 0;
    let mut layers = Vec::new();
    // The j-th layer contains precisely the roots of multiplicity >= j.
    // Dividing by the existing primitive squarefree part lowers every
    // remaining multiplicity by one, so summing distinct counts is exact.
    while p.len() > 1 {
        let distinct_real_roots = exact_real_roots::count_real_roots_prs_bigint_coeffs(&p);
        count += distinct_real_roots;
        let squarefree = exact_real_roots::squarefree_part_primitive_bigint_coeffs(&p);
        let quotient = divide_exact(&p, &squarefree);
        layers.push(RealRootLayer {
            polynomial: p,
            distinct_real_roots,
        });
        p = primitive(clear_pair(&quotient, &[]).0);
    }
    (count, layers)
}
