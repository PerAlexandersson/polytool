use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};

pub type Q = BigRational;

pub fn trim_q(mut p: Vec<Q>) -> Vec<Q> {
    while p.last().is_some_and(Zero::is_zero) {
        p.pop();
    }
    p
}

pub fn is_one_q(p: &[Q]) -> bool {
    p.len() == 1 && p[0].is_one()
}

pub fn derivative_q(p: &[Q]) -> Vec<Q> {
    if p.len() <= 1 {
        return Vec::new();
    }
    trim_q(
        p.iter()
            .enumerate()
            .skip(1)
            .map(|(i, c)| c * BigInt::from(i))
            .collect(),
    )
}

pub fn add_q(a: &[Q], b: &[Q]) -> Vec<Q> {
    let mut result = vec![Q::zero(); a.len().max(b.len())];
    for (i, value) in a.iter().enumerate() {
        result[i] += value;
    }
    for (i, value) in b.iter().enumerate() {
        result[i] += value;
    }
    trim_q(result)
}

pub fn sub_q(a: &[Q], b: &[Q]) -> Vec<Q> {
    let mut result = vec![Q::zero(); a.len().max(b.len())];
    for (i, value) in a.iter().enumerate() {
        result[i] += value;
    }
    for (i, value) in b.iter().enumerate() {
        result[i] -= value;
    }
    trim_q(result)
}

pub fn mul_q(a: &[Q], b: &[Q]) -> Vec<Q> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut result = vec![Q::zero(); a.len() + b.len() - 1];
    for (i, left) in a.iter().enumerate() {
        for (j, right) in b.iter().enumerate() {
            result[i + j] += left * right;
        }
    }
    trim_q(result)
}

pub fn scale_q(a: &[Q], scalar: &Q) -> Vec<Q> {
    trim_q(a.iter().map(|value| value * scalar).collect())
}

pub fn div_rem_q(dividend: &[Q], divisor: &[Q]) -> (Vec<Q>, Vec<Q>) {
    let divisor = trim_q(divisor.to_vec());
    assert!(!divisor.is_empty(), "division by zero polynomial");
    let mut remainder = trim_q(dividend.to_vec());
    if remainder.len() < divisor.len() {
        return (Vec::new(), remainder);
    }
    let mut quotient = vec![Q::zero(); remainder.len() - divisor.len() + 1];
    let divisor_leading = divisor.last().expect("nonempty divisor").clone();
    while remainder.len() >= divisor.len() && !remainder.is_empty() {
        let shift = remainder.len() - divisor.len();
        let multiple = remainder.last().expect("nonempty remainder").clone() / &divisor_leading;
        quotient[shift] = multiple.clone();
        for (i, value) in divisor.iter().enumerate() {
            remainder[i + shift] -= &multiple * value;
        }
        remainder = trim_q(remainder);
    }
    (trim_q(quotient), remainder)
}

pub fn exact_div_q(dividend: &[Q], divisor: &[Q]) -> Vec<Q> {
    let (quotient, remainder) = div_rem_q(dividend, divisor);
    assert!(remainder.is_empty(), "polynomial division was not exact");
    quotient
}

pub fn monic_q(p: &[Q]) -> Vec<Q> {
    let p = trim_q(p.to_vec());
    if p.is_empty() {
        return p;
    }
    let leading = p.last().expect("nonempty polynomial").clone();
    scale_q(&p, &(&Q::one() / leading))
}

pub fn gcd_q(a: &[Q], b: &[Q]) -> Vec<Q> {
    let mut left = monic_q(a);
    let mut right = monic_q(b);
    while !right.is_empty() {
        let (_, remainder) = div_rem_q(&left, &right);
        left = right;
        right = remainder;
    }
    monic_q(&left)
}

pub fn evaluate_q(p: &[Q], x: &Q) -> Q {
    p.iter()
        .rev()
        .fold(Q::zero(), |value, coefficient| value * x + coefficient)
}

pub fn primitive_positive(p: &[Q]) -> (Vec<BigInt>, Q) {
    let p = trim_q(p.to_vec());
    assert!(!p.is_empty(), "a factor cannot be zero");
    let denominator_lcm = p.iter().fold(BigInt::one(), |lcm, coefficient| {
        lcm.lcm(coefficient.denom())
    });
    let mut integer_coeffs: Vec<BigInt> = p
        .iter()
        .map(|coefficient| coefficient.numer() * (&denominator_lcm / coefficient.denom()))
        .collect();
    let content = integer_coeffs
        .iter()
        .filter(|coefficient| !coefficient.is_zero())
        .fold(BigInt::zero(), |content, coefficient| {
            content.gcd(&coefficient.abs())
        });
    if !content.is_zero() {
        for coefficient in &mut integer_coeffs {
            *coefficient /= &content;
        }
    }
    let mut scalar = Q::new(content, denominator_lcm);
    if integer_coeffs
        .last()
        .expect("nonempty polynomial")
        .is_negative()
    {
        for coefficient in &mut integer_coeffs {
            *coefficient = -coefficient.clone();
        }
        scalar = -scalar;
    }
    (integer_coeffs, scalar)
}

pub fn mod_bigint(value: &BigInt, modulus: &BigInt) -> BigInt {
    let residue = value % modulus;
    if residue.is_negative() {
        residue + modulus
    } else {
        residue
    }
}

pub fn centered_mod_bigint(value: &BigInt, modulus: &BigInt) -> BigInt {
    let residue = mod_bigint(value, modulus);
    let half = modulus / 2u32;
    if residue > half {
        residue - modulus
    } else {
        residue
    }
}

pub fn mod_inverse_bigint(value: &BigInt, modulus: &BigInt) -> BigInt {
    let egcd = value.extended_gcd(modulus);
    assert!(egcd.gcd.is_one(), "modular inverse does not exist");
    mod_bigint(&egcd.x, modulus)
}
