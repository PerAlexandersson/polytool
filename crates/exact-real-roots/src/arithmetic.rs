//! Integer coefficient arithmetic and primitive square-free parts.
//! Coefficients are in ascending degree order; helpers support the root backends.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

pub(super) fn trim(mut p: Vec<BigInt>) -> Vec<BigInt> {
    while p.last().is_some_and(|c| c.is_zero()) {
        p.pop();
    }
    p
}

pub(super) fn trim_slice(p: &[BigInt]) -> Vec<BigInt> {
    trim(p.to_vec())
}

pub(super) fn degree(p: &[BigInt]) -> Option<usize> {
    p.iter().rposition(|c| !c.is_zero())
}

pub(super) fn leading(p: &[BigInt]) -> BigInt {
    degree(p).map(|d| p[d].clone()).unwrap_or_else(BigInt::zero)
}

pub(super) fn sign_i8(x: &BigInt) -> i8 {
    if x.is_positive() {
        1
    } else if x.is_negative() {
        -1
    } else {
        0
    }
}

pub(super) fn content(p: &[BigInt]) -> BigInt {
    let mut g = BigInt::zero();
    for c in p {
        if !c.is_zero() {
            let a = c.abs();
            g = if g.is_zero() { a } else { g.gcd(&a) };
        }
    }
    g
}

pub(super) fn primitive_keep_sign(p: Vec<BigInt>) -> Vec<BigInt> {
    let p = trim(p);
    let g = content(&p);
    if g.is_zero() || g.is_one() {
        p
    } else {
        trim(p.into_iter().map(|c| c / &g).collect())
    }
}

pub(super) fn primitive_positive(p: Vec<BigInt>) -> Vec<BigInt> {
    let mut p = primitive_keep_sign(p);
    if leading(&p).is_negative() {
        for c in &mut p {
            *c = -c.clone();
        }
    }
    p
}

pub(super) fn poly_neg(p: &[BigInt]) -> Vec<BigInt> {
    p.iter().map(|c| -c).collect()
}

pub(super) fn poly_derivative(p: &[BigInt]) -> Vec<BigInt> {
    if p.len() <= 1 {
        return vec![];
    }
    trim(
        p.iter()
            .enumerate()
            .skip(1)
            .map(|(i, c)| c * BigInt::from(i))
            .collect(),
    )
}

pub(super) fn poly_mul(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r = vec![BigInt::zero(); a.len() + b.len() - 1];
    for (i, ca) in a.iter().enumerate() {
        if ca.is_zero() {
            continue;
        }
        for (j, cb) in b.iter().enumerate() {
            if !cb.is_zero() {
                r[i + j] += ca * cb;
            }
        }
    }
    trim(r)
}

/// Pseudo-remainder with a positive multiplier.
///
/// If `b_pos` is `b` multiplied by a sign so that its leading coefficient is
/// positive, this returns `lc(b_pos)^k * rem(a, b_pos)` for a nonnegative `k`.
/// The scalar is positive, so it does not change Sturm sign variations.
pub(super) fn pseudo_remainder_positive_multiplier(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let mut output = Vec::new();
    pseudo_remainder_positive_multiplier_into(&mut output, a, b);
    output
}

/// In-place-buffer variant of [`pseudo_remainder_positive_multiplier`].
///
/// Existing coefficient slots are updated with `clone_from`, allowing their
/// limb allocations to be reused by a streaming remainder sequence.
pub(super) fn pseudo_remainder_positive_multiplier_into(
    output: &mut Vec<BigInt>,
    a: &[BigInt],
    b: &[BigInt],
) {
    let da = match degree(a) {
        Some(d) => d,
        None => {
            output.clear();
            return;
        }
    };
    let db = degree(b).expect("pseudo_remainder: division by zero polynomial");
    if da < db {
        output.resize_with(a.len(), BigInt::zero);
        for (target, source) in output.iter_mut().zip(a) {
            target.clone_from(source);
        }
        output.truncate(da + 1);
        return;
    }

    let b = &b[..=db];
    let divisor_is_negative = b[db].is_negative();
    let lc = b[db].abs();
    debug_assert!(lc.is_positive());

    output.resize_with(da + 1, BigInt::zero);
    for (target, source) in output.iter_mut().zip(&a[..=da]) {
        target.clone_from(source);
    }
    output.truncate(da + 1);

    let mut e = da - db + 1;
    while let Some(dr) = degree(output) {
        if dr < db {
            break;
        }
        // The leading term cancels exactly in
        // `lc(b) * output - lc(output) * b`. Move its coefficient out rather
        // than cloning it, and leave that slot zero instead of performing the
        // two large-integer operations whose results would cancel.
        let c = std::mem::take(&mut output[dr]);
        let shift = dr - db;
        for coefficient in &mut output[..dr] {
            *coefficient *= &lc;
        }
        for (j, coefficient) in b[..db].iter().enumerate() {
            if !coefficient.is_zero() {
                if divisor_is_negative {
                    output[j + shift] += &c * coefficient;
                } else {
                    output[j + shift] -= &c * coefficient;
                }
            }
        }
        while output.last().is_some_and(BigInt::is_zero) {
            output.pop();
        }
        e -= 1;
    }
    if e > 0 && !output.is_empty() {
        let remaining_scale = lc.pow(e as u32);
        for coefficient in output.iter_mut() {
            *coefficient *= &remaining_scale;
        }
    }
}

fn primitive_gcd(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let mut r0 = primitive_positive(trim_slice(a));
    let mut r1 = primitive_positive(trim_slice(b));
    if r0.is_empty() {
        return r1;
    }
    if r1.is_empty() {
        return r0;
    }

    while !r1.is_empty() {
        let rem = primitive_positive(pseudo_remainder_positive_multiplier(&r0, &r1));
        r0 = r1;
        r1 = rem;
    }
    primitive_positive(r0)
}

pub(super) fn exact_div(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let da = match degree(a) {
        Some(d) => d,
        None => return vec![],
    };
    let db = degree(b).expect("exact_div: division by zero polynomial");
    if da < db {
        return vec![];
    }

    let mut rem = trim_slice(a);
    let mut q = vec![BigInt::zero(); da - db + 1];
    let lc_b = leading(b);
    for k in (0..=da - db).rev() {
        let coeff = rem[db + k].clone();
        if coeff.is_zero() {
            continue;
        }
        assert!(
            (&coeff % &lc_b).is_zero(),
            "exact_div: non-exact leading quotient"
        );
        let qk = coeff / &lc_b;
        q[k] = qk.clone();
        for (j, bj) in b.iter().enumerate().take(db + 1) {
            rem[j + k] -= &qk * bj;
        }
    }
    debug_assert!(trim(rem).is_empty(), "exact_div: nonzero remainder");
    trim(q)
}

pub(super) fn strip_initial_zeros_bigint(coeffs: &[BigInt]) -> &[BigInt] {
    match coeffs.iter().position(|c| !c.is_zero()) {
        Some(i) => &coeffs[i..],
        None => &[],
    }
}

/// Return the primitive, positive-leading square-free part.
///
/// This lower-level entry point is exposed so polynomial libraries can reuse
/// the same fraction-free integer GCD as the root-counting backends.
#[doc(hidden)]
pub fn squarefree_part_primitive_bigint_coeffs(coeffs: &[BigInt]) -> Vec<BigInt> {
    let p = primitive_positive(trim_slice(coeffs));
    if degree(&p).unwrap_or(0) == 0 {
        return p;
    }
    let dp = poly_derivative(&p);
    if dp.is_empty() {
        return p;
    }
    let g = primitive_gcd(&p, &dp);
    if degree(&g).unwrap_or(0) == 0 {
        p
    } else {
        primitive_positive(exact_div(&p, &g))
    }
}

/// Try to reduce a BigInt vector to i64 coefficients.
pub fn bigint_coeffs_to_i64(coeffs: &[BigInt]) -> Option<Vec<i64>> {
    coeffs.iter().map(ToPrimitive::to_i64).collect()
}
