//! Exact coefficient filters, low-degree formulas and adaptive real-rootedness.
//! The general fallback delegates to the signed PRS backend.

use crate::arithmetic::{degree, strip_initial_zeros_bigint, trim, trim_slice};
use crate::prs::{derivative_sturm_summary, is_real_rooted_prs_bigint_coeffs};
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

pub(super) fn alternating_neg_argument(coeffs: &[BigInt]) -> Vec<BigInt> {
    trim(
        coeffs
            .iter()
            .enumerate()
            .map(|(i, c)| if i % 2 == 0 { c.clone() } else { -c })
            .collect(),
    )
}

fn all_one_sign(coeffs: &[BigInt]) -> bool {
    let mut seen_pos = false;
    let mut seen_neg = false;
    for c in coeffs {
        if c.is_positive() {
            seen_pos = true;
        } else if c.is_negative() {
            seen_neg = true;
        }
    }
    !(seen_pos && seen_neg)
}

/// Return whether a polynomial of degree at most four has only real roots.
///
/// The discriminant formulas below are evaluated only up to a positive
/// factor, since only their sign is needed.  The quartic invariants also cover
/// zero-discriminant cases, so repeated roots do not require a GCD fallback.
pub(super) fn low_degree_is_real_rooted(coeffs: &[BigInt]) -> Option<bool> {
    let p = trim_slice(coeffs);
    let Some(d) = degree(&p) else {
        return Some(true);
    };

    match d {
        0 | 1 => Some(true),
        2 => {
            // b^2 - 4ac >= 0.
            let discriminant = &p[1] * &p[1] - BigInt::from(4) * &p[2] * &p[0];
            Some(!discriminant.is_negative())
        }
        3 => {
            // For ax^3 + bx^2 + cx + d, the usual discriminant has the
            // same sign as 4 Delta_0^3 - Delta_1^2, since their ratio is
            // the positive quantity 27a^2.
            let delta_0 = &p[2] * &p[2] - BigInt::from(3) * &p[3] * &p[1];
            let delta_1 = BigInt::from(2) * &p[2] * &p[2] * &p[2]
                - BigInt::from(9) * &p[3] * &p[2] * &p[1]
                + BigInt::from(27) * &p[3] * &p[3] * &p[0];
            let discriminant_numerator =
                BigInt::from(4) * &delta_0 * &delta_0 * &delta_0 - &delta_1 * &delta_1;
            Some(!discriminant_numerator.is_negative())
        }
        4 => {
            // Binary-quartic invariants.  The usual discriminant has the
            // same sign as 4 I^3 - J^2 (the ratio is the positive constant
            // 27).  When it is positive, the quartic has four real roots
            // exactly when P < 0 and D < 0.
            let i =
                BigInt::from(12) * &p[4] * &p[0] - BigInt::from(3) * &p[3] * &p[1] + &p[2] * &p[2];
            let j = BigInt::from(72) * &p[4] * &p[2] * &p[0]
                + BigInt::from(9) * &p[3] * &p[2] * &p[1]
                - BigInt::from(27) * &p[4] * &p[1] * &p[1]
                - BigInt::from(27) * &p[3] * &p[3] * &p[0]
                - BigInt::from(2) * &p[2] * &p[2] * &p[2];
            let discriminant_numerator = BigInt::from(4) * &i * &i * &i - &j * &j;

            if discriminant_numerator.is_negative() {
                return Some(false);
            }

            let quartic_p = BigInt::from(8) * &p[4] * &p[2] - BigInt::from(3) * &p[3] * &p[3];
            if quartic_p.is_positive()
                || (discriminant_numerator.is_positive() && quartic_p.is_zero())
            {
                return Some(false);
            }

            let quartic_d = BigInt::from(64) * &p[4] * &p[4] * &p[4] * &p[0]
                - BigInt::from(16) * &p[4] * &p[4] * &p[2] * &p[2]
                + BigInt::from(16) * &p[4] * &p[3] * &p[3] * &p[2]
                - BigInt::from(16) * &p[4] * &p[4] * &p[3] * &p[1]
                - BigInt::from(3) * &p[3] * &p[3] * &p[3] * &p[3];

            if discriminant_numerator.is_positive() {
                Some(quartic_d.is_negative())
            } else {
                // Strict inequalities give one double and two simple real
                // roots. Equality includes two double roots and the triple-
                // or quadruple-root cases.
                Some(!quartic_p.is_positive() && !quartic_d.is_positive())
            }
        }
        _ => None,
    }
}

fn make_nonnegative(mut coeffs: Vec<BigInt>) -> Vec<BigInt> {
    if coeffs.iter().any(|c| c.is_negative()) {
        for c in &mut coeffs {
            *c = -c.clone();
        }
    }
    coeffs
}

/// Check Newton's inequalities
/// `a_k^2 binom(d,k-1) binom(d,k+1) >= a_{k-1} a_{k+1} binom(d,k)^2`.
///
/// They hold for every real-rooted polynomial with real coefficients, of any
/// signs (the products `a_{k-1} a_{k+1}` and `a_k^2` are sign-invariant under
/// `t -> -t` and scaling), so this is a cheap necessary condition for
/// real-rootedness.  It is exact and works with arbitrary-size coefficients.
pub fn satisfies_newton_inequalities_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let p = trim_slice(strip_initial_zeros_bigint(coeffs));
    let d = match degree(&p) {
        Some(d) if d >= 2 => d,
        _ => return true,
    };

    let mut binom = vec![BigInt::one(); d + 1];
    for k in 1..=d {
        binom[k] = &binom[k - 1] * BigInt::from(d - k + 1) / BigInt::from(k);
    }

    for k in 1..d {
        let lhs = p[k].pow(2) * &binom[k - 1] * &binom[k + 1];
        let rhs = &p[k - 1] * &p[k + 1] * &binom[k] * &binom[k];
        if lhs < rhs {
            return false;
        }
    }
    true
}

/// Deprecated name of [`satisfies_newton_inequalities_bigint_coeffs`].
#[deprecated(
    since = "0.2.0",
    note = "renamed to `satisfies_newton_inequalities_bigint_coeffs`"
)]
pub fn satisfies_newton_inequalities_bigint(coeffs: &[BigInt]) -> bool {
    satisfies_newton_inequalities_bigint_coeffs(coeffs)
}

/// Check Kurtz's sufficient condition for distinct real roots.
///
/// If `f(t)=sum a_i t^i` has positive coefficients and
/// `a_i^2 > 4 a_{i-1} a_{i+1}` for every `1 <= i < d`, then all roots of `f`
/// are distinct and real.  Since all coefficients are positive, those roots
/// are automatically negative.  The condition is sufficient, not necessary.
///
/// This returns `false` for degree `0` and `1`, where the criterion has no
/// inequalities to check.
pub fn satisfies_kurtz_condition_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let p = trim_slice(strip_initial_zeros_bigint(coeffs));
    let d = match degree(&p) {
        Some(d) if d >= 2 => d,
        _ => return false,
    };
    if p.iter().take(d + 1).any(|c| !c.is_positive()) {
        return false;
    }

    for i in 1..d {
        if p[i].pow(2) <= BigInt::from(4) * &p[i - 1] * &p[i + 1] {
            return false;
        }
    }
    true
}

/// Deprecated name of [`satisfies_kurtz_condition_bigint_coeffs`].
#[deprecated(
    since = "0.2.0",
    note = "renamed to `satisfies_kurtz_condition_bigint_coeffs`"
)]
pub fn satisfies_kurtz_condition_bigint(coeffs: &[BigInt]) -> bool {
    satisfies_kurtz_condition_bigint_coeffs(coeffs)
}

/// Exact real-rootedness optimized for one-signed coefficient polynomials.
///
/// Returns `None` if the nonzero coefficients do not all have the same sign.
/// Powers of `t` are removed first, since zero is already a real root.  After
/// that, polynomials of degree at most 1 always give `Some(true)`, whatever the
/// signs.
pub fn is_real_rooted_one_signed_bigint_coeffs(coeffs: &[BigInt]) -> Option<bool> {
    let p = trim_slice(strip_initial_zeros_bigint(coeffs));
    let d = match degree(&p) {
        Some(d) => d,
        None => return Some(true),
    };
    if d <= 1 {
        return Some(true);
    }
    if !all_one_sign(&p) {
        return None;
    }
    if let Some(real_rooted) = low_degree_is_real_rooted(&p) {
        return Some(real_rooted);
    }

    let p = make_nonnegative(p);
    if satisfies_kurtz_condition_bigint_coeffs(&p) {
        return Some(true);
    }
    if !satisfies_newton_inequalities_bigint_coeffs(&p) {
        return Some(false);
    }

    // After stripping zero roots, a one-signed polynomial has only negative
    // real roots, so counting all distinct real roots (signs at +-infinity
    // only) is equivalent to, and cheaper than, counting positive roots of
    // p(-t). The subresultant PRS beats the Uspensky/Descartes path on every
    // benchmarked one-signed family (Eulerian, type B Eulerian, Touchard,
    // Narayana, linear-factor products; see `polytool bench real-rootedness`).
    let summary = derivative_sturm_summary(&p, true, false);
    Some(summary.real_roots == summary.squarefree_degree)
}

/// Exact real-rootedness with an adaptive one-signed path and PRS fallback.
///
/// The one-signed path first applies the Kurtz (sufficient) and Newton
/// (necessary) coefficient filters, then counts distinct real roots with the
/// subresultant PRS; mixed-sign inputs use PRS directly. (An earlier
/// heuristic switched large palindromic inputs to Uspensky; after the
/// Collins--Brown PRS it was about 10--26 times slower on those families.) This avoids constructing Bézout/PSD
/// matrices entirely and is the backend used by the public
/// `real_rootedness::is_real_rooted_bigint_coeffs` default.
pub fn is_real_rooted_fast_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    if let Some(real_rooted) = low_degree_is_real_rooted(coeffs) {
        return real_rooted;
    }
    if let Some(rr) = is_real_rooted_one_signed_bigint_coeffs(coeffs) {
        rr
    } else {
        is_real_rooted_prs_bigint_coeffs(coeffs)
    }
}

/// Convenience wrapper for `i64` coefficient vectors.
pub fn is_real_rooted_fast_i64(coeffs: &[i64]) -> bool {
    let coeffs: Vec<BigInt> = coeffs.iter().map(|&c| BigInt::from(c)).collect();
    is_real_rooted_fast_bigint_coeffs(&coeffs)
}
