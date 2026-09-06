//! Real-rootedness testing, interlacing, log-concavity, and ultra-log-concavity.
//!
//! All functions accept polynomials as `i64` coefficient vectors in ascending degree order:
//! `coeffs[i]` is the coefficient of t^i.
//!
//! # Default algorithm: primitive integer root counting
//!
//! The default real-rootedness functions use the primitive integer PRS path in
//! [`crate::root_count`].  For one-signed coefficient polynomials, this first
//! tries cheap coefficient tests and then counts positive roots after the
//! transformation `t ↦ -t`; otherwise it counts all real roots of the
//! squarefree part directly.  Interlacing uses Bézout matrices, which are
//! dramatically faster than older rational Sturm chains for that task.
//!
//! ## Bézout matrix for interlacing
//!
//! Given polynomials f (degree d) and g (degree d−1), both with positive leading
//! coefficients, the **Bézout matrix** B(f, g) is the d×d symmetric matrix
//! whose (i, j) entry is the coefficient of x^i y^j in the bivariate polynomial
//!
//! ```text
//! (f(x)g(y) - f(y)g(x)) / (x - y).
//! ```
//!
//! There are two standard oriented versions.
//!
//! - If `deg(f) = deg(g) + 1`, then `g` interlaces `f` iff `B(f, g)` is
//!   positive semidefinite; in the coprime/strict case this is positive
//!   definite. See Kummer--Naldi--Plaumann, Thm. 2.13, citing
//!   Krein--Naimark, §2.2.
//! - If `deg(f) = deg(g)`, the same-degree alternation criterion uses the
//!   same Bézoutian with the orientation fixed by the argument order. Fisk,
//!   §9.21, Cor. 9.145 gives the positive definite/no-common-root form.
//!
//! In semidefinite/common-root cases, the common factor must be handled
//! separately: a positive semidefinite Bézoutian does not by itself certify
//! that a shared factor has only real roots. This is why
//! [`check_weak_interlacing`] divides out the gcd and explicitly verifies that
//! the gcd is real-rooted.
//!
//! This reduces interlacing to a single exact matrix definiteness check,
//! avoiding root isolation entirely.
//!
//! ## Bézout matrix for real-rootedness
//!
//! A polynomial f of degree d is real-rooted if and only if f' interlaces f. Since
//! f may have repeated roots, B(f, f') is only positive **semi**-definite for real-rooted f.
//!
//! ## Shared roots and same-degree interlacing
//!
//! The Bézout matrix detects **strict** interlacing only (deg difference exactly 1).
//! [`check_weak_interlacing`] extends this to handle:
//!
//! - **Shared roots**: divides out gcd(f, g) over ℚ, verifies the GCD is real-rooted,
//!   then checks strict interlacing of the reduced polynomials.
//! - **Same degree**: `check_weak_interlacing(f, g)` tests f ≪ g (f on the LEFT)
//!   by extending f with a root far to the right (via the Cauchy bound), reducing
//!   to the deg+1 case. If all coefficients of f are positive (all roots negative),
//!   multiplying by t suffices.
//! - **Directed**: `check_weak_interlacing(f, g)` tests `f ≪ g` (f interlaces g
//!   from the left). This is **not symmetric**: `f ≪ g` and `g ≪ f` are different
//!   for same-degree polynomials. Requires deg(f) ≤ deg(g) ≤ deg(f) + 1.
//!
//! ## Sturm chain fallback
//!
//! The original Sturm-chain implementations are available as [`is_real_rooted_sturm`],
//! [`check_interlacing_sturm`], and [`real_roots`] for cases where actual root locations
//! are needed.
//!
//! References:
//!
//! - S. Fisk, *Polynomials, roots, and interlacing*, arXiv:math/0612833,
//!   §9.21, Cor. 9.145.
//! - M. Kummer, S. Naldi, and D. Plaumann, *Spectrahedral representations of
//!   plane hyperbolic curves*, arXiv:1807.10901, Thm. 2.13.
//! - M. G. Krein and M. A. Naimark, *The method of symmetric and Hermitian
//!   forms in the theory of the separation of the roots of algebraic
//!   equations*, Linear and Multilinear Algebra 10 (1981), §2.2.
//! - MathOverflow discussion of the common-factor caveat:
//!   <https://mathoverflow.net/questions/403708/bezout-matrices-and-interlacing-roots>

use crate::sturm::SturmChain;
use crate::Polynomial;
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{One, Signed, ToPrimitive, Zero};

type Q = Ratio<BigInt>;

// ---------------------------------------------------------------------------
// Display utility for &[i64] coefficient vectors
// ---------------------------------------------------------------------------

/// Format an `i64` coefficient vector as a human-readable polynomial string.
///
/// Uses the variable name `t` by default. Coefficients of 1 or -1 are displayed
/// without the leading "1". Examples:
///
/// - `[1, 2, 1]` → `"1 + 2t + t^2"`
/// - `[0, -1, 0, 3]` → `"-t + 3t^3"`
/// - `[]` or `[0]` → `"0"`
pub fn format_poly(coeffs: &[i64]) -> String {
    format_poly_var(coeffs, "t")
}

/// Format an `i64` coefficient vector with a custom variable name.
pub fn format_poly_var(coeffs: &[i64], var: &str) -> String {
    let mut terms = Vec::new();
    for (i, &c) in coeffs.iter().enumerate() {
        if c == 0 {
            continue;
        }
        let term = match (c, i) {
            (_, 0) => format!("{}", c),
            (1, 1) => var.to_string(),
            (-1, 1) => format!("-{}", var),
            (_, 1) => format!("{}{}", c, var),
            (1, e) => format!("{}^{}", var, e),
            (-1, e) => format!("-{}^{}", var, e),
            (_, e) => format!("{}{}^{}", c, var, e),
        };
        terms.push(term);
    }
    if terms.is_empty() {
        return "0".to_string();
    }
    let mut result = terms[0].clone();
    for term in &terms[1..] {
        if let Some(rest) = term.strip_prefix('-') {
            result.push_str(" - ");
            result.push_str(rest);
        } else {
            result.push_str(" + ");
            result.push_str(term);
        }
    }
    result
}

/// Format a `BigInt` coefficient vector as a human-readable polynomial string.
pub fn format_poly_bigint_coeffs(coeffs: &[BigInt]) -> String {
    format_poly_bigint_var(coeffs, "t")
}

/// Format a `BigInt` coefficient vector with a custom variable name.
pub fn format_poly_bigint_var(coeffs: &[BigInt], var: &str) -> String {
    let one = BigInt::from(1);
    let minus_one = BigInt::from(-1);
    let mut terms = Vec::new();
    for (i, c) in coeffs.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let term = match i {
            0 => c.to_string(),
            1 if c == &one => var.to_string(),
            1 if c == &minus_one => format!("-{}", var),
            1 => format!("{}{}", c, var),
            e if c == &one => format!("{}^{}", var, e),
            e if c == &minus_one => format!("-{}^{}", var, e),
            e => format!("{}{}^{}", c, var, e),
        };
        terms.push(term);
    }
    if terms.is_empty() {
        return "0".to_string();
    }
    let mut result = terms[0].clone();
    for term in &terms[1..] {
        if let Some(rest) = term.strip_prefix('-') {
            result.push_str(" - ");
            result.push_str(rest);
        } else {
            result.push_str(" + ");
            result.push_str(term);
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Real-rootedness and interlacing
// ---------------------------------------------------------------------------

/// Check if a polynomial is real-rooted using Sturm chains with exact arithmetic.
///
/// This is the slower rational Sturm-chain method. Prefer
/// [`is_real_rooted`] unless you need the older Sturm-chain route explicitly.
///
/// Returns true for the zero polynomial and for constant/linear polynomials.
pub fn is_real_rooted_sturm(coeffs: &[i64]) -> bool {
    let Some(start) = coeffs.iter().position(|&c| c != 0) else {
        return true;
    };
    let end = coeffs
        .iter()
        .rposition(|&c| c != 0)
        .expect("nonzero coefficient must have a last position");
    let trimmed = &coeffs[start..=end];
    if trimmed.len() <= 2 {
        return true;
    }

    let sc = SturmChain::from_i64_coeffs(trimmed);
    let sf_degree = sc.square_free_degree();
    sc.count_real_roots() == sf_degree
}

/// Check if a `BigInt`-coefficient polynomial is real-rooted using Sturm
/// chains with exact rational arithmetic.
pub fn is_real_rooted_sturm_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let Some(start) = coeffs.iter().position(|c| !c.is_zero()) else {
        return true;
    };
    let end = coeffs
        .iter()
        .rposition(|c| !c.is_zero())
        .expect("nonzero coefficient must have a last position");
    let trimmed = &coeffs[start..=end];
    if trimmed.len() <= 2 {
        return true;
    }

    let sc = SturmChain::from_bigint_coeffs(trimmed);
    let sf_degree = sc.square_free_degree();
    sc.count_real_roots() == sf_degree
}

/// Check if a coefficient sequence is log-concave.
///
/// A sequence a_0, a_1, ..., a_d is log-concave if a_k^2 >= a_{k-1} * a_{k+1}
/// for all 1 <= k <= d-1.
pub fn is_log_concave(coeffs: &[i64]) -> bool {
    if coeffs.len() <= 2 {
        return true;
    }
    for i in 1..coeffs.len() - 1 {
        let lhs = i128::from(coeffs[i]) * i128::from(coeffs[i]);
        let rhs = i128::from(coeffs[i - 1]) * i128::from(coeffs[i + 1]);
        if lhs < rhs {
            return false;
        }
    }
    true
}

/// Check if a `BigInt` coefficient sequence is log-concave.
pub fn is_log_concave_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    if coeffs.len() <= 2 {
        return true;
    }
    for i in 1..coeffs.len() - 1 {
        let lhs = &coeffs[i] * &coeffs[i];
        let rhs = &coeffs[i - 1] * &coeffs[i + 1];
        if lhs < rhs {
            return false;
        }
    }
    true
}

/// Check if a coefficient sequence is unimodal.
///
/// A sequence is unimodal if it is weakly increasing up to some index and
/// weakly decreasing afterwards.  Plateaus are allowed.
pub fn is_unimodal(coeffs: &[i64]) -> bool {
    if coeffs.len() <= 2 {
        return true;
    }

    let mut i = 1;
    while i < coeffs.len() && coeffs[i - 1] <= coeffs[i] {
        i += 1;
    }
    while i < coeffs.len() && coeffs[i - 1] >= coeffs[i] {
        i += 1;
    }
    i == coeffs.len()
}

/// Check if a `BigInt` coefficient sequence is unimodal.
pub fn is_unimodal_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    if coeffs.len() <= 2 {
        return true;
    }

    let mut i = 1;
    while i < coeffs.len() && coeffs[i - 1] <= coeffs[i] {
        i += 1;
    }
    while i < coeffs.len() && coeffs[i - 1] >= coeffs[i] {
        i += 1;
    }
    i == coeffs.len()
}

/// Check if a coefficient sequence is ultra-log-concave (satisfies Newton's inequalities).
///
/// A polynomial a_0 + a_1 t + ... + a_d t^d of degree d is ultra-log-concave if
/// a_k / C(d,k) is log-concave, where C(d,k) is the binomial coefficient.
///
/// This is implied by real-rootedness but is strictly weaker.
pub fn is_ultra_log_concave(coeffs: &[i64]) -> bool {
    let d = match coeffs.iter().rposition(|&c| c != 0) {
        Some(d) => d,
        None => return true,
    };
    if d <= 1 {
        return true;
    }
    let mut binom = vec![BigInt::from(1); d + 1];
    for k in 1..=d {
        binom[k] = &binom[k - 1] * BigInt::from(d - k + 1) / BigInt::from(k);
    }
    for k in 1..d {
        let lhs = BigInt::from(coeffs[k]).pow(2) * &binom[k - 1] * &binom[k + 1];
        let rhs = BigInt::from(coeffs[k - 1]) * BigInt::from(coeffs[k + 1]) * &binom[k] * &binom[k];
        if lhs < rhs {
            return false;
        }
    }
    true
}

/// Check if a `BigInt` coefficient sequence is ultra-log-concave.
pub fn is_ultra_log_concave_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let d = match coeffs.iter().rposition(|c| !c.is_zero()) {
        Some(d) => d,
        None => return true,
    };
    if d <= 1 {
        return true;
    }
    let mut binom = vec![BigInt::from(1); d + 1];
    for k in 1..=d {
        binom[k] = &binom[k - 1] * BigInt::from(d - k + 1) / BigInt::from(k);
    }
    for k in 1..d {
        let lhs = coeffs[k].pow(2) * &binom[k - 1] * &binom[k + 1];
        let rhs = &coeffs[k - 1] * &coeffs[k + 1] * &binom[k] * &binom[k];
        if lhs < rhs {
            return false;
        }
    }
    true
}

/// Check if a polynomial has palindromic (symmetric) coefficients.
///
/// A polynomial a_0 + a_1 t + ... + a_d t^d is palindromic if a_i = a_{d-i}
/// for all i, equivalently t^d p(1/t) = p(t).
///
/// Returns true for the zero polynomial and constant polynomials.
/// Trailing zeros are trimmed before checking.
pub fn is_palindromic(coeffs: &[i64]) -> bool {
    let d = match coeffs.iter().rposition(|&c| c != 0) {
        Some(d) => d,
        None => return true,
    };
    for i in 0..=d / 2 {
        if coeffs[i] != coeffs[d - i] {
            return false;
        }
    }
    true
}

/// Check if a `BigInt` coefficient polynomial is palindromic.
pub fn is_palindromic_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let d = match coeffs.iter().rposition(|c| !c.is_zero()) {
        Some(d) => d,
        None => return true,
    };
    for i in 0..=d / 2 {
        if coeffs[i] != coeffs[d - i] {
            return false;
        }
    }
    true
}

/// Compute the gamma coefficients of a palindromic polynomial.
///
/// For a palindromic polynomial p(t) of degree d, there is a unique expansion
///
/// ```text
/// p(t) = sum_{i=0}^{floor(d/2)} gamma_i * t^i * (1+t)^{d-2i}
/// ```
///
/// Returns `None` if the polynomial is not palindromic or if an exact gamma
/// coefficient does not fit in `i64`; use
/// [`gamma_coefficients_bigint_coeffs`] for the no-overflow result.
/// Returns `Some(gamma_coefficients)` where `gamma[i]` is the coefficient
/// of t^i (1+t)^{d-2i} in the expansion.
///
/// The polynomial is gamma-positive iff all returned coefficients are non-negative.
pub fn gamma_coefficients(coeffs: &[i64]) -> Option<Vec<i64>> {
    let coeffs = coeffs.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
    gamma_coefficients_bigint_coeffs(&coeffs)?
        .into_iter()
        .map(|gamma| gamma.to_i64())
        .collect()
}

/// Compute gamma coefficients for a palindromic `BigInt` coefficient polynomial.
pub fn gamma_coefficients_bigint_coeffs(coeffs: &[BigInt]) -> Option<Vec<BigInt>> {
    let d = match coeffs.iter().rposition(|c| !c.is_zero()) {
        Some(d) => d,
        None => return Some(vec![]),
    };
    if !is_palindromic_bigint_coeffs(coeffs) {
        return None;
    }
    let half = d / 2;

    let mut binomials = vec![vec![BigInt::from(0); d + 1]; d + 1];
    for n in 0..=d {
        binomials[n][0] = BigInt::from(1);
        for k in 1..=n {
            binomials[n][k] = &binomials[n - 1][k - 1] + &binomials[n - 1][k];
        }
    }

    let mut gamma = vec![BigInt::from(0); half + 1];
    for i in 0..=half {
        let mut val = coeffs[i].clone();
        for j in 0..i {
            val -= &gamma[j] * &binomials[d - 2 * j][i - j];
        }
        gamma[i] = val;
    }

    Some(gamma)
}

/// Check if a polynomial is gamma-positive.
///
/// A palindromic polynomial is gamma-positive if all its gamma coefficients
/// (in the basis {t^i (1+t)^{d-2i}}) are non-negative.
///
/// Returns false for non-palindromic polynomials.
pub fn is_gamma_positive(coeffs: &[i64]) -> bool {
    match gamma_coefficients(coeffs) {
        None => false,
        Some(gammas) => gammas.iter().all(|&g| g >= 0),
    }
}

/// Check if a `BigInt` coefficient polynomial is gamma-positive.
pub fn is_gamma_positive_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    match gamma_coefficients_bigint_coeffs(coeffs) {
        None => false,
        Some(gammas) => gammas.iter().all(|g| g >= &BigInt::zero()),
    }
}

/// Remove all initial zeros, i.e. divide out the largest power of `t`.
pub fn strip_initial_zeros(coeffs: &[i64]) -> &[i64] {
    let start = match coeffs.iter().position(|&c| c != 0) {
        Some(start) => start,
        None => return &[],
    };
    &coeffs[start..]
}

/// Remove all initial zeros from a `BigInt` coefficient vector.
pub fn strip_initial_zeros_bigint_coeffs(coeffs: &[BigInt]) -> &[BigInt] {
    let start = match coeffs.iter().position(|c| !c.is_zero()) {
        Some(start) => start,
        None => return &[],
    };
    &coeffs[start..]
}

/// Check palindromicity after dividing out all factors of `t`.
pub fn is_palindromic_ignoring_initial_zeros(coeffs: &[i64]) -> bool {
    is_palindromic(strip_initial_zeros(coeffs))
}

/// Check palindromicity after dividing out all factors of `t`.
pub fn is_palindromic_ignoring_initial_zeros_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    is_palindromic_bigint_coeffs(strip_initial_zeros_bigint_coeffs(coeffs))
}

/// Compute gamma-coefficients after dividing out all factors of `t`.
pub fn gamma_coefficients_ignoring_initial_zeros(coeffs: &[i64]) -> Option<Vec<i64>> {
    gamma_coefficients(strip_initial_zeros(coeffs))
}

/// Compute gamma coefficients after dividing out all factors of `t`.
pub fn gamma_coefficients_ignoring_initial_zeros_bigint_coeffs(
    coeffs: &[BigInt],
) -> Option<Vec<BigInt>> {
    gamma_coefficients_bigint_coeffs(strip_initial_zeros_bigint_coeffs(coeffs))
}

/// Check gamma-positivity after dividing out all factors of `t`.
pub fn is_gamma_positive_ignoring_initial_zeros(coeffs: &[i64]) -> bool {
    is_gamma_positive(strip_initial_zeros(coeffs))
}

/// Check gamma-positivity after dividing out all factors of `t`.
pub fn is_gamma_positive_ignoring_initial_zeros_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    is_gamma_positive_bigint_coeffs(strip_initial_zeros_bigint_coeffs(coeffs))
}

/// Compute the Stapledon decomposition with respect to a degree bound `n`.
///
/// For a polynomial `p(x)` of degree at most `n`, there is a unique decomposition
///
/// ```text
/// p(x) = a(x) + x b(x),
/// ```
///
/// where `a(x)` is symmetric with center `n/2` and `b(x)` is symmetric with center
/// `(n-1)/2`. This returns the coefficient vectors of `(a(x), b(x))` in ascending
/// degree order, or `None` if `deg(p) > n`.
pub fn stapledon_decomposition(coeffs: &[i64], n: usize) -> Option<(Vec<i64>, Vec<i64>)> {
    let coeffs = coeffs.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
    stapledon_decomposition_bigint_coeffs(&coeffs, n).map(|(a, b)| {
        (
            a.iter()
                .map(|c| i64::try_from(c).expect("Stapledon coefficient too large for i64"))
                .collect(),
            b.iter()
                .map(|c| i64::try_from(c).expect("Stapledon coefficient too large for i64"))
                .collect(),
        )
    })
}

/// Compute the Stapledon decomposition for a `BigInt` coefficient polynomial.
pub fn stapledon_decomposition_bigint_coeffs(
    coeffs: &[BigInt],
    n: usize,
) -> Option<(Vec<BigInt>, Vec<BigInt>)> {
    crate::polynomial::Polynomial::<BigInt>::new(coeffs.to_vec())
        .stapledon_decomposition(n)
        .map(|(a, b)| (a.coeffs().to_vec(), b.coeffs().to_vec()))
}

/// Return the Hermite--Biehler decomposition `(E, O)` where
///
/// ```text
/// p(t) = E(t^2) + t O(t^2).
/// ```
pub fn hermite_biehler_parts(coeffs: &[i64]) -> (Vec<i64>, Vec<i64>) {
    let p = Polynomial::<i64>::from_i64_coeffs(coeffs);
    let (even, odd) = p.hermite_biehler_decomposition();
    (even.coeffs().to_vec(), odd.coeffs().to_vec())
}

/// Check whether a polynomial has only simple roots.
///
/// Returns `false` for the zero polynomial and `true` for nonzero constants.
pub fn has_simple_roots(coeffs: &[i64]) -> bool {
    let p = Polynomial::<Q>::new(
        coeffs
            .iter()
            .map(|&c| Q::from_integer(BigInt::from(c)))
            .collect(),
    );
    p.has_simple_roots()
}

/// Check whether a BigInt-coefficient polynomial has only simple roots.
///
/// Returns `false` for the zero polynomial and `true` for nonzero constants.
pub fn has_simple_roots_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let p = Polynomial::<Q>::new(coeffs.iter().map(|c| Q::from_integer(c.clone())).collect());
    p.has_simple_roots()
}

/// Find all real roots of a polynomial as rational interval midpoints.
///
/// Returns `None` if the polynomial is not real-rooted.
/// Returns `Some(vec![])` for constant/zero polynomials.
pub fn real_roots(coeffs: &[i64]) -> Option<Vec<Q>> {
    let Some(start) = coeffs.iter().position(|&c| c != 0) else {
        return Some(vec![]);
    };
    let end = coeffs
        .iter()
        .rposition(|&c| c != 0)
        .expect("nonzero polynomial has a final nonzero coefficient");
    let trimmed = &coeffs[start..=end];
    let mut roots = if start > 0 {
        vec![Q::zero()]
    } else {
        Vec::new()
    };
    if trimmed.len() <= 1 {
        return Some(roots);
    }
    if trimmed.len() == 2 {
        let a = Q::from_integer(BigInt::from(trimmed[0]));
        let b = Q::from_integer(BigInt::from(trimmed[1]));
        roots.push(-a / b);
        roots.sort();
        return Some(roots);
    }

    let sc = SturmChain::from_i64_coeffs(trimmed);
    let sf_degree = sc.square_free_degree();
    let eps = Q::new(BigInt::from(1), BigInt::from(1_000_000));
    let intervals = sc.isolate_roots(&eps);

    if intervals.len() == sf_degree {
        roots.extend(
            intervals
                .into_iter()
                .map(|(lo, hi)| (&lo + &hi) / Q::from_integer(BigInt::from(2))),
        );
        roots.sort();
        Some(roots)
    } else {
        None
    }
}

/// Check if roots of p interlace roots of q, using Sturm chains.
///
/// This is the slower Sturm-chain method that isolates all roots. Prefer
/// [`check_interlacing`] (Bézout-based) for better performance.
///
/// Convention: `check_interlacing_sturm(f, g)` checks f ≪ g (f has smaller degree).
///
/// Returns `None` if either polynomial is not real-rooted.
/// Returns `Some(false)` if deg(g) ≠ deg(f) + 1.
pub fn check_interlacing_sturm(f: &[i64], g: &[i64]) -> Option<bool> {
    if !is_real_rooted_sturm(f) || !is_real_rooted_sturm(g) {
        return None;
    }

    let (Some(df), Some(dg)) = (poly_degree_trimmed(f), poly_degree_trimmed(g)) else {
        return Some(false);
    };
    if dg != df + 1 {
        return Some(false);
    }
    if !has_simple_roots(f) || !has_simple_roots(g) || resultant(f, g).is_zero() {
        return Some(false);
    }

    let f_chain = SturmChain::from_i64_coeffs(f);
    let g_chain = SturmChain::from_i64_coeffs(g);
    let mut epsilon = Q::one();
    loop {
        let f_intervals = f_chain.isolate_roots(&epsilon);
        let g_intervals = g_chain.isolate_roots(&epsilon);
        debug_assert_eq!(f_intervals.len(), df);
        debug_assert_eq!(g_intervals.len(), dg);

        let certified = f_intervals
            .iter()
            .enumerate()
            .all(|(i, (f_lo, f_hi))| g_intervals[i].1 < *f_lo && *f_hi < g_intervals[i + 1].0);
        if certified {
            return Some(true);
        }

        let definitely_wrong = f_intervals
            .iter()
            .enumerate()
            .any(|(i, (f_lo, f_hi))| *f_hi < g_intervals[i].0 || g_intervals[i + 1].1 < *f_lo);
        if definitely_wrong {
            return Some(false);
        }
        epsilon /= Q::from_integer(BigInt::from(2));
    }
}

// ---------------------------------------------------------------------------
// Bézout matrix approach
// ---------------------------------------------------------------------------

/// Compute the Bézout matrix B(f, g) where deg(f) = deg(g) + 1.
///
/// This is the degree-difference-one Bézoutian used in the interlacing
/// criterion of Kummer--Naldi--Plaumann, Thm. 2.13, after Krein--Naimark,
/// §2.2.  Same-degree interlacing is handled by [`check_weak_interlacing`]
/// via a degree-extension reduction; see Fisk, §9.21, Cor. 9.145 for the
/// same-degree positive definite form.
///
/// The (i, j) entry (0-indexed) is the coefficient of x^i y^j in
///
/// ```text
/// (f(x)g(y) - f(y)g(x)) / (x - y)
/// ```
///
/// Returns a d x d matrix (as `i64`) where d = deg(f), or `None` if the degree
/// constraint is not satisfied or an exact entry does not fit in `i64`.
///
/// The division is exact because h(x,y) = f(x)g(y) - f(y)g(x) vanishes on x = y.
/// Expanding: for each pair (k, l) with k > l, the contribution is
///
/// ```text
/// (f_k g_l - f_l g_k) * sum_{m=0}^{k-l-1} x^{l+m} y^{k-1-m}
/// ```
///
/// which follows from factoring x^k y^l - x^l y^k = x^l y^l (x^{k-l} - y^{k-l})
/// and using the geometric sum (x^n - y^n)/(x - y) = sum x^i y^{n-1-i}.
pub fn bezout_matrix(f: &[i64], g: &[i64]) -> Option<Vec<Vec<i64>>> {
    bezout_matrix_bigint(f, g)?
        .into_iter()
        .map(|row| row.into_iter().map(|entry| entry.to_i64()).collect())
        .collect()
}

/// Compute the Bézout matrix with BigInt entries (no overflow).
fn bezout_matrix_bigint(f: &[i64], g: &[i64]) -> Option<Vec<Vec<BigInt>>> {
    let f_big: Vec<BigInt> = f.iter().map(|&c| BigInt::from(c)).collect();
    let g_big: Vec<BigInt> = g.iter().map(|&c| BigInt::from(c)).collect();
    bezout_matrix_bigint_coeffs(&f_big, &g_big)
}

/// Compute the Bézout matrix from BigInt coefficient vectors.
///
/// The coefficient vectors are in ascending degree order.  This is the exact
/// no-overflow variant of [`bezout_matrix`], useful for benchmarking and for
/// downstream linear-algebra experiments.
pub fn bezout_matrix_bigint_coeffs(f: &[BigInt], g: &[BigInt]) -> Option<Vec<Vec<BigInt>>> {
    let zero_big = BigInt::from(0);
    let df = {
        let mut d = None;
        for i in (0..f.len()).rev() {
            if f[i] != zero_big {
                d = Some(i);
                break;
            }
        }
        d
    }?;
    let dg = {
        let mut d = None;
        for i in (0..g.len()).rev() {
            if g[i] != zero_big {
                d = Some(i);
                break;
            }
        }
        d
    }?;
    if df != dg + 1 {
        return None;
    }
    let d = df;
    let zero = BigInt::from(0);
    let mut b = vec![vec![zero.clone(); d]; d];

    for k in 0..=df {
        for l in 0..k {
            let fk = if k < f.len() { &f[k] } else { &zero };
            let fl = if l < f.len() { &f[l] } else { &zero };
            let gk = if k < g.len() { &g[k] } else { &zero };
            let gl = if l < g.len() { &g[l] } else { &zero };
            let c = fk * gl - fl * gk;
            if c == zero {
                continue;
            }
            for m in 0..=(k - l - 1) {
                let xi = l + m;
                let yj = k - 1 - m;
                if xi < d && yj < d {
                    b[xi][yj] = &b[xi][yj] + &c;
                }
            }
        }
    }
    Some(b)
}

/// Check positive definiteness of a symmetric BigInt matrix via Gaussian elimination.
fn is_positive_definite_bigint(mat: &[Vec<BigInt>]) -> bool {
    crate::linalg::is_positive_definite(mat)
}

fn is_positive_semidefinite_bigint(mat: &[Vec<BigInt>]) -> bool {
    crate::linalg::is_positive_semidefinite(mat)
}

/// Check **strict** interlacing f ≪ g via the Bézout matrix.
///
/// Convention: `check_interlacing(f, g)` checks f ≪ g, where f has degree d−1
/// and g has degree d. The caller must pass f (smaller degree) first.
///
/// Two real-rooted polynomials with deg(g) = deg(f) + 1 strictly interlace
/// if and only if the Bézout matrix B(g, f) is positive definite. This is the
/// coprime/strict form of the degree-difference-one criterion in
/// Kummer--Naldi--Plaumann, Thm. 2.13, citing Krein--Naimark, §2.2.
///
/// This is 100–400× faster than Sturm chains at degree 15+ because it avoids
/// root isolation entirely — just one exact Gaussian elimination on a d×d matrix.
///
/// Returns `Some(true)` if f ≪ g, `Some(false)` if not,
/// or `None` if deg(g) ≠ deg(f) + 1.
///
/// # Limitations
///
/// - Requires deg(g) = deg(f) + 1 exactly; returns `None` otherwise.
/// - The Bézout matrix is singular when f and g share a root, so this function
///   returns `Some(false)` in that case.
///
/// For same-degree interlacing, shared roots, or either argument order, use
/// [`check_weak_interlacing`].
pub fn check_interlacing(f: &[i64], g: &[i64]) -> Option<bool> {
    let f_big: Vec<BigInt> = f.iter().map(|&c| BigInt::from(c)).collect();
    let g_big: Vec<BigInt> = g.iter().map(|&c| BigInt::from(c)).collect();
    check_interlacing_bigint_coeffs(&f_big, &g_big)
}

/// Check **strict** interlacing f ≪ g for `BigInt` coefficient vectors.
///
/// This has the same convention as [`check_interlacing`], but avoids all
/// coefficient overflow. The leading coefficients are normalized to be
/// positive before the Bézout matrix is formed, so the result is invariant
/// under multiplying either input polynomial by a nonzero scalar.
pub fn check_interlacing_bigint_coeffs(f: &[BigInt], g: &[BigInt]) -> Option<bool> {
    let df = poly_degree_trimmed_bigint(f);
    let dg = poly_degree_trimmed_bigint(g);

    // 0 ≪ g for any real-rooted g (vacuously true).
    if df.is_none() {
        return Some(is_real_rooted_bigint_coeffs(g));
    }

    let (df, dg) = match (df, dg) {
        (Some(df), Some(dg)) if dg == df + 1 => (df, dg),
        _ => return None,
    };

    let f = normalize_positive_leading_bigint(f, df);
    let g = normalize_positive_leading_bigint(g, dg);
    let mat = bezout_matrix_bigint_coeffs(&g, &f)?;
    Some(is_positive_definite_bigint(&mat))
}

fn normalize_positive_leading_bigint(coeffs: &[BigInt], degree: usize) -> Vec<BigInt> {
    let leading = coeffs[degree].clone();
    debug_assert!(!leading.is_zero());
    let mut normalized: Vec<BigInt> = coeffs.iter().take(degree + 1).cloned().collect();
    if leading.is_negative() {
        for coeff in &mut normalized {
            *coeff = -coeff.clone();
        }
    }
    normalized
}

/// Check **directed weak** interlacing `p ≪ q` via the Bézout matrix.
///
/// Tests whether `p` interlaces `q` **from the left**: the roots of `p` and `q`
/// alternate with p's roots weakly to the left of q's roots:
///
/// ```text
/// ... ≤ a₃ ≤ b₃ ≤ a₂ ≤ b₂ ≤ a₁ ≤ b₁ ≤ 0
/// ```
///
/// where aᵢ are roots of p and bᵢ are roots of q. Shared roots are allowed.
/// Requires `deg(p) = deg(q)` or `deg(p) = deg(q) - 1` (equivalently,
/// `deg(q) ∈ {deg(p), deg(p) + 1}`). The zero polynomial interlaces
/// any real-rooted polynomial from the left (vacuously).
///
/// **This function is NOT symmetric**: `check_weak_interlacing(p, q)` tests
/// `p ≪ q`, which is different from `q ≪ p` when `deg(p) = deg(q)`.
///
/// Handles:
/// - **deg(q) = deg(p) + 1**: divides out gcd(p,q) if nontrivial (verifying
///   the shared roots are real), then checks strict interlacing of the reduced pair.
/// - **deg(p) = deg(q)** (same degree): extends p with a root far to the right
///   (via the Cauchy bound), giving p_ext of degree deg(p)+1, then checks q ≪ p_ext.
///   **Fast path**: if all coefficients of p are positive (all roots negative),
///   multiplying by t (root at 0) suffices.
///
/// References: Kummer--Naldi--Plaumann, Thm. 2.13 for the
/// degree-difference-one Bézoutian criterion; Fisk, §9.21, Cor. 9.145 for the
/// same-degree oriented alternation criterion. The explicit gcd real-rootedness
/// check below guards the semidefinite/common-factor caveat.
///
/// Returns `Some(true)` if `p ≪ q`, `Some(false)` if not, or `None`
/// if the degree relationship is invalid (deg(p) > deg(q) or deg(q) > deg(p) + 1).
pub fn check_weak_interlacing(p: &[i64], q: &[i64]) -> Option<bool> {
    let p_big: Vec<BigInt> = p.iter().map(|&c| BigInt::from(c)).collect();
    let q_big: Vec<BigInt> = q.iter().map(|&c| BigInt::from(c)).collect();
    check_weak_interlacing_bigint_coeffs(&p_big, &q_big)
}

/// Check **directed weak** interlacing `p ≪ q` for `BigInt` coefficients.
///
/// This has the same convention as [`check_weak_interlacing`], but avoids
/// overflow in the same-degree Cauchy-bound reduction and in the reduced
/// Bézout matrix after common factors are removed.
pub fn check_weak_interlacing_bigint_coeffs(p: &[BigInt], q: &[BigInt]) -> Option<bool> {
    // 0 ≪ q for any real-rooted q (vacuously true).
    if poly_degree_trimmed_bigint(p).is_none() {
        return Some(is_real_rooted_bigint_coeffs(q));
    }
    // p ≪ 0 only if p is also zero (already handled above).
    if poly_degree_trimmed_bigint(q).is_none() {
        return Some(false);
    }

    let dp = poly_degree_trimmed_bigint(p).unwrap();
    let dq = poly_degree_trimmed_bigint(q).unwrap();

    // Same-degree case: p ≪ q means p's roots are to the LEFT of q's roots.
    // The alternation pattern is: p_d ≤ q_d ≤ ... ≤ p_1 ≤ q_1 ≤ 0.
    //
    // Reduction to the deg+1 case: extend p (the LEFT polynomial) with a root
    // far to the RIGHT. Then p_ext has degree d+1 with roots {p_d,...,p_1,R},
    // and we check q ≪ p_ext, which gives:
    //   p_d ≤ q_d ≤ p_{d-1} ≤ ... ≤ q_1 ≤ R.
    // Since R is far right, q_1 ≤ R holds automatically, and the rest
    // recovers the original same-degree interlacing p_i ≤ q_i.
    if dp == dq {
        let p_ext = if has_nonzero_one_signed_trimmed_coefficients_bigint(p)
            && has_nonzero_one_signed_trimmed_coefficients_bigint(q)
        {
            // Both polynomials have no nonnegative real roots, so root at 0 is
            // far right enough if the pair is real-rooted.
            poly_mul_linear_factor_bigint(p, &BigInt::zero())
        } else {
            let r = cauchy_root_bound_bigint(p).max(cauchy_root_bound_bigint(q)) + BigInt::from(1);
            poly_mul_linear_factor_bigint(p, &(-r))
        };
        return check_weak_interlacing_impl(q, &p_ext);
    }

    // deg(q) = deg(p) + 1: standard case, p has smaller degree.
    if dq == dp + 1 {
        return check_weak_interlacing_impl(p, q);
    }

    // deg(p) > deg(q): invalid for p ≪ q.
    None
}

/// Cauchy bound for root radius: all roots of `coeffs` lie in |z| < bound.
///
/// For f(t) = a_n t^n + ... + a_0, the Cauchy bound is
/// 1 + max(|a_0|, ..., |a_{n-1}|) / |a_n|.
/// Returns 2 for constant or zero polynomials.
fn cauchy_root_bound_bigint(coeffs: &[BigInt]) -> BigInt {
    let deg = match poly_degree_trimmed_bigint(coeffs) {
        Some(d) if d > 0 => d,
        _ => return BigInt::from(2),
    };
    let lc = coeffs[deg].abs();
    let max_other = coeffs[..deg]
        .iter()
        .map(|c| c.abs())
        .max()
        .unwrap_or_else(BigInt::zero);
    BigInt::from(1) + (&max_other + &lc - BigInt::from(1)) / lc
}

/// Multiply polynomial by (t + r) in BigInt coefficients.
fn poly_mul_linear_factor_bigint(coeffs: &[BigInt], r: &BigInt) -> Vec<BigInt> {
    // (t + r) * (a_0 + a_1 t + ... + a_n t^n) = r*a_0 + (a_0 + r*a_1)t + ... + a_n t^{n+1}
    let n = coeffs.len();
    let mut result = vec![BigInt::zero(); n + 1];
    for i in 0..n {
        result[i] = &result[i] + r * &coeffs[i];
        result[i + 1] = &result[i + 1] + &coeffs[i];
    }
    result
}

fn has_nonzero_one_signed_trimmed_coefficients_bigint(coeffs: &[BigInt]) -> bool {
    let Some(degree) = poly_degree_trimmed_bigint(coeffs) else {
        return false;
    };
    if coeffs.first().is_none_or(BigInt::is_zero) {
        return false;
    }
    let active = &coeffs[..=degree];
    active.iter().all(|c| c >= &BigInt::zero()) || active.iter().all(|c| c <= &BigInt::zero())
}

/// Core weak interlacing check for deg(g) = deg(f) + 1.
fn check_weak_interlacing_impl(f: &[BigInt], g: &[BigInt]) -> Option<bool> {
    let to_q = |coeffs: &[BigInt]| -> Vec<Q> {
        coeffs.iter().map(|c| Q::from_integer(c.clone())).collect()
    };

    let pq = to_q(f);
    let qq = to_q(g);

    let gcd = poly_gcd_q(&pq, &qq);
    let gcd_deg = q_degree(&gcd);

    // If GCD is trivial (constant), just do strict interlacing
    if gcd_deg == 0 {
        return check_interlacing_bigint_coeffs(f, g);
    }

    // The shared roots (given by the GCD) must themselves be real. PSD of the
    // Bézoutian does not certify this: e.g. f=t(t^2+1), g=t^2+1 has a non-real
    // common factor, so the original pair cannot weakly interlace even though
    // the reduced pair is harmless.
    let gcd_big = q_poly_to_primitive_bigint(&gcd);
    if !is_real_rooted_bigint_coeffs(&gcd_big) {
        return Some(false);
    }

    // Divide out the GCD
    let f_red = poly_exact_div_q(&pq, &gcd);
    let g_red = poly_exact_div_q(&qq, &gcd);

    let df = q_degree(&f_red);
    let dg = q_degree(&g_red);

    if df == 0 && dg == 0 {
        // Both reduced to constants: all roots are shared.
        return Some(true);
    }

    // After dividing out GCD, f_red should have smaller degree than g_red
    if dg != df + 1 {
        return None;
    }

    let f_big = q_poly_to_primitive_bigint(&f_red);
    let g_big = q_poly_to_primitive_bigint(&g_red);
    check_interlacing_bigint_coeffs(&f_big, &g_big)
}

// ---------------------------------------------------------------------------
// Polynomial arithmetic over Q (for GCD computation)
// ---------------------------------------------------------------------------

fn q_degree(p: &[Q]) -> usize {
    let zero = Q::from_integer(BigInt::from(0));
    p.iter().rposition(|c| *c != zero).unwrap_or(0)
}

fn poly_gcd_q(a: &[Q], b: &[Q]) -> Vec<Q> {
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

fn poly_exact_div_q(a: &[Q], b: &[Q]) -> Vec<Q> {
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
fn q_poly_to_primitive_bigint(p: &[Q]) -> Vec<BigInt> {
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

fn derivative_bigint_coeffs(coeffs: &[BigInt]) -> Vec<BigInt> {
    let Some(degree) = poly_degree_trimmed_bigint(coeffs) else {
        return vec![BigInt::zero()];
    };
    if degree == 0 {
        return vec![BigInt::zero()];
    }
    let mut derivative = Vec::with_capacity(degree);
    for k in 0..degree {
        derivative.push(&coeffs[k + 1] * BigInt::from((k + 1) as u64));
    }
    derivative
}

/// Return the squarefree part `p / gcd(p, p')`, with primitive integer
/// coefficients.
///
/// This is useful for boolean real-rootedness checks: a polynomial is weakly
/// real-rooted if and only if its squarefree part is strictly real-rooted.
/// The zero polynomial is returned as `0`.
pub fn squarefree_part_bigint_coeffs(coeffs: &[BigInt]) -> Vec<BigInt> {
    let Some(degree) = poly_degree_trimmed_bigint(coeffs) else {
        return vec![BigInt::zero()];
    };
    let trimmed: Vec<BigInt> = coeffs.iter().take(degree + 1).cloned().collect();
    if degree <= 1 {
        return q_poly_to_primitive_bigint(
            &trimmed
                .iter()
                .map(|c| Q::from_integer(c.clone()))
                .collect::<Vec<_>>(),
        );
    }

    let p_q: Vec<Q> = trimmed.iter().map(|c| Q::from_integer(c.clone())).collect();
    let dp_q: Vec<Q> = derivative_bigint_coeffs(&trimmed)
        .iter()
        .map(|c| Q::from_integer(c.clone()))
        .collect();
    let gcd = poly_gcd_q(&p_q, &dp_q);
    let squarefree = poly_exact_div_q(&p_q, &gcd);
    q_poly_to_primitive_bigint(&squarefree)
}

fn newton_sums_rational_bigint_coeffs(coeffs: &[BigInt], max_power: usize) -> Option<Vec<Q>> {
    let degree = poly_degree_trimmed_bigint(coeffs)?;
    let lc = coeffs[degree].clone();
    if lc.is_zero() {
        return None;
    }

    let lc_q = Q::from_integer(lc);
    let monic_coeffs: Vec<Q> = coeffs
        .iter()
        .take(degree)
        .map(|c| Q::from_integer(c.clone()) / lc_q.clone())
        .collect();

    let mut sums = vec![Q::from_integer(BigInt::zero()); max_power + 1];
    sums[0] = Q::from_integer(BigInt::from(degree as u64));

    for k in 1..=max_power {
        let mut total = Q::from_integer(BigInt::zero());
        let upper = (k - 1).min(degree);
        for i in 1..=upper {
            total += monic_coeffs[degree - i].clone() * sums[k - i].clone();
        }
        if k <= degree {
            total += Q::from_integer(BigInt::from(k as u64)) * monic_coeffs[degree - k].clone();
        }
        sums[k] = -total;
    }

    Some(sums)
}

/// Build the Hermite/Newton-sum matrix for a polynomial.
///
/// For a degree `d` polynomial with roots `alpha_i`, the rational Hermite
/// matrix has entry `s_{i+j}`, where `s_k=sum_i alpha_i^k`.  This function
/// clears denominators by a single positive scalar and returns an integer
/// matrix with the same definiteness.
pub fn hermite_matrix_bigint_coeffs(coeffs: &[BigInt]) -> Option<Vec<Vec<BigInt>>> {
    let degree = poly_degree_trimmed_bigint(coeffs)?;
    if degree == 0 {
        return Some(Vec::new());
    }

    let sums = newton_sums_rational_bigint_coeffs(coeffs, 2 * degree - 2)?;
    let zero = Q::from_integer(BigInt::zero());
    let denominator_lcm = sums
        .iter()
        .take(2 * degree - 1)
        .filter(|s| **s != zero)
        .fold(BigInt::from(1), |acc, s| acc.lcm(s.denom()));
    let scale = Q::from_integer(denominator_lcm);

    let mut matrix = vec![vec![BigInt::zero(); degree]; degree];
    for i in 0..degree {
        for j in 0..degree {
            let scaled = sums[i + j].clone() * scale.clone();
            debug_assert_eq!(scaled.denom(), &BigInt::from(1));
            matrix[i][j] = scaled.to_integer();
        }
    }
    Some(matrix)
}

/// Check real-rootedness by first taking the squarefree part, then checking
/// strict real-rootedness with the Bézout positive-definite criterion.
///
/// This is kept separate from [`is_real_rooted_bezout_bigint_coeffs`] so the
/// old semidefinite Bézout implementation remains available for comparison.
pub fn is_real_rooted_bezout_squarefree_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let squarefree = squarefree_part_bigint_coeffs(coeffs);
    is_strictly_real_rooted_bezout_bigint_coeffs(&squarefree)
}

/// Check strict real-rootedness with the Bézout positive-definite criterion,
/// without removing repeated factors first.
pub fn is_strictly_real_rooted_bezout_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let Some(degree) = poly_degree_trimmed_bigint(coeffs) else {
        return true;
    };
    if degree <= 1 {
        return true;
    }
    let derivative = derivative_bigint_coeffs(coeffs);
    check_interlacing_bigint_coeffs(&derivative, coeffs).unwrap_or(false)
}

/// Check real-rootedness by first taking the squarefree part, then checking
/// positive definiteness of the Hermite/Newton-sum matrix.
pub fn is_real_rooted_hermite_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let squarefree = squarefree_part_bigint_coeffs(coeffs);
    is_strictly_real_rooted_hermite_bigint_coeffs(&squarefree)
}

/// Check strict real-rootedness with the Hermite/Newton-sum positive-definite
/// criterion, without removing repeated factors first.
pub fn is_strictly_real_rooted_hermite_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let Some(degree) = poly_degree_trimmed_bigint(coeffs) else {
        return true;
    };
    if degree <= 1 {
        return true;
    }
    let Some(matrix) = hermite_matrix_bigint_coeffs(coeffs) else {
        return false;
    };
    is_positive_definite_bigint(&matrix)
}

/// Explicit squarefree+Bézout real-rootedness check for `i64` coefficients.
pub fn is_real_rooted_bezout_squarefree(coeffs: &[i64]) -> bool {
    let big: Vec<BigInt> = coeffs.iter().map(|&c| BigInt::from(c)).collect();
    is_real_rooted_bezout_squarefree_bigint_coeffs(&big)
}

/// Explicit squarefree+Hermite/Newton real-rootedness check for `i64`
/// coefficients.
pub fn is_real_rooted_hermite(coeffs: &[i64]) -> bool {
    let big: Vec<BigInt> = coeffs.iter().map(|&c| BigInt::from(c)).collect();
    is_real_rooted_hermite_bigint_coeffs(&big)
}

/// Explicit strict Bézout real-rootedness check for `i64` coefficients.
pub fn is_strictly_real_rooted_bezout(coeffs: &[i64]) -> bool {
    let big: Vec<BigInt> = coeffs.iter().map(|&c| BigInt::from(c)).collect();
    is_strictly_real_rooted_bezout_bigint_coeffs(&big)
}

/// Explicit strict Hermite/Newton real-rootedness check for `i64`
/// coefficients.
pub fn is_strictly_real_rooted_hermite(coeffs: &[i64]) -> bool {
    let big: Vec<BigInt> = coeffs.iter().map(|&c| BigInt::from(c)).collect();
    is_strictly_real_rooted_hermite_bigint_coeffs(&big)
}

/// Check if a polynomial is real-rooted.
///
/// This is the public default exact path.  It delegates to the primitive
/// integer PRS root-counting implementation in [`crate::root_count`].  Use
/// [`is_real_rooted_bezout`], [`is_real_rooted_bezout_squarefree`], or
/// [`is_real_rooted_hermite`] when an explicit matrix-based comparison or
/// certificate is desired.
pub fn is_real_rooted(coeffs: &[i64]) -> bool {
    crate::root_count::is_real_rooted_fast_i64(coeffs)
}

fn is_real_rooted_bezout_i64_impl(coeffs: &[i64]) -> bool {
    let d = match poly_degree_trimmed(coeffs) {
        Some(d) => d,
        None => return true,
    };
    if d <= 1 {
        return true;
    }

    // Compute f'(t) in BigInt to avoid i64 overflow.
    // f' = Σ (k+1) a_{k+1} t^k, degree d-1.
    let mut fp: Vec<BigInt> = Vec::with_capacity(d);
    let mut overflowed = false;
    for k in 0..d {
        let mult = (k as i64) + 1;
        let ak1 = coeff_i64(coeffs, k + 1);
        match mult.checked_mul(ak1) {
            Some(v) => fp.push(BigInt::from(v)),
            None => {
                fp.push(BigInt::from(mult) * BigInt::from(ak1));
                overflowed = true;
            }
        }
    }

    // Ensure leading coefficients have the same sign.
    let lc_f = BigInt::from(coeff_i64(coeffs, d));
    let lc_fp = fp.last().cloned().unwrap_or_else(|| BigInt::from(0));
    let zero_big = BigInt::from(0);
    if (lc_f > zero_big) != (lc_fp > zero_big) {
        for c in fp.iter_mut() {
            *c = -c.clone();
        }
    }

    if overflowed {
        // Use fully BigInt Bézout path
        let f_big: Vec<BigInt> = coeffs.iter().map(|&c| BigInt::from(c)).collect();
        let mat = match bezout_matrix_bigint_coeffs(&f_big, &fp) {
            Some(m) => m,
            None => return false,
        };
        is_positive_semidefinite_bigint(&mat)
    } else {
        // Convert back to i64 for the existing path
        let fp_i64: Vec<i64> = fp
            .iter()
            .map(|b| {
                use num_traits::ToPrimitive;
                b.to_i64().unwrap()
            })
            .collect();
        let mat = match bezout_matrix_bigint(coeffs, &fp_i64) {
            Some(m) => m,
            None => return false,
        };
        is_positive_semidefinite_bigint(&mat)
    }
}

/// Check if a polynomial with `BigInt` coefficients is real-rooted.
///
/// This is the public default exact path.  It uses the primitive integer PRS
/// root-counting implementation in [`crate::root_count`], with a specialized
/// one-signed shortcut when applicable.
///
/// Coefficients are given in ascending degree order.
pub fn is_real_rooted_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    crate::root_count::is_real_rooted_fast_bigint_coeffs(coeffs)
}

/// Check if a polynomial with `BigInt` coefficients is real-rooted using the
/// old semidefinite Bézout-matrix criterion.
///
/// This is the old exact BigInt implementation and is kept public for
/// benchmarking and cases where the matrix certificate itself is desired.
pub fn is_real_rooted_bezout_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let d = match poly_degree_trimmed_bigint(coeffs) {
        Some(d) => d,
        None => return true,
    };
    if d <= 1 {
        return true;
    }

    let mut gcd = BigInt::zero();
    for c in coeffs.iter().take(d + 1) {
        if !c.is_zero() {
            gcd = if gcd.is_zero() {
                c.abs()
            } else {
                gcd.gcd(&c.abs())
            };
        }
    }
    if !gcd.is_zero() {
        let reduced: Vec<BigInt> = coeffs.iter().take(d + 1).map(|c| c / &gcd).collect();
        if reduced.iter().all(|c| c.to_i64().is_some()) {
            let reduced_i64: Vec<i64> = reduced
                .iter()
                .map(|c| c.to_i64().expect("checked above"))
                .collect();
            return is_real_rooted_bezout(&reduced_i64);
        }
    }

    let mut fp: Vec<BigInt> = Vec::with_capacity(d);
    for k in 0..d {
        fp.push(&coeffs[k + 1] * BigInt::from((k + 1) as u64));
    }

    let lc_f = coeffs[d].clone();
    let lc_fp = fp.last().cloned().unwrap_or_else(|| BigInt::from(0));
    let zero_big = BigInt::from(0);
    if (lc_f > zero_big) != (lc_fp > zero_big) {
        for c in &mut fp {
            *c = -c.clone();
        }
    }

    let mat = match bezout_matrix_bigint_coeffs(coeffs, &fp) {
        Some(m) => m,
        None => return false,
    };
    is_positive_semidefinite_bigint(&mat)
}

// ---------------------------------------------------------------------------
// Aliases for explicit algorithm selection
// ---------------------------------------------------------------------------

/// Explicit Bézout-matrix real-rootedness check for `i64` coefficients.
pub fn is_real_rooted_bezout(coeffs: &[i64]) -> bool {
    is_real_rooted_bezout_i64_impl(coeffs)
}

/// Alias for [`check_interlacing`] — Bézout matrix method (the default).
pub fn check_interlacing_bezout(p: &[i64], q: &[i64]) -> Option<bool> {
    check_interlacing(p, q)
}

/// Alias for [`check_weak_interlacing`] — Bézout matrix with GCD factoring.
pub fn check_weak_interlacing_bezout(p: &[i64], q: &[i64]) -> Option<bool> {
    check_weak_interlacing(p, q)
}

// ---------------------------------------------------------------------------
// Resultant and discriminant
// ---------------------------------------------------------------------------

/// Compute the Sylvester matrix of f and g.
///
/// For f of degree n and g of degree m, the Sylvester matrix is (n+m) x (n+m).
/// The first m rows are shifted copies of f's coefficients (in descending order),
/// and the last n rows are shifted copies of g's coefficients.
///
/// Returns `None` if either polynomial is zero.
pub fn sylvester_matrix(f: &[i64], g: &[i64]) -> Option<Vec<Vec<BigInt>>> {
    let n = poly_degree_trimmed(f)?;
    let m = poly_degree_trimmed(g)?;
    let size = n + m;
    let zero = BigInt::from(0);
    let mut mat = vec![vec![zero.clone(); size]; size];

    // First m rows: coefficients of f, shifted right by row index.
    // Row i has f_n, f_{n-1}, ..., f_0 starting at column i.
    for i in 0..m {
        for k in 0..=n {
            mat[i][i + k] = BigInt::from(coeff_i64(f, n - k));
        }
    }
    // Last n rows: coefficients of g, shifted right by row index.
    for i in 0..n {
        for k in 0..=m {
            mat[m + i][i + k] = BigInt::from(coeff_i64(g, m - k));
        }
    }
    Some(mat)
}

fn sylvester_matrix_bigint_coeffs(f: &[BigInt], g: &[BigInt]) -> Option<Vec<Vec<BigInt>>> {
    let n = poly_degree_trimmed_bigint(f)?;
    let m = poly_degree_trimmed_bigint(g)?;
    let size = n + m;
    let zero = BigInt::from(0);
    let mut mat = vec![vec![zero; size]; size];

    for i in 0..m {
        for k in 0..=n {
            mat[i][i + k] = f[n - k].clone();
        }
    }
    for i in 0..n {
        for k in 0..=m {
            mat[m + i][i + k] = g[m - k].clone();
        }
    }
    Some(mat)
}

fn determinant_bigint(mat: &[Vec<BigInt>]) -> BigInt {
    crate::linalg::determinant(mat)
}

/// Compute the resultant Res(f, g) of two polynomials.
///
/// The resultant is zero if and only if f and g share a common root (over the
/// algebraic closure) or one of them is zero. Computed via the Sylvester matrix
/// determinant.
///
/// Returns the resultant as a `BigInt`.
pub fn resultant(f: &[i64], g: &[i64]) -> BigInt {
    let f_big = f.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
    let g_big = g.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
    resultant_bigint_coeffs(&f_big, &g_big)
}

/// Compute the resultant for `BigInt` coefficient polynomials.
pub fn resultant_bigint_coeffs(f: &[BigInt], g: &[BigInt]) -> BigInt {
    match sylvester_matrix_bigint_coeffs(f, g) {
        Some(mat) => determinant_bigint(&mat),
        None => BigInt::from(0),
    }
}

/// Compute the discriminant of a polynomial.
///
/// ```text
/// disc(f) = (-1)^{n(n-1)/2} * Res(f, f') / lc(f)
/// ```
///
/// where n = deg(f) and f' is the formal derivative.
///
/// The discriminant is zero iff f has a repeated root. For a polynomial with
/// roots r_1, ..., r_n, disc(f) = lc(f)^{2n-2} * prod_{i<j} (r_i - r_j)^2
/// (up to the sign convention).
pub fn discriminant(f: &[i64]) -> BigInt {
    let f_big = f.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
    discriminant_bigint_coeffs(&f_big)
}

/// Compute the discriminant for a `BigInt` coefficient polynomial.
pub fn discriminant_bigint_coeffs(f: &[BigInt]) -> BigInt {
    let n = match poly_degree_trimmed_bigint(f) {
        Some(n) => n,
        None => return BigInt::from(0),
    };
    if n <= 1 {
        return BigInt::from(1);
    }

    let mut fp: Vec<BigInt> = Vec::with_capacity(n);
    for k in 0..n {
        fp.push(BigInt::from(k + 1) * &f[k + 1]);
    }

    let res = match sylvester_matrix_bigint_coeffs(f, &fp) {
        Some(mat) => determinant_bigint(&mat),
        None => BigInt::from(0),
    };
    let lc = f[n].clone();
    let sign_exp = n * (n - 1) / 2;
    let sign = if sign_exp % 2 == 0 {
        BigInt::from(1)
    } else {
        BigInt::from(-1)
    };

    sign * res / lc
}

// ---------------------------------------------------------------------------
// Ehrhart polynomial <-> h*-vector conversion
// ---------------------------------------------------------------------------

/// Convert an h\*-vector to an Ehrhart polynomial.
///
/// Given h\* = (h\*\_0, ..., h\*\_d), the Ehrhart polynomial is
///
/// ```text
/// L(t) = sum_{i=0}^{d} h*_i * C(t + d - i, d)
/// ```
///
/// where C(t+a, d) = (t+a)(t+a-1)...(t+a-d+1)/d! is a polynomial of degree d in t.
///
/// Returns the polynomial coefficients as rationals in ascending degree order.
/// For a valid h\*-vector of a lattice polytope, L(t) has degree d = len(hstar) - 1.
pub fn hstar_to_ehrhart(hstar: &[i64]) -> Vec<Q> {
    let hstar = hstar.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
    hstar_to_ehrhart_bigint_coeffs(&hstar)
}

/// Convert an arbitrary-size integer h\*-vector to an Ehrhart polynomial.
pub fn hstar_to_ehrhart_bigint_coeffs(hstar: &[BigInt]) -> Vec<Q> {
    if hstar.is_empty() {
        return vec![];
    }
    let d = hstar.len() - 1;
    if d == 0 {
        return vec![Q::from_integer(hstar[0].clone())];
    }

    let zero = Q::from_integer(BigInt::from(0));
    let one = Q::from_integer(BigInt::from(1));
    let mut result = vec![zero.clone(); d + 1];

    // Compute d! for normalization
    let mut d_fact = Q::from_integer(BigInt::from(1));
    for j in 2..=d {
        d_fact *= Q::from_integer(BigInt::from(j as i64));
    }

    for (i, hi) in hstar.iter().enumerate() {
        if hi == &BigInt::from(0) {
            continue;
        }
        // Build C(t + d - i, d) as a polynomial in t.
        // C(t + a, d) = product_{j=0}^{d-1} (t + a - j) / d!  where a = d - i.
        let a = (d - i) as i64;

        // Multiply linear factors (t + a)(t + a - 1)...(t + a - d + 1)
        // Start with polynomial = 1
        let mut poly = vec![one.clone()]; // degree 0
        for j in 0..d {
            let shift = a - j as i64; // constant term of the linear factor (t + shift)
            let shift_q = Q::from_integer(BigInt::from(shift));
            // Multiply poly by (t + shift): new[k] = poly[k-1] + shift * poly[k]
            let mut new_poly = vec![zero.clone(); poly.len() + 1];
            for (k, c) in poly.iter().enumerate() {
                new_poly[k] = new_poly[k].clone() + shift_q.clone() * c.clone();
                new_poly[k + 1] = new_poly[k + 1].clone() + c.clone();
            }
            poly = new_poly;
        }

        // Divide by d! and scale by h*_i
        let scale = Q::from_integer(hi.clone()) / d_fact.clone();
        for (k, c) in poly.iter().enumerate() {
            if k <= d {
                result[k] = result[k].clone() + scale.clone() * c.clone();
            }
        }
    }

    // Trim trailing zeros
    while result.last().is_some_and(|c| *c == zero) {
        result.pop();
    }
    result
}

/// Convert an Ehrhart polynomial to an h\*-vector.
///
/// Given L(t) of degree d with rational coefficients, the h\*-vector entries are
///
/// ```text
/// h*_i = sum_{k=0}^{i} (-1)^k * C(d+1, k) * L(i-k)    for i = 0, ..., d
/// ```
///
/// where L(j) is the Ehrhart polynomial evaluated at the non-negative integer j.
///
/// The h\*-vector entries are always integers for Ehrhart polynomials of lattice
/// polytopes. Nonintegral results and entries outside `i64` are reported as
/// errors rather than rounded or panicked on.
pub fn ehrhart_to_hstar(ehrhart_coeffs: &[Q]) -> Result<Vec<i64>, EhrhartToHstarError> {
    ehrhart_to_hstar_bigint(ehrhart_coeffs)?
        .into_iter()
        .enumerate()
        .map(|(index, value)| {
            i64::try_from(&value)
                .map_err(|_| EhrhartToHstarError::EntryOutOfI64Range { index, value })
        })
        .collect()
}

/// Errors returned when exact Ehrhart-to-h\* conversion is not defined for the
/// supplied representation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EhrhartToHstarError {
    ZeroDenominator,
    NonIntegralEntry { index: usize, value: Q },
    EntryOutOfI64Range { index: usize, value: BigInt },
}

impl std::fmt::Display for EhrhartToHstarError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ZeroDenominator => write!(f, "Ehrhart coefficient denominator must be nonzero"),
            Self::NonIntegralEntry { index, value } => {
                write!(f, "h*-vector entry {index} is not an integer: {value}")
            }
            Self::EntryOutOfI64Range { index, value } => {
                write!(f, "h*-vector entry {index} does not fit in i64: {value}")
            }
        }
    }
}

impl std::error::Error for EhrhartToHstarError {}

/// Convert an Ehrhart polynomial to an arbitrary-size integer h\*-vector.
pub fn ehrhart_to_hstar_bigint(ehrhart_coeffs: &[Q]) -> Result<Vec<BigInt>, EhrhartToHstarError> {
    if ehrhart_coeffs.is_empty() {
        return Ok(vec![]);
    }
    let zero = Q::from_integer(BigInt::from(0));
    // Find degree
    let d = ehrhart_coeffs.iter().rposition(|c| *c != zero).unwrap_or(0);

    // Evaluate L(t) at t = 0, 1, ..., d
    let mut values = Vec::with_capacity(d + 1);
    for t in 0..=d {
        let t_q = Q::from_integer(BigInt::from(t as i64));
        let mut val = zero.clone();
        let mut t_pow = Q::from_integer(BigInt::from(1));
        for c in ehrhart_coeffs.iter() {
            val += c.clone() * t_pow.clone();
            t_pow *= t_q.clone();
        }
        values.push(val);
    }

    // Binomial coefficients C(d+1, k) for k = 0, ..., d+1
    let mut binom = vec![BigInt::from(1); d + 2];
    for k in 1..=d + 1 {
        binom[k] = &binom[k - 1] * BigInt::from(d + 1 - k + 1) / BigInt::from(k);
    }

    // h*_i = sum_{k=0}^{i} (-1)^k * C(d+1, k) * L(i-k)
    let mut hstar = Vec::with_capacity(d + 1);
    for i in 0..=d {
        let mut val = zero.clone();
        for k in 0..=i {
            let signed_binom = if k % 2 == 0 {
                binom[k].clone()
            } else {
                -binom[k].clone()
            };
            val += Q::from_integer(signed_binom) * values[i - k].clone();
        }
        if val.denom() != &BigInt::from(1) {
            return Err(EhrhartToHstarError::NonIntegralEntry {
                index: i,
                value: val,
            });
        }
        hstar.push(val.to_integer());
    }

    Ok(hstar)
}

/// Convenience: convert Ehrhart polynomial given as `i64` coefficients of the
/// **numerator** polynomial when written over a common denominator.
///
/// That is, if L(t) = (a_0 + a_1 t + ... + a_d t^d) / denom, pass `(numerator_coeffs, denom)`.
/// This is useful when the Ehrhart polynomial is stored in integer-denominator form.
pub fn ehrhart_to_hstar_with_denom(
    numerator_coeffs: &[i64],
    denom: i64,
) -> Result<Vec<i64>, EhrhartToHstarError> {
    if denom == 0 {
        return Err(EhrhartToHstarError::ZeroDenominator);
    }
    let denom_q = Q::from_integer(BigInt::from(denom));
    let coeffs: Vec<Q> = numerator_coeffs
        .iter()
        .map(|&c| Q::from_integer(BigInt::from(c)) / denom_q.clone())
        .collect();
    ehrhart_to_hstar(&coeffs)
}

fn poly_degree_trimmed(coeffs: &[i64]) -> Option<usize> {
    coeffs.iter().rposition(|&c| c != 0)
}

fn poly_degree_trimmed_bigint(coeffs: &[BigInt]) -> Option<usize> {
    coeffs.iter().rposition(|c| *c != BigInt::from(0))
}

fn coeff_i64(coeffs: &[i64], k: usize) -> i64 {
    if k < coeffs.len() {
        coeffs[k]
    } else {
        0
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_poly() {
        assert_eq!(format_poly(&[1, 2, 1]), "1 + 2t + t^2");
        assert_eq!(format_poly(&[0, -1, 0, 3]), "-t + 3t^3");
        assert_eq!(format_poly(&[]), "0");
        assert_eq!(format_poly(&[0]), "0");
        assert_eq!(format_poly(&[5]), "5");
        assert_eq!(format_poly(&[0, 0, -1]), "-t^2");
        assert_eq!(format_poly_var(&[1, 1], "q"), "1 + q");
    }

    #[test]
    fn test_real_rooted() {
        assert!(is_real_rooted(&[1, 2, 1])); // (1+t)^2
        assert!(is_real_rooted(&[1, 11, 11, 1])); // Eulerian A_4(t)
        assert!(!is_real_rooted(&[1, 43, 196, 168, 23, 1])); // counterexample
        assert!(is_real_rooted(&[0])); // zero
        assert!(is_real_rooted(&[5])); // constant
        assert!(is_real_rooted(&[1, 1])); // linear
    }

    #[test]
    fn test_sturm_zero_polynomial_edge_cases() {
        assert!(is_real_rooted_sturm(&[]));
        assert!(is_real_rooted_sturm(&[0, 0, 0]));

        let empty: Vec<BigInt> = Vec::new();
        assert!(is_real_rooted_sturm_bigint_coeffs(&empty));
        assert!(is_real_rooted_sturm_bigint_coeffs(&[
            BigInt::from(0),
            BigInt::from(0),
        ]));
    }

    #[test]
    fn test_log_concave() {
        assert!(is_log_concave(&[1, 2, 1]));
        assert!(is_log_concave(&[1, 3, 3, 1]));
        assert!(!is_log_concave(&[1, 1, 3]));
    }

    #[test]
    fn test_unimodal() {
        assert!(is_unimodal(&[]));
        assert!(is_unimodal(&[1]));
        assert!(is_unimodal(&[1, 2, 2, 1]));
        assert!(is_unimodal(&[3, 2, 1]));
        assert!(is_unimodal(&[1, 2, 3]));
        assert!(!is_unimodal(&[1, 3, 2, 4]));
    }

    #[test]
    fn test_log_concave_uses_wide_products() {
        assert!(is_log_concave(&[i64::MAX, i64::MAX, i64::MAX]));
        assert!(!is_log_concave(&[i64::MAX, 0, i64::MAX]));
    }

    #[test]
    fn test_ultra_log_concave() {
        assert!(is_ultra_log_concave(&[1, 3, 3, 1])); // (1+t)^3
        assert!(is_ultra_log_concave(&[1, 4, 6, 4, 1])); // (1+t)^4
    }

    #[test]
    fn test_bigint_shape_properties() {
        let scale = BigInt::from(10u64).pow(40);
        let coeffs = [1_i64, 3, 3, 1]
            .into_iter()
            .map(|c| BigInt::from(c) * &scale)
            .collect::<Vec<_>>();
        assert!(is_palindromic_bigint_coeffs(&coeffs));
        assert!(is_unimodal_bigint_coeffs(&coeffs));
        assert!(is_log_concave_bigint_coeffs(&coeffs));
        assert!(is_ultra_log_concave_bigint_coeffs(&coeffs));
    }

    #[test]
    fn test_hermite_biehler_parts() {
        let (even, odd) = hermite_biehler_parts(&[1, 2, 3, 4, 5]);
        assert_eq!(even, vec![1, 3, 5]);
        assert_eq!(odd, vec![2, 4]);
    }

    #[test]
    fn test_has_simple_roots_helpers() {
        assert!(!has_simple_roots(&[1, 2, 1])); // (1+t)^2
        assert!(has_simple_roots(&[1, 0, -1])); // 1 - t^2
        assert!(has_simple_roots(&[7])); // constant
        assert!(!has_simple_roots(&[0])); // zero

        let repeated: Vec<BigInt> = [1_i64, 2, 1].into_iter().map(BigInt::from).collect();
        assert!(!has_simple_roots_bigint_coeffs(&repeated));

        let simple: Vec<BigInt> = [1_i64, 0, -1].into_iter().map(BigInt::from).collect();
        assert!(has_simple_roots_bigint_coeffs(&simple));
    }

    #[test]
    fn test_real_roots() {
        // (t-1)(t-2) = 2 - 3t + t^2
        let roots = real_roots(&[2, -3, 1]).unwrap();
        assert_eq!(roots.len(), 2);

        // t^2 + 1: not real-rooted
        assert!(real_roots(&[1, 0, 1]).is_none());

        // Initial zero coefficients represent a root at zero, not coefficients
        // that may simply be discarded from the returned root set.
        assert_eq!(real_roots(&[0, 1]), Some(vec![Q::zero()]));
        assert_eq!(real_roots(&[0, 0, 1]), Some(vec![Q::zero()]));
    }

    #[test]
    fn test_interlacing_sturm() {
        // (t-2) ≪ (t-1)(t-3): roots 1 < 2 < 3
        assert_eq!(check_interlacing_sturm(&[-2, 1], &[3, -4, 1]), Some(true));

        // Wrong order: big first → Some(false) since deg diff is -1
        assert_eq!(check_interlacing_sturm(&[3, -4, 1], &[-2, 1]), Some(false));

        // Same-degree: not interlacing (requires deg diff = 1)
        assert_eq!(
            check_interlacing_sturm(&[-1, 0, 1], &[-4, 0, 1]),
            Some(false)
        );

        // Not real-rooted: t^2+1 → None
        assert_eq!(check_interlacing_sturm(&[1, 1], &[1, 0, 1]), None);

        // The linear root 333333334/1000000000 lies just to the right of
        // 1/3, the smaller root of 2 - 7t + 3t^2. Fixed-width midpoint
        // approximations used to reverse this ordering.
        assert_eq!(
            check_interlacing_sturm(&[-333_333_334, 1_000_000_000], &[2, -7, 3]),
            Some(true)
        );
    }

    #[test]
    fn test_interlacing_bezout() {
        // (t-2) ≪ (t-1)(t-3): interlacing, roots 1 < 2 < 3
        assert_eq!(check_interlacing(&[-2, 1], &[3, -4, 1]), Some(true));

        // Wrong order: big first → None
        assert_eq!(check_interlacing(&[3, -4, 1], &[-2, 1]), None);

        // Same-degree polynomials: returns None
        assert_eq!(check_interlacing(&[-1, 0, 1], &[-4, 0, 1]), None);

        // Not real-rooted: t+1 and t^2+1. Degree diff = 1 but B is not pos def.
        assert_eq!(check_interlacing(&[1, 1], &[1, 0, 1]), Some(false));
    }

    #[test]
    fn test_interlacing_is_invariant_under_scalar_signs() {
        // (t-2) interlaces (t-1)(t-3), independent of multiplying either
        // polynomial by -1.
        assert_eq!(check_interlacing(&[-2, 1], &[3, -4, 1]), Some(true));
        assert_eq!(check_interlacing(&[2, -1], &[3, -4, 1]), Some(true));
        assert_eq!(check_interlacing(&[-2, 1], &[-3, 4, -1]), Some(true));
        assert_eq!(check_interlacing(&[2, -1], &[-3, 4, -1]), Some(true));
    }

    #[test]
    fn test_interlacing_bigint_coeffs() {
        let scale = BigInt::from(10u64).pow(40);
        let f: Vec<BigInt> = [2_i64, -1]
            .into_iter()
            .map(|c| BigInt::from(c) * &scale)
            .collect();
        let g: Vec<BigInt> = [-3_i64, 4, -1]
            .into_iter()
            .map(|c| BigInt::from(c) * &scale)
            .collect();

        assert_eq!(check_interlacing_bigint_coeffs(&f, &g), Some(true));
    }

    // -- Bézout matrix tests --

    #[test]
    fn test_bezout_matrix_simple() {
        // f = (t-1)(t-3) = 3 - 4t + t^2, g = t - 2
        // B should be 2x2... wait deg(f)=2, deg(g)=1, so d=2.
        let b = bezout_matrix(&[3, -4, 1], &[-2, 1]).unwrap();
        assert_eq!(b.len(), 2);
        // Should be positive definite since g interlaces f (roots: 1 < 2 < 3).
        let b_big = bezout_matrix_bigint(&[3, -4, 1], &[-2, 1]).unwrap();
        assert!(is_positive_definite_bigint(&b_big));
    }

    #[test]
    fn test_bezout_matrix_checked_narrowing() {
        let max = i64::MAX;
        let min = i64::MIN;
        assert!(bezout_matrix(&[max, min, max], &[min, max]).is_none());
    }

    #[test]
    fn test_bezout_interlacing() {
        // (t-2) ≪ (t-1)(t-3): interlacing, roots 1 < 2 < 3
        assert_eq!(check_interlacing_bezout(&[-2, 1], &[3, -4, 1]), Some(true));

        // (t-4) ≪ (t-1)(t-3): NOT interlacing (4 not between 1 and 3)
        assert_eq!(check_interlacing_bezout(&[-4, 1], &[3, -4, 1]), Some(false));
    }

    #[test]
    fn test_bezout_interlacing_degree3() {
        // (t-2)(t-4) ≪ (t-1)(t-3)(t-5): interlacing, 1 < 2 < 3 < 4 < 5
        assert_eq!(
            check_interlacing_bezout(&[8, -6, 1], &[-15, 23, -9, 1]),
            Some(true)
        );

        // (t-2)(t-6) ≪ (t-1)(t-3)(t-5): NOT interlacing (no root of g between 3 and 5)
        assert_eq!(
            check_interlacing_bezout(&[12, -8, 1], &[-15, 23, -9, 1]),
            Some(false)
        );
    }

    #[test]
    fn test_bezout_real_rooted() {
        assert!(is_real_rooted_bezout(&[1, 2, 1])); // (1+t)^2
        assert!(is_real_rooted_bezout(&[1, 3, 3, 1])); // (1+t)^3
        assert!(is_real_rooted_bezout(&[1, 11, 11, 1])); // Eulerian A_4(t)
        assert!(!is_real_rooted_bezout(&[1, 0, 1])); // t^2 + 1
        assert!(is_real_rooted_bezout(&[1])); // constant
        assert!(is_real_rooted_bezout(&[1, 1])); // linear
    }

    #[test]
    fn test_default_real_rooted_bigint_coeffs() {
        let rr: Vec<BigInt> = [1_i64, 3, 3, 1].into_iter().map(BigInt::from).collect();
        assert!(is_real_rooted_bigint_coeffs(&rr));

        let not_rr: Vec<BigInt> = [1_i64, 0, 1].into_iter().map(BigInt::from).collect();
        assert!(!is_real_rooted_bigint_coeffs(&not_rr));

        let huge_scale = BigInt::from(10u64).pow(40);
        let scaled: Vec<BigInt> = [1_i64, 2, 1]
            .into_iter()
            .map(|c| BigInt::from(c) * &huge_scale)
            .collect();
        assert!(is_real_rooted_bigint_coeffs(&scaled));
    }

    #[test]
    fn test_default_real_rootedness_uses_prs_for_mixed_signs() {
        let cases: Vec<(&[i64], bool)> = vec![
            (&[-6, 11, -6, 1], true),  // (t-1)(t-2)(t-3)
            (&[-2, 1, -2, 1], false),  // (t-2)(t^2+1)
            (&[0, 0, 1, -2, 1], true), // t^2(t-1)^2
            (&[1, 0, 2, 0, 1], false), // (t^2+1)^2
        ];

        for (coeffs, expected) in cases {
            let big: Vec<BigInt> = coeffs.iter().copied().map(BigInt::from).collect();
            assert_eq!(
                is_real_rooted_bigint_coeffs(&big),
                crate::root_count::is_real_rooted_fast_bigint_coeffs(&big),
                "default should delegate to PRS backend for {coeffs:?}"
            );
            assert_eq!(
                is_real_rooted(coeffs),
                expected,
                "i64 default disagrees on {coeffs:?}"
            );
            assert_eq!(
                is_real_rooted_bigint_coeffs(&big),
                expected,
                "BigInt default disagrees on {coeffs:?}"
            );
        }
    }

    #[test]
    fn test_squarefree_part_bigint_coeffs() {
        let repeated: Vec<BigInt> = [0_i64, 0, 1, -2, 1].into_iter().map(BigInt::from).collect();
        assert_eq!(
            squarefree_part_bigint_coeffs(&repeated),
            vec![BigInt::from(0), BigInt::from(-1), BigInt::from(1)]
        );

        let complex_repeated: Vec<BigInt> =
            [1_i64, 0, 2, 0, 1].into_iter().map(BigInt::from).collect();
        assert_eq!(
            squarefree_part_bigint_coeffs(&complex_repeated),
            vec![BigInt::from(1), BigInt::from(0), BigInt::from(1)]
        );
    }

    #[test]
    fn test_hermite_matrix_bigint_coeffs() {
        // (t - 1)(t - 2) has Newton sums s_0=2, s_1=3, s_2=5.
        let p: Vec<BigInt> = [2_i64, -3, 1].into_iter().map(BigInt::from).collect();
        assert_eq!(
            hermite_matrix_bigint_coeffs(&p),
            Some(vec![
                vec![BigInt::from(2), BigInt::from(3)],
                vec![BigInt::from(3), BigInt::from(5)],
            ])
        );

        // Non-monic example: roots are still 1 and 2, and denominators are
        // cleared by a positive scalar.
        let scaled: Vec<BigInt> = [4_i64, -6, 2].into_iter().map(BigInt::from).collect();
        assert_eq!(
            hermite_matrix_bigint_coeffs(&scaled),
            Some(vec![
                vec![BigInt::from(2), BigInt::from(3)],
                vec![BigInt::from(3), BigInt::from(5)],
            ])
        );
    }

    #[test]
    fn test_squarefree_real_rootedness_paths() {
        let cases: Vec<(&[i64], bool)> = vec![
            (&[1, 2, 1], true),        // repeated real root
            (&[4, -4, 1], true),       // non-monic repeated real root
            (&[0, 0, 1, -2, 1], true), // t^2(t-1)^2
            (&[-2, 3, -1], true),      // negative scalar times (t-1)(t-2)
            (&[1, 0, 1], false),       // t^2+1
            (&[1, 0, 2, 0, 1], false), // (t^2+1)^2
            (&[1, 43, 196, 168, 23, 1], false),
            (&[1, 11, 11, 1], true),
        ];

        for (coeffs, expected) in cases {
            assert_eq!(
                is_real_rooted_bezout_squarefree(coeffs),
                expected,
                "squarefree Bézout disagrees on {coeffs:?}"
            );
            assert_eq!(
                is_real_rooted_hermite(coeffs),
                expected,
                "Hermite/Newton disagrees on {coeffs:?}"
            );
            assert_eq!(
                is_real_rooted_bezout(coeffs),
                expected,
                "old Bézout disagrees on {coeffs:?}"
            );
        }

        assert!(!is_strictly_real_rooted_bezout(&[1, 2, 1]));
        assert!(!is_strictly_real_rooted_hermite(&[1, 2, 1]));
        assert!(is_strictly_real_rooted_bezout(&[2, -3, 1]));
        assert!(is_strictly_real_rooted_hermite(&[2, -3, 1]));
    }

    #[test]
    fn test_weak_interlacing_shared_roots() {
        // g = (t-1)(t-2) = 2 - 3t + t^2, roots {1, 2} (deg 2, small)
        // f = (t-1)^2 (t-3) = -3 + 7t - 5t^2 + t^3, roots {1, 1, 3} (deg 3, big)
        // Shared root at 1. After dividing out (t-1):
        //   g/(t-1) = (t-2) = -2 + t   (deg 1, small)
        //   f/(t-1) = (t-1)(t-3) = 3 - 4t + t^2  (deg 2, big)
        // These strictly interlace: 1 < 2 < 3.
        assert_eq!(
            check_weak_interlacing_bezout(&[2, -3, 1], &[-3, 7, -5, 1]),
            Some(true)
        );

        // f = (t-1)(t-3) = 3-4t+t^2, g = (t-1)(t-4) = 4-5t+t^2
        // Same degree after GCD removal: f/(t-1) = t-3, g/(t-1) = t-4
        // These have degree diff = 0 with 1 root each. Strict check doesn't apply.
        // But they do weakly interlace if we allow same-degree.
        // Actually with deg diff = 0, check_interlacing returns None.
        // So weak interlacing for same-degree reduced polys returns None.
        // Same degree. Roots of f: {1,3}, roots of g: {1,4}.
        // After removing shared root 1, reduced roots {3} and {4}: same degree.
        // The Cauchy extension puts a root far right: 3 < 4 < R, so interlacing holds.
        let result = check_weak_interlacing_bezout(&[3, -4, 1], &[4, -5, 1]);
        assert_eq!(result, Some(true));

        // Same-degree interlacing that DOES hold:
        // f = (t-1)(t-3), roots {1,3}; g = (t-2)(t-4), roots {2,4}
        // Roots alternate: 1 < 2 < 3 < 4.
        assert_eq!(check_weak_interlacing(&[3, -4, 1], &[8, -6, 1]), Some(true));

        // The check is directed, so reversing the order changes the answer.
        assert_eq!(
            check_weak_interlacing(&[8, -6, 1], &[3, -4, 1]),
            Some(false)
        );

        // Same-degree, NOT interlacing: f = (t-1)(t-4), g = (t-2)(t-3)
        // Roots nested: 1 < 2 < 3 < 4 but pattern is f,g,g,f — not alternating.
        assert_eq!(
            check_weak_interlacing(&[4, -5, 1], &[6, -5, 1]),
            Some(false)
        );

        // Same-degree, positive coefficients (fast path via *t):
        // f = 1+t, g = 1+2t. Roots: -1 < -1/2. Alternating.
        assert_eq!(check_weak_interlacing(&[1, 1], &[1, 2]), Some(true));

        // Same-degree, positive left polynomial but right polynomial has a
        // positive root. The fast root-at-0 reduction is not valid here; the
        // Cauchy-bound reduction must put the added root to the right of 1.
        assert_eq!(check_weak_interlacing(&[1, 1], &[-1, 1]), Some(true));
        assert_eq!(check_weak_interlacing(&[-1, 1], &[1, 1]), Some(false));
    }

    #[test]
    fn test_weak_interlacing_identical() {
        // f = g = (t-1)(t-2): all roots shared, reduced to constants.
        assert_eq!(
            check_weak_interlacing_bezout(&[2, -3, 1], &[2, -3, 1]),
            Some(true)
        );
    }

    #[test]
    fn test_weak_interlacing_bigint_coeffs() {
        let scale = BigInt::from(10u64).pow(40);
        let p: Vec<BigInt> = [1_i64, 1]
            .into_iter()
            .map(|c| BigInt::from(c) * &scale)
            .collect();
        let q: Vec<BigInt> = [-1_i64, 1]
            .into_iter()
            .map(|c| BigInt::from(c) * &scale)
            .collect();

        assert_eq!(check_weak_interlacing_bigint_coeffs(&p, &q), Some(true));
        assert_eq!(check_weak_interlacing_bigint_coeffs(&q, &p), Some(false));
    }

    #[test]
    fn test_weak_interlacing_complex_shared_roots() {
        // f = (t^2+1)(t-1) = -1 + t - t^2 + t^3, complex roots at ±i, real root at 1
        // g = (t^2+1)(t-2) = -2 + t - 2t^2 + t^3, complex roots at ±i, real root at 2
        // GCD = t^2+1, which is NOT real-rooted.
        // Should return Some(false) — cannot weakly interlace with shared complex roots.
        assert_eq!(
            check_weak_interlacing(&[-1, 1, -1, 1], &[-2, 1, -2, 1]),
            Some(false)
        );
    }

    #[test]
    fn test_bezout_agrees_with_sturm() {
        // Test on several polynomials that both methods agree
        let cases: Vec<(&[i64], bool)> = vec![
            (&[1, 2, 1], true),
            (&[1, 3, 3, 1], true),
            (&[1, 4, 6, 4, 1], true),
            (&[1, 11, 11, 1], true),
            (&[1, 26, 66, 26, 1], true),
            (&[1, 0, 1], false),
            (&[1, 43, 196, 168, 23, 1], false),
            (&[-6, 11, -6, 1], true), // (t-1)(t-2)(t-3)
        ];
        for (coeffs, expected) in cases {
            let sturm = is_real_rooted_sturm(coeffs);
            let bezout = is_real_rooted(coeffs);
            assert_eq!(sturm, expected, "Sturm disagrees on {:?}", coeffs);
            assert_eq!(bezout, expected, "Bézout disagrees on {:?}", coeffs);
        }
    }

    // -- Palindromic / gamma tests --

    #[test]
    fn test_is_palindromic() {
        assert!(is_palindromic(&[1, 2, 1]));
        assert!(is_palindromic(&[1, 11, 11, 1]));
        assert!(is_palindromic(&[1, 4, 6, 4, 1]));
        assert!(is_palindromic(&[1]));
        assert!(is_palindromic(&[]));
        assert!(!is_palindromic(&[0, 1, 1]));
        assert!(!is_palindromic(&[0, 0, 1, 2, 1]));
        assert!(!is_palindromic(&[1, 2, 3]));
        // Trailing zeros trimmed: [1, 2, 1, 0] is palindromic
        assert!(is_palindromic(&[1, 2, 1, 0]));
        assert!(!is_palindromic(&[0, 1, 2, 1, 0]));
    }

    #[test]
    fn test_gamma_coefficients() {
        // (1+t)^3 = [1,3,3,1]: gamma = [1, 0]
        assert_eq!(gamma_coefficients(&[1, 3, 3, 1]), Some(vec![1, 0]));

        // A_4(t) = [1,11,11,1]: gamma_0=1, gamma_1=11-3=8
        assert_eq!(gamma_coefficients(&[1, 11, 11, 1]), Some(vec![1, 8]));

        // (1+t)^4 = [1,4,6,4,1]: gamma_0=1, gamma_1=4-C(4,1)=0, gamma_2=6-C(4,2)=0
        assert_eq!(gamma_coefficients(&[1, 4, 6, 4, 1]), Some(vec![1, 0, 0]));

        // Not palindromic
        assert_eq!(gamma_coefficients(&[1, 2, 3]), None);

        // Zero/constant
        assert_eq!(gamma_coefficients(&[]), Some(vec![]));
        assert_eq!(gamma_coefficients(&[5]), Some(vec![5]));
        assert_eq!(gamma_coefficients(&[0, 1, 1]), None);
        assert_eq!(gamma_coefficients(&[0, 0, 1, 3, 3, 1]), None);

        // [1,2,2,1] degree 3: gamma_0=1, gamma_1=2-3=-1
        assert_eq!(gamma_coefficients(&[1, 2, 2, 1]), Some(vec![1, -1]));

        // The exact gamma expansion exists, but later coefficients do not fit
        // in i64. The convenience wrapper must fail cleanly without overflowing
        // its former i128 Pascal triangle.
        let mut high_degree = vec![0; 201];
        high_degree[0] = 1;
        high_degree[200] = 1;
        assert_eq!(gamma_coefficients(&high_degree), None);
    }

    #[test]
    fn test_bigint_gamma_coefficients() {
        let scale = BigInt::from(10u64).pow(40);
        let coeffs = [1_i64, 11, 11, 1]
            .into_iter()
            .map(|c| BigInt::from(c) * &scale)
            .collect::<Vec<_>>();
        assert_eq!(
            gamma_coefficients_bigint_coeffs(&coeffs),
            Some(vec![scale.clone(), BigInt::from(8) * scale])
        );
    }

    #[test]
    fn test_is_gamma_positive() {
        assert!(is_gamma_positive(&[1, 3, 3, 1]));
        assert!(is_gamma_positive(&[1, 11, 11, 1]));
        assert!(is_gamma_positive(&[1, 4, 6, 4, 1]));
        assert!(!is_gamma_positive(&[0, 1, 1]));
        // [1,2,2,1] has gamma_1 = -1, NOT gamma-positive
        assert!(!is_gamma_positive(&[1, 2, 2, 1]));
        // Non-palindromic: not gamma-positive
        assert!(!is_gamma_positive(&[1, 2, 3]));
    }

    #[test]
    fn test_strip_initial_zero_wrappers() {
        let empty: &[i64] = &[];
        assert_eq!(strip_initial_zeros(empty), empty);
        assert_eq!(strip_initial_zeros(&[0, 0, 1, 2, 1]), &[1, 2, 1]);
        assert!(is_palindromic_ignoring_initial_zeros(&[0, 1, 1]));
        assert!(is_palindromic_ignoring_initial_zeros(&[0, 0, 1, 2, 1]));
        assert_eq!(
            gamma_coefficients_ignoring_initial_zeros(&[0, 1, 1]),
            Some(vec![1])
        );
        assert_eq!(
            gamma_coefficients_ignoring_initial_zeros(&[0, 0, 1, 3, 3, 1]),
            Some(vec![1, 0])
        );
        assert!(is_gamma_positive_ignoring_initial_zeros(&[0, 1, 1]));
    }

    #[test]
    fn test_stapledon_decomposition_basic() {
        let (a, b) = stapledon_decomposition(&[1, 2, 3], 2).unwrap();
        assert_eq!(a, vec![1, 0, 1]);
        assert_eq!(b, vec![2, 2]);
    }

    #[test]
    fn test_stapledon_decomposition_palindromic() {
        let (a, b) = stapledon_decomposition(&[1, 11, 11, 1], 3).unwrap();
        assert_eq!(a, vec![1, 11, 11, 1]);
        assert_eq!(b, Vec::<i64>::new());
    }

    #[test]
    fn test_stapledon_decomposition_with_larger_bound() {
        let (a, b) = stapledon_decomposition(&[1, 4, 1], 4).unwrap();
        assert_eq!(a, vec![1, 5, 6, 5, 1]);
        assert_eq!(b, vec![-1, -5, -5, -1]);
    }

    #[test]
    fn test_stapledon_decomposition_rejects_large_degree() {
        assert_eq!(stapledon_decomposition(&[1, 2, 3], 1), None);
    }

    // -- Resultant / discriminant tests --

    #[test]
    fn test_resultant() {
        // Res((t-1)(t-2), (t-3)) = (1-3)(2-3) = (-2)(-1) = 2
        assert_eq!(resultant(&[2, -3, 1], &[-3, 1]), BigInt::from(2));

        // Res(t^2+1, t-1) = 1^2 + 1 = 2
        assert_eq!(resultant(&[1, 0, 1], &[-1, 1]), BigInt::from(2));

        // Res(t-1, t-1) = 0 (shared root)
        assert_eq!(resultant(&[-1, 1], &[-1, 1]), BigInt::from(0));
    }

    #[test]
    fn test_discriminant() {
        // disc(t^2 - 1) = disc((t-1)(t+1)). Roots 1,-1.
        // disc = (-1)^{2*1/2} * Res(t^2-1, 2t) / 1
        // = (-1)^1 * Res([−1,0,1],[0,2]) / 1
        // Sylvester of t^2-1 (deg 2) and 2t (deg 1): 3x3 matrix.
        //   row0 (f shifted): [1, 0, -1]
        //   row1 (g shifted): [2, 0, 0]
        //   row2 (g shifted): [0, 2, 0]
        // det = 1*(0*0 - 0*2) - 0*(2*0 - 0*0) + (-1)*(2*2 - 0*0) = -4
        // disc = (-1)^1 * (-4) / 1 = 4
        assert_eq!(discriminant(&[-1, 0, 1]), BigInt::from(4));

        // disc((t-1)(t-2)) = (1-2)^2 = 1. f = [2,-3,1], lc=1, n=2.
        // disc = (-1)^1 * Res(f, f') / lc
        // f' = [-3, 2]. Sylvester 3x3.
        assert_eq!(discriminant(&[2, -3, 1]), BigInt::from(1));

        // disc(t^2 + 1) = -4 (no real roots, negative discriminant)
        assert_eq!(discriminant(&[1, 0, 1]), BigInt::from(-4));
    }

    #[test]
    fn test_discriminant_derivative_uses_bigint() {
        // f(t) = i64::MAX * t^2 + 1 has derivative 2*i64::MAX*t,
        // which does not fit in i64.
        assert_eq!(
            discriminant(&[1, 0, i64::MAX]),
            BigInt::from(-4) * BigInt::from(i64::MAX)
        );
    }

    #[test]
    fn test_discriminant_degree3() {
        // f = (t-1)(t-2)(t-3) = -6 + 11t - 6t^2 + t^3
        // disc = product (r_i - r_j)^2 = (1-2)^2(1-3)^2(2-3)^2 = 1*4*1 = 4
        // But with the sign: (-1)^{3*2/2} = (-1)^3 = -1, and
        // disc = (-1)^3 * Res(f,f') / lc(f).
        // Let's just compute and check the absolute value.
        let d = discriminant(&[-6, 11, -6, 1]);
        // |disc| = 4 for monic with roots 1,2,3
        assert_eq!(d.clone() * d.clone().signum(), BigInt::from(4));
    }

    // -- Ehrhart / h*-vector tests --

    fn q(n: i64, d: i64) -> Q {
        Q::new(BigInt::from(n), BigInt::from(d))
    }

    #[test]
    fn test_hstar_to_ehrhart_simplex() {
        // Standard d-simplex has h* = [1, 0, ..., 0].
        // L(t) = C(t+d, d).

        // 1-simplex [0,1]: h* = [1, 0], L(t) = t + 1
        let ehrhart = hstar_to_ehrhart(&[1, 0]);
        assert_eq!(ehrhart.len(), 2);
        assert_eq!(ehrhart[0], q(1, 1)); // constant term = 1
        assert_eq!(ehrhart[1], q(1, 1)); // t coefficient = 1

        // 2-simplex: h* = [1, 0, 0], L(t) = C(t+2,2) = (t+1)(t+2)/2 = 1 + 3t/2 + t^2/2
        let ehrhart = hstar_to_ehrhart(&[1, 0, 0]);
        assert_eq!(ehrhart.len(), 3);
        assert_eq!(ehrhart[0], q(1, 1)); // 1
        assert_eq!(ehrhart[1], q(3, 2)); // 3/2
        assert_eq!(ehrhart[2], q(1, 2)); // 1/2
    }

    #[test]
    fn test_hstar_to_ehrhart_square() {
        // Unit square: h* = [1, 1], L(t) = (t+1)^2 = 1 + 2t + t^2
        let ehrhart = hstar_to_ehrhart(&[1, 1]);
        assert_eq!(ehrhart.len(), 2);
        // Wait, h* = [1, 1] means d = 1, but the unit square is dimension 2...
        // Actually h* has d+1 entries for a d-dimensional polytope.
        // Unit square: d=2, h* = [1, 1, 0]
        let ehrhart = hstar_to_ehrhart(&[1, 1, 0]);
        assert_eq!(ehrhart.len(), 3);
        assert_eq!(ehrhart[0], q(1, 1));
        assert_eq!(ehrhart[1], q(2, 1));
        assert_eq!(ehrhart[2], q(1, 1));
    }

    #[test]
    fn test_ehrhart_to_hstar_simplex() {
        // 1-simplex: L(t) = t + 1 = [1, 1], h* = [1, 0]
        let ehrhart = vec![q(1, 1), q(1, 1)];
        assert_eq!(ehrhart_to_hstar(&ehrhart).unwrap(), vec![1, 0]);

        // 2-simplex: L(t) = 1 + 3t/2 + t^2/2, h* = [1, 0, 0]
        let ehrhart = vec![q(1, 1), q(3, 2), q(1, 2)];
        assert_eq!(ehrhart_to_hstar(&ehrhart).unwrap(), vec![1, 0, 0]);
    }

    #[test]
    fn test_ehrhart_to_hstar_square() {
        // Unit square: L(t) = (t+1)^2 = 1 + 2t + t^2, h* = [1, 1, 0]
        let ehrhart = vec![q(1, 1), q(2, 1), q(1, 1)];
        assert_eq!(ehrhart_to_hstar(&ehrhart).unwrap(), vec![1, 1, 0]);
    }

    #[test]
    fn test_roundtrip_hstar_ehrhart() {
        // h* -> Ehrhart -> h* should be identity
        let hstar = vec![1, 4, 6, 4, 1];
        let ehrhart = hstar_to_ehrhart(&hstar);
        let recovered = ehrhart_to_hstar(&ehrhart).unwrap();
        assert_eq!(recovered, hstar);
    }

    #[test]
    fn test_ehrhart_to_hstar_with_denom() {
        // 2-simplex: L(t) = (t^2 + 3t + 2) / 2
        // numerator = [2, 3, 1], denom = 2
        assert_eq!(
            ehrhart_to_hstar_with_denom(&[2, 3, 1], 2).unwrap(),
            vec![1, 0, 0]
        );
    }

    #[test]
    fn test_ehrhart_to_hstar_reports_invalid_inputs() {
        assert_eq!(
            ehrhart_to_hstar_bigint(&[q(1, 2)]),
            Err(EhrhartToHstarError::NonIntegralEntry {
                index: 0,
                value: q(1, 2),
            })
        );
        assert_eq!(
            ehrhart_to_hstar_with_denom(&[1], 0),
            Err(EhrhartToHstarError::ZeroDenominator)
        );
    }
}
