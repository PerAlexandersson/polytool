//! Interlacing checks: strict Bézout interlacing, weak (shared-root and
//! same-degree) interlacing, and the undirected `interlaces_either_order`.
//!
//! Directed functions return `Option<bool>`; `None` means the degree relation
//! does not admit that directed test.  Shared factors are divided out over the
//! rationals with the helpers in the sibling `qpoly` module.

use super::bezout::{bezout_matrix_bigint_coeffs, is_positive_definite_bigint};
use super::qpoly::{poly_exact_div_q, poly_gcd_q, q_degree, q_poly_to_primitive_bigint};
use super::real_rooted::is_real_rooted_bigint_coeffs;
use super::Q;
use crate::vec_poly;
use num_bigint::BigInt;
use num_traits::{Signed, Zero};

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
    let df = vec_poly::degree(f);
    let dg = vec_poly::degree(g);

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

/// Check **directed weak** interlacing `p ≪ q` ("`p` interlaces `q`") exactly.
///
/// # Definition
///
/// Let `p` and `q` be real-rooted, with roots `a_1 ≥ a_2 ≥ ... ≥ a_m` and
/// `b_1 ≥ b_2 ≥ ...` listed with multiplicity.  Then `p ≪ q` means:
///
/// - **`deg q = deg p + 1`**: the roots of `p` separate those of `q`,
///
///   ```text
///   b_1 ≥ a_1 ≥ b_2 ≥ a_2 ≥ ... ≥ b_m ≥ a_m ≥ b_(m+1);
///   ```
///
/// - **`deg q = deg p`** (same degree): every root of `p` lies weakly to the
///   left of the matching root of `q`,
///
///   ```text
///   b_1 ≥ a_1 ≥ b_2 ≥ a_2 ≥ ... ≥ b_m ≥ a_m.
///   ```
///
/// All inequalities are weak, so shared (and repeated) roots are allowed.  The
/// relation only depends on the roots: multiplying either polynomial by a
/// nonzero constant (in particular, flipping the sign of its leading
/// coefficient) does not change the answer.
///
/// **This relation is directed.**  For same-degree pairs `p ≪ q` and `q ≪ p`
/// are different conditions (both hold only when `p` and `q` have the same
/// roots), and for `deg q = deg p + 1` only `p ≪ q` is meaningful.  Use
/// [`interlaces_either_order`] for the undirected question "do these two
/// polynomials interlace in some order?".
///
/// # Zero and constant polynomials
///
/// - `p = 0`: returns `Some(is_real_rooted(q))`; the zero polynomial
///   interlaces every real-rooted polynomial, including `0` itself.
/// - `q = 0` and `p ≠ 0`: returns `Some(false)`.
/// - Nonzero constants have no roots, so `c ≪ q` holds for every nonzero
///   constant `c` and every real-rooted `q` of degree 0 or 1.
///
/// Some research conventions instead treat a pair with a zero member as
/// vacuously interlacing in *both* positions, whatever the other polynomial
/// is.  No separate function is needed for that; write
///
/// ```text
/// vec_poly::is_zero(p) || vec_poly::is_zero(q) || check_weak_interlacing(p, q) == Some(true)
/// ```
///
/// (and use [`interlaces_either_order`] for the undirected question, where a
/// zero argument interlaces exactly the real-rooted polynomials).
///
/// # Return value
///
/// - `Some(true)`: the degrees are compatible and `p ≪ q` holds (so both
///   polynomials are real-rooted).
/// - `Some(false)`: the degrees are compatible but `p ≪ q` fails, including
///   when `p` or `q` is not real-rooted.
/// - `None`: both polynomials are nonzero and `deg q ∉ {deg p, deg p + 1}`;
///   the directed relation is not defined for this argument order.  In
///   particular, `None` is returned when the arguments are passed in the
///   wrong order for a degree-difference-one pair.
///
/// # Method
///
/// - **`deg q = deg p + 1`**: divides out `gcd(p, q)` over ℚ if it is
///   nontrivial (verifying that the shared factor is real-rooted), then checks
///   strict interlacing of the reduced pair with the Bézout matrix.
/// - **`deg q = deg p`**: extends `p` with a root to the right of all roots
///   (via the Cauchy bound), giving `p_ext` of degree `deg p + 1`, then checks
///   `q ≪ p_ext`.  **Fast path**: if both `p` and `q` have one-signed
///   coefficients with nonzero constant terms (so no nonnegative roots),
///   multiplying `p` by `t` suffices.
///
/// References: Kummer--Naldi--Plaumann, Thm. 2.13 for the
/// degree-difference-one Bézoutian criterion; Fisk, §9.21, Cor. 9.145 for the
/// same-degree oriented alternation criterion. The explicit gcd real-rootedness
/// check guards the semidefinite/common-factor caveat.
///
/// ```
/// use polytool::check_weak_interlacing;
///
/// // t ≪ t^2 - 1: the root 0 separates -1 and 1.
/// assert_eq!(check_weak_interlacing(&[0, 1], &[-1, 0, 1]), Some(true));
/// // Wrong argument order for a degree-difference-one pair.
/// assert_eq!(check_weak_interlacing(&[-1, 0, 1], &[0, 1]), None);
/// // Same degree: t + 1 (root -1) lies to the left of t (root 0).
/// assert_eq!(check_weak_interlacing(&[1, 1], &[0, 1]), Some(true));
/// assert_eq!(check_weak_interlacing(&[0, 1], &[1, 1]), Some(false));
/// ```
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
    if vec_poly::degree(p).is_none() {
        return Some(is_real_rooted_bigint_coeffs(q));
    }
    // p ≪ 0 only if p is also zero (already handled above).
    if vec_poly::degree(q).is_none() {
        return Some(false);
    }

    let dp = vec_poly::degree(p).unwrap();
    let dq = vec_poly::degree(q).unwrap();

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
            vec_poly::shift(p)
        } else {
            let r = cauchy_root_bound_bigint(p).max(cauchy_root_bound_bigint(q)) + BigInt::from(1);
            vec_poly::mul(p, &[-r, BigInt::from(1)])
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

/// Undirected weak interlacing: whether `f` and `g` interlace in some order.
///
/// Returns `true` iff one of the directed relations of
/// [`check_weak_interlacing`] holds in an order allowed by the degrees:
///
/// | degrees | condition |
/// |---|---|
/// | `deg g = deg f + 1` | `f ≪ g` |
/// | `deg f = deg g + 1` | `g ≪ f` |
/// | `deg f = deg g` | `f ≪ g` **or** `g ≪ f` |
/// | `|deg f - deg g| ≥ 2` | never (`false`) |
///
/// The zero polynomial interlaces every real-rooted polynomial (in both
/// orders), so `interlaces_either_order(0, g) == is_real_rooted(g)` and
/// `interlaces_either_order(0, 0) == true`; a polynomial that is not
/// real-rooted interlaces nothing.  Leading-coefficient signs do not matter,
/// and nonzero constants interlace each other and every linear polynomial.
///
/// This is the predicate that research code usually wants when checking
/// that consecutive members of a family interlace, without tracking which
/// one has the larger degree or which same-degree orientation applies.
///
/// ```
/// use polytool::interlaces_either_order;
///
/// // Either argument order works.
/// assert!(interlaces_either_order(&[0, 1], &[-1, 0, 1]));
/// assert!(interlaces_either_order(&[-1, 0, 1], &[0, 1]));
/// // Same degree, roots {-1} and {0}.
/// assert!(interlaces_either_order(&[0, 1], &[1, 1]));
/// // Zero interlaces real-rooted polynomials only.
/// assert!(interlaces_either_order(&[0], &[2, 3, 1]));
/// assert!(!interlaces_either_order(&[], &[1, 0, 1]));
/// // Degree difference two.
/// assert!(!interlaces_either_order(&[1], &[-1, 0, 1]));
/// ```
pub fn interlaces_either_order(f: &[i64], g: &[i64]) -> bool {
    let f_big: Vec<BigInt> = f.iter().map(|&c| BigInt::from(c)).collect();
    let g_big: Vec<BigInt> = g.iter().map(|&c| BigInt::from(c)).collect();
    interlaces_either_order_bigint_coeffs(&f_big, &g_big)
}

/// Undirected weak interlacing for `BigInt` coefficients.
///
/// Same convention as [`interlaces_either_order`].
pub fn interlaces_either_order_bigint_coeffs(f: &[BigInt], g: &[BigInt]) -> bool {
    let holds =
        |p: &[BigInt], q: &[BigInt]| check_weak_interlacing_bigint_coeffs(p, q) == Some(true);
    match (vec_poly::degree(f), vec_poly::degree(g)) {
        (None, None) => true,
        (None, Some(_)) => is_real_rooted_bigint_coeffs(g),
        (Some(_), None) => is_real_rooted_bigint_coeffs(f),
        (Some(df), Some(dg)) if dg == df + 1 => holds(f, g),
        (Some(df), Some(dg)) if df == dg + 1 => holds(g, f),
        (Some(df), Some(dg)) if df == dg => holds(f, g) || holds(g, f),
        _ => false,
    }
}

/// Cauchy bound for root radius: all roots of `coeffs` lie in |z| < bound.
///
/// For f(t) = a_n t^n + ... + a_0, the Cauchy bound is
/// 1 + max(|a_0|, ..., |a_{n-1}|) / |a_n|.
/// Returns 2 for constant or zero polynomials.
fn cauchy_root_bound_bigint(coeffs: &[BigInt]) -> BigInt {
    let deg = match vec_poly::degree(coeffs) {
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

fn has_nonzero_one_signed_trimmed_coefficients_bigint(coeffs: &[BigInt]) -> bool {
    let Some(degree) = vec_poly::degree(coeffs) else {
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

/// Alias for [`check_interlacing`] — Bézout matrix method (the default).
pub fn check_interlacing_bezout(p: &[i64], q: &[i64]) -> Option<bool> {
    check_interlacing(p, q)
}

/// Alias for [`check_weak_interlacing`] — Bézout matrix with GCD factoring.
pub fn check_weak_interlacing_bezout(p: &[i64], q: &[i64]) -> Option<bool> {
    check_weak_interlacing(p, q)
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::real_rootedness::{
        real_rooted::is_real_rooted, sturm_fallback::is_real_rooted_sturm,
    };

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
    fn test_a046741_degrees_seven_eight_interlacing_regression() {
        let degree_seven = [1, 19, 137, 473, 815, 667, 223, 21];
        let degree_eight = [1, 22, 191, 838, 1982, 2504, 1577, 424, 34];
        assert_eq!(check_interlacing(&degree_seven, &degree_eight), Some(true));
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

    // -- Undirected interlacing --

    /// Polynomial `sign * prod (t - r)` with the given integer roots.
    fn from_roots(roots: &[i64], scale: i64) -> Vec<i64> {
        roots
            .iter()
            .fold(vec![scale], |acc, &r| vec_poly::mul(&acc, &[-r, 1]))
    }

    /// Multisets of size `len` from `values` (nondecreasing sequences).
    fn multisets(values: &[i64], len: usize) -> Vec<Vec<i64>> {
        if len == 0 {
            return vec![vec![]];
        }
        let mut out = Vec::new();
        for (i, &v) in values.iter().enumerate() {
            for mut rest in multisets(&values[i..], len - 1) {
                rest.insert(0, v);
                out.push(rest);
            }
        }
        out
    }

    /// Root-based oracle for the directed relation `p ≪ q`, with roots given
    /// in any order.
    fn directed_oracle(p_roots: &[i64], q_roots: &[i64]) -> Option<bool> {
        let mut a = p_roots.to_vec();
        let mut b = q_roots.to_vec();
        a.sort_unstable_by(|x, y| y.cmp(x));
        b.sort_unstable_by(|x, y| y.cmp(x));
        if b.len() != a.len() && b.len() != a.len() + 1 {
            return None;
        }
        // Interleave b_1, a_1, b_2, a_2, ... and require it to be nonincreasing.
        let mut merged = Vec::new();
        for i in 0..b.len() {
            merged.push(b[i]);
            if i < a.len() {
                merged.push(a[i]);
            }
        }
        Some(merged.windows(2).all(|w| w[0] >= w[1]))
    }

    enum TestPoly {
        Zero,
        Rooted(Vec<i64>, Vec<i64>),
        NotRealRooted(Vec<i64>),
    }

    impl TestPoly {
        fn coeffs(&self) -> &[i64] {
            match self {
                TestPoly::Zero => &[0, 0],
                TestPoly::Rooted(_, c) | TestPoly::NotRealRooted(c) => c,
            }
        }
    }

    fn interlacing_test_family() -> Vec<TestPoly> {
        let values = [-2, -1, 0, 1, 2];
        let mut family = vec![TestPoly::Zero];
        for len in 0..=3 {
            for roots in multisets(&values, len) {
                for scale in [1, -3] {
                    let coeffs = from_roots(&roots, scale);
                    family.push(TestPoly::Rooted(roots.clone(), coeffs));
                }
            }
        }
        for coeffs in [
            vec![1, 0, 1],      // t^2 + 1
            vec![-1, 1, -1],    // -(t^2 - t + 1)
            vec![-1, 1, -1, 1], // (t^2 + 1)(t - 1)
            vec![0, 1, 0, 1],   // t (t^2 + 1)
            vec![-2, 0, 0, 1],  // t^3 - 2
        ] {
            family.push(TestPoly::NotRealRooted(coeffs));
        }
        family
    }

    #[test]
    fn zero_polynomial_conventions() {
        assert!(is_real_rooted(&[0]));
        assert!(is_real_rooted(&[]));
        assert_eq!(check_weak_interlacing(&[], &[]), Some(true));
        assert_eq!(check_weak_interlacing(&[0], &[1, 0, 1]), Some(false));
        assert_eq!(check_weak_interlacing(&[1], &[0]), Some(false));
        assert!(interlaces_either_order(&[], &[0, 0]));
        assert!(interlaces_either_order(&[5], &[0]));
        assert!(interlaces_either_order(&[0], &[1, 2, 1]));
        assert!(!interlaces_either_order(&[1, 0, 1], &[0]));
    }

    #[test]
    fn weak_interlacing_matches_root_oracle_exhaustively() {
        let family = interlacing_test_family();
        let mut checked_true = 0usize;
        for f in &family {
            for g in &family {
                let (fc, gc) = (f.coeffs(), g.coeffs());
                let directed = check_weak_interlacing(fc, gc);
                // Directed relation against the oracle.
                let expected_directed = match (f, g) {
                    (TestPoly::Zero, _) => Some(is_real_rooted(gc)),
                    (_, TestPoly::Zero) => Some(false),
                    (TestPoly::Rooted(a, _), TestPoly::Rooted(b, _)) => directed_oracle(a, b),
                    _ => {
                        let (df, dg) = (vec_poly::degree(fc), vec_poly::degree(gc));
                        let (df, dg) = (df.unwrap(), dg.unwrap());
                        (dg == df || dg == df + 1).then_some(false)
                    }
                };
                assert_eq!(
                    directed, expected_directed,
                    "directed f = {fc:?}, g = {gc:?}"
                );

                // Undirected relation against the oracle.
                let expected = match (f, g) {
                    (TestPoly::Zero, TestPoly::Zero) => true,
                    (TestPoly::Zero, other) | (other, TestPoly::Zero) => {
                        !matches!(other, TestPoly::NotRealRooted(_))
                    }
                    (TestPoly::Rooted(a, _), TestPoly::Rooted(b, _)) => {
                        directed_oracle(a, b) == Some(true) || directed_oracle(b, a) == Some(true)
                    }
                    _ => false,
                };
                let actual = interlaces_either_order(fc, gc);
                assert_eq!(actual, expected, "undirected f = {fc:?}, g = {gc:?}");
                assert_eq!(actual, interlaces_either_order(gc, fc), "symmetry");
                // Consistency with the directed checker.
                assert_eq!(
                    actual,
                    check_weak_interlacing(fc, gc) == Some(true)
                        || check_weak_interlacing(gc, fc) == Some(true)
                );
                let f_big: Vec<BigInt> = fc.iter().map(|&c| BigInt::from(c)).collect();
                let g_big: Vec<BigInt> = gc.iter().map(|&c| BigInt::from(c)).collect();
                assert_eq!(
                    interlaces_either_order_bigint_coeffs(&f_big, &g_big),
                    actual
                );
                checked_true += usize::from(actual);
            }
        }
        assert!(
            checked_true > 1000,
            "family should contain many interlacing pairs"
        );
    }

    #[test]
    fn interlaces_either_order_bigint_large_coefficients() {
        // Roots c - 1, c + 1 and c with c = 10^20.
        let c = BigInt::from(10).pow(20);
        let f = vec![-&c, BigInt::from(1)];
        let g = vec![
            (&c - 1u32) * (&c + 1u32),
            -BigInt::from(2) * &c,
            BigInt::from(1),
        ];
        assert!(interlaces_either_order_bigint_coeffs(&f, &g));
        assert!(interlaces_either_order_bigint_coeffs(&g, &f));
        // Shift f's root outside [c - 1, c + 1].
        let f_far = vec![-(&c + 2u32), BigInt::from(1)];
        assert!(!interlaces_either_order_bigint_coeffs(&f_far, &g));
    }
}
