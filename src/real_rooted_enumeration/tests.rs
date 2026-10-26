//! Tests for the real-rooted enumeration and its kernel.
//!
//! The oracles below are independent of the kernel: they use the Sturm-chain
//! real-rootedness backend and coefficient boxes derived from Newton's
//! identities and Maclaurin's inequality, and never use Lemma 2.5.

use super::*;
use crate::numeric::binomial_bigint;
use crate::real_rootedness::{discriminant_bigint_coeffs, is_real_rooted_sturm_bigint_coeffs};

fn big(values: &[i64]) -> Vec<BigInt> {
    values.iter().map(|&v| BigInt::from(v)).collect()
}

fn range(lower: i64, upper: i64) -> Option<InclusiveIntegerRange> {
    Some(InclusiveIntegerRange {
        lower: BigInt::from(lower),
        upper: BigInt::from(upper),
    })
}

fn kernel(q: &[i64]) -> Option<InclusiveIntegerRange> {
    real_rooted_antiderivative_integer_constants_bigint(&big(q)).unwrap()
}

fn eval_rational_rational(p: &[BigRational], x: &BigRational) -> BigRational {
    let mut acc = BigRational::zero();
    for c in p.iter().rev() {
        acc = acc * x + c;
    }
    acc
}

fn ceil_div_rational(value: &BigRational, modulus: &BigInt) -> BigInt {
    (value / BigRational::from_integer(modulus.clone()))
        .ceil()
        .to_integer()
}

fn floor_div_rational(value: &BigRational, modulus: &BigInt) -> BigInt {
    (value / BigRational::from_integer(modulus.clone()))
        .floor()
        .to_integer()
}

fn mul(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let mut out = vec![BigInt::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    out
}

fn derivative_k(p: &[BigInt], k: usize) -> Vec<BigInt> {
    let mut p = p.to_vec();
    for _ in 0..k {
        p = derivative(&p);
    }
    p
}

// ---------------------------------------------------------------------------
// Antiderivatives
// ---------------------------------------------------------------------------

#[test]
fn exact_antiderivative_divides_and_reports_failures() {
    assert_eq!(
        antiderivative_exact_bigint(&big(&[2, 4, 3])).unwrap(),
        big(&[0, 2, 2, 1])
    );
    assert_eq!(antiderivative_exact_bigint(&big(&[])).unwrap(), big(&[0]));
    assert_eq!(
        antiderivative_exact_bigint(&big(&[0, 0])).unwrap(),
        big(&[0])
    );
    assert_eq!(
        antiderivative_exact_bigint(&big(&[5, 0, 0])).unwrap(),
        big(&[0, 5])
    );
    let error = antiderivative_exact_bigint(&big(&[1, 2, 4])).unwrap_err();
    assert_eq!(error.exponent, 3);
    assert_eq!(error.coefficient, BigInt::from(4));
    assert!(error.to_string().contains("not divisible by 3"));
}

#[test]
fn rational_antiderivative_matches_exact_one() {
    let q = big(&[-3, 4, 9, -8]);
    let rational = antiderivative_rational_bigint(&q);
    let exact = antiderivative_exact_bigint(&q).unwrap();
    assert_eq!(rational.len(), exact.len());
    for (r, e) in rational.iter().zip(&exact) {
        assert_eq!(r, &BigRational::from_integer(e.clone()));
    }
    assert_eq!(
        antiderivative_rational_bigint(&[]),
        vec![BigRational::zero()]
    );
}

// ---------------------------------------------------------------------------
// Kernel: worked examples
// ---------------------------------------------------------------------------

#[test]
fn kernel_cubic_example_with_integral_repeated_root_endpoints() {
    // P = x^3 - 3x: endpoints -2, 2 give (x-1)^2 (x+2) and (x+1)^2 (x-2).
    assert_eq!(kernel(&[-3, 0, 3]), range(-2, 2));
    for c in [-2, 2] {
        let shifted = big(&[c, -3, 0, 1]);
        assert!(discriminant_bigint_coeffs(&shifted).is_zero());
        assert!(is_real_rooted_sturm_bigint_coeffs(&shifted));
    }
}

#[test]
fn kernel_nonintegral_endpoints() {
    // P = x^3/3 - x: interval [-2/3, 2/3].
    assert_eq!(kernel(&[-1, 0, 1]), range(0, 0));
    // q = 3x^2 - 5, P = x^3 - 5x, critical values
    // ±(10/3) sqrt(5/3) ≈ ±4.30, so the integers are -4..=4.
    assert_eq!(kernel(&[-5, 0, 3]), range(-4, 4));
    // Asymmetric: q = 3x^2 - 6x, P = x^3 - 3x^2, interval [0, 4].
    assert_eq!(kernel(&[0, -6, 3]), range(0, 4));
}

#[test]
fn kernel_negative_leading_coefficient_negates_interval() {
    assert_eq!(kernel(&[3, 0, -3]), range(-2, 2));
    assert_eq!(kernel(&[0, 6, -3]), range(-4, 0));
    assert_eq!(kernel(&[1, 0, -1]), range(0, 0));
}

#[test]
fn kernel_repeated_roots_force_a_single_constant() {
    // q = x^2: P = x^3/3, only C = 0.
    assert_eq!(kernel(&[0, 0, 1]), range(0, 0));
    // q = 4 (x-1)^3: P = (x-1)^4 - 1, only C = 1.
    assert_eq!(kernel(&[-4, 12, -12, 4]), range(1, 1));
    // q = (x-1)^2: P = x^3/3 - x^2 + x, forced C = -1/3 is not an integer.
    assert_eq!(kernel(&[1, -2, 1]), None);
    // q = (x-1)^2 (x+2): the double root at 1 forces C = -P(1) = -3/4.
    assert_eq!(kernel(&[2, -3, 0, 1]), None);
    // q = 12 (x-1)^2 (x+2): P(1) = 9, forced C = -9, and P - 9 = 3(x-1)^3 (x+3).
    assert_eq!(kernel(&[24, -36, 0, 12]), range(-9, -9));
}

#[test]
fn kernel_several_distinct_multiple_roots() {
    // q = 30 x^2 (x-1)^2: P = 6x^5 - 15x^4 + 10x^3, P(0) = 0, P(1) = 1:
    // the two forced values disagree.
    assert_eq!(kernel(&[0, 0, 30, -60, 30]), None);
    // q = 15 (x^2 - 1)^2, P = 3x^5 - 10x^3 + 15x, P(1) = 8, P(-1) = -8:
    // the forced values disagree again.
    assert_eq!(kernel(&[15, 0, -30, 0, 15]), None);
    // q = x^2 (x^2 - 1) has a double root at 0 and simple roots ±1:
    // P = x^5/5 - x^3/3, forced C = 0, and P = x^3 (3x^2 - 5)/15 is real-rooted.
    assert_eq!(kernel(&[0, 0, -1, 0, 1]), range(0, 0));
}

#[test]
fn kernel_coincident_critical_values() {
    // P = x^4 - 2x^2: -P(-1) = -P(1) = 1 (two minima), -P(0) = 0.
    assert_eq!(kernel(&[0, -4, 0, 4]), range(0, 1));
    // q = 5x^4 - 15x^2 + 4 has simple roots and a symmetric pair of equal
    // critical values; compare against the brute-force reference.
    let q = big(&[4, 0, -15, 0, 5]);
    assert_eq!(
        real_rooted_antiderivative_integer_constants_bigint(&q).unwrap(),
        reference_constants(&q, &BigInt::one(), 3)
    );
}

#[test]
fn kernel_empty_interval_for_non_real_rooted_derivative() {
    assert_eq!(kernel(&[1, 0, 1]), None);
    assert_eq!(kernel(&[1, 0, 0, 1]), None);
    // Real-rooted q with an empty integer interval: q = 3x^2 - 1,
    // P = x^3 - x, critical values ±2/(3 sqrt 3) ≈ ±0.385: only C = 0.
    assert_eq!(kernel(&[-1, 0, 3]), range(0, 0));
    // q = 12x^2 - 1: P = 4x^3 - x has critical values ±1/(3 sqrt 3), and
    // with modulus 7 the only multiple is 0.
    let q = big(&[-1, 0, 12]);
    assert_eq!(
        real_rooted_antiderivative_constant_multiples_bigint(&q, &BigInt::from(7)).unwrap(),
        range(0, 0)
    );
    // q = 6x^2 - 6x: P = 2x^3 - 3x^2 has admissible interval [0, 1];
    // modulus 2 keeps only C = 0.
    let q = big(&[0, -6, 6]);
    assert_eq!(
        real_rooted_antiderivative_constant_multiples_bigint(&q, &BigInt::from(2)).unwrap(),
        range(0, 0)
    );
    // q = 3x^2 - 3x: P = x^3 - 3x^2/2, interval [0, 1/2].
    assert_eq!(kernel(&[0, -3, 3]), range(0, 0));
    // A real-rooted q whose admissible interval contains no integer:
    // q = 24x^2 - 24x + 5, P = 8x^3 - 12x^2 + 5x = 8y^3 - y + 1/2 with
    // y = x - 1/2, so P + C is real-rooted iff C + 1/2 lies in
    // [-1/(3 sqrt 6), 1/(3 sqrt 6)] ≈ [-0.136, 0.136].
    assert_eq!(kernel(&[5, -24, 24]), None);
    assert_eq!(kernel(&[-5, 24, -24]), None);
    // P = x^3 - 6x^2 + 11x = (x-1)(x-2)(x-3) + 6: only C = -6.
    assert_eq!(kernel(&[11, -12, 3]), range(-6, -6));
}

#[test]
fn kernel_coefficients_beyond_i64() {
    let scale = BigInt::from(10).pow(30);
    let q: Vec<BigInt> = big(&[-3, 0, 3]).iter().map(|c| c * &scale).collect();
    let r = real_rooted_antiderivative_integer_constants_bigint(&q)
        .unwrap()
        .unwrap();
    assert_eq!(r.lower, -BigInt::from(2) * &scale);
    assert_eq!(r.upper, BigInt::from(2) * &scale);

    // q = 3x^2 - 3B^2 with B = 10^20: interval [-2B^3, 2B^3].
    let b = BigInt::from(10).pow(20);
    let q = vec![-BigInt::from(3) * &b * &b, BigInt::zero(), BigInt::from(3)];
    let r = real_rooted_antiderivative_integer_constants_bigint(&q)
        .unwrap()
        .unwrap();
    let b3 = &b * &b * &b;
    assert_eq!(r.lower, -BigInt::from(2) * &b3);
    assert_eq!(r.upper, BigInt::from(2) * &b3);
    assert_eq!(r.count(), BigInt::from(4) * &b3 + 1);
}

#[test]
fn kernel_rejects_unsupported_input() {
    assert_eq!(
        real_rooted_antiderivative_integer_constants_bigint(&big(&[0, 0])),
        Err(AntiderivativeIntervalError::ZeroPolynomial)
    );
    assert_eq!(
        real_rooted_antiderivative_integer_constants_bigint(&big(&[1, 2, 0])),
        Err(AntiderivativeIntervalError::UnboundedDegree { degree: 1 })
    );
    assert_eq!(
        real_rooted_antiderivative_integer_constants_bigint(&big(&[5])),
        Err(AntiderivativeIntervalError::UnboundedDegree { degree: 0 })
    );
    assert_eq!(
        real_rooted_antiderivative_constant_multiples_bigint(&big(&[-3, 0, 3]), &BigInt::zero()),
        Err(AntiderivativeIntervalError::NonPositiveModulus)
    );
}

#[test]
fn inclusive_range_helpers() {
    let r = range(-2, 1).unwrap();
    assert_eq!(r.count(), BigInt::from(4));
    assert!(r.contains(&BigInt::from(-2)) && !r.contains(&BigInt::from(2)));
    assert_eq!(r.iter().collect::<Vec<_>>(), big(&[-2, -1, 0, 1]));
    assert_eq!(range(3, 3).unwrap().iter().count(), 1);
}

// ---------------------------------------------------------------------------
// Kernel: comparison with a brute-force reference
// ---------------------------------------------------------------------------

/// Independent reference: test every multiple of `modulus` in `[lo, hi]`,
/// which must contain all critical values `-P(µ)` (the admissible set lies
/// in their hull), and check that the admissible multiples are consecutive.
fn reference_multiples(
    q: &[BigInt],
    modulus: &BigInt,
    lo: &BigRational,
    hi: &BigRational,
) -> Option<InclusiveIntegerRange> {
    let p = antiderivative_rational_bigint(q);
    let denominator = p.iter().fold(BigInt::one(), |acc, c| acc.lcm(c.denom()));
    let scaled: Vec<BigInt> = p
        .iter()
        .map(|c| (c * BigRational::from_integer(denominator.clone())).to_integer())
        .collect();
    let mut found: Vec<BigInt> = Vec::new();
    let mut k = ceil_div_rational(lo, modulus);
    let k_last = floor_div_rational(hi, modulus);
    while k <= k_last {
        let mut shifted = scaled.clone();
        shifted[0] += &denominator * &k * modulus;
        if is_real_rooted_sturm_bigint_coeffs(&shifted) {
            found.push(k.clone());
        }
        k += 1;
    }
    let (first, last) = (found.first()?.clone(), found.last()?.clone());
    assert_eq!(
        BigInt::from(found.len()),
        &last - &first + 1,
        "admissible set of {q:?} is not an interval"
    );
    Some(InclusiveIntegerRange {
        lower: first,
        upper: last,
    })
}

/// Reference with the coarse bound `|P(µ)| <= Σ |P_k| R^k`, where
/// `root_bound = R` bounds the absolute values of the roots of `q`.
fn reference_constants(
    q: &[BigInt],
    modulus: &BigInt,
    root_bound: i64,
) -> Option<InclusiveIntegerRange> {
    let p = antiderivative_rational_bigint(q);
    let r = BigRational::from_integer(BigInt::from(root_bound));
    let abs: Vec<BigRational> = p.iter().map(Signed::abs).collect();
    let bound = eval_rational_rational(&abs, &r);
    reference_multiples(q, modulus, &-bound.clone(), &bound)
}

/// Reference for `q` with known rational roots: scan the exact hull of the
/// critical values.
fn reference_from_roots(
    q: &[BigInt],
    modulus: &BigInt,
    roots: &[BigRational],
) -> Option<InclusiveIntegerRange> {
    let p = antiderivative_rational_bigint(q);
    let values: Vec<BigRational> = roots
        .iter()
        .map(|r| -eval_rational_rational(&p, r))
        .collect();
    let lo = values.iter().min().unwrap();
    let hi = values.iter().max().unwrap();
    reference_multiples(q, modulus, lo, hi)
}

/// Products of linear factors `a x - b` with small integers, including
/// repeated factors, and their negatives.
fn real_rooted_fixture_derivatives() -> Vec<(Vec<BigInt>, Vec<BigRational>)> {
    let factors: [(i64, i64); 7] = [(1, 0), (1, 1), (1, -1), (1, 2), (2, 1), (3, -2), (1, -3)];
    let mut out = Vec::new();
    let mut push = |indices: &[usize], lead: i64| {
        let mut q = big(&[lead]);
        let mut roots = Vec::new();
        for &i in indices {
            let (a, b) = factors[i];
            q = mul(&q, &big(&[-b, a]));
            roots.push(BigRational::new(BigInt::from(b), BigInt::from(a)));
        }
        out.push((q, roots));
    };
    for i in 0..factors.len() {
        for j in i..factors.len() {
            push(&[i, j], 1);
            push(&[i, j], -2);
            for k in j..factors.len() {
                push(&[i, j, k], 1);
                if (i + j + k) % 3 == 0 {
                    push(&[i, j, k], -3);
                }
            }
        }
    }
    // A few degree-four products, including double and triple roots.
    for indices in [
        [0, 1, 2, 3],
        [0, 0, 1, 2],
        [1, 1, 1, 2],
        [2, 2, 3, 3],
        [4, 5, 5, 6],
        [0, 0, 0, 0],
    ] {
        push(&indices, 1);
        push(&indices, 4);
    }
    out
}

#[test]
fn kernel_agrees_with_reference_on_real_rooted_fixtures() {
    for (q, roots) in real_rooted_fixture_derivatives() {
        for modulus in [1, 2, 6] {
            let modulus = BigInt::from(modulus);
            let fast = real_rooted_antiderivative_constant_multiples_bigint(&q, &modulus).unwrap();
            let slow = reference_from_roots(&q, &modulus, &roots);
            assert_eq!(fast, slow, "q = {q:?}, modulus = {modulus}");
        }
    }
}

#[test]
fn kernel_agrees_with_reference_on_small_coefficient_grid() {
    // Every cubic and quadratic with coefficients in -2..=2 (nonzero lead).
    let mut checked = 0;
    for degree in 2..=3usize {
        let total = 5usize.pow(degree as u32 + 1);
        for mut code in 0..total {
            let mut q = Vec::with_capacity(degree + 1);
            for _ in 0..=degree {
                q.push(BigInt::from((code % 5) as i64 - 2));
                code /= 5;
            }
            if q[degree].is_zero() {
                continue;
            }
            // Cauchy bound 1 + max |q_k / q_d| <= 3.
            let fast = real_rooted_antiderivative_integer_constants_bigint(&q).unwrap();
            let slow = if is_real_rooted_sturm_bigint_coeffs(&q) {
                reference_constants(&q, &BigInt::one(), 3)
            } else {
                None
            };
            assert_eq!(fast, slow, "q = {q:?}");
            checked += 1;
        }
    }
    assert_eq!(checked, 100 + 500);
}

// ---------------------------------------------------------------------------
// Enumerator: brute-force oracle
// ---------------------------------------------------------------------------

/// Coefficient box from Proposition 2.8: for real roots `λ` of
/// `p = a_0 Π (x - λ_j)`, `Σ λ_j^2 = N / a_0^2` with `N = a_1^2 - 2 a_0 a_2`,
/// and Maclaurin's inequality together with the power-mean inequality give
/// `|a_i| <= |a_0| C(n,i) (N / (n a_0^2))^(i/2)`.  Returns `None` if `N < 0`.
fn coefficient_box(n: usize, prefix: &[BigInt]) -> Option<Vec<BigInt>> {
    let (a0, a1, a2) = (&prefix[0], &prefix[1], &prefix[2]);
    let norm = a1 * a1 - BigInt::from(2) * a0 * a2;
    if norm.is_negative() {
        return None;
    }
    let a0_sq = a0 * a0;
    Some(
        (0..=n)
            .map(|i| {
                let c = binomial_bigint(n, i);
                let numerator = &a0_sq * &c * &c * norm.pow(i as u32);
                let denominator = BigInt::from(n).pow(i as u32) * a0_sq.pow(i as u32);
                (numerator / denominator).sqrt()
            })
            .collect(),
    )
}

/// Ascending coefficients of `a_0 x^n + ... + a_i x^(n-i)` divided by
/// `x^(n-i)` after differentiating `n - i` times, i.e. `p_i`.
fn derivative_of_partial(n: usize, descending: &[BigInt]) -> Vec<BigInt> {
    let mut ascending = vec![BigInt::zero(); n + 1];
    for (j, a) in descending.iter().enumerate() {
        ascending[n - j] = a.clone();
    }
    derivative_k(&ascending, n + 1 - descending.len())
}

/// Brute force: every coefficient vector in the box, pruned only by the
/// necessary condition that each derivative `p_i` (which depends on
/// `a_0..a_i` only) is real-rooted, then filtered by the Sturm test.
fn brute_force(n: usize, prefix: &[BigInt]) -> Vec<Vec<BigInt>> {
    let mut out = Vec::new();
    let Some(bounds) = coefficient_box(n, prefix) else {
        return out;
    };
    for (i, a) in prefix.iter().enumerate() {
        if a.abs() > bounds[i] {
            return out;
        }
    }
    fn recurse(
        n: usize,
        bounds: &[BigInt],
        descending: &mut Vec<BigInt>,
        out: &mut Vec<Vec<BigInt>>,
    ) {
        if !is_real_rooted_sturm_bigint_coeffs(&derivative_of_partial(n, descending)) {
            return;
        }
        if descending.len() == n + 1 {
            out.push(descending.iter().rev().cloned().collect());
            return;
        }
        let b = &bounds[descending.len()];
        let mut a = -b.clone();
        while &a <= b {
            descending.push(a.clone());
            recurse(n, bounds, descending, out);
            descending.pop();
            a += 1;
        }
    }
    let mut descending = prefix.to_vec();
    recurse(n, &bounds, &mut descending, &mut out);
    out
}

fn enumerate(n: usize, prefix: &[i64]) -> RealRootedEnumerationResult {
    enumerate_real_rooted_integer_polynomials_bigint(
        n,
        &big(prefix),
        &RealRootedEnumerationOptions::default(),
    )
    .unwrap()
}

/// Check the documented properties of every emitted polynomial.
fn check_outputs(n: usize, prefix: &[BigInt], result: &RealRootedEnumerationResult) {
    for p in &result.polynomials {
        assert_eq!(p.len(), n + 1, "degree of {p:?}");
        for (j, a) in prefix.iter().enumerate() {
            assert_eq!(&p[n - j], a, "prefix of {p:?}");
        }
        assert!(is_real_rooted_sturm_bigint_coeffs(p), "{p:?}");
        for k in 0..=n {
            let tower = derivative_k(p, k);
            assert!(
                is_real_rooted_sturm_bigint_coeffs(&tower),
                "{p:?} level {k}"
            );
            // Constant term of p_(n-k) is k! a_(n-k).
            assert_eq!(tower[0], &p[k] * factorial_bigint(k));
        }
    }
    for pair in result.polynomials.windows(2) {
        let order = pair[0].iter().rev().cmp(pair[1].iter().rev());
        assert_eq!(order, Ordering::Less, "lexicographic order");
    }
}

fn assert_matches_oracle(n: usize, prefix: &[i64]) -> usize {
    let prefix_big = big(prefix);
    let result = enumerate(n, prefix);
    assert_eq!(result.status, RealRootedEnumerationStatus::Complete);
    check_outputs(n, &prefix_big, &result);
    let oracle = brute_force(n, &prefix_big);
    assert_eq!(result.polynomials, oracle, "n = {n}, prefix = {prefix:?}");
    assert_eq!(result.stats.results_emitted, oracle.len());
    assert_eq!(
        *result.stats.polynomials_by_level.last().unwrap(),
        oracle.len() as u64
    );
    oracle.len()
}

#[test]
fn enumeration_matches_brute_force_degree_three() {
    // x^3 - 3x + c: c in -2..=2.
    assert_eq!(assert_matches_oracle(3, &[1, 0, -3]), 5);
    assert_eq!(assert_matches_oracle(3, &[1, 0, -1]), 1);
    // x^3 + c is real-rooted only for c = 0.
    assert_eq!(assert_matches_oracle(3, &[1, 0, 0]), 1);
    // Newton bound violated: a_1^2 < 2 a_0 a_2.
    assert_eq!(assert_matches_oracle(3, &[1, 0, 1]), 0);
    // (x - 2)^3 = x^3 - 6x^2 + 12x - 8 is the only completion.
    assert_eq!(assert_matches_oracle(3, &[1, -6, 12]), 1);
    for prefix in [
        [2, 1, -3],
        [3, 1, -2],
        [-1, 0, 3],
        [-2, 3, 1],
        [1, 1, -7],
        [4, -4, -9],
    ] {
        assert_matches_oracle(3, &prefix);
    }
}

#[test]
fn enumeration_matches_brute_force_degree_four_and_five() {
    let mut total = 0;
    for prefix in [
        [1, 0, -4],
        [1, 0, -6],
        [1, -4, 6],
        [2, 0, -5],
        [-1, 0, 5],
        [1, 2, -3],
        [3, -1, -4],
        [1, 1, 1],
    ] {
        total += assert_matches_oracle(4, &prefix);
    }
    for prefix in [[1, 0, -5], [1, -5, 10], [1, 0, -8], [-2, 1, 4]] {
        total += assert_matches_oracle(5, &prefix);
    }
    assert!(total > 0);
}

#[test]
#[cfg_attr(
    debug_assertions,
    ignore = "slow oracle in debug builds; run with --release"
)]
fn enumeration_matches_brute_force_degree_six() {
    assert_matches_oracle(6, &[1, 0, -6]);
}

#[test]
fn enumeration_with_longer_prefixes() {
    // t = 4.
    assert_matches_oracle(5, &[1, 0, -5, 0]);
    assert_matches_oracle(4, &[1, 0, -5, 2]);
    // t = n + 1: the input itself, if real-rooted.
    let r = enumerate(4, &[1, 0, -5, 0, 4]); // (x^2 - 1)(x^2 - 4)
    assert_eq!(r.status, RealRootedEnumerationStatus::Complete);
    assert_eq!(r.polynomials, vec![big(&[4, 0, -5, 0, 1])]);
    let r = enumerate(4, &[1, -4, 6, -4, 1]); // (x - 1)^4
    assert_eq!(r.polynomials, vec![big(&[1, -4, 6, -4, 1])]);
    let r = enumerate(4, &[1, 0, 1, 0, 0]); // x^2 (x^2 + 1)
    assert_eq!(r.status, RealRootedEnumerationStatus::Complete);
    assert!(r.polynomials.is_empty());
    assert_eq!(r.stats.derivative_levels, 0);
}

#[test]
fn enumeration_includes_reducible_and_repeated_root_polynomials() {
    let r = enumerate(4, &[1, 0, -4]);
    // (x^2 - 2)^2 has two double roots; x^2 (x - 2)(x + 2) is reducible.
    assert!(r.polynomials.contains(&big(&[4, 0, -4, 0, 1])));
    assert!(r.polynomials.contains(&big(&[0, 0, -4, 0, 1])));
    // (x - a)^n for several a: prefix (1, -n a, C(n,2) a^2).
    for a in [-2i64, -1, 0, 1, 3] {
        for n in [3usize, 4, 5] {
            let ni = n as i64;
            let prefix = [1, -ni * a, ni * (ni - 1) / 2 * a * a];
            let r = enumerate(n, &prefix);
            let power: Vec<BigInt> = (0..=n)
                .map(|k| binomial_bigint(n, k) * BigInt::from(-a).pow((n - k) as u32))
                .collect();
            // The prefix fixes the root mean and variance (zero), so the
            // power is the only completion.
            assert_eq!(r.polynomials, vec![power], "a = {a}, n = {n}");
        }
    }
}

#[test]
fn seidel_trace_polynomial_counts() {
    // Monic, a_1 = 0, a_2 = -C(n,2): the trace polynomials T_n of the paper
    // (|T_3| = 5: x^3 - 3x + c with |c| <= 2).  The counts for n <= 5 are
    // checked against the brute-force oracle.
    assert_eq!(assert_matches_oracle(3, &[1, 0, -3]), 5);
    assert_eq!(assert_matches_oracle(4, &[1, 0, -6]), 64);
}

#[test]
#[cfg_attr(
    debug_assertions,
    ignore = "slow oracle in debug builds; run with --release"
)]
fn seidel_trace_polynomial_count_degree_five() {
    assert_eq!(assert_matches_oracle(5, &[1, 0, -10]), 3116);
}

#[test]
#[cfg_attr(debug_assertions, ignore = "slow in debug builds; run with --release")]
fn seidel_trace_polynomials_degree_six_regression() {
    // Regression value (no oracle at this size): 110190 examined nodes.
    let r = enumerate(6, &[1, 0, -15]);
    assert_eq!(r.status, RealRootedEnumerationStatus::Complete);
    assert_eq!(r.polynomials.len(), 588_724);
    assert_eq!(
        r.stats.polynomials_by_level,
        vec![1, 81, 4288, 105_820, 588_724]
    );
    // Spot-check the emitted polynomials (every 997th).
    let sample = RealRootedEnumerationResult {
        polynomials: r.polynomials.iter().step_by(997).cloned().collect(),
        ..r.clone()
    };
    check_outputs(6, &big(&[1, 0, -15]), &sample);
}

#[test]
fn invalid_input_is_rejected() {
    let opts = RealRootedEnumerationOptions::default();
    let e = |n: usize, prefix: &[i64]| {
        enumerate_real_rooted_integer_polynomials_bigint(n, &big(prefix), &opts).unwrap_err()
    };
    assert_eq!(
        e(1, &[1, 0, 0]),
        RealRootedEnumerationError::DegreeTooSmall { degree: 1 }
    );
    assert_eq!(
        e(4, &[1, 0]),
        RealRootedEnumerationError::PrefixLength {
            prefix_len: 2,
            degree: 4
        }
    );
    assert_eq!(
        e(2, &[1, 0, 0, 0]),
        RealRootedEnumerationError::PrefixLength {
            prefix_len: 4,
            degree: 2
        }
    );
    assert_eq!(
        e(4, &[0, 1, 1]),
        RealRootedEnumerationError::ZeroLeadingCoefficient
    );
    assert!(e(4, &[0, 1, 1]).to_string().contains("nonzero"));
}

#[test]
fn degree_two_with_full_prefix() {
    let r = enumerate(2, &[1, 0, -1]);
    assert_eq!(r.polynomials, vec![big(&[-1, 0, 1])]);
    let r = enumerate(2, &[1, 0, 1]);
    assert!(r.polynomials.is_empty());
    assert_eq!(r.status, RealRootedEnumerationStatus::Complete);
}

#[test]
fn coefficients_beyond_i64_in_enumeration() {
    // x^3 - 3 B^2 x + c with B = 10^7: c in [-2B^3, 2B^3] has 4*10^21 + 1
    // members; take a small result budget and check the first ones.
    let b = BigInt::from(10).pow(7);
    let prefix = vec![BigInt::one(), BigInt::zero(), -BigInt::from(3) * &b * &b];
    let options = RealRootedEnumerationOptions {
        max_results: Some(3),
        ..RealRootedEnumerationOptions::default()
    };
    let r = enumerate_real_rooted_integer_polynomials_bigint(3, &prefix, &options).unwrap();
    assert_eq!(
        r.status,
        RealRootedEnumerationStatus::ResourceLimit(RealRootedEnumerationLimit::Results)
    );
    let b3 = &b * &b * &b;
    let first: Vec<BigInt> = r.polynomials.iter().map(|p| p[0].clone()).collect();
    assert_eq!(
        first,
        vec![
            -BigInt::from(2) * &b3,
            -BigInt::from(2) * &b3 + 1,
            -BigInt::from(2) * &b3 + 2
        ]
    );
    assert!(first[0].bits() > 63);
}

// ---------------------------------------------------------------------------
// Budgets and determinism
// ---------------------------------------------------------------------------

#[test]
fn tiny_budgets_report_resource_limits() {
    let prefix = big(&[1, 0, -6]);
    let full = enumerate(4, &[1, 0, -6]);
    assert!(!full.polynomials.is_empty());
    let run = |options: RealRootedEnumerationOptions| {
        enumerate_real_rooted_integer_polynomials_bigint(4, &prefix, &options).unwrap()
    };

    let r = run(RealRootedEnumerationOptions {
        max_nodes: Some(0),
        ..Default::default()
    });
    assert_eq!(
        r.status,
        RealRootedEnumerationStatus::ResourceLimit(RealRootedEnumerationLimit::Nodes)
    );
    assert!(r.polynomials.is_empty());

    let r = run(RealRootedEnumerationOptions {
        max_results: Some(0),
        ..Default::default()
    });
    assert_eq!(
        r.status,
        RealRootedEnumerationStatus::ResourceLimit(RealRootedEnumerationLimit::Results)
    );
    assert!(r.polynomials.is_empty());

    let r = run(RealRootedEnumerationOptions {
        max_coefficient_bits: Some(2),
        ..Default::default()
    });
    assert_eq!(
        r.status,
        RealRootedEnumerationStatus::ResourceLimit(RealRootedEnumerationLimit::CoefficientBits)
    );
    assert_eq!(r.status.as_str(), "resource_limit");
    assert_eq!(
        RealRootedEnumerationLimit::CoefficientBits.as_str(),
        "max_coefficient_bits"
    );

    // A result budget equal to the full count is complete.
    let r = run(RealRootedEnumerationOptions {
        max_results: Some(full.polynomials.len()),
        ..Default::default()
    });
    assert_eq!(r.status, RealRootedEnumerationStatus::Complete);
    assert_eq!(r.polynomials, full.polynomials);
}

#[test]
fn partial_results_are_prefixes_and_budgets_are_monotone() {
    let prefix = big(&[1, 0, -6]);
    let full = enumerate(4, &[1, 0, -6]);
    assert_eq!(full.status, RealRootedEnumerationStatus::Complete);
    let nodes = full.stats.nodes_examined;
    for budget in 0..=nodes + 2 {
        let r = enumerate_real_rooted_integer_polynomials_bigint(
            4,
            &prefix,
            &RealRootedEnumerationOptions {
                max_nodes: Some(budget),
                ..Default::default()
            },
        )
        .unwrap();
        assert!(full.polynomials.starts_with(&r.polynomials));
        if budget >= nodes {
            assert_eq!(r, full);
        } else {
            assert_eq!(
                r.status,
                RealRootedEnumerationStatus::ResourceLimit(RealRootedEnumerationLimit::Nodes)
            );
        }
    }
    for budget in 0..=full.polynomials.len() {
        let r = enumerate_real_rooted_integer_polynomials_bigint(
            4,
            &prefix,
            &RealRootedEnumerationOptions {
                max_results: Some(budget),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(r.polynomials, full.polynomials[..budget]);
    }
}

#[test]
fn enumeration_is_deterministic_and_counts_are_consistent() {
    let a = enumerate(5, &[1, 0, -10]);
    let b = enumerate(5, &[1, 0, -10]);
    assert_eq!(a, b);
    let s = &a.stats;
    assert_eq!(s.derivative_levels, 3);
    assert_eq!(s.polynomials_by_level.len(), 4);
    assert_eq!(s.polynomials_by_level[0], 1);
    // Every non-leaf tower polynomial is examined exactly once.
    let non_leaf: u64 = s.polynomials_by_level[..3].iter().sum();
    assert_eq!(s.nodes_examined, non_leaf);
    assert_eq!(s.exact_root_isolations, s.nodes_examined);
    assert!(s.nodes_pruned <= s.nodes_examined);
    assert_eq!(s.final_verifications, a.polynomials.len() as u64);
}

#[test]
fn certifying_path_has_no_floating_point() {
    let source = include_str!("../real_rooted_enumeration.rs");
    for token in ["f64", "f32", "as_f", "to_f"] {
        assert!(!source.contains(token), "found {token}");
    }
}
