//! Charge generating functions of cylindric staircase tableaux.
//!
//! The staircase `(1 0^z)^y` of width `w` is the cylindric diagram with `y` rows of `w` cells,
//! each row shifted `z` steps to the left of the row above. Its cells are `(i, d)` with
//! `i in Z/y` (row) and `d = 1..=w` (position). The cover relations are
//! `(i, d) -> (i, d + 1)` and `(i, d) -> (i + 1 mod y, d + z)`. For `z = 1` this is the
//! cylindrical rectangle `a^b` with `a = w`, `b = y`.
//!
//! A standard cylindric tableau `T` (a linear extension, `n = yw` labels) is recorded by its
//! **hop word** `u_1 ... u_n`, where `u_r in {0, ..., y-1}` is the row of label `r`. The hop
//! words are exactly the words with each letter `w` times whose prefix counts `h` satisfy
//! `h_{i+1} <= h_i + z` for all `i`, cyclically.
//!
//! The **charge** of the row reading word (rows `y-1, ..., 0`, each left to right) is
//! `sum_{r in A(u)} (n - r)`, where `A(u) = { r in [n-1] : u_{r+1} <= u_r }` (weak descents
//! of the hop word). The **cocharge** is `binom(n, 2) - charge`.
//!
//! All polynomials are dense coefficient vectors in ascending degree with `BigInt`
//! coefficients and no trailing zeros. Exponents of the formula-based routines are computed
//! with checked `i64` arithmetic; out-of-range parameters panic rather than overflow.

use num_bigint::BigInt;
use num_traits::{One, Zero};
use std::collections::{BTreeMap, HashMap};

fn trim(mut v: Vec<BigInt>) -> Vec<BigInt> {
    while v.last().is_some_and(|c| c.is_zero()) {
        v.pop();
    }
    v
}

fn check_shape(shift: usize, rows: usize, width: usize) {
    assert!(shift >= 1, "the staircase shift z must be at least 1");
    assert!(rows >= 1, "the number of rows y must be at least 1");
    assert!(width >= 1, "the width w must be at least 1");
}

/// Dynamic programme over (prefix counts, last letter); returns `charge` coefficients.
fn hop_word_charge_dp(shift: usize, rows: usize, width: usize) -> Vec<BigInt> {
    let n = rows * width;
    let max_deg = n * (n - 1) / 2;
    // state: (h, last letter) -> polynomial
    let mut cur: HashMap<(Vec<usize>, usize), Vec<BigInt>> = HashMap::new();
    // place the first letter separately (no charge contribution)
    for i in 0..rows {
        let mut h = vec![0; rows];
        h[i] = 1;
        if admissible(&h, shift) {
            let mut p = vec![BigInt::zero(); max_deg + 1];
            p[0] = BigInt::one();
            cur.insert((h, i), p);
        }
    }
    for r in 1..n {
        // placing label r + 1; the pair (r, r + 1) contributes n - r on a weak descent
        let mut nxt: HashMap<(Vec<usize>, usize), Vec<BigInt>> = HashMap::new();
        for ((h, last), poly) in &cur {
            for i in 0..rows {
                if h[i] == width {
                    continue;
                }
                let mut h2 = h.clone();
                h2[i] += 1;
                if !admissible(&h2, shift) {
                    continue;
                }
                let add = if i <= *last { n - r } else { 0 };
                let tgt = nxt
                    .entry((h2, i))
                    .or_insert_with(|| vec![BigInt::zero(); max_deg + 1]);
                for (e, c) in poly.iter().enumerate() {
                    if !c.is_zero() {
                        tgt[e + add] += c;
                    }
                }
            }
        }
        cur = nxt;
    }
    let mut tot = vec![BigInt::zero(); max_deg + 1];
    for poly in cur.values() {
        for (e, c) in poly.iter().enumerate() {
            tot[e] += c;
        }
    }
    trim(tot)
}

fn admissible(h: &[usize], shift: usize) -> bool {
    let y = h.len();
    (0..y).all(|i| h[(i + 1) % y] <= h[i] + shift)
}

/// The charge polynomial `X(q) = sum_T q^{charge(rw T)}` of the staircase `(1 0^z)^y` of
/// width `w` (see the module documentation for the conventions).
///
/// Computed exactly by a dynamic programme over hop-word prefixes; the cost is
/// `O(#states * y * n^2)` big-integer additions, where the states are the admissible prefix
/// count vectors. Panics if `z`, `y` or `w` is zero.
///
/// For `z = 1` this is `X(a, b; q)` of the cylindrical rectangle `a^b`; for example
/// `X(1, 3; q) = 1 + 2q + 2q^2 + q^3`.
pub fn staircase_charge_polynomial_bigint(shift: usize, rows: usize, width: usize) -> Vec<BigInt> {
    check_shape(shift, rows, width);
    hop_word_charge_dp(shift, rows, width)
}

/// The cocharge polynomial `CC(q) = q^{binom(n,2)} X(1/q) = sum_T q^{cocharge(rw T)}`.
pub fn staircase_cocharge_polynomial_bigint(
    shift: usize,
    rows: usize,
    width: usize,
) -> Vec<BigInt> {
    let x = staircase_charge_polynomial_bigint(shift, rows, width);
    let n = rows * width;
    let top = n * (n - 1) / 2;
    let mut cc = vec![BigInt::zero(); top + 1];
    for (e, c) in x.iter().enumerate() {
        cc[top - e] += c;
    }
    trim(cc)
}

/// The number of standard cylindric tableaux of the staircase `(1 0^z)^y` of width `w`
/// (the number of hop words). Conventions used by
/// [`staircase_cocharge_root_of_unity_prediction`]: width `0` gives `1` (the empty tableau),
/// and `z = 0` with `y >= 2` and `w >= 1` gives `0`.
pub fn staircase_tableau_count_bigint(shift: usize, rows: usize, width: usize) -> BigInt {
    if width == 0 {
        return BigInt::one();
    }
    if shift == 0 {
        return if rows == 1 {
            BigInt::one()
        } else {
            BigInt::zero()
        };
    }
    // count-only dynamic programme
    let mut cur: HashMap<Vec<usize>, BigInt> = HashMap::new();
    cur.insert(vec![0; rows], BigInt::one());
    for _ in 0..rows * width {
        let mut nxt: HashMap<Vec<usize>, BigInt> = HashMap::new();
        for (h, c) in &cur {
            for i in 0..rows {
                if h[i] == width {
                    continue;
                }
                let mut h2 = h.clone();
                h2[i] += 1;
                if admissible(&h2, shift) {
                    *nxt.entry(h2).or_insert_with(BigInt::zero) += c;
                }
            }
        }
        cur = nxt;
    }
    cur.values().fold(BigInt::zero(), |a, b| a + b)
}

/// The conjectured value of the cocharge polynomial at a primitive `d`-th root of unity
/// (conjecture A of the project notes): for `d | n`,
/// `CC(zeta_d) = |CSYT((1 0^{z'})^y, width n/(d y))|` with `z' + 1 = (z+1)/gcd(d, z+1)`,
/// read as `0` when `z' = 0` or `y` does not divide `n/d`. Returns `None` if `d` does not
/// divide `n`.
pub fn staircase_cocharge_root_of_unity_prediction(
    shift: usize,
    rows: usize,
    width: usize,
    d: usize,
) -> Option<BigInt> {
    let n = rows * width;
    if d == 0 || !n.is_multiple_of(d) {
        return None;
    }
    let m = n / d;
    if !m.is_multiple_of(rows) {
        return Some(BigInt::zero());
    }
    let p = shift + 1;
    let zp = p / gcd(d, p) - 1;
    if zp == 0 {
        return Some(if m == 0 {
            BigInt::one()
        } else {
            BigInt::zero()
        });
    }
    Some(staircase_tableau_count_bigint(zp, rows, m / rows))
}

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

/// Gaussian binomial `[n; k]_q` (ascending coefficients), `0` outside `0 <= k <= n`.
fn q_binomial(n: usize, k: usize, cache: &mut HashMap<(usize, usize), Vec<BigInt>>) -> Vec<BigInt> {
    if k > n {
        return vec![];
    }
    if k == 0 || k == n {
        return vec![BigInt::one()];
    }
    if let Some(v) = cache.get(&(n, k)) {
        return v.clone();
    }
    // [n;k] = [n-1;k-1] + q^k [n-1;k]
    let a = q_binomial(n - 1, k - 1, cache);
    let b = q_binomial(n - 1, k, cache);
    let mut r = vec![BigInt::zero(); (k * (n - k)) + 1];
    for (e, c) in a.iter().enumerate() {
        r[e] += c;
    }
    for (e, c) in b.iter().enumerate() {
        r[e + k] += c;
    }
    let r = trim(r);
    cache.insert((n, k), r.clone());
    r
}

fn poly_mul(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r = vec![BigInt::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            r[i + j] += x * y;
        }
    }
    r
}

/// The `q`-multinomial `[n; c_0, ..., c_{k-1}]_q` (`n = sum c_i`), empty if some part is
/// negative.
fn q_multinomial(parts: &[i64], cache: &mut HashMap<(usize, usize), Vec<BigInt>>) -> Vec<BigInt> {
    if parts.iter().any(|&c| c < 0) {
        return vec![];
    }
    let mut rest: usize = parts.iter().map(|&c| c as usize).sum();
    let mut acc = vec![BigInt::one()];
    for &c in parts {
        let c = c as usize;
        acc = poly_mul(&acc, &q_binomial(rest, c, cache));
        rest -= c;
    }
    acc
}

fn add_shifted(acc: &mut BTreeMap<i64, BigInt>, p: &[BigInt], shift: i64, sign: i64) {
    for (e, c) in p.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let k = shift + e as i64;
        let v = acc.entry(k).or_insert_with(BigInt::zero);
        if sign > 0 {
            *v += c;
        } else {
            *v -= c;
        }
    }
}

fn laurent_to_polynomial(acc: BTreeMap<i64, BigInt>, what: &str) -> Vec<BigInt> {
    let acc: BTreeMap<i64, BigInt> = acc.into_iter().filter(|(_, c)| !c.is_zero()).collect();
    let Some((&lo, _)) = acc.iter().next() else {
        return vec![];
    };
    assert!(
        lo >= 0,
        "{what}: negative exponent {lo} survived cancellation"
    );
    let hi = *acc.keys().last().unwrap();
    let mut v = vec![BigInt::zero(); (hi + 1) as usize];
    for (k, c) in acc {
        v[k as usize] = c;
    }
    v
}

fn permutations(n: usize) -> Vec<(Vec<usize>, i64)> {
    // all permutations of 0..n with their signs
    fn rec(
        prefix: &mut Vec<usize>,
        used: &mut Vec<bool>,
        n: usize,
        out: &mut Vec<(Vec<usize>, i64)>,
    ) {
        if prefix.len() == n {
            let inv = (0..n)
                .flat_map(|i| (i + 1..n).map(move |j| (i, j)))
                .filter(|&(i, j)| prefix[i] > prefix[j])
                .count();
            out.push((prefix.clone(), if inv % 2 == 0 { 1 } else { -1 }));
            return;
        }
        for v in 0..n {
            if !used[v] {
                used[v] = true;
                prefix.push(v);
                rec(prefix, used, n, out);
                prefix.pop();
                used[v] = false;
            }
        }
    }
    let mut out = vec![];
    rec(&mut vec![], &mut vec![false; n], n, &mut out);
    out
}

/// All `alpha in Z^y` with `sum alpha = 0` and `|alpha_i| <= bound`.
fn root_lattice_box(y: usize, bound: i64) -> Vec<Vec<i64>> {
    let mut out = vec![];
    let mut cur = vec![0i64; y];
    fn rec(i: usize, y: usize, bound: i64, cur: &mut Vec<i64>, out: &mut Vec<Vec<i64>>) {
        if i == y - 1 {
            let s: i64 = cur[..y - 1].iter().sum();
            if s.abs() <= bound {
                cur[y - 1] = -s;
                out.push(cur.clone());
            }
            return;
        }
        for v in -bound..=bound {
            cur[i] = v;
            rec(i + 1, y, bound, cur, out);
        }
    }
    rec(0, y, bound, &mut cur, &mut out);
    out
}

/// The bosonic (affine alternating-sum) formula for the staircase charge polynomial:
///
/// `X(q) = q^{n(w-1)/2} sum_{tau, alpha} sign(tau) q^{E(eps)} [n; w + (z+1) eps]_q`,
///
/// where `eps = eps_tau + y alpha`, `eps_tau = (i - tau(i))_i`, `alpha` runs over the root
/// lattice `{alpha in Z^y : sum alpha = 0}`, `K = (z+1) y`, and
/// `E(eps) = (z+1) ((K-1)|eps|^2 + 2 sum_i i' eps_i) / (2y)` with `i' = i` for `i >= 1` and
/// `0' = y`.
///
/// This is the exact formula reported by the worker analysis of the project (derived from the
/// Schilling--Shimozono bosonic formula for level-restricted paths); this routine evaluates it
/// so that it can be compared with [`staircase_charge_polynomial_bigint`]. Panics if a
/// fractional exponent or a surviving negative exponent occurs.
pub fn staircase_bosonic_charge_polynomial_bigint(
    shift: usize,
    rows: usize,
    width: usize,
) -> Vec<BigInt> {
    check_shape(shift, rows, width);
    let (y, w, p) = (rows as i64, width as i64, (shift + 1) as i64);
    let n = y * w;
    let big_k = p * y;
    let bound = w / p + 2;
    let mut cache = HashMap::new();
    let mut acc = BTreeMap::new();
    let perms = permutations(rows);
    for alpha in root_lattice_box(rows, bound) {
        for (tau, sign) in &perms {
            let eps: Vec<i64> = (0..rows)
                .map(|i| i as i64 - tau[i] as i64 + y * alpha[i])
                .collect();
            let parts: Vec<i64> = eps.iter().map(|e| w + p * e).collect();
            if parts.iter().any(|&c| c < 0) {
                continue;
            }
            let norm: i64 = eps.iter().map(|e| e * e).sum();
            let lin: i64 = eps
                .iter()
                .enumerate()
                .map(|(i, e)| if i == 0 { y * e } else { i as i64 * e })
                .sum();
            let num = p * ((big_k - 1) * norm + 2 * lin);
            assert!(
                num % (2 * y) == 0,
                "fractional exponent in the bosonic formula"
            );
            let base = n * (w - 1);
            assert!(base % 2 == 0);
            let expo = base / 2 + num / (2 * y);
            let mult = q_multinomial(&parts, &mut cache);
            add_shifted(&mut acc, &mult, expo, *sign);
        }
    }
    laurent_to_polynomial(acc, "staircase bosonic formula")
}

/// The cocharge form of the Krattenthaler--Mohanty formula for the staircase `(1 0^z)^y` of
/// width `w` (`p = z + 1`, `n = yw`):
///
/// `CC(q) = sum_{eps in O} sign(tau) q^{T_p(eps)} [n; w - p eps_0, ..., w - p eps_{y-1}]_q`,
///
/// where `O` is the disjoint union over `tau in S_y` of `eps_tau + yQ`,
/// `eps_tau = (i - tau(i))_i`, `Q = {alpha in Z^y : sum alpha = 0}`, and
/// `T_p(eps) = (p / y) (|eps|^2 / 2 - sum_i i eps_i)` (an integer on `O`). This is
/// Krattenthaler--Mohanty, Discrete Math. 126 (1994), Theorem B, with `tau` the reversal,
/// `mu_i = z (y - 1 - i)` and `K = zy`, written in the variables `eps`.
pub fn staircase_cocharge_bosonic_polynomial_bigint(
    shift: usize,
    rows: usize,
    width: usize,
) -> Vec<BigInt> {
    check_shape(shift, rows, width);
    let (y, w, p) = (rows as i64, width as i64, (shift + 1) as i64);
    let bound = w / p + 2;
    let mut cache = HashMap::new();
    let mut acc = BTreeMap::new();
    let perms = permutations(rows);
    for alpha in root_lattice_box(rows, bound) {
        for (tau, sign) in &perms {
            let eps: Vec<i64> = (0..rows)
                .map(|i| i as i64 - tau[i] as i64 + y * alpha[i])
                .collect();
            let parts: Vec<i64> = eps.iter().map(|e| w - p * e).collect();
            if parts.iter().any(|&c| c < 0) {
                continue;
            }
            let norm: i64 = eps.iter().map(|e| e * e).sum();
            let lin: i64 = eps.iter().enumerate().map(|(i, e)| i as i64 * e).sum();
            let num = p * (norm - 2 * lin);
            assert!(num % (2 * y) == 0, "T_p(eps) is not an integer");
            let mult = q_multinomial(&parts, &mut cache);
            add_shifted(&mut acc, &mult, num / (2 * y), *sign);
        }
    }
    laurent_to_polynomial(acc, "staircase cocharge bosonic formula")
}

/// Krattenthaler--Mohanty, *q-generalization of a ballot problem*, Discrete Math. 126 (1994),
/// Theorem B: the generating function `sum_p q^{maj_tau p}` over lattice paths in `Z^n` from
/// `mu` to `lambda` (unit positive steps) that never touch the hyperplanes
/// `x_i - x_{i+1} = -1` (`1 <= i < n`) and `x_n - x_1 = -1 - K`, with
/// `maj_tau p = maj(tau(p))` for a permutation `tau` of the letters.
///
/// The determinant of the theorem is expanded over permutations `sigma`, giving
/// `sum_{u, sigma} sign(sigma) q^{sum_s T_tau(s, sigma(s), u, mu)} [|lambda - mu|; m]_q` with
/// `m_s = lambda_s - s - mu_{sigma(s)} + sigma(s) - (K + n) u_s`. Letters, `mu`, `lambda` and
/// `tau` are one-indexed in the paper; here the slices are zero-indexed and `tau` is a
/// permutation of `0..n`. Requires `mu, lambda` in the chamber
/// `x_1 >= ... >= x_n >= x_1 - K` and `lambda >= mu`.
pub fn krattenthaler_mohanty_major_polynomial_bigint(
    big_k: i64,
    mu: &[i64],
    lambda: &[i64],
    tau: &[usize],
) -> Vec<BigInt> {
    let n = mu.len();
    assert!(n >= 1 && lambda.len() == n && tau.len() == n);
    let ni = n as i64;
    // c(j) = [tau(j) < tau(j+1)], c(n) = [tau(n) < tau(1)], j one-indexed
    let c: Vec<i64> = (0..n)
        .map(|j| i64::from(tau[j] < tau[(j + 1) % n]))
        .collect();
    let csum: i64 = c.iter().sum();
    // S(k, l) = sum_{j=k}^{l-1} c(j) with the paper's sign convention (one-indexed)
    let ssum = |k: i64, l: i64| -> i64 {
        if k <= l {
            (k..l).map(|j| c[(j - 1) as usize]).sum()
        } else {
            -(l..k).map(|j| c[(j - 1) as usize]).sum::<i64>()
        }
    };
    let total: i64 = lambda.iter().sum::<i64>() - mu.iter().sum::<i64>();
    assert!(total >= 0);
    let bound = total / (big_k + ni) + 2;
    let mut cache = HashMap::new();
    let mut acc = BTreeMap::new();
    let perms = permutations(n);
    for u in root_lattice_box(n, bound) {
        for (sigma, sign) in &perms {
            let mut parts = Vec::with_capacity(n);
            let mut twice_expo: i64 = 0;
            for s in 1..=ni {
                let t = sigma[(s - 1) as usize] as i64 + 1;
                let us = u[(s - 1) as usize];
                let mut_ = mu[(t - 1) as usize];
                parts.push(lambda[(s - 1) as usize] - s - mut_ + t - (big_k + ni) * us);
                // 2 T_tau(s, t, u, mu)
                twice_expo += 2 * (mut_ - t) * ssum(t, s)
                    + 2 * us * (mut_ - t) * csum
                    + (big_k + ni) * us * us * csum
                    + 2 * (big_k + ni) * us * ssum(1, s);
            }
            if parts.iter().any(|&m| m < 0) {
                continue;
            }
            assert!(twice_expo % 2 == 0, "half-integral exponent in Theorem B");
            let mult = q_multinomial(&parts, &mut cache);
            add_shifted(&mut acc, &mult, twice_expo / 2, *sign);
        }
    }
    laurent_to_polynomial(acc, "Krattenthaler-Mohanty Theorem B")
}

/// Brute-force `sum_p q^{maj_tau p}` over the lattice paths of
/// [`krattenthaler_mohanty_major_polynomial_bigint`] (an independent check; exponential).
pub fn chamber_paths_major_polynomial_bruteforce(
    big_k: i64,
    mu: &[i64],
    lambda: &[i64],
    tau: &[usize],
) -> Vec<BigInt> {
    let n = mu.len();
    let total: usize = (lambda.iter().sum::<i64>() - mu.iter().sum::<i64>()) as usize;
    let mut out = vec![BigInt::zero(); total * total / 2 + 2];
    let mut word = Vec::with_capacity(total);
    let mut pos = mu.to_vec();
    fn ok(x: &[i64], k: i64) -> bool {
        let n = x.len();
        (0..n - 1).all(|i| x[i] - x[i + 1] >= 0) && x[n - 1] - x[0] >= -k
    }
    #[allow(clippy::too_many_arguments)]
    fn rec(
        pos: &mut Vec<i64>,
        word: &mut Vec<usize>,
        lambda: &[i64],
        tau: &[usize],
        k: i64,
        total: usize,
        out: &mut Vec<BigInt>,
    ) {
        if word.len() == total {
            let maj: usize = (0..total - 1)
                .filter(|&i| tau[word[i]] > tau[word[i + 1]])
                .map(|i| i + 1)
                .sum();
            out[maj] += 1;
            return;
        }
        for i in 0..pos.len() {
            if pos[i] < lambda[i] {
                pos[i] += 1;
                if ok(pos, k) {
                    word.push(i);
                    rec(pos, word, lambda, tau, k, total, out);
                    word.pop();
                }
                pos[i] -= 1;
            }
        }
    }
    let _ = n;
    if total == 0 {
        return vec![BigInt::one()];
    }
    rec(&mut pos, &mut word, lambda, tau, big_k, total, &mut out);
    trim(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ints(v: &[i64]) -> Vec<BigInt> {
        v.iter().map(|&x| BigInt::from(x)).collect()
    }

    #[test]
    fn small_charge_polynomials() {
        assert_eq!(
            staircase_charge_polynomial_bigint(1, 3, 1),
            ints(&[1, 2, 2, 1])
        );
        assert_eq!(
            staircase_charge_polynomial_bigint(1, 4, 1),
            ints(&[1, 3, 5, 6, 5, 3, 1])
        );
        // two rows: X(a,2;q) = q^{a(a-1)} prod_{j<=a} (1+q^j)
        for a in 1..=6usize {
            let mut p = vec![BigInt::one()];
            for j in 1..=a {
                let mut f = vec![BigInt::zero(); j + 1];
                f[0] = BigInt::one();
                f[j] = BigInt::one();
                p = poly_mul(&p, &f);
            }
            let mut want = vec![BigInt::zero(); a * (a - 1)];
            want.extend(p);
            assert_eq!(staircase_charge_polynomial_bigint(1, 2, a), want, "a = {a}");
        }
    }

    #[test]
    fn counts_match_known_values() {
        // |CSYT(a^3)| = 6 * 8^{a-1}; |CSYT(k^4)| = 24, 1088, 50688
        for a in 1..=5u32 {
            assert_eq!(
                staircase_tableau_count_bigint(1, 3, a as usize),
                BigInt::from(6) * BigInt::from(8).pow(a - 1)
            );
        }
        let c4: Vec<BigInt> = (1..=3)
            .map(|a| staircase_tableau_count_bigint(1, 4, a))
            .collect();
        assert_eq!(c4, ints(&[24, 1088, 50688]));
        // (100)^2: 2 * 3^{w-1}
        for w in 1..=6u32 {
            assert_eq!(
                staircase_tableau_count_bigint(2, 2, w as usize),
                BigInt::from(2) * BigInt::from(3).pow(w - 1)
            );
        }
        for (z, y, w) in [(1, 3, 3), (2, 3, 2), (3, 2, 4)] {
            let x = staircase_charge_polynomial_bigint(z, y, w);
            let s = x.iter().fold(BigInt::zero(), |a, b| a + b);
            assert_eq!(s, staircase_tableau_count_bigint(z, y, w));
        }
    }

    #[test]
    fn bosonic_formula_matches_dp() {
        for (z, y, w) in [
            (1, 2, 5),
            (1, 3, 3),
            (1, 4, 2),
            (1, 5, 1),
            (2, 2, 4),
            (2, 3, 2),
            (3, 2, 3),
            (2, 4, 1),
        ] {
            assert_eq!(
                staircase_bosonic_charge_polynomial_bigint(z, y, w),
                staircase_charge_polynomial_bigint(z, y, w),
                "z={z} y={y} w={w}"
            );
        }
    }

    #[test]
    fn krattenthaler_mohanty_theorem_b_matches_bruteforce() {
        for (k, mu, lam) in [
            (2, vec![1, 0], vec![3, 2]),
            (3, vec![2, 1, 0], vec![4, 3, 2]),
            (3, vec![1, 1, 0], vec![3, 2, 2]),
            (4, vec![3, 2, 1, 0], vec![4, 3, 2, 1]),
        ] {
            let n = mu.len();
            let id: Vec<usize> = (0..n).collect();
            let rev: Vec<usize> = (0..n).rev().collect();
            for tau in [id, rev] {
                assert_eq!(
                    krattenthaler_mohanty_major_polynomial_bigint(k, &mu, &lam, &tau),
                    chamber_paths_major_polynomial_bruteforce(k, &mu, &lam, &tau),
                    "K={k} mu={mu:?} lambda={lam:?} tau={tau:?}"
                );
            }
        }
    }

    #[test]
    fn rectangle_cocharge_is_krattenthaler_mohanty_lesser_index() {
        // hop words of a^b are the chamber paths from mu = (b-1, ..., 0) to mu + (a, ..., a)
        // with K = b; reversing and relabelling i -> b-1-i turns cocharge into maj_{tau_0}.
        for (b, a) in [(2usize, 4usize), (3, 2), (3, 3), (4, 2), (5, 1)] {
            let mu: Vec<i64> = (0..b).map(|i| (b - 1 - i) as i64).collect();
            let lam: Vec<i64> = mu.iter().map(|m| m + a as i64).collect();
            let rev: Vec<usize> = (0..b).rev().collect();
            assert_eq!(
                krattenthaler_mohanty_major_polynomial_bigint(b as i64, &mu, &lam, &rev),
                staircase_cocharge_polynomial_bigint(1, b, a),
                "a={a} b={b}"
            );
        }
    }

    #[test]
    fn staircase_cocharge_is_krattenthaler_mohanty_lesser_index() {
        // staircase hop words are the chamber paths from mu_i = z(y-1-i) to mu + w, K = zy
        for (z, y, w) in [
            (2usize, 2usize, 4usize),
            (2, 3, 2),
            (3, 2, 3),
            (3, 3, 2),
            (2, 4, 1),
        ] {
            let mu: Vec<i64> = (0..y).map(|i| (z * (y - 1 - i)) as i64).collect();
            let lam: Vec<i64> = mu.iter().map(|m| m + w as i64).collect();
            let rev: Vec<usize> = (0..y).rev().collect();
            assert_eq!(
                krattenthaler_mohanty_major_polynomial_bigint((z * y) as i64, &mu, &lam, &rev),
                staircase_cocharge_polynomial_bigint(z, y, w),
                "z={z} y={y} w={w}"
            );
        }
    }

    #[test]
    fn cocharge_bosonic_formula_matches_dp() {
        for (z, y, w) in [
            (1, 3, 3),
            (1, 4, 2),
            (2, 2, 5),
            (2, 3, 2),
            (3, 2, 4),
            (2, 4, 1),
        ] {
            assert_eq!(
                staircase_cocharge_bosonic_polynomial_bigint(z, y, w),
                staircase_cocharge_polynomial_bigint(z, y, w),
                "z={z} y={y} w={w}"
            );
        }
    }

    #[test]
    fn cocharge_values_at_roots_of_unity() {
        use polytool::cyclic_sieving::root_of_unity_evaluation_bigint;
        for (z, y, w) in [(1, 3, 3), (2, 2, 6), (2, 3, 4), (3, 2, 6)] {
            let cc = staircase_cocharge_polynomial_bigint(z, y, w);
            let n = y * w;
            for d in (2..=n).filter(|d| n % d == 0) {
                let ev = root_of_unity_evaluation_bigint(&cc, d, 1);
                assert_eq!(
                    ev.integer_value,
                    staircase_cocharge_root_of_unity_prediction(z, y, w, d),
                    "z={z} y={y} w={w} d={d}"
                );
            }
        }
    }
}
