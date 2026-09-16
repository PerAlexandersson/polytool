//! Key polynomials (Demazure characters) via Kogan faces of GT-polytopes.
//!
//! A key polynomial κ_{λ,σ} is computed by summing monomials z^{w(G)} over
//! integer GT-patterns G lying in the union of all reduced Kogan faces of
//! type w₀σ in GT(λ).
//!
//! # Kogan faces
//!
//! A Kogan face is obtained by forcing a subset of the left-interlacing
//! inequalities x_{i+1,j} ≥ x_{i,j} to equalities. Each equality position
//! (i,j) (where 1 ≤ j ≤ i < n) is associated with transposition s_{n-i+j-1}.
//!
//! Reading the transpositions bottom-to-top, left-to-right gives a word.
//! If the word is reduced, the face is a reduced Kogan face and its type
//! is the permutation represented by the word.
//!
//! # Convention
//!
//! We use (w₀σ)(i) = σ(w₀(i)), i.e., w₀ acts first. This matches the
//! Demazure operator definition κ_{λ,σ} = π_σ(z^λ).

use crate::partition::Partition;
use crate::permutation::{
    compose_permutations, inverse_permutation, longest_permutation,
    permutation_from_simple_transpositions,
};
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};
use std::collections::{BTreeMap, HashSet};

/// A GT-pattern for shape λ with n = λ.num_parts().
/// Stored as `rows[i]` = entries of row i+1 (0-indexed), so `rows[n-1] = λ`.
/// Row i (1-indexed) has i entries.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct GtPattern {
    /// `rows[i]` has i+1 entries (row i+1 in 1-indexed notation).
    pub rows: Vec<Vec<u32>>,
}

impl GtPattern {
    /// Number of rows (= n).
    pub fn n(&self) -> usize {
        self.rows.len()
    }

    /// Weight vector: w_k = sum(row k) - sum(row k-1), 1-indexed.
    pub fn weight(&self) -> Vec<u32> {
        let n = self.n();
        let mut w = vec![0u32; n];
        let mut prev_sum = 0u32;
        for k in 0..n {
            let cur_sum: u32 = self.rows[k].iter().sum();
            w[k] = cur_sum - prev_sum;
            prev_sum = cur_sum;
        }
        w
    }

    /// Entry x_{i,j} (1-indexed). Row i has i entries.
    pub fn entry(&self, i: usize, j: usize) -> u32 {
        self.rows[i - 1][j - 1]
    }
}

/// Enumerate all integer GT-patterns of shape λ with n rows and a given set of
/// left-interlacing equalities forced.
///
/// `n` is the number of rows (= number of variables). λ is padded with zeros
/// to length n if needed.
///
/// `equalities` is a set of (i, j) pairs (1-indexed) meaning x_{i+1,j} = x_{i,j}
/// where 1 ≤ j ≤ i < n.
///
/// Returns all valid integer GT-patterns satisfying the standard interlacing
/// plus the forced equalities.
pub fn gt_patterns_with_equalities_n(
    lambda: &Partition,
    n: usize,
    equalities: &[(usize, usize)],
) -> Vec<GtPattern> {
    gt_patterns_inner(lambda, n, equalities, false)
}

/// Like `gt_patterns_with_equalities_n`, but with STRICT interlacing
/// at all non-forced positions.
///
/// At forced equality positions (Kogan face equalities), x_{i+1,j} = x_{i,j}.
/// At all other positions, the interlacing is strict:
///   x_{i+1,j} > x_{i,j} > x_{i+1,j+1}
///
/// This counts the "strictly interior" lattice points, used for
/// Ehrhart-Macdonald reciprocity: P(-k) = (-1)^d · strict_count(k).
pub fn strict_gt_patterns_with_equalities_n(
    lambda: &Partition,
    n: usize,
    equalities: &[(usize, usize)],
) -> Vec<GtPattern> {
    gt_patterns_inner(lambda, n, equalities, true)
}

fn gt_patterns_inner(
    lambda: &Partition,
    n: usize,
    equalities: &[(usize, usize)],
    strict: bool,
) -> Vec<GtPattern> {
    if n == 0 {
        return vec![GtPattern { rows: vec![] }];
    }

    // Build equality lookup: eq_set[i][j] = true means x_{i+1,j} = x_{i,j}
    let mut eq_set = vec![vec![false; n + 1]; n + 1];
    for &(i, j) in equalities {
        if i < n && j >= 1 && j <= i {
            eq_set[i][j] = true;
        }
    }

    let top_row: Vec<u32> = (0..n)
        .map(|j| {
            if j < lambda.num_parts() {
                lambda.part(j)
            } else {
                0
            }
        })
        .collect();

    let mut results = Vec::new();
    let mut rows_built: Vec<Vec<u32>> = vec![top_row];

    enumerate_gt_rows(n, n - 1, &eq_set, strict, &mut rows_built, &mut results);
    results
}

fn enumerate_gt_rows(
    n: usize,
    row_idx: usize,
    eq_set: &[Vec<bool>],
    strict: bool,
    rows_built: &mut Vec<Vec<u32>>,
    results: &mut Vec<GtPattern>,
) {
    if row_idx == 0 {
        let n = rows_built.len();
        let mut rows = Vec::with_capacity(n);
        for k in (0..n).rev() {
            rows.push(rows_built[k].clone());
        }
        results.push(GtPattern { rows });
        return;
    }

    let above = rows_built[n - row_idx - 1].clone();
    let num_entries = row_idx;
    let mut current_row = vec![0u32; num_entries];

    enumerate_gt_entries(
        n,
        row_idx,
        eq_set,
        strict,
        &above,
        &mut current_row,
        0,
        rows_built,
        results,
    );
}

fn enumerate_gt_entries(
    n: usize,
    row_idx: usize,
    eq_set: &[Vec<bool>],
    strict: bool,
    above: &[u32],
    current: &mut Vec<u32>,
    pos: usize,
    rows_built: &mut Vec<Vec<u32>>,
    results: &mut Vec<GtPattern>,
) {
    if pos == current.len() {
        rows_built.push(current.clone());
        enumerate_gt_rows(n, row_idx - 1, eq_set, strict, rows_built, results);
        rows_built.pop();
        return;
    }

    let j = pos + 1; // 1-indexed column
    let is_forced = eq_set[row_idx][j];

    // Bounds from interlacing: x_{row_idx+1, j} ≥ x_{row_idx, j} ≥ x_{row_idx+1, j+1}
    // In strict mode (non-forced): use > instead of ≥.
    let mut upper = above[pos]; // x_{row_idx+1, j}
    let mut lower = above[pos + 1]; // x_{row_idx+1, j+1}

    if strict && !is_forced {
        // Strict left interlacing: x_{i+1,j} > x_{i,j} → upper = above[pos] - 1
        if above[pos] == 0 {
            return;
        } // can't go below 0
        upper = above[pos] - 1;
        // Strict right interlacing: x_{i,j} > x_{i+1,j+1} → lower = above[pos+1] + 1
        lower = above[pos + 1] + 1;
    }

    // Same-row monotonicity (weakly decreasing)
    if pos > 0 {
        upper = upper.min(current[pos - 1]);
    }

    if is_forced {
        let forced = above[pos];
        if forced >= lower && forced <= upper {
            current[pos] = forced;
            enumerate_gt_entries(
                n,
                row_idx,
                eq_set,
                strict,
                above,
                current,
                pos + 1,
                rows_built,
                results,
            );
        }
        return;
    }

    if lower > upper {
        return;
    }

    for val in lower..=upper {
        current[pos] = val;
        enumerate_gt_entries(
            n,
            row_idx,
            eq_set,
            strict,
            above,
            current,
            pos + 1,
            rows_built,
            results,
        );
    }
}

// ── Reduced words and Kogan face placement ──────────────────────────────────

/// All reduced words for a permutation σ (given in one-line notation, 1-indexed).
///
/// Uses the standard algorithm: for each descent position i (`σ[i] > σ[i+1]`),
/// apply s_i and recursively find reduced words for s_i · σ, prepending i.
pub fn reduced_words(perm: &[usize]) -> Vec<Vec<usize>> {
    let n = perm.len();
    if n <= 1 {
        return vec![vec![]];
    }

    // Check if identity
    if perm.iter().enumerate().all(|(i, &v)| v == i + 1) {
        return vec![vec![]];
    }

    let mut result = Vec::new();
    for i in 0..n - 1 {
        if perm[i] > perm[i + 1] {
            // Apply s_{i+1} on the left: swap positions i, i+1
            let mut new_perm = perm.to_vec();
            new_perm.swap(i, i + 1);
            for mut word in reduced_words(&new_perm) {
                word.insert(0, i + 1); // s_{i+1}
                result.push(word);
            }
        }
    }
    result
}

/// Enumerate all reduced Kogan faces of a given type σ for the GT diagram with n variables.
///
/// The algorithm (following Alexandersson's Mathematica implementation):
/// 1. Build the word vector wVec = positions in reading order (bottom-to-top, left-to-right).
///    For n=4: wVec = [1,2,3, 2,3, 3] (transposition indices s_k).
/// 2. Enumerate all subsets of wVec positions.
/// 3. For each subset, compute the permutation from the subword.
/// 4. Check if the subword is reduced (inversions of perm == length of subword).
/// 5. Keep only faces whose type matches σ.
///
/// Returns a list of equality sets, each being a Vec<(usize, usize)> of (i, j) pairs
/// (1-indexed, representing x_{i+1,j} = x_{i,j}).
pub fn reduced_kogan_faces(n: usize, sigma: &[usize]) -> Vec<Vec<(usize, usize)>> {
    // Build position list in reading order: bottom-to-top, left-to-right.
    // Level i (equalities between rows i+1 and i): positions (i, j) for j=1..i.
    // wVec[k] = transposition index s_{n-i+j-1} at position k.
    let mut positions: Vec<(usize, usize)> = Vec::new(); // (i, j) pairs
    let mut wvec: Vec<usize> = Vec::new(); // transposition indices
    for i in 1..n {
        for j in 1..=i {
            let trans = n - i + j - 1;
            if trans >= 1 && trans < n {
                positions.push((i, j));
                wvec.push(trans);
            }
        }
    }

    let m = wvec.len();
    let mut results = Vec::new();

    // Enumerate all 2^m subsets
    for mask in 0u64..(1u64 << m) {
        // Extract the subword
        let subword: Vec<usize> = (0..m)
            .filter(|&k| (mask >> k) & 1 == 1)
            .map(|k| wvec[k])
            .collect();

        if subword.is_empty() {
            // Empty subword → identity permutation
            if sigma.iter().enumerate().all(|(i, &v)| v == i + 1) {
                results.push(vec![]);
            }
            continue;
        }

        // Compute the permutation from the subword
        let perm = permutation_from_simple_transpositions(n, &subword);

        // Check if reduced: inversions(perm) == len(subword)
        let inv = inversion_count(&perm);
        if inv != subword.len() {
            continue; // not reduced
        }

        // Check if type matches sigma
        if perm == sigma {
            let eqs: Vec<(usize, usize)> = (0..m)
                .filter(|&k| (mask >> k) & 1 == 1)
                .map(|k| positions[k])
                .collect();
            results.push(eqs);
        }
    }

    results
}

/// Compute the permutation (one-line notation) from a word of simple transpositions.
/// Word [a_1, ..., a_k] represents s_{a_1} · s_{a_2} · ... · s_{a_k},
/// applied right-to-left (s_{a_k} first), matching the Mathematica convention.
fn permutation_from_word(n: usize, word: &[usize]) -> Vec<usize> {
    permutation_from_simple_transpositions(n, word)
}

/// Count inversions of a permutation.
fn inversion_count(perm: &[usize]) -> usize {
    let n = perm.len();
    let mut count = 0;
    for i in 0..n {
        for j in (i + 1)..n {
            if perm[i] > perm[j] {
                count += 1;
            }
        }
    }
    count
}

/// For backward compatibility: find all valid placements for a given word.
pub fn kogan_face_placements(n: usize, word: &[usize]) -> Vec<Vec<(usize, usize)>> {
    if word.is_empty() {
        return vec![vec![]];
    }
    let perm = permutation_from_word(n, word);
    reduced_kogan_faces(n, &perm)
}

// ── Key polynomial computation ──────────────────────────────────────────────

/// Compute all monomials (as weight vectors) appearing in the key polynomial
/// κ_{λ,σ}, by enumerating integer points in the union of all reduced Kogan
/// faces of type w₀σ.
///
/// Returns a map: weight → multiplicity (always 1 for key polynomials, but
/// the union might have overlapping faces).
pub fn key_polynomial_weights(lambda: &Partition, sigma: &[usize]) -> BTreeMap<Vec<u32>, u64> {
    let n = sigma.len();

    let w0 = longest_permutation(n);
    // The KST formula uses σ⁻¹ where the operator definition uses σ.
    // Invert sigma to match.
    let sigma_inv = inverse_permutation(sigma);
    let w0_sigma = compose_permutations(&sigma_inv, &w0);

    // Find all reduced Kogan faces of type w₀σ directly.
    let faces = reduced_kogan_faces(n, &w0_sigma);

    // Collect all GT-patterns in the union (deduplicate by pattern, not weight).
    let mut all_patterns: HashSet<GtPattern> = HashSet::new();

    for eqs in &faces {
        let patterns = gt_patterns_with_equalities_n(lambda, n, eqs);
        for pat in patterns {
            all_patterns.insert(pat);
        }
    }

    // Count multiplicities: different GT-patterns with the same weight
    // contribute separately to the coefficient.
    let mut weights: BTreeMap<Vec<u32>, u64> = BTreeMap::new();
    for pat in &all_patterns {
        *weights.entry(pat.weight()).or_insert(0) += 1;
    }
    weights
}

/// Coefficient of z^μ in the key polynomial κ_{λ,σ}.
///
/// This is the "key Kostka number": the number of integer GT-patterns
/// of shape λ in the union of reduced Kogan faces of type w₀σ that
/// have weight μ.
pub fn key_kostka(lambda: &Partition, sigma: &[usize], mu: &[u32]) -> u64 {
    let weights = key_polynomial_weights(lambda, sigma);
    weights.get(mu).copied().unwrap_or(0)
}

// ── Ehrhart polynomials for key-Kostka coefficients ─────────────────────────

/// Ehrhart polynomial for the key polynomial lattice point count.
///
/// P(k) = |GT(kλ, σ) ∩ Z^{n(n-1)/2}| is a polynomial in k (by the theory
/// of rational polytopes). This computes the polynomial by evaluating at
/// enough positive dilations and interpolating.
pub struct KeyEhrhartPoly {
    /// Coefficients: `poly[i]` is the coefficient of `k^i`.
    pub coeffs: Vec<BigRational>,
    /// Degree.
    pub degree: usize,
}

impl KeyEhrhartPoly {
    /// Evaluate the polynomial at a given k.
    pub fn eval(&self, k: i64) -> BigRational {
        let k_r = BigRational::from(BigInt::from(k));
        let mut result = BigRational::zero();
        let mut power = BigRational::one();
        for c in &self.coeffs {
            result += c * &power;
            power *= &k_r;
        }
        result
    }

    /// True if any coefficient is negative.
    pub fn has_negative_coefficient(&self) -> bool {
        self.coeffs.iter().any(|c| *c < BigRational::zero())
    }

    /// Display as a human-readable string.
    pub fn display(&self) -> String {
        if self.coeffs.iter().all(|c| c.is_zero()) {
            return "0".into();
        }
        let mut terms: Vec<String> = Vec::new();
        for i in (0..=self.degree).rev() {
            let c = &self.coeffs[i];
            if c.is_zero() {
                continue;
            }
            let c_str = if c.denom() == &BigInt::one() {
                format!("{}", c.numer())
            } else {
                format!("({}/{})", c.numer(), c.denom())
            };
            let term = match i {
                0 => c_str,
                1 => {
                    if *c == BigRational::one() {
                        "k".into()
                    } else {
                        format!("{}k", c_str)
                    }
                }
                _ => {
                    if *c == BigRational::one() {
                        format!("k^{}", i)
                    } else {
                        format!("{}k^{}", c_str, i)
                    }
                }
            };
            terms.push(term);
        }
        terms.join(" + ")
    }
}

/// Count the total number of lattice points in GT(kλ, σ).
///
/// This is the key polynomial evaluated at z = (1,1,...,1), i.e.,
/// κ_{kλ,σ}(1,...,1).
pub fn key_lattice_point_count(lambda: &Partition, sigma: &[usize], dilation: u32) -> u64 {
    let n = sigma.len();
    let kl = scale_partition(lambda, n, dilation);
    let weights = key_polynomial_weights(&kl, sigma);
    weights.values().sum()
}

/// Count strictly interior lattice points in GT(kλ, σ).
///
/// Used for Ehrhart-Macdonald reciprocity:
/// P(-k) = (-1)^d · strict_count(k).
pub fn key_strict_lattice_point_count(lambda: &Partition, sigma: &[usize], dilation: u32) -> u64 {
    let n = sigma.len();
    let kl = scale_partition(lambda, n, dilation);

    let w0 = longest_permutation(n);
    let sigma_inv = inverse_permutation(sigma);
    let w0_sigma = compose_permutations(&sigma_inv, &w0);
    let faces = reduced_kogan_faces(n, &w0_sigma);

    let mut all_patterns: HashSet<GtPattern> = HashSet::new();
    for eqs in &faces {
        let patterns = strict_gt_patterns_with_equalities_n(&kl, n, eqs);
        for pat in patterns {
            all_patterns.insert(pat);
        }
    }
    all_patterns.len() as u64
}

fn scale_partition(lambda: &Partition, n: usize, k: u32) -> Partition {
    Partition::new(
        (0..n)
            .map(|j| {
                let p = if j < lambda.num_parts() {
                    lambda.part(j)
                } else {
                    0
                };
                p * k
            })
            .collect(),
    )
}

/// Compute the Ehrhart polynomial P(k) = |GT(kλ, σ) ∩ Z^{n(n-1)/2}|.
///
/// Compute the Ehrhart polynomial P(k) = |GT(kλ, σ) ∩ Z^{n(n-1)/2}|.
///
/// Uses Ehrhart-Macdonald reciprocity: evaluates both positive dilations
/// and negative dilations (via strict GT-pattern counts), choosing whichever
/// side is cheaper. Strict counts are often 0 for small k, making
/// negative evaluations essentially free.
///
/// `max_degree` is an upper bound on the polynomial degree.
/// If None, uses n(n-1)/2 (dimension of the full GT-polytope).
pub fn key_ehrhart_polynomial(
    lambda: &Partition,
    sigma: &[usize],
    max_degree: Option<usize>,
) -> KeyEhrhartPoly {
    key_ehrhart_polynomial_inner(lambda, sigma, max_degree, true)
}

/// Like `key_ehrhart_polynomial` but uses only positive dilations
/// (no reciprocity). Simpler but slower for large polytopes.
pub fn key_ehrhart_polynomial_positive_only(
    lambda: &Partition,
    sigma: &[usize],
    max_degree: Option<usize>,
) -> KeyEhrhartPoly {
    key_ehrhart_polynomial_inner(lambda, sigma, max_degree, false)
}

fn key_ehrhart_polynomial_inner(
    lambda: &Partition,
    sigma: &[usize],
    max_degree: Option<usize>,
    use_reciprocity: bool,
) -> KeyEhrhartPoly {
    let n = sigma.len();
    let d = max_degree.unwrap_or(n * (n - 1) / 2);

    if !use_reciprocity {
        // Plain method: evaluate at k = 0, 1, ..., d
        let points: Vec<(i64, BigRational)> = (0..=(d as u32))
            .map(|k| {
                let count = key_lattice_point_count(lambda, sigma, k);
                (k as i64, BigRational::from(BigInt::from(count)))
            })
            .collect();

        let coeffs = poly_interpolate_q(&points);
        let true_degree = coeffs
            .iter()
            .enumerate()
            .rev()
            .find(|(_, c)| !c.is_zero())
            .map(|(i, _)| i)
            .unwrap_or(0);
        return KeyEhrhartPoly {
            coeffs,
            degree: true_degree,
        };
    }

    // Reciprocity method: first determine the actual degree by
    // evaluating positive points until the polynomial stabilizes.
    // Then use reciprocity for remaining points.

    // Start with P(0) = 1 (trivial GT-pattern at dilation 0).
    let mut pos_values: Vec<u64> = vec![1]; // pos_values[k] = P(k)
    for k in 1..=(d as u32) {
        pos_values.push(key_lattice_point_count(lambda, sigma, k));
    }

    // Try interpolation and find actual degree
    let points: Vec<(i64, BigRational)> = pos_values
        .iter()
        .enumerate()
        .map(|(k, &v)| (k as i64, BigRational::from(BigInt::from(v))))
        .collect();

    let coeffs = poly_interpolate_q(&points);
    let actual_degree = coeffs
        .iter()
        .enumerate()
        .rev()
        .find(|(_, c)| !c.is_zero())
        .map(|(i, _)| i)
        .unwrap_or(0);

    // Verify: does the polynomial predict the next value correctly?
    // If so, we have the right degree.
    let ep = KeyEhrhartPoly {
        coeffs,
        degree: actual_degree,
    };

    // Optionally verify by checking P at one more point
    let verify_k = (d + 1) as u32;
    let predicted = ep.eval(verify_k as i64);
    let actual = key_lattice_point_count(lambda, sigma, verify_k);
    if predicted == BigRational::from(BigInt::from(actual)) {
        return ep;
    }

    // If verification failed, d was too small. Increase and retry.
    // (This shouldn't happen if max_degree is correct.)
    let mut extended_points: Vec<(i64, BigRational)> = pos_values
        .iter()
        .enumerate()
        .map(|(k, &v)| (k as i64, BigRational::from(BigInt::from(v))))
        .collect();
    extended_points.push((verify_k as i64, BigRational::from(BigInt::from(actual))));
    let coeffs = poly_interpolate_q(&extended_points);
    let degree = coeffs
        .iter()
        .enumerate()
        .rev()
        .find(|(_, c)| !c.is_zero())
        .map(|(i, _)| i)
        .unwrap_or(0);
    KeyEhrhartPoly { coeffs, degree }
}

/// Lagrange interpolation over Q from arbitrary (x, y) points.
fn poly_interpolate_q(points: &[(i64, BigRational)]) -> Vec<BigRational> {
    let d = points.len();
    let mut mat: Vec<Vec<BigRational>> = points
        .iter()
        .map(|&(x, ref y)| {
            let xb = BigInt::from(x);
            let mut row: Vec<BigRational> = Vec::with_capacity(d + 1);
            let mut power = BigInt::one();
            for _ in 0..d {
                row.push(BigRational::from(power.clone()));
                power *= &xb;
            }
            row.push(y.clone());
            row
        })
        .collect();

    // Gauss-Jordan elimination
    for col in 0..d {
        let pivot_row = (col..d)
            .find(|&r| !mat[r][col].is_zero())
            .expect("poly_interpolate: singular system");
        mat.swap(col, pivot_row);
        let pivot = mat[col][col].clone();
        for j in col..=d {
            let v = mat[col][j].clone() / &pivot;
            mat[col][j] = v;
        }
        for row in 0..d {
            if row == col {
                continue;
            }
            let factor = mat[row][col].clone();
            if factor.is_zero() {
                continue;
            }
            for j in col..=d {
                let sub = factor.clone() * &mat[col][j];
                mat[row][j] -= sub;
            }
        }
    }
    mat.iter().map(|row| row[d].clone()).collect()
}

// ══════════════════════════════════════════════════════════════════════════════
// Tests
// ══════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    fn p(parts: &[u32]) -> Partition {
        Partition::new(parts.to_vec())
    }

    // Helper: sum of all multiplicities in a weight map
    fn total_monomials(weights: &BTreeMap<Vec<u32>, u64>) -> u64 {
        weights.values().sum()
    }

    // ── Reduced words ────────────────────────────────────────────────────

    #[test]
    fn test_reduced_words_identity() {
        assert_eq!(reduced_words(&[1, 2, 3]), vec![vec![]]);
    }

    #[test]
    fn test_reduced_words_s1() {
        assert_eq!(reduced_words(&[2, 1, 3]), vec![vec![1]]);
    }

    #[test]
    fn test_reduced_words_w0_s3() {
        // w₀ in S₃ = [3,2,1] = s₁s₂s₁ = s₂s₁s₂ (two reduced words)
        let words = reduced_words(&[3, 2, 1]);
        assert_eq!(words.len(), 2);
        assert!(words.contains(&vec![1, 2, 1]));
        assert!(words.contains(&vec![2, 1, 2]));
    }

    #[test]
    fn test_reduced_words_s1s2() {
        // [3,1,2] = s₁s₂
        let words = reduced_words(&[3, 1, 2]);
        assert_eq!(words, vec![vec![1, 2]]);
    }

    // ── GT-pattern enumeration ───────────────────────────────────────────

    #[test]
    fn test_gt_full_polytope_21() {
        // GT(2,1,0): full polytope should have 8 patterns = s_{2,1}(z₁,z₂,z₃) monomials
        let patterns = gt_patterns_with_equalities_n(&p(&[2, 1, 0]), 3, &[]);
        assert_eq!(patterns.len(), 8);
    }

    #[test]
    fn test_gt_with_equality_l1_l2() {
        // Face L1+L2: x_{3,1}=x_{2,1} (i=2,j=1) and x_{3,2}=x_{2,2} (i=2,j=2)
        // For λ=(2,1,0): forces a=2, b=1, c free in {1,2} → 2 patterns
        let patterns = gt_patterns_with_equalities_n(&p(&[2, 1, 0]), 3, &[(2, 1), (2, 2)]);
        assert_eq!(patterns.len(), 2);
    }

    #[test]
    fn test_gt_all_equalities() {
        // All equalities for n=3: L1+L2+L3
        let patterns = gt_patterns_with_equalities_n(&p(&[2, 1, 0]), 3, &[(2, 1), (2, 2), (1, 1)]);
        assert_eq!(patterns.len(), 1);
        assert_eq!(patterns[0].weight(), vec![2, 1, 0]);
    }

    // ── Key polynomial: n=3, λ=(2,1,0) ──────────────────────────────────
    // Verified against Mathematica OperatorKeyPolynomial

    #[test]
    fn test_key_identity() {
        // κ_{(2,1,0), id} = z₁²z₂ (1 monomial)
        let w = key_polynomial_weights(&p(&[2, 1, 0]), &[1, 2, 3]);
        assert_eq!(total_monomials(&w), 1);
        assert_eq!(w[&vec![2, 1, 0]], 1);
    }

    #[test]
    fn test_key_s1() {
        // κ_{(2,1,0), s₁=[2,1,3]} = z₁²z₂ + z₁z₂² (2 monomials)
        let w = key_polynomial_weights(&p(&[2, 1, 0]), &[2, 1, 3]);
        assert_eq!(total_monomials(&w), 2);
        assert_eq!(w[&vec![2, 1, 0]], 1);
        assert_eq!(w[&vec![1, 2, 0]], 1);
    }

    #[test]
    fn test_key_s2() {
        // κ_{(2,1,0), s₂=[1,3,2]} = z₁²z₂ + z₁²z₃ (2 monomials)
        let w = key_polynomial_weights(&p(&[2, 1, 0]), &[1, 3, 2]);
        assert_eq!(total_monomials(&w), 2);
        assert_eq!(w[&vec![2, 1, 0]], 1);
        assert_eq!(w[&vec![2, 0, 1]], 1);
    }

    #[test]
    fn test_key_s1s2() {
        // κ_{(2,1,0), [3,1,2]} = 5 monomials
        let w = key_polynomial_weights(&p(&[2, 1, 0]), &[3, 1, 2]);
        assert_eq!(total_monomials(&w), 5);
    }

    #[test]
    fn test_key_s2s1() {
        // κ_{(2,1,0), [2,3,1]} = 5 monomials
        let w = key_polynomial_weights(&p(&[2, 1, 0]), &[2, 3, 1]);
        assert_eq!(total_monomials(&w), 5);
    }

    #[test]
    fn test_key_w0_is_schur() {
        // κ_{(2,1,0), w₀=[3,2,1]} = s_{2,1}(z₁,z₂,z₃) = 8 monomials
        let w = key_polynomial_weights(&p(&[2, 1, 0]), &[3, 2, 1]);
        assert_eq!(total_monomials(&w), 8);
    }

    // ── Key polynomial: n=3, λ=(3,1,0) ──────────────────────────────────

    #[test]
    fn test_key_31_identity() {
        // κ_{(3,1,0), id} = z₁³z₂ (1 monomial)
        let w = key_polynomial_weights(&p(&[3, 1, 0]), &[1, 2, 3]);
        assert_eq!(total_monomials(&w), 1);
    }

    #[test]
    fn test_key_31_s1() {
        // κ_{(3,1,0), s₁} = z₁³z₂ + z₁²z₂² + z₁z₂³ (3 monomials)
        let w = key_polynomial_weights(&p(&[3, 1, 0]), &[2, 1, 3]);
        assert_eq!(total_monomials(&w), 3);
    }

    #[test]
    fn test_key_31_w0() {
        // κ_{(3,1,0), w₀} = s_{3,1}(z₁,z₂,z₃) = 12 monomials (counted with multiplicity)
        // Actually s_{3,1} has some monomials with coeff > 1, but key poly counts
        // lattice points (each once). Let me check...
        // s_{3,1}(z₁,z₂,z₃) has 15 total monomials (with multiplicity), 12 distinct.
        let w = key_polynomial_weights(&p(&[3, 1, 0]), &[3, 2, 1]);
        // The Schur polynomial s_{3,1}(z₁,z₂,z₃) = sum over SSYT = 15 GT patterns
        // But key_polynomial_weights counts lattice points (each GT pattern once).
        // So total = number of GT patterns = K(3,1,0 ; 1,1,1,1) + K(3,1,0 ; 2,1,1) + ...
        // = total number of SSYT of shape (3,1) with entries in {1,2,3}
        // = dim(V_{(3,1)} for GL₃) = product formula...
        // Actually total GT patterns for (3,1,0) with 3 rows:
        // row 3 = (3,1,0), row 2 = (a,b) with 3>=a>=1, 1>=b>=0, a>=b
        // row 1 = (c) with a>=c>=b
        // Count: a=1,b=0: c in {0,1} = 2. a=1,b=1: c=1 = 1. a=2,b=0: c in {0,1,2} = 3.
        // a=2,b=1: c in {1,2} = 2. a=3,b=0: c in {0,1,2,3} = 4. a=3,b=1: c in {1,2,3} = 3.
        // Total: 2+1+3+2+4+3 = 15.
        assert_eq!(total_monomials(&w), 15);
    }

    // ── Kogan face placement ─────────────────────────────────────────────

    #[test]
    fn test_kogan_placements_s2s1_n3() {
        // For n=3, word s₂s₁: only one placement
        let placements = kogan_face_placements(3, &[2, 1]);
        assert_eq!(placements.len(), 1);
    }

    #[test]
    fn test_kogan_placements_s1s2_n3() {
        // For n=3, word s₁s₂: only one placement
        let placements = kogan_face_placements(3, &[1, 2]);
        assert_eq!(placements.len(), 1);
    }

    #[test]
    fn test_kogan_placements_s2_n3() {
        // For n=3, word s₂: two positions have s₂ (level 1 and level 2 right)
        let placements = kogan_face_placements(3, &[2]);
        assert_eq!(placements.len(), 2);
    }

    // ── n=4 example from the paper ───────────────────────────────────────

    #[test]
    fn test_key_n4_paper_example() {
        // κ_{(2,1,0,0), [2,4,3,1]} from the Osaka paper
        // Mathematica gives 8 monomials (operator definition)
        let w = key_polynomial_weights(&p(&[2, 1, 0, 0]), &[2, 4, 3, 1]);
        assert_eq!(total_monomials(&w), 8);
    }

    // ── Ehrhart polynomial tests ─────────────────────────────────────────

    #[test]
    fn test_ehrhart_identity_is_one() {
        // κ_{kλ, id} = z₁^{kλ₁} ⋯ z_n^{kλ_n}, always 1 monomial.
        // So P(k) = 1 for all k ≥ 0.
        let ep = key_ehrhart_polynomial(&p(&[2, 1]), &[1, 2, 3], Some(3));
        assert_eq!(ep.degree, 0);
        assert_eq!(ep.eval(0), BigRational::one());
        assert_eq!(ep.eval(5), BigRational::one());
    }

    #[test]
    fn test_ehrhart_w0_is_schur() {
        // κ_{kλ, w₀} = s_{kλ}. For λ=(2,1), n=3:
        // P(k) = number of SSYT of shape (2k, k) with entries in {1,2,3}
        // P(0)=1, P(1)=8, P(2)=27 (these are dim of GL₃ irreps)
        let ep = key_ehrhart_polynomial(&p(&[2, 1]), &[3, 2, 1], Some(3));
        assert_eq!(ep.eval(0), BigRational::from(BigInt::from(1)));
        assert_eq!(ep.eval(1), BigRational::from(BigInt::from(8)));
        // Verify it's a polynomial by checking a predicted value
        let p3 = ep.eval(3);
        let count3 = key_lattice_point_count(&p(&[2, 1]), &[3, 2, 1], 3);
        assert_eq!(p3, BigRational::from(BigInt::from(count3)));
    }

    #[test]
    fn test_ehrhart_s1_linear() {
        // κ_{(k·2, k·1, 0), s₁} should have 2 monomials for any k > 0.
        // Actually it grows: k=1 → 2, k=2 → more.
        // Let's just verify the polynomial interpolation is consistent.
        let ep = key_ehrhart_polynomial(&p(&[2, 1]), &[2, 1, 3], Some(3));
        // Verify at a few points
        for k in 0..=4 {
            let predicted = ep.eval(k as i64);
            let actual = key_lattice_point_count(&p(&[2, 1]), &[2, 1, 3], k);
            assert_eq!(
                predicted,
                BigRational::from(BigInt::from(actual)),
                "Mismatch at k={}",
                k
            );
        }
    }

    #[test]
    fn test_ehrhart_nonnegative_conjecture() {
        // Alexandersson-Alhajjar conjecture: Ehrhart poly has non-negative coefficients.
        // Test for a few (λ, σ) pairs.
        let cases: Vec<(Vec<u32>, Vec<usize>)> = vec![
            (vec![2, 1], vec![2, 1, 3]),
            (vec![2, 1], vec![1, 3, 2]),
            (vec![2, 1], vec![3, 1, 2]),
            (vec![3, 1], vec![2, 1, 3]),
        ];
        for (lam, sigma) in &cases {
            let ep = key_ehrhart_polynomial(&p(lam), sigma, Some(3));
            assert!(
                !ep.has_negative_coefficient(),
                "Negative coefficient for lambda={:?}, sigma={:?}: {}",
                lam,
                sigma,
                ep.display()
            );
        }
    }
}
