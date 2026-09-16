//! QSym power sum bases.
//!
//! This module implements the type 1 and type 2 quasisymmetric power sums
//! indexed by compositions, as defined in:
//!
//! > Ballantine, Daugherty, Hicks, Mason, Niese.
//! > *Quasisymmetric Power Sums.* J. Combin. Theory Ser. A (2020).
//! > <https://doi.org/10.1016/j.jcta.2020.105273>
//!
//! Both are defined via the same structural formula over coarsenings of α,
//! differing only in the block-weight function:
//!
//! ```text
//! Ψ_α = z_α  Σ_{β coarsens α}  1 / Π_B π(B)   · M_{sum(β)}
//! Φ_α = z_α  Σ_{β coarsens α}  1 / Π_B sp(B)  · M_{sum(β)}
//! ```
//!
//! where:
//! - β ranges over all coarsenings of α (merging adjacent parts)
//! - π(B) = Π_{i=1}^m (b_1 + ... + b_i) — product of partial sums
//! - sp(B) = m! · Π_{i=1}^m b_i — factorial of block length times product of parts
//! - z_α = Π i^{m_i} · m_i! — standard z-coefficient
//!
//! Both Ψ ↔ M and Φ ↔ M conversions involve rational coefficients.
//!
//! It also implements the combinatorial power sums `p_alpha` and reverse
//! combinatorial power sums `p^r_alpha` of:
//!
//! > Aliniaeifard, Wang, van Willigenburg.
//! > *P-partition power sums.* European J. Combin. (2023).
//! > <https://doi.org/10.1016/j.ejc.2023.103688>
//!
//! Their monomial expansions are integral and nonnegative.  We use
//! Theorem 5.12 and Remark 5.18 of that paper, counting the stated matrices
//! directly.

use std::collections::BTreeMap;

use num_bigint::BigInt;
use num_rational::BigRational;
use sym_poly_core::{Composition, Ring};

use crate::basis::QSymBasis;
use crate::qsym_function::QSymFunction;

// =========================================================================
// Public API
// =========================================================================

/// Compute Ψ_α (type 1 power sum) expanded in the monomial basis.
///
/// Requires rational coefficients since the expansion involves 1/π(B) factors.
pub fn psi_in_monomial_basis<C: Ring>(alpha: &Composition) -> QSymFunction<C> {
    power_sum_in_monomial_basis::<C>(alpha, BlockWeight::Pi)
}

/// Compute Ψ̃_α = Ψ_α / z_α (normalized type 1 power sum) in the monomial basis.
///
/// Ψ̃_α = Σ_{β coarsens α} (1 / Π_B π(B)) · M_{sum(β)}
///
/// The normalized basis has the property that for naturally labeled posets P,
/// ω(Γ(P)) is Ψ̃-positive with non-negative integer coefficients.
pub fn psi_normalized_in_monomial_basis<C: Ring>(alpha: &Composition) -> QSymFunction<C> {
    normalized_power_sum_in_monomial_basis::<C>(alpha, BlockWeight::Pi)
}

/// Compute Φ̃_α = Φ_α / z_α (normalized type 2 power sum) in the monomial basis.
pub fn phi_normalized_in_monomial_basis<C: Ring>(alpha: &Composition) -> QSymFunction<C> {
    normalized_power_sum_in_monomial_basis::<C>(alpha, BlockWeight::Sp)
}

/// Compute Φ_α (type 2 power sum) expanded in the monomial basis.
///
/// Requires rational coefficients since the expansion involves 1/sp(B) factors.
pub fn phi_in_monomial_basis<C: Ring>(alpha: &Composition) -> QSymFunction<C> {
    power_sum_in_monomial_basis::<C>(alpha, BlockWeight::Sp)
}

/// Convert a QSymFunction from Ψ basis to M basis.
pub fn psi_to_monomial<C: Ring>(f: &QSymFunction<C>) -> QSymFunction<C> {
    expand_to_monomial(f, BlockWeight::Pi)
}

/// Convert a QSymFunction from Φ basis to M basis.
pub fn phi_to_monomial<C: Ring>(f: &QSymFunction<C>) -> QSymFunction<C> {
    expand_to_monomial(f, BlockWeight::Sp)
}

/// Convert a QSymFunction from M basis to Ψ basis.
///
/// Uses per-degree rational matrix inversion.
/// Requires rational coefficients (`Ratio<BigInt>`, `Ratio<i64>`).
pub fn monomial_to_psi<C: Ring>(f: &QSymFunction<C>) -> QSymFunction<C> {
    monomial_to_power_sum(f, QSymBasis::PowerSumPsi, BlockWeight::Pi)
}

/// Convert a QSymFunction from M basis to Φ basis.
///
/// Uses per-degree rational matrix inversion.
/// Requires rational coefficients (`Ratio<BigInt>`, `Ratio<i64>`).
pub fn monomial_to_phi<C: Ring>(f: &QSymFunction<C>) -> QSymFunction<C> {
    monomial_to_power_sum(f, QSymBasis::PowerSumPhi, BlockWeight::Sp)
}

/// Compute the combinatorial power sum `p_alpha` expanded in the monomial basis.
pub fn combinatorial_power_sum_in_monomial_basis<C: Ring>(alpha: &Composition) -> QSymFunction<C> {
    combinatorial_power_sum_in_monomial_basis_with_reverse::<C>(alpha, false)
}

/// Compute the reverse combinatorial power sum `p^r_alpha` expanded in the
/// monomial basis.
pub fn reverse_combinatorial_power_sum_in_monomial_basis<C: Ring>(
    alpha: &Composition,
) -> QSymFunction<C> {
    combinatorial_power_sum_in_monomial_basis_with_reverse::<C>(alpha, true)
}

/// Convert a QSymFunction from the combinatorial power-sum basis to M.
pub fn combinatorial_power_sum_to_monomial<C: Ring>(f: &QSymFunction<C>) -> QSymFunction<C> {
    expand_combinatorial_power_sum_to_monomial(f, false)
}

/// Convert a QSymFunction from the reverse combinatorial power-sum basis to M.
pub fn reverse_combinatorial_power_sum_to_monomial<C: Ring>(
    f: &QSymFunction<C>,
) -> QSymFunction<C> {
    expand_combinatorial_power_sum_to_monomial(f, true)
}

/// Convert from M to the combinatorial power-sum basis.
///
/// The inverse transition matrix may have rational entries; use a rational
/// coefficient ring such as `Ratio<i64>` unless exact divisibility is known.
pub fn monomial_to_combinatorial_power_sum<C: Ring>(f: &QSymFunction<C>) -> QSymFunction<C> {
    monomial_to_combinatorial_power_sum_basis(f, QSymBasis::CombinatorialPowerSum, false)
}

/// Convert from M to the reverse combinatorial power-sum basis.
///
/// The inverse transition matrix may have rational entries; use a rational
/// coefficient ring such as `Ratio<i64>` unless exact divisibility is known.
pub fn monomial_to_reverse_combinatorial_power_sum<C: Ring>(
    f: &QSymFunction<C>,
) -> QSymFunction<C> {
    monomial_to_combinatorial_power_sum_basis(f, QSymBasis::ReverseCombinatorialPowerSum, true)
}

// =========================================================================
// Block-weight functions
// =========================================================================

#[derive(Clone, Copy)]
enum BlockWeight {
    Pi,
    Sp,
}

impl BlockWeight {
    /// Return a factorization of the block weight into individually bounded
    /// factors.  Keeping the factorization avoids narrowing their product to
    /// `i64` in generic arbitrary-precision computations.
    fn factors(self, block: &[u32]) -> Vec<i64> {
        match self {
            BlockWeight::Pi => {
                let mut sum = 0u32;
                block
                    .iter()
                    .map(|&part| {
                        sum = sum.checked_add(part).expect("pi block weight overflow");
                        i64::from(sum)
                    })
                    .collect()
            }
            BlockWeight::Sp => {
                let mut factors: Vec<i64> = (1..=block.len())
                    .map(|factor| i64::try_from(factor).expect("block length overflow"))
                    .collect();
                factors.extend(block.iter().map(|&part| i64::from(part)));
                factors
            }
        }
    }

    #[cfg(test)]
    fn value_i64(self, block: &[u32]) -> i64 {
        self.factors(block).into_iter().fold(1i64, |value, factor| {
            value.checked_mul(factor).unwrap_or_else(|| match self {
                BlockWeight::Pi => panic!("pi block weight overflow"),
                BlockWeight::Sp => panic!("sp block weight overflow"),
            })
        })
    }
}

/// π(B) = product of partial sums of block B.
///
/// For B = (b_1, ..., b_m): π(B) = b_1 · (b_1+b_2) · ... · (b_1+...+b_m).
#[cfg(test)]
fn pi_value(block: &[u32]) -> i64 {
    BlockWeight::Pi.value_i64(block)
}

/// sp(B) = |B|! · Π b_i — factorial of block length times product of parts.
///
/// For B = (b_1, ..., b_m): sp(B) = m! · b_1 · b_2 · ... · b_m.
#[cfg(test)]
fn sp_value(block: &[u32]) -> i64 {
    BlockWeight::Sp.value_i64(block)
}

// =========================================================================
// Generic implementation (shared by Ψ and Φ)
// =========================================================================

/// Normalized power sum: Σ_{β coarsens α} 1/Π_B w(B) · M_{sum(β)} (no z_α factor).
fn normalized_power_sum_in_monomial_basis<C: Ring>(
    alpha: &Composition,
    block_weight: BlockWeight,
) -> QSymFunction<C> {
    if alpha.is_empty() {
        return QSymFunction::basis_element(QSymBasis::Monomial, Composition::empty());
    }

    let parts = alpha.parts();
    let mut terms: BTreeMap<Composition, C> = BTreeMap::new();

    for coarsening in composition_coarsenings(parts) {
        let mut result_parts = Vec::new();
        let mut denominator_factors = Vec::new();

        for block in &coarsening {
            result_parts.push(checked_block_sum(block));
            denominator_factors.extend(block_weight.factors(block));
        }

        let comp = Composition::new(result_parts);
        let coeff = coefficient_from_factors::<C>(Vec::new(), denominator_factors);

        let entry = terms.entry(comp).or_insert_with(C::zero);
        *entry = entry.clone() + coeff;
    }

    QSymFunction::from_terms(QSymBasis::Monomial, terms)
}

/// Generic power sum expansion: z_α · Σ_{β coarsens α} 1/Π_B w(B) · M_{sum(β)}.
fn power_sum_in_monomial_basis<C: Ring>(
    alpha: &Composition,
    block_weight: BlockWeight,
) -> QSymFunction<C> {
    if alpha.is_empty() {
        return QSymFunction::basis_element(QSymBasis::Monomial, Composition::empty());
    }

    let parts = alpha.parts();
    let numerator_factors = z_factors(parts);

    let mut terms: BTreeMap<Composition, C> = BTreeMap::new();

    for coarsening in composition_coarsenings(parts) {
        let mut result_parts = Vec::new();
        let mut denominator_factors = Vec::new();

        for block in &coarsening {
            result_parts.push(checked_block_sum(block));
            denominator_factors.extend(block_weight.factors(block));
        }

        let comp = Composition::new(result_parts);
        let coeff = coefficient_from_factors::<C>(numerator_factors.clone(), denominator_factors);

        let entry = terms.entry(comp).or_insert_with(C::zero);
        *entry = entry.clone() + coeff;
    }

    QSymFunction::from_terms(QSymBasis::Monomial, terms)
}

/// Expand from a power sum basis to monomial basis.
fn expand_to_monomial<C: Ring>(f: &QSymFunction<C>, block_weight: BlockWeight) -> QSymFunction<C> {
    let mut result_terms: BTreeMap<Composition, C> = BTreeMap::new();

    for (alpha, coeff) in f.terms() {
        let expansion = power_sum_in_monomial_basis::<C>(alpha, block_weight);
        for (beta, c) in expansion.terms() {
            let entry = result_terms.entry(beta.clone()).or_insert_with(C::zero);
            *entry = entry.clone() + coeff.clone() * c.clone();
        }
    }

    QSymFunction::from_terms(QSymBasis::Monomial, result_terms)
}

/// Convert from monomial to a power sum basis via per-degree matrix inversion.
fn monomial_to_power_sum<C: Ring>(
    f: &QSymFunction<C>,
    target_basis: QSymBasis,
    block_weight: BlockWeight,
) -> QSymFunction<C> {
    if f.is_zero() {
        return QSymFunction::zero(target_basis);
    }

    let mut by_degree: BTreeMap<u32, Vec<(Composition, C)>> = BTreeMap::new();
    for (comp, coeff) in f.terms() {
        by_degree
            .entry(comp.size())
            .or_default()
            .push((comp.clone(), coeff.clone()));
    }

    let mut result_terms: BTreeMap<Composition, C> = BTreeMap::new();

    for (deg, terms) in by_degree {
        let compositions = Composition::integer_compositions(deg);
        let k = compositions.len();
        let comp_index: BTreeMap<&Composition, usize> = compositions
            .iter()
            .enumerate()
            .map(|(i, c)| (c, i))
            .collect();

        // Build power-sum → M transition matrix as (num, den) pairs
        let ps_to_m = build_transition_matrix(&compositions, &comp_index, block_weight);

        // Invert to get M → power-sum matrix
        let m_to_ps = invert_rational_matrix(&ps_to_m);

        // Build input vector
        let mut input = vec![C::zero(); k];
        for (comp, coeff) in &terms {
            if let Some(&idx) = comp_index.get(comp) {
                input[idx] = coeff.clone();
            }
        }

        // Matrix-vector multiply
        for i in 0..k {
            let mut sum = C::zero();
            for j in 0..k {
                if !input[j].is_zero() {
                    let coefficient = &m_to_ps[j][i];
                    if *coefficient.numer() != BigInt::from(0) {
                        let term = input[j].clone() * C::from_bigint(coefficient.numer());
                        sum = sum + term.exact_div_bigint(coefficient.denom());
                    }
                }
            }
            if !sum.is_zero() {
                result_terms.insert(compositions[i].clone(), sum);
            }
        }
    }

    QSymFunction::from_terms(target_basis, result_terms)
}

fn combinatorial_power_sum_in_monomial_basis_with_reverse<C: Ring>(
    alpha: &Composition,
    reverse: bool,
) -> QSymFunction<C> {
    if alpha.is_empty() {
        return QSymFunction::basis_element(QSymBasis::Monomial, Composition::empty());
    }

    let mut terms = BTreeMap::new();
    for beta in Composition::integer_compositions(alpha.size()) {
        let coeff = combinatorial_power_sum_monomial_coefficient(alpha, &beta, reverse);
        if coeff != 0 {
            terms.insert(beta, C::from_i64(coeff));
        }
    }

    QSymFunction::from_terms(QSymBasis::Monomial, terms)
}

fn expand_combinatorial_power_sum_to_monomial<C: Ring>(
    f: &QSymFunction<C>,
    reverse: bool,
) -> QSymFunction<C> {
    let mut result_terms: BTreeMap<Composition, C> = BTreeMap::new();

    for (alpha, coeff) in f.terms() {
        let expansion = combinatorial_power_sum_in_monomial_basis_with_reverse::<C>(alpha, reverse);
        for (beta, c) in expansion.terms() {
            let entry = result_terms.entry(beta.clone()).or_insert_with(C::zero);
            *entry = entry.clone() + coeff.clone() * c.clone();
        }
    }

    QSymFunction::from_terms(QSymBasis::Monomial, result_terms)
}

fn monomial_to_combinatorial_power_sum_basis<C: Ring>(
    f: &QSymFunction<C>,
    target_basis: QSymBasis,
    reverse: bool,
) -> QSymFunction<C> {
    if f.is_zero() {
        return QSymFunction::zero(target_basis);
    }

    let mut by_degree: BTreeMap<u32, Vec<(Composition, C)>> = BTreeMap::new();
    for (comp, coeff) in f.terms() {
        by_degree
            .entry(comp.size())
            .or_default()
            .push((comp.clone(), coeff.clone()));
    }

    let mut result_terms: BTreeMap<Composition, C> = BTreeMap::new();

    for (deg, terms) in by_degree {
        let compositions = Composition::integer_compositions(deg);
        let k = compositions.len();
        let comp_index: BTreeMap<&Composition, usize> = compositions
            .iter()
            .enumerate()
            .map(|(i, c)| (c, i))
            .collect();

        let ps_to_m =
            build_combinatorial_power_sum_transition_matrix(&compositions, &comp_index, reverse);
        let m_to_ps = invert_rational_matrix(&ps_to_m);

        let mut input = vec![C::zero(); k];
        for (comp, coeff) in &terms {
            if let Some(&idx) = comp_index.get(comp) {
                input[idx] = coeff.clone();
            }
        }

        for i in 0..k {
            let mut sum = C::zero();
            for j in 0..k {
                if !input[j].is_zero() {
                    let coefficient = &m_to_ps[j][i];
                    if *coefficient.numer() != BigInt::from(0) {
                        let term = input[j].clone() * C::from_bigint(coefficient.numer());
                        sum = sum + term.exact_div_bigint(coefficient.denom());
                    }
                }
            }
            if !sum.is_zero() {
                result_terms.insert(compositions[i].clone(), sum);
            }
        }
    }

    QSymFunction::from_terms(target_basis, result_terms)
}

// =========================================================================
// Combinatorial helpers
// =========================================================================

/// z_α = Π i^{m_i} · m_i! where m_i = #{j : α_j = i}.
#[cfg(test)]
fn z_coefficient(parts: &[u32]) -> i64 {
    let mut counts: BTreeMap<u32, u32> = BTreeMap::new();
    for &p in parts {
        let count = counts.entry(p).or_insert(0);
        *count = count.checked_add(1).expect("z coefficient overflow");
    }
    let mut z: i64 = 1;
    for (&val, &mult) in &counts {
        for _ in 0..mult {
            z = z
                .checked_mul(i64::from(val))
                .expect("z coefficient overflow");
        }
        for k in 1..=i64::from(mult) {
            z = z.checked_mul(k).expect("z coefficient overflow");
        }
    }
    z
}

fn z_factors(parts: &[u32]) -> Vec<i64> {
    let mut counts: BTreeMap<u32, usize> = BTreeMap::new();
    for &part in parts {
        *counts.entry(part).or_insert(0) += 1;
    }
    let mut factors = Vec::new();
    for (part, multiplicity) in counts {
        factors.extend(std::iter::repeat(i64::from(part)).take(multiplicity));
        factors.extend(
            (1..=multiplicity)
                .map(|factor| i64::try_from(factor).expect("z coefficient multiplicity overflow")),
        );
    }
    factors
}

fn coefficient_from_factors<C: Ring>(
    mut numerator_factors: Vec<i64>,
    mut denominator_factors: Vec<i64>,
) -> C {
    for denominator in &mut denominator_factors {
        for numerator in &mut numerator_factors {
            let divisor = gcd(numerator.unsigned_abs(), denominator.unsigned_abs()) as i64;
            *numerator /= divisor;
            *denominator /= divisor;
            if *denominator == 1 {
                break;
            }
        }
    }

    let mut coefficient = numerator_factors
        .into_iter()
        .fold(C::one(), |value, factor| value * C::from_i64(factor));
    for denominator in denominator_factors {
        if denominator != 1 {
            coefficient = coefficient.exact_div_i64(denominator);
        }
    }
    coefficient
}

fn checked_block_sum(block: &[u32]) -> u32 {
    block
        .iter()
        .try_fold(0u32, |total, &part| total.checked_add(part))
        .expect("composition block sum overflow")
}

/// All coarsenings of a composition (merging adjacent parts).
///
/// For k parts there are 2^{k-1} coarsenings, one for each subset of
/// the k-1 gaps between consecutive parts.
fn composition_coarsenings(parts: &[u32]) -> Vec<Vec<Vec<u32>>> {
    let k = parts.len();
    if k == 0 {
        return vec![vec![]];
    }
    let num_gaps = k - 1;
    assert!(
        num_gaps < usize::BITS as usize,
        "too many composition gaps to enumerate"
    );
    let mut result = Vec::with_capacity(1usize << num_gaps);

    for mask in 0..(1usize << num_gaps) {
        let mut coarsening = Vec::new();
        let mut block = vec![parts[0]];
        for i in 0..num_gaps {
            if mask & (1usize << i) != 0 {
                block.push(parts[i + 1]);
            } else {
                coarsening.push(block);
                block = vec![parts[i + 1]];
            }
        }
        coarsening.push(block);
        result.push(coarsening);
    }

    result
}

fn combinatorial_power_sum_monomial_coefficient(
    alpha: &Composition,
    beta: &Composition,
    reverse: bool,
) -> i64 {
    if alpha.size() != beta.size() {
        return 0;
    }
    if alpha.is_empty() {
        return if beta.is_empty() { 1 } else { 0 };
    }

    if reverse {
        let alpha_rev = reverse_composition(alpha);
        let beta_rev = reverse_composition(beta);
        return combinatorial_power_sum_monomial_coefficient(&alpha_rev, &beta_rev, false);
    }

    let mut columns = alpha.parts().to_vec();
    columns.sort_by(|left, right| right.cmp(left));

    let mut remaining = beta.parts().to_vec();
    let mut assignment = vec![0usize; columns.len()];
    let mut count = 0i64;
    count_combinatorial_power_sum_matrices(
        alpha.parts(),
        &columns,
        &mut remaining,
        &mut assignment,
        0,
        &mut count,
    );
    count
}

fn count_combinatorial_power_sum_matrices(
    target_word: &[u32],
    columns: &[u32],
    remaining_row_sums: &mut [u32],
    assignment: &mut [usize],
    col: usize,
    count: &mut i64,
) {
    if col == columns.len() {
        if remaining_row_sums.iter().any(|&sum| sum != 0) {
            return;
        }
        if matrix_reading_word(columns, assignment, remaining_row_sums.len()) == target_word {
            *count = count
                .checked_add(1)
                .expect("combinatorial power-sum coefficient overflow");
        }
        return;
    }

    let part = columns[col];
    for row in 0..remaining_row_sums.len() {
        if remaining_row_sums[row] < part {
            continue;
        }
        remaining_row_sums[row] -= part;
        assignment[col] = row;
        count_combinatorial_power_sum_matrices(
            target_word,
            columns,
            remaining_row_sums,
            assignment,
            col + 1,
            count,
        );
        remaining_row_sums[row] += part;
    }
}

fn matrix_reading_word(columns: &[u32], assignment: &[usize], rows: usize) -> Vec<u32> {
    let mut word = Vec::with_capacity(columns.len());
    for row in 0..rows {
        for (col, &part) in columns.iter().enumerate() {
            if assignment[col] == row {
                word.push(part);
            }
        }
    }
    word
}

fn reverse_composition(alpha: &Composition) -> Composition {
    Composition::new(alpha.parts().iter().rev().copied().collect())
}

// =========================================================================
// Rational matrix utilities
// =========================================================================

/// Build the power-sum → M transition matrix as (num, den) pairs.
fn build_transition_matrix(
    compositions: &[Composition],
    comp_index: &BTreeMap<&Composition, usize>,
    block_weight: BlockWeight,
) -> Vec<Vec<BigRational>> {
    let k = compositions.len();
    let mut matrix = vec![vec![BigRational::from_integer(BigInt::from(0)); k]; k];

    for (j, alpha) in compositions.iter().enumerate() {
        let parts = alpha.parts();
        let z = z_factors(parts)
            .into_iter()
            .fold(BigInt::from(1), |value, factor| value * factor);

        for coarsening in composition_coarsenings(parts) {
            let mut result_parts = Vec::new();
            let mut denominator = BigInt::from(1);

            for block in &coarsening {
                result_parts.push(checked_block_sum(block));
                for factor in block_weight.factors(block) {
                    denominator *= factor;
                }
            }

            let comp = Composition::new(result_parts);
            if let Some(&i) = comp_index.get(&comp) {
                matrix[j][i] += BigRational::new(z.clone(), denominator);
            }
        }
    }

    matrix
}

fn build_combinatorial_power_sum_transition_matrix(
    compositions: &[Composition],
    comp_index: &BTreeMap<&Composition, usize>,
    reverse: bool,
) -> Vec<Vec<BigRational>> {
    let k = compositions.len();
    let mut matrix = vec![vec![BigRational::from_integer(BigInt::from(0)); k]; k];

    for (j, alpha) in compositions.iter().enumerate() {
        for beta in Composition::integer_compositions(alpha.size()) {
            let coeff = combinatorial_power_sum_monomial_coefficient(alpha, &beta, reverse);
            if coeff == 0 {
                continue;
            }
            if let Some(&i) = comp_index.get(&beta) {
                matrix[j][i] = BigRational::from_integer(BigInt::from(coeff));
            }
        }
    }

    matrix
}

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

/// Invert a rational matrix via exact Gaussian elimination.
fn invert_rational_matrix(m: &[Vec<BigRational>]) -> Vec<Vec<BigRational>> {
    let n = m.len();
    if n == 0 {
        return vec![];
    }

    let mut aug: Vec<Vec<BigRational>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut row = Vec::with_capacity(2 * n);
        for j in 0..n {
            row.push(m[i][j].clone());
        }
        for j in 0..n {
            row.push(if i == j {
                BigRational::from_integer(BigInt::from(1))
            } else {
                BigRational::from_integer(BigInt::from(0))
            });
        }
        aug.push(row);
    }

    for col in 0..n {
        let mut pivot = None;
        for row in col..n {
            if aug[row][col] != BigRational::from_integer(BigInt::from(0)) {
                pivot = Some(row);
                break;
            }
        }
        let pivot = pivot.expect("power-sum → M matrix is singular");
        aug.swap(col, pivot);

        let diag = aug[col][col].clone();
        for j in 0..2 * n {
            aug[col][j] = aug[col][j].clone() / diag.clone();
        }

        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = aug[row][col].clone();
            for j in 0..2 * n {
                let val = aug[col][j].clone() * factor.clone();
                aug[row][j] = aug[row][j].clone() - val;
            }
        }
    }

    let mut inv = vec![vec![BigRational::from_integer(BigInt::from(0)); n]; n];
    for i in 0..n {
        for j in 0..n {
            inv[i][j] = aug[i][n + j].clone();
        }
    }
    inv
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;
    use num_rational::Ratio;
    type Q = Ratio<i64>;

    // -- Helper tests --

    #[test]
    fn test_pi_value() {
        assert_eq!(pi_value(&[2, 1]), 6); // 2 * 3
        assert_eq!(pi_value(&[3]), 3);
        assert_eq!(pi_value(&[1, 1, 1]), 6); // 1 * 2 * 3
    }

    #[test]
    #[should_panic(expected = "pi block weight overflow")]
    fn test_pi_value_rejects_overflow() {
        let _ = pi_value(&[u32::MAX, 1]);
    }

    #[test]
    fn test_sp_value() {
        assert_eq!(sp_value(&[2, 1]), 2 * 2 * 1); // 2! * 2 * 1 = 4
        assert_eq!(sp_value(&[3]), 1 * 3); // 1! * 3 = 3
        assert_eq!(sp_value(&[1, 1, 1]), 6 * 1 * 1 * 1); // 3! * 1 * 1 * 1 = 6
    }

    #[test]
    #[should_panic(expected = "sp block weight overflow")]
    fn test_sp_value_rejects_overflow() {
        let _ = sp_value(&vec![1; 21]);
    }

    #[test]
    fn test_z_coefficient() {
        assert_eq!(z_coefficient(&[2, 1]), 2);
        assert_eq!(z_coefficient(&[1, 1, 1]), 6);
        assert_eq!(z_coefficient(&[3]), 3);
        assert_eq!(z_coefficient(&[2, 2]), 8);
    }

    #[test]
    #[should_panic(expected = "z coefficient overflow")]
    fn test_z_coefficient_rejects_overflow() {
        let _ = z_coefficient(&vec![1; 21]);
    }

    #[test]
    fn test_factorized_power_sum_coefficient_avoids_i64_overflow() {
        type BigQ = Ratio<BigInt>;

        let parts = vec![1; 21];
        let z = z_factors(&parts);
        let psi = coefficient_from_factors::<BigQ>(z.clone(), BlockWeight::Pi.factors(&parts));
        let phi = coefficient_from_factors::<BigQ>(z, BlockWeight::Sp.factors(&parts));
        assert_eq!(psi, BigQ::from_integer(BigInt::from(1)));
        assert_eq!(phi, BigQ::from_integer(BigInt::from(1)));
    }

    #[test]
    fn test_rational_matrix_inverse_keeps_bigint_coefficients() {
        type BigQ = Ratio<BigInt>;

        let large = BigInt::from(i64::MAX) + BigInt::from(1);
        let matrix = vec![vec![BigQ::new(BigInt::from(1), large.clone())]];
        let inverse = invert_rational_matrix(&matrix);
        assert_eq!(inverse[0][0], BigQ::from_integer(large));
    }

    #[test]
    fn test_coarsenings_count() {
        assert_eq!(composition_coarsenings(&[2, 1]).len(), 2);
        assert_eq!(composition_coarsenings(&[1, 1, 1]).len(), 4);
        assert_eq!(composition_coarsenings(&[3]).len(), 1);
    }

    #[test]
    #[should_panic(expected = "too many composition gaps to enumerate")]
    fn test_coarsenings_reject_too_many_gaps() {
        let parts = vec![1; usize::BITS as usize + 1];

        let _ = composition_coarsenings(&parts);
    }

    fn comp(parts: &[u32]) -> Composition {
        Composition::new(parts.to_vec())
    }

    fn monomial_from_i64_terms(terms: &[(&[u32], i64)]) -> QSymFunction<i64> {
        let mut map = BTreeMap::new();
        for &(parts, coeff) in terms {
            map.insert(comp(parts), coeff);
        }
        QSymFunction::from_terms(QSymBasis::Monomial, map)
    }

    // -- P-partition combinatorial power-sum tests --

    #[test]
    fn test_combinatorial_power_sum_paper_degree_four_table() {
        // Aliniaeifard--Wang--van Willigenburg, Example 5.2.
        let examples: Vec<(&[u32], Vec<(&[u32], i64)>)> = vec![
            (&[4], vec![(&[4], 1)]),
            (&[3, 1], vec![(&[4], 1), (&[3, 1], 1)]),
            (&[1, 3], vec![(&[1, 3], 1)]),
            (
                &[2, 1, 1],
                vec![(&[4], 1), (&[3, 1], 2), (&[2, 2], 1), (&[2, 1, 1], 2)],
            ),
            (&[2, 2], vec![(&[4], 1), (&[2, 2], 2)]),
            (&[1, 2, 1], vec![(&[1, 3], 2), (&[1, 2, 1], 2)]),
            (&[1, 1, 2], vec![(&[2, 2], 1), (&[1, 1, 2], 2)]),
            (
                &[1, 1, 1, 1],
                vec![
                    (&[4], 1),
                    (&[3, 1], 4),
                    (&[1, 3], 4),
                    (&[2, 2], 6),
                    (&[2, 1, 1], 12),
                    (&[1, 2, 1], 12),
                    (&[1, 1, 2], 12),
                    (&[1, 1, 1, 1], 24),
                ],
            ),
        ];

        for (alpha, terms) in examples {
            let actual = combinatorial_power_sum_in_monomial_basis::<i64>(&comp(alpha));
            let expected = monomial_from_i64_terms(&terms);
            assert_eq!(actual, expected, "p_{alpha:?} monomial expansion");
        }
    }

    #[test]
    fn test_reverse_combinatorial_power_sum_monomial_formula() {
        // Remark 5.18 gives [M_beta] p^r_alpha = R_{alpha^r,beta^r}.
        let actual = reverse_combinatorial_power_sum_in_monomial_basis::<i64>(&comp(&[1, 1, 2]));
        let expected =
            monomial_from_i64_terms(&[(&[4], 1), (&[1, 3], 2), (&[2, 2], 1), (&[1, 1, 2], 2)]);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_combinatorial_power_sum_partition_refinement_example() {
        // AWW Example 5.5: the symmetric p_(2,1,1) refines as the sum over
        // all rearrangements of (2,1,1).
        let sum = combinatorial_power_sum_in_monomial_basis::<i64>(&comp(&[2, 1, 1]))
            + combinatorial_power_sum_in_monomial_basis::<i64>(&comp(&[1, 2, 1]))
            + combinatorial_power_sum_in_monomial_basis::<i64>(&comp(&[1, 1, 2]));
        let expected = monomial_from_i64_terms(&[
            (&[4], 1),
            (&[3, 1], 2),
            (&[1, 3], 2),
            (&[2, 2], 2),
            (&[2, 1, 1], 2),
            (&[1, 2, 1], 2),
            (&[1, 1, 2], 2),
        ]);
        assert_eq!(sum, expected);
    }

    #[test]
    fn test_combinatorial_power_sum_product_example_57() {
        // AWW Example 5.7:
        // p_(1,2) p_(1) = (1/2) p_(1,2,1) + p_(1,1,2).
        let lhs = QSymFunction::<Q>::combinatorial_power_sum(comp(&[1, 2]))
            .multiply(&QSymFunction::<Q>::combinatorial_power_sum(comp(&[1])));
        let rhs = QSymFunction::scaled_basis_element(
            QSymBasis::CombinatorialPowerSum,
            comp(&[1, 2, 1]),
            Q::new(1, 2),
        ) + QSymFunction::combinatorial_power_sum(comp(&[1, 1, 2]));
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn test_combinatorial_power_sum_involution_examples_511() {
        // AWW Example 5.11:
        // psi(p_112) = -p^r_112, rho(p_112) = p^r_211,
        // omega(p_112) = -p_211.
        let p112 = combinatorial_power_sum_in_monomial_basis::<i64>(&comp(&[1, 1, 2]));
        let pr112 = reverse_combinatorial_power_sum_in_monomial_basis::<i64>(&comp(&[1, 1, 2]));
        let pr211 = reverse_combinatorial_power_sum_in_monomial_basis::<i64>(&comp(&[2, 1, 1]));
        let p211 = combinatorial_power_sum_in_monomial_basis::<i64>(&comp(&[2, 1, 1]));

        assert_eq!(
            p112.clone().psi_involution().to_monomial_basis(),
            pr112.scale(&-1)
        );
        assert_eq!(p112.clone().rho_involution().to_monomial_basis(), pr211);
        assert_eq!(p112.omega_involution().to_monomial_basis(), p211.scale(&-1));
    }

    #[test]
    fn test_combinatorial_power_sum_basis_roundtrips_degree_four() {
        for alpha in Composition::integer_compositions(4) {
            let p: QSymFunction<Q> =
                QSymFunction::basis_element(QSymBasis::CombinatorialPowerSum, alpha.clone());
            let back = p.to_monomial_basis().to_combinatorial_power_sum_basis();
            assert_eq!(back, p, "roundtrip failed for p_{}", alpha);

            let pr: QSymFunction<Q> =
                QSymFunction::basis_element(QSymBasis::ReverseCombinatorialPowerSum, alpha.clone());
            let back = pr
                .to_monomial_basis()
                .to_reverse_combinatorial_power_sum_basis();
            assert_eq!(back, pr, "roundtrip failed for pr_{}", alpha);
        }
    }

    // -- Ψ tests --

    #[test]
    fn test_psi_single_part() {
        // Ψ_(n) = M_(n)
        let psi3: QSymFunction<Q> = psi_in_monomial_basis(&Composition::new(vec![3]));
        assert_eq!(
            psi3.coefficient(&Composition::new(vec![3])),
            Q::from_integer(1)
        );
        assert_eq!(psi3.terms().len(), 1);
    }

    #[test]
    fn test_psi_21_ne_psi_12() {
        let psi21: QSymFunction<Q> = psi_in_monomial_basis(&Composition::new(vec![2, 1]));
        let psi12: QSymFunction<Q> = psi_in_monomial_basis(&Composition::new(vec![1, 2]));
        assert_ne!(psi21, psi12);

        // Ψ_(2,1) = M_(2,1) + 1/3 M_(3)
        assert_eq!(
            psi21.coefficient(&Composition::new(vec![2, 1])),
            Q::from_integer(1)
        );
        assert_eq!(psi21.coefficient(&Composition::new(vec![3])), Q::new(1, 3));

        // Ψ_(1,2) = M_(1,2) + 2/3 M_(3)
        assert_eq!(
            psi12.coefficient(&Composition::new(vec![1, 2])),
            Q::from_integer(1)
        );
        assert_eq!(psi12.coefficient(&Composition::new(vec![3])), Q::new(2, 3));
    }

    #[test]
    fn test_psi_11() {
        // Ψ_(1,1) = 2*M_(1,1) + M_(2)
        let psi11: QSymFunction<Q> = psi_in_monomial_basis(&Composition::new(vec![1, 1]));
        assert_eq!(
            psi11.coefficient(&Composition::new(vec![1, 1])),
            Q::from_integer(2)
        );
        assert_eq!(
            psi11.coefficient(&Composition::new(vec![2])),
            Q::from_integer(1)
        );
    }

    #[test]
    fn test_psi_111() {
        // Ψ_(1,1,1) = 6*M_(1,1,1) + 3*M_(2,1) + 3*M_(1,2) + M_(3)
        let psi: QSymFunction<Q> = psi_in_monomial_basis(&Composition::new(vec![1, 1, 1]));
        assert_eq!(
            psi.coefficient(&Composition::new(vec![1, 1, 1])),
            Q::from_integer(6)
        );
        assert_eq!(
            psi.coefficient(&Composition::new(vec![2, 1])),
            Q::from_integer(3)
        );
        assert_eq!(
            psi.coefficient(&Composition::new(vec![1, 2])),
            Q::from_integer(3)
        );
        assert_eq!(
            psi.coefficient(&Composition::new(vec![3])),
            Q::from_integer(1)
        );
    }

    #[test]
    fn test_psi_roundtrip_degree3() {
        for alpha in Composition::integer_compositions(3) {
            let m: QSymFunction<Q> = QSymFunction::monomial_qsym(alpha.clone());
            let in_psi = monomial_to_psi(&m);
            let back = psi_to_monomial(&in_psi);
            assert_eq!(
                back.coefficient(&alpha),
                Q::from_integer(1),
                "Ψ roundtrip failed for M_{}",
                alpha
            );
            assert_eq!(back.terms().len(), 1, "Ψ extra terms for M_{}", alpha);
        }
    }

    // -- Φ tests --

    #[test]
    fn test_phi_single_part() {
        // Φ_(n) should also equal M_(n) since the only coarsening of (n) is itself,
        // and both π({n}) = n and sp({n}) = 1!*n = n.
        let phi3: QSymFunction<Q> = phi_in_monomial_basis(&Composition::new(vec![3]));
        assert_eq!(
            phi3.coefficient(&Composition::new(vec![3])),
            Q::from_integer(1)
        );
        assert_eq!(phi3.terms().len(), 1);
    }

    #[test]
    fn test_phi_21_ne_phi_12() {
        let phi21: QSymFunction<Q> = phi_in_monomial_basis(&Composition::new(vec![2, 1]));
        let phi12: QSymFunction<Q> = phi_in_monomial_basis(&Composition::new(vec![1, 2]));
        assert_ne!(phi21, phi12);
    }

    #[test]
    fn test_phi_21() {
        // Φ_(2,1): coarsenings of (2,1):
        //   {(2),(1)}: sp({2})=2, sp({1})=1. denom=2. coeff = z/2 = 2/2 = 1. comp=(2,1).
        //   {(2,1)}:   sp({2,1})=2!*2*1=4. denom=4. coeff = z/4 = 2/4 = 1/2. comp=(3).
        // Φ_(2,1) = M_(2,1) + 1/2 M_(3)
        let phi21: QSymFunction<Q> = phi_in_monomial_basis(&Composition::new(vec![2, 1]));
        assert_eq!(
            phi21.coefficient(&Composition::new(vec![2, 1])),
            Q::from_integer(1)
        );
        assert_eq!(phi21.coefficient(&Composition::new(vec![3])), Q::new(1, 2));
    }

    #[test]
    fn test_phi_12() {
        // Φ_(1,2): coarsenings of (1,2):
        //   {(1),(2)}: sp({1})=1, sp({2})=2. denom=2. coeff = 2/2 = 1. comp=(1,2).
        //   {(1,2)}:   sp({1,2})=2!*1*2=4. denom=4. coeff = 2/4 = 1/2. comp=(3).
        // Φ_(1,2) = M_(1,2) + 1/2 M_(3)
        let phi12: QSymFunction<Q> = phi_in_monomial_basis(&Composition::new(vec![1, 2]));
        assert_eq!(
            phi12.coefficient(&Composition::new(vec![1, 2])),
            Q::from_integer(1)
        );
        assert_eq!(phi12.coefficient(&Composition::new(vec![3])), Q::new(1, 2));
    }

    #[test]
    fn test_phi_roundtrip_degree3() {
        for alpha in Composition::integer_compositions(3) {
            let m: QSymFunction<Q> = QSymFunction::monomial_qsym(alpha.clone());
            let in_phi = monomial_to_phi(&m);
            let back = phi_to_monomial(&in_phi);
            assert_eq!(
                back.coefficient(&alpha),
                Q::from_integer(1),
                "Φ roundtrip failed for M_{}",
                alpha
            );
            assert_eq!(back.terms().len(), 1, "Φ extra terms for M_{}", alpha);
        }
    }

    // -- Cross-check: Ψ ≠ Φ --

    #[test]
    fn test_psi_ne_phi() {
        let psi: QSymFunction<Q> = psi_in_monomial_basis(&Composition::new(vec![2, 1]));
        let phi: QSymFunction<Q> = phi_in_monomial_basis(&Composition::new(vec![2, 1]));
        assert_ne!(psi, phi);
    }

    // -- Omega involution --

    #[test]
    fn test_omega_involution() {
        // ω² = id on a fundamental basis element
        let f: QSymFunction<Q> = QSymFunction::fundamental_qsym(Composition::new(vec![2, 1]));
        let omega_f = f.omega_involution();
        let back = omega_f.omega_involution();
        assert_eq!(back, f, "ω² should be identity");
    }

    #[test]
    fn test_omega_on_psi() {
        // ω(Ψ_α) = (-1)^{n-ℓ(α)} Ψ_{α^r}
        // For α = (2,1), n=3, ℓ=2: sign = (-1)^{3-2} = -1, α^r = (1,2)
        let psi21: QSymFunction<Q> = psi_in_monomial_basis(&Composition::new(vec![2, 1]));
        let omega_psi = psi21.omega_involution();
        let psi12: QSymFunction<Q> = psi_in_monomial_basis(&Composition::new(vec![1, 2]));
        let expected = psi12.scale(&Q::from_integer(-1));
        assert_eq!(
            omega_psi.to_monomial_basis(),
            expected.to_monomial_basis(),
            "ω(Ψ_(2,1)) should be -Ψ_(1,2)"
        );
    }

    #[test]
    fn test_omega_psi_identity_degree3() {
        // Verify ω(Ψ_α) = (-1)^{n-ℓ(α)} Ψ_{α^r} for all compositions of 3
        for alpha in Composition::integer_compositions(3) {
            let psi_alpha: QSymFunction<Q> = psi_in_monomial_basis(&alpha);
            let omega_psi = psi_alpha.omega_involution().to_monomial_basis();

            let n = alpha.size();
            let ell = alpha.num_parts() as u32;
            let sign = if (n - ell) % 2 == 0 {
                Q::from_integer(1)
            } else {
                Q::from_integer(-1)
            };
            let reversed_parts: Vec<u32> = alpha.parts().iter().rev().copied().collect();
            let alpha_r = Composition::new(reversed_parts);
            let psi_rev: QSymFunction<Q> = psi_in_monomial_basis(&alpha_r);
            let expected = psi_rev.scale(&sign).to_monomial_basis();

            assert_eq!(
                omega_psi, expected,
                "ω(Ψ_{}) should be ({})Ψ_{}",
                alpha, sign, alpha_r
            );
        }
    }

    // -- Normalized Ψ̃ = Ψ/z --

    #[test]
    fn test_psi_normalized_21() {
        // Ψ̃_(2,1) = Ψ_(2,1)/z_(2,1) = (M_(2,1) + 1/3 M_(3)) / 2
        //          = 1/2 M_(2,1) + 1/6 M_(3)
        let psi_n: QSymFunction<Q> =
            psi_normalized_in_monomial_basis(&Composition::new(vec![2, 1]));
        assert_eq!(
            psi_n.coefficient(&Composition::new(vec![2, 1])),
            Q::new(1, 2)
        );
        assert_eq!(psi_n.coefficient(&Composition::new(vec![3])), Q::new(1, 6));
    }

    #[test]
    fn test_p_partition_psi_normalized_positive() {
        // Γ(P) is Ψ̃-positive for naturally labeled posets:
        //   Γ(P) = Σ c_α · Ψ̃_α  with c_α ∈ ℤ_{≥0}
        // i.e., the Ψ_α coefficient of Γ(P) times z_α is a non-negative integer.
        let posets: Vec<(usize, Vec<(usize, usize)>)> = vec![
            (3, vec![(0, 1), (1, 2)]),         // chain
            (3, vec![]),                       // antichain
            (3, vec![(0, 2), (1, 2)]),         // V-poset
            (4, vec![(0, 1), (1, 2), (2, 3)]), // chain 4
            (4, vec![(0, 2), (1, 3)]),         // two disjoint edges
            (4, vec![(0, 1), (0, 2), (0, 3)]), // star
        ];
        for (n, covers) in &posets {
            let gamma: QSymFunction<Q> =
                crate::p_partition::p_partition_generating_function(*n, covers);
            let in_psi = gamma.to_monomial_basis().to_basis(QSymBasis::PowerSumPsi);
            for (alpha, coeff) in in_psi.terms() {
                let z = z_coefficient(alpha.parts());
                let scaled = coeff.clone() * Q::from_integer(z);
                assert!(
                    *scaled.denom() == 1,
                    "poset {:?}: z·coeff of Ψ_{} not integer: {}",
                    covers,
                    alpha,
                    scaled
                );
                assert!(
                    *scaled.numer() >= 0,
                    "poset {:?}: z·coeff of Ψ_{} negative: {}",
                    covers,
                    alpha,
                    scaled
                );
            }
        }
    }
}
