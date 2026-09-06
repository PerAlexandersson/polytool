//! Linear algebra over exact rationals (BigInt/BigRational).
//!
//! Provides fraction-free Bareiss elimination, Gaussian elimination, positive
//! definiteness checks, determinants, linear system solving, and total
//! non-negativity (TNN) checking.
//!
//! # Total non-negativity
//!
//! A matrix is *totally non-negative* (TNN) if every minor (determinant of every
//! square submatrix) is non-negative. The brute-force check
//! ([`check_total_positivity`]) enumerates all submatrices up to a given size,
//! which is exponential. The Neville elimination method ([`check_tnn_neville`])
//! checks TNN for *all* minors in O(n²m) time by performing adjacent-row
//! subtraction and verifying that no entry becomes negative.

#![allow(clippy::needless_range_loop)]

use crate::Polynomial;
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, ToPrimitive, Zero};

type Q = Ratio<BigInt>;

/// Matrix dimension at which the default positive-definiteness check switches
/// from BigInt Bareiss to modular CRT reconstruction.
pub const MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD: usize = 30;

/// Check if a symmetric BigInt matrix is positive definite.
///
/// Uses Sylvester's criterion.  Small matrices use fraction-free BigInt
/// Bareiss elimination; matrices of size
/// [`MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD`] or larger try the modular
/// CRT path first and fall back to Bareiss if the CRT bounds do not certify in
/// time. The explicit [`is_positive_definite_bareiss`] and
/// [`is_positive_definite_modular`] entry points are kept for benchmarking.
pub fn is_positive_definite(mat: &[Vec<BigInt>]) -> bool {
    if !is_symmetric_bigint_matrix(mat) {
        return false;
    }
    if mat.len() >= MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD {
        is_positive_definite_modular(mat).unwrap_or_else(|| is_positive_definite_bareiss(mat))
    } else {
        is_positive_definite_bareiss(mat)
    }
}

/// Check positive definiteness via BigInt Bareiss leading principal minors.
pub fn is_positive_definite_bareiss(mat: &[Vec<BigInt>]) -> bool {
    if !is_symmetric_bigint_matrix(mat) {
        return false;
    }
    bareiss_leading_principal_minors_bigint(mat)
        .map(|minors| minors.into_iter().all(|minor| minor > BigInt::zero()))
        .unwrap_or(false)
}

/// Try to check positive definiteness via modular CRT.
///
/// This is a second implementation of Sylvester's criterion.  It reconstructs
/// the leading principal minors from determinants computed over prime fields,
/// using Hadamard bounds to certify when CRT reconstruction is exact.
///
/// Returns `None` when the matrix is not square or when the accumulated CRT
/// modulus does not certify all leading minors within the prime budget.
pub fn is_positive_definite_modular(mat: &[Vec<BigInt>]) -> Option<bool> {
    if !is_square_bigint_matrix(mat) {
        return None;
    }
    if !is_symmetric_bigint_matrix(mat) {
        return Some(false);
    }
    modular_leading_principal_minors_bigint(mat)
        .map(|minors| minors.into_iter().all(|minor| minor > BigInt::zero()))
}

/// Check if a symmetric BigInt matrix is positive semi-definite.
///
/// A symmetric matrix is positive semi-definite iff all diagonal pivots
/// are non-negative, and zero pivots have zero entries in their column below.
pub fn is_positive_semidefinite(mat: &[Vec<BigInt>]) -> bool {
    let n = mat.len();
    if !is_symmetric_bigint_matrix(mat) {
        return false;
    }
    let qmat = bigint_to_q(mat);
    let mut a = qmat;

    for k in 0..n {
        if a[k][k] < Q::zero() {
            return false;
        }
        if a[k][k].is_zero() {
            // All entries below must also be zero
            for i in (k + 1)..n {
                if !a[i][k].is_zero() {
                    return false;
                }
            }
            continue;
        }
        let pivot = a[k][k].clone();
        for i in (k + 1)..n {
            if a[i][k].is_zero() {
                continue;
            }
            let factor = a[i][k].clone() / pivot.clone();
            for j in (k + 1)..n {
                let sub = factor.clone() * a[k][j].clone();
                a[i][j] -= sub;
            }
            a[i][k] = Q::zero();
        }
    }
    true
}

/// Compute the determinant of a BigInt matrix using fraction-free elimination.
pub fn determinant(mat: &[Vec<BigInt>]) -> BigInt {
    assert!(
        is_square_bigint_matrix(mat),
        "determinant requires a square matrix"
    );
    let n = mat.len();
    if n == 0 {
        return BigInt::one();
    }

    let mut a = mat.to_vec();
    let mut denom = BigInt::one();
    let mut sign_is_negative = false;

    for k in 0..(n - 1) {
        if a[k][k].is_zero() {
            let Some(pivot_row) = ((k + 1)..n).find(|&i| !a[i][k].is_zero()) else {
                return BigInt::zero();
            };
            a.swap(k, pivot_row);
            sign_is_negative = !sign_is_negative;
        }

        let pivot = a[k][k].clone();
        let mut next_entries = Vec::with_capacity((n - k - 1) * (n - k - 1));
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                let numerator = &a[i][j] * &pivot - &a[i][k] * &a[k][j];
                debug_assert!(
                    (&numerator % &denom).is_zero(),
                    "Bareiss division should be exact"
                );
                next_entries.push((i, j, numerator / &denom));
            }
        }
        for (i, j, value) in next_entries {
            a[i][j] = value;
        }
        denom = pivot;
    }

    if sign_is_negative {
        -a[n - 1][n - 1].clone()
    } else {
        a[n - 1][n - 1].clone()
    }
}

fn is_square_bigint_matrix(mat: &[Vec<BigInt>]) -> bool {
    let n = mat.len();
    mat.iter().all(|row| row.len() == n)
}

fn is_symmetric_bigint_matrix(mat: &[Vec<BigInt>]) -> bool {
    let n = mat.len();
    if !is_square_bigint_matrix(mat) {
        return false;
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if mat[i][j] != mat[j][i] {
                return false;
            }
        }
    }
    true
}

/// Compute the leading principal determinants by fraction-free Bareiss
/// elimination over `BigInt`, without row swaps.
///
/// If every Bareiss pivot is nonzero, the `k`th returned value is the
/// determinant of the leading `(k+1) x (k+1)` principal submatrix.  This is
/// the exact integer version of the fixed-order elimination used in
/// Sylvester's criterion.
///
/// Returns `None` if the matrix is not square, a nonfinal required pivot is
/// zero, or an exact Bareiss division unexpectedly fails.  A nonfinal zero
/// pivot does not imply that the whole matrix is singular; it only means this
/// no-pivot leading-principal computation cannot continue.
pub fn bareiss_leading_principal_minors_bigint(mat: &[Vec<BigInt>]) -> Option<Vec<BigInt>> {
    let n = mat.len();
    if mat.iter().any(|row| row.len() != n) {
        return None;
    }
    if n == 0 {
        return Some(Vec::new());
    }

    let mut a = mat.to_vec();
    let mut denom = BigInt::one();
    let mut minors = Vec::with_capacity(n);

    for k in 0..n {
        let pivot = a[k][k].clone();
        if pivot.is_zero() {
            if k == n - 1 {
                minors.push(pivot);
                break;
            }
            return None;
        }
        minors.push(pivot.clone());
        if k == n - 1 {
            break;
        }

        let mut next_entries = Vec::with_capacity((n - k - 1) * (n - k - 1));
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                let numerator = &a[i][j] * &pivot - &a[i][k] * &a[k][j];
                if &numerator % &denom != BigInt::zero() {
                    return None;
                }
                next_entries.push((i, j, numerator / &denom));
            }
        }
        for (i, j, value) in next_entries {
            a[i][j] = value;
        }
        denom = pivot;
    }

    Some(minors)
}

/// Compute a determinant by no-pivot fraction-free Bareiss elimination.
///
/// This returns `None` in the same cases as
/// [`bareiss_leading_principal_minors_bigint`].  Use [`determinant`] when a
/// row-swapping determinant routine is needed.
pub fn bareiss_determinant_bigint(mat: &[Vec<BigInt>]) -> Option<BigInt> {
    if mat.is_empty() {
        return Some(BigInt::one());
    }
    bareiss_leading_principal_minors_bigint(mat).and_then(|mut minors| minors.pop())
}

/// Compute leading principal determinants by modular determinant computation
/// and Chinese remaindering.
///
/// For a square integer matrix, this returns the determinant of each leading
/// principal submatrix.  Each determinant is computed modulo a sequence of
/// large primes and reconstructed by CRT once the accumulated modulus exceeds
/// twice a Hadamard bound for every leading minor.
///
/// Unlike [`bareiss_leading_principal_minors_bigint`], this method can handle
/// zero Bareiss pivots because each leading determinant is computed separately
/// with modular pivoting.  It is intended as a comparison/large-entry path, not
/// as a replacement for the simpler Bareiss implementation.
pub fn modular_leading_principal_minors_bigint(mat: &[Vec<BigInt>]) -> Option<Vec<BigInt>> {
    const MAX_MODULAR_PRIMES: usize = 256;

    let n = mat.len();
    if mat.iter().any(|row| row.len() != n) {
        return None;
    }
    if n == 0 {
        return Some(Vec::new());
    }

    let hadamard_squared_bounds = leading_principal_hadamard_squared_bounds(mat)?;
    let mut residues = vec![BigInt::zero(); n];
    let mut modulus = BigInt::one();
    let mut prime_search_start = (1u64 << 61) - 1;

    for _ in 0..MAX_MODULAR_PRIMES {
        let prime = previous_prime_at_or_below(prime_search_start)?;
        prime_search_start = prime.saturating_sub(2);
        let modular_minors = leading_principal_minors_mod_prime(mat, prime);
        crt_update_residues(&mut residues, &mut modulus, &modular_minors, prime)?;

        if crt_modulus_certifies_bounds(&modulus, &hadamard_squared_bounds) {
            return Some(
                residues
                    .into_iter()
                    .map(|r| symmetric_residue(r, &modulus))
                    .collect(),
            );
        }
    }

    None
}

fn leading_principal_hadamard_squared_bounds(mat: &[Vec<BigInt>]) -> Option<Vec<BigInt>> {
    let n = mat.len();
    let mut bounds = Vec::with_capacity(n);
    for k in 1..=n {
        let mut product = BigInt::one();
        for i in 0..k {
            let mut row_square_sum = BigInt::zero();
            for j in 0..k {
                row_square_sum += &mat[i][j] * &mat[i][j];
            }
            product *= row_square_sum;
        }
        bounds.push(product);
    }
    Some(bounds)
}

fn crt_modulus_certifies_bounds(modulus: &BigInt, squared_bounds: &[BigInt]) -> bool {
    let modulus_squared = modulus * modulus;
    squared_bounds.iter().all(|bound| {
        let threshold = BigInt::from(4) * bound;
        modulus_squared > threshold
    })
}

fn symmetric_residue(residue: BigInt, modulus: &BigInt) -> BigInt {
    if &residue * BigInt::from(2) > *modulus {
        residue - modulus
    } else {
        residue
    }
}

fn crt_update_residues(
    residues: &mut [BigInt],
    modulus: &mut BigInt,
    modular_values: &[u64],
    prime: u64,
) -> Option<()> {
    let modulus_mod_prime = bigint_mod_u64(modulus, prime);
    let inverse = mod_inverse_prime(modulus_mod_prime, prime)?;
    let old_modulus = modulus.clone();

    for (residue, &value_mod_prime) in residues.iter_mut().zip(modular_values.iter()) {
        let residue_mod_prime = bigint_mod_u64(residue, prime);
        let correction = if value_mod_prime >= residue_mod_prime {
            value_mod_prime - residue_mod_prime
        } else {
            prime - (residue_mod_prime - value_mod_prime)
        };
        let t = mul_mod_u64(correction, inverse, prime);
        *residue += &old_modulus * BigInt::from(t);
    }

    *modulus = old_modulus * BigInt::from(prime);
    Some(())
}

fn bigint_mod_u64(value: &BigInt, modulus: u64) -> u64 {
    let modulus_big = BigInt::from(modulus);
    let mut reduced = value % &modulus_big;
    if reduced.is_negative() {
        reduced += &modulus_big;
    }
    reduced.to_u64().expect("reduced residue should fit in u64")
}

fn leading_principal_minors_mod_prime(mat: &[Vec<BigInt>], prime: u64) -> Vec<u64> {
    let reduced: Vec<Vec<u64>> = mat
        .iter()
        .map(|row| {
            row.iter()
                .map(|entry| bigint_mod_u64(entry, prime))
                .collect()
        })
        .collect();
    if let Some(minors) = leading_principal_minors_mod_prime_bareiss(&reduced, prime) {
        return minors;
    }

    let n = mat.len();
    (1..=n)
        .map(|prefix_size| determinant_mod_prime_prefix(&reduced, prefix_size, prime))
        .collect()
}

fn leading_principal_minors_mod_prime_bareiss(mat: &[Vec<u64>], prime: u64) -> Option<Vec<u64>> {
    let n = mat.len();
    let mut a = mat.to_vec();
    let mut denom = 1u64;
    let mut minors = Vec::with_capacity(n);

    for k in 0..n {
        let pivot = a[k][k];
        if pivot == 0 {
            if k == n - 1 {
                minors.push(0);
                break;
            }
            return None;
        }
        minors.push(pivot);
        if k == n - 1 {
            break;
        }

        let denom_inverse = mod_inverse_prime(denom, prime)?;
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                let left = mul_mod_u64(a[i][j], pivot, prime);
                let right = mul_mod_u64(a[i][k], a[k][j], prime);
                let numerator = sub_mod_u64(left, right, prime);
                a[i][j] = mul_mod_u64(numerator, denom_inverse, prime);
            }
        }
        denom = pivot;
    }

    Some(minors)
}

fn determinant_mod_prime_prefix(mat: &[Vec<u64>], size: usize, prime: u64) -> u64 {
    let mut a: Vec<Vec<u64>> = mat
        .iter()
        .take(size)
        .map(|row| row.iter().take(size).copied().collect())
        .collect();
    let mut det = 1u64;
    let mut sign_is_negative = false;

    for k in 0..size {
        let Some(pivot_row) = (k..size).find(|&i| a[i][k] != 0) else {
            return 0;
        };
        if pivot_row != k {
            a.swap(k, pivot_row);
            sign_is_negative = !sign_is_negative;
        }
        let pivot = a[k][k];
        det = mul_mod_u64(det, pivot, prime);
        let inverse = mod_inverse_prime(pivot, prime).expect("nonzero pivot modulo prime");
        for i in (k + 1)..size {
            if a[i][k] == 0 {
                continue;
            }
            let factor = mul_mod_u64(a[i][k], inverse, prime);
            for j in (k + 1)..size {
                let sub = mul_mod_u64(factor, a[k][j], prime);
                a[i][j] = sub_mod_u64(a[i][j], sub, prime);
            }
        }
    }

    if sign_is_negative && det != 0 {
        prime - det
    } else {
        det
    }
}

fn mul_mod_u64(a: u64, b: u64, modulus: u64) -> u64 {
    if modulus <= u32::MAX as u64 {
        return (a * b) % modulus;
    }
    ((a as u128 * b as u128) % modulus as u128) as u64
}

fn sub_mod_u64(a: u64, b: u64, modulus: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        modulus - (b - a)
    }
}

fn add_mod_u64(a: u64, b: u64, modulus: u64) -> u64 {
    let (sum, overflow) = a.overflowing_add(b);
    if overflow {
        ((a as u128 + b as u128) % modulus as u128) as u64
    } else if sum >= modulus {
        sum - modulus
    } else {
        sum
    }
}

fn pow_mod_u64(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut result = 1u64;
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = mul_mod_u64(result, base, modulus);
        }
        base = mul_mod_u64(base, base, modulus);
        exponent >>= 1;
    }
    result
}

fn mod_inverse_prime(value: u64, prime: u64) -> Option<u64> {
    (value != 0).then(|| pow_mod_u64(value, prime - 2, prime))
}

fn signed_mod_u64(value: i64, prime: u64) -> u64 {
    (i128::from(value).rem_euclid(i128::from(prime))) as u64
}

/// Errors raised while constructing or updating a sparse modular row.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SparseModRowError {
    ZeroModulus,
}

impl std::fmt::Display for SparseModRowError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ZeroModulus => write!(f, "sparse modular row modulus must be nonzero"),
        }
    }
}

impl std::error::Error for SparseModRowError {}

/// Sparse augmented row over the prime field `F_p`.
///
/// The row represents
///
/// ```text
///   sum_i entries[i].1 * x_{entries[i].0} = rhs   (mod p).
/// ```
///
/// Entries are stored sorted by column and normalized into `0..p`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SparseModRow {
    entries: Vec<(usize, u64)>,
    rhs: u64,
}

impl SparseModRow {
    /// Create an empty sparse row with right-hand side `rhs mod prime`.
    ///
    /// The `prime` argument is not stored; it is used only to normalize input.
    pub fn new(rhs: u64, prime: u64) -> Result<Self, SparseModRowError> {
        if prime == 0 {
            return Err(SparseModRowError::ZeroModulus);
        }
        Ok(Self {
            entries: Vec::new(),
            rhs: rhs % prime,
        })
    }

    pub(crate) fn with_capacity(rhs: u64, prime: u64, capacity: usize) -> Self {
        debug_assert_ne!(prime, 0);
        Self {
            entries: Vec::with_capacity(capacity),
            rhs: rhs % prime,
        }
    }

    /// Create a sparse row from unsigned `(column, value)` entries.
    ///
    /// Duplicate columns are combined modulo `prime`, and zero coefficients are
    /// removed.
    pub fn from_entries<I>(entries: I, rhs: u64, prime: u64) -> Result<Self, SparseModRowError>
    where
        I: IntoIterator<Item = (usize, u64)>,
    {
        let mut row = Self::new(rhs, prime)?;
        for (col, value) in entries {
            row.add_entry_unchecked(col, value, prime);
        }
        Ok(row)
    }

    /// Create an empty sparse row with a signed right-hand side.
    pub fn from_signed_rhs(rhs: i64, prime: u64) -> Result<Self, SparseModRowError> {
        if prime == 0 {
            return Err(SparseModRowError::ZeroModulus);
        }
        Ok(Self {
            entries: Vec::new(),
            rhs: signed_mod_u64(rhs, prime),
        })
    }

    /// Create a sparse row from signed `(column, value)` entries.
    ///
    /// Duplicate columns are combined modulo `prime`, and zero coefficients are
    /// removed.
    pub fn from_signed_entries<I>(
        entries: I,
        rhs: i64,
        prime: u64,
    ) -> Result<Self, SparseModRowError>
    where
        I: IntoIterator<Item = (usize, i64)>,
    {
        let mut row = Self::from_signed_rhs(rhs, prime)?;
        for (col, value) in entries {
            row.add_entry_unchecked(col, signed_mod_u64(value, prime), prime);
        }
        Ok(row)
    }

    /// Return the nonzero column/value entries, sorted by column.
    pub fn entries(&self) -> &[(usize, u64)] {
        &self.entries
    }

    /// Return the right-hand side.
    pub fn rhs(&self) -> u64 {
        self.rhs
    }

    /// Return true if `solution` satisfies this augmented row over `F_p`.
    pub fn is_satisfied_by(&self, solution: &[u64], prime: u64) -> Result<bool, SparseModRowError> {
        if prime == 0 {
            return Err(SparseModRowError::ZeroModulus);
        }
        let mut lhs = 0;
        for &(col, coeff) in &self.entries {
            let Some(&value) = solution.get(col) else {
                return Ok(false);
            };
            lhs = add_mod_u64(lhs, mul_mod_u64(coeff, value, prime), prime);
        }
        Ok(lhs == self.rhs)
    }

    fn is_satisfied_by_validated_solution(&self, solution: &[u64], prime: u64) -> bool {
        let mut lhs = 0;
        for &(col, coeff) in &self.entries {
            debug_assert!(col < solution.len());
            lhs = add_mod_u64(lhs, mul_mod_u64(coeff, solution[col], prime), prime);
        }
        lhs == self.rhs
    }

    /// Add `value` to the coefficient of `col`, reducing modulo `prime`.
    pub fn add_entry(
        &mut self,
        col: usize,
        value: u64,
        prime: u64,
    ) -> Result<(), SparseModRowError> {
        if prime == 0 {
            return Err(SparseModRowError::ZeroModulus);
        }
        self.add_entry_unchecked(col, value, prime);
        Ok(())
    }

    fn add_entry_unchecked(&mut self, col: usize, value: u64, prime: u64) {
        debug_assert_ne!(prime, 0);
        let value = value % prime;
        if value == 0 {
            return;
        }
        match self.entries.binary_search_by_key(&col, |&(col, _)| col) {
            Ok(index) => {
                let new_value = add_mod_u64(self.entries[index].1, value, prime);
                if new_value == 0 {
                    self.entries.remove(index);
                } else {
                    self.entries[index].1 = new_value;
                }
            }
            Err(index) => self.entries.insert(index, (col, value)),
        }
    }

    /// Add a signed value to the coefficient of `col`.
    pub fn add_signed_entry(
        &mut self,
        col: usize,
        value: i64,
        prime: u64,
    ) -> Result<(), SparseModRowError> {
        if prime == 0 {
            return Err(SparseModRowError::ZeroModulus);
        }
        self.add_entry_unchecked(col, signed_mod_u64(value, prime), prime);
        Ok(())
    }

    pub(crate) fn add_reduced_entry_sorted(&mut self, col: usize, value: u64, prime: u64) {
        debug_assert!(value < prime);
        if value == 0 {
            return;
        }
        let Some(&(last_col, _)) = self.entries.last() else {
            self.entries.push((col, value));
            return;
        };
        if col > last_col {
            self.entries.push((col, value));
        } else {
            self.add_entry_unchecked(col, value, prime);
        }
    }

    /// Return true if the row is the tautology `0 = 0`.
    pub fn is_tautology(&self) -> bool {
        self.entries.is_empty() && self.rhs == 0
    }

    /// Return true if the row is the contradiction `0 = c` with `c != 0`.
    pub fn is_contradiction(&self) -> bool {
        self.entries.is_empty() && self.rhs != 0
    }

    fn leading_entry(&self) -> Option<(usize, u64)> {
        self.entries.first().copied()
    }

    fn coefficient_at(&self, col: usize) -> Option<u64> {
        self.entries
            .binary_search_by_key(&col, |&(entry_col, _)| entry_col)
            .ok()
            .map(|index| self.entries[index].1)
    }

    fn normalize_pivot(&mut self, prime: u64) {
        let Some((_, pivot)) = self.leading_entry() else {
            return;
        };
        let pivot_inv = mod_inverse_prime(pivot, prime)
            .expect("nonzero element modulo prime should be invertible");
        for (_, value) in &mut self.entries {
            *value = mul_mod_u64(*value, pivot_inv, prime);
        }
        self.rhs = mul_mod_u64(self.rhs, pivot_inv, prime);
    }

    fn normalize_pivot_at(&mut self, pivot_col: usize, prime: u64) {
        let Some(pivot) = self.coefficient_at(pivot_col) else {
            return;
        };
        let pivot_inv = mod_inverse_prime(pivot, prime)
            .expect("nonzero element modulo prime should be invertible");
        for (_, value) in &mut self.entries {
            *value = mul_mod_u64(*value, pivot_inv, prime);
        }
        self.rhs = mul_mod_u64(self.rhs, pivot_inv, prime);
    }

    fn subtract_entry(&mut self, col: usize, value: u64, prime: u64) {
        let value = value % prime;
        if value == 0 {
            return;
        }
        match self.entries.binary_search_by_key(&col, |&(col, _)| col) {
            Ok(index) => {
                let new_value = sub_mod_u64(self.entries[index].1, value, prime);
                if new_value == 0 {
                    self.entries.remove(index);
                } else {
                    self.entries[index].1 = new_value;
                }
            }
            Err(index) => {
                let value = sub_mod_u64(0, value, prime);
                self.entries.insert(index, (col, value));
            }
        }
    }

    fn subtract_scaled_small_pivot(&mut self, pivot: &SparseModRow, factor: u64, prime: u64) {
        if self
            .entries
            .first()
            .is_some_and(|&(col, value)| col == pivot.entries[0].0 && value == factor)
        {
            self.entries.remove(0);
            for &(col, value) in pivot.entries.iter().skip(1) {
                self.subtract_entry(col, mul_mod_u64(factor, value, prime), prime);
            }
        } else {
            for &(col, value) in &pivot.entries {
                self.subtract_entry(col, mul_mod_u64(factor, value, prime), prime);
            }
        }
        self.rhs = sub_mod_u64(self.rhs, mul_mod_u64(factor, pivot.rhs, prime), prime);
    }

    fn subtract_scaled_leading_pivot(&mut self, pivot: &SparseModRow, factor: u64, prime: u64) {
        let mut merged = Vec::with_capacity(
            self.entries
                .len()
                .saturating_add(pivot.entries.len().saturating_sub(2)),
        );
        let mut lhs = 1;
        let mut rhs = 1;
        while lhs < self.entries.len() || rhs < pivot.entries.len() {
            let next_col = match (self.entries.get(lhs), pivot.entries.get(rhs)) {
                (Some(&(lhs_col, _)), Some(&(rhs_col, _))) => lhs_col.min(rhs_col),
                (Some(&(lhs_col, _)), None) => lhs_col,
                (None, Some(&(rhs_col, _))) => rhs_col,
                (None, None) => break,
            };

            let mut value = 0;
            if lhs < self.entries.len() && self.entries[lhs].0 == next_col {
                value = self.entries[lhs].1;
                lhs += 1;
            }
            if rhs < pivot.entries.len() && pivot.entries[rhs].0 == next_col {
                value = sub_mod_u64(
                    value,
                    mul_mod_u64(factor, pivot.entries[rhs].1, prime),
                    prime,
                );
                rhs += 1;
            }
            if value != 0 {
                merged.push((next_col, value));
            }
        }

        self.entries = merged;
        self.rhs = sub_mod_u64(self.rhs, mul_mod_u64(factor, pivot.rhs, prime), prime);
    }

    fn subtract_scaled(&mut self, pivot: &SparseModRow, factor: u64, prime: u64) {
        let factor = factor % prime;
        if factor == 0 {
            return;
        }
        if self
            .entries
            .first()
            .is_some_and(|&(col, value)| col == pivot.entries[0].0 && value == factor)
        {
            self.subtract_scaled_leading_pivot(pivot, factor, prime);
            return;
        }
        if pivot.entries.len() <= 4 {
            self.subtract_scaled_small_pivot(pivot, factor, prime);
            return;
        }

        let mut merged = Vec::with_capacity(self.entries.len() + pivot.entries.len());
        let mut lhs = 0;
        let mut rhs = 0;
        while lhs < self.entries.len() || rhs < pivot.entries.len() {
            let next_col = match (self.entries.get(lhs), pivot.entries.get(rhs)) {
                (Some(&(lhs_col, _)), Some(&(rhs_col, _))) => lhs_col.min(rhs_col),
                (Some(&(lhs_col, _)), None) => lhs_col,
                (None, Some(&(rhs_col, _))) => rhs_col,
                (None, None) => break,
            };

            let mut value = 0;
            if lhs < self.entries.len() && self.entries[lhs].0 == next_col {
                value = self.entries[lhs].1;
                lhs += 1;
            }
            if rhs < pivot.entries.len() && pivot.entries[rhs].0 == next_col {
                value = sub_mod_u64(
                    value,
                    mul_mod_u64(factor, pivot.entries[rhs].1, prime),
                    prime,
                );
                rhs += 1;
            }
            if value != 0 {
                merged.push((next_col, value));
            }
        }

        self.entries = merged;
        self.rhs = sub_mod_u64(self.rhs, mul_mod_u64(factor, pivot.rhs, prime), prime);
    }

    fn subtract_scaled_general(&mut self, pivot: &SparseModRow, factor: u64, prime: u64) {
        let factor = factor % prime;
        if factor == 0 {
            return;
        }

        let mut merged = Vec::with_capacity(self.entries.len() + pivot.entries.len());
        let mut lhs = 0;
        let mut rhs = 0;
        while lhs < self.entries.len() || rhs < pivot.entries.len() {
            let next_col = match (self.entries.get(lhs), pivot.entries.get(rhs)) {
                (Some(&(lhs_col, _)), Some(&(rhs_col, _))) => lhs_col.min(rhs_col),
                (Some(&(lhs_col, _)), None) => lhs_col,
                (None, Some(&(rhs_col, _))) => rhs_col,
                (None, None) => break,
            };

            let mut value = 0;
            if lhs < self.entries.len() && self.entries[lhs].0 == next_col {
                value = self.entries[lhs].1;
                lhs += 1;
            }
            if rhs < pivot.entries.len() && pivot.entries[rhs].0 == next_col {
                value = sub_mod_u64(
                    value,
                    mul_mod_u64(factor, pivot.entries[rhs].1, prime),
                    prime,
                );
                rhs += 1;
            }
            if value != 0 {
                merged.push((next_col, value));
            }
        }

        self.entries = merged;
        self.rhs = sub_mod_u64(self.rhs, mul_mod_u64(factor, pivot.rhs, prime), prime);
    }
}

/// Row ordering strategy for sparse modular elimination.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum SparseModRowOrder {
    /// Process rows in input order.  This streams rows without collecting them.
    #[default]
    Input,
    /// Process rows by increasing number of nonzero coefficients.
    ///
    /// This is a cheap fill-in heuristic.  The result still reports original
    /// zero-based input row indices in `pivot_rows` and `inconsistent_row`.
    IncreasingNonzeros,
    /// Process rows by an initial Markowitz-style fill-in score.
    ///
    /// The score is `(row_nnz - 1) * (leading_column_nnz - 1)`, where
    /// `leading_column_nnz` is computed from the unreduced input system.  This
    /// is cheaper than dynamic Markowitz pivoting but tends to prefer sparse
    /// rows whose leading pivot column is also sparse.
    IncreasingMarkowitzCost,
    /// Repeatedly choose a pivot entry with low current Markowitz fill score.
    ///
    /// This collects the matrix and updates all active rows after each pivot,
    /// so it has higher overhead than the streaming backends. It can preserve
    /// sparsity much better on systems where the leading-column order produces
    /// unnecessary fill-in.
    DynamicMarkowitz,
}

/// Options for sparse modular elimination.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SparseModEliminationOptions {
    /// Row ordering strategy.
    pub row_order: SparseModRowOrder,
    /// Whether to compute one solution when the system is consistent.
    pub compute_solution: bool,
}

impl Default for SparseModEliminationOptions {
    fn default() -> Self {
        Self {
            row_order: SparseModRowOrder::Input,
            compute_solution: true,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum SparseModSolutionMode {
    None,
    AnyConsistent,
    FullRank,
}

impl SparseModSolutionMode {
    fn from_compute_solution(compute_solution: bool) -> Self {
        if compute_solution {
            Self::AnyConsistent
        } else {
            Self::None
        }
    }
}

/// Lightweight statistics from sparse modular elimination.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct SparseModEliminationStats {
    /// Number of input rows seen by the solver.
    pub input_rows: usize,
    /// Number of nonzero coefficients in the input rows.
    pub input_nonzeros: usize,
    /// Number of tautological rows reduced to `0 = 0`.
    pub zero_rows: usize,
    /// Number of sparse row reductions against existing pivot rows.
    pub row_reductions: usize,
    /// Number of rows checked by evaluating the unique solution after full
    /// column rank was reached.
    pub full_rank_checks: usize,
    /// Sum of nonzero coefficients in stored pivot rows.
    pub pivot_nonzeros: usize,
    /// Maximum number of nonzero coefficients in any active row during
    /// elimination.
    pub max_active_nonzeros: usize,
}

/// Summary of sparse Gaussian elimination over a prime field.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SparseModEliminationResult {
    /// Whether the input system is consistent.
    pub consistent: bool,
    /// Number of pivot rows found before termination.
    pub rank: usize,
    /// Zero-based input row indices that introduced new pivots.
    pub pivot_rows: Vec<usize>,
    /// Pivot columns in the order they were inserted.
    pub pivot_columns: Vec<usize>,
    /// Zero-based input row index witnessing inconsistency, if any.
    pub inconsistent_row: Option<usize>,
    /// One solution with all free variables set to zero, if the system is
    /// consistent and solution computation was requested.
    pub solution: Option<Vec<u64>>,
    /// Sparse elimination statistics useful for tuning row-order and fill-in
    /// heuristics.
    pub stats: SparseModEliminationStats,
}

/// Errors for sparse modular linear-system elimination.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SparseModEliminationError {
    /// The modulus is not prime.  The elimination routine divides by nonzero
    /// pivots, so it must run over a field.
    ModulusNotPrime { modulus: u64 },
    /// A row contains a column outside `0..num_vars`.
    ColumnOutOfRange {
        row: usize,
        column: usize,
        num_vars: usize,
    },
}

impl std::fmt::Display for SparseModEliminationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ModulusNotPrime { modulus } => {
                write!(f, "modulus {modulus} is not prime")
            }
            Self::ColumnOutOfRange {
                row,
                column,
                num_vars,
            } => write!(
                f,
                "row {row} contains column {column}, but the system has {num_vars} variables"
            ),
        }
    }
}

impl std::error::Error for SparseModEliminationError {}

fn sparse_modular_solution_from_pivots(
    pivots: &[Option<SparseModRow>],
    num_vars: usize,
    prime: u64,
) -> Vec<u64> {
    let mut solution = vec![0; num_vars];
    for pivot_col in (0..num_vars).rev() {
        let Some(row) = &pivots[pivot_col] else {
            continue;
        };
        let mut value = row.rhs;
        for &(col, coeff) in row.entries.iter().skip(1) {
            value = sub_mod_u64(value, mul_mod_u64(coeff, solution[col], prime), prime);
        }
        solution[pivot_col] = value;
    }
    solution
}

fn sparse_modular_solution_from_ordered_pivots(
    pivot_rows: &[SparseModRow],
    pivot_columns: &[usize],
    num_vars: usize,
    prime: u64,
) -> Vec<u64> {
    let mut solution = vec![0; num_vars];
    for (row, &pivot_col) in pivot_rows.iter().zip(pivot_columns.iter()).rev() {
        let mut value = row.rhs;
        for &(col, coeff) in &row.entries {
            if col != pivot_col {
                value = sub_mod_u64(value, mul_mod_u64(coeff, solution[col], prime), prime);
            }
        }
        solution[pivot_col] = value;
    }
    solution
}

fn sparse_modular_row_matches_solution(row: &SparseModRow, solution: &[u64], prime: u64) -> bool {
    row.is_satisfied_by_validated_solution(solution, prime)
}

fn sparse_modular_column_counts(
    indexed_rows: &[(usize, SparseModRow)],
    num_vars: usize,
) -> Vec<usize> {
    let mut column_counts = vec![0; num_vars];
    for (_, row) in indexed_rows {
        for &(col, _) in &row.entries {
            if col < num_vars {
                column_counts[col] += 1;
            }
        }
    }
    column_counts
}

fn sparse_modular_initial_markowitz_cost(
    row: &SparseModRow,
    column_counts: &[usize],
) -> (usize, usize, usize) {
    let Some(&(col, _)) = row.entries.first() else {
        return (0, 0, 0);
    };
    let row_fill = row.entries.len().saturating_sub(1);
    let column_fill = column_counts
        .get(col)
        .copied()
        .unwrap_or(usize::MAX)
        .saturating_sub(1);
    (
        row_fill.saturating_mul(column_fill),
        row.entries.len(),
        column_fill,
    )
}

fn sparse_modular_dynamic_markowitz_pivot(
    active_rows: &[(usize, SparseModRow)],
    num_vars: usize,
) -> Option<(usize, usize)> {
    let column_counts = sparse_modular_column_counts(active_rows, num_vars);
    let mut best: Option<(usize, usize, usize, usize, usize, usize)> = None;

    for (row_pos, (row_index, row)) in active_rows.iter().enumerate() {
        let row_fill = row.entries.len().saturating_sub(1);
        for &(col, _) in &row.entries {
            let column_fill = column_counts[col].saturating_sub(1);
            let score = row_fill.saturating_mul(column_fill);
            let candidate = (
                score,
                row.entries.len(),
                column_counts[col],
                *row_index,
                col,
                row_pos,
            );
            if best.is_none_or(|best| candidate < best) {
                best = Some(candidate);
            }
        }
    }

    best.map(|(_, _, _, _, col, row_pos)| (row_pos, col))
}

fn sparse_modular_linear_system_consistency_indexed<I>(
    rows: I,
    num_vars: usize,
    prime: u64,
    solution_mode: SparseModSolutionMode,
) -> Result<SparseModEliminationResult, SparseModEliminationError>
where
    I: IntoIterator<Item = (usize, SparseModRow)>,
{
    let mut pivots = vec![None; num_vars];
    let mut pivot_rows = Vec::new();
    let mut pivot_columns = Vec::new();
    let mut stats = SparseModEliminationStats::default();
    let mut full_rank_solution: Option<Vec<u64>> = None;

    for (row_index, mut row) in rows {
        stats.input_rows += 1;
        stats.input_nonzeros += row.entries.len();
        stats.max_active_nonzeros = stats.max_active_nonzeros.max(row.entries.len());
        for &(column, _) in &row.entries {
            if column >= num_vars {
                return Err(SparseModEliminationError::ColumnOutOfRange {
                    row: row_index,
                    column,
                    num_vars,
                });
            }
        }

        if let Some(solution) = &full_rank_solution {
            stats.full_rank_checks += 1;
            if !sparse_modular_row_matches_solution(&row, solution, prime) {
                return Ok(SparseModEliminationResult {
                    consistent: false,
                    rank: pivot_columns.len(),
                    pivot_rows,
                    pivot_columns,
                    inconsistent_row: Some(row_index),
                    solution: None,
                    stats,
                });
            }
            continue;
        }

        loop {
            let Some((pivot_col, pivot_coeff)) = row.leading_entry() else {
                if row.rhs != 0 {
                    return Ok(SparseModEliminationResult {
                        consistent: false,
                        rank: pivot_columns.len(),
                        pivot_rows,
                        pivot_columns,
                        inconsistent_row: Some(row_index),
                        solution: None,
                        stats,
                    });
                }
                stats.zero_rows += 1;
                break;
            };

            let Some(pivot) = &pivots[pivot_col] else {
                row.normalize_pivot(prime);
                stats.pivot_nonzeros += row.entries.len();
                pivots[pivot_col] = Some(row);
                pivot_rows.push(row_index);
                pivot_columns.push(pivot_col);
                if pivot_columns.len() == num_vars {
                    full_rank_solution = Some(sparse_modular_solution_from_pivots(
                        &pivots, num_vars, prime,
                    ));
                }
                break;
            };
            row.subtract_scaled(pivot, pivot_coeff, prime);
            stats.row_reductions += 1;
            stats.max_active_nonzeros = stats.max_active_nonzeros.max(row.entries.len());
        }
    }
    let solution = match solution_mode {
        SparseModSolutionMode::None => None,
        SparseModSolutionMode::AnyConsistent => full_rank_solution.or_else(|| {
            Some(sparse_modular_solution_from_pivots(
                &pivots, num_vars, prime,
            ))
        }),
        SparseModSolutionMode::FullRank => full_rank_solution,
    };

    Ok(SparseModEliminationResult {
        consistent: true,
        rank: pivot_columns.len(),
        pivot_rows,
        pivot_columns,
        inconsistent_row: None,
        solution,
        stats,
    })
}

fn sparse_modular_linear_system_consistency_dynamic_markowitz(
    rows: Vec<(usize, SparseModRow)>,
    num_vars: usize,
    prime: u64,
    solution_mode: SparseModSolutionMode,
) -> Result<SparseModEliminationResult, SparseModEliminationError> {
    let mut active_rows = Vec::with_capacity(rows.len());
    let mut pivot_rows = Vec::new();
    let mut pivot_columns = Vec::new();
    let mut ordered_pivots = Vec::new();
    let mut stats = SparseModEliminationStats::default();

    for (row_index, row) in rows {
        stats.input_rows += 1;
        stats.input_nonzeros += row.entries.len();
        stats.max_active_nonzeros = stats.max_active_nonzeros.max(row.entries.len());
        for &(column, _) in &row.entries {
            if column >= num_vars {
                return Err(SparseModEliminationError::ColumnOutOfRange {
                    row: row_index,
                    column,
                    num_vars,
                });
            }
        }
        if row.is_contradiction() {
            return Ok(SparseModEliminationResult {
                consistent: false,
                rank: 0,
                pivot_rows,
                pivot_columns,
                inconsistent_row: Some(row_index),
                solution: None,
                stats,
            });
        }
        if row.is_tautology() {
            stats.zero_rows += 1;
        } else {
            active_rows.push((row_index, row));
        }
    }

    loop {
        let mut row_pos = 0;
        while row_pos < active_rows.len() {
            if active_rows[row_pos].1.is_contradiction() {
                return Ok(SparseModEliminationResult {
                    consistent: false,
                    rank: pivot_columns.len(),
                    pivot_rows,
                    pivot_columns,
                    inconsistent_row: Some(active_rows[row_pos].0),
                    solution: None,
                    stats,
                });
            }
            if active_rows[row_pos].1.is_tautology() {
                stats.zero_rows += 1;
                active_rows.swap_remove(row_pos);
            } else {
                row_pos += 1;
            }
        }

        let Some((pivot_pos, pivot_col)) =
            sparse_modular_dynamic_markowitz_pivot(&active_rows, num_vars)
        else {
            break;
        };
        let (pivot_row_index, mut pivot) = active_rows.swap_remove(pivot_pos);
        pivot.normalize_pivot_at(pivot_col, prime);
        stats.pivot_nonzeros += pivot.entries.len();
        pivot_rows.push(pivot_row_index);
        pivot_columns.push(pivot_col);
        ordered_pivots.push(pivot);

        if pivot_columns.len() == num_vars {
            let solution = sparse_modular_solution_from_ordered_pivots(
                &ordered_pivots,
                &pivot_columns,
                num_vars,
                prime,
            );
            stats.full_rank_checks += active_rows.len();
            for (row_index, row) in &active_rows {
                if !sparse_modular_row_matches_solution(row, &solution, prime) {
                    return Ok(SparseModEliminationResult {
                        consistent: false,
                        rank: pivot_columns.len(),
                        pivot_rows,
                        pivot_columns,
                        inconsistent_row: Some(*row_index),
                        solution: None,
                        stats,
                    });
                }
            }
            return Ok(SparseModEliminationResult {
                consistent: true,
                rank: pivot_columns.len(),
                pivot_rows,
                pivot_columns,
                inconsistent_row: None,
                solution: match solution_mode {
                    SparseModSolutionMode::None => None,
                    SparseModSolutionMode::AnyConsistent | SparseModSolutionMode::FullRank => {
                        Some(solution)
                    }
                },
                stats,
            });
        }

        let pivot = ordered_pivots
            .last()
            .expect("pivot was just pushed before active-row reduction");
        for (_, row) in &mut active_rows {
            let Some(factor) = row.coefficient_at(pivot_col) else {
                continue;
            };
            row.subtract_scaled_general(pivot, factor, prime);
            stats.row_reductions += 1;
            stats.max_active_nonzeros = stats.max_active_nonzeros.max(row.entries.len());
        }
    }

    let solution = match solution_mode {
        SparseModSolutionMode::None | SparseModSolutionMode::FullRank => None,
        SparseModSolutionMode::AnyConsistent => Some(sparse_modular_solution_from_ordered_pivots(
            &ordered_pivots,
            &pivot_columns,
            num_vars,
            prime,
        )),
    };

    Ok(SparseModEliminationResult {
        consistent: true,
        rank: pivot_columns.len(),
        pivot_rows,
        pivot_columns,
        inconsistent_row: None,
        solution,
        stats,
    })
}

/// Check consistency of a sparse linear system over the prime field `F_p`.
///
/// The rows are augmented rows `A_i x = b_i`.  The implementation performs
/// incremental sparse Gaussian elimination, normalizing each stored pivot row.
/// It is intended for sparse systems where one mostly needs rank/consistency
/// information or a single solution over `F_p` before a more expensive exact
/// solve.
pub fn sparse_modular_linear_system_consistency<I>(
    rows: I,
    num_vars: usize,
    prime: u64,
) -> Result<SparseModEliminationResult, SparseModEliminationError>
where
    I: IntoIterator<Item = SparseModRow>,
{
    sparse_modular_linear_system_consistency_with_options(
        rows,
        num_vars,
        prime,
        SparseModEliminationOptions::default(),
    )
}

/// Check consistency of a sparse linear system over `F_p` with explicit
/// elimination options.
pub fn sparse_modular_linear_system_consistency_with_options<I>(
    rows: I,
    num_vars: usize,
    prime: u64,
    options: SparseModEliminationOptions,
) -> Result<SparseModEliminationResult, SparseModEliminationError>
where
    I: IntoIterator<Item = SparseModRow>,
{
    if !is_prime_u64(prime) {
        return Err(SparseModEliminationError::ModulusNotPrime { modulus: prime });
    }

    sparse_modular_linear_system_consistency_prime_unchecked_with_options(
        rows, num_vars, prime, options,
    )
}

pub(crate) fn sparse_modular_linear_system_consistency_prime_unchecked_with_options<I>(
    rows: I,
    num_vars: usize,
    prime: u64,
    options: SparseModEliminationOptions,
) -> Result<SparseModEliminationResult, SparseModEliminationError>
where
    I: IntoIterator<Item = SparseModRow>,
{
    sparse_modular_linear_system_consistency_prime_unchecked_with_solution_mode(
        rows,
        num_vars,
        prime,
        options.row_order,
        SparseModSolutionMode::from_compute_solution(options.compute_solution),
    )
}

pub(crate) fn sparse_modular_linear_system_consistency_prime_unchecked_with_solution_mode<I>(
    rows: I,
    num_vars: usize,
    prime: u64,
    row_order: SparseModRowOrder,
    solution_mode: SparseModSolutionMode,
) -> Result<SparseModEliminationResult, SparseModEliminationError>
where
    I: IntoIterator<Item = SparseModRow>,
{
    match row_order {
        SparseModRowOrder::Input => sparse_modular_linear_system_consistency_indexed(
            rows.into_iter().enumerate(),
            num_vars,
            prime,
            solution_mode,
        ),
        SparseModRowOrder::IncreasingNonzeros => {
            let mut indexed_rows = rows.into_iter().enumerate().collect::<Vec<_>>();
            indexed_rows
                .sort_unstable_by_key(|(index, row)| (row.entries.len(), row.rhs == 0, *index));
            sparse_modular_linear_system_consistency_indexed(
                indexed_rows,
                num_vars,
                prime,
                solution_mode,
            )
        }
        SparseModRowOrder::IncreasingMarkowitzCost => {
            let mut indexed_rows = rows.into_iter().enumerate().collect::<Vec<_>>();
            let column_counts = sparse_modular_column_counts(&indexed_rows, num_vars);
            indexed_rows.sort_unstable_by_key(|(index, row)| {
                (
                    sparse_modular_initial_markowitz_cost(row, &column_counts),
                    row.rhs == 0,
                    *index,
                )
            });
            sparse_modular_linear_system_consistency_indexed(
                indexed_rows,
                num_vars,
                prime,
                solution_mode,
            )
        }
        SparseModRowOrder::DynamicMarkowitz => {
            let indexed_rows = rows.into_iter().enumerate().collect::<Vec<_>>();
            sparse_modular_linear_system_consistency_dynamic_markowitz(
                indexed_rows,
                num_vars,
                prime,
                solution_mode,
            )
        }
    }
}

/// Convenience boolean wrapper around
/// [`sparse_modular_linear_system_consistency`].
pub fn sparse_modular_linear_system_consistent<I>(
    rows: I,
    num_vars: usize,
    prime: u64,
) -> Result<bool, SparseModEliminationError>
where
    I: IntoIterator<Item = SparseModRow>,
{
    sparse_modular_linear_system_consistency_with_options(
        rows,
        num_vars,
        prime,
        SparseModEliminationOptions {
            compute_solution: false,
            ..SparseModEliminationOptions::default()
        },
    )
    .map(|result| result.consistent)
}

/// Return one solution of a sparse linear system over `F_p`, if the system is
/// consistent.
///
/// Free variables are set to zero.  Use
/// [`sparse_modular_linear_system_consistency`] when rank, pivot row indices,
/// pivot columns, or an inconsistency witness are also needed.
pub fn sparse_modular_linear_system_solution<I>(
    rows: I,
    num_vars: usize,
    prime: u64,
) -> Result<Option<Vec<u64>>, SparseModEliminationError>
where
    I: IntoIterator<Item = SparseModRow>,
{
    sparse_modular_linear_system_consistency(rows, num_vars, prime).map(|result| result.solution)
}

fn previous_prime_at_or_below(mut n: u64) -> Option<u64> {
    if n < 2 {
        return None;
    }
    if n == 2 {
        return Some(2);
    }
    if n.is_multiple_of(2) {
        n -= 1;
    }
    while n >= 3 {
        if is_prime_u64(n) {
            return Some(n);
        }
        n = n.saturating_sub(2);
    }
    None
}

fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for p in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n == p {
            return true;
        }
        if n.is_multiple_of(p) {
            return false;
        }
    }

    let mut d = n - 1;
    let mut s = 0;
    while d.is_multiple_of(2) {
        d /= 2;
        s += 1;
    }

    for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if a >= n {
            continue;
        }
        let mut x = pow_mod_u64(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        let mut probably_prime_for_base = false;
        for _ in 1..s {
            x = mul_mod_u64(x, x, n);
            if x == n - 1 {
                probably_prime_for_base = true;
                break;
            }
        }
        if !probably_prime_for_base {
            return false;
        }
    }
    true
}

/// Compute the leading principal determinants by fraction-free Bareiss
/// elimination over matrices with entries in `Z[t]`.
///
/// The entries are [`Polynomial<BigInt>`] in ascending degree order.  All
/// divisions are checked as exact polynomial divisions over `Z[t]`.  This is
/// useful when a Bézout or Sylvester matrix depends polynomially on a parameter
/// and one wants exact leading principal minors without passing through
/// rational functions.
///
/// Returns `None` if the matrix is not square, a nonfinal required pivot is
/// zero, or an exact polynomial division fails.
pub fn bareiss_leading_principal_minors_polynomial_bigint(
    mat: &[Vec<Polynomial<BigInt>>],
) -> Option<Vec<Polynomial<BigInt>>> {
    let n = mat.len();
    if mat.iter().any(|row| row.len() != n) {
        return None;
    }
    if n == 0 {
        return Some(Vec::new());
    }

    let mut a = mat.to_vec();
    let mut denom = Polynomial::<BigInt>::one();
    let mut minors = Vec::with_capacity(n);

    for k in 0..n {
        let pivot = a[k][k].clone();
        if pivot.is_zero() {
            if k == n - 1 {
                minors.push(pivot);
                break;
            }
            return None;
        }
        minors.push(pivot.clone());
        if k == n - 1 {
            break;
        }

        let mut next_entries = Vec::with_capacity((n - k - 1) * (n - k - 1));
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                let numerator = a[i][j].clone() * pivot.clone() - a[i][k].clone() * a[k][j].clone();
                let value = exact_div_polynomial_bigint(&numerator, &denom)?;
                next_entries.push((i, j, value));
            }
        }
        for (i, j, value) in next_entries {
            a[i][j] = value;
        }
        denom = pivot;
    }

    Some(minors)
}

/// Compute a determinant in `Z[t]` by no-pivot fraction-free Bareiss
/// elimination.
///
/// Use this when the fixed pivot order is part of the certificate.  The return
/// value is `None` under the same conditions as
/// [`bareiss_leading_principal_minors_polynomial_bigint`].
pub fn bareiss_determinant_polynomial_bigint(
    mat: &[Vec<Polynomial<BigInt>>],
) -> Option<Polynomial<BigInt>> {
    if mat.is_empty() {
        return Some(Polynomial::one());
    }
    bareiss_leading_principal_minors_polynomial_bigint(mat).and_then(|mut minors| minors.pop())
}

fn exact_div_polynomial_bigint(
    numerator: &Polynomial<BigInt>,
    denominator: &Polynomial<BigInt>,
) -> Option<Polynomial<BigInt>> {
    if denominator.is_zero() {
        return None;
    }
    if numerator.is_zero() {
        return Some(Polynomial::zero());
    }

    let den = denominator.coeffs();
    let den_degree = den.len() - 1;
    let den_lc = den.last().expect("nonzero denominator");
    let mut rem = numerator.coeffs().to_vec();

    if rem.len() < den.len() {
        return if rem.iter().all(BigInt::is_zero) {
            Some(Polynomial::zero())
        } else {
            None
        };
    }

    let mut quotient = vec![BigInt::zero(); rem.len() - den_degree];
    trim_bigint_coeffs(&mut rem);
    while rem.len() >= den.len() && !rem.is_empty() {
        let shift = rem.len() - den.len();
        let rem_lc = rem.last().expect("nonempty remainder").clone();
        if (&rem_lc % den_lc) != BigInt::zero() {
            return None;
        }
        let q = rem_lc / den_lc;
        quotient[shift] = q.clone();
        for (i, c) in den.iter().enumerate() {
            rem[shift + i] -= &q * c;
        }
        trim_bigint_coeffs(&mut rem);
    }

    if rem.is_empty() {
        Some(Polynomial::new(quotient))
    } else {
        None
    }
}

/// A particular solution of a consistent exact linear system, together with
/// the rank information needed to decide whether it is identifiable.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LinearSystemSolution {
    /// One solution, with every free variable set to zero.
    pub particular: Vec<Q>,
    /// Rank of the coefficient matrix.
    pub rank: usize,
    /// Dimension of the solution space of the associated homogeneous system.
    pub nullity: usize,
}

/// Solve `Ax = b` over ℚ and report the rank and nullity of `A`.
///
/// Returns `None` if the dimensions are invalid or the system is
/// inconsistent. For an underdetermined consistent system, free variables in
/// [`LinearSystemSolution::particular`] are set to zero.
pub fn solve_linear_system_with_rank(a: &[Vec<Q>], b: &[Q]) -> Option<LinearSystemSolution> {
    let num_rows = a.len();
    if b.len() != num_rows {
        return None;
    }
    if num_rows == 0 {
        return Some(LinearSystemSolution {
            particular: vec![],
            rank: 0,
            nullity: 0,
        });
    }
    let num_cols = a[0].len();
    if a.iter().any(|row| row.len() != num_cols) {
        return None;
    }

    // Build augmented matrix [A | b].
    let mut aug: Vec<Vec<Q>> = a
        .iter()
        .zip(b.iter())
        .map(|(row, bi)| {
            let mut r = row.clone();
            r.push(bi.clone());
            r
        })
        .collect();

    // Forward + upward elimination (reduced row-echelon form).
    let mut pivot_cols: Vec<(usize, usize)> = Vec::new();
    let mut pivot_row = 0;

    for col in 0..num_cols {
        let Some(pr) = (pivot_row..num_rows).find(|&r| !aug[r][col].is_zero()) else {
            continue;
        };

        aug.swap(pivot_row, pr);
        pivot_cols.push((pivot_row, col));

        let pivot_val = aug[pivot_row][col].clone();
        let pivot_snapshot: Vec<_> = aug[pivot_row][col..=num_cols].to_vec();

        // Eliminate all other rows in this column.
        for row in 0..num_rows {
            if row == pivot_row || aug[row][col].is_zero() {
                continue;
            }
            let factor = aug[row][col].clone() / pivot_val.clone();
            for (aug_j, pivot_j) in aug[row][col..=num_cols]
                .iter_mut()
                .zip(pivot_snapshot.iter())
            {
                let sub = pivot_j.clone() * &factor;
                *aug_j -= sub;
            }
        }

        pivot_row += 1;
    }

    // Consistency check: any row [0 ... 0 | nonzero] means no solution.
    for aug_row in &aug[pivot_row..num_rows] {
        if !aug_row[num_cols].is_zero() {
            return None;
        }
    }

    // Extract solution (free variables = 0).
    let mut x = vec![Q::zero(); num_cols];
    for &(pr, pc) in &pivot_cols {
        x[pc] = aug[pr][num_cols].clone() / aug[pr][pc].clone();
    }
    let rank = pivot_cols.len();
    Some(LinearSystemSolution {
        particular: x,
        rank,
        nullity: num_cols - rank,
    })
}

/// Solve Ax = b via Gaussian elimination over ℚ.
///
/// Returns `None` if the system is inconsistent. For underdetermined systems,
/// free variables are set to zero. Use [`solve_linear_system_with_rank`] when
/// the caller needs to distinguish a unique solution from an arbitrary
/// particular solution.
pub fn solve_linear_system(a: &[Vec<Q>], b: &[Q]) -> Option<Vec<Q>> {
    solve_linear_system_with_rank(a, b).map(|solution| solution.particular)
}

/// Solve a nonsingular square system `Ax = b` over `Q`.
///
/// This uses ordinary forward elimination followed by back-substitution.  It is
/// faster than [`solve_linear_system`] when the caller already knows that the
/// chosen rows form a full-rank square subsystem, because it does not compute a
/// full reduced row-echelon form.  Returns `None` if the input is not square or
/// if a zero pivot proves that the matrix is singular.
pub fn solve_full_rank_square_linear_system(a: &[Vec<Q>], b: &[Q]) -> Option<Vec<Q>> {
    let n = a.len();
    if b.len() != n || a.iter().any(|row| row.len() != n) {
        return None;
    }
    if n == 0 {
        return Some(vec![]);
    }

    if let Some(integer_aug) = square_rational_system_as_bigint_augmented(a, b) {
        if let Some(solution) = solve_full_rank_square_integer_augmented(integer_aug) {
            return Some(solution);
        }
    }

    let mut aug: Vec<Vec<Q>> = a
        .iter()
        .zip(b.iter())
        .map(|(row, bi)| {
            let mut r = row.clone();
            r.push(bi.clone());
            r
        })
        .collect();

    for col in 0..n {
        let pivot_row = (col..n).find(|&row| !aug[row][col].is_zero())?;
        aug.swap(col, pivot_row);

        let pivot_val = aug[col][col].clone();
        let (pivot_block, lower_rows) = aug.split_at_mut(col + 1);
        let pivot_row = &pivot_block[col];
        for row in lower_rows {
            if row[col].is_zero() {
                continue;
            }
            let factor = row[col].clone() / &pivot_val;
            row[col] = Q::zero();
            for (aug_j, pivot_j) in row[col + 1..=n]
                .iter_mut()
                .zip(pivot_row[col + 1..=n].iter())
            {
                let sub = pivot_j.clone() * &factor;
                *aug_j -= sub;
            }
        }
    }

    let mut x = vec![Q::zero(); n];
    for row in (0..n).rev() {
        let mut value = aug[row][n].clone();
        for (col, coeff) in aug[row].iter().enumerate().take(n).skip(row + 1) {
            if !coeff.is_zero() {
                value -= coeff.clone() * &x[col];
            }
        }
        x[row] = value / &aug[row][row];
    }

    Some(x)
}

fn square_rational_system_as_bigint_augmented(a: &[Vec<Q>], b: &[Q]) -> Option<Vec<Vec<BigInt>>> {
    a.iter()
        .zip(b.iter())
        .map(|row| {
            let (row, rhs) = row;
            let mut augmented = row
                .iter()
                .map(|entry| entry.denom().is_one().then(|| entry.numer().clone()))
                .collect::<Option<Vec<_>>>()?;
            augmented.push(rhs.denom().is_one().then(|| rhs.numer().clone())?);
            Some(augmented)
        })
        .collect()
}

fn solve_full_rank_square_integer_augmented(mut aug: Vec<Vec<BigInt>>) -> Option<Vec<Q>> {
    let n = aug.len();
    if aug.iter().any(|row| row.len() != n + 1) {
        return None;
    }
    if n == 0 {
        return Some(vec![]);
    }

    let mut denom = BigInt::one();

    for col in 0..n {
        let pivot_row = (col..n).find(|&row| !aug[row][col].is_zero())?;
        aug.swap(col, pivot_row);
        let pivot = aug[col][col].clone();
        if pivot.is_zero() {
            return None;
        }
        if col == n - 1 {
            break;
        }

        let pivot_tail = aug[col][col + 1..=n].to_vec();
        for row in (col + 1)..n {
            let row_pivot = aug[row][col].clone();
            aug[row][col] = BigInt::zero();
            for (offset, target_col) in ((col + 1)..=n).enumerate() {
                let numerator = if row_pivot.is_zero() {
                    &aug[row][target_col] * &pivot
                } else {
                    &aug[row][target_col] * &pivot - &row_pivot * &pivot_tail[offset]
                };
                if &numerator % &denom != BigInt::zero() {
                    return None;
                }
                aug[row][target_col] = numerator / &denom;
            }
        }
        denom = pivot;
    }

    let mut x = vec![Q::zero(); n];
    for row in (0..n).rev() {
        let diagonal = aug[row][row].clone();
        if diagonal.is_zero() {
            return None;
        }
        let mut value = Q::from_integer(aug[row][n].clone());
        for (col, coeff) in aug[row].iter().enumerate().take(n).skip(row + 1) {
            if !coeff.is_zero() {
                value -= Q::from_integer(coeff.clone()) * &x[col];
            }
        }
        x[row] = value / Q::from_integer(diagonal);
    }

    Some(x)
}

// ---------------------------------------------------------------------------
// Total positivity
// ---------------------------------------------------------------------------

/// Exact failure reported by brute-force total-positivity checking.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TotalPositivityError {
    RaggedMatrix {
        row: usize,
        expected_columns: usize,
        found_columns: usize,
    },
    NegativeMinor {
        size: usize,
        rows: Vec<usize>,
        columns: Vec<usize>,
        determinant: BigInt,
    },
    ZeroMinor {
        size: usize,
        rows: Vec<usize>,
        columns: Vec<usize>,
    },
}

impl std::fmt::Display for TotalPositivityError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::RaggedMatrix {
                row,
                expected_columns,
                found_columns,
            } => write!(
                f,
                "matrix row {row} has {found_columns} columns; expected {expected_columns}"
            ),
            Self::NegativeMinor {
                size,
                rows,
                columns,
                determinant,
            } => write!(
                f,
                "{size}x{size} minor rows {rows:?} cols {columns:?} has det = {determinant} < 0"
            ),
            Self::ZeroMinor {
                size,
                rows,
                columns,
            } => write!(
                f,
                "{size}x{size} minor rows {rows:?} cols {columns:?} has det = 0 (strict mode)"
            ),
        }
    }
}

impl std::error::Error for TotalPositivityError {}

/// Check if an integer matrix is totally positive (all minors strictly positive,
/// except those that are trivially zero due to zero rows/columns).
///
/// More precisely, checks that all minors are **non-negative** (totally non-negative, TNN).
/// Set `strict` to require all non-trivially-zero minors to be strictly positive.
///
/// `max_minor_size` limits the largest minor checked. For a full check, pass
/// `min(nrows, ncols)`. Returns `Ok(())` if all minors pass, or a structured
/// [`TotalPositivityError`] describing ragged input or the first violation.
///
/// # Example
/// ```
/// use polytool::linalg::check_total_positivity;
/// // Pascal's triangle is totally positive
/// let mat = vec![
///     vec![1, 0, 0],
///     vec![1, 1, 0],
///     vec![1, 2, 1],
///     vec![1, 3, 3],
///     vec![1, 4, 6],
/// ];
/// assert!(check_total_positivity(&mat, 3, false).is_ok());
/// ```
pub fn check_total_positivity(
    mat: &[Vec<i64>],
    max_minor_size: usize,
    strict: bool,
) -> Result<(), TotalPositivityError> {
    let nrows = mat.len();
    if nrows == 0 {
        return Ok(());
    }
    let ncols = mat[0].len();
    if let Some((row, found_columns)) = mat
        .iter()
        .enumerate()
        .find_map(|(row, values)| (values.len() != ncols).then_some((row, values.len())))
    {
        return Err(TotalPositivityError::RaggedMatrix {
            row,
            expected_columns: ncols,
            found_columns,
        });
    }
    let max_k = max_minor_size.min(nrows).min(ncols);

    for k in 1..=max_k {
        for rows in CombinationIterator::new(nrows, k) {
            for cols in CombinationIterator::new(ncols, k) {
                let sub = extract_submatrix_i64(mat, &rows, &cols);
                let det = determinant(&sub);
                if det < BigInt::zero() {
                    return Err(TotalPositivityError::NegativeMinor {
                        size: k,
                        rows,
                        columns: cols,
                        determinant: det,
                    });
                }
                if strict && det.is_zero() {
                    // Check if this is a non-trivial zero
                    let all_zero = rows.iter().all(|&r| cols.iter().all(|&c| mat[r][c] == 0));
                    if !all_zero {
                        return Err(TotalPositivityError::ZeroMinor {
                            size: k,
                            rows,
                            columns: cols,
                        });
                    }
                }
            }
        }
    }
    Ok(())
}

/// Check total non-negativity: all minors >= 0.
/// Convenience wrapper around [`check_total_positivity`].
pub fn is_totally_nonnegative(mat: &[Vec<i64>], max_minor_size: usize) -> bool {
    check_total_positivity(mat, max_minor_size, false).is_ok()
}

/// Check total non-negativity via Neville elimination (adjacent-row subtraction).
///
/// This is *much* faster than [`check_total_positivity`]: O(n²m) instead of
/// exponential in `min(n,m)`. It checks ALL minors, not just up to a given size.
///
/// A matrix is TNN iff Neville elimination completes with all multipliers ≥ 0
/// and no negative entries arise during the process.
///
/// Returns `Ok(())` if TNN, or `Err(msg)` describing the failure.
///
/// # Example
/// ```
/// use polytool::linalg::check_tnn_neville;
/// let pascal = vec![
///     vec![1, 0, 0],
///     vec![1, 1, 0],
///     vec![1, 2, 1],
///     vec![1, 3, 3],
///     vec![1, 4, 6],
/// ];
/// assert!(check_tnn_neville(&pascal).is_ok());
/// ```
fn check_tnn_neville_rational_rows(a: &mut [Vec<Q>]) -> Result<(), String> {
    let nrows = a.len();
    let ncols = a.first().map_or(0, Vec::len);
    let min_dim = nrows.min(ncols);
    for k in 0..min_dim {
        for i in (k + 1..nrows).rev() {
            if a[i][k].is_zero() {
                continue;
            }
            if a[i - 1][k].is_zero() {
                // A positive entry below a zero pivot gives a negative 2x2
                // minor if the pivot row has any later positive entry. A
                // zero suffix row can be moved down without changing a minor.
                if ((k + 1)..ncols).any(|j| !a[i - 1][j].is_zero()) {
                    return Err(format!(
                        "Neville elimination: positive entry [{},{}] below zero pivot \
                         with nonzero pivot-row suffix",
                        i, k
                    ));
                }
                a.swap(i - 1, i);
                continue;
            }
            let multiplier = a[i][k].clone() / a[i - 1][k].clone();
            debug_assert!(multiplier >= Q::zero());

            for j in k..ncols {
                let sub = multiplier.clone() * a[i - 1][j].clone();
                a[i][j] -= sub;
            }
            a[i][k] = Q::zero();

            for j in (k + 1)..ncols {
                if a[i][j] < Q::zero() {
                    return Err(format!(
                        "Neville elimination: entry [{},{}] became negative \
                         (multiplier at column {}, rows [{},{}])",
                        i,
                        j,
                        k,
                        i - 1,
                        i
                    ));
                }
            }
        }
    }
    Ok(())
}

pub fn check_tnn_neville(mat: &[Vec<i64>]) -> Result<(), String> {
    let nrows = mat.len();
    if nrows == 0 {
        return Ok(());
    }
    let ncols = mat[0].len();
    if mat.iter().any(|row| row.len() != ncols) {
        return Err("matrix rows have inconsistent lengths".to_string());
    }
    if ncols == 0 {
        return Ok(());
    }

    // Convert to rationals for exact arithmetic
    let mut a: Vec<Vec<Q>> = mat
        .iter()
        .map(|row| {
            row.iter()
                .map(|&v| Q::from_integer(BigInt::from(v)))
                .collect()
        })
        .collect();

    // Check all entries are non-negative
    for i in 0..nrows {
        for j in 0..ncols {
            if a[i][j] < Q::zero() {
                return Err(format!("Entry [{},{}] = {} < 0", i, j, mat[i][j]));
            }
        }
    }

    check_tnn_neville_rational_rows(&mut a)?;

    // Row elimination detects the initial minors involving the first column;
    // applying the same test to the transpose checks the complementary family
    // of initial minors. Both families are required in the Neville criterion.
    let mut transpose = (0..ncols)
        .map(|j| {
            (0..nrows)
                .map(|i| Q::from_integer(BigInt::from(mat[i][j])))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    check_tnn_neville_rational_rows(&mut transpose)
}

/// Fast TNN check via Neville elimination. Returns bool.
/// Equivalent to `check_tnn_neville(mat).is_ok()`.
pub fn is_tnn(mat: &[Vec<i64>]) -> bool {
    check_tnn_neville(mat).is_ok()
}

/// Check TNN via Neville elimination on a BigInt matrix.
/// Use this when entries may overflow i64.
pub fn check_tnn_neville_bigint(mat: &[Vec<BigInt>]) -> Result<(), String> {
    let nrows = mat.len();
    if nrows == 0 {
        return Ok(());
    }
    let ncols = mat[0].len();
    if mat.iter().any(|row| row.len() != ncols) {
        return Err("matrix rows have inconsistent lengths".to_string());
    }
    if ncols == 0 {
        return Ok(());
    }

    let mut a: Vec<Vec<Q>> = mat
        .iter()
        .map(|row| row.iter().map(|v| Q::from_integer(v.clone())).collect())
        .collect();

    for i in 0..nrows {
        for j in 0..ncols {
            if a[i][j] < Q::zero() {
                return Err(format!("Entry [{},{}] < 0", i, j));
            }
        }
    }

    check_tnn_neville_rational_rows(&mut a)?;
    let mut transpose = (0..ncols)
        .map(|j| {
            (0..nrows)
                .map(|i| Q::from_integer(mat[i][j].clone()))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    check_tnn_neville_rational_rows(&mut transpose)
}

fn extract_submatrix_i64(mat: &[Vec<i64>], rows: &[usize], cols: &[usize]) -> Vec<Vec<BigInt>> {
    rows.iter()
        .map(|&r| cols.iter().map(|&c| BigInt::from(mat[r][c])).collect())
        .collect()
}

struct CombinationIterator {
    n: usize,
    current: Option<Vec<usize>>,
}

impl CombinationIterator {
    fn new(n: usize, k: usize) -> Self {
        let current = (k <= n).then(|| (0..k).collect());
        Self { n, current }
    }
}

impl Iterator for CombinationIterator {
    type Item = Vec<usize>;

    fn next(&mut self) -> Option<Self::Item> {
        let result = self.current.clone()?;
        let k = result.len();
        if k == 0 {
            self.current = None;
            return Some(result);
        }

        let current = self.current.as_mut().expect("combination exists");
        let Some(index) = (0..k).rev().find(|&i| current[i] < self.n - k + i) else {
            self.current = None;
            return Some(result);
        };
        current[index] += 1;
        for i in index + 1..k {
            current[i] = current[i - 1] + 1;
        }
        Some(result)
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn bigint_to_q(mat: &[Vec<BigInt>]) -> Vec<Vec<Q>> {
    mat.iter()
        .map(|row| row.iter().map(|v| Q::from_integer(v.clone())).collect())
        .collect()
}

fn trim_bigint_coeffs(coeffs: &mut Vec<BigInt>) {
    while coeffs.last().is_some_and(BigInt::is_zero) {
        coeffs.pop();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bi(v: i64) -> BigInt {
        BigInt::from(v)
    }

    fn bi_mat(rows: &[&[i64]]) -> Vec<Vec<BigInt>> {
        rows.iter()
            .map(|row| row.iter().map(|&v| bi(v)).collect())
            .collect()
    }

    fn poly(coeffs: &[i64]) -> Polynomial<BigInt> {
        Polynomial::from_i64_coeffs(coeffs)
    }

    // -----------------------------------------------------------------------
    // Determinant
    // -----------------------------------------------------------------------

    #[test]
    fn test_determinant_empty() {
        assert_eq!(determinant(&[]), BigInt::one());
    }

    #[test]
    fn test_determinant_1x1() {
        assert_eq!(determinant(&bi_mat(&[&[7]])), bi(7));
    }

    #[test]
    fn test_determinant_2x2() {
        // [[1, 2], [3, 4]] -> det = 1*4 - 2*3 = -2
        let m = bi_mat(&[&[1, 2], &[3, 4]]);
        assert_eq!(determinant(&m), bi(-2));
    }

    #[test]
    fn test_determinant_3x3() {
        // [[1, 2, 3], [4, 5, 6], [7, 8, 10]] -> det = 1*(50-48) - 2*(40-42) + 3*(32-35) = 2+4-9 = -3
        let m = bi_mat(&[&[1, 2, 3], &[4, 5, 6], &[7, 8, 10]]);
        assert_eq!(determinant(&m), bi(-3));
    }

    #[test]
    fn test_determinant_row_pivot_after_first_step() {
        let m = bi_mat(&[&[1, 1, 0], &[1, 1, 1], &[0, 1, 1]]);
        assert_eq!(determinant(&m), bi(-1));
    }

    #[test]
    fn test_determinant_singular() {
        // [[1, 2], [2, 4]] -> det = 0
        let m = bi_mat(&[&[1, 2], &[2, 4]]);
        assert_eq!(determinant(&m), BigInt::zero());
    }

    #[test]
    #[should_panic(expected = "determinant requires a square matrix")]
    fn test_determinant_rejects_nonsquare_matrix() {
        let m = bi_mat(&[&[1, 2, 3], &[4, 5, 6]]);

        let _ = determinant(&m);
    }

    #[test]
    fn test_determinant_identity() {
        let m = bi_mat(&[&[1, 0, 0], &[0, 1, 0], &[0, 0, 1]]);
        assert_eq!(determinant(&m), bi(1));
    }

    #[test]
    fn test_bareiss_leading_principal_minors_bigint() {
        let m = bi_mat(&[&[2, 1], &[1, 2]]);
        assert_eq!(
            bareiss_leading_principal_minors_bigint(&m),
            Some(vec![bi(2), bi(3)])
        );
        assert_eq!(bareiss_determinant_bigint(&m), Some(bi(3)));
    }

    #[test]
    fn test_bareiss_bigint_no_pivot_zero_pivot() {
        let m = bi_mat(&[&[0, 1], &[1, 0]]);
        assert_eq!(determinant(&m), bi(-1));
        assert_eq!(bareiss_leading_principal_minors_bigint(&m), None);
    }

    #[test]
    fn test_bareiss_bigint_final_zero_pivot() {
        let m = bi_mat(&[&[1, 1], &[1, 1]]);
        assert_eq!(
            bareiss_leading_principal_minors_bigint(&m),
            Some(vec![bi(1), bi(0)])
        );
        assert_eq!(bareiss_determinant_bigint(&m), Some(bi(0)));
    }

    #[test]
    fn test_modular_leading_principal_minors_match_bareiss() {
        let m = bi_mat(&[&[4, 2, 1], &[2, 5, 3], &[1, 3, 6]]);
        let bareiss = bareiss_leading_principal_minors_bigint(&m).unwrap();
        let modular = modular_leading_principal_minors_bigint(&m).unwrap();
        assert_eq!(modular, bareiss);
    }

    #[test]
    fn test_modular_leading_principal_minors_handles_zero_bareiss_pivot() {
        let m = bi_mat(&[&[0, 1], &[1, 0]]);
        assert_eq!(bareiss_leading_principal_minors_bigint(&m), None);
        assert_eq!(
            modular_leading_principal_minors_bigint(&m),
            Some(vec![bi(0), bi(-1)])
        );
    }

    fn shifted_hilbert_like_matrix(size: usize) -> Vec<Vec<BigInt>> {
        let mut mat = vec![vec![BigInt::zero(); size]; size];
        for i in 0..size {
            for j in 0..size {
                mat[i][j] = BigInt::from((i + j + 1) as i64);
            }
            mat[i][i] += BigInt::from((size * size) as i64);
        }
        mat
    }

    fn large_positive_diagonal_matrix(size: usize) -> Vec<Vec<BigInt>> {
        let diagonal = BigInt::one() << 700usize;
        let mut mat = vec![vec![BigInt::zero(); size]; size];
        for (i, row) in mat.iter_mut().enumerate() {
            row[i] = diagonal.clone();
        }
        mat
    }

    #[test]
    fn test_bareiss_polynomial_bigint_2x2() {
        // [[1+z, z], [z, 1+z]] has leading minors 1+z and 1+2z.
        let m = vec![
            vec![poly(&[1, 1]), poly(&[0, 1])],
            vec![poly(&[0, 1]), poly(&[1, 1])],
        ];
        assert_eq!(
            bareiss_leading_principal_minors_polynomial_bigint(&m),
            Some(vec![poly(&[1, 1]), poly(&[1, 2])])
        );
        assert_eq!(
            bareiss_determinant_polynomial_bigint(&m),
            Some(poly(&[1, 2]))
        );
    }

    #[test]
    fn test_bareiss_polynomial_bigint_final_zero_pivot() {
        let m = vec![
            vec![poly(&[1, 1]), poly(&[1, 1])],
            vec![poly(&[1, 1]), poly(&[1, 1])],
        ];
        assert_eq!(
            bareiss_leading_principal_minors_polynomial_bigint(&m),
            Some(vec![poly(&[1, 1]), Polynomial::zero()])
        );
        assert_eq!(
            bareiss_determinant_polynomial_bigint(&m),
            Some(Polynomial::zero())
        );
    }

    #[test]
    fn test_bareiss_polynomial_bigint_nonconstant_exact_division() {
        // The final stage divides by the nonconstant first pivot 1+z.
        let m = vec![
            vec![poly(&[1, 1]), poly(&[0, 1]), Polynomial::zero()],
            vec![poly(&[0, 1]), poly(&[1, 1]), Polynomial::zero()],
            vec![Polynomial::zero(), Polynomial::zero(), poly(&[1, 1])],
        ];
        assert_eq!(
            bareiss_leading_principal_minors_polynomial_bigint(&m),
            Some(vec![poly(&[1, 1]), poly(&[1, 2]), poly(&[1, 3, 2])])
        );
        assert_eq!(
            bareiss_determinant_polynomial_bigint(&m),
            Some(poly(&[1, 3, 2]))
        );
    }

    // -----------------------------------------------------------------------
    // Positive definiteness
    // -----------------------------------------------------------------------

    #[test]
    fn test_pd_identity() {
        let m = bi_mat(&[&[1, 0, 0], &[0, 1, 0], &[0, 0, 1]]);
        assert!(is_positive_definite(&m));
    }

    #[test]
    fn test_pd_2x2() {
        // [[2, 1], [1, 2]] -> pivots: 2, 2 - 1/2 = 3/2 > 0
        let m = bi_mat(&[&[2, 1], &[1, 2]]);
        assert!(is_positive_definite(&m));
        assert!(is_positive_definite_bareiss(&m));
        assert_eq!(is_positive_definite_modular(&m), Some(true));
    }

    #[test]
    fn test_pd_rejects_nonsymmetric_matrix() {
        let m = bi_mat(&[&[1, 2], &[0, 1]]);

        assert!(!is_positive_definite(&m));
        assert!(!is_positive_definite_bareiss(&m));
        assert_eq!(is_positive_definite_modular(&m), Some(false));
    }

    #[test]
    fn test_pd_default_agrees_with_paths_across_threshold() {
        for size in [
            MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD - 1,
            MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD,
        ] {
            let m = shifted_hilbert_like_matrix(size);
            let bareiss = is_positive_definite_bareiss(&m);
            let modular = is_positive_definite_modular(&m).unwrap();
            assert_eq!(bareiss, modular);
            assert_eq!(is_positive_definite(&m), modular);
        }
    }

    #[test]
    fn test_pd_default_falls_back_when_modular_bounds_do_not_certify() {
        let m = large_positive_diagonal_matrix(MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD);

        assert!(modular_leading_principal_minors_bigint(&m).is_none());
        assert_eq!(is_positive_definite_modular(&m), None);
        assert!(is_positive_definite_bareiss(&m));
        assert!(is_positive_definite(&m));
    }

    #[test]
    fn test_pd_negative_definite() {
        // [[-1, 0], [0, -1]] -> first pivot = -1 < 0
        let m = bi_mat(&[&[-1, 0], &[0, -1]]);
        assert!(!is_positive_definite(&m));
        assert_eq!(is_positive_definite_modular(&m), Some(false));
    }

    #[test]
    fn test_pd_zero_matrix() {
        let m = bi_mat(&[&[0, 0], &[0, 0]]);
        assert!(!is_positive_definite(&m)); // PD requires strict positivity
    }

    #[test]
    fn test_pd_indefinite() {
        // [[1, 0], [0, -1]]
        let m = bi_mat(&[&[1, 0], &[0, -1]]);
        assert!(!is_positive_definite(&m));
        assert_eq!(is_positive_definite_modular(&m), Some(false));
    }

    // -----------------------------------------------------------------------
    // Positive semi-definiteness
    // -----------------------------------------------------------------------

    #[test]
    fn test_psd_identity() {
        let m = bi_mat(&[&[1, 0], &[0, 1]]);
        assert!(is_positive_semidefinite(&m));
    }

    #[test]
    fn test_psd_zero_matrix() {
        let m = bi_mat(&[&[0, 0], &[0, 0]]);
        assert!(is_positive_semidefinite(&m));
    }

    #[test]
    fn test_psd_rank_deficient() {
        // [[1, 1], [1, 1]] -> eigenvalues 0, 2 -> PSD
        let m = bi_mat(&[&[1, 1], &[1, 1]]);
        assert!(is_positive_semidefinite(&m));
    }

    #[test]
    fn test_psd_indefinite() {
        // [[1, 2], [2, 1]] -> eigenvalues 3, -1 -> not PSD
        let m = bi_mat(&[&[1, 2], &[2, 1]]);
        assert!(!is_positive_semidefinite(&m));
    }

    #[test]
    fn test_psd_rejects_nonsquare_matrix() {
        let m = bi_mat(&[&[1, 0, 0], &[0, 1, 0]]);

        assert!(!is_positive_semidefinite(&m));
    }

    #[test]
    fn test_psd_rejects_nonsymmetric_matrix() {
        let m = bi_mat(&[&[1, 0], &[2, 1]]);

        assert!(!is_positive_semidefinite(&m));
    }

    // -----------------------------------------------------------------------
    // Solve linear system
    // -----------------------------------------------------------------------

    #[test]
    fn linear_system_reports_rank_and_nullity() {
        let a = vec![
            vec![Q::one(), Q::one(), Q::zero()],
            vec![
                Q::from_integer(BigInt::from(2)),
                Q::from_integer(BigInt::from(2)),
                Q::zero(),
            ],
        ];
        let b = vec![Q::one(), Q::from_integer(BigInt::from(2))];
        let solution = solve_linear_system_with_rank(&a, &b).unwrap();

        assert_eq!(solution.rank, 1);
        assert_eq!(solution.nullity, 2);
        assert_eq!(solution.particular, vec![Q::one(), Q::zero(), Q::zero()]);
    }

    #[test]
    fn test_solve_identity() {
        let a = vec![vec![Q::one(), Q::zero()], vec![Q::zero(), Q::one()]];
        let b = vec![Q::from_integer(bi(3)), Q::from_integer(bi(5))];
        let x = solve_linear_system(&a, &b).unwrap();
        assert_eq!(x[0], Q::from_integer(bi(3)));
        assert_eq!(x[1], Q::from_integer(bi(5)));
    }

    #[test]
    fn test_solve_2x2() {
        // x + 2y = 5, 3x + 4y = 11 -> x = 1, y = 2
        let a = vec![
            vec![Q::from_integer(bi(1)), Q::from_integer(bi(2))],
            vec![Q::from_integer(bi(3)), Q::from_integer(bi(4))],
        ];
        let b = vec![Q::from_integer(bi(5)), Q::from_integer(bi(11))];
        let x = solve_linear_system(&a, &b).unwrap();
        assert_eq!(x[0], Q::from_integer(bi(1)));
        assert_eq!(x[1], Q::from_integer(bi(2)));
    }

    #[test]
    fn test_solve_full_rank_square_2x2() {
        // x + 2y = 5, 3x + 4y = 11 -> x = 1, y = 2
        let a = vec![
            vec![Q::from_integer(bi(1)), Q::from_integer(bi(2))],
            vec![Q::from_integer(bi(3)), Q::from_integer(bi(4))],
        ];
        let b = vec![Q::from_integer(bi(5)), Q::from_integer(bi(11))];
        let x = solve_full_rank_square_linear_system(&a, &b).unwrap();
        assert_eq!(x[0], Q::from_integer(bi(1)));
        assert_eq!(x[1], Q::from_integer(bi(2)));
    }

    #[test]
    fn test_solve_full_rank_square_integer_row_swap() {
        // Solution is (3, -1, 2); the first pivot requires a row swap.
        let a = vec![
            vec![
                Q::from_integer(bi(0)),
                Q::from_integer(bi(2)),
                Q::from_integer(bi(1)),
            ],
            vec![
                Q::from_integer(bi(1)),
                Q::from_integer(bi(1)),
                Q::from_integer(bi(0)),
            ],
            vec![
                Q::from_integer(bi(2)),
                Q::from_integer(bi(3)),
                Q::from_integer(bi(1)),
            ],
        ];
        let b = vec![
            Q::from_integer(bi(0)),
            Q::from_integer(bi(2)),
            Q::from_integer(bi(5)),
        ];
        let x = solve_full_rank_square_linear_system(&a, &b).unwrap();
        assert_eq!(
            x,
            vec![
                Q::from_integer(bi(3)),
                Q::from_integer(bi(-1)),
                Q::from_integer(bi(2)),
            ]
        );
    }

    #[test]
    fn test_solve_full_rank_square_rational_fallback() {
        // x - y = 0, x/2 + y = 2 -> x = y = 4/3.
        let half = Q::new(bi(1), bi(2));
        let a = vec![
            vec![half, Q::one()],
            vec![Q::one(), Q::from_integer(bi(-1))],
        ];
        let b = vec![Q::from_integer(bi(2)), Q::zero()];
        let x = solve_full_rank_square_linear_system(&a, &b).unwrap();
        assert_eq!(x, vec![Q::new(bi(4), bi(3)), Q::new(bi(4), bi(3))]);
    }

    #[test]
    fn test_solve_full_rank_square_rejects_singular() {
        let a = vec![vec![Q::one(), Q::one()], vec![Q::one(), Q::one()]];
        let b = vec![Q::one(), Q::one()];
        assert!(solve_full_rank_square_linear_system(&a, &b).is_none());
    }

    #[test]
    fn test_solve_inconsistent() {
        // x + y = 1, x + y = 2 -> no solution
        let a = vec![vec![Q::one(), Q::one()], vec![Q::one(), Q::one()]];
        let b = vec![Q::one(), Q::from_integer(bi(2))];
        assert!(solve_linear_system(&a, &b).is_none());
    }

    #[test]
    fn test_solve_rejects_mismatched_shapes() {
        assert!(solve_linear_system(&[vec![Q::one()]], &[]).is_none());
        assert!(solve_linear_system(&[], &[Q::one()]).is_none());
        assert!(solve_linear_system(
            &[vec![Q::one(), Q::zero()], vec![Q::one()]],
            &[Q::one(), Q::one()]
        )
        .is_none());
    }

    // -----------------------------------------------------------------------
    // Sparse modular linear systems
    // -----------------------------------------------------------------------

    fn dense_modular_linear_system_consistent(
        mut aug: Vec<Vec<u64>>,
        num_vars: usize,
        prime: u64,
    ) -> bool {
        let num_rows = aug.len();
        let mut pivot_row = 0;

        for col in 0..num_vars {
            let Some(pr) = (pivot_row..num_rows).find(|&row| aug[row][col] != 0) else {
                continue;
            };
            aug.swap(pivot_row, pr);
            let pivot_inv = mod_inverse_prime(aug[pivot_row][col], prime)
                .expect("nonzero element modulo prime is invertible");
            let pivot_snapshot = aug[pivot_row][col..=num_vars].to_vec();

            for aug_row in aug.iter_mut().take(num_rows).skip(pivot_row + 1) {
                if aug_row[col] == 0 {
                    continue;
                }
                let factor = mul_mod_u64(aug_row[col], pivot_inv, prime);
                for (entry, pivot_entry) in aug_row[col..=num_vars]
                    .iter_mut()
                    .zip(pivot_snapshot.iter())
                {
                    *entry = sub_mod_u64(*entry, mul_mod_u64(factor, *pivot_entry, prime), prime);
                }
            }

            pivot_row += 1;
            if pivot_row == num_rows {
                return true;
            }
        }

        aug[pivot_row..]
            .iter()
            .all(|row| row[..num_vars].iter().any(|&entry| entry != 0) || row[num_vars] == 0)
    }

    fn sparse_rows_from_dense_aug(
        aug: &[Vec<u64>],
        num_vars: usize,
        prime: u64,
    ) -> Vec<SparseModRow> {
        aug.iter()
            .map(|dense_row| {
                let mut row = SparseModRow::new(dense_row[num_vars], prime).unwrap();
                for (col, &value) in dense_row.iter().take(num_vars).enumerate() {
                    row.add_entry(col, value, prime).unwrap();
                }
                row
            })
            .filter(|row| !row.is_tautology())
            .collect()
    }

    fn solution_satisfies_dense_aug(
        aug: &[Vec<u64>],
        num_vars: usize,
        prime: u64,
        solution: &[u64],
    ) -> bool {
        aug.iter().all(|row| {
            let lhs = row
                .iter()
                .take(num_vars)
                .zip(solution.iter())
                .fold(0, |acc, (&coeff, &value)| {
                    add_mod_u64(acc, mul_mod_u64(coeff, value, prime), prime)
                });
            lhs == row[num_vars] % prime
        })
    }

    #[test]
    fn test_sparse_mod_row_normalizes_and_cancels_entries() {
        let prime = 101;
        let mut row = SparseModRow::from_signed_rhs(-3, prime).unwrap();
        row.add_signed_entry(2, -2, prime).unwrap();
        row.add_entry(2, 103, prime).unwrap();
        row.add_signed_entry(1, -1, prime).unwrap();

        assert_eq!(row.rhs(), 98);
        assert_eq!(row.entries(), &[(1, 100)]);

        let row = SparseModRow::from_entries([(2, 5), (1, 3), (2, 96)], 204, prime).unwrap();
        assert_eq!(row.rhs(), 2);
        assert_eq!(row.entries(), &[(1, 3)]);
        assert!(!row.is_tautology());
        assert!(!row.is_contradiction());

        let row = SparseModRow::new(7, prime).unwrap();
        assert!(row.is_contradiction());
    }

    #[test]
    fn test_sparse_mod_sorted_add_fast_path_and_fallback() {
        let prime = 101;
        let mut row = SparseModRow::new(0, prime).unwrap();

        row.add_reduced_entry_sorted(1, 4, prime);
        row.add_reduced_entry_sorted(3, 8, prime);
        row.add_reduced_entry_sorted(2, 5, prime);
        row.add_reduced_entry_sorted(3, 93, prime);

        assert_eq!(row.entries(), &[(1, 4), (2, 5)]);
    }

    #[test]
    fn test_sparse_mod_small_pivot_reduction_cancels_in_place() {
        let prime = 101;
        let pivot = SparseModRow::from_entries([(0, 1), (2, 4)], 7, prime).unwrap();
        let mut row = SparseModRow::from_entries([(0, 5), (1, 3), (2, 20)], 9, prime).unwrap();

        row.subtract_scaled(&pivot, 5, prime);

        assert_eq!(row.rhs(), 75);
        assert_eq!(row.entries(), &[(1, 3)]);
    }

    #[test]
    fn test_sparse_modular_solver_matches_dense_consistency() {
        let prime = 101;
        let num_vars = 3;
        let aug = vec![
            vec![1, 2, 0, 5],
            vec![0, 3, 4, 7],
            vec![2, 7, 4, 17],
            vec![0, 0, 0, 0],
        ];
        let sparse_rows = sparse_rows_from_dense_aug(&aug, num_vars, prime);
        let result =
            sparse_modular_linear_system_consistency(sparse_rows, num_vars, prime).unwrap();

        assert_eq!(
            result.consistent,
            dense_modular_linear_system_consistent(aug.clone(), num_vars, prime)
        );
        assert_eq!(result.rank, 2);
        assert_eq!(result.pivot_rows, vec![0, 1]);
        assert_eq!(result.pivot_columns, vec![0, 1]);
        assert_eq!(result.inconsistent_row, None);
        let solution = result.solution.as_ref().expect("consistent system");
        assert_eq!(solution.len(), num_vars);
        assert!(solution_satisfies_dense_aug(
            &aug, num_vars, prime, solution
        ));

        let direct_solution = sparse_modular_linear_system_solution(
            sparse_rows_from_dense_aug(&aug, num_vars, prime),
            num_vars,
            prime,
        )
        .unwrap()
        .expect("consistent system");
        assert!(solution_satisfies_dense_aug(
            &aug,
            num_vars,
            prime,
            &direct_solution
        ));
    }

    #[test]
    fn test_sparse_modular_solver_can_process_sparsest_rows_first() {
        let prime = 101;
        let num_vars = 3;
        let aug = vec![vec![1, 1, 1, 1], vec![0, 1, 0, 2], vec![0, 0, 1, 3]];
        let result = sparse_modular_linear_system_consistency_with_options(
            sparse_rows_from_dense_aug(&aug, num_vars, prime),
            num_vars,
            prime,
            SparseModEliminationOptions {
                row_order: SparseModRowOrder::IncreasingNonzeros,
                compute_solution: true,
            },
        )
        .unwrap();

        assert!(result.consistent);
        assert_eq!(result.rank, 3);
        assert_eq!(result.pivot_rows, vec![1, 2, 0]);
        assert_eq!(result.pivot_columns, vec![1, 2, 0]);
        assert_eq!(result.stats.input_rows, 3);
        assert_eq!(result.stats.input_nonzeros, 5);
        assert_eq!(result.stats.row_reductions, 0);
        assert!(solution_satisfies_dense_aug(
            &aug,
            num_vars,
            prime,
            result.solution.as_ref().expect("consistent system")
        ));

        let result = sparse_modular_linear_system_consistency_with_options(
            sparse_rows_from_dense_aug(&aug, num_vars, prime),
            num_vars,
            prime,
            SparseModEliminationOptions {
                row_order: SparseModRowOrder::IncreasingNonzeros,
                compute_solution: false,
            },
        )
        .unwrap();
        assert!(result.consistent);
        assert_eq!(result.solution, None);
    }

    #[test]
    fn test_sparse_modular_solver_checks_tail_rows_after_full_rank() {
        let prime = 101;
        let num_vars = 2;
        let aug = vec![vec![1, 0, 2], vec![0, 1, 3], vec![1, 1, 5], vec![2, 3, 13]];
        let result = sparse_modular_linear_system_consistency_with_options(
            sparse_rows_from_dense_aug(&aug, num_vars, prime),
            num_vars,
            prime,
            SparseModEliminationOptions {
                row_order: SparseModRowOrder::Input,
                compute_solution: false,
            },
        )
        .unwrap();

        assert!(result.consistent);
        assert_eq!(result.rank, 2);
        assert_eq!(result.pivot_rows, vec![0, 1]);
        assert_eq!(result.stats.row_reductions, 0);
        assert_eq!(result.stats.full_rank_checks, 2);
        assert_eq!(result.solution, None);

        let inconsistent_aug = vec![vec![1, 0, 2], vec![0, 1, 3], vec![1, 1, 6]];
        let result = sparse_modular_linear_system_consistency_with_options(
            sparse_rows_from_dense_aug(&inconsistent_aug, num_vars, prime),
            num_vars,
            prime,
            SparseModEliminationOptions {
                row_order: SparseModRowOrder::Input,
                compute_solution: false,
            },
        )
        .unwrap();

        assert!(!result.consistent);
        assert_eq!(result.rank, 2);
        assert_eq!(result.inconsistent_row, Some(2));
        assert_eq!(result.stats.full_rank_checks, 1);
    }

    #[test]
    fn test_sparse_modular_solver_full_rank_solution_mode() {
        let prime = 101;
        let num_vars = 3;
        let rank_deficient_aug = vec![vec![1, 0, 0, 2], vec![0, 1, 0, 3]];
        let full_rank_only =
            sparse_modular_linear_system_consistency_prime_unchecked_with_solution_mode(
                sparse_rows_from_dense_aug(&rank_deficient_aug, num_vars, prime),
                num_vars,
                prime,
                SparseModRowOrder::IncreasingNonzeros,
                SparseModSolutionMode::FullRank,
            )
            .unwrap();
        assert!(full_rank_only.consistent);
        assert_eq!(full_rank_only.rank, 2);
        assert_eq!(full_rank_only.solution, None);

        let any_consistent =
            sparse_modular_linear_system_consistency_prime_unchecked_with_solution_mode(
                sparse_rows_from_dense_aug(&rank_deficient_aug, num_vars, prime),
                num_vars,
                prime,
                SparseModRowOrder::IncreasingNonzeros,
                SparseModSolutionMode::AnyConsistent,
            )
            .unwrap();
        assert_eq!(any_consistent.rank, 2);
        assert!(solution_satisfies_dense_aug(
            &rank_deficient_aug,
            num_vars,
            prime,
            any_consistent
                .solution
                .as_ref()
                .expect("rank-deficient solution")
        ));

        let full_rank_aug = vec![vec![1, 0, 0, 2], vec![0, 1, 0, 3], vec![0, 0, 1, 4]];
        let full_rank =
            sparse_modular_linear_system_consistency_prime_unchecked_with_solution_mode(
                sparse_rows_from_dense_aug(&full_rank_aug, num_vars, prime),
                num_vars,
                prime,
                SparseModRowOrder::IncreasingNonzeros,
                SparseModSolutionMode::FullRank,
            )
            .unwrap();
        assert_eq!(full_rank.rank, 3);
        assert!(solution_satisfies_dense_aug(
            &full_rank_aug,
            num_vars,
            prime,
            full_rank.solution.as_ref().expect("full-rank solution")
        ));
    }

    #[test]
    fn test_sparse_modular_solver_can_use_initial_markowitz_cost() {
        let prime = 101;
        let num_vars = 5;
        let rows = vec![
            SparseModRow::from_entries([(0, 1), (3, 1)], 1, prime).unwrap(),
            SparseModRow::from_entries([(0, 1), (4, 1)], 2, prime).unwrap(),
            SparseModRow::from_entries([(1, 1), (3, 1)], 3, prime).unwrap(),
            SparseModRow::from_entries([(1, 1), (4, 1)], 4, prime).unwrap(),
            SparseModRow::from_entries([(2, 1), (3, 1)], 5, prime).unwrap(),
        ];
        let result = sparse_modular_linear_system_consistency_with_options(
            rows,
            num_vars,
            prime,
            SparseModEliminationOptions {
                row_order: SparseModRowOrder::IncreasingMarkowitzCost,
                compute_solution: false,
            },
        )
        .unwrap();

        assert!(result.consistent);
        assert_eq!(result.pivot_rows[0], 4);
        assert_eq!(result.pivot_columns[0], 2);
        assert_eq!(result.solution, None);
    }

    #[test]
    fn test_sparse_modular_solver_dynamic_markowitz_matches_dense() {
        let prime = 101;
        let num_vars = 4;
        let aug = vec![
            vec![1, 1, 1, 1, 10],
            vec![0, 0, 1, 1, 7],
            vec![0, 1, 0, 1, 8],
            vec![1, 0, 0, 0, 2],
            vec![2, 1, 1, 2, 19],
            vec![0, 0, 0, 0, 0],
        ];
        let result = sparse_modular_linear_system_consistency_with_options(
            sparse_rows_from_dense_aug(&aug, num_vars, prime),
            num_vars,
            prime,
            SparseModEliminationOptions {
                row_order: SparseModRowOrder::DynamicMarkowitz,
                compute_solution: true,
            },
        )
        .unwrap();

        assert_eq!(
            result.consistent,
            dense_modular_linear_system_consistent(aug.clone(), num_vars, prime)
        );
        assert_eq!(result.rank, 4);
        assert_eq!(result.solution.as_ref().map(Vec::len), Some(num_vars));
        assert!(solution_satisfies_dense_aug(
            &aug,
            num_vars,
            prime,
            result.solution.as_ref().expect("consistent system")
        ));
        assert_eq!(result.stats.zero_rows, 0);
        assert!(result.stats.row_reductions > 0);
    }

    #[test]
    fn test_sparse_modular_solver_dynamic_markowitz_detects_inconsistency() {
        let prime = 101;
        let num_vars = 3;
        let aug = vec![vec![1, 0, 0, 2], vec![0, 1, 0, 3], vec![1, 1, 0, 6]];
        let result = sparse_modular_linear_system_consistency_with_options(
            sparse_rows_from_dense_aug(&aug, num_vars, prime),
            num_vars,
            prime,
            SparseModEliminationOptions {
                row_order: SparseModRowOrder::DynamicMarkowitz,
                compute_solution: true,
            },
        )
        .unwrap();

        assert!(!result.consistent);
        assert_eq!(result.rank, 2);
        assert_eq!(result.inconsistent_row, Some(2));
        assert_eq!(result.solution, None);
    }

    #[test]
    fn test_sparse_modular_solver_detects_inconsistency() {
        let prime = 101;
        let num_vars = 2;
        let aug = vec![vec![1, 1, 1], vec![2, 2, 3]];
        let sparse_rows = sparse_rows_from_dense_aug(&aug, num_vars, prime);
        let result =
            sparse_modular_linear_system_consistency(sparse_rows, num_vars, prime).unwrap();

        assert!(!result.consistent);
        assert_eq!(result.rank, 1);
        assert_eq!(result.pivot_rows, vec![0]);
        assert_eq!(result.pivot_columns, vec![0]);
        assert_eq!(result.inconsistent_row, Some(1));
        assert_eq!(result.solution, None);
    }

    #[test]
    fn test_sparse_modular_solver_rejects_bad_inputs() {
        let row = {
            let mut row = SparseModRow::new(0, 101).unwrap();
            row.add_entry(3, 1, 101).unwrap();
            row
        };
        assert_eq!(
            sparse_modular_linear_system_consistency(vec![row], 3, 101),
            Err(SparseModEliminationError::ColumnOutOfRange {
                row: 0,
                column: 3,
                num_vars: 3
            })
        );
        assert_eq!(
            sparse_modular_linear_system_consistency(Vec::new(), 3, 100),
            Err(SparseModEliminationError::ModulusNotPrime { modulus: 100 })
        );
        assert_eq!(SparseModRow::new(1, 0), Err(SparseModRowError::ZeroModulus));
        assert_eq!(
            SparseModRow::from_signed_rhs(-1, 0),
            Err(SparseModRowError::ZeroModulus)
        );
        let mut row = SparseModRow::new(0, 101).unwrap();
        assert_eq!(row.add_entry(0, 1, 0), Err(SparseModRowError::ZeroModulus));
    }

    // -----------------------------------------------------------------------
    // Total non-negativity
    // -----------------------------------------------------------------------

    #[test]
    fn test_tnn_pascal() {
        let pascal = vec![
            vec![1, 0, 0],
            vec![1, 1, 0],
            vec![1, 2, 1],
            vec![1, 3, 3],
            vec![1, 4, 6],
        ];
        assert!(is_tnn(&pascal));
        assert!(check_tnn_neville(&pascal).is_ok());
        assert!(check_total_positivity(&pascal, 3, false).is_ok());
    }

    #[test]
    fn test_tnn_identity() {
        let id = vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]];
        assert!(is_tnn(&id));
    }

    #[test]
    fn test_tnn_negative_entry() {
        let m = vec![vec![1, -1], vec![0, 1]];
        assert!(!is_tnn(&m));
    }

    #[test]
    fn test_tnn_rejects_row_swap_false_positive() {
        let permutation = vec![vec![0, 1], vec![1, 0]];
        assert!(check_tnn_neville(&permutation).is_err());

        let permutation_bigint = vec![vec![bi(0), bi(1)], vec![bi(1), bi(0)]];
        assert!(check_tnn_neville_bigint(&permutation_bigint).is_err());

        // Moving an all-zero suffix row remains safe and avoids rejecting a
        // genuinely TNN rank-deficient matrix.
        assert!(check_tnn_neville(&[vec![0, 0], vec![1, 1]]).is_ok());
    }

    #[test]
    fn test_tnn_rejects_ragged_matrices() {
        assert_eq!(
            check_total_positivity(&[vec![1, 0], vec![1]], 2, false),
            Err(TotalPositivityError::RaggedMatrix {
                row: 1,
                expected_columns: 2,
                found_columns: 1,
            })
        );
        assert!(check_tnn_neville(&[vec![1, 0], vec![1]]).is_err());
        assert!(check_tnn_neville_bigint(&[vec![bi(1), bi(0)], vec![bi(1)]]).is_err());
    }

    #[test]
    fn lazy_combinations_match_eager_reference() {
        fn eager(n: usize, k: usize) -> Vec<Vec<usize>> {
            fn generate(
                start: usize,
                remaining: usize,
                current: &mut Vec<usize>,
                output: &mut Vec<Vec<usize>>,
                n: usize,
            ) {
                if remaining == 0 {
                    output.push(current.clone());
                    return;
                }
                for value in start..=n - remaining {
                    current.push(value);
                    generate(value + 1, remaining - 1, current, output, n);
                    current.pop();
                }
            }

            if k > n {
                return Vec::new();
            }
            let mut output = Vec::new();
            generate(0, k, &mut Vec::new(), &mut output, n);
            output
        }

        for n in 0..=9 {
            for k in 0..=n + 1 {
                assert_eq!(
                    CombinationIterator::new(n, k).collect::<Vec<_>>(),
                    eager(n, k)
                );
            }
        }
    }

    #[test]
    fn test_tnn_not_tnn_positive_entries() {
        // All entries ≥ 0 but 2x2 minor det = 1*1 - 2*2 = -3 < 0
        let m = vec![vec![1, 2], vec![2, 1]];
        assert!(!is_tnn(&m));
    }

    #[test]
    fn test_tnn_path_matrix() {
        // Path matrix for P3: entry (i,j) = 1 if i <= j, 0 otherwise
        // These are always TNN
        let m = vec![vec![1, 1, 1], vec![0, 1, 1], vec![0, 0, 1]];
        assert!(is_tnn(&m));
    }

    #[test]
    fn test_tnn_neville_vs_brute() {
        // Verify Neville and brute-force agree on a 3x3 matrix
        let m = vec![vec![1, 1, 0], vec![0, 1, 1], vec![0, 0, 1]];
        let neville = check_tnn_neville(&m).is_ok();
        let brute = check_total_positivity(&m, 3, false).is_ok();
        assert_eq!(neville, brute);
    }

    #[test]
    fn test_tnn_neville_matches_brute_force_on_binary_3x3_matrices() {
        for bits in 0_u16..(1 << 9) {
            let mut matrix = vec![vec![0; 3]; 3];
            for entry in 0..9 {
                matrix[entry / 3][entry % 3] = i64::from((bits >> entry) & 1);
            }
            assert_eq!(
                check_tnn_neville(&matrix).is_ok(),
                check_total_positivity(&matrix, 3, false).is_ok(),
                "Neville and brute-force checks disagree for {matrix:?}"
            );
        }
    }

    #[test]
    fn test_tnn_bigint() {
        let m = vec![vec![bi(1), bi(1)], vec![bi(0), bi(1)]];
        assert!(check_tnn_neville_bigint(&m).is_ok());
    }
}
