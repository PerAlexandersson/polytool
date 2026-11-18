//! Exact symmetric inertia and replayable, symmetrically pivoted block LDLᵀ.
//!
//! [`symmetric_ldlt`] uses rational congruences with 1×1 and 2×2 pivots.
//! Zero diagonal entries never require a tolerance or an eigenvalue approximation.

use num_rational::BigRational as Q;
use num_traits::{One, Signed, Zero};
use std::fmt;

/// Numbers of positive, negative and zero eigenvalues, with multiplicity.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Inertia {
    /// Number of positive eigenvalues.
    pub positive: usize,
    /// Number of negative eigenvalues.
    pub negative: usize,
    /// Number of zero eigenvalues.
    pub zero: usize,
}

/// A diagonal block of an exact LDLᵀ decomposition.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DiagonalBlock {
    /// A 1×1 block, possibly zero.
    Single(Q),
    /// The block `[[0, b], [b, 0]]`, where `b` is nonzero.
    ///
    /// Such a block is used only when the remaining diagonal is zero.
    Pair(Q),
}

/// Exact certificate for `A[permutation, permutation] = L D Lᵀ`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SymmetricLdlt {
    /// Position `i` contains the original zero-based row/column index.
    pub permutation: Vec<usize>,
    /// Dense unit lower triangular matrix `L`, in permuted coordinates.
    pub lower: Vec<Vec<Q>>,
    /// Successive diagonal blocks of `D`.
    pub blocks: Vec<DiagonalBlock>,
}

/// Invalid input to symmetric rational decomposition.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum InertiaError {
    /// The matrix has a row whose length differs from its number of rows.
    NotSquare,
    /// An off-diagonal entry differs from its transpose.
    NotSymmetric,
}

impl fmt::Display for InertiaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::NotSquare => "matrix must be square",
            Self::NotSymmetric => "matrix must be symmetric",
        })
    }
}

impl std::error::Error for InertiaError {}

fn validate(matrix: &[Vec<Q>]) -> Result<(), InertiaError> {
    let n = matrix.len();
    if matrix.iter().any(|row| row.len() != n) {
        return Err(InertiaError::NotSquare);
    }
    for (i, row) in matrix.iter().enumerate() {
        for (j, entry) in row.iter().enumerate().take(i) {
            if entry != &matrix[j][i] {
                return Err(InertiaError::NotSymmetric);
            }
        }
    }
    Ok(())
}

impl SymmetricLdlt {
    /// Read the exact inertia from the diagonal blocks.
    ///
    /// For a caller-constructed certificate, call [`Self::verify`] first.
    pub fn inertia(&self) -> Inertia {
        let mut result = Inertia::default();
        for block in &self.blocks {
            match block {
                DiagonalBlock::Single(value) if value.is_positive() => result.positive += 1,
                DiagonalBlock::Single(value) if value.is_negative() => result.negative += 1,
                DiagonalBlock::Single(_) => result.zero += 1,
                DiagonalBlock::Pair(value) if value.is_zero() => result.zero += 2,
                DiagonalBlock::Pair(_) => {
                    result.positive += 1;
                    result.negative += 1;
                }
            }
        }
        result
    }

    /// Check shape, permutation, unit triangularity and `P A Pᵀ = L D Lᵀ`.
    ///
    /// Malformed or altered certificates return `false` without indexing panics.
    pub fn verify(&self, matrix: &[Vec<Q>]) -> bool {
        let n = matrix.len();
        if validate(matrix).is_err()
            || self.permutation.len() != n
            || self.lower.len() != n
            || self.lower.iter().any(|row| row.len() != n)
        {
            return false;
        }
        let mut seen = vec![false; n];
        for &index in &self.permutation {
            if index >= n || seen[index] {
                return false;
            }
            seen[index] = true;
        }
        for (i, row) in self.lower.iter().enumerate() {
            if !row[i].is_one() || row[i + 1..].iter().any(|x| !x.is_zero()) {
                return false;
            }
        }
        let mut reconstructed = vec![vec![Q::zero(); n]; n];
        let mut k = 0;
        for block in &self.blocks {
            let width = match block {
                DiagonalBlock::Single(_) => 1,
                DiagonalBlock::Pair(b) if !b.is_zero() => 2,
                DiagonalBlock::Pair(_) => return false,
            };
            if width > n.saturating_sub(k) {
                return false;
            }
            for (i, row) in reconstructed.iter_mut().enumerate() {
                for (j, entry) in row.iter_mut().enumerate() {
                    *entry += match block {
                        DiagonalBlock::Single(d) => &self.lower[i][k] * d * &self.lower[j][k],
                        DiagonalBlock::Pair(b) => {
                            b * (&self.lower[i][k] * &self.lower[j][k + 1]
                                + &self.lower[i][k + 1] * &self.lower[j][k])
                        }
                    };
                }
            }
            k += width;
        }
        k == n
            && reconstructed.iter().enumerate().all(|(i, row)| {
                row.iter()
                    .enumerate()
                    .all(|(j, entry)| entry == &matrix[self.permutation[i]][self.permutation[j]])
            })
    }
}

fn swap_basis(a: &mut [Vec<Q>], certificate: &mut SymmetricLdlt, i: usize, j: usize, k: usize) {
    a.swap(i, j);
    for row in a {
        row.swap(i, j);
    }
    certificate.permutation.swap(i, j);
    for column in 0..k {
        let value = certificate.lower[i][column].clone();
        certificate.lower[i][column] = certificate.lower[j][column].clone();
        certificate.lower[j][column] = value;
    }
}

/// Compute a block LDLᵀ certificate using exact rational arithmetic.
///
/// Returns a typed error for nonsquare or nonsymmetric input. The empty matrix
/// has an empty certificate. Cost: O(n³) rational operations and O(n²) storage.
// Paired row/column updates preserve symmetry in place; indexed traversal
// keeps these congruence updates together without cloning the trailing block.
#[allow(clippy::needless_range_loop)]
pub fn symmetric_ldlt(matrix: &[Vec<Q>]) -> Result<SymmetricLdlt, InertiaError> {
    validate(matrix)?;
    let n = matrix.len();
    let mut a = matrix.to_vec();
    let mut certificate = SymmetricLdlt {
        permutation: (0..n).collect(),
        lower: vec![vec![Q::zero(); n]; n],
        blocks: Vec::new(),
    };
    for i in 0..n {
        certificate.lower[i][i] = Q::one();
    }
    let mut k = 0;
    while k < n {
        if let Some(pivot) = (k..n).find(|&i| !a[i][i].is_zero()) {
            swap_basis(&mut a, &mut certificate, k, pivot, k);
            let d = a[k][k].clone();
            for i in k + 1..n {
                certificate.lower[i][k] = &a[i][k] / &d;
            }
            for i in k + 1..n {
                for j in i..n {
                    let updated = &a[i][j] - &certificate.lower[i][k] * &a[j][k];
                    a[i][j] = updated.clone();
                    a[j][i] = updated;
                }
            }
            certificate.blocks.push(DiagonalBlock::Single(d));
            k += 1;
        } else if let Some((i, j)) =
            (k..n).find_map(|i| (i + 1..n).find(|&j| !a[i][j].is_zero()).map(|j| (i, j)))
        {
            swap_basis(&mut a, &mut certificate, k, i, k);
            swap_basis(&mut a, &mut certificate, k + 1, j, k);
            let b = a[k][k + 1].clone();
            for i in k + 2..n {
                certificate.lower[i][k] = &a[i][k + 1] / &b;
                certificate.lower[i][k + 1] = &a[i][k] / &b;
            }
            for i in k + 2..n {
                for j in i..n {
                    let updated = &a[i][j]
                        - &certificate.lower[i][k] * &a[j][k]
                        - &certificate.lower[i][k + 1] * &a[j][k + 1];
                    a[i][j] = updated.clone();
                    a[j][i] = updated;
                }
            }
            certificate.blocks.push(DiagonalBlock::Pair(b));
            k += 2;
        } else {
            certificate
                .blocks
                .extend((k..n).map(|_| DiagonalBlock::Single(Q::zero())));
            break;
        }
    }
    Ok(certificate)
}

/// Count eigenvalue signs exactly, including singular and empty matrices.
pub fn symmetric_inertia(matrix: &[Vec<Q>]) -> Result<Inertia, InertiaError> {
    Ok(symmetric_ldlt(matrix)?.inertia())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix(rows: &[&[i64]]) -> Vec<Vec<Q>> {
        rows.iter()
            .map(|row| row.iter().map(|&x| Q::from_integer(x.into())).collect())
            .collect()
    }

    #[test]
    fn pivots_singular_shapes_and_replay() {
        for (a, expected) in [
            (matrix(&[]), (0, 0, 0)),
            (matrix(&[&[0, 2], &[2, 0]]), (1, 1, 0)),
            (matrix(&[&[0, 0, 2], &[0, 0, 0], &[2, 0, 0]]), (1, 1, 1)),
            (matrix(&[&[0, 1, 2], &[1, 3, 4], &[2, 4, 5]]), (2, 1, 0)),
            (matrix(&[&[1, 1], &[1, 1]]), (1, 0, 1)),
        ] {
            let mut certificate = symmetric_ldlt(&a).unwrap();
            assert!(certificate.verify(&a));
            let inertia = certificate.inertia();
            assert_eq!((inertia.positive, inertia.negative, inertia.zero), expected);
            certificate.blocks.push(DiagonalBlock::Single(Q::zero()));
            assert!(!certificate.verify(&a));
        }
        assert_eq!(
            symmetric_ldlt(&matrix(&[&[1, 2]])),
            Err(InertiaError::NotSquare)
        );
        assert_eq!(
            symmetric_ldlt(&matrix(&[&[1, 2], &[3, 4]])),
            Err(InertiaError::NotSymmetric)
        );
    }

    #[test]
    fn congruences_preserve_known_inertia() {
        // A = T diag(d) Tᵀ for invertible unit triangular T.
        for seed in 0..40i64 {
            let n = 7;
            let mut t = vec![vec![Q::zero(); n]; n];
            for (i, row) in t.iter_mut().enumerate() {
                row[i] = Q::one();
                for (j, entry) in row.iter_mut().enumerate().take(i) {
                    *entry = Q::new(((seed + 3 * i as i64 + j as i64) % 9 - 4).into(), 3.into());
                }
            }
            let d = [-2, 0, 3, -1, 0, 4, 5].map(|x| Q::from_integer(x.into()));
            let a: Vec<Vec<Q>> = (0..n)
                .map(|i| {
                    (0..n)
                        .map(|j| (0..n).map(|k| &t[i][k] * &d[k] * &t[j][k]).sum())
                        .collect()
                })
                .collect();
            let certificate = symmetric_ldlt(&a).unwrap();
            assert!(certificate.verify(&a));
            assert_eq!(
                certificate.inertia(),
                Inertia {
                    positive: 3,
                    negative: 2,
                    zero: 2
                }
            );
        }
    }
}
