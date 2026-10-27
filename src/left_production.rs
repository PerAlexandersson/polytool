//! Left production matrices: one certificate for a triangle, its reversal and
//! the Pólya frequency property of every row (Chen--Fu--Ruan).
//!
//! Let `A = [a_{n,k}]` be lower triangular with nonzero diagonal, so row `n`
//! is a polynomial `P_n(t) = sum_k a_{n,k} t^k` of degree exactly `n`.  Its
//! *left production matrix* (Liang et al.; Chen--Fu--Ruan, arXiv:2512.08369)
//! is `Q(A) = A (1 ⊕ A^{-1})`, equivalently the lower-triangular matrix with
//!
//! ```text
//! a_{n,0} = q_{n,0},    a_{n,k} = sum_{i=k}^{n} q_{n,i} a_{i-1,k-1}   (k >= 1),
//! ```
//!
//! i.e. `A = Q (1 ⊕ A) = Q [1 ⊕ Q] [I_2 ⊕ Q] ...`.  Row `n` of `A` depends
//! only on the leading block `Q_n = Q[0..=n, 0..=n]`, and conversely.
//!
//! # The theorem
//!
//! Chen, Fu and Ruan (*A unified planar network approach to total positivity
//! of combinatorial matrices and real-rootedness of polynomials*,
//! arXiv:2512.08369, Theorems 3 and 4) build one planar network from copies
//! of a network for `Q`.  With different sources and sinks its path matrices
//! are `A`, the reversal `[a_{n,n-k}]`, and the transposed Toeplitz sections
//! of each row:
//!
//! ```text
//! T_r(α_n)^T = M_{n,r}[n..=n+r; 0..=r],
//! M_{n,r} = (Q_n ⊕ I_r)(I_1 ⊕ Q_n ⊕ I_{r-1}) ... (I_r ⊕ Q_n),
//! ```
//!
//! where `α_n = (a_{n,0}, ..., a_{n,n}, 0, ...)`.  Hence: **if `Q` is totally
//! nonnegative (TNN), then `A` and its reversal are TNN and every row is a PF
//! sequence**, i.e. real-rooted with nonpositive zeros.
//!
//! # What a finite certificate proves
//!
//! Given rows `P_0, ..., P_N`, [`left_production_matrix`] computes `Q_N`
//! exactly, and it is uniquely determined.  If `Q_N` is TNN, then extending it
//! by an identity block, `Q_N ⊕ I`, gives an infinite TNN lower-triangular
//! matrix whose triangle agrees with `A` in rows `0..=N` (row `n` only uses
//! `Q_n`).  So the theorem applies and proves, for the supplied rows only:
//! `A_N` and its reversal are TNN, and each of `P_0, ..., P_N` is PF.  The
//! Toeplitz identity is checked for every row with `r = n` (a window that
//! contains every coefficient); it is an algebraic identity, verified here as
//! a consistency check, and the PF conclusion comes from TNN of `Q_N`.
//!
//! Failure of `Q_N` to be TNN proves nothing about the rows: the criterion is
//! sufficient, not necessary.  The Eulerian triangle is an example: its rows
//! are real-rooted, but `q_{4,2} = -5`.  (Total nonnegativity of the Eulerian
//! triangle itself is Brenti's open conjecture.)
//!
//! The families covered by the paper are all instances: exponential Riordan
//! arrays `R[f', f]` with `f'` PF (Stirling numbers of both kinds, Lah
//! numbers), iteration matrices of multiplier sequences (idempotent numbers),
//! `r`-Whitney numbers, and Brenti's `n`-recursive matrices
//! `t_{n,k} = a_n t_{n-1,k-1} + b_n t_{n-1,k} + c_n t_{n-2,k-1}` with
//! nonnegative coefficients (Pascal, Delannoy, derangement triangles).  The
//! certificate is computed from the rows, so no family needs to be named.
//!
//! When `Q_N` is TNN, [`check_left_production_matrix`] also returns an explicit
//! planar network for it: `Q_N = L D` with `D` its (positive) diagonal and `L`
//! lower-unitriangular TNN, and `L` is the path matrix of the canonical network
//! of [`crate::tnn_network`].
//!
//! Conventions: triangles are in triangular row format (row `n` has `n + 1`
//! entries), zero-based, exact.  For a search over height-invariant networks
//! see [`crate::pf_network`].

use crate::tnn_network::{
    evaluate_path_matrix, find_negative_minor, reconstruct_lambda_from_q_rows, BigRational,
    CanonicalPlanarNetwork,
};
use crate::{check_tnn_neville_bigint, Polynomial};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::fmt;

type Q = BigRational;

/// Invalid input or a failed certificate check.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LeftProductionError {
    /// No rows were supplied.
    EmptyInput,
    /// `first_row_index` was not 0 or 1.
    UnsupportedFirstRowIndex {
        /// The rejected index.
        found: usize,
    },
    /// Row `row` does not have degree exactly `row` (the triangle needs a
    /// nonzero diagonal).
    DegreeMismatch {
        /// Zero-based row index in the triangle.
        row: usize,
        /// Its degree (`None` for the zero polynomial).
        found: Option<usize>,
    },
    /// A matrix row has the wrong length for triangular row format.
    InvalidRowLength {
        /// Zero-based row index.
        row: usize,
        /// Expected length `row + 1`.
        expected: usize,
        /// Actual length.
        found: usize,
    },
    /// A diagonal entry is zero.
    ZeroDiagonal {
        /// Zero-based row index.
        row: usize,
    },
    /// The triangle regenerated from the production matrix differs from the
    /// stored triangle.
    RegeneratedMismatch {
        /// First differing row.
        row: usize,
    },
    /// The production matrix is not totally nonnegative.
    NotTnn {
        /// The Neville-elimination diagnosis.
        reason: String,
    },
    /// A Toeplitz section of a row differs from the unified-network window.
    ToeplitzMismatch {
        /// The row `n`.
        row: usize,
        /// Row index of the first differing entry in the window.
        i: usize,
        /// Column index of the first differing entry in the window.
        j: usize,
    },
}

impl fmt::Display for LeftProductionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyInput => write!(f, "at least one row is required"),
            Self::UnsupportedFirstRowIndex { found } => {
                write!(f, "first row index must be 0 or 1, found {}", found)
            }
            Self::DegreeMismatch { row, found } => match found {
                Some(degree) => write!(
                    f,
                    "row {} has degree {}; a triangle with nonzero diagonal needs degree {}",
                    row, degree, row
                ),
                None => write!(f, "row {} is the zero polynomial", row),
            },
            Self::InvalidRowLength {
                row,
                expected,
                found,
            } => write!(
                f,
                "row {} has {} entries, expected {}",
                row, found, expected
            ),
            Self::ZeroDiagonal { row } => write!(f, "diagonal entry of row {} is zero", row),
            Self::RegeneratedMismatch { row } => write!(
                f,
                "the triangle regenerated from the production matrix differs in row {}",
                row
            ),
            Self::NotTnn { reason } => {
                write!(
                    f,
                    "production matrix is not totally nonnegative: {}",
                    reason
                )
            }
            Self::ToeplitzMismatch { row, i, j } => write!(
                f,
                "row {}: unified-network window differs from the Toeplitz section at ({}, {})",
                row, i, j
            ),
        }
    }
}

impl std::error::Error for LeftProductionError {}

/// Coefficient rows `P_0, ..., P_N` in triangular row format, checking that
/// `deg P_n = n`.  With `first_row_index == 1`, `P_0 = 1` is prepended.
pub fn triangle_from_polynomials(
    polys: &[Polynomial<BigInt>],
    first_row_index: usize,
) -> Result<Vec<Vec<BigInt>>, LeftProductionError> {
    if polys.is_empty() {
        return Err(LeftProductionError::EmptyInput);
    }
    if first_row_index > 1 {
        return Err(LeftProductionError::UnsupportedFirstRowIndex {
            found: first_row_index,
        });
    }
    let mut rows = Vec::with_capacity(polys.len() + first_row_index);
    if first_row_index == 1 {
        rows.push(vec![BigInt::one()]);
    }
    for poly in polys {
        let row = rows.len();
        if poly.degree() != Some(row) {
            return Err(LeftProductionError::DegreeMismatch {
                row,
                found: poly.degree(),
            });
        }
        rows.push(poly.coeffs().to_vec());
    }
    Ok(rows)
}

fn validate_triangle<T: Zero>(rows: &[Vec<T>]) -> Result<(), LeftProductionError> {
    if rows.is_empty() {
        return Err(LeftProductionError::EmptyInput);
    }
    for (n, row) in rows.iter().enumerate() {
        if row.len() != n + 1 {
            return Err(LeftProductionError::InvalidRowLength {
                row: n,
                expected: n + 1,
                found: row.len(),
            });
        }
        if row[n].is_zero() {
            return Err(LeftProductionError::ZeroDiagonal { row: n });
        }
    }
    Ok(())
}

/// The leading block `Q_N` of the left production matrix of the triangle
/// `rows` (triangular row format, nonzero diagonal), in triangular row format.
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::left_production::left_production_matrix;
/// use polytool::tnn_network::BigRational;
///
/// // Pascal's triangle has the all-ones lower-triangular production matrix.
/// let pascal: Vec<Vec<BigInt>> = vec![vec![1], vec![1, 1], vec![1, 2, 1], vec![1, 3, 3, 1]]
///     .into_iter()
///     .map(|r: Vec<i64>| r.into_iter().map(BigInt::from).collect())
///     .collect();
/// let q = left_production_matrix(&pascal).unwrap();
/// assert!(q.iter().flatten().all(|v| *v == BigRational::from_integer(BigInt::from(1))));
/// ```
pub fn left_production_matrix(
    rows: &[Vec<BigInt>],
) -> Result<Vec<Vec<BigRational>>, LeftProductionError> {
    validate_triangle(rows)?;
    let a: Vec<Vec<Q>> = rows
        .iter()
        .map(|row| row.iter().map(|c| Q::from_integer(c.clone())).collect())
        .collect();
    let mut q: Vec<Vec<Q>> = Vec::with_capacity(a.len());
    for n in 0..a.len() {
        let mut row = vec![Q::zero(); n + 1];
        row[0] = a[n][0].clone();
        for k in (1..=n).rev() {
            let mut value = a[n][k].clone();
            for i in (k + 1)..=n {
                value -= &row[i] * &a[i - 1][k - 1];
            }
            row[k] = value / &a[k - 1][k - 1];
        }
        q.push(row);
    }
    Ok(q)
}

/// Regenerate the triangle from a left production matrix:
/// `a_{n,0} = q_{n,0}` and `a_{n,k} = sum_{i=k}^{n} q_{n,i} a_{i-1,k-1}`.
/// `production` is in triangular row format; missing entries count as zero.
pub fn triangle_from_left_production_matrix(
    production: &[Vec<BigRational>],
) -> Vec<Vec<BigRational>> {
    let entry = |n: usize, i: usize| production[n].get(i).cloned().unwrap_or_else(Q::zero);
    let mut a: Vec<Vec<Q>> = Vec::with_capacity(production.len());
    for n in 0..production.len() {
        let mut row = vec![Q::zero(); n + 1];
        row[0] = entry(n, 0);
        for k in 1..=n {
            let mut value = Q::zero();
            for i in k..=n {
                value += entry(n, i) * &a[i - 1][k - 1];
            }
            row[k] = value;
        }
        a.push(row);
    }
    a
}

/// The window `M_{n,r}[n..=n+r; 0..=r]` of the unified Chen--Fu--Ruan
/// network, where `M_{n,r} = prod_{i=0}^{r} (I_i ⊕ Q_n ⊕ I_{r-i})` and `Q_n`
/// is the leading `(n+1) x (n+1)` block of `production` (triangular row
/// format, at least `n + 1` rows).  By Theorem 4 of arXiv:2512.08369 it equals
/// the transposed Toeplitz section `[a_{n, j-i}]_{0 <= i, j <= r}` of row `n`.
///
/// # Panics
///
/// Panics if `production` has at most `n` rows.
pub fn unified_network_toeplitz_window(
    production: &[Vec<BigRational>],
    n: usize,
    r: usize,
) -> Vec<Vec<BigRational>> {
    assert!(production.len() > n, "production matrix has too few rows");
    let entry = |a: usize, c: usize| production[a].get(c).cloned().unwrap_or_else(Q::zero);
    let size = n + r + 1;
    (n..=n + r)
        .map(|start| {
            let mut v = vec![Q::zero(); size];
            v[start] = Q::one();
            for shift in 0..=r {
                // v <- v (I_shift ⊕ Q_n ⊕ I_{r-shift}); only columns
                // shift..=shift+n change.
                let block: Vec<Q> = (0..=n)
                    .map(|c| {
                        let mut total = Q::zero();
                        for a in c..=n {
                            let x = &v[shift + a];
                            if !x.is_zero() {
                                total += x * entry(a, c);
                            }
                        }
                        total
                    })
                    .collect();
                for (c, value) in block.into_iter().enumerate() {
                    v[shift + c] = value;
                }
            }
            v.truncate(r + 1);
            v
        })
        .collect()
}

/// A negative minor of a rational matrix: zero-based increasing `rows` and
/// `cols` of the square (zero-padded) matrix, and its exact value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RationalNegativeMinor {
    /// Row indices (strictly increasing).
    pub rows: Vec<usize>,
    /// Column indices (strictly increasing, same length as `rows`).
    pub cols: Vec<usize>,
    /// The exact (negative) value.
    pub value: BigRational,
}

/// Scale a rational triangular matrix to a square integer matrix by the
/// positive lcm of its denominators.
fn scaled_square(production: &[Vec<Q>]) -> (BigInt, Vec<Vec<BigInt>>) {
    let size = production.len();
    let scale = production
        .iter()
        .flatten()
        .fold(BigInt::one(), |acc, v| acc.lcm(v.denom()));
    let matrix = (0..size)
        .map(|n| {
            (0..size)
                .map(|k| {
                    production[n].get(k).map_or_else(BigInt::zero, |v| {
                        (v * Q::from_integer(scale.clone())).to_integer()
                    })
                })
                .collect()
        })
        .collect();
    (scale, matrix)
}

/// Why a production matrix is not totally nonnegative.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProductionTnnFailure {
    /// The Neville-elimination diagnosis.
    pub reason: String,
    /// An explicit negative minor, when the bounded search finds one.
    pub witness: Option<Box<RationalNegativeMinor>>,
}

/// Decide exactly whether a rational lower-triangular matrix (triangular row
/// format) is totally nonnegative.  The matrix is scaled by the positive lcm
/// of its denominators and tested by exact Neville elimination
/// ([`check_tnn_neville_bigint`]); on failure a bounded search looks for an
/// explicit negative minor, whose value is reported for the unscaled matrix.
pub fn check_production_tnn(production: &[Vec<BigRational>]) -> Result<(), ProductionTnnFailure> {
    let (scale, matrix) = scaled_square(production);
    check_tnn_neville_bigint(&matrix).map_err(|error| {
        let reason = error.to_string();
        let witness = find_negative_minor(&matrix).map(|minor| {
            let size = minor.rows.len() as u32;
            Box::new(RationalNegativeMinor {
                value: Q::new(minor.value, scale.pow(size)),
                rows: minor.rows,
                cols: minor.cols,
            })
        });
        ProductionTnnFailure { reason, witness }
    })
}

/// A claimed left-production certificate: the triangle and its production
/// matrix.  [`LeftProductionCertificate::verify`] re-checks everything.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LeftProductionCertificate {
    /// The triangle `A_N` in triangular row format.
    pub triangle: Vec<Vec<BigInt>>,
    /// The production matrix `Q_N` in triangular row format.
    pub production_matrix: Vec<Vec<BigRational>>,
}

impl LeftProductionCertificate {
    /// Check that the triangle has a nonzero diagonal, that it regenerates
    /// from the production matrix, that the production matrix is TNN, and
    /// that for every row `n` the unified-network window
    /// [`unified_network_toeplitz_window`] with `r = n` equals the transposed
    /// Toeplitz section of the row.
    pub fn verify(&self) -> Result<(), LeftProductionError> {
        validate_triangle(&self.triangle)?;
        if self.production_matrix.len() != self.triangle.len() {
            return Err(LeftProductionError::RegeneratedMismatch {
                row: self.production_matrix.len().min(self.triangle.len()),
            });
        }
        let regenerated = triangle_from_left_production_matrix(&self.production_matrix);
        for (n, (got, want)) in regenerated.iter().zip(&self.triangle).enumerate() {
            let want: Vec<Q> = want.iter().map(|c| Q::from_integer(c.clone())).collect();
            if *got != want {
                return Err(LeftProductionError::RegeneratedMismatch { row: n });
            }
        }
        check_production_tnn(&self.production_matrix).map_err(|failure| {
            LeftProductionError::NotTnn {
                reason: failure.reason,
            }
        })?;
        for (n, row) in self.triangle.iter().enumerate() {
            let window = unified_network_toeplitz_window(&self.production_matrix, n, n);
            for (i, got) in window.iter().enumerate() {
                for (j, value) in got.iter().enumerate() {
                    let want = j
                        .checked_sub(i)
                        .and_then(|d| row.get(d))
                        .map_or_else(Q::zero, |c| Q::from_integer(c.clone()));
                    if *value != want {
                        return Err(LeftProductionError::ToeplitzMismatch { row: n, i, j });
                    }
                }
            }
        }
        Ok(())
    }
}

/// An explicit planar network for a TNN production matrix: `Q_N = L D`, where
/// `L` is the path matrix of the canonical network `network` (see
/// [`crate::tnn_network::CanonicalPlanarNetwork`]) and `D = diag(column_scales)`
/// multiplies the weight of every path ending at sink `k` by
/// `column_scales[k] = q_{k,k} > 0`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProductionNetwork {
    /// Canonical network of the unit-normalized production matrix `L`.
    pub network: CanonicalPlanarNetwork<BigRational>,
    /// The diagonal of `Q_N`.
    pub column_scales: Vec<BigRational>,
}

/// Result of [`check_left_production_matrix`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LeftProductionReport {
    /// The triangle `A_N` (triangular row format; `P_0 = 1` included when it
    /// was prepended).
    pub triangle: Vec<Vec<BigInt>>,
    /// The leading block `Q_N` of the left production matrix.
    pub production_matrix: Vec<Vec<BigRational>>,
    /// Whether `Q_N` is totally nonnegative.
    pub tnn: bool,
    /// When not TNN: the Neville-elimination diagnosis.
    pub failure_reason: Option<String>,
    /// When not TNN and found: an explicit negative minor of `Q_N`.
    pub negative_minor: Option<RationalNegativeMinor>,
    /// When TNN: whether [`LeftProductionCertificate::verify`] passed.
    pub certificate_verified: Option<bool>,
    /// When TNN: an explicit planar network for `Q_N`, re-evaluated exactly.
    pub production_network: Option<ProductionNetwork>,
}

impl LeftProductionReport {
    /// Whether the report certifies (Chen--Fu--Ruan) that every row
    /// `P_0, ..., P_N` is PF and that `A_N` and its reversal are TNN.
    pub fn rows_certified_pf(&self) -> bool {
        self.tnn && self.certificate_verified == Some(true)
    }
}

impl fmt::Display for LeftProductionReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let last = self.triangle.len() - 1;
        if self.rows_certified_pf() {
            write!(
                f,
                "certified: Q_{} is TNN, so rows 0..={} are PF (real-rooted, nonpositive zeros) \
                 and the triangle and its reversal are TNN (Chen--Fu--Ruan); this prefix only",
                last, last
            )
        } else if self.tnn {
            write!(
                f,
                "Q_{} is TNN but the certificate failed re-verification",
                last
            )
        } else {
            write!(
                f,
                "not certified: Q_{} is not TNN ({}); this says nothing about the rows",
                last,
                self.failure_reason.as_deref().unwrap_or("unknown reason")
            )
        }
    }
}

fn production_network(production: &[Vec<Q>]) -> Option<ProductionNetwork> {
    let column_scales: Vec<Q> = production
        .iter()
        .enumerate()
        .map(|(k, row)| row[k].clone())
        .collect();
    if column_scales.iter().any(|d| !d.is_positive()) {
        return None;
    }
    let unit: Vec<Vec<Q>> = production
        .iter()
        .map(|row| row.iter().zip(&column_scales).map(|(v, d)| v / d).collect())
        .collect();
    // Zero-tail condition of the canonical network (holds for TNN input).
    if let Some(start) = (1..unit.len()).find(|&n| unit[n][0].is_zero()) {
        if unit[start..].iter().any(|row| !row[0].is_zero()) {
            return None;
        }
    }
    let lambda = reconstruct_lambda_from_q_rows(&unit).ok()?;
    let network = CanonicalPlanarNetwork { lambda };
    (evaluate_path_matrix(&network) == unit).then_some(ProductionNetwork {
        network,
        column_scales,
    })
}

/// Compute the left production matrix of the rows `P_0, ..., P_N` (with
/// `first_row_index == 1`, the input is `P_1, ..., P_N` and `P_0 = 1` is
/// prepended) and decide whether it certifies PF rows.
///
/// Requires `deg P_n = n`.  See the module documentation for exactly what a
/// positive answer proves.
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::left_production::check_left_production_matrix;
/// use polytool::Polynomial;
///
/// // Bell polynomials sum_k S(n,k) t^k (Stirling numbers of the second kind).
/// let rows: Vec<Polynomial<BigInt>> = [vec![0, 1], vec![0, 1, 1], vec![0, 1, 3, 1]]
///     .iter()
///     .map(|r| Polynomial::from_i64_coeffs(r))
///     .collect();
/// let report = check_left_production_matrix(&rows, 1).unwrap();
/// assert!(report.rows_certified_pf());
///
/// // 1 + t + t^2 is not certified.
/// let rows = vec![Polynomial::from_i64_coeffs(&[1, 1]), Polynomial::from_i64_coeffs(&[1, 1, 1])];
/// assert!(!check_left_production_matrix(&rows, 1).unwrap().rows_certified_pf());
/// ```
pub fn check_left_production_matrix(
    polys: &[Polynomial<BigInt>],
    first_row_index: usize,
) -> Result<LeftProductionReport, LeftProductionError> {
    let triangle = triangle_from_polynomials(polys, first_row_index)?;
    let production_matrix = left_production_matrix(&triangle)?;
    let verdict = check_production_tnn(&production_matrix);
    let tnn = verdict.is_ok();
    let (failure_reason, negative_minor) = match verdict {
        Ok(()) => (None, None),
        Err(failure) => (Some(failure.reason), failure.witness.map(|minor| *minor)),
    };
    let (certificate_verified, network) = if tnn {
        let certificate = LeftProductionCertificate {
            triangle: triangle.clone(),
            production_matrix: production_matrix.clone(),
        };
        (
            Some(certificate.verify().is_ok()),
            production_network(&production_matrix),
        )
    } else {
        (None, None)
    };
    Ok(LeftProductionReport {
        triangle,
        production_matrix,
        tnn,
        failure_reason,
        negative_minor,
        certificate_verified,
        production_network: network,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numeric::binomial_bigint;

    fn q(v: i64) -> Q {
        Q::from_integer(BigInt::from(v))
    }

    fn triangle(rows: &[&[i64]]) -> Vec<Vec<BigInt>> {
        rows.iter()
            .map(|r| r.iter().map(|&c| BigInt::from(c)).collect())
            .collect()
    }

    fn polys(rows: &[Vec<BigInt>]) -> Vec<Polynomial<BigInt>> {
        rows.iter().map(|r| Polynomial::new(r.clone())).collect()
    }

    fn pascal(last: usize) -> Vec<Vec<BigInt>> {
        (0..=last)
            .map(|n| (0..=n).map(|k| binomial_bigint(n, k)).collect())
            .collect()
    }

    fn eulerian(last: usize) -> Vec<Vec<BigInt>> {
        let mut rows: Vec<Vec<BigInt>> = vec![vec![BigInt::one()]];
        for m in 1..=last {
            let prev = &rows[m - 1];
            let row = (0..=m)
                .map(|k| {
                    let left = if k >= 1 {
                        BigInt::from(m - k + 1) * &prev[k - 1]
                    } else {
                        BigInt::zero()
                    };
                    left + prev
                        .get(k)
                        .map_or_else(BigInt::zero, |c| BigInt::from(k + 1) * c)
                })
                .collect();
            rows.push(row);
        }
        rows
    }

    fn stirling2(last: usize) -> Vec<Vec<BigInt>> {
        let mut rows: Vec<Vec<BigInt>> = vec![vec![BigInt::one()]];
        for m in 1..=last {
            let prev = &rows[m - 1];
            let row = (0..=m)
                .map(|k| {
                    if k == 0 {
                        return BigInt::zero();
                    }
                    prev[k - 1].clone()
                        + prev
                            .get(k)
                            .map_or_else(BigInt::zero, |c| BigInt::from(k) * c)
                })
                .collect();
            rows.push(row);
        }
        rows
    }

    fn narayana(last: usize) -> Vec<Vec<BigInt>> {
        // N(n+1, k+1) = C(n+1, k+1) C(n+1, k) / (n+1).
        (0..=last)
            .map(|n| {
                (0..=n)
                    .map(|k| {
                        binomial_bigint(n + 1, k + 1) * binomial_bigint(n + 1, k)
                            / BigInt::from(n + 1)
                    })
                    .collect()
            })
            .collect()
    }

    #[test]
    fn pascal_production_and_certificate() {
        let a = pascal(8);
        let production = left_production_matrix(&a).unwrap();
        assert!(production.iter().flatten().all(|v| *v == q(1)));
        let report = check_left_production_matrix(&polys(&a), 0).unwrap();
        assert!(report.rows_certified_pf(), "{}", report);
        let network = report.production_network.unwrap();
        assert!(network.column_scales.iter().all(|d| *d == q(1)));
    }

    #[test]
    fn stirling_rows_certified() {
        let a = stirling2(9);
        let report = check_left_production_matrix(&polys(&a), 0).unwrap();
        assert!(report.rows_certified_pf(), "{}", report);
        // Q = 1 ⊕ Pascal.
        assert_eq!(report.production_matrix[3], vec![q(0), q(1), q(2), q(1)]);
    }

    #[test]
    fn eulerian_and_narayana_not_certified() {
        let report = check_left_production_matrix(&polys(&eulerian(3)), 0).unwrap();
        assert!(report.rows_certified_pf(), "{}", report);
        let report = check_left_production_matrix(&polys(&eulerian(6)), 0).unwrap();
        assert!(!report.tnn);
        assert_eq!(report.production_matrix[4][2], q(-5));
        let minor = report.negative_minor.expect("a witness");
        assert!(minor.value < Q::zero());

        let report = check_left_production_matrix(&polys(&narayana(3)), 0).unwrap();
        assert!(report.rows_certified_pf(), "{}", report);
        // Q_4 = [.., [1, 3, 2, 4, 1]] has the minor 2 * 2 - 3 * 3 < 0, and
        // q_{5,1} = -1.
        let report = check_left_production_matrix(&polys(&narayana(6)), 0).unwrap();
        assert!(!report.tnn);
        assert_eq!(
            report.production_matrix[4],
            vec![q(1), q(3), q(2), q(4), q(1)]
        );
        assert_eq!(report.production_matrix[5][1], q(-1));
    }

    #[test]
    fn one_plus_t_plus_t_squared_never_certified() {
        let a = triangle(&[&[1], &[1, 1], &[1, 1, 1]]);
        let report = check_left_production_matrix(&polys(&a), 0).unwrap();
        assert!(!report.rows_certified_pf());
        assert_eq!(report.production_matrix[2], vec![q(1), q(0), q(1)]);
        // Any first row: [c, 1] followed by 1 + t + t^2.
        for c in 0..6 {
            let a = triangle(&[&[1], &[c, 1], &[1, 1, 1]]);
            let report = check_left_production_matrix(&polys(&a), 0).unwrap();
            assert!(!report.rows_certified_pf(), "c = {}", c);
        }
    }

    #[test]
    fn toeplitz_identity_holds_for_any_triangle() {
        // Theorem 4 of arXiv:2512.08369 is algebraic: check it on non-TNN data.
        for a in [
            eulerian(5),
            narayana(5),
            triangle(&[&[2], &[-1, 3], &[4, 0, 5]]),
        ] {
            let production = left_production_matrix(&a).unwrap();
            assert_eq!(
                triangle_from_left_production_matrix(&production),
                a.iter()
                    .map(|r| r
                        .iter()
                        .map(|c| Q::from_integer(c.clone()))
                        .collect::<Vec<_>>())
                    .collect::<Vec<_>>()
            );
            for (n, row) in a.iter().enumerate() {
                for r in 0..=n + 1 {
                    let window = unified_network_toeplitz_window(&production, n, r);
                    for (i, w) in window.iter().enumerate() {
                        for (j, value) in w.iter().enumerate() {
                            let want = j
                                .checked_sub(i)
                                .and_then(|d| row.get(d))
                                .map_or_else(Q::zero, |c| Q::from_integer(c.clone()));
                            assert_eq!(*value, want, "n {} r {} ({}, {})", n, r, i, j);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn three_term_n_recursive_triangle_certified() {
        // Delannoy: d(n,k) = d(n-1,k-1) + d(n-1,k) + d(n-2,k-1).
        let mut rows: Vec<Vec<BigInt>> = vec![vec![BigInt::one()]];
        for n in 1..=7usize {
            let row = (0..=n)
                .map(|k| {
                    let get = |m: usize, j: Option<usize>| {
                        j.and_then(|j| rows.get(m).and_then(|r: &Vec<BigInt>| r.get(j)))
                            .cloned()
                            .unwrap_or_else(BigInt::zero)
                    };
                    get(n - 1, k.checked_sub(1))
                        + get(n - 1, Some(k))
                        + n.checked_sub(2)
                            .map_or_else(BigInt::zero, |m| get(m, k.checked_sub(1)))
                })
                .collect();
            rows.push(row);
        }
        let report = check_left_production_matrix(&polys(&rows), 0).unwrap();
        assert!(report.rows_certified_pf(), "{}", report);
    }

    #[test]
    fn input_errors() {
        assert_eq!(
            check_left_production_matrix(&[], 1),
            Err(LeftProductionError::EmptyInput)
        );
        let rows = vec![Polynomial::from_i64_coeffs(&[1, 1, 1])];
        assert_eq!(
            check_left_production_matrix(&rows, 1),
            Err(LeftProductionError::DegreeMismatch {
                row: 1,
                found: Some(2)
            })
        );
        assert!(matches!(
            check_left_production_matrix(&rows, 2),
            Err(LeftProductionError::UnsupportedFirstRowIndex { found: 2 })
        ));
        let tampered = LeftProductionCertificate {
            triangle: pascal(2),
            production_matrix: vec![vec![q(1)], vec![q(1), q(1)], vec![q(1), q(2), q(1)]],
        };
        assert_eq!(
            tampered.verify(),
            Err(LeftProductionError::RegeneratedMismatch { row: 2 })
        );
    }
}
