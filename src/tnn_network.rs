//! Canonical planar-network certificates for lower-unitriangular TNN matrices.
//!
//! This module implements the constructive converse to the LGV lemma in the
//! lower-unitriangular setting: if a finite lower-unitriangular matrix is totally
//! nonnegative, then it is the path matrix of a unique triangular planar network
//! `Gamma_N(lambda)` satisfying the usual zero-propagation condition.
//!
//! Conventions: matrices are exact `BigInt` matrices in *triangular row format*
//! (row `n` has `n + 1` entries, `rows[n][k]` with `0 <= k <= n`), all indices are
//! zero-based, and polynomial coefficient vectors are ascending.  Verdicts are exact
//! (Neville elimination over `BigInt`, rational path-matrix evaluation) and
//! failures are returned as [`NetworkError`] values rather than panics.
//! For the underlying TNN test see [`crate::linalg::check_totally_nonnegative_neville_bigint`];
//! for planar strip digraphs with height-dependent transitions see
//! [`crate::brenti_sequence`].
//!
//! The main intended workflow is:
//!
//! 1. start from a monic polynomial sequence,
//! 2. extract its coefficient matrix,
//! 3. check total non-negativity,
//! 4. reconstruct the canonical network,
//! 5. verify the certificate by re-evaluating the path matrix.
//!
//! When a matrix is rejected, [`NetworkError`] says why: a malformed shape
//! (row length, nonzero entry above the diagonal, non-unit diagonal), a
//! negative entry, or a failure of total non-negativity.  In the last case
//! [`check_lower_unitriangular_totally_nonnegative`] also searches for an explicit
//! [`NegativeMinor`] whose exact value can be recomputed independently.
//!
//! # Chain polynomials
//!
//! For a lower-unitriangular matrix `R = (r_{n,k})`, the *chain polynomials*
//! are `p_0 = 1` and `p_n(t) = t * sum_{k<n} r_{n,k} p_k(t)`
//! ([`chain_polynomials_bigint`], [`chain_polynomials_rational`]).
//! Brändén and Saud Maia Leite (arXiv:2412.06595, Theorem 3.7) prove that if
//! `R` is totally nonnegative, then every `p_n` is
//! real-rooted with all zeros in `[-1, 0]` and `p_n` interlaces `p_{n+1}`.
//! [`check_chain_polynomials_bigint`] reports the TNN verdict and the
//! chain-polynomial properties separately.
//!
//! # Guardrail
//!
//! A network certificate proves TNN of the *matrix*.  It does **not** prove that
//! the row polynomials `sum_k r_{n,k} t^k` are real-rooted: the matrix
//! `[[1,0,0],[1,1,0],[1,1,1]]` is lower-unitriangular and TNN, yet its last row
//! is `1 + t + t^2`.  Real-rootedness of a single row is the
//! Aissen--Schoenberg--Whitney criterion (its Toeplitz matrix is TNN), which is
//! a different matrix.  The theorem above concerns the chain polynomials, not
//! the rows.

use crate::linalg::determinant;
use crate::real_rootedness::{check_weak_interlacing_bigint_coeffs, is_real_rooted_bigint_coeffs};
use crate::root_count::count_root_signs_prs_bigint_coeffs;
use crate::{check_totally_nonnegative_neville_bigint, vec_poly, Polynomial};
use num_bigint::BigInt;
use num_traits::{One, Zero};
use std::fmt;

/// Exact rational type used for canonical network weights.
pub use crate::BigRational;

/// Canonical triangular planar network `Gamma_N(lambda)`.
///
/// The network has `num_rows = lambda.len() + 1` source rows. Horizontal edges
/// have weight `1`, and `lambda[n][k]` is the weight of the vertical edge
/// `(n + 1, k) -> (n, k)`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CanonicalPlanarNetwork<T> {
    /// Vertical weights `lambda[n][k]` for `0 <= k <= n < num_rows - 1`.
    pub lambda: Vec<Vec<T>>,
}

impl<T> CanonicalPlanarNetwork<T> {
    /// Number of source rows in the associated triangular network.
    pub fn num_rows(&self) -> usize {
        self.lambda.len() + 1
    }

    /// Access a vertical weight.
    pub fn vertical_weight(&self, n: usize, k: usize) -> Option<&T> {
        self.lambda.get(n).and_then(|row| row.get(k))
    }
}

/// Canonical certificate consisting of a coefficient matrix and its network.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CanonicalTotallyNonnegativeProof {
    /// The lower-unitriangular coefficient matrix in triangular row format:
    /// `coefficient_matrix[n][k]` is the coefficient of `t^k` in the `n`-th
    /// polynomial of the input sequence.
    pub coefficient_matrix: Vec<Vec<BigInt>>,
    /// The canonical network whose path matrix equals `coefficient_matrix`.
    pub network: CanonicalPlanarNetwork<BigRational>,
}

/// Deprecated name of [`CanonicalTotallyNonnegativeProof`].
#[deprecated(
    since = "0.6.0",
    note = "renamed to `CanonicalTotallyNonnegativeProof`"
)]
pub type CanonicalTnnProof = CanonicalTotallyNonnegativeProof;

impl CanonicalTotallyNonnegativeProof {
    /// Recompute the path matrix of the stored network.
    ///
    /// Returns `None` for a malformed network (some `lambda[n]` without exactly
    /// `n + 1` entries; the fields are public, so a certificate can be built by
    /// hand), the same condition [`verify`](Self::verify) rejects.
    pub fn evaluated_path_matrix(&self) -> Option<Vec<Vec<BigRational>>> {
        is_well_formed_network(&self.network).then(|| evaluate_path_matrix(&self.network))
    }

    /// Check that the stored network evaluates back to the stored matrix.
    pub fn verify(&self) -> bool {
        verify_path_matrix_certificate(&self.coefficient_matrix, &self.network)
    }
}

/// An explicit minor with a negative value: a checkable witness that a matrix
/// is not totally nonnegative.
///
/// `rows` and `cols` are increasing, zero-based indices into the square
/// (zero-padded) form of the lower-triangular matrix, and `value` is the exact
/// determinant of the corresponding submatrix.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NegativeMinor {
    /// Row indices of the minor (strictly increasing, zero-based).
    pub rows: Vec<usize>,
    /// Column indices of the minor (strictly increasing, zero-based; same
    /// length as `rows`).
    pub cols: Vec<usize>,
    /// The exact (negative) value of the minor.
    pub value: BigInt,
}

impl NegativeMinor {
    /// Recompute the minor in `matrix` (rows of any length; missing entries
    /// count as zero) and check that it equals the stored negative value.
    pub fn verify(&self, matrix: &[Vec<BigInt>]) -> bool {
        if self.rows.len() != self.cols.len() || self.value >= BigInt::zero() {
            return false;
        }
        if self.rows.iter().any(|&r| r >= matrix.len()) {
            return false;
        }
        determinant(&submatrix(matrix, &self.rows, &self.cols)) == self.value
    }
}

impl fmt::Display for NegativeMinor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "minor with rows {:?} and columns {:?} equals {}",
            self.rows, self.cols, self.value
        )
    }
}

/// Errors returned while building or verifying a TNN network certificate.
///
/// The shape variants (`InvalidRowLength`, `NotSquare`,
/// `NonzeroAboveDiagonal`, `NonUnitDiagonal`) and `NegativeEntry` name the
/// offending entry.  `NotTotallyNonnegative` carries the Neville-elimination diagnosis and,
/// when the bounded search finds one, an explicit [`NegativeMinor`].
/// `NegativeReducedEntry` names the failing Whitney-reduction step.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NetworkError {
    /// The matrix has no rows.
    EmptyMatrix,
    /// The polynomial sequence is empty.
    EmptyPolynomialSequence,
    /// A square input matrix has a row of the wrong length.
    NotSquare {
        /// Zero-based index of the offending row.
        row: usize,
        /// Required row length (the number of rows).
        expected: usize,
        /// Actual row length.
        found: usize,
    },
    /// A square input matrix has a nonzero entry above the diagonal.
    NonzeroAboveDiagonal {
        /// Row of the first offending entry.
        row: usize,
        /// Column of the first offending entry (`col > row`).
        col: usize,
        /// The nonzero value found there.
        value: BigInt,
    },
    /// The first polynomial of a monic sequence has degree other than 0 or 1.
    UnsupportedStartingDegree {
        /// Degree of the first supplied polynomial.
        found: usize,
    },
    /// The polynomial at a position of the sequence does not have degree equal
    /// to its position (the sequence must have degrees `0, 1, 2, ...`).
    UnexpectedPolynomialDegree {
        /// Position in the sequence (counting the implicit `P_0 = 1` when the
        /// input started at degree 1).
        index: usize,
        /// Required degree.
        expected: usize,
        /// Actual degree, `None` for the zero polynomial.
        found: Option<usize>,
    },
    /// A polynomial of the sequence is not monic.
    NonMonicPolynomial {
        /// Position in the sequence (same convention as
        /// [`NetworkError::UnexpectedPolynomialDegree`]).
        index: usize,
        /// Its leading coefficient (`None` for the zero polynomial).
        leading: Option<BigInt>,
    },
    /// A triangular-format row does not have `row + 1` entries.
    InvalidRowLength {
        /// Zero-based row index.
        row: usize,
        /// Required length, `row + 1`.
        expected: usize,
        /// Actual length.
        found: usize,
    },
    /// A matrix entry is negative.
    NegativeEntry {
        /// Row of the entry.
        row: usize,
        /// Column of the entry.
        col: usize,
        /// The negative value.
        value: BigInt,
    },
    /// A diagonal entry is not 1.  For rational input the value is the integer
    /// part of the entry (see [`chain_polynomials_rational`]).
    NonUnitDiagonal {
        /// Row (and column) of the diagonal entry.
        row: usize,
        /// Its value.
        value: BigInt,
    },
    /// Neville elimination proved that the matrix is not totally
    /// nonnegative.  `reason` describes the failing elimination step;
    /// `witness` is an explicit negative minor when the bounded search finds
    /// one (it always does for matrices with at most
    /// [`EXHAUSTIVE_MINOR_SEARCH_MAX_DIMENSION`] rows).
    NotTotallyNonnegative {
        /// Message of the failing Neville-elimination step.
        reason: String,
        /// An explicit negative minor, when the bounded search finds one.
        witness: Option<NegativeMinor>,
    },
    /// Once the first column of a TNN lower-unitriangular matrix vanishes it
    /// stays zero; this row has a nonzero first-column entry after a zero.
    ZeroTailViolation {
        /// Row whose first entry is nonzero after an earlier zero.
        row: usize,
        /// That nonzero first-column entry.
        value: BigInt,
    },
    /// A Whitney reduction step produced a negative weight (cannot occur for a
    /// totally nonnegative input).
    NegativeReducedEntry {
        /// Row of the reduced matrix.
        row: usize,
        /// Column of the reduced matrix.
        col: usize,
        /// The negative reduced value.
        value: BigRational,
    },
    /// The reconstructed network did not evaluate back to the input matrix
    /// (internal consistency failure).
    ReconstructionMismatch,
}

impl fmt::Display for NetworkError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyMatrix => write!(f, "matrix must contain at least one row"),
            Self::EmptyPolynomialSequence => {
                write!(
                    f,
                    "polynomial sequence must contain at least one polynomial"
                )
            }
            Self::UnsupportedStartingDegree { found } => write!(
                f,
                "expected the supplied monic sequence to start at degree 0 or 1, found degree {}",
                found
            ),
            Self::UnexpectedPolynomialDegree {
                index,
                expected,
                found,
            } => write!(
                f,
                "polynomial P_{} should have degree {}, found {:?}",
                index, expected, found
            ),
            Self::NonMonicPolynomial { index, leading } => write!(
                f,
                "polynomial P_{} is not monic; leading coefficient is {:?}",
                index, leading
            ),
            Self::InvalidRowLength {
                row,
                expected,
                found,
            } => write!(
                f,
                "row {} should have length {}, found {}",
                row, expected, found
            ),
            Self::NegativeEntry { row, col, value } => {
                write!(f, "matrix entry ({},{}) = {} is negative", row, col, value)
            }
            Self::NonUnitDiagonal { row, value } => write!(
                f,
                "diagonal entry ({},{}) should be 1, found {}",
                row, row, value
            ),
            Self::NotSquare {
                row,
                expected,
                found,
            } => write!(
                f,
                "matrix is not square: row {} has length {}, expected {}",
                row, found, expected
            ),
            Self::NonzeroAboveDiagonal { row, col, value } => write!(
                f,
                "matrix is not lower triangular: entry ({},{}) = {} lies above the diagonal",
                row, col, value
            ),
            Self::NotTotallyNonnegative { reason, witness } => {
                write!(f, "matrix is not totally nonnegative: {}", reason)?;
                if let Some(witness) = witness {
                    write!(f, "; witness: {}", witness)?;
                }
                Ok(())
            }
            Self::ZeroTailViolation { row, value } => write!(
                f,
                "first-column zero tail violated at row {}: expected 0, found {}",
                row, value
            ),
            Self::NegativeReducedEntry { row, col, value } => write!(
                f,
                "Whitney reduction produced a negative entry at ({},{}): {}",
                row, col, value
            ),
            Self::ReconstructionMismatch => {
                write!(
                    f,
                    "reconstructed network does not evaluate back to the input matrix"
                )
            }
        }
    }
}

impl std::error::Error for NetworkError {}

/// Build the lower-unitriangular coefficient matrix from a monic polynomial sequence.
///
/// The input may start at `P_0 = 1` or at `P_1`. In the latter case this
/// function prepends `P_0 = 1` automatically.
pub fn coefficient_matrix_from_monic_polynomials(
    polys: &[Polynomial<BigInt>],
) -> Result<Vec<Vec<BigInt>>, NetworkError> {
    if polys.is_empty() {
        return Err(NetworkError::EmptyPolynomialSequence);
    }

    let first_degree = polys[0]
        .degree()
        .ok_or(NetworkError::UnexpectedPolynomialDegree {
            index: 0,
            expected: 0,
            found: None,
        })?;

    let offset = match first_degree {
        0 => 0,
        1 => 1,
        d => return Err(NetworkError::UnsupportedStartingDegree { found: d }),
    };

    let mut rows = Vec::with_capacity(polys.len() + offset);
    if offset == 1 {
        rows.push(vec![BigInt::one()]);
    }

    for (supplied_index, poly) in polys.iter().enumerate() {
        let sequence_index = supplied_index + offset;
        let found_degree = poly.degree();
        if found_degree != Some(sequence_index) {
            return Err(NetworkError::UnexpectedPolynomialDegree {
                index: sequence_index,
                expected: sequence_index,
                found: found_degree,
            });
        }

        let leading = poly.leading_coefficient();
        if leading.as_ref() != Some(&BigInt::one()) {
            return Err(NetworkError::NonMonicPolynomial {
                index: sequence_index,
                leading,
            });
        }

        rows.push(poly.coeffs().to_vec());
    }

    Ok(rows)
}

/// Validate the lower-unitriangular / nonnegative shape constraints.
pub fn validate_lower_unitriangular_nonnegative_matrix(
    rows: &[Vec<BigInt>],
) -> Result<(), NetworkError> {
    if rows.is_empty() {
        return Err(NetworkError::EmptyMatrix);
    }

    for (n, row) in rows.iter().enumerate() {
        let expected_len = n + 1;
        if row.len() != expected_len {
            return Err(NetworkError::InvalidRowLength {
                row: n,
                expected: expected_len,
                found: row.len(),
            });
        }

        for (k, value) in row.iter().enumerate() {
            if value < &BigInt::zero() {
                return Err(NetworkError::NegativeEntry {
                    row: n,
                    col: k,
                    value: value.clone(),
                });
            }
        }

        if row[n] != BigInt::one() {
            return Err(NetworkError::NonUnitDiagonal {
                row: n,
                value: row[n].clone(),
            });
        }
    }

    Ok(())
}

/// Convert a square lower-triangular matrix into the triangular row format
/// used by this module (row `n` keeps its first `n + 1` entries).
///
/// Errors name the first row of the wrong length ([`NetworkError::NotSquare`])
/// or the first nonzero entry above the diagonal
/// ([`NetworkError::NonzeroAboveDiagonal`]).  Nonnegativity and the unit
/// diagonal are checked separately by
/// [`validate_lower_unitriangular_nonnegative_matrix`].
pub fn lower_triangular_rows_from_square(
    matrix: &[Vec<BigInt>],
) -> Result<Vec<Vec<BigInt>>, NetworkError> {
    if matrix.is_empty() {
        return Err(NetworkError::EmptyMatrix);
    }
    let size = matrix.len();
    for (row, entries) in matrix.iter().enumerate() {
        if entries.len() != size {
            return Err(NetworkError::NotSquare {
                row,
                expected: size,
                found: entries.len(),
            });
        }
        if let Some((col, value)) = entries
            .iter()
            .enumerate()
            .skip(row + 1)
            .find(|(_, value)| !value.is_zero())
        {
            return Err(NetworkError::NonzeroAboveDiagonal {
                row,
                col,
                value: value.clone(),
            });
        }
    }
    Ok(matrix
        .iter()
        .enumerate()
        .map(|(row, entries)| entries[..=row].to_vec())
        .collect())
}

/// Matrices with at most this many rows get an exhaustive search for a
/// negative minor when Neville elimination rejects them.
pub const EXHAUSTIVE_MINOR_SEARCH_MAX_DIMENSION: usize = 10;

/// Upper bound on the number of minors evaluated by the witness search for
/// larger matrices.
const MINOR_SEARCH_BUDGET: usize = 200_000;

/// Decide whether a lower-unitriangular matrix (triangular row format) is
/// totally nonnegative, with a structured diagnosis on failure.
///
/// The decision is made by exact Neville elimination
/// ([`check_totally_nonnegative_neville_bigint`]) on the zero-padded square matrix.  On
/// failure, [`NetworkError::NotTotallyNonnegative`] carries the elimination message and an
/// explicit [`NegativeMinor`] found by a search that first tries structured
/// families (contiguous rows against contiguous or initial columns) and then
/// all minors, exhaustively for matrices with at most
/// [`EXHAUSTIVE_MINOR_SEARCH_MAX_DIMENSION`] rows and within a fixed budget
/// otherwise.  Shape problems and negative entries are reported first, by
/// [`validate_lower_unitriangular_nonnegative_matrix`].
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::tnn_network::{check_lower_unitriangular_totally_nonnegative, NetworkError};
///
/// let rows: Vec<Vec<BigInt>> = vec![vec![1], vec![1, 1], vec![2, 0, 1]]
///     .into_iter()
///     .map(|r: Vec<i64>| r.into_iter().map(BigInt::from).collect())
///     .collect();
/// match check_lower_unitriangular_totally_nonnegative(&rows) {
///     Err(NetworkError::NotTotallyNonnegative { witness: Some(minor), .. }) => {
///         // rows {1,2}, columns {0,1}: 1*0 - 1*2 = -2
///         assert_eq!(minor.rows, vec![1, 2]);
///         assert_eq!(minor.cols, vec![0, 1]);
///         assert_eq!(minor.value, BigInt::from(-2));
///     }
///     other => panic!("unexpected {:?}", other),
/// }
/// ```
pub fn check_lower_unitriangular_totally_nonnegative(
    rows: &[Vec<BigInt>],
) -> Result<(), NetworkError> {
    validate_lower_unitriangular_nonnegative_matrix(rows)?;
    let padded = pad_to_square_matrix(rows);
    check_totally_nonnegative_neville_bigint(&padded).map_err(|error| {
        NetworkError::NotTotallyNonnegative {
            reason: error.to_string(),
            witness: find_negative_minor(&padded),
        }
    })
}

/// Deprecated name of [`check_lower_unitriangular_totally_nonnegative`].
#[deprecated(
    since = "0.6.0",
    note = "renamed to `check_lower_unitriangular_totally_nonnegative`"
)]
pub fn check_lower_unitriangular_tnn(rows: &[Vec<BigInt>]) -> Result<(), NetworkError> {
    check_lower_unitriangular_totally_nonnegative(rows)
}

/// Reconstruct the canonical network certificate from a lower-unitriangular TNN matrix.
///
/// Errors are those of [`check_lower_unitriangular_totally_nonnegative`], plus the Whitney
/// reduction diagnostics [`NetworkError::ZeroTailViolation`] and
/// [`NetworkError::NegativeReducedEntry`] (which cannot occur for a TNN input).
pub fn reconstruct_canonical_totally_nonnegative_network(
    rows: &[Vec<BigInt>],
) -> Result<CanonicalPlanarNetwork<BigRational>, NetworkError> {
    check_lower_unitriangular_totally_nonnegative(rows)?;

    let first_zero = rows
        .iter()
        .enumerate()
        .skip(1)
        .find(|(_, row)| row[0].is_zero())
        .map(|(idx, _)| idx);
    if let Some(start) = first_zero {
        for (row_idx, row) in rows.iter().enumerate().skip(start) {
            if !row[0].is_zero() {
                return Err(NetworkError::ZeroTailViolation {
                    row: row_idx,
                    value: row[0].clone(),
                });
            }
        }
    }

    let qrows = bigint_triangular_to_q(rows);
    let lambda = reconstruct_lambda_from_q_rows(&qrows)?;
    let network = CanonicalPlanarNetwork { lambda };
    if !verify_path_matrix_certificate(rows, &network) {
        return Err(NetworkError::ReconstructionMismatch);
    }

    Ok(network)
}

/// Deprecated name of [`reconstruct_canonical_totally_nonnegative_network`].
#[deprecated(
    since = "0.6.0",
    note = "renamed to `reconstruct_canonical_totally_nonnegative_network`"
)]
pub fn reconstruct_canonical_tnn_network(
    rows: &[Vec<BigInt>],
) -> Result<CanonicalPlanarNetwork<BigRational>, NetworkError> {
    reconstruct_canonical_totally_nonnegative_network(rows)
}

/// Build a full canonical TNN certificate from a monic polynomial sequence.
pub fn build_totally_nonnegative_certificate_from_monic_polynomials(
    polys: &[Polynomial<BigInt>],
) -> Result<CanonicalTotallyNonnegativeProof, NetworkError> {
    let coefficient_matrix = coefficient_matrix_from_monic_polynomials(polys)?;
    let network = reconstruct_canonical_totally_nonnegative_network(&coefficient_matrix)?;
    Ok(CanonicalTotallyNonnegativeProof {
        coefficient_matrix,
        network,
    })
}

/// Deprecated name of [`build_totally_nonnegative_certificate_from_monic_polynomials`].
#[deprecated(
    since = "0.6.0",
    note = "renamed to `build_totally_nonnegative_certificate_from_monic_polynomials`"
)]
pub fn build_tnn_certificate_from_monic_polynomials(
    polys: &[Polynomial<BigInt>],
) -> Result<CanonicalTotallyNonnegativeProof, NetworkError> {
    build_totally_nonnegative_certificate_from_monic_polynomials(polys)
}

/// Evaluate the triangular path matrix of `Gamma_N(lambda)`.
///
/// # Panics
///
/// Panics if `network` is malformed, that is, if `lambda[n]` does not have
/// exactly `n + 1` entries for every row `n`.  Networks built by
/// [`reconstruct_canonical_totally_nonnegative_network`] always have this shape;
/// [`verify_path_matrix_certificate`] checks it instead of panicking.
pub fn evaluate_path_matrix(
    network: &CanonicalPlanarNetwork<BigRational>,
) -> Vec<Vec<BigRational>> {
    let num_rows = network.num_rows();
    let mut matrix = Vec::with_capacity(num_rows);

    for source in 0..num_rows {
        let mut dp: Vec<Vec<BigRational>> = (0..=source)
            .map(|i| vec![BigRational::zero(); i + 1])
            .collect();
        dp[source][0] = BigRational::one();

        for i in (0..=source).rev() {
            for j in 0..=i {
                let current = dp[i][j].clone();
                if current.is_zero() {
                    continue;
                }

                if j < i {
                    dp[i][j + 1] += current.clone();
                    dp[i - 1][j] += current * network.lambda[i - 1][j].clone();
                }
            }
        }

        matrix.push((0..=source).map(|k| dp[k][k].clone()).collect());
    }

    matrix
}

/// Check whether a network evaluates back to a given integer matrix.
///
/// A malformed network (some `lambda[n]` without exactly `n + 1` entries)
/// evaluates to nothing and is rejected with `false`.
pub fn verify_path_matrix_certificate(
    rows: &[Vec<BigInt>],
    network: &CanonicalPlanarNetwork<BigRational>,
) -> bool {
    is_well_formed_network(network) && evaluate_path_matrix(network) == bigint_triangular_to_q(rows)
}

/// Whether `lambda[n]` has exactly `n + 1` entries for every row `n`, the shape
/// [`evaluate_path_matrix`] requires.
fn is_well_formed_network(network: &CanonicalPlanarNetwork<BigRational>) -> bool {
    network
        .lambda
        .iter()
        .enumerate()
        .all(|(n, row)| row.len() == n + 1)
}

/// Chain polynomials of a lower-unitriangular integer matrix.
///
/// Returns `p_0, ..., p_N` for the `N + 1` rows of `rows` (triangular row
/// format: row `n` has `n + 1` entries and `rows[n][n] = 1`), where `p_0 = 1`
/// and `p_n(t) = t * sum_{k<n} r_{n,k} p_k(t)`.  Coefficient vectors are in
/// ascending degree order, trimmed, with the zero polynomial as `[0]`.
///
/// Only the shape and unit diagonal are required here; signs are not
/// checked, so the generator also serves non-TNN matrices.  Use
/// [`check_chain_polynomials_bigint`] for the Brändén--Saud Maia Leite
/// properties.
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::tnn_network::chain_polynomials_bigint;
///
/// // Pascal's triangle: p_n = sum_k k! S(n,k) t^k.
/// let pascal: Vec<Vec<BigInt>> = vec![vec![1], vec![1, 1], vec![1, 2, 1], vec![1, 3, 3, 1]]
///     .into_iter()
///     .map(|r: Vec<i64>| r.into_iter().map(BigInt::from).collect())
///     .collect();
/// let p = chain_polynomials_bigint(&pascal).unwrap();
/// let expected: Vec<BigInt> = [0, 1, 6, 6].iter().map(|&c| BigInt::from(c)).collect();
/// assert_eq!(p[3], expected);
/// ```
pub fn chain_polynomials_bigint(rows: &[Vec<BigInt>]) -> Result<Vec<Vec<BigInt>>, NetworkError> {
    validate_unitriangular_shape(rows, |row, value| NetworkError::NonUnitDiagonal {
        row,
        value: value.clone(),
    })?;
    Ok(chain_polynomials_unchecked(rows))
}

/// Chain polynomials of a lower-unitriangular rational matrix.
///
/// Same conventions as [`chain_polynomials_bigint`].  A diagonal entry other
/// than `1` is reported as [`NetworkError::NonUnitDiagonal`] with the value's
/// integer part (exact rational diagonals equal to `1` are required).
pub fn chain_polynomials_rational(
    rows: &[Vec<BigRational>],
) -> Result<Vec<Vec<BigRational>>, NetworkError> {
    validate_unitriangular_shape(rows, |row, value| NetworkError::NonUnitDiagonal {
        row,
        value: value.to_integer(),
    })?;
    Ok(chain_polynomials_unchecked(rows))
}

fn validate_unitriangular_shape<C: crate::CoeffRing>(
    rows: &[Vec<C>],
    non_unit: impl Fn(usize, &C) -> NetworkError,
) -> Result<(), NetworkError> {
    if rows.is_empty() {
        return Err(NetworkError::EmptyMatrix);
    }
    for (n, row) in rows.iter().enumerate() {
        if row.len() != n + 1 {
            return Err(NetworkError::InvalidRowLength {
                row: n,
                expected: n + 1,
                found: row.len(),
            });
        }
        if row[n] != C::one() {
            return Err(non_unit(n, &row[n]));
        }
    }
    Ok(())
}

fn chain_polynomials_unchecked<C: crate::CoeffRing>(rows: &[Vec<C>]) -> Vec<Vec<C>> {
    let mut chain: Vec<Vec<C>> = Vec::with_capacity(rows.len());
    for (n, row) in rows.iter().enumerate() {
        if n == 0 {
            chain.push(vec![C::one()]);
            continue;
        }
        let mut sum = vec![C::zero()];
        for (k, entry) in row.iter().take(n).enumerate() {
            if !entry.is_zero() {
                sum = vec_poly::add(&sum, &vec_poly::scale(&chain[k], entry));
            }
        }
        chain.push(vec_poly::trim(&vec_poly::shift(&sum)));
    }
    chain
}

/// Report of [`check_chain_polynomials_bigint`].
///
/// The TNN verdict and the chain-polynomial properties are computed
/// independently, so the report also shows what happens for non-TNN
/// matrices (where the theorem says nothing).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChainPolynomialReport {
    /// `p_0, ..., p_N`, ascending coefficients, zero polynomial as `[0]`.
    pub chain_polynomials: Vec<Vec<BigInt>>,
    /// TNN verdict of the input matrix, with diagnostics on failure
    /// (see [`check_lower_unitriangular_totally_nonnegative`]).
    pub totally_nonnegative_verdict: Result<(), NetworkError>,
    /// Whether `p_n` is real-rooted (the zero polynomial counts as real-rooted).
    pub real_rooted: Vec<bool>,
    /// Whether `p_n` is real-rooted with all zeros in `[-1, 0]`.
    pub roots_in_minus_one_zero: Vec<bool>,
    /// Entry `n` says whether `p_n` interlaces `p_{n+1}` (zeros of `p_n`
    /// weakly to the left, in the directed sense of
    /// [`crate::check_weak_interlacing_bigint_coeffs`]); the zero polynomial
    /// interlaces every real-rooted polynomial in both directions.
    pub interlacing: Vec<bool>,
}

impl ChainPolynomialReport {
    /// Whether the matrix is totally nonnegative.
    pub fn is_totally_nonnegative(&self) -> bool {
        self.totally_nonnegative_verdict.is_ok()
    }

    /// Deprecated name of [`ChainPolynomialReport::is_totally_nonnegative`].
    #[deprecated(since = "0.6.0", note = "renamed to `is_totally_nonnegative`")]
    pub fn is_tnn(&self) -> bool {
        self.is_totally_nonnegative()
    }

    /// Whether every chain polynomial is real-rooted in `[-1, 0]` and
    /// consecutive ones interlace (the conclusion of the theorem).
    pub fn chain_conclusion_holds(&self) -> bool {
        self.roots_in_minus_one_zero.iter().all(|&ok| ok) && self.interlacing.iter().all(|&ok| ok)
    }

    /// `true` would contradict Brändén--Saud Maia Leite, Theorem 3.7: the
    /// matrix is TNN but the chain-polynomial conclusion fails.
    pub fn contradicts_theorem(&self) -> bool {
        self.is_totally_nonnegative() && !self.chain_conclusion_holds()
    }
}

/// Generate the chain polynomials of a lower-unitriangular integer matrix
/// and check them against Brändén--Saud Maia Leite (arXiv:2412.06595,
/// Theorem 3.7).
///
/// The report contains the TNN verdict of the matrix and, separately, exact
/// checks of real-rootedness, location of the zeros in `[-1, 0]`, and
/// interlacing `p_n ≺ p_{n+1}`.  Only shape errors (row lengths, unit
/// diagonal) are returned as `Err`; a negative entry or a failing minor
/// appears in [`ChainPolynomialReport::totally_nonnegative_verdict`].
///
/// The theorem is about the chain polynomials, **not** about the row
/// polynomials `sum_k r_{n,k} t^k`: `[[1,0,0],[1,1,0],[1,1,1]]` is TNN, its
/// chain polynomials `1, t, t + t^2` behave as predicted, but its last row
/// `1 + t + t^2` has no real zeros.
pub fn check_chain_polynomials_bigint(
    rows: &[Vec<BigInt>],
) -> Result<ChainPolynomialReport, NetworkError> {
    let chain_polynomials = chain_polynomials_bigint(rows)?;
    let totally_nonnegative_verdict = check_lower_unitriangular_totally_nonnegative(rows);
    let is_zero = |p: &[BigInt]| vec_poly::degree(p).is_none();
    let real_rooted: Vec<bool> = chain_polynomials
        .iter()
        .map(|p| is_zero(p) || is_real_rooted_bigint_coeffs(p))
        .collect();
    let minus_one_zero_test = [BigInt::zero(), BigInt::one(), BigInt::one()];
    let roots_in_minus_one_zero = chain_polynomials
        .iter()
        .zip(&real_rooted)
        .map(|(p, &rr)| {
            rr && (is_zero(p)
                || count_root_signs_prs_bigint_coeffs(p, &minus_one_zero_test).positive == 0)
        })
        .collect();
    let interlacing = chain_polynomials
        .windows(2)
        .zip(real_rooted.windows(2))
        .map(|(pair, rr)| match (is_zero(&pair[0]), is_zero(&pair[1])) {
            (true, true) => true,
            (true, false) => rr[1],
            (false, true) => rr[0],
            (false, false) => {
                check_weak_interlacing_bigint_coeffs(&pair[0], &pair[1]) == Some(true)
            }
        })
        .collect();
    Ok(ChainPolynomialReport {
        chain_polynomials,
        totally_nonnegative_verdict,
        real_rooted,
        roots_in_minus_one_zero,
        interlacing,
    })
}

fn submatrix(matrix: &[Vec<BigInt>], rows: &[usize], cols: &[usize]) -> Vec<Vec<BigInt>> {
    rows.iter()
        .map(|&r| {
            cols.iter()
                .map(|&c| matrix[r].get(c).cloned().unwrap_or_else(BigInt::zero))
                .collect()
        })
        .collect()
}

/// Search a square matrix for a negative minor, structured families first.
pub(crate) fn find_negative_minor(matrix: &[Vec<BigInt>]) -> Option<NegativeMinor> {
    let n = matrix.len();
    let mut budget = MINOR_SEARCH_BUDGET;
    let test = |rows: &[usize], cols: &[usize], budget: &mut usize| -> Option<NegativeMinor> {
        *budget = budget.saturating_sub(1);
        let value = determinant(&submatrix(matrix, rows, cols));
        (value < BigInt::zero()).then(|| NegativeMinor {
            rows: rows.to_vec(),
            cols: cols.to_vec(),
            value,
        })
    };

    // Entries (1x1 minors).
    for (r, row) in matrix.iter().enumerate() {
        if let Some((c, value)) = row.iter().enumerate().find(|(_, v)| **v < BigInt::zero()) {
            return Some(NegativeMinor {
                rows: vec![r],
                cols: vec![c],
                value: value.clone(),
            });
        }
    }

    // Contiguous rows against contiguous columns, and against
    // `{0, ..., size-2} ∪ {c}` (the minors behind Whitney/Neville steps).
    for size in 2..=n {
        for r0 in 0..=(n - size) {
            let rows: Vec<usize> = (r0..r0 + size).collect();
            for c0 in 0..=(n - size) {
                let cols: Vec<usize> = (c0..c0 + size).collect();
                if let Some(w) = test(&rows, &cols, &mut budget) {
                    return Some(w);
                }
            }
            for c in size..n {
                let mut cols: Vec<usize> = (0..size - 1).collect();
                cols.push(c);
                if let Some(w) = test(&rows, &cols, &mut budget) {
                    return Some(w);
                }
            }
        }
    }

    // All minors, smallest first, within the budget (exhaustive for small n).
    let exhaustive = n <= EXHAUSTIVE_MINOR_SEARCH_MAX_DIMENSION;
    for size in 2..=n {
        let row_sets = combinations(n, size);
        for rows in &row_sets {
            for cols in &row_sets {
                if !exhaustive && budget == 0 {
                    return None;
                }
                if let Some(w) = test(rows, cols, &mut budget) {
                    return Some(w);
                }
            }
        }
    }
    None
}

fn combinations(n: usize, k: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut current: Vec<usize> = (0..k).collect();
    if k > n {
        return out;
    }
    loop {
        out.push(current.clone());
        let Some(i) = (0..k).rev().find(|&i| current[i] < n - k + i) else {
            return out;
        };
        current[i] += 1;
        for j in i + 1..k {
            current[j] = current[j - 1] + 1;
        }
    }
}

pub(crate) fn reconstruct_lambda_from_q_rows(
    rows: &[Vec<BigRational>],
) -> Result<Vec<Vec<BigRational>>, NetworkError> {
    let num_rows = rows.len();
    if num_rows <= 1 {
        return Ok(vec![]);
    }

    let first_zero = rows
        .iter()
        .enumerate()
        .skip(1)
        .find(|(_, row)| row[0].is_zero())
        .map(|(idx, _)| idx)
        .unwrap_or(num_rows);

    let mut mu = Vec::with_capacity(num_rows - 1);
    for n in 0..(num_rows - 1) {
        if n < first_zero {
            mu.push(rows[n + 1][0].clone() / rows[n][0].clone());
        } else {
            mu.push(BigRational::zero());
        }
    }

    let mut reduced = Vec::with_capacity(num_rows - 1);
    for n in 0..(num_rows - 1) {
        let mut row = Vec::with_capacity(n + 1);
        for k in 0..=n {
            let entry = rows[n + 1][k + 1].clone() - mu[n].clone() * qcoeff(&rows[n], k + 1);
            if entry < BigRational::zero() {
                return Err(NetworkError::NegativeReducedEntry {
                    row: n,
                    col: k,
                    value: entry,
                });
            }
            row.push(entry);
        }
        reduced.push(row);
    }

    let gamma = reconstruct_lambda_from_q_rows(&reduced)?;
    let mut lambda = Vec::with_capacity(num_rows - 1);
    for n in 0..(num_rows - 1) {
        let mut row = Vec::with_capacity(n + 1);
        row.push(mu[n].clone());
        if n > 0 {
            row.extend(gamma[n - 1].iter().cloned());
        }
        lambda.push(row);
    }

    Ok(lambda)
}

fn pad_to_square_matrix(rows: &[Vec<BigInt>]) -> Vec<Vec<BigInt>> {
    let n = rows.len();
    let mut padded = vec![vec![BigInt::zero(); n]; n];
    for (i, row) in rows.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            padded[i][j] = value.clone();
        }
    }
    padded
}

fn bigint_triangular_to_q(rows: &[Vec<BigInt>]) -> Vec<Vec<BigRational>> {
    rows.iter()
        .map(|row| {
            row.iter()
                .map(|value| BigRational::from_integer(value.clone()))
                .collect()
        })
        .collect()
}

fn qcoeff(row: &[BigRational], index: usize) -> BigRational {
    row.get(index).cloned().unwrap_or_else(BigRational::zero)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bi(v: i64) -> BigInt {
        BigInt::from(v)
    }

    fn br(v: i64) -> BigRational {
        BigRational::from_integer(bi(v))
    }

    fn bi_poly(coeffs: &[i64]) -> Polynomial<BigInt> {
        Polynomial::new(coeffs.iter().map(|&v| bi(v)).collect())
    }

    fn bi_rows(rows: &[&[i64]]) -> Vec<Vec<BigInt>> {
        rows.iter()
            .map(|row| row.iter().map(|&v| bi(v)).collect())
            .collect()
    }

    #[test]
    fn test_coefficient_matrix_from_monic_polynomials_starting_at_zero() {
        let polys = vec![bi_poly(&[1]), bi_poly(&[1, 1]), bi_poly(&[1, 3, 1])];
        let rows = coefficient_matrix_from_monic_polynomials(&polys).unwrap();
        assert_eq!(rows, bi_rows(&[&[1], &[1, 1], &[1, 3, 1]]));
    }

    #[test]
    fn test_coefficient_matrix_from_monic_polynomials_prepends_p0() {
        let polys = vec![
            bi_poly(&[1, 1]),
            bi_poly(&[1, 2, 1]),
            bi_poly(&[1, 3, 3, 1]),
        ];
        let rows = coefficient_matrix_from_monic_polynomials(&polys).unwrap();
        assert_eq!(rows, bi_rows(&[&[1], &[1, 1], &[1, 2, 1], &[1, 3, 3, 1]]));
    }

    #[test]
    fn test_evaluate_pascal_network() {
        let network = CanonicalPlanarNetwork {
            lambda: vec![vec![br(1)], vec![br(1), br(1)], vec![br(1), br(1), br(1)]],
        };
        let matrix = evaluate_path_matrix(&network);
        assert_eq!(
            matrix,
            vec![
                vec![br(1)],
                vec![br(1), br(1)],
                vec![br(1), br(2), br(1)],
                vec![br(1), br(3), br(3), br(1)],
            ]
        );
    }

    #[test]
    fn test_vertical_weight_and_certificate_path_matrix() {
        // Pascal's triangle (OEIS A007318) is the path matrix of the all-ones network.
        let polys = vec![
            bi_poly(&[1, 1]),
            bi_poly(&[1, 2, 1]),
            bi_poly(&[1, 3, 3, 1]),
            bi_poly(&[1, 4, 6, 4, 1]),
        ];
        let proof = build_totally_nonnegative_certificate_from_monic_polynomials(&polys).unwrap();
        assert_eq!(proof.network.num_rows(), 5);
        // `lambda[n][k]` is defined exactly for 0 <= k <= n < num_rows - 1.
        for n in 0..6 {
            for k in 0..6 {
                let weight = proof.network.vertical_weight(n, k);
                if k <= n && n < 4 {
                    assert_eq!(weight, Some(&br(1)), "weight ({n}, {k})");
                } else {
                    assert_eq!(weight, None, "weight ({n}, {k})");
                }
            }
        }
        // The stored network evaluates back to the stored integer matrix.
        let evaluated = proof.evaluated_path_matrix().unwrap();
        assert_eq!(evaluated.len(), proof.coefficient_matrix.len());
        for (row, expected) in evaluated.iter().zip(&proof.coefficient_matrix) {
            let integers: Vec<BigInt> = row.iter().map(|q| q.to_integer()).collect();
            assert!(row.iter().all(|q| q.is_integer()));
            assert_eq!(&integers, expected);
        }
        assert!(proof.verify());
        // Tampering with a weight is detected by `verify` and changes the path matrix.
        let mut tampered = proof.clone();
        tampered.network.lambda[3][1] = br(2);
        assert!(!tampered.verify());
        assert_ne!(tampered.evaluated_path_matrix().unwrap(), evaluated);
    }

    #[test]
    fn evaluated_path_matrix_of_a_malformed_network_is_none() {
        let polys = vec![bi_poly(&[1]), bi_poly(&[1, 1]), bi_poly(&[1, 2, 1])];
        let proof = build_totally_nonnegative_certificate_from_monic_polynomials(&polys).unwrap();
        assert!(proof.evaluated_path_matrix().is_some());
        for malformed in [
            vec![vec![]],
            vec![vec![br(1), br(1)]],
            vec![vec![br(1)], vec![br(1)]],
        ] {
            let mut tampered = proof.clone();
            tampered.network.lambda = malformed;
            assert_eq!(tampered.evaluated_path_matrix(), None);
            assert!(!tampered.verify());
        }
    }

    #[test]
    fn test_reconstruct_pascal_network() {
        let rows = bi_rows(&[&[1], &[1, 1], &[1, 2, 1], &[1, 3, 3, 1]]);
        let network = reconstruct_canonical_totally_nonnegative_network(&rows).unwrap();
        assert_eq!(
            network.lambda,
            vec![vec![br(1)], vec![br(1), br(1)], vec![br(1), br(1), br(1)],]
        );
        assert!(verify_path_matrix_certificate(&rows, &network));
    }

    #[test]
    fn test_reconstruct_zero_tail_network() {
        let network = CanonicalPlanarNetwork {
            lambda: vec![vec![br(2)], vec![br(0), br(3)]],
        };
        let rows_q = evaluate_path_matrix(&network);
        let rows: Vec<Vec<BigInt>> = rows_q
            .iter()
            .map(|row| row.iter().map(|q| q.to_integer()).collect())
            .collect();
        assert_eq!(rows, bi_rows(&[&[1], &[2, 1], &[0, 3, 1]]));

        let reconstructed = reconstruct_canonical_totally_nonnegative_network(&rows).unwrap();
        assert_eq!(reconstructed, network);
    }

    #[test]
    fn test_reconstruct_rejects_non_tnn_matrix() {
        let rows = bi_rows(&[&[1], &[1, 1], &[2, 0, 1]]);
        match reconstruct_canonical_totally_nonnegative_network(&rows) {
            Err(NetworkError::NotTotallyNonnegative {
                witness: Some(witness),
                ..
            }) => {
                assert!(witness.verify(&rows));
                assert_eq!(witness.rows, vec![1, 2]);
                assert_eq!(witness.cols, vec![0, 1]);
                assert_eq!(witness.value, bi(-2));
            }
            other => panic!(
                "expected NotTotallyNonnegative with witness, got {:?}",
                other
            ),
        }
    }

    #[test]
    fn test_not_tnn_witness_needs_a_3x3_minor() {
        // The lower-triangular Toeplitz matrix of 1 + t + t^2: entries and
        // 2x2 minors are nonnegative (the sequence is log-concave), but the
        // sequence is not PF, so a 3x3 minor is negative.
        let rows = bi_rows(&[&[1], &[1, 1], &[1, 1, 1], &[0, 1, 1, 1]]);
        let error = check_lower_unitriangular_totally_nonnegative(&rows).unwrap_err();
        match error {
            NetworkError::NotTotallyNonnegative {
                witness: Some(witness),
                ..
            } => {
                assert!(witness.verify(&rows));
                assert_eq!(witness.rows, vec![1, 2, 3]);
                assert_eq!(witness.cols, vec![0, 1, 2]);
                assert_eq!(witness.value, bi(-1));
            }
            other => panic!(
                "expected NotTotallyNonnegative with witness, got {:?}",
                other
            ),
        }
    }

    #[test]
    fn test_witness_search_agrees_with_neville_on_small_matrices() {
        // Every lower-unitriangular 3x3 matrix with entries in 0..=2 below
        // the diagonal: Neville rejects exactly when a witness exists.
        for a in 0..=2 {
            for b in 0..=2 {
                for c in 0..=2 {
                    let rows = bi_rows(&[&[1], &[a, 1], &[b, c, 1]]);
                    match check_lower_unitriangular_totally_nonnegative(&rows) {
                        Ok(()) => {
                            let padded = pad_to_square_matrix(&rows);
                            assert_eq!(find_negative_minor(&padded), None, "{:?}", rows);
                        }
                        Err(NetworkError::NotTotallyNonnegative {
                            witness: Some(w), ..
                        }) => assert!(w.verify(&rows)),
                        other => panic!("unexpected {:?} for {:?}", other, rows),
                    }
                }
            }
        }
    }

    #[test]
    fn test_shape_diagnostics() {
        let square = bi_rows(&[&[1, 0], &[1, 1]]);
        assert_eq!(
            lower_triangular_rows_from_square(&square).unwrap(),
            bi_rows(&[&[1], &[1, 1]])
        );
        assert_eq!(
            lower_triangular_rows_from_square(&bi_rows(&[&[1, 2], &[1, 1]])),
            Err(NetworkError::NonzeroAboveDiagonal {
                row: 0,
                col: 1,
                value: bi(2)
            })
        );
        assert_eq!(
            lower_triangular_rows_from_square(&bi_rows(&[&[1, 0], &[1]])),
            Err(NetworkError::NotSquare {
                row: 1,
                expected: 2,
                found: 1
            })
        );
        assert_eq!(
            check_lower_unitriangular_totally_nonnegative(&bi_rows(&[&[1], &[-1, 1]])),
            Err(NetworkError::NegativeEntry {
                row: 1,
                col: 0,
                value: bi(-1)
            })
        );
        assert_eq!(
            check_lower_unitriangular_totally_nonnegative(&bi_rows(&[&[1], &[1, 2]])),
            Err(NetworkError::NonUnitDiagonal {
                row: 1,
                value: bi(2)
            })
        );
    }

    #[test]
    fn test_chain_polynomials_pascal() {
        // p_n = sum_k k! S(n,k) t^k (ordered set partitions / Fubini).
        let rows = bi_rows(&[&[1], &[1, 1], &[1, 2, 1], &[1, 3, 3, 1], &[1, 4, 6, 4, 1]]);
        let chain = chain_polynomials_bigint(&rows).unwrap();
        let expected: Vec<Vec<BigInt>> = [
            vec![1],
            vec![0, 1],
            vec![0, 1, 2],
            vec![0, 1, 6, 6],
            vec![0, 1, 14, 36, 24],
        ]
        .iter()
        .map(|r| r.iter().map(|&v| bi(v)).collect())
        .collect();
        assert_eq!(chain, expected);
        let report = check_chain_polynomials_bigint(&rows).unwrap();
        assert!(report.is_totally_nonnegative());
        assert!(report.chain_conclusion_holds());
        assert!(!report.contradicts_theorem());
    }

    #[test]
    fn test_chain_polynomials_rational_matches_integer() {
        let rows = bi_rows(&[&[1], &[2, 1], &[1, 3, 1]]);
        let qrows = bigint_triangular_to_q(&rows);
        let chain_q = chain_polynomials_rational(&qrows).unwrap();
        let chain_z = chain_polynomials_bigint(&rows).unwrap();
        let lifted: Vec<Vec<BigRational>> = chain_z
            .iter()
            .map(|p| {
                p.iter()
                    .map(|c| BigRational::from_integer(c.clone()))
                    .collect()
            })
            .collect();
        assert_eq!(chain_q, lifted);
    }

    #[test]
    fn test_chain_polynomials_identity_are_zero() {
        let rows = bi_rows(&[&[1], &[0, 1], &[0, 0, 1]]);
        let report = check_chain_polynomials_bigint(&rows).unwrap();
        assert_eq!(report.chain_polynomials[1], vec![bi(0)]);
        assert!(report.is_totally_nonnegative());
        assert!(report.chain_conclusion_holds());
    }

    #[test]
    fn test_build_certificate_from_polynomials() {
        let polys = vec![
            bi_poly(&[1, 1]),
            bi_poly(&[1, 2, 1]),
            bi_poly(&[1, 3, 3, 1]),
        ];
        let proof = build_totally_nonnegative_certificate_from_monic_polynomials(&polys).unwrap();
        assert!(proof.verify());
        assert_eq!(
            proof.coefficient_matrix,
            bi_rows(&[&[1], &[1, 1], &[1, 2, 1], &[1, 3, 3, 1]])
        );
    }

    #[test]
    fn malformed_networks_fail_verification_instead_of_panicking() {
        let rows = bi_rows(&[&[1], &[1, 1]]);
        let ragged = CanonicalPlanarNetwork {
            lambda: vec![vec![BigRational::one(), BigRational::one()]],
        };
        assert!(!verify_path_matrix_certificate(&rows, &ragged));
        let wrong_row = CanonicalPlanarNetwork {
            lambda: vec![vec![BigRational::one()], vec![BigRational::one()]],
        };
        assert!(!verify_path_matrix_certificate(&rows, &wrong_row));
        let good = reconstruct_canonical_totally_nonnegative_network(&rows).unwrap();
        assert!(verify_path_matrix_certificate(&rows, &good));
    }
}
