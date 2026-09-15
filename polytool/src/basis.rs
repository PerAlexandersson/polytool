//! Exact coordinate computations in polynomial bases.
//!
//! This module provides:
//!
//! - a general exact routine for expressing a polynomial in a prescribed basis of
//!   `Pol_{\le n}`,
//! - the standard Bernstein basis
//!   `{ binomial(n,j) t^j (1-t)^(n-j) }_{j=0}^n`,
//! - the "magic basis"
//!   `{ t^i (1+t)^{d-i} }_{i=0}^d`,
//! - convenience helpers for checking magic positivity.

use crate::polynomial::{CoeffRing, FieldRing};
use crate::Polynomial;
use num_bigint::BigInt;
use num_rational::Ratio;
use std::fmt;

/// Exact rational type used for basis coordinates from integer input.
pub type BigRational = Ratio<BigInt>;

/// Errors produced while expanding a polynomial in a basis.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BasisError {
    EmptyBasis,
    TargetDegreeTooLarge {
        degree: usize,
        basis_degree_bound: usize,
    },
    BasisPolynomialTooLarge {
        index: usize,
        degree: usize,
        basis_degree_bound: usize,
    },
    SingularBasis,
}

impl fmt::Display for BasisError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyBasis => write!(f, "basis must contain at least one polynomial"),
            Self::TargetDegreeTooLarge {
                degree,
                basis_degree_bound,
            } => write!(
                f,
                "target polynomial has degree {}, but the basis is for degree at most {}",
                degree, basis_degree_bound
            ),
            Self::BasisPolynomialTooLarge {
                index,
                degree,
                basis_degree_bound,
            } => write!(
                f,
                "basis polynomial {} has degree {}, exceeding the degree bound {}",
                index, degree, basis_degree_bound
            ),
            Self::SingularBasis => write!(
                f,
                "the supplied basis is singular and does not form a basis of the ambient space"
            ),
        }
    }
}

impl std::error::Error for BasisError {}

/// Express `target` in the supplied basis of `Pol_{\le n}`.
///
/// The basis must have length `n + 1`, where `n = basis.len() - 1`, and each
/// basis polynomial must have degree at most `n`. Coefficients are computed
/// exactly over the same field as the polynomial coefficients.
pub fn coordinates_in_basis<C: FieldRing>(
    target: &Polynomial<C>,
    basis: &[Polynomial<C>],
) -> Result<Vec<C>, BasisError> {
    if basis.is_empty() {
        return Err(BasisError::EmptyBasis);
    }

    let degree_bound = basis.len() - 1;
    if let Some(degree) = target.degree() {
        if degree > degree_bound {
            return Err(BasisError::TargetDegreeTooLarge {
                degree,
                basis_degree_bound: degree_bound,
            });
        }
    }

    for (index, basis_poly) in basis.iter().enumerate() {
        if let Some(degree) = basis_poly.degree() {
            if degree > degree_bound {
                return Err(BasisError::BasisPolynomialTooLarge {
                    index,
                    degree,
                    basis_degree_bound: degree_bound,
                });
            }
        }
    }

    let matrix: Vec<Vec<C>> = (0..=degree_bound)
        .map(|row| basis.iter().map(|poly| poly.coeff(row)).collect())
        .collect();
    let rhs: Vec<C> = (0..=degree_bound).map(|row| target.coeff(row)).collect();

    solve_square_linear_system(&matrix, &rhs).ok_or(BasisError::SingularBasis)
}

/// Express an `i64`-coefficient polynomial in an integral basis exactly over `Q`.
pub fn coordinates_in_basis_i64(
    target: &Polynomial<i64>,
    basis: &[Polynomial<i64>],
) -> Result<Vec<BigRational>, BasisError> {
    let target_q = polynomial_i64_to_q(target);
    let basis_q: Vec<_> = basis.iter().map(polynomial_i64_to_q).collect();
    coordinates_in_basis(&target_q, &basis_q)
}

/// Express a `BigInt`-coefficient polynomial in an integral basis exactly over `Q`.
pub fn coordinates_in_basis_bigint(
    target: &Polynomial<BigInt>,
    basis: &[Polynomial<BigInt>],
) -> Result<Vec<BigRational>, BasisError> {
    let target_q = polynomial_bigint_to_q(target);
    let basis_q: Vec<_> = basis.iter().map(polynomial_bigint_to_q).collect();
    coordinates_in_basis(&target_q, &basis_q)
}

/// Build the standard Bernstein basis of degree `degree`:
///
/// ```text
/// binomial(degree,j) t^j (1-t)^(degree-j),  0 <= j <= degree.
/// ```
pub fn bernstein_basis_bigint(degree: usize) -> Vec<Polynomial<BigInt>> {
    let t = Polynomial::variable();
    let one_minus_t = Polynomial::from_i64_coeffs(&[1, -1]);
    (0..=degree)
        .map(|j| {
            let polynomial = poly_pow(&t, j) * poly_pow(&one_minus_t, degree - j);
            polynomial.scale(&binomial_bigint(degree, j))
        })
        .collect()
}

/// Express a rational polynomial in the standard Bernstein basis of a fixed
/// ambient degree.
///
/// If `target = sum_i a_i t^i`, the returned coordinates are
///
/// ```text
/// b_j = sum_{i=0}^j a_i binomial(j,i) / binomial(degree,i).
/// ```
///
/// The ambient degree may exceed the degree of the target; this is the usual
/// Bernstein degree-elevation operation.
pub fn bernstein_basis_coordinates(
    target: &Polynomial<BigRational>,
    degree: usize,
) -> Result<Vec<BigRational>, BasisError> {
    if let Some(target_degree) = target.degree() {
        if target_degree > degree {
            return Err(BasisError::TargetDegreeTooLarge {
                degree: target_degree,
                basis_degree_bound: degree,
            });
        }
    }

    Ok((0..=degree)
        .map(|j| {
            (0..=j.min(target.coeffs().len().saturating_sub(1)))
                .map(|i| {
                    let numerator = binomial_bigint(j, i);
                    let denominator = binomial_bigint(degree, i);
                    target.coeff(i) * BigRational::new(numerator, denominator)
                })
                .fold(BigRational::from_integer(BigInt::from(0)), |sum, term| {
                    sum + term
                })
        })
        .collect())
}

/// Express an integer polynomial in the standard Bernstein basis exactly.
pub fn bernstein_basis_coordinates_bigint(
    coeffs: &[BigInt],
    degree: usize,
) -> Result<Vec<BigRational>, BasisError> {
    bernstein_basis_coordinates(
        &polynomial_bigint_to_q(&Polynomial::new(coeffs.to_vec())),
        degree,
    )
}

/// Express an `i64` polynomial in the standard Bernstein basis exactly.
pub fn bernstein_basis_coordinates_i64(
    coeffs: &[i64],
    degree: usize,
) -> Result<Vec<BigRational>, BasisError> {
    bernstein_basis_coordinates(
        &polynomial_i64_to_q(&Polynomial::from_i64_coeffs(coeffs)),
        degree,
    )
}

/// Build the magic basis of `Pol_{\le degree}`:
///
/// ```text
/// (1+t)^degree, t(1+t)^(degree-1), ..., t^degree.
/// ```
pub fn magic_basis<C: CoeffRing>(degree: usize) -> Vec<Polynomial<C>> {
    let one_plus_t = Polynomial::from_i64_coeffs(&[1, 1]);
    (0..=degree)
        .map(|i| poly_pow(&Polynomial::variable(), i) * poly_pow(&one_plus_t, degree - i))
        .collect()
}

/// Exact coordinates in the magic basis of `Pol_{\le degree}` for an `i64`-coefficient polynomial.
pub fn magic_basis_coordinates_i64(
    coeffs: &[i64],
    degree: usize,
) -> Result<Vec<BigRational>, BasisError> {
    let target = Polynomial::<i64>::from_i64_coeffs(coeffs);
    let basis = magic_basis::<i64>(degree);
    coordinates_in_basis_i64(&target, &basis)
}

/// Exact coordinates in the magic basis of `Pol_{\le degree}` for a `BigInt`-coefficient polynomial.
pub fn magic_basis_coordinates_bigint(
    coeffs: &[BigInt],
    degree: usize,
) -> Result<Vec<BigRational>, BasisError> {
    let target = Polynomial::<BigInt>::new(coeffs.to_vec());
    let basis = magic_basis::<BigInt>(degree);
    coordinates_in_basis_bigint(&target, &basis)
}

/// Check whether an `i64`-coefficient polynomial is magic positive with respect
/// to the degree bound `degree`.
pub fn is_magic_positive_i64(coeffs: &[i64], degree: usize) -> Result<bool, BasisError> {
    let zero = BigRational::from_integer(BigInt::from(0));
    Ok(magic_basis_coordinates_i64(coeffs, degree)?
        .iter()
        .all(|c| c >= &zero))
}

/// Check whether a `BigInt`-coefficient polynomial is magic positive with
/// respect to the degree bound `degree`.
pub fn is_magic_positive_bigint(coeffs: &[BigInt], degree: usize) -> Result<bool, BasisError> {
    let zero = BigRational::from_integer(BigInt::from(0));
    Ok(magic_basis_coordinates_bigint(coeffs, degree)?
        .iter()
        .all(|c| c >= &zero))
}

/// Theorem 2.13-style data attached to a magic-basis expansion.
///
/// For coordinates `c_0, ..., c_d` in the basis `{x^j (1+x)^(d-j)}`, this stores:
///
/// - the coordinates themselves,
/// - the partial sums `c_0 + ... + c_j`,
/// - the reversed partial sums `c_d + ... + c_{d-j}`,
/// - whether all coordinates are nonnegative,
/// - whether the inequalities
///   `c_0 + ... + c_j <= c_d + ... + c_{d-j}` hold for all `0 <= j <= floor(d/2)`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MagicBasisAnalysis<R> {
    pub coordinates: Vec<R>,
    pub left_partial_sums: Vec<R>,
    pub right_partial_sums: Vec<R>,
    pub all_nonnegative: bool,
    pub left_leq_right: bool,
}

/// Analyze the magic-basis coordinates of an `i64`-coefficient polynomial.
///
/// This is convenient when using the Brandén-Solus criterion from Theorem 2.13:
/// if the input is the auxiliary polynomial
/// `i(x) = sum_j c_j x^j (1+x)^(d-j)`, then `all_nonnegative && left_leq_right`
/// is exactly the condition checked there.
pub fn analyze_magic_basis_i64(
    coeffs: &[i64],
    degree: usize,
) -> Result<MagicBasisAnalysis<BigRational>, BasisError> {
    let coordinates = magic_basis_coordinates_i64(coeffs, degree)?;
    Ok(analyze_magic_coordinates(coordinates))
}

/// Analyze the magic-basis coordinates of a `BigInt`-coefficient polynomial.
pub fn analyze_magic_basis_bigint(
    coeffs: &[BigInt],
    degree: usize,
) -> Result<MagicBasisAnalysis<BigRational>, BasisError> {
    let coordinates = magic_basis_coordinates_bigint(coeffs, degree)?;
    Ok(analyze_magic_coordinates(coordinates))
}

fn analyze_magic_coordinates(coordinates: Vec<BigRational>) -> MagicBasisAnalysis<BigRational> {
    let degree = coordinates.len().saturating_sub(1);
    let zero = BigRational::from_integer(BigInt::from(0));
    let all_nonnegative = coordinates.iter().all(|c| c >= &zero);

    let mut left_running = zero.clone();
    let mut right_running = zero.clone();
    let mut left_partial_sums = Vec::new();
    let mut right_partial_sums = Vec::new();
    let mut left_leq_right = true;

    for j in 0..=degree / 2 {
        left_running += coordinates[j].clone();
        right_running += coordinates[degree - j].clone();
        left_partial_sums.push(left_running.clone());
        right_partial_sums.push(right_running.clone());
        if left_running > right_running {
            left_leq_right = false;
        }
    }

    MagicBasisAnalysis {
        coordinates,
        left_partial_sums,
        right_partial_sums,
        all_nonnegative,
        left_leq_right,
    }
}

fn polynomial_i64_to_q(poly: &Polynomial<i64>) -> Polynomial<BigRational> {
    Polynomial::new(
        poly.coeffs()
            .iter()
            .map(|&c| BigRational::from_integer(BigInt::from(c)))
            .collect(),
    )
}

fn polynomial_bigint_to_q(poly: &Polynomial<BigInt>) -> Polynomial<BigRational> {
    Polynomial::new(
        poly.coeffs()
            .iter()
            .map(|c| BigRational::from_integer(c.clone()))
            .collect(),
    )
}

fn poly_pow<C: CoeffRing>(base: &Polynomial<C>, exp: usize) -> Polynomial<C> {
    if exp == 0 {
        return Polynomial::one();
    }
    let mut result = Polynomial::one();
    let mut power = base.clone();
    let mut e = exp;
    while e > 0 {
        if e & 1 == 1 {
            result = result * power.clone();
        }
        e >>= 1;
        if e > 0 {
            power = power.clone() * power;
        }
    }
    result
}

fn binomial_bigint(n: usize, k: usize) -> BigInt {
    if k > n {
        return BigInt::from(0);
    }
    let k = k.min(n - k);
    let mut value = BigInt::from(1);
    for j in 0..k {
        value = value * BigInt::from(n - j) / BigInt::from(j + 1);
    }
    value
}

fn solve_square_linear_system<C: FieldRing>(a: &[Vec<C>], b: &[C]) -> Option<Vec<C>> {
    let n = a.len();
    if n == 0 || b.len() != n || a.iter().any(|row| row.len() != n) {
        return None;
    }

    let mut aug: Vec<Vec<C>> = a
        .iter()
        .zip(b.iter())
        .map(|(row, rhs)| {
            let mut r = row.clone();
            r.push(rhs.clone());
            r
        })
        .collect();

    for col in 0..n {
        let pivot = (col..n).find(|&row| !aug[row][col].is_zero())?;
        if pivot != col {
            aug.swap(pivot, col);
        }

        let pivot_value = aug[col][col].clone();
        for entry in aug[col].iter_mut().take(n + 1).skip(col) {
            *entry = entry.clone().field_div(pivot_value.clone());
        }
        let pivot_row = aug[col].clone();

        for (row, row_entries) in aug.iter_mut().enumerate().take(n) {
            if row == col || row_entries[col].is_zero() {
                continue;
            }
            let factor = row_entries[col].clone();
            for (j, pivot_entry) in pivot_row.iter().enumerate().take(n + 1).skip(col) {
                row_entries[j] = row_entries[j].clone() - factor.clone() * pivot_entry.clone();
            }
        }
    }

    Some((0..n).map(|i| aug[i][n].clone()).collect())
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

    #[test]
    fn test_coordinates_in_monomial_basis() {
        let target = Polynomial::<BigRational>::new(vec![br(2), br(3), br(5)]);
        let basis = vec![
            Polynomial::<BigRational>::from_i64_coeffs(&[1]),
            Polynomial::<BigRational>::from_i64_coeffs(&[0, 1]),
            Polynomial::<BigRational>::from_i64_coeffs(&[0, 0, 1]),
        ];
        let coords = coordinates_in_basis(&target, &basis).unwrap();
        assert_eq!(coords, vec![br(2), br(3), br(5)]);
    }

    #[test]
    fn test_magic_basis_degree_two() {
        let basis = magic_basis::<i64>(2);
        assert_eq!(basis[0], Polynomial::from_i64_coeffs(&[1, 2, 1]));
        assert_eq!(basis[1], Polynomial::from_i64_coeffs(&[0, 1, 1]));
        assert_eq!(basis[2], Polynomial::from_i64_coeffs(&[0, 0, 1]));
    }

    #[test]
    fn test_bernstein_basis_degree_two() {
        let basis = bernstein_basis_bigint(2);
        assert_eq!(basis[0], Polynomial::from_i64_coeffs(&[1, -2, 1]));
        assert_eq!(basis[1], Polynomial::from_i64_coeffs(&[0, 2, -2]));
        assert_eq!(basis[2], Polynomial::from_i64_coeffs(&[0, 0, 1]));
    }

    #[test]
    fn test_bernstein_coordinates_and_degree_elevation() {
        assert_eq!(
            bernstein_basis_coordinates_i64(&[1, 2, 3], 2).unwrap(),
            vec![br(1), br(2), br(6)]
        );
        assert_eq!(
            bernstein_basis_coordinates_i64(&[1, 2, 3], 3).unwrap(),
            vec![
                br(1),
                BigRational::new(bi(5), bi(3)),
                BigRational::new(bi(10), bi(3)),
                br(6),
            ]
        );
        assert_eq!(
            bernstein_basis_coordinates_i64(&[0, 0, 1], 4).unwrap(),
            vec![
                br(0),
                br(0),
                BigRational::new(bi(1), bi(6)),
                BigRational::new(bi(1), bi(2)),
                br(1),
            ]
        );
    }

    #[test]
    fn test_bernstein_coordinates_reject_small_ambient_degree() {
        assert_eq!(
            bernstein_basis_coordinates_i64(&[1, 2, 3], 1),
            Err(BasisError::TargetDegreeTooLarge {
                degree: 2,
                basis_degree_bound: 1,
            })
        );
    }

    #[test]
    fn test_magic_basis_coordinates_for_one_plus_two_t_cubed() {
        // (1 + 2t)^3 = ((1+t) + t)^3
        let coords = magic_basis_coordinates_i64(&[1, 6, 12, 8], 3).unwrap();
        assert_eq!(coords, vec![br(1), br(3), br(3), br(1)]);
        assert!(is_magic_positive_i64(&[1, 6, 12, 8], 3).unwrap());
    }

    #[test]
    fn test_magic_basis_coordinates_nonpositive_example() {
        let coords = magic_basis_coordinates_i64(&[1, 3, 1], 2).unwrap();
        assert_eq!(
            coords,
            vec![br(1), br(1), BigRational::from_integer(bi(-1))]
        );
        assert!(!is_magic_positive_i64(&[1, 3, 1], 2).unwrap());
    }

    #[test]
    fn test_magic_basis_analysis_tracks_partial_sums() {
        let analysis = analyze_magic_basis_i64(&[1, 6, 12, 8], 3).unwrap();
        assert_eq!(analysis.coordinates, vec![br(1), br(3), br(3), br(1)]);
        assert_eq!(analysis.left_partial_sums, vec![br(1), br(4)]);
        assert_eq!(analysis.right_partial_sums, vec![br(1), br(4)]);
        assert!(analysis.all_nonnegative);
        assert!(analysis.left_leq_right);
    }

    #[test]
    fn test_magic_basis_analysis_detects_failed_partial_sum_inequality() {
        let basis = magic_basis::<i64>(4);
        let poly = basis[0].clone().scale(&1) + basis[1].clone().scale(&10) + basis[4].clone();
        let analysis = analyze_magic_basis_i64(poly.coeffs(), 4).unwrap();
        assert_eq!(
            analysis.coordinates,
            vec![br(1), br(10), br(0), br(0), br(1)]
        );
        assert_eq!(analysis.left_partial_sums, vec![br(1), br(11), br(11)]);
        assert_eq!(analysis.right_partial_sums, vec![br(1), br(1), br(1)]);
        assert!(analysis.all_nonnegative);
        assert!(!analysis.left_leq_right);
    }

    #[test]
    fn test_singular_basis_rejected() {
        let target = Polynomial::<BigRational>::from_i64_coeffs(&[1, 2]);
        let basis = vec![
            Polynomial::<BigRational>::from_i64_coeffs(&[1, 0]),
            Polynomial::<BigRational>::from_i64_coeffs(&[2, 0]),
        ];
        match coordinates_in_basis(&target, &basis) {
            Err(BasisError::SingularBasis) => {}
            other => panic!("expected SingularBasis, got {:?}", other),
        }
    }
}
