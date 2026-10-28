//! Exact factorization over `Q` for univariate and monic bivariate polynomials.
//!
//! The univariate entry points use squarefree decomposition followed by the
//! Zassenhaus modular factorization algorithm.  The bivariate entry points use
//! a squarefree specialization and Hensel lifting in the translated variable.
//! Coefficients are always exact, and coefficient vectors use ascending degree
//! order.  Related arithmetic lives in [`crate::polynomial`] and
//! [`crate::vec_poly`].  The implementation follows the presentation of
//! Zassenhaus factorization in von zur Gathen and Gerhard, *Modern Computer
//! Algebra*.
//!
//! The public functions are [`irreducible_factorization_bigint_coeffs`],
//! [`irreducible_factorization_rational_coeffs`], and
//! [`irreducible_factorization_monic_bivariate`].  The `_with_budget` variants
//! provide a deterministic bound for the subset trials used during
//! recombination.

mod arithmetic;
mod bivariate;
mod modular;
mod univariate;

use std::fmt;

use num_bigint::BigInt;

use crate::BigRational;

/// Errors of the factorization routines.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PolynomialFactorizationError {
    /// The input is the zero polynomial.
    ZeroPolynomial,
    /// A bivariate input is not monic in its first variable (the coefficient
    /// row of the highest power of the first variable is not the constant 1).
    NotMonicInFirstVariable,
    /// A bivariate input is not squarefree as a polynomial in its first
    /// variable over `Q(x)`, or no squarefree specialization was found.
    NotSquarefree,
    /// The recombination step needed more subset trials than allowed.
    RecombinationBudgetExceeded {
        /// The subset-trial limit that was hit.
        limit: usize,
    },
}

impl fmt::Display for PolynomialFactorizationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroPolynomial => write!(formatter, "cannot factor the zero polynomial"),
            Self::NotMonicInFirstVariable => {
                write!(
                    formatter,
                    "the bivariate polynomial is not monic in its first variable"
                )
            }
            Self::NotSquarefree => write!(
                formatter,
                "the polynomial is not squarefree in its first variable"
            ),
            Self::RecombinationBudgetExceeded { limit } => {
                write!(
                    formatter,
                    "factor recombination exceeded the subset limit {limit}"
                )
            }
        }
    }
}

impl std::error::Error for PolynomialFactorizationError {}

/// Default limit on the number of factor subsets tried during recombination.
pub const DEFAULT_RECOMBINATION_SUBSET_LIMIT: usize = 1 << 16;

/// Factorization of a univariate polynomial over `Q`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RationalPolynomialFactorization {
    /// Rational constant `c` with `input = c * prod factor^multiplicity`.
    pub constant: BigRational,
    /// Irreducible factors over `Q`: primitive integer coefficient vectors in
    /// ascending degree order with positive leading coefficient, each of
    /// degree at least one, with multiplicity at least one, sorted by degree
    /// and then coefficient vector.
    pub factors: Vec<(Vec<BigInt>, usize)>,
}

/// One irreducible factor of a bivariate polynomial that is monic in its first
/// variable.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MonicBivariateFactor {
    /// `coefficients[i][j]` is the coefficient of `s^i x^j`; the factor is
    /// monic in `s`: the last row is `[1]`.  Rows have trailing zeros removed
    /// (a zero row is the empty vector).
    pub coefficients: Vec<Vec<BigRational>>,
    /// Multiplicity of this factor in the input.
    pub multiplicity: usize,
}

/// Factor an integer-coefficient polynomial over `Q`.
///
/// Coefficients are in ascending degree order.  The result lists the
/// irreducible factors with multiplicities and the rational constant with
/// `input = constant * prod factor^multiplicity`.  Recombination uses
/// [`DEFAULT_RECOMBINATION_SUBSET_LIMIT`].
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::factorization::irreducible_factorization_bigint_coeffs;
///
/// // 2 x^5 + 2 x = 2 x (x^4 + 1); x^4 + 1 is irreducible over Q.
/// let coeffs: Vec<BigInt> = [0, 2, 0, 0, 0, 2].into_iter().map(BigInt::from).collect();
/// let factorization = irreducible_factorization_bigint_coeffs(&coeffs).unwrap();
/// assert_eq!(factorization.constant, polytool::BigRational::from_integer(2.into()));
/// let degrees: Vec<usize> = factorization.factors.iter().map(|(f, _)| f.len() - 1).collect();
/// assert_eq!(degrees, vec![1, 4]);
/// ```
pub fn irreducible_factorization_bigint_coeffs(
    coeffs: &[BigInt],
) -> Result<RationalPolynomialFactorization, PolynomialFactorizationError> {
    univariate::factor_bigint(coeffs, DEFAULT_RECOMBINATION_SUBSET_LIMIT)
}

/// Factor a rational-coefficient polynomial over `Q`.
///
/// Same output conventions as [`irreducible_factorization_bigint_coeffs`];
/// the constant absorbs the denominators.
pub fn irreducible_factorization_rational_coeffs(
    coeffs: &[BigRational],
) -> Result<RationalPolynomialFactorization, PolynomialFactorizationError> {
    univariate::factor_rational(coeffs, DEFAULT_RECOMBINATION_SUBSET_LIMIT)
}

/// Factor a rational-coefficient polynomial over `Q` with a recombination
/// subset-trial limit.
pub fn irreducible_factorization_rational_coeffs_with_budget(
    coeffs: &[BigRational],
    max_recombination_subsets: usize,
) -> Result<RationalPolynomialFactorization, PolynomialFactorizationError> {
    univariate::factor_rational(coeffs, max_recombination_subsets)
}

/// Factor a monic-in-`s` bivariate polynomial over `Q`.
///
/// `coefficients[i][j]` is the coefficient of `s^i x^j`; the row of the top
/// power of `s` must be the constant 1.  The input must be squarefree as a
/// polynomial in `s` over `Q(x)`; repeated factors give
/// [`PolynomialFactorizationError::NotSquarefree`].  The factors are sorted by
/// `s`-degree, then `x`-degree, then coefficients, and their product is the
/// input.
///
/// ```
/// use polytool::factorization::irreducible_factorization_monic_bivariate;
/// use polytool::BigRational;
///
/// let q = |n: i64| BigRational::from_integer(n.into());
/// // (s - 1)^2 - x^2 = (s - 1 - x)(s - 1 + x)
/// let input = vec![vec![q(1), q(0), q(-1)], vec![q(-2)], vec![q(1)]];
/// let factors = irreducible_factorization_monic_bivariate(&input).unwrap();
/// assert_eq!(factors.len(), 2);
/// assert!(factors.iter().all(|factor| factor.multiplicity == 1));
/// ```
pub fn irreducible_factorization_monic_bivariate(
    coefficients: &[Vec<BigRational>],
) -> Result<Vec<MonicBivariateFactor>, PolynomialFactorizationError> {
    bivariate::factor(coefficients, DEFAULT_RECOMBINATION_SUBSET_LIMIT)
}

/// Factor a monic-in-`s` bivariate polynomial with a recombination
/// subset-trial limit.
pub fn irreducible_factorization_monic_bivariate_with_budget(
    coefficients: &[Vec<BigRational>],
    max_recombination_subsets: usize,
) -> Result<Vec<MonicBivariateFactor>, PolynomialFactorizationError> {
    bivariate::factor(coefficients, max_recombination_subsets)
}
