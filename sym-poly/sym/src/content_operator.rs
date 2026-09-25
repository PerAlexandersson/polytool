//! The Schur-diagonal content operator `C_n` and its rational inverse.
//!
//! With rows and columns numbered from one, the content of `(i,j)` is `j-i`.
//! Thus `C_n(s_lambda) = product_(i,j in lambda)(n+j-i) * s_lambda`.
//! The empty product is one. For nonnegative `n`, the eigenvalue vanishes
//! exactly when `lambda` has more than `n` rows. On symmetric polynomials in
//! `n` variables only partitions with at most `n` rows survive, so the operator
//! is invertible over the rationals there.
//!
//! Results are returned in the Schur basis. The library represents formal
//! symmetric functions, not a finite-variable quotient: [`c_n_inverse`]
//! rejects a nonzero Schur term longer than `n` rather than silently discarding
//! it. To invert a finite-variable specialization, first convert to Schur and
//! retain only terms with at most `n` rows.
//!
//! Content products and rational division are arbitrary precision. Conversion
//! from other bases still uses the existing i64 transition tables and inherits
//! their degree/coefficient limits. Schur input avoids those conversions.

use std::error::Error;
use std::fmt;

use num_bigint::BigInt;
use num_rational::BigRational;
use sym_poly_core::{Partition, Ring};

use crate::{Basis, SymmetricFunction};

/// A nonzero Schur term lies in the kernel of `C_n`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ContentOperatorError {
    /// Index of a nonzero term having more than `num_vars` rows.
    pub partition: Partition,
    /// Number of variables supplied to the operator.
    pub num_vars: usize,
}

impl fmt::Display for ContentOperatorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "C_{} has zero eigenvalue on s_{} ({} rows)",
            self.num_vars,
            self.partition,
            self.partition.num_parts()
        )
    }
}

impl Error for ContentOperatorError {}

/// Return the integer eigenvalue `product_(i,j in lambda)(num_vars+j-i)`.
///
/// Rows and columns in this formula are one-based. Returns one for the empty
/// partition, and zero exactly when `lambda.num_parts() > num_vars`.
/// Arithmetic does not narrow `num_vars` to a signed machine integer.
pub fn schur_content_product(lambda: &Partition, num_vars: usize) -> BigInt {
    if lambda.num_parts() > num_vars {
        return BigInt::zero();
    }
    let mut result = BigInt::one();
    for (row, &length) in lambda.parts().iter().enumerate() {
        let first = BigInt::from(num_vars - row);
        for column in 0..length {
            result *= &first + BigInt::from(column);
        }
    }
    result
}

/// Apply `C_n` to a symmetric function in any classical basis.
///
/// Returns a Schur expansion with the same coefficient ring. Schur terms with
/// more than `num_vars` rows are annihilated. Use `BigInt` or `BigRational` to
/// avoid fixed-width coefficient overflow; other rings must support the
/// eigenvalues through [`Ring::from_bigint`].
pub fn c_n<C: Ring>(value: &SymmetricFunction<C>, num_vars: usize) -> SymmetricFunction<C> {
    SymmetricFunction::from_terms(
        Basis::Schur,
        value
            .to_schur_basis()
            .into_terms()
            .into_iter()
            .filter_map(|(partition, coefficient)| {
                let factor = schur_content_product(&partition, num_vars);
                (!factor.is_zero()).then(|| (partition, coefficient * C::from_bigint(&factor)))
            })
            .collect(),
    )
}

/// Apply `C_n^{-1}`, returning exact rational coefficients in the Schur basis.
///
/// Accepts `BigInt` or `BigRational` coefficients (or any `Ring` convertible
/// into `BigRational`). In particular, integer graph polynomials need no
/// manual rational lift. Coefficients are lifted before basis conversion.
///
/// # Errors
///
/// Returns [`ContentOperatorError`] if the Schur expansion contains a nonzero
/// term with more than `num_vars` rows. Such a formal symmetric function is not
/// in the image of `C_n`. For finite-variable input, discard those terms after
/// Schur conversion before calling this function. Zero is accepted for every
/// `num_vars`, and constants are unchanged, including when `num_vars = 0`.
///
/// # Examples
///
/// ```
/// use num_bigint::BigInt;
/// use num_rational::BigRational;
/// use sym_poly_core::Partition;
/// use sym_poly_sym::{c_n, c_n_inverse, SymmetricFunction};
///
/// let lambda = Partition::new(vec![2, 1]);
/// let s = SymmetricFunction::<BigInt>::schur_symmetric(lambda.clone());
/// assert_eq!(c_n(&s, 3).coefficient(&lambda), BigInt::from(24));
/// assert_eq!(
///     c_n_inverse(&s, 3)?.coefficient(&lambda),
///     BigRational::new(BigInt::from(1), BigInt::from(24)),
/// );
/// # Ok::<(), sym_poly_sym::ContentOperatorError>(())
/// ```
pub fn c_n_inverse<C>(
    value: &SymmetricFunction<C>,
    num_vars: usize,
) -> Result<SymmetricFunction<BigRational>, ContentOperatorError>
where
    C: Ring + Into<BigRational>,
{
    let rational: SymmetricFunction<BigRational> = SymmetricFunction::from_terms(
        value.basis(),
        value
            .terms()
            .iter()
            .map(|(partition, coefficient)| (partition.clone(), coefficient.clone().into()))
            .collect(),
    );
    let mut terms = rational.to_schur_basis().into_terms();
    for (partition, coefficient) in &mut terms {
        let factor = schur_content_product(partition, num_vars);
        if factor.is_zero() {
            return Err(ContentOperatorError {
                partition: partition.clone(),
                num_vars,
            });
        }
        *coefficient /= BigRational::from_integer(factor);
    }
    Ok(SymmetricFunction::from_terms(Basis::Schur, terms))
}
