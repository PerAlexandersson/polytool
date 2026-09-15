//! Symmetric-decomposition utilities inspired by Brandén--Solus.
//!
//! This module ties together:
//!
//! - the `I_d`-decomposition `p = a + x b`,
//! - the corresponding `f`-polynomial and `R_d`-decomposition,
//! - alternatingly increasing checks,
//! - magic-basis coordinates and the partial-sum inequalities from
//!   Theorem 2.13 of Brandén--Solus,
//!
//! so the web frontend can inspect the full package of decomposition data for a
//! given polynomial in one shot.

use crate::basis::{
    analyze_magic_basis_bigint, analyze_magic_basis_i64, BasisError, MagicBasisAnalysis,
};
use crate::polynomial::CoeffRing;
use crate::polynomial::Polynomial;
use crate::real_rootedness::{
    check_weak_interlacing, check_weak_interlacing_bigint_coeffs, is_real_rooted,
    is_real_rooted_bigint_coeffs,
};
use num_bigint::BigInt;
use num_rational::Ratio;

/// Exact rational type used in magic-basis analysis.
pub type BigRational = Ratio<BigInt>;

/// A structured view of the `I_d` / `R_d` decomposition data of an `i64` polynomial.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SymmetricDecompositionAnalysis {
    pub degree: usize,
    pub reciprocal: Vec<i64>,
    pub a: Vec<i64>,
    pub b: Vec<i64>,
    pub a_real_rooted: bool,
    pub b_real_rooted: bool,
    pub b_interlaces_a: Option<bool>,
    pub reciprocal_interlaces_input: Option<bool>,
    pub alternatingly_increasing: bool,
    pub f_polynomial: Vec<i64>,
    pub r_transform_of_f: Vec<i64>,
    pub r_a: Vec<i64>,
    pub r_b: Vec<i64>,
    pub r_interlaces_f: Option<bool>,
    pub magic: MagicBasisAnalysis<BigRational>,
}

/// A structured view of the `I_d` / `R_d` decomposition data of a `BigInt`
/// polynomial.
///
/// This is the arbitrary-precision counterpart of
/// [`SymmetricDecompositionAnalysis`]. No coefficient is narrowed to `i64`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SymmetricDecompositionAnalysisBigInt {
    pub degree: usize,
    pub reciprocal: Vec<BigInt>,
    pub a: Vec<BigInt>,
    pub b: Vec<BigInt>,
    pub a_real_rooted: bool,
    pub b_real_rooted: bool,
    pub b_interlaces_a: Option<bool>,
    pub reciprocal_interlaces_input: Option<bool>,
    pub alternatingly_increasing: bool,
    pub f_polynomial: Vec<BigInt>,
    pub r_transform_of_f: Vec<BigInt>,
    pub r_a: Vec<BigInt>,
    pub r_b: Vec<BigInt>,
    pub r_interlaces_f: Option<bool>,
    pub magic: MagicBasisAnalysis<BigRational>,
}

/// Compute the `f`-polynomial
///
/// ```text
/// f(h; x) = (1 + x)^d h(x / (1 + x))
///         = sum_{k=0}^d h_k x^k (1 + x)^(d-k)
/// ```
///
/// with respect to the degree bound `degree`.
pub fn f_polynomial<C: CoeffRing>(h: &Polynomial<C>, degree: usize) -> Option<Polynomial<C>> {
    match h.degree() {
        Some(d) if d > degree => return None,
        _ => {}
    }

    let one_plus_x = Polynomial::<C>::from_i64_coeffs(&[1, 1]);
    let mut result = Polynomial::<C>::zero();

    for k in 0..=degree {
        let coeff = h.coeff(k);
        if coeff.is_zero() {
            continue;
        }
        let term = Polynomial::<C>::monomial(C::one(), k)
            * poly_pow(&one_plus_x, degree.saturating_sub(k));
        result = result + term.scale(&coeff);
    }

    Some(result)
}

/// Convenience wrapper for `i64` coefficients.
pub fn f_polynomial_i64(coeffs: &[i64], degree: usize) -> Option<Vec<i64>> {
    let h = Polynomial::<i64>::from_i64_coeffs(coeffs);
    f_polynomial(&h, degree).map(|p| p.coeffs().to_vec())
}

/// Compute the reflection
///
/// ```text
/// R_d(p)(x) = (-1)^d p(-1 - x)
/// ```
///
/// with respect to the degree bound `degree`.
pub fn r_transform<C: CoeffRing>(p: &Polynomial<C>, degree: usize) -> Option<Polynomial<C>> {
    match p.degree() {
        Some(d) if d > degree => return None,
        _ => {}
    }

    let minus_one = C::from_i64(-1);
    let sign = if degree.is_multiple_of(2) {
        C::one()
    } else {
        C::from_i64(-1)
    };

    Some(p.shift(&minus_one).dilate(&minus_one).scale(&sign))
}

/// Convenience wrapper for `i64` coefficients.
pub fn r_transform_i64(coeffs: &[i64], degree: usize) -> Option<Vec<i64>> {
    let p = Polynomial::<i64>::from_i64_coeffs(coeffs);
    r_transform(&p, degree).map(|q| q.coeffs().to_vec())
}

/// Compute the `R_d`-decomposition `p = \tilde a + x \tilde b`.
pub fn r_decomposition<C: CoeffRing>(
    p: &Polynomial<C>,
    degree: usize,
) -> Option<(Polynomial<C>, Polynomial<C>)> {
    let reflected = r_transform(p, degree)?;
    let one_plus_x = Polynomial::<C>::from_i64_coeffs(&[1, 1]);
    let x = Polynomial::<C>::variable();

    let a_tilde = one_plus_x * p.clone() - x.clone() * reflected.clone();
    let b_tilde = reflected - p.clone();
    Some((a_tilde, b_tilde))
}

/// Convenience wrapper for `i64` coefficients.
pub fn r_decomposition_i64(coeffs: &[i64], degree: usize) -> Option<(Vec<i64>, Vec<i64>)> {
    let p = Polynomial::<i64>::from_i64_coeffs(coeffs);
    r_decomposition(&p, degree).map(|(a, b)| (a.coeffs().to_vec(), b.coeffs().to_vec()))
}

/// Check whether the coefficients are alternatingly increasing.
///
/// That is, for `p_0 + p_1 x + ... + p_d x^d`, we test
///
/// ```text
/// 0 <= p_0 <= p_d <= p_1 <= p_{d-1} <= ...
/// ```
pub fn is_alternatingly_increasing(coeffs: &[i64]) -> bool {
    let d = match coeffs.iter().rposition(|&c| c != 0) {
        Some(d) => d,
        None => return true,
    };

    let mut zigzag = Vec::with_capacity(d + 1);
    for j in 0..=d / 2 {
        zigzag.push(coeffs[j]);
        let mirror = d - j;
        if mirror != j {
            zigzag.push(coeffs[mirror]);
        }
    }

    zigzag.windows(2).all(|w| w[0] <= w[1])
}

/// Check whether arbitrary-precision integer coefficients are alternatingly
/// increasing.
pub fn is_alternatingly_increasing_bigint(coeffs: &[BigInt]) -> bool {
    let d = match coeffs.iter().rposition(|c| c != &BigInt::from(0)) {
        Some(d) => d,
        None => return true,
    };

    let mut zigzag = Vec::with_capacity(d + 1);
    for j in 0..=d / 2 {
        zigzag.push(&coeffs[j]);
        let mirror = d - j;
        if mirror != j {
            zigzag.push(&coeffs[mirror]);
        }
    }

    zigzag.windows(2).all(|w| w[0] <= w[1])
}

/// Analyze the decomposition data of an `i64`-coefficient polynomial with
/// respect to its actual degree.
pub fn analyze_symmetric_decomposition_i64(
    coeffs: &[i64],
) -> Result<SymmetricDecompositionAnalysis, BasisError> {
    let p = Polynomial::<i64>::from_i64_coeffs(coeffs);
    let degree = p.degree().unwrap_or(0);
    let input = p.coeffs().to_vec();

    let reciprocal = p
        .reverse_with_degree(degree)
        .expect("actual degree should always be a valid bound");
    let (a, b) = p
        .stapledon_decomposition(degree)
        .expect("actual degree should always be a valid bound");

    let f = f_polynomial(&p, degree).expect("actual degree should always be a valid bound");
    let r_of_f = r_transform(&f, degree).expect("actual degree should always be a valid bound");
    let (r_a, r_b) =
        r_decomposition(&f, degree).expect("actual degree should always be a valid bound");

    let magic = analyze_magic_basis_i64(p.coeffs(), degree)?;

    Ok(SymmetricDecompositionAnalysis {
        degree,
        reciprocal: reciprocal.coeffs().to_vec(),
        a: a.coeffs().to_vec(),
        b: b.coeffs().to_vec(),
        a_real_rooted: is_real_rooted(a.coeffs()),
        b_real_rooted: is_real_rooted(b.coeffs()),
        b_interlaces_a: check_weak_interlacing(b.coeffs(), a.coeffs()),
        reciprocal_interlaces_input: check_weak_interlacing(reciprocal.coeffs(), &input),
        alternatingly_increasing: is_alternatingly_increasing(&input),
        f_polynomial: f.coeffs().to_vec(),
        r_transform_of_f: r_of_f.coeffs().to_vec(),
        r_a: r_a.coeffs().to_vec(),
        r_b: r_b.coeffs().to_vec(),
        r_interlaces_f: check_weak_interlacing(r_of_f.coeffs(), f.coeffs()),
        magic,
    })
}

/// Analyze the decomposition data of a `BigInt`-coefficient polynomial with
/// respect to its actual degree, without narrowing any coefficient.
pub fn analyze_symmetric_decomposition_bigint(
    coeffs: &[BigInt],
) -> Result<SymmetricDecompositionAnalysisBigInt, BasisError> {
    let p = Polynomial::<BigInt>::new(coeffs.to_vec());
    let degree = p.degree().unwrap_or(0);
    let input = p.coeffs().to_vec();

    let reciprocal = p
        .reverse_with_degree(degree)
        .expect("actual degree should always be a valid bound");
    let (a, b) = p
        .stapledon_decomposition(degree)
        .expect("actual degree should always be a valid bound");

    let f = f_polynomial(&p, degree).expect("actual degree should always be a valid bound");
    let r_of_f = r_transform(&f, degree).expect("actual degree should always be a valid bound");
    let (r_a, r_b) =
        r_decomposition(&f, degree).expect("actual degree should always be a valid bound");

    let magic = analyze_magic_basis_bigint(p.coeffs(), degree)?;

    Ok(SymmetricDecompositionAnalysisBigInt {
        degree,
        reciprocal: reciprocal.coeffs().to_vec(),
        a: a.coeffs().to_vec(),
        b: b.coeffs().to_vec(),
        a_real_rooted: is_real_rooted_bigint_coeffs(a.coeffs()),
        b_real_rooted: is_real_rooted_bigint_coeffs(b.coeffs()),
        b_interlaces_a: check_weak_interlacing_bigint_coeffs(b.coeffs(), a.coeffs()),
        reciprocal_interlaces_input: check_weak_interlacing_bigint_coeffs(
            reciprocal.coeffs(),
            &input,
        ),
        alternatingly_increasing: is_alternatingly_increasing_bigint(&input),
        f_polynomial: f.coeffs().to_vec(),
        r_transform_of_f: r_of_f.coeffs().to_vec(),
        r_a: r_a.coeffs().to_vec(),
        r_b: r_b.coeffs().to_vec(),
        r_interlaces_f: check_weak_interlacing_bigint_coeffs(r_of_f.coeffs(), f.coeffs()),
        magic,
    })
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_f_polynomial_basic() {
        let f = f_polynomial_i64(&[1, 2, 1], 2).unwrap();
        assert_eq!(f, vec![1, 4, 4]);
    }

    #[test]
    fn test_r_decomposition_matches_paper_example() {
        let f = vec![1, 4, 6];
        let reflected = r_transform_i64(&f, 2).unwrap();
        assert_eq!(reflected, vec![3, 8, 6]);

        let (a_tilde, b_tilde) = r_decomposition_i64(&f, 2).unwrap();
        assert_eq!(a_tilde, vec![1, 2, 2]);
        assert_eq!(b_tilde, vec![2, 4]);
    }

    #[test]
    fn test_alternatingly_increasing_examples() {
        assert!(is_alternatingly_increasing(&[
            1, 1018, 10678, 14498, 2933, 32
        ]));
        assert!(is_alternatingly_increasing(&[1, 3, 1]));
        assert!(!is_alternatingly_increasing(&[1, 2, 3]));
    }

    #[test]
    fn test_symmetric_decomposition_analysis_matches_small_example() {
        let h = vec![1, 4, 2];
        let analysis = analyze_symmetric_decomposition_i64(&h).unwrap();

        assert_eq!(analysis.degree, 2);
        assert_eq!(analysis.reciprocal, vec![2, 4, 1]);
        assert_eq!(analysis.a, vec![1, 3, 1]);
        assert_eq!(analysis.b, vec![1, 1]);
        assert_eq!(analysis.f_polynomial, vec![1, 6, 7]);
        assert_eq!(analysis.r_transform_of_f, vec![2, 8, 7]);
        assert_eq!(analysis.r_a, vec![1, 5, 5]);
        assert_eq!(analysis.r_b, vec![1, 2]);
        assert!(analysis.a_real_rooted);
        assert!(analysis.b_real_rooted);
        assert_eq!(analysis.b_interlaces_a, Some(true));
        assert_eq!(analysis.reciprocal_interlaces_input, Some(true));
        assert_eq!(analysis.r_interlaces_f, Some(true));
        assert!(analysis.alternatingly_increasing);
        assert_eq!(
            analysis.magic.coordinates,
            vec![
                BigRational::from_integer(BigInt::from(1)),
                BigRational::from_integer(BigInt::from(2)),
                BigRational::from_integer(BigInt::from(-1)),
            ]
        );
        assert!(!analysis.magic.left_leq_right);
    }

    #[test]
    fn test_bigint_symmetric_decomposition_preserves_huge_coefficients() {
        let scale = BigInt::from(10u8).pow(40);
        let h = [1_i64, 4, 2]
            .into_iter()
            .map(|coefficient| BigInt::from(coefficient) * &scale)
            .collect::<Vec<_>>();
        let analysis = analyze_symmetric_decomposition_bigint(&h).unwrap();

        assert_eq!(analysis.degree, 2);
        assert_eq!(
            analysis.reciprocal,
            [2_i64, 4, 1]
                .into_iter()
                .map(|coefficient| BigInt::from(coefficient) * &scale)
                .collect::<Vec<_>>()
        );
        assert_eq!(
            analysis.a,
            [1_i64, 3, 1]
                .into_iter()
                .map(|coefficient| BigInt::from(coefficient) * &scale)
                .collect::<Vec<_>>()
        );
        assert_eq!(
            analysis.b,
            [1_i64, 1]
                .into_iter()
                .map(|coefficient| BigInt::from(coefficient) * &scale)
                .collect::<Vec<_>>()
        );
        assert_eq!(
            analysis.magic.coordinates[0],
            BigRational::from_integer(scale)
        );
        assert!(analysis.a_real_rooted);
        assert!(analysis.b_real_rooted);
        assert_eq!(analysis.b_interlaces_a, Some(true));
        assert!(analysis.alternatingly_increasing);
    }
}
