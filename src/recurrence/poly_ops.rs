//! Exact univariate and bivariate coefficient-vector helpers shared by the
//! recurrence finders, verifiers and generators.

use super::*;
use num_traits::One;
use num_traits::Zero;

/// Convert an `i64` coefficient vector to exact rational coefficients.
pub(super) fn i64_poly_to_rational(coeffs: &[i64]) -> Vec<BigRational> {
    coeffs
        .iter()
        .map(|&c| BigRational::from_integer(BigInt::from(c)))
        .collect()
}

pub(super) fn i64_polys_to_rational(polys: &[Vec<i64>]) -> Vec<Vec<BigRational>> {
    polys.iter().map(|p| i64_poly_to_rational(p)).collect()
}

pub(super) fn rational_zero_poly() -> Vec<BigRational> {
    vec![BigRational::zero()]
}

/// d-th derivative of a polynomial given as a rational coefficient vector.
pub(super) fn poly_nth_derivative_rational(coeffs: &[BigRational], d: usize) -> Vec<BigRational> {
    let mut result = coeffs.to_vec();
    for _ in 0..d {
        result = poly_derivative_once_rational(&result);
    }
    if result.is_empty() {
        rational_zero_poly()
    } else {
        result
    }
}

pub(super) fn poly_derivative_once_rational(coeffs: &[BigRational]) -> Vec<BigRational> {
    if coeffs.len() <= 1 {
        return rational_zero_poly();
    }
    coeffs[1..]
        .iter()
        .enumerate()
        .map(|(i, c)| {
            let factor = BigInt::from(i + 1);
            if c.denom().is_one() {
                BigRational::from_integer(c.numer() * factor)
            } else if factor.is_one() {
                c.clone()
            } else {
                BigRational::new(c.numer() * factor, c.denom().clone())
            }
        })
        .collect()
}

/// Coefficient of t^k (0 when k is out of range).
pub(super) fn poly_coeff_rational(coeffs: &[BigRational], k: usize) -> BigRational {
    if k < coeffs.len() {
        coeffs[k].clone()
    } else {
        BigRational::zero()
    }
}

pub(super) fn poly_coeff_nonzero_rational(coeffs: &[BigRational], k: usize) -> bool {
    coeffs.get(k).is_some_and(|coeff| !coeff.is_zero())
}

/// Degree of a polynomial (returns 0 for the zero polynomial).
pub(super) fn poly_degree_rational(coeffs: &[BigRational]) -> usize {
    coeffs.iter().rposition(|c| !c.is_zero()).unwrap_or(0)
}

pub(super) fn trim_poly_rational(mut coeffs: Vec<BigRational>) -> Vec<BigRational> {
    while coeffs.len() > 1 && coeffs.last().is_some_and(|c| c.is_zero()) {
        coeffs.pop();
    }
    if coeffs.is_empty() {
        rational_zero_poly()
    } else {
        coeffs
    }
}

pub(super) fn poly_is_zero_rational(coeffs: &[BigRational]) -> bool {
    coeffs.iter().all(|c| c.is_zero())
}

pub(super) fn poly_add_scaled_assign(
    target: &mut Vec<BigRational>,
    source: &[BigRational],
    scale: &BigRational,
) {
    if scale.is_zero() || poly_is_zero_rational(source) {
        return;
    }
    if target.len() < source.len() {
        target.resize(source.len(), BigRational::zero());
    }
    for (idx, coeff) in source.iter().enumerate() {
        if !coeff.is_zero() {
            target[idx] += coeff.clone() * scale.clone();
        }
    }
}

pub(super) fn poly_add_signed_assign(
    target: &mut Vec<BigRational>,
    source: &[BigRational],
    subtract: bool,
) {
    if poly_is_zero_rational(source) {
        return;
    }
    if target.len() < source.len() {
        target.resize(source.len(), BigRational::zero());
    }
    for (idx, coeff) in source.iter().enumerate() {
        if coeff.is_zero() {
            continue;
        }
        if subtract {
            target[idx] -= coeff.clone();
        } else {
            target[idx] += coeff.clone();
        }
    }
}

pub(super) fn poly_add_bivar_eval_signed_assign(
    target: &mut Vec<BigRational>,
    poly: &BivarPoly,
    n_powers: &[BigRational],
    subtract: bool,
) {
    let width = poly.coeffs.iter().map(Vec::len).max().unwrap_or(0);
    if width == 0 {
        return;
    }
    if target.len() < width {
        target.resize(width, BigRational::zero());
    }
    for (i, row) in poly.coeffs.iter().enumerate() {
        let n_power = n_powers
            .get(i)
            .expect("bivariate evaluation requires enough powers of n");
        if n_power.is_zero() {
            continue;
        }
        for (j, coeff) in row.iter().enumerate() {
            if coeff.is_zero() {
                continue;
            }
            let value = coeff.clone() * n_power.clone();
            if subtract {
                target[j] -= value;
            } else {
                target[j] += value;
            }
        }
    }
}

pub(super) fn poly_add_bivar_eval_times_poly_signed_assign(
    target: &mut Vec<BigRational>,
    coeff_poly: &BivarPoly,
    n_powers: &[BigRational],
    source: &[BigRational],
    subtract: bool,
) {
    if poly_is_zero_rational(source) {
        return;
    }
    let width = coeff_poly.coeffs.iter().map(Vec::len).max().unwrap_or(0);
    if width == 0 {
        return;
    }
    let needed_len = width + source.len() - 1;
    if target.len() < needed_len {
        target.resize(needed_len, BigRational::zero());
    }

    for (i, row) in coeff_poly.coeffs.iter().enumerate() {
        let n_power = n_powers
            .get(i)
            .expect("bivariate evaluation requires enough powers of n");
        if n_power.is_zero() {
            continue;
        }
        for (j, coeff) in row.iter().enumerate() {
            if coeff.is_zero() {
                continue;
            }
            let scale = coeff.clone() * n_power.clone();
            if scale.is_zero() {
                continue;
            }
            for (k, source_coeff) in source.iter().enumerate() {
                if source_coeff.is_zero() {
                    continue;
                }
                let value = scale.clone() * source_coeff.clone();
                if subtract {
                    target[j + k] -= value;
                } else {
                    target[j + k] += value;
                }
            }
        }
    }
}

pub(super) fn poly_mul_rational(lhs: &[BigRational], rhs: &[BigRational]) -> Vec<BigRational> {
    if poly_is_zero_rational(lhs) || poly_is_zero_rational(rhs) {
        return rational_zero_poly();
    }
    let mut result = vec![BigRational::zero(); lhs.len() + rhs.len() - 1];
    for (i, a) in lhs.iter().enumerate() {
        if a.is_zero() {
            continue;
        }
        for (j, b) in rhs.iter().enumerate() {
            if !b.is_zero() {
                result[i + j] += a.clone() * b.clone();
            }
        }
    }
    trim_poly_rational(result)
}

pub(super) fn poly_div_exact_rational(
    numerator: &[BigRational],
    denominator: &[BigRational],
) -> Result<Vec<BigRational>, RecurrenceEvaluationError> {
    let denominator = trim_poly_rational(denominator.to_vec());
    if poly_is_zero_rational(&denominator) {
        return Err(RecurrenceEvaluationError::ZeroDenominator);
    }

    let mut remainder = trim_poly_rational(numerator.to_vec());
    if poly_is_zero_rational(&remainder) {
        return Ok(rational_zero_poly());
    }

    let denom_deg = poly_degree_rational(&denominator);
    let denom_lc = denominator[denom_deg].clone();
    let mut quotient =
        vec![BigRational::zero(); poly_degree_rational(&remainder).saturating_sub(denom_deg) + 1];

    while !poly_is_zero_rational(&remainder) && poly_degree_rational(&remainder) >= denom_deg {
        let rem_deg = poly_degree_rational(&remainder);
        let shift = rem_deg - denom_deg;
        let factor = remainder[rem_deg].clone() / denom_lc.clone();
        quotient[shift] += factor.clone();
        for (idx, denom_coeff) in denominator.iter().enumerate().take(denom_deg + 1) {
            let target = idx + shift;
            remainder[target] -= factor.clone() * denom_coeff.clone();
        }
        remainder = trim_poly_rational(remainder);
    }

    if poly_is_zero_rational(&remainder) {
        Ok(trim_poly_rational(quotient))
    } else {
        Err(RecurrenceEvaluationError::NonPolynomialQuotient)
    }
}

pub(super) fn bivar_eval_n(poly: &BivarPoly, n: usize) -> Vec<BigRational> {
    bivar_eval_n_with_powers(poly, &rational_index_powers(n, bivar_n_degree(poly)))
}

fn bivar_eval_n_with_powers(poly: &BivarPoly, n_powers: &[BigRational]) -> Vec<BigRational> {
    let width = poly.coeffs.iter().map(Vec::len).max().unwrap_or(0);
    if width == 0 {
        return rational_zero_poly();
    }
    let mut result = vec![BigRational::zero(); width];
    for (i, row) in poly.coeffs.iter().enumerate() {
        let n_power = n_powers
            .get(i)
            .expect("bivariate evaluation requires enough powers of n");
        for (j, coeff) in row.iter().enumerate() {
            if !coeff.is_zero() {
                result[j] += coeff.clone() * n_power.clone();
            }
        }
    }
    trim_poly_rational(result)
}

fn bivar_n_degree(poly: &BivarPoly) -> usize {
    poly.coeffs.len().saturating_sub(1)
}

pub(super) fn recurrence_n_degree(rec: &Recurrence) -> usize {
    let mut degree = rec.denominator.as_ref().map_or(0, bivar_n_degree);
    for term in &rec.terms {
        degree = degree.max(bivar_n_degree(&term.coeff));
    }
    if let Some(inhomogeneous) = &rec.inhomogeneous {
        degree = degree.max(bivar_n_degree(inhomogeneous));
    }
    degree
}
