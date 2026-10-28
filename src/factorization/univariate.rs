use std::collections::BTreeMap;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};

use super::arithmetic::{
    derivative_q, div_rem_q, exact_div_q, gcd_q, is_one_q, monic_q, primitive_positive, scale_q,
    trim_q, Q,
};
use super::modular::{candidate_primes, factor_mod_prime, hensel_lift};
use super::{PolynomialFactorizationError, RationalPolynomialFactorization};

pub(super) fn factor_bigint(
    coefficients: &[BigInt],
    budget: usize,
) -> Result<RationalPolynomialFactorization, PolynomialFactorizationError> {
    let rational_coefficients = coefficients
        .iter()
        .cloned()
        .map(BigRational::from_integer)
        .collect::<Vec<_>>();
    factor_rational(&rational_coefficients, budget)
}

pub(super) fn factor_rational(
    coefficients: &[Q],
    budget: usize,
) -> Result<RationalPolynomialFactorization, PolynomialFactorizationError> {
    let coefficients = trim_q(coefficients.to_vec());
    if coefficients.is_empty() {
        return Err(PolynomialFactorizationError::ZeroPolynomial);
    }
    if coefficients.len() == 1 {
        return Ok(RationalPolynomialFactorization {
            constant: coefficients[0].clone(),
            factors: Vec::new(),
        });
    }
    let leading = coefficients.last().expect("nonconstant polynomial").clone();
    let monic = scale_q(&coefficients, &(Q::one() / &leading));
    let mut constant = leading;
    let mut factor_map: BTreeMap<Vec<BigInt>, usize> = BTreeMap::new();
    let mut subset_trials = 0usize;
    for (squarefree_factor, multiplicity) in squarefree_decomposition(&monic) {
        if squarefree_factor.len() <= 1 {
            continue;
        }
        for factor in factor_squarefree_monic(&squarefree_factor, budget, &mut subset_trials)? {
            let (primitive, scalar) = primitive_positive(&factor);
            constant *= scalar_pow(&scalar, multiplicity);
            *factor_map.entry(primitive).or_default() += multiplicity;
        }
    }
    Ok(RationalPolynomialFactorization {
        constant,
        factors: factor_map.into_iter().collect(),
    })
}

fn scalar_pow(scalar: &Q, exponent: usize) -> Q {
    (0..exponent).fold(Q::one(), |product, _| product * scalar)
}

fn squarefree_decomposition(polynomial: &[Q]) -> Vec<(Vec<Q>, usize)> {
    let derivative = derivative_q(polynomial);
    let mut common = gcd_q(polynomial, &derivative);
    let mut squarefree_part = exact_div_q(polynomial, &common);
    let mut multiplicity = 1;
    let mut result = Vec::new();
    while !is_one_q(&squarefree_part) {
        let next_common = gcd_q(&squarefree_part, &common);
        let factor = exact_div_q(&squarefree_part, &next_common);
        if !is_one_q(&factor) {
            result.push((monic_q(&factor), multiplicity));
        }
        squarefree_part = next_common;
        common = exact_div_q(&common, &squarefree_part);
        multiplicity += 1;
    }
    result
}

fn factor_squarefree_monic(
    polynomial: &[Q],
    budget: usize,
    subset_trials: &mut usize,
) -> Result<Vec<Vec<Q>>, PolynomialFactorizationError> {
    let polynomial = monic_q(polynomial);
    if polynomial.len() <= 2 {
        return Ok(vec![polynomial]);
    }
    let (integer_polynomial, _) = primitive_positive(&polynomial);
    let leading = integer_polynomial
        .last()
        .expect("nonconstant polynomial")
        .clone();
    let degree = integer_polynomial.len() - 1;
    let coefficient_sum = integer_polynomial
        .iter()
        .map(|coefficient| coefficient.abs())
        .fold(BigInt::zero(), |sum, coefficient| sum + coefficient);
    // Landau--Mignotte: every integer factor `h` of the primitive polynomial
    // has `|h|_inf <= 2^degree * |f|_1`.  Recombination reconstructs
    // `leading * g` for a monic factor `g`, which is an integer multiple of
    // such an `h` by at most `|leading|`, so the bound includes that factor.
    let bound = (BigInt::one() << (degree + 2)) * coefficient_sum * leading.abs() + BigInt::one();
    for prime in candidate_primes().filter(|prime| (&leading % *prime) != BigInt::zero()) {
        let Some(factors_mod_prime) = factor_mod_prime(&integer_polynomial, prime) else {
            continue;
        };
        // The prime does not divide the leading coefficient, so a polynomial
        // that stays irreducible modulo the prime is irreducible over Q.
        if factors_mod_prime.len() <= 1 {
            return Ok(vec![polynomial]);
        }
        let (lifted_factors, modulus) = hensel_lift(
            &integer_polynomial,
            &leading,
            prime,
            &factors_mod_prime,
            &(bound.clone() * 2u32),
        );
        let mut indices = Vec::new();
        for subset_size in 1..=lifted_factors.len() / 2 {
            if let Some((candidate, quotient)) = find_divisor_subset(
                &lifted_factors,
                &modulus,
                &leading,
                &polynomial,
                subset_size,
                0,
                &mut indices,
                budget,
                subset_trials,
            )? {
                let mut result = factor_squarefree_monic(&candidate, budget, subset_trials)?;
                result.extend(factor_squarefree_monic(&quotient, budget, subset_trials)?);
                return Ok(result);
            }
        }
        return Ok(vec![polynomial]);
    }
    // The ascending prime iterator is unbounded in the practical domain; this
    // fallback is only reachable if that iterator itself is exhausted.
    Ok(vec![polynomial])
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn find_divisor_subset(
    factors: &[Vec<BigInt>],
    modulus: &BigInt,
    leading: &BigInt,
    polynomial: &[Q],
    subset_size: usize,
    start: usize,
    selected: &mut Vec<usize>,
    budget: usize,
    subset_trials: &mut usize,
) -> Result<Option<(Vec<Q>, Vec<Q>)>, PolynomialFactorizationError> {
    if selected.len() == subset_size {
        if *subset_trials >= budget {
            return Err(PolynomialFactorizationError::RecombinationBudgetExceeded {
                limit: budget,
            });
        }
        *subset_trials += 1;
        let mut integer_product = vec![BigInt::one()];
        for index in selected.iter().copied() {
            let factor = &factors[index];
            let mut product = vec![BigInt::zero(); integer_product.len() + factor.len() - 1];
            for (left_degree, left) in integer_product.iter().enumerate() {
                for (right_degree, right) in factor.iter().enumerate() {
                    product[left_degree + right_degree] += left * right;
                }
            }
            for coefficient in &mut product {
                *coefficient = super::arithmetic::mod_bigint(coefficient, modulus);
            }
            integer_product = product;
        }
        let candidate = trim_q(
            integer_product
                .iter()
                .map(|coefficient| {
                    let numerator =
                        super::arithmetic::centered_mod_bigint(&(coefficient * leading), modulus);
                    BigRational::new(numerator, leading.clone())
                })
                .collect(),
        );
        let (quotient, remainder) = div_rem_q(polynomial, &candidate);
        if remainder.is_empty() && quotient.len() > 1 {
            return Ok(Some((monic_q(&candidate), monic_q(&quotient))));
        }
        return Ok(None);
    }
    let remaining = subset_size - selected.len();
    let last_start = factors.len().saturating_sub(remaining);
    for index in start..=last_start {
        selected.push(index);
        if let Some(found) = find_divisor_subset(
            factors,
            modulus,
            leading,
            polynomial,
            subset_size,
            index + 1,
            selected,
            budget,
            subset_trials,
        )? {
            selected.pop();
            return Ok(Some(found));
        }
        selected.pop();
    }
    Ok(None)
}
