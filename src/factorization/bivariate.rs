use num_bigint::BigInt;
use num_traits::{One, Zero};

use super::arithmetic::{add_q, div_rem_q, evaluate_q, mul_q, scale_q, sub_q, trim_q, Q};
use super::univariate::factor_rational;
use super::{MonicBivariateFactor, PolynomialFactorizationError};

type Bivariate = Vec<Vec<Q>>;
type Series = Vec<Vec<Q>>;

pub(super) fn factor(
    coefficients: &[Vec<Q>],
    budget: usize,
) -> Result<Vec<MonicBivariateFactor>, PolynomialFactorizationError> {
    let input = normalize_input(coefficients)?;
    if input.len() == 1 {
        return Ok(Vec::new());
    }
    let x_degree = input
        .iter()
        .map(|row| row.len().saturating_sub(1))
        .max()
        .unwrap_or(0);
    let mut selected_specialization = None;
    for specialization in specialization_points().take(64) {
        let value = evaluate_at_x(&input, &specialization);
        let factorization = factor_rational(&value, budget)?;
        if factorization
            .factors
            .iter()
            .any(|(_, multiplicity)| *multiplicity != 1)
        {
            continue;
        }
        selected_specialization = Some((specialization, factorization));
        break;
    }
    let (x0, factorization) =
        selected_specialization.ok_or(PolynomialFactorizationError::NotSquarefree)?;
    let base_factors = factorization
        .factors
        .iter()
        .map(|(factor, _)| {
            let leading = factor
                .last()
                .expect("nonconstant univariate factor")
                .clone();
            factor
                .iter()
                .map(|coefficient| Q::new(coefficient.clone(), leading.clone()))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let target = shifted_target(&input, &x0, x_degree);
    let lifted = lift_factorization(&base_factors, &target, x_degree);
    recombine(&input, &lifted, &x0, x_degree, budget)
}

fn normalize_input(coefficients: &[Vec<Q>]) -> Result<Bivariate, PolynomialFactorizationError> {
    let mut result = coefficients
        .iter()
        .map(|row| trim_q(row.clone()))
        .collect::<Vec<_>>();
    while result.last().is_some_and(Vec::is_empty) {
        result.pop();
    }
    if result.is_empty() {
        return Err(PolynomialFactorizationError::ZeroPolynomial);
    }
    if result.last() != Some(&vec![Q::one()]) {
        return Err(PolynomialFactorizationError::NotMonicInFirstVariable);
    }
    Ok(result)
}

fn specialization_points() -> impl Iterator<Item = BigInt> {
    std::iter::once(BigInt::zero()).chain((1usize..).flat_map(|magnitude| {
        let magnitude = BigInt::from(magnitude);
        [magnitude.clone(), -magnitude].into_iter()
    }))
}

fn evaluate_at_x(input: &[Vec<Q>], x: &BigInt) -> Vec<Q> {
    let x = Q::from_integer(x.clone());
    input.iter().map(|row| evaluate_q(row, &x)).collect()
}

fn binomial(n: usize, k: usize) -> BigInt {
    if k > n {
        return BigInt::zero();
    }
    (0..k).fold(BigInt::one(), |result, index| {
        result * BigInt::from(n - index) / BigInt::from(index + 1)
    })
}

fn shifted_target(input: &[Vec<Q>], x0: &BigInt, degree: usize) -> Vec<Vec<Q>> {
    let x0 = Q::from_integer(x0.clone());
    let mut target = vec![Vec::new(); degree + 1];
    for (s_degree, row) in input.iter().enumerate() {
        for (t_degree, target_row) in target.iter_mut().enumerate().take(degree + 1) {
            let coefficient = row
                .iter()
                .enumerate()
                .skip(t_degree)
                .map(|(x_degree, value)| {
                    value
                        * binomial(x_degree, t_degree)
                        * x0.clone().pow((x_degree - t_degree) as i32)
                })
                .fold(Q::zero(), |sum, value| sum + value);
            if !coefficient.is_zero() {
                if target_row.len() <= s_degree {
                    target_row.resize(s_degree + 1, Q::zero());
                }
                target_row[s_degree] = coefficient;
            }
        }
    }
    target.into_iter().map(trim_q).collect()
}

fn inverse_mod_q(value: &[Q], modulus: &[Q]) -> Vec<Q> {
    let mut old_r = trim_q(value.to_vec());
    let mut r = trim_q(modulus.to_vec());
    let mut old_s = vec![Q::one()];
    let mut s = Vec::new();
    while !r.is_empty() {
        let (quotient, remainder) = div_rem_q(&old_r, &r);
        old_r = r;
        r = remainder;
        let next_s = sub_q(&old_s, &mul_q(&quotient, &s));
        old_s = s;
        s = next_s;
    }
    assert_eq!(old_r.len(), 1, "specialized factors are not coprime");
    scale_q(&old_s, &(Q::one() / &old_r[0]))
}

fn remainder_mod_q(value: &[Q], modulus: &[Q]) -> Vec<Q> {
    div_rem_q(value, modulus).1
}

fn series_product_coefficient(factors: &[Series], degree: usize) -> Vec<Q> {
    let mut result = vec![Vec::new(); degree + 1];
    result[0] = vec![Q::one()];
    for factor in factors {
        let mut next = vec![Vec::new(); degree + 1];
        for left_degree in 0..=degree {
            for right_degree in 0..=degree - left_degree {
                let product = mul_q(&result[left_degree], &factor[right_degree]);
                next[left_degree + right_degree] =
                    add_q(&next[left_degree + right_degree], &product);
            }
        }
        result = next;
    }
    result[degree].clone()
}

fn lift_factorization(base: &[Vec<Q>], target: &[Vec<Q>], degree: usize) -> Vec<Series> {
    let inverses = base
        .iter()
        .enumerate()
        .map(|(index, factor)| {
            let cofactor = base
                .iter()
                .enumerate()
                .filter(|(other, _)| *other != index)
                .fold(vec![Q::one()], |product, (_, other)| mul_q(&product, other));
            inverse_mod_q(&cofactor, factor)
        })
        .collect::<Vec<_>>();
    let mut lifted = base
        .iter()
        .map(|factor| {
            let mut series = vec![Vec::new(); degree + 1];
            series[0] = factor.clone();
            series
        })
        .collect::<Vec<_>>();
    for t_degree in 1..=degree {
        let current = series_product_coefficient(&lifted, t_degree);
        let error = sub_q(&target[t_degree], &current);
        for index in 0..lifted.len() {
            let correction = remainder_mod_q(&mul_q(&error, &inverses[index]), &base[index]);
            lifted[index][t_degree] = correction;
        }
    }
    lifted
}

fn series_product(factors: &[Series], degree: usize) -> Series {
    (0..=degree)
        .map(|t_degree| series_product_coefficient(factors, t_degree))
        .collect()
}

fn series_to_bivariate(series: &[Vec<Q>], x0: &BigInt) -> Bivariate {
    let x0 = Q::from_integer(x0.clone());
    let s_degree = series.iter().map(|row| row.len()).max().unwrap_or(0);
    let mut result = vec![Vec::new(); s_degree];
    for (t_degree, row) in series.iter().enumerate() {
        for (s_degree, coefficient) in row.iter().enumerate() {
            if coefficient.is_zero() {
                continue;
            }
            let result_row = &mut result[s_degree];
            if result_row.len() < t_degree + 1 {
                result_row.resize(t_degree + 1, Q::zero());
            }
            for (x_degree, result_coefficient) in
                result_row.iter_mut().enumerate().take(t_degree + 1)
            {
                let term = coefficient
                    * binomial(t_degree, x_degree)
                    * (-x0.clone()).pow((t_degree - x_degree) as i32);
                *result_coefficient += term;
            }
        }
    }
    result.into_iter().map(trim_q).collect()
}

fn bivariate_div_rem(dividend: &[Vec<Q>], divisor: &[Vec<Q>]) -> (Bivariate, Bivariate) {
    let divisor_degree = divisor.len() - 1;
    let mut remainder = dividend.to_vec();
    let mut quotient = vec![Vec::new(); dividend.len().saturating_sub(divisor_degree)];
    while remainder.len() >= divisor.len() && !remainder.is_empty() {
        let shift = remainder.len() - divisor.len();
        let coefficient = remainder.last().cloned().unwrap_or_default();
        if quotient.len() <= shift {
            quotient.resize(shift + 1, Vec::new());
        }
        quotient[shift] = add_q(&quotient[shift], &coefficient);
        for (index, divisor_row) in divisor.iter().enumerate() {
            let product = mul_q(&coefficient, divisor_row);
            remainder[index + shift] = sub_q(&remainder[index + shift], &product);
        }
        while remainder.last().is_some_and(Vec::is_empty) {
            remainder.pop();
        }
    }
    (quotient.into_iter().map(trim_q).collect(), remainder)
}

#[allow(clippy::too_many_arguments)]
fn find_bivariate_divisor_subset(
    factors: &[Series],
    remainder: &[Vec<Q>],
    degree: usize,
    subset_size: usize,
    start: usize,
    selected: &mut Vec<usize>,
    x0: &BigInt,
    budget: usize,
    trials: &mut usize,
) -> Result<Option<(Series, Bivariate, Vec<usize>)>, PolynomialFactorizationError> {
    if selected.len() == subset_size {
        if *trials >= budget {
            return Err(PolynomialFactorizationError::RecombinationBudgetExceeded {
                limit: budget,
            });
        }
        *trials += 1;
        let chosen = selected
            .iter()
            .map(|index| factors[*index].clone())
            .collect::<Vec<_>>();
        let candidate_series = series_product(&chosen, degree);
        let candidate = series_to_bivariate(&candidate_series, x0);
        let (quotient, remainder_of_division) = bivariate_div_rem(remainder, &candidate);
        if remainder_of_division.is_empty() {
            return Ok(Some((candidate_series, quotient, selected.clone())));
        }
        return Ok(None);
    }
    let needed = subset_size - selected.len();
    let last_start = factors.len().saturating_sub(needed);
    for index in start..=last_start {
        selected.push(index);
        if let Some(found) = find_bivariate_divisor_subset(
            factors,
            remainder,
            degree,
            subset_size,
            index + 1,
            selected,
            x0,
            budget,
            trials,
        )? {
            selected.pop();
            return Ok(Some(found));
        }
        selected.pop();
    }
    Ok(None)
}

fn recombine(
    input: &[Vec<Q>],
    lifted: &[Series],
    x0: &BigInt,
    degree: usize,
    budget: usize,
) -> Result<Vec<MonicBivariateFactor>, PolynomialFactorizationError> {
    let mut active = lifted.to_vec();
    let mut remainder = input.to_vec();
    let mut trials = 0usize;
    let mut factors = Vec::new();
    while !active.is_empty() {
        let mut found = None;
        let mut selected = Vec::new();
        for subset_size in 1..=active.len() / 2 {
            if let Some(candidate) = find_bivariate_divisor_subset(
                &active,
                &remainder,
                degree,
                subset_size,
                0,
                &mut selected,
                x0,
                budget,
                &mut trials,
            )? {
                found = Some(candidate);
                break;
            }
        }
        if let Some((candidate_series, quotient, selected_indices)) = found {
            let candidate = series_to_bivariate(&candidate_series, x0);
            factors.push(MonicBivariateFactor {
                coefficients: candidate,
                multiplicity: 1,
            });
            for index in selected_indices.into_iter().rev() {
                active.remove(index);
            }
            remainder = quotient;
        } else {
            factors.push(MonicBivariateFactor {
                coefficients: remainder.clone(),
                multiplicity: 1,
            });
            active.clear();
        }
    }
    factors.sort_by(|left, right| {
        let left_key = (
            left.coefficients.len().saturating_sub(1),
            left.coefficients
                .iter()
                .map(|row| row.len().saturating_sub(1))
                .max()
                .unwrap_or(0),
            &left.coefficients,
        );
        let right_key = (
            right.coefficients.len().saturating_sub(1),
            right
                .coefficients
                .iter()
                .map(|row| row.len().saturating_sub(1))
                .max()
                .unwrap_or(0),
            &right.coefficients,
        );
        left_key.cmp(&right_key)
    });
    Ok(factors)
}
