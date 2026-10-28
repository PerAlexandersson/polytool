use num_bigint::BigInt;
use num_traits::{One, ToPrimitive, Zero};

use super::arithmetic::{centered_mod_bigint, mod_bigint};
type Fp = Vec<u64>;

fn trim(mut p: Fp) -> Fp {
    while p.last() == Some(&0) {
        p.pop();
    }
    p
}

fn sub(a: &[u64], b: &[u64], p: u64) -> Fp {
    let mut result = vec![0; a.len().max(b.len())];
    for (i, value) in result.iter_mut().enumerate() {
        let left = a.get(i).copied().unwrap_or(0);
        let right = b.get(i).copied().unwrap_or(0);
        *value = (left + p - right) % p;
    }
    trim(result)
}

fn mul(a: &[u64], b: &[u64], p: u64) -> Fp {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut result = vec![0; a.len() + b.len() - 1];
    for (i, left) in a.iter().copied().enumerate() {
        for (j, right) in b.iter().copied().enumerate() {
            result[i + j] = (result[i + j] + left * right) % p;
        }
    }
    trim(result)
}

fn scale(a: &[u64], scalar: u64, p: u64) -> Fp {
    trim(a.iter().map(|value| value * scalar % p).collect())
}

fn div_rem(dividend: &[u64], divisor: &[u64], p: u64) -> (Fp, Fp) {
    assert!(!divisor.is_empty(), "division by zero polynomial");
    let divisor_leading = *divisor.last().expect("nonempty divisor");
    let divisor_inverse = inverse_scalar(divisor_leading, p);
    let mut remainder = trim(dividend.to_vec());
    if remainder.len() < divisor.len() {
        return (Vec::new(), remainder);
    }
    let mut quotient = vec![0; remainder.len() - divisor.len() + 1];
    while remainder.len() >= divisor.len() && !remainder.is_empty() {
        let shift = remainder.len() - divisor.len();
        let coefficient = remainder.last().copied().unwrap() * divisor_inverse % p;
        quotient[shift] = coefficient;
        for (i, divisor_coefficient) in divisor.iter().copied().enumerate() {
            let subtraction = coefficient * divisor_coefficient % p;
            remainder[i + shift] = (remainder[i + shift] + p - subtraction) % p;
        }
        remainder = trim(remainder);
    }
    (trim(quotient), remainder)
}

fn exact_div(dividend: &[u64], divisor: &[u64], p: u64) -> Fp {
    let (quotient, remainder) = div_rem(dividend, divisor, p);
    assert!(remainder.is_empty(), "finite-field division was not exact");
    quotient
}

fn gcd(a: &[u64], b: &[u64], p: u64) -> Fp {
    let mut left = monic(a, p);
    let mut right = monic(b, p);
    while !right.is_empty() {
        let (_, remainder) = div_rem(&left, &right, p);
        left = right;
        right = remainder;
    }
    monic(&left, p)
}

fn monic(p: &[u64], prime: u64) -> Fp {
    let p = trim(p.to_vec());
    if p.is_empty() {
        return p;
    }
    scale(&p, inverse_scalar(*p.last().unwrap(), prime), prime)
}

fn derivative(p: &[u64], prime: u64) -> Fp {
    if p.len() <= 1 {
        return Vec::new();
    }
    trim(
        p.iter()
            .enumerate()
            .skip(1)
            .map(|(index, value)| value * (index as u64) % prime)
            .collect(),
    )
}

fn inverse_scalar(value: u64, prime: u64) -> u64 {
    assert!(value != 0 && value < prime);
    let mut old_r = value as i128;
    let mut r = prime as i128;
    let mut old_s = 1i128;
    let mut s = 0i128;
    while r != 0 {
        let quotient = old_r / r;
        (old_r, r) = (r, old_r - quotient * r);
        (old_s, s) = (s, old_s - quotient * s);
    }
    (old_s.rem_euclid(prime as i128)) as u64
}

fn pow_mod(base: &[u64], mut exponent: u64, modulus: &[u64], prime: u64) -> Fp {
    let mut result = vec![1];
    let mut base = div_rem(base, modulus, prime).1;
    while exponent != 0 {
        if exponent & 1 == 1 {
            result = div_rem(&mul(&result, &base, prime), modulus, prime).1;
        }
        exponent >>= 1;
        if exponent != 0 {
            base = div_rem(&mul(&base, &base, prime), modulus, prime).1;
        }
    }
    result
}

fn extended_gcd(a: &[u64], b: &[u64], prime: u64) -> (Fp, Fp, Fp) {
    let mut old_r = trim(a.to_vec());
    let mut r = trim(b.to_vec());
    let mut old_s = vec![1];
    let mut s = Vec::new();
    let mut old_t = Vec::new();
    let mut t = vec![1];
    while !r.is_empty() {
        let (quotient, remainder) = div_rem(&old_r, &r, prime);
        old_r = r;
        r = remainder;
        let next_s = sub(&old_s, &mul(&quotient, &s, prime), prime);
        old_s = s;
        s = next_s;
        let next_t = sub(&old_t, &mul(&quotient, &t, prime), prime);
        old_t = t;
        t = next_t;
    }
    let leading_inverse = inverse_scalar(*old_r.last().unwrap(), prime);
    (
        scale(&old_r, leading_inverse, prime),
        scale(&old_s, leading_inverse, prime),
        scale(&old_t, leading_inverse, prime),
    )
}

fn inverse_mod(value: &[u64], modulus: &[u64], prime: u64) -> Fp {
    let (gcd_value, coefficient, _) = extended_gcd(value, modulus, prime);
    assert_eq!(gcd_value, vec![1], "polynomial inverse does not exist");
    div_rem(&coefficient, modulus, prime).1
}

fn nullspace(mut matrix: Vec<Vec<u64>>, prime: u64) -> Vec<Fp> {
    let rows = matrix.len();
    let columns = matrix.first().map_or(0, Vec::len);
    let mut pivot_columns = Vec::new();
    let mut pivot_row = 0;
    for column in 0..columns {
        let Some(source_row) = (pivot_row..rows).find(|row| matrix[*row][column] != 0) else {
            continue;
        };
        matrix.swap(pivot_row, source_row);
        let inverse = inverse_scalar(matrix[pivot_row][column], prime);
        for value in &mut matrix[pivot_row] {
            *value = *value * inverse % prime;
        }
        for row in 0..rows {
            if row == pivot_row || matrix[row][column] == 0 {
                continue;
            }
            let multiple = matrix[row][column];
            let pivot_tail = matrix[pivot_row][column..].to_vec();
            for (offset, value) in matrix[row][column..].iter_mut().enumerate() {
                *value = (*value + prime - multiple * pivot_tail[offset] % prime) % prime;
            }
        }
        pivot_columns.push(column);
        pivot_row += 1;
        if pivot_row == rows {
            break;
        }
    }
    let is_pivot = (0..columns)
        .map(|column| pivot_columns.contains(&column))
        .collect::<Vec<_>>();
    let mut result = Vec::new();
    for free_column in 0..columns {
        if is_pivot[free_column] {
            continue;
        }
        let mut vector = vec![0; columns];
        vector[free_column] = 1;
        for (row, pivot_column) in pivot_columns.iter().copied().enumerate() {
            vector[pivot_column] = (prime - matrix[row][free_column]) % prime;
        }
        result.push(trim(vector));
    }
    result
}

fn berlekamp_basis(polynomial: &[u64], prime: u64) -> Vec<Fp> {
    let degree = polynomial.len() - 1;
    let x = vec![0, 1];
    let mut matrix = vec![vec![0; degree]; degree];
    let powers = (0..degree)
        .map(|column| pow_mod(&x, prime * column as u64, polynomial, prime))
        .collect::<Vec<_>>();
    for (column, power) in powers.iter().enumerate() {
        for (row, matrix_row) in matrix.iter_mut().enumerate() {
            matrix_row[column] = power.get(row).copied().unwrap_or(0);
        }
    }
    for (column, matrix_row) in matrix.iter_mut().enumerate() {
        matrix_row[column] = (matrix_row[column] + prime - 1) % prime;
    }
    nullspace(matrix, prime)
}

fn berlekamp_recursive(polynomial: Fp, prime: u64, output: &mut Vec<Fp>) {
    if polynomial.len() <= 2 {
        output.push(polynomial);
        return;
    }
    let basis = berlekamp_basis(&polynomial, prime);
    if basis.len() <= 1 {
        output.push(polynomial);
        return;
    }
    for basis_element in &basis {
        for scalar in 0..prime {
            let mut shifted = basis_element.clone();
            if shifted.is_empty() {
                shifted.push((prime - scalar) % prime);
            } else {
                shifted[0] = (shifted[0] + prime - scalar) % prime;
            }
            let divisor = gcd(&polynomial, &shifted, prime);
            if divisor.len() > 1 && divisor.len() < polynomial.len() {
                let quotient = exact_div(&polynomial, &divisor, prime);
                berlekamp_recursive(divisor, prime, output);
                berlekamp_recursive(quotient, prime, output);
                return;
            }
        }
    }
    panic!("Berlekamp basis did not split a reducible squarefree polynomial");
}

/// Factor a squarefree polynomial over a prime field into monic factors.
pub(super) fn factor_mod_prime(coefficients: &[BigInt], prime: u64) -> Option<Vec<Fp>> {
    let polynomial: Fp = coefficients
        .iter()
        .map(|coefficient| {
            mod_bigint(coefficient, &BigInt::from(prime))
                .to_u64()
                .unwrap()
        })
        .collect();
    let polynomial = trim(polynomial);
    if polynomial.len() <= 1 || polynomial.last() == Some(&0) {
        return None;
    }
    if gcd(&polynomial, &derivative(&polynomial, prime), prime).len() != 1 {
        return None;
    }
    let mut factors = Vec::new();
    berlekamp_recursive(monic(&polynomial, prime), prime, &mut factors);
    factors.sort_by_key(Vec::len);
    Some(factors)
}

fn product_bigint(factors: &[Vec<BigInt>]) -> Vec<BigInt> {
    factors.iter().fold(vec![BigInt::one()], |product, factor| {
        let mut result = vec![BigInt::zero(); product.len() + factor.len() - 1];
        for (i, left) in product.iter().enumerate() {
            for (j, right) in factor.iter().enumerate() {
                result[i + j] += left * right;
            }
        }
        result
    })
}

fn product_mod_prime_except(factors: &[Vec<BigInt>], skip: usize, prime: u64) -> Fp {
    factors
        .iter()
        .enumerate()
        .filter(|(index, _)| *index != skip)
        .fold(vec![1], |product, (_, factor)| {
            let reduced: Fp = factor
                .iter()
                .map(|coefficient| {
                    mod_bigint(coefficient, &BigInt::from(prime))
                        .to_u64()
                        .expect("small prime residue")
                })
                .collect();
            mul(&product, &reduced, prime)
        })
}

/// Multifactor Hensel lifting for a primitive integer polynomial whose leading
/// coefficient is invertible modulo `prime`.
pub(super) fn hensel_lift(
    coefficients: &[BigInt],
    leading_coefficient: &BigInt,
    prime: u64,
    factors_mod_prime: &[Fp],
    required_modulus: &BigInt,
) -> (Vec<Vec<BigInt>>, BigInt) {
    let prime_bigint = BigInt::from(prime);
    let mut modulus = prime_bigint.clone();
    let mut factors: Vec<Vec<BigInt>> = factors_mod_prime
        .iter()
        .map(|factor| factor.iter().map(|value| BigInt::from(*value)).collect())
        .collect();
    while modulus <= *required_modulus {
        let next_modulus = &modulus * &prime_bigint;
        let leading_inverse =
            super::arithmetic::mod_inverse_bigint(leading_coefficient, &next_modulus);
        let target: Vec<BigInt> = coefficients
            .iter()
            .map(|coefficient| mod_bigint(&(coefficient * &leading_inverse), &next_modulus))
            .collect();
        let product = product_bigint(&factors);
        let error: Fp = target
            .iter()
            .enumerate()
            .map(|(index, target_coefficient)| {
                let product_coefficient = product.get(index).cloned().unwrap_or_default();
                let quotient = (target_coefficient - product_coefficient) / &modulus;
                mod_bigint(&quotient, &prime_bigint)
                    .to_u64()
                    .expect("small prime residue")
            })
            .collect();
        for index in 0..factors.len() {
            let factor_mod_prime: Fp = factors[index]
                .iter()
                .map(|coefficient| {
                    mod_bigint(coefficient, &prime_bigint)
                        .to_u64()
                        .expect("small prime residue")
                })
                .collect();
            let cofactor = product_mod_prime_except(&factors, index, prime);
            let inverse = inverse_mod(&cofactor, &factor_mod_prime, prime);
            let correction = div_rem(&mul(&error, &inverse, prime), &factor_mod_prime, prime).1;
            for (degree, value) in correction.iter().copied().enumerate() {
                factors[index][degree] += &modulus * value;
            }
        }
        modulus = next_modulus;
    }
    for factor in &mut factors {
        for coefficient in factor {
            *coefficient = centered_mod_bigint(coefficient, &modulus);
        }
    }
    (factors, modulus)
}

/// Deterministic primes in ascending order, as used by Zassenhaus selection.
pub(super) fn candidate_primes() -> impl Iterator<Item = u64> {
    (2u64..).filter(|candidate| is_prime(*candidate))
}

fn is_prime(candidate: u64) -> bool {
    if candidate < 2 {
        return false;
    }
    if candidate == 2 {
        return true;
    }
    if candidate % 2 == 0 {
        return false;
    }
    let mut divisor = 3u64;
    while divisor <= candidate / divisor {
        if candidate % divisor == 0 {
            return false;
        }
        divisor += 2;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn product_mod_prime(factors: &[Fp], prime: u64) -> Fp {
        factors
            .iter()
            .fold(vec![1], |product, factor| mul(&product, factor, prime))
    }

    #[test]
    fn berlekamp_factors_a_squarefree_polynomial() {
        let polynomial = vec![1, 0, 0, 0, 1];
        let factors = factor_mod_prime(
            &polynomial
                .iter()
                .map(|value| BigInt::from(*value))
                .collect::<Vec<_>>(),
            3,
        )
        .expect("x^4+1 is squarefree modulo 3");
        assert!(factors.iter().all(|factor| factor.last() == Some(&1)));
        assert_eq!(product_mod_prime(&factors, 3), polynomial);
    }

    #[test]
    fn multifactor_hensel_lifting_preserves_the_target() {
        let coefficients = vec![BigInt::from(-1), BigInt::from(0), BigInt::from(1)];
        let factors_mod_prime = factor_mod_prime(&coefficients, 3).unwrap();
        let (factors, modulus) = hensel_lift(
            &coefficients,
            &BigInt::one(),
            3,
            &factors_mod_prime,
            &BigInt::from(50),
        );
        let product = product_bigint(&factors);
        for (index, coefficient) in product.iter().enumerate() {
            assert_eq!(
                mod_bigint(coefficient, &modulus),
                mod_bigint(&coefficients[index], &modulus)
            );
        }
    }

    #[test]
    fn multifactor_hensel_lifting_handles_nonmonic_integer_input() {
        let coefficients = vec![
            BigInt::from(3),
            BigInt::from(11),
            BigInt::from(-11),
            BigInt::from(2),
            BigInt::from(9),
            BigInt::from(-10),
            BigInt::from(8),
        ];
        let factors_mod_prime = factor_mod_prime(&coefficients, 3).unwrap();
        let (factors, modulus) = hensel_lift(
            &coefficients,
            &BigInt::from(8),
            3,
            &factors_mod_prime,
            &BigInt::from(50),
        );
        let product = product_bigint(&factors);
        let inverse = super::super::arithmetic::mod_inverse_bigint(&BigInt::from(8), &modulus);
        for (index, coefficient) in product.iter().enumerate() {
            let target = &coefficients[index] * &inverse;
            assert_eq!(
                mod_bigint(coefficient, &modulus),
                mod_bigint(&target, &modulus)
            );
        }
    }
}
