//! Exact arithmetic behind the Veronese-section lifts: coefficient sections,
//! bivariate polynomials in the shift variable `L`, the Veronese norm and
//! monic division in `L`.
//!
//! A bivariate polynomial is a [`Polynomial`] in `L` whose coefficients are
//! polynomials in a second variable (`x` for a parent, `y = x^r` for a
//! section); the public coefficient layout `[i][j]` is the coefficient of
//! `L^i` times the `j`-th power of the second variable.

use super::*;
use crate::polynomial::Polynomial;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

/// Exact univariate polynomial with rational coefficients.
pub(super) type UnivariatePolynomial = Polynomial<BigRational>;
/// Polynomial in `L` with coefficients in `Q[x]` (or `Q[y]`).
pub(super) type BivariatePolynomial = Polynomial<UnivariatePolynomial>;

/// Build a bivariate polynomial from `[i][j]` = coefficient of `L^i v^j`.
pub(super) fn bivariate_from_coefficients(rows: &[Vec<BigRational>]) -> BivariatePolynomial {
    Polynomial::new(
        rows.iter()
            .map(|row| Polynomial::new(row.clone()))
            .collect(),
    )
}

/// The `[i][j]` coefficient layout, each row trimmed of trailing zeros (a
/// zero row is empty).
pub(super) fn bivariate_to_coefficients(polynomial: &BivariatePolynomial) -> Vec<Vec<BigRational>> {
    polynomial
        .coeffs()
        .iter()
        .map(|row| row.coeffs().to_vec())
        .collect()
}

/// Degree in `L` (0 for the zero polynomial).
pub(super) fn degree_in_shift(polynomial: &BivariatePolynomial) -> usize {
    polynomial.degree().unwrap_or(0)
}

/// Largest degree in the second variable over all coefficients of `L^i`.
pub(super) fn degree_in_second_variable(polynomial: &BivariatePolynomial) -> usize {
    polynomial
        .coeffs()
        .iter()
        .filter_map(Polynomial::degree)
        .max()
        .unwrap_or(0)
}

/// Whether the leading coefficient in `L` is the constant 1.
pub(super) fn is_monic_in_shift(polynomial: &BivariatePolynomial) -> bool {
    polynomial
        .leading_coefficient()
        .is_some_and(|leading| leading == Polynomial::one())
}

/// `p(v^power)` for a univariate `p(v)`.
pub(super) fn substitute_power_univariate(
    polynomial: &UnivariatePolynomial,
    power: usize,
) -> UnivariatePolynomial {
    if polynomial.is_zero() {
        return Polynomial::zero();
    }
    let degree = polynomial.degree().unwrap_or(0);
    let mut coefficients = vec![BigRational::zero(); degree * power + 1];
    for (index, coefficient) in polynomial.coeffs().iter().enumerate() {
        coefficients[index * power] = coefficient.clone();
    }
    Polynomial::new(coefficients)
}

/// Substitute `v^power` for the second variable in every coefficient of `L^i`.
pub(super) fn substitute_power_bivariate(
    polynomial: &BivariatePolynomial,
    power: usize,
) -> BivariatePolynomial {
    Polynomial::new(
        polynomial
            .coeffs()
            .iter()
            .map(|coefficient| substitute_power_univariate(coefficient, power))
            .collect(),
    )
}

/// The coefficient section `V_{r,j}(p)(y) = sum_k [x^(r k + j)] p(x) y^k`.
pub(super) fn coefficient_section(
    coefficients: &[BigRational],
    modulus: usize,
    residue: usize,
) -> UnivariatePolynomial {
    Polynomial::new(
        coefficients
            .iter()
            .skip(residue)
            .step_by(modulus)
            .cloned()
            .collect(),
    )
}

/// `sum_a x^a s_a(x^r)` from the sections `s_0, ..., s_(r-1)`.
pub(super) fn assemble_from_sections(sections: &[UnivariatePolynomial]) -> UnivariatePolynomial {
    let modulus = sections.len();
    let length = sections
        .iter()
        .enumerate()
        .filter_map(|(residue, section)| section.degree().map(|d| d * modulus + residue + 1))
        .max()
        .unwrap_or(0);
    let mut coefficients = vec![BigRational::zero(); length];
    for (residue, section) in sections.iter().enumerate() {
        for (index, coefficient) in section.coeffs().iter().enumerate() {
            coefficients[index * modulus + residue] = coefficient.clone();
        }
    }
    Polynomial::new(coefficients)
}

/// The Veronese norm `prod_{w^r = 1} G(L, w x)` written in `y = x^r`.
///
/// Computed as the determinant of the `r x r` matrix of multiplication by `G`
/// on the basis `1, x, ..., x^(r-1)` of `Q[L, y][x] / (x^r - y)`: the entry in
/// row `c` and column `b` is the coefficient of `x^c` in `G x^b`, namely
/// `G_(c-b)(L, y)` for `c >= b` and `y G_(c-b+r)(L, y)` otherwise, where
/// `G_a(L, y)` is the section of `G` at residue `a` (taken coefficientwise in
/// `L`).  Since `y` is transcendental, `x^r - y` is irreducible, and the
/// determinant of multiplication is the field norm, which is the product
/// over the `r` conjugates `x -> w x`.
pub(super) fn veronese_norm(parent: &BivariatePolynomial, modulus: usize) -> BivariatePolynomial {
    let sections: Vec<BivariatePolynomial> = (0..modulus)
        .map(|residue| {
            Polynomial::new(
                parent
                    .coeffs()
                    .iter()
                    .map(|coefficient| coefficient_section(coefficient.coeffs(), modulus, residue))
                    .collect(),
            )
        })
        .collect();
    let y = Polynomial::constant(Polynomial::variable());
    let matrix: Vec<Vec<BivariatePolynomial>> = (0..modulus)
        .map(|row| {
            (0..modulus)
                .map(|column| {
                    if row >= column {
                        sections[row - column].clone()
                    } else {
                        y.clone() * sections[row + modulus - column].clone()
                    }
                })
                .collect()
        })
        .collect();
    let rows: Vec<usize> = (0..modulus).collect();
    laplace_determinant(&matrix, &rows, 0)
}

/// Determinant of the minor with the given rows and the columns `column..`,
/// by expansion along its first column (the matrices here are at most 6 x 6).
fn laplace_determinant(
    matrix: &[Vec<BivariatePolynomial>],
    rows: &[usize],
    column: usize,
) -> BivariatePolynomial {
    if rows.is_empty() {
        return Polynomial::one();
    }
    let mut determinant = BivariatePolynomial::zero();
    for (position, &row) in rows.iter().enumerate() {
        let entry = &matrix[row][column];
        if entry.is_zero() {
            continue;
        }
        let remaining: Vec<usize> = rows
            .iter()
            .enumerate()
            .filter(|&(other, _)| other != position)
            .map(|(_, &index)| index)
            .collect();
        let term = entry.clone() * laplace_determinant(matrix, &remaining, column + 1);
        determinant = if position % 2 == 0 {
            determinant + term
        } else {
            determinant - term
        };
    }
    determinant
}

/// Division in `Q[v][L]` by a divisor monic in `L`: `(quotient, remainder)`.
///
/// Returns `None` when the divisor is not monic in `L`.
pub(super) fn divide_monic_in_shift(
    dividend: &BivariatePolynomial,
    divisor: &BivariatePolynomial,
) -> Option<(BivariatePolynomial, BivariatePolynomial)> {
    if !is_monic_in_shift(divisor) {
        return None;
    }
    let divisor_degree = degree_in_shift(divisor);
    let mut remainder: Vec<UnivariatePolynomial> = dividend.coeffs().to_vec();
    if remainder.len() <= divisor_degree {
        return Some((Polynomial::zero(), dividend.clone()));
    }
    let mut quotient = vec![UnivariatePolynomial::zero(); remainder.len() - divisor_degree];
    for top in (divisor_degree..remainder.len()).rev() {
        let factor = remainder[top].clone();
        if factor.is_zero() {
            continue;
        }
        let shift = top - divisor_degree;
        for (index, coefficient) in divisor.coeffs().iter().enumerate() {
            let product = factor.clone() * coefficient.clone();
            remainder[shift + index] = remainder[shift + index].clone() - product;
        }
        quotient[shift] = factor;
    }
    Some((Polynomial::new(quotient), Polynomial::new(remainder)))
}

/// Exact quotient `dividend / divisor` in `Q[v]`, or `None` when the division
/// leaves a remainder (or the divisor is zero).
pub(super) fn exact_quotient(
    dividend: &UnivariatePolynomial,
    divisor: &UnivariatePolynomial,
) -> Option<UnivariatePolynomial> {
    if divisor.is_zero() {
        return None;
    }
    let (quotient, remainder) = dividend.div_rem(divisor);
    remainder.is_zero().then_some(quotient)
}

/// Exact sign-and-roots summary of one polynomial row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct RowShape {
    /// The row is the zero polynomial.
    pub(super) is_zero: bool,
    /// All nonzero coefficients have the same sign.
    pub(super) has_one_sign: bool,
    /// The row is nonzero and has only real roots.
    pub(super) is_real_rooted: bool,
}

impl RowShape {
    /// Nonzero, coefficients of one sign and real-rooted: all roots in `(-inf, 0]`.
    pub(super) fn has_nonpositive_real_roots(self) -> bool {
        !self.is_zero && self.has_one_sign && self.is_real_rooted
    }
}

/// Classify a row exactly: zero, one coefficient sign, real-rootedness (after
/// clearing denominators, through [`crate::is_real_rooted_bigint_coeffs`]).
pub(super) fn classify_row(coefficients: &[BigRational]) -> RowShape {
    if coefficients.iter().all(Zero::is_zero) {
        return RowShape {
            is_zero: true,
            has_one_sign: true,
            is_real_rooted: false,
        };
    }
    let has_positive = coefficients.iter().any(Signed::is_positive);
    let has_negative = coefficients.iter().any(Signed::is_negative);
    let common_denominator = coefficients
        .iter()
        .fold(BigInt::one(), |accumulator, coefficient| {
            accumulator.lcm(coefficient.denom())
        });
    let integer_coefficients: Vec<BigInt> = coefficients
        .iter()
        .map(|coefficient| coefficient.numer() * (&common_denominator / coefficient.denom()))
        .collect();
    RowShape {
        is_zero: false,
        has_one_sign: !(has_positive && has_negative),
        is_real_rooted: crate::is_real_rooted_bigint_coeffs(&integer_coefficients),
    }
}

/// Plain-text form of a bivariate polynomial with `[i][j]` = coefficient of
/// `L^i v^j`, highest power of `L` first, for example `L^2 - L - x`.
pub(super) fn format_shift_polynomial(coefficients: &[Vec<BigRational>], variable: &str) -> String {
    let mut terms: Vec<String> = Vec::new();
    for (shift_power, row) in coefficients.iter().enumerate().rev() {
        for (variable_power, coefficient) in row.iter().enumerate() {
            if coefficient.is_zero() {
                continue;
            }
            let mut monomial = String::new();
            match shift_power {
                0 => {}
                1 => monomial.push('L'),
                power => monomial.push_str(&format!("L^{power}")),
            }
            match variable_power {
                0 => {}
                1 => monomial.push_str(variable),
                power => monomial.push_str(&format!("{variable}^{power}")),
            }
            let text = if monomial.is_empty() {
                format_rational_coeff(coefficient)
            } else if coefficient.is_one() {
                monomial
            } else if *coefficient == -BigRational::one() {
                format!("-{monomial}")
            } else {
                format!("{}{monomial}", format_rational_coeff(coefficient))
            };
            terms.push(text);
        }
    }
    let Some((first, rest)) = terms.split_first() else {
        return "0".into();
    };
    let mut result = first.clone();
    for term in rest {
        if let Some(positive) = term.strip_prefix('-') {
            result.push_str(" - ");
            result.push_str(positive);
        } else {
            result.push_str(" + ");
            result.push_str(term);
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(value: i64) -> BigRational {
        BigRational::from_integer(BigInt::from(value))
    }

    fn rows(values: &[&[i64]]) -> Vec<Vec<BigRational>> {
        values
            .iter()
            .map(|row| row.iter().map(|&value| q(value)).collect())
            .collect()
    }

    #[test]
    fn norm_of_linear_parent_is_product_of_conjugates() {
        // G = L - 1 - x: the r = 2 norm is (L - 1)^2 - y.
        let parent = bivariate_from_coefficients(&rows(&[&[-1, -1], &[1]]));
        let norm = veronese_norm(&parent, 2);
        assert_eq!(
            bivariate_to_coefficients(&norm),
            rows(&[&[1, -1], &[-2], &[1]])
        );
        // r = 3: (L - 1)^3 - y.
        let norm = veronese_norm(&parent, 3);
        assert_eq!(
            bivariate_to_coefficients(&norm),
            rows(&[&[-1, -1], &[3], &[-3], &[1]])
        );
        assert_eq!(
            format_shift_polynomial(&rows(&[&[-1, -1], &[1]]), "x"),
            "L - 1 - x"
        );
    }

    #[test]
    fn norm_of_fibonacci_parent() {
        // G = L^2 - L - x: G(L, x) G(L, -x) = (L^2 - L)^2 - y.
        let parent = bivariate_from_coefficients(&rows(&[&[0, -1], &[-1], &[1]]));
        let norm = veronese_norm(&parent, 2);
        assert_eq!(
            bivariate_to_coefficients(&norm),
            rows(&[&[0, -1], &[], &[1], &[-2], &[1]])
        );
        assert_eq!(
            format_shift_polynomial(&bivariate_to_coefficients(&norm), "y"),
            "L^4 - 2L^3 + L^2 - y"
        );
    }

    #[test]
    fn monic_division_and_sections() {
        let dividend = bivariate_from_coefficients(&rows(&[&[1, -1], &[-2], &[1]]));
        let divisor = bivariate_from_coefficients(&rows(&[&[-1], &[1]]));
        let (quotient, remainder) = divide_monic_in_shift(&dividend, &divisor).unwrap();
        assert_eq!(bivariate_to_coefficients(&quotient), rows(&[&[-1], &[1]]));
        assert_eq!(bivariate_to_coefficients(&remainder), rows(&[&[0, -1]]));
        let polynomial = rows(&[&[1, 2, 3, 4, 5]]).remove(0);
        assert_eq!(
            coefficient_section(&polynomial, 2, 1).coeffs(),
            &[q(2), q(4)]
        );
        let sections = vec![
            coefficient_section(&polynomial, 3, 0),
            coefficient_section(&polynomial, 3, 1),
            coefficient_section(&polynomial, 3, 2),
        ];
        assert_eq!(assemble_from_sections(&sections).coeffs(), &polynomial[..]);
    }

    #[test]
    fn row_classification() {
        assert!(classify_row(&rows(&[&[1, 2, 1]])[0]).has_nonpositive_real_roots());
        assert!(classify_row(&rows(&[&[3]])[0]).has_nonpositive_real_roots());
        assert!(classify_row(&rows(&[&[-1, -1]])[0]).has_nonpositive_real_roots());
        let mixed = classify_row(&rows(&[&[1, -1]])[0]);
        assert!(!mixed.has_one_sign && mixed.is_real_rooted);
        let complex = classify_row(&rows(&[&[1, 1, 1]])[0]);
        assert!(complex.has_one_sign && !complex.is_real_rooted);
        assert!(classify_row(&rows(&[&[0]])[0]).is_zero);
        let half = vec![BigRational::new(BigInt::from(1), BigInt::from(2)), q(1)];
        assert!(classify_row(&half).has_nonpositive_real_roots());
    }
}
