//! Integration tests for exact rational polynomial factorization.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};
use polytool::factorization::{
    irreducible_factorization_bigint_coeffs, irreducible_factorization_monic_bivariate,
    irreducible_factorization_rational_coeffs,
    irreducible_factorization_rational_coeffs_with_budget, MonicBivariateFactor,
    PolynomialFactorizationError,
};

type Q = BigRational;

fn q(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

fn trim(mut polynomial: Vec<Q>) -> Vec<Q> {
    while polynomial.last().is_some_and(Zero::is_zero) {
        polynomial.pop();
    }
    polynomial
}

fn mul(a: &[Q], b: &[Q]) -> Vec<Q> {
    let mut result = vec![Q::zero(); a.len() + b.len() - 1];
    for (i, left) in a.iter().enumerate() {
        for (j, right) in b.iter().enumerate() {
            result[i + j] += left * right;
        }
    }
    trim(result)
}

fn reconstruct_univariate(constant: &Q, factors: &[(Vec<BigInt>, usize)]) -> Vec<Q> {
    let mut result = vec![constant.clone()];
    for (factor, multiplicity) in factors {
        let factor = factor
            .iter()
            .cloned()
            .map(Q::from_integer)
            .collect::<Vec<_>>();
        for _ in 0..*multiplicity {
            result = mul(&result, &factor);
        }
    }
    result
}

fn bivariate_mul(a: &[Vec<Q>], b: &[Vec<Q>]) -> Vec<Vec<Q>> {
    let mut result = vec![Vec::new(); a.len() + b.len() - 1];
    for (i, left) in a.iter().enumerate() {
        for (j, right) in b.iter().enumerate() {
            let product = mul(left, right);
            if result[i + j].len() < product.len() {
                result[i + j].resize(product.len(), Q::zero());
            }
            for (degree, coefficient) in product.iter().enumerate() {
                result[i + j][degree] += coefficient;
            }
        }
    }
    result.into_iter().map(trim).collect()
}

fn reconstruct_bivariate(factors: &[MonicBivariateFactor]) -> Vec<Vec<Q>> {
    let mut result = vec![vec![Q::one()]];
    for factor in factors {
        for _ in 0..factor.multiplicity {
            result = bivariate_mul(&result, &factor.coefficients);
        }
    }
    while result.last().is_some_and(Vec::is_empty) {
        result.pop();
    }
    result
}

fn assert_univariate_reconstructs(coefficients: &[i64]) {
    let input = coefficients.iter().copied().map(q).collect::<Vec<_>>();
    let factorization = irreducible_factorization_rational_coeffs(&input).unwrap();
    assert_eq!(
        reconstruct_univariate(&factorization.constant, &factorization.factors),
        input
    );
}

#[test]
fn factors_repeated_content_and_rational_coefficients_exactly() {
    let repeated = irreducible_factorization_bigint_coeffs(&[
        BigInt::from(-1),
        BigInt::from(3),
        BigInt::from(-3),
        BigInt::from(1),
    ])
    .unwrap();
    assert_eq!(
        repeated.factors,
        vec![(vec![BigInt::from(-1), BigInt::one()], 3)]
    );
    assert_eq!(repeated.constant, Q::one());

    let rational = vec![
        q(1),
        q(0),
        Q::new(BigInt::from(-3), BigInt::from(2)),
        q(1) / q(2),
    ];
    let factorization = irreducible_factorization_rational_coeffs(&rational).unwrap();
    assert_eq!(
        reconstruct_univariate(&factorization.constant, &factorization.factors),
        rational
    );
    assert_eq!(factorization.factors.len(), 2);
}

#[test]
fn cyclotomic_products_and_special_irreducibles() {
    for degree in 1..=30 {
        let mut polynomial = vec![BigInt::zero(); degree + 1];
        polynomial[0] = -BigInt::one();
        polynomial[degree] = BigInt::one();
        let factorization = irreducible_factorization_bigint_coeffs(&polynomial).unwrap();
        let reconstructed = reconstruct_univariate(&factorization.constant, &factorization.factors);
        let expected = polynomial
            .into_iter()
            .map(Q::from_integer)
            .collect::<Vec<_>>();
        assert_eq!(reconstructed, expected, "x^{degree} - 1");
    }
    let cyclotomic = irreducible_factorization_bigint_coeffs(
        &[1i64, 0, 0, 0, 1]
            .iter()
            .copied()
            .map(BigInt::from)
            .collect::<Vec<_>>(),
    )
    .unwrap();
    assert_eq!(cyclotomic.factors.len(), 1);
    assert_univariate_reconstructs(&[-16, 0, 0, 0, 0, 0, 0, 0, 1]);
    let swinnerton_dyer = irreducible_factorization_bigint_coeffs(
        &[576i64, 0, -960, 0, 352, 0, -40, 0, 1]
            .iter()
            .copied()
            .map(BigInt::from)
            .collect::<Vec<_>>(),
    )
    .unwrap();
    assert_eq!(swinnerton_dyer.factors.len(), 1);
}

#[test]
fn bivariate_examples_reconstruct_and_sort() {
    let examples = [
        (vec![vec![q(1), q(0), q(-1)], vec![q(-2)], vec![q(1)]], 2),
        (
            vec![
                vec![q(1), q(0), q(0), q(0), q(-1)],
                vec![q(-4)],
                vec![q(6)],
                vec![q(-4)],
                vec![q(1)],
            ],
            3,
        ),
        (
            vec![
                vec![q(-1), q(0), q(0), q(-1)],
                vec![q(3)],
                vec![q(-3)],
                vec![q(1)],
            ],
            2,
        ),
        (
            vec![
                vec![q(0), q(0), q(-1)],
                vec![],
                vec![q(1)],
                vec![q(-2)],
                vec![q(1)],
            ],
            2,
        ),
    ];
    for (input, expected_count) in examples {
        let factors = irreducible_factorization_monic_bivariate(&input).unwrap();
        assert_eq!(factors.len(), expected_count);
        assert_eq!(reconstruct_bivariate(&factors), input);
        assert!(factors.windows(2).all(|pair| {
            (pair[0].coefficients.len(), &pair[0].coefficients)
                <= (pair[1].coefficients.len(), &pair[1].coefficients)
        }));
    }

    let first = vec![vec![q(-1), q(-3), q(-2)], vec![q(1)]];
    let second = vec![
        vec![q(1), q(-3), q(7), q(-6), q(4)],
        vec![q(-2), q(3), q(2)],
        vec![q(1)],
    ];
    let product = bivariate_mul(&first, &second);
    let factors = irreducible_factorization_monic_bivariate(&product).unwrap();
    assert_eq!(factors.len(), 2);
    assert_eq!(reconstruct_bivariate(&factors), product);
}

#[test]
fn bivariate_irreducible_and_errors() {
    let irreducible = vec![vec![q(0), q(-1)], vec![], vec![q(1)]];
    let factors = irreducible_factorization_monic_bivariate(&irreducible).unwrap();
    assert_eq!(factors.len(), 1);
    assert_eq!(reconstruct_bivariate(&factors), irreducible);

    let second_irreducible = vec![vec![q(2), q(0), q(-1)], vec![q(-2)], vec![q(1)]];
    let factors = irreducible_factorization_monic_bivariate(&second_irreducible).unwrap();
    assert_eq!(factors.len(), 1);
    assert_eq!(reconstruct_bivariate(&factors), second_irreducible);

    assert_eq!(
        irreducible_factorization_bigint_coeffs(&[]),
        Err(PolynomialFactorizationError::ZeroPolynomial)
    );
    assert_eq!(
        irreducible_factorization_monic_bivariate(&[vec![q(0)]]),
        Err(PolynomialFactorizationError::ZeroPolynomial)
    );
    assert_eq!(
        irreducible_factorization_monic_bivariate(&[vec![q(2)], vec![q(2)]]),
        Err(PolynomialFactorizationError::NotMonicInFirstVariable)
    );
    let repeated = vec![vec![q(0), q(0), q(1)], vec![q(0), q(-2)], vec![q(1)]];
    assert_eq!(
        irreducible_factorization_monic_bivariate(&repeated),
        Err(PolynomialFactorizationError::NotSquarefree)
    );
}

#[test]
fn recombination_budget_is_enforced() {
    let input = vec![q(-1), q(0), q(1)];
    assert_eq!(
        irreducible_factorization_rational_coeffs_with_budget(&input, 0),
        Err(PolynomialFactorizationError::RecombinationBudgetExceeded { limit: 0 })
    );
}

#[test]
fn large_leading_coefficients_recombine_completely() {
    let factors: [&[i64]; 4] = [&[-3, 7, 1000], &[-2, 0, 0, 999], &[5, 1], &[-1, 12345]];
    let mut product = vec![q(1)];
    for factor in factors {
        product = mul(&product, &factor.iter().copied().map(q).collect::<Vec<_>>());
    }
    let factorization = irreducible_factorization_rational_coeffs(&product).unwrap();
    assert_eq!(factorization.factors.len(), 4);
    assert_eq!(
        reconstruct_univariate(&factorization.constant, &factorization.factors),
        trim(product)
    );
}
