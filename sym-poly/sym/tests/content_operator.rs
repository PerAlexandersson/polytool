use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};
use sym_poly_core::Partition;
use sym_poly_sym::{
    c_n, c_n_inverse, partition_graph_symmetric, schur_content_product, Basis,
    ContentOperatorError, SymmetricFunction,
};

fn part(parts: &[u32]) -> Partition {
    Partition::new(parts.to_vec())
}

fn rational(numerator: i64, denominator: i64) -> BigRational {
    BigRational::new(numerator.into(), denominator.into())
}

fn lift(value: &SymmetricFunction<BigInt>) -> SymmetricFunction<BigRational> {
    SymmetricFunction::from_terms(
        value.basis(),
        value
            .terms()
            .iter()
            .map(|(p, c)| (p.clone(), BigRational::from_integer(c.clone())))
            .collect(),
    )
}

#[test]
fn content_convention_and_small_products() {
    assert_eq!(schur_content_product(&part(&[3]), 3), 60.into());
    assert_eq!(schur_content_product(&part(&[2, 1]), 3), 24.into());
    assert_eq!(schur_content_product(&part(&[1, 1, 1]), 3), 6.into());
    for degree in 0..=8 {
        for lambda in Partition::all_of_size(degree) {
            for n in 0..=10 {
                // Independent signed arithmetic, including boxes after a zero factor.
                let expected: BigInt = lambda
                    .diagram_boxes()
                    .iter()
                    .map(|&(i, j)| BigInt::from(n as i64 + j as i64 - i as i64))
                    .product();
                let actual = schur_content_product(&lambda, n);
                assert_eq!(actual, expected, "lambda={lambda}, n={n}");
                assert_eq!(actual.is_zero(), lambda.num_parts() > n);
            }
        }
    }
}

#[test]
fn large_products_and_full_usize_variable_count() {
    let n = BigInt::from(usize::MAX);
    assert_eq!(
        schur_content_product(&part(&[2]), usize::MAX),
        &n * (&n + 1)
    );
    let lambda = part(&[40]);
    let factor: BigInt = (100..140).map(BigInt::from).product();
    assert!(factor > BigInt::from(i64::MAX));
    assert_eq!(schur_content_product(&lambda, 100), factor);
    let s = SymmetricFunction::<BigInt>::schur_symmetric(lambda.clone());
    assert_eq!(c_n(&s, 100).coefficient(&lambda), factor);
    let inverse = c_n_inverse(&s, 100).unwrap();
    assert_eq!(
        inverse.coefficient(&lambda),
        BigRational::new(BigInt::one(), factor)
    );
    assert_eq!(c_n(&inverse, 100), lift(&s));
}

#[test]
fn zero_empty_and_zero_variables() {
    for n in [0, 1, 5] {
        let zero = SymmetricFunction::<BigInt>::zero(Basis::Forgotten);
        assert!(c_n(&zero, n).is_zero());
        assert!(c_n_inverse(&zero, n).unwrap().is_zero());
        let constant = SymmetricFunction::scaled_basis_element(
            Basis::Monomial,
            Partition::empty(),
            BigInt::from(-7),
        );
        assert_eq!(c_n(&constant, n), constant.to_schur_basis());
        assert_eq!(
            c_n_inverse(&constant, n).unwrap(),
            lift(&constant).to_schur_basis()
        );
    }
    let s1 = SymmetricFunction::<BigInt>::schur_symmetric(part(&[1]));
    assert!(c_n(&s1, 0).is_zero());
    assert_eq!(
        c_n_inverse(&s1, 0),
        Err(ContentOperatorError {
            partition: part(&[1]),
            num_vars: 0
        })
    );
}

#[test]
fn inverse_domain_is_checked_after_schur_conversion() {
    // p_2 = s_2 - s_11: its indexing partition has one row, but it is not
    // in the image of C_1 as a formal symmetric function.
    let p2 = SymmetricFunction::<BigInt>::power_sum_symmetric(part(&[2]));
    let error = c_n_inverse(&p2, 1).unwrap_err();
    assert_eq!(
        error,
        ContentOperatorError {
            partition: part(&[1, 1]),
            num_vars: 1
        }
    );
    assert!(error.to_string().contains("zero eigenvalue"));
    let projected = SymmetricFunction::<BigInt>::schur_symmetric(part(&[2]));
    assert_eq!(c_n_inverse(&c_n(&p2, 1), 1).unwrap(), lift(&projected));
    // At length = n, the eigenvalue is still nonzero.
    assert!(c_n_inverse(&p2, 2).is_ok());
}

#[test]
fn integer_input_has_fractional_inverse() {
    let s21 = SymmetricFunction::<BigInt>::schur_symmetric(part(&[2, 1]));
    assert_eq!(
        c_n_inverse(&s21, 3).unwrap().coefficient(&part(&[2, 1])),
        rational(1, 24)
    );
    // Generic forward operator also works in the existing machine-integer API.
    let s21_i64 = SymmetricFunction::<i64>::schur_symmetric(part(&[2, 1]));
    assert_eq!(c_n(&s21_i64, 3).coefficient(&part(&[2, 1])), 24);
}

#[test]
fn all_bases_and_both_roundtrip_directions() {
    let schur = SymmetricFunction::from_terms(
        Basis::Schur,
        [
            (part(&[4]), rational(2, 3)),
            (part(&[2, 1, 1]), rational(-7, 5)),
        ]
        .into_iter()
        .collect(),
    );
    for basis in [
        Basis::Schur,
        Basis::Monomial,
        Basis::Elementary,
        Basis::CompleteH,
        Basis::PowerSum,
        Basis::Forgotten,
    ] {
        let value = schur.to_basis(basis);
        assert_eq!(c_n(&value, 4), c_n(&schur, 4));
        assert_eq!(c_n_inverse(&c_n(&value, 4), 4).unwrap(), schur);
        assert_eq!(c_n(&c_n_inverse(&value, 4).unwrap(), 4), schur);
    }
}

#[test]
fn nonhomogeneous_linearity() {
    let a =
        SymmetricFunction::scaled_basis_element(Basis::Schur, Partition::empty(), rational(3, 7));
    let b = SymmetricFunction::scaled_basis_element(Basis::Schur, part(&[2]), rational(-4, 9));
    let sum = a.clone() + b.clone();
    assert_eq!(c_n(&sum, 3), c_n(&a, 3) + c_n(&b, 3));
    assert_eq!(
        c_n_inverse(&sum, 3).unwrap(),
        c_n_inverse(&a, 3).unwrap() + c_n_inverse(&b, 3).unwrap()
    );
    assert_eq!(c_n(&c_n_inverse(&sum, 3).unwrap(), 3), sum);
}

#[test]
fn graph_basis_integration() {
    for degree in 2..=7 {
        let n = 2 * degree as usize;
        for lambda in Partition::all_of_size(degree)
            .into_iter()
            .filter(|p| p.parts().iter().all(|&a| a >= 2))
        {
            let graph = partition_graph_symmetric(&lambda, n).unwrap();
            let inverse = c_n_inverse(&graph, n).unwrap();
            assert_eq!(c_n(&inverse, n), lift(&graph).to_schur_basis());
            assert_eq!(
                c_n_inverse(&c_n(&graph, n), n).unwrap(),
                lift(&graph).to_schur_basis()
            );
        }
    }
    // b_(2)^(4) = 6 s_2 - 10 s_11; eigenvalues 20 and 12.
    let graph = partition_graph_symmetric(&part(&[2]), 4).unwrap();
    let inverse = c_n_inverse(&graph, 4).unwrap();
    assert_eq!(inverse.coefficient(&part(&[2])), rational(3, 10));
    assert_eq!(inverse.coefficient(&part(&[1, 1])), rational(-5, 6));
}
