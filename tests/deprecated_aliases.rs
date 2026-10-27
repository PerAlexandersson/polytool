//! The deprecated pre-0.6.0 names still compile and agree with their replacements.
#![allow(deprecated)]

use num_bigint::BigInt;
use polytool::sequences::deco;
use polytool::*;

fn big(values: &[i64]) -> Vec<BigInt> {
    values.iter().copied().map(BigInt::from).collect()
}

#[test]
fn renamed_real_rootedness_helpers_forward() {
    assert_eq!(real_roots(&[2, -3, 1]), approximate_real_roots(&[2, -3, 1]));
    assert_eq!(
        ehrhart_to_hstar_with_denom(&[2, 3, 1], 2),
        real_rootedness::ehrhart_to_hstar_with_denominator(&[2, 3, 1], 2)
    );
    let coefficients: Vec<BigRational> = [0, -2, 1]
        .iter()
        .map(|&c| BigRational::from_integer(BigInt::from(c)))
        .collect();
    assert_eq!(
        real_rootedness::isolate_real_roots_exact(&coefficients),
        real_rootedness::isolate_real_roots_rational_coeffs(&coefficients)
    );
}

#[test]
fn renamed_root_count_filters_forward() {
    let polynomial = big(&[1, 5, 6, 1]);
    assert_eq!(
        satisfies_newton_inequalities_bigint(&polynomial),
        satisfies_newton_inequalities_bigint_coeffs(&polynomial)
    );
    assert_eq!(
        root_count::satisfies_kurtz_condition_bigint(&polynomial),
        satisfies_kurtz_condition_bigint_coeffs(&polynomial)
    );
}

#[test]
fn renamed_total_nonnegativity_helpers_forward() {
    let pascal = vec![vec![1, 0], vec![1, 1], vec![1, 2]];
    assert!(is_tnn(&pascal));
    assert_eq!(
        is_totally_nonnegative(&pascal, 2),
        is_totally_nonnegative_by_minors(&pascal, 2)
    );
    let matrix = vec![vec![vec![1, 1]]];
    assert_eq!(
        is_lace_totally_nonnegative_i64(&matrix, 4, 5, 4),
        is_lace_totally_nonnegative_by_minors_i64(&matrix, 4, 5, 4)
    );
    let polynomials = vec![vec![2, 1], vec![8, 6, 1], vec![3, 4, 1]];
    assert_eq!(
        is_lace_sequence_totally_nonnegative_i64(&polynomials, 1, 3, 3),
        is_lace_sequence_totally_nonnegative_by_minors_i64(&polynomials, 1, 3, 3)
    );
    assert!(check_lace_total_nonnegative_i64(&matrix, 4, 5, 4).is_ok());
    assert!(check_lace_sequence_total_nonnegative_i64(&polynomials, 1, 3, 3).is_err());
}

#[test]
fn renamed_sequence_and_sieving_helpers_forward() {
    assert_eq!(deco::polynomials(4, 1), deco::deco_polynomials(4, 1));
    assert_eq!(
        deco::polynomials_bigint(4, &BigInt::from(1)),
        deco::deco_polynomials_bigint(4, &BigInt::from(1))
    );
    let polynomial = big(&[1, 1, 1]);
    let message = cyclic_sieving_check_orbit_sizes(&polynomial, 3, &[3]).map(|check| check.holds);
    assert_eq!(
        message,
        check_cyclic_sieving_orbit_sizes(&polynomial, 3, &[3])
            .map(|check| check.holds)
            .map_err(|error| error.to_string())
    );
    assert!(matches!(
        cyclic_sieving_check_orbit_sizes(&polynomial, 0, &[]),
        Err(message) if message == "cyclic group order must be positive"
    ));
}

#[test]
fn deprecated_decomposition_alias_is_the_generic_type() {
    let analysis: SymmetricDecompositionAnalysisBigInt =
        analyze_symmetric_decomposition_bigint(&big(&[1, 4, 1])).unwrap();
    let generic: SymmetricDecompositionAnalysis<BigInt> = analysis.clone();
    assert_eq!(analysis, generic);
    assert_eq!(
        analyze_symmetric_decomposition_i64(&[1, 4, 1])
            .unwrap()
            .degree,
        2
    );
}
