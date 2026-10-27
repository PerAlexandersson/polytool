//! Chain polynomials of lower-unitriangular matrices (Brändén--Saud Maia
//! Leite, arXiv:2412.06595, Theorem 3.7) and structured TNN diagnostics.

use num_bigint::BigInt;
use polytool::numeric::binomial_bigint;
use polytool::{
    check_chain_polynomials_bigint, check_lower_unitriangular_totally_nonnegative,
    evaluate_path_matrix, is_real_rooted_bigint_coeffs, lower_triangular_rows_from_square,
    reconstruct_canonical_totally_nonnegative_network, BigRational, CanonicalPlanarNetwork,
    NetworkError,
};

fn bi(v: i64) -> BigInt {
    BigInt::from(v)
}

fn rows_i64(rows: &[&[i64]]) -> Vec<Vec<BigInt>> {
    rows.iter()
        .map(|r| r.iter().map(|&v| bi(v)).collect())
        .collect()
}

fn pascal(size: usize) -> Vec<Vec<BigInt>> {
    (0..size)
        .map(|n| (0..=n).map(|k| binomial_bigint(n, k)).collect())
        .collect()
}

fn stirling2(size: usize) -> Vec<Vec<BigInt>> {
    let mut rows: Vec<Vec<BigInt>> = vec![vec![bi(1)]];
    for n in 1..size {
        let prev = rows[n - 1].clone();
        rows.push(
            (0..=n)
                .map(|k| {
                    let stay = prev.get(k).cloned().unwrap_or_default() * bi(k as i64);
                    let up = if k > 0 { prev[k - 1].clone() } else { bi(0) };
                    stay + up
                })
                .collect(),
        );
    }
    rows
}

/// Path matrix of the canonical network with integer weights `lambda`.
fn network_matrix(lambda: Vec<Vec<i64>>) -> Vec<Vec<BigInt>> {
    let network = CanonicalPlanarNetwork {
        lambda: lambda
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|v| BigRational::from_integer(bi(v)))
                    .collect()
            })
            .collect(),
    };
    evaluate_path_matrix(&network)
        .into_iter()
        .map(|row| row.into_iter().map(|q| q.to_integer()).collect())
        .collect()
}

fn assert_theorem_holds(rows: &[Vec<BigInt>]) {
    let report = check_chain_polynomials_bigint(rows).unwrap();
    assert!(
        report.is_totally_nonnegative(),
        "{:?}",
        report.totally_nonnegative_verdict
    );
    assert!(report.chain_conclusion_holds(), "{:?}", report);
    assert!(!report.contradicts_theorem());
    for (n, p) in report.chain_polynomials.iter().enumerate() {
        assert!(p.len() <= n + 1);
    }
}

#[test]
fn pascal_chain_polynomials_are_fubini_polynomials() {
    let report = check_chain_polynomials_bigint(&pascal(7)).unwrap();
    // p_n = sum_k k! S(n,k) t^k.
    assert_eq!(
        report.chain_polynomials[5],
        rows_i64(&[&[0, 1, 30, 150, 240, 120]])[0]
    );
    assert_theorem_holds(&pascal(9));
}

#[test]
fn tnn_matrices_give_real_rooted_interlacing_chain_polynomials() {
    assert_theorem_holds(&stirling2(8));
    assert_theorem_holds(&network_matrix(vec![
        vec![2],
        vec![1, 3],
        vec![0, 1, 2],
        vec![0, 4, 0, 1],
        vec![0, 1, 1, 5, 2],
    ]));
    assert_theorem_holds(&network_matrix(vec![
        vec![1],
        vec![3, 1],
        vec![2, 2, 7],
        vec![1, 0, 0, 3],
    ]));
    // Lower-triangular Toeplitz matrix of the PF sequence (1 + t)^2 = 1, 2, 1.
    assert_theorem_holds(&rows_i64(&[
        &[1],
        &[2, 1],
        &[1, 2, 1],
        &[0, 1, 2, 1],
        &[0, 0, 1, 2, 1],
    ]));
}

#[test]
fn all_ones_triangle_is_tnn_but_its_last_row_is_not_real_rooted() {
    let rows = rows_i64(&[&[1], &[1, 1], &[1, 1, 1]]);
    let report = check_chain_polynomials_bigint(&rows).unwrap();
    assert!(report.is_totally_nonnegative());
    assert_eq!(
        report.chain_polynomials,
        rows_i64(&[&[1], &[0, 1], &[0, 1, 1]])
    );
    assert!(report.chain_conclusion_holds());
    // The guardrail: the raw row 1 + t + t^2 is not real-rooted.
    assert!(!is_real_rooted_bigint_coeffs(&rows[2]));
}

#[test]
fn non_tnn_matrices_report_a_negative_minor() {
    // Lower-triangular Toeplitz matrix of 1 + t + t^2, which is not PF.
    let rows = rows_i64(&[&[1], &[1, 1], &[1, 1, 1], &[0, 1, 1, 1]]);
    let report = check_chain_polynomials_bigint(&rows).unwrap();
    match &report.totally_nonnegative_verdict {
        Err(NetworkError::NotTotallyNonnegative {
            witness: Some(witness),
            ..
        }) => {
            assert!(witness.verify(&rows));
            assert_eq!(witness.value, bi(-1));
        }
        other => panic!("expected a witness, got {:?}", other),
    }
    assert!(!report.contradicts_theorem());
    // The canonical-network reconstruction reports the same diagnosis.
    assert!(matches!(
        reconstruct_canonical_totally_nonnegative_network(&rows),
        Err(NetworkError::NotTotallyNonnegative {
            witness: Some(_),
            ..
        })
    ));
}

#[test]
fn chain_report_keeps_tnn_and_shape_errors_apart() {
    // Negative entry: the chain polynomials are still generated.
    let rows = rows_i64(&[&[1], &[-1, 1]]);
    let report = check_chain_polynomials_bigint(&rows).unwrap();
    assert_eq!(report.chain_polynomials[1], vec![bi(0), bi(-1)]);
    assert!(matches!(
        report.totally_nonnegative_verdict,
        Err(NetworkError::NegativeEntry { row: 1, col: 0, .. })
    ));
    // Without TNN the conclusion can fail: p_2 = t (1 - 3t) has the zero 1/3.
    let rows = rows_i64(&[&[1], &[1, 1], &[1, -3, 1]]);
    let report = check_chain_polynomials_bigint(&rows).unwrap();
    assert_eq!(report.chain_polynomials[2], vec![bi(0), bi(1), bi(-3)]);
    assert!(report.real_rooted[2]);
    assert!(!report.roots_in_minus_one_zero[2]);
    assert!(!report.is_totally_nonnegative() && !report.chain_conclusion_holds());
    // Shape errors are returned directly.
    assert_eq!(
        check_chain_polynomials_bigint(&rows_i64(&[&[1], &[1, 2]])),
        Err(NetworkError::NonUnitDiagonal {
            row: 1,
            value: bi(2)
        })
    );
    assert!(matches!(
        check_chain_polynomials_bigint(&rows_i64(&[&[1], &[1]])),
        Err(NetworkError::InvalidRowLength { row: 1, .. })
    ));
}

#[test]
fn square_input_is_converted_with_diagnostics() {
    let square = rows_i64(&[&[1, 0, 0], &[2, 1, 0], &[1, 2, 1]]);
    let rows = lower_triangular_rows_from_square(&square).unwrap();
    assert!(check_lower_unitriangular_totally_nonnegative(&rows).is_ok());
    let upper = rows_i64(&[&[1, 0, 0], &[2, 1, 3], &[1, 2, 1]]);
    assert_eq!(
        lower_triangular_rows_from_square(&upper),
        Err(NetworkError::NonzeroAboveDiagonal {
            row: 1,
            col: 2,
            value: bi(3)
        })
    );
}
