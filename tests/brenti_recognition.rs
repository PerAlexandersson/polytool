//! Recognition of polynomial rows as rows of a planar strip digraph.

use num_bigint::BigInt;
use num_rational::Ratio;
use polytool::numeric::binomial_bigint;
use polytool::sequences::{
    binomial_eulerian_polynomials_bigint, eulerian_polynomials_bigint, narayana_polynomials_bigint,
    type_b_eulerian_polynomials_bigint,
};
use polytool::{
    check_totally_nonnegative_neville_bigint, is_real_rooted_bigint_coeffs,
    recognize_brenti_strip_digraph, BigRational, BrentiEdge, BrentiRecognition,
    BrentiRecognitionError, BrentiRecognitionOptions, BrentiStripDigraph, BrentiWeightDomain,
    Polynomial,
};

fn bi(v: i64) -> BigInt {
    BigInt::from(v)
}

fn br(v: i64) -> BigRational {
    BigRational::from_integer(bi(v))
}

fn rows_i64(rows: &[&[i64]]) -> Vec<Vec<BigInt>> {
    rows.iter()
        .map(|r| r.iter().map(|&v| bi(v)).collect())
        .collect()
}

fn polys(rows: &[Vec<BigInt>]) -> Vec<Polynomial<BigInt>> {
    rows.iter().map(|r| Polynomial::new(r.clone())).collect()
}

fn recognize(rows: &[Vec<BigInt>], options: &BrentiRecognitionOptions) -> BrentiRecognition {
    recognize_brenti_strip_digraph(&polys(rows), options).unwrap()
}

fn rational() -> BrentiRecognitionOptions {
    BrentiRecognitionOptions::default()
}

fn integer() -> BrentiRecognitionOptions {
    BrentiRecognitionOptions {
        weights: BrentiWeightDomain::Integer,
        ..Default::default()
    }
}

fn certified_digraph(result: BrentiRecognition, rows: &[Vec<BigInt>]) -> BrentiStripDigraph {
    match result {
        BrentiRecognition::Certified {
            digraph,
            certificate,
        } => {
            assert!(digraph.validate().is_ok());
            assert_eq!(certificate.digraph, digraph);
            certificate
                .verify_against_bigint_polynomials(&polys(rows))
                .unwrap();
            digraph
        }
        other => panic!("expected a certificate, got {}", other),
    }
}

fn assert_impossible(result: &BrentiRecognition) -> usize {
    match result {
        BrentiRecognition::ImpossibleWithinBounds {
            supports_checked, ..
        } => *supports_checked,
        other => panic!("expected impossibility, got {}", other),
    }
}

fn edge(to: usize, weight: BigRational) -> BrentiEdge {
    BrentiEdge { to, weight }
}

fn pascal_rows(max_n: usize) -> Vec<Vec<BigInt>> {
    (1..=max_n)
        .map(|n| (0..=n).map(|k| binomial_bigint(n, k)).collect())
        .collect()
}

/// Rows `sum_k S(n,k) t^k` for `n = 1..=max_n` (Stirling numbers of the
/// second kind).
fn stirling2_rows(max_n: usize) -> Vec<Vec<BigInt>> {
    let mut rows = vec![vec![bi(1)]];
    for n in 1..=max_n {
        let prev = &rows[n - 1];
        let row: Vec<BigInt> = (0..=n)
            .map(|k| {
                let stay = prev.get(k).cloned().unwrap_or_default() * bi(k as i64);
                let up = if k > 0 {
                    prev.get(k - 1).cloned().unwrap_or_default()
                } else {
                    bi(0)
                };
                stay + up
            })
            .collect();
        rows.push(row);
    }
    rows.split_off(1)
}

#[test]
fn pascal_rows_give_the_binomial_strip() {
    let rows = pascal_rows(6);
    for options in [rational(), integer()] {
        let digraph = certified_digraph(recognize(&rows, &options), &rows);
        // Below the top height: y -> y and y -> y + 1 with weight 1.  The top
        // height is never left within the prefix, so it keeps no edges.
        let mut expected: Vec<Vec<BrentiEdge>> = (0..6)
            .map(|y| vec![edge(y, br(1)), edge(y + 1, br(1))])
            .collect();
        expected.push(vec![]);
        assert_eq!(digraph.outgoing, expected);
    }
}

#[test]
fn stirling_rows_give_weights_equal_to_the_height() {
    let rows = stirling2_rows(6);
    let digraph = certified_digraph(recognize(&rows, &integer()), &rows);
    // S(n+1,k) = S(n,k-1) + k S(n,k): edges y -> y (weight y) and y -> y + 1.
    assert_eq!(digraph.outgoing[0], vec![edge(1, br(1))]);
    for y in 1..6 {
        assert_eq!(
            digraph.outgoing[y],
            vec![edge(y, br(y as i64)), edge(y + 1, br(1))]
        );
    }
}

#[test]
fn standard_real_rooted_families_are_not_strip_digraph_rows() {
    // Each of these real-rooted families violates the forced transitions
    // already without planarity (a column of W would need a negative weight).
    let eulerian = eulerian_polynomials_bigint(7);
    let narayana = narayana_polynomials_bigint(7);
    let type_b = type_b_eulerian_polynomials_bigint(6);
    let binomial_eulerian = binomial_eulerian_polynomials_bigint(6);
    let first_zero = BrentiRecognitionOptions {
        first_row_index: 0,
        ..Default::default()
    };
    // A_1 = 1, A_2, ... read as rows 0, 1, ...; and A_2, A_3, ... as rows 1, 2, ...
    assert_eq!(
        assert_impossible(&recognize(&eulerian[1..], &first_zero)),
        0
    );
    assert_eq!(
        assert_impossible(&recognize(&eulerian[2..], &rational())),
        0
    );
    assert_eq!(
        assert_impossible(&recognize(&narayana[1..], &first_zero)),
        0
    );
    assert_eq!(
        assert_impossible(&recognize(&narayana[2..], &rational())),
        0
    );
    assert_eq!(assert_impossible(&recognize(&type_b[1..], &rational())), 0);
    assert_eq!(
        assert_impossible(&recognize(&binomial_eulerian[1..], &rational())),
        0
    );
    assert_eq!(
        assert_impossible(&recognize(&binomial_eulerian, &first_zero)),
        0
    );
}

#[test]
fn certificate_does_not_imply_real_rooted_rows() {
    // The rows of the TNN matrix [[1,0,0],[1,1,0],[1,1,1]]: recognized, yet
    // 1 + t + t^2 has no real zeros.  The certificate is about TNN of the
    // coefficient matrix, not about the rows.
    let rows = rows_i64(&[&[1, 1], &[1, 1, 1]]);
    let digraph = certified_digraph(recognize(&rows, &integer()), &rows);
    assert_eq!(
        digraph.outgoing,
        vec![
            vec![edge(0, br(1)), edge(1, br(1))],
            vec![edge(2, br(1))],
            vec![]
        ]
    );
    assert!(!is_real_rooted_bigint_coeffs(&rows[1]));
}

#[test]
fn planarity_rules_out_a_square() {
    // [1 + t + t^2, (1 + t + t^2)^2]: every column has nonnegative solutions,
    // but target 1 needs an edge from height 1 or 2 while height 0 already
    // reaches height 2.
    let rows = rows_i64(&[&[1, 1, 1], &[1, 2, 3, 2, 1]]);
    let checked = assert_impossible(&recognize(&rows, &rational()));
    assert!(checked > 0);
    assert_impossible(&recognize(&rows, &integer()));
}

#[test]
fn underdetermined_rows_use_the_linear_program() {
    // v_1 = (1,1,1); W_1 + W_2 = (0,0,1,2,1) has a one-parameter family of
    // planar solutions.
    let rows = rows_i64(&[&[1, 1, 1], &[1, 1, 2, 2, 1]]);
    certified_digraph(recognize(&rows, &rational()), &rows);
    certified_digraph(recognize(&rows, &integer()), &rows);
}

#[test]
fn rational_weights_can_be_necessary() {
    // [2t, 2t + t^2] forces the weight 1/2 on the edge 1 -> 2.
    let rows = rows_i64(&[&[0, 2], &[0, 2, 1]]);
    let digraph = certified_digraph(recognize(&rows, &rational()), &rows);
    assert_eq!(
        digraph.outgoing[1],
        vec![
            edge(1, br(1)),
            edge(2, BigRational::from(Ratio::new(bi(1), bi(2))))
        ]
    );
    assert_impossible(&recognize(&rows, &integer()));
}

#[test]
fn target_span_bound_is_respected() {
    let rows = pascal_rows(4);
    let narrow = BrentiRecognitionOptions {
        max_target_span: Some(1),
        ..Default::default()
    };
    assert!(assert_impossible(&recognize(&rows, &narrow)) > 0);
    let two = BrentiRecognitionOptions {
        max_target_span: Some(2),
        ..Default::default()
    };
    certified_digraph(recognize(&rows, &two), &rows);
}

#[test]
fn obstructions_and_input_errors() {
    let negative = rows_i64(&[&[1, -1]]);
    assert_eq!(assert_impossible(&recognize(&negative, &rational())), 0);

    // A zero row followed by a nonzero row.
    let dead = rows_i64(&[&[0], &[1]]);
    assert_impossible(&recognize(&dead, &rational()));

    assert_eq!(
        recognize_brenti_strip_digraph(&[], &rational()),
        Err(BrentiRecognitionError::EmptyInput)
    );
    let bad_offset = BrentiRecognitionOptions {
        first_row_index: 2,
        ..Default::default()
    };
    assert_eq!(
        recognize_brenti_strip_digraph(&polys(&pascal_rows(2)), &bad_offset),
        Err(BrentiRecognitionError::UnsupportedFirstRowIndex { found: 2 })
    );
    let first_zero = BrentiRecognitionOptions {
        first_row_index: 0,
        ..Default::default()
    };
    assert!(matches!(
        recognize_brenti_strip_digraph(&polys(&pascal_rows(2)), &first_zero),
        Err(BrentiRecognitionError::FirstRowNotOne { .. })
    ));
    let zero_span = BrentiRecognitionOptions {
        max_target_span: Some(0),
        ..Default::default()
    };
    assert_eq!(
        recognize_brenti_strip_digraph(&polys(&pascal_rows(2)), &zero_span),
        Err(BrentiRecognitionError::ZeroTargetSpan)
    );

    // P_0 = 1 alone, read as row 0: the empty digraph.
    let only_one = rows_i64(&[&[1]]);
    certified_digraph(recognize(&only_one, &first_zero), &only_one);
}

#[test]
fn tiny_budget_is_inconclusive() {
    let rows = rows_i64(&[&[1, 1, 1], &[1, 2, 3, 2, 1]]);
    let options = BrentiRecognitionOptions {
        max_search_nodes: 1,
        ..Default::default()
    };
    assert!(matches!(
        recognize(&rows, &options),
        BrentiRecognition::Inconclusive { .. }
    ));
}

struct Lcg(u64);

impl Lcg {
    fn next(&mut self, modulus: u64) -> usize {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 33) % modulus) as usize
    }
}

/// A random planar weight matrix: nondecreasing target intervals with small
/// integer weights.
fn random_planar_weights(rng: &mut Lcg, size: usize) -> Vec<Vec<BigRational>> {
    let mut cuts: Vec<usize> = (0..=size).map(|_| rng.next(size as u64)).collect();
    cuts.sort_unstable();
    (0..size)
        .map(|y| {
            (0..size)
                .map(|z| {
                    if cuts[y] <= z && z <= cuts[y + 1] {
                        br(rng.next(3) as i64)
                    } else {
                        br(0)
                    }
                })
                .collect()
        })
        .collect()
}

#[test]
fn random_planar_digraphs_are_recovered() {
    let mut rng = Lcg(0x9e37_79b9_7f4a_7c15);
    let mut recovered = 0;
    for _ in 0..40 {
        let size = 2 + rng.next(5);
        let num_rows = 1 + rng.next(4);
        let source = BrentiStripDigraph::from_weight_matrix(&random_planar_weights(&mut rng, size));
        assert!(source.validate().is_ok());
        let rows: Vec<Vec<BigInt>> = (1..=num_rows)
            .map(|n| {
                let counts = source.row_counts(n).unwrap();
                let mut row: Vec<BigInt> = counts.iter().map(|c| c.to_integer()).collect();
                while row.len() > 1 && row.last() == Some(&bi(0)) {
                    row.pop();
                }
                row
            })
            .collect();
        for options in [rational(), integer()] {
            let digraph = certified_digraph(recognize(&rows, &options), &rows);
            // LGV: the coefficient matrix (rows P_0..P_N, padded) is TNN.
            let width = 1 + rows.iter().map(Vec::len).max().unwrap();
            let mut matrix = vec![vec![bi(0); width]; rows.len() + 1];
            matrix[0][0] = bi(1);
            for (n, row) in rows.iter().enumerate() {
                for (k, c) in row.iter().enumerate() {
                    matrix[n + 1][k] = c.clone();
                }
            }
            assert!(
                check_totally_nonnegative_neville_bigint(&matrix).is_ok(),
                "{:?}",
                matrix
            );
            assert!(digraph.validate().is_ok());
        }
        recovered += 1;
    }
    assert_eq!(recovered, 40);
}
