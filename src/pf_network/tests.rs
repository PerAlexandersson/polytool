//! Tests for the enclosing module.

use super::recognition::{SearchOutcome, SupportSearch};
use super::*;
use crate::numeric::binomial_bigint;

fn q(v: i64) -> Q {
    Q::from_integer(BigInt::from(v))
}

fn step(length: usize, rise: usize, weight: i64) -> PfStep {
    PfStep {
        length,
        rise,
        weight: q(weight),
    }
}

fn polys(rows: &[&[i64]]) -> Vec<Polynomial<BigInt>> {
    rows.iter()
        .map(|r| Polynomial::from_i64_coeffs(r))
        .collect()
}

fn pascal_rows(last: usize) -> Vec<Polynomial<BigInt>> {
    (1..=last)
        .map(|n| Polynomial::new((0..=n).map(|k| binomial_bigint(n, k)).collect()))
        .collect()
}

fn eulerian_rows(last: usize) -> Vec<Polynomial<BigInt>> {
    // A(n, k): permutations of [n + 1] with k descents, rows n = 1..=last.
    let mut rows: Vec<Vec<BigInt>> = vec![vec![BigInt::one()]];
    for m in 1..=last {
        let prev = &rows[m - 1];
        let row = (0..=m)
            .map(|k| {
                let left = if k >= 1 {
                    BigInt::from(m - k + 1) * &prev[k - 1]
                } else {
                    BigInt::zero()
                };
                let right = prev
                    .get(k)
                    .map_or_else(BigInt::zero, |c| BigInt::from(k + 1) * c);
                left + right
            })
            .collect();
        rows.push(row);
    }
    rows.into_iter().skip(1).map(Polynomial::new).collect()
}

#[test]
fn segment_conflicts() {
    // Proper crossing.
    assert!(segments_conflict((0, 0), (1, 2), (0, 1), (1, 1)));
    // Common endpoint only.
    assert!(!segments_conflict((0, 0), (1, 1), (0, 0), (1, 0)));
    assert!(!segments_conflict((0, 0), (2, 1), (1, 0), (2, 1)));
    // An endpoint in the interior of the other segment.
    assert!(segments_conflict((0, 0), (2, 2), (1, 1), (2, 1)));
    // Collinear overlap versus collinear touching.
    assert!(segments_conflict((0, 0), (2, 0), (1, 0), (3, 0)));
    assert!(!segments_conflict((0, 0), (1, 0), (1, 0), (2, 0)));
    // Parallel, disjoint.
    assert!(!segments_conflict((0, 0), (1, 1), (0, 1), (1, 2)));
}

#[test]
fn compatibility_of_step_types() {
    let t = |length, rise| PfStepType { length, rise };
    // Brenti's recursive-matrix pattern is planar.
    for (a, b) in [(t(1, 0), t(1, 1)), (t(1, 0), t(2, 1)), (t(1, 1), t(2, 1))] {
        assert!(step_types_compatible(a, b), "{:?} {:?}", a, b);
        assert!(step_types_compatible(b, a), "{:?} {:?}", b, a);
    }
    // Rises 0 and 2 in one column cross.
    assert!(!step_types_compatible(t(1, 0), t(1, 2)));
    // Non-primitive steps are never allowed.
    assert!(!step_types_compatible(t(2, 2), t(1, 0)));
    assert!(!step_types_compatible(t(2, 0), t(2, 0)));
    // Every primitive type is compatible with itself.
    for s in primitive_step_types(4, 4) {
        assert!(step_types_compatible(s, s), "{:?}", s);
    }
}

#[test]
fn pascal_network_rows_and_toeplitz() {
    let network = HeightInvariantNetwork {
        steps: (0..4).map(|_| vec![step(1, 0, 1), step(1, 1, 1)]).collect(),
    };
    network.validate().unwrap();
    assert_eq!(
        network.row_coefficients(4).unwrap(),
        vec![q(1), q(4), q(6), q(4), q(1)]
    );
    let matrix = network.toeplitz_path_matrix(3, 5).unwrap();
    let row = network.row_coefficients(3).unwrap();
    assert_eq!(matrix, toeplitz_window(&row, 5));
    let cert = build_pf_rows_certificate(network, &[1, 2, 3, 4], &pascal_rows(4)).unwrap();
    cert.verify_against_bigint_polynomials(&pascal_rows(4))
        .unwrap();
}

#[test]
fn row_polynomials_follow_the_three_term_recurrence() {
    // Network for P_n = (1 + t) P_{n-1} + 2 t P_{n-2}, P_0 = 1: steps of length 1
    // with rises 0 and 1 (weight 1) and, from column 2 on, a step of length 2,
    // rise 1 and weight 2.  The oracle is the recurrence on coefficient vectors.
    let last = 7;
    let network = HeightInvariantNetwork {
        steps: (1..=last)
            .map(|n| {
                let mut column = vec![step(1, 0, 1), step(1, 1, 1)];
                if n >= 2 {
                    column.push(step(2, 1, 2));
                }
                column
            })
            .collect(),
    };
    network.validate().unwrap();
    let mut rows: Vec<Vec<i64>> = vec![vec![1]];
    for n in 1..=last {
        let mut next = vec![0i64; n + 1];
        for (k, c) in rows[n - 1].iter().enumerate() {
            next[k] += c;
            next[k + 1] += c;
        }
        if n >= 2 {
            for (k, c) in rows[n - 2].iter().enumerate() {
                next[k + 1] += 2 * c;
            }
        }
        while next.last() == Some(&0) {
            next.pop();
        }
        rows.push(next);
    }
    for (n, expected) in rows.iter().enumerate() {
        let expected_poly = Polynomial::new(expected.iter().map(|&v| q(v)).collect());
        assert_eq!(network.row_polynomial(n).unwrap(), expected_poly, "row {n}");
    }
    let prefix = network.row_polynomials_prefix(last).unwrap();
    assert_eq!(prefix.len(), last + 1);
    for (n, poly) in prefix.iter().enumerate() {
        assert_eq!(poly, &network.row_polynomial(n).unwrap());
    }
    assert!(matches!(
        network.row_polynomial(last + 1),
        Err(PfNetworkError::RowOutOfRange { .. })
    ));
}

#[test]
fn invalid_networks_rejected() {
    let crossing = HeightInvariantNetwork {
        steps: vec![vec![step(1, 0, 1), step(1, 2, 1)]],
    };
    assert!(matches!(
        crossing.validate(),
        Err(PfNetworkError::Crossing { .. })
    ));
    // The crossing network would produce the row 1 + t^2, which is not PF.
    let through = HeightInvariantNetwork {
        steps: vec![vec![], vec![step(2, 2, 1)]],
    };
    assert!(matches!(
        through.validate(),
        Err(PfNetworkError::StepThroughVertex { .. })
    ));
    let early = HeightInvariantNetwork {
        steps: vec![vec![step(2, 1, 1)]],
    };
    assert!(matches!(
        early.validate(),
        Err(PfNetworkError::StepBeforeSource { .. })
    ));
    let negative = HeightInvariantNetwork {
        steps: vec![vec![step(1, 0, -1)]],
    };
    assert!(matches!(
        negative.validate(),
        Err(PfNetworkError::NegativeWeight { .. })
    ));
    let duplicate = HeightInvariantNetwork {
        steps: vec![vec![step(1, 0, 1), step(1, 0, 2)]],
    };
    assert!(matches!(
        duplicate.validate(),
        Err(PfNetworkError::DuplicateStep { .. })
    ));
    // Zero-weight edges are ignored by the planarity check.
    let zero = HeightInvariantNetwork {
        steps: vec![vec![step(1, 0, 1), step(1, 2, 0)]],
    };
    zero.validate().unwrap();
    // Crossing between columns: a long step over a short one.
    let long = HeightInvariantNetwork {
        steps: vec![vec![step(1, 0, 1)], vec![step(1, 2, 1), step(2, 1, 1)]],
    };
    assert!(matches!(
        long.validate(),
        Err(PfNetworkError::Crossing { .. })
    ));
}

#[test]
fn finite_toeplitz_window_alone_proves_nothing() {
    let to_int = |m: Vec<Vec<Q>>| -> Vec<Vec<BigInt>> {
        m.into_iter()
            .map(|r| r.into_iter().map(|v| v.to_integer()).collect())
            .collect()
    };
    let row = vec![q(1), q(1), q(1)];
    assert!(
        crate::check_totally_nonnegative_neville_bigint(&to_int(toeplitz_window(&row, 3))).is_ok()
    );
    assert!(
        crate::check_totally_nonnegative_neville_bigint(&to_int(toeplitz_window(&row, 4))).is_err()
    );
}

#[test]
fn recognizes_pascal_rows() {
    for weights in [BrentiWeightDomain::Rational, BrentiWeightDomain::Integer] {
        let options = PfNetworkOptions {
            weights,
            ..PfNetworkOptions::default()
        };
        let rows = pascal_rows(6);
        match recognize_pf_network(&rows, &options).unwrap() {
            PfNetworkRecognition::Certified { certificate, .. } => {
                certificate
                    .verify_against_bigint_polynomials(&rows)
                    .unwrap();
                assert_eq!(certificate.row_indices, (1..=6).collect::<Vec<_>>());
            }
            other => panic!("expected a certificate, got {}", other),
        }
    }
    // With P_0 given explicitly.
    let mut rows = vec![Polynomial::from_i64_coeffs(&[1])];
    rows.extend(pascal_rows(3));
    let options = PfNetworkOptions {
        first_row_index: 0,
        ..PfNetworkOptions::default()
    };
    assert!(recognize_pf_network(&rows, &options)
        .unwrap()
        .is_certified());
}

#[test]
fn recognizes_three_term_recurrence_rows() {
    // P_n = (1 + t) P_{n-1} + 2 t P_{n-2}: Brenti's recursive pattern.
    let mut rows: Vec<Vec<i64>> = vec![vec![1], vec![1, 1]];
    for n in 2..7 {
        let mut next = vec![0i64; n + 1];
        for (k, c) in rows[n - 1].iter().enumerate() {
            next[k] += c;
            next[k + 1] += c;
        }
        for (k, c) in rows[n - 2].iter().enumerate() {
            next[k + 1] += 2 * c;
        }
        rows.push(next);
    }
    let input: Vec<Polynomial<BigInt>> = rows[1..]
        .iter()
        .map(|r| Polynomial::from_i64_coeffs(r))
        .collect();
    let result = recognize_pf_network(&input, &PfNetworkOptions::default()).unwrap();
    assert!(result.is_certified(), "{}", result);
    // Only single steps: the rows are not products of linear factors.
    let single = PfNetworkOptions {
        max_step_length: 1,
        ..PfNetworkOptions::default()
    };
    let result = recognize_pf_network(&input, &single).unwrap();
    assert!(
        matches!(result, PfNetworkRecognition::ImpossibleWithinBounds { .. }),
        "{}",
        result
    );
}

#[test]
fn never_certifies_one_plus_t_plus_t_squared() {
    for rows in [
        polys(&[&[1, 1, 1]]),
        polys(&[&[1, 1], &[1, 1, 1]]),
        polys(&[&[1, 2], &[1, 1, 1]]),
    ] {
        let result = recognize_pf_network(&rows, &PfNetworkOptions::default()).unwrap();
        assert!(
            matches!(result, PfNetworkRecognition::ImpossibleWithinBounds { .. }),
            "{}",
            result
        );
    }
    // The support search alone (without the real-rootedness shortcut)
    // also finds nothing, even with long steps and large rises.
    let inputs: [&[&[i64]]; 3] = [
        &[&[1, 1, 1]],
        &[&[1, 1], &[1, 1, 1]],
        &[&[1, 2], &[1, 1, 1]],
    ];
    for rows in inputs {
        let mut vectors: Vec<Vec<BigInt>> = vec![vec![BigInt::one()]];
        vectors.extend(
            rows.iter()
                .map(|r| r.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>()),
        );
        let options = PfNetworkOptions {
            max_step_length: 4,
            max_rise: 4,
            ..PfNetworkOptions::default()
        };
        let mut search = SupportSearch::new(&vectors, &options);
        assert!(matches!(search.run(), SearchOutcome::Impossible(_)));
    }
}

#[test]
fn eulerian_and_narayana_rows() {
    // A_1..A_3 and the Narayana rows N_1..N_3 fit Brenti's recursive
    // pattern (1,0), (1,1), (2,1); with the fourth rows no network has
    // steps of length <= 4 and rise <= 4 (exact, within these bounds).
    for (rows, name) in [
        (eulerian_rows(4), "eulerian"),
        (narayana_rows(4), "narayana"),
    ] {
        let result = recognize_pf_network(&rows[..3], &PfNetworkOptions::default()).unwrap();
        match &result {
            PfNetworkRecognition::Certified { certificate } => {
                certificate
                    .verify_against_bigint_polynomials(&rows[..3])
                    .unwrap();
            }
            other => panic!("{}: {}", name, other),
        }
        let options = PfNetworkOptions {
            max_step_length: 4,
            max_rise: 4,
            ..PfNetworkOptions::default()
        };
        let result = recognize_pf_network(&rows, &options).unwrap();
        assert!(
            matches!(result, PfNetworkRecognition::ImpossibleWithinBounds { .. }),
            "{}: {}",
            name,
            result
        );
    }
}

fn narayana_rows(last: usize) -> Vec<Polynomial<BigInt>> {
    (1..=last)
        .map(|n| {
            Polynomial::new(
                (0..=n)
                    .map(|k| {
                        binomial_bigint(n + 1, k + 1) * binomial_bigint(n + 1, k)
                            / BigInt::from(n + 1)
                    })
                    .collect(),
            )
        })
        .collect()
}

#[test]
fn negative_coefficient_and_bad_options() {
    let result = recognize_pf_network(&polys(&[&[1, -1]]), &PfNetworkOptions::default()).unwrap();
    assert!(matches!(
        result,
        PfNetworkRecognition::ImpossibleWithinBounds { .. }
    ));
    assert_eq!(
        recognize_pf_network(&[], &PfNetworkOptions::default()),
        Err(PfNetworkRecognitionError::EmptyInput)
    );
    let options = PfNetworkOptions {
        first_row_index: 0,
        ..PfNetworkOptions::default()
    };
    assert!(matches!(
        recognize_pf_network(&polys(&[&[1, 1]]), &options),
        Err(PfNetworkRecognitionError::FirstRowNotOne { .. })
    ));
    let options = PfNetworkOptions {
        max_step_length: 0,
        ..PfNetworkOptions::default()
    };
    assert_eq!(
        recognize_pf_network(&polys(&[&[1, 1]]), &options),
        Err(PfNetworkRecognitionError::ZeroStepLength)
    );
}

#[test]
fn search_budget_is_inconclusive() {
    let options = PfNetworkOptions {
        max_search_nodes: 2,
        ..PfNetworkOptions::default()
    };
    let result = recognize_pf_network(&pascal_rows(5), &options).unwrap();
    assert!(
        matches!(result, PfNetworkRecognition::Inconclusive { .. }),
        "{}",
        result
    );
}

#[test]
fn search_finds_every_random_valid_network() {
    // Completeness check against the independent validator: build random
    // networks, keep the planar ones, and recognize their rows.
    let mut state: u64 = 0x9e37_79b9_7f4a_7c15;
    let mut next = |modulus: u64| -> u64 {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state >> 33) % modulus
    };
    let types = primitive_step_types(3, 3);
    let (mut valid, mut tried) = (0, 0);
    while valid < 60 && tried < 20_000 {
        tried += 1;
        let last = 1 + next(4) as usize;
        let mut steps: Vec<Vec<PfStep>> = Vec::with_capacity(last);
        for n in 1..=last {
            let mut column = Vec::new();
            for t in types.iter().filter(|t| t.length <= n) {
                if next(3) == 0 {
                    column.push(step(t.length, t.rise, 1 + next(3) as i64));
                }
            }
            steps.push(column);
        }
        let network = HeightInvariantNetwork { steps };
        if network.validate().is_err() {
            continue;
        }
        valid += 1;
        let rows: Vec<Polynomial<BigInt>> = network
            .row_polynomials_prefix(last)
            .unwrap()
            .into_iter()
            .skip(1)
            .map(|p| Polynomial::new(p.coeffs().iter().map(|c| c.to_integer()).collect()))
            .collect();
        for c in network.row_polynomials_prefix(last).unwrap() {
            let coeffs: Vec<BigInt> = c.coeffs().iter().map(|v| v.to_integer()).collect();
            assert!(coeffs.is_empty() || is_real_rooted_bigint_coeffs(&coeffs));
        }
        for weights in [BrentiWeightDomain::Rational, BrentiWeightDomain::Integer] {
            let options = PfNetworkOptions {
                weights,
                ..PfNetworkOptions::default()
            };
            let result = recognize_pf_network(&rows, &options).unwrap();
            assert!(result.is_certified(), "{:?}: {}", network, result);
        }
    }
    assert!(
        valid >= 60,
        "only {} valid networks in {} tries",
        valid,
        tried
    );
}
