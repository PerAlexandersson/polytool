//! Benchmark exact square solves with integer and fractional inputs.
//!
//! Run with `--release -- --size 48 --repeat 5`. JSON Lines contain the exact
//! solution so separate executable builds can be compared without timing fields.

use num_bigint::BigInt;
use num_traits::Zero;
use polytool::linalg::solve_full_rank_square_linear_system;
use polytool::BigRational as Q;
use std::hint::black_box;
use std::time::Instant;

fn main() {
    let mut size = 40usize;
    let mut repeat = 5usize;
    let mut only = String::new();
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--size" => size = args.next().unwrap().parse().unwrap(),
            "--repeat" => repeat = args.next().unwrap().parse().unwrap(),
            "--only" => only = args.next().unwrap(),
            _ => panic!("use --size N --repeat N --only SUBSTRING"),
        }
    }
    assert!((1..=256).contains(&size) && repeat > 0);
    let denominators = [2i64, 3, 5, 7, 11, 13];
    for case in [
        "integer",
        "row_scaled",
        "mixed",
        "rhs_only",
        "hilbert",
        "singular_early",
        "singular_late",
        "bidiagonal",
        "block_diagonal",
    ] {
        if !case.contains(&only) {
            continue;
        }
        // Strict diagonal dominance ensures nonsingularity in the first four
        // cases. Hilbert matrices are positive definite; the sparse cases are
        // nonsingular triangular or strictly diagonally dominant block matrices.
        // The singular cases deliberately destroy full rank. These are matrix
        // workloads, not polynomial sequences.
        let mut matrix: Vec<Vec<Q>> = (0..size)
            .map(|i| {
                (0..size)
                    .map(|j| {
                        if case == "hilbert" {
                            return Q::new(BigInt::from(1), BigInt::from(i + j + 1));
                        }
                        if case == "bidiagonal" {
                            return Q::from_integer(BigInt::from(if i == j {
                                2
                            } else if j == i + 1 {
                                -1
                            } else {
                                0
                            }));
                        }
                        if case == "block_diagonal" && i / 8 != j / 8 {
                            return Q::zero();
                        }
                        let value = if i == j {
                            100 * size as i64 + 1
                        } else {
                            ((13 * i + 7 * j + 3 * i * j) % 11) as i64 - 5
                        };
                        let denominator = match case {
                            "row_scaled" => denominators[i % denominators.len()],
                            "mixed" => denominators[(i + 2 * j) % denominators.len()],
                            _ => 1,
                        };
                        Q::new(BigInt::from(value), BigInt::from(denominator))
                    })
                    .collect()
            })
            .collect();
        if case == "singular_early" {
            // The first coefficient column vanishes, so rank is at most n-1.
            for row in &mut matrix {
                row[0] = Q::zero();
            }
        } else if case == "singular_late" {
            // A dense system whose dependence is exposed at the last pivot.
            matrix[size - 1] = if size == 1 {
                vec![Q::zero()]
            } else {
                matrix[0].clone()
            };
        }
        let witness: Vec<Q> = (0..size)
            .map(|j| {
                Q::new(
                    BigInt::from((j % 5) as i64 - 2),
                    BigInt::from(if case == "rhs_only" { 3 } else { 1 }),
                )
            })
            .collect();
        let rhs: Vec<Q> = matrix
            .iter()
            .map(|row| row.iter().zip(&witness).map(|(a, x)| a * x).sum())
            .collect();
        let expected = if case.starts_with("singular_") {
            None
        } else {
            Some(witness)
        };
        let input_bits = matrix
            .iter()
            .flatten()
            .chain(&rhs)
            .map(|x| x.numer().bits().max(x.denom().bits()))
            .max()
            .unwrap();
        // One warm-up is excluded. Matrix generation and exact comparison are
        // outside the solve timer, which includes conversion and elimination.
        for run in 0..=repeat {
            let start = Instant::now();
            let actual = solve_full_rank_square_linear_system(black_box(&matrix), black_box(&rhs));
            let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
            assert_eq!(actual, expected);
            if run > 0 {
                println!(
                    "{}",
                    serde_json::json!({
                        "case": case, "size": size, "run": run,
                        "elapsed_ms": elapsed_ms, "input_bits": input_bits,
                        "solution": actual.as_ref().map(|solution| solution.iter().map(ToString::to_string).collect::<Vec<_>>()),
                    })
                );
            }
        }
    }
}
