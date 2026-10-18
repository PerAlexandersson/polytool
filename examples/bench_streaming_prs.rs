//! Reproduce the deterministic corpus used to benchmark exact root counting.
//!
//! Run through degree 100 with:
//!
//! ```text
//! cargo run --release -p polytool --example bench_streaming_prs -- 100
//! ```
//!
//! To compare implementations, run this command at each commit under the same
//! build profile and machine load. The timed region contains only the public
//! real-rootedness call; deterministic corpus construction and square-free
//! validation happen beforehand.

use num_bigint::BigInt;
use num_traits::{One, Zero};
use polytool::{is_real_rooted_bigint_coeffs, squarefree_degree_bigint_coeffs};
use std::hint::black_box;
use std::time::{Duration, Instant};

fn multiply(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let mut result = vec![BigInt::zero(); a.len() + b.len() - 1];
    for (i, left) in a.iter().enumerate() {
        for (j, right) in b.iter().enumerate() {
            result[i + j] += left * right;
        }
    }
    result
}

fn main() {
    let maximum_degree = std::env::args()
        .nth(1)
        .and_then(|value| value.parse().ok())
        .unwrap_or(20usize);
    let mut state = 0x243f_6a88_85a3_08d3u64;

    println!("degree\tsamples\tns_per_poly\trepetitions\tchecksum");
    for degree in 2usize..=maximum_degree {
        let mut corpus = Vec::new();
        for sample in 0..64usize {
            let polynomial = if sample < 48 {
                let mut coefficients = Vec::with_capacity(degree + 1);
                for _ in 0..=degree {
                    state = state
                        .wrapping_mul(6_364_136_223_846_793_005)
                        .wrapping_add(1_442_695_040_888_963_407);
                    coefficients.push(BigInt::from(((state >> 32) % 21) as i64 - 10));
                }
                if coefficients[degree].is_zero() {
                    coefficients[degree] = BigInt::one();
                }
                coefficients
            } else {
                let mut polynomial = vec![BigInt::one()];
                for root_index in 0..degree {
                    let root = BigInt::from(1 + ((sample - 48) + 3 * root_index) % 23);
                    polynomial = multiply(&polynomial, &[root, BigInt::one()]);
                }
                polynomial
            };
            let mut polynomial = polynomial;
            while squarefree_degree_bigint_coeffs(&polynomial) != degree {
                polynomial[0] += BigInt::one();
            }
            corpus.push(polynomial);
        }

        let mut checksum = 0usize;
        let mut repetitions = 0u32;
        let start = Instant::now();
        while repetitions == 0
            || (repetitions < 100_000 && start.elapsed() < Duration::from_millis(50))
        {
            for polynomial in &corpus {
                checksum += usize::from(black_box(is_real_rooted_bigint_coeffs(black_box(
                    polynomial,
                ))));
            }
            repetitions += 1;
        }
        let nanoseconds =
            start.elapsed().as_nanos() as f64 / (f64::from(repetitions) * corpus.len() as f64);
        println!(
            "{degree}\t{}\t{nanoseconds:.1}\t{repetitions}\t{checksum}",
            corpus.len()
        );
    }
}
