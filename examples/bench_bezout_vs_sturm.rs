//! Benchmark: Bézout matrix vs Sturm chain for real-rootedness and interlacing.
//!
//! Run with: cargo run --release --example bench_bezout_vs_sturm

use polytool::real_rootedness::*;
use std::time::Instant;

/// Largest `n` for which the Sturm interlacing check is timed.
const STURM_INTERLACING_MAX_N: usize = 10;

fn main() {
    // A_21 has a coefficient above i64::MAX, so the i64 benchmark stops at A_20.
    let max_n = 20;
    // polys[i] = A_{i+1}(t).
    let polys = polytool::sequences::eulerian_polynomials(max_n)[1..].to_vec();

    println!("=== Real-rootedness: Sturm vs explicit Bézout vs default ===");
    println!(
        "{:>4} {:>8} {:>12} {:>12} {:>12}",
        "n", "deg", "Sturm (μs)", "Bézout (μs)", "Default (μs)"
    );
    println!("{}", "-".repeat(57));

    for (i, p) in polys.iter().enumerate() {
        let n = i + 1;
        if n < 3 {
            continue;
        }

        let iters = if n < 10 {
            100
        } else if n < 20 {
            10
        } else {
            3
        };

        let t0 = Instant::now();
        for _ in 0..iters {
            let _ = is_real_rooted_sturm(p);
        }
        let sturm_us = t0.elapsed().as_micros() as f64 / iters as f64;

        let t0 = Instant::now();
        for _ in 0..iters {
            let _ = is_real_rooted_bezout(p);
        }
        let bezout_us = t0.elapsed().as_micros() as f64 / iters as f64;

        let t0 = Instant::now();
        for _ in 0..iters {
            let _ = is_real_rooted(p);
        }
        let default_us = t0.elapsed().as_micros() as f64 / iters as f64;

        let deg = p.len() - 1;
        println!(
            "{:>4} {:>8} {:>12.1} {:>12.1} {:>12.1}  {}",
            n,
            deg,
            sturm_us,
            bezout_us,
            default_us,
            if default_us < bezout_us && default_us < sturm_us {
                "← default"
            } else if bezout_us < sturm_us {
                "← explicit Bézout"
            } else {
                "← Sturm"
            }
        );
    }

    println!();
    println!("=== Interlacing: Sturm vs Bézout (consecutive Eulerian polynomials) ===");
    println!(
        "{:>4} {:>8} {:>12} {:>12}",
        "n", "deg", "Sturm (μs)", "Bézout (μs)"
    );
    println!("{}", "-".repeat(42));

    for i in 1..polys.len() {
        let n = i + 1;
        if n < 4 {
            continue;
        }
        // check_interlacing(f, g) tests f ≪ g with deg g = deg f + 1.
        let f = &polys[i - 1]; // A_{n-1}(t), degree n-2
        let g = &polys[i]; // A_{n}(t), degree n-1

        let iters = if n < 10 {
            100
        } else if n < 20 {
            10
        } else {
            3
        };

        // The Sturm interlacing check roughly doubles in cost per degree
        // (about 4 s per call at n = 12), so it is timed only up to n = 10.
        let sturm_us = if n <= STURM_INTERLACING_MAX_N {
            let t0 = Instant::now();
            for _ in 0..iters {
                let _ = check_interlacing_sturm(f, g);
            }
            t0.elapsed().as_micros() as f64 / iters as f64
        } else {
            f64::INFINITY
        };

        let t0 = Instant::now();
        for _ in 0..iters {
            let _ = check_interlacing(f, g);
        }
        let bezout_us = t0.elapsed().as_micros() as f64 / iters as f64;

        let deg = g.len() - 1;
        let result = check_interlacing(f, g);
        println!(
            "{:>4} {:>8} {:>12.1} {:>12.1}  {} (interlaces: {:?})",
            n,
            deg,
            sturm_us,
            bezout_us,
            if bezout_us < sturm_us {
                "← Bézout"
            } else {
                "← Sturm"
            },
            result
        );
    }
}
