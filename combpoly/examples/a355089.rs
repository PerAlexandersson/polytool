//! Generate exact terms of OEIS A355089 with the count-only parity-word DP.
//!
//! Usage: `cargo run --release -p combpoly --example a355089 -- 40`

use combpoly::parity_pattern::parity_preserving_123_avoiding_counts;
use std::time::Instant;

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .as_deref()
        .unwrap_or("40")
        .parse()
        .expect("maximum n must be an integer in 0..=63");

    let started = Instant::now();
    let counts = parity_preserving_123_avoiding_counts(max_n);
    for (n, count) in counts.iter().enumerate() {
        println!("{n}\t{count}");
    }
    eprintln!("computed n=0..={max_n} in {:.3?}", started.elapsed());
}
