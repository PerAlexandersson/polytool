//! Generate exact statistic-polynomial rows for OEIS A355089.
//!
//! Usage: `cargo run --release -p combpoly --example a355089_statistics -- 20 des`

use combpoly::parity_pattern::{parity_preserving_123_statistic_polynomials, Parity123Statistic};
use std::time::Instant;

fn main() {
    let mut args = std::env::args().skip(1);
    let max_n: u8 = args
        .next()
        .as_deref()
        .unwrap_or("20")
        .parse()
        .expect("maximum n must be an integer in 0..=63");
    let statistic_name = args.next().unwrap_or_else(|| "des".to_owned());
    let statistic = match statistic_name.as_str() {
        "des" => Parity123Statistic::Descents,
        "asc" => Parity123Statistic::Ascents,
        "maj" => Parity123Statistic::MajorIndex,
        "amaj" => Parity123Statistic::AscentMajorIndex,
        "lrmin" => Parity123Statistic::LeftToRightMinima,
        _ => panic!("statistic must be des, asc, maj, amaj, or lrmin"),
    };

    let started = Instant::now();
    let rows = parity_preserving_123_statistic_polynomials(max_n, statistic);
    for row in rows {
        let coefficients: Vec<_> = row.iter().map(ToString::to_string).collect();
        println!("[{}]", coefficients.join(", "));
    }
    eprintln!(
        "computed statistic={statistic_name} n=0..={max_n} in {:.3?}",
        started.elapsed()
    );
}
