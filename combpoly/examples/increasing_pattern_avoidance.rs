//! Print exact descent and peak polynomials for `Av_n(12...m)`.

use combpoly::monotone_pattern::{
    increasing_pattern_avoiding_descent_polynomials, increasing_pattern_avoiding_peak_polynomials,
};

fn main() {
    let mut arguments = std::env::args().skip(1);
    let max_n = arguments
        .next()
        .map(|value| value.parse::<usize>().expect("max_n must be an integer"))
        .unwrap_or(10);
    let pattern_length = arguments
        .next()
        .map(|value| {
            value
                .parse::<usize>()
                .expect("pattern_length must be an integer")
        })
        .unwrap_or(4);
    assert!(arguments.next().is_none(), "usage: max_n [pattern_length]");

    let descents = increasing_pattern_avoiding_descent_polynomials(max_n, pattern_length);
    let peaks = increasing_pattern_avoiding_peak_polynomials(max_n, pattern_length);
    println!("avoiding the increasing pattern of length {pattern_length}");
    for n in 0..=max_n {
        println!("n={n}: des={:?}; peak={:?}", descents[n], peaks[n]);
    }
}
