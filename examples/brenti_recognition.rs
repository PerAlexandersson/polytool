//! Recognize Pascal rows as rows of a planar strip digraph and check the
//! chain polynomials of a lower-unitriangular TNN matrix.
//!
//! Run with `cargo run --example brenti_recognition`.

use num_bigint::BigInt;
use polytool::numeric::binomial_bigint;
use polytool::{
    check_chain_polynomials_bigint, format_poly_bigint_coeffs, recognize_brenti_strip_digraph,
    BrentiRecognition, BrentiRecognitionOptions, Polynomial,
};

fn main() {
    // Rows (1+t)^n for n = 1..=5, with P_0 = 1 implicit.
    let rows: Vec<Polynomial<BigInt>> = (1..=5)
        .map(|n| Polynomial::new((0..=n).map(|k| binomial_bigint(n, k)).collect()))
        .collect();
    let result = recognize_brenti_strip_digraph(&rows, &BrentiRecognitionOptions::default())
        .expect("valid input");
    println!("{}", result);
    if let BrentiRecognition::Certified { digraph, .. } = &result {
        for (height, edges) in digraph.outgoing.iter().enumerate() {
            let targets: Vec<String> = edges
                .iter()
                .map(|edge| format!("{} (weight {})", edge.to, edge.weight))
                .collect();
            println!("  height {} -> {}", height, targets.join(", "));
        }
    }

    // [[1,0,0],[1,1,0],[1,1,1]] is TNN: its chain polynomials are real-rooted
    // in [-1, 0] and interlace, although its last row 1 + t + t^2 is not
    // real-rooted.
    let matrix: Vec<Vec<BigInt>> = vec![vec![1], vec![1, 1], vec![1, 1, 1]]
        .into_iter()
        .map(|row: Vec<i64>| row.into_iter().map(BigInt::from).collect())
        .collect();
    let report = check_chain_polynomials_bigint(&matrix).expect("lower-unitriangular input");
    println!("\nmatrix TNN: {}", report.is_totally_nonnegative());
    for (n, p) in report.chain_polynomials.iter().enumerate() {
        println!(
            "  p_{} = {}  (real-rooted in [-1,0]: {})",
            n,
            format_poly_bigint_coeffs(p),
            report.roots_in_minus_one_zero[n]
        );
    }
    println!(
        "  consecutive interlacing: {:?}; theorem conclusion holds: {}",
        report.interlacing,
        report.chain_conclusion_holds()
    );
}
