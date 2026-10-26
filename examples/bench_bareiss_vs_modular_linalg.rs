//! Compare BigInt Bareiss and modular/CRT leading principal minors.
//!
//! Run with:
//!
//! ```text
//! cargo run --release -p polytool --example bench_bareiss_vs_modular_linalg
//! ```

use num_bigint::BigInt;
use num_traits::{One, Zero};
use polytool::vec_poly::{derivative, trim};
use polytool::{
    bareiss_leading_principal_minors_bigint, bezout_matrix_bigint_coeffs,
    modular_leading_principal_minors_bigint,
};
use std::time::Instant;

fn product_linear(degree: usize) -> Vec<BigInt> {
    let mut p = vec![BigInt::one()];
    for a in 1..=degree {
        let mut next = vec![BigInt::zero(); p.len() + 1];
        for (i, c) in p.iter().enumerate() {
            next[i] += BigInt::from(a) * c;
            next[i + 1] += c;
        }
        p = trim(&next);
    }
    p
}

fn bench_degree(degree: usize) {
    let f = product_linear(degree);
    let fp = derivative(&f);
    let bezout = bezout_matrix_bigint_coeffs(&f, &fp).expect("degree difference one");

    let start = Instant::now();
    let bareiss = bareiss_leading_principal_minors_bigint(&bezout)
        .expect("Bareiss should succeed for this positive definite matrix");
    let bareiss_elapsed = start.elapsed();

    let start = Instant::now();
    let modular = modular_leading_principal_minors_bigint(&bezout)
        .expect("modular CRT reconstruction should certify this matrix");
    let modular_elapsed = start.elapsed();

    assert_eq!(bareiss, modular);
    println!(
        "degree {:>3}: matrix {:>3}x{:<3}  bareiss {:>9.3?}  modular {:>9.3?}",
        degree,
        bezout.len(),
        bezout.len(),
        bareiss_elapsed,
        modular_elapsed
    );
}

fn main() {
    println!("Bezout leading-principal-minor benchmark");
    println!("f(t)=prod_(a=1)^d (t+a), comparing B(f,f')");
    for degree in [12, 20, 30, 40] {
        bench_degree(degree);
    }
}
