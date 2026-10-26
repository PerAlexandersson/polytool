//! Benchmark exact real-rootedness checks on representative polynomial families.
//!
//! Run with:
//!
//! ```text
//! cargo run --release -p polytool --example bench_positive_real_rooted
//! ```

use num_bigint::BigInt;
use polytool::sequences::{
    chebyshev_polynomials_t_bigint, chebyshev_polynomials_u_bigint, eulerian_polynomials_bigint,
    hermite_polynomials_bigint, narayana_polynomials_bigint, type_b_eulerian_polynomials_bigint,
};
use polytool::vec_poly::{mul, trim};
use polytool::{
    is_real_rooted_bezout_bigint_coeffs, is_real_rooted_bigint_coeffs,
    is_real_rooted_prs_bigint_coeffs, is_real_rooted_uspensky_bigint_coeffs,
    primitive_sturm_max_coefficient_bits,
};
use std::hint::black_box;
use std::time::{Duration, Instant};

fn product_linear(deg: usize) -> Vec<BigInt> {
    let mut p = vec![BigInt::from(1)];
    for a in 1..=deg {
        p = mul(&p, &[BigInt::from(a), BigInt::from(1)]);
    }
    p
}

fn sequence_entry(
    degree: usize,
    make_sequence: impl FnOnce(usize) -> Vec<Vec<BigInt>>,
) -> Vec<BigInt> {
    make_sequence(degree + 1)
        .into_iter()
        .find(|p| p.len().saturating_sub(1) == degree)
        .expect("sequence does not contain the requested degree")
}

fn eulerian_polynomial(degree: usize) -> Vec<BigInt> {
    sequence_entry(degree, eulerian_polynomials_bigint)
}

fn narayana_polynomial(degree: usize) -> Vec<BigInt> {
    sequence_entry(degree, narayana_polynomials_bigint)
}

fn type_b_eulerian_polynomial(degree: usize) -> Vec<BigInt> {
    sequence_entry(degree, type_b_eulerian_polynomials_bigint)
}

fn chebyshev_t_polynomial(degree: usize) -> Vec<BigInt> {
    sequence_entry(degree, chebyshev_polynomials_t_bigint)
}

fn chebyshev_u_polynomial(degree: usize) -> Vec<BigInt> {
    sequence_entry(degree, chebyshev_polynomials_u_bigint)
}

fn hermite_polynomial(degree: usize) -> Vec<BigInt> {
    sequence_entry(degree, hermite_polynomials_bigint)
}

fn touchard_polynomial(degree: usize) -> Vec<BigInt> {
    let mut p = vec![BigInt::from(1)];
    for _ in 0..degree {
        let mut next = vec![BigInt::from(0); p.len() + 1];
        for (k, coefficient) in p.iter().enumerate() {
            next[k + 1] += coefficient;
            if k > 0 {
                next[k] += coefficient * BigInt::from(k);
            }
        }
        p = trim(&next);
    }
    p
}

fn time_it(name: &str, f: impl FnOnce() -> bool) -> (bool, Duration) {
    let t0 = Instant::now();
    let result = black_box(f());
    let elapsed = t0.elapsed();
    println!(
        "  {:<18} {:<5} {:>10.3} ms",
        name,
        result,
        elapsed.as_secs_f64() * 1000.0
    );
    (result, elapsed)
}

fn time_silent(f: impl FnOnce() -> bool) -> (bool, Duration) {
    let t0 = Instant::now();
    let result = black_box(f());
    (result, t0.elapsed())
}

fn bench_case(name: &str, p: &[BigInt], run_bezout: bool) {
    let degree = p.len().saturating_sub(1);
    let coeff_bits = p.iter().map(|c| c.bits()).max().unwrap_or(0);
    let prs_bits = (degree <= 40).then(|| primitive_sturm_max_coefficient_bits(p));
    println!();
    let prs_bits_display = prs_bits
        .map(|bits| bits.to_string())
        .unwrap_or_else(|| "skipped".to_string());
    println!(
        "{name}: degree={degree}, max_coeff_bits={coeff_bits}, max_prs_bits={prs_bits_display}"
    );
    let (fast, _) = time_it("default-fast", || is_real_rooted_bigint_coeffs(p));
    let (prs, _) = time_it("primitive-prs", || is_real_rooted_prs_bigint_coeffs(p));
    let (uspensky, _) = time_it("uspensky", || is_real_rooted_uspensky_bigint_coeffs(p));
    assert_eq!(fast, prs);
    assert_eq!(fast, uspensky);

    if run_bezout {
        let (bezout, _) = time_it("bezout-bigint", || is_real_rooted_bezout_bigint_coeffs(p));
        assert_eq!(fast, bezout);
    } else {
        println!("  {:<18} skipped", "bezout-bigint");
    }
}

fn cutoff_sweep(name: &str, degrees: &[usize], make_poly: impl Fn(usize) -> Vec<BigInt>) {
    println!();
    println!("{name} cutoff sweep");
    println!(
        "{:>6} {:>8} {:>14} {:>14} {:>14} {:>14} {:>10}",
        "deg", "bits", "default ms", "prs ms", "uspensky ms", "bezout ms", "winner"
    );
    println!("{}", "-".repeat(91));

    for &degree in degrees {
        let p = make_poly(degree);
        let bits = p.iter().map(|c| c.bits()).max().unwrap_or(0);
        let (default_rr, default_t) = time_silent(|| is_real_rooted_bigint_coeffs(&p));
        let (prs_rr, prs_t) = time_silent(|| is_real_rooted_prs_bigint_coeffs(&p));
        let (uspensky_rr, uspensky_t) = time_silent(|| is_real_rooted_uspensky_bigint_coeffs(&p));
        assert_eq!(default_rr, prs_rr);
        assert_eq!(default_rr, uspensky_rr);

        let run_bezout = degree <= 30;
        let (bezout_rr, bezout_t) = if run_bezout {
            let r = time_silent(|| is_real_rooted_bezout_bigint_coeffs(&p));
            assert_eq!(default_rr, r.0);
            (Some(r.0), Some(r.1))
        } else {
            (None, None)
        };
        black_box(bezout_rr);

        let default_ms = default_t.as_secs_f64() * 1000.0;
        let prs_ms = prs_t.as_secs_f64() * 1000.0;
        let uspensky_ms = uspensky_t.as_secs_f64() * 1000.0;
        let bezout_ms = bezout_t.map(|t| t.as_secs_f64() * 1000.0);
        let mut winner = ("default", default_ms);
        for candidate in [("prs", prs_ms), ("uspensky", uspensky_ms)] {
            if candidate.1 < winner.1 {
                winner = candidate;
            }
        }
        if let Some(bezout_ms) = bezout_ms {
            if bezout_ms < winner.1 {
                winner = ("bezout", bezout_ms);
            }
        }

        let bezout_display = bezout_ms
            .map(|ms| format!("{ms:14.3}"))
            .unwrap_or_else(|| format!("{:>14}", "skipped"));
        println!(
            "{degree:>6} {bits:>8} {default_ms:>14.3} {prs_ms:>14.3} \
             {uspensky_ms:>14.3} {bezout_display} {:>10}",
            winner.0,
        );
    }
}

fn family_sweep(name: &str, degrees: &[usize], make_poly: impl Fn(usize) -> Vec<BigInt>) {
    println!();
    println!("{name} PRS/Uspensky sweep");
    println!(
        "{:>6} {:>8} {:>8} {:>14} {:>14} {:>14} {:>10}",
        "deg", "bits", "pal?", "default ms", "prs ms", "uspensky ms", "winner"
    );
    println!("{}", "-".repeat(85));

    for &degree in degrees {
        let p = make_poly(degree);
        let bits = p.iter().map(|c| c.bits()).max().unwrap_or(0);
        let first_nonzero = p
            .iter()
            .position(|c| c != &BigInt::from(0))
            .unwrap_or(p.len());
        let core = &p[first_nonzero..];
        let palindromic = core.iter().eq(core.iter().rev());
        let (default_rr, default_t) = time_silent(|| is_real_rooted_bigint_coeffs(&p));
        let (prs_rr, prs_t) = time_silent(|| is_real_rooted_prs_bigint_coeffs(&p));
        let (uspensky_rr, uspensky_t) = time_silent(|| is_real_rooted_uspensky_bigint_coeffs(&p));
        assert_eq!(default_rr, prs_rr);
        assert_eq!(prs_rr, uspensky_rr);

        let default_ms = default_t.as_secs_f64() * 1000.0;
        let prs_ms = prs_t.as_secs_f64() * 1000.0;
        let uspensky_ms = uspensky_t.as_secs_f64() * 1000.0;
        let winner = if uspensky_ms < prs_ms {
            "uspensky"
        } else {
            "prs"
        };
        println!(
            "{degree:>6} {bits:>8} {palindromic:>8} {default_ms:>14.3} {prs_ms:>14.3} \
             {uspensky_ms:>14.3} {winner:>10}"
        );
    }
}

fn main() {
    println!("Exact real-rootedness benchmark on representative families");
    bench_case("prod_{a=1}^{30} (x+a)", &product_linear(30), true);
    bench_case("prod_{a=1}^{80} (x+a)", &product_linear(80), false);
    bench_case("Eulerian degree 35", &eulerian_polynomial(35), true);
    bench_case("Eulerian degree 79", &eulerian_polynomial(79), false);
    bench_case(
        "small non-real-rooted",
        &[
            BigInt::from(1),
            BigInt::from(43),
            BigInt::from(196),
            BigInt::from(168),
            BigInt::from(23),
            BigInt::from(1),
        ],
        true,
    );

    cutoff_sweep(
        "prod_{a=1}^{d} (x+a)",
        &[5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80],
        product_linear,
    );
    cutoff_sweep(
        "Eulerian degree d",
        &[5, 10, 15, 20, 25, 30, 35, 40],
        eulerian_polynomial,
    );

    let family_degrees = [10, 20, 25, 30, 35, 40];
    family_sweep("Narayana", &family_degrees, narayana_polynomial);
    family_sweep(
        "type-B Eulerian",
        &family_degrees,
        type_b_eulerian_polynomial,
    );
    family_sweep("Touchard", &family_degrees, touchard_polynomial);
    family_sweep("Chebyshev T", &family_degrees, chebyshev_t_polynomial);
    family_sweep("Chebyshev U", &family_degrees, chebyshev_u_polynomial);
    family_sweep("Hermite", &family_degrees, hermite_polynomial);
}
