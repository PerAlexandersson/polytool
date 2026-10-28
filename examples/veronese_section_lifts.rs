//! Recognize `Q_(n+2) = 2 Q_(n+1) + (y - 1) Q_n`, `Q_0 = Q_1 = 1`, as the even
//! coefficient sections of `(1 + x)^n`, which proves that every `Q_n` is
//! real-rooted with nonpositive roots.
//!
//! Run with `cargo run --example veronese_section_lifts`.

use num_bigint::BigInt;
use polytool::recurrence::{
    find_veronese_section_lifts, format_characteristic_polynomial, BigRational, BivarPoly,
    Recurrence, RecurrenceError, RecurrenceTerm, VeroneseSectionLiftOptions,
};

fn q(value: i64) -> BigRational {
    BigRational::from_integer(BigInt::from(value))
}

fn main() -> Result<(), RecurrenceError> {
    let recurrence = Recurrence {
        terms: vec![
            RecurrenceTerm::new(
                1,
                0,
                BivarPoly {
                    coeffs: vec![vec![q(2)]],
                },
            ),
            RecurrenceTerm::new(
                2,
                0,
                BivarPoly {
                    coeffs: vec![vec![q(-1), q(1)]],
                },
            ),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let initial_rows = vec![vec![q(1)], vec![q(1)]];
    let report = find_veronese_section_lifts(
        &recurrence,
        &initial_rows,
        0,
        &VeroneseSectionLiftOptions::default(),
    )?;
    println!(
        "chi(L, y) = {}",
        format_characteristic_polynomial(&report.input_characteristic_polynomial, "y")
    );
    for lift in &report.lifts {
        println!(
            "modulus {}, residue {}: G(L, x) = {}, parent {} with seeds {:?} -> {}",
            lift.modulus,
            lift.residue,
            format_characteristic_polynomial(&lift.parent_characteristic_polynomial, "x"),
            lift.parent_recurrence,
            lift.parent_recurrence_json().initial_polynomials,
            lift.status.label()
        );
    }
    for skipped in &report.skipped_moduli {
        println!("skipped modulus {}: {}", skipped.modulus, skipped.reason);
    }
    Ok(())
}
