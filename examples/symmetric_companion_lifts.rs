//! Recognize the symmetric parts of the `1/2`-Eulerian polynomials, given
//! only their second-order scalar recurrence and the seeds `A_2 = 1 + x`,
//! `A_3 = 1 + 7x + x^2`, as the symmetric part of the parent
//! `P_(n+1) = (1 + 2n x) P_n + 2x (1 - x) P_n'`, `P_2 = 1 + 2x`, and prove
//! real-rootedness of every `A_n` through the preservation theorem.
//!
//! Run with `cargo run --example symmetric_companion_lifts`.

use num_bigint::BigInt;
use polytool::recurrence::{
    find_symmetric_companion_lifts, verify_symmetric_companion_certificate, BigRational, BivarPoly,
    Recurrence, RecurrenceError, RecurrenceTerm, SymmetricCompanionLiftOptions,
};

fn q(value: i64) -> BigRational {
    BigRational::from_integer(BigInt::from(value))
}

fn term(offset: usize, derivative: usize, rows: &[&[i64]]) -> RecurrenceTerm {
    let coeffs = rows
        .iter()
        .map(|row| row.iter().map(|&c| q(c)).collect())
        .collect();
    RecurrenceTerm::new(
        offset,
        derivative,
        BivarPoly::from_coefficients_normalized(coeffs),
    )
}

fn main() -> Result<(), RecurrenceError> {
    // k = 2, target index m, coefficient rows [m^i][x^j]:
    // A_m = [1 + (4m - 7) x] A_(m-1) + 4x(1 - x) A'_(m-1)
    //       - [(6m - 15) x + (4m^2 - 22m + 30) x^2] A_(m-2)
    //       - x(1 - x)[6 + (8m - 26) x] A'_(m-2) - 4x^2 (1 - x)^2 A''_(m-2).
    let recurrence = Recurrence {
        terms: vec![
            term(1, 0, &[&[1, -7], &[0, 4]]),
            term(1, 1, &[&[0, 4, -4]]),
            term(2, 0, &[&[0, 15, -30], &[0, -6, 22], &[0, 0, -4]]),
            term(2, 1, &[&[0, -6, 32, -26], &[0, 0, -8, 8]]),
            term(2, 2, &[&[0, 0, -4, 8, -4]]),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let seeds = vec![vec![q(1), q(1)], vec![q(1), q(7), q(1)]];
    let report = find_symmetric_companion_lifts(
        &recurrence,
        &seeds,
        2,
        &SymmetricCompanionLiftOptions::default(),
    )?;
    println!("schedule: {}", report.symmetry_schedule.rule);
    for attempt in &report.search.attempts {
        println!(
            "attempt (degree {}): {:?}, {}",
            attempt.parameter_degree_in_n, attempt.outcome, attempt.detail
        );
    }
    for lift in &report.lifts {
        println!("companion: {}", lift.companion_formula);
        println!("parent:    {}", lift.parent_rule);
        println!(
            "status:    identity {:?}, real-rootedness {:?}, Lean {}",
            lift.verification.identity_status,
            lift.verification.real_rootedness_status,
            lift.verification.lean_replay
        );
        for check in lift
            .verification
            .identities
            .iter()
            .chain(&lift.verification.hypotheses)
        {
            println!("  {:<32} {:?}", check.name, check.status);
        }
        println!("  {}", lift.verification.conclusion);
        let replay = verify_symmetric_companion_certificate(&lift.certificate, 10)?;
        println!("independent replay certified: {}", replay.is_certified());
    }
    Ok(())
}
