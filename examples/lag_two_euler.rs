//! Certify simple negative roots and strict interlacing for a positive
//! lag-two Euler triangle, using its canonical catalog recurrence.

use polytool::recurrence::{check_lag_two_euler_theorem, RowInterlacing};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // A157637: a(n)=b(n)=1, c(n)=n/2, P_0=1, P_1=1+x.
    let entry = polytool::oeis::by_id("A157637").expect("bundled family");
    let (rule, seeds) = entry.recurrence_parts()?;
    let report = check_lag_two_euler_theorem(
        &rule,
        entry.recurrence_first_index,
        &seeds,
        &Default::default(),
    )?;
    assert!(report.is_proved());
    report.verify()?;
    let conclusion = report.conclusion.as_ref().expect("proved conclusion");
    assert!(conclusion.simple_roots);
    assert_eq!(conclusion.interlacing, Some(RowInterlacing::Strict));
    println!("{}", conclusion.statement);
    Ok(())
}
