//! Discover the sharp scale in a third-order recurrence and replay the proof.

use polytool::recurrence::{
    discover_interlacing_refinements, verify_interlacing_refinements, BigRational, BivarPoly,
    InterlacingRefinementRequest, Recurrence, RecurrenceJson, RecurrenceTerm,
};

fn main() {
    let q = |n: i64| BigRational::from_integer(n.into());
    // P_(n+2)=(1+x)P_(n+1)+4x P_n+4x^2 P_(n-1).
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(
                1,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![q(1), q(1)]]),
            ),
            RecurrenceTerm::new(
                2,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(4)]]),
            ),
            RecurrenceTerm::new(
                3,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(0), q(4)]]),
            ),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let record = RecurrenceJson::from_recurrence_rational(
        &rule,
        0,
        &[vec![q(1)], vec![q(1)], vec![q(1), q(2)]],
        None,
    );
    let request: InterlacingRefinementRequest = serde_json::from_value(serde_json::json!({
        "recurrence": record,
        "limits": {"maximum_candidates": 2}
    }))
    .unwrap();
    let report = discover_interlacing_refinements(&request).unwrap();
    assert_eq!(report.status, "proved");
    assert_eq!(
        report
            .scaled_interval
            .as_ref()
            .unwrap()
            .rational_scale
            .as_deref(),
        Some("2")
    );
    verify_interlacing_refinements(&report).unwrap();
    println!("{}", serde_json::to_string_pretty(&report).unwrap());
}
