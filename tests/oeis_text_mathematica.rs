//! Evaluates the generated Mathematica programs (`%t`) with `wolframscript`
//! and compares their flattened rows with polytool's rows.
//!
//! Opt-in: the test runs only with `POLYTOOL_TEST_WOLFRAMSCRIPT=1` and a
//! `wolframscript` on the `PATH`; otherwise it returns early.

use num_bigint::BigInt;
use polytool::oeis::{
    BFileRowOrientation, BFileZeroPolicy, OeisRecurrencePresentation, OeisTriangleReading,
};
use polytool::recurrence::RecurrenceJson;
use std::process::Command;

const THREE_TERM: &str = r#"{"schema": "polytool.recurrence.v1", "first_index": 0,
  "initial_polynomials": [["1"], ["1", "1"]],
  "recurrence": {"terms": [
    {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["2", "1"]]}},
    {"offset": 2, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["1", "1"]]}}]}}"#;

const DERIVATIVE: &str = r#"{"schema": "polytool.recurrence.v1", "first_index": 0,
  "initial_polynomials": [["1"]],
  "recurrence": {"terms": [
    {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["1", "-2"], ["0", "2"]]}},
    {"offset": 1, "deriv_order": 1, "sign": "none", "coeff": {"coeffs": [["0", "2", "-2"]]}}]}}"#;

const INDEX_DEPENDENT: &str = r#"{"schema": "polytool.recurrence.v1", "first_index": 0,
  "initial_polynomials": [["1"]],
  "recurrence": {"terms": [
    {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["0", "1"], ["1"]]}}]}}"#;

const DENOMINATOR: &str = r#"{"schema": "polytool.recurrence.v1", "first_index": 1,
  "initial_polynomials": [["1"], ["1", "1"]],
  "recurrence": {"terms": [
    {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["-1", "-1"], ["2", "2"]]}},
    {"offset": 2, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["2", "-4", "2"], ["-1", "2", "-1"]]}}],
   "denominator": {"coeffs": [["1"], ["1"]]}}}"#;

/// `P(n) = x (1 + x) P(n-1) + (-1)^n x^2 P(n-2)`: low zeros to drop.
const LOW_ZEROS: &str = r#"{"schema": "polytool.recurrence.v1", "first_index": 0,
  "initial_polynomials": [["1"], ["0", "1"]],
  "recurrence": {"terms": [
    {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["0", "1", "1"]]}},
    {"offset": 2, "deriv_order": 0, "sign": "alternating_n", "coeff": {"coeffs": [["0", "0", "1"]]}}]}}"#;

fn wolframscript_enabled() -> bool {
    std::env::var("POLYTOOL_TEST_WOLFRAMSCRIPT").as_deref() == Ok("1")
        && Command::new("wolframscript")
            .arg("-version")
            .output()
            .is_ok_and(|output| output.status.success())
}

fn evaluate(program: &str) -> Vec<BigInt> {
    let output = Command::new("wolframscript")
        .arg("-code")
        .arg(program)
        .output()
        .expect("run wolframscript");
    assert!(output.status.success(), "wolframscript failed: {output:?}");
    let text = String::from_utf8(output.stdout).expect("utf-8 output");
    let inner = text
        .trim()
        .strip_prefix('{')
        .and_then(|rest| rest.strip_suffix('}'))
        .unwrap_or_else(|| panic!("unexpected output {text:?} for\n{program}"));
    inner
        .split(',')
        .map(|term| {
            term.trim()
                .parse()
                .unwrap_or_else(|_| panic!("non-integer term {term:?}"))
        })
        .collect()
}

#[test]
fn mathematica_programs_reproduce_the_rows() {
    if !wolframscript_enabled() {
        eprintln!("skipped: set POLYTOOL_TEST_WOLFRAMSCRIPT=1 with wolframscript installed");
        return;
    }
    let reversed_regular = OeisTriangleReading {
        orientation: BFileRowOrientation::Reversed,
        regular_triangle: true,
        ..OeisTriangleReading::default()
    };
    let offset_one = OeisTriangleReading {
        first_row_index: 1,
        ..OeisTriangleReading::default()
    };
    let trimmed_reversed = OeisTriangleReading {
        first_row_index: 2,
        orientation: BFileRowOrientation::Reversed,
        zero_policy: BFileZeroPolicy::TrimLowOrder,
        regular_triangle: false,
    };
    let cases = [
        (THREE_TERM, None, OeisTriangleReading::default()),
        (DERIVATIVE, Some(1), reversed_regular),
        (INDEX_DEPENDENT, None, offset_one),
        (DENOMINATOR, None, offset_one),
        (LOW_ZEROS, Some(1), trimmed_reversed),
    ];
    for (json, first_row, reading) in cases {
        let record: RecurrenceJson = serde_json::from_str(json).unwrap();
        let presentation =
            OeisRecurrencePresentation::from_recurrence_json(&record, first_row, reading).unwrap();
        let text = presentation.oeis_text(9).unwrap();
        let expected: Vec<BigInt> = presentation
            .triangle_rows(9)
            .unwrap()
            .into_iter()
            .flatten()
            .collect();
        assert_eq!(
            evaluate(&text.mathematica),
            expected,
            "program:\n{}",
            text.mathematica
        );
    }
}
