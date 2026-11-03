use super::*;

fn record(json: &str) -> RecurrenceJson {
    serde_json::from_str(json).expect("valid recurrence record")
}

/// `P(n) = (2 + x) P(n-1) + (1 + x) P(n-2)`, `P(0) = 1`, `P(1) = 1 + x`.
pub(super) const THREE_TERM: &str = r#"{
  "schema": "polytool.recurrence.v1", "first_index": 0,
  "initial_polynomials": [["1"], ["1", "1"]],
  "recurrence": {"terms": [
    {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["2", "1"]]}},
    {"offset": 2, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["1", "1"]]}}]}
}"#;

/// `P(n) = (1 - 2x + 2nx) P(n-1) + (2x - 2x^2) P'(n-1)`, `P(0) = 1`: the rows
/// of A156919 reversed, from the row `n = 1` on.
pub(super) const DERIVATIVE: &str = r#"{
  "schema": "polytool.recurrence.v1", "first_index": 0,
  "initial_polynomials": [["1"]],
  "recurrence": {"terms": [
    {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["1", "-2"], ["0", "2"]]}},
    {"offset": 1, "deriv_order": 1, "sign": "none", "coeff": {"coeffs": [["0", "2", "-2"]]}}]}
}"#;

/// `P(n) = (n + x) P(n-1)`, `P(0) = 1` (unsigned Stirling numbers of the
/// first kind).
pub(super) const INDEX_DEPENDENT: &str = r#"{
  "schema": "polytool.recurrence.v1", "first_index": 0,
  "initial_polynomials": [["1"]],
  "recurrence": {"terms": [
    {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["0", "1"], ["1"]]}}]}
}"#;

/// Narayana polynomials: `(n + 1) N(n) = (2n - 1)(1 + x) N(n-1) - (n - 2)(1 - x)^2 N(n-2)`,
/// `N(1) = 1`, `N(2) = 1 + x`.
pub(super) const DENOMINATOR: &str = r#"{
  "schema": "polytool.recurrence.v1", "first_index": 1,
  "initial_polynomials": [["1"], ["1", "1"]],
  "recurrence": {"terms": [
    {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["-1", "-1"], ["2", "2"]]}},
    {"offset": 2, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["2", "-4", "2"], ["-1", "2", "-1"]]}}],
   "denominator": {"coeffs": [["1"], ["1"]]}}
}"#;

fn int_rows(rows: &[&[i64]]) -> Vec<Vec<BigInt>> {
    rows.iter()
        .map(|row| row.iter().map(|&c| BigInt::from(c)).collect())
        .collect()
}

#[test]
fn three_term_constant_coefficients() {
    let presentation = OeisRecurrencePresentation::from_recurrence_json(
        &record(THREE_TERM),
        None,
        OeisTriangleReading::default(),
    )
    .unwrap();
    let text = presentation.oeis_text(5).unwrap();
    assert_eq!(
        text.recurrence_sentence,
        "The row polynomials satisfy P(n,x) = (2+x)*P(n-1,x) + (1+x)*P(n-2,x) for n >= 2, \
         with P(0,x) = 1, P(1,x) = 1 + x."
    );
    assert_eq!(
        text.triangle_example,
        "Triangle begins:\n  0:  1;\n  1:  1,  1;\n  2:  3,  4,  1;\n  3:  7, 13,  7,  1;\n  \
         4: 17, 40, 32, 10, 1;\n"
    );
    assert_eq!(
        text.mathematica,
        "P[0, x_] := 1;\nP[1, x_] := 1 + x;\n\
         P[n_, x_] := P[n, x] = Expand[(2 + x)*P[n - 1, x] + (1 + x)*P[n - 2, x]];\n\
         Table[CoefficientList[P[n, x], x], {n, 0, 4}] // Flatten"
    );
    assert_eq!(
        text.pari,
        "P(N) = {my(v = vector(N + 1)); v[1] = 1; v[2] = 1 + x; \
         for(n = 2, N, my(i = n + 1); v[i] = (2 + x)*v[i - 1] + (1 + x)*v[i - 2]); v};\n\
         row(n, p) = {my(c = Vecrev(p)); c};\n\
         T(N) = {my(v = P(N)); concat(vector(#v, i, row(i - 1, v[i])))};\n\
         T(4)"
    );
}

#[test]
fn derivative_reversed_with_index_shift() {
    // A156919: offset 0, row n of the entry is the reversed node row n + 1.
    let reading = OeisTriangleReading {
        first_row_index: 0,
        orientation: BFileRowOrientation::Reversed,
        zero_policy: BFileZeroPolicy::KeepAll,
        regular_triangle: true,
    };
    let presentation =
        OeisRecurrencePresentation::from_recurrence_json(&record(DERIVATIVE), Some(1), reading)
            .unwrap();
    let text = presentation.oeis_text(4).unwrap();
    assert_eq!(
        text.recurrence_sentence,
        "Row n lists the coefficients of P(n,x) in descending powers of x, where \
         P(n,x) = (1+2*n*x)*P(n-1,x) + (2*x-2*x^2)*d/dx P(n-1,x) for n >= 1, with P(0,x) = 1."
    );
    assert_eq!(
        text.triangle_example,
        "Triangle begins:\n  0: 1;\n  1: 2,  1;\n  2: 4, 10,  1;\n  3: 8, 60, 36, 1;\n"
    );
    assert_eq!(
        text.mathematica,
        "P[0, x_] := 1;\n\
         P[n_, x_] := P[n, x] = Expand[(1 + 2*n*x)*P[n - 1, x] + (2*x - 2*x^2)*D[P[n - 1, x], x]];\n\
         Table[Reverse[PadRight[CoefficientList[P[n, x], x], n + 1]], {n, 0, 3}] // Flatten"
    );
    assert_eq!(
        text.pari,
        "P(N) = {my(v = vector(N + 1)); v[1] = 1; \
         for(n = 1, N, my(i = n + 1); v[i] = (1 + 2*n*x)*v[i - 1] + (2*x - 2*x^2)*deriv(v[i - 1], x)); v};\n\
         row(n, p) = {my(c = Vecrev(p)); c = concat(c, vector(n + 1 - #c)); c = Vecrev(c); c};\n\
         T(N) = {my(v = P(N)); concat(vector(#v, i, row(i - 1, v[i])))};\n\
         T(3)"
    );
}

#[test]
fn index_dependent_coefficient_with_offset_one() {
    let reading = OeisTriangleReading {
        first_row_index: 1,
        ..OeisTriangleReading::default()
    };
    let presentation =
        OeisRecurrencePresentation::from_recurrence_json(&record(INDEX_DEPENDENT), None, reading)
            .unwrap();
    let text = presentation.oeis_text(4).unwrap();
    assert_eq!(
        text.recurrence_sentence,
        "The row polynomials satisfy P(n,x) = (n-1+x)*P(n-1,x) for n >= 2, with P(1,x) = 1."
    );
    assert_eq!(
        text.triangle_example,
        "Triangle begins:\n  1: 1;\n  2: 1,  1;\n  3: 2,  3, 1;\n  4: 6, 11, 6, 1;\n"
    );
    assert_eq!(
        text.mathematica,
        "P[1, x_] := 1;\n\
         P[n_, x_] := P[n, x] = Expand[(n - 1 + x)*P[n - 1, x]];\n\
         Table[CoefficientList[P[n, x], x], {n, 1, 4}] // Flatten"
    );
    assert_eq!(
        text.pari,
        "P(N) = {my(v = vector(N)); v[1] = 1; for(n = 2, N, my(i = n); v[i] = (n - 1 + x)*v[i - 1]); v};\n\
         row(n, p) = {my(c = Vecrev(p)); c};\n\
         T(N) = {my(v = P(N)); concat(vector(#v, i, row(i, v[i])))};\n\
         T(4)"
    );
}

#[test]
fn denominator_is_divided_out_exactly() {
    let reading = OeisTriangleReading {
        first_row_index: 1,
        ..OeisTriangleReading::default()
    };
    let presentation =
        OeisRecurrencePresentation::from_recurrence_json(&record(DENOMINATOR), None, reading)
            .unwrap();
    let text = presentation.oeis_text(5).unwrap();
    assert_eq!(
        text.recurrence_sentence,
        "The row polynomials satisfy (n+1)*P(n,x) = (2*n-1+(2*n-1)*x)*P(n-1,x) \
         - (n-2-(2*n-4)*x+(n-2)*x^2)*P(n-2,x) for n >= 3, with P(1,x) = 1, P(2,x) = 1 + x."
    );
    assert_eq!(
        text.triangle_example,
        "Triangle begins:\n  1: 1;\n  2: 1,  1;\n  3: 1,  3,  1;\n  4: 1,  6,  6,  1;\n  \
         5: 1, 10, 20, 10, 1;\n"
    );
    assert_eq!(
        text.mathematica,
        "P[1, x_] := 1;\nP[2, x_] := 1 + x;\n\
         P[n_, x_] := P[n, x] = Expand[Together[((2*n - 1 + (2*n - 1)*x)*P[n - 1, x] \
         - (n - 2 - (2*n - 4)*x + (n - 2)*x^2)*P[n - 2, x])/(n + 1)]];\n\
         Table[CoefficientList[P[n, x], x], {n, 1, 5}] // Flatten"
    );
    assert!(text
        .pari
        .contains("v[i] = ((2*n - 1 + (2*n - 1)*x)*v[i - 1]"));
    assert!(text.pari.contains(")/(n + 1)); v};"));
}

#[test]
fn low_zero_trimming_and_rows() {
    // P(n) = x P(n-1), P(0) = 1: rows x^n, read without the low zeros.
    let json = r#"{
      "schema": "polytool.recurrence.v1", "first_index": 0,
      "initial_polynomials": [["1"]],
      "recurrence": {"terms": [
        {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["0", "1"]]}}]}
    }"#;
    let reading = OeisTriangleReading {
        zero_policy: BFileZeroPolicy::TrimLowOrder,
        ..OeisTriangleReading::default()
    };
    let presentation =
        OeisRecurrencePresentation::from_recurrence_json(&record(json), None, reading).unwrap();
    assert_eq!(
        presentation.triangle_rows(3).unwrap(),
        int_rows(&[&[1], &[1], &[1]])
    );
    let text = presentation.oeis_text(3).unwrap();
    assert!(text
        .recurrence_sentence
        .starts_with("Row n lists the coefficients of P(n,x) in ascending powers of x, starting at the lowest nonzero coefficient, where P(n,x) = x*P(n-1,x)"));
    assert!(text
        .mathematica
        .contains("CoefficientList[Expand[P[n, x]/x^Exponent[P[n, x], x, Min]], x]"));
    assert!(text.pari.contains("Vecrev(p/x^valuation(p, x))"));
}

#[test]
fn regular_triangle_pads_and_stops() {
    let rows = int_rows(&[&[1], &[1, 1], &[2], &[1, 2, 3, 4, 5]]);
    let reading = OeisTriangleReading {
        regular_triangle: true,
        ..OeisTriangleReading::default()
    };
    assert_eq!(
        reading.read_rows(&rows),
        int_rows(&[&[1], &[1, 1], &[2, 0, 0]])
    );
}

#[test]
fn alternating_sign_and_inhomogeneous_terms() {
    // P(n) = (-1)^n P(n-1) + n x, P(0) = 1.
    let json = r#"{
      "schema": "polytool.recurrence.v1", "first_index": 0,
      "initial_polynomials": [["1"]],
      "recurrence": {"terms": [
        {"offset": 1, "deriv_order": 0, "sign": "alternating_n", "coeff": {"coeffs": [["1"]]}}],
       "inhomogeneous": {"coeffs": [["0", "0"], ["0", "1"]]}}
    }"#;
    let presentation = OeisRecurrencePresentation::from_recurrence_json(
        &record(json),
        None,
        OeisTriangleReading::default(),
    )
    .unwrap();
    assert_eq!(
        presentation.format_recurrence_sentence(),
        "The row polynomials satisfy P(n,x) = (-1)^n*P(n-1,x) + n*x for n >= 1, with P(0,x) = 1."
    );
    assert_eq!(
        presentation.triangle_rows(3).unwrap(),
        int_rows(&[&[1], &[-1, 1], &[-1, 3]])
    );
}

#[test]
fn first_row_before_record_is_rejected() {
    let error = OeisRecurrencePresentation::from_recurrence_json(
        &record(DENOMINATOR),
        Some(0),
        OeisTriangleReading::default(),
    )
    .unwrap_err();
    assert_eq!(
        error,
        OeisTextError::FirstRowBeforeRecurrence {
            first_row: 0,
            first_index: 1
        }
    );
}

#[test]
fn rational_rows_are_reported() {
    // P(n) = P(n-1)/2.
    let json = r#"{
      "schema": "polytool.recurrence.v1", "first_index": 0,
      "initial_polynomials": [["1"]],
      "recurrence": {"terms": [
        {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["1/2"]]}}]}
    }"#;
    let presentation = OeisRecurrencePresentation::from_recurrence_json(
        &record(json),
        None,
        OeisTriangleReading::default(),
    )
    .unwrap();
    assert_eq!(
        presentation.triangle_rows(3).unwrap_err(),
        OeisTextError::NonIntegralCoefficient { row_index: 1 }
    );
    assert_eq!(
        presentation.format_recurrence_sentence(),
        "The row polynomials satisfy P(n,x) = 1/2*P(n-1,x) for n >= 1, with P(0,x) = 1."
    );
}
