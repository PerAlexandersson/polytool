//! LaTeX definitions preserve exact coefficients, indices and recurrence data.

use polytool::recurrence::RecurrenceJson;
use serde_json::json;

#[test]
fn complete_definition_groups_signs_fractions_and_derivatives() {
    let record: RecurrenceJson = serde_json::from_value(json!({
        "schema": "polytool.recurrence.v1", "first_index": 5,
        "initial_polynomials": [["1"], ["1", "1"], ["0"]],
        "recurrence": {
            "denominator": {"coeffs": [["2"], ["1"]]},
            "terms": [
                {"offset": 1, "deriv_order": 0, "sign": "none",
                 "coeff": {"coeffs": [["1"], ["1"]]}},
                {"offset": 2, "deriv_order": 3, "sign": "none",
                 "coeff": {"coeffs": [["0", "0", "-1/2"]]}},
                {"offset": 1, "deriv_order": 2, "sign": "alternating_n",
                 "coeff": {"coeffs": [["-1"]]}},
                {"offset": 1, "deriv_order": 0, "sign": "none",
                 "coeff": {"coeffs": [["0"]]}}
            ],
            "inhomogeneous": {"coeffs": [["0", "1"], ["-1"]]}
        }
    }))
    .unwrap();
    assert_eq!(
        record.to_latex().unwrap(),
        concat!(
            "$P_{5} = 1$\n$P_{6} = 1 + t$\n$P_{7} = 0$\n",
            "$\\bigl(n + 2\\bigr)\\, P_{n} = \\bigl(n + 1\\bigr)\\, P_{n-1}",
            " - \\frac{1}{2}t^{2}\\, P_{n-2}^{(3)} - (-1)^n P_{n-1}''",
            " - \\bigl(n - t\\bigr) \\quad (n \\geq 8)$"
        )
    );
    let original = record.to_recurrence_parts().unwrap().0.to_string();
    let decorated = record.clone().with_latex().unwrap();
    let decoded: RecurrenceJson =
        serde_json::from_str(&serde_json::to_string(&decorated).unwrap()).unwrap();
    assert_eq!(
        decoded.to_recurrence_parts().unwrap().0.to_string(),
        original
    );
    assert_eq!(decoded.latex, Some(record.to_latex().unwrap()));
    let shifted = decoded.reindexed(0).unwrap();
    assert_eq!(shifted.latex, Some(shifted.to_latex().unwrap()));
    assert!(shifted.latex.unwrap().ends_with("(n \\geq 3)$"));
}

#[test]
fn zero_rhs_negative_seed_and_invalid_domains() {
    let mut record: RecurrenceJson = serde_json::from_value(json!({
        "schema": "polytool.recurrence.v1", "first_index": 0,
        "initial_polynomials": [["1", "-1/2"]],
        "recurrence": {"terms": [{"offset": 1, "deriv_order": 1,
            "sign": "none", "coeff": {"coeffs": [["0"]]}}],
            "inhomogeneous": {"coeffs": [["0"]]}}
    }))
    .unwrap();
    assert_eq!(
        record.to_latex().unwrap(),
        "$P_{0} = 1 - \\frac{1}{2}t$\n$P_{n} = 0 \\quad (n \\geq 1)$"
    );
    record.first_index = usize::MAX;
    assert!(record.to_latex().is_err());
    record.first_index = 0;
    record.recurrence.terms[0].offset = 2;
    assert!(record.to_latex().is_err());
}
