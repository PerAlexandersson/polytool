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
        record.to_latex_expanded().unwrap(),
        concat!(
            "$P_{5} = 1$\n$P_{6} = 1 + t$\n$P_{7} = 0$\n",
            "$\\bigl(2n + 4\\bigr)\\, P_{n} = \\bigl(2n + 2\\bigr)\\, P_{n-1}",
            " - t^{2}\\, P_{n-2}^{(3)} - 2\\, (-1)^n P_{n-1}''",
            " - \\bigl(2n - 2t\\bigr) \\quad (n \\geq 8)$"
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

fn erik_record(extra: serde_json::Value) -> RecurrenceJson {
    let mut record = json!({
        "schema": "polytool.recurrence.v1", "first_index": 0,
        "initial_polynomials": [["1"], ["1", "1"]],
        "recurrence": {
            "denominator": {"coeffs": [["2"], ["1"]]},
            "terms": [
                {"offset": 1, "deriv_order": 0, "sign": "none",
                 "coeff": {"coeffs": [["1", "1"], ["2", "2"]]}},
                {"offset": 2, "deriv_order": 0, "sign": "none",
                 "coeff": {"coeffs": [["1", "-2", "1"], ["-1", "2", "-1"]]}}
            ]
        }
    });
    if let Some(object) = extra.as_object() {
        for (key, value) in object {
            record[key] = value.clone();
        }
    }
    serde_json::from_value(record).unwrap()
}

#[test]
fn coefficients_are_factored_over_the_rationals() {
    let record = erik_record(json!({}));
    assert_eq!(
        record.to_latex().unwrap(),
        concat!(
            "$P_{0} = 1$\n$P_{1} = 1 + t$\n",
            "$(n + 2) P_{n} = (2n + 1)(1 + t) P_{n-1} - (n - 1)(1 - t)^{2} P_{n-2}",
            " \\quad (n \\geq 2)$"
        )
    );
    // The expanded presentation is still available and is a different text.
    assert!(record
        .to_latex_expanded()
        .unwrap()
        .contains("\\bigl(n + 2\\bigr)"));
}

#[test]
fn factored_text_re_expands_to_the_stored_coefficients() {
    use polytool::recurrence::BivarPoly;
    let record = erik_record(json!({}));
    let (rule, _, _) = record.to_recurrence_parts().unwrap();
    let mut coefficients: Vec<BivarPoly> =
        rule.terms.iter().map(|term| term.coeff.clone()).collect();
    coefficients.extend(rule.denominator.clone());
    for coefficient in coefficients {
        let factorization = coefficient.irreducible_factorization().unwrap();
        assert_eq!(factorization.expand(), coefficient);
        assert!(factorization
            .factors
            .iter()
            .all(|(_, multiplicity)| *multiplicity >= 1));
    }
}

#[test]
fn search_records_keep_their_candidate_status_line() {
    let record = erik_record(json!({
        "search": {
            "recurrence_text": "(n + 2) P(n) = ...", "source_rows": 2, "skip_prefix": 0,
            "unknowns": 8, "weighted_unknowns": 8, "equations": 8, "fit_polynomials": 2,
            "verification_polynomials": 0, "candidates_tried": 1,
            "options": {
                "var_deg": 2, "idx_deg": 1, "diff_deg": 0, "rec_len": 2,
                "homogeneous": true, "inhomo_var_deg": 1, "inhomo_idx_deg": 1,
                "denom_var_deg": 0, "denom_idx_deg": 1, "alternating_sign": false
            }
        }
    }));
    let latex = record.to_latex().unwrap();
    assert!(
        latex.starts_with("% Candidate fitted to finite rows; validity for all n is unproved.\n")
    );
    assert!(latex.contains("(2n + 1)(1 + t) P_{n-1}"));
}
