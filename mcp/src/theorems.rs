//! Conditional theorem matches (`polytool::renewal`) for MCP recurrence tools.

use crate::parsing::{invalid_params, validate_mcp_text};
use polytool::recurrence::{
    AdaptiveSearchOptions, AdaptiveSearchResult, BigRational, RecurrenceJson,
};
use polytool::renewal::{match_recurrence_theorems, RecurrenceTheoremReport};
use rmcp::ErrorData as McpError;
use std::collections::BTreeMap;

/// Additive `theorem_matches` field of a found recurrence: rows after
/// `skip_prefix` are labeled from 1, as in its `recurrence_json`.  `None` when
/// no theorem matches, so the field is omitted.
pub(crate) fn search_theorem_matches(
    result: &AdaptiveSearchResult,
    polynomials: &[Vec<BigRational>],
    search: &AdaptiveSearchOptions,
) -> Option<Vec<serde_json::Value>> {
    let searched = polynomials.get(search.skip_prefix..).unwrap_or(&[]);
    let matches = match_recurrence_theorems(&result.recurrence, 1, searched).theorem_matches;
    (!matches.is_empty()).then(|| {
        matches
            .iter()
            .map(|m| serde_json::to_value(m).expect("serialize theorem match"))
            .collect()
    })
}

/// `polytool.recurrence-theorems.v1` report for a supplied recurrence record.
pub(crate) fn supplied_theorem_report(
    recurrence_json: &str,
) -> Result<BTreeMap<String, serde_json::Value>, McpError> {
    validate_mcp_text(recurrence_json, "recurrence JSON").map_err(invalid_params)?;
    let record: RecurrenceJson = serde_json::from_str(recurrence_json)
        .map_err(|error| invalid_params(format!("failed to parse recurrence JSON: {error}")))?;
    let report = RecurrenceTheoremReport::for_supplied(&record)
        .map_err(|error| invalid_params(format!("invalid recurrence JSON: {error}")))?;
    let serde_json::Value::Object(fields) =
        serde_json::to_value(report).expect("serialize theorem report")
    else {
        unreachable!("theorem report is an object")
    };
    Ok(fields.into_iter().collect())
}

#[cfg(test)]
mod tests {
    use crate::requests::{FindRecurrenceRequest, RationalCoefficientInput};
    use crate::{PolynomialToolsServer, RecurrenceTheoremsRequest};
    use rmcp::{handler::server::wrapper::Parameters, Json};

    fn rows(text: &[&[i64]]) -> Vec<Vec<RationalCoefficientInput>> {
        text.iter()
            .map(|row| row.iter().map(|&c| c.into()).collect())
            .collect()
    }

    #[test]
    fn find_recurrence_reports_conditional_renewal_match() {
        let server = PolynomialToolsServer::new();
        // 1/(1 - x(1 + 3tx + 2t^2x^2)): R = (1+z)(1+2z), ell = 1.
        let Json(found) = server
            .find_recurrence(Parameters(FindRecurrenceRequest {
                coefficients: Some(rows(&[
                    &[1],
                    &[1],
                    &[1, 3],
                    &[1, 6, 2],
                    &[1, 9, 13],
                    &[1, 12, 33, 12],
                    &[1, 15, 62, 63, 4],
                    &[1, 18, 100, 180, 66],
                    &[1, 21, 147, 390, 321, 36],
                    &[1, 24, 203, 720, 985, 360, 8],
                ])),
                max_var_deg: Some(2),
                max_idx_deg: Some(0),
                max_rec_len: Some(3),
                ..Default::default()
            }))
            .unwrap();
        let matches = found.theorem_matches.expect("theorem match");
        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0]["theorem"], "weighted_renewal_window");
        assert_eq!(matches[0]["status"], "conditional_theorem_match");
        assert_eq!(matches[0]["window"]["max_step"], 1);
        assert_eq!(matches[0]["boundary"]["status"], "canonical_after_relabel");
        assert_eq!(matches[0]["proof_status"], "paper_proof_not_lean_checked");
        // The finite-fit fields are unchanged and separate.
        assert!(found.fit_polynomials.is_some());

        let Json(report) = server
            .recurrence_theorems(Parameters(RecurrenceTheoremsRequest {
                recurrence_json: found.recurrence_json.unwrap(),
            }))
            .unwrap();
        assert_eq!(report["schema"], "polytool.recurrence-theorems.v1");
        assert_eq!(report["finite_evidence"]["source"], "supplied_recurrence");
        assert_eq!(report["theorem_matches"][0]["window"]["max_step"], 1);
        assert_eq!(
            report["rejected_theorems"][0]["theorem"],
            "fixed_delay_predecessor_window"
        );
    }

    #[test]
    fn non_matching_and_invalid_records() {
        let server = PolynomialToolsServer::new();
        let Json(fibonacci) = server
            .find_recurrence(Parameters(FindRecurrenceRequest {
                coefficients: Some(rows(&[&[1], &[1], &[2], &[3], &[5], &[8], &[13]])),
                ..Default::default()
            }))
            .unwrap();
        assert!(fibonacci.found);
        assert!(fibonacci.theorem_matches.is_none());
        let value = serde_json::to_value(&fibonacci).unwrap();
        assert!(value.get("theorem_matches").is_none());
        let Json(report) = server
            .recurrence_theorems(Parameters(RecurrenceTheoremsRequest {
                recurrence_json: fibonacci.recurrence_json.unwrap(),
            }))
            .unwrap();
        assert_eq!(report["theorem_matches"], serde_json::json!([]));
        assert_eq!(report["rejected_theorems"].as_array().unwrap().len(), 2);
        assert!(server
            .recurrence_theorems(Parameters(RecurrenceTheoremsRequest {
                recurrence_json: "{}".into(),
            }))
            .is_err());
    }
}
