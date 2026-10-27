//! WASM exports for coordinate transformations, implicit relations, exact
//! root enclosures and conditional theorem matches.
//!
//! Every export takes and returns JSON text.  Reports are the library's
//! versioned schemas unchanged (`polytool.recurrence-transform.v1`,
//! `polytool.recurrence-implicit.v1`, `polytool.root-enclosures.v1`,
//! `polytool.recurrence-theorems.v1`), so a report produced in the browser
//! replays with the CLI and MCP verifiers and vice versa.  Failures are
//! returned as `{"error": "..."}`.
//!
//! CPU use is bounded by the library limit types: requests above the default
//! service ceilings of [`TransformLimits`] are rejected, as in the MCP
//! server, and the root-domain section is validated by the library.  The
//! browser additionally caps request size at 1 MiB and the refinement width
//! of root enclosures at `2^-256`.  The page runs these exports in a Web
//! Worker that the user can cancel.

use num_bigint::BigInt;
use polytool::recurrence::{
    ImplicitTransformReport, ImplicitTransformRequest, RecurrenceJson, RecurrenceTransformReport,
    RecurrenceTransformRequest, RootEnclosureReport, RootEnclosureRequest, TransformLimits,
};
use polytool::renewal::RecurrenceTheoremReport;
use wasm_bindgen::prelude::*;

/// Largest JSON request or report accepted by these exports.
const MAX_JSON_BYTES: usize = 1_048_576;

/// Root enclosures are refined to widths of at least `2^-MAX_WIDTH_BITS`.
const MAX_WIDTH_BITS: u64 = 256;

const TRANSFORM_SCHEMA: &str = "polytool.recurrence-transform.v1";
const IMPLICIT_SCHEMA: &str = "polytool.recurrence-implicit.v1";
const ROOT_ENCLOSURES_SCHEMA: &str = "polytool.root-enclosures.v1";
const THEOREMS_SCHEMA: &str = "polytool.recurrence-theorems.v1";
const SPACES_SCHEMA: &str = "polytool.recurrence-spaces.v1";

fn respond(run: impl FnOnce() -> Result<serde_json::Value, String>) -> String {
    match run() {
        Ok(value) => value.to_string(),
        Err(error) => serde_json::json!({ "error": error }).to_string(),
    }
}

fn bounded(json: &str, what: &str) -> Result<(), String> {
    if json.len() > MAX_JSON_BYTES {
        return Err(format!("{what} exceeds 1 MiB"));
    }
    Ok(())
}

fn to_value(value: impl serde::Serialize) -> Result<serde_json::Value, String> {
    serde_json::to_value(value).map_err(|e| e.to_string())
}

fn check_schema(found: &str, expected: &str) -> Result<(), String> {
    if found != expected {
        return Err(format!("expected schema {expected}, found {found}"));
    }
    Ok(())
}

fn check_limits(limits: &TransformLimits) -> Result<(), String> {
    limits.check_ceiling(&TransformLimits::default())
}

/// Reject refinement widths below `2^-MAX_WIDTH_BITS`; the library rejects
/// non-positive and malformed widths.
fn check_enclosure_request(request: &RootEnclosureRequest) -> Result<(), String> {
    let Some(width) = &request.max_width else {
        return Ok(());
    };
    let width = polytool::recurrence::parse_rational_coeff(width)?;
    let numerator: BigInt = width.numer().clone();
    if numerator > BigInt::from(0) && width.denom().bits() > numerator.bits() + MAX_WIDTH_BITS {
        return Err(format!(
            "max_width must be at least 2^-{MAX_WIDTH_BITS} in the browser"
        ));
    }
    Ok(())
}

/// Index of the bundled OEIS families: `[{id, name, first_row}]`.
#[wasm_bindgen]
pub fn oeis_catalog_index() -> String {
    let entries: Vec<serde_json::Value> = polytool::oeis::catalog()
        .iter()
        .map(|entry| {
            serde_json::json!({
                "id": entry.id,
                "name": entry.name,
                "first_row": entry.first_row,
            })
        })
        .collect();
    serde_json::Value::Array(entries).to_string()
}

/// Same `polytool.recurrence-transform.v1` contract as the CLI
/// `recurrence-transform` and the MCP tool `transform_recurrence`.
#[wasm_bindgen]
pub fn transform_recurrence(request_json: &str) -> String {
    respond(|| {
        bounded(request_json, "recurrence-transform request")?;
        let request: RecurrenceTransformRequest =
            serde_json::from_str(request_json).map_err(|e| e.to_string())?;
        check_limits(&request.limits)?;
        to_value(polytool::recurrence::transform_recurrence(&request)?)
    })
}

fn verify_transform_value(report_json: &str) -> Result<serde_json::Value, String> {
    bounded(report_json, "recurrence-transform report")?;
    let report: RecurrenceTransformReport =
        serde_json::from_str(report_json).map_err(|e| e.to_string())?;
    check_schema(&report.schema, TRANSFORM_SCHEMA)?;
    check_limits(&report.request.limits)?;
    polytool::recurrence::verify_recurrence_transform(&report)?;
    Ok(serde_json::json!({
        "schema": TRANSFORM_SCHEMA,
        "verified": true,
        "transported": report.transported.is_some(),
        "independent_search": report.independent_search.is_some(),
        "root_domain": report.root_domain.is_some(),
    }))
}

/// Replay a transform report exactly from its embedded request.
#[wasm_bindgen]
pub fn verify_recurrence_transform(report_json: &str) -> String {
    respond(|| verify_transform_value(report_json))
}

/// Same `polytool.recurrence-implicit.v1` contract as the CLI
/// `recurrence-implicit` and the MCP tool `transform_recurrence_implicit`.
#[wasm_bindgen]
pub fn transform_recurrence_implicit(request_json: &str) -> String {
    respond(|| {
        bounded(request_json, "recurrence-implicit request")?;
        let request: ImplicitTransformRequest =
            serde_json::from_str(request_json).map_err(|e| e.to_string())?;
        check_limits(&request.limits)?;
        to_value(polytool::recurrence::implicit_transform_recurrence(
            &request,
        )?)
    })
}

fn verify_implicit_value(report_json: &str) -> Result<serde_json::Value, String> {
    bounded(report_json, "recurrence-implicit report")?;
    let report: ImplicitTransformReport =
        serde_json::from_str(report_json).map_err(|e| e.to_string())?;
    check_schema(&report.schema, IMPLICIT_SCHEMA)?;
    check_limits(&report.request.limits)?;
    polytool::recurrence::verify_recurrence_implicit(&report)?;
    Ok(serde_json::json!({
        "schema": IMPLICIT_SCHEMA,
        "verified": true,
        "assessment": report.assessment.status,
        "independent_search": report.independent_search.map(|search| search.evidence),
    }))
}

/// Replay an implicit-transport report exactly, re-checking every
/// reconstruction certificate.
#[wasm_bindgen]
pub fn verify_recurrence_implicit(report_json: &str) -> String {
    respond(|| verify_implicit_value(report_json))
}

/// Same `polytool.root-enclosures.v1` contract as the CLI `root-enclosures`
/// and the MCP tool `exact_root_enclosures`.
#[wasm_bindgen]
pub fn root_enclosures(request_json: &str) -> String {
    respond(|| {
        bounded(request_json, "root-enclosures request")?;
        let request: RootEnclosureRequest =
            serde_json::from_str(request_json).map_err(|e| e.to_string())?;
        check_enclosure_request(&request)?;
        to_value(polytool::recurrence::root_enclosures(&request)?)
    })
}

fn verify_enclosures_value(report_json: &str) -> Result<serde_json::Value, String> {
    bounded(report_json, "root-enclosures report")?;
    let report: RootEnclosureReport =
        serde_json::from_str(report_json).map_err(|e| e.to_string())?;
    check_schema(&report.schema, ROOT_ENCLOSURES_SCHEMA)?;
    check_enclosure_request(&report.request)?;
    let recomputed = polytool::recurrence::root_enclosures(&report.request)?;
    if recomputed != report {
        return Err("root enclosures differ from the exact recomputation".into());
    }
    Ok(serde_json::json!({
        "schema": ROOT_ENCLOSURES_SCHEMA,
        "verified": true,
        "polynomials": report.items.len(),
    }))
}

/// Recompute every enclosure and interval verdict from the embedded request
/// and require equality.
#[wasm_bindgen]
pub fn verify_root_enclosures(report_json: &str) -> String {
    respond(|| verify_enclosures_value(report_json))
}

fn parse_record(recurrence_json: &str) -> Result<RecurrenceJson, String> {
    bounded(recurrence_json, "recurrence JSON")?;
    serde_json::from_str(recurrence_json).map_err(|e| format!("invalid recurrence JSON: {e}"))
}

/// Same `polytool.recurrence-theorems.v1` contract as the CLI
/// `recurrence-theorems` and the MCP tool `recurrence_theorems`, for a
/// supplied `polytool.recurrence.v1` record.
#[wasm_bindgen]
pub fn recurrence_theorems(recurrence_json: &str) -> String {
    respond(|| {
        let record = parse_record(recurrence_json)?;
        to_value(
            RecurrenceTheoremReport::for_supplied(&record)
                .map_err(|e| format!("invalid recurrence JSON: {e}"))?,
        )
    })
}

fn verify_theorems_value(report_json: &str) -> Result<serde_json::Value, String> {
    bounded(report_json, "recurrence-theorems report")?;
    let report: RecurrenceTheoremReport =
        serde_json::from_str(report_json).map_err(|e| e.to_string())?;
    check_schema(&report.schema, THEOREMS_SCHEMA)?;
    if report.finite_evidence.source != "supplied_recurrence" {
        return Err(format!(
            "only supplied_recurrence reports can be replayed; a {} report needs the fitted rows",
            report.finite_evidence.source
        ));
    }
    let recomputed = RecurrenceTheoremReport::for_supplied(&report.recurrence_json)?;
    if recomputed != report {
        return Err("theorem report differs from the exact recomputation".into());
    }
    Ok(serde_json::json!({
        "schema": THEOREMS_SCHEMA,
        "verified": true,
        "theorem_matches": report.theorem_matches.len(),
    }))
}

/// Recompute the matches of a supplied-recurrence theorem report from its
/// embedded record and require equality.
#[wasm_bindgen]
pub fn verify_recurrence_theorems(report_json: &str) -> String {
    respond(|| verify_theorems_value(report_json))
}

/// Plain text and LaTeX of a `polytool.recurrence.v1` record:
/// `{"text": ..., "latex": ...}`.
#[wasm_bindgen]
pub fn recurrence_latex(recurrence_json: &str) -> String {
    respond(|| {
        let record = parse_record(recurrence_json)?;
        let (recurrence, _, _) = record.to_recurrence_parts()?;
        Ok(serde_json::json!({
            "text": recurrence.to_string(),
            "latex": recurrence.to_latex(),
        }))
    })
}

/// Replay any supported report, dispatching on its `schema` field.
///
/// Returns `{"schema", "verified": true, ...}` or `{"error"}`; an unknown
/// schema is an error, never a silent success.
#[wasm_bindgen]
pub fn verify_report(report_json: &str) -> String {
    respond(|| {
        bounded(report_json, "report")?;
        let value: serde_json::Value =
            serde_json::from_str(report_json).map_err(|e| e.to_string())?;
        let schema = value
            .get("schema")
            .and_then(serde_json::Value::as_str)
            .ok_or("the report has no schema field")?;
        match schema {
            TRANSFORM_SCHEMA => verify_transform_value(report_json),
            IMPLICIT_SCHEMA => verify_implicit_value(report_json),
            ROOT_ENCLOSURES_SCHEMA => verify_enclosures_value(report_json),
            THEOREMS_SCHEMA => verify_theorems_value(report_json),
            SPACES_SCHEMA => {
                let checked: serde_json::Value =
                    serde_json::from_str(&crate::verify_recurrence_spaces(report_json))
                        .map_err(|e| e.to_string())?;
                if let Some(error) = checked.get("error").and_then(serde_json::Value::as_str) {
                    return Err(error.to_string());
                }
                Ok(serde_json::json!({
                    "schema": SPACES_SCHEMA,
                    "verified": true,
                    "incomplete_shapes": checked["incomplete_shapes"],
                }))
            }
            other => Err(format!("no browser replay for schema {other}")),
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::{json, Value};

    fn parse(text: &str) -> Value {
        serde_json::from_str(text).expect("export returns JSON")
    }

    fn eulerian_f_request() -> Value {
        json!({
            "oeis_id": "A008292", "row_count": 8, "search": {}, "root_domain": {},
            "transformation": {"kind": "f_polynomial",
                "schedule": {"first_index": 1, "first_degree": 0, "step": 1}}
        })
    }

    #[test]
    fn transform_report_round_trips_and_rejects_tampering() {
        let report = transform_recurrence(&eulerian_f_request().to_string());
        let value = parse(&report);
        assert_eq!(value["schema"], TRANSFORM_SCHEMA);
        assert_eq!(value["transported"]["evidence"], "transported_identity");
        assert_eq!(value["independent_search"]["evidence"], "finite_fit");
        let claims = value["root_domain"]["claims"].as_array().unwrap();
        assert!(claims
            .iter()
            .any(|c| c["evidence"] == "certified_family_invariant"));
        let checked = parse(&verify_recurrence_transform(&report));
        assert_eq!(checked["verified"], true, "{checked}");
        assert_eq!(parse(&verify_report(&report))["verified"], true);

        let mut tampered = value.clone();
        tampered["transformed_states"][3][0][1] = json!("9");
        assert!(parse(&verify_report(&tampered.to_string()))["error"].is_string());
        let mut relabelled = value;
        relabelled["root_domain"]["claims"][0]["evidence"] = json!("certified_family_invariant");
        assert!(parse(&verify_recurrence_transform(&relabelled.to_string()))["error"].is_string());
    }

    #[test]
    fn transform_limits_and_size_are_bounded() {
        let mut request = eulerian_f_request();
        request["limits"] = json!({"max_rows": 1_000_000});
        let error = parse(&transform_recurrence(&request.to_string()));
        assert!(error["error"].as_str().unwrap().contains("service ceiling"));
        let huge = " ".repeat(MAX_JSON_BYTES + 1);
        assert!(parse(&transform_recurrence(&huge))["error"]
            .as_str()
            .unwrap()
            .contains("1 MiB"));
    }

    #[test]
    fn implicit_report_round_trips_and_rejects_tampering() {
        let request = json!({
            "oeis_id": "A007318", "row_count": 8, "search": {},
            "transformation": {"kind": "diagonal",
                "weight": {"kind": "factorial_power", "power": 1}}
        });
        let report = transform_recurrence_implicit(&request.to_string());
        let value = parse(&report);
        assert_eq!(value["schema"], IMPLICIT_SCHEMA);
        assert_eq!(value["assessment"]["status"], "equivalent_description");
        assert_eq!(parse(&verify_report(&report))["verified"], true);
        let mut tampered = value;
        tampered["assessment"]["status"] = json!("certified_consequence");
        assert!(parse(&verify_recurrence_implicit(&tampered.to_string()))["error"].is_string());
    }

    #[test]
    fn root_enclosures_replay_and_width_floor() {
        let request = json!({
            "polynomials": [["0", "0", "1", "3", "3", "1"], ["1", "3", "1"]],
            "max_width": "1/1000",
            "interval": {"lower": {"value": "-1"}, "upper": {"value": "0"}}
        });
        let report = root_enclosures(&request.to_string());
        let value = parse(&report);
        assert_eq!(value["items"][0]["real_root_count"], 5);
        assert_eq!(value["items"][1]["real_rooted_in_interval"], false);
        assert_eq!(parse(&verify_report(&report))["verified"], true);
        let mut tampered = value;
        tampered["items"][1]["roots"][0]["lower"] = json!("-3");
        assert!(parse(&verify_root_enclosures(&tampered.to_string()))["error"].is_string());

        let mut tiny = request;
        tiny["max_width"] = json!(format!("1/{}", BigInt::from(2).pow(300)));
        assert!(parse(&root_enclosures(&tiny.to_string()))["error"]
            .as_str()
            .unwrap()
            .contains("2^-256"));
    }

    #[test]
    fn theorem_report_replays_supplied_records_only() {
        // P_n = (1 + t) P_(n-1): a weighted renewal recurrence with ell = 1.
        let record = json!({
            "schema": "polytool.recurrence.v1", "first_index": 1,
            "initial_polynomials": [["1", "1"]],
            "recurrence": {"terms": [{"offset": 1, "deriv_order": 0, "sign": "none",
                "coeff": {"coeffs": [["1", "1"]]}}]}
        });
        let report = recurrence_theorems(&record.to_string());
        let value = parse(&report);
        assert_eq!(value["schema"], THEOREMS_SCHEMA, "{value}");
        assert_eq!(value["finite_evidence"]["source"], "supplied_recurrence");
        assert_eq!(parse(&verify_report(&report))["verified"], true);
        let mut tampered = value.clone();
        tampered["finite_evidence"]["rows_compared"] = json!(99);
        assert!(parse(&verify_recurrence_theorems(&tampered.to_string()))["error"].is_string());
        let mut searched = value;
        searched["finite_evidence"]["source"] = json!("adaptive_search_fit");
        assert!(
            parse(&verify_recurrence_theorems(&searched.to_string()))["error"]
                .as_str()
                .unwrap()
                .contains("fitted rows")
        );
        assert!(parse(&recurrence_theorems("{}"))["error"].is_string());

        let latex = parse(&recurrence_latex(&record.to_string()));
        assert!(latex["latex"].as_str().unwrap().contains("P(n)"));
    }

    #[test]
    fn verify_report_rejects_unknown_schemas_and_catalog_lists_families() {
        assert!(
            parse(&verify_report(r#"{"schema":"polytool.unknown.v1"}"#))["error"]
                .as_str()
                .unwrap()
                .contains("no browser replay")
        );
        assert!(parse(&verify_report("{}"))["error"].is_string());
        let catalog = parse(&oeis_catalog_index());
        let entries = catalog.as_array().unwrap();
        assert!(entries.iter().any(|e| e["id"] == "A008292"));
        assert!(entries.len() > 100);
    }
}
