//! Shared-parser stability adapter and compact report storage.

use crate::{limits::*, parsing::invalid_params, report_store::ReportStore, requests::Detail};
use polytool::stability::{
    check_univariate_stability_with_convention,
    parse::{parse_univariate, ParseLimits},
    wire,
};
use rmcp::ErrorData;
use schemars::{JsonSchema, Schema, SchemaGenerator};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{
    borrow::Cow,
    time::{Duration, Instant},
};

/// Exact stability batch input. Coefficient arrays use ascending degree order.
#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct StabilityUnivariateRequest {
    /// Expressions or ascending Gaussian coefficient arrays; 1 through 512 items.
    pub polynomials: Vec<Value>,
    /// Named region, default upper; exclusive with other selectors.
    pub region: Option<String>,
    /// Four Gaussian entries of a nonsingular fractional linear map.
    pub mobius: Option<Vec<Value>>,
    /// Exact Gaussian disk center; requires radius.
    pub disk_center: Option<Value>,
    /// Strictly positive exact rational radius string.
    pub radius: Option<String>,
    /// Override the region's default: avoid or inside.
    pub convention: Option<String>,
    /// Summary by default; full returns the entire replay certificate.
    #[serde(default)]
    pub detail: Detail,
}

impl JsonSchema for StabilityUnivariateRequest {
    fn schema_name() -> Cow<'static, str> {
        "StabilityUnivariateRequest".into()
    }
    fn json_schema(_: &mut SchemaGenerator) -> Schema {
        Schema::try_from(json!({
            "type":"object","additionalProperties":false,"required":["polynomials"],
            "properties":{
                "polynomials":{"type":"array","minItems":1,"maxItems":512,"items":{"oneOf":[{"type":"string"},{"type":"array","items":{"$ref":"#/$defs/gaussian"}}]}},
                "region":{"type":"string","enum":["upper","lower","left","right","disk","exterior"],"default":"upper"},
                "mobius":{"type":"array","minItems":4,"maxItems":4,"items":{"$ref":"#/$defs/gaussian"}},
                "disk_center":{"$ref":"#/$defs/gaussian"},"radius":{"type":"string"},
                "convention":{"type":"string","enum":["avoid","inside"]},
                "detail":{"type":"string","enum":["summary","full"],"default":"summary"}
            },
            "$defs":{"gaussian":{"oneOf":[{"type":"integer"},{"type":"string"},{"type":"object","additionalProperties":false,"required":["real","imaginary"],"properties":{"real":{"type":"string"},"imaginary":{"type":"string"}}}]}}
        })).expect("valid stability schema")
    }
}

pub(crate) const HELP: &str = "stability_univariate accepts expressions with integer/rational parts and I, or ascending Gaussian coefficient arrays. Gaussian values are JSON integers, exact strings, or {real,imaginary} rational strings; JSON floats fail. Choose region (upper default), mobius [a,b,c,d], or disk_center and positive rational radius, exclusively. Upper/lower/mobius default to avoid; other regions default to inside. Counts include multiplicities. Zero has null degree/counts and passes only avoid; constants pass both. Limits: 512 polynomials, degree 500, 1000000 bits per numerator/denominator, existing 1 MiB text, 16384 bytes per coefficient literal, 65536 total coefficients. A 50-second elapsed-time guard rejects the whole response at computation checkpoints; it does not interrupt an individual exact arithmetic operation. No partial counts are returned. detail:full includes canonical input/options and certificates in polytool.stability.v1. If configured, reports are stored; use report_id with read_report to inspect evidence. CLI: polytool stability verify --input report.json replays the complete report.";

pub(crate) fn response(
    input: StabilityUnivariateRequest,
    store: &ReportStore,
) -> Result<Value, ErrorData> {
    let start = Instant::now();
    compute(input, store, || start.elapsed() >= Duration::from_secs(50))
}

fn compute(
    input: StabilityUnivariateRequest,
    store: &ReportStore,
    expired: impl Fn() -> bool,
) -> Result<Value, ErrorData> {
    let limits = ParseLimits {
        max_degree: MCP_MAX_STABILITY_DEGREE,
        max_coefficient_bits: MCP_MAX_STABILITY_COEFFICIENT_BITS,
        max_literal_bytes: MCP_MAX_COEFFICIENT_TEXT_BYTES,
    };
    let invalid = |e: polytool::stability::parse::ParseError| invalid_params(e.to_string());
    if input.polynomials.is_empty() || input.polynomials.len() > MCP_MAX_BATCH_POLYNOMIALS {
        return Err(invalid_params("expected 1 through 512 polynomials"));
    }
    let text = serde_json::to_vec(&input).map_err(|e| invalid_params(e.to_string()))?;
    if text.len() > MCP_MAX_TEXT_BYTES {
        return Err(invalid_params("stability request exceeds 1 MiB text limit"));
    }
    let mobius = input
        .mobius
        .as_ref()
        .map(|m| {
            m.iter()
                .map(|c| wire::parse_gaussian_json(c, limits))
                .collect::<Result<Vec<_>, _>>()
        })
        .transpose()
        .map_err(invalid)?;
    let center = input
        .disk_center
        .as_ref()
        .map(|c| wire::parse_gaussian_json(c, limits))
        .transpose()
        .map_err(invalid)?;
    let center_text = center
        .as_ref()
        .map(|c| format!("({})+({})*I", c.real, c.imaginary));
    let region = wire::parse_region(
        input.region.as_deref(),
        mobius.as_deref(),
        center_text.as_deref(),
        input.radius.as_deref(),
        limits,
    )
    .map_err(invalid)?;
    let convention =
        wire::parse_convention(input.convention.as_deref(), &region).map_err(invalid)?;
    let mut parsed = Vec::new();
    let mut total = 0usize;
    for value in &input.polynomials {
        if expired() {
            return Err(invalid_params(
                "stability time guard exceeded; no exact result returned",
            ));
        }
        let coefficients = match value {
            Value::String(s) => parse_univariate(s, limits).map_err(invalid)?,
            Value::Array(a) => {
                if a.len() > limits.max_degree + 1 {
                    return Err(invalid_params("coefficient array exceeds degree 500"));
                }
                for c in a {
                    if c.to_string().len() > MCP_MAX_COEFFICIENT_TEXT_BYTES {
                        return Err(invalid_params("coefficient literal exceeds text limit"));
                    }
                }
                a.iter()
                    .map(|c| wire::parse_gaussian_json(c, limits))
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(invalid)?
            }
            _ => {
                return Err(invalid_params(
                    "polynomial must be expression or ascending coefficient array",
                ))
            }
        };
        limits.validate(&coefficients).map_err(invalid)?;
        total += coefficients.len();
        if total > MCP_MAX_TOTAL_COEFFICIENTS {
            return Err(invalid_params(
                "stability batch exceeds aggregate coefficient limit",
            ));
        }
        parsed.push(coefficients);
    }
    let mut items = Vec::new();
    let mut full = Vec::new();
    for (index, coefficients) in parsed.iter().enumerate() {
        if expired() {
            return Err(invalid_params(
                "stability time guard exceeded; no exact result returned",
            ));
        }
        let report = check_univariate_stability_with_convention(coefficients, &region, convention)
            .map_err(|e| invalid_params(e.to_string()))?;
        if expired() {
            return Err(invalid_params(
                "stability time guard exceeded; no exact result returned",
            ));
        }
        items.push(wire::summary(index, &report));
        full.push(wire::full_item(index, coefficients, &report));
    }
    let full = wire::full_report(full);
    let bytes = serde_json::to_vec(&full).expect("JSON report");
    if expired() {
        return Err(invalid_params(
            "stability time guard exceeded; no exact result returned",
        ));
    }
    let handle = store.store(&bytes);
    let mut result = if input.detail == Detail::Full {
        full
    } else {
        json!({"items":items})
    };
    if let Some(handle) = handle {
        if let Some(id) = handle.get("id") {
            result["report_id"] = id.clone();
        }
        result["report"] = handle;
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn request(value: Value) -> StabilityUnivariateRequest {
        serde_json::from_value(value).unwrap()
    }
    #[test]
    fn summary_full_and_budgets() {
        let store = ReportStore::new(Default::default());
        let r = response(
            request(json!({"polynomials":["(t+1)^3","t-I",["I",1],"0",[]],"region":"left"})),
            &store,
        )
        .unwrap();
        assert_eq!(
            r["items"][0]["counts"],
            json!({"inside":3,"boundary":0,"outside":0})
        );
        assert_eq!(r["items"][3]["degree"], Value::Null);
        assert_eq!(r["items"][3]["counts"], Value::Null);
        for item in r["items"].as_array().unwrap() {
            assert!(serde_json::to_vec(item).unwrap().len() <= 200);
            assert!(item.get("certificate").is_none());
        }
        let full = response(
            request(json!({"polynomials":["t^2+(1+2I)t-I"],"detail":"full"})),
            &store,
        )
        .unwrap();
        let limits = ParseLimits {
            max_degree: 500,
            max_coefficient_bits: 1_000_000,
            max_literal_bytes: 16384,
        };
        assert!(wire::verify_report(&full, limits).unwrap());
        let many = response(request(json!({"polynomials":vec!["t+1";512]})), &store).unwrap();
        assert!(serde_json::to_vec(&many).unwrap().len() < 512 * 200 + 32);
    }
    #[test]
    fn validates_selectors_exactness_and_caps() {
        let store = ReportStore::new(Default::default());
        for value in [
            json!({"polynomials":[]}),
            json!({"polynomials":vec!["1";513]}),
            json!({"polynomials":[[1.0]]}),
            json!({"polynomials":["t^501"]}),
            json!({"polynomials":["x+y"]}),
            json!({"polynomials":["1/(t+1)"]}),
            json!({"polynomials":["1.1"]}),
            json!({"polynomials":["t"],"radius":"1"}),
            json!({"polynomials":["t"],"disk_center":0,"radius":"0"}),
            json!({"polynomials":["t"],"region":"upper","mobius":[1,0,0,1]}),
            json!({"polynomials":["t"],"mobius":[1,1,1,1]}),
            json!({"polynomials":["t"],"convention":"other"}),
            json!({"polynomials":[["1".repeat(16385)]]}),
            json!({"polynomials":[format!("{}+t","1".repeat(16385))]}),
            json!({"polynomials":vec![vec![1;501];131]}),
            json!({"polynomials":[" ".repeat(MCP_MAX_TEXT_BYTES+1)]}),
        ] {
            assert!(response(request(value.clone()), &store).is_err(), "{value}");
        }
        assert!(serde_json::from_value::<StabilityUnivariateRequest>(
            json!({"polynomials":["1"],"unknown":1})
        )
        .is_err());
        assert!(compute(request(json!({"polynomials":["t"]})), &store, || true).is_err());
    }
    #[test]
    fn expiry_after_a_completed_item_leaves_no_partial_report() {
        let temp = crate::report_store::tests::TempDirectory::new();
        let store = ReportStore::new(temp.config());
        let checkpoints = std::cell::Cell::new(0);
        let result = compute(
            request(json!({"polynomials":["t+1","t+2"]})),
            &store,
            || {
                checkpoints.set(checkpoints.get() + 1);
                checkpoints.get() >= 5
            },
        );
        assert!(result.is_err());
        assert_eq!(std::fs::read_dir(&temp.0).unwrap().count(), 0);
    }
    #[test]
    fn stored_report_and_schema() {
        let temp = crate::report_store::tests::TempDirectory::new();
        let store = ReportStore::new(temp.config());
        let r = response(
            request(json!({"polynomials":["t-I"],"mobius":[1,0,0,1]})),
            &store,
        )
        .unwrap();
        let bytes = store.read(r["report_id"].as_str().unwrap()).unwrap();
        let full: Value = serde_json::from_slice(&bytes).unwrap();
        assert!(full["items"][0].get("certificate").is_some());
        assert!(serde_json::to_vec(&r).unwrap().len() < 700);
        let schema =
            serde_json::to_value(schemars::schema_for!(StabilityUnivariateRequest)).unwrap();
        assert_eq!(schema["additionalProperties"], false);
        assert_eq!(schema["properties"]["polynomials"]["maxItems"], 512);
        assert_eq!(schema["$defs"]["gaussian"]["oneOf"][0]["type"], "integer");
    }
}
