//! Bounded family catalog adapters; no filesystem or network writes.
use crate::parsing::invalid_params;
use polytool::families::*;
use rmcp::ErrorData;
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{json, Value};

#[derive(Clone, Debug, Default, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ListFamiliesRequest {
    /// Case-insensitive ID/title substring, at most 256 bytes.
    pub query: Option<String>,
    /// Real-rootedness status: unknown, conjectured, finite_evidence,
    /// lean_proved, theorem_proved, known_not_real_rooted.
    pub status: Option<String>,
    /// Page offset, default 0.
    pub offset: Option<usize>,
    /// Page size 1..100, default 50. Summaries omit full generators.
    pub limit: Option<usize>,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GetFamilyRequest {
    /// Bundled A###### or PT-F#### ID; case-insensitive.
    pub id: String,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenerateFamilyRowsRequest {
    /// Bundled A###### or PT-F#### ID; case-insensitive.
    pub id: String,
    /// Complete rows from first_row, 1..64, default 8 (12 for export).
    pub count: Option<usize>,
}
fn find(id: &str) -> Result<FamilyRecord, ErrorData> {
    if id.len() > 16 {
        return Err(invalid_params("family ID too long"));
    }
    catalog()
        .map_err(|e| invalid_params(e.to_string()))?
        .into_iter()
        .find(|r| r.id.eq_ignore_ascii_case(id))
        .ok_or_else(|| invalid_params("unknown family ID"))
}
fn bounded(value: Value) -> Result<Value, ErrorData> {
    if value.to_string().len() > 1024 * 1024 {
        Err(invalid_params(
            "family response exceeds 1 MiB; request fewer rows",
        ))
    } else {
        Ok(value)
    }
}
fn count(input: &GenerateFamilyRowsRequest, default: usize) -> Result<usize, ErrorData> {
    let n = input.count.unwrap_or(default);
    if !(1..=64).contains(&n) {
        Err(invalid_params("count must be 1..64"))
    } else {
        Ok(n)
    }
}
pub(crate) fn list(input: ListFamiliesRequest) -> Result<Value, ErrorData> {
    let query = input.query.unwrap_or_default().to_lowercase();
    if query.len() > 256 {
        return Err(invalid_params("query exceeds 256 bytes"));
    }
    if let Some(status) = &input.status {
        serde_json::from_value::<FamilyProofStatus>(json!(status))
            .map_err(|_| invalid_params("unknown family proof status"))?;
    }
    let limit = input.limit.unwrap_or(50);
    if !(1..=100).contains(&limit) {
        return Err(invalid_params("limit must be 1..100"));
    }
    let records = catalog().map_err(|e| invalid_params(e.to_string()))?;
    let entries: Vec<_> = records
        .iter()
        .filter(|r| r.id.to_lowercase().contains(&query) || r.title.to_lowercase().contains(&query))
        .filter(|r| {
            input
                .status
                .as_ref()
                .is_none_or(|s| json!(r.real_rootedness.status) == *s)
        })
        .collect();
    let total = entries.len();
    let offset = input.offset.unwrap_or(0);
    let page:Vec<_>=entries.into_iter().skip(offset).take(limit).map(|r|json!({"id":r.id,"title":r.title,"status":r.real_rootedness.status,"interlacing":r.interlacing.status,"degree_pattern":r.degree_pattern,"lean":r.real_rootedness.lean.as_ref().map(|l|&l.permalink),"models":r.models})).collect();
    let end = offset.saturating_add(page.len());
    Ok(
        json!({"schema":"polytool.families.v1","families":page,"total":total,"next_offset":if end<total {Some(end)} else {None}}),
    )
}
pub(crate) fn show(input: GetFamilyRequest) -> Result<Value, ErrorData> {
    bounded(json!({"schema":"polytool.family.v1","family":find(&input.id)?}))
}
pub(crate) fn generate(input: GenerateFamilyRowsRequest) -> Result<Value, ErrorData> {
    let n = count(&input, 8)?;
    let record = find(&input.id)?;
    let rows = record
        .generate_rows(n)
        .map_err(|e| invalid_params(e.to_string()))?;
    bounded(
        json!({"schema":"polytool.family-rows.v1","id":record.id,"first_row":record.first_row,"rows":rows}),
    )
}
pub(crate) fn export(input: GenerateFamilyRowsRequest) -> Result<Value, ErrorData> {
    let n = count(&input, 12)?;
    bounded(
        serde_json::to_value(
            export_oeis(&find(&input.id)?, n).map_err(|e| invalid_params(e.to_string()))?,
        )
        .unwrap(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use rmcp::handler::server::wrapper::Parameters;
    #[test]
    fn tools_share_generators_and_reject_excessive_requests() {
        let server = crate::PolynomialToolsServer::new();
        let page = server
            .list_families(Parameters(ListFamiliesRequest {
                query: Some("PT-F".into()),
                limit: Some(1),
                ..Default::default()
            }))
            .unwrap()
            .0;
        assert_eq!(page["total"], 2);
        assert_eq!(page["next_offset"], 1);
        let row = server
            .generate_family_rows(Parameters(GenerateFamilyRowsRequest {
                id: "pt-f0001".into(),
                count: Some(5),
            }))
            .unwrap()
            .0;
        assert_eq!(row["rows"][4], json!(["1", "53", "3"]));
        assert!(server
            .get_family(Parameters(GetFamilyRequest {
                id: "A000000".into()
            }))
            .is_err());
        assert!(server
            .generate_family_rows(Parameters(GenerateFamilyRowsRequest {
                id: "PT-F0001".into(),
                count: Some(65)
            }))
            .is_err());
        let export = server
            .export_family_oeis(Parameters(GenerateFamilyRowsRequest {
                id: "PT-F0002".into(),
                count: Some(6),
            }))
            .unwrap()
            .0;
        assert!(export["blocks"].as_str().unwrap().contains("%S A000000"));
        assert!(server
            .list_families(Parameters(ListFamiliesRequest {
                limit: Some(101),
                ..Default::default()
            }))
            .is_err());
        assert!(list(ListFamiliesRequest {
            status: Some("proved".into()),
            ..Default::default()
        })
        .is_err());
    }
}
