//! Shared bounded recurrence input for replay and theorem matching.
//!
//! Extracting a recurrence from a report does not verify the report's claims.

use crate::input::read_limited_text;
use polytool::recurrence::RecurrenceJson;
use serde_json::Value;
use std::fs;

pub(crate) fn read_recurrence_json(
    path: &str,
    candidate: Option<usize>,
) -> Result<RecurrenceJson, String> {
    let input = if path == "-" {
        crate::input::read_command_input("recurrence JSON from stdin")?
    } else {
        let file = fs::File::open(path)
            .map_err(|e| format!("failed to open recurrence JSON `{path}`: {e}"))?;
        read_limited_text(file, &format!("recurrence JSON `{path}`"))?
    };
    let value: Value = serde_json::from_str(&input)
        .map_err(|e| format!("failed to parse recurrence JSON: {e}"))?;
    let record = match value.get("schema").and_then(Value::as_str) {
        Some("polytool.recurrence-representations.v1") => {
            let items = value["representations"].as_array()
                .ok_or("recurrence menu needs a representations array")?;
            let index = match (candidate, items.len()) {
                (Some(index), _) => index,
                (None, 1) => 0,
                (None, 0) => return Err("recurrence menu contains no candidates".into()),
                (None, n) => return Err(format!(
                    "recurrence menu contains {n} candidates; choose --candidate <index> (0-based)"
                )),
            };
            items.get(index).ok_or_else(|| format!(
                "candidate {index} is out of range for {} candidates", items.len()
            ))?.get("recurrence_json").ok_or("candidate has no recurrence_json")?.clone()
        }
        Some("polytool.recurrence-theorems.v1") => {
            if candidate.is_some() {
                return Err("--candidate is only valid for a recurrence representation menu".into());
            }
            value.get("recurrence_json").ok_or("theorem report has no recurrence_json")?.clone()
        }
        Some("polytool.recurrence.v1") => {
            if candidate.is_some() {
                return Err("--candidate is only valid for a recurrence representation menu".into());
            }
            value
        }
        _ => return Err("expected polytool.recurrence.v1, polytool.recurrence-theorems.v1 or polytool.recurrence-representations.v1 JSON".into()),
    };
    if record.get("schema").and_then(Value::as_str) != Some("polytool.recurrence.v1") {
        return Err("selected record must have schema polytool.recurrence.v1".into());
    }
    serde_json::from_value(record).map_err(|e| format!("failed to parse recurrence JSON: {e}"))
}
