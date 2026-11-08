//! `polytool recurrence-transform`: exact coordinate changes of a polynomial
//! sequence, transport of its defining recurrence, and exact replay.

use crate::input::error_text;
use crate::input::read_stdin_or_exit;
use crate::json_out::print_json_result;
use polytool::recurrence::{
    transform_recurrence, verify_recurrence_transform, RecurrenceTransformReport,
    RecurrenceTransformRequest,
};
use serde_json::json;

pub(crate) fn cmd_recurrence_transform(args: &[String]) {
    let mut verify = false;
    let mut compact = false;
    for arg in args {
        match arg.as_str() {
            "--verify" => verify = true,
            "--compact" => compact = true,
            "--json" => {}
            other => {
                eprintln!("unknown option: {other}; bounds belong in the JSON request");
                std::process::exit(2);
            }
        }
    }
    let input = read_stdin_or_exit();
    let outcome = if verify {
        serde_json::from_str::<RecurrenceTransformReport>(&input)
            .map_err(error_text)
            .and_then(|report| {
                verify_recurrence_transform(&report).map_err(error_text)?;
                Ok(json!({
                    "verified": true,
                    "transported": report.transported.is_some(),
                    "independent_search": report.independent_search.is_some(),
                }))
            })
    } else {
        serde_json::from_str::<RecurrenceTransformRequest>(&input)
            .map_err(error_text)
            .and_then(|request| transform_recurrence(&request).map_err(error_text))
            .and_then(|report| serde_json::to_value(report).map_err(error_text))
    };
    print_json_result(outcome, compact);
}
