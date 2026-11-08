//! `polytool recurrence-search-representations`: the bounded queue search
//! over equivalent representations (coordinate changes, transported and
//! certified relations, fits) with exact replay.

use crate::input::error_text;
use crate::input::read_stdin_or_exit;
use crate::json_out::print_json_result;
use polytool::recurrence::{
    search_representations, verify_representation_search, RepresentationSearchReport,
    RepresentationSearchRequest, EQUIVALENT_DESCRIPTION_LABEL,
};
use serde_json::json;

pub(crate) fn cmd_recurrence_search_representations(args: &[String]) {
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
        serde_json::from_str::<RepresentationSearchReport>(&input)
            .map_err(error_text)
            .and_then(|report| {
                verify_representation_search(&report).map_err(error_text)?;
                let equivalent = report
                    .representations
                    .iter()
                    .flat_map(|n| &n.descriptions)
                    .filter(|d| d.label == EQUIVALENT_DESCRIPTION_LABEL)
                    .count();
                Ok(json!({
                    "verified": true,
                    "status": report.status,
                    "representations": report.representations.len(),
                    "equivalent_descriptions": equivalent,
                }))
            })
    } else {
        serde_json::from_str::<RepresentationSearchRequest>(&input)
            .map_err(error_text)
            .and_then(|request| search_representations(&request).map_err(error_text))
            .and_then(|report| serde_json::to_value(report).map_err(error_text))
    };
    print_json_result(outcome, compact);
}
