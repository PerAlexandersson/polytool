//! `polytool recurrence-implicit`: diagonal maps and derivative antecedents,
//! implicit differential relations with reconstruction data, uniqueness
//! arguments, reconstruction certificates and exact replay.

use crate::input::error_text;
use crate::input::read_stdin_or_exit;
use crate::json_out::print_json_result;
use polytool::recurrence::{
    assess_generating_description, implicit_transform_recurrence_with_root_domain,
    parse_implicit_transform_input, search_implicit_relations, verify_implicit_relation_search,
    verify_recurrence_implicit, GeneratingDescription, ImplicitRelationSearchReport,
    ImplicitRelationSearchRequest, ImplicitTransformReport,
};
use serde_json::json;

pub(crate) fn cmd_recurrence_implicit(args: &[String]) {
    let mut verify = false;
    let mut assess = false;
    let mut search = false;
    let mut verify_search = false;
    let mut compact = false;
    for arg in args {
        match arg.as_str() {
            "--verify" => verify = true,
            "--assess" => assess = true,
            "--search" => search = true,
            "--verify-search" => verify_search = true,
            "--compact" => compact = true,
            "--json" => {}
            other => {
                eprintln!("unknown option: {other}; bounds belong in the JSON request");
                std::process::exit(2);
            }
        }
    }
    if [verify, assess, search, verify_search]
        .iter()
        .filter(|flag| **flag)
        .count()
        > 1
    {
        eprintln!("give at most one of --verify, --assess, --search and --verify-search");
        std::process::exit(2);
    }
    let input = read_stdin_or_exit();
    let outcome = if search {
        serde_json::from_str::<ImplicitRelationSearchRequest>(&input)
            .map_err(error_text)
            .and_then(|request| search_implicit_relations(&request).map_err(error_text))
            .and_then(|report| serde_json::to_value(report).map_err(error_text))
    } else if verify_search {
        serde_json::from_str::<ImplicitRelationSearchReport>(&input)
            .map_err(error_text)
            .and_then(|report| {
                verify_implicit_relation_search(&report).map_err(error_text)?;
                Ok(json!({
                    "verified": true,
                    "results": report
                        .results
                        .iter()
                        .map(|r| json!({"status": r.status, "evidence": r.evidence}))
                        .collect::<Vec<_>>(),
                }))
            })
    } else if verify {
        serde_json::from_str::<ImplicitTransformReport>(&input)
            .map_err(error_text)
            .and_then(|report| {
                verify_recurrence_implicit(&report).map_err(error_text)?;
                Ok(json!({
                    "verified": true,
                    "assessment": report.assessment.status,
                    "independent_search": report
                        .independent_search
                        .as_ref()
                        .map(|s| s.evidence.clone()),
                }))
            })
    } else if assess {
        serde_json::from_str::<GeneratingDescription>(&input)
            .map_err(error_text)
            .and_then(|description| assess_generating_description(&description).map_err(error_text))
            .and_then(|assessment| serde_json::to_value(assessment).map_err(error_text))
    } else {
        serde_json::from_str::<serde_json::Value>(&input)
            .map_err(error_text)
            .and_then(|value| parse_implicit_transform_input(value).map_err(error_text))
            .and_then(|(request, root_domain)| {
                implicit_transform_recurrence_with_root_domain(&request, root_domain.as_ref())
                    .map_err(error_text)
            })
            .and_then(|report| serde_json::to_value(report).map_err(error_text))
    };
    print_json_result(outcome, compact);
}
