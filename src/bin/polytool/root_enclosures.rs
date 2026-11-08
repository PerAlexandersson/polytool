//! `polytool root-enclosures`: exact rational root enclosures with
//! multiplicities and exact interval tests.

use crate::input::error_text;
use crate::input::read_stdin_or_exit;
use crate::json_out::print_json_result;
use polytool::recurrence::{
    root_enclosures, verify_root_enclosures, RootEnclosureReport, RootEnclosureRequest,
    ROOT_ENCLOSURES_SCHEMA,
};

pub(crate) fn cmd_root_enclosures(args: &[String]) {
    let mut compact = false;
    let mut verify = false;
    for arg in args {
        match arg.as_str() {
            "--compact" => compact = true,
            "--verify" => verify = true,
            "--json" => {}
            other => {
                eprintln!("unknown option: {other}; parameters belong in the JSON request");
                std::process::exit(2);
            }
        }
    }
    let input = read_stdin_or_exit();
    if verify {
        let outcome = serde_json::from_str::<RootEnclosureReport>(&input)
            .map_err(error_text)
            .and_then(|report| {
                verify_root_enclosures(&report).map_err(error_text)?;
                Ok(serde_json::json!({
                    "schema": ROOT_ENCLOSURES_SCHEMA,
                    "verified": true,
                    "polynomials": report.items.len(),
                }))
            });
        print_json_result(outcome, compact);
        return;
    }
    let outcome = serde_json::from_str::<RootEnclosureRequest>(&input)
        .map_err(error_text)
        .and_then(|request| root_enclosures(&request).map_err(error_text));
    print_json_result(outcome, compact);
}
