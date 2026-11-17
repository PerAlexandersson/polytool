//! JSON adapter for stationary interlacing-matrix orbit certificates.

use crate::input::{error_text, read_stdin_or_exit};
use polytool::recurrence::{
    check_matrix_orbit, verify_matrix_orbit, MatrixOrbitReport, MatrixOrbitRequest,
};

pub(crate) fn cmd_matrix_orbit(args: &[String]) {
    let verify = args.iter().any(|a| a == "--verify");
    let compact = args.iter().any(|a| a == "--compact");
    if args
        .iter()
        .any(|a| a != "--verify" && a != "--json" && a != "--compact")
    {
        eprintln!("unknown matrix-orbit option");
        std::process::exit(2);
    }
    let input = read_stdin_or_exit();
    let outcome: Result<(serde_json::Value, bool), String> = if verify {
        serde_json::from_str::<MatrixOrbitReport>(&input).map_err(error_text).and_then(|report| {
            verify_matrix_orbit(&report).map_err(error_text)?;
            Ok((serde_json::json!({"verified":true,"status":report.status,"lean_replay":report.lean_replay}),true))
        })
    } else {
        serde_json::from_str::<MatrixOrbitRequest>(&input)
            .map_err(error_text)
            .and_then(|request| {
                let report = check_matrix_orbit(&request).map_err(error_text)?;
                let proved = report.status == "proved";
                Ok((serde_json::to_value(report).map_err(error_text)?, proved))
            })
    };
    match outcome {
        Ok((value, proved)) => {
            let text = if compact {
                serde_json::to_string(&value)
            } else {
                serde_json::to_string_pretty(&value)
            };
            println!("{}", text.expect("serialize matrix-orbit output"));
            if !proved {
                std::process::exit(1);
            }
        }
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(2);
        }
    }
}
