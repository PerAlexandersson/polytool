//! Exact BKW and all-row theorem certificates, with deterministic replay.

use crate::input::{error_text, parse_usize_option, read_stdin_or_exit};
use crate::recurrence_input::read_recurrence_json;
use polytool::bkw::certified::{certify_recurrence, CertificationOptions, CertificationReport};

pub(crate) fn cmd_bkw(args: &[String]) {
    if args.first().map(String::as_str) != Some("certify") {
        crate::input::usage_error("bkw", "expected certify");
    }
    match run(&args[1..]) {
        Ok(()) => {}
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(2);
        }
    }
}

fn run(args: &[String]) -> Result<(), String> {
    let mut options = CertificationOptions::default();
    let (mut path, mut candidate, mut json, mut verify) = (None, None, false, false);
    let mut supplied_options = false;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--recurrence" => {
                i += 1;
                path = Some(
                    args.get(i)
                        .ok_or("--recurrence requires FILE or -")?
                        .as_str(),
                );
            }
            "--candidate" => candidate = Some(parse_usize_option(args, &mut i, "--candidate")?),
            "--max-row-index" => {
                options.maximum_row_index = parse_usize_option(args, &mut i, "--max-row-index")?;
                supplied_options = true;
            }
            "--max-pencil-degree" => {
                options.maximum_pencil_degree =
                    parse_usize_option(args, &mut i, "--max-pencil-degree")?;
                supplied_options = true;
            }
            "--radius" => {
                i += 1;
                options.circle_radius =
                    Some(args.get(i).ok_or("--radius requires a rational")?.clone());
                supplied_options = true;
            }
            "--json" => json = true,
            "--verify" => verify = true,
            x => return Err(format!("unknown bkw certify option: {x}")),
        }
        i += 1;
    }
    if verify {
        if path.is_some() || candidate.is_some() || supplied_options {
            return Err("--verify takes only a certificate on stdin".into());
        }
        let report: CertificationReport =
            serde_json::from_str(&read_stdin_or_exit()).map_err(error_text)?;
        report.verify().map_err(error_text)?;
        if json {
            println!(
                "{}",
                serde_json::json!({"verified":true,"result":report.result})
            );
        } else {
            println!("BKW certificate verified: {}", report.result);
        }
        return Ok(());
    }
    let record = read_recurrence_json(path.unwrap_or("-"), candidate)?;
    let report = certify_recurrence(&record, &options).map_err(error_text)?;
    if json {
        println!(
            "{}",
            serde_json::to_string_pretty(&report).map_err(error_text)?
        );
    } else {
        println!(
            "{}; limit set: {:?}\n{}",
            report.result, report.limit_set_status, report.reason
        );
        if let Some(row) = &report.first_failing_row {
            println!("First failing row {}: {:?}", row.index, row.coefficients);
        }
        if let Some(circle) = &report.circle {
            println!("Circle radius: {}; {}", circle.radius, circle.conclusion);
        }
    }
    Ok(())
}
