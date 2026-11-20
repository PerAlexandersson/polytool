//! Thin CLI adapter for effective dominant-pair certificates and replay.

use crate::input::{error_text, parse_usize_option, read_stdin_or_exit};
use crate::recurrence_input::read_recurrence_json;
use polytool::bkw::effective::{
    certify_effective_bkw, certify_effective_bkw_recurrence, EffectiveBkwReport,
    EffectiveBkwRequest,
};

pub(crate) fn cmd_bkw_effective(args: &[String]) {
    if let Err(error) = run(args) {
        eprintln!("{error}");
        std::process::exit(2);
    }
}
fn run(args: &[String]) -> Result<(), String> {
    let (mut request_path, mut recurrence_path, mut candidate, mut json, mut verify) =
        (None, None, None, false, false);
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--request" | "--recurrence" => {
                let flag = args[i].as_str();
                i += 1;
                let path = args
                    .get(i)
                    .ok_or_else(|| format!("{flag} requires FILE or -"))?
                    .as_str();
                if flag == "--request" {
                    request_path = Some(path);
                } else {
                    recurrence_path = Some(path);
                }
            }
            "--candidate" => candidate = Some(parse_usize_option(args, &mut i, "--candidate")?),
            "--verify" => verify = true,
            "--json" => json = true,
            s => return Err(format!("unknown bkw effective option: {s}")),
        }
        i += 1;
    }
    if verify {
        if request_path.is_some() || recurrence_path.is_some() || candidate.is_some() {
            return Err("--verify accepts only a report on stdin".into());
        }
        let report: EffectiveBkwReport =
            serde_json::from_str(&read_stdin_or_exit()).map_err(error_text)?;
        report.verify().map_err(error_text)?;
        if json {
            println!(
                "{}",
                serde_json::json!({"verified":true,"status":report.status})
            );
        } else {
            println!("Effective-BKW certificate verified: {}", report.status);
        }
        return Ok(());
    }
    if request_path.is_some() && recurrence_path.is_some() {
        return Err("choose --request or --recurrence".into());
    }
    let report = if let Some(path) = recurrence_path {
        let record = read_recurrence_json(path, candidate)?;
        certify_effective_bkw_recurrence(&record).map_err(error_text)?
    } else {
        if candidate.is_some() {
            return Err("--candidate requires --recurrence".into());
        }
        let text = match request_path.unwrap_or("-") {
            "-" => read_stdin_or_exit(),
            path => std::fs::read_to_string(path).map_err(error_text)?,
        };
        let request: EffectiveBkwRequest = serde_json::from_str(&text).map_err(error_text)?;
        certify_effective_bkw(&request).map_err(error_text)?
    };
    if json {
        println!(
            "{}",
            serde_json::to_string_pretty(&report).map_err(error_text)?
        );
    } else {
        println!("{}: {}", report.status, report.reason);
        if let Some(c) = report.certificate {
            println!(
                "Local cutoff: {}; finite rows: {}; zero count: {}*n + {:?}",
                c.cutoff,
                c.finite_rows.len(),
                c.zero_count.slope,
                c.zero_count.offsets
            );
        }
    }
    Ok(())
}
