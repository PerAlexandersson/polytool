//! `polytool recurrence-degree`: the degree law and the signs of the leading
//! coefficients of a recurrence family for every index.

use crate::input::{error_text, parse_usize_option, read_stdin_or_exit};
use crate::recurrence_input::read_recurrence_json;
use polytool::recurrence::{
    certify_leading_coefficient, LeadingCoefficientCertificate, LeadingCoefficientOptions,
    LiftHypothesisStatus,
};

pub(crate) fn cmd_recurrence_degree(args: &[String]) {
    match run(args) {
        Ok(true) => {}
        Ok(false) => std::process::exit(1),
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(2);
        }
    }
}

fn run(args: &[String]) -> Result<bool, String> {
    let (mut path, mut candidate) = (None, None);
    let (mut rows, mut period) = (None, None);
    let (mut json, mut verify) = (false, false);
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--recurrence" => {
                i += 1;
                path = Some(
                    args.get(i)
                        .ok_or("--recurrence needs a file or -")?
                        .as_str(),
                );
            }
            "--candidate" => candidate = Some(parse_usize_option(args, &mut i, "--candidate")?),
            "--count" => rows = Some(parse_usize_option(args, &mut i, "--count")?),
            "--max-period" => period = Some(parse_usize_option(args, &mut i, "--max-period")?),
            "--json" => json = true,
            "--verify" => verify = true,
            other => return Err(format!("unknown recurrence-degree option: {other}")),
        }
        i += 1;
    }
    if verify {
        if path.is_some() || candidate.is_some() || rows.is_some() || period.is_some() {
            return Err("--verify cannot be combined with search flags".into());
        }
        let certificate: LeadingCoefficientCertificate =
            serde_json::from_str(&read_stdin_or_exit()).map_err(error_text)?;
        check_limits(&certificate.options)?;
        certificate.verify().map_err(error_text)?;
        if json {
            println!(
                "{}",
                serde_json::json!({"verified": true, "status": certificate.status})
            );
        } else {
            println!("Certificate verified; status: {:?}.", certificate.status);
        }
        return Ok(true);
    }
    let options = LeadingCoefficientOptions {
        sample_rows: rows.unwrap_or(LeadingCoefficientOptions::default().sample_rows),
        maximum_period: period.unwrap_or(LeadingCoefficientOptions::default().maximum_period),
    };
    check_limits(&options)?;
    let record = read_recurrence_json(
        path.ok_or("recurrence-degree needs --recurrence <file|->")?,
        candidate,
    )?;
    let (rule, first_index, seeds) = record.to_recurrence_parts().map_err(error_text)?;
    let certificate =
        certify_leading_coefficient(&rule, first_index, &seeds, &options).map_err(error_text)?;
    if json {
        println!(
            "{}",
            serde_json::to_string_pretty(&certificate).map_err(error_text)?
        );
    } else {
        print_text(&certificate);
    }
    Ok(certificate.is_proved())
}

/// Bounds on the search options (also for a certificate read with
/// `--verify`, whose replay may rerun the search).
fn check_limits(options: &LeadingCoefficientOptions) -> Result<(), String> {
    if options.sample_rows > 400 {
        return Err("--count must be at most 400".into());
    }
    if options.maximum_period == 0 || options.maximum_period > 24 {
        return Err("--max-period must be between 1 and 24".into());
    }
    Ok(())
}

fn print_text(certificate: &LeadingCoefficientCertificate) {
    println!("Degree law: {:?}", certificate.status);
    println!("Rule: {}", certificate.normal_form);
    if let Some(conclusion) = &certificate.conclusion {
        println!("{conclusion}");
    }
    for hypothesis in &certificate.hypotheses {
        let status = match &hypothesis.status {
            LiftHypothesisStatus::Verified => "verified".to_string(),
            LiftHypothesisStatus::Failed { witness: Some(n) } => format!("failed at n = {n}"),
            LiftHypothesisStatus::Failed { witness: None } => "failed".to_string(),
            LiftHypothesisStatus::Unresolved => "unresolved".to_string(),
        };
        println!("  {} [{status}]: {}", hypothesis.name, hypothesis.evidence);
    }
    for attempt in &certificate.attempts {
        println!("  attempt: {attempt}");
    }
    if !certificate.lean_declarations.is_empty() {
        println!(
            "Related Lean: route {} ({}); replay not attempted",
            certificate.lean_route.as_deref().unwrap_or("-"),
            certificate.lean_declarations.join(", ")
        );
    }
}
