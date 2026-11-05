//! Discover and replay general recurrence proofs through the library registry.

use crate::input::{error_text, is_help_arg, parse_usize_option, read_stdin_or_exit};
use crate::recurrence_input::read_recurrence_json;
use polytool::recurrence::{
    discover_recurrence_proofs, RealRootednessTheoremOptions, RecurrenceProofDiscoveryReport,
};

pub(crate) fn print_recurrence_prove_help() {
    println!("Usage: polytool recurrence-prove --recurrence <file|-> [--candidate N] [--json]");
    println!("       polytool recurrence-prove --verify [--json] < proof.json");
    println!("       polytool recurrence-prove --help");
    println!();
    println!("Discover general proofs for a recurrence and seeds; no named-sequence lookup.");
    println!("Tries product induction, Favard/three-term, Liu-Wang, derivative-lag,");
    println!("Ma-Wang interval and Liu-Wang three-term criteria, with exact all-index");
    println!("hypothesis checks. Reports every attempt, its obligations and conclusion.");
    println!("A fit still needs an all-index identity for an independently defined family.");
    println!("Rust certificates can be replayed; Lean replay is not attempted.");
    println!();
    println!("  --recurrence <file|->  recurrence.v1, theorem wrapper or representation menu");
    println!(
        "  --candidate N          Explicit 0-based menu choice (required for multiple entries)"
    );
    println!("  --max-start-delay N    Induction-start retries, default 3 (maximum 16)");
    println!("  --json                 Emit polytool.recurrence-proof-discovery.v1");
    println!("  --verify               Replay the complete report from stdin; no search flags");
    println!("Exit status: 0 proved (or valid replay), 1 no proof, 2 invalid input/replay.");
    println!();
    println!("Example: polytool recurrence --json < rows.txt | polytool recurrence-prove --recurrence - --json");
}

pub(crate) fn cmd_recurrence_prove(args: &[String]) {
    if args.iter().any(|s| is_help_arg(s)) {
        print_recurrence_prove_help();
        return;
    }
    match run(args) {
        Ok(proved) => {
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

fn run(args: &[String]) -> Result<bool, String> {
    let (mut path, mut candidate, mut delay) = (None, None, None);
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
            "--max-start-delay" => {
                delay = Some(parse_usize_option(args, &mut i, "--max-start-delay")?)
            }
            "--json" => json = true,
            "--verify" => verify = true,
            other => return Err(format!("unknown recurrence-prove option: {other}")),
        }
        i += 1;
    }
    if verify {
        if path.is_some() || candidate.is_some() || delay.is_some() {
            return Err("--verify cannot be combined with search flags".into());
        }
        let report: RecurrenceProofDiscoveryReport =
            serde_json::from_str(&read_stdin_or_exit()).map_err(error_text)?;
        if report.maximum_start_delay > 16 {
            return Err("maximum_start_delay exceeds 16".into());
        }
        report.verify().map_err(error_text)?;
        if json {
            println!(
                "{}",
                serde_json::json!({"verified": true, "status": report.status, "lean_replay": report.lean_replay})
            );
        } else {
            println!(
                "Report verified; discovery status: {}. Lean replay: {}.",
                report.status, report.lean_replay
            );
        }
        return Ok(true);
    }
    let maximum_start_delay = delay.unwrap_or(3);
    if maximum_start_delay > 16 {
        return Err("--max-start-delay must be at most 16".into());
    }
    let record = read_recurrence_json(
        path.ok_or("recurrence-prove needs --recurrence <file|->")?,
        candidate,
    )?;
    let report = discover_recurrence_proofs(
        &record,
        &RealRootednessTheoremOptions {
            maximum_start_delay,
            ..Default::default()
        },
    )
    .map_err(error_text)?;
    if json {
        println!(
            "{}",
            serde_json::to_string_pretty(&report).map_err(error_text)?
        );
    } else {
        println!("Proof discovery: {}", report.status);
        println!("{}", report.statement);
        for attempt in &report.attempts {
            println!("\n{:?}: {:?}", attempt.theorem, attempt.status);
            if let Some(conclusion) = &attempt.conclusion {
                println!("  {}", conclusion.statement);
            } else {
                println!("  {}", attempt.normal_form);
            }
            for hypothesis in &attempt.hypotheses {
                println!(
                    "  {} [{:?}]: {}",
                    hypothesis.name, hypothesis.status, hypothesis.evidence
                );
            }
            println!("  Reference: {}", attempt.reference);
            if !attempt.lean_declarations.is_empty() {
                println!(
                    "  Related Lean: {} (replay not attempted)",
                    attempt.lean_declarations.join(", ")
                );
            }
        }
    }
    Ok(report.status == "proved")
}
