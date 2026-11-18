//! Thin exact stability CLI adapter.

use crate::input::{next_option_value, read_command_input, usage_error};
use polytool::stability::{
    check_univariate_stability_with_convention,
    parse::{parse_gaussian, parse_univariate, ParseLimits},
    wire,
};
use serde_json::json;

pub(crate) fn cmd_stability(args: &[String]) {
    if let Err(error) = run(args) {
        usage_error("stability", &error);
    }
}

fn run(args: &[String]) -> Result<(), String> {
    let mut limits = ParseLimits {
        max_degree: 2000,
        max_coefficient_bits: 1_000_000,
        max_literal_bytes: 16 * 1024 * 1024,
    };
    if args.first().is_some_and(|s| s == "verify") {
        let text = read_command_input("stability report")?;
        let value = serde_json::from_str(&text).map_err(|e| e.to_string())?;
        let verified = wire::verify_report(&value, limits).map_err(|e| e.to_string())?;
        println!(
            "{}",
            json!({"schema":"polytool.stability-verification.v1","verified":verified})
        );
        if !verified {
            std::process::exit(1);
        }
        return Ok(());
    }
    let (mut region, mut mobius, mut center, mut radius, mut convention, mut certificate) =
        (None, None, None, None, None, None);
    let mut format = "summary";
    let mut expressions = Vec::new();
    let mut i = 0;
    while i < args.len() {
        let flag = args[i].as_str();
        match flag {
            "--json" => format = "json",
            "--format" => format = next_option_value(args, &mut i, flag)?,
            "--region" => region = Some(next_option_value(args, &mut i, flag)?),
            "--mobius" => mobius = Some(next_option_value(args, &mut i, flag)?),
            "--disk-center" => center = Some(next_option_value(args, &mut i, flag)?),
            "--radius" => radius = Some(next_option_value(args, &mut i, flag)?),
            "--convention" => convention = Some(next_option_value(args, &mut i, flag)?),
            "--certificate" => certificate = Some(next_option_value(args, &mut i, flag)?),
            "--max-degree" => {
                limits.max_degree = next_option_value(args, &mut i, flag)?
                    .parse()
                    .map_err(|_| "invalid degree cap")?
            }
            "--max-coefficient-bits" => {
                limits.max_coefficient_bits = next_option_value(args, &mut i, flag)?
                    .parse()
                    .map_err(|_| "invalid coefficient bit cap")?
            }
            s if s.starts_with("--") => return Err(format!("unknown option {s}")),
            s => expressions.push(s.to_string()),
        }
        i += 1;
    }
    if limits.max_degree > 2000 || limits.max_coefficient_bits > 1_000_000 {
        return Err("caps cannot exceed degree 2000 / coefficient-part 1000000 bits".into());
    }
    let mobius = mobius
        .map(|s| {
            s.split(',')
                .map(|c| parse_gaussian(c, limits))
                .collect::<Result<Vec<_>, _>>()
        })
        .transpose()
        .map_err(|e| e.to_string())?;
    let region = wire::parse_region(region, mobius.as_deref(), center, radius, limits)
        .map_err(|e| e.to_string())?;
    let convention = wire::parse_convention(convention, &region).map_err(|e| e.to_string())?;
    if expressions.is_empty() {
        expressions = read_command_input("stability polynomials")?
            .lines()
            .map(str::trim)
            .filter(|s| !s.is_empty() && !s.starts_with('#'))
            .map(str::to_string)
            .collect();
    }
    if expressions.is_empty() {
        return Err("expected at least one polynomial".into());
    }
    if expressions.iter().map(String::len).sum::<usize>() > 16 * 1024 * 1024 {
        return Err("polynomial text exceeds 16 MiB".into());
    }
    let mut summaries = Vec::new();
    let mut full = Vec::new();
    for (index, expression) in expressions.iter().enumerate() {
        let coefficients =
            parse_univariate(expression, limits).map_err(|e| format!("item {index}: {e}"))?;
        let report = check_univariate_stability_with_convention(&coefficients, &region, convention)
            .map_err(|e| e.to_string())?;
        summaries.push(wire::summary(index, &report));
        if format == "full" || certificate.is_some() {
            full.push(wire::full_item(index, &coefficients, &report));
        }
    }
    let full = wire::full_report(full);
    if let Some(path) = certificate {
        std::fs::write(path, format!("{full}\n"))
            .map_err(|e| format!("cannot write certificate: {e}"))?;
    }
    match format {
        "json" => println!("{}", json!({"items":summaries})),
        "full" => println!("{full}"),
        "summary" => {
            for item in summaries {
                println!(
                    "{}: degree={} counts={} verdict={} meaning={} method={}",
                    item["index"],
                    item["degree"],
                    item["counts"],
                    item["verdict"],
                    item["meaning"].as_str().unwrap(),
                    item["method"].as_str().unwrap()
                );
            }
        }
        _ => return Err("invalid output format".into()),
    }
    Ok(())
}
