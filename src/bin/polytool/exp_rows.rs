//! Generate exponential-family rows and run exact checks in process.
use crate::input::{next_option_value, usage_error};
use polytool::exponential_rows::*;

fn spec(value: &str) -> Result<BlockWeightSpec, String> {
    if value.starts_with('{') {
        return serde_json::from_str(value).map_err(|e| e.to_string());
    }
    if value == "none" {
        return Ok(BlockWeightSpec::Sizes { sizes: vec![] });
    }
    if let Some(s) = value.strip_prefix("sizes:") {
        return Ok(BlockWeightSpec::Sizes {
            sizes: if s.is_empty() {
                vec![]
            } else {
                s.split(',')
                    .map(|r| r.parse().map_err(|_| "invalid size".to_string()))
                    .collect::<Result<_, _>>()?
            },
        });
    }
    if let Some(s) = value.strip_prefix("tail:") {
        return Ok(BlockWeightSpec::Tail {
            minimum: s.parse().map_err(|_| "invalid tail minimum")?,
        });
    }
    if let Some(s) = value.strip_prefix("weights:") {
        return Ok(BlockWeightSpec::Finite {
            weights: s.split(',').map(str::to_string).collect(),
        });
    }
    Err("expected none, sizes:1,2, weights:1,1/2, tail:2, or a JSON block specification".into())
}
fn rows(value: &str) -> Result<Vec<usize>, String> {
    let parse = |s: &str| {
        s.parse::<usize>()
            .map_err(|_| "invalid row index".to_string())
    };
    if let Some((lo, hi)) = value.split_once("..=") {
        let (lo, hi) = (parse(lo)?, parse(hi)?);
        if lo > hi || hi == usize::MAX || hi - lo > 1_000_000 {
            return Err("invalid or excessive row range".into());
        }
        return Ok((lo..=hi).collect());
    }
    value.split(',').map(parse).collect()
}
pub(crate) fn cmd_exp_rows(args: &[String]) {
    let fail = |e: &str| -> ! { usage_error("exp-rows", e) };
    let mut request = None;
    let mut unmarked = None;
    let mut marked = None;
    let mut indices = None;
    let mut check = ExponentialRowCheck::None;
    let mut summary = false;
    let mut json = false;
    let mut verify = false;
    let mut i = 0;
    while i < args.len() {
        let flag = args[i].as_str();
        match flag {
            "--json" => json = true,
            "--summary" => summary = true,
            "--verify" => verify = true,
            "--format" | "--input" | "--unmarked" | "--marked" | "--rows" | "--check" => {
                let value = next_option_value(args, &mut i, flag).unwrap_or_else(|e| fail(&e));
                match flag {
                    "--format" => json = value == "json",
                    "--input" => {
                        let text = if value == "-" {
                            std::io::read_to_string(std::io::stdin())
                                .unwrap_or_else(|e| fail(&e.to_string()))
                        } else {
                            std::fs::read_to_string(value).unwrap_or_else(|e| fail(&e.to_string()))
                        };
                        request = Some(text);
                    }
                    "--unmarked" => unmarked = Some(spec(value).unwrap_or_else(|e| fail(&e))),
                    "--marked" => marked = Some(spec(value).unwrap_or_else(|e| fail(&e))),
                    "--rows" => indices = Some(rows(value).unwrap_or_else(|e| fail(&e))),
                    "--check" => {
                        check = match value {
                            "none" => ExponentialRowCheck::None,
                            "real-rooted" => ExponentialRowCheck::RealRooted,
                            "interlacing" => ExponentialRowCheck::Interlacing,
                            "first-failure" => ExponentialRowCheck::FirstFailure,
                            _ => fail("unknown check"),
                        }
                    }
                    _ => unreachable!(),
                }
            }
            _ => fail("unknown option"),
        }
        i += 1;
    }
    if verify {
        let saved: ExponentialRowsReport = serde_json::from_str(
            request
                .as_deref()
                .unwrap_or_else(|| fail("--verify requires --input")),
        )
        .unwrap_or_else(|e| fail(&e.to_string()));
        let replay =
            check_exponential_rows(&saved.request).unwrap_or_else(|e| fail(&e.to_string()));
        if replay != saved {
            fail("report does not replay");
        }
        println!(
            "{}",
            if json {
                "{\"verified\":true}"
            } else {
                "Exact report replay verified"
            }
        );
        return;
    }
    let request = if let Some(text) = request {
        if unmarked.is_some()
            || marked.is_some()
            || indices.is_some()
            || check != ExponentialRowCheck::None
            || summary
        {
            fail("--input is exclusive with family, rows, check and summary options");
        }
        serde_json::from_str(&text).unwrap_or_else(|e| fail(&e.to_string()))
    } else {
        ExponentialRowsRequest {
            family: ExponentialFamily {
                unmarked: unmarked.unwrap_or(BlockWeightSpec::Sizes { sizes: vec![] }),
                marked: marked.unwrap_or_else(|| fail("--marked is required")),
            },
            rows: indices.unwrap_or_else(|| fail("--rows is required")),
            check,
            summary,
            finite_rank: vec![],
        }
    };
    let report = check_exponential_rows(&request).unwrap_or_else(|e| fail(&e.to_string()));
    if json {
        println!(
            "{}",
            serde_json::to_string(&report).expect("serializable report")
        );
    } else {
        for row in &report.rows {
            println!(
                "n={}: {} degree={} real_rooted={}",
                row.index,
                row.coefficients
                    .as_ref()
                    .map(|r| r.join(", "))
                    .unwrap_or_default(),
                degree_label(row.degree),
                verdict_label(row.real_rooted, "not checked")
            );
        }
        for pair in &report.interlacing {
            println!(
                "{} -> {}: forward={}, reverse={}",
                pair.left_index,
                pair.right_index,
                verdict_label(pair.forward, "not applicable"),
                verdict_label(pair.reverse, "not applicable")
            );
        }
        for route in &report.routes {
            println!(
                "route {} ({}): {}",
                route.route,
                evidence_label(route.evidence),
                route.statement
            );
        }
        if let Some(f) = &report.first_failure {
            println!(
                "First failure n={}: {}, degree={}, distinct_real_roots={}, squarefree_degree={}",
                f.index,
                f.property,
                degree_label(f.degree),
                f.distinct_real_roots,
                f.squarefree_degree
            );
        }
        for cert in &report.finite_rank {
            println!(
                "rank {}: {} certified={}",
                cert.request.index, cert.route, cert.certified
            );
        }
    }
    if report.first_failure.is_some() {
        std::process::exit(1);
    }
}

fn degree_label(degree: Option<usize>) -> String {
    degree.map_or_else(|| "zero polynomial".into(), |d| d.to_string())
}
fn verdict_label(verdict: Option<bool>, absent: &'static str) -> &'static str {
    match verdict {
        Some(true) => "yes",
        Some(false) => "no",
        None => absent,
    }
}
fn evidence_label(level: ExponentialEvidenceLevel) -> &'static str {
    match level {
        ExponentialEvidenceLevel::TheoremBacked => "theorem backed",
        ExponentialEvidenceLevel::PaperProof => "paper proof",
        ExponentialEvidenceLevel::UnderReview => "under review",
        ExponentialEvidenceLevel::FiniteRank => "finite rank",
    }
}
