//! Read triangle rows from a local OEIS checkout.

use polytool::oeis::{infer_oeis_rows, OeisBFileRejection, OeisRowLengthRule};
use serde_json::json;
use std::path::Path;

#[derive(Debug)]
enum LocalRowReadError {
    MissingSequence {
        relative_path: String,
    },
    ReadSequence {
        relative_path: String,
        source: std::io::Error,
    },
}
impl std::fmt::Display for LocalRowReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MissingSequence { relative_path } => write!(f, "missing OEIS sequence file: expected --data/{relative_path} (pattern --data/seq/Annn/Annnnnn.seq)"),
            Self::ReadSequence { relative_path, source } => write!(f, "cannot read --data/{relative_path}: {source}"),
        }
    }
}
impl std::error::Error for LocalRowReadError {}

pub(super) fn run(args: &[String]) {
    let mut id = None;
    let mut data = None;
    let mut json_output = false;
    let mut iter = args.iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--json" => json_output = true,
            "--data" => data = iter.next(),
            value if !value.starts_with('-') && id.is_none() => id = Some(value),
            _ => crate::input::usage_error("oeis rows", "invalid argument"),
        }
    }
    let Some(id) = id.filter(|id| {
        id.len() == 7 && id.starts_with('A') && id[1..].bytes().all(|b| b.is_ascii_digit())
    }) else {
        crate::input::usage_error("oeis rows", "expected an A-number with six digits");
    };
    let Some(data) = data else {
        crate::input::usage_error("oeis rows", "--data DIR is required");
    };
    let root = Path::new(data);
    let relative = format!("seq/{}/{id}.seq", &id[..4]);
    let text = std::fs::read_to_string(root.join(&relative))
        .map_err(|source| {
            if source.kind() == std::io::ErrorKind::NotFound {
                LocalRowReadError::MissingSequence {
                    relative_path: relative.clone(),
                }
            } else {
                LocalRowReadError::ReadSequence {
                    relative_path: relative.clone(),
                    source,
                }
            }
        })
        .unwrap_or_else(|error| {
            eprintln!("{error}");
            std::process::exit(2)
        });
    let bpath = root
        .join("files")
        .join(&id[..4])
        .join(format!("b{}.txt", &id[1..]));
    let (bfile, unreadable) = match std::fs::read_to_string(&bpath) {
        Ok(text) => (Some(text), None),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => (None, None),
        Err(error) => (None, Some(error.to_string())),
    };
    let mut report = infer_oeis_rows(&text, id, bfile.as_deref()).unwrap_or_else(|error| {
        eprintln!("{error}");
        std::process::exit(2)
    });
    if unreadable.is_some() {
        report.bfile_rejection = Some(OeisBFileRejection::Unreadable);
    }
    let rule = report.row_length_rule.as_ref().map(|rule| match rule {
        OeisRowLengthRule::ConstantIncrease(step) => format!("+{step} per row"),
        OeisRowLengthRule::PeriodicIncrease { period, phase } => {
            format!("+1 every {period} rows (transition phase {phase})")
        }
    });
    let rule_details = report.row_length_rule.as_ref().map(|rule| match rule {
        OeisRowLengthRule::ConstantIncrease(step) => {
            json!({"kind":"constant_step", "step":step, "period":1, "phase":0})
        }
        OeisRowLengthRule::PeriodicIncrease { period, phase } => {
            json!({"kind":"periodic_step", "step":1, "period":period, "phase":phase})
        }
    });
    if json_output {
        let rows: Vec<Vec<String>> = report
            .rows
            .iter()
            .map(|row| row.iter().map(ToString::to_string).collect())
            .collect();
        println!("{}", serde_json::to_string_pretty(&json!({"id": id, "name": report.name, "keywords": report.keywords, "term_source": report.term_source.as_str(), "row_length_rule": rule, "row_length_rule_details": rule_details, "rows": rows, "status": if report.reason.is_some() { "refused" } else { "ok" }, "reason": report.reason.map(|r| r.message()), "reason_code": report.reason.map(|r| r.as_str()), "bfile_rejection": report.bfile_rejection.map(|r| r.as_str()), "bfile_rejection_reason": report.bfile_rejection.map(|r| r.message()), "bfile_rejection_detail": unreadable})).unwrap());
    } else if let Some(reason) = report.reason {
        if let Some(rejection) = report.bfile_rejection {
            eprintln!(
                "{id}: {}{}",
                rejection.message(),
                unreadable
                    .as_ref()
                    .map(|s| format!(": {s}"))
                    .unwrap_or_default()
            );
        }
        eprintln!("{id}: {reason} (term source: {})", report.term_source);
    } else {
        if let Some(rejection) = report.bfile_rejection {
            eprintln!(
                "{id}: {}{}",
                rejection.message(),
                unreadable
                    .as_ref()
                    .map(|s| format!(": {s}"))
                    .unwrap_or_default()
            );
        }
        eprintln!("{id}: term source: {}", report.term_source);
        for row in report.rows {
            println!(
                "{}",
                row.iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join(" ")
            );
        }
    }
}
