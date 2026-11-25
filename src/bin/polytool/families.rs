//! Thin CLI adapter over the family catalog library.
use crate::input::{next_option_value, usage_error};
use polytool::families::*;
use serde_json::json;
use std::io::Write;

pub(crate) fn cmd_families(args: &[String]) {
    let fail = |e: &str| -> ! { usage_error("families", e) };
    let Some((action, rest)) = args.split_first() else {
        fail("expected list, show, generate, export-oeis or validate")
    };
    let mut id = None;
    let mut count = if action == "export-oeis" { 12 } else { 8 };
    let mut format = if action == "show" { "json" } else { "rows" };
    let mut query = String::new();
    let mut status = None;
    let mut bfile = None;
    let mut lean_export = None;
    let mut i = 0;
    while i < rest.len() {
        let flag = rest[i].as_str();
        match flag {
            "--json" => format = "json",
            "--format" | "--count" | "--query" | "--status" | "--bfile" | "--lean-export" => {
                let value = next_option_value(rest, &mut i, flag).unwrap_or_else(|e| fail(&e));
                match flag {
                    "--format" => format = value,
                    "--count" => {
                        count = value
                            .parse::<usize>()
                            .ok()
                            .filter(|n| (1..=512).contains(n))
                            .unwrap_or_else(|| fail("count must be 1..512"))
                    }
                    "--query" => query = value.to_lowercase(),
                    "--status" => status = Some(value),
                    "--bfile" => bfile = Some(value),
                    "--lean-export" => lean_export = Some(value),
                    _ => unreachable!(),
                }
            }
            _ if !flag.starts_with('-') && id.is_none() => id = Some(flag.to_uppercase()),
            _ => fail("unexpected argument"),
        }
        i += 1;
    }
    if let Some(status) = status {
        serde_json::from_value::<FamilyProofStatus>(json!(status))
            .unwrap_or_else(|_| fail("unknown proof status"));
    }
    let records = catalog().unwrap_or_else(|e| fail(&e.to_string()));
    let output = match action.as_str() {
        "validate" => {
            validate_catalog(&records).unwrap_or_else(|e| fail(&e.to_string()));
            if let Some(path) = lean_export {
                let text = std::fs::read_to_string(path).unwrap_or_else(|e| fail(&e.to_string()));
                let declarations = serde_json::from_str::<Vec<FamilyLeanDeclaration>>(&text)
                    .unwrap_or_else(|e| fail(&e.to_string()));
                validate_lean_declarations(&records, &declarations)
                    .unwrap_or_else(|e| fail(&e.to_string()));
            }
            json!({"schema":"polytool.families-validation.v1","valid":true,"count":records.len()})
        }
        "list" => {
            let entries:Vec<_>=records.iter().filter(|r|query.is_empty() || r.id.to_lowercase().contains(&query) || r.title.to_lowercase().contains(&query))
                .filter(|r|status.is_none_or(|s|serde_json::to_value(r.real_rootedness.status).unwrap()==s))
                .map(|r|json!({"id":r.id,"title":r.title,"status":r.real_rootedness.status,"interlacing":r.interlacing.status,"degree_pattern":r.degree_pattern,"lean":r.real_rootedness.lean,"models":r.models})).collect();
            if format != "json" {
                for r in &entries {
                    println!(
                        "{}\t{}\t{}",
                        r["id"].as_str().unwrap(),
                        r["status"].as_str().unwrap(),
                        r["title"].as_str().unwrap()
                    );
                }
                return;
            }
            json!({"schema":"polytool.families.v1","count":entries.len(),"families":entries})
        }
        "show" | "generate" | "export-oeis" => {
            let id = id.unwrap_or_else(|| fail("missing family ID"));
            let r = records
                .iter()
                .find(|r| r.id == id)
                .unwrap_or_else(|| fail("unknown family ID"));
            match action.as_str() {
                "show" => json!({"schema":"polytool.family.v1","family":r}),
                "generate" => {
                    let rows = r
                        .generate_rows(count)
                        .unwrap_or_else(|e| fail(&e.to_string()));
                    if format != "json" {
                        if format == "terms" {
                            println!(
                                "{}",
                                rows.into_iter().flatten().collect::<Vec<_>>().join(", ")
                            );
                        } else {
                            for row in rows {
                                println!("{}", row.join(", "));
                            }
                        }
                        return;
                    }
                    json!({"schema":"polytool.family-rows.v1","id":r.id,"first_row":r.first_row,"rows":rows})
                }
                _ => {
                    let export = export_oeis(r, count).unwrap_or_else(|e| fail(&e.to_string()));
                    if let Some(path) = bfile {
                        let mut file = std::fs::OpenOptions::new()
                            .write(true)
                            .create_new(true)
                            .open(path)
                            .unwrap_or_else(|e| fail(&e.to_string()));
                        file.write_all(export.bfile.as_bytes())
                            .unwrap_or_else(|e| fail(&e.to_string()));
                    }
                    if format != "json" {
                        println!(
                            "{}\n# oeis-contributions TOML skeleton\n{}",
                            export.blocks, export.proposal_toml
                        );
                        return;
                    }
                    serde_json::to_value(export).unwrap()
                }
            }
        }
        _ => fail("unknown families command"),
    };
    println!("{}", serde_json::to_string(&output).unwrap());
}
