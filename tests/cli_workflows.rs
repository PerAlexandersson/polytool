//! CLI error contracts, nested help and recurrence pipelines.

use serde_json::{json, Value};
use std::io::Write;
use std::process::{Command, Output, Stdio};

fn run(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_polytool"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    if let Err(error) = child.stdin.take().unwrap().write_all(input.as_bytes()) {
        assert_eq!(error.kind(), std::io::ErrorKind::BrokenPipe);
    }
    child.wait_with_output().unwrap()
}

fn success(output: Output) -> String {
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

#[test]
fn invalid_options_and_missing_values_fail() {
    for command in [
        "real-rooted",
        "properties",
        "interlacing",
        "interlacing-profile",
        "gamma-expansion",
        "coefficient-tests",
        "family-check",
        "recurrence",
        "recurrence-generate",
        "resultant",
        "discriminant",
        "bkw-scout",
        "cyclic-sieving",
        "cyclic-sieving-sequence",
        "bicyclic-sieving",
        "hstar-inequalities",
        "hstar-to-ehrhart",
        "ehrhart-to-hstar",
        "pf-pencil",
        "sequence",
    ] {
        let output = run(&[command, "--typo"], "");
        assert!(
            !output.status.success(),
            "{command} ignored an invalid option"
        );
        assert!(!output.stderr.is_empty(), "{command} needs a diagnostic");
    }
    for args in [
        vec!["recurrence", "--max-candidates", "bogus"],
        vec!["recurrence-generate", "--rows", "bogus"],
        vec!["recurrence-generate", "--additional"],
        vec!["sequence", "eulerian", "5", "--format"],
    ] {
        let output = run(&args, "");
        assert_eq!(output.status.code(), Some(2), "{args:?}");
    }
}

#[test]
fn nested_help_is_detailed_and_does_not_read_input() {
    for (subcommand, marker) in [
        ("add", "--ids-file"),
        ("check", "--oeis-data"),
        ("evidence", "Usage:"),
    ] {
        let direct = success(run(&["oeis", subcommand, "--help"], ""));
        let prefixed = success(run(&["help", "oeis", subcommand], ""));
        assert_eq!(direct, prefixed);
        assert!(direct.contains(marker));
        assert!(direct.contains(&format!("polytool oeis {subcommand}")));
    }
    assert!(!run(&["help", "oeis", "unknown"], "").status.success());
}

#[test]
fn row_generation_search_wrapper_and_replay_compose() {
    let rows = success(run(&["sequence", "eulerian", "12", "--format", "rows"], ""));
    assert!(rows.starts_with("1\n1, 1\n1, 4, 1\n"));
    let report = success(run(
        &[
            "recurrence",
            "--json",
            "--theorem-matches",
            "--max-rec-len",
            "1",
            "--max-diff-deg",
            "1",
            "--max-candidates",
            "500",
        ],
        &rows,
    ));
    let value: Value = serde_json::from_str(&report).unwrap();
    assert_eq!(value["schema"], "polytool.recurrence-theorems.v1");
    let replay = success(run(
        &[
            "recurrence-generate",
            "--recurrence",
            "-",
            "--rows",
            "12",
            "--format",
            "rows",
        ],
        &report,
    ));
    assert_eq!(rows, replay);
}

#[test]
fn menu_selection_is_explicit_and_preserves_the_selected_rule() {
    let rule = |factor: &str| {
        json!({
            "schema":"polytool.recurrence.v1", "first_index":0,
            "initial_polynomials":[["1"]],
            "recurrence":{"terms":[{"offset":1,"deriv_order":0,"sign":"none",
                "coeff":{"coeffs":[[factor]]}}]}
        })
    };
    let menu = json!({"schema":"polytool.recurrence-representations.v1", "representations":[
        {"recurrence_json":rule("2")}, {"recurrence_json":rule("3")}
    ]})
    .to_string();
    let args = ["recurrence-generate", "--recurrence", "-", "--rows", "4"];
    let ambiguous = run(&args, &menu);
    assert_eq!(ambiguous.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&ambiguous.stderr).contains("--candidate"));
    for (index, expected) in [("0", "1\n2\n4\n8\n"), ("1", "1\n3\n9\n27\n")] {
        let selected = [&args[..], &["--candidate", index]].concat();
        assert_eq!(success(run(&selected, &menu)), expected);
    }
    let invalid = [&args[..], &["--candidate", "2"]].concat();
    assert_eq!(run(&invalid, &menu).status.code(), Some(2));
    assert_eq!(run(&args, "{}").status.code(), Some(2));
    let selected = [&args[..], &["--candidate", "0"]].concat();
    assert_eq!(
        run(&selected, &rule("2").to_string()).status.code(),
        Some(2)
    );
}

#[test]
fn exact_scalar_json_and_resultant_arity() {
    let discriminants: Value = serde_json::from_str(&success(run(
        &["discriminant", "--json"],
        "2,-3,1\n1,0,1\n",
    )))
    .unwrap();
    assert_eq!(discriminants["items"][0]["discriminant"], "1");
    assert_eq!(discriminants["items"][1]["discriminant"], "-4");
    let resultant: Value =
        serde_json::from_str(&success(run(&["resultant", "--json"], "-1,1\n-3,1\n"))).unwrap();
    assert_eq!(resultant["resultant"], "-2");
    for input in ["1,1\n", "1,1\n1,2\n1,3\n"] {
        assert_eq!(run(&["resultant"], input).status.code(), Some(2));
    }
}

#[test]
fn registry_help_and_completion_cover_every_command() {
    let metadata: Value = serde_json::from_str(&success(run(&["help", "--json"], ""))).unwrap();
    assert_eq!(metadata["schema"], "polytool.commands.v1");
    let commands = metadata["commands"].as_array().unwrap();
    assert!(commands.len() >= 35);
    let top = success(run(&["--help"], ""));
    let shells = ["bash", "zsh", "fish"].map(|shell| success(run(&["completions", shell], "")));
    for command in commands {
        let name = command["name"].as_str().unwrap();
        assert!(top.contains(name));
        let help = success(run(&[name, "--help"], ""));
        assert!(help.contains(name), "missing detailed help: {name}");
        for alias in command["aliases"].as_array().unwrap() {
            assert_eq!(success(run(&[alias.as_str().unwrap(), "--help"], "")), help);
        }
        assert!(shells.iter().all(|script| script.contains(name)));
    }
    assert_eq!(run(&["completions", "unknown"], "").status.code(), Some(2));
}
