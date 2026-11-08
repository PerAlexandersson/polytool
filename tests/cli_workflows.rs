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
    let metadata: Value = serde_json::from_str(&success(run(&["help", "--json"], ""))).unwrap();
    for command in metadata["commands"].as_array().unwrap() {
        let command = command["name"].as_str().unwrap();
        let output = run(&[command, "--typo"], "");
        assert!(
            !output.status.success(),
            "{command} ignored an invalid option"
        );
        assert!(!output.stderr.is_empty(), "{command} needs a diagnostic");
        assert!(
            output.stdout.is_empty(),
            "{command} polluted stdout on an argument error"
        );
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
fn missing_arguments_and_invalid_help_paths_leave_stdout_empty() {
    for args in [
        vec!["sequence"],
        vec!["sequence", "unknown", "2"],
        vec!["oeis"],
        vec!["bench"],
        vec!["bench", "compare", "--typo"],
        vec!["bench", "interlacing", "--typo"],
        vec!["bench", "recurrence-fixtures", "--typo"],
        vec!["bench", "real-rootedness", "--typo"],
        vec!["pf-pencil"],
        vec!["pf-pencil", "--all-family-h"],
        vec!["help", "recurrence", "typo"],
        vec!["help", "oeis", "add", "typo"],
        vec!["help", "oeis", "unknown"],
        vec!["help", "unknown"],
        vec!["--version", "typo"],
        vec!["--help", "typo"],
    ] {
        let output = run(&args, "");
        assert_eq!(output.status.code(), Some(2), "{args:?}");
        assert!(output.stdout.is_empty(), "{args:?}");
        assert!(!output.stderr.is_empty(), "{args:?}");
    }
}

#[test]
fn input_errors_locate_physical_lines_without_changing_row_indices() {
    let input = "# input\r\n1,2\r\n\r\n  # comment\r\n1,,3\r\n";
    for command in ["real-rooted", "properties", "interlacing", "recurrence"] {
        let output = run(&[command, "--json"], input);
        assert_eq!(output.status.code(), Some(2), "{command}");
        assert!(output.stdout.is_empty(), "{command}");
        let message = String::from_utf8(output.stderr).unwrap();
        assert!(message.contains("line 5"), "{command}: {message}");
        if command == "recurrence" {
            assert!(message.contains("coefficient 2"), "{message}");
        }
    }
    let output = run(&["family-check", "--json"], input);
    let report: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["parse_errors"][0]["index"], 1);
    assert!(report["parse_errors"][0]["error"]
        .as_str()
        .unwrap()
        .contains("line 5"));

    for command in ["recurrence", "bernstein-expansion"] {
        let output = run(&[command], "# input\n1, 3/0\n");
        assert_eq!(output.status.code(), Some(2));
        assert!(output.stdout.is_empty());
        let message = String::from_utf8(output.stderr).unwrap();
        assert!(
            message.contains("line 2") && message.contains("coefficient 2"),
            "{message}"
        );
        assert!(message.contains("zero denominator"), "{message}");
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

#[test]
fn json_only_commands_accept_explicit_json_without_changing_reports() {
    let jobs = [
        ("root-enclosures", r#"{"polynomials":[["-1","0","1"]]}"#),
        (
            "recurrence-transform",
            r#"{"oeis_id":"A008292","row_count":4,"transformation":{"kind":"f_polynomial","schedule":{"first_index":1,"first_degree":0,"step":1}}}"#,
        ),
        (
            "recurrence-implicit",
            r#"{"oeis_id":"A007318","row_count":4,"transformation":{"kind":"diagonal","weight":{"kind":"factorial_power","power":1}}}"#,
        ),
        (
            "recurrence-search-representations",
            r#"{"oeis_id":"A007318","row_count":4,"limits":{"max_depth":1},"transformations":[{"kind":"diagonal","weight":{"kind":"factorial_power","power":1}}]}"#,
        ),
    ];
    for (command, request) in jobs {
        let default = success(run(&[command], request));
        assert_eq!(success(run(&[command, "--json"], request)), default);
        let compact = success(run(&[command, "--json", "--compact"], request));
        assert_eq!(compact.lines().count(), 1);
        assert_eq!(
            serde_json::from_str::<Value>(&default).unwrap(),
            serde_json::from_str::<Value>(&compact).unwrap()
        );
    }
}

#[test]
fn enclosure_replay_rejects_tampering_in_library_and_cli() {
    use polytool::recurrence::{verify_root_enclosures, RootEnclosureReport};
    let request =
        r#"{"polynomials":[["0","0","1","2","1"],["1","0","1"],[],["3/2"]],"max_width":"1/100"}"#;
    let report: Value = serde_json::from_str(&success(run(&["root-enclosures"], request))).unwrap();
    let decoded: RootEnclosureReport = serde_json::from_value(report.clone()).unwrap();
    verify_root_enclosures(&decoded).unwrap();
    let verified: Value = serde_json::from_str(&success(run(
        &["root-enclosures", "--verify", "--json", "--compact"],
        &report.to_string(),
    )))
    .unwrap();
    assert_eq!(verified["verified"], true);
    assert_eq!(verified["polynomials"], 4);
    for (pointer, value) in [
        ("/schema", json!("unknown.schema")),
        ("/items/0/roots/0/multiplicity", json!(7)),
        ("/items/0/roots/0/lower", json!("-99")),
        ("/items/1/real_rooted", json!(true)),
        ("/request/polynomials/0/0", json!("1")),
        ("/request/max_width", json!("-1")),
    ] {
        let mut tampered = report.clone();
        *tampered.pointer_mut(pointer).unwrap() = value;
        let decoded: RootEnclosureReport = serde_json::from_value(tampered.clone()).unwrap();
        assert!(verify_root_enclosures(&decoded).is_err(), "{pointer}");
        let output = run(&["root-enclosures", "--verify"], &tampered.to_string());
        assert_eq!(output.status.code(), Some(2), "{pointer}");
        assert!(output.stdout.is_empty());
    }
}

#[test]
fn veronese_search_reuses_report_wrappers_and_explicit_menu_selection() {
    let record = json!({"schema":"polytool.recurrence.v1", "first_index":0,
    "initial_polynomials":[["1"],["1"]], "recurrence":{"terms":[
        {"offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["2"]]}},
        {"offset":2,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["-1","1"]]}}
    ]}});
    let args = ["recurrence-veronese-lifts", "--json", "--parent-rows", "4"];
    let expected = success(run(&args, &record.to_string()));
    let wrapper = json!({"schema":"polytool.recurrence-theorems.v1", "recurrence_json":record});
    assert_eq!(success(run(&args, &wrapper.to_string())), expected);
    let menu = json!({"schema":"polytool.recurrence-representations.v1", "representations":[
        {"recurrence_json":record}, {"recurrence_json":record}
    ]});
    assert_eq!(run(&args, &menu.to_string()).status.code(), Some(2));
    let selected = [&args[..], &["--candidate", "1"]].concat();
    assert_eq!(success(run(&selected, &menu.to_string())), expected);
}

#[test]
fn scout_rejects_invalid_options_before_reading_a_symbol() {
    for args in [
        vec!["--format", "bogus"],
        vec!["--box", "NaN", "1", "-1", "1"],
        vec!["--tol", "inf"],
        vec!["--min-imag", "-inf"],
        vec!["--min-imag", "-1"],
        vec!["--tol", "0"],
        vec!["--max-iter", "0"],
        vec!["--box", "-1e308", "1e308", "-1", "1"],
    ] {
        let args = [&["bkw-scout"][..], &args[..]].concat();
        let output = run(&args, "not a symbol");
        assert_eq!(output.status.code(), Some(2), "{args:?}");
        assert!(output.stdout.is_empty());
        let error = String::from_utf8(output.stderr).unwrap();
        assert!(
            !error.contains("parse BKW symbol"),
            "late validation: {args:?}: {error}"
        );
    }
}

#[test]
fn scout_json_with_mathematica_is_one_valid_document() {
    let args = [
        "bkw-scout",
        "--grid",
        "1",
        "--include-real-axis",
        "--no-refine",
        "--json",
    ];
    let plain: Value = serde_json::from_str(&success(run(&args, "1;0;1"))).unwrap();
    let with_formula = [&args[..], &["--mathematica"]].concat();
    let mut augmented: Value = serde_json::from_str(&success(run(&with_formula, "1;0;1"))).unwrap();
    assert!(augmented["mathematica"]
        .as_str()
        .unwrap()
        .contains("Reduce"));
    augmented.as_object_mut().unwrap().remove("mathematica");
    assert_eq!(augmented, plain);
    assert!(!plain["candidates"].as_array().unwrap().is_empty());
}
