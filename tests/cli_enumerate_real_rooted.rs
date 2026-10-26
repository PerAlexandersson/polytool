//! CLI tests for `polytool enumerate-real-rooted`.

use std::process::{Command, Output};

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_polytool"))
        .arg("enumerate-real-rooted")
        .args(args)
        .output()
        .expect("run polytool")
}

fn stdout(output: &Output) -> String {
    String::from_utf8(output.stdout.clone()).expect("stdout is utf8")
}

fn json(output: &Output) -> serde_json::Value {
    serde_json::from_str(&stdout(output)).expect("valid JSON")
}

#[test]
fn text_output_states_conventions_and_lists_ascending_rows() {
    let output = run(&["--degree", "3", "--leading-desc", "1,0,-3"]);
    assert!(output.status.success());
    let text = stdout(&output);
    assert!(text.contains("descending a_0, a_1"));
    assert!(text.contains("output: ascending coefficients"));
    assert!(text.contains("# status: complete (5 polynomials)"));
    assert!(text.contains("max_nodes=1000000 max_results=10000 max_coefficient_bits=none"));
    let rows: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(
        rows,
        vec![
            "-2, -3, 0, 1",
            "-1, -3, 0, 1",
            "0, -3, 0, 1",
            "1, -3, 0, 1",
            "2, -3, 0, 1"
        ]
    );
}

#[test]
fn json_output_is_exact_and_reports_statistics() {
    let output = run(&[
        "--degree",
        "4",
        "--leading-desc",
        "1 0 -4",
        "--format",
        "json",
    ]);
    assert!(output.status.success());
    let value = json(&output);
    assert_eq!(value["status"], "complete");
    assert_eq!(value["limit_reached"], serde_json::Value::Null);
    assert_eq!(value["output_coefficient_order"], "ascending");
    assert_eq!(
        value["leading_coefficients_descending"],
        serde_json::json!(["1", "0", "-4"])
    );
    assert_eq!(value["limits"]["max_nodes"], 1_000_000);
    assert_eq!(
        value["limits"]["max_coefficient_bits"],
        serde_json::Value::Null
    );
    let polys = value["polynomials"].as_array().unwrap();
    assert_eq!(value["count"], polys.len());
    assert!(polys.contains(&serde_json::json!(["4", "0", "-4", "0", "1"])));
    assert_eq!(value["stats"]["results_emitted"], polys.len());
    assert_eq!(value["stats"]["derivative_levels"], 2);
}

#[test]
fn json_keeps_integers_beyond_double_precision_exact() {
    // x^3 - 3 * 10^20 x + c is real-rooted for |c| <= 2 * 10^30.
    let output = run(&[
        "--degree",
        "3",
        "--leading-desc",
        "1,0,-300000000000000000000",
        "--max-results",
        "2",
        "--json",
    ]);
    assert!(output.status.success());
    let value = json(&output);
    assert_eq!(value["status"], "resource_limit");
    assert_eq!(value["limit_reached"], "max_results");
    assert_eq!(
        value["polynomials"],
        serde_json::json!([
            [
                "-2000000000000000000000000000000",
                "-300000000000000000000",
                "0",
                "1"
            ],
            [
                "-1999999999999999999999999999999",
                "-300000000000000000000",
                "0",
                "1"
            ]
        ])
    );
}

#[test]
fn empty_complete_result_is_distinct_from_limits_and_errors() {
    // a_1^2 < 2 a_0 a_2: no real-rooted completion.
    let output = run(&["--degree", "4", "--leading-desc", "1,0,1", "--json"]);
    assert!(output.status.success());
    let value = json(&output);
    assert_eq!(value["status"], "complete");
    assert_eq!(value["count"], 0);

    let output = run(&[
        "--degree",
        "4",
        "--leading-desc",
        "1,0,-4",
        "--max-nodes",
        "0",
        "--json",
    ]);
    assert!(output.status.success());
    let value = json(&output);
    assert_eq!(value["status"], "resource_limit");
    assert_eq!(value["limit_reached"], "max_nodes");
    assert_eq!(value["count"], 0);

    let output = run(&[
        "--degree",
        "4",
        "--leading-desc",
        "1,0",
        "--max-results",
        "none",
        "--json",
    ]);
    assert_eq!(output.status.code(), Some(2));
    let value = json(&output);
    assert_eq!(value["status"], "invalid_input");
    assert!(value["error"]
        .as_str()
        .unwrap()
        .contains("leading coefficients"));
}

#[test]
fn invalid_arguments_are_rejected() {
    for args in [
        vec!["--leading-desc", "1,0,-3"],
        vec!["--degree", "3"],
        vec!["--degree", "3", "--leading-desc", "1,x,-3"],
        vec!["--degree", "3", "--leading-desc", "0,0,-3"],
        vec!["--degree", "3", "--leading-desc", "1,0,-3", "--bogus"],
        vec![
            "--degree",
            "3",
            "--leading-desc",
            "1,0,-3",
            "--max-nodes",
            "-1",
        ],
    ] {
        let output = run(&args);
        assert_eq!(output.status.code(), Some(2), "{args:?}");
        assert!(!output.stderr.is_empty());
    }
}

#[test]
fn help_mentions_coefficient_orders() {
    let output = run(&["--help"]);
    assert!(output.status.success());
    let text = stdout(&output);
    assert!(text.contains("DESCENDING"));
    assert!(text.contains("ASCENDING"));
    assert!(text.contains("resource_limit"));
}
