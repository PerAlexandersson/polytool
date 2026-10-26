//! Regression tests for `polytool bench real-rootedness`.
//!
//! The suite carries an exact expected verdict for every check, so running it
//! in `--quick` mode doubles as a correctness regression test for exact
//! real-rootedness, strict interlacing, distinct-real-root counting and the
//! Ehrhart/h* conversion on positive-coefficient families.

use serde_json::Value;
use std::collections::BTreeSet;
use std::fs;
use std::process::{Command, Output};

fn run_bench(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_polytool"))
        .arg("bench")
        .arg("real-rootedness")
        .args(args)
        .output()
        .expect("run polytool bench real-rootedness")
}

fn stdout_of(output: &Output) -> String {
    String::from_utf8(output.stdout.clone()).expect("stdout is utf8")
}

fn stderr_of(output: &Output) -> String {
    String::from_utf8(output.stderr.clone()).expect("stderr is utf8")
}

#[test]
fn quick_suite_has_no_verdict_failures() {
    let output = run_bench(&["--quick"]);
    assert!(output.status.success(), "{}", stderr_of(&output));
    assert!(stderr_of(&output).contains("0 failures"));

    let stdout = stdout_of(&output);
    let mut lines = stdout.lines();
    let header: Vec<&str> = lines.next().expect("header").split('\t').collect();
    assert_eq!(
        header,
        [
            "case",
            "family",
            "kind",
            "degree",
            "coeff_bits",
            "expected",
            "verdict",
            "status",
            "repeat",
            "min_us",
            "median_us"
        ]
    );
    let mut families = BTreeSet::new();
    let mut kinds = BTreeSet::new();
    let mut rows = 0;
    for line in lines {
        let cells: Vec<&str> = line.split('\t').collect();
        assert_eq!(cells.len(), header.len(), "row: {line}");
        assert_eq!(
            cells[5], cells[6],
            "verdict differs from expectation: {line}"
        );
        assert_eq!(cells[7], "ok", "row: {line}");
        families.insert(cells[1].to_string());
        kinds.insert(cells[2].to_string());
        rows += 1;
    }
    assert!(rows >= 40, "expected the full quick suite, got {rows} rows");
    for family in ["eulerian", "narayana", "hstar", "counterexample"] {
        assert!(families.contains(family), "missing family {family}");
    }
    for kind in ["real_rooted", "interlacing", "root_count", "hstar_identity"] {
        assert!(kinds.contains(kind), "missing check kind {kind}");
    }
}

#[test]
fn quick_suite_covers_both_verdicts_of_real_rootedness() {
    let output = run_bench(&["--quick", "--format", "json"]);
    assert!(output.status.success(), "{}", stderr_of(&output));
    let report: Value = serde_json::from_str(&stdout_of(&output)).expect("json report");
    assert_eq!(report["schema"], "polytool.bench.real-rootedness.v1");
    assert_eq!(report["failures"], 0);
    let records = report["records"].as_array().expect("records");
    let verdicts: BTreeSet<&str> = records
        .iter()
        .filter(|r| r["kind"] == "real_rooted")
        .map(|r| r["verdict"].as_str().expect("verdict"))
        .collect();
    assert_eq!(verdicts, BTreeSet::from(["true", "false"]));
    // Interlacing covers true, false and "degrees do not admit the test".
    let interlacing: BTreeSet<&str> = records
        .iter()
        .filter(|r| r["kind"] == "interlacing")
        .map(|r| r["verdict"].as_str().expect("verdict"))
        .collect();
    assert_eq!(interlacing, BTreeSet::from(["true", "false", "none"]));
}

#[test]
fn degrees_option_scales_the_families() {
    let output = run_bench(&["--degrees", "20", "--only", "narayana", "--repeat", "1"]);
    assert!(output.status.success(), "{}", stderr_of(&output));
    let stdout = stdout_of(&output);
    assert!(stdout.contains("narayana_20\tnarayana\treal_rooted\t20\t"));
    assert!(stdout.contains("narayana_20\tnarayana\troot_count\t20\t"));
}

#[test]
fn baseline_comparison_detects_verdict_changes() {
    let dir = std::env::temp_dir();
    let id = std::process::id();
    let baseline_path = dir.join(format!("polytool-rr-baseline-{id}.json"));
    let tampered_path = dir.join(format!("polytool-rr-tampered-{id}.json"));

    let baseline = run_bench(&["--quick", "--json"]);
    assert!(baseline.status.success(), "{}", stderr_of(&baseline));
    fs::write(&baseline_path, &baseline.stdout).expect("write baseline");

    // The same run compared with itself passes and reports ratios.
    let path = baseline_path.to_str().expect("utf8 path");
    let same = run_bench(&["--quick", "--baseline", path]);
    assert!(same.status.success(), "{}", stderr_of(&same));
    assert!(stdout_of(&same)
        .lines()
        .next()
        .expect("header")
        .ends_with("ratio_vs_baseline"));

    // Flipping one saved verdict makes the comparison fail.
    let mut report: Value = serde_json::from_slice(&baseline.stdout).expect("json");
    report["records"][1]["verdict"] = Value::String("flipped".to_string());
    fs::write(&tampered_path, report.to_string()).expect("write tampered baseline");
    let tampered_arg = tampered_path.to_str().expect("utf8 path");
    let failed = run_bench(&["--quick", "--baseline", tampered_arg]);
    assert_eq!(failed.status.code(), Some(1));
    assert!(stderr_of(&failed).contains("differs from baseline"));

    fs::remove_file(&baseline_path).ok();
    fs::remove_file(&tampered_path).ok();
}

#[test]
fn invalid_options_are_rejected() {
    for args in [
        vec!["--degrees", "x"],
        vec!["--degrees", "2"],
        vec!["--max-slowdown", "2"],
        vec!["--baseline", "/nonexistent/polytool-baseline.json"],
        vec!["--unknown"],
    ] {
        let output = run_bench(&args);
        assert_eq!(output.status.code(), Some(1), "args: {args:?}");
        assert!(!stderr_of(&output).is_empty(), "args: {args:?}");
        assert!(
            !stdout_of(&output).contains("\tverdict\t"),
            "no report is printed for bad arguments: {args:?}"
        );
    }
}
