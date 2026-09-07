use std::io::Write;
use std::process::{Command, Output, Stdio};

fn run_polytool(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_polytool"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn polytool");

    child
        .stdin
        .as_mut()
        .expect("polytool stdin")
        .write_all(input.as_bytes())
        .expect("write polytool stdin");
    child.wait_with_output().expect("wait for polytool")
}

#[test]
fn version_flag_prints_the_exact_build_version() {
    let output = run_polytool(&["--version"], "");

    assert!(output.status.success());
    assert!(output.stderr.is_empty());
    assert_eq!(
        String::from_utf8(output.stdout).expect("version output is UTF-8"),
        format!("{}\n", polytool::version::build_version())
    );
}

const FIBONACCI_ROWS: &str = "1\n1\n2\n3\n5\n";
const GEOMETRIC_ROWS: &str = "1\n2\n4\n8\n16\n";

fn stdout_json(output: &Output) -> serde_json::Value {
    serde_json::from_slice(&output.stdout).expect("polytool stdout is valid JSON")
}

#[test]
fn zero_candidate_budget_reports_exhaustion_without_work() {
    let output = run_polytool(
        &["recurrence", "--max-candidates", "0", "--format", "json"],
        FIBONACCI_ROWS,
    );

    assert_eq!(output.status.code(), Some(3));
    let json = stdout_json(&output);
    assert_eq!(json["status"], "budget_exhausted");
    assert_eq!(json["candidates_considered"], 0);
    assert_eq!(json["candidates_tried"], 0);
    assert_eq!(json["max_candidates"], 0);
}

#[test]
fn small_candidate_budget_reports_exhaustion_after_exact_count() {
    let output = run_polytool(
        &["recurrence", "--max-candidates", "1", "--format", "json"],
        FIBONACCI_ROWS,
    );

    assert_eq!(output.status.code(), Some(3));
    let json = stdout_json(&output);
    assert_eq!(json["status"], "budget_exhausted");
    assert_eq!(json["candidates_considered"], 1);
    assert_eq!(json["max_candidates"], 1);
}

#[test]
fn exact_budget_boundary_reports_complete_search_space_failure() {
    let output = run_polytool(
        &[
            "recurrence",
            "--max-candidates",
            "1",
            "--min-rec-len",
            "1",
            "--max-rec-len",
            "1",
            "--max-var-deg",
            "0",
            "--max-idx-deg",
            "0",
            "--max-diff-deg",
            "0",
            "--format",
            "json",
        ],
        FIBONACCI_ROWS,
    );

    assert!(output.status.success());
    let json = stdout_json(&output);
    assert_eq!(json["status"], "not_found");
    assert_eq!(json["candidates_considered"], 1);
    assert_eq!(json["max_candidates"], 1);
}

#[test]
fn recurrence_can_succeed_on_the_last_budgeted_candidate() {
    let output = run_polytool(&["recurrence", "--max-candidates", "1"], GEOMETRIC_ROWS);

    assert!(output.status.success());
    assert_eq!(output.stdout, b"P(n) = 2 P(n-1)\n");
    assert!(String::from_utf8(output.stderr)
        .expect("stderr is UTF-8")
        .contains("1 candidates tried"));
}

#[test]
fn omitted_candidate_budget_preserves_unbounded_cli_behavior() {
    let output = run_polytool(&["recurrence"], GEOMETRIC_ROWS);

    assert!(output.status.success());
    assert_eq!(output.stdout, b"P(n) = 2 P(n-1)\n");
    assert!(!String::from_utf8(output.stderr)
        .expect("stderr is UTF-8")
        .contains("budget"));
}
