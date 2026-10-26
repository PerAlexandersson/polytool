//! CLI: conditional theorem matches stay separate from finite-fit evidence.
use std::io::Write;
use std::process::{Command, Output, Stdio};

/// Rows of `1/(1 - x R(t x))` with `R = 1 + 3z + 2z^2` (labels 0..9).
const RENEWAL_ROWS: &[u8] = b"1\n1\n1,3\n1,6,2\n1,9,13\n1,12,33,12\n1,15,62,63,4\n\
1,18,100,180,66\n1,21,147,390,321,36\n1,24,203,720,985,360,8\n";

/// Rows of `Q_n = Q_(n-1) + n t Q_(max(0, n-2))`.
const PREDECESSOR_ROWS: &[u8] =
    b"1\n1,1\n1,3\n1,6,3\n1,10,15\n1,15,45,15\n1,21,105,105\n1,28,210,420,105\n";

fn run(args: &[&str], stdin: &[u8]) -> Output {
    let mut process = Command::new(env!("CARGO_BIN_EXE_polytool"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    process.stdin.take().unwrap().write_all(stdin).unwrap();
    process.wait_with_output().unwrap()
}

const SEARCH: [&str; 7] = [
    "recurrence",
    "--max-var-deg",
    "2",
    "--max-idx-deg",
    "0",
    "--max-rec-len",
    "3",
];

#[test]
fn default_json_is_unchanged_and_match_is_noted_on_stderr() {
    let output = run(&[&SEARCH[..], &["--json"]].concat(), RENEWAL_ROWS);
    assert!(output.status.success());
    let record: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(record["schema"], "polytool.recurrence.v1");
    assert!(record.get("theorem_matches").is_none());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("conditional theorem match weighted_renewal_window"));

    // The supplied record is checked by recurrence-theorems.
    let checked = run(
        &["recurrence-theorems", "--recurrence", "-", "--json"],
        &output.stdout,
    );
    assert!(checked.status.success());
    let report: serde_json::Value = serde_json::from_slice(&checked.stdout).unwrap();
    assert_eq!(report["schema"], "polytool.recurrence-theorems.v1");
    assert_eq!(report["finite_evidence"]["source"], "supplied_recurrence");
    let theorem = &report["theorem_matches"][0];
    assert_eq!(theorem["theorem"], "weighted_renewal_window");
    assert_eq!(theorem["window"]["max_step"], 1);
    assert_eq!(theorem["boundary"]["status"], "canonical_after_relabel");
}

#[test]
fn theorem_matches_flag_emits_separate_report() {
    let output = run(
        &[&SEARCH[..], &["--json", "--theorem-matches"]].concat(),
        RENEWAL_ROWS,
    );
    assert!(output.status.success());
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["schema"], "polytool.recurrence-theorems.v1");
    let plain = run(&[&SEARCH[..], &["--json"]].concat(), RENEWAL_ROWS);
    let plain: serde_json::Value = serde_json::from_slice(&plain.stdout).unwrap();
    assert_eq!(report["recurrence_json"], plain);
    assert_eq!(report["finite_evidence"]["source"], "adaptive_search_fit");
    assert_eq!(report["finite_evidence"]["status"], "finite_evidence_only");
    assert_eq!(report["finite_evidence"]["rows_compared"], 10);
    let theorem = &report["theorem_matches"][0];
    assert_eq!(theorem["status"], "conditional_theorem_match");
    assert_eq!(theorem["weighted_renewal"]["seed"][2], "2");
    assert_eq!(theorem["proof_status"], "paper_proof_not_lean_checked");

    let output = run(&["recurrence", "--theorem-matches"], PREDECESSOR_ROWS);
    assert!(output.status.success());
    let text = String::from_utf8(output.stdout).unwrap();
    assert!(text.contains("Conditional theorem match: fixed_delay_predecessor_window"));
    assert!(text.contains("Not matched weighted_renewal_window"));
}

#[test]
fn recurrence_theorems_exit_status() {
    let fibonacci = run(&["recurrence", "--json"], b"1\n1\n2\n3\n5\n8\n13\n");
    let none = run(
        &["recurrence-theorems", "--recurrence", "-"],
        &fibonacci.stdout,
    );
    assert_eq!(none.status.code(), Some(1));
    assert!(String::from_utf8(none.stdout)
        .unwrap()
        .contains("No conditional theorem match."));
    let invalid = run(&["recurrence-theorems", "--recurrence", "-"], b"{}");
    assert_eq!(invalid.status.code(), Some(2));
}
