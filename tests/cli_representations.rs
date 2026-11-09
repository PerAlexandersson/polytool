//! CLI tests for recurrence representation search.

use std::io::Write;
use std::process::{Command, Stdio};

mod common;
use common::run_polytool as run;

#[test]
fn unsupported_menu_flags_fail_before_reading_rows() {
    for flag in ["--theorem-matches", "--python"] {
        let output = run(
            &["recurrence", "--representations", flag],
            "not a polynomial",
        );
        assert_eq!(output.status.code(), Some(2));
        assert!(output.stdout.is_empty());
        let message = String::from_utf8(output.stderr).unwrap();
        assert!(message.contains("--representations"));
        assert!(!message.contains("Invalid polynomial"));
        if flag == "--theorem-matches" {
            assert!(message.contains("recurrence-theorems --recurrence menu.json --candidate N"));
        }
    }
}

#[test]
fn text_menu_indices_reuse_the_same_json_candidates() {
    let input = run(
        &[
            "oeis", "generate", "A171608", "--rows", "40", "--format", "rows",
        ],
        "",
    );
    assert!(input.status.success());
    let input = String::from_utf8(input.stdout).unwrap();
    let args = [
        "recurrence",
        "--representations",
        "--max-rec-len",
        "3",
        "--max-diff-deg",
        "2",
        "--max-var-deg",
        "3",
        "--max-idx-deg",
        "1",
        "--denominator",
        "--max-denom-var-deg",
        "0",
        "--max-denom-idx-deg",
        "1",
        "--max-candidates",
        "2000",
    ];
    let text = run(&args, &input);
    assert!(text.status.success());
    let text = String::from_utf8(text.stdout).unwrap();
    let json = run(&[&args[..], &["--json"]].concat(), &input);
    assert!(json.status.success());
    let menu = String::from_utf8(json.stdout).unwrap();
    let report: serde_json::Value = serde_json::from_str(&menu).unwrap();
    let candidates = report["representations"].as_array().unwrap();
    assert!(candidates.len() >= 3);
    assert!(text.contains("lag  deriv") && text.contains("fit  check  eqs"));
    assert!(text.contains("not all-rank proofs"));
    for (index, candidate) in candidates.iter().enumerate() {
        assert!(text.contains(&format!(
            "[{index}] {}",
            candidate["recurrence"].as_str().unwrap()
        )));
        let checked = run(
            &[
                "recurrence-theorems",
                "--recurrence",
                "-",
                "--candidate",
                &index.to_string(),
                "--json",
            ],
            &menu,
        );
        // A missing theorem match is status 1; both outcomes must select the same record.
        assert!(matches!(checked.status.code(), Some(0 | 1)));
        let checked: serde_json::Value = serde_json::from_slice(&checked.stdout).unwrap();
        assert_eq!(checked["recurrence_json"], candidate["recurrence_json"]);
    }

    let empty = run(&[&args[..], &["--max-candidates", "0"]].concat(), &input);
    assert_eq!(empty.status.code(), Some(3));
    let text = String::from_utf8(empty.stdout).unwrap();
    assert!(text.contains("BudgetExhausted"));
    assert!(text.contains("No representations found within the inspected bounds"));
}

#[test]
fn cli_menu_json_has_replayable_candidates_and_partial_exit_status() {
    for limit in ["0", "1"] {
        let mut process = Command::new(env!("CARGO_BIN_EXE_polytool"))
            .args([
                "recurrence",
                "--representations",
                "--json",
                "--max-rec-len",
                "1",
                "--max-diff-deg",
                "0",
                "--max-var-deg",
                "0",
                "--max-idx-deg",
                "0",
                "--max-candidates",
                limit,
            ])
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .spawn()
            .unwrap();
        process
            .stdin
            .take()
            .unwrap()
            .write_all(b"1\n2\n4\n8\n16\n32\n")
            .unwrap();
        let output = process.wait_with_output().unwrap();
        let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
        assert_eq!(report["schema"], "polytool.recurrence-representations.v1");
        if limit == "0" {
            assert_eq!(output.status.code(), Some(3));
            assert_eq!(report["status"], "budget_exhausted");
        } else {
            assert!(output.status.success());
            let record: polytool::recurrence::RecurrenceJson =
                serde_json::from_value(report["representations"][0]["recurrence_json"].clone())
                    .unwrap();
            let (rec, first, initial) = record.to_recurrence_parts().unwrap();
            assert_eq!(
                rec.generate_rows_rational(&initial, first, 8).unwrap()[7][0].to_string(),
                "128"
            );
        }
    }
}
