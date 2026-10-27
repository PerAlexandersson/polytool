//! CLI tests for recurrence representation search.

use std::io::Write;
use std::process::{Command, Stdio};

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
