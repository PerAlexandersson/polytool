//! A rule found with a search window longer than its lag keeps every row
//! before the window as an initial row, so the exported JSON regenerates
//! the input rows.

use std::io::Write;
use std::process::{Command, Stdio};

fn polytool() -> Command {
    Command::new(env!("CARGO_BIN_EXE_polytool"))
}

#[test]
fn exported_rule_keeps_rows_before_the_search_window() {
    // Row 0 is 1; from row 1 on, P(n) = (1 + t) P(n-1) with P(1) = t - 1.
    // The rule fails at n = 1, so it is found with window 2 and must keep
    // both P(0) and P(1).
    let rows = "1\n-1 1\n-1 0 1\n-1 -1 1 1\n-1 -2 0 2 1\n-1 -3 -2 2 3 1\n-1 -4 -5 0 5 4 1\n-1 -5 -9 -5 5 9 5 1\n";
    let mut child = polytool()
        .args(["recurrence", "--json"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .expect("polytool runs");
    child
        .stdin
        .take()
        .unwrap()
        .write_all(rows.as_bytes())
        .unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success());
    let json: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let initial = json["initial_polynomials"].as_array().unwrap();
    assert!(initial.len() >= 2, "{json}");
    assert_eq!(initial[1], serde_json::json!(["-1", "1"]));

    let dir = std::env::temp_dir().join(format!("polytool-export-seeds-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let file = dir.join("rule.json");
    std::fs::write(&file, serde_json::to_vec(&json).unwrap()).unwrap();
    let generated = polytool()
        .args(["recurrence-generate", "--recurrence"])
        .arg(&file)
        .args(["--rows", "8"])
        .output()
        .unwrap();
    assert!(generated.status.success());
    let text = String::from_utf8(generated.stdout).unwrap();
    let normalize = |s: &str| {
        s.lines()
            .map(|l| {
                l.replace(',', " ")
                    .split_whitespace()
                    .collect::<Vec<_>>()
                    .join(" ")
            })
            .filter(|l| !l.is_empty())
            .collect::<Vec<_>>()
    };
    assert_eq!(normalize(&text), normalize(rows));
    std::fs::remove_dir_all(&dir).ok();
}
