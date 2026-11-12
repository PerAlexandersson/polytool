//! CLI discovery, explicit menu selection, certificate replay and usage failures.

mod common;
use common::run_polytool as run;

fn input() -> serde_json::Value {
    serde_json::json!({
        "schema": "polytool.recurrence.v1", "first_index": 0,
        "initial_polynomials": [["1"]],
        "recurrence": {"terms": [{"offset": 1, "deriv_order": 0, "sign": "none",
            "coeff": {"coeffs": [["1", "1"], ["1"]]}}]}
    })
}

#[test]
fn discover_and_replay_then_reject_tampering() {
    let result = run(
        &["recurrence-prove", "--recurrence", "-", "--json"],
        &input().to_string(),
    );
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let mut report: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(report["status"], "proved");
    assert_eq!(report["attempts"][0]["theorem"], "product_recurrence");
    assert_eq!(report["lean_replay"], "not_attempted");
    let mut old_report = report.clone();
    old_report
        .as_object_mut()
        .unwrap()
        .remove("trivially_real_rooted_if_rule_holds");
    assert!(run(
        &["recurrence-prove", "--verify", "--json"],
        &old_report.to_string()
    )
    .status
    .success());
    let mut altered = report.clone();
    altered["trivially_real_rooted_if_rule_holds"] = serde_json::json!(false);
    assert_eq!(
        run(&["recurrence-prove", "--verify"], &altered.to_string())
            .status
            .code(),
        Some(2)
    );
    let verified = run(
        &["recurrence-prove", "--verify", "--json"],
        &report.to_string(),
    );
    assert!(verified.status.success());
    report["attempts"][0]["conclusion"]["simple_roots"] = serde_json::json!(true);
    let bad = run(&["recurrence-prove", "--verify"], &report.to_string());
    assert_eq!(bad.status.code(), Some(2));
    assert!(bad.stdout.is_empty());
}

#[test]
fn explicit_menu_selection_and_no_proof_status() {
    let mut bad = input();
    bad["recurrence"]["terms"][0]["coeff"]["coeffs"] = serde_json::json!([["1", "0", "1"]]);
    let menu = serde_json::json!({"schema": "polytool.recurrence-representations.v1",
        "representations": [{"recurrence_json": input()}, {"recurrence_json": bad}]});
    let ambiguous = run(
        &["recurrence-prove", "--recurrence", "-"],
        &menu.to_string(),
    );
    assert_eq!(ambiguous.status.code(), Some(2));
    let yes = run(
        &["recurrence-prove", "--recurrence", "-", "--candidate", "0"],
        &menu.to_string(),
    );
    assert!(yes.status.success());
    assert!(String::from_utf8_lossy(&yes.stdout).contains("ProductRecurrence: Proved"));
    let no = run(
        &[
            "recurrence-prove",
            "--recurrence",
            "-",
            "--candidate",
            "1",
            "--json",
        ],
        &menu.to_string(),
    );
    assert_eq!(no.status.code(), Some(1));
    let report: serde_json::Value = serde_json::from_slice(&no.stdout).unwrap();
    assert_eq!(report["status"], "not_proved");
    assert!(run(&["recurrence-prove", "--verify"], &report.to_string())
        .status
        .success());
}

#[test]
fn help_and_invalid_flags_do_not_consume_input() {
    assert!(run(&["help", "recurrence-prove"], "").status.success());
    for args in [
        vec!["recurrence-prove", "--unknown"],
        vec!["recurrence-prove", "--verify", "--recurrence", "-"],
        vec![
            "recurrence-prove",
            "--max-start-delay",
            "17",
            "--recurrence",
            "-",
        ],
    ] {
        let result = run(&args, "invalid JSON");
        assert_eq!(result.status.code(), Some(2));
        assert!(result.stdout.is_empty());
        assert!(!String::from_utf8_lossy(&result.stderr).contains("parse recurrence JSON"));
    }
}

#[test]
fn the_record_defaults_to_the_command_input() {
    // Without --recurrence the record is read from stdin, or from the file
    // named by the universal --input option, for every command whose main
    // input is a recurrence record.
    let record = input().to_string();
    let path = std::env::temp_dir().join(format!(
        "polytool-cli-recurrence-input-{}.json",
        std::process::id()
    ));
    std::fs::write(&path, &record).unwrap();
    let file = path.to_str().unwrap();
    for (command, extra) in [
        ("recurrence-prove", vec!["--json"]),
        ("recurrence-degree", vec!["--json"]),
        ("recurrence-theorems", vec!["--json"]),
        ("recurrence-generate", vec!["--count", "3"]),
    ] {
        let from_stdin = run(&[&[command][..], &extra].concat(), &record);
        let from_file = run(&[&[command, "--input", file][..], &extra].concat(), "");
        let explicit = run(&[&[command, "--recurrence", file][..], &extra].concat(), "");
        assert!(!from_stdin.stdout.is_empty(), "{command}");
        assert_eq!(from_stdin.stdout, from_file.stdout, "{command}");
        assert_eq!(from_stdin.stdout, explicit.stdout, "{command}");
        assert_eq!(
            from_stdin.status.code(),
            from_file.status.code(),
            "{command}"
        );
    }
    std::fs::remove_file(&path).unwrap();
}

#[test]
fn zero_residue_classes_are_reported_by_the_cli() {
    // P_n = x^2 P_(n-2), P_0 = 1, P_1 = 0: every odd row is zero.
    let record = serde_json::json!({
        "schema": "polytool.recurrence.v1", "first_index": 0,
        "initial_polynomials": [["1"], ["0"]],
        "recurrence": {"terms": [{"offset": 2, "deriv_order": 0, "sign": "none",
            "coeff": {"coeffs": [["0", "0", "1"]]}}]}
    });
    let result = run(&["recurrence-prove", "--json"], &record.to_string());
    assert!(result.status.success());
    let report: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(report["status"], "proved");
    assert!(report.get("real_rooted_from_index").is_none());
    assert_eq!(
        report["real_rooted_outside_zero_residue_classes"],
        serde_json::json!({"from_index": 0,
            "zero_residue_classes": {"period": 2, "residues": [1]}})
    );
    let text = run(&["recurrence-prove"], &record.to_string());
    let text = String::from_utf8_lossy(&text.stdout);
    assert!(
        !text.contains("Every row P_n with n >= 0 is nonzero"),
        "{text}"
    );
    assert!(
        text.contains("Zero residue classes: every row P_n with n >= 0 and n mod 2 in [1] is zero")
    );
}

#[test]
fn conditional_product_classification_is_visible_in_json() {
    let supplied = run(
        &["recurrence-prove", "--recurrence", "-", "--json"],
        &input().to_string(),
    );
    let supplied: serde_json::Value = serde_json::from_slice(&supplied.stdout).unwrap();
    assert_eq!(supplied["trivially_real_rooted"], true);
    assert_eq!(supplied["trivially_real_rooted_if_rule_holds"], true);
    let mut nonreal = input();
    nonreal["recurrence"]["terms"][0]["coeff"]["coeffs"] = serde_json::json!([["1", "0", "1"]]);
    let result = run(
        &["recurrence-prove", "--recurrence", "-", "--json"],
        &nonreal.to_string(),
    );
    let result: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(result["trivially_real_rooted"], false);
    assert_eq!(result["trivially_real_rooted_if_rule_holds"], false);
    let fit = run(
        &["recurrence", "--json"],
        "1\n1 1\n1 2 1\n1 3 3 1\n1 4 6 4 1\n1 5 10 10 5 1\n",
    );
    assert!(fit.status.success());
    let result = run(
        &["recurrence-prove", "--recurrence", "-", "--json"],
        std::str::from_utf8(&fit.stdout).unwrap(),
    );
    let result: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(result["trivially_real_rooted"], false);
    assert_eq!(result["trivially_real_rooted_if_rule_holds"], true);
    for command in ["recurrence-prove", "recurrence-theorems"] {
        let output = run(
            &[command, "--recurrence", "-"],
            std::str::from_utf8(&fit.stdout).unwrap(),
        );
        assert!(String::from_utf8_lossy(&output.stdout)
            .contains("Trivially real-rooted if the rule holds: yes."));
        let plain = run(&[command, "--recurrence", "-"], &input().to_string());
        assert!(!String::from_utf8_lossy(&plain.stdout)
            .contains("Trivially real-rooted if the rule holds:"));
    }
    let theorem = run(
        &["recurrence", "--json", "--theorem-matches"],
        "1\n1 1\n1 2 1\n1 3 3 1\n1 4 6 4 1\n1 5 10 10 5 1\n",
    );
    let theorem: serde_json::Value = serde_json::from_slice(&theorem.stdout).unwrap();
    assert!(theorem.get("trivially_real_rooted").is_none());
    assert_eq!(theorem["trivially_real_rooted_if_rule_holds"], true);
}

#[test]
fn structural_gamma_chain_text_json_and_replay() {
    let entry = polytool::oeis::by_id("A008303").unwrap();
    let (rule, seeds) = entry.recurrence_parts().unwrap();
    let record = polytool::recurrence::RecurrenceJson::from_recurrence_rational(
        &rule,
        entry.recurrence_first_index,
        &seeds,
        None,
    );
    let input = serde_json::to_string(&record).unwrap();
    let output = run(&["recurrence-prove", "--recurrence", "-", "--json"], &input);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let mut report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["closure_lifts"]["attempts"][0]["stage"], "gamma");
    assert!(run(
        &["recurrence-prove", "--verify", "--json"],
        &report.to_string()
    )
    .status
    .success());
    let text = run(&["recurrence-prove", "--recurrence", "-"], &input);
    let text = String::from_utf8(text.stdout).unwrap();
    for step in [
        "All-index identity",
        "Computed parent proof",
        "Closure theorem",
        "Nonzero rows",
    ] {
        assert!(text.contains(step), "{text}");
    }
    report["closure_lifts"]["attempts"][0]["steps"][0]["statement"] = "finite rows only".into();
    assert_eq!(
        run(&["recurrence-prove", "--verify"], &report.to_string())
            .status
            .code(),
        Some(2)
    );
}
