use std::fs;
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
fn real_rooted_accepts_bigint_coefficients() {
    let output = run_polytool(&["real-rooted"], "-1000000000000000000000000000000,1\n");

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("stdout is utf8");
    assert!(stdout.contains("-1000000000000000000000000000000 + t: real-rooted"));
}

#[test]
fn malformed_polynomial_input_is_fatal() {
    let output = run_polytool(&["real-rooted"], "1,1\nnot-a-polynomial\n1,2,1\n");

    assert!(!output.status.success());
    assert_eq!(output.status.code(), Some(2));
    assert!(output.stdout.is_empty());
    let stderr = String::from_utf8(output.stderr).expect("stderr is utf8");
    assert!(stderr.contains("Invalid polynomial input"));
}

#[test]
fn huge_exponents_are_rejected_without_panicking() {
    let output = run_polytool(&["real-rooted"], &format!("t^{}\n", usize::MAX));

    assert_eq!(output.status.code(), Some(2));
    let stderr = String::from_utf8(output.stderr).expect("stderr is utf8");
    assert!(stderr.contains("too large"));
}

#[test]
fn cli_rejects_oversized_standard_input() {
    const LIMIT: usize = 16 * 1024 * 1024;
    let mut input = "1,".repeat(LIMIT / 2);
    input.push_str("1\n");
    let output = run_polytool(&["real-rooted"], &input);

    assert_eq!(output.status.code(), Some(2));
    let stderr = String::from_utf8(output.stderr).expect("stderr is utf8");
    assert!(stderr.contains("CLI input limit"));
}

#[test]
fn nonintegral_ehrhart_conversion_is_a_clean_error() {
    let output = run_polytool(&["ehrhart-to-hstar"], "1/2\n");

    assert_eq!(output.status.code(), Some(2));
    assert!(output.stdout.is_empty());
    let stderr = String::from_utf8(output.stderr).expect("stderr is utf8");
    assert!(stderr.contains("not an integer"));
}

#[test]
fn explicit_hstar_dimension_reports_degree_violation() {
    let output = run_polytool(
        &["hstar-inequalities", "--dimension", "1", "--json"],
        "1,0,1\n",
    );

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("stdout is utf8");
    let json: serde_json::Value = serde_json::from_str(&stdout).expect("valid hstar JSON");
    assert_eq!(
        json["items"][0]["hstar"],
        serde_json::json!(["1", "0", "1"])
    );
    assert_eq!(json["items"][0]["degree"], 2);
    assert_eq!(json["items"][0]["all_applicable_hold"], false);
}

#[test]
fn gamma_expansion_json_accepts_bigint_coefficients() {
    let output = run_polytool(
        &["gamma-expansion", "--json"],
        "1,100000000000000000000,1\n",
    );

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("stdout is utf8");
    assert!(stdout.contains("\"coefficients\":[\"1\",\"100000000000000000000\",\"1\"]"));
    assert!(stdout.contains("\"gamma\":[\"1\",\"99999999999999999998\"]"));
    assert!(stdout.contains("\"expansion\":\"(1+t)^2 + 99999999999999999998 t\""));
}

#[test]
fn sequence_json_generates_bigint_coefficients() {
    let output = run_polytool(&["sequence", "chebyshev-t", "64", "--json"], "");

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("stdout is utf8");
    assert!(stdout.contains("\"degree\":64"));
    assert!(stdout.contains("\"9223372036854775808\""));
}

#[test]
fn resultant_and_discriminant_accept_bigint_coefficients() {
    let resultant = run_polytool(
        &["resultant"],
        "-100000000000000000000,1\n-100000000000000000001,1\n",
    );
    assert!(resultant.status.success());
    let stdout = String::from_utf8(resultant.stdout).expect("stdout is utf8");
    assert!(stdout.contains("= -1"));

    let discriminant = run_polytool(&["discriminant"], "1,0,1000000000000000000000000000000\n");
    assert!(discriminant.status.success());
    let stdout = String::from_utf8(discriminant.stdout).expect("stdout is utf8");
    assert!(stdout.contains("-4000000000000000000000000000000"));
}

#[test]
fn pf_pencil_a036969_degree_five_reports_real_rooted_checks() {
    let output = run_polytool(&["pf-pencil", "--case", "A036969", "--degree", "5"], "");

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("stdout is utf8");
    assert!(stdout.contains("A036969 d=5"));
    assert!(stdout.contains("J_alpha,d(t) = "));
    assert!(stdout.contains("t J_beta,d(t) = "));
    assert!(stdout.contains("lambda=0: "));
    assert!(stdout.contains("lambda=1: "));
    assert!(stdout.contains("lambda=10: "));
    assert!(stdout.contains("lambda=100: "));
    assert!(stdout.contains("finite evidence ok: true"));
    assert!(!stdout.contains("NOT real-rooted"));
}

#[test]
fn pf_pencil_all_family_h_json_reports_all_cases() {
    let output = run_polytool(
        &["pf-pencil", "--all-family-h", "--max-degree", "3", "--json"],
        "",
    );

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("stdout is utf8");
    let json: serde_json::Value = serde_json::from_str(&stdout).expect("valid pf-pencil JSON");
    assert_eq!(json["schema"], "polytool.pf-pencil.v1");
    assert_eq!(json["overall_finite_evidence_ok"], true);
    let items = json["items"].as_array().expect("items is an array");
    assert_eq!(items.len(), 45);

    let cases = items
        .iter()
        .map(|item| item["case"].as_str().expect("case string"))
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(cases.len(), 15);
    assert!(cases.contains("A036969"));
    assert!(cases.contains("A198204"));
}

#[test]
fn pf_pencil_denominator_clearing_cases_do_not_panic() {
    let a080248 = run_polytool(&["pf-pencil", "--case", "A080248", "--degree", "5"], "");
    assert!(a080248.status.success());
    let stdout = String::from_utf8(a080248.stdout).expect("stdout is utf8");
    assert!(stdout.contains("common denominator cleared: 2"));
    assert!(stdout.contains("finite evidence ok: true"));

    let a198204 = run_polytool(&["pf-pencil", "--case", "A198204", "--degree", "2"], "");
    assert!(a198204.status.success());
    let stdout = String::from_utf8(a198204.stdout).expect("stdout is utf8");
    assert!(stdout.contains("common denominator cleared: 2"));
    assert!(stdout.contains("finite evidence ok: true"));

    let a198204_d0 = run_polytool(&["pf-pencil", "--case", "A198204", "--degree", "0"], "");
    assert!(a198204_d0.status.success());
    let stdout = String::from_utf8(a198204_d0.stdout).expect("stdout is utf8");
    assert!(stdout.contains("unsupported"));
}

#[test]
fn bench_interlacing_reports_tsv() {
    let output = run_polytool(
        &[
            "bench",
            "interlacing",
            "--sequence",
            "eulerian",
            "--max-n",
            "4",
            "--repeat",
            "1",
        ],
        "",
    );

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("stdout is utf8");
    assert!(
        stdout.starts_with("sequence\tleft_index\tright_index\tdegree\trepeat\tavg_us\tresult\n")
    );
    assert!(stdout.contains("eulerian\t"));
}

#[test]
fn bench_recurrence_fixture_reports_tsv() {
    let report_path = std::env::temp_dir().join(format!(
        "polytool-recurrence-bench-report-{}.md",
        std::process::id()
    ));
    let report_arg = report_path.to_string_lossy().into_owned();
    let output = run_polytool(
        &[
            "bench",
            "recurrence-fixtures",
            "--only",
            "01_scalar_geometric",
            "--repeat",
            "1",
            "--summary",
            "--report",
            &report_arg,
        ],
        "",
    );

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("stdout is utf8");
    assert!(stdout.starts_with(
        "slug\trun\tfound\telapsed_ms\tcandidates\tunknowns\tweighted\tfit_rows\tverify_rows\trecurrence\n"
    ));
    assert!(stdout.contains("01_scalar_geometric\t1\ttrue\t"));
    assert!(stdout.contains("P(n) = 2 P(n-1)"));
    assert!(stdout.contains("# fixture_summary\n"));
    assert!(stdout.contains("# category_summary\n"));

    let report = fs::read_to_string(&report_path).expect("read benchmark report");
    assert!(report.contains("# Recurrence Fixture Benchmark Report"));
    assert!(report.contains("| synthetic | 1 | 1 | 1 |"));
    assert!(report.contains("| 01_scalar_geometric |"));
    let _ = fs::remove_file(report_path);
}

#[test]
fn bench_recurrence_fixture_reports_json_and_compare() {
    let output = run_polytool(
        &[
            "bench",
            "recurrence-fixtures",
            "--only",
            "01_scalar_geometric",
            "--repeat",
            "1",
            "--format",
            "json",
        ],
        "",
    );

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("stdout is utf8");
    assert!(stdout.contains("\"schema\": \"polytool.bench.recurrence-fixtures.v1\""));
    assert!(stdout.contains("\"fixture_summaries\""));
    assert!(stdout.contains("\"diagnostics\""));
    assert!(stdout.contains("\"generated_candidates\""));

    let old_path = std::env::temp_dir().join(format!(
        "polytool-recurrence-bench-old-{}.json",
        std::process::id()
    ));
    let new_path = std::env::temp_dir().join(format!(
        "polytool-recurrence-bench-new-{}.json",
        std::process::id()
    ));
    fs::write(&old_path, &stdout).expect("write old benchmark JSON");
    fs::write(&new_path, &stdout).expect("write new benchmark JSON");
    let old_arg = old_path.to_string_lossy().into_owned();
    let new_arg = new_path.to_string_lossy().into_owned();

    let compare = run_polytool(&["bench", "compare", &old_arg, &new_arg, "--top", "1"], "");
    assert!(compare.status.success());
    let compare_stdout = String::from_utf8(compare.stdout).expect("compare stdout is utf8");
    assert!(compare_stdout.contains("# fixture_compare\n"));
    assert!(compare_stdout.contains("01_scalar_geometric"));

    let compare_json = run_polytool(
        &[
            "bench", "compare", &old_arg, &new_arg, "--top", "1", "--format", "json",
        ],
        "",
    );
    assert!(compare_json.status.success());
    let compare_json_stdout =
        String::from_utf8(compare_json.stdout).expect("compare JSON stdout is utf8");
    assert!(compare_json_stdout.contains("\"schema\": \"polytool.bench.compare.v1\""));
    assert!(compare_json_stdout.contains("\"worst_regressions\""));

    let _ = fs::remove_file(old_path);
    let _ = fs::remove_file(new_path);
}

#[test]
fn hstar_inequalities_json_reports_named_failure() {
    let output = run_polytool(
        &["hstar-inequalities", "--dimension", "3", "--json"],
        "1,20,1,0\n",
    );

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("stdout is utf8");
    let json: serde_json::Value = serde_json::from_str(&stdout).expect("valid hstar JSON");
    assert_eq!(json["items"][0]["all_applicable_hold"], false);
    let checks = json["items"][0]["checks"].as_array().expect("checks array");
    assert!(checks.iter().any(|check| {
        check["family"] == "Balletti-Higashitani"
            && check["applicable"] == true
            && check["holds"] == false
            && check["reference"].as_str().unwrap().contains("Theorem 1.4")
    }));
}

#[test]
fn coefficient_tests_json_reports_kurtz() {
    let output = run_polytool(
        &["coefficient-tests", "--json"],
        "1000,1110,111,1\n1,4,6,4,1\n",
    );

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("stdout is utf8");
    let json: serde_json::Value = serde_json::from_str(&stdout).expect("valid coefficient JSON");
    assert_eq!(json["items"][0]["kurtz"]["holds"], true);
    assert_eq!(json["items"][0]["kurtz"]["implies_real_rooted"], true);
    assert_eq!(json["items"][1]["kurtz"]["holds"], false);
}

#[test]
fn cyclic_sieving_json_checks_order_two() {
    let output = run_polytool(
        &[
            "cyclic-sieving",
            "--order",
            "2",
            "--fixed-counts",
            "2,0",
            "--json",
        ],
        "1,1\n",
    );

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("stdout is utf8");
    let json: serde_json::Value = serde_json::from_str(&stdout).expect("valid csp JSON");
    assert_eq!(json["items"][0]["report"]["holds"], true);
    assert_eq!(
        json["items"][0]["report"]["evaluations"][1]["integer_value"],
        "0"
    );
}

#[test]
fn cyclic_sieving_sequence_uses_default_offsets() {
    let output = run_polytool(
        &["cyclic-sieving-sequence", "--first-index", "2", "--json"],
        "1,1\n1,1,1\n",
    );

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("stdout is utf8");
    let json: serde_json::Value = serde_json::from_str(&stdout).expect("valid csp sequence JSON");
    let orders = json["items"][0]["candidate_orders"]
        .as_array()
        .expect("candidate orders");
    assert!(orders.iter().any(|item| item["order"] == 2));
    assert!(orders.iter().any(|item| item["order"] == 5));
}
