//! Shared CLI seed selection, normalization, Lean output and refusal reasons.

use crate::common::{run_polytool as run, stderr, stdout};
use polytool::recurrence::RecurrenceJson;
use serde_json::Value;

#[test]
fn lean_json_fields_require_an_explicit_name_in_records_and_menus() {
    let input = "1\n2\n4\n8\n16\n32\n64\n128";
    for mode in [None, Some("--representations"), Some("--theorem-matches")] {
        for name in [None, Some("Rows")] {
            let mut args = vec![
                "recurrence",
                "--json",
                "--no-laurent",
                "--max-recurrence-length",
                "1",
                "--max-variable-degree",
                "0",
                "--max-index-degree",
                "0",
                "--max-derivative-order",
                "0",
            ];
            if let Some(mode) = mode {
                args.push(mode);
            }
            if let Some(name) = name {
                args.extend(["--lean-name", name]);
            }
            let output = run(&args, input);
            assert!(output.status.success(), "{}", stderr(&output));
            let report: Value = serde_json::from_slice(&output.stdout).unwrap();
            let record = match mode {
                Some("--representations") => {
                    assert!(!report["representations"].as_array().unwrap().is_empty());
                    &report["representations"][0]["recurrence_json"]
                }
                Some("--theorem-matches") => &report["recurrence_json"],
                _ => &report,
            };
            assert!(record["presentation"].is_object());
            assert!(record["normalized_recurrence"].is_object());
            assert_eq!(record.get("lean").is_some(), name.is_some());
            assert!(record.get("lean_reason").is_none());
            if let Some(name) = name {
                let name = if mode == Some("--representations") {
                    format!("{name}_0")
                } else {
                    name.to_string()
                };
                assert!(record["lean"]
                    .as_str()
                    .unwrap()
                    .contains(&format!("def {name} :")));
            }
        }
    }
    let output = run(&["recurrence", "--format", "lean"], input);
    assert!(output.status.success(), "{}", stderr(&output));
    assert!(stdout(&output).contains("def P : ℕ → ℝ[X]"));
}

fn catalog_rows(id: &str, count: usize) -> String {
    polytool::oeis::by_id(id)
        .unwrap()
        .generate_rows(count)
        .unwrap()
        .iter()
        .map(|row| {
            row.iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(",")
        })
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn earliest_start_is_two_for_weighted_eulerian_rows_in_all_formats() {
    let rows = catalog_rows("A065826", 12);
    let output = run(
        &[
            "recurrence",
            "--first-index",
            "0",
            "--json",
            "--lean-name",
            "A065826",
        ],
        &rows,
    );
    assert!(output.status.success(), "{}", stderr(&output));
    let value: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["initial_polynomials"].as_array().unwrap().len(), 2);
    assert_eq!(value["presentation"]["start_index"], 2);
    assert_eq!(value["presentation"]["earliest_valid_index"], 2);
    assert!(value["latex"].as_str().unwrap().contains("(n \\geq 2)"));
    let text = run(&["recurrence", "--first-index", "0"], &rows);
    assert!(stdout(&text).contains("Range: n >= 2"));
    assert!(stdout(&text).contains("P(1) = [1, 2]"));
    let latex = run(&["recurrence", "--first-index", "0", "--latex"], &rows);
    assert_eq!(stdout(&latex).trim(), value["latex"].as_str().unwrap());
    let lean = run(
        &[
            "recurrence",
            "--first-index",
            "0",
            "--format",
            "lean",
            "--lean-name",
            "A065826",
        ],
        &rows,
    );
    assert!(lean.status.success(), "{}", stderr(&lean));
    assert_eq!(stdout(&lean).trim(), value["lean"].as_str().unwrap());
    assert!(stdout(&lean).contains("| n + 2 =>"));
}

#[test]
fn singularity_seeds_shifted_indices_normalization_and_lean_examples() {
    let rows = catalog_rows("A216118", 10);
    for first in ["0", "7"] {
        let output = run(
            &[
                "recurrence",
                "--first-index",
                first,
                "--json",
                "--lean-name",
                "A216118",
            ],
            &rows,
        );
        assert!(output.status.success(), "{}", stderr(&output));
        let value: Value = serde_json::from_slice(&output.stdout).unwrap();
        let first: usize = first.parse().unwrap();
        assert_eq!(
            value["presentation"]["left_hand_zeros"],
            serde_json::json!([(first + 1).to_string()])
        );
        assert_eq!(value["presentation"]["start_index"], first + 2);
        assert_eq!(value["initial_polynomials"].as_array().unwrap().len(), 2);
        let record: RecurrenceJson = serde_json::from_value(value.clone()).unwrap();
        let normalized = value["normalized_recurrence"].clone();
        assert!(
            normalized["denominator"]["coeffs"][1][0]
                .as_str()
                .unwrap()
                .parse::<i64>()
                .unwrap()
                > 0
        );
        assert!(value["lean"]
            .as_str()
            .unwrap()
            .contains("def A216118 : ℕ → ℝ[X]"));
        assert_eq!(
            value["lean"].as_str().unwrap().matches("example :").count(),
            4
        );
        assert_eq!(record.initial_polynomials.len(), 2);
    }
}

#[test]
fn invalid_names_and_formats_print_no_definition() {
    for args in [
        vec!["recurrence", "--format", "lean", "--lean-name", "n"],
        vec!["recurrence", "--format", "lean", "--lean-name", "Bad Name"],
        vec!["recurrence", "--format", "unknown"],
    ] {
        let output = run(&args, "1\n2\n4\n8\n16\n32");
        assert_eq!(output.status.code(), Some(2));
        assert!(output.stdout.is_empty());
    }
}

fn generated_input(record: Value, count: usize) -> String {
    let record: RecurrenceJson = serde_json::from_value(record).unwrap();
    let (rule, first, seeds) = record.to_recurrence_parts().unwrap();
    rule.generate_rows_rational(&seeds, first, count)
        .unwrap()
        .iter()
        .map(|row| {
            row.iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(",")
        })
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn future_singularity_has_an_explicit_reason_and_no_defining_range() {
    let input = generated_input(
        serde_json::json!({
            "schema":"polytool.recurrence.v1", "first_index":0, "initial_polynomials":[["1"]],
            "recurrence":{"denominator":{"coeffs":[["-20"],["1"]]},
                "terms":[{"offset":1,"deriv_order":0,"sign":"none", "coeff":{"coeffs":[["1"],["1"]]}}]}
        }),
        12,
    );
    let args = [
        "recurrence",
        "--first-index",
        "0",
        "--no-laurent",
        "--try-denominator",
        "--max-recurrence-length",
        "1",
        "--max-variable-degree",
        "0",
        "--max-index-degree",
        "1",
        "--max-derivative-order",
        "0",
        "--max-denominator-variable-degree",
        "0",
        "--max-denominator-index-degree",
        "1",
    ];
    let mut json_args = args.to_vec();
    json_args.push("--json");
    let default_output = run(&json_args, &input);
    assert!(
        default_output.status.success(),
        "{}",
        stderr(&default_output)
    );
    let default_record: Value = serde_json::from_slice(&default_output.stdout).unwrap();
    assert!(default_record.get("lean").is_none());
    assert!(default_record.get("lean_reason").is_none());
    json_args.extend(["--lean-name", "P"]);
    let output = run(&json_args, &input);
    assert!(output.status.success(), "{}", stderr(&output));
    let value: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(
        value["presentation"]["left_hand_zeros"],
        serde_json::json!(["20"])
    );
    assert_eq!(value["presentation"]["status"], "insufficient_initial_rows");
    assert!(value["presentation"]["start_index"].is_null());
    assert_eq!(value["initial_polynomials"].as_array().unwrap().len(), 12);
    assert!(value.get("lean").is_none());
    assert!(value["lean_reason"]
        .as_str()
        .unwrap()
        .contains("through index 20"));
    assert!(!value["latex"].as_str().unwrap().contains("\\geq"));
    let text = run(&args, &input);
    assert!(stdout(&text).contains("No defining range:"));
    assert!(!stdout(&text).contains("Range: n >="));
    let mut lean_args = args.to_vec();
    lean_args.extend(["--format", "lean"]);
    let lean = run(&lean_args, &input);
    assert_eq!(lean.status.code(), Some(2));
    assert!(lean.stdout.is_empty());
}

#[test]
fn t_dependent_leader_refuses_lean_but_preserves_the_json_report() {
    let input = generated_input(
        serde_json::json!({
            "schema":"polytool.recurrence.v1", "first_index":0, "initial_polynomials":[["0","1"]],
            "recurrence":{"denominator":{"coeffs":[["-1","1"],["1"]]},
                "terms":[{"offset":1,"deriv_order":0,"sign":"none", "coeff":{"coeffs":[["0","1"],["1"]]}}]}
        }),
        12,
    );
    let args = [
        "recurrence",
        "--first-index",
        "0",
        "--no-laurent",
        "--try-denominator",
        "--max-recurrence-length",
        "1",
        "--max-variable-degree",
        "1",
        "--max-index-degree",
        "1",
        "--max-derivative-order",
        "0",
        "--max-denominator-variable-degree",
        "1",
        "--max-denominator-index-degree",
        "1",
    ];
    let mut json_args = args.to_vec();
    json_args.push("--json");
    let default_output = run(&json_args, &input);
    assert!(
        default_output.status.success(),
        "{}",
        stderr(&default_output)
    );
    let default_record: Value = serde_json::from_slice(&default_output.stdout).unwrap();
    assert!(default_record.get("lean").is_none());
    assert!(default_record.get("lean_reason").is_none());
    json_args.extend(["--lean-name", "P"]);
    let output = run(&json_args, &input);
    assert!(output.status.success(), "{}", stderr(&output));
    let value: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(
        value["presentation"]["left_hand_zeros"],
        serde_json::json!([])
    );
    assert!(value.get("lean").is_none());
    assert!(value["lean_reason"]
        .as_str()
        .unwrap()
        .contains("depends on t"));
    let mut lean_args = args.to_vec();
    lean_args.extend(["--format", "lean"]);
    let lean = run(&lean_args, &input);
    assert_eq!(lean.status.code(), Some(2));
    assert!(lean.stdout.is_empty());
    assert!(stderr(&lean).contains("depends on t"));
}

#[test]
fn fitted_text_scaling_scalar_latex_and_nonempty_unquoted_zero_lines() {
    let input = "1\n1/2\n1/4\n1/8\n1/16\n1/32\n1/64\n1/128\n";
    let text = run(&["recurrence", "--no-laurent"], input);
    assert!(text.status.success(), "{}", stderr(&text));
    assert!(stdout(&text).starts_with("P(n) = 1/2 P(n-1)\n"));
    assert!(!stdout(&text).contains("Left-hand zeros"));
    let latex = run(&["recurrence", "--no-laurent", "--latex"], input);
    assert!(
        stdout(&latex).contains("$2 P_{n} = P_{n-1}"),
        "{}",
        stdout(&latex)
    );
    assert!(!stdout(&latex).contains("(2) P_{n}"));
    assert!(!stdout(&latex).contains("Left-hand zeros"));
    let singular = run(
        &["recurrence", "--first-index", "0"],
        &catalog_rows("A216118", 10),
    );
    assert!(stdout(&singular).contains("Left-hand zeros (nonnegative integers): [1]"));
    assert!(!stdout(&singular).contains("[\"1\"]"));
}

#[test]
fn python_and_oeis_programs_use_the_checked_seed_prefix() {
    let input = "1\n-1 1\n-1 0 1\n-1 -1 1 1\n-1 -2 0 2 1\n-1 -3 -2 2 3 1\n-1 -4 -5 0 5 4 1\n-1 -5 -9 -5 5 9 5 1\n";
    let searched = run(&["recurrence", "--json"], input);
    assert!(searched.status.success(), "{}", stderr(&searched));
    let value: Value = serde_json::from_slice(&searched.stdout).unwrap();
    let count = value["initial_polynomials"].as_array().unwrap().len();
    assert_eq!(count, 2);
    let python = run(&["recurrence", "--python"], input);
    let python = stdout(&python);
    assert!(python.contains("    1: [1],"), "{python}");
    assert!(python.contains("    2: [-1, 1],"), "{python}");
    assert!(!python.contains("    3: ["));
    let programs = run(
        &[
            "recurrence-generate",
            "--recurrence",
            "-",
            "--format",
            "oeis-json",
            "--count",
            "4",
        ],
        &value.to_string(),
    );
    assert!(programs.status.success(), "{}", stderr(&programs));
    let programs: Value = serde_json::from_slice(&programs.stdout).unwrap();
    let mathematica = programs["mathematica"].as_str().unwrap();
    assert!(mathematica.contains("P[1, x_] := 1;"), "{mathematica}");
    assert!(mathematica.contains("P[2, x_] := -1 + x;"), "{mathematica}");
    assert!(!mathematica.contains("P[3, x_] :="));
    let pari = programs["pari"].as_str().unwrap();
    assert!(pari.contains("v[1] = 1; v[2] = -1 + x;"), "{pari}");
    assert!(!pari.contains("v[3] ="));
    assert!(pari.contains("for(n = 3, N"), "{pari}");
}

#[test]
fn lean_keywords_and_collisions_are_rejected_with_a_reason_by_the_cli() {
    for (name, reason) in [
        ("have", "keyword"),
        ("Type", "keyword"),
        ("at", "keyword"),
        ("derivative", "clashes"),
        ("Polynomial", "clashes"),
        ("Real", "clashes"),
    ] {
        let output = run(
            &["recurrence", "--format", "lean", "--lean-name", name],
            "1\n2\n4\n8\n16",
        );
        assert_eq!(output.status.code(), Some(2));
        assert!(output.stdout.is_empty());
        assert!(stderr(&output).contains(reason), "{}", stderr(&output));
    }
}

#[test]
fn lean_menus_keep_renderable_candidates_and_explain_the_others() {
    let input = generated_input(
        serde_json::json!({
            "schema":"polytool.recurrence.v1","first_index":0,"initial_polynomials":[["0","1"]],
            "recurrence":{"denominator":{"coeffs":[["-1","1"],["1"]]},"terms":[
                {"offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["0","1"],["1"]]}}]}
        }),
        12,
    );
    let args = [
        "recurrence",
        "--representations",
        "--format",
        "lean",
        "--lean-name",
        "Menu",
        "--no-laurent",
        "--try-denominator",
        "--max-recurrence-length",
        "2",
        "--max-variable-degree",
        "1",
        "--max-index-degree",
        "1",
        "--max-derivative-order",
        "0",
        "--max-denominator-variable-degree",
        "1",
        "--max-denominator-index-degree",
        "1",
    ];
    let mixed = run(&args, &input);
    assert!(
        mixed.status.success(),
        "{}\n{}",
        stdout(&mixed),
        stderr(&mixed)
    );
    assert!(stdout(&mixed).contains("def Menu_"), "{}", stdout(&mixed));
    assert!(
        stdout(&mixed).contains("-- Representation"),
        "{}",
        stdout(&mixed)
    );
    assert!(
        stdout(&mixed).contains("depends on t"),
        "{}",
        stdout(&mixed)
    );
    let mut only_unsupported = args;
    let lag = only_unsupported
        .iter()
        .position(|&s| s == "--max-recurrence-length")
        .unwrap();
    only_unsupported[lag + 1] = "1";
    let none = run(&only_unsupported, &input);
    assert_eq!(none.status.code(), Some(2));
    assert!(!stdout(&none).contains("def Menu_"));
    assert!(stdout(&none).contains("depends on t"));
    let partial = run(
        &[
            "recurrence",
            "--representations",
            "--format",
            "lean",
            "--max-candidates",
            "1",
        ],
        "1\n2\n4\n8\n16\n32\n64\n128",
    );
    assert!(partial.status.success(), "{}", stderr(&partial));
    assert!(stdout(&partial).contains("def P_0"));
    assert!(stdout(&partial).contains("menu is incomplete"));
}

#[test]
fn default_output_changes_are_documented_separately_from_additions() {
    let changelog = include_str!("../../CHANGELOG.md");
    let changed = changelog
        .split("### Changed\n")
        .nth(1)
        .unwrap()
        .split("### Added\n")
        .next()
        .unwrap();
    for subject in [
        "seed count",
        "stated range",
        "extra left-hand zero line",
        "Recurrence LaTeX",
        "MCP search/family `latex`",
    ] {
        assert!(
            changed.contains(subject),
            "missing default-output change: {subject}"
        );
    }
}

#[test]
fn lean_examples_cover_alternation_forcing_rational_constants_and_order_three() {
    let cases = [
        (
            "AlternatingRows",
            serde_json::json!({"schema":"polytool.recurrence.v1","first_index":0,
            "initial_polynomials":[["1"]],"recurrence":{"terms":[{"offset":1,"deriv_order":0,
                "sign":"alternating_n","coeff":{"coeffs":[["2","1"]]}}]}}),
            vec![
                "--max-recurrence-length",
                "1",
                "--max-variable-degree",
                "1",
                "--try-alternating-sign",
            ],
            "(-1 : ℝ) ^ (n + 1)",
        ),
        (
            "InhomogeneousRows",
            serde_json::json!({"schema":"polytool.recurrence.v1","first_index":0,
            "initial_polynomials":[["1"]],"recurrence":{"terms":[{"offset":1,"deriv_order":0,
                "sign":"none","coeff":{"coeffs":[["0","1"]]}}],"inhomogeneous":{"coeffs":[["1"]]}}}),
            vec![
                "--max-recurrence-length",
                "1",
                "--max-variable-degree",
                "1",
                "--try-inhomogeneous",
                "--max-inhomogeneous-variable-degree",
                "0",
                "--max-inhomogeneous-index-degree",
                "0",
            ],
            "+ 1",
        ),
        (
            "RationalRows",
            serde_json::json!({"schema":"polytool.recurrence.v1","first_index":0,
            "initial_polynomials":[["1"]],"recurrence":{"terms":[{"offset":1,"deriv_order":0,
                "sign":"none","coeff":{"coeffs":[["1/2"]]}}]}}),
            vec!["--max-recurrence-length", "1", "--max-variable-degree", "0"],
            "C (1 / 2)",
        ),
        (
            "OrderThreeRows",
            serde_json::json!({"schema":"polytool.recurrence.v1","first_index":0,
            "initial_polynomials":[["1"],["2"],["4"]],"recurrence":{"terms":[
                {"offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["1"]]}},
                {"offset":2,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["1"]]}},
                {"offset":3,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["1"]]}}]}}),
            vec![
                "--min-recurrence-length",
                "3",
                "--max-recurrence-length",
                "3",
                "--max-variable-degree",
                "0",
            ],
            "| n + 3 =>",
        ),
    ];
    for (name, record, options, feature) in cases {
        let input = generated_input(record, 12);
        let mut args = vec![
            "recurrence",
            "--first-index",
            "0",
            "--no-laurent",
            "--format",
            "lean",
            "--lean-name",
            name,
            "--max-index-degree",
            "0",
            "--max-derivative-order",
            "0",
        ];
        args.extend(options);
        let output = run(&args, &input);
        assert!(output.status.success(), "{name}: {}", stderr(&output));
        let text = stdout(&output);
        assert!(text.contains(feature), "{name}: {text}");
        assert!(text.lines().all(|line| line.chars().count() <= 100));
    }
}
