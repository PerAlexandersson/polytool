//! Root-only JSON, held-out evidence, and index-preserving recurrence exports.

use polytool::recurrence::RecurrenceJson;
use serde_json::{json, Value};
use std::process::Output;

use crate::common::run_polytool as run;

fn output_json(output: Output) -> Value {
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).unwrap()
}

fn source() -> RecurrenceJson {
    serde_json::from_value(json!({
        "schema":"polytool.recurrence.v1", "first_index":0,
        "initial_polynomials":[["1"]],
        "recurrence":{"terms":[{"offset":1,"deriv_order":0,"sign":"none",
            "coeff":{"coeffs":[["0","1"],["1"]]}}]}
    }))
    .unwrap()
}

fn source_rows() -> String {
    let (rule, first, seeds) = source().to_recurrence_parts().unwrap();
    rule.generate_rows_rational(&seeds, first, 12)
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

fn search(extra: &[&str]) -> Value {
    let mut args = vec![
        "recurrence",
        "--json",
        "--no-laurent",
        "--max-rec-len",
        "1",
        "--max-diff-deg",
        "0",
        "--max-idx-deg",
        "1",
        "--max-var-deg",
        "1",
    ];
    args.extend(extra);
    output_json(run(&args, &source_rows()))
}

#[test]
fn root_json_is_minimal_exact_and_handles_zero() {
    let result = output_json(run(&["real-rooted", "--json"],
        "0,0\n7\n1,2,1,0\n1,0,1\n1000000000000000000000000000000,2000000000000000000000000000000,1000000000000000000000000000000\n"));
    let items = result["items"].as_array().unwrap();
    assert_eq!(items.len(), 5);
    assert_eq!(items[0]["degree"], Value::Null);
    for (i, (degree, real)) in [(0, true), (2, true), (2, false), (2, true)]
        .iter()
        .enumerate()
    {
        assert_eq!(
            items[i + 1],
            json!({"index":i+1,"degree":degree,"real_rooted":real})
        );
    }
    assert_eq!(
        output_json(run(&["real-rooted", "--json"], "")),
        json!({"items":[]})
    );
    assert_eq!(run(&["real-rooted", "--bad"], "1").status.code(), Some(2));
    assert_eq!(
        run(&["real-rooted", "--json"], "not a polynomial!")
            .status
            .code(),
        Some(2)
    );
    assert_eq!(
        String::from_utf8(run(&["real-rooted"], "1,2,1").stdout).unwrap(),
        "1 + 2t + t^2: real-rooted\n"
    );
}

#[test]
fn search_index_and_skip_replay_the_same_rows() {
    for (index, skip) in [("0", "0"), ("5", "0"), ("0", "2"), ("5", "2")] {
        let value = search(&["--first-index", index, "--skip-prefix", skip]);
        let skipped: usize = skip.parse().unwrap();
        let first = index.parse::<usize>().unwrap() + skipped;
        assert_eq!(value["first_index"], first);
        let generated = output_json(run(
            &[
                "recurrence-generate",
                "--recurrence",
                "-",
                "--json",
                "--rows",
                &(12 - skipped).to_string(),
            ],
            &value.to_string(),
        ));
        let expected = output_json(run(
            &[
                "recurrence-generate",
                "--recurrence",
                "-",
                "--json",
                "--rows",
                "12",
            ],
            &serde_json::to_string(&source()).unwrap(),
        ));
        assert_eq!(
            generated["polynomials"],
            json!(&expected["polynomials"].as_array().unwrap()[skipped..])
        );
    }
    assert_eq!(search(&[])["first_index"], 1);
    assert_eq!(search(&["--skip-prefix", "2"])["first_index"], 1);
    assert_eq!(
        run(&["recurrence", "--first-index", "-1"], "1\n1")
            .status
            .code(),
        Some(2)
    );
}

#[test]
fn heldout_count_uses_actual_polynomial_degrees() {
    let value = search(&["--first-index", "0"]);
    let fit = value["search"]["fit_polynomials"].as_u64().unwrap();
    let expected: u64 = (fit..12).map(|n| n + 1).sum();
    assert!(expected > 0);
    assert_eq!(value["search"]["verification_equations"], expected);
    let no_verify = search(&["--no-verify"]);
    assert_eq!(no_verify["search"]["verification_equations"], 0);
    assert_eq!(no_verify["search"]["verification_polynomials"], 0);
}

#[test]
fn representation_indices_and_theorem_record_agree() {
    let menu = search(&["--first-index", "5", "--representations"]);
    for item in menu["representations"].as_array().unwrap() {
        assert_eq!(item["recurrence_json"]["first_index"], 5);
        assert_eq!(item["last_checked_index"], 16);
        assert!(item["first_checked_index"].as_u64().unwrap() >= 6);
        assert!(item["verification_equations"].as_u64().unwrap() > 0);
    }
    assert!(!menu["representations"].as_array().unwrap().is_empty());
    let zero_menu = search(&["--first-index", "0", "--representations"]);
    assert_eq!(
        zero_menu["representations"][0]["complexity"]["nonzero_coefficients"],
        3
    );
    assert_eq!(
        menu["representations"][0]["complexity"]["nonzero_coefficients"],
        4
    );
    let theorem = search(&["--first-index", "0", "--theorem-matches"]);
    assert_eq!(theorem["recurrence_json"]["first_index"], 0);
    assert_eq!(theorem["finite_evidence"]["first_index"], 0);
}

#[test]
fn generation_reindex_transports_all_terms_and_parity() {
    let record: RecurrenceJson = serde_json::from_value(json!({
        "schema":"polytool.recurrence.v1", "first_index":2,
        "initial_polynomials":[["1","1"]],
        "recurrence":{
            "denominator":{"coeffs":[["1"],["1"]]},
            "inhomogeneous":{"coeffs":[["0"],["1"]]},
            "terms":[{"offset":1,"deriv_order":1,"sign":"alternating_n",
                "coeff":{"coeffs":[["0","1"],["1"]]}}]
        }
    }))
    .unwrap();
    let input = serde_json::to_string(&record).unwrap();
    let baseline = output_json(run(
        &[
            "recurrence-generate",
            "--recurrence",
            "-",
            "--json",
            "--rows",
            "8",
        ],
        &input,
    ));
    for first in ["0", "9"] {
        let generated = output_json(run(
            &[
                "recurrence-generate",
                "--recurrence",
                "-",
                "--json",
                "--rows",
                "8",
                "--first-index",
                first,
            ],
            &input,
        ));
        assert_eq!(generated["first_index"], first.parse::<usize>().unwrap());
        assert_eq!(generated["polynomials"], baseline["polynomials"]);
    }
    let mut old = serde_json::to_value(source()).unwrap();
    // Old search metadata remains accepted without the new count.
    old["search"] = search(&[])["search"].clone();
    old["search"]
        .as_object_mut()
        .unwrap()
        .remove("verification_equations");
    assert!(serde_json::from_value::<RecurrenceJson>(old).is_ok());
}

#[test]
fn latex_cli_and_json_share_the_indexed_definition() {
    let value = search(&["--first-index", "5", "--skip-prefix", "2"]);
    let latex = value["latex"].as_str().unwrap();
    assert!(latex.starts_with("% Candidate fitted to finite rows;"));
    assert!(latex.contains("$P_{7} = "));
    let seed_count = value["initial_polynomials"].as_array().unwrap().len();
    assert!(latex.ends_with(&format!("(n \\geq {})$", 7 + seed_count)));
    assert!(!latex.contains("P(n"));

    for option in [vec!["--latex"], vec!["--format", "latex"]] {
        let mut args = vec![
            "recurrence",
            "--first-index",
            "5",
            "--skip-prefix",
            "2",
            "--no-laurent",
            "--max-rec-len",
            "1",
            "--max-diff-deg",
            "0",
            "--max-idx-deg",
            "1",
            "--max-var-deg",
            "1",
        ];
        args.extend(option);
        let output = run(&args, &source_rows());
        assert!(output.status.success());
        assert_eq!(String::from_utf8(output.stdout).unwrap().trim_end(), latex);
    }
    let report = search(&["--theorem-matches", "--first-index", "0"]);
    assert!(report["recurrence_json"]["latex"]
        .as_str()
        .unwrap()
        .contains("$P_{0} = 1$"));
    let report = search(&["--representations", "--first-index", "0"]);
    for item in report["representations"].as_array().unwrap() {
        assert!(item["latex"].as_str().unwrap().contains("P_{n}"));
        assert!(item["recurrence_json"]["latex"]
            .as_str()
            .unwrap()
            .contains("$P_{0} = 1$"));
    }
}

#[test]
fn latex_factors_index_dependent_coefficients() {
    // Rows of (n + 2) P_n = (2n + 1)(1 + t) P_{n-1} - (n - 1)(1 - t)^2 P_{n-2}.
    let rows = "1\n1,1\n1,3,1\n1,6,6,1\n1,10,20,10,1\n1,15,50,50,15,1\n\
                1,21,105,175,105,21,1\n1,28,196,490,490,196,28,1\n\
                1,36,336,1176,1764,1176,336,36,1\n1,45,540,2520,5292,5292,2520,540,45,1\n\
                1,55,825,4950,13860,19404,13860,4950,825,55,1\n\
                1,66,1210,9075,32670,60984,60984,32670,9075,1210,66,1";
    let output = run(&["recurrence", "--latex", "--first-index", "0"], rows);
    assert!(output.status.success());
    let text = String::from_utf8(output.stdout).unwrap();
    assert!(text.contains(
        "$(n + 2) P_{n} = (2n + 1)(1 + t) P_{n-1} - (n - 1)(1 - t)^{2} P_{n-2} \\quad (n \\geq 2)$"
    ));
    assert!(text.starts_with("% Candidate fitted to finite rows;"));
}
