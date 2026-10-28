//! CLI tests for coefficients beyond `i64` in `polytool family-check` and
//! `polytool properties`: the recurrence search and the gamma expansion.

use num_bigint::BigInt;
use polytool::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};
use serde_json::Value;
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
    // The CLI may exit on an argument error before reading its input.
    if let Err(error) = child
        .stdin
        .take()
        .expect("polytool stdin")
        .write_all(input.as_bytes())
    {
        assert_eq!(error.kind(), std::io::ErrorKind::BrokenPipe);
    }
    child.wait_with_output().expect("wait for polytool")
}

fn json(output: &Output) -> Value {
    serde_json::from_slice(&output.stdout).expect("stdout is JSON")
}

/// Rows `(1 + scale x)^n` for `n = 0..=last`, as comma-separated lines and as
/// `BigInt` coefficient vectors.
fn scaled_binomial_rows(scale: i64, last: usize) -> (String, Vec<Vec<BigInt>>) {
    let mut rows: Vec<Vec<BigInt>> = vec![vec![BigInt::from(1)]];
    for _ in 0..last {
        let previous = rows.last().unwrap().clone();
        let mut next = vec![BigInt::from(0); previous.len() + 1];
        for (k, value) in previous.iter().enumerate() {
            next[k] += value;
            next[k + 1] += value * scale;
        }
        rows.push(next);
    }
    let text = rows
        .iter()
        .map(|row| {
            row.iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(",")
        })
        .collect::<Vec<_>>()
        .join("\n");
    (text + "\n", rows)
}

#[test]
fn family_check_recurrence_accepts_coefficients_beyond_i64() {
    let (text, rows) = scaled_binomial_rows(10_000_000_000, 6);
    assert!(rows[6].iter().any(|value| *value > BigInt::from(i64::MAX)));

    let output = run_polytool(&["family-check", "--recurrence", "--json"], &text);
    assert!(output.status.success(), "{output:?}");
    let value = json(&output);
    assert_eq!(value["all_required_checks_passed"], true);
    assert_eq!(value["item_count"], 7);
    let recurrence = &value["recurrence"];
    assert_eq!(recurrence["found"], true, "{recurrence}");
    assert_eq!(recurrence["error"], Value::Null);
    let found = recurrence["recurrence"].as_str().expect("recurrence text");
    assert!(!found.is_empty());

    let text_output = run_polytool(&["family-check", "--recurrence"], &text);
    assert!(text_output.status.success(), "{text_output:?}");
    let stdout = String::from_utf8(text_output.stdout).unwrap();
    assert!(stdout.contains("- Found: "), "{stdout}");
    assert!(!stdout.contains("fit in i64"), "{stdout}");
}

#[test]
fn family_check_recurrence_is_unchanged_for_rows_that_fit_in_i64() {
    let (text, rows) = scaled_binomial_rows(1, 7);
    let rows_i64: Vec<Vec<i64>> = rows
        .iter()
        .map(|row| {
            row.iter()
                .map(|value| i64::try_from(value).unwrap())
                .collect()
        })
        .collect();
    let expected = find_recurrence_adaptive(&rows_i64, &AdaptiveSearchOptions::default())
        .expect("Pascal rows satisfy a recurrence");

    let output = run_polytool(&["family-check", "--recurrence", "--json"], &text);
    assert!(output.status.success(), "{output:?}");
    let value = json(&output);
    let recurrence = &value["recurrence"];
    assert_eq!(recurrence["found"], true);
    assert_eq!(recurrence["recurrence"], expected.recurrence.to_string());
    assert_eq!(recurrence["unknowns"], expected.num_unknowns);
    assert_eq!(recurrence["equations"], expected.num_equations);
    assert_eq!(recurrence["candidates_tried"], expected.candidates_tried);
}

#[test]
fn family_check_recurrence_reports_a_missing_recurrence_for_large_rows() {
    // Two rows: fewer than the three needed, with entries beyond i64.
    let (text, _) = scaled_binomial_rows(10_000_000_000, 2);
    let two_rows: String = text.lines().take(2).collect::<Vec<_>>().join("\n") + "\n";
    let output = run_polytool(&["family-check", "--recurrence", "--json"], &two_rows);
    let value = json(&output);
    assert_eq!(value["recurrence"]["found"], false);
    assert_eq!(value["recurrence"]["error"], "need at least 3 polynomials");
}

/// Gamma coefficients `(-1)^i d/(d-i) C(d-i, i)` of `1 + t^d`.
fn gamma_of_one_plus_power(d: usize) -> Vec<String> {
    let mut binomial = vec![vec![0u128; d + 1]; d + 1];
    for n in 0..=d {
        binomial[n][0] = 1;
        for k in 1..=n {
            binomial[n][k] = binomial[n - 1][k - 1] + binomial[n - 1][k];
        }
    }
    (0..=d / 2)
        .map(|i| {
            let magnitude = if i == 0 {
                1
            } else {
                binomial[d - i][i] + binomial[d - i - 1][i - 1]
            };
            if i % 2 == 0 {
                magnitude.to_string()
            } else {
                format!("-{magnitude}")
            }
        })
        .collect()
}

fn one_plus_power_row(d: usize) -> String {
    let mut coefficients = vec!["0"; d + 1];
    coefficients[0] = "1";
    coefficients[d] = "1";
    coefficients.join(",") + "\n"
}

#[test]
fn gamma_coefficients_are_exact_when_they_overflow_i64_but_the_row_fits() {
    // 1 + t^100 has coefficients 0 and 1, but its gamma coefficients reach 10^20.
    let expected = gamma_of_one_plus_power(100);
    assert!(expected
        .iter()
        .any(|value| value.trim_start_matches('-').parse::<u128>().unwrap() > i64::MAX as u128));
    let input = one_plus_power_row(100);

    let properties = json(&run_polytool(&["properties", "--json"], &input));
    let item = &properties["items"][0];
    assert_eq!(item["palindromic"], true);
    assert_eq!(item["gamma_positive"], false);
    assert_eq!(
        item["gamma_coefficients"],
        serde_json::to_value(&expected).unwrap()
    );

    let family = json(&run_polytool(&["family-check", "--json"], &input));
    assert_eq!(
        family["items"][0]["gamma_coefficients"],
        serde_json::to_value(&expected).unwrap()
    );
}

#[test]
fn gamma_coefficients_are_unchanged_for_small_rows() {
    let properties = json(&run_polytool(&["properties", "--json"], "1,4,1\n1,2,0\n"));
    assert_eq!(
        properties["items"][0]["gamma_coefficients"],
        serde_json::json!(["1", "2"])
    );
    assert_eq!(properties["items"][0]["gamma_positive"], true);
    // Not palindromic: no gamma expansion.
    assert_eq!(properties["items"][1]["gamma_coefficients"], Value::Null);
    assert_eq!(properties["items"][1]["gamma_positive"], false);
}
