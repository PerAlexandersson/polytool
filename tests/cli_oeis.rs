use std::process::{Command, Output};

fn run_polytool(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_polytool"))
        .args(args)
        .output()
        .expect("spawn polytool")
}

#[test]
fn oeis_list_contains_verified_catalog() {
    let output = run_polytool(&["oeis", "list", "--json"]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["schema"], "polytool.oeis-catalog.v1");
    assert_eq!(value["count"], 755);
    assert!(value["sequences"]
        .as_array()
        .unwrap()
        .iter()
        .any(|entry| entry["id"] == "A008292"));
    assert!(value["sequences"]
        .as_array()
        .unwrap()
        .iter()
        .any(|entry| entry["id"] == "A390883" && entry["status"] == "validated"));

    let all = run_polytool(&["oeis", "list", "--json", "--include-experimental"]);
    assert!(all.status.success());
    let all_value: serde_json::Value = serde_json::from_slice(&all.stdout).unwrap();
    assert_eq!(all_value["count"], 785);
}

#[test]
fn oeis_triangle_and_polynomial_formats_preserve_indices() {
    let triangle = run_polytool(&[
        "oeis", "generate", "A008292", "--rows", "4", "--format", "triangle",
    ]);
    assert!(triangle.status.success());
    assert_eq!(
        String::from_utf8(triangle.stdout).unwrap(),
        "1: 1\n2: 1 1\n3: 1 4 1\n4: 1 11 11 1\n"
    );

    let polynomial = run_polytool(&[
        "oeis",
        "generate",
        "A008292",
        "--start-row",
        "3",
        "--rows",
        "2",
    ]);
    assert!(polynomial.status.success());
    assert_eq!(
        String::from_utf8(polynomial.stdout).unwrap(),
        "P_3(t) = 1 + 4t + t^2\nP_4(t) = 1 + 11t + 11t^2 + t^3\n"
    );
}

#[test]
fn oeis_bfile_uses_flattened_indices_and_whole_rows() {
    let output = run_polytool(&[
        "oeis",
        "generate",
        "A008292",
        "--rows",
        "5",
        "--format",
        "bfile",
        "--max-terms",
        "7",
    ]);
    assert!(output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        "1 1\n2 1\n3 1\n4 1\n5 4\n6 1\n"
    );
}

#[test]
fn oeis_bfile_preserves_irregular_row_widths() {
    let output = run_polytool(&[
        "oeis", "generate", "A390883", "--rows", "6", "--format", "bfile",
    ]);
    assert!(output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        "0 1\n1 1\n2 1\n3 2\n4 1\n5 10\n6 1\n7 37\n8 10\n9 1\n10 126\n11 105\n"
    );
}

#[test]
fn oeis_bfile_rejects_unverified_prefix_mapping() {
    let output = run_polytool(&[
        "oeis", "generate", "A095704", "--rows", "3", "--format", "bfile",
    ]);
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8(output.stderr)
        .unwrap()
        .contains("prefix mapping is not verified"));
}

#[test]
fn oeis_info_exports_recurrence_dialects() {
    let output = run_polytool(&["oeis", "info", "A008292", "--json"]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["id"], "A008292");
    assert!(value["recurrence"].as_str().unwrap().contains("P(n-1)"));
    assert!(value["latex"].as_str().unwrap().contains("P(n-1)"));
    assert!(value["mathematica"].as_str().unwrap().contains("P["));
    assert!(value["sage"].as_str().unwrap().contains("def P"));
    assert!(value["python"].as_str().unwrap().contains("def P"));
}

#[test]
fn experimental_oeis_entries_require_an_explicit_flag() {
    let rejected = run_polytool(&["oeis", "generate", "A035469", "--rows", "2"]);
    assert_eq!(rejected.status.code(), Some(2));
    assert!(String::from_utf8(rejected.stderr)
        .unwrap()
        .contains("is experimental"));

    let accepted = run_polytool(&[
        "oeis",
        "generate",
        "A035469",
        "--rows",
        "2",
        "--include-experimental",
    ]);
    assert!(accepted.status.success());
}
