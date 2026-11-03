//! `polytool recurrence-generate` OEIS text formats (`oeis-example`,
//! `oeis-formula`, `mathematica`, `pari`, `oeis-json`) and the b-file with
//! `--first-row` and `--regular-triangle`.

use std::path::PathBuf;
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

/// `P(n) = (1 - 2x + 2nx) P(n-1) + (2x - 2x^2) P'(n-1)`, `P(0) = 1`; A156919
/// is its triangle reversed from the row `n = 1` on, with offset 0.
const DERIVATIVE: &str = r#"{"schema": "polytool.recurrence.v1", "first_index": 0,
  "initial_polynomials": [["1"]],
  "recurrence": {"terms": [
    {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["1", "-2"], ["0", "2"]]}},
    {"offset": 1, "deriv_order": 1, "sign": "none", "coeff": {"coeffs": [["0", "2", "-2"]]}}]}}"#;

/// A fresh record file per call: the tests run in parallel threads.
fn record_file() -> PathBuf {
    static CALLS: AtomicUsize = AtomicUsize::new(0);
    let dir = std::env::temp_dir().join(format!("polytool-oeis-text-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let file = dir.join(format!(
        "derivative-{}.json",
        CALLS.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::write(&file, DERIVATIVE).unwrap();
    file
}

fn generate(extra: &[&str]) -> Output {
    let file = record_file();
    Command::new(env!("CARGO_BIN_EXE_polytool"))
        .args([
            "recurrence-generate",
            "--recurrence",
            file.to_str().unwrap(),
        ])
        .args(extra)
        .output()
        .expect("spawn polytool")
}

fn stdout(output: &Output) -> String {
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout.clone()).unwrap()
}

const A156919: &[&str] = &[
    "--first-row",
    "1",
    "--offset",
    "0",
    "--reversed",
    "--regular-triangle",
    "--rows",
    "5",
];

#[test]
fn example_formula_and_programs_follow_the_oeis_alignment() {
    let with = |format: &str| {
        let mut args = vec!["--format", format];
        args.extend_from_slice(A156919);
        stdout(&generate(&args))
    };
    assert_eq!(
        with("oeis-example"),
        "Triangle begins:\n  0:  1;\n  1:  2,   1;\n  2:  4,  10,   1;\n  \
         3:  8,  60,  36,   1;\n  4: 16, 296, 516, 116, 1;\n"
    );
    assert_eq!(
        with("oeis-formula"),
        "Row n lists the coefficients of P(n,x) in descending powers of x, where \
         P(n,x) = (1+2*n*x)*P(n-1,x) + (2*x-2*x^2)*d/dx P(n-1,x) for n >= 1, with P(0,x) = 1.\n"
    );
    assert!(with("mathematica").ends_with(
        "Table[Reverse[PadRight[CoefficientList[P[n, x], x], n + 1]], {n, 0, 4}] // Flatten\n"
    ));
    assert!(with("pari").ends_with("T(4)\n"));
    let json: serde_json::Value = serde_json::from_str(&with("oeis-json")).unwrap();
    assert_eq!(json["schema"], "polytool.oeis-text.v1");
    assert_eq!(json["rows_shown"], 5);
    assert_eq!(json["triangle_example"], with("oeis-example"));
}

#[test]
fn without_offset_the_rows_keep_their_own_index() {
    let formula = stdout(&generate(&["--format", "oeis-formula"]));
    assert_eq!(
        formula,
        "The row polynomials satisfy P(n,x) = (1+(2*n-2)*x)*P(n-1,x) \
         + (2*x-2*x^2)*d/dx P(n-1,x) for n >= 1, with P(0,x) = 1.\n"
    );
}

#[test]
fn bfile_from_the_first_row_matches_the_example() {
    let bfile = stdout(&generate(&[
        "--format",
        "bfile",
        "--first-row",
        "1",
        "--offset",
        "0",
        "--reversed",
        "--regular-triangle",
        "--max-terms",
        "15",
    ]));
    let terms: Vec<&str> = bfile
        .lines()
        .map(|line| line.split(' ').nth(1).unwrap())
        .collect();
    assert_eq!(
        terms,
        ["1", "2", "1", "4", "10", "1", "8", "60", "36", "1", "16", "296", "516", "116", "1"]
    );
    assert!(bfile.starts_with("0 1\n1 2\n"));
}

#[test]
fn reading_options_need_a_reading_format() {
    let output = generate(&["--rows", "3", "--first-row", "1"]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("--first-row"));
    let output = generate(&["--format", "pari", "--max-terms", "3"]);
    assert!(!output.status.success());
}
