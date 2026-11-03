//! `polytool recurrence-generate --format bfile`: b-files from any recurrence,
//! and the `oeis generate --format bfile` output through the same formatter.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn run_polytool(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_polytool"))
        .args(args)
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

/// The recurrence JSON of a catalog entry, written to a file of its own.
fn recurrence_file(id: &str) -> PathBuf {
    let info = run_polytool(&["oeis", "info", id, "--json"]);
    let value: serde_json::Value = serde_json::from_slice(&info.stdout).unwrap();
    let dir = std::env::temp_dir().join(format!("polytool-bfile-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let file = dir.join(format!("{id}.json"));
    std::fs::write(
        &file,
        serde_json::to_vec(&value["recurrence_data"]).unwrap(),
    )
    .unwrap();
    file
}

fn generate(file: &Path, extra: &[&str]) -> Output {
    let mut args = vec![
        "recurrence-generate",
        "--recurrence",
        file.to_str().unwrap(),
    ];
    args.extend_from_slice(extra);
    run_polytool(&args)
}

#[test]
fn eulerian_bfile_matches_the_catalog_bfile_for_every_stopping_cap() {
    let file = recurrence_file("A008292");
    // Row n has n terms, so 45 terms are exactly 9 rows.
    for cap in ["1", "2", "5", "6", "44", "45", "46"] {
        let from_recurrence = stdout(&generate(
            &file,
            &["--format", "bfile", "--offset", "1", "--max-terms", cap],
        ));
        let from_catalog = stdout(&run_polytool(&[
            "oeis",
            "generate",
            "A008292",
            "--rows",
            "12",
            "--format",
            "bfile",
            "--max-terms",
            cap,
        ]));
        assert_eq!(from_recurrence, from_catalog, "cap {cap}");
        assert!(!from_recurrence.is_empty());
    }
}

#[test]
fn first_terms_and_row_counts_agree_with_the_catalog() {
    let file = recurrence_file("A008292");
    let text = stdout(&generate(
        &file,
        &["--format", "bfile", "--offset", "1", "--rows", "4"],
    ));
    assert_eq!(
        text,
        "1 1\n2 1\n3 1\n4 1\n5 4\n6 1\n7 1\n8 11\n9 11\n10 1\n"
    );
    let catalog = stdout(&run_polytool(&[
        "oeis", "generate", "A008292", "--rows", "4", "--format", "bfile",
    ]));
    assert_eq!(text, catalog);
    // `--additional` keeps the initial rows and appends: one seed row here.
    let additional = stdout(&generate(
        &file,
        &["--format", "bfile", "--offset", "1", "--additional", "3"],
    ));
    assert_eq!(additional, text);
}

#[test]
fn reversed_rows_and_the_nearest_row_rule() {
    let file = recurrence_file("A008292");
    let reversed = stdout(&generate(
        &file,
        &[
            "--format",
            "bfile",
            "--offset",
            "0",
            "--rows",
            "4",
            "--reversed",
        ],
    ));
    // Eulerian rows are palindromic, so reversal changes nothing.
    assert_eq!(
        reversed
            .lines()
            .map(|l| l.split(' ').nth(1).unwrap())
            .collect::<Vec<_>>(),
        ["1", "1", "1", "1", "4", "1", "1", "11", "11", "1"]
    );
    // Prefix sizes 1, 3, 6, 10: a cap of 5 is nearer to 6 than to 3.
    let at_most = stdout(&generate(&file, &["--format", "bfile", "--max-terms", "5"]));
    assert_eq!(at_most.lines().count(), 3);
    let nearest = stdout(&generate(
        &file,
        &["--format", "bfile", "--max-terms", "5", "--nearest-row"],
    ));
    assert_eq!(nearest.lines().count(), 6);
    assert!(nearest.starts_with("0 1\n1 1\n"));
}

#[test]
fn a_nonpalindromic_triangle_shows_the_orientation() {
    // A hand-written recurrence P(n) = (1 + t) P(n-1) with P(0) = 1 + 2t, so
    // the rows are 1 2 / 1 3 2 / 1 4 5 2 (not palindromic).
    let json = serde_json::json!({
        "schema": "polytool.recurrence.v1",
        "first_index": 0,
        "initial_polynomials": [["1", "2"]],
        "recurrence": {
            "terms": [{
                "offset": 1, "deriv_order": 0, "sign": "none",
                "coeff": {"coeffs": [["1", "1"]]}
            }]
        }
    });
    let dir = std::env::temp_dir().join(format!("polytool-bfile-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let file = dir.join("shift.json");
    std::fs::write(&file, serde_json::to_vec(&json).unwrap()).unwrap();
    let values = |args: &[&str]| -> Vec<String> {
        let mut all = vec!["--format", "bfile", "--rows", "3"];
        all.extend_from_slice(args);
        stdout(&generate(&file, &all))
            .lines()
            .map(|line| line.split(' ').nth(1).unwrap().to_string())
            .collect()
    };
    assert_eq!(values(&[]), ["1", "2", "1", "3", "2", "1", "4", "5", "2"]);
    assert_eq!(
        values(&["--reversed"]),
        ["2", "1", "2", "3", "1", "2", "5", "4", "1"]
    );
}

#[test]
fn bfile_options_are_validated() {
    let file = recurrence_file("A008292");
    // A cap smaller than the first row, with whole rows only.
    let tiny = generate(
        &file,
        &["--format", "bfile", "--rows", "3", "--max-terms", "0"],
    );
    assert_eq!(tiny.status.code(), Some(2));
    // Neither a row count nor a cap.
    let unbounded = generate(&file, &["--format", "bfile"]);
    assert_eq!(unbounded.status.code(), Some(2));
    // B-file options without the format.
    let stray = generate(&file, &["--rows", "3", "--offset", "1"]);
    assert_eq!(stray.status.code(), Some(2));
    // The default text format is unchanged.
    assert_eq!(stdout(&generate(&file, &["--rows", "2"])), "1\n1, 1\n");
}

#[test]
fn oeis_generate_bfile_is_byte_identical_to_whole_row_flattening() {
    // The format `oeis generate --format bfile` had before it used the library
    // formatter: `index value` lines, rows kept whole, at most the cap.
    for (id, rows, cap) in [
        ("A008292", 20usize, None),
        ("A008292", 20, Some(37usize)),
        ("A390883", 12, None),
        ("A390883", 12, Some(25)),
    ] {
        let entry = polytool::oeis::by_id(id).unwrap();
        let mut expected = String::new();
        let mut index = entry.flattened_offset;
        let mut terms = 0usize;
        for row in entry.generate_rows(rows).unwrap() {
            if cap.is_some_and(|limit| terms + row.len() > limit) {
                break;
            }
            for value in &row {
                expected.push_str(&format!("{index} {value}\n"));
                index += 1;
            }
            terms += row.len();
        }
        let rows_text = rows.to_string();
        let mut args = vec![
            "oeis", "generate", id, "--rows", &rows_text, "--format", "bfile",
        ];
        let cap_text = cap.map(|limit| limit.to_string());
        if let Some(cap_text) = &cap_text {
            args.extend_from_slice(&["--max-terms", cap_text]);
        }
        assert_eq!(stdout(&run_polytool(&args)), expected, "{id} cap {cap:?}");
    }
}
