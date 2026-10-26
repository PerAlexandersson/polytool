//! `polytool oeis add` and `polytool oeis check` against a tiny OEIS fixture
//! (`tests/fixtures/oeis-data`) and a scratch catalog data file.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn fixture_data() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/oeis-data")
}

fn run_polytool(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_polytool"))
        .args(args)
        .env_remove("OEIS_DATA_DIR")
        .output()
        .expect("spawn polytool")
}

fn stdout(output: &Output) -> String {
    String::from_utf8(output.stdout.clone()).unwrap()
}

/// A scratch catalog holding only the first record of the bundled data file.
fn scratch_catalog(name: &str) -> PathBuf {
    let directory = Path::new(env!("CARGO_TARGET_TMPDIR")).join(name);
    std::fs::create_dir_all(&directory).unwrap();
    let data = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("data/oeis-catalog.jsonl"),
    )
    .unwrap();
    let first = data.lines().next().unwrap();
    assert!(first.contains("\"id\":\"A001263\""));
    let path = directory.join("oeis-catalog.jsonl");
    std::fs::write(&path, format!("{first}\n")).unwrap();
    path
}

#[test]
fn oeis_add_records_a_verified_entry_with_provenance() {
    let catalog = scratch_catalog("oeis-add-records");
    let catalog_arg = catalog.to_str().unwrap();
    let data = fixture_data();
    let data_arg = data.to_str().unwrap();

    let output = run_polytool(&[
        "oeis",
        "add",
        "A008292",
        "--oeis-data",
        data_arg,
        "--catalog",
        catalog_arg,
    ]);
    assert!(output.status.success(), "{:?}", output);
    assert!(stdout(&output).starts_with("A008292\tadded\tP(n) = (1 - t + nt) P(n-1)"));

    let text = std::fs::read_to_string(&catalog).unwrap();
    let lines = text.lines().collect::<Vec<_>>();
    assert_eq!(lines.len(), 2);
    assert!(lines[0].contains("\"id\":\"A001263\""));
    let record: serde_json::Value = serde_json::from_str(lines[1]).unwrap();
    assert_eq!(record["id"], "A008292");
    assert_eq!(record["first_row"], 1);
    let provenance = &record["provenance"];
    assert_eq!(provenance["source"], "polytool_recurrence_search");
    assert_eq!(provenance["status"], "verified");
    assert_eq!(provenance["source_rows"], 10);
    assert_eq!(provenance["fit_rows"], 8);
    assert_eq!(provenance["verification_rows"], 2);
    assert_eq!(provenance["checked_terms"], 55);
    assert_eq!(provenance["terms_source"], "oeis_internal_format");
    assert_eq!(provenance["search_options"]["max_rec_len"], 5);

    let again = run_polytool(&[
        "oeis",
        "add",
        "A008292",
        "--oeis-data",
        data_arg,
        "--catalog",
        catalog_arg,
    ]);
    assert!(again.status.success());
    assert!(stdout(&again).starts_with("A008292\talready_present\t"));

    let check = run_polytool(&[
        "oeis",
        "check",
        "--catalog",
        catalog_arg,
        "--oeis-data",
        data_arg,
        "--json",
    ]);
    assert!(check.status.success());
    let report: serde_json::Value = serde_json::from_slice(&check.stdout).unwrap();
    assert_eq!(report["schema"], "polytool.oeis-check.v1");
    assert_eq!(report["checked"], 1);
    assert_eq!(report["unchecked"], 1);
    assert_eq!(report["results"][1]["id"], "A008292");
    assert_eq!(report["results"][1]["checked_terms"], 55);
}

#[test]
fn oeis_add_reports_batch_outcomes_without_writing_on_dry_run() {
    let catalog = scratch_catalog("oeis-add-batch");
    let before = std::fs::read_to_string(&catalog).unwrap();
    let ids = catalog.with_file_name("ids.txt");
    std::fs::write(
        &ids,
        "# batch\nA001263\nA000045 Fibonacci\nA008292\nA999999\n",
    )
    .unwrap();
    let data = fixture_data();
    let output = run_polytool(&[
        "oeis",
        "add",
        "--ids-file",
        ids.to_str().unwrap(),
        "--oeis-data",
        data.to_str().unwrap(),
        "--catalog",
        catalog.to_str().unwrap(),
        "--dry-run",
        "--json",
    ]);
    assert_eq!(output.status.code(), Some(1));
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["schema"], "polytool.oeis-add.v1");
    let outcomes = report["results"]
        .as_array()
        .unwrap()
        .iter()
        .map(|result| {
            (
                result["id"].as_str().unwrap().to_string(),
                result["outcome"].as_str().unwrap().to_string(),
            )
        })
        .collect::<Vec<_>>();
    let expected = [
        ("A001263", "already_present"),
        ("A000045", "layout_unsupported"),
        ("A008292", "added"),
        ("A999999", "unavailable"),
    ];
    assert_eq!(
        outcomes,
        expected
            .iter()
            .map(|(id, outcome)| (id.to_string(), outcome.to_string()))
            .collect::<Vec<_>>()
    );
    assert_eq!(report["results"][2]["record"]["id"], "A008292");
    assert_eq!(std::fs::read_to_string(&catalog).unwrap(), before);
}

#[test]
fn oeis_add_accepts_explicit_row_width_rules_and_bfiles() {
    let catalog = scratch_catalog("oeis-add-rules");
    let catalog_arg = catalog.to_str().unwrap();
    let data = fixture_data();
    let data_arg = data.to_str().unwrap();

    let table = run_polytool(&[
        "oeis",
        "add",
        "A055151",
        "--layout",
        "table",
        "--row-widths",
        "1+1/2",
        "--oeis-data",
        data_arg,
        "--catalog",
        catalog_arg,
        "--dry-run",
    ]);
    assert!(table.status.success(), "{:?}", table);
    let text = stdout(&table);
    assert!(text.starts_with("A055151\tadded\t"), "{text}");
    let record: serde_json::Value = serde_json::from_str(text.lines().nth(1).unwrap()).unwrap();
    assert_eq!(record["layout"], "table");

    let bfile = catalog.with_file_name("b008292.txt");
    let terms = [
        1, 1, 1, 1, 4, 1, 1, 11, 11, 1, 1, 26, 66, 26, 1, 1, 57, 302, 302, 57, 1, 1, 120, 1191,
        2416, 1191, 120, 1, 1, 247, 4293, 15619, 15619, 4293, 247, 1,
    ];
    let mut lines = String::from("# A008292 excerpt\n");
    for (index, term) in terms.iter().enumerate() {
        lines.push_str(&format!("{} {term}\n", index + 1));
    }
    std::fs::write(&bfile, lines).unwrap();
    let explicit = run_polytool(&[
        "oeis",
        "add",
        "A008292",
        "--bfile",
        bfile.to_str().unwrap(),
        "--name",
        "Eulerian triangle",
        "--layout",
        "triangle",
        "--catalog",
        catalog_arg,
        "--dry-run",
        "--json",
    ]);
    assert!(explicit.status.success(), "{:?}", explicit);
    let report: serde_json::Value = serde_json::from_slice(&explicit.stdout).unwrap();
    let record = &report["results"][0]["record"];
    assert_eq!(record["name"], "Eulerian triangle");
    assert_eq!(record["provenance"]["terms_source"], "bfile");
    assert_eq!(record["provenance"]["checked_terms"], 36);

    let missing_rule = run_polytool(&[
        "oeis",
        "add",
        "A055151",
        "--layout",
        "table",
        "--oeis-data",
        data_arg,
        "--dry-run",
    ]);
    assert_eq!(missing_rule.status.code(), Some(2));
}

#[test]
fn oeis_check_validates_bundled_entries_against_local_records() {
    let data = fixture_data();
    let output = run_polytool(&[
        "oeis",
        "check",
        "A008292",
        "A055151",
        "A001263",
        "--oeis-data",
        data.to_str().unwrap(),
    ]);
    assert!(output.status.success(), "{:?}", output);
    let text = stdout(&output);
    assert!(text.contains("A008292\tok\t55 terms match"), "{text}");
    assert!(text.contains("A055151\tok\t"), "{text}");
    assert!(text.contains("A001263\tunchecked\t"), "{text}");
    assert!(text.ends_with("summary: 2 ok, 0 mismatched, 1 unchecked (3 entries)\n"));
}
