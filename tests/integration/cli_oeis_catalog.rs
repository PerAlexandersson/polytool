//! `polytool oeis add` and `polytool oeis check` against a tiny OEIS fixture
//! (`tests/fixtures/oeis-data`) and a scratch catalog data file.

use crate::common::{run_polytool_catalog as run_polytool, stdout};
use std::path::{Path, PathBuf};

fn fixture_data() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/oeis-data")
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

/// The Eulerian rule `P_n = (1 - t + n t) P_(n-1) + (t - t^2) P'_(n-1)`,
/// `P_1 = 1`, as a `polytool.recurrence.v1` object; `last` is the coefficient
/// of `n t` in the first term.
fn eulerian_recurrence_json(last: &str) -> String {
    format!(
        r#"{{"schema":"polytool.recurrence.v1","first_index":1,"initial_polynomials":[["1"]],
        "recurrence":{{"terms":[
          {{"offset":1,"deriv_order":0,"sign":"none","coeff":{{"coeffs":[["1","-1"],["0","{last}"]]}}}},
          {{"offset":1,"deriv_order":1,"sign":"none","coeff":{{"coeffs":[["0","1","-1"]]}}}}]}}}}"#
    )
}

#[test]
fn oeis_add_validates_a_supplied_recurrence_and_rejects_a_wrong_one() {
    let catalog = scratch_catalog("oeis-add-supplied");
    let catalog_arg = catalog.to_str().unwrap();
    let data = fixture_data();
    let data_arg = data.to_str().unwrap();
    let directory = catalog.parent().unwrap();

    // A wrong coefficient: the second row would be 1 + 3t.
    let wrong = directory.join("wrong.json");
    std::fs::write(&wrong, eulerian_recurrence_json("2")).unwrap();
    let before = std::fs::read_to_string(&catalog).unwrap();
    let output = run_polytool(&[
        "oeis",
        "add",
        "A008292",
        "--oeis-data",
        data_arg,
        "--catalog",
        catalog_arg,
        "--recurrence-json",
        wrong.to_str().unwrap(),
    ]);
    assert_eq!(output.status.code(), Some(1), "{output:?}");
    assert!(
        stdout(&output).starts_with("A008292\tmismatch\tOEIS term a(3) = 1 but"),
        "{}",
        stdout(&output)
    );
    assert_eq!(std::fs::read_to_string(&catalog).unwrap(), before);

    // Search options do not apply to a supplied recurrence.
    let right = directory.join("right.json");
    std::fs::write(&right, eulerian_recurrence_json("1")).unwrap();
    let refused = run_polytool(&[
        "oeis",
        "add",
        "A008292",
        "--recurrence-json",
        right.to_str().unwrap(),
        "--max-rec-len",
        "2",
        "--oeis-data",
        data_arg,
        "--dry-run",
    ]);
    assert_eq!(refused.status.code(), Some(2));

    let output = run_polytool(&[
        "oeis",
        "add",
        "A008292",
        "--oeis-data",
        data_arg,
        "--catalog",
        catalog_arg,
        "--recurrence-json",
        right.to_str().unwrap(),
    ]);
    assert!(output.status.success(), "{output:?}");
    assert!(stdout(&output)
        .starts_with("A008292\tadded\tP(n) = (1 - t + nt) P(n-1) + (t - t^2) P'(n-1) [supplied;"));
    let text = std::fs::read_to_string(&catalog).unwrap();
    let lines = text.lines().collect::<Vec<_>>();
    assert_eq!(lines.len(), 2);
    let record: serde_json::Value = serde_json::from_str(lines[1]).unwrap();
    assert_eq!(record["id"], "A008292");
    assert_eq!(record["recurrence_first_index"], 1);
    assert_eq!(record["initial_rows"], serde_json::json!([["1"]]));
    let provenance = &record["provenance"];
    assert_eq!(provenance["source"], "supplied_recurrence");
    assert_eq!(provenance["status"], "validated");
    assert_eq!(provenance["verification_rows"], 0);
    assert_eq!(provenance["source_rows"], 10);
    assert_eq!(provenance["checked_terms"], 55);
    assert!(provenance.get("fit_rows").is_none());
    assert!(provenance.get("search_options").is_none());
    assert!(record.get("real_rootedness").is_none());

    let check = run_polytool(&[
        "oeis",
        "check",
        "A008292",
        "--catalog",
        catalog_arg,
        "--oeis-data",
        data_arg,
    ]);
    assert!(check.status.success(), "{check:?}");
    assert!(stdout(&check).contains("summary: 1 ok, 0 mismatched"));

    // The first row as an authoritative prefix row, the recurrence from row 2.
    let shifted = directory.join("shifted.json");
    std::fs::write(
        &shifted,
        eulerian_recurrence_json("1")
            .replace("\"first_index\":1", "\"first_index\":2")
            .replace("[[\"1\"]]", "[[\"1\",\"1\"]]"),
    )
    .unwrap();
    let fresh = scratch_catalog("oeis-add-supplied-prefix");
    let output = run_polytool(&[
        "oeis",
        "add",
        "A008292",
        "--oeis-data",
        data_arg,
        "--catalog",
        fresh.to_str().unwrap(),
        "--recurrence-json",
        shifted.to_str().unwrap(),
        "--prefix-rows",
        "1",
        "--dry-run",
        "--json",
    ]);
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let result = &report["results"][0];
    assert_eq!(result["outcome"], "added", "{report:#}");
    assert_eq!(result["record"]["prefix_rows"], serde_json::json!([["1"]]));
    assert_eq!(result["record"]["recurrence_first_width"], 2);
}

#[test]
fn oeis_info_reports_the_real_rootedness_flag() {
    let output = run_polytool(&["oeis", "info", "A001263", "--json"]);
    assert!(output.status.success());
    let info: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(info["real_rootedness"]["status"], "lean_proved");
    assert_eq!(info["real_rootedness"]["as_of"], "2026-10-03");
    let output = run_polytool(&["oeis", "info", "A157153"]);
    assert!(stdout(&output).contains("real-rootedness: unknown"));
    let output = run_polytool(&["oeis", "info", "A001263"]);
    assert!(stdout(&output).contains(
        "real-rootedness: lean_proved as of 2026-10-03 (recorded flag, not proved by polytool)"
    ));
}

/// A scratch catalog of three bundled records (which carry no evidence key).
fn scratch_catalog_of_three(name: &str) -> (PathBuf, Vec<String>) {
    let directory = Path::new(env!("CARGO_TARGET_TMPDIR")).join(name);
    std::fs::create_dir_all(&directory).unwrap();
    let data = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("data/oeis-catalog.jsonl"),
    )
    .unwrap();
    let mut lines = Vec::new();
    for id in ["A001263", "A008292", "A122431"] {
        let prefix = format!("{{\"id\":\"{id}\"");
        let line = data.lines().find(|line| line.starts_with(&prefix)).unwrap();
        assert!(
            !line.contains("\"evidence\""),
            "unexpected evidence in {id}"
        );
        lines.push(line.to_string());
    }
    let path = directory.join("oeis-catalog.jsonl");
    std::fs::write(&path, format!("{}\n", lines.join("\n"))).unwrap();
    (path, lines)
}

#[test]
fn oeis_evidence_write_changes_only_the_evidence_key() {
    let (catalog, lines) = scratch_catalog_of_three("oeis-evidence-write");
    let catalog_arg = catalog.to_str().unwrap();

    // Without --write nothing changes.
    let output = run_polytool(&[
        "oeis",
        "evidence",
        "A008292",
        "A122431",
        "--rows",
        "12",
        "--catalog",
        catalog_arg,
        "--as-of",
        "2026-11-02",
    ]);
    assert!(output.status.success(), "{output:?}");
    let text = stdout(&output);
    assert!(text.contains("A008292\treal_rooted\trows=12\t-"), "{text}");
    assert!(
        text.contains("A122431\tlog_concave\trows=12\treal_rooted@2,ultra_log_concave@2"),
        "{text}"
    );
    assert!(text.contains("lean_proved with evidence weaker than real_rooted: none"));
    assert_eq!(
        std::fs::read_to_string(&catalog).unwrap(),
        format!("{}\n", lines.join("\n"))
    );

    // --write adds the key after real_rootedness to the computed entries only.
    let output = run_polytool(&[
        "oeis",
        "evidence",
        "A008292",
        "A122431",
        "--rows",
        "12",
        "--catalog",
        catalog_arg,
        "--as-of",
        "2026-11-02",
        "--write",
    ]);
    assert!(output.status.success(), "{output:?}");
    let written = std::fs::read_to_string(&catalog).unwrap();
    let new_lines = written.lines().collect::<Vec<_>>();
    assert_eq!(new_lines.len(), 3);
    assert_eq!(new_lines[0], lines[0], "A001263 was not selected");
    let a008292 = format!(
        "{},\"evidence\":{{\"strongest\":\"real_rooted\",\"rows_checked\":12,\
         \"first_failure\":{{}},\"as_of\":\"2026-11-02\"}}}}",
        lines[1].strip_suffix('}').unwrap()
    );
    assert_eq!(new_lines[1], a008292);
    let a122431 = format!(
        "{},\"evidence\":{{\"strongest\":\"log_concave\",\"rows_checked\":12,\
         \"first_failure\":{{\"real_rooted\":2,\"ultra_log_concave\":2}},\
         \"as_of\":\"2026-11-02\"}}}}",
        lines[2].strip_suffix('}').unwrap()
    );
    assert_eq!(new_lines[2], a122431);
    assert!(new_lines[2].contains("\"real_rootedness\":{\"status\":\"known_not_real_rooted\""));

    // A second run with other rows replaces the key, again touching nothing else.
    let output = run_polytool(&[
        "oeis",
        "evidence",
        "A122431",
        "--rows",
        "5",
        "--catalog",
        catalog_arg,
        "--as-of",
        "2026-11-03",
        "--write",
        "--json",
    ]);
    assert!(output.status.success(), "{output:?}");
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["schema"], "polytool.oeis-evidence.v1");
    assert_eq!(report["written"], true);
    assert_eq!(report["distribution"]["log_concave"], 1);
    assert_eq!(
        report["entries"][0]["evidence"]["first_failure"],
        serde_json::json!({"real_rooted": 2, "ultra_log_concave": 2})
    );
    let rewritten = std::fs::read_to_string(&catalog).unwrap();
    let again = rewritten.lines().collect::<Vec<_>>();
    assert_eq!(again[0], lines[0]);
    assert_eq!(again[1], new_lines[1]);
    assert!(again[2].starts_with(lines[2].strip_suffix('}').unwrap()));
    assert!(again[2].ends_with(
        "\"evidence\":{\"strongest\":\"log_concave\",\"rows_checked\":5,\
         \"first_failure\":{\"real_rooted\":2,\"ultra_log_concave\":2},\
         \"as_of\":\"2026-11-03\"}}"
    ));
}

#[test]
fn oeis_evidence_rejects_bad_requests_and_noncanonical_files() {
    let (catalog, _) = scratch_catalog_of_three("oeis-evidence-errors");
    let catalog_arg = catalog.to_str().unwrap();
    for args in [
        vec!["oeis", "evidence", "A008292"],
        vec!["oeis", "evidence", "--rows", "5"],
        vec!["oeis", "evidence", "A008292", "--all", "--rows", "5"],
        vec!["oeis", "evidence", "A008292", "--rows", "0"],
        vec!["oeis", "evidence", "A999999", "--rows", "5"],
    ] {
        let mut full = args.clone();
        full.extend(["--catalog", catalog_arg]);
        let output = run_polytool(&full);
        assert_eq!(output.status.code(), Some(2), "{full:?}");
    }
    // A file that is not in canonical form (here: a space) is never rewritten.
    let text = std::fs::read_to_string(&catalog).unwrap();
    let spaced = text.replacen("\"id\":", "\"id\": ", 1);
    std::fs::write(&catalog, &spaced).unwrap();
    let output = run_polytool(&[
        "oeis",
        "evidence",
        "A008292",
        "--rows",
        "5",
        "--catalog",
        catalog_arg,
        "--write",
    ]);
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&output.stderr).contains("not in canonical form"));
    assert_eq!(std::fs::read_to_string(&catalog).unwrap(), spaced);
}

#[test]
fn oeis_info_reports_missing_evidence() {
    let output = run_polytool(&["oeis", "info", "A122431", "--json"]);
    assert!(output.status.success());
    let info: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert!(info["evidence"].is_null());
    let output = run_polytool(&["oeis", "info", "A122431"]);
    assert!(stdout(&output).contains("evidence: none recorded"));
}
