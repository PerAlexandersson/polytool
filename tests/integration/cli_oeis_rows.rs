//! Local OEIS row extraction.

#[test]
fn separate_bfile_directory_extends_rows_and_preserves_validation() {
    use polytool::oeis::{read_oeis_bfile, OeisSourceError};

    let root = std::env::temp_dir().join(format!("polytool-rows-bfile-dir-{}", std::process::id()));
    let cache = root.join("bfiles");
    std::fs::create_dir_all(root.join("seq/A000")).unwrap();
    std::fs::create_dir_all(root.join("files/A000")).unwrap();
    std::fs::create_dir_all(cache.join("A000")).unwrap();
    std::fs::write(
        root.join("seq/A000/A000001.seq"),
        "%N A000001 Example\n%O A000001 0,1\n%K A000001 tabl\n%S A000001 1,2,3\n",
    )
    .unwrap();
    let fallback = root.join("files/A000/b000001.txt");
    let override_file = cache.join("A000/b000001.txt");
    let pointer = "version https://git-lfs.github.com/spec/v1\noid sha256:0000000000000000000000000000000000000000000000000000000000000000\nsize 24\n";
    std::fs::write(&fallback, pointer).unwrap();
    let run = |with_cache: bool| {
        let mut command = crate::common::polytool_command();
        command.args([
            "oeis",
            "rows",
            "A000001",
            "--data",
            root.to_str().unwrap(),
            "--json",
        ]);
        if with_cache {
            command.args(["--bfile-dir", cache.to_str().unwrap()]);
        }
        let output = command.output().unwrap();
        assert!(output.status.success(), "{:?}", output);
        serde_json::from_slice::<serde_json::Value>(&output.stdout).unwrap()
    };
    let baseline = run(false);
    assert_eq!(run(true), baseline);
    assert_eq!(baseline["bfile_rejection"], "invalid_bfile");
    assert!(matches!(
        polytool::oeis::parse_bfile(
            &read_oeis_bfile(&root, "A000001", Some(&cache))
                .unwrap()
                .unwrap()
        ),
        Err(OeisSourceError::GitLfsPointer)
    ));
    let extension = "0 1\n1 2\n2 3\n3 4\n4 5\n5 6\n";
    std::fs::write(&override_file, extension).unwrap();
    assert_eq!(
        read_oeis_bfile(&root, "A000001", Some(&cache))
            .unwrap()
            .as_deref(),
        Some(extension)
    );
    assert_eq!(
        run(true)["rows"],
        serde_json::json!([["1"], ["2", "3"], ["4", "5", "6"]])
    );
    assert_eq!(run(true)["term_source"], "bfile");
    assert!(run(true)["bfile_rejection"].is_null());
    // A present override wins even when the fallback is usable. Rejection
    // retains the internal terms instead of silently trying another source.
    std::fs::write(&fallback, extension).unwrap();
    for (text, reason) in [
        ("0 9\n1 2\n2 3\n", "prefix_mismatch"),
        ("0 1\n2 2\n3 3\n", "non_consecutive_indices"),
    ] {
        std::fs::write(&override_file, text).unwrap();
        let report = run(true);
        assert_eq!(report["bfile_rejection"], reason);
        assert_eq!(report["term_source"], "oeis_internal_format");
        assert_eq!(report["rows"].as_array().unwrap().len(), 2);
    }
    std::fs::write(&override_file, [0xff]).unwrap();
    let report = run(true);
    assert_eq!(report["bfile_rejection"], "unreadable");
    assert!(report["bfile_rejection_detail"].is_string());
    std::fs::remove_file(&override_file).unwrap();
    assert_eq!(run(true)["term_source"], "bfile");
    assert_eq!(run(true)["rows"].as_array().unwrap().len(), 3);
    std::fs::remove_file(&fallback).unwrap();
    assert!(read_oeis_bfile(&root, "A000001", Some(&cache))
        .unwrap()
        .is_none());
    assert_eq!(run(true), run(false));
    for id in ["", "A12", "A00000/", "a000001", "é00000"] {
        assert_eq!(
            read_oeis_bfile(&root, id, Some(&cache)).unwrap_err().kind(),
            std::io::ErrorKind::InvalidInput
        );
    }
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn local_rows_and_refusal() {
    let root = std::env::temp_dir().join(format!("polytool-rows-{}", std::process::id()));
    std::fs::create_dir_all(root.join("seq/A000")).unwrap();
    let path = root.join("seq/A000/A000001.seq");
    std::fs::write(
        &path,
        "%N A000001 Example\n%O A000001 0,1\n%K A000001 tabl\n%S A000001 1,2,3,4,5,6\n",
    )
    .unwrap();
    let run = |extra: &[&str]| {
        crate::common::polytool_command()
            .args(["oeis", "rows", "A000001", "--data", root.to_str().unwrap()])
            .args(extra)
            .output()
            .unwrap()
    };
    let result = run(&[]);
    assert!(result.status.success());
    assert_eq!(result.stdout, b"1\n2 3\n4 5 6\n");
    let json: serde_json::Value = serde_json::from_slice(&run(&["--json"]).stdout).unwrap();
    assert_eq!(json["term_source"], "oeis_internal_format");
    assert_eq!(json["status"], "ok");
    std::fs::write(
        &path,
        "%N A000001 Example\n%O A000001 0,1\n%K A000001 nonn\n%S A000001 1,2\n",
    )
    .unwrap();
    let json: serde_json::Value = serde_json::from_slice(&run(&["--json"]).stdout).unwrap();
    assert_eq!(json["reason"], "not tabl or tabf");
    assert_eq!(json["reason_code"], "not_triangle_keyword");
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn local_table_rules_and_bfile_validation() {
    let root = std::env::temp_dir().join(format!("polytool-rows-table-{}", std::process::id()));
    std::fs::create_dir_all(root.join("seq/A000")).unwrap();
    std::fs::create_dir_all(root.join("files/A000")).unwrap();
    let path = root.join("seq/A000/A000001.seq");
    let bpath = root.join("files/A000/b000001.txt");
    let run = || {
        let output = crate::common::polytool_command()
            .args([
                "oeis",
                "rows",
                "A000001",
                "--data",
                root.to_str().unwrap(),
                "--json",
            ])
            .output()
            .unwrap();
        assert!(output.status.success());
        serde_json::from_slice::<serde_json::Value>(&output.stdout).unwrap()
    };
    let write = |terms: &str, examples: &str| {
        std::fs::write(&path, format!("%N A000001 Example\n%O A000001 0,1\n%K A000001 tabf\n%S A000001 {terms}\n{examples}")).unwrap()
    };
    write(
        "1,2,3,4,5,6,7,8,9",
        "%e A000001 n=0: 1\n%e A000001 n=1: 2,3,4\n%e A000001 n=2: 5,6,7,8,9\n",
    );
    assert_eq!(run()["row_length_rule"], "+2 per row");
    write(
        "1,2,3,4,5,6,7,8,9",
        "%e A000001 1\n%e A000001 2\n%e A000001 3,4\n%e A000001 5,6\n%e A000001 7,8,9\n",
    );
    assert_eq!(
        run()["row_length_rule"],
        "+1 every 2 rows (transition phase 1)"
    );
    assert_eq!(
        run()["row_length_rule_details"],
        serde_json::json!({"kind":"periodic_step","step":1,"period":2,"phase":1})
    );
    write(
        "1,2,3,4,5,6,7",
        "%e A000001 1\n%e A000001 2,3\n%e A000001 4,5,6,7\n",
    );
    assert_eq!(run()["reason"], "irregular lengths");
    assert_eq!(run()["reason_code"], "irregular_lengths");
    std::fs::write(
        &path,
        "%N A000001 Example\n%O A000001 0,1\n%K A000001 tabl\n%S A000001 1,2,3\n",
    )
    .unwrap();
    std::fs::write(&bpath, "0 1\n1 2\n2 3\n3 4\n").unwrap();
    assert_eq!(run()["term_source"], "bfile");
    std::fs::write(&bpath, "0 1\n2 2\n3 3\n").unwrap();
    assert_eq!(run()["term_source"], "oeis_internal_format");
    std::fs::write(&bpath, "0 9\n1 2\n2 3\n").unwrap();
    assert_eq!(run()["term_source"], "oeis_internal_format");
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn diagnostics_labels_and_structured_rules() {
    let root =
        std::env::temp_dir().join(format!("polytool-rows-diagnostics-{}", std::process::id()));
    std::fs::create_dir_all(root.join("seq/A000")).unwrap();
    std::fs::create_dir_all(root.join("files/A000")).unwrap();
    let seq = root.join("seq/A000/A000001.seq");
    let bfile = root.join("files/A000/b000001.txt");
    let run = |flags: &[&str]| {
        crate::common::polytool_command()
            .args(["oeis", "rows", "A000001", "--data", root.to_str().unwrap()])
            .args(flags)
            .output()
            .unwrap()
    };
    let json = || serde_json::from_slice::<serde_json::Value>(&run(&["--json"]).stdout).unwrap();
    let write = |terms: &str, table: &str| {
        std::fs::write(
            &seq,
            format!(
                "%N A000001 Example\n%O A000001 0,1\n%K A000001 tabf\n%S A000001 {terms}\n{table}"
            ),
        )
        .unwrap()
    };
    let missing = run(&[]);
    assert_eq!(missing.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&missing.stderr).contains("--data/seq/Annn/Annnnnn.seq"));
    write(
        "1,2,3,4,5,6",
        "%e A000001 [0] 1;\n%e A000001 1 | 2,3;\n%e A000001 n=2: 4,5,6.\n",
    );
    assert_eq!(
        json()["row_length_rule_details"],
        serde_json::json!({"kind":"constant_step","step":1,"period":1,"phase":0})
    );
    write("1,2,3,4,5,6", "%e A000001 1\n%e A000001 2,3\n");
    assert_eq!(json()["reason_code"], "too_few_example_rows");
    write("1,2,3,4,5,6", "%e A000001 2,3\n%e A000001 4,5,6\n");
    assert_eq!(json()["reason_code"], "table_does_not_start_at_first_term");
    write(
        "1,2,3,4,5,6",
        "%e A000001 1\n%e A000001 2,3\n%e A000001 4,9,6\n",
    );
    assert_eq!(json()["reason_code"], "terms_disagree_with_table");
    write("1,2,3,4,5,6", "");
    assert_eq!(json()["reason_code"], "no_example_table");
    std::fs::write(&bfile, "0 1\r\n  \r\n1 2\r\n2 3\r\n3 4\r\n4 5\r\n5 6\r\n").unwrap();
    assert_eq!(json()["term_source"], "bfile");
    std::fs::write(&bfile, "version https://git-lfs.github.com/spec/v1\n").unwrap();
    assert_eq!(json()["bfile_rejection"], "invalid_bfile");
    assert!(String::from_utf8_lossy(&run(&[]).stderr).contains("not a valid b-file"));
    std::fs::write(&bfile, "0 1\n2 2\n").unwrap();
    assert_eq!(json()["bfile_rejection"], "non_consecutive_indices");
    std::fs::write(&bfile, "0 9\n1 2\n").unwrap();
    assert_eq!(json()["bfile_rejection"], "prefix_mismatch");
    std::fs::write(&bfile, [0xff]).unwrap();
    assert_eq!(json()["bfile_rejection"], "unreadable");
    assert!(String::from_utf8_lossy(&run(&[]).stderr).contains("not UTF-8"));
    std::fs::remove_dir_all(root).unwrap();
}
