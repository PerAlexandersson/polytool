//! Family CLI contracts: canonical generation, validation and OEIS draft export.
use crate::common::{run_polytool_snapshot_with_deadline as run, stdout};
use serde_json::{json, Value};

#[test]
fn family_commands_share_exact_generators() {
    let output = run(
        &["families", "generate", "PT-F0001", "--count", "5", "--json"],
        "",
    );
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let rows: Value = serde_json::from_str(&stdout(&output)).unwrap();
    assert_eq!(rows["rows"][4], json!(["1", "53", "3"]));
    let output = run(&["families", "list", "--query", "PT-F", "--json"], "");
    let list: Value = serde_json::from_str(&stdout(&output)).unwrap();
    assert_eq!(list["count"], 2);
    let output = run(
        &[
            "families",
            "export-oeis",
            "PT-F0002",
            "--count",
            "9",
            "--json",
        ],
        "",
    );
    let draft: Value = serde_json::from_str(&stdout(&output)).unwrap();
    assert!(draft["blocks"].as_str().unwrap().contains("%S A000000"));
    assert!(draft["proposal_toml"].as_str().unwrap().contains("[data]"));
    assert!(draft["bfile"].as_str().unwrap().starts_with("0 1\n1 1\n"));
    for arguments in [
        vec!["families", "show", "PT-F9999"],
        vec!["families", "generate", "PT-F0001", "--count", "513"],
        vec!["families", "generate", "PT-F0001", "--format", "invalid"],
    ] {
        assert_eq!(run(&arguments, "").status.code(), Some(2));
    }
    let output = run(&["families", "validate", "--json"], "");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn bfile_creation_refuses_overwrites_outside_the_fixture_directory() {
    let path = std::env::temp_dir().join(format!(
        "polytool-family-export-{}-{}.txt",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    let arguments = [
        "families",
        "export-oeis",
        "PT-F0001",
        "--count",
        "6",
        "--bfile",
        path.to_str().unwrap(),
        "--json",
    ];
    let first = run(&arguments, "");
    assert!(
        first.status.success(),
        "{}",
        String::from_utf8_lossy(&first.stderr)
    );
    let before = std::fs::read(&path).unwrap();
    let second = run(&arguments, "");
    assert_eq!(second.status.code(), Some(2));
    assert_eq!(std::fs::read(&path).unwrap(), before);
    std::fs::remove_file(path).unwrap();
}
