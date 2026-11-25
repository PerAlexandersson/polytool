//! Exponential row CLI schemas, sparse extraction, failure status and replay.
use crate::common::run_polytool as run;
use serde_json::{json, Value};
#[test]
fn sparse_row_and_counterexample_replay() {
    let result = run(
        &[
            "exp-rows",
            "--unmarked",
            "sizes:1,2",
            "--marked",
            "sizes:7",
            "--rows",
            "1000",
            "--json",
        ],
        "",
    );
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let value: Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(value["rows"][0]["degree"], 142);
    let result = run(
        &[
            "exp-rows",
            "--unmarked",
            "sizes:1",
            "--marked",
            "weights:0,2,0,0,0,720",
            "--rows",
            "0..=8",
            "--check",
            "first-failure",
            "--summary",
            "--json",
        ],
        "",
    );
    assert_eq!(result.status.code(), Some(1));
    let mut value: Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(value["first_failure"]["index"], 6);
    assert!(value["rows"][0]["coefficients"].is_null());
    let saved = serde_json::to_string(&value).unwrap();
    let replay = run(&["exp-rows", "--verify", "--input", "-", "--json"], &saved);
    assert!(
        replay.status.success(),
        "{}",
        String::from_utf8_lossy(&replay.stderr)
    );
    value["first_failure"]["distinct_real_roots"] = json!(3);
    let bad = run(
        &["exp-rows", "--verify", "--input", "-"],
        &value.to_string(),
    );
    assert_eq!(bad.status.code(), Some(2));
}
#[test]
fn json_input_pencil_and_invalid_requests() {
    let request = json!({"family":{"unmarked":{"kind":"sizes","sizes":[1]},"marked":{"kind":"sizes","sizes":[1]}},"rows":[2],"finite_rank":[{"index":2,"reference":["1","2","1"],"pencil":{"a":"1","b":"1","factor":["1","1"]}}]});
    let result = run(
        &["exp-rows", "--input", "-", "--json"],
        &request.to_string(),
    );
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let value: Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(value["finite_rank"][0]["certified"], true);
    for args in [
        vec!["exp-rows", "--marked", "sizes:0", "--rows", "1"],
        vec!["exp-rows", "--marked", "tail:2", "--rows", "3,1"],
        vec!["exp-rows", "--marked", "weights:1/0", "--rows", "1"],
    ] {
        assert_eq!(run(&args, "").status.code(), Some(2));
    }
    let result = run(
        &[
            "exp-rows",
            "--marked",
            "tail:1",
            "--rows",
            "1,3",
            "--check",
            "interlacing",
            "--json",
        ],
        "",
    );
    assert!(result.status.success());
    let v: Value = serde_json::from_slice(&result.stdout).unwrap();
    assert!(v["interlacing"].as_array().unwrap().is_empty());
}
