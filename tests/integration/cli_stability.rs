//! CLI stability syntax, exact counts, replay and resource bounds.

use crate::common::run_polytool as run;
use serde_json::{json, Value};

fn value(args: &[&str], input: &str) -> Value {
    let output = run(args, input);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).unwrap()
}

#[test]
fn regions_gaussian_batch_and_zero() {
    let r = value(&["stability", "(1+t)^3", "--region", "left", "--json"], "");
    assert_eq!(
        r["items"][0]["counts"],
        json!({"inside":3,"boundary":0,"outside":0})
    );
    for (region, counts) in [
        ("upper", [0, 3, 0]),
        ("lower", [0, 3, 0]),
        ("right", [0, 0, 3]),
        ("disk", [0, 3, 0]),
        ("exterior", [0, 3, 0]),
    ] {
        let r = value(&["stability", "(1+t)^3", "--region", region, "--json"], "");
        assert_eq!(
            r["items"][0]["counts"],
            json!({"inside":counts[0],"boundary":counts[1],"outside":counts[2]})
        );
    }
    let r = value(&["stability", "--json"], "# comment\nt-I\nt+I\n0\n5\n");
    assert_eq!(r["items"][0]["counts"]["inside"], 1);
    assert_eq!(r["items"][1]["counts"]["inside"], 0);
    assert_eq!(r["items"][2]["degree"], Value::Null);
    assert_eq!(r["items"][2]["counts"], Value::Null);
    assert_eq!(r["items"][2]["verdict"], true);
    assert_eq!(r["items"][3]["verdict"], true);
    let disk = value(
        &[
            "stability",
            "t-I",
            "--disk-center",
            "I",
            "--radius",
            "1",
            "--json",
        ],
        "",
    );
    assert_eq!(disk["items"][0]["verdict"], true);
    let mobius = value(
        &[
            "stability",
            "t-I",
            "--mobius",
            "1,0,0,1",
            "--convention",
            "inside",
            "--json",
        ],
        "",
    );
    assert_eq!(mobius["items"][0]["verdict"], true);
}

#[test]
fn replay_and_tampering() {
    let r = value(&["stability", "t^2+(1+2I)t-I", "--format", "full"], "");
    assert_eq!(r["schema"], "polytool.stability.v1");
    assert_eq!(
        value(&["stability", "verify"], &r.to_string())["verified"],
        true
    );
    let mut altered = r.clone();
    altered["version"] = json!("0.6.0");
    assert_eq!(
        value(&["stability", "verify"], &altered.to_string())["verified"],
        true
    );
    altered["version"] = json!("99.0.0");
    assert_eq!(
        value(&["stability", "verify"], &altered.to_string())["verified"],
        true
    );
    altered["schema"] = json!("polytool.stability.v2");
    assert_eq!(
        run(&["stability", "verify"], &altered.to_string())
            .status
            .code(),
        Some(2)
    );
    for version in [Value::Null, json!(1), json!("")] {
        let mut malformed = r.clone();
        malformed["version"] = version;
        assert_eq!(
            run(&["stability", "verify"], &malformed.to_string())
                .status
                .code(),
            Some(2)
        );
    }
    let mut altered = r.clone();
    altered["items"][0]["counts"]["inside"] = json!(99);
    let output = run(&["stability", "verify"], &altered.to_string());
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        serde_json::from_slice::<Value>(&output.stdout).unwrap()["verified"],
        false
    );
    let mut altered = r;
    altered["items"][0]["certificate"]["real_part"] = json!(["3"]);
    assert_eq!(
        run(&["stability", "verify"], &altered.to_string())
            .status
            .code(),
        Some(1)
    );
}

#[test]
fn invalid_inputs_are_atomic_usage_errors() {
    for args in [
        vec!["stability", "t", "--region", "upper", "--mobius", "1,0,0,1"],
        vec!["stability", "t", "--disk-center", "I"],
        vec!["stability", "t", "--disk-center", "I", "--radius", "-1"],
        vec!["stability", "t", "--mobius", "1,1,1,1"],
        vec!["stability", "t^2001"],
        vec!["stability", "t", "--max-degree", "2001"],
        vec!["stability", "8", "--max-coefficient-bits", "3"],
        vec!["stability", "t", "--max-coefficient-bits", "1000001"],
        vec!["stability", "x+y"],
        vec!["stability", "1+i"],
        vec!["stability", "1+j"],
        vec!["stability", "1.2"],
        vec!["stability", "1/(t+1)"],
        vec!["stability", "t+1", "1/0"],
    ] {
        let output = run(&args, "");
        assert_eq!(
            output.status.code(),
            Some(2),
            "{args:?}: {}",
            String::from_utf8_lossy(&output.stdout)
        );
        assert!(output.stdout.is_empty());
    }
}

#[test]
fn input_files_and_certificate_files() {
    let dir = std::env::temp_dir().join(format!("stability-cli-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let input = dir.join("input.txt");
    let certificate = dir.join("certificate.json");
    std::fs::write(&input, "t+1\n0\n").unwrap();
    let r = value(
        &[
            "stability",
            "--input",
            input.to_str().unwrap(),
            "--certificate",
            certificate.to_str().unwrap(),
            "--region",
            "left",
            "--json",
        ],
        "",
    );
    assert_eq!(r["items"].as_array().unwrap().len(), 2);
    assert_eq!(
        value(
            &[
                "stability",
                "verify",
                "--input",
                certificate.to_str().unwrap()
            ],
            ""
        )["verified"],
        true
    );
    std::fs::remove_dir_all(dir).unwrap();
}
