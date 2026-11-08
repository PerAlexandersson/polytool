//! Frozen CLI contract: invocations used by external scripts keep their exact
//! stdout and exit status, and the help pages match their reviewed snapshots.
//!
//! Golden files live in `tests/snapshots/cli/`:
//!
//! - `frozen/<name>.out` is the stdout of one frozen invocation. These files are
//!   a compatibility promise; re-bless them (`POLYTOOL_BLESS_CLI_FROZEN=1`) only
//!   for an intentional, documented output change.
//! - `help/<page>.txt` is a help page. Re-bless them with
//!   `POLYTOOL_BLESS_CLI_HELP=1` after reviewing a deliberate help change.
//!
//! Standard error is not compared: renamed spellings may add a deprecation
//! warning there.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

fn snapshot_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/snapshots/cli")
}

fn run(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_polytool"))
        .args(args)
        .current_dir(snapshot_dir().join("inputs"))
        .env_remove("POLYTOOL_DEPRECATION_WARNINGS")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    if let Err(error) = child.stdin.take().unwrap().write_all(input.as_bytes()) {
        assert_eq!(error.kind(), std::io::ErrorKind::BrokenPipe);
    }
    child.wait_with_output().unwrap()
}

/// Standard input of a frozen invocation.
enum Stdin {
    None,
    /// A file in `tests/snapshots/cli/inputs/`.
    File(&'static str),
    /// The golden stdout of an earlier frozen invocation.
    Golden(&'static str),
}

struct Frozen {
    name: &'static str,
    args: &'static [&'static str],
    stdin: Stdin,
    exit: i32,
}

const fn frozen(
    name: &'static str,
    args: &'static [&'static str],
    stdin: Stdin,
    exit: i32,
) -> Frozen {
    Frozen {
        name,
        args,
        stdin,
        exit,
    }
}

/// The invocations listed as frozen in the CLI usability plan.
const FROZEN: &[Frozen] = &[
    frozen(
        "real-rooted",
        &["real-rooted"],
        Stdin::File("mixed-roots.txt"),
        0,
    ),
    frozen(
        "real-rooted-json",
        &["real-rooted", "--json"],
        Stdin::File("mixed-roots.txt"),
        0,
    ),
    frozen(
        "interlacing",
        &["interlacing"],
        Stdin::File("eulerian.txt"),
        0,
    ),
    frozen(
        "interlacing-json",
        &["interlacing", "--json"],
        Stdin::File("eulerian.txt"),
        0,
    ),
    frozen(
        "interlacing-profile",
        &["interlacing-profile"],
        Stdin::File("eulerian.txt"),
        0,
    ),
    frozen(
        "interlacing-profile-json",
        &["interlacing-profile", "--json"],
        Stdin::File("eulerian.txt"),
        0,
    ),
    frozen(
        "properties",
        &["properties"],
        Stdin::File("mixed-roots.txt"),
        0,
    ),
    frozen(
        "properties-json",
        &["properties", "--json"],
        Stdin::File("mixed-roots.txt"),
        0,
    ),
    frozen(
        "family-check-json",
        &["family-check", "--json"],
        Stdin::File("eulerian.txt"),
        0,
    ),
    frozen(
        "family-check-require-pass",
        &[
            "family-check",
            "--require-real-rooted",
            "--require-weak-interlacing",
        ],
        Stdin::File("eulerian.txt"),
        0,
    ),
    frozen(
        "family-check-require-fail",
        &[
            "family-check",
            "--require-real-rooted",
            "--require-weak-interlacing",
        ],
        Stdin::File("not-interlacing.txt"),
        1,
    ),
    frozen(
        "recurrence-search-text",
        &[
            "recurrence",
            "--max-rec-len",
            "2",
            "--max-var-deg",
            "2",
            "--max-idx-deg",
            "2",
            "--skip-prefix",
            "1",
            "--denominator",
            "--fit-extra-rows",
            "1",
        ],
        Stdin::File("eulerian.txt"),
        0,
    ),
    frozen(
        "recurrence-search-json",
        &[
            "recurrence",
            "--json",
            "--min-rec-len",
            "1",
            "--max-rec-len",
            "1",
            "--max-var-deg",
            "2",
            "--max-idx-deg",
            "1",
        ],
        Stdin::File("eulerian.txt"),
        0,
    ),
    frozen(
        "recurrence-search-not-found",
        &[
            "recurrence",
            "--json",
            "--max-rec-len",
            "1",
            "--max-var-deg",
            "1",
            "--max-idx-deg",
            "1",
        ],
        Stdin::File("eulerian.txt"),
        0,
    ),
    frozen(
        "recurrence-representations-json",
        &[
            "recurrence",
            "--representations",
            "--json",
            "--max-rec-len",
            "1",
            "--max-var-deg",
            "2",
            "--max-idx-deg",
            "1",
        ],
        Stdin::File("eulerian.txt"),
        0,
    ),
    frozen(
        "recurrence-space",
        &["recurrence", "--space"],
        Stdin::File("spaces-request.json"),
        0,
    ),
    frozen(
        "recurrence-verify-space",
        &["recurrence", "--verify-space"],
        Stdin::Golden("recurrence-space"),
        0,
    ),
    frozen(
        "recurrence-generate-rows",
        &[
            "recurrence-generate",
            "--recurrence",
            "eulerian-recurrence.json",
            "--rows",
            "6",
        ],
        Stdin::None,
        0,
    ),
    frozen(
        "recurrence-generate-stdin-json",
        &[
            "recurrence-generate",
            "--recurrence",
            "-",
            "--rows",
            "5",
            "--json",
        ],
        Stdin::File("eulerian-recurrence.json"),
        0,
    ),
    frozen(
        "oeis-generate-rows",
        &[
            "oeis", "generate", "A008292", "--rows", "5", "--format", "rows",
        ],
        Stdin::None,
        0,
    ),
    frozen(
        "oeis-generate-json",
        &[
            "oeis", "generate", "A008292", "--rows", "4", "--format", "json",
        ],
        Stdin::None,
        0,
    ),
    frozen("sequence", &["sequence", "eulerian", "5"], Stdin::None, 0),
    frozen(
        "sequence-rows",
        &["sequence", "eulerian", "5", "--format", "rows"],
        Stdin::None,
        0,
    ),
    frozen(
        "pf-pencil",
        &["pf-pencil", "--case", "A036969", "--degree", "3"],
        Stdin::None,
        0,
    ),
    frozen(
        "pf-pencil-json",
        &[
            "pf-pencil",
            "--case",
            "A036969",
            "--min-degree",
            "1",
            "--max-degree",
            "2",
            "--json",
        ],
        Stdin::None,
        0,
    ),
    frozen(
        "gamma-expansion-json",
        &["gamma-expansion", "--json"],
        Stdin::File("eulerian.txt"),
        0,
    ),
    frozen("gamma-alias", &["gamma"], Stdin::File("eulerian.txt"), 0),
    frozen(
        "bernstein-expansion-json",
        &["bernstein-expansion", "--degree", "8", "--json"],
        Stdin::File("eulerian.txt"),
        0,
    ),
    frozen(
        "bernstein-alias",
        &["bernstein", "--degree", "8"],
        Stdin::File("eulerian.txt"),
        0,
    ),
];

fn bless(variable: &str) -> bool {
    std::env::var_os(variable).is_some_and(|value| value == "1")
}

fn compare_or_bless(path: &Path, actual: &str, variable: &str, failures: &mut Vec<String>) {
    if bless(variable) {
        std::fs::write(path, actual).unwrap();
        return;
    }
    match std::fs::read_to_string(path) {
        Ok(expected) if expected == actual => {}
        Ok(expected) => failures.push(format!(
            "{} differs ({} vs {} bytes); first differing line:\n  expected: {:?}\n  actual:   {:?}",
            path.display(),
            expected.len(),
            actual.len(),
            expected
                .lines()
                .zip(actual.lines())
                .find(|(a, b)| a != b)
                .map_or("<length>", |(a, _)| a),
            expected
                .lines()
                .zip(actual.lines())
                .find(|(a, b)| a != b)
                .map_or("<length>", |(_, b)| b),
        )),
        Err(_) => failures.push(format!(
            "{} is missing; set {variable}=1 to create it",
            path.display()
        )),
    }
}

/// JSON keys added next to an old key under the additive-JSON policy.
const ADDITIVE_KEYS: &[&str] = &["first_index"];

/// `actual` is the golden JSON document plus keys from [`ADDITIVE_KEYS`].
fn extends(golden: &serde_json::Value, actual: &serde_json::Value) -> bool {
    use serde_json::Value;
    match (golden, actual) {
        (Value::Object(old), Value::Object(new)) => {
            old.iter()
                .all(|(key, value)| new.get(key).is_some_and(|other| extends(value, other)))
                && new
                    .keys()
                    .all(|key| old.contains_key(key) || ADDITIVE_KEYS.contains(&key.as_str()))
        }
        (Value::Array(old), Value::Array(new)) => {
            old.len() == new.len() && old.iter().zip(new).all(|(a, b)| extends(a, b))
        }
        _ => golden == actual,
    }
}

/// Whether a JSON report differs from its golden file only by additive keys.
fn additive_json_extension(golden: &Path, actual: &str) -> bool {
    let Ok(expected) = std::fs::read_to_string(golden) else {
        return false;
    };
    match (
        serde_json::from_str::<serde_json::Value>(&expected),
        serde_json::from_str::<serde_json::Value>(actual),
    ) {
        (Ok(expected), Ok(actual)) => expected != actual && extends(&expected, &actual),
        _ => false,
    }
}

#[test]
fn frozen_invocations_keep_stdout_and_exit_status() {
    let dir = snapshot_dir();
    let mut failures = Vec::new();
    for case in FROZEN {
        let input = match case.stdin {
            Stdin::None => String::new(),
            Stdin::File(name) => std::fs::read_to_string(dir.join("inputs").join(name)).unwrap(),
            Stdin::Golden(name) => {
                std::fs::read_to_string(dir.join("frozen").join(format!("{name}.out"))).unwrap()
            }
        };
        let output = run(case.args, &input);
        assert_eq!(
            output.status.code(),
            Some(case.exit),
            "{}: {:?}\n{}",
            case.name,
            case.args,
            String::from_utf8_lossy(&output.stderr)
        );
        let stdout = String::from_utf8(output.stdout).unwrap();
        let golden = dir.join("frozen").join(format!("{}.out", case.name));
        if !bless("POLYTOOL_BLESS_CLI_FROZEN") && additive_json_extension(&golden, &stdout) {
            continue;
        }
        compare_or_bless(&golden, &stdout, "POLYTOOL_BLESS_CLI_FROZEN", &mut failures);
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn version_first_line_is_stable() {
    let output = run(&["--version"], "");
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let first = stdout.lines().next().unwrap();
    let prefix = format!("polytool {} (", env!("CARGO_PKG_VERSION"));
    assert!(first.starts_with(&prefix), "{first}");
    assert!(first.ends_with(')'), "{first}");
    assert_eq!(run(&["-V"], "").stdout, stdout.as_bytes());
}

#[test]
fn scout_keeps_its_json_report_and_exit_status() {
    // Floating-point values may differ across platforms, so only the shape is frozen.
    let output = run(
        &[
            "bkw-scout",
            "--symbol",
            "1; -x; 1",
            "--grid",
            "5",
            "--top",
            "2",
            "--json",
        ],
        "",
    );
    assert!(output.status.success());
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    for key in ["symbol", "box", "grid", "candidates"] {
        assert!(report.get(key).is_some(), "{key}");
    }
}

/// Every help page reachable from the command index, as `(file stem, args)`.
fn help_pages() -> Vec<(String, Vec<String>)> {
    let output = run(&["help", "--json"], "");
    assert!(output.status.success());
    let index: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let mut pages = vec![("polytool".to_string(), vec!["--help".to_string()])];
    for command in index["commands"].as_array().unwrap() {
        let name = command["name"].as_str().unwrap().to_string();
        pages.push((name.clone(), vec![name.clone(), "--help".to_string()]));
        for sub in command["subcommands"].as_array().into_iter().flatten() {
            let sub = sub["name"].as_str().unwrap().to_string();
            pages.push((
                format!("{name}-{sub}"),
                vec![name.clone(), sub, "--help".to_string()],
            ));
        }
    }
    for topic in index["topics"].as_array().into_iter().flatten() {
        let topic = topic["name"].as_str().unwrap().to_string();
        pages.push((
            format!("topic-{topic}"),
            vec!["help".to_string(), topic.clone()],
        ));
    }
    pages
}

#[test]
fn help_pages_match_snapshots() {
    let dir = snapshot_dir().join("help");
    let mut failures = Vec::new();
    let pages = help_pages();
    // The build commit changes with every commit; keep it out of the snapshots.
    let version = String::from_utf8(run(&["--version"], "").stdout).unwrap();
    let version = version.lines().next().unwrap().to_string();
    for (stem, args) in &pages {
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        let output = run(&args, "");
        assert!(output.status.success(), "{args:?}");
        let stdout = String::from_utf8(output.stdout)
            .unwrap()
            .replace(&version, "polytool <version>");
        compare_or_bless(
            &dir.join(format!("{stem}.txt")),
            &stdout,
            "POLYTOOL_BLESS_CLI_HELP",
            &mut failures,
        );
    }
    let index = run(&["help", "--json"], "");
    compare_or_bless(
        &dir.join("help-json.json"),
        &String::from_utf8(index.stdout).unwrap(),
        "POLYTOOL_BLESS_CLI_HELP",
        &mut failures,
    );
    if !bless("POLYTOOL_BLESS_CLI_HELP") {
        let expected: std::collections::BTreeSet<String> = pages
            .iter()
            .map(|(stem, _)| format!("{stem}.txt"))
            .collect();
        for entry in std::fs::read_dir(&dir).unwrap() {
            let name = entry.unwrap().file_name().into_string().unwrap();
            if name.ends_with(".txt") && !expected.contains(&name) {
                failures.push(format!("stale help snapshot {name}"));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
