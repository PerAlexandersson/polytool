//! The clap-driven CLI contract: help pages and topics, `help --json`,
//! did-you-mean errors, universal options (`--flag=value`, `--json`/`--format`,
//! `--input`) and the old-versus-new spelling table (identical stdout, a
//! deprecation warning on stderr for the old spelling only).

use serde_json::Value;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
use std::time::{Duration, Instant};

fn inputs() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/snapshots/cli/inputs")
}

/// Run polytool with a deadline (a stuck command fails the test instead of hanging).
fn run(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_polytool"))
        .args(args)
        .current_dir(inputs())
        .env_remove("POLYTOOL_DEPRECATION_WARNINGS")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    if let Err(error) = child.stdin.take().unwrap().write_all(input.as_bytes()) {
        assert_eq!(error.kind(), std::io::ErrorKind::BrokenPipe);
    }
    // Drain both pipes while waiting, so a large output cannot block the child.
    let drain = |mut pipe: Box<dyn std::io::Read + Send>| {
        std::thread::spawn(move || {
            let mut bytes = Vec::new();
            pipe.read_to_end(&mut bytes).unwrap();
            bytes
        })
    };
    let out = drain(Box::new(child.stdout.take().unwrap()));
    let err = drain(Box::new(child.stderr.take().unwrap()));
    let start = Instant::now();
    let status = loop {
        if let Some(status) = child.try_wait().unwrap() {
            break status;
        }
        if start.elapsed() > Duration::from_secs(60) {
            child.kill().unwrap();
            panic!("polytool {args:?} did not finish within 60 s");
        }
        std::thread::sleep(Duration::from_millis(5));
    };
    Output {
        status,
        stdout: out.join().unwrap(),
        stderr: err.join().unwrap(),
    }
}

mod common;
use common::{stderr, stdout};

fn input(name: &str) -> String {
    std::fs::read_to_string(inputs().join(name)).unwrap()
}

fn index() -> Value {
    let output = run(&["help", "--json"], "");
    assert!(output.status.success());
    serde_json::from_slice(&output.stdout).unwrap()
}

/// Every command and subcommand path in the index.
fn command_paths(index: &Value) -> Vec<(Vec<String>, Value)> {
    let mut paths = Vec::new();
    for command in index["commands"].as_array().unwrap() {
        let name = command["name"].as_str().unwrap().to_string();
        paths.push((vec![name.clone()], command.clone()));
        for sub in command["subcommands"].as_array().unwrap() {
            let sub_name = sub["name"].as_str().unwrap().to_string();
            paths.push((vec![name.clone(), sub_name], sub.clone()));
        }
    }
    paths
}

#[test]
fn every_page_has_usage_examples_and_exit_status() {
    let index = index();
    for (path, spec) in command_paths(&index) {
        let words: Vec<&str> = path.iter().map(String::as_str).collect();
        let page = stdout(&run(&[&words[..], &["--help"]].concat(), ""));
        for section in ["Usage:", "Examples:", "Exit status:"] {
            assert!(page.contains(section), "{path:?} lacks {section}\n{page}");
        }
        assert!(!spec["examples"].as_array().unwrap().is_empty(), "{path:?}");
        let via_help = run(&[&["help"][..], &words[..]].concat(), "");
        assert!(via_help.status.success(), "{path:?}");
        assert_eq!(stdout(&via_help), page, "{path:?}");
        assert_eq!(stdout(&run(&[&words[..], &["-h"]].concat(), "")), page);
    }
}

#[test]
fn topics_are_listed_and_non_empty() {
    let index = index();
    let top = stdout(&run(&["--help"], ""));
    let topics = index["topics"].as_array().unwrap();
    for name in [
        "getting-started",
        "input",
        "coefficient-order",
        "indexing",
        "output-formats",
        "json-schemas",
        "exit-codes",
        "budgets",
        "recurrence-json",
        "compatibility",
    ] {
        assert!(topics.iter().any(|topic| topic["name"] == name), "{name}");
        assert!(top.contains(name), "{name} missing from the overview");
        let output = run(&["help", name], "");
        assert!(output.status.success(), "{name}");
        assert!(stdout(&output).lines().count() > 5, "{name}");
    }
    let compatibility = stdout(&run(&["help", "compatibility"], ""));
    assert!(compatibility.contains("--max-rec-len -> --max-recurrence-length"));
    assert!(compatibility.contains("--rows -> --count"));
    let exit_codes = stdout(&run(&["help", "exit-codes"], ""));
    assert!(exit_codes.contains("family-check"));
}

#[test]
fn help_json_lists_flags_aliases_and_topics() {
    let index = index();
    assert_eq!(index["schema"], "polytool.commands.v1");
    let recurrence = index["commands"]
        .as_array()
        .unwrap()
        .iter()
        .find(|command| command["name"] == "recurrence")
        .unwrap();
    let flags = recurrence["flags"].as_array().unwrap();
    let flag = flags
        .iter()
        .find(|flag| flag["name"] == "max-recurrence-length")
        .unwrap();
    assert_eq!(flag["deprecated_aliases"][0], "--max-rec-len");
    assert_eq!(flag["default"], "5");
    assert!(flags.iter().any(|flag| flag["name"] == "json"));
    let gamma = index["commands"]
        .as_array()
        .unwrap()
        .iter()
        .find(|command| command["name"] == "gamma-expansion")
        .unwrap();
    assert_eq!(gamma["aliases"][0], "gamma");
}

#[test]
fn usage_errors_point_to_help_and_suggest_corrections() {
    for (args, expected) in [
        (vec!["real-rooted", "--jsn"], "--json"),
        (vec!["recurrence", "--max-idx"], "--max-index-degree"),
        (vec!["real-roote"], "real-rooted"),
        (
            vec!["recurrence", "generate"],
            "polytool recurrence-generate",
        ),
        (vec!["real-rooted", "--format", "latex"], "possible values"),
    ] {
        let output = run(&args, "");
        assert_eq!(output.status.code(), Some(2), "{args:?}");
        assert!(output.stdout.is_empty(), "{args:?}");
        let message = stderr(&output);
        assert!(message.contains(expected), "{args:?}: {message}");
        assert!(
            message
                .trim_end()
                .lines()
                .last()
                .unwrap()
                .starts_with("Run `polytool"),
            "{args:?}: {message}"
        );
    }
    let output = run(&["help", "recurence"], "");
    assert_eq!(output.status.code(), Some(2));
    assert!(stderr(&output).contains("polytool help recurrence"));
    let output = run(&["help", "oeis", "generat"], "");
    assert_eq!(output.status.code(), Some(2));
    assert!(stderr(&output).contains("polytool help oeis generate"));

    let output = run(&[], "");
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    assert!(stderr(&output).contains("Start here"));
}

#[test]
fn completion_covers_subcommands_flags_and_topics() {
    for shell in ["bash", "zsh", "fish"] {
        let script = stdout(&run(&["completions", shell], ""));
        for word in ["real-rooted", "generate", "max-recurrence-length"] {
            assert!(script.contains(word), "{shell}: {word}");
        }
        if shell != "fish" {
            assert!(script.contains("getting-started"), "{shell}: topics");
        }
    }
}

/// Values for flags and positionals when probing that every documented
/// option is accepted by its handler.
fn probe_args(path: &[String], spec: &Value) -> Vec<Vec<String>> {
    let mut base: Vec<String> = path.to_vec();
    for positional in spec["positionals"].as_array().unwrap() {
        if positional["required"] == true {
            let value = match positional["name"].as_str().unwrap() {
                "NAME" => "eulerian",
                "SHELL" => "bash",
                name if name.contains("A-NUMBER") || name == "ID" => "A008292",
                _ => "3",
            };
            base.push(value.to_string());
        }
    }
    let mut runs = Vec::new();
    for flag in spec["flags"].as_array().unwrap() {
        let mut args = base.clone();
        args.push(flag["long"].as_str().unwrap().to_string());
        let arity = flag["arity"].as_u64().unwrap() as usize;
        let choices = flag["choices"].as_array().unwrap();
        for _ in 0..arity {
            let value = choices
                .first()
                .and_then(Value::as_str)
                .unwrap_or("1")
                .to_string();
            args.push(value);
        }
        runs.push(args);
    }
    runs
}

#[test]
fn every_documented_flag_reaches_its_handler() {
    let index = index();
    for (path, spec) in command_paths(&index) {
        // These maintain catalog files, read external OEIS data or time
        // benchmarks; their flags are covered by their own tests.
        if matches!(path[0].as_str(), "bench")
            || (path[0] == "oeis"
                && matches!(
                    path.get(1).map(String::as_str),
                    Some("add" | "check" | "evidence")
                ))
            || !spec["subcommands"].as_array().unwrap().is_empty()
        {
            continue;
        }
        for args in probe_args(&path, &spec) {
            let words: Vec<&str> = args.iter().map(String::as_str).collect();
            let output = run(&words, "");
            let message = stderr(&output).to_lowercase();
            assert!(
                !message.contains("unexpected argument")
                    && !message.contains("unknown option")
                    && !message.contains(&format!("unknown {} option", path[0])),
                "{words:?} was not accepted: {message}"
            );
        }
    }
}

#[test]
fn universal_options_are_equivalent() {
    let rows = input("mixed-roots.txt");
    let json = stdout(&run(&["real-rooted", "--json"], &rows));
    assert_eq!(
        stdout(&run(&["real-rooted", "--format", "json"], &rows)),
        json
    );
    assert_eq!(stdout(&run(&["real-rooted", "--format=json"], &rows)), json);
    assert_eq!(
        stdout(&run(
            &["real-rooted", "--json", "--input", "mixed-roots.txt"],
            ""
        )),
        json
    );
    assert_eq!(
        stdout(&run(&["real-rooted", "--json", "--input=-"], &rows)),
        json
    );
    let text = stdout(&run(&["real-rooted"], &rows));
    assert_eq!(
        stdout(&run(&["real-rooted", "--json", "--format", "text"], &rows)),
        text
    );
    let euler = input("eulerian.txt");
    let spaced = stdout(&run(
        &["recurrence", "--max-recurrence-length", "2"],
        &euler,
    ));
    assert_eq!(
        stdout(&run(&["recurrence", "--max-recurrence-length=2"], &euler)),
        spaced
    );
    // Values that start with a dash are values, not options.
    let output = run(
        &[
            "bkw-scout",
            "--symbol=1; -x; 1",
            "--grid",
            "5",
            "--top",
            "1",
        ],
        "",
    );
    assert!(output.status.success(), "{}", stderr(&output));
    let missing = run(&["real-rooted", "--input", "no-such-file.txt"], "");
    assert_eq!(missing.status.code(), Some(2));
    assert!(stderr(&missing).contains("no-such-file.txt"));
    // JSON-only commands accept --format json and reject --format text.
    let request = r#"{"polynomials":[["-1","0","1"]]}"#;
    let default = stdout(&run(&["root-enclosures"], request));
    assert_eq!(
        stdout(&run(&["root-enclosures", "--format", "json"], request)),
        default
    );
    assert_eq!(
        run(&["root-enclosures", "--format", "text"], request)
            .status
            .code(),
        Some(2)
    );
}

/// `(old arguments, new arguments, stdin file)`: identical stdout and exit
/// status; the old spelling warns on stderr, the new one does not.
const RENAMES: &[(&[&str], &[&str], &str)] = &[
    (
        &[
            "recurrence",
            "--json",
            "--min-rec-len",
            "1",
            "--max-rec-len",
            "2",
            "--min-var-deg",
            "0",
            "--max-var-deg",
            "2",
            "--min-idx-deg",
            "0",
            "--max-idx-deg",
            "1",
            "--min-diff-deg",
            "0",
            "--max-diff-deg",
            "1",
            "--min-margin",
            "1",
            "--fit-extra-rows",
            "1",
        ],
        &[
            "recurrence",
            "--json",
            "--min-recurrence-length",
            "1",
            "--max-recurrence-length",
            "2",
            "--min-variable-degree",
            "0",
            "--max-variable-degree",
            "2",
            "--min-index-degree",
            "0",
            "--max-index-degree",
            "1",
            "--min-derivative-order",
            "0",
            "--max-derivative-order",
            "1",
            "--min-equation-margin",
            "1",
            "--fit-extra-rows",
            "1",
        ],
        "eulerian.txt",
    ),
    (
        &[
            "recurrence",
            "--json",
            "--max-rec-len",
            "1",
            "--max-var-deg",
            "1",
            "--max-idx-deg",
            "1",
            "--inhomogeneous",
            "--min-inhomo-var-deg",
            "0",
            "--max-inhomo-var-deg",
            "1",
            "--min-inhomo-idx-deg",
            "0",
            "--max-inhomo-idx-deg",
            "1",
            "--alternating-sign",
            "--max-candidates",
            "200",
        ],
        &[
            "recurrence",
            "--json",
            "--max-recurrence-length",
            "1",
            "--max-variable-degree",
            "1",
            "--max-index-degree",
            "1",
            "--try-inhomogeneous",
            "--min-inhomogeneous-variable-degree",
            "0",
            "--max-inhomogeneous-variable-degree",
            "1",
            "--min-inhomogeneous-index-degree",
            "0",
            "--max-inhomogeneous-index-degree",
            "1",
            "--try-alternating-sign",
            "--max-candidates",
            "200",
        ],
        "eulerian.txt",
    ),
    (
        &[
            "recurrence",
            "--json",
            "--max-rec-len",
            "1",
            "--denominator",
            "--max-denom-var-deg",
            "1",
            "--max-denom-idx-deg",
            "1",
            "--no-verify",
            "--max-laurent-x-power",
            "1",
            "--max-candidates",
            "200",
        ],
        &[
            "recurrence",
            "--json",
            "--max-recurrence-length",
            "1",
            "--try-denominator",
            "--max-denominator-variable-degree",
            "1",
            "--max-denominator-index-degree",
            "1",
            "--no-held-out-rows",
            "--max-laurent-power",
            "1",
            "--max-candidates",
            "200",
        ],
        "eulerian.txt",
    ),
    (
        &[
            "recurrence-generate",
            "--recurrence",
            "eulerian-recurrence.json",
            "--rows",
            "5",
        ],
        &[
            "recurrence-generate",
            "--recurrence",
            "eulerian-recurrence.json",
            "--count",
            "5",
        ],
        "",
    ),
    (
        &[
            "recurrence-generate",
            "--recurrence",
            "-",
            "--additional",
            "2",
            "--json",
        ],
        &[
            "recurrence-generate",
            "--recurrence",
            "-",
            "--additional-count",
            "2",
            "--json",
        ],
        "eulerian-recurrence.json",
    ),
    (
        &[
            "recurrence-generate",
            "--recurrence",
            "eulerian-recurrence.json",
            "--rows",
            "6",
            "--format",
            "bfile",
            "--first-row",
            "2",
            "--offset",
            "0",
            "--reversed",
        ],
        &[
            "recurrence-generate",
            "--recurrence",
            "eulerian-recurrence.json",
            "--count",
            "6",
            "--format",
            "bfile",
            "--from-index",
            "2",
            "--first-term-index",
            "0",
            "--descending",
        ],
        "",
    ),
    (
        &[
            "recurrence-degree",
            "--recurrence",
            "eulerian-recurrence.json",
            "--rows",
            "12",
        ],
        &[
            "recurrence-degree",
            "--recurrence",
            "eulerian-recurrence.json",
            "--count",
            "12",
        ],
        "",
    ),
    (
        &[
            "oeis",
            "generate",
            "A008292",
            "--start-row",
            "3",
            "--rows",
            "3",
            "--format",
            "json",
        ],
        &[
            "oeis",
            "generate",
            "A008292",
            "--from-index",
            "3",
            "--count",
            "3",
            "--format",
            "json",
        ],
        "",
    ),
    (
        &["family-check", "--recurrence", "--json"],
        &["family-check", "--find-recurrence", "--json"],
        "eulerian.txt",
    ),
    (
        &[
            "enumerate-real-rooted",
            "--degree",
            "4",
            "--leading-desc",
            "1,0,-6",
        ],
        &[
            "enumerate-real-rooted",
            "--degree",
            "4",
            "--leading-coefficients-descending",
            "1,0,-6",
        ],
        "",
    ),
    (
        &[
            "brenti-recognize",
            "--first-row",
            "1",
            "--max-integer-nodes",
            "1000",
            "--max-span",
            "2",
            "--json",
        ],
        &[
            "brenti-recognize",
            "--first-index",
            "1",
            "--max-integer-search-nodes",
            "1000",
            "--max-target-span",
            "2",
            "--json",
        ],
        "eulerian.txt",
    ),
    (
        &[
            "pf-network",
            "--first-row",
            "1",
            "--max-integer-nodes",
            "1000",
            "--max-length",
            "2",
            "--json",
        ],
        &[
            "pf-network",
            "--first-index",
            "1",
            "--max-integer-search-nodes",
            "1000",
            "--max-step-length",
            "2",
            "--json",
        ],
        "eulerian.txt",
    ),
    (
        &["left-production", "--first-row", "0", "--json"],
        &["left-production", "--first-index", "0", "--json"],
        "eulerian.txt",
    ),
    (
        &["cyclic-sieving-sequence", "--offsets", "0,1", "--json"],
        &[
            "cyclic-sieving-sequence",
            "--order-offsets",
            "0,1",
            "--json",
        ],
        "eulerian.txt",
    ),
    (
        &[
            "bkw-scout",
            "--symbol",
            "1; -x; 1",
            "--grid",
            "5",
            "--grid-re",
            "5",
            "--grid-im",
            "7",
            "--tol",
            "1e-10",
            "--max-iter",
            "50",
            "--min-imag",
            "0.01",
            "--top",
            "2",
        ],
        &[
            "bkw-scout",
            "--symbol",
            "1; -x; 1",
            "--grid",
            "5",
            "--real-grid-points",
            "5",
            "--imaginary-grid-points",
            "7",
            "--tolerance",
            "1e-10",
            "--max-iterations",
            "50",
            "--min-imaginary-part",
            "0.01",
            "--top",
            "2",
        ],
        "",
    ),
];

#[test]
fn old_and_new_spellings_give_identical_output() {
    for (old, new, file) in RENAMES {
        let stdin = if file.is_empty() {
            String::new()
        } else {
            input(file)
        };
        let old_output = run(old, &stdin);
        let new_output = run(new, &stdin);
        assert_eq!(
            old_output.status.code(),
            new_output.status.code(),
            "{new:?}"
        );
        assert!(
            new_output.status.code().is_some_and(|code| code != 2),
            "{new:?}: {}",
            stderr(&new_output)
        );
        assert_eq!(
            stdout(&old_output),
            stdout(&new_output),
            "{old:?} vs {new:?}"
        );
        assert!(
            stderr(&old_output).contains("is deprecated"),
            "{old:?} did not warn: {}",
            stderr(&old_output)
        );
        assert!(
            !stderr(&new_output).contains("is deprecated"),
            "{new:?} warned: {}",
            stderr(&new_output)
        );
    }
}

#[test]
fn deprecation_warnings_can_be_silenced() {
    let output = Command::new(env!("CARGO_BIN_EXE_polytool"))
        .args(["oeis", "generate", "A008292", "--rows", "2"])
        .env("POLYTOOL_DEPRECATION_WARNINGS", "0")
        .output()
        .unwrap();
    assert!(output.status.success());
    assert!(output.stderr.is_empty(), "{}", stderr(&output));
}

#[test]
fn renamed_json_keys_are_emitted_under_both_names() {
    let output = run(
        &[
            "oeis",
            "generate",
            "A008292",
            "--count",
            "2",
            "--from-index",
            "2",
            "--json",
        ],
        "",
    );
    let report: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["first_row"], 2);
    assert_eq!(report["first_index"], report["first_row"]);
    for command in ["brenti-recognize", "pf-network"] {
        let output = run(
            &[command, "--json", "--first-index", "1"],
            &input("eulerian.txt"),
        );
        let report: Value = serde_json::from_slice(&output.stdout).unwrap();
        assert_eq!(report["options"]["first_index"], 1, "{command}");
        assert_eq!(report["options"]["first_row_index"], 1, "{command}");
    }
    let output = run(
        &["left-production", "--json", "--first-index", "0"],
        &input("eulerian.txt"),
    );
    let report: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["first_index"], report["first_row_index"]);
}

#[cfg(unix)]
#[test]
fn getting_started_commands_run() {
    let page = stdout(&run(&["help", "getting-started"], ""));
    let binary = PathBuf::from(env!("CARGO_BIN_EXE_polytool"));
    let path = format!(
        "{}:{}",
        binary.parent().unwrap().display(),
        std::env::var("PATH").unwrap_or_default()
    );
    let dir = std::env::temp_dir().join(format!("polytool-getting-started-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let commands: Vec<&str> = page
        .lines()
        .filter_map(|line| line.strip_prefix("  "))
        .filter(|line| line.starts_with("polytool "))
        .collect();
    assert!(commands.len() >= 6);
    for command in commands {
        let output = Command::new("sh")
            .args(["-c", command])
            .current_dir(&dir)
            .env("PATH", &path)
            .stdin(Stdio::null())
            .output()
            .unwrap();
        assert!(output.status.success(), "{command}: {}", stderr(&output));
    }
    std::fs::remove_dir_all(&dir).unwrap();
}
