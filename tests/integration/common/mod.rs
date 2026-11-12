//! Shared CLI runners, output readers and exact-rational constructors for the
//! single integration-test binary. Child-process settings never change the
//! test process's environment or working directory.

use num_bigint::BigInt;
use polytool::BigRational;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Output, Stdio};
use std::time::{Duration, Instant};

/// A command for the tested CLI binary, ready for per-invocation settings.
pub fn polytool_command() -> Command {
    Command::new(env!("CARGO_BIN_EXE_polytool"))
}

fn spawn_with_input(mut command: Command, input: &[u8]) -> Child {
    let mut child = command
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn polytool");
    // Help and argument errors can close the pipe before reading the input.
    if let Err(error) = child.stdin.take().expect("polytool stdin").write_all(input) {
        assert_eq!(error.kind(), std::io::ErrorKind::BrokenPipe);
    }
    child
}

/// Run `polytool args...` with text on standard input.
pub fn run_polytool(args: &[&str], input: &str) -> Output {
    run_polytool_bytes(args, input.as_bytes())
}

/// [`run_polytool`] with raw bytes on standard input.
pub fn run_polytool_bytes(args: &[&str], input: &[u8]) -> Output {
    let mut command = polytool_command();
    command.args(args);
    spawn_with_input(command, input)
        .wait_with_output()
        .expect("wait for polytool")
}

/// [`run_polytool`] for `polytool command args...`.
pub fn run_polytool_subcommand(command: &str, args: &[&str], input: &str) -> Output {
    let mut all_args = vec![command];
    all_args.extend_from_slice(args);
    run_polytool(&all_args, input)
}

/// Run `polytool args...` without standard input.
pub fn run_polytool_without_input(args: &[&str]) -> Output {
    polytool_command()
        .args(args)
        .output()
        .expect("spawn polytool")
}

/// Run a catalog command without inheriting an external OEIS data directory.
pub fn run_polytool_catalog(args: &[&str]) -> Output {
    polytool_command()
        .args(args)
        .env_remove("OEIS_DATA_DIR")
        .output()
        .expect("spawn polytool")
}

/// Run the real-rootedness benchmark CLI.
pub fn run_real_rootedness_benchmark(args: &[&str]) -> Output {
    polytool_command()
        .arg("bench")
        .arg("real-rootedness")
        .args(args)
        .output()
        .expect("run polytool bench real-rootedness")
}

/// Run the real-rooted integer-polynomial enumeration CLI.
pub fn run_real_rooted_enumeration(args: &[&str]) -> Output {
    polytool_command()
        .arg("enumerate-real-rooted")
        .args(args)
        .output()
        .expect("run polytool")
}

/// Directory of frozen CLI outputs and help snapshots.
pub fn snapshot_directory() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/snapshots/cli")
}

/// Directory of inputs used by the frozen CLI invocations.
pub fn snapshot_inputs() -> PathBuf {
    snapshot_directory().join("inputs")
}

fn snapshot_command(args: &[&str]) -> Command {
    let mut command = polytool_command();
    command
        .args(args)
        .current_dir(snapshot_inputs())
        .env_remove("POLYTOOL_DEPRECATION_WARNINGS");
    command
}

/// Run a frozen invocation in its fixture directory.
pub fn run_polytool_snapshot(args: &[&str], input: &str) -> Output {
    spawn_with_input(snapshot_command(args), input.as_bytes())
        .wait_with_output()
        .unwrap()
}

/// Run a help/alias invocation with its existing 60-second deadline.
pub fn run_polytool_snapshot_with_deadline(args: &[&str], input: &str) -> Output {
    let mut child = spawn_with_input(snapshot_command(args), input.as_bytes());
    // Drain both pipes while waiting, so large output cannot block the child.
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

/// Standard output as text.
pub fn stdout(output: &Output) -> String {
    String::from_utf8(output.stdout.clone()).expect("stdout is utf8")
}

/// Standard output as text, asserting successful completion first.
pub fn successful_stdout(output: &Output) -> String {
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    stdout(output)
}

/// Standard error as text.
pub fn stderr(output: &Output) -> String {
    String::from_utf8(output.stderr.clone()).expect("stderr is utf8")
}

/// Standard output parsed as JSON.
pub fn json(output: &Output) -> serde_json::Value {
    serde_json::from_slice(&output.stdout).expect("stdout is JSON")
}

/// The integer `value` as an exact rational.
pub fn q(value: i64) -> BigRational {
    BigRational::from_integer(BigInt::from(value))
}

/// A rational string parsed with the recurrence input convention.
pub fn q_from_str(value: &str) -> BigRational {
    polytool::recurrence::parse_rational_coeff(value).unwrap()
}

/// Integer rows as exact rationals.
pub fn rows(values: &[&[i64]]) -> Vec<Vec<BigRational>> {
    values
        .iter()
        .map(|row| row.iter().map(|&value| q(value)).collect())
        .collect()
}
