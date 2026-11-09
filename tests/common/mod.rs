//! Helpers shared by the integration tests: running the `polytool` binary and
//! reading its output, and small exact-rational constructors.
//!
//! Include with `mod common;` and import what the test needs; a test file that
//! uses only some of the helpers is fine (`dead_code` is allowed here).

#![allow(dead_code)]

use num_bigint::BigInt;
use polytool::BigRational;
use std::io::Write;
use std::process::{Command, Output, Stdio};

/// Run `polytool args...` with `input` on standard input and collect its
/// output.  The CLI may exit before reading its input (`--help`, an argument
/// error); a closed pipe is then expected, not a test failure.
pub fn run_polytool(args: &[&str], input: &str) -> Output {
    run_polytool_bytes(args, input.as_bytes())
}

/// [`run_polytool`] with raw bytes on standard input.
pub fn run_polytool_bytes(args: &[&str], input: &[u8]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_polytool"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn polytool");
    if let Err(error) = child.stdin.take().expect("polytool stdin").write_all(input) {
        assert_eq!(error.kind(), std::io::ErrorKind::BrokenPipe);
    }
    child.wait_with_output().expect("wait for polytool")
}

/// [`run_polytool`] for `polytool command args...`.
pub fn run_polytool_subcommand(command: &str, args: &[&str], input: &str) -> Output {
    let mut all_args = vec![command];
    all_args.extend_from_slice(args);
    run_polytool(&all_args, input)
}

/// Run `polytool args...` without standard input.
pub fn run_polytool_without_input(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_polytool"))
        .args(args)
        .output()
        .expect("spawn polytool")
}

/// Standard output as text.
pub fn stdout(output: &Output) -> String {
    String::from_utf8(output.stdout.clone()).expect("stdout is utf8")
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

/// Integer rows as exact rationals.
pub fn rows(values: &[&[i64]]) -> Vec<Vec<BigRational>> {
    values
        .iter()
        .map(|row| row.iter().map(|&value| q(value)).collect())
        .collect()
}
