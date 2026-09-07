use std::io::Write;
use std::process::{Command, Output, Stdio};

fn run_polytool(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_polytool"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn polytool");

    child
        .stdin
        .as_mut()
        .expect("polytool stdin")
        .write_all(input.as_bytes())
        .expect("write polytool stdin");
    child.wait_with_output().expect("wait for polytool")
}

#[test]
fn version_flag_prints_the_exact_build_version() {
    let output = run_polytool(&["--version"], "");

    assert!(output.status.success());
    assert!(output.stderr.is_empty());
    assert_eq!(
        String::from_utf8(output.stdout).expect("version output is UTF-8"),
        format!("{}\n", polytool::version::build_version())
    );
}
