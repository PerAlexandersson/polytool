use std::env;
use std::path::Path;
use std::process::Command;

const ABBREVIATED_COMMIT_LENGTH: usize = 12;

fn normalized_commit(value: &str) -> Option<String> {
    let value = value.trim();
    if value.len() < ABBREVIATED_COMMIT_LENGTH
        || !value.bytes().all(|byte| byte.is_ascii_hexdigit())
    {
        return None;
    }
    Some(value[..ABBREVIATED_COMMIT_LENGTH].to_ascii_lowercase())
}

fn git_output(manifest_dir: &Path, args: &[&str]) -> Option<String> {
    let output = Command::new("git")
        .arg("-C")
        .arg(manifest_dir)
        .args(args)
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    String::from_utf8(output.stdout).ok()
}

fn git_commit(manifest_dir: &Path) -> Option<String> {
    git_output(manifest_dir, &["ls-files", "--error-unmatch", "Cargo.toml"])?;
    normalized_commit(&git_output(
        manifest_dir,
        &["rev-parse", "--verify", "HEAD"],
    )?)
}

fn track_git_ref(manifest_dir: &Path) {
    let Some(head_path) = git_output(manifest_dir, &["rev-parse", "--git-path", "HEAD"]) else {
        return;
    };
    println!("cargo:rerun-if-changed={}", head_path.trim());

    if let Some(reference) =
        git_output(manifest_dir, &["symbolic-ref", "-q", "HEAD"]).and_then(|reference| {
            git_output(manifest_dir, &["rev-parse", "--git-path", reference.trim()])
        })
    {
        println!("cargo:rerun-if-changed={}", reference.trim());
    }

    if let Some(packed_refs) = git_output(manifest_dir, &["rev-parse", "--git-path", "packed-refs"])
    {
        println!("cargo:rerun-if-changed={}", packed_refs.trim());
    }
}

fn main() {
    println!("cargo:rerun-if-env-changed=POLYTOOL_GIT_COMMIT");
    let manifest_dir = env::var_os("CARGO_MANIFEST_DIR")
        .map(std::path::PathBuf::from)
        .expect("Cargo sets CARGO_MANIFEST_DIR for build scripts");

    let commit = match env::var("POLYTOOL_GIT_COMMIT") {
        Ok(value) => normalized_commit(&value),
        Err(_) => {
            track_git_ref(&manifest_dir);
            git_commit(&manifest_dir)
        }
    };
    println!(
        "cargo:rustc-env=POLYTOOL_GIT_COMMIT={}",
        commit.as_deref().unwrap_or("")
    );
}
