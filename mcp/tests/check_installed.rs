//! Installed-binary freshness checks without builds, installation or Git writes.

#![cfg(unix)]

use std::os::unix::fs::PermissionsExt;
use std::path::PathBuf;
use std::process::{Command, Output};

struct Fixture(PathBuf);

impl Fixture {
    fn new() -> Self {
        let directory = std::env::temp_dir().join(format!(
            "polytool-installed-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir(&directory).unwrap();
        std::fs::write(
            directory.join("Cargo.toml"),
            format!("[package]\nversion = \"{}\"\n", env!("CARGO_PKG_VERSION")),
        )
        .unwrap();
        std::fs::write(
            directory.join("check-installed.sh"),
            include_str!("../check-installed.sh"),
        )
        .unwrap();
        let fixture = Self(directory);
        // Mock only the read-only Git queries used by the check script.
        fixture.executable("git", "#!/bin/sh\nprintf '%s\\n' abcdef0123456789\n");
        fixture
    }

    fn executable(&self, name: &str, content: &str) {
        let path = self.0.join(name);
        std::fs::write(&path, content).unwrap();
        std::fs::set_permissions(path, std::fs::Permissions::from_mode(0o755)).unwrap();
    }

    fn binary_version(&self, version: &str) {
        self.executable(
            "polytool-mcp",
            &format!("#!/bin/sh\nprintf '%s\\n' '{version}'\n"),
        );
    }

    fn check(&self, explicit_binary: bool) -> Output {
        let mut command = Command::new("bash");
        command
            .arg(self.0.join("check-installed.sh"))
            .env("BIN_DIR", &self.0)
            .env(
                "PATH",
                std::env::join_paths(
                    std::iter::once(self.0.clone())
                        .chain(std::env::split_paths(&std::env::var_os("PATH").unwrap())),
                )
                .unwrap(),
            );
        if explicit_binary {
            command.arg(self.0.join("polytool-mcp"));
        }
        command.output().unwrap()
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).unwrap();
    }
}

#[test]
fn checks_version_and_commit_and_rejects_unverifiable_builds() {
    let fixture = Fixture::new();
    let version = env!("CARGO_PKG_VERSION");
    fixture.binary_version(&format!("polytool-mcp {version} (git abcdef012345)"));
    for explicit_binary in [false, true] {
        let output = fixture.check(explicit_binary);
        assert!(output.status.success(), "{output:?}");
        assert!(String::from_utf8(output.stdout)
            .unwrap()
            .starts_with("current:"));
    }
    for stale in [
        format!("polytool-mcp {version} (git 000000000000)"),
        "polytool-mcp 0.7.0 (git abcdef012345)".into(),
        "polytool-mcp 0.7.0".into(),
        format!("polytool-mcp {version} (git unknown)"),
    ] {
        fixture.binary_version(&stale);
        let output = fixture.check(true);
        assert_eq!(output.status.code(), Some(1));
        assert!(String::from_utf8(output.stderr)
            .unwrap()
            .contains("stale: rebuild with install.sh"));
    }
    fixture.executable("polytool-mcp", "#!/bin/sh\nexit 1\n");
    assert_eq!(fixture.check(true).status.code(), Some(1));
    std::fs::remove_file(fixture.0.join("polytool-mcp")).unwrap();
    assert_eq!(fixture.check(true).status.code(), Some(1));
    fixture.executable("git", "#!/bin/sh\nexit 1\n");
    let output = fixture.check(true);
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8(output.stderr)
        .unwrap()
        .contains("unknown"));
}

#[test]
fn unrelated_repository_head_does_not_make_subtree_build_stale() {
    let fixture = Fixture::new();
    fixture.executable(
        "git",
        r#"#!/bin/sh
case "$*" in
  *ls-files*) exit 0;;
  *'log -1 --format=%H -- '*) printf '%s\n' abcdef0123456789;;
  *) printf '%s\n' 9999999999999999;;
esac
"#,
    );
    fixture.binary_version(&format!(
        "polytool-mcp {} (git abcdef012345)",
        env!("CARGO_PKG_VERSION")
    ));
    assert!(fixture.check(true).status.success());
}
