//! Build-version formatting shared by the CLI and tests.

/// Length of the abbreviated hexadecimal Git commit in version output.
pub const ABBREVIATED_COMMIT_LENGTH: usize = 12;

/// Git commit captured by `build.rs`, or `None` when trustworthy metadata was
/// unavailable at build time.
pub fn build_git_commit() -> Option<&'static str> {
    let commit = env!("POLYTOOL_GIT_COMMIT");
    (!commit.is_empty()).then_some(commit)
}

fn valid_abbreviated_commit(commit: &str) -> bool {
    commit.len() == ABBREVIATED_COMMIT_LENGTH && commit.bytes().all(|byte| byte.is_ascii_hexdigit())
}

/// Format a deterministic Polytool version line.
///
/// Invalid or missing metadata is reported as unavailable instead of being
/// displayed as a Git revision.
pub fn format_version(version: &str, git_commit: Option<&str>) -> String {
    match git_commit.filter(|commit| valid_abbreviated_commit(commit)) {
        Some(commit) => format!("polytool {version} (git {})", commit.to_ascii_lowercase()),
        None => format!("polytool {version} (git unavailable)"),
    }
}

/// Return the version line for the currently built Polytool binary.
pub fn build_version() -> String {
    format_version(env!("CARGO_PKG_VERSION"), build_git_commit())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formats_normal_build_metadata_deterministically() {
        assert_eq!(
            format_version("1.2.3", Some("ABCDEF012345")),
            "polytool 1.2.3 (git abcdef012345)"
        );
    }

    #[test]
    fn reports_missing_or_invalid_metadata_honestly() {
        assert_eq!(
            format_version("1.2.3", None),
            "polytool 1.2.3 (git unavailable)"
        );
        assert_eq!(
            format_version("1.2.3", Some("not-a-commit")),
            "polytool 1.2.3 (git unavailable)"
        );
    }
}
