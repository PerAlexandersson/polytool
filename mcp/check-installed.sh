#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: ./mcp/check-installed.sh [BINARY]

Compare an installed polytool-mcp version/build commit with this checkout.
Default BINARY: $BIN_DIR/polytool-mcp, or ${PREFIX:-$HOME/.local}/bin/polytool-mcp.
Exit 0: current; 1: stale, missing or unusable binary; 2: unknown checkout commit.
This check never builds or installs anything. Uncommitted edits are not detected.
EOF
}

if [[ ${1:-} == --help || ${1:-} == -h ]]; then
  usage
  exit 0
fi
if [[ $# -gt 1 || ${1:-} == -* ]]; then
  usage >&2
  exit 2
fi

script_dir="$(CDPATH='' cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
binary="${1:-${BIN_DIR:-${PREFIX:-$HOME/.local}/bin}/polytool-mcp}"
version=$(sed -n 's/^version = "\([^"]*\)"$/\1/p' "$script_dir/Cargo.toml")
commit=unknown
# Match the library build script: only trust Git if the manifest is tracked.
if git -C "$script_dir" ls-files --error-unmatch Cargo.toml >/dev/null 2>&1; then
  revision=$(git -C "$script_dir" log -1 --format=%H -- "$script_dir/.." 2>/dev/null || true)
  if [[ $revision =~ ^[0-9a-fA-F]{12,}$ ]]; then
    commit=$(printf '%.12s' "$revision" | tr '[:upper:]' '[:lower:]')
  fi
fi
if [[ $commit == unknown ]]; then
  printf 'cannot verify: checkout build commit is unknown (Git metadata required)\n' >&2
  exit 2
fi

expected="polytool-mcp $version (git $commit)"
if [[ ! -x $binary ]]; then
  printf 'stale: rebuild with install.sh (binary missing or not executable: %s)\n' "$binary" >&2
  exit 1
fi
if ! installed=$("$binary" --version 2>/dev/null); then
  printf 'stale: rebuild with install.sh (cannot read binary version: %s)\n' "$binary" >&2
  exit 1
fi
if [[ $installed != "$expected" ]]; then
  printf 'stale: rebuild with install.sh\n  installed: %s\n  checkout:  %s\n' "$installed" "$expected" >&2
  exit 1
fi
printf 'current: %s\n' "$installed"
