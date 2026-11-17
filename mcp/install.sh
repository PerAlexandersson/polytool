#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Install the polytool MCP server.

Usage:
  ./mcp/install.sh [options]

Options:
  --prefix DIR      Install under DIR. Default: $PREFIX or $HOME/.local
  --bin-dir DIR     Install binary directly into DIR. Default: PREFIX/bin
  --build-dir DIR   Cargo target directory. Default: external cache (see below)
  --cargo PATH      Cargo executable. Default: $CARGO or cargo
  --skip-build      Copy an already-built binary from BUILD_DIR/release
  -h, --help        Show this help

Environment:
  PREFIX            Default install prefix
  BIN_DIR           Default binary directory
  CARGO             Default cargo executable
  CARGO_TARGET_DIR  Default build directory, otherwise XDG_CACHE_HOME/polytool/target
  XDG_CACHE_HOME    Cache root. Default: $HOME/.cache
  RUSTC/RUSTDOC     Passed through to cargo when set

Examples:
  ./mcp/install.sh
  ./mcp/install.sh --prefix "$HOME/.local"
  BIN_DIR="$HOME/bin" ./mcp/install.sh
EOF
}

die() {
  printf 'error: %s\n' "$*" >&2
  exit 1
}

json_escape() {
  printf '%s' "$1" | sed 's/\\/\\\\/g; s/"/\\"/g'
}

script_dir="$(CDPATH='' cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

prefix="${PREFIX:-${HOME:-}/.local}"
if [[ -z "${BIN_DIR:-}" ]]; then
  bin_dir="$prefix/bin"
  bin_dir_explicit=0
else
  bin_dir="$BIN_DIR"
  bin_dir_explicit=1
fi
build_dir="${CARGO_TARGET_DIR:-${XDG_CACHE_HOME:-$HOME/.cache}/polytool/target}"
cargo_bin="${CARGO:-cargo}"
skip_build=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix)
      [[ $# -ge 2 ]] || die "--prefix needs a directory"
      prefix="$2"
      if [[ "$bin_dir_explicit" -eq 0 ]]; then
        bin_dir="$prefix/bin"
      fi
      shift 2
      ;;
    --bin-dir)
      [[ $# -ge 2 ]] || die "--bin-dir needs a directory"
      bin_dir="$2"
      bin_dir_explicit=1
      shift 2
      ;;
    --build-dir)
      [[ $# -ge 2 ]] || die "--build-dir needs a directory"
      build_dir="$2"
      shift 2
      ;;
    --cargo)
      [[ $# -ge 2 ]] || die "--cargo needs an executable path"
      cargo_bin="$2"
      shift 2
      ;;
    --skip-build)
      skip_build=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "unknown option: $1"
      ;;
  esac
done

[[ -n "$bin_dir" ]] || die "could not determine install bin directory"
[[ -n "$build_dir" ]] || die "could not determine cargo build directory"

manifest="$script_dir/Cargo.toml"
binary="$build_dir/release/polytool-mcp"

if [[ "$skip_build" -eq 0 ]]; then
  command -v "$cargo_bin" >/dev/null 2>&1 || die "cargo not found: $cargo_bin"
  "$cargo_bin" build \
    --release \
    --manifest-path "$manifest" \
    --target-dir "$build_dir" \
    --bin polytool-mcp
fi

[[ -x "$binary" ]] || die "built binary not found or not executable: $binary"

mkdir -p "$bin_dir"
dest="$bin_dir/polytool-mcp"
install -m 0755 "$binary" "$dest"

dest_escaped="$(json_escape "$dest")"

cat <<EOF
Installed polytool-mcp to:
  $dest

MCP client configuration:
{
  "mcpServers": {
    "polytool": {
      "command": "$dest_escaped"
    }
  }
}

If $bin_dir is not on PATH, use the absolute command above in your MCP client.
EOF
