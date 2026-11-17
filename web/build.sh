#!/usr/bin/env bash
set -euo pipefail

script_dir="$(CDPATH='' cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# An external cache works in standalone checkouts too. Explicit paths may be
# relative to the caller; normalize them before changing the working directory.
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-${XDG_CACHE_HOME:-$HOME/.cache}/polytool/target}"
[[ $CARGO_TARGET_DIR == /* ]] || export CARGO_TARGET_DIR="$PWD/$CARGO_TARGET_DIR"
output_dir="${POLYTOOL_WASM_DIR:-$CARGO_TARGET_DIR/web-pkg}"
[[ $output_dir == /* ]] || output_dir="$PWD/$output_dir"
wasm_pack="${WASM_PACK:-wasm-pack}"
cd "$script_dir"
"$wasm_pack" build --target no-modules --release --out-dir "$output_dir"
printf '\nBuild complete: %s\n' "$output_dir"
printf 'Stage this bundle as pkg/ beside index.html to serve the application.\n'
