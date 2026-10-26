# Polytool contributor guide

## Standalone build

- The crate must build standalone.  Its manifest, tests, examples, build script,
  and published sources must not use path dependencies outside this directory
  (the in-tree `crates/`, `mcp/`, and `web/` packages are fine).
- Validate with `cargo metadata --locked --format-version 1 --no-deps`.
- Keep build output (`target/`, `pkg/`, and similar) out of the source tree,
  for example by setting `CARGO_TARGET_DIR`.

## Mathematical and API conventions

- Preserve exact arithmetic.  Do not replace exact PRS, Bézout, rational, or
  `BigInt` routines with floating-point heuristics.  Use the documented `BigInt`
  APIs when `i64` coefficients can overflow, and keep explicit overflow and
  shape contracts.
- Coefficient vectors are in ascending degree order: entry `i` is the coefficient
  of `t^i`.  Preserve the documented directed-degree and zero-polynomial
  conventions, and add focused tests whenever behavior changes.

## CLI and documentation

- Keep CLI help detailed and AI-friendly.  Machine-readable commands must retain
  stable JSON schemas; serialize exact integers as strings where the JSON number
  model could lose precision.
- Update `README.md` and focused tests in the same change when public API, CLI
  input/output, exactness guarantees, or standalone build behavior changes.

## Verification

- Run focused library and CLI tests before broader checks.
- Run `python3 scripts/build_oeis_catalog.py --check-bundled`; it needs no
  external OEIS or Lean sources.  It verifies the committed
  `src/oeis_catalog_generated.rs` against `src/oeis_catalog_generated.sha256`.
- Verify that `cargo package --list` excludes agent notes, CI configuration,
  benchmark output, and generator tooling, while retaining the recurrence
  fixtures used by the installed `bench recurrence-fixtures` command.
- Preserve unrelated changes when editing.
