# Polytool contributor guide

## Standalone boundary

- `polytool/` is exported as the standalone public `main` branch.  Its manifest,
  tests, examples, build script, and published sources must not use sibling path
  dependencies or import workspace-only crates.  Put cross-crate research drivers
  in `../experiments/` instead.
- Before changing release-boundary files, validate the standalone subtree with
  `cargo metadata --locked --manifest-path polytool/Cargo.toml --format-version 1
  --no-deps`.  Keep build output outside this Dropbox checkout.

## Mathematical and API conventions

- Preserve exact arithmetic.  Do not replace exact PRS, Bézout, rational, or
  `BigInt` routines with floating-point heuristics.  Use the documented `BigInt`
  APIs when `i64` coefficients can overflow.
- Coefficient vectors are in ascending degree order: entry `i` is the coefficient
  of `t^i`.  Preserve the documented directed-degree and zero-polynomial
  conventions, and add focused tests whenever behavior changes.

## CLI and documentation

- Keep CLI help detailed and AI-friendly.  Machine-readable commands must retain
  stable JSON schemas; serialize exact integers as strings where the JSON number
  model could lose precision.
- Update `README.md` and focused tests in the same change when public API, CLI
  input/output, exactness guarantees, or standalone build behavior changes.
