# Polytool Handoff

Full chronological history through 2026-09-16 is archived in
[`doc/HANDOFF_ARCHIVE_THROUGH_2026-09-16.md`](doc/HANDOFF_ARCHIVE_THROUGH_2026-09-16.md).

## Current state

`master` in the parent Rust repository is canonical. The public `main` branch
is a generated standalone subtree and must be updated only with
`../scripts/sync-polytool-main.sh` after the parent checkpoint is committed,
pushed, clean for eligible paths, and verified.

The repository-health batch owns Polytool manifest/package metadata, the
self-contained generated-OEIS checksum check, and the completed subset-DP
permanent optimization in `src/vec_poly.rs`. Source ownership is released. No
release or deployment has been authorized.

## Verification expectations

- Keep exact arithmetic and explicit overflow/shape contracts.
- Run focused library/CLI tests before broader checks.
- Run `python3 scripts/build_oeis_catalog.py --check-bundled` without requiring
  external OEIS or Lean sources.
- Verify `cargo package --list` excludes internal agent notes, handoffs, CI,
  benchmark output, and generator tooling while retaining the recurrence
  fixtures used by the installed benchmark CLI.
- Keep targets outside the source checkout and preserve unrelated changes.
