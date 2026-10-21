# Polytool Handoff

## Verified: complete recurrence spaces — 2026-10-21

The isolated recurrence-representations lane adds `recurrence --space` and
`--verify-space`, shared library/MCP/WASM APIs, and a **Complete spaces** browser
view. Reports contain full bounded kernels, successor projections, source-
conditional certificates and finite-row absence duals. Original indices,
source consistency, singular ranks and resource shortages remain explicit.
Legacy finite-fit APIs are unchanged. Exact scope and all passing checks are
in the parent handoff; ownership is released. No installed MCP replacement,
shared-checkout edit, merge, subtree publication or live deployment occurred.
RealRooted #1978 tracks the remaining Lean checker/tactic work. Spectral and
unrestricted-order algorithms remain future stages.

## Verified recurrence representation checkpoint — 2026-10-20

The isolated `feat/recurrence-representations-20261020` lane adds the library,
CLI `--representations`, MCP tool and browser comparison view. Source ownership
is released after the checkpoint; exact scope/checks are in the root handoff.
The shared checkout and installed MCP remain untouched. All outputs are finite
exact fits with explicit domains and coverage; complete-kernel/all-rank source
certification remains future work. The optimized browser bundle passes all
48 Node checks. Website staging is separate; no live deployment or merge.

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
