# Polytool handoff

## Active GitHub issues #1 and #2 (2026-09-06)

The coding worker owns the issue implementation in the isolated monorepo
worktree `/tmp/polytool-issues-1-2-20260906` on branch
`fix/polytool-issues-1-2-20260906`, starting at freshly fetched
`origin/master` commit `6486f93`.  No open pull request existed when work
started.

Owned files are `build.rs`, `src/version.rs`, `src/lib.rs`,
`src/recurrence.rs`, `src/bin/polytool.rs`, focused new or existing tests under
`tests/`, `README.md`, `mcp/src/lib.rs`, `mcp/README.md`, this handoff, and the
root `HANDOFF.md`.  Ownership covers only GitHub #1 (`polytool --version` with
honest reproducible Git metadata) and #2 (an exact recurrence-candidate
budget with distinct exhaustion).  Standalone `main` remains a generated
subtree projection and will not be edited directly.

Checkpoint `92f7612` implements the crate-version plus lowercase 12-digit Git
commit line, with build-time ref tracking, an explicit reproducible-build
override, and `(git unavailable)` for missing or invalid metadata.  Checkpoint
`87ef7f8` adds `AdaptiveSearchBudget`, `AdaptiveSearchOutcome`, exact
candidate counting before all filters, `--max-candidates`, JSON/exit-status
termination reporting, and matching MCP options/status/counters.  Existing
unbounded entry points remain unchanged wrappers.

Focused tests cover normal/fallback version formatting and the real CLI, zero
and one-candidate budgets, success on the exact boundary, full search-space
failure at the boundary, distinct exhaustion, MCP parity, and unchanged
unbounded behavior.  Established verification passes: 318 non-OEIS library
tests, all 61 recurrence tests, 19 CLI BigInt tests, 7 CLI OEIS tests, 6 new
CLI tests, 2 overfit fixtures, 5 interlacing API tests, 22 MCP tests plus its
binary/docs targets, 5 Polytool doctests, focused imported-OEIS validation,
strict Polytool/MCP Clippy, Cargo metadata, formatting, and
`git diff --check`.  All Rust commands used external
`CARGO_TARGET_DIR=/cargo-target/ai-projects`, `timeout 60s`, and `nice -n 10`.

The task branch is published as PR #3.  The repository has no GitHub Actions
workflow and unprotected `master` reports no status checks.  A local
projected-subtree preflight at `a7bcd0a` passes all six CLI tests and its
actual version is `polytool 0.2.1-rc.5 (git a7bcd0abf054)`.

On 2026-09-07 the user clarified that the green-CI wording came from
supervisor caution rather than an explicit user constraint and accepted the
complete local monorepo and standalone checks as the merge gate.  The earlier
blocker wording was therefore incorrect.  PR #3 publication, standalone
projection, issue closure, and ownership release are proceeding.

## Final monorepo integration (2026-09-06)

The completed review-fix history was integrated with clean local monorepo
`master` by merge commit `a801eb1`. The merge preserves parent `41d57bd` with
all Polytool fixes and parent `29c075e` with all audited monorepo changes. It
merged without textual conflicts, so no LLT or handoff content was discarded.

Post-merge verification used external
`CARGO_TARGET_DIR=/cargo-target/ai-projects`, `timeout 60s`, and `nice -n 10`:

- the Polytool non-OEIS library suite passed 312 tests;
- focused recurrence, linear-algebra, and parser coverage is included, with 57
  recurrence and 62 linear-algebra tests passing independently;
- CLI BigInt, recurrence-overfit, interlacing API, documentation, and the
  focused imported-OEIS validation replay passed;
- Polytool MCP passed all library, binary, and documentation targets;
- strict Clippy for Polytool and Polytool MCP passed;
- the other changed monorepo packages and tracked experiment binaries passed
  their proportional checks;
- Cargo metadata and `git diff --check` passed.

This handoff is part of the user-authorized final monorepo publication and
standalone `polytool/` subtree projection. Integration ownership is released
after remote-ref verification. Ehrcalc is outside this operation and remains
untouched.

## Completed review fixes (2026-09-06)

The review-fix worker used the regular isolated worktree
`/tmp/polytool-review-fixes-20260906` on branch
`fix/polytool-review-20260906`. Ownership covered `src/parse.rs`, `src/lib.rs`,
`src/real_rootedness.rs`, `src/recurrence.rs`, `src/linalg.rs`,
`src/bin/polytool.rs`, `mcp/src/lib.rs`, `tests/cli_bigint.rs`, focused inline
tests, and this handoff. The implementation is complete and that ownership is
released after the final handoff commit.

History and implementation checkpoints:

- `76e3b0c` merges `origin/master` into local master `d0a1ca3`; its two parents
  preserve both the six local-only and nine remote-only commits;
- `6ca87f6` records initial isolated-worktree ownership;
- `ffa4ad4` implements all confirmed review fixes and both requested lazy
  iteration improvements.

The implementation makes modular recurrence rejection exact by requiring a
full-column-rank modular certificate and otherwise falling back to rational
solving. It checks exponent and recurrence-index arithmetic, bounds dense
parser allocation and CLI/MCP input, caps MCP sequence, recurrence, and finite
Lace workloads, replaces Ehrhart and modular-row panics with structured errors,
and validates ragged total-positivity matrices. Bivariate recurrence
coefficients tolerate the existing public ragged representation without
panicking or dropping terms and are normalized at the JSON boundary. Total
positivity combinations and score-ordered recurrence candidates are now lazy;
equivalence tests compare both iterators with their previous eager order.

Compatibility choices: `BivarPoly::coeffs` remains public, so existing struct
literals continue to compile. Ragged rows are interpreted with missing exact
zeros and JSON round-trips become rectangular. The Ehrhart conversion and
public `SparseModRow` construction/update APIs now return typed `Result`s; this
source-level change is intentional because those APIs previously panicked on
invalid exact input or modulus zero. Existing valid arithmetic and recurrence
ordering are unchanged.

Verification used external `CARGO_TARGET_DIR=/cargo-target/ai-projects` and
60-second, reduced-priority Rust commands:

```text
cargo test -q -p polytool --lib -- --skip oeis                 312 passed
cargo test -q -p polytool --lib recurrence::tests               57 passed
cargo test -q -p polytool --lib linalg::tests                   62 passed
cargo test -q -p polytool --lib parse::tests                    17 passed
cargo test -q -p polytool --lib real_rootedness::tests::test_ehrhart
                                                                  4 passed
cargo test -q -p polytool --test cli_bigint                     19 passed
cargo test -q -p polytool --test recurrence_overfit_fixtures     2 passed
cargo test -q -p polytool --test interlacing_api                 5 passed
cargo test -q -p polytool --doc                                  5 passed
cargo test -q -p polytool-mcp                         21 + 2 + 1 passed
cargo test -q -p polytool --lib \
  oeis::tests::every_imported_lean_definition_reproduces_its_validation_row
                                                                  1 passed
cargo clippy -q -p polytool -p polytool-mcp --all-targets --
  -D warnings -A clippy::manual-is-multiple-of
  -A clippy::needless-range-loop -A clippy::bool-assert-comparison
                                                                  passed
git diff --check                                                 passed
```

Two unchanged exhaustive fixture replays exceeded the required 60-second cap:
`oeis::tests::every_sparse_definition_reproduces_its_fixture_rows` and
`recurrence_json_fixtures_regenerate_raw_rows`. Each was terminated by
`timeout` with status 124 and emitted no failure before termination. Their
focused recurrence paths and the other 312 library tests pass.

That checkpoint itself was not pushed and left the canonical checkout at
`d0a1ca3`; the later final-integration section above supersedes that historical
state while preserving the original verification record.

## OEIS recurrence catalog

The host Codex supervisor completed the OEIS catalog expansion at the user's
request.  No Rust worker is active and catalog file ownership is released.
Unrelated dirty Rust-workspace files remain untouched.  Files changed by the
completed task are:

- `src/oeis.rs` (new);
- `src/lib.rs`;
- `src/bin/polytool.rs`;
- `mcp/src/lib.rs`;
- `tests/cli_oeis.rs` (new);
- `scripts/build_oeis_catalog.py` and generated catalog files (new);
- additive documentation in `README.md` and this `HANDOFF.md`.

The pre-existing coupled-recurrence changes in `src/recurrence.rs` and
`src/linalg.rs` are adopted as dependencies but will not be mixed into the
catalog implementation without a separate verified checkpoint.  The active
`real-rooted-oeis` worker owns its checkout; this task reads recurrence and row
data there and in `projects/OEIS-polynomials` but does not edit either project.

Implementation status on 2026-09-04:

- `polytool::oeis` contains 785 recurrence-backed A-number functions, exact
  sparse recurrence decoding, dynamic lookup, row-range generation, and
  structural-zero restoration for OEIS output;
- 125 holdout-backed entries and 630 OEIS-prefix-validated entries are
  available by default; 30 entries require `--include-experimental`;
- 767 entries have a locally validated OEIS flattened-prefix mapping and
  support strict b-file output; the remaining 18 still support rows, triangles,
  polynomials, JSON, JSONL, and CSV;
- CLI commands are `polytool oeis list`, `polytool oeis info`, and
  `polytool oeis generate`; MCP tools are `list_oeis_sequences`,
  `get_oeis_sequence`, and `generate_oeis_rows`;
- `scripts/build_oeis_catalog.py` imports the 73 machine recurrence fixtures,
  the plain-file OEIS queue, and all 728 generated recurrence definitions from
  `projects/real-rooted-oeis-proofs/ProofsOeis`.  Every Lean definition fits
  the supported canonical grammar (maximum lag 5 and derivative order 2).
  The importer clears rational recurrence denominators, converts index
  conventions exactly, validates generated rows against local OEIS data, emits
  one embedded Rust replay row for each of the 630 new validated entries, and
  supports `--check` drift detection.

The unsafe plain-file queue recurrences remain rejected.  Where the Lean proof
repository contains an independently generated definition, it is imported and
validated on its own merits; for example, its A102413 recurrence does match the
current OEIS prefix.  Ten newly imported definitions lack a safe row-layout
alignment and therefore remain experimental: A099040, A103451, A105278,
A144217, A145677, A158821, A185740, A185911, A225117, and A258993.

## Web example catalogue (2026-09-02)

The web UI presents one flat menu of 34 distinct OEIS-labelled polynomial
families. It absorbs the former eight example buttons, removes the four
families duplicated between those buttons and the 30-entry OEIS menu, and
removes the category groups.

All 30 imported OEIS entries have recursive definitions in the corresponding
`ProofsOeis/A*.lean` files. The four additional distinct built-in families
(derangement excedances, Fibonacci matchings, Touchard polynomials, and Simsun
descents) also carry explicit recurrences. This does not imply that the web
recurrence search will rediscover every formula under its default bounds.

Long polynomial previews now parse and rejoin term signs before inserting the
ellipsis, so positive terms no longer render as `+ +` and a negative final term
renders with `-` rather than `+ -`.

Verification:

```text
python3 scripts/build_oeis_catalog.py --check                 passed
cargo test -q -p polytool --lib                              308 passed
cargo test -q -p polytool --test cli_oeis                      7 passed
cargo test -q -p polytool-mcp                                 22 passed
cargo test -q -p polytool --doc                                5 passed
cargo test -q -p polytool --test cli_bigint                   16 passed
cargo test -q -p polytool --test interlacing_api               5 passed
cargo clippy -q -p polytool -p polytool-mcp --all-targets --
  -D warnings -A clippy::manual-is-multiple-of
  -A clippy::needless-range-loop -A clippy::bool-assert-comparison
                                                               passed
```

The full library replay takes about 75 seconds because it regenerates every
row of all 73 benchmark recurrences.  A combined legacy recurrence-fixture
integration run reached its 180-second cap without reporting a failure; the
new independent sparse-definition replay completed successfully.

## Coupled Weyl recurrences

The main Rust worker owns these files for the current task:

- `src/recurrence.rs`
- `src/linalg.rs`
- `README.md`
- `HANDOFF.md`

The first library slice is implemented in `recurrence.rs`:

- `WeylOperator` stores `sum_d c_d(n,x) D_x^d` in normal order;
- `VectorRecurrence` represents
  `F_(n+1) = M(n,x,D_x) F_n + G(n,x)`;
- `find_vector_recurrence[_rational]` performs fixed-bound exact fitting;
- each output row is solved separately and reports its exact rank/nullity;
- non-identifiable rows are rejected by default (`require_unique = true`);
- final complete transitions are excluded from fitting and verified exactly;
- affine polynomial forcing is supported distinctly from the matrix state;
- `find_companion_vector_recurrence[_rational]` handles larger index lags,
  fitting only the final block row and inserting shift identities directly.

The convention is explicit: `states[k] = F_(first_index+k)`, and coefficients
in the transition `F_n -> F_(n+1)` are evaluated at the source index `n`.

Exact regression fixtures recover:

- the even/odd up-down-run Eulerian pair from Ma--Ma--Yeh--Yeh, Discrete Math.
  345 (2022), 112716;
- the lag-two type-B `1/k`-Eulerian system from Ma et al., EJC 27(3) (2020),
  P3.27, specialized to `k=1`;
- the affine q-integer recurrence `[n+1]_x = x[n]_x + 1`.

Verification on 2026-08-18:

```text
cargo test -q -p polytool --lib                 302 passed
cargo test -q -p polytool --lib recurrence::tests
                                                   51 passed
cargo test -q -p polytool --doc                   5 passed
cargo clippy -q -p polytool --lib -- -D warnings \
  -A clippy::manual-is-multiple-of \
  -A clippy::needless-range-loop                  passed
```

Potential follow-ups are an adaptive bound search, modular prefiltering for the
vector systems, JSON/CLI support, and optional nullspace-basis output for
exploratory non-identifiable fits. None is required for the fixed-bound exact
library API.

External review status: a read-only Claude Code review was attempted on
2026-08-18 at 08:34 UTC, but Claude exited before reading the diff because the
account session limit was reached (reported reset: 10:40 UTC). No review edits
were made.

## Web catalogue verification and deployment history

```text
inline JavaScript parse                                      passed
flat picker count, uniqueness, and key resolution           34/34
loadExample textarea/status checks                          34/34
recurrence provenance audit                                 34/34
positive/negative abbreviation regression cases              passed
desktop and 390px-wide headless-Chrome inspection            passed
standalone wasm-pack release build                            passed
git diff --check                                              passed
```

Commit `108916e` was pushed to monorepo `master`, projected to standalone
`main` as `fe697ec`, and deployed to `poly.symmetricfunctions.com` on
2026-09-02. Cache-busted public fetches matched the staged HTML, JavaScript,
and WASM byte for byte; a live headless-Chrome load reported no WASM or
JavaScript error.

The private proof-repository name was subsequently removed from all 26
imported-sequence comments before publication. The public examples now retain
only their mathematical descriptions and OEIS URLs; they do not mention the
author's related projects.

The recurrence-result page no longer generates or displays the bulky standalone
Python export. It keeps the compact recurrence JSON and links directly to the
documented `polytool recurrence-generate` command, which reconstructs exact
rows from that JSON. The library and CLI Python exporter remain available for
backward compatibility; only the browser payload and UI were simplified.

The browser's adaptive bounds are now recurrence depth 10, `t`-degree 5,
`n`-degree 5, and derivative order 5. In particular, this includes A059427's
cubic derivative coefficient `t - t^3`; a regression test confirms that degree
two fails and degree three finds the recurrence. The adaptive-mode tooltip
states that there is no elapsed-time or candidate-count cutoff: the search
stops only on a match, user cancellation, or exhaustion of its finite bounds,
and a failed exhaustive search may therefore take a long time.

Recurrence-option tooltips now mark their formulas with `data-tex` and render
them through the already loaded KaTeX runtime. Rendering is applied both to the
original tooltip nodes and to the floating tooltip layer; the alternating-sign
label and tooltip therefore display `(-1)^n` with an actual superscript.
The two focused web-crate tests, standalone WASM build, JavaScript parse,
code-card checks, and a headless-browser KaTeX/WASM load all pass.
Commit `1e8b865` was pushed to monorepo `master`, projected to standalone
`main` as `931a4e1`, and deployed on 2026-09-02. Cache-busted public files
matched the staged bundle byte for byte, and the live browser check passed.

A166073 was removed from the example picker because adaptive recurrence search
had to explore too much of the enlarged search space before finding its more
complicated recurrence. It was replaced by A008288, the Delannoy array read by
antidiagonals, whose row polynomials satisfy
`P(n) = (1 + t) P(n-1) + t P(n-2)`. A focused browser-crate regression test
checks that the full adaptive configuration finds this recurrence within ten
candidates. All remaining `Related project` comments were removed at the same
time. Commit `ccc9bf1` was pushed to monorepo `master`, projected to standalone
`main` as `d9f7361`, and deployed on 2026-09-02. Cache-busted public HTML,
JavaScript, and WASM matched the staged bundle byte for byte; a live
headless-Chrome load also passed.

## Uspensky/Descartes comparison

The main Rust worker owns these files:

- `src/root_count.rs`
- `src/lib.rs`
- `examples/bench_positive_real_rooted.rs`
- `README.md`
- `HANDOFF.md`

The exact Uspensky/Descartes comparison path is implemented and public. It
uses a strict Fujiwara bound, dyadic magnitude bands, reciprocal reduction
below `1`, homographic subdivision, and exact `BigInt` sign variations.  It
uses no finite fields or floating point.

The default one-signed real-rootedness path is now adaptive. It runs the cheap
Kurtz/Newton filters first, uses PRS generally, and selects Uspensky only when:

- the degree after stripping zero roots is at least 35;
- the two endpoint coefficients are equal; and
- an interior coefficient is at least `4^degree` times an endpoint.

The last two exact tests are invariant under scalar multiplication. This
conservative signature was chosen because degree or palindromicity alone
regressed important families. Square-free preprocessing is now shared with the
selected counter, avoiding a duplicate polynomial GCD in both the adaptive and
explicit PRS real-rootedness paths.

Release benchmark highlights from 2026-08-17 (single-process runs, so treat
sub-millisecond differences as noise):

```text
family                              primitive PRS    Uspensky
prod_{a=1}^{30} (x+a)                   0.88 ms       2.61 ms
prod_{a=1}^{80} (x+a)                  15.14 ms      82.43 ms
Eulerian (degree 35)                   45.82 ms      39.18 ms
Eulerian (degree 79)                    9.83 s        5.15 s
Narayana (degree 40)                    1.24 ms       9.41 ms
type-B Eulerian (degree 40)           182.37 ms     151.77 ms
Touchard (degree 40)                   55.54 ms      45.80 ms
Chebyshev T (degree 40)                 0.43 ms      61.05 ms
Chebyshev U (degree 40)                 0.43 ms      36.03 ms
Hermite (degree 40)                     0.33 ms      10.46 ms
```

These fair timings include the shared square-free preprocessing only once.
Eulerian, type-B Eulerian, and Touchard cross over around degree 35, while
Narayana, Chebyshev, Hermite, and evenly spaced linear-factor products strongly
favor PRS in the tested range. The benchmark example now covers all of these
families; Touchard is generated locally by its standard recurrence.

Verification:

```text
cargo test -q -p polytool --lib                 296 passed
cargo test -q -p polytool --test cli_bigint      16 passed
cargo test -q -p polytool --test interlacing_api  5 passed
cargo test -q -p polytool --doc                   5 passed
cargo clippy -q -p polytool --lib --examples -- \
  -D warnings -A clippy::manual-is-multiple-of \
  -A clippy::needless-range-loop                  passed
```

The two allowed Clippy lints are pre-existing in `cyclic_sieving.rs` and
`hstar_inequalities.rs`; neither file is part of this work.

## Repository ownership

- `/workspace/rust` branch `master` is the canonical monorepo history.
- All polytool changes are made under `polytool/` on `master`.
- Repository branch `main` is a generated standalone projection of this
  directory. Do not commit to it directly.
- After tested polytool changes have been committed and pushed on `master`, run
  `./scripts/sync-polytool-main.sh` from the monorepo root.

## Current state

The former independently rooted standalone history was backed up before
`main` was replaced by a `git subtree split --prefix=polytool` projection.
The projected tree includes the crate rename to `polytool` and all correctness
fixes from the August 2026 Rust review.
