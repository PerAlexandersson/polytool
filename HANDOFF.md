# Polytool handoff

## OEIS recurrence catalog

The host Codex supervisor completed the OEIS catalog increment at the user's
request.  No Rust worker is active and no catalog file ownership remains.
Files changed by the completed task are:

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

- `polytool::oeis` contains 145 recurrence-backed A-number functions, exact
  sparse recurrence decoding, dynamic lookup, row-range generation, and
  structural-zero restoration for OEIS output;
- 125 holdout-backed entries are available by default and 20 short-prefix
  entries require `--include-experimental`;
- 137 entries have a locally verified OEIS flattened-prefix mapping and support
  strict b-file output; the remaining 8 still support rows, triangles,
  polynomials, JSON, JSONL, and CSV;
- CLI commands are `polytool oeis list`, `polytool oeis info`, and
  `polytool oeis generate`; MCP tools are `list_oeis_sequences`,
  `get_oeis_sequence`, and `generate_oeis_rows`;
- `scripts/build_oeis_catalog.py` imports the 73 machine recurrence fixtures
  and the plain-file OEIS queue, validates queue recurrences against cached
  rows, emits exact sparse Rust definitions, and supports `--check` drift
  detection.

Seven unsafe queue entries are excluded.  `A123125` was already tagged
`invalid_recurrence`.  Independent replay additionally found that `A102413`,
`A153520`, `A153521`, `A201701`, `A271704`, and `A285066` do not reproduce their
cached rows; the first three currently have stale holdout-verification tags in
the source queue.

Verification:

```text
python3 scripts/build_oeis_catalog.py --check                 passed
cargo test -q -p polytool --lib                              307 passed
cargo test -q -p polytool --test cli_oeis                      6 passed
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
