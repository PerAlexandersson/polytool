# Handoff

## Completed: forgotten-basis investigation — 2026-09-25

All-degree partition graphs have strictly positive forgotten coefficients
after multiplying by `(-1)^(|lambda|-length(lambda))`, using `f=omega(m)`.
`docs/partition-graph-forgotten.md` records the star-expansion and set-partition
inclusion-exclusion proof, normalization, examples, and exact evidence.
All 41 two-partitions through degree 10 were checked at three variable counts
(123 cases), with complete support, the predicted all-ones coefficient,
monomial roundtrips, and an independent Schur-conjugation conversion route.
Only this handoff and the mathematical note were changed; existing library
generators and conversions were reused by a bounded temporary Rust driver.
These note files are released. The prior no-push boundary below still applies.

## Completed: all-degree partition-graph basis — 2026-09-25

`sym-poly-sym::partition_graph_symmetric(&lambda, n)` implements the
all-degree basis from Alexandersson--Shapiro, arXiv:1104.0589 Section 2.2.
The user explicitly selected this family; the even-degree square-graph basis
is not included. The generator uses disjoint outward stars, BigInt coefficients,
explicit variable count, and Definition 1's full permutation sum. The reusable
`symmetrized_graph_monomial` expands via MultiPoly and exponent stabilizers,
without enumerating S_n; Schur conversion reuses the existing Sym API.
The checked `partition_graph_schur` example lists degrees 2--6 at n=2*d,
or accepts `degree [num_vars]`. README and library-guide routing are updated.

Eight new tests pass: literal permutation-sum agreement (including isolates,
repeated edges and disconnected stars), known even/odd Schur expansions,
all 2-partition diagonal coefficients and Schur roundtrips through degree 8,
translation invariance, n=30 BigInt coefficients, boundary/errors, orientation,
and relabeling. Full package tests pass (223 passed, 3 pre-existing ignored);
the subsequently strengthened nonzero translation probes also pass. The
example, strict Rustdoc, formatting, and diff checks pass. Focused Clippy has
no findings in the new module/tests/example; existing unrelated warnings remain.
Logs are `/tmp/partition-graph-{sym-tests,doc,clippy}.log` and
`/tmp/partition-graph-clippy.jsonl`.

All task files are released after the local checkpoint on
`feat/translation-invariant-graph-bases`, based on `62d6e6f`. No push was made:
this branch includes 13 pre-existing unpublished master commits, including a
prior handoff's explicitly unpublished work, which this task does not publish.
The i64 limitation in existing Schur transition data is documented and a
separate exactness follow-up is proposed in
`suggestions/sym-transition-exactness-2026-09-25.md`.
Paper-cache lookup found no copy; its arXiv fetch failed with HTTP 406 from
the metadata endpoint. The PDF was read through the browser instead.
No Lean work or additional worker was involved.

## Completed: colored exchange library — 2026-09-22

The user-authorized 3-Matroids integration adds `combinatoric-core::colored_exchange`:
validated fixed-content families, deterministic symmetric cycle witnesses,
explicit search-budget exhaustion, exact BigUint marginals, and verified
product/binary-port/private-block certificates. The API has no 128-word cap.
Its example demonstrates witnesses, obstructions, weighted coefficients,
and decompositions. README and library-guide routing are updated.

Claude implemented the module, export, 11 integration tests, and example.
Root reviewed and integrated it, fixed actual traversal-budget accounting,
removed finite sentinels from unbounded searches, compressed arbitrary raw
color labels, rejected unused private-certificate data, and guarded corrupt
binary count overflow. A twelfth integration test covers these defects.
All eight unit tests and twelve integration tests pass. Independent brute
force agrees on all 4,095 quota-(2,1,1) supports (405 successes); tests also
cover larger-than-128-word families and arbitrary-precision coefficients.
The example, strict Rustdoc, and focused formatting pass. Clippy reports
pre-existing findings outside the owned files (44 in the lib-test target,
38 duplicated from the library); none concern this addition.

All task files are released after the local checkpoint. No nested Claude
process remains active. No push/publication was authorized or performed.
The project consumer and further exact classification replay are recorded in
`/workspace/projects/3-Matroids/RUST_INTEGRATION.md` and its handoff.
No analytic log-concavity conclusion is built into the certificate API.

## Completed: A355089 statistic refinements — 2026-09-18

`combpoly::parity_pattern` now uses one generic maximum-and-gap traversal for
exact counts and five exact polynomial refinements: descents, ascents, major
index, ascent-major index, and left-to-right minima. The public
`parity_preserving_123_statistic_polynomials` API returns ascending-degree
`BigUint` rows, and the checked `a355089_statistics` example emits rows in
Polytool's input format. Independent permutation generation verifies every
refinement through `n=8`; specialization at `t=1` verifies all five through
`n=20`. Reverse-complement and the disjoint record-set decomposition prove
that the left-to-right-minimum polynomial is palindromic at every even length;
this has a focused regression through `n=20`. All 200 Combpoly library tests,
both examples, strict Rustdoc,
formatting, and focused strict Clippy pass.

Project-local Polytool searches found no scalar or even/odd-vector recurrence
in the recorded bounded spaces. The ignored
`experiments/src/bin/a355089_stat_recurrence.rs` remains a use-once search
driver, not maintained library functionality. No generic permutation or
Polytool source was changed. All files from this batch are released.

## Completed: parity-alternating 123-avoidance DP — 2026-09-18

`combpoly::parity_pattern` now provides an exact `BigUint`, count-only memoized
DP for both parity orientations of `321`-avoidance and the canonical A355089
generator obtained by value complement. The checked `a355089` example
reproduces generation. Its independent brute-force test covers `n <= 9`, all
published A355089 terms through `n=30` are regression-tested, and the preserving
orientation agrees with A354208 through `n=20`. All 197 Combpoly library tests,
the example target, strict Rustdoc, formatting, and focused Clippy pass. Strict
whole-crate Clippy initially reported only 14 pre-existing findings outside the
owned files; allowing precisely those existing lint classes exposed no new
finding. The ignored `experiments/src/bin/parity_123.rs` and generic
`combpoly/src/permutation.rs` were profiled/read but not modified. All files
from this task are released.

Full chronology through the start of the current repository-health batch is in
[`docs/HANDOFF_ARCHIVE_THROUGH_2026-09-16.md`](docs/HANDOFF_ARCHIVE_THROUGH_2026-09-16.md).

## Active: repository health remediation — 2026-09-16

`/root` (2026-09-17) has released this handoff and the ignored project
consumer `experiments/src/bin/chordal_clawfree_boundary_scan.rs` after a
bounded higher-rank transfer run.  No maintained library source was owned or
changed.

The chordal-claw-free Schur project has released all Rust files.  Reusable
work now includes `Graph::join`, the monomial-basis universal-vertex operator
`chromatic_symmetric_add_universal_vertex`, and cached degree-dependent
permutations, cycle indices, and character tables in the exact clique-boundary
routine.  The universal-vertex API documents and tests its exact positive
action on every two-column Schur function.  The family-specific scan remains in the ignored
`experiments/src/bin/chordal_clawfree_boundary_scan.rs` consumer.  Its new
`positive-compression` mode refutes the naive model obtained by leaving the
new universal label unused in the old harmonic module; that model cannot
produce the lower outer two-row shapes present in the true transfer.  The
`algebra-cone` mode implements exact insertion of the universal choice in
each PEO transposition code and shows that applying this map to the naive
orthogonal primitive-idempotent section also gives the wrong refined kernel.
The further `orbit-compression` mode tests a diagonal projector onto
`c`-subsets containing exactly `c-b` boundary labels.  It matches all `29`
endpoint columns through order five and all `18` endpoint columns in the
order-six spanning samples, but fails for mixed boundary content already at
`(n,k,c,b)=(5,2,2,1)`: at `C=1` it gives `(1,9)/10`, rather than the exact
`(7,9)/16`.  This is project-specific finite evidence and a refuted general
model, not shared-library functionality.  During this checkpoint `/root`
temporarily owned only this handoff and the ignored experiment driver; both
are released.  The failed compression models are project-specific refuted
routes.  Its
`cone-transform` mode reconstructs exact rational candidate transfer matrices:
all fourteen matrices through order five have full input rank, and every
negative entry has an outer source shape of first row at least three.  The
same statement now holds for the full-rank order-seven spanning
reconstructions at boundary sizes two, three, and four (`28`, `19`, and `10`
input coordinates, respectively).  The ignored driver permits these bounded
order-seven runs; this project-specific cap change is not maintained-library
functionality.  The
order-six unlabeled universal-clique audit covered 17 relevant graphs and 361
positive boundary states without a failure; the attempted order-seven audit
was interrupted before completion and has no claimed result.

Verification passed: 269 `combinatoric-core` tests, 135 `sym-poly-sym` unit
tests plus all integration and documentation tests, strict Rustdoc, the
checked universal-vertex example, and the project-specific exact transform
sweeps.  Focused strict Clippy reports only the 38 pre-existing
`combinatoric-core` findings and 36 pre-existing `sym-poly-sym` findings; the
changed source files add none.

The chordal-claw-free Schur project has again released its ignored experiment
driver.  `experiments/src/bin/chordal_clawfree_boundary_scan.rs` now has a
`detail` mode for individual channel tables and a `universal` mode for the
universal-vertex preservation probe.  The latter finds the unrestricted claim
false (the edgeless three-vertex graph cones to a claw), but finds no failure
among 3,181 boundary-positive chordal claw-free inputs through order five when
the cone also remains claw-free.  These are project-specific proof-search
checks and remain outside the maintained workspace.

The chordal-claw-free Schur project has released all Rust file ownership after
promoting its reusable core.  `combinatoric-core` now provides the complete
`3`-sun clique-blowup generator and clique-rooted perfect elimination orders;
`sym-poly-sym` now provides exact cleared chordal clique-boundary Schur
channels through order nine.  `chromatic_symmetric` uses one stable-set-
partition traversal shared by all monomial types.  The checked example and
library guide expose the new boundary API.  The family-specific scan remains
an ignored driver at
`experiments/src/bin/chordal_clawfree_boundary_scan.rs`, outside the maintained
workspace.

Verification passed: all 268 `combinatoric-core` tests, all 134
`sym-poly-sym` unit tests plus its integration and documentation tests, strict
Rustdoc for `sym-poly-sym`, the checked boundary example, and the explicit
negative-intermediate/positive-completion cross-check.  Focused Clippy adds no
warning in the new module; its 36 reports are pre-existing elsewhere in
`sym-poly-sym`.  The project scan found no negative outer-private-clique
channel through total order eight (28 blowups, 243 nonzero coefficients),
whereas the overlapping central boundary is already negative for the bare
`3`-sun.  This sharply identifies the outer private clique as the viable
interface for the next proof.
The shared-library checkpoint is committed locally.  Its routine push to the
existing SSH remote failed because this environment has no usable public-key
authentication; no remote or credential configuration was changed.

The host supervisor owns final integration, CI observation, documentation,
suggestions, and this handoff. All bounded source workers have completed and
released their files; there is no active source ownership in this repository.
Commits `4a897bf`, `e777c57`, and `f7d136f` are pushed on `master`.

Kostka retirement is complete in standalone commit `95cc655` and workspace
commits `b4f1174` and `06340f3`; GitHub Actions run `35095705694` passed. The
historical `kostka` submodule and database-backed `KTT-search` application remain in Git
for reproducibility, but are excluded from the maintained root workspace and
have no dedicated CI job. New GT/Ehrhart work routes to Ehrcalc;
LR/Kostka/Schur compatibility work routes to lrcalc-rs. The Kostka repository
README and GitHub description identify both successors. The tracked-experiments
CI job still checks out the frozen submodule because four retained GT
experiment binaries consume it.

The final performance follow-up streams Combpoly permutation statistics,
flagged-tableau weights, QSym quasi-shuffle multiplicities, and both general
and chordal graph sink-polynomial independent sets. Noncrossing and nonnesting
matchings now prune invalid branches directly. Full affected-package tests pass:
266 combinatoric-core, 194 Combpoly library plus 2 binary, 204 multipoly, and
136 QSym tests, together with all affected examples.

The 961 ignored experiment binaries remain local and untracked. Broad builds
must not compete with the active eight-process KTT computation; focused checks
use external targets, reduced priority, and separate CPU affinity when safe.

The root workspace now excludes the standalone `experiments` and retired
`KTT-search` workspaces, each with a committed lockfile. The Combpoly inner Git
metadata was preserved as a verified bundle and recoverable metadata directory
under `/mnt/2TB-Babel/ai-storage/source/repo-metadata-backups/` before removal
from the working tree. The nested Kostka metadata-only checkpoint is
`a4f01fe`, pushed to `origin/main`.

The remaining deliberate debt is narrow: add durable performance benchmarks,
resolve the 38 pre-existing strict-Clippy findings in combinatoric-core, and
upgrade the isolated experiments/KTT MySQL chain when its upstream
future-incompatibility warning is resolved. See
[`suggestions/rust-repository-health-2026-09-16.md`](suggestions/rust-repository-health-2026-09-16.md).

The complete workspace documentation builds with `RUSTDOCFLAGS="-D warnings"`;
the root CI enforces this alongside the declared Polytool MSRVs. Kostka's two
documentation repairs are pushed separately as `ec1383c`.

The first MSRV run found post-1.82 `is_multiple_of` calls in Polytool; all
library, CLI, and example occurrences are now 1.82-compatible. Strict Clippy
also passes for Polytool and its MCP companion. The full 356-test Polytool
library suite passes after those changes, as does an all-target check on the
exact Rust 1.82.0 toolchain for Polytool and its web companion.

Polytool's `oeis info <A-number> --json` now embeds `recurrence_data`, the
standard `polytool.recurrence.v1` object with exact coefficients and initial
polynomials, while preserving the existing human-readable recurrence and
language exports. Focused CLI tests, recurrence JSON fixture replay, formatting,
and strict binary Clippy pass.

The OEIS catalog's former 30 experimental entries are now reconciled. Nineteen
queue recurrences with complete OEIS-row alignment are correctly `validated`;
the corrected A099040 import is also validated. Ten entries that cannot align
to complete OEIS rows are no longer bundled: A103328 has reversed structural
zeros and already has its valid recurrence in OEIS, seven are stale row-boundary
artifacts, and A266178/A266298 lost cellular-automaton zero positions. The
generator now rejects unaligned queue and Lean imports. Regeneration also drops
the corrected-away A062154 entry. The catalog contains 774 entries (125
verified, 649 validated, zero experimental). Generator checks, 9 focused OEIS
library tests, 9 CLI tests, all 31 MCP tests, formatting, and strict all-target
Clippy pass.

## Durable policy

- `master` is the canonical monorepo branch; `polytool/main` is generated only
  by `scripts/sync-polytool-main.sh` after its prerequisites pass.
- Use root Cargo commands for maintained crates and explicit `--manifest-path`
  commands for experiments or legacy KTT.
- Keep ignored experiments disposable. Promote only stable reusable kernels
  into documented and tested library APIs.
- Keep build output in `/cargo-target/ai-projects`, not this source tree.
- Preserve unrelated worker changes and record new ownership here before
  adopting a dirty file.
