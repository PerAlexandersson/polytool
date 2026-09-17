# Handoff

Full chronology through the start of the current repository-health batch is in
[`docs/HANDOFF_ARCHIVE_THROUGH_2026-09-16.md`](docs/HANDOFF_ARCHIVE_THROUGH_2026-09-16.md).

## Active: repository health remediation — 2026-09-16

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
Both are project-specific refuted routes.  Its
`cone-transform` mode reconstructs exact rational candidate transfer matrices:
all fourteen matrices through order five have full input rank, and every
negative entry has an outer source shape of first row at least three.  The
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
