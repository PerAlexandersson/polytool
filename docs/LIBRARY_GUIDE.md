# Library capability guide

Use this guide to select an existing API before writing code. It is a routing
map, not a replacement for each crate's README, module documentation, tests,
or examples.

## Choose a crate by task

| Task | Start with | Public entry points / checked examples |
| --- | --- | --- |
| Partitions, compositions, graphs, posets, permutations, integral linear algebra | `combinatoric-core` | `Partition`, `Composition`, `Graph`, `Poset`, `smith_normal_form`; `combinatoric-core/examples/optimist_sorting.rs` |
| Shared tableaux, SSAFs, exact finite-field/rational-function tools, formal sums | `sym-poly-core` | `Tableau`, `SkewTableau`, `Ssaf`, `FormalSum`, `TransitionCache`; `sym-poly/core/examples/hamel_goulden_site_example.rs` |
| Symmetric functions, chromatic and bond-lattice symmetric functions, LLT/Macdonald-related constructions | `sym-poly-sym` | `SymmetricFunction`, `Basis`, `chromatic_symmetric`, `graph_mobius_symmetric_function`, `chromatic_mobius_symmetric_function`, `unicellular_llt`, `delta_modified_macdonald`; `sym-poly/sym/examples/*_site_example.rs` |
| Quasisymmetric functions, P-partitions, QSym Schur variants, QSym chromatics | `sym-poly-qsym` | `QSymFunction`, `QSymBasis`, `p_partition_generating_function`, `qsym_schur`, `young_qsym_schur`, `chromatic_qsym`; `sym-poly/qsym/examples/qsym_schur_degree4.rs`, `p_partition_site_example.rs` |
| Sparse multivariate and nonsymmetric polynomial bases/operators | `sym-poly-multipoly` | `MultiPoly`, `key_polynomial`, `atom_polynomial`, `schubert_polynomial`, `flagged_schur`, `nonsymmetric_macdonald_filling_formula`; `sym-poly/multipoly/examples/key_site_example.rs`, `schubert_site_example.rs` |
| Exact univariate root/interlacing, recurrence, h* and standard sequences | `polytool` | `Polynomial`, `is_real_rooted`, `check_weak_interlacing`, `find_recurrence_adaptive`, `eulerian_polynomials_bigint`; `polytool/examples/bench_bezout_vs_sturm.rs` |
| Permutation/word/parking-function generating polynomials and CLI scans | `combpoly` | `build_generating_polynomial` and the `combpoly` CLI; see `combpoly/README.md` |
| Structured evidence, counterexamples, and proof-search records | `polynomial-lab` | `PolynomialFamilyRegistry`, `real_rooted_evidence_id`, and the `poly-lab` CLI; see `polynomial-lab/tests/fixtures/minimal_lab/` |
| Historical GT/Ehrhart, Kostka, and LR reproducibility | `kostka` | its CLI and modules `kostka_dp`, `ehrhart`, `lr`; read `kostka/README.md` first because new feature work belongs in Ehrcalc |

`experiments/` is for bounded, use-once research probes. Its large ignored
binary collection is intentional. Do not attempt to track, rename, test, or
catalogue it wholesale. When a result becomes reusable, promote its stable core
to the owning library with tests; add a tracked example only when it is a useful
canonical demonstration. It is a standalone Cargo workspace: use
`--manifest-path experiments/Cargo.toml` for intentional runs. The retired
database-backed `KTT-search` application is likewise isolated behind its own
manifest and lockfile; new GT/Ehrhart work belongs in Ehrcalc.

## Architecture and names that must not be conflated

The library stack is `combinatoric-core` → `sym-poly-core` →
`sym-poly-multipoly`/`sym-poly-sym` → `sym-poly-qsym`. `sym-poly-core`
re-exports canonical foundation types; it does not supply them to
`combinatoric-core`.

- `Partition` and `Composition`: use `combinatoric_core::{Partition,
  Composition}` as the canonical types, or the compatible re-exports from
  `sym_poly_core`. `kostka::Partition` is a legacy, separate type: convert at
  its boundary rather than mixing the APIs.
- `Ring`: `combinatoric_core::Ring` (also re-exported by `sym_poly_core`) is
  for the symmetric-polynomial stack. `polytool::{CoeffRing, FieldRing}` are
  its separate univariate coefficient traits.
- `Polynomial`: `polytool::Polynomial<C>` is a dense univariate type;
  `sym_poly_core::UnivariatePolynomial` is a distinct shared-algebra type;
  `sym_poly_multipoly::MultiPoly<C>` is sparse multivariate. Do not assume
  conversions or coefficient order without reading the target API.
- `key_polynomial`: `sym_poly_multipoly::key_polynomial` constructs a
  polynomial from a composition. `combinatoric_core::key_polynomial` contains
  legacy combinatorial weight routines with different inputs and outputs.
- Chromatic APIs have different codomains: `sym_poly_sym::chromatic_symmetric`
  returns a symmetric function, while `sym_poly_qsym::chromatic_qsym` returns
  a quasisymmetric function. Pick the theorem's target algebra first.

## Exact arithmetic and result interpretation

- Coefficient vectors accepted by `polytool` convenience functions are in
  ascending degree order. Use its `*_bigint_coeffs` APIs when `i64` is not a
  proved bound. `Option<bool>` interlacing results can be `None` for an invalid
  directed degree relation, rather than a failed mathematical assertion.
- The symmetric-function crates are generic over their `Ring`; QSym power-sum
  bases require rational coefficients. Select `BigInt` or `BigRational` when
  growth or division demands it, and do not add lossy casts for convenience.
- Some legacy and convenience APIs deliberately expose `i64` results and may
  reject values that do not fit. Keep an exact primary implementation whenever
  a combinatorial quantity can grow, and add checked ergonomic wrappers only
  where the owning crate's convention supports them.
- Root checks, finite enumerations, and recurrence fits are evidence or exact
  computations for their stated inputs; they are not proofs beyond those
  inputs. Record method and checked range in `polynomial-lab` when appropriate.

## Find before write

1. Run `cargo metadata --no-deps --format-version 1` and locate the owning
   crate; inspect its local guide and README.
2. Search `src/`, `tests/`, and `examples/` with `rg` for the object, theorem,
   and likely operation. Read the candidate's signature, coefficient type,
   indexing convention, and tests before designing a replacement.
3. Reuse or extend the lowest suitable owner. Do not add a dependency upward
   through the symmetric-polynomial stack, duplicate a general routine in an
   experiment, or use a legacy `kostka` type in the canonical stack.
4. For a reusable feature, add focused tests and documentation; add a small
   tracked example when it clarifies a stable public use. Leave truly local
   probes in `experiments/` without turning them into maintenance obligations.
5. Verify with the narrowest relevant command, preserve concurrent changes,
   and update project handoff/documentation according to the active task.
