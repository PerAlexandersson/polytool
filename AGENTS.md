# Workspace: Rust combinatorics projects

## Layout

```
rust/
  Cargo.toml                  ← workspace root (shared Cargo.lock)
  sym-poly/core/              ← Ring, Partition, Composition, BasisIndex, FormalSum
  sym-poly/sym/               ← Symmetric functions (6 classical bases)
  sym-poly/qsym/              ← Quasisymmetric functions (M, F, Ψ, Φ bases)
  sym-poly/multipoly/         ← Multivariate polynomials, divided differences, key polys
  polytool/           ← Polynomial toolkit (will be public)
  combinatoric-core/          ← Graphs, posets, chromatic; re-exports sym-poly
  combpoly/                   ← Permutation statistics + combpoly CLI
  kostka/                     ← Public crate: Kostka/LR coefficients, GT polytopes
  experiments/                ← temporary working directory for research binaries
```

## Operational guidance

- Run Rust commands with both `nice` and a 60 second timeout.
  Preferred pattern: `timeout 60s nice -n 10 cargo ...`
- Assume work in several projects may be happening simultaneously.
  Avoid broad cleanup, long-running background jobs, or edits outside the
  current Rust task unless they are clearly intended.

### sym-poly/ (multi-crate)
Modular workspace for symmetric polynomial algebras. Four crates, all generic
over `Ring` coefficients (i64, BigInt, Ratio<BigInt>, Ratio<i64>).

**sym-poly/core/** — Shared foundation:
- `ring` — Ring trait and implementations
- `partition`, `composition` — Indexing types with full combinatorial utilities
- `index` — `BasisIndex` trait (generic over Partition/Composition)
- `formal_sum` — `FormalSum<I, C>` generic linear combination
- `matrix` — Integer matrix utilities (multiply, transpose, invert)
- `transition_cache` — `TransitionCache<B>` per-algebra cached transition matrices

**sym-poly/sym/** — Symmetric functions:
- `SymmetricFunction<C>` with 6 bases (m, e, h, p, s, f), all 30 conversions
- Kostka coefficients, S_n characters, omega involution, multiplication
- Hall inner product, skew Schur (Jacobi-Trudi), plethysm, specializations

**sym-poly/qsym/** — Quasisymmetric functions:
- `QSymFunction<C>` with 4 bases: M (monomial), F (fundamental), Ψ (type 1), Φ (type 2)
- Power sum bases Ψ, Φ per Ballantine--Daugherty--Hicks--Mason--Niese (JCTA 2020)
- Normalized Ψ̃ = Ψ/z, Φ̃ = Φ/z utilities
- Omega involution (ω² = id, ω(Ψ_α) = (-1)^{n-ℓ(α)} Ψ_{α^r})
- Quasi-shuffle multiplication, Sym↔QSym maps
- (P,w)-partitions (Ψ̃-positive for naturally labeled posets)
- Chromatic quasisymmetric functions

**sym-poly/multipoly/** — Multivariate polynomials:
- `MultiPoly<C>` sparse multivariate polynomial type
- Divided difference operators ∂_i and π_i (Demazure)
- Key, atom, slide/glide, finite lock, Kohnert, Lascoux, and Schubert
  polynomial routines
- Nonsymmetric Hall-Littlewood `E_alpha(x; 0,t)` by operators, and full `q,t`
  permuted-basement nonsymmetric Macdonald polynomials by the SSAF filling
  formula

### polytool/
Univariate polynomial toolkit for combinatorial research. Uses primitive
integer PRS/root counting as the default exact real-rootedness backend, with
Bézout matrices as the main backend for interlacing and explicit matrix
certificates.

**Modules:**
- `polynomial` — `Polynomial<C>` with `CoeffRing`/`FieldRing` traits, arithmetic,
  derivative, evaluate, shift, reverse, dilate, GCD, division, Lagrange interpolation
- `real_rootedness` — public real-rootedness wrappers, explicit
  Bézout/Hermite comparison checks, strict/weak interlacing, log-concavity,
  ultra-log-concavity, palindromic, gamma-positivity, resultant,
  discriminant, Ehrhart h*-vector conversion, format_poly
- `root_count` — primitive integer PRS root counting and the default exact
  real-rootedness backend
- `sturm` — Sturm chains for exact root isolation
- `recurrence` — Adaptive recurrence search for polynomial sequences
- `sequences` — Eulerian, Narayana, type B Eulerian, Chebyshev T/U, Hermite

### combinatoric-core/
Facade crate: re-exports sym-poly-core and sym-poly-sym for backward compatibility.
Retains local modules for combinatorial structures not in the function algebra hierarchy.

**Re-exported** (from sym-poly): `Ring`, `Partition`, `Composition`, `Basis`,
`SymmetricFunction`, `kostka`, `transition` — all existing `use combinatoric_core::*`
paths continue to work.

**Local modules:**
- `graph` — Graph type with generators (complete, path, cycle, bipartite,
  multipartite, Ferrers board, unit interval, Petersen), operations (complement,
  induced subgraph, delete/contract edge), predicates (connected, bipartite,
  claw-free), matchings, independence number, graph6 parser
- `poset` — Poset type (Hasse diagram). Constructors: chain, antichain, fence,
  k-alternating, from Young/skew diagram, dual. Algorithms: linear extensions,
  order-preserving maps (backtracking + frontier DP), order polytope Ehrhart
  polynomial and h*-vector (BigRational, no overflow), P-Eulerian polynomial,
  natural relabeling. Frontier DP gives 100-18000x speedup over backtracking.
- `chromatic` — Chromatic symmetric function (Sym version, via deletion-contraction)
- `key_polynomial` — Key polynomials κ_{λ,σ} via Kogan faces of GT-polytopes,
  Ehrhart polynomials for key-Kostka coefficients

**Numeric policy / upgrade path:**
- Prefer one exact implementation plus one ergonomic wrapper, not two separate algorithms.
- For integer-valued combinatorial functions whose coefficients may grow, the target pattern is:
  `foo_bigint(...) -> Vec<BigInt>` as the canonical implementation, and
  `foo(...) -> Vec<i64>` as a checked convenience wrapper calling `to_i64().expect(...)`.
- For genuinely rational algorithms (for example Ehrhart interpolation), keep the exact
  `BigRational` version as the primary API.
- Avoid making every public API fully generic over coefficient type unless the generic
  version clearly improves reuse; generic helper layers are preferred over genericizing
  every algorithm.

### combpoly/
Permutation statistics and generating polynomials. Provides the `combpoly` CLI.

**Library modules:**
- `permutation` — Generation, pattern avoidance, backtrack image, filtered generation
- `statistics` — 18 statistics: des, exc, peak, valley, inv, maj, lrmax, fix, cyc, ...
- `polynomial_builder` — `build_generating_polynomial(objects, stat)`
- `order` — Bruhat/weak order ideals
- `word` — Words on Ferrers boards
- `parking` — Parking functions
- `catalan` — Catalan/Dyck path utilities
- `cayley` — Cayley permutations

### kostka/
Public crate. Kostka coefficients (DP), LR coefficients (GT polytope),
Ehrhart polynomials and h*-vectors of Gelfand-Tsetlin polytopes. Has CLI.

### experiments/
Temporary working directory for research exploration binaries.
Depends on combpoly and polytool. Flat layout with prefix-based grouping.

Notes:
- This folder is disposable and may be gitignored.
- Do not assume binaries here are stable, curated, or meant for long-term tracking.

## Naming conventions

- Use verbose, Mathematica-style names: `conjugate_partition`, `kostka_coefficient`,
  `to_monomial_basis`, `schur_symmetric`, `build_generating_polynomial`.
- CLI should be AI-friendly: use `clap` with detailed `--help`, support `--format json` output.

## Mathematica reference packages

Located in `~/AI-projects/combinatoric-tools/mathematica-packages/` (shared, read-only):
- `SymmetricFunctions.m` (~2000 lines) — full symmetric function package
- `CombinatoricTools.m` (~1600 lines) — partitions, compositions, tableaux, charge
- `NewTableaux.m` — SSYT/SYT generation, cylindric tableaux

Also: `~/Dropbox/AI-projects/docs/research/INTERLACING.md` has Mathematica code for the Bézout
interlacing algorithm, ready to add to the symmetric functions package.

## Non-nesting rook polynomial paper

Paper: `~/Dropbox/AI-projects/projects/Real-rooted-non-nesting-rooks/non-nesting-rooks.tex`
Overleaf companion: `Rook-Eulerian-Polynomials-and-permutation-ideals/`

### Status (2026-04-07)
- **Theorem**: R_μ(t) is real-rooted for every partition μ (non-nesting rook polynomial)
- **Proof approach**: Row-by-row induction, column-strip adjacent interlacing
- **Gap identified**: Original proof incorrectly applies Brändén Cor. 8.7 (requires pairwise
  interlacing, but input sequence is only adjacent-interlacing; verified computationally).
- **Revised approach** (Section 5 of tex): Uses shift lemma + Wagner cone instead.
  - The reduction `D_c/t ≪ G(c+1, s)  ⟹  G(c+1, s) ≪ G(c, s)` is clean (Lemma 5.6).
  - Shift lemma gives `D_c/t ≪ f_c` (proved, Lemma 5.7).
  - **Remaining gap**: Need `D_c/t ≪ f_{c+1}` (base case, verified 693/693)
    and `D_c/t ≪ G(c+1, s)` (full claim, verified 1026/1026).
    The inductive step s → s+1 via left cone requires `f_{s+1} ≪ D_c/t`, which fails.
  - ψ_c = (f_c − f_{c+1})/t has nonneg coefficients (693/693) and is real-rooted (693/693).
  - q-deformation R(q,t) = Σ q^{nest} t^k is real-rooted in t for all q ∈ [0,1] (verified n≤10).

### Experiment binaries
- `nn_rook_proof_check.rs` — First round: adjacent vs pairwise, sub-lemma, matrix N step
- `nn_rook_proof_check2.rs` — Partial sums G(c,s), diagonal interlacing, (A+tR)≪(B+tR)
- `nn_rook_proof_check3.rs` — 4-polynomial diagram, summand compatibility
- `nn_rook_qdeform.rs` — q-deformation (q^{nest} weight), column-strip for various q
- `nn_rook_qdeform2.rs` — Standard rook (q=1) properties, q^{crossing} weight
- `nn_rook_common.rs` — Common interlacing approach (fails for degree-gap reasons)
- `nn_rook_deform.rs` — Continuous deformation α∈[0,1], D/t properties (all pass)
- `nn_rook_proof_final.rs` — Final: all 6 properties P1–P6 checked

## What's next

### sym-poly/sym — Tier 2-4
- [ ] SSYT generation
- [ ] Ribbon Schur functions
- [ ] Cylindric Schur functions
- [ ] Hall-Littlewood P, T (needs `Polynomial<BigInt>` coefficients)
- [ ] Jack P, J (needs rational parameter)
- [ ] Schur Q, P (shifted symmetric functions)
- [ ] Macdonald P, H, J
- [ ] LLT polynomials, Delta/Nabla operators

### sym-poly/qsym
- [ ] Young quasisymmetric Schur (YQS) basis
- [ ] General (P,w)-partitions (non-natural labelings)
- [ ] q-chromatic quasisymmetric function (Shareshian-Wachs, polynomial coefficients in q)
- [ ] NCSym and Hopf algebra duality with QSym

### sym-poly/multipoly
- [ ] Operator recursion for full `q,t` nonsymmetric Macdonald polynomials
      (the filling formula is implemented)
- [ ] Transition-basis support for glide/lock/Kohnert/Lascoux families where
      the combinatorial transition matrices are tractable
- [ ] Larger checked examples for permuted-basement Macdonalds, glide
      polynomials, and Lascoux polynomials

## References
- Stanley, EC2 for symmetric function foundations
- Fisk, "Polynomials, roots, and interlacing" for Bézout matrix theory
- Postnikov (2005) for cylindric Schur functions
- Ballantine--Daugherty--Hicks--Mason--Niese, JCTA 2020, for QSym power sums Ψ and Φ
