# Project Guide: combpoly

## Purpose

Library and CLI for generating polynomials from combinatorial objects
(permutations, words, parking functions) by counting with a statistic.

Core operation: given a set of objects and a statistic,
`build_generating_polynomial` returns the coefficient vector where
`coeffs[k] = #{w : stat(w) = k}`.

## CLI tool: `combpoly`

```bash
cargo run --release -- poly --perms 7 --avoiding 312 --stat des --real-rooted
cargo run --release -- scan --size 7 --avoiding 312 --ideal bruhat --stat exc
cargo run --release -- list --perms 5 --avoiding 312
cargo run --release -- recurrence --perms --max-n 12 --avoiding 312 --stat des --auto
```

## Library modules (`src/`)

| Module | Contents |
|--------|----------|
| `permutation.rs` | Generation, pattern avoidance, backtrack_image, filtered_permutations |
| `statistics.rs` | 18 statistics: des, exc, peak, valley, inv, maj, lrmax, fix, cyc, ... |
| `polynomial_builder.rs` | `build_generating_polynomial(objects, stat)` |
| `order.rs` | Bruhat/weak order ideals |
| `word.rs` | Words on Ferrers boards |
| `parking.rs` | Parking functions, run-sorted variants |
| `catalan.rs` | Catalan/Dyck path utilities |
| `cayley.rs` | Cayley permutations |
| `fixed_descent.rs` | Fixed-descent insertion and transfer data |
| `lattice_path_matroid.rs` | Lattice-path-matroid h* computations |
| `rook_placements.rs` | Ordinary and non-nesting rook generators and packets |

## Dependencies

- `polytool` owns real-rootedness, interlacing, coefficient-vector analysis,
  formatting, and recurrence search. Reuse it instead of adding a second
  polynomial implementation here.
- `combinatoric-core` owns general foundational graph, poset, partition, and
  permutation structures. Do not create a reverse dependency from that crate
  to Combpoly.
- `sym-poly-*` owns symmetric, quasisymmetric, and multivariate function
  algebras. Combpoly should produce combinatorial data rather than duplicate
  those algebra types.

## Exploration binaries

Research binaries are in the `experiments/` crate (workspace sibling), not
here. Most are deliberately ignored and use-once. Promote a generally useful
generator into this library only with a documented contract, exact tests, and
an independent small-instance check; do not move the surrounding scan loop.

## Implementation conventions

- Search this crate and `../docs/LIBRARY_GUIDE.md` before adding enumeration,
  pattern, statistic, rook, polynomial, or interpolation helpers.
- State whether permutations are zero- or one-indexed and use the shared
  pattern APIs. Do not encode a named pattern through an undocumented custom
  inequality.
- Coefficient vectors use ascending degree. Trim only when trailing zeros are
  mathematically inessential; h* dimension padding is an important exception.
- For coefficients that can grow, implement a `BigInt` API as the canonical
  routine and make any `i64` API a checked convenience wrapper.
- Prefer iterators or callbacks for large object families. Do not require full
  materialization merely to compute a distribution.
- Public APIs need rustdoc covering the mathematical definition, empty input,
  invalid input, indexing, arithmetic limits, and complexity.
- Add known examples plus an independent implementation or exhaustive small
  comparison. A research observation alone is not a regression oracle.
- Update `README.md`, `../docs/LIBRARY_GUIDE.md`, and a runnable example when
  adding a named public family.

## Paper: Backtrack permutations

Associated research paper in `paper/Backtrack-permutations.tex`.
See the paper directory for current status.

## OEIS connections

| Triangle | OEIS | Context |
|----------|------|---------|
| Narayana (des on Av_132 etc.) | A001263 | h-vector of associahedron |
| Av_321 + des | A091156 | Dyck paths by long ascents |
| Av_132 + des | A048994 | Unsigned Stirling 1st kind |
| Peak on Av_132 (Class A) | A091894 | Dyck paths by ddu's |
| Peak on Av_312 (Class B) | A236406 | 123-avoiding by peaks |
| Exc on Av_312 | NOT IN OEIS | New, fails real-rootedness |
| Exc on Av_231 | NOT IN OEIS | New, real-rooted, no recurrence |

## Style notes

- Paper uses amsart, explicit `\ref` (no cleveref).
- LaTeX macros: `\backtrack`, `\symS`, `\Av`, `\des`, `\exc`, `\peak`, `\ltr`, `\oeis{Annnnnn}`.
