# sym-poly workspace guide

## Dependency direction and crate boundaries

The four crates are deliberately layered:

```text
combinatoric-core
        |
  sym-poly-core -> sym-poly-multipoly -> sym-poly-sym -> sym-poly-qsym
```

`sym-poly-sym` also directly uses `sym-poly-core` and
`sym-poly-multipoly`; `sym-poly-qsym` uses core, sym, and combinatoric-core.
Keep shared indexing, coefficient abstractions, tableaux, fillings, and exact
linear algebra in `core`. Keep sparse multivariate algebra and nonsymmetric
families in `multipoly`, symmetric functions and symmetric-function families
in `sym`, and quasisymmetric families in `qsym`. Do not move an algorithm
upstream just to shorten an import, and do not add a reverse dependency.

The root `README.md` is a roadmap/breadcrumb for website work. The crate
READMEs are the concise capability maps and should name actual modules,
examples, and conventions.

## Shared conventions and exactness

- `Partition` is normalized weakly decreasing with zeros removed;
  `Composition` strips trailing zeros; `WeakComposition` preserves length.
  Partitions index `Sym`, compositions index `QSym`, and weak compositions
  encode multivariate exponent/shape vectors.
- `Ring` implementations include `i64`, `BigInt`, `Ratio<i64>`, and
  `Ratio<BigInt>` (with polynomial/rational-function layers built on top).
  Integer-ring exact division panics when a quotient is not integral; rational
  rings represent the quotient.
- Generic coefficients do not imply every algorithm is arbitrary precision.
  Transition caches and several conversion/Kostka paths store `Vec<Vec<i64>>`
  or use i64 helpers. A `BigInt` coefficient choice therefore does not upgrade
  current i64 transition/Kostka data. Multivariate exponents and homogeneous
  degrees are `u32` and checked for overflow.
- Respect basis order, reading-word, basement, variable-index, and permutation
  conventions in the module rustdoc. Use explicit constructors rather than
  silently changing normalization or padding.

## Development rules

Search with `rg` across all workspace crates before adding a generator,
conversion, transition matrix, filling rule, or helper. Reuse a canonical
implementation and make wrappers call it; do not create duplicate algorithms
in a sibling crate. One-use experiments are intentionally not in this repo and
are not a substitute for a public API test or a checked example.

Every new public type/function needs rustdoc, focused unit tests for invariants
and at least one small source/paper example when applicable. Keep examples in
the owning crate's `examples/` directory and cite them from website/TeX work
with a short `Related Rust` breadcrumb. Run focused tests first, then the
owning package's test/check command; keep Cargo build output in the external
cache rather than the Dropbox checkout. Do not edit generated `.lake`, `target`,
or experiment output.

When a change crosses crate boundaries, update the owning README and relevant
example in the same checkpoint. Verify package names and dependency edges with
`cargo metadata --no-deps --format-version 1`; do not infer them from old notes.
