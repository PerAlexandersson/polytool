# Rust combinatorics workspace

This workspace contains exact combinatorics libraries, polynomial tools,
research applications, and checked examples for `symmetricfunctions.com`.

Start with [docs/LIBRARY_GUIDE.md](docs/LIBRARY_GUIDE.md): it maps mathematical
tasks to the owning crate, public symbols, examples, coefficient conventions,
and the required find-before-write workflow. Operational rules for contributors
and agents live in [AGENTS.md](AGENTS.md).

## Crate map

- `combinatoric-core`: foundational partitions, compositions, graphs, posets,
  permutations, exact linear algebra, and the shared `Ring` trait.
- `sym-poly/core`: shared algebra/tableau infrastructure that depends on and
  re-exports the canonical foundation types.
- `sym-poly/multipoly`, `sym-poly/sym`, `sym-poly/qsym`: respectively
  multivariate/nonsymmetric, symmetric, and quasisymmetric polynomial layers.
- `polytool`: exact univariate polynomial algorithms and CLI; `web` and `mcp`
  are adapters around that library.
- `combpoly`: permutation, word, parking-function, and related generating
  polynomial exploration on top of `polytool`.
- `polynomial-lab`: structured evidence ledger for real-rootedness and
  interlacing projects, also backed by `polytool`.
- `kostka`: legacy reproducibility implementation; see its README for the
  maintained successor and new-work policy. `KTT-search` consumes it.
- `flagged-lorentzian`: targeted research crate using `sym-poly-core`.
- `experiments`: deliberately disposable research binaries consuming the
  reusable crates. Most are intentionally ignored; promote only stable library
  cores and selected canonical examples.

The dependency direction begins with `combinatoric-core`, then
`sym-poly-core`, then `sym-poly-multipoly`/`sym-poly-sym`/`sym-poly-qsym`.
`combinatoric-core` is not a facade over `sym-poly`.

## Polytool standalone branch

`master` is the canonical history for this workspace. The `main` branch is the
generated standalone view of `polytool/`; do not develop directly on it. After
committing and pushing eligible `polytool` work on `master`, publish it with:

```sh
./scripts/sync-polytool-main.sh
```

The script requires a clean `polytool/` tree and refuses to overwrite a
diverged `main` branch.
