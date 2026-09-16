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
- `kostka`: retired historical submodule, excluded from the maintained Cargo
  workspace and CI. Use [Ehrcalc](https://github.com/PerAlexandersson/ehrcalc)
  for GT/Ehrhart and related exact workflows, or
  [lrcalc-rs](https://github.com/PerAlexandersson/lrcalc-rs) for the
  LR/Kostka/Schur compatibility surface. The old database-backed `KTT-search`
  application is retained only for reproducibility.
- `flagged-lorentzian`: targeted research crate using `sym-poly-core`.
- `experiments`: deliberately disposable research binaries consuming the
  reusable crates. Most are intentionally ignored; promote only stable library
  cores and selected canonical examples. It is a separate Cargo workspace so
  local ignored probes cannot alter the maintained root target graph.

The dependency direction begins with `combinatoric-core`, then
`sym-poly-core`, then `sym-poly-multipoly`/`sym-poly-sym`/`sym-poly-qsym`.
`combinatoric-core` is not a facade over `sym-poly`.

## Cargo and CI boundaries

The root workspace contains only maintained libraries and adapters. Ordinary
root commands do not discover the large local experiment forest, the retired
Kostka submodule, or the retired KTT application's MySQL dependency graph. Use
an explicit standalone manifest only when reproducing a historical result:

```sh
cargo test --locked --manifest-path experiments/Cargo.toml --bin <name>
cargo test --locked --manifest-path kostka/Cargo.toml --lib
cargo check --locked --manifest-path KTT-search/Cargo.toml
```

The root GitHub workflow tests each maintained package separately and verifies
the tracked experiment artifacts in a clean checkout. There is no dedicated
Kostka/KTT job; the experiment job alone checks out the frozen Kostka submodule
because a few retained GT experiments consume it for reproducibility. The
workflow also checks the declared Polytool library/web MSRV (Rust 1.82) and MCP
MSRV (Rust 1.88), and treats Rust documentation warnings as errors. Standalone
workspace lockfiles are committed at `experiments/Cargo.lock`,
`kostka/Cargo.lock`, and `KTT-search/Cargo.lock`.

## Polytool standalone branch

`master` is the canonical history for this workspace. The `main` branch is the
generated standalone view of `polytool/`; do not develop directly on it. After
committing and pushing eligible `polytool` work on `master`, publish it with:

```sh
./scripts/sync-polytool-main.sh
```

The script requires a clean `polytool/` tree and refuses to overwrite a
diverged `main` branch.
