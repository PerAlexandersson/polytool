# Rust combinatorics workspace

This is a multi-crate Cargo workspace for exact combinatorics, polynomial
experiments, and checked `symmetricfunctions.com` examples. Read the root
workspace guide first, then any deeper `AGENTS.md` before changing a crate.
For crate selection, public entry points, arithmetic, and naming distinctions,
read [docs/LIBRARY_GUIDE.md](docs/LIBRARY_GUIDE.md).

## Before editing

1. Check `git status --short`, the project `HANDOFF.md` when present, and
   current ownership. Preserve unrelated edits and the deliberate ignored
   experiment area; never clean it broadly.
2. Use `cargo metadata --no-deps --format-version 1` and `rg` to confirm the
   workspace member, dependency direction, existing public symbol, tests, and
   examples. Do this before adding a parallel implementation.
3. Work in the narrowest owning crate. Coordinate before changing a shared
   crate or a file another worker owns. Record durable work in that project's
   handoff unless a higher-priority task explicitly assigns it elsewhere.

## Dependency direction and ownership

The canonical foundations point upward; do not reverse these dependencies or
recreate their shared types:

```text
combinatoric-core
  └─ sym-poly-core
       ├─ sym-poly-multipoly ─┐
       └────────────────────── sym-poly-sym
                                └─ sym-poly-qsym
```

- `combinatoric-core` owns foundational combinatorics, including `Partition`,
  `Composition`, `Graph`, `Poset`, and `Ring`.
- `sym-poly-core` depends on and re-exports those foundation types; it owns
  shared algebra/tableau infrastructure. It is not a dependency of
  `combinatoric-core`.
- `sym-poly-multipoly`, `sym-poly-sym`, and `sym-poly-qsym` own the respective
  polynomial/function-algebra layers. `sym-poly-sym` uses `multipoly`, and
  `sym-poly-qsym` uses `sym-poly-sym`.
- `polytool` and `combinatoric-core` are separate maintained library roots.
  `combpoly` and `polynomial-lab` consume `polytool`; their MCP/web siblings
  are adapters, not duplicate libraries. The `kostka` submodule and standalone
  `KTT-search` application are retired reproducibility artifacts, excluded from
  the root workspace and routine CI. Route new GT/Ehrhart work to Ehrcalc and
  LR/Kostka/Schur compatibility work to lrcalc-rs. `flagged-lorentzian`
  consumes `sym-poly-core`.
- `experiments` is a standalone consumer workspace for the reusable crates. It
  is intentionally a disposable, often ignored staging area, not a catalog to
  complete or clean. Its ignored local binaries must not enter the maintained
  root target graph.

## Implementation conventions

- Reuse the canonical type and algorithm in the owning crate. New reusable
  sequence families need a documented, tested library generator; a one-off
  investigation may remain an ignored experiment. Promote only stable cores
  and selectively useful, tracked examples.
- Keep exact arithmetic exact. Prefer `BigInt`/`BigRational` when bounds can
  outgrow `i64`; do not silently truncate. Follow each API's coefficient order
  and conversion semantics.
- Add a small focused test and, for a user-facing mathematical feature, a
  checked example in the owning crate when it clarifies use. Do not treat a
  successful finite experiment as a proof.
- Use descriptive, domain-specific public names. Search for an existing name
  and distinguish similarly named APIs by crate and return type.
- Keep generated artifacts out of the synced checkout. In Docker set
  `CARGO_TARGET_DIR=/cargo-target/ai-projects`; do not retain a workspace-local
  `target/` tree.

## Verification and collaboration

- Run the narrowest relevant formatter, test, or Cargo check from the owning
  crate/workspace. For a new potentially unbounded executable, start with
  `timeout 60s nice -n 10 cargo ...`; use a proportionate monitored timeout
  for established Cargo work.
- Root Cargo commands cover maintained members. For a deliberate experiment or
  historical reproduction, pass the explicit `experiments/Cargo.toml`,
  `kostka/Cargo.toml`, or `KTT-search/Cargo.toml` manifest and preserve its
  committed standalone lockfile. CI tests tracked experiments only from a clean
  checkout; that job checks out the frozen Kostka submodule only because a few
  retained GT experiments still consume it. There is no dedicated Kostka/KTT
  job, and CI does not treat the ignored local forest as a maintained suite.
- Stage and commit only task-owned files after relevant checks; do not absorb
  generated files or another worker's changes. Follow the workspace's normal
  checkpoint/push policy unless the task says otherwise.
- The workspace `master` branch is canonical. `polytool/main` is a generated
  standalone subtree view; update it only with `scripts/sync-polytool-main.sh`
  after the documented prerequisites are met.
