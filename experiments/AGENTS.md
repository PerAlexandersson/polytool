# Experiments Agent Guide

## Purpose and persistence

This crate is a scratch surface for bounded research computations. Most files
under `src/bin/` are deliberately ignored and use-once. Do not force-add them,
catalogue the whole local forest, or treat an ignored path as a durable API.

Tracked experiment binaries are exceptional. Track one only when it is a
reproducibility artifact, an independent oracle for a maintained library, or a
durable integration check that cannot live in a standalone crate. Add a narrow
`.gitignore` exception and state the reason in the Rust handoff.

## Before writing an experiment

1. Search `../docs/LIBRARY_GUIDE.md` and existing APIs with `rg`.
2. Reuse `combinatoric-core`, `combpoly`, `sym-poly-*`, `polytool`, or `kostka`
   instead of copying enumeration, coefficient-vector, interpolation, graph,
   or exact-polynomial helpers.
3. Check `../docs/EXPERIMENT_PROMOTION.md` before extending a candidate already
   marked for promotion.
4. Keep the research driver separate from its mathematical kernel. A driver
   may print tables and scan parameter grids; a reusable kernel needs a named
   library API and tests.

`src/lib.rs` exists only to share plumbing among local binaries. It is not a
second general-purpose library. When a shared helper acquires multiple durable
consumers or a mathematical contract, promote it to the owning crate.

Use the checked pattern APIs in `combpoly::permutation` and state the pattern
convention explicitly. Do not add another handwritten triple inequality for a
named permutation pattern.

## Promotion threshold

Promote a kernel when at least one of these holds:

- it is the canonical generator for a named polynomial or sequence;
- two independent projects need it;
- several binaries have copied it or depend on a shared experiment helper;
- it implements a generally useful exact combinatorial structure;
- a website, paper, OEIS entry, CLI, or MCP result should cite it durably.

A promotion must include:

- a precise mathematical definition and target module;
- exact coefficient and indexing conventions, including the empty object;
- invalid-input and overflow behavior;
- focused unit tests and an independent small-instance oracle;
- one documented usage example;
- updated capability and promotion guides;
- migration of durable consumers without absorbing the exploratory scan loop.

Do not promote code merely because it is large. Failed conjecture probes,
parameter grids, manuscript-specific assertions, and report formatting normally
stay local and ignored.

## Computation and output

- Prefer exact arithmetic. A floating-point scan may suggest a conjecture but
  cannot certify it.
- Make expensive searches bounded, resumable when practical, and explicit
  about the tested range.
- Keep generated tables, checkpoints, database dumps, and build output out of
  Git and outside the source tree when practical.
- Never embed credentials. Load configured services through the established
  environment and do not print secrets.
- Use the external Cargo target directory and the workspace's bounded,
  reduced-priority command convention for newly written searches.

For a tracked experiment, record the exact command, bounds, source commit, and
result location in the relevant research project rather than growing this
guide into a chronological ledger.
