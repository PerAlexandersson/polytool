# Polytool handoff

## Web example catalogue (2026-09-02)

The web UI presents one flat menu of 34 distinct OEIS-labelled polynomial
families. It absorbs the former eight example buttons, removes the four
families duplicated between those buttons and the 30-entry OEIS menu, and
removes the category groups.

All 30 imported OEIS entries have recursive definitions in the corresponding
`ProofsOeis/A*.lean` files. The four additional distinct built-in families
(derangement excedances, Fibonacci matchings, Touchard polynomials, and Simsun
descents) also carry explicit recurrences. This does not imply that the web
recurrence search will rediscover every formula under its default bounds.

Long polynomial previews now parse and rejoin term signs before inserting the
ellipsis, so positive terms no longer render as `+ +` and a negative final term
renders with `-` rather than `+ -`.

Verification:

```text
inline JavaScript parse                                      passed
flat picker count, uniqueness, and key resolution           34/34
loadExample textarea/status checks                          34/34
recurrence provenance audit                                 34/34
positive/negative abbreviation regression cases              passed
desktop and 390px-wide headless-Chrome inspection            passed
standalone wasm-pack release build                            passed
git diff --check                                              passed
```

Commit `108916e` was pushed to monorepo `master`, projected to standalone
`main` as `fe697ec`, and deployed to `poly.symmetricfunctions.com` on
2026-09-02. Cache-busted public fetches matched the staged HTML, JavaScript,
and WASM byte for byte; a live headless-Chrome load reported no WASM or
JavaScript error.

The private proof-repository name was subsequently removed from all 26
imported-sequence comments before publication. The public examples now retain
only their mathematical descriptions and OEIS URLs; they do not mention the
author's related projects.

The recurrence-result page no longer generates or displays the bulky standalone
Python export. It keeps the compact recurrence JSON and links directly to the
documented `polytool recurrence-generate` command, which reconstructs exact
rows from that JSON. The library and CLI Python exporter remain available for
backward compatibility; only the browser payload and UI were simplified.

The browser's adaptive bounds are now recurrence depth 10, `t`-degree 5,
`n`-degree 5, and derivative order 5. In particular, this includes A059427's
cubic derivative coefficient `t - t^3`; a regression test confirms that degree
two fails and degree three finds the recurrence. The adaptive-mode tooltip
states that there is no elapsed-time or candidate-count cutoff: the search
stops only on a match, user cancellation, or exhaustion of its finite bounds,
and a failed exhaustive search may therefore take a long time.

Recurrence-option tooltips now mark their formulas with `data-tex` and render
them through the already loaded KaTeX runtime. Rendering is applied both to the
original tooltip nodes and to the floating tooltip layer; the alternating-sign
label and tooltip therefore display `(-1)^n` with an actual superscript.
The two focused web-crate tests, standalone WASM build, JavaScript parse,
code-card checks, and a headless-browser KaTeX/WASM load all pass.
Commit `1e8b865` was pushed to monorepo `master`, projected to standalone
`main` as `931a4e1`, and deployed on 2026-09-02. Cache-busted public files
matched the staged bundle byte for byte, and the live browser check passed.

A166073 was removed from the example picker because adaptive recurrence search
had to explore too much of the enlarged search space before finding its more
complicated recurrence. It was replaced by A008288, the Delannoy array read by
antidiagonals, whose row polynomials satisfy
`P(n) = (1 + t) P(n-1) + t P(n-2)`. A focused browser-crate regression test
checks that the full adaptive configuration finds this recurrence within ten
candidates. All remaining `Related project` comments were removed at the same
time. Commit `ccc9bf1` was pushed to monorepo `master`, projected to standalone
`main` as `d9f7361`, and deployed on 2026-09-02. Cache-busted public HTML,
JavaScript, and WASM matched the staged bundle byte for byte; a live
headless-Chrome load also passed.

## Uspensky/Descartes comparison

The main Rust worker owns these files:

- `src/root_count.rs`
- `src/lib.rs`
- `examples/bench_positive_real_rooted.rs`
- `README.md`
- `HANDOFF.md`

The exact Uspensky/Descartes comparison path is implemented and public. It
uses a strict Fujiwara bound, dyadic magnitude bands, reciprocal reduction
below `1`, homographic subdivision, and exact `BigInt` sign variations.  It
uses no finite fields or floating point.

The default one-signed real-rootedness path is now adaptive. It runs the cheap
Kurtz/Newton filters first, uses PRS generally, and selects Uspensky only when:

- the degree after stripping zero roots is at least 35;
- the two endpoint coefficients are equal; and
- an interior coefficient is at least `4^degree` times an endpoint.

The last two exact tests are invariant under scalar multiplication. This
conservative signature was chosen because degree or palindromicity alone
regressed important families. Square-free preprocessing is now shared with the
selected counter, avoiding a duplicate polynomial GCD in both the adaptive and
explicit PRS real-rootedness paths.

Release benchmark highlights from 2026-08-17 (single-process runs, so treat
sub-millisecond differences as noise):

```text
family                              primitive PRS    Uspensky
prod_{a=1}^{30} (x+a)                   0.88 ms       2.61 ms
prod_{a=1}^{80} (x+a)                  15.14 ms      82.43 ms
Eulerian (degree 35)                   45.82 ms      39.18 ms
Eulerian (degree 79)                    9.83 s        5.15 s
Narayana (degree 40)                    1.24 ms       9.41 ms
type-B Eulerian (degree 40)           182.37 ms     151.77 ms
Touchard (degree 40)                   55.54 ms      45.80 ms
Chebyshev T (degree 40)                 0.43 ms      61.05 ms
Chebyshev U (degree 40)                 0.43 ms      36.03 ms
Hermite (degree 40)                     0.33 ms      10.46 ms
```

These fair timings include the shared square-free preprocessing only once.
Eulerian, type-B Eulerian, and Touchard cross over around degree 35, while
Narayana, Chebyshev, Hermite, and evenly spaced linear-factor products strongly
favor PRS in the tested range. The benchmark example now covers all of these
families; Touchard is generated locally by its standard recurrence.

Verification:

```text
cargo test -q -p polytool --lib                 296 passed
cargo test -q -p polytool --test cli_bigint      16 passed
cargo test -q -p polytool --test interlacing_api  5 passed
cargo test -q -p polytool --doc                   5 passed
cargo clippy -q -p polytool --lib --examples -- \
  -D warnings -A clippy::manual-is-multiple-of \
  -A clippy::needless-range-loop                  passed
```

The two allowed Clippy lints are pre-existing in `cyclic_sieving.rs` and
`hstar_inequalities.rs`; neither file is part of this work.

## Repository ownership

- `/workspace/rust` branch `master` is the canonical monorepo history.
- All polytool changes are made under `polytool/` on `master`.
- Repository branch `main` is a generated standalone projection of this
  directory. Do not commit to it directly.
- After tested polytool changes have been committed and pushed on `master`, run
  `./scripts/sync-polytool-main.sh` from the monorepo root.

## Current state

The former independently rooted standalone history was backed up before
`main` was replaced by a `git subtree split --prefix=polytool` projection.
The projected tree includes the crate rename to `polytool` and all correctness
fixes from the August 2026 Rust review.
