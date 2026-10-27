# OEIS catalog commands

Part of the [CLI reference](cli.md).

## Generate recurrence-backed OEIS families

The bundled OEIS catalog exposes one Rust function per A-number, for example
`polytool::oeis::A008292()`, together with dynamic lookup by identifier.  Each
entry is generated from a sparse exact recurrence description and its initial
rows; installed binaries do not need an OEIS database, Python, or SymPy.

```sh
polytool oeis list
polytool oeis info A008292 --json
polytool oeis generate A008292 --rows 8
polytool oeis generate A008292 --rows 8 --format triangle
polytool oeis generate A008292 --rows 40 --format bfile --max-terms 1000
```

Generation formats are `rows`, `triangle`, `polynomial`, `json`, `jsonl`,
`csv`, and `bfile`.  `--start-row` selects a later displayed OEIS row.  Strict
b-file output is enabled only when the bundled complete-row prefix has been
matched against OEIS data; `--max-terms` stops before a row that would cross
the cap.  The JSON returned by `oeis info` keeps the human-readable
`recurrence` string and also includes `recurrence_data`, a directly consumable
`polytool.recurrence.v1` object with exact rational coefficients and initial
polynomials.  Catalog entries do not fabricate recurrence-search diagnostics,
so this embedded object omits the optional `search` member.  The catalog
distinguishes holdout-`verified` recurrences from
`validated` recurrences whose generated rows match the current OEIS prefix but
whose original fitting/holdout provenance is unavailable.  Both are enabled by
default. The original generator rejected imported queue and Lean definitions
that could not be aligned with complete OEIS rows; in particular, the catalog
does not publish a recurrence fitted to a row fragment or to coefficient data
with positional zeros removed. The `experimental` status and `--include-experimental` option
remain available for explicitly curated future entries, but the current
bundled catalog contains none.

### Catalog data, provenance, and adding sequences

The catalog's single source of truth is the tracked data file
`data/oeis-catalog.jsonl`: one JSON object per sequence, sorted by A-number,
with a fixed key order so that diffs stay small.  The build script
(`build.rs`) validates it and generates the static catalog, the accessor
functions and their rustdoc titles; no external data or tooling is needed to
build.  The record schema is public as `polytool::oeis::OeisCatalogRecord`.
Each record holds the executable recurrence (`layout`, `first_row`,
`flattened_offset`, `bfile_prefix_verified`, `prefix_rows`,
`recurrence_first_index`, `recurrence_first_width`, `initial_rows`,
`recurrence`, with recurrence coefficients as sparse `[n_degree, t_degree,
"value"]` monomials) and a `provenance` object with the `status`, the
`source_rows` / `verification_rows` counts and the `rows_sha256` fingerprint.
The `source` field records where the recurrence came from:

| `source` | Entries | Origin |
| --- | --- | --- |
| `benchmark_fixture` | 73 | `polytool recurrence --json` on `fixtures/recurrence-benchmarks` (held-out verified); also records `fixture_slug` |
| `lean_definition` | 630 | Fitted by an earlier version of the polytool recurrence search and transcribed into generated Lean definitions; validated against the OEIS data prefix (16 rows); records the last compared row as `validation_row` |
| `curated_queue` | 71 | Fitted by recurrence search on database-cached rows; `verified` when checked on held-out rows |
| `polytool_recurrence_search` | new | `polytool oeis add` (below); records `fit_rows`, `checked_terms`, `search_options` and the polytool version |

The first three tiers were produced by an earlier offline generator, which
combined the benchmark fixtures, a curated recurrence queue, and generated
Lean definitions with a local OEIS copy; its output was migrated to the data
file unchanged (a test compares every entry's definition and first 20 rows
with fingerprints captured before the migration).

To add sequences, run from the source directory with a local copy of the OEIS
in its internal format (`seq/A123/A123456.seq` files):

```sh
polytool oeis add A008275 --oeis-data ~/oeis --dry-run   # print the entry only
polytool oeis add A008275 A001497 --oeis-data ~/oeis      # write data/oeis-catalog.jsonl
polytool oeis add --ids-file ids.txt --max-candidates 20000
polytool oeis add A055151 --layout table --row-widths 1+1/2 --replace
polytool oeis add A123456 --bfile b123456.txt --name "..." --layout triangle
polytool oeis check --oeis-data ~/oeis                    # re-validate every entry
```

`oeis add` reads the name, offset, keywords and terms (or an explicit b-file),
cuts the terms into complete rows (regular triangles for `tabl`; any other
sequence needs `--layout triangle` or `--layout table --row-widths W+S[/R]`),
runs the adaptive recurrence search with the options of `polytool recurrence`
on a prefix of the rows (`--fit-rows`, by default all but `max(2, rows/4)`),
and then replays the candidate entry through the catalog's own row generator
and compares it with every available term, held-out rows and a trailing
partial row included.  The status is `verified` when held-out rows were
checked and `validated` otherwise.  Each A-number gets one outcome: `added`,
`updated` (with `--replace`), `already_present`, `no_recurrence`,
`layout_unsupported`, `insufficient_data`, `mismatch`, or `unavailable`.
Rebuild polytool to bundle new entries.  `oeis check` regenerates entries
(bundled, or from `--catalog FILE`) and compares them term by term with the
local OEIS data, listing entries it could not check.  A finite match is
evidence, not a proof that the recurrence holds for every row.

Release packages omit repository-only agent notes, handoffs, CI configuration,
and benchmark output; they include the catalog data file, from which the build
script generates the catalog. The recurrence fixtures stay
in the package deliberately: the installed CLI's `bench recurrence-fixtures`
command uses them as its default self-contained corpus.
