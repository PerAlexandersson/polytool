# OEIS catalog commands

Part of the [CLI reference](cli.md).

## Read rows from local OEIS data

`polytool oeis rows A001263 --data DIR [--json]` reads local `.seq` and
matching consecutive b-file terms without the catalog or network. It prints
complete coefficient rows, or a refusal reason. For `tabf`, the consecutive
`%e` table must start at the first term: earlier row boundaries cannot be
inferred safely. At least three example rows are required, with two observed
increases for periodic growth. JSON includes the existing display string and
`row_length_rule_details`: `kind` is `constant_step` or `periodic_step`,
`step` is the increase, `period` counts transitions per increase, and `phase`
is the zero-based transition index modulo the period. It also includes a
stable `reason_code`, human `reason`, exact decimal-string rows, term source,
and `bfile_rejection`, `bfile_rejection_reason` and optional I/O detail when
a present b-file cannot be used. Text output gives b-file rejection on stderr.

## Generate recurrence-backed OEIS families

The bundled OEIS catalog exposes one Rust function per A-number, for example
`polytool::oeis::A008292()`, together with dynamic lookup by identifier.  Each
entry is generated from a sparse exact recurrence description and its initial
rows; installed binaries do not need an OEIS database, Python, or SymPy.

```sh
polytool oeis list
polytool oeis info A008292 --json
polytool oeis generate A008292 --count 8
polytool oeis generate A008292 --count 8 --format triangle
polytool oeis generate A008292 --count 40 --format bfile --max-terms 1000
```

Generation formats are `rows`, `triangle`, `polynomial`, `json`, `jsonl`,
`csv`, and `bfile`.  `--from-index` selects a later displayed OEIS row.  Strict
b-file output is enabled only when the bundled complete-row prefix has been
matched against OEIS data; `--max-terms` stops before a row that would cross
the cap.  The same formatter writes b-files for any recurrence with
`polytool recurrence-generate --format bfile` (see `doc/cli-recurrence.md`).  The JSON returned by `oeis info` keeps the human-readable
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
| `polytool_recurrence_search` | 54 | `polytool oeis add` (below); records `fit_rows`, `checked_terms`, `search_options` and the polytool version |
| `supplied_recurrence` | 91 | `polytool oeis add --recurrence-json` (below): a supplied recurrence, validated against every OEIS term of the snapshot; records `checked_terms`, `terms_source`, the last compared row and the polytool version |

The first three tiers were produced by an earlier offline generator, which
combined the benchmark fixtures, a curated recurrence queue, and generated
Lean definitions with a local OEIS copy; its output was migrated to the data
file unchanged (a test compares every entry's definition and first 20 rows
with fingerprints captured before the migration).

The catalog contains 919 entries. The [November 2026 additions](oeis-catalog-additions-2026-11.md)
record the source revision, validation and conservative proof flags for 59
additional families.

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

### Supplied recurrences

`oeis add <A-number> --recurrence-json FILE` records a known recurrence instead
of searching for one.  `FILE` is a `polytool.recurrence.v1` object (the
`recurrence_data` member of `oeis info --json` has this form): the rule, the
recurrence index `first_index` of its first initial row, and the initial rows.
The first initial row is matched with the first OEIS row after the
`--prefix-rows P` complete OEIS rows (default 0), which are stored as
authoritative prefix rows.  `tabl` sequences are regular triangles; another
sequence needs `--layout triangle` or `--layout table` (with a supplied
recurrence a table needs no `--row-widths`: the generated rows are the table
rows, and `--row-widths` is needed only to cut prefix rows).  The candidate
entry is replayed through the catalog's row generator and compared term by term
with every OEIS term (the same comparison as `oeis check`); it is written only
when every term matches, at least five complete rows and two rows computed by
the recurrence are covered, and the status is `validated` (no rows are held out
because nothing was fitted).  A wrong recurrence gives the outcome `mismatch`
with the first differing term.  Search options and `--fit-rows` do not apply.

```sh
polytool oeis info A008292 --json | jq .recurrence_data > a008292.json
polytool oeis add A008292 --recurrence-json a008292.json --replace --dry-run
polytool oeis add A051142 --recurrence-json a051142.json --prefix-rows 2
```

As with the search, a finite match is evidence that the supplied rule produces
the OEIS rows, not a proof for every row.

### Real-rootedness flags

A record may carry the optional key `real_rootedness`, the last key of its
line:

```json
"real_rootedness":{"status":"lean_proved","as_of":"2026-10-03"}
```

| `status` | Meaning |
| --- | --- |
| `lean_proved` | Every row of the entry is zero or real-rooted, by a machine-checked Lean proof (only the axioms `propext`, `Classical.choice` and `Quot.sound`) of a definition whose rows equal the catalog rows for every index.  The equality is shown exactly for all `n`: after the index shift that aligns the rows, both recurrences have the same normal form (every coefficient divided by the common polynomial gcd, leader normalized), the catalog leader does not vanish on the computed rows, the Lean rule divides by no expression that vanishes, and the rows agree exactly up to the index from which both rules apply.  Leading catalog rows before the Lean definition's first row (at most three) were checked exactly instead. |
| `known_not_real_rooted` | Some row has a non-real root (an exact check of a generated row). |
| absent | Unknown. |

`as_of` (optional, `YYYY-MM-DD`) dates the proof-status snapshot or the check.
The flag is not a polytool proof: `oeis info` prints it as a recorded flag, the
JSON of `oeis info` and `oeis list --json` reports it as `real_rootedness`
(`null` when unknown), and the Veronese parent provers use a `lean_proved`
entry only as an attestation (verdict `conditional_on_attestation`, label
`oeis_catalog:<A-number>:lean_proved`), see the library documentation.  The
bundled catalog flags 795 entries `lean_proved` and one (A122431, whose row
`3 + 3x + x^2` is not real-rooted) `known_not_real_rooted`.  A replaced entry
(`oeis add --replace`) loses its flag, because the flag certifies the rows of
the recurrence it was established for.

### Evidence records

A record may carry the optional key `evidence`, written after `real_rootedness`
(and so last):

```json
"evidence":{"strongest":"log_concave","rows_checked":50,"first_failure":{"real_rooted":2,"ultra_log_concave":2},"as_of":"2026-11-02"}
```

It states which coefficient property held on **every** one of the first
`rows_checked` rows of the entry (consecutive rows from `first_row`, none
skipped, all checked exactly by polytool), and where stronger properties fail.
It is **evidence, not a proof**: it says nothing about later rows and is
separate from the `real_rootedness` flag, which certifies all rows.

Properties, strongest first, for a row with coefficients `a_0, ..., a_d`
(ascending powers, trailing zeros removed; the zero row satisfies all of them):

| Property | Meaning for a row |
| --- | --- |
| `real_rooted` | `is_real_rooted_bigint_coeffs`: all roots real (exact, any coefficient signs) |
| `ultra_log_concave` | admissible shape and `is_ultra_log_concave_bigint_coeffs` (Newton's inequalities at the nominal degree `d`) |
| `log_concave` | admissible shape and `is_log_concave_bigint_coeffs` (`a_k^2 >= a_{k-1} a_{k+1}`) |
| `unimodal` | admissible shape and `is_unimodal_bigint_coeffs` |
| `none` | (value of `strongest` only) no property holds on every checked row |

The three coefficient properties require an *admissible shape*: all
coefficients nonnegative and no zero between two nonzero coefficients (leading
and trailing zeros are allowed).  A row with a negative coefficient or an
internal zero therefore fails them, and for such a row only `real_rooted` can
hold; this keeps "log-concave" from holding by accident, for example for
`1 + t^3`, which satisfies `a_k^2 >= a_{k-1} a_{k+1}` pointwise.  On admissible
rows ultra-log-concave implies log-concave implies unimodal.

`strongest` is the first property in this order that holds on every checked
row, or `none`.  `first_failure` has one entry, the first failing row, for each
property stronger than `strongest`; the row is the *displayed OEIS row index*
(the indexing of `first_row` in `oeis info`, so row 2 of A122431 is
`3 + 3x + x^2`), and lies among the checked rows.  A real-rooted family with
mixed-sign rows has `strongest = real_rooted` and no failures recorded for the
weaker properties.  `rows_checked` is smaller than the requested count only
when the entry cannot generate more rows.

```sh
polytool oeis evidence A008292 A122431 --count 50            # print only
polytool oeis evidence --all --count 50 --write              # record every entry in the data file
polytool oeis evidence A001263 --count 30 --json
polytool oeis info A122431                                  # shows the record, labelled evidence
                                                            # ("none recorded" until one is written)
```

`oeis evidence` regenerates the rows from the catalog's own recurrence
(`--catalog FILE`, default `data/oeis-catalog.jsonl` in the current directory
when present, else the bundled catalog) and reports one line per entry
(`A-number`, `strongest`, `rows=N`, failing rows as `property@row`), the
distribution of `strongest`, and every `lean_proved` entry whose evidence is
weaker than `real_rooted`, which would contradict the flag.  `--write` replaces
the `evidence` key of each computed entry in the data file and nothing else:
it refuses to rewrite a file that is not in canonical form, and it checks that
every line with its evidence removed is unchanged and that records it did not
touch keep their bytes.  `--per-entry-timeout-seconds S` (default 60, 0 for no
limit) is checked before each row; an entry over it gets no new record (an
existing one is kept) and is listed as `timeout`.  `--as-of DATE` overrides the
recorded date (default: today, UTC).  `oeis add --replace` drops the entry's
flag and evidence, because both describe the replaced recurrence's rows.

Release packages omit repository-only agent notes, handoffs, CI configuration,
and benchmark output; they include the catalog data file, from which the build
script generates the catalog. The recurrence fixtures stay
in the package deliberately: the installed CLI's `bench recurrence-fixtures`
command uses them as its default self-contained corpus.
