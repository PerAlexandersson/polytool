# Curated sequence families

The family catalog combines every bundled OEIS family with curated sequences
that are candidates for OEIS submission. Its scope is one polynomial row per
integer index. Larger classes indexed by arbitrary graphs, matroids or shapes
belong in the external open-problems catalog. An automated atlas is a discovery
census; promotion into this catalog adds a curated definition and provenance.

```sh
polytool families list --query packing
polytool families list --status theorem_proved --json
polytool families show PT-F0001
polytool families generate PT-F0001 --count 9 --format rows
polytool families generate A035469 --count 6 --format terms
polytool families validate --json
polytool families export-oeis PT-F0002 --count 40 --bfile packing-cap5.txt
polytool families export-oeis PT-F0001 --json
```

`list` emits compact summaries. `show` emits the complete record and generator
as JSON. `generate` emits complete coefficient rows in ascending powers of `t`,
flattened comma-separated `terms`, or JSON with exact integer strings and the
first displayed row index. `--count` counts complete rows, including the first;
it is bounded to 1–512. Exponential-family ranks are bounded to 1000.
No root check is implied by generating rows. Use `exp-rows` for exact finite
root and interlacing checks of an exponential generator.

`export-oeis` prints internal-format draft blocks and a proposal TOML skeleton
matching [oeis-contributions](https://github.com/PerAlexandersson/oeis-contributions).
JSON additionally contains the b-file text. `--bfile` creates a new file and
refuses to overwrite an existing one. Data fields are generated from the same
canonical rows as `generate`: at most three `%S/%T/%U` lines, each at most 80
bytes including the tag and A-number, and at most 240 data payload bytes in
total. Integers are never split. B-files include up to 10,000 terms or 1 MiB
from the requested complete rows; `bfile_truncated` reports either cap. To
obtain longer b-files request more rows. `A000000` is a placeholder for an
unassigned A-number. The second `%O` component is inferred from the generated
prefix and needs review. Keywords describe the prefix only.

The proposal is deliberately a skeleton: fill in the author, check novelty
again, and add a repository-local generator adapter that reproduces its
**bounded data prefix**, not the entire longer b-file. Its `data.terms` field
is the exact `%S/%T/%U` prefix, at most 240 characters. The adapter command is
included in a comment and `%o`. Exporting does not submit or modify OEIS or
the proposal repository. Existing A-number exports use a `comment` proposal
kind for review rather than claiming to be new sequences.

## Metadata and validation

`data/families/*.json` contains strict version-1 sidecars keyed by `A######`
or `PT-F####`; the filename must equal the ID. The Rust serde schema in
`src/families/schema.rs` is authoritative. Unknown fields, invalid enum values,
invalid calendar dates, unresolved generator IDs, mismatched OEIS layout or
indices, duplicate IDs/reference keys, dangling theorem/related references,
incorrect initial rows and inconsistent submission bookkeeping are rejected.
`families validate` generates at least six rows per record and compares every
recorded initial row. This does not prove the interpretation of a generator.
The test suite and documentation build both run this validator.

The fields are `schema_version`, `id`, `title`, `definition`, `formula`,
`generating_function`, `generator`, `first_row`, `flattened_offset`, `layout`
(`tabl` or `tabf`), `degree_pattern`, `initial_rows`, `references`,
`real_rootedness`, `interlacing`, `models`, `related`, `oeis`, `submission`,
`links`, and `generator_status`. Exact initial coefficients are integer
strings in ascending powers; OEIS structural trailing zeros are preserved.
Empty curated strings and lists mean “not curated”, not an assertion.

Generators have exactly one tagged form:

```json
{"kind":"oeis","id":"A035469"}
```

```json
{"kind":"exp_rows","family":{"unmarked":{"kind":"sizes","sizes":[1]},"marked":{"kind":"finite","weights":["0","1","4","31","125"]}}}
```

The third form is `{"kind":"recurrence","record":...}` containing a complete
`polytool.recurrence.v1` object. Its first generated recurrence row is displayed
as the family's `first_row`; subsequent rows advance both indices by one.
The exponential generator's `first_row` is its rank in `n![z^n]exp(A+tB)`.
Transfer matrices are outside phase 1.

All OEIS records are derived at runtime; uncurated entries need no sidecar.
Their generator validation status and root attestation (`status`, `as_of`,
`statement` from the audit) always come from `data/oeis-catalog.jsonl`.
Those fields in an OEIS sidecar are placeholders, not another source of truth.
Exact finite evidence is derived where available. Curated theorem keys,
references and public Lean links supplement the canonical attestation.

Properties use the existing `lean_proved`, `theorem_proved` and
`known_not_real_rooted` flags, plus `unknown`, `conjectured` and
`finite_evidence`. Each property may record `as_of`, a local `theorem` reference
key, `statement`, `lean`, `paper_proof`, `finite_range` and positive `lag`.
A finite range is inclusive in displayed indices; for lag-one interlacing it
indexes the **left** rows. Reported external evidence is labeled as such and
is not an all-index proof. A model has `family`, JSON-object `parameters`, and
`replay_status` (`pending`, `passed`, `failed`). A related entry has an `id` and
precise `relation`. References carry `key`, `citation`, optional `doi`/`arxiv`,
`locator` and public `url`.

`oeis.id` is an A-number or JSON `null`. Non-OEIS entries require `checked`
(`YYYY-MM-DD`) and nonempty `queries`. A stable PT ID remains after assignment;
update `oeis.id` and `submission.a_number`. Submission states are
`not_started`, `draft`, `ready`, `submitted`, `accepted`, `withdrawn`;
submitted/accepted records require `date`, and accepted records require a
consistent A-number. External `links` point one way only.

## Public Lean mappings

Manual links currently allow only the public RealRooted repository. Each has
`repository`, `declaration`, a full lowercase 40-hex `commit`, and a source
`permalink` pinned to that commit. Generic theorem links do not assert that the
catalog identification was checked in Lean. Older canonical `lean_proved`
flags without a curated public mapping keep their status and have no link.
Private proof repositories are never linked.

The proving repository will eventually own the family tags and export the
mapping for links in both directions. This phase does not add tags to another
repository or invent reverse links. An optional hook already accepts a local
JSON array exported at pinned commits:

```json
[{"family_id":"A143494","property":"real_rootedness","repository":"https://github.com/PerAlexandersson/RealRooted","commit":"28fa02e98eedb2e351ff7aa203cbf0f8181668f0","declaration":"RealRooted.DistinctValues.poissonPolynomial_isPFPolynomial"}]
```

`families validate --lean-export declarations.json` requires every manual
mapping to occur in the supplied list, including its family, property and
commit. An empty list rejects all mapped declarations. The caller supplies a
trusted proof-owner export; this hook does not verify Lean or fetch the network.
The docs builder accepts the same `--lean-export` option. By default the build
checks link shape and public allowlist only; declaration existence is manually
attested, and later renames are isolated by commit pinning.

## Library and MCP

`polytool::families::{catalog, parse_family, validate_catalog,
validate_lean_declarations, export_oeis}` and `FamilyRecord::generate_rows`
share the same data and typed `FamilyCatalogError`. The catalog is loaded
without filesystem access and ships in the published crate.

MCP tools mirror the CLI: `list_families`, `get_family`, `generate_family_rows`,
`export_family_oeis`. For example:

```json
{"id":"PT-F0001","count":9}
```

Use that request for generation or export. `get_family` takes only `id`;
`list_families` accepts `query`, `status`, `offset`, and `limit` (1–100, default
50). It returns `total` and `next_offset`. MCP generation/export counts are
bounded to 1–64 and their serialized responses to 1 MiB. No tool writes files
or posts submissions. CLI JSON envelopes are `polytool.families.v1`,
`polytool.family.v1`, `polytool.family-rows.v1`,
`polytool.families-validation.v1`, `polytool.family-oeis-export.v1`.

The [documentation builder](development.md#documentation-site) produces a
searchable/status-filterable family index, one page per family, generated
first rows, copy buttons, b-files, proposal skeletons and validated metadata.
All generation runs through the built CLI; Python only renders the result.
