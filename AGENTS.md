# Polytool contributor guide

## Standalone build

- The crate must build standalone.  Its manifest, tests, examples, build script,
  and published sources must not use path dependencies outside this directory
  (the in-tree `crates/`, `mcp/`, and `web/` packages are fine).
- Validate with `cargo metadata --locked --format-version 1 --no-deps`.
- Keep build output (`target/`, `pkg/`, and similar) out of the source tree,
  for example by setting `CARGO_TARGET_DIR`.

## Layout and where new code goes

- `src/lib.rs` documents every public module and re-exports the main entry
  points; keep its module list in sync when a module is added, renamed or
  removed.
- Module map (library crate):
  - polynomial types and arithmetic: `polynomial` (`Polynomial<C>`, `CoeffRing`),
    `vec_poly` (arithmetic on plain `&[C]` coefficient vectors), `numeric`
    (binomials, factorials, q-analogues), `parse` (text to coefficients);
  - real-rootedness, interlacing, log-concavity, gamma-positivity, resultants,
    Ehrhart/h* conversion: `real_rootedness`; its default exact backend is the
    separate `exact-real-roots` crate (`crates/exact-real-roots`), re-exported as
    `root_count`; `sturm` and `sturm_cf` hold the Sturm-chain routes;
    `real_rooted_enumeration` enumerates real-rooted integer polynomials with
    a prescribed leading prefix (exact antiderivative-shift kernel);
  - linear algebra, total non-negativity, determinants: `linalg`;
  - recurrence search and renewal recognition: `recurrence`, `renewal`; bundled
    OEIS families: `oeis` (data in `data/oeis-catalog.jsonl`, generated into
    `$OUT_DIR` by `build.rs`; record schema in `src/oeis_catalog_schema.rs`);
  - structured diagnostics: `coefficient_tests`, `hstar_inequalities`,
    `cyclic_sieving`, `decomposition`, `basis`, `bkw`;
  - planar-network and TNN certificates: `tnn_network`, `brenti_sequence`,
    `brenti_recognition`, `interlacing_matrix`; Pólya frequency certificates
    for rows: `pf_network` (height-invariant networks), `left_production`
    (Chen--Fu--Ruan left production matrices);
  - standard sequences: `sequences` (with `deco` and `literature`).
- The CLI is `src/bin/polytool/`: `main.rs` dispatches, each subcommand (or
  close group) has its own module, and `input`, `json_out` and `help` hold the
  shared parsing, JSON fragments and help text.  `mcp/` and `web/` are adapters
  over the library, not separate implementations.  Put mathematics in the library
  and keep adapters thin.
- `recurrence.rs`, `real_rootedness.rs` and `linalg.rs` are facades: their
  private submodules (under `src/recurrence/`, `src/real_rootedness/`,
  `src/linalg/`) are split by concern and re-exported, so public paths stay
  `polytool::recurrence::*` etc. Add new code to the matching submodule and
  re-export it from the facade. `src/recurrence/` is grouped one level deeper
  (`base`, `fit`, `transform`, `implicit`, `root_domain`, `representation`);
  each group has a `<group>.rs` facade that re-exports its members, and a
  cross-group helper is `pub(in crate::recurrence)`. The MCP crate (`mcp/src/`)
  is split the same way (requests, responses, parsing, server, tests).
- Reuse before adding: look for an existing helper in `vec_poly`, `numeric` and
  `real_rootedness` before writing a new coefficient-vector routine, and build
  new real-rootedness checks on `is_real_rooted*` or `exact-real-roots` rather
  than on floating-point root finders.
- `data/oeis-catalog.jsonl` is the single source of the OEIS catalog.  Add or
  replace entries with `polytool oeis add` (it writes canonical, sorted lines
  with provenance) rather than by hand; never change an existing entry's
  recurrence without a reason recorded in its provenance, and re-bless
  `tests/fixtures/oeis-catalog-fingerprints.tsv` only for intentional changes
  (`POLYTOOL_BLESS_OEIS_FINGERPRINTS=1`).  `src/oeis_catalog_schema.rs` is
  included by both `build.rs` and `src/oeis.rs`, so it has no inner
  attributes or `use` items.

## Mathematical and API conventions

- Preserve exact arithmetic.  Do not replace exact PRS, Bézout, rational, or
  `BigInt` routines with floating-point heuristics.  Use the documented `BigInt`
  APIs when `i64` coefficients can overflow, and keep explicit overflow and
  shape contracts.
- Coefficient vectors are in ascending degree order: entry `i` is the coefficient
  of `t^i`.  Preserve the documented directed-degree and zero-polynomial
  conventions, and add focused tests whenever behavior changes.
- Zero-polynomial representation differs by layer: `vec_poly` arithmetic returns
  `[0]` (use `vec_poly::trim`), while `Polynomial::coeffs()` and
  `vec_poly::strip_trailing_zeros` give the empty vector.  Parsers do not trim.
- `i64` entry points are convenience wrappers.  Generic `Polynomial<i64>` and
  `vec_poly` arithmetic panics on overflow (checked, in release builds too);
  `*_bigint` / `*_bigint_coeffs` functions never overflow.  Document any panic in
  a `# Panics` section and prefer returning `Result`/`Option` for invalid input.
- Interlacing functions return `Option<bool>`: `None` means the degrees do not
  admit that directed test, not "false".  `check_weak_interlacing(f, g)` tests
  `f << g`; use `interlaces_either_order` for the undirected question.
- Diagnostic report types (`*Report`, `*Check`, `*Error`) are serialized by the
  CLI and MCP server.  Document every public field (meaning, indexing base,
  coefficient order, units), and treat renaming or reordering fields as a schema
  change.
- A passing finite check (coefficient tests, h* inequalities, TNN of a finite
  truncation, a floating-point BKW scout) is evidence, not a proof; say so in the
  rustdoc of such functions.

## CLI and documentation

- Keep CLI help detailed and AI-friendly.  Machine-readable commands must retain
  stable JSON schemas; serialize exact integers as strings where the JSON number
  model could lose precision.
- Update the documentation and focused tests in the same change when public API,
  CLI input/output, exactness guarantees, or standalone build behavior changes.
  `README.md` is a short landing page; the detailed material lives in `doc/`
  (`cli*.md` for the CLI reference, `library.md`, `mcp.md` and `mcp/README.md`,
  `web.md`, `development.md`, `algorithms.md`) and in `CHANGELOG.md`.  Put CLI
  details in the matching `doc/cli*.md` page and release notes in `CHANGELOG.md`.
- Every module starts with `//!` documentation (purpose, entry points,
  conventions, related modules).  Public items need rustdoc; check with
  `cargo rustdoc -q --lib -p polytool -- -W missing_docs` (it reports none;
  keep it that way, including in the generated OEIS catalog, whose generator
  emits accessor docs).  Doc examples must compile as doctests and
  must not use APIs newer than the declared `rust-version` (1.82).
- This directory is published as a standalone repository: do not put local
  paths, private project names or workflow notes into code, comments, docs or
  fixtures.

## Verification

- Run focused library and CLI tests before broader checks.
- Before finishing a change run `cargo fmt --all -- --check`,
  `cargo clippy --locked -p polytool -p exact-real-roots --all-targets -- -D warnings`,
  `cargo test --locked -p polytool -p exact-real-roots --doc`,
  `RUSTDOCFLAGS="-D warnings" cargo doc --locked --no-deps -p polytool -p exact-real-roots`
  and `cargo +1.82 check --locked -p polytool -p exact-real-roots`.
- `cargo test --lib oeis::tests` checks the catalog data file (canonical lines,
  provenance tiers, recorded fingerprints) without external data.  With a
  local OEIS copy, `polytool oeis check --oeis-data DIR` re-validates every
  entry term by term.
- Verify that `cargo package --list` excludes agent notes, CI configuration,
  and benchmark output, while retaining `data/oeis-catalog.jsonl` and the
  recurrence fixtures used by the installed `bench recurrence-fixtures` command.
- Preserve unrelated changes when editing.
