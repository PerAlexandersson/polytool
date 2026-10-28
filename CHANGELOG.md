# Changelog

## Unreleased

This release is `polytool` 0.6.0 and `exact-real-roots` 0.2.0.  Existing CLI
commands and flags, MCP tool names and the web export names are unchanged; new
commands and tools are listed under Additions.  Some MCP responses now give
coefficients as decimal strings (see Breaking changes).

Breaking changes:

- MCP: no tool limits polynomial input to 64-bit integers any more.  The
  `polynomials` items of `parse_polynomials`, `real_roots`,
  `analyze_decomposition`, `check_polynomial_family`, `find_recurrence` and
  `find_recurrence_representations` take integers or exact decimal strings
  (the `PolynomialInput` and `PolynomialBatchInput` schemas, and the
  `polytool-mcp` re-exports `PolynomialInput`, `PolynomialBatchInput` and
  `NormalizedPolynomial`, are gone); their `expressions` and `text` forms are
  parsed exactly, and bare JSON integers up to `u64::MAX` in the rational
  `coefficients` rows are read exactly.  Response coefficients are decimal
  strings: `items[].coefficients` of `parse_polynomials`,
  `items[].polynomial.coefficients` of `real_roots`, every coefficient vector of
  `analyze_decomposition`, and `items[].coefficients`,
  `items[].gamma_coefficients` and `lace.matrix` of `check_polynomial_family`.
  `check_polynomial_family` with `sequence` no longer rejects rows beyond
  64 bits, and `find_recurrence_representations` names the rows that failed to
  parse.
- `cyclic_sieving::cyclic_sieving_sequence_reports_bigint` returns
  `Result<Vec<CyclicSievingSequenceItem>, CyclicSievingError>` (it panicked or
  silently skipped candidates); `CyclicSievingError` gains the variants
  `OrderTooLarge` and `SequenceIndexOverflow`, and
  `RecurrenceError` gains `Context`.
- `tnn_network::CanonicalTotallyNonnegativeProof::evaluated_path_matrix` returns
  `Option<Vec<Vec<BigRational>>>` (`None` for a malformed network).
- `sequences::literature::hoggatt_polynomials_bigint` now follows the
  `rows[n] = H_n` rule of the other `sequences` generators: it returns
  `max_n + 1` rows for `n = 0..=max_n`, and `rows[0]` is the zero polynomial
  `[0]` (it returned `max_n` rows with `rows[i] = H_{i+1}`, and no rows for
  `max_n = 0`).  Use `hoggatt_polynomials_bigint(n, m, q)[1..]` for the old
  result.  The `literature_sequence_rows hoggatt` example prints the extra row.
- Typed errors replace `String` errors; each `Display` text is the former
  message (checked by tests), and each type implements `std::error::Error`:
  - `parse::parse_polynomial`, `parse_polynomial_bigint`, `parse_polynomials`
    and `parse_polynomials_bigint` return `ParsePolynomialError` (new,
    re-exported at the crate root);
  - the `cyclic_sieving` functions that returned `Result<_, String>` return
    `CyclicSievingError` (new, re-exported at the crate root);
  - `linalg::check_tnn_neville` and `check_tnn_neville_bigint` return the
    existing `TotalPositivityError`, which gains the variants
    `InconsistentRowLengths`, `NegativeEntry`, `PositiveEntryBelowZeroPivot`
    and `NegativeEliminationEntry` for the Neville failures;
  - the `interlacing_matrix::check_lace_*` functions return
    `InterlacingMatrixError`, which gains the variants
    `NotTotallyNonnegative(TotalPositivityError)` and `SizeOverflow`;
  - `oeis::OeisCatalogRecord::{parse_jsonl, validate, validate_all}` return
    `oeis::OeisCatalogError` (new, one variant per rejected condition);
  - `renewal::{WeightedRenewal::new, recognize_weighted_renewal,
    recognize_fixed_delay_predecessor}`,
    `VariablePredecessorRecurrence::{certify, certify_all_indices,
    verify_certificate, rows}` and `RecurrenceTheoremReport::for_supplied`
    return `renewal::RenewalError` (new);
  - `recurrence::parse_rational_coeff` returns `recurrence::ParseRationalError`
    (new), and every other public `recurrence` function that returned
    `Result<_, String>` returns `recurrence::RecurrenceError` (new): the
    `Display` text is the former message, and the variant says whether the input
    is invalid (`InvalidInput`), unsupported, over a service ceiling
    (`LimitExceeded`), out of computation budget (`BudgetExhausted`, the search
    is incomplete rather than invalid), or a failed certificate replay
    (`VerificationFailed`); `Transform`, `Oeis` and `Rational` wrap
    `TransformError`, `OeisSequenceError` and `ParseRationalError`.  The
    functions are `RecurrenceJson::to_recurrence_parts`,
    `find_recurrence_representations{,_rational}`, `find_recurrence_spaces`,
    `verify_recurrence_spaces`, `RecurrenceCoefficientRescaling::{validate,
    rescale_coefficients, unscale_coefficients}`,
    `ConstantLagSource::{generate_rows, from_recurrence_json}`,
    `IndexedLagSource::generate_rows`, `RecurrenceSpaceLimits::check_ceiling`,
    `RecurrenceSpaceRepresentative::to_recurrence`,
    `RecurrenceLeaderCoverage::to_recurrence`,
    `GammaTransport::relation_system`,
    `TransformVectorRecurrenceJson::{text, to_recurrence}`,
    `TransformLimits::check_ceiling`, `transform_recurrence`,
    `verify_recurrence_transform`, `RepresentationSearchLimits::check_ceiling`,
    `search_representations`, `verify_representation_search`,
    `lossy_reconstruction`, `ImplicitRelation::{from_recurrence, from_json}`,
    `CoefficientCondition::from_json`, `assess_generating_description`,
    `ReconstructionCertificate::{verify_identity, replay_on_rows}`,
    `implicit_transform_recurrence{,_with_root_domain}`,
    `parse_implicit_transform_input`, `verify_recurrence_implicit`,
    `ImplicitSearchSupport::monomials`, `search_implicit_relations`,
    `verify_implicit_relation_search`, `RootInterval::at`,
    `RootDomainRequest::validate`, `root_implication`, `root_enclosures`,
    `jensen_fixed_cap_check` and `analyze_implicit_root_domain`.  The
    `String` fields of the serialized report structs, and so every JSON schema,
    are unchanged.
  Code that only prints the error, uses `?` or calls `.to_string()` is
  unaffected; code that compares the `String` or calls `str` methods on it must
  call `.to_string()` first.  The deprecated aliases that keep a `String` error
  are listed below.
- Spelled-out `tnn` names: public fields and enum variants cannot keep an alias,
  so `ChainPolynomialReport::tnn` is now `totally_nonnegative_verdict` (a
  `Result<(), NetworkError>`), `LeftProductionReport::tnn` is
  `totally_nonnegative`, and the variants `NetworkError::NotTnn` and
  `LeftProductionError::NotTnn` are `NotTotallyNonnegative`.  The methods and
  functions in the table below keep deprecated aliases.  CLI commands and flags
  (`tnn-network`), MCP tool names (`check_tnn_network`), JSON field names
  (`tnn`, `not_tnn`) and the `tnn_network` module path are unchanged.  In
  `polytool-mcp` the request, response and item types are now
  `TotallyNonnegativeNetworkRequest`, `TotallyNonnegativeNetworkResponse` and
  `TotallyNonnegativeFailureItem` (deprecated aliases keep the old names; the
  schema name `TnnFailureItem` and every JSON field are unchanged), and the fields named `tnn` are
  `totally_nonnegative` in Rust (still `tnn` in JSON).
- `decomposition::SymmetricDecompositionAnalysis` is generic over the
  coefficient type, `SymmetricDecompositionAnalysis<C = i64>`.  The unparameterized
  name is unchanged, and `SymmetricDecompositionAnalysisBigInt` is now a
  deprecated alias of `SymmetricDecompositionAnalysis<BigInt>`.
- `BigRational` is one definition, `pub use num_rational::BigRational` at the
  crate root; `basis`, `decomposition`, `recurrence` and `tnn_network` re-export
  it, so every existing path keeps working (they were four identical type
  aliases, and `polytool::BigRational` came from `tnn_network`).  Public
  signatures that used a private `Q` alias now name `BigRational`.
- `exact-real-roots` 0.2.0 and `polytool` 0.6.0 (lockstep with `polytool-mcp`
  and `polytool-web`); `polytool` requires `exact-real-roots` 0.2.

Renamed (deprecated aliases kept for one release):

Every old name still compiles and forwards to the new one with a
`#[deprecated]` warning; the aliases will be removed in the next release.

| Old name | New name |
|---|---|
| `exact_real_roots::satisfies_kurtz_condition_bigint` (also `polytool::` and `polytool::root_count::`) | `satisfies_kurtz_condition_bigint_coeffs` |
| `exact_real_roots::satisfies_newton_inequalities_bigint` (same re-exports) | `satisfies_newton_inequalities_bigint_coeffs` |
| `real_rootedness::isolate_real_roots_exact` | `isolate_real_roots_rational_coeffs` |
| `real_rootedness::real_roots` (also `polytool::real_roots`) | `approximate_real_roots` (returns interval midpoints, not exact roots) |
| `linalg::sparse_modular_linear_system_consistent` | `is_sparse_modular_linear_system_consistent` |
| `linalg::is_tnn` | `is_totally_nonnegative_neville` |
| `linalg::check_tnn_neville` (also `polytool::`) | `check_totally_nonnegative_neville` |
| `linalg::check_tnn_neville_bigint` (also `polytool::`) | `check_totally_nonnegative_neville_bigint` |
| `interlacing_matrix::check_lace_tnn_neville_i64` (also `polytool::`) | `check_lace_totally_nonnegative_neville_i64` |
| `interlacing_matrix::check_lace_tnn_neville_bigint` (also `polytool::`) | `check_lace_totally_nonnegative_neville_bigint` |
| `tnn_network::check_lower_unitriangular_tnn` (also `polytool::`) | `check_lower_unitriangular_totally_nonnegative` |
| `tnn_network::reconstruct_canonical_tnn_network` (also `polytool::`) | `reconstruct_canonical_totally_nonnegative_network` |
| `tnn_network::build_tnn_certificate_from_monic_polynomials` (also `polytool::`) | `build_totally_nonnegative_certificate_from_monic_polynomials` |
| `tnn_network::CanonicalTnnProof` (also `polytool::`) | `CanonicalTotallyNonnegativeProof` |
| `tnn_network::ChainPolynomialReport::is_tnn` | `ChainPolynomialReport::is_totally_nonnegative` |
| `left_production::check_production_tnn` | `check_production_totally_nonnegative` |
| `left_production::ProductionTnnFailure` | `ProductionTotallyNonnegativeFailure` |
| `linalg::is_totally_nonnegative` | `is_totally_nonnegative_by_minors` |
| `interlacing_matrix::check_lace_total_nonnegative_i64` | `check_lace_totally_nonnegative_by_minors_i64` |
| `interlacing_matrix::check_lace_sequence_total_nonnegative_i64` | `check_lace_sequence_totally_nonnegative_by_minors_i64` |
| `interlacing_matrix::is_lace_totally_nonnegative_i64` | `is_lace_totally_nonnegative_by_minors_i64` |
| `interlacing_matrix::is_lace_sequence_totally_nonnegative_i64` | `is_lace_sequence_totally_nonnegative_by_minors_i64` |
| `sequences::deco::polynomials` | `sequences::deco::deco_polynomials` |
| `sequences::deco::polynomials_bigint` | `sequences::deco::deco_polynomials_bigint` |
| `real_rootedness::ehrhart_to_hstar_with_denom` | `ehrhart_to_hstar_with_denominator` |
| `cyclic_sieving::cyclic_sieving_check_orbit_sizes` | `check_cyclic_sieving_orbit_sizes` |
| `decomposition::SymmetricDecompositionAnalysisBigInt` | `SymmetricDecompositionAnalysis<BigInt>` |
| `recurrence::SharedDenominatorVectorRecurrenceOptions` | `VectorRecurrenceDenominatorOptions` (unused alias, deprecated for removal) |
| `recurrence::VectorRecurrenceWithDenominator` | `SharedDenominatorVectorRecurrence` (unused alias, deprecated for removal) |
| `recurrence::VectorRecurrenceWithDenominatorFit` | `SharedDenominatorVectorRecurrenceFit` (unused alias, deprecated for removal) |

The deprecated `interlacing_matrix::check_lace_total_nonnegative_i64`,
`check_lace_sequence_total_nonnegative_i64` and
`cyclic_sieving::cyclic_sieving_check_orbit_sizes` keep their former `String`
error type (they are the only public functions that do); the new names return
the typed errors above.  `exact-real-roots` has no `Result<_, String>` API.

Fixes:

- Gamma coefficients beyond 64 bits: `polynomial_properties`,
  `check_polynomial_family` and the CLI `properties` / `family-check`
  reported no gamma coefficients (and `gamma_positive = false`) for a
  palindromic row whose gamma coefficients do not fit in `i64`; they now use
  the exact `BigInt` routine.
- CLI `family-check --recurrence` accepts coefficients of any size (it used the
  `i64` recurrence search).
- `hstar_inequality_report_bigint` no longer panics (debug) or wraps
  (release) when the supplied dimension is more than one below the degree of
  `h*`; the report has `codegree: None` and a failing `dimension bound` check.
- `tnn_network::verify_path_matrix_certificate` (and
  `CanonicalTotallyNonnegativeProof::verify`) returns `false` for a malformed
  network (a row `lambda[n]` without `n + 1` entries) instead of panicking.
- `interlacing_matrix::lace_matrix*` return `InterlacingMatrixError::SizeOverflow`
  for truncations whose size overflows, instead of panicking.  The size check
  now counts the outer vector of rows as well as the entries (so `block_cols = 0`
  no longer admits any `block_rows`, which panicked with a capacity overflow or
  aborted), rejects more than the new `MAX_LACE_DIMENSION` (`u32::MAX`) rows or
  columns, and the block loops no longer iterate over rows beyond `block_cols`.
  The MCP `check_polynomial_family` `lace` check bounds the rows and the block
  columns each by 16384 (the request fields' descriptions say so) and no longer
  overflows when `block_rows` is `usize::MAX` and `block_cols` is omitted.
- `cyclic_sieving_sequence_reports_bigint` no longer panics when
  `first_index + row` overflows `isize` and no longer allocates a report for
  a huge candidate order: it returns `CyclicSievingError::SequenceIndexOverflow`
  or `CyclicSievingError::OrderTooLarge` (new variants), and
  `cyclic_sieving_report_bigint` rejects an order above the new
  `CYCLIC_SIEVING_MAX_ORDER` (1000000) before allocating.  Candidate orders that
  overflow `isize` were skipped before; they are now `OrderTooLarge`.  The MCP
  `cyclic_sieving_sequence` tool and the CLI report these as errors; the MCP tool
  also requires every positive candidate order to be at most 10000, like
  `cyclic_sieving` (the `offsets` field's description says so).
- `hstar_inequality_report_bigint` no longer overflows `dimension + 1` or
  allocates a vector of the requested length for an absurd `dimension`: a
  dimension above the new `HSTAR_INEQUALITIES_MAX_DIMENSION` (4096), also when
  inferred from a longer input, gives a report with the single failing
  `dimension limit` check.  The MCP `hstar_inequalities` tool (`dimension`) and
  `polytool hstar-inequalities --dimension` reject larger values.
- `recurrence`: budget exhaustion while replaying a certificate
  (`verify_implicit_relation_search`, `verify_representation_search`) is
  `RecurrenceError::BudgetExhausted` instead of `VerificationFailed`;
  `check_kernel_certificate` failures of a tampered certificate (null vector,
  basis witness, rank minor inverse, kernel dimensions) and the replay checks of
  vector descriptions are `VerificationFailed` instead of `InvalidInput`;
  exceeding `max_space_shapes` is `LimitExceeded` and an unknown recurrence
  `sign` is `InvalidInput`.  `Display` texts are unchanged.
- `recurrence`: a `ParseRationalError` or `TransformError` raised while
  mapping or parsing is no longer flattened into an `InvalidInput` string:
  `RecurrenceError::with_context` wraps it in the new
  `RecurrenceError::Context { context, source }` variant, with the same
  `Display` text, and `Error::source` returns the original typed error.
- `tnn_network::CanonicalTotallyNonnegativeProof::evaluated_path_matrix` no
  longer panics on a hand-built network with a malformed `lambda`; it returns
  `None` (see Breaking changes).
- `recurrence_holds_rational` and `recurrence_holds_from_rational` accept
  `rec_len = usize::MAX`.
- Rustdoc: every public function that panics on caller-reachable invalid input
  (zero divisors, non-square matrices, a non-monic or inexact division, an
  `i64` convenience wrapper whose result does not fit) documents it under a
  `# Panics` section.

Maintenance:

- The recurrence implementation is grouped into private submodules by
  concern; all `polytool::recurrence` paths are unchanged.
- Lints are configured in each package's `[lints]` table (`missing_docs`,
  `rust_2018_idioms`, `forbid(unsafe_code)`).
- The README is a landing page; the CLI, library, MCP, web and algorithm
  references moved to `doc/`, and release notes to this file.
- The web page is split into `index.html`, `app.css` and `js/*.js`; serve the
  whole `web/` directory. All web tests live in `web/tests/`
  (`node --test web/tests/*.test.cjs`).
- The `generate_recurrence_benchmarks` example takes an output directory and
  no longer writes into `fixtures/recurrence-benchmarks/`, whose manifest it
  would have truncated to the synthetic fixtures.
- Examples: `bench_bezout_vs_sturm` no longer overflows and passes the
  interlacing arguments in the documented order.
- The MCP server uses `rmcp` 3.8 (was 0.16) and needs Rust 1.89. Tool names,
  descriptions and results are unchanged; input/output schemas encode
  optional fields as `"type": [..., "null"]` instead of `nullable: true` and
  drop the root `title`/`description`; arguments that fail to deserialize are
  reported in-band (`isError: true`) instead of as JSON-RPC errors.

Performance:

- The default exact real-rootedness test no longer switches large
  one-signed palindromic inputs to the Uspensky path and counts real roots
  directly (signs at infinity only); degree-79 Eulerian, type-B Eulerian and
  Touchard checks are about 10, 15 and 26 times faster (6.9 s to 0.72 s,
  14.5 s to 0.94 s, 4.9 s to 0.19 s).

Additions:

- Symmetric-companion lifts: `recurrence::find_symmetric_companion_lifts`
  takes a palindromic sequence `A_n` (scalar recurrence, derivatives allowed,
  with seeds and a symmetry schedule `d_n = d_0 + n - n_0`, supplied or
  inferred) and searches parameters `a_n, b_n, c_n` (degree at most 1 in `n`)
  such that the reconstructed companion `B_n = [A_(n+1) - U_n A_n]/(epsilon_n x)`
  makes `P_n = A_n + x B_n` satisfy `P_(n+1) = (a_n + b_n x) P_n + c_n x (1-x) P_n'`.
  Certificates (`SymmetricCompanionCertificate`, replayed by
  `verify_symmetric_companion_certificate`) prove the identity for every `n`
  by exact operator identities; the interlacing-preservation route
  (Brändén--Solus) is reported separately with its sign and seed hypotheses
  (`a_n >= 0`, `0 <= delta_n <= c_n`, nonnegative interlacing seed pair).
  Schema `polytool.recurrence-symmetric-companion.v1`.
- `sequences::literature::one_over_k_eulerian_polynomials_rational` and
  `one_over_k_eulerian_symmetric_parts_rational`: the `1/k`-Eulerian
  polynomials (`A_(n+1) = (1 + k n x) A_n + k x (1-x) A_n'`, exact rational
  `k`) and their symmetric decompositions.
- `approximate_real_roots_bigint_coeffs` (the `i64` `approximate_real_roots`
  now calls it) and `linalg::check_total_positivity_bigint` (bounded-minor
  total positivity with `BigInt` entries, sharing one kernel with the `i64`
  version); both are re-exported at the crate root.
- `factorization`: exact irreducible factorization over Q.
  `irreducible_factorization_bigint_coeffs` / `_rational_coeffs`
  (`_with_budget`) return a `RationalPolynomialFactorization` (constant and
  primitive integer factors with multiplicities; Yun squarefree decomposition,
  Berlekamp modulo a small prime, Hensel lifting, Zassenhaus recombination);
  `irreducible_factorization_monic_bivariate` (`_with_budget`) factors a
  polynomial monic in its first variable, with multiplicities (interpolated
  squarefree decomposition, a squarefree specialization and Hensel lifting).
  Errors are `PolynomialFactorizationError`.
- Veronese-section lifts: `recurrence::find_veronese_section_lifts` takes a
  constant-coefficient, derivative-free recurrence for `Q_n(y)` with its
  initial rows and returns parent recurrences `P_n(x)` with an exact identity
  `Q_n = V_(r,j)(P_n)` (the coefficients of `x^(rk+j)`), found by factoring
  `chi(L, x^r)` for the input characteristic polynomial `chi` and solving for
  the parent seeds exactly.  If every parent row is real-rooted with
  coefficients of one sign, every `Q_n` is zero or real-rooted; each lift
  reports exact checks of its first parent rows, and parents
  `P_(n+k) = c(x) P_n` (including every order-one parent) with real-rooted,
  one-sign `c` and seeds are proved outright.  Report types
  `VeroneseSectionLiftReport`, `VeroneseSectionLift`,
  `VeroneseSectionLiftStatus`, `VeroneseSectionLiftRejection` (JSON schema
  `polytool.recurrence-veronese-lifts.v1`); CLI `recurrence-veronese-lifts`;
  MCP tool `find_veronese_section_lifts`.
- Lazy row iterators.  `recurrence::RecurrenceRowIterator` (made with
  `Recurrence::generate_rows_rational_iter(initial_rows, first_index)`) yields the seed
  rows and then one exactly evaluated row per `next`, with
  `Item = Result<Vec<BigRational>, RecurrenceEvaluationError>`.  It evaluates
  nothing until a row is requested, keeps only the last `max_offset()` rows
  (plus unyielded seeds), can be paused and continued (`by_ref().take(k)`),
  yields a step error once and then ends (`FusedIterator`), and exposes
  `first_index`, `next_index`, `rows_yielded`, `steps_evaluated`,
  `history_limit` and `retained_row_count`.
  `Recurrence::generate_rows_rational` now collects a prefix of it (same rows,
  same errors).
- `oeis::OeisSequenceDefinition::{generate_rows_iter, generate_rows_from_iter}` and
  `oeis::OeisRowIterator` (`Item = Result<Vec<BigInt>, OeisSequenceError>`):
  lazy catalog rows over the shared kernel, with the authoritative prefix,
  displayed row indices, typed nonintegral errors and triangle padding of
  `generate_rows`, which now collects from it (`generate_rows_from` and
  `flattened_offset_for_row` no longer hold the skipped rows).  An error is
  reported when its row is reached, so with several failures the earliest row
  wins (the batch form used to report a step failure before an earlier
  nonintegral row).
- Lazy named families over the same kernel, yielding member `n` as the `n`-th
  item: `sequences::{eulerian_polynomials_iter, type_b_eulerian_polynomials_iter,
  chebyshev_polynomials_t_iter, chebyshev_polynomials_u_iter,
  hermite_polynomials_iter}` and
  `sequences::literature::{delannoy_polynomials_iter,
  delannoy_square_polynomials_iter}` (`sequences::PolynomialFamilyRowIterator`).
  The batch generators are unchanged.
- Example `lazy_rows`, and the integration tests in `tests/lazy_rows.rs`.

- `Display` and `std::error::Error` for `linalg::RationalKernelError`.
- New error types `ParsePolynomialError` and `CyclicSievingError` (see above).
- New module `pf_network`: height-shift-invariant planar networks
  (`HeightInvariantNetwork`, `PfStep`, `PfStepType`) with an exact planarity
  test, `build_pf_rows_certificate` / `PfRowsCertificate` (rows plus
  explicitly rebuilt Toeplitz windows), and `recognize_pf_network`, an exact
  column-by-column search over minimal planar supports with step-length,
  rise and budget bounds.  A certificate proves that every matched row is a
  Pólya frequency sequence (Brenti 1995): real-rooted with nonpositive zeros.
  CLI `pf-network`, MCP tool `recognize_pf_network`.
- New module `left_production` (Chen--Fu--Ruan, arXiv:2512.08369):
  `left_production_matrix`, `triangle_from_left_production_matrix`,
  `unified_network_toeplitz_window`, `check_production_tnn`,
  `LeftProductionCertificate::verify` and `check_left_production_matrix`.  A
  TNN left production matrix certifies PF rows and a TNN triangle and
  reversal for the supplied prefix; the report includes a negative-minor
  witness on failure and an explicit planar network on success.  CLI
  `left-production`, MCP tool `check_left_production_matrix`.

- New module `real_rooted_enumeration`:
  `enumerate_real_rooted_integer_polynomials_bigint` enumerates all
  real-rooted integer polynomials of fixed degree with at least three
  prescribed leading coefficients (Greaves--Syatriadi `AllRealRooted`, with
  certified exact endpoints instead of floating point), with node, result and
  coefficient-bit budgets; the kernel
  `real_rooted_antiderivative_integer_constants_bigint` (and
  `..._constant_multiples_bigint`) returns the integer shifts making an
  antiderivative real-rooted.  CLI command `enumerate-real-rooted`.

- `recurrence::find_recurrence_representations_rational`: bounded menus of
  candidate recurrence representations (lags, derivatives, polynomial
  coefficients).
- `recurrence::{find_recurrence_spaces, verify_recurrence_spaces,
  ConstantLagSource, IndexedLagSource}`: complete bounded recurrence spaces
  with exact certificates and replay, automatic right-hand-side degree
  bounds, collective leaders and polynomial-index propagation.
- CLI `--space` / `--verify-space`, matching MCP tools, and a guided
  full-search mode in the web app.
- Recurrence spaces: optional request field `coefficient_rescaling`
  (`recurrence::RecurrenceCoefficientRescaling`, kind `factorial_power`)
  certifies complete spaces for `Q[n] = T_m(P[n])`, `T_m(X^k) = X^k/(k!)^m`,
  from `r(U+m(A+J)+1)` residuals of a constant-index source; results carry a
  replayable `rescaling_bound` (`RecurrenceRescalingBound`).
- `recurrence::{RecurrenceTransformation, DegreeSchedule, transform_recurrence,
  verify_recurrence_transform}`: exact coordinate changes (`x^r`
  multiplication/removal, reversal, affine and Möbius maps, `f`-polynomial,
  gamma coordinates, symmetric decomposition) with inverses, image conditions
  and affine degree schedules on explicit index domains; exact transport of
  derivative recurrences (a vector system for the symmetric decomposition);
  replayable `polytool.recurrence-transform.v1` reports keeping the
  transported rule and an independent fit separately.  CLI
  `recurrence-transform` and MCP tools `transform_recurrence` /
  `verify_recurrence_transform`.  `decomposition::f_polynomial` now uses a
  homogeneous Horner scheme (same results and degree-bound semantics).
- Optional root-domain information for transformed results: request field
  `root_domain` (`recurrence::RootDomainRequest`) adds a `root_domain` report
  section (`RootDomainReport`) with theorem adapters (`root_implication`:
  monomial factors, real affine/Möbius maps, nonnegative-coefficient window,
  gamma quarter bound, Hermite--Kakeya--Obreschkoff symmetric pairs) recording
  direction, hypotheses, degree/sign/endpoint/zero conventions and failing
  converses; evidence-labelled claims; a Rolle recurrence match
  `Q_m = s(m) D^j(g Q_(m-1))` closing inductions (canonical Eulerian rows on
  `[-1, 0]`); and exact root observations.  Replayed by
  `verify_recurrence_transform`; reports without the field are unchanged.
- `real_rootedness::{isolate_real_roots_exact, RealInterval,
  RealRootIsolation, RealRootEnclosure}`: exact rational root enclosures with
  multiplicities, degree and valuation, and exact interval/ray tests;
  `recurrence::root_enclosures` with CLI `root-enclosures` and MCP tool
  `exact_root_enclosures`.  `recurrence::jensen_fixed_cap_check`: fixed-cap
  Jensen-polynomial test of diagonal multipliers (`fixed_cap` evidence).
- Implicit differential relations: `recurrence::{ImplicitRelation,
  CoefficientCondition, GeneratingDescription, assess_generating_description,
  DiagonalWeight, find_reconstruction_certificate, ReconstructionCertificate,
  implicit_transform_recurrence, verify_recurrence_implicit}` and
  `RecurrenceTransformation::{transport_implicit, reconstruction_operator}`.
  New registry maps `diagonal` (`T_n(x^k) = lambda(n,k) x^k`: factorial,
  `W_m`, inverse `W_m`, products of polynomial factors) and
  `derivative_antecedent` (`P_n = D^m(g Q_n)/m!`, exact back substitution).
  Diagonal transport uses one common clearing multiplier with certified
  exceptional-index conditions; uncertified sets stay visible. Generating
  descriptions record seeds, exceptional-index and boundary conditions, and a
  uniqueness argument separate from annihilation. Reconstruction certificates
  `L o T = sum H_nu S^nu M` record their differentiations and are replayed.
  `polytool.recurrence-implicit.v1` reports; CLI `recurrence-implicit`
  (`--verify`, `--assess`); MCP tools `transform_recurrence_implicit`,
  `verify_recurrence_implicit`, `assess_generating_description`.
- Bounded search over equivalent representations:
  `recurrence::{search_representations, verify_representation_search,
  RepresentationSearchRequest, RepresentationSearchLimits,
  RepresentationSearchReport, representation_label}`. A breadth-first queue
  applies a request-specified set of registry maps, deduplicates maps by
  certified equality of normal forms (never by matching rows), reconstructs
  the source exactly, transports certified relations through embedded
  replayable transform/implicit steps, certifies complete spaces after
  factorial maps, fits the new coordinates and feeds certified relations to
  later steps. Descriptions are labelled `finite_fit`,
  `certified_annihilating_relation` or `equivalent_generating_description`
  and ranked; limits produce explicit budget events. New evidence label
  `certified_complete_space` (`CERTIFIED_COMPLETE_SPACE_EVIDENCE`).
  `polytool.recurrence-representation-search.v1` reports; CLI
  `recurrence-search-representations` (`--verify`); MCP tools
  `search_recurrence_representations`,
  `verify_recurrence_representation_search`. Existing commands and defaults
  are unchanged.
- Direct search for implicit relations:
  `recurrence::{search_implicit_relations, verify_implicit_relation_search,
  ImplicitRelationSearchRequest, ImplicitSearchSupport,
  ImplicitRelationSearchReport}` fit `sum c(n,x) D^j Q[n+i] = 0` with
  derivatives on every row (roles `successor`, `implicit_successor`,
  `row_equation`) as one exact kernel per support with a rank certificate,
  reusing the complete-space columns. A constant-index monic pure-lag source
  under the identity, `factorial_power` or `derivative` map certifies the
  complete space (`certified_complete_space`); other fits are `finite_fit` or
  `certified_by_division` (new evidence `CERTIFIED_BY_DIVISION_EVIDENCE`,
  `M o T = sum H_nu S^nu L`). Coefficient conditions for exceptional indices
  are searched with a truncation theorem, and every representative carries a
  generating description and assessment.
  `polytool.recurrence-implicit-search.v1` reports; CLI
  `recurrence-implicit --search` / `--verify-search`; MCP tools
  `search_implicit_relations`, `verify_implicit_relation_search`.
- Lossy registry map `derivative` (`Q_n = D^m P_n`):
  `RecurrenceTransformation::{is_lossy, lost_data,
  inverse_row_with_lost_data}` (`inverse_row` rejects it), implicit transport
  `M = D^(m+r) L D^(-m)` with a recorded identity, pullbacks
  (`pull_back_relation_through_derivative`,
  `pull_back_condition_through_derivative`,
  `push_condition_through_derivative`, `relation_coefficient_condition`) and
  `lossy_reconstruction` / `LossyReconstruction`, which assesses the
  description of the source through `Q` with and without the lost data.
- Representation search: certified implicit-only relations of a parent are
  transported through diagonal maps, derivative antecedents and the
  `derivative` map as `implicit_relation` steps (`ImplicitRelationStep`) with
  carried coefficient conditions; lossy chains record `lost_data`
  (`LostDataRecord`) and are inverted with it. Replay re-transports every
  such step and recomputes the lost data. Reports without these maps are
  unchanged apart from one added convention line.
- `vec_poly::{exact_div_by_t, exact_div_by_t_minus_one, exact_div_by_t_minus,
  sum, prefix_sums, mul_linear, mul_linear_pow, has_nonnegative_coefficients,
  evaluate_homogeneous, evaluate_at_ratio}`: exact division by `t` and
  `t - a` (`None` on a nonzero remainder), sums and exclusive prefix sums,
  products with powers of `a + b t`, a coefficient sign test, and
  division-free evaluation at fractions.
- New module `hurwitz` with `hurwitz_matrix_bigint`,
  `hurwitz_leading_minors_bigint` and `is_hurwitz_stable_bigint` (exact
  Routh--Hurwitz criterion).
- `bezoutian_bigint_coeffs`: Bézoutian for arbitrary degree pairs; agrees
  with `bezout_matrix_bigint_coeffs` when `deg f = deg g + 1`.
- `sequences::shifted_legendre_polynomials_bigint`: shifted Legendre
  polynomials on `[0, 1]` with `P~_n(1) = 1`.
- The `check_weak_interlacing` documentation now spells out the idiom for
  treating a zero polynomial as vacuously interlacing.
- CLI `tnn-network`: TNN verdict for the coefficient triangle of monic rows
  (or, with `--matrix`, of a lower-triangular matrix), with the failing shape
  condition, negative entry or negative-minor witness; canonical network
  weights `lambda` and a reconstruct → evaluate → compare certificate;
  `--chain` adds the Brändén--Saud Maia Leite chain polynomials with exact
  real-rootedness, `[-1, 0]` and interlacing checks.  TNN of the triangle is
  reported separately from, and does not imply, real-rootedness of the rows.
- CLI `brenti-recognize`: exact planar strip-digraph search with
  `--first-row`, `--max-span`, `--weights`, `--max-search-nodes` and
  `--max-integer-nodes`; verdict `certified` (exit 0),
  `impossible_within_bounds` (exit 1) or `inconclusive` (exit 3).  A
  certificate proves TNN of the coefficient matrix, not real-rooted or PF rows.
- MCP tools `check_tnn_network` and `recognize_brenti_strip_digraph` with the
  same structured results (at most 128 rows for `check_tnn_network`; at most
  64 rows of degree at most 64 and bounded search budgets for
  `recognize_brenti_strip_digraph`).
- CLI `bench real-rootedness`: deterministic, verdict-checked regression suite
  (real-rootedness, strict interlacing, distinct-real-root counts and h\*
  conversion round trips) on Eulerian, Narayana, type B Eulerian, h\*-vector
  and positive-coefficient counterexample families of degree 20--80, with a
  `--quick` mode run by `cargo test`, TSV/JSON output and `--baseline`
  comparison.
- Internal: `linalg::sparse_modular` is split into the row type, a private
  `solver` submodule and a `tests` submodule; public paths are unchanged.
- `renewal::{VariablePredecessorRecurrence, PredecessorRule,
  PredecessorWeight, VariablePredecessorCertificate,
  recognize_fixed_delay_predecessor}`: symbolic variable-predecessor
  recurrences `Q_n = Q_(n-1) + w(n) t Q_(p(n))` with exact, replayable
  hypotheses certificates (fixed delay; polynomial weights proven nonnegative
  on all integers `n >= 1` by a Taylor-shift certificate, periodic weights by
  residue) and a canonical row generator; finite lists certify only a prefix.
- `renewal::{match_recurrence_theorems, RecurrenceTheoremReport}`: explicit
  conditional theorem matches (`weighted_renewal_window`,
  `fixed_delay_predecessor_window`) for found or supplied recurrences, schema
  `polytool.recurrence-theorems.v1`, kept separate from finite fit evidence.
- CLI `recurrence --theorem-matches` and `recurrence-theorems --recurrence
  <file|->`; default `recurrence` output is unchanged apart from a stderr note.
  MCP: additive `theorem_matches` field of `find_recurrence` (and the
  `check_polynomial_family` recurrence) and the tool `recurrence_theorems`.
- Web app: a **Transform** tab for coordinate transformations (with optional
  root location), implicit relations, exact root enclosures and conditional
  theorem matches. Every result carries a visually distinct evidence label
  (finite fit, checked on supplied rows, transported identity, certified for
  all indices, conditional theorem match). Reports can be copied, downloaded,
  imported and replayed exactly. New WASM exports `transform_recurrence`,
  `verify_recurrence_transform`, `transform_recurrence_implicit`,
  `verify_recurrence_implicit`, `root_enclosures`, `verify_root_enclosures`,
  `recurrence_theorems`, `verify_recurrence_theorems`, `verify_report`,
  `recurrence_latex` and `oeis_catalog_index`. The guided full search can
  rescale coefficients by `(k!)^m`.
- `recurrence::transport_gamma` and gamma transport in
  `RecurrenceTransformation::transport`: derivative rules are transported into
  gamma coordinates by clearing powers of `1 + x` and splitting the identity
  into its symmetric and antisymmetric parts; proportional parts give the
  scalar gamma rule (Eulerian, Narayana and type B Eulerian are covered by
  tests that regenerate the gamma rows), independent parts are reported with
  both relations. The gamma symmetry condition is also certified for rules
  with a left factor that is provably nonzero at every target.
- Root domain: the Rolle induction accepts factors `g(x)(c(m) + l(m) x)`
  (new optional field `index_linear_factor`) and index-dependent candidate
  intervals, checked symbolically (`I_(m-1) ⊆ I_m`, Taylor-shift sign
  certificates); index-dependent intervals map exactly through affine maps
  and through their index-domain hull otherwise; gamma inductions transfer
  to `P_n` only inside `(-inf, 1/4]`. New hypothesis status `violated`.
- Implicit reports: optional `root_domain` (`ImplicitRootDomainReport`,
  `analyze_implicit_root_domain`, `implicit_transform_recurrence_with_root_domain`,
  `parse_implicit_transform_input`) with multiplier-sequence and
  Rolle-antecedent adapters, Rolle inductions on the defining rule and on a
  certified fit, assumptions and a fixed-cap Jensen test; accepted by CLI
  `recurrence-implicit` and MCP `transform_recurrence_implicit`, replayed by
  `--verify` / `verify_recurrence_implicit`.
- Representation search: vector representations (after
  `symmetric_decomposition`) are transformed further. Scalar bijective maps
  act componentwise (`RecurrenceTransformation::{map_state, inverse_state,
  with_degree_offset, component_degree_offsets, is_mobius_class}`; a schedule
  acts on `b_n` at `d_n - 1`), map normal forms keep per-component suffixes,
  certified vector systems are transported exactly through Möbius-class maps
  (`transport_vector_componentwise`, steps `vector_transport` /
  `VectorTransportStep`, degree conditions certified from a certified
  decomposition), states get first-order vector fits (`vector_fit`), and
  vector systems carry `vector_assessment` (`VectorSystemAssessment`:
  uniqueness from a left-factor coefficient without integer roots), so a
  certified vector system can be an equivalent generating description and
  equal vector representations are detected. New optional report fields
  `RepresentationDescription::vector_assessment`,
  `RepresentationStep::vector_transport` and
  `RepresentationNode::implicit_spaces`; `TransformVectorRecurrenceJson::text`.
- Rank-two gamma transport: `GammaTransport::relation_system`
  (`GammaRelationSystem`) states the even part and its odd companion as
  implicit relations; `polytool.recurrence-transform.v1` reports carry them
  as `transported.system` (`GammaRelationSystemJson`) instead of a transport
  error, and the representation search lists them as
  `transported_gamma_system`, certified only with a certified symmetry
  condition (from the same or another step of the edge).
- Complete spaces after further maps: `search_implicit_relations` certifies
  complete spaces after Möbius-class maps (`multiply_power`, `divide_power`,
  `affine`, and `reversal` / `mobius` / `f_polynomial` with an unbounded
  schedule; transported-source theorem, `k (max(u+s)+1)` residuals) and
  diagonal maps with polynomial weights times `(k!)^(-m)` (`binomial`,
  `product` with positive exponents; weight theorem,
  `k (max(u+s)+w+1)` or `k (U+m(A+J)+w+1)` residuals), with coefficient
  conditions where their truncation theorem applies. The representation
  search embeds these as `implicit_spaces` (default menu or `space_shapes`).
  `ImplicitRelation::to_successor_rule` converts a successor relation back to
  a rule.

Real-rootedness contracts (found by a differential review; no default
real-rootedness or interlacing verdict changed):

- `satisfies_newton_inequalities_bigint_coeffs` now checks Newton's
  inequalities for coefficients of any sign, so it is a true necessary
  condition (it returned `false` for real-rooted polynomials with a negative
  coefficient, e.g. `(t - 1)(t - 2)`).
- `check_interlacing_sturm(0, g)` now returns `Some(true)` for real-rooted `g`,
  the documented convention of `check_interlacing`; both docs state it.
- `approximate_real_roots` docs: only a root `0` and the root of a linear
  input are exact.
- `tests/real_rootedness_differential.rs` compares every exact method against
  ground truth on about 84,000 cases (ignored by default; run with
  `--ignored`).

## 0.5.0

Breaking changes:

- `tnn_network::NetworkError::NotTnn` carries an optional explicit negative
  minor (`witness: Option<NegativeMinor>`), and `NetworkError` has the new
  variants `NotSquare` and `NonzeroAboveDiagonal`.
- The `brenti_sequence` documentation no longer claims that a valid strip
  digraph makes its row polynomials real-rooted.  Its transitions depend on the
  height (the digraph is invariant in `x`), so by Lindström--Gessel--Viennot it
  certifies that the coefficient matrix is totally nonnegative; heights
  `0 -> {0, 1}`, `1 -> {2}` give the row `1 + t + t^2`.

Additions:

- `tnn_network::{check_lower_unitriangular_tnn, NegativeMinor,
  lower_triangular_rows_from_square, chain_polynomials_bigint,
  chain_polynomials_rational, check_chain_polynomials_bigint}`: TNN diagnostics
  with explicit failing minors, and the chain polynomials of a
  lower-unitriangular matrix (real-rooted in `[-1, 0]` and interlacing when the
  matrix is TNN; Brändén--Saud Maia Leite, Theorem 3.7).
- `brenti_recognition::recognize_brenti_strip_digraph`: exact search for a
  planar strip digraph reproducing a finite prefix of rows (certificate,
  impossibility within the stated bounds, or inconclusive).
- `BrentiStripDigraph::{from_weight_matrix, weight_matrix}`.
- `exact-real-roots` 0.1.1: `isolate_real_roots_bigint_coeffs` with exact
  dyadic isolating intervals.

## 0.4.0

Breaking changes:

- `sequences::eulerian_polynomials(_bigint)` and
  `sequences::narayana_polynomials(_bigint)` now follow the module-wide
  `rows[n]` convention, include the row `n = 0` (`A_0 = N_0 = [1]`) and return
  `max_n + 1` rows; previously they returned `max_n` rows with
  `rows[i] = A_{i+1}` (resp. `N_{i+1}`).  Migrate with `eulerian_polynomials(n)[1..]`
  to get the old rows.  The CLI `sequence` command and the MCP `generate_sequence`
  inputs still list `A_1..A_n` and `N_1..N_n`.

Additions:

- `sequences::{derangement_descent_polynomials_bigint,
  fixed_point_descent_polynomials_bigint, binomial_eulerian_polynomials_bigint,
  refined_eulerian_polynomials_bigint, s_eulerian_polynomial_bigint}` and an
  "Indexing conventions" section in the `sequences` module documentation.
- `linalg::{characteristic_polynomial_bigint, closed_walk_counts_bigint}`.
- `numeric::{q_integer, q_factorial, q_binomial, q_catalan, q_narayana}`.
- `is_ultra_log_concave_infinite` / `_bigint_coeffs` (ULC(∞)).

## 0.3.0

Breaking changes:

- `vec_poly` is generic over `CoeffRing` instead of fixed to `i64`.
  Scalars are now passed by reference (`scale(&p, &c)`, `evaluate(&p, &x)`),
  every arithmetic result is trimmed (`neg` and `shift` included), and calls
  whose arguments are only integer literals need a type annotation such as
  `&[1i64, 1]`.
- `i64` (and the new `i128`) coefficient arithmetic in `vec_poly` and in
  `Polynomial` arithmetic panics on overflow instead of wrapping in release
  builds; the same holds for `Ratio<i64>`.
- `CoeffRing` gains borrowing arithmetic hooks (`add_ref`, `mul_ref`,
  `add_mul_assign`, ...) with default implementations; a `CoeffRing` impl for
  `i128` is added, which makes integer-literal type inference ambiguous where
  it previously defaulted to `i64`.

Additions:

- `vec_poly::{derivative, reverse, reciprocal, compose, pow, wronskian,
  shift_by, strip_trailing_zeros}` and `Polynomial::pow`.
- `numeric`: `binomial_bigint`, `factorial_bigint`, `falling_factorial_bigint`,
  `multinomial_bigint`, `binomial_u128`.
- `interlaces_either_order` / `interlaces_either_order_bigint_coeffs`
  (undirected, zero-tolerant weak interlacing) and
  `is_ultra_log_concave_of_degree` / `_bigint_coeffs` (ULC at a prescribed
  nominal degree).
- Expanded documentation of the directed and same-degree semantics of
  `check_weak_interlacing`.
