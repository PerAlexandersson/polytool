# Changelog

## Unreleased

- Certified rowwise products of two rule-defined polynomial families:
  `recurrence::{certify_polynomial_family_product, apply_polynomial_family_product}`
  support ordinary, coefficientwise, factorial-weighted and inverse-factorial
  products, with explicit parent index shifts. The new independent
  `PolynomialFamilyProductCertificate` schema records both parent proofs,
  exact annihilator closure, support and coefficient-anchor checks, and
  certified equality with the target recurrence. PF sign hypotheses are
  proved for all indices; ordinary products and Schur's factorial product
  also admit the appropriate whole-line parent proofs. Replay rejects
  altered certificates. Existing lift schemas and discovery routes are unchanged.

- Behaviour-preserving refactoring (no output changes). New public API:
  `vec_poly::clear_denominators` and `vec_poly::clear_denominators_with_scale`
  (integer coefficients of `lcm(denominators) * p`, with the scale), replacing
  about a dozen private copies; `oeis::parse_oeis_internal_format` and
  `oeis::parse_bfile` (with `OeisSourceRecord`, `OeisTermsSource`,
  `ParsedBFile` and the typed `OeisSourceError`), promoted from the CLI, which
  keeps its messages. Deprecated, as unused: `OeisRecurrencePresentation::displayed_recurrence`
  and `LeadingCoefficientCertificate::degree_steps_from`. New evidence gate
  `tests/discovery_fingerprints.rs` (fingerprints of the discovery report and
  the parent proof of every bundled catalog family; the full sweep is
  `#[ignore]`d) and shared integration-test helpers in `tests/common/`.
- Factor gauge with a bounded residual (RealRooted row routes `powerForm`,
  `closedForm` and `factorPower`). New theorem
  `recurrence::check_factor_gauge_theorem` (theorem `factor_gauge`,
  certificate `recurrence::FactorGaugeCertificate` with
  `FactorGaugeResidualShape` and `FactorGaugeDominantRoot`): for a rule
  `f(n) P_n = sum_j kappa_j(n) q(x)^j P_(n-j)` of order at least two with one
  fixed monic `q`, `P_n = q^(n-e) g S_n` with `S_n` in the fixed span of the
  gauged seed rows; `S_n` is shown real-rooted for every index by its degree
  (at most one), the discriminant (degree two), the formal cubic discriminant
  (degree three) or the palindromic-quartic criterion in `u = x + 1/x`, each
  a sign of a polynomial in `n` decided exactly when the coefficients of
  `S_n` are polynomials in `n` (proved by an identity; otherwise
  `not_applicable`), and nonzero by a coefficient polynomial, a dominant
  rational characteristic root or the leading-coefficient certificate.
  Proves OEIS A110200, A124376, A128541, A133109, A141611, A157169 and
  A157171 (catalog sweep 688 -> 695 proved, none lost).
  `RealRootednessTheoremReport` gains the optional field `factor_gauge`.
  `recurrence-prove` / MCP `discover_recurrence_proofs` record the attempt
  only when the rule has the theorem's form (other reports are unchanged,
  and saved tail or alternative attempt lists without it still replay); new
  parent route 8f (`VeroneseParentProofRoute::FactorGauge`) runs only for a
  parent left open and is recorded only when the rule has the form.
- Soundness fix: residue-class proofs with zero classes.
  `check_residue_class_subsequence_theorem` admits classes whose seed rows
  are all zero (infinitely many zero rows, as for OEIS A003982, A132918 and
  A286100), but its conclusion's `first_index` read as "every row from here
  on is nonzero and real-rooted", and `recurrence-prove` reported
  `real_rooted_from_index` for such families. The conclusion now records
  the zero classes in the new optional field
  `RealRootednessTheoremConclusion::zero_residue_classes`
  (`ZeroResidueClasses { period, residues }`, residues of the index `n`;
  helpers `has_zero_residue_classes` and
  `proves_row_nonzero_and_real_rooted`), and the discovery report omits
  `real_rooted_from_index` for a proof with zero classes and sets the new
  `real_rooted_outside_zero_residue_classes` (`from_index` and the classes;
  Python `ProofReport.real_rooted_outside_zero_residue_classes`, CLI text
  line); the status stays `proved`. Both fields are omitted otherwise, so
  every other report keeps its exact JSON; saved reports of the three
  families above no longer replay. A rule whose classes are all zero is not
  applicable. The parent provers keep route 8d for such parents, whose
  statement already allows zero rows.
- `recurrence-prove`, `recurrence-degree`, `recurrence-theorems` and
  `recurrence-generate` read the record from standard input (or
  `--input FILE`) when `--recurrence` is omitted, as
  `recurrence-veronese-lifts` already did; the `recurrence-prove` help lists
  the RealRooted row routes.
- Index-section lift (RealRooted row route `pathPower`, the inverse of
  `subseq`). New theorem `recurrence::check_index_section_lift_theorem`
  (theorem `index_section_lift`, certificate
  `recurrence::IndexSectionCertificate` with `IndexSectionTransform`,
  constant `MAX_INDEX_SECTION_PERIOD`): for a rule with multipliers depending
  on `x` only, the rows (possibly reversed, after an exceptional prefix row)
  are every `p`-th row (`p <= 5`) of a constant-coefficient parent found by
  factoring `chi_T(w^p)`, with `Res_w(psi(w), z - w^p) = chi_T(z)` and
  polynomial parent seeds; the identity holds for every index by an
  equal-order recurrence argument, and the parent is proved by the other
  theorems. Proves OEIS A109955, A109960, A116088, A116089, A127893, A193636
  and A206294; A114123, A184883 and A216182 are reported with an open
  parent. `RealRootednessTheoremReport` gains the optional field
  `index_section`. `recurrence-prove` / MCP `discover_recurrence_proofs`
  record the attempt only when a section identity is found (other reports
  are unchanged, and saved tail or alternative attempt lists without it
  still replay); new parent route 8e
  (`VeroneseParentProofRoute::IndexSectionLift`) runs only for a parent left
  open and is recorded only when a section identity is found.
- Proof discovery (issues #42, #45, #46 step 1). `recurrence-prove`
  (`discover_recurrence_proofs`, MCP `discover_recurrence_proofs`, Python
  `prove_real_rooted`) runs the iterated-preserver and hypergeometric-row
  criteria after the six earlier ones, in the order of the parent cascade
  (`DISCOVERY_THEOREMS`). When no criterion proves the record it tries the
  tail after rows that are zero or not real-rooted (`tail`,
  `real_rooted_from_index`) and alternative representations from the
  recurrence finder (`alternative_representations`), each used only after an
  all-index equality certificate (new `certify_recurrence_equality`,
  `RecurrenceEqualityCertificate`, schema `polytool.recurrence-equality.v1`).
  New flag `--max-representation-candidates` (MCP
  `max_representation_candidates`), new
  `discover_recurrence_proofs_with_options` and
  `AlternativeRepresentationOptions`. Inhomogeneous rules get the status
  `not_applicable` with `not_applicable_reason`. All report fields are
  additive. `RecurrenceProofDiscoveryReport::verify` now checks a saved
  report's claims instead of byte equality: the recorded attempts must be a
  subsequence of the recomputed ones, each equal or replaying on its own, and
  the summary must follow from them, so reports saved before a criterion was
  added keep replaying.
- `check_liu_wang_three_term_theorem` proves the degenerate case of
  proportional rows (`A_n = a_n g`, `B_n = b_n g^2`, rows `c_n g^k P_(N-2)`)
  by a positive scalar sequence or a hypergeometric closed form for `c_n`
  (which allows an isolated zero of `a_n`), tried only when nothing else
  proves the rule.
- Veronese-section lifts: supplied rows that do not follow the recurrence are
  dropped and the search runs on the tail (`tail_first_index`) instead of
  failing; an open parent is passed to the shared parent provers, with the
  new status `proved_parent_route` (`VeroneseSectionLiftStatus::ProvedParentRoute`).
- RealRooted row routes `variableDelay` and `subseq` (issue #46, first
  step).  New theorems `recurrence::check_variable_delay_theorem` (theorem
  `variable_delay`, certificate `recurrence::VariableDelayCertificate`):
  fixed-delay rules `f(n) P_n = A(n) P_(n-1) + B(n) x P_(n-L)`, `L >= 2`,
  with `A/f > 0` and `B/f >= 0` for every index, a constant first row and
  the seed rows as exact window steps, proving nonnegative coefficients,
  every root in `(-inf, 0)` and `P_A << P_B` for `A <= B <= A + L`
  (RealRooted `strictInterl_of_scaled_fixedDelay_recurrence`); and
  `recurrence::check_residue_class_subsequence_theorem` (theorem
  `residue_class_subsequences`, certificate
  `recurrence::ResidueClassCertificate`): rules whose lags are all multiples
  of `p >= 2` are proved one residue class at a time by the other theorems,
  each class report replayed, zero classes allowed.  `RealRootednessTheoremReport`
  gains the optional fields `variable_delay` and `residue_classes` (omitted
  from the JSON otherwise).  `recurrence-prove` / MCP
  `discover_recurrence_proofs` record the two attempts only when the rule
  has their form (other reports are unchanged); parent routes 8c and 8d
  (`VeroneseParentProofRoute::VariableDelay`, `ResidueClassSubsequences`)
  run only for a parent the earlier routes leave open and are recorded only
  where they apply.  Catalog sweep (`recurrence-prove`, 60 s per entry):
  519 -> 540 of 787 proved, no regressions; newly proved A003982, A099605,
  A102547, A109446, A109447, A118394, A118931, A118933, A132918, A137374,
  A171608, A201013, A211956, A223168, A286100, A305540, A317496, A317497,
  A318772, A318773, A321434.  **Breaking:** `RealRootednessTheorem` and
  `VeroneseParentProofRoute` gain two variants each; exhaustive matches need
  arms.  Discovery reports of rules with lags of gcd `>= 2` or of the
  fixed-delay form now carry one more attempt, so such reports written
  before this change still replay (see the discovery verification
  entry above).

- CLI usability (issue #41, phases 0-3). The command line is parsed by `clap`
  from one declarative specification (`src/bin/polytool/cli/spec.rs`,
  `commands/`) that also drives every help page, `polytool help --json`, shell
  completion and the deprecated-spelling tables.
  - Help: `polytool help <command> [<subcommand>]` at any depth (each `oeis`
    and `bench` subcommand has its own page), topic pages (`polytool help
    getting-started|input|coefficient-order|indexing|output-formats|
    json-schemas|exit-codes|budgets|recurrence-json|compatibility`), examples
    and exit statuses on every page, a grouped overview with a "Start here"
    block, and `polytool help --json` with flags, defaults, choices,
    positionals, subcommands, examples, exit codes and topics (additive keys
    in `polytool.commands.v1`). `polytool` without arguments prints the
    overview to stderr (exit status 1, as before).
  - Usage errors suggest the nearest flag, command or topic and end with
    ``Run `polytool help <command>` for usage.`` (stderr, exit status 2).
    Usage-error messages now come from clap and their wording changed.
  - `polytool completions bash|zsh|fish` completes subcommands, flags, flag
    values and help topics, not only command names.
  - Uniform options: `--flag=value` everywhere; `--format <...>` and `--json`
    on every command with structured output (`--format json` on JSON-only
    commands); `--input <file|->` on every command that reads input.
  - **Renamed flags (breaking-change notes).** Abbreviated and inconsistent
    flag names have spelled-out canonical names; the old spellings still work
    unchanged on stdout, but each use prints a deprecation warning on stderr
    (`POLYTOOL_DEPRECATION_WARNINGS=0` silences it). Scripts that require an
    empty stderr must switch to the new names or set the variable. The full
    table is in `doc/cli.md` and `polytool help compatibility`: the recurrence
    search flags (`--max-rec-len` -> `--max-recurrence-length`, `--max-var-deg`
    -> `--max-variable-degree`, `--max-idx-deg` -> `--max-index-degree`,
    `--max-diff-deg` -> `--max-derivative-order`, `--inhomogeneous` ->
    `--try-inhomogeneous`, `--denominator` -> `--try-denominator`,
    `--min-margin` -> `--min-equation-margin`, `--no-verify` ->
    `--no-held-out-rows`, ...); `--rows` -> `--count`; `--additional` ->
    `--additional-count`; `--first-row K`/`--start-row K` -> `--from-index K`
    (output selection) and `--first-row <0|1>` -> `--first-index <0|1>`
    (labelling the first input row); `--offset` -> `--first-term-index`;
    `--reversed` -> `--descending`; `--leading-desc` ->
    `--leading-coefficients-descending`; `family-check --recurrence` ->
    `--find-recurrence`; `--offsets` -> `--order-offsets`; the network
    budgets `--max-integer-nodes`/`--max-span`/`--max-length` ->
    `--max-integer-search-nodes`/`--max-target-span`/`--max-step-length`;
    the `bkw-scout` flags `--max-iter`/`--tol`/`--min-imag`/`--grid-re`/
    `--grid-im` -> `--max-iterations`/`--tolerance`/`--min-imaginary-part`/
    `--real-grid-points`/`--imaginary-grid-points`; `bench interlacing
    --max-n` -> `--max-index`. Hand-written per-command help functions and
    the old `help` module were removed.
  - Additive JSON keys: `first_index` next to `first_row` (`oeis generate
    --json`, OEIS catalog JSON) and next to `first_row_index`
    (`brenti-recognize`, `pf-network`, `left-production`).
  - MCP: request fields renamed to the same spelled-out words, with the old
    names accepted as serde aliases (`min_rec_len` -> `min_recurrence_length`
    and the rest of the search-option family, `first_row_index` ->
    `first_index`, `generate_oeis_rows` `rows`/`first_row` ->
    `count`/`from_index`, `generate_recurrence_rows` `rows`/`additional` ->
    `count`/`additional_count`, `maximum_start_delay` -> `max_start_delay`,
    `offsets` -> `order_offsets`); the published tool schemas show the new
    names. Responses add `first_index` next to `first_row` and
    `first_row_index`.
  - Tests: golden snapshots of every help page and topic
    (`tests/snapshots/cli/help`, re-bless with `POLYTOOL_BLESS_CLI_HELP=1`),
    frozen external invocations whose stdout and exit status must not change
    (`tests/cli_compatibility.rs`), and an old-versus-new spelling table,
    universal-option and flag-coverage tests (`tests/cli_help_and_aliases.rs`).
  - Dependencies: `clap` and `clap_complete` 4.5 (MSRV 1.74, compatible with
    the crate's 1.82).

- Python bindings (issue #40): the new crate `python/` (PyO3 0.30 and
  maturin, one `abi3` wheel for CPython 3.10 and newer, Rust 1.85) builds the
  package `polytool`, a thin adapter over the library.  Coefficients are exact
  Python `int` and `fractions.Fraction` values in ascending degree order;
  floats are refused.  The first surface: the `Polynomial` class (arithmetic,
  polynomial division over Q, evaluation, composition, shift, dilation,
  reversal, Bernstein, magic and arbitrary basis coordinates);
  `is_real_rooted`, `count_real_roots`, `isolate_real_roots` (exact
  enclosures with multiplicities), `interlaces` (directed, `None` when the
  degrees do not admit the test), `interlaces_either_order`, the coefficient
  predicates, `gamma_coefficients`, `discriminant`, `resultant` and
  `properties` (an item of `polytool properties --json`);
  `find_recurrence(rows, **options)` with the CLI search options under
  spelled-out names, returning a `Recurrence` (`generate`, lazy `iterate`,
  `latex`, `to_json`/`from_json`, `reindexed`, `match_theorems`) whose
  `to_json(indent=2)` equals `polytool recurrence --json`;
  `prove_real_rooted` (`ProofReport` with `status`, `route`, `certificate`,
  `attempts` and an exact `verify()` replay, the report of `recurrence-prove`)
  and `certify_leading_coefficient` (the certificate of `recurrence-degree`);
  and `polytool.oeis` (`ids`, `info`, `generate`, `recurrence`).  Reports
  convert to `dict` through their existing JSON schemas; notebooks get
  `_repr_latex_`.  A pytest suite checks big integers, fractions, round trips,
  the docstring and tutorial examples, and parity with the CLI JSON; the
  workflow `.github/workflows/python.yml` builds wheels on Linux, macOS and
  Windows and runs the suite.  Publishing to PyPI stays disabled until the
  owner enables it.  Guide: `doc/python.md`.

- New module `recurrence::leading_coefficient` (issue #37):
  `certify_leading_coefficient` decides a degree law `deg P_n = D(n)`,
  periodic by residue class (`RecurrenceDegreeLaw`: affine laws, half growth,
  Veronese sections, two-step products), and the sign of every leading
  coefficient of a recurrence family for all `n`, exactly.  The law is read
  off computed rows; the rows before the induction start are recorded
  exactly and every later index is covered by a route of the RealRooted
  tactic `rr_row_leadingCoeff_pos` (`cone`, `recurrence` with a rational
  growth-ratio invariant, `subseq`, `twoStep`, first-order `cancel`
  (optionally with the coefficient-ratio invariant `[x^(D-1)] P_n >= rho
  lead P_n`, which proves A101343), and
  `coefficients` with an affine degree envelope), each all-`n` condition
  decided by `decide_index_rational_function_sign` on the residue classes.
  The replayable `LeadingCoefficientCertificate`
  (`polytool.recurrence-leading-coefficient.v1`) names the route and its
  RealRooted declarations, reports a constant leading coefficient (monic
  rows) when the top-coefficient recurrence holds and fixes it (never on a
  class where the `cancel` route's top terms cancel), and a not-proved
  certificate lists its attempts (not a disproof).  `verify` replays the
  recorded route and refuses laws outside the searched range.  On the 860
  catalog entries it proves 784 laws (653 `cone`, 97 `recurrence`, 14
  `twoStep`, 12 `cancel`, 8 `subseq`), 390 of them with a constant leading
  coefficient.  New command `polytool recurrence-degree` (with `--verify`).
  The derivative-lag check takes its degree law from the module where its
  one-signed degree steps fail or an x-degree exceeds its bounds; the
  Liu--Wang check admits higher x-degrees when the module proves growth one;
  the Liu--Wang three-term check uses it in the identity coordinates where
  its degree-step state machine is unresolved (hypothesis `degree_law` in
  each); the hypergeometric-row check proves its first coefficient with it
  (the scalar recurrence of `U_0` as a family of constant rows) where the
  hypergeometric fit and the positive-coefficient test fail.  Newly proved
  theorem checks: derivative-lag for A046716, A054142, A082985, A084534,
  A104684, A121314, A126216, A172431, A209414, A360282 (a negative lag
  multiplier with a ratio invariant) and hypergeometric-row for A134346;
  A367270 (Liu--Wang three-term) is proved at the supplied seeds instead of
  one row later.  The parent cascade proves the same 585 catalog parents as
  before; A121314, A209414 and A360282 are now reached by the derivative-lag
  route (strict interlacing) before the Liu--Wang three-term route.  Every
  other proved report is unchanged; not-proved derivative-lag reports whose
  degree hypotheses failed now carry the shared hypothesis.  All earlier
  proved theorem reports and parent proofs replay
  (`check_recorded_parent_proof`).

- `polytool recurrence --max-work N`: a deterministic work budget for the
  adaptive search (`AdaptiveSearchBudget::max_work`, `with_max_work`).  Work
  units are computed from the dimensions of the exact and modular solves and
  checks, never from the clock, so the same input and `N` give the same output
  on any machine; the search stops with `budget_exhausted` (exit status 3).
  JSON reports `work_units` and `max_work` (and `search.work_units` for a found
  record); `AdaptiveSearchDiagnostics::work_units`.  Combines with
  `--max-candidates`; rejected with `--representations`.

- "Trivially real-rooted" proof classification:
  `recurrence::is_trivially_real_rooted(rule, first_index, seeds)` is true when
  the defining recurrence is `B(n) P_n = A(x,n) P_(n-k)` and the product theorem
  proves it.  New `trivially_real_rooted` field in
  `RealRootednessTheoremReport`, `RecurrenceProofDiscoveryReport`
  (`recurrence-prove`) and `RecurrenceTheoremReport` (`recurrence-theorems`,
  supplied records only; absent for fitted ones).  `real-rooted` sees only
  polynomials and has no such flag.

- `polytool recurrence --latex` (and the `latex` field of the JSON records, also
  with `--representations` and `--theorem-matches`) now factors every
  coefficient over `Q`: `(n + 2) P_{n} = (2n + 1)(1 + t) P_{n-1} - (n - 1)(1 - t)^{2}
  P_{n-2}`.  New `BivarPoly::irreducible_factorization` /
  `irreducible_factorization_with_budget` (univariate and monic fast paths,
  otherwise the shear `t = v + c n` onto the monic bivariate factorization;
  factors are primitive integer polynomials with a canonical sign; the product is
  re-expanded and compared exactly) returning `BivarFactorization` or a
  `BivarFactorizationError` (an exhausted budget is `Incomplete`, never a partial
  result); `Recurrence::to_latex_factored{,_with_budget}`.  A coefficient whose
  factorization does not finish stays expanded and is noted by a `%` comment
  line.  **Behaviour change:** `RecurrenceJson::to_latex` and `with_latex` now
  factor; `RecurrenceJson::to_latex_expanded` keeps the previous expanded text.
  `Recurrence::to_latex` is unchanged.

- New theorem `recurrence::check_hypergeometric_row_theorem` and parent route
  8b (`VeroneseParentProofRoute::HypergeometricRow`, tried after the iterated
  preserver, recorded only when the rule has a row relation of hypergeometric
  shape): rows whose coefficient ratio `[x^k] P_n / [x^(k-1)] P_n` factors
  into terms linear in `k` and affine in `n`.  The row relation
  `p(n, k) [x^k] P_n + q(n, k) [x^(k-1)] P_n = 0` is fitted exactly and
  proved for every `n` by an invariant row-annihilator system of the rule
  (the Hadamard-factor certificate); in the identity coordinates (valuation
  `v(n)` affine) or after the reversal at the degree, the ratio is matched
  with the image of `(1 + x)^M(n)` under a product of nonnegative multiplier
  sequences (complements `M'!/(M' - k)!` by Pólya--Schur, Theorem 7.1 of
  Brändén, arXiv:1410.6601; polynomials with roots `<= 0` by Laguerre's
  theorem; `1/(alpha)_k`, `alpha > 0`; a positive scale), every condition
  decided for all `n` by the sign checker, with the support of the rows from
  valuation and degree bounds and the first coefficient from its scalar
  recurrence.  It proves every row nonzero and real-rooted with every root
  in `(-inf, 0)` (`(-inf, 0]` with a valuation).  Of the 860 catalog
  entries, 28 parents are newly proved: A056939, A089503, A094305, A106246,
  A108426, A110608, A119308, A120986, A123160, A124848, A124860, A124932,
  A125092, A128540, A132812, A135065, A145596, A181543, A197208, A253283,
  A286785, A330798, A347056, A352880, A357613, A359363, A376827, A391045.
  Every other parent proof is unchanged (the route runs only when the earlier
  routes fail) and recorded proofs replay.  **Breaking:**
  `RealRootednessTheorem` and `VeroneseParentProofRoute` gain the variant
  `HypergeometricRow`; exhaustive matches need an arm.  The general proof
  discovery (`recurrence-prove`, MCP `discover_recurrence_proofs`) is
  unchanged.

- The iterated-preserver theorem (`recurrence::check_iterated_preserver_theorem`,
  parent route 8) recognizes Brändén's Weyl step `S_n = A + B D + C D^2`
  (P. Brändén, Trans. Amer. Math. Soc. 358 (2006), Corollary 3.5): the
  discriminant `B^2 - 4AC >= 0` on the line, `A` strictly interlacing or
  strictly alternating left of `B` with leading coefficients of the same
  sign, and `S_n(x^k)` of degree `deg A + k` with leading coefficients of one
  sign for every `k >= 0`, decided exactly for fixed operators and by the
  sign checker for index-dependent ones; coefficients of one sign keep every
  root in `(-inf, 0]`.  It is tried after the other menu operators, so the
  half-line Euler step still proves its rows (A036969).  A step free of `n`
  with checked alternating start rows also gives weak interlacing
  (Corollary 2.5).  The new `recurrence::check_weyl_step_preserver` decides a
  fixed operator, or refutes it as a global preserver at a rational point
  with a negative discriminant (Borcea--Brändén 2010, Theorem 1.3).  The
  theorem also accepts single-lag rules `P_n = T_n(P_(n-L))` (induction in
  each residue class, degree step decided class by class) and runs, after the
  identity, in the coordinates `recurrence::IteratedPreserverCoordinates`
  (power gauge `P_n / x^e`, translation `P_n(x + r)`, reversal at the nominal
  degree), each an exact Ore-algebra identity replayed by
  `IteratedPreserverCertificate::verify_operator_identity`.  A proof in the
  identity coordinates at a later induction start is preferred, so every
  report proved before is unchanged.  Of the 860 catalog entries, 31 parents
  are newly proved (most were settled before only by a `lean_proved`
  catalog flag): A008957, A048854, A079621, A086915, A088729, A102537,
  A104978, A105278, A130850, A133081, A137597, A143497, A151775, A176230,
  A191935, A241171, A243662, A269945, A286724, A290596, A290598, A292219,
  A304330, A304334, A305540, A321434, A331333, A345013, A370258, A371080,
  A397048.  The parent of A026386, proved by the product induction before, is
  now settled by route 8; its recorded proofs replay by their own route.

- `recurrence-prove` and MCP `discover_recurrence_proofs` try six general
  recurrence criteria with replayable exact certificates and explicit
  obligations; no named-sequence lookup or external attestation is used.
  `--verify` / `verify_recurrence_proofs` replays the complete report.
- General product induction supports index-dependent linear/quadratic and
  separable multipliers, plus bounded factorization of higher-degree
  multipliers with constant leading coefficient. It is also available to
  parent transfers when the required negative-root domain is proved.
  Rust proof certificates remain distinct from Lean replay.

- Recurrence alternatives have an aligned text comparison table, stable
  candidate numbers, verification counts and explicit checked singularities.
  The unsupported `--representations --theorem-matches` combination now fails
  before reading input and gives the existing per-candidate theorem command.
  Search results and JSON/LaTeX formats are unchanged.
- Polynomial input errors include physical line numbers; rational-list errors
  also identify the coefficient. Rational inputs reuse the library parser.
  Report row indices and accepted input syntax are unchanged.
- Invalid benchmark/PF-pencil arguments consistently exit with status 2.
  Usage hints go to stderr, keeping stdout clear on argument errors; explicit
  help/version requests reject extra arguments. Bash completion handles
  `help oeis` and stops suggesting commands in unrelated argument positions.
- `root-enclosures --verify` replays saved reports through a shared exact
  library verifier, also used by the browser. JSON-only enclosure, transform,
  implicit and representation-search commands accept explicit `--json`.
- Veronese lift searches accept theorem wrappers and representation menus
  through the same reader as recurrence generation and theorem matching.
- BKW scout validates formats and finite numerical options before searching.
  `--json --mathematica` emits one JSON document with a formula string, and
  non-finite numerical results fail before writing JSON.
- CLI help groups commands by task; `help --json` exports their names, aliases
  and summaries. OEIS maintenance subcommands now reach their detailed help.
  `completions bash|zsh|fish` emits command-name shell completion scripts.
- `sequence --format rows` supplies coefficient lists for pipelines;
  `recurrence-generate --format rows` names its existing default format.
  Recurrence generation and theorem matching accept theorem wrappers and
  representation menus, with explicit zero-based `--candidate` selection
  required for menus containing multiple entries.
- Invalid CLI options and input no longer silently succeed in the affected
  checking, recurrence, BKW and sieving commands. `resultant` requires exactly
  two polynomials. Resultants and discriminants support exact decimal-string
  `--json` output; unsupported options are rejected.
- Interlacing reports reuse successful strict tests for the weak result.
  Shared-root, same-degree, failed and incompatible cases retain their checks.
- A static documentation builder combines task recipes, an executable-derived
  CLI reference and searchable rustdoc. Its GitHub workflow builds artifacts
  and supports separately requested manual Pages publication.

- Square elimination delays scaling rows whose pivot-column entry is zero,
  using the same exact row-level invariant as the general solver. Sparse
  systems avoid unnecessary updates; pivot choice and exact results are unchanged.

- Square solves return immediately when exact elimination proves singularity,
  instead of repeating the solve with rational arithmetic. Non-exact division
  still requests the existing fallback.

- Exact recurrence fitting builds integer matrices directly when source rows
  and derivatives are integral. Fractional inputs retain the rational path;
  fitted coefficients may still be rational, and exact verification is unchanged.

- Rational square systems clear denominators once per equation, including the
  right-hand side, and use the existing fraction-free integer solver. Integral
  inputs retain their direct conversion path; exact results are unchanged.

- Exact linear solves stop as soon as a contradictory row proves inconsistency;
  consistent solves and complete kernel/certificate reductions retain their
  existing rank and budget contracts.
- Integer elimination combines delayed row scaling with its update and skips
  the products that cancel in the pivot column. Integer square solves also use
  integer back substitution, avoiding rational reductions in each dot product.
  Exact outputs and certificate budgets are unchanged.

This release is `polytool` 0.6.0, `exact-real-roots` 0.2.0 and the new
`ore-algebra` 0.1.0.  Existing CLI
commands and flags, MCP tool names and the web export names are unchanged; new
commands and tools are listed under Additions.  Some MCP responses now give
coefficients as decimal strings (see Breaking changes).

Breaking changes:

- `RecurrenceJson` has a new optional `latex` presentation field; Rust struct
  literals must add `latex: None`. Older JSON records remain readable. The
  recurrence generator ignores this field, and `reindexed` refreshes it when
  present. `Recurrence::to_latex` now uses subscripts and cleaner signs, while
  plain-text output keeps `P(n)`.

- `AdaptiveSearchResult`, `RecurrenceJsonSearch` and `RecurrenceRepresentation`
  have a new `verification_equations` field; Rust struct literals must add it.
  Older serialized records without this field remain readable (default 0).
- `LiftOperation::BasisMap::basis_family` is boxed to keep the enum compact
  as recurrence metadata grows. Its serialized JSON is unchanged.

- `recurrence::IteratedPreserverCertificate` has the new public fields `lag`
  (default 1) and `coordinates` (`IteratedPreserverCoordinates`, default the
  identity), both omitted from the JSON at their defaults, so existing
  certificates load and serialize unchanged; struct literals must add them.
  The certificate's `preserver` may now be a
  `PreserverOperator::WeylOperator` (Brändén's Weyl step).
- `recurrence::RealRootednessTheorem` and
  `recurrence::VeroneseParentProofRoute` have the new variant
  `LiuWangThreeTerm`; exhaustive matches must add it.
- `recurrence::RealRootednessTheorem` has the new variant
  `IteratedPreserver`, `recurrence::VeroneseParentProofRoute` the new variant
  `IteratedPreserver`, and `recurrence::RealRootednessTheoremReport` the new
  public field `iterated_preserver` (`Option<IteratedPreserverCertificate>`,
  omitted from the JSON when absent, so existing reports load and serialize
  unchanged); exhaustive matches and struct literals must add them.  The
  parent provers run the new theorem as route 8 (after Ma--Wang) and record
  its attempt only for rules of its form, so the recorded parent proofs of
  other rules are unchanged; the attestation and catalog-flag routes are now
  numbered 9 and 10 (route 7b, Liu--Wang three-term, runs before route 8).
  A first-order parent that only an attestation or a catalog flag settled
  before (without the root-domain route) is now proved; records made before
  keep replaying through the new `is_supported_by` checks (Additions).
- `recurrence::RootEndpoint` has the new public field `exact`
  (`Option<ExactQuadraticEndpoint>`); struct literals must add it (or use the
  new `RootEndpoint::rational`).  Serialized endpoints without the field still
  load, and an endpoint without an exact form serializes byte-identically.
  Request intervals (`RootDomainRequest` candidate intervals and assumptions,
  `RootEnclosureRequest::interval`) reject exact forms.
- The basis-map lift summaries (`root_transfer.input_root_region`, from
  `AlgebraicInterval::to_root_interval`) record an irrational endpoint's
  exact value next to its rational outer bound, and an open irrational end
  stays open (it was closed before; the outer bound lies strictly outside it,
  so both are correct).
- The real-rooted atlas engine (the `polytool::atlas` module and the
  `polytool atlas` command, added during this development cycle and never
  released) is removed again: polytool provides the transforms, lifts,
  certificates and the OEIS catalog, and the atlas is a separate application
  built on the public API (see Additions for the pieces made public).
- `recurrence::RecurrenceOptions` has the new public field `leader_x_power`
  (default 0) and `recurrence::AdaptiveSearchOptions` the new fields
  `try_laurent` (default false in the library, true in `polytool recurrence`
  and the MCP recurrence tools) and `max_laurent_x_power` (default 2); struct
  literals without `..Default::default()` must add them.  Serialized
  `AdaptiveSearchOptions` without the new fields still load; a found rule
  keeps `x^e` in its denominator, so `polytool.recurrence.v1` is unchanged.
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

- The Euler step of the preserver menu no longer attributes its theorem to
  Brändén's Corollary 3.5 (Trans. Amer. Math. Soc. 358 (2006)), whose
  hypothesis (i) fails for these steps; the theorem text now names the
  half-line RealRooted theorem `IsPFPolynomial.generalStep_of_natDegree_le`
  as what is applied.
- Exported adaptive-search rules (`polytool recurrence --json`,
  `recurrence-theorems`, `oeis add`, MCP `find_recurrence`) keep every row
  before the search window as an initial row.  A rule found with a window
  longer than its own lag (vanishing top coefficients) was exported with
  only `max_offset()` initial rows, so `recurrence-generate` produced wrong
  early rows.  New `AdaptiveSearchResult::export_initial_count`.

- `find_preserver_image_lifts` no longer aborts when one candidate's
  certificate fails its self-check (for example A061896 and A199478 stopped
  with "the pseudo-division identity does not hold"): such a candidate is
  rejected with a note in `PreserverImageReport::notes`, and an assembled
  identity or certificate that does not replay leaves the identity
  supported by rows only.  The cause was a missing factor: after a division
  remainder was removed by row annihilators, the index content `c(n)` of
  the annihilation multiplier entered the recorded multiplier `w(n)` but
  not the recorded quotients.  The Euler-step candidates of these entries
  are now certified identities (their side conditions fail, so nothing new
  is proved).

- `check_three_term_recurrence_theorem` and `check_liu_wang_theorem` normalize
  seed rows (trailing zero coefficients removed) before building the report,
  so reports from padded seeds, such as some catalog records, replay with
  `verify()`; previously such reports were rejected as differing from their
  recomputation (the verdicts were unaffected).

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

- The coordinate maps of the transform registry (`multiply_power`,
  `divide_power`, `affine`, `reversal`, `f_polynomial`, `mobius`) transport
  scalar rules through the forward rules of `ore_algebra::Substitution`
  instead of a separate chain-rule implementation; the transported rules are
  identical (checked against the former implementation on every rule of the
  bundled catalog).  The preserver-image coordinate maps (affine, reversal,
  Möbius) use the same substitutions for their row action and their
  dividend (`ore_algebra::pull_back_substitution`).

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

- The lift modules use the spelled-out `ore-algebra` API; the crate's hidden
  short method names (`mul`, `sub`, `eval_n`, `dy`, `x_valuation`, ...) are
  removed.  The lifts' own copies of the modular kernel, the operator
  pseudo-division and the gauge, affine and reversal conjugations are
  replaced by `exact_kernel_vectors_by_group`, `right_pseudo_division`
  (`right_pseudo_division_of_row` for systems), `pull_back` and
  `push_forward`.  Certificates change accordingly: pseudo-division
  multipliers carry the shifted parent leader of every level (possibly more
  exceptional bases, checked exactly), and a preserver dividend is divided
  by its common power of `x` (for a reversal `y^N` uses `N = b + h R`).

Performance:

- `ore-algebra`: integer operands skip the gcd reductions of the rational
  type (sums and products of integers are formed directly), polynomial
  translates use the binomial expansion, compositions translate each
  coefficient once per shift, operators are updated in place, the fits of
  `find_annihilator_system` and the ansatz equations are built without
  rational products, and the modular kernel caches its primes, reduces
  residues from the machine digits and eliminates only over the support of
  each pivot row.  Every result is unchanged (the new
  `tests/ore_algebra_workloads.rs` compares fingerprints of the complete
  outputs of a fixed workload set with recorded ones).  Release timings on
  that set (fastest of five runs): image recurrences of the preserver menu
  shapes 16 to 29 times faster on Eulerian, Narayana and Delannoy rows
  (Eulerian `1 + D` 1.21 s to 0.049 s, Narayana `theta (theta + 1)` 0.38 s
  to 0.020 s) and 4 to 7 times on binomial rows (`1/k!` 27 ms to 3.7 ms),
  annihilator systems 4 to 23 times, difference image annihilators
  (Touchard `Delta`) 5 to 8 times, right pseudo-division of large dividends
  26 to 33 times; the whole set 16.6 s to 2.7 s.
- `linalg`: the exact linear systems and kernels
  (`solve_linear_system`, `solve_linear_system_with_rank`,
  `rational_kernel_with_certificate`) eliminate on integer rows, with the
  fraction-free update `(v * row - t * pivot_row) / previous_pivot` (exact
  division), instead of `Ratio<BigInt>` entries that were reduced by a gcd
  after every operation.  Solving uses forward elimination with the smallest
  pivot row and back substitution on integers (the solution, rank, nullity
  and consistency do not depend on the pivot rows); the kernel keeps the
  Gauss-Jordan pivot order, so the basis, `independent_rows`, the inverse
  minor and the point at which a `RationalKernelLimits` budget is exceeded
  (operation count and coefficient bits, checked on the reduced rational
  values) are exactly those of the rational elimination.  The recurrence fit
  builds its equations without rational products for integer entries.
  Every result is unchanged (the new `tests/linalg_workloads.rs` compares
  fingerprints of the complete outputs of 24 workloads with recorded ones
  from the previous implementation, and unit tests compare against the
  rational reference).  Release timings on that set (best of two runs of
  three, 28.0 s to 1.6 s in all): `solve_linear_system_with_rank` on
  integer, sparse and rational matrices 49 to 67 times faster, kernels 11 to
  28 times, adaptive searches 4 to 38 times (Narayana with the default
  options 14.4 s to 0.38 s, shifted Legendre with `try_laurent` 4.7 s to
  0.28 s), the Hadamard lifts of the 321- and 123-avoiding families 2.7 and
  5.2 times (1.09 s to 0.41 s, 0.95 s to 0.18 s; the rest is the lift's own
  annihilator-system search).
- The default exact real-rootedness test no longer switches large
  one-signed palindromic inputs to the Uspensky path and counts real roots
  directly (signs at infinity only); degree-79 Eulerian, type-B Eulerian and
  Touchard checks are about 10, 15 and 26 times faster (6.9 s to 0.72 s,
  14.5 s to 0.94 s, 4.9 s to 0.19 s).

Additions:

- `polytool recurrence --latex` (alias `--format latex`) prints initial rows
  and the recurrence's starting index, one formula per line. Search JSON
  includes the same presentation in `latex`; coefficients remain expanded.
  `--representations --latex` prints the candidate menu as LaTeX definitions.
  `RecurrenceJson::{to_latex, with_latex}` expose the complete renderer in Rust.
- LaTeX output groups polynomial coefficients, omits unit multipliers, puts
  negative signs outside rational fractions, and handles zero right sides,
  forcing terms and derivatives with subscripted row indices.

- `polytool real-rooted --json` returns only the exact real-rootedness result,
  zero-based input index and degree for each polynomial; the zero polynomial
  has degree `null`. Plain-text output is unchanged.
- Recurrence search reports `verification_equations`, the dense coefficient
  slots checked on held-out rows before cancellation, not an independence count.
- `polytool recurrence --first-index N` uses the requested numbering in every
  output, with skipped input rows still counted. `recurrence-generate` accepts
  the same flag to relabel a record while preserving its generated rows.
  `RecurrenceJson::reindexed` exposes this exact index substitution in Rust.
  `RecurrenceRepresentations::reindexed` also updates complexity and the frontier.

- `recurrence::check_liu_wang_three_term_theorem`: Liu and Wang's Theorem 2.1
  and Corollary 2.4 (Adv. Appl. Math. 38 (2007) 542-560) for lag rules
  `P_n = A_n(x) P_(n-1) [+ V_n(x) P'_(n-1)] + B_n(x) P_(n-2)` with
  `x`-dependent coefficients, every hypothesis decided for all indices:
  degree steps `0` or `1` (half growth), `B_n <= 0` on the real line, on
  `(-inf, 0]` with an exact root-location induction, or on a root window,
  in the identity, divided-by-`x^m`, reversed and gauged coordinates; weak
  interlacing and the certified root region in a
  `RealRootednessTheoremReport` (`theorem: "liu_wang_three_term"`) with
  `verify()`.  It is route 7b of `prove_parent_real_rootedness`
  (`VeroneseParentProofRoute::LiuWangThreeTerm`, lift-summary route
  `liu_wang_three_term`), which settles the parent when the region lies in
  `(-inf, 0]`.
- Iterated real-rootedness preservers: `recurrence::check_iterated_preserver_theorem`
  (theorem `iterated_preserver`, with `recurrence::IteratedPreserverCertificate`
  in the report) proves a first-order rule `P_n = T_n(P_(n-1))` real-rooted
  from the last seed row on when its step is recognized exactly, as an
  identity of Ore operators, as `T_n = a_n(x) S_n` with `S_n` the identity or
  an operator of the preserver menu (`phi(n, D)`, `f(n, theta)`, the Euler
  step, `x^m D^j`, `D^j x^m`) and the multiplier `a_n` real-rooted (fixed
  factors exactly, index-dependent factors of degree `<= 2` by the sign
  checker), with every root in `(-inf, 0]`, `[0, inf)` or on the line, an
  exact degree step, and every earlier seed row checked exactly (zero seed
  rows allowed; the conclusion starts after the last one).  Rows that do not
  follow the rule, as in the 3-Pascal triangle `P_n = (1 + x) P_(n-1)` after
  `P_2 = 1 + 3x + x^2` (OEIS A028262), are covered this way, without a
  special product rule.  A linear multiplier also gives weak interlacing of
  consecutive rows.  An Euler step is applied in its degree-local form
  (`u(n) + v(n) k > 0` for `k <= deg P_(n-1)` instead of `v >= 0`).
  `verify()` replays the operator identity; the parent
  provers use the theorem as route 8 (region `(-inf, 0]`).  Of the 389
  first-order catalog rules, 202 are proved by it.
- Recorded certifications stay valid when the parent provers grow stronger:
  `VeroneseParentProof::is_supported_by` and
  `recurrence::check_recorded_parent_proof` accept a recorded parent proof that
  claims nothing beyond today's result (`not_proved`, or conditional on an
  attestation where the statement is now proved outright; a proved record
  needs the same route, statement and interlacing; attempts may differ), and
  reject one that claims more; `VeroneseSectionCertification::check_against_recomputation`
  and `LiftSummary::is_supported_by` (used by `BasisMapLiftCertification::verify`)
  apply the same rule to whole certifications.  A record proved by a route
  that the current cascade no longer reaches first (for example the paper
  theorem or the root domain, now preceded by the product induction) is
  replayed by that route alone and accepted only when the replay proves the
  parent with the recorded statement and interlacing (also for the
  coefficient family of a basis-map lift).  Recorded Veronese and basis-map
  certificates are kept as regression fixtures
  (`tests/recorded_certifications.rs`).
- Exact quadratic endpoints of root regions: `recurrence::ExactQuadraticEndpoint`
  (`a + b sqrt(d)` as rational strings, with `to_quadratic_number`,
  `from_quadratic_number` and the plain-ASCII `to_plain_text`, e.g.
  `-3-2*sqrt(2)`), the optional `RootEndpoint::exact` field (its rational
  `value` is a strict outer bound, `closed` refers to the exact end), checked
  by `RootInterval::check_exact_endpoints`, with `has_exact_endpoints`,
  `without_exact_endpoints`, `agrees_with_recomputed` (a record without exact
  forms agrees with an exact recomputation), `RootEndpoint::to_exact_text` and
  `format_root_interval_exact` (`[-3-2*sqrt(2), -3+2*sqrt(2)]`; the `Display`
  of `RootInterval` keeps the rational form).  The A008288 (Delannoy) basis-map
  lifts now record `[-3 - 2 sqrt(2), -3 + 2 sqrt(2)]` exactly, and
  `BasisMapLiftCertification::verify` accepts certificates written with or
  without the exact form and rejects a wrong one.  Exact affine images of
  index-dependent intervals keep the exact form; other interval transports
  work on the rational outer interval and drop it.
- OEIS-ready text from any recurrence: `oeis::OeisRecurrencePresentation`
  re-indexes a recurrence (or a `polytool.recurrence.v1` record) to the row
  numbers of an OEIS triangle, read as `oeis::OeisTriangleReading` says (first
  row number, ascending or descending powers, zero trimming, regular-triangle
  padding), and formats the `%e` "Triangle begins:" block
  (`oeis::format_triangle_example`), the `%F` recurrence sentence in plain
  ASCII, and `%t` Mathematica and `%o` PARI/GP programs that list the same
  rows (`OeisText`, schema `polytool.oeis-text.v1`; errors `OeisTextError`).
  Golden tests cover four recurrence shapes; an opt-in test
  (`POLYTOOL_TEST_WOLFRAMSCRIPT=1`) evaluates the Mathematica programs with
  `wolframscript`.  The PARI/GP programs are covered by golden tests only.
  `polytool recurrence-generate` gains `--format
  oeis-example|oeis-formula|mathematica|pari|oeis-json`, and `--first-row` and
  `--regular-triangle` (for these formats and `--format bfile`).  The web
  interface exports `oeis_text_for_recurrence(request)` and
  `generate_oeis_rows(request)` (OEIS-aligned rows and their b-file text,
  `polytool.oeis-rows.v1`, at most 1000 rows or 200000 terms).

- B-files from any recurrence: `oeis::flatten_rows_to_bfile` (with
  `oeis::BFileOptions`, `BFileRowOrientation`, `BFileZeroPolicy`,
  `BFileStoppingRule`, `BFileText`, `BFileError` and `oeis::bfile_row_terms`)
  flattens integer coefficient rows into b-file text with an offset, an optional
  reversed row orientation, a zero-trimming policy and a term cap that never
  splits a row (at most the cap, or the whole row nearest to it).
  `polytool recurrence-generate` gains `--format text|json|bfile` with
  `--offset`, `--max-terms` (rows are generated until the cap is reached when
  neither `--rows` nor `--additional` is given), `--nearest-row`, `--reversed`
  and `--trim-zeros`.  `polytool oeis generate --format bfile` now uses the same
  library function; its output is unchanged.

- Public building blocks of the lifts, for applications that chain lifts
  over their own parent records (previously crate-internal):
  `recurrence::prepare_lift_rule` with `recurrence::PreparedLiftRule`
  (a validated rule with its identifying seeds),
  `recurrence::prove_parent_real_rootedness` with
  `recurrence::ParentRecurrenceRecord` (the parent provers),
  `recurrence::find_bounded_root_region_report` (a replayable theorem report
  for a bounded root interval) and
  `PreserverParentLibrary::add_parent_record` (one named parent record).

- Computed gamma and Veronese image rules:
  `recurrence::GammaLiftParent::gamma_polynomial_recurrence` gives a runnable
  rule with seeds for the gamma polynomials of a palindromic parent, from the
  gamma transport of its rule (`transport_gamma`), and
  `recurrence::VeroneseSectionCandidate::section_recurrence` (options
  `VeroneseSectionRecurrenceOptions`) one for `lambda(n) V_(r,j)(P_(n+s))`,
  from the transported section system and
  `ore_algebra::find_component_annihilator`, or
  `RecurrenceError::BudgetExhausted`.  Both are normalized like the image
  rules: rows at bases where the relation is not proved (below the parent's
  identifying rows, roots of the division multiplier or of the content) are
  checked exactly and kept as seeds where the relation fails there.
  Callers that previously fitted rules to the lifted rows can use them
  instead; on known gamma and Veronese lifts the computed rules give the
  same rows, and most Veronese sections that timed out with fitted rules now
  finish in well under a second.

- OEIS catalog: 43 more entries (860 in all), each a supplied recurrence
  validated term by term against the OEIS data. Among them are 10 entries
  whose recurrences the adaptive search now finds only through its
  left-hand-side variation, with a left-hand factor in `n` alone (`n`,
  `n ± 2n^2`) and zero constant term: A085880, A088617, A101919, A102537,
  A120986, A243662, A280580, A281000, A390725, A391045.

- OEIS catalog real-rootedness flags: a record may carry
  `"real_rootedness": {"status": ..., "as_of": "YYYY-MM-DD"}` (schema
  `oeis::OeisCatalogRealRootedness`, `OeisRealRootednessStatus::{LeanProved,
  KnownNotRealRooted}`; absent means unknown), read with
  `OeisSequenceDefinition::real_rootedness()` and listed by
  `oeis::entries_with_real_rootedness_flag()`; `oeis info` and `oeis list
  --json` report it.  `lean_proved` records that every row is zero or
  real-rooted by a machine-checked Lean proof made outside polytool, for a
  definition whose rows were shown equal to the catalog rows for every index
  (same recurrence normal form after the aligning index shift, nonvanishing
  leaders, exact agreement of the initial rows).  The bundled catalog flags 753
  entries `lean_proved` and A122431 `known_not_real_rooted` (its row
  `3 + 3x + x^2`).  The flag is never reported as a polytool proof.
- OEIS catalog evidence records: a record may carry, after `real_rootedness`,
  `"evidence": {"strongest": ..., "rows_checked": N, "first_failure": {...},
  "as_of": "YYYY-MM-DD"}` (schema `oeis::OeisCatalogEvidence`,
  `OeisEvidenceProperty::{RealRooted, UltraLogConcave, LogConcave, Unimodal,
  NoProperty}`, `OeisEvidenceFirstFailures`): the strongest property that holds
  on every one of the first `rows_checked` rows (exact, none skipped) and the
  first failing row of each stronger property.  Read with
  `OeisSequenceDefinition::evidence()` (`oeis::OeisEvidence`), computed with
  `OeisSequenceDefinition::compute_evidence`, `oeis::OeisEvidenceAccumulator`
  and `oeis::row_has_evidence_property`; `oeis info` and `oeis list --json`
  report it.  New CLI `polytool oeis evidence [A-number...|--all] --rows N
  [--write] [--per-entry-timeout-seconds S]` computes and records it.  It is
  evidence about finitely many rows, not a proof, and separate from the
  `real_rootedness` flag.  `OeisCatalogRecord` has the new field `evidence` and
  `OeisCatalogError` the variants `InvalidEvidenceDate` and
  `InconsistentEvidence`.  The bundled data file holds no evidence records yet.
- Parent provers: a Veronese, Hadamard-factor, row-factor or gamma lift parent
  whose rule, first index and initial rows equal a `lean_proved` catalog entry
  exactly gets the verdict `conditional_on_attestation` with the label
  `oeis_catalog:<A-number>:lean_proved` (the caller-attestation mechanism),
  provided polytool proves that every parent row has coefficients of one sign
  (nonnegative seeds, positive `x`-free leader, and nonnegative combined
  coefficient weights under an inductive degree bound, shown exactly in a
  falling-factorial basis).  Used directly as parents, 138 flagged entries
  that no polytool route proves become conditional.
- `polytool oeis add --recurrence-json FILE [--prefix-rows P]` records a
  supplied `polytool.recurrence.v1` recurrence without a search, after
  checking every OEIS term exactly through the catalog's row generator; new
  provenance source `supplied_recurrence` (status `validated`, with
  `checked_terms`, `terms_source` and the last compared row).  The catalog
  gains 43 such families (817 entries).  `OeisCatalogError` gains
  `InvalidSuppliedProvenance` and `InvalidRealRootednessDate`, and
  `OeisProvenanceSource` gains `SuppliedRecurrence`.
- Hadamard-factor lifts: `recurrence::find_hadamard_factor_lifts` writes the
  input rows as `Q_n = lambda(n) sum_k u(n, k) p(n, k) x^k` with `u` from a
  menu (`binom(n + a, k + b)`, `k!`, `1/k!`, `(n - k)!`, `1/(n - k)!`,
  `(alpha)_k`, `1/(alpha)_k`, optionally products of two) and a parent with a
  smaller rule, certifies the identity for every index by conjugating the
  input relation and dividing it by the parent relation
  (`HadamardFactorCertificate`, replayable), transfers real-rootedness by the
  Pólya--Schur characterization of multiplier sequences and runs the parent
  provers; report `polytool.recurrence-hadamard-factor-lifts.v1` with a
  `lift_summary` per lift (new `LiftOperation::HadamardFactor`).  The Narayana
  and type-B Narayana triangles are proved end to end.
- Gamma lifts: `recurrence::find_gamma_lifts` and `certify_gamma_lift`
  recognise `Q_n = lambda(n) Gamma_(n+s)(c z)` with `Gamma_m` the gamma
  polynomial of a palindromic parent (menu `GammaLiftParent::standard_menu`:
  Eulerian, type-B Eulerian, Narayana; also
  `GammaLiftParent::one_over_k_eulerian_symmetric_part`), certify the parent
  symmetry for every row by the reversal induction and the identity by
  pseudo-division of the input relation by the transported gamma rule
  (`GammaLiftCertificate`, replayable), and transfer real-rootedness with
  roots `<= 0` by the gamma criterion; reports
  `polytool.recurrence-gamma-lift.v1` and `polytool.recurrence-gamma-lifts.v1`
  (new `LiftOperation::GammaPolynomial`).  The peak (and valley) polynomials
  A008303 are proved from the Eulerian polynomials.
- Adaptive recurrence search: left-hand side variation.  After the two
  historical passes (whose results and candidate counts are unchanged), a
  third complexity-ordered stage tries left-hand factors `x^e f(n, t)`
  (`RecurrenceOptions::leader_x_power`, `e <= max_laurent_x_power`) and
  normalizations of `f` on a coefficient other than its constant term, so
  factors such as `n` (shifted Legendre) are reached without
  `try_denominator`.  `AdaptiveSearchOptions::with_laurent` toggles it; the
  CLI has `--no-laurent` and `--max-laurent-x-power`, the MCP recurrence
  tools `try_laurent` and `max_laurent_x_power`.
- Left-hand factors `x^e f(n)` in the theorem checks and lifts: a rule whose
  left-hand factor is `x^e f(n)` (`e >= 1`, `f` free of `x`, as the adaptive
  search stores it) and whose `x^e` divides every right-hand coefficient is
  normalized exactly to `f(n) P_n = sum (c / x^e) D^d P_(n-r)` (the same rows)
  before the three-term, Liu--Wang, derivative-lag and Ma--Wang checks, the
  parent provers of the lifts, the shared lift rule preparation (Hadamard,
  row-factor, basis-map, gamma and preserver-image lifts), the Veronese-section
  preparation and the Rolle root-domain matcher.  The new record
  `recurrence::LeftHandFactorNormalization` (`x_power`, `original_rule`,
  `normalized_rule`, `statement`; `from_recurrence`, `normalized_recurrence`,
  `verify`) appears in the new optional fields
  `RealRootednessTheoremReport::left_hand_factor_normalization` and
  `VeroneseParentProof::left_hand_factor_normalization` (omitted from the JSON
  when no normalization was applied, so existing reports are unchanged);
  `verify()` re-derives it.  When `x^e` does not divide the right-hand side
  the division happens only on the rows, and the checks report
  `not_applicable` (the lifts `Unsupported`) with that reason.
- Preserver-image lifts: `recurrence::find_preserver_image_lifts` and
  `certify_preserver_image_lift` recognise `Q_n = lambda(n) T_n(P_(n+s))` for
  a proved parent `P` (`PreserverParentLibrary`: named families, optionally
  OEIS catalog entries) and an operator `T_n` of a menu (`PreserverOperator`:
  `f(n, theta)`, `phi(n, D)`, `D^j o x^m`, `x^m D^j`, affine, reversal,
  Möbius, hypergeometric multipliers, general Weyl operators without a
  theorem); the identity is certified by pseudo-division of the gauged input
  relation composed with (or conjugated through) the operator by the shifted
  parent relation, with row annihilators for a remainder, or by the Hadamard
  conjugation (`PreserverImageCertificate`, replayed from the recorded rules
  and seeds); side conditions are decided exactly for every `n`; reports
  `polytool.recurrence-preserver-image.v1` and
  `polytool.recurrence-preserver-image-lifts.v1` (new
  `LiftOperation::PreserverImage`), with the Lean route and RealRooted
  declarations of each theorem and an interlacing claim
  (`PreserverInterlacingClaim`, Hermite--Kakeya--Obreschkoff transfer of the
  parent's proved interlacing for operators independent of `n`).
- Ma--Wang interval theorem: `recurrence::check_ma_wang_interval_theorem`
  (theorem `ma_wang_interval`, report
  `polytool.recurrence-real-rootedness-theorem.v1`) for first-order rules
  `f(n) P_n = c_0(n, x) P_(n-1) + gamma(n) (x - L)(x - U) P'_(n-1)`,
  `deg_x c_0 <= 1`, `L < U` rational: every row real-rooted with every root in
  `[L, U]` and weak interlacing of consecutive rows, by the exact two-endpoint
  criterion of Ma and Wang (EJC 15 (2008) R17) in the closed form of RealRooted
  (`affineEndpointStep_closed`), with the endpoint multiplicities tracked
  exactly through the induction (sign checker between endpoint hits, and the
  lattice argument for affine thresholds with an infinite hit sequence); a
  failing threshold names the row that leaves the interval.  It is a new parent
  route (`VeroneseParentProofRoute::MaWangInterval`, lift-summary route
  `ma_wang_interval`) when `U <= 0`, run after the derivative-lag theorem.
- Preserver-image lifts with interval parents: `PreserverParentRegion::Interval`
  and `PreserverParent::interval_proof` (the Ma--Wang interval theorem, or the
  three-term theorem with a root region), and the certification field
  `PreserverImageCertification::interval_parent_proof` (replayed by
  `verify()`).  Affine, reversal and Möbius maps and the Rolle operators give
  the exact image of the parent interval (lift-summary routes
  `ma_wang_interval` and `three_term_interval`); multipliers and `phi(D)` use
  their ray theorems only.  The named Chebyshev and shifted Legendre parents
  are now interval parents.
- Preserver-image operator `PreserverOperator::EulerStep` and shape
  `PreserverOperatorShape::EulerStep`: the Euler step
  `kappa (theta + a)(theta + b) + x (u + v theta)` with the side conditions of
  RealRooted `IsPFPolynomial.generalStep_of_natDegree_le` (`kappa, a, b > 0`,
  `u > 0`, `v >= 0`, negative comparison defect), decided exactly for every
  `n`; region `(-inf, 0]`, no interlacing claim.
- New standalone crate `ore-algebra` (`crates/ore-algebra`, dependencies
  `num-*` only): the exact Ore algebra of shifts and derivatives
  (`BivariatePolynomial`, `OreOperator`, the Euler form), moved out of
  `recurrence::implicit` with unchanged behaviour, plus right
  pseudo-division, maps between sequences with the transport of annihilators
  (`RowMapStep`, `pull_back`, `push_forward`: Weyl operators, gauges, affine
  maps, reversals, coefficient multipliers), jointly invariant systems of
  additional annihilators of recurrence-defined rows
  (`find_annihilator_system`), and runnable image annihilators found by
  push-forward or as common left multiples from an increasing ansatz with an
  exact modular kernel (`find_image_annihilator`), each with an identity
  certificate that replays on its own.
- `ore-algebra`: the variable shift `E f(x) = f(x + 1)`
  (`DifferentialDifferenceOperator`, pseudo-division by a recurrence and in
  `E`), the falling-factorial basis change with its conjugation
  `x -> x (1 + D)`, `E -> 1 + D` and Euler form, maps with `E` and certified
  image annihilators (`DifferenceRowMapStep`, `pull_back_difference`,
  `find_difference_image_annihilator`), Laurent coefficients in `x`
  (`LaurentOreOperator`: clearing, reversal, power gauge and dilation as
  conjugations), and substitutions as first-class operations
  (`Substitution`: affine, reversal, weighted Möbius, `f`-polynomial, power
  gauge, with row maps, inverses, pullbacks and forward rules).
  certificate that replays on its own; right pseudo-division of a row of
  operators by a square system with leader `g(n) S^k I`
  (`right_pseudo_division_of_row`, `system_index_only_leader`) and runnable
  scalar annihilators of one component of such a system
  (`find_component_annihilator`, `ComponentAnnihilatorCertificate`, the
  identity `w(n) (M e_c) = sum_b H_b o Phi_b`).
- Image recurrences: `recurrence::image_recurrence_of_operator` computes a
  runnable rule (top shift free of derivatives, leader decided nonzero after
  the seeds) for `Q_n = lambda(n) T_n(P_(n+s))` from a runnable parent rule,
  for a Weyl operator, an affine, reversal or Möbius map, or a coefficient
  multiplier `x^k -> u(k) x^k` (`ImageOperatorSpecification`, `ImageRowMap`,
  `ImageGauge`), with seeds computed as `T_n(P_n)` and an
  `ImageRecurrenceCertificate` (the identity `w E = H o L + Z o A`, the
  invariance of auxiliary annihilators, exceptional bases checked on rows,
  replay) that `verify()` recomputes; schema
  `polytool.image_recurrence.v1`.  When no runnable rule is found within the
  budget the result is `RecurrenceError::BudgetExhausted`.  For example
  `A_n + A_n'` (Eulerian) gets an order-2 rule, `theta (theta + 1) N_n`
  (Narayana) and the `1/k!` image of the binomial rows (Laguerre) order-1
  rules.  Preserver-image candidates convert to image specifications
  (`ImageRowMap::try_from(&PreserverOperator)`, every operator except the
  multipliers that depend on `n`;
  `PreserverImageCandidate::image_operator_specification` and
  `image_recurrence`), and with `PreserverImageOptions::image_rule` a
  certified preserver image carries its computed rule
  (`PreserverImageCertification::image_rule`, replayed by `verify()`).
- Basis maps and basis-map lifts: `recurrence::BasisMap` (the map
  `T_P(u^k) = P_(i0 + k)` of a recurrence-defined family as a serializable
  value, with `apply`, the exact triangular inverse `expand`, `tail` and the
  intertwining relation with multiplication by `u`,
  `BasisMapIntertwining`); `check_basis_map_preserves_pf` decides exactly the
  hypotheses of the basis-map theorems (geometric, eventually geometric,
  constant- and affine-current order-two, ray and interval differential,
  Hermite-type and index-dependent scalar families) and reports source and
  target intervals with exact quadratic endpoints (`QuadraticNumber`,
  `AlgebraicInterval`; for example `[-1 - sqrt(2), 0]` for A035607), exact
  counterexamples for maps that do not preserve, and RealRooted declarations
  (schema `polytool.recurrence-basis-map-theorem.v1`);
  `certify_basis_map_lift` certifies `Q_n = T_P(R_n)` for all `n` by
  transporting the coefficient rule through the companion form of `T_P` and
  pseudo-dividing the input rule (`BasisMapDivisionCertificate`, replayed by
  `verify()`), and proves the coefficient roots in the source with the
  parent provers on affinely transported families
  (`polytool.recurrence-basis-map-lift.v1`, new `LiftOperation::BasisMap`);
  `find_basis_map_lifts` searches a `BasisMapLibrary` (by default the admitted
  maps of the OEIS catalog; `polytool.recurrence-basis-map-lifts.v1`), and
  `apply_basis_map_to_family` builds and certifies the rule of an image
  `T_P(R)` (`polytool.recurrence-basis-map-image.v1`): transported for a map
  of order one, computed from the transported system for higher order
  (`BasisMapImageConstruction::Eliminated`), fitted only when that finds
  nothing within its budget.
- `sequences::literature::peak_polynomials_bigint` and
  `peak_polynomials_iter` (permutations by peaks, A008303), and library
  recurrences returned as `(Recurrence, seeds, first_index)`
  (`literature::RecurrenceWithSeeds`): `peak_polynomials_recurrence`,
  `one_over_k_eulerian_symmetric_part_a_recurrence` and
  `one_over_k_eulerian_symmetric_part_b_recurrence` (second-order rules of the
  symmetric decompositions of the `1/k`-Eulerian polynomials, derived in the
  rustdoc).  With the `a_n` rule the symmetric-companion lift proves the
  symmetric parts real-rooted for `k = 2, 3`.
- Row-factor lifts: `recurrence::find_row_factor_lifts` detects
  `Q_n = F(x)^(alpha n + beta) S_n` from exact row factorizations, certifies it
  by substituting into the input relation (`RowFactorCertificate`, replayable),
  checks `F` once and proves the cofactors with the parent provers; report
  `polytool.recurrence-row-factor-lifts.v1` (new `LiftOperation::RowFactor`).
- Exact sign decisions for an integer index (`integer_range_sign`,
  re-exported at the crate root): `decide_index_polynomial_sign`,
  `decide_index_rational_function_sign` and `..._with_budget` decide
  `f(n) > 0`, `>= 0`, `< 0`, `<= 0`, `!= 0` or `= 0` for every integer `n` of
  an `IntegerRange` (`[first, last]` or `[first, +inf)`), `f = p / q` with
  exact rational coefficients.  The result (`IndexSignDecision`) is `holds`
  with a replayable `IndexSignCertificate` (Taylor-shift coefficients of one
  sign, or exact signs at sample indices between Sturm-isolated integer root
  windows; `verify()` replays it), `fails_at` with the smallest failing index
  and the exact value (`null` where the denominator vanishes), or `undecided`
  only over the `IndexSignBudget` (degree 64, 100000 bisection steps by
  default).  Huge roots are located by bisection, never by enumeration.
  Serialized with exact integer and rational strings.
- Common lift summary: `recurrence::LiftSummary` (operation and parameters,
  parent rule and seeds with gauge and index shift, identity status, root
  transfer with the parent and input root regions as `RootInterval`s, the
  transfer hypotheses with their exact sign decisions, parent proof status and
  route, the verdict of the shared `lift_verdict` (`proved`,
  `conditional_on_attestation`, `open`, or `route_failed`, which says nothing
  against the input), and `lean_replay`).  It is
  an additive field `lift_summary` of each lift in
  `polytool.recurrence-veronese-lifts.v1` (computed by
  `VeroneseSectionLift::lift_summary`), of
  `polytool.recurrence-veronese-certificate.v1`
  (`VeroneseSectionCertification::lift_summary`, an `Option` so that older
  records still deserialize; `verify()` checks it against
  `recompute_lift_summary`) and of each lift in
  `polytool.recurrence-symmetric-companion.v1`
  (`SymmetricCompanionLift::lift_summary`, likewise optional).  No existing
  field is renamed or removed, so the schema tags are unchanged; the new
  public field makes struct literals of `VeroneseSectionCertification` and
  `SymmetricCompanionLift` need `lift_summary`.
- Interlacing as a proved lift property: `LiftSummary::input_interlacing`
  (`LiftInterlacingClaim`: `proved` with a `LiftInterlacing` of kind, an
  `InterlacingOrientation`, `from_index` and a
  `LiftInterlacingJustification` with the theorem and its premises, or
  `not_proved` with the reason), set with
  `LiftSummary::with_input_interlacing`; the verdict does not depend on it.
  The parent provers record the interlacing their route concludes in the new
  field `VeroneseParentProof::interlacing` (window theorems, three-term,
  Liu--Wang, derivative-lag and Ma--Wang checks, all `P_(m-1) << P_m`).
  Transfer rules: the Hermite--Kakeya--Obreschkoff rule for operations free
  of `n` that preserve real-rootedness on the whole line (preserver images,
  Hadamard multiplier sequences such as `1/k!`, basis maps with a whole-line
  map theorem; weak, `unordered`); Athanasiadis--Wagner (arXiv:2404.12989,
  Corollary 5.6(b)) for the Veronese certifications (weak, directed); a
  common real-rooted factor for row factors with a constant exponent.
  `n`-dependent operators, gamma lifts and symmetric companions claim
  nothing.  The `lift_summary` JSON gains the required field
  `input_interlacing`, so summaries serialized before this change no longer
  deserialize (`lift_summary` itself is new in this release).
- The sign hypotheses of the lifts use the exact checker: the
  symmetric-companion hypotheses `a_n >= 0`, `0 <= delta_n <= c_n` (formerly
  unresolved for degree above one unless the Taylor shift had one-signed
  coefficients) and the coupling `epsilon_n != 0` are decided for every
  degree, a failure naming the smallest failing index; the Veronese gauge
  check and the root-domain helpers behind the Rolle and gamma certificates
  (index-polynomial signs, nonvanishing leaders, now via the gcd of the
  leader's coefficients) decide signs at every integer index instead of using
  a sufficient Taylor-shift test.
- `recurrence::certify_veronese_section_lift` certifies one candidate
  Veronese-section lift `Q_n = lambda(n) V_(r,j)(P_(n+s))` for inputs with
  index-dependent coefficients and derivatives: the parent rule is transported
  to the section matrices (`x -> X`, `D -> DD`, an algebra homomorphism), the
  input relation is pseudo-divided by it, and a zero remainder with exact
  checks at the finitely many exceptional bases proves the identity for every
  `n` (`VeroneseForwardDivisionCertificate`, replayed by `verify()`).  The
  parent is then proved by the periodic-product rule, a canonical paper
  theorem, an unconditional root-domain claim (e.g. Eulerian through the
  `f`-polynomial), or conditionally by a caller-supplied attestation.  Verdicts
  keep the identity, the parent proof and Lean replay (`not_attempted`) apart;
  schema `polytool.recurrence-veronese-certificate.v1`; MCP tool
  `certify_veronese_section_lift`.  Example: the even sections of the Eulerian
  polynomials are proved real-rooted this way.
- Real-rootedness theorems for recurrences with every hypothesis checked
  exactly (`recurrence::check_three_term_recurrence_theorem`,
  `recurrence::check_liu_wang_theorem`, report
  `RealRootednessTheoremReport`, schema
  `polytool.recurrence-real-rootedness-theorem.v1`, replayed by `verify()`).
  The three-term theorem covers `P_n = (a_n x + b_n) P_(n-1) - c_n P_(n-2)`
  (orthogonal-polynomial recurrences: simple real roots and strict
  interlacing when `a_n != 0` and `c_n / (a_n a_(n-1)) > 0`, with interlacing
  seeds), also after the coordinate change `y = 1/(x + v_0)` when the last
  coefficient is `c_n (x + v_0)^2`, and certifies an optional root region by
  a positive-ratio (chain-sequence) induction.  The Liu--Wang theorem (Adv.
  Appl. Math. 2007, Theorem 2.1 with `g = f'`) covers
  `P_n = a_n(x) P_(n-1) + b_n(x) P'_(n-1)` with every root kept in
  `(-inf, 0)`.  Hypotheses on rational functions of `n` go through the index
  sign checker; a failing hypothesis names the smallest failing index.  Both
  are new parent routes of `certify_veronese_section_lift`
  (`VeroneseParentProofRoute::{ThreeTermRecurrence, LiuWang}`, lift-summary
  routes `three_term_recurrence` and `liu_wang`), run after the root-domain
  route: the even sections of the Narayana polynomials (A001263, order-two
  parent, through `y = 1/(x - 1)`) and of the type-B Eulerian polynomials
  (A060187) are now proved real-rooted, and the Eulerian parent is proved by
  Liu--Wang when the root-domain route is switched off.
- The three-term and Liu--Wang checks move their induction start up to
  `RealRootednessTheoremOptions::maximum_start_delay` (default 3) rows past
  the supplied seeds when the theorem fails at the seeds; the skipped rows
  enter as seed rows checked exactly, the report records the delay in the new
  field `start_delay` (omitted from the JSON when zero) and `verify()`
  replays it.  The parent provers' root-domain route moves the parent's first
  index the same way.  The Eulerian and 1/k-Eulerian rules from `A_0 = 1`
  (where `A_1 = A_0`) are now proved without shifting the start by hand.
- Liu--Wang derivative-lag theorem:
  `recurrence::check_liu_wang_derivative_lag_theorem` (theorem
  `liu_wang_derivative_lag`) for second-order rules
  `f(n) P_n = c_0 P_(n-1) + c_1 P'_(n-1) + c_2 P_(n-2)`: strict interlacing of
  consecutive rows and every root `< 0` when the derivative and lag
  multipliers vanish at `0`, are nonpositive left of `0` (one of them
  negative there), `A_n(0) > 0`, the degree grows by one, and the two seed
  rows of the induction strictly interlace, each decided exactly for all
  `n`.  It mirrors the RealRooted row routes `halfGrowth`
  (`derivLag_strictInterl_of_nonnegCoeffs`) and `sign` (three-term rules),
  with nonnegative coefficients replaced by the root-location induction.
  When `x^e` divides the two seed rows it proves the rows divided by `x^e`
  (roots in `(-inf, 0]`, weak interlacing).  Proved: Delannoy A008288, deco
  A144438 (and the deco rows with `w = 2`), A046739 and derangement
  excedances A271697 (start moved two rows, gauge `x`); binomial Eulerian A046802 fails the lag sign at `n = 2`.  New
  report fields `lean_route` and `lean_declarations` (omitted when empty)
  record the mirrored RealRooted route; new parent route
  `VeroneseParentProofRoute::LiuWangDerivativeLag` (lift-summary route
  `liu_wang_derivative_lag`), run after Liu--Wang.  A theorem conclusion no
  longer extends the interlacing claim back over a zero seed row.
- Hadamard-factor lifts remove a division remainder with an invariant row
  annihilator of any first-order parent (an operator `A(n)` with an exact
  identity `nu A(n + 1) T(n) = B A(n)` under the parent step and exact row
  checks; new optional certificate fields `row_annihilator_step_multiplier`,
  `row_annihilator_step_quotient`, `row_annihilator_start_base`,
  `row_annihilator_restart_bases`), not only of product parents.  The Baxter
  triangle A056939 (Hoggatt rows with `m = 3`) is now proved end to end:
  factor `binom(n + 2, k + 1)`, hypergeometric parent proved by Liu--Wang.
- Hadamard-factor lifts also remove a division remainder for parents of order
  two or more, `g(n) P_(n+k) = sum_(i < k) T_i(n) P_(n+i)`, with an invariant
  row-annihilator system on the windows `(P_n, ..., P_(n+k-1))`: generators
  `A^(l)` of shift below `k` with exact identities
  `nu_l N(A^(l)) = sum_j B_(l,j) A^(j)` under the companion step
  `N(A) = g S A - A_(k-1)(n + 1) L`, exact base and restart row checks, and
  `mu r = sum_l Z_l A^(l)`.  New optional certificate field
  `row_annihilator_system` (types `HadamardAnnihilatorSystem`,
  `HadamardAnnihilatorGenerator`); the schema stays v1 and the first-order
  path is unchanged.  The generators are guessed from the parent rows by
  modular kernel computations and checked exactly.  This certifies the
  `binom(n + 1, k)` identities of the descent polynomials of 321-avoiding
  permutations and the peak polynomials of 123-avoiding permutations, and the
  Baxter triangle is now also proved through its second-order (Jacobi)
  parent.
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
