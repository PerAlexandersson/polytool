# Python bindings

The `polytool` Python package calls the exact polynomial, recurrence and
proof routines of the Rust library from Python and Jupyter (SageMath
included), without shelling out to the command-line tool.  It is a thin
adapter, like the CLI, the [MCP server](mcp.md) and the [web app](web.md): the
mathematics and the JSON schemas are those of the [library](library.md).

- Exact numbers only: coefficients and results are Python `int` and
  `fractions.Fraction`.  Floats are refused with a `TypeError`, never rounded.
- Polynomials are coefficient lists in ascending degree order (entry `i` is
  the coefficient of `t^i`), or `polytool.Polynomial` objects.
- Sequence positions are called *indices*: `first_index` labels the first
  polynomial `P_first_index`, and `from_index` selects where output starts.
- Reports convert to `dict` through the same stable JSON schemas as the CLI's
  `--json` output (tested for equality, see [Tests](#tests-and-ci)).

## Installation

The package is built with [maturin](https://www.maturin.rs/) from the
`python/` directory of the repository.  It needs Rust 1.85 or newer and
CPython 3.10 or newer; one `abi3` wheel serves every CPython version from 3.10
on.  Prebuilt wheels will be published on PyPI as `polytool-exact`
(`pip install polytool-exact`; the import name is `polytool`); until then,
install from source inside a virtual environment:

```sh
python -m venv ~/.venvs/polytool && . ~/.venvs/polytool/bin/activate
pip install "git+https://github.com/PerAlexandersson/polytool#subdirectory=python"
```

or, from a checkout:

```sh
cd python
pip install maturin
maturin build --release --out /tmp/polytool-wheels
pip install /tmp/polytool-wheels/polytool_exact-*.whl
```

`maturin develop --release` installs an editable build into the active
environment instead (it copies the compiled module next to the package
sources; that file is ignored by Git).  Set `CARGO_TARGET_DIR` to keep the
Rust build output out of the checkout.

## Tutorial

The examples below are a notebook session; they are checked by the test
suite (`python/tests/test_documentation.py`).

### Polynomials

```pycon
>>> import polytool
>>> from fractions import Fraction
>>> p = polytool.Polynomial([1, 11, 11, 1])
>>> p
Polynomial([1, 11, 11, 1])
>>> print(p)
1 + 11t + 11t^2 + t^3
>>> p.degree, p.coefficients
(3, [1, 11, 11, 1])
>>> q = polytool.Polynomial("1 + 2t")          # integer polynomial text
>>> q**2, p // q, p % q
(Polynomial([1, 4, 4]), Polynomial([Fraction(23, 8), Fraction(21, 4), Fraction(1, 2)]), Polynomial([Fraction(-15, 8)]))
>>> p(Fraction(1, 2)), p(q)                    # evaluate, compose
(Fraction(75, 8), Polynomial([24, 72, 56, 8]))
>>> p.shift(-1), p.dilate(2), p.reverse(degree=4)
(Polynomial([0, -8, 8, 1]), Polynomial([1, 22, 44, 8]), Polynomial([0, 1, 11, 11, 1]))
>>> polytool.Polynomial([10**30, Fraction(1, 3)])
Polynomial([1000000000000000000000000000000, Fraction(1, 3)])

```

In a notebook, polynomials and recurrences render as LaTeX (`_repr_latex_`).
Arithmetic accepts polynomials, `int` and `Fraction`; `/` divides by a
number, and `//`, `%`, `divmod` and `gcd` divide polynomials over the
rationals.  Further transforms: `derivative(order)`, `compose`,
`squarefree_part`, and coordinates in a basis (`bernstein_coordinates`,
`magic_coordinates` for `t^i (1 + t)^(d - i)`, `coordinates_in_basis`).

### Real roots and interlacing

Every function takes a `Polynomial`, a coefficient list or polynomial text,
and has a method of the same name on `Polynomial`.

```pycon
>>> polytool.is_real_rooted([1, 11, 11, 1])
True
>>> r = polytool.Polynomial([0, -2, -4, -1, 2, 1])   # t (t + 1)^2 (t^2 - 2)
>>> r.count_real_roots(), r.count_real_roots(distinct=True)
(5, 4)
>>> roots = r.isolate_real_roots()
>>> [(root.lower, root.upper, root.multiplicity) for root in roots[:3]]
[(Fraction(-46345, 32768), Fraction(-92655, 65536), 1), (-1, -1, 2), (0, 0, 1)]
>>> polytool.interlaces([8, -6, 1], [-15, 23, -9, 1])     # (t-2)(t-4) << (t-1)(t-3)(t-5)
True
>>> polytool.interlaces([-15, 23, -9, 1], [8, -6, 1]) is None
True
>>> polytool.interlaces_either_order([-15, 23, -9, 1], [8, -6, 1])
True

```

`isolate_real_roots` returns exact rational enclosures (width at most `1/1024`
by default, or `max_width`); an enclosure with `lower == upper` is an exact
root.  `interlaces(f, g)` is the directed weak test `f << g`
(`strict=True` for strict interlacing); `None` means that the degrees do not
admit the directed test in this order, not "does not interlace".

### Coefficient properties

```pycon
>>> polytool.gamma_coefficients([1, 11, 11, 1])
[1, 8]
>>> polytool.is_log_concave([1, 4, 1]), polytool.is_unimodal([1, 0, 1])
(True, False)
>>> polytool.properties([1, 4, 1])
{'polynomial': '1 + 4t + t^2', 'coefficients': [1, 4, 1], 'degree': 2, 'real_rooted': True, 'simple_roots': True, 'palindromic': True, 'gamma_positive': True, 'gamma_coefficients': [1, 2], 'unimodal': True, 'log_concave': True, 'ultra_log_concave': True}

```

`properties` returns an item of `polytool properties --json` (without
`index`), with exact numbers where the CLI writes decimal strings.  The
predicates `is_log_concave`, `is_ultra_log_concave`, `is_unimodal`,
`is_palindromic`, `is_gamma_positive` and `has_simple_roots` are exact;
rational input is decided on its smallest positive integer multiple, which
none of them distinguishes from the input.  `discriminant` and `resultant`
need integer coefficients.

### Recurrence search

`find_recurrence(rows, **options)` runs the adaptive search of
`polytool recurrence` (same defaults) and returns a `Recurrence`, or `None`
when the configured search space contains none.

```pycon
>>> rows = polytool.oeis.generate("A001263", 8)      # Narayana triangle
>>> rows[:4]
[[1], [1, 1], [1, 3, 1], [1, 6, 6, 1]]
>>> narayana = polytool.find_recurrence(rows)
>>> narayana
<Recurrence (1 + n) P(n) = (-1 + 2n - t + 2nt) P(n-1) + (2 - n - 4t + 2nt + 2t^2 - nt^2) P(n-2) with 2 initial polynomials from index 1, fitted>
>>> print(narayana.latex())
% Candidate fitted to finite rows; validity for all n is unproved.
$P_{1} = 1$
$P_{2} = 1 + t$
$(n + 1) P_{n} = (2n - 1)(1 + t) P_{n-1} - (n - 2)(1 - t)^{2} P_{n-2} \quad (n \geq 3)$
>>> narayana.generate(2, from_index=9)
[Polynomial([1, 36, 336, 1176, 1764, 1176, 336, 36, 1]), Polynomial([1, 45, 540, 2520, 5292, 5292, 2520, 540, 45, 1])]
>>> rows_lazily = narayana.iterate()
>>> next(rows_lazily), next(rows_lazily)
(Polynomial([1]), Polynomial([1, 1]))

```

The keyword options are the CLI search flags with spelled-out names:

| Python keyword | CLI flag |
| --- | --- |
| `first_index` (default 1), `skip_prefix` | `--first-index`, `--skip-prefix` |
| `min_recurrence_length`, `max_recurrence_length` | `--min-recurrence-length`, `--max-recurrence-length` |
| `min_variable_degree`, `max_variable_degree` | `--min-variable-degree`, `--max-variable-degree` |
| `min_index_degree`, `max_index_degree` | `--min-index-degree`, `--max-index-degree` |
| `min_derivative_order`, `max_derivative_order` | `--min-derivative-order`, `--max-derivative-order` |
| `inhomogeneous`, `min_/max_inhomogeneous_variable_degree`, `min_/max_inhomogeneous_index_degree` | `--try-inhomogeneous`, `--min/max-inhomogeneous-variable-degree`, `--min/max-inhomogeneous-index-degree` |
| `denominator`, `max_denominator_variable_degree`, `max_denominator_index_degree` | `--try-denominator`, `--max-denominator-variable-degree`, `--max-denominator-index-degree` |
| `alternating_sign` | `--try-alternating-sign` |
| `laurent` (default `True`), `max_laurent_power` | `--no-laurent`, `--max-laurent-power` |
| `min_margin`, `fit_extra_rows`, `verify` (default `True`), `modular_prefilter` | `--min-equation-margin`, `--fit-extra-rows`, `--no-held-out-rows`, `--(no-)modular-prefilter` |
| `max_candidates`, `max_work` | `--max-candidates`, `--max-work` |
| `lean_name` (default `None`) | `--lean-name` |

A budget that runs out raises `polytool.BudgetExhaustedError` (with the
attributes `candidates_tried`, `candidates_considered` and `work_units`).
Long computations release the GIL but cannot be interrupted with Ctrl-C; use
`max_work` for a bounded run.

A found recurrence is a fit to finitely many rows (its `search` metadata
records the fit).  Proofs below concern the sequence that the recurrence and
its initial polynomials define.

### JSON records

A `Recurrence` is one record of the schema `polytool.recurrence.v1`, the
output of `polytool recurrence --json`; `to_json(indent=2)` reproduces that
output byte for byte when the search options match. Fitted records include
the checked `presentation`, `normalized_recurrence` and the minimal initial-row
prefix, extended through any left-hand singularity. All supplied rows are used
to check the boundary, including retained rows before the fitting window.
Supplied coefficient lists retain their zero padding in the JSON seed rows,
matching CLI input; `Polynomial` objects already have trimmed coefficients.

Lean metadata is opt-in: `find_recurrence(rows, lean_name="Rows")` adds `lean`
on success or `lean_reason` when the rule cannot be rendered. With the default
`lean_name=None`, neither field is present. This matches CLI JSON with and
without `--lean-name Rows`, respectively. Invalid names raise `ValueError`.
Lean output is checked by exact semantic read-back; it is not a proof that a
finite fit holds for every index. JSON/dictionary round trips preserve all
requested metadata; reindexing clears Lean output because its coefficients
and index shifts must be rendered again.

```pycon
>>> named = polytool.find_recurrence(rows, lean_name="Rows")
>>> "def Rows : ℕ → ℝ[X]" in named.to_dict()["lean"]
True
>>> "lean" in narayana.to_dict(), "lean_reason" in narayana.to_dict()
(False, False)

```

```pycon
>>> record = narayana.to_dict()
>>> record["schema"], record["first_index"], record["initial_polynomials"]
('polytool.recurrence.v1', 1, [['1'], ['1', '1']])
>>> polytool.Recurrence.from_json(narayana.to_json()) == narayana
True
>>> legendre = polytool.Recurrence.from_json({
...     "schema": "polytool.recurrence.v1", "first_index": 0,
...     "initial_polynomials": [["1"], ["0", "1"]],
...     "recurrence": {
...         "terms": [
...             {"offset": 1, "deriv_order": 0, "sign": "none",
...              "coeff": {"coeffs": [["0", "-1"], ["0", "2"]]}},
...             {"offset": 2, "deriv_order": 0, "sign": "none",
...              "coeff": {"coeffs": [["1"], ["-1"]]}}],
...         "denominator": {"coeffs": [["0"], ["1"]]}}})
>>> legendre.generate(3)
[Polynomial([1]), Polynomial([0, 1]), Polynomial([Fraction(-1, 2), 0, Fraction(3, 2)])]

```

In a coefficient array `coeffs[i][j]` is the coefficient of `n^i t^j`, an
exact decimal string; see [cli-recurrence.md](cli-recurrence.md) for the
schema.  `reindexed(first_index)` relabels the rows, and `match_theorems()`
returns the conditional theorem report of `polytool recurrence-theorems`.

### Real-rootedness proofs

`prove_real_rooted` (also a `Recurrence` method) runs the general theorem
search of `polytool recurrence-prove` ([proof-discovery.md](proof-discovery.md)):
product induction, three-term (Favard), Liu--Wang, Liu--Wang derivative lag,
Ma--Wang interval, Liu--Wang three-term, iterated-preserver and
hypergeometric-row criteria, each hypothesis decided exactly for every `n`,
then the tail after rows that are zero or not real-rooted, certified
alternative representations and catalog lifts to a family of the bundled
OEIS catalog (`max_lift_candidates=6`, `0` disables them), and last the
criteria in exact affine coordinates `x = a u + b` derived from the rule
(`max_affine_coordinates=8`, `0` disables them; `report.affine_coordinate` is
the conclusion mapped back to `x`, and `report.route` is then
`"affine_coordinate:<theorem>"`).

```pycon
>>> report = narayana.prove_real_rooted()
>>> report
<ProofReport proved by three_term_recurrence>
>>> report.status, report.route
('proved', 'three_term_recurrence')
>>> print(report.certificate["conclusion"]["statement"])
every P_n, n >= 1, is nonzero and real-rooted, with simple roots from n = 1 on; Q_(n-1) strictly interlaces Q_n in y for n >= 2
>>> report.verify()
True
>>> polytool.prove_real_rooted(legendre, root_region=(-1, 1)).certificate["conclusion"]["root_region"]
{'lower': {'value': '-1', 'closed': False}, 'upper': {'value': '1', 'closed': False}}

```

`report.attempts` lists every theorem tried, with its hypotheses;
`certificate` is the first proved one.  `real_rooted_from_index` is the
first index `N` with every row `P_n`, `n >= N`, proved nonzero and
real-rooted.  A residue-class proof may leave whole residue classes of zero
rows (for `P_n = t^2 P_(n-2)` with `P_1 = 0`, every odd row): then
`real_rooted_from_index` is `None` and
`real_rooted_outside_zero_residue_classes` is a `dict` with `from_index` and
`zero_residue_classes` (`period`, `residues` of `n`).  `status ==
"not_proved"` is not a disproof.  `to_dict()` is the `polytool.recurrence-proof-discovery.v1` report
of `polytool recurrence-prove --json`, `ProofReport.from_json` reads one
back, and `verify()` replays the whole search, raising
`polytool.VerificationError` on any mismatch.  A record proved by a catalog
lift has `route` `"catalog_lift:<route>"` and its conclusion as `certificate`
and `catalog_lift`; a lift resting only on the catalog's `lean_proved` flag
gives `status == "conditional_on_attestation"` (`is_proved` is `False`) with
the conditional claim in `catalog_lift`.

### Degree laws and leading coefficients

```pycon
>>> floor_third = polytool.Recurrence.from_json({
...     "schema": "polytool.recurrence.v1", "first_index": 0,
...     "initial_polynomials": [["1"], ["1"], ["1"]],
...     "recurrence": {"terms": [
...         {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["1"]]}},
...         {"offset": 3, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["0", "1"]]}}]}})
>>> degree = floor_third.certify_leading_coefficient()
>>> degree.is_proved, degree.route, degree.degree_at(31), degree.sign_at(31)
(True, 'cone', 10, 'positive')
>>> degree.verify()
True

```

`certify_leading_coefficient` is `polytool recurrence-degree`: a degree law
`deg P_n = D(n)`, periodic by residue class, with the sign of every leading
coefficient, certified for all `n` by a route of the RealRooted tactic
`rr_row_leadingCoeff_pos` (`degree_law`, `holds_from`, `conclusion`,
`exact_leading_coefficient`; schema `polytool.recurrence-leading-coefficient.v1`).

### The OEIS catalog

```pycon
>>> "A008292" in polytool.oeis.ids()
True
>>> info = polytool.oeis.info("A008292")
>>> info["name"]
'Triangle of Eulerian numbers T(n,k) (n >= 1, 1 <= k <= n) read by rows.'
>>> info["first_index"], info["real_rootedness"]["status"]
(1, 'lean_proved')
>>> polytool.oeis.generate("A008292", 3, from_index=4)
[[1, 11, 11, 1], [1, 26, 66, 26, 1], [1, 57, 302, 302, 57, 1]]
>>> eulerian = polytool.oeis.recurrence("A008292")
>>> print(eulerian)
P(n) = (1 - t + nt) P(n-1) + (t - t^2) P'(n-1)
>>> eulerian.prove_real_rooted().route
'liu_wang'

```

`generate` returns rows exactly as the OEIS lists them (lists of `int`,
including padding zeros of a triangle); `from_index` is the displayed OEIS
row number.  `info` has the keys of `polytool oeis info --json`, with
`first_index` for the CLI's `first_row`, and `recurrence` returns the bundled
recurrence indexed by the displayed rows.  See [cli-oeis.md](cli-oeis.md) for
the catalog's statuses, real-rootedness flags and evidence records.

## API summary

| Name | Returns |
| --- | --- |
| `Polynomial(coefficients)` | exact polynomial; arithmetic, evaluation, transforms |
| `is_real_rooted`, `count_real_roots`, `isolate_real_roots` | `bool`, `int`, list of `RealRootEnclosure` |
| `interlaces`, `interlaces_either_order` | `bool` or `None`, `bool` |
| `is_log_concave`, `is_ultra_log_concave`, `is_unimodal`, `is_palindromic`, `is_gamma_positive`, `has_simple_roots` | `bool` |
| `gamma_coefficients`, `discriminant`, `resultant`, `properties` | exact values, `dict` |
| `find_recurrence` | `Recurrence` or `None` |
| `Recurrence`: `generate`, `iterate`, `latex`, `to_json`, `to_dict`, `from_json`, `reindexed`, `match_theorems` | rows, LaTeX, JSON |
| `prove_real_rooted` | `ProofReport` (`status`, `route`, `certificate`, `attempts`, `real_rooted_from_index`, `real_rooted_outside_zero_residue_classes`, `catalog_lift`, `verify`) |
| `certify_leading_coefficient` | `LeadingCoefficientCertificate` (`degree_at`, `sign_at`, `verify`) |
| `oeis.ids`, `oeis.info`, `oeis.generate`, `oeis.recurrence` | catalog access |
| `PolytoolError` (a `ValueError`), `VerificationError`, `BudgetExhaustedError` | exceptions |

`help(polytool.find_recurrence)` and the other docstrings give the full
contracts.

## Tests and CI

The pytest suite in `python/tests/` covers big integers and fractions in both
directions, JSON and pickle round trips, the docstring and tutorial examples,
and parity with the CLI: with `POLYTOOL_CLI` set to a `polytool` binary built
from the same source, `test_cli_parity.py` compares `properties`,
`real-rooted`, `recurrence`, `recurrence-prove`, `recurrence-degree` and
`oeis` JSON output with the Python results.

```sh
cd python
maturin build --release --out /tmp/polytool-wheels
pip install --force-reinstall /tmp/polytool-wheels/polytool_exact-*.whl pytest
cargo build --release --bin polytool --manifest-path ../Cargo.toml
POLYTOOL_CLI=<target directory>/release/polytool pytest tests
```

The GitHub workflow `.github/workflows/python.yml` builds the wheels on Linux,
macOS and Windows, runs the suite on Linux with the oldest and newest
supported CPython, and checks the crate (formatting, Clippy, Rust unit tests).
Publishing to PyPI is a separate, manually approved step reserved for the
repository owner.
