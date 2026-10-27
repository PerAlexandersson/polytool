# Recurrence transformations and implicit relations

Part of the [CLI reference](cli.md); continues
[recurrence search](cli-recurrence.md).

### Coordinate transformations and recurrence transport

`polytool recurrence-transform` applies an exact change of polynomial
coordinates to a sequence, transports its defining recurrence into the new
coordinates, and can search independently there. It reads one JSON request on
stdin and prints a `polytool.recurrence-transform.v1` report; `--verify`
replays a report exactly from its embedded request (exit status 2 on a
mismatch or invalid input), and `--compact` prints compact JSON.

| `kind` | Row map | Precondition / inverse image |
| --- | --- | --- |
| `multiply_power {power}` | `Q_n = x^r P_n` | inverse needs `x^r` to divide `Q_n` |
| `divide_power {power}` | `Q_n = P_n / x^r` | `x^r` divides `P_n` |
| `reversal {schedule}` | `Q_n = x^(d_n) P_n(1/x)` | `deg P_n <= d_n` |
| `affine {a, b}` | `Q_n = P_n(a x + b)`, `a != 0` | none |
| `mobius {a, b, c, e, schedule}` | `Q_n = (c x + e)^(d_n) P_n((a x + b)/(c x + e))`, `a e - b c != 0` | `deg P_n <= d_n` |
| `f_polynomial {schedule}` | `Q_n = (1 + x)^(d_n) P_n(x/(1 + x))` | `deg P_n <= d_n` |
| `gamma {schedule}` | `P_n = sum_i Q_(n,i) x^i (1 + x)^(d_n - 2i)` | `x^(d_n) P_n(1/x) = P_n`; inverse needs `deg Q_n <= d_n/2` |
| `symmetric_decomposition {schedule}` | `P_n = a_n + x b_n`, state `(a_n, b_n)` | `deg P_n <= d_n`; inverse needs `a_n`, `b_n` symmetric at `d_n`, `d_n - 1` |

Rational parameters are strings such as `"3/2"`. A schedule
`{"first_index": n0, "first_degree": d0, "step": h, "last_index": n1}` declares
`d_n = d0 + h (n - n0)` on the index domain `n0 <= n (<= n1)`, so
`d_(n+1) = d_n + h`; a decreasing schedule needs `last_index`, and negative
degrees, empty domains and rows outside the domain are rejected. The declared
bound is used even when a row has lower degree (`f_polynomial` calls
`decomposition::f_polynomial` with `d_n`, never the observed degree), so
missing top coefficients and endpoint factors are preserved exactly.

The source is a defining `polytool.recurrence.v1` rule with seeds
(`source_recurrence`), a bundled OEIS family (`oeis_id`), or explicit `rows`
with `row_origin`; with a source, `row_count` rows are generated (default 10).
Every row is mapped, the inverse reconstructs it exactly
(`reconstruction_verified`), and a row violating a precondition rejects the
request. Each condition is reported as `certified_all_indices` when the source
rule proves it for the whole domain (the transported rule is denominator-free
and the seeds qualify; for gamma, the reversal-transported rule equals the
source rule) and otherwise as `checked_supplied_rows`, which is finite
evidence.

A rule `f(n,x) P_n = sum s(n) c(n,x) D^k P_(n-r) + g(n,x)` is transported
exactly. For a Möbius map with `phi = (a x + b)/ell`, `ell = c x + e`,
`Delta = a e - b c`, a first-order rule `R P_(n+1) = A P_n + B P_n'` becomes

```text
R Q_(n+1) = ell^h (A(phi) - d_n c ell B(phi)/Delta) Q_n + ell^(h+2) B(phi)/Delta Q_n'
```

Higher derivatives and lags follow from the product and chain rules; the
identity is multiplied by the power of `ell` that clears every denominator and
a common power of `ell` is then divided out exactly, never evaluated at the
pole. An additive term `g` is transported only when `ell^(d_n)` does not
depend on `n` (`h = 0`, or constant `ell` with `ell^h = 1`). The symmetric
decomposition turns a first-order rule without `(-1)^n` factors into the
complete two-component shared-denominator vector system for `(a_n, b_n)`.

Gamma coordinates `P_n = (1 + x)^(d_n) G_n(z)`, `z = x/(1 + x)^2`, are a
quadratic substitution and are transported differently
(`recurrence::transport_gamma`). With `u = 1 + x` and `z' = (1 - x)/u^3`,
`D^k [u^d G(z)] = sum_i u^(d - k - 2i) T_(k,i) G^(i)(z)` with `T_(0,0) = 1` and
`T_(k+1,i) = (d - k - 2i) T_(k,i) + u T_(k,i)' + (1 - x) T_(k,i-1)`. The rule
is divided by `u^(d_m)` and multiplied by the power of `u` that clears every
denominator; then each coefficient splits at an index `N` into its symmetric
and antisymmetric parts, `E = u^N A(z) + (1 - x) u^(N-1) B(z)`. Since
`w = (1 - x)/(1 + x)` satisfies `w^2 = 1 - 4z` and `x -> 1/x` fixes `z` and
negates `w`, the identity is equivalent to two relations: the even part (all
`A`) and the odd part (all `B`). When they are proportional (for a rule that
is invariant under reversal at `d_n`) the part with a nonzero leader is the
transported scalar rule; otherwise the transport error names both parts (the
odd part is a companion constraint, and no single scalar rule carries the
source rule). The symmetry condition is certified for every index when the
reversal-transported rule defines the same row map as the source rule, its
left factor is provably nonzero at every target, and the seeds are symmetric.
For the canonical families (row index 1 of degree 0, `d_n = n - 1`):

| Family | Gamma rule |
| --- | --- |
| Eulerian A008292 | `G_m = (1 + 2(m - 2) z) G_(m-1) + z (1 - 4z) G'_(m-1)` |
| Narayana A001263 | `(m + 1) G_m = (2m - 1) G_(m-1) - (m - 2)(1 - 4z) G_(m-2)` |
| type B Eulerian A060187 | `G_m = (1 + 4(m - 2) z) G_(m-1) + 2z (1 - 4z) G'_(m-1)` |

A rule such as `P_m = (1 + 2x) P_(m-1)` has the independent parts
`G_m = (3/2) G_(m-1)` and `0 = -(1/2) G_(m-1)` and is not transported.

On symmetric rows the independent pair is still an exact description: the
report states it as `transported.system`, the even part and its odd companion
as implicit relations of the gamma rows (written `Q[n+i]`, base index `n`, see
"Implicit differential relations" below), each with its first base and the
evidence `transported_identity` (conditional on the symmetry condition in
`conditions`). Together they are equivalent to the source rule on symmetric
rows; a rule with an additive term or `(-1)^n` factors has no homogeneous
form and keeps the transport error. For `(1 + x)^n` generated by the rule
`P_m = (1 + 2x - m x) P_(m-1) + (x + x^2) P'_(m-1)`, which is not invariant
under reversal at `d_n = n`, the gamma rows are `G_n = 1` and the system is
`Q[n+1] - Q[n] = 0` with the row equation `-x D Q[n] = 0`
(`GammaTransport::relation_system`).

The transported rule (`evidence: "transported_identity"`, a replayable
`polytool.recurrence.v1` record with transformed seeds) is an exact
consequence of the supplied source rule and is kept separately from the
optional independent search (`"search": {}` or explicit bounds,
`evidence: "finite_fit"`). Canonical Eulerian rows A008292 (`E_1 = 1`, degree
`n - 1`) under the `f`-polynomial map give `Q_1 = 1`,
`Q_(n+1) = (1 + 2x) Q_n + x (1 + x) Q_n'`:

```sh
echo '{"oeis_id":"A008292","row_count":8,"search":{},
  "transformation":{"kind":"f_polynomial",
    "schedule":{"first_index":1,"first_degree":0,"step":1}}}' \
  | polytool recurrence-transform > report.json
polytool recurrence-transform --verify < report.json
```

In Rust, `RecurrenceTransformation::{map_row, inverse_row, transport,
record}` act on single rows and rules, and
`recurrence::{transform_recurrence, verify_recurrence_transform}` build and
replay reports. The MCP tools `transform_recurrence` and
`verify_recurrence_transform` take `{"request": {...}}` and
`{"report": {...}}`. Root-location information is opt-in (next section).

### Root-domain information for transformed results

Algebraic discovery never needs an interval or a real-rootedness proof. An
optional `"root_domain": {}` (or explicit fields) in a `recurrence-transform`
request adds a `root_domain` section to the report; requests and reports
without it are unchanged. The section keeps apart

- `implications`: theorem adapters for the map, each with its `direction`
  (`equivalent`, `source_to_target`, `target_to_source`), hypotheses with a
  status (`certified_all_indices`, `checked_supplied_rows`, `failed_on_rows`,
  ...), degree bounds, coefficient-sign hypotheses, endpoint conventions, zero
  cases, index domain and, for one-directional maps, why the converse fails;
- `claims` about the source rows `P_n` or the transformed rows, each with an
  `evidence` label: `finite_exact_observation` (exact root isolation of the
  listed rows only), `supplied_assumption` (from the request, checked on the
  rows), `proposed_interval` (a candidate window), `completed_recurrence_induction`
  (a closed induction over the transported rule for every index) and
  `certified_family_invariant` (transferred to every index by a map theorem
  whose conditions are certified); unproved dependencies are listed in
  `conditional_on`;
- `recurrence_theorem`: the match of the transported rule against Rolle's
  theorem for derivatives of products, `Q_m = s(m) D^j (g Q_(m-1))`, with every
  interval attempt (bounded by `max_candidate_intervals`);
- `observations`: exact rational root enclosures with multiplicities for the
  first `max_observed_rows` rows.

| Adapter | Statement | Direction |
| --- | --- | --- |
| `monomial_factor` | `Q_n = x^r P_n`; only the multiplicity of the root 0 changes | equivalent |
| `real_mobius_equivalence` (`reversal`, `mobius`, `f_polynomial`), `real_affine_equivalence` | real Möbius maps with `a e - b c != 0` and `deg P_n <= d_n`; the pole `-e/c` is a root of multiplicity `d_n - deg P_n`, roots at `a/c` lower `deg Q_n`; intervals and rays map exactly | equivalent |
| `nonnegative_coefficient_window` | nonnegative coefficients put the real roots of `P_n` in `(-inf, 0]`; for the `f`-polynomial every real root of `Q_n` lies in `[-1, 0]` (`-1` records missing top degree, `0` the valuation) | proposes a window only |
| `gamma_quarter_bound` | `P_n` is real-rooted iff every root of the gamma polynomial is real and at most `1/4` (`G = 1 - z` at `d = 2` gives `1 + x + x^2`) | equivalent |
| `hermite_kakeya_obreschkoff_symmetric_pair` | real-rooted interlacing `a_n`, `x b_n` give real-rooted `P_n = a_n + x b_n` | components to `P_n` only |
| `rolle_derivative_of_product` | roots of `g` and `Q_(m-1)` in an interval `I` put the roots of `s(m) D^j(g Q_(m-1))` in `I` (zero allowed); also `g(x)(c(m) + l(m) x)` and intervals `I_n` with `I_(m-1) ⊆ I_m` | recurrence step |

Real-rooted components alone do not suffice for the symmetric pair
(`a = (1 + x)^3`, `b = x` give `1 + 3x + 4x^2 + x^3`), and its converse fails
(`(1 + 2x)^2` has `a = 1 + x + x^2`). Invertible maps need not preserve
real-rootedness both ways: for `P = (xQ)'`, real-rooted `Q` gives real-rooted
`P`, but `(1 + x)^2` has the antecedent `1 + x + x^2/3`.

The induction tries supplied `candidate_intervals`, then the default window,
then the hull of the roots of `g` and the seeds, checking `g` and every seed
exactly. Canonical Eulerian rows A008292 under the `f`-polynomial map become
`Q_(n+1) = (x(1 + x) Q_n)'`, `Q_1 = 1`; the induction closes on `[-1, 0]`
(`completed_recurrence_induction`), and since `deg P_n <= d_n` is certified
and `deg Q_n = d_n` follows from the degree formula, the Möbius equivalence
gives real-rooted Eulerian polynomials with roots in `(-inf, 0]`
(`certified_family_invariant`):

```sh
echo '{"oeis_id":"A008292","row_count":8,"root_domain":{},
  "transformation":{"kind":"f_polynomial",
    "schedule":{"first_index":1,"first_degree":0,"step":1}}}' \
  | polytool recurrence-transform > report.json
polytool recurrence-transform --verify < report.json
```

`root_domain` fields: `candidate_intervals`, `assumptions`, `default_window`
(default `true`), `max_candidate_intervals` (8), `max_observed_rows` (32),
`max_observed_degree` (128). An interval is `{"lower": {...}, "upper": {...}}`
with endpoints `{"value": "-1", "per_index": "0", "closed": true}`; an omitted
end is infinite (a ray) and `value + per_index * n` depends on the index (for
`P_n = x + n` the roots `-n` escape while the `f`-polynomial roots stay in
`(-1, 0]`). An assumption is `{"subject": "source" | "transformed",
"component": 0, "property": "real_rooted" | "real_rooted_in" | "real_roots_in",
"interval": ..., "first_index": ..., "last_index": ..., "note": ...}`; it is
checked exactly on the observed rows (`assumption_checks`), transferred through
equivalences with its label in `conditional_on`, and replayed with the report.
Finite exact checks are never promoted to statements about every index.

Index-dependent candidate intervals `[l_0 + l_1 n, u_0 + u_1 n]` are tried in
the induction. The step needs `I_(m-1) ⊆ I_m` (`l_1 <= 0`, `u_1 >= 0`), `g`
real-rooted in the first `I_m`, and, for a factor `g(x)(c(m) + l(m) x)` with
an index-dependent linear part, a certified sign of `l(m)` and the root
`-c(m)/l(m)` in `I_m` for every `m`; these polynomial inequalities in `m` are
certified when every Taylor coefficient at the first target has the right
sign, and an attempt that cannot be certified fails. Under `Q_n(x) = P_n(-x)`
the rule `Q_m = (x (x + m) Q_(m-1))'` closes on `[-n, 0]`, and the affine map
gives `P_n` real-rooted in `[0, n]`. Index-dependent intervals map exactly
through affine maps; through a Möbius map with a pole their hull over the
index domain is mapped and the claim says so. In gamma coordinates a completed
induction for `G_n` transfers to `P_n` only when every interval lies in
`(-inf, 1/4]`: `P_m = (1 + x^2) P_(m-1)` gives `G_m = (1 - 2z) G_(m-1)`,
proved real-rooted at `1/2`, and no conclusion for `P_n`.

`polytool root-enclosures` (MCP `exact_root_enclosures`) is the exact interval
interface: for `{"polynomials": [["0", "0", "1", "3", "3", "1"]],
"max_width": "1/1000", "interval": {...}}` it returns enclosures with
multiplicities (`lower == upper` is an exact rational root; otherwise exactly
one root lies strictly inside and neither endpoint is a root), degree,
valuation, `real_rooted` and exact interval verdicts. The MCP `real_roots`
midpoints are representatives only. In Rust use
`real_rootedness::{isolate_real_roots_rational_coeffs, RealInterval}`; fixed-cap diagonal
multiplier tests use `recurrence::jensen_fixed_cap_check` (Jensen polynomial
`T[(1 + x)^N]`, Schur--Malo composition; evidence `fixed_cap`, never a proof
for larger degrees).

### Implicit differential relations, diagonal maps and reconstruction

Some coordinate changes turn a successor rule into an implicit relation with
derivatives on the newest row, or into a differential equation within one row.
An implicit relation is `sum c(n,x) D^d Q[n+i] = 0` at base index `n`
(JSON `{"terms": [{"shift": i, "deriv_order": d, "coeff": {"coeffs": ...}}]}`),
and a coefficient condition is `sum c(n) [x^k] Q[n+i] = 0`
(`{"terms": [{"shift": i, "coefficient_index": k, "coeff": [...]}]}`).
`polytool recurrence-implicit` maps a sequence with a defining rule
(`source_recurrence` or `oeis_id`) by one of two registry maps, which
`recurrence-transform` also accepts for row mapping:

| `kind` | Row map | Inverse |
| --- | --- | --- |
| `diagonal {weight}` | `Q_n = T_n P_n`, `T_n(x^k) = lambda(n,k) x^k` | divide coefficient `k` by `lambda(n,k)` |
| `derivative_antecedent {g}` | `Q_n` with `P_n = D^m(g Q_n)/m!`, `m = deg g` | `P_n = D^m(g Q_n)/m!` |

Weights are `factorial_power {power}` (`1/(k!)^m`), `binomial {m}`
(`binom(k+m,m)`, so `Q = W_m P` with `W_m = D^m x^m/m!`), `inverse_binomial {m}`
(`P = W_m Q`) and `product {scale, factorial_exponent, factors}` for
`scale (k!)^p prod f_j(n,k)^(e_j)`. A weight factor with a certified zero
(for example `k - 2`) is rejected; a factor whose zero set depends on `n`
is reported as `checked_supplied_rows` in `map_conditions`, and the
transported components then carry the same finite evidence. The antecedent
map is triangular with diagonal `lc(g) binom(k+m,m)`, so antecedents are found
by exact back substitution.

With `theta = x D`, a relation multiplied by `x^J` is written as
`sum x^s a(n,theta) S^i` and transported by
`x^s a(n,theta) S^i -> x^s a(n,theta) lambda(n,theta+s)/lambda(n+i,theta) S^i`.
One common multiplier `H(n,theta) = lcm b_(i,s)(n, theta - s)` clears the
denominators of every column (`H(theta) x^s = x^s H(theta+s)`). For the
original residual `Z_n`, `Z_n = 0` is equivalent to `H(n,theta) Z_n = 0`
together with `[x^k] Z_n = 0` at every `k >= 0` with `H(n,k) = 0`. These
exceptional indices are certified uniformly in `n` (exact nonnegative integer
roots of an index-free factor, or a one-signed shifted factor); otherwise the
clearing record says `checked_supplied_rows`, lists the factor and the zeros
seen on the supplied rows, and makes no equivalence claim. The derivative
antecedent transports a relation by substitution, `L T_g Q = 0`.

The report keeps evidence apart. `transports` holds each transported relation
(`transported_identity`), its clearing certificate and, when `P_n = T Q_n` with
a differential operator `T` (`W_m`, `T_g`), a reconstruction certificate
`L o T = sum_nu H_nu o S^nu o M` whose multipliers `H_nu = sum_d h_d(n,x) D^d`
record their differentiations. `description` is a generating description:
seeds, relations, `exceptional_index_conditions` and supplied
`boundary_conditions`. `assessment` gives the status (`equivalent_description`,
`certified_consequence`, `finite_evidence`, `inconsistent_seeds`) and a
uniqueness argument that is separate from annihilation. The argument picks a
relation of positive leader shift and compares coefficients of a difference
of two solutions. In the degree grading the top coefficient is multiplied by
`psi(n,k)`. In the valuation grading the lowest coefficient is multiplied by
`phi(n,k)`, and every certified exceptional index needs a condition that pins
`[x^k]` of the new row. Singular leaders and uncovered indices are listed as
reasons. `additional_relations` stores supplied or transported row equations
without treating them as successor rules. `"search": {}` fits a successor rule
to `Q` (`finite_fit`). The fit becomes `certified_by_reconstruction` when the
identity `L o T = sum H_nu S^nu M` holds, the source rule determines `P` from
its seeds, the generated `T Q*` agrees with the source seeds, and `T` is
injective.

Canonical binomial rows A007318 divided coefficientwise by `k!` give
`Q'_(n+1) = Q'_n + Q_n` with the exceptional condition
`Q_(n+1)(0) = Q_n(0)` and `Q_0 = 1`, an equivalent description; without the
constant condition (`--assess` on the edited description) the status drops to
`certified_consequence`. The row equation `x Q'' + (1 + x) Q' - n Q = 0` can be
stored through `target_relations` or transported from `(1 + x) P' - n P = 0`
through `source_relations`. Canonical A065826 divided by `k + 1` gives the
Eulerian rows A008292, and the fitted Eulerian rule is proved by
reconstruction:

```sh
echo '{"oeis_id":"A065826","row_count":9,"search":{},
  "transformation":{"kind":"diagonal","weight":{"kind":"inverse_binomial","m":1}}}' \
  | polytool recurrence-implicit > report.json
polytool recurrence-implicit --verify < report.json
```

An optional `"root_domain": {}` (or the fields of `recurrence-transform`)
adds a `root_domain` section with its own request. Adapters, each with its
direction and exact hypotheses: `polya_schur_multiplier_sequence` for an
index-free weight `scale (k!)^p prod f_j(k)^(e_j)` (`P_n` real-rooted gives
`Q_n` real-rooted when `p <= 0`, every `e_j > 0` and every `f_j` has only real
roots, all `<= 0`: Laguerre's theorem, `1/k!` and products; the reverse when
`1/lambda` qualifies, as for `inverse_binomial`), and
`rolle_derivative_antecedent` (`Q_n` real-rooted with `g` real-rooted gives
`P_n = D^m(g Q_n)/m!` real-rooted with roots in the hull of the roots of `g`
and the interval of `Q_n`; never the converse). The defining rule of `P_n`
and a successor rule of `Q_n` certified by reconstruction are matched against
Rolle's theorem (`source_rule_theorem`, `fitted_rule_theorem`), and a
completed induction is transferred only in the adapter's direction. An
index-free weight is also tested by `jensen_fixed_cap` at the largest
observed degree (`fixed_cap` evidence). Binomial rows divided by `k!` thus
give real-rooted `sum_k binom(n,k) x^k/k!` for every `n`.

`--verify` replays the report exactly and re-checks every reconstruction
certificate symbolically and on the rows. In Rust, see
`recurrence::{ImplicitRelation, CoefficientCondition, GeneratingDescription,
assess_generating_description, DiagonalWeight,
RecurrenceTransformation::transport_implicit, find_reconstruction_certificate,
ReconstructionCertificate, implicit_transform_recurrence,
verify_recurrence_implicit}`. The MCP tools are `transform_recurrence_implicit`,
`verify_recurrence_implicit` and `assess_generating_description`.

### Searching implicit relations directly

`polytool recurrence-implicit --search` fits implicit relations
`sum c(n,x) D^j Q[n+i] = 0` directly to exact rows, with derivatives allowed
on every row. The rows come from a defining rule (`source_recurrence`), a
bundled family (`oeis_id`) or `rows` (finite fits only), optionally after a
scalar registry map (`transformation`). Each entry of `supports`
(`role`, `max_shift`, `derivative_order`, `index_degree`, `variable_degree`)
is one ansatz box:

| `role` | Leader | Shift |
| --- | --- | --- |
| `successor` | no derivative on the row `n + R` | `max_shift >= 1` |
| `implicit_successor` | derivatives on every row | `max_shift >= 1` |
| `row_equation` | one row | `max_shift = 0` |

The columns and residual matrix are those of the complete-space search, and
every support gives one exact kernel with a rank certificate. Basis vectors
with independent leader projections (every basis vector of a row equation)
are the representatives, each in canonical form with its role, evidence,
generating description and assessment. A support whose kernel has no leader
is `absent_in_support` (`rows_only_obstruction`).

Completeness is certified when the source is a constant-index monic
pure-lag rule `chi = S^k - sum a_i(x) S^i` valid from its origin and the map
is the identity, `factorial_power` or `derivative`. Then every column
`x^v n^u D^s Q[n+j]` is annihilated by a power of `chi`: `h = max(u+s) + 1`
for the identity, `max(u+s) + m + 1` for `Q = D^m P` and `U + m(A+J) + 1`
for the factorial rescaling. Vanishing at `k h` consecutive bases then
implies vanishing for all bases, with derivatives on the leader and for row
equations alike, and the kernel is the complete space
(`certified_complete_space`).

Two further finite theorems cover other maps; `source.completeness_basis`
names the one used and `results[i].annihilator_power` its exponent `h`:

| Map | Residuals | Why |
| --- | --- | --- |
| `multiply_power`, `divide_power`, `affine`; `reversal`, `mobius`, `f_polynomial` with an unbounded schedule | `k (max(u+s) + 1)` | transported source: `Q_(n+k) = sum_i a~_i(x) Q_(n+i)` with `a~_i = a_i(phi) ell^(h(k-i))` index free in `Q(x)`, so the identity argument applies to `Q` itself |
| `diagonal` with `lambda = c (k!)^(-m) prod f_j(n,k)^(e_j)`, all `e_j > 0` (`binomial`, `product`) | `k (max(u+s) + w + 1)` for `m = 0`, `k (U + m(A+J) + w + 1)` for `m >= 1` | weight theorem: with `theta = x D`, `x^J` and `(k!)^m`-rescaling turn a column into `n^u x^s theta^(j falling) prod (theta+h)^m Lambda(n+i, theta) P_(n+i)`, `w` the total `(n,k)` degree of `Lambda = prod f_j^(e_j)` |

In the first theorem the lemma of the complete-space module holds verbatim
over the differential field `Q(x)`: a monic, index-free `chi~` commutes with
`S`, `x` and its coefficient derivatives. The theorems are proved in the
module documentation of `src/recurrence/implicit/search.rs`. Weights with negative exponents
(`inverse_binomial`), `(k!)^p` with `p > 0`, derivative antecedents, gamma
and bounded schedules are not covered and stay finite fits; a request whose
rows cannot reach the required count is `incomplete_budget`, never certified
on fewer residuals. Coefficient conditions use the truncation theorem with
`K' = k`, the index degree of `Lambda` added to `U_c` for weights; after a
Möbius-class map they are complete only when every `a~_i` is regular at
`x = 0`. Without a finite theorem, representatives are `finite_fit`,
unless an identity `M o T = sum H_nu S^nu L` with the defining rule `L`
(right division; `T` the identity, `D^m` or a polynomial diagonal weight)
proves them at every base where `L` holds (`certified_by_division`).

When a uniqueness argument fails only at certified exceptional indices of the
valuation grading, coefficient conditions `sum c(n) [x^j] Q[n+i] = 0` that pin
`[x^k]` of the new row are searched as well (`conditions`,
`condition_index_degree`). Under the same source the truncations
`P_n mod x^(K+1)` satisfy a constant linear recurrence of order `k (K+1)`,
so the condition space is complete from `k (K'+1)(U_c+1)` bases (`K' = k`, or
`k + m` after `D^m`). For A007318 divided by `k!` the search recovers
`D Q[n+1] - D Q[n] - Q[n] = 0` with the condition
`[x^0] Q[n+1] - [x^0] Q[n] = 0` (an equivalent description, both complete)
and the row equation `x D^2 Q[n] + (1 + x) D Q[n] - n Q[n] = 0`
(a consequence: a row equation does not determine later rows):

```sh
echo '{"oeis_id":"A007318",
  "transformation":{"kind":"diagonal","weight":{"kind":"factorial_power","power":1}},
  "supports":[{"role":"implicit_successor","max_shift":1,"derivative_order":1},
    {"role":"row_equation","derivative_order":2,"index_degree":1,"variable_degree":1}]}' \
  | polytool recurrence-implicit --search > search.json
polytool recurrence-implicit --verify-search < search.json
```

`--verify-search` reruns the search, requires equality, and re-checks every
kernel and condition certificate against matrices rebuilt from the reported
rows, every relation and condition on the rows, every division certificate
and every assessment. The report is `polytool.recurrence-implicit-search.v1`
(`source` with `constant_index_order` and `completeness_basis`, `rows`,
`results[i]` with `status`, `evidence`, `monomials`, `annihilator_power`,
`residual_count`, `certificate`, `representatives[j]` with `relation`,
`division_certificate`, `condition_searches`, `description`, `assessment`).
In Rust, see `recurrence::{search_implicit_relations,
verify_implicit_relation_search, ImplicitRelationSearchRequest,
ImplicitSearchSupport}`; the MCP tools are `search_implicit_relations` and
`verify_implicit_relation_search`.

### The lossy derivative map

The registry map `{"kind": "derivative", "order": m}` sends `P_n` to
`Q_n = D^m P_n` and forgets `[x^0..x^(m-1)] P_n`
(`RecurrenceTransformation::lost_data`). With these coefficients
`P_n = I^m Q_n + L_n` is rebuilt exactly
(`inverse_row_with_lost_data`); `inverse_row`, and therefore
`recurrence-transform`, rejects the map. A relation `L` of `P` is transported
by `transport_implicit` to `M = D^(m+r) L D^(-m)`, where
`r = max(0, max (deg_x c - d))` over the terms `c(n,x) D^d` of `L` makes
`D^(m+r) L` right divisible by `D^m`. The identity `M o D^m = D^(m+r) o L` is
recorded as a replayable certificate; `M` is a consequence of `L`. Conversely
a relation `M` of `Q` pulls back to `M o D^m` and a condition of `Q` to one of
`P` (`[x^k] Q = (k+m)!/k! [x^(k+m)] P`). The leader of `M o D^m` kills
`x^0..x^(m-1)`, so a description of `P` through `Q` is equivalent only
together with conditions pinning the lost coefficients
(`lossy_reconstruction`, `LossyReconstruction`). For the Eulerian rows
A008292 and `m = 1` the transported rule is
`D Q[n+1] + (-x + x^2) D^2 Q[n] + (-3 + 4 x - n x) D Q[n] + (2 - 2 n) Q[n] = 0`;
with `[x^0] P[n+1] - [x^0] P[n] = 0` and
`[x^1] P[n+1] - 2 [x^1] P[n] - n [x^0] P[n] = 0` the pulled-back description
of `P` is equivalent, and without these lost-data conditions it is a
certified consequence.
