# Searching equivalent representations

Part of the [CLI reference](cli.md); builds on
[recurrence search](cli-recurrence.md) and
[recurrence transformations](cli-recurrence-transforms.md).

### Searching equivalent representations

`polytool recurrence-search-representations` composes the coordinate changes,
transports, certified spaces and fits of the [recurrence search](cli-recurrence.md)
and [transformation](cli-recurrence-transforms.md) pages into one bounded queue search
over representations of a sequence. The source is a defining rule with seeds
and origin (`source_recurrence`), a bundled family (`oeis_id`), or rows only
(`rows`, `row_origin`; every result is then a finite fit). `transformations`
is the allowed set of registry maps. Representation 0 is the source; for each
representation below `limits.max_depth` every allowed map is applied in turn:

1. The composite map from the source is reduced to a normal form.
   Möbius-class maps (`multiply_power`, `divide_power`, `affine`, `reversal`,
   `mobius`, `f_polynomial`) compose to
   `Q_n(x) = K_0 K_1^n prod (x + s)^(u + v n) P_n((a x + b)/(c x + e))`,
   consecutive diagonal maps multiply their weights, and identity segments
   are removed. A map whose normal form is certified equal to an existing
   one is a `duplicate_map` edge. Rows are never compared: reversal twice at
   the same degree bound is the identity and deduplicates, while bounds `n`
   and `n + 1` give `x P_n`, a new representation. The palindromic rows of
   A007318 equal their reversal, but that map is a separate representation.
   It is marked `equal_sequence_as` only because both carry identical
   equivalent generating descriptions.
2. The rows are mapped exactly, and the source rows are rebuilt by inverting
   the whole chain (`reconstruction_verified`). The original source stays
   the reference.
3. Known relations are transported. Every certified scalar description of
   the parent (at most `max_sources_per_edge`, simplest first) becomes the
   defining rule of an embedded, replayable step. Möbius-class, gamma and
   symmetric-decomposition maps use a `polytool.recurrence-transform.v1`
   step. Diagonal maps and derivative antecedents use a
   `polytool.recurrence-implicit.v1` step, whose successor fit can be proved
   by reconstruction. A parent with certified implicit relations but no
   certified scalar rule feeds them to diagonal maps, derivative antecedents
   and the lossy `derivative` map as `implicit_relation` steps
   (`transport_implicit`), carrying its coefficient conditions where the map
   allows; the `derivative` map always takes this route. Without a certified
   parent relation the step only maps rows.
4. After a `factorial_power` diagonal map of a constant-index monic
   certified rule, complete spaces are certified on a bounded menu of
   supports. Lags run over `1..=max_rec_len` and derivative orders over
   `0..=max_diff_deg`. Uniform degrees increase up to the search bounds, the
   menu stops at the first success, and supports dominated by a success are
   skipped. `space_shapes` replaces the menu. After a Möbius-class map, a
   polynomial diagonal weight or `derivative`, the same menu runs the direct
   implicit search with one `successor` support per embedded
   `polytool.recurrence-implicit-search.v1` report (`implicit_spaces`;
   transported-source and weight theorems of
   [the transformation page](cli-recurrence-transforms.md#implicit-differential-relations-diagonal-maps-and-reconstruction)).
5. The new coordinates are fitted by the bounded menu (`search`, default
   lags 2, `x`-degree 2, index degree 1, derivative order 1, left factors in
   `n`, 2000 candidates).

Certified descriptions of a new representation feed its children.
Descriptions of the same relation are merged; the strongest label is kept
and the other provenances are listed in `also_found_by`. Each description has
one of three labels:

| `label` | Meaning |
| --- | --- |
| `finite_fit` | fitted to finitely many rows, or resting on a condition checked on rows only |
| `certified_annihilating_relation` | annihilation proved for every index of the domain, uniqueness not certified |
| `equivalent_generating_description` | certified annihilation, sufficient initial or boundary data, and a certified uniqueness argument (`assess_generating_description`) |

Certification is relative to the defining source, the maps actually
enumerated and the finite supports searched. Descriptions are ranked by shift
length, derivative order, coefficient degree, index degree, component count
and reconstruction cost (normal-form segments plus coefficient conditions),
then by label. Every limit that cuts work (`max_representations`,
`max_steps`, `max_space_shapes`, fit or elimination budgets) is listed in
`budget`, and the status is then `budget_exhausted`. The depth limit is a
requested bound: unexpanded representations at `max_depth` are listed in
`depth_frontier`. An optional `root_domain` request is attached to every
coordinate-change step, with the evidence labels of
[root-domain information](cli-recurrence-transforms.md#root-domain-information-for-transformed-results).
Discovery never depends on it.

Vector representations, such as `(a_n, b_n)` after `symmetric_decomposition`,
are expanded as well. Every scalar bijective map acts componentwise; a degree
schedule acts on `b_n` at `d_n - 1`, its symmetry index, so `gamma` with the
decomposition's schedule gives the gamma vectors of both components (their
symmetry is then certified by construction). Lossy and vector-producing maps
are rejected on vector states. The chain is inverted componentwise, and the
map normal form keeps one reduced componentwise suffix per component (`map`
segments `componentwise on a: ...`), so a reversal applied twice at the same
bounds is a `duplicate_map` on vectors too. A certified vector system of the
parent is transported exactly through Möbius-class maps (`steps[k].vector_transport`,
`transport_vector_componentwise`): it is `transported_identity` when the
parent system and every degree condition are certified, the degree bounds
`deg a_n <= d_n`, `deg b_n <= d_n - 1` being inherited from a certified
decomposition. States are fitted by a first-order vector system with a
shared left factor (`vector_fit`, finite evidence). Every vector system
carries `vector_assessment`: with its initial state it determines the states
when its left factor `q(n,x)` is nonzero at every index, certified by an
`x`-coefficient of `q` without integer roots in the domain; a certified system
starting at the row origin is then an `equivalent_generating_description`,
and two vector representations with identical such systems and initial
states are the same sequence (`equal_sequence_as`). For Eulerian A008292
split at `d_n = n`, the componentwise `f`-polynomials satisfy the certified
system `a(n+1) = (1 + 2t) a(n) + (t + t^2) a'(n) + (2t + 2t^2) b(n)`,
`b(n+1) = -a(n) + (t + t^2) b'(n)`.

A rank-two gamma transport (independent even and odd parts) is a
`transported_gamma_system` description: seeds and the two relations, with
component count 2 in the ranking. It is certified only where the symmetry of
the rows is certified, by its own step or by another step of the same edge
whose source rule is invariant under reversal; otherwise its evidence is
`checked_supplied_rows` and its label `finite_fit`.

A `derivative` map records the coefficients it forgets in the node's
`lost_data` (`{chain_position, order, rows}`); the chain is inverted with
them, and replay recomputes them from the source rows. Its
`implicit_relation` steps carry a `lossy` record: the parent's description
through `Q` with the conditions that pin the lost coefficients, and without
every condition that refers to them, each assessed. For A008292 the first is
an `equivalent_description` and the second a `certified_consequence`.

```sh
echo '{"oeis_id":"A007318","row_count":10,"limits":{"max_depth":1},
  "transformations":[{"kind":"diagonal","weight":{"kind":"factorial_power","power":1}}]}' \
  | polytool recurrence-search-representations > report.json
polytool recurrence-search-representations --verify < report.json
```

For A007318 divided coefficientwise by `k!`, representation 1 has three
equivalent descriptions. The first is the transported implicit rule
`Q'_(n+1) = Q'_n + Q_n` with its exceptional constant condition. The other
two are certified complete-space members: the first-order derivative rule
`n Q_n = (n + x) Q_(n-1) + x Q'_(n-1)` and the three-term rule
`n Q_n = (2n - 1 + x) Q_(n-1) - (n - 1) Q_(n-2)`. Without the constant
condition the implicit rule is only a `certified_annihilating_relation`.
For A065826 with `inverse_binomial {m: 1}` and `f_polynomial` (schedule
`d_n = n - 1` from index 1), the first map gives the Eulerian rows with a
reconstruction-certified rule. The second map transports that rule to
`Q_(n+1) = (x (1 + x) Q_n)'`. An excerpt of a report:

```json
{"schema": "polytool.recurrence-representation-search.v1", "status": "complete",
 "source": {"kind": "oeis", "row_origin": 0, "row_count": 10,
            "rule": {"schema": "polytool.recurrence.v1", "...": "..."}},
 "representations": [
  {"id": 0, "depth": 0, "chain": [], "map": {"identity": true, "segments": [],
   "index_domain": "n >= 0"}, "components": ["P"], "...": "..."},
  {"id": 1, "depth": 1, "parent": 0, "transformation": 0, "chain": [0],
   "map": {"identity": false, "index_domain": "n >= 0", "segments":
     ["Q_n(x) = sum_k lambda(n,k) [x^k]P_n x^k with lambda = (1/(k!)^1)"]},
   "components": ["Q"], "row_origin": 0,
   "states": [[["1"]], [["1", "1"]], [["1", "2", "1/2"]], "..."],
   "reconstruction_verified": true, "expanded": false,
   "descriptions": [
    {"label": "equivalent_generating_description", "kind": "certified_space",
     "evidence": "certified_complete_space",
     "text": "(n) P(n) = (n + t) P(n-1) + t P'(n-1)",
     "scalar": {"schema": "polytool.recurrence.v1", "first_index": 0,
                "initial_polynomials": [["1"]], "recurrence": {"...": "..."}},
     "description": {"first_index": 0, "seeds": [["1"]], "relations": ["..."],
                     "exceptional_index_conditions": [], "boundary_conditions": []},
     "assessment": {"status": "equivalent_description", "...": "..."},
     "cost": {"shift_length": 1, "derivative_order": 1, "coefficient_degree": 1,
              "index_degree": 1, "component_count": 1, "reconstruction_cost": 1},
     "provenance": "spaces[4].results[0].representatives[0]",
     "also_found_by": ["fits"]}, "..."],
   "steps": [{"source_description": 0, "implicit": {"schema":
     "polytool.recurrence-implicit.v1", "...": "..."}}],
   "spaces": [{"schema": "polytool.recurrence-spaces.v1", "...": "..."}],
   "fits": {"status": "complete", "configurations_considered": 42, "fits_found": 3}}],
 "edges": [{"parent": 0, "transformation": 0, "outcome": "new_representation",
            "representation": 1}],
 "ranked": [{"representation": 0, "description": 0,
             "label": "equivalent_generating_description",
             "text": "P(n) = (1 + t) P(n-1)", "cost": {"...": "..."}}, "..."],
 "steps_used": 10, "conventions": ["..."]}
```

The schema, for renderers:

| Path | Content |
| --- | --- |
| `schema`, `status` | `polytool.recurrence-representation-search.v1`; `complete` or `budget_exhausted` |
| `request` | the request, retained for replay |
| `source` | `kind` (`defining_rule`, `oeis`, `rows_only`), `rule` (with seeds), `row_origin`, `row_count` |
| `representations[i]` | `id`, `depth`, `parent`, `transformation` (index into `request.transformations`), `chain` (indices in application order), `map` (`identity`, `segments`, `index_domain`), `components`, `row_origin`, `states[row][component]` (ascending rational strings), `reconstruction_verified`, `map_conditions`, `descriptions`, `steps`, `spaces`, `fits`, `equal_sequence_as`, `expanded`, `notes` |
| `representations[i].descriptions[j]` | `label`, `kind` (`defining_rule`, `transported_rule`, `transported_vector_system`, `transported_gamma_system`, `implicit_transport`, `reconstruction_certified_fit`, `certified_space`, `finite_fit`, `vector_fit`), `evidence`, `text`, optional `scalar` (`polytool.recurrence.v1` with seeds, original indices), `vector`, `description` (generating description), `assessment`, `cost`, `provenance`, `also_found_by`, `note`, `vector_assessment` (`status`, `denominator_coefficient`, `statement`) |
| `representations[i].steps[k]` | `source_description` (parent description used as the defining rule, absent for rows only), one of `transform` / `implicit` (embedded reports) / `implicit_relation` (`source_relation`, `relation`, `clearing`, `certificate`, `dropped_conditions`, `description`, `assessment`, `lossy`, `holds_on_rows`, `checked_bases`) / `vector_transport` (`component_maps`, `transported`, `text`, `evidence`, `holds_on_rows`, `checked_transitions`), or `error` |
| `representations[i].implicit_spaces[k]` | embedded `polytool.recurrence-implicit-search.v1` reports (one `successor` support each) after Möbius-class maps, polynomial diagonal weights and `derivative` |
| `representations[i].lost_data[k]` | for each lossy map of the chain: `chain_position`, `order`, `rows` (the forgotten coefficients of the rows entering the map) |
| `edges[k]` | `parent`, `transformation`, `outcome` (`new_representation`, `duplicate_map`, `rejected`, `skipped_budget`), `representation`, `detail` |
| `ranked[k]` | `representation`, `description`, `label`, `text`, `cost` |
| `depth_frontier`, `budget`, `steps_used`, `conventions` | unexpanded representations at `max_depth`; budget events `{limit, representation, detail}`; backend calls; conventions |

Successor-rule `text` uses the `polytool.recurrence.v1` display (`P(n)`,
variable `t`, target index), and implicit relations use `Q[n+i]` with base
index `n`. A renderer should distinguish `finite_fit` from the two certified
labels and show `evidence` and `assessment.uniqueness` for certified results.
`--verify` reruns the search from the embedded request and requires equality.
It then replays every embedded step and space certificate with its own
verifier and checks that each step's defining rule is a certified
description of its parent. Every representation is re-derived from the
source rows (map and exact inverse), and every description is re-checked on
its rows. In Rust, see `recurrence::{RepresentationSearchRequest,
RepresentationSearchLimits, search_representations, RepresentationSearchReport,
verify_representation_search, representation_label}`; the MCP tools are
`search_recurrence_representations` and
`verify_recurrence_representation_search`. Not yet covered: implicit
transport through Möbius-class maps, transport of vector systems through
gamma and diagonal maps (their rows are mapped and fitted only), and complete
spaces after gamma maps, derivative antecedents and inverse polynomial
weights.
