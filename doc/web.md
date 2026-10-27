# Web arbitrary-precision contract and build

The recurrence panel also offers **Compare representations**. It uses the
visible length, coefficient-degree and derivative-order bounds, and a global
configuration limit (default 2000). A table lets the user select a formula,
copy its text/LaTeX/JSON, and extend the sequence with that selected formula.
It displays checked starting indices, held-out row counts, and partial coverage.
The interface preserves the task-oriented UI and row-extension feature deployed
from `878113b`, which had not yet reached the parent integration branch.

The **Full search** mode uses the entered polynomial rows by default. Choose
the maximum number of previous terms and derivative order; coefficient bounds
are in a separate expandable section. The first-row index is explicit. The
browser searches every length/order combination within these bounds, with up
to 16 combinations in the guided form. Examples use defining recurrences for
A171608, Hermite polynomials and automatic polynomial-degree selection.

Results lead with typeset formulas and a copy button, with additional formulas
and search bounds expandable. Short labels distinguish finite-row matches,
identities following from a defining recurrence, bounded absence and unfinished
searches. Certificate exports and linear-algebra diagnostics are kept out of
the webpage; the CLI and library retain their complete reports and replay APIs.

**Advanced: request JSON** supports custom defining recurrences, sparse
supports and resource limits. Editing the JSON switches to custom mode and
disables the guided form; the main input is then unused. Any relevant input or
settings change cancels active work and clears stale results. The existing
WASM export `parse_and_format(input)` supplies exact decimal-string coefficient
rows without computing root properties. The form rejects the entire input if
any row is invalid. Rational rows remain
available through the custom JSON request.
The WASM exports `find_recurrence_spaces(request_json)` and
`verify_recurrence_spaces(report_json)` share the Rust backend.
For the sources `A171608` and the binomial example, the guided form can
rescale coefficients first (`coefficient_rescaling`, divide the coefficient of
`x^k` by `k!` or `(k!)^2`). Formulas then name the rescaled rows `Q`.

The **Transform** tab (polynomial-sequence workflow) applies the coordinate
changes of `recurrence-transform` and the diagonal and antecedent maps of
`recurrence-implicit` in the browser. The source is a bundled OEIS family
(searchable by A-number), a `polytool.recurrence.v1` record with initial rows
(for example the one found in the Recurrence tab), or the rows entered above,
which have no defining recurrence. The form builds the same JSON request as the
CLI. **Advanced: request JSON** shows the request and accepts a custom one.
Root location (`root_domain`) and an independent fit (`search`) are optional.
The root-location checkbox applies to coordinate changes and to the diagonal
and derivative-antecedent maps alike: the request carries a top-level
`root_domain`, which `transform_recurrence_implicit` passes to the library
(`parse_implicit_transform_input` and
`implicit_transform_recurrence_with_root_domain`), so implicit reports carry a
`root_domain` section that replays with the report. A candidate interval end
is a rational or a linear function of the index (`-1`, `-1-n`, `n/2`), shown
as written by the report (for example `[-1-n, 0]`). Gamma coordinates are
transported when the symmetric and antisymmetric parts of the rule are
proportional (Eulerian, Narayana, type B Eulerian), with the rule and the rows
written in `z = x/(1 + x)^2`; other rules are refused with the library's
explanation. Examples cover the Eulerian `f`-polynomial, gamma transport,
an index-dependent root window and symmetric decomposition, the Narayana
reversal, Pascal rows divided by `k!` (with and without root location), and
A065826 divided by `k + 1`.

Results are grouped by report section. Each claim carries one evidence label,
with its own colour, border and symbol:

| Label | Codes | Meaning |
| --- | --- | --- |
| Certified (✔) | `certified_all_indices`, `completed_recurrence_induction`, `certified_family_invariant`, `certified_by_reconstruction`, `equivalent_description` | proved for every index |
| Transported (⇒) | `transported_identity`, `conditional_on_source`, `certified_consequence` | exact consequence of the defining recurrence |
| Theorem (§) | `conditional_theorem_match`, `paper_proof_not_lean_checked` | hypotheses checked exactly; the remaining obligations are listed |
| Checked (◇) | `checked_supplied_rows`, `finite_exact_observation`, `fixed_cap` | exact on the listed rows only |
| Finite fit (≈) | `finite_fit`, `exact_fit`, `finite_evidence_only`, `finite_evidence` | fitted to finitely many rows |
| Assumption (?) | `supplied_assumption`, `proposed_interval` | supplied or proposed, not proved |

The Jensen test of an implicit root-domain section (`fixed_cap`, labelled
"Fixed degree cap only") sits in its own dotted block under a warning that it is
a finite check, never a proof for larger degrees; certified claims are listed
before it and keep their solid style. An unknown label is shown verbatim in a
neutral style. A report schema without
a dedicated view is rendered by listing every object that carries an
`evidence` field. To add a view, add an entry to `REPORT_RENDERERS` in
`index.html`. Each report card can check the source, transported or fitted
scalar recurrence against the interlacing theorems
(`polytool.recurrence-theorems.v1`). It can also compute exact root enclosures
of the source or transformed rows (`polytool.root-enclosures.v1`), copy or
download the report JSON, and replay it. **Import a report and replay it**
accepts any of these schemas or a `polytool.recurrence-spaces.v1` report.
Replay recomputes the report from its embedded request. A report that fails
replay is marked as failed, and its content is shown only on request. Replay
does not strengthen any label.

WASM exports (JSON text in, JSON text out, `{"error": ...}` on failure):
`transform_recurrence(request)` and `verify_recurrence_transform(report)`;
`transform_recurrence_implicit(request)` and
`verify_recurrence_implicit(report)` (the implicit request may carry
`root_domain`, and a report with that section is replayed with it);
`root_enclosures(request)` and
`verify_root_enclosures(report)` (recompute and compare);
`recurrence_theorems(recurrence_json)` and `verify_recurrence_theorems(report)`
(supplied-recurrence reports only); `verify_report(report)`, which dispatches
on `schema`; `recurrence_latex(recurrence_json)`; and `oeis_catalog_index()`.
The reports are the library schemas unchanged, so browser reports replay with
`polytool recurrence-transform --verify` and the MCP verifiers, and CLI
reports replay in the browser. Requests are limited to 1 MiB. Transformation
limits above the default `TransformLimits` service ceilings are rejected, as
in the MCP server. The guided form generates at most 64 rows, and enclosure
widths must be at least `2^-256`. Computations run in a Web Worker that
**Cancel** terminates.

`search_recurrence_representations(request)` and
`verify_recurrence_representation_search(report)` run and replay the bounded
representation search (`polytool.recurrence-representation-search.v1`);
`verify_report` dispatches that schema as well. Limits above
`RepresentationSearchLimits::service_ceiling()` are rejected, as in the MCP
server, and a report larger than 1 MiB is refused at the source because its
replay would be. In the Transform tab, the mode **Search representations**
offers a checklist of allowed maps (shared degree-bound and diagonal-order
fields), depth 1 to 3, at most 24 representations and at most 24 rows. The
dedicated view shows the status and every budget event (an incomplete search
is labelled as such), the tree of maps with duplicate, rejected and skipped
applications and `equal_sequence_as` marked by symbol and border as well as
colour, the simplest descriptions overall, and for each representation its map,
reconstruction check, ranked descriptions with their label, evidence and
uniqueness argument, and its rows. The labels `equivalent_generating_description`
and `certified_complete_space` use the certified style,
`certified_annihilating_relation` and `defining` the transported style, and
`finite_fit` the fit style, so the six-class scheme is unchanged. Small searches
(at most six representations) open every representation; larger ones open from
the tree links.

The browser wrapper parses polynomial input with `parse_polynomials_bigint`.
Every integer coefficient returned across its string-based WASM/JSON boundary
is a canonical decimal string, including input, gamma, decomposition, and
recurrence-initial coefficients. Resultants and discriminants are decimal
strings as well. Browser code must keep these values as strings or convert
them to JavaScript `BigInt` for exact arithmetic; converting them to `Number`
can round values above `2^53`.

Build the no-modules bundle with `wasm-pack`:

```sh
wasm-pack build web --target no-modules --release --out-dir pkg
```

`web/build.sh` runs the same build from inside `web/`. The page `web/index.html`
loads `pkg/polytool_web.js` and `pkg/polytool_web_bg.wasm`; serve the `web/`
directory (for example with `python3 -m http.server`) to try it locally.
