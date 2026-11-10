# MCP server

The repository also contains a local Model Context Protocol server in
`mcp/`. It exposes the exact polynomial routines to MCP clients over stdio.

Build, install and client-configuration instructions (including
`polytool-mcp --help` and `--version`), the tools-only transport and the output
conventions are documented in [`../mcp/README.md`](../mcp/README.md) (sections
"Build", "Install", "MCP client configuration" and "Output format").

Available MCP tools:

- `parse_polynomials`
- `check_polynomial_family`
- `polynomial_properties`
- `check_interlacing_pair`
- `check_interlacing_sequence`
- `check_interlacing_profile`
- `real_roots`
- `find_recurrence`
- `generate_recurrence_rows`
- `resultant`
- `discriminant`
- `ehrhart_hstar`
- `analyze_decomposition`
- `generate_sequence`
- `check_tnn_network`
- `recognize_brenti_strip_digraph`
- `recognize_pf_network`
- `check_left_production_matrix`
- `transform_recurrence`
- `verify_recurrence_transform`
- `transform_recurrence_implicit`
- `verify_recurrence_implicit`
- `assess_generating_description`
- `search_recurrence_representations`
- `verify_recurrence_representation_search`
- `search_implicit_relations`
- `verify_implicit_relation_search`
- `enumerate_real_rooted`
- `find_veronese_section_lifts`
- `certify_veronese_section_lift`

`enumerate_real_rooted` takes `{"degree", "leading_coefficients_descending",
"max_results", "max_nodes", "max_coefficient_bits"}` (the options of
`polytool enumerate-real-rooted`; the prefix `a_0, a_1, ...` is descending,
`leading_desc` is accepted as an alias) and returns the fields of its JSON
output with ascending exact integer strings.  Unlike the CLI, the MCP server
bounds the work: degree at most 32, at most 33 prefix coefficients of at most
256 bits, `max_nodes` default 5000 and at most 25000, `max_results` default 500
and at most 2000, `max_coefficient_bits` default 1024 and at most 2048.  A
larger value is rejected and a limit cannot be lifted.  A budget that stops the
search yields `status: "resource_limit"` with `limit_reached`, never an empty
"complete" answer.

`check_tnn_network` and `recognize_brenti_strip_digraph` return the same
fields as the JSON output of `polytool tnn-network` and
`polytool brenti-recognize` (without the `schema` tag), including the fixed
guardrail `note`: TNN of the coefficient matrix, and a strip-digraph
certificate, say nothing about real-rootedness of the row polynomials.
`recognize_pf_network` and `check_left_production_matrix` mirror
`polytool pf-network` and `polytool left-production`; their certificates do
prove PF (real-rooted) rows, and their notes say what they do not prove (at
most 64 rows of degree at most 64 and steps of length and rise at most 8 for
`recognize_pf_network`; at most 32 rows for `check_left_production_matrix`).

`find_veronese_section_lifts` takes `{"recurrence_json", "max_modulus",
"parent_rows", "max_candidates"}` (the options of
`polytool recurrence-veronese-lifts`) and returns its
`polytool.recurrence-veronese-lifts.v1` report: parents `P_n(x)` whose
coefficient sections are the supplied sequence, each with its
`parent_recurrence` record (usable with `generate_recurrence_rows`) and a
status, of which `proved_product_parent` and `proved_parent_route` are
proofs (`tail_first_index` names the first used row when supplied rows that
do not follow the recurrence were dropped).  The MCP server
bounds the work: recurrence order at most 12, `max_modulus` at most 6,
`parent_rows` default 12 and at most 64, `max_candidates` default 256 and at
most 4096.

`certify_veronese_section_lift` takes `{"recurrence_json", "candidate",
"checked_rows", "attestations", "try_root_domain"}` and returns the
`polytool.recurrence-veronese-certificate.v1` report of
`recurrence::certify_veronese_section_lift` for one candidate lift
`Q_n = lambda(n) V_(r,j)(P_(n+s))`: the identity status (`certified` with a
replayable forward-division certificate, `supported_by_rows` or `refuted`),
the parent proof and the verdict, of which only `proved` is a proof.  The MCP
server bounds the work: input and parent order at most 12, modulus at most 4,
`|index_shift|` at most 4, derivative order at most 6, `checked_rows` default
16 and at most 64, at most 8 attestations, and the library's budget of 10000
operator terms.

Every MCP tool that takes polynomial coefficients accepts arbitrary-size exact
integers: `parse_polynomials`, `check_polynomial_family`,
`polynomial_properties`, `coefficient_tests`, the interlacing tools,
`real_roots`, `analyze_decomposition`, `resultant`, `discriminant`,
`ehrhart_hstar`, `hstar_inequalities`, the cyclic-sieving tools, the network
tools, `find_recurrence` and `find_recurrence_representations` (including
`polynomials`, `expressions` and `text` rows), and sequence generation.  A
coefficient may be a JSON integer (read exactly within the 64-bit range, up to
`u64::MAX`) or a quoted decimal string, which is the form to use beyond that
range and beyond the client's safe-integer range; a floating-point JSON number
such as `1e30` is rejected.  Recurrence
`coefficients` rows may also be exact rational strings such as `"3/7"`.
Responses return integers as decimal strings, including the coefficients of
`parse_polynomials`, `real_roots`, `analyze_decomposition` and
`check_polynomial_family` (item coefficients, gamma coefficients and the Lace
matrix).  The remaining integer-typed arguments are indices and sizes, not
coefficients.

`find_recurrence` exposes all adaptive-search controls as top-level MCP
arguments. For example, a compact exact search can use:

```json
{
  "coefficients": [[1], [2], [4], [8], [16]],
  "max_rec_len": 1,
  "max_var_deg": 0,
  "max_idx_deg": 0,
  "max_diff_deg": 0,
  "include_code": false
}
```

The legacy nested `options` object remains accepted, with top-level values
taking precedence. Omitting `include_code` preserves the full generated-code
response.

See [`mcp/README.md`](../mcp/README.md) for request schemas, examples, and
development notes.

## General recurrence proof discovery

`discover_recurrence_proofs` takes `recurrence_json` and optional
`max_start_delay` (default 3, at most 16), `max_representation_candidates`
(default 300, at most 2000, 0 disables the alternative-representation stage)
`max_lift_candidates` (default 6, at most 100, 0 disables the catalog-lift
stage) and `max_affine_coordinates` (default 8, at most 64, 0 disables the
affine-coordinate stage).
It tries the eight general theorem families of the parent cascade (product,
three-term, Liu--Wang, derivative lag, Ma--Wang interval, Liu--Wang
three-term, iterated preserver, hypergeometric row) without catalogue
lookup, then the tail after rows that are zero or not real-rooted, then
alternative rules from the recurrence finder, each used only after an
all-index equality certificate, and last the catalog lifts: a normalized
identification with a bundled OEIS catalog family or a backward basis-map
lift, whose parent polytool proves by the stages above (`catalog_lifts`;
a parent that is only flagged `lean_proved` gives the status
`conditional_on_attestation`, never `proved`), then the affine coordinates:
the criteria on the record in exact coordinates `x = a u + b` derived from
the rule, a proof mapped back to `x` (`affine_coordinates`).
`verify_recurrence_proofs` takes
`report_json` and replays the recorded routes (every attempt, the tail and
every equality certificate; the attempts must be of exactly the criteria
the current search records, so reports saved before a criterion was added
are recomputed). A proved
result concerns the recurrence-defined sequence; finite fits require a
separate all-index identity. Lean replay is not attempted. See
[the proof-discovery guide](proof-discovery.md).

## Renamed request fields

After 0.6.0, polytool renamed request fields so that a tool field, the matching CLI flag
and the JSON key use the same spelled-out words (`max_index_degree` is
`--max-index-degree`). The old names remain accepted as aliases, so existing
clients keep working; the published tool schemas show the new names only.

| Old field | New field | Tools |
|-----------|-----------|-------|
| `min_rec_len`, `max_rec_len` | `min_recurrence_length`, `max_recurrence_length` | find_recurrence, find_recurrence_representations, recurrence options |
| `min_var_deg`, `max_var_deg` | `min_variable_degree`, `max_variable_degree` | same |
| `min_idx_deg`, `max_idx_deg` | `min_index_degree`, `max_index_degree` | same |
| `min_diff_deg`, `max_diff_deg` | `min_derivative_order`, `max_derivative_order` | same |
| `min/max_inhomo_var_deg`, `min/max_inhomo_idx_deg` | `min/max_inhomogeneous_variable_degree`, `min/max_inhomogeneous_index_degree` | same |
| `max_denom_var_deg`, `max_denom_idx_deg` | `max_denominator_variable_degree`, `max_denominator_index_degree` | same |
| `max_laurent_x_power`, `min_margin`, `no_verify` | `max_laurent_power`, `min_equation_margin`, `no_held_out_rows` | same |
| `first_row_index` | `first_index` | recognize_brenti_strip_digraph, recognize_pf_network, check_left_production_matrix |
| `rows`, `first_row` | `count`, `from_index` | generate_oeis_rows |
| `rows`, `additional` | `count`, `additional_count` | generate_recurrence_rows |
| `maximum_start_delay` | `max_start_delay` | discover_recurrence_proofs |
| `offsets` | `order_offsets` | cyclic_sieving_sequence |

Responses are additive: `first_index` is reported next to `first_row`
(OEIS tools) and next to `first_row_index` (planar-network tools).

