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

The MCP property, interlacing, resultant, discriminant, Ehrhart conversion,
sequence-generation, and recurrence tools accept arbitrary-size integer
coefficients. JSON integers may be provided directly when they fit the client
stack, and larger exact integers should be sent as strings.

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
