# Family sidecars

One strict JSON record per curated family. The authoritative typed schema and
validators are in `src/families/`; see [the guide](../../doc/cli-families.md).
Uncurated OEIS records are derived directly and need no sidecar. OEIS generators,
status flags, audit strings and dates always come from `oeis-catalog.jsonl`.
Run `polytool families validate` before building the site or exporting drafts.
