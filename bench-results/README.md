# Benchmark Results

This directory is the conventional home for longer `polytool bench` outputs.
Keep short, named Markdown reports here when a run is useful for comparing
solver or search changes.

Recommended recurrence fixture command:

```sh
polytool bench recurrence-fixtures --repeat 3 \
  --summary --report bench-results/recurrence-fixtures/all.md
```

Recommended focused OEIS command:

```sh
polytool bench recurrence-fixtures --only oeis --repeat 3 \
  --summary --report bench-results/recurrence-fixtures/oeis.md
```

Machine-readable focused OEIS command:

```sh
polytool bench recurrence-fixtures --only oeis --repeat 3 --format json \
  > bench-results/recurrence-fixtures/oeis.json
```

Compare two JSON runs:

```sh
polytool bench compare old.json new.json --top 10
polytool bench compare old.json new.json --format json > comparison.json
```

The per-run TSV printed to stdout can still be redirected separately when raw
timings are needed. The Markdown report stores category summaries and fixture
summaries, which is usually enough for code review and handoff notes.
The JSON output additionally stores adaptive search diagnostics, so it is the
preferred format before and after search or solver changes. Those diagnostics
include both search rejection counters and cumulative timing buckets for the
main recurrence-finder stages, including modular lifting when a prime-field
candidate is exact-verified before the rational solver fallback.

Real-rootedness regression suite (verdicts are checked against exact expected
values; the baseline is a saved JSON run, not a checked-in file):

```sh
polytool bench real-rootedness --repeat 3 --format json \
  > bench-results/real-rootedness/full.json
polytool bench real-rootedness --baseline bench-results/real-rootedness/full.json
```
