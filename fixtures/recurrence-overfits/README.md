# Recurrence Overfit Fixtures

These are negative recurrence-search fixtures from curated OEIS
real-rootedness research datasets.  They are finite row prefixes that produced
implausible high-coefficient recurrences under broad or unverified recurrence
searches.

Unlike `fixtures/recurrence-benchmarks/`, these are not expected-success
fixtures.  They are intended for testing recurrence-confidence heuristics,
held-out verification, and search-space changes that should reduce spurious
fits.

## Fixtures

| id | rows | recurrence bytes | note |
|---|---:|---:|---|
| `A177970` | 98 | 2537 | High-rational-coefficient derivative recurrence; compare with the compact true fixture `30_oeis_a177970_order8_closed_form` in `recurrence-benchmarks`. |
| `A181000` | 19 | 16256 | Extreme 144-unknown no-verification fit with very large rational coefficients. |

## Manual Stress Commands

The A181000 bad fit can be reproduced deliberately with:

```sh
polytool recurrence \
  --min-rec-len 4 --max-rec-len 4 \
  --min-var-deg 3 --max-var-deg 3 \
  --min-idx-deg 2 --max-idx-deg 2 \
  --min-diff-deg 2 --max-diff-deg 2 \
  --no-verify \
  < fixtures/recurrence-overfits/rows/A181000.txt
```

The important diagnostic is that the result uses all 19 rows for fitting and
has zero held-out verification rows.  Treat such output as a stress artifact,
not as recurrence evidence.

The stored recurrence text files and row files were derived from curated
research datasets of OEIS recurrences.
