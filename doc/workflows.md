# CLI workflows

## Check a polynomial or a whole family

```sh
printf '1,11,11,1\n' | polytool real-rooted --json
polytool sequence eulerian 12 --format rows | polytool family-check --json
```

Coefficients are in ascending degree order. Use `real-rooted` when only
the root test is needed; `properties` computes several additional properties.
`family-check` also checks consecutive interlacing. These are exact tests of
the supplied rows, not proofs for an infinite family.

## Find a recurrence and extend the sequence

```sh
polytool sequence eulerian 12 --format rows > rows.txt
polytool recurrence --json --first-index 1 < rows.txt > recurrence.json
polytool recurrence-generate --recurrence recurrence.json --rows 20 --format rows
```

Eulerian and Narayana output starts at index 1; the other standard families
start at 0. Row files do not carry indices: set `--first-index` explicitly
when needed. `recurrence` defaults to 1. The generator reads the index stored
in the recurrence. `--rows` counts total rows, including the initial rows.
Recurrence discovery from finitely many rows gives a candidate identity;
held-out checks are evidence, not an all-indices proof.

`sequence`, `oeis generate` and `recurrence-generate` all accept
`--format rows`. Human-readable labelled output remains the default for the
first two. OEIS `--format csv` is a long `n,k,value` table, not coefficient
rows suitable for recurrence input.

## Reuse theorem reports and choose a formula

```sh
polytool recurrence --json --theorem-matches < rows.txt > report.json
polytool recurrence-generate --recurrence report.json --rows 20

polytool recurrence --representations --json < rows.txt > menu.json
polytool recurrence-generate --recurrence menu.json --candidate 0 --rows 20
polytool recurrence-theorems --recurrence menu.json --candidate 0 --json
polytool recurrence-veronese-lifts --recurrence menu.json --candidate 0 --json
```

All three readers accept `polytool.recurrence.v1`, the embedded `recurrence_json`
in `polytool.recurrence-theorems.v1`, and entries of the `representations`
array in `polytool.recurrence-representations.v1`. Candidate indices are
zero-based positions in that array, not positions in its `frontier` array.
A menu with multiple entries requires `--candidate`; one entry is selected
automatically. Other report schemas need their own replay commands. Extracting
a recurrence does not verify a report or strengthen its evidence labels.
Use `--recurrence -` to read any of these supported inputs from stdin.

## Replay exact root enclosures

```sh
printf '{"polynomials":[["-2","0","1"]],"max_width":"1/100"}\n' |
  polytool root-enclosures --json > enclosures.json
polytool root-enclosures --verify < enclosures.json
```

Replay recomputes the embedded request and checks the entire report, including
root multiplicities and interval verdicts. It verifies the supplied finite
computation; it does not prove real-rootedness for an infinite family. The
browser verifier and CLI use the same library replay function.

`root-enclosures`, `recurrence-transform`, `recurrence-implicit` and
`recurrence-search-representations` always produce JSON. They accept `--json`
for consistency with other commands and `--compact` for one-line output.

## Get formulas and exact scalar results

```sh
polytool recurrence --latex < rows.txt
printf '2,-3,1\n' | polytool discriminant --json
printf '%s\n' '-1,1' '-3,1' | polytool resultant --json
```

LaTeX includes initial conditions and the starting index. General coefficient
factorization is a separate feature. Scalar JSON stores exact integers as
decimal strings: discriminants use `{"items":[{"index":0,"discriminant":"1"}]}`;
resultants use `{"resultant":"-2"}` and require exactly two input polynomials.

## Help, completion and scripting

```sh
polytool --help
polytool help recurrence
polytool oeis add --help
polytool help --json
polytool completions bash > polytool-completion.bash
source polytool-completion.bash
```

Help groups commands by task. `help --json` emits the `polytool.commands.v1`
index, with names, aliases, groups and summaries. Completion scripts cover
command names for Bash, Zsh and Fish (plus OEIS subcommands in Bash);
they do not complete every option value.

Diagnostics go to stderr; data goes to stdout. Invalid options and malformed
input fail with a nonzero status (normally 2). Some commands use 1 for a
negative verification result, and bounded searches use 3 for unfinished work.
An exact root test reporting `false` is still a successful computation.
A completed recurrence search finding nothing also remains successful, with
`status: "not_found"` in JSON. Inspect result fields as well as exit status.
For shell pipelines, enable `set -o pipefail` in Bash to catch upstream failures.
