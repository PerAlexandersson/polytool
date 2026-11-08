//! Specification of the `oeis` command and its subcommands.

use crate::cli::spec::{
    option, positional, switch, CommandSpec, InputKind, OutputKind, COMMAND, EXIT_ANSWERED,
};

const ID_HELP: &str = "OEIS A-number (A followed by six digits)";

const LIST: CommandSpec = CommandSpec {
    name: "list",
    summary: "List the bundled recurrence-backed OEIS families",
    description: r#"List the catalog of OEIS polynomial families bundled with polytool. The text
output has one tab-separated line per family: A-number, status, whether a
verified b-file prefix exists (bfile or no-bfile) and the sequence name. Entries
with the status experimental are hidden unless --include-experimental is given."#,
    input: InputKind::None,
    output: OutputKind::JsonSwitch,
    flags: &[switch(
        "include-experimental",
        "Also list experimental catalog entries",
    )],
    details: r#"JSON output is polytool.oeis-catalog.v1: {"count": ..., "sequences": [...]}."#,
    examples: &[
        "polytool oeis list",
        "polytool oeis list --json --include-experimental",
    ],
    exit_codes: EXIT_ANSWERED,
    see_also: &["oeis info", "oeis generate"],
    ..COMMAND
};

const INFO: CommandSpec = CommandSpec {
    name: "info",
    summary: "Show one bundled OEIS family and its recurrence",
    description: r#"Show the catalog record of one bundled family: name, status, row layout and
indexing, evidence, and the exact recurrence that generates the rows. With --json
the report also carries the recurrence as polytool.recurrence.v1 data (which
`polytool recurrence-generate` reads), its LaTeX form and a Mathematica form."#,
    input: InputKind::None,
    output: OutputKind::JsonSwitch,
    positionals: &[positional("A-NUMBER", ID_HELP)],
    examples: &[
        "polytool oeis info A008292",
        "polytool oeis info A008292 --json",
    ],
    exit_codes: &[
        (0, "found"),
        (
            2,
            "usage error, or the A-number is not in the bundled catalog",
        ),
    ],
    see_also: &["oeis list", "oeis generate", "topic:recurrence-json"],
    ..COMMAND
};

const GENERATE: CommandSpec = CommandSpec {
    name: "generate",
    summary: "Generate rows of a bundled OEIS family",
    usage: &["oeis generate <A-NUMBER> --count <N> [OPTIONS]"],
    description: r#"Generate OEIS polynomial rows from the bundled exact recurrence of one family.
--count gives the number of complete rows and is required. --from-index selects
the first displayed OEIS row (default: the catalog's start row); strict b-file
output must begin at the catalog's start row.

Formats (several have aliases):
  polynomial   expanded polynomials, one per row (default; poly, polynomials)
  rows         coefficient lists, ascending degree (coefficients, text)
  triangle     aligned triangle (table)
  json         rows as JSON
  jsonl        one JSON row per line (json-lines)
  csv          comma-separated coefficients
  bfile        OEIS b-file lines `index value` (b-file); --max-terms caps the
               number of terms and never splits a row

An experimental catalog entry needs --include-experimental."#,
    input: InputKind::None,
    output: OutputKind::Formats {
        values: &[
            "polynomial",
            "polynomials",
            "poly",
            "rows",
            "coefficients",
            "text",
            "triangle",
            "table",
            "json",
            "jsonl",
            "json-lines",
            "csv",
            "bfile",
            "b-file",
        ],
        default: "polynomial",
    },
    flags: &[
        option("count", "N", "Number of complete rows to emit (required)").deprecated(&["rows"]),
        option(
            "from-index",
            "K",
            "First displayed OEIS row (default: catalog start)",
        )
        .deprecated(&["start-row"]),
        option("max-terms", "N", "B-file cap on terms; never splits a row"),
        switch(
            "include-experimental",
            "Permit an experimental catalog entry",
        ),
    ],
    positionals: &[positional("A-NUMBER", ID_HELP)],
    examples: &[
        "polytool oeis generate A008292 --count 6 --format triangle",
        "polytool oeis generate A008292 --count 30 --format bfile",
        "polytool oeis generate A008292 --count 6 --format rows | polytool real-rooted",
    ],
    exit_codes: &[
        (0, "rows written"),
        (
            2,
            "usage error, unknown or experimental entry, or --max-terms smaller than the first row",
        ),
    ],
    see_also: &["oeis info", "recurrence-generate", "topic:output-formats"],
    ..COMMAND
};

const ADD: CommandSpec = CommandSpec {
    name: "add",
    summary: "Search a recurrence for OEIS rows and record it in the catalog",
    usage: &[
        "oeis add <A-NUMBER>... [OPTIONS]",
        "oeis add --ids-file <FILE> [OPTIONS]",
        "oeis add <A-NUMBER> --bfile <FILE> --name <TITLE> --layout <LAYOUT> [OPTIONS]",
        "oeis add <A-NUMBER> --recurrence-json <FILE> [--prefix-rows <P>] [OPTIONS]",
    ],
    description: r#"Read OEIS terms, cut them into complete rows, run the adaptive recurrence search
(the same search and options as `polytool recurrence`) on a prefix of the rows,
and check the recurrence exactly against every remaining term through the
catalog's own row generator. A passing entry is written to the catalog data file
with provenance `polytool_recurrence_search`; its status is `verified` when
held-out rows were checked, else `validated`. Rebuild polytool to bundle new
entries. This is a maintenance command for a source checkout.

Rows: OEIS `tabl` sequences are read as regular triangles (rows of length
1, 2, 3, ...). Other sequences need an explicit rule: --layout triangle (also
regular_triangle, tabl) treats the terms as a regular triangle; --layout table
(also tabf) is an irregular table and requires --row-widths. --fit-rows gives
the number of rows to the search (default: all but max(2, rows/4), which are
held out).

Supplied recurrence (no search): --recurrence-json reads a polytool.recurrence.v1
object (rule, first_index and initial polynomials; one A-number only). Its first
initial row is matched with the first OEIS row after the prefix, and the rows it
generates are compared exactly with every OEIS term. Provenance is
`supplied_recurrence`, status `validated`. A table layout may omit --row-widths:
the generated rows give the widths. Search options and --fit-rows do not apply.

Outcomes per A-number: added, updated, already_present, no_recurrence (nothing
within the bounds or budget), layout_unsupported, insufficient_data, mismatch (a
held-out term differs) and unavailable (no record)."#,
    input: InputKind::None,
    output: OutputKind::JsonSwitch,
    flags: &[
        option(
            "oeis-data",
            "DIR",
            "Local OEIS internal-format copy (default: $OEIS_DATA_DIR)",
        )
        .heading("Input"),
        option("bfile", "FILE", "Explicit b-file (one A-number only)").heading("Input"),
        option("name", "TITLE", "Sequence title (with --bfile)").heading("Input"),
        option(
            "first-term-index",
            "K",
            "Index of the first term and displayed row (with --bfile)",
        )
        .heading("Input")
        .deprecated(&["offset"]),
        option(
            "ids-file",
            "FILE",
            "A-numbers, one per line (# comments allowed)",
        )
        .heading("Input"),
        option("layout", "LAYOUT", "How terms form rows: triangle or table").heading("Rows"),
        option(
            "row-widths",
            "W+S[/R]",
            "Row i has W + S*floor(i/R) terms (R defaults to 1)",
        )
        .heading("Rows"),
        option(
            "fit-rows",
            "N",
            "Rows given to the search (rest are held out)",
        )
        .heading("Rows"),
        option(
            "recurrence-json",
            "FILE",
            "Use a supplied polytool.recurrence.v1 record (no search)",
        )
        .heading("Supplied recurrence"),
        option(
            "prefix-rows",
            "P",
            "Store the first P OEIS rows as authoritative prefix rows",
        )
        .default("0")
        .heading("Supplied recurrence"),
        option(
            "catalog",
            "FILE",
            "Catalog data file (default: ./data/oeis-catalog.jsonl if present)",
        )
        .heading("Catalog"),
        switch(
            "dry-run",
            "Do not write; print each new entry as a JSONL line",
        )
        .heading("Catalog"),
        switch("replace", "Search again and update entries already present").heading("Catalog"),
        option(
            "max-candidates",
            "N",
            "Inspect at most N candidate configurations",
        )
        .heading("Recurrence search"),
        switch("verbose", "Print candidate search details").heading("Recurrence search"),
        switch(
            "modular-prefilter",
            "Enable the modular rejection (the default)",
        )
        .hidden(),
    ],
    shared: &[&crate::commands::search::RECURRENCE_SEARCH],
    positionals: &[positional("A-NUMBER", ID_HELP).optional().multiple()],
    details: r#"Row widths: --row-widths 1+2 gives rows of 1, 3, 5, ... terms; 1+1/2 gives
1, 1, 2, 2, ... The triangle layout means --row-widths 1+1.

--replace: a replaced entry loses its real_rootedness flag and evidence record.

--bfile needs --name and --layout, and applies to exactly one A-number.
Otherwise OEIS data is required: --oeis-data <dir> or $OEIS_DATA_DIR. It reads
<dir>/seq/A123/A123456.seq (or <dir>/A123/...), and a real b-file
<dir>/files/A123/b123456.txt when it extends the %S/%T/%U terms.

JSON output is polytool.oeis-add.v1."#,
    examples: &[
        "polytool oeis add A008292 --oeis-data ~/oeis --dry-run",
        "polytool oeis add --ids-file ids.txt --max-candidates 20000",
        "polytool oeis add A008292 --oeis-data ~/oeis --recurrence-json rec.json --dry-run",
    ],
    exit_codes: &[
        (0, "every outcome is added, updated or already_present"),
        (1, "some A-number had another outcome"),
        (2, "usage error or unreadable data"),
    ],
    see_also: &["oeis check", "recurrence", "oeis evidence"],
    ..COMMAND
};

const CHECK: CommandSpec = CommandSpec {
    name: "check",
    summary: "Re-validate catalog entries against a local OEIS copy",
    usage: &["oeis check [A-NUMBER...] [--oeis-data <DIR>] [OPTIONS]"],
    description: r#"Regenerate catalog entries and compare their flattened rows with every term of
the local OEIS record (the %S/%T/%U data, or a real b-file when it extends
them). Entries whose flattened alignment is not verified, or whose record is
missing locally, are reported as unchecked."#,
    input: InputKind::None,
    output: OutputKind::JsonSwitch,
    flags: &[
        option(
            "oeis-data",
            "DIR",
            "Local OEIS internal-format copy (default: $OEIS_DATA_DIR)",
        ),
        option(
            "catalog",
            "FILE",
            "Check the entries of a data file instead of the bundled ones",
        ),
        option("ids-file", "FILE", "Restrict to the listed A-numbers"),
    ],
    positionals: &[positional("A-NUMBER", ID_HELP).optional().multiple()],
    details: r#"--catalog is useful after `oeis add`, before rebuilding polytool.
JSON output is polytool.oeis-check.v1."#,
    examples: &[
        "polytool oeis check --oeis-data ~/oeis",
        "polytool oeis check A008292 A008293 --oeis-data ~/oeis --json",
    ],
    exit_codes: &[
        (0, "no entry mismatches"),
        (1, "some entry mismatches"),
        (2, "usage error or unreadable data"),
    ],
    see_also: &["oeis add", "oeis evidence"],
    ..COMMAND
};

const EVIDENCE: CommandSpec = CommandSpec {
    name: "evidence",
    summary: "Strongest coefficient property on the first rows of each entry",
    usage: &["oeis evidence [A-NUMBER...|--all] --count <N> [OPTIONS]"],
    description: r#"Compute, for each entry, the strongest coefficient property that holds on every
one of its first N rows (fewer if the entry cannot generate that many),
exactly, from the catalog's own rows. Properties, strongest first: real_rooted,
ultra_log_concave, log_concave, unimodal, none. real_rooted is the exact
polytool check (any coefficient signs); the other three need a row with
nonnegative coefficients and no zero between two nonzero ones, and then use
polytool's own log-concavity, ultra-log-concavity and unimodality checks
(ultra-log-concavity at the row's degree). The record lists the first failing
row (displayed OEIS row index) of each stronger property.

This is evidence about finitely many rows, not a proof, and it is separate from
the real_rootedness flag."#,
    input: InputKind::None,
    output: OutputKind::JsonSwitch,
    flags: &[
        option(
            "count",
            "N",
            "Number of rows to check (required, at least 1)",
        )
        .deprecated(&["rows"]),
        switch("all", "Every catalog entry (instead of listing A-numbers)"),
        option("ids-file", "FILE", "Restrict to the listed A-numbers"),
        switch("write", "Store the records in the catalog data file"),
        option(
            "per-entry-timeout-seconds",
            "S",
            "Limit per entry, checked before each row (0 = none)",
        )
        .default("60"),
        option(
            "as-of",
            "YYYY-MM-DD",
            "Date to record (default: today, UTC)",
        ),
        option(
            "catalog",
            "FILE",
            "Data file to read (and to write with --write)",
        ),
    ],
    positionals: &[positional("A-NUMBER", ID_HELP).optional().multiple()],
    details: r#"Give either A-numbers (or --ids-file) or --all, not both.

--write changes only the `evidence` key of each computed entry. It needs a data
file: the default is ./data/oeis-catalog.jsonl if present, else the bundled
catalog, and then --write is not possible.

An entry over --per-entry-timeout-seconds gets no record (an existing record is
kept) and is listed.

Besides the per-entry lines the report gives the distribution of `strongest` and
lists every lean_proved entry whose evidence is weaker than real_rooted (a
contradiction to investigate). JSON output is polytool.oeis-evidence.v1."#,
    examples: &[
        "polytool oeis evidence A008292 --count 20",
        "polytool oeis evidence --all --count 12 --catalog data/oeis-catalog.jsonl --write",
    ],
    see_also: &["oeis check", "real-rooted", "properties"],
    ..COMMAND
};

pub(crate) const OEIS: CommandSpec = CommandSpec {
    name: "oeis",
    group: "OEIS families",
    summary: "List and generate recurrence-backed OEIS families",
    description: r#"Work with the OEIS polynomial families bundled with polytool, each generated by
an exact recurrence. `list`, `info` and `generate` use the bundled catalog;
`add`, `check` and `evidence` maintain the catalog data file in a source
checkout."#,
    subcommands: &[LIST, INFO, GENERATE, ADD, CHECK, EVIDENCE],
    examples: &[
        "polytool oeis list",
        "polytool oeis generate A008292 --count 6 --format triangle",
        "polytool oeis add A008292 --oeis-data ~/oeis --dry-run",
    ],
    see_also: &["recurrence-generate", "sequence", "topic:output-formats"],
    run: Some(crate::cmd_oeis),
    ..COMMAND
};
