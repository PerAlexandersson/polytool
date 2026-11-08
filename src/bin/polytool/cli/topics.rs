//! `polytool help <topic>` pages: concepts shared by many commands.
//!
//! The `exit-codes` and `compatibility` pages end with tables generated from
//! the command specification, so they cannot drift from the commands.

use crate::cli::spec::{TopicSpec, EXIT_ANSWERED};
use crate::commands::COMMANDS;

pub(crate) const TOPICS: &[TopicSpec] = &[
    TopicSpec {
        name: "getting-started",
        title: "A five-minute tour of polytool",
        body: GETTING_STARTED,
    },
    TopicSpec {
        name: "input",
        title: "Input formats: polynomials, rationals, JSON, --input",
        body: INPUT,
    },
    TopicSpec {
        name: "coefficient-order",
        title: "Coefficient order: ascending, and the three exceptions",
        body: COEFFICIENT_ORDER,
    },
    TopicSpec {
        name: "indexing",
        title: "Indices: index, --first-index and --from-index",
        body: INDEXING,
    },
    TopicSpec {
        name: "output-formats",
        title: "Output: --format, --json, --compact, stdout and stderr",
        body: OUTPUT_FORMATS,
    },
    TopicSpec {
        name: "json-schemas",
        title: "JSON schema identifiers and their versioning policy",
        body: JSON_SCHEMAS,
    },
    TopicSpec {
        name: "exit-codes",
        title: "Exit statuses of every command",
        body: EXIT_CODES,
    },
    TopicSpec {
        name: "budgets",
        title: "Search budgets and what happens when one is exhausted",
        body: BUDGETS,
    },
    TopicSpec {
        name: "recurrence-json",
        title: "The polytool.recurrence.v1 record, wrappers and menus",
        body: RECURRENCE_JSON,
    },
    TopicSpec {
        name: "compatibility",
        title: "Renamed commands and flags, and their deprecated spellings",
        body: COMPATIBILITY,
    },
];

pub(crate) fn find(name: &str) -> Option<&'static TopicSpec> {
    TOPICS.iter().find(|topic| topic.name == name)
}

pub(crate) fn render(topic: &TopicSpec) -> String {
    let mut text = format!(
        "{}\n\n{}\n",
        topic.title,
        crate::cli::wrap_paragraphs(topic.body.trim_end())
    );
    match topic.name {
        "exit-codes" => text.push_str(&exit_code_table()),
        "compatibility" => text.push_str(&deprecated_spelling_table()),
        _ => {}
    }
    text
}

fn exit_code_table() -> String {
    let mut text = String::from("\nPer command:\n");
    for command in COMMANDS {
        let mut pages = vec![(command.name.to_string(), command)];
        for sub in command.subcommands {
            pages.push((format!("{} {}", command.name, sub.name), sub));
        }
        for (name, spec) in pages {
            let codes = if spec.exit_codes.is_empty() {
                EXIT_ANSWERED
            } else {
                spec.exit_codes
            };
            let codes: Vec<String> = codes
                .iter()
                .map(|(code, meaning)| format!("{code} {meaning}"))
                .collect();
            text.push_str(&format!("  {name}\n      {}\n", codes.join("\n      ")));
        }
    }
    text
}

fn deprecated_spelling_table() -> String {
    let mut text = String::from("\nDeprecated spellings (old -> new):\n");
    for command in COMMANDS {
        let mut pages = vec![(command.name.to_string(), command)];
        for sub in command.subcommands {
            pages.push((format!("{} {}", command.name, sub.name), sub));
        }
        for (name, spec) in pages {
            for alias in spec.deprecated_aliases {
                text.push_str(&format!("  {alias} -> {name}\n"));
            }
            let mut flags: Vec<String> = Vec::new();
            for (flag, _) in spec.all_flags() {
                for old in flag.deprecated {
                    flags.push(format!("--{old} -> --{}", flag.name));
                }
            }
            if !flags.is_empty() {
                text.push_str(&format!("  {name}:\n"));
                for line in flags {
                    text.push_str(&format!("      {line}\n"));
                }
            }
        }
    }
    text
}

const GETTING_STARTED: &str = "\
Every command reads its input from stdin (or from --input FILE) and writes its \
answer to stdout; progress notes, warnings and errors go to stderr. Polynomials \
are dense coefficient lists in ascending degree order, one per line.

Generate a family, check it, and look for a recurrence:
  polytool sequence eulerian 8 --format rows > eulerian.txt
  polytool real-rooted --input eulerian.txt
  polytool properties --json < eulerian.txt
  polytool family-check --require-real-rooted --require-weak-interlacing < eulerian.txt
  polytool recurrence --json < eulerian.txt > recurrence.json
  polytool recurrence-generate --recurrence recurrence.json --count 10

Work with the bundled OEIS families:
  polytool oeis list
  polytool oeis info A008292
  polytool oeis generate A008292 --count 8 --format rows | polytool interlacing

Pipe JSON between commands: --json (or --format json) gives a stable, versioned \
JSON report (see `polytool help json-schemas`); exact integers are JSON strings.

Next steps: `polytool help <command>` for any command, `polytool help input`, \
`polytool help indexing` and `polytool help exit-codes`.";

const INPUT: &str = "\
Polynomial input (most commands): one polynomial per line on stdin, or in the \
file named by --input FILE (`--input -` is stdin).
  1, 3, 2              dense coefficients a_0, a_1, ..., a_d (ascending)
  [1, 3, 2]            brackets ([], {}, ()) are allowed
  1 3 2                whitespace-separated coefficients
  1 + 3t + 2t^2        expanded notation (integer-coefficient commands)
  # a comment          blank lines and lines starting with # are skipped

Coefficients are arbitrary-size integers. Commands that read rational \
coefficient lists (recurrence, bernstein-expansion, ehrhart-to-hstar) also \
accept exact rationals such as 3/7, but not expanded notation.

Errors name the physical input line (`line 5: coefficient 2: ...`), while \
reported indices count only polynomial lines. The input is limited to 16 MiB \
and a polynomial to the library's coefficient limit.

JSON input: root-enclosures, recurrence-transform, recurrence-implicit, \
recurrence-search-representations, bicyclic-sieving and `recurrence --space` \
read one JSON document; their help pages list the request fields. Commands \
that consume a recurrence take --recurrence FILE (or `-` for stdin), which may \
be a polytool.recurrence.v1 record, a theorem wrapper or a representation menu \
(see `polytool help recurrence-json`).";

const COEFFICIENT_ORDER: &str = "\
Coefficient lists are in ascending degree order everywhere: entry i is the \
coefficient of t^i, so `1, 3, 2` is 1 + 3t + 2t^2. JSON reports use the same \
order. The exceptions are explicit in their names:

  enumerate-real-rooted --leading-coefficients-descending A0,A1,...
      gives the leading coefficients first (descending); its output is ascending.
  recurrence-generate --descending
      prints each generated row highest degree first.
  bkw-scout --symbol '1; -x; 1'
      lists the coefficients of z^0, z^1, ... (ascending in z); each is a
      polynomial in x.

Triangles (oeis generate --format triangle, b-files) list row n, then row n+1, \
each row in ascending order.";

const INDEXING: &str = "\
An index is the position of a polynomial P_n in a sequence; a row is a row of a \
triangle or table. Two options set indices, with different meanings:

  --first-index K   labels the first input row: it is P_K. Defaults differ by
                    command and are shown in each help page: recurrence
                    numbers the retained rows from 1 (rows dropped by
                    --skip-prefix still count when --first-index is given),
                    cyclic-sieving-sequence starts at 0, brenti-recognize,
                    pf-network and left-production use 1 (P_0 = 1 implicit).
  --from-index K    selects where the output starts (recurrence-generate,
                    oeis generate); earlier indices are computed but not shown.

Reports: a zero-based `index` is a position in the input list (real-rooted, \
properties, ... JSON items); `first_index` is the sequence index of the first \
polynomial (recurrence JSON, transform schedules, root-enclosures). OEIS row \
JSON carries both `first_row` (old name) and `first_index` with equal values.";

const OUTPUT_FORMATS: &str = "\
Every command with structured output accepts --format FORMAT and the shorthand \
--json (= --format json); `polytool help <command>` lists its formats and the \
default. Text output is for reading; JSON output is for scripts and is \
versioned (see `polytool help json-schemas`).

JSON-only commands (root-enclosures, recurrence-transform, recurrence-implicit, \
recurrence-search-representations, recurrence --space/--verify-space) print \
pretty JSON by default and accept --compact for one line; --json is accepted \
there as a no-op.

Exact values: integers and rationals that can exceed 2^53 are JSON strings \
(\"12345678901234567890\", \"3/7\"); small counts and indices are numbers.

Streams: the answer goes to stdout; progress notes, deprecation warnings and \
errors go to stderr, so `polytool ... --json | jq` stays valid. A closed pipe \
(`polytool ... | head`) ends the program quietly with status 141.";

const JSON_SCHEMAS: &str = "\
Machine-readable reports carry a `schema` key `polytool.<name>.v<N>`. Within a \
version, keys are only added, never removed or renamed; a renamed key keeps \
being emitted under its old name next to the new one. A breaking change bumps \
the version. Simple per-polynomial reports (real-rooted, interlacing, \
properties, gamma-expansion, ...) have the shape {\"items\": [...]} without a \
schema key.

  polytool.commands.v1                     polytool help --json
  polytool.recurrence.v1                   recurrence --json (a found recurrence)
  polytool.recurrence-representations.v1   recurrence --representations --json
  polytool.recurrence-spaces.v1            recurrence --space / --verify-space
  polytool.recurrence-theorems.v1          recurrence --theorem-matches --json,
                                           recurrence-theorems
  polytool.recurrence-leading-coefficient.v1  recurrence-degree
  polytool.recurrence-proof-discovery.v1   recurrence-prove
  polytool.recurrence-real-rootedness-theorem.v1  recurrence-prove reports
  polytool.recurrence-veronese-lifts.v1    recurrence-veronese-lifts
  polytool.recurrence-veronese-certificate.v1  recurrence-veronese-lifts
  polytool.recurrence-transform.v1         recurrence-transform
  polytool.recurrence-implicit.v1          recurrence-implicit
  polytool.recurrence-implicit-search.v1   recurrence-implicit --search
  polytool.recurrence-representation-search.v1  recurrence-search-representations
  polytool.image_recurrence.v1             implicit image rules
  polytool.recurrence-*-lift(s).v1         lift searches and certificates
  polytool.root-enclosures.v1              root-enclosures
  polytool.oeis-catalog.v1                 oeis list --json, oeis info --json
  polytool.oeis-rows.v1                    oeis generate --json
  polytool.oeis-text.v1                    recurrence-generate --format oeis-json
  polytool.oeis-add.v1, oeis-check.v1, oeis-evidence.v1  oeis add/check/evidence
  polytool.brenti-recognize.v1             brenti-recognize --json
  polytool.pf-network.v1                   pf-network --json
  polytool.left-production.v1              left-production --json
  polytool.tnn-network.v1                  tnn-network --json
  polytool.pf-pencil.v1                    pf-pencil --json
  polytool.bench.real-rootedness.v1        bench real-rootedness --json";

const EXIT_CODES: &str = "\
General rule:
  0  the command answered, also when the mathematical answer is negative
     (a polynomial that is not real-rooted, a search that found nothing)
  1  a decided negative where the command documents it (family-check with
     --require-*, failed certificates and replays, impossible networks)
  2  usage error or invalid input; stdout stays empty
  3  inconclusive: a search budget was exhausted before a decision
  141 stdout was closed early (broken pipe)";

const BUDGETS: &str = "\
Searches are bounded by deterministic budgets: the same budget gives the same \
answer on every machine (time limits are only used by `oeis evidence`).

  --max-candidates N        recurrence, oeis add, recurrence-veronese-lifts:
                            candidate configurations inspected
  --max-work N              recurrence: deterministic work units (solve sizes)
  --max-search-nodes N      brenti-recognize, pf-network: search-tree nodes
  --max-integer-search-nodes N  brenti-recognize, pf-network: integer search
                            per column (old spelling --max-integer-nodes)
  --max-results, --max-nodes, --max-coefficient-bits
                            enumerate-real-rooted
  --max-iterations N        bkw-scout: Durand--Kerner iterations
  --max-start-delay N       recurrence-prove
  --per-entry-timeout-seconds S  oeis evidence (wall clock)
  limits {...}              JSON requests of recurrence --space,
                            recurrence-search-representations, recurrence-implicit

When a budget is exhausted the answer is never a negative: the report says so \
(`budget_exhausted`, `inconclusive`, `incomplete_budget`, `resource_limit`). \
recurrence, recurrence --representations, recurrence --space, brenti-recognize \
and pf-network then exit with status 3; enumerate-real-rooted exits 0 with status \
resource_limit and an initial segment of the sorted list. Each command's help \
page lists its exit statuses.";

const RECURRENCE_JSON: &str = "\
`polytool recurrence --json` prints a polytool.recurrence.v1 record:
  schema                 \"polytool.recurrence.v1\"
  first_index            index n of the first initial polynomial
  initial_polynomials    seed rows (ascending coefficient strings)
  recurrence.terms       one term per (offset, derivative order, sign) with an
                         exact coefficient polynomial in t and n
  search                 options and diagnostics of the search that found it

Commands that consume a recurrence (recurrence-generate, recurrence-degree, \
recurrence-prove, recurrence-theorems, recurrence-veronese-lifts) take \
--recurrence FILE or --recurrence - (stdin) and also accept two wrappers: a \
polytool.recurrence-theorems.v1 report (its recurrence_json is used) and a \
polytool.recurrence-representations.v1 menu, from which --candidate N \
(zero-based) selects one entry; a menu with several entries requires \
--candidate. Extracting a recurrence from a report does not verify the \
report's claims.

A record fitted to rows is a candidate: it reproduces the rows it was fitted \
and verified on, which is evidence, not a proof for all n.";

const COMPATIBILITY: &str = "\
Old spellings keep working. Each use prints one warning on stderr (stdout is \
unchanged); set POLYTOOL_DEPRECATION_WARNINGS=0 to silence the warnings. The \
spellings used by published scripts (the recurrence search flags, --rows, \
--first-row, --start-row, --offset and the recurrence-* commands) are frozen: \
they will not be removed in the 0.x series.

Naming rules: names are spelled out; `index` is a position in a sequence and \
`row` a row of a triangle; --first-index labels the first input row and \
--from-index selects where output starts; --count is the number of rows or \
polynomials to produce; options that widen a search start with --try-; limits \
start with --max-. A renamed JSON key is emitted under both names, and the MCP \
server accepts both request field names.";
