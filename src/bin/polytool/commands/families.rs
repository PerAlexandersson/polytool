//! Specifications for the curated family catalog.
use crate::cli::spec::{option, positional, CommandSpec, OutputKind, COMMAND};
const LIST: CommandSpec = CommandSpec {
    name: "list",
    summary: "List OEIS and submission-bound sequence families",
    flags: &[
        option("query", "TEXT", "Case-insensitive ID/title filter"),
        option("status", "STATUS", "Real-rootedness status filter"),
    ],
    output: OutputKind::JsonSwitch,
    examples: &["polytool families list --status theorem_proved --json"],
    ..COMMAND
};
const SHOW: CommandSpec = CommandSpec {
    name: "show",
    summary: "Show complete curated metadata and canonical generator",
    positionals: &[positional("ID", "A###### or PT-F####")],
    output: OutputKind::JsonOnly,
    examples: &["polytool families show PT-F0001"],
    ..COMMAND
};
const GENERATE: CommandSpec = CommandSpec {
    name: "generate",
    summary: "Generate the first complete exact coefficient rows",
    positionals: &[positional("ID", "A###### or PT-F####")],
    flags: &[option("count", "N", "Number of complete rows, 1..512").default("8")],
    output: OutputKind::Formats {
        values: &["rows", "terms", "json"],
        default: "rows",
    },
    examples: &["polytool families generate PT-F0001 --count 9 --format rows"],
    ..COMMAND
};
const EXPORT: CommandSpec = CommandSpec {
    name: "export-oeis", summary: "Print copyable OEIS draft blocks and a proposal TOML skeleton",
    description: "Generate a conservative 240-byte data prefix in at most three 80-byte OEIS lines, plus a b-file bounded to 10,000 terms and 1 MiB. Review authors, search date, keywords and generator adapter in the proposal skeleton. This command does not submit anything.",
    positionals: &[positional("ID","A###### or PT-F####")],
    flags: &[option("count","N","Number of complete rows to generate, 1..512").default("12"), option("bfile","FILE","Write b-file here (create new file; never overwrite)")],
    output: OutputKind::JsonSwitch,
    examples: &["polytool families export-oeis PT-F0001 --count 40 --bfile packing.txt", "polytool families export-oeis PT-F0002 --json"], ..COMMAND
};
const VALIDATE: CommandSpec = CommandSpec {
    name: "validate",
    summary: "Validate all bundled sidecars, references and reproduced rows",
    flags: &[option(
        "lean-export",
        "FILE",
        "Optional public proof-owner JSON declaration list at pinned commits",
    )],
    output: OutputKind::JsonSwitch,
    examples: &["polytool families validate --json"],
    ..COMMAND
};
pub(crate) const FAMILIES: CommandSpec = CommandSpec {
    name:"families", group:"OEIS families", summary:"Curated sequence families, proof attestations and OEIS submission drafts",
    description:"One overview of bundled OEIS families and curated submission-bound families. Canonical OEIS generators/statuses are inherited, with curated JSON metadata and non-OEIS exp-rows or recurrence generators. This catalog is separate from an automated discovery census.",
    examples: &["polytool families list --query packing", "polytool families generate PT-F0001 --count 9", "polytool families export-oeis PT-F0002 --json"],
    subcommands:&[LIST,SHOW,GENERATE,EXPORT,VALIDATE], run:Some(crate::families::cmd_families),
    see_also:&["oeis","exp-rows"], ..COMMAND
};
