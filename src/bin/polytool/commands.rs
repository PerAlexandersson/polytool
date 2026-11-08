//! The command registry: every command's specification, in overview order.
//!
//! Each group of commands has a child module with its `CommandSpec`s (flags,
//! help text, examples, exit statuses); `search` holds the recurrence-search
//! flags shared by several commands. The types are in `cli::spec`.

mod developer;
mod expansions;
mod inspect;
mod networks;
mod oeis;
mod recurrence;
pub(crate) mod search;
mod sieving;

use crate::cli::spec::{positional, CommandSpec, OutputKind, COMMAND};

/// `polytool completions <shell>`.
const COMPLETIONS: CommandSpec = CommandSpec {
    name: "completions",
    group: "Shell support",
    summary: "Generate shell completion (commands, subcommands, flags, topics)",
    description: "Print a completion script for bash, zsh or fish. It completes command \
names, nested subcommands, every flag, flag values with a fixed set of choices, and \
`polytool help <command|topic>`.",
    output: OutputKind::None,
    positionals: &[positional("SHELL", "Shell to generate for").choices(&["bash", "zsh", "fish"])],
    examples: &[
        "polytool completions bash > ~/.local/share/bash-completion/completions/polytool",
        "polytool completions zsh > \"${fpath[1]}/_polytool\"",
        "polytool completions fish > ~/.config/fish/completions/polytool.fish",
    ],
    exit_codes: &[(0, "script printed"), (2, "unknown shell")],
    run: Some(crate::cli::print_completions),
    ..COMMAND
};

pub(crate) const COMMANDS: &[CommandSpec] = &[
    inspect::REAL_ROOTED,
    inspect::INTERLACING,
    inspect::INTERLACING_PROFILE,
    inspect::PROPERTIES,
    inspect::FAMILY_CHECK,
    inspect::COEFFICIENT_TESTS,
    inspect::ENUMERATE_REAL_ROOTED,
    inspect::ROOT_ENCLOSURES,
    inspect::RESULTANT,
    inspect::DISCRIMINANT,
    recurrence::SEQUENCE,
    recurrence::RECURRENCE,
    recurrence::RECURRENCE_GENERATE,
    recurrence::RECURRENCE_TRANSFORM,
    recurrence::RECURRENCE_DEGREE,
    recurrence::RECURRENCE_PROVE,
    recurrence::RECURRENCE_THEOREMS,
    recurrence::RECURRENCE_VERONESE_LIFTS,
    recurrence::RECURRENCE_IMPLICIT,
    recurrence::RECURRENCE_SEARCH_REPRESENTATIONS,
    oeis::OEIS,
    expansions::GAMMA_EXPANSION,
    expansions::BERNSTEIN_EXPANSION,
    expansions::HSTAR_TO_EHRHART,
    expansions::EHRHART_TO_HSTAR,
    expansions::HSTAR_INEQUALITIES,
    expansions::STAPLEDON,
    sieving::CYCLIC_SIEVING,
    sieving::CYCLIC_SIEVING_SEQUENCE,
    sieving::BICYCLIC_SIEVING,
    networks::TNN_NETWORK,
    networks::BRENTI_RECOGNIZE,
    networks::PF_NETWORK,
    networks::LEFT_PRODUCTION,
    developer::BENCH,
    developer::BKW_SCOUT,
    developer::PF_PENCIL,
    COMPLETIONS,
];

pub(crate) fn find(name: &str) -> Option<&'static CommandSpec> {
    COMMANDS.iter().find(|cmd| {
        cmd.name == name || cmd.aliases.contains(&name) || cmd.deprecated_aliases.contains(&name)
    })
}
