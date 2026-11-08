//! Declarative command, flag and topic specification.
//!
//! One [`CommandSpec`] per command (and per nested subcommand) is the single
//! source for the clap parser, every help page, shell completion, the
//! machine-readable `polytool help --json` index and the documentation site.
//! The tables live in `commands.rs` and its children; this module defines the
//! types and the `const` builder helpers that keep those tables short.

/// What a command reads from standard input (or from `--input <file|->`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum InputKind {
    /// Reads nothing; `--input` is not offered.
    None,
    /// One polynomial per line: integer coefficient lists or expanded notation.
    Polynomials,
    /// One exact rational coefficient list per line.
    RationalLists,
    /// A JSON document whose shape the command page describes.
    Json,
    /// Another text format described on the command page.
    Text,
}

impl InputKind {
    pub(crate) fn reads_input(self) -> bool {
        self != Self::None
    }

    pub(crate) fn json_name(self) -> Option<&'static str> {
        match self {
            Self::None => None,
            Self::Polynomials => Some("polynomials"),
            Self::RationalLists => Some("rational-lists"),
            Self::Json => Some("json"),
            Self::Text => Some("text"),
        }
    }
}

/// How a command selects its output format.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum OutputKind {
    /// No structured output; neither `--json` nor `--format` is offered.
    None,
    /// Text by default; `--json` (or `--format json`) switches to JSON.
    JsonSwitch,
    /// JSON only (the default); `--json` and `--format json` are accepted no-ops.
    JsonOnly,
    /// The command's own `--format` values; `--json` means `--format json`.
    Formats {
        values: &'static [&'static str],
        default: &'static str,
    },
}

impl OutputKind {
    /// The values accepted by `--format`, if the command has the option.
    pub(crate) fn format_values(self) -> &'static [&'static str] {
        match self {
            Self::None => &[],
            Self::JsonSwitch => &["text", "json"],
            Self::JsonOnly => &["json"],
            Self::Formats { values, .. } => values,
        }
    }

    pub(crate) fn default_format(self) -> Option<&'static str> {
        match self {
            Self::None => None,
            Self::JsonSwitch => Some("text"),
            Self::JsonOnly => Some("json"),
            Self::Formats { default, .. } => Some(default),
        }
    }
}

/// One `--long` option or switch.
#[derive(Clone, Copy, Debug)]
pub(crate) struct FlagSpec {
    /// Canonical long name without the leading dashes.
    pub name: &'static str,
    /// Value placeholder (`N`, `FILE`); `None` for a switch.
    pub value: Option<&'static str>,
    /// Number of values one occurrence takes (1 for ordinary options).
    pub arity: usize,
    pub help: &'static str,
    /// Documented default, shown in help (never injected into the arguments).
    pub default: Option<&'static str>,
    /// Accepted values (validated by the parser, offered by completion).
    pub choices: &'static [&'static str],
    /// Old spellings that still work and print a deprecation warning.
    pub deprecated: &'static [&'static str],
    /// Help section heading; `None` puts the flag under "Options".
    pub heading: Option<&'static str>,
    /// Accepted but not shown in help.
    pub hidden: bool,
}

/// A switch without a value.
pub(crate) const fn switch(name: &'static str, help: &'static str) -> FlagSpec {
    FlagSpec {
        name,
        value: None,
        arity: 0,
        help,
        default: None,
        choices: &[],
        deprecated: &[],
        heading: None,
        hidden: false,
    }
}

/// An option taking one value.
pub(crate) const fn option(
    name: &'static str,
    value: &'static str,
    help: &'static str,
) -> FlagSpec {
    FlagSpec {
        name,
        value: Some(value),
        arity: 1,
        help,
        default: None,
        choices: &[],
        deprecated: &[],
        heading: None,
        hidden: false,
    }
}

impl FlagSpec {
    pub(crate) const fn default(self, default: &'static str) -> Self {
        Self {
            default: Some(default),
            ..self
        }
    }

    pub(crate) const fn choices(self, choices: &'static [&'static str]) -> Self {
        Self { choices, ..self }
    }

    pub(crate) const fn deprecated(self, deprecated: &'static [&'static str]) -> Self {
        Self { deprecated, ..self }
    }

    pub(crate) const fn heading(self, heading: &'static str) -> Self {
        Self {
            heading: Some(heading),
            ..self
        }
    }

    pub(crate) const fn hidden(self) -> Self {
        Self {
            hidden: true,
            ..self
        }
    }

    /// Take `arity` values per occurrence (`--box R0 R1 I0 I1`).
    pub(crate) const fn values(self, arity: usize) -> Self {
        Self { arity, ..self }
    }
}

/// A positional argument.
#[derive(Clone, Copy, Debug)]
pub(crate) struct PositionalSpec {
    /// Placeholder shown in usage (`NAME`, `A-NUMBER`).
    pub name: &'static str,
    pub help: &'static str,
    pub required: bool,
    /// Accepts several values.
    pub multiple: bool,
    pub choices: &'static [&'static str],
}

pub(crate) const fn positional(name: &'static str, help: &'static str) -> PositionalSpec {
    PositionalSpec {
        name,
        help,
        required: true,
        multiple: false,
        choices: &[],
    }
}

impl PositionalSpec {
    pub(crate) const fn optional(self) -> Self {
        Self {
            required: false,
            ..self
        }
    }

    pub(crate) const fn multiple(self) -> Self {
        Self {
            multiple: true,
            ..self
        }
    }

    pub(crate) const fn choices(self, choices: &'static [&'static str]) -> Self {
        Self { choices, ..self }
    }
}

/// A named set of flags shared by several commands (one help section).
#[derive(Debug)]
pub(crate) struct FlagGroup {
    pub heading: &'static str,
    pub flags: &'static [FlagSpec],
}

/// One documented exit status.
pub(crate) type ExitCode = (i32, &'static str);

/// A command or nested subcommand.
pub(crate) struct CommandSpec {
    pub name: &'static str,
    /// Equivalent names shown in help and completion (`gamma`).
    pub aliases: &'static [&'static str],
    /// Old command names that still work and print a deprecation warning.
    pub deprecated_aliases: &'static [&'static str],
    /// Overview heading (top-level commands only).
    pub group: &'static str,
    /// One line, no final period.
    pub summary: &'static str,
    /// Synopsis lines without the leading `polytool `; empty: generated.
    pub usage: &'static [&'static str],
    /// Paragraphs shown after the summary on the full help page.
    pub description: &'static str,
    pub input: InputKind,
    pub output: OutputKind,
    /// Command-specific flags (section "Options" unless a flag sets a heading).
    pub flags: &'static [FlagSpec],
    /// Shared flag groups, each shown under its own heading.
    pub shared: &'static [&'static FlagGroup],
    pub positionals: &'static [PositionalSpec],
    pub subcommands: &'static [CommandSpec],
    /// Long reference text (request fields, output notes) for the full page.
    pub details: &'static str,
    /// Example command lines (shell syntax, without a prompt).
    pub examples: &'static [&'static str],
    pub exit_codes: &'static [ExitCode],
    /// Related commands and topics (`recurrence-generate`, `topic:indexing`).
    pub see_also: &'static [&'static str],
    /// Handler; receives the normalized arguments after the command name.
    /// `None` for nested subcommands, which the parent's handler dispatches.
    pub run: Option<fn(&[String])>,
}

/// Defaults for the optional fields of a [`CommandSpec`].
pub(crate) const COMMAND: CommandSpec = CommandSpec {
    name: "",
    aliases: &[],
    deprecated_aliases: &[],
    group: "",
    summary: "",
    usage: &[],
    description: "",
    input: InputKind::None,
    output: OutputKind::None,
    flags: &[],
    shared: &[],
    positionals: &[],
    subcommands: &[],
    details: "",
    examples: &[],
    exit_codes: &[],
    see_also: &[],
    run: None,
};

impl CommandSpec {
    /// Every flag of the command, shared groups last, with its heading.
    pub(crate) fn all_flags(
        &self,
    ) -> impl Iterator<Item = (&'static FlagSpec, Option<&'static str>)> + 'static {
        self.flags
            .iter()
            .map(|flag| (flag, flag.heading))
            .chain(self.shared.iter().flat_map(|group| {
                group
                    .flags
                    .iter()
                    .map(move |flag| (flag, Some(flag.heading.unwrap_or(group.heading))))
            }))
    }

    pub(crate) fn find_subcommand(&self, name: &str) -> Option<&'static CommandSpec> {
        self.subcommands
            .iter()
            .find(|sub| sub.name == name || sub.aliases.contains(&name))
    }
}

/// Exit statuses shared by commands that always answer.
pub(crate) const EXIT_ANSWERED: &[ExitCode] = &[
    (
        0,
        "answered (also when the mathematical answer is negative)",
    ),
    (2, "usage error or invalid input"),
];

/// A `polytool help <topic>` page.
pub(crate) struct TopicSpec {
    pub name: &'static str,
    pub title: &'static str,
    pub body: &'static str,
}
