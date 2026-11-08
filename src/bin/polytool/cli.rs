//! The clap front end: builds the parser from the command specification,
//! normalizes the parsed arguments for the command handlers, prints help
//! pages and topics, the `help --json` index and shell completions.
//!
//! Parsing contract:
//!
//! - `--flag value` and `--flag=value` are equivalent; values may start with
//!   `-` (negative numbers, BKW symbols).
//! - A deprecated spelling is accepted, mapped to its canonical name and
//!   reported once on stderr (set `POLYTOOL_DEPRECATION_WARNINGS=0` to silence).
//! - Handlers receive canonical tokens in the original order
//!   (`--name value...` or positional values), so the per-command semantics
//!   (last occurrence wins, repeated options accumulate) are unchanged.
//! - `--input <file|->` replaces standard input for every reading command;
//!   `--json` and `--format <...>` are offered by every command with
//!   structured output and normalized to what the handler understands.
//! - Usage errors go to stderr with exit status 2 and end with a pointer to
//!   `polytool help <command>`; stdout stays empty.

pub(crate) mod spec;
pub(crate) mod topics;

use crate::cli::spec::{CommandSpec, FlagSpec, InputKind, OutputKind};
use crate::commands::COMMANDS;
use clap::builder::PossibleValuesParser;
use clap::error::{ContextKind, ContextValue, ErrorKind};
use clap::parser::ValueSource;
use clap::{Arg, ArgAction, ArgMatches, Command};
use serde_json::{json, Value};
use std::io::Write;

const DEPRECATION_WARNINGS_VARIABLE: &str = "POLYTOOL_DEPRECATION_WARNINGS";
const START_HERE: &[&str] = &[
    "real-rooted",
    "properties",
    "family-check",
    "recurrence",
    "oeis",
];

fn exit_usage() -> ! {
    std::process::exit(2)
}

// ---------------------------------------------------------------------------
// Building the clap commands.

fn flag_help(flag: &FlagSpec) -> String {
    let mut text = flag.help.to_string();
    if let Some(default) = flag.default {
        text.push_str(&format!(" [default: {default}]"));
    }
    protect(text)
}

fn value_arg(arg: Arg, flag: &'static FlagSpec) -> Arg {
    let names: Vec<&'static str> = flag.value.unwrap_or("VALUE").split_whitespace().collect();
    let mut arg = arg
        .num_args(flag.arity)
        .action(ArgAction::Append)
        .allow_hyphen_values(true);
    arg = if names.len() > 1 {
        arg.value_names(names)
    } else {
        arg.value_name(names[0])
    };
    if !flag.choices.is_empty() {
        arg = arg.value_parser(PossibleValuesParser::new(flag.choices));
    }
    arg
}

fn flag_args(flag: &'static FlagSpec, heading: Option<&'static str>) -> Vec<Arg> {
    let mut args = Vec::new();
    for (index, id) in std::iter::once(flag.name)
        .chain(flag.deprecated.iter().copied())
        .enumerate()
    {
        let mut arg = Arg::new(id).long(id);
        if index == 0 {
            arg = arg.help(flag_help(flag)).hide(flag.hidden);
            if let Some(heading) = heading {
                arg = arg.help_heading(heading);
            }
        } else {
            arg = arg.hide(true);
        }
        arg = if flag.value.is_some() {
            value_arg(arg, flag)
        } else {
            arg.action(ArgAction::SetTrue)
        };
        args.push(arg);
    }
    args
}

/// The `--input`, `--json` and `--format` options a command gets from its
/// input and output kinds.
fn universal_flags(spec: &CommandSpec) -> Vec<Arg> {
    let mut args = Vec::new();
    let declares = |name: &str| spec.all_flags().any(|(flag, _)| flag.name == name);
    if spec.input.reads_input() && !declares("input") {
        args.push(
            Arg::new("input")
                .long("input")
                .value_name("FILE")
                .num_args(1)
                .action(ArgAction::Append)
                .allow_hyphen_values(true)
                .help("Read the input from FILE instead of stdin (`-` is stdin)")
                .help_heading("Input and output"),
        );
    }
    let values = spec.output.format_values();
    if !values.is_empty() && !declares("format") {
        let default = spec.output.default_format().unwrap_or("text");
        args.push(
            Arg::new("format")
                .long("format")
                .value_name("FORMAT")
                .num_args(1)
                .action(ArgAction::Append)
                .allow_hyphen_values(true)
                .value_parser(PossibleValuesParser::new(values))
                .help(format!("Output format [default: {default}]"))
                .help_heading("Input and output"),
        );
    }
    if !values.is_empty() && !declares("json") {
        let help = match spec.output {
            OutputKind::JsonOnly => "Emit JSON (the default; accepted for uniformity)",
            _ => "Emit machine-readable JSON (same as --format json)",
        };
        args.push(
            Arg::new("json")
                .long("json")
                .action(ArgAction::SetTrue)
                .help(help)
                .help_heading("Input and output"),
        );
    }
    args
}

fn exit_status_text(spec: &CommandSpec) -> String {
    let codes = if spec.exit_codes.is_empty() {
        spec::EXIT_ANSWERED
    } else {
        spec.exit_codes
    };
    let mut text = String::from("Exit status:\n");
    for (code, meaning) in codes {
        text.push_str(&format!("  {code}  {meaning}\n"));
    }
    text
}

fn examples_text(spec: &CommandSpec) -> String {
    if spec.examples.is_empty() {
        return String::new();
    }
    let mut text = String::from("Examples:\n");
    for example in spec.examples {
        for (index, line) in example.lines().enumerate() {
            let prefix = if index == 0 { "  " } else { "    " };
            text.push_str(&format!("{prefix}{line}\n"));
        }
    }
    text
}

fn input_text(spec: &CommandSpec) -> &'static str {
    match spec.input {
        InputKind::None => "",
        InputKind::Polynomials => {
            "Input:\n  One polynomial per line on stdin (or --input FILE): ascending integer\n  \
             coefficient lists such as `1, 3, 2`, or expanded notation such as\n  \
             `1 + 3t + 2t^2`. See `polytool help input`.\n"
        }
        InputKind::RationalLists => {
            "Input:\n  One exact coefficient list per line on stdin (or --input FILE), in\n  \
             ascending degree order; entries are integers or rationals such as 3/7.\n  \
             See `polytool help input`.\n"
        }
        InputKind::Json => {
            "Input:\n  A JSON document on stdin (or --input FILE), described on this page.\n"
        }
        InputKind::Text => {
            "Input:\n  Text on stdin (or --input FILE); see the description above.\n"
        }
    }
}

fn deprecated_text(spec: &CommandSpec) -> String {
    let mut lines = Vec::new();
    for name in spec.deprecated_aliases {
        lines.push(format!("  {name} -> {}", spec.name));
    }
    for (flag, _) in spec.all_flags() {
        for old in flag.deprecated {
            lines.push(format!("  --{old} -> --{}", flag.name));
        }
    }
    if lines.is_empty() {
        return String::new();
    }
    format!(
        "Deprecated spellings (still accepted; they print a warning on stderr):\n{}\n",
        lines.join("\n")
    )
}

fn see_also_text(spec: &CommandSpec) -> String {
    if spec.see_also.is_empty() {
        return String::new();
    }
    let entries: Vec<String> = spec
        .see_also
        .iter()
        .map(|entry| match entry.strip_prefix("topic:") {
            Some(topic) => format!("polytool help {topic}"),
            None => format!("polytool {entry}"),
        })
        .collect();
    format!("See also:\n  {}\n", entries.join(", "))
}

fn join_sections(sections: &[String]) -> String {
    sections
        .iter()
        .filter(|section| !section.is_empty())
        .map(|section| section.trim_end().to_string())
        .collect::<Vec<_>>()
        .join("\n\n")
}

/// Reflow the paragraphs of `text` (runs of unindented lines) and wrap them at
/// 80 columns; indented lines are preformatted (formulas, tables) and kept.
pub(crate) fn wrap_paragraphs(text: &str) -> String {
    const WIDTH: usize = 80;
    fn flush(paragraph: &mut Vec<&str>, out: &mut Vec<String>) {
        if paragraph.is_empty() {
            return;
        }
        let mut current = String::new();
        for word in paragraph.join(" ").split_whitespace() {
            if !current.is_empty() && current.chars().count() + 1 + word.chars().count() > WIDTH {
                out.push(std::mem::take(&mut current));
            }
            if !current.is_empty() {
                current.push(' ');
            }
            current.push_str(word);
        }
        out.push(current);
        paragraph.clear();
    }
    let mut out = Vec::new();
    let mut paragraph = Vec::new();
    for line in text.lines() {
        if line.trim().is_empty() || line.starts_with(' ') {
            flush(&mut paragraph, &mut out);
            out.push(line.to_string());
        } else {
            paragraph.push(line);
        }
    }
    flush(&mut paragraph, &mut out);
    out.join("\n")
}

/// Page layout: summary, description, usage, options, then the sections of
/// [`after_text`]. The summary alone (`about`) also labels the command in its
/// parent's subcommand list.
const HELP_TEMPLATE: &str =
    "{about-with-newline}\n{before-help}{usage-heading} {usage}\n\n{all-args}{after-help}";

/// The description paragraphs shown between the summary and the usage line.
fn description_text(spec: &CommandSpec) -> String {
    if spec.description.is_empty() {
        String::new()
    } else {
        wrap_paragraphs(spec.description.trim_end())
    }
}

/// The sections after the option list.
fn after_text(spec: &CommandSpec) -> String {
    join_sections(&[
        input_text(spec).to_string(),
        if spec.details.is_empty() {
            String::new()
        } else {
            format!(
                "Details:\n{}",
                indent(&wrap_paragraphs(spec.details.trim_end()))
            )
        },
        examples_text(spec),
        exit_status_text(spec),
        deprecated_text(spec),
        see_also_text(spec),
    ])
}

/// clap turns the text `{n}` into a line break in help strings; specification
/// text (LaTeX such as `P_{n}`) is protected with this placeholder and
/// restored by [`restore_braces`] when a page is printed.
const BRACE_N_PLACEHOLDER: char = '\u{E000}';

fn protect(text: String) -> String {
    text.replace("{n}", &BRACE_N_PLACEHOLDER.to_string())
}

pub(crate) fn restore_braces(text: &str) -> String {
    text.replace(BRACE_N_PLACEHOLDER, "{n}")
}

fn build_command(spec: &'static CommandSpec) -> Command {
    // One page for `-h`, `--help` and `polytool help <command>`: no long
    // variants are set, so clap renders the compact layout everywhere.
    let mut command = Command::new(spec.name)
        .about(protect(format!("{}.", spec.summary)))
        .before_help(protect(description_text(spec)))
        .help_template(HELP_TEMPLATE)
        .after_help(protect(after_text(spec)))
        .args_override_self(true)
        .disable_help_subcommand(true)
        .disable_version_flag(true);
    for alias in spec.aliases {
        command = command.visible_alias(*alias);
    }
    for alias in spec.deprecated_aliases {
        command = command.alias(*alias);
    }
    if !spec.usage.is_empty() {
        let usage: Vec<String> = spec
            .usage
            .iter()
            .map(|line| format!("polytool {line}"))
            .collect();
        command = command.override_usage(usage.join("\n       "));
    }
    for (flag, heading) in spec.all_flags() {
        for arg in flag_args(flag, heading) {
            command = command.arg(arg);
        }
    }
    for arg in universal_flags(spec) {
        command = command.arg(arg);
    }
    for (index, positional) in spec.positionals.iter().enumerate() {
        let mut arg = Arg::new(positional.name)
            .value_name(positional.name)
            .index(index + 1)
            .required(positional.required)
            .help(positional.help)
            .action(if positional.multiple {
                ArgAction::Append
            } else {
                ArgAction::Set
            })
            .num_args(if positional.multiple {
                clap::builder::ValueRange::new(1..)
            } else {
                clap::builder::ValueRange::new(1..=1)
            });
        if !positional.choices.is_empty() {
            arg = arg.value_parser(PossibleValuesParser::new(positional.choices));
        }
        command = command.arg(arg);
    }
    if !spec.subcommands.is_empty() {
        command = command.subcommand_required(true);
        for sub in spec.subcommands {
            command = command.subcommand(build_command(sub));
        }
    }
    command
}

fn indent(text: &str) -> String {
    text.lines()
        .map(|line| {
            if line.is_empty() {
                String::new()
            } else {
                format!("  {line}")
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
}

/// The full clap command tree (used for parsing, help pages and completion).
pub(crate) fn root_command() -> Command {
    let mut root = Command::new("polytool")
        .bin_name("polytool")
        .about("Exact univariate polynomial toolkit for combinatorial research")
        .disable_help_subcommand(true)
        .disable_version_flag(true)
        .subcommand_required(true);
    for spec in COMMANDS {
        root = root.subcommand(build_command(spec));
    }
    root
}

// ---------------------------------------------------------------------------
// Normalizing parsed arguments for the handlers.

fn warn_deprecated(command: &str, old: &str, new: &str) {
    if std::env::var(DEPRECATION_WARNINGS_VARIABLE).is_ok_and(|value| value == "0") {
        return;
    }
    eprintln!(
        "polytool: warning: `{old}` is deprecated in `polytool {command}`; use `{new}` \
         (the old spelling keeps working; see `polytool help compatibility`)"
    );
}

fn command_line(matches: &ArgMatches, id: &str) -> bool {
    matches.value_source(id) == Some(ValueSource::CommandLine)
}

fn first_index(matches: &ArgMatches, id: &str) -> usize {
    matches
        .indices_of(id)
        .and_then(|mut indices| indices.next())
        .unwrap_or(0)
}

fn last_index(matches: &ArgMatches, id: &str) -> Option<usize> {
    command_line(matches, id)
        .then(|| matches.indices_of(id).and_then(Iterator::last))
        .flatten()
}

/// Tokens of every occurrence of an option, keyed by position.
fn option_occurrences(
    matches: &ArgMatches,
    id: &str,
    canonical: &str,
    arity: usize,
) -> Vec<(usize, Vec<String>)> {
    let indices: Vec<usize> = matches.indices_of(id).into_iter().flatten().collect();
    let values: Vec<String> = matches
        .get_many::<String>(id)
        .into_iter()
        .flatten()
        .cloned()
        .collect();
    let arity = arity.max(1);
    values
        .chunks(arity)
        .enumerate()
        .map(|(occurrence, chunk)| {
            let position = indices.get(occurrence * arity).copied().unwrap_or(0);
            let mut tokens = vec![format!("--{canonical}")];
            tokens.extend(chunk.iter().cloned());
            (position, tokens)
        })
        .collect()
}

fn normalize(spec: &'static CommandSpec, matches: &ArgMatches, path: &str) -> Vec<String> {
    let mut items: Vec<(usize, Vec<String>)> = Vec::new();
    for (flag, _) in spec.all_flags() {
        let spellings = std::iter::once((flag.name, false))
            .chain(flag.deprecated.iter().map(|old| (*old, true)));
        for (id, deprecated) in spellings {
            if !command_line(matches, id) {
                continue;
            }
            if deprecated {
                warn_deprecated(path, &format!("--{id}"), &format!("--{}", flag.name));
            }
            if flag.value.is_some() {
                items.extend(option_occurrences(matches, id, flag.name, flag.arity));
            } else {
                items.push((first_index(matches, id), vec![format!("--{}", flag.name)]));
            }
        }
    }
    if spec.input.reads_input() && command_line(matches, "input") {
        if let Some(path) = matches
            .get_many::<String>("input")
            .and_then(|mut values| values.next_back())
        {
            crate::input::set_input_path(path);
        }
    }
    normalize_output(spec, matches, &mut items);
    for positional in spec.positionals {
        let indices: Vec<usize> = matches
            .indices_of(positional.name)
            .into_iter()
            .flatten()
            .collect();
        for (offset, value) in matches
            .get_many::<String>(positional.name)
            .into_iter()
            .flatten()
            .enumerate()
        {
            let position = indices.get(offset).copied().unwrap_or(0);
            items.push((position, vec![value.clone()]));
        }
    }
    items.sort_by_key(|(position, _)| *position);
    let mut tokens: Vec<String> = items.into_iter().flat_map(|(_, tokens)| tokens).collect();
    if let Some((name, sub_matches)) = matches.subcommand() {
        let sub = spec
            .find_subcommand(name)
            .expect("clap only returns declared subcommands");
        tokens.push(sub.name.to_string());
        tokens.extend(normalize(sub, sub_matches, &format!("{path} {}", sub.name)));
    }
    tokens
}

/// Translate the universal `--json` / `--format` options into the spelling
/// the handler understands.
fn normalize_output(
    spec: &CommandSpec,
    matches: &ArgMatches,
    items: &mut Vec<(usize, Vec<String>)>,
) {
    if spec.output == OutputKind::None {
        return;
    }
    let declares = |name: &str| spec.all_flags().any(|(flag, _)| flag.name == name);
    let json = (!declares("json"))
        .then(|| last_index(matches, "json"))
        .flatten();
    let format = (!declares("format"))
        .then(|| {
            let position = last_index(matches, "format")?;
            let value = matches.get_many::<String>("format")?.next_back()?.clone();
            Some((position, value))
        })
        .flatten();
    match spec.output {
        OutputKind::None => {}
        OutputKind::Formats { .. } => {
            // The handler parses its own `--format` and `--json`.
            if let Some(position) = json {
                items.push((position, vec!["--json".to_string()]));
            }
            if command_line(matches, "format") && !declares("format") {
                items.extend(option_occurrences(matches, "format", "format", 1));
            }
        }
        OutputKind::JsonSwitch | OutputKind::JsonOnly => {
            let wants_json = match (json, format) {
                (Some(json_position), Some((format_position, value))) => {
                    if json_position > format_position {
                        Some(json_position)
                    } else {
                        (value == "json").then_some(format_position)
                    }
                }
                (Some(position), None) => Some(position),
                (None, Some((position, value))) => (value == "json").then_some(position),
                (None, None) => None,
            };
            if let Some(position) = wants_json {
                items.push((position, vec!["--json".to_string()]));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Errors.

/// The command path (`oeis generate`) an argument list addresses.
fn command_path(args: &[String]) -> String {
    let Some(first) = args.first() else {
        return String::new();
    };
    let Some(spec) = crate::commands::find(first) else {
        return String::new();
    };
    let mut path = spec.name.to_string();
    let mut current = spec;
    for arg in &args[1..] {
        if arg.starts_with('-') {
            continue;
        }
        match current.find_subcommand(arg) {
            Some(sub) => {
                path.push(' ');
                path.push_str(sub.name);
                current = sub;
            }
            None => break,
        }
    }
    path
}

fn spec_for_path(path: &str) -> Option<&'static CommandSpec> {
    let mut words = path.split(' ');
    let mut spec = crate::commands::find(words.next()?)?;
    for word in words {
        spec = spec.find_subcommand(word)?;
    }
    Some(spec)
}

fn report_parse_error(error: clap::Error, args: &[String]) -> ! {
    match error.kind() {
        ErrorKind::DisplayHelp | ErrorKind::DisplayVersion => {
            print!("{}", restore_braces(&error.render().to_string()));
            let _ = std::io::stdout().flush();
            std::process::exit(0);
        }
        _ => {}
    }
    let path = command_path(args);
    let mut message = error.render().to_string();
    let mut tip = String::new();
    if let Some(spec) = spec_for_path(&path) {
        // Suggest the canonical name rather than a deprecated spelling.
        for (flag, _) in spec.all_flags() {
            for old in flag.deprecated {
                message = message.replace(&format!("'--{old}'"), &format!("'--{}'", flag.name));
            }
        }
        if let Some(ContextValue::String(word)) = error.get(ContextKind::InvalidArg) {
            if !word.starts_with('-') {
                let candidate = format!("{}-{word}", spec.name);
                if crate::commands::find(&candidate).is_some() {
                    tip = format!(
                        "  tip: `{word}` is not a subcommand of `{}`; did you mean \
                         `polytool {candidate}`?\n\n",
                        spec.name
                    );
                }
            }
        }
    }
    let pointer = if path.is_empty() {
        "Run `polytool --help` for the list of commands.".to_string()
    } else {
        format!("Run `polytool help {path}` for usage.")
    };
    let message = message.replace(
        "For more information, try '--help'.",
        &format!("{tip}{pointer}"),
    );
    let message = if message.contains(&pointer) {
        message
    } else {
        format!("{}\n{pointer}\n", message.trim_end())
    };
    eprint!("{message}");
    exit_usage();
}

// ---------------------------------------------------------------------------
// Dispatch.

/// Parse `args` (without the program name) and run the command.
pub(crate) fn run(args: &[String]) {
    let mut root = root_command();
    let argv = std::iter::once("polytool".to_string()).chain(args.iter().cloned());
    let matches = match root.try_get_matches_from_mut(argv) {
        Ok(matches) => matches,
        Err(error) => report_parse_error(error, args),
    };
    let (name, sub_matches) = matches.subcommand().expect("a subcommand is required");
    let spec = crate::commands::find(name).expect("clap only returns declared commands");
    if spec.name != args[0] && spec.deprecated_aliases.contains(&args[0].as_str()) {
        warn_deprecated(&args[0], &args[0], &format!("polytool {}", spec.name));
    }
    let tokens = normalize(spec, sub_matches, spec.name);
    let run = spec.run.expect("top-level commands have a handler");
    run(&tokens);
}

// ---------------------------------------------------------------------------
// Help.

fn levenshtein(a: &str, b: &str) -> usize {
    let b: Vec<char> = b.chars().collect();
    let mut previous: Vec<usize> = (0..=b.len()).collect();
    for (i, ca) in a.chars().enumerate() {
        let mut current = vec![i + 1];
        for (j, cb) in b.iter().enumerate() {
            let substitution = previous[j] + usize::from(ca != *cb);
            current.push(substitution.min(previous[j + 1] + 1).min(current[j] + 1));
        }
        previous = current;
    }
    previous[b.len()]
}

/// The closest of `candidates` to `word`, if it is close enough to be a typo.
pub(crate) fn did_you_mean<'a>(
    word: &str,
    candidates: impl IntoIterator<Item = &'a str>,
) -> Option<&'a str> {
    candidates
        .into_iter()
        .map(|candidate| (levenshtein(word, candidate), candidate))
        .filter(|(distance, candidate)| {
            *distance <= 2.max(word.len() / 3) || candidate.starts_with(word)
        })
        .min_by_key(|(distance, _)| *distance)
        .map(|(_, candidate)| candidate)
}

fn help_usage_error(message: &str, suggestion: Option<&str>) -> ! {
    eprintln!("error: {message}");
    if let Some(suggestion) = suggestion {
        eprintln!("\n  tip: did you mean `polytool help {suggestion}`?");
    }
    eprintln!("\nRun `polytool help` for the list of commands and topics.");
    exit_usage();
}

/// `polytool help [<command> [<subcommand>] | <topic> | --json]`.
pub(crate) fn run_help(words: &[String]) {
    match words {
        [] => print_top_level_help(&mut std::io::stdout()),
        [flag] if crate::input::is_help_arg(flag) => print_top_level_help(&mut std::io::stdout()),
        [flag] if flag == "--json" => print_help_json(),
        [first, rest @ ..] => {
            if let Some(topic) = topics::find(first) {
                if !rest.is_empty() {
                    help_usage_error(&format!("the topic `{}` has no subpages", topic.name), None);
                }
                print!("{}", topics::render(topic));
                return;
            }
            let Some(mut spec) = crate::commands::find(first) else {
                let candidates = COMMANDS
                    .iter()
                    .map(|command| command.name)
                    .chain(topics::TOPICS.iter().map(|topic| topic.name));
                help_usage_error(
                    &format!("unknown command or topic `{first}`"),
                    did_you_mean(first, candidates),
                );
            };
            let mut path = vec![spec.name];
            for word in rest {
                match spec.find_subcommand(word) {
                    Some(sub) => {
                        spec = sub;
                        path.push(sub.name);
                    }
                    None => {
                        let candidates = spec.subcommands.iter().map(|sub| sub.name);
                        let suggestion = did_you_mean(word, candidates)
                            .map(|sub| format!("{} {sub}", path.join(" ")));
                        let message = if spec.subcommands.is_empty() {
                            format!("`{}` has no subcommands (got `{word}`)", path.join(" "))
                        } else {
                            format!("unknown subcommand `{word}` of `{}`", path.join(" "))
                        };
                        help_usage_error(&message, suggestion.as_deref());
                    }
                }
            }
            print!("{}", render_long_help(&path));
        }
    }
}

/// The full help page of a command path such as `["oeis", "generate"]`.
pub(crate) fn render_long_help(path: &[&str]) -> String {
    let mut root = root_command();
    root.build();
    let mut command = &mut root;
    for name in path {
        command = command
            .find_subcommand_mut(name)
            .expect("help paths are resolved against the specification");
    }
    restore_braces(&command.render_help().to_string())
}

pub(crate) fn print_top_level_help(out: &mut dyn Write) {
    let mut text = String::new();
    text.push_str(&format!("{}\n", polytool::version::build_version()));
    text.push_str("Exact univariate polynomial toolkit for combinatorial research.\n\n");
    text.push_str("Usage:\n");
    text.push_str("  polytool <command> [options]           Run a command (input on stdin)\n");
    text.push_str("  polytool help <command> [<subcommand>] Full help for a command\n");
    text.push_str("  polytool help <topic>                  Read a topic (see Topics)\n");
    text.push_str("  polytool help --json                   Machine-readable index\n\n");
    text.push_str("Start here:\n");
    for name in START_HERE {
        let spec = crate::commands::find(name).expect("start-here commands exist");
        text.push_str(&format!("  {:<34} {}\n", spec.name, spec.summary));
    }
    let mut group = "";
    for command in COMMANDS {
        if command.group != group {
            group = command.group;
            text.push_str(&format!("\n{group}:\n"));
        }
        text.push_str(&format!("  {:<34} {}\n", command.name, command.summary));
    }
    text.push_str("\nTopics (polytool help <topic>):\n");
    for topic in topics::TOPICS {
        text.push_str(&format!("  {:<34} {}\n", topic.name, topic.title));
    }
    text.push_str("\nOptions:\n");
    text.push_str("  -h, --help     Print help\n");
    text.push_str("  -V, --version  Print version and build commit\n\n");
    text.push_str("Run `polytool help <command>` for command-specific help.\n");
    let _ = out.write_all(text.as_bytes());
}

fn flag_json(flag: &FlagSpec, heading: Option<&str>) -> Value {
    json!({
        "name": flag.name,
        "long": format!("--{}", flag.name),
        "value": flag.value,
        "arity": flag.arity,
        "help": flag.help,
        "default": flag.default,
        "choices": flag.choices,
        "deprecated_aliases": flag
            .deprecated
            .iter()
            .map(|old| format!("--{old}"))
            .collect::<Vec<_>>(),
        "heading": heading.unwrap_or("Options"),
        "hidden": flag.hidden,
    })
}

fn command_json(spec: &CommandSpec) -> Value {
    let mut flags: Vec<Value> = spec
        .all_flags()
        .map(|(flag, heading)| flag_json(flag, heading))
        .collect();
    for arg in universal_flags(spec) {
        let name = arg.get_id().as_str().to_string();
        let value = arg
            .get_value_names()
            .and_then(|names| names.first())
            .map(|name| name.as_str().to_string());
        let choices: Vec<String> = arg
            .get_possible_values()
            .iter()
            .map(|value| value.get_name().to_string())
            .collect();
        flags.push(json!({
            "name": name,
            "long": format!("--{name}"),
            "value": value,
            "arity": usize::from(value.is_some()),
            "help": arg.get_help().map(ToString::to_string),
            "default": Value::Null,
            "choices": choices,
            "deprecated_aliases": Vec::<String>::new(),
            "heading": "Input and output",
            "hidden": false,
        }));
    }
    json!({
        "name": spec.name,
        "aliases": spec.aliases,
        "deprecated_aliases": spec.deprecated_aliases,
        "group": spec.group,
        "summary": spec.summary,
        "usage": spec.usage.iter().map(|line| format!("polytool {line}")).collect::<Vec<_>>(),
        "description": spec.description,
        "input": spec.input.json_name(),
        "output": {
            "formats": spec.output.format_values(),
            "default": spec.output.default_format(),
        },
        "flags": flags,
        "positionals": spec.positionals.iter().map(|positional| json!({
            "name": positional.name,
            "help": positional.help,
            "required": positional.required,
            "multiple": positional.multiple,
            "choices": positional.choices,
        })).collect::<Vec<_>>(),
        "subcommands": spec.subcommands.iter().map(command_json).collect::<Vec<_>>(),
        "details": spec.details,
        "examples": spec.examples,
        "exit_codes": (if spec.exit_codes.is_empty() { spec::EXIT_ANSWERED } else { spec.exit_codes })
            .iter()
            .map(|(code, meaning)| json!({"code": code, "meaning": meaning}))
            .collect::<Vec<_>>(),
        "see_also": spec.see_also,
    })
}

/// `polytool help --json`: the machine-readable command, flag and topic index.
pub(crate) fn print_help_json() {
    let index = json!({
        "schema": "polytool.commands.v1",
        "commands": COMMANDS.iter().map(command_json).collect::<Vec<_>>(),
        "topics": topics::TOPICS.iter().map(|topic| json!({
            "name": topic.name,
            "title": topic.title,
            "body": topic.body,
        })).collect::<Vec<_>>(),
        "environment": [{
            "name": DEPRECATION_WARNINGS_VARIABLE,
            "help": "Set to 0 to silence deprecation warnings for old spellings",
        }],
    });
    println!(
        "{}",
        serde_json::to_string_pretty(&index).expect("JSON values serialize")
    );
}

// ---------------------------------------------------------------------------
// Completion.

/// `polytool completions <bash|zsh|fish>`.
pub(crate) fn print_completions(args: &[String]) {
    let shell = match args {
        [shell] if shell == "bash" => clap_complete::Shell::Bash,
        [shell] if shell == "zsh" => clap_complete::Shell::Zsh,
        [shell] if shell == "fish" => clap_complete::Shell::Fish,
        _ => crate::input::usage_error("completions", "expected bash, zsh or fish"),
    };
    let mut names: Vec<&'static str> = COMMANDS.iter().map(|command| command.name).collect();
    names.extend(topics::TOPICS.iter().map(|topic| topic.name));
    let help = Command::new("help")
        .about("Print help for a command or topic")
        .arg(
            Arg::new("topic")
                .value_name("COMMAND_OR_TOPIC")
                .num_args(0..)
                .value_parser(PossibleValuesParser::new(names)),
        )
        .arg(Arg::new("json").long("json").action(ArgAction::SetTrue));
    let mut root = root_command().subcommand(help).arg(
        Arg::new("version")
            .short('V')
            .long("version")
            .action(ArgAction::SetTrue)
            .help("Print version and build commit"),
    );
    let mut script = Vec::new();
    clap_complete::generate(shell, &mut root, "polytool", &mut script);
    print!("{}", restore_braces(&String::from_utf8_lossy(&script)));
}
