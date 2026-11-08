//! CLI for polytool: real-rootedness, interlacing, recurrence search.
//!
//! Reads polynomials from stdin as comma-separated integer coefficients
//! in ascending degree order, one polynomial per line.
//!
//! Each subcommand (or closely related group) lives in its own module;
//! `input`, `json_out` and `help` hold the shared input parsing, JSON
//! fragments and help text.

mod bench;
mod bench_real_rootedness;
mod bicyclic_sieving;
mod bkw_scout;
mod brenti_recognize;
mod coefficient_tests;
mod commands;
mod cyclic_sieving;
mod cyclic_sieving_sequence;
mod ehrhart;
mod enumerate_real_rooted;
mod expansion;
mod family_check;
mod help;
mod input;
mod interlacing;
mod json_out;
mod left_production;
mod oeis;
mod pf_network;
mod pf_pencil;
mod properties;
mod recurrence;
mod recurrence_degree;
mod recurrence_implicit;
mod recurrence_input;
mod recurrence_prove;
mod recurrence_search_representations;
mod recurrence_theorems;
mod recurrence_transform;
mod recurrence_veronese_lifts;
mod resultant;
mod root_enclosures;
mod sequence;
mod stapledon;
mod tnn_network;

use crate::bench::cmd_bench;
use crate::bicyclic_sieving::cmd_bicyclic_sieving;
use crate::bkw_scout::cmd_bkw_scout;
use crate::brenti_recognize::cmd_brenti_recognize;
use crate::coefficient_tests::cmd_coefficient_tests;
use crate::cyclic_sieving::cmd_cyclic_sieving;
use crate::cyclic_sieving_sequence::cmd_cyclic_sieving_sequence;
use crate::ehrhart::{cmd_ehrhart_to_hstar, cmd_hstar_inequalities, cmd_hstar_to_ehrhart};
use crate::enumerate_real_rooted::cmd_enumerate_real_rooted;
use crate::expansion::{cmd_bernstein_expansion, cmd_gamma_expansion};
use crate::family_check::cmd_family_check;
use crate::help::{print_command_help, print_top_level_help};
use crate::input::{is_help_arg, usage_error};
use crate::interlacing::{cmd_interlacing, cmd_interlacing_profile};
use crate::left_production::cmd_left_production;
use crate::oeis::cmd_oeis;
use crate::pf_network::cmd_pf_network;
use crate::pf_pencil::cmd_pf_pencil;
use crate::properties::{cmd_properties, cmd_real_rooted};
use crate::recurrence::{cmd_recurrence, cmd_recurrence_generate};
use crate::recurrence_implicit::cmd_recurrence_implicit;
use crate::recurrence_search_representations::cmd_recurrence_search_representations;
use crate::recurrence_theorems::cmd_recurrence_theorems;
use crate::recurrence_transform::cmd_recurrence_transform;
use crate::recurrence_veronese_lifts::cmd_recurrence_veronese_lifts;
use crate::resultant::{cmd_discriminant, cmd_resultant};
use crate::root_enclosures::cmd_root_enclosures;
use crate::sequence::cmd_sequence;
use crate::stapledon::cmd_stapledon;
use crate::tnn_network::cmd_tnn_network;

/// Exit quietly when stdout is closed early (e.g. `polytool ... | head`).
///
/// Rust ignores `SIGPIPE`, so `println!` panics with "failed printing to
/// stdout: Broken pipe". Treat exactly that panic as the conventional
/// SIGPIPE exit status 141 and keep the default report for every other panic.
fn exit_quietly_on_broken_pipe() {
    let default_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        let message = info
            .payload()
            .downcast_ref::<String>()
            .map(String::as_str)
            .or_else(|| info.payload().downcast_ref::<&str>().copied())
            .unwrap_or("");
        if message.contains("failed printing to stdout") && message.contains("Broken pipe") {
            std::process::exit(141);
        }
        default_hook(info);
    }));
}

fn print_help_path(path: &[String]) -> bool {
    if path.first().is_some_and(|name| name == "oeis") {
        return crate::oeis::print_oeis_subcommand_help(path.get(1).map(String::as_str));
    }
    path.first().is_some_and(|name| print_command_help(name))
}

fn main() {
    exit_quietly_on_broken_pipe();
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        print_top_level_help();
        std::process::exit(1);
    }

    let cmd = &args[1];
    let rest = &args[2..];
    if matches!(cmd.as_str(), "-V" | "--version") {
        if !rest.is_empty() {
            usage_error("", "--version does not accept additional arguments");
        }
        println!("{}", polytool::version::build_version());
        return;
    }
    if cmd == "help" && rest == ["--json"] {
        crate::commands::print_metadata();
        return;
    }
    if is_help_arg(cmd) || cmd == "help" {
        let valid_path = if cmd == "help" {
            rest.len() <= 1 || (rest[0] == "oeis" && rest.len() == 2)
        } else {
            rest.is_empty()
        };
        if !valid_path {
            usage_error(
                "",
                "expected `polytool help <command>` or `polytool help oeis <subcommand>`",
            );
        }
        if cmd == "help" && !rest.is_empty() {
            if is_help_arg(&rest[0]) {
                print_top_level_help();
            } else if !print_help_path(rest) {
                usage_error("", &format!("unknown help topic: {}", rest.join(" ")));
            }
        } else {
            print_top_level_help();
        }
        return;
    }
    if rest.iter().any(|arg| is_help_arg(arg)) {
        if !print_help_path(&args[1..]) {
            usage_error("", &format!("unknown help topic: {}", args[1..].join(" ")));
        }
        return;
    }

    if let Some(command) = crate::commands::find(cmd) {
        (command.run)(rest);
    } else {
        eprintln!("Unknown command: {cmd}. Run `polytool --help` for available commands.");
        std::process::exit(2);
    }
}
