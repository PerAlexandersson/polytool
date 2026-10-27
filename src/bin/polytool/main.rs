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
mod recurrence_implicit;
mod recurrence_search_representations;
mod recurrence_theorems;
mod recurrence_transform;
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
use crate::input::is_help_arg;
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
        println!("{}", polytool::version::build_version());
        return;
    }
    if is_help_arg(cmd) {
        if cmd == "help" && !rest.is_empty() {
            if is_help_arg(&rest[0]) {
                print_top_level_help();
            } else if !print_command_help(&rest[0]) {
                eprintln!("Unknown command: {}", rest[0]);
                std::process::exit(1);
            }
        } else {
            print_top_level_help();
        }
        return;
    }
    if rest.iter().any(|arg| is_help_arg(arg)) {
        if !print_command_help(cmd) {
            eprintln!("Unknown command: {}", cmd);
            std::process::exit(1);
        }
        return;
    }

    match cmd.as_str() {
        "real-rooted" => cmd_real_rooted(),
        "interlacing" => cmd_interlacing(rest),
        "interlacing-profile" => cmd_interlacing_profile(rest),
        "properties" => cmd_properties(rest),
        "gamma-expansion" | "gamma" => cmd_gamma_expansion(rest),
        "bernstein-expansion" | "bernstein" => cmd_bernstein_expansion(rest),
        "family-check" => cmd_family_check(rest),
        "sequence" => cmd_sequence(rest),
        "oeis" => cmd_oeis(rest),
        "pf-pencil" | "pf-bidiagonal-pencil" => cmd_pf_pencil(rest),
        "recurrence" => cmd_recurrence(rest),
        "recurrence-generate" => cmd_recurrence_generate(rest),
        "recurrence-transform" => cmd_recurrence_transform(rest),
        "recurrence-theorems" => cmd_recurrence_theorems(rest),
        "root-enclosures" => cmd_root_enclosures(rest),
        "recurrence-implicit" => cmd_recurrence_implicit(rest),
        "recurrence-search-representations" => cmd_recurrence_search_representations(rest),
        "bench" => cmd_bench(rest),
        "bkw-scout" => cmd_bkw_scout(rest),
        "resultant" => cmd_resultant(),
        "discriminant" => cmd_discriminant(),
        "coefficient-tests" => cmd_coefficient_tests(rest),
        "hstar-to-ehrhart" => cmd_hstar_to_ehrhart(rest),
        "ehrhart-to-hstar" => cmd_ehrhart_to_hstar(rest),
        "hstar-inequalities" => cmd_hstar_inequalities(rest),
        "cyclic-sieving" => cmd_cyclic_sieving(rest),
        "cyclic-sieving-sequence" => cmd_cyclic_sieving_sequence(rest),
        "bicyclic-sieving" => cmd_bicyclic_sieving(rest),
        "stapledon" => cmd_stapledon(rest),
        "tnn-network" => cmd_tnn_network(rest),
        "brenti-recognize" => cmd_brenti_recognize(rest),
        "enumerate-real-rooted" => cmd_enumerate_real_rooted(rest),
        "pf-network" => cmd_pf_network(rest),
        "left-production" => cmd_left_production(rest),
        _ => {
            eprintln!("Unknown command: {}", cmd);
            std::process::exit(1);
        }
    }
}
