//! CLI for polytool: real-rootedness, interlacing, recurrence search.
//!
//! Reads polynomials from stdin as comma-separated integer coefficients
//! in ascending degree order, one polynomial per line.
//!
//! Each subcommand (or closely related group) lives in its own module;
//! `input`, `json_out` and `help` hold the shared input parsing, JSON
//! fragments and help text.

mod bench;
mod bicyclic_sieving;
mod bkw_scout;
mod brenti_recognize;
mod coefficient_tests;
mod cyclic_sieving;
mod cyclic_sieving_sequence;
mod ehrhart;
mod expansion;
mod family_check;
mod help;
mod input;
mod interlacing;
mod json_out;
mod oeis;
mod pf_pencil;
mod properties;
mod recurrence;
mod resultant;
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
use crate::expansion::{cmd_bernstein_expansion, cmd_gamma_expansion};
use crate::family_check::cmd_family_check;
use crate::help::{print_command_help, print_top_level_help};
use crate::input::is_help_arg;
use crate::interlacing::{cmd_interlacing, cmd_interlacing_profile};
use crate::oeis::cmd_oeis;
use crate::pf_pencil::cmd_pf_pencil;
use crate::properties::{cmd_properties, cmd_real_rooted};
use crate::recurrence::{cmd_recurrence, cmd_recurrence_generate};
use crate::resultant::{cmd_discriminant, cmd_resultant};
use crate::sequence::cmd_sequence;
use crate::stapledon::cmd_stapledon;
use crate::tnn_network::cmd_tnn_network;

fn main() {
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
        _ => {
            eprintln!("Unknown command: {}", cmd);
            std::process::exit(1);
        }
    }
}
