//! Top-level help, shared input descriptions and the `help <command>`
//! dispatcher.  Command-specific help lives with each command.

use crate::bench::print_bench_help;
use crate::bicyclic_sieving::print_bicyclic_sieving_help;
use crate::bkw_scout::bkw_scout_usage;
use crate::brenti_recognize::print_brenti_recognize_help;
use crate::coefficient_tests::print_coefficient_tests_help;
use crate::cyclic_sieving::print_cyclic_sieving_help;
use crate::cyclic_sieving_sequence::print_cyclic_sieving_sequence_help;
use crate::ehrhart::print_hstar_inequalities_help;
use crate::enumerate_real_rooted::print_enumerate_real_rooted_help;
use crate::family_check::print_family_check_help;
use crate::left_production::print_left_production_help;
use crate::oeis::print_oeis_help;
use crate::pf_network::print_pf_network_help;
use crate::pf_pencil::print_pf_pencil_help;
use crate::recurrence::{print_recurrence_generate_help, print_recurrence_help};
use crate::recurrence_theorems::print_recurrence_theorems_help;
use crate::recurrence_transform::print_recurrence_transform_help;
use crate::sequence::print_sequence_help;
use crate::stapledon::print_stapledon_help;
use crate::tnn_network::print_tnn_network_help;

pub(crate) fn print_coefficient_input_help() {
    println!("Input:");
    println!("  Read one polynomial per line from stdin.");
    println!("  Dense coefficient lists are in ascending degree order:");
    println!("    a_0, a_1, ..., a_d  represents  a_0 + a_1 t + ... + a_d t^d");
    println!("  Brackets, whitespace-separated coefficients, and expanded notation");
    println!("  such as 1 + 2t + 3t^2 are also accepted.");
    println!("  Blank lines and lines starting with # are ignored.");
}

pub(crate) fn print_rational_coefficient_input_help() {
    println!("Input:");
    println!("  Read dense coefficient lists from stdin, one polynomial per line.");
    println!("  Coefficients are in ascending degree order and may be arbitrary-size");
    println!("  integers or exact rationals such as 42, -17, or 3/7.");
    println!("  Brackets and whitespace-separated coefficients are accepted.");
    println!("  Expanded polynomial expressions are not accepted by recurrence.");
    println!("  Blank lines and lines starting with # are ignored.");
}

pub(crate) fn print_top_level_help() {
    println!("{}", polytool::version::build_version());
    println!("Dense univariate polytool for combinatorial research.");
    println!();
    println!("Usage:");
    println!("  polytool <command> [options]");
    println!("  polytool --help");
    println!("  polytool help <command>");
    println!();
    println!("Commands:");
    println!("  real-rooted       Check real-rootedness of each polynomial");
    println!("  interlacing       Check interlacing of consecutive polynomial pairs");
    println!("  interlacing-profile");
    println!("                    Count consecutive previous interlacings until first fail");
    println!("  properties        Show real-rootedness, unimodality, and related properties");
    println!("  gamma-expansion   Expand palindromic polynomials in the gamma basis");
    println!("  bernstein-expansion");
    println!("                    Expand polynomials in a fixed-degree Bernstein basis");
    println!("  family-check      Check properties and consecutive interlacing together");
    println!("  sequence          Generate standard polynomial sequences");
    println!("  oeis              List and generate recurrence-backed OEIS families");
    println!("  pf-pencil         Check built-in PF/Jensen endpoint pencils");
    println!("  recurrence        Search for a polynomial recurrence");
    println!("  recurrence-generate");
    println!("                    Generate coefficient rows from recurrence JSON");
    println!("  recurrence-transform");
    println!("                    Change polynomial coordinates and transport a recurrence");
    println!("  recurrence-theorems");
    println!("                    Conditional interlacing theorem matches for recurrence JSON");
    println!("  bench             Run built-in timing suites");
    println!("  bkw-scout         Scout BKW equal-modulus loci for a recurrence symbol");
    println!("  resultant         Compute the resultant of two polynomials");
    println!("  discriminant      Compute the discriminant of each polynomial");
    println!("  coefficient-tests Check Newton and Kurtz coefficient criteria");
    println!("  hstar-to-ehrhart  Convert h*-vectors to Ehrhart polynomials");
    println!("  ehrhart-to-hstar  Convert Ehrhart polynomials to h*-vectors");
    println!("  hstar-inequalities");
    println!("                    Check named Ehrhart h*-vector inequalities");
    println!("  cyclic-sieving    Check/profile one cyclic sieving order");
    println!("  cyclic-sieving-sequence");
    println!("                    Profile candidate cyclic sieving orders for a sequence");
    println!("  bicyclic-sieving  Check bicyclic sieving for C_n1 x C_n2 and f(q,t)");
    println!("  stapledon         Compute a Stapledon decomposition");
    println!("  tnn-network       Check TNN of a coefficient triangle; canonical planar network");
    println!("  brenti-recognize  Search for a planar strip digraph reproducing the rows");
    println!("  enumerate-real-rooted");
    println!("                    Enumerate real-rooted integer polynomials by leading prefix");
    println!("  pf-network        Search for a height-invariant planar network: PF rows");
    println!("  left-production   Chen--Fu--Ruan left production matrix: TNN triangle, PF rows");
    println!();
    println!("Options:");
    println!("  -h, --help        Print help text");
    println!("  -V, --version     Print version and build commit");
    println!();
    println!("Run `polytool help <command>` for command-specific help.");
}

fn print_stdin_command_help(command: &str, summary: &str, extra: &[&str]) {
    println!("Usage:");
    println!("  polytool {command}");
    println!("  polytool {command} --help");
    println!();
    println!("{summary}");
    println!();
    print_coefficient_input_help();
    if !extra.is_empty() {
        println!();
        for line in extra {
            println!("{line}");
        }
    }
}

pub(crate) fn print_command_help(command: &str) -> bool {
    match command {
        "real-rooted" => print_stdin_command_help(
            "real-rooted",
            "Check whether each input polynomial has only real roots.",
            &[
                "Note:",
                "  Accepts arbitrary-size integer coefficients.",
                "",
                "Example:",
                "  echo '1, 3, 2' | polytool real-rooted",
            ],
        ),
        "interlacing" => print_stdin_command_help(
            "interlacing",
            "Check strict and weak interlacing for consecutive input pairs.",
            &[
                "Options:",
                "  --json    Emit machine-readable JSON",
                "",
                "Note:",
                "  Accepts arbitrary-size integer coefficients.",
                "",
                "Example:",
                "  printf '2,-3,1\\n-1,1\\n' | polytool interlacing",
            ],
        ),
        "interlacing-profile" => print_stdin_command_help(
            "interlacing-profile",
            "For each row, count backward consecutive previous interlacings until first failure.",
            &[
                "Options:",
                "  --json    Emit machine-readable JSON, including checked previous-pair reports",
                "",
                "Note:",
                "  Accepts arbitrary-size integer coefficients.",
                "",
                "Example:",
                "  printf '1\\n1,1\\n1,2,1\\n' | polytool interlacing-profile",
            ],
        ),
        "properties" => print_stdin_command_help(
            "properties",
            concat!(
                "Report real-rootedness, palindromicity, gamma-positivity, ",
                "unimodality, log-concavity, and ultra-log-concavity."
            ),
            &[
                "Options:",
                "  --json    Emit machine-readable JSON",
                "",
                "Note:",
                "  Accepts arbitrary-size integer coefficients.",
                "",
                "Example:",
                "  echo '1, 11, 11, 1' | polytool properties",
            ],
        ),
        "gamma-expansion" | "gamma" => print_stdin_command_help(
            "gamma-expansion",
            concat!(
                "Expand each palindromic input polynomial in the gamma basis ",
                "t^i (1+t)^(d-2i)."
            ),
            &[
                "Options:",
                "  --json    Emit machine-readable JSON",
                "",
                "Note:",
                "  Accepts arbitrary-size integer coefficients.",
                "",
                "Aliases:",
                "  polytool gamma",
                "",
                "Example:",
                "  echo '1, 11, 11, 1' | polytool gamma-expansion",
            ],
        ),
        "bernstein-expansion" | "bernstein" => {
            println!("Usage:");
            println!("  polytool bernstein-expansion [--degree <n>] [--json]");
            println!("  polytool bernstein --degree <n>");
            println!();
            println!("Expand each input polynomial in the standard Bernstein basis");
            println!("  binomial(n,j) t^j (1-t)^(n-j),  0 <= j <= n.");
            println!();
            print_rational_coefficient_input_help();
            println!();
            println!("Options:");
            println!("  --degree <n>  Use ambient degree n for every row");
            println!("                (default: the degree of each input polynomial)");
            println!("  --json        Emit machine-readable JSON with rational strings");
            println!();
            println!("Example:");
            println!("  echo '1, 2, 3' | polytool bernstein --degree 3");
        }
        "family-check" => print_family_check_help(),
        "sequence" => print_sequence_help(),
        "oeis" => print_oeis_help(),
        "pf-pencil" | "pf-bidiagonal-pencil" => print_pf_pencil_help(),
        "recurrence" => print_recurrence_help(),
        "recurrence-generate" => print_recurrence_generate_help(),
        "recurrence-transform" => print_recurrence_transform_help(),
        "recurrence-theorems" => print_recurrence_theorems_help(),
        "bench" => print_bench_help(),
        "bkw-scout" => bkw_scout_usage(),
        "resultant" => print_stdin_command_help(
            "resultant",
            "Compute the resultant of the first two input polynomials.",
            &[
                "Note:",
                "  Accepts arbitrary-size integer coefficients.",
                "",
                "Example:",
                "  printf '1,0,1\\n-1,1\\n' | polytool resultant",
            ],
        ),
        "discriminant" => print_stdin_command_help(
            "discriminant",
            "Compute the discriminant of each input polynomial.",
            &[
                "Note:",
                "  Accepts arbitrary-size integer coefficients.",
                "",
                "Example:",
                "  echo '1, 0, 1' | polytool discriminant",
            ],
        ),
        "coefficient-tests" => print_coefficient_tests_help(),
        "hstar-to-ehrhart" => print_stdin_command_help(
            "hstar-to-ehrhart",
            "Convert each h*-vector into Ehrhart polynomial coefficients.",
            &[
                "Options:",
                "  --json    Emit machine-readable JSON",
                "",
                "Note:",
                "  Accepts arbitrary-size integer h*-vector entries.",
                "",
                "Example:",
                "  echo '1, 2, 1' | polytool hstar-to-ehrhart",
            ],
        ),
        "ehrhart-to-hstar" => print_stdin_command_help(
            "ehrhart-to-hstar",
            "Convert each Ehrhart polynomial into an h*-vector.",
            &[
                "Input coefficients may be integers or exact rationals.",
                "",
                "Options:",
                "  --json    Emit machine-readable JSON",
                "",
                "Note:",
                "  Returns arbitrary-size integer h*-vector entries.",
                "",
                "Example:",
                "  echo '1, 2, 2' | polytool ehrhart-to-hstar",
            ],
        ),
        "hstar-inequalities" => print_hstar_inequalities_help(),
        "cyclic-sieving" => print_cyclic_sieving_help(),
        "cyclic-sieving-sequence" => print_cyclic_sieving_sequence_help(),
        "bicyclic-sieving" => print_bicyclic_sieving_help(),
        "stapledon" => print_stapledon_help(),
        "tnn-network" => print_tnn_network_help(),
        "brenti-recognize" => print_brenti_recognize_help(),
        "enumerate-real-rooted" => print_enumerate_real_rooted_help(),
        "pf-network" => print_pf_network_help(),
        "left-production" => print_left_production_help(),
        _ => return false,
    }
    true
}
