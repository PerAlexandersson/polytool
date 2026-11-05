//! Command metadata shared by dispatch, grouped help, completion and documentation.

use serde::Serialize;

#[derive(Serialize)]
pub(crate) struct CommandSpec {
    pub name: &'static str,
    pub aliases: &'static [&'static str],
    pub group: &'static str,
    pub summary: &'static str,
    #[serde(skip)]
    pub run: fn(&[String]),
}

pub(crate) const COMMANDS: &[CommandSpec] = &[
    CommandSpec {
        name: "real-rooted",
        aliases: &[],
        group: "Inspect polynomials",
        summary: "Check real-rootedness of each polynomial",
        run: crate::cmd_real_rooted,
    },
    CommandSpec {
        name: "interlacing",
        aliases: &[],
        group: "Inspect polynomials",
        summary: "Check interlacing of consecutive polynomial pairs",
        run: crate::cmd_interlacing,
    },
    CommandSpec {
        name: "interlacing-profile",
        aliases: &[],
        group: "Inspect polynomials",
        summary: "Count consecutive previous interlacings until first fail",
        run: crate::cmd_interlacing_profile,
    },
    CommandSpec {
        name: "properties",
        aliases: &[],
        group: "Inspect polynomials",
        summary: "Show real-rootedness, unimodality, and related properties",
        run: crate::cmd_properties,
    },
    CommandSpec {
        name: "family-check",
        aliases: &[],
        group: "Inspect polynomials",
        summary: "Check properties and consecutive interlacing together",
        run: crate::cmd_family_check,
    },
    CommandSpec {
        name: "root-enclosures",
        aliases: &[],
        group: "Inspect polynomials",
        summary: "Exact rational root enclosures with multiplicities",
        run: crate::cmd_root_enclosures,
    },
    CommandSpec {
        name: "resultant",
        aliases: &[],
        group: "Inspect polynomials",
        summary: "Compute the resultant of two polynomials",
        run: crate::cmd_resultant,
    },
    CommandSpec {
        name: "discriminant",
        aliases: &[],
        group: "Inspect polynomials",
        summary: "Compute the discriminant of each polynomial",
        run: crate::cmd_discriminant,
    },
    CommandSpec {
        name: "coefficient-tests",
        aliases: &[],
        group: "Inspect polynomials",
        summary: "Check Newton and Kurtz coefficient criteria",
        run: crate::cmd_coefficient_tests,
    },
    CommandSpec {
        name: "enumerate-real-rooted",
        aliases: &[],
        group: "Inspect polynomials",
        summary: "Enumerate real-rooted integer polynomials by leading prefix",
        run: crate::cmd_enumerate_real_rooted,
    },
    CommandSpec {
        name: "sequence",
        aliases: &[],
        group: "Generate and find recurrences",
        summary: "Generate standard polynomial sequences",
        run: crate::cmd_sequence,
    },
    CommandSpec {
        name: "oeis",
        aliases: &[],
        group: "Generate and find recurrences",
        summary: "List and generate recurrence-backed OEIS families",
        run: crate::cmd_oeis,
    },
    CommandSpec {
        name: "recurrence",
        aliases: &[],
        group: "Generate and find recurrences",
        summary: "Search for a polynomial recurrence",
        run: crate::cmd_recurrence,
    },
    CommandSpec {
        name: "recurrence-generate",
        aliases: &[],
        group: "Generate and find recurrences",
        summary: "Generate coefficient rows from recurrence JSON",
        run: crate::cmd_recurrence_generate,
    },
    CommandSpec {
        name: "recurrence-transform",
        aliases: &[],
        group: "Generate and find recurrences",
        summary: "Change polynomial coordinates and transport a recurrence",
        run: crate::cmd_recurrence_transform,
    },
    CommandSpec {
        name: "recurrence-prove",
        aliases: &[],
        group: "Generate and find recurrences",
        summary: "Discover and replay general recurrence proofs",
        run: crate::recurrence_prove::cmd_recurrence_prove,
    },
    CommandSpec {
        name: "recurrence-theorems",
        aliases: &[],
        group: "Generate and find recurrences",
        summary: "Conditional interlacing theorem matches for recurrence JSON",
        run: crate::cmd_recurrence_theorems,
    },
    CommandSpec {
        name: "recurrence-veronese-lifts",
        aliases: &[],
        group: "Generate and find recurrences",
        summary: "Parent sequences whose coefficient sections are the sequence",
        run: crate::cmd_recurrence_veronese_lifts,
    },
    CommandSpec {
        name: "recurrence-implicit",
        aliases: &[],
        group: "Generate and find recurrences",
        summary: "Diagonal maps and antecedents: implicit relations, uniqueness",
        run: crate::cmd_recurrence_implicit,
    },
    CommandSpec {
        name: "recurrence-search-representations",
        aliases: &[],
        group: "Generate and find recurrences",
        summary: "Bounded search over equivalent representations and relations",
        run: crate::cmd_recurrence_search_representations,
    },
    CommandSpec {
        name: "gamma-expansion",
        aliases: &["gamma"],
        group: "Bases and combinatorics",
        summary: "Expand palindromic polynomials in the gamma basis",
        run: crate::cmd_gamma_expansion,
    },
    CommandSpec {
        name: "bernstein-expansion",
        aliases: &["bernstein"],
        group: "Bases and combinatorics",
        summary: "Expand polynomials in a fixed-degree Bernstein basis",
        run: crate::cmd_bernstein_expansion,
    },
    CommandSpec {
        name: "hstar-to-ehrhart",
        aliases: &[],
        group: "Bases and combinatorics",
        summary: "Convert h*-vectors to Ehrhart polynomials",
        run: crate::cmd_hstar_to_ehrhart,
    },
    CommandSpec {
        name: "ehrhart-to-hstar",
        aliases: &[],
        group: "Bases and combinatorics",
        summary: "Convert Ehrhart polynomials to h*-vectors",
        run: crate::cmd_ehrhart_to_hstar,
    },
    CommandSpec {
        name: "hstar-inequalities",
        aliases: &[],
        group: "Bases and combinatorics",
        summary: "Check named Ehrhart h*-vector inequalities",
        run: crate::cmd_hstar_inequalities,
    },
    CommandSpec {
        name: "cyclic-sieving",
        aliases: &[],
        group: "Bases and combinatorics",
        summary: "Check/profile one cyclic sieving order",
        run: crate::cmd_cyclic_sieving,
    },
    CommandSpec {
        name: "cyclic-sieving-sequence",
        aliases: &[],
        group: "Bases and combinatorics",
        summary: "Profile candidate cyclic sieving orders for a sequence",
        run: crate::cmd_cyclic_sieving_sequence,
    },
    CommandSpec {
        name: "bicyclic-sieving",
        aliases: &[],
        group: "Bases and combinatorics",
        summary: "Check bicyclic sieving for C_n1 x C_n2 and f(q,t)",
        run: crate::cmd_bicyclic_sieving,
    },
    CommandSpec {
        name: "stapledon",
        aliases: &[],
        group: "Bases and combinatorics",
        summary: "Compute a Stapledon decomposition",
        run: crate::cmd_stapledon,
    },
    CommandSpec {
        name: "pf-pencil",
        aliases: &["pf-bidiagonal-pencil"],
        group: "Networks and exploration",
        summary: "Check built-in PF/Jensen endpoint pencils",
        run: crate::cmd_pf_pencil,
    },
    CommandSpec {
        name: "bench",
        aliases: &[],
        group: "Networks and exploration",
        summary: "Run built-in timing suites",
        run: crate::cmd_bench,
    },
    CommandSpec {
        name: "bkw-scout",
        aliases: &[],
        group: "Networks and exploration",
        summary: "Scout BKW equal-modulus loci for a recurrence symbol",
        run: crate::cmd_bkw_scout,
    },
    CommandSpec {
        name: "tnn-network",
        aliases: &[],
        group: "Networks and exploration",
        summary: "Check TNN of a coefficient triangle; canonical planar network",
        run: crate::cmd_tnn_network,
    },
    CommandSpec {
        name: "brenti-recognize",
        aliases: &[],
        group: "Networks and exploration",
        summary: "Search for a planar strip digraph reproducing the rows",
        run: crate::cmd_brenti_recognize,
    },
    CommandSpec {
        name: "pf-network",
        aliases: &[],
        group: "Networks and exploration",
        summary: "Search for a height-invariant planar network: PF rows",
        run: crate::cmd_pf_network,
    },
    CommandSpec {
        name: "left-production",
        aliases: &[],
        group: "Networks and exploration",
        summary: "Chen--Fu--Ruan left production matrix: TNN triangle, PF rows",
        run: crate::cmd_left_production,
    },
    CommandSpec {
        name: "completions",
        aliases: &[],
        group: "Shell support",
        summary: "Generate command-name completion for bash, zsh or fish",
        run: print_completions,
    },
];

pub(crate) fn find(name: &str) -> Option<&'static CommandSpec> {
    COMMANDS
        .iter()
        .find(|cmd| cmd.name == name || cmd.aliases.contains(&name))
}

pub(crate) fn print_metadata() {
    println!(
        "{}",
        serde_json::json!({"schema":"polytool.commands.v1", "commands":COMMANDS})
    );
}

fn print_completions(args: &[String]) {
    let names = COMMANDS
        .iter()
        .flat_map(|cmd| std::iter::once(cmd.name).chain(cmd.aliases.iter().copied()))
        .chain(["help", "--help", "--version"])
        .collect::<Vec<_>>()
        .join(" ");
    match args.first().map(String::as_str).filter(|_| args.len() == 1) {
        Some("bash") => println!(
            r#"_polytool() {{
  local words='' reply
  COMPREPLY=()
  if (( COMP_CWORD == 1 )) || {{ [[ ${{COMP_WORDS[1]}} == help ]] && (( COMP_CWORD == 2 )); }}; then
    words='{names}'
  elif {{ [[ ${{COMP_WORDS[1]}} == oeis ]] && (( COMP_CWORD == 2 )); }} ||
       {{ [[ ${{COMP_WORDS[1]}} == help && ${{COMP_WORDS[2]}} == oeis ]] && (( COMP_CWORD == 3 )); }}; then
    words='list info generate add check evidence --help'
  fi
  if [[ -n $words ]]; then
    while IFS= read -r reply; do COMPREPLY+=("$reply"); done < <(compgen -W "$words" -- "${{COMP_WORDS[COMP_CWORD]}}")
  fi
}}
complete -o default -F _polytool polytool"#
        ),
        Some("zsh") => {
            println!("#compdef polytool\n_arguments '1:command:({names})' '*:file:_files'")
        }
        Some("fish") => println!("complete -c polytool -n '__fish_use_subcommand' -a '{names}'"),
        _ => {
            eprintln!("usage: polytool completions <bash|zsh|fish>");
            std::process::exit(2);
        }
    }
}
