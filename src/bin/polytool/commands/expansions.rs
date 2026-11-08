//! Specifications of the "Expansions and conversions" commands.

use crate::cli::spec::{
    option, positional, CommandSpec, InputKind, OutputKind, COMMAND, EXIT_ANSWERED,
};

const GROUP: &str = "Expansions and conversions";

pub(crate) const GAMMA_EXPANSION: CommandSpec = CommandSpec {
    name: "gamma-expansion",
    aliases: &["gamma"],
    group: GROUP,
    summary: "Expand palindromic polynomials in the gamma basis",
    description: "Expand each palindromic input polynomial in the gamma basis \
t^i (1+t)^(d-2i), where d is the degree. A polynomial that is not palindromic has no \
such expansion and is reported as an error. Coefficients are arbitrary-size integers.",
    input: InputKind::Polynomials,
    output: OutputKind::JsonSwitch,
    examples: &[
        "echo '1, 11, 11, 1' | polytool gamma-expansion",
        "polytool sequence eulerian 6 --format rows | polytool gamma --json",
    ],
    exit_codes: &[
        (0, "every row has a gamma expansion"),
        (1, "some row is not palindromic"),
        (2, "usage error or invalid input"),
    ],
    see_also: &["properties", "topic:coefficient-order"],
    run: Some(crate::cmd_gamma_expansion),
    ..COMMAND
};

pub(crate) const BERNSTEIN_EXPANSION: CommandSpec = CommandSpec {
    name: "bernstein-expansion",
    aliases: &["bernstein"],
    group: GROUP,
    summary: "Expand polynomials in a fixed-degree Bernstein basis",
    description: "Expand each input polynomial in the standard Bernstein basis\n\n\
\x20\x20binomial(n,j) t^j (1-t)^(n-j),  0 <= j <= n.\n\n\
Coefficients may be arbitrary-size integers or exact rationals such as 42, -17 or 3/7. \
The command expects coefficient lists, not expanded expressions. With --json the \
coordinates are rational strings. A row whose degree exceeds the ambient degree is \
reported as an error.",
    input: InputKind::RationalLists,
    output: OutputKind::JsonSwitch,
    flags: &[option(
        "degree",
        "N",
        "Use ambient degree N for every row (default: degree of each row)",
    )],
    examples: &[
        "echo '1, 2, 3' | polytool bernstein-expansion --degree 3",
        "echo '1, 1/2' | polytool bernstein --json",
    ],
    exit_codes: &[
        (0, "every row was expanded"),
        (1, "some row could not be expanded in the requested degree"),
        (2, "usage error or invalid input"),
    ],
    see_also: &["gamma-expansion", "topic:coefficient-order"],
    run: Some(crate::cmd_bernstein_expansion),
    ..COMMAND
};

pub(crate) const HSTAR_TO_EHRHART: CommandSpec = CommandSpec {
    name: "hstar-to-ehrhart",
    group: GROUP,
    summary: "Convert h*-vectors to Ehrhart polynomials",
    description: "Convert each h*-vector into the coefficients of the Ehrhart polynomial \
L(n). Entries are arbitrary-size integers; the output coefficients are exact rationals.",
    input: InputKind::Polynomials,
    output: OutputKind::JsonSwitch,
    examples: &["echo '1, 2, 1' | polytool hstar-to-ehrhart"],
    exit_codes: EXIT_ANSWERED,
    see_also: &["ehrhart-to-hstar", "hstar-inequalities"],
    run: Some(crate::cmd_hstar_to_ehrhart),
    ..COMMAND
};

pub(crate) const EHRHART_TO_HSTAR: CommandSpec = CommandSpec {
    name: "ehrhart-to-hstar",
    group: GROUP,
    summary: "Convert Ehrhart polynomials to h*-vectors",
    description: "Convert each Ehrhart polynomial into its h*-vector. Input coefficients \
may be integers or exact rationals; the h*-vector entries are arbitrary-size integers.",
    input: InputKind::RationalLists,
    output: OutputKind::JsonSwitch,
    examples: &["echo '1, 2, 2' | polytool ehrhart-to-hstar"],
    exit_codes: EXIT_ANSWERED,
    see_also: &["hstar-to-ehrhart", "hstar-inequalities"],
    run: Some(crate::cmd_ehrhart_to_hstar),
    ..COMMAND
};

pub(crate) const HSTAR_INEQUALITIES: CommandSpec = CommandSpec {
    name: "hstar-inequalities",
    group: GROUP,
    summary: "Check named Ehrhart h*-vector inequalities",
    description: "Check named Ehrhart h*-vector inequalities with exact integer \
arithmetic. Each input polynomial is an h*-vector in ascending degree order.",
    input: InputKind::Polynomials,
    output: OutputKind::JsonSwitch,
    flags: &[option(
        "dimension",
        "D",
        "Use dimension D instead of len(h*)-1 (at most 4096)",
    )],
    examples: &["echo '1, 4, 1, 0' | polytool hstar-inequalities --dimension 3"],
    exit_codes: EXIT_ANSWERED,
    see_also: &["hstar-to-ehrhart", "ehrhart-to-hstar"],
    run: Some(crate::cmd_hstar_inequalities),
    ..COMMAND
};

pub(crate) const STAPLEDON: CommandSpec = CommandSpec {
    name: "stapledon",
    group: GROUP,
    summary: "Compute a Stapledon decomposition",
    description: "Compute the Stapledon decomposition of each input polynomial with \
respect to the degree bound N, printed as `f = a + x (b)`. A polynomial whose degree \
exceeds N is reported on standard error and skipped. Coefficients are arbitrary-size \
integers.",
    input: InputKind::Polynomials,
    output: OutputKind::None,
    positionals: &[positional("N", "Degree bound (nonnegative integer)")],
    examples: &["echo '1, 4, 1' | polytool stapledon 3"],
    exit_codes: &[
        (
            0,
            "finished (rows beyond the bound are reported but do not change the status)",
        ),
        (2, "usage error or invalid input"),
    ],
    see_also: &["hstar-inequalities"],
    run: Some(crate::cmd_stapledon),
    ..COMMAND
};
