//! Specifications of the "Check polynomials" and "Roots and algebra" commands.

use crate::cli::spec::{option, switch, CommandSpec, InputKind, OutputKind, COMMAND};

const GROUP_CHECK: &str = "Check polynomials";
const GROUP_ROOTS: &str = "Roots and algebra";

pub(crate) const REAL_ROOTED: CommandSpec = CommandSpec {
    name: "real-rooted",
    group: GROUP_CHECK,
    summary: "Check real-rootedness of each polynomial",
    description: "Decide exactly whether each input polynomial has only real roots. \
Coefficients are arbitrary-size integers.\n\n\
With --json, each item has a zero-based `index` (position in the input), `degree` \
and `real_rooted`. The zero polynomial has degree null; constants have degree 0.",
    input: InputKind::Polynomials,
    output: OutputKind::JsonSwitch,
    examples: &[
        "echo '1, 3, 2' | polytool real-rooted",
        "polytool sequence eulerian 6 --format rows | polytool real-rooted --json",
    ],
    see_also: &[
        "properties",
        "family-check",
        "root-enclosures",
        "topic:input",
    ],
    run: Some(crate::cmd_real_rooted),
    ..COMMAND
};

pub(crate) const INTERLACING: CommandSpec = CommandSpec {
    name: "interlacing",
    group: GROUP_CHECK,
    summary: "Check interlacing of consecutive polynomial pairs",
    description: "Check strict and weak interlacing for each consecutive pair of input \
polynomials (rows 0-1, 1-2, ...). Coefficients are arbitrary-size integers.",
    input: InputKind::Polynomials,
    output: OutputKind::JsonSwitch,
    examples: &["printf '2,-3,1\\n-1,1\\n' | polytool interlacing"],
    see_also: &["interlacing-profile", "family-check"],
    run: Some(crate::cmd_interlacing),
    ..COMMAND
};

pub(crate) const INTERLACING_PROFILE: CommandSpec = CommandSpec {
    name: "interlacing-profile",
    group: GROUP_CHECK,
    summary: "Count consecutive previous interlacings until first fail",
    description: "For each row, count how many consecutive previous rows it interlaces \
with, walking backwards until the first failure. With --json the report includes \
the checked previous-pair reports. Coefficients are arbitrary-size integers.",
    input: InputKind::Polynomials,
    output: OutputKind::JsonSwitch,
    examples: &["printf '1\\n1,1\\n1,2,1\\n' | polytool interlacing-profile"],
    see_also: &["interlacing", "family-check"],
    run: Some(crate::cmd_interlacing_profile),
    ..COMMAND
};

pub(crate) const PROPERTIES: CommandSpec = CommandSpec {
    name: "properties",
    group: GROUP_CHECK,
    summary: "Show real-rootedness, unimodality, and related properties",
    description: "Report real-rootedness, palindromicity, gamma-positivity, unimodality, \
log-concavity and ultra-log-concavity of each input polynomial. Coefficients are \
arbitrary-size integers.",
    input: InputKind::Polynomials,
    output: OutputKind::JsonSwitch,
    examples: &["echo '1, 11, 11, 1' | polytool properties"],
    see_also: &["real-rooted", "family-check", "coefficient-tests"],
    run: Some(crate::cmd_properties),
    ..COMMAND
};

pub(crate) const FAMILY_CHECK: CommandSpec = CommandSpec {
    name: "family-check",
    group: GROUP_CHECK,
    summary: "Check properties and consecutive interlacing together",
    description: "Check a polynomial family in one pass: the properties of every row, \
weak and strict interlacing of consecutive rows and, with --find-recurrence, an \
adaptive recurrence search. The --require-* options turn a property into a \
requirement: the command then exits 1 when some row fails it. All checks accept \
arbitrary-size integer coefficients.",
    input: InputKind::Polynomials,
    output: OutputKind::JsonSwitch,
    flags: &[
        switch("find-recurrence", "Also search for an adaptive recurrence")
            .deprecated(&["recurrence"]),
        switch(
            "require-real-rooted",
            "Exit 1 if any row is not real-rooted",
        )
        .heading("Requirements"),
        switch(
            "require-palindromic",
            "Exit 1 if any row is not palindromic",
        )
        .heading("Requirements"),
        switch(
            "require-gamma-positive",
            "Exit 1 if any row is not gamma-positive",
        )
        .heading("Requirements"),
        switch("require-unimodal", "Exit 1 if any row is not unimodal").heading("Requirements"),
        switch(
            "require-log-concave",
            "Exit 1 if any row is not log-concave",
        )
        .heading("Requirements"),
        switch(
            "require-ultra-log-concave",
            "Exit 1 if any row is not ultra-log-concave",
        )
        .heading("Requirements"),
        switch(
            "require-weak-interlacing",
            "Exit 1 if some consecutive rows do not weakly interlace",
        )
        .heading("Requirements"),
    ],
    examples: &[
        "polytool sequence eulerian 8 --format rows | polytool family-check",
        "polytool family-check --require-real-rooted --require-weak-interlacing --input rows.txt",
    ],
    exit_codes: &[
        (0, "checked; every --require-* condition holds"),
        (1, "some row fails a --require-* condition"),
        (2, "usage error or unreadable input"),
    ],
    see_also: &["properties", "interlacing", "recurrence"],
    run: Some(crate::cmd_family_check),
    ..COMMAND
};

pub(crate) const COEFFICIENT_TESTS: CommandSpec = CommandSpec {
    name: "coefficient-tests",
    group: GROUP_CHECK,
    summary: "Check Newton and Kurtz coefficient criteria",
    description: "Check the Newton inequalities (necessary for real-rootedness of a \
polynomial with positive coefficients) and Kurtz's sufficient condition. A passing \
necessary test is evidence, not a proof.",
    input: InputKind::Polynomials,
    output: OutputKind::JsonSwitch,
    examples: &["echo '1, 4, 1' | polytool coefficient-tests --json"],
    see_also: &["properties", "real-rooted"],
    run: Some(crate::cmd_coefficient_tests),
    ..COMMAND
};

pub(crate) const ENUMERATE_REAL_ROOTED: CommandSpec = CommandSpec {
    name: "enumerate-real-rooted",
    group: GROUP_CHECK,
    summary: "Enumerate real-rooted integer polynomials by leading prefix",
    usage: &[
        "enumerate-real-rooted --degree <N> --leading-coefficients-descending <A0,A1,...> [OPTIONS]",
    ],
    description: "Enumerate every real-rooted integer polynomial
  p(x) = a_0 x^n + a_1 x^(n-1) + ... + a_n
of degree n whose first t leading coefficients a_0, ..., a_(t-1) are given, with
3 <= t <= n + 1 and a_0 != 0 (either sign). Reducible polynomials and polynomials
with repeated roots are included. The search is the Greaves--Syatriadi
derivative-tower algorithm (arXiv:2504.09241) with exact, certified endpoints: no
floating point is used. Nothing is read from stdin.

Coefficient order: --leading-coefficients-descending is DESCENDING (a_0 first);
each output polynomial is ASCENDING, as in all other polytool commands
(c_0, c_1, ..., c_n with c_k = a_(n-k)), so output lines can be piped into
`polytool real-rooted`. Results are sorted lexicographically by (a_t, ..., a_n).",
    output: OutputKind::Formats {
        values: &["text", "json"],
        default: "text",
    },
    flags: &[
        option("degree", "N", "Degree n >= 2 (required)"),
        option(
            "leading-coefficients-descending",
            "LIST",
            "Comma- or space-separated a_0, ..., a_(t-1), leading first (required)",
        )
        .deprecated(&["leading-desc"]),
        option("max-results", "K|none", "Stop after K results")
            .default("10000")
            .heading("Limits"),
        option("max-nodes", "K|none", "Stop after examining K search nodes")
            .default("1000000")
            .heading("Limits"),
        option(
            "max-coefficient-bits",
            "B|none",
            "Stop if a derivative-tower coefficient needs more than B bits",
        )
        .default("none")
        .heading("Limits"),
    ],
    details: "Status values:
  complete        the list is the full set (it may be empty)
  resource_limit  a budget was hit; the list is an initial segment of the full
                  sorted set, and limit_reached names the budget
  invalid_input   (exit status 2) the arguments were rejected
  internal_error  (exit status 1) an internal certification check failed
Errors are reported on stderr; with JSON output a status object is also printed
on stdout. JSON coefficients are strings, so they stay exact beyond 2^53.",
    examples: &[
        "polytool enumerate-real-rooted --degree 6 --leading-coefficients-descending 1,0,-15 --json",
    ],
    exit_codes: &[
        (0, "complete, or stopped by a budget (status resource_limit)"),
        (1, "internal certification error"),
        (2, "usage error or invalid arguments"),
    ],
    see_also: &["real-rooted", "topic:coefficient-order", "topic:budgets"],
    run: Some(crate::cmd_enumerate_real_rooted),
    ..COMMAND
};

pub(crate) const ROOT_ENCLOSURES: CommandSpec = CommandSpec {
    name: "root-enclosures",
    group: GROUP_ROOTS,
    summary: "Exact rational root enclosures with multiplicities",
    usage: &[
        "root-enclosures [--compact] < request.json",
        "root-enclosures --verify [--compact] < report.json",
    ],
    description: "Isolate the real roots of each polynomial exactly. Every distinct real \
root gets a rational enclosure [lower, upper] and its multiplicity: lower == upper is \
an exact rational root, otherwise exactly one root lies strictly inside and neither \
endpoint is a root. Midpoints are representatives only, never certified bounds. \
Reads one JSON request and prints a polytool.root-enclosures.v1 report.",
    input: InputKind::Json,
    output: OutputKind::JsonOnly,
    flags: &[
        switch(
            "verify",
            "Replay a report: recompute its embedded request and compare",
        ),
        switch("compact", "Print compact (single-line) JSON").heading("Input and output"),
    ],
    details: "Request fields:
  polynomials   [[\"c0\", \"c1\", ...], ...] ascending rational strings
  max_width     optional rational: refine enclosures to at most this width
                (default refinement: width <= 1/1024)
  interval      optional {lower?, upper?} with endpoints
                {value, per_index?, closed?}; an omitted end is infinite and
                value + per_index * n uses n = first_index + position
  first_index   index of the first polynomial (default 0)

Each item reports degree, valuation (multiplicity of 0), real_rooted,
real_root_count (with multiplicity), the enclosures and, with an interval,
real_roots_in_interval and real_rooted_in_interval (decided exactly). The zero
polynomial has no degree and satisfies interval tests vacuously.",
    examples: &["echo '{\"polynomials\":[[\"0\",\"0\",\"1\",\"3\",\"3\",\"1\"]],
\"interval\":{\"lower\":{\"value\":\"-1\"},\"upper\":{\"value\":\"0\"}}}' | polytool root-enclosures"],
    exit_codes: &[
        (0, "success (or a valid replay with --verify)"),
        (2, "invalid input, or a replay that does not match"),
    ],
    see_also: &["real-rooted", "topic:json-schemas"],
    run: Some(crate::cmd_root_enclosures),
    ..COMMAND
};

pub(crate) const RESULTANT: CommandSpec = CommandSpec {
    name: "resultant",
    group: GROUP_ROOTS,
    summary: "Compute the resultant of two polynomials",
    description: "Compute the resultant of exactly two input polynomials. Coefficients \
are arbitrary-size integers; with --json exact values are decimal strings.",
    input: InputKind::Polynomials,
    output: OutputKind::JsonSwitch,
    examples: &["printf '1,0,1\\n-1,1\\n' | polytool resultant"],
    see_also: &["discriminant"],
    run: Some(crate::cmd_resultant),
    ..COMMAND
};

pub(crate) const DISCRIMINANT: CommandSpec = CommandSpec {
    name: "discriminant",
    group: GROUP_ROOTS,
    summary: "Compute the discriminant of each polynomial",
    description: "Compute the discriminant of each input polynomial. Coefficients are \
arbitrary-size integers; with --json exact values are decimal strings.",
    input: InputKind::Polynomials,
    output: OutputKind::JsonSwitch,
    examples: &["echo '1, 0, 1' | polytool discriminant"],
    see_also: &["resultant"],
    run: Some(crate::cmd_discriminant),
    ..COMMAND
};
