//! Exact univariate stability command specifications.

use crate::cli::spec::{option, positional, CommandSpec, InputKind, OutputKind, COMMAND};

pub(crate) const STABILITY: CommandSpec = CommandSpec {
    name: "stability",
    group: "Check polynomials",
    summary: "Count roots exactly in an open circular region",
    description: "Read Gaussian-rational univariate expressions or ascending comma-separated \
coefficient lists, from positional arguments or one per line on stdin/--input FILE. \
Use I for the imaginary unit; floating inputs and nonconstant denominators fail. \
Counts include multiplicities and distinguish inside, boundary and outside. Upper/lower \
default to avoidance; left/right/disk/exterior default to all-inside. Zero has null \
degree/counts, passes avoidance by convention and fails all-inside. Constants pass both.",
    input: InputKind::Text,
    output: OutputKind::Formats { values: &["summary", "json", "full"], default: "summary" },
    positionals: &[positional("POLYNOMIAL", "Exact univariate polynomial (repeat for a batch)").optional().multiple()],
    flags: &[
        option("region", "REGION", "Named open region; exclusive with mobius and centered disk").choices(&["upper","lower","left","right","disk","exterior"]).default("upper"),
        option("mobius", "A,B,C,D", "Image of the upper half-plane under (Az+B)/(Cz+D)"),
        option("disk-center", "CENTER", "Exact Gaussian center; requires --radius"),
        option("radius", "R", "Positive rational radius; requires --disk-center"),
        option("convention", "CONVENTION", "Override the region's default verdict meaning").choices(&["avoid","inside"]),
        option("certificate", "FILE", "Write the full replayable polytool.stability.v1 report"),
        option("max-degree", "N", "Degree cap, at most 2000").default("2000"),
        option("max-coefficient-bits", "N", "Numerator/denominator bit cap, at most 1000000").default("1000000"),
    ],
    subcommands: &[VERIFY],
    details: "--format json (or --json) returns compact items; full includes canonical input, \
options and exact certificate entries as strings. --certificate writes the same full report. \
Selectors are exclusive; disk center/radius are required together. Möbius maps must be nonsingular. \
Use -- before a positional expression beginning with a minus sign.",
    examples: &[
        "polytool stability '(1+t)^3' --region left",
        "polytool stability 't^2+(1+2I)t-I' --region upper --json",
        "polytool stability 't-I' --disk-center I --radius 1 --certificate report.json",
        "polytool stability verify --input report.json",
    ],
    run: Some(crate::stability::cmd_stability),
    ..COMMAND
};

const VERIFY: CommandSpec = CommandSpec {
    name: "verify",
    summary: "Replay every exact stability certificate",
    description: "Read a full polytool.stability.v1 report from --input FILE or stdin and \
recompute every field exactly. Altered or malformed reports fail. Uses degree 2000 and \
coefficient-part 1000000-bit caps.",
    input: InputKind::Json,
    output: OutputKind::JsonOnly,
    examples: &["polytool stability verify --input report.json"],
    exit_codes: &[
        (0, "report verified"),
        (1, "report does not verify"),
        (2, "invalid input"),
    ],
    ..COMMAND
};
