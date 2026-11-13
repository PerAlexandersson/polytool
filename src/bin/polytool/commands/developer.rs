//! Specifications of the "Developer tools" commands.

use crate::cli::spec::{
    option, positional, switch, CommandSpec, InputKind, OutputKind, COMMAND, EXIT_ANSWERED,
};

const GROUP: &str = "Developer tools";

pub(crate) const BKW: CommandSpec = CommandSpec {
    name: "bkw",
    group: GROUP,
    summary: "Exact BKW certificates, theorem proofs and named conditional claims",
    description: "Certify stationary recurrences with rational polynomial coefficients and supplied seeds. Use bkw certify; the numerical bkw-scout remains evidence only.",
    subcommands: &[BKW_CERTIFY],
    examples: &["polytool bkw certify --recurrence recurrence.json --json"],
    run: Some(crate::bkw_certify::cmd_bkw),
    ..COMMAND
};

const BKW_CERTIFY: CommandSpec = CommandSpec {
    name: "certify",
    summary: "Decide sparse tie pencils, check the circle criterion, and search exact rows",
    description: "Read existing recurrence JSON with seeds. Reduce its generating function before BKW checks. Public standard quadratic product proofs run first. Exact Toeplitz circle nonnegativity gives proved: every row has simple real zeros and consecutive rows strictly interlace. The theorem source is circle criterion for Hessenberg Toeplitz symbols (unpublished manuscript), recorded in circle.provenance. Versioned Schmidt-Spitzer regions for standard Hessenberg tails of bandwidth at most 3 certify limit-spectrum geometry exactly; a real limit spectrum can give a separate conditional all-row claim under hessenberg_real_limit_spectrum, recorded in conditional_on. Outside the region bandwidth range no conditional verdict is made; the circle theorem can still give proved. --verify replays an embedded certificate exactly.",
    input: InputKind::Json,
    output: OutputKind::JsonSwitch,
    flags: &[
        option("recurrence","FILE","Recurrence JSON, theorem report or representation menu; - reads stdin"),
        option("candidate","N","Zero-based representation-menu candidate"),
        option("max-row-index","N","Largest absolute row index searched, inclusive").default("200"),
        option("max-pencil-degree","N","Degree cap for exact symbolic elimination").default("24"),
        option("radius","R","Positive rational circle radius to try first"),
        switch("verify","Replay a certification report from stdin"),
    ],
    examples: &["polytool bkw certify --recurrence recurrence.json --json", "polytool bkw certify --verify --json < certificate.json"],
    exit_codes: &[(0,"exact report or successful replay"),(2,"invalid input or failed replay")],
    see_also: &["bkw-scout","recurrence-prove"],
    ..COMMAND
};

const BENCH_FORMATS: OutputKind = OutputKind::Formats {
    values: &["tsv", "json", "text"],
    default: "tsv",
};

const EXIT_BENCH: &[(i32, &str)] = &[
    (0, "finished"),
    (1, "a fixture, file or check failed"),
    (2, "usage error"),
];

pub(crate) const BENCH: CommandSpec = CommandSpec {
    name: "bench",
    group: GROUP,
    summary: "Run built-in timing suites",
    description: "Run built-in timing suites and print tab-separated results. Choose a \
suite with a subcommand; run `polytool help bench <subcommand>` for its options.",
    subcommands: &[
        BENCH_RECURRENCE_FIXTURES,
        BENCH_COMPARE,
        BENCH_INTERLACING,
        BENCH_REAL_ROOTEDNESS,
    ],
    examples: &[
        "polytool bench interlacing --max-index 10",
        "polytool bench real-rootedness --quick",
    ],
    exit_codes: EXIT_BENCH,
    run: Some(crate::cmd_bench),
    ..COMMAND
};

const BENCH_RECURRENCE_FIXTURES: CommandSpec = CommandSpec {
    name: "recurrence-fixtures",
    aliases: &["recurrence"],
    summary: "Time adaptive recurrence search on fixture rows",
    description: "Time the adaptive recurrence search on the fixture rows listed in the \
manifest.tsv of a fixture directory. The recurrence-search options of each fixture come \
from its manifest entry. JSON output includes adaptive counters and timing buckets.",
    output: BENCH_FORMATS,
    flags: &[
        option("base", "DIR", "Fixture directory containing manifest.tsv"),
        option("only", "TEXT", "Run only fixture slugs containing TEXT"),
        option("repeat", "N", "Repeat each fixture N times").default("1"),
        switch(
            "no-modular-prefilter",
            "Disable the default modular recurrence prefilter",
        ),
        switch("summary", "Append fixture and category summaries to stdout"),
        option("report", "FILE", "Write a Markdown benchmark report"),
    ],
    examples: &[
        "polytool bench recurrence-fixtures --base fixtures --only eulerian --repeat 3",
        "polytool bench recurrence-fixtures --base fixtures --format json > run.json",
    ],
    exit_codes: EXIT_BENCH,
    see_also: &["bench compare", "recurrence"],
    ..COMMAND
};

const BENCH_COMPARE: CommandSpec = CommandSpec {
    name: "compare",
    summary: "Compare two recurrence-fixture JSON benchmark runs",
    description: "Compare two JSON runs written by `bench recurrence-fixtures --format \
json` and show the worst regressions of the second run relative to the first.",
    output: BENCH_FORMATS,
    flags: &[option("top", "N", "Number of worst regressions to show").default("10")],
    positionals: &[
        positional("OLD", "Earlier benchmark run (JSON)"),
        positional("NEW", "Later benchmark run (JSON)"),
    ],
    examples: &["polytool bench compare --top 5 old.json new.json"],
    exit_codes: EXIT_BENCH,
    see_also: &["bench recurrence-fixtures"],
    ..COMMAND
};

const BENCH_INTERLACING: CommandSpec = CommandSpec {
    name: "interlacing",
    summary: "Time consecutive interlacing checks on a sequence",
    description: "Time the interlacing check of each consecutive pair of polynomials of \
a built-in sequence and print one tab-separated line per pair.",
    flags: &[
        option("sequence", "NAME", "Sequence name").default("eulerian"),
        option("max-index", "N", "Largest sequence index to generate")
            .default("20")
            .deprecated(&["max-n"]),
        option("repeat", "N", "Repeat each pair N times").default("3"),
    ],
    examples: &["polytool bench interlacing --sequence eulerian --max-index 12 --repeat 5"],
    exit_codes: EXIT_BENCH,
    see_also: &["interlacing", "sequence"],
    ..COMMAND
};

const BENCH_REAL_ROOTEDNESS: CommandSpec = CommandSpec {
    name: "real-rootedness",
    aliases: &["real-rooted"],
    summary: "Verdict-checked timing suite on degree 20-80 families",
    description: "Run a timing suite over families of polynomials of degree 20-80. Each \
check has an exact expected verdict; any mismatch (or a verdict that differs from the \
baseline) is reported as a failure and exits with status 1.\n\n\
Check kinds: real_rooted, interlacing (strict, f << g), root_count (distinct real \
roots), hstar_identity (Ehrhart <-> h* conversion round trips).",
    output: BENCH_FORMATS,
    flags: &[
        switch("quick", "Small degrees (about 8-14); fast enough for CI"),
        option(
            "degrees",
            "D1,D2,...",
            "Degrees of the scaling families (each 4..=200)",
        )
        .default("20,40,60,80; quick: 8,14"),
        option("only", "TEXT", "Run only cases whose name contains TEXT"),
        option("repeat", "N", "Time each check N times").default("1"),
        option(
            "baseline",
            "FILE",
            "Compare verdicts and timings with a saved JSON run",
        ),
        option(
            "max-slowdown",
            "X",
            "With --baseline, fail when a check is more than X times slower",
        ),
    ],
    details: "--max-slowdown applies to checks whose baseline median is at least 1 ms, \
and requires --baseline. The full default degree list takes a few minutes.",
    examples: &[
        "polytool bench real-rootedness --quick",
        "polytool bench real-rootedness --quick --json > baseline.json",
        "polytool bench real-rootedness --quick --baseline baseline.json --max-slowdown 3",
    ],
    exit_codes: EXIT_BENCH,
    see_also: &["real-rooted", "bench interlacing"],
    ..COMMAND
};

pub(crate) const BKW_SCOUT: CommandSpec = CommandSpec {
    name: "bkw-scout",
    group: GROUP,
    summary: "Scout BKW equal-modulus loci for a recurrence symbol",
    description: "Scout, numerically, the equal-modulus loci of the roots of a \
recurrence symbol on a grid of complex values of x. The result is floating-point \
evidence, not a proof; it comes with a Mathematica Reduce skeleton for exact follow-up.\n\n\
The symbol is given as z-coefficient polynomials in x in ascending z-degree, separated \
by semicolons. For example, F(x,z) = 1 - xz + z^2 is `1; -x; 1`. Without --symbol the \
symbol is read from standard input. All floating-point options must be finite; \
--min-imaginary-part must be nonnegative.",
    input: InputKind::Text,
    output: OutputKind::JsonSwitch,
    flags: &[
        option(
            "symbol",
            "TEXT",
            "Symbol coefficients, for example '1; -x; 1'",
        ),
        option("box", "R0 R1 I0 I1", "Complex x rectangle")
            .values(4)
            .default("-3 3 -3 3"),
        option("grid", "N", "Use an N by N grid").default("61"),
        option("real-grid-points", "N", "Number of real-axis grid samples")
            .deprecated(&["grid-re"]),
        option(
            "imaginary-grid-points",
            "N",
            "Number of imaginary-axis grid samples",
        )
        .deprecated(&["grid-im"]),
        option("top", "N", "Number of candidates to print").default("20"),
        switch("include-real-axis", "Include samples with Im(x) = 0"),
        option(
            "min-imaginary-part",
            "EPS",
            "Minimum |Im(x)| unless the real axis is included",
        )
        .deprecated(&["min-imag"]),
        option("refine-steps", "N", "Local coordinate-refinement steps").default("8"),
        switch("no-refine", "Disable local refinement"),
        option("tolerance", "EPS", "Positive Durand--Kerner root tolerance").deprecated(&["tol"]),
        option(
            "max-iterations",
            "N",
            "Positive Durand--Kerner iteration limit",
        )
        .deprecated(&["max-iter"]),
        switch(
            "mathematica",
            "Also print an exact Mathematica Reduce skeleton",
        ),
    ],
    details: "With --json, --mathematica adds a `mathematica` string field.",
    examples: &[
        "polytool bkw-scout --symbol '1; -x; 1' --top 5",
        "echo '1; -x; 1' | polytool bkw-scout --box -2 2 -2 2 --grid 41 --json",
    ],
    exit_codes: EXIT_ANSWERED,
    see_also: &["recurrence-theorems"],
    run: Some(crate::cmd_bkw_scout),
    ..COMMAND
};

pub(crate) const PF_PENCIL: CommandSpec = CommandSpec {
    name: "pf-pencil",
    aliases: &["pf-bidiagonal-pencil"],
    group: GROUP,
    summary: "Check built-in PF/Jensen endpoint pencils",
    usage: &[
        "pf-pencil --case <OEIS-ID> --degree <D> [OPTIONS]",
        "pf-pencil --case <OEIS-ID> --min-degree <A> --max-degree <B> [OPTIONS]",
        "pf-pencil --all-family-h --max-degree <B> [OPTIONS]",
    ],
    description: "Check built-in Family H PF-bidiagonal/Jensen endpoint pencils. Choose \
exactly one of --case and --all-family-h, and either --degree or --max-degree \
(optionally with --min-degree); --degree cannot be combined with the range options.",
    output: OutputKind::JsonSwitch,
    flags: &[
        option("case", "OEIS-ID", "Built-in case such as A036969"),
        option("degree", "D", "Check one old-row degree"),
        option("min-degree", "D", "First degree of a range").default("1"),
        option("max-degree", "D", "Last degree of a range"),
        switch("all-family-h", "Check all built-in Family H cases"),
        option(
            "lambdas",
            "LIST",
            "Comma-separated nonnegative integer lambdas",
        )
        .default("0,1,10,100"),
    ],
    details: "Supported cases:\n\n\
\x20\x20A036969 A071951 A080248 A156289 A160562 A269945\n\
\x20\x20A166960 A166961 A166962 A166972 A191935\n\
\x20\x20A371081 A371259 A390433 A198204",
    examples: &[
        "polytool pf-pencil --case A036969 --degree 4",
        "polytool pf-pencil --all-family-h --max-degree 6 --json",
    ],
    exit_codes: EXIT_ANSWERED,
    see_also: &["oeis"],
    run: Some(crate::cmd_pf_pencil),
    ..COMMAND
};
