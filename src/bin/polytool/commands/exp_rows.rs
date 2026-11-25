//! Specification of exact exponential-family row generation and checking.
use crate::cli::spec::{option, switch, CommandSpec, OutputKind, COMMAND};
pub(crate) const EXP_ROWS:CommandSpec=CommandSpec {
    name:"exp-rows",group:"Sequences and recurrences",
    summary:"Generate exact labelled-block exponential rows and check roots in process",
    description:"Generate F_n(t)=n![z^n]exp(A(z)+t B(z)) using factorial-normalized block weights. Dense ranges stream the insertion recurrence; sparse rows use scalar coefficient extraction and a closed form for a single marked size. Root checks reuse the exact PRS engine. Symbolic routes distinguish theorem-backed seeds, paper obstructions, under-review k=8 and finite-rank certificates. Finite success is not an all-rank proof.",
    output:OutputKind::JsonSwitch,
    flags:&[
        option("unmarked","SPEC","none, sizes:1,2, weights:1,1/2, tail:2, exact expression, or JSON block spec").default("none"),
        option("marked","SPEC","Marked block specification in the same syntax"),
        option("rows","INDICES","Increasing ranks: 0..=80 or 304,305,1000"),
        option("check","PROPERTY","Finite check").choices(&["none","real-rooted","interlacing","first-failure"]).default("none"),
        option("input","FILE","Read shared request JSON (or saved report with --verify); '-' reads stdin"),
        switch("summary","Omit passing coefficient rows, retaining failure certificates"),
        switch("verify","Replay complete saved JSON report from --input"),
    ],
    details:"JSON schema polytool.exp-rows.v1. Weights start at size one; expressions use ordinary powers and coeffs:[...] includes the constant slot. Functions: exp, log, inv, compose, trunc, derivative. All rational numbers are exact strings. Zero rows and constants are real-rooted. Interlacing records both orientations for adjacent selected ranks only; either orientation passes. A failure includes integer coefficients, square-free part and distinct PRS root count. First-failure stops at the first selected failure. JSON requests can include finite_rank reference rows with optional pencil factorizations. See doc/cli-exp-rows.md.",
    examples:&["polytool exp-rows --unmarked sizes:1,2 --marked sizes:7 --rows 304,305,1000 --check real-rooted --json","polytool exp-rows --marked tail:1 --rows 0..=20 --check interlacing","polytool exp-rows --input request.json --json","polytool exp-rows --verify --input report.json --json"],
    exit_codes:&[(0,"generation, passing finite checks, or successful replay"),(1,"a requested finite check failed"),(2,"invalid input or replay mismatch")],
    see_also:&["sequence","family-check","recurrence-generate"],
    run:Some(crate::exp_rows::cmd_exp_rows),..COMMAND
};
