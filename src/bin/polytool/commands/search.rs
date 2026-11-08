//! The recurrence-search flags shared by `recurrence`, `oeis add` and
//! `bench recurrence-fixtures` (one parser per command, one specification).

use crate::cli::spec::{option, switch, FlagGroup};

/// Shape of the recurrences to try and how rows are split into fit and
/// held-out rows.
pub(crate) const RECURRENCE_SEARCH: FlagGroup = FlagGroup {
    heading: "Recurrence search",
    flags: &[
        option("skip-prefix", "N", "Ignore the first N input polynomials").default("0"),
        option(
            "min-recurrence-length",
            "N",
            "Minimum number of previous rows a recurrence uses",
        )
        .default("1")
        .deprecated(&["min-rec-len"]),
        option(
            "max-recurrence-length",
            "N",
            "Maximum number of previous rows a recurrence uses",
        )
        .default("5")
        .deprecated(&["max-rec-len"]),
        option(
            "min-variable-degree",
            "D",
            "Minimum degree in t of the recurrence coefficients",
        )
        .default("0")
        .deprecated(&["min-var-deg"]),
        option(
            "max-variable-degree",
            "D",
            "Maximum degree in t of the recurrence coefficients",
        )
        .default("3")
        .deprecated(&["max-var-deg"]),
        option(
            "min-index-degree",
            "D",
            "Minimum degree in the index n of the recurrence coefficients",
        )
        .default("0")
        .deprecated(&["min-idx-deg"]),
        option(
            "max-index-degree",
            "D",
            "Maximum degree in the index n of the recurrence coefficients",
        )
        .default("3")
        .deprecated(&["max-idx-deg"]),
        option(
            "min-derivative-order",
            "D",
            "Minimum order of the derivatives P'(n-k), P''(n-k), ... used",
        )
        .default("0")
        .deprecated(&["min-diff-deg"]),
        option(
            "max-derivative-order",
            "D",
            "Maximum order of the derivatives used",
        )
        .default("2")
        .deprecated(&["max-diff-deg"]),
        switch("try-inhomogeneous", "Also try an inhomogeneous term")
            .deprecated(&["inhomogeneous"]),
        option(
            "min-inhomogeneous-variable-degree",
            "D",
            "Minimum degree in t of the inhomogeneous term",
        )
        .default("0")
        .deprecated(&["min-inhomo-var-deg"]),
        option(
            "max-inhomogeneous-variable-degree",
            "D",
            "Maximum degree in t of the inhomogeneous term",
        )
        .default("3")
        .deprecated(&["max-inhomo-var-deg"]),
        option(
            "min-inhomogeneous-index-degree",
            "D",
            "Minimum degree in n of the inhomogeneous term",
        )
        .default("0")
        .deprecated(&["min-inhomo-idx-deg"]),
        option(
            "max-inhomogeneous-index-degree",
            "D",
            "Maximum degree in n of the inhomogeneous term",
        )
        .default("3")
        .deprecated(&["max-inhomo-idx-deg"]),
        switch(
            "try-denominator",
            "Also try a left-hand-side denominator (rational recurrences)",
        )
        .deprecated(&["denominator"]),
        option(
            "max-denominator-variable-degree",
            "D",
            "Maximum degree in t of the denominator",
        )
        .default("2")
        .deprecated(&["max-denom-var-deg"]),
        option(
            "max-denominator-index-degree",
            "D",
            "Maximum degree in n of the denominator",
        )
        .default("2")
        .deprecated(&["max-denom-idx-deg"]),
        switch("try-alternating-sign", "Also try alternating-sign terms")
            .deprecated(&["alternating-sign"]),
        option(
            "min-equation-margin",
            "N",
            "Require at least N more equations than unknowns",
        )
        .default("1")
        .deprecated(&["min-margin"]),
        option(
            "fit-extra-rows",
            "N",
            "Fit on N rows beyond the first solvable prefix",
        )
        .default("1"),
        switch(
            "no-held-out-rows",
            "Fit on every input row (no held-out verification rows)",
        )
        .deprecated(&["no-verify"]),
        switch(
            "no-modular-prefilter",
            "Disable the default modular rejection of candidates",
        ),
    ],
};
