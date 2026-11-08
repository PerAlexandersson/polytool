//! Specifications of the "Planar networks" commands.

use crate::cli::spec::{option, switch, CommandSpec, ExitCode, InputKind, OutputKind, COMMAND};

const GROUP: &str = "Planar networks";

const EXIT_CERTIFICATE: &[ExitCode] = &[
    (0, "certified"),
    (1, "impossible within the bounds"),
    (2, "invalid input or options"),
    (3, "inconclusive (a budget was exhausted)"),
];

pub(crate) const TNN_NETWORK: CommandSpec = CommandSpec {
    name: "tnn-network",
    group: GROUP,
    summary: "Check TNN of a coefficient triangle; canonical planar network",
    description: "Decide exactly whether the coefficient triangle R = (r_{n,k}) of the \
input rows is a lower-unitriangular totally nonnegative (TNN) matrix. On success, print \
the weights lambda[n][k] of its canonical planar network Gamma_N(lambda) (lambda[n][k] \
is the weight of the vertical edge (n+1,k) -> (n,k); horizontal edges have weight 1) and \
re-evaluate the network's path matrix against R (reconstruct, evaluate, compare). On \
failure, say why: a shape problem (not lower-unitriangular), a negative entry, or an \
explicit negative minor (rows, columns, exact value; zero-based indices).\n\n\
By default, row n is a monic polynomial P_n of degree n, with r_{n,k} = [t^k] P_n. The \
rows may start at P_0 = 1 or at P_1, in which case P_0 = 1 is prepended. With --matrix, \
each line is a matrix row, either a square lower-triangular matrix (every row as long \
as the number of rows) or triangular rows (row n has n + 1 entries); shape errors then \
name the offending entry.",
    input: InputKind::Polynomials,
    output: OutputKind::JsonSwitch,
    flags: &[
        switch("matrix", "Read matrix rows instead of monic polynomial rows"),
        switch("chain", "Also compute and check the chain polynomials"),
    ],
    details: "--chain computes the chain polynomials p_0 = 1, p_n = t * sum_{k<n} \
r_{n,k} p_k, and checks that each is real-rooted with zeros in [-1, 0] and that p_n \
interlaces p_{n+1}.\n\n\
Guardrails:\n\n\
\x20\x20TNN of R does NOT make the row polynomials real-rooted: [[1,0,0],[1,1,0],[1,1,1]] is\n\
\x20\x20TNN, yet its last row is 1 + t + t^2. Branden--Saud Maia Leite (arXiv:2412.06595,\n\
\x20\x20Thm 3.7): if R is TNN then the chain polynomials are real-rooted in [-1, 0] and\n\
\x20\x20interlace. The --chain report keeps the TNN verdict separate from these polynomial\n\
\x20\x20checks. A finite triangle certifies only the supplied prefix.\n\n\
JSON output has the schema polytool.tnn-network.v1.",
    examples: &[
        "printf '1,1\\n1,2,1\\n1,3,3,1\\n' | polytool tnn-network --chain",
        "printf '1\\n1,1\\n1,1,1\\n' | polytool tnn-network",
        "printf '1,0,0\\n1,1,0\\n2,0,1\\n' | polytool tnn-network --matrix --json",
    ],
    exit_codes: &[
        (0, "R is TNN and the network certificate re-evaluates to R"),
        (
            1,
            "R is not lower-unitriangular TNN (the report says why), or a check the theory guarantees failed",
        ),
        (2, "invalid input or options"),
    ],
    see_also: &["brenti-recognize", "pf-network", "left-production"],
    run: Some(crate::cmd_tnn_network),
    ..COMMAND
};

pub(crate) const BRENTI_RECOGNIZE: CommandSpec = CommandSpec {
    name: "brenti-recognize",
    group: GROUP,
    summary: "Search for a planar strip digraph reproducing the rows",
    description: "Search exactly for a planar strip digraph D (heights 0..H, edges from \
column x to column x + 1 that depend only on the height, monotone targets) whose \
weighted path counts from (0,0) to (n,k) are the coefficients [t^k] P_n of the input \
rows. By default the input is P_1, P_2, ... with P_0 = 1 implicit.\n\n\
The verdict is one of:\n\n\
\x20\x20certified                 a digraph, re-verified against every input row\n\
\x20\x20impossible_within_bounds  an exact proof that no digraph within the bounds\n\
\x20                           reproduces the rows\n\
\x20\x20inconclusive              a search budget ran out",
    input: InputKind::Polynomials,
    output: OutputKind::JsonSwitch,
    flags: &[
        option(
            "first-index",
            "0|1",
            "Index of the first input row (0: it must be P_0 = 1)",
        )
        .default("1")
        .deprecated(&["first-row"]),
        option(
            "max-target-span",
            "K",
            "At most K consecutive target heights per height",
        )
        .default("no bound")
        .deprecated(&["max-span"]),
        option("weights", "KIND", "Nonnegative rational or integer weights")
            .default("rational")
            .choices(&["rational", "integer"]),
        option("max-search-nodes", "N", "Budget on partial supports")
            .default("1000000")
            .heading("Limits"),
        option(
            "max-integer-search-nodes",
            "N",
            "Budget on the integer search per column",
        )
        .default("10000")
        .deprecated(&["max-integer-nodes"])
        .heading("Limits"),
    ],
    details: "Guardrails:\n\n\
\x20\x20A certificate proves that the coefficient matrix of P_0, ..., P_N is totally\n\
\x20\x20nonnegative (Lindstrom--Gessel--Viennot). It does NOT prove that the rows are\n\
\x20\x20real-rooted or PF: [1 + t, 1 + t + t^2] is certified. The transitions here depend\n\
\x20\x20on the height, so the digraph is invariant under shifts in x, not in y; Brenti's\n\
\x20\x20PF theorem needs height-shift invariance. Only the supplied prefix is covered.\n\
\x20\x20For PF certificates of the rows themselves use pf-network or left-production.\n\n\
JSON output has the schema polytool.brenti-recognize.v1.",
    examples: &[
        "printf '1,1\\n1,2,1\\n1,3,3,1\\n' | polytool brenti-recognize",
        "printf '0,1\\n0,1,1\\n0,1,3,1\\n0,1,7,6,1\\n' | polytool brenti-recognize --weights integer",
        "printf '1,1\\n1,4,1\\n1,11,11,1\\n' | polytool brenti-recognize --json",
    ],
    exit_codes: EXIT_CERTIFICATE,
    see_also: &["tnn-network", "pf-network", "left-production", "topic:budgets"],
    run: Some(crate::cmd_brenti_recognize),
    ..COMMAND
};

pub(crate) const PF_NETWORK: CommandSpec = CommandSpec {
    name: "pf-network",
    group: GROUP,
    summary: "Search for a height-invariant planar network: PF rows",
    description: "Search exactly for a planar network D on {0..N} x Z that is invariant \
under height shifts: for each column n it has steps (length, rise, weight) standing \
for the edges (n - length, y) -> (n, y + rise) at every height y, so\n\n\
\x20\x20P_n = sum over steps into n of weight * t^rise * P_(n - length),  P_0 = 1.\n\n\
Edges are straight segments; the network is valid when every positive-weight step is \
primitive (gcd(length, rise) = 1) and no two edges meet away from a common endpoint. \
Then the Toeplitz matrix of every row is a path matrix of fully compatible sources \
(0, i) and sinks (n, j), hence totally nonnegative, so every row is a Polya frequency \
(PF) sequence: real-rooted with only nonpositive zeros. The certificate is re-verified \
exactly: structure, planarity, every row, and the explicitly rebuilt (deg+1) x (deg+1) \
path matrix against the Toeplitz window of the row.\n\n\
The search tries, column by column, every inclusion-minimal set of steps that is planar \
with the edges already chosen and admits nonnegative weights for that row (exact linear \
feasibility). By default the input is P_1, P_2, ... with P_0 = 1 implicit.\n\n\
The verdict is one of:\n\n\
\x20\x20certified                 a network, re-verified against every input row\n\
\x20\x20impossible_within_bounds  an exact proof that no network within the bounds\n\
\x20                           reproduces the rows (a row that is not real-rooted\n\
\x20                           rules out every network)\n\
\x20\x20inconclusive              a search budget ran out",
    input: InputKind::Polynomials,
    output: OutputKind::JsonSwitch,
    flags: &[
        option(
            "first-index",
            "0|1",
            "Index of the first input row (0: it must be P_0 = 1)",
        )
        .default("1")
        .deprecated(&["first-row"]),
        option("max-step-length", "L", "Largest step length")
            .default("3")
            .deprecated(&["max-length"]),
        option("max-rise", "R", "Largest rise").default("3"),
        option("weights", "KIND", "Nonnegative rational or integer weights")
            .default("rational")
            .choices(&["rational", "integer"]),
        option("max-search-nodes", "N", "Budget on column supports")
            .default("1000000")
            .heading("Limits"),
        option(
            "max-integer-search-nodes",
            "N",
            "Budget on the integer search per row",
        )
        .default("10000")
        .deprecated(&["max-integer-nodes"])
        .heading("Limits"),
    ],
    details: "Guardrails:\n\n\
\x20\x20A certificate proves that every supplied row is PF (real-rooted, nonpositive\n\
\x20\x20zeros). It does NOT prove that the coefficient matrix of the rows is TNN, and it\n\
\x20\x20says nothing about rows beyond the input. Conversely, tnn-network and\n\
\x20\x20brenti-recognize certify TNN of the coefficient matrix, not PF rows.\n\
\x20\x20impossible_within_bounds is relative to --max-step-length and --max-rise unless\n\
\x20\x20the obstruction is a row that is not real-rooted. See also left-production for\n\
\x20\x20the Chen--Fu--Ruan certificate of a whole triangle.\n\n\
JSON output has the schema polytool.pf-network.v1.",
    examples: &[
        "printf '1,1\\n1,2,1\\n1,3,3,1\\n' | polytool pf-network",
        "printf '1,1\\n1,4,1\\n1,11,11,1\\n' | polytool pf-network --json",
        "printf '1,1,1\\n' | polytool pf-network",
    ],
    exit_codes: EXIT_CERTIFICATE,
    see_also: &[
        "left-production",
        "tnn-network",
        "brenti-recognize",
        "topic:budgets",
    ],
    run: Some(crate::cmd_pf_network),
    ..COMMAND
};

pub(crate) const LEFT_PRODUCTION: CommandSpec = CommandSpec {
    name: "left-production",
    group: GROUP,
    summary: "Chen--Fu--Ruan left production matrix: TNN triangle, PF rows",
    description: "Read rows P_0, ..., P_N with deg P_n = n (a lower-triangular matrix A \
with nonzero diagonal, a_{n,k} = [t^k] P_n), compute exactly the leading block Q_N of \
its left production matrix Q = A (1 ⊕ A^{-1}), that is\n\n\
\x20\x20a_{n,0} = q_{n,0},  a_{n,k} = sum_{i=k..n} q_{n,i} a_{i-1,k-1},\n\n\
and decide whether Q_N is totally nonnegative (TNN). By Chen--Fu--Ruan \
(arXiv:2512.08369) a TNN Q_N certifies, for the supplied rows: every row is a Polya \
frequency sequence (real-rooted, nonpositive zeros), and A_N and its reversal \
[a_{n,n-k}] are TNN. On success the certificate is re-verified: A_N is regenerated \
from Q_N, Q_N is checked TNN, and for every row n the unified-network window \
M_{n,n}[n..2n; 0..n] is compared with the transposed Toeplitz section of P_n \
(Theorem 4). An explicit planar network for Q_N is printed: Q_N = L D with L the \
canonical network Gamma_N(lambda) and D the diagonal of Q_N. On failure an explicit \
negative minor of Q_N is given when found.\n\n\
By default the input is P_1, P_2, ... with P_0 = 1 prepended.",
    input: InputKind::Polynomials,
    output: OutputKind::JsonSwitch,
    flags: &[option(
        "first-index",
        "0|1",
        "Index of the first input row (1: P_0 = 1 is prepended)",
    )
    .default("1")
    .deprecated(&["first-row"])],
    details: "Guardrails:\n\n\
\x20\x20The criterion is sufficient, not necessary: if Q_N is not TNN nothing follows\n\
\x20\x20about the rows (the Eulerian rows are real-rooted, but q_{4,2} = -5). A\n\
\x20\x20certificate covers only the supplied rows. Covered families include exponential\n\
\x20\x20Riordan arrays R[f', f] with f' PF (Stirling numbers of both kinds, Lah numbers),\n\
\x20\x20r-Whitney numbers, and Brenti's n-recursive triangles with nonnegative\n\
\x20\x20coefficients (Pascal, Delannoy).\n\n\
JSON output has the schema polytool.left-production.v1.",
    examples: &[
        "printf '0,1\\n0,1,1\\n0,1,3,1\\n0,1,7,6,1\\n' | polytool left-production",
        "printf '1,1\\n1,4,1\\n1,11,11,1\\n1,26,66,26,1\\n' | polytool left-production --json",
    ],
    exit_codes: &[
        (
            0,
            "Q_N is TNN and the certificate re-verifies (rows certified PF)",
        ),
        (1, "not certified"),
        (2, "invalid input or options (for example deg P_n != n)"),
    ],
    see_also: &["pf-network", "tnn-network", "brenti-recognize"],
    run: Some(crate::cmd_left_production),
    ..COMMAND
};
