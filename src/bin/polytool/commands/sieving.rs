//! Specifications of the "Cyclic sieving" commands.

use crate::cli::spec::{option, CommandSpec, InputKind, OutputKind, COMMAND, EXIT_ANSWERED};

const GROUP: &str = "Cyclic sieving";

pub(crate) const CYCLIC_SIEVING: CommandSpec = CommandSpec {
    name: "cyclic-sieving",
    group: GROUP,
    summary: "Check/profile one cyclic sieving order",
    usage: &["cyclic-sieving --order <N> [OPTIONS]"],
    description: "Profile or check P(q) at powers of a primitive N-th root of unity. \
Without fixed counts or orbit sizes, only the exact root-of-unity evaluations are \
profiled.\n\n\
Every run also reports the orbit decomposition of the polynomial (exact integer \
arithmetic; Alexandersson-Amini, 'The cone of cyclic sieving phenomena', Discrete \
Math. 342 (2019), arXiv:1804.01447):\n\n\
\x20\x20- reduce f modulo q^n-1 to residues c_0..c_{n-1}; f(zeta^d) is an integer for\n\
\x20   all d iff c_j = c_gcd(j,n) (Desarmenien); violating gcd classes are listed;\n\
\x20\x20- if integral, f = sum_(r|n) a_r h_r mod q^n-1 with h_r = 1+q^r+...+q^(n-r),\n\
\x20   a_r = sum_(e|r) mu(r/e) c_e = number of orbits of size n/r a CSP action needs;\n\
\x20   S_k = k*a_(n/k) = sum_(j|k) mu(k/j) f(zeta^j) elements lie in orbits of size k;\n\
\x20\x20- cone condition: all a_r >= 0 (equivalently all S_k >= 0). For f in N[q] this\n\
\x20   is necessary and sufficient for a C_n-action making (X, C_n, f) a CSP.\n\n\
The hypothesis (natural coefficients) and the cone condition are reported \
separately. With --orbit-sizes the predicted orbit counts are compared with the \
actual ones. JSON adds `orbit_decomposition` and `orbit_sizes_check` next to the \
unchanged `report` of every item; exact integers are strings.",
    input: InputKind::Polynomials,
    output: OutputKind::JsonSwitch,
    flags: &[
        option("order", "N", "Cyclic group order (positive integer)"),
        option(
            "fixed-counts",
            "LIST",
            "Fixed-point counts fix(g^0),...,fix(g^(N-1)), comma-separated",
        ),
        option(
            "orbit-sizes",
            "LIST",
            "Alternative to --fixed-counts: orbit sizes (divisors of N), comma-separated",
        ),
    ],
    examples: &[
        "echo '1, 1, 1' | polytool cyclic-sieving --order 3",
        "echo '1, 1, 1' | polytool cyclic-sieving --order 3 --orbit-sizes 3 --json",
    ],
    exit_codes: EXIT_ANSWERED,
    see_also: &["cyclic-sieving-sequence", "bicyclic-sieving"],
    run: Some(crate::cmd_cyclic_sieving),
    ..COMMAND
};

pub(crate) const CYCLIC_SIEVING_SEQUENCE: CommandSpec = CommandSpec {
    name: "cyclic-sieving-sequence",
    group: GROUP,
    summary: "Profile candidate cyclic sieving orders for a sequence",
    description: "For a polynomial sequence P_n(q), profile candidate cyclic group \
orders n+offset for every input polynomial; by default the offsets are \
-2,-1,0,1,2,3.\n\n\
The optional fixed-counts file has one entry per line, in the form \
`<index> <order>: c_0,...,c_{order-1}`; blank lines and lines starting with # are \
ignored.",
    input: InputKind::Polynomials,
    output: OutputKind::JsonSwitch,
    flags: &[
        option("first-index", "N", "Index of the first input polynomial").default("0"),
        option(
            "order-offsets",
            "LIST",
            "Candidate order offsets, comma-separated",
        )
        .default("-2,-1,0,1,2,3")
        .deprecated(&["offsets"]),
        option(
            "fixed-counts-file",
            "FILE",
            "Lines `<index> <order>: c_0,...,c_{order-1}` of fixed-point counts",
        ),
    ],
    examples: &[
        "polytool sequence eulerian 6 --format rows | polytool cyclic-sieving-sequence",
        "polytool sequence eulerian 6 --format rows | \
polytool cyclic-sieving-sequence --first-index 1 --order-offsets 0,1 --json",
    ],
    exit_codes: EXIT_ANSWERED,
    see_also: &["cyclic-sieving", "topic:indexing"],
    run: Some(crate::cmd_cyclic_sieving_sequence),
    ..COMMAND
};

pub(crate) const BICYCLIC_SIEVING: CommandSpec = CommandSpec {
    name: "bicyclic-sieving",
    group: GROUP,
    summary: "Check bicyclic sieving for C_n1 x C_n2 and f(q,t)",
    usage: &["bicyclic-sieving --order1 <N1> --order2 <N2> [OPTIONS]"],
    description: "Check bicyclic sieving (Barcelo-Reiner-Stanton) of f(q,t) for \
G = C_n1 x C_n2.\n\n\
Input is a JSON coefficient matrix: f[i][j] is the coefficient of q^i t^j, for \
example\n\n\
\x20\x20[[1,0],[0,1]]   represents   1 + q t\n\n\
Entries are integers or decimal strings (use strings beyond 64 bits).\n\n\
The polynomial is reduced modulo (q^n1-1, t^n2-1) and evaluated exactly at all \
points (zeta_n1^d1, zeta_n2^d2); the limit is n1*n2 <= 2500. If every coefficient \
c_(j1,j2) equals c_(gcd(j1,n1),gcd(j2,n2)), the polynomial lies in the product \
basis h_r1(q) h_r2(t) (Alexandersson-Amini, multivariate basis theorem) and its \
coordinates a_(r1,r2) (two-dimensional Moebius inversion) are reported, with \
f(zeta^d1, omega^d2) = sum over (n1/r1)|d1, (n2/r2)|d2 of (n1/r1)(n2/r2) a_(r1,r2).\n\n\
These coordinates are NOT orbit counts and need not be nonnegative: for the \
diagonal action of C_2 x C_2 on two points,
\x20 f = 1 + qt  has coordinates  1, -1, -1, 2.

Integrality at all points is reported separately from product form (they differ in \
general when gcd(n1,n2) > 1). With --fixed-table the command checks \
f(zeta^d1, omega^d2) = fix(g1^d1 g2^d2) for all d1, d2.",
    input: InputKind::Json,
    output: OutputKind::JsonSwitch,
    flags: &[
        option("order1", "N1", "Order of the cyclic group acting through q"),
        option("order2", "N2", "Order of the cyclic group acting through t"),
        option(
            "fixed-table",
            "JSON",
            "N1 x N2 JSON matrix of fix(g1^d1 g2^d2), d1<N1, d2<N2",
        ),
    ],
    examples: &["echo '[[1,0],[0,1]]' | polytool bicyclic-sieving --order1 2 --order2 2"],
    exit_codes: EXIT_ANSWERED,
    see_also: &["cyclic-sieving"],
    run: Some(crate::cmd_bicyclic_sieving),
    ..COMMAND
};
