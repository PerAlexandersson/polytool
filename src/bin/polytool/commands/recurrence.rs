//! Specifications of the "Sequences and recurrences" commands.

use crate::cli::spec::{option, positional, switch, CommandSpec, InputKind, OutputKind, COMMAND};

const GROUP: &str = "Sequences and recurrences";

/// Exit statuses of the commands that prove or replay a statement.
const EXIT_PROVED: &[(i32, &str)] = &[
    (0, "proved (or the replayed certificate is valid)"),
    (1, "not proved"),
    (2, "usage error, invalid input or a failed replay"),
];

/// Exit statuses of `recurrence-prove`.
const EXIT_PROVE: &[(i32, &str)] = &[
    (0, "proved by polytool (or the replayed report is valid)"),
    (
        1,
        "not proved by polytool (also conditional_on_attestation)",
    ),
    (2, "usage error, invalid input or a failed replay"),
];

/// Exit statuses of the JSON-request commands.
const EXIT_REQUEST: &[(i32, &str)] = &[
    (0, "success (or a replay that verified)"),
    (2, "invalid input or a failed replay"),
];

pub(crate) const SEQUENCE: CommandSpec = CommandSpec {
    name: "sequence",
    group: GROUP,
    summary: "Generate standard polynomial sequences",
    description: r#"Generate the first polynomials of a standard family. Generated coefficients use
arbitrary-size integers.

Names:
  eulerian
  narayana
  type-b-eulerian  (aliases: type-b, type_b_eulerian, type_b)
  chebyshev-t      (aliases: chebyshev_t, t)
  chebyshev-u      (aliases: chebyshev_u, u)
  hermite

The Eulerian and Narayana sequences are labelled A_1, ..., A_N and N_1, ..., N_N;
the other families start at index 0. With --format rows the output is one
ascending coefficient list per line, ready for the other commands."#,
    output: OutputKind::Formats {
        values: &["text", "rows", "json"],
        default: "text",
    },
    positionals: &[
        positional("NAME", "Family name (see the list above)"),
        positional("MAX-N", "Number of polynomials (the largest index)"),
    ],
    examples: &[
        "polytool sequence eulerian 5",
        "polytool sequence eulerian 6 --format rows | polytool real-rooted",
    ],
    see_also: &["oeis", "recurrence", "topic:output-formats"],
    run: Some(crate::cmd_sequence),
    ..COMMAND
};

pub(crate) const RECURRENCE: CommandSpec = CommandSpec {
    name: "recurrence",
    group: GROUP,
    summary: "Search for a polynomial recurrence",
    usage: &[
        "recurrence [OPTIONS]",
        "recurrence --space [--json] < request.json",
        "recurrence --verify-space [--json] < report.json",
    ],
    description: r#"Search adaptively for a linear differential-polynomial recurrence among the input
polynomial sequence. The recurrence relates P_n to earlier rows P_(n-k) and
their derivatives with coefficients that are polynomials in t and in the index
n. A recurrence found this way is a finite fit (evidence on the supplied rows),
not a proof for all n; held-out rows are used to check it unless
--no-held-out-rows is given.

Input is read as dense coefficient lists, one polynomial per line, ascending
degree, with arbitrary-size integers or exact rationals such as 42, -17 or 3/7.
This command expects coefficient lists, not expanded expressions.

--first-index labels the first input row (nonnegative); with --skip-prefix the
discarded rows still count. By default the retained rows are renumbered from 1.

Output formats: text (default), json (the recurrence JSON with initial
conditions), latex (every coefficient is factored over Q and the product
verified by exact re-expansion, e.g.
(n + 2) P_{n} = (2n + 1)(1 + t) P_{n-1} - (n - 1)(1 - t)^{2} P_{n-2};
a coefficient whose factorization exhausts its budget stays expanded, noted by
a % comment line) and python (exact Python code). --latex and --python are
shorthands for the matching --format value.

--space and --verify-space are separate modes: they read a JSON request (or a
replayed report) from standard input, emit exact JSON and accept only an
optional --json; bounds belong in the request. Source-backed positives are
conditional on the supplied source and seeds, rows-only positives are fits, and
dual negatives prove absence within the shape. The opt-in request field
"coefficient_rescaling": {"kind":"factorial_power","power":m} searches
Q[n] = T(P[n]) with T(X^k) = X^k/(k!)^m, m >= 1, and certifies complete spaces
from r(U+m(A+J)+1) residuals (constant-index source)."#,
    input: InputKind::RationalLists,
    output: OutputKind::Formats {
        values: &["text", "json", "latex", "python"],
        default: "text",
    },
    flags: &[
        option(
            "first-index",
            "N",
            "Index of the first input row (nonnegative)",
        ),
        option(
            "max-candidates",
            "N",
            "Inspect at most N candidate configurations",
        )
        .heading("Limits"),
        option("max-work", "N", "Stop after N deterministic work units").heading("Limits"),
        switch(
            "no-laurent",
            "Do not vary the left-hand side (x^e, n-factors, charts)",
        ),
        option(
            "max-laurent-power",
            "E",
            "Largest left-hand power x^E (default 2)",
        )
        .deprecated(&["max-laurent-x-power"]),
        switch(
            "representations",
            "Return a menu of length/derivative tradeoffs",
        ),
        switch("theorem-matches", "Also report conditional theorem matches"),
        switch("latex", "Emit LaTeX (same as --format latex)"),
        switch("python", "Emit exact Python code (same as --format python)"),
        switch("verbose", "Print candidate search details"),
        switch(
            "space",
            "Complete bounded spaces: read a space request (JSON)",
        ),
        switch(
            "verify-space",
            "Replay a recurrence-space report (JSON), no search",
        ),
        switch(
            "modular-prefilter",
            "Enable the modular rejection (the default)",
        )
        .hidden(),
    ],
    shared: &[&crate::commands::search::RECURRENCE_SEARCH],
    details: r#"Budgets:
  --max-work counts deterministic work units over the whole search (modular and
  exact solves, held-out checks, precomputation) from the dimensions of the
  systems, never from the clock: the same N gives the same output on every
  machine. A step runs only if its whole cost fits, so the units used never
  exceed N. One solve of r equations in c unknowns costs r*c*(min(r,c)+1) units.
  Reaching --max-work or --max-candidates ends the search as budget_exhausted
  (exit status 3). JSON reports work_units and max_work (a found record has
  search.work_units). --max-work is not supported with --representations.

--representations:
  Text output is a comparison table and numbered formulas; JSON is the full
  menu. Save the menu with --json and select an entry with --candidate N in
  recurrence-generate or recurrence-theorems. It cannot be combined with
  --theorem-matches or --format python.

--theorem-matches:
  Also report conditional theorem matches (weighted renewal, fixed-delay
  predecessor), kept separate from the finite fit evidence. With --json the
  output is polytool.recurrence-theorems.v1 (the usual record is its
  recurrence_json). Without the flag a match is only noted on stderr.

--space / --verify-space (JSON):
  The reply is a space report; the JSON output lists structured candidates and
  coverage. Exit status 3 means some requested shape is not completed."#,
    examples: &[
        "printf '1\\n1\\n2\\n3\\n5\\n8\\n' | polytool recurrence",
        "polytool sequence eulerian 8 --format rows | polytool recurrence --json > recurrence.json",
        "polytool recurrence --input rows.txt --representations",
    ],
    exit_codes: &[
        (
            0,
            "finished (a recurrence was found, or none within the bounds)",
        ),
        (2, "usage error or invalid input"),
        (
            3,
            "the search budget was exhausted before the space was covered",
        ),
    ],
    see_also: &[
        "recurrence-generate",
        "recurrence-theorems",
        "recurrence-prove",
        "topic:budgets",
        "topic:recurrence-json",
        "topic:indexing",
    ],
    run: Some(crate::cmd_recurrence),
    ..COMMAND
};

pub(crate) const RECURRENCE_GENERATE: CommandSpec = CommandSpec {
    name: "recurrence-generate",
    group: GROUP,
    summary: "Generate coefficient rows from recurrence JSON",
    usage: &[
        "recurrence-generate [--recurrence <FILE|->] --count <N> [OPTIONS]",
        "recurrence-generate [--recurrence <FILE|->] --additional-count <N> [OPTIONS]",
        "recurrence-generate [--recurrence <FILE|->] --format bfile --max-terms <N> [OPTIONS]",
    ],
    description: r#"Generate coefficient rows from recurrence JSON emitted by `polytool recurrence
--json`. The input also accepts theorem wrappers and representation menus (pick
a menu entry with --candidate, 0-based). Without --recurrence (or with
`--recurrence -`) the record is read from standard input, or from --input FILE.

Row selection: --count N outputs N rows in total, including the initial rows;
--additional-count N outputs the initial rows plus N generated rows. Use only
one of them. --first-index relabels the first row and shifts the rule
accordingly (default: keep the JSON record's first_index).

Formats: text (the default; `rows` is an alias) and json print the generated
rows. The OEIS text formats read the rows from --from-index on as an OEIS
triangle, with the recurrence index n shifted to the OEIS row number, and show
--count rows (default 10):
  oeis-example   %e block "Triangle begins:" with aligned columns
  oeis-formula   %F sentence: the recurrence in plain ASCII with its initial
                 rows, e.g. P(n,x) = (2+x)*P(n-1,x) + ...
  mathematica    %t program: memoized P[n_, x_] and the flattened rows
  pari           %o PARI/GP program listing the same rows
  oeis-json      all four as JSON (schema polytool.oeis-text.v1)
Without --first-term-index the first row keeps its own index as its OEIS row
number.

B-file output (--format bfile, alias b-file): the rows are flattened into OEIS
b-file lines `index value`, coefficients in ascending powers within each row,
and every coefficient must be an integer. A term cap never splits a row."#,
    input: InputKind::Json,
    output: OutputKind::Formats {
        values: &[
            "text",
            "rows",
            "json",
            "bfile",
            "b-file",
            "oeis-example",
            "oeis-formula",
            "mathematica",
            "pari",
            "oeis-json",
        ],
        default: "text",
    },
    flags: &[
        option(
            "recurrence",
            "FILE|-",
            "Read recurrence JSON from FILE (default -: stdin or --input FILE)",
        ),
        option("candidate", "N", "Select a representation menu entry (0-based)"),
        option(
            "first-index",
            "N",
            "Relabel the first row and shift the rule accordingly",
        ),
        option("count", "N", "Output N rows in total, including initial rows")
            .deprecated(&["rows"]),
        option("additional-count", "N", "Output initial rows plus N generated rows")
            .deprecated(&["additional"]),
        option(
            "first-term-index",
            "K",
            "Index of the first term (the OEIS offset; default 0)",
        )
        .heading("B-file and OEIS text output")
        .deprecated(&["offset"]),
        option(
            "from-index",
            "K",
            "Start at the row with absolute index K (default: first index)",
        )
        .heading("B-file and OEIS text output")
        .deprecated(&["first-row"]),
        switch("descending", "Write each row from its highest coefficient down")
            .heading("B-file and OEIS text output")
            .deprecated(&["reversed"]),
        option(
            "trim-zeros",
            "WHERE",
            "Drop zeros before --descending is applied (default none)",
        )
        .choices(&["none", "low", "high", "both"])
        .heading("B-file and OEIS text output"),
        switch(
            "regular-triangle",
            "Row i has exactly L0 + i entries: pad with high zeros",
        )
        .heading("B-file and OEIS text output"),
        option("max-terms", "N", "Write at most N terms, whole rows only")
            .heading("B-file output"),
        switch(
            "nearest-row",
            "With --max-terms, stop at the whole row nearest to N terms",
        )
        .heading("B-file output"),
    ],
    details: r#"B-file and OEIS text options:
  --trim-zeros     none keeps everything; low drops zeros at the constant end,
                   high at the leading end, both at both ends of each row.
  --regular-triangle  Row i (counted from the first) has exactly L0 + i entries,
                   where L0 is the length of the first row: shorter rows are
                   padded with high zeros and generation stops before a longer
                   row.
  --first-term-index, --from-index, --descending, --trim-zeros and
                   --regular-triangle need --format bfile or an OEIS text format.
  --max-terms and --nearest-row need --format bfile. Without --count or
                   --additional-count, rows are generated until the cap is
                   reached. --nearest-row stops at the whole row nearest to N
                   terms (it may exceed N; ties go to fewer rows) instead of the
                   longest prefix with at most N terms."#,
    examples: &[
        "polytool recurrence --input rows.txt --json > recurrence.json",
        "polytool recurrence-generate --recurrence recurrence.json --count 50",
        "polytool recurrence-generate --recurrence eulerian.json --format bfile \\\n    --first-term-index 1 --max-terms 1000 > b008292.txt",
        "polytool recurrence-generate --recurrence rec.json --format mathematica \\\n    --from-index 1 --first-term-index 0 --descending --count 11",
    ],
    see_also: &[
        "recurrence",
        "oeis",
        "topic:recurrence-json",
        "topic:output-formats",
        "topic:indexing",
    ],
    run: Some(crate::cmd_recurrence_generate),
    ..COMMAND
};

pub(crate) const RECURRENCE_TRANSFORM: CommandSpec = CommandSpec {
    name: "recurrence-transform",
    group: GROUP,
    summary: "Change polynomial coordinates and transport a recurrence",
    usage: &[
        "recurrence-transform [--compact] < request.json",
        "recurrence-transform --verify < report.json",
    ],
    description: r#"Apply an exact change of polynomial coordinates to a sequence, transport its
defining recurrence into the new coordinates, and optionally search
independently there. Reads one JSON request from standard input and prints a
polytool.recurrence-transform.v1 report. Bounds belong in the request, not on
the command line."#,
    input: InputKind::Json,
    output: OutputKind::JsonOnly,
    flags: &[
        switch("verify", "Replay a report exactly from its embedded request"),
        switch("compact", "Print compact JSON"),
    ],
    details: r#"Request fields:
  transformation   {"kind": ...} with kind one of
                     multiply_power {power}      Q_n = x^r P_n
                     divide_power {power}        Q_n = P_n / x^r (needs x^r | P_n)
                     reversal {schedule}         Q_n = x^(d_n) P_n(1/x)
                     affine {a, b}               Q_n = P_n(a x + b), a != 0
                     mobius {a, b, c, e, schedule}
                                     Q_n = (c x + e)^(d_n) P_n((a x + b)/(c x + e))
                     f_polynomial {schedule}     Q_n = (1 + x)^(d_n) P_n(x/(1 + x))
                     gamma {schedule}            gamma coordinates at symmetry index d_n
                                                 (rules transported when the even and odd
                                                 parts of the identity are proportional)
                                                 (otherwise transported.system: the even
                                                 part and its odd companion, two relations)
                     symmetric_decomposition {schedule}
                                     (a_n, b_n) with P_n = a_n + x b_n
                   a, b, c, e are exact rationals as strings ("3/2").
  schedule         {first_index, first_degree, step, last_index?}:
                   d_n = first_degree + step*(n - first_index) on
                   first_index <= n (<= last_index); decreasing needs last_index
  source_recurrence  defining polytool.recurrence.v1 rule with seeds, or
  oeis_id          a bundled OEIS family (e.g. "A008292"), or
  rows             explicit rows (ascending rational strings) with row_origin
  row_count        rows generated from the source (default 10)
  search           {} or bounds {max_rec_len, max_var_deg, max_idx_deg,
                   max_diff_deg, try_denominator, max_candidates, ...}:
                   also search independently in the new coordinates
  limits           {max_rows, max_degree, ...} resource limits
  root_domain      optional {} or {candidate_intervals, assumptions,
                   default_window, max_candidate_intervals, max_observed_rows,
                   max_observed_degree}: attach root-domain information
                   (never needed for the algebraic results). An interval is
                   {lower?, upper?} with endpoints {value, per_index?, closed?};
                   an assumption is {subject: source|transformed, component?,
                   property: real_rooted|real_rooted_in|real_roots_in, interval?,
                   first_index?, last_index?, note?}

Report:
  transformation      exact row map, inverse, image conditions, domain
  transformed_states  [row][component][coefficient]
  reconstruction_verified  inverse(map(P_n)) == P_n for every row
  conditions          certified_all_indices (from the source rule) or
                      checked_supplied_rows (finite evidence)
  transported         evidence transported_identity: exact consequence of
                      the source rule, with transformed seeds
  transported.system  rank-two gamma case: [even part, odd companion] as
                      implicit relations of the gamma rows, exact consequences
                      where the rows are symmetric (see conditions)
  independent_search  evidence finite_fit, kept separately
  root_domain         implications (direction, hypotheses with status,
                      converse), claims labelled finite_exact_observation,
                      supplied_assumption, proposed_interval,
                      completed_recurrence_induction or
                      certified_family_invariant, the Rolle recurrence match
                      with its interval attempts, and exact root enclosures"#,
    examples: &[
        "echo '{\"oeis_id\":\"A008292\",\"row_count\":8,\"search\":{},\n  \"transformation\":{\"kind\":\"f_polynomial\",\n  \"schedule\":{\"first_index\":1,\"first_degree\":0,\"step\":1}}}' \\\n  | polytool recurrence-transform",
        "polytool recurrence-transform --input request.json > report.json\npolytool recurrence-transform --verify --input report.json",
    ],
    exit_codes: EXIT_REQUEST,
    see_also: &[
        "recurrence-implicit",
        "recurrence-search-representations",
        "topic:json-schemas",
    ],
    run: Some(crate::cmd_recurrence_transform),
    ..COMMAND
};

pub(crate) const RECURRENCE_DEGREE: CommandSpec = CommandSpec {
    name: "recurrence-degree",
    group: GROUP,
    summary: "Degree law and leading-coefficient signs of a recurrence for all n",
    usage: &[
        "recurrence-degree [--recurrence <FILE|->] [--candidate N] [--count N] [--max-period P]",
        "recurrence-degree --verify < certificate.json",
    ],
    description: r#"Decide the degree law deg P_n = D(n) and the sign of the leading coefficient of
every row P_n of a recurrence family (rule, first index and seed rows), for all
n, exactly. D(n) is periodic by residue class (affine laws, half growth and
other periodic laws); it is read off exactly computed rows and then proved for
every index by one of the routes of the RealRooted tactic
rr_row_leadingCoeff_pos: cone (nonnegative top-coefficient multipliers),
recurrence (a growth-ratio invariant rho), subseq (separate residue classes),
twoStep (P_n = q P_(n-2)), cancel (first-order rules whose top terms cancel,
optionally with a ratio bound [x^(D-1)] P_n >= rho lead P_n) and coefficients
(nonnegative coefficients in an affine degree envelope; the degree is then only
bounded). Every statement about all n is decided by the exact integer-range sign
checker; the rows before the induction start are computed exactly. A route that
fails is not a disproof.

The recurrence input is recurrence.v1 JSON, a theorem wrapper or a
representation menu; with a menu, --candidate is required unless it has one
entry. --verify replays a certificate read from standard input and takes no
search flags."#,
    input: InputKind::Json,
    output: OutputKind::JsonSwitch,
    flags: &[
        option(
            "recurrence",
            "FILE|-",
            "Read recurrence JSON from FILE (default -: stdin or --input FILE)",
        ),
        option("candidate", "N", "Explicit 0-based menu choice"),
        option("count", "N", "Rows computed to read off the law")
            .default("36")
            .deprecated(&["rows"]),
        option("max-period", "P", "Largest period of the law").default("6"),
        switch("verify", "Replay a certificate from stdin; no search flags"),
    ],
    details: r#"JSON output is polytool.recurrence-leading-coefficient.v1. With --verify and --json
the reply is {"verified": true, "status": ...}."#,
    examples: &[
        "polytool recurrence --json < rows.txt | polytool recurrence-degree --recurrence -",
        "polytool recurrence-degree --recurrence rec.json --json > certificate.json\npolytool recurrence-degree --verify < certificate.json",
    ],
    exit_codes: EXIT_PROVED,
    see_also: &["recurrence-prove", "recurrence", "topic:recurrence-json"],
    run: Some(crate::recurrence_degree::cmd_recurrence_degree),
    ..COMMAND
};

pub(crate) const RECURRENCE_PROVE: CommandSpec = CommandSpec {
    name: "recurrence-prove",
    group: GROUP,
    summary: "Discover and replay general recurrence proofs",
    usage: &[
        "recurrence-prove [--recurrence <FILE|->] [--candidate N] [--max-start-delay N] \
         [--max-representation-candidates N] [--max-lift-candidates N] \
         [--max-affine-coordinates N]",
        "recurrence-prove --verify < proof.json",
    ],
    description: r#"Discover general proofs for a recurrence and seeds; no named-sequence lookup.
Tries, in the order of the parent cascade, product induction, Favard/three-term,
Liu-Wang, derivative-lag, Ma-Wang interval, Liu-Wang three-term, iterated
preserver and hypergeometric-row criteria, then the RealRooted row routes
fixed-delay window (variable_delay: f(n) P_n = A(n) P_(n-1) + B(n) x P_(n-L)
with a delay L of at least 2) and residue classes (residue_class_subsequences:
every lag a multiple of a period p of at least 2, each class P_(first + j + p m)
proved by a criterion above), the last two recorded only when the rule has their
form; every hypothesis is checked exactly for all indices. Reports every attempt, its obligations and
conclusion. A fit still needs
an all-index identity for an independently defined family. Rust certificates can
be replayed; Lean replay is not attempted.

Iterated preserver also searches bounded exact factorizations of lag-one
second-order steps over Q(n)[x]. Each first-order factor uses the existing
preserver menu or the shared Liu-Wang conditions on its actual input degree.
The certificate records application order and replays the full Ore identity;
intermediate degree and sign conditions hold for every rule index. This can
certify a closed ray even when the composite fails a whole-line preserver
test. No simplicity or interlacing is inferred from the factorization alone.

Index-section lift (recorded only when it applies): the rows, possibly
reversed, are every p-th row (p <= 5) of a constant-coefficient parent found by
factoring, an identity for every index; the parent is proved by the criteria.

Factor gauge (recorded only when it applies): every lag-j multiplier is
kappa_j(n) q(x)^j for one fixed q, so P_n = q^(n-e) g S_n with S_n in the span
of the gauged seed rows; S_n is real-rooted for every n by its degree (<= 1),
a discriminant (degree 2 or 3) or the palindromic-quartic test (u = x + 1/x),
signs of polynomials in n decided exactly.

Positive lag-two Euler (lag_two_euler, recorded only when its form matches):
with theta = x D and m = n - first_index, the operators are
[(a+b theta)+x(a+b(m-1)-b theta)] P_(n-1) and
c x(theta+1)(m-1-theta) P_(n-2). The full operator identity is exact;
a > 0, b >= a, c >= 0 and a nonzero denominator are checked for every index.
Positive constant/linear seeds and exact later seed steps give simple negative
roots, adjacent strict interlacing and auxiliary interlacing, with a replayable
certificate and RealRooted lagTwoEuler_recurrence_certificate correspondence.
Negative lag corrections are excluded. Tried after the existing routes.

Binomial-perturbed Euler orbit (binomial_perturbed_euler_orbit, recorded only
when its operators match): the exact lag-two theta recurrence of RealRooted
issue #2513, s >= 1, kappa >= s+4, Q_0=1, Q_1=1+2u,
Q_2=1+kappa*u+kappa*u^2, starting at degree m >= 3. Integer degree shifts,
every seed and nonzero scalar denominators are checked exactly. Also recognizes
Q_m(u)=(1+u)^m P_n(u/(1+u)) by exact Mobius transport with an inductive degree
cap. Gives simple roots in (-1,0) for Q, or (-inf,0) for P, and strict adjacent
interlacing with the direction preserved. A certificate replays the transformed
rule and seed identities. Lean: strictInterl_of_orbitPerturbation_recurrence.

Rational PF Riordan (rational_pf_riordan, recorded only when it applies):
for an index-free homogeneous rule without derivatives, build the exact reduced
row generating function from every prescribed seed. Recognize
F(x,z)=d(z)/(1-x z d(z) K(z)), with d(0)>0 and nonzero rational PF series d,K.
The rational Edrei criterion is decided by exact gcd layers, PRS root signs and
root isolation: numerator zeros nonpositive, denominator poles positive,
positive initial scalar. Chen--Fu--Ruan Theorems 3(3) and 4 give PF rows for all
indices, including K(0)=0 and noninvertible arrays. Also tries padded reversal
x^(a*m+b) P_(first_index+m)(1/x) with a certified degree bound (a=1,2,3).
Certificates replay the generating-function identity and exact PF evidence.
Zero rows are reported exactly; no simplicity or interlacing is claimed.
The route is available to index-section parents and the affine marker shifts.

When no criterion proves the record:
  - tail: if one of the first rows is zero or not real-rooted, the criteria run
    on the record starting after the last such row; a proof then says "every
    row from index N on is real-rooted" (real_rooted_from_index);
  - alternative representations: the recurrence finder's menu on generated rows
    proposes other rules (first-order derivative rules with a factor f(n), lag
    at most two); a rule is used only after an all-index equality certificate
    (Ore-algebra division and invariant row annihilators, exact), then the
    criteria run on it. --max-representation-candidates bounds each finder pass
    (deterministic count; 0 disables the stage);
  - catalog lifts: the record is connected to a family of the bundled OEIS
    catalog by an identity certified for every index: a normalized
    identification T_n = c(n) x^(a n + b) S_(n+s) (candidates from a
    fingerprint index of the catalog rows, the scalar proved nonzero), or a
    backward basis-map lift Q_n = T_P(R_n) whose map theorem carries the roots
    of R into a target interval. The parent (S, or R) is proved by the stages
    above (never by another lift); R may also be identified with a catalog
    family. When the parent is only flagged lean_proved in the catalog (a proof
    made outside polytool), the status is conditional_on_attestation, never
    proved. --max-lift-candidates bounds the catalog families certified
    (deterministic order; 0 disables the stage). Not run when a root region is
    requested;
  - affine coordinates: when still nothing proves the record, the criteria run
    on it in exact coordinates x = a u + b (Q_n(u) = P_n(a u + b): every
    coefficient substituted, D_u = a D_x on derivative terms, the seeds
    mapped), with b a rational root of a factor free of n of a coefficient (a
    = 1, then a = -1) and the reflection x = -u. A proof maps back: same rows
    nonzero and real-rooted, root region a R + b (ends swapped for a < 0),
    interlacing kept for a > 0 and mirrored for a < 0. No other stage runs in
    the new coordinate. --max-affine-coordinates bounds the coordinates tried
    (deterministic menu; 0 disables the stage).
A rule with an inhomogeneous term matches no criterion: status not_applicable.

The recurrence input is recurrence.v1 JSON, a theorem wrapper or a
representation menu; with a menu, --candidate is required unless it has one
entry. --verify replays the recorded routes of a report from standard input
(every attempt, the tail, every equality certificate and the catalog-lift
certificates with their parent reports, the affine coordinates with their
transformed records) and takes no search flags."#,
    input: InputKind::Json,
    output: OutputKind::JsonSwitch,
    flags: &[
        option(
            "recurrence",
            "FILE|-",
            "Read recurrence JSON from FILE (default -: stdin or --input FILE)",
        ),
        option("candidate", "N", "Explicit 0-based menu choice"),
        option(
            "max-start-delay",
            "N",
            "Induction-start retries (maximum 16)",
        )
        .default("3"),
        option(
            "max-representation-candidates",
            "N",
            "Finder configurations per alternative-representation pass (0: off)",
        )
        .default("300"),
        option(
            "max-lift-candidates",
            "N",
            "Catalog families certified by the catalog-lift stage (0: off)",
        )
        .default("6"),
        option(
            "max-affine-coordinates",
            "N",
            "Coordinates x = a u + b tried by the affine-coordinate stage (0: off)",
        )
        .default("8"),
        switch("verify", "Replay the recorded routes of a report from stdin; no search flags"),
    ],
    details: r#"JSON output is polytool.recurrence-proof-discovery.v1. Its trivially_real_rooted is
true when the supplied (not fitted) rule is B(n) P_n = A(x,n) P_(n-k) and
product induction proves it. real_rooted_from_index is the first index N with
every row P_n, n >= N, proved nonzero and real-rooted. A residue-class proof
may leave whole residue classes of zero rows (conclusion.zero_residue_classes,
residues of n modulo the period); real_rooted_from_index is then omitted and
real_rooted_outside_zero_residue_classes gives from_index and the zero classes
(status proved). tail and alternative_representations hold the later stages, each
alternative with its polytool.recurrence-equality.v1 certificate. catalog_lifts
holds the catalog-lift stage: identifications (catalog id, normalized-equality
certificate, the parent's own discovery report, the transfer), basis_maps (the
backward basis-map search and coefficient identifications) and conclusion
(route, catalog_family or basis_map, parent_proof polytool or
catalog_attestation, conditional_on, from_index, root_region). Status
conditional_on_attestation: the only lift rests on a catalog flag;
real_rooted_from_index is then omitted and the claim is in
catalog_lifts.conclusion. affine_coordinates holds the affine-coordinate stage:
candidates_available, coordinates (slope, intercept, reason, the transformed
record and its attempts, a transformed tail) and conclusion (coordinate,
theorem, from_index, root_region mapped to x, interlacing with
interlacing_orientation preserved or mirrored). Replay needs
the attempts of exactly the criteria the current search records; recompute
reports saved before a criterion was added."#,
    examples: &[
        "polytool recurrence --json < rows.txt | polytool recurrence-prove --recurrence - --json",
        "polytool recurrence-prove --recurrence rec.json --json > proof.json\npolytool recurrence-prove --verify < proof.json",
    ],
    exit_codes: EXIT_PROVE,
    see_also: &["recurrence-degree", "recurrence-theorems", "topic:recurrence-json"],
    run: Some(crate::recurrence_prove::cmd_recurrence_prove),
    ..COMMAND
};

pub(crate) const RECURRENCE_THEOREMS: CommandSpec = CommandSpec {
    name: "recurrence-theorems",
    group: GROUP,
    summary: "Conditional interlacing theorem matches for recurrence JSON",
    usage: &["recurrence-theorems [--recurrence <FILE|->] [--candidate N] [--json]"],
    description: r#"Check a supplied recurrence JSON record (polytool.recurrence.v1, e.g. from
`polytool recurrence --json`) against the supported interlacing theorems and
report explicit conditional theorem matches. Also accepts theorem wrappers and
recurrence representation menus (select an entry with --candidate, 0-based).

Hypotheses are checked exactly. The record's initial rows are compared with the
canonical boundary rows (finite check; `canonical_after_relabel` when the first
supplied row plays the role of index 0). A match is conditional: it never shows
that the recurrence holds for the intended sequence for all n."#,
    input: InputKind::Json,
    output: OutputKind::JsonSwitch,
    flags: &[
        option(
            "recurrence",
            "FILE|-",
            "Read recurrence JSON from FILE (default -: stdin or --input FILE)",
        ),
        option(
            "candidate",
            "N",
            "Select a representation menu entry (0-based)",
        ),
    ],
    details: r#"Theorems (paper proofs, not Lean checked):
  weighted_renewal_window          P_n = sum_j r_j t^j P_(n-1-ell*j), r_0 = 1,
                                   r_j >= 0, R(z) = sum r_j z^j real-rooted;
                                   P_n << P_(n+s) for d*s <= d*ell+1, d = deg R.
  fixed_delay_predecessor_window   P_n = P_(n-1) + w(n) t P_(n-L), w polynomial in n
                                   with w(n) >= 0 for all integers n >= 1 (exact
                                   certificate, not sampling); Q_A << Q_B for
                                   B - A <= L with Q_(max(0, n-L)) boundary.

JSON output is the polytool.recurrence-theorems.v1 report; its
trivially_real_rooted (supplied, not fitted, records only) is true when the rule
is B(n) P_n = A(x,n) P_(n-k) and product induction proves it."#,
    examples: &[
        "polytool recurrence --json < rows.txt | polytool recurrence-theorems --recurrence -",
        "polytool recurrence-theorems --recurrence menu.json --candidate 0 --json",
    ],
    exit_codes: &[
        (0, "at least one theorem matches"),
        (1, "no theorem matches"),
        (2, "usage error or invalid input"),
    ],
    see_also: &["recurrence", "recurrence-prove", "topic:recurrence-json"],
    run: Some(crate::cmd_recurrence_theorems),
    ..COMMAND
};

pub(crate) const RECURRENCE_VERONESE_LIFTS: CommandSpec = CommandSpec {
    name: "recurrence-veronese-lifts",
    group: GROUP,
    summary: "Parent sequences whose coefficient sections are the sequence",
    usage: &["recurrence-veronese-lifts [--recurrence <FILE|->] [OPTIONS]"],
    description: r#"Search for parent sequences P_n(x) whose Veronese coefficient section is the
sequence Q_n(y) given by a polytool.recurrence.v1 record (for example from
`polytool recurrence --json`):

  V_(r,j)(P)(y) = sum_k [x^(r k + j)] P(x) y^k,   Q_n = V_(r,j)(P_n) for all n >= n0.

If P_n is real-rooted with coefficients of one sign, every section is zero or
real-rooted with nonpositive roots (Polya frequency sequences,
Aissen-Schoenberg-Whitney-Edrei), so a lift transfers real-rootedness to Q_n.

Supported input: derivative-free terms without (-1)^n factors and with
coefficients independent of n, a constant (or no) denominator, no inhomogeneous
term. Anything else exits with status 2. Supplied rows that do not follow the
recurrence (zero or other prefix rows) are dropped: the search runs on the
shortest tail of rows that does, and the report names its first index
(tail_first_index). Theorem wrappers and representation menus are accepted
(select an entry with --candidate, 0-based). Without --recurrence the record is
read from standard input.

The search is bounded: no lift found does not mean that none exists."#,
    input: InputKind::Json,
    output: OutputKind::JsonSwitch,
    flags: &[
        option("recurrence", "FILE|-", "Read recurrence JSON from FILE (default: stdin)"),
        option("candidate", "N", "Select a representation menu entry (0-based)"),
        option("max-modulus", "R", "Largest modulus, 2..=6").default("4"),
        option("parent-rows", "N", "Parent rows checked exactly").default("12"),
        option("max-candidates", "N", "Candidate parents G tried per modulus").default("256"),
    ],
    details: r#"Method, for each modulus r = 2..max-modulus (at most the order d):
  chi(L, y) = L^d - sum_k c_k(y) L^(d-k) is the characteristic polynomial of
  Q_m = sum_k c_k(y) Q_(m-k). chi(L, x^r) is factored over Q; each monic divisor
  G(L, x) with r deg_L G <= d is a candidate parent recurrence. G is kept when
  N = prod_(w^r=1) G(L, w x) divides chi; for each residue j the parent seeds
  are solved exactly (polynomial seeds required) and the first d rows are
  compared. N | chi, chi monic and d matching rows prove the identity for all n.

Status of a lift:
  proved_product_parent        order-one parent P_n = c^(n-n0) P_(n0) with c and
                               P_(n0) real-rooted of one sign: every Q_n (n >= n0)
                               is zero or real-rooted with nonpositive roots
  proved_parent_route          another parent, proved zero or real-rooted of one
                               sign for every n by the shared parent provers
                               (theorem routes of recurrence-prove): same
                               conclusion for every Q_n
  parent_real_rootedness_open  all checked parent rows are real-rooted of one
                               sign (finite evidence, not a proof)
  parent_row_zero, parent_row_mixed_signs, parent_row_not_real_rooted
                               a checked parent row fails: no proof from this lift
Lifts are listed proved first, then open, then failing. Rejected candidates give
a reason: norm_does_not_divide, seeds_not_determined, seeds_not_polynomial or
row_mismatch.

JSON output is the polytool.recurrence-veronese-lifts.v1 report (exact rationals
as strings; each parent as a polytool.recurrence.v1 record in
parent_recurrence)."#,
    examples: &[
        "echo '{\"schema\":\"polytool.recurrence.v1\",\"first_index\":0,\n  \"initial_polynomials\":[[\"1\"],[\"1\"]],\"recurrence\":{\"terms\":[\n  {\"offset\":1,\"deriv_order\":0,\"sign\":\"none\",\"coeff\":{\"coeffs\":[[\"2\"]]}},\n  {\"offset\":2,\"deriv_order\":0,\"sign\":\"none\",\"coeff\":{\"coeffs\":[[\"-1\",\"1\"]]}}]}}' \\\n  | polytool recurrence-veronese-lifts",
        "polytool recurrence-veronese-lifts --recurrence rec.json --max-modulus 3 --json",
    ],
    exit_codes: &[
        (0, "at least one lift found"),
        (1, "no lift found"),
        (2, "invalid or unsupported input"),
    ],
    see_also: &["recurrence", "recurrence-transform", "topic:recurrence-json"],
    run: Some(crate::cmd_recurrence_veronese_lifts),
    ..COMMAND
};

pub(crate) const RECURRENCE_IMPLICIT: CommandSpec = CommandSpec {
    name: "recurrence-implicit",
    group: GROUP,
    summary: "Diagonal maps and antecedents: implicit relations, uniqueness",
    usage: &[
        "recurrence-implicit [--compact] < request.json",
        "recurrence-implicit --verify < report.json",
        "recurrence-implicit --assess < description.json",
        "recurrence-implicit --search [--compact] < search.json",
        "recurrence-implicit --verify-search < search-report.json",
    ],
    description: r#"Map a sequence with a defining rule by a diagonal coefficient map or a derivative
antecedent, transport the rule into an implicit differential relation
(derivatives on any row) with exceptional-index conditions, and assess the
resulting generating description: annihilation evidence and a uniqueness
argument are reported separately. Prints a polytool.recurrence-implicit.v1
report. All bounds belong in the JSON request on standard input.

The modes are exclusive: transport (default), --verify, --assess, --search and
--verify-search."#,
    input: InputKind::Json,
    output: OutputKind::JsonOnly,
    flags: &[
        switch("verify", "Replay a report exactly and re-check every certificate"),
        switch("assess", "Assess a generating description and print the assessment"),
        switch("search", "Search implicit relations directly (JSON search request)"),
        switch("verify-search", "Replay a search report and re-check every kernel certificate"),
        switch("compact", "Print compact JSON"),
    ],
    details: r#"Request fields:
  transformation   {"kind": "diagonal", "weight": W} with W one of
                     {"kind":"factorial_power","power":m}  lambda = 1/(k!)^m
                     {"kind":"binomial","m":m}             lambda = binom(k+m,m), Q = W_m P
                     {"kind":"inverse_binomial","m":m}     lambda = 1/binom(k+m,m), P = W_m Q
                     {"kind":"product","scale":"c","factorial_exponent":p,
                      "factors":[{"coeffs":{"coeffs":[[..]]},"exponent":e}]}
                                     lambda = c (k!)^p prod f(n,k)^e
                   or {"kind": "derivative_antecedent", "g": ["g0", ...]}:
                     Q_n with P_n = D^m(g Q_n)/m!, m = deg g (exact back substitution)
                   ({"kind": "derivative", "order": m}, Q_n = D^m P_n, is lossy:
                   it is searched and transported by --search and the representation
                   search, which record the lost coefficients)
  source_recurrence  defining polytool.recurrence.v1 rule with seeds, or
  oeis_id          a bundled OEIS family (e.g. "A065826")
  row_count        rows generated from the source (default 10)
  source_relations [{relation, first_base}] further relations of P, transported
  target_relations [{relation, first_base}] relations of Q to store (checked on rows)
  boundary_conditions [{condition, first_base}] coefficient conditions of Q
  search           {} or bounds: fit a successor rule to Q and try to prove it
                   by a reconstruction certificate L o T = sum H_nu S^nu M
  limits           {max_rows, max_degree, ...} resource limits
  root_domain      optional {} or {candidate_intervals, assumptions, ...} as in
                   recurrence-transform: multiplier-sequence (diagonal) and Rolle
                   (derivative antecedent) adapters with their direction, Rolle
                   inductions on the defining rule of P and on a certified fit of
                   Q, assumptions, a fixed-cap Jensen test and exact root enclosures

A relation is {"terms": [{"shift": i, "deriv_order": d, "coeff": {"coeffs": ...}}]}
meaning sum c(n,x) D^d Q[n+i] = 0 (coeffs[a][b] multiplies n^a x^b); a condition
is {"terms": [{"shift": i, "coefficient_index": k, "coeff": [..]}]} meaning
sum c(n) [x^k] Q[n+i] = 0 (coeff ascending in n).

Report:
  map_conditions     weights defined and nonzero: certified_all_indices or
                     checked_supplied_rows (uncertified, visible)
  transports         transported relation, clearing certificate (common
                     multiplier H, exceptional indices), reconstruction certificate
  description        seeds, relations, exceptional_index_conditions,
                     boundary_conditions
  assessment         equivalent_description | certified_consequence |
                     finite_evidence | inconsistent_seeds, with the uniqueness argument
  additional_relations  stored relations outside the description (row equations)
  independent_search finite_fit, or certified_by_reconstruction after an upgrade
  root_domain        when requested: its request, implications, claims (with
                     evidence labels), Rolle theorems, jensen_fixed_cap (fixed_cap
                     evidence only), assumption checks and observations

Modes:
  --verify         Replay a report exactly and re-check every certificate
  --assess         Assess a generating description (for example after removing a
                   condition) and print the assessment
  --search         Search implicit relations directly
                   (polytool.recurrence-implicit-search.v1)
  --verify-search  Replay a search report and re-check every kernel certificate

Direct search (--search) request fields:
  source_recurrence | oeis_id | rows (+ row_origin)   the rows P (rows only: finite fits)
  transformation   optional scalar map Q = T(P) applied first (diagonal,
                   derivative_antecedent, derivative {order}, Moebius-class maps, ...)
  supports         [{role, max_shift, derivative_order, index_degree, variable_degree}]
                   role successor (no derivative on row n+R) | implicit_successor
                   (derivatives on every row) | row_equation (max_shift 0); default:
                   implicit_successor (1,1,1,1) and row_equation (0,2,1,1)
  row_count, domain_start, conditions (default true), condition_index_degree (1),
  limits           complete-space limits {max_rows, max_columns, max_degree, ...}
Each support is one exact kernel with a rank certificate. With a constant-index
monic pure-lag source and the identity, factorial_power or derivative map a
finite theorem makes the kernel the complete space (certified_complete_space);
otherwise relations are finite_fit, or certified_by_division when
M o T = sum H S^nu L for the defining rule L. Moebius-class maps (multiply_power,
divide_power, affine; reversal, mobius, f_polynomial with an unbounded schedule)
need k(max(u+s)+1) residuals (transported-source theorem), and diagonal weights
c (k!)^(-m) prod f_j(n,k)^(e_j), e_j > 0 (binomial, product), need
k(max(u+s)+w+1) for m = 0 or k(U+m(A+J)+w+1), w the total (n,k) degree (weight
theorem); completeness_basis names the theorem. Exceptional indices of a leader
get a coefficient-condition search; every representative carries a generating
description and assessment."#,
    examples: &[
        "echo '{\"oeis_id\":\"A065826\",\"row_count\":9,\"search\":{},\n  \"transformation\":{\"kind\":\"diagonal\",\n  \"weight\":{\"kind\":\"inverse_binomial\",\"m\":1}}}' \\\n  | polytool recurrence-implicit",
        "echo '{\"oeis_id\":\"A007318\",\"transformation\":{\"kind\":\"diagonal\",\n  \"weight\":{\"kind\":\"factorial_power\",\"power\":1}},\n  \"supports\":[{\"role\":\"implicit_successor\",\"max_shift\":1,\"derivative_order\":1}]}' \\\n  | polytool recurrence-implicit --search",
    ],
    exit_codes: EXIT_REQUEST,
    see_also: &[
        "recurrence-transform",
        "recurrence-search-representations",
        "topic:json-schemas",
    ],
    run: Some(crate::cmd_recurrence_implicit),
    ..COMMAND
};

pub(crate) const RECURRENCE_SEARCH_REPRESENTATIONS: CommandSpec = CommandSpec {
    name: "recurrence-search-representations",
    group: GROUP,
    summary: "Bounded search over equivalent representations and relations",
    usage: &[
        "recurrence-search-representations [--compact] < request.json",
        "recurrence-search-representations --verify < report.json",
    ],
    description: r#"Search a bounded queue of equivalent representations of one polynomial sequence:
apply the allowed coordinate changes, transport the known certified relations,
certify complete spaces after factorial maps, fit the new coordinates, and feed
certified relations into later steps. Prints a
polytool.recurrence-representation-search.v1 report. All bounds belong in the
JSON request on standard input."#,
    input: InputKind::Json,
    output: OutputKind::JsonOnly,
    flags: &[
        switch(
            "verify",
            "Replay a report: rerun the search and every step and certificate",
        ),
        switch("compact", "Print compact JSON"),
    ],
    details: r#"Request fields:
  source_recurrence  defining polytool.recurrence.v1 rule with seeds, or
  oeis_id            a bundled OEIS family (e.g. "A007318"), or
  rows               rows only (exploratory: every description is a finite fit)
  row_origin         original index of rows[0] (default: source first index, else 0)
  row_count          rows generated from a source (default 12)
  transformations    the allowed maps (recurrence-transform / recurrence-implicit
                     registry: multiply_power, divide_power, reversal, affine,
                     mobius, f_polynomial, gamma, symmetric_decomposition,
                     diagonal, derivative_antecedent)
                     and the lossy derivative {"kind":"derivative","order":m}:
                     its lost coefficients are recorded (lost_data) and the parent
                     is assessed through Q with and without them
  search             fit bounds {max_rec_len, max_var_deg, max_idx_deg,
                     max_diff_deg, try_denominator, max_candidates}
                     (default 2, 2, 1, 1, true, 2000)
  fits               run the fit menu in every representation (default true)
  certified_spaces   certify complete spaces after factorial maps (default true),
                     and after Moebius-class maps, polynomial diagonal weights and
                     derivative maps (implicit_spaces)
  space_shapes       explicit complete-space supports (default: a bounded menu)
  root_domain        optional root-domain analysis of every coordinate change
  limits             {max_depth (2), max_representations (16), max_steps (160),
                     max_sources_per_edge (2), max_space_shapes (16), max_rows,
                     max_degree, max_search_candidates, max_search_degree,
                     max_operations, max_coefficient_bits}

Report:
  representations    the source (0) and every new representation: chain of maps,
                     certified map normal form, rows, exact reconstruction of
                     the source, embedded replayable steps and complete spaces,
                     and ranked descriptions labelled
                       finite_fit | certified_annihilating_relation |
                       equivalent_generating_description
  edges              every map application: new_representation | duplicate_map
                     (certified equal normal forms; rows are never compared) |
                     rejected | skipped_budget
  ranked             all descriptions by shift length, derivative order,
                     coefficient degree, index degree, components, reconstruction
  budget, status     limits that cut work; status complete | budget_exhausted
  steps[k].implicit_relation  a certified (also implicit-only) parent relation
                     transported through diagonal, derivative_antecedent or
                     derivative maps, with carried conditions and, for lossy maps,
                     the parent's description with and without the lost data
  steps[k].vector_transport  a certified vector system of a vector parent
                     transported exactly through a componentwise Moebius-class map
  implicit_spaces    polytool.recurrence-implicit-search.v1 complete spaces after
                     Moebius-class maps, polynomial diagonal weights and derivatives

Vector representations: after symmetric_decomposition the pair (a_n, b_n) is
mapped componentwise by every scalar bijective map (a schedule acts on b at
d_n - 1, its symmetry index), fitted by first-order vector systems (vector_fit),
and its certified systems are transported through Moebius-class maps. Vector
descriptions carry vector_assessment (uniqueness: the shared left factor is a
nonzero polynomial at every index). A rank-two gamma transport is reported as
transported_gamma_system (even part and odd companion), certified only where the
symmetry of the rows is certified."#,
    examples: &[
        "echo '{\"oeis_id\":\"A007318\",\"limits\":{\"max_depth\":1},\n  \"transformations\":[{\"kind\":\"diagonal\",\n  \"weight\":{\"kind\":\"factorial_power\",\"power\":1}}]}' \\\n  | polytool recurrence-search-representations",
        "polytool recurrence-search-representations --verify --input report.json",
    ],
    exit_codes: EXIT_REQUEST,
    see_also: &[
        "recurrence-transform",
        "recurrence-implicit",
        "topic:json-schemas",
    ],
    run: Some(crate::cmd_recurrence_search_representations),
    ..COMMAND
};
