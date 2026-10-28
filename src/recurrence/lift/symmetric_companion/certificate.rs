//! The certificate of a symmetric-companion lift, its construction from a
//! candidate, and [`verify_symmetric_companion_certificate`], which replays it
//! from the recorded data alone.

use super::operators::{
    apply_to_row, divide_by_x, evaluate_at_index, index_only, is_symmetric_at, scalar,
    CompanionParameters,
};
use super::source::PreparedSource;
use crate::real_rootedness::check_weak_interlacing_bigint_coeffs;
use crate::recurrence::implicit::algebra::{poly_text, roots_in_domain, Biv, DOp};
use crate::recurrence::implicit::biv_json;
use crate::recurrence::*;
use crate::vec_poly;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};

type Q = BigRational;

/// The Lean replay status of every certificate of this milestone: no Lean
/// theorem is replayed.
pub const SYMMETRIC_COMPANION_LEAN_REPLAY_NOT_ATTEMPTED: &str = "not_attempted";

/// The parameters of a lift as polynomials in the index `n`: entry `i` of each
/// vector is the coefficient of `n^i`, an exact integer or `num/den` string.
///
/// `delta = c - epsilon` and `b = a + c d + delta`, with the schedule
/// `d_n = d0 + (n - n0)`; the verifier rejects a record violating either.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymmetricCompanionParameters {
    /// `a_n`, the constant coefficient of the parent multiplier.
    pub a: Vec<String>,
    /// `b_n`, the linear coefficient of the parent multiplier.
    pub b: Vec<String>,
    /// `c_n`, the derivative coefficient of the parent rule.
    pub c: Vec<String>,
    /// `delta_n = b_n - a_n - c_n d_n`.
    pub delta: Vec<String>,
    /// `epsilon_n = c_n - delta_n`, the coupling of the companion.
    pub epsilon: Vec<String>,
}

/// A replayable certificate: the source, the schedule, the parameters, the
/// eliminated operator with its derivation from the source rule, and the
/// collapsed parent.
///
/// Everything needed by [`verify_symmetric_companion_certificate`] is
/// recorded; the verifier recomputes every derived item and rejects a record
/// that disagrees with its own parameters.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymmetricCompanionCertificate {
    /// The target sequence `A_n`: its scalar recurrence and seed rows
    /// (`polytool.recurrence.v1`); `first_index` is `n0`, the first index of
    /// the reconstruction domain `n >= n0`.
    pub source: RecurrenceJson,
    /// `d_(n0)`; the schedule is `d_n = d_(n0) + (n - n0)`.
    pub symmetry_degree_at_first_index: usize,
    /// The parameter witnesses.
    pub parameters: SymmetricCompanionParameters,
    /// The eliminated operator `E` of the module documentation, base `n`.
    pub eliminated_relation: ImplicitRelationJson,
    /// `g(n) E = sum_nu H_nu S^nu M` with `M` the source rule in forward form:
    /// `source_relation` is `g(n) E`, the map is the identity.  `None` when no
    /// derivation was found (the identity then stays a candidate).
    pub source_consequence: Option<ReconstructionCertificate>,
    /// The scalar factor `g(n)` of `source_consequence` (coefficients of
    /// `n^i`).
    pub source_consequence_factor: Vec<String>,
    /// Bases `n >= n0` where `(E A)_n = 0` is not covered by
    /// `source_consequence` (the source rule does not hold there, or
    /// `g(n) = 0`); each is checked on exact rows.
    pub exceptional_bases: Vec<usize>,
    /// `B_(n0) = R_(n0) / (epsilon_(n0) x)`, ascending coefficients.
    pub companion_seed: Vec<String>,
    /// The parent `P_(n+1) = (a_n + b_n x) P_n + c_n x (1 - x) D P_n` with seed
    /// `P_(n0) = A_(n0) + x B_(n0)`, as a `polytool.recurrence.v1` record.
    pub parent: RecurrenceJson,
}

/// Status of one obligation or hypothesis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SymmetricCompanionCheckStatus {
    /// Established for every index of the domain.
    Verified,
    /// Refuted exactly.
    Failed,
    /// Neither established nor refuted by the implemented arguments.
    Unresolved,
    /// Outside the supported class (for example a vanishing coupling).
    Unsupported,
}

/// One named obligation or hypothesis with its status and evidence.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymmetricCompanionCheck {
    /// Stable snake-case name.
    pub name: String,
    /// What is claimed, for every `n >= n0` unless stated otherwise.
    pub statement: String,
    /// The status.
    pub status: SymmetricCompanionCheckStatus,
    /// How the status was obtained.
    pub evidence: String,
}

/// Status of the all-index algebraic reconstruction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SymmetricCompanionIdentityStatus {
    /// Every identity obligation is verified for all `n >= n0`.
    Certified,
    /// No obligation fails, but some are unresolved.
    Candidate,
    /// The candidate leaves the supported class (for example `epsilon_n = 0`
    /// at some index); nothing is cancelled away.
    Unsupported,
    /// Some obligation fails exactly; this is not a companion of `A_n`.
    Rejected,
}

/// Status of the real-rootedness route through the preservation theorem.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SymmetricCompanionRealRootednessStatus {
    /// The identity is certified and every hypothesis of the preservation
    /// theorem holds: every `A_n`, `B_n` and `P_n` (`n >= n0`) is zero or
    /// real-rooted with nonnegative coefficients, and `B_n <= A_n`.
    TheoremApplies,
    /// The identity is certified but some hypothesis is not verified.  The
    /// route fails; this says nothing against real-rootedness of `A_n`.
    HypothesesNotVerified,
    /// The identity is not certified, so the route is not evaluated.
    NotEvaluated,
}

/// The result of replaying a certificate.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymmetricCompanionVerification {
    /// The all-index reconstruction status.
    pub identity_status: SymmetricCompanionIdentityStatus,
    /// The real-rootedness route status (paper theorem).
    pub real_rootedness_status: SymmetricCompanionRealRootednessStatus,
    /// Lean replay: always [`SYMMETRIC_COMPANION_LEAN_REPLAY_NOT_ATTEMPTED`].
    pub lean_replay: String,
    /// The identity obligations.
    pub identities: Vec<SymmetricCompanionCheck>,
    /// The hypotheses of the preservation theorem.
    pub hypotheses: Vec<SymmetricCompanionCheck>,
    /// Names and statements of the obligations and hypotheses that are not
    /// verified.
    pub unresolved_obligations: Vec<String>,
    /// Number of indices covered by the exact row self-check.
    pub self_checked_rows: usize,
    /// One-paragraph conclusion.
    pub conclusion: String,
}

impl SymmetricCompanionVerification {
    /// Whether the all-index identity is certified.
    pub fn is_certified(&self) -> bool {
        self.identity_status == SymmetricCompanionIdentityStatus::Certified
    }

    /// Whether the preservation theorem applies, so every `A_n` is zero or
    /// real-rooted.
    pub fn theorem_applies(&self) -> bool {
        self.real_rootedness_status == SymmetricCompanionRealRootednessStatus::TheoremApplies
    }
}

// ---------------------------------------------------------------------------
// Encoding helpers
// ---------------------------------------------------------------------------

pub(super) fn encode(p: &[Q]) -> Vec<String> {
    vec_poly::trim(p)
        .iter()
        .map(format_rational_coeff)
        .collect()
}

fn decode(p: &[String], what: &str) -> Result<Vec<Q>, RecurrenceError> {
    if p.is_empty() {
        return Err(RecurrenceError::InvalidInput(format!(
            "{what} has no coefficients"
        )));
    }
    let values = p
        .iter()
        .map(|c| {
            parse_rational_coeff(c)
                .map_err(|e| RecurrenceError::from(e).with_context(format!("invalid {what}")))
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(vec_poly::trim(&values))
}

impl SymmetricCompanionParameters {
    pub(super) fn from_parameters(parameters: &CompanionParameters) -> Self {
        Self {
            a: encode(&parameters.a),
            b: encode(&parameters.b()),
            c: encode(&parameters.c),
            delta: encode(&parameters.delta()),
            epsilon: encode(&parameters.epsilon),
        }
    }

    /// The exact values `(a, b, c, delta, epsilon)` as polynomials in `n`.
    #[allow(clippy::type_complexity)]
    pub fn values(&self) -> Result<(Vec<Q>, Vec<Q>, Vec<Q>, Vec<Q>, Vec<Q>), RecurrenceError> {
        Ok((
            decode(&self.a, "parameter a")?,
            decode(&self.b, "parameter b")?,
            decode(&self.c, "parameter c")?,
            decode(&self.delta, "parameter delta")?,
            decode(&self.epsilon, "parameter epsilon")?,
        ))
    }
}

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

/// The identity map as a shift-free relation.
fn identity_relation() -> ImplicitRelation {
    ImplicitRelation::from_op(&DOp::single(0, 0, Biv::one()))
}

/// `g(n) E = sum_nu H_nu S^nu M`, with `g = 1` (right division by the source
/// rule, which needs a constant leader) or, for an order-two source,
/// `g = f(n + 2)` the forward leader with `f E = epsilon_n M`.
pub(super) fn source_consequence(
    source: &PreparedSource,
    eliminated: &DOp,
    epsilon: &[Q],
) -> Option<(ReconstructionCertificate, Vec<Q>)> {
    let eliminated_relation = ImplicitRelation::from_op(eliminated);
    if let Some(certificate) = find_reconstruction_certificate(
        &eliminated_relation,
        &identity_relation(),
        &source.relation,
    ) {
        return Some((certificate, vec![Q::one()]));
    }
    if source.order != 2 {
        return None;
    }
    let factor = source.forward_leader();
    let left = scalar(&factor).compose(eliminated);
    let right = scalar(epsilon).compose(&source.relation.op());
    if left.is_zero() || left != right {
        return None;
    }
    let operator_term = |coeff: &Biv| OperatorTermJson {
        deriv_order: 0,
        coeff: biv_json(coeff),
    };
    let certificate = ReconstructionCertificate {
        source_relation: ImplicitRelation::from_op(&left).to_json(),
        map: vec![operator_term(&Biv::one())],
        target_relation: source.relation.to_json(),
        multipliers: vec![ReconstructionMultiplier {
            shift: 0,
            differentiations: vec![operator_term(&index_only(epsilon))],
        }],
        statement: "f(n + 2) E = epsilon_n M (both rules have a leader free of x and D): every \
                    sequence satisfying M at base n satisfies f(n + 2) (E A)_n = 0"
            .into(),
    };
    certificate.verify_identity().ok()?;
    Some((certificate, factor))
}

/// The bases `n >= n0` not covered by the source consequence: those where the
/// source rule is not known to hold (`n < n0 + K - R`) and the roots of
/// `g(n)` beyond them.
pub(super) fn exceptional_bases(source: &PreparedSource, factor: &[Q]) -> Vec<usize> {
    let regular = source.first_regular_base();
    let mut bases: Vec<usize> = (source.first_index..regular).collect();
    bases.extend(roots_in_domain(factor, regular, None));
    bases.sort_unstable();
    bases.dedup();
    bases
}

/// `R_n = A_(n+1) - U_n A_n`.
fn remainder(parameters: &CompanionParameters, rows: &[Vec<Q>], n: usize) -> Vec<Q> {
    let ops = parameters.operators();
    let offset = n - parameters.first_index;
    vec_poly::sub(&rows[offset + 1], &apply_to_row(&ops.u, &rows[offset], n))
}

/// `B_n = R_n / (epsilon_n x)`, when `epsilon_n != 0` and `x | R_n`.
fn companion_row(parameters: &CompanionParameters, rows: &[Vec<Q>], n: usize) -> Option<Vec<Q>> {
    let epsilon = evaluate_at_index(&parameters.epsilon, n);
    if epsilon.is_zero() {
        return None;
    }
    let quotient = divide_by_x(&remainder(parameters, rows, n))?;
    Some(vec_poly::trim(&vec_poly::scale(
        &quotient,
        &epsilon.recip(),
    )))
}

/// The parent rule `P_m = (a_(m-1) + b_(m-1) x) P_(m-1) + c_(m-1) x (1 - x) D P_(m-1)`.
fn parent_recurrence(parameters: &CompanionParameters) -> Recurrence {
    let ops = parameters.operators();
    let mut terms = Vec::new();
    for ((_, derivative), coefficient) in &ops.parent.0 {
        terms.push(RecurrenceTerm::new(
            1,
            *derivative,
            coefficient.shift_n(-1).to_bivar(),
        ));
    }
    Recurrence {
        terms,
        denominator: None,
        inhomogeneous: None,
    }
}

/// Build the certificate of a candidate; `Err` with a reason when the
/// companion seed `B_(n0)` does not exist.
pub(super) fn build_certificate(
    source: &PreparedSource,
    parameters: &CompanionParameters,
) -> Result<SymmetricCompanionCertificate, String> {
    let n0 = source.first_index;
    let ops = parameters.operators();
    let (consequence, factor) =
        match source_consequence(source, &ops.eliminated, &parameters.epsilon) {
            Some((certificate, factor)) => (Some(certificate), factor),
            None => (None, vec![Q::one()]),
        };
    let exceptional = exceptional_bases(source, &factor);
    let rows = source.rows(2).map_err(|e| e.to_string())?;
    if evaluate_at_index(&parameters.epsilon, n0).is_zero() {
        return Err(format!(
            "epsilon vanishes at the first index n0 = {n0}, so B_(n0) is not defined"
        ));
    }
    let companion = companion_row(parameters, &rows, n0)
        .ok_or_else(|| format!("R_(n0)(0) != 0 at n0 = {n0}: A_(n0+1)(0) != a_(n0) A_(n0)(0)"))?;
    let parent_seed = vec_poly::add(&rows[0], &vec_poly::shift(&companion));
    let parent = RecurrenceJson::from_recurrence_rational(
        &parent_recurrence(parameters),
        n0,
        &[vec_poly::trim(&parent_seed)],
        None,
    );
    Ok(SymmetricCompanionCertificate {
        source: source.to_json(),
        symmetry_degree_at_first_index: parameters.symmetry_degree_at_first_index,
        parameters: SymmetricCompanionParameters::from_parameters(parameters),
        eliminated_relation: ImplicitRelation::from_op(&ops.eliminated).to_json(),
        source_consequence: consequence,
        source_consequence_factor: encode(&factor),
        exceptional_bases: exceptional,
        companion_seed: encode(&companion),
        parent,
    })
}

// ---------------------------------------------------------------------------
// Verification
// ---------------------------------------------------------------------------

fn check(
    name: &str,
    statement: impl Into<String>,
    status: SymmetricCompanionCheckStatus,
    evidence: impl Into<String>,
) -> SymmetricCompanionCheck {
    SymmetricCompanionCheck {
        name: name.into(),
        statement: statement.into(),
        status,
        evidence: evidence.into(),
    }
}

fn mismatch(what: &str) -> RecurrenceError {
    RecurrenceError::VerificationFailed(format!(
        "the recorded {what} does not match the one recomputed from the parameters"
    ))
}

/// `p(n) >= 0` for every integer `n >= n0`: verified when `p(n0 + u)` has
/// nonnegative coefficients; for degree at most one the converse holds, so a
/// failure is exact; otherwise unresolved.
fn nonnegative_on_domain(p: &[Q], n0: usize) -> (SymmetricCompanionCheckStatus, String) {
    let shifted = vec_poly::trim(&vec_poly::compose(
        p,
        &[Q::from_integer(n0.into()), Q::one()],
    ));
    let text = poly_text(&shifted, "u");
    if shifted.iter().all(|c| !c.is_negative()) {
        return (
            SymmetricCompanionCheckStatus::Verified,
            format!("with n = {n0} + u the polynomial is {text}, with nonnegative coefficients"),
        );
    }
    if vec_poly::degree(&shifted).unwrap_or(0) <= 1 {
        return (
            SymmetricCompanionCheckStatus::Failed,
            format!(
                "with n = {n0} + u the polynomial is {text}, of degree at most one with a \
                 negative coefficient, so it is negative for some integer n >= {n0}"
            ),
        );
    }
    (
        SymmetricCompanionCheckStatus::Unresolved,
        format!("with n = {n0} + u the polynomial is {text}; its sign is not decided"),
    )
}

/// Positive integer multiple of a rational row (same roots, same signs).
fn integer_row(p: &[Q]) -> Vec<BigInt> {
    let denominator = p.iter().fold(BigInt::one(), |acc, c| acc.lcm(c.denom()));
    p.iter()
        .map(|c| (c * Q::from_integer(denominator.clone())).to_integer())
        .collect()
}

/// `b <= a` in the directed weak-interlacing convention of the preservation
/// theorem: nonzero rows interlace with `a` holding the rightmost root, which
/// is `check_weak_interlacing(b, a) == Some(true)`; a zero member is allowed
/// when the other is real-rooted (and `(0, 0)` holds).
fn companion_interlaces(b: &[Q], a: &[Q]) -> (bool, String) {
    let b_int = integer_row(b);
    let a_int = integer_row(a);
    match (vec_poly::is_zero(b), vec_poly::is_zero(a)) {
        (true, true) => (true, "both members are zero".into()),
        (true, false) | (false, true) => {
            let other = if vec_poly::is_zero(b) { &a_int } else { &b_int };
            let real = crate::real_rootedness::is_real_rooted_bigint_coeffs(other);
            (
                real,
                format!(
                    "one member is zero; the other is {}real-rooted (exact)",
                    if real { "" } else { "not " }
                ),
            )
        }
        (false, false) => {
            let verdict = check_weak_interlacing_bigint_coeffs(&b_int, &a_int);
            (
                verdict == Some(true),
                format!(
                    "check_weak_interlacing(B, A) = {verdict:?}: B << A means the roots weakly \
                     alternate with A holding the rightmost root"
                ),
            )
        }
    }
}

struct Outcome {
    identities: Vec<SymmetricCompanionCheck>,
    hypotheses: Vec<SymmetricCompanionCheck>,
    self_checked_rows: usize,
}

/// Replay a certificate from its recorded data, independently of the search.
///
/// The verifier recomputes the eliminated operator, the companion seed and the
/// parent from the parameters (a disagreement is a
/// [`RecurrenceError::VerificationFailed`]) and then decides every obligation
/// of the module documentation:
///
/// 1. `eliminated_annihilates_source`: `source_consequence` passes
///    [`ReconstructionCertificate::verify_identity`], its relations are
///    `g(n) E` and the source rule, and `(E A)_n = 0` holds exactly at every
///    exceptional base;
/// 2. `coupling_nonvanishing`: `epsilon_n != 0` for every `n >= n0` (exact
///    integer roots); a vanishing coupling is unsupported;
/// 3. `constant_term_divisibility`: the derivative terms of `E` vanish at
///    `x = 0`, the projection factors as `(epsilon_n S - g_n) o (S - a_n)`, and
///    `A_(n0+1)(0) = a_(n0) A_(n0)(0)`;
/// 4. `second_coupled_equation`: the elimination identity holds as an
///    operator identity, which with 1--3 gives the second equation;
/// 5. `symmetry_schedule`: `d_(n0) >= 1` and the two seeds are symmetric;
/// 6. `parent_collapse`: the collapse identities hold as operator identities;
/// 7. `row_self_check`: the coupled system, symmetry and parent rows agree on
///    `self_check_rows` exact rows (a defensive check, not part of the proof).
///
/// The hypotheses of the preservation theorem are then checked exactly:
/// `a_n >= 0`, `0 <= delta_n <= c_n` on the integer domain, `d_(n0) >= 1`,
/// nonnegative coefficients of `A_(n0)`, `B_(n0)`, and `B_(n0) <= A_(n0)`.
pub fn verify_symmetric_companion_certificate(
    certificate: &SymmetricCompanionCertificate,
    self_check_rows: usize,
) -> Result<SymmetricCompanionVerification, RecurrenceError> {
    let (recurrence, first_index, seeds) = certificate.source.to_recurrence_parts()?;
    let source = PreparedSource::new(&recurrence, &seeds, first_index)?;
    let (a, b, c, delta, epsilon) = certificate.parameters.values()?;
    let parameters = CompanionParameters::new(
        a,
        c,
        epsilon,
        first_index,
        certificate.symmetry_degree_at_first_index,
    );
    if parameters.delta() != delta {
        return Err(RecurrenceError::VerificationFailed(
            "the recorded delta is not c - epsilon".into(),
        ));
    }
    if parameters.b() != b {
        return Err(RecurrenceError::VerificationFailed(
            "the recorded b is not a + c d + delta".into(),
        ));
    }
    let ops = parameters.operators();
    let recorded = ImplicitRelation::from_json(&certificate.eliminated_relation)?.op();
    if recorded != ops.eliminated {
        return Err(mismatch("eliminated operator"));
    }
    let factor = decode(
        &certificate.source_consequence_factor,
        "source consequence factor",
    )?;
    if exceptional_bases(&source, &factor) != certificate.exceptional_bases {
        return Err(mismatch("list of exceptional bases"));
    }
    let outcome = verify_obligations(certificate, &source, &parameters, &factor, self_check_rows)?;
    Ok(assemble(outcome))
}

fn verify_obligations(
    certificate: &SymmetricCompanionCertificate,
    source: &PreparedSource,
    parameters: &CompanionParameters,
    factor: &[Q],
    self_check_rows: usize,
) -> Result<Outcome, RecurrenceError> {
    use SymmetricCompanionCheckStatus::*;
    let n0 = source.first_index;
    let d0 = parameters.symmetry_degree_at_first_index;
    let ops = parameters.operators();
    let last_exceptional = certificate.exceptional_bases.last().copied().unwrap_or(n0);
    let row_count = (self_check_rows.max(2) + 2).max(last_exceptional - n0 + 3);
    let rows = source.rows(row_count)?;
    let mut identities = Vec::new();

    // 1. E annihilates A for every n >= n0.
    let annihilation_statement =
        "(E A)_n = 0 for every n >= n0, where E is the eliminated operator";
    identities.push(match &certificate.source_consequence {
        None => check(
            "eliminated_annihilates_source",
            annihilation_statement,
            Unresolved,
            "no derivation of g(n) E from the source rule was found (right division needs a \
             constant leader; proportionality needs an order-two source)",
        ),
        Some(consequence) => {
            let source_relation = ImplicitRelation::from_json(&consequence.source_relation)?.op();
            let target_relation = ImplicitRelation::from_json(&consequence.target_relation)?.op();
            let identity_map = consequence.map.len() == 1
                && consequence.map[0].deriv_order == 0
                && biv_from_json(&consequence.map[0].coeff)? == Biv::one();
            if source_relation != scalar(factor).compose(&ops.eliminated)
                || target_relation != source.relation.op()
                || !identity_map
            {
                return Err(mismatch("source consequence"));
            }
            consequence.verify_identity()?;
            let failing = certificate.exceptional_bases.iter().copied().find(|&base| {
                PreparedSource::apply(&ops.eliminated, &rows, n0, base)
                    .is_none_or(|r| !vec_poly::is_zero(&r))
            });
            match failing {
                Some(base) => check(
                    "eliminated_annihilates_source",
                    annihilation_statement,
                    Failed,
                    format!("(E A)_n != 0 at the exceptional base n = {base} (exact rows)"),
                ),
                None => check(
                    "eliminated_annihilates_source",
                    annihilation_statement,
                    Verified,
                    format!(
                        "g(n) E = sum_nu H_nu S^nu M is an exact Ore-algebra identity with \
                         g(n) = {}; M holds at every base n >= {} (the rows from index {} on \
                         are generated by the source rule), so (E A)_n = 0 wherever g(n) != 0; \
                         the remaining bases {:?} (earlier bases and roots of g) are checked \
                         exactly on rows",
                        poly_text(factor, "n"),
                        source.first_regular_base(),
                        n0 + source.seeds.len(),
                        certificate.exceptional_bases
                    ),
                ),
            }
        }
    });

    // 2. The coupling does not vanish.
    let coupling_statement = "epsilon_n != 0 for every n >= n0";
    let epsilon_roots = roots_in_domain(&parameters.epsilon, n0, None);
    identities.push(if vec_poly::is_zero(&parameters.epsilon) {
        check(
            "coupling_nonvanishing",
            coupling_statement,
            Unsupported,
            "epsilon is identically zero: the companion is not determined by A",
        )
    } else if let Some(root) = epsilon_roots.first() {
        check(
            "coupling_nonvanishing",
            coupling_statement,
            Unsupported,
            format!(
                "epsilon_n = {} vanishes at n = {root}; vanishing couplings are not supported \
                 (the undivided equations are not cancelled)",
                poly_text(&parameters.epsilon, "n")
            ),
        )
    } else {
        check(
            "coupling_nonvanishing",
            coupling_statement,
            Verified,
            format!(
                "epsilon_n = {} has no integer root n >= {n0} (exact)",
                poly_text(&parameters.epsilon, "n")
            ),
        )
    });

    // 3. R_n(0) = 0 for every n >= n0.
    let divisibility_statement = "R_n(0) = 0, i.e. A_(n+1)(0) = a_n A_n(0), for every n >= n0";
    let w0 = &rows[1][0] - evaluate_at_index(&parameters.a, n0) * &rows[0][0];
    // The induction needs (E A)_n = 0 and epsilon_n != 0 for every n.
    let annihilation_and_coupling = identities.iter().all(|c| c.status == Verified);
    identities.push(match ops.constant_term_projection() {
        None => check(
            "constant_term_divisibility",
            divisibility_statement,
            Unresolved,
            "a derivative term of E does not vanish at x = 0",
        ),
        Some(projection) if !ops.constant_term_factorization_holds(&projection) => check(
            "constant_term_divisibility",
            divisibility_statement,
            Unresolved,
            "the constant-term projection of E does not factor through S - a_n",
        ),
        Some(_) if !w0.is_zero() => check(
            "constant_term_divisibility",
            divisibility_statement,
            Failed,
            format!("A_(n0+1)(0) - a_(n0) A_(n0)(0) = {w0} != 0 at n0 = {n0}"),
        ),
        Some(_) => check(
            "constant_term_divisibility",
            divisibility_statement,
            if annihilation_and_coupling {
                Verified
            } else {
                Unresolved
            },
            "E at x = 0 equals (epsilon_n S - g_n) o (S - a_n) with \
             g_n = epsilon_(n+1) (a_n - epsilon_n) (exact operator identity), so \
             w_n = A_(n+1)(0) - a_n A_n(0) satisfies epsilon_n w_(n+1) = g_n w_n wherever \
             (E A)_n(0) = 0; with w_(n0) = 0 (exact) and epsilon_n != 0, w_n = 0 for all n \
             (unresolved while the two obligations above are not verified)",
        ),
    });

    // 4. The second coupled equation, by elimination.
    let second_statement =
        "B_(n+1) = delta_n A_n + V_n B_n for every n >= n0, with B_n = R_n / (epsilon_n x)";
    let prerequisites = identities.iter().all(|c| c.status == Verified);
    identities.push(if !ops.elimination_identity_holds() {
        check(
            "second_coupled_equation",
            second_statement,
            Failed,
            "the elimination identity fails as an operator identity",
        )
    } else if prerequisites {
        check(
            "second_coupled_equation",
            second_statement,
            Verified,
            "exact operator identity E A = eps_n E_(n+1) + eps_n eps_(n+1) x F_n - \
             eps_(n+1) Vt_n E_n (both components; the B component is Vt_n o x = x o V_n); \
             E_n = 0 by the definition of B_n, (E A)_n = 0 and eps_n eps_(n+1) != 0, so F_n = 0",
        )
    } else {
        check(
            "second_coupled_equation",
            second_statement,
            Unresolved,
            "the elimination identity holds as an operator identity, but it gives the second \
             equation only together with the obligations above",
        )
    });

    // 5. Symmetry for every n.
    let symmetry_statement = "A_n is symmetric at d_n and B_n at d_n - 1 for every n >= n0";
    let seeds_symmetric = is_symmetric_at(&rows[0], d0) && is_symmetric_at(&rows[1], d0 + 1);
    let coupled_system = identities.iter().all(|c| c.status == Verified);
    identities.push(if d0 == 0 {
        check(
            "symmetry_schedule",
            symmetry_statement,
            Unsupported,
            "d_(n0) = 0: the regular domain needs d_n >= 1; start the domain at the next index",
        )
    } else if !seeds_symmetric {
        check(
            "symmetry_schedule",
            symmetry_statement,
            Failed,
            format!(
                "the seeds A_(n0), A_(n0+1) are not symmetric at {d0}, {}",
                d0 + 1
            ),
        )
    } else {
        check(
            "symmetry_schedule",
            symmetry_statement,
            if coupled_system { Verified } else { Unresolved },
            "seeds symmetric at d_(n0), d_(n0) + 1 (exact); H(d, a, c) maps symmetric-at-d to \
             symmetric-at-(d+1) and x maps symmetric-at-d to symmetric-at-(d+2), so induction \
             over the coupled system carries the symmetry (degree deficits included; \
             unresolved while the coupled system is not established)",
        )
    });

    // 6. The parent collapse.
    let parent_statement = "P_n = A_n + x B_n satisfies P_(n+1) = (a_n + b_n x) P_n \
                            + c_n x (1 - x) D P_n for every n >= n0";
    let collapse = ops.collapse_identity_holds();
    let rows_divisible = !evaluate_at_index(&parameters.epsilon, n0).is_zero()
        && companion_row(parameters, &rows, n0).is_some();
    if rows_divisible {
        let companion = companion_row(parameters, &rows, n0).unwrap_or_default();
        if encode(&companion) != certificate.companion_seed {
            return Err(mismatch("companion seed B_(n0)"));
        }
        let parent_seed = vec_poly::trim(&vec_poly::add(&rows[0], &vec_poly::shift(&companion)));
        let parent = RecurrenceJson::from_recurrence_rational(
            &parent_recurrence(parameters),
            n0,
            &[parent_seed],
            None,
        );
        if parent != certificate.parent {
            return Err(mismatch("parent"));
        }
    }
    identities.push(if !collapse {
        check(
            "parent_collapse",
            parent_statement,
            Failed,
            "a collapse identity fails as an operator identity",
        )
    } else {
        check(
            "parent_collapse",
            parent_statement,
            if identities.iter().all(|c| c.status == Verified) {
                Verified
            } else {
                Unresolved
            },
            "exact operator identities L_n = U_n + delta_n x and L_n o x = epsilon_n x + x o V_n; \
             they give the parent rule from the two coupled equations",
        )
    });

    // 7. Exact rows (defensive self-check).
    let (self_checked_rows, failure) = row_self_check(parameters, &rows, certificate)?;
    identities.push(match failure {
        Some(reason) => check(
            "row_self_check",
            "the coupled system, symmetry and parent agree on exact rows",
            Failed,
            reason,
        ),
        None => check(
            "row_self_check",
            "the coupled system, symmetry and parent agree on exact rows",
            Verified,
            format!("{self_checked_rows} indices from n0 = {n0} checked exactly (not a proof)"),
        ),
    });

    let hypotheses = theorem_hypotheses(parameters, &rows);
    Ok(Outcome {
        identities,
        hypotheses,
        self_checked_rows,
    })
}

/// Check the coupled system on exact rows; `(indices checked, first failure)`.
fn row_self_check(
    parameters: &CompanionParameters,
    rows: &[Vec<Q>],
    certificate: &SymmetricCompanionCertificate,
) -> Result<(usize, Option<String>), RecurrenceError> {
    let n0 = parameters.first_index;
    let d0 = parameters.symmetry_degree_at_first_index;
    let ops = parameters.operators();
    let (parent_recurrence, parent_first, parent_seeds) =
        certificate.parent.to_recurrence_parts()?;
    let parents = parent_recurrence
        .generate_rows_rational(&parent_seeds, parent_first, rows.len() - 1)
        .map_err(|e| RecurrenceError::VerificationFailed(format!("parent replay failed: {e}")))?;
    // Stop before the first vanishing coupling (reported separately as
    // unsupported); a nondivisible remainder is a failure.
    let mut count = rows.len() - 2;
    let mut companions = Vec::new();
    for offset in 0..=count {
        let n = n0 + offset;
        if evaluate_at_index(&parameters.epsilon, n).is_zero() {
            count = offset.saturating_sub(1);
            break;
        }
        match companion_row(parameters, rows, n) {
            Some(companion) => companions.push(companion),
            None => return Ok((offset, Some(format!("R_n(0) != 0 at n = {n}")))),
        }
    }
    for offset in 0..count {
        let n = n0 + offset;
        let d = d0 + offset;
        if !is_symmetric_at(&rows[offset], d) {
            return Ok((offset, Some(format!("A_{n} is not symmetric at {d}"))));
        }
        if d >= 1 && !is_symmetric_at(&companions[offset], d - 1) {
            return Ok((offset, Some(format!("B_{n} is not symmetric at {}", d - 1))));
        }
        let delta = evaluate_at_index(&parameters.delta(), n);
        let expected = vec_poly::add(
            &vec_poly::scale(&rows[offset], &delta),
            &apply_to_row(&ops.v, &companions[offset], n),
        );
        if vec_poly::trim(&expected) != companions[offset + 1] {
            return Ok((
                offset,
                Some(format!("the second coupled equation fails at n = {n}")),
            ));
        }
        let collapsed = vec_poly::add(&rows[offset], &vec_poly::shift(&companions[offset]));
        if vec_poly::trim(&collapsed) != vec_poly::trim(&parents[offset]) {
            return Ok((offset, Some(format!("P_{n} != A_{n} + x B_{n}"))));
        }
    }
    Ok((count, None))
}

fn theorem_hypotheses(
    parameters: &CompanionParameters,
    rows: &[Vec<Q>],
) -> Vec<SymmetricCompanionCheck> {
    use SymmetricCompanionCheckStatus::*;
    let n0 = parameters.first_index;
    let mut hypotheses = Vec::new();
    for (name, statement, p) in [
        (
            "a_nonnegative",
            "a_n >= 0 for every n >= n0",
            parameters.a.clone(),
        ),
        (
            "delta_nonnegative",
            "delta_n >= 0 for every n >= n0",
            parameters.delta(),
        ),
        (
            "delta_at_most_c",
            "delta_n <= c_n (epsilon_n >= 0) for every n >= n0",
            parameters.epsilon.clone(),
        ),
    ] {
        let (status, evidence) = nonnegative_on_domain(&p, n0);
        hypotheses.push(check(name, statement, status, evidence));
    }
    let d0 = parameters.symmetry_degree_at_first_index;
    hypotheses.push(check(
        "symmetry_index_positive",
        "d_n >= 1 for every n >= n0",
        if d0 >= 1 { Verified } else { Failed },
        format!("d_(n0) = {d0} and d_n increases by one per step"),
    ));
    let companion = companion_row(parameters, rows, n0);
    match companion {
        None => {
            hypotheses.push(check(
                "initial_pair_nonnegative",
                "A_(n0) and B_(n0) have nonnegative coefficients",
                Unresolved,
                "B_(n0) is not defined",
            ));
            hypotheses.push(check(
                "initial_pair_interlacing",
                "B_(n0) <= A_(n0) (directed weak interlacing)",
                Unresolved,
                "B_(n0) is not defined",
            ));
        }
        Some(companion) => {
            let nonnegative = vec_poly::has_nonnegative_coefficients(&rows[0])
                && vec_poly::has_nonnegative_coefficients(&companion);
            hypotheses.push(check(
                "initial_pair_nonnegative",
                "A_(n0) and B_(n0) have nonnegative coefficients",
                if nonnegative { Verified } else { Failed },
                format!(
                    "A_(n0) = {}, B_(n0) = {}",
                    poly_text(&rows[0], "x"),
                    poly_text(&companion, "x")
                ),
            ));
            let (interlacing, evidence) = companion_interlaces(&companion, &rows[0]);
            hypotheses.push(check(
                "initial_pair_interlacing",
                "B_(n0) <= A_(n0) (directed weak interlacing)",
                if interlacing { Verified } else { Failed },
                evidence,
            ));
        }
    }
    hypotheses
}

fn assemble(outcome: Outcome) -> SymmetricCompanionVerification {
    use SymmetricCompanionCheckStatus::*;
    let statuses = || outcome.identities.iter().map(|c| c.status);
    let identity_status = if statuses().any(|s| s == Failed) {
        SymmetricCompanionIdentityStatus::Rejected
    } else if statuses().any(|s| s == Unsupported) {
        SymmetricCompanionIdentityStatus::Unsupported
    } else if statuses().any(|s| s == Unresolved) {
        SymmetricCompanionIdentityStatus::Candidate
    } else {
        SymmetricCompanionIdentityStatus::Certified
    };
    let hypotheses_hold = outcome.hypotheses.iter().all(|c| c.status == Verified);
    let real_rootedness_status = match identity_status {
        SymmetricCompanionIdentityStatus::Certified if hypotheses_hold => {
            SymmetricCompanionRealRootednessStatus::TheoremApplies
        }
        SymmetricCompanionIdentityStatus::Certified => {
            SymmetricCompanionRealRootednessStatus::HypothesesNotVerified
        }
        _ => SymmetricCompanionRealRootednessStatus::NotEvaluated,
    };
    let mut unresolved_obligations = Vec::new();
    for item in &outcome.identities {
        if item.status != Verified {
            unresolved_obligations.push(format!("{}: {}", item.name, item.statement));
        }
    }
    if identity_status == SymmetricCompanionIdentityStatus::Certified {
        for item in &outcome.hypotheses {
            if item.status != Verified {
                unresolved_obligations.push(format!("{}: {}", item.name, item.statement));
            }
        }
    }
    let conclusion = match real_rootedness_status {
        SymmetricCompanionRealRootednessStatus::TheoremApplies => {
            "identity certified for every n >= n0 and every hypothesis of the preservation \
             theorem holds: every A_n, B_n and P_n (n >= n0) is zero or real-rooted with \
             nonnegative coefficients and B_n <= A_n; Lean replay not attempted"
        }
        SymmetricCompanionRealRootednessStatus::HypothesesNotVerified => {
            "identity certified for every n >= n0, but a hypothesis of the preservation theorem \
             is not verified: this route gives no real-rootedness proof (which says nothing \
             against real-rootedness of A_n)"
        }
        SymmetricCompanionRealRootednessStatus::NotEvaluated => match identity_status {
            SymmetricCompanionIdentityStatus::Rejected => {
                "an obligation fails exactly: this is not a symmetric companion of A_n"
            }
            SymmetricCompanionIdentityStatus::Unsupported => {
                "the candidate leaves the supported class; nothing is claimed"
            }
            _ => "candidate only: some identity obligation is unresolved; nothing is claimed",
        },
    }
    .to_string();
    SymmetricCompanionVerification {
        identity_status,
        real_rootedness_status,
        lean_replay: SYMMETRIC_COMPANION_LEAN_REPLAY_NOT_ATTEMPTED.into(),
        identities: outcome.identities,
        hypotheses: outcome.hypotheses,
        unresolved_obligations,
        self_checked_rows: outcome.self_checked_rows,
        conclusion,
    }
}
