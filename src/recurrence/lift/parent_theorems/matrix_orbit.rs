//! Stationary matrix orbits preserving ordered interlacing vectors.
//!
//! Affine entries are decided by the quantified pencil criterion, including
//! repeated indices and zero coordinates. Scalar recurrence identities are
//! checked as matrix polynomial identities, rather than fitted on rows.

use super::*;
use crate::Polynomial;
use std::collections::BTreeSet;

pub(in crate::recurrence) mod affine;
mod realization;
#[cfg(test)]
mod tests;

type Poly = Polynomial<Q>;
type Matrix = Vec<Vec<Poly>>;
type RationalMatrix = Vec<Vec<Q>>;

/// Version of the matrix-orbit request/report protocol.
pub const MATRIX_ORBIT_SCHEMA: &str = "polytool.matrix-orbit.v1";
const REFERENCE: &str = "Branden, Unimodality, log-concavity, real-rootedness and beyond, \
    Theorem 8.5 and Lemma 8.4 (arXiv:1410.6601); affine quantified pencils";
const MAXIMUM_DIMENSION: usize = 5;
const MAXIMUM_BASIS_CANDIDATES: usize = 4096;

/// A stationary state system, with ascending rational-string coefficients.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MatrixOrbitStateSystem {
    /// Square polynomial matrix, with dimension between one and five.
    pub matrix: Vec<Vec<Vec<String>>>,
    /// Polynomial column vector at the recurrence record's first index.
    pub start_vector: Vec<Vec<String>>,
    /// Constant rational row vector selecting the scalar target.
    pub projection: Vec<String>,
}

/// A recurrence and optional exact stationary state realization.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MatrixOrbitRequest {
    /// Defining scalar recurrence, index and every prescribed seed.
    pub recurrence: RecurrenceJson,
    /// Supplied realization; absence selects bounded automatic recognition.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub state_system: Option<MatrixOrbitStateSystem>,
}

/// Exact stationary points checked to decide quadratic copositivity.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MatrixOrbitCopositivityCheck {
    /// Symmetric rational matrix, in row order.
    pub matrix: Vec<Vec<String>>,
    /// Minimum on the nonnegative unit simplex, as a rational string.
    pub minimum: String,
}

/// Universal affine pencil check for one weakly ordered pair of index pairs.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MatrixOrbitPencilCheck {
    /// Zero-based row indices, in weakly increasing order.
    pub rows: [usize; 2],
    /// Zero-based column indices, in weakly increasing order.
    pub columns: [usize; 2],
    /// True when one of the two pencils is identically zero.
    pub zero_side: bool,
    /// Rational constants A, C, K, D0, D1, D2 of the Wronskian reduction.
    pub wronskian_parameters: Vec<String>,
    /// The dimension-two and dimension-three copositivity checks, when nonzero.
    pub copositivity: Vec<MatrixOrbitCopositivityCheck>,
}

/// A successful state realization and its exact all-index proof data.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MatrixOrbitCertificate {
    /// Original state realization, before reordering and changing basis.
    pub state_system: MatrixOrbitStateSystem,
    /// Coordinate order: new position `i` takes old coordinate `permutation[i]`.
    pub permutation: Vec<usize>,
    /// Nonnegative rational unit triangular basis, after reordering.
    pub basis: Vec<Vec<String>>,
    /// Transformed matrix, seed and projection, checked by exact identities.
    pub transformed_system: MatrixOrbitStateSystem,
    /// Exact checks for every weakly ordered row/column index pair.
    pub pencils: Vec<MatrixOrbitPencilCheck>,
    /// Boolean support masks of successive orbit vectors until the first repeat.
    /// Evaluation at one identifies nonzero nonnegative polynomials exactly.
    pub nonzero_supports: Vec<usize>,
    /// Provenance of the realization shape, without sequence identifiers.
    pub realization: String,
    /// Matching Lean preservation declarations, when their shape applies.
    pub lean_declarations: Vec<String>,
}

/// Replayable matrix-orbit report. A negative result is bounded search failure.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MatrixOrbitReport {
    /// [`MATRIX_ORBIT_SCHEMA`].
    pub schema: String,
    /// Embedded recurrence and optional supplied state system.
    pub request: MatrixOrbitRequest,
    /// `proved` or `not_proved`; no negative all-index assertion is inferred.
    pub status: String,
    /// Exact theorem reference.
    pub reference: String,
    /// Description of the proved conclusion or the search obstruction.
    pub reason: String,
    /// Successful certificate, present exactly when the status is `proved`.
    pub certificate: Option<MatrixOrbitCertificate>,
    /// Lean replay status, always `not_attempted`.
    pub lean_replay: String,
}

fn invalid(message: impl Into<String>) -> RecurrenceError {
    RecurrenceError::InvalidInput(message.into())
}
fn parse_numbers(values: &[String]) -> Result<Vec<Q>, RecurrenceError> {
    values
        .iter()
        .map(|x| {
            x.parse::<Q>()
                .map_err(|_| invalid("invalid rational coefficient"))
        })
        .collect()
}
fn strings(values: &[Q]) -> Vec<String> {
    values.iter().map(ToString::to_string).collect()
}
fn encode(matrix: &Matrix, seed: &[Poly], projection: &[Q]) -> MatrixOrbitStateSystem {
    MatrixOrbitStateSystem {
        matrix: matrix
            .iter()
            .map(|r| {
                r.iter()
                    .map(|p| strings(&vec_poly::trim(p.coeffs())))
                    .collect()
            })
            .collect(),
        start_vector: seed
            .iter()
            .map(|p| strings(&vec_poly::trim(p.coeffs())))
            .collect(),
        projection: strings(projection),
    }
}
fn decode(system: &MatrixOrbitStateSystem) -> Result<(Matrix, Vec<Poly>, Vec<Q>), RecurrenceError> {
    let n = system.matrix.len();
    if n == 0
        || n > MAXIMUM_DIMENSION
        || system.matrix.iter().any(|r| r.len() != n)
        || system.start_vector.len() != n
        || system.projection.len() != n
    {
        return Err(invalid(
            "state system must be square, of dimension 1..=5, with matching seed and projection",
        ));
    }
    let matrix = system
        .matrix
        .iter()
        .map(|r| r.iter().map(|p| parse_numbers(p).map(Poly::new)).collect())
        .collect::<Result<Matrix, _>>()?;
    let seed = system
        .start_vector
        .iter()
        .map(|p| parse_numbers(p).map(Poly::new))
        .collect::<Result<Vec<_>, _>>()?;
    if matrix
        .iter()
        .flatten()
        .chain(seed.iter())
        .any(|p| p.degree().is_some_and(|d| d > 24))
    {
        return Err(invalid("state polynomial degree exceeds 24"));
    }
    Ok((matrix, seed, parse_numbers(&system.projection)?))
}
fn identity(n: usize) -> Matrix {
    (0..n)
        .map(|i| {
            (0..n)
                .map(|j| Poly::new(vec![Q::from_integer((i == j).into())]))
                .collect()
        })
        .collect()
}
fn multiply(a: &Matrix, b: &Matrix) -> Matrix {
    let n = a.len();
    (0..n)
        .map(|i| {
            (0..n)
                .map(|j| (0..n).fold(Poly::zero(), |v, k| v + a[i][k].clone() * b[k][j].clone()))
                .collect()
        })
        .collect()
}
fn action(matrix: &Matrix, seed: &[Poly]) -> Vec<Poly> {
    matrix
        .iter()
        .map(|r| {
            r.iter()
                .zip(seed)
                .fold(Poly::zero(), |s, (a, b)| s + a.clone() * b.clone())
        })
        .collect()
}
fn project(projection: &[Q], seed: &[Poly]) -> Poly {
    projection.iter().zip(seed).fold(Poly::zero(), |v, (a, b)| {
        v + Poly::new(vec![a.clone()]) * b.clone()
    })
}
fn nonnegative(p: &Poly) -> bool {
    p.coeffs().iter().all(|c| !c.is_negative())
}
fn interlacing_seed(seed: &[Poly]) -> bool {
    let rows: Vec<Vec<BigInt>> = seed
        .iter()
        .map(|p| integer_coefficients(p.coeffs()))
        .collect();
    seed.iter().all(nonnegative)
        && rows.iter().all(|p| crate::is_real_rooted_bigint_coeffs(p))
        && (0..seed.len()).all(|i| {
            (i + 1..seed.len()).all(|j| {
                seed[i].is_zero()
                    || seed[j].is_zero()
                    || crate::check_weak_interlacing_bigint_coeffs(&rows[i], &rows[j]) == Some(true)
            })
        })
}
fn nonzero_supports(matrix: &Matrix, seed: &[Poly], projection: &[Q]) -> Option<Vec<usize>> {
    let mut state = seed
        .iter()
        .enumerate()
        .fold(0, |s, (i, p)| s | (usize::from(!p.is_zero()) << i));
    let target = projection
        .iter()
        .enumerate()
        .fold(0, |s, (i, c)| s | (usize::from(!c.is_zero()) << i));
    let mut seen = BTreeSet::new();
    let mut supports = Vec::new();
    while seen.insert(state) {
        if state & target == 0 {
            return None;
        }
        supports.push(state);
        state = matrix.iter().enumerate().fold(0, |s, (i, r)| {
            s | (usize::from(
                r.iter()
                    .enumerate()
                    .any(|(j, p)| state & (1 << j) != 0 && !p.is_zero()),
            ) << i)
        });
    }
    Some(supports)
}
fn identity_and_seeds(
    matrix: &Matrix,
    seed: &[Poly],
    projection: &[Q],
    multipliers: &[Poly],
    rows: &[Vec<Q>],
) -> bool {
    let n = matrix.len();
    let mut power = identity(n);
    let mut rhs = vec![vec![Poly::zero(); n]; n];
    let order = multipliers.len();
    for k in 0..order {
        let coefficient = &multipliers[order - k - 1];
        for i in 0..n {
            for j in 0..n {
                rhs[i][j] = rhs[i][j].clone() + coefficient.clone() * power[i][j].clone();
            }
        }
        power = multiply(matrix, &power);
    }
    if power != rhs {
        return false;
    }
    let mut state = seed.to_vec();
    for row in rows {
        if project(projection, &state) != Poly::new(row.clone()) {
            return false;
        }
        state = action(matrix, &state);
    }
    true
}

/// Decide the supplied recurrence by bounded stationary affine-matrix recognition.
///
/// Supplied systems use every state ordering and a finite menu of nonnegative
/// triangular bases; automatic recognition uses companion forms and threshold
/// and separated-characteristic shapes. A proof checks a matrix annihilating
/// identity, every prescribed scalar seed, the transformed ordered PF seed,
/// nonnegative projection, all quantified pencils, and exact nonvanishing.
/// Only real-rootedness of projected rows is claimed, without temporal interlacing.
///
/// # Errors
/// Returns [`RecurrenceError`] for malformed rational data, dimensions or records.
pub fn check_matrix_orbit(
    request: &MatrixOrbitRequest,
) -> Result<MatrixOrbitReport, RecurrenceError> {
    let (rule, _, rows) = request.recurrence.to_recurrence_parts()?;
    let supplied = request.state_system.as_ref().map(decode).transpose()?;
    let mut report = MatrixOrbitReport {
        schema: MATRIX_ORBIT_SCHEMA.into(),
        request: request.clone(),
        status: "not_proved".into(),
        reference: REFERENCE.into(),
        reason: "no admissible realization within the bounded search".into(),
        certificate: None,
        lean_replay: "not_attempted".into(),
    };
    if rule.terms.iter().any(|term| term.offset > 8) {
        report.reason = "recurrence order exceeds eight".into();
        return Ok(report);
    }
    let Ok(multipliers) = degree_law_reversal::constant_coefficient_multipliers(&rule) else {
        report.reason = "requires a stationary homogeneous recurrence without derivatives".into();
        return Ok(report);
    };
    if multipliers.len() > 8 || rows.len() < multipliers.len() || rows.len() > 32 {
        report.reason =
            "requires order at most eight, enough seeds, and at most 32 prescribed rows".into();
        return Ok(report);
    }
    let certificate = if let Some((matrix, seed, projection)) = supplied {
        if !identity_and_seeds(&matrix, &seed, &projection, &multipliers, &rows) {
            report.reason =
                "supplied system fails the matrix recurrence identity or a prescribed scalar seed"
                    .into();
            return Ok(report);
        }
        let mut remaining = MAXIMUM_BASIS_CANDIDATES;
        let declarations = realization::lean_declarations(&matrix, "supplied");
        realization::search_bases(
            &matrix,
            &seed,
            &projection,
            "supplied",
            &declarations,
            &mut remaining,
        )
    } else {
        realization::recognize(&multipliers, &rows)
    };
    if let Some(certificate) = certificate {
        report.status = "proved".into();
        report.reason = "every projected row from the first index is nonzero and real-rooted; no consecutive interlacing or simplicity is asserted".into();
        report.certificate = Some(certificate);
    }
    Ok(report)
}

/// Replay a matrix-orbit report from its recorded realization, without searching.
///
/// # Errors
/// Rejects altered identities, seeds, projections, pencil data or conclusions.
pub fn verify_matrix_orbit(report: &MatrixOrbitReport) -> Result<(), RecurrenceError> {
    let fail =
        || RecurrenceError::VerificationFailed("matrix-orbit certificate does not replay".into());
    if report.schema != MATRIX_ORBIT_SCHEMA {
        return Err(fail());
    }
    if let Some(c) = &report.certificate {
        let mut request = report.request.clone();
        request.state_system = Some(c.state_system.clone());
        let (rule, _, rows) = request.recurrence.to_recurrence_parts()?;
        if rule.terms.iter().any(|term| term.offset > 8) {
            return Err(fail());
        }
        let multipliers =
            degree_law_reversal::constant_coefficient_multipliers(&rule).map_err(|_| fail())?;
        let (matrix, seed, projection) = decode(&c.state_system)?;
        if multipliers.len() > 8
            || rows.len() < multipliers.len()
            || rows.len() > 32
            || !realization::valid_provenance(
                &matrix,
                &c.realization,
                report.request.state_system.is_some(),
            )
        {
            return Err(fail());
        }
        if !identity_and_seeds(&matrix, &seed, &projection, &multipliers, &rows) {
            return Err(fail());
        }
        if let Some(system) = &report.request.state_system {
            if decode(system)? != (matrix.clone(), seed.clone(), projection.clone()) {
                return Err(fail());
            }
        }
        let basis = c
            .basis
            .iter()
            .map(|r| parse_numbers(r))
            .collect::<Result<Vec<_>, _>>()?;
        let expected = realization::certify_basis(
            &matrix,
            &seed,
            &projection,
            &c.permutation,
            &basis,
            &c.realization,
            &c.lean_declarations,
        )
        .ok_or_else(fail)?;
        if &expected != c || report.status!="proved" || report.reference!=REFERENCE
            || report.lean_replay!="not_attempted" || report.reason!="every projected row from the first index is nonzero and real-rooted; no consecutive interlacing or simplicity is asserted" {
            return Err(fail());
        }
        // Metadata is derived from the actual original matrix, never trusted.
        if c.lean_declarations != realization::lean_declarations(&matrix, &c.realization) {
            return Err(fail());
        }
    } else if check_matrix_orbit(&report.request)? != *report {
        return Err(fail());
    }
    Ok(())
}

/// Prove a scalar recurrence by automatic interlacing-matrix orbit recognition.
///
/// # Errors
/// Returns [`RecurrenceError`] for an invalid recurrence record or seed data.
pub fn check_matrix_orbit_theorem(
    rule: &Recurrence,
    first: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    check_real_rootedness_theorem(
        RealRootednessTheorem::MatrixOrbit,
        rule,
        first,
        seeds,
        options,
    )
}
pub(super) const ROUTE: RoutePolicy = RoutePolicy {
    reference: REFERENCE,
    phase: DiscoveryPhase::Direct,
    after_direct_eligible: None,
    replay: ParentReplayPolicy::RecordedRoute,
    discovery_replay: DiscoveryReplayPolicy::Direct,
    check_at_supplied_seeds: theorem_report,
    start_policy: StartPolicy::OwnStarts,
    is_preferred: any_proof,
    certificate_name: Some("matrix-orbit"),
    discovery: DiscoveryPolicy::WhenStronger,
    in_residue_class_menu: false,
    in_index_section_menu: false,
    in_degree_law_reversal_menu: false,
    parent: ParentPolicy {
        stage: ParentStage::WhileUnproved,
        // The PF-parent menu requests a root ray. This splitting-only route
        // declines root-region requests rather than certifying an open ray.
        region: ParentRegion::OpenNonpositiveRay,
        record_only_if_applicable: true,
    },
};
fn theorem_report(
    rule: &Recurrence,
    first: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    let mut report = RealRootednessTheoremReport::new(
        RealRootednessTheorem::MatrixOrbit,
        rule,
        first,
        seeds,
        options,
    );
    if options.root_region.is_some() {
        return Ok(
            report.not_applicable("matrix orbit does not certify requested open root regions")
        );
    }
    let checked = check_matrix_orbit(&MatrixOrbitRequest {
        recurrence: report.recurrence.clone(),
        state_system: None,
    })?;
    let Some(c) = checked.certificate.as_ref() else {
        return Ok(report.not_applicable(checked.reason));
    };
    report.normal_form = format!(
        "P_(first+n) = projection * M(x)^n * start; {}",
        c.realization
    );
    report.lean_declarations = c.lean_declarations.clone();
    report.hypotheses.push(exact_hypothesis(
        "matrix_orbit_identity",
        "matrix annihilation and every prescribed scalar seed agree exactly",
        Ok("identity over Q[x] and recurrence uniqueness".into()),
    ));
    report.hypotheses.push(exact_hypothesis("interlacing_matrix", "nonnegative affine entries satisfy every quantified pencil, the seed is ordered PF, and projection weights are nonnegative", Ok(format!("{} exact universal pencil checks",c.pencils.len()))));
    report.hypotheses.push(exact_hypothesis(
        "nonzero_projection",
        "every projected row is nonzero",
        Ok(format!(
            "{} exact support states until repetition",
            c.nonzero_supports.len()
        )),
    ));
    report.certificate = Some(RouteCertificate::MatrixOrbit(Box::new(checked)));
    Ok(report.finish(RealRootednessTheoremConclusion {
        first_index: first, zero_residue_classes: None, interlacing_from: None, interlacing: None,
        simple_roots: false, root_region: None,
        statement: format!("every P_n, n >= {first}, is nonzero and real-rooted by a stationary interlacing-matrix orbit; no consecutive interlacing or simplicity is asserted"),
    }))
}
pub(super) fn verify_certificate_in_report(
    c: &MatrixOrbitReport,
    report: &RealRootednessTheoremReport,
) -> Result<(), RecurrenceError> {
    let (mut rule, first, seeds) = report.recurrence.to_recurrence_parts()?;
    if let Some(normalization) = &report.left_hand_factor_normalization {
        normalization.verify()?;
        rule = normalization.normalized_recurrence()?;
    }
    let expected = RecurrenceJson::from_recurrence_rational(&rule, first, &seeds, None);
    if c.request.recurrence != expected || c.request.state_system.is_some() || c.status != "proved"
    {
        return Err(RecurrenceError::VerificationFailed(
            "matrix-orbit recurrence mismatch".into(),
        ));
    }
    verify_matrix_orbit(c)
}
