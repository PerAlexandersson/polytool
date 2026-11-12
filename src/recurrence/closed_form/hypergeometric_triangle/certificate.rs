//! Hypergeometric operator certificates, replay and derivation.

use super::*;

/// Prove `E P_n = 0` for every `n >= n_0` (module documentation) and build
/// the recurrence record of `E`.
///
/// The operator must have a single term at its largest shift `I >= 1`, with
/// `variable_power = euler_power = 0` (a leader free of `t` and `theta`).
///
/// # Errors
///
/// [`RecurrenceError::InvalidInput`] for a malformed triangle or operator,
/// [`RecurrenceError::Unsupported`] when the definition cannot be checked
/// (a ratio factor that may vanish on the support, rows too narrow), and
/// [`RecurrenceError::VerificationFailed`] when an identity fails (with the
/// failing block or row).
pub fn certify_hypergeometric_row_operator(
    triangle: &HypergeometricTriangle,
    operator: &[EulerOperatorTerm],
) -> Result<HypergeometricRowRecurrenceCertificate, RecurrenceError> {
    certify_with_limit(
        triangle,
        operator,
        HypergeometricRowRecurrenceBudget::default().max_checked_row,
    )
}

fn certify_with_limit(
    triangle: &HypergeometricTriangle,
    operator: &[EulerOperatorTerm],
    max_checked_row: usize,
) -> Result<HypergeometricRowRecurrenceCertificate, RecurrenceError> {
    let terms = triangle.parse()?;
    let n0 = triangle.first_index()?;
    let parsed = parse_operator(operator)?;
    let order = parsed
        .keys()
        .map(|(i, _, _)| *i)
        .max()
        .ok_or_else(|| RecurrenceError::InvalidInput("the operator is zero".into()))?;
    let leader_keys: Vec<_> = parsed.keys().filter(|(i, _, _)| *i == order).collect();
    if order == 0 || leader_keys != vec![&(order, 0, 0)] {
        return Err(RecurrenceError::InvalidInput(
            "the operator needs order >= 1 and a leader c(n) S^I free of t and theta".into(),
        ));
    }
    // Rows: before the separation index, and seeds past the leader's roots.
    let leader = parsed[&(order, 0, 0)].clone();
    let roots = crate::recurrence::implicit::algebra::roots_in_domain(&leader, n0, None);
    let seeds = order.max(roots.last().map_or(0, |r| r + order + 1 - n0));
    let proof = prove_annihilation(&terms, n0, n0, &parsed, seeds, max_checked_row)?;
    let derivative = derivative_operator(&parsed);
    let rows = proof.rows;
    let strips = proof.strips;
    let checked = proof.checked_bases;
    let (rule, shift) = ImplicitRelation::from_op(&derivative)
        .to_successor_rule()
        .ok_or_else(|| {
            RecurrenceError::VerificationFailed("the operator is not runnable".into())
        })?;
    debug_assert_eq!(shift, order);
    let record = RecurrenceJson::from_recurrence_rational(&rule, n0, &rows[..seeds], None);
    let text = operator_text(&parsed);
    let statement = format!(
        "For every n >= {n0}: {text} with P_n(t) = sum_k T(n,k) t^k (theta = t d/dt); the \
         record with the seed rows {n0}..{} generates exactly the rows P_n.",
        n0 + seeds - 1
    );
    Ok(HypergeometricRowRecurrenceCertificate {
        schema: HYPERGEOMETRIC_ROW_RECURRENCE_SCHEMA.into(),
        triangle: triangle.clone(),
        operator: operator_terms(&parsed),
        operator_text: text,
        boundary_strips: strips,
        checked_bases: checked,
        record,
        statement,
    })
}

/// JSON schema tag of [`HypergeometricRecordCertificate`].
pub const HYPERGEOMETRIC_RECORD_SCHEMA: &str = "polytool.hypergeometric-record.v1";

/// A proof that a recurrence record generates the rows of a hypergeometric
/// triangle for every index ([`certify_hypergeometric_row_record`]),
/// schema `polytool.hypergeometric-record.v1`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HypergeometricRecordCertificate {
    /// [`HYPERGEOMETRIC_RECORD_SCHEMA`].
    pub schema: String,
    /// The definition.
    pub triangle: HypergeometricTriangle,
    /// The record (search metadata and presentation are ignored).
    pub record: RecurrenceJson,
    /// The rule in forward Euler form at base `n = m - R` (`R` the largest
    /// lag of the record), multiplied on the left by
    /// `t^variable_power_multiplier`, which clears
    /// `D^d = t^(-d) theta (theta - 1) ... (theta - d + 1)`.
    pub operator: Vec<EulerOperatorTerm>,
    /// The power of `t` the rule was multiplied by.
    pub variable_power_multiplier: usize,
    /// The proved strip identities and separation indices.
    pub boundary_strips: Vec<HypergeometricBoundaryStrips>,
    /// Bases at which the rule was checked exactly on rows.
    pub checked_bases: Vec<usize>,
    /// What the certificate establishes.
    pub statement: String,
}

impl HypergeometricRecordCertificate {
    /// Replay: recompute the certificate from the recorded triangle and
    /// record and compare.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::VerificationFailed`] when the recomputation fails
    /// or differs.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        if self.schema != HYPERGEOMETRIC_RECORD_SCHEMA {
            return Err(RecurrenceError::VerificationFailed(format!(
                "unsupported schema `{}`",
                self.schema
            )));
        }
        let replay = certify_hypergeometric_row_record(&self.triangle, &self.record)
            .map_err(|error| RecurrenceError::VerificationFailed(error.to_string()))?;
        if &replay != self {
            return Err(RecurrenceError::VerificationFailed(
                "the hypergeometric record certificate differs from its recomputation".into(),
            ));
        }
        Ok(())
    }
}

/// Prove that `record` generates the rows `P_m` of `triangle` for every
/// `m >= record.first_index` (which must be at least the triangle's first
/// index).
///
/// The rule `f(m, t) P_m = sum c(m, t) D^d P_(m-r)` is written in forward
/// form at base `n = m - R` and multiplied by `t^J` (`J` its largest
/// derivative order), which gives an operator in `t` and `theta`; the
/// module's argument proves that it annihilates the triangle's rows at
/// every base `n >= record.first_index + (seed count) - R` (and the
/// identity transfers back, since `t^J` is not a zero divisor).  The seeds
/// are compared with the rows exactly, and the left-hand factor `f(m, .)`
/// is a nonzero polynomial at every row the record computes (its content
/// in `m` has no integer root there), so by induction the record's rows are
/// the triangle's rows.
///
/// # Errors
///
/// As [`certify_hypergeometric_row_operator`]; a seed that differs from the
/// triangle's row, a vanishing left-hand factor, an inhomogeneous or
/// sign-alternating rule ([`RecurrenceError::Unsupported`]).
pub fn certify_hypergeometric_row_record(
    triangle: &HypergeometricTriangle,
    record: &RecurrenceJson,
) -> Result<HypergeometricRecordCertificate, RecurrenceError> {
    let terms = triangle.parse()?;
    let n0 = triangle.first_index()?;
    let (rule, first, seeds) = record.to_recurrence_parts()?;
    if first < n0 {
        return Err(RecurrenceError::InvalidInput(format!(
            "the record starts at {first}, before the triangle's first index {n0}"
        )));
    }
    let (relation, lag) = ImplicitRelation::from_recurrence(&rule)
        .map_err(|error| RecurrenceError::Unsupported(error.to_string()))?;
    let op = relation.op();
    let computed_from = first + seeds.len();
    if computed_from < first + lag {
        return Err(RecurrenceError::InvalidInput(
            "the record has fewer seeds than its largest lag".into(),
        ));
    }
    let first_base = computed_from - lag;
    // The left-hand factor at base n is the coefficient of S^R.
    let leader =
        op.0.get(&(lag, 0))
            .cloned()
            .ok_or_else(|| RecurrenceError::InvalidInput("the rule has no leader".into()))?;
    let (content, _) = crate::recurrence::implicit::algebra::content_and_primitive(&leader);
    let roots = crate::recurrence::implicit::algebra::roots_in_domain(&content, first_base, None);
    if let Some(root) = roots.first() {
        return Err(RecurrenceError::Unsupported(format!(
            "the left-hand factor vanishes at the row {}",
            root + lag
        )));
    }
    let (multiplier, theta) = crate::recurrence::implicit::algebra::to_theta(&op);
    let mut parsed: EulerMap = BTreeMap::new();
    for ((i, j), c) in &theta {
        let width = c.variable_degree().map_or(0, |d| d + 1);
        for l in 0..width {
            let column = vec_poly::trim(&c.variable_coefficient(l));
            if !vec_poly::is_zero(&column) {
                parsed.insert((*i, *j, l), column);
            }
        }
    }
    let rows_needed = computed_from - n0;
    let proof = prove_annihilation(&terms, n0, first_base, &parsed, rows_needed, 400)?;
    for (offset, seed) in seeds.iter().enumerate() {
        let m = first + offset;
        if vec_poly::trim(seed) != proof.rows[m - n0] {
            return Err(RecurrenceError::VerificationFailed(format!(
                "the seed row {m} differs from the triangle's row"
            )));
        }
    }
    let text = operator_text(&parsed);
    let statement = format!(
        "For every m >= {first} the record's row equals P_m(t) = sum_k T(n,k) t^k: the seeds \
         {first}..{} are rows of the triangle, t^{multiplier} times the rule is {text} at \
         base n = m - {lag}, proved for every n >= {first_base}, and the left-hand factor does \
         not vanish at a computed row.",
        computed_from - 1
    );
    let mut record = record.clone();
    record.search = None;
    record.latex = None;
    Ok(HypergeometricRecordCertificate {
        schema: HYPERGEOMETRIC_RECORD_SCHEMA.into(),
        triangle: triangle.clone(),
        record,
        operator: operator_terms(&parsed),
        variable_power_multiplier: multiplier,
        boundary_strips: proof.strips,
        checked_bases: proof.checked_bases,
        statement,
    })
}

/// The proof of `E P_n = 0` for every `n >= first_base` (module
/// documentation, steps 2 to 4) and the rows it evaluated.
struct AnnihilationProof {
    strips: Vec<HypergeometricBoundaryStrips>,
    checked_bases: Vec<usize>,
    /// The rows `P_(n_0), P_(n_0 + 1), ...` (at least `minimum_rows`).
    rows: Vec<Vec<Q>>,
}

fn prove_annihilation(
    terms: &[Term],
    n0: usize,
    first_base: usize,
    parsed: &EulerMap,
    minimum_rows: usize,
    max_checked_row: usize,
) -> Result<AnnihilationProof, RecurrenceError> {
    let order = parsed.keys().map(|(i, _, _)| *i).max().unwrap_or(0);
    let pairs: BTreeSet<(usize, usize)> = parsed.keys().map(|(i, j, _)| (*i, *j)).collect();
    let mut strips = Vec::new();
    let mut separation = first_base;
    for (position, term) in terms.iter().enumerate() {
        for residue in 0..term.modulus {
            let (blocks, mut info) = blocks_for(term, residue, &pairs)?;
            for (b, block) in blocks.iter().enumerate() {
                if !block_value(block, parsed).is_zero() {
                    let what = if b == 0 {
                        "the interior identity".to_string()
                    } else if b <= info.lower_strips {
                        format!("the boundary strip k = lower(n) + {}", b - 1)
                    } else {
                        format!(
                            "the boundary strip k = upper(n) + {}",
                            b - info.lower_strips
                        )
                    };
                    return Err(RecurrenceError::VerificationFailed(format!(
                        "term {position}, residue {residue} mod {}: {what} fails",
                        term.modulus
                    )));
                }
            }
            info.term = position;
            separation = separation.max(info.separation_index);
            strips.push(info);
        }
    }
    let needed = (separation - n0 + order).max(minimum_rows);
    if n0 + needed > max_checked_row {
        return Err(RecurrenceError::BudgetExhausted(format!(
            "the certificate needs rows up to {} (limit {max_checked_row})",
            n0 + needed
        )));
    }
    let rows = total_rows(terms, needed)?;
    let derivative = derivative_operator(parsed);
    let mut checked = Vec::new();
    for n in first_base..separation {
        let value = derivative.apply(&rows, n0, n).ok_or_else(|| {
            RecurrenceError::VerificationFailed(format!("row check at base {n}: rows missing"))
        })?;
        if !vec_poly::is_zero(&value) {
            return Err(RecurrenceError::VerificationFailed(format!(
                "E P_n != 0 at the base n = {n} (exact rows)"
            )));
        }
        checked.push(n);
    }
    Ok(AnnihilationProof {
        strips,
        checked_bases: checked,
        rows,
    })
}

// ---------------------------------------------------------------------------
// Search
// ---------------------------------------------------------------------------

/// Find and certify an annihilating operator of the rows of `triangle`
/// (module documentation, "The search").
///
/// # Errors
///
/// As [`certify_hypergeometric_row_operator`] for the definition;
/// [`RecurrenceError::BudgetExhausted`] when no operator exists within the
/// budget.
pub fn derive_hypergeometric_row_recurrence(
    triangle: &HypergeometricTriangle,
    budget: &HypergeometricRowRecurrenceBudget,
) -> Result<HypergeometricRowRecurrenceCertificate, RecurrenceError> {
    let terms = triangle.parse()?;
    for order in 1..=budget.max_order {
        for euler in 0..=budget.max_euler_power {
            for variable in 0..=budget.max_variable_power {
                let unknowns =
                    |degree: usize| (order * (variable + 1) * (euler + 1) + 1) * (degree + 1);
                let top = (0..=budget.max_index_degree)
                    .rev()
                    .find(|d| unknowns(*d) <= budget.max_unknowns);
                let Some(top) = top else { continue };
                if solve(&terms, order, variable, euler, top)?.is_none() {
                    continue;
                }
                for degree in 0..=top {
                    if let Some(operator) = solve(&terms, order, variable, euler, degree)? {
                        let terms = operator_terms(&operator);
                        return certify_with_limit(triangle, &terms, budget.max_checked_row);
                    }
                }
            }
        }
    }
    Err(RecurrenceError::BudgetExhausted(
        "no annihilating operator with a free leader within the budget".into(),
    ))
}

/// The kernel vector with a nonzero leader of the ansatz, made primitive.
fn solve(
    terms: &[Term],
    order: usize,
    variable: usize,
    euler: usize,
    degree: usize,
) -> Result<Option<EulerMap>, RecurrenceError> {
    // Unknown columns: (i, j, l, d), the leader block last.
    let mut columns = Vec::new();
    for i in 0..order {
        for j in 0..=variable {
            for l in 0..=euler {
                for d in 0..=degree {
                    columns.push((i, j, l, d));
                }
            }
        }
    }
    let leader_start = columns.len();
    for d in 0..=degree {
        columns.push((order, 0, 0, d));
    }
    let mut pairs: BTreeSet<(usize, usize)> = BTreeSet::new();
    for i in 0..order {
        for j in 0..=variable {
            pairs.insert((i, j));
        }
    }
    pairs.insert((order, 0));
    // Equations: one row per (block, monomial).
    let mut equations: BTreeMap<(usize, usize, usize), Vec<(usize, Q)>> = BTreeMap::new();
    let mut block_id = 0;
    for term in terms {
        for residue in 0..term.modulus {
            let (blocks, _) = blocks_for(term, residue, &pairs)?;
            for block in &blocks {
                for (column, (i, j, l, d)) in columns.iter().enumerate() {
                    let Some(numerator) = block.numerators.get(&(*i, *j)) else {
                        continue;
                    };
                    let value = Biv::monomial(*d, 0, Q::one())
                        .multiply(&column_power(&block.column_value, *j, *l))
                        .multiply(numerator);
                    for (a, row) in value.rows().iter().enumerate() {
                        for (b, c) in row.iter().enumerate() {
                            if !c.is_zero() {
                                equations
                                    .entry((block_id, a, b))
                                    .or_default()
                                    .push((column, c.clone()));
                            }
                        }
                    }
                }
                block_id += 1;
            }
        }
    }
    let width = columns.len();
    let matrix: Vec<Vec<BigInt>> = equations
        .values()
        .map(|entries| {
            let (_, integers) =
                vec_poly::clear_denominators_with_scale_iter(entries.iter().map(|(_, c)| c));
            let mut row = vec![BigInt::zero(); width];
            for ((column, _), integer) in entries.iter().zip(integers) {
                row[*column] += integer;
            }
            row
        })
        .collect();
    let group: Vec<usize> = (0..width).map(|c| usize::from(c >= leader_start)).collect();
    let vectors = ore_algebra::exact_kernel_vectors_by_group(&matrix, &group, 2);
    let Some(Some(vector)) = vectors.into_iter().nth(1) else {
        return Ok(None);
    };
    let mut operator: EulerMap = BTreeMap::new();
    for ((i, j, l, d), value) in columns.iter().zip(vector) {
        if value.is_zero() {
            continue;
        }
        let entry = operator.entry((*i, *j, *l)).or_default();
        if entry.len() <= *d {
            entry.resize(d + 1, Q::zero());
        }
        entry[*d] = value;
    }
    operator.retain(|_, c| !vec_poly::is_zero(c));
    if !operator.contains_key(&(order, 0, 0)) {
        return Ok(None);
    }
    Ok(Some(primitive(operator)))
}

/// Divide by the gcd of the coefficients in `Q[n]` and scale to coprime
/// integers with a positive leading coefficient of the leader.
fn primitive(mut operator: EulerMap) -> EulerMap {
    let mut gcd = crate::polynomial::Polynomial::<Q>::zero();
    for c in operator.values() {
        gcd = gcd.gcd(&crate::polynomial::Polynomial::new(c.clone()));
    }
    if vec_poly::degree(gcd.coeffs()).unwrap_or(0) > 0 {
        for c in operator.values_mut() {
            let (quotient, _) = crate::polynomial::Polynomial::new(c.clone()).div_rem(&gcd);
            *c = vec_poly::trim(quotient.coeffs());
        }
    }
    let all: Vec<Q> = operator.values().flatten().cloned().collect();
    let (denominator, integers) = vec_poly::clear_denominators_with_scale(&all);
    let numerator = integers.iter().fold(BigInt::zero(), |acc, c| acc.gcd(c));
    let leader_sign = operator
        .iter()
        .next_back()
        .and_then(|(_, c)| c.iter().rev().find(|v| !v.is_zero()))
        .is_some_and(Signed::is_negative);
    let mut scale = Q::new(denominator, numerator);
    if leader_sign {
        scale = -scale;
    }
    for c in operator.values_mut() {
        *c = vec_poly::trim(&c.iter().map(|v| v * &scale).collect::<Vec<_>>());
    }
    operator
}
