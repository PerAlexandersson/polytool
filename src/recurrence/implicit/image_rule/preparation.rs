//! Parent records and prepared image operators.

use super::*;

// ---------------------------------------------------------------------------
// Parent
// ---------------------------------------------------------------------------

/// A validated parent: its relation `L` (base form), first index, seeds and
/// the first base from which `L` holds.
pub(super) struct Parent {
    pub(super) rule: Recurrence,
    pub(super) first_index: usize,
    pub(super) seeds: Vec<Vec<Q>>,
    pub(super) relation: DOp,
    pub(super) first_base: usize,
}

impl Parent {
    pub(super) fn new(json: &RecurrenceJson) -> Result<Self, RecurrenceError> {
        let (rule, first_index, rows) = json.to_recurrence_parts()?;
        if rule.terms.is_empty() {
            return Err(RecurrenceError::Unsupported(
                "the parent recurrence has no terms".into(),
            ));
        }
        if rule.terms.iter().any(|t| t.offset == 0) {
            return Err(RecurrenceError::InvalidInput(
                "parent recurrence offsets must be positive".into(),
            ));
        }
        if rule.terms.iter().any(RecurrenceTerm::is_alternating) {
            return Err(RecurrenceError::Unsupported(
                "the parent recurrence has a (-1)^n factor".into(),
            ));
        }
        if rule.inhomogeneous.as_ref().is_some_and(|g| !g.is_zero()) {
            return Err(RecurrenceError::Unsupported(
                "the parent recurrence has an inhomogeneous term".into(),
            ));
        }
        let leader = rule.denominator.clone().unwrap_or(BivarPoly {
            coeffs: vec![vec![Q::one()]],
        });
        if leader.is_zero() {
            return Err(RecurrenceError::InvalidInput(
                "the parent recurrence denominator is zero".into(),
            ));
        }
        if !Biv::from_bivar(&leader).is_variable_free() {
            return Err(RecurrenceError::Unsupported(
                "the parent leader f(n, x) must be free of x".into(),
            ));
        }
        let order = rule.max_offset();
        if rows.len() < order {
            return Err(RecurrenceError::InvalidInput(format!(
                "the parent recurrence has order {order} but only {} initial rows",
                rows.len()
            )));
        }
        let count = rule.generation_initial_count(first_index, rows.len());
        if count > MAX_IMAGE_RECURRENCE_ROWS {
            return Err(RecurrenceError::LimitExceeded(format!(
                "the parent needs {count} seeds, more than {MAX_IMAGE_RECURRENCE_ROWS}"
            )));
        }
        let rows: Vec<Vec<Q>> = rows.iter().map(|row| vec_poly::trim(row)).collect();
        let seeds = rows[..count].to_vec();
        let generated = rule
            .generate_rows_rational(&seeds, first_index, rows.len())
            .map_err(|error| {
                RecurrenceError::InvalidInput(format!(
                    "the parent recurrence cannot be replayed: {error}"
                ))
            })?;
        for (offset, (supplied, generated)) in rows.iter().zip(&generated).enumerate() {
            if supplied != &vec_poly::trim(generated) {
                return Err(RecurrenceError::InvalidInput(format!(
                    "supplied parent row {} disagrees with its recurrence",
                    first_index + offset
                )));
            }
        }
        if !bivariate_leader_nonvanishing(&leader, first_index + count) {
            return Err(RecurrenceError::Unsupported(format!(
                "the parent leader is not proved nonzero at every target index >= {}",
                first_index + count
            )));
        }
        let (relation, _) = ImplicitRelation::from_recurrence(&rule)?;
        Ok(Self {
            first_base: first_index + count - order,
            rule,
            first_index,
            seeds,
            relation: relation.op(),
        })
    }

    /// Rows `P_(start), ..., P_(start + count - 1)`, `start >= p0`.
    pub(super) fn rows(&self, start: usize, count: usize) -> Result<Vec<Vec<Q>>, RecurrenceError> {
        let offset = start - self.first_index;
        let total = (offset + count).max(self.seeds.len());
        let rows = self
            .rule
            .generate_rows_rational(&self.seeds, self.first_index, total)
            .map_err(|error| {
                RecurrenceError::InvalidInput(format!("the parent cannot be replayed: {error}"))
            })?;
        Ok(rows[offset..offset + count]
            .iter()
            .map(|row| vec_poly::trim(row))
            .collect())
    }
}

// ---------------------------------------------------------------------------
// Specification
// ---------------------------------------------------------------------------

fn parse_value(text: &str, what: &str) -> Result<Q, RecurrenceError> {
    parse_rational_coeff(text).map_err(|error| RecurrenceError::from(error).with_context(what))
}

fn parse_polynomial(values: &[String], what: &str) -> Result<Vec<Q>, RecurrenceError> {
    if values.len() > 64 {
        return Err(RecurrenceError::LimitExceeded(format!(
            "the {what} has more than 64 coefficients"
        )));
    }
    let parsed = values
        .iter()
        .map(|value| parse_value(value, what))
        .collect::<Result<Vec<_>, _>>()?;
    let trimmed = vec_poly::trim(&parsed);
    if vec_poly::is_zero(&trimmed) {
        return Err(RecurrenceError::InvalidInput(format!(
            "the {what} is the zero polynomial"
        )));
    }
    Ok(trimmed)
}

/// `p(n) != 0` for every integer `n >= from`, decided by the sign checker.
fn require_nonzero_from(p: &[Q], from: usize, what: &str) -> Result<(), RecurrenceError> {
    match decide_index_polynomial_sign(
        p,
        &IntegerRange::starting_at(from),
        IndexSignCondition::Nonzero,
    ) {
        IndexSignDecision::Holds { .. } => Ok(()),
        IndexSignDecision::FailsAt { index, .. } => Err(RecurrenceError::InvalidInput(format!(
            "the {what} vanishes at n = {index}; it must be nonzero for every n >= {from}"
        ))),
        IndexSignDecision::Undecided { reason } => Err(RecurrenceError::BudgetExhausted(format!(
            "the {what} could not be decided nonzero for every n >= {from}: {reason}"
        ))),
    }
}

/// `p(k) != 0` for every integer `k >= from` (`from` may be negative),
/// decided by exact integer roots.
pub(super) fn require_no_integer_roots_from(
    p: &[Q],
    from: i64,
    what: &str,
) -> Result<(), RecurrenceError> {
    let shifted = vec_poly::compose(p, &[Q::from_integer(from.into()), Q::one()]);
    if let Some(root) = roots_in_domain(&shifted, 0, None).first() {
        return Err(RecurrenceError::Unsupported(format!(
            "the {what} vanishes at k = {}; it must be nonzero for every integer k >= {from}",
            *root as i64 + from
        )));
    }
    Ok(())
}

/// The specification as crate steps with its gauge values.
pub(super) struct PreparedImage {
    pub(super) shift: i64,
    /// The map steps (without the gauge).
    map_steps: Vec<RowMapStep>,
    /// `mu(n)/mu(n-1) = alpha(n)/beta(n)`, `mu = lambda / delta`.
    alpha: Vec<Q>,
    beta: Vec<Q>,
    /// `mu(n0)`.
    initial_gauge: Q,
    /// The coefficient multiplier `(A, B)`, if any.
    pub(super) multiplier: Option<(Vec<Q>, Vec<Q>)>,
}

impl PreparedImage {
    pub(super) fn new(
        specification: &ImageOperatorSpecification,
        parent: &Parent,
        n0: usize,
    ) -> Result<Self, RecurrenceError> {
        let shift = specification.index_shift;
        if shift.unsigned_abs() > MAX_IMAGE_RECURRENCE_INDEX_SHIFT as u64 {
            return Err(RecurrenceError::LimitExceeded(format!(
                "|index_shift| exceeds {MAX_IMAGE_RECURRENCE_INDEX_SHIFT}"
            )));
        }
        let start = i64::try_from(n0)
            .ok()
            .and_then(|n0| n0.checked_add(shift))
            .and_then(|start| usize::try_from(start).ok());
        if start.is_none_or(|start| start < parent.first_index) {
            return Err(RecurrenceError::InvalidInput(format!(
                "the image uses parent row n0 + s = {n0} + {shift}, before the parent's first \
                 index {}",
                parent.first_index
            )));
        }
        let (map_steps, delta, multiplier) = map_steps(&specification.map, n0)?;
        require_nonzero_from(&delta, n0, "Weyl-operator denominator")?;
        let at = |p: &[Q], n: usize| vec_poly::evaluate(p, &Q::from_integer(n.into()));
        // lambda in ratio form.
        let (a, b, lambda0) = match &specification.gauge {
            ImageGauge::None => (vec![Q::one()], vec![Q::one()], Q::one()),
            ImageGauge::Ratio {
                numerator,
                denominator,
                initial_value,
            } => {
                let a = parse_polynomial(numerator, "gauge numerator")?;
                let b = parse_polynomial(denominator, "gauge denominator")?;
                let initial = parse_value(initial_value, "gauge initial value")?;
                if initial.is_zero() {
                    return Err(RecurrenceError::InvalidInput(
                        "the gauge initial value lambda(n0) must be nonzero".into(),
                    ));
                }
                require_nonzero_from(&a, n0 + 1, "gauge ratio numerator")?;
                require_nonzero_from(&b, n0 + 1, "gauge ratio denominator")?;
                (a, b, initial)
            }
            ImageGauge::Closed {
                numerator,
                denominator,
            } => {
                let p = parse_polynomial(numerator, "gauge numerator")?;
                let q = parse_polynomial(denominator, "gauge denominator")?;
                require_nonzero_from(&p, n0, "gauge numerator")?;
                require_nonzero_from(&q, n0, "gauge denominator")?;
                let previous = |f: &[Q]| vec_poly::compose(f, &[-Q::one(), Q::one()]);
                (
                    vec_poly::mul(&p, &previous(&q)),
                    vec_poly::mul(&q, &previous(&p)),
                    at(&p, n0) / at(&q, n0),
                )
            }
        };
        // mu = lambda / delta: alpha = a(n) delta(n - 1), beta = b(n) delta(n).
        let previous_delta = vec_poly::compose(&delta, &[-Q::one(), Q::one()]);
        Ok(Self {
            shift,
            map_steps,
            alpha: vec_poly::trim(&vec_poly::mul(&a, &previous_delta)),
            beta: vec_poly::trim(&vec_poly::mul(&b, &delta)),
            initial_gauge: lambda0 / at(&delta, n0),
            multiplier,
        })
    }

    /// Every step, the gauge last (omitted when `mu` is constant).
    pub(super) fn steps(&self) -> Vec<RowMapStep> {
        let mut steps = self.map_steps.clone();
        if self.alpha != self.beta {
            steps.push(RowMapStep::Gauge {
                numerator: self.alpha.clone(),
                denominator: self.beta.clone(),
            });
        }
        steps
    }

    /// The rows `Q_n0, ..., Q_(n0 + count - 1)`.
    pub(super) fn image_rows(
        &self,
        parent: &Parent,
        n0: usize,
        count: usize,
    ) -> Result<Vec<Vec<Q>>, RecurrenceError> {
        let start = (n0 as i64 + self.shift) as usize;
        let parent_rows = parent.rows(start, count)?;
        let mut mu = self.initial_gauge.clone();
        let mut out = Vec::with_capacity(count);
        for (offset, row) in parent_rows.iter().enumerate() {
            let n = n0 + offset;
            if offset > 0 {
                let index = Q::from_integer(n.into());
                mu = mu * vec_poly::evaluate(&self.alpha, &index)
                    / vec_poly::evaluate(&self.beta, &index);
            }
            let mut current = row.clone();
            for step in &self.map_steps {
                current = step.apply_row(&current, n).ok_or_else(|| {
                    RecurrenceError::Unsupported(format!(
                        "the row map is undefined at n = {n} (for a reversal: the parent row \
                         has degree above d_n)"
                    ))
                })?;
            }
            out.push(vec_poly::trim(&vec_poly::scale(&current, &mu)));
        }
        Ok(out)
    }

    /// The shifted parent rows `P'_n0, ...`.
    pub(super) fn parent_rows(
        &self,
        parent: &Parent,
        n0: usize,
        count: usize,
    ) -> Result<Vec<Vec<Q>>, RecurrenceError> {
        parent.rows((n0 as i64 + self.shift) as usize, count)
    }
}

/// The map steps, the Weyl denominator `delta` and the multiplier `(A, B)`.
#[allow(clippy::type_complexity)]
fn map_steps(
    map: &ImageRowMap,
    n0: usize,
) -> Result<(Vec<RowMapStep>, Vec<Q>, Option<(Vec<Q>, Vec<Q>)>), RecurrenceError> {
    let reversal_degree = |schedule: &DegreeSchedule| -> Result<Vec<Q>, RecurrenceError> {
        schedule.validate()?;
        if schedule.last_index.is_some() || schedule.step < 0 || schedule.first_index > n0 {
            return Err(RecurrenceError::InvalidInput(
                "the reversal schedule must be unbounded, with a nonnegative step, and start at \
                 or before n0"
                    .into(),
            ));
        }
        let step = Q::from_integer(schedule.step.into());
        let first = Q::from_integer(schedule.first_index.into());
        let degree = Q::from_integer(schedule.first_degree.into());
        Ok(vec_poly::trim(&[degree - &step * first, step]))
    };
    let one = vec![Q::one()];
    Ok(match map {
        ImageRowMap::Identity => (Vec::new(), one, None),
        ImageRowMap::WeylOperator { terms, denominator } => {
            if terms.is_empty() || terms.len() > 64 {
                return Err(RecurrenceError::InvalidInput(
                    "the Weyl operator needs between 1 and 64 terms".into(),
                ));
            }
            let mut op = DOp::default();
            for term in terms {
                if term.deriv_order > 16 {
                    return Err(RecurrenceError::LimitExceeded(
                        "a Weyl-operator derivative order exceeds 16".into(),
                    ));
                }
                op.add_term((0, term.deriv_order), biv_from_json(&term.coeff)?);
            }
            if op.is_zero() {
                return Err(RecurrenceError::InvalidInput(
                    "the Weyl operator is zero".into(),
                ));
            }
            let delta = parse_polynomial(denominator, "Weyl-operator denominator")?;
            (vec![RowMapStep::WeylOperator(op)], delta, None)
        }
        ImageRowMap::Affine { slope, intercept } => {
            let slope = parse_value(slope, "affine slope")?;
            if slope.is_zero() {
                return Err(RecurrenceError::InvalidInput(
                    "the affine slope must be nonzero".into(),
                ));
            }
            let intercept = parse_value(intercept, "affine intercept")?;
            (vec![RowMapStep::Affine { slope, intercept }], one, None)
        }
        ImageRowMap::Reversal { schedule } => (
            vec![RowMapStep::Reversal {
                degree: reversal_degree(schedule)?,
            }],
            one,
            None,
        ),
        ImageRowMap::Mobius {
            a,
            b,
            c,
            e,
            schedule,
        } => {
            let a = parse_value(a, "Mobius a")?;
            let b = parse_value(b, "Mobius b")?;
            let c = parse_value(c, "Mobius c")?;
            let e = parse_value(e, "Mobius e")?;
            if c.is_zero() {
                return Err(RecurrenceError::InvalidInput(
                    "a Mobius map needs c != 0 (use the affine map otherwise)".into(),
                ));
            }
            let determinant = &a * &e - &b * &c;
            if determinant.is_zero() {
                return Err(RecurrenceError::InvalidInput(
                    "a Mobius map needs a e - b c != 0".into(),
                ));
            }
            // (c x + e)^d P((a x + b)/(c x + e)) = A2(R(A1(P))) with
            // A1: F(y) -> F(a/c + kappa y), kappa = (b c - a e)/c,
            // R: G(u) -> u^d G(1/u), A2: K(u) -> K(c x + e).
            let kappa = -determinant / &c;
            (
                vec![
                    RowMapStep::Affine {
                        slope: kappa,
                        intercept: a / &c,
                    },
                    RowMapStep::Reversal {
                        degree: reversal_degree(schedule)?,
                    },
                    RowMapStep::Affine {
                        slope: c,
                        intercept: e,
                    },
                ],
                one,
                None,
            )
        }
        ImageRowMap::CoefficientMultiplier {
            numerator,
            denominator,
        } => {
            let a = parse_polynomial(numerator, "multiplier numerator")?;
            let b = parse_polynomial(denominator, "multiplier denominator")?;
            require_no_integer_roots_from(&a, 0, "multiplier numerator A(k)")?;
            require_no_integer_roots_from(&b, 0, "multiplier denominator B(k)")?;
            (
                vec![RowMapStep::CoefficientMultiplier {
                    numerator: a.clone(),
                    denominator: b.clone(),
                }],
                one,
                Some((a, b)),
            )
        }
    })
}
