//! Implicit-search sources and finite-theorem context.

use super::*;

// ---------------------------------------------------------------------------
// Context
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum Theorem {
    Identity,
    Factorial(usize),
    Derivative(usize),
    /// A Möbius-class map: `Q` satisfies a monic constant-index recurrence
    /// of the source order over `Q(x)`; `regular_at_zero` when its
    /// coefficients have no pole at `x = 0` (needed for conditions).
    Transported {
        regular_at_zero: bool,
    },
    /// A diagonal weight `scale (k!)^(-m) prod f_j(n,k)^(e_j)` with every
    /// `e_j > 0`: `total` and `index` are the total `(n,k)` degree and the
    /// `n` degree of the polynomial part.
    Weight {
        m: usize,
        total: usize,
        index: usize,
    },
}

/// `(total degree, n degree)` of the nonzero entries of `coeffs[a][b]`
/// (`n^a k^b`); `None` for a zero or unparsable polynomial.
fn weight_factor_degrees(coeffs: &BivarPolyJson) -> Option<(usize, usize)> {
    let mut degrees = None;
    for (a, row) in coeffs.coeffs.iter().enumerate() {
        for (b, c) in row.iter().enumerate() {
            if !parse_rational_coeff(c).ok()?.is_zero() {
                let (t, i) = degrees.unwrap_or((0, 0));
                degrees = Some((t.max(a + b), i.max(a)));
            }
        }
    }
    degrees
}

/// The weight theorem of a diagonal map whose weight is a polynomial in
/// `(n, k)` times `(k!)^(-m)`.
fn weight_theorem(weight: &DiagonalWeight) -> Option<Theorem> {
    match weight {
        DiagonalWeight::Binomial { m } => Some(Theorem::Weight {
            m: 0,
            total: *m,
            index: 0,
        }),
        DiagonalWeight::Product {
            factorial_exponent,
            factors,
            ..
        } => {
            if *factorial_exponent > 0 || factors.iter().any(|f| f.exponent <= 0) {
                return None;
            }
            let m = usize::try_from(factorial_exponent.unsigned_abs()).ok()?;
            let (mut total, mut index) = (0usize, 0usize);
            for f in factors {
                let e = usize::try_from(f.exponent).ok()?;
                let (t, i) = weight_factor_degrees(&f.coeffs)?;
                total = total.checked_add(t.checked_mul(e)?)?;
                index = index.checked_add(i.checked_mul(e)?)?;
            }
            Some(Theorem::Weight { m, total, index })
        }
        _ => None,
    }
}

/// Whether the coefficients `a~_i = c_i / f` of the transported source rule
/// `f Q_m = sum_i c_i Q_(m-i)` are regular at `x = 0`; `None` when the
/// transported rule is not index free (it always is for a Möbius-class map
/// of a constant-index pure-lag source).
pub(super) fn transported_regular_at_zero(
    t: &RecurrenceTransformation,
    json: &RecurrenceJson,
) -> Option<bool> {
    let (rule, _, _) = json.to_recurrence_parts().ok()?;
    let TransportedRecurrence::Scalar(transported) = t.transport(&rule).ok()? else {
        return None;
    };
    let valuation = |p: &BivarPoly| -> Option<Option<usize>> {
        if p.coeffs.iter().skip(1).flatten().any(|c| !c.is_zero()) {
            return None;
        }
        Some(
            p.coeffs
                .first()
                .and_then(|row| row.iter().position(|c| !c.is_zero())),
        )
    };
    let leader = match &transported.denominator {
        Some(f) => valuation(f)?.ok_or(()).ok()?,
        None => 0,
    };
    for term in &transported.terms {
        if term.deriv_order != 0 || term.is_alternating() {
            return None;
        }
        if valuation(&term.coeff)?.is_some_and(|v| v < leader) {
            return Some(false);
        }
    }
    Some(true)
}

/// Whether a finite verification theorem of this module covers the map `t`
/// for the source rule `json` (a constant-index monic pure-lag rule valid
/// from its first index).
pub(in crate::recurrence) fn finite_theorem_applies(
    t: &RecurrenceTransformation,
    json: &RecurrenceJson,
) -> bool {
    constant_order(json).is_some() && map_theorem(Some(t), json).is_some()
}

/// The finite theorem for the map `t` (`None`: the identity) of a
/// constant-index monic pure-lag source `json`, if one is proved.
pub(super) fn map_theorem(
    t: Option<&RecurrenceTransformation>,
    json: &RecurrenceJson,
) -> Option<Theorem> {
    use RecurrenceTransformation as T;
    let Some(t) = t else {
        return Some(Theorem::Identity);
    };
    match t {
        T::Diagonal {
            weight: DiagonalWeight::FactorialPower { power },
        } => Some(Theorem::Factorial(*power)),
        T::Diagonal { weight } => weight_theorem(weight),
        T::Derivative { order } => Some(Theorem::Derivative(*order)),
        T::MultiplyPower { .. } | T::DividePower { .. } | T::Affine { .. } => {
            Some(Theorem::Transported {
                regular_at_zero: true,
            })
        }
        T::Reversal { schedule } | T::Mobius { schedule, .. } | T::FPolynomial { schedule } => {
            if schedule.last_index.is_some() {
                return None;
            }
            Some(Theorem::Transported {
                regular_at_zero: transported_regular_at_zero(t, json)?,
            })
        }
        _ => None,
    }
}

pub(super) struct Context<'a> {
    pub(super) request: &'a ImplicitRelationSearchRequest,
    pub(super) origin: usize,
    pub(super) start: usize,
    pub(super) rule: Option<(Recurrence, Vec<Vec<Q>>)>,
    pub(super) order: Option<usize>,
    pub(super) theorem: Option<Theorem>,
    source_rows: Vec<Vec<Q>>,
    pub(super) rows: Vec<Vec<Q>>,
    pub(super) defining: Option<(ImplicitRelation, usize)>,
    pub(super) forward: Option<ImplicitRelation>,
    pub(super) symbol: &'static str,
}

impl Context<'_> {
    /// Make at least `count` rows available (from the origin) when a rule
    /// generates them; rows-only contexts keep their rows.
    pub(super) fn ensure(&mut self, count: usize) -> Result<(), RecurrenceError> {
        if count <= self.rows.len() {
            return Ok(());
        }
        let Some((rule, seeds)) = &self.rule else {
            return Ok(());
        };
        if count > self.request.limits.max_rows {
            return Err(RecurrenceError::BudgetExhausted("row budget".into()));
        }
        let generated = rule
            .generate_rows_rational(seeds, self.origin, count.max(seeds.len()))
            .map_err(|e| {
                RecurrenceError::InvalidInput(format!("source rule generation failed: {e}"))
            })?;
        for row in generated.iter().skip(self.source_rows.len()) {
            let row = vec_poly::trim(row);
            if vec_poly::degree(&row).unwrap_or(0) > self.request.limits.max_degree {
                return Err(RecurrenceError::BudgetExhausted("row degree budget".into()));
            }
            let index = self.origin + self.source_rows.len();
            let mapped = map_row(self.request.transformation.as_ref(), index, &row)?;
            self.source_rows.push(row);
            self.rows.push(mapped);
        }
        Ok(())
    }

    pub(super) fn polys(&self) -> Vec<Poly> {
        self.rows.iter().map(|r| Poly::new(r.clone())).collect()
    }
}

pub(super) fn map_row(
    t: Option<&RecurrenceTransformation>,
    index: usize,
    row: &[Q],
) -> Result<Vec<Q>, RecurrenceError> {
    let Some(t) = t else {
        return Ok(vec_poly::trim(row));
    };
    let mut state = t
        .map_row(index, row)
        .map_err(|e| RecurrenceError::from(e).with_context(t.kind_name()))?;
    if state.len() != 1 {
        return Err(RecurrenceError::InvalidInput(format!(
            "{} gives {} components; the search needs scalar rows",
            t.kind_name(),
            state.len()
        )));
    }
    Ok(vec_poly::trim(&state.remove(0)))
}

/// `T` with `Q_n = T P_n` as a shift-free differential operator, when the map
/// is one (identity, `derivative`, a diagonal map with polynomial weight).
fn forward_operator(t: Option<&RecurrenceTransformation>) -> Option<ImplicitRelation> {
    let op = match t {
        None => DOp::single(0, 0, Biv::one()),
        Some(RecurrenceTransformation::Derivative { order }) => DOp::single(0, *order, Biv::one()),
        Some(RecurrenceTransformation::Diagonal { weight }) => {
            let reciprocal = match weight {
                DiagonalWeight::Binomial { m } => DiagonalWeight::InverseBinomial { m: *m },
                DiagonalWeight::Product {
                    scale,
                    factorial_exponent: 0,
                    factors,
                } if factors.iter().all(|f| f.exponent > 0) => {
                    let scale = parse_rational_coeff(scale).ok()?;
                    if scale.is_zero() {
                        return None;
                    }
                    DiagonalWeight::Product {
                        scale: format_rational_coeff(&scale.recip()),
                        factorial_exponent: 0,
                        factors: factors
                            .iter()
                            .map(|f| DiagonalWeightFactor {
                                coeffs: f.coeffs.clone(),
                                exponent: -f.exponent,
                            })
                            .collect(),
                    }
                }
                _ => return None,
            };
            return reciprocal.inverse_operator();
        }
        Some(_) => return None,
    };
    Some(ImplicitRelation::from_op(&op))
}

/// Whether the rule, its seeds and its extra initial rows form a
/// constant-index monic pure-lag source valid from the first index.
fn constant_order(json: &RecurrenceJson) -> Option<usize> {
    let (source, rows) = ConstantLagSource::from_recurrence_json(json).ok()?;
    let k = source.coefficients.len();
    let a = source
        .coefficients
        .iter()
        .map(|c| decode_poly(c))
        .collect::<Result<Vec<_>, _>>()
        .ok()?;
    let rows = rows
        .iter()
        .map(|r| decode_poly(r))
        .collect::<Result<Vec<_>, _>>()
        .ok()?;
    for j in k..rows.len() {
        let mut value = vec![Q::zero()];
        for (i, c) in a.iter().enumerate() {
            value = vec_poly::add(&value, &vec_poly::mul(c, &rows[j - k + i]));
        }
        if vec_poly::trim(&value) != vec_poly::trim(&rows[j]) {
            return None;
        }
    }
    Some(k)
}

pub(super) fn resolve<'a>(
    request: &'a ImplicitRelationSearchRequest,
) -> Result<(Context<'a>, ImplicitSearchSource, Vec<String>), RecurrenceError> {
    let limits = &request.limits;
    let mut notes = Vec::new();
    if let Some(t) = &request.transformation {
        t.validate()
            .map_err(|e| RecurrenceError::from(e).with_context("transformation"))?;
        if t.component_count() != 1 {
            return Err(RecurrenceError::InvalidInput(format!(
                "{} gives vector rows; the search needs a scalar map",
                t.kind_name()
            )));
        }
    }
    let explicit = request
        .rows
        .as_ref()
        .map(|rows| {
            rows.iter()
                .map(|r| decode_poly(r))
                .collect::<Result<Vec<_>, _>>()
        })
        .transpose()?;
    let defining = match (&request.source_recurrence, &request.oeis_id) {
        (Some(_), Some(_)) => {
            return Err(RecurrenceError::InvalidInput(
                "give at most one of source_recurrence and oeis_id".into(),
            ))
        }
        (Some(json), None) => {
            let (rule, first, seeds) = json.to_recurrence_parts()?;
            Some(("defining_rule", rule, first, seeds))
        }
        (None, Some(id)) => {
            let definition = crate::oeis::by_id(id).ok_or_else(|| {
                RecurrenceError::InvalidInput(format!("unknown OEIS family `{id}`"))
            })?;
            let (rule, seeds, first_row) = definition.recurrence_export_parts()?;
            let first = usize::try_from(first_row).map_err(|_| {
                RecurrenceError::InvalidInput(format!(
                    "OEIS family `{id}` starts at a negative row index"
                ))
            })?;
            Some(("oeis", rule, first, seeds))
        }
        (None, None) => None,
    };
    let row_count = request
        .row_count
        .unwrap_or(DEFAULT_IMPLICIT_SEARCH_ROW_COUNT);
    if row_count > limits.max_rows {
        return Err(RecurrenceError::LimitExceeded(format!(
            "{row_count} rows exceed max_rows = {}",
            limits.max_rows
        )));
    }
    let (kind, rule, origin, json) = match defining {
        Some((kind, rule, first, seeds)) => {
            if rule.terms.iter().any(|t| t.offset == 0) {
                return Err(RecurrenceError::InvalidInput(
                    "source rule terms need positive offsets".into(),
                ));
            }
            if seeds.is_empty() || rule.max_offset() > seeds.len() {
                return Err(RecurrenceError::InvalidInput(
                    "source rule needs at least max-offset seeds".into(),
                ));
            }
            if request.row_origin.is_some_and(|o| o != first) {
                return Err(RecurrenceError::InvalidInput(format!(
                    "row_origin differs from the source first index {first}"
                )));
            }
            let json = RecurrenceJson::from_recurrence_rational(&rule, first, &seeds, None);
            (kind, Some((rule, seeds)), first, Some(json))
        }
        None => ("rows_only", None, request.row_origin.unwrap_or(0), None),
    };
    let start = request.domain_start.unwrap_or(origin);
    if start < origin {
        return Err(RecurrenceError::InvalidInput(
            "domain_start precedes the row origin".into(),
        ));
    }
    let order = json.as_ref().and_then(constant_order);
    let theorem = match (&json, order) {
        (Some(json), Some(_)) => map_theorem(request.transformation.as_ref(), json),
        _ => None,
    };
    let completeness_basis = match (theorem, order, &json) {
        (Some(Theorem::Identity), Some(k), _) => format!(
            "constant-index monic pure-lag source of order {k}: chi^(u+s+1) annihilates every \
             column x^v n^u D^s P[n+j]"
        ),
        (Some(Theorem::Factorial(m)), Some(k), _) => format!(
            "constant-index monic pure-lag source of order {k} under T(x^j) = x^j/(j!)^{m}: \
             H = U + {m}(A+J) + 1 (rescaling theorem)"
        ),
        (Some(Theorem::Derivative(m)), Some(k), _) => format!(
            "constant-index monic pure-lag source of order {k} under Q = D^{m} P: \
             chi^(u+s+{m}+1) annihilates every column x^v n^u D^s Q[n+j]"
        ),
        (Some(Theorem::Transported { regular_at_zero }), Some(k), _) => format!(
            "constant-index monic pure-lag source of order {k} under a Möbius-class map: Q \
             satisfies the transported monic rule Q[n+{k}] = sum a~_i(x) Q[n+i] with \
             coefficients in Q(x) independent of n, so chi~^(u+s+1) annihilates every column \
             x^v n^u D^s Q[n+j] (transported-source theorem){}",
            if regular_at_zero {
                "; the a~_i are regular at x = 0, so coefficient conditions are complete too"
            } else {
                "; an a~_i has a pole at x = 0, so coefficient conditions are finite fits"
            }
        ),
        (Some(Theorem::Weight { m, total, index }), Some(k), _) => format!(
            "constant-index monic pure-lag source of order {k} under a diagonal weight \
             (k!)^(-{m}) times a polynomial of total degree {total} (index degree {index}) in \
             (n, k): H = {} + {total} + 1 (weight theorem)",
            if m == 0 {
                "max(u+s)".to_string()
            } else {
                format!("U + {m}(A+J)")
            }
        ),
        (None, Some(_), _) => "no finite theorem for this map: results are finite fits unless \
                               certified by division"
            .to_string(),
        (_, None, Some(_)) => "the rule is not a constant-index monic pure-lag rule valid from \
                               its origin: results are finite fits unless certified by division"
            .to_string(),
        _ => "rows only: every relation is a finite fit; absence in a support is still \
              certified by the rows"
            .to_string(),
    };
    let defining_relation = match &rule {
        Some((r, seeds)) => match ImplicitRelation::from_recurrence(r) {
            Ok((relation, shift)) => {
                let first_base = (origin + seeds.len()).checked_sub(shift).ok_or(
                    RecurrenceError::InvalidInput("the rule starts before its leader shift".into()),
                )?;
                Some((relation, first_base))
            }
            Err(error) => {
                notes.push(format!(
                    "the defining rule has no implicit form ({error}); no division certificates"
                ));
                None
            }
        },
        None => None,
    };
    let forward = forward_operator(request.transformation.as_ref());
    let mut context = Context {
        request,
        origin,
        start,
        rule,
        order,
        theorem,
        source_rows: Vec::new(),
        rows: Vec::new(),
        defining: defining_relation,
        forward,
        symbol: if request.transformation.is_some() {
            "Q"
        } else {
            "P"
        },
    };
    match &explicit {
        Some(rows) if context.rule.is_none() => {
            if rows.is_empty() {
                return Err(RecurrenceError::InvalidInput(
                    "at least one row is required".into(),
                ));
            }
            if rows.len() > limits.max_rows {
                return Err(RecurrenceError::LimitExceeded(format!(
                    "{} rows exceed max_rows = {}",
                    rows.len(),
                    limits.max_rows
                )));
            }
            for (i, row) in rows.iter().enumerate() {
                if vec_poly::degree(row).unwrap_or(0) > limits.max_degree {
                    return Err(RecurrenceError::LimitExceeded(
                        "row degree exceeds max_degree".into(),
                    ));
                }
                let mapped = map_row(request.transformation.as_ref(), origin + i, row)?;
                context.source_rows.push(vec_poly::trim(row));
                context.rows.push(mapped);
            }
        }
        Some(rows) => {
            context.ensure(row_count.max(rows.len()))?;
            if context.source_rows[..rows.len()]
                .iter()
                .zip(rows)
                .any(|(a, b)| *a != vec_poly::trim(b))
            {
                return Err(RecurrenceError::InvalidInput(
                    "rows disagree with the rows generated by the source rule".into(),
                ));
            }
        }
        None if context.rule.is_none() => {
            return Err(RecurrenceError::InvalidInput(
                "give source_recurrence, oeis_id or rows".into(),
            ))
        }
        None => context.ensure(row_count)?,
    }
    let source = ImplicitSearchSource {
        kind: kind.to_string(),
        rule: json,
        row_origin: origin,
        constant_index_order: order,
        completeness_basis,
    };
    Ok((context, source, notes))
}
