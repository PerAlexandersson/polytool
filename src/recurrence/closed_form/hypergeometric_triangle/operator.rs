//! Euler operators and hypergeometric ratio blocks.

use super::*;

// ---------------------------------------------------------------------------
// Identities
// ---------------------------------------------------------------------------

/// One identity block: the column index `k` of the coefficient
/// (`k` itself in the interior, `lower(n) + m` or `upper(n) + m` in a strip)
/// and, for every pair `(i, j)` whose point is in the support, the
/// numerator `N_(i,j)` of `T(n + i, k - j) / base` over the block's common
/// denominator.  The identity is `sum c_(i,j,l)(n) (k - j)^l N_(i,j) = 0`.
pub(super) struct Block {
    pub(super) column_value: Linear,
    pub(super) numerators: BTreeMap<(usize, usize), Biv>,
}

/// The ratio `T(n + i, k') / base` along a path, as a factored function.
struct PathBuilder<'a> {
    term: &'a Term,
    residue: usize,
}

impl PathBuilder<'_> {
    fn residue_after(&self, s: usize) -> usize {
        (self.residue + s) % self.term.modulus
    }

    /// `prod_(s < i) rho_n(n + s, column)`.
    fn up(
        &self,
        from_shift: usize,
        steps: usize,
        column: &Linear,
    ) -> Result<Factored, RecurrenceError> {
        let mut out = Factored::one();
        for s in from_shift..from_shift + steps {
            let ratio = self.term.index_ratio(self.residue_after(s));
            out = out.multiply(&ratio.compose(s as i64, column)?);
        }
        Ok(out)
    }

    /// Along row `n + shift` from `column` by `steps` (negative: backwards).
    fn along(
        &self,
        shift: usize,
        column: &Linear,
        steps: i64,
    ) -> Result<Factored, RecurrenceError> {
        let ratio = self.term.column_ratio(self.residue_after(shift));
        let mut out = Factored::one();
        if steps >= 0 {
            for u in 0..steps {
                let at = column.add(&Linear::constant(Q::from_integer(u.into())));
                out = out.multiply(&ratio.compose(shift as i64, &at)?);
            }
        } else {
            for u in 1..=-steps {
                let at = column.add(&Linear::constant(Q::from_integer((-u).into())));
                out = out.multiply(&ratio.compose(shift as i64, &at)?.inverse()?);
            }
        }
        Ok(out)
    }
}

/// The blocks of one term and residue class for the pairs `pairs`, with
/// the strip counts and the separation index (module documentation, steps
/// 2 and 3).
pub(super) fn blocks_for(
    term: &Term,
    residue: usize,
    pairs: &BTreeSet<(usize, usize)>,
) -> Result<(Vec<Block>, HypergeometricBoundaryStrips), RecurrenceError> {
    let order = pairs.iter().map(|(i, _)| *i).max().unwrap_or(0);
    let width = pairs.iter().map(|(_, j)| *j).max().unwrap_or(0);
    let path = PathBuilder { term, residue };
    let modulus = term.modulus;
    let lower = |s: usize| {
        term.lower
            .affine((residue + s) % modulus)
            .compose(s as i64, &Linear::constant(Q::zero()))
    };
    let upper = |s: usize| {
        term.upper
            .affine((residue + s) % modulus)
            .compose(s as i64, &Linear::constant(Q::zero()))
    };
    let difference = |a: &Linear, b: &Linear| -> Result<i64, RecurrenceError> {
        let d = a.add(&b.scale(&-Q::one()));
        if !d.index.is_zero() || !d.constant.is_integer() {
            return Err(RecurrenceError::VerificationFailed(
                "support offsets are not integer constants".into(),
            ));
        }
        d.constant
            .to_integer()
            .to_i64()
            .ok_or_else(|| RecurrenceError::LimitExceeded("support offset".into()))
    };
    let lower_steps = (0..=order)
        .map(|i| difference(&lower(i), &lower(0)))
        .collect::<Result<Vec<_>, _>>()?;
    let upper_steps = (0..=order)
        .map(|i| difference(&upper(i), &upper(0)))
        .collect::<Result<Vec<_>, _>>()?;
    let width_i64 = width as i64;
    let lower_strips = lower_steps.iter().copied().max().unwrap_or(0) + width_i64;
    let upper_strips = upper_steps.iter().copied().max().unwrap_or(0) + width_i64;

    let mut blocks = Vec::new();
    // Interior: base (n, k).
    let k = Linear::column_variable();
    let mut interior = BTreeMap::new();
    for &(i, j) in pairs {
        let ratio = path
            .up(0, i, &k)?
            .multiply(&path.along(i, &k, -(j as i64))?);
        let at = k.add(&Linear::constant(-Q::from_integer(j.into())));
        let factor = factored::compose_polynomial(&term.multiplier, i as i64, &at);
        interior.insert((i, j), (ratio, factor));
    }
    blocks.push(make_block(k, interior));
    // Lower strips: base (n, lower(n)).
    let base_low = lower(0);
    for m in 0..lower_strips {
        let mut ratios = BTreeMap::new();
        for &(i, j) in pairs {
            let offset = m - j as i64;
            if offset < lower_steps[i] {
                continue;
            }
            let column = base_low.add(&Linear::constant(Q::from_integer(offset.into())));
            let ratio = path
                .along(0, &base_low, offset)?
                .multiply(&path.up(0, i, &column)?);
            let factor = factored::compose_polynomial(&term.multiplier, i as i64, &column);
            ratios.insert((i, j), (ratio, factor));
        }
        let value = base_low.add(&Linear::constant(Q::from_integer(m.into())));
        blocks.push(make_block(value, ratios));
    }
    // Upper strips: base (n, upper(n)).
    let base_high = upper(0);
    for m in 1..=upper_strips {
        let mut ratios = BTreeMap::new();
        for &(i, j) in pairs {
            let offset = m - j as i64;
            if offset > upper_steps[i] {
                continue;
            }
            let ratio = path
                .up(0, i, &base_high)?
                .multiply(&path.along(i, &base_high, offset)?);
            let column = base_high.add(&Linear::constant(Q::from_integer(offset.into())));
            let factor = factored::compose_polynomial(&term.multiplier, i as i64, &column);
            ratios.insert((i, j), (ratio, factor));
        }
        let value = base_high.add(&Linear::constant(Q::from_integer(m.into())));
        blocks.push(make_block(value, ratios));
    }
    // Separation index: lower(n) + lower_strips - 1 <= upper(n),
    // upper(n) + 1 - width >= lower(n + order) and upper(n) >= lower(n + order).
    let conditions = [
        upper(0)
            .add(&lower(0).scale(&-Q::one()))
            .add(&Linear::constant(Q::from_integer(
                (1 - lower_strips).into(),
            ))),
        upper(0)
            .add(&lower(order).scale(&-Q::one()))
            .add(&Linear::constant(Q::from_integer((1 - width_i64).into()))),
        upper(0).add(&lower(order).scale(&-Q::one())),
    ];
    let start = term.first_in_class(residue);
    let mut separation = start;
    for condition in &conditions {
        separation = separation.max(first_nonnegative(condition, start, modulus)?);
    }
    let strips = HypergeometricBoundaryStrips {
        term: 0,
        residue,
        modulus,
        lower_strips: lower_strips as usize,
        upper_strips: upper_strips as usize,
        separation_index: separation,
    };
    Ok((blocks, strips))
}

/// The first `n >= start` in the residue class of `start` modulo `modulus`
/// with `condition(n) >= 0`, for a condition nondecreasing in `n`.
fn first_nonnegative(
    condition: &Linear,
    start: usize,
    modulus: usize,
) -> Result<usize, RecurrenceError> {
    let zero = Q::zero();
    let value = condition.evaluate(&Q::from_integer(start.into()), &zero);
    if !value.is_negative() {
        return Ok(start);
    }
    if !condition.index.is_positive() {
        return Err(RecurrenceError::Unsupported(
            "the rows are too narrow for the boundary strips to separate".into(),
        ));
    }
    let steps = (-value / (&condition.index * Q::from_integer(modulus.into()))).ceil();
    let steps = steps
        .to_integer()
        .to_usize()
        .ok_or_else(|| RecurrenceError::LimitExceeded("separation index".into()))?;
    Ok(start + steps * modulus)
}

fn make_block(column_value: Linear, ratios: BTreeMap<(usize, usize), (Factored, Biv)>) -> Block {
    let mut denominator: BTreeMap<Linear, i64> = BTreeMap::new();
    for (ratio, _) in ratios.values() {
        for (phi, e) in ratio.denominator_exponents() {
            let entry = denominator.entry(phi).or_insert(0);
            *entry = (*entry).max(e);
        }
    }
    let numerators = ratios
        .into_iter()
        .map(|(key, (ratio, factor))| (key, ratio.numerator_over(&denominator).multiply(&factor)))
        .collect();
    Block {
        column_value,
        numerators,
    }
}

/// `(column_value - j)^l` as a polynomial in `(n, k)`.
pub(super) fn column_power(column_value: &Linear, j: usize, l: usize) -> Biv {
    column_value
        .add(&Linear::constant(-Q::from_integer(j.into())))
        .to_biv()
        .power(l)
}

/// `sum c (k - j)^l N_(i,j)` of a block for an operator.
pub(super) fn block_value(block: &Block, operator: &EulerMap) -> Biv {
    let mut total = Biv::zero();
    for ((i, j, l), c) in operator {
        let Some(numerator) = block.numerators.get(&(*i, *j)) else {
            continue;
        };
        let term = Biv::from_index_polynomial(c)
            .multiply(&column_power(&block.column_value, *j, *l))
            .multiply(numerator);
        total = total.add(&term);
    }
    total
}

// ---------------------------------------------------------------------------
// Certification
// ---------------------------------------------------------------------------

/// Parse and normalize operator terms: `(i, j, l) -> c(n)`.
pub(super) fn parse_operator(operator: &[EulerOperatorTerm]) -> Result<EulerMap, RecurrenceError> {
    let mut out: EulerMap = BTreeMap::new();
    for term in operator {
        let c = term
            .coefficient
            .iter()
            .map(|text| parse_rational_coeff(text))
            .collect::<Result<Vec<_>, _>>()?;
        let key = (term.shift, term.variable_power, term.euler_power);
        let sum = vec_poly::add(out.get(&key).map(Vec::as_slice).unwrap_or(&[]), &c);
        out.insert(key, sum);
    }
    out.retain(|_, c| !vec_poly::is_zero(c));
    Ok(out)
}

pub(super) fn operator_terms(operator: &EulerMap) -> Vec<EulerOperatorTerm> {
    operator
        .iter()
        .map(|((i, j, l), c)| EulerOperatorTerm {
            shift: *i,
            variable_power: *j,
            euler_power: *l,
            coefficient: vec_poly::trim(c)
                .iter()
                .map(format_rational_coeff)
                .collect(),
        })
        .collect()
}

pub(super) fn operator_text(operator: &EulerMap) -> String {
    let mut parts = Vec::new();
    for ((i, j, l), c) in operator.iter().rev() {
        let coefficient = crate::recurrence::implicit::algebra::poly_text(c, "n");
        let t = match j {
            0 => String::new(),
            1 => " t".into(),
            _ => format!(" t^{j}"),
        };
        let theta = match l {
            0 => String::new(),
            1 => " theta".into(),
            _ => format!(" theta^{l}"),
        };
        let row = match i {
            0 => "P(n)".to_string(),
            _ => format!("P(n+{i})"),
        };
        parts.push(format!("({coefficient}){t}{theta} {row}"));
    }
    format!("{} = 0", parts.join(" + "))
}

/// The derivative form of `E` as an Ore operator.
pub(super) fn derivative_operator(operator: &EulerMap) -> DOp {
    let mut theta = ThetaOp::new();
    for ((i, j, l), c) in operator {
        let mut rows: Vec<Vec<Q>> = c.iter().map(|_| vec![Q::zero(); l + 1]).collect();
        for (a, value) in c.iter().enumerate() {
            rows[a][*l] = value.clone();
        }
        theta_add(&mut theta, (*i, *j), Biv::new(rows));
    }
    from_theta(&theta)
}
