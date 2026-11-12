//! Hypergeometric terms and their support bounds.

use super::*;

// ---------------------------------------------------------------------------
// Parsed terms
// ---------------------------------------------------------------------------

/// `ceil((a n + b) / d)` (lower) or `floor((a n + b) / d)` (upper).
#[derive(Clone, Debug)]
pub(super) struct Bound {
    a: BigInt,
    b: BigInt,
    d: BigInt,
    round_up: bool,
}

impl Bound {
    fn parse(text: &str, round_up: bool) -> Result<Self, RecurrenceError> {
        let form = parse_linear(text)?;
        if !form.column.is_zero() {
            return Err(RecurrenceError::InvalidInput(format!(
                "support bound `{text}` must not depend on k"
            )));
        }
        let (d, integers) =
            vec_poly::clear_denominators_with_scale_iter([&form.index, &form.constant].into_iter());
        let [a, b]: [BigInt; 2] = integers.try_into().expect("two coefficients");
        if a.is_negative() {
            return Err(RecurrenceError::Unsupported(format!(
                "support bound `{text}` decreases in n"
            )));
        }
        Ok(Self { a, b, d, round_up })
    }

    fn at(&self, n: i64) -> i64 {
        let value = &self.a * BigInt::from(n) + &self.b;
        let rounded = if self.round_up {
            value.div_ceil(&self.d)
        } else {
            value.div_floor(&self.d)
        };
        rounded.to_i64().unwrap_or(i64::MAX)
    }

    /// The bound on the residue class `n = r mod modulus` as an exact
    /// linear form in `n` (`d` divides `modulus`).
    pub(super) fn affine(&self, residue: usize) -> Linear {
        let r = BigInt::from(residue);
        let value = &self.a * &r + &self.b;
        let rounded = if self.round_up {
            value.div_ceil(&self.d)
        } else {
            value.div_floor(&self.d)
        };
        let offset = &rounded * &self.d - &value;
        Linear {
            index: Q::new(self.a.clone(), self.d.clone()),
            column: Q::zero(),
            constant: Q::new(&self.b + offset, self.d.clone()),
        }
    }
}

#[derive(Clone, Debug)]
pub(super) struct Term {
    pub(super) first_index: usize,
    pub(super) column: Vec<Factored>,
    pub(super) index: Vec<Factored>,
    pub(super) lower: Bound,
    pub(super) upper: Bound,
    pub(super) first_value: Q,
    /// The polynomial factor `M(n, k)` (one when absent).
    pub(super) multiplier: Biv,
    /// The period of ratios and supports.
    pub(super) modulus: usize,
}

impl Term {
    pub(super) fn parse(
        definition: &HypergeometricTermDefinition,
    ) -> Result<Self, RecurrenceError> {
        let period = definition.column_ratios.len();
        if period == 0 || definition.index_ratios.len() != period {
            return Err(RecurrenceError::InvalidInput(
                "column_ratios and index_ratios need the same positive length".into(),
            ));
        }
        let column = definition
            .column_ratios
            .iter()
            .map(|text| parse_factored(text))
            .collect::<Result<Vec<_>, _>>()?;
        let index = definition
            .index_ratios
            .iter()
            .map(|text| parse_factored(text))
            .collect::<Result<Vec<_>, _>>()?;
        let lower = Bound::parse(&definition.lower_column, true)?;
        let upper = Bound::parse(&definition.upper_column, false)?;
        let first_value = parse_rational_coeff(&definition.first_value)?;
        if first_value.is_zero() {
            return Err(RecurrenceError::InvalidInput(
                "the first value must be nonzero".into(),
            ));
        }
        let lcm = |x: usize, y: &BigInt| -> Result<usize, RecurrenceError> {
            let y = y.to_usize().ok_or_else(|| {
                RecurrenceError::LimitExceeded("support denominator too large".into())
            })?;
            Ok(x.lcm(&y))
        };
        let modulus = lcm(lcm(period, &lower.d)?, &upper.d)?;
        if modulus > 64 {
            return Err(RecurrenceError::LimitExceeded(format!(
                "period {modulus} of the ratios and supports exceeds 64"
            )));
        }
        let multiplier = match &definition.polynomial_factor {
            Some(text) => factored::parse_polynomial(text)?,
            None => Biv::one(),
        };
        Ok(Self {
            first_index: definition.first_index,
            column,
            index,
            lower,
            upper,
            first_value,
            multiplier,
            modulus,
        })
    }

    fn period(&self) -> usize {
        self.column.len()
    }

    pub(super) fn column_ratio(&self, n_residue: usize) -> &Factored {
        &self.column[n_residue % self.period()]
    }

    pub(super) fn index_ratio(&self, n_residue: usize) -> &Factored {
        &self.index[n_residue % self.period()]
    }

    /// The first `n >= first_index` with `n = residue mod modulus`.
    pub(super) fn first_in_class(&self, residue: usize) -> usize {
        let n0 = self.first_index;
        n0 + (residue + self.modulus - n0 % self.modulus) % self.modulus
    }

    /// Step 1 of the module documentation.
    pub(super) fn check_definition(&self) -> Result<(), RecurrenceError> {
        // Consistency of the two ratios.
        let k = Linear::column_variable();
        let k_plus_one = k.add(&Linear::constant(Q::one()));
        for r in 0..self.period() {
            let left = self
                .column_ratio(r + 1)
                .compose(1, &k)?
                .multiply(self.index_ratio(r));
            let right = self
                .index_ratio(r)
                .compose(0, &k_plus_one)?
                .multiply(self.column_ratio(r));
            if left != right {
                return Err(RecurrenceError::Unsupported(format!(
                    "the ratios are not consistent (residue {r}): rho_k(n+1,k) rho_n(n,k) = {} \
                     but rho_n(n,k+1) rho_k(n,k) = {}",
                    left.text(),
                    right.text()
                )));
            }
        }
        for residue in 0..self.modulus {
            let lower = self.lower.affine(residue);
            let upper = self.upper.affine(residue);
            let next_residue = (residue + 1) % self.modulus;
            let lower_next = self
                .lower
                .affine(next_residue)
                .compose(1, &Linear::constant(Q::zero()));
            let start = self.first_in_class(residue);
            let n = Q::from_integer(BigInt::from(start));
            let zero = Q::zero();
            // Nonempty, monotone supports with lower(n) >= 0 and lower(n+1) <= upper(n).
            let width = upper.add(&lower.scale(&-Q::one()));
            let step = upper.add(&lower_next.scale(&-Q::one()));
            if lower.evaluate(&n, &zero).is_negative()
                || width.evaluate(&n, &zero).is_negative()
                || width.index.is_negative()
                || step.evaluate(&n, &zero).is_negative()
                || step.index.is_negative()
            {
                return Err(RecurrenceError::Unsupported(format!(
                    "the support must satisfy 0 <= lower(n) <= upper(n) and lower(n+1) <= \
                     upper(n) for every n >= {} (residue {residue} mod {})",
                    self.first_index, self.modulus
                )));
            }
            let column_region = (lower.clone(), upper.add(&Linear::constant(-Q::one())));
            let index_region = (lower_next, upper.clone());
            for (ratio, region, name) in [
                (self.column_ratio(residue), column_region, "T(n,k+1)/T(n,k)"),
                (self.index_ratio(residue), index_region, "T(n+1,k)/T(n,k)"),
            ] {
                for phi in ratio.factors.keys() {
                    if !nonzero_on_region(phi, residue, self.modulus, start, &region) {
                        return Err(RecurrenceError::Unsupported(format!(
                            "the factor {} of {name} may vanish on the support (residue \
                             {residue} mod {}); the term must be nonzero on its support",
                            phi.text(),
                            self.modulus
                        )));
                    }
                }
            }
        }
        Ok(())
    }

    /// The rows of this term as `(lower(n), values)` for `count` rows.
    pub(super) fn rows(&self, count: usize) -> Result<Vec<(i64, Vec<Q>)>, RecurrenceError> {
        let mut out: Vec<(i64, Vec<Q>)> = Vec::with_capacity(count);
        let n0 = i64::try_from(self.first_index)
            .map_err(|_| RecurrenceError::LimitExceeded("first index too large".into()))?;
        let evaluate = |ratio: &Factored, n: i64, k: i64| -> Result<Q, RecurrenceError> {
            ratio
                .evaluate(&Q::from_integer(n.into()), &Q::from_integer(k.into()))
                .filter(|v| !v.is_zero())
                .ok_or_else(|| {
                    RecurrenceError::Unsupported(format!(
                        "a ratio vanishes or is undefined at (n, k) = ({n}, {k})"
                    ))
                })
        };
        for step in 0..count {
            let n = n0 + step as i64;
            let residue = (n as usize) % self.period();
            let lower = self.lower.at(n);
            let upper = self.upper.at(n);
            let first = if step == 0 {
                self.first_value.clone()
            } else {
                let (previous_lower, previous) = &out[step - 1];
                let position = usize::try_from(lower - previous_lower).map_err(|_| {
                    RecurrenceError::Unsupported("the lower support bound decreases".into())
                })?;
                let previous_residue = ((n - 1) as usize) % self.period();
                let value = previous.get(position).ok_or_else(|| {
                    RecurrenceError::Unsupported(format!("row {n} starts after row {} ends", n - 1))
                })?;
                value * evaluate(self.index_ratio(previous_residue), n - 1, lower)?
            };
            let mut values = vec![first];
            for k in lower..upper {
                let next =
                    values.last().expect("nonempty") * evaluate(self.column_ratio(residue), n, k)?;
                values.push(next);
            }
            out.push((lower, values));
        }
        Ok(out)
    }

    /// The rows of `T = h M` (module documentation).
    pub(super) fn rows_with_factor(
        &self,
        count: usize,
    ) -> Result<Vec<(i64, Vec<Q>)>, RecurrenceError> {
        let mut rows = self.rows(count)?;
        if self.multiplier == Biv::one() {
            return Ok(rows);
        }
        for (step, (lower, values)) in rows.iter_mut().enumerate() {
            let n = Q::from_integer(BigInt::from(self.first_index + step));
            let at_n = self.multiplier.evaluate_index(&n);
            for (offset, value) in values.iter_mut().enumerate() {
                let k = Q::from_integer(BigInt::from(*lower + offset as i64));
                *value *= vec_poly::evaluate(&at_n, &k);
            }
        }
        Ok(rows)
    }
}

/// Whether the linear form `phi` is nonzero at every lattice point
/// `(n, k)` with `n >= start`, `n = residue mod modulus` and
/// `region.0(n) <= k <= region.1(n)` (module documentation, step 1).
fn nonzero_on_region(
    phi: &Linear,
    residue: usize,
    modulus: usize,
    start: usize,
    region: &(Linear, Linear),
) -> bool {
    let (low, high) = region;
    let zero = Q::zero();
    let slope = &high.index - &low.index;
    let mut n = start;
    let gap = |n: usize| {
        let n = Q::from_integer(n.into());
        high.evaluate(&n, &zero) - low.evaluate(&n, &zero)
    };
    if gap(n).is_negative() {
        if !slope.is_positive() {
            return true; // empty region
        }
        let needed = (-gap(n) / (&slope * Q::from_integer(modulus.into()))).ceil();
        let Some(steps) = needed.to_integer().to_usize() else {
            return false;
        };
        n += steps * modulus;
        if gap(n).is_negative() {
            return false;
        }
    }
    let nq = Q::from_integer(n.into());
    let corners = [
        phi.evaluate(&nq, &low.evaluate(&nq, &zero)),
        phi.evaluate(&nq, &high.evaluate(&nq, &zero)),
    ];
    let directions = [
        &phi.index + &phi.column * &low.index,
        &phi.index + &phi.column * &high.index,
    ];
    if corners.iter().all(Signed::is_positive) && directions.iter().all(|d| !d.is_negative()) {
        return true;
    }
    if corners.iter().all(Signed::is_negative) && directions.iter().all(|d| !d.is_positive()) {
        return true;
    }
    // No integer zero at all: A modulus m + B k + (A r + C) with
    // gcd(A modulus, B) not dividing A r + C.
    let (_, integers) = vec_poly::clear_denominators_with_scale_iter(
        [&phi.index, &phi.column, &phi.constant].into_iter(),
    );
    let [index, b, constant]: [BigInt; 3] = integers.try_into().expect("three coefficients");
    let a = &index * BigInt::from(modulus);
    let c = index * BigInt::from(residue) + constant;
    let g = a.gcd(&b);
    !g.is_zero() && !c.is_multiple_of(&g)
}
