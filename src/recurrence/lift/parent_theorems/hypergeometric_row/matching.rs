//! Hypergeometric row relations and linear-factor matching.

use super::*;

/// `a_0 + a_1 n`.
pub(super) type Affine = (Q, Q);

pub(super) fn affine_function(a: &Affine) -> IndexFunction {
    IndexFunction::polynomial(vec![a.0.clone(), a.1.clone()])
}

pub(super) fn affine_add(a: &Affine, b: &Affine) -> Affine {
    (&a.0 + &b.0, &a.1 + &b.1)
}

pub(super) fn affine_neg(a: &Affine) -> Affine {
    (-a.0.clone(), -a.1.clone())
}

pub(super) fn affine_constant(c: i64) -> Affine {
    (Q::from_integer(c.into()), Q::zero())
}

pub(super) fn affine_is_integral(a: &Affine) -> bool {
    a.0.is_integer() && a.1.is_integer()
}

pub(super) fn affine_text(a: &Affine) -> String {
    IndexFunction::polynomial(vec![a.0.clone(), a.1.clone()]).text()
}

/// `j + a` as text.
pub(super) fn factor_text(a: &Affine) -> String {
    let text = affine_text(a);
    if text == "0" {
        "j".into()
    } else if let Some(rest) = text.strip_prefix('-') {
        format!("j - {rest}")
    } else {
        format!("j + {text}")
    }
}

pub(super) fn times(a: &IndexFunction, b: &IndexFunction) -> IndexFunction {
    IndexFunction::new(
        vec_poly::mul(&a.numerator, &b.numerator),
        vec_poly::mul(&a.denominator, &b.denominator),
    )
}

fn holds(function: &IndexFunction, first: usize, condition: IndexSignCondition) -> bool {
    matches!(
        function.decide(first, condition),
        IndexSignDecision::Holds { .. }
    )
}

/// `p(n, k) = scalar(n) prod_i (k + roots[i](n))`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct LinearFactors {
    pub(super) scalar: Vec<Q>,
    pub(super) roots: Vec<Affine>,
}

impl LinearFactors {
    pub(super) fn biv(&self) -> Biv {
        let mut out = Biv::from_index_polynomial(&self.scalar);
        for root in &self.roots {
            let factor = Biv::new(vec![vec![root.0.clone(), Q::one()], vec![root.1.clone()]]);
            out = out.multiply(&factor);
        }
        out
    }

    /// `p(n, k0)` as a polynomial in `n`.
    pub(super) fn at(&self, k0: &Affine) -> Vec<Q> {
        let mut out = self.scalar.clone();
        for root in &self.roots {
            let value = affine_add(root, k0);
            out = vec_poly::mul(&out, &[value.0, value.1]);
        }
        vec_poly::trim(&out)
    }

    pub(super) fn text(&self, variable: &str) -> String {
        let mut parts = vec![format!(
            "({})",
            IndexFunction::polynomial(self.scalar.clone()).text()
        )];
        for root in &self.roots {
            parts.push(format!(
                "({})",
                factor_text(root).replacen('j', variable, 1)
            ));
        }
        parts.join(" ")
    }
}

/// The factorization of `p(n, k)` (`Biv` in `(n, k)`) into a scalar in `n`
/// and linear factors `k + u(n)`, `u` affine; `Err` naming the shape
/// otherwise.
pub(super) fn linear_factors(p: &Biv) -> Result<LinearFactors, String> {
    let Some(degree) = p.variable_degree() else {
        return Err("the relation has a zero coefficient".into());
    };
    if degree == 0 {
        return Ok(LinearFactors {
            scalar: vec_poly::trim(&p.variable_coefficient(0)),
            roots: Vec::new(),
        });
    }
    let (content, primitive) = content_and_primitive(p);
    let leading = vec_poly::trim(&primitive.variable_coefficient(degree));
    if leading.len() != 1 {
        return Err(
            "the leading coefficient in k of the relation depends on n after removing \
                    the content"
                .into(),
        );
    }
    let lead = leading[0].clone();
    let coefficients: Vec<Vec<Q>> = (0..=degree)
        .map(|i| {
            vec_poly::trim(&vec_poly::scale(
                &primitive.variable_coefficient(i),
                &lead.recip(),
            ))
        })
        .collect();
    let factors =
        irreducible_factorization_monic_bivariate_with_budget(&coefficients, FACTORIZATION_BUDGET)
            .map_err(|error| format!("the relation was not factored: {error}"))?;
    let mut roots = Vec::new();
    for factor in factors {
        if factor.coefficients.len() != 2 || factor.coefficients[0].len() > 2 {
            return Err("a factor of the relation is not linear in k and affine in n".into());
        }
        let root = (
            factor.coefficients[0]
                .first()
                .cloned()
                .unwrap_or_else(Q::zero),
            factor.coefficients[0]
                .get(1)
                .cloned()
                .unwrap_or_else(Q::zero),
        );
        for _ in 0..factor.multiplicity {
            roots.push(root.clone());
        }
    }
    roots.sort();
    Ok(LinearFactors {
        scalar: vec_poly::trim(&vec_poly::scale(&content, &lead)),
        roots,
    })
}

/// The relation `p(n, k) T_k + q(n, k) T_(k-1) = 0` of smallest total degree
/// (at most [`MAX_RELATION_DEGREE`]) satisfied by every coefficient of the
/// rows `rows[i] = P_(first + i)` with index `>= from`; `None` when there is
/// none.  Checked exactly on every row used.
pub(super) fn fit_row_relation(rows: &[Vec<Q>], first: usize, from: usize) -> Option<(Biv, Biv)> {
    let used: Vec<(usize, &Vec<Q>)> = rows
        .iter()
        .enumerate()
        .map(|(i, row)| (first + i, row))
        .filter(|(n, _)| *n >= from)
        .collect();
    let coefficient = |row: &Vec<Q>, k: i64| -> Q {
        if k < 0 {
            Q::zero()
        } else {
            row.get(k as usize).cloned().unwrap_or_else(Q::zero)
        }
    };
    for degree in 1..=MAX_RELATION_DEGREE {
        let monomials: Vec<(usize, usize)> = (0..=degree)
            .flat_map(|a| (0..=degree - a).map(move |b| (a, b)))
            .collect();
        let width = 2 * monomials.len();
        let mut equations: Vec<Vec<Q>> = Vec::new();
        for (n, row) in &used {
            let n_value = Q::from_integer(BigInt::from(*n));
            for k in 0..=row.len() as i64 {
                let k_value = Q::from_integer(BigInt::from(k));
                let (current, previous) = (coefficient(row, k), coefficient(row, k - 1));
                if current.is_zero() && previous.is_zero() {
                    continue;
                }
                let mut equation = Vec::with_capacity(width);
                for value in [&current, &previous] {
                    for (a, b) in &monomials {
                        equation.push(value * n_value.pow(*a as i32) * k_value.pow(*b as i32));
                    }
                }
                equations.push(equation);
            }
        }
        if equations.len() < width + 8 {
            return None;
        }
        let integral: Vec<Vec<BigInt>> = equations
            .iter()
            .map(|equation| vec_poly::clear_denominators(equation))
            .collect();
        let Some(vector) = exact_kernel_vectors_by_group(&integral, &vec![0; width], 1)
            .into_iter()
            .next()
            .flatten()
        else {
            continue;
        };
        let solves = equations.iter().all(|equation| {
            equation
                .iter()
                .zip(&vector)
                .fold(Q::zero(), |acc, (a, b)| acc + a * b)
                .is_zero()
        });
        if !solves {
            continue;
        }
        let build = |offset: usize| {
            let mut rows = vec![vec![Q::zero(); degree + 1]; degree + 1];
            for (index, (a, b)) in monomials.iter().enumerate() {
                rows[*a][*b] = vector[offset + index].clone();
            }
            Biv::new(rows)
        };
        let (p, q) = (build(0), build(monomials.len()));
        if p.is_zero() || q.is_zero() {
            continue;
        }
        return Some((p, q));
    }
    None
}

/// `H(n) = p(n, theta) + x q(n, theta + 1)`: `[x^k] H P = p(n, k) [x^k] P +
/// q(n, k) [x^(k-1)] P`.
pub(super) fn relation_operator(p: &Biv, q: &Biv) -> DOp {
    let mut op = ThetaOp::new();
    theta_add(&mut op, (0, 0), p.clone());
    theta_add(&mut op, (0, 1), q.shift(&Q::zero(), &Q::one()));
    from_theta(&op)
}

/// The coordinates of [`check_hypergeometric_row_theorem`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Coordinates {
    Identity,
    Reversal,
}

/// The support facts (module documentation) decided from `first`: the
/// valuation bound when the valuation is not zero, and for the reversal
/// the degree bound.
pub(super) fn support_hypotheses(
    shape: &RecurrenceShape,
    p: &LinearFactors,
    q: &LinearFactors,
    valuation: &Affine,
    degree: Option<&Affine>,
    (rows, first_index, rule_start): (&[Vec<Q>], usize, usize),
    first: usize,
) -> Vec<LiftHypothesis> {
    // Rows checked exactly before the bounds hold by induction on the rule.
    let max_lag = shape.terms.keys().map(|(lag, _)| *lag).max().unwrap_or(1);
    let earlier: Vec<(usize, &Vec<Q>)> = (first_index.max(first.saturating_sub(max_lag))
        ..rule_start.max(first))
        .filter_map(|n| rows.get(n - first_index).map(|row| (n, row)))
        .collect();
    let mut out = Vec::new();
    if *valuation != affine_constant(0) {
        out.push(lower_support_hypothesis(
            shape, p, valuation, &earlier, first,
        ));
    }
    if let Some(degree) = degree {
        out.extend(upper_support_hypotheses(
            shape, p, q, degree, &earlier, first,
        ));
    }
    out
}

/// A matching of the ratio factors with multiplier sequences (module
/// documentation): `numerators[i]` is `j + num[i]`, `denominators[m]` is
/// `j + den[m]`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Matching {
    pub(super) base_numerator: usize,
    pub(super) complements: Vec<usize>,
    pub(super) pairs: Vec<(usize, usize)>,
    pub(super) reciprocals: Vec<usize>,
}

/// The ratio `R(n, j) = C(n) prod (j + num) / prod (j + den)` in the chosen
/// coordinates.
#[derive(Clone)]
pub(super) struct Ratio {
    pub(super) scale: IndexFunction,
    pub(super) numerators: Vec<Affine>,
    pub(super) denominators: Vec<Affine>,
}

impl Ratio {
    /// `R(n, t)` for a fixed integer `t >= 1`.
    pub(super) fn at(&self, t: usize) -> IndexFunction {
        let t = affine_constant(t as i64);
        let mut out = self.scale.clone();
        for a in &self.numerators {
            out = times(&out, &affine_function(&affine_add(a, &t)));
        }
        for b in &self.denominators {
            out = out.div(&affine_function(&affine_add(b, &t)));
        }
        out
    }

    pub(super) fn text(&self) -> String {
        let product = |factors: &[Affine]| {
            if factors.is_empty() {
                "1".to_string()
            } else {
                factors
                    .iter()
                    .map(|a| format!("({})", factor_text(a)))
                    .collect::<Vec<_>>()
                    .join(" ")
            }
        };
        format!(
            "({}) {} / ({})",
            self.scale.text(),
            product(&self.numerators),
            product(&self.denominators)
        )
    }
}

/// `M` of a numerator factor `j - M - 1` with integer coefficients.
pub(super) fn complement_size(a: &Affine) -> Option<Affine> {
    let size = affine_add(&affine_neg(a), &affine_constant(-1));
    affine_is_integral(&size).then_some(size)
}

/// The first matching (in a fixed search order) whose conditions all hold
/// for `n >= first`; `upper` is `d(n) - v` for the reversal.
pub(super) fn find_matching(
    ratio: &Ratio,
    first: usize,
    upper: Option<&Affine>,
) -> Option<Matching> {
    let zero = affine_constant(0);
    let one = affine_constant(1);
    let nonnegative =
        |a: &Affine| holds(&affine_function(a), first, IndexSignCondition::Nonnegative);
    let positive = |a: &Affine| holds(&affine_function(a), first, IndexSignCondition::Positive);
    for (base_numerator, a) in ratio.numerators.iter().enumerate() {
        let Some(size) = complement_size(a) else {
            continue;
        };
        if !nonnegative(&size) {
            continue;
        }
        if let Some(upper) = upper {
            if !nonnegative(&affine_add(upper, &affine_neg(&size))) {
                continue;
            }
        }
        let Some(base_denominator) = ratio.denominators.iter().position(|b| *b == zero) else {
            continue;
        };
        let mut used_denominators = vec![false; ratio.denominators.len()];
        used_denominators[base_denominator] = true;
        let mut complements = Vec::new();
        let mut pairs = Vec::new();
        let mut ok = true;
        for (index, a) in ratio.numerators.iter().enumerate() {
            if index == base_numerator {
                continue;
            }
            if let Some(other) = complement_size(a) {
                if nonnegative(&affine_add(&other, &affine_neg(&size))) {
                    complements.push(index);
                    continue;
                }
            }
            let partner = ratio.denominators.iter().enumerate().position(|(m, b)| {
                !used_denominators[m] && {
                    let difference = affine_add(a, &affine_neg(b));
                    affine_is_integral(&difference)
                        && nonnegative(&difference)
                        && positive(&affine_add(b, &one))
                }
            });
            match partner {
                Some(m) => {
                    used_denominators[m] = true;
                    pairs.push((index, m));
                }
                None => {
                    ok = false;
                    break;
                }
            }
        }
        if !ok {
            continue;
        }
        let reciprocals: Vec<usize> = (0..ratio.denominators.len())
            .filter(|m| !used_denominators[*m])
            .collect();
        if reciprocals
            .iter()
            .any(|m| !positive(&affine_add(&ratio.denominators[*m], &one)))
        {
            continue;
        }
        let sign = if (complements.len() + 1) % 2 == 0 {
            Q::one()
        } else {
            -Q::one()
        };
        if !holds(
            &ratio.scale.scale(&sign),
            first,
            IndexSignCondition::Positive,
        ) {
            continue;
        }
        return Some(Matching {
            base_numerator,
            complements,
            pairs,
            reciprocals,
        });
    }
    None
}
