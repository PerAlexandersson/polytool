//! Exact coefficient-relation closure, including support and singularity checks.

use super::*;
use crate::recurrence::implicit::algebra::{from_theta, theta_add, ThetaOp};
use crate::recurrence::implicit::biv_json;
use crate::recurrence::lift::hadamard_factor::annihilator::vanishing_bases;
use crate::recurrence::lift::hadamard_factor::annihilator_system::find_invariant_annihilator_system;
use crate::recurrence::lift::hadamard_factor::certify::{check_system_rows, system_json};
use ore_algebra::exact_kernel_vectors_by_group;

/// An all-index coefficient relation for one defining rule, with its anchor
/// and degree bound. Coefficient polynomials use `(n,k)` in place of `(n,x)`.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PolynomialFamilyCoefficientRelation {
    /// `p(n,k)` in `p(n,k) a(n,k) + q(n,k) a(n,k-1) = 0`.
    pub current: BivarPolyJson,
    /// `q(n,k)` in the same relation.
    pub previous: BivarPolyJson,
    /// The invariant annihilator identities proving the relation.
    pub system: HadamardAnnihilatorSystem,
    /// All row checks before the invariant-system start and at zeros of its
    /// multiplier, including the identifying prefix.
    pub checked_indices: Vec<usize>,
    /// Exact inductive upper bound `deg P_n <= n + degree_offset`.
    pub degree_offset: i64,
    /// Exact constant coefficient, fixed for every index by the rule.
    pub constant_coefficient: String,
    /// Exact positive or negative coefficient-cone test proving
    /// `p(n,k) != 0` for `1 <= k <= n + degree_offset`.
    pub leader_decision: String,
}

fn relation_operator(p: &Biv, q: &Biv) -> DOp {
    let mut theta = ThetaOp::new();
    theta_add(&mut theta, (0, 0), p.clone());
    theta_add(&mut theta, (0, 1), q.shift(&Q::zero(), &Q::one()));
    from_theta(&theta)
}

/// A proposed relation, from exact equations only. The invariance certificate
/// below, rather than this fit, proves it for all indices.
fn fit_relation(parent: &PreparedLiftRule, count: usize) -> Result<(Biv, Biv), RecurrenceError> {
    let rows = parent.rows(count)?;
    for degree in 1..=5 {
        let monomials: Vec<_> = (0..=degree)
            .flat_map(|a| (0..=degree - a).map(move |b| (a, b)))
            .collect();
        let width = monomials.len() * 2;
        let mut equations = Vec::new();
        for (i, row) in rows.iter().enumerate() {
            let n = Q::from_integer((parent.first_index + i).into());
            for k in 0..=row.len() {
                let current = row.get(k).cloned().unwrap_or_else(Q::zero);
                let previous = k
                    .checked_sub(1)
                    .and_then(|j| row.get(j))
                    .cloned()
                    .unwrap_or_else(Q::zero);
                if current.is_zero() && previous.is_zero() {
                    continue;
                }
                let index = Q::from_integer(k.into());
                let mut equation = Vec::with_capacity(width);
                for value in [&current, &previous] {
                    equation.extend(
                        monomials
                            .iter()
                            .map(|&(a, b)| value * n.pow(a as i32) * index.pow(b as i32)),
                    );
                }
                equations.push(vec_poly::clear_denominators(&equation));
            }
        }
        if equations.len() < width + 8 {
            continue;
        }
        let groups: Vec<_> = (0..width).collect();
        for vector in exact_kernel_vectors_by_group(&equations, &groups, width)
            .into_iter()
            .flatten()
        {
            let build = |offset: usize| {
                let mut result = vec![vec![Q::zero(); degree + 1]; degree + 1];
                for (i, &(a, b)) in monomials.iter().enumerate() {
                    result[a][b] = vector[offset + i].clone();
                }
                Biv::new(result)
            };
            let (p, q) = (build(0), build(monomials.len()));
            if !p.is_zero() && !q.is_zero() && leader_decision(&p, degree_offset(parent)?).is_ok() {
                return Ok((p, q));
            }
        }
    }
    Err(unsupported(
        "no coefficient relation within total degree five",
    ))
}

/// The rule proves the envelope by induction: x^e D^d P_(n-r) has degree
/// at most n+c+e-d-r. Initial rows establish the induction window.
fn degree_offset(parent: &PreparedLiftRule) -> Result<i64, RecurrenceError> {
    if parent.rule.terms.iter().any(|term| {
        term.coeff.coeffs.iter().any(|row| {
            row.iter().enumerate().any(|(e, value)| {
                !value.is_zero() && e as i64 - term.deriv_order as i64 - term.offset as i64 > 0
            })
        })
    }) {
        return Err(unsupported(
            "the rule does not preserve a unit-slope degree envelope",
        ));
    }
    parent
        .seeds
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let row = vec_poly::trim(row);
            row.len() as i64 - 1 - (parent.first_index + i) as i64
        })
        .max()
        .ok_or_else(|| unsupported("a coefficient relation needs identifying rows"))
}

fn power(base: &Biv, exponent: usize) -> Biv {
    (0..exponent).fold(Biv::one(), |out, _| out.multiply(base))
}

/// Substitute k=1+u, n=1-c+u+v. On the whole quadrant u,v>=0 a
/// polynomial with coefficients of one sign and nonzero constant cannot
/// vanish. Every 1<=k<=n+c is in that quadrant. No numerical root test.
fn leader_decision(p: &Biv, offset: i64) -> Result<String, RecurrenceError> {
    let n = Biv::new(vec![
        vec![Q::from_integer((1 - offset).into()), Q::one()],
        vec![Q::one()],
    ]);
    let k = Biv::new(vec![vec![Q::one()], vec![Q::one()]]);
    let mut transformed = Biv::zero();
    for (a, row) in p.rows().iter().enumerate() {
        for (b, value) in row.iter().enumerate() {
            if !value.is_zero() {
                transformed = transformed.add(&power(&n, a).multiply(&power(&k, b)).scale(value));
            }
        }
    }
    let constant = transformed
        .rows()
        .first()
        .and_then(|r| r.first())
        .cloned()
        .unwrap_or_else(Q::zero);
    let values: Vec<_> = transformed.rows().iter().flatten().collect();
    if (constant.is_positive() && values.iter().all(|v| !v.is_negative()))
        || (constant.is_negative() && values.iter().all(|v| !v.is_positive()))
    {
        Ok(format!(
            "p(1-({offset})+u+v,1+u) has coefficients of one sign and a nonzero constant; u,v >= 0"
        ))
    } else {
        Err(unsupported(
            "the coefficient leader's all-index nonvanishing is undecided",
        ))
    }
}

fn certify_relation(
    parent: &PreparedLiftRule,
    p: &Biv,
    q: &Biv,
    offset: i64,
) -> Result<PolynomialFamilyCoefficientRelation, RecurrenceError> {
    let anchor = constant_coefficient(parent).ok_or_else(|| {
        unsupported("the row constant coefficient is not proved fixed by the rule")
    })?;
    let leader_decision = leader_decision(p, offset)?;
    let operator = relation_operator(p, q);
    let from = parent.valid_from_base().max(parent.first_index);
    // Check a bounded prefix first to reject incorrect proposed relations cheaply.
    let prefix = parent.rows((from - parent.first_index + 1).max(12))?;
    if prefix.iter().enumerate().any(|(i, row)| {
        !operator
            .apply_row(row, parent.first_index + i)
            .iter()
            .all(Zero::is_zero)
    }) {
        return Err(unsupported(
            "the proposed coefficient relation fails on a row",
        ));
    }
    let decomposition =
        find_invariant_annihilator_system(&operator, parent, from, MAX_HADAMARD_OPERATOR_TERMS)?
            .ok_or_else(|| {
                unsupported("no invariant annihilator system certifies the coefficient relation")
            })?;
    if decomposition.multiplier.is_zero() {
        return Err(unsupported("zero annihilator multiplier"));
    }
    let restarts = decomposition.system.restart_bases(from);
    check_system_rows(&decomposition.system.generators, parent, from, &restarts)?;
    let mut checked_indices: Vec<_> = (parent.first_index..from).collect();
    checked_indices.extend(vanishing_bases(&decomposition.multiplier, from));
    checked_indices.sort_unstable();
    checked_indices.dedup();
    let last = checked_indices
        .last()
        .copied()
        .unwrap_or(parent.first_index);
    if last - parent.first_index > 64 {
        return Err(unsupported(
            "coefficient relation exceptional row exceeds the check bound",
        ));
    }
    let rows = parent.rows(last - parent.first_index + 1)?;
    if checked_indices.iter().any(|&n| {
        !operator
            .apply_row(&rows[n - parent.first_index], n)
            .iter()
            .all(Zero::is_zero)
    }) {
        return Err(unsupported(
            "the coefficient relation fails at an exceptional index",
        ));
    }
    let mut system = system_json(
        &decomposition.system,
        &decomposition.quotients,
        from,
        restarts,
    );
    system.row_multiplier = biv_json(&decomposition.multiplier);
    Ok(PolynomialFamilyCoefficientRelation {
        current: biv_json(p),
        previous: biv_json(q),
        system,
        checked_indices,
        degree_offset: offset,
        constant_coefficient: format_rational_coeff(&anchor),
        leader_decision,
    })
}

pub(super) fn certify(
    image: &PreparedLiftRule,
    left: &PreparedLiftRule,
    right: &PreparedLiftRule,
    operation: PolynomialFamilyProductOperation,
    options: &PolynomialFamilyProductOptions,
) -> Result<PolynomialFamilyProductClosure, RecurrenceError> {
    let (pa, qa) = fit_relation(left, options.search_rows)?;
    let (pb, qb) = fit_relation(right, options.search_rows)?;
    let left_offset = degree_offset(left)?;
    let right_offset = degree_offset(right)?;
    let offset = left_offset.min(right_offset);
    if degree_offset(image)? > offset {
        return Err(unsupported("the image degree envelope is too large"));
    }
    let left_relation = certify_relation(left, &pa, &qa, left_offset)?;
    let right_relation = certify_relation(right, &pb, &qb, right_offset)?;
    if let Ok(closure) = multiplier_closure(
        image,
        right,
        left,
        &right_relation,
        &pb,
        &qb,
        true,
        operation,
    ) {
        return Ok(closure);
    }
    if let Ok(closure) = multiplier_closure(
        image,
        left,
        right,
        &left_relation,
        &pa,
        &qa,
        false,
        operation,
    ) {
        return Ok(closure);
    }
    let mut p = pa.multiply(&pb);
    let mut q = qa.multiply(&qb).scale(&(-Q::one()));
    let k = Biv::from_variable_polynomial(&[Q::zero(), Q::one()]);
    match operation {
        PolynomialFamilyProductOperation::FactorialHadamard => q = q.multiply(&k),
        PolynomialFamilyProductOperation::InverseFactorialHadamard => p = p.multiply(&k),
        _ => (),
    }
    if constant_coefficient(image)
        != constant_coefficient(left)
            .zip(constant_coefficient(right))
            .map(|(a, b)| a * b)
    {
        return Err(unsupported(
            "image constant coefficient is not the product of the parent anchors",
        ));
    }
    let image_relation = certify_relation(image, &p, &q, offset)?;
    Ok(PolynomialFamilyProductClosure::CoefficientRelations {
        left: Box::new(left_relation),
        right: Box::new(right_relation),
        image: Box::new(image_relation),
    })
}

pub(super) fn binomial_closure(
    image: &PreparedLiftRule,
    left: &PreparedLiftRule,
    right: &PreparedLiftRule,
    operation: PolynomialFamilyProductOperation,
) -> Result<PolynomialFamilyProductClosure, RecurrenceError> {
    use crate::recurrence::lift::hadamard_factor::search::certify_candidate;
    let options = HadamardFactorLiftOptions {
        checked_rows: 16,
        ..Default::default()
    };
    for right_is_binomial in [true, false] {
        let (binomial, other) = if right_is_binomial {
            (right, left)
        } else {
            (left, right)
        };
        // The bounded offset is discovered from the degree; equality, not that
        // sample, proves the second family's entire definition is binomial.
        let Ok(offset) = degree_offset(binomial) else {
            continue;
        };
        if !(-4..=4).contains(&offset) {
            continue;
        }
        let rule = Recurrence {
            terms: vec![RecurrenceTerm::new(
                1,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![Q::one(), Q::one()]]),
            )],
            denominator: None,
            inhomogeneous: None,
        };
        let first_row = binomial.rows(1)?;
        let canonical =
            RecurrenceJson::from_recurrence_rational(&rule, binomial.first_index, &first_row, None);
        let equality = certify_recurrence_equality(&record(binomial), &canonical)?;
        if !equality.is_certified() {
            continue;
        }
        // The seed has to be exactly (1+x)^(first+offset), not an arbitrary
        // real-rooted factor carried along by the same first-order rule.
        let degree = binomial.first_index as i64 + offset;
        if degree < 0 {
            continue;
        }
        let expected: Vec<Q> = (0..=degree as usize)
            .map(|k| Q::from_integer(crate::numeric::binomial_bigint(degree as usize, k)))
            .collect();
        if first_row[0] != expected {
            continue;
        }
        let mut terms = vec![HadamardFactorTerm::Binomial {
            index_offset: offset,
            coefficient_offset: 0,
        }];
        match operation {
            PolynomialFamilyProductOperation::FactorialHadamard => {
                terms.push(HadamardFactorTerm::Factorial)
            }
            PolynomialFamilyProductOperation::InverseFactorialHadamard => {
                terms.push(HadamardFactorTerm::InverseFactorial)
            }
            _ => (),
        }
        let (status, _) = certify_candidate(
            image,
            other,
            &HadamardFactor { terms },
            &VeroneseSectionGauge::None,
            &options,
        )?;
        if let HadamardFactorIdentityStatus::Certified { certificate } = status {
            return Ok(PolynomialFamilyProductClosure::BinomialFactor {
                right_is_binomial,
                binomial_equality: Box::new(equality),
                certificate,
            });
        }
    }
    Err(unsupported(
        "neither source admits certified binomial conjugation",
    ))
}

#[allow(clippy::too_many_arguments)]
fn multiplier_closure(
    image: &PreparedLiftRule,
    source: &PreparedLiftRule,
    other: &PreparedLiftRule,
    relation: &PolynomialFamilyCoefficientRelation,
    p: &Biv,
    q: &Biv,
    right_is_multiplier: bool,
    operation: PolynomialFamilyProductOperation,
) -> Result<PolynomialFamilyProductClosure, RecurrenceError> {
    use crate::recurrence::lift::hadamard_factor::search::certify_candidate;
    let anchor =
        constant_coefficient(source).ok_or_else(|| unsupported("no fixed multiplier anchor"))?;
    let mut singles = Vec::new();
    for a in -2..=4 {
        for b in 0..=2 {
            if source.first_index as i64 + a >= b {
                singles.push((a, b));
            }
        }
    }
    let mut menus: Vec<Vec<(i64, i64)>> = singles.iter().map(|&entry| vec![entry]).collect();
    if operation == PolynomialFamilyProductOperation::Hadamard {
        for i in 0..singles.len() {
            for j in i..singles.len() {
                menus.push(vec![singles[i], singles[j]]);
            }
        }
    }
    for menu in menus {
        if menu.iter().map(|(a, b)| a - b).min() != Some(relation.degree_offset) {
            continue;
        }
        let mut numerator = Biv::one();
        let mut denominator = Biv::one();
        let mut initial = vec![Q::one()];
        for &(a, b) in &menu {
            numerator = numerator.multiply(&Biv::new(vec![
                vec![Q::from_integer((a - b + 1).into()), -Q::one()],
                vec![Q::one()],
            ]));
            denominator = denominator.multiply(&Biv::from_variable_polynomial(&[
                Q::from_integer(b.into()),
                Q::one(),
            ]));
            for j in 0..b {
                initial = vec_poly::mul(&initial, &[Q::from_integer((a - j).into()), Q::one()]);
                initial = vec_poly::scale(&initial, &Q::from_integer((j + 1).into()).recip());
            }
        }
        if !p
            .multiply(&numerator)
            .add(&q.multiply(&denominator))
            .is_zero()
        {
            continue;
        }
        // On 1<=k<=n+c, every binomial coefficient and both ratio
        // denominators are strictly positive. Outside that support the
        // binomial product is zero. The certified row relation and fixed
        // anchor therefore identify the source for every index.
        let initial_biv = Biv::from_index_polynomial(&initial);
        let at_first = vec_poly::evaluate(&initial, &Q::from_integer(source.first_index.into()));
        if at_first.is_zero() {
            continue;
        }
        let gauge = VeroneseSectionGauge::Ratio {
            numerator: initial_biv
                .shift_index(-1)
                .variable_coefficient(0)
                .iter()
                .map(format_rational_coeff)
                .collect(),
            denominator: initial.iter().map(format_rational_coeff).collect(),
            initial_value: format_rational_coeff(&(&anchor / at_first)),
        };
        let mut terms: Vec<_> = menu
            .iter()
            .map(|&(a, b)| HadamardFactorTerm::Binomial {
                index_offset: a,
                coefficient_offset: b,
            })
            .collect();
        match operation {
            PolynomialFamilyProductOperation::FactorialHadamard => {
                terms.push(HadamardFactorTerm::Factorial)
            }
            PolynomialFamilyProductOperation::InverseFactorialHadamard => {
                terms.push(HadamardFactorTerm::InverseFactorial)
            }
            _ => (),
        }
        let (identity, _) = certify_candidate(
            image,
            other,
            &HadamardFactor { terms },
            &gauge,
            &HadamardFactorLiftOptions::default(),
        )?;
        if let HadamardFactorIdentityStatus::Certified { certificate } = identity {
            return Ok(PolynomialFamilyProductClosure::CoefficientMultiplier { right_is_multiplier,
                source:Box::new(relation.clone()),ratio_identity:format!("p(n,k) numerator(n,k) + q(n,k) denominator(n,k) = 0; binomial factors {menu:?}; anchor / multiplier(n,0)"),certificate });
        }
    }
    Err(unsupported(
        "no certified coefficient multiplier representation",
    ))
}
