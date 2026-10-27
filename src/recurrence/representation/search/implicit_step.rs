//! Implicit-relation steps: application of a map chain to an implicit relation and
//! transport of its conditions.

use super::*;

/// `c / lambda(n + shift, k)` for a diagonal weight, as `(numerator(n),
/// denominator(n))` with `lambda = K (k!)^p prod f_j(n,k)^(e_j)`.
fn inverse_weight(weight: &DiagonalWeight, shift: usize, k: usize) -> Option<(Vec<Q>, Vec<Q>)> {
    let kf = crate::recurrence::implicit::algebra::factorial(k);
    let k_value = Q::from_integer(BigInt::from(k));
    let constant = |c: Q| Some((vec![c], vec![Q::from_integer(1.into())]));
    match weight {
        DiagonalWeight::FactorialPower { power } => {
            let mut value = Q::from_integer(1.into());
            for _ in 0..*power {
                value *= &kf;
            }
            constant(value)
        }
        DiagonalWeight::Binomial { m } => constant(binomial(k + m, *m).recip()),
        DiagonalWeight::InverseBinomial { m } => constant(binomial(k + m, *m)),
        DiagonalWeight::Product {
            scale,
            factorial_exponent,
            factors,
        } => {
            let mut value = parse_rational_coeff(scale).ok()?.recip();
            for _ in 0..factorial_exponent.unsigned_abs() {
                if *factorial_exponent > 0 {
                    value /= &kf;
                } else {
                    value *= &kf;
                }
            }
            let mut numerator = vec![value];
            let mut denominator = vec![Q::from_integer(1.into())];
            for factor in factors {
                let f = biv_from_json(&factor.coeffs).ok()?;
                let at = f
                    .shift(&Q::from_integer(BigInt::from(shift)), &Q::zero())
                    .eval_y(&k_value);
                let power = vec_poly::pow(&at, factor.exponent.unsigned_abs() as usize);
                if factor.exponent > 0 {
                    denominator = vec_poly::mul(&denominator, &power);
                } else {
                    numerator = vec_poly::mul(&numerator, &power);
                }
            }
            Some((numerator, denominator))
        }
    }
}

/// A coefficient condition of the parent rows `P` rewritten for `Q = T(P)`,
/// or `None` when the map forgets a coefficient it refers to (or is not an
/// implicit map).
pub(super) fn transport_condition(
    t: &RecurrenceTransformation,
    condition: &CoefficientCondition,
) -> Option<CoefficientCondition> {
    let condition = condition.canonical();
    let transported = match t {
        RecurrenceTransformation::Derivative { order } => {
            push_condition_through_derivative(&condition, *order)?
        }
        RecurrenceTransformation::DerivativeAntecedent { g } => {
            // [x^e] P = binom(e+m,m) sum_j g_j [x^(e+m-j)] Q.
            let (g, m) = derivative_antecedent_parameters(g).ok()?;
            let mut terms = Vec::new();
            for term in &condition.terms {
                let e = term.coefficient_index;
                for (j, gj) in g.iter().enumerate() {
                    if gj.is_zero() || e + m < j {
                        continue;
                    }
                    terms.push(ConditionTerm {
                        shift: term.shift,
                        coefficient_index: e + m - j,
                        coeff: vec_poly::scale(&term.coeff, &(binomial(e + m, m) * gj)),
                    });
                }
            }
            CoefficientCondition { terms }.canonical()
        }
        RecurrenceTransformation::Diagonal { weight } => {
            // [x^k] P_(n+i) = [x^k] Q_(n+i) / lambda(n+i, k); clear the
            // polynomial denominators by their lcm (nonzero on the domain by
            // the weight hypothesis).
            let mut raw = Vec::new();
            for term in &condition.terms {
                let (numerator, denominator) =
                    inverse_weight(weight, term.shift, term.coefficient_index)?;
                raw.push((term, vec_poly::mul(&term.coeff, &numerator), denominator));
            }
            let lcm = raw.iter().fold(
                Polynomial::new(vec![Q::from_integer(1.into())]),
                |acc, (_, _, d)| {
                    let d = Polynomial::new(d.clone());
                    let g = acc.gcd(&d);
                    let (quotient, _) = (acc.clone() * d).div_rem(&g);
                    quotient
                },
            );
            let terms = raw
                .into_iter()
                .map(|(term, numerator, denominator)| {
                    let (factor, _) = lcm.div_rem(&Polynomial::new(denominator));
                    ConditionTerm {
                        shift: term.shift,
                        coefficient_index: term.coefficient_index,
                        coeff: vec_poly::trim(&vec_poly::mul(&numerator, factor.coeffs())),
                    }
                })
                .collect();
            CoefficientCondition { terms }.canonical()
        }
        _ => return None,
    };
    // Scale the leading term's leading coefficient to one.
    let lead = transported
        .terms
        .last()
        .and_then(|t| vec_poly::degree(&t.coeff).map(|d| t.coeff[d].clone()));
    Some(match lead {
        Some(lead) if !lead.is_zero() => CoefficientCondition {
            terms: transported
                .terms
                .into_iter()
                .map(|t| ConditionTerm {
                    coeff: vec_poly::scale(&t.coeff, &lead.recip()),
                    ..t
                })
                .collect(),
        },
        _ => transported,
    })
}

fn push_unique(list: &mut Vec<DescribedCondition>, record: DescribedCondition) {
    if !list.iter().any(|c| c.condition == record.condition) {
        list.push(record);
    }
}

/// Transport the single relation of a certified parent description (and its
/// conditions) through an implicit map.
pub(super) fn build_implicit_relation_step(
    t: &RecurrenceTransformation,
    parent: &RepresentationDescription,
    parent_rows: &[Vec<Q>],
    rows: &[Vec<Q>],
    origin: usize,
) -> Result<(ImplicitRelationStep, Vec<TransformConditionStatus>), RecurrenceError> {
    let description = parent
        .description
        .as_ref()
        .ok_or(RecurrenceError::InvalidInput(
            "the parent description has no generating description".into(),
        ))?;
    let [source] = description.relations.as_slice() else {
        return Err(RecurrenceError::InvalidInput(
            "the parent description needs exactly one relation".into(),
        ));
    };
    let relation = ImplicitRelation::from_json(&source.relation)?;
    let map_conditions = match t {
        RecurrenceTransformation::Diagonal { weight } => weight_regularity(weight, origin)?,
        RecurrenceTransformation::DerivativeAntecedent { .. } => vec![TransformConditionStatus {
            condition: "T_g bijective".into(),
            status: CERTIFIED_ALL_INDICES.into(),
            detail: "the monomial matrix of T_g is triangular with diagonal lc(g) binom(k+m,m) \
                     != 0"
                .into(),
        }],
        RecurrenceTransformation::Derivative { order } => vec![TransformConditionStatus {
            condition: format!(
                "Q_n = D^{order} P_n with lost data {} P_n",
                crate::recurrence::implicit::lost_range(*order)
            ),
            status: CERTIFIED_ALL_INDICES.into(),
            detail: "defined for every row; the lost coefficients are recorded as boundary data"
                .into(),
        }],
        other => {
            return Err(RecurrenceError::InvalidInput(format!(
                "implicit relations are not transported through {}",
                other.kind_name()
            )))
        }
    };
    let map_certified = map_conditions
        .iter()
        .all(|c| c.status == CERTIFIED_ALL_INDICES);
    let certified = evidence_is_certified(&source.evidence) && map_certified;
    let evidence = if certified {
        TRANSPORTED_IDENTITY_EVIDENCE
    } else {
        CHECKED_SUPPLIED_ROWS
    };
    let first_base = source.first_base;
    let transport = t.transport_implicit(&relation, first_base)?;
    let mut exceptional = Vec::new();
    for (k, condition) in &transport.exceptional_conditions {
        push_unique(
            &mut exceptional,
            DescribedCondition::new(condition, first_base, evidence, Some(*k), "Q").with_note(
                format!("[x^{k}] of the uncleared residual: H(n,{k}) = 0 at the exceptional index"),
            ),
        );
    }
    let mut boundary = Vec::new();
    let mut dropped = Vec::new();
    let carried_evidence = |record: &DescribedCondition| {
        if evidence_is_certified(&record.evidence) && map_certified {
            record.evidence.clone()
        } else {
            CHECKED_SUPPLIED_ROWS.to_string()
        }
    };
    let carried_index = |k: Option<usize>| match t {
        RecurrenceTransformation::Diagonal { .. } => k,
        RecurrenceTransformation::Derivative { order } => k.and_then(|k| k.checked_sub(*order)),
        _ => None,
    };
    for (list, records, target) in [
        (
            "exceptional_index_conditions",
            &description.exceptional_index_conditions,
            true,
        ),
        (
            "boundary_conditions",
            &description.boundary_conditions,
            false,
        ),
    ] {
        for (i, record) in records.iter().enumerate() {
            let condition = CoefficientCondition::from_json(&record.condition)?;
            match transport_condition(t, &condition) {
                Some(carried) if !carried.is_vacuous() => {
                    let described = DescribedCondition::new(
                        &carried,
                        record.first_base,
                        &carried_evidence(record),
                        carried_index(record.coefficient_index),
                        "Q",
                    )
                    .with_note(format!("the parent's {list}[{i}] carried through the map"));
                    push_unique(
                        if target {
                            &mut exceptional
                        } else {
                            &mut boundary
                        },
                        described,
                    );
                }
                Some(_) => {}
                None => dropped.push(format!("{list}[{i}]")),
            }
        }
    }
    let mut lost_conditions = Vec::new();
    if let RecurrenceTransformation::Derivative { order } = t {
        let pulled = pull_back_relation_through_derivative(&transport.relation, *order);
        let indices = valuation_exceptional_indices(&pulled, first_base)
            .unwrap_or_default()
            .into_iter()
            .filter(|k| *k <= MAX_CONDITION_COEFFICIENT_INDEX)
            .collect::<Vec<_>>();
        lost_conditions = coefficient_conditions_of(source, &indices, "P")?;
        for record in &lost_conditions {
            let condition = CoefficientCondition::from_json(&record.condition)?;
            if let Some(pushed) = push_condition_through_derivative(&condition, *order) {
                if !pushed.is_vacuous() {
                    push_unique(
                        &mut boundary,
                        DescribedCondition::new(
                            &pushed,
                            record.first_base,
                            &carried_evidence(record),
                            carried_index(record.coefficient_index),
                            "Q",
                        )
                        .with_note("[x^k] of the parent relation, free of lost coefficients"),
                    );
                }
            }
        }
        lost_conditions.extend(description.exceptional_index_conditions.iter().cloned());
        lost_conditions.extend(description.boundary_conditions.iter().cloned());
    }
    let r = transport.relation.max_shift();
    let seed_end = (first_base + r).max(origin + 1);
    if rows.len() < seed_end - origin {
        return Err(RecurrenceError::InvalidInput(
            "too few rows for the seeds of the transported relation".into(),
        ));
    }
    let note = match t {
        RecurrenceTransformation::Derivative { .. } => {
            "transported through the lossy derivative map (a consequence of the parent relation)"
        }
        _ if transport.clearing.is_some() => {
            "transported through the diagonal map and cleared by one common multiplier"
        }
        _ => "transported by substitution P = T_g Q",
    };
    let described = DescribedRelation::new(&transport.relation, first_base, evidence, "Q")
        .with_note(format!("{note}; source: {}", source.text));
    let q_description = GeneratingDescription {
        first_index: origin,
        seeds: rows[..seed_end - origin]
            .iter()
            .map(|r| encode_poly(r))
            .collect(),
        relations: vec![described.clone()],
        exceptional_index_conditions: exceptional,
        boundary_conditions: boundary,
    };
    let assessment = assess_generating_description(&q_description)?;
    let lossy = match t {
        RecurrenceTransformation::Derivative { order } => Some(lossy_reconstruction(
            *order,
            &q_description,
            origin,
            parent_rows,
            &lost_conditions,
        )?),
        _ => None,
    };
    let (holds, checked) = description_holds(&q_description, rows, origin)?;
    Ok((
        ImplicitRelationStep {
            source_relation: source.clone(),
            relation: described,
            clearing: transport.clearing,
            certificate: transport.reconstruction,
            dropped_conditions: dropped,
            description: q_description,
            assessment,
            lossy,
            holds_on_rows: holds,
            checked_bases: checked,
        },
        map_conditions,
    ))
}

/// Whether every relation and condition of a description holds on the rows
/// at each base of its domain, and the number of relation bases checked.
fn description_holds(
    description: &GeneratingDescription,
    rows: &[Vec<Q>],
    origin: usize,
) -> Result<(bool, usize), RecurrenceError> {
    let mut holds = true;
    let mut checked = 0;
    for record in &description.relations {
        let relation = ImplicitRelation::from_json(&record.relation)?;
        let end = (origin + rows.len()).saturating_sub(relation.max_shift());
        for base in record.first_base.max(origin)..end {
            holds &= relation.holds_at(rows, origin, base) != Some(false);
            checked += 1;
        }
    }
    for record in description
        .exceptional_index_conditions
        .iter()
        .chain(&description.boundary_conditions)
    {
        let condition = CoefficientCondition::from_json(&record.condition)?;
        let end = (origin + rows.len()).saturating_sub(condition.max_shift());
        for base in record.first_base.max(origin)..end {
            holds &= condition.holds_at(rows, origin, base) != Some(false);
        }
    }
    Ok((holds, checked))
}

/// Independent checks of an implicit-relation step on the new rows.
pub(super) fn verify_implicit_relation_step(
    t: &RecurrenceTransformation,
    record: &ImplicitRelationStep,
    rows: &[Vec<Q>],
    origin: usize,
) -> Result<(), RecurrenceError> {
    let source = ImplicitRelation::from_json(&record.source_relation.relation)?;
    let transport = t.transport_implicit(&source, record.source_relation.first_base)?;
    let relation = ImplicitRelation::from_json(&record.relation.relation)?;
    if !relation.same_relation(&transport.relation) {
        return Err(RecurrenceError::VerificationFailed(
            "the transported relation does not replay".into(),
        ));
    }
    if let Some(certificate) = &record.certificate {
        certificate.verify_identity()?;
    }
    if !record.holds_on_rows || !description_holds(&record.description, rows, origin)?.0 {
        return Err(RecurrenceError::VerificationFailed(
            "the transported description fails on the rows".into(),
        ));
    }
    if assess_generating_description(&record.description)? != record.assessment {
        return Err(RecurrenceError::VerificationFailed(
            "the assessment does not replay".into(),
        ));
    }
    if let Some(lossy) = &record.lossy {
        if assess_generating_description(&lossy.with_lost_data)? != lossy.with_lost_data_assessment
            || assess_generating_description(&lossy.without_lost_data)?
                != lossy.without_lost_data_assessment
        {
            return Err(RecurrenceError::VerificationFailed(
                "the lossy reconstruction assessments do not replay".into(),
            ));
        }
    }
    Ok(())
}
