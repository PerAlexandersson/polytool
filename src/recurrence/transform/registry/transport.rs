//! Recurrence transport, canonicalization and identity checks.

use super::*;

/// Transport a scalar rule through a coordinate map by the forward rule of
/// its substitution (module documentation, "Transport"): every piece of
/// the rule is pushed forward with the common bounds of the whole rule,
/// the additive term is mapped by the substitution and the left factor
/// `q`, powers of `lambda` are cleared and divided out, and the scale is
/// that of `lambda^e hom(f)`.
pub(super) fn transport_scalar(
    rule: &Recurrence,
    substitution: &Substitution,
    gauge: &Gauge,
) -> Result<Recurrence, TransformError> {
    let top = rule.max_offset();
    let index_shift = i64::try_from(top).map_err(|_| TransformError::Limit("offset".into()))?;
    let to_base = |p: &BivarPoly| Biv::from_bivar(p).shift_index(index_shift);
    // The pieces f(n + R) S^R and c(n + R) D^k S^(R - r) of the relation at
    // base n = m - R (the sign of a (-1)^m term is kept in its key).
    let lhs = rule
        .denominator
        .as_ref()
        .map(to_base)
        .unwrap_or_else(Biv::one);
    let mut pieces = vec![(PieceKey::Lhs, DOp::single(top, 0, lhs))];
    for term in &rule.terms {
        let coefficient = to_base(&term.coeff);
        if coefficient.is_zero() {
            continue;
        }
        let key = PieceKey::Term {
            offset: term.offset,
            derivative: term.deriv_order,
            alternating: term.sign == RecurrenceSign::AlternatingN,
        };
        pieces.push((
            key,
            DOp::single(top - term.offset, term.deriv_order, coefficient),
        ));
    }
    let bounds = pieces
        .iter()
        .map(|(_, piece)| PullbackBounds::of(piece))
        .fold(PullbackBounds::of(&DOp::identity()), |a, b| {
            PullbackBounds {
                top_shift: a.top_shift.max(b.top_shift),
                derivative_order: a.derivative_order.max(b.derivative_order),
                variable_degree: a.variable_degree.max(b.variable_degree),
            }
        });
    // push_forward_substitution with the bounds of the whole rule, so that
    // the left factor q is common to every piece.
    let inverse = substitution.inverse();
    let push = |op: &DOp| pull_back_substitution(op, &inverse, bounds).map_err(substitution_error);
    let from_base = |c: &Biv| bp_trim(c.shift_index(-index_shift).rows().to_vec());
    let mut combined: BTreeMap<PieceKey, Bp> = BTreeMap::new();
    for (key, piece) in &pieces {
        for ((shift, derivative), c) in &push(piece)?.0 {
            let image_key = match key {
                PieceKey::Lhs if (*shift, *derivative) == (top, 0) => PieceKey::Lhs,
                PieceKey::Term {
                    offset,
                    alternating,
                    ..
                } if *shift == top - offset => PieceKey::Term {
                    offset: *offset,
                    derivative: *derivative,
                    alternating: *alternating,
                },
                _ => {
                    return Err(TransformError::TransportUnsupported(
                        "the forward rule moved a term to another shift".into(),
                    ))
                }
            };
            let entry = combined.entry(image_key).or_default();
            *entry = bp_add(entry, &from_base(c));
        }
    }
    combined.retain(|key, poly| *key == PieceKey::Lhs || !bp_is_zero(poly));

    // The additive term: (P Q)_n = q C_n(g(n + R)), with
    // C_n(g) = lambda^(w + h n - tau deg g) hom(g) (index free by the check).
    let mut inhomogeneous: Option<(Bp, i64)> = None;
    if let Some(g) = &rule.inhomogeneous {
        let g = bp_from(g);
        if !bp_is_zero(&g) {
            let index_independent = gauge.step == 0
                || (gauge.lambda_is_constant()
                    && rational_power(&gauge.lambda[0], gauge.step).is_one());
            if !index_independent {
                return Err(TransformError::InhomogeneousNotTransportable);
            }
            let left_factor = push(&DOp::identity())?;
            let q_factor = left_factor
                .0
                .get(&(0, 0))
                .map(|c| bp_trim(c.rows().to_vec()))
                .unwrap_or_default();
            let (hom, degree) = bp_homogenize(&g, &gauge.numerator, &gauge.denominator);
            let degree = i64::try_from(degree)
                .ok()
                .and_then(|d| d.checked_mul(gauge.tau));
            let exponent = checked(degree.and_then(|d| gauge.weight_constant.checked_sub(d)))?;
            inhomogeneous = Some((bp_mul(&hom, &q_factor), exponent));
        }
    }
    let lambda_poly = vec_poly::trim(&gauge.lambda);
    if gauge.lambda_is_constant() {
        if let Some((g, exponent)) = inhomogeneous {
            let g = bp_scale(&g, &rational_power(&gauge.lambda[0], exponent));
            combined.insert(PieceKey::Inhomogeneous, g);
        }
    } else {
        if let Some((g, exponent)) = inhomogeneous {
            let deficit = usize::try_from((-exponent).max(0))
                .map_err(|_| TransformError::Limit("exponent".into()))?;
            if deficit > 0 {
                let factor = vec_poly::pow(&lambda_poly, deficit);
                for poly in combined.values_mut() {
                    *poly = bp_mul_x(poly, &factor);
                }
            }
            let surplus = usize::try_from(exponent.max(0))
                .map_err(|_| TransformError::Limit("exponent".into()))?;
            combined.insert(
                PieceKey::Inhomogeneous,
                bp_mul_x(&g, &vec_poly::pow(&lambda_poly, surplus)),
            );
        }
        // Divide out the largest common power of lambda (exact division).
        loop {
            let divided = combined
                .iter()
                .map(|(key, poly)| bp_div_linear(poly, &gauge.lambda).map(|p| (key.clone(), p)))
                .collect::<Option<BTreeMap<_, _>>>();
            match divided {
                Some(next) if next.get(&PieceKey::Lhs).is_some_and(|p| !bp_is_zero(p)) => {
                    combined = next;
                }
                _ => break,
            }
        }
    }
    normalize_left_factor(&mut combined, rule, gauge);
    Ok(assemble_rule(combined))
}

/// Scale the transported pieces so that the left factor is
/// `lambda^e hom(f)` (`hom(f) = den^(deg f) f(num/den)`, `e` an integer;
/// `hom(f) lambda^(-tau deg f)` for a constant `lambda`), the scale of the
/// chain-rule transport: the forward rule of a substitution fixes the left
/// factor only up to a nonzero constant.  Unchanged when the left factor is
/// not of that form (it always is).
fn normalize_left_factor(combined: &mut BTreeMap<PieceKey, Bp>, rule: &Recurrence, gauge: &Gauge) {
    let f = rule
        .denominator
        .as_ref()
        .map(bp_from)
        .unwrap_or_else(|| bp_constant(Q::one()));
    let (hom, degree) = bp_homogenize(&f, &gauge.numerator, &gauge.denominator);
    let Some(lhs) = combined.get(&PieceKey::Lhs).cloned() else {
        return;
    };
    let (expected, candidate) = if gauge.lambda_is_constant() {
        let exponent = -(degree as i64) * gauge.tau;
        (
            bp_scale(&hom, &rational_power(&gauge.lambda[0], exponent)),
            lhs,
        )
    } else {
        let lambda_poly = vec_poly::trim(&gauge.lambda);
        let difference =
            bp_x_degree(&lhs).unwrap_or(0) as i64 - bp_x_degree(&hom).unwrap_or(0) as i64;
        let power = |e: i64| vec_poly::pow(&lambda_poly, e.unsigned_abs() as usize);
        (
            bp_mul_x(&hom, &power(difference.max(0))),
            bp_mul_x(&lhs, &power(difference.min(0))),
        )
    };
    let position = expected
        .iter()
        .enumerate()
        .find_map(|(i, row)| row.iter().position(|c| !c.is_zero()).map(|j| (i, j)));
    let Some((i, j)) = position else {
        return;
    };
    let Some(scale) = candidate
        .get(i)
        .and_then(|row| row.get(j))
        .filter(|c| !c.is_zero())
        .map(|c| c / &expected[i][j])
    else {
        return;
    };
    if bp_scale(&expected, &scale) != candidate {
        return;
    }
    let inverse = scale.recip();
    for poly in combined.values_mut() {
        *poly = bp_scale(poly, &inverse);
    }
}

pub(in crate::recurrence::transform) fn assemble_rule(
    mut combined: BTreeMap<PieceKey, Bp>,
) -> Recurrence {
    let mut lhs = combined.remove(&PieceKey::Lhs).unwrap_or_default();
    if let Some(constant) = bp_is_constant(&lhs) {
        let inverse = constant.recip();
        for poly in combined.values_mut() {
            *poly = bp_scale(poly, &inverse);
        }
        lhs = bp_constant(Q::one());
    }
    let inhomogeneous = combined
        .remove(&PieceKey::Inhomogeneous)
        .filter(|g| !bp_is_zero(g))
        .map(bp_to_bivar);
    let terms = combined
        .into_iter()
        .filter(|(_, poly)| !bp_is_zero(poly))
        .map(|(key, poly)| match key {
            PieceKey::Term {
                offset,
                derivative,
                alternating,
            } => RecurrenceTerm::new(offset, derivative, bp_to_bivar(poly)).with_sign(
                if alternating {
                    RecurrenceSign::AlternatingN
                } else {
                    RecurrenceSign::None
                },
            ),
            _ => unreachable!("only term pieces remain"),
        })
        .collect();
    Recurrence {
        terms,
        denominator: (bp_is_constant(&lhs) != Some(Q::one())).then(|| bp_to_bivar(lhs)),
        inhomogeneous,
    }
}

pub(in crate::recurrence::transform) fn canonical_pieces(
    rule: &Recurrence,
) -> BTreeMap<PieceKey, Bp> {
    let mut combined: BTreeMap<PieceKey, Bp> = BTreeMap::new();
    combined.insert(
        PieceKey::Lhs,
        rule.denominator
            .as_ref()
            .map(bp_from)
            .unwrap_or_else(|| bp_constant(Q::one())),
    );
    for term in &rule.terms {
        let key = PieceKey::Term {
            offset: term.offset,
            derivative: term.deriv_order,
            alternating: term.sign == RecurrenceSign::AlternatingN,
        };
        let entry = combined.entry(key).or_default();
        *entry = bp_add(entry, &bp_from(&term.coeff));
    }
    if let Some(g) = &rule.inhomogeneous {
        combined.insert(PieceKey::Inhomogeneous, bp_from(g));
    }
    if let Some(constant) = combined.get(&PieceKey::Lhs).and_then(bp_is_constant) {
        let inverse = constant.recip();
        for poly in combined.values_mut() {
            *poly = bp_scale(poly, &inverse);
        }
    }
    combined.retain(|key, poly| *key == PieceKey::Lhs || !bp_is_zero(poly));
    combined
}

/// Canonical form of a rule: terms with equal offset, derivative order and
/// sign merged, zero terms dropped, coefficients trimmed and a constant left
/// factor divided out.
pub fn canonical_recurrence(rule: &Recurrence) -> Recurrence {
    assemble_rule(canonical_pieces(rule))
}

/// Whether two rules have the same canonical form (identical recursions).
pub fn same_canonical_recurrence(left: &Recurrence, right: &Recurrence) -> bool {
    canonical_pieces(left) == canonical_pieces(right)
}

// ---------------------------------------------------------------------------
// Symmetric decomposition transport
// ---------------------------------------------------------------------------

pub(super) fn transport_symmetric_decomposition(
    rule: &Recurrence,
    schedule: &DegreeSchedule,
) -> Result<SharedDenominatorVectorRecurrence, TransformError> {
    if rule.terms.iter().any(|t| t.offset != 1) {
        return Err(TransformError::TransportUnsupported(
            "the symmetric decomposition transports first-order rules (all offsets 1)".into(),
        ));
    }
    if rule.terms.iter().any(RecurrenceTerm::is_alternating) {
        return Err(TransformError::TransportUnsupported(
            "vector systems carry no (-1)^n factors".into(),
        ));
    }
    let reversal = RecurrenceTransformation::Reversal {
        schedule: schedule.clone(),
    };
    let reversed = match reversal.transport(rule)? {
        TransportedRecurrence::Scalar(r) => r,
        TransportedRecurrence::Vector(_) => unreachable!("reversal is scalar"),
    };
    let source = canonical_recurrence(rule);
    let one = bp_constant(Q::one());
    let f = source
        .denominator
        .as_ref()
        .map(bp_from)
        .unwrap_or_else(|| one.clone());
    let f_bar = reversed
        .denominator
        .as_ref()
        .map(bp_from)
        .unwrap_or_else(|| one.clone());
    let one_minus_x = [q(1), q(-1)];
    let x = [q(0), q(1)];
    let (denominator, u, u_bar) = if f == f_bar {
        (bp_mul_x(&f, &one_minus_x), one.clone(), one)
    } else {
        (bp_mul_x(&bp_mul(&f, &f_bar), &one_minus_x), f_bar, f)
    };
    let coefficients = |r: &Recurrence| -> Vec<Bp> {
        let order = r.terms.iter().map(|t| t.deriv_order).max().unwrap_or(0);
        let mut c = vec![Vec::new(); order + 2];
        for term in &r.terms {
            c[term.deriv_order] = bp_add(&c[term.deriv_order], &bp_from(&term.coeff));
        }
        c
    };
    let mut c = coefficients(&source);
    let mut c_bar = coefficients(&reversed);
    let width = c.len().max(c_bar.len());
    c.resize(width, Vec::new());
    c_bar.resize(width, Vec::new());
    let order_count = width - 1;
    let mut matrix = vec![vec![vec![Vec::new(); order_count]; 2]; 2];
    for k in 0..order_count {
        let uc = bp_mul(&u, &c[k]);
        let ubc = bp_mul(&u_bar, &c_bar[k]);
        let lowered = bp_scale(&bp_mul(&u, &c[k + 1]), &q((k + 1) as i64));
        let x_uc = bp_mul_x(&uc, &x);
        let x_ubc = bp_mul_x(&ubc, &x);
        matrix[0][0][k] = bp_add(&uc, &bp_scale(&x_ubc, &q(-1)));
        matrix[0][1][k] = bp_add(&bp_add(&x_uc, &bp_scale(&x_ubc, &q(-1))), &lowered);
        matrix[1][0][k] = bp_add(&ubc, &bp_scale(&uc, &q(-1)));
        matrix[1][1][k] = bp_add(
            &bp_add(&ubc, &bp_scale(&x_uc, &q(-1))),
            &bp_scale(&lowered, &q(-1)),
        );
    }
    let g = source
        .inhomogeneous
        .as_ref()
        .map(bp_from)
        .unwrap_or_default();
    let g_bar = reversed
        .inhomogeneous
        .as_ref()
        .map(bp_from)
        .unwrap_or_default();
    let forcing = (!bp_is_zero(&g) || !bp_is_zero(&g_bar)).then(|| {
        let ug = bp_mul(&u, &g);
        let ubg = bp_mul(&u_bar, &g_bar);
        vec![
            bp_add(&ug, &bp_scale(&bp_mul_x(&ubg, &x), &q(-1))),
            bp_add(&ubg, &bp_scale(&ug, &q(-1))),
        ]
    });
    // Coefficients are in the target index m = n + 1; the vector convention uses n.
    let to_source_index = |p: Bp| bp_to_bivar(p).shift_index(1);
    Ok(SharedDenominatorVectorRecurrence {
        denominator: to_source_index(denominator),
        matrix: matrix
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|entry| {
                        WeylOperator::new(entry.into_iter().map(to_source_index).collect())
                    })
                    .collect()
            })
            .collect(),
        forcing: forcing.map(|components| components.into_iter().map(to_source_index).collect()),
    })
}

// ---------------------------------------------------------------------------
// Residual checks
// ---------------------------------------------------------------------------

/// Whether the scalar rule holds exactly at target index `target`, where
/// `rows[i]` is the row with original index `origin + i`. `None` when a
/// referenced row is unavailable.
pub fn transformed_rule_holds_at(
    rule: &Recurrence,
    rows: &[Vec<BigRational>],
    origin: usize,
    target: usize,
) -> Option<bool> {
    let position = target.checked_sub(origin)?;
    let current = rows.get(position)?;
    let lhs = match &rule.denominator {
        Some(f) => poly_mul_rational(&bivar_eval_n(f, target), current),
        None => current.clone(),
    };
    let mut rhs = rational_zero_poly();
    for term in &rule.terms {
        let source = rows.get(position.checked_sub(term.offset)?)?;
        let derivative = poly_nth_derivative_rational(source, term.deriv_order);
        let product = poly_mul_rational(&bivar_eval_n(&term.coeff, target), &derivative);
        let sign = if term.is_alternating() && target % 2 == 1 {
            q(-1)
        } else {
            Q::one()
        };
        poly_add_scaled_assign(&mut rhs, &product, &sign);
    }
    if let Some(g) = &rule.inhomogeneous {
        poly_add_scaled_assign(&mut rhs, &bivar_eval_n(g, target), &Q::one());
    }
    Some(vec_poly::trim(&lhs) == vec_poly::trim(&rhs))
}
