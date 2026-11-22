//! Supplied-candidate checks, exact fitting and bounded coordinate proposals.
use super::*;
use preservation::{check_operator_preserver, directed};

fn validate(request: &InterlacingRefinementRequest) -> Result<(), RecurrenceError> {
    let l = &request.limits;
    if l.maximum_dimension == 0
        || l.maximum_dimension > 5
        || l.maximum_candidates == 0
        || l.maximum_candidates > 256
        || l.variable_degree > 3
        || l.index_degree > 4
        || l.derivative_order > 2
        || l.maximum_fit_cells == 0
        || l.maximum_fit_cells > 2_000_000
        || l.maximum_supports == 0
        || l.maximum_supports > 64
        || request.candidates.len() > l.maximum_candidates
        || request.states.len() > 128
        || request.target_rows.len() > 128
        || request.target_rows.iter().any(|p| p.len() > 65)
    {
        return Err(invalid("refinement search bounds exceed supported limits"));
    }
    if usize::from(request.source.is_some())
        + usize::from(request.recurrence.is_some())
        + usize::from(!request.states.is_empty())
        != 1
    {
        return Err(invalid(
            "supply exactly one of recurrence, source or finite states",
        ));
    }
    if request.first_index > usize::MAX - 256 {
        return Err(invalid("first index too large"));
    }
    if request.states.is_empty() && !request.target_rows.is_empty() {
        return Err(invalid("target_rows are supported with finite states only"));
    }
    if !request.states.is_empty() && !request.source_projection.is_empty() {
        return Err(invalid("source_projection requires a defining source"));
    }
    Ok(())
}
fn source_weights(
    request: &InterlacingRefinementRequest,
    n: usize,
) -> Result<Vec<Q>, RecurrenceError> {
    if request.source_projection.is_empty() {
        let mut w = vec![Q::zero(); n];
        w[n - 1] = Q::one();
        Ok(w)
    } else {
        let w = numbers(&request.source_projection)?;
        if w.len() != n || w.iter().all(Zero::is_zero) {
            return Err(invalid("source projection dimensions differ"));
        }
        Ok(w)
    }
}
fn weights(
    candidate: &InterlacingRefinementCandidate,
    n: usize,
) -> Result<Vec<Q>, RecurrenceError> {
    let w = if candidate.projection.is_empty() {
        vec![Q::one(); n]
    } else {
        numbers(&candidate.projection)?
    };
    if w.len() != n || w.iter().any(Signed::is_negative) || w.iter().all(Zero::is_zero) {
        return Err(invalid(
            "target projection must have matching nonnegative scalar weights and be nonzero",
        ));
    }
    Ok(w)
}
fn ordered(state: &[Vec<Q>]) -> bool {
    state.iter().all(|p| {
        p.iter().all(|c| !c.is_negative())
            && crate::is_real_rooted_bigint_coeffs(&vp::clear_denominators(p))
    }) && (0..state.len()).all(|i| (i + 1..state.len()).all(|j| directed(&state[i], &state[j])))
}
fn degrees(state: &[Vec<Q>]) -> (Vec<Option<usize>>, bool) {
    let d: Vec<_> = state.iter().map(|p| vp::degree(p)).collect();
    let valid = (0..d.len()).all(|i| {
        (i + 1..d.len()).all(|j| match (d[i], d[j]) {
            (Some(a), Some(b)) => a == b || a + 1 == b,
            _ => true,
        })
    });
    (d, valid)
}
fn sum(state: &[Vec<Q>], weights: &[Q]) -> Vec<Q> {
    state
        .iter()
        .zip(weights)
        .fold(vec![Q::zero()], |s, (p, w)| vp::add(&s, &vp::scale(p, w)))
}

fn check(
    candidate: &InterlacingRefinementCandidate,
    source: Option<&TransformVectorRecurrenceJson>,
    states: &[Vec<Vec<Q>>],
    targets: &[Vec<Q>],
    source_projection: &[Q],
    first: usize,
    limits: &RefinementSearchLimits,
) -> Result<InterlacingRefinementCheck, RecurrenceError> {
    let input_dimension = states[0].len();
    let dimension = candidate
        .coordinates
        .as_ref()
        .map_or(input_dimension, Vec::len);
    if dimension == 0 || dimension > limits.maximum_dimension {
        return Err(invalid("refinement dimension exceeds limits"));
    }
    let t = match &candidate.coordinates {
        Some(t) => decode_coordinates(t, dimension, input_dimension)?,
        None => identity(input_dimension),
    };
    let refined: Vec<_> = states
        .iter()
        .enumerate()
        .map(|(i, s)| apply(&t, s, first + i))
        .collect();
    let projection = weights(candidate, dimension)?;
    let (initial_degrees, degree_conditions) = degrees(&refined[0]);
    let mut check = InterlacingRefinementCheck {
        candidate: candidate.clone(),
        identity: "unresolved".into(),
        target_recovery: "unresolved".into(),
        initial_invariant: ordered(&refined[0]),
        initial_degrees,
        degree_conditions,
        transition: None,
        preservation: None,
        status: "unresolved".into(),
        reason: String::new(),
    };
    let recovered = (0..input_dimension).all(|j| {
        let op = (0..dimension).fold(Op::default(), |s, i| s.add(&t[i][j].scale(&projection[i])));
        op == Op::identity().scale(&source_projection[j])
    });
    if source.is_some() {
        check.target_recovery = if recovered { "exact" } else { "failed" }.into();
    } else if !targets.is_empty() {
        check.target_recovery = if targets.len() == refined.len()
            && refined
                .iter()
                .zip(targets)
                .all(|(s, p)| vp::trim(&sum(s, &projection)) == vp::trim(p))
        {
            "finite_fit"
        } else {
            "failed"
        }
        .into();
    } else {
        check.target_recovery = "unresolved".into();
    }
    let mut denominator = Biv::one();
    let matrix = if let Some(transition) = &candidate.transition {
        let (m, q, seed) = decode_system(transition)?;
        if m.len() != dimension
            || transition.first_index != first
            || seed.iter().map(|p| vp::trim(p)).collect::<Vec<_>>()
                != refined[0].iter().map(|p| vp::trim(p)).collect::<Vec<_>>()
        {
            check.status = "rejected".into();
            check.reason =
                "supplied transfer has wrong dimensions, first index or initial vector".into();
            return Ok(check);
        }
        denominator = q;
        Some(m)
    } else if let Some(source) = source {
        let (c, q, _) = decode_system(source)?;
        denominator = q;
        match fit_symbolic(&t, &c, limits) {
            Ok(matrix) => matrix,
            Err(RecurrenceError::InvalidInput(reason))
                if reason.contains("fit exceeds maximum_fit_cells") =>
            {
                check.reason = reason;
                return Ok(check);
            }
            Err(error) => return Err(error),
        }
    } else {
        let n = dimension;
        let unknowns = n * (limits.variable_degree + 1) * (limits.index_degree + 1);
        let max_width =
            refined.iter().flatten().map(Vec::len).max().unwrap_or(1) + limits.variable_degree;
        if unknowns
            .saturating_mul(max_width)
            .saturating_mul(refined.len())
            > limits.maximum_fit_cells
        {
            check.reason = "finite fit exceeds maximum_fit_cells".into();
            return Ok(check);
        }
        let opts = VectorRecurrenceOptions {
            var_deg: limits.variable_degree,
            idx_deg: limits.index_degree,
            diff_deg: 0,
            held_out_transitions: usize::from(refined.len() > 3),
            require_unique: false,
            ..Default::default()
        };
        find_vector_recurrence_rational(&refined, first, &opts)
            .ok()
            .map(|fit| {
                fit.recurrence
                    .matrix
                    .iter()
                    .map(|row| {
                        row.iter()
                            .map(|op| {
                                let mut out = Op::default();
                                for (d, c) in op.coefficients.iter().enumerate() {
                                    out.add_term((0, d), Biv::from_bivar(c));
                                }
                                out
                            })
                            .collect()
                    })
                    .collect()
            })
    };
    let Some(matrix) = matrix else {
        if !check.initial_invariant || check.target_recovery == "failed" {
            check.status = "rejected".into();
        }
        check.reason = "no exact transfer found under the degree and work caps".into();
        return Ok(check);
    };
    let transition = encode_system(&matrix, &denominator, &refined[0], first);
    check.transition = Some(transition.clone());
    if let Some(source) = source {
        let (c, q, _) = decode_system(source)?;
        let shifted = t
            .iter()
            .map(|r| r.iter().map(|op| op.shift_index(1)).collect())
            .collect();
        let left = multiply(&shifted, &c);
        let right = multiply(&matrix, &t);
        let exact = left
            .iter()
            .flatten()
            .zip(right.iter().flatten())
            .all(|(l, r)| l.left_multiply(&denominator) == r.left_multiply(&q));
        check.identity = if exact { "exact" } else { "failed" }.into();
    } else {
        let (system, _, _) = transition.to_recurrence()?;
        check.identity = if refined
            .windows(2)
            .enumerate()
            .all(|(i, w)| system.holds_transition_rational(&w[0], &w[1], first + i))
        {
            "finite_fit"
        } else {
            "failed"
        }
        .into();
    }
    let mut preservation = check_operator_preserver(&matrix, first)?;
    if !denominator_positive(&denominator, first) {
        preservation.status = "unresolved".into();
        preservation.reason =
            "positive common denominator on the entire integer ray has not been certified".into();
    }
    if let Some(source) = source {
        let (_, q, _) = decode_system(source)?;
        if !denominator_positive(&q, first) {
            preservation.status = "unresolved".into();
            preservation.reason =
                "positive source denominator has not been certified throughout the integer ray"
                    .into();
        }
    }
    let preserved = preservation.status == "proved";
    let rejected = preservation.status == "rejected"
        || check.identity == "failed"
        || check.target_recovery == "failed"
        || !check.initial_invariant;
    check.reason=if rejected{"candidate fails an explicit algebraic, initial-order or preservation obligation"}else if preserved && check.identity=="exact" && check.target_recovery=="exact"{"exact intertwiner, scalar target recovery, ordered PF seed and all-index preservation; target is real-rooted or zero from the recorded first index"}else{"finite fitting, unresolved target identity or unresolved preservation remains; no infinite-family conclusion"}.into();
    check.status = if rejected {
        "rejected"
    } else if preserved && check.identity == "exact" && check.target_recovery == "exact" {
        "proved"
    } else if check.identity == "finite_fit" {
        "finite_fit"
    } else {
        "unresolved"
    }
    .into();
    check.preservation = Some(preservation);
    Ok(check)
}

fn generated(t: OperatorMatrix, name: String, projection: &[Q]) -> InterlacingRefinementCandidate {
    InterlacingRefinementCandidate {
        name,
        coordinates: Some(encode_coordinates(&t)),
        transition: None,
        projection: projection.iter().map(ToString::to_string).collect(),
    }
}
fn scaled_template(
    source: &TransformVectorRecurrenceJson,
) -> Result<Option<(ScaledRefinementInterval, InterlacingRefinementCandidate)>, RecurrenceError> {
    let (m, q, _) = decode_system(source)?;
    if m.len() != 3 || q != Biv::one() {
        return Ok(None);
    }
    let get = |j: usize, x: usize| {
        m[2][j]
            .0
            .get(&(0, 0))
            .map(|p| {
                p.rows()
                    .first()
                    .and_then(|r| r.get(x))
                    .cloned()
                    .unwrap_or_default()
            })
            .unwrap_or_default()
    };
    let p = [
        get(2, 0).to_string(),
        get(2, 1).to_string(),
        get(0, 2).to_string(),
        get(1, 1).to_string(),
    ];
    let expected = vec![
        vec![Op::default(), Op::identity(), Op::default()],
        vec![Op::default(), Op::default(), Op::identity()],
        vec![
            Op::single(0, 0, Biv::new(vec![vec![Q::zero(), Q::zero(), get(0, 2)]])),
            Op::single(0, 0, Biv::new(vec![vec![Q::zero(), get(1, 1)]])),
            Op::single(0, 0, Biv::new(vec![vec![get(2, 0), get(2, 1)]])),
        ],
    ];
    if m != expected {
        return Ok(None);
    }
    let interval = solve_scaled_refinement_interval(&p)?;
    let Some(k) = &interval.rational_scale else {
        return Ok(Some((
            interval,
            generated(
                identity(3),
                "infeasible scaled lag-sum template".into(),
                &[Q::zero(), Q::zero(), Q::one()],
            ),
        )));
    };
    let k = parse_rational_coeff(k)?;
    let t = vec![
        vec![Op::default(), Op::identity(), Op::default()],
        vec![
            Op::single(0, 0, Biv::new(vec![vec![Q::zero(), k.clone()]])),
            Op::identity(),
            Op::default(),
        ],
        vec![Op::default(), Op::default(), Op::identity()],
    ];
    Ok(Some((
        interval,
        generated(
            t,
            format!("scaled lag-sum kappa={k}"),
            &[Q::zero(), Q::zero(), Q::one()],
        ),
    )))
}

pub(super) fn run_search(
    request: &InterlacingRefinementRequest,
) -> Result<InterlacingRefinementReport, RecurrenceError> {
    validate(request)?;
    let source = match (&request.source, &request.recurrence) {
        (Some(s), _) => Some(s.clone()),
        (_, Some(r)) => Some(companion(r, request.limits.maximum_dimension)?),
        _ => None,
    };
    let (first, states) = if let Some(s) = &source {
        let (_, _, seed) = decode_system(s)?;
        (s.first_index, vec![seed])
    } else {
        let n = request.states[0].len();
        if n == 0
            || n > request.limits.maximum_dimension
            || request
                .states
                .iter()
                .any(|s| s.len() != n || s.iter().any(|p| p.len() > 65))
        {
            return Err(invalid(
                "invalid finite state dimensions or polynomial degrees",
            ));
        }
        (
            request.first_index,
            request
                .states
                .iter()
                .map(|s| s.iter().map(|p| numbers(p)).collect())
                .collect::<Result<Vec<_>, _>>()?,
        )
    };
    if first > usize::MAX - 256 {
        return Err(invalid("first index too large"));
    }
    let n = states[0].len();
    if n > request.limits.maximum_dimension {
        return Err(invalid("source dimension exceeds limits"));
    }
    let weights = if source.is_none() {
        vec![Q::one(); n]
    } else {
        source_weights(request, n)?
    };
    let targets = request
        .target_rows
        .iter()
        .map(|p| numbers(p))
        .collect::<Result<Vec<_>, _>>()?;
    let mut candidates = request.candidates.clone();
    let mut scale = None;
    if request.derive_candidates {
        if let Some(source) = &source {
            if let Some((interval, c)) = scaled_template(source)? {
                if interval.feasible {
                    candidates.push(c);
                }
                scale = Some(interval);
            }
        }
        candidates.push(generated(
            identity(n),
            "companion coordinates".into(),
            &weights,
        ));
        // Refine the newest row as its difference from the preceding row and
        // that preceding row. The sum recovers the target identically.
        if n >= 2 && weights[n - 1].is_one() && weights[..n - 1].iter().all(Zero::is_zero) {
            for reverse in [false, true] {
                let mut t = vec![vec![Op::default(); n]; 2];
                t[0][n - 2] = Op::identity();
                t[1][n - 1] = Op::identity();
                t[1][n - 2] = Op::identity().scale(&-Q::one());
                if reverse {
                    t.reverse();
                }
                candidates.push(generated(
                    t,
                    format!("last-step difference reverse={reverse}"),
                    &[Q::one(), Q::one()],
                ));
            }
        }
        // Preserve a designated target coordinate while adding derivative
        // coordinates; their target weights are zero, so recovery is exact.
        if request.limits.derivative_order > 0 && n < request.limits.maximum_dimension {
            let mut t = identity(n);
            let mut derivative = vec![Op::default(); n];
            derivative[n - 1] = Op::single(0, 1, Biv::one());
            t.insert(0, derivative);
            let mut w = vec![Q::zero()];
            w.extend(weights.clone());
            candidates.push(generated(
                t,
                "derivative plus companion coordinates".into(),
                &w,
            ));
        }
        // Small sparse polynomial shears; rectangular/singular proposals above
        // are retained rather than requiring invertibility.
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                for k in [
                    Q::new(1.into(), 2.into()),
                    Q::one(),
                    Q::from_integer(2.into()),
                ] {
                    let mut t = identity(n);
                    t[i][j] = Op::single(0, 0, Biv::new(vec![vec![Q::zero(), k.clone()]]));
                    candidates.push(generated(
                        t,
                        format!("sparse x shear ({i},{j})={k}"),
                        &weights,
                    ));
                }
            }
        }
    }
    let capped = candidates.len() > request.limits.maximum_candidates;
    candidates.truncate(request.limits.maximum_candidates);
    let mut checks = Vec::new();
    for c in candidates {
        checks.push(check(
            &c,
            source.as_ref(),
            &states,
            &targets,
            &weights,
            first,
            &request.limits,
        )?);
    }
    let selected = checks
        .iter()
        .enumerate()
        .filter(|(_, c)| c.status == "proved")
        .min_by_key(|(_, c)| {
            let m = c.transition.as_ref().unwrap();
            let nonzero = m
                .matrix
                .iter()
                .flatten()
                .flatten()
                .flat_map(|p| p.coeffs.iter().flatten())
                .filter(|s| s.as_str() != "0")
                .count();
            (m.matrix.len(), nonzero)
        })
        .map(|(i, _)| i);
    let status = if selected.is_some() {
        "proved"
    } else if checks.iter().any(|c| c.status == "finite_fit") {
        "finite_fit"
    } else {
        "not_proved"
    };
    Ok(InterlacingRefinementReport{schema:INTERLACING_REFINEMENT_SCHEMA.into(),request:request.clone(),status:status.into(),first_index:first,candidates:checks,selected,candidate_limit_reached:capped,scaled_interval:scale,reference:"Branden, arXiv:1410.6601, Theorem 8.5, zero-aware form with every weakly ordered row/column pair (RealRooted.PreservesInterlacingSeq0.iff_nonneg_and_2x2). A nonnegative scalar collapse is real-rooted or zero. Weyl Euler-lag transitions use a separate Euler/Rolle and common-interleaver argument.".into(),lean_replay:"not_attempted".into()})
}
