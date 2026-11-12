//! Row observations, assumptions and root-domain theorem checks.

use super::*;

type IsolationKey = (RootSubject, usize, usize);

fn property_holds(
    isolation: &RealRootIsolation,
    property: RootProperty,
    interval: Option<&RealInterval>,
) -> bool {
    match (property, interval) {
        (RootProperty::RealRooted, _) => isolation.is_real_rooted(),
        (RootProperty::RealRootedIn, Some(i)) => isolation.real_rooted_in(i),
        (RootProperty::RealRootsIn, Some(i)) => isolation.all_real_roots_in(i),
        _ => false,
    }
}

fn subject_name(subject: RootSubject, component: usize, names: &[String]) -> String {
    match subject {
        RootSubject::Source => "P".into(),
        RootSubject::Transformed => names.get(component).cloned().unwrap_or_else(|| "Q".into()),
    }
}

fn property_text(property: RootProperty, interval: Option<&RootInterval>) -> String {
    match (property, interval) {
        (RootProperty::RealRooted, _) => "real-rooted".into(),
        (RootProperty::RealRootedIn, Some(i)) => format!("real-rooted in {i}"),
        (RootProperty::RealRootsIn, Some(i)) => format!("every real root in {i}"),
        _ => "malformed".into(),
    }
}

/// Image of `interval` (possibly index-dependent, on the indices
/// `first <= n (<= last)`) under a scalar equivalence.  Index-independent
/// intervals map exactly; index-dependent ones map exactly through affine maps
/// (endpoints stay affine in `n`) and through monomial factors that add no
/// root outside them; otherwise (a Möbius map with a pole, whose image
/// endpoints are rational in `n`) the hull over the index domain, an
/// index-independent superset, is mapped and `exact` is `false`.
fn interval_image(
    geometry: &Geometry,
    to_target: bool,
    interval: &RootInterval,
    full_degree: bool,
    first: usize,
    last: Option<usize>,
) -> Option<(RootInterval, bool)> {
    let map = |real: &RealInterval| {
        if to_target {
            geometry.forward(real)
        } else {
            geometry.backward(real, full_degree)
        }
    };
    if interval.is_index_independent() {
        let image = map(&interval.at(0).ok()?)?;
        return Some((RootInterval::from_real(&image), true));
    }
    match geometry {
        Geometry::Monomial { power, multiply } => {
            let adds_zero = *power > 0 && *multiply == to_target;
            if !adds_zero || contains_point_throughout(interval, &Q::zero(), first, last).ok()? {
                return Some((interval.clone(), true));
            }
        }
        Geometry::Mobius { a, b, c, e, .. } if c.is_zero() => {
            // Forward y -> (e y - b)/a, backward y -> (a y + b)/e.
            let (alpha, beta) = if to_target {
                (e / a, -(b / a))
            } else {
                (a / e, b / e)
            };
            return affine_interval_image(interval, &alpha, &beta)
                .ok()
                .map(|image| (image, true));
        }
        _ => {}
    }
    let image = map(&domain_hull(interval, first, last).ok()?)?;
    Some((RootInterval::from_real(&image), false))
}

/// Transfer `property` on `interval` through a scalar equivalence; the flag
/// is `false` when an index-dependent interval was replaced by its hull.
fn transfer(
    geometry: &Geometry,
    to_target: bool,
    property: RootProperty,
    interval: Option<&RootInterval>,
    full_degree: bool,
    domain: (usize, Option<usize>),
) -> Option<(RootProperty, Option<RootInterval>, bool)> {
    if !geometry.is_scalar_equivalence() {
        return None;
    }
    let image = interval
        .and_then(|i| interval_image(geometry, to_target, i, full_degree, domain.0, domain.1));
    match (property, image) {
        (RootProperty::RealRooted, _) => Some((RootProperty::RealRooted, None, true)),
        (p, Some((image, exact))) => Some((p, Some(image), exact)),
        (RootProperty::RealRootedIn, None) => Some((RootProperty::RealRooted, None, true)),
        (RootProperty::RealRootsIn, None) => None,
    }
}

const HULL_NOTE: &str = "; the index-dependent interval is mapped through its hull over the \
                         index domain (an index-independent superset)";

/// Exact isolation of the first `max_observed_rows` rows: the source row and
/// every component of the transformed state at each index.
pub(in crate::recurrence::root_domain) fn observe_rows(
    request: &RootDomainRequest,
    origin: usize,
    rows: &[Vec<Q>],
    states: &[Vec<Vec<Q>>],
    notes: &mut Vec<String>,
) -> (
    BTreeMap<IsolationKey, RealRootIsolation>,
    Vec<RowRootObservation>,
) {
    // Exact isolation of the observed rows.
    let mut isolations: BTreeMap<IsolationKey, RealRootIsolation> = BTreeMap::new();
    let mut observations = Vec::new();
    let mut skipped = Vec::new();
    for (i, row) in rows.iter().enumerate().take(request.max_observed_rows) {
        let index = origin + i;
        let state = &states[i];
        let too_large = std::iter::once(row)
            .chain(state.iter())
            .any(|p| vec_poly::degree(p).is_some_and(|d| d > request.max_observed_degree));
        if too_large {
            skipped.push(index);
            continue;
        }
        let rows = std::iter::once((RootSubject::Source, 0, row)).chain(
            state
                .iter()
                .enumerate()
                .map(|(c, p)| (RootSubject::Transformed, c, p)),
        );
        for (subject, component, polynomial) in rows {
            let isolation = isolate_real_roots_rational_coeffs(polynomial);
            observations.push(RowRootObservation {
                subject,
                component,
                index,
                degree: isolation.degree,
                valuation: isolation.valuation,
                real_rooted: isolation.is_real_rooted(),
                real_root_count: isolation.real_root_count,
                roots: RootEnclosureJson::from_isolation(&isolation),
            });
            isolations.insert((subject, component, index), isolation);
        }
    }
    if rows.len() > request.max_observed_rows {
        notes.push(format!(
            "roots isolated on the first {} of {} rows (max_observed_rows)",
            request.max_observed_rows,
            rows.len()
        ));
    }
    if !skipped.is_empty() {
        notes.push(format!(
            "rows {skipped:?} exceed max_observed_degree = {} and were not isolated",
            request.max_observed_degree
        ));
    }

    (isolations, observations)
}

/// Finite exact observation claims, one per subject and component.
pub(in crate::recurrence::root_domain) fn observation_claims(
    isolations: &BTreeMap<IsolationKey, RealRootIsolation>,
    components: usize,
    names: &[String],
    claims: &mut Vec<RootDomainClaim>,
    notes: &mut Vec<String>,
) {
    let subjects = std::iter::once((RootSubject::Source, 0))
        .chain((0..components).map(|c| (RootSubject::Transformed, c)));
    for (subject, component) in subjects {
        let observed: Vec<(&usize, &RealRootIsolation)> = isolations
            .iter()
            .filter(|((s, c, _), _)| *s == subject && *c == component)
            .map(|((_, _, index), isolation)| (index, isolation))
            .collect();
        let (Some(first), Some(last)) = (observed.first(), observed.last()) else {
            continue;
        };
        let failing: Vec<usize> = observed
            .iter()
            .filter(|(_, iso)| !iso.is_real_rooted())
            .map(|(index, _)| **index)
            .collect();
        let hull = observed
            .iter()
            .filter_map(|(_, iso)| iso.enclosure_hull())
            .reduce(|a, b| a.hull(&b));
        let name = subject_name(subject, component, names);
        let (property, interval) = match (&hull, failing.is_empty()) {
            (Some(h), true) => (RootProperty::RealRootedIn, Some(RootInterval::from_real(h))),
            (None, true) => (RootProperty::RealRooted, None),
            (Some(h), false) => (RootProperty::RealRootsIn, Some(RootInterval::from_real(h))),
            (None, false) => {
                notes.push(format!(
                    "{name}_n has no real roots and is not real-rooted at the observed indices \
                     {failing:?}"
                ));
                continue;
            }
        };
        let mut derivation = format!(
            "exact root isolation of {name}_n at the {} observed indices; the interval is the \
             closed hull of the rational isolating enclosures",
            observed.len()
        );
        if !failing.is_empty() {
            derivation.push_str(&format!("; not real-rooted at indices {failing:?}"));
        }
        claims.push(RootDomainClaim {
            subject,
            component,
            property,
            interval,
            first_index: *first.0,
            last_index: Some(*last.0),
            evidence: RootEvidence::FiniteExactObservation,
            conditional_on: Vec::new(),
            derivation,
        });
    }
}

/// Check supplied assumption `k` exactly on the observed rows, record the
/// check and the supplied claim, and return its first index and label, or
/// `None` when the rows refute it.
pub(in crate::recurrence::root_domain) fn check_assumption(
    k: usize,
    assumption: &RootAssumption,
    isolations: &BTreeMap<IsolationKey, RealRootIsolation>,
    domain_first: usize,
    names: &[String],
    sinks: (
        &mut Vec<AssumptionCheck>,
        &mut Vec<RootDomainClaim>,
        &mut Vec<String>,
    ),
) -> Result<Option<(usize, String)>, RecurrenceError> {
    let (assumption_checks, claims, notes) = sinks;
    let first = assumption.first_index.unwrap_or(domain_first);
    let mut checked = Vec::new();
    let mut refuted_at = None;
    for ((subject, component, index), isolation) in isolations {
        if *subject != assumption.subject
            || *component != assumption.component
            || *index < first
            || assumption.last_index.is_some_and(|last| *index > last)
        {
            continue;
        }
        let interval = assumption
            .interval
            .as_ref()
            .map(|i| i.at(*index))
            .transpose()?;
        checked.push(*index);
        if refuted_at.is_none()
            && !property_holds(isolation, assumption.property, interval.as_ref())
        {
            refuted_at = Some(*index);
        }
    }
    let status = match (checked.is_empty(), refuted_at) {
        (_, Some(_)) => "refuted_on_rows",
        (true, None) => "not_checked",
        (false, None) => "consistent_with_rows",
    };
    assumption_checks.push(AssumptionCheck {
        assumption: k,
        status: status.into(),
        checked_indices: checked,
        refuted_at,
    });
    let name = subject_name(assumption.subject, assumption.component, names);
    let label = format!(
        "assumption {k}: {name}_n {} for n >= {first}{}",
        property_text(assumption.property, assumption.interval.as_ref()),
        assumption
            .last_index
            .map_or(String::new(), |last| format!(" and n <= {last}"))
    );
    claims.push(RootDomainClaim {
        subject: assumption.subject,
        component: assumption.component,
        property: assumption.property,
        interval: assumption.interval.clone(),
        first_index: first,
        last_index: assumption.last_index,
        evidence: RootEvidence::SuppliedAssumption,
        conditional_on: vec![label.clone()],
        derivation: format!(
            "supplied in the request{}",
            assumption
                .note
                .as_ref()
                .map_or(String::new(), |note| format!(" ({note})"))
        ),
    });
    if refuted_at.is_some() {
        notes.push(format!(
            "{label} is refuted on the rows; nothing is transferred"
        ));
        return Ok(None);
    }
    Ok(Some((first, label)))
}

pub(crate) fn analyze_root_domain(
    inputs: &RootDomainInputs<'_>,
) -> Result<RootDomainReport, RecurrenceError> {
    let request = inputs.request;
    request.validate()?;
    let transformation = inputs.transformation;
    let geometry = Geometry::of(transformation)?;
    let record = transformation.record();
    let names: Vec<String> = record.components.iter().map(|c| c.name.clone()).collect();
    let components = transformation.component_count();
    let last_index = transformation.schedule().and_then(|s| s.last_index);
    let domain_first = transformation
        .schedule()
        .map_or(inputs.origin, |s| s.first_index.max(inputs.origin));
    let mut notes = Vec::new();
    let mut claims = Vec::new();

    let (isolations, observations) = observe_rows(
        request,
        inputs.origin,
        inputs.rows,
        inputs.states,
        &mut notes,
    );
    observation_claims(&isolations, components, &names, &mut claims, &mut notes);

    // The map's implication with hypothesis statuses.
    let mut implication = root_implication(transformation)?;
    let base_hypotheses = condition_hypotheses(inputs.conditions);
    let extra = implication
        .hypotheses
        .split_off(record.image_conditions.len());
    implication.hypotheses = base_hypotheses;
    for mut hypothesis in extra {
        if matches!(geometry, Geometry::Symmetric) {
            let mut checked = 0;
            let mut failed = Vec::new();
            let observed = inputs
                .states
                .iter()
                .enumerate()
                .take(request.max_observed_rows)
                .filter(|(_, state)| {
                    state.iter().all(|p| {
                        vec_poly::degree(p).is_none_or(|d| d <= request.max_observed_degree)
                    })
                });
            for (i, state) in observed {
                let a = crate::vec_poly::clear_denominators(&state[0]);
                let xb = crate::vec_poly::clear_denominators(&vec_poly::shift(&state[1]));
                checked += 1;
                if !interlaces_either_order_bigint_coeffs(&a, &xb) {
                    failed.push(inputs.origin + i);
                }
            }
            (hypothesis.status, hypothesis.detail) = if failed.is_empty() {
                (
                    HypothesisStatus::CheckedSuppliedRows,
                    format!("exact on the {checked} observed rows; finite evidence only"),
                )
            } else {
                (
                    HypothesisStatus::FailedOnRows,
                    format!("fails at indices {failed:?}"),
                )
            };
        }
        implication.hypotheses.push(hypothesis);
    }
    let conditional = uncertified(inputs.conditions);
    let all_certified = conditional.is_empty();
    let mut implications = vec![implication];

    // The nonnegative-coefficient window.
    let mut candidates: Vec<(RootInterval, String)> = request
        .candidate_intervals
        .iter()
        .map(|interval| (interval.clone(), "supplied".to_string()))
        .collect();
    if request.default_window {
        if let Some(window) = geometry.forward(&RealInterval::at_most(Q::zero())) {
            let negative: Vec<usize> = inputs
                .rows
                .iter()
                .enumerate()
                .filter(|(_, row)| row.iter().any(Signed::is_negative))
                .map(|(i, _)| inputs.origin + i)
                .collect();
            let window_json = RootInterval::from_real(&window);
            let (status, detail) = if negative.is_empty() {
                (
                    HypothesisStatus::CheckedSuppliedRows,
                    format!(
                        "exact on all {} supplied rows; finite evidence only",
                        inputs.rows.len()
                    ),
                )
            } else {
                (
                    HypothesisStatus::FailedOnRows,
                    format!("negative coefficients at indices {negative:?}"),
                )
            };
            let mut hypotheses = condition_hypotheses(inputs.conditions);
            hypotheses.push(RootHypothesis {
                condition: "P_n has nonnegative coefficients".into(),
                status,
                detail,
            });
            let f_polynomial = matches!(
                &geometry,
                Geometry::Mobius { a, b, c, e, .. }
                    if a.is_one() && b.is_zero() && c.is_one() && e.is_one()
            );
            implications.push(RootImplication {
                theorem: "nonnegative_coefficient_window".into(),
                direction: RootImplicationDirection::SourceToTarget,
                source_property: "P_n nonzero with nonnegative coefficients (every real root \
                                  in (-inf, 0])"
                    .into(),
                target_property: format!("every real root of Q_n in {window_json}"),
                statement: if f_polynomial {
                    "Q = (1 + x)^d P(x/(1 + x)) = sum_j p_j x^j (1 + x)^(d - j) and \
                     (-1)^d Q(-1 - x) = sum_j p_j (1 + x)^j x^(d - j) have nonnegative \
                     coefficients and deg Q <= d, so every real root of Q lies in [-1, 0]; \
                     P is real-rooted iff Q is.  The window does not prove that Q splits."
                        .into()
                } else {
                    format!(
                        "the real roots of P_n lie in (-inf, 0] and the map sends them into \
                         {window_json}; the window does not prove real-rootedness"
                    )
                },
                hypotheses,
                degree_bounds: record.degree_schedule.clone(),
                coefficient_signs: Some("p_(n,j) >= 0 for every j".into()),
                endpoint_conventions: implications[0].endpoint_conventions.clone(),
                zero_cases: ZERO_CASES.to_string(),
                index_domain: record.index_domain.clone(),
                converse: Some(
                    "fails: P = 1 + z + z^2 has nonnegative coefficients, the window holds \
                     vacuously for Q = 1 + 3x + 3x^2, and neither is real-rooted"
                        .into(),
                ),
            });
            if negative.is_empty() && transformation.component_count() == 1 {
                let mut conditional_on = vec![
                    "P_n has nonnegative coefficients for every index (checked on the supplied \
                     rows only)"
                        .to_string(),
                ];
                conditional_on.extend(conditional.iter().cloned());
                claims.push(RootDomainClaim {
                    subject: RootSubject::Transformed,
                    component: 0,
                    property: RootProperty::RealRootsIn,
                    interval: Some(window_json.clone()),
                    first_index: domain_first,
                    last_index,
                    evidence: RootEvidence::ProposedInterval,
                    conditional_on,
                    derivation: "nonnegative_coefficient_window; a candidate window for \
                                 proofs, not a proof of real-rootedness"
                        .into(),
                });
                candidates.push((window_json.clone(), "default_window".into()));
            }
        }
    }

    // Supplied assumptions.
    let mut assumption_checks = Vec::new();
    for (k, assumption) in request.assumptions.iter().enumerate() {
        if assumption.subject == RootSubject::Transformed && assumption.component >= components {
            return Err(RecurrenceError::InvalidInput(format!(
                "assumption {k}: component {} does not exist",
                assumption.component
            )));
        }
        let Some((first, label)) = check_assumption(
            k,
            assumption,
            &isolations,
            domain_first,
            &names,
            (&mut assumption_checks, &mut claims, &mut notes),
        )?
        else {
            continue;
        };
        let to_target = assumption.subject == RootSubject::Source;
        if let Some((property, interval, exact)) = transfer(
            &geometry,
            to_target,
            assumption.property,
            assumption.interval.as_ref(),
            false,
            (first, assumption.last_index),
        ) {
            let mut conditional_on = vec![label];
            conditional_on.extend(conditional.iter().cloned());
            claims.push(RootDomainClaim {
                subject: if to_target {
                    RootSubject::Transformed
                } else {
                    RootSubject::Source
                },
                component: 0,
                property,
                interval,
                first_index: first,
                last_index: assumption.last_index,
                evidence: RootEvidence::SuppliedAssumption,
                conditional_on,
                derivation: format!(
                    "transferred from assumption {k} by {}{}",
                    geometry.theorem(),
                    if exact { "" } else { HULL_NOTE }
                ),
            });
        } else if matches!(geometry, Geometry::Gamma) {
            let quarter = RealInterval::at_most(Q::new(BigInt::one(), BigInt::from(4)));
            let transferred = match (to_target, assumption.property) {
                (true, RootProperty::RealRooted | RootProperty::RealRootedIn) => Some((
                    RootSubject::Transformed,
                    RootProperty::RealRootedIn,
                    Some(RootInterval::from_real(&quarter)),
                )),
                (false, RootProperty::RealRootedIn)
                    if assumption.interval.as_ref().is_some_and(|i| {
                        domain_hull(i, first, assumption.last_index)
                            .is_ok_and(|hull| hull.is_subset_of(&quarter))
                    }) =>
                {
                    Some((RootSubject::Source, RootProperty::RealRooted, None))
                }
                _ => None,
            };
            if let Some((subject, property, interval)) = transferred {
                let mut conditional_on = vec![label];
                conditional_on.extend(conditional.iter().cloned());
                claims.push(RootDomainClaim {
                    subject,
                    component: 0,
                    property,
                    interval,
                    first_index: first,
                    last_index: assumption.last_index,
                    evidence: RootEvidence::SuppliedAssumption,
                    conditional_on,
                    derivation: format!("transferred from assumption {k} by gamma_quarter_bound"),
                });
            }
        }
    }

    // The recurrence theorem match of the transported rule.
    let recurrence_theorem = inputs.transported.as_ref().map(|transported| {
        let context = RolleContext {
            geometry: &geometry,
            supplied: &candidates,
            all_certified,
            conditional: &conditional,
        };
        rolle_induction(inputs, transported, &context, &mut claims, &mut notes)
    });

    Ok(RootDomainReport {
        conventions: conventions(),
        implications,
        claims,
        recurrence_theorem,
        assumption_checks,
        observations,
        notes,
    })
}

/// What the analysis established before the recurrence theorem is matched.
struct RolleContext<'a> {
    geometry: &'a Geometry,
    /// Candidate intervals supplied with the request.
    supplied: &'a [(RootInterval, String)],
    /// Whether every map condition is certified for all indices.
    all_certified: bool,
    /// Labels of the assumptions the claims depend on.
    conditional: &'a [String],
}

fn rolle_induction(
    inputs: &RootDomainInputs<'_>,
    transported: &TransportedInput<'_>,
    context: &RolleContext<'_>,
    claims: &mut Vec<RootDomainClaim>,
    notes: &mut Vec<String>,
) -> RecurrenceRootTheorem {
    let RolleContext {
        geometry,
        supplied,
        all_certified,
        conditional,
    } = *context;
    let rule = match transported.rule {
        TransportedRecurrence::Scalar(rule) => rule,
        _ => return unmatched_theorem("only scalar transported rules are matched"),
    };
    let name = inputs
        .transformation
        .record()
        .components
        .into_iter()
        .map(|c| c.name)
        .next()
        .unwrap_or_else(|| "Q".into());
    let first = transported.first_index;
    let valid_from = transported.valid_from;
    let last_index = inputs.transformation.schedule().and_then(|s| s.last_index);
    let seeds: Vec<(usize, Vec<Q>)> = (first..valid_from)
        .filter_map(|k| {
            k.checked_sub(inputs.origin)
                .and_then(|i| inputs.states.get(i))
                .map(|state| (k, state[0].clone()))
        })
        .collect();
    let outcome = rolle_core(
        &RolleSetup {
            rule,
            name: name.clone(),
            first,
            valid_from,
            seeds,
            candidates: supplied.to_vec(),
            limit: inputs.request.max_candidate_intervals,
            schedule: inputs.transformation.schedule(),
        },
        notes,
    );
    let Some(interval) = outcome.proven.clone() else {
        return outcome.theorem;
    };
    claims.push(RootDomainClaim {
        subject: RootSubject::Transformed,
        component: 0,
        property: RootProperty::RealRootedIn,
        interval: Some(interval.clone()),
        first_index: first,
        last_index,
        evidence: RootEvidence::CompletedRecurrenceInduction,
        conditional_on: conditional.to_vec(),
        derivation: format!(
            "rolle_derivative_of_product on the transported rule (an exact consequence of the \
             defining source rule) with g and the seeds checked exactly on {interval}; the rows \
             it generates are the transformed rows"
        ),
    });
    if !all_certified {
        notes.push(
            "the induction is not transferred to P_n: some map conditions are only checked on \
             the supplied rows"
                .into(),
        );
        return outcome.theorem;
    }
    if matches!(geometry, Geometry::Gamma) {
        gamma_transfer(&interval, (first, last_index), &name, claims, notes);
        return outcome.theorem;
    }
    match transfer(
        geometry,
        false,
        RootProperty::RealRootedIn,
        Some(&interval),
        outcome.full_degree,
        (first, last_index),
    ) {
        Some((property, source_interval, exact)) => claims.push(RootDomainClaim {
            subject: RootSubject::Source,
            component: 0,
            property,
            interval: source_interval,
            first_index: first,
            last_index,
            evidence: RootEvidence::CertifiedFamilyInvariant,
            conditional_on: Vec::new(),
            derivation: format!(
                "the completed induction for {name}_n and {} (target to source; every map \
                 condition certified for all indices){}{}",
                geometry.theorem(),
                if outcome.full_degree {
                    "; deg Q_n = d_n by the degree formula, so no root of P_n is added at the \
                     image of infinity"
                } else {
                    ""
                },
                if exact { "" } else { HULL_NOTE }
            ),
        }),
        None => notes.push(format!(
            "the induction is not transferred to P_n: {} gives no scalar equivalence",
            geometry.theorem()
        )),
    }
    outcome.theorem
}

/// Transfer a completed induction for the gamma rows `G_n` to `P_n`: only
/// when every interval `I_n` lies in `(-inf, 1/4]` (`gamma_quarter_bound`);
/// real-rootedness of `G_n` alone proves nothing about `P_n`.
fn gamma_transfer(
    interval: &RootInterval,
    domain: (usize, Option<usize>),
    name: &str,
    claims: &mut Vec<RootDomainClaim>,
    notes: &mut Vec<String>,
) {
    let quarter = RealInterval::at_most(Q::new(BigInt::one(), BigInt::from(4)));
    let hull = match domain_hull(interval, domain.0, domain.1) {
        Ok(hull) if hull.is_subset_of(&quarter) => hull,
        _ => {
            notes.push(format!(
                "the induction proves {name}_n real-rooted in {interval}, which is not contained \
                 in (-inf, 1/4] for every index: no conclusion for P_n (gamma_quarter_bound needs \
                 every root of the gamma polynomial at most 1/4; real-rootedness alone is \
                 insufficient)"
            ));
            return;
        }
    };
    let nonpositive = hull.is_subset_of(&RealInterval::at_most(Q::zero()));
    claims.push(RootDomainClaim {
        subject: RootSubject::Source,
        component: 0,
        property: if nonpositive {
            RootProperty::RealRootedIn
        } else {
            RootProperty::RealRooted
        },
        interval: nonpositive.then(|| RootInterval::from_real(&RealInterval::at_most(Q::zero()))),
        first_index: domain.0,
        last_index: domain.1,
        evidence: RootEvidence::CertifiedFamilyInvariant,
        conditional_on: Vec::new(),
        derivation: format!(
            "the completed induction puts every root of {name}_n in {interval}, contained in \
             (-inf, 1/4] for every index (hull {hull}); gamma_quarter_bound with the symmetry \
             certified for all indices gives real-rooted P_n{}",
            if nonpositive {
                "; roots of the gamma polynomial at most 0 put every root of P_n in (-inf, 0]"
            } else {
                ""
            }
        ),
    });
}
