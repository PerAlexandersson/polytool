//! Theorem and root-domain routes for parent proofs.

use super::*;

/// Route 2, with the interlacing of the first successful match (module
/// documentation, "Interlacing").
pub(super) fn theorem_route(
    record: &ParentRecurrenceRecord<'_>,
    rows: &[Vec<Q>],
) -> (Vec<VeroneseParentProofAttempt>, Option<LiftInterlacing>) {
    let matches = match_recurrence_theorems(record.rule, record.first_index, rows);
    let mut attempts = Vec::new();
    let mut interlacing = None;
    for theorem_match in matches.theorem_matches {
        let boundary_ok = matches!(
            theorem_match.boundary.status,
            TheoremBoundaryStatus::Canonical | TheoremBoundaryStatus::CanonicalAfterRelabel
        );
        let scope_ok = theorem_match
            .variable_predecessor
            .as_ref()
            .is_none_or(|certificate| certificate.scope == CertificateScope::AllIndices);
        let success = boundary_ok && scope_ok;
        if success && interlacing.is_none() {
            interlacing = paper_theorem_interlacing(&theorem_match, record.first_index);
        }
        let detail = if success {
            format!(
                "{} ({}): {}; the identifying seeds equal the canonical rows ({:?}), so every \
                 parent row is a canonical row",
                theorem_match.family,
                theorem_match.proof_status,
                theorem_match.conclusion,
                theorem_match.boundary.status
            )
        } else if !boundary_ok {
            format!(
                "the parent rows are not the canonical boundary rows ({:?})",
                theorem_match.boundary.status
            )
        } else {
            "the theorem certificate covers finitely many indices only".into()
        };
        attempts.push(VeroneseParentProofAttempt {
            route: VeroneseParentProofRoute::PaperTheorem,
            variant: theorem_match.theorem.clone(),
            success,
            detail,
        });
    }
    if attempts.is_empty() {
        let reasons: Vec<String> = matches
            .rejected_theorems
            .iter()
            .map(|r| format!("{}: {}", r.theorem, r.reason))
            .collect();
        attempts.push(VeroneseParentProofAttempt {
            route: VeroneseParentProofRoute::PaperTheorem,
            variant: "match_recurrence_theorems".into(),
            success: false,
            detail: format!("no theorem matched ({})", reasons.join("; ")),
        });
    }
    (attempts, interlacing)
}

/// The interlacing `P_(m-1) << P_m`, `m >= first_index + 1`, of a matched
/// window theorem whose window contains the step `1` (every canonical row
/// `P_n` interlaces `P_(n+s)` from the left for `0 <= s <= max_step`).
fn paper_theorem_interlacing(
    theorem_match: &RecurrenceTheoremMatch,
    first_index: usize,
) -> Option<LiftInterlacing> {
    if theorem_match.window.max_step.is_some_and(|step| step == 0) {
        return None;
    }
    let mut premises: Vec<String> = theorem_match
        .hypotheses
        .iter()
        .map(|h| format!("{} ({})", h.statement, h.check))
        .collect();
    premises.push(format!(
        "window: {} ({})",
        theorem_match.window.statement, theorem_match.window.inequality
    ));
    premises.push(format!(
        "the identifying parent seeds equal the canonical rows ({:?}), so every parent row is a \
         canonical row and consecutive parent rows are consecutive canonical rows",
        theorem_match.boundary.status
    ));
    Some(LiftInterlacing {
        kind: RowInterlacing::Weak,
        orientation: InterlacingOrientation::PreviousInterlacesNext,
        from_index: first_index.checked_add(1)?,
        justification: LiftInterlacingJustification {
            theorem: format!(
                "{} ({}; {}): {}",
                theorem_match.theorem,
                theorem_match.proof_status,
                theorem_match.reference,
                theorem_match.conclusion
            ),
            premises,
        },
    })
}

/// `d_n = d_(n0) + step (n - n0)` fitted to the degrees of the first rows,
/// when every observed row is nonzero and the degrees are affine.
fn fitted_schedule(first_index: usize, rows: &[Vec<Q>]) -> Option<DegreeSchedule> {
    let degrees: Vec<i64> = rows
        .iter()
        .map(|row| vec_poly::degree(row).and_then(|d| i64::try_from(d).ok()))
        .collect::<Option<_>>()?;
    let first = *degrees.first()?;
    let step = degrees.get(1).map_or(0, |second| second - first);
    let affine = degrees
        .iter()
        .enumerate()
        .all(|(i, d)| i64::try_from(i).is_ok_and(|i| *d == first + step * i));
    affine.then(|| DegreeSchedule::new(first_index, first as usize, step))
}

fn accepted_source_claim(claim: &RootDomainClaim, first_index: usize) -> bool {
    let evidence_ok = matches!(
        claim.evidence,
        RootEvidence::CompletedRecurrenceInduction | RootEvidence::CertifiedFamilyInvariant
    );
    let upper_ok = claim.interval.as_ref().is_some_and(|interval| {
        interval.is_index_independent()
            && interval.upper.as_ref().is_some_and(|upper| {
                parse_rational_coeff(&upper.value).is_ok_and(|value| value <= Q::zero())
            })
    });
    claim.subject == RootSubject::Source
        && claim.component == 0
        && claim.property == RootProperty::RealRootedIn
        && upper_ok
        && evidence_ok
        && claim.conditional_on.is_empty()
        && claim.last_index.is_none()
        && claim.first_index <= first_index
}

pub(super) fn root_domain_route(
    record: &ParentRecurrenceRecord<'_>,
    rows: &[Vec<Q>],
) -> Vec<VeroneseParentProofAttempt> {
    let mut attempts = Vec::new();
    let seed_count = record.seeds.len();
    for delay in 0..=MAXIMUM_ROOT_DOMAIN_START_DELAY {
        // The moved record needs its seeds among the exactly checked rows.
        let Some(seeds) = rows.get(delay..delay + seed_count) else {
            break;
        };
        let moved = ParentRecurrenceRecord {
            rule: record.rule,
            seeds,
            first_index: record.first_index + delay,
        };
        let (moved_attempts, success) = root_domain_menu(&moved, &rows[delay..], delay);
        attempts.extend(moved_attempts);
        if success {
            break;
        }
    }
    attempts
}

/// The root-domain menu for `record` (moved `delay` rows past the parent's
/// own seeds; `rows` start at the record's first index).
fn root_domain_menu(
    record: &ParentRecurrenceRecord<'_>,
    rows: &[Vec<Q>],
    delay: usize,
) -> (Vec<VeroneseParentProofAttempt>, bool) {
    let source = RecurrenceJson::from_recurrence_rational(
        record.rule,
        record.first_index,
        record.seeds,
        None,
    );
    let schedule = fitted_schedule(record.first_index, rows);
    let mut menu = vec![(
        "identity (affine 1, 0)".to_string(),
        Some(RecurrenceTransformation::Affine {
            a: "1".into(),
            b: "0".into(),
        }),
    )];
    for (name, make) in [
        (
            "f_polynomial",
            (|s: DegreeSchedule| RecurrenceTransformation::FPolynomial { schedule: s })
                as fn(DegreeSchedule) -> RecurrenceTransformation,
        ),
        ("reversal", |s| RecurrenceTransformation::Reversal {
            schedule: s,
        }),
        ("mobius (1, 0, -1, 1)", |s| {
            RecurrenceTransformation::Mobius {
                a: "1".into(),
                b: "0".into(),
                c: "-1".into(),
                e: "1".into(),
                schedule: s,
            }
        }),
    ] {
        menu.push((name.to_string(), schedule.clone().map(make)));
    }
    let variant = |name: String| {
        if delay == 0 {
            name
        } else {
            format!("{name}, induction start moved {delay} row(s)")
        }
    };
    let moved_note = if delay == 0 {
        String::new()
    } else {
        format!(
            " (source moved {delay} row(s) past the parent seeds; the earlier parent rows are \
             among the rows checked exactly to be zero or real-rooted with coefficients of one \
             sign)"
        )
    };
    let mut attempts = Vec::new();
    for (name, transformation) in menu {
        let Some(transformation) = transformation else {
            attempts.push(VeroneseParentProofAttempt {
                route: VeroneseParentProofRoute::RootDomain,
                variant: variant(name),
                success: false,
                detail: "no affine degree schedule fits the first parent rows".into(),
            });
            continue;
        };
        let mut request = RecurrenceTransformRequest::from_source(
            transformation,
            source.clone(),
            ROOT_DOMAIN_ROWS,
        );
        request.root_domain = Some(RootDomainRequest::default());
        let attempt = match transform_recurrence(&request) {
            Err(error) => VeroneseParentProofAttempt {
                route: VeroneseParentProofRoute::RootDomain,
                variant: variant(name),
                success: false,
                detail: format!("the transform failed: {error}"),
            },
            Ok(report) => {
                let claim = report.root_domain.as_ref().and_then(|root| {
                    root.claims
                        .iter()
                        .find(|claim| accepted_source_claim(claim, record.first_index))
                });
                match claim {
                    Some(claim) => VeroneseParentProofAttempt {
                        route: VeroneseParentProofRoute::RootDomain,
                        variant: variant(name),
                        success: true,
                        detail: format!(
                            "source claim: every P_n, n >= {}, is real-rooted in {} ({:?}): \
                             {}{moved_note}",
                            claim.first_index,
                            claim
                                .interval
                                .as_ref()
                                .map_or_else(String::new, ToString::to_string),
                            claim.evidence,
                            claim.derivation
                        ),
                    },
                    None => VeroneseParentProofAttempt {
                        route: VeroneseParentProofRoute::RootDomain,
                        variant: variant(name),
                        success: false,
                        detail: "no unconditional all-index source claim with roots in \
                                 (-inf, 0]"
                            .into(),
                    },
                }
            }
        };
        let success = attempt.success;
        attempts.push(attempt);
        if success {
            return (attempts, true);
        }
    }
    (attempts, false)
}

/// A theorem route of the registry ([`TheoremRoute::policy`]): the theorem
/// with the requested region of its [`ParentRegion`] policy, proved for the
/// parent only when the certified region lies in `(-inf, 0]` (every parent
/// row then has coefficients of one sign).  The last value is whether the
/// theorem applies (status other than `not_applicable`).
pub(super) fn theorem_attempt(
    record: &ParentRecurrenceRecord<'_>,
    route: &TheoremRoute,
) -> Result<(VeroneseParentProofAttempt, Option<LiftInterlacing>, bool), RecurrenceError> {
    let check = |options: &RealRootednessTheoremOptions| {
        check_real_rootedness_theorem(
            route.theorem,
            record.rule,
            record.first_index,
            record.seeds,
            options,
        )
    };
    let ray = RealRootednessTheoremOptions {
        root_region: Some(RealInterval {
            lower: None,
            lower_closed: false,
            upper: Some(Q::zero()),
            upper_closed: false,
        }),
        ..RealRootednessTheoremOptions::default()
    };
    let unrestricted = RealRootednessTheoremOptions::default();
    let report = match route.policy.parent.region {
        ParentRegion::OpenNonpositiveRay => check(&ray)?,
        ParentRegion::Unrestricted => check(&unrestricted)?,
        ParentRegion::OpenNonpositiveRayThenUnrestricted => {
            let requested = check(&ray)?;
            if requested.is_proved()
                || requested.status == RealRootednessTheoremStatus::NotApplicable
            {
                requested
            } else {
                check(&unrestricted)?
            }
        }
    };
    let route_name = VeroneseParentProofRoute::of_theorem(route.theorem);
    let variant = route_name.as_str();
    // A proof whose certified region is not contained in (-inf, 0] (the
    // Ma--Wang interval, the Liu--Wang three-term or iterated-preserver
    // region, ...) gives real-rooted rows but not one-signed coefficients.
    let outside_ray = report.is_proved()
        && !report
            .conclusion
            .as_ref()
            .and_then(|c| c.root_region.as_ref())
            .and_then(|region| region.at(0).ok())
            .and_then(|region| region.upper)
            .is_some_and(|upper| !upper.is_positive());
    let detail = match (&report.status, &report.conclusion) {
        (RealRootednessTheoremStatus::Proved, Some(conclusion)) if outside_ray => format!(
            "{}: {}; the {} not contained in (-inf, 0], so the rows are real-rooted but this \
             route does not give coefficients of one sign",
            report.normal_form,
            conclusion.statement,
            if conclusion.root_region.is_some() {
                "certified region is"
            } else {
                "theorem certifies no root region, so the roots are"
            }
        ),
        (RealRootednessTheoremStatus::Proved, Some(conclusion)) if report.start_delay > 0 => {
            format!(
                "induction start moved {} row(s) past the seeds: {}: {}",
                report.start_delay, report.normal_form, conclusion.statement
            )
        }
        (RealRootednessTheoremStatus::Proved, Some(conclusion)) => {
            format!("{}: {}", report.normal_form, conclusion.statement)
        }
        (RealRootednessTheoremStatus::NotApplicable, _) => {
            format!("not applicable: {}", report.normal_form)
        }
        (status, _) => format!(
            "{status:?}: {}; open hypotheses: {}",
            report.normal_form,
            report.open_hypotheses().join("; ")
        ),
    };
    let success = report.is_proved() && !outside_ray;
    let interlacing = success
        .then(|| theorem_report_interlacing(&report, variant))
        .flatten();
    let applicable = report.status != RealRootednessTheoremStatus::NotApplicable;
    Ok((
        VeroneseParentProofAttempt {
            route: route_name,
            variant: variant.into(),
            success,
            detail,
        },
        interlacing,
        applicable,
    ))
}

/// The interlacing conclusion of a proved theorem report: `P_(m-1)`
/// interlaces `P_m` from the left (`P_m` carries the largest zero, the
/// orientation of [`RowInterlacing`]) for every `m >= interlacing_from`.
pub(in crate::recurrence::lift) fn theorem_report_interlacing(
    report: &RealRootednessTheoremReport,
    variant: &str,
) -> Option<LiftInterlacing> {
    let conclusion = report.conclusion.as_ref().filter(|_| report.is_proved())?;
    let (kind, from_index) = (conclusion.interlacing?, conclusion.interlacing_from?);
    let mut premises = vec![format!(
        "every hypothesis of the theorem is verified exactly for every index: {}",
        report.normal_form
    )];
    premises.push(conclusion.statement.clone());
    Some(LiftInterlacing {
        kind,
        orientation: InterlacingOrientation::PreviousInterlacesNext,
        from_index,
        justification: LiftInterlacingJustification {
            theorem: format!("{variant}: {}", report.reference),
            premises,
        },
    })
}

/// Whether two records define the same parent: same first index, the same
/// relation term by term, the same number of identifying seeds and equal
/// seeds (the supplied rows of `other` must agree with its rule).
pub(super) fn same_record(
    record: &ParentRecurrenceRecord<'_>,
    other: &RecurrenceJson,
) -> Result<bool, RecurrenceError> {
    let (rule, first_index, rows) = other.to_recurrence_parts()?;
    if first_index != record.first_index || rows.len() < rule.max_offset() {
        return Ok(false);
    }
    let (left, left_shift) = ImplicitRelation::from_recurrence(record.rule)?;
    let (right, right_shift) = ImplicitRelation::from_recurrence(&rule)?;
    if left_shift != right_shift || left.op() != right.op() {
        return Ok(false);
    }
    let count = rule.generation_initial_count(first_index, rows.len());
    if count != record.seeds.len() {
        return Ok(false);
    }
    let generated = match rule.generate_rows_rational(&rows[..count], first_index, rows.len()) {
        Ok(generated) => generated,
        Err(_) => return Ok(false),
    };
    let equal = |a: &[Q], b: &[Q]| vec_poly::trim(a) == vec_poly::trim(b);
    Ok(rows.iter().zip(&generated).all(|(a, b)| equal(a, b))
        && rows[..count]
            .iter()
            .zip(record.seeds)
            .all(|(a, b)| equal(a, b)))
}
