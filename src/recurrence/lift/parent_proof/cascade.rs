//! Ordered parent-proof cascade and public entry points.

use super::*;

/// Run the parent provers (module documentation): try to prove that every
/// row `P_n`, `n >= record.first_index`, of the parent is zero or
/// real-rooted with coefficients of one sign, for every `n`.
///
/// `needed_from` is the first parent index the caller uses (at least
/// `first_index`); `attestations` are the caller-supplied statements of
/// route 9 (a verdict resting on one is conditional, never proved);
/// `try_root_domain` enables route 3.  A refuted or unproved statement is
/// an `Ok` proof record with that status.
///
/// # Errors
///
/// [`RecurrenceError::InvalidInput`] when the parent rule cannot be
/// replayed from its seeds, and errors of the theorem checks on malformed
/// records.
pub fn prove_parent_real_rootedness(
    record: &ParentRecurrenceRecord<'_>,
    needed_from: usize,
    attestations: &[VeroneseParentAttestation],
    try_root_domain: bool,
) -> Result<VeroneseParentProof, RecurrenceError> {
    prove_parent_with_routes(record, needed_from, attestations, try_root_domain, None)
}

/// Internal closure policy: the structural mode never enters either attestation route.
pub(in crate::recurrence) fn prove_parent_with_policy(
    record: &ParentRecurrenceRecord<'_>,
    needed_from: usize,
    attestations: &[VeroneseParentAttestation],
    try_root_domain: bool,
    allow_attestations: bool,
) -> Result<VeroneseParentProof, RecurrenceError> {
    if allow_attestations {
        return prove_parent_real_rootedness(record, needed_from, attestations, try_root_domain);
    }
    let mut routes: Vec<_> = RealRootednessTheorem::ALL
        .iter()
        .copied()
        .map(VeroneseParentProofRoute::of_theorem)
        .collect();
    routes.extend([
        VeroneseParentProofRoute::PeriodicProduct,
        VeroneseParentProofRoute::PaperTheorem,
    ]);
    if try_root_domain {
        routes.push(VeroneseParentProofRoute::RootDomain);
    }
    prove_parent_by_routes(record, needed_from, &routes)
}

/// Replay of the routes `routes` alone (module documentation, "Replaying
/// recorded proofs"): the provers with every other route skipped (no
/// attempt recorded for it), the root-domain route enabled when it is
/// among `routes`, and no attestations.  The result is proved only by one
/// of `routes`; it is a correct proof whenever it is proved, like every
/// result of the provers.
pub(in crate::recurrence) fn prove_parent_by_routes(
    record: &ParentRecurrenceRecord<'_>,
    needed_from: usize,
    routes: &[VeroneseParentProofRoute],
) -> Result<VeroneseParentProof, RecurrenceError> {
    let try_root_domain = routes.contains(&VeroneseParentProofRoute::RootDomain);
    prove_parent_with_routes(record, needed_from, &[], try_root_domain, Some(routes))
}

/// Run the theorem routes of `stage` in the registry order
/// ([`TheoremRoute::policy`]): each route runs while the stage allows it
/// (unproved, or for [`ParentStage::ReplacesEarly`] proved only by an early
/// route), its attempt is recorded (unless it does not apply and the route
/// records only applicable attempts), and its proof sets `proved` and the
/// interlacing.
pub(super) fn run_theorem_routes(
    record: &ParentRecurrenceRecord<'_>,
    stage: ParentStage,
    allowed: &dyn Fn(VeroneseParentProofRoute) -> bool,
    proved: &mut Option<VeroneseParentProofRoute>,
    interlacing: &mut Option<LiftInterlacing>,
    attempts: &mut Vec<VeroneseParentProofAttempt>,
) -> Result<(), RecurrenceError> {
    for route in ROUTES
        .iter()
        .filter(|route| route.policy.parent.stage == stage)
    {
        let proved_early = proved
            .and_then(VeroneseParentProofRoute::theorem)
            .is_some_and(|theorem| {
                ROUTES.iter().any(|early| {
                    early.theorem == theorem && early.policy.parent.stage == ParentStage::Early
                })
            });
        let runs = match stage {
            ParentStage::Early | ParentStage::WhileUnproved => proved.is_none(),
            ParentStage::ReplacesEarly => proved.is_none() || proved_early,
        };
        let route_name = VeroneseParentProofRoute::of_theorem(route.theorem);
        if !runs || !allowed(route_name) {
            continue;
        }
        let (attempt, theorem_interlacing, applicable) = theorem_attempt(record, route)?;
        if route.policy.parent.record_only_if_applicable && !applicable {
            continue;
        }
        let success = attempt.success;
        attempts.push(attempt);
        if success {
            *proved = Some(route_name);
            *interlacing = theorem_interlacing;
        }
    }
    Ok(())
}

/// The parent provers, restricted to `only` when given (attestation and
/// catalog-flag routes are then skipped too).
fn prove_parent_with_routes(
    record: &ParentRecurrenceRecord<'_>,
    needed_from: usize,
    attestations: &[VeroneseParentAttestation],
    try_root_domain: bool,
    only: Option<&[VeroneseParentProofRoute]>,
) -> Result<VeroneseParentProof, RecurrenceError> {
    let allowed = |route: VeroneseParentProofRoute| only.is_none_or(|r| r.contains(&route));
    let rows = record
        .rule
        .generate_rows_rational(
            record.seeds,
            record.first_index,
            PARENT_ROWS_CHECKED.max(record.seeds.len()),
        )
        .map_err(|error| {
            RecurrenceError::InvalidInput(format!(
                "the parent recurrence cannot be replayed: {error}"
            ))
        })?;
    // Routes 1 to 8 and the one-signed test of route 10 run on the rule with
    // a left-hand factor free of x when the normalization applies; a rule
    // that does not normalize is kept, and the theorem routes name the reason.
    let (normalized_rule, normalization) = rule_with_x_free_left_hand_factor(record.rule)
        .unwrap_or_else(|_| (record.rule.clone(), None));
    let normalized = ParentRecurrenceRecord {
        rule: &normalized_rule,
        seeds: record.seeds,
        first_index: record.first_index,
    };
    let mut proof = VeroneseParentProof {
        status: VeroneseParentProofStatus::NotProved,
        route: None,
        statement: String::new(),
        conditional_on: Vec::new(),
        attempts: Vec::new(),
        failing_row_index: first_failing_row(&rows, record.first_index),
        rows_checked: rows.len(),
        interlacing: None,
        left_hand_factor_normalization: normalization,
    };
    let conclusion = format!(
        "every parent row P_n, n >= {}, is zero or real-rooted with coefficients of one sign",
        record.first_index
    );
    if let Some(index) = proof.failing_row_index {
        proof.status = VeroneseParentProofStatus::Refuted;
        proof.statement = format!(
            "parent row P_{index} is nonzero and not real-rooted with coefficients of one sign; \
             this parent gives no proof"
        );
        return Ok(proof);
    }

    // Route 1.
    let mut proved = None;
    if allowed(VeroneseParentProofRoute::PeriodicProduct) {
        let periodic = periodic_product_parent(normalized.rule, normalized.seeds);
        proof.attempts.push(VeroneseParentProofAttempt {
            route: VeroneseParentProofRoute::PeriodicProduct,
            variant: "P_(n+k) = c(x) P_n".into(),
            success: periodic.is_ok(),
            detail: match &periodic {
                Ok(()) => "c and the seeds are nonzero, real-rooted and of one sign, so \
                           P_(n0 + k m + i) = c^m P_(n0 + i) is too"
                    .into(),
                Err(reason) => reason.clone(),
            },
        });
        proved = periodic
            .is_ok()
            .then_some(VeroneseParentProofRoute::PeriodicProduct);
    }
    let mut interlacing = None;

    // The early theorem routes (product induction): real-rootedness only;
    // a later theorem route may replace the proof.
    run_theorem_routes(
        &normalized,
        ParentStage::Early,
        &allowed,
        &mut proved,
        &mut interlacing,
        &mut proof.attempts,
    )?;

    // Route 2.
    if proved.is_none() && allowed(VeroneseParentProofRoute::PaperTheorem) {
        let (attempts, window_interlacing) = theorem_route(&normalized, &rows);
        if attempts.iter().any(|a| a.success) {
            proved = Some(VeroneseParentProofRoute::PaperTheorem);
            interlacing = window_interlacing;
        }
        proof.attempts.extend(attempts);
    }

    // Route 3.
    if proved.is_none() && try_root_domain && allowed(VeroneseParentProofRoute::RootDomain) {
        let attempts = root_domain_route(&normalized, &rows);
        if attempts.iter().any(|a| a.success) {
            proved = Some(VeroneseParentProofRoute::RootDomain);
        }
        proof.attempts.extend(attempts);
    }

    // Routes 4 to 8b: product induction already settles real-rootedness,
    // but these routes may supply the stronger interlacing conclusion.
    // Then routes 8c to 8f, for a parent left open.
    for stage in [ParentStage::ReplacesEarly, ParentStage::WhileUnproved] {
        run_theorem_routes(
            &normalized,
            stage,
            &allowed,
            &mut proved,
            &mut interlacing,
            &mut proof.attempts,
        )?;
    }

    if let Some(route) = proved {
        proof.status = VeroneseParentProofStatus::Proved;
        proof.route = Some(route);
        proof.statement = conclusion;
        proof.interlacing = interlacing;
        return Ok(proof);
    }

    if only.is_some() {
        proof.statement = format!("open: the replayed routes do not prove that {conclusion}");
        return Ok(proof);
    }
    // Route 9.
    for attestation in attestations {
        let matches = same_record(record, &attestation.recurrence)?
            || (proof.left_hand_factor_normalization.is_some()
                && same_record(&normalized, &attestation.recurrence)?);
        let covers = attestation.first_index <= needed_from;
        let success = matches && covers;
        proof.attempts.push(VeroneseParentProofAttempt {
            route: VeroneseParentProofRoute::Attestation,
            variant: attestation.label.clone(),
            success,
            detail: if success {
                format!("{} (source: {})", attestation.statement, attestation.source)
            } else if !matches {
                "the attested record differs from the parent".into()
            } else {
                format!(
                    "the attestation starts at index {}, after the first used parent index \
                     {needed_from}",
                    attestation.first_index
                )
            },
        });
        if success {
            proof.status = VeroneseParentProofStatus::ConditionalOnAttestation;
            proof.route = Some(VeroneseParentProofRoute::Attestation);
            proof.statement = format!("{conclusion} (attested, not proved here)");
            proof.conditional_on = vec![attestation.label.clone()];
            return Ok(proof);
        }
    }
    // Route 10.
    let candidates = catalog_flag_attestations(record);
    let mut differing = Vec::new();
    for attestation in candidates {
        let matches = same_record(record, &attestation.recurrence)?
            || (proof.left_hand_factor_normalization.is_some()
                && same_record(&normalized, &attestation.recurrence)?);
        if !matches {
            differing.push(attestation.label);
            continue;
        }
        let covers = attestation.first_index <= needed_from;
        let one_signed = covers.then(|| {
            one_signed_rows_for_all_indices(normalized.rule, record.first_index, record.seeds)
        });
        let success = matches!(one_signed, Some(Ok(_)));
        proof.attempts.push(VeroneseParentProofAttempt {
            route: VeroneseParentProofRoute::Attestation,
            variant: attestation.label.clone(),
            success,
            detail: match &one_signed {
                Some(Ok(signs)) => format!("{}; {signs}", attestation.statement),
                Some(Err(reason)) => format!(
                    "the parent equals the flagged entry, but the flag gives real-rootedness \
                     only and one-signed coefficients are not shown: {reason}"
                ),
                None => format!(
                    "the flagged entry starts at index {}, after the first used parent index \
                     {needed_from}",
                    attestation.first_index
                ),
            },
        });
        if success {
            proof.status = VeroneseParentProofStatus::ConditionalOnAttestation;
            proof.route = Some(VeroneseParentProofRoute::Attestation);
            proof.statement = format!(
                "{conclusion} (real-rootedness from the catalog flag, not proved here; \
                 one-signed coefficients proved here)"
            );
            proof.conditional_on = vec![attestation.label.clone()];
            return Ok(proof);
        }
    }
    if !differing.is_empty() {
        proof.attempts.push(VeroneseParentProofAttempt {
            route: VeroneseParentProofRoute::Attestation,
            variant: "oeis_catalog".into(),
            success: false,
            detail: format!(
                "no lean_proved catalog entry equals the parent; {} with the same first index \
                 and seeds differ in the rule or in later initial rows ({})",
                differing.len(),
                differing.join(", ")
            ),
        });
    }
    proof.statement = format!(
        "open: no route proved that {conclusion}; the first {} parent rows pass the exact \
         check (evidence only)",
        proof.rows_checked
    );
    Ok(proof)
}
