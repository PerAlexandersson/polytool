//! Image identity certification and parent proof checks.

use super::*;

// ---------------------------------------------------------------------------
// Certification
// ---------------------------------------------------------------------------

/// The outcome of the identity part of a certification.
pub(super) struct IdentityOutcome {
    pub(super) identity: PreserverImageIdentityStatus,
    pub(super) rows_compared: usize,
    pub(super) notes: Vec<String>,
}

fn supported(rows: usize, reason: String) -> PreserverImageIdentityStatus {
    PreserverImageIdentityStatus::SupportedByRows { rows, reason }
}

/// The first input row that differs from the candidate (an undefined
/// candidate row counts as a difference).
fn first_mismatch(input: &[Vec<Q>], candidate: &[Option<Vec<Q>>]) -> Option<(usize, bool)> {
    input
        .iter()
        .zip(candidate)
        .enumerate()
        .find_map(|(offset, (a, b))| match b {
            None => Some((offset, true)),
            Some(b) if vec_poly::trim(a) != *b => Some((offset, false)),
            Some(_) => None,
        })
}

fn refuted(n0: usize, offset: usize, undefined: bool) -> PreserverImageIdentityStatus {
    let index = n0 + offset;
    PreserverImageIdentityStatus::Refuted {
        index,
        reason: if undefined {
            format!("the operator is undefined on the parent row used for n = {index}")
        } else {
            format!("Q_{index} differs from the candidate lambda(n) T_n(P_(n+s))")
        },
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn certify_identity(
    input: &PreparedLiftRule,
    parent: &PreparedLiftRule,
    shifted: &PreparedLiftRule,
    candidate: &PreserverImageCandidate,
    form: &OperatorForm,
    gauge: &RatioGauge,
    options: &PreserverImageOptions,
) -> Result<IdentityOutcome, RecurrenceError> {
    let n0 = input.first_index;
    let mut notes = Vec::new();
    let mut count = options.checked_rows.max(input.seeds.len());
    let compare =
        |count: usize| -> Result<Option<PreserverImageIdentityStatus>, RecurrenceError> {
            let inputs = input.rows(count)?;
            let rows = candidate_rows(shifted, form, gauge, n0, count)?;
            Ok(first_mismatch(&inputs, &rows)
                .map(|(offset, undefined)| refuted(n0, offset, undefined)))
        };
    if let Some(status) = compare(count)? {
        return Ok(IdentityOutcome {
            identity: status,
            rows_compared: count,
            notes,
        });
    }
    if let Some(decision) = denominator_decision(form, n0) {
        if !decision.holds() {
            return Ok(IdentityOutcome {
                identity: supported(
                    count,
                    format!(
                        "the operator coefficients are not proved finite at every n >= {n0}: {}",
                        decision.describe("delta(n)")
                    ),
                ),
                rows_compared: count,
                notes,
            });
        }
    }

    // Hypergeometric multipliers: the Hadamard conjugation.
    if let OperatorForm::Hadamard { factor } = form {
        let hadamard_options = HadamardFactorLiftOptions {
            checked_rows: options
                .checked_rows
                .clamp(1, MAX_PRESERVER_IMAGE_CHECKED_ROWS),
            max_operator_terms: options
                .max_operator_terms
                .min(crate::recurrence::MAX_HADAMARD_OPERATOR_TERMS),
            ..HadamardFactorLiftOptions::default()
        };
        let (status, rows) =
            certify_candidate(input, shifted, factor, &candidate.gauge, &hadamard_options)?;
        count = count.max(rows);
        let identity = match status {
            HadamardFactorIdentityStatus::Certified { certificate } => {
                let wrapped = PreserverHadamardCertificate {
                    input: record(input),
                    parent: record(parent),
                    index_shift: candidate.index_shift,
                    factor: factor.clone(),
                    hadamard: *certificate,
                };
                wrapped.verify()?;
                PreserverImageIdentityStatus::Certified {
                    certificate: Box::new(PreserverImageCertificate::HadamardConjugation(
                        Box::new(wrapped),
                    )),
                }
            }
            HadamardFactorIdentityStatus::SupportedByRows { rows, reason } => supported(
                rows,
                format!("the Hadamard conjugation does not apply: {reason}"),
            ),
            HadamardFactorIdentityStatus::Refuted { index } => {
                PreserverImageIdentityStatus::Refuted {
                    index,
                    reason: format!("Q_{index} differs from the candidate Hadamard image"),
                }
            }
        };
        return Ok(IdentityOutcome {
            identity,
            rows_compared: count,
            notes,
        });
    }

    // Weyl operators and coordinate maps: pseudo-division of the dividend.
    let (k, _) = match x_free_leader(&shifted.relation) {
        Ok(found) => found,
        Err(error) => {
            return Ok(IdentityOutcome {
                identity: supported(count, format!("the parent relation is not usable: {error}")),
                rows_compared: count,
                notes,
            })
        }
    };
    let delta = form
        .index_denominator()
        .map(<[Q]>::to_vec)
        .unwrap_or_else(|| vec![Q::one()]);
    let (alpha, beta) = effective_ratio(gauge, &delta);
    let e =
        dividend(&input.relation, input.leader_shift, &alpha, &beta, form).ok_or_else(|| {
            RecurrenceError::VerificationFailed("the operator form has no dividend".into())
        })?;
    let division = match pseudo_divide_scalar(&e, &shifted.relation, options.max_operator_terms) {
        Ok(division) => division,
        Err(RecurrenceError::BudgetExhausted(reason) | RecurrenceError::LimitExceeded(reason)) => {
            return Ok(IdentityOutcome {
                identity: supported(count, format!("the division ran out of budget: {reason}")),
                rows_compared: count,
                notes,
            })
        }
        Err(error) => return Err(error),
    };
    let input_from = input.valid_from_base();
    let parent_from = shifted.valid_from_base().max(n0);
    let from = input_from.max(parent_from);
    let mut multiplier = division.multiplier.clone();
    let mut row_multiplier = Biv::one();
    let mut quotient_multiplier = Biv::one();
    let mut annihilators: Vec<(DOp, DOp)> = Vec::new();
    let mut annihilation_json = None;
    let mut extra_exceptional = Vec::new();
    if !division.remainder_is_zero() {
        let hadamard_options = HadamardFactorLiftOptions {
            max_operator_terms: options
                .max_operator_terms
                .min(crate::recurrence::MAX_HADAMARD_OPERATOR_TERMS),
            ..HadamardFactorLiftOptions::default()
        };
        if !remainder_is_small(&division.remainder[0]) {
            return Ok(IdentityOutcome {
                identity: supported(
                    count,
                    "the division leaves a nonzero remainder too large for the row-annihilator \
                     search"
                        .into(),
                ),
                rows_compared: count,
                notes,
            });
        }
        let annihilation =
            match remove_remainder(shifted, &division.remainder[0], from, &hadamard_options)? {
                Ok(annihilation) => annihilation,
                Err(reason) => {
                    notes.push(format!(
                        "the pseudo-division of the dividend by the shifted parent relation \
                         leaves a nonzero remainder of shift < {k}"
                    ));
                    return Ok(IdentityOutcome {
                        identity: supported(count, reason),
                        rows_compared: count,
                        notes,
                    });
                }
            };
        // The annihilation is `mu r = sum Z A` with `mu = c(n) p(x) mu_2(n, x)`.
        // From `w E = sum H S^nu L' + r`, multiplying by `mu` gives
        // `p mu_2 (c w) E = sum (mu H) S^nu L' + sum Z A`: the index content
        // `c(n)` goes into the multiplier `w` *and* into the quotients
        // (it was missing from the quotients, so the recorded identity
        // failed whenever `c != 1`).
        let index_content = n_polynomial(&annihilation.index_multiplier);
        multiplier = multiplier.multiply(&index_content);
        let (record, full, pairs) = annihilation_record(annihilation)?;
        if let PreserverRowAnnihilation::System { system } = &record {
            extra_exceptional = vanishing_bases(&biv_from_json(&system.row_multiplier)?, from);
        }
        quotient_multiplier = full.multiply(&index_content);
        row_multiplier = full;
        annihilators = pairs;
        annihilation_json = Some(record);
        notes.push(
            "the division remainder is removed by row annihilators of the shifted parent".into(),
        );
    }
    let quotients: BTreeMap<usize, DOp> = division
        .quotients
        .iter()
        .map(|(nu, row)| {
            let mut h = DOp::default();
            for (key, c) in &row[0].0 {
                h.add_term(*key, quotient_multiplier.multiply(c));
            }
            (*nu, h)
        })
        .filter(|(_, h)| !h.is_zero())
        .collect();
    if !identity_holds(
        &e,
        &shifted.relation,
        &row_multiplier.multiply(&multiplier),
        &quotients,
        &annihilators,
    ) {
        // A failed self-check rejects this candidate only (fail-safe: the
        // identity is then not certified).
        notes.push("the assembled division identity failed its exact re-check".into());
        return Ok(IdentityOutcome {
            identity: supported(
                count,
                "the pseudo-division identity does not hold exactly".into(),
            ),
            rows_compared: count,
            notes,
        });
    }
    let mut exceptional = exceptional_bases(&multiplier, input_from, parent_from);
    exceptional.extend(extra_exceptional);
    exceptional.sort_unstable();
    exceptional.dedup();
    let needed = rows_needed(&exceptional, input.leader_shift, n0, input.seeds.len());
    if needed > MAX_PRESERVER_IMAGE_CHECKED_ROWS {
        return Ok(IdentityOutcome {
            identity: supported(
                count,
                format!("the exceptional bases need {needed} compared rows"),
            ),
            rows_compared: count,
            notes,
        });
    }
    if needed > count {
        count = needed;
        if let Some(status) = compare(count)? {
            return Ok(IdentityOutcome {
                identity: status,
                rows_compared: count,
                notes,
            });
        }
    }
    let rows = candidate_rows(shifted, form, gauge, n0, needed)?;
    if let Err(RecurrenceError::VerificationFailed(reason)) = check_rows(input, &rows, &exceptional)
    {
        return Ok(IdentityOutcome {
            identity: supported(count, reason),
            rows_compared: count,
            notes,
        });
    }
    let identifying = input.seeds.len();
    let s = candidate.index_shift;
    let certificate = PreserverDivisionCertificate {
        input: record(input),
        parent: record(parent),
        index_shift: s,
        gauge: candidate.gauge.clone(),
        operator: candidate.operator.clone(),
        shifted_parent_relation: ImplicitRelation::from_op(&shifted.relation).to_json(),
        parent_leader_shift: k,
        dividend_relation: ImplicitRelation::from_op(&e).to_json(),
        multiplier: encode_poly(&n_coefficients(&multiplier)),
        quotients: quotients
            .iter()
            .map(|(nu, h)| PreserverDivisionQuotient {
                shift: *nu,
                differentiations: shift_free_json(h),
            })
            .collect(),
        row_annihilation: annihilation_json,
        input_valid_from_base: input_from,
        parent_valid_from_base: parent_from,
        exceptional_bases: exceptional,
        identifying_indices: (n0..n0 + identifying).collect(),
        statement: format!(
            "w(n) E = sum_nu H_nu S^nu L' (with any remainder removed by row annihilators of \
             the parent), where L' is the parent rule at \
             P'_n = P_(n{s:+}) (valid from base {parent_from}) and E the input rule gauged by \
             lambda(n)/delta(n) and composed with (or conjugated through) T = {}; so the \
             candidate rows satisfy the input rule at every base n >= {input_from} (the \
             exceptional bases checked exactly), and agreement on the {identifying} \
             identifying rows gives Q_n = lambda(n) T_n(P_(n{s:+})) for every n >= {n0}",
            candidate.operator.text()
        ),
    };
    match certificate.verify() {
        Ok(()) => {}
        Err(RecurrenceError::VerificationFailed(reason)) => {
            notes.push("the assembled certificate failed its replay".into());
            return Ok(IdentityOutcome {
                identity: supported(count, format!("the certificate does not replay: {reason}")),
                rows_compared: count,
                notes,
            });
        }
        Err(error) => return Err(error),
    }
    Ok(IdentityOutcome {
        identity: PreserverImageIdentityStatus::Certified {
            certificate: Box::new(PreserverImageCertificate::OperatorDivision(Box::new(
                certificate,
            ))),
        },
        rows_compared: count,
        notes,
    })
}

/// The parent provers: `prove_parent_with_policy` (rows zero or real-rooted with
/// coefficients of one sign) and, when that is not proved and the
/// operator's theorem needs only real-rootedness, the three-term theorem
/// without a root region.
pub(super) fn prove_parents(
    parent: &PreparedLiftRule,
    needed_from: usize,
    options: &PreserverImageOptions,
    whole_line_suffices: bool,
) -> Result<(VeroneseParentProof, Option<RealRootednessTheoremReport>), RecurrenceError> {
    let proof = prove_parent_with_policy(
        &ParentRecurrenceRecord {
            rule: &parent.rule,
            seeds: &parent.seeds,
            first_index: parent.first_index,
        },
        needed_from,
        &options.attestations,
        options.try_root_domain,
        options.allow_attestations,
    )?;
    if proof.status == VeroneseParentProofStatus::Proved || !whole_line_suffices {
        return Ok((proof, None));
    }
    let report = check_three_term_recurrence_theorem(
        &parent.rule,
        parent.first_index,
        &parent.seeds,
        &RealRootednessTheoremOptions {
            root_region: None,
            ..RealRootednessTheoremOptions::default()
        },
    )?;
    Ok((proof, report.is_proved().then_some(report)))
}

/// The parent theorem report carrying the parent's interlacing: the
/// whole-line report, else (for a parent proved by the parent provers) the
/// first of the three-term and Liu--Wang checks (roots in `(-inf, 0)`) and
/// the derivative-lag check that proves the parent with an interlacing
/// conclusion, else the interval report.
pub(super) fn parent_interlacing_report(
    parent: &PreparedLiftRule,
    proof: &VeroneseParentProof,
    whole_line: Option<&RealRootednessTheoremReport>,
    interval: Option<&RealRootednessTheoremReport>,
) -> Result<Option<RealRootednessTheoremReport>, RecurrenceError> {
    if let Some(report) = whole_line {
        return Ok(Some(report.clone()));
    }
    let interval = interval
        .filter(|report| {
            report.is_proved()
                && report
                    .conclusion
                    .as_ref()
                    .is_some_and(|c| c.interlacing.is_some())
        })
        .cloned();
    if proof.status != VeroneseParentProofStatus::Proved {
        return Ok(interval);
    }
    let negative = RealRootednessTheoremOptions {
        root_region: Some(crate::real_rootedness::RealInterval {
            lower: None,
            lower_closed: false,
            upper: Some(Q::zero()),
            upper_closed: false,
        }),
        ..RealRootednessTheoremOptions::default()
    };
    let (rule, p0, seeds) = (&parent.rule, parent.first_index, &parent.seeds);
    if proof.route == Some(VeroneseParentProofRoute::Attestation) {
        return Ok(interval);
    }
    let reports = [
        check_three_term_recurrence_theorem(rule, p0, seeds, &negative)?,
        check_liu_wang_theorem(rule, p0, seeds, &negative)?,
        check_liu_wang_derivative_lag_theorem(
            rule,
            p0,
            seeds,
            &RealRootednessTheoremOptions::default(),
        )?,
    ];
    Ok(reports
        .into_iter()
        .find(|report| {
            report.is_proved()
                && report
                    .conclusion
                    .as_ref()
                    .is_some_and(|c| c.interlacing.is_some())
        })
        .or(interval))
}
