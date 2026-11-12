//! Certificate replay and exact supporting checks.

use super::*;

// ---------------------------------------------------------------------------
// Replay
// ---------------------------------------------------------------------------

/// Replay an annihilation against the shifted parent from base `from`:
/// returns the row multiplier `p(x) mu_2(n, x)`, the pairs `(A, Z)` and the
/// bases where `mu_2` vanishes (exceptional).
fn replay_annihilation(
    annihilation: &PreserverRowAnnihilation,
    shifted: &PreparedLiftRule,
    from: usize,
) -> Result<(Biv, AnnihilatorTerms, Vec<usize>), RecurrenceError> {
    let failed = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
    match annihilation {
        PreserverRowAnnihilation::Single {
            row_multiplier,
            annihilator,
            quotient,
            step_multiplier,
            step_quotient,
            start_base,
            restart_bases,
        } => {
            let p = Biv::from_variable_polynomial(&decode_poly(row_multiplier)?);
            if p.is_zero() {
                return failed("the row multiplier is zero");
            }
            let a = shift_free_from_json(annihilator)?;
            match (step_multiplier, start_base) {
                (Some(nu), Some(start)) => {
                    let invariance = RowAnnihilatorInvariance {
                        step_multiplier: biv_from_json(nu)?,
                        step_quotient: shift_free_from_json(step_quotient)?,
                        start_base: *start,
                        restart_bases: restart_bases.clone(),
                    };
                    verify_invariance(&a, &invariance, shifted, from)?;
                }
                (None, None) if step_quotient.is_empty() && restart_bases.is_empty() => {
                    let Some(product) = ProductParent::recognise(shifted) else {
                        return failed("the shifted parent is not a product parent");
                    };
                    let expected = product.annihilator();
                    if expected != a || !product.annihilator_is_invariant(&a) {
                        return failed("the row annihilator does not match the parent");
                    }
                }
                _ => return failed("incomplete row-annihilator invariance fields"),
            }
            Ok((p, vec![(a, shift_free_from_json(quotient)?)], Vec::new()))
        }
        PreserverRowAnnihilation::System { system } => {
            let (decoded, terms) = system_from_json(system)?;
            verify_system(
                &decoded,
                system.start_base,
                &system.restart_bases,
                shifted,
                from,
            )?;
            let second = biv_from_json(&system.row_multiplier)?;
            if second.is_zero() {
                return failed("the row-annihilator system multiplier is zero");
            }
            let bases = vanishing_bases(&second, from);
            Ok((second, terms, bases))
        }
    }
}

/// The JSON form of an annihilation found by the Hadamard machinery, with
/// the full row multiplier and the pairs `(A, Z)`.
pub(super) fn annihilation_record(
    annihilation: RowAnnihilation,
) -> Result<(PreserverRowAnnihilation, Biv, AnnihilatorTerms), RecurrenceError> {
    let p = Biv::from_variable_polynomial(&annihilation.row_multiplier);
    match annihilation.kind {
        RowAnnihilationKind::Single {
            annihilator,
            quotient,
            invariance,
        } => {
            let record = PreserverRowAnnihilation::Single {
                row_multiplier: encode_poly(&annihilation.row_multiplier),
                annihilator: shift_free_json(&annihilator),
                quotient: shift_free_json(&quotient),
                step_multiplier: invariance.as_ref().map(|i| biv_json(&i.step_multiplier)),
                step_quotient: invariance
                    .as_ref()
                    .map_or_else(Vec::new, |i| shift_free_json(&i.step_quotient)),
                start_base: invariance.as_ref().map(|i| i.start_base),
                restart_bases: invariance.map_or_else(Vec::new, |i| i.restart_bases),
            };
            Ok((record, p, vec![(annihilator, quotient)]))
        }
        RowAnnihilationKind::System {
            system,
            row_multiplier,
            quotients,
            start_base,
            restart_bases,
        } => {
            let mut json = system_json(&system, &quotients, start_base, restart_bases);
            json.row_multiplier = biv_json(&row_multiplier);
            let (_, terms) = system_from_json(&json)?;
            Ok((
                PreserverRowAnnihilation::System { system: json },
                p.multiply(&row_multiplier),
                terms,
            ))
        }
    }
}

impl PreserverDivisionCertificate {
    /// Replay the certificate from its recorded input, parent, shift, gauge
    /// and operator: regenerate the input and parent rows from the recorded
    /// rules and seeds (and re-validate both rules), recompute `L'` and `E`,
    /// check `w E = sum_nu H_nu S^nu L'` exactly, check that every integer
    /// root of `w` from the valid bases on is listed as exceptional, and
    /// re-check the identifying rows and the input relation at every
    /// exceptional base on freshly computed candidate rows.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        let failed = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
        let input = prepare_record(&self.input, "input")?;
        let n0 = input.first_index;
        if self.identifying_indices != (n0..n0 + input.seeds.len()).collect::<Vec<_>>() {
            return failed("the identifying indices do not match the input");
        }
        let parent = prepare_record(&self.parent, "parent")?;
        let shifted = shifted_parent(&parent, n0, self.index_shift)?;
        if ImplicitRelation::from_op(&shifted.relation).to_json() != self.shifted_parent_relation {
            return failed("the shifted parent relation does not match its recomputation");
        }
        let (k, leader) = x_free_leader(&shifted.relation)?;
        if k != self.parent_leader_shift {
            return failed("the parent leader shift does not match");
        }
        check_schedule_domain(&self.operator, n0)?;
        let form = self.operator.form()?;
        if matches!(form, OperatorForm::Hadamard { .. }) {
            return failed("a Hadamard multiplier needs the Hadamard conjugation certificate");
        }
        if let Some(decision) = denominator_decision(&form, n0) {
            if !decision.holds() {
                return failed("the operator denominator is not proved nonzero on the domain");
            }
        }
        let gauge = self.gauge.ratio_form(n0)?;
        let delta = form.index_denominator().unwrap_or(&[]).to_vec();
        let delta = if delta.is_empty() {
            vec![Q::one()]
        } else {
            delta
        };
        let (alpha, beta) = effective_ratio(&gauge, &delta);
        let Some(e) = dividend(&input.relation, input.leader_shift, &alpha, &beta, &form) else {
            return failed("the operator has no dividend");
        };
        if ImplicitRelation::from_op(&e).to_json() != self.dividend_relation {
            return failed("the dividend does not match its recomputation");
        }
        let multiplier = n_polynomial(&decode_poly(&self.multiplier)?);
        if multiplier.is_zero() {
            return failed("the multiplier is zero");
        }
        let mut quotients: BTreeMap<usize, DOp> = BTreeMap::new();
        for entry in &self.quotients {
            let h = shift_free_from_json(&entry.differentiations)?;
            let slot = quotients.entry(entry.shift).or_default();
            *slot = slot.add(&h);
        }
        let _ = leader;
        let input_from = input.valid_from_base();
        let parent_from = shifted.valid_from_base().max(n0);
        if input_from != self.input_valid_from_base || parent_from != self.parent_valid_from_base {
            return failed("the recorded valid bases do not match");
        }
        let from = input_from.max(parent_from);
        let (row_multiplier, annihilators, mut expected) = match &self.row_annihilation {
            None => (Biv::one(), Vec::new(), Vec::new()),
            Some(annihilation) => replay_annihilation(annihilation, &shifted, from)?,
        };
        if !identity_holds(
            &e,
            &shifted.relation,
            &row_multiplier.multiply(&multiplier),
            &quotients,
            &annihilators,
        ) {
            return failed("the identity w E = sum H_nu S^nu L' (+ Z A) fails");
        }
        expected.extend(exceptional_bases(&multiplier, input_from, parent_from));
        if expected
            .iter()
            .any(|base| !self.exceptional_bases.contains(base))
        {
            return failed("an exceptional base is missing");
        }
        let needed = rows_needed(
            &self.exceptional_bases,
            input.leader_shift,
            n0,
            input.seeds.len(),
        );
        if needed > MAX_PRESERVER_IMAGE_CHECKED_ROWS {
            return Err(RecurrenceError::LimitExceeded(format!(
                "the certificate needs {needed} candidate rows, more than \
                 {MAX_PRESERVER_IMAGE_CHECKED_ROWS}"
            )));
        }
        let rows = candidate_rows(&shifted, &form, &gauge, n0, needed)?;
        check_rows(&input, &rows, &self.exceptional_bases)
    }
}

impl PreserverHadamardCertificate {
    /// Replay: [`HadamardFactorCertificate::verify`], the agreement of its
    /// factor, gauge, input relation and parent with the records here (the
    /// parent being the shifted parent `P'` recomputed from the parent
    /// record), and the identifying rows and exceptional bases re-checked
    /// on candidate rows computed from the regenerated parent rows.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        let failed = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
        self.hadamard.verify()?;
        let input = prepare_record(&self.input, "input")?;
        let n0 = input.first_index;
        let parent = prepare_record(&self.parent, "parent")?;
        let shifted = shifted_parent(&parent, n0, self.index_shift)?;
        if self.hadamard.factor != self.factor {
            return failed("the Hadamard certificate has another factor");
        }
        if self.hadamard.parent != record(&shifted) {
            return failed("the Hadamard certificate's parent is not the shifted parent");
        }
        if self.hadamard.input_relation != ImplicitRelation::from_op(&input.relation).to_json()
            || self.hadamard.input_leader_shift != input.leader_shift
        {
            return failed("the Hadamard certificate's input relation differs from the input");
        }
        if self.hadamard.identifying_indices != (n0..n0 + input.seeds.len()).collect::<Vec<_>>() {
            return failed("the identifying indices do not match the input");
        }
        if self.hadamard.input_valid_from_base != input.valid_from_base()
            || self.hadamard.parent_valid_from_base != shifted.valid_from_base().max(n0)
        {
            return failed("the recorded valid bases do not match");
        }
        let gauge = self.hadamard.gauge.ratio_form(n0)?;
        let form = OperatorForm::Hadamard {
            factor: self.factor.clone(),
        };
        let needed = rows_needed(
            &self.hadamard.exceptional_bases,
            input.leader_shift,
            n0,
            input.seeds.len(),
        );
        if needed > MAX_PRESERVER_IMAGE_CHECKED_ROWS {
            return Err(RecurrenceError::LimitExceeded(format!(
                "the certificate needs {needed} candidate rows"
            )));
        }
        let rows = candidate_rows(&shifted, &form, &gauge, n0, needed)?;
        check_rows(&input, &rows, &self.hadamard.exceptional_bases)
    }
}

impl PreserverImageCertificate {
    /// Replay the certificate ([`PreserverDivisionCertificate::verify`] or
    /// [`PreserverHadamardCertificate::verify`]).
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        match self {
            Self::OperatorDivision(certificate) => certificate.verify(),
            Self::HadamardConjugation(certificate) => certificate.verify(),
        }
    }

    /// The recorded input, parent, shift, gauge and operator.
    fn parts(
        &self,
    ) -> (
        &RecurrenceJson,
        &RecurrenceJson,
        i64,
        &VeroneseSectionGauge,
        PreserverOperator,
    ) {
        match self {
            Self::OperatorDivision(c) => (
                &c.input,
                &c.parent,
                c.index_shift,
                &c.gauge,
                c.operator.clone(),
            ),
            Self::HadamardConjugation(c) => (
                &c.input,
                &c.parent,
                c.index_shift,
                &c.hadamard.gauge,
                PreserverOperator::HadamardMultiplier {
                    factor: c.factor.clone(),
                },
            ),
        }
    }
}

impl PreserverImageCertification {
    /// Whether every `Q_n`, `n >= n0`, is proved zero or real-rooted in the
    /// summary's input root region.
    pub fn is_proved(&self) -> bool {
        self.lift_summary.is_proved()
    }

    /// The [`LiftSummary`] computed from the record.
    pub fn recompute_lift_summary(&self) -> LiftSummary {
        preserver_image_summary(PreserverImageSummaryParts {
            candidate: &self.candidate,
            input_first_index: self.input_first_index,
            identity: &self.identity,
            parent_proof: &self.parent_proof,
            whole_line_parent_proof: self.whole_line_parent_proof.as_ref(),
            interval_parent_proof: self.interval_parent_proof.as_ref(),
            interlacing: &self.interlacing,
        })
    }

    /// Replay the certificate (when certified), check that it belongs to the
    /// candidate, replay the whole-line parent theorem report (when present),
    /// and check that the summary equals the recomputed one.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        let failed = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
        if self.schema != PRESERVER_IMAGE_CERTIFICATION_SCHEMA {
            return Err(RecurrenceError::InvalidInput(format!(
                "unsupported schema `{}`",
                self.schema
            )));
        }
        if let PreserverImageIdentityStatus::Certified { certificate } = &self.identity {
            certificate.verify()?;
            let (input, parent, shift, gauge, operator) = certificate.parts();
            let candidate_parent = record(&prepare_record(&self.candidate.parent, "parent")?);
            if *parent != candidate_parent
                || shift != self.candidate.index_shift
                || *gauge != self.candidate.gauge
                || operator != self.candidate.operator
                || input.first_index != self.input_first_index
            {
                return failed("the certificate does not belong to the candidate");
            }
        }
        if let Some(report) = &self.whole_line_parent_proof {
            report.verify()?;
            let candidate_parent = record(&prepare_record(&self.candidate.parent, "parent")?);
            if report.recurrence != candidate_parent || report.requested_root_region.is_some() {
                return failed("the whole-line parent report does not belong to the parent");
            }
        }
        if let Some(report) = &self.interval_parent_proof {
            report.verify()?;
            let candidate_parent = record(&prepare_record(&self.candidate.parent, "parent")?);
            let bounded = report
                .conclusion
                .as_ref()
                .and_then(|c| c.root_region.as_ref())
                .is_some_and(|r| r.lower.is_some() && r.upper.is_some());
            if report.recurrence != candidate_parent || !report.is_proved() || !bounded {
                return failed("the interval parent report does not belong to the parent");
            }
        }
        if let Some(report) = &self.parent_interlacing_report {
            report.verify()?;
            let candidate_parent = record(&prepare_record(&self.candidate.parent, "parent")?);
            if report.recurrence != candidate_parent || !report.is_proved() {
                return failed("the parent interlacing report does not belong to the parent");
            }
        }
        self.verify_image_rule()?;
        let claim = interlacing_claim(
            &self.candidate,
            self.input_first_index,
            &self.identity,
            self.parent_interlacing_report.as_ref(),
        );
        if claim != self.interlacing {
            return failed("the interlacing claim does not match the one recomputed");
        }
        let summary = self.recompute_lift_summary();
        if summary != self.lift_summary || !summary.is_consistent() {
            return failed("the lift summary does not match the one recomputed from the record");
        }
        Ok(())
    }
}

impl PreserverImageCertification {
    /// Replay the computed image rule (when present): its certificate, that
    /// it belongs to the candidate and the input's first index, that the
    /// identity is certified, and that its first rows are the input rows.
    fn verify_image_rule(&self) -> Result<(), RecurrenceError> {
        let Some(image) = &self.image_rule else {
            return Ok(());
        };
        let failed = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
        let PreserverImageIdentityStatus::Certified { certificate } = &self.identity else {
            return failed("an image rule is recorded without a certified identity");
        };
        image.verify()?;
        if image.certificate.parent != self.candidate.parent
            || image.certificate.first_index != self.input_first_index
            || image.certificate.specification != self.candidate.image_operator_specification()?
        {
            return failed("the image rule does not belong to the candidate");
        }
        let (input, ..) = certificate.parts();
        let input = prepare_record(input, "input")?;
        let count = input.seeds.len() + image.order + 1;
        if !image_rows_match(image, &input, count)? {
            return failed("the image rule disagrees with the input rows");
        }
        Ok(())
    }
}
