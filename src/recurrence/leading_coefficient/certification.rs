//! Leading-coefficient route certification and exact leading values.

use super::*;

impl Context {
    /// The certificate of `route` for `law`.
    pub(super) fn certify(
        &self,
        law: &RecurrenceDegreeLaw,
        route: LeadingCoefficientRoute,
        ratio: Option<&Q>,
        rows: &[Vec<Q>],
    ) -> Result<LeadingCoefficientCertificate, RecurrenceError> {
        let start = self.induction_start(law);
        let mut certificate = self.not_applicable(self.normal_form.clone());
        certificate.status = LeadingCoefficientStatus::NotProved;
        certificate.route = Some(route);
        certificate.degree_law = Some(law.clone());
        certificate.degree_is_exact = route != LeadingCoefficientRoute::Coefficients;
        certificate.ratio = ratio.cloned();
        certificate.induction_start = start;
        if law.start < self.first_index
            || (self.shape.alternating && law.period % 2 == 1)
            || rows.len() < start - self.first_index
        {
            return Err(RecurrenceError::InvalidInput(
                "the degree law starts before the first index, has an odd period for a rule \
                 with (-1)^n factors, or the rows are missing"
                    .into(),
            ));
        }
        let initial = &rows[..start - self.first_index];
        certificate.initial_rows = initial
            .iter()
            .enumerate()
            .map(|(i, row)| LeadingCoefficientRow {
                index: self.first_index + i,
                degree: row_degree(row),
                sign: leading_coefficient(row).and_then(LeadingCoefficientSign::of),
            })
            .collect();
        let mut hypotheses = Vec::new();
        hypotheses.push(self.initial_rows_hypothesis(law, route, ratio, initial));
        if self.shape.left_factor.len() > 1 {
            hypotheses.push(sign_hypothesis(
                "left_factor_nonzero",
                &format!("f(n) != 0 for every n >= {start}"),
                "f(n)",
                on_class(&self.shape.left_factor, start, 1),
                vec![Q::one()],
                start,
                1,
                IndexSignCondition::Nonzero,
            ));
        }
        let classes = self.classes(law);
        match route {
            LeadingCoefficientRoute::Cone
            | LeadingCoefficientRoute::Subseq
            | LeadingCoefficientRoute::TwoStep => {
                for class in &classes {
                    hypotheses.push(self.degree_bound_hypothesis(law, class, 0));
                    hypotheses.extend(self.cone_hypotheses(law, class));
                }
            }
            LeadingCoefficientRoute::Recurrence => {
                let Some(rho) = ratio.filter(|rho| rho.is_positive()) else {
                    return Err(RecurrenceError::InvalidInput(
                        "the recurrence route needs a positive ratio".into(),
                    ));
                };
                for class in &classes {
                    hypotheses.push(self.degree_bound_hypothesis(law, class, 0));
                    hypotheses.extend(self.ratio_hypotheses(law, class, rho));
                }
            }
            LeadingCoefficientRoute::Cancel => {
                if ratio.is_some_and(|rho| !rho.is_positive()) {
                    return Err(RecurrenceError::InvalidInput(
                        "the ratio of the cancel route must be positive".into(),
                    ));
                }
                hypotheses.extend(self.cancel_hypotheses(law, &classes, ratio));
            }
            LeadingCoefficientRoute::Coefficients => {
                for class in &classes {
                    hypotheses.push(self.degree_bound_hypothesis(law, class, 0));
                    hypotheses.extend(self.coefficient_hypotheses(law, class, true));
                }
            }
        }
        // The route label must fit the rule.
        let label_fits = match route {
            LeadingCoefficientRoute::Subseq => {
                law.period > 1 && self.shape.terms.iter().all(|t| t.lag % law.period == 0)
            }
            LeadingCoefficientRoute::TwoStep => {
                self.is_two_step() && law.has_constant_sign(LeadingCoefficientSign::Positive)
            }
            // Nonnegative coefficients: at most first derivatives, no
            // `(-1)^n` factor, positive leading coefficients (and first
            // order for `cancel`).
            LeadingCoefficientRoute::Cancel | LeadingCoefficientRoute::Coefficients => {
                !self.shape.alternating
                    && self.shape.terms.iter().all(|t| t.derivative_order <= 1)
                    && law.has_constant_sign(LeadingCoefficientSign::Positive)
                    && (route == LeadingCoefficientRoute::Coefficients || self.shape.max_lag == 1)
            }
            _ => true,
        };
        if !label_fits {
            hypotheses.push(LiftHypothesis {
                name: "route_shape".into(),
                statement: format!(
                    "the rule has the shape of the route {}",
                    route.lean_route_name()
                ),
                status: LiftHypothesisStatus::Failed { witness: None },
                evidence: "the rule does not have this shape".into(),
                sign_decision: None,
            });
        }
        certificate.hypotheses = hypotheses;
        if certificate
            .hypotheses
            .iter()
            .all(|h| h.status == LiftHypothesisStatus::Verified)
        {
            certificate.status = LeadingCoefficientStatus::Proved;
            let route = self.refined_route(law, route);
            certificate.route = Some(route);
            certificate.lean_route = Some(route.lean_route_name().into());
            certificate.lean_declarations = lean_declarations(route, law, ratio.is_some());
            let holds_from = certificate
                .initial_rows
                .iter()
                .rposition(|row| row.degree.is_none())
                .map_or(self.first_index, |i| self.first_index + i + 1);
            certificate.holds_from = Some(holds_from);
            certificate.conclusion = Some(self.conclusion(law, route, holds_from, &certificate));
            if certificate.degree_is_exact {
                if let Some((value, hypotheses)) = self.exact_leading_coefficient(law, rows) {
                    if let Some(conclusion) = certificate.conclusion.as_mut() {
                        conclusion.push_str(&format!(
                            "; lead P_n = {} for every n >= {}",
                            format_rational_coeff(&value),
                            start - self.shape.max_lag
                        ));
                    }
                    certificate.exact_leading_coefficient = Some(value);
                    certificate.exact_leading_coefficient_hypotheses = hypotheses;
                }
            }
        }
        Ok(certificate)
    }

    /// The constant leading coefficient `c` (see
    /// [`LeadingCoefficientCertificate::exact_leading_coefficient`]) with
    /// its verified hypotheses, when every one holds.
    pub(super) fn exact_leading_coefficient(
        &self,
        law: &RecurrenceDegreeLaw,
        rows: &[Vec<Q>],
    ) -> Option<(Q, Vec<LiftHypothesis>)> {
        if !(law.has_constant_sign(LeadingCoefficientSign::Positive)
            || law.has_constant_sign(LeadingCoefficientSign::Negative))
        {
            return None;
        }
        let start = self.induction_start(law);
        let window = start - self.shape.max_lag;
        let value = leading_coefficient(&rows[window - self.first_index])?.clone();
        if (window..start).any(|n| leading_coefficient(&rows[n - self.first_index]) != Some(&value))
        {
            return None;
        }
        let mut hypotheses = Vec::new();
        for class in &self.classes(law) {
            // The top-coefficient recurrence needs the degree bound on every
            // class (it fails where the `cancel` route's top terms cancel:
            // there the coefficient also involves `[x^(D-1)] P_(n-1)`).
            let bound = self.degree_bound_hypothesis(law, class, 0);
            if bound.status != LiftHypothesisStatus::Verified {
                return None;
            }
            hypotheses.push(bound);
            let denominator = self.left_factor_on(class, law);
            let total = self
                .signed_multipliers(law, class)
                .iter()
                .fold(Vec::new(), |sum, m| vec_poly::add(&sum, m));
            let hypothesis = sign_hypothesis(
                &format!("top_multiplier_sum_is_one_class_{}", class.c),
                &format!(
                    "sum_j M_j(n) = f(n) for n = {} + {} m, so lead P_n = lead P_(n-j) for \
                     every lag",
                    class.c, law.period
                ),
                "sum_j M_j - f",
                vec_poly::sub(&total, &denominator).stripped(),
                vec![Q::one()],
                class.c,
                law.period,
                IndexSignCondition::Zero,
            );
            if hypothesis.status != LiftHypothesisStatus::Verified {
                return None;
            }
            hypotheses.push(hypothesis);
        }
        Some((value, hypotheses))
    }

    fn is_two_step(&self) -> bool {
        matches!(self.shape.terms.as_slice(), [term] if term.lag == 2
            && term.derivative_order == 0
            && !term.alternating
            && term.coefficients.iter().skip(1).all(|row| row.iter().all(Zero::is_zero)))
            && self.shape.left_factor.len() == 1
    }

    /// `cone` becomes `twoStep` or `subseq` when the rule has that shape.
    fn refined_route(
        &self,
        law: &RecurrenceDegreeLaw,
        route: LeadingCoefficientRoute,
    ) -> LeadingCoefficientRoute {
        match route {
            LeadingCoefficientRoute::Cone
                if self.is_two_step()
                    && law.has_constant_sign(LeadingCoefficientSign::Positive) =>
            {
                LeadingCoefficientRoute::TwoStep
            }
            LeadingCoefficientRoute::Cone
                if law.period > 1 && self.shape.terms.iter().all(|t| t.lag % law.period == 0) =>
            {
                LeadingCoefficientRoute::Subseq
            }
            other => other,
        }
    }

    pub(super) fn conclusion(
        &self,
        law: &RecurrenceDegreeLaw,
        route: LeadingCoefficientRoute,
        holds_from: usize,
        certificate: &LeadingCoefficientCertificate,
    ) -> String {
        let relation = if certificate.degree_is_exact {
            "the degree law"
        } else {
            "the degree envelope (deg P_n <= D(n)), nonnegative coefficients and P_n(0) > 0, with"
        };
        let mut text = format!(
            "{relation} {} (route {}, induction from n = {})",
            law.text(),
            route.lean_route_name(),
            certificate.induction_start
        );
        if law.start > self.first_index {
            let earlier: Vec<String> = certificate
                .initial_rows
                .iter()
                .filter(|row| row.index < law.start)
                .map(|row| match (row.degree, row.sign) {
                    (Some(degree), Some(sign)) => {
                        format!("P_{}: degree {degree}, {}", row.index, sign.text())
                    }
                    _ => format!("P_{} = 0", row.index),
                })
                .collect();
            text.push_str(&format!("; earlier rows (exact): {}", earlier.join(", ")));
        }
        if holds_from > self.first_index {
            text.push_str(&format!("; every row from P_{holds_from} on is nonzero"));
        } else {
            text.push_str("; every row is nonzero");
        }
        text
    }

    // -- hypotheses -------------------------------------------------------

    fn initial_rows_hypothesis(
        &self,
        law: &RecurrenceDegreeLaw,
        route: LeadingCoefficientRoute,
        ratio: Option<&Q>,
        initial: &[Vec<Q>],
    ) -> LiftHypothesis {
        let start = self.induction_start(law);
        let window = start - self.shape.max_lag;
        let exact = route != LeadingCoefficientRoute::Coefficients;
        let needs_nonnegative = matches!(
            route,
            LeadingCoefficientRoute::Coefficients | LeadingCoefficientRoute::Cancel
        );
        let statement = format!(
            "the rows P_{}, ..., P_{} computed exactly: every row P_n with n >= {} has {} D(n) and \
             the sign of the law{}{}",
            self.first_index,
            start - 1,
            law.start,
            if exact { "degree" } else { "degree at most" },
            if needs_nonnegative {
                format!("; the rows from P_{window} on have nonnegative coefficients")
            } else {
                String::new()
            },
            if route == LeadingCoefficientRoute::Coefficients {
                " and a positive constant term"
            } else {
                ""
            },
        );
        let mut failure: Option<(usize, String)> = None;
        for (i, row) in initial.iter().enumerate() {
            let n = self.first_index + i;
            let Some(expected) = law.degree_at(n) else {
                continue;
            };
            let degree = row_degree(row);
            let sign = leading_coefficient(row).and_then(LeadingCoefficientSign::of);
            let degree_fits = match degree {
                None => false,
                Some(d) if exact => d == expected,
                Some(d) => d <= expected,
            };
            if !degree_fits || sign != law.sign_at(n) {
                failure = Some((n, format!("P_{n} has degree {degree:?} and sign {sign:?}")));
                break;
            }
            if needs_nonnegative && n >= window {
                if row.iter().any(Signed::is_negative) {
                    failure = Some((n, format!("P_{n} has a negative coefficient")));
                    break;
                }
                if route == LeadingCoefficientRoute::Coefficients
                    && !row.first().is_some_and(Signed::is_positive)
                {
                    failure = Some((n, format!("P_{n}(0) is not positive")));
                    break;
                }
            }
        }
        if failure.is_none() {
            if let (LeadingCoefficientRoute::Cancel, Some(rho)) = (route, ratio) {
                let n = start - 1;
                let row = &initial[n - self.first_index];
                let holds = row_degree(row)
                    .filter(|&d| d >= 1)
                    .is_some_and(|d| rho * &row[d] <= row[d - 1]);
                if !holds {
                    failure = Some((
                        n,
                        format!("the ratio base fails: [x^(D-1)] P_{n} < rho lead P_{n}"),
                    ));
                }
            }
            if let (LeadingCoefficientRoute::Recurrence, Some(rho)) = (route, ratio) {
                let values: Vec<Q> = (window..start)
                    .map(|n| {
                        let row = &initial[n - self.first_index];
                        leading_coefficient(row).cloned().unwrap_or_else(Q::zero)
                            * law.sign_at(n).map_or_else(Q::zero, |s| s.factor())
                    })
                    .collect();
                if let Some(t) = values.windows(2).position(|w| rho * &w[0] > w[1]) {
                    failure = Some((
                        window + t + 1,
                        format!(
                            "the ratio base fails: rho mu({}) > mu({}) for mu(n) = sign * lead P_n",
                            window + t,
                            window + t + 1
                        ),
                    ));
                }
            }
        }
        let (status, evidence) = match failure {
            None => (
                LiftHypothesisStatus::Verified,
                format!(
                    "checked exactly on {} rows{}",
                    initial.len(),
                    match (route, ratio) {
                        (LeadingCoefficientRoute::Recurrence, Some(rho)) => format!(
                            "; rho mu(t) <= mu(t + 1) for t = {window}, ..., {} with rho = {}",
                            start.saturating_sub(2),
                            format_rational_coeff(rho)
                        ),
                        (LeadingCoefficientRoute::Cancel, Some(rho)) => format!(
                            "; [x^(D-1)] P_{0} >= rho lead P_{0} with rho = {1}",
                            start - 1,
                            format_rational_coeff(rho)
                        ),
                        _ => String::new(),
                    }
                ),
            ),
            Some((n, reason)) => (
                LiftHypothesisStatus::Failed {
                    witness: Some(BigInt::from(n)),
                },
                format!("{reason} (exact)"),
            ),
        };
        LiftHypothesis {
            name: "initial_rows".into(),
            statement,
            status,
            evidence,
            sign_decision: None,
        }
    }
}
