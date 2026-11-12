//! Degree, sign, cone, ratio and cancellation hypotheses.

use super::*;

impl Context {
    /// The degree bound `deg_x c_t <= e_j + d_t + extra` and `e_j >= 0` on
    /// a class.
    pub(super) fn degree_bound_hypothesis(
        &self,
        law: &RecurrenceDegreeLaw,
        class: &Class,
        extra: i64,
    ) -> LiftHypothesis {
        let c = class.c;
        let name = format!("degree_bound_class_{c}");
        let statement = format!(
            "for n = {c} + {} m: D(n - j) <= D(n) for every lag j, and deg_x c_t(n, x) <= \
             D(n) - D(n - j_t) + d_t{} for every summand",
            law.period,
            if extra > 0 { " + 1" } else { "" }
        );
        let mut failure = None;
        for (lag, order, coefficients) in self.class_terms(class) {
            let gap = class.gaps[&lag];
            if gap < 0 {
                failure = Some((None, format!("D({c}) < D({}) (lag {lag})", c - lag)));
                break;
            }
            let allowed = gap + order as i64 + extra;
            let top = x_degree(&coefficients).unwrap_or(0) as i64;
            if top > allowed {
                // The coefficient of x^top is a nonzero polynomial in n: it
                // is nonzero at some n of the class.
                let excess = column(&coefficients, top as usize);
                let witness = (0..=excess.len()).find_map(|m| {
                    let n = c + law.period * m;
                    (!vec_poly::evaluate(&excess, &integer(n as i64)).is_zero()).then_some(n)
                });
                failure = Some((
                    witness,
                    format!(
                        "the summand of lag {lag} and derivative order {order} has degree \
                         {top} > {allowed} in x"
                    ),
                ));
                break;
            }
        }
        let (status, evidence) = match failure {
            None => (
                LiftHypothesisStatus::Verified,
                format!(
                    "gaps {} (exact)",
                    class
                        .gaps
                        .iter()
                        .map(|(j, e)| format!("e_{j} = {e}"))
                        .collect::<Vec<_>>()
                        .join(", ")
                ),
            ),
            Some((witness, reason)) => (
                LiftHypothesisStatus::Failed {
                    witness: witness.map(BigInt::from),
                },
                format!("{reason} (exact)"),
            ),
        };
        LiftHypothesis {
            name,
            statement,
            status,
            evidence,
            sign_decision: None,
        }
    }

    /// The numerators (over `f(c + period m)`) of the signed multipliers
    /// `m_j` on a class, by lag `1..=K`.
    pub(super) fn signed_multipliers(
        &self,
        law: &RecurrenceDegreeLaw,
        class: &Class,
    ) -> Vec<Vec<Q>> {
        let c = class.c;
        let mut multipliers = vec![Vec::new(); self.shape.max_lag + 1];
        for (lag, order, coefficients) in self.class_terms(class) {
            let gap = class.gaps[&lag];
            if gap < 0 {
                continue;
            }
            let top = column(&coefficients, gap as usize + order);
            let factor = falling_factorial_on_class(
                law.degree_at(c - lag).unwrap(),
                law.growth_per_period,
                order,
            );
            let sign = law.sign_at(c).unwrap().factor() * law.sign_at(c - lag).unwrap().factor();
            let summand = vec_poly::scale(
                &vec_poly::mul(&on_class(&top, c, law.period), &factor),
                &sign,
            );
            multipliers[lag] = vec_poly::add(&multipliers[lag], &summand).stripped();
        }
        multipliers
    }

    pub(super) fn cone_hypotheses(
        &self,
        law: &RecurrenceDegreeLaw,
        class: &Class,
    ) -> Vec<LiftHypothesis> {
        let c = class.c;
        let denominator = self.left_factor_on(class, law);
        let multipliers = self.signed_multipliers(law, class);
        let mut hypotheses = Vec::new();
        let mut total = Vec::new();
        for (lag, numerator) in multipliers.iter().enumerate().skip(1) {
            if !class.gaps.contains_key(&lag) {
                continue;
            }
            total = vec_poly::add(&total, numerator);
            hypotheses.push(sign_hypothesis(
                &format!("top_multiplier_lag_{lag}_nonnegative_class_{c}"),
                &format!(
                    "the signed top-coefficient multiplier m_{lag}(n) >= 0 for n = {c} + {} m",
                    law.period
                ),
                &format!("m_{lag}"),
                numerator.clone(),
                denominator.clone(),
                c,
                law.period,
                IndexSignCondition::Nonnegative,
            ));
        }
        hypotheses.push(sign_hypothesis(
            &format!("top_multiplier_sum_positive_class_{c}"),
            &format!(
                "sum_j m_j(n) > 0 for n = {c} + {} m (the top coefficients stay positive)",
                law.period
            ),
            "sum_j m_j",
            total.stripped(),
            denominator,
            c,
            law.period,
            IndexSignCondition::Positive,
        ));
        hypotheses
    }

    pub(super) fn ratio_hypotheses(
        &self,
        law: &RecurrenceDegreeLaw,
        class: &Class,
        rho: &Q,
    ) -> Vec<LiftHypothesis> {
        let c = class.c;
        let k = self.shape.max_lag;
        let denominator = self.left_factor_on(class, law);
        let multipliers = self.signed_multipliers(law, class);
        let power = |e: usize| (0..e).fold(Q::one(), |acc, _| acc * rho);
        (0..k)
            .map(|i| {
                // C_i = sum_(t = i)^(K-1) m_(K-t) rho^(t-i) - rho^(K-i), over f.
                let mut numerator = vec_poly::scale(&denominator, &-power(k - i));
                for t in i..k {
                    numerator = vec_poly::add(
                        &numerator,
                        &vec_poly::scale(&multipliers[k - t], &power(t - i)),
                    );
                }
                sign_hypothesis(
                    &format!("ratio_step_{i}_class_{c}"),
                    &format!(
                        "C_{i}(n) = sum_(t = {i})^{} m_(K-t)(n) rho^(t-{i}) - rho^{} >= 0 for \
                         n = {c} + {} m, rho = {} (the ratio invariant mu(n) >= rho mu(n - 1) \
                         propagates)",
                        k - 1,
                        k - i,
                        law.period,
                        format_rational_coeff(rho)
                    ),
                    &format!("C_{i}"),
                    numerator.stripped(),
                    denominator.clone(),
                    c,
                    law.period,
                    IndexSignCondition::Nonnegative,
                )
            })
            .collect()
    }

    /// The pairs `(B_j, A_j)` of non-derivative and first-derivative
    /// coefficients by lag, on a class.
    fn paired_terms(&self, class: &Class) -> BTreeMap<usize, (Bivariate, Bivariate)> {
        let mut pairs: BTreeMap<usize, (Bivariate, Bivariate)> = BTreeMap::new();
        for (lag, order, coefficients) in self.class_terms(class) {
            let entry = pairs.entry(lag).or_default();
            let target = if order == 0 {
                &mut entry.0
            } else {
                &mut entry.1
            };
            *target = add_bivariate(target, &coefficients);
        }
        pairs
    }

    /// Nonnegative coefficients: `[x^0] A_j >= 0` and `B_(j,i) + r A_(j,i+1)
    /// >= 0` at `r = 0` and `r = D(n - j)`; with `witness`, also
    /// `sum_j [x^0] B_j > 0`.
    pub(super) fn coefficient_hypotheses(
        &self,
        law: &RecurrenceDegreeLaw,
        class: &Class,
        witness: bool,
    ) -> Vec<LiftHypothesis> {
        let c = class.c;
        let p = law.period;
        let denominator = self.left_factor_on(class, law);
        let mut hypotheses = Vec::new();
        let mut constant_terms = Vec::new();
        for (lag, (b, a)) in self.paired_terms(class) {
            let on = |poly: Vec<Q>| on_class(&poly, c, p);
            let previous = vec![
                integer(law.degree_at(c - lag).unwrap() as i64),
                integer(law.growth_per_period as i64),
            ];
            let a0 = on(column(&a, 0));
            if !a0.is_empty() {
                hypotheses.push(sign_hypothesis(
                    &format!("derivative_constant_nonnegative_lag_{lag}_class_{c}"),
                    &format!("[x^0] A_{lag}(n) >= 0 for n = {c} + {p} m"),
                    &format!("[x^0] A_{lag}"),
                    a0,
                    denominator.clone(),
                    c,
                    p,
                    IndexSignCondition::Nonnegative,
                ));
            }
            let top = x_degree(&b).max(x_degree(&a).map(|d| d.saturating_sub(1)));
            for i in 0..=top.unwrap_or(0) {
                let bi = on(column(&b, i));
                let ai = on(column(&a, i + 1));
                if bi.is_empty() && ai.is_empty() {
                    continue;
                }
                if i == 0 {
                    constant_terms = vec_poly::add(&constant_terms, &bi);
                }
                if !bi.is_empty() {
                    hypotheses.push(sign_hypothesis(
                        &format!("coefficient_multiplier_{i}_low_lag_{lag}_class_{c}"),
                        &format!("[x^{i}] B_{lag}(n) >= 0 for n = {c} + {p} m (r = 0)"),
                        &format!("[x^{i}] B_{lag}"),
                        bi.clone(),
                        denominator.clone(),
                        c,
                        p,
                        IndexSignCondition::Nonnegative,
                    ));
                }
                if !ai.is_empty() {
                    let upper = vec_poly::add(&bi, &vec_poly::mul(&ai, &previous));
                    hypotheses.push(sign_hypothesis(
                        &format!("coefficient_multiplier_{i}_high_lag_{lag}_class_{c}"),
                        &format!(
                            "[x^{i}] B_{lag}(n) + D(n - {lag}) [x^{}] A_{lag}(n) >= 0 for n = {c} \
                             + {p} m (r = D(n - {lag}))",
                            i + 1
                        ),
                        &format!("[x^{i}] B_{lag} + D(n - {lag}) [x^{}] A_{lag}", i + 1),
                        upper.stripped(),
                        denominator.clone(),
                        c,
                        p,
                        IndexSignCondition::Nonnegative,
                    ));
                }
            }
        }
        if witness {
            hypotheses.push(sign_hypothesis(
                &format!("constant_term_multiplier_positive_class_{c}"),
                &format!(
                    "sum_j [x^0] B_j(n) > 0 for n = {c} + {p} m (so P_n(0) >= sum_j [x^0] B_j(n) \
                     P_(n-j)(0) > 0)"
                ),
                "sum_j [x^0] B_j",
                constant_terms.stripped(),
                denominator,
                c,
                p,
                IndexSignCondition::Positive,
            ));
        }
        hypotheses
    }

    /// The `cancel` route on every class (first-order rules), with the
    /// coefficient-ratio invariant `[x^(D-1)] P_n >= rho lead P_n` when
    /// `ratio` is given (module documentation).
    pub(super) fn cancel_hypotheses(
        &self,
        law: &RecurrenceDegreeLaw,
        classes: &[Class],
        ratio: Option<&Q>,
    ) -> Vec<LiftHypothesis> {
        let mut hypotheses = Vec::new();
        let Some(first_class) = classes.first() else {
            return hypotheses;
        };
        let pairs = self.paired_terms(first_class);
        let (b, a) = pairs.get(&1).cloned().unwrap_or_default();
        let m = x_degree(&b)
            .unwrap_or(0)
            .max(x_degree(&a).map_or(0, |d| d.saturating_sub(1)));
        let shape_ok = m >= 1;
        hypotheses.push(LiftHypothesis {
            name: "cancel_shape".into(),
            statement:
                "P_n = A_n(x) P'_(n-1) + B_n(x) P_(n-1) with deg_x A <= m + 1, deg_x B <= m \
                        for some m >= 1"
                    .into(),
            status: if shape_ok {
                LiftHypothesisStatus::Verified
            } else {
                LiftHypothesisStatus::Failed { witness: None }
            },
            evidence: format!("m = {m} (exact)"),
            sign_decision: None,
        });
        if !shape_ok {
            return hypotheses;
        }
        if ratio.is_some() {
            // The invariant needs a coefficient below the leading one.
            let lowest = classes
                .iter()
                .map(|class| (class.c - 1, law.degree_at(class.c - 1).unwrap()))
                .min_by_key(|(_, degree)| *degree);
            let positive = lowest.is_some_and(|(_, degree)| degree >= 1);
            hypotheses.push(LiftHypothesis {
                name: "cancel_positive_degree".into(),
                statement: "D(n - 1) >= 1 for every n from the induction start on (the ratio \
                            invariant compares [x^(D-1)] P_(n-1) with lead P_(n-1))"
                    .into(),
                status: if positive {
                    LiftHypothesisStatus::Verified
                } else {
                    LiftHypothesisStatus::Failed {
                        witness: lowest.map(|(n, _)| BigInt::from(n)),
                    }
                },
                evidence: format!(
                    "the law is nondecreasing; smallest D(n - 1) on the first classes: {:?} \
                     (exact)",
                    lowest.map(|(_, degree)| degree)
                ),
                sign_decision: None,
            });
        }
        for class in classes {
            let c = class.c;
            let p = law.period;
            let denominator = self.left_factor_on(class, law);
            let on = |poly: Vec<Q>| on_class(&poly, c, p);
            let previous = vec![
                integer(law.degree_at(c - 1).unwrap() as i64),
                integer(law.growth_per_period as i64),
            ];
            // E_i(D(n-1) - shift) = B_i + (D(n-1) - shift) A_(i+1), with
            // B_(-1) = 0: the coefficient of P_(n-1)[x^(D - shift)] in
            // [x^(D - shift + i)] P_n (with the factor 1 / f(n)).
            let level = |i: i64, shift: i64| -> Vec<Q> {
                let r = vec_poly::add(&previous, &[integer(-shift)]).stripped();
                let b_part = if i >= 0 {
                    on(column(&b, i as usize))
                } else {
                    Vec::new()
                };
                vec_poly::add(
                    &b_part,
                    &vec_poly::mul(&on(column(&a, (i + 1) as usize)), &r),
                )
                .stripped()
            };
            let a_top = on(column(&a, m + 1));
            let m_signed = m as i64;
            let top = level(m_signed, 0);
            let next = level(m_signed - 1, 0);
            let step = class.gaps[&1];
            let sign =
                |name: &str, statement: String, function: &str, numerator: Vec<Q>, condition| {
                    sign_hypothesis(
                        &format!("{name}_class_{c}"),
                        &format!("{statement} for n = {c} + {p} m'"),
                        function,
                        numerator,
                        denominator.clone(),
                        c,
                        p,
                        condition,
                    )
                };
            // `alpha p_D + beta p_(D-1) >= 0` on the cone `p_D > 0`,
            // `p_(D-1) >= rho p_D`: `beta >= 0` and `alpha + rho beta >= 0`.
            let invariant = |alpha: Vec<Q>, beta: Vec<Q>, rho: &Q, text: &str| {
                vec![
                    sign(
                        "ratio_invariant_slope",
                        format!(
                            "the coefficient of [x^(D-1)] P_(n-1) in rho^-1 [x^(D'-1)] P_n - \
                             lead P_n is >= 0 ({text})"
                        ),
                        "beta",
                        beta.clone(),
                        IndexSignCondition::Nonnegative,
                    ),
                    sign(
                        "ratio_invariant",
                        format!(
                            "alpha + rho beta >= 0 with rho = {}, so [x^(D'-1)] P_n >= rho lead \
                             P_n from [x^(D-1)] P_(n-1) >= rho lead P_(n-1) ({text})",
                            format_rational_coeff(rho)
                        ),
                        "alpha + rho beta",
                        vec_poly::add(&alpha, &vec_poly::scale(&beta, rho)).stripped(),
                        IndexSignCondition::Nonnegative,
                    ),
                ]
            };
            if step == m_signed {
                hypotheses.push(sign(
                    "top_grows",
                    format!("A_(m+1)(n) D(n-1) + B_m(n) > 0 (step m = {m})"),
                    "A_(m+1) D + B_m",
                    top.clone(),
                    IndexSignCondition::Positive,
                ));
                if let Some(rho) = ratio {
                    // lead P_n = E_m(D) p_D, [x^(D'-1)] P_n = E_(m-1)(D) p_D +
                    // E_m(D-1) p_(D-1).
                    let alpha = vec_poly::sub(&next, &vec_poly::scale(&top, rho));
                    hypotheses.extend(invariant(
                        alpha,
                        level(m_signed, 1),
                        rho,
                        "the top terms grow",
                    ));
                }
            } else if step + 1 == m_signed {
                hypotheses.push(sign(
                    "top_cancels",
                    "A_(m+1)(n) D(n-1) + B_m(n) = 0 (step m - 1)".into(),
                    "A_(m+1) D + B_m",
                    top,
                    IndexSignCondition::Zero,
                ));
                hypotheses.push(sign(
                    "top_derivative_nonpositive",
                    "A_(m+1)(n) <= 0".into(),
                    "A_(m+1)",
                    a_top.clone(),
                    IndexSignCondition::Nonpositive,
                ));
                match ratio {
                    None => hypotheses.push(sign(
                        "next_level_positive",
                        "A_m(n) D(n-1) + B_(m-1)(n) > 0".into(),
                        "A_m D + B_(m-1)",
                        next,
                        IndexSignCondition::Positive,
                    )),
                    Some(rho) => {
                        // lead P_n = E_(m-1)(D) p_D - A_(m+1) p_(D-1) > 0 on
                        // the cone; [x^(D'-1)] P_n = E_(m-2)(D) p_D +
                        // E_(m-1)(D-1) p_(D-1) - 2 A_(m+1) p_(D-2), the last
                        // term >= 0.
                        hypotheses.push(sign(
                            "next_level_with_ratio_positive",
                            format!(
                                "A_m(n) D(n-1) + B_(m-1)(n) - rho A_(m+1)(n) > 0 with rho = {}",
                                format_rational_coeff(rho)
                            ),
                            "A_m D + B_(m-1) - rho A_(m+1)",
                            vec_poly::sub(&next, &vec_poly::scale(&a_top, rho)).stripped(),
                            IndexSignCondition::Positive,
                        ));
                        let alpha =
                            vec_poly::sub(&level(m_signed - 2, 0), &vec_poly::scale(&next, rho));
                        let beta =
                            vec_poly::add(&level(m_signed - 1, 1), &vec_poly::scale(&a_top, rho));
                        hypotheses.extend(invariant(alpha, beta, rho, "the top terms cancel"));
                    }
                }
            } else {
                hypotheses.push(LiftHypothesis {
                    name: format!("cancel_step_class_{c}"),
                    statement: format!("the degree step on the class of {c} is m or m - 1"),
                    status: LiftHypothesisStatus::Failed {
                        witness: Some(BigInt::from(c)),
                    },
                    evidence: format!("the step is {step}, m = {m} (exact)"),
                    sign_decision: None,
                });
            }
            hypotheses.extend(self.coefficient_hypotheses(law, class, false));
        }
        hypotheses
    }
}
