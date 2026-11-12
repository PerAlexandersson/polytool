//! Differential and scalar basis-map preservation criteria.

use super::*;

// ---------------------------------------------------------------------------
// First-order differential families
// ---------------------------------------------------------------------------

pub(super) fn differential(
    map: &BasisMap,
    g: &[Q],
    multiplier: &[Q],
    v: &[Q],
    _options: &BasisMapTheoremOptions,
) -> Result<Analysis, RecurrenceError> {
    let _ = map;
    let a = vec_poly::coeff(multiplier, 0);
    let b = vec_poly::coeff(multiplier, 1);
    let v = vec_poly::trim(v);
    let degree = vec_poly::degree(&v).unwrap_or(0);
    let g_text = poly_text(g);
    let positive_leading = leading(g).is_positive();
    match degree {
        0 => {
            let kappa = -vec_poly::coeff(&v, 0);
            let mut analysis = Analysis {
                shape: BasisMapShape::HermiteDifferential {
                    initial_row: g_text,
                    a: text(&a),
                    b: text(&b),
                    kappa: text(&kappa),
                },
                status: BasisMapTheoremStatus::Proved,
                theorem: "Hermite-type operators L f = (a + b t) f - kappa f', b > 0, kappa >= \
                          0: every shift L - rho keeps real-rootedness (the derivative \
                          coefficient is <= 0 everywhere: Ma and Wang, EJC 15 (2008) R17, \
                          Theorem 2), and T(p) = c prod_j (L - rho_j) g for p = c prod_j (u - \
                          rho_j); so every real-rooted input has a real-rooted image"
                    .into(),
                hypotheses: Vec::new(),
                conclusions: Vec::new(),
                counterexample: None,
                lean: Vec::new(),
            };
            let checks = [
                (
                    "multiplier_slope_positive",
                    "b > 0".to_string(),
                    b.is_positive(),
                ),
                (
                    "kappa_positive",
                    "kappa > 0".to_string(),
                    kappa.is_positive(),
                ),
            ];
            for (name, statement, holds) in checks {
                analysis
                    .hypotheses
                    .push(hypothesis(name, statement, holds, "exact".into()));
            }
            let real = !vec_poly::is_zero(g) && is_real_rooted(g);
            analysis.hypotheses.push(hypothesis(
                "initial_row_real_rooted",
                "g = P_0 is nonzero and real-rooted".into(),
                real,
                "exact root isolation".into(),
            ));
            if !vec_poly::is_zero(g) && !real {
                analysis.status = BasisMapTheoremStatus::NotPreserving;
                analysis.counterexample = Some(counterexample(
                    &[Q::one()],
                    g,
                    "T(1) = g is not real-rooted".into(),
                ));
                return Ok(analysis);
            }
            if analysis
                .hypotheses
                .iter()
                .all(|h| h.status == LiftHypothesisStatus::Verified)
            {
                analysis.push_conclusion(BasisMapConclusion {
                    input_interval: AlgebraicInterval::whole_line(),
                    output_interval: AlgebraicInterval::whole_line(),
                    nonnegative_coefficients: false,
                    sharp: true,
                    statement: "every real-rooted input has a real-rooted image".into(),
                    hypotheses: Vec::new(),
                });
            }
            Ok(analysis)
        }
        1 => {
            let kappa = vec_poly::coeff(&v, 1);
            let h = -vec_poly::coeff(&v, 0) / &kappa;
            let mut analysis = Analysis {
                shape: BasisMapShape::RayDifferential {
                    initial_row: g_text,
                    a: text(&a),
                    b: text(&b),
                    kappa: text(&kappa),
                    root: text(&h),
                },
                status: BasisMapTheoremStatus::Proved,
                theorem: "ray basis maps: for L f = (a + b t) f + kappa (t - h) f', b > 0, kappa \
                          >= 0, and g real-rooted with every root <= h, of multiplicity m at h, \
                          T(u^k) = L^k g sends every real polynomial with all roots in S to one \
                          with all roots in (-inf, h] if and only if S lies in (-inf, a + b h + \
                          kappa m] (RealRooted rayBasisMap_maps_realRootedIn_iff; the step is \
                          Ma and Wang, EJC 15 (2008) R17, Theorem 2)"
                    .into(),
                hypotheses: Vec::new(),
                conclusions: Vec::new(),
                counterexample: None,
                lean: vec!["RealRooted.rayBasisMap_maps_realRootedIn_iff"],
            };
            analysis.hypotheses.push(hypothesis(
                "multiplier_slope_positive",
                "b > 0".into(),
                b.is_positive(),
                "exact".into(),
            ));
            analysis.hypotheses.push(hypothesis(
                "kappa_nonnegative",
                "kappa >= 0".into(),
                !kappa.is_negative(),
                "exact".into(),
            ));
            let nonzero = !vec_poly::is_zero(g);
            let rooted =
                nonzero && AlgebraicInterval::at_most(rational(&h)).contains_all_roots_of(g);
            analysis.hypotheses.push(hypothesis(
                "initial_row_roots",
                format!(
                    "g = P_0 is nonzero and real-rooted with every root <= {}",
                    text(&h)
                ),
                rooted,
                "exact root isolation".into(),
            ));
            if nonzero && !is_real_rooted(g) {
                analysis.status = BasisMapTheoremStatus::NotPreserving;
                analysis.counterexample = Some(counterexample(
                    &[Q::one()],
                    g,
                    "T(1) = g is not real-rooted".into(),
                ));
                return Ok(analysis);
            }
            if analysis
                .hypotheses
                .iter()
                .all(|x| x.status == LiftHypothesisStatus::Verified)
            {
                let m = multiplicity(g, &h);
                let upper = &a + &b * &h + &kappa * integer(m as i64);
                let output = AlgebraicInterval::at_most(rational(&h));
                analysis.push_conclusion(BasisMapConclusion {
                    input_interval: AlgebraicInterval::at_most(rational(&upper)),
                    nonnegative_coefficients: positive_leading && !h.is_positive(),
                    output_interval: output,
                    sharp: true,
                    statement: format!(
                        "source (-inf, a + b h + kappa m] = (-inf, {}] with m = {m}, target \
                         (-inf, {}]",
                        text(&upper),
                        text(&h)
                    ),
                    hypotheses: Vec::new(),
                });
            }
            Ok(analysis)
        }
        2 => {
            let kappa = vec_poly::coeff(&v, 2);
            let p = vec_poly::coeff(&v, 1) / &kappa;
            let q0 = vec_poly::coeff(&v, 0) / &kappa;
            // Roots of t^2 + p t + q0.
            let discriminant = &p * &p - integer(4) * &q0;
            let root = QuadraticNumber::square_root(discriminant.clone()).ok();
            let rational_roots = discriminant.is_positive()
                && root.as_ref().is_some_and(QuadraticNumber::is_rational);
            if !rational_roots {
                return Ok(Analysis::not_applicable(
                    "the derivative coefficient is quadratic without two distinct rational \
                     roots"
                        .into(),
                ));
            }
            let s = root
                .expect("checked")
                .as_rational()
                .cloned()
                .expect("rational");
            let half = Q::new(BigInt::one(), BigInt::from(2));
            let lower = (-&p - &s) * &half;
            let upper = (-&p + &s) * &half;
            let mut analysis = Analysis {
                shape: BasisMapShape::IntervalDifferential {
                    initial_row: g_text,
                    a: text(&a),
                    b: text(&b),
                    kappa: text(&kappa),
                    lower_root: text(&lower),
                    upper_root: text(&upper),
                },
                status: BasisMapTheoremStatus::Proved,
                theorem: "interval basis maps: for L f = (a + b t) f + kappa (t - l)(t - h) f', \
                          l < h, kappa > 0, b + kappa deg g > 0, and g with positive leading \
                          coefficient, real-rooted in [l, h] with endpoint multiplicities m_l, \
                          m_h, T(u^k) = L^k g sends every real polynomial with all roots in S \
                          to one with all roots in [l, h] if and only if S lies in [a + b l - \
                          kappa (h - l) m_l, a + b h + kappa (h - l) m_h] (RealRooted \
                          intervalBasisMap_maps_realRootedIn_iff; Ma and Wang, EJC 15 (2008) \
                          R17, Theorem 2 and Corollary 3); a negative leading coefficient of g \
                          only changes the sign of T"
                    .into(),
                hypotheses: Vec::new(),
                conclusions: Vec::new(),
                counterexample: None,
                lean: vec!["RealRooted.intervalBasisMap_maps_realRootedIn_iff"],
            };
            let deg_g = vec_poly::degree(g).unwrap_or(0);
            let top = &b + &kappa * integer(deg_g as i64);
            analysis.hypotheses.push(hypothesis(
                "kappa_positive",
                "kappa > 0".into(),
                kappa.is_positive(),
                "exact".into(),
            ));
            analysis.hypotheses.push(hypothesis(
                "top_coefficient_positive",
                format!("b + kappa deg g = {} > 0", text(&top)),
                top.is_positive(),
                "exact".into(),
            ));
            let interval = AlgebraicInterval::closed(rational(&lower), rational(&upper));
            let nonzero = !vec_poly::is_zero(g);
            let rooted = nonzero && interval.contains_all_roots_of(g);
            analysis.hypotheses.push(hypothesis(
                "initial_row_roots",
                format!(
                    "g = P_0 is nonzero and real-rooted with every root in [{}, {}]",
                    text(&lower),
                    text(&upper)
                ),
                rooted,
                "exact root isolation".into(),
            ));
            if nonzero && !is_real_rooted(g) {
                analysis.status = BasisMapTheoremStatus::NotPreserving;
                analysis.counterexample = Some(counterexample(
                    &[Q::one()],
                    g,
                    "T(1) = g is not real-rooted".into(),
                ));
                return Ok(analysis);
            }
            if analysis
                .hypotheses
                .iter()
                .all(|x| x.status == LiftHypothesisStatus::Verified)
            {
                let m_lower = integer(multiplicity(g, &lower) as i64);
                let m_upper = integer(multiplicity(g, &upper) as i64);
                let width = &upper - &lower;
                let source_lower = &a + &b * &lower - &kappa * &width * &m_lower;
                let source_upper = &a + &b * &upper + &kappa * &width * &m_upper;
                analysis.push_conclusion(BasisMapConclusion {
                    input_interval: AlgebraicInterval::closed(
                        rational(&source_lower),
                        rational(&source_upper),
                    ),
                    output_interval: interval,
                    nonnegative_coefficients: positive_leading && !upper.is_positive(),
                    sharp: true,
                    statement: format!(
                        "source [{}, {}], target [{}, {}]",
                        text(&source_lower),
                        text(&source_upper),
                        text(&lower),
                        text(&upper)
                    ),
                    hypotheses: Vec::new(),
                });
            }
            Ok(analysis)
        }
        _ => Ok(Analysis::not_applicable(
            "the derivative coefficient has degree above two".into(),
        )),
    }
}

// ---------------------------------------------------------------------------
// Index-dependent scalar term
// ---------------------------------------------------------------------------

pub(super) fn scalar_index(
    map: &BasisMap,
    options: &BasisMapTheoremOptions,
) -> Result<Analysis, RecurrenceError> {
    let rule = map.rule();
    let not_applicable = |reason: &str| Ok(Analysis::not_applicable(reason.to_string()));
    if rule.max_offset() != 1 || map.seeds().len() != 1 {
        return not_applicable("an index-dependent rule other than first order without a prefix");
    }
    let leader = match &rule.denominator {
        None => Q::one(),
        Some(f) => match Biv::from_bivar(f).as_constant() {
            Some(c) => c,
            None => return not_applicable("an index-dependent leader"),
        },
    };
    let mut multiplier = Biv::zero();
    let mut v = Biv::zero();
    for term in &rule.terms {
        let c = Biv::from_bivar(&term.coeff).scale(&leader.recip());
        match term.deriv_order {
            0 => multiplier = multiplier.add(&c),
            1 => v = v.add(&c),
            _ => return not_applicable("an index-dependent rule with higher derivatives"),
        }
    }
    if !v.is_index_free() {
        return not_applicable("an index-dependent derivative coefficient");
    }
    // Only the t^0 coefficient of the multiplier may depend on n.
    let rows = multiplier.rows();
    if rows.iter().skip(1).any(|row| row.len() > 1) {
        return not_applicable("an index-dependent coefficient of t");
    }
    let alpha_n: Vec<Q> = vec_poly::trim(
        &rows
            .iter()
            .map(|row| row.first().cloned().unwrap_or_else(Q::zero))
            .collect::<Vec<_>>(),
    );
    let t_part = t_polynomial(&multiplier);
    let b = vec_poly::coeff(&t_part, 1);
    if vec_poly::degree(&t_part).unwrap_or(0) > 1 {
        return not_applicable("an index-dependent rule whose multiplier has degree above one");
    }
    let v = t_polynomial(&v);
    if vec_poly::degree(&v).unwrap_or(0) > 1 {
        return not_applicable("an index-dependent rule with a nonlinear derivative coefficient");
    }
    // alpha as a polynomial in the basis index k: n = i0 + 1 + k.
    let shift = integer(map.first_index() as i64 + 1);
    let alpha_k = vec_poly::compose(&alpha_n, &[shift, Q::one()]);
    let kappa = vec_poly::coeff(&v, 1);
    let a = vec_poly::coeff(&alpha_k, 0);
    // Without a derivative term the ray root h is free: the smallest h >= 0
    // with a + b h >= 0.
    let h = if !kappa.is_zero() {
        -vec_poly::coeff(&v, 0) / &kappa
    } else if a.is_negative() && b.is_positive() {
        -(&a / &b)
    } else {
        Q::zero()
    };
    let g = &map.seeds()[0];
    let slope = vec_poly::coeff(&alpha_k, 1);
    let affine = vec_poly::degree(&alpha_k).unwrap_or(0) <= 1;
    let mut analysis = Analysis {
        shape: BasisMapShape::ScalarIndex {
            initial_row: poly_text(g),
            scalar: poly_text(&alpha_k),
            b: text(&b),
            kappa: text(&kappa),
            root: text(&h),
        },
        status: BasisMapTheoremStatus::Proved,
        theorem: "index-dependent scalar term: P_(k+1) = (L + k s) P_k with L f = (a + b t) f + \
                  kappa (t - h) f' gives T = U_g o R_s (RealRooted \
                  basisMap_factorization_of_scalarIndex_recurrence), where R_s(u^k) = u (u + s) \
                  ... (u + (k-1) s) keeps every root <= 0 for s >= 0 (Su, Yang and Zhang, \
                  Australas. J. Combin. 56 (2013), Theorem 4; Brenti, Mem. AMS 413 (1989), \
                  Theorem 2.4.3; RealRooted NewtonParametersAffine.preservesNonposRooted) and \
                  U_g(u^k) = L^k g sends (-inf, a + b h + kappa m] into (-inf, h] (RealRooted \
                  rayBasisMap_maps_realRootedIn_iff)"
            .into(),
        hypotheses: Vec::new(),
        conclusions: Vec::new(),
        counterexample: None,
        lean: vec![
            "RealRooted.basisMap_factorization_of_scalarIndex_recurrence",
            "RealRooted.NewtonParametersAffine.preservesNonposRooted",
            "RealRooted.rayBasisMap_maps_realRootedIn_iff",
        ],
    };
    let pure_newton = kappa.is_zero() && is_constant(g) && !vec_poly::is_zero(g) && b.is_positive();
    let slope_ok = affine && !slope.is_negative();
    analysis.hypotheses.push(hypothesis(
        "scalar_affine_nondecreasing",
        format!(
            "alpha(k) = {} is affine in k with slope >= 0",
            crate::recurrence::implicit::algebra::poly_text(&alpha_k, "k")
        ),
        slope_ok,
        "exact".into(),
    ));
    if !slope_ok && pure_newton {
        // Newton basis prod_(j < k) (b t + alpha(j)): iff theorem.
        analysis.status = BasisMapTheoremStatus::NotPreserving;
        analysis.theorem = "Newton bases P_k = prod_(j < k) (b t + c_j), b > 0: T preserves \
                            real-rootedness of the inputs with every root <= 0 if and only if \
                            c_j is an arithmetic progression with nonnegative step (RealRooted \
                            newtonPreservesNonposRooted_iff_parametersAffine; Su, Yang and Zhang, \
                            Theorem 4)"
            .into();
        analysis.lean = vec!["RealRooted.newtonPreservesNonposRooted_iff_parametersAffine"];
        let limit = options.witness_search_limit.clamp(1, 12);
        'search: for r in [1i64, 2, 4, 8, 16] {
            for n in 2..=limit {
                let input = vec_poly::pow(&[integer(r), Q::one()], n);
                let image = map.apply(&input)?;
                if !is_real_rooted(&image) {
                    analysis.counterexample = Some(counterexample(
                        &input,
                        &image,
                        format!("T((u + {r})^{n}) is not real-rooted (exact)"),
                    ));
                    break 'search;
                }
            }
        }
        return Ok(analysis);
    }
    analysis.hypotheses.push(hypothesis(
        "multiplier_slope_positive",
        "b > 0".into(),
        b.is_positive(),
        "exact".into(),
    ));
    analysis.hypotheses.push(hypothesis(
        "kappa_nonnegative",
        "kappa >= 0".into(),
        !kappa.is_negative(),
        "exact".into(),
    ));
    let nonzero = !vec_poly::is_zero(g);
    let rooted = nonzero && AlgebraicInterval::at_most(rational(&h)).contains_all_roots_of(g);
    analysis.hypotheses.push(hypothesis(
        "initial_row_roots",
        format!(
            "g = P_0 is nonzero and real-rooted with every root <= {}",
            text(&h)
        ),
        rooted,
        "exact root isolation".into(),
    ));
    let m = if rooted { multiplicity(g, &h) } else { 0 };
    let source_end = &a + &b * &h + &kappa * integer(m as i64);
    analysis.hypotheses.push(hypothesis(
        "ray_source_contains_nonpositive_ray",
        format!(
            "a + b h + kappa m = {} >= 0 (the source of U_g contains (-inf, 0])",
            text(&source_end)
        ),
        !source_end.is_negative(),
        "exact".into(),
    ));
    if analysis
        .hypotheses
        .iter()
        .all(|x| x.status == LiftHypothesisStatus::Verified)
    {
        analysis.push_conclusion(BasisMapConclusion {
            input_interval: AlgebraicInterval::nonpositive_ray(),
            output_interval: AlgebraicInterval::at_most(rational(&h)),
            nonnegative_coefficients: leading(g).is_positive() && !h.is_positive(),
            sharp: false,
            statement: format!(
                "(-inf, 0] -> (-inf, {}] through R_s, s = {}, and U_g",
                text(&h),
                text(&slope)
            ),
            hypotheses: Vec::new(),
        });
    }
    Ok(analysis)
}

// ---------------------------------------------------------------------------
// Order two
// ---------------------------------------------------------------------------

pub(super) fn order_two(
    map: &BasisMap,
    current: &[Q],
    lag: &[Q],
    options: &BasisMapTheoremOptions,
) -> Result<Analysis, RecurrenceError> {
    let seeds = map.seeds();
    let p0 = &seeds[0];
    let p1 = &seeds[1];
    if !is_constant(p0) || vec_poly::is_zero(p0) || vec_poly::degree(p1).unwrap_or(0) > 1 {
        return Ok(Analysis::not_applicable(
            "an order-two rule whose first row is not a nonzero constant or whose second row \
             has degree above one"
                .into(),
        ));
    }
    let normalization = p0[0].clone();
    let l = vec_poly::coeff(p1, 0) / &normalization;
    let m = vec_poly::coeff(p1, 1) / &normalization;
    let a = vec_poly::coeff(current, 0);
    let d = vec_poly::coeff(current, 1);
    let c = vec_poly::coeff(lag, 0);
    let b = vec_poly::coeff(lag, 1);
    if d.is_zero() {
        if !m.is_zero() {
            return Ok(Analysis::not_applicable(
                "a constant current multiplier with a nonconstant second row".into(),
            ));
        }
        return Ok(constant_current(&normalization, &a, &b, &c, &l, options));
    }
    if d.is_negative() {
        return Ok(Analysis::not_applicable(
            "an affine current multiplier with negative slope".into(),
        ));
    }
    Ok(affine_current(&normalization, &a, &d, &b, &c, &l, &m))
}

/// The scalar criterion `v_n >= 0` for `v_0 = 1`, `v_1 = s`,
/// `v_(n+2) = a v_(n+1) + c v_n` with `a, s >= 0` (RealRooted
/// `orderTwoConstants_nonneg_iff`): `a^2 + 4c >= 0` and
/// (`2s >= a` or `c >= s^2 - a s`).
fn scalar_criterion(a: &Q, c: &Q, s: &Q) -> (bool, String) {
    let disc = a * a + integer(4) * c;
    let first = !disc.is_negative();
    let second = integer(2) * s >= *a || *c >= s * s - a * s;
    (
        first && second,
        format!(
            "a^2 + 4c = {} >= 0: {first}; 2s >= a or c >= s^2 - a s: {second}",
            text(&disc)
        ),
    )
}

fn constant_current(
    normalization: &Q,
    a: &Q,
    b: &Q,
    c: &Q,
    s: &Q,
    options: &BasisMapTheoremOptions,
) -> Analysis {
    let mut analysis = Analysis {
        shape: BasisMapShape::ConstantCurrent {
            normalization: text(normalization),
            a: text(a),
            b: text(b),
            c: text(c),
            initial_constant: text(s),
        },
        status: BasisMapTheoremStatus::Proved,
        theorem: "constant current multiplier: P_0 = 1, P_1 = s, P_(k+2) = a P_(k+1) + (b t + c) \
                  P_k with a, s >= 0, b > 0.  P_k(t) = Q_k(b t + c) for the rows Q_k with b = 1, \
                  c = 0, whose map sends PF inputs to PF outputs (an interlacing-preserving 2 x \
                  2 matrix, Branden, arXiv:1410.6601, Theorem 8.5), so every PF input has a \
                  real-rooted image with roots <= -c/b; PF inputs (also those with roots in \
                  any [-R, 0]) have PF images if and only if every P_k(0) >= 0, if and only if \
                  a^2 + 4c >= 0 and (2s >= a or c >= s^2 - a s) (RealRooted \
                  orderTwoBasis_pf_equivalences, basisTransform_orderTwoBasis_preserves_bounded_\
                  pf_iff)"
            .into(),
        hypotheses: Vec::new(),
        conclusions: Vec::new(),
        counterexample: None,
        lean: vec![
            "RealRooted.basisTransform_orderTwoBasis_eq_comp",
            "RealRooted.orderTwoZeroTransform_isPFPolynomial",
            "RealRooted.orderTwoBasis_pf_equivalences",
            "RealRooted.basisTransform_orderTwoBasis_preserves_bounded_pf_iff",
        ],
    };
    for (name, statement, holds) in [
        ("current_nonnegative", "a >= 0", !a.is_negative()),
        ("initial_constant_nonnegative", "s >= 0", !s.is_negative()),
        ("lag_slope_positive", "b > 0", b.is_positive()),
    ] {
        analysis
            .hypotheses
            .push(hypothesis(name, statement.into(), holds, "exact".into()));
    }
    if analysis
        .hypotheses
        .iter()
        .any(|h| h.status != LiftHypothesisStatus::Verified)
    {
        return analysis;
    }
    let root_bound = -(c / b);
    analysis.push_conclusion(BasisMapConclusion {
        input_interval: AlgebraicInterval::nonpositive_ray(),
        output_interval: AlgebraicInterval::at_most(rational(&root_bound)),
        nonnegative_coefficients: false,
        sharp: false,
        statement: format!(
            "every PF input has a real-rooted image with every root <= -c/b = {}",
            text(&root_bound)
        ),
        hypotheses: Vec::new(),
    });
    let (criterion, evidence) = scalar_criterion(a, c, s);
    analysis.hypotheses.push(hypothesis(
        "constant_terms_nonnegative",
        "a^2 + 4c >= 0 and (2s >= a or c >= s^2 - a s), equivalently every P_k(0) >= 0".into(),
        criterion,
        evidence,
    ));
    if criterion {
        let bound = if root_bound.is_negative() {
            root_bound
        } else {
            Q::zero()
        };
        analysis.push_conclusion(BasisMapConclusion {
            input_interval: AlgebraicInterval::nonpositive_ray(),
            output_interval: AlgebraicInterval::at_most(rational(&bound)),
            nonnegative_coefficients: normalization.is_positive(),
            sharp: true,
            statement: "every PF input has a PF image (zero allowed)".into(),
            hypotheses: Vec::new(),
        });
    } else {
        // The first monomial with a negative constant term.
        let (mut previous, mut current) = (Q::one(), s.clone());
        for n in 2..=options.witness_search_limit {
            let next = a * &current + c * &previous;
            previous = current;
            current = next;
            if current.is_negative() {
                analysis.counterexample = Some(BasisMapCounterexample {
                    input: {
                        let mut monomial = vec!["0".to_string(); n];
                        monomial.push("1".into());
                        monomial
                    },
                    image: Vec::new(),
                    failure: format!(
                        "P_{n}(0) / P_0 = {} < 0 while P_{n} has a positive leading \
                         coefficient: the PF input u^{n} has an image with a positive root (the \
                         image is omitted)",
                        text(&current)
                    ),
                });
                break;
            }
        }
    }
    analysis
}

fn affine_current(normalization: &Q, a: &Q, d: &Q, b: &Q, c: &Q, l: &Q, m: &Q) -> Analysis {
    let beta = b / d;
    let k = &beta * &beta + a * &beta - c;
    let eta = m / d;
    let s = l - &eta * a;
    let mut analysis = Analysis {
        shape: BasisMapShape::AffineCurrent {
            normalization: text(normalization),
            a: text(a),
            d: text(d),
            b: text(b),
            c: text(c),
            initial_constant: text(l),
            initial_slope: text(m),
        },
        status: BasisMapTheoremStatus::Proved,
        theorem: "affine current multiplier with an affine initial row (module documentation, \
                  \"The affine-current theorem\"): with beta = b/d, K = beta^2 + a beta - c > 0 \
                  and eta = m/d >= 0, every real polynomial with all roots in E = [-beta - \
                  sqrt(K), -beta + sqrt(K)] has a real-rooted image, and no interval of positive \
                  length outside E works (Theorem 1); for eta > 0 the images lie in I = [(-a - \
                  2 beta - 2 sqrt(K))/d, (-a - 2 beta + 2 sqrt(K))/d] if and only if eta >= 1/2 \
                  and (s - (2 eta - 1) beta)^2 <= (2 eta - 1)^2 K, s = l - eta a (Theorem 2); \
                  with a >= 0, K >= beta^2, a + 2 beta >= 0, a^2 + 4c >= 0, l >= 0 and the \
                  scalar criterion, PF inputs with roots in [-R, 0], R = beta + sqrt(K), have \
                  PF images, sharply (the common-interlacer proof; the Chebyshev case is \
                  Iserles and Saff, Math. Proc. Cambridge Philos. Soc. 105 (1989), Theorem 1; \
                  the scalar criterion is RealRooted orderTwoConstants_nonneg_iff)"
            .into(),
        hypotheses: Vec::new(),
        conclusions: Vec::new(),
        counterexample: None,
        lean: vec!["RealRooted.orderTwoConstants_nonneg_iff"],
    };
    analysis.hypotheses.push(hypothesis(
        "characteristic_parameter_positive",
        format!("K = beta^2 + a beta - c = {} > 0", text(&k)),
        k.is_positive(),
        "exact".into(),
    ));
    analysis.hypotheses.push(hypothesis(
        "initial_slope_nonnegative",
        format!("eta = m/d = {} >= 0", text(&eta)),
        !eta.is_negative(),
        "exact".into(),
    ));
    if analysis
        .hypotheses
        .iter()
        .any(|h| h.status != LiftHypothesisStatus::Verified)
    {
        return analysis;
    }
    let root_k = QuadraticNumber::square_root(k.clone()).expect("K > 0");
    let source = AlgebraicInterval::closed(
        root_k.negated().add_rational(&-beta.clone()),
        root_k.add_rational(&-beta.clone()),
    );
    analysis.push_conclusion(BasisMapConclusion {
        input_interval: source.clone(),
        output_interval: AlgebraicInterval::whole_line(),
        nonnegative_coefficients: false,
        sharp: false,
        statement: format!(
            "every real polynomial with all roots in E = {source} has a real-rooted image"
        ),
        hypotheses: Vec::new(),
    });
    // Target interval I in t: ((-a - 2 beta) -+ 2 sqrt(K))/d.
    let center = (-a - integer(2) * &beta) / d;
    let half_width = root_k.scale(&(integer(2) / d));
    let target = AlgebraicInterval::closed(
        half_width.negated().add_rational(&center),
        half_width.add_rational(&center),
    );
    let mut target_proved = false;
    if eta.is_positive() {
        let two_eta = integer(2) * &eta - Q::one();
        let deviation = &s - &two_eta * &beta;
        let holds = !two_eta.is_negative() && &deviation * &deviation <= &two_eta * &two_eta * &k;
        analysis.hypotheses.push(hypothesis(
            "target_interval",
            format!(
                "2 eta - 1 = {} >= 0 and (s - (2 eta - 1) beta)^2 <= (2 eta - 1)^2 K, s = {}",
                text(&two_eta),
                text(&s)
            ),
            holds,
            "exact".into(),
        ));
        if holds {
            target_proved = true;
            analysis.push_conclusion(BasisMapConclusion {
                input_interval: source.clone(),
                output_interval: target.clone(),
                nonnegative_coefficients: false,
                sharp: false,
                statement: format!(
                    "every real polynomial with all roots in E = {source} has an image with all \
                     roots in I = {target}"
                ),
                hypotheses: Vec::new(),
            });
        }
    }
    // PF on [-R, 0].
    let radius = root_k.add_rational(&beta);
    let checks: Vec<(&str, String, bool)> = vec![
        (
            "current_constant_nonnegative",
            format!("a = {} >= 0", text(a)),
            !a.is_negative(),
        ),
        (
            "zero_in_source",
            format!("K >= beta^2 ({} >= {})", text(&k), text(&(&beta * &beta))),
            k >= &beta * &beta,
        ),
        (
            "radius_positive",
            format!("R = beta + sqrt(K) = {radius} > 0"),
            radius.signum() == Ordering::Greater,
        ),
        (
            "target_upper_end_nonpositive",
            format!(
                "a + 2 beta = {} >= 0 and a^2 + 4c = {} >= 0",
                text(&(a + integer(2) * &beta)),
                text(&(a * a + integer(4) * c))
            ),
            !(a + integer(2) * &beta).is_negative() && !(a * a + integer(4) * c).is_negative(),
        ),
        (
            "initial_constant_nonnegative",
            format!("l = {} >= 0", text(l)),
            !l.is_negative(),
        ),
    ];
    for (name, statement, holds) in checks {
        analysis
            .hypotheses
            .push(hypothesis(name, statement, holds, "exact".into()));
    }
    let (criterion, evidence) = scalar_criterion(a, c, l);
    analysis.hypotheses.push(hypothesis(
        "constant_terms_nonnegative",
        "2l >= a or c >= l^2 - a l (with a^2 + 4c >= 0): every P_k(0) >= 0".into(),
        criterion,
        evidence,
    ));
    let pf_names = [
        "current_constant_nonnegative",
        "zero_in_source",
        "radius_positive",
        "target_upper_end_nonpositive",
        "initial_constant_nonnegative",
        "constant_terms_nonnegative",
    ];
    let pf = analysis
        .hypotheses
        .iter()
        .filter(|h| pf_names.contains(&h.name.as_str()))
        .all(|h| h.status == LiftHypothesisStatus::Verified);
    if pf {
        let input =
            AlgebraicInterval::closed(radius.negated(), QuadraticNumber::rational(Q::zero()));
        let output = if target_proved {
            target.intersection(&AlgebraicInterval::nonpositive_ray())
        } else {
            AlgebraicInterval::nonpositive_ray()
        };
        analysis.push_conclusion(BasisMapConclusion {
            input_interval: input.clone(),
            output_interval: output,
            nonnegative_coefficients: normalization.is_positive(),
            sharp: false,
            statement: format!(
                "every PF input with all roots in {input} has a PF image; no larger [-R, 0] works"
            ),
            hypotheses: Vec::new(),
        });
    }
    analysis
}
