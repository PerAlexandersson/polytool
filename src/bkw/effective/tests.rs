//! Regression examples, negative controls and independent replay mutations.

use super::*;

fn request(n: &[&[i64]], d: &[&[i64]]) -> EffectiveBkwRequest {
    let encode = |rows: &[&[i64]]| {
        rows.iter()
            .map(|r| r.iter().map(ToString::to_string).collect())
            .collect()
    };
    EffectiveBkwRequest {
        generating_function: RationalGeneratingFunction {
            schema: crate::recurrence::RATIONAL_GENERATING_FUNCTION_SCHEMA.into(),
            first_index: 0,
            numerator: encode(n),
            denominator: encode(d),
        },
        interval: None,
        target_zero_count: None,
        uniform_constants: None,
        pair_factorization: None,
    }
}
fn chebyshev() -> EffectiveBkwRequest {
    request(&[&[1]], &[&[1], &[0, -2], &[1]])
}
fn tiling() -> EffectiveBkwRequest {
    request(
        &[&[1], &[0, -1]],
        &[&[1], &[-1, -1], &[0, -2, -1], &[0, 0, 1, 1]],
    )
}
fn path_cube() -> EffectiveBkwRequest {
    let mut r = request(
        &[&[1], &[0, -2], &[0, -2], &[0, 0, -1, 1]],
        &[
            &[1],
            &[-1, -3],
            &[0, -2, 2],
            &[0, -2, 1, 1],
            &[0, 0, -1, 2, -1],
        ],
    );
    r.uniform_constants = Some(EffectiveBkwUniformConstants {
        modulus_ratio: "4/5".into(),
        amplitude_lower_bound: "13/100".into(),
        numerator_weights: ["1", "8/5", "3/5", "1/2"].map(str::to_string).to_vec(),
        cofactor_weights: None,
    });
    r
}
#[test]
fn chebyshev_no_remainder_and_replay() {
    let r = certify_effective_bkw(&chebyshev()).unwrap();
    assert_eq!(r.status, "proved", "{}", r.reason);
    let c = r.certificate.as_ref().unwrap();
    assert_eq!(c.cutoff, 0);
    assert_eq!(c.domination_cutoff, Some(0));
    assert_eq!(c.root_region, EffectiveBkwRootRegion::WholeRealLine);
    assert_eq!(c.zero_count.offsets, [0]);
    r.verify().unwrap();
    let mut altered = r.clone();
    altered.certificate.as_mut().unwrap().domination_cutoff = Some(1);
    assert!(altered.verify().is_err());
    let mut altered = r;
    altered.certificate.as_mut().unwrap().signs[0].resultant[0] = "123".into();
    assert!(altered.verify().is_err());
}

#[test]
fn supplied_quartic_uniform_supplier_records_its_replay_mode() {
    let mut input = chebyshev();
    let h = B::new(vec![vec![q(1)], vec![-q(3) / q(10)], vec![q(1) / q(50)]]);
    input.generating_function.denominator = encode_b(
        &decode_b(&input.generating_function.denominator)
            .unwrap()
            .multiply(&h),
    );
    input.pair_factorization = Some(EffectiveBkwPairFactorization {
        real_part: vec!["0".into(), "1".into()],
        norm_squared: vec!["1".into()],
        cofactor: encode_b(&h),
    });
    input.uniform_constants = Some(EffectiveBkwUniformConstants {
        modulus_ratio: "1/2".into(),
        amplitude_lower_bound: "1/2".into(),
        numerator_weights: vec!["1".into()],
        cofactor_weights: Some(vec!["3/10".into(), "1/50".into()]),
    });
    let report = certify_effective_bkw(&input).unwrap();
    assert_eq!(report.status, "proved", "{}", report.reason);
    let c = report.certificate.as_ref().unwrap();
    assert!(c.domination.starts_with("uniform_compact_tail"));
    assert_eq!(c.domination_cutoff, Some(c.cutoff));
    assert!(c.branch_selection.starts_with("r=0"));
    report.verify().unwrap();
}

#[test]
fn order_two_rejects_uniform_constants_in_both_charts() {
    for mut input in [chebyshev(), request(&[&[1]], &[&[1], &[-1], &[0, -1]])] {
        input.uniform_constants = Some(EffectiveBkwUniformConstants {
            modulus_ratio: "1/2".into(),
            amplitude_lower_bound: "1".into(),
            numerator_weights: vec!["1".into()],
            cofactor_weights: None,
        });
        assert!(
            matches!(certify_effective_bkw(&input), Err(RecurrenceError::InvalidInput(message)) if message.contains("order-two"))
        );
    }
}

#[test]
fn phase_period_and_index_arithmetic_is_checked() {
    let base = certify_effective_bkw(&chebyshev())
        .unwrap()
        .certificate
        .unwrap();
    let mut out = base.clone();
    out.zero_count.offsets.clear();
    assert!(degree::close_phase(&mut out).is_err());
    let mut out = base.clone();
    out.zero_count.divisor = 0;
    assert!(degree::close_phase(&mut out).is_err());
    let mut out = base.clone();
    out.zero_count.divisor = usize::MAX;
    out.zero_count.slope = usize::MAX;
    assert!(degree::close_phase(&mut out)
        .unwrap_err()
        .to_string()
        .contains("overflow"));
    let mut out = base;
    out.cutoff = usize::MAX;
    assert!(degree::close_phase(&mut out)
        .unwrap_err()
        .to_string()
        .contains("overflow"));
    assert!(degree::value(
        &EffectiveBkwZeroCount {
            slope: usize::MAX,
            offsets: vec![0],
            divisor: 1,
        },
        usize::MAX
    )
    .is_err());
}

#[test]
fn denominator_norm_reuse_matches_direct_certificates() {
    let domain = Domain::select(&P::new(vec![q(-4), q(0), q(4)]), 2, &None).unwrap();
    let p = B::new(vec![vec![q(-2), q(0), q(-1)], vec![], vec![q(1)]]);
    let a = R::new(c(1), variable().add(&B::constant(q(3) / q(2)))).modulo(&p);
    let h = vec![
        R::of(c(1)),
        a.clone(),
        a.square(),
        a.square().square(),
        R::new(c(1), variable().add(&c(3))).modulo(&p),
    ];
    let direct: Vec<_> = std::iter::once(&a)
        .chain(&h)
        .enumerate()
        .filter(|(_, r)| r.d != B::one())
        .map(|(j, r)| {
            check_sign(
                &format!("pair/cofactor denominator {j} is nonzero"),
                &p,
                &r.d.power(2),
                &domain,
                Branch::LargestPositive,
                1,
            )
            .unwrap()
        })
        .collect();
    assert!(direct.iter().any(|c| !c.off_branch_boxes.is_empty()));
    assert_eq!(
        compact::check_pair_denominators(&p, &a, &h, &domain, Branch::LargestPositive).unwrap(),
        direct
    );
    let bad = R::new(c(1), p.clone());
    assert!(
        compact::check_pair_denominators(&p, &bad, &[], &domain, Branch::LargestPositive).is_err()
    );
}

#[test]
#[ignore = "exact irreducible quintic elimination; run explicitly with --release --ignored"]
fn automatic_irreducible_quintic_proved_and_replays() {
    let mut input = chebyshev();
    let mut h = B::one();
    for j in [1, 2, 3] {
        h = h.multiply(&B::new(vec![vec![q(1)], vec![-q(j) / q(20)]]));
    }
    let d = decode_b(&input.generating_function.denominator)
        .unwrap()
        .multiply(&h)
        .add(&B::monomial(5, 0, q(1) / q(1000000)));
    input.generating_function.denominator = encode_b(&d);
    let d: Vec<_> = d.rows().iter().map(|row| P::new(row.clone())).collect();
    let constant = P::new(d.iter().rev().map(|p| p.coeff(0)).collect());
    let coefficient = P::new(d.iter().rev().map(|p| p.coeff(1)).collect());
    // Linear in t and primitive in Q[lambda][t], hence irreducible.
    assert_eq!(constant.gcd(&coefficient).degree(), Some(0));
    assert!(higher::polynomial_pairs(&d).unwrap().is_empty());
    let report = certify_effective_bkw(&input).unwrap();
    assert_eq!(report.status, "proved", "{}", report.reason);
    assert!(report
        .certificate
        .as_ref()
        .unwrap()
        .branch_selection
        .contains("largest positive pair-product"));
    report.verify().unwrap();
}
#[test]
fn tiling_single_remainder_at_infinity() {
    let r = certify_effective_bkw(&tiling()).unwrap();
    assert_eq!(r.status, "proved", "{}", r.reason);
    let c = r.certificate.as_ref().unwrap();
    assert_eq!(c.zero_count.offsets, [0, -1]);
    assert_eq!(c.root_region, EffectiveBkwRootRegion::NegativeRay);
    assert_eq!(c.cutoff, 0);
    assert!(matches!(
        c.interval.lower,
        EffectiveBkwEndpoint::NegativeInfinity
    ));
    assert!(c.phase.iter().any(|p| p.contains("s=1/t")));
    assert!(c.signs.iter().all(|s| s.off_branch_boxes.is_empty()));
    r.verify().unwrap();
    let mut wrong_region = r.clone();
    wrong_region.certificate.as_mut().unwrap().root_region = EffectiveBkwRootRegion::PositiveRay;
    assert!(wrong_region.verify().is_err());
    let mut altered = r;
    altered.certificate.as_mut().unwrap().phase[0].push('!');
    assert!(altered.verify().is_err());
}

#[test]
fn positive_infinity_edge_uses_exact_reflection() {
    let input = reflect_request(&tiling()).unwrap();
    let report = certify_effective_bkw(&input).unwrap();
    assert_eq!(report.status, "proved", "{}", report.reason);
    assert_eq!(
        report.certificate.as_ref().unwrap().coordinate_transport,
        EffectiveBkwCoordinateTransport::ParameterReflectionWithAlternatingRows
    );
    assert_eq!(
        report.certificate.as_ref().unwrap().root_region,
        EffectiveBkwRootRegion::PositiveRay
    );
    report.verify().unwrap();
}
#[test]
#[cfg_attr(
    not(feature = "full-tests"),
    ignore = "exact quartic elimination; run with --release --features full-tests"
)]
fn path_cube_uniform_compact_certificate() {
    let r = certify_effective_bkw(&path_cube()).unwrap();
    assert_eq!(r.status, "proved", "{}", r.reason);
    let c = r.certificate.as_ref().unwrap();
    assert_eq!(c.cutoff, 39);
    assert_eq!(c.root_region, EffectiveBkwRootRegion::NegativeRay);
    assert_eq!(c.finite_rows.len(), 39);
    assert_eq!(c.finite_rows[1].coefficients, ["1", "1"]);
    assert_eq!(c.finite_rows[2].coefficients, ["1", "4", "1"]);
    assert_eq!(c.finite_rows[3].coefficients, ["1", "11", "11", "1"]);
    assert_eq!(
        c.tail_at_cutoff,
        "1060190164870166829055606784/9094947017729282379150390625"
    );
    assert!(c.signs.iter().any(|s| !s.off_branch_boxes.is_empty()));
    r.verify().unwrap();

    let mut automatic = path_cube();
    automatic.uniform_constants = None;
    let rule = crate::recurrence::Recurrence {
        terms: automatic
            .generating_function
            .denominator
            .iter()
            .enumerate()
            .skip(1)
            .map(|(offset, row)| {
                crate::recurrence::RecurrenceTerm::new(
                    offset,
                    0,
                    crate::recurrence::BivarPoly::from_coefficients_normalized(vec![
                        vec_poly::neg(decode(row).unwrap().coeffs()),
                    ]),
                )
            })
            .collect(),
        denominator: None,
        inhomogeneous: None,
    };
    let seeds: Vec<_> = c
        .finite_rows
        .iter()
        .take(4)
        .map(|row| decode(&row.coefficients).unwrap().coeffs().to_vec())
        .collect();
    let record =
        crate::recurrence::RecurrenceJson::from_recurrence_rational(&rule, 0, &seeds, None);
    let discovered = certify_effective_bkw_recurrence(&record).unwrap();
    assert_eq!(discovered.request, automatic);
    assert_eq!(discovered.status, "proved", "{}", discovered.reason);
    let certificate = discovered.certificate.as_ref().unwrap();
    eprintln!(
        "automatic compact cutoff {}: {}",
        certificate.cutoff, certificate.domination
    );
    assert_eq!(certificate.cutoff, 0);
    assert_eq!(certificate.domination_cutoff, Some(0));
    assert!(certificate
        .domination
        .starts_with("pointwise_two_real_remainders"));
    assert_eq!(certificate.finite_rows.len(), certificate.cutoff);
    assert_eq!(certificate.root_region, EffectiveBkwRootRegion::NegativeRay);
    discovered.verify().unwrap();
}
#[test]
fn rejects_nonreal_rows_and_incorrect_targets() {
    // A genuine non-real-rooted family: P_2=(1+t)^2+1.
    let bad = request(&[&[1]], &[&[1], &[-1, -1], &[-1]]);
    let r = certify_effective_bkw(&bad).unwrap();
    assert_eq!(r.status, "not_proved");
    assert!(!r.reason.is_empty());
    r.verify().unwrap();
    let mut wrong = chebyshev();
    wrong.target_zero_count = Some(EffectiveBkwZeroCount {
        slope: 1,
        offsets: vec![-1],
        divisor: 1,
    });
    assert!(certify_effective_bkw(&wrong)
        .unwrap()
        .reason
        .contains("degree budget"));
    // A compact pi-to-zero pair phase grows by one crossing per row;
    // this recurrence has the incompatible degree cap floor(n/5).
    let sublinear = request(&[&[1]], &[&[1], &[-1], &[], &[], &[], &[0, -1]]);
    let report = certify_effective_bkw(&sublinear).unwrap();
    assert_eq!(report.status, "not_proved");
    assert!(report.reason.contains("5 equal limiting modes"));
    assert!(report.reason.contains("no global domination cutoff"));
    report.verify().unwrap();
    let four = decode_b(
        &request(&[&[1]], &[&[1], &[-1], &[], &[], &[0, -1]])
            .generating_function
            .denominator,
    )
    .unwrap();
    assert_eq!(degree::equimodular_obstruction(&B::one(), &four), None);
}
#[test]
fn invalid_compact_constants_are_rejected() {
    let mut r = path_cube();
    r.uniform_constants.as_mut().unwrap().modulus_ratio = "1".into();
    assert!(certify_effective_bkw(&r).is_err());
}
#[test]
#[cfg_attr(
    not(feature = "full-tests"),
    ignore = "exact quartic negative controls; run with --release --features full-tests"
)]
fn path_cube_negative_controls() {
    for (q, m) in [("39/50", "13/100"), ("4/5", "7/50")] {
        let mut r = path_cube();
        let c = r.uniform_constants.as_mut().unwrap();
        c.modulus_ratio = q.into();
        c.amplitude_lower_bound = m.into();
        let report = certify_effective_bkw(&r).unwrap();
        assert_eq!(report.status, "not_proved", "{q} {m}");
        assert!(
            report.reason.contains("uniform modulus gap") || report.reason.contains("amplitude"),
            "{}",
            report.reason
        );
    }
}
#[test]
fn algebraic_refinement_keeps_the_requested_root() {
    let polynomial = P::new(vec![q(-2), q(0), q(1)]);
    // This interval contains +sqrt(2) only. A coarse enclosure of -sqrt(2)
    // nevertheless overlaps its left edge, so overlap cannot select the root.
    let mut endpoint = End {
        polynomial: polynomial.clone(),
        bounds: I::new(Q::new((-7071).into(), 5000.into()), q(2)),
    };
    let chain = crate::sturm::SturmChain::from_bigint_coeffs(&vec_poly::clear_denominators(
        polynomial.coeffs(),
    ));
    assert_eq!(
        chain.count_roots_in(&endpoint.bounds.lo, &endpoint.bounds.hi),
        1
    );
    endpoint.refine(&Q::new(1.into(), 256.into()));
    assert!(endpoint.bounds.lo > q(1));
    assert!(endpoint.bounds.hi < q(2));
    let saved = endpoint.bounds.wire();
    endpoint.refine(&Q::new(1.into(), 16.into()));
    assert_eq!(
        endpoint.bounds.wire(),
        saved,
        "coarser requests must not widen an enclosure"
    );
    endpoint.refine(&Q::new(1.into(), (1u64 << 32).into()));
    assert!(polynomial.evaluate(&endpoint.bounds.lo).is_negative());
    assert!(polynomial.evaluate(&endpoint.bounds.hi).is_positive());
}

#[test]
fn off_branch_resultant_zero_has_a_box() {
    let p = B::new(vec![vec![q(0), q(-1)], vec![q(1)], vec![], vec![q(1)]]); // r^3+r-t
    let f = variable().power(2).add(&c(1));
    let domain = Domain {
        wire: EffectiveBkwInterval {
            lower: EffectiveBkwEndpoint::NegativeInfinity,
            upper: EffectiveBkwEndpoint::PositiveInfinity,
        },
        lower: None,
        upper: None,
        sample: Q::zero(),
    };
    let cert = check_sign("r^2+1", &p, &f, &domain, Branch::UniqueReal, 1).unwrap();
    assert_eq!(cert.off_branch_boxes.len(), 1);
    let weak = check_sign(
        "r squared",
        &p,
        &variable().power(2),
        &domain,
        Branch::UniqueReal,
        0,
    )
    .unwrap();
    assert_eq!(weak.nonnegative_cell_samples.len(), 2);
    assert!(check_sign(
        "r squared strict",
        &p,
        &variable().power(2),
        &domain,
        Branch::UniqueReal,
        1
    )
    .is_err());
    assert!(check_sign(
        "r-1",
        &p,
        &variable().subtract(&c(1)),
        &domain,
        Branch::UniqueReal,
        -1
    )
    .is_err());
}
#[test]
fn higher_order_uniform_tail_is_exact() {
    let weights = [
        q(1),
        Q::new(8.into(), 5.into()),
        Q::new(3.into(), 5.into()),
        Q::new(1.into(), 2.into()),
    ];
    let (n, _) = uniform_tail_cutoff(
        4,
        &Q::new(4.into(), 5.into()),
        &Q::new(13.into(), 100.into()),
        &weights,
    )
    .unwrap();
    assert_eq!(n, 39);
    assert!(
        uniform_tail_cutoff(
            5,
            &Q::new(1.into(), 2.into()),
            &Q::new(1.into(), 2.into()),
            &[q(1)]
        )
        .unwrap()
        .0 > 1
    );
}

#[test]
fn supplied_higher_order_pair_factorization() {
    let mut input = chebyshev();
    let pair = decode_b(&input.generating_function.denominator).unwrap();
    let h = B::new(vec![
        vec![q(1)],
        vec![Q::new((-3).into(), 10.into())],
        vec![Q::new(3.into(), 100.into())],
        vec![Q::new((-1).into(), 1000.into())],
    ]);
    input.generating_function.denominator = encode_b(&pair.multiply(&h));
    input.pair_factorization = Some(EffectiveBkwPairFactorization {
        real_part: vec!["0".into(), "1".into()],
        norm_squared: vec!["1".into()],
        cofactor: encode_b(&h),
    });
    input.uniform_constants = Some(EffectiveBkwUniformConstants {
        modulus_ratio: "1/2".into(),
        amplitude_lower_bound: "1/2".into(),
        numerator_weights: vec!["1".into()],
        cofactor_weights: Some(["3/10", "3/100", "1/1000"].map(str::to_string).to_vec()),
    });
    let report = certify_effective_bkw(&input).unwrap();
    assert_eq!(report.status, "proved", "{}", report.reason);
    report.verify().unwrap();
    let mut automatic = input.clone();
    automatic.pair_factorization = None;
    automatic.uniform_constants = None;
    let report = certify_effective_bkw(&automatic).unwrap();
    assert_eq!(report.status, "proved", "{}", report.reason);
    report.verify().unwrap();
    input.pair_factorization.as_mut().unwrap().norm_squared[0] = "2".into();
    assert_eq!(certify_effective_bkw(&input).unwrap().status, "not_proved");

    // A canceled cofactor can vanish identically after substituting the
    // endpoint double-root formula. Reject the unsupported expansion without
    // dividing by that zero symbolic expression.
    let mut canceled = chebyshev();
    let h = decode_b(&[vec!["1".into()], vec!["0".into(), "-1".into()]]).unwrap();
    let pair = decode_b(&canceled.generating_function.denominator).unwrap();
    canceled.generating_function.numerator = encode_b(&h);
    canceled.generating_function.denominator = encode_b(&pair.multiply(&h));
    canceled.pair_factorization = Some(EffectiveBkwPairFactorization {
        real_part: vec!["0".into(), "1".into()],
        norm_squared: vec!["1".into()],
        cofactor: encode_b(&h),
    });
    let report = certify_effective_bkw(&canceled).unwrap();
    assert_eq!(report.status, "not_proved");
    assert!(report.reason.contains("endpoint pair cofactor"));
}

#[test]
fn newton_degrees_and_quadratic_ray_replay() {
    // Fibonacci polynomials: deg P_n=floor(n/2), including both residues.
    let input = request(&[&[1]], &[&[1], &[-1], &[0, -1]]);
    let report = certify_effective_bkw(&input).unwrap();
    assert_eq!(report.status, "proved", "{}", report.reason);
    let certificate = report.certificate.as_ref().unwrap();
    assert_eq!(certificate.zero_count.divisor, 2);
    assert_eq!(certificate.zero_count.offsets, [0, -1]);
    assert_eq!(certificate.phase_limits.left_eigenvalue_argument, "1/2");
    assert_eq!(certificate.root_region, EffectiveBkwRootRegion::NegativeRay);
    assert_eq!(
        certificate.degree_certificate.leading_numerator,
        ["1", "1", "-1"]
    );
    assert_eq!(
        certificate.degree_certificate.leading_denominator,
        ["1", "0", "-2", "0", "1"]
    );
    report.verify().unwrap();
    let mut tampered = report.clone();
    tampered.certificate.as_mut().unwrap().zero_count.divisor = 1;
    assert!(tampered.verify().is_err());
    let mut tampered = report;
    tampered
        .certificate
        .as_mut()
        .unwrap()
        .degree_certificate
        .newton
        .as_mut()
        .unwrap()
        .first_layers[1] += 1;
    assert!(tampered.verify().is_err());

    // A second Newton slope; degree evidence alone is not a root proof.
    let input = request(&[&[1]], &[&[1], &[-1], &[], &[0, 0, -1]]);
    let budget = degree_budget(
        &decode_b(&input.generating_function.numerator).unwrap(),
        &decode_b(&input.generating_function.denominator).unwrap(),
    )
    .unwrap();
    assert_eq!(budget.count.slope, 2);
    assert_eq!(budget.count.divisor, 3);
    for n in 0..30 {
        assert_eq!(degree::value(&budget.count, n).unwrap(), 2 * (n / 3));
    }
    let input = request(&[&[1]], &[&[1], &[-1], &[0, -1], &[0, 0, -1]]);
    let budget = degree_budget(
        &decode_b(&input.generating_function.numerator).unwrap(),
        &decode_b(&input.generating_function.denominator).unwrap(),
    )
    .unwrap();
    for n in 0..30 {
        assert_eq!(degree::value(&budget.count, n).unwrap(), 2 * n / 3);
    }
    let input = request(&[&[1], &[0, 1]], &[&[1], &[-1], &[0, -1]]);
    let report = certify_effective_bkw(&input).unwrap();
    assert_eq!(report.status, "proved", "{}", report.reason);
    report.verify().unwrap();
}

#[test]
fn general_leading_coefficients_and_quadratic_chart_selection() {
    for input in [
        request(&[&[2]], &[&[1], &[0, -2], &[1]]),
        request(&[&[10, 1]], &[&[1], &[0, -2], &[1]]),
        request(&[&[-4, 0, 1]], &[&[1], &[0, -2], &[1]]),
        request(&[&[1]], &[&[1], &[-4, -2], &[3, 3, 1]]),
    ] {
        let report = certify_effective_bkw(&input).unwrap();
        assert_eq!(report.status, "proved", "{}", report.reason);
        report.verify().unwrap();
    }
    // A fixed nonreal factor cannot be absorbed by exterior sign counts.
    let bad = request(&[&[1, 0, 1]], &[&[1], &[0, -2], &[1]]);
    assert_eq!(certify_effective_bkw(&bad).unwrap().status, "not_proved");
    // The diagonal top layer is the polynomial -1. After that finite
    // exception, the degree is n-1 and the leading coefficient is positive.
    let input = request(&[&[-1], &[3, 2]], &[&[1], &[0, -2], &[1]]);
    let report = certify_effective_bkw(&input).unwrap();
    assert_eq!(report.status, "proved", "{}", report.reason);
    let certificate = report.certificate.as_ref().unwrap();
    assert_eq!(certificate.zero_count.offsets, [-1]);
    assert_eq!(
        certificate.degree_certificate.leading_numerator,
        ["-1", "5"]
    );
    assert_eq!(
        certificate.degree_certificate.leading_denominator,
        ["1", "-2"]
    );
    report.verify().unwrap();
}

#[test]
fn recurrence_zero_prefix_is_shifted_and_bound_to_original() {
    let input: crate::recurrence::RecurrenceJson = serde_json::from_str(
        r#"{
        "schema":"polytool.recurrence.v1","first_index":7,
        "initial_polynomials":[["0"],["1"]],"recurrence":{"terms":[
        {"offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["1"]]}},
        {"offset":2,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["0","1"]]}}
    ]}}"#,
    )
    .unwrap();
    let report = certify_effective_bkw_recurrence(&input).unwrap();
    assert_eq!(report.status, "proved", "{}", report.reason);
    assert_eq!(report.request.generating_function.first_index, 8);
    assert_eq!(report.skipped_prefix.len(), 1);
    assert_eq!(report.skipped_prefix[0].degree, None);
    report.verify().unwrap();
    let mut tampered = report.clone();
    tampered
        .source_recurrence
        .as_mut()
        .unwrap()
        .initial_polynomials[0] = vec!["1".into()];
    assert!(tampered.verify().is_err());
    let mut tampered = report;
    tampered.skipped_prefix[0].index += 1;
    assert!(tampered.verify().is_err());

    // The order-two rule starts only after all three prescribed seeds.
    // Keeping that exception gives P_2=1 and the Chebyshev tail from there.
    let extra: crate::recurrence::RecurrenceJson = serde_json::from_str(
        r#"{
        "schema":"polytool.recurrence.v1","first_index":0,
        "initial_polynomials":[["2"],["0"],["1"]],"recurrence":{"terms":[
        {"offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["0","2"]]}},
        {"offset":2,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["-1"]]}}
    ]}}"#,
    )
    .unwrap();
    let report = certify_effective_bkw_recurrence(&extra).unwrap();
    assert_eq!(report.status, "proved", "{}", report.reason);
    assert_eq!(report.request.generating_function.first_index, 2);
    assert_eq!(report.skipped_prefix.len(), 2);
    report.verify().unwrap();
}

#[test]
fn newton_cubic_ray_with_a_small_real_remainder() {
    let mut input = request(&[&[1]], &[&[1], &[-1], &[0, -1]]);
    let pair = decode_b(&input.generating_function.denominator).unwrap();
    let h = B::new(vec![vec![q(1)], vec![q(-1) / q(10)]]);
    input.generating_function.denominator = encode_b(&pair.multiply(&h));
    let report = certify_effective_bkw(&input).unwrap();
    assert_eq!(report.status, "proved", "{}", report.reason);
    assert_eq!(report.certificate.as_ref().unwrap().zero_count.divisor, 2);
    report.verify().unwrap();
}

#[test]
#[cfg_attr(
    not(feature = "full-tests"),
    ignore = "automatic factorizations through order eight"
)]
fn automatic_polynomial_pairs_through_order_eight() {
    for order in 5..=8 {
        let mut input = chebyshev();
        let mut h = B::one();
        for j in 1..=order - 2 {
            h = h.multiply(&B::new(vec![
                vec![q(1)],
                vec![-Q::new(1.into(), (100 * j).into())],
            ]));
        }
        input.generating_function.denominator = encode_b(
            &decode_b(&input.generating_function.denominator)
                .unwrap()
                .multiply(&h),
        );
        let report = certify_effective_bkw(&input).unwrap();
        assert_eq!(report.status, "proved", "order {order}: {}", report.reason);
        report.verify().unwrap();
    }
    // Two nonreal pairs: the constant quadratic cofactor is not dominant.
    // Automatic search must continue to the other quadratic factor.
    let mut input = chebyshev();
    let mut h = B::new(vec![vec![q(1)], vec![], vec![q(1) / q(10000)]]);
    for j in [1, 2] {
        h = h.multiply(&B::new(vec![vec![q(1)], vec![-q(j) / q(100)]]));
    }
    input.generating_function.denominator = encode_b(
        &decode_b(&input.generating_function.denominator)
            .unwrap()
            .multiply(&h),
    );
    let report = certify_effective_bkw(&input).unwrap();
    assert_eq!(report.status, "proved", "{}", report.reason);
    report.verify().unwrap();
}

#[test]
#[ignore = "bounded diagnostic for the higher-order algebraic subresultant"]
fn algebraic_pair_subresultant_factorization_identity() {
    let mut input = chebyshev();
    let mut h = B::one();
    for j in [1, 3, 7] {
        h = h.multiply(&B::new(vec![
            vec![q(1)],
            vec![-Q::new(j.into(), 100.into())],
        ]));
    }
    let mut d = decode_b(&input.generating_function.denominator)
        .unwrap()
        .multiply(&h);
    d = d.add(&B::monomial(5, 0, q(1) / q(100000000)));
    input.generating_function.denominator = encode_b(&d);
    let d: Vec<_> = d.rows().iter().map(|row| P::new(row.clone())).collect();
    // A primitive polynomial of degree one in t is irreducible in Q[t,lambda].
    let constant = P::new(d.iter().rev().map(|p| p.coeff(0)).collect());
    let coefficient = P::new(d.iter().rev().map(|p| p.coeff(1)).collect());
    assert!(d.iter().all(|p| p.degree().unwrap_or(0) <= 1));
    assert_eq!(constant.gcd(&coefficient).degree(), Some(0));
    let (p, a, b, h, _) = higher::algebraic_pair(&d).unwrap();
    let pair = [R::of(c(1)), a.scale(&q(-2)), b];
    let mut product = vec![R::of(c(0)); 6];
    for (i, u) in pair.iter().enumerate() {
        for (j, v) in h.iter().enumerate() {
            product[i + j] = product[i + j].add(&u.mul(v));
        }
    }
    for (i, v) in product.iter().enumerate() {
        assert!(
            rem(&v.sub(&R::of(t(&d[i]))).n, &p).is_zero(),
            "coefficient {i}"
        );
    }
}

#[test]
fn a_cofactor_collision_cannot_supply_a_pair_phase_endpoint() {
    let interval = EffectiveBkwInterval {
        lower: EffectiveBkwEndpoint::Rational { value: "-2".into() },
        upper: EffectiveBkwEndpoint::Rational {
            value: "1/2".into(),
        },
    };
    let domain = Domain::select(&P::new(vec![q(-1)]), 5, &Some(interval)).unwrap();
    // The dominant pair has |lambda|=2, Re(lambda)=t. Two small positive
    // cofactor roots collide at t=1/2, x=2/5; that collision is not the pair.
    let pair = B::new(vec![vec![q(1)], vec![q(0), q(-2)], vec![q(4)]]);
    let real = B::new(vec![vec![q(1)], vec![-q(1) / q(100)]]);
    let small = B::new(vec![
        vec![q(1)],
        vec![-q(4) / q(5)],
        vec![q(31) / q(200), q(1) / q(100)],
    ]);
    let actual = real.multiply(&small);
    let wrong = real.multiply(&pair);
    let h = |p: &B| {
        p.rows()
            .iter()
            .map(|row| R::of(t(&P::new(row.clone()))))
            .collect::<Vec<_>>()
    };
    let gamma = [q(1) / q(2), q(1) / q(16), q(1) / q(4096)];
    let ratio = q(3) / q(4);
    assert!(
        gamma
            .iter()
            .enumerate()
            .map(|(j, g)| g / ratio.pow((j + 1) as i32))
            .sum::<Q>()
            < q(1)
    );
    let actual_h = h(&actual);
    for (j, (coefficient, weight)) in actual_h.iter().skip(1).zip(&gamma).enumerate() {
        let bound = R::of(c(4))
            .pow(j + 1)
            .scale(&(weight * weight))
            .sub(&coefficient.square());
        check_sign(
            "actual cofactor Cauchy bound",
            &variable(),
            &bound.signed_polynomial(&variable()),
            &domain,
            Branch::Zero,
            0,
        )
        .unwrap();
    }
    compact::check_collision_cofactor(&domain, &R::of(c(-2)), &actual_h, &gamma, true).unwrap();
    let parameter = q(1) / q(2);
    let d = P::new(pair.multiply(&actual).evaluate_variable(&parameter));
    let collision = P::new(vec![q(1), -q(4) / q(5), q(4) / q(25)]);
    assert_eq!(d, collision * P::new(wrong.evaluate_variable(&parameter)));
    let rejected = compact::check_collision_cofactor(
        &domain,
        &R::of(c(2)).scale(&(q(1) / q(5))),
        &h(&wrong),
        &gamma,
        false,
    );
    assert!(
        rejected.is_err(),
        "a cofactor double root must not fix the pair phase limit"
    );
    // Equality in an endpoint weight is allowed; the ratio budget is strict.
    compact::check_collision_cofactor(
        &domain,
        &R::of(c(-2)),
        &[R::of(c(1)), R::of(c(1))],
        &[q(1) / q(2)],
        true,
    )
    .unwrap();
}

#[test]
fn rational_pi_angles_use_complete_cyclotomic_tests() {
    for (square, angle) in [
        (q(0), q(1) / q(2)),
        (q(1) / q(4), q(1) / q(3)),
        (q(1) / q(2), q(1) / q(4)),
        (q(3) / q(4), q(1) / q(6)),
        (q(1), q(0)),
    ] {
        assert_eq!(
            quadratic::angle_from_cosine_squared(&square, 1).unwrap(),
            angle
        );
        assert_eq!(
            quadratic::angle_from_cosine_squared(&square, -1).unwrap(),
            q(1) - angle
        );
    }
    for (square, angle) in [
        (q(0), q(0)),
        (q(1) / q(3), q(1) / q(6)),
        (q(1), q(1) / q(4)),
        (q(3), q(1) / q(3)),
    ] {
        assert_eq!(
            quadratic::angle_from_tangent_squared(&square, -1).unwrap(),
            -angle
        );
    }
    for denominator in 2..=24 {
        for numerator in 0..=denominator {
            let value = q(numerator) / q(denominator);
            let twice_cosine = q(4) * &value - q(2);
            // Roots of unity are algebraic integers. A rational trace is
            // integral; conversely this quadratic with integral trace and
            // all roots on the unit circle is cyclotomic here.
            let expected = twice_cosine.denom() == &num_bigint::BigInt::one();
            assert_eq!(
                quadratic::angle_from_cosine_squared(&value, 1).is_ok(),
                expected
            );
        }
    }
    assert!(quadratic::angle_from_cosine_squared(&q(-1), 1).is_err());
    assert!(quadratic::angle_from_cosine_squared(&q(2), 1).is_err());
    assert!(quadratic::angle_from_tangent_squared(&q(2), 1).is_err());
    assert!(quadratic::angle_from_tangent_squared(&q(-1), 1).is_err());
}

#[test]
fn binomial_edge_dilation_uses_general_bezout_exponents() {
    for (degree, k) in [(3, 1), (3, 2), (3, 4), (3, 5), (4, 3), (5, 2), (7, 3)] {
        let p = B::monomial(degree, 0, q(1)).add(&B::monomial(
            0,
            k,
            q(if k % 2 == 0 { -2 } else { 2 }),
        ));
        let (eigen, parameter) = cubic::binomial_newton_dilation(&p, k, degree);
        assert_eq!(eigen.pow(degree as i32), q(2) * parameter.pow(k as i32));
        assert_ne!((eigen, parameter), (q(1), q(1)));
    }
    // For the cubic integer-slope edge, rational dilation exists only when
    // c is already a cube; the original exact linear-factor search handles it.
    let p = B::monomial(3, 0, q(1)).add(&B::monomial(0, 3, q(2)));
    assert_eq!(cubic::binomial_newton_dilation(&p, 1, 1), (q(1), q(1)));
}

#[test]
fn newton_budget_rejects_a_vanishing_residue_class_quickly() {
    // P_0 = 1 + 2t, P_1 = 0, P_n = (2t^2 - t - 1) P_(n-2): every odd row is zero.
    // The layer search would never find an odd positive layer; reject at once.
    let numerator = B::new(vec![vec![q(1), q(2)]]);
    let denominator = B::new(vec![vec![q(1)], vec![], vec![q(1), q(1), q(-2)]]);
    let started = std::time::Instant::now();
    let error = degree::newton_budget(&numerator, &denominator)
        .err()
        .unwrap()
        .to_string();
    assert!(error.contains("vanish"), "{error}");
    assert!(started.elapsed() < std::time::Duration::from_secs(5));
}

#[test]
fn local_newton_unit_check_and_branch_denominator_guards() {
    let numerator = B::one();
    let denominator = B::new(vec![vec![q(1), q(1)], vec![q(-1)], vec![q(0), q(-1)]]);
    assert!(degree::newton_budget(&numerator, &denominator)
        .err()
        .unwrap()
        .to_string()
        .contains("D(t,0)=1"));
    let interval = EffectiveBkwInterval {
        lower: EffectiveBkwEndpoint::Rational { value: "-1".into() },
        upper: EffectiveBkwEndpoint::Rational { value: "1".into() },
    };
    let domain = Domain::select(&P::new(vec![q(-1)]), 2, &Some(interval)).unwrap();
    // This pole is away from the sample. Later zero-numerator arithmetic
    // could erase it, so the fractions are checked before that arithmetic.
    let mut fraction = R::of(c(0));
    fraction.d = t(&P::new(vec![-q(1) / q(2), q(1)]));
    assert!(
        compact::check_pair_denominators(&variable(), &fraction, &[], &domain, Branch::Zero)
            .is_err()
    );
    assert!(compact::check_pair_denominators(
        &variable(),
        &R::of(c(0)),
        &[fraction],
        &domain,
        Branch::Zero
    )
    .is_err());
    let mut regular = R::of(c(0));
    regular.d = t(&P::new(vec![q(1), q(0), q(1)]));
    assert_eq!(
        compact::check_pair_denominators(&variable(), &regular, &[], &domain, Branch::Zero)
            .unwrap()
            .len(),
        1
    );
}

#[test]
fn automatic_cauchy_weights_allow_identically_zero_coefficients() {
    let mut input = chebyshev();
    let small_pair = B::new(vec![vec![q(1)], vec![], vec![q(1) / q(10000)]]);
    let h = small_pair.power(2);
    // N=H makes F=N/H identically one, isolating the Cauchy supplier
    // from the separate restricted amplitude-winding test.
    input.generating_function.numerator = encode_b(&h);
    input.pair_factorization = Some(EffectiveBkwPairFactorization {
        real_part: vec!["0".into(), "1".into()],
        norm_squared: vec!["1".into()],
        cofactor: encode_b(&h),
    });
    input.generating_function.denominator = encode_b(
        &decode_b(&input.generating_function.denominator)
            .unwrap()
            .multiply(&h),
    );
    let report = certify_effective_bkw(&input).unwrap();
    assert_eq!(report.status, "proved", "{}", report.reason);
    let constants = report
        .certificate
        .as_ref()
        .unwrap()
        .automatic_uniform_constants
        .as_ref()
        .unwrap();
    let gamma = constants.cofactor_weights.as_ref().unwrap();
    assert_eq!(gamma[0], "0");
    assert_eq!(gamma[2], "0");
    report.verify().unwrap();
}

#[test]
fn generic_triple_edge_and_equidominant_newton_face() {
    // [z^n] 1/((1-z)^3-t*z^3) = sum_k binom(n+2,3k+2)*t^k.
    let request = request(&[&[1]], &[&[1], &[-3], &[3], &[-1, -1]]);
    let report = certify_effective_bkw(&request).unwrap();
    assert_eq!(report.status, "proved", "{}", report.reason);
    let c = report.certificate.as_ref().unwrap();
    assert_eq!(c.domination_cutoff, Some(0));
    assert_eq!(c.phase_limits.left_eigenvalue_argument, "1/3");
    assert_eq!(c.phase_limits.right_amplitude_argument, "-2/3");
    assert_eq!(c.zero_count.divisor, 3);
    report.verify().unwrap();
    let mut altered = report.clone();
    altered
        .certificate
        .as_mut()
        .unwrap()
        .phase_limits
        .right_amplitude_argument = "-1/2".into();
    assert!(altered.verify().is_err());
    let mut altered = report;
    altered.certificate.as_mut().unwrap().domination_cutoff = Some(1);
    assert!(altered.verify().is_err());
}

#[test]
#[cfg_attr(
    not(feature = "full-tests"),
    ignore = "exact quartic Newton elimination"
)]
fn quartic_fractional_ray_with_two_fixed_real_modes() {
    // D=((1-2z)^2-t*z^2)*(1-z^2). The pair tends to angle pi/2
    // at infinity, with two fixed real cofactor modes +1 and -1.
    let r = request(&[&[1]], &[&[1], &[-4], &[3, -1], &[4], &[-4, 1]]);
    let report = certify_effective_bkw(&r).unwrap();
    assert_eq!(report.status, "proved", "{}", report.reason);
    let c = report.certificate.as_ref().unwrap();
    assert_eq!(c.phase_limits.left_eigenvalue_argument, "1/2");
    assert_eq!(c.domination_cutoff, Some(0));
    assert_eq!(c.zero_count.divisor, 2);
    report.verify().unwrap();
}

#[test]
fn legacy_uniform_and_single_mode_reports_require_regeneration() {
    for input in [
        chebyshev(),
        tiling(),
        request(&[&[1]], &[&[1], &[-1], &[0, -1], &[0, -1]]),
    ] {
        let report = certify_with_mode(&input, true).unwrap();
        assert_eq!(report.status, "proved", "{}", report.reason);
        assert!(report
            .certificate
            .as_ref()
            .unwrap()
            .domination_cutoff
            .is_none());
        assert!(report
            .verify()
            .unwrap_err()
            .to_string()
            .contains("legacy certificate"));
        certify_effective_bkw(&input).unwrap().verify().unwrap();
        let mut bad = report;
        bad.certificate.as_mut().unwrap().cutoff += 1;
        assert!(bad.verify().is_err());
    }
}

#[test]
fn pointwise_conjugate_cofactor_and_real_amplitude_lift() {
    let mut r = chebyshev();
    let h = B::new(vec![vec![q(1)], vec![], vec![q(1) / q(100)]]);
    r.generating_function.denominator = encode_b(
        &decode_b(&r.generating_function.denominator)
            .unwrap()
            .multiply(&h),
    );
    r.pair_factorization = Some(EffectiveBkwPairFactorization {
        real_part: vec!["0".into(), "1".into()],
        norm_squared: vec!["1".into()],
        cofactor: encode_b(&h),
    });
    let report = certify_effective_bkw(&r).unwrap();
    assert_eq!(report.status, "proved", "{}", report.reason);
    let c = report.certificate.as_ref().unwrap();
    assert_eq!(c.domination_cutoff, Some(0));
    assert!(c.domination.starts_with("pointwise_conjugate_remainder"));
    assert!(c
        .signs
        .iter()
        .any(|s| s.name == "positive real part of pair amplitude"));
    report.verify().unwrap();
    let mut bad = report;
    bad.certificate.as_mut().unwrap().expansion_data.cofactor[2].numerator[0][0] = "1/101".into();
    assert!(bad.verify().is_err());
}

#[test]
fn collision_factor_preconditioning_preserves_all_resultant_roots() {
    let left = P::new(vec![q(2), q(1)]);
    let right = P::new(vec![q(-3), q(1)]);
    let interior = P::new(vec![q(0), q(1)]);
    let outside = P::new(vec![q(-5), q(1)]);
    let nonreal = P::new(vec![q(1), q(0), q(1)]);
    let domain = Domain::select(
        &(left.clone() * right.clone()),
        2,
        &Some(EffectiveBkwInterval {
            lower: EffectiveBkwEndpoint::Rational { value: "-2".into() },
            upper: EffectiveBkwEndpoint::Rational { value: "3".into() },
        }),
    )
    .unwrap();
    let resultant = left.pow(12) * right.pow(7) * interior.pow(4) * outside.pow(2) * nonreal.pow(3);
    let reduced = signs::critical_polynomial(&resultant, &domain);
    assert_eq!(reduced.degree(), Some(6));
    assert!(
        reduced
            .gcd(&(left * right * interior * outside * nonreal))
            .degree()
            == Some(6)
    );
    assert_eq!(domain.critical_roots(&reduced).unwrap().len(), 1);
}

#[test]
fn degree_zero_exceptions_do_not_bypass_the_phase_slope() {
    let report = certify_effective_bkw(&chebyshev()).unwrap();
    let mut c = report.certificate.unwrap();
    c.zero_count = EffectiveBkwZeroCount {
        slope: 0,
        offsets: vec![0],
        divisor: 1,
    };
    c.phase_limits.left_eigenvalue_argument = "1/2".into();
    c.phase_limits.left_amplitude_argument = "-10".into();
    c.phase_limits.right_amplitude_argument = "0".into();
    assert!(degree::close_phase(&mut c).is_err());
}

#[test]
fn parity_upper_bounds_do_not_supply_generic_exterior_signs() {
    let report = certify_effective_bkw(&tiling()).unwrap();
    let mut c = report.certificate.as_ref().unwrap().clone();
    assert!(c.degree_certificate.newton.is_none());
    assert_eq!(c.zero_count.offsets, [0, -1]);
    assert!(degree::close_phase(&mut c).is_err());
    // Its original parity phase argument is independent of this closure.
    report.verify().unwrap();
}

#[test]
fn five_mode_zero_edge_obstructs_global_pair_domination() {
    let input = request(&[&[1]], &[&[1], &[], &[0, -2], &[], &[], &[0, -7]]);
    let report = certify_effective_bkw(&input).unwrap();
    assert_eq!(report.status, "not_proved");
    assert!(report.reason.contains("zero-end Newton face has 5"));
    assert!(report
        .reason
        .contains("excluding global pair domination near t=0"));
    report.verify().unwrap();
}

#[test]
fn direct_witness_checker_consumes_samples_and_factorization() {
    let report = certify_effective_bkw(&chebyshev()).unwrap();
    replay::witnesses(&report).unwrap();
    let mut alternative = report.clone();
    // A different exact sample of the same strict discriminant sign is valid
    // evidence even though deterministic producer comparison would reject it.
    alternative.certificate.as_mut().unwrap().signs[0].sample_parameter = "1/2".into();
    replay::witnesses(&alternative).unwrap();
    alternative.certificate.as_mut().unwrap().signs[0].sample_parameter = "2".into();
    assert!(replay::witnesses(&alternative).is_err());
    let mut bad = report.clone();
    bad.certificate
        .as_mut()
        .unwrap()
        .expansion_data
        .pair_norm_squared
        .numerator = vec![vec!["2".into()]];
    assert!(replay::witnesses(&bad).is_err());
    let mut bad = report;
    bad.certificate.as_mut().unwrap().domination_cutoff = None;
    assert!(bad
        .verify()
        .unwrap_err()
        .to_string()
        .contains("legacy certificate"));
}
