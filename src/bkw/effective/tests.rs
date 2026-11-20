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
    assert_eq!(c.cutoff, 1);
    assert_eq!(c.root_region, EffectiveBkwRootRegion::WholeRealLine);
    assert_eq!(c.zero_count.offsets, [0]);
    r.verify().unwrap();
    let mut altered = r.clone();
    altered.certificate.as_mut().unwrap().finite_rows[0].coefficients[0] = "2".into();
    assert!(altered.verify().is_err());
    let mut altered = r;
    altered.certificate.as_mut().unwrap().signs[0].resultant[0] = "123".into();
    assert!(altered.verify().is_err());
}
#[test]
fn tiling_single_remainder_at_infinity() {
    let r = certify_effective_bkw(&tiling()).unwrap();
    assert_eq!(r.status, "proved", "{}", r.reason);
    let c = r.certificate.as_ref().unwrap();
    assert_eq!(c.zero_count.offsets, [0, -1]);
    assert_eq!(c.root_region, EffectiveBkwRootRegion::NegativeRay);
    assert_eq!(c.cutoff, 1);
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
    assert!(certificate.cutoff <= 59);
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
    });
    assert!(certify_effective_bkw(&wrong)
        .unwrap()
        .reason
        .contains("degree budget"));
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
