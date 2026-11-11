//! All-index rational PF Riordan certificates and rejection/regression cases.

use num_traits::Zero;
use polytool::recurrence::*;
use polytool::vec_poly;

fn q(n: i64) -> BigRational {
    BigRational::from_integer(n.into())
}
fn p(v: &[i64]) -> Vec<BigRational> {
    v.iter().map(|&n| q(n)).collect()
}

// Canonical rational Riordan generator from d=dn/dd, K=kn/kd. The recurrence
// is the denominator identity, with all numerator corrections included as seeds.
fn family(dn: &[i64], dd: &[i64], kn: &[i64], kd: &[i64]) -> (Recurrence, Vec<Vec<BigRational>>) {
    let n = vec_poly::mul(&p(dn), &p(kd));
    let d0 = vec_poly::mul(&p(dd), &p(kd));
    let l = vec_poly::shift_by(&vec_poly::mul(&p(dn), &p(kn)), 1);
    let order = d0.len().max(l.len()) - 1;
    let mut terms = Vec::new();
    for lag in 1..=order {
        let coeff = vec![
            -d0.get(lag).cloned().unwrap_or_else(BigRational::zero),
            l.get(lag).cloned().unwrap_or_else(BigRational::zero),
        ];
        if !vec_poly::is_zero(&coeff) {
            terms.push(RecurrenceTerm::new(
                lag,
                0,
                BivarPoly::from_coefficients_normalized(vec![coeff]),
            ));
        }
    }
    let rule = Recurrence {
        terms,
        denominator: None,
        inhomogeneous: None,
    };
    let mut seeds: Vec<Vec<BigRational>> = Vec::new();
    for i in 0..order.max(n.len()) {
        let mut row = vec![n.get(i).cloned().unwrap_or_else(BigRational::zero)];
        for lag in 1..=i.min(order) {
            let c = vec![
                -d0.get(lag).cloned().unwrap_or_else(BigRational::zero),
                l.get(lag).cloned().unwrap_or_else(BigRational::zero),
            ];
            row = vec_poly::add(&row, &vec_poly::mul(&c, &seeds[i - lag]));
        }
        seeds.push(vec_poly::trim(&row));
    }
    (rule, seeds)
}

fn check(rule: &Recurrence, seeds: &[Vec<BigRational>]) -> RealRootednessTheoremReport {
    check_rational_pf_riordan_theorem(rule, 0, seeds, &Default::default()).unwrap()
}

#[test]
fn exact_audit_pairs_and_algebraic_factors_have_replayable_proofs() {
    for (dn, dd, kn, kd) in [
        (&[1][..], &[1, -1][..], &[1, 2, 1][..], &[1, -1][..]), // C_n
        (&[1][..], &[1, -1][..], &[0, 0, 1][..], &[1, -1][..]), // R_m
        (&[1][..], &[1, -2, 1][..], &[0, 0, 1][..], &[1][..]),  // H_m
        (&[1, 1][..], &[1, -2, 1][..], &[1, 1][..], &[1][..]),  // A216182
        (&[1, 3, 1][..], &[1, -3, 1][..], &[0, 0, 2][..], &[1][..]), // irrational factors
    ] {
        let (rule, seeds) = family(dn, dd, kn, kd);
        let report = check(&rule, &seeds);
        assert!(report.is_proved(), "{report:?}");
        assert!(!report.conclusion.as_ref().unwrap().simple_roots);
        assert!(report.conclusion.as_ref().unwrap().interlacing.is_none());
        report.verify().unwrap();
        report
            .rational_pf_riordan_certificate()
            .unwrap()
            .verify_identity(&rule, 0, &seeds)
            .unwrap();
        let rows = rule.generate_rows_rational(&seeds, 0, 30).unwrap();
        for row in rows {
            assert!(polytool::is_real_rooted_bigint_coeffs(
                &row.iter().map(|c| c.to_integer()).collect::<Vec<_>>()
            ));
        }
    }
}

#[test]
fn rational_pf_decision_checks_signs_cancellation_and_multiplicities() {
    let series = RationalGeneratingSeries::new(&p(&[0, 0, 1, 2, 1]), &p(&[1, -2, 1])).unwrap();
    let evidence = check_rational_polya_frequency(&series).unwrap();
    assert!(evidence.is_polya_frequency);
    assert_eq!(evidence.valuation, Some(2));
    assert_eq!(evidence.numerator_layers.len(), 2);
    assert_eq!(evidence.denominator_layers.len(), 2);
    for (n, d) in [
        (&[-1][..], &[1][..]),
        (&[1, -1][..], &[1][..]),
        (&[1, 0, 1][..], &[1][..]),
        (&[1][..], &[1, 1][..]),
        (&[0][..], &[1][..]),
    ] {
        assert!(
            !check_rational_polya_frequency(&RationalGeneratingSeries::new(&p(n), &p(d)).unwrap())
                .unwrap()
                .is_polya_frequency
        );
    }
    let cancelled = RationalGeneratingSeries::new(&p(&[1, -1]), &p(&[1, -2, 1])).unwrap();
    assert!(
        check_rational_polya_frequency(&cancelled)
            .unwrap()
            .is_polya_frequency
    );
}

#[test]
fn noninvertible_array_reports_periodic_and_finite_zero_rows_exactly() {
    // F=1/(1-xz^3), all but multiples of three vanish.
    let (rule, seeds) = family(&[1], &[1], &[0, 0, 1], &[1]);
    let report = check(&rule, &seeds);
    assert!(report.is_proved());
    let conclusion = report.conclusion.as_ref().unwrap();
    assert_eq!(conclusion.first_index, 0);
    assert_eq!(
        conclusion.zero_residue_classes,
        Some(ZeroResidueClasses {
            period: 3,
            residues: vec![1, 2]
        })
    );
    report.verify().unwrap();
    // d=1, h=z^3+z^4: supports [3k,4k], gaps 1,2,5 only.
    let (rule, seeds) = family(&[1], &[1], &[0, 0, 1, 1], &[1]);
    let report = check(&rule, &seeds);
    assert!(report.is_proved());
    assert_eq!(report.conclusion.as_ref().unwrap().first_index, 6);
    assert_eq!(
        report
            .rational_pf_riordan_certificate()
            .unwrap()
            .finite_zero_indices,
        vec![1, 2, 5]
    );
    report.verify().unwrap();
    let record = RecurrenceJson::from_recurrence_rational(&rule, 0, &seeds, None);
    let discovery = discover_recurrence_proofs(&record, &Default::default()).unwrap();
    discovery.verify().unwrap();
}

#[test]
fn geometric_row_counterexample_nonreal_seeds_and_bad_pf_catalog_entries_fail() {
    // d=1/(1-z), h=z, hence K=1-z, geometric-sum rows.
    let (rule, seeds) = family(&[1], &[1, -1], &[1, -1], &[1]);
    assert!(!check(&rule, &seeds).is_proved());
    let rows = rule.generate_rows_rational(&seeds, 0, 4).unwrap();
    assert_eq!(rows[2], p(&[1, 1, 1]));
    assert!(!polytool::is_real_rooted_bigint_coeffs(&[
        1.into(),
        1.into(),
        1.into()
    ]));
    // In contrast 1/(1-z-xz) has K=1 and binomial rows: it must pass.
    let (rule, seeds) = family(&[1], &[1, -1], &[1], &[1]);
    assert!(check(&rule, &seeds).is_proved());
    let (rule, mut seeds) = family(&[1], &[1, -1], &[0, 0, 1], &[1]);
    seeds[1] = p(&[1, 4]);
    seeds[2] = p(&[1, 2, 1]);
    assert!(!check(&rule, &seeds).is_proved());
    for id in [
        "A109954", "A156308", "A158454", "A158909", "A181327", "A210039", "A275442", "A277513",
        "A334009", "A046741", "A082649", "A092879", "A245963", "A177970",
    ] {
        eprintln!("checking {id}");
        let entry = polytool::oeis::by_id(id).unwrap();
        let (r, seeds) = entry.recurrence_parts().unwrap();
        let report = check_rational_pf_riordan_theorem(
            &r,
            entry.recurrence_first_index,
            &seeds,
            &Default::default(),
        )
        .unwrap();
        assert!(!report.is_proved(), "{id}");
    }
}

#[test]
fn tampered_identity_factorization_and_zero_claims_are_rejected() {
    let (rule, seeds) = family(&[1, 1], &[1, -2, 1], &[1, 1], &[1]);
    let report = check(&rule, &seeds);
    let original = serde_json::to_value(&report).unwrap();
    for pointer in [
        "/certificate/rational_pf_riordan/generating_function/numerator/0/0",
        "/certificate/rational_pf_riordan/decomposition/k/numerator/0",
        "/certificate/rational_pf_riordan/d_evidence/numerator_layers/0/polynomial/0",
    ] {
        let mut value = original.clone();
        *value.pointer_mut(pointer).unwrap() = serde_json::json!("99");
        let forged: RealRootednessTheoremReport = serde_json::from_value(value).unwrap();
        assert!(forged.verify().is_err());
        assert!(forged
            .rational_pf_riordan_certificate()
            .unwrap()
            .verify_identity(&rule, 0, &seeds)
            .is_err());
    }
    let mut forged = report.clone();
    if let Some(RouteCertificate::RationalPfRiordan(c)) = &mut forged.certificate {
        c.finite_zero_indices.push(7);
    }
    assert!(forged.verify().is_err());
    let mut forged = report;
    forged.certificate = None;
    assert!(forged.verify().is_err());
}

#[test]
fn target_routes_reversals_sections_and_marker_shift_compose() {
    for id in [
        "A098172", "A098173", "A114123", "A118424", "A118884", "A184883", "A190088", "A190152",
        "A216182", "A277627",
    ] {
        eprintln!("checking {id}");
        let entry = polytool::oeis::by_id(id).unwrap();
        let (rule, seeds) = entry.recurrence_parts().unwrap();
        let record = RecurrenceJson::from_recurrence_rational(
            &rule,
            entry.recurrence_first_index,
            &seeds,
            None,
        );
        let report = discover_recurrence_proofs(&record, &Default::default()).unwrap();
        assert_eq!(
            report.status,
            "proved",
            "{id}: {}",
            serde_json::to_string(&report).unwrap()
        );
        report.verify().unwrap();
    }
}

#[test]
fn requested_open_ray_does_not_admit_zero_roots() {
    let (rule, seeds) = family(&[1], &[1], &[1], &[1]);
    let report = check_rational_pf_riordan_theorem(
        &rule,
        0,
        &seeds,
        &RealRootednessTheoremOptions {
            root_region: Some(polytool::real_rootedness::RealInterval {
                lower: None,
                lower_closed: false,
                upper: Some(BigRational::zero()),
                upper_closed: false,
            }),
            ..Default::default()
        },
    )
    .unwrap();
    assert!(!report.is_proved());
    let closed = check(&rule, &seeds);
    assert!(closed.is_proved());
    assert!(
        closed
            .conclusion
            .as_ref()
            .unwrap()
            .root_region
            .as_ref()
            .unwrap()
            .at(0)
            .unwrap()
            .upper_closed
    );
    closed.verify().unwrap();
    let requested = check_rational_pf_riordan_theorem(
        &rule,
        0,
        &seeds,
        &RealRootednessTheoremOptions {
            root_region: Some(polytool::real_rootedness::RealInterval::at_most(
                BigRational::zero(),
            )),
            ..Default::default()
        },
    )
    .unwrap();
    assert!(!requested.is_proved());
    requested.verify().unwrap();
}
