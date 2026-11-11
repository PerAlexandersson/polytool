//! Exact products of first-order preservers, including the insertion operator
//! for Jacobi--Stirling descents (Ma--Wang, arXiv:2610.03111, Section 5).
//! Fixtures are composed from their defining operators, never fitted to rows.

use num_bigint::BigInt;
use num_traits::{One, Zero};
use ore_algebra::{BivariatePolynomial as Biv, OreOperator};
use polytool::recurrence::{
    check_iterated_preserver_theorem, BivarPoly, RealRootednessTheoremOptions, Recurrence,
    RecurrenceTerm, RouteCertificate,
};
use polytool::BigRational as Q;

fn q(n: i64) -> Q {
    Q::from_integer(BigInt::from(n))
}

fn insertion(slope: i64, intercept: i64) -> OreOperator {
    let mut op = OreOperator::single(
        0,
        0,
        Biv::new(vec![vec![q(1), q(intercept)], vec![q(0), q(slope)]]),
    );
    op.add_term((0, 1), Biv::from_variable_polynomial(&[q(0), q(1), q(-1)]));
    op
}

fn rule(op: &OreOperator) -> Recurrence {
    Recurrence {
        terms: op
            .0
            .iter()
            .map(|((_, d), c)| {
                RecurrenceTerm::new(
                    1,
                    *d,
                    BivarPoly::from_coefficients_normalized(c.rows().to_vec()),
                )
            })
            .collect(),
        denominator: None,
        inhomogeneous: None,
    }
}

fn options() -> RealRootednessTheoremOptions {
    RealRootednessTheoremOptions {
        maximum_start_delay: 0,
        ..Default::default()
    }
}

#[test]
fn jacobi_stirling_insertion_product_is_proved_and_replayed() {
    let op = insertion(3, -2).compose(&insertion(3, -3));
    let r = rule(&op);
    let fixture: polytool::recurrence::RecurrenceJson = serde_json::from_str(include_str!(
        "fixtures/operator-factorization/jacobi-stirling-i0.json"
    ))
    .unwrap();
    let fixture_rule = fixture.to_recurrence_parts().unwrap().0;
    let mut fixture_op = OreOperator::default();
    for term in fixture_rule.terms {
        fixture_op.add_term((0, term.deriv_order), Biv::new(term.coeff.coeffs));
    }
    assert_eq!(fixture_op, op);
    // At n=1 the nominal degree-two leading term cancels. The standard
    // moved-start check generates A_1 exactly, then applies the product.
    let from_zero = check_iterated_preserver_theorem(
        &r,
        0,
        &[vec![q(1)]],
        &RealRootednessTheoremOptions::default(),
    )
    .unwrap();
    assert!(from_zero.is_proved(), "{from_zero:#?}");
    from_zero.verify().unwrap();
    let report = check_iterated_preserver_theorem(&r, 1, &[vec![q(1), q(1)]], &options()).unwrap();
    assert!(report.is_proved(), "{report:#?}");
    let certificate = report.iterated_preserver_certificate().unwrap();
    assert_eq!(certificate.factors.len(), 2);
    assert_eq!(certificate.degree_step, 2);
    assert_eq!(
        report
            .conclusion
            .as_ref()
            .unwrap()
            .root_region
            .as_ref()
            .unwrap()
            .to_string(),
        "(-inf, 0]"
    );
    assert!(!report.conclusion.as_ref().unwrap().simple_roots);
    assert!(report.conclusion.as_ref().unwrap().interlacing.is_none());
    certificate.verify_operator_identity(&r).unwrap();
    report.verify().unwrap();
    let json = serde_json::to_string(&report).unwrap();
    let decoded: polytool::recurrence::RealRootednessTheoremReport =
        serde_json::from_str(&json).unwrap();
    decoded.verify().unwrap();
    // Two independently published rows from Table 1.
    let row2 = op.apply_row(&[q(1), q(1)], 2);
    assert_eq!(row2, vec![q(1), q(14), q(21), q(4)]);
    assert_eq!(
        op.apply_row(&row2, 3),
        vec![q(1), q(75), q(603), q(1065), q(460), q(36)]
    );
    let mut tampered = report.clone();
    let Some(RouteCertificate::IteratedPreserver(c)) = tampered.certificate.as_mut() else {
        panic!()
    };
    c.factors[0].multiplication[0].numerator[0] = "37".into();
    assert!(c.verify_operator_identity(&r).is_err());
    assert!(tampered.verify().is_err());
    let mut reordered = certificate.clone();
    reordered.factors.swap(0, 1);
    assert!(reordered.verify_operator_identity(&r).is_err());
    let mut different_rule = r.clone();
    different_rule.terms[0].coeff = BivarPoly::from_coefficients_normalized(vec![vec![q(2)]]);
    assert!(certificate
        .verify_operator_identity(&different_rule)
        .is_err());
}

#[test]
fn synthetic_affine_insertions_and_closed_endpoint() {
    for (a, b, c, d) in [(3, 0, 4, 1), (2, 3, 2, 2), (5, -1, 3, 0)] {
        let r = rule(&insertion(a, b).compose(&insertion(c, d)));
        // Including an endpoint root, handled by continuity of the closed ray.
        for seed in [vec![Q::one()], vec![Q::zero(), Q::one()]] {
            let report = check_iterated_preserver_theorem(&r, 0, &[seed], &options()).unwrap();
            assert!(report.is_proved(), "{a} {b} {c} {d}: {report:#?}");
            report.verify().unwrap();
        }
    }
}

#[test]
fn bad_factor_and_irreducible_airy_step_are_not_proofs() {
    // (1+x-D) o (1+x+D): applied to 1 gives x^2+2x,
    // but to that gives x^4+4x^3+4x^2-2, which has nonreal roots.
    let a = OreOperator::single(0, 0, Biv::from_variable_polynomial(&[q(1), q(1)]));
    let derivative = OreOperator::single(0, 1, Biv::one());
    let r = rule(&a.subtract(&derivative).compose(&a.add(&derivative)));
    let report = check_iterated_preserver_theorem(&r, 0, &[vec![q(1)]], &options()).unwrap();
    assert!(!report.is_proved(), "{report:#?}");
    let airy = OreOperator::single(0, 2, Biv::one()).add(&OreOperator::single(
        0,
        0,
        Biv::from_variable_polynomial(&[q(0), q(1)]),
    ));
    let report =
        check_iterated_preserver_theorem(&rule(&airy), 0, &[vec![q(1)]], &options()).unwrap();
    assert!(!report.is_proved());
    assert!(report
        .iterated_preserver_certificate()
        .is_none_or(|c| c.factors.is_empty()));
}

#[test]
fn distinct_derivative_coefficients_index_scalars_and_outer_multiplier() {
    let mut right = OreOperator::single(0, 0, Biv::new(vec![vec![q(2), q(2)], vec![q(0), q(10)]]));
    right.add_term((0, 1), Biv::from_variable_polynomial(&[q(0), q(2), q(-2)]));
    right = right.left_multiply(&Biv::from_index_polynomial(&[q(1), q(1)]));
    let mut left = OreOperator::single(0, 0, Biv::new(vec![vec![q(3), q(5)], vec![q(0), q(25)]]));
    left.add_term((0, 1), Biv::from_variable_polynomial(&[q(0), q(3), q(-6)]));
    for multiplier in [
        Biv::from_variable_polynomial(&[q(1), q(1)]),
        Biv::new(vec![vec![q(1), q(1)], vec![q(0), q(1)]]),
    ] {
        let op = left.compose(&right).left_multiply(&multiplier);
        let mut r = rule(&op);
        r.denominator = Some(BivarPoly::from_coefficients_normalized(vec![
            vec![q(2)],
            vec![q(1)],
        ]));
        let report = check_iterated_preserver_theorem(&r, 0, &[vec![q(1)]], &options()).unwrap();
        assert!(report.is_proved(), "{report:#?}");
        assert_eq!(
            report
                .iterated_preserver_certificate()
                .unwrap()
                .factors
                .len(),
            2
        );
        assert_eq!(
            report.iterated_preserver_certificate().unwrap().degree_step,
            3
        );
        report.verify().unwrap();
    }
}

#[test]
fn reflection_keeps_the_positive_ray_and_bad_seed_is_rejected() {
    let mut op = insertion(3, -2).compose(&insertion(3, -3));
    for ((_, d), c) in &mut op.0 {
        let mut rows = c.rows().to_vec();
        for row in &mut rows {
            for (j, v) in row.iter_mut().enumerate() {
                if (j + *d) % 2 == 1 {
                    *v = -v.clone();
                }
            }
        }
        *c = Biv::new(rows);
    }
    let r = rule(&op);
    let report = check_iterated_preserver_theorem(&r, 1, &[vec![q(1), q(-1)]], &options()).unwrap();
    assert!(report.is_proved(), "{report:#?}");
    assert_eq!(
        report
            .conclusion
            .as_ref()
            .unwrap()
            .root_region
            .as_ref()
            .unwrap()
            .to_string(),
        "[0, +inf)"
    );
    report.verify().unwrap();
    let bad =
        check_iterated_preserver_theorem(&r, 1, &[vec![q(1), q(0), q(1)]], &options()).unwrap();
    assert!(!bad.is_proved());
}

#[test]
fn rational_index_coefficients_are_defined_on_the_entire_rule_range() {
    let mut right = OreOperator::single(
        0,
        0,
        Biv::new(vec![vec![q(1)], vec![q(0), q(4)], vec![q(0), q(4)]]),
    );
    right.add_term(
        (0, 1),
        Biv::new(vec![vec![q(0), q(1), q(-1)], vec![q(0), q(1), q(-1)]]),
    );
    let mut r = rule(&insertion(3, 1).compose(&right));
    let report = check_iterated_preserver_theorem(&r, 0, &[vec![q(1)]], &options()).unwrap();
    assert!(report.is_proved(), "{report:#?}");
    let c = report.iterated_preserver_certificate().unwrap();
    assert_eq!(c.factors.len(), 2);
    assert!(c
        .factors
        .iter()
        .flat_map(|f| f.multiplication.iter().chain(&f.derivative))
        .any(|v| v.denominator.len() > 1));
    report.verify().unwrap();
    r.denominator = Some(BivarPoly::from_coefficients_normalized(vec![
        vec![q(-1)],
        vec![q(1)],
    ]));
    let bad = check_iterated_preserver_theorem(&r, 0, &[vec![q(1)]], &options()).unwrap();
    assert!(!bad.is_proved());
}
