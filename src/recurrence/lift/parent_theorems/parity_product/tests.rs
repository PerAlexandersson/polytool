//! Exact cubic, allocation, replay and negative tests.
use super::*;

fn p(c: &[i64]) -> Poly {
    Poly::new(c.iter().map(|v| Q::from_integer((*v).into())).collect())
}
fn rule(mu: &Poly, beta: &Poly) -> Recurrence {
    let a = beta.clone() + mu.clone();
    Recurrence {
        terms: [a.clone(), -(a * mu.clone()), mu.pow(3)]
            .into_iter()
            .enumerate()
            .filter(|(_, p)| !p.is_zero())
            .map(|(i, p)| {
                RecurrenceTerm::new(
                    i + 1,
                    0,
                    BivarPoly::from_coefficients_normalized(vec![p.coeffs().to_vec()]),
                )
            })
            .collect(),
        denominator: None,
        inhomogeneous: None,
    }
}
fn rows(ps: &[Poly]) -> Vec<Vec<Q>> {
    ps.iter().map(|p| vec_poly::trim(p.coeffs())).collect()
}
fn example() -> (Recurrence, Vec<Vec<Q>>) {
    (
        rule(&p(&[-1]), &p(&[2, 1])),
        rows(&[p(&[1]), p(&[1, 1]), p(&[2, 3, 1])]),
    )
}
#[test]
fn exact_cubic_recognition_including_zero_trace_and_wrong_cofactor() {
    for (mu, beta) in [
        (p(&[-1]), p(&[2, 1])),
        (p(&[0, 1]), p(&[1, 2])),
        (p(&[2, 3, 1]), p(&[-2, -3, -1])),
        (
            p(&[0, 1]).scale(&Q::new(1.into(), 2.into())),
            p(&[1, 2]).scale(&Q::new(1.into(), 3.into())),
        ),
        (
            p(&[1, 1]).scale(&Q::new(1.into(), 2.into())),
            p(&[-1, -1]).scale(&Q::new(1.into(), 2.into())),
        ),
    ] {
        let r = rule(&mu, &beta);
        let m = constant_coefficient_multipliers(&r).unwrap();
        assert_eq!(recognize(&m), Some((mu, beta)));
        let mut wrong = m;
        wrong[1] = wrong[1].clone() + p(&[1]);
        assert!(recognize(&wrong).is_none());
    }
    assert!(recognize(&[p(&[0]), p(&[0]), p(&[2])]).is_none());
}
#[test]
fn proves_products_and_replays_without_claiming_strictness() {
    let (r, seeds) = example();
    let report = check_parity_product_theorem(&r, 0, &seeds, &Default::default()).unwrap();
    assert!(report.is_proved(), "{report:?}");
    report.verify().unwrap();
    let c = report.parity_product_certificate().unwrap();
    c.verify(&report).unwrap();
    assert_eq!(c.factor_proofs.len(), 2);
    assert_eq!(c.proof_depth, 1);
    let conclusion = report.conclusion.as_ref().unwrap();
    assert_eq!(conclusion.interlacing, Some(RowInterlacing::Weak));
    assert_eq!(conclusion.interlacing_from, Some(1));
    assert!(!conclusion.simple_roots);
    let mut changed = report.clone();
    if let Some(RouteCertificate::ParityProduct(c)) = &mut changed.certificate {
        c.mu = vec!["1".into()];
    }
    assert!(changed.verify().is_err());
    let generated = r.generate_rows_rational(&seeds, 0, 12).unwrap();
    let (sr, _, ss) = c.factor_proofs[0].recurrence.to_recurrence_parts().unwrap();
    let (tr, _, ts) = c.factor_proofs[1].recurrence.to_recurrence_parts().unwrap();
    let s = sr.generate_rows_rational(&ss, 0, 7).unwrap();
    let t = tr.generate_rows_rational(&ts, 0, 7).unwrap();
    for m in 0..6 {
        assert_eq!(
            Poly::new(generated[2 * m].clone()),
            Poly::new(s[m].clone()) * Poly::new(t[m].clone())
        );
        assert_eq!(
            Poly::new(generated[2 * m + 1].clone()),
            Poly::new(s[m + 1].clone()) * Poly::new(t[m].clone())
        );
    }
}
#[test]
fn nontrivial_gcd_allocation_is_found() {
    // Start two rows later: gcd(P2,P3)=x+2, whose constant divisor fails.
    let (r, seeds) = example();
    let generated = r.generate_rows_rational(&seeds, 0, 5).unwrap();
    let report = check_parity_product_theorem(&r, 2, &generated[2..], &Default::default()).unwrap();
    assert!(report.is_proved());
    let c = report.parity_product_certificate().unwrap();
    assert_eq!(c.initial_common_divisor, vec!["2", "1"]);
    assert_eq!(c.divisors_examined, 2);
    assert_eq!(c.divisor_lattice_size, 2);
    report.verify().unwrap();
}
#[test]
fn casoratian_failure_and_two_nonreal_factors_do_not_prove() {
    let (r, mut seeds) = example();
    seeds[2] = p(&[3, 3, 1]).coeffs().to_vec();
    let report = check_parity_product_theorem(&r, 0, &seeds, &Default::default()).unwrap();
    assert!(!report.is_proved());
    assert!(report.normal_form.contains("Casoratian"));
    // The cubic and Casoratian are exact, but s1=x^2+1 and t1=x^2+2 are nonreal.
    let beta = p(&[2, 0, 1]);
    let r = rule(&p(&[-1]), &beta);
    let seeds = rows(&[p(&[1]), p(&[1, 0, 1]), p(&[1, 0, 1]) * beta]);
    let report = check_parity_product_theorem(
        &r,
        0,
        &seeds,
        &RealRootednessTheoremOptions {
            maximum_start_delay: 0,
            root_region: None,
        },
    )
    .unwrap();
    assert!(!report.is_proved());
}
#[test]
fn reversal_and_constant_row_normalization_are_exact() {
    let (r, seeds) = example();
    let record = DegreeLawTransform::Reversal {
        degree_slope: 1,
        degree_offset: 0,
    }
    .apply(&r, 0, &seeds)
    .unwrap();
    let (r, first, seeds) = record.to_recurrence_parts().unwrap();
    let report = check_parity_product_theorem(&r, first, &seeds, &Default::default()).unwrap();
    assert!(report.is_proved());
    report.verify().unwrap();
    let mut scaled = r.clone();
    scaled.denominator = Some(BivarPoly::from_coefficients_normalized(vec![vec![
        Q::from_integer(7.into()),
    ]]));
    for t in &mut scaled.terms {
        for row in &mut t.coeff.coeffs {
            for c in row {
                *c *= Q::from_integer(7.into());
            }
        }
    }
    assert!(
        check_parity_product_theorem(&scaled, first, &seeds, &Default::default())
            .unwrap()
            .is_proved()
    );
}

#[test]
fn recording_preserves_existing_conclusions_and_rejects_excessive_bounds() {
    let (r, seeds) = example();
    let report = check_parity_product_theorem(&r, 0, &seeds, &Default::default()).unwrap();
    assert!(!DiscoveryPolicy::WhenStronger.records_after(&report, std::slice::from_ref(&report)));
    let mut weaker = report.clone();
    weaker.conclusion.as_mut().unwrap().interlacing = None;
    weaker.conclusion.as_mut().unwrap().interlacing_from = None;
    assert!(DiscoveryPolicy::WhenStronger.records_after(&report, &[weaker]));
    assert!(check_parity_product_theorem(
        &r,
        0,
        &seeds,
        &RealRootednessTheoremOptions {
            maximum_start_delay: 17,
            root_region: None
        }
    )
    .is_err());
    let invalid = rows(&[p(&[1, 0, 1]), p(&[1]), p(&[1]), p(&[1, 1]), p(&[2, 3, 1])]);
    assert!(
        !check_parity_product_theorem(&r, 0, &invalid, &Default::default())
            .unwrap()
            .is_proved()
    );
}
