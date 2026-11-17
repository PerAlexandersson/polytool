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

#[test]
fn independent_sections_prove_catalog_models_and_replay_allocations() {
    for id in ["A082649", "A156308", "A158454", "A334009"] {
        let entry = crate::oeis::by_id(id).unwrap();
        let (rule, seeds) = entry.recurrence_parts().unwrap();
        let report = check_parity_product_theorem(
            &rule,
            entry.recurrence_first_index,
            &seeds,
            &Default::default(),
        )
        .unwrap();
        assert!(report.is_proved(), "{id}: {report:?}");
        let certificate = report.parity_product_certificate().unwrap();
        assert_eq!(certificate.factor_proofs.len(), 4);
        assert_eq!(certificate.independent_sections.as_ref().unwrap().len(), 2);
        assert!(report.conclusion.as_ref().unwrap().interlacing.is_none());
        assert!(!report.conclusion.as_ref().unwrap().simple_roots);
        report.verify().unwrap();
        certificate.verify(&report).unwrap();
        for field in 0..5 {
            let mut changed = report.clone();
            let Some(RouteCertificate::ParityProduct(c)) = changed.certificate.as_mut() else {
                unreachable!()
            };
            let section = &mut c.independent_sections.as_mut().unwrap()[0];
            match field {
                0 => section.initial_scalar_ratio = "1234567/2".into(),
                1 => section.seed_products[2] = vec!["123".into()],
                2 => section.subsequence.first_source_index += 1,
                3 => section.initial_divisors[0] = vec!["123".into()],
                _ => section.product_characteristic_polynomial[0] = vec!["123".into()],
            }
            assert!(changed.verify().is_err(), "{id}, altered field {field}");
        }
    }
}

#[test]
fn independent_sections_use_whole_line_liu_wang_factors() {
    let (r, seeds) = (
        rule(&p(&[1]), &p(&[0, 2])),
        rows(&[p(&[1]), p(&[1, 1]), p(&[0, 1, 2])]),
    );
    let report = check_parity_product_theorem(&r, 0, &seeds, &Default::default()).unwrap();
    assert!(report.is_proved(), "{report:?}");
    assert!(report
        .parity_product_certificate()
        .unwrap()
        .factor_proofs
        .iter()
        .any(|r| r.theorem == RealRootednessTheorem::LiuWangThreeTerm));
    report.verify().unwrap();
    // This proved family itself refutes an adjacent-interlacing transfer:
    // P1=x+1 has its root outside the two roots of P2=x*(2x+1).
    assert_eq!(
        crate::check_weak_interlacing_bigint_coeffs(
            &[1.into(), 1.into()],
            &[0.into(), 1.into(), 2.into()],
        ),
        Some(false)
    );
    let original = r.generate_rows_rational(&seeds, 0, 16).unwrap();
    let c = report.parity_product_certificate().unwrap();
    for (residue, section) in c.independent_sections.as_ref().unwrap().iter().enumerate() {
        let generated: Vec<_> = c.factor_proofs[2 * residue..2 * residue + 2]
            .iter()
            .map(|proof| {
                let (r, first, seeds) = proof.recurrence.to_recurrence_parts().unwrap();
                r.generate_rows_rational(&seeds, first, 7).unwrap()
            })
            .collect();
        for (m, (left, right)) in generated[0].iter().zip(&generated[1]).enumerate() {
            let source = section.subsequence.first_source_index + 2 * (section.skipped_rows + m);
            assert_eq!(
                Poly::new(original[source].clone()),
                Poly::new(left.clone()) * Poly::new(right.clone())
            );
        }
    }
}

#[test]
fn failed_casoratian_with_zero_seed_replays_for_each_search_bound() {
    let r = rule(&p(&[1]), &p(&[-2, 0, 4]));
    let seeds = rows(&[p(&[0]), p(&[1, 0, -1]), p(&[0, 0, 4, 0, -4])]);
    for bound in [0, 1, 3] {
        let report = check_parity_product_theorem(
            &r,
            0,
            &seeds,
            &RealRootednessTheoremOptions {
                maximum_start_delay: bound,
                root_region: None,
            },
        )
        .unwrap();
        assert!(!report.is_proved());
        report.verify().unwrap();
    }
    // Starting after the zero row permits the existing index-section lift
    // to prove the factors of the quadratic marker substitution.
    let generated = r.generate_rows_rational(&seeds, 0, 4).unwrap();
    let report = check_parity_product_theorem(&r, 1, &generated[1..], &Default::default()).unwrap();
    assert!(report.is_proved(), "{report:?}");
    assert!(report
        .parity_product_certificate()
        .unwrap()
        .factor_proofs
        .iter()
        .any(|p| p.theorem == RealRootednessTheorem::IndexSectionLift));
    report.verify().unwrap();
}

#[test]
fn shared_success_keeps_legacy_certificate_serialization() {
    let (r, seeds) = example();
    let report = check_parity_product_theorem(&r, 0, &seeds, &Default::default()).unwrap();
    let c = report.parity_product_certificate().unwrap();
    assert!(c.independent_sections.is_none());
    assert!(!serde_json::to_value(c)
        .unwrap()
        .as_object()
        .unwrap()
        .contains_key("independent_sections"));
    report.verify().unwrap();
}

#[test]
fn independent_sections_allow_repeated_roots_and_scalar_normalization() {
    let r = rule(&p(&[1]), &p(&[2]));
    let scalar_seeds = rows(&[p(&[1]), p(&[4]), p(&[9])]);
    let scalar = check_parity_product_theorem(&r, 0, &scalar_seeds, &Default::default()).unwrap();
    assert!(!scalar.is_proved());
    scalar.verify().unwrap();
    // Constant beta and repeated characteristic roots remain supported for
    // polynomial seed families: P_n=(x+n+1)^2 has characteristic (z-1)^3.
    let seeds = rows(&[p(&[1, 2, 1]), p(&[4, 4, 1]), p(&[9, 6, 1])]);
    let report = check_parity_product_theorem(&r, 0, &seeds, &Default::default()).unwrap();
    assert!(report.is_proved(), "{report:?}");
    assert!(report
        .parity_product_certificate()
        .unwrap()
        .independent_sections
        .is_some());
    report.verify().unwrap();

    let entry = crate::oeis::by_id("A082649").unwrap();
    let (r, seeds) = entry.recurrence_parts().unwrap();
    let record = DegreeLawTransform::Reversal {
        degree_slope: 1,
        degree_offset: 0,
    }
    .apply(&r, entry.recurrence_first_index, &seeds)
    .unwrap();
    let (mut r, first, seeds) = record.to_recurrence_parts().unwrap();
    let seeds: Vec<_> = seeds
        .iter()
        .map(|p| vec_poly::scale(p, &Q::new((-3).into(), 2.into())))
        .collect();
    r.denominator = Some(BivarPoly::from_coefficients_normalized(vec![vec![
        Q::from_integer(7.into()),
    ]]));
    for term in &mut r.terms {
        for row in &mut term.coeff.coeffs {
            for c in row {
                *c *= Q::from_integer(7.into());
            }
        }
    }
    let report = check_parity_product_theorem(&r, first, &seeds, &Default::default()).unwrap();
    assert!(report.is_proved(), "{report:?}");
    report.verify().unwrap();
}
