//! Regression tests for exact two-family products.
use super::*;

fn catalog(id: &str) -> RecurrenceJson {
    let entry = crate::oeis::by_id(id).unwrap();
    let (rule, seeds) = entry.recurrence_parts().unwrap();
    RecurrenceJson::from_recurrence_rational(&rule, entry.recurrence_first_index, &seeds, None)
}

#[test]
fn narayana_hadamard_a350512_proves_a174148() {
    // The catalog's excluded finite prefix is PF as well: 1 and 1+x.
    assert_eq!(
        crate::oeis::by_id("A174148").unwrap().prefix_rows,
        &[&["1"][..], &["1", "1"][..]]
    );
    let candidate = PolynomialFamilyProductCandidate {
        left: catalog("A001263"),
        right: catalog("A350512"),
        left_index_shift: 1,
        right_index_shift: 1,
        operation: PolynomialFamilyProductOperation::Hadamard,
    };
    let proof =
        certify_polynomial_family_product(&catalog("A174148"), &candidate, &Default::default())
            .unwrap();
    assert_eq!(proof.conclusion, "pf");
    proof.verify().unwrap();
    let mut tampered = proof;
    let PolynomialFamilyProductClosure::CoefficientMultiplier { certificate, .. } =
        &mut tampered.closure
    else {
        panic!("expected the certified Narayana multiplier closure");
    };
    certificate.gauge = VeroneseSectionGauge::None;
    assert!(tampered.verify().is_err());
}

#[test]
fn consecutive_fibonacci_row_product_proves_a092879() {
    let candidate = PolynomialFamilyProductCandidate {
        left: catalog("A011973"),
        right: catalog("A011973"),
        left_index_shift: 0,
        right_index_shift: 1,
        operation: PolynomialFamilyProductOperation::Ordinary,
    };
    let proof =
        certify_polynomial_family_product(&catalog("A092879"), &candidate, &Default::default())
            .unwrap();
    proof.verify().unwrap();
}

fn named(name: &str) -> RecurrenceJson {
    PreserverParentLibrary::named_family_records()
        .into_iter()
        .find(|(key, _)| key == name)
        .unwrap()
        .1
}

#[test]
fn binomial_pair_all_weights_forward_and_tamper_rejection() {
    for operation in [
        PolynomialFamilyProductOperation::Hadamard,
        PolynomialFamilyProductOperation::FactorialHadamard,
        PolynomialFamilyProductOperation::InverseFactorialHadamard,
    ] {
        let candidate = PolynomialFamilyProductCandidate {
            left: named("binomial"),
            right: named("binomial"),
            left_index_shift: 0,
            right_index_shift: 0,
            operation,
        };
        let proof = apply_polynomial_family_product(&candidate, 0, &Default::default()).unwrap();
        assert_eq!(proof.conclusion, "pf");
        let json = serde_json::to_string(&proof).unwrap();
        let decoded: PolynomialFamilyProductCertificate = serde_json::from_str(&json).unwrap();
        decoded.verify().unwrap();
        let mut changed = decoded.clone();
        changed.theorem.push_str(" tampered");
        assert!(changed.verify().is_err());
        let mut changed = decoded;
        changed.candidate.left_index_shift = 1;
        assert!(changed.verify().is_err());
    }
}

#[test]
fn growing_eulerian_binomial_needs_a_runnable_image_rule() {
    let candidate = PolynomialFamilyProductCandidate {
        left: named("eulerian"),
        right: named("binomial"),
        left_index_shift: 0,
        right_index_shift: 0,
        operation: PolynomialFamilyProductOperation::Hadamard,
    };
    let error = apply_polynomial_family_product(&candidate, 1, &Default::default()).unwrap_err();
    assert!(error.to_string().contains("no product image rule"));
}

#[test]
fn schur_accepts_whole_line_left_but_unweighted_requires_pf() {
    let mut candidate = PolynomialFamilyProductCandidate {
        left: named("hermite_he"),
        right: named("binomial"),
        left_index_shift: 0,
        right_index_shift: 0,
        operation: PolynomialFamilyProductOperation::Hadamard,
    };
    assert!(apply_polynomial_family_product(&candidate, 0, &Default::default()).is_err());
    candidate.operation = PolynomialFamilyProductOperation::FactorialHadamard;
    let proof = apply_polynomial_family_product(&candidate, 0, &Default::default()).unwrap();
    assert_eq!(proof.conclusion, "zero_or_real_rooted");
    proof.verify().unwrap();
    candidate.operation = PolynomialFamilyProductOperation::InverseFactorialHadamard;
    assert!(apply_polynomial_family_product(&candidate, 0, &Default::default()).is_err());
}

#[test]
fn negative_sign_is_not_a_pf_certificate() {
    let (rule, first, mut seeds) = named("binomial").to_recurrence_parts().unwrap();
    seeds[0][0] = -Q::one();
    let candidate = PolynomialFamilyProductCandidate {
        left: RecurrenceJson::from_recurrence_rational(&rule, first, &seeds, None),
        right: named("binomial"),
        left_index_shift: 0,
        right_index_shift: 0,
        operation: PolynomialFamilyProductOperation::Hadamard,
    };
    assert!(apply_polynomial_family_product(&candidate, 0, &Default::default()).is_err());
}

fn scalar(value: i64) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(vec![vec![Q::from_integer(value.into())]])
}

#[test]
fn fixed_eulerian_polynomial_hadamard_binomial_forward() {
    // The fixed PF family A_n=A_3=1+4x+x^2; the other family is (1+x)^n.
    let rule = Recurrence {
        terms: vec![RecurrenceTerm::new(1, 0, scalar(1))],
        denominator: None,
        inhomogeneous: None,
    };
    let rows = vec![vec![Q::one(), Q::from_integer(4.into()), Q::one()]];
    let candidate = PolynomialFamilyProductCandidate {
        left: RecurrenceJson::from_recurrence_rational(&rule, 0, &rows, None),
        right: named("binomial"),
        left_index_shift: 0,
        right_index_shift: 0,
        operation: PolynomialFamilyProductOperation::Hadamard,
    };
    let proof = apply_polynomial_family_product(&candidate, 0, &Default::default()).unwrap();
    proof.verify().unwrap();
    assert_eq!(proof.conclusion, "pf");
}

#[test]
fn two_non_menu_families_coefficient_closure_and_tampering() {
    // A_n=B_n=1+(n+1)x, so Q_n=1+(n+1)^2 x. The source and image
    // are defined by second and third finite differences respectively.
    let rows: Vec<Vec<Q>> = (1..=3)
        .map(|n| vec![Q::one(), Q::from_integer((n * n).into())])
        .collect();
    let source_rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, scalar(2)),
            RecurrenceTerm::new(2, 0, scalar(-1)),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let source_rows: Vec<Vec<Q>> = (1..=2)
        .map(|n| vec![Q::one(), Q::from_integer(n.into())])
        .collect();
    let source = RecurrenceJson::from_recurrence_rational(&source_rule, 0, &source_rows, None);
    let image_rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, scalar(3)),
            RecurrenceTerm::new(2, 0, scalar(-3)),
            RecurrenceTerm::new(3, 0, scalar(1)),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let target = RecurrenceJson::from_recurrence_rational(&image_rule, 0, &rows, None);
    let candidate = PolynomialFamilyProductCandidate {
        left: source.clone(),
        right: source,
        left_index_shift: 0,
        right_index_shift: 0,
        operation: PolynomialFamilyProductOperation::Hadamard,
    };
    let proof =
        certify_polynomial_family_product(&target, &candidate, &Default::default()).unwrap();
    assert!(matches!(
        proof.closure,
        PolynomialFamilyProductClosure::CoefficientRelations { .. }
    ));
    proof.verify().unwrap();
    let mut changed = proof.clone();
    if let PolynomialFamilyProductClosure::CoefficientRelations { image, .. } = &mut changed.closure
    {
        image.system.generators[0].step_quotients.clear();
    }
    assert!(changed.verify().is_err());
    let mut changed = proof;
    changed.conclusion = "strictly_interlacing".into();
    assert!(changed.verify().is_err());
}

#[test]
fn forward_fibonacci_product_equals_target_by_recurrence_equality() {
    let candidate = PolynomialFamilyProductCandidate {
        left: catalog("A011973"),
        right: catalog("A011973"),
        left_index_shift: 0,
        right_index_shift: 1,
        operation: PolynomialFamilyProductOperation::Ordinary,
    };
    let proof = apply_polynomial_family_product(&candidate, 0, &Default::default()).unwrap();
    let equality = certify_recurrence_equality(&catalog("A092879"), &proof.image).unwrap();
    assert!(equality.is_certified());
    equality.verify().unwrap();
    let mut changed = proof;
    if let PolynomialFamilyProductClosure::Ordinary { certificate } = &mut changed.closure {
        certificate.quotients.clear();
    }
    assert!(changed.verify().is_err());
}

#[test]
fn forward_narayana_product_equals_a174148() {
    let candidate = PolynomialFamilyProductCandidate {
        left: catalog("A001263"),
        right: catalog("A350512"),
        left_index_shift: 1,
        right_index_shift: 1,
        operation: PolynomialFamilyProductOperation::Hadamard,
    };
    let proof = apply_polynomial_family_product(&candidate, 1, &Default::default()).unwrap();
    let equality = certify_recurrence_equality(&catalog("A174148"), &proof.image).unwrap();
    assert!(equality.is_certified());
    proof.verify().unwrap();
    equality.verify().unwrap();
}

#[test]
fn pf_hypothesis_cannot_be_inferred_from_a_long_positive_prefix() {
    // All sampled rows are PF, but at n=201 the root becomes positive.
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, scalar(2)),
            RecurrenceTerm::new(2, 0, scalar(-1)),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let rows: Vec<Vec<Q>> = [200, 199]
        .into_iter()
        .map(|c| vec![Q::one(), Q::from_integer(c.into())])
        .collect();
    let candidate = PolynomialFamilyProductCandidate {
        left: RecurrenceJson::from_recurrence_rational(&rule, 0, &rows, None),
        right: named("binomial"),
        left_index_shift: 0,
        right_index_shift: 0,
        operation: PolynomialFamilyProductOperation::Hadamard,
    };
    assert!(apply_polynomial_family_product(&candidate, 0, &Default::default()).is_err());
}
