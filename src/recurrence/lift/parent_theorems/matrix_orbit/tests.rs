//! Exact matrix-orbit recognition and replay regressions.
use super::*;

fn request(id: &str) -> MatrixOrbitRequest {
    let entry = crate::oeis::by_id(id).unwrap();
    let (rule, seeds) = entry.recurrence_parts().unwrap();
    MatrixOrbitRequest {
        recurrence: RecurrenceJson::from_recurrence_rational(
            &rule,
            entry.recurrence_first_index,
            &seeds,
            None,
        ),
        state_system: None,
    }
}
#[test]
fn automatically_recognizes_requested_models() {
    for id in ["A120987", "A304429", "A330509", "A340620", "A245963"] {
        let report = check_matrix_orbit(&request(id)).unwrap();
        assert_eq!(report.status, "proved", "{id}: {}", report.reason);
        verify_matrix_orbit(&report).unwrap();
    }
}
#[test]
fn rejects_tampered_pencil_and_projection() {
    let mut report = check_matrix_orbit(&request("A245963")).unwrap();
    report.certificate.as_mut().unwrap().pencils.pop();
    assert!(verify_matrix_orbit(&report).is_err());
    let mut report = check_matrix_orbit(&request("A304429")).unwrap();
    report
        .certificate
        .as_mut()
        .unwrap()
        .transformed_system
        .projection[0] = "-1".into();
    assert!(verify_matrix_orbit(&report).is_err());
}

#[test]
fn support_cycles_certify_nonvanishing_and_detect_later_zeros() {
    let one = Poly::one();
    let zero = Poly::zero();
    let cycle = vec![
        vec![zero.clone(), one.clone()],
        vec![one.clone(), zero.clone()],
    ];
    let projection = vec![Q::one(), Q::zero()];
    assert_eq!(
        nonzero_supports(&cycle, &[one.clone(), one.clone()], &projection),
        Some(vec![3])
    );
    assert_eq!(
        nonzero_supports(&cycle, &[one.clone(), zero.clone()], &projection),
        None
    );
    let dies = vec![
        vec![one.clone(), zero.clone()],
        vec![zero.clone(), zero.clone()],
    ];
    assert!(nonzero_supports(&dies, &[one, Poly::one()], &[Q::zero(), Q::one()]).is_none());
    assert!(affine::certify(&vec![vec![zero]]).is_some());
}

#[test]
fn basis_recovery_preserves_target_and_requires_nonnegative_weights() {
    let matrix = vec![
        vec![Poly::one(), Poly::zero()],
        vec![Poly::zero(), Poly::one()],
    ];
    let seed = vec![Poly::one(), Poly::one()];
    let basis = vec![vec![Q::one(), Q::one()], vec![Q::zero(), Q::one()]];
    assert!(realization::certify_basis(
        &matrix,
        &seed,
        &[Q::one(), Q::zero()],
        &[0, 1],
        &basis,
        "supplied",
        &[]
    )
    .is_none());
    let certificate = realization::certify_basis(
        &matrix,
        &seed,
        &[Q::one(), Q::one()],
        &[0, 1],
        &basis,
        "supplied",
        &[],
    )
    .unwrap();
    assert_eq!(certificate.transformed_system.projection, vec!["1", "0"]);
    assert_eq!(
        certificate.transformed_system.start_vector,
        vec![vec!["2"], vec!["1"]]
    );
}

#[test]
fn seed_checks_retain_order_signs_and_coordinate_real_rootedness() {
    let p = |c: &[i64]| Poly::new(c.iter().map(|&v| Q::from_integer(v.into())).collect());
    assert!(interlacing_seed(&[p(&[2, 1]), p(&[1, 1]), Poly::zero()]));
    assert!(!interlacing_seed(&[p(&[1, 1]), p(&[2, 1])]));
    assert!(!interlacing_seed(&[Poly::zero(), p(&[1, 0, 1])]));
    assert!(!interlacing_seed(&[p(&[-1, 1])]));
    assert!(affine::certify(&vec![vec![p(&[1, 0, 1])]]).is_none());
}

#[test]
fn theorem_replay_retains_shared_left_hand_factor_normalization() {
    let (mut rule, first, seeds) = request("A245963").recurrence.to_recurrence_parts().unwrap();
    for term in &mut rule.terms {
        term.coeff = BivarPoly::from_coefficients_normalized(
            term.coeff
                .normalized_coefficients()
                .iter()
                .map(|row| vec_poly::shift_by(row, 2))
                .collect(),
        );
    }
    rule.denominator = Some(BivarPoly::from_coefficients_normalized(vec![vec![
        Q::zero(),
        Q::zero(),
        Q::one(),
    ]]));
    let report = check_matrix_orbit_theorem(&rule, first, &seeds, &Default::default()).unwrap();
    assert!(report.is_proved());
    assert!(report.left_hand_factor_normalization.is_some());
    report.verify().unwrap();
}
