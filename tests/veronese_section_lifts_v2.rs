//! End-to-end tests of `recurrence::certify_veronese_section_lift`: forward
//! certificates for index-dependent and gauged inputs, the parent provers,
//! near misses, budgets and certificate replay.

use polytool::recurrence::{
    certify_veronese_section_lift, find_veronese_section_lifts, BigRational, BivarPoly,
    LiftParentProof, LiftVerdict, Recurrence, RecurrenceError, RecurrenceJson, RecurrenceTerm,
    VeroneseParentAttestation, VeroneseParentProofRoute, VeroneseParentProofStatus,
    VeroneseSectionCandidate, VeroneseSectionCertification, VeroneseSectionCertifyOptions,
    VeroneseSectionGauge, VeroneseSectionIdentityStatus, VeroneseSectionLiftOptions,
    VeroneseSectionRecurrenceOptions, VeroneseSectionVerdict,
};

mod common;
use common::q;

fn row(values: &[i64]) -> Vec<BigRational> {
    values.iter().map(|&value| q(value)).collect()
}

fn rows(values: &[&[i64]]) -> Vec<Vec<BigRational>> {
    values.iter().map(|values| row(values)).collect()
}

/// `coeffs[i][j]` = coefficient of `n^i t^j`.
fn bivar(coeffs: &[&[i64]]) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(coeffs.iter().map(|r| row(r)).collect())
}

/// `f(n, t) P_n = sum c_k(n, t) D^d P_(n-k)` from `(k, d, c_k)`.
fn rule(terms: &[(usize, usize, &[&[i64]])], denominator: Option<&[&[i64]]>) -> Recurrence {
    Recurrence {
        terms: terms
            .iter()
            .map(|(offset, order, coeffs)| RecurrenceTerm::new(*offset, *order, bivar(coeffs)))
            .collect(),
        denominator: denominator.map(bivar),
        inhomogeneous: None,
    }
}

/// `P_n = c(x) P_(n-1)`, `P_0 = 1`.
fn product_parent(c: &[i64]) -> RecurrenceJson {
    RecurrenceJson::from_recurrence_rational(&rule(&[(1, 0, &[c])], None), 0, &rows(&[&[1]]), None)
}

/// `Q_(n+2) = 2 Q_(n+1) + (y - 1) Q_n`.
fn worked() -> Recurrence {
    rule(&[(1, 0, &[&[2]]), (2, 0, &[&[-1, 1]])], None)
}

fn candidate(modulus: usize, residue: usize, parent: RecurrenceJson) -> VeroneseSectionCandidate {
    VeroneseSectionCandidate {
        modulus,
        residue,
        parent,
        index_shift: 0,
        gauge: VeroneseSectionGauge::None,
    }
}

fn certify(
    input: &Recurrence,
    seeds: &[Vec<BigRational>],
    candidate: &VeroneseSectionCandidate,
) -> VeroneseSectionCertification {
    let certification = certify_veronese_section_lift(
        input,
        seeds,
        0,
        candidate,
        &VeroneseSectionCertifyOptions::default(),
    )
    .expect("certification runs");
    certification.verify().expect("the certification replays");
    certification
}

fn assert_replays(certification: &VeroneseSectionCertification) {
    let json = serde_json::to_string(certification).unwrap();
    let back: VeroneseSectionCertification = serde_json::from_str(&json).unwrap();
    assert_eq!(&back, certification);
    back.verify().unwrap();
    let value: serde_json::Value = serde_json::from_str(&json).unwrap();
    assert_eq!(
        value["schema"],
        "polytool.recurrence-veronese-certificate.v1"
    );
    assert_eq!(value["lean_replay"], "not_attempted");
}

#[test]
fn worked_example_is_proved_and_agrees_with_v1() {
    let seeds = rows(&[&[1], &[1]]);
    let certification = certify(&worked(), &seeds, &candidate(2, 0, product_parent(&[1, 1])));
    assert_eq!(certification.verdict, VeroneseSectionVerdict::Proved);
    assert_eq!(
        certification.parent_proof.route,
        Some(VeroneseParentProofRoute::PeriodicProduct)
    );
    let VeroneseSectionIdentityStatus::Certified { certificate, .. } = &certification.identity
    else {
        panic!("not certified: {:?}", certification.identity);
    };
    assert_eq!(certificate.multiplier, vec!["1".to_string()]);
    assert!(certificate.exceptional_bases.is_empty());
    assert_eq!(certificate.identifying_indices, vec![0, 1]);
    assert_replays(&certification);

    // The v1 search finds the same parent; its record certifies the same way.
    let report =
        find_veronese_section_lifts(&worked(), &seeds, 0, &VeroneseSectionLiftOptions::default())
            .unwrap();
    let lift = &report.lifts[0];
    assert_eq!((lift.modulus, lift.residue), (2, 0));
    let parent = RecurrenceJson::from_recurrence_rational(
        &lift.parent_recurrence,
        lift.parent_first_index,
        &lift.parent_initial_rows,
        None,
    );
    let from_v1 = certify(&worked(), &seeds, &candidate(2, 0, parent));
    assert_eq!(from_v1.verdict, VeroneseSectionVerdict::Proved);
}

/// `V_(2,0)((1+x)^n) / n!`: `m (m-1) Q_m = 2 (m-1) Q_(m-1) + (y - 1) Q_(m-2)`.
#[test]
fn factorial_gauge_with_index_dependent_input() {
    let input = rule(
        &[(1, 0, &[&[-2], &[2]]), (2, 0, &[&[-1, 1]])],
        Some(&[&[0], &[-1], &[1]]),
    );
    let mut candidate = candidate(2, 0, product_parent(&[1, 1]));
    candidate.gauge = VeroneseSectionGauge::Ratio {
        numerator: vec!["1".into()],
        denominator: vec!["0".into(), "1".into()],
        initial_value: "1".into(),
    };
    let certification = certify(&input, &rows(&[&[1], &[1]]), &candidate);
    assert_eq!(
        certification.verdict,
        VeroneseSectionVerdict::Proved,
        "{certification:#?}"
    );
    let VeroneseSectionIdentityStatus::Certified { certificate, .. } = &certification.identity
    else {
        panic!("not certified");
    };
    // w(n) = (n + 1)(n + 2).
    assert_eq!(certificate.multiplier, vec!["2", "3", "1"]);
    assert_replays(&certification);

    // The closed form lambda(n) = 1/n! is not rational; the same gauge with a
    // wrong initial value is refuted at the first row.
    candidate.gauge = VeroneseSectionGauge::Ratio {
        numerator: vec!["1".into()],
        denominator: vec!["0".into(), "1".into()],
        initial_value: "2".into(),
    };
    let refuted = certify(&input, &rows(&[&[1], &[1]]), &candidate);
    assert_eq!(
        refuted.identity,
        VeroneseSectionIdentityStatus::Refuted { index: 0 }
    );
    assert_eq!(refuted.verdict, VeroneseSectionVerdict::IdentityRefuted);
}

/// `2^n V_(3,1)((1+x)^(n+2))`: `Q_(n+3) = 6 Q_(n+2) - 12 Q_(n+1) + 8 (1 + y) Q_n`.
#[test]
fn shifted_and_gauged_modulus_three_section() {
    let input = rule(
        &[(1, 0, &[&[6]]), (2, 0, &[&[-12]]), (3, 0, &[&[8, 8]])],
        None,
    );
    let seeds = rows(&[&[2], &[6], &[16, 4]]);
    let mut candidate = candidate(3, 1, product_parent(&[1, 1]));
    candidate.index_shift = 2;
    candidate.gauge = VeroneseSectionGauge::Ratio {
        numerator: vec!["2".into()],
        denominator: vec!["1".into()],
        initial_value: "1".into(),
    };
    let certification = certify(&input, &seeds, &candidate);
    assert_eq!(
        certification.verdict,
        VeroneseSectionVerdict::Proved,
        "{certification:#?}"
    );
    assert_replays(&certification);

    // The other residues and shifts do not match.
    for (residue, shift) in [(0, 2), (2, 2), (1, 1)] {
        let mut other = candidate.clone();
        other.residue = residue;
        other.index_shift = shift;
        let refuted = certify(&input, &seeds, &other);
        assert!(
            matches!(
                refuted.identity,
                VeroneseSectionIdentityStatus::Refuted { .. }
            ),
            "residue {residue}, shift {shift}: {:?}",
            refuted.identity
        );
    }
}

#[test]
fn same_annihilator_with_different_seeds_is_refuted() {
    let certification = certify(
        &worked(),
        &rows(&[&[1], &[2]]),
        &candidate(2, 0, product_parent(&[1, 1])),
    );
    assert_eq!(
        certification.identity,
        VeroneseSectionIdentityStatus::Refuted { index: 1 }
    );
    assert!(!certification.verdict.is_proved());
}

#[test]
fn perturbed_input_coefficient_is_refuted() {
    // Q_(n+2) = 2 Q_(n+1) + (y - 2) Q_n: Q_2 = y, while V_(2,0)((1+x)^2) = 1 + y.
    let input = rule(&[(1, 0, &[&[2]]), (2, 0, &[&[-2, 1]])], None);
    let certification = certify(
        &input,
        &rows(&[&[1], &[1]]),
        &candidate(2, 0, product_parent(&[1, 1])),
    );
    assert_eq!(
        certification.identity,
        VeroneseSectionIdentityStatus::Refuted { index: 2 }
    );
}

#[test]
fn parents_with_mixed_signs_or_complex_roots_are_not_proved() {
    // Even sections of (1 - x)^n equal those of (1 + x)^n: the identity is
    // certified, but P_1 = 1 - x has mixed signs.
    let certification = certify(
        &worked(),
        &rows(&[&[1], &[1]]),
        &candidate(2, 0, product_parent(&[1, -1])),
    );
    assert!(matches!(
        certification.identity,
        VeroneseSectionIdentityStatus::Certified { .. }
    ));
    assert_eq!(
        certification.verdict,
        VeroneseSectionVerdict::IdentityCertifiedParentRefuted
    );
    assert_eq!(certification.parent_proof.failing_row_index, Some(1));

    // E_(n+2) = 2 (1 + y) E_(n+1) - (1 + y + y^2) E_n: the even sections of
    // (1 + x + x^2)^n, whose parent rows have non-real roots.
    let input = rule(&[(1, 0, &[&[2, 2]]), (2, 0, &[&[-1, -1, -1]])], None);
    let certification = certify(
        &input,
        &rows(&[&[1], &[1, 1]]),
        &candidate(2, 0, product_parent(&[1, 1, 1])),
    );
    assert!(matches!(
        certification.identity,
        VeroneseSectionIdentityStatus::Certified { .. }
    ));
    assert_eq!(
        certification.parent_proof.status,
        VeroneseParentProofStatus::Refuted
    );
    assert!(!certification.verdict.is_proved());
}

#[test]
fn budgets_and_limits_are_typed_errors() {
    let seeds = rows(&[&[1], &[1]]);
    let options = VeroneseSectionCertifyOptions {
        max_operator_terms: 1,
        ..VeroneseSectionCertifyOptions::default()
    };
    let error = certify_veronese_section_lift(
        &worked(),
        &seeds,
        0,
        &candidate(2, 0, product_parent(&[1, 1])),
        &options,
    )
    .unwrap_err();
    assert!(
        matches!(error, RecurrenceError::BudgetExhausted(_)),
        "{error}"
    );

    for (modulus, residue, shift, checked_rows) in [(5, 0, 0, 16), (2, 0, 5, 16), (2, 0, 0, 65)] {
        let mut candidate = candidate(modulus, residue, product_parent(&[1, 1]));
        candidate.index_shift = shift;
        let options = VeroneseSectionCertifyOptions {
            checked_rows,
            ..VeroneseSectionCertifyOptions::default()
        };
        let error =
            certify_veronese_section_lift(&worked(), &seeds, 0, &candidate, &options).unwrap_err();
        assert!(
            matches!(error, RecurrenceError::LimitExceeded(_)),
            "{error}"
        );
    }
    let error = certify_veronese_section_lift(
        &worked(),
        &seeds,
        0,
        &candidate(2, 2, product_parent(&[1, 1])),
        &VeroneseSectionCertifyOptions::default(),
    )
    .unwrap_err();
    assert!(matches!(error, RecurrenceError::InvalidInput(_)), "{error}");
    // A negative shift before the parent's first index.
    let mut early = candidate(2, 0, product_parent(&[1, 1]));
    early.index_shift = -1;
    let error = certify_veronese_section_lift(
        &worked(),
        &seeds,
        0,
        &early,
        &VeroneseSectionCertifyOptions::default(),
    )
    .unwrap_err();
    assert!(matches!(error, RecurrenceError::InvalidInput(_)), "{error}");
    // A first used parent row far beyond the parent's first index.
    let error = certify_veronese_section_lift(
        &worked(),
        &seeds,
        2000,
        &candidate(2, 0, product_parent(&[1, 1])),
        &VeroneseSectionCertifyOptions::default(),
    )
    .unwrap_err();
    assert!(
        matches!(error, RecurrenceError::LimitExceeded(_)),
        "{error}"
    );
}

#[test]
fn tampered_certificates_fail_to_replay() {
    let certification = certify(
        &worked(),
        &rows(&[&[1], &[1]]),
        &candidate(2, 0, product_parent(&[1, 1])),
    );
    let mut tampered = certification.clone();
    let VeroneseSectionIdentityStatus::Certified { certificate, .. } = &mut tampered.identity
    else {
        panic!("not certified");
    };
    certificate.multiplier = vec!["2".into()];
    assert!(tampered.verify().is_err());

    let mut tampered = certification.clone();
    tampered.verdict = VeroneseSectionVerdict::IdentityCertifiedParentOpen;
    assert!(tampered.verify().is_err());

    let mut tampered = certification;
    let VeroneseSectionIdentityStatus::Certified { certificate, .. } = &mut tampered.identity
    else {
        panic!("not certified");
    };
    certificate.residue = 1;
    assert!(tampered.verify().is_err());
}

/// `Q_m = 2 m Q_(m-1) + m (m - 1) (y - 1) Q_(m-2)`, `Q_0 = Q_1 = 1`: the even
/// sections of `P_m = m (1 + x) P_(m-1) = m! (1 + x)^m`.
fn factorial_binomial_input() -> Recurrence {
    rule(
        &[
            (1, 0, &[&[0], &[2]]),
            (2, 0, &[&[0, 0], &[1, -1], &[-1, 1]]),
        ],
        None,
    )
}

fn factorial_binomial_parent() -> RecurrenceJson {
    RecurrenceJson::from_recurrence_rational(
        &rule(&[(1, 0, &[&[0, 0], &[1, 1]])], None),
        0,
        &rows(&[&[1]]),
        None,
    )
}

/// `P_m = (x + m)^3 P_(m-1)`, `P_0 = 1`: the factor depends on the index,
/// so it is not a periodic product; it has degree three, so neither the
/// three-term nor the Liu--Wang route applies; and its index-dependent cubic
/// factor is beyond the multiplier test of the iterated-preserver route
/// (fixed factors and index-dependent factors of degree `<= 2`).  (The
/// earlier example `m (1 + x)^2` is now proved by that route.)
fn cube_product_parent() -> RecurrenceJson {
    RecurrenceJson::from_recurrence_rational(
        &rule(&[(1, 0, &[&[0, 0, 0, 1], &[0, 0, 3], &[0, 3], &[1]])], None),
        0,
        &rows(&[&[1]]),
        None,
    )
}

#[test]
fn index_dependent_parent_and_attestations() {
    let run_with = |input: &Recurrence,
                    seeds: &[Vec<BigRational>],
                    parent: &RecurrenceJson,
                    options: &VeroneseSectionCertifyOptions| {
        let certification = certify_veronese_section_lift(
            input,
            seeds,
            0,
            &candidate(2, 0, parent.clone()),
            options,
        )
        .expect("certification runs");
        certification.verify().expect("replays");
        certification
    };
    // The input is the computed rule of the even sections of the parent.
    let parent = cube_product_parent();
    let section = candidate(2, 0, parent.clone())
        .section_recurrence(0, &VeroneseSectionRecurrenceOptions::default())
        .expect("the section rule is computed");
    let (input, first_index, seeds) = section.to_recurrence_parts().unwrap();
    assert_eq!(first_index, 0);
    let run = |options: &VeroneseSectionCertifyOptions| run_with(&input, &seeds, &parent, options);
    let attestation = |record: RecurrenceJson, first_index: usize| VeroneseParentAttestation {
        label: "test attestation".into(),
        recurrence: record,
        first_index,
        statement: "every row is real-rooted with nonnegative coefficients".into(),
        source: "test".into(),
    };
    let without_root_domain = |attestations| VeroneseSectionCertifyOptions {
        attestations,
        try_root_domain: false,
        ..VeroneseSectionCertifyOptions::default()
    };

    // The parent factor (x + m)^3 depends on the index; the general product
    // rule proves it without an attestation or root-domain transformation
    // search.
    let open = run(&without_root_domain(Vec::new()));
    assert!(matches!(
        open.identity,
        VeroneseSectionIdentityStatus::Certified { .. }
    ));
    assert_eq!(open.verdict, VeroneseSectionVerdict::Proved);

    assert_eq!(
        open.parent_proof.route,
        Some(VeroneseParentProofRoute::ProductRecurrence)
    );
    // Attestations, whether matching, different or late, cannot downgrade
    // the unconditional product proof to a conditional verdict.
    let different = run(&without_root_domain(vec![attestation(
        product_parent(&[1, 1]),
        0,
    )]));
    assert_eq!(different.verdict, VeroneseSectionVerdict::Proved);
    let late = run(&without_root_domain(vec![attestation(parent.clone(), 1)]));
    assert_eq!(late.verdict, VeroneseSectionVerdict::Proved);
    let applied = run(&without_root_domain(vec![attestation(parent.clone(), 0)]));
    assert_eq!(applied.verdict, VeroneseSectionVerdict::Proved);
    assert!(applied.parent_proof.conditional_on.is_empty());
    assert_replays(&applied);

    // For the first-order parent P_m = m (1 + x) P_(m-1) of the even sections
    // of m! (1 + x)^m, preserve Liu--Wang's stronger interlacing conclusion
    // even though product induction already proves real-rootedness.
    let input = factorial_binomial_input();
    let seeds = rows(&[&[1], &[1]]);
    let first_order = factorial_binomial_parent();
    let proved = run_with(
        &input,
        &seeds,
        &first_order,
        &VeroneseSectionCertifyOptions::default(),
    );
    assert_eq!(proved.verdict, VeroneseSectionVerdict::Proved);
    assert_eq!(
        proved.parent_proof.route,
        Some(VeroneseParentProofRoute::LiuWang)
    );
    let proved = run_with(
        &input,
        &seeds,
        &first_order,
        &without_root_domain(Vec::new()),
    );
    assert_eq!(proved.verdict, VeroneseSectionVerdict::Proved);
    assert_eq!(
        proved.parent_proof.route,
        Some(VeroneseParentProofRoute::LiuWang)
    );
    assert_replays(&proved);
}

/// The even sections `E_m` of the Narayana polynomials (A001263, displayed
/// rows `N_1 = 1`, `N_2 = 1 + x`, ...) satisfy an order-4 recurrence with
/// index-dependent coefficients and a leader depending on `y`, found by an
/// exact ansatz on 40 rows (bidegree 4 in `n`, 3 in `y`; unique there).  The
/// forward certificate proves `E_m = V_(2,0)(N_m)` for every `m >= 1`, and the
/// order-two Narayana parent is proved by the three-term route (Jacobi-type
/// recurrence in `y = 1/(x - 1)`, roots in `(-1, 0)`, so every `N_m` has its
/// roots in `(-inf, 0)`), even with the root-domain route switched off.
#[test]
fn narayana_even_sections_have_a_forward_certificate() {
    let input = rule(
        &[
            (
                1,
                0,
                &[&[0, 0], &[-96, -24], &[356, 92], &[-288, -72], &[64, 16]],
            ),
            (
                2,
                0,
                &[
                    &[-144, 108, 36],
                    &[858, -180, -102],
                    &[-1158, 124, 106],
                    &[576, -48, -48],
                    &[-96, 8, 8],
                ],
            ),
            (
                3,
                0,
                &[
                    &[252, 936, 252],
                    &[-1152, -3792, -912],
                    &[1220, 3968, 924],
                    &[-480, -1560, -360],
                    &[64, 208, 48],
                ],
            ),
            (
                4,
                0,
                &[
                    &[-108, 180, -36, -36],
                    &[447, -777, 213, 117],
                    &[-425, 743, -211, -107],
                    &[144, -252, 72, 36],
                    &[-16, 28, -8, -4],
                ],
            ),
        ],
        Some(&[&[0, 0], &[57, 15], &[-7, -1], &[-48, -12], &[16, 4]]),
    );
    let (parent_rule, parent_seeds, first_row) = polytool::oeis::by_id("A001263")
        .unwrap()
        .recurrence_export_parts()
        .unwrap();
    assert_eq!(first_row, 1, "displayed rows start at N_1 = 1");
    let parent_rows = parent_rule
        .generate_rows_rational(&parent_seeds, 1, 4)
        .unwrap();
    let seeds: Vec<Vec<BigRational>> = parent_rows
        .iter()
        .map(|row| row.iter().step_by(2).cloned().collect())
        .collect();
    let parent = RecurrenceJson::from_recurrence_rational(&parent_rule, 1, &parent_seeds, None);
    let options = VeroneseSectionCertifyOptions {
        try_root_domain: false,
        ..VeroneseSectionCertifyOptions::default()
    };
    let certification = certify_veronese_section_lift(
        &input,
        &seeds,
        1,
        &candidate(2, 0, parent.clone()),
        &options,
    )
    .unwrap();
    certification.verify().unwrap();
    let VeroneseSectionIdentityStatus::Certified { certificate, .. } = &certification.identity
    else {
        panic!("not certified: {certification:#?}");
    };
    assert_eq!(certificate.identifying_indices, vec![1, 2, 3, 4]);
    // w(n) = (n + 3)(n + 4)(n + 5) has no integer roots in the domain.
    assert_eq!(certificate.multiplier, vec!["60", "47", "12", "1"]);
    assert!(certificate.exceptional_bases.is_empty());
    assert_eq!(
        certification.verdict,
        VeroneseSectionVerdict::Proved,
        "{certification:#?}"
    );
    assert_eq!(
        certification.parent_proof.route,
        Some(VeroneseParentProofRoute::ThreeTermRecurrence)
    );
    let summary = certification.lift_summary.as_ref().expect("summary");
    assert!(matches!(
        &summary.parent_proof,
        LiftParentProof::Proved { route, .. } if route == "three_term_recurrence"
    ));
    assert_eq!(summary.verdict, LiftVerdict::Proved);
    assert_replays(&certification);

    // The odd sections are a different sequence: refuted at the first row.
    let refuted =
        certify_veronese_section_lift(&input, &seeds, 1, &candidate(2, 1, parent), &options)
            .unwrap();
    assert_eq!(
        refuted.identity,
        VeroneseSectionIdentityStatus::Refuted { index: 1 }
    );
}

/// The even sections of the Eulerian polynomials (A008292, `A_1 = 1`,
/// `A_2 = 1 + x`, ...) satisfy a monic order-3 successor recurrence with
/// derivatives up to order 2 on the earlier rows, found by an exact ansatz on
/// 36 rows (unique at bidegree 3 in `n` and `y`).  The forward certificate
/// proves the identity and the root-domain route (f-polynomial map) proves
/// the parent: every even Eulerian section is real-rooted.
#[test]
fn eulerian_even_sections_are_proved() {
    let input = rule(
        &[
            (1, 0, &[&[-1], &[1]]),
            (1, 1, &[&[0, 4]]),
            (2, 0, &[&[2, 4], &[-1, -4], &[0, 1]]),
            (2, 1, &[&[0, 0, 12], &[0, -4, -4]]),
            (2, 2, &[&[0, 0, -4, 4]]),
            (3, 0, &[&[0, 12], &[0, -19], &[0, 8], &[0, -1]]),
            (3, 1, &[&[0, -6, 18], &[0, 6, -22], &[0, 0, 4]]),
            (3, 2, &[&[0, 0, -4, 4], &[0, 0, 4, -4]]),
        ],
        None,
    );
    let (parent_rule, parent_seeds, first_row) = polytool::oeis::by_id("A008292")
        .unwrap()
        .recurrence_export_parts()
        .unwrap();
    assert_eq!(first_row, 1);
    let even = |count: usize| -> Vec<Vec<BigRational>> {
        parent_rule
            .generate_rows_rational(&parent_seeds, 1, count)
            .unwrap()
            .iter()
            .map(|row| row.iter().step_by(2).cloned().collect())
            .collect()
    };
    let parent = RecurrenceJson::from_recurrence_rational(&parent_rule, 1, &parent_seeds, None);
    let certification = certify_veronese_section_lift(
        &input,
        &even(3),
        1,
        &candidate(2, 0, parent.clone()),
        &VeroneseSectionCertifyOptions::default(),
    )
    .unwrap();
    assert!(
        matches!(
            certification.identity,
            VeroneseSectionIdentityStatus::Certified { .. }
        ),
        "{certification:#?}"
    );
    assert_eq!(certification.verdict, VeroneseSectionVerdict::Proved);
    assert_eq!(
        certification.parent_proof.route,
        Some(VeroneseParentProofRoute::RootDomain)
    );
    assert_replays(&certification);

    // Extra supplied rows must agree with the recurrence; a wrong seed is refuted.
    let mut wrong = even(3);
    wrong[2] = row(&[1, 2]);
    let refuted = certify_veronese_section_lift(
        &input,
        &wrong,
        1,
        &candidate(2, 0, parent),
        &VeroneseSectionCertifyOptions::default(),
    )
    .unwrap();
    assert_eq!(
        refuted.identity,
        VeroneseSectionIdentityStatus::Refuted { index: 3 }
    );
}

/// The even sections of the type-B Eulerian polynomials (A060187) satisfy a
/// monic order-3 successor recurrence with second derivatives on the earlier
/// rows (exact ansatz on 36 rows, unique at bidegree 3).  The identity is
/// certified and the parent
/// `B_n = (1 + (2n - 3) x) B_(n-1) + 2 x (1 - x) B'_(n-1)` is proved by the
/// Liu--Wang route: every even type-B Eulerian section is real-rooted.
#[test]
fn type_b_eulerian_even_sections_are_certified() {
    let input = rule(
        &[
            (1, 0, &[&[-6], &[2]]),
            (1, 1, &[&[0, 8]]),
            (2, 0, &[&[3, 25], &[0, -20], &[0, 4]]),
            (2, 1, &[&[0, 24, 56], &[0, -16, -16]]),
            (2, 2, &[&[0, 0, -16, 16]]),
            (3, 0, &[&[4, 252], &[-2, -254], &[0, 80], &[0, -8]]),
            (3, 1, &[&[0, -64, 320], &[0, 32, -224], &[0, 0, 32]]),
            (3, 2, &[&[0, 0, -64, 64], &[0, 0, 32, -32]]),
        ],
        None,
    );
    let (parent_rule, parent_seeds, first_row) = polytool::oeis::by_id("A060187")
        .unwrap()
        .recurrence_export_parts()
        .unwrap();
    assert_eq!(first_row, 1);
    let seeds: Vec<Vec<BigRational>> = parent_rule
        .generate_rows_rational(&parent_seeds, 1, 3)
        .unwrap()
        .iter()
        .map(|row| row.iter().step_by(2).cloned().collect())
        .collect();
    let parent = RecurrenceJson::from_recurrence_rational(&parent_rule, 1, &parent_seeds, None);
    let certification = certify_veronese_section_lift(
        &input,
        &seeds,
        1,
        &candidate(2, 0, parent),
        &VeroneseSectionCertifyOptions::default(),
    )
    .unwrap();
    assert!(
        matches!(
            certification.identity,
            VeroneseSectionIdentityStatus::Certified { .. }
        ),
        "{certification:#?}"
    );
    // The root-domain menu does not prove the type-B Eulerian parent; the
    // Liu--Wang route does.
    assert_eq!(certification.verdict, VeroneseSectionVerdict::Proved);
    assert_eq!(
        certification.parent_proof.route,
        Some(VeroneseParentProofRoute::LiuWang)
    );
    let root_domain_attempts = certification
        .parent_proof
        .attempts
        .iter()
        .filter(|a| a.route == VeroneseParentProofRoute::RootDomain)
        .count();
    assert!(root_domain_attempts > 0);
    assert!(certification
        .parent_proof
        .attempts
        .iter()
        .filter(|a| a.route == VeroneseParentProofRoute::RootDomain)
        .all(|a| !a.success));
    assert_eq!(
        certification.lift_summary.as_ref().unwrap().verdict,
        LiftVerdict::Proved
    );
    assert_replays(&certification);
}

/// The section rows `lambda(n) V_(r, j)(P_(n+s))`, `n = n0, ..., n0 + count - 1`,
/// with `lambda(n) = 2^(n - n0)` when `doubling`.
fn section_rows(
    parent: &RecurrenceJson,
    modulus: usize,
    residue: usize,
    first_index: usize,
    shift: usize,
    doubling: bool,
    count: usize,
) -> Vec<Vec<BigRational>> {
    let (rule, p0, seeds) = parent.to_recurrence_parts().unwrap();
    let start = first_index + shift - p0;
    let rows = rule
        .generate_rows_rational(&seeds, p0, start + count)
        .unwrap();
    let mut factor = q(1);
    rows[start..]
        .iter()
        .map(|row| {
            let mut section: Vec<BigRational> = row
                .iter()
                .skip(residue)
                .step_by(modulus)
                .map(|c| c * &factor)
                .collect();
            while section.last().is_some_and(|c| *c == q(0)) {
                section.pop();
            }
            if doubling {
                factor *= q(2);
            }
            section
        })
        .collect()
}

fn trimmed(rows: Vec<Vec<BigRational>>) -> Vec<Vec<BigRational>> {
    rows.into_iter()
        .map(|mut row| {
            while row.last().is_some_and(|c| *c == q(0)) {
                row.pop();
            }
            row
        })
        .collect()
}

/// `T_n = x T_(n-1) + x T'_(n-1)`, `T_0 = 1` (Touchard rows).
fn touchard() -> RecurrenceJson {
    RecurrenceJson::from_recurrence_rational(
        &rule(&[(1, 0, &[&[0, 1]]), (1, 1, &[&[0, 1]])], None),
        0,
        &rows(&[&[1]]),
        None,
    )
}

/// Section rules computed from the parent rule (not fitted) generate the
/// sections and are certified for every index: a derivative-free parent of
/// order two with an index-dependent leader (Narayana) and first-order
/// parents with derivatives (Touchard, Eulerian).  The modulus-3 sections
/// of the Eulerian rows need a rule of order 5 with derivatives up to order
/// 6, beyond the default budget, so a larger one is passed there.
#[test]
fn computed_section_rules_generate_the_sections_and_are_certified() {
    let large = VeroneseSectionRecurrenceOptions {
        max_order_increase: Some(4),
        max_derivative_order: Some(6),
        max_unknowns: 2000,
        max_shapes: 500,
        ..VeroneseSectionRecurrenceOptions::default()
    };
    let default = VeroneseSectionRecurrenceOptions::default();
    let eulerian = polytool::recurrence::GammaLiftParent::eulerian().recurrence;
    let cases = [
        (
            "narayana",
            polytool::recurrence::GammaLiftParent::narayana().recurrence,
            2,
            &default,
        ),
        (
            "narayana",
            polytool::recurrence::GammaLiftParent::narayana().recurrence,
            3,
            &default,
        ),
        ("touchard", touchard(), 2, &default),
        ("touchard", touchard(), 3, &default),
        ("eulerian", eulerian.clone(), 2, &default),
        ("eulerian", eulerian, 3, &large),
    ];
    for (name, parent, modulus, options) in cases {
        let residues = if modulus == 3 && name == "eulerian" {
            0..1
        } else {
            0..modulus
        };
        for residue in residues {
            let candidate = VeroneseSectionCandidate {
                modulus,
                residue,
                parent: parent.clone(),
                index_shift: 0,
                gauge: VeroneseSectionGauge::None,
            };
            let p0 = parent.first_index;
            let record = candidate
                .section_recurrence(p0, options)
                .unwrap_or_else(|error| panic!("{name} r={modulus} j={residue}: {error}"));
            let (input, first_index, seeds) = record.to_recurrence_parts().unwrap();
            assert_eq!(first_index, p0);
            let generated = input.generate_rows_rational(&seeds, p0, 30).unwrap();
            assert_eq!(
                trimmed(generated),
                section_rows(&parent, modulus, residue, p0, 0, false, 30),
                "{name} r={modulus} j={residue}"
            );
            let certification = certify_veronese_section_lift(
                &input,
                &seeds,
                first_index,
                &candidate,
                &VeroneseSectionCertifyOptions::default(),
            )
            .unwrap();
            certification.verify().unwrap();
            assert!(
                matches!(
                    certification.identity,
                    VeroneseSectionIdentityStatus::Certified { .. }
                ),
                "{name} r={modulus} j={residue}: {certification:#?}"
            );
        }
    }
}

/// A shifted and gauged candidate `Q_n = 2^n V_(2, 1)(N_(n+1))`.
#[test]
fn computed_section_rule_of_a_shifted_gauged_candidate_is_certified() {
    let parent = polytool::recurrence::GammaLiftParent::narayana().recurrence;
    let gauge = VeroneseSectionGauge::Ratio {
        numerator: vec!["2".into()],
        denominator: vec!["1".into()],
        initial_value: "1".into(),
    };
    let candidate = VeroneseSectionCandidate {
        modulus: 2,
        residue: 1,
        parent: parent.clone(),
        index_shift: 1,
        gauge,
    };
    let record = candidate
        .section_recurrence(1, &VeroneseSectionRecurrenceOptions::default())
        .unwrap();
    let (input, first_index, seeds) = record.to_recurrence_parts().unwrap();
    let generated = input
        .generate_rows_rational(&seeds, first_index, 25)
        .unwrap();
    assert_eq!(
        trimmed(generated),
        section_rows(&parent, 2, 1, 1, 1, true, 25)
    );
    let certification = certify_veronese_section_lift(
        &input,
        &seeds,
        first_index,
        &candidate,
        &VeroneseSectionCertifyOptions::default(),
    )
    .unwrap();
    assert!(
        matches!(
            certification.identity,
            VeroneseSectionIdentityStatus::Certified { .. }
        ),
        "{certification:#?}"
    );
}

/// With a budget too small for the elimination the search says so.
#[test]
fn computed_section_rule_reports_an_exhausted_budget() {
    let candidate = VeroneseSectionCandidate {
        modulus: 3,
        residue: 0,
        parent: polytool::recurrence::GammaLiftParent::narayana().recurrence,
        index_shift: 0,
        gauge: VeroneseSectionGauge::None,
    };
    let options = VeroneseSectionRecurrenceOptions {
        max_order_increase: Some(0),
        max_coefficient_degree: 1,
        ..VeroneseSectionRecurrenceOptions::default()
    };
    assert!(matches!(
        candidate.section_recurrence(1, &options),
        Err(RecurrenceError::BudgetExhausted(_))
    ));
}
