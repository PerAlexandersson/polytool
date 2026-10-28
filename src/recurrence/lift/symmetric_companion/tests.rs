//! Unit tests of the symmetric-companion lifts: the worked `1/k`-Eulerian
//! fixture, the operator identities, and near misses.

use super::certificate::build_certificate;
use super::operators::{h_operator, is_symmetric_at, CompanionParameters};
use super::search::infer_schedule;
use super::source::PreparedSource;
use super::*;
use crate::recurrence::implicit::algebra::Biv;
use crate::recurrence::*;
use crate::sequences::literature::one_over_k_eulerian_symmetric_parts_rational;
use crate::vec_poly;
use num_bigint::BigInt;
use num_traits::Zero;

type Q = BigRational;

fn q(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

fn qr(numerator: i64, denominator: i64) -> Q {
    Q::new(BigInt::from(numerator), BigInt::from(denominator))
}

fn row(values: &[Q]) -> Vec<Q> {
    vec_poly::trim(values)
}

/// The expanded eliminated recurrence of the `1/k`-Eulerian symmetric parts,
/// in the target index `m = n + 2`, written out by hand from the issue:
///
/// ```text
/// A_m = (U_(m-1) + Vt_(m-2)) A_(m-1) - (Vt_(m-2) o U_(m-2) - (k - 1) x) A_(m-2),
/// U_n = [1 + (1 + k (n - 1)) x] + k x (1 - x) D,   Vt_n = k n x + k x (1 - x) D.
/// ```
fn eliminated_one_over_k_recurrence(k: &Q) -> Recurrence {
    let m = Biv::new(vec![vec![], vec![q(1)]]);
    let x = Biv::from_y(&[q(0), q(1)]);
    let one_minus_x = Biv::from_y(&[q(1), q(-1)]);
    let x_one_minus_x = x.mul(&one_minus_x);
    let kk = Biv::constant(k.clone());
    let constant = |c: i64| Biv::constant(q(c));
    let n_big = m.add(&constant(-2));
    let alpha = Biv::one();
    let beta = Biv::one().add(&kk.mul(&m.add(&constant(-3))));
    // (U_(m-1) + Vt_(m-2)): 1 + (1 + 2k(m - 2)) x and 2k x (1 - x) D.
    let first_multiplication =
        Biv::one().add(&x.mul(&Biv::one().add(&kk.scale(&q(2)).mul(&n_big))));
    let first_derivative = kk.scale(&q(2)).mul(&x_one_minus_x);
    // Vt o U - (k - 1) x.
    let alpha_beta_x = alpha.add(&beta.mul(&x));
    let second_multiplication = kk
        .mul(&n_big)
        .mul(&x)
        .mul(&alpha_beta_x)
        .add(&kk.mul(&beta).mul(&x_one_minus_x))
        .add(&x.mul(&kk.add(&constant(-1))).scale(&q(-1)));
    let one_minus_two_x = Biv::from_y(&[q(1), q(-2)]);
    let second_derivative = kk
        .mul(&kk)
        .mul(&n_big)
        .mul(&x)
        .mul(&x_one_minus_x)
        .add(&kk.mul(&x_one_minus_x).mul(&alpha_beta_x))
        .add(&kk.mul(&kk).mul(&x_one_minus_x).mul(&one_minus_two_x));
    let second_second = kk.mul(&kk).mul(&x_one_minus_x).mul(&x_one_minus_x);
    let minus = |b: &Biv| b.scale(&q(-1)).to_bivar();
    Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, first_multiplication.to_bivar()),
            RecurrenceTerm::new(1, 1, first_derivative.to_bivar()),
            RecurrenceTerm::new(2, 0, minus(&second_multiplication)),
            RecurrenceTerm::new(2, 1, minus(&second_derivative)),
            RecurrenceTerm::new(2, 2, minus(&second_second)),
        ],
        denominator: None,
        inhomogeneous: None,
    }
}

fn one_over_k_seeds(k: &Q) -> Vec<Vec<Q>> {
    vec![row(&[q(1), q(1)]), row(&[q(1), k * q(3) + q(1), q(1)])]
}

fn ks() -> Vec<Q> {
    vec![q(1), q(2), qr(3, 2)]
}

fn find(
    recurrence: &Recurrence,
    seeds: &[Vec<Q>],
    first_index: usize,
) -> SymmetricCompanionLiftReport {
    find_symmetric_companion_lifts(
        recurrence,
        seeds,
        first_index,
        &SymmetricCompanionLiftOptions::default(),
    )
    .unwrap()
}

/// Altered inputs may be accepted only through a separately valid certificate:
/// every certified lift replays, and its parent's symmetric parts are the
/// input rows.
fn assert_certified_lifts_are_genuine(
    report: &SymmetricCompanionLiftReport,
    recurrence: &Recurrence,
    seeds: &[Vec<Q>],
    first_index: usize,
) {
    let source = PreparedSource::new(recurrence, seeds, first_index).unwrap();
    let rows = source.rows(12).unwrap();
    for lift in report.certified_lifts() {
        let replay = verify_symmetric_companion_certificate(&lift.certificate, 8).unwrap();
        assert_eq!(replay, lift.verification);
        let (parent, parent_first, parent_seeds) =
            lift.certificate.parent.to_recurrence_parts().unwrap();
        let parents = parent
            .generate_rows_rational(&parent_seeds, parent_first, rows.len())
            .unwrap();
        let d0 = lift.certificate.symmetry_degree_at_first_index;
        for (offset, (a, p)) in rows.iter().zip(&parents).enumerate() {
            let (symmetric, _) = crate::polynomial::Polynomial::new(p.clone())
                .stapledon_decomposition(d0 + offset)
                .unwrap();
            assert_eq!(vec_poly::trim(symmetric.coeffs()), *a);
        }
    }
}

#[test]
fn hand_recurrence_reproduces_the_canonical_symmetric_parts() {
    for k in ks() {
        let source = PreparedSource::new(
            &eliminated_one_over_k_recurrence(&k),
            &one_over_k_seeds(&k),
            2,
        )
        .unwrap();
        let rows = source.rows(12).unwrap();
        let parts = one_over_k_eulerian_symmetric_parts_rational(13, &k);
        for (offset, generated) in rows.iter().enumerate() {
            assert_eq!(
                *generated,
                parts[2 + offset].0,
                "k = {k}, n = {}",
                2 + offset
            );
        }
        assert_eq!(infer_schedule(&rows), Ok(1));
    }
}

#[test]
fn worked_fixture_recovers_the_parent_and_the_theorem_applies() {
    for k in ks() {
        let recurrence = eliminated_one_over_k_recurrence(&k);
        let report = find(&recurrence, &one_over_k_seeds(&k), 2);
        let lift = report
            .proving_lift()
            .unwrap_or_else(|| panic!("k = {k}: {report:#?}"));
        assert_eq!(
            report.symmetry_schedule.symmetry_degree_at_first_index,
            Some(1)
        );
        let (a, b, c, delta, epsilon) = lift.certificate.parameters.values().unwrap();
        assert_eq!(a, vec![q(1)]);
        assert_eq!(b, vec![q(0), k.clone()]);
        assert_eq!(c, vec![k.clone()]);
        assert_eq!(delta, row(&[&k - q(1)]));
        assert_eq!(epsilon, vec![q(1)]);
        let expected_companion = format_rational_coeff(&(&k - q(1)));
        assert_eq!(lift.certificate.companion_seed, vec![expected_companion]);
        let (parent, first, seeds) = lift.certificate.parent.to_recurrence_parts().unwrap();
        assert_eq!(first, 2);
        assert_eq!(seeds, vec![row(&[q(1), k.clone()])]);
        // The parent rows are the 1/k-Eulerian polynomials.
        let parents = parent.generate_rows_rational(&seeds, 2, 10).unwrap();
        let canonical =
            crate::sequences::literature::one_over_k_eulerian_polynomials_rational(11, &k);
        for (offset, p) in parents.iter().enumerate() {
            assert_eq!(vec_poly::trim(p), canonical[2 + offset]);
        }
        assert!(lift.verification.is_certified());
        assert!(lift
            .verification
            .identities
            .iter()
            .all(|c| c.status == SymmetricCompanionCheckStatus::Verified));
        assert_eq!(lift.verification.lean_replay, "not_attempted");
        assert!(lift.verification.unresolved_obligations.is_empty());
        // Independent replay, also after a JSON round trip.
        let replay = verify_symmetric_companion_certificate(&lift.certificate, 10).unwrap();
        assert!(replay.theorem_applies());
        let json = serde_json::to_string(&lift.certificate).unwrap();
        let decoded: SymmetricCompanionCertificate = serde_json::from_str(&json).unwrap();
        assert_eq!(decoded, lift.certificate);
        assert!(verify_symmetric_companion_certificate(&decoded, 6)
            .unwrap()
            .theorem_applies());
        let report_json = report.to_json_string().unwrap();
        assert!(report_json.contains(SYMMETRIC_COMPANION_LIFT_SCHEMA));
        let decoded: SymmetricCompanionLiftReport = serde_json::from_str(&report_json).unwrap();
        assert_eq!(decoded, report);
    }
}

#[test]
fn zero_companion_for_k_equal_to_one() {
    let k = q(1);
    let report = find(
        &eliminated_one_over_k_recurrence(&k),
        &one_over_k_seeds(&k),
        2,
    );
    let lift = report.proving_lift().unwrap();
    assert_eq!(lift.certificate.companion_seed, vec!["0".to_string()]);
    let (_, _, _, delta, _) = lift.certificate.parameters.values().unwrap();
    assert_eq!(delta, vec![q(0)]);
}

#[test]
fn operator_identities_hold_for_arbitrary_parameters() {
    let parameters = CompanionParameters::new(
        vec![q(2), qr(1, 3)],
        vec![qr(5, 2), q(-1)],
        vec![q(3), q(7)],
        4,
        6,
    );
    let ops = parameters.operators();
    assert!(ops.elimination_identity_holds());
    assert!(ops.collapse_identity_holds());
    let projection = ops.constant_term_projection().unwrap();
    assert!(ops.constant_term_factorization_holds(&projection));
    // H(d, a, c) maps symmetric-at-d to symmetric-at-(d+1), deficits included.
    let h = h_operator(&[q(3)], &[q(2)], &[qr(1, 2)]);
    for f in [
        row(&[q(1), q(4), q(4), q(1)]),
        row(&[q(0), q(2), q(2)]),
        row(&[q(0)]),
    ] {
        assert!(is_symmetric_at(&f, 3));
        let image = super::operators::apply_to_row(&h, &f, 0);
        assert!(is_symmetric_at(&image, 4), "{image:?}");
    }
}

#[test]
fn nonconstant_leader_uses_the_proportional_derivation() {
    // Multiply the fixture rule by (m + 1): the leader is no longer constant.
    let k = q(2);
    let mut recurrence = eliminated_one_over_k_recurrence(&k);
    let factor = Biv::new(vec![vec![q(1)], vec![q(1)]]);
    for term in &mut recurrence.terms {
        term.coeff = Biv::from_bivar(&term.coeff).mul(&factor).to_bivar();
    }
    recurrence.denominator = Some(factor.to_bivar());
    let report = find(&recurrence, &one_over_k_seeds(&k), 2);
    let lift = report.proving_lift().unwrap();
    assert_eq!(
        lift.certificate.source_consequence_factor,
        vec!["3".to_string(), "1".to_string()]
    );
    assert!(lift.certificate.exceptional_bases.is_empty());
}

#[test]
fn higher_order_source_without_derivation_stays_a_candidate() {
    // M = (S + 1) o E has order three; E is not a left multiple of M.
    let k = q(2);
    let parameters = CompanionParameters::new(vec![q(1)], vec![k.clone()], vec![q(1)], 2, 1);
    let e = parameters.operators().eliminated;
    let s_plus_one = super::operators::shift_operator().add(&super::operators::scalar(&[q(1)]));
    let m = ImplicitRelation::from_op(&s_plus_one.compose(&e));
    let (recurrence, order) = m.to_successor_rule().unwrap();
    assert_eq!(order, 3);
    let fixture = PreparedSource::new(
        &eliminated_one_over_k_recurrence(&k),
        &one_over_k_seeds(&k),
        2,
    )
    .unwrap();
    let seeds = fixture.rows(3).unwrap();
    let report = find(&recurrence, &seeds, 2);
    assert!(report.certified_lifts().next().is_none());
    let lift = &report.lifts[0];
    assert_eq!(
        lift.verification.identity_status,
        SymmetricCompanionIdentityStatus::Candidate
    );
    assert_eq!(
        lift.verification.real_rootedness_status,
        SymmetricCompanionRealRootednessStatus::NotEvaluated
    );
    assert!(lift
        .verification
        .unresolved_obligations
        .iter()
        .any(|o| o.starts_with("eliminated_annihilates_source")));
}

#[test]
fn changed_seeds_and_coefficients_are_not_accepted_without_a_valid_certificate() {
    let k = q(2);
    let recurrence = eliminated_one_over_k_recurrence(&k);
    // Changed seed: still symmetric, different interior coefficient.
    let seeds = vec![row(&[q(1), q(1)]), row(&[q(1), q(8), q(1)])];
    let report = find(&recurrence, &seeds, 2);
    assert_certified_lifts_are_genuine(&report, &recurrence, &seeds, 2);
    // Here a separately valid certificate exists: the same rule with the
    // companion seed B_2 = 2, i.e. the parent seed P_2 = 1 + 3x.
    let lift = report.proving_lift().unwrap();
    assert_eq!(lift.certificate.companion_seed, vec!["2".to_string()]);
    // Changed recurrence coefficient.
    let mut changed = recurrence.clone();
    changed.terms[0].coeff.coeffs[0][1] += q(1);
    let seeds = one_over_k_seeds(&k);
    let report = find(&changed, &seeds, 2);
    assert_certified_lifts_are_genuine(&report, &changed, &seeds, 2);
    assert!(report.lifts.is_empty());
    assert!(!report.obstructions.is_empty());
}

#[test]
fn wrong_schedule_is_rejected() {
    let k = q(2);
    let options = SymmetricCompanionLiftOptions {
        symmetry_degree_at_first_index: Some(2),
        ..SymmetricCompanionLiftOptions::default()
    };
    let report = find_symmetric_companion_lifts(
        &eliminated_one_over_k_recurrence(&k),
        &one_over_k_seeds(&k),
        2,
        &options,
    )
    .unwrap();
    assert!(report.certified_lifts().next().is_none());
    for lift in &report.lifts {
        assert_eq!(
            lift.verification.identity_status,
            SymmetricCompanionIdentityStatus::Rejected
        );
    }
    // A certificate with the true parameters but the wrong schedule is
    // rejected by the seed symmetry check.
    let k_parameters = CompanionParameters::new(vec![q(1)], vec![k.clone()], vec![q(1)], 2, 2);
    let source = PreparedSource::new(
        &eliminated_one_over_k_recurrence(&k),
        &one_over_k_seeds(&k),
        2,
    )
    .unwrap();
    if let Ok(certificate) = build_certificate(&source, &k_parameters) {
        let verification = verify_symmetric_companion_certificate(&certificate, 6).unwrap();
        assert_ne!(
            verification.identity_status,
            SymmetricCompanionIdentityStatus::Certified
        );
    }
}

#[test]
fn constant_terms_that_break_divisibility_are_an_obstruction() {
    let k = q(2);
    let recurrence = eliminated_one_over_k_recurrence(&k);
    let seeds = vec![row(&[q(1), q(1)]), row(&[q(2), q(7), q(2)])];
    let report = find(&recurrence, &seeds, 2);
    assert!(report.lifts.is_empty());
    assert!(report
        .search
        .attempts
        .iter()
        .all(|a| a.outcome == SymmetricCompanionAttemptOutcome::ConstantTermsIncompatible));
    // The verifier refutes R_(n0)(0) = 0 on a certificate whose seeds were altered.
    let good = find(&recurrence, &one_over_k_seeds(&k), 2);
    let mut certificate = good.proving_lift().unwrap().certificate.clone();
    certificate.source.initial_polynomials[1] = vec!["2".into(), "7".into(), "2".into()];
    let verification = verify_symmetric_companion_certificate(&certificate, 6).unwrap();
    assert_eq!(
        verification.identity_status,
        SymmetricCompanionIdentityStatus::Rejected
    );
    let divisibility = verification
        .identities
        .iter()
        .find(|c| c.name == "constant_term_divisibility")
        .unwrap();
    assert_eq!(divisibility.status, SymmetricCompanionCheckStatus::Failed);
}

/// Rows `A_(n0..)` of the coupled system with the given parameters and
/// seeds, and the eliminated operator as a successor rule.
fn system_fixture(
    parameters: &CompanionParameters,
    a_seed: Vec<Q>,
    b_seed: Vec<Q>,
    count: usize,
) -> (Recurrence, Vec<Vec<Q>>) {
    let ops = parameters.operators();
    let mut a = vec![a_seed];
    let mut b = vec![b_seed];
    for offset in 0..count {
        let n = parameters.first_index + offset;
        let epsilon = super::operators::evaluate_at_index(&parameters.epsilon, n);
        let delta = super::operators::evaluate_at_index(&parameters.delta(), n);
        let next_a = vec_poly::add(
            &super::operators::apply_to_row(&ops.u, &a[offset], n),
            &vec_poly::scale(&vec_poly::shift(&b[offset]), &epsilon),
        );
        let next_b = vec_poly::add(
            &vec_poly::scale(&a[offset], &delta),
            &super::operators::apply_to_row(&ops.v, &b[offset], n),
        );
        a.push(vec_poly::trim(&next_a));
        b.push(vec_poly::trim(&next_b));
    }
    let (recurrence, _) = ImplicitRelation::from_op(&ops.eliminated)
        .to_successor_rule()
        .unwrap();
    (recurrence, a)
}

#[test]
fn vanishing_coupling_is_unsupported() {
    // a = 1, c = 1, epsilon_n = n - 4 (vanishes at n = 4), delta_n = 5 - n.
    let parameters = CompanionParameters::new(vec![q(1)], vec![q(1)], vec![q(-4), q(1)], 2, 1);
    let (recurrence, rows) = system_fixture(&parameters, row(&[q(1), q(1)]), vec![q(3)], 8);
    // The leader epsilon vanishes at base 4 (target 6): rows through 6 are seeds.
    let source = PreparedSource::new(&recurrence, &rows, 2).unwrap();
    assert!(source.seeds.len() >= 5);
    let certificate = build_certificate(&source, &parameters).unwrap();
    let verification = verify_symmetric_companion_certificate(&certificate, 6).unwrap();
    assert_eq!(
        verification.identity_status,
        SymmetricCompanionIdentityStatus::Unsupported
    );
    let coupling = verification
        .identities
        .iter()
        .find(|c| c.name == "coupling_nonvanishing")
        .unwrap();
    assert_eq!(coupling.status, SymmetricCompanionCheckStatus::Unsupported);
    // The search finds the same candidate and reports it as unsupported.
    let report = find(&recurrence, &rows, 2);
    assert!(report.certified_lifts().next().is_none());
    assert_eq!(
        report.lifts[0].verification.identity_status,
        SymmetricCompanionIdentityStatus::Unsupported
    );
    assert_certified_lifts_are_genuine(&report, &recurrence, &rows, 2);
}

#[test]
fn degree_deficit_with_zero_constant_terms() {
    // x A_n with d_n + 2: a' = a - c = 1 - k, the same c and epsilon.
    for (k, theorem) in [(q(1), true), (q(2), false)] {
        let parameters =
            CompanionParameters::new(vec![q(1) - &k], vec![k.clone()], vec![q(1)], 2, 3);
        let (recurrence, rows) = system_fixture(
            &parameters,
            row(&[q(0), q(1), q(1)]),
            row(&[q(0), &k - q(1)]),
            8,
        );
        assert!(rows.iter().all(|r| r[0].is_zero()));
        let report = find(&recurrence, &rows[..2], 2);
        assert_eq!(
            report.symmetry_schedule.symmetry_degree_at_first_index,
            Some(3)
        );
        let lift = report
            .certified_lifts()
            .next()
            .unwrap_or_else(|| panic!("{report:#?}"));
        assert!(lift.verification.self_checked_rows > 0);
        assert_eq!(lift.verification.theorem_applies(), theorem, "k = {k}");
        if !theorem {
            assert_eq!(
                lift.verification.real_rootedness_status,
                SymmetricCompanionRealRootednessStatus::HypothesesNotVerified
            );
            let a = lift
                .verification
                .hypotheses
                .iter()
                .find(|c| c.name == "a_nonnegative")
                .unwrap();
            assert_eq!(a.status, SymmetricCompanionCheckStatus::Failed);
        }
        assert_certified_lifts_are_genuine(&report, &recurrence, &rows[..2], 2);
    }
}

#[test]
fn exhausted_budget_is_reported() {
    let k = q(2);
    let options = SymmetricCompanionLiftOptions {
        max_linearized_monomials: 3,
        ..SymmetricCompanionLiftOptions::default()
    };
    let report = find_symmetric_companion_lifts(
        &eliminated_one_over_k_recurrence(&k),
        &one_over_k_seeds(&k),
        2,
        &options,
    )
    .unwrap();
    assert!(report.lifts.is_empty());
    assert!(report.search.exhausted_budget);
    assert!(report
        .search
        .attempts
        .iter()
        .all(|a| a.outcome == SymmetricCompanionAttemptOutcome::BudgetExceeded));
}

#[test]
fn tampered_certificates_are_refused() {
    let k = qr(3, 2);
    let report = find(
        &eliminated_one_over_k_recurrence(&k),
        &one_over_k_seeds(&k),
        2,
    );
    let certificate = report.proving_lift().unwrap().certificate.clone();
    let mut changed = certificate.clone();
    changed.parameters.c = vec!["2".into()];
    assert!(verify_symmetric_companion_certificate(&changed, 6).is_err());
    let mut changed = certificate.clone();
    changed.companion_seed = vec!["1".into()];
    assert!(verify_symmetric_companion_certificate(&changed, 6).is_err());
    let mut changed = certificate.clone();
    changed.parent.initial_polynomials = vec![vec!["1".into(), "2".into()]];
    assert!(verify_symmetric_companion_certificate(&changed, 6).is_err());
    let mut changed = certificate;
    changed.exceptional_bases = vec![2];
    assert!(verify_symmetric_companion_certificate(&changed, 6).is_err());
}

#[test]
fn unsupported_inputs_are_errors() {
    let k = q(2);
    let mut recurrence = eliminated_one_over_k_recurrence(&k);
    recurrence.denominator = Some(BivarPoly::from_coefficients_normalized(vec![vec![
        q(1),
        q(1),
    ]]));
    assert!(matches!(
        find_symmetric_companion_lifts(
            &recurrence,
            &one_over_k_seeds(&k),
            2,
            &SymmetricCompanionLiftOptions::default()
        ),
        Err(RecurrenceError::Unsupported(_))
    ));
    let options = SymmetricCompanionLiftOptions {
        max_parameter_degree_in_n: 2,
        ..SymmetricCompanionLiftOptions::default()
    };
    assert!(find_symmetric_companion_lifts(
        &eliminated_one_over_k_recurrence(&k),
        &one_over_k_seeds(&k),
        2,
        &options
    )
    .is_err());
}
