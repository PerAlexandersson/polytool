//! End-to-end tests of `recurrence::find_symmetric_companion_lifts` and
//! `recurrence::verify_symmetric_companion_certificate` on the worked
//! `1/k`-Eulerian fixture (rust-workspace issue 9), through the public API.

use num_bigint::BigInt;
use polytool::recurrence::{
    find_symmetric_companion_lifts, verify_symmetric_companion_certificate, BigRational, BivarPoly,
    Recurrence, RecurrenceTerm, SymmetricCompanionAttemptOutcome, SymmetricCompanionCertificate,
    SymmetricCompanionCheckStatus, SymmetricCompanionIdentityStatus, SymmetricCompanionLiftOptions,
    SymmetricCompanionLiftReport, SymmetricCompanionRealRootednessStatus,
    SYMMETRIC_COMPANION_LIFT_SCHEMA,
};
use polytool::sequences::literature::{
    one_over_k_eulerian_polynomials_rational, one_over_k_eulerian_symmetric_parts_rational,
};

type Q = BigRational;

fn q(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

fn qr(numerator: i64, denominator: i64) -> Q {
    Q::new(BigInt::from(numerator), BigInt::from(denominator))
}

fn bivariate(rows: Vec<Vec<Q>>) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(rows)
}

/// The expanded eliminated recurrence of the `1/k`-Eulerian symmetric parts
/// (target index `m`, `coeffs[i][j]` multiplies `m^i x^j`):
///
/// ```text
/// A_m = [1 + (1 + 2k(m - 2)) x] A_(m-1) + 2k x (1 - x) A'_(m-1)
///       - [(1 + k(m - 2) + k^2 (m - 3)) x + k(m - 3)(1 + k(m - 3)) x^2] A_(m-2)
///       - x (1 - x) [(k + k^2) + (k + k^2 (2m - 7)) x] A'_(m-2)
///       - k^2 x^2 (1 - x)^2 A''_(m-2).
/// ```
fn eliminated_recurrence(k: &Q) -> Recurrence {
    let k2 = k * k;
    let zero = q(0);
    let t10 = vec![vec![q(1), q(1) - k * q(4)], vec![zero.clone(), k * q(2)]];
    let t11 = vec![vec![zero.clone(), k * q(2), -(k * q(2))]];
    let t20 = vec![
        vec![
            zero.clone(),
            -(q(1) - k * q(2) - &k2 * q(3)),
            -(-(k * q(3)) + &k2 * q(9)),
        ],
        vec![zero.clone(), -(k + &k2), -(k - &k2 * q(6))],
        vec![zero.clone(), zero.clone(), -k2.clone()],
    ];
    let t21 = vec![
        vec![zero.clone(), -(k + &k2), &k2 * q(8), k - &k2 * q(7)],
        vec![zero.clone(), zero.clone(), -(&k2 * q(2)), &k2 * q(2)],
    ];
    let t22 = vec![vec![
        zero.clone(),
        zero,
        -k2.clone(),
        &k2 * q(2),
        -k2.clone(),
    ]];
    Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, bivariate(t10)),
            RecurrenceTerm::new(1, 1, bivariate(t11)),
            RecurrenceTerm::new(2, 0, bivariate(t20)),
            RecurrenceTerm::new(2, 1, bivariate(t21)),
            RecurrenceTerm::new(2, 2, bivariate(t22)),
        ],
        denominator: None,
        inhomogeneous: None,
    }
}

/// `A_2 = 1 + x`, `A_3 = 1 + (3k + 1) x + x^2`.
fn seeds(k: &Q) -> Vec<Vec<Q>> {
    vec![vec![q(1), q(1)], vec![q(1), k * q(3) + q(1), q(1)]]
}

fn search(
    k: &Q,
    seeds: &[Vec<Q>],
    options: &SymmetricCompanionLiftOptions,
) -> SymmetricCompanionLiftReport {
    find_symmetric_companion_lifts(&eliminated_recurrence(k), seeds, 2, options).unwrap()
}

#[test]
fn the_eliminated_recurrence_generates_the_canonical_symmetric_parts() {
    for k in [q(1), q(2), qr(3, 2)] {
        let rows = eliminated_recurrence(&k)
            .generate_rows_rational(&seeds(&k), 2, 14)
            .unwrap();
        let parts = one_over_k_eulerian_symmetric_parts_rational(15, &k);
        for (offset, row) in rows.iter().enumerate() {
            assert_eq!(
                polytool::vec_poly::trim(row),
                parts[2 + offset].0,
                "k = {k}, n = {}",
                2 + offset
            );
        }
    }
}

#[test]
fn worked_fixture_is_certified_and_the_preservation_theorem_applies() {
    for k in [q(1), q(2), qr(3, 2)] {
        let report = search(&k, &seeds(&k), &SymmetricCompanionLiftOptions::default());
        assert_eq!(report.schema, SYMMETRIC_COMPANION_LIFT_SCHEMA);
        let lift = report.proving_lift().expect("a proving lift");
        let verification = &lift.verification;
        assert_eq!(
            verification.identity_status,
            SymmetricCompanionIdentityStatus::Certified
        );
        assert_eq!(
            verification.real_rootedness_status,
            SymmetricCompanionRealRootednessStatus::TheoremApplies
        );
        assert_eq!(verification.lean_replay, "not_attempted");
        assert!(verification
            .hypotheses
            .iter()
            .all(|h| h.status == SymmetricCompanionCheckStatus::Verified));
        // a_n = 1, b_n = k n, c_n = k, delta_n = k - 1, epsilon_n = 1.
        let (a, b, c, delta, epsilon) = lift.certificate.parameters.values().unwrap();
        assert_eq!(a, vec![q(1)]);
        assert_eq!(b, vec![q(0), k.clone()]);
        assert_eq!(c, vec![k.clone()]);
        assert_eq!(delta, polytool::vec_poly::trim(&[&k - q(1)]));
        assert_eq!(epsilon, vec![q(1)]);
        // B_2 = k - 1 and P_2 = 1 + k x; the parent rows are the 1/k-Eulerian
        // polynomials.
        let (parent, first, parent_seeds) = lift.certificate.parent.to_recurrence_parts().unwrap();
        assert_eq!(first, 2);
        assert_eq!(
            parent_seeds,
            vec![polytool::vec_poly::trim(&[q(1), k.clone()])]
        );
        let canonical = one_over_k_eulerian_polynomials_rational(14, &k);
        let parents = parent.generate_rows_rational(&parent_seeds, 2, 12).unwrap();
        for (offset, row) in parents.iter().enumerate() {
            assert_eq!(polytool::vec_poly::trim(row), canonical[2 + offset]);
        }
        // Independent replay from JSON.
        let json = serde_json::to_string(&lift.certificate).unwrap();
        let certificate: SymmetricCompanionCertificate = serde_json::from_str(&json).unwrap();
        let replay = verify_symmetric_companion_certificate(&certificate, 12).unwrap();
        assert!(replay.is_certified() && replay.theorem_applies());
    }
}

#[test]
fn report_json_round_trips() {
    let k = qr(3, 2);
    let report = search(&k, &seeds(&k), &SymmetricCompanionLiftOptions::default());
    let json = report.to_json_string().unwrap();
    assert!(json.contains("\"schema\": \"polytool.recurrence-symmetric-companion.v1\""));
    assert!(json.contains("\"identity_status\": \"certified\""));
    assert!(json.contains("\"real_rootedness_status\": \"theorem_applies\""));
    let decoded: SymmetricCompanionLiftReport = serde_json::from_str(&json).unwrap();
    assert_eq!(decoded, report);
}

#[test]
fn near_misses_are_not_certified_without_a_valid_certificate() {
    let k = q(2);
    // Wrong schedule.
    let options = SymmetricCompanionLiftOptions {
        symmetry_degree_at_first_index: Some(3),
        ..SymmetricCompanionLiftOptions::default()
    };
    let report = search(&k, &seeds(&k), &options);
    assert!(report.certified_lifts().next().is_none());
    // Constant terms that break R_n(0) = 0.
    let altered = vec![vec![q(1), q(1)], vec![q(3), q(7), q(3)]];
    let report = search(&k, &altered, &SymmetricCompanionLiftOptions::default());
    assert!(report.lifts.is_empty());
    assert!(report
        .search
        .attempts
        .iter()
        .all(|a| a.outcome == SymmetricCompanionAttemptOutcome::ConstantTermsIncompatible));
    // Exhausted budget.
    let options = SymmetricCompanionLiftOptions {
        max_linearized_monomials: 2,
        ..SymmetricCompanionLiftOptions::default()
    };
    let report = search(&k, &seeds(&k), &options);
    assert!(report.lifts.is_empty() && report.search.exhausted_budget);
    // A tampered certificate is refused.
    let report = search(&k, &seeds(&k), &SymmetricCompanionLiftOptions::default());
    let mut certificate = report.proving_lift().unwrap().certificate.clone();
    certificate.parameters.epsilon = vec!["2".into()];
    assert!(verify_symmetric_companion_certificate(&certificate, 6).is_err());
}
