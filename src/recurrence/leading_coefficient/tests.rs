//! Tests of the degree-law and leading-coefficient certificates.

use super::*;
use crate::recurrence::{BivarPoly, RecurrenceTerm};

fn q(value: i64) -> Q {
    integer(value)
}

/// `c(n, x)` from `(n power, x power, value)` triples.
fn bivariate(entries: &[(usize, usize, i64)]) -> BivarPoly {
    let rows = entries.iter().map(|e| e.0).max().unwrap_or(0) + 1;
    let width = entries.iter().map(|e| e.1).max().unwrap_or(0) + 1;
    let mut coefficients = vec![vec![Q::zero(); width]; rows];
    for &(i, j, value) in entries {
        coefficients[i][j] += q(value);
    }
    BivarPoly::from_coefficients_normalized(coefficients)
}

/// `(n power, x power, value)` entries of a coefficient.
type Entries<'a> = &'a [(usize, usize, i64)];

/// A rule from `(lag, derivative order, coefficient entries)` terms.
fn rule(terms: &[(usize, usize, Entries<'_>)]) -> Recurrence {
    Recurrence {
        terms: terms
            .iter()
            .map(|(lag, order, entries)| RecurrenceTerm::new(*lag, *order, bivariate(entries)))
            .collect(),
        denominator: None,
        inhomogeneous: None,
    }
}

fn rows(values: &[&[i64]]) -> Vec<Vec<Q>> {
    values
        .iter()
        .map(|row| row.iter().map(|&v| q(v)).collect())
        .collect()
}

fn certify(rule: &Recurrence, first: usize, seeds: &[Vec<Q>]) -> LeadingCoefficientCertificate {
    let certificate =
        certify_leading_coefficient(rule, first, seeds, &LeadingCoefficientOptions::default())
            .unwrap();
    certificate.verify().unwrap();
    let json = serde_json::to_string(&certificate).unwrap();
    let back: LeadingCoefficientCertificate = serde_json::from_str(&json).unwrap();
    assert_eq!(back, certificate);
    certificate
}

fn catalog(id: &str) -> (Recurrence, usize, Vec<Vec<Q>>) {
    let entry = crate::oeis::by_id(id).unwrap();
    let (rule, seeds) = entry.recurrence_parts().unwrap();
    (rule, entry.recurrence_first_index, seeds)
}

/// The degrees and signs of the certificate agree with computed rows.
fn assert_matches_rows(certificate: &LeadingCoefficientCertificate, rule: &Recurrence) {
    let (_, first, seeds) = certificate.recurrence.to_recurrence_parts().unwrap();
    let computed = rule.generate_rows_rational(&seeds, first, 60).unwrap();
    let from = certificate.holds_from.unwrap();
    for (i, row) in computed.iter().enumerate() {
        let n = first + i;
        if n < from {
            continue;
        }
        let expected_sign = leading_coefficient(row).and_then(LeadingCoefficientSign::of);
        assert_eq!(certificate.sign_at(n), expected_sign, "sign of P_{n}");
        if certificate.degree_is_exact {
            assert_eq!(certificate.degree_at(n), row_degree(row), "degree of P_{n}");
        }
    }
}

#[test]
fn half_growth_associated_stirling_numbers() {
    let (rule, first, seeds) = catalog("A008299");
    let certificate = certify(&rule, first, &seeds);
    assert!(certificate.is_proved(), "{certificate:#?}");
    let law = certificate.degree_law.as_ref().unwrap();
    assert_eq!((law.period, law.growth_per_period), (2, 1));
    assert_matches_rows(&certificate, &rule);
}

#[test]
fn periodic_degree_law_of_period_three() {
    // A102547: P_n = P_(n-1) + x P_(n-3), deg P_n = floor(n / 3).
    let (rule, first, seeds) = catalog("A102547");
    let certificate = certify(&rule, first, &seeds);
    assert!(certificate.is_proved());
    let law = certificate.degree_law.as_ref().unwrap();
    assert_eq!((law.period, law.growth_per_period), (3, 1));
    assert_eq!(certificate.route, Some(LeadingCoefficientRoute::Cone));
    assert_eq!(
        certificate.lean_declarations,
        ["RealRooted.LinRec.natDegree_eq_and_leadingCoeff_pos_of_degreeLaw"]
    );
    assert_eq!(certificate.degree_at(30), Some(10));
    assert_matches_rows(&certificate, &rule);
}

#[test]
fn veronese_section_residue_classes() {
    // P_n = (1 + x) P_(n-2) with P_0 = 1, P_1 = 1 + 2x + 3x^2: the classes
    // evolve separately, deg P_(2k) = k, deg P_(2k+1) = k + 2.
    let rule = rule(&[(2, 0, &[(0, 0, 1), (0, 1, 1)])]);
    let certificate = certify(&rule, 0, &rows(&[&[1], &[1, 2, 3]]));
    assert!(certificate.is_proved());
    assert_eq!(certificate.route, Some(LeadingCoefficientRoute::TwoStep));
    let law = certificate.degree_law.as_ref().unwrap();
    assert_eq!(law.period, 2);
    assert_eq!(certificate.degree_at(11), Some(7));
    assert_matches_rows(&certificate, &rule);

    // With an index-dependent factor the two-step theorem does not apply,
    // but the residue classes still evolve separately.
    let rule = super::tests::rule(&[(2, 0, &[(0, 0, 1), (1, 1, 1)])]);
    let certificate = certify(&rule, 0, &rows(&[&[1], &[1, 2, 3]]));
    assert!(certificate.is_proved());
    assert_eq!(certificate.route, Some(LeadingCoefficientRoute::Subseq));
    assert_matches_rows(&certificate, &rule);
}

#[test]
fn two_step_product_with_negative_factor_alternates_in_pairs() {
    // P_n = (1 - x) P_(n-2): the signs repeat with period four.
    let rule = rule(&[(2, 0, &[(0, 0, 1), (0, 1, -1)])]);
    let certificate = certify(&rule, 0, &rows(&[&[1], &[0, 1]]));
    assert!(certificate.is_proved());
    let law = certificate.degree_law.as_ref().unwrap();
    assert_eq!(law.period, 4);
    assert!(certificate
        .lean_declarations
        .contains(&"RealRooted.hasPosLeadingCoeff_C_mul".to_string()));
    assert_matches_rows(&certificate, &rule);
}

#[test]
fn ratio_route_for_a_negative_multiplier() {
    // The RealRooted example: P_n = 2x P_(n-1) - x^3 P_(n-3), top
    // coefficients c(n) = 2 c(n-1) - c(n-3).
    let rule = rule(&[(1, 0, &[(0, 1, 2)]), (3, 0, &[(0, 3, -1)])]);
    let certificate = certify(&rule, 0, &rows(&[&[1], &[0, 2], &[0, 0, 3]]));
    assert!(certificate.is_proved(), "{:#?}", certificate.attempts);
    assert_eq!(certificate.route, Some(LeadingCoefficientRoute::Recurrence));
    assert!(certificate.ratio.as_ref().unwrap().is_positive());
    assert_matches_rows(&certificate, &rule);
}

#[test]
fn cancellation_of_the_top_term() {
    // P_n = (x - x^2) P'_(n-1) + (1 + x) P_(n-1), P_0 = 1 + x: the x^2
    // terms cancel, every row has degree one.
    let rule = rule(&[
        (1, 1, &[(0, 1, 1), (0, 2, -1)]),
        (1, 0, &[(0, 0, 1), (0, 1, 1)]),
    ]);
    let certificate = certify(&rule, 0, &rows(&[&[1, 1]]));
    assert!(certificate.is_proved(), "{:#?}", certificate.attempts);
    assert_eq!(certificate.route, Some(LeadingCoefficientRoute::Cancel));
    assert_eq!(certificate.degree_at(40), Some(1));
    // The top-coefficient recurrence does not hold where the top terms
    // cancel: no constant leading coefficient is claimed (P_0 = 1 + x,
    // P_1 = 1 + 3x).
    assert_eq!(certificate.exact_leading_coefficient, None);
    assert_matches_rows(&certificate, &rule);
}

#[test]
fn nonnegative_coefficients_within_an_envelope() {
    // Eulerian polynomials A_n = (1 + (n - 1) x) A_(n-1) + x (1 - x) A'_(n-1)
    // (first index 1), by the coefficients route with an affine envelope.
    let rule = rule(&[
        (1, 0, &[(0, 0, 1), (1, 1, 1), (0, 1, -1)]),
        (1, 1, &[(0, 1, 1), (0, 2, -1)]),
    ]);
    let seeds = rows(&[&[1]]);
    let law = RecurrenceDegreeLaw {
        start: 1,
        period: 1,
        base_degrees: vec![0],
        growth_per_period: 1,
        leading_coefficient_signs: vec![LeadingCoefficientSign::Positive],
    };
    let certificate = certify_leading_coefficient_route(
        &rule,
        1,
        &seeds,
        &LeadingCoefficientOptions::default(),
        &law,
        LeadingCoefficientRoute::Coefficients,
        None,
    )
    .unwrap();
    assert!(certificate.is_proved(), "{:#?}", certificate.hypotheses);
    assert!(!certificate.degree_is_exact);
    assert_eq!(certificate.degree_at(5), None);
    assert_eq!(
        certificate.sign_at(5),
        Some(LeadingCoefficientSign::Positive)
    );
    certificate.verify().unwrap();
    // The search proves the exact law by the cone.
    let searched = certify(&rule, 1, &seeds);
    assert_eq!(searched.route, Some(LeadingCoefficientRoute::Cone));
    assert_matches_rows(&searched, &rule);
}

#[test]
fn paired_finite_lag_coefficient_bounds() {
    // P_n = (1 + x) P_(n-1) + x^2 P'_(n-2), P_0 = 1, P_1 = 1 + x: lags one
    // and two, a first derivative at lag two; the envelope D(n) = n.
    let rule = rule(&[(1, 0, &[(0, 0, 1), (0, 1, 1)]), (2, 1, &[(0, 2, 1)])]);
    let seeds = rows(&[&[1], &[1, 1]]);
    let law = RecurrenceDegreeLaw {
        start: 0,
        period: 1,
        base_degrees: vec![0],
        growth_per_period: 1,
        leading_coefficient_signs: vec![LeadingCoefficientSign::Positive],
    };
    let certificate = certify_leading_coefficient_route(
        &rule,
        0,
        &seeds,
        &LeadingCoefficientOptions::default(),
        &law,
        LeadingCoefficientRoute::Coefficients,
        None,
    )
    .unwrap();
    assert!(certificate.is_proved(), "{:#?}", certificate.hypotheses);
    assert_eq!(certificate.induction_start, 2);
    assert!(certificate
        .hypotheses
        .iter()
        .any(|h| h.name == "coefficient_multiplier_1_high_lag_2_class_2"));
    assert_eq!(
        certificate.lean_declarations,
        [
            "RealRooted.LinRec.natDegree_le_and_hasNonnegCoeffs",
            "RealRooted.HasNonnegCoeffs.pos_leadingCoeff"
        ]
    );
    certificate.verify().unwrap();
    for n in 0..30 {
        assert_eq!(
            certificate.sign_at(n),
            Some(LeadingCoefficientSign::Positive)
        );
    }
    // An envelope that the seeds exceed fails the initial rows.
    let narrow = RecurrenceDegreeLaw {
        growth_per_period: 0,
        ..law.clone()
    };
    let certificate = certify_leading_coefficient_route(
        &rule,
        0,
        &seeds,
        &LeadingCoefficientOptions::default(),
        &narrow,
        LeadingCoefficientRoute::Coefficients,
        None,
    )
    .unwrap();
    assert!(!certificate.is_proved());
    // A negative constant term multiplier is refused (r = 0 end).
    let negative = super::tests::rule(&[(1, 0, &[(0, 0, -1), (0, 1, 1)]), (2, 1, &[(0, 2, 1)])]);
    let certificate = certify_leading_coefficient_route(
        &negative,
        0,
        &seeds,
        &LeadingCoefficientOptions::default(),
        &law,
        LeadingCoefficientRoute::Coefficients,
        None,
    )
    .unwrap();
    assert!(!certificate.is_proved());
}

#[test]
fn a_vanishing_left_hand_factor_is_not_proved() {
    // (n - 50) P_n = x P_(n-1): the rows exist up to n = 49, but the
    // left-hand factor vanishes at n = 50.
    let mut rule = rule(&[(1, 0, &[(0, 1, 1)])]);
    rule.denominator = Some(bivariate(&[(1, 0, 1), (0, 0, -50)]));
    let certificate = certify(&rule, 0, &rows(&[&[1]]));
    assert_eq!(certificate.status, LeadingCoefficientStatus::NotProved);
    let hypothesis = certificate
        .hypotheses
        .iter()
        .find(|h| h.name == "left_factor_nonzero")
        .unwrap();
    assert_eq!(
        hypothesis.status,
        LiftHypothesisStatus::Failed {
            witness: Some(BigInt::from(50))
        }
    );
    // A factor that vanishes inside the sampled rows: not applicable.
    let mut rule = super::tests::rule(&[(1, 0, &[(0, 1, 1)])]);
    rule.denominator = Some(bivariate(&[(1, 0, 1), (0, 0, -5)]));
    let certificate = certify(&rule, 0, &rows(&[&[1]]));
    assert_eq!(certificate.status, LeadingCoefficientStatus::NotApplicable);
}

#[test]
fn a_vanishing_leader_is_not_proved() {
    // P_n = ((n - 50) x + 1) P_(n-1): the leader vanishes at n = 50.
    let rule = rule(&[(1, 0, &[(0, 0, 1), (1, 1, 1), (0, 1, -50)])]);
    let certificate = certify(&rule, 0, &rows(&[&[1]]));
    assert_eq!(certificate.status, LeadingCoefficientStatus::NotProved);
    assert!(certificate.conclusion.is_none());
    let failure = certificate
        .hypotheses
        .iter()
        .filter_map(|h| match &h.status {
            LiftHypothesisStatus::Failed { witness } => witness.clone(),
            _ => None,
        })
        .min()
        .unwrap();
    assert_eq!(failure, BigInt::from(50));
    assert!(!certificate.attempts.is_empty());
}

#[test]
fn a_zero_row_moves_the_start() {
    // A zero seed: P_0 = 0, P_1 = 1, P_n = x P_(n-1) + P_(n-2).
    let rule = rule(&[(1, 0, &[(0, 1, 1)]), (2, 0, &[(0, 0, 1)])]);
    let certificate = certify(&rule, 0, &rows(&[&[], &[1]]));
    assert!(certificate.is_proved());
    assert_eq!(certificate.holds_from, Some(1));
    assert_eq!(certificate.degree_at(0), None);
    assert_eq!(certificate.degree_at(7), Some(6));
}

#[test]
fn tampered_certificates_fail_to_verify() {
    let (rule, first, seeds) = catalog("A102547");
    let certificate = certify(&rule, first, &seeds);
    let mut tampered = certificate.clone();
    tampered.degree_law.as_mut().unwrap().growth_per_period = 2;
    assert!(tampered.verify().is_err());
    let mut tampered = certificate.clone();
    tampered.route = Some(LeadingCoefficientRoute::Recurrence);
    assert!(tampered.verify().is_err());
    let mut tampered = certificate.clone();
    tampered.hypotheses[1].evidence.push('!');
    assert!(tampered.verify().is_err());
    let mut tampered = certificate.clone();
    tampered.initial_rows[0].degree = Some(3);
    assert!(tampered.verify().is_err());
    let mut tampered = certificate.clone();
    tampered.status = LeadingCoefficientStatus::NotProved;
    assert!(tampered.verify().is_err());
    // A recorded law outside the searched range is refused before replay.
    let mut tampered = certificate.clone();
    tampered.degree_law.as_mut().unwrap().base_degrees[0] = usize::MAX;
    assert!(tampered.verify().is_err());
    let mut tampered = certificate;
    tampered.degree_law.as_mut().unwrap().start = 1 << 40;
    assert!(tampered.verify().is_err());

    // A not-proved certificate claiming a proof is refused.
    let rule = super::tests::rule(&[(1, 0, &[(0, 0, 1), (1, 1, 1), (0, 1, -50)])]);
    let mut forged = certify(&rule, 0, &rows(&[&[1]]));
    forged.status = LeadingCoefficientStatus::Proved;
    assert!(forged.verify().is_err());
}

#[test]
fn unsupported_rules_are_not_applicable() {
    let mut inhomogeneous = rule(&[(1, 0, &[(0, 1, 1)])]);
    inhomogeneous.inhomogeneous = Some(bivariate(&[(0, 0, 1)]));
    let certificate = certify(&inhomogeneous, 0, &rows(&[&[1]]));
    assert_eq!(certificate.status, LeadingCoefficientStatus::NotApplicable);
}

/// The regression families of the issue: degree laws, signs and the exact
/// leading coefficient where the rows are monic.
#[test]
fn catalog_regression_families() {
    for (id, period, growth) in [
        ("A019538", 1, 1),
        ("A062145", 1, 1),
        ("A062196", 1, 1),
        ("A085478", 1, 1),
        ("A281260", 1, 1),
    ] {
        let (rule, first, seeds) = catalog(id);
        let certificate = certify(&rule, first, &seeds);
        assert!(certificate.is_proved(), "{id}: {:#?}", certificate.attempts);
        let law = certificate.degree_law.as_ref().unwrap();
        assert_eq!(
            (law.period, law.growth_per_period),
            (period, growth),
            "{id}"
        );
        assert_matches_rows(&certificate, &rule);
    }
    let (rule, first, seeds) = catalog("A085478");
    let certificate = certify(&rule, first, &seeds);
    assert_eq!(certificate.exact_leading_coefficient, Some(Q::one()));
}

/// A101343: first-order rule whose top terms cancel at every other step with
/// a vanishing next level (`deg P_n = floor((n + 1) / 2)`); the `cancel`
/// route with the coefficient-ratio invariant `[x^(D-1)] P_n >= rho lead P_n`
/// proves it, and without the ratio it is not proved.
#[test]
fn cancellation_with_a_coefficient_ratio_invariant() {
    let (rule, first, seeds) = catalog("A101343");
    let certificate = certify(&rule, first, &seeds);
    assert!(certificate.is_proved(), "{:#?}", certificate.attempts);
    assert_eq!(certificate.route, Some(LeadingCoefficientRoute::Cancel));
    assert!(certificate.ratio.as_ref().unwrap().is_positive());
    assert_eq!(
        certificate.lean_declarations,
        ["RealRooted.LinRec.natDegree_eq_and_leadingCoeff_pos_of_cancel"]
    );
    let law = certificate.degree_law.as_ref().unwrap();
    assert_eq!((law.period, law.growth_per_period), (2, 1));
    assert_eq!(certificate.exact_leading_coefficient, None);
    assert_matches_rows(&certificate, &rule);
    let without_ratio = certify_leading_coefficient_route(
        &rule,
        first,
        &seeds,
        &certificate.options,
        law,
        LeadingCoefficientRoute::Cancel,
        None,
    )
    .unwrap();
    assert!(!without_ratio.is_proved());
    // A larger ratio fails the base row or the invariant; a tampered ratio
    // does not replay.
    let mut tampered = certificate.clone();
    tampered.ratio = Some(certificate.ratio.as_ref().unwrap() * q(64));
    assert!(tampered.verify().is_err());
    assert!(certify_leading_coefficient_route(
        &rule,
        first,
        &seeds,
        &certificate.options,
        law,
        LeadingCoefficientRoute::Cancel,
        Some(&q(-1)),
    )
    .is_err());
}
