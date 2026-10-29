//! Tests of `recurrence::check_three_term_recurrence_theorem` and
//! `recurrence::check_liu_wang_theorem`: classical orthogonal families,
//! derivative families, near misses with the failing index as witness,
//! wrong seeds, and replay.

use num_bigint::BigInt;
use polytool::real_rootedness::RealInterval;
use polytool::recurrence::{
    check_liu_wang_theorem, check_three_term_recurrence_theorem, BigRational, BivarPoly,
    LiftHypothesisStatus, RealRootednessTheoremOptions, RealRootednessTheoremReport,
    RealRootednessTheoremStatus, Recurrence, RecurrenceTerm, RowInterlacing,
};

fn q(value: i64) -> BigRational {
    BigRational::from_integer(BigInt::from(value))
}

fn rows(values: &[&[i64]]) -> Vec<Vec<BigRational>> {
    values
        .iter()
        .map(|row| row.iter().map(|&v| q(v)).collect())
        .collect()
}

/// `coeffs[i][j]` = coefficient of `n^i x^j`.
fn bivar(coeffs: &[&[i64]]) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(rows(coeffs))
}

/// `f(n) P_n = sum c_k(n, x) D^d P_(n-k)` from `(k, d, c_k)`.
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

fn region(lower: Option<i64>, upper: Option<i64>) -> RealRootednessTheoremOptions {
    RealRootednessTheoremOptions {
        root_region: Some(RealInterval {
            lower: lower.map(q),
            lower_closed: false,
            upper: upper.map(q),
            upper_closed: false,
        }),
        ..RealRootednessTheoremOptions::default()
    }
}

fn no_region() -> RealRootednessTheoremOptions {
    RealRootednessTheoremOptions::default()
}

fn assert_replays(report: &RealRootednessTheoremReport) {
    report.verify().expect("the report replays");
    let json = serde_json::to_string(report).unwrap();
    let back: RealRootednessTheoremReport = serde_json::from_str(&json).unwrap();
    assert_eq!(&back, report);
    back.verify().unwrap();
    let value: serde_json::Value = serde_json::from_str(&json).unwrap();
    assert_eq!(
        value["schema"],
        "polytool.recurrence-real-rootedness-theorem.v1"
    );
    assert_eq!(value["lean_replay"], "not_attempted");
    // A tampered status or hypothesis is refused.
    let mut tampered = report.clone();
    tampered.status = if report.is_proved() {
        RealRootednessTheoremStatus::Unresolved
    } else {
        RealRootednessTheoremStatus::Proved
    };
    assert!(tampered.verify().is_err());
    if let Some(first) = report.hypotheses.first() {
        let mut tampered = report.clone();
        tampered.hypotheses[0].status = match first.status {
            LiftHypothesisStatus::Verified => LiftHypothesisStatus::Unresolved,
            _ => LiftHypothesisStatus::Verified,
        };
        assert!(tampered.verify().is_err());
    }
}

/// The witness of the failed hypothesis `name` (the report must have failed).
fn failure(report: &RealRootednessTheoremReport, name: &str) -> Option<BigInt> {
    assert_eq!(
        report.status,
        RealRootednessTheoremStatus::HypothesisFailed,
        "{report:#?}"
    );
    let hypothesis = report
        .hypotheses
        .iter()
        .find(|h| h.name == name)
        .unwrap_or_else(|| panic!("no hypothesis {name}: {report:#?}"));
    match &hypothesis.status {
        LiftHypothesisStatus::Failed { witness } => witness.clone(),
        other => panic!("{name} did not fail: {other:?}"),
    }
}

/// The name of the first failed hypothesis.
fn first_failure(report: &RealRootednessTheoremReport) -> String {
    report
        .hypotheses
        .iter()
        .find(|h| matches!(h.status, LiftHypothesisStatus::Failed { .. }))
        .map(|h| h.name.clone())
        .expect("a failed hypothesis")
}

fn three_term(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[&[i64]],
    options: &RealRootednessTheoremOptions,
) -> RealRootednessTheoremReport {
    let report = check_three_term_recurrence_theorem(rule, first_index, &rows(seeds), options)
        .expect("the check runs");
    assert_replays(&report);
    report
}

fn liu_wang(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[&[i64]],
    options: &RealRootednessTheoremOptions,
) -> RealRootednessTheoremReport {
    let report =
        check_liu_wang_theorem(rule, first_index, &rows(seeds), options).expect("the check runs");
    assert_replays(&report);
    report
}

/// Hermite: `H_n = 2 x H_(n-1) - 2 (n - 1) H_(n-2)`, `H_0 = 1`, `H_1 = 2 x`.
fn hermite() -> Recurrence {
    rule(&[(1, 0, &[&[0, 2]]), (2, 0, &[&[2], &[-2]])], None)
}

/// Laguerre: `n L_n = (2n - 1 - x) L_(n-1) - (n - 1) L_(n-2)`, `L_0 = 1`,
/// `L_1 = 1 - x`; with `sign = -1` the family `L_n(-x)`.
fn laguerre(sign: i64) -> Recurrence {
    rule(
        &[(1, 0, &[&[-1, -sign], &[2]]), (2, 0, &[&[1], &[-1]])],
        Some(&[&[0], &[1]]),
    )
}

#[test]
fn hermite_is_real_rooted_with_strict_interlacing() {
    let report = three_term(&hermite(), 0, &[&[1], &[0, 2]], &no_region());
    assert_eq!(report.status, RealRootednessTheoremStatus::Proved);
    let conclusion = report.conclusion.as_ref().unwrap();
    assert_eq!(conclusion.interlacing, Some(RowInterlacing::Strict));
    assert_eq!(conclusion.interlacing_from, Some(1));
    assert!(conclusion.simple_roots);
    assert!(conclusion.root_region.is_none());

    // Hermite rows have positive roots: the region (-inf, 0) fails at the
    // first seed with a root >= 0 (H_1 = 2 x).
    let report = three_term(&hermite(), 0, &[&[1], &[0, 2]], &region(None, Some(0)));
    assert_eq!(failure(&report, "root_upper_bound"), Some(1.into()));
    assert_eq!(first_failure(&report), "root_upper_bound");
    // ... and the region (-inf, 1) at the first computed row with a root >= 1.
    let report = three_term(&hermite(), 0, &[&[1], &[0, 2]], &region(None, Some(1)));
    // H_2 = 4x^2 - 2 has roots +-0.71, H_3 roots 0, +-1.22.
    assert_eq!(failure(&report, "root_upper_bound"), Some(3.into()));
}

/// Seed rows with trailing zero coefficients (as in some catalog records)
/// give the same verdict as the trimmed rows, and the report replays.
#[test]
fn seeds_with_trailing_zeros_replay() {
    let padded = three_term(&hermite(), 0, &[&[1, 0], &[0, 2, 0]], &no_region());
    let trimmed = three_term(&hermite(), 0, &[&[1], &[0, 2]], &no_region());
    assert_eq!(padded.status, RealRootednessTheoremStatus::Proved);
    assert_eq!(padded, trimmed);
    let padded = liu_wang(&type_b(&[1]), 0, &[&[1, 0, 0]], &no_region());
    let trimmed = liu_wang(&type_b(&[1]), 0, &[&[1]], &no_region());
    assert_eq!(padded, trimmed);
}

#[test]
fn laguerre_roots_are_positive_and_the_reflection_negative() {
    let report = three_term(&laguerre(1), 0, &[&[1], &[1, -1]], &region(Some(0), None));
    assert_eq!(
        report.status,
        RealRootednessTheoremStatus::Proved,
        "{report:#?}"
    );
    let bound = report
        .hypotheses
        .iter()
        .find(|h| h.name == "root_lower_bound")
        .unwrap();
    assert!(bound.evidence.contains("exact ratio"), "{}", bound.evidence);

    let report = three_term(&laguerre(-1), 0, &[&[1], &[1, 1]], &region(None, Some(0)));
    assert_eq!(
        report.status,
        RealRootednessTheoremStatus::Proved,
        "{report:#?}"
    );
    assert_eq!(
        report.conclusion.unwrap().root_region.unwrap().to_string(),
        "(-inf, 0)"
    );
}

#[test]
fn chebyshev_and_legendre_roots_lie_in_the_open_unit_interval() {
    // Chebyshev T: T_n = 2 x T_(n-1) - T_(n-2), T_0 = 1, T_1 = x.
    let chebyshev = rule(&[(1, 0, &[&[0, 2]]), (2, 0, &[&[-1]])], None);
    let report = three_term(&chebyshev, 0, &[&[1], &[0, 1]], &region(Some(-1), Some(1)));
    assert_eq!(
        report.status,
        RealRootednessTheoremStatus::Proved,
        "{report:#?}"
    );
    // Legendre: n P_n = (2n - 1) x P_(n-1) - (n - 1) P_(n-2).
    let legendre = rule(
        &[(1, 0, &[&[0, -1], &[0, 2]]), (2, 0, &[&[1], &[-1]])],
        Some(&[&[0], &[1]]),
    );
    let report = three_term(&legendre, 0, &[&[1], &[0, 1]], &region(Some(-1), Some(1)));
    assert_eq!(
        report.status,
        RealRootednessTheoremStatus::Proved,
        "{report:#?}"
    );
    // A smaller interval fails at an exact row: the largest root of P_3 is
    // sqrt(3/5) > 3/4 ... but the region (-1/2, 1/2) already fails at P_2
    // (roots +-1/sqrt(3)).
    let half = RealRootednessTheoremOptions {
        root_region: Some(RealInterval {
            lower: Some(BigRational::new((-1).into(), 2.into())),
            lower_closed: false,
            upper: Some(BigRational::new(1.into(), 2.into())),
            upper_closed: false,
        }),
        ..RealRootednessTheoremOptions::default()
    };
    let report = three_term(&legendre, 0, &[&[1], &[0, 1]], &half);
    assert_eq!(failure(&report, "root_upper_bound"), Some(2.into()));
    assert_eq!(first_failure(&report), "root_upper_bound");
}

#[test]
fn narayana_is_proved_after_the_reciprocal_coordinate_change() {
    let (rule, seeds, first_row) = polytool::oeis::by_id("A001263")
        .unwrap()
        .recurrence_export_parts()
        .unwrap();
    let first_index = usize::try_from(first_row).unwrap();
    let report =
        check_three_term_recurrence_theorem(&rule, first_index, &seeds, &region(None, Some(0)))
            .unwrap();
    assert_eq!(
        report.status,
        RealRootednessTheoremStatus::Proved,
        "{report:#?}"
    );
    assert!(report.normal_form.contains("Q_n(y) = y^(d_n) P_n(1/y + 1)"));
    let conclusion = report.conclusion.as_ref().unwrap();
    assert_eq!(conclusion.interlacing, Some(RowInterlacing::Strict));
    assert_eq!(
        conclusion.root_region.as_ref().unwrap().to_string(),
        "(-inf, 0)"
    );
    assert_replays(&report);

    // Without a region only Q interlaces; the region containing the pole
    // x = 1 cannot be pulled back.
    let report =
        check_three_term_recurrence_theorem(&rule, first_index, &seeds, &no_region()).unwrap();
    assert_eq!(report.status, RealRootednessTheoremStatus::Proved);
    assert_eq!(report.conclusion.as_ref().unwrap().interlacing, None);
    let report =
        check_three_term_recurrence_theorem(&rule, first_index, &seeds, &region(None, Some(2)))
            .unwrap();
    assert_eq!(report.status, RealRootednessTheoremStatus::Unresolved);
    assert_replays(&report);
}

#[test]
fn three_term_near_misses_name_the_failing_index() {
    // c_n = 2 (7 - n) vanishes at n = 7.
    let vanishing = rule(&[(1, 0, &[&[0, 2]]), (2, 0, &[&[-14], &[2]])], None);
    let report = three_term(&vanishing, 0, &[&[1], &[0, 2]], &no_region());
    assert_eq!(
        failure(&report, "recurrence_coefficient_positive"),
        Some(7.into())
    );
    assert_eq!(first_failure(&report), "recurrence_coefficient_positive");

    // a_n = n - 5 vanishes at n = 5.
    let degenerate = rule(&[(1, 0, &[&[0, -5], &[0, 1]]), (2, 0, &[&[-1]])], None);
    let report = three_term(&degenerate, 0, &[&[1], &[0, 1]], &no_region());
    assert_eq!(
        failure(&report, "leading_coefficient_nonzero"),
        Some(5.into())
    );
    assert_eq!(first_failure(&report), "leading_coefficient_nonzero");

    // Seeds x and (x - 1)(x - 3): the root 0 does not separate 1 and 3.
    let report = three_term(&hermite(), 0, &[&[0, 1], &[3, -4, 1]], &no_region());
    assert_eq!(failure(&report, "seeds_interlace"), Some(1.into()));
    assert_eq!(first_failure(&report), "seeds_interlace");
    // Seeds of degrees 0 and 2.
    let report = three_term(&hermite(), 0, &[&[1], &[-1, 0, 1]], &no_region());
    assert_eq!(failure(&report, "seeds_interlace"), Some(1.into()));
    assert_eq!(first_failure(&report), "seeds_interlace");
    // Seeds that interlace, with higher degrees: x and x^2 - 1.
    let report = three_term(&hermite(), 0, &[&[0, 1], &[-1, 0, 1]], &no_region());
    assert_eq!(report.status, RealRootednessTheoremStatus::Proved);

    // An earlier seed that is not real-rooted (three seeds, start at 1).
    let report = three_term(&hermite(), 0, &[&[1, 0, 1], &[1], &[0, 2]], &no_region());
    assert_eq!(failure(&report, "earlier_seeds"), Some(0.into()));
    assert_eq!(first_failure(&report), "earlier_seeds");
}

#[test]
fn three_term_shapes_that_do_not_apply() {
    let first_order = rule(&[(1, 0, &[&[1, 1]])], None);
    let report = three_term(&first_order, 0, &[&[1], &[1, 1]], &no_region());
    assert_eq!(report.status, RealRootednessTheoremStatus::NotApplicable);
    // c_2 = -(1 + x): not gamma(n) times a constant or a square.
    let linear_tail = rule(&[(1, 0, &[&[0, 2]]), (2, 0, &[&[-1, -1]])], None);
    let report = three_term(&linear_tail, 0, &[&[1], &[0, 2]], &no_region());
    assert_eq!(report.status, RealRootednessTheoremStatus::NotApplicable);
    let report = three_term(&hermite(), 0, &[&[1]], &no_region());
    assert_eq!(report.status, RealRootednessTheoremStatus::NotApplicable);
}

/// `B_n = (1 + (2n - 1) x) B_(n-1) + 2 x (1 - x) B'_(n-1)`, `B_0 = 1`, with
/// `a_n(0) = constant(n)`.
fn type_b(constant: &[i64]) -> Recurrence {
    let mut a: Vec<Vec<i64>> = vec![vec![0, -1], vec![0, 2]];
    for (i, c) in constant.iter().enumerate() {
        if a.len() <= i {
            a.push(vec![0, 0]);
        }
        a[i][0] = *c;
    }
    let a: Vec<&[i64]> = a.iter().map(Vec::as_slice).collect();
    rule(&[(1, 0, &a), (1, 1, &[&[0, 2, -2]])], None)
}

#[test]
fn liu_wang_proves_eulerian_type_b_and_factorial_families() {
    let report = liu_wang(&type_b(&[1]), 0, &[&[1]], &no_region());
    assert_eq!(
        report.status,
        RealRootednessTheoremStatus::Proved,
        "{report:#?}"
    );
    let conclusion = report.conclusion.as_ref().unwrap();
    assert_eq!(conclusion.interlacing, Some(RowInterlacing::Weak));
    assert_eq!(
        conclusion.root_region.as_ref().unwrap().to_string(),
        "(-inf, 0)"
    );

    for id in ["A060187", "A008292"] {
        let (rule, seeds, first_row) = polytool::oeis::by_id(id)
            .unwrap()
            .recurrence_export_parts()
            .unwrap();
        let report = check_liu_wang_theorem(
            &rule,
            usize::try_from(first_row).unwrap(),
            &seeds,
            &region(None, Some(0)),
        )
        .unwrap();
        assert_eq!(
            report.status,
            RealRootednessTheoremStatus::Proved,
            "{id}: {report:#?}"
        );
        assert_replays(&report);
    }

    // P_n = n (1 + x) P_(n-1): no derivative term.
    let factorial = rule(&[(1, 0, &[&[0, 0], &[1, 1]])], None);
    let report = liu_wang(&factorial, 0, &[&[1]], &no_region());
    assert_eq!(
        report.status,
        RealRootednessTheoremStatus::Proved,
        "{report:#?}"
    );
}

#[test]
fn liu_wang_near_misses_name_the_failing_index() {
    // a_n(0) = 6 - n vanishes at n = 6.
    let report = liu_wang(&type_b(&[6, -1]), 0, &[&[1]], &no_region());
    assert_eq!(failure(&report, "constant_term_positive"), Some(6.into()));
    assert_eq!(first_failure(&report), "constant_term_positive");

    // b_n = 2 (4 - n) x (1 - x) is positive on (-inf, 0) from n = 5 on; the
    // leading part a_n = 1 + 10 n^2 x keeps lambda_n > 0.
    let late_sign = rule(
        &[
            (1, 0, &[&[1, 0], &[0, 0], &[0, 10]]),
            (1, 1, &[&[0, 8, -8], &[0, -2, 2]]),
        ],
        None,
    );
    let report = liu_wang(&late_sign, 0, &[&[1]], &no_region());
    assert_eq!(
        failure(&report, "derivative_coefficient_nonpositive"),
        Some(5.into())
    );
    assert_eq!(first_failure(&report), "derivative_coefficient_nonpositive");

    // lambda_n = (2n - 1) - 2 (n - 1) - [n >= 3] ... : a_n = 1 + (4 - n) x with
    // no derivative term has lambda_n = 4 - n, failing at n = 4.
    let shrinking = rule(&[(1, 0, &[&[1, 4], &[0, -1]])], None);
    let report = liu_wang(&shrinking, 0, &[&[1]], &no_region());
    assert_eq!(failure(&report, "degree_step_positive"), Some(4.into()));
    assert_eq!(first_failure(&report), "degree_step_positive");

    // b_n(0) = 1 != 0 (lambda_n = 3 - 2 n fails too, at n = 2).
    let shifted = rule(&[(1, 0, &[&[1, 1]]), (1, 1, &[&[1, 2, -2]])], None);
    let report = liu_wang(&shifted, 0, &[&[1]], &no_region());
    assert_eq!(
        failure(&report, "derivative_coefficient_vanishes_at_zero"),
        Some(1.into())
    );

    // A seed with a positive root.
    let report = liu_wang(&type_b(&[1]), 0, &[&[1, -1]], &no_region());
    assert_eq!(failure(&report, "seed_roots_negative"), Some(0.into()));
    assert_eq!(first_failure(&report), "seed_roots_negative");

    // A requested region that does not contain (-inf, 0).
    let report = liu_wang(&type_b(&[1]), 0, &[&[1]], &region(Some(-1), Some(0)));
    assert_eq!(report.status, RealRootednessTheoremStatus::Unresolved);
}

#[test]
fn liu_wang_shapes_that_do_not_apply() {
    let report = liu_wang(&hermite(), 0, &[&[1], &[0, 2]], &no_region());
    assert_eq!(report.status, RealRootednessTheoremStatus::NotApplicable);
    let second_derivative = rule(&[(1, 0, &[&[1, 1]]), (1, 2, &[&[0, 1]])], None);
    let report = liu_wang(&second_derivative, 0, &[&[1]], &no_region());
    assert_eq!(report.status, RealRootednessTheoremStatus::NotApplicable);
    let cubic = rule(&[(1, 0, &[&[1, 1]]), (1, 1, &[&[0, 1, 0, 1]])], None);
    let report = liu_wang(&cubic, 0, &[&[1]], &no_region());
    assert_eq!(report.status, RealRootednessTheoremStatus::NotApplicable);
}

/// `A_n = (1 + (n - 1) x) A_(n-1) + x (1 - x) A'_(n-1)` (Eulerian).
fn eulerian() -> Recurrence {
    rule(
        &[(1, 0, &[&[1, -1], &[0, 1]]), (1, 1, &[&[0, 1, -1]])],
        None,
    )
}

fn without_start_delay() -> RealRootednessTheoremOptions {
    RealRootednessTheoremOptions {
        maximum_start_delay: 0,
        ..RealRootednessTheoremOptions::default()
    }
}

#[test]
fn liu_wang_moves_the_induction_start_past_a_degree_plateau() {
    // From A_0 = 1 the degree does not grow at n = 1 (A_1 = 1).
    let fixed = liu_wang(&eulerian(), 0, &[&[1]], &without_start_delay());
    assert_eq!(failure(&fixed, "degree_step_positive"), Some(1.into()));
    assert_eq!(fixed.start_delay, 0);

    let report = liu_wang(&eulerian(), 0, &[&[1]], &no_region());
    assert_eq!(
        report.status,
        RealRootednessTheoremStatus::Proved,
        "{report:#?}"
    );
    assert_eq!(report.start_delay, 1);
    // The record keeps the supplied seed; the theorem starts at A_1 and
    // checks A_0 exactly.
    assert_eq!(report.recurrence.initial_polynomials.len(), 1);
    assert!(report.normal_form.contains("for n >= 2"));
    assert!(report
        .hypotheses
        .iter()
        .any(|h| h.name == "earlier_seeds" && h.status == LiftHypothesisStatus::Verified));
    let conclusion = report.conclusion.as_ref().unwrap();
    assert_eq!(conclusion.first_index, 0);
    // A_0 = 1 and A_1 = 1 interlace (constants), checked exactly.
    assert_eq!(conclusion.interlacing_from, Some(1));
    assert!(conclusion.statement.contains("start moved 1 row(s)"));
    let json = serde_json::to_value(&report).unwrap();
    assert_eq!(json["start_delay"], 1);
    // The recorded delay is replayed; a different one is refused.
    let mut tampered = report.clone();
    tampered.start_delay = 2;
    assert!(tampered.verify().is_err());
    // Reports without a delay do not serialize the field.
    let json = serde_json::to_value(&fixed).unwrap();
    assert!(json.get("start_delay").is_none());
}

#[test]
fn three_term_moves_the_induction_start_past_equal_degree_seeds() {
    // P_n = 2 x P_(n-1) - P_(n-2), P_0 = P_1 = 1: P_n = V_(n-1) (Chebyshev
    // polynomials of the third kind) from n = 1 on.  The seeds have equal
    // degrees, so the theorem starts at P_1, P_2; P_0 is checked exactly.
    let chebyshev = rule(&[(1, 0, &[&[0, 2]]), (2, 0, &[&[-1]])], None);
    let fixed = three_term(&chebyshev, 0, &[&[1], &[1]], &without_start_delay());
    assert_eq!(failure(&fixed, "seeds_interlace"), Some(1.into()));

    let report = three_term(&chebyshev, 0, &[&[1], &[1]], &region(Some(-1), Some(1)));
    assert_eq!(
        report.status,
        RealRootednessTheoremStatus::Proved,
        "{report:#?}"
    );
    assert_eq!(report.start_delay, 1);
    let conclusion = report.conclusion.as_ref().unwrap();
    assert_eq!(conclusion.first_index, 0);
    // P_0 = 1 does not strictly interlace P_1 = 1 (equal degrees), so the
    // interlacing claim stays at the moved start.
    assert_eq!(conclusion.interlacing_from, Some(2));
    assert_eq!(
        conclusion.root_region.as_ref().unwrap().to_string(),
        "(-1, 1)"
    );
}

#[test]
fn a_later_start_does_not_help_when_an_early_row_fails() {
    // The seed 1 + x + x^2 is not real-rooted; every moved start keeps it as
    // a checked seed row, so the report at the supplied seed is returned.
    let report = liu_wang(&eulerian(), 0, &[&[1, 1, 1]], &no_region());
    assert_eq!(report.start_delay, 0);
    assert_eq!(failure(&report, "seed_roots_negative"), Some(0.into()));

    // A hypothesis failing for every large n: a_n(0) = 6 - n.
    let report = liu_wang(&type_b(&[6, -1]), 0, &[&[1]], &no_region());
    assert_eq!(report.start_delay, 0);
    assert_eq!(failure(&report, "constant_term_positive"), Some(6.into()));
}

/// The example of `doc/library.md`.
#[test]
fn library_documentation_example() -> Result<(), polytool::recurrence::RecurrenceError> {
    let q = |n: i64| BigRational::from_integer(n.into());
    // Legendre: n P_n = (2n - 1) x P_(n-1) - (n - 1) P_(n-2), P_0 = 1, P_1 = x.
    let legendre = Recurrence {
        terms: vec![
            RecurrenceTerm::new(
                1,
                0,
                BivarPoly {
                    coeffs: vec![vec![q(0), q(-1)], vec![q(0), q(2)]],
                },
            ),
            RecurrenceTerm::new(
                2,
                0,
                BivarPoly {
                    coeffs: vec![vec![q(1)], vec![q(-1)]],
                },
            ),
        ],
        denominator: Some(BivarPoly {
            coeffs: vec![vec![q(0)], vec![q(1)]],
        }),
        inhomogeneous: None,
    };
    let options = RealRootednessTheoremOptions {
        root_region: Some(RealInterval::closed(q(-1), q(1))),
        ..RealRootednessTheoremOptions::default()
    };
    let report = check_three_term_recurrence_theorem(
        &legendre,
        0,
        &[vec![q(1)], vec![q(0), q(1)]],
        &options,
    )?;
    assert!(report.is_proved());
    let conclusion = report.conclusion.as_ref().unwrap();
    assert_eq!(conclusion.interlacing, Some(RowInterlacing::Strict));
    assert_eq!(
        conclusion.root_region.as_ref().unwrap().to_string(),
        "(-1, 1)"
    );
    report.verify()?;
    Ok(())
}
