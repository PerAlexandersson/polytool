//! Rules with a left-hand factor `x^e f(n)` (as the adaptive search stores
//! them with `try_laurent`): the theorem checks and the lifts normalize them
//! exactly to `f(n) P_n = sum (c / x^e) D^d P_(n-r)` when `x^e` divides the
//! right-hand side, reach the same conclusion as the `x`-free rule, record
//! the normalization and replay it; a rule whose division by `x^e` happens
//! only on the rows is rejected with the reason.

use polytool::oeis;
use polytool::real_rootedness::RealInterval;
use polytool::recurrence::{
    certify_veronese_section_lift, check_liu_wang_derivative_lag_theorem, check_liu_wang_theorem,
    check_ma_wang_interval_theorem, check_three_term_recurrence_theorem, find_row_factor_lifts,
    BigRational, BivarPoly, LeftHandFactorNormalization, RealRootednessTheoremOptions,
    RealRootednessTheoremReport, RealRootednessTheoremStatus, Recurrence, RecurrenceError,
    RecurrenceJson, RecurrenceTerm, RowFactorLiftOptions, VeroneseSectionCandidate,
    VeroneseSectionCertifyOptions, VeroneseSectionGauge, VeroneseSectionVerdict,
};

type Q = BigRational;

mod common;
use common::{q, rows};

fn bivar(coeffs: &[&[i64]]) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(rows(coeffs))
}

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

/// `x^e c(n, x)`.
fn times_x_power(polynomial: &BivarPoly, e: usize) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(
        polynomial
            .coeffs
            .iter()
            .map(|row| {
                let mut shifted = vec![q(0); e];
                shifted.extend(row.iter().cloned());
                shifted
            })
            .collect(),
    )
}

/// Both sides of `rule` multiplied by `x^e`: the left-hand factor becomes
/// `x^e f(n)`, stored in `denominator` as the adaptive search stores it.
fn multiplied_by_x_power(rule: &Recurrence, e: usize) -> Recurrence {
    let one = bivar(&[&[1]]);
    Recurrence {
        terms: rule
            .terms
            .iter()
            .map(|term| RecurrenceTerm {
                coeff: times_x_power(&term.coeff, e),
                ..term.clone()
            })
            .collect(),
        denominator: Some(times_x_power(rule.denominator.as_ref().unwrap_or(&one), e)),
        inhomogeneous: rule.inhomogeneous.as_ref().map(|g| times_x_power(g, e)),
    }
}

fn no_region() -> RealRootednessTheoremOptions {
    RealRootednessTheoremOptions::default()
}

fn open_region(lower: Option<i64>, upper: Option<i64>) -> RealRootednessTheoremOptions {
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

type Check = fn(
    &Recurrence,
    usize,
    &[Vec<Q>],
    &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError>;

/// Run `check` on `base` and on `base` times `x^e`; both are proved with the
/// same conclusion, the second records the normalization, and both replay.
fn assert_same_conclusion(
    check: Check,
    base: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
    e: usize,
) -> RealRootednessTheoremReport {
    let plain = check(base, first_index, seeds, options).unwrap();
    assert!(plain.is_proved(), "{plain:#?}");
    assert!(plain.left_hand_factor_normalization.is_none());
    let lifted_rule = multiplied_by_x_power(base, e);
    // The two rules give the same rows.
    assert_eq!(
        base.generate_rows_rational(seeds, first_index, 10).unwrap(),
        lifted_rule
            .generate_rows_rational(seeds, first_index, 10)
            .unwrap()
    );
    let report = check(&lifted_rule, first_index, seeds, options).unwrap();
    assert!(report.is_proved(), "{report:#?}");
    report.verify().unwrap();
    let json = serde_json::to_string(&report).unwrap();
    let back: RealRootednessTheoremReport = serde_json::from_str(&json).unwrap();
    assert_eq!(back, report);
    back.verify().unwrap();

    // Same conclusion, hypotheses, normal form and delay as the x-free rule.
    let (got, want) = (
        report.conclusion.as_ref().unwrap(),
        plain.conclusion.as_ref().unwrap(),
    );
    assert_eq!(got.first_index, want.first_index);
    assert_eq!(got.interlacing, want.interlacing);
    assert_eq!(got.interlacing_from, want.interlacing_from);
    assert_eq!(got.simple_roots, want.simple_roots);
    assert_eq!(got.root_region, want.root_region);
    assert!(got.statement.starts_with(&want.statement));
    assert!(got.statement.contains("applied to the normalized rule"));
    assert_eq!(report.hypotheses, plain.hypotheses);
    assert_eq!(report.normal_form, plain.normal_form);
    assert_eq!(report.start_delay, plain.start_delay);
    assert_eq!(report.lean_route, plain.lean_route);

    // The supplied rule stays the recorded rule; the record names both.
    assert_eq!(
        report.recurrence.recurrence,
        RecurrenceJson::from_recurrence_rational(&lifted_rule, first_index, seeds, None).recurrence
    );
    assert_eq!(
        report.recurrence,
        plain_record_with_rule(&plain, &lifted_rule)
    );
    let normalization = report.left_hand_factor_normalization.as_ref().unwrap();
    assert_eq!(normalization.x_power, e);
    assert_eq!(normalization.original_rule, report.recurrence.recurrence);
    normalization.verify().unwrap();
    assert_eq!(
        LeftHandFactorNormalization::from_recurrence(&lifted_rule).unwrap(),
        Some(normalization.clone())
    );
    let normalized = normalization.normalized_recurrence().unwrap();
    assert_eq!(
        normalized
            .generate_rows_rational(seeds, first_index, 10)
            .unwrap(),
        base.generate_rows_rational(seeds, first_index, 10).unwrap()
    );
    assert_tampering_is_refused(&report);
    report
}

/// The record of `plain` (seeds in their normal form) with `rule` instead.
fn plain_record_with_rule(
    plain: &RealRootednessTheoremReport,
    rule: &Recurrence,
) -> RecurrenceJson {
    let mut record = plain.recurrence.clone();
    record.recurrence = RecurrenceJson::from_recurrence_rational(rule, 0, &[], None).recurrence;
    record
}

fn assert_tampering_is_refused(report: &RealRootednessTheoremReport) {
    let normalization = report.left_hand_factor_normalization.as_ref().unwrap();

    let mut dropped = report.clone();
    dropped.left_hand_factor_normalization = None;
    assert!(dropped.verify().is_err());

    let mut wrong_power = report.clone();
    wrong_power
        .left_hand_factor_normalization
        .as_mut()
        .unwrap()
        .x_power += 1;
    assert!(wrong_power.verify().is_err());
    let mut record = normalization.clone();
    record.x_power += 1;
    assert!(record.verify().is_err());

    let mut wrong_rule = normalization.clone();
    wrong_rule.normalized_rule.terms[0].coeff.coeffs[0][0] = "17".into();
    assert!(wrong_rule.verify().is_err());
    let mut tampered = report.clone();
    tampered.left_hand_factor_normalization = Some(wrong_rule);
    assert!(tampered.verify().is_err());

    let mut wrong_statement = normalization.clone();
    wrong_statement.statement.push('.');
    assert!(wrong_statement.verify().is_err());

    // A record whose original rule is the x-free rule does not normalize.
    let mut not_original = normalization.clone();
    not_original.original_rule = not_original.normalized_rule.clone();
    assert!(not_original.verify().is_err());
}

fn hermite() -> Recurrence {
    rule(&[(1, 0, &[&[0, 2]]), (2, 0, &[&[2], &[-2]])], None)
}

fn laguerre_reflected() -> Recurrence {
    // n L_n = (2n - 1 + x) L_(n-1) - (n - 1) L_(n-2): the rows L_n(-x).
    rule(
        &[(1, 0, &[&[-1, 1], &[2]]), (2, 0, &[&[1], &[-1]])],
        Some(&[&[0], &[1]]),
    )
}

#[test]
fn three_term_rules_times_x_power_are_proved_as_before() {
    assert_same_conclusion(
        check_three_term_recurrence_theorem,
        &hermite(),
        0,
        &rows(&[&[1], &[0, 2]]),
        &no_region(),
        1,
    );
    // Chebyshev U with the region (-1, 1), times x^2.
    let chebyshev = rule(&[(1, 0, &[&[0, 2]]), (2, 0, &[&[-1]])], None);
    assert_same_conclusion(
        check_three_term_recurrence_theorem,
        &chebyshev,
        0,
        &rows(&[&[1], &[0, 2]]),
        &open_region(Some(-1), Some(1)),
        2,
    );
    // A left-hand factor x n (f(n) = n is not constant).
    let report = assert_same_conclusion(
        check_three_term_recurrence_theorem,
        &laguerre_reflected(),
        0,
        &rows(&[&[1], &[1, 1]]),
        &open_region(None, Some(0)),
        1,
    );
    let normalization = report.left_hand_factor_normalization.as_ref().unwrap();
    assert!(normalization.statement.contains("f(n) = n"));
    assert!(normalization.normalized_rule.denominator.is_some());
}

#[test]
fn liu_wang_rules_times_x_power_are_proved_as_before() {
    // Type B Eulerian: B_n = (1 + (2n - 1) x) B_(n-1) + 2 x (1 - x) B'_(n-1).
    let type_b = rule(
        &[(1, 0, &[&[1, -1], &[0, 2]]), (1, 1, &[&[0, 2, -2]])],
        None,
    );
    let report = assert_same_conclusion(
        check_liu_wang_theorem,
        &type_b,
        0,
        &rows(&[&[1]]),
        &no_region(),
        1,
    );
    // f = 1 after the division: the normalized rule has no left-hand factor.
    let normalization = report.left_hand_factor_normalization.unwrap();
    assert!(normalization.normalized_rule.denominator.is_none());

    // Eulerian A_n from A_0 = A_1 = 1: the induction start moves one row.
    let eulerian = rule(
        &[(1, 0, &[&[1, -1], &[0, 1]]), (1, 1, &[&[0, 1, -1]])],
        None,
    );
    let report = assert_same_conclusion(
        check_liu_wang_theorem,
        &eulerian,
        0,
        &rows(&[&[1]]),
        &open_region(None, Some(0)),
        2,
    );
    assert!(report.start_delay > 0);
}

#[test]
fn derivative_lag_rules_times_x_power_are_proved_as_before() {
    for id in ["A144438", "A008288", "A271697"] {
        let entry = oeis::by_id(id).unwrap();
        let (base, seeds) = entry.recurrence_parts().unwrap();
        let report = assert_same_conclusion(
            check_liu_wang_derivative_lag_theorem,
            &base,
            entry.recurrence_first_index,
            &seeds,
            &no_region(),
            1,
        );
        assert!(report.lean_route.is_some(), "{id}");
    }
}

#[test]
fn ma_wang_rules_times_x_power_are_proved_as_before() {
    // Fubini F_n = x F_(n-1) + x (1 + x) F'_(n-1).
    let fubini = rule(&[(1, 0, &[&[0, 1]]), (1, 1, &[&[0, 1, 1]])], None);
    let report = assert_same_conclusion(
        check_ma_wang_interval_theorem,
        &fubini,
        0,
        &rows(&[&[1]]),
        &no_region(),
        1,
    );
    assert_eq!(
        report.conclusion.unwrap().root_region.unwrap().to_string(),
        "(-1, 0]"
    );
}

/// `x P_n = (x^2 + n - 1) P_(n-1) - x P'_(n-1)` from `P_0 = 1`: the rows
/// are `x^n` (the terms `(n - 1) x^(n-1) - x (n - 1) x^(n-2)` cancel), so
/// every row divides exactly, but `x` does not divide the coefficient of
/// `P_(n-1)`.
fn row_only_division() -> Recurrence {
    rule(
        &[(1, 0, &[&[-1, 0, 1], &[1]]), (1, 1, &[&[0, -1]])],
        Some(&[&[0, 1]]),
    )
}

#[test]
fn a_division_on_the_rows_only_is_not_applicable_with_the_reason() {
    let rule = row_only_division();
    let seeds = rows(&[&[1]]);
    let generated = rule.generate_rows_rational(&seeds, 0, 5).unwrap();
    assert_eq!(generated[3], rows(&[&[0, 0, 0, 1]])[0]);
    let error = LeftHandFactorNormalization::from_recurrence(&rule).unwrap_err();
    assert!(matches!(error, RecurrenceError::Unsupported(_)));
    assert!(error.to_string().contains("x^1 does not divide"), "{error}");
    let checks: [Check; 4] = [
        check_three_term_recurrence_theorem,
        check_liu_wang_theorem,
        check_liu_wang_derivative_lag_theorem,
        check_ma_wang_interval_theorem,
    ];
    for check in checks {
        let report = check(&rule, 0, &seeds, &no_region()).unwrap();
        assert_eq!(report.status, RealRootednessTheoremStatus::NotApplicable);
        assert!(
            report
                .normal_form
                .contains("happens only on the rows, so the rule is not exactly equivalent"),
            "{}",
            report.normal_form
        );
        assert!(report.hypotheses.is_empty());
        assert!(report.left_hand_factor_normalization.is_none());
        report.verify().unwrap();
    }
    // A left-hand factor with two powers of x is left as it is.
    let mut mixed = rule.clone();
    mixed.denominator = Some(bivar(&[&[0, 1, 1]]));
    assert_eq!(
        LeftHandFactorNormalization::from_recurrence(&mixed).unwrap(),
        None
    );
}

#[test]
fn terms_with_the_same_lag_are_added_before_dividing() {
    // x P_n = (x^2 + 1) P_(n-1) - P_(n-1): the sum x^2 is divisible.
    let split = rule(&[(1, 0, &[&[1, 0, 1]]), (1, 0, &[&[-1]])], Some(&[&[0, 1]]));
    let normalization = LeftHandFactorNormalization::from_recurrence(&split)
        .unwrap()
        .unwrap();
    let normalized = normalization.normalized_recurrence().unwrap();
    assert_eq!(normalized.terms.len(), 1);
    assert_eq!(normalized.terms[0].coeff, bivar(&[&[0, 1]]));
    assert_eq!(normalized.denominator, None);
    normalization.verify().unwrap();
}

#[test]
fn a_veronese_parent_with_an_x_power_left_factor_is_certified() {
    // Q_(n+2) = 2 Q_(n+1) + (y - 1) Q_n, Q_0 = Q_1 = 1: the even sections of
    // (1 + x)^n, with parent P_n = (1 + x) P_(n-1) given as x P_n = (x + x^2) P_(n-1).
    let input = rule(&[(1, 0, &[&[2]]), (2, 0, &[&[-1, 1]])], None);
    let parent = rule(&[(1, 0, &[&[1, 1]])], None);
    let certify = |parent: &Recurrence, input: &Recurrence| {
        let candidate = VeroneseSectionCandidate {
            modulus: 2,
            residue: 0,
            parent: RecurrenceJson::from_recurrence_rational(parent, 0, &[vec![q(1)]], None),
            index_shift: 0,
            gauge: VeroneseSectionGauge::None,
        };
        certify_veronese_section_lift(
            input,
            &rows(&[&[1], &[1]]),
            0,
            &candidate,
            &VeroneseSectionCertifyOptions::default(),
        )
    };
    let certification = certify(&multiplied_by_x_power(&parent, 1), &input).unwrap();
    assert_eq!(certification.verdict, VeroneseSectionVerdict::Proved);
    certification.verify().unwrap();
    // The input may carry the factor too.
    let certification = certify(&parent, &multiplied_by_x_power(&input, 2)).unwrap();
    assert_eq!(certification.verdict, VeroneseSectionVerdict::Proved);
    certification.verify().unwrap();
    // A parent whose division by x happens only on the rows is rejected.
    let error = certify(&row_only_division(), &input).unwrap_err();
    assert!(matches!(error, RecurrenceError::Unsupported(_)), "{error}");
    assert!(error.to_string().contains("parent recurrence"), "{error}");
}

#[test]
fn a_row_factor_lift_accepts_an_x_power_left_factor() {
    // Q_n = (1 + x)^(n-1) A_n (A_n Eulerian, n >= 1), fitted, then times x.
    let eulerian = oeis::by_id("A008292").unwrap();
    let mut power = vec![q(1)];
    let multiply = |a: &[Q], b: &[Q]| {
        let mut out = vec![q(0); a.len() + b.len() - 1];
        for (i, x) in a.iter().enumerate() {
            for (j, y) in b.iter().enumerate() {
                out[i + j] += x * y;
            }
        }
        out
    };
    let family: Vec<Vec<Q>> = eulerian
        .generate_rows(14)
        .unwrap()
        .iter()
        .map(|row| {
            let row: Vec<Q> = row.iter().cloned().map(Q::from_integer).collect();
            let out = multiply(&row, &power);
            power = multiply(&power, &[q(1), q(1)]);
            out
        })
        .collect();
    let fit = polytool::recurrence::find_recurrence_adaptive_rational(
        &family,
        &polytool::recurrence::AdaptiveSearchOptions::default(),
    )
    .unwrap();
    let seeds = family[..fit.recurrence.max_offset()].to_vec();
    let lifted = multiplied_by_x_power(&fit.recurrence, 1);
    let report =
        find_row_factor_lifts(&lifted, &seeds, 1, &RowFactorLiftOptions::default()).unwrap();
    report.verify().unwrap();
    assert!(report.proved_lifts().next().is_some());
}
