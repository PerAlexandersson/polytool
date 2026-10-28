//! Row-factor lifts `Q_n = F^(e(n)) S_n`: a `(1 + x)^n` factor on a
//! differential family, a fixed quadratic factor, near misses and
//! certificate replay.

use polytool::oeis;
use polytool::recurrence::*;

fn q(value: i64) -> BigRational {
    BigRational::from_integer(value.into())
}

fn multiply(a: &[BigRational], b: &[BigRational]) -> Vec<BigRational> {
    let mut out = vec![q(0); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    out
}

fn certificate(lift: &RowFactorLift) -> &RowFactorCertificate {
    match &lift.identity {
        RowFactorIdentityStatus::Certified { certificate } => certificate,
        other => panic!("not certified: {other:?}"),
    }
}

/// `Q_n = (1 + x)^(n-1) A_n` with `A_n` the Eulerian polynomials (A008292,
/// `n >= 1`), given by a fitted first-order differential rule.
fn eulerian_times_power() -> (Recurrence, Vec<Vec<BigRational>>) {
    let eulerian = oeis::by_id("A008292").unwrap();
    let mut power = vec![q(1)];
    let rows: Vec<Vec<BigRational>> = eulerian
        .generate_rows(14)
        .unwrap()
        .iter()
        .map(|row| {
            let row: Vec<BigRational> =
                row.iter().cloned().map(BigRational::from_integer).collect();
            let out = multiply(&row, &power);
            power = multiply(&power, &[q(1), q(1)]);
            out
        })
        .collect();
    // The fitted rule indexes the first row as 1, which is the first index.
    let fit = find_recurrence_adaptive_rational(&rows, &AdaptiveSearchOptions::default()).unwrap();
    assert!(fit.recurrence.terms.iter().any(|t| t.deriv_order > 0));
    let seeds = rows[..fit.recurrence.max_offset()].to_vec();
    (fit.recurrence, seeds)
}

#[test]
fn eulerian_times_power_of_one_plus_x_is_proved() {
    let (rule, seeds) = eulerian_times_power();
    let report = find_row_factor_lifts(&rule, &seeds, 1, &RowFactorLiftOptions::default()).unwrap();
    report.verify().unwrap();
    let lift = report.proved_lifts().next().expect("a proved lift");
    assert_eq!(lift.factor, vec!["1".to_string(), "1".to_string()]);
    assert_eq!((lift.exponent_slope, lift.exponent_offset), (1, -1));
    assert!(lift.factor_roots_nonpositive);
    let certificate = certificate(lift);
    assert_eq!(certificate.removed_factor_power, 2);
    certificate.verify().unwrap();
    assert_eq!(
        lift.lift_summary.operation,
        LiftOperation::RowFactor {
            factor: vec!["1".into(), "1".into()],
            exponent_slope: 1,
            exponent_offset: -1,
        }
    );

    let json = serde_json::to_value(&report).unwrap();
    assert_eq!(json["schema"], ROW_FACTOR_LIFTS_SCHEMA);
    assert_eq!(
        json["lifts"][0]["lift_summary"]["operation"]["kind"],
        "row_factor"
    );
    let back: RowFactorLiftReport = serde_json::from_value(json).unwrap();
    assert_eq!(back, report);
    back.verify().unwrap();

    // Tampering is caught by the replay.
    let mut wrong_power = certificate.clone();
    wrong_power.removed_factor_power += 1;
    assert!(wrong_power.verify().is_err());
    let mut wrong_factor = certificate.clone();
    wrong_factor.factor = vec!["1".into(), "2".into()];
    assert!(wrong_factor.verify().is_err());
    let mut wrong_exponent = certificate.clone();
    wrong_exponent.exponent_offset = 0;
    assert!(wrong_exponent.verify().is_err());
}

#[test]
fn a_fixed_quadratic_factor_is_found() {
    // Q_n = (1 + 3x + x^2) (1 + x)^n.
    let rule = Recurrence {
        terms: vec![RecurrenceTerm::new(
            1,
            0,
            BivarPoly {
                coeffs: vec![vec![q(1), q(1)]],
            },
        )],
        denominator: None,
        inhomogeneous: None,
    };
    let report = find_row_factor_lifts(
        &rule,
        &[vec![q(1), q(3), q(1)]],
        0,
        &RowFactorLiftOptions::default(),
    )
    .unwrap();
    report.verify().unwrap();
    let quadratic = report
        .lifts
        .iter()
        .find(|lift| lift.factor == ["1", "3", "1"])
        .expect("the quadratic factor");
    assert_eq!(
        (quadratic.exponent_slope, quadratic.exponent_offset),
        (0, 1)
    );
    assert!(quadratic.lift_summary.is_proved());
    let linear = report
        .lifts
        .iter()
        .find(|lift| lift.factor == ["1", "1"])
        .expect("the linear factor");
    assert_eq!((linear.exponent_slope, linear.exponent_offset), (1, 0));
    assert!(linear.lift_summary.is_proved());
}

#[test]
fn a_non_real_rooted_factor_fails_the_route() {
    // Q_n = (1 + x + x^2) (1 + x)^n.
    let rule = Recurrence {
        terms: vec![RecurrenceTerm::new(
            1,
            0,
            BivarPoly {
                coeffs: vec![vec![q(1), q(1)]],
            },
        )],
        denominator: None,
        inhomogeneous: None,
    };
    let report = find_row_factor_lifts(
        &rule,
        &[vec![q(1), q(1), q(1)]],
        0,
        &RowFactorLiftOptions::default(),
    )
    .unwrap();
    report.verify().unwrap();
    let lift = report
        .lifts
        .iter()
        .find(|lift| lift.factor == ["1", "1", "1"])
        .expect("the quadratic factor");
    certificate(lift).verify().unwrap();
    assert!(!lift.factor_real_rooted);
    assert_eq!(lift.lift_summary.verdict, LiftVerdict::RouteFailed);
    assert_eq!(report.proved_lifts().count(), 0);
}

#[test]
fn a_factor_of_finitely_many_rows_is_rejected() {
    // Q_n = (1 + x)^2 + (n - 1)(n - 2)(n - 3) x for n >= 1: (1 + x)^2 divides
    // the first three rows only.
    let rule = Recurrence {
        terms: [(1, 4), (2, -6), (3, 4), (4, -1)]
            .into_iter()
            .map(|(offset, c)| {
                RecurrenceTerm::new(
                    offset,
                    0,
                    BivarPoly {
                        coeffs: vec![vec![q(c)]],
                    },
                )
            })
            .collect(),
        denominator: None,
        inhomogeneous: None,
    };
    let seeds: Vec<Vec<BigRational>> = (1..=4i64)
        .map(|n| vec![q(1), q(2 + (n - 1) * (n - 2) * (n - 3)), q(1)])
        .collect();
    let report = find_row_factor_lifts(
        &rule,
        &seeds,
        1,
        &RowFactorLiftOptions {
            detection_rows: 3,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(report.lifts.is_empty());
    let rejection = &report.rejected[0];
    assert_eq!(rejection.factor, ["1", "1"]);
    assert_eq!(rejection.failing_index, Some(4));
}
