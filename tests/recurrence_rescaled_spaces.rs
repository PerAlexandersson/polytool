//! Complete recurrence spaces after factorial coefficient rescaling
//! `T_m(X^k) = X^k/(k!)^m`.
//!
//! The source rows come from `ConstantLagSource::generate_rows`; an
//! independent direct formula `Q[n] = sum_k binom(n,k) X^k/(k!)^m` and a
//! direct residual evaluator check certified identities far beyond the
//! certified prefix.

use polytool::numeric::{binomial_bigint, factorial_bigint};
use polytool::recurrence::*;
use polytool::Polynomial;
use serde_json::{json, Value};

type Poly = Polynomial<BigRational>;

/// Canonical A007318 source `P[n+1] = (1+X) P[n]`, seed `P[origin] = (1+X)^origin`.
fn binomial_source(origin: usize) -> Value {
    let seed: Vec<String> = (0..=origin)
        .map(|k| binomial_bigint(origin, k).to_string())
        .collect();
    json!({"coefficients":[["1","1"]],"seeds":[seed],"origin":origin})
}

fn rescaled_request(
    source: Value,
    m: usize,
    domain_start: usize,
    shapes: Value,
) -> RecurrenceSpaceRequest {
    serde_json::from_value(json!({
        "source": source,
        "coefficient_rescaling": {"kind": "factorial_power", "power": m},
        "domain_start": domain_start,
        "shapes": shapes
    }))
    .unwrap()
}

fn box_shape(length: usize, derivative: usize, index: usize, variable: usize) -> Value {
    json!({"kind":"box","length":length,"derivative_order":derivative,
        "index_degree":index,"variable_degree":variable,"leading_index_degree":index})
}

fn run(req: &RecurrenceSpaceRequest) -> RecurrenceSpaces {
    let report = find_recurrence_spaces(req).unwrap();
    verify_recurrence_spaces(req, &report).unwrap();
    // Serialization round trip, then replay against the embedded request.
    let text = serde_json::to_string(&report).unwrap();
    let back: RecurrenceSpaces = serde_json::from_str(&text).unwrap();
    assert_eq!(back, report);
    verify_recurrence_spaces(&back.request, &back).unwrap();
    report
}

fn q(s: &str) -> BigRational {
    parse_rational_coeff(s).unwrap()
}

/// Independent rows `Q[n] = sum_k binom(n,k) X^k/(k!)^m`.
fn direct_rows(m: u32, count: usize) -> Vec<Poly> {
    (0..count)
        .map(|n| {
            Poly::new(
                (0..=n)
                    .map(|k| BigRational::new(binomial_bigint(n, k), factorial_bigint(k).pow(m)))
                    .collect(),
            )
        })
        .collect()
}

/// Residual `sum_c sign(c) v_c n^u X^v D^s Q[n+lag]` at base index `n`.
fn residual(result: &RecurrenceSpaceResult, vector: &[String], rows: &[Poly], n: usize) -> Poly {
    let mut total = Poly::new(vec![]);
    for (c, v) in result.monomials.iter().zip(vector) {
        let sign = if c.lag == result.target_lag { 1 } else { -1 };
        let factor =
            q(v) * BigRational::from_integer((sign * (n as i64).pow(c.index_degree as u32)).into());
        let term = rows[n + c.lag]
            .nth_derivative(c.derivative_order)
            .scale(&factor)
            * Poly::monomial(BigRational::from_integer(1.into()), c.variable_degree);
        total = total + term;
    }
    total
}

fn assert_kernel_holds_beyond_prefix(result: &RecurrenceSpaceResult, rows: &[Poly], from: usize) {
    let basis = &result.certificate.as_ref().unwrap().basis;
    for vector in basis {
        for n in from..rows.len() - result.target_lag {
            assert!(
                residual(result, vector, rows, n).is_zero(),
                "fails at n = {n}"
            );
        }
    }
}

fn kernel_contains(result: &RecurrenceSpaceResult, expected: &[(RecurrenceMonomial, i64)]) {
    let basis = &result.certificate.as_ref().unwrap().basis;
    assert_eq!(basis.len(), 1, "expected a one-dimensional kernel");
    let pivot = result
        .monomials
        .iter()
        .position(|c| *c == expected[0].0)
        .unwrap();
    let scale = q(&basis[0][pivot]) / BigRational::from_integer(expected[0].1.into());
    for (i, c) in result.monomials.iter().enumerate() {
        let want = expected.iter().find(|(m, _)| m == c).map_or(0, |(_, v)| *v);
        assert_eq!(
            q(&basis[0][i]),
            scale.clone() * BigRational::from_integer(want.into()),
            "column {c:?}"
        );
    }
}

fn mono(
    lag: usize,
    derivative_order: usize,
    index_degree: usize,
    variable_degree: usize,
) -> RecurrenceMonomial {
    RecurrenceMonomial {
        lag,
        derivative_order,
        index_degree,
        variable_degree,
    }
}

#[test]
fn a007318_first_order_derivative_recurrence_is_certified() {
    // (n+1) Q[n+1] = (n+1+X) Q[n] + X Q'[n], n >= 0, with Q[0] = 1.
    let req = rescaled_request(binomial_source(0), 1, 0, json!([box_shape(1, 1, 1, 1)]));
    let report = run(&req);
    let r = &report.results[0];
    assert_eq!(r.status, "successor_exists");
    assert_eq!(r.evidence, "conditional_on_source");
    // U = A = J = 1, m = 1: H = 1 + 1(1+1) + 1 = 4 residuals for order r = 1.
    assert_eq!((r.h, r.residual_count), (4, 4));
    let bound = r.rescaling_bound.as_ref().unwrap();
    assert_eq!(bound.annihilator_power, 4);
    assert_eq!(bound.source_order, 1);
    assert_eq!(bound.initial_rows, vec![vec!["1".to_string()]]);
    assert!(bound.statement.contains("complete space"));
    kernel_contains(
        r,
        &[
            (mono(1, 0, 0, 0), 1),
            (mono(1, 0, 1, 0), 1),
            (mono(0, 0, 0, 0), 1),
            (mono(0, 0, 1, 0), 1),
            (mono(0, 0, 0, 1), 1),
            (mono(0, 1, 0, 1), 1),
        ],
    );
    let rep = &r.representatives[0];
    assert_eq!(
        rep.formula,
        "(1 + n) Q[n+1] = (1 + n + X) Q[n+0] + (X) D^1(Q[n+0])"
    );
    assert_eq!(rep.singular_ranks, Some(vec![]));
    let rows = direct_rows(1, 40);
    assert_kernel_holds_beyond_prefix(r, &rows, 0);
    // The existing backward evaluator regenerates Q from Q[0] = 1.
    let rec = rep.to_recurrence(r).unwrap();
    let exact: Vec<Vec<BigRational>> = rows.iter().map(|p| p.coeffs().to_vec()).collect();
    for target in 1..exact.len() {
        assert_eq!(
            Poly::new(rec.evaluate_next_rational(&exact[..target], 0).unwrap()),
            rows[target]
        );
    }
}

#[test]
fn a007318_three_term_recurrence_is_certified() {
    // (n+1) Q[n+1] = (2n+1+X) Q[n] - n Q[n-1] for n >= 1, written from base
    // index n' = n - 1 >= 0: (n'+2) Q[n'+2] = (2n'+3+X) Q[n'+1] - (n'+1) Q[n'].
    let req = rescaled_request(binomial_source(0), 1, 0, json!([box_shape(2, 0, 1, 1)]));
    let report = run(&req);
    let r = &report.results[0];
    assert_eq!(r.status, "successor_exists");
    assert_eq!(r.residual_count, 3);
    kernel_contains(
        r,
        &[
            (mono(2, 0, 0, 0), 2),
            (mono(2, 0, 1, 0), 1),
            (mono(1, 0, 0, 0), 3),
            (mono(1, 0, 1, 0), 2),
            (mono(1, 0, 0, 1), 1),
            (mono(0, 0, 0, 0), -1),
            (mono(0, 0, 1, 0), -1),
        ],
    );
    assert_eq!(
        r.rescaling_bound.as_ref().unwrap().initial_rows,
        vec![
            vec!["1".to_string()],
            vec!["1".to_string(), "1".to_string()]
        ]
    );
    assert_kernel_holds_beyond_prefix(r, &direct_rows(1, 40), 0);
}

#[test]
fn source_recurrence_input_and_source_coordinate_consistency_rows() {
    let rec = Recurrence {
        terms: vec![RecurrenceTerm::new(
            1,
            0,
            BivarPoly::from_coefficients_normalized(vec![vec![q("1"), q("1")]]),
        )],
        denominator: None,
        inhomogeneous: None,
    };
    let json = RecurrenceJson::from_recurrence_rational(&rec, 0, &[vec![q("1")]], None);
    let req: RecurrenceSpaceRequest = serde_json::from_value(json!({
        "source_recurrence": json,
        "coefficient_rescaling": {"kind": "factorial_power", "power": 1},
        "shapes": [box_shape(1, 1, 1, 1)]
    }))
    .unwrap();
    assert_eq!(run(&req).results[0].status, "successor_exists");
    // Consistency rows are rows of the source P, not of Q.
    let mut source_rows =
        rescaled_request(binomial_source(0), 1, 0, json!([box_shape(1, 1, 1, 1)]));
    source_rows.rows = Some(vec![
        vec!["1".into()],
        vec!["1".into(), "1".into()],
        vec!["1".into(), "2".into(), "1".into()],
    ]);
    run(&source_rows);
    let mut rescaled_rows = source_rows.clone();
    rescaled_rows.rows.as_mut().unwrap()[2] = vec!["1".into(), "2".into(), "1/2".into()];
    assert!(find_recurrence_spaces(&rescaled_rows)
        .unwrap_err()
        .contains("disagree"));
}

#[test]
fn second_power_with_derivatives_on_shifted_rows() {
    let rows = direct_rows(2, 40);
    // Pure-shift order-3 recurrence for m = 2.
    let req = rescaled_request(binomial_source(0), 2, 0, json!([box_shape(3, 0, 2, 1)]));
    let report = run(&req);
    let r = &report.results[0];
    assert_eq!(r.status, "successor_exists");
    // U = 2, A = 1, J = 0, m = 2: H = 2 + 2 + 1 = 5.
    assert_eq!(r.residual_count, 5);
    assert_eq!(
        r.representatives[0].formula,
        "(9 + 6n + n^2) Q[n+3] = (2 + 3n + n^2) Q[n+0] + (-12 - 12n - 3n^2) Q[n+1] \
         + (19 + 15n + 3n^2 + X) Q[n+2]"
    );
    assert_kernel_holds_beyond_prefix(r, &rows, 0);
    // Derivatives on the shifted row Q[n+1]: theta^2 (Q[n+1] - Q[n]) = X Q[n]
    // gives X Q''[n+1] + Q'[n+1] - X Q''[n] - Q'[n] - Q[n] = 0, a zero-leader
    // identity, so the shape has a kernel but no successor.
    let req = rescaled_request(binomial_source(0), 2, 0, json!([box_shape(2, 2, 1, 1)]));
    let report = run(&req);
    let r = &report.results[0];
    assert_eq!(r.status, "absent_in_shape");
    // U = 1, A = 1, J = 2, m = 2: H = 1 + 2(1+2) + 1 = 8.
    assert_eq!(r.residual_count, 8);
    let basis = &r.certificate.as_ref().unwrap().basis;
    assert!(!basis.is_empty());
    assert_kernel_holds_beyond_prefix(r, &rows, 0);
    let identity: Vec<String> = r
        .monomials
        .iter()
        .map(|c| {
            let v = match (c.lag, c.derivative_order, c.index_degree, c.variable_degree) {
                (1, 2, 0, 1) | (1, 1, 0, 0) => -1,
                (0, 2, 0, 1) | (0, 1, 0, 0) | (0, 0, 0, 0) => 1,
                _ => 0,
            };
            v.to_string()
        })
        .collect();
    for n in 0..rows.len() - 2 {
        assert!(residual(r, &identity, &rows, n).is_zero());
    }
    // The identity lies in the certified kernel: it is fixed by its free coordinates.
    let free = &r.certificate.as_ref().unwrap().free_columns;
    let combination =
        basis
            .iter()
            .zip(free)
            .fold(vec![q("0"); r.monomials.len()], |mut acc, (v, &f)| {
                for (a, x) in acc.iter_mut().zip(v) {
                    *a += q(x) * q(&identity[f]);
                }
                acc
            });
    assert_eq!(
        combination,
        identity.iter().map(|s| q(s)).collect::<Vec<_>>()
    );
}

#[test]
fn changed_origin_keeps_original_indices() {
    // P[n] = (1+X)^n from origin 3, residuals from n = 4.
    let mut req = rescaled_request(binomial_source(3), 1, 4, json!([box_shape(1, 1, 1, 1)]));
    req.include_leader_coverage = true;
    let report = run(&req);
    let r = &report.results[0];
    assert_eq!(r.status, "successor_exists");
    assert_eq!(r.last_required_source_row, Some(8));
    let coverage = r.leader_coverage.as_ref().unwrap();
    assert!(coverage.formula.contains("Q[n+1]") && !coverage.formula.contains('P'));
    assert_eq!(
        r.representatives[0].formula,
        "(1 + n) Q[n+1] = (1 + n + X) Q[n+0] + (X) D^1(Q[n+0])"
    );
    assert_eq!(
        r.rescaling_bound.as_ref().unwrap().initial_rows,
        vec![["1", "4", "3", "2/3", "1/24"].map(String::from).to_vec()]
    );
    let mut moved = req.clone();
    moved.domain_start = 3;
    assert!(verify_recurrence_spaces(&moved, &report).is_err());
}

#[test]
fn empty_and_full_kernels() {
    let req = rescaled_request(binomial_source(0), 2, 0, json!([box_shape(1, 2, 1, 1)]));
    let report = run(&req);
    let r = &report.results[0];
    assert_eq!(r.status, "absent_in_shape");
    assert_eq!(r.evidence, "rows_only_obstruction");
    let cert = r.certificate.as_ref().unwrap();
    assert!(cert.basis.is_empty());
    assert_eq!(cert.absence_duals.len(), 2);
    // A zero sequence: every identity holds, so the kernel is everything.
    let zero = json!({"coefficients":[["1","1"]],"seeds":[["0"]]});
    let req = rescaled_request(zero, 1, 0, json!([box_shape(1, 1, 1, 1)]));
    let report = run(&req);
    let r = &report.results[0];
    assert!(r.all_checked_rows_zero);
    assert_eq!(
        r.certificate.as_ref().unwrap().basis.len(),
        r.monomials.len()
    );
    assert_eq!(r.status, "successor_exists");
}

#[test]
fn theorem_hypotheses_are_validated() {
    let shapes = json!([box_shape(1, 1, 1, 1)]);
    let zero_power = rescaled_request(binomial_source(0), 0, 0, shapes.clone());
    assert!(find_recurrence_spaces(&zero_power)
        .unwrap_err()
        .contains("m >= 1"));
    let mut huge = rescaled_request(binomial_source(0), 600, 0, shapes.clone());
    assert!(find_recurrence_spaces(&huge).is_err());
    huge.limits.max_degree = 1000;
    huge.coefficient_rescaling =
        Some(RecurrenceCoefficientRescaling::FactorialPower { power: 600 });
    let report = find_recurrence_spaces(&huge).unwrap();
    assert_eq!(report.results[0].status, "incomplete_budget");
    let rows_only: Result<RecurrenceSpaceRequest, _> = serde_json::from_value(json!({
        "rows": [["1"], ["1", "1"], ["1", "2", "1"]],
        "coefficient_rescaling": {"kind": "factorial_power", "power": 1},
        "shapes": shapes
    }));
    assert!(find_recurrence_spaces(&rows_only.unwrap())
        .unwrap_err()
        .contains("constant-index monic source"));
    let indexed: RecurrenceSpaceRequest = serde_json::from_value(json!({
        "indexed_source": {"coefficients": [[["0","1"],["0","1"]]], "seeds": [["1"]]},
        "coefficient_rescaling": {"kind": "factorial_power", "power": 1},
        "shapes": shapes
    }))
    .unwrap();
    assert!(find_recurrence_spaces(&indexed).is_err());
    let automatic = rescaled_request(
        binomial_source(0),
        1,
        0,
        json!([{"kind":"automatic_variable_degree","length":1,"derivative_order":1,"index_degree":1}]),
    );
    assert!(find_recurrence_spaces(&automatic)
        .unwrap_err()
        .contains("automatic_variable_degree"));
    let before_origin = rescaled_request(binomial_source(2), 1, 1, shapes.clone());
    assert!(find_recurrence_spaces(&before_origin).is_err());
    for bad in [
        json!({"kind": "factorial_power"}),
        json!({"kind": "factorial_power", "power": 1, "shift": 1}),
        json!({"kind": "falling_factorial", "power": 1}),
    ] {
        let parsed: Result<RecurrenceSpaceRequest, _> = serde_json::from_value(json!({
            "source": binomial_source(0), "coefficient_rescaling": bad, "shapes": shapes
        }));
        assert!(parsed.is_err());
    }
}

#[test]
fn tampered_source_transform_support_and_certificates_are_rejected() {
    let req = rescaled_request(binomial_source(0), 1, 0, json!([box_shape(1, 1, 1, 1)]));
    let report = run(&req);
    // Changed transform in the expected request.
    let mut other_power = req.clone();
    other_power.coefficient_rescaling =
        Some(RecurrenceCoefficientRescaling::FactorialPower { power: 2 });
    assert!(verify_recurrence_spaces(&other_power, &report).is_err());
    // A report whose embedded transform was altered (or removed) no longer replays.
    let mut altered = report.clone();
    altered.request.coefficient_rescaling =
        Some(RecurrenceCoefficientRescaling::FactorialPower { power: 2 });
    assert!(verify_recurrence_spaces(&altered.request, &altered).is_err());
    let mut removed = report.clone();
    removed.request.coefficient_rescaling = None;
    assert!(verify_recurrence_spaces(&removed.request, &removed).is_err());
    // Changed source seeds.
    let mut other_seed = req.clone();
    other_seed.source.as_mut().unwrap().seeds = vec![vec!["2".into()]];
    assert!(verify_recurrence_spaces(&other_seed, &report).is_err());
    let mut seeded = report.clone();
    seeded.request = other_seed.clone();
    assert!(verify_recurrence_spaces(&other_seed, &seeded).is_err());
    // Changed support.
    let mut support = report.clone();
    support.results[0].monomials[0].variable_degree = 2;
    assert!(verify_recurrence_spaces(&req, &support).is_err());
    // Changed theorem data.
    let mut bound = report.clone();
    bound.results[0]
        .rescaling_bound
        .as_mut()
        .unwrap()
        .annihilator_power = 2;
    assert!(verify_recurrence_spaces(&req, &bound).is_err());
    let mut initial = report.clone();
    initial.results[0]
        .rescaling_bound
        .as_mut()
        .unwrap()
        .initial_rows = vec![vec!["2".into()]];
    assert!(verify_recurrence_spaces(&req, &initial).is_err());
    let mut count = report.clone();
    count.results[0].residual_count = 2;
    assert!(verify_recurrence_spaces(&req, &count).is_err());
    // Changed kernel vector.
    let mut basis = report.clone();
    basis.results[0].certificate.as_mut().unwrap().basis[0][0] = "2".into();
    assert!(verify_recurrence_spaces(&req, &basis).is_err());
}

#[test]
fn budget_outcomes_are_explicit() {
    let shapes = json!([box_shape(1, 1, 1, 1)]);
    let mut rows = rescaled_request(binomial_source(0), 1, 0, shapes.clone());
    rows.limits.max_rows = 3;
    let report = run(&rows);
    assert_eq!(report.results[0].status, "incomplete_budget");
    assert_eq!(report.results[0].evidence, "none");
    assert!(report.results[0]
        .reason
        .as_ref()
        .unwrap()
        .contains("budget"));
    let mut bits = rescaled_request(binomial_source(0), 3, 0, shapes);
    bits.limits.max_coefficient_bits = 8;
    let report = run(&bits);
    assert_eq!(report.results[0].status, "incomplete_budget");
    assert!(report.results[0]
        .reason
        .as_ref()
        .unwrap()
        .contains("rescaling coefficient bit budget"));
}
