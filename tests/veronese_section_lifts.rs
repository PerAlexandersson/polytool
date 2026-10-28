//! End-to-end tests of `recurrence::find_veronese_section_lifts` and the
//! `polytool recurrence-veronese-lifts` command.
//!
//! The searches factor `chi(L, x^r)` with `polytool::factorization`; tests
//! that need a successful factorization are marked
//! is unit-tested in the library.

use num_bigint::BigInt;
use polytool::recurrence::{
    find_veronese_section_lifts, format_characteristic_polynomial, BigRational, BivarPoly,
    Recurrence, RecurrenceError, RecurrenceJson, RecurrenceSign, RecurrenceTerm,
    VeroneseSectionLift, VeroneseSectionLiftOptions, VeroneseSectionLiftRejectionReason,
    VeroneseSectionLiftReport, VeroneseSectionLiftStatus,
};
use std::io::Write;
use std::process::{Command, Output, Stdio};

fn q(value: i64) -> BigRational {
    BigRational::from_integer(BigInt::from(value))
}

fn row(values: &[i64]) -> Vec<BigRational> {
    values.iter().map(|&value| q(value)).collect()
}

fn rows(values: &[&[i64]]) -> Vec<Vec<BigRational>> {
    values.iter().map(|values| row(values)).collect()
}

/// `Q_n = sum_k c_k(y) Q_(n-k)` from `(k, ascending coefficients of c_k)`.
fn recurrence(terms: &[(usize, &[i64])]) -> Recurrence {
    Recurrence {
        terms: terms
            .iter()
            .map(|(offset, coefficients)| {
                RecurrenceTerm::new(
                    *offset,
                    0,
                    BivarPoly::from_coefficients_normalized(vec![row(coefficients)]),
                )
            })
            .collect(),
        denominator: None,
        inhomogeneous: None,
    }
}

/// `Q_(n+2) = 2 Q_(n+1) + (y - 1) Q_n`.
fn worked_recurrence() -> Recurrence {
    recurrence(&[(1, &[2]), (2, &[-1, 1])])
}

fn lifts(
    rule: &Recurrence,
    seeds: &[&[i64]],
    options: &VeroneseSectionLiftOptions,
) -> VeroneseSectionLiftReport {
    find_veronese_section_lifts(rule, &rows(seeds), 0, options).expect("search runs")
}

fn parent_text(lift: &VeroneseSectionLift) -> String {
    format_characteristic_polynomial(&lift.parent_characteristic_polynomial, "x")
}

fn find<'a>(
    report: &'a VeroneseSectionLiftReport,
    modulus: usize,
    residue: usize,
    parent: &str,
) -> Option<&'a VeroneseSectionLift> {
    report.lifts.iter().find(|lift| {
        lift.modulus == modulus && lift.residue == residue && parent_text(lift) == parent
    })
}

fn binomial(n: u64, k: u64) -> i64 {
    if k > n {
        return 0;
    }
    (0..k).fold(1i64, |acc, i| acc * (n - i) as i64 / (i + 1) as i64)
}

/// `V_(r, j)((a + b x)^n)` for `n` in `0..count`, as integer rows.
fn product_sections(a: i64, b: i64, modulus: u64, residue: u64, count: u64) -> Vec<Vec<i64>> {
    (0..count)
        .map(|n| {
            let mut coefficients = Vec::new();
            let mut power = residue;
            while power <= n {
                coefficients
                    .push(binomial(n, power) * a.pow((n - power) as u32) * b.pow(power as u32));
                power += modulus;
            }
            if coefficients.is_empty() {
                coefficients.push(0);
            }
            coefficients
        })
        .collect()
}

fn as_slices(rows: &[Vec<i64>]) -> Vec<&[i64]> {
    rows.iter().map(Vec::as_slice).collect()
}

#[test]
fn worked_example_is_the_even_section_of_one_plus_x() {
    let report = lifts(&worked_recurrence(), &[&[1], &[1]], &Default::default());
    let best = &report.lifts[0];
    assert_eq!((best.modulus, best.residue), (2, 0));
    assert_eq!(parent_text(best), "L - 1 - x");
    assert_eq!(best.parent_initial_rows, rows(&[&[1]]));
    assert_eq!(best.status, VeroneseSectionLiftStatus::ProvedProductParent);
    let conjugate = find(&report, 2, 0, "L - 1 + x").expect("conjugate parent listed");
    assert!(conjugate.status.is_failing());
    assert_eq!(
        report.lifts.last().unwrap().status.label(),
        "parent_row_mixed_signs"
    );
    assert!(report.skipped_moduli.is_empty());
    assert!(!report.exhausted_budget);
}

#[test]
fn odd_section_of_one_plus_x() {
    let report = lifts(&worked_recurrence(), &[&[0], &[1]], &Default::default());
    let best = &report.lifts[0];
    assert_eq!((best.modulus, best.residue), (2, 1));
    assert_eq!(parent_text(best), "L - 1 - x");
    assert_eq!(best.parent_initial_rows, rows(&[&[1]]));
    assert!(best.status.is_proved());
}

#[test]
fn modulus_three_sections_of_one_plus_x() {
    // Q_(n+3) = 3 Q_(n+2) - 3 Q_(n+1) + (1 + y) Q_n.
    let rule = recurrence(&[(1, &[3]), (2, &[-3]), (3, &[1, 1])]);
    for residue in 0..3u64 {
        let seeds = product_sections(1, 1, 3, residue, 3);
        let report = lifts(&rule, &as_slices(&seeds), &Default::default());
        let lift = find(&report, 3, residue as usize, "L - 1 - x").expect("modulus-3 lift");
        assert_eq!(lift.parent_initial_rows, rows(&[&[1]]));
        assert!(lift.status.is_proved());
        assert!(report.lifts[0].status.is_proved());
    }
    let report = lifts(&rule, &[&[1], &[1], &[1]], &Default::default());
    assert_eq!((report.lifts[0].modulus, report.lifts[0].residue), (3, 0));
}

#[test]
fn sections_of_two_plus_three_x_and_modulus_four() {
    // (L - 2)^2 - 9y: Q_(n+2) = 4 Q_(n+1) + (9y - 4) Q_n.
    let rule = recurrence(&[(1, &[4]), (2, &[-4, 9])]);
    for residue in 0..2u64 {
        let seeds = product_sections(2, 3, 2, residue, 2);
        let report = lifts(&rule, &as_slices(&seeds), &Default::default());
        let best = &report.lifts[0];
        assert_eq!((best.modulus, best.residue), (2, residue as usize));
        assert_eq!(parent_text(best), "L - 2 - 3x");
        assert_eq!(best.parent_initial_rows, rows(&[&[1]]));
        assert!(best.status.is_proved());
    }
    // (L - 1)^4 - y: sections of (1 + x)^n modulo 4.
    let rule = recurrence(&[(1, &[4]), (2, &[-6]), (3, &[4]), (4, &[-1, 1])]);
    let seeds = product_sections(1, 1, 4, 0, 4);
    let report = lifts(&rule, &as_slices(&seeds), &Default::default());
    let best = &report.lifts[0];
    assert_eq!((best.modulus, best.residue), (4, 0));
    assert_eq!(parent_text(best), "L - 1 - x");
    assert!(best.status.is_proved());
    // The order-two parent sum_k C(n, 2k) x^k at modulus 2 stays open.
    let open = find(&report, 2, 0, "L^2 - 2L + 1 - x").expect("order-two parent");
    assert_eq!(
        open.status,
        VeroneseSectionLiftStatus::ParentRealRootednessOpen
    );
}

#[test]
fn even_sections_of_fibonacci_polynomials() {
    // chi = L^4 - 2L^3 + L^2 - y.
    let rule = recurrence(&[(1, &[2]), (2, &[-1]), (4, &[0, 1])]);
    let report = lifts(&rule, &[&[1], &[1], &[1], &[1]], &Default::default());
    assert_eq!(
        format_characteristic_polynomial(&report.input_characteristic_polynomial, "y"),
        "L^4 - 2L^3 + L^2 - y"
    );
    let best = &report.lifts[0];
    assert_eq!((best.modulus, best.residue), (2, 0));
    assert_eq!(parent_text(best), "L^2 - L - x");
    assert_eq!(best.parent_initial_rows, rows(&[&[1], &[1]]));
    assert_eq!(
        best.status,
        VeroneseSectionLiftStatus::ParentRealRootednessOpen
    );
    assert_eq!(best.parent_row_check.rows_checked, 12);
    assert_eq!(best.parent_row_check.first_not_real_rooted_row_index, None);
}

#[test]
fn non_minimal_input_recurrence() {
    // chi = (L - 2)((L - 1)^2 - y).
    let rule = recurrence(&[(1, &[4]), (2, &[-5, 1]), (3, &[2, -2])]);
    let report = lifts(&rule, &[&[1], &[1], &[1, 1]], &Default::default());
    let best = &report.lifts[0];
    assert_eq!((best.modulus, best.residue), (2, 0));
    assert_eq!(parent_text(best), "L - 1 - x");
    assert_eq!(best.cofactor, rows(&[&[-2], &[1]]));
    assert!(best.status.is_proved());
    assert!(report.rejected_candidates.iter().any(|rejection| {
        rejection.reason
            == VeroneseSectionLiftRejectionReason::NormDoesNotDivideCharacteristicPolynomial
    }));
}

#[test]
fn near_misses_are_not_accepted() {
    // Changed seeds: residue 0 would need the seed section 1/y.
    let report = lifts(&worked_recurrence(), &[&[1], &[2]], &Default::default());
    assert!(find(&report, 2, 0, "L - 1 - x").is_none());
    assert!(report
        .rejected_candidates
        .iter()
        .any(|rejection| rejection.residue == Some(0)
            && rejection.reason == VeroneseSectionLiftRejectionReason::SeedsNotPolynomial));
    // (These rows are the odd sections of (1 + x)^(n+1), a different parent instance.)
    let shifted = find(&report, 2, 1, "L - 1 - x").expect("shifted parent");
    assert_eq!(shifted.parent_initial_rows, rows(&[&[1, 1]]));

    // Changed coefficient: chi(L, x^2) = (L - 1)^2 + 1 - x^2 is irreducible.
    let changed = recurrence(&[(1, &[2]), (2, &[-2, 1])]);
    let report = lifts(&changed, &[&[1], &[1]], &Default::default());
    assert!(report.lifts.is_empty());

    // Exhausted candidate budget.
    let options = VeroneseSectionLiftOptions {
        max_parent_candidates: 0,
        ..Default::default()
    };
    let report = lifts(&worked_recurrence(), &[&[1], &[1]], &options);
    assert!(report.exhausted_budget);
    assert!(report.lifts.is_empty());
}

#[test]
fn unsupported_recurrences_are_rejected_before_factoring() {
    let seeds = rows(&[&[1], &[1]]);
    let unsupported = |rule: Recurrence| {
        matches!(
            find_veronese_section_lifts(&rule, &seeds, 0, &Default::default()),
            Err(RecurrenceError::Unsupported(_))
        )
    };
    let mut rule = worked_recurrence();
    rule.terms[1].coeff =
        BivarPoly::from_coefficients_normalized(vec![row(&[-1, 1]), row(&[0, 1])]);
    assert!(unsupported(rule), "n-dependent coefficient");
    let mut rule = worked_recurrence();
    rule.terms[1].deriv_order = 1;
    assert!(unsupported(rule), "derivative term");
    let mut rule = worked_recurrence();
    rule.inhomogeneous = Some(BivarPoly::from_coefficients_normalized(vec![row(&[1])]));
    assert!(unsupported(rule), "inhomogeneous term");
    let mut rule = worked_recurrence();
    rule.terms[0] = rule.terms[0]
        .clone()
        .with_sign(RecurrenceSign::AlternatingN);
    assert!(unsupported(rule), "alternating sign");
    assert!(matches!(
        find_veronese_section_lifts(&worked_recurrence(), &rows(&[&[1]]), 0, &Default::default()),
        Err(RecurrenceError::InvalidInput(_))
    ));
}

#[test]
fn parent_recurrence_round_trips_through_recurrence_json() {
    let rule = recurrence(&[(1, &[2]), (2, &[-1]), (4, &[0, 1])]);
    let report = lifts(&rule, &[&[1], &[1], &[1], &[1]], &Default::default());
    let lift = &report.lifts[0];
    let encoded = serde_json::to_string(&lift.parent_recurrence_json()).unwrap();
    let decoded: RecurrenceJson = serde_json::from_str(&encoded).unwrap();
    let (parent, first_index, seeds) = decoded.to_recurrence_parts().unwrap();
    let replayed = parent
        .generate_rows_rational(&seeds, first_index, 10)
        .unwrap();
    let input = rule
        .generate_rows_rational(&rows(&[&[1], &[1], &[1], &[1]]), 0, 10)
        .unwrap();
    for (parent_row, input_row) in replayed.iter().zip(&input) {
        let section: Vec<BigRational> = parent_row.iter().step_by(2).cloned().collect();
        let trim = |row: &[BigRational]| {
            let mut row = row.to_vec();
            while row.last().is_some_and(|c| *c == q(0)) {
                row.pop();
            }
            row
        };
        assert_eq!(trim(&section), trim(input_row));
    }
    // The report itself serializes with exact strings.
    let json = serde_json::to_value(&report).unwrap();
    assert_eq!(json["schema"], "polytool.recurrence-veronese-lifts.v1");
    assert_eq!(
        json["lifts"][0]["parent_recurrence"],
        serde_json::to_value(lift.parent_recurrence_json()).unwrap()
    );
}

const WORKED_JSON: &str = r#"{"schema":"polytool.recurrence.v1","first_index":0,
"initial_polynomials":[["1"],["1"]],"recurrence":{"terms":[
{"offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["2"]]}},
{"offset":2,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["-1","1"]]}}]}}"#;

fn run(args: &[&str], stdin: &[u8]) -> Output {
    let mut process = Command::new(env!("CARGO_BIN_EXE_polytool"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    process.stdin.take().unwrap().write_all(stdin).unwrap();
    process.wait_with_output().unwrap()
}

#[test]
fn cli_help_and_input_errors() {
    let help = run(&["recurrence-veronese-lifts", "--help"], b"");
    assert!(help.status.success());
    let text = String::from_utf8(help.stdout).unwrap();
    assert!(text.contains("V_(r,j)(P)(y)"));
    assert!(text.contains("proved_product_parent"));
    let top = run(&["--help"], b"");
    assert!(String::from_utf8(top.stdout)
        .unwrap()
        .contains("recurrence-veronese-lifts"));

    let alternating = WORKED_JSON.replacen("\"none\"", "\"alternating_n\"", 1);
    let output = run(&["recurrence-veronese-lifts"], alternating.as_bytes());
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8(output.stderr)
        .unwrap()
        .contains("alternating"));
    let output = run(
        &["recurrence-veronese-lifts", "--max-modulus", "7"],
        WORKED_JSON.as_bytes(),
    );
    assert_eq!(output.status.code(), Some(2));
    let output = run(
        &["recurrence-veronese-lifts", "--bogus"],
        WORKED_JSON.as_bytes(),
    );
    assert_eq!(output.status.code(), Some(2));
}

#[test]
fn cli_finds_the_worked_example() {
    let output = run(
        &["recurrence-veronese-lifts", "--recurrence", "-", "--json"],
        WORKED_JSON.as_bytes(),
    );
    assert!(output.status.success());
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["schema"], "polytool.recurrence-veronese-lifts.v1");
    assert_eq!(report["lifts"][0]["status"], "proved_product_parent");
    assert_eq!(
        report["lifts"][0]["parent_characteristic_polynomial_text"],
        "L - 1 - x"
    );
    assert_eq!(
        report["lifts"][0]["parent_recurrence"]["initial_polynomials"],
        serde_json::json!([["1"]])
    );
    let text = run(
        &["recurrence-veronese-lifts", "--max-modulus", "2"],
        WORKED_JSON.as_bytes(),
    );
    assert!(text.status.success());
    let text = String::from_utf8(text.stdout).unwrap();
    assert!(text.contains("Lift 1: modulus 2, residue 0 [proved_product_parent]"));
    assert!(text.contains("parent G(L, x) = L - 1 - x"));
}

#[test]
fn binomial_rows_of_one_parity_have_a_proved_period_two_parent() {
    // A109446: C(n, k) with n - k even.  Its residue alternates with n for the
    // parent (1 + x)^n, so the fixed-residue parent has P_(n+2) = (1 + x)^2 P_n
    // with seeds P_0 = x, P_1 = 1 + x.
    let definition = polytool::oeis::catalog()
        .iter()
        .find(|definition| definition.id == "A109446")
        .expect("A109446 is in the catalog");
    let (rule, initial_rows, first_index) = definition.recurrence_export_parts().unwrap();
    let report = find_veronese_section_lifts(
        &rule,
        &initial_rows,
        usize::try_from(first_index).unwrap(),
        &VeroneseSectionLiftOptions::default(),
    )
    .unwrap();
    let lift = find(&report, 2, 1, "L^2 - 1 - 2x - x^2").expect("period-two parent");
    assert_eq!(lift.status, VeroneseSectionLiftStatus::ProvedProductParent);
    assert_eq!(lift.parent_initial_rows, rows(&[&[0, 1], &[1, 1]]));
}
