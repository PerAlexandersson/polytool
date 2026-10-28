//! Boundary and error-path tests: the documented input limits, option
//! validation and the invalid-input branches of the tool handlers.
//!
//! The oracles are the published limit constants (a value at the limit is
//! accepted and one more is rejected), known mathematical values (unit-square
//! Ehrhart data, Eulerian and Chebyshev rows) and the exact option mapping.

// `rmcp::Json` has no `Debug`, so `expect_err` is unavailable for tool results.
#![allow(clippy::err_expect)]

use super::*;
use crate::{
    family::family_source_and_polynomials,
    limits::{
        MCP_MAX_BATCH_POLYNOMIALS, MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL, MCP_MAX_RECURRENCE_DEGREE,
        MCP_MAX_RECURRENCE_INPUT_ROWS, MCP_MAX_RECURRENCE_LENGTH, MCP_MAX_SEQUENCE_N,
        MCP_MAX_TEXT_BYTES, MCP_MAX_TOTAL_COEFFICIENTS,
    },
    parsing::{
        parse_bigint_batch, parse_bigint_matrix, parse_recurrence_batch_rational,
        parse_required_bigint_polynomial,
    },
    recurrence_options::{apply_recurrence_options, validate_recurrence_for_generation},
};
use polytool::recurrence::{BigRational, BivarPoly, Recurrence, RecurrenceTerm};
use rmcp::{handler::server::wrapper::Parameters, ErrorData as McpError, Json};
use serde_json::json;

/// One polynomial input parsed exactly, with the error as text.
fn parse_one(input: &BigIntPolynomialInput) -> Result<Vec<num_bigint::BigInt>, String> {
    parse_required_bigint_polynomial(input, "p").map_err(message)
}

fn message(error: McpError) -> String {
    error.message.to_string()
}

fn coeffs(values: &[i64]) -> BigIntPolynomialInput {
    BigIntPolynomialInput {
        coefficients: Some(
            values
                .iter()
                .copied()
                .map(BigIntCoefficientInput::from)
                .collect(),
        ),
        expression: None,
    }
}

fn strings(values: &[i64]) -> Vec<String> {
    values.iter().map(ToString::to_string).collect()
}

fn bigint_coeffs(values: &[&str]) -> BigIntPolynomialInput {
    BigIntPolynomialInput {
        coefficients: Some(
            values
                .iter()
                .map(|value| BigIntCoefficientInput::Text((*value).to_string()))
                .collect(),
        ),
        expression: None,
    }
}

fn family_request(polynomials: &[&[i64]]) -> CheckPolynomialFamilyRequest {
    CheckPolynomialFamilyRequest {
        polynomials: Some(polynomials.iter().map(|values| coeffs(values)).collect()),
        text: None,
        sequence: None,
        max_n: None,
        options: None,
        recurrence_options: None,
        lace: None,
    }
}

fn family_options() -> FamilyCheckOptions {
    FamilyCheckOptions {
        require_real_rooted: None,
        require_simple_roots: None,
        require_palindromic: None,
        require_gamma_positive: None,
        require_unimodal: None,
        require_log_concave: None,
        require_ultra_log_concave: None,
        check_consecutive_interlacing: None,
        require_consecutive_weak_interlacing: None,
        find_recurrence: None,
    }
}

fn first_failure(
    polynomials: &[&[i64]],
    configure: impl FnOnce(&mut FamilyCheckOptions),
) -> Option<String> {
    let mut options = family_options();
    configure(&mut options);
    let mut request = family_request(polynomials);
    request.options = Some(options);
    let Json(response) = PolynomialToolsServer::new()
        .check_polynomial_family(Parameters(request))
        .unwrap();
    assert_eq!(
        response.all_required_checks_passed,
        response.first_failure.is_none()
    );
    response.first_failure
}

// -- check_polynomial_family ------------------------------------------------

#[test]
fn family_sequence_source_matches_classical_rows() {
    let server = PolynomialToolsServer::new();
    let mut request = family_request(&[]);
    request.polynomials = None;
    request.sequence = Some(SequenceKind::Eulerian);
    request.max_n = Some(4);
    let Json(response) = server.check_polynomial_family(Parameters(request)).unwrap();
    // Eulerian polynomials A_1..A_4 (OEIS A008292); the real-rootedness and
    // interlacing of consecutive rows is a classical theorem.
    assert_eq!(response.source, "sequence:Eulerian, max_n=4");
    assert_eq!(response.item_count, 4);
    let rows: Vec<Vec<i64>> = response
        .items
        .iter()
        .map(|item| {
            item.polynomial
                .coefficients
                .iter()
                .map(|c| c.parse().unwrap())
                .collect()
        })
        .collect();
    assert_eq!(
        rows,
        vec![vec![1], vec![1, 1], vec![1, 4, 1], vec![1, 11, 11, 1]]
    );
    assert!(response.all_required_checks_passed);
    assert_eq!(response.consecutive_pairs.len(), 3);
    assert!(response.consecutive_pairs.iter().all(|pair| pair
        .result
        .as_ref()
        .and_then(|result| result.weak)
        == Some(true)));
    assert!(response.items.iter().all(|item| item.gamma_positive));
    assert!(response.markdown.contains("Consecutive interlacing"));

    // The other sequence kinds follow their defining recurrences.
    let rows_of = |sequence: SequenceKind, max_n: usize| -> Vec<Vec<i64>> {
        let mut request = family_request(&[]);
        request.polynomials = None;
        request.sequence = Some(sequence);
        request.max_n = Some(max_n);
        let Json(response) = server.check_polynomial_family(Parameters(request)).unwrap();
        response
            .items
            .iter()
            .map(|item| {
                item.polynomial
                    .coefficients
                    .iter()
                    .map(|c| c.parse().unwrap())
                    .collect()
            })
            .collect()
    };
    let add_scaled = |target: &mut Vec<i64>, source: &[i64], shift: usize, factor: i64| {
        if target.len() < source.len() + shift {
            target.resize(source.len() + shift, 0);
        }
        for (index, value) in source.iter().enumerate() {
            target[index + shift] += factor * value;
        }
    };
    // Chebyshev T_{n+1} = 2 t T_n - T_{n-1} (T_0 = 1, T_1 = t) and
    // U_{n+1} = 2 t U_n - U_{n-1} (U_0 = 1, U_1 = 2t).
    for (kind, second) in [
        (SequenceKind::ChebyshevT, vec![0, 1]),
        (SequenceKind::ChebyshevU, vec![0, 2]),
    ] {
        let rows = rows_of(kind, 6);
        assert_eq!(rows[0], vec![1]);
        assert_eq!(rows[1], second);
        for n in 1..rows.len() - 1 {
            let mut expected = Vec::new();
            add_scaled(&mut expected, &rows[n], 1, 2);
            add_scaled(&mut expected, &rows[n - 1], 0, -1);
            assert_eq!(rows[n + 1], expected, "recurrence at n = {n}");
        }
    }
    // Hermite (probabilists', OEIS A066325): He_{n+1} = t He_n - n He_{n-1}.
    let hermite = rows_of(SequenceKind::Hermite, 6);
    for n in 1..hermite.len() - 1 {
        let mut expected = Vec::new();
        add_scaled(&mut expected, &hermite[n], 1, 1);
        add_scaled(&mut expected, &hermite[n - 1], 0, -(n as i64));
        assert_eq!(hermite[n + 1], expected, "recurrence at n = {n}");
    }
    // Narayana polynomials N_n = sum_k N(n, k) t^(k-1), with N(n, k) = C(n,k) C(n,k-1)/n
    // (OEIS A001263): the rows are palindromic and real-rooted.
    let narayana = rows_of(SequenceKind::Narayana, 5);
    assert_eq!(narayana.last().unwrap(), &vec![1, 10, 20, 10, 1]);
    // Type B Eulerian polynomials B_1 = 1 + t, B_2 = 1 + 6t + t^2 (OEIS A060187).
    let type_b = rows_of(SequenceKind::TypeBEulerian, 3);
    assert!(type_b.contains(&vec![1, 6, 1]), "{type_b:?}");
}

#[test]
fn family_source_validation_errors() {
    let server = PolynomialToolsServer::new();
    let error = |request: CheckPolynomialFamilyRequest| {
        message(
            server
                .check_polynomial_family(Parameters(request))
                .err()
                .expect("expected an error"),
        )
    };
    // Exactly one source.
    let mut none = family_request(&[]);
    none.polynomials = None;
    assert!(error(none).contains("expected exactly one source"));
    let mut two = family_request(&[&[1, 1]]);
    two.text = Some("1,1".to_string());
    assert!(error(two).contains("expected exactly one source"));
    // `sequence` needs `max_n`, within the limit; `max_n` needs `sequence`.
    let mut sequence = family_request(&[]);
    sequence.polynomials = None;
    sequence.sequence = Some(SequenceKind::ChebyshevT);
    assert!(error(sequence.clone()).contains("`max_n` is required"));
    sequence.max_n = Some(MCP_MAX_SEQUENCE_N + 1);
    assert!(error(sequence.clone()).contains("max_n must be at most 200"));
    let mut stray = family_request(&[&[1, 1]]);
    stray.max_n = Some(3);
    assert!(error(stray).contains("only valid when `sequence` is provided"));
    // Rows beyond i64 are accepted: the Eulerian polynomial A_30 has coefficients
    // up to 10^30, and the sequence source keeps them exactly (the per-row
    // checks of such high degree are too slow for a debug-profile test).
    let mut large = family_request(&[]);
    large.polynomials = None;
    large.sequence = Some(SequenceKind::Eulerian);
    large.max_n = Some(30);
    let (source, parsed) = family_source_and_polynomials(&large).unwrap();
    assert_eq!(source, "sequence:Eulerian, max_n=30");
    let parsed = parsed.expect("generated rows have no parse errors");
    assert_eq!(parsed.len(), 30);
    let widest = parsed[29]
        .coefficients
        .iter()
        .max()
        .expect("nonempty row")
        .clone();
    assert!(widest > num_bigint::BigInt::from(i64::MAX), "{widest}");
}

#[test]
fn family_parse_errors_are_reported_with_markdown() {
    let server = PolynomialToolsServer::new();
    let mut request = family_request(&[]);
    request.polynomials = None;
    request.text = Some("1,1\nbad + x\n1,2,1\nalso bad".to_string());
    let Json(response) = server.check_polynomial_family(Parameters(request)).unwrap();
    assert_eq!(response.source, "text");
    assert!(!response.all_required_checks_passed);
    assert_eq!(
        response.first_failure.as_deref(),
        Some("one or more polynomials failed to parse")
    );
    assert_eq!(response.item_count, 0);
    assert!(response.items.is_empty() && response.consecutive_pairs.is_empty());
    assert_eq!(
        response
            .parse_errors
            .iter()
            .map(|e| e.index)
            .collect::<Vec<_>>(),
        vec![1, 3]
    );
    assert!(response.markdown.contains("## Parse errors"));
    assert!(response.markdown.contains("- polynomial 1: "));
    assert!(response.markdown.contains("- polynomial 3: "));
}

#[test]
fn family_required_checks_report_the_first_violated_property() {
    // Defaults: real-rootedness is required, everything else is not.
    assert_eq!(first_failure(&[&[1, 1], &[1, 2, 1]], |_| {}), None);
    assert_eq!(
        first_failure(&[&[1, 1], &[1, 0, 1]], |_| {}).as_deref(),
        Some("polynomial 1 is not real-rooted")
    );
    assert_eq!(
        first_failure(&[&[1, 0, 1]], |o| o.require_real_rooted = Some(false)),
        None
    );
    // (1 + t)^2 is real-rooted with a double root.
    assert_eq!(
        first_failure(&[&[1, 2, 1]], |o| o.require_simple_roots = Some(true)).as_deref(),
        Some("polynomial 0 does not have simple roots")
    );
    // 2 + 3t + t^2 = (1 + t)(2 + t) is not palindromic and hence not gamma-positive.
    assert_eq!(
        first_failure(&[&[2, 3, 1]], |o| o.require_palindromic = Some(true)).as_deref(),
        Some("polynomial 0 is not palindromic")
    );
    assert_eq!(
        first_failure(&[&[2, 3, 1]], |o| o.require_gamma_positive = Some(true)).as_deref(),
        Some("polynomial 0 is not gamma-positive")
    );
    // 2 + t + 2t^2 decreases and then increases; it is not real-rooted either.
    assert_eq!(
        first_failure(&[&[2, 1, 2]], |o| {
            o.require_real_rooted = Some(false);
            o.require_unimodal = Some(true);
        })
        .as_deref(),
        Some("polynomial 0 is not unimodal")
    );
    // 1 + t + 2t^2: unimodal but 1^2 < 1 * 2, so not log-concave.
    assert_eq!(
        first_failure(&[&[1, 1, 2]], |o| {
            o.require_real_rooted = Some(false);
            o.require_unimodal = Some(true);
            o.require_log_concave = Some(true);
        })
        .as_deref(),
        Some("polynomial 0 is not log-concave")
    );
    // 1 + t + t^2: log-concave, but (1/2)^2 < 1 * 1 after dividing by C(2, k).
    assert_eq!(
        first_failure(&[&[1, 1, 1]], |o| {
            o.require_real_rooted = Some(false);
            o.require_log_concave = Some(true);
            o.require_ultra_log_concave = Some(true);
        })
        .as_deref(),
        Some("polynomial 0 is not ultra-log-concave")
    );
    // Consecutive weak interlacing: the roots of 1 + t^2 are not real.
    assert_eq!(
        first_failure(&[&[1, 1], &[1, 0, 1]], |o| {
            o.require_real_rooted = Some(false);
            o.require_consecutive_weak_interlacing = Some(true);
        })
        .as_deref(),
        Some("consecutive pair 0-1 is not weakly interlacing")
    );
    assert_eq!(
        first_failure(&[&[1, 1], &[1, 2, 1]], |o| {
            o.require_consecutive_weak_interlacing = Some(true);
        }),
        None
    );
}

#[test]
fn family_recurrence_search_failures_are_reported() {
    let server = PolynomialToolsServer::new();
    let mut request = family_request(&[&[1], &[2]]);
    let mut options = family_options();
    options.find_recurrence = Some(true);
    request.options = Some(options.clone());
    let Json(response) = server.check_polynomial_family(Parameters(request)).unwrap();
    let recurrence = response.recurrence.as_ref().unwrap();
    assert!(!recurrence.found);
    assert_eq!(recurrence.status, RecurrenceSearchStatus::InvalidInput);
    assert_eq!(
        response.first_failure.as_deref(),
        Some("recurrence search did not find a recurrence")
    );
    assert!(response
        .markdown
        .contains("- Not found: need at least 3 polynomials"));

    // Constants 1, 3, 7, 2, 11, 5 satisfy no length-1 constant-coefficient recurrence.
    let mut request = family_request(&[&[1], &[3], &[7], &[2], &[11], &[5]]);
    request.options = Some(options);
    request.recurrence_options = Some(RecurrenceSearchOptionsInput {
        max_rec_len: Some(1),
        max_var_deg: Some(0),
        max_idx_deg: Some(0),
        max_diff_deg: Some(0),
        ..Default::default()
    });
    let Json(response) = server.check_polynomial_family(Parameters(request)).unwrap();
    let recurrence = response.recurrence.as_ref().unwrap();
    assert_eq!(recurrence.status, RecurrenceSearchStatus::NotFound);
    assert!(response
        .markdown
        .contains("- Not found: no recurrence found within the search bounds"));
}

#[test]
fn family_lace_limits_and_shape_errors() {
    let server = PolynomialToolsServer::new();
    let lace_response = |request: LaceCheckRequest, polynomials: &[&[i64]]| {
        let mut family = family_request(polynomials);
        family.lace = Some(request);
        let Json(response) = server.check_polynomial_family(Parameters(family)).unwrap();
        response.lace.unwrap()
    };
    let lace = |block_rows, block_cols, max_minor_size| LaceCheckRequest {
        block_rows,
        block_cols,
        max_minor_size,
        include_matrix: Some(true),
    };

    // The banded Toeplitz matrix of the real-rooted polynomial (1 + t)^3 is totally
    // nonnegative (Aissen--Schoenberg--Whitney); the defaults are one block row per
    // polynomial and `block_rows + degree + 1` columns.
    let pascal: [&[i64]; 1] = [&[1, 3, 3, 1]];
    let ok = lace_response(lace(Some(4), None, Some(3)), &pascal);
    assert!(ok.totally_nonnegative, "{:?}", ok.error);
    assert_eq!(ok.block_rows, 4);
    assert_eq!(ok.block_cols, 4 + 3 + 1);
    assert_eq!(ok.rows, 4);
    assert_eq!(ok.columns, ok.block_cols);
    let matrix = ok.matrix.as_ref().unwrap();
    assert_eq!(matrix[0], strings(&[1, 3, 3, 1, 0, 0, 0, 0]));
    assert_eq!(matrix[1], strings(&[0, 1, 3, 3, 1, 0, 0, 0]));
    let defaults = lace_response(lace(None, None, None), &pascal);
    assert_eq!((defaults.block_rows, defaults.block_cols), (1, 1 + 3 + 1));
    assert!(defaults.totally_nonnegative);

    // Too many cells: rows * columns above the cap is refused before any computation.
    let large = lace_response(lace(Some(200), Some(200), Some(1)), &pascal);
    assert!(!large.totally_nonnegative);
    assert!(large.matrix.is_none());
    assert!(large
        .error
        .as_deref()
        .unwrap()
        .contains("exceeds MCP limits"));
    // Few cells but too many minors: sum_k C(10,k) C(40,k) = C(50,10) > 250000.
    let many_minors = lace_response(lace(Some(10), Some(40), Some(10)), &pascal);
    assert!(many_minors
        .error
        .as_deref()
        .unwrap()
        .contains("exceeds MCP limits"));
    // Row count overflow is reported rather than panicking.
    let overflow = lace_response(
        lace(Some(usize::MAX), Some(2), Some(1)),
        &[&[1, 1], &[1, 2, 1]],
    );
    assert!(overflow
        .error
        .as_deref()
        .unwrap()
        .contains("row count overflow"));
    // Zero columns must not hide an enormous row count (the row allocation alone
    // would abort the process), and an overflowing default column count must not
    // panic.
    for (block_rows, block_cols) in [
        (100_000_000_000_000_000, Some(0)),
        (usize::MAX / 2, Some(0)),
        (usize::MAX, None),
        (1_000_000, Some(1)),
        (1, Some(1_000_000_000)),
    ] {
        let huge = lace_response(lace(Some(block_rows), block_cols, Some(1)), &pascal);
        assert!(!huge.totally_nonnegative, "{block_rows} {block_cols:?}");
        assert!(huge.matrix.is_none());
        assert!(
            huge.error
                .as_deref()
                .unwrap()
                .contains("exceeds MCP limits")
                || huge.error.as_deref().unwrap().contains("overflow"),
            "{:?}",
            huge.error
        );
    }
    // Truncating to one column keeps just the first column of each block row.
    let narrow = lace_response(lace(Some(2), Some(1), Some(1)), &pascal);
    assert!(narrow.totally_nonnegative);
    assert_eq!(narrow.matrix, Some(vec![strings(&[1]), strings(&[0])]));
    // An empty family has no Lace matrix: reported as an error, not as TNN.
    let empty = lace_response(lace(None, None, None), &[]);
    assert!(!empty.totally_nonnegative);
    assert!(empty.error.is_some());
    assert!(empty.matrix.is_none());
}

// -- recurrence options ------------------------------------------------------

#[test]
fn recurrence_options_map_every_field_and_imply_their_switches() {
    // Defaults.
    let (default, budget) = apply_recurrence_options(None).unwrap();
    assert!(!default.verbose);
    assert_eq!(budget.max_candidates, None);

    let input = RecurrenceSearchOptionsInput {
        skip_prefix: Some(1),
        min_rec_len: Some(2),
        max_rec_len: Some(3),
        min_var_deg: Some(0),
        max_var_deg: Some(2),
        min_idx_deg: Some(0),
        max_idx_deg: Some(1),
        min_diff_deg: Some(0),
        max_diff_deg: Some(1),
        try_inhomogeneous: Some(false),
        min_inhomo_var_deg: Some(0),
        max_inhomo_var_deg: Some(1),
        min_inhomo_idx_deg: Some(0),
        max_inhomo_idx_deg: Some(1),
        try_denominator: Some(false),
        try_alternating_sign: Some(true),
        max_denom_var_deg: Some(1),
        max_denom_idx_deg: Some(1),
        min_margin: Some(2),
        no_verify: Some(true),
        fit_extra_rows: Some(3),
        modular_prefilter: Some(false),
        max_candidates: Some(1_000),
    };
    let (options, budget) = apply_recurrence_options(Some(input)).unwrap();
    assert_eq!(budget.max_candidates, Some(1_000));
    assert_eq!(options.skip_prefix, 1);
    assert_eq!((options.min_rec_len, options.max_rec_len), (2, 3));
    assert_eq!((options.min_var_deg, options.max_var_deg), (0, 2));
    assert_eq!((options.min_idx_deg, options.max_idx_deg), (0, 1));
    assert_eq!((options.min_diff_deg, options.max_diff_deg), (0, 1));
    assert_eq!(
        (options.min_inhomo_var_deg, options.max_inhomo_var_deg),
        (0, 1)
    );
    assert_eq!(
        (options.min_inhomo_idx_deg, options.max_inhomo_idx_deg),
        (0, 1)
    );
    assert_eq!(
        (options.max_denom_var_deg, options.max_denom_idx_deg),
        (1, 1)
    );
    // Giving an inhomogeneous or denominator degree switches that search on, even
    // when `try_*` was given as false earlier in the same request.
    assert!(options.try_inhomogeneous && options.try_denominator);
    assert!(options.try_alternating_sign);
    assert_eq!(options.min_margin, 2);
    assert!(options.no_verify);
    assert_eq!(options.fit_extra_rows, 3);
    assert!(!options.modular_prefilter);
    assert!(!options.verbose);

    // Without the degree controls the switches keep their explicit values.
    let (options, _) = apply_recurrence_options(Some(RecurrenceSearchOptionsInput {
        try_inhomogeneous: Some(true),
        try_denominator: Some(true),
        ..Default::default()
    }))
    .unwrap();
    assert!(options.try_inhomogeneous && options.try_denominator);
    for (field, value) in [
        (
            RecurrenceSearchOptionsInput {
                min_inhomo_idx_deg: Some(1),
                ..Default::default()
            },
            "inhomogeneous",
        ),
        (
            RecurrenceSearchOptionsInput {
                max_inhomo_var_deg: Some(1),
                ..Default::default()
            },
            "inhomogeneous",
        ),
        (
            RecurrenceSearchOptionsInput {
                max_denom_idx_deg: Some(1),
                ..Default::default()
            },
            "denominator",
        ),
    ] {
        let (options, _) = apply_recurrence_options(Some(field)).unwrap();
        let switched = if value == "denominator" {
            options.try_denominator
        } else {
            options.try_inhomogeneous
        };
        assert!(switched, "{value}");
    }
}

#[test]
fn recurrence_options_are_validated_against_the_mcp_limits() {
    let error = |options: RecurrenceSearchOptionsInput| {
        message(
            apply_recurrence_options(Some(options))
                .err()
                .expect("expected an error"),
        )
    };
    // Empty ranges.
    assert!(error(RecurrenceSearchOptionsInput {
        min_rec_len: Some(5),
        max_rec_len: Some(3),
        ..Default::default()
    })
    .contains("minimum recurrence length 5 exceeds maximum 3"));
    assert!(error(RecurrenceSearchOptionsInput {
        min_var_deg: Some(4),
        max_var_deg: Some(1),
        ..Default::default()
    })
    .contains("minimum variable degree 4 exceeds maximum 1"));
    assert!(error(RecurrenceSearchOptionsInput {
        min_inhomo_var_deg: Some(9),
        ..Default::default()
    })
    .contains("minimum inhomogeneous variable degree 9 exceeds maximum 3"));
    // Limits are inclusive.
    let at_length_limit = apply_recurrence_options(Some(RecurrenceSearchOptionsInput {
        max_rec_len: Some(MCP_MAX_RECURRENCE_LENGTH),
        max_var_deg: Some(0),
        max_idx_deg: Some(0),
        max_diff_deg: Some(0),
        ..Default::default()
    }));
    assert!(at_length_limit.is_ok());
    assert!(error(RecurrenceSearchOptionsInput {
        max_rec_len: Some(MCP_MAX_RECURRENCE_LENGTH + 1),
        ..Default::default()
    })
    .contains("exceeds the MCP limit 32"));
    for (options, name) in [
        (
            RecurrenceSearchOptionsInput {
                max_var_deg: Some(MCP_MAX_RECURRENCE_DEGREE + 1),
                ..Default::default()
            },
            "variable",
        ),
        (
            RecurrenceSearchOptionsInput {
                max_diff_deg: Some(MCP_MAX_RECURRENCE_DEGREE + 1),
                ..Default::default()
            },
            "derivative",
        ),
        (
            RecurrenceSearchOptionsInput {
                max_denom_idx_deg: Some(MCP_MAX_RECURRENCE_DEGREE + 1),
                ..Default::default()
            },
            "denominator index",
        ),
    ] {
        assert!(
            error(options).contains(&format!("maximum {name} degree 17 exceeds")),
            "{name}"
        );
    }
    assert!(error(RecurrenceSearchOptionsInput {
        skip_prefix: Some(MCP_MAX_RECURRENCE_INPUT_ROWS + 1),
        ..Default::default()
    })
    .contains("skip_prefix and fit_extra_rows must be at most 256"));
    assert!(error(RecurrenceSearchOptionsInput {
        fit_extra_rows: Some(MCP_MAX_RECURRENCE_INPUT_ROWS + 1),
        ..Default::default()
    })
    .contains("skip_prefix and fit_extra_rows must be at most 256"));
    assert!(error(RecurrenceSearchOptionsInput {
        min_margin: Some(20_001),
        ..Default::default()
    })
    .contains("min_margin 20001 exceeds"));

    // A huge search space is refused by candidate count, unless a candidate budget
    // caps it, in which case the unknown-count limit still applies.
    let huge = RecurrenceSearchOptionsInput {
        max_rec_len: Some(32),
        max_var_deg: Some(16),
        max_idx_deg: Some(16),
        max_diff_deg: Some(16),
        ..Default::default()
    };
    assert!(error(huge.clone()).contains("may inspect"));
    assert!(error(RecurrenceSearchOptionsInput {
        max_candidates: Some(10),
        ..huge
    })
    .contains("unknowns; the MCP limit is 20000"));
}

fn bivar_of_size(rows: usize, columns: usize) -> BivarPoly {
    BivarPoly {
        coeffs: vec![vec![BigRational::from_integer(1.into()); columns]; rows],
    }
}

#[test]
fn recurrence_generation_validation_enforces_each_limit() {
    let one = |value: i64| vec![BigRational::from_integer(value.into())];
    let rule = |terms: Vec<RecurrenceTerm>| Recurrence {
        terms,
        denominator: None,
        inhomogeneous: None,
    };
    let term =
        |offset: usize, deriv: usize, coeff: BivarPoly| RecurrenceTerm::new(offset, deriv, coeff);
    let ok_rule = rule(vec![term(1, 0, bivar_of_size(1, 1))]);
    assert!(validate_recurrence_for_generation(&ok_rule, &[one(1)]).is_ok());

    let error = |recurrence: &Recurrence, rows: &[Vec<BigRational>]| {
        message(
            validate_recurrence_for_generation(recurrence, rows)
                .err()
                .expect("expected an error"),
        )
    };
    // Number of terms.
    let many_terms = rule(
        (0..=MCP_MAX_RECURRENCE_LENGTH)
            .map(|_| term(1, 0, bivar_of_size(1, 1)))
            .collect(),
    );
    assert!(error(&many_terms, &[one(1)]).contains("has 33 terms"));
    // Rows and their lengths.
    let rows = vec![one(1); MCP_MAX_RECURRENCE_INPUT_ROWS + 1];
    assert!(error(&ok_rule, &rows).contains("initial rows"));
    let long_row =
        vec![BigRational::from_integer(1.into()); MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL + 1];
    assert!(error(&ok_rule, &[long_row]).contains("initial polynomial exceeds"));
    // Term offset, derivative order and bivariate degrees.
    let far = rule(vec![term(
        MCP_MAX_RECURRENCE_LENGTH + 1,
        0,
        bivar_of_size(1, 1),
    )]);
    assert!(error(&far, &[one(1)]).contains("has offset 33"));
    let deep = rule(vec![term(
        1,
        MCP_MAX_RECURRENCE_DEGREE + 1,
        bivar_of_size(1, 1),
    )]);
    assert!(error(&deep, &[one(1)]).contains("derivative order 17"));
    let wide = rule(vec![term(
        1,
        0,
        bivar_of_size(1, MCP_MAX_RECURRENCE_DEGREE + 2),
    )]);
    assert!(error(&wide, &[one(1)]).contains("recurrence term 0 exceeds the MCP bivariate degree"));
    let tall = Recurrence {
        denominator: Some(bivar_of_size(MCP_MAX_RECURRENCE_DEGREE + 2, 1)),
        ..ok_rule.clone()
    };
    assert!(error(&tall, &[one(1)]).contains("recurrence denominator exceeds"));
    let inhomogeneous = Recurrence {
        inhomogeneous: Some(bivar_of_size(1, MCP_MAX_RECURRENCE_DEGREE + 2)),
        ..ok_rule.clone()
    };
    assert!(error(&inhomogeneous, &[one(1)]).contains("inhomogeneous term exceeds"));
    // Total stored coefficients: 256 rows of 4096 coefficients exceed 65536.
    let big_rows =
        vec![
            vec![BigRational::from_integer(1.into()); MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL];
            MCP_MAX_RECURRENCE_INPUT_ROWS
        ];
    assert!(error(&ok_rule, &big_rows).contains("stored coefficients"));
}

// -- parsing limits -----------------------------------------------------------

#[test]
fn polynomial_input_forms_and_limits() {
    // Exactly one of `coefficients` and `expression`.
    let both = BigIntPolynomialInput {
        coefficients: Some(vec![1.into()]),
        expression: Some("1".to_string()),
    };
    assert!(parse_one(&both)
        .err()
        .expect("expected an error")
        .contains("got both"));
    let neither = BigIntPolynomialInput {
        coefficients: None,
        expression: None,
    };
    assert!(parse_one(&neither)
        .err()
        .expect("expected an error")
        .contains("expected exactly one"));
    // Length limit (inclusive) for coefficients and for the expression text.
    let at_limit = vec![1; MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL];
    assert!(parse_one(&coeffs(&at_limit)).is_ok());
    let over = vec![1; MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL + 1];
    assert!(parse_one(&coeffs(&over))
        .err()
        .expect("expected an error")
        .contains("the MCP limit is 4096"));
    let long_text = BigIntPolynomialInput {
        coefficients: None,
        expression: Some("1".repeat(MCP_MAX_TEXT_BYTES + 1)),
    };
    assert!(parse_one(&long_text)
        .err()
        .expect("expected an error")
        .contains("polynomial expression is 1048577 bytes"));
    let degree_expression = BigIntPolynomialInput {
        coefficients: None,
        expression: Some(format!("t^{}", MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL)),
    };
    assert!(parse_one(&degree_expression)
        .err()
        .expect("expected an error")
        .contains("4097 coefficients"));

    // BigInt coefficients: text length, syntax and per-coefficient context.
    let bad = bigint_coeffs(&["1", "x1"]);
    let error = message(
        parse_required_bigint_polynomial(&bad, "p")
            .err()
            .expect("expected an error"),
    );
    assert!(
        error.starts_with("p: coefficient 1: invalid integer 'x1'"),
        "{error}"
    );
    let huge = bigint_coeffs(&[&"9".repeat(16_385)]);
    assert!(message(
        parse_required_bigint_polynomial(&huge, "q")
            .err()
            .expect("expected an error")
    )
    .contains("integer coefficient is 16385 bytes"));
    let both = BigIntPolynomialInput {
        coefficients: Some(vec![1.into()]),
        expression: Some("1".to_string()),
    };
    assert!(message(
        parse_required_bigint_polynomial(&both, "p")
            .err()
            .expect("expected an error")
    )
    .contains("got both"));
    let neither = BigIntPolynomialInput {
        coefficients: None,
        expression: None,
    };
    assert!(message(
        parse_required_bigint_polynomial(&neither, "p")
            .err()
            .expect("expected an error")
    )
    .contains("expected exactly one"));
    let expression = BigIntPolynomialInput {
        coefficients: None,
        expression: Some("2t^2 + 100000000000000000000".to_string()),
    };
    assert_eq!(
        parse_required_bigint_polynomial(&expression, "p")
            .unwrap()
            .last(),
        Some(&2.into())
    );
}

#[test]
fn batch_parsing_enforces_size_limits() {
    let error = |input: BigIntPolynomialBatchInput| {
        message(parse_bigint_batch(&input).err().expect("expected an error"))
    };
    // Exactly one of `polynomials` and `text`.
    assert!(error(BigIntPolynomialBatchInput {
        polynomials: Some(vec![coeffs(&[1])]),
        text: Some("1".to_string()),
    })
    .contains("got both"));
    assert!(error(BigIntPolynomialBatchInput {
        polynomials: None,
        text: None,
    })
    .contains("expected exactly one"));
    // Number of polynomials (inclusive limit); blank and `#` comment lines do not count.
    let at_limit = vec![coeffs(&[1]); MCP_MAX_BATCH_POLYNOMIALS];
    assert_eq!(
        parse_bigint_batch(&BigIntPolynomialBatchInput {
            polynomials: Some(at_limit),
            text: None,
        })
        .unwrap()
        .len(),
        MCP_MAX_BATCH_POLYNOMIALS
    );
    assert!(error(BigIntPolynomialBatchInput {
        polynomials: Some(vec![coeffs(&[1]); MCP_MAX_BATCH_POLYNOMIALS + 1]),
        text: None,
    })
    .contains("batch has 513 polynomials"));
    let text_lines = |count: usize| {
        let mut text = String::from("# header\n\n");
        text.push_str(&"1,1\n".repeat(count));
        text
    };
    assert_eq!(
        parse_bigint_batch(&BigIntPolynomialBatchInput {
            polynomials: None,
            text: Some(text_lines(MCP_MAX_BATCH_POLYNOMIALS)),
        })
        .unwrap()
        .len(),
        MCP_MAX_BATCH_POLYNOMIALS
    );
    assert!(error(BigIntPolynomialBatchInput {
        polynomials: None,
        text: Some(text_lines(MCP_MAX_BATCH_POLYNOMIALS + 1)),
    })
    .contains("batch has 513 polynomials"));
    assert!(error(BigIntPolynomialBatchInput {
        polynomials: None,
        text: Some(" ".repeat(MCP_MAX_TEXT_BYTES + 1)),
    })
    .contains("polynomial batch is 1048577 bytes"));
    // One over-long polynomial is an item error and does not poison its neighbors.
    let batch = parse_bigint_batch(&BigIntPolynomialBatchInput {
        polynomials: Some(vec![
            coeffs(&[1, 1]),
            coeffs(&vec![1; MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL + 1]),
        ]),
        text: None,
    })
    .unwrap();
    assert!(batch[0].is_ok());
    assert!(batch[1]
        .as_ref()
        .err()
        .expect("expected an error")
        .contains("coefficients"));
    let failed: Vec<usize> = batch
        .iter()
        .enumerate()
        .filter_map(|(index, item)| item.is_err().then_some(index))
        .collect();
    assert_eq!(failed, vec![1]);
    // The total coefficient count is limited across the batch: 17 * 4096 > 65536.
    assert!(error(BigIntPolynomialBatchInput {
        polynomials: Some(vec![
            coeffs(&vec![1; MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL]);
            MCP_MAX_TOTAL_COEFFICIENTS
                / MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL
                + 1
        ]),
        text: None,
    })
    .contains("batch has 69632 coefficients"));

    // The BigInt batch parser has the same rules.
    let bigint_error = |input: BigIntPolynomialBatchInput| {
        message(parse_bigint_batch(&input).err().expect("expected an error"))
    };
    assert!(bigint_error(BigIntPolynomialBatchInput {
        polynomials: Some(vec![bigint_coeffs(&["1"])]),
        text: Some("1".to_string()),
    })
    .contains("got both"));
    assert!(bigint_error(BigIntPolynomialBatchInput {
        polynomials: None,
        text: None,
    })
    .contains("expected exactly one"));
    assert!(bigint_error(BigIntPolynomialBatchInput {
        polynomials: Some(vec![bigint_coeffs(&["1"]); MCP_MAX_BATCH_POLYNOMIALS + 1]),
        text: None,
    })
    .contains("batch has 513 polynomials"));
    assert!(bigint_error(BigIntPolynomialBatchInput {
        polynomials: None,
        text: Some(text_lines(MCP_MAX_BATCH_POLYNOMIALS + 1)),
    })
    .contains("batch has 513 polynomials"));
    assert!(bigint_error(BigIntPolynomialBatchInput {
        polynomials: None,
        text: Some(" ".repeat(MCP_MAX_TEXT_BYTES + 1)),
    })
    .contains("polynomial batch is 1048577 bytes"));
    let long_row = "1,".repeat(MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL + 1);
    let batch = parse_bigint_batch(&BigIntPolynomialBatchInput {
        polynomials: None,
        text: Some(format!("{}1\n1,1", long_row)),
    })
    .unwrap();
    assert!(batch[0]
        .as_ref()
        .err()
        .expect("expected an error")
        .contains("coefficients"));
    assert!(batch[1].is_ok());
    let wide = bigint_coeffs(&vec!["1"; MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL]);
    assert!(bigint_error(BigIntPolynomialBatchInput {
        polynomials: Some(vec![
            wide;
            MCP_MAX_TOTAL_COEFFICIENTS
                / MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL
                + 1
        ]),
        text: None,
    })
    .contains("batch has 69632 coefficients"));
}

#[test]
fn recurrence_input_parsing_enforces_its_own_limits() {
    let error = |input: FindRecurrenceRequest| {
        message(
            parse_recurrence_batch_rational(&input)
                .err()
                .expect("expected an error"),
        )
    };
    let rows = |count: usize| -> Vec<Vec<RationalCoefficientInput>> { vec![vec![1.into()]; count] };
    assert!(error(FindRecurrenceRequest::default()).contains("expected exactly one of"));
    assert!(error(FindRecurrenceRequest {
        coefficients: Some(rows(1)),
        text: Some("1".to_string()),
        ..Default::default()
    })
    .contains("expected exactly one of"));

    // Row limits for each source form.
    let too_many = MCP_MAX_RECURRENCE_INPUT_ROWS + 1;
    assert!(error(FindRecurrenceRequest {
        polynomials: Some(vec![coeffs(&[1]); too_many]),
        ..Default::default()
    })
    .contains("recurrence input has 257 rows"));
    assert!(error(FindRecurrenceRequest {
        coefficients: Some(rows(too_many)),
        ..Default::default()
    })
    .contains("recurrence input has 257 rows"));
    assert!(error(FindRecurrenceRequest {
        expressions: Some(vec!["1".to_string(); too_many]),
        ..Default::default()
    })
    .contains("recurrence input has 257 rows"));
    assert!(error(FindRecurrenceRequest {
        text: Some("1\n".repeat(too_many)),
        ..Default::default()
    })
    .contains("recurrence input has 257 rows"));
    assert!(error(FindRecurrenceRequest {
        text: Some(" ".repeat(MCP_MAX_TEXT_BYTES + 1)),
        ..Default::default()
    })
    .contains("recurrence polynomial batch is 1048577 bytes"));
    // At the row limit every form is accepted.
    let at_limit = MCP_MAX_RECURRENCE_INPUT_ROWS;
    for request in [
        FindRecurrenceRequest {
            polynomials: Some(vec![coeffs(&[1]); at_limit]),
            ..Default::default()
        },
        FindRecurrenceRequest {
            coefficients: Some(rows(at_limit)),
            ..Default::default()
        },
        FindRecurrenceRequest {
            expressions: Some(vec!["1".to_string(); at_limit]),
            ..Default::default()
        },
        FindRecurrenceRequest {
            text: Some("1\n".repeat(at_limit)),
            ..Default::default()
        },
    ] {
        assert_eq!(
            parse_recurrence_batch_rational(&request).unwrap().len(),
            at_limit
        );
    }
    // Total coefficient limit: 17 rows of 4096 coefficients.
    let wide_row = vec![RationalCoefficientInput::Integer(1); MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL];
    assert!(error(FindRecurrenceRequest {
        coefficients: Some(vec![
            wide_row.clone();
            MCP_MAX_TOTAL_COEFFICIENTS
                / MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL
                + 1
        ]),
        ..Default::default()
    })
    .contains("recurrence input has 69632 coefficients"));

    // Rational coefficients: exact parsing, item-level errors with indices.
    let parsed = parse_recurrence_batch_rational(&FindRecurrenceRequest {
        coefficients: Some(vec![
            vec![
                RationalCoefficientInput::Text("3/6".to_string()),
                2.into(),
                RationalCoefficientInput::Integer(0),
            ],
            vec![RationalCoefficientInput::Text("not a number".to_string())],
            vec![RationalCoefficientInput::Text("1".repeat(16_385))],
            vec![RationalCoefficientInput::Integer(1); MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL + 1],
        ]),
        ..Default::default()
    })
    .unwrap();
    // 3/6 reduces to 1/2 and the trailing zero is dropped.
    let first = parsed[0].as_ref().unwrap();
    assert_eq!(first.len(), 2);
    assert_eq!(first[0].to_string(), "1/2");
    assert!(parsed[1]
        .as_ref()
        .err()
        .expect("expected an error")
        .starts_with("coefficient 0: invalid rational"));
    assert!(parsed[2]
        .as_ref()
        .err()
        .expect("expected an error")
        .contains("rational coefficient is 16385 bytes"));
    assert!(parsed[3]
        .as_ref()
        .err()
        .expect("expected an error")
        .contains("the MCP limit is 4096"));
}

#[test]
fn bigint_matrix_parsing_enforces_limits() {
    let entry = |value: i64| BigIntCoefficientInput::Integer(value);
    let parsed =
        parse_bigint_matrix(&[vec![entry(1), entry(2)], vec![entry(3)]], "matrix").unwrap();
    assert_eq!(parsed.len(), 2);
    assert_eq!(parsed[1], vec![3.into()]);
    let too_many_rows = vec![vec![entry(1)]; MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL + 1];
    assert!(message(
        parse_bigint_matrix(&too_many_rows, "matrix")
            .err()
            .expect("expected an error")
    )
    .contains("`matrix` has 4097 rows"));
    let too_many_cells = vec![vec![entry(1); 4096]; 17];
    assert!(message(
        parse_bigint_matrix(&too_many_cells, "fixed_table")
            .err()
            .expect("expected an error")
    )
    .contains("`fixed_table` has 69632 entries"));
    let long_row = vec![vec![entry(1); MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL + 1]];
    assert!(message(
        parse_bigint_matrix(&long_row, "matrix")
            .err()
            .expect("expected an error")
    )
    .contains("polynomial has 4097 coefficients"));
}

// -- handler error branches ----------------------------------------------------

#[test]
fn ehrhart_hstar_known_values_and_input_errors() {
    let server = PolynomialToolsServer::new();
    let request = |mode| EhrhartHstarRequest {
        mode,
        hstar: None,
        ehrhart_coefficients: None,
        numerator_coefficients: None,
        denominator: None,
    };
    let error = |request: EhrhartHstarRequest| {
        message(
            server
                .ehrhart_hstar(Parameters(request))
                .err()
                .expect("expected an error"),
        )
    };

    // The unit square has Ehrhart polynomial (t + 1)^2 = 1 + 2t + t^2 and
    // h*-polynomial 1 + t; (1 + 2t)^3 with h* = (1, 4, 1) is the unit cube's.
    let Json(cube) = server
        .ehrhart_hstar(Parameters(EhrhartHstarRequest {
            ehrhart_coefficients: Some(
                ["1", "3", "3", "1"].iter().map(|s| s.to_string()).collect(),
            ),
            ..request(EhrhartHstarMode::EhrhartToHstar)
        }))
        .unwrap();
    assert_eq!(
        cube.hstar,
        Some(vec!["1".into(), "4".into(), "1".into(), "0".into()])
    );
    // The same polynomial with a common denominator: (2 + 6t + 6t^2 + 2t^3) / 2.
    let Json(scaled) = server
        .ehrhart_hstar(Parameters(EhrhartHstarRequest {
            numerator_coefficients: Some(vec![2.into(), 6.into(), 6.into(), 2.into()]),
            denominator: Some(2.into()),
            ..request(EhrhartHstarMode::EhrhartToHstar)
        }))
        .unwrap();
    assert_eq!(scaled.hstar, cube.hstar);

    assert!(error(request(EhrhartHstarMode::HstarToEhrhart)).contains("`hstar` is required"));
    assert!(error(request(EhrhartHstarMode::EhrhartToHstar))
        .contains("provide either `ehrhart_coefficients` or both"));
    assert!(error(EhrhartHstarRequest {
        ehrhart_coefficients: Some(vec!["1".into()]),
        numerator_coefficients: Some(vec![1.into()]),
        denominator: Some(1.into()),
        ..request(EhrhartHstarMode::EhrhartToHstar)
    })
    .contains("provide either"));
    assert!(error(EhrhartHstarRequest {
        numerator_coefficients: Some(vec![1.into()]),
        denominator: Some(0.into()),
        ..request(EhrhartHstarMode::EhrhartToHstar)
    })
    .contains("`denominator` must be nonzero"));
    assert!(error(EhrhartHstarRequest {
        ehrhart_coefficients: Some(vec!["1/x".into()]),
        ..request(EhrhartHstarMode::EhrhartToHstar)
    })
    .contains("invalid rational"));
    assert!(error(EhrhartHstarRequest {
        ehrhart_coefficients: Some(vec!["1".into(); MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL + 1]),
        ..request(EhrhartHstarMode::EhrhartToHstar)
    })
    .contains("Ehrhart polynomial has 4097 coefficients"));
    assert!(error(EhrhartHstarRequest {
        ehrhart_coefficients: Some(vec!["1".repeat(16_385)]),
        ..request(EhrhartHstarMode::EhrhartToHstar)
    })
    .contains("Ehrhart coefficient is 16385 bytes"));
    // A non-integral h*-vector is rejected cleanly: L(t) = 1/2 + t has h* = (1/2, 1/2).
    assert!(!error(EhrhartHstarRequest {
        ehrhart_coefficients: Some(vec!["1/2".into(), "1".into()]),
        ..request(EhrhartHstarMode::EhrhartToHstar)
    })
    .is_empty());
}

#[test]
fn sieving_tools_reject_invalid_parameters() {
    let server = PolynomialToolsServer::new();
    let cyclic_error = |polynomial: &[&str], order: usize, orbit_sizes: Option<Vec<usize>>| {
        message(
            server
                .cyclic_sieving(Parameters(CyclicSievingRequest {
                    polynomial: bigint_coeffs(polynomial),
                    order,
                    fixed_counts: None,
                    orbit_sizes,
                }))
                .err()
                .expect("expected an error"),
        )
    };
    assert!(cyclic_error(&["1", "1"], 0, None).contains("`order` must be positive"));
    assert!(cyclic_error(&["1", "1"], 10_001, None).contains("exceeds the MCP limit 10000"));
    assert!(cyclic_error(
        &["1", "1"],
        2,
        Some(vec![1; MCP_MAX_TOTAL_COEFFICIENTS + 1])
    )
    .contains("`orbit_sizes` may have at most 65536 entries"));
    assert!(cyclic_error(&["1", "x"], 2, None).starts_with("polynomial: "));

    let bicyclic = |order1, order2, matrix: Vec<Vec<i64>>| {
        server.bicyclic_sieving(Parameters(BicyclicSievingRequest {
            matrix: matrix
                .into_iter()
                .map(|row| row.into_iter().map(BigIntCoefficientInput::from).collect())
                .collect(),
            order1,
            order2,
            fixed_table: None,
        }))
    };
    assert!(message(
        bicyclic(0, 2, vec![vec![1]])
            .err()
            .expect("expected an error")
    )
    .contains("`order1` and `order2` must be positive"));
    assert!(message(
        bicyclic(2, 0, vec![vec![1]])
            .err()
            .expect("expected an error")
    )
    .contains("must be positive"));
    let table_error = server.bicyclic_sieving(Parameters(BicyclicSievingRequest {
        matrix: vec![vec![1.into(), 1.into()], vec![1.into(), 1.into()]],
        order1: 2,
        order2: 2,
        fixed_table: Some(vec![vec![BigIntCoefficientInput::Text("z".to_string())]]),
    }));
    assert!(table_error.is_err());

    let sequence = |request: CyclicSievingSequenceRequest| {
        message(
            server
                .cyclic_sieving_sequence(Parameters(request))
                .err()
                .expect("expected an error"),
        )
    };
    let base = || CyclicSievingSequenceRequest {
        polynomials: Some(vec![bigint_coeffs(&["1", "1"])]),
        text: None,
        first_index: None,
        offsets: None,
        fixed_counts: None,
    };
    assert!(sequence(CyclicSievingSequenceRequest {
        offsets: Some(Vec::new()),
        ..base()
    })
    .contains("`offsets` may not be empty"));
    // An overflowing first index or a huge candidate order is a parameter error,
    // not a panic or an enormous allocation.
    for (first_index, offsets) in [
        (isize::MAX, vec![0]),
        (isize::MAX - 1, vec![0]),
        (0, vec![isize::MAX]),
        (0, vec![10_001]),
        (1_000_000_000_000, vec![0, 1]),
    ] {
        let error = sequence(CyclicSievingSequenceRequest {
            polynomials: Some(vec![
                bigint_coeffs(&["1", "1"]),
                bigint_coeffs(&["1", "1", "1"]),
            ]),
            first_index: Some(first_index),
            offsets: Some(offsets.clone()),
            ..base()
        });
        assert!(
            error.contains("overflows") || error.contains("exceeds the MCP limit"),
            "{first_index} {offsets:?}: {error}"
        );
    }
    assert!(sequence(CyclicSievingSequenceRequest {
        text: Some("1,1\nbad + x".to_string()),
        polynomials: None,
        ..base()
    })
    .contains("polynomial 1: "));
    assert!(sequence(CyclicSievingSequenceRequest {
        fixed_counts: Some(vec![CyclicSievingFixedCountsInput {
            index: 0,
            order: 0,
            counts: Vec::new(),
        }]),
        ..base()
    })
    .contains("fixed-count order must be positive"));
    assert!(sequence(CyclicSievingSequenceRequest {
        fixed_counts: Some(vec![CyclicSievingFixedCountsInput {
            index: 0,
            order: 3,
            counts: vec![1.into(), 1.into()],
        }]),
        ..base()
    })
    .contains("expected 3 counts, got 2"));
}

#[test]
fn oeis_tools_reject_unknown_and_experimental_requests() {
    let server = PolynomialToolsServer::new();
    let unknown = message(
        server
            .get_oeis_sequence(Parameters(GetOeisSequenceRequest {
                id: "A000000".to_string(),
            }))
            .err()
            .expect("expected an error"),
    );
    assert!(unknown.contains("unknown bundled OEIS sequence: A000000"));
    let generate = |id: &str, rows: usize, first_row: Option<i64>, experimental: Option<bool>| {
        server.generate_oeis_rows(Parameters(GenerateOeisRowsRequest {
            id: id.to_string(),
            rows,
            first_row,
            include_experimental: experimental,
        }))
    };
    assert!(message(
        generate("A008292", 201, None, None)
            .err()
            .expect("expected an error")
    )
    .contains("rows must be at most 200"));
    assert!(message(
        generate("A000000", 1, None, None)
            .err()
            .expect("expected an error")
    )
    .contains("unknown bundled OEIS sequence"));

    // Eulerian rows A008292: row n = 4 is 1 + 11t + 11t^2 + t^3; `first_row` shifts the
    // window.  An experimental entry needs the explicit opt-in.
    let Json(rows) = generate("A008292", 2, Some(3), None).unwrap();
    assert_eq!(rows.first_row, 3);
    assert_eq!(rows.rows[0].n, 3);
    assert_eq!(rows.rows[1].n, 4);
    assert_eq!(rows.rows[0].coefficients, vec!["1", "4", "1"]);
    assert_eq!(rows.rows[1].coefficients, vec!["1", "11", "11", "1"]);
    let Json(listing) = server
        .list_oeis_sequences(Parameters(ListOeisSequencesRequest {
            include_experimental: Some(true),
        }))
        .unwrap();
    let Json(stable) = server
        .list_oeis_sequences(Parameters(ListOeisSequencesRequest::default()))
        .unwrap();
    assert_eq!(listing.count, listing.sequences.len());
    assert!(stable.count <= listing.count);
    if let Some(experimental) = listing
        .sequences
        .iter()
        .find(|sequence| !stable.sequences.iter().any(|s| s.id == sequence.id))
    {
        let id = experimental.id.clone();
        assert!(message(
            generate(&id, 1, None, None)
                .err()
                .expect("expected an error")
        )
        .contains("is experimental; set include_experimental to true"));
        assert!(generate(&id, 1, None, Some(true)).is_ok());
    }
}

#[test]
fn json_request_tools_reject_malformed_and_oversized_payloads() {
    let server = PolynomialToolsServer::new();
    let malformed = json!([1, 2, 3]);
    let oversized = json!({ "padding": "x".repeat(MCP_MAX_TEXT_BYTES + 1) });
    use crate::requests::{
        ImplicitRelationSearchToolRequest, VerifyImplicitRelationSearchToolRequest,
    };
    let call_spaces = |value: serde_json::Value| {
        message(
            server
                .find_recurrence_spaces(Parameters(RecurrenceSpacesToolRequest { request: value }))
                .err()
                .expect("expected an error"),
        )
    };
    assert!(!call_spaces(malformed.clone()).is_empty());
    assert!(call_spaces(oversized.clone()).contains("recurrence-space request is"));

    let verify_spaces = |report: serde_json::Value, expected: Option<serde_json::Value>| {
        message(
            server
                .verify_recurrence_spaces(Parameters(VerifyRecurrenceSpacesToolRequest {
                    report,
                    expected_request: expected,
                }))
                .err()
                .expect("expected an error"),
        )
    };
    assert!(!verify_spaces(malformed.clone(), None).is_empty());
    assert!(verify_spaces(oversized.clone(), None).contains("recurrence-space report is"));
    // A well-formed report object that is not a report.
    assert!(!verify_spaces(json!({"request": {}}), Some(oversized.clone())).is_empty());

    macro_rules! request_tool {
        ($method:ident, $request:ident, $label:expr) => {{
            assert!(
                !message(
                    server
                        .$method(Parameters($request {
                            request: malformed.clone()
                        }))
                        .err()
                        .expect("expected an error")
                )
                .is_empty(),
                "{}",
                $label
            );
            let error = message(
                server
                    .$method(Parameters($request {
                        request: oversized.clone(),
                    }))
                    .err()
                    .expect("expected an error"),
            );
            assert!(
                error.contains("the MCP limit is 1048576"),
                "{}: {error}",
                $label
            );
        }};
    }
    macro_rules! report_tool {
        ($method:ident, $request:ident, $label:expr) => {{
            assert!(
                !message(
                    server
                        .$method(Parameters($request {
                            report: malformed.clone()
                        }))
                        .err()
                        .expect("expected an error")
                )
                .is_empty(),
                "{}",
                $label
            );
            let error = message(
                server
                    .$method(Parameters($request {
                        report: oversized.clone(),
                    }))
                    .err()
                    .expect("expected an error"),
            );
            assert!(
                error.contains("the MCP limit is 1048576"),
                "{}: {error}",
                $label
            );
        }};
    }
    request_tool!(
        transform_recurrence,
        RecurrenceTransformToolRequest,
        "transform"
    );
    request_tool!(
        exact_root_enclosures,
        RootEnclosuresToolRequest,
        "enclosures"
    );
    request_tool!(
        transform_recurrence_implicit,
        RecurrenceImplicitToolRequest,
        "implicit"
    );
    request_tool!(
        search_recurrence_representations,
        RecurrenceRepresentationSearchToolRequest,
        "representations"
    );
    request_tool!(
        search_implicit_relations,
        ImplicitRelationSearchToolRequest,
        "implicit search"
    );
    report_tool!(
        verify_recurrence_transform,
        VerifyRecurrenceTransformToolRequest,
        "verify transform"
    );
    report_tool!(
        verify_recurrence_implicit,
        VerifyRecurrenceImplicitToolRequest,
        "verify implicit"
    );
    report_tool!(
        verify_recurrence_representation_search,
        VerifyRecurrenceRepresentationSearchToolRequest,
        "verify representations"
    );
    report_tool!(
        verify_implicit_relation_search,
        VerifyImplicitRelationSearchToolRequest,
        "verify implicit search"
    );

    // Generating-description assessment.
    assert!(!message(
        server
            .assess_generating_description(Parameters(AssessGeneratingDescriptionToolRequest {
                description: malformed.clone(),
            }))
            .err()
            .expect("expected an error")
    )
    .is_empty());
    assert!(message(
        server
            .assess_generating_description(Parameters(AssessGeneratingDescriptionToolRequest {
                description: oversized,
            }))
            .err()
            .expect("expected an error")
    )
    .contains("generating description is"));
}

#[test]
fn recurrence_row_generation_and_search_reject_invalid_requests() {
    let server = PolynomialToolsServer::new();
    // P(n) = P(n-1) + P(n-2) on constants, seeds 1, 1 (Fibonacci, OEIS A000045).
    let Json(found) = server
        .find_recurrence(Parameters(FindRecurrenceRequest {
            coefficients: Some(
                [1, 1, 2, 3, 5, 8, 13]
                    .iter()
                    .map(|&c| vec![RationalCoefficientInput::Integer(c)])
                    .collect(),
            ),
            ..Default::default()
        }))
        .unwrap();
    let recurrence_json = found.recurrence_json.clone().expect("recurrence json");
    let generate = |rows: Option<usize>, additional: Option<usize>| {
        server.generate_recurrence_rows(Parameters(GenerateRecurrenceRowsRequest {
            recurrence_json: recurrence_json.clone(),
            rows,
            additional,
        }))
    };
    // Row counts: exactly one of `rows` and `additional`; at most 200 rows.
    assert!(
        message(generate(None, None).err().expect("expected an error"))
            .contains("exactly one of rows or additional")
    );
    assert!(
        message(generate(Some(3), Some(3)).err().expect("expected an error"))
            .contains("exactly one of rows or additional")
    );
    assert!(
        message(generate(Some(201), None).err().expect("expected an error"))
            .contains("requested 201 recurrence rows; the MCP limit is 200")
    );
    let Json(more) = generate(None, Some(3)).unwrap();
    assert!(more.row_count >= 3);
    // The generated rows continue the Fibonacci numbers.
    let Json(rows) = generate(Some(10), None).unwrap();
    let values: Vec<String> = rows.polynomials.iter().map(|row| row[0].clone()).collect();
    assert_eq!(
        values,
        ["1", "1", "2", "3", "5", "8", "13", "21", "34", "55"]
    );
    assert!(message(
        server
            .generate_recurrence_rows(Parameters(GenerateRecurrenceRowsRequest {
                recurrence_json: "not json".to_string(),
                rows: Some(3),
                additional: None,
            }))
            .err()
            .expect("expected an error")
    )
    .contains("failed to parse recurrence JSON"));
    assert!(message(
        server
            .generate_recurrence_rows(Parameters(GenerateRecurrenceRowsRequest {
                recurrence_json: "x".repeat(MCP_MAX_TEXT_BYTES + 1),
                rows: Some(3),
                additional: None,
            }))
            .err()
            .expect("expected an error")
    )
    .contains("recurrence JSON is"));

    // find_recurrence reports unparsable rows and a too-short input as invalid input
    // (not as protocol errors), and representations reject unparsable rows.
    let Json(bad) = server
        .find_recurrence(Parameters(FindRecurrenceRequest {
            text: Some("1\nbad + x\n2\n3".to_string()),
            ..Default::default()
        }))
        .unwrap();
    assert_eq!(bad.status, RecurrenceSearchStatus::InvalidInput);
    assert_eq!(bad.parse_errors.len(), 1);
    assert_eq!(bad.parse_errors[0].index, 1);
    let Json(short) = server
        .find_recurrence(Parameters(FindRecurrenceRequest {
            text: Some("1\n2".to_string()),
            ..Default::default()
        }))
        .unwrap();
    assert_eq!(short.status, RecurrenceSearchStatus::InvalidInput);
    assert_eq!(short.error.as_deref(), Some("need at least 3 polynomials"));
    assert!(message(
        server
            .find_recurrence_representations(Parameters(FindRecurrenceRequest {
                text: Some("1\nbad + x\n2\n3".to_string()),
                ..Default::default()
            }))
            .err()
            .expect("expected an error")
    )
    .contains("one or more polynomials failed to parse"));
}

#[test]
fn server_reports_its_identity_and_default_constructs() {
    use rmcp::ServerHandler;
    let from_default = PolynomialToolsServer::default();
    let info = from_default.get_info();
    assert_eq!(info.server_info.name, "polytool");
    assert_eq!(info.server_info.version, env!("CARGO_PKG_VERSION"));
    assert!(info.capabilities.tools.is_some());
    assert!(info
        .instructions
        .as_deref()
        .is_some_and(|text| text.contains("polytool")));
}

// -- planar-network tools ------------------------------------------------------

fn tnn_request() -> TotallyNonnegativeNetworkRequest {
    TotallyNonnegativeNetworkRequest {
        polynomials: None,
        text: None,
        matrix: None,
        chain: None,
    }
}

fn matrix_input(rows: &[&[i64]]) -> Option<Vec<Vec<BigIntCoefficientInput>>> {
    Some(
        rows.iter()
            .map(|row| {
                row.iter()
                    .map(|&v| BigIntCoefficientInput::Integer(v))
                    .collect()
            })
            .collect(),
    )
}

#[test]
fn tnn_network_matrix_input_and_failure_kinds() {
    let server = PolynomialToolsServer::new();
    let run = |request: TotallyNonnegativeNetworkRequest| {
        let Json(response) = server.check_tnn_network(Parameters(request)).unwrap();
        response
    };
    // Pascal's triangle (OEIS A007318) as a square zero-padded matrix and as triangular
    // rows: totally nonnegative, canonical network with all vertical weights 1.
    for rows in [
        &[&[1, 0, 0][..], &[1, 1, 0], &[1, 2, 1]][..],
        &[&[1][..], &[1, 1], &[1, 2, 1]][..],
    ] {
        let response = run(TotallyNonnegativeNetworkRequest {
            matrix: matrix_input(rows),
            ..tnn_request()
        });
        assert_eq!(response.input_mode, "matrix");
        assert_eq!(response.input_rows, 3);
        assert!(!response.prepended_p0);
        assert!(response.totally_nonnegative, "{:?}", response.failure);
        assert_eq!(
            response.coefficient_matrix,
            Some(vec![
                vec!["1".to_string()],
                vec!["1".to_string(), "1".to_string()],
                vec!["1".to_string(), "2".to_string(), "1".to_string()],
            ])
        );
        assert_eq!(response.certificate_verified, Some(true));
        let lambda = response.lambda.as_ref().unwrap();
        assert!(lambda.iter().flatten().all(|weight| weight == "1"));
        assert!(response.note.contains("does not imply"));
    }

    // Each structural failure names its kind and the offending entry.
    let failure = |rows: &[&[i64]]| {
        let response = run(TotallyNonnegativeNetworkRequest {
            matrix: matrix_input(rows),
            ..tnn_request()
        });
        assert!(!response.totally_nonnegative);
        assert!(response.lambda.is_none() && response.certificate_verified.is_none());
        response.failure.expect("failure")
    };
    let above = failure(&[&[1, 1], &[0, 1]]);
    assert_eq!(above.kind, "nonzero_above_diagonal");
    assert_eq!(
        (above.row, above.col, above.value.as_deref()),
        (Some(0), Some(1), Some("1"))
    );
    let diagonal = failure(&[&[2, 0], &[1, 1]]);
    assert_eq!(diagonal.kind, "non_unit_diagonal");
    assert_eq!(
        (diagonal.row, diagonal.col, diagonal.value.as_deref()),
        (Some(0), Some(0), Some("2"))
    );
    let negative = failure(&[&[1], &[-1, 1]]);
    assert_eq!(negative.kind, "negative_entry");
    assert_eq!(
        (negative.row, negative.col, negative.value.as_deref()),
        (Some(1), Some(0), Some("-1"))
    );
    let short_row = failure(&[&[1, 0], &[1]]);
    assert_eq!(short_row.kind, "invalid_row_length");
    assert_eq!(short_row.row, Some(0));
    // [[1],[0,1],[1,0,1]] has the negative 2 x 2 minor on rows {1, 2}, columns {0, 1}:
    // det [[0, 1], [1, 0]] = -1, reported as an explicit witness.
    let not_tnn = failure(&[&[1], &[0, 1], &[1, 0, 1]]);
    assert_eq!(not_tnn.kind, "not_tnn");
    let witness = not_tnn.witness.expect("explicit negative minor");
    assert_eq!(witness.value, "-1");
    assert_eq!(witness.rows.len(), witness.cols.len());

    // Polynomial input: the monic sequence P_1, P_2, ... gets P_0 = 1 prepended; the
    // chain polynomials of Pascal's rows are reported on request.
    let response = run(TotallyNonnegativeNetworkRequest {
        text: Some("1,1\n1,2,1\n1,3,3,1".to_string()),
        chain: Some(true),
        ..tnn_request()
    });
    assert_eq!(response.input_mode, "polynomials");
    assert!(response.prepended_p0 && response.totally_nonnegative);
    assert_eq!(response.coefficient_matrix.as_ref().unwrap().len(), 4);
    let chain = response.chain.expect("chain requested");
    assert!(chain.error.is_none());
    assert_eq!(chain.polynomials.len(), 4);
    assert_eq!(chain.interlacing.len(), 3);
    assert!(chain.contradicts_theorem == Some(false));
    // A non-monic polynomial is a structured failure, also for the requested chain.
    let response = run(TotallyNonnegativeNetworkRequest {
        text: Some("1,2".to_string()),
        chain: Some(true),
        ..tnn_request()
    });
    let failure = response.failure.expect("non-monic");
    assert_eq!(failure.kind, "non_monic_polynomial");
    assert_eq!(failure.value.as_deref(), Some("2"));
    assert!(response.chain.unwrap().error.is_some());
    let response = run(TotallyNonnegativeNetworkRequest {
        text: Some("1\n1,1,1".to_string()),
        ..tnn_request()
    });
    assert_eq!(
        response.failure.unwrap().kind,
        "unexpected_polynomial_degree"
    );
}

#[test]
fn network_tools_enforce_their_row_degree_and_search_limits() {
    use crate::limits::{
        MCP_MAX_BRENTI_INTEGER_NODES, MCP_MAX_BRENTI_ROWS, MCP_MAX_BRENTI_SEARCH_NODES,
        MCP_MAX_LEFT_PRODUCTION_ROWS, MCP_MAX_PF_NETWORK_DEGREE, MCP_MAX_PF_NETWORK_ROWS,
        MCP_MAX_PF_NETWORK_STEP, MCP_MAX_TNN_NETWORK_ROWS,
    };
    let server = PolynomialToolsServer::new();
    let rows = |count: usize, row: &str| vec![row; count].join("\n");

    // check_tnn_network: exactly one source, at least one row, at most 128 rows.
    let tnn_error = |request: TotallyNonnegativeNetworkRequest| {
        message(
            server
                .check_tnn_network(Parameters(request))
                .err()
                .expect("error"),
        )
    };
    assert!(tnn_error(tnn_request())
        .contains("expected exactly one of `polynomials`, `text` or `matrix`"));
    assert!(tnn_error(TotallyNonnegativeNetworkRequest {
        text: Some("1,1".to_string()),
        matrix: matrix_input(&[&[1]]),
        ..tnn_request()
    })
    .contains("expected exactly one"));
    assert!(tnn_error(TotallyNonnegativeNetworkRequest {
        polynomials: Some(Vec::new()),
        ..tnn_request()
    })
    .contains("at least one row is required"));
    assert!(tnn_error(TotallyNonnegativeNetworkRequest {
        text: Some(rows(MCP_MAX_TNN_NETWORK_ROWS + 1, "1,1")),
        ..tnn_request()
    })
    .contains("the MCP limit for check_tnn_network is 128"));
    assert!(tnn_error(TotallyNonnegativeNetworkRequest {
        text: Some("1,1\nbad + x".to_string()),
        ..tnn_request()
    })
    .contains("polynomial 1: "));

    // recognize_brenti_strip_digraph.
    let brenti = |text: String, nodes: Option<usize>, integer_nodes: Option<usize>| {
        server.recognize_brenti_strip_digraph(Parameters(BrentiRecognitionRequest {
            polynomials: None,
            text: Some(text),
            first_row_index: None,
            max_target_span: None,
            weights: None,
            max_search_nodes: nodes,
            max_integer_search_nodes: integer_nodes,
        }))
    };
    assert!(message(
        brenti(rows(MCP_MAX_BRENTI_ROWS + 1, "1,1"), None, None)
            .err()
            .expect("error")
    )
    .contains("the MCP limit for recognize_brenti_strip_digraph is 64"));
    assert!(message(
        brenti(
            "1,1".to_string(),
            None,
            Some(MCP_MAX_BRENTI_INTEGER_NODES + 1)
        )
        .err()
        .expect("error")
    )
    .contains("`max_integer_search_nodes` 100001 exceeds the MCP limit 100000"));
    assert!(message(
        brenti(
            "1,1".to_string(),
            Some(MCP_MAX_BRENTI_SEARCH_NODES + 1),
            None
        )
        .err()
        .expect("error")
    )
    .contains("`max_search_nodes` 1000001 exceeds the MCP limit 1000000"));

    // recognize_pf_network.
    let pf = |text: String,
              step: Option<usize>,
              rise: Option<usize>,
              nodes: Option<usize>,
              integer_nodes: Option<usize>| {
        server.recognize_pf_network(Parameters(PfNetworkRequest {
            polynomials: None,
            text: Some(text),
            first_row_index: None,
            max_step_length: step,
            max_rise: rise,
            weights: None,
            max_search_nodes: nodes,
            max_integer_search_nodes: integer_nodes,
        }))
    };
    let pf_error =
        |result: Result<Json<PfNetworkResponse>, McpError>| message(result.err().expect("error"));
    assert!(pf_error(pf(
        rows(MCP_MAX_PF_NETWORK_ROWS + 1, "1,1"),
        None,
        None,
        None,
        None
    ))
    .contains("the MCP limit for recognize_pf_network is 64"));
    assert!(pf_error(pf(
        vec!["1"; MCP_MAX_PF_NETWORK_DEGREE + 2].join(","),
        None,
        None,
        None,
        None
    ))
    .contains("polynomial 0 has degree 65"));
    assert!(pf_error(pf(
        "1,1".to_string(),
        Some(MCP_MAX_PF_NETWORK_STEP + 1),
        None,
        None,
        None
    ))
    .contains("`max_step_length` 9 exceeds the MCP limit 8"));
    assert!(pf_error(pf(
        "1,1".to_string(),
        None,
        Some(MCP_MAX_PF_NETWORK_STEP + 1),
        None,
        None
    ))
    .contains("`max_rise` 9 exceeds the MCP limit 8"));
    assert!(pf_error(pf(
        "1,1".to_string(),
        None,
        None,
        Some(MCP_MAX_BRENTI_SEARCH_NODES + 1),
        None
    ))
    .contains("`max_search_nodes` 1000001 exceeds"));
    assert!(pf_error(pf(
        "1,1".to_string(),
        None,
        None,
        None,
        Some(MCP_MAX_BRENTI_INTEGER_NODES + 1)
    ))
    .contains("`max_integer_search_nodes` 100001 exceeds"));
    // Binomial rows (OEIS A007318) are recognized by a height-invariant network with the
    // two single steps (length 1, rise 0 and 1), weight 1 each.
    let Json(pascal) = pf("1,1\n1,2,1\n1,3,3,1".to_string(), None, None, None, None).unwrap();
    assert!(pascal.certified, "{:?}", pascal.verdict);
    let network = pascal.network.expect("network");
    assert!(network
        .steps
        .iter()
        .flatten()
        .all(|step| step.length == 1 && step.weight == "1"));

    // check_left_production_matrix.
    let left = |text: String| {
        server.check_left_production_matrix(Parameters(LeftProductionRequest {
            polynomials: None,
            text: Some(text),
            first_row_index: None,
        }))
    };
    assert!(message(
        left(rows(MCP_MAX_LEFT_PRODUCTION_ROWS + 1, "1,1"))
            .err()
            .expect("error")
    )
    .contains("the MCP limit for check_left_production_matrix is 32"));
    assert!(
        message(left("1,1\nbad + x".to_string()).err().expect("error")).contains("polynomial 1: ")
    );
}

// -- error items and empty inputs ----------------------------------------------

#[test]
fn batch_tools_report_malformed_entries_in_place() {
    let server = PolynomialToolsServer::new();
    let text = || Some("1,2,1\nbad + x\n1,0,1".to_string());
    let batch = || BigIntPolynomialBatchInput {
        polynomials: None,
        text: text(),
    };

    let Json(properties) = server.polynomial_properties(Parameters(batch())).unwrap();
    assert!(properties.items[0].ok && properties.items[2].ok);
    assert!(!properties.items[1].ok);
    assert!(properties.items[1].error.is_some());
    assert!(properties.items[1].real_rooted.is_none());
    // (1 + t)^2 is real-rooted, 1 + t^2 is not.
    assert_eq!(properties.items[0].real_rooted, Some(true));
    assert_eq!(properties.items[2].real_rooted, Some(false));

    let Json(tests) = server.coefficient_tests(Parameters(batch())).unwrap();
    assert!(tests.items[0].ok && !tests.items[1].ok);
    assert!(tests.items[1].newton.is_none() && tests.items[1].kurtz.is_none());

    // Discriminants: disc((1 + t)^2) = 0 (double root), disc(1 + t^2) = -4.
    let Json(discriminants) = server.discriminant(Parameters(batch())).unwrap();
    assert_eq!(discriminants.items[0].discriminant.as_deref(), Some("0"));
    assert!(!discriminants.items[1].ok);
    assert!(discriminants.items[1].discriminant.is_none());
    assert_eq!(discriminants.items[2].discriminant.as_deref(), Some("-4"));

    let Json(hstar) = server
        .hstar_inequalities(Parameters(HStarInequalitiesRequest {
            polynomials: None,
            text: text(),
            dimension: None,
        }))
        .unwrap();
    assert!(hstar.items[0].ok && !hstar.items[1].ok);
    for dimension in [4097, 1 << 40, usize::MAX] {
        let refused = server.hstar_inequalities(Parameters(HStarInequalitiesRequest {
            polynomials: Some(vec![bigint_coeffs(&["1", "4", "1"])]),
            text: None,
            dimension: Some(dimension),
        }));
        assert!(
            message(refused.err().expect("dimension limit")).contains("at most 4096"),
            "{dimension}"
        );
    }
    assert!(hstar.items[1].error.is_some());
    assert!(hstar.items[1].checks.is_empty() && hstar.items[1].dimension.is_none());

    // Profile: a malformed entry has no polynomial, and the pairs around it carry its error.
    let Json(profile) = server
        .check_interlacing_profile(Parameters(batch()))
        .unwrap();
    assert!(!profile.profile[1].ok);
    assert_eq!(profile.profile[1].previous_count, 1);
    assert_eq!(profile.profile[1].checked_previous_count, 0);
    assert!(profile.profile[1].error.is_some());
    assert!(profile.profile[2].ok);
    assert_eq!(profile.profile[2].checked_previous_count, 1);
    assert_eq!(profile.profile[2].interlacing_previous_count, 0);
    assert!(!profile.profile[2].previous[0].ok);
}

#[test]
fn empty_coefficient_lists_normalize_to_the_zero_polynomial() {
    let server = PolynomialToolsServer::new();
    let Json(parsed) = server
        .parse_polynomials(Parameters(BigIntPolynomialBatchInput {
            polynomials: Some(vec![coeffs(&[]), coeffs(&[0, 0])]),
            text: None,
        }))
        .unwrap();
    for item in &parsed.items {
        assert_eq!(item.coefficients, Some(strings(&[0])));
        assert_eq!(item.degree, Some(0));
    }
    // Exact integers beyond i64 round-trip through the BigInt normalization, and
    // trailing zeros are removed there as well.
    let Json(properties) = server
        .polynomial_properties(Parameters(BigIntPolynomialBatchInput {
            polynomials: Some(vec![
                BigIntPolynomialInput {
                    coefficients: Some(Vec::new()),
                    expression: None,
                },
                bigint_coeffs(&["1", "0", "0"]),
            ]),
            text: None,
        }))
        .unwrap();
    let zero = properties.items[0].polynomial.as_ref().unwrap();
    assert_eq!(zero.coefficients, vec!["0"]);
    let one = properties.items[1].polynomial.as_ref().unwrap();
    assert_eq!(
        (one.coefficients.clone(), one.degree),
        (vec!["1".to_string()], 0)
    );
    // Rational recurrence rows drop trailing zeros too (the rows 1, 1, 2, 3 stay constant).
    let rows = parse_recurrence_batch_rational(&FindRecurrenceRequest {
        coefficients: Some(vec![
            Vec::new(),
            vec![
                RationalCoefficientInput::Integer(0),
                RationalCoefficientInput::Integer(0),
            ],
            vec![
                RationalCoefficientInput::Text("2".to_string()),
                RationalCoefficientInput::Integer(0),
            ],
        ]),
        ..Default::default()
    })
    .unwrap();
    assert_eq!(rows[0].as_ref().unwrap().len(), 1);
    assert_eq!(rows[1].as_ref().unwrap().len(), 1);
    assert_eq!(rows[2].as_ref().unwrap().len(), 1);
    assert_eq!(rows[2].as_ref().unwrap()[0].to_string(), "2");
}

#[test]
fn bigint_coefficient_json_accepts_integers_strings_and_rejects_other_types() {
    let parse = |value: serde_json::Value| {
        serde_json::from_value::<BigIntPolynomialInput>(json!({ "coefficients": [value] }))
    };
    let exact = "123456789012345678901234567890";
    let input = parse(json!(exact)).unwrap();
    let parsed = parse_required_bigint_polynomial(&input, "p").unwrap();
    assert_eq!(parsed[0].to_string(), exact);
    // Plain JSON integers keep their value; u64 values beyond i64 are kept exactly.
    let parsed = parse_required_bigint_polynomial(&parse(json!(-7)).unwrap(), "p").unwrap();
    assert_eq!(parsed[0], (-7).into());
    let parsed = parse_required_bigint_polynomial(&parse(json!(u64::MAX)).unwrap(), "p").unwrap();
    assert_eq!(parsed[0].to_string(), u64::MAX.to_string());
    // Booleans and nulls are rejected with a message naming the accepted forms.
    let error = parse(json!(true))
        .err()
        .expect("boolean rejected")
        .to_string();
    assert!(
        error.contains("an integer or an exact integer encoded as a decimal string"),
        "{error}"
    );
    assert!(parse(json!(null)).is_err());
}

#[test]
fn polynomial_properties_reports_exact_gamma_coefficients_for_rows_that_fit_in_i64() {
    // 1 + t^100 has coefficients 0 and 1, so it takes the i64 routes, but its
    // gamma coefficients (-1)^i 100/(100-i) C(100-i, i) reach 10^20 and exceed i64.
    // They used to be dropped (null), now they are exact.
    let server = PolynomialToolsServer::new();
    let mut row = vec!["0"; 101];
    row[0] = "1";
    row[100] = "1";
    let Json(response) = server
        .polynomial_properties(Parameters(BigIntPolynomialBatchInput {
            polynomials: Some(vec![bigint_coeffs(&row), bigint_coeffs(&["1", "4", "1"])]),
            text: None,
        }))
        .unwrap();
    let item = &response.items[0];
    assert_eq!(item.palindromic, Some(true));
    assert_eq!(item.gamma_positive, Some(false));
    let gamma = item
        .gamma_coefficients
        .as_ref()
        .expect("exact gamma coefficients");
    assert_eq!(gamma.len(), 51);
    assert_eq!(&gamma[..3], ["1", "-100", "4850"]);
    assert_eq!(gamma[50], "2");
    assert!(gamma
        .iter()
        .any(|value| value.trim_start_matches('-').len() >= 20));

    // Small rows are unchanged: 1 + 4t + t^2 = 1 (1+t)^2 + 2 t.
    let small = &response.items[1];
    assert_eq!(small.gamma_positive, Some(true));
    assert_eq!(
        small.gamma_coefficients,
        Some(vec!["1".to_string(), "2".to_string()])
    );
}

// -- coefficients beyond i64 ----------------------------------------------------

/// Coefficients of `(1 + scale x)^n`, `n = 0..=last`, as decimal strings.
fn scaled_binomial_strings(scale: u64, last: usize) -> Vec<Vec<String>> {
    let scale = num_bigint::BigInt::from(scale);
    (0..=last)
        .map(|n| {
            let mut binomial = num_bigint::BigInt::from(1);
            let mut power = num_bigint::BigInt::from(1);
            (0..=n)
                .map(|k| {
                    let value = (&binomial * &power).to_string();
                    binomial = &binomial * num_bigint::BigInt::from(n - k)
                        / num_bigint::BigInt::from(k + 1);
                    power *= &scale;
                    value
                })
                .collect()
        })
        .collect()
}

const TWO_TO_64: &str = "18446744073709551616";
const TEN_TO_20: &str = "100000000000000000000";

#[test]
fn parse_polynomials_reports_exact_coefficients_beyond_i64() {
    let server = PolynomialToolsServer::new();
    let Json(response) = server
        .parse_polynomials(Parameters(BigIntPolynomialBatchInput {
            polynomials: Some(vec![
                bigint_coeffs(&[TWO_TO_64, "-3", "0"]),
                coeffs(&[1, 2, 3]),
                BigIntPolynomialInput {
                    coefficients: None,
                    expression: Some(format!("{TEN_TO_20}t^2 + 1")),
                },
            ]),
            text: None,
        }))
        .unwrap();
    let first = &response.items[0];
    assert_eq!(
        first.coefficients,
        Some(vec![TWO_TO_64.to_string(), "-3".to_string()])
    );
    assert_eq!(first.degree, Some(1));
    assert_eq!(
        first.polynomial.as_deref(),
        Some(format!("{TWO_TO_64} - 3t").as_str())
    );
    // Small inputs keep their values (now as strings).
    assert_eq!(response.items[1].coefficients, Some(strings(&[1, 2, 3])));
    assert_eq!(
        response.items[2].coefficients,
        Some(vec![
            "1".to_string(),
            "0".to_string(),
            TEN_TO_20.to_string()
        ])
    );

    // Bare JSON integers above i64::MAX are read exactly; text input works too.
    let input: BigIntPolynomialBatchInput = serde_json::from_value(json!({
        "polynomials": [{ "coefficients": [u64::MAX, 1] }]
    }))
    .unwrap();
    let Json(response) = server.parse_polynomials(Parameters(input)).unwrap();
    assert_eq!(
        response.items[0].coefficients,
        Some(vec![u64::MAX.to_string(), "1".to_string()])
    );
    let Json(response) = server
        .parse_polynomials(Parameters(BigIntPolynomialBatchInput {
            polynomials: None,
            text: Some(format!("1,{TEN_TO_20}\nbad + x")),
        }))
        .unwrap();
    assert_eq!(
        response.items[0].coefficients,
        Some(vec!["1".to_string(), TEN_TO_20.to_string()])
    );
    assert!(!response.items[1].ok);
    let value = serde_json::to_value(&response.items[0]).unwrap();
    assert_eq!(value["coefficients"][1], TEN_TO_20);
}

#[test]
fn real_roots_accepts_coefficients_beyond_i64() {
    let server = PolynomialToolsServer::new();
    let one_plus_ten20 = format!("-{}1", &TEN_TO_20[..TEN_TO_20.len() - 1]);
    let Json(response) = server
        .real_roots(Parameters(BigIntPolynomialBatchInput {
            polynomials: Some(vec![
                // (t - 1)(t - 10^20) = 10^20 - (10^20 + 1) t + t^2.
                bigint_coeffs(&[TEN_TO_20, &one_plus_ten20, "1"]),
                // 1 + 10^20 t^2 has no real root.
                bigint_coeffs(&["1", "0", TEN_TO_20]),
                // 1 + 10^20 t has the exact root -10^-20.
                bigint_coeffs(&["1", TEN_TO_20]),
                coeffs(&[2, -3, 1]),
            ]),
            text: None,
        }))
        .unwrap();
    let tolerance: BigRational = "1/1000000".parse().unwrap();
    let close = |text: &str, expected: &str| {
        let root: BigRational = text.parse().unwrap();
        let expected: BigRational = expected.parse().unwrap();
        let difference = &root - &expected;
        let magnitude = if difference < BigRational::from_integer(0.into()) {
            -difference
        } else {
            difference
        };
        magnitude < tolerance
    };
    let big = &response.items[0];
    assert_eq!(big.real_rooted, Some(true));
    let roots = big.roots.as_ref().unwrap();
    assert_eq!(roots.len(), 2);
    assert!(
        close(&roots[0], "1") && close(&roots[1], TEN_TO_20),
        "{roots:?}"
    );
    assert_eq!(big.polynomial.as_ref().unwrap().coefficients[0], TEN_TO_20);
    assert_eq!(response.items[1].real_rooted, Some(false));
    assert_eq!(response.items[1].roots, None);
    assert_eq!(
        response.items[2].roots,
        Some(vec![format!("-1/{TEN_TO_20}")])
    );
    // Small inputs: roots 1 and 2.
    let small = response.items[3].roots.as_ref().unwrap();
    assert!(close(&small[0], "1") && close(&small[1], "2"));
}

#[test]
fn analyze_decomposition_accepts_coefficients_beyond_i64() {
    let server = PolynomialToolsServer::new();
    // 10^20 (1 + 4t + 2t^2): every coordinate of the small example scales by 10^20.
    let scaled = |values: &[i64]| -> Vec<String> {
        values
            .iter()
            .map(|value| format!("{value}{}", &TEN_TO_20[1..]))
            .collect()
    };
    let row = scaled(&[1, 4, 2]);
    let Json(response) = server
        .analyze_decomposition(Parameters(BigIntPolynomialBatchInput {
            polynomials: Some(vec![bigint_coeffs(
                &row.iter().map(String::as_str).collect::<Vec<_>>(),
            )]),
            text: None,
        }))
        .unwrap();
    let item = &response.items[0];
    assert!(item.ok, "{:?}", item.error);
    let report = item.report.as_ref().unwrap();
    assert_eq!(report.polynomial.coefficients, row);
    assert_eq!(report.polynomial.degree, 2);
    assert_eq!(report.reciprocal.coefficients, scaled(&[2, 4, 1]));
    assert_eq!(report.a.coefficients, scaled(&[1, 3, 1]));
    assert_eq!(report.b.coefficients, scaled(&[1, 1]));
    assert_eq!(report.f_polynomial.coefficients, scaled(&[1, 6, 7]));
    assert_eq!(report.r_transform_of_f.coefficients, scaled(&[2, 8, 7]));
    assert_eq!(report.magic.coordinates, scaled(&[1, 2, -1]));
    assert_eq!(report.magic.partial_sum_checks, vec![false, false]);
    assert!(report.a_real_rooted && report.b_real_rooted);
    assert_eq!(report.b_interlaces_a, Some(true));
}

#[test]
fn find_recurrence_accepts_every_input_form_beyond_i64() {
    let server = PolynomialToolsServer::new();
    // (1 + 10^10 x)^n for n = 0..=5: coefficients up to 10^50.
    let rows = scaled_binomial_strings(10_000_000_000, 5);
    assert!(rows[5][2].len() > 19);
    let search = |request: FindRecurrenceRequest| {
        let Json(response) = server.find_recurrence(Parameters(request)).unwrap();
        assert!(response.found, "{:?}", response.error);
        response.recurrence.unwrap()
    };
    let polynomials = serde_json::from_value::<FindRecurrenceRequest>(json!({
        "polynomials": rows.iter().map(|row| json!({ "coefficients": row })).collect::<Vec<_>>(),
        "include_code": false,
    }))
    .unwrap();
    let from_polynomials = search(polynomials);
    assert!(from_polynomials.contains("P(n-1)"), "{from_polynomials}");
    let text = rows
        .iter()
        .map(|row| row.join(","))
        .collect::<Vec<_>>()
        .join("\n");
    let by_text = search(FindRecurrenceRequest {
        text: Some(text),
        include_code: Some(false),
        ..Default::default()
    });
    assert_eq!(by_text, from_polynomials);
    let by_coefficients = search(
        serde_json::from_value(json!({ "coefficients": rows, "include_code": false })).unwrap(),
    );
    assert_eq!(by_coefficients, from_polynomials);
    let expressions = rows
        .iter()
        .map(|row| {
            row.iter()
                .enumerate()
                .map(|(power, value)| format!("{value}t^{power}"))
                .collect::<Vec<_>>()
                .join(" + ")
        })
        .collect();
    let by_expression = search(FindRecurrenceRequest {
        expressions: Some(expressions),
        include_code: Some(false),
        ..Default::default()
    });
    assert_eq!(by_expression, from_polynomials);

    // Bare JSON integers above i64::MAX are read exactly in `coefficients` and `polynomials`.
    let u64_rows = json!([[u64::MAX], [u64::MAX], [u64::MAX], [u64::MAX]]);
    let request: FindRecurrenceRequest = serde_json::from_value(json!({
        "coefficients": u64_rows, "include_code": false
    }))
    .unwrap();
    let rows = parse_recurrence_batch_rational(&request).unwrap();
    assert_eq!(
        rows[0].as_ref().unwrap()[0].to_string(),
        u64::MAX.to_string()
    );
    let request: FindRecurrenceRequest = serde_json::from_value(json!({
        "polynomials": [{ "coefficients": [u64::MAX, 1] }]
    }))
    .unwrap();
    let rows = parse_recurrence_batch_rational(&request).unwrap();
    assert_eq!(
        rows[0].as_ref().unwrap()[0].to_string(),
        u64::MAX.to_string()
    );
    // Floating-point JSON numbers are still refused, with the actionable message.
    let error = serde_json::from_value::<FindRecurrenceRequest>(json!({
        "coefficients": [[1.5e30]]
    }))
    .unwrap_err()
    .to_string();
    assert!(error.contains("floating-point JSON number"), "{error}");
}

#[test]
fn find_recurrence_representations_handles_big_rows_and_reports_row_errors() {
    let server = PolynomialToolsServer::new();
    let rows = scaled_binomial_strings(10_000_000_000, 4);
    let request: FindRecurrenceRequest = serde_json::from_value(json!({
        "polynomials": rows.iter().map(|row| json!({ "coefficients": row })).collect::<Vec<_>>(),
        "max_rec_len": 1, "max_diff_deg": 0, "max_idx_deg": 0, "max_var_deg": 1,
    }))
    .unwrap();
    let Json(report) = server
        .find_recurrence_representations(Parameters(request))
        .unwrap();
    assert_eq!(report["schema"], "polytool.recurrence-representations.v1");
    assert_eq!(report["status"], "complete");
    assert_eq!(report["representations"][0]["evidence"], "exact_fit");
    assert_eq!(
        report["representations"][0]["recurrence_json"]["initial_polynomials"][0][0],
        "1"
    );

    // Per-row parse errors are named in the error.
    let bad: FindRecurrenceRequest = serde_json::from_value(json!({
        "text": "1,1\nbad + x\n1,2,1\nalso bad\n1"
    }))
    .unwrap();
    let error = message(
        server
            .find_recurrence_representations(Parameters(bad))
            .err()
            .expect("expected an error"),
    );
    assert!(
        error.starts_with("one or more polynomials failed to parse"),
        "{error}"
    );
    assert!(
        error.contains("polynomial 1:") && error.contains("polynomial 3:"),
        "{error}"
    );
}

#[test]
fn check_polynomial_family_accepts_coefficients_beyond_i64() {
    let server = PolynomialToolsServer::new();
    // (1 + 10^10 x)^n for n = 1..=6: all rows are real-rooted; the last has coefficients up to 10^60.
    let rows = scaled_binomial_strings(10_000_000_000, 6);
    let request = CheckPolynomialFamilyRequest {
        polynomials: Some(
            rows[1..]
                .iter()
                .map(|row| bigint_coeffs(&row.iter().map(String::as_str).collect::<Vec<_>>()))
                .collect(),
        ),
        text: None,
        sequence: None,
        max_n: None,
        options: Some(FamilyCheckOptions {
            find_recurrence: Some(true),
            ..family_options()
        }),
        recurrence_options: None,
        lace: Some(LaceCheckRequest {
            block_rows: Some(2),
            block_cols: Some(4),
            max_minor_size: Some(2),
            include_matrix: Some(true),
        }),
    };
    let Json(response) = server.check_polynomial_family(Parameters(request)).unwrap();
    assert_eq!(response.source, "polynomials");
    assert_eq!(response.item_count, 6);
    assert!(response.parse_errors.is_empty());
    assert!(response.items.iter().all(|item| item.real_rooted));
    let last = response.items.last().unwrap();
    assert_eq!(last.polynomial.coefficients, rows[6]);
    assert_eq!(last.polynomial.degree, 6);
    assert!(!last.palindromic && !last.simple_roots);
    assert!(last.unimodal && last.log_concave);
    assert_eq!(response.consecutive_pairs.len(), 5);
    assert!(response
        .consecutive_pairs
        .iter()
        .all(|pair| pair.ok && pair.result.is_some()));
    let recurrence = response.recurrence.as_ref().expect("recurrence searched");
    assert!(recurrence.found, "{:?}", recurrence.error);
    let lace = response.lace.as_ref().expect("lace checked");
    assert_eq!((lace.rows, lace.columns), (12, 4));
    let matrix = lace.matrix.as_ref().expect("matrix requested");
    assert_eq!(matrix[0], vec!["1", "10000000000", "0", "0"]);
    assert!(matrix.iter().flatten().any(|value| value.len() > 19));
    // The wire format: strings, not numbers.
    let value = serde_json::to_value(&response).unwrap();
    assert_eq!(value["items"][5]["coefficients"][6], rows[6][6]);
    assert_eq!(value["lace"]["matrix"][0][1], "10000000000");
    assert!(
        response.markdown.contains(&rows[6][6])
            || response.markdown.contains("Polynomial family check")
    );

    // A palindromic row whose gamma coefficients exceed i64 reports them exactly.
    let mut palindrome = vec!["0"; 101];
    palindrome[0] = "1";
    palindrome[100] = "1";
    let Json(response) = server
        .check_polynomial_family(Parameters(CheckPolynomialFamilyRequest {
            polynomials: Some(vec![bigint_coeffs(&palindrome)]),
            options: Some(FamilyCheckOptions {
                check_consecutive_interlacing: Some(false),
                ..family_options()
            }),
            ..family_request(&[])
        }))
        .unwrap();
    let item = &response.items[0];
    assert!(item.palindromic && !item.gamma_positive);
    let gamma = item
        .gamma_coefficients
        .as_ref()
        .expect("exact gamma coefficients");
    assert_eq!(&gamma[..3], ["1", "-100", "4850"]);

    // Text input accepts big values too, and reports parse errors per item.
    let Json(response) = server
        .check_polynomial_family(Parameters(CheckPolynomialFamilyRequest {
            polynomials: None,
            text: Some(format!("1,{TWO_TO_64}\n1,{TEN_TO_20},{TEN_TO_20}")),
            ..family_request(&[])
        }))
        .unwrap();
    assert_eq!(response.source, "text");
    assert_eq!(response.items[0].polynomial.coefficients[1], TWO_TO_64);
}

#[test]
fn rational_coefficient_json_keeps_large_unsigned_integers() {
    let parsed: Vec<RationalCoefficientInput> =
        serde_json::from_value(json!([1, -2, u64::MAX, "3/7"])).unwrap();
    assert!(matches!(parsed[0], RationalCoefficientInput::Integer(1)));
    assert!(matches!(parsed[1], RationalCoefficientInput::Integer(-2)));
    assert!(
        matches!(&parsed[2], RationalCoefficientInput::Text(text) if text == &u64::MAX.to_string())
    );
    assert!(matches!(&parsed[3], RationalCoefficientInput::Text(text) if text == "3/7"));
    let error = serde_json::from_value::<RationalCoefficientInput>(json!(true))
        .unwrap_err()
        .to_string();
    assert!(error.contains("an integer or an exact rational"), "{error}");
}
