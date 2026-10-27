//! Tool-level tests: call the MCP tool methods directly and check the
//! structured responses and error reporting.

use super::*;
use crate::{
    family::check_family_lace,
    limits::{
        MCP_MAX_BATCH_POLYNOMIALS, MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL, MCP_MAX_LACE_MATRIX_CELLS,
        MCP_MAX_RECURRENCE_DEGREE, MCP_MAX_RECURRENCE_INPUT_ROWS, MCP_MAX_RECURRENCE_LENGTH,
        MCP_MAX_RECURRENCE_ROWS, MCP_MAX_SEQUENCE_N,
    },
    parsing::{
        normalize_polynomial, parse_batch, parse_polynomial_input, parse_recurrence_batch_rational,
    },
    recurrence_options::apply_recurrence_options,
    report_items::parse_items,
};
use rmcp::{handler::server::wrapper::Parameters, Json};
use serde_json::json;

fn coeffs(values: &[i64]) -> PolynomialInput {
    PolynomialInput {
        coefficients: Some(values.to_vec()),
        expression: None,
    }
}

fn expr(value: &str) -> PolynomialInput {
    PolynomialInput {
        coefficients: None,
        expression: Some(value.to_string()),
    }
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

#[test]
fn parses_coefficients_and_expression() {
    assert_eq!(
        normalize_polynomial(parse_polynomial_input(&coeffs(&[1, 2, 3, 0])).unwrap()),
        NormalizedPolynomial {
            polynomial: "1 + 2t + 3t^2".to_string(),
            coefficients: vec![1, 2, 3],
            degree: 2,
        }
    );
    assert_eq!(
        normalize_polynomial(parse_polynomial_input(&expr("3t^2 + 2t + 1")).unwrap()).coefficients,
        vec![1, 2, 3]
    );
}

#[test]
fn preserves_batch_parse_errors() {
    let batch = parse_batch(&PolynomialBatchInput {
        polynomials: None,
        text: Some("1,2,1\nbad + x\n".to_string()),
    })
    .unwrap();
    let items = parse_items(&batch);
    assert_eq!(items.len(), 2);
    assert!(items[0].ok);
    assert!(!items[1].ok);
}

#[test]
fn checks_eulerian_properties() {
    let server = PolynomialToolsServer::new();
    let Json(response) = server
        .polynomial_properties(Parameters(BigIntPolynomialBatchInput {
            polynomials: Some(vec![bigint_coeffs(&["1", "11", "11", "1"])]),
            text: None,
        }))
        .unwrap();
    let item = &response.items[0];
    assert_eq!(item.real_rooted, Some(true));
    assert_eq!(item.gamma_positive, Some(true));
    assert_eq!(
        item.gamma_coefficients,
        Some(vec!["1".to_string(), "8".to_string()])
    );
    assert_eq!(item.unimodal, Some(true));
}

#[test]
fn polynomial_properties_accepts_bigint_coefficients() {
    let server = PolynomialToolsServer::new();
    let large = "1000000000000000000000000000000000000000000";
    let Json(response) = server
        .polynomial_properties(Parameters(BigIntPolynomialBatchInput {
            polynomials: Some(vec![bigint_coeffs(&[large, "2", "1"])]),
            text: None,
        }))
        .unwrap();
    let item = &response.items[0];
    assert!(item.ok);
    assert_eq!(
        item.polynomial
            .as_ref()
            .and_then(|polynomial| polynomial.coefficients.first()),
        Some(&large.to_string())
    );
    let value = serde_json::to_value(response).unwrap();
    assert_eq!(value["items"][0]["polynomial"]["coefficients"][0], large);
}

#[test]
fn coefficient_tests_report_kurtz() {
    let server = PolynomialToolsServer::new();
    let Json(response) = server
        .coefficient_tests(Parameters(BigIntPolynomialBatchInput {
            polynomials: Some(vec![bigint_coeffs(&["1000", "1110", "111", "1"])]),
            text: None,
        }))
        .unwrap();
    let item = &response.items[0];
    assert!(item.ok);
    assert!(item.kurtz.as_ref().unwrap().holds);
    assert!(item.kurtz.as_ref().unwrap().implies_real_rooted);
}

#[test]
fn hstar_inequalities_report_named_failures() {
    let server = PolynomialToolsServer::new();
    let Json(response) = server
        .hstar_inequalities(Parameters(HStarInequalitiesRequest {
            polynomials: Some(vec![bigint_coeffs(&["1", "20", "1", "0"])]),
            text: None,
            dimension: Some(3),
        }))
        .unwrap();
    let item = &response.items[0];
    assert!(item.ok);
    assert_eq!(item.all_applicable_hold, Some(false));
    assert!(item.checks.iter().any(|check| {
        check.family == "Balletti-Higashitani" && check.applicable && !check.holds
    }));
}

#[test]
fn cyclic_sieving_checks_order_two() {
    let server = PolynomialToolsServer::new();
    let Json(response) = server
        .cyclic_sieving(Parameters(CyclicSievingRequest {
            polynomial: bigint_coeffs(&["1", "1"]),
            order: 2,
            fixed_counts: Some(vec![1.into(), 0.into()]),
            orbit_sizes: None,
        }))
        .unwrap();
    assert_eq!(response.report.holds, Some(false));

    let Json(response) = server
        .cyclic_sieving(Parameters(CyclicSievingRequest {
            polynomial: bigint_coeffs(&["1", "1"]),
            order: 2,
            fixed_counts: Some(vec![2.into(), 0.into()]),
            orbit_sizes: None,
        }))
        .unwrap();
    assert_eq!(response.report.holds, Some(true));
}

fn cyclic_request(
    coefficients: &[&str],
    order: usize,
    orbit_sizes: Option<Vec<usize>>,
) -> CyclicSievingRequest {
    CyclicSievingRequest {
        polynomial: bigint_coeffs(coefficients),
        order,
        fixed_counts: None,
        orbit_sizes,
    }
}

#[test]
fn cyclic_sieving_reports_orbit_decomposition_and_cone_condition() {
    let server = PolynomialToolsServer::new();
    // Alexandersson-Amini counterexample: integral, natural, but S_3 = -3.
    let Json(response) = server
        .cyclic_sieving(Parameters(cyclic_request(
            &["10", "1", "0", "3", "0", "1"],
            6,
            None,
        )))
        .unwrap();
    let dec = &response.orbit_decomposition;
    assert!(dec.integral && dec.natural_coefficients);
    assert_eq!(dec.cone_condition, Some(false));
    assert!(!dec.realizable);
    assert_eq!(
        dec.fixed_points.as_ref().unwrap(),
        &["15", "8", "12", "5", "12", "8"]
    );
    let size3 = dec
        .orbits
        .as_ref()
        .unwrap()
        .iter()
        .find(|entry| entry.orbit_size == 3)
        .unwrap();
    assert_eq!(size3.orbit_count, "-1");
    assert_eq!(size3.element_count, "-3");
    assert!(response.orbit_sizes_check.is_none());

    // non-integral
    let Json(response) = server
        .cyclic_sieving(Parameters(cyclic_request(&["1", "1"], 3, None)))
        .unwrap();
    assert!(!response.orbit_decomposition.integral);
    assert_eq!(response.orbit_decomposition.violating_classes, vec![1]);
    assert!(response.orbit_decomposition.orbits.is_none());
    assert!(!response.orbit_decomposition.realizable);
}

#[test]
fn cyclic_sieving_orbit_sizes_alternative_to_fixed_counts() {
    let server = PolynomialToolsServer::new();
    let coefficients = ["1", "1", "2", "1", "1"];
    let Json(response) = server
        .cyclic_sieving(Parameters(cyclic_request(
            &coefficients,
            4,
            Some(vec![4, 2]),
        )))
        .unwrap();
    assert_eq!(response.report.holds, Some(true));
    assert!(response.orbit_sizes_check.as_ref().unwrap().holds);
    let Json(response) = server
        .cyclic_sieving(Parameters(cyclic_request(
            &coefficients,
            4,
            Some(vec![4, 1, 1]),
        )))
        .unwrap();
    assert_eq!(response.report.holds, Some(false));
    let check = response.orbit_sizes_check.unwrap();
    assert!(!check.holds);
    assert_eq!(check.mismatched_orbit_sizes, vec![1, 2]);

    let mut both = cyclic_request(&coefficients, 4, Some(vec![4, 2]));
    both.fixed_counts = Some(vec![6.into(), 0.into(), 2.into(), 0.into()]);
    assert!(server.cyclic_sieving(Parameters(both)).is_err());
    assert!(server
        .cyclic_sieving(Parameters(cyclic_request(&coefficients, 4, Some(vec![3]))))
        .is_err());
}

#[test]
fn bicyclic_sieving_one_plus_qt() {
    let server = PolynomialToolsServer::new();
    let matrix = |rows: &[&[i64]]| -> Vec<Vec<BigIntCoefficientInput>> {
        rows.iter()
            .map(|row| {
                row.iter()
                    .map(|&v| BigIntCoefficientInput::Integer(v))
                    .collect()
            })
            .collect()
    };
    let Json(response) = server
        .bicyclic_sieving(Parameters(BicyclicSievingRequest {
            matrix: matrix(&[&[1, 0], &[0, 1]]),
            order1: 2,
            order2: 2,
            fixed_table: Some(matrix(&[&[2, 0], &[0, 2]])),
        }))
        .unwrap();
    assert_eq!(response.holds, Some(true));
    assert!(response.integral && response.product_basis_form);
    assert_eq!(
        response.product_coordinates.unwrap(),
        vec![vec!["1", "-1"], vec!["-1", "2"]]
    );
    let Json(response) = server
        .bicyclic_sieving(Parameters(BicyclicSievingRequest {
            matrix: matrix(&[&[1, 0], &[0, 1]]),
            order1: 2,
            order2: 2,
            fixed_table: Some(matrix(&[&[2, 0], &[0, 0]])),
        }))
        .unwrap();
    assert_eq!(response.holds, Some(false));
    assert_eq!(response.mismatches, vec![[1, 1]]);
    assert!(server
        .bicyclic_sieving(Parameters(BicyclicSievingRequest {
            matrix: matrix(&[&[1]]),
            order1: 0,
            order2: 2,
            fixed_table: None,
        }))
        .is_err());
}

#[test]
fn cyclic_sieving_sequence_profiles_default_offsets() {
    let server = PolynomialToolsServer::new();
    let Json(response) = server
        .cyclic_sieving_sequence(Parameters(CyclicSievingSequenceRequest {
            polynomials: Some(vec![
                bigint_coeffs(&["1", "1"]),
                bigint_coeffs(&["1", "1", "1"]),
            ]),
            text: None,
            first_index: Some(2),
            offsets: None,
            fixed_counts: None,
        }))
        .unwrap();
    assert_eq!(response.items[0].index, 2);
    assert!(response.items[0]
        .candidate_orders
        .iter()
        .any(|report| report.order == 2));
    assert!(response.items[0]
        .candidate_orders
        .iter()
        .any(|report| report.order == 5));
}

#[test]
fn checks_strict_and_weak_interlacing() {
    let server = PolynomialToolsServer::new();
    let Json(strict) = server
        .check_interlacing_pair(Parameters(BigIntPolynomialPairRequest {
            p: bigint_coeffs(&["-2", "1"]),
            q: bigint_coeffs(&["3", "-4", "1"]),
        }))
        .unwrap();
    assert_eq!(strict.strict, Some(true));

    let Json(weak) = server
        .check_interlacing_pair(Parameters(BigIntPolynomialPairRequest {
            p: bigint_coeffs(&["-1", "1"]),
            q: bigint_coeffs(&["2", "-3", "1"]),
        }))
        .unwrap();
    assert_eq!(weak.weak, Some(true));
}

#[test]
fn interlacing_accepts_bigint_coefficients() {
    let server = PolynomialToolsServer::new();
    let center = "100000000000000000000";
    let Json(result) = server
        .check_interlacing_pair(Parameters(BigIntPolynomialPairRequest {
            p: bigint_coeffs(&["-100000000000000000000", "1"]),
            q: bigint_coeffs(&[
                "9999999999999999999999999999999999999999",
                "-200000000000000000000",
                "1",
            ]),
        }))
        .unwrap();

    assert_eq!(result.strict, Some(true));
    assert_eq!(result.p.coefficients[0], format!("-{center}"));
    assert_eq!(result.q.degree, 2);
}

#[test]
fn interlacing_profile_stops_at_first_failure() {
    let server = PolynomialToolsServer::new();
    let Json(response) = server
        .check_interlacing_profile(Parameters(BigIntPolynomialBatchInput {
            polynomials: Some(vec![
                bigint_coeffs(&["1"]),
                bigint_coeffs(&["1", "1"]),
                bigint_coeffs(&["1", "2", "1"]),
                bigint_coeffs(&["1", "0", "1"]),
                bigint_coeffs(&["1", "3", "3", "1"]),
            ]),
            text: None,
        }))
        .unwrap();

    assert_eq!(response.profile[2].previous_count, 2);
    assert_eq!(response.profile[2].checked_previous_count, 2);
    assert_eq!(response.profile[2].interlacing_previous_count, 1);
    assert_eq!(response.profile[2].previous_interlacing_indices, vec![1]);

    assert_eq!(response.profile[3].previous_count, 3);
    assert_eq!(response.profile[3].checked_previous_count, 1);
    assert_eq!(response.profile[3].interlacing_previous_count, 0);
}

#[test]
fn generate_sequence_returns_bigint_coefficients() {
    let server = PolynomialToolsServer::new();
    let Json(response) = server
        .generate_sequence(Parameters(GenerateSequenceRequest {
            sequence: SequenceKind::ChebyshevT,
            max_n: 64,
        }))
        .unwrap();

    let last = response.polynomials.last().expect("T_64");
    assert_eq!(last.degree, 64);
    assert_eq!(
        last.coefficients.last().map(String::as_str),
        Some("9223372036854775808")
    );
}

#[test]
fn lists_and_generates_oeis_catalog_rows() {
    let server = PolynomialToolsServer::new();
    let Json(list) = server
        .list_oeis_sequences(Parameters(ListOeisSequencesRequest::default()))
        .unwrap();
    assert_eq!(list.count, 774);
    assert!(list.sequences.iter().any(|entry| entry.id == "A008292"));
    assert!(list
        .sequences
        .iter()
        .any(|entry| entry.id == "A390883" && entry.status == "validated"));

    let Json(response) = server
        .generate_oeis_rows(Parameters(GenerateOeisRowsRequest {
            id: "a008292".to_string(),
            rows: 3,
            first_row: None,
            include_experimental: None,
        }))
        .unwrap();
    assert_eq!(response.first_row, 1);
    assert_eq!(response.rows[2].n, 3);
    assert_eq!(response.rows[2].coefficients, vec!["1", "4", "1"]);

    let Json(export) = server
        .get_oeis_sequence(Parameters(GetOeisSequenceRequest {
            id: "A166345".to_string(),
        }))
        .unwrap();
    assert_eq!(
        export.recurrence,
        "P(n) = (1 - t + nt) P(n-1) + (t - t^2) P'(n-1)"
    );
    assert!(export.python.contains("    3: [1, 2, 1],"));
}

#[test]
fn computes_resultant_and_discriminant() {
    let server = PolynomialToolsServer::new();
    let Json(resultant) = server
        .resultant(Parameters(BigIntPolynomialPairRequest {
            p: bigint_coeffs(&["2", "-3", "1"]),
            q: bigint_coeffs(&["-3", "1"]),
        }))
        .unwrap();
    assert_eq!(resultant.resultant, "2");

    let Json(discriminant) = server
        .discriminant(Parameters(BigIntPolynomialBatchInput {
            polynomials: Some(vec![bigint_coeffs(&["-1", "0", "1"])]),
            text: None,
        }))
        .unwrap();
    assert_eq!(discriminant.items[0].discriminant.as_deref(), Some("4"));

    let large = "1000000000000000000000000000000";
    let Json(large_discriminant) = server
        .discriminant(Parameters(BigIntPolynomialBatchInput {
            polynomials: Some(vec![bigint_coeffs(&["1", "0", large])]),
            text: None,
        }))
        .unwrap();
    assert_eq!(
        large_discriminant.items[0].discriminant.as_deref(),
        Some("-4000000000000000000000000000000")
    );
}

#[test]
fn roundtrips_hstar_and_ehrhart() {
    let server = PolynomialToolsServer::new();
    let Json(to_ehrhart) = server
        .ehrhart_hstar(Parameters(EhrhartHstarRequest {
            mode: EhrhartHstarMode::HstarToEhrhart,
            hstar: Some(vec![1.into(), 1.into(), 0.into()]),
            ehrhart_coefficients: None,
            numerator_coefficients: None,
            denominator: None,
        }))
        .unwrap();
    assert_eq!(
        to_ehrhart.ehrhart_coefficients.as_ref().unwrap(),
        &vec!["1".to_string(), "2".to_string(), "1".to_string()]
    );

    let Json(to_hstar) = server
        .ehrhart_hstar(Parameters(EhrhartHstarRequest {
            mode: EhrhartHstarMode::EhrhartToHstar,
            hstar: None,
            ehrhart_coefficients: Some(vec!["1".to_string(), "2".to_string(), "1".to_string()]),
            numerator_coefficients: None,
            denominator: None,
        }))
        .unwrap();
    assert_eq!(
        to_hstar.hstar,
        Some(vec!["1".to_string(), "1".to_string(), "0".to_string()])
    );

    let Json(common_denominator) = server
        .ehrhart_hstar(Parameters(EhrhartHstarRequest {
            mode: EhrhartHstarMode::EhrhartToHstar,
            hstar: None,
            ehrhart_coefficients: None,
            numerator_coefficients: Some(vec![2.into(), 3.into(), 1.into()]),
            denominator: Some(2.into()),
        }))
        .unwrap();
    assert_eq!(
        common_denominator.hstar,
        Some(vec!["1".to_string(), "0".to_string(), "0".to_string()])
    );
}

#[test]
fn find_recurrence_schema_is_flat_and_documents_every_control() {
    let schema = serde_json::to_value(schemars::schema_for!(FindRecurrenceRequest)).unwrap();
    assert_eq!(schema["type"], "object");
    assert!(schema.get("oneOf").is_none());
    assert_eq!(schema["additionalProperties"], false);
    let properties = schema["properties"].as_object().expect("schema properties");
    for name in [
        "polynomials",
        "coefficients",
        "expressions",
        "text",
        "options",
        "include_code",
        "skip_prefix",
        "min_rec_len",
        "max_rec_len",
        "min_var_deg",
        "max_var_deg",
        "min_idx_deg",
        "max_idx_deg",
        "min_diff_deg",
        "max_diff_deg",
        "try_inhomogeneous",
        "min_inhomo_var_deg",
        "max_inhomo_var_deg",
        "min_inhomo_idx_deg",
        "max_inhomo_idx_deg",
        "try_denominator",
        "try_alternating_sign",
        "max_denom_var_deg",
        "max_denom_idx_deg",
        "min_margin",
        "no_verify",
        "fit_extra_rows",
        "modular_prefilter",
        "max_candidates",
    ] {
        let property = properties
            .get(name)
            .unwrap_or_else(|| panic!("missing {name}"));
        assert!(
            property["description"]
                .as_str()
                .is_some_and(|text| !text.is_empty()),
            "{name} lacks a description: {property}"
        );
    }
}

#[test]
fn flat_and_legacy_recurrence_controls_work_with_top_level_precedence() {
    let server = PolynomialToolsServer::new();
    let flat: FindRecurrenceRequest = serde_json::from_value(json!({
        "coefficients": [[1], [2], [4], [8], [16]],
        "min_rec_len": 1,
        "max_rec_len": 1,
        "max_var_deg": 0,
        "max_idx_deg": 0,
        "max_diff_deg": 0,
        "max_candidates": 1
    }))
    .unwrap();
    let Json(flat_result) = server.find_recurrence(Parameters(flat)).unwrap();
    assert_eq!(flat_result.status, RecurrenceSearchStatus::Found);
    assert_eq!(flat_result.candidates_considered, Some(1));

    let legacy: FindRecurrenceRequest = serde_json::from_value(json!({
        "coefficients": [[1], [2], [4], [8], [16]],
        "options": {
            "min_rec_len": 1,
            "max_rec_len": 1,
            "max_var_deg": 0,
            "max_idx_deg": 0,
            "max_diff_deg": 0,
            "max_candidates": 1
        }
    }))
    .unwrap();
    let Json(legacy_result) = server.find_recurrence(Parameters(legacy)).unwrap();
    assert_eq!(legacy_result.status, RecurrenceSearchStatus::Found);

    let conflict: FindRecurrenceRequest = serde_json::from_value(json!({
        "coefficients": [[1], [2], [4], [8], [16]],
        "max_candidates": 1,
        "max_rec_len": 1,
        "options": { "max_candidates": 0, "max_rec_len": 7 }
    }))
    .unwrap();
    let effective = conflict.effective_options().unwrap();
    assert_eq!(effective.max_candidates, Some(1));
    assert_eq!(effective.max_rec_len, Some(1));
    let Json(conflict_result) = server.find_recurrence(Parameters(conflict)).unwrap();
    assert_eq!(conflict_result.status, RecurrenceSearchStatus::Found);
    assert_eq!(conflict_result.candidates_considered, Some(1));
}

#[test]
fn compact_recurrence_result_omits_code_and_default_result_remains_full() {
    let server = PolynomialToolsServer::new();
    let request = || {
        json!({
            "coefficients": [[1], [2], [4], [8], [16]],
            "max_rec_len": 1,
            "max_var_deg": 0,
            "max_idx_deg": 0,
            "max_diff_deg": 0
        })
    };

    let full: FindRecurrenceRequest = serde_json::from_value(request()).unwrap();
    let Json(full) = server.find_recurrence(Parameters(full)).unwrap();
    assert!(full.recurrence.is_some());
    assert!(full.latex.is_some());
    assert!(full.mathematica.is_some());
    assert!(full.python.is_some());
    assert!(full.sage.is_some());
    assert!(full.recurrence_json.is_some());

    let mut compact_request = request();
    compact_request["include_code"] = json!(false);
    let compact: FindRecurrenceRequest = serde_json::from_value(compact_request).unwrap();
    let Json(compact) = server.find_recurrence(Parameters(compact)).unwrap();
    assert!(compact.recurrence.is_some());
    assert!(compact.latex.is_some());
    assert!(compact.unknowns.is_some());
    assert!(compact.candidates_considered.is_some());
    assert!(compact.mathematica.is_none());
    assert!(compact.python.is_none());
    assert!(compact.sage.is_none());
    assert!(compact.recurrence_json.is_none());
    let serialized = serde_json::to_value(compact).unwrap();
    for absent in ["mathematica", "python", "sage", "recurrence_json"] {
        assert!(serialized.get(absent).is_none(), "unexpected {absent}");
    }
}

#[test]
fn recurrence_spaces_tool_certifies_and_replays() {
    let server = PolynomialToolsServer::new();
    let request = json!({"source":{"coefficients":[["2"]],"seeds":[["9007199254740993"]]},"shapes":[{"kind":"box","length":1}]});
    let Json(report) = server
        .find_recurrence_spaces(Parameters(RecurrenceSpacesToolRequest {
            request: request.clone(),
        }))
        .unwrap();
    assert_eq!(report["schema"], "polytool.recurrence-spaces.v1");
    let Json(checked) = server
        .verify_recurrence_spaces(Parameters(VerifyRecurrenceSpacesToolRequest {
            report: json!(report),
            expected_request: Some(request.clone()),
        }))
        .unwrap();
    assert_eq!(checked["verified_completed_claims"], true);
    let mut wrong = request.clone();
    wrong["domain_start"] = json!(1);
    assert!(server
        .verify_recurrence_spaces(Parameters(VerifyRecurrenceSpacesToolRequest {
            report: json!(report),
            expected_request: Some(wrong)
        }))
        .is_err());
    let mut excessive = request;
    excessive["limits"] = json!({"max_columns":1000000});
    assert!(server
        .find_recurrence_spaces(Parameters(RecurrenceSpacesToolRequest {
            request: excessive
        }))
        .is_err());
}

#[test]
fn recurrence_spaces_indexed_automatic_and_leader_coverage() {
    let server = PolynomialToolsServer::new();
    let request = json!({"indexed_source":{"coefficients":[[["0","1"],["0","1"]]],"seeds":[["1"]]},
        "include_leader_coverage":true,"shapes":[{"kind":"automatic_variable_degree","length":1,"index_degree":1}]});
    let Json(report) = server
        .find_recurrence_spaces(Parameters(RecurrenceSpacesToolRequest {
            request: request.clone(),
        }))
        .unwrap();
    assert!(report["results"][0]["propagation"].is_object());
    assert!(report["results"][0]["automatic_degree_bound"].is_object());
    assert!(report["results"][0]["leader_coverage"].is_object());
    let Json(checked) = server
        .verify_recurrence_spaces(Parameters(VerifyRecurrenceSpacesToolRequest {
            report: json!(report),
            expected_request: Some(request),
        }))
        .unwrap();
    assert_eq!(checked["verified_completed_claims"], true);
}

#[test]
fn recurrence_spaces_factorial_rescaling_certifies_and_binds_power() {
    let server = PolynomialToolsServer::new();
    let request = json!({"source":{"coefficients":[["1","1"]],"seeds":[["1"]]},
        "coefficient_rescaling":{"kind":"factorial_power","power":1},
        "shapes":[{"kind":"box","length":2,"index_degree":1,"variable_degree":1,"leading_index_degree":1}]});
    let Json(report) = server
        .find_recurrence_spaces(Parameters(RecurrenceSpacesToolRequest {
            request: request.clone(),
        }))
        .unwrap();
    let result = &report["results"][0];
    assert_eq!(result["status"], "successor_exists");
    assert_eq!(result["rescaling_bound"]["annihilator_power"], 3);
    assert_eq!(
        result["representatives"][0]["formula"],
        "(2 + n) Q[n+2] = (-1 - n) Q[n+0] + (3 + 2n + X) Q[n+1]"
    );
    let Json(checked) = server
        .verify_recurrence_spaces(Parameters(VerifyRecurrenceSpacesToolRequest {
            report: json!(report),
            expected_request: Some(request.clone()),
        }))
        .unwrap();
    assert_eq!(checked["verified_completed_claims"], true);
    let mut other = request;
    other["coefficient_rescaling"]["power"] = json!(2);
    assert!(server
        .verify_recurrence_spaces(Parameters(VerifyRecurrenceSpacesToolRequest {
            report: json!(report),
            expected_request: Some(other),
        }))
        .is_err());
}

#[test]
fn recurrence_representations_tool_uses_shared_schema_and_global_budget() {
    let server = PolynomialToolsServer::new();
    for limit in [0, 1] {
        let input: FindRecurrenceRequest = serde_json::from_value(json!({
            "coefficients": [["9007199254740993"], ["18014398509481986"],
                ["36028797018963972"], ["72057594037927944"], ["144115188075855888"]],
            "max_rec_len": 1, "max_diff_deg": 0, "max_idx_deg": 0,
            "max_var_deg": 0, "max_candidates": limit
        }))
        .unwrap();
        let Json(report) = server
            .find_recurrence_representations(Parameters(input))
            .unwrap();
        assert_eq!(report["schema"], "polytool.recurrence-representations.v1");
        assert_eq!(report["configurations_considered"], limit);
        if limit == 0 {
            assert_eq!(report["status"], "budget_exhausted");
            assert_eq!(report["coverage"][0]["status"], "incomplete");
        } else {
            assert_eq!(report["status"], "complete");
            assert_eq!(report["representations"][0]["evidence"], "exact_fit");
            assert_eq!(
                report["representations"][0]["recurrence_json"]["initial_polynomials"][0][0],
                "9007199254740993"
            );
        }
    }
}

#[test]
fn recurrence_input_exclusivity_remains_runtime_validated() {
    let empty = FindRecurrenceRequest::default();
    assert!(
        format!("{:?}", parse_recurrence_batch_rational(&empty).unwrap_err())
            .contains("expected exactly one")
    );

    let conflict: FindRecurrenceRequest = serde_json::from_value(json!({
        "coefficients": [[1], [2], [4]],
        "expressions": ["1", "2", "4"]
    }))
    .unwrap();
    assert!(format!(
        "{:?}",
        parse_recurrence_batch_rational(&conflict).unwrap_err()
    )
    .contains("expected exactly one"));
}

#[test]
fn recurrence_search_examples() {
    let server = PolynomialToolsServer::new();
    let Json(geometric) = server
        .find_recurrence(Parameters(FindRecurrenceRequest {
            polynomials: Some(vec![
                coeffs(&[1]),
                coeffs(&[2]),
                coeffs(&[4]),
                coeffs(&[8]),
                coeffs(&[16]),
            ]),
            coefficients: None,
            expressions: None,
            text: None,
            options: None,
            ..Default::default()
        }))
        .unwrap();
    assert_eq!(geometric.recurrence.as_deref(), Some("P(n) = 2 P(n-1)"));
    assert!(geometric
        .python
        .as_deref()
        .is_some_and(|python| python.contains("def P(n):")));
    let recurrence_json = geometric
        .recurrence_json
        .clone()
        .expect("recurrence JSON output");
    let Json(generated) = server
        .generate_recurrence_rows(Parameters(GenerateRecurrenceRowsRequest {
            recurrence_json,
            rows: Some(6),
            additional: None,
        }))
        .unwrap();
    assert_eq!(
        generated.polynomials,
        vec![
            vec!["1".to_string()],
            vec!["2".to_string()],
            vec!["4".to_string()],
            vec!["8".to_string()],
            vec!["16".to_string()],
            vec!["32".to_string()],
        ]
    );

    let Json(fibonacci) = server
        .find_recurrence(Parameters(FindRecurrenceRequest {
            polynomials: None,
            coefficients: Some(vec![
                vec![1.into()],
                vec![1.into()],
                vec![2.into()],
                vec![3.into()],
                vec![5.into()],
                vec![8.into()],
            ]),
            expressions: None,
            text: None,
            options: None,
            ..Default::default()
        }))
        .unwrap();
    assert_eq!(
        fibonacci.recurrence.as_deref(),
        Some("P(n) = P(n-1) + P(n-2)")
    );

    let Json(expression_geometric) = server
        .find_recurrence(Parameters(FindRecurrenceRequest {
            polynomials: None,
            coefficients: None,
            expressions: Some(vec![
                "1".to_string(),
                "2".to_string(),
                "4".to_string(),
                "8".to_string(),
                "16".to_string(),
            ]),
            text: None,
            options: None,
            ..Default::default()
        }))
        .unwrap();
    assert_eq!(
        expression_geometric.recurrence.as_deref(),
        Some("P(n) = 2 P(n-1)")
    );

    let Json(rational_geometric) = server
        .find_recurrence(Parameters(FindRecurrenceRequest {
            polynomials: None,
            coefficients: Some(vec![
                vec![RationalCoefficientInput::Text("1/2".to_string())],
                vec![RationalCoefficientInput::Text("1/4".to_string())],
                vec![RationalCoefficientInput::Text("1/8".to_string())],
                vec![RationalCoefficientInput::Text("1/16".to_string())],
                vec![RationalCoefficientInput::Text("1/32".to_string())],
            ]),
            expressions: None,
            text: None,
            options: None,
            ..Default::default()
        }))
        .unwrap();
    assert_eq!(
        rational_geometric.recurrence.as_deref(),
        Some("P(n) = 1/2 P(n-1)")
    );

    let large = "1267650600228229401496703205376";
    let double_large = "2535301200456458802993406410752";
    let Json(large_geometric) = server
        .find_recurrence(Parameters(FindRecurrenceRequest {
            polynomials: None,
            coefficients: Some(vec![
                vec![RationalCoefficientInput::Text(large.to_string())],
                vec![RationalCoefficientInput::Text(double_large.to_string())],
                vec![RationalCoefficientInput::Text(
                    "5070602400912917605986812821504".to_string(),
                )],
                vec![RationalCoefficientInput::Text(
                    "10141204801825835211973625643008".to_string(),
                )],
                vec![RationalCoefficientInput::Text(
                    "20282409603651670423947251286016".to_string(),
                )],
            ]),
            expressions: None,
            text: None,
            options: None,
            ..Default::default()
        }))
        .unwrap();
    assert_eq!(
        large_geometric.recurrence.as_deref(),
        Some("P(n) = 2 P(n-1)")
    );
    assert!(large_geometric
        .sage
        .as_deref()
        .is_some_and(|sage| sage.contains(large)));

    let Json(alternating) = server
        .find_recurrence(Parameters(FindRecurrenceRequest {
            polynomials: None,
            coefficients: Some(vec![
                vec![1.into()],
                vec![1.into()],
                vec![(-1).into()],
                vec![(-1).into()],
                vec![1.into()],
                vec![1.into()],
                vec![(-1).into()],
                vec![(-1).into()],
            ]),
            expressions: None,
            text: None,
            options: Some(
                serde_json::from_value(json!({
                    "try_alternating_sign": true,
                    "max_rec_len": 1,
                    "max_var_deg": 0,
                    "max_idx_deg": 0,
                    "max_diff_deg": 0
                }))
                .unwrap(),
            ),
            ..Default::default()
        }))
        .unwrap();
    assert_eq!(
        alternating.recurrence.as_deref(),
        Some("P(n) = (-1)^n P(n-1)")
    );
    assert!(alternating
        .sage
        .as_deref()
        .is_some_and(|sage| sage.contains("(-1)**n*P(n - 1)")));

    let eulerian = vec![
        coeffs(&[1]),
        coeffs(&[1]),
        coeffs(&[1, 1]),
        coeffs(&[1, 4, 1]),
        coeffs(&[1, 11, 11, 1]),
        coeffs(&[1, 26, 66, 26, 1]),
    ];
    let Json(eulerian_result) = server
        .find_recurrence(Parameters(FindRecurrenceRequest {
            polynomials: Some(eulerian),
            coefficients: None,
            expressions: None,
            text: None,
            options: None,
            ..Default::default()
        }))
        .unwrap();
    assert!(eulerian_result.found);

    let conflict = parse_recurrence_batch_rational(&FindRecurrenceRequest {
        polynomials: Some(vec![coeffs(&[1]), coeffs(&[2]), coeffs(&[4])]),
        coefficients: Some(vec![vec![1.into()], vec![2.into()], vec![4.into()]]),
        expressions: None,
        text: None,
        options: None,
        ..Default::default()
    })
    .unwrap_err();
    assert!(format!("{conflict:?}").contains("expected exactly one"));
}

#[test]
fn recurrence_candidate_budget_has_distinct_mcp_outcomes() {
    let server = PolynomialToolsServer::new();
    let fibonacci = || FindRecurrenceRequest {
        polynomials: None,
        coefficients: Some(vec![
            vec![1.into()],
            vec![1.into()],
            vec![2.into()],
            vec![3.into()],
            vec![5.into()],
        ]),
        expressions: None,
        text: None,
        options: None,
        ..Default::default()
    };

    let mut zero_budget = fibonacci();
    zero_budget.options = Some(
        serde_json::from_value(json!({ "max_candidates": 0 })).expect("valid recurrence options"),
    );
    let Json(exhausted) = server.find_recurrence(Parameters(zero_budget)).unwrap();
    assert_eq!(exhausted.status, RecurrenceSearchStatus::BudgetExhausted);
    assert!(!exhausted.found);
    assert_eq!(exhausted.candidates_considered, Some(0));
    assert_eq!(exhausted.candidates_tried, Some(0));

    let mut exact_space = fibonacci();
    exact_space.options = Some(
        serde_json::from_value(json!({
            "max_candidates": 1,
            "min_rec_len": 1,
            "max_rec_len": 1,
            "max_var_deg": 0,
            "max_idx_deg": 0,
            "max_diff_deg": 0
        }))
        .expect("valid recurrence options"),
    );
    let Json(not_found) = server.find_recurrence(Parameters(exact_space)).unwrap();
    assert_eq!(not_found.status, RecurrenceSearchStatus::NotFound);
    assert_eq!(not_found.candidates_considered, Some(1));

    let Json(found) = server
        .find_recurrence(Parameters(FindRecurrenceRequest {
            polynomials: None,
            coefficients: Some(vec![
                vec![1.into()],
                vec![2.into()],
                vec![4.into()],
                vec![8.into()],
                vec![16.into()],
            ]),
            expressions: None,
            text: None,
            options: Some(
                serde_json::from_value(json!({ "max_candidates": 1 }))
                    .expect("valid recurrence options"),
            ),
            ..Default::default()
        }))
        .unwrap();
    assert_eq!(found.status, RecurrenceSearchStatus::Found);
    assert_eq!(found.candidates_considered, Some(1));
}

#[test]
fn family_check_reports_properties_and_recurrence() {
    let server = PolynomialToolsServer::new();
    let Json(response) = server
        .check_polynomial_family(Parameters(CheckPolynomialFamilyRequest {
            polynomials: Some(vec![
                coeffs(&[1]),
                coeffs(&[2]),
                coeffs(&[4]),
                coeffs(&[8]),
                coeffs(&[16]),
            ]),
            text: None,
            sequence: None,
            max_n: None,
            options: Some(FamilyCheckOptions {
                require_real_rooted: Some(true),
                require_simple_roots: None,
                require_palindromic: None,
                require_gamma_positive: None,
                require_unimodal: None,
                require_log_concave: None,
                require_ultra_log_concave: None,
                check_consecutive_interlacing: Some(false),
                require_consecutive_weak_interlacing: None,
                find_recurrence: Some(true),
            }),
            recurrence_options: None,
            lace: None,
        }))
        .unwrap();

    assert!(response.all_required_checks_passed);
    assert_eq!(response.item_count, 5);
    assert_eq!(
        response
            .recurrence
            .as_ref()
            .and_then(|r| r.recurrence.as_deref()),
        Some("P(n) = 2 P(n-1)")
    );
    assert!(response.markdown.contains("Polynomial properties"));
}

#[test]
fn family_check_can_require_finite_lace_tnn() {
    let server = PolynomialToolsServer::new();
    let Json(response) = server
        .check_polynomial_family(Parameters(CheckPolynomialFamilyRequest {
            polynomials: Some(vec![
                coeffs(&[2, 1]),
                coeffs(&[8, 6, 1]),
                coeffs(&[3, 4, 1]),
            ]),
            text: None,
            sequence: None,
            max_n: None,
            options: None,
            recurrence_options: None,
            lace: Some(LaceCheckRequest {
                block_rows: Some(1),
                block_cols: Some(3),
                max_minor_size: Some(3),
                include_matrix: Some(true),
            }),
        }))
        .unwrap();

    assert!(!response.all_required_checks_passed);
    assert_eq!(
        response.first_failure.as_deref(),
        Some("finite Lace truncation is not TNN up to minors of size 3")
    );
    assert_eq!(response.lace.as_ref().map(|l| l.tnn), Some(false));
    assert!(response
        .lace
        .as_ref()
        .and_then(|l| l.error.as_deref())
        .is_some_and(|error| error.contains("det = -1")));
}

#[test]
fn serializes_bigint_and_bigrational_as_strings() {
    let server = PolynomialToolsServer::new();
    let Json(resultant) = server
        .resultant(Parameters(BigIntPolynomialPairRequest {
            p: bigint_coeffs(&["1", "0", "1"]),
            q: bigint_coeffs(&["-1", "1"]),
        }))
        .unwrap();
    let value = serde_json::to_value(resultant).unwrap();
    assert_eq!(value["resultant"], "2");
    assert_eq!(value["p"]["coefficients"][0], "1");

    let Json(ehrhart) = server
        .ehrhart_hstar(Parameters(EhrhartHstarRequest {
            mode: EhrhartHstarMode::HstarToEhrhart,
            hstar: Some(vec![1.into(), 0.into(), 0.into()]),
            ehrhart_coefficients: None,
            numerator_coefficients: None,
            denominator: None,
        }))
        .unwrap();
    let value = serde_json::to_value(ehrhart).unwrap();
    assert_eq!(value["ehrhart_coefficients"][1], "3/2");
}

#[test]
fn mcp_rejects_oversized_polynomial_and_recurrence_inputs() {
    let oversized_expression = PolynomialInput {
        coefficients: None,
        expression: Some(format!("x^{MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL}")),
    };
    assert!(parse_polynomial_input(&oversized_expression)
        .unwrap_err()
        .contains("MCP limit"));

    let oversized_batch = PolynomialBatchInput {
        polynomials: Some(vec![coeffs(&[1]); MCP_MAX_BATCH_POLYNOMIALS + 1]),
        text: None,
    };
    assert!(format!("{:?}", parse_batch(&oversized_batch).unwrap_err()).contains("MCP limit"));

    let oversized_recurrence_batch = FindRecurrenceRequest {
        polynomials: None,
        coefficients: Some(vec![
            vec![RationalCoefficientInput::from(1)];
            MCP_MAX_RECURRENCE_INPUT_ROWS + 1
        ]),
        expressions: None,
        text: None,
        options: None,
        ..Default::default()
    };
    assert!(format!(
        "{:?}",
        parse_recurrence_batch_rational(&oversized_recurrence_batch).unwrap_err()
    )
    .contains("MCP limit"));

    let excessive_length: RecurrenceSearchOptionsInput =
        serde_json::from_value(json!({ "max_rec_len": MCP_MAX_RECURRENCE_LENGTH + 1 })).unwrap();
    assert!(format!(
        "{:?}",
        apply_recurrence_options(Some(excessive_length)).unwrap_err()
    )
    .contains("maximum recurrence length"));

    let excessive_grid: RecurrenceSearchOptionsInput = serde_json::from_value(json!({
        "max_rec_len": MCP_MAX_RECURRENCE_LENGTH,
        "max_var_deg": MCP_MAX_RECURRENCE_DEGREE,
        "max_idx_deg": MCP_MAX_RECURRENCE_DEGREE,
        "max_diff_deg": MCP_MAX_RECURRENCE_DEGREE
    }))
    .unwrap();
    assert!(format!(
        "{:?}",
        apply_recurrence_options(Some(excessive_grid)).unwrap_err()
    )
    .contains("candidates"));
}

#[test]
fn mcp_bounds_generation_and_reports_exact_conversion_errors() {
    let server = PolynomialToolsServer::new();
    assert!(server
        .generate_sequence(Parameters(GenerateSequenceRequest {
            sequence: SequenceKind::Eulerian,
            max_n: MCP_MAX_SEQUENCE_N + 1,
        }))
        .is_err());

    let nonintegral = server.ehrhart_hstar(Parameters(EhrhartHstarRequest {
        mode: EhrhartHstarMode::EhrhartToHstar,
        hstar: None,
        ehrhart_coefficients: Some(vec!["1/2".to_string()]),
        numerator_coefficients: None,
        denominator: None,
    }));
    assert!(
        format!("{:?}", nonintegral.err().expect("nonintegral error")).contains("not an integer")
    );

    let Json(found) = server
        .find_recurrence(Parameters(FindRecurrenceRequest {
            polynomials: Some(vec![coeffs(&[1]), coeffs(&[2]), coeffs(&[4]), coeffs(&[8])]),
            coefficients: None,
            expressions: None,
            text: None,
            options: None,
            ..Default::default()
        }))
        .unwrap();
    let recurrence_json = found.recurrence_json.expect("small recurrence found");
    let too_many_rows =
        server.generate_recurrence_rows(Parameters(GenerateRecurrenceRowsRequest {
            recurrence_json,
            rows: Some(MCP_MAX_RECURRENCE_ROWS + 1),
            additional: None,
        }));
    assert!(format!("{:?}", too_many_rows.err().expect("row limit error")).contains("MCP limit"));

    let oversized_lace = check_family_lace(
        &vec![normalize_polynomial(vec![1]); 4],
        &LaceCheckRequest {
            block_rows: Some(MCP_MAX_LACE_MATRIX_CELLS),
            block_cols: Some(2),
            max_minor_size: Some(1),
            include_matrix: Some(false),
        },
    );
    assert!(!oversized_lace.tnn);
    assert!(oversized_lace
        .error
        .is_some_and(|error| error.contains("exceeds MCP limits")));
}

fn tnn_request(text: &str, chain: bool) -> TnnNetworkRequest {
    TnnNetworkRequest {
        polynomials: None,
        text: Some(text.to_string()),
        matrix: None,
        chain: Some(chain),
    }
}

fn int_matrix(rows: &[&[i64]]) -> Vec<Vec<BigIntCoefficientInput>> {
    rows.iter()
        .map(|row| {
            row.iter()
                .map(|&v| BigIntCoefficientInput::from(v))
                .collect()
        })
        .collect()
}

fn brenti_request(text: &str) -> BrentiRecognitionRequest {
    BrentiRecognitionRequest {
        polynomials: None,
        text: Some(text.to_string()),
        first_row_index: None,
        max_target_span: None,
        weights: None,
        max_search_nodes: None,
        max_integer_search_nodes: None,
    }
}

#[test]
fn tnn_network_certifies_pascal_with_chain_polynomials() {
    let server = PolynomialToolsServer::new();
    let Json(response) = server
        .check_tnn_network(Parameters(tnn_request(
            "1,1\n1,2,1\n1,3,3,1\n1,4,6,4,1",
            true,
        )))
        .unwrap();
    assert_eq!(response.input_mode, "polynomials");
    assert!(response.prepended_p0);
    assert!(response.tnn);
    assert_eq!(response.failure, None);
    assert_eq!(response.certificate_verified, Some(true));
    let lambda = response.lambda.as_ref().unwrap();
    assert_eq!(lambda.len(), 4);
    for (n, row) in lambda.iter().enumerate() {
        assert_eq!(row, &vec!["1".to_string(); n + 1]);
    }
    let chain = response.chain.as_ref().unwrap();
    assert_eq!(chain.error, None);
    assert_eq!(chain.conclusion_holds, Some(true));
    assert_eq!(chain.contradicts_theorem, Some(false));
    assert_eq!(chain.polynomials[3].polynomial, "t + 6t^2 + 6t^3");
    assert!(chain.polynomials.iter().all(|p| p.roots_in_minus_one_zero));
    assert_eq!(chain.interlacing.len(), 4);
    assert!(chain.interlacing.iter().all(|pair| pair.interlaces));

    let value = serde_json::to_value(&response).unwrap();
    assert_eq!(value["coefficient_matrix"][2], json!(["1", "2", "1"]));
    assert_eq!(value["lambda"][0], json!(["1"]));
}

#[test]
fn tnn_network_guardrail_example_makes_no_row_claim() {
    let server = PolynomialToolsServer::new();
    let Json(response) = server
        .check_tnn_network(Parameters(TnnNetworkRequest {
            polynomials: None,
            text: None,
            matrix: Some(int_matrix(&[&[1, 0, 0], &[1, 1, 0], &[1, 1, 1]])),
            chain: None,
        }))
        .unwrap();
    assert_eq!(response.input_mode, "matrix");
    assert!(response.tnn);
    assert_eq!(response.certificate_verified, Some(true));
    assert_eq!(
        response.lambda,
        Some(vec![
            vec!["1".to_string()],
            vec!["1".to_string(), "0".to_string()]
        ])
    );
    assert_eq!(response.chain, None);
    assert!(response
        .note
        .contains("does not imply that the row polynomials are real-rooted"));
    let value = serde_json::to_value(&response).unwrap();
    assert!(value
        .as_object()
        .unwrap()
        .keys()
        .all(|key| !key.contains("real_rooted")));
}

#[test]
fn tnn_network_reports_negative_minor_and_shape_failures() {
    let server = PolynomialToolsServer::new();
    let Json(response) = server
        .check_tnn_network(Parameters(TnnNetworkRequest {
            polynomials: None,
            text: None,
            matrix: Some(int_matrix(&[&[1], &[1, 1], &[2, 0, 1]])),
            chain: Some(true),
        }))
        .unwrap();
    assert!(!response.tnn);
    assert_eq!(response.lambda, None);
    let failure = response.failure.as_ref().unwrap();
    assert_eq!(failure.kind, "not_tnn");
    assert_eq!(
        failure.witness,
        Some(NegativeMinorItem {
            rows: vec![1, 2],
            cols: vec![0, 1],
            value: "-2".to_string(),
        })
    );
    // The chain polynomials are still reported; the theorem does not apply.
    assert_eq!(
        response.chain.as_ref().unwrap().contradicts_theorem,
        Some(false)
    );

    let Json(non_monic) = server
        .check_tnn_network(Parameters(tnn_request("1,2\n1,2,1", false)))
        .unwrap();
    let failure = non_monic.failure.unwrap();
    assert_eq!(failure.kind, "non_monic_polynomial");
    assert_eq!(
        (failure.row, failure.value.as_deref()),
        (Some(1), Some("2"))
    );
    assert_eq!(non_monic.coefficient_matrix, None);

    let Json(above) = server
        .check_tnn_network(Parameters(TnnNetworkRequest {
            polynomials: None,
            text: None,
            matrix: Some(int_matrix(&[&[1, 2], &[1, 1]])),
            chain: Some(true),
        }))
        .unwrap();
    assert_eq!(above.failure.unwrap().kind, "nonzero_above_diagonal");
    assert!(above.chain.unwrap().error.is_some());
}

#[test]
fn tnn_network_validates_input_shape_and_limits() {
    let server = PolynomialToolsServer::new();
    let both = server.check_tnn_network(Parameters(TnnNetworkRequest {
        polynomials: None,
        text: Some("1".to_string()),
        matrix: Some(int_matrix(&[&[1]])),
        chain: None,
    }));
    assert!(format!("{:?}", both.err().unwrap()).contains("exactly one"));
    let bad = server.check_tnn_network(Parameters(tnn_request("1\n1,x", false)));
    assert!(format!("{:?}", bad.err().unwrap()).contains("polynomial 1"));
    let rows = crate::limits::MCP_MAX_TNN_NETWORK_ROWS + 1;
    let too_many = server.check_tnn_network(Parameters(TnnNetworkRequest {
        polynomials: None,
        text: None,
        matrix: Some(vec![vec![BigIntCoefficientInput::from(1)]; rows]),
        chain: None,
    }));
    assert!(format!("{:?}", too_many.err().unwrap()).contains("MCP limit"));
}

#[test]
fn brenti_recognition_certifies_stirling_and_rejects_eulerian() {
    let server = PolynomialToolsServer::new();
    let mut request = brenti_request("0,1\n0,1,1\n0,1,3,1\n0,1,7,6,1\n0,1,15,25,10,1");
    request.weights = Some(BrentiWeightsInput::Integer);
    let Json(stirling) = server
        .recognize_brenti_strip_digraph(Parameters(request))
        .unwrap();
    assert_eq!(stirling.verdict, "certified");
    assert!(stirling.certified);
    assert_eq!(stirling.certificate_verified, Some(true));
    assert_eq!(stirling.row_indices, Some(vec![1, 2, 3, 4, 5]));
    assert_eq!(stirling.options.weights, "integer");
    let digraph = stirling.digraph.unwrap();
    let edge = |to: usize, weight: usize| StripDigraphEdgeItem {
        to,
        weight: weight.to_string(),
    };
    assert_eq!(digraph.outgoing[0], vec![edge(1, 1)]);
    for y in 1..5 {
        assert_eq!(digraph.outgoing[y], vec![edge(y, y), edge(y + 1, 1)]);
    }

    let Json(eulerian) = server
        .recognize_brenti_strip_digraph(Parameters(brenti_request(
            "1,1\n1,4,1\n1,11,11,1\n1,26,66,26,1",
        )))
        .unwrap();
    assert_eq!(eulerian.verdict, "impossible_within_bounds");
    assert!(!eulerian.certified);
    assert_eq!(eulerian.digraph, None);
    assert_eq!(eulerian.supports_checked, Some(0));
    assert!(eulerian.obstruction.is_some());
}

#[test]
fn brenti_recognition_guardrail_budgets_and_limits() {
    let server = PolynomialToolsServer::new();
    let Json(guardrail) = server
        .recognize_brenti_strip_digraph(Parameters(brenti_request("1,1\n1,1,1")))
        .unwrap();
    assert_eq!(guardrail.verdict, "certified");
    assert!(guardrail
        .note
        .contains("does not show that the rows are real-rooted"));

    let mut tiny = brenti_request("1,1,1\n1,2,3,2,1");
    tiny.max_search_nodes = Some(1);
    let Json(inconclusive) = server
        .recognize_brenti_strip_digraph(Parameters(tiny))
        .unwrap();
    assert_eq!(inconclusive.verdict, "inconclusive");
    assert!(inconclusive.reason.is_some());

    let mut huge = brenti_request("1,1");
    huge.max_search_nodes = Some(crate::limits::MCP_MAX_BRENTI_SEARCH_NODES + 1);
    assert!(server
        .recognize_brenti_strip_digraph(Parameters(huge))
        .is_err());
    let mut bad_offset = brenti_request("1,1");
    bad_offset.first_row_index = Some(2);
    assert!(format!(
        "{:?}",
        server
            .recognize_brenti_strip_digraph(Parameters(bad_offset))
            .err()
            .unwrap()
    )
    .contains("first row index"));
    let degree = crate::limits::MCP_MAX_BRENTI_DEGREE + 1;
    let high = brenti_request(&format!("t^{degree}"));
    assert!(format!(
        "{:?}",
        server
            .recognize_brenti_strip_digraph(Parameters(high))
            .err()
            .unwrap()
    )
    .contains("degree"));
}

fn pf_network_request(text: &str) -> PfNetworkRequest {
    PfNetworkRequest {
        polynomials: None,
        text: Some(text.to_string()),
        first_row_index: None,
        max_step_length: None,
        max_rise: None,
        weights: None,
        max_search_nodes: None,
        max_integer_search_nodes: None,
    }
}

#[test]
fn pf_network_certifies_pascal_and_never_one_plus_t_plus_t_squared() {
    let server = PolynomialToolsServer::new();
    let Json(pascal) = server
        .recognize_pf_network(Parameters(pf_network_request("1,1\n1,2,1\n1,3,3,1")))
        .unwrap();
    assert_eq!(pascal.verdict, "certified");
    assert_eq!(pascal.certificate_verified, Some(true));
    assert_eq!(pascal.row_indices, Some(vec![1, 2, 3]));
    let network = pascal.network.unwrap();
    assert_eq!(network.last_column, 3);
    let step = |length: usize, rise: usize| PfNetworkStepItem {
        length,
        rise,
        weight: "1".to_string(),
    };
    assert_eq!(network.steps[0], vec![step(1, 0), step(1, 1)]);
    assert!(pascal.note.contains("Pólya frequency"));
    assert!(pascal
        .note
        .contains("does not certify that the coefficient matrix"));

    let Json(bad) = server
        .recognize_pf_network(Parameters(pf_network_request("1,1\n1,1,1")))
        .unwrap();
    assert_eq!(bad.verdict, "impossible_within_bounds");
    assert!(!bad.certified);
    assert!(bad.obstruction.unwrap().contains("not real-rooted"));

    let Json(eulerian) = server
        .recognize_pf_network(Parameters(pf_network_request(
            "1,1\n1,4,1\n1,11,11,1\n1,26,66,26,1",
        )))
        .unwrap();
    assert_eq!(eulerian.verdict, "impossible_within_bounds");

    let mut tiny = pf_network_request("1,1\n1,2,1\n1,3,3,1");
    tiny.max_search_nodes = Some(1);
    let Json(inconclusive) = server.recognize_pf_network(Parameters(tiny)).unwrap();
    assert_eq!(inconclusive.verdict, "inconclusive");

    let mut long = pf_network_request("1,1");
    long.max_step_length = Some(crate::limits::MCP_MAX_PF_NETWORK_STEP + 1);
    assert!(format!(
        "{:?}",
        server.recognize_pf_network(Parameters(long)).err().unwrap()
    )
    .contains("MCP limit"));
}

fn left_production_request(text: &str, first_row_index: Option<usize>) -> LeftProductionRequest {
    LeftProductionRequest {
        polynomials: None,
        text: Some(text.to_string()),
        first_row_index,
    }
}

#[test]
fn left_production_certifies_stirling_and_rejects_eulerian() {
    let server = PolynomialToolsServer::new();
    let Json(stirling) = server
        .check_left_production_matrix(Parameters(left_production_request(
            "0,1\n0,1,1\n0,1,3,1\n0,1,7,6,1",
            None,
        )))
        .unwrap();
    assert!(stirling.tnn);
    assert!(stirling.rows_certified_pf);
    assert_eq!(stirling.certificate_verified, Some(true));
    assert_eq!(stirling.production_matrix[4], vec!["0", "1", "3", "3", "1"]);
    assert!(stirling.production_network.is_some());

    let Json(eulerian) = server
        .check_left_production_matrix(Parameters(left_production_request(
            "1,1\n1,4,1\n1,11,11,1\n1,26,66,26,1",
            None,
        )))
        .unwrap();
    assert!(!eulerian.tnn);
    assert!(!eulerian.rows_certified_pf);
    assert_eq!(eulerian.production_matrix[4][2], "-5");
    assert_eq!(eulerian.negative_minor.unwrap().value, "-5");
    assert!(eulerian.note.contains("nothing follows about the rows"));

    let Json(ones) = server
        .check_left_production_matrix(Parameters(left_production_request(
            "1\n1,1\n1,1,1",
            Some(0),
        )))
        .unwrap();
    assert!(!ones.rows_certified_pf);

    let bad =
        server.check_left_production_matrix(Parameters(left_production_request("1,1,1", None)));
    assert!(format!("{:?}", bad.err().unwrap()).contains("degree"));
}

#[test]
fn recurrence_transform_tool_attaches_root_domain_and_enclosures() {
    let server = PolynomialToolsServer::new();
    let request = json!({"oeis_id":"A008292","row_count":6,"root_domain":{},
        "transformation":{"kind":"f_polynomial",
            "schedule":{"first_index":1,"first_degree":0,"step":1}}});
    let Json(report) = server
        .transform_recurrence(Parameters(RecurrenceTransformToolRequest { request }))
        .unwrap();
    let root = &report["root_domain"];
    assert_eq!(root["recurrence_theorem"]["matched"], true);
    assert!(root["claims"]
        .as_array()
        .unwrap()
        .iter()
        .any(|c| c["evidence"] == "completed_recurrence_induction"));
    let Json(checked) = server
        .verify_recurrence_transform(Parameters(VerifyRecurrenceTransformToolRequest {
            report: json!(report),
        }))
        .unwrap();
    assert_eq!(checked["verified"], true);

    let Json(enclosures) = server
        .exact_root_enclosures(Parameters(RootEnclosuresToolRequest {
            request: json!({"polynomials": [["0", "-2", "0", "1"]], "max_width": "1/100",
                "interval": {"lower": {"value": "-2"}, "upper": {"value": "2"}}}),
        }))
        .unwrap();
    let item = &enclosures["items"][0];
    assert_eq!(item["real_root_count"], 3);
    assert_eq!(item["roots"][1]["lower"], "0");
    assert_eq!(item["real_rooted_in_interval"], true);
    assert!(server
        .exact_root_enclosures(Parameters(RootEnclosuresToolRequest {
            request: json!({"polynomials": [["1"]], "unknown": 1}),
        }))
        .is_err());
}

#[test]
fn recurrence_transform_tool_transports_and_replays() {
    let server = PolynomialToolsServer::new();
    let request = json!({"oeis_id":"A008292","row_count":6,"search":{},
        "transformation":{"kind":"f_polynomial",
            "schedule":{"first_index":1,"first_degree":0,"step":1}}});
    let Json(report) = server
        .transform_recurrence(Parameters(RecurrenceTransformToolRequest {
            request: request.clone(),
        }))
        .unwrap();
    assert_eq!(report["schema"], "polytool.recurrence-transform.v1");
    assert_eq!(report["transported"]["evidence"], "transported_identity");
    assert_eq!(report["independent_search"]["evidence"], "finite_fit");
    let Json(checked) = server
        .verify_recurrence_transform(Parameters(VerifyRecurrenceTransformToolRequest {
            report: json!(report),
        }))
        .unwrap();
    assert_eq!(checked["verified"], true);
    let mut tampered = json!(report);
    tampered["row_origin"] = json!(0);
    assert!(server
        .verify_recurrence_transform(Parameters(VerifyRecurrenceTransformToolRequest {
            report: tampered,
        }))
        .is_err());
    let mut excessive = request;
    excessive["limits"] = json!({"max_rows": 1_000_000});
    assert!(server
        .transform_recurrence(Parameters(RecurrenceTransformToolRequest {
            request: excessive
        }))
        .is_err());
}

#[test]
fn recurrence_implicit_tools_transport_assess_and_replay() {
    let server = PolynomialToolsServer::new();
    let request = json!({"oeis_id":"A007318","row_count":6,
        "transformation":{"kind":"diagonal",
            "weight":{"kind":"factorial_power","power":1}}});
    let Json(report) = server
        .transform_recurrence_implicit(Parameters(RecurrenceImplicitToolRequest {
            request: request.clone(),
        }))
        .unwrap();
    assert_eq!(report["schema"], "polytool.recurrence-implicit.v1");
    assert_eq!(report["assessment"]["status"], "equivalent_description");
    assert_eq!(
        report["transports"][0]["relation"]["text"],
        "D Q[n+1] - D Q[n] - Q[n] = 0"
    );
    let Json(checked) = server
        .verify_recurrence_implicit(Parameters(VerifyRecurrenceImplicitToolRequest {
            report: json!(report),
        }))
        .unwrap();
    assert_eq!(checked["verified"], true);
    let mut description = report["description"].clone();
    description["exceptional_index_conditions"] = json!([]);
    let Json(assessment) = server
        .assess_generating_description(Parameters(AssessGeneratingDescriptionToolRequest {
            description,
        }))
        .unwrap();
    assert_eq!(assessment["status"], "certified_consequence");
    let mut tampered = json!(report);
    tampered["row_origin"] = json!(1);
    assert!(server
        .verify_recurrence_implicit(Parameters(VerifyRecurrenceImplicitToolRequest {
            report: tampered,
        }))
        .is_err());
    let mut excessive = request;
    excessive["limits"] = json!({"max_rows": 1_000_000});
    assert!(server
        .transform_recurrence_implicit(Parameters(RecurrenceImplicitToolRequest {
            request: excessive
        }))
        .is_err());
}

#[test]
fn representation_search_tools_round_trip_and_replay() {
    let server = PolynomialToolsServer::new();
    let request = json!({"oeis_id": "A007318", "row_count": 8,
        "limits": {"max_depth": 1},
        "transformations": [{"kind": "diagonal",
            "weight": {"kind": "factorial_power", "power": 1}}]});
    let Json(report) = server
        .search_recurrence_representations(Parameters(RecurrenceRepresentationSearchToolRequest {
            request: request.clone(),
        }))
        .unwrap();
    assert_eq!(
        report["schema"],
        "polytool.recurrence-representation-search.v1"
    );
    assert_eq!(report["status"], "complete");
    let descriptions = report["representations"][1]["descriptions"]
        .as_array()
        .unwrap();
    assert!(descriptions.iter().any(
        |d| d["label"] == "equivalent_generating_description" && d["kind"] == "certified_space"
    ));
    // The MCP JSON is the crate report.
    let parsed: polytool::recurrence::RepresentationSearchReport =
        serde_json::from_value(json!(report)).unwrap();
    assert_eq!(json!(parsed), json!(report));
    let Json(checked) = server
        .verify_recurrence_representation_search(Parameters(
            VerifyRecurrenceRepresentationSearchToolRequest {
                report: json!(report),
            },
        ))
        .unwrap();
    assert_eq!(checked["verified"], true);
    let mut tampered = json!(report);
    tampered["representations"][1]["states"][2][0][1] = json!("3");
    assert!(server
        .verify_recurrence_representation_search(Parameters(
            VerifyRecurrenceRepresentationSearchToolRequest { report: tampered }
        ))
        .is_err());
    let mut excessive = request;
    excessive["limits"] = json!({"max_depth": 100});
    assert!(server
        .search_recurrence_representations(Parameters(RecurrenceRepresentationSearchToolRequest {
            request: excessive
        }))
        .is_err());
}

fn enumerate_request(value: serde_json::Value) -> EnumerateRealRootedRequest {
    serde_json::from_value(value).expect("valid enumerate_real_rooted request")
}

#[test]
fn enumerate_real_rooted_lists_the_cubic_family_in_ascending_order() {
    let server = PolynomialToolsServer::new();
    let Json(response) = server
        .enumerate_real_rooted(Parameters(enumerate_request(json!({
            "degree": 3, "leading_coefficients_descending": [1, 0, -3]
        }))))
        .unwrap();
    assert_eq!(response.command, "enumerate-real-rooted");
    assert_eq!(response.status, "complete");
    assert_eq!(response.limit_reached, None);
    assert_eq!(response.output_coefficient_order, "ascending");
    assert_eq!(response.leading_coefficients_descending, ["1", "0", "-3"]);
    // x^3 - 3x + c is real-rooted exactly for c = -2, ..., 2.
    let expected: Vec<Vec<String>> = (-2..=2)
        .map(|c: i32| vec![c.to_string(), "-3".into(), "0".into(), "1".into()])
        .collect();
    assert_eq!(response.polynomials, expected);
    assert_eq!(response.count, 5);
    assert_eq!(response.stats.results_emitted, 5);
    assert_eq!(response.limits.max_nodes, 5_000);
    assert_eq!(response.limits.max_results, 500);
    assert_eq!(response.limits.max_coefficient_bits, 1_024);
    // The JSON field names are those of the CLI.
    let value = json!(response);
    for field in [
        "command",
        "degree",
        "leading_coefficients_descending",
        "output_coefficient_order",
        "status",
        "limit_reached",
        "limits",
        "count",
        "polynomials",
        "stats",
    ] {
        assert!(value.get(field).is_some(), "missing field {field}");
    }
    assert!(value["limit_reached"].is_null());
    for field in [
        "derivative_levels",
        "nodes_examined",
        "polynomials_by_level",
    ] {
        assert!(value["stats"].get(field).is_some(), "missing stat {field}");
    }
}

#[test]
fn enumerate_real_rooted_accepts_string_coefficients_and_the_cli_alias() {
    let server = PolynomialToolsServer::new();
    let Json(response) = server
        .enumerate_real_rooted(Parameters(enumerate_request(json!({
            "degree": 4, "leading_desc": ["1", "0", "-4"]
        }))))
        .unwrap();
    assert_eq!(response.status, "complete");
    // (x^2 - 2)^2 and x^2 (x^2 - 4) are real-rooted members.
    for expected in [["4", "0", "-4", "0", "1"], ["0", "0", "-4", "0", "1"]] {
        assert!(response
            .polynomials
            .iter()
            .any(|p| p.iter().map(String::as_str).eq(expected)));
    }
    // A negative leading coefficient is allowed.
    let Json(negative) = server
        .enumerate_real_rooted(Parameters(enumerate_request(json!({
            "degree": 3, "leading_coefficients_descending": ["-1", "0", "3"]
        }))))
        .unwrap();
    assert_eq!(negative.count, 5);
}

#[test]
fn enumerate_real_rooted_reports_each_budget_as_a_resource_limit() {
    let server = PolynomialToolsServer::new();
    let Json(by_results) = server
        .enumerate_real_rooted(Parameters(enumerate_request(json!({
            "degree": 3, "leading_coefficients_descending": [1, 0, -3], "max_results": 2
        }))))
        .unwrap();
    assert_eq!(by_results.status, "resource_limit");
    assert_eq!(by_results.limit_reached.as_deref(), Some("max_results"));
    assert_eq!(by_results.count, 2);
    // An initial segment of the complete sorted list.
    assert_eq!(by_results.polynomials[0][0], "-2");
    assert_eq!(by_results.polynomials[1][0], "-1");

    let Json(by_nodes) = server
        .enumerate_real_rooted(Parameters(enumerate_request(json!({
            "degree": 10, "leading_coefficients_descending": [1, 0, -15], "max_nodes": 3
        }))))
        .unwrap();
    assert_eq!(by_nodes.status, "resource_limit");
    assert_eq!(by_nodes.limit_reached.as_deref(), Some("max_nodes"));
    assert!(by_nodes.stats.nodes_examined <= 3);

    let Json(by_bits) = server
        .enumerate_real_rooted(Parameters(enumerate_request(json!({
            "degree": 8, "leading_coefficients_descending": [1, 0, -15],
            "max_coefficient_bits": 4
        }))))
        .unwrap();
    assert_eq!(by_bits.status, "resource_limit");
    assert_eq!(
        by_bits.limit_reached.as_deref(),
        Some("max_coefficient_bits")
    );
}

#[test]
fn enumerate_real_rooted_rejects_invalid_and_excessive_requests() {
    let server = PolynomialToolsServer::new();
    let call = |value: serde_json::Value| {
        server.enumerate_real_rooted(Parameters(enumerate_request(value)))
    };
    let message = |value: serde_json::Value| {
        call(value)
            .err()
            .expect("request is rejected")
            .message
            .to_string()
    };
    assert!(
        message(json!({"degree": 1, "leading_coefficients_descending": [1, 0, 1]}))
            .contains("degree must be at least 2")
    );
    assert!(
        message(json!({"degree": 3, "leading_coefficients_descending": [1, 0]}))
            .contains("leading coefficients")
    );
    assert!(
        message(json!({"degree": 3, "leading_coefficients_descending": [0, 0, 1]}))
            .contains("nonzero")
    );
    assert!(
        message(json!({"degree": 3, "leading_coefficients_descending": [1, 0, "x"]}))
            .contains("invalid integer")
    );
    // Ceilings on the degree, the budgets and the prefix size.
    assert!(
        message(json!({"degree": 33, "leading_coefficients_descending": [1, 0, -3]}))
            .contains("degree = 33 exceeds the MCP limit 32")
    );
    for (name, value) in [
        ("max_nodes", 25_001u64),
        ("max_results", 2_001),
        ("max_coefficient_bits", 2_049),
    ] {
        let text = message(json!({"degree": 3,
            "leading_coefficients_descending": [1, 0, -3], name: value}));
        assert!(
            text.contains(name) && text.contains("exceeds the MCP limit"),
            "{text}"
        );
    }
    let huge = format!("1{}", "0".repeat(100));
    assert!(
        message(json!({"degree": 3, "leading_coefficients_descending": [1, 0, huge]}))
            .contains("bits")
    );
    assert!(
        message(json!({"degree": 32, "leading_coefficients_descending": vec![1; 34]}))
            .contains("at most 33 leading coefficients")
    );
    // Unknown fields are refused by the schema.
    assert!(serde_json::from_value::<EnumerateRealRootedRequest>(
        json!({"degree": 3, "leading_coefficients_descending": [1, 0, -3], "unbounded": true})
    )
    .is_err());
}

#[test]
fn enumerate_real_rooted_schema_documents_the_cli_fields() {
    let schema = serde_json::to_value(schemars::schema_for!(EnumerateRealRootedRequest)).unwrap();
    assert_eq!(schema["additionalProperties"], false);
    let properties = schema["properties"].as_object().expect("schema properties");
    for name in [
        "degree",
        "leading_coefficients_descending",
        "max_results",
        "max_nodes",
        "max_coefficient_bits",
    ] {
        assert!(properties.contains_key(name), "missing {name}");
    }
    let required = schema["required"].as_array().expect("required fields");
    assert!(required.contains(&json!("degree")));
    assert!(required.contains(&json!("leading_coefficients_descending")));
    assert!(!required.contains(&json!("max_nodes")));
}
