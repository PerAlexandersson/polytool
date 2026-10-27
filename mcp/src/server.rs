//! The MCP server: the tool router with one method per tool, and the
//! `ServerHandler` metadata.

use crate::family::{
    check_family_lace, family_markdown, family_recurrence_response, family_report,
    family_source_and_polynomials, first_family_failure, generated_sequence_polynomials_bigint,
};
use crate::limits::{
    MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL, MCP_MAX_COEFFICIENT_TEXT_BYTES, MCP_MAX_CYCLIC_ORDER,
    MCP_MAX_RECURRENCE_ROWS, MCP_MAX_SEQUENCE_N, MCP_MAX_TOTAL_COEFFICIENTS,
};
use crate::parsing::{
    bigint_strings, collect_rational_polynomials_or_errors, invalid_params,
    normalize_bigint_polynomial, parse_batch, parse_bigint_batch, parse_bigint_coefficient,
    parse_bigint_coefficients, parse_bigint_matrix, parse_rational,
    parse_recurrence_batch_rational, parse_required_bigint_polynomial, validate_mcp_text,
};
use crate::planar_network::{
    brenti_recognition_response, left_production_response, pf_network_response,
    tnn_network_response,
};
use crate::recurrence_options::{
    apply_recurrence_options, recurrence_failure_response, recurrence_json_string,
    validate_recurrence_for_generation,
};
use crate::report_items::{
    bigint_matrix_strings, coefficient_tests_item, cyclic_sequence_item_response,
    cyclic_sieving_report_item, displayed, format_rational, hstar_inequalities_item,
    interlacing_pair_has_interlacing, interlacing_pair_item, interlacing_result,
    interlacing_result_i64, oeis_sequence_summary, orbit_decomposition_item,
    orbit_sizes_check_item, parse_bigint_items, parse_items, polynomial_properties_item_bigint,
};
use crate::requests::{
    AssessGeneratingDescriptionToolRequest, RecurrenceImplicitToolRequest,
    VerifyRecurrenceImplicitToolRequest,
};
use crate::requests::{
    BicyclicSievingRequest, BigIntPolynomialBatchInput, BigIntPolynomialPairRequest,
    BrentiRecognitionRequest, CheckPolynomialFamilyRequest, CyclicSievingRequest,
    CyclicSievingSequenceRequest, EhrhartHstarMode, EhrhartHstarRequest, FindRecurrenceRequest,
    GenerateOeisRowsRequest, GenerateRecurrenceRowsRequest, GenerateSequenceRequest,
    GetOeisSequenceRequest, HStarInequalitiesRequest, ListOeisSequencesRequest,
    PolynomialBatchInput, RecurrenceSpacesToolRequest, TnnNetworkRequest,
    VerifyRecurrenceSpacesToolRequest,
};
use crate::requests::{ImplicitRelationSearchToolRequest, VerifyImplicitRelationSearchToolRequest};
use crate::requests::{
    RecurrenceRepresentationSearchToolRequest, VerifyRecurrenceRepresentationSearchToolRequest,
};
use crate::requests::{
    RecurrenceTransformToolRequest, RootEnclosuresToolRequest, VerifyRecurrenceTransformToolRequest,
};
use crate::responses::{
    BicyclicSievingResponse, BrentiRecognitionResponse, CheckPolynomialFamilyResponse,
    CoefficientTestsResponse, CyclicSievingResponse, CyclicSievingSequenceResponse,
    DecompositionItem, DecompositionReport, DecompositionResponse, DiscriminantItem,
    DiscriminantResponse, EhrhartHstarResponse, FindRecurrenceResponse, GenerateOeisRowsResponse,
    GenerateRecurrenceRowsResponse, GenerateSequenceResponse, GetOeisSequenceResponse,
    HStarInequalitiesResponse, InterlacingPairItem, InterlacingPairResult, InterlacingProfileItem,
    InterlacingProfileResponse, InterlacingSequenceResponse, ListOeisSequencesResponse,
    MagicBasisReport, OeisPolynomialRow, ParsePolynomialsResponse, PolynomialPropertiesItem,
    PolynomialPropertiesResponse, RealRootsItem, RealRootsResponse, RecurrenceSearchStatus,
    ResultantResponse, TnnNetworkResponse,
};
use num_bigint::BigInt;
use polytool::recurrence::{
    find_recurrence_adaptive_rational_with_budget, find_recurrence_representations_rational,
    format_rational_coeff, AdaptiveSearchOutcome, BigRational, RecurrenceJson,
};
use polytool::*;
use rmcp::{
    handler::server::{router::tool::ToolRouter, wrapper::Parameters},
    model::{Implementation, ServerCapabilities, ServerInfo},
    tool, tool_handler, tool_router, ErrorData as McpError, Json, ServerHandler,
};
use serde_json::json;
use std::collections::BTreeMap;

#[derive(Debug, Clone)]
pub struct PolynomialToolsServer {
    tool_router: ToolRouter<Self>,
}

#[tool_router(router = tool_router)]
impl PolynomialToolsServer {
    pub fn new() -> Self {
        Self {
            tool_router: Self::tool_router(),
        }
    }

    #[tool(description = "Parse and format one or more dense univariate polynomials.")]
    pub fn parse_polynomials(
        &self,
        Parameters(input): Parameters<PolynomialBatchInput>,
    ) -> Result<Json<ParsePolynomialsResponse>, McpError> {
        let batch = parse_batch(&input)?;
        Ok(Json(ParsePolynomialsResponse {
            items: parse_items(&batch),
        }))
    }

    #[tool(
        description = "Check a polynomial family in one exact batch: properties, consecutive interlacing, optional finite Lace TNN, and optional recurrence search."
    )]
    pub fn check_polynomial_family(
        &self,
        Parameters(input): Parameters<CheckPolynomialFamilyRequest>,
    ) -> Result<Json<CheckPolynomialFamilyResponse>, McpError> {
        let (source, polynomials) = family_source_and_polynomials(&input)?;
        let polynomials = match polynomials {
            Ok(polynomials) => polynomials,
            Err(parse_errors) => {
                let first_failure = Some("one or more polynomials failed to parse".to_string());
                let mut markdown = String::new();
                markdown.push_str("# Polynomial family check\n\n");
                markdown.push_str(&format!("- Source: `{source}`\n"));
                markdown.push_str("- Required checks: failed\n");
                markdown.push_str("- First failure: one or more polynomials failed to parse\n\n");
                markdown.push_str("## Parse errors\n\n");
                for item in &parse_errors {
                    markdown.push_str(&format!(
                        "- polynomial {}: {}\n",
                        item.index,
                        item.error.as_deref().unwrap_or("parse error")
                    ));
                }
                return Ok(Json(CheckPolynomialFamilyResponse {
                    source,
                    item_count: 0,
                    all_required_checks_passed: false,
                    first_failure,
                    parse_errors,
                    items: Vec::new(),
                    consecutive_pairs: Vec::new(),
                    lace: None,
                    recurrence: None,
                    markdown,
                }));
            }
        };

        let options = input.options.as_ref();
        let find_recurrence = options.and_then(|o| o.find_recurrence).unwrap_or(false);
        let recurrence_search = if find_recurrence {
            Some(apply_recurrence_options(input.recurrence_options.clone())?)
        } else {
            None
        };
        let check_consecutive_interlacing = options
            .and_then(|o| o.check_consecutive_interlacing)
            .unwrap_or(true);

        let items: Vec<_> = polynomials
            .iter()
            .cloned()
            .enumerate()
            .map(|(index, polynomial)| family_report(index, polynomial))
            .collect();

        let consecutive_pairs = if check_consecutive_interlacing {
            polynomials
                .windows(2)
                .enumerate()
                .map(|(pair_index, pair)| InterlacingPairItem {
                    pair_index,
                    left_index: pair_index,
                    right_index: pair_index + 1,
                    ok: true,
                    result: Some(interlacing_result_i64(pair[0].clone(), pair[1].clone())),
                    error: None,
                })
                .collect()
        } else {
            Vec::new()
        };

        let lace = input
            .lace
            .as_ref()
            .map(|request| check_family_lace(&polynomials, request));

        let coefficients: Vec<Vec<i64>> =
            polynomials.iter().map(|p| p.coefficients.clone()).collect();
        let recurrence = recurrence_search
            .as_ref()
            .map(|(search, budget)| family_recurrence_response(&coefficients, search, *budget));

        let first_failure = first_family_failure(
            &items,
            &consecutive_pairs,
            lace.as_ref(),
            recurrence.as_ref(),
            options,
        );
        let all_required_checks_passed = first_failure.is_none();
        let markdown = family_markdown(
            &source,
            all_required_checks_passed,
            first_failure.as_deref(),
            &items,
            &consecutive_pairs,
            lace.as_ref(),
            recurrence.as_ref(),
        );

        Ok(Json(CheckPolynomialFamilyResponse {
            source,
            item_count: items.len(),
            all_required_checks_passed,
            first_failure,
            parse_errors: Vec::new(),
            items,
            consecutive_pairs,
            lace,
            recurrence,
            markdown,
        }))
    }

    #[tool(
        description = "Check real-rootedness, gamma-positivity, unimodality, and concavity properties."
    )]
    pub fn polynomial_properties(
        &self,
        Parameters(input): Parameters<BigIntPolynomialBatchInput>,
    ) -> Result<Json<PolynomialPropertiesResponse>, McpError> {
        let batch = parse_bigint_batch(&input)?;
        let items = batch
            .into_iter()
            .enumerate()
            .map(|(index, item)| match item {
                Ok(polynomial) => polynomial_properties_item_bigint(index, polynomial),
                Err(error) => PolynomialPropertiesItem {
                    index,
                    ok: false,
                    error: Some(error),
                    polynomial: None,
                    real_rooted: None,
                    simple_roots: None,
                    palindromic: None,
                    gamma_positive: None,
                    gamma_coefficients: None,
                    unimodal: None,
                    log_concave: None,
                    ultra_log_concave: None,
                },
            })
            .collect();
        Ok(Json(PolynomialPropertiesResponse { items }))
    }

    #[tool(
        description = "Check Newton inequalities and Kurtz's sufficient real-rootedness criterion for each polynomial."
    )]
    pub fn coefficient_tests(
        &self,
        Parameters(input): Parameters<BigIntPolynomialBatchInput>,
    ) -> Result<Json<CoefficientTestsResponse>, McpError> {
        let batch = parse_bigint_batch(&input)?;
        let items = batch
            .into_iter()
            .enumerate()
            .map(|(index, item)| coefficient_tests_item(index, item))
            .collect();
        Ok(Json(CoefficientTestsResponse { items }))
    }

    #[tool(
        description = "Check named Ehrhart h*-vector inequalities, including Stanley, Hibi, Balletti-Higashitani, and Stapledon-derived checks."
    )]
    pub fn hstar_inequalities(
        &self,
        Parameters(input): Parameters<HStarInequalitiesRequest>,
    ) -> Result<Json<HStarInequalitiesResponse>, McpError> {
        let batch = parse_bigint_batch(&BigIntPolynomialBatchInput {
            polynomials: input.polynomials,
            text: input.text,
        })?;
        let items = batch
            .into_iter()
            .enumerate()
            .map(|(index, item)| hstar_inequalities_item(index, item, input.dimension))
            .collect();
        Ok(Json(HStarInequalitiesResponse { items }))
    }

    #[tool(
        description = "Check strict and weak interlacing for a pair of polynomials with arbitrary-size integer coefficients."
    )]
    pub fn check_interlacing_pair(
        &self,
        Parameters(input): Parameters<BigIntPolynomialPairRequest>,
    ) -> Result<Json<InterlacingPairResult>, McpError> {
        let p = normalize_bigint_polynomial(parse_required_bigint_polynomial(&input.p, "p")?);
        let q = normalize_bigint_polynomial(parse_required_bigint_polynomial(&input.q, "q")?);
        Ok(Json(interlacing_result(p, q)))
    }

    #[tool(
        description = "Check strict and weak interlacing for consecutive polynomial pairs with arbitrary-size integer coefficients."
    )]
    pub fn check_interlacing_sequence(
        &self,
        Parameters(input): Parameters<BigIntPolynomialBatchInput>,
    ) -> Result<Json<InterlacingSequenceResponse>, McpError> {
        let batch = parse_bigint_batch(&input)?;
        let items = parse_bigint_items(&batch);
        let pairs = batch
            .windows(2)
            .enumerate()
            .map(|(pair_index, pair)| {
                interlacing_pair_item(pair_index, pair_index, pair_index + 1, &pair[0], &pair[1])
            })
            .collect();
        Ok(Json(InterlacingSequenceResponse { items, pairs }))
    }

    #[tool(
        description = "For each polynomial, count backward consecutive previous interlacings until the first failure. Accepts arbitrary-size integer coefficients."
    )]
    pub fn check_interlacing_profile(
        &self,
        Parameters(input): Parameters<BigIntPolynomialBatchInput>,
    ) -> Result<Json<InterlacingProfileResponse>, McpError> {
        let batch = parse_bigint_batch(&input)?;
        let items = parse_bigint_items(&batch);
        let mut pair_index = 0;
        let profile = batch
            .iter()
            .enumerate()
            .map(|(index, item)| match item {
                Err(error) => InterlacingProfileItem {
                    index,
                    ok: false,
                    polynomial: None,
                    previous_count: index,
                    checked_previous_count: 0,
                    interlacing_previous_count: 0,
                    strict_previous_count: 0,
                    weak_previous_count: 0,
                    previous_interlacing_indices: Vec::new(),
                    previous: Vec::new(),
                    error: Some(error.clone()),
                },
                Ok(polynomial) => {
                    let mut previous = Vec::new();
                    for left_index in (0..index).rev() {
                        let pair = interlacing_pair_item(
                            pair_index,
                            left_index,
                            index,
                            &batch[left_index],
                            item,
                        );
                        pair_index += 1;
                        let interlaces = interlacing_pair_has_interlacing(&pair);
                        previous.push(pair);
                        if !interlaces {
                            break;
                        }
                    }
                    let previous_interlacing_indices = previous
                        .iter()
                        .filter(|pair| interlacing_pair_has_interlacing(pair))
                        .map(|pair| pair.left_index)
                        .collect::<Vec<_>>();
                    InterlacingProfileItem {
                        index,
                        ok: true,
                        polynomial: Some(polynomial.display.clone()),
                        previous_count: index,
                        checked_previous_count: previous.len(),
                        interlacing_previous_count: previous_interlacing_indices.len(),
                        strict_previous_count: previous
                            .iter()
                            .filter(|pair| {
                                pair.result.as_ref().and_then(|result| result.strict) == Some(true)
                            })
                            .count(),
                        weak_previous_count: previous
                            .iter()
                            .filter(|pair| {
                                pair.result.as_ref().and_then(|result| result.weak) == Some(true)
                            })
                            .count(),
                        previous_interlacing_indices,
                        previous,
                        error: None,
                    }
                }
            })
            .collect();
        Ok(Json(InterlacingProfileResponse { items, profile }))
    }

    #[tool(
        description = "Return rational midpoint representatives for isolated real roots, or real_rooted=false. Midpoints are neither roots nor certified bounds; use exact_root_enclosures for certified enclosures with multiplicities."
    )]
    pub fn real_roots(
        &self,
        Parameters(input): Parameters<PolynomialBatchInput>,
    ) -> Result<Json<RealRootsResponse>, McpError> {
        let batch = parse_batch(&input)?;
        let items = batch
            .into_iter()
            .enumerate()
            .map(|(index, item)| match item {
                Ok(polynomial) => match polytool::real_roots(&polynomial.coefficients) {
                    Some(roots) => RealRootsItem {
                        index,
                        ok: true,
                        error: None,
                        polynomial: Some(polynomial),
                        real_rooted: Some(true),
                        roots: Some(roots.into_iter().map(format_rational).collect()),
                    },
                    None => RealRootsItem {
                        index,
                        ok: true,
                        error: None,
                        polynomial: Some(polynomial),
                        real_rooted: Some(false),
                        roots: None,
                    },
                },
                Err(error) => RealRootsItem {
                    index,
                    ok: false,
                    error: Some(error),
                    polynomial: None,
                    real_rooted: None,
                    roots: None,
                },
            })
            .collect();
        Ok(Json(RealRootsResponse { items }))
    }

    #[tool(
        description = "Search adaptively for a polynomial recurrence. Provide exactly one of `polynomials`, `coefficients`, `expressions`, or `text`. Pass search controls such as `min_rec_len`, `max_rec_len`, degree bounds, `try_denominator`, `fit_extra_rows`, `modular_prefilter`, and `max_candidates` as top-level arguments. Legacy nested `options` is accepted; a top-level control wins when both are supplied. Set `include_code` to false for a compact result without Mathematica, Python, Sage, or recurrence JSON. When the found recurrence matches a supported interlacing theorem (weighted renewal, fixed-delay predecessor), the additive field theorem_matches lists conditional matches: exact hypotheses, conclusion, window and a finite boundary-row check. They are separate from the finite-fit evidence and never prove the recurrence for all n."
    )]
    pub fn find_recurrence(
        &self,
        Parameters(input): Parameters<FindRecurrenceRequest>,
    ) -> Result<Json<FindRecurrenceResponse>, McpError> {
        let recurrence_options = input.effective_options();
        let include_code = input.include_code.unwrap_or(true);
        let batch = parse_recurrence_batch_rational(&input)?;
        let polynomials = match collect_rational_polynomials_or_errors(batch) {
            Ok(polynomials) => polynomials,
            Err(parse_errors) => {
                return Ok(Json(recurrence_failure_response(
                    RecurrenceSearchStatus::InvalidInput,
                    "one or more polynomials failed to parse",
                    None,
                    parse_errors,
                )));
            }
        };
        if polynomials.len() < 3 {
            return Ok(Json(recurrence_failure_response(
                RecurrenceSearchStatus::InvalidInput,
                "need at least 3 polynomials",
                None,
                Vec::new(),
            )));
        }
        let (search, budget) = apply_recurrence_options(recurrence_options)?;
        match find_recurrence_adaptive_rational_with_budget(&polynomials, &search, budget) {
            AdaptiveSearchOutcome::Found(result) => Ok(Json(FindRecurrenceResponse {
                status: RecurrenceSearchStatus::Found,
                found: true,
                recurrence: Some(format!("{}", result.recurrence)),
                latex: Some(result.recurrence.to_latex()),
                mathematica: include_code.then(|| {
                    result
                        .recurrence
                        .to_mathematica_definition_rational(&polynomials)
                }),
                sage: include_code
                    .then(|| result.recurrence.to_sage_definition_rational(&polynomials)),
                python: include_code.then(|| {
                    result
                        .recurrence
                        .to_python_definition_rational(&polynomials)
                }),
                recurrence_json: include_code.then(|| {
                    recurrence_json_string(&result, &polynomials, &search, polynomials.len())
                }),
                unknowns: Some(result.num_unknowns),
                weighted_unknowns: Some(result.weighted_unknowns),
                equations: Some(result.num_equations),
                fit_polynomials: Some(result.fit_polynomials),
                verification_polynomials: Some(result.verification_polynomials),
                candidates_tried: Some(result.candidates_tried),
                candidates_considered: Some(result.diagnostics.considered_candidates),
                error: None,
                parse_errors: Vec::new(),
                theorem_matches: crate::theorems::search_theorem_matches(
                    &result,
                    &polynomials,
                    &search,
                ),
            })),
            AdaptiveSearchOutcome::NoRecurrence(summary) => Ok(Json(recurrence_failure_response(
                RecurrenceSearchStatus::NotFound,
                "no recurrence found within the search bounds",
                Some(summary),
                Vec::new(),
            ))),
            AdaptiveSearchOutcome::BudgetExhausted(summary) => {
                let error = format!(
                    "recurrence candidate budget exhausted after {} candidates; unsearched candidates remain",
                    summary.diagnostics.considered_candidates
                );
                Ok(Json(recurrence_failure_response(
                    RecurrenceSearchStatus::BudgetExhausted,
                    error,
                    Some(summary),
                    Vec::new(),
                )))
            }
        }
    }

    #[tool(
        description = "Find several recurrence representations trading lag length against derivative order. Uses the same inputs and bounds as find_recurrence, but keeps one preferred exact fit per length/order bound pair. Returns polytool.recurrence-representations.v1: candidates with reusable recurrence_json objects, actual complexity, checked starting indices, finite-fit evidence, a nondominated frontier, and per-pair coverage. The candidate budget is global; incomplete coverage is explicit. No automatic denominator escalation. These are finite-data candidates, not all-rank certificates, a complete annihilator basis, or nonexistence proofs. include_code is accepted for request compatibility; recurrence_json is always included and language code is not generated."
    )]
    pub fn find_recurrence_representations(
        &self,
        Parameters(input): Parameters<FindRecurrenceRequest>,
    ) -> Result<Json<BTreeMap<String, serde_json::Value>>, McpError> {
        let (search, budget) = apply_recurrence_options(input.effective_options())?;
        let batch = parse_recurrence_batch_rational(&input)?;
        let polynomials = collect_rational_polynomials_or_errors(batch)
            .map_err(|_| invalid_params("one or more polynomials failed to parse"))?;
        let report = find_recurrence_representations_rational(&polynomials, &search, budget)
            .map_err(invalid_params)?;
        let serde_json::Value::Object(fields) =
            serde_json::to_value(report).expect("serialize recurrence menu")
        else {
            unreachable!("recurrence report is an object")
        };
        Ok(Json(fields.into_iter().collect()))
    }

    #[tool(
        description = "Compute complete bounded recurrence spaces with exact kernel/rank certificates. request is a JSON object with shapes (kind: box with length, derivative_order, index_degree, variable_degree, leading_index_degree; or kind: explicit with target_lag and monomials), domain_start (original zero-based index by default), and source {coefficients: ascending rational-string polynomial rows, seeds: exactly k rows, origin}, source_recurrence (existing recurrence.v1), or rows plus row_origin. Alternatively indexed_source uses coefficients[i][n_power][X_power] for monic pure-lag sources polynomial in n, with a replayable jet-state propagator. Shape kind automatic_variable_degree removes the RHS X-degree cap while retaining index/leader bounds. include_leader_coverage returns a sum-of-squares leader combination. Optional coefficient_rescaling {kind: factorial_power, power: m >= 1} (constant-index source only, no automatic shapes) searches Q[n] = T(P[n]) with T(X^k) = X^k/(k!)^m; each result then carries rescaling_bound, certifying the complete space from r(U+m(A+J)+1) residuals. Positive identities are conditional_on_source; rows-only positives are exact_fit. Dual negatives prove absent_in_shape from rows alone. Returns polytool.recurrence-spaces.v1, including full kernels, successors, bounded frontier, singular-rank diagnostics and explicit incomplete results. No spectral or unrestricted-order claim."
    )]
    pub fn find_recurrence_spaces(
        &self,
        Parameters(input): Parameters<RecurrenceSpacesToolRequest>,
    ) -> Result<Json<BTreeMap<String, serde_json::Value>>, McpError> {
        validate_mcp_text(&input.request.to_string(), "recurrence-space request")
            .map_err(invalid_params)?;
        let request: polytool::recurrence::RecurrenceSpaceRequest =
            serde_json::from_value(input.request).map_err(|e| invalid_params(e.to_string()))?;
        request
            .limits
            .check_ceiling(&Default::default())
            .map_err(invalid_params)?;
        let report =
            polytool::recurrence::find_recurrence_spaces(&request).map_err(invalid_params)?;
        let serde_json::Value::Object(fields) =
            serde_json::to_value(report).expect("serialize recurrence spaces")
        else {
            unreachable!()
        };
        Ok(Json(fields.into_iter().collect()))
    }

    #[tool(
        description = "Replay exact recurrence-space kernel, rank and dual certificates without rerunning the full search. Supply report from find_recurrence_spaces and optionally expected_request to bind independently supplied assumptions. Without expected_request this verifies the report's internal source-conditional claims, not the truth of a guessed source. Incomplete shapes remain explicitly incomplete."
    )]
    pub fn verify_recurrence_spaces(
        &self,
        Parameters(input): Parameters<VerifyRecurrenceSpacesToolRequest>,
    ) -> Result<Json<BTreeMap<String, serde_json::Value>>, McpError> {
        validate_mcp_text(&input.report.to_string(), "recurrence-space report")
            .map_err(invalid_params)?;
        let report: polytool::recurrence::RecurrenceSpaces =
            serde_json::from_value(input.report).map_err(|e| invalid_params(e.to_string()))?;
        let request = match input.expected_request {
            Some(value) => {
                validate_mcp_text(&value.to_string(), "expected request")
                    .map_err(invalid_params)?;
                serde_json::from_value(value).map_err(|e| invalid_params(e.to_string()))?
            }
            None => report.request.clone(),
        };
        request
            .limits
            .check_ceiling(&Default::default())
            .map_err(invalid_params)?;
        polytool::recurrence::verify_recurrence_spaces(&request, &report)
            .map_err(invalid_params)?;
        Ok(Json(BTreeMap::from([
            ("verified_completed_claims".into(), json!(true)),
            (
                "incomplete_shapes".into(),
                json!(report
                    .results
                    .iter()
                    .filter(|r| r.certificate.is_none())
                    .count()),
            ),
        ])))
    }

    #[tool(
        description = "Apply an exact change of polynomial coordinates and transport a defining recurrence. request: {transformation: {kind: multiply_power {power} | divide_power {power} (needs x^r | P_n) | reversal {schedule} (x^(d_n) P_n(1/x)) | affine {a, b} (P_n(a x + b)) | mobius {a, b, c, e, schedule} ((c x + e)^(d_n) P_n((a x + b)/(c x + e)), ae - bc != 0) | f_polynomial {schedule} | gamma {schedule} (symmetry index d_n) | symmetric_decomposition {schedule} (P_n = a_n + x b_n)}, rationals as strings; schedule {first_index, first_degree, step, last_index?} gives d_n = first_degree + step*(n - first_index) on an explicit index domain. Source: source_recurrence (polytool.recurrence.v1 with seeds), oeis_id (bundled family) or rows with row_origin; row_count rows are generated from a source. Optional search ({} for defaults) fits independently in the new coordinates. Returns polytool.recurrence-transform.v1: exact row map, inverse and image conditions, transformed states, exact reconstruction, condition statuses (certified_all_indices from the source rule, or checked_supplied_rows), the transported rule (evidence transported_identity, with transformed seeds) and the independent fit (evidence finite_fit), kept separately. Gamma coordinates P_n = (1+x)^(d_n) G_n(x/(1+x)^2) transport a rule when the symmetric (even) and antisymmetric (odd) parts of the cleared identity are proportional, otherwise the transport error names both parts. Optional root_domain ({} for defaults, or {candidate_intervals, assumptions, default_window, max_candidate_intervals, max_observed_rows, max_observed_degree}; interval {lower?, upper?} with endpoints {value, per_index?, closed?}, an omitted end is infinite; assumption {subject: source|transformed, component?, property: real_rooted|real_rooted_in|real_roots_in, interval?, first_index?, last_index?, note?}) attaches root implications (direction, hypotheses with status, converse), claims labelled finite_exact_observation, supplied_assumption, proposed_interval, completed_recurrence_induction or certified_family_invariant, the Rolle recurrence match Q_m = s(m) D^j(g Q_(m-1)) (also g(x)(c(m) + l(m) x) and index-dependent intervals, checked symbolically; gamma transfers only inside (-inf, 1/4]) with its interval attempts, and exact root enclosures; the algebraic results never depend on it."
    )]
    pub fn transform_recurrence(
        &self,
        Parameters(input): Parameters<RecurrenceTransformToolRequest>,
    ) -> Result<Json<BTreeMap<String, serde_json::Value>>, McpError> {
        validate_mcp_text(&input.request.to_string(), "recurrence-transform request")
            .map_err(invalid_params)?;
        let request: polytool::recurrence::RecurrenceTransformRequest =
            serde_json::from_value(input.request).map_err(|e| invalid_params(e.to_string()))?;
        request
            .limits
            .check_ceiling(&Default::default())
            .map_err(invalid_params)?;
        let report =
            polytool::recurrence::transform_recurrence(&request).map_err(invalid_params)?;
        let serde_json::Value::Object(fields) =
            serde_json::to_value(report).expect("serialize recurrence transform")
        else {
            unreachable!()
        };
        Ok(Json(fields.into_iter().collect()))
    }

    #[tool(
        description = "Replay a polytool.recurrence-transform.v1 report from transform_recurrence exactly: recompute the map, reconstruction, transported rule and independent search from the embedded request and require equality. The transported rule remains conditional on the supplied source rule."
    )]
    pub fn verify_recurrence_transform(
        &self,
        Parameters(input): Parameters<VerifyRecurrenceTransformToolRequest>,
    ) -> Result<Json<BTreeMap<String, serde_json::Value>>, McpError> {
        validate_mcp_text(&input.report.to_string(), "recurrence-transform report")
            .map_err(invalid_params)?;
        let report: polytool::recurrence::RecurrenceTransformReport =
            serde_json::from_value(input.report).map_err(|e| invalid_params(e.to_string()))?;
        report
            .request
            .limits
            .check_ceiling(&Default::default())
            .map_err(invalid_params)?;
        polytool::recurrence::verify_recurrence_transform(&report).map_err(invalid_params)?;
        Ok(Json(BTreeMap::from([
            ("verified".into(), json!(true)),
            ("transported".into(), json!(report.transported.is_some())),
            (
                "independent_search".into(),
                json!(report.independent_search.is_some()),
            ),
        ])))
    }

    #[tool(
        description = "Exact real-root enclosures with multiplicities. request: {polynomials: [[ascending rational strings]], max_width?: rational (default refinement width 1/1024), interval?: {lower?, upper?} with endpoints {value, per_index?, closed?} (omitted end infinite; value + per_index * n with n = first_index + position), first_index?}. Returns polytool.root-enclosures.v1 items with degree, valuation, real_rooted, real_root_count (with multiplicity), enclosures {lower, upper, multiplicity} (lower == upper is an exact rational root; otherwise exactly one root lies strictly inside and neither endpoint is a root) and exact interval verdicts real_roots_in_interval and real_rooted_in_interval. Each verdict is exact for the given polynomial only."
    )]
    pub fn exact_root_enclosures(
        &self,
        Parameters(input): Parameters<RootEnclosuresToolRequest>,
    ) -> Result<Json<BTreeMap<String, serde_json::Value>>, McpError> {
        validate_mcp_text(&input.request.to_string(), "root-enclosures request")
            .map_err(invalid_params)?;
        let request: polytool::recurrence::RootEnclosureRequest =
            serde_json::from_value(input.request).map_err(|e| invalid_params(e.to_string()))?;
        let report = polytool::recurrence::root_enclosures(&request).map_err(invalid_params)?;
        let serde_json::Value::Object(fields) =
            serde_json::to_value(report).expect("serialize root enclosures")
        else {
            unreachable!()
        };
        Ok(Json(fields.into_iter().collect()))
    }

    #[tool(
        description = "Map a sequence with a defining rule by a diagonal coefficient map or a derivative antecedent and transport the rule into an implicit differential relation (derivatives on any row). request: {transformation: {kind: diagonal, weight: {kind: factorial_power {power} (lambda = 1/(k!)^m) | binomial {m} (lambda = binom(k+m,m), Q = W_m P) | inverse_binomial {m} (P = W_m Q) | product {scale, factorial_exponent, factors: [{coeffs: {coeffs}, exponent}]}} | {kind: derivative_antecedent, g: [..]} (P_n = D^m(g Q_n)/m!, m = deg g)}, source_recurrence (polytool.recurrence.v1 with seeds) or oeis_id, row_count, optional source_relations / target_relations [{relation: {terms: [{shift, deriv_order, coeff: {coeffs}}]}, first_base}], boundary_conditions [{condition: {terms: [{shift, coefficient_index, coeff}]}, first_base}], search ({} fits a successor rule and tries a reconstruction certificate). Returns polytool.recurrence-implicit.v1: map conditions (certified_all_indices or checked_supplied_rows), transported relations with clearing certificates (one common multiplier H, certified exceptional indices and their coefficient conditions; uncertified sets stay visible), reconstruction certificates L o T = sum H_nu S^nu M, a generating description with an assessment (equivalent_description, certified_consequence, finite_evidence, inconsistent_seeds) whose uniqueness argument is separate from annihilation, additional row relations, and the search (finite_fit or certified_by_reconstruction). Optional root_domain (as in transform_recurrence) adds root_domain: multiplier-sequence (diagonal; Laguerre f(k) with real roots <= 0, 1/k!, products) or Rolle (derivative antecedent, Q => P only) adapters with direction and exact hypotheses, Rolle inductions on the defining rule and on a certified fit, transfers only in the certified direction, assumptions, a jensen_fixed_cap test (fixed_cap evidence only) and exact root enclosures."
    )]
    pub fn transform_recurrence_implicit(
        &self,
        Parameters(input): Parameters<RecurrenceImplicitToolRequest>,
    ) -> Result<Json<BTreeMap<String, serde_json::Value>>, McpError> {
        validate_mcp_text(&input.request.to_string(), "recurrence-implicit request")
            .map_err(invalid_params)?;
        let (request, root_domain) =
            polytool::recurrence::parse_implicit_transform_input(input.request)
                .map_err(invalid_params)?;
        request
            .limits
            .check_ceiling(&Default::default())
            .map_err(invalid_params)?;
        let report = polytool::recurrence::implicit_transform_recurrence_with_root_domain(
            &request,
            root_domain.as_ref(),
        )
        .map_err(invalid_params)?;
        let serde_json::Value::Object(fields) =
            serde_json::to_value(report).expect("serialize implicit transport")
        else {
            unreachable!()
        };
        Ok(Json(fields.into_iter().collect()))
    }

    #[tool(
        description = "Replay a polytool.recurrence-implicit.v1 report from transform_recurrence_implicit exactly from its embedded request, re-check every reconstruction certificate symbolically and by deriving its identities on the transformed rows, and re-assess the description."
    )]
    pub fn verify_recurrence_implicit(
        &self,
        Parameters(input): Parameters<VerifyRecurrenceImplicitToolRequest>,
    ) -> Result<Json<BTreeMap<String, serde_json::Value>>, McpError> {
        validate_mcp_text(&input.report.to_string(), "recurrence-implicit report")
            .map_err(invalid_params)?;
        let report: polytool::recurrence::ImplicitTransformReport =
            serde_json::from_value(input.report).map_err(|e| invalid_params(e.to_string()))?;
        report
            .request
            .limits
            .check_ceiling(&Default::default())
            .map_err(invalid_params)?;
        polytool::recurrence::verify_recurrence_implicit(&report).map_err(invalid_params)?;
        Ok(Json(BTreeMap::from([
            ("verified".into(), json!(true)),
            ("assessment".into(), json!(report.assessment.status)),
            (
                "independent_search".into(),
                json!(report.independent_search.map(|s| s.evidence)),
            ),
        ])))
    }

    #[tool(
        description = "Assess a generating description {first_index, seeds, relations: [{relation, first_base, role, evidence, text}], exceptional_index_conditions, boundary_conditions} (for example one from transform_recurrence_implicit with a condition removed): seed consistency, annihilation evidence and a uniqueness argument, reported separately."
    )]
    pub fn assess_generating_description(
        &self,
        Parameters(input): Parameters<AssessGeneratingDescriptionToolRequest>,
    ) -> Result<Json<BTreeMap<String, serde_json::Value>>, McpError> {
        validate_mcp_text(&input.description.to_string(), "generating description")
            .map_err(invalid_params)?;
        let description: polytool::recurrence::GeneratingDescription =
            serde_json::from_value(input.description).map_err(|e| invalid_params(e.to_string()))?;
        if description.seeds.len() > crate::limits::MCP_MAX_RECURRENCE_INPUT_ROWS {
            return Err(invalid_params("too many seed rows".to_string()));
        }
        let assessment = polytool::recurrence::assess_generating_description(&description)
            .map_err(invalid_params)?;
        let serde_json::Value::Object(fields) =
            serde_json::to_value(assessment).expect("serialize assessment")
        else {
            unreachable!()
        };
        Ok(Json(fields.into_iter().collect()))
    }

    #[tool(
        description = "Generate or extend coefficient rows from recurrence JSON returned by find_recurrence or one recurrence_json object from find_recurrence_representations."
    )]
    pub fn generate_recurrence_rows(
        &self,
        Parameters(input): Parameters<GenerateRecurrenceRowsRequest>,
    ) -> Result<Json<GenerateRecurrenceRowsResponse>, McpError> {
        if input.rows.is_some() == input.additional.is_some() {
            return Err(invalid_params("provide exactly one of rows or additional"));
        }
        validate_mcp_text(&input.recurrence_json, "recurrence JSON").map_err(invalid_params)?;
        let recurrence_json: RecurrenceJson = serde_json::from_str(&input.recurrence_json)
            .map_err(|error| invalid_params(format!("failed to parse recurrence JSON: {error}")))?;
        let (recurrence, first_index, initial_polys) = recurrence_json
            .to_recurrence_parts()
            .map_err(|error| invalid_params(format!("invalid recurrence JSON: {error}")))?;
        let row_count = match (input.rows, input.additional) {
            (Some(rows), None) => rows,
            (None, Some(additional)) => initial_polys
                .len()
                .checked_add(additional)
                .ok_or_else(|| invalid_params("recurrence row count overflow"))?,
            _ => unreachable!("validated exactly one row-count option"),
        };
        if row_count > MCP_MAX_RECURRENCE_ROWS {
            return Err(invalid_params(format!(
                "requested {row_count} recurrence rows; the MCP limit is {MCP_MAX_RECURRENCE_ROWS}"
            )));
        }
        first_index
            .checked_add(row_count)
            .ok_or_else(|| invalid_params("recurrence index range overflow"))?;
        validate_recurrence_for_generation(&recurrence, &initial_polys)?;
        let generated = recurrence
            .generate_rows_rational(&initial_polys, first_index, row_count)
            .map_err(|error| invalid_params(format!("failed to generate rows: {error}")))?;
        let polynomials = generated
            .iter()
            .map(|row| row.iter().map(format_rational_coeff).collect())
            .collect::<Vec<Vec<_>>>();
        Ok(Json(GenerateRecurrenceRowsResponse {
            first_index,
            row_count: generated.len(),
            polynomials,
        }))
    }

    #[tool(
        description = "Check a supplied recurrence JSON record (polytool.recurrence.v1, e.g. recurrence_json from find_recurrence) against the supported interlacing theorems: weighted_renewal_window (P_n = sum_j r_j t^j P_(n-1-ell*j), real-rooted nonnegative seed, window d*s <= d*ell+1) and fixed_delay_predecessor_window (P_n = P_(n-1) + w(n) t P_(n-L), polynomial weight certified nonnegative on all integers n >= 1, window B - A <= L). Returns polytool.recurrence-theorems.v1 with finite_evidence, theorem_matches (exact hypotheses, conclusion, window, boundary check of the record's initial rows, remaining obligations, proof_status paper_proof_not_lean_checked) and rejected_theorems with reasons. Matches are conditional: the caller must still establish the recurrence for the intended sequence for all n."
    )]
    pub fn recurrence_theorems(
        &self,
        Parameters(input): Parameters<crate::requests::RecurrenceTheoremsRequest>,
    ) -> Result<Json<BTreeMap<String, serde_json::Value>>, McpError> {
        crate::theorems::supplied_theorem_report(&input.recurrence_json).map(Json)
    }

    #[tool(description = "Compute the exact resultant of two polynomials.")]
    pub fn resultant(
        &self,
        Parameters(input): Parameters<BigIntPolynomialPairRequest>,
    ) -> Result<Json<ResultantResponse>, McpError> {
        let p = normalize_bigint_polynomial(parse_required_bigint_polynomial(&input.p, "p")?);
        let q = normalize_bigint_polynomial(parse_required_bigint_polynomial(&input.q, "q")?);
        let resultant =
            polytool::resultant_bigint_coeffs(&p.coefficients, &q.coefficients).to_string();
        Ok(Json(ResultantResponse {
            p: p.display,
            q: q.display,
            resultant,
        }))
    }

    #[tool(description = "Compute exact discriminants for one or more polynomials.")]
    pub fn discriminant(
        &self,
        Parameters(input): Parameters<BigIntPolynomialBatchInput>,
    ) -> Result<Json<DiscriminantResponse>, McpError> {
        let batch = parse_bigint_batch(&input)?;
        let items = batch
            .into_iter()
            .enumerate()
            .map(|(index, item)| match item {
                Ok(polynomial) => DiscriminantItem {
                    index,
                    ok: true,
                    error: None,
                    discriminant: Some(
                        polytool::discriminant_bigint_coeffs(&polynomial.coefficients).to_string(),
                    ),
                    polynomial: Some(polynomial.display),
                },
                Err(error) => DiscriminantItem {
                    index,
                    ok: false,
                    error: Some(error),
                    polynomial: None,
                    discriminant: None,
                },
            })
            .collect();
        Ok(Json(DiscriminantResponse { items }))
    }

    #[tool(
        description = "Convert h*-vectors to Ehrhart polynomials or Ehrhart polynomials to h*-vectors."
    )]
    pub fn ehrhart_hstar(
        &self,
        Parameters(input): Parameters<EhrhartHstarRequest>,
    ) -> Result<Json<EhrhartHstarResponse>, McpError> {
        match input.mode {
            EhrhartHstarMode::HstarToEhrhart => {
                let hstar_input = input
                    .hstar
                    .ok_or_else(|| invalid_params("`hstar` is required for hstar_to_ehrhart"))?;
                let hstar = parse_bigint_coefficients(&hstar_input)?;
                let ehrhart = hstar_to_ehrhart_bigint_coeffs(&hstar)
                    .into_iter()
                    .map(format_rational)
                    .collect();
                Ok(Json(EhrhartHstarResponse {
                    mode: EhrhartHstarMode::HstarToEhrhart,
                    hstar: Some(bigint_strings(&hstar)),
                    ehrhart_coefficients: Some(ehrhart),
                }))
            }
            EhrhartHstarMode::EhrhartToHstar => {
                let hstar = match (
                    input.ehrhart_coefficients,
                    input.numerator_coefficients,
                    input.denominator,
                ) {
                    (Some(coefficients), None, None) => {
                        if coefficients.len() > MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL {
                            return Err(invalid_params(format!(
                                "Ehrhart polynomial has {} coefficients; the MCP limit is {MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL}",
                                coefficients.len()
                            )));
                        }
                        for coefficient in &coefficients {
                            if coefficient.len() > MCP_MAX_COEFFICIENT_TEXT_BYTES {
                                return Err(invalid_params(format!(
                                    "Ehrhart coefficient is {} bytes; the MCP limit is {MCP_MAX_COEFFICIENT_TEXT_BYTES}",
                                    coefficient.len()
                                )));
                            }
                        }
                        let coefficients: Result<Vec<_>, _> =
                            coefficients.iter().map(|c| parse_rational(c)).collect();
                        ehrhart_to_hstar_bigint(&coefficients.map_err(invalid_params)?)
                            .map_err(|error| invalid_params(error.to_string()))?
                    }
                    (None, Some(numerator_coefficients), Some(denominator)) => {
                        let numerator_coefficients =
                            parse_bigint_coefficients(&numerator_coefficients)?;
                        let denominator =
                            parse_bigint_coefficient(&denominator).map_err(invalid_params)?;
                        if denominator == BigInt::from(0) {
                            return Err(invalid_params("`denominator` must be nonzero"));
                        }
                        let coefficients = numerator_coefficients
                            .into_iter()
                            .map(|numerator| BigRational::new(numerator, denominator.clone()))
                            .collect::<Vec<_>>();
                        ehrhart_to_hstar_bigint(&coefficients)
                            .map_err(|error| invalid_params(error.to_string()))?
                    }
                    _ => {
                        return Err(invalid_params(
                            "for ehrhart_to_hstar, provide either `ehrhart_coefficients` or both `numerator_coefficients` and `denominator`",
                        ));
                    }
                };
                Ok(Json(EhrhartHstarResponse {
                    mode: EhrhartHstarMode::EhrhartToHstar,
                    hstar: Some(bigint_strings(&hstar)),
                    ehrhart_coefficients: None,
                }))
            }
        }
    }

    #[tool(
        description = "Profile or check a cyclic sieving polynomial for one cyclic group order using exact cyclotomic arithmetic."
    )]
    pub fn cyclic_sieving(
        &self,
        Parameters(input): Parameters<CyclicSievingRequest>,
    ) -> Result<Json<CyclicSievingResponse>, McpError> {
        if input.order == 0 {
            return Err(invalid_params("`order` must be positive"));
        }
        if input.order > MCP_MAX_CYCLIC_ORDER {
            return Err(invalid_params(format!(
                "`order` {} exceeds the MCP limit {MCP_MAX_CYCLIC_ORDER}",
                input.order
            )));
        }
        if input.fixed_counts.is_some() && input.orbit_sizes.is_some() {
            return Err(invalid_params(
                "`fixed_counts` and `orbit_sizes` are alternatives; give only one",
            ));
        }
        if input
            .orbit_sizes
            .as_ref()
            .is_some_and(|sizes| sizes.len() > MCP_MAX_TOTAL_COEFFICIENTS)
        {
            return Err(invalid_params(format!(
                "`orbit_sizes` may have at most {MCP_MAX_TOTAL_COEFFICIENTS} entries"
            )));
        }
        let polynomial = normalize_bigint_polynomial(parse_required_bigint_polynomial(
            &input.polynomial,
            "polynomial",
        )?);
        let orbit_sizes_check = input
            .orbit_sizes
            .as_ref()
            .map(|sizes| {
                cyclic_sieving_check_orbit_sizes(&polynomial.coefficients, input.order, sizes)
            })
            .transpose()
            .map_err(invalid_params)?;
        let fixed_counts = match &orbit_sizes_check {
            Some(check) => Some(check.actual_fixed_points.clone()),
            None => input
                .fixed_counts
                .as_ref()
                .map(|counts| parse_bigint_coefficients(counts))
                .transpose()?,
        };
        let report = cyclic_sieving_report_bigint(
            &polynomial.coefficients,
            input.order,
            fixed_counts.as_deref(),
        )
        .map_err(invalid_params)?;
        let decomposition =
            cyclic_sieving_orbit_decomposition_bigint(&polynomial.coefficients, input.order)
                .map_err(invalid_params)?;
        Ok(Json(CyclicSievingResponse {
            report: cyclic_sieving_report_item(&report),
            orbit_decomposition: orbit_decomposition_item(&decomposition),
            orbit_sizes_check: orbit_sizes_check.as_ref().map(orbit_sizes_check_item),
        }))
    }

    #[tool(
        description = "Check bicyclic sieving for G = C_n1 x C_n2 and f(q,t) given as a coefficient matrix (matrix[i][j] = coeff of q^i t^j). Reduces modulo (q^n1-1, t^n2-1), evaluates exactly at all (zeta_n1^d1, zeta_n2^d2), reports integrality, and, when c_(j1,j2)=c_(gcd(j1,n1),gcd(j2,n2)), the coordinates in the product basis h_r1(q)h_r2(t) (not orbit counts; may be negative: 1+qt gives 1,-1,-1,2). With fixed_table (n1 x n2 table of fix(g1^d1 g2^d2)) it checks the CSP equalities."
    )]
    pub fn bicyclic_sieving(
        &self,
        Parameters(input): Parameters<BicyclicSievingRequest>,
    ) -> Result<Json<BicyclicSievingResponse>, McpError> {
        if input.order1 == 0 || input.order2 == 0 {
            return Err(invalid_params("`order1` and `order2` must be positive"));
        }
        let matrix = parse_bigint_matrix(&input.matrix, "matrix")?;
        let table = input
            .fixed_table
            .as_ref()
            .map(|table| parse_bigint_matrix(table, "fixed_table"))
            .transpose()?;
        let report =
            bicyclic_sieving_report_bigint(&matrix, input.order1, input.order2, table.as_deref())
                .map_err(invalid_params)?;
        Ok(Json(BicyclicSievingResponse {
            order1: report.order1,
            order2: report.order2,
            reduced_matrix: bigint_matrix_strings(&report.reduced_matrix),
            evaluations: report
                .evaluations
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|cell| cell.as_ref().map(ToString::to_string))
                        .collect()
                })
                .collect(),
            integral: report.integral,
            natural_coefficients: report.natural_coefficients,
            divisors1: report.divisors1.clone(),
            divisors2: report.divisors2.clone(),
            product_basis_form: report.product_basis_form,
            product_coordinates: report
                .product_coordinates
                .as_deref()
                .map(bigint_matrix_strings),
            product_formula_evaluations: report
                .product_formula_evaluations
                .as_deref()
                .map(bigint_matrix_strings),
            fixed_table: report.fixed_table.as_deref().map(bigint_matrix_strings),
            mismatches: report.mismatches.iter().map(|&(a, b)| [a, b]).collect(),
            holds: report.holds,
        }))
    }

    #[tool(
        description = "For a polynomial sequence P_n(q), profile candidate cyclic sieving group orders n-2,n-1,n,n+1,n+2,n+3 by default."
    )]
    pub fn cyclic_sieving_sequence(
        &self,
        Parameters(input): Parameters<CyclicSievingSequenceRequest>,
    ) -> Result<Json<CyclicSievingSequenceResponse>, McpError> {
        let first_index = input.first_index.unwrap_or(0);
        let offsets = input.offsets.unwrap_or_else(|| vec![-2, -1, 0, 1, 2, 3]);
        if offsets.is_empty() {
            return Err(invalid_params("`offsets` may not be empty"));
        }
        let batch = parse_bigint_batch(&BigIntPolynomialBatchInput {
            polynomials: input.polynomials,
            text: input.text,
        })?;
        let mut polynomials = Vec::new();
        for (index, item) in batch.into_iter().enumerate() {
            match item {
                Ok(polynomial) => polynomials.push(polynomial.coefficients),
                Err(error) => {
                    return Err(invalid_params(format!("polynomial {index}: {error}")));
                }
            }
        }
        let mut fixed_counts = BTreeMap::new();
        for item in input.fixed_counts.unwrap_or_default() {
            if item.order == 0 {
                return Err(invalid_params("fixed-count order must be positive"));
            }
            let counts = parse_bigint_coefficients(&item.counts)?;
            if counts.len() != item.order {
                return Err(invalid_params(format!(
                    "fixed counts for index {} order {}: expected {} counts, got {}",
                    item.index,
                    item.order,
                    item.order,
                    counts.len()
                )));
            }
            fixed_counts.insert((item.index, item.order), counts);
        }
        let reports = cyclic_sieving_sequence_reports_bigint(
            &polynomials,
            first_index,
            &offsets,
            &fixed_counts,
        );
        Ok(Json(CyclicSievingSequenceResponse {
            first_index,
            offsets,
            items: reports.iter().map(cyclic_sequence_item_response).collect(),
        }))
    }

    #[tool(
        description = "Analyze symmetric decomposition, R-transform, interlacing checks, and magic-basis coordinates."
    )]
    pub fn analyze_decomposition(
        &self,
        Parameters(input): Parameters<PolynomialBatchInput>,
    ) -> Result<Json<DecompositionResponse>, McpError> {
        let batch = parse_batch(&input)?;
        let items = batch
            .into_iter()
            .enumerate()
            .map(|(index, item)| match item {
                Ok(polynomial) => {
                    match analyze_symmetric_decomposition_i64(&polynomial.coefficients) {
                        Ok(analysis) => {
                            let partial_sum_checks = analysis
                                .magic
                                .left_partial_sums
                                .iter()
                                .zip(analysis.magic.right_partial_sums.iter())
                                .map(|(left, right)| left <= right)
                                .collect();
                            DecompositionItem {
                                index,
                                ok: true,
                                report: Some(DecompositionReport {
                                    polynomial,
                                    reciprocal: displayed(&analysis.reciprocal),
                                    a: displayed(&analysis.a),
                                    b: displayed(&analysis.b),
                                    a_real_rooted: analysis.a_real_rooted,
                                    b_real_rooted: analysis.b_real_rooted,
                                    b_interlaces_a: analysis.b_interlaces_a,
                                    reciprocal_interlaces_input: analysis
                                        .reciprocal_interlaces_input,
                                    alternatingly_increasing: analysis.alternatingly_increasing,
                                    f_polynomial: displayed(&analysis.f_polynomial),
                                    r_transform_of_f: displayed(&analysis.r_transform_of_f),
                                    r_a: displayed(&analysis.r_a),
                                    r_b: displayed(&analysis.r_b),
                                    r_interlaces_f: analysis.r_interlaces_f,
                                    magic: MagicBasisReport {
                                        coordinates: analysis
                                            .magic
                                            .coordinates
                                            .iter()
                                            .map(ToString::to_string)
                                            .collect(),
                                        left_partial_sums: analysis
                                            .magic
                                            .left_partial_sums
                                            .iter()
                                            .map(ToString::to_string)
                                            .collect(),
                                        right_partial_sums: analysis
                                            .magic
                                            .right_partial_sums
                                            .iter()
                                            .map(ToString::to_string)
                                            .collect(),
                                        partial_sum_checks,
                                        all_nonnegative: analysis.magic.all_nonnegative,
                                        left_leq_right: analysis.magic.left_leq_right,
                                    },
                                }),
                                error: None,
                            }
                        }
                        Err(error) => DecompositionItem {
                            index,
                            ok: false,
                            report: None,
                            error: Some(error.to_string()),
                        },
                    }
                }
                Err(error) => DecompositionItem {
                    index,
                    ok: false,
                    report: None,
                    error: Some(error),
                },
            })
            .collect();
        Ok(Json(DecompositionResponse { items }))
    }

    #[tool(description = "List bundled recurrence-backed OEIS polynomial families.")]
    pub fn list_oeis_sequences(
        &self,
        Parameters(input): Parameters<ListOeisSequencesRequest>,
    ) -> Result<Json<ListOeisSequencesResponse>, McpError> {
        let include_experimental = input.include_experimental.unwrap_or(false);
        let sequences = polytool::oeis::catalog()
            .iter()
            .filter(|entry| {
                include_experimental
                    || entry.status != polytool::oeis::OeisSequenceStatus::Experimental
            })
            .map(oeis_sequence_summary)
            .collect::<Vec<_>>();
        Ok(Json(ListOeisSequencesResponse {
            count: sequences.len(),
            sequences,
        }))
    }

    #[tool(
        description = "Get metadata and exact recurrence exports for one bundled OEIS polynomial family."
    )]
    pub fn get_oeis_sequence(
        &self,
        Parameters(input): Parameters<GetOeisSequenceRequest>,
    ) -> Result<Json<GetOeisSequenceResponse>, McpError> {
        let entry = polytool::oeis::by_id(&input.id).ok_or_else(|| {
            invalid_params(format!("unknown bundled OEIS sequence: {}", input.id))
        })?;
        let (recurrence, initial_rows, first_index) = entry
            .recurrence_export_parts()
            .map_err(|error| invalid_params(error.to_string()))?;
        Ok(Json(GetOeisSequenceResponse {
            sequence: oeis_sequence_summary(entry),
            recurrence: recurrence.to_string(),
            latex: recurrence.to_latex(),
            mathematica: recurrence
                .to_mathematica_definition_rational_indexed(&initial_rows, first_index),
            sage: recurrence.to_sage_definition_rational_indexed(&initial_rows, first_index),
            python: recurrence.to_python_definition_rational_indexed(&initial_rows, first_index),
        }))
    }

    #[tool(
        description = "Generate exact polynomial coefficient rows for one bundled recurrence-backed OEIS family."
    )]
    pub fn generate_oeis_rows(
        &self,
        Parameters(input): Parameters<GenerateOeisRowsRequest>,
    ) -> Result<Json<GenerateOeisRowsResponse>, McpError> {
        const MAX_ROWS: usize = 200;
        if input.rows > MAX_ROWS {
            return Err(invalid_params(format!(
                "rows must be at most {MAX_ROWS} for MCP generation"
            )));
        }
        let entry = polytool::oeis::by_id(&input.id).ok_or_else(|| {
            invalid_params(format!("unknown bundled OEIS sequence: {}", input.id))
        })?;
        if entry.status == polytool::oeis::OeisSequenceStatus::Experimental
            && !input.include_experimental.unwrap_or(false)
        {
            return Err(invalid_params(format!(
                "{} is experimental; set include_experimental to true",
                entry.id
            )));
        }
        let first_row = input.first_row.unwrap_or(entry.first_row);
        let generated = entry
            .generate_rows_from(first_row, input.rows)
            .map_err(|error| invalid_params(error.to_string()))?;
        let rows = generated
            .into_iter()
            .enumerate()
            .map(|(offset, coefficients)| {
                let offset = i64::try_from(offset)
                    .map_err(|_| invalid_params("OEIS row offset does not fit in i64"))?;
                let n = first_row
                    .checked_add(offset)
                    .ok_or_else(|| invalid_params("OEIS row index overflow"))?;
                Ok(OeisPolynomialRow {
                    n,
                    polynomial: format_poly_bigint_coeffs(&coefficients),
                    coefficients: coefficients.iter().map(ToString::to_string).collect(),
                })
            })
            .collect::<Result<Vec<_>, McpError>>()?;
        Ok(Json(GenerateOeisRowsResponse {
            id: entry.id.to_string(),
            first_row,
            row_count: rows.len(),
            rows,
        }))
    }

    #[tool(description = "Generate standard polynomial sequences.")]
    pub fn generate_sequence(
        &self,
        Parameters(input): Parameters<GenerateSequenceRequest>,
    ) -> Result<Json<GenerateSequenceResponse>, McpError> {
        if input.max_n > MCP_MAX_SEQUENCE_N {
            return Err(invalid_params(format!(
                "max_n must be at most {MCP_MAX_SEQUENCE_N} for MCP sequence generation"
            )));
        }
        let polynomials = generated_sequence_polynomials_bigint(&input.sequence, input.max_n)
            .into_iter()
            .map(|coefficients| normalize_bigint_polynomial(coefficients).display)
            .collect();

        Ok(Json(GenerateSequenceResponse {
            sequence: input.sequence,
            max_n: input.max_n,
            polynomials,
        }))
    }

    #[tool(
        description = "Decide exactly whether the coefficient triangle r_{n,k} = [t^k] P_n of monic rows P_0 (or P_1), P_1, ... with deg P_n = n (or given `matrix` rows, square lower-triangular or triangular) is lower-unitriangular and totally nonnegative (TNN). On success returns the canonical planar network weights lambda[n][k] (vertical edge (n+1,k) -> (n,k)) and whether the network's path matrix re-evaluates to the input; on failure the reason (shape, negative entry, or an explicit negative minor with zero-based rows/cols and exact value). With chain=true also returns the chain polynomials p_0 = 1, p_n = t sum_{k<n} r_{n,k} p_k with checks of real-rootedness in [-1,0] and consecutive interlacing (Branden--Saud Maia Leite, arXiv:2412.06595, Thm 3.7: these hold when the matrix is TNN). Guardrail: TNN of the matrix does NOT make the row polynomials real-rooted ([[1,0,0],[1,1,0],[1,1,1]] is TNN with row 1+t+t^2)."
    )]
    pub fn check_tnn_network(
        &self,
        Parameters(input): Parameters<TnnNetworkRequest>,
    ) -> Result<Json<TnnNetworkResponse>, McpError> {
        Ok(Json(tnn_network_response(&input)?))
    }

    #[tool(
        description = "Search exactly for a planar strip digraph (edges from column x to x+1 depending only on the height, monotone targets) whose path counts from (0,0) to (n,k) are the coefficients [t^k] P_n of the input rows (P_1, P_2, ... with P_0 = 1 implicit, or first_row_index=0). Bounds: max_target_span, weights (rational|integer), max_search_nodes, max_integer_search_nodes. Verdict: certified (digraph re-verified against every row), impossible_within_bounds (exact), or inconclusive (budget exhausted). Guardrail: a certificate proves that the coefficient matrix of P_0..P_N is TNN (Lindstrom--Gessel--Viennot), NOT that the rows are real-rooted or PF ([1+t, 1+t+t^2] is certified); it covers only the supplied prefix."
    )]
    pub fn recognize_brenti_strip_digraph(
        &self,
        Parameters(input): Parameters<BrentiRecognitionRequest>,
    ) -> Result<Json<BrentiRecognitionResponse>, McpError> {
        Ok(Json(brenti_recognition_response(&input)?))
    }

    #[tool(
        description = "Search exactly for a planar network on {0..N} x Z that is invariant under height shifts (steps (length, rise, weight) into each column n stand for the edges (n-length, y) -> (n, y+rise) at every height y, so P_n = sum weight * t^rise * P_(n-length), P_0 = 1) whose rows are the input polynomials (P_1, P_2, ... or first_row_index=0). Valid networks have primitive steps and no crossing edges; then the Toeplitz matrix of every row is a path matrix of fully compatible sources and sinks, hence TNN (Lindstrom--Gessel--Viennot, Brenti 1995), so every supplied row is a Polya frequency sequence: real-rooted with nonpositive zeros. Bounds: max_step_length, max_rise, weights (rational|integer), max_search_nodes, max_integer_search_nodes. Verdict: certified (re-verified: planarity, rows, rebuilt Toeplitz windows), impossible_within_bounds (exact; a non-real-rooted row rules out every network), or inconclusive. Guardrail: certifies PF rows, NOT TNN of the coefficient matrix (that is check_tnn_network / recognize_brenti_strip_digraph); covers only the supplied rows."
    )]
    pub fn recognize_pf_network(
        &self,
        Parameters(input): Parameters<crate::requests::PfNetworkRequest>,
    ) -> Result<Json<crate::responses::PfNetworkResponse>, McpError> {
        Ok(Json(pf_network_response(&input)?))
    }

    #[tool(
        description = "Chen--Fu--Ruan (arXiv:2512.08369) certificate for a triangle: rows P_0 (or P_1), P_1, ... with deg P_n = n form a lower-triangular A with nonzero diagonal; compute exactly the leading block Q_N of its left production matrix Q = A (1 + A^{-1}) (a_{n,0} = q_{n,0}, a_{n,k} = sum_{i=k..n} q_{n,i} a_{i-1,k-1}) and decide whether Q_N is TNN. If so (certificate re-verified: regeneration, TNN, unified-network Toeplitz windows), rows P_0..P_N are Polya frequency sequences (real-rooted, nonpositive zeros) and A_N and its reversal are TNN; an explicit planar network Q_N = L D is returned. Otherwise returns the failure and a negative minor when found. Guardrail: sufficient, not necessary (Eulerian rows are real-rooted, yet q_{4,2} = -5); covers only the supplied rows."
    )]
    pub fn check_left_production_matrix(
        &self,
        Parameters(input): Parameters<crate::requests::LeftProductionRequest>,
    ) -> Result<Json<crate::responses::LeftProductionResponse>, McpError> {
        Ok(Json(left_production_response(&input)?))
    }

    #[tool(
        description = "Bounded queue search over equivalent representations of one polynomial sequence. request: {source_recurrence (polytool.recurrence.v1 with seeds and first_index) | oeis_id | rows (exploratory fits only), row_origin, row_count (default 12), transformations: [registry maps as in transform_recurrence / transform_recurrence_implicit: multiply_power, divide_power, reversal, affine, mobius, f_polynomial, gamma, symmetric_decomposition, diagonal, derivative_antecedent, and the lossy derivative {order} whose lost coefficients are recorded as lost_data], search (fit bounds, default max_rec_len 2, max_var_deg 2, max_idx_deg 1, max_diff_deg 1, try_denominator true, max_candidates 2000), fits (default true), certified_spaces (default true: complete spaces after factorial diagonal maps), space_shapes, root_domain (optional, attached to coordinate-change steps), limits {max_depth 2, max_representations 16, max_steps 160, max_sources_per_edge 2, max_space_shapes 16, ...}}. Each representation applies every allowed map, deduplicates maps by certified equality of normal forms (never by matching rows), reconstructs the source exactly, transports its certified relations (replayable polytool.recurrence-transform.v1 / polytool.recurrence-implicit.v1 steps), certifies complete spaces, fits the new coordinates, and feeds certified relations to later steps (implicit-only relations through diagonal, derivative_antecedent and derivative maps as implicit_relation steps; for the lossy map the parent's description is assessed with and without the lost data). Returns polytool.recurrence-representation-search.v1: representations with ranked descriptions labelled finite_fit, certified_annihilating_relation or equivalent_generating_description (with seeds, conditions and a uniqueness assessment), edges (new_representation, duplicate_map, rejected, skipped_budget), a global ranking (shift length, derivative order, coefficient degree, index degree, components, reconstruction cost), budget events and status complete | budget_exhausted."
    )]
    pub fn search_recurrence_representations(
        &self,
        Parameters(input): Parameters<RecurrenceRepresentationSearchToolRequest>,
    ) -> Result<Json<BTreeMap<String, serde_json::Value>>, McpError> {
        validate_mcp_text(&input.request.to_string(), "representation-search request")
            .map_err(invalid_params)?;
        let request: polytool::recurrence::RepresentationSearchRequest =
            serde_json::from_value(input.request).map_err(|e| invalid_params(e.to_string()))?;
        request
            .limits
            .check_ceiling(&polytool::recurrence::RepresentationSearchLimits::service_ceiling())
            .map_err(invalid_params)?;
        let report =
            polytool::recurrence::search_representations(&request).map_err(invalid_params)?;
        let serde_json::Value::Object(fields) =
            serde_json::to_value(report).expect("serialize representation search")
        else {
            unreachable!()
        };
        Ok(Json(fields.into_iter().collect()))
    }

    #[tool(
        description = "Replay a polytool.recurrence-representation-search.v1 report from search_recurrence_representations: rerun the search from its embedded request and require equality, replay every embedded transform/implicit step and complete-space certificate with its own verifier, check that each step's defining rule is a certified description of its parent, re-derive every representation from the original source rows (map and exact inverse), and re-check every description on its rows."
    )]
    pub fn verify_recurrence_representation_search(
        &self,
        Parameters(input): Parameters<VerifyRecurrenceRepresentationSearchToolRequest>,
    ) -> Result<Json<BTreeMap<String, serde_json::Value>>, McpError> {
        validate_mcp_text(&input.report.to_string(), "representation-search report")
            .map_err(invalid_params)?;
        let report: polytool::recurrence::RepresentationSearchReport =
            serde_json::from_value(input.report).map_err(|e| invalid_params(e.to_string()))?;
        report
            .request
            .limits
            .check_ceiling(&polytool::recurrence::RepresentationSearchLimits::service_ceiling())
            .map_err(invalid_params)?;
        polytool::recurrence::verify_representation_search(&report).map_err(invalid_params)?;
        Ok(Json(BTreeMap::from([
            ("verified".into(), json!(true)),
            ("status".into(), json!(report.status)),
            (
                "representations".into(),
                json!(report.representations.len()),
            ),
        ])))
    }

    #[tool(
        description = "List EVERY real-rooted integer polynomial p(x) = a_0 x^n + a_1 x^(n-1) + ... + a_n of degree n whose first t leading coefficients are prescribed (3 <= t <= n + 1, a_0 != 0 of either sign); reducible polynomials and repeated roots are included. Exact Greaves--Syatriadi derivative-tower search (arXiv:2504.09241) with certified endpoints; no floating point. Arguments: degree (2..=32), leading_coefficients_descending (a_0, a_1, ..., a_(t-1) in DESCENDING order; integers or exact decimal strings of at most 256 bits; alias leading_desc), and optional budgets max_results (default 500, maximum 2000), max_nodes (default 5000, maximum 25000), max_coefficient_bits (default 1024, maximum 2048); a larger value is rejected. Returns the fields of `polytool enumerate-real-rooted --json`: polynomials as ASCENDING exact integer strings c_0..c_n (c_k = a_(n-k)) sorted lexicographically by (a_t, ..., a_n), count, status complete (the full set, possibly empty) or resource_limit (an initial segment of the full sorted set; limit_reached names max_nodes, max_results or max_coefficient_bits), the applied limits, and deterministic search stats. A complete enumeration is a certified finite computation for that one prefix, not a statement about a family; with only two prescribed coefficients the set can be infinite, which is why t >= 3 is required."
    )]
    pub fn enumerate_real_rooted(
        &self,
        Parameters(input): Parameters<crate::requests::EnumerateRealRootedRequest>,
    ) -> Result<Json<crate::responses::EnumerateRealRootedResponse>, McpError> {
        Ok(Json(crate::enumeration::enumerate_real_rooted_response(
            &input,
        )?))
    }

    #[tool(
        description = "Search implicit differential relations sum c(n,x) D^j Q[n+i] = 0 (derivatives on any row, the newest included) directly from exact rows. request: {source_recurrence (polytool.recurrence.v1 with seeds) | oeis_id | rows + row_origin (finite fits only), row_count (default 12), transformation (optional scalar registry map Q = T(P), e.g. diagonal factorial_power, derivative_antecedent, or the lossy {kind: derivative, order: m}), domain_start, supports: [{role: successor | implicit_successor | row_equation, max_shift, derivative_order, index_degree, variable_degree}] (default implicit_successor (1,1,1,1) and row_equation (0,2,1,1)), conditions (default true), condition_index_degree (default 1), limits (complete-space limits)}. Each support is one exact kernel with a rank certificate. With a constant-index monic pure-lag source and the identity, factorial_power or derivative map a finite theorem makes the kernel the complete space of all-index relations (certified_complete_space); otherwise representatives are finite_fit, or certified_by_division when M o T = sum H_nu S^nu L for the defining rule L. Empty leader projections are absent_in_support (rows_only_obstruction). Exceptional indices of a leader trigger a coefficient-condition search (complete under the same source by a truncation theorem). Returns polytool.recurrence-implicit-search.v1 with representatives (role, evidence, generating description and assessment: equivalent_description, certified_consequence, finite_evidence)."
    )]
    pub fn search_implicit_relations(
        &self,
        Parameters(input): Parameters<ImplicitRelationSearchToolRequest>,
    ) -> Result<Json<BTreeMap<String, serde_json::Value>>, McpError> {
        validate_mcp_text(&input.request.to_string(), "implicit-search request")
            .map_err(invalid_params)?;
        let request: polytool::recurrence::ImplicitRelationSearchRequest =
            serde_json::from_value(input.request).map_err(|e| invalid_params(e.to_string()))?;
        request
            .limits
            .check_ceiling(&Default::default())
            .map_err(invalid_params)?;
        let report =
            polytool::recurrence::search_implicit_relations(&request).map_err(invalid_params)?;
        let serde_json::Value::Object(fields) =
            serde_json::to_value(report).expect("serialize implicit search")
        else {
            unreachable!()
        };
        Ok(Json(fields.into_iter().collect()))
    }

    #[tool(
        description = "Replay a polytool.recurrence-implicit-search.v1 report from search_implicit_relations: rerun the search from its embedded request and require equality, re-check every kernel and condition certificate against matrices rebuilt from the reported rows, every relation and condition on the rows, every division certificate and every assessment."
    )]
    pub fn verify_implicit_relation_search(
        &self,
        Parameters(input): Parameters<VerifyImplicitRelationSearchToolRequest>,
    ) -> Result<Json<BTreeMap<String, serde_json::Value>>, McpError> {
        validate_mcp_text(&input.report.to_string(), "implicit-search report")
            .map_err(invalid_params)?;
        let report: polytool::recurrence::ImplicitRelationSearchReport =
            serde_json::from_value(input.report).map_err(|e| invalid_params(e.to_string()))?;
        report
            .request
            .limits
            .check_ceiling(&Default::default())
            .map_err(invalid_params)?;
        polytool::recurrence::verify_implicit_relation_search(&report).map_err(invalid_params)?;
        Ok(Json(BTreeMap::from([
            ("verified".into(), json!(true)),
            (
                "results".into(),
                json!(report
                    .results
                    .iter()
                    .map(|r| json!({"status": r.status, "evidence": r.evidence}))
                    .collect::<Vec<_>>()),
            ),
        ])))
    }
}

#[tool_handler(router = self.tool_router)]
impl ServerHandler for PolynomialToolsServer {
    fn get_info(&self) -> ServerInfo {
        ServerInfo {
            protocol_version: Default::default(),
            capabilities: ServerCapabilities::builder().enable_tools().build(),
            server_info: Implementation {
                name: "polytool".to_string(),
                title: Some("polytool".to_string()),
                version: env!("CARGO_PKG_VERSION").to_string(),
                description: Some(
                    "Exact univariate polytool for combinatorial research.".to_string(),
                ),
                icons: None,
                website_url: None,
            },
            instructions: Some("Exact univariate polytool for combinatorial research.".to_string()),
        }
    }
}

impl Default for PolynomialToolsServer {
    fn default() -> Self {
        Self::new()
    }
}
