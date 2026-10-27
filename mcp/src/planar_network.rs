//! The planar-network tools: `check_tnn_network` (TNN verdict, canonical
//! network and chain polynomials of a coefficient triangle),
//! `recognize_brenti_strip_digraph` (exact strip-digraph search; TNN of the
//! coefficient matrix), `recognize_pf_network` (height-invariant networks;
//! PF rows) and `check_left_production_matrix` (Chen--Fu--Ruan; TNN triangle,
//! reversal and PF rows).
//!
//! All are thin adapters over `polytool::tnn_network`,
//! `polytool::brenti_recognition`, `polytool::pf_network` and
//! `polytool::left_production`; they only parse, bound and serialize.

use crate::limits::{
    MCP_MAX_BRENTI_DEGREE, MCP_MAX_BRENTI_INTEGER_NODES, MCP_MAX_BRENTI_ROWS,
    MCP_MAX_BRENTI_SEARCH_NODES, MCP_MAX_LEFT_PRODUCTION_ROWS, MCP_MAX_PF_NETWORK_DEGREE,
    MCP_MAX_PF_NETWORK_ROWS, MCP_MAX_PF_NETWORK_STEP, MCP_MAX_TNN_NETWORK_ROWS,
};
use crate::parsing::{bigint_strings, invalid_params, parse_bigint_batch, parse_bigint_matrix};
use crate::requests::{
    BigIntCoefficientInput, BigIntPolynomialBatchInput, BigIntPolynomialInput,
    BrentiRecognitionRequest, BrentiWeightsInput, LeftProductionRequest, PfNetworkRequest,
    TotallyNonnegativeNetworkRequest,
};
use crate::responses::{
    BrentiRecognitionOptionsItem, BrentiRecognitionResponse, ChainInterlacingItem,
    ChainPolynomialItem, ChainPolynomialsItem, NegativeMinorItem, StripDigraphEdgeItem,
    StripDigraphItem, TotallyNonnegativeFailureItem, TotallyNonnegativeNetworkResponse,
};
use crate::responses::{
    LeftProductionResponse, PfNetworkItem, PfNetworkOptionsItem, PfNetworkResponse,
    PfNetworkStepItem, ProductionNetworkItem, RationalNegativeMinorItem,
};
use num_bigint::BigInt;
use polytool::*;
use rmcp::ErrorData as McpError;

/// Guardrail returned with every `check_tnn_network` response.
pub(crate) const TNN_NETWORK_NOTE: &str = "TNN of the coefficient matrix does not imply \
     that the row polynomials are real-rooted: [[1,0,0],[1,1,0],[1,1,1]] is TNN and has \
     the row 1 + t + t^2. Brändén--Saud Maia Leite (arXiv:2412.06595, Theorem 3.7) \
     concerns the chain polynomials, not the rows. A finite triangle certifies only \
     this finite prefix.";

/// Guardrail returned with every `recognize_brenti_strip_digraph` response.
pub(crate) const BRENTI_RECOGNIZE_NOTE: &str = "A certificate shows that the supplied rows \
     are rows of the path matrix of a planar strip digraph with height-dependent \
     transitions, so by Lindström--Gessel--Viennot the coefficient matrix of the rows \
     P_0, ..., P_N is totally nonnegative. It does not show that the rows are real-rooted \
     or Pólya frequency sequences ([1 + t, 1 + t + t^2] is recognized), and it covers \
     only the supplied prefix.";

const CHAIN_THEOREM: &str = "Brändén--Saud Maia Leite, arXiv:2412.06595, Theorem 3.7";

/// Parse a `polynomials`/`text` batch, rejecting the request on the first
/// unparsable row (these tools need every row).
fn parse_rows(
    polynomials: &Option<Vec<BigIntPolynomialInput>>,
    text: &Option<String>,
) -> Result<Vec<Vec<BigInt>>, McpError> {
    let batch = parse_bigint_batch(&BigIntPolynomialBatchInput {
        polynomials: polynomials.clone(),
        text: text.clone(),
    })?;
    batch
        .into_iter()
        .enumerate()
        .map(|(index, item)| {
            item.map(|p| p.coefficients)
                .map_err(|error| invalid_params(format!("polynomial {index}: {error}")))
        })
        .collect()
}

fn parse_tnn_input(
    input: &TotallyNonnegativeNetworkRequest,
) -> Result<(bool, Vec<Vec<BigInt>>), McpError> {
    let given = [
        input.polynomials.is_some(),
        input.text.is_some(),
        input.matrix.is_some(),
    ]
    .iter()
    .filter(|&&present| present)
    .count();
    if given != 1 {
        return Err(invalid_params(
            "expected exactly one of `polynomials`, `text` or `matrix`",
        ));
    }
    let (matrix_mode, rows) = match &input.matrix {
        Some(matrix) => (true, parse_matrix(matrix)?),
        None => (false, parse_rows(&input.polynomials, &input.text)?),
    };
    if rows.is_empty() {
        return Err(invalid_params("at least one row is required"));
    }
    if rows.len() > MCP_MAX_TNN_NETWORK_ROWS {
        return Err(invalid_params(format!(
            "input has {} rows; the MCP limit for check_tnn_network is {MCP_MAX_TNN_NETWORK_ROWS}",
            rows.len()
        )));
    }
    Ok((matrix_mode, rows))
}

fn parse_matrix(matrix: &[Vec<BigIntCoefficientInput>]) -> Result<Vec<Vec<BigInt>>, McpError> {
    parse_bigint_matrix(matrix, "matrix")
}

/// Stable snake-case name of a [`NetworkError`] variant.
fn network_error_kind(error: &NetworkError) -> &'static str {
    match error {
        NetworkError::EmptyMatrix => "empty_matrix",
        NetworkError::EmptyPolynomialSequence => "empty_polynomial_sequence",
        NetworkError::NotSquare { .. } => "not_square",
        NetworkError::NonzeroAboveDiagonal { .. } => "nonzero_above_diagonal",
        NetworkError::UnsupportedStartingDegree { .. } => "unsupported_starting_degree",
        NetworkError::UnexpectedPolynomialDegree { .. } => "unexpected_polynomial_degree",
        NetworkError::NonMonicPolynomial { .. } => "non_monic_polynomial",
        NetworkError::InvalidRowLength { .. } => "invalid_row_length",
        NetworkError::NegativeEntry { .. } => "negative_entry",
        NetworkError::NonUnitDiagonal { .. } => "non_unit_diagonal",
        NetworkError::NotTotallyNonnegative { .. } => "not_tnn",
        NetworkError::ZeroTailViolation { .. } => "zero_tail_violation",
        NetworkError::NegativeReducedEntry { .. } => "negative_reduced_entry",
        NetworkError::ReconstructionMismatch => "reconstruction_mismatch",
    }
}

pub(crate) fn tnn_failure_item(error: &NetworkError) -> TotallyNonnegativeFailureItem {
    let (row, col, value) = match error {
        NetworkError::NotSquare { row, .. } | NetworkError::InvalidRowLength { row, .. } => {
            (Some(*row), None, None)
        }
        NetworkError::UnexpectedPolynomialDegree { index, .. } => (Some(*index), None, None),
        NetworkError::NonMonicPolynomial { index, leading } => (
            Some(*index),
            None,
            Some(
                leading
                    .as_ref()
                    .map_or("0".to_string(), ToString::to_string),
            ),
        ),
        NetworkError::NonzeroAboveDiagonal { row, col, value }
        | NetworkError::NegativeEntry { row, col, value } => {
            (Some(*row), Some(*col), Some(value.to_string()))
        }
        NetworkError::NonUnitDiagonal { row, value } => {
            (Some(*row), Some(*row), Some(value.to_string()))
        }
        NetworkError::ZeroTailViolation { row, value } => {
            (Some(*row), Some(0), Some(value.to_string()))
        }
        NetworkError::NegativeReducedEntry { row, col, value } => {
            (Some(*row), Some(*col), Some(value.to_string()))
        }
        _ => (None, None, None),
    };
    let witness = match error {
        NetworkError::NotTotallyNonnegative {
            witness: Some(minor),
            ..
        } => Some(NegativeMinorItem {
            rows: minor.rows.clone(),
            cols: minor.cols.clone(),
            value: minor.value.to_string(),
        }),
        _ => None,
    };
    TotallyNonnegativeFailureItem {
        kind: network_error_kind(error).to_string(),
        message: error.to_string(),
        row,
        col,
        value,
        witness,
    }
}

fn chain_item(chain: Result<ChainPolynomialReport, NetworkError>) -> ChainPolynomialsItem {
    let report = match chain {
        Ok(report) => report,
        Err(error) => {
            return ChainPolynomialsItem {
                theorem: CHAIN_THEOREM.to_string(),
                error: Some(error.to_string()),
                polynomials: Vec::new(),
                interlacing: Vec::new(),
                conclusion_holds: None,
                contradicts_theorem: None,
            }
        }
    };
    let polynomials = report
        .chain_polynomials
        .iter()
        .enumerate()
        .map(|(n, p)| ChainPolynomialItem {
            n,
            polynomial: format_poly_bigint_coeffs(p),
            coefficients: bigint_strings(p),
            real_rooted: report.real_rooted[n],
            roots_in_minus_one_zero: report.roots_in_minus_one_zero[n],
        })
        .collect();
    let interlacing = report
        .interlacing
        .iter()
        .enumerate()
        .map(|(n, &interlaces)| ChainInterlacingItem {
            left: n,
            right: n + 1,
            interlaces,
        })
        .collect();
    ChainPolynomialsItem {
        theorem: CHAIN_THEOREM.to_string(),
        error: None,
        polynomials,
        interlacing,
        conclusion_holds: Some(report.chain_conclusion_holds()),
        contradicts_theorem: Some(report.contradicts_theorem()),
    }
}

fn rational_rows(rows: &[Vec<BigRational>]) -> Vec<Vec<String>> {
    rows.iter()
        .map(|row| row.iter().map(ToString::to_string).collect())
        .collect()
}

/// Run `check_tnn_network`.
pub(crate) fn tnn_network_response(
    input: &TotallyNonnegativeNetworkRequest,
) -> Result<TotallyNonnegativeNetworkResponse, McpError> {
    let (matrix_mode, lines) = parse_tnn_input(input)?;
    let input_rows = lines.len();
    let rows = if matrix_mode {
        if lines.iter().all(|row| row.len() == input_rows) {
            lower_triangular_rows_from_square(&lines)
        } else {
            Ok(lines)
        }
    } else {
        let polys: Vec<Polynomial<BigInt>> = lines.into_iter().map(Polynomial::new).collect();
        coefficient_matrix_from_monic_polynomials(&polys)
    };
    let prepended_p0 = !matrix_mode && rows.as_ref().is_ok_and(|r| r.len() == input_rows + 1);
    let verdict = rows
        .as_ref()
        .map_err(Clone::clone)
        .and_then(|r| check_lower_unitriangular_totally_nonnegative(r));
    let (lambda, certificate_verified) = match (&rows, &verdict) {
        (Ok(r), Ok(())) => match reconstruct_canonical_totally_nonnegative_network(r) {
            Ok(network) => (
                Some(rational_rows(&network.lambda)),
                Some(verify_path_matrix_certificate(r, &network)),
            ),
            Err(_) => (None, Some(false)),
        },
        _ => (None, None),
    };
    let chain = input.chain.unwrap_or(false).then(|| {
        chain_item(
            rows.as_ref()
                .map_err(Clone::clone)
                .and_then(|r| check_chain_polynomials_bigint(r)),
        )
    });
    Ok(TotallyNonnegativeNetworkResponse {
        input_mode: if matrix_mode { "matrix" } else { "polynomials" }.to_string(),
        input_rows,
        prepended_p0,
        coefficient_matrix: rows
            .as_ref()
            .ok()
            .map(|r| r.iter().map(|row| bigint_strings(row)).collect()),
        totally_nonnegative: verdict.is_ok(),
        failure: verdict.as_ref().err().map(tnn_failure_item),
        lambda,
        certificate_verified,
        chain,
        note: TNN_NETWORK_NOTE.to_string(),
    })
}

fn brenti_options(input: &BrentiRecognitionRequest) -> Result<BrentiRecognitionOptions, McpError> {
    let defaults = BrentiRecognitionOptions::default();
    let max_search_nodes = input.max_search_nodes.unwrap_or(defaults.max_search_nodes);
    if max_search_nodes > MCP_MAX_BRENTI_SEARCH_NODES {
        return Err(invalid_params(format!(
            "`max_search_nodes` {max_search_nodes} exceeds the MCP limit {MCP_MAX_BRENTI_SEARCH_NODES}"
        )));
    }
    let max_integer_search_nodes = input
        .max_integer_search_nodes
        .unwrap_or(defaults.max_integer_search_nodes);
    if max_integer_search_nodes > MCP_MAX_BRENTI_INTEGER_NODES {
        return Err(invalid_params(format!(
            "`max_integer_search_nodes` {max_integer_search_nodes} exceeds the MCP limit \
             {MCP_MAX_BRENTI_INTEGER_NODES}"
        )));
    }
    Ok(BrentiRecognitionOptions {
        first_row_index: input.first_row_index.unwrap_or(defaults.first_row_index),
        max_target_span: input.max_target_span,
        weights: match input.weights.unwrap_or_default() {
            BrentiWeightsInput::Rational => BrentiWeightDomain::Rational,
            BrentiWeightsInput::Integer => BrentiWeightDomain::Integer,
        },
        max_search_nodes,
        max_integer_search_nodes,
    })
}

fn digraph_item(digraph: &BrentiStripDigraph) -> StripDigraphItem {
    StripDigraphItem {
        num_heights: digraph.num_heights(),
        outgoing: digraph
            .outgoing
            .iter()
            .map(|edges| {
                edges
                    .iter()
                    .map(|edge| StripDigraphEdgeItem {
                        to: edge.to,
                        weight: edge.weight.to_string(),
                    })
                    .collect()
            })
            .collect(),
    }
}

/// Run `recognize_brenti_strip_digraph`.
pub(crate) fn brenti_recognition_response(
    input: &BrentiRecognitionRequest,
) -> Result<BrentiRecognitionResponse, McpError> {
    let options = brenti_options(input)?;
    let rows = parse_rows(&input.polynomials, &input.text)?;
    if rows.len() > MCP_MAX_BRENTI_ROWS {
        return Err(invalid_params(format!(
            "input has {} rows; the MCP limit for recognize_brenti_strip_digraph is \
             {MCP_MAX_BRENTI_ROWS}",
            rows.len()
        )));
    }
    let polys: Vec<Polynomial<BigInt>> = rows.into_iter().map(Polynomial::new).collect();
    if let Some((index, degree)) = polys
        .iter()
        .enumerate()
        .filter_map(|(index, p)| p.degree().map(|d| (index, d)))
        .find(|&(_, degree)| degree > MCP_MAX_BRENTI_DEGREE)
    {
        return Err(invalid_params(format!(
            "polynomial {index} has degree {degree}; the MCP limit for \
             recognize_brenti_strip_digraph is {MCP_MAX_BRENTI_DEGREE}"
        )));
    }
    let result = recognize_brenti_strip_digraph(&polys, &options)
        .map_err(|error| invalid_params(error.to_string()))?;
    let options_item = BrentiRecognitionOptionsItem {
        first_row_index: options.first_row_index,
        max_target_span: options.max_target_span,
        weights: match options.weights {
            BrentiWeightDomain::Rational => "rational",
            BrentiWeightDomain::Integer => "integer",
        }
        .to_string(),
        max_search_nodes: options.max_search_nodes,
        max_integer_search_nodes: options.max_integer_search_nodes,
    };
    let mut response = BrentiRecognitionResponse {
        options: options_item,
        rows: polys.len(),
        verdict: String::new(),
        certified: result.is_certified(),
        digraph: None,
        row_indices: None,
        certificate_verified: None,
        supports_checked: None,
        obstruction: None,
        reason: None,
        message: result.to_string(),
        note: BRENTI_RECOGNIZE_NOTE.to_string(),
    };
    match &result {
        BrentiRecognition::Certified {
            digraph,
            certificate,
        } => {
            response.verdict = "certified".to_string();
            response.digraph = Some(digraph_item(digraph));
            response.row_indices = Some(certificate.row_indices.clone());
            response.certificate_verified = Some(
                certificate
                    .verify_against_bigint_polynomials(&polys)
                    .is_ok(),
            );
        }
        BrentiRecognition::ImpossibleWithinBounds {
            supports_checked,
            obstruction,
        } => {
            response.verdict = "impossible_within_bounds".to_string();
            response.supports_checked = Some(*supports_checked);
            response.obstruction = Some(obstruction.clone());
        }
        BrentiRecognition::Inconclusive { reason } => {
            response.verdict = "inconclusive".to_string();
            response.reason = Some(reason.clone());
        }
    }
    Ok(response)
}

/// Guardrail returned with every `recognize_pf_network` response.
pub(crate) const PF_NETWORK_NOTE: &str = "A certificate shows that each supplied row is a \
     row of a planar network that is invariant under height shifts, so the Toeplitz matrix \
     of each row is the path matrix of fully compatible sources and sinks and is totally \
     nonnegative (Lindström--Gessel--Viennot; Brenti 1995): every supplied row is a Pólya \
     frequency sequence, i.e. real-rooted with only nonpositive zeros \
     (Aissen--Schoenberg--Whitney). It covers only the supplied rows, and it does not \
     certify that the coefficient matrix of the rows is TNN (check_tnn_network and \
     recognize_brenti_strip_digraph concern that matrix, not the rows).";

/// Guardrail returned with every `check_left_production_matrix` response.
pub(crate) const LEFT_PRODUCTION_NOTE: &str = "If the left production matrix Q_N is TNN, \
     Chen--Fu--Ruan (arXiv:2512.08369, Theorems 3 and 4) give one planar network whose \
     path matrices are the triangle, its reversal and the Toeplitz matrix of each row, so \
     the rows P_0, ..., P_N are Pólya frequency sequences (real-rooted, nonpositive zeros) \
     and the triangle and its reversal are TNN, for this prefix only. If Q_N is not TNN, \
     nothing follows about the rows (the Eulerian rows are real-rooted, yet q_{4,2} = -5).";

fn bounded_option(
    value: Option<usize>,
    default: usize,
    limit: usize,
    name: &str,
) -> Result<usize, McpError> {
    let value = value.unwrap_or(default);
    if value > limit {
        return Err(invalid_params(format!(
            "`{name}` {value} exceeds the MCP limit {limit}"
        )));
    }
    Ok(value)
}

fn pf_network_options(input: &PfNetworkRequest) -> Result<PfNetworkOptions, McpError> {
    let defaults = PfNetworkOptions::default();
    Ok(PfNetworkOptions {
        first_row_index: input.first_row_index.unwrap_or(defaults.first_row_index),
        max_step_length: bounded_option(
            input.max_step_length,
            defaults.max_step_length,
            MCP_MAX_PF_NETWORK_STEP,
            "max_step_length",
        )?,
        max_rise: bounded_option(
            input.max_rise,
            defaults.max_rise,
            MCP_MAX_PF_NETWORK_STEP,
            "max_rise",
        )?,
        weights: match input.weights.unwrap_or_default() {
            BrentiWeightsInput::Rational => BrentiWeightDomain::Rational,
            BrentiWeightsInput::Integer => BrentiWeightDomain::Integer,
        },
        max_search_nodes: bounded_option(
            input.max_search_nodes,
            defaults.max_search_nodes,
            MCP_MAX_BRENTI_SEARCH_NODES,
            "max_search_nodes",
        )?,
        max_integer_search_nodes: bounded_option(
            input.max_integer_search_nodes,
            defaults.max_integer_search_nodes,
            MCP_MAX_BRENTI_INTEGER_NODES,
            "max_integer_search_nodes",
        )?,
    })
}

fn pf_network_item(network: &HeightInvariantNetwork) -> PfNetworkItem {
    PfNetworkItem {
        last_column: network.last_column(),
        steps: network
            .steps
            .iter()
            .map(|column| {
                column
                    .iter()
                    .map(|step| PfNetworkStepItem {
                        length: step.length,
                        rise: step.rise,
                        weight: step.weight.to_string(),
                    })
                    .collect()
            })
            .collect(),
    }
}

/// Run `recognize_pf_network`.
pub(crate) fn pf_network_response(input: &PfNetworkRequest) -> Result<PfNetworkResponse, McpError> {
    let options = pf_network_options(input)?;
    let rows = parse_rows(&input.polynomials, &input.text)?;
    if rows.len() > MCP_MAX_PF_NETWORK_ROWS {
        return Err(invalid_params(format!(
            "input has {} rows; the MCP limit for recognize_pf_network is \
             {MCP_MAX_PF_NETWORK_ROWS}",
            rows.len()
        )));
    }
    let polys: Vec<Polynomial<BigInt>> = rows.into_iter().map(Polynomial::new).collect();
    if let Some((index, degree)) = polys
        .iter()
        .enumerate()
        .filter_map(|(index, p)| p.degree().map(|d| (index, d)))
        .find(|&(_, degree)| degree > MCP_MAX_PF_NETWORK_DEGREE)
    {
        return Err(invalid_params(format!(
            "polynomial {index} has degree {degree}; the MCP limit for recognize_pf_network \
             is {MCP_MAX_PF_NETWORK_DEGREE}"
        )));
    }
    let result = recognize_pf_network(&polys, &options)
        .map_err(|error| invalid_params(error.to_string()))?;
    let mut response = PfNetworkResponse {
        options: PfNetworkOptionsItem {
            first_row_index: options.first_row_index,
            max_step_length: options.max_step_length,
            max_rise: options.max_rise,
            weights: match options.weights {
                BrentiWeightDomain::Rational => "rational",
                BrentiWeightDomain::Integer => "integer",
            }
            .to_string(),
            max_search_nodes: options.max_search_nodes,
            max_integer_search_nodes: options.max_integer_search_nodes,
        },
        rows: polys.len(),
        verdict: String::new(),
        certified: result.is_certified(),
        network: None,
        row_indices: None,
        certificate_verified: None,
        supports_checked: None,
        obstruction: None,
        reason: None,
        message: result.to_string(),
        note: PF_NETWORK_NOTE.to_string(),
    };
    match &result {
        PfNetworkRecognition::Certified { certificate } => {
            response.verdict = "certified".to_string();
            response.network = Some(pf_network_item(&certificate.network));
            response.row_indices = Some(certificate.row_indices.clone());
            response.certificate_verified = Some(
                certificate
                    .verify_against_bigint_polynomials(&polys)
                    .is_ok(),
            );
        }
        PfNetworkRecognition::ImpossibleWithinBounds {
            supports_checked,
            obstruction,
        } => {
            response.verdict = "impossible_within_bounds".to_string();
            response.supports_checked = Some(*supports_checked);
            response.obstruction = Some(obstruction.clone());
        }
        PfNetworkRecognition::Inconclusive { reason } => {
            response.verdict = "inconclusive".to_string();
            response.reason = Some(reason.clone());
        }
    }
    Ok(response)
}

/// Run `check_left_production_matrix`.
pub(crate) fn left_production_response(
    input: &LeftProductionRequest,
) -> Result<LeftProductionResponse, McpError> {
    let first_row_index = input.first_row_index.unwrap_or(1);
    let rows = parse_rows(&input.polynomials, &input.text)?;
    if rows.len() > MCP_MAX_LEFT_PRODUCTION_ROWS {
        return Err(invalid_params(format!(
            "input has {} rows; the MCP limit for check_left_production_matrix is \
             {MCP_MAX_LEFT_PRODUCTION_ROWS}",
            rows.len()
        )));
    }
    let polys: Vec<Polynomial<BigInt>> = rows.into_iter().map(Polynomial::new).collect();
    let report = check_left_production_matrix(&polys, first_row_index)
        .map_err(|error| invalid_params(error.to_string()))?;
    Ok(LeftProductionResponse {
        first_row_index,
        rows: polys.len(),
        triangle: report
            .triangle
            .iter()
            .map(|row| bigint_strings(row))
            .collect(),
        production_matrix: rational_rows(&report.production_matrix),
        totally_nonnegative: report.totally_nonnegative,
        failure_reason: report.failure_reason.clone(),
        negative_minor: report
            .negative_minor
            .as_ref()
            .map(|minor| RationalNegativeMinorItem {
                rows: minor.rows.clone(),
                cols: minor.cols.clone(),
                value: minor.value.to_string(),
            }),
        certificate_verified: report.certificate_verified,
        rows_certified_pf: report.rows_certified_pf(),
        production_network: report.production_network.as_ref().map(|network| {
            ProductionNetworkItem {
                lambda: rational_rows(&network.network.lambda),
                column_scales: network
                    .column_scales
                    .iter()
                    .map(ToString::to_string)
                    .collect(),
            }
        }),
        message: report.to_string(),
        note: LEFT_PRODUCTION_NOTE.to_string(),
    })
}
