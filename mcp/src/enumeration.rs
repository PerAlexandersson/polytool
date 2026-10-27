//! `enumerate_real_rooted`: validation, bounded execution and response of the
//! real-rooted integer polynomial enumeration.

use crate::limits::{
    MCP_DEFAULT_ENUMERATE_COEFFICIENT_BITS, MCP_DEFAULT_ENUMERATE_NODES,
    MCP_DEFAULT_ENUMERATE_RESULTS, MCP_MAX_ENUMERATE_COEFFICIENT_BITS, MCP_MAX_ENUMERATE_DEGREE,
    MCP_MAX_ENUMERATE_NODES, MCP_MAX_ENUMERATE_PREFIX_BITS, MCP_MAX_ENUMERATE_RESULTS,
};
use crate::parsing::{invalid_params, parse_bigint_coefficient};
use crate::requests::EnumerateRealRootedRequest;
use crate::responses::{
    EnumerateRealRootedLimits, EnumerateRealRootedResponse, EnumerateRealRootedStats,
};
use num_bigint::BigInt;
use polytool::real_rooted_enumeration::{
    enumerate_real_rooted_integer_polynomials_bigint, RealRootedEnumerationError,
    RealRootedEnumerationOptions, RealRootedEnumerationStatus,
};
use rmcp::ErrorData as McpError;

/// A requested budget: the default when absent, an error above the ceiling.
fn bounded_limit<T: PartialOrd + std::fmt::Display + Copy>(
    name: &str,
    requested: Option<T>,
    default: T,
    ceiling: T,
) -> Result<T, McpError> {
    let value = requested.unwrap_or(default);
    if value > ceiling {
        return Err(invalid_params(format!(
            "{name} = {value} exceeds the MCP limit {ceiling}"
        )));
    }
    Ok(value)
}

pub(crate) fn enumerate_real_rooted_response(
    input: &EnumerateRealRootedRequest,
) -> Result<EnumerateRealRootedResponse, McpError> {
    if input.degree > MCP_MAX_ENUMERATE_DEGREE {
        return Err(invalid_params(format!(
            "degree = {} exceeds the MCP limit {MCP_MAX_ENUMERATE_DEGREE}",
            input.degree
        )));
    }
    // The prefix has at most `degree + 1` entries once the library accepts it;
    // reject longer input before parsing any coefficient.
    if input.leading_coefficients_descending.len() > MCP_MAX_ENUMERATE_DEGREE + 1 {
        return Err(invalid_params(format!(
            "at most {} leading coefficients are accepted",
            MCP_MAX_ENUMERATE_DEGREE + 1
        )));
    }
    let mut prefix = Vec::with_capacity(input.leading_coefficients_descending.len());
    for coefficient in &input.leading_coefficients_descending {
        let value: BigInt = parse_bigint_coefficient(coefficient).map_err(invalid_params)?;
        if value.bits() > MCP_MAX_ENUMERATE_PREFIX_BITS {
            return Err(invalid_params(format!(
                "a leading coefficient has {} bits; the MCP limit is {MCP_MAX_ENUMERATE_PREFIX_BITS}",
                value.bits()
            )));
        }
        prefix.push(value);
    }
    let limits = EnumerateRealRootedLimits {
        max_nodes: bounded_limit(
            "max_nodes",
            input.max_nodes,
            MCP_DEFAULT_ENUMERATE_NODES,
            MCP_MAX_ENUMERATE_NODES,
        )?,
        max_results: bounded_limit(
            "max_results",
            input.max_results,
            MCP_DEFAULT_ENUMERATE_RESULTS,
            MCP_MAX_ENUMERATE_RESULTS,
        )?,
        max_coefficient_bits: bounded_limit(
            "max_coefficient_bits",
            input.max_coefficient_bits,
            MCP_DEFAULT_ENUMERATE_COEFFICIENT_BITS,
            MCP_MAX_ENUMERATE_COEFFICIENT_BITS,
        )?,
    };
    let options = RealRootedEnumerationOptions {
        max_nodes: Some(limits.max_nodes),
        max_results: Some(limits.max_results),
        max_coefficient_bits: Some(limits.max_coefficient_bits),
    };
    let result = enumerate_real_rooted_integer_polynomials_bigint(input.degree, &prefix, &options)
        .map_err(|error| match error {
            RealRootedEnumerationError::CertificationFailure { .. } => {
                McpError::internal_error(error.to_string(), None)
            }
            other => invalid_params(other.to_string()),
        })?;
    let stats = result.stats;
    Ok(EnumerateRealRootedResponse {
        command: "enumerate-real-rooted".to_string(),
        degree: input.degree,
        leading_coefficients_descending: prefix.iter().map(ToString::to_string).collect(),
        output_coefficient_order: "ascending".to_string(),
        status: result.status.as_str().to_string(),
        limit_reached: match result.status {
            RealRootedEnumerationStatus::Complete => None,
            RealRootedEnumerationStatus::ResourceLimit(limit) => Some(limit.as_str().to_string()),
        },
        limits,
        count: result.polynomials.len(),
        polynomials: result
            .polynomials
            .iter()
            .map(|polynomial| polynomial.iter().map(ToString::to_string).collect())
            .collect(),
        stats: EnumerateRealRootedStats {
            derivative_levels: stats.derivative_levels,
            nodes_examined: stats.nodes_examined,
            nodes_pruned: stats.nodes_pruned,
            exact_root_isolations: stats.exact_root_isolations,
            refinement_steps: stats.refinement_steps,
            endpoint_tests: stats.endpoint_tests,
            final_verifications: stats.final_verifications,
            results_emitted: stats.results_emitted,
            polynomials_by_level: stats.polynomials_by_level,
        },
    })
}
