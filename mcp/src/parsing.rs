//! Input parsing, normalization and size validation shared by the tools.
//!
//! Parse failures are reported per polynomial where a batch allows it, and as
//! `invalid_params` errors otherwise.

use crate::limits::{
    MCP_MAX_BATCH_POLYNOMIALS, MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL, MCP_MAX_COEFFICIENT_TEXT_BYTES,
    MCP_MAX_RECURRENCE_INPUT_ROWS, MCP_MAX_TEXT_BYTES, MCP_MAX_TOTAL_COEFFICIENTS,
};
use crate::requests::{
    BigIntCoefficientInput, BigIntPolynomialBatchInput, BigIntPolynomialInput,
    FindRecurrenceRequest, PolynomialBatchInput, PolynomialInput, RationalCoefficientInput,
};
use crate::responses::{BigIntNormalizedPolynomial, NormalizedPolynomial, ParsePolynomialItem};
use num_bigint::BigInt;
use polytool::recurrence::{parse_rational_coeff, BigRational};
use polytool::*;
use rmcp::ErrorData as McpError;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct NormalizedBigIntPolynomial {
    pub(crate) display: BigIntNormalizedPolynomial,
    pub(crate) coefficients: Vec<BigInt>,
}

pub(crate) type ParsedBatch = Vec<Result<NormalizedPolynomial, String>>;
pub(crate) type BigIntParsedBatch = Vec<Result<NormalizedBigIntPolynomial, String>>;
type RationalParsedBatch = Vec<Result<Vec<BigRational>, String>>;
pub(crate) type FamilySourceParse = (
    String,
    Result<Vec<NormalizedPolynomial>, Vec<ParsePolynomialItem>>,
);

pub(crate) fn invalid_params(message: impl std::fmt::Display) -> McpError {
    McpError::invalid_params(message.to_string(), None)
}

pub(crate) fn validate_mcp_text(text: &str, description: &str) -> Result<(), String> {
    if text.len() > MCP_MAX_TEXT_BYTES {
        return Err(format!(
            "{description} is {} bytes; the MCP limit is {MCP_MAX_TEXT_BYTES}",
            text.len()
        ));
    }
    Ok(())
}

fn validate_mcp_polynomial_length<T>(coefficients: Vec<T>) -> Result<Vec<T>, String> {
    if coefficients.len() > MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL {
        return Err(format!(
            "polynomial has {} coefficients; the MCP limit is {MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL}",
            coefficients.len()
        ));
    }
    Ok(coefficients)
}

fn validate_mcp_batch_size(item_count: usize) -> Result<(), McpError> {
    if item_count > MCP_MAX_BATCH_POLYNOMIALS {
        return Err(invalid_params(format!(
            "batch has {item_count} polynomials; the MCP limit is {MCP_MAX_BATCH_POLYNOMIALS}"
        )));
    }
    Ok(())
}

pub(crate) fn parse_polynomial_input(input: &PolynomialInput) -> Result<Vec<i64>, String> {
    match (&input.coefficients, &input.expression) {
        (Some(coefficients), None) => validate_mcp_polynomial_length(coefficients.clone()),
        (None, Some(expression)) => {
            validate_mcp_text(expression, "polynomial expression")?;
            validate_mcp_polynomial_length(
                parse_polynomial(expression).map_err(|error| error.to_string())?,
            )
        }
        (Some(_), Some(_)) => {
            Err("expected exactly one of `coefficients` or `expression`, got both".to_string())
        }
        (None, None) => Err("expected exactly one of `coefficients` or `expression`".to_string()),
    }
}

pub(crate) fn parse_bigint_coefficient(input: &BigIntCoefficientInput) -> Result<BigInt, String> {
    match input {
        BigIntCoefficientInput::Integer(value) => Ok(BigInt::from(*value)),
        BigIntCoefficientInput::Text(value) => {
            if value.len() > MCP_MAX_COEFFICIENT_TEXT_BYTES {
                return Err(format!(
                    "integer coefficient is {} bytes; the MCP limit is {MCP_MAX_COEFFICIENT_TEXT_BYTES}",
                    value.len()
                ));
            }
            value
                .parse::<BigInt>()
                .map_err(|e| format!("invalid integer '{}': {}", value, e))
        }
    }
}

pub(crate) fn parse_bigint_coefficients(
    input: &[BigIntCoefficientInput],
) -> Result<Vec<BigInt>, McpError> {
    if input.len() > MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL {
        return Err(invalid_params(format!(
            "polynomial has {} coefficients; the MCP limit is {MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL}",
            input.len()
        )));
    }
    input
        .iter()
        .map(parse_bigint_coefficient)
        .collect::<Result<Vec<_>, _>>()
        .map_err(invalid_params)
}

fn parse_bigint_polynomial_input(input: &BigIntPolynomialInput) -> Result<Vec<BigInt>, String> {
    match (&input.coefficients, &input.expression) {
        (Some(coefficients), None) => validate_mcp_polynomial_length(
            coefficients
                .iter()
                .enumerate()
                .map(|(index, coefficient)| {
                    parse_bigint_coefficient(coefficient)
                        .map_err(|error| format!("coefficient {index}: {error}"))
                })
                .collect::<Result<Vec<_>, _>>()?,
        ),
        (None, Some(expression)) => {
            validate_mcp_text(expression, "polynomial expression")?;
            validate_mcp_polynomial_length(
                parse_polynomial_bigint(expression).map_err(|error| error.to_string())?,
            )
        }
        (Some(_), Some(_)) => {
            Err("expected exactly one of `coefficients` or `expression`, got both".to_string())
        }
        (None, None) => Err("expected exactly one of `coefficients` or `expression`".to_string()),
    }
}

pub(crate) fn parse_required_bigint_polynomial(
    input: &BigIntPolynomialInput,
    name: &str,
) -> Result<Vec<BigInt>, McpError> {
    parse_bigint_polynomial_input(input).map_err(|e| invalid_params(format!("{name}: {e}")))
}

fn normalize_coefficients(mut coefficients: Vec<i64>) -> Vec<i64> {
    while coefficients.len() > 1 && coefficients.last() == Some(&0) {
        coefficients.pop();
    }
    if coefficients.is_empty() {
        coefficients.push(0);
    }
    coefficients
}

fn normalize_rational_coefficients(mut coefficients: Vec<BigRational>) -> Vec<BigRational> {
    let zero = parse_rational("0").expect("0 parses as BigRational");
    while coefficients.len() > 1 && coefficients.last() == Some(&zero) {
        coefficients.pop();
    }
    if coefficients.is_empty() {
        coefficients.push(zero);
    }
    coefficients
}

fn normalize_bigint_coefficients(mut coefficients: Vec<BigInt>) -> Vec<BigInt> {
    while coefficients.len() > 1 && coefficients.last() == Some(&BigInt::from(0)) {
        coefficients.pop();
    }
    if coefficients.is_empty() {
        coefficients.push(BigInt::from(0));
    }
    coefficients
}

fn degree(coefficients: &[i64]) -> usize {
    coefficients.iter().rposition(|&c| c != 0).unwrap_or(0)
}

fn bigint_degree(coefficients: &[BigInt]) -> usize {
    coefficients
        .iter()
        .rposition(|c| c != &BigInt::from(0))
        .unwrap_or(0)
}

pub(crate) fn bigint_strings(coefficients: &[BigInt]) -> Vec<String> {
    coefficients.iter().map(ToString::to_string).collect()
}

pub(crate) fn normalize_polynomial(coefficients: Vec<i64>) -> NormalizedPolynomial {
    let coefficients = normalize_coefficients(coefficients);
    NormalizedPolynomial {
        polynomial: format_poly(&coefficients),
        degree: degree(&coefficients),
        coefficients,
    }
}

pub(crate) fn normalize_bigint_polynomial(coefficients: Vec<BigInt>) -> NormalizedBigIntPolynomial {
    let coefficients = normalize_bigint_coefficients(coefficients);
    NormalizedBigIntPolynomial {
        display: BigIntNormalizedPolynomial {
            polynomial: format_poly_bigint_coeffs(&coefficients),
            coefficients: bigint_strings(&coefficients),
            degree: bigint_degree(&coefficients),
        },
        coefficients,
    }
}

pub(crate) fn parse_batch(input: &PolynomialBatchInput) -> Result<ParsedBatch, McpError> {
    let parsed: ParsedBatch = match (&input.polynomials, &input.text) {
        (Some(polynomials), None) => {
            validate_mcp_batch_size(polynomials.len())?;
            polynomials
                .iter()
                .map(|p| parse_polynomial_input(p).map(normalize_polynomial))
                .collect()
        }
        (None, Some(text)) => {
            validate_mcp_text(text, "polynomial batch").map_err(invalid_params)?;
            validate_mcp_batch_size(
                text.lines()
                    .filter(|line| {
                        let line = line.trim();
                        !line.is_empty() && !line.starts_with('#')
                    })
                    .count(),
            )?;
            parse_polynomials(text)
                .into_iter()
                .map(|result| {
                    result
                        .map_err(|error| error.to_string())
                        .and_then(validate_mcp_polynomial_length)
                        .map(normalize_polynomial)
                })
                .collect()
        }
        (Some(_), Some(_)) => Err(invalid_params(
            "expected exactly one of `polynomials` or `text`, got both",
        ))?,
        (None, None) => Err(invalid_params(
            "expected exactly one of `polynomials` or `text`",
        ))?,
    };
    let total = parsed
        .iter()
        .filter_map(|item| item.as_ref().ok())
        .try_fold(0usize, |total, polynomial| {
            total.checked_add(polynomial.coefficients.len())
        });
    let total = total.ok_or_else(|| invalid_params("total coefficient count overflow"))?;
    if total > MCP_MAX_TOTAL_COEFFICIENTS {
        return Err(invalid_params(format!(
            "batch has {total} coefficients; the MCP limit is {MCP_MAX_TOTAL_COEFFICIENTS}"
        )));
    }
    Ok(parsed)
}

pub(crate) fn parse_bigint_batch(
    input: &BigIntPolynomialBatchInput,
) -> Result<BigIntParsedBatch, McpError> {
    let parsed: BigIntParsedBatch = match (&input.polynomials, &input.text) {
        (Some(polynomials), None) => {
            validate_mcp_batch_size(polynomials.len())?;
            polynomials
                .iter()
                .map(|p| parse_bigint_polynomial_input(p).map(normalize_bigint_polynomial))
                .collect()
        }
        (None, Some(text)) => {
            validate_mcp_text(text, "polynomial batch").map_err(invalid_params)?;
            validate_mcp_batch_size(
                text.lines()
                    .filter(|line| {
                        let line = line.trim();
                        !line.is_empty() && !line.starts_with('#')
                    })
                    .count(),
            )?;
            parse_polynomials_bigint(text)
                .into_iter()
                .map(|result| {
                    result
                        .map_err(|error| error.to_string())
                        .and_then(validate_mcp_polynomial_length)
                        .map(normalize_bigint_polynomial)
                })
                .collect()
        }
        (Some(_), Some(_)) => Err(invalid_params(
            "expected exactly one of `polynomials` or `text`, got both",
        ))?,
        (None, None) => Err(invalid_params(
            "expected exactly one of `polynomials` or `text`",
        ))?,
    };
    let total = parsed
        .iter()
        .filter_map(|item| item.as_ref().ok())
        .try_fold(0usize, |total, polynomial| {
            total.checked_add(polynomial.coefficients.len())
        });
    let total = total.ok_or_else(|| invalid_params("total coefficient count overflow"))?;
    if total > MCP_MAX_TOTAL_COEFFICIENTS {
        return Err(invalid_params(format!(
            "batch has {total} coefficients; the MCP limit is {MCP_MAX_TOTAL_COEFFICIENTS}"
        )));
    }
    Ok(parsed)
}

fn i64_coefficients_to_rational(coefficients: &[i64]) -> Vec<BigRational> {
    coefficients
        .iter()
        .map(|coeff| parse_rational(&coeff.to_string()).expect("i64 parses as BigRational"))
        .collect()
}

fn parse_rational_coefficient(input: &RationalCoefficientInput) -> Result<BigRational, String> {
    match input {
        RationalCoefficientInput::Integer(value) => parse_rational(&value.to_string()),
        RationalCoefficientInput::Text(value) => {
            if value.len() > MCP_MAX_COEFFICIENT_TEXT_BYTES {
                return Err(format!(
                    "rational coefficient is {} bytes; the MCP limit is {MCP_MAX_COEFFICIENT_TEXT_BYTES}",
                    value.len()
                ));
            }
            parse_rational(value)
        }
    }
}

fn parse_rational_coefficients(
    coefficients: &[RationalCoefficientInput],
) -> Result<Vec<BigRational>, String> {
    if coefficients.len() > MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL {
        return Err(format!(
            "polynomial has {} coefficients; the MCP limit is {MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL}",
            coefficients.len()
        ));
    }
    coefficients
        .iter()
        .enumerate()
        .map(|(index, coefficient)| {
            parse_rational_coefficient(coefficient)
                .map_err(|error| format!("coefficient {index}: {error}"))
        })
        .collect::<Result<Vec<_>, _>>()
        .map(normalize_rational_coefficients)
}

fn validate_recurrence_batch(batch: RationalParsedBatch) -> Result<RationalParsedBatch, McpError> {
    if batch.len() > MCP_MAX_RECURRENCE_INPUT_ROWS {
        return Err(invalid_params(format!(
            "recurrence input has {} rows; the MCP limit is {MCP_MAX_RECURRENCE_INPUT_ROWS}",
            batch.len()
        )));
    }
    let total = batch
        .iter()
        .filter_map(|item| item.as_ref().ok())
        .try_fold(0usize, |total, polynomial| {
            total.checked_add(polynomial.len())
        })
        .ok_or_else(|| invalid_params("total recurrence coefficient count overflow"))?;
    if total > MCP_MAX_TOTAL_COEFFICIENTS {
        return Err(invalid_params(format!(
            "recurrence input has {total} coefficients; the MCP limit is {MCP_MAX_TOTAL_COEFFICIENTS}"
        )));
    }
    Ok(batch)
}

pub(crate) fn parse_recurrence_batch_rational(
    input: &FindRecurrenceRequest,
) -> Result<RationalParsedBatch, McpError> {
    let source_count = [
        input.polynomials.is_some(),
        input.coefficients.is_some(),
        input.expressions.is_some(),
        input.text.is_some(),
    ]
    .into_iter()
    .filter(|present| *present)
    .count();

    if source_count != 1 {
        return Err(invalid_params(
            "expected exactly one of `polynomials`, `coefficients`, `expressions`, or `text`",
        ));
    }

    if let Some(polynomials) = &input.polynomials {
        if polynomials.len() > MCP_MAX_RECURRENCE_INPUT_ROWS {
            return Err(invalid_params(format!(
                "recurrence input has {} rows; the MCP limit is {MCP_MAX_RECURRENCE_INPUT_ROWS}",
                polynomials.len()
            )));
        }
        let batch = parse_batch(&PolynomialBatchInput {
            polynomials: Some(polynomials.clone()),
            text: None,
        })?;
        return validate_recurrence_batch(
            batch
                .into_iter()
                .map(|item| {
                    item.map(|polynomial| i64_coefficients_to_rational(&polynomial.coefficients))
                })
                .collect(),
        );
    }

    if let Some(coefficients) = &input.coefficients {
        if coefficients.len() > MCP_MAX_RECURRENCE_INPUT_ROWS {
            return Err(invalid_params(format!(
                "recurrence input has {} rows; the MCP limit is {MCP_MAX_RECURRENCE_INPUT_ROWS}",
                coefficients.len()
            )));
        }
        let total = coefficients
            .iter()
            .try_fold(0usize, |total, row| total.checked_add(row.len()))
            .ok_or_else(|| invalid_params("total recurrence coefficient count overflow"))?;
        if total > MCP_MAX_TOTAL_COEFFICIENTS {
            return Err(invalid_params(format!(
                "recurrence input has {total} coefficients; the MCP limit is {MCP_MAX_TOTAL_COEFFICIENTS}"
            )));
        }
        return validate_recurrence_batch(
            coefficients
                .iter()
                .map(|coefficients| parse_rational_coefficients(coefficients))
                .collect(),
        );
    }

    if let Some(expressions) = &input.expressions {
        if expressions.len() > MCP_MAX_RECURRENCE_INPUT_ROWS {
            return Err(invalid_params(format!(
                "recurrence input has {} rows; the MCP limit is {MCP_MAX_RECURRENCE_INPUT_ROWS}",
                expressions.len()
            )));
        }
        let polynomials = expressions
            .iter()
            .cloned()
            .map(|expression| PolynomialInput {
                coefficients: None,
                expression: Some(expression),
            })
            .collect();
        let batch = parse_batch(&PolynomialBatchInput {
            polynomials: Some(polynomials),
            text: None,
        })?;
        return validate_recurrence_batch(
            batch
                .into_iter()
                .map(|item| {
                    item.map(|polynomial| i64_coefficients_to_rational(&polynomial.coefficients))
                })
                .collect(),
        );
    }

    if let Some(text) = &input.text {
        validate_mcp_text(text, "recurrence polynomial batch").map_err(invalid_params)?;
        let row_count = text
            .lines()
            .filter(|line| {
                let line = line.trim();
                !line.is_empty() && !line.starts_with('#')
            })
            .count();
        if row_count > MCP_MAX_RECURRENCE_INPUT_ROWS {
            return Err(invalid_params(format!(
                "recurrence input has {row_count} rows; the MCP limit is {MCP_MAX_RECURRENCE_INPUT_ROWS}"
            )));
        }
    }
    let batch = parse_batch(&PolynomialBatchInput {
        polynomials: None,
        text: input.text.clone(),
    })?;
    validate_recurrence_batch(
        batch
            .into_iter()
            .map(|item| {
                item.map(|polynomial| i64_coefficients_to_rational(&polynomial.coefficients))
            })
            .collect(),
    )
}

pub(crate) fn parse_bigint_matrix(
    input: &[Vec<BigIntCoefficientInput>],
    name: &str,
) -> Result<Vec<Vec<BigInt>>, McpError> {
    if input.len() > MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL {
        return Err(invalid_params(format!(
            "`{name}` has {} rows; the MCP limit is {MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL}",
            input.len()
        )));
    }
    let cells: usize = input.iter().map(Vec::len).sum();
    if cells > MCP_MAX_TOTAL_COEFFICIENTS {
        return Err(invalid_params(format!(
            "`{name}` has {cells} entries; the MCP limit is {MCP_MAX_TOTAL_COEFFICIENTS}"
        )));
    }
    input
        .iter()
        .map(|row| parse_bigint_coefficients(row))
        .collect()
}

pub(crate) fn collect_polynomials_or_errors(
    batch: ParsedBatch,
) -> Result<Vec<NormalizedPolynomial>, Vec<ParsePolynomialItem>> {
    let errors: Vec<ParsePolynomialItem> = batch
        .iter()
        .enumerate()
        .filter_map(|(index, item)| match item {
            Ok(_) => None,
            Err(error) => Some(ParsePolynomialItem {
                index,
                ok: false,
                polynomial: None,
                coefficients: None,
                degree: None,
                error: Some(error.clone()),
            }),
        })
        .collect();
    if errors.is_empty() {
        Ok(batch.into_iter().map(Result::unwrap).collect())
    } else {
        Err(errors)
    }
}

pub(crate) fn collect_rational_polynomials_or_errors(
    batch: RationalParsedBatch,
) -> Result<Vec<Vec<BigRational>>, Vec<ParsePolynomialItem>> {
    let errors: Vec<ParsePolynomialItem> = batch
        .iter()
        .enumerate()
        .filter_map(|(index, item)| match item {
            Ok(_) => None,
            Err(error) => Some(ParsePolynomialItem {
                index,
                ok: false,
                polynomial: None,
                coefficients: None,
                degree: None,
                error: Some(error.clone()),
            }),
        })
        .collect();
    if errors.is_empty() {
        Ok(batch.into_iter().map(Result::unwrap).collect())
    } else {
        Err(errors)
    }
}

pub(crate) fn integer_polys_to_rational(coefficients: &[Vec<i64>]) -> Vec<Vec<BigRational>> {
    coefficients
        .iter()
        .map(|row| {
            row.iter()
                .map(|coeff| {
                    parse_rational_coeff(&coeff.to_string())
                        .expect("integer coefficients parse as rationals")
                })
                .collect()
        })
        .collect()
}

pub(crate) fn parse_rational(input: &str) -> Result<BigRational, String> {
    input
        .parse::<BigRational>()
        .map_err(|e| format!("invalid rational `{input}`: {e}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    use serde_json::json;

    #[test]
    fn bigint_json_rejects_floating_point_with_actionable_message() {
        let error = serde_json::from_value::<BigIntPolynomialInput>(json!({
            "coefficients": [1.2345678901234568e20]
        }))
        .unwrap_err()
        .to_string();

        assert!(error.contains("floating-point JSON number"), "{error}");
        assert!(error.contains("quoted decimal strings"), "{error}");
        assert!(error.contains("9007199254740991"), "{error}");
    }

    #[test]
    fn bigint_json_preserves_unsigned_integers_beyond_i64() {
        let input = serde_json::from_value::<BigIntPolynomialInput>(json!({
            "coefficients": [u64::MAX]
        }))
        .unwrap();
        let parsed = parse_bigint_polynomial_input(&input).unwrap();

        assert_eq!(parsed, vec![BigInt::from(u64::MAX)]);
    }
}
