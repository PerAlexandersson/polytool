//! The `check_polynomial_family` workflow: family sources, per-polynomial
//! reports, finite Lace TNN checks, recurrence search and the Markdown summary.

use crate::limits::{MCP_MAX_LACE_MATRIX_CELLS, MCP_MAX_LACE_MINORS, MCP_MAX_SEQUENCE_N};
use crate::parsing::{
    bigint_strings, invalid_params, normalize_bigint_polynomial, parse_bigint_batch,
    BigIntParsedBatch, NormalizedBigIntPolynomial,
};
use crate::recurrence_options::{recurrence_failure_response, recurrence_json_string};
use crate::report_items::{bigint_matrix_strings, interlacing_result, polynomial_property_values};
use crate::requests::{
    BigIntPolynomialBatchInput, CheckPolynomialFamilyRequest, FamilyCheckOptions, LaceCheckRequest,
    SequenceKind,
};
use crate::responses::{
    FamilyPolynomialReport, FindRecurrenceResponse, InterlacingPairItem, LaceCheckResponse,
    ParsePolynomialItem, RecurrenceSearchStatus,
};
use num_bigint::BigInt;
use polytool::linalg::check_total_positivity_bigint;
use polytool::recurrence::{
    find_recurrence_adaptive_rational_with_budget, AdaptiveSearchBudget, AdaptiveSearchOptions,
    AdaptiveSearchOutcome, BigRational,
};
use polytool::sequences::{
    chebyshev_polynomials_t_bigint, chebyshev_polynomials_u_bigint, eulerian_polynomials_bigint,
    hermite_polynomials_bigint, narayana_polynomials_bigint, type_b_eulerian_polynomials_bigint,
};
use polytool::*;
use rmcp::ErrorData as McpError;

/// Exact property verdicts of one family member.
pub(crate) fn family_report(
    index: usize,
    polynomial: NormalizedBigIntPolynomial,
) -> FamilyPolynomialReport {
    let values = polynomial_property_values(&polynomial.coefficients);
    FamilyPolynomialReport {
        index,
        polynomial: polynomial.display,
        real_rooted: values.real_rooted,
        simple_roots: values.simple_roots,
        palindromic: values.palindromic,
        gamma_positive: values.gamma_positive,
        gamma_coefficients: values
            .gamma_coefficients
            .map(|gamma| bigint_strings(&gamma)),
        unimodal: values.unimodal,
        log_concave: values.log_concave,
        ultra_log_concave: values.ultra_log_concave,
    }
}

/// Consecutive interlacing results of the family.
pub(crate) fn family_consecutive_pairs(
    polynomials: &[NormalizedBigIntPolynomial],
) -> Vec<InterlacingPairItem> {
    polynomials
        .windows(2)
        .enumerate()
        .map(|(pair_index, pair)| InterlacingPairItem {
            pair_index,
            left_index: pair_index,
            right_index: pair_index + 1,
            ok: true,
            result: Some(interlacing_result(pair[0].clone(), pair[1].clone())),
            error: None,
        })
        .collect()
}

/// Where the family comes from, after the one-source and `max_n` rules.
enum FamilySource {
    Sequence {
        sequence: SequenceKind,
        max_n: usize,
    },
    Batch,
}

fn select_family_source(
    has_polynomials: bool,
    has_text: bool,
    sequence: Option<&SequenceKind>,
    max_n: Option<usize>,
) -> Result<FamilySource, McpError> {
    let explicit_sources =
        usize::from(has_polynomials) + usize::from(has_text) + usize::from(sequence.is_some());
    if explicit_sources != 1 {
        return Err(invalid_params(
            "expected exactly one source: `polynomials`, `text`, or `sequence`",
        ));
    }

    if let Some(sequence) = sequence {
        let max_n = max_n.ok_or_else(|| invalid_params("`max_n` is required with `sequence`"))?;
        if max_n > MCP_MAX_SEQUENCE_N {
            return Err(invalid_params(format!(
                "max_n must be at most {MCP_MAX_SEQUENCE_N} for MCP sequence generation"
            )));
        }
        return Ok(FamilySource::Sequence {
            sequence: sequence.clone(),
            max_n,
        });
    }

    if max_n.is_some() {
        return Err(invalid_params(
            "`max_n` is only valid when `sequence` is provided",
        ));
    }
    Ok(FamilySource::Batch)
}

fn family_source_label(sequence: &SequenceKind, max_n: usize) -> String {
    format!("sequence:{sequence:?}, max_n={max_n}")
}

/// Family source with exact `BigInt` coefficients: a parsed `polynomials` or
/// `text` batch, or generated rows of a named sequence.
pub(crate) type FamilySourceParse = (
    String,
    Result<Vec<NormalizedBigIntPolynomial>, Vec<ParsePolynomialItem>>,
);

fn collect_bigint_polynomials_or_errors(
    batch: BigIntParsedBatch,
) -> Result<Vec<NormalizedBigIntPolynomial>, Vec<ParsePolynomialItem>> {
    let errors: Vec<ParsePolynomialItem> = batch
        .iter()
        .enumerate()
        .filter_map(|(index, item)| {
            item.as_ref().err().map(|error| ParsePolynomialItem {
                index,
                ok: false,
                polynomial: None,
                coefficients: None,
                degree: None,
                error: Some(error.clone()),
            })
        })
        .collect();
    if errors.is_empty() {
        Ok(batch.into_iter().flatten().collect())
    } else {
        Err(errors)
    }
}

pub(crate) fn family_source_and_polynomials(
    input: &CheckPolynomialFamilyRequest,
) -> Result<FamilySourceParse, McpError> {
    family_source_and_polynomials_from(
        &BigIntPolynomialBatchInput {
            polynomials: input.polynomials.clone(),
            text: input.text.clone(),
        },
        input.sequence.as_ref(),
        input.max_n,
    )
}

fn family_source_and_polynomials_from(
    batch: &BigIntPolynomialBatchInput,
    sequence: Option<&SequenceKind>,
    max_n: Option<usize>,
) -> Result<FamilySourceParse, McpError> {
    match select_family_source(
        batch.polynomials.is_some(),
        batch.text.is_some(),
        sequence,
        max_n,
    )? {
        FamilySource::Sequence { sequence, max_n } => {
            let polynomials = generated_sequence_polynomials_bigint(&sequence, max_n)
                .into_iter()
                .map(normalize_bigint_polynomial)
                .collect();
            Ok((family_source_label(&sequence, max_n), Ok(polynomials)))
        }
        FamilySource::Batch => {
            let parsed = parse_bigint_batch(batch)?;
            let source = if batch.text.is_some() {
                "text".to_string()
            } else {
                "polynomials".to_string()
            };
            Ok((source, collect_bigint_polynomials_or_errors(parsed)))
        }
    }
}

pub(crate) fn generated_sequence_polynomials_bigint(
    sequence: &SequenceKind,
    max_n: usize,
) -> Vec<Vec<BigInt>> {
    match sequence {
        // The library families start at A_0 = N_0 = 1 (equal to A_1, N_1); the
        // tool keeps the classical rows A_1..A_max_n and N_1..N_max_n.
        SequenceKind::Eulerian => eulerian_polynomials_bigint(max_n).split_off(1),
        SequenceKind::Narayana => narayana_polynomials_bigint(max_n).split_off(1),
        SequenceKind::TypeBEulerian => type_b_eulerian_polynomials_bigint(max_n),
        SequenceKind::ChebyshevT => chebyshev_polynomials_t_bigint(max_n),
        SequenceKind::ChebyshevU => chebyshev_polynomials_u_bigint(max_n),
        SequenceKind::Hermite => hermite_polynomials_bigint(max_n),
    }
}

/// The Lace truncation to check, or the reason the request is refused.
enum LacePlan {
    Check {
        block_rows: usize,
        block_cols: usize,
        max_minor_size: usize,
    },
    Refused(LaceCheckResponse),
}

/// Apply the defaults and the MCP size limits to a finite Lace request.
fn plan_lace_check(
    polynomial_count: usize,
    max_degree: usize,
    request: &LaceCheckRequest,
) -> LacePlan {
    let block_rows = request
        .block_rows
        .unwrap_or_else(|| polynomial_count.max(1));
    let block_cols = request
        .block_cols
        .unwrap_or_else(|| block_rows.saturating_add(max_degree).saturating_add(1));
    let refused = |max_minor_size, rows, error: String| {
        LacePlan::Refused(LaceCheckResponse {
            block_rows,
            block_cols,
            max_minor_size,
            rows,
            columns: block_cols,
            totally_nonnegative: false,
            error: Some(error),
            matrix: None,
        })
    };
    let rows = match polynomial_count.checked_mul(block_rows) {
        Some(rows) => rows,
        None => {
            return refused(
                request.max_minor_size.unwrap_or(0),
                0,
                "finite Lace matrix row count overflow".to_string(),
            );
        }
    };
    let max_minor_size = request
        .max_minor_size
        .unwrap_or_else(|| rows.min(block_cols).min(4));
    // Each dimension is bounded on its own: the product bound alone lets
    // `block_cols = 0` through for any number of rows.
    let cells = rows.checked_mul(block_cols);
    let within_cells = rows <= MCP_MAX_LACE_MATRIX_CELLS
        && block_cols <= MCP_MAX_LACE_MATRIX_CELLS
        && cells.is_some_and(|cells| cells <= MCP_MAX_LACE_MATRIX_CELLS);
    // Evaluated only for a matrix within the cell bound, so the number of sizes
    // summed over is small.
    let within_minors = within_cells
        && (1..=max_minor_size.min(rows).min(block_cols))
            .try_fold(0usize, |total, size| {
                binomial_capped(rows, size, MCP_MAX_LACE_MINORS)
                    .checked_mul(binomial_capped(block_cols, size, MCP_MAX_LACE_MINORS))
                    .and_then(|count| total.checked_add(count))
            })
            .is_some_and(|count| count <= MCP_MAX_LACE_MINORS);
    if !within_cells || !within_minors {
        return refused(
            max_minor_size,
            rows,
            format!(
                "finite Lace request exceeds MCP limits ({MCP_MAX_LACE_MATRIX_CELLS} matrix cells and {MCP_MAX_LACE_MINORS} checked minors)"
            ),
        );
    }
    LacePlan::Check {
        block_rows,
        block_cols,
        max_minor_size,
    }
}

/// Finite Lace check: builds the truncation with
/// `lace_matrix_sequence_bigint` and checks all minors up to `max_minor_size`
/// with `check_total_positivity_bigint`.
pub(crate) fn check_family_lace(
    polynomials: &[NormalizedBigIntPolynomial],
    request: &LaceCheckRequest,
) -> LaceCheckResponse {
    let coefficients: Vec<Vec<BigInt>> =
        polynomials.iter().map(|p| p.coefficients.clone()).collect();
    let max_degree = polynomials
        .iter()
        .map(|p| p.display.degree)
        .max()
        .unwrap_or(0);
    let (block_rows, block_cols, max_minor_size) =
        match plan_lace_check(polynomials.len(), max_degree, request) {
            LacePlan::Check {
                block_rows,
                block_cols,
                max_minor_size,
            } => (block_rows, block_cols, max_minor_size),
            LacePlan::Refused(response) => return response,
        };

    match lace_matrix_sequence_bigint(&coefficients, block_rows, block_cols) {
        Ok(matrix) => {
            let rows = matrix.len();
            let columns = matrix.first().map(Vec::len).unwrap_or(0);
            let check = check_total_positivity_bigint(&matrix, max_minor_size, false)
                .map_err(InterlacingMatrixError::NotTotallyNonnegative);
            LaceCheckResponse {
                block_rows,
                block_cols,
                max_minor_size,
                rows,
                columns,
                totally_nonnegative: check.is_ok(),
                error: check.err().map(|err| err.to_string()),
                matrix: request
                    .include_matrix
                    .unwrap_or(false)
                    .then(|| bigint_matrix_strings(&matrix)),
            }
        }
        Err(error) => LaceCheckResponse {
            block_rows,
            block_cols,
            max_minor_size: request.max_minor_size.unwrap_or(0),
            rows: 0,
            columns: 0,
            totally_nonnegative: false,
            error: Some(error.to_string()),
            matrix: None,
        },
    }
}

/// `C(n, k)`, saturating at `cap + 1` so that callers can compare against `cap`.
fn binomial_capped(n: usize, k: usize, cap: usize) -> usize {
    let saturated = cap.saturating_add(1);
    polytool::numeric::binomial_u128(n, k)
        .and_then(|value| usize::try_from(value).ok())
        .map_or(saturated, |value| value.min(saturated))
}

/// Adaptive recurrence search over exact rational rows.
fn family_recurrence_response_rational(
    polynomials: &[Vec<BigRational>],
    search: &AdaptiveSearchOptions,
    budget: AdaptiveSearchBudget,
) -> FindRecurrenceResponse {
    if polynomials.len() < 3 {
        return recurrence_failure_response(
            RecurrenceSearchStatus::InvalidInput,
            "need at least 3 polynomials",
            None,
            Vec::new(),
        );
    }

    match find_recurrence_adaptive_rational_with_budget(polynomials, search, budget) {
        AdaptiveSearchOutcome::Found(result) => FindRecurrenceResponse {
            status: RecurrenceSearchStatus::Found,
            found: true,
            recurrence: Some(format!("{}", result.recurrence)),
            latex: Some(result.recurrence.to_latex()),
            mathematica: Some(
                result
                    .recurrence
                    .to_mathematica_definition_rational(polynomials),
            ),
            sage: Some(result.recurrence.to_sage_definition_rational(polynomials)),
            python: Some(result.recurrence.to_python_definition_rational(polynomials)),
            recurrence_json: Some(recurrence_json_string(
                &result,
                polynomials,
                search,
                polynomials.len(),
            )),
            unknowns: Some(result.num_unknowns),
            weighted_unknowns: Some(result.weighted_unknowns),
            equations: Some(result.num_equations),
            fit_polynomials: Some(result.fit_polynomials),
            verification_polynomials: Some(result.verification_polynomials),
            candidates_tried: Some(result.candidates_tried),
            candidates_considered: Some(result.diagnostics.considered_candidates),
            error: None,
            parse_errors: Vec::new(),
            theorem_matches: crate::theorems::search_theorem_matches(&result, polynomials, search),
        },
        AdaptiveSearchOutcome::NoRecurrence(summary) => recurrence_failure_response(
            RecurrenceSearchStatus::NotFound,
            "no recurrence found within the search bounds",
            Some(summary),
            Vec::new(),
        ),
        AdaptiveSearchOutcome::BudgetExhausted(summary) => recurrence_failure_response(
            RecurrenceSearchStatus::BudgetExhausted,
            format!(
                "recurrence candidate budget exhausted after {} candidates; unsearched candidates remain",
                summary.diagnostics.considered_candidates
            ),
            Some(summary),
            Vec::new(),
        ),
    }
}

/// Recurrence search for the family rows (integer rows of any size, converted
/// with `BigRational::from_integer`).
pub(crate) fn family_recurrence_response(
    polynomials: &[NormalizedBigIntPolynomial],
    search: &AdaptiveSearchOptions,
    budget: AdaptiveSearchBudget,
) -> FindRecurrenceResponse {
    let rows: Vec<Vec<BigRational>> = polynomials
        .iter()
        .map(|polynomial| {
            polynomial
                .coefficients
                .iter()
                .cloned()
                .map(BigRational::from_integer)
                .collect()
        })
        .collect();
    family_recurrence_response_rational(&rows, search, budget)
}

pub(crate) fn first_family_failure(
    items: &[FamilyPolynomialReport],
    pairs: &[InterlacingPairItem],
    lace: Option<&LaceCheckResponse>,
    recurrence: Option<&FindRecurrenceResponse>,
    options: Option<&FamilyCheckOptions>,
) -> Option<String> {
    let require_real_rooted = options.and_then(|o| o.require_real_rooted).unwrap_or(true);
    let require_simple_roots = options
        .and_then(|o| o.require_simple_roots)
        .unwrap_or(false);
    let require_palindromic = options.and_then(|o| o.require_palindromic).unwrap_or(false);
    let require_gamma_positive = options
        .and_then(|o| o.require_gamma_positive)
        .unwrap_or(false);
    let require_unimodal = options.and_then(|o| o.require_unimodal).unwrap_or(false);
    let require_log_concave = options.and_then(|o| o.require_log_concave).unwrap_or(false);
    let require_ultra_log_concave = options
        .and_then(|o| o.require_ultra_log_concave)
        .unwrap_or(false);
    let require_interlacing = options
        .and_then(|o| o.require_consecutive_weak_interlacing)
        .unwrap_or(false);

    for item in items {
        if require_real_rooted && !item.real_rooted {
            return Some(format!("polynomial {} is not real-rooted", item.index));
        }
        if require_simple_roots && !item.simple_roots {
            return Some(format!(
                "polynomial {} does not have simple roots",
                item.index
            ));
        }
        if require_palindromic && !item.palindromic {
            return Some(format!("polynomial {} is not palindromic", item.index));
        }
        if require_gamma_positive && !item.gamma_positive {
            return Some(format!("polynomial {} is not gamma-positive", item.index));
        }
        if require_unimodal && !item.unimodal {
            return Some(format!("polynomial {} is not unimodal", item.index));
        }
        if require_log_concave && !item.log_concave {
            return Some(format!("polynomial {} is not log-concave", item.index));
        }
        if require_ultra_log_concave && !item.ultra_log_concave {
            return Some(format!(
                "polynomial {} is not ultra-log-concave",
                item.index
            ));
        }
    }

    if require_interlacing {
        for pair in pairs {
            if pair.result.as_ref().and_then(|r| r.weak) != Some(true) {
                return Some(format!(
                    "consecutive pair {}-{} is not weakly interlacing",
                    pair.left_index, pair.right_index
                ));
            }
        }
    }

    if let Some(lace) = lace {
        if !lace.totally_nonnegative {
            return Some(format!(
                "finite Lace truncation is not TNN up to minors of size {}",
                lace.max_minor_size
            ));
        }
    }

    if let Some(recurrence) = recurrence {
        if !recurrence.found {
            return Some("recurrence search did not find a recurrence".to_string());
        }
    }

    None
}

fn yes_no(value: bool) -> &'static str {
    if value {
        "yes"
    } else {
        "no"
    }
}

fn opt_yes_no(value: Option<bool>) -> &'static str {
    match value {
        Some(true) => "yes",
        Some(false) => "no",
        None => "n/a",
    }
}

pub(crate) fn family_markdown(
    source: &str,
    all_required_checks_passed: bool,
    first_failure: Option<&str>,
    items: &[FamilyPolynomialReport],
    pairs: &[InterlacingPairItem],
    lace: Option<&LaceCheckResponse>,
    recurrence: Option<&FindRecurrenceResponse>,
) -> String {
    let mut out = String::new();
    out.push_str("# Polynomial family check\n\n");
    out.push_str(&format!("- Source: `{source}`\n"));
    out.push_str(&format!("- Polynomials: {}\n", items.len()));
    out.push_str(&format!(
        "- Required checks: {}\n",
        if all_required_checks_passed {
            "passed"
        } else {
            "failed"
        }
    ));
    if let Some(first_failure) = first_failure {
        out.push_str(&format!("- First failure: {first_failure}\n"));
    }

    out.push_str("\n## Polynomial properties\n\n");
    out.push_str("| i | polynomial | degree | RR | gamma+ | unimodal | LC | ULC | pal |\n");
    out.push_str("|---:|---|---:|:---:|:---:|:---:|:---:|:---:|:---:|\n");
    for item in items {
        out.push_str(&format!(
            "| {} | `{}` | {} | {} | {} | {} | {} | {} | {} |\n",
            item.index,
            item.polynomial.polynomial,
            item.polynomial.degree,
            yes_no(item.real_rooted),
            yes_no(item.gamma_positive),
            yes_no(item.unimodal),
            yes_no(item.log_concave),
            yes_no(item.ultra_log_concave),
            yes_no(item.palindromic)
        ));
    }

    if !pairs.is_empty() {
        out.push_str("\n## Consecutive interlacing\n\n");
        out.push_str("| pair | weak | strict | status |\n");
        out.push_str("|---:|:---:|:---:|---|\n");
        for pair in pairs {
            if let Some(result) = &pair.result {
                out.push_str(&format!(
                    "| {}-{} | {} | {} | `{}` |\n",
                    pair.left_index,
                    pair.right_index,
                    opt_yes_no(result.weak),
                    opt_yes_no(result.strict),
                    result.status
                ));
            } else {
                out.push_str(&format!(
                    "| {}-{} | n/a | n/a | `{}` |\n",
                    pair.left_index,
                    pair.right_index,
                    pair.error.as_deref().unwrap_or("parse error")
                ));
            }
        }
    }

    if let Some(lace) = lace {
        out.push_str("\n## Finite Lace check\n\n");
        out.push_str(&format!(
            "- Truncation: block rows `{}`, block cols `{}`; matrix size `{} x {}`.\n",
            lace.block_rows, lace.block_cols, lace.rows, lace.columns
        ));
        out.push_str(&format!(
            "- Checked minors up to size `{}`: {}.\n",
            lace.max_minor_size,
            if lace.totally_nonnegative {
                "TNN"
            } else {
                "not TNN"
            }
        ));
        if let Some(error) = &lace.error {
            out.push_str(&format!("- Failure: `{error}`\n"));
        }
    }

    if let Some(recurrence) = recurrence {
        out.push_str("\n## Recurrence search\n\n");
        if let Some(found) = &recurrence.recurrence {
            out.push_str(&format!("- Found: `{found}`\n"));
        } else if let Some(error) = &recurrence.error {
            out.push_str(&format!("- Not found: {error}\n"));
        } else {
            out.push_str("- Not found.\n");
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence_options::apply_recurrence_options;
    use crate::requests::BigIntPolynomialInput;

    /// Rows `(1 + scale x)^n`, `n = first..=last`, as `BigInt` coefficient vectors.
    fn scaled_binomial_rows(scale: i64, first: usize, last: usize) -> Vec<Vec<BigInt>> {
        (first..=last)
            .map(|n| {
                let mut binomial = BigInt::from(1);
                let mut power = BigInt::from(1);
                let mut row = Vec::new();
                for k in 0..=n {
                    row.push(&binomial * &power);
                    binomial = binomial * BigInt::from(n - k) / BigInt::from(k + 1);
                    power *= scale;
                }
                row
            })
            .collect()
    }

    fn normalized(rows: &[Vec<BigInt>]) -> Vec<NormalizedBigIntPolynomial> {
        rows.iter()
            .cloned()
            .map(normalize_bigint_polynomial)
            .collect()
    }

    fn text_batch(rows: &[Vec<BigInt>]) -> BigIntPolynomialBatchInput {
        let text = rows
            .iter()
            .map(|row| {
                row.iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            })
            .collect::<Vec<_>>()
            .join("\n");
        BigIntPolynomialBatchInput {
            polynomials: None,
            text: Some(text),
        }
    }

    fn i64_rows(rows: &[Vec<BigInt>]) -> Vec<Vec<i64>> {
        rows.iter()
            .map(|row| row.iter().map(|v| i64::try_from(v).unwrap()).collect())
            .collect()
    }

    #[test]
    fn bigint_source_accepts_rows_beyond_i64_from_every_source() {
        // (1 + x)^70 has the middle coefficient C(70, 35) > i64::MAX.
        let rows = scaled_binomial_rows(1, 68, 70);
        assert!(rows[2][35] > BigInt::from(i64::MAX));

        let (source, parsed) =
            family_source_and_polynomials_from(&text_batch(&rows), None, None).unwrap();
        assert_eq!(source, "text");
        assert_eq!(parsed.unwrap(), normalized(&rows));

        let listed = BigIntPolynomialBatchInput {
            polynomials: Some(
                rows.iter()
                    .map(|row| BigIntPolynomialInput {
                        coefficients: Some(
                            row.iter()
                                .map(|value| {
                                    crate::requests::BigIntCoefficientInput::Text(value.to_string())
                                })
                                .collect(),
                        ),
                        expression: None,
                    })
                    .collect(),
            ),
            text: None,
        };
        let (source, parsed) = family_source_and_polynomials_from(&listed, None, None).unwrap();
        assert_eq!(source, "polynomials");
        assert_eq!(parsed.unwrap(), normalized(&rows));

        // Eulerian A_30 has leading coefficient 30! > i64::MAX; the sequence
        // source no longer narrows to i64.
        let empty = BigIntPolynomialBatchInput {
            polynomials: None,
            text: None,
        };
        let (source, parsed) =
            family_source_and_polynomials_from(&empty, Some(&SequenceKind::Eulerian), Some(30))
                .unwrap();
        assert_eq!(source, "sequence:Eulerian, max_n=30");
        let parsed = parsed.unwrap();
        assert_eq!(parsed.len(), 30);
        assert!(parsed[29]
            .coefficients
            .iter()
            .any(|value| *value > BigInt::from(i64::MAX)));
    }

    #[test]
    fn bigint_source_keeps_the_source_rules_and_parse_errors() {
        let error = |batch: &BigIntPolynomialBatchInput, sequence, max_n| {
            family_source_and_polynomials_from(batch, sequence, max_n)
                .expect_err("rejected")
                .message
                .to_string()
        };
        let empty = BigIntPolynomialBatchInput {
            polynomials: None,
            text: None,
        };
        assert!(error(&empty, None, None).contains("expected exactly one source"));
        let both = BigIntPolynomialBatchInput {
            polynomials: None,
            text: Some("1,1".to_string()),
        };
        assert!(error(&both, Some(&SequenceKind::Hermite), Some(3))
            .contains("expected exactly one source"));
        assert!(error(&empty, Some(&SequenceKind::Hermite), None).contains("`max_n` is required"));
        assert!(error(
            &empty,
            Some(&SequenceKind::Hermite),
            Some(MCP_MAX_SEQUENCE_N + 1)
        )
        .contains("max_n must be at most"));
        assert!(error(&both, None, Some(3)).contains("only valid when `sequence` is provided"));

        let (_, parsed) = family_source_and_polynomials_from(
            &BigIntPolynomialBatchInput {
                polynomials: None,
                text: Some("1,1\nbad + x\n1,2,1\nalso bad".to_string()),
            },
            None,
            None,
        )
        .unwrap();
        let errors = parsed.expect_err("parse errors");
        assert_eq!(
            errors.iter().map(|item| item.index).collect::<Vec<_>>(),
            [1, 3]
        );
        assert!(errors.iter().all(|item| !item.ok && item.error.is_some()));
    }

    #[test]
    fn reports_agree_with_the_i64_library_routines_for_small_rows() {
        let rows: Vec<Vec<i64>> = vec![
            vec![1, 4, 1],
            vec![1, 11, 11, 1],
            vec![0, 1, 3, 1],
            vec![1, 2, 1],
            vec![1, 0, 1],
            vec![2, 3, 5],
        ];
        for (index, row) in rows.iter().enumerate() {
            let big: Vec<BigInt> = row.iter().copied().map(BigInt::from).collect();
            let report = family_report(index, normalize_bigint_polynomial(big));
            assert_eq!(report.index, index);
            assert_eq!(report.polynomial.polynomial, format_poly(row));
            assert_eq!(
                (
                    report.real_rooted,
                    report.simple_roots,
                    report.palindromic,
                    report.gamma_positive
                ),
                (
                    is_real_rooted(row),
                    has_simple_roots(row),
                    is_palindromic_ignoring_initial_zeros(row),
                    is_gamma_positive_ignoring_initial_zeros(row)
                ),
                "{row:?}"
            );
            assert_eq!(
                (
                    report.unimodal,
                    report.log_concave,
                    report.ultra_log_concave
                ),
                (
                    is_unimodal(row),
                    is_log_concave(row),
                    is_ultra_log_concave(row)
                )
            );
            assert_eq!(
                report.gamma_coefficients,
                gamma_coefficients_ignoring_initial_zeros(row)
                    .map(|gamma| gamma.iter().map(ToString::to_string).collect()),
                "{row:?}"
            );
        }
    }

    #[test]
    fn bigint_reports_and_pairs_are_exact_beyond_i64() {
        let rows = scaled_binomial_rows(1, 67, 68);
        let polynomials = normalized(&rows);
        let reports: Vec<_> = polynomials
            .iter()
            .cloned()
            .enumerate()
            .map(|(index, polynomial)| family_report(index, polynomial))
            .collect();
        for (index, report) in reports.iter().enumerate() {
            let degree = 67 + index;
            assert!(report.real_rooted && report.palindromic && report.gamma_positive);
            assert!(!report.simple_roots);
            // (1 + x)^d = 1 * (1 + x)^d: gamma = (1, 0, ..., 0).
            let gamma = report.gamma_coefficients.as_ref().unwrap();
            assert_eq!(gamma.len(), degree / 2 + 1);
            assert_eq!(gamma[0], "1");
            assert!(gamma[1..].iter().all(|value| value == "0"));
        }
        assert_eq!(
            reports[1].polynomial.coefficients[34],
            rows[1][34].to_string()
        );

        let pairs = family_consecutive_pairs(&polynomials);
        assert_eq!(pairs.len(), 1);
        assert!(pairs.iter().all(|pair| pair.ok));
        // Scaling x by 10^10 does not change the interlacing verdicts.
        let small_pairs = family_consecutive_pairs(&normalized(&scaled_binomial_rows(1, 3, 5)));
        let scaled = scaled_binomial_rows(10_000_000_000, 3, 5);
        assert!(scaled[2]
            .iter()
            .any(|value| *value > BigInt::from(i64::MAX)));
        let scaled_pairs = family_consecutive_pairs(&normalized(&scaled));
        for (pair, scaled_pair) in small_pairs.iter().zip(&scaled_pairs) {
            let (pair, scaled_pair) = (
                pair.result.as_ref().unwrap(),
                scaled_pair.result.as_ref().unwrap(),
            );
            assert_eq!(
                (pair.strict, pair.weak),
                (scaled_pair.strict, scaled_pair.weak)
            );
            assert_eq!(pair.status, scaled_pair.status);
        }
    }

    #[test]
    fn bigint_gamma_coefficients_are_exact_beyond_i64() {
        // 1 + t^100 has coefficients 0 and 1, but its gamma coefficients
        // (-1)^i 100/(100-i) C(100-i, i) reach 10^20; the former i64 route
        // reported none.
        let mut row = vec![BigInt::from(0); 101];
        row[0] = BigInt::from(1);
        row[100] = BigInt::from(1);
        let report = family_report(0, normalize_bigint_polynomial(row));
        assert!(report.palindromic && !report.gamma_positive && !report.real_rooted);
        let gamma = report.gamma_coefficients.unwrap();
        assert_eq!(gamma.len(), 51);
        assert_eq!(&gamma[..3], ["1", "-100", "4850"]);
        assert_eq!(gamma[50], "2");
        assert!(gamma.iter().any(
            |value| value.trim_start_matches('-').parse::<BigInt>().unwrap()
                > BigInt::from(i64::MAX)
        ));

        // A_4 = 1 + 11 t + 11 t^2 + t^3 has gamma coefficients 1, 8 (from the
        // sequence source, whose rows are exact beyond A_21).
        let empty = BigIntPolynomialBatchInput {
            polynomials: None,
            text: None,
        };
        let (_, parsed) =
            family_source_and_polynomials_from(&empty, Some(&SequenceKind::Eulerian), Some(30))
                .unwrap();
        let report = family_report(3, parsed.unwrap()[3].clone());
        assert!(report.real_rooted && report.palindromic && report.gamma_positive);
        assert_eq!(
            report.gamma_coefficients,
            Some(vec!["1".to_string(), "8".to_string()])
        );
    }

    /// The finite Lace check through the `i64` library routines.
    fn lace_expected_i64(
        rows: &[Vec<BigInt>],
        block_rows: usize,
        block_cols: usize,
        size: usize,
    ) -> LaceCheckResponse {
        let small = i64_rows(rows);
        let matrix = lace_matrix_sequence_i64(&small, block_rows, block_cols).unwrap();
        let check = check_lace_sequence_totally_nonnegative_by_minors_i64(
            &small, block_rows, block_cols, size,
        );
        LaceCheckResponse {
            block_rows,
            block_cols,
            max_minor_size: size,
            rows: matrix.len(),
            columns: matrix[0].len(),
            totally_nonnegative: check.is_ok(),
            error: check.err().map(|err| err.to_string()),
            matrix: Some(
                matrix
                    .iter()
                    .map(|row| row.iter().map(ToString::to_string).collect())
                    .collect(),
            ),
        }
    }

    #[test]
    fn lace_check_agrees_with_the_i64_library_routines_for_small_rows() {
        let rows = scaled_binomial_rows(1, 1, 3);
        let request = |block_rows, block_cols, max_minor_size| LaceCheckRequest {
            block_rows: Some(block_rows),
            block_cols: Some(block_cols),
            max_minor_size: Some(max_minor_size),
            include_matrix: Some(true),
        };
        for (block_rows, block_cols, size) in [(3, 7, 3), (2, 5, 2), (1, 4, 1)] {
            assert_eq!(
                check_family_lace(&normalized(&rows), &request(block_rows, block_cols, size)),
                lace_expected_i64(&rows, block_rows, block_cols, size),
            );
        }
        // Defaults: one block row per polynomial, `block_rows + max_degree + 1` columns.
        let defaults = check_family_lace(
            &normalized(&rows),
            &LaceCheckRequest {
                block_rows: None,
                block_cols: None,
                max_minor_size: None,
                include_matrix: None,
            },
        );
        assert_eq!((defaults.block_rows, defaults.block_cols), (3, 7));
        assert_eq!(defaults.max_minor_size, 4);
        assert_eq!(defaults.matrix, None);
        // Over the MCP limits: refused before any matrix is built.
        let refused = check_family_lace(&normalized(&rows), &request(1_000_000, 3, 1));
        assert!(!refused.totally_nonnegative);
        assert!(refused.error.unwrap().contains("exceeds MCP limits"));
        // A non-TNN truncation reports the library's failure text.
        let bad = vec![
            vec![BigInt::from(1), BigInt::from(-1)],
            vec![BigInt::from(1)],
        ];
        let response = check_family_lace(&normalized(&bad), &request(2, 3, 2));
        assert!(!response.totally_nonnegative);
        assert_eq!(response, lace_expected_i64(&bad, 2, 3, 2));
    }

    #[test]
    fn bigint_lace_check_handles_entries_beyond_i64() {
        // Scaling x by 10^10 multiplies rows and columns of the Lace matrix by
        // positive factors, so TNN verdicts agree with the unscaled rows.
        let unscaled = scaled_binomial_rows(1, 1, 3);
        let scaled = scaled_binomial_rows(10_000_000_000, 1, 3);
        assert!(scaled[2]
            .iter()
            .any(|value| *value > BigInt::from(i64::MAX)));
        let request = LaceCheckRequest {
            block_rows: None,
            block_cols: None,
            max_minor_size: Some(3),
            include_matrix: Some(true),
        };
        let expected = check_family_lace(&normalized(&unscaled), &request);
        let response = check_family_lace(&normalized(&scaled), &request);
        // Both truncations are non-TNN, with a 2 x 2 minor that scales by 10^30.
        assert!(!expected.totally_nonnegative && !response.totally_nonnegative);
        assert!(expected.error.unwrap().contains("has det = -1 < 0"));
        assert!(response
            .error
            .unwrap()
            .contains("has det = -1000000000000000000000000000000 < 0"));
        assert_eq!(
            (response.rows, response.columns),
            (expected.rows, expected.columns)
        );
        let matrix = response.matrix.unwrap();
        assert!(matrix
            .iter()
            .flatten()
            .any(|value| value.parse::<BigInt>().unwrap() > BigInt::from(i64::MAX)));

        // The Toeplitz truncation of the single polynomial (1 + 10^10 x)^3 is TNN.
        let single = normalized(&scaled[2..]);
        let response = check_family_lace(&single, &request);
        assert_eq!(response.error, None, "{response:?}");
        assert!(response.totally_nonnegative);

        // A huge negative entry makes a 1 x 1 minor negative.
        let mut negative = scaled[2..].to_vec();
        negative[0][2] = -negative[0][2].clone();
        let response = check_family_lace(&normalized(&negative), &request);
        assert!(!response.totally_nonnegative);
        assert!(response.error.unwrap().contains("< 0"));
    }

    #[test]
    fn bigint_recurrence_search_handles_rows_beyond_i64() {
        let (search, budget) = apply_recurrence_options(None).unwrap();
        let scaled = scaled_binomial_rows(10_000_000_000, 0, 6);
        assert!(scaled[6]
            .iter()
            .any(|value| *value > BigInt::from(i64::MAX)));
        let response = family_recurrence_response(&normalized(&scaled), &search, budget);
        assert!(response.found, "{:?}", response.error);
        assert_eq!(response.status, RecurrenceSearchStatus::Found);
        assert!(response.recurrence_json.is_some() && response.python.is_some());
        let recurrence = response.recurrence.unwrap();
        // P_n = (1 + 10^10 x) P_{n-1}.
        assert!(
            recurrence.contains("P(n-1)") && recurrence.contains("10000000000"),
            "{recurrence}"
        );

        let few = family_recurrence_response(&normalized(&scaled[..2]), &search, budget);
        assert_eq!(few.status, RecurrenceSearchStatus::InvalidInput);
    }

    #[test]
    fn recurrence_response_agrees_with_the_i64_library_search_for_small_rows() {
        let (search, budget) = apply_recurrence_options(None).unwrap();
        let rows = scaled_binomial_rows(1, 0, 7);
        let expected = polytool::recurrence::find_recurrence_adaptive(&i64_rows(&rows), &search)
            .expect("Pascal rows satisfy a recurrence");
        let response = family_recurrence_response(&normalized(&rows), &search, budget);
        assert!(response.found);
        assert_eq!(response.recurrence, Some(expected.recurrence.to_string()));
        assert_eq!(response.unknowns, Some(expected.num_unknowns));
        assert_eq!(response.equations, Some(expected.num_equations));
        assert_eq!(response.candidates_tried, Some(expected.candidates_tried));
    }
    #[test]
    fn binomial_capped_saturates_at_cap_plus_one() {
        assert_eq!(binomial_capped(10, 3, 1000), 120);
        assert_eq!(binomial_capped(10, 3, 120), 120);
        assert_eq!(binomial_capped(10, 3, 119), 120);
        assert_eq!(binomial_capped(10, 3, 10), 11);
        assert_eq!(binomial_capped(3, 5, 10), 0);
        assert_eq!(binomial_capped(200, 100, 250_000), 250_001);
        assert_eq!(binomial_capped(usize::MAX, 3, usize::MAX), usize::MAX);
    }
}
