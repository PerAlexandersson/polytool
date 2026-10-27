//! The `check_polynomial_family` workflow: family sources, per-polynomial
//! reports, finite Lace TNN checks, recurrence search and the Markdown summary.

use crate::limits::{MCP_MAX_LACE_MATRIX_CELLS, MCP_MAX_LACE_MINORS, MCP_MAX_SEQUENCE_N};
use crate::parsing::{
    collect_polynomials_or_errors, integer_polys_to_rational, invalid_params, normalize_polynomial,
    parse_batch, FamilySourceParse,
};
use crate::recurrence_options::{recurrence_failure_response, recurrence_json_string};
use crate::requests::{
    CheckPolynomialFamilyRequest, FamilyCheckOptions, LaceCheckRequest, PolynomialBatchInput,
    SequenceKind,
};
use crate::responses::{
    FamilyPolynomialReport, FindRecurrenceResponse, InterlacingPairItem, LaceCheckResponse,
    NormalizedPolynomial, RecurrenceSearchStatus,
};
use num_bigint::BigInt;
use polytool::recurrence::{
    find_recurrence_adaptive_with_budget, AdaptiveSearchBudget, AdaptiveSearchOptions,
    AdaptiveSearchOutcome,
};
use polytool::sequences::{
    chebyshev_polynomials_t_bigint, chebyshev_polynomials_u_bigint, eulerian_polynomials_bigint,
    hermite_polynomials_bigint, narayana_polynomials_bigint, type_b_eulerian_polynomials_bigint,
};
use polytool::*;
use rmcp::ErrorData as McpError;

pub(crate) fn family_report(
    index: usize,
    polynomial: NormalizedPolynomial,
) -> FamilyPolynomialReport {
    let coefficients = &polynomial.coefficients;
    FamilyPolynomialReport {
        index,
        real_rooted: is_real_rooted(coefficients),
        simple_roots: has_simple_roots(coefficients),
        palindromic: is_palindromic_ignoring_initial_zeros(coefficients),
        gamma_positive: is_gamma_positive_ignoring_initial_zeros(coefficients),
        gamma_coefficients: gamma_coefficients_ignoring_initial_zeros(coefficients),
        unimodal: is_unimodal(coefficients),
        log_concave: is_log_concave(coefficients),
        ultra_log_concave: is_ultra_log_concave(coefficients),
        polynomial,
    }
}

pub(crate) fn family_source_and_polynomials(
    input: &CheckPolynomialFamilyRequest,
) -> Result<FamilySourceParse, McpError> {
    let explicit_sources = usize::from(input.polynomials.is_some())
        + usize::from(input.text.is_some())
        + usize::from(input.sequence.is_some());
    if explicit_sources != 1 {
        return Err(invalid_params(
            "expected exactly one source: `polynomials`, `text`, or `sequence`",
        ));
    }

    if let Some(sequence) = input.sequence.clone() {
        let max_n = input
            .max_n
            .ok_or_else(|| invalid_params("`max_n` is required with `sequence`"))?;
        if max_n > MCP_MAX_SEQUENCE_N {
            return Err(invalid_params(format!(
                "max_n must be at most {MCP_MAX_SEQUENCE_N} for MCP sequence generation"
            )));
        }
        let bigint_polynomials = generated_sequence_polynomials_bigint(&sequence, max_n);
        let polynomials = bigint_polynomials
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|coefficient| {
                        i64::try_from(&coefficient).map_err(|_| {
                            invalid_params(format!(
                                "sequence {sequence:?} at max_n={max_n} exceeds the i64 coefficient range supported by check_polynomial_family"
                            ))
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()
                    .map(normalize_polynomial)
            })
            .collect::<Result<Vec<_>, _>>()?;
        return Ok((
            format!("sequence:{sequence:?}, max_n={max_n}"),
            Ok(polynomials),
        ));
    }

    if input.max_n.is_some() {
        return Err(invalid_params(
            "`max_n` is only valid when `sequence` is provided",
        ));
    }

    let batch = parse_batch(&PolynomialBatchInput {
        polynomials: input.polynomials.clone(),
        text: input.text.clone(),
    })?;
    let source = if input.text.is_some() {
        "text".to_string()
    } else {
        "polynomials".to_string()
    };
    Ok((source, collect_polynomials_or_errors(batch)))
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

pub(crate) fn check_family_lace(
    polynomials: &[NormalizedPolynomial],
    request: &LaceCheckRequest,
) -> LaceCheckResponse {
    let coefficients: Vec<Vec<i64>> = polynomials.iter().map(|p| p.coefficients.clone()).collect();
    let max_degree = polynomials.iter().map(|p| p.degree).max().unwrap_or(0);
    let block_rows = request
        .block_rows
        .unwrap_or_else(|| polynomials.len().max(1));
    let block_cols = request
        .block_cols
        .unwrap_or_else(|| block_rows + max_degree + 1);
    let rows = match polynomials.len().checked_mul(block_rows) {
        Some(rows) => rows,
        None => {
            return LaceCheckResponse {
                block_rows,
                block_cols,
                max_minor_size: request.max_minor_size.unwrap_or(0),
                rows: 0,
                columns: block_cols,
                totally_nonnegative: false,
                error: Some("finite Lace matrix row count overflow".to_string()),
                matrix: None,
            };
        }
    };
    let max_minor_size = request
        .max_minor_size
        .unwrap_or_else(|| rows.min(block_cols).min(4));
    let cells = rows.checked_mul(block_cols);
    let minor_count =
        (1..=max_minor_size.min(rows).min(block_cols)).try_fold(0usize, |total, size| {
            binomial_capped(rows, size, MCP_MAX_LACE_MINORS)
                .checked_mul(binomial_capped(block_cols, size, MCP_MAX_LACE_MINORS))
                .and_then(|count| total.checked_add(count))
        });
    if cells.is_none_or(|cells| cells > MCP_MAX_LACE_MATRIX_CELLS)
        || minor_count.is_none_or(|count| count > MCP_MAX_LACE_MINORS)
    {
        return LaceCheckResponse {
            block_rows,
            block_cols,
            max_minor_size,
            rows,
            columns: block_cols,
            totally_nonnegative: false,
            error: Some(format!(
                "finite Lace request exceeds MCP limits ({MCP_MAX_LACE_MATRIX_CELLS} matrix cells and {MCP_MAX_LACE_MINORS} checked minors)"
            )),
            matrix: None,
        };
    }

    match lace_matrix_sequence_i64(&coefficients, block_rows, block_cols) {
        Ok(matrix) => {
            let rows = matrix.len();
            let columns = matrix.first().map(Vec::len).unwrap_or(0);
            let check = check_lace_sequence_totally_nonnegative_by_minors_i64(
                &coefficients,
                block_rows,
                block_cols,
                max_minor_size,
            );
            LaceCheckResponse {
                block_rows,
                block_cols,
                max_minor_size,
                rows,
                columns,
                totally_nonnegative: check.is_ok(),
                error: check.err().map(|err| err.to_string()),
                matrix: request.include_matrix.unwrap_or(false).then_some(matrix),
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

pub(crate) fn family_recurrence_response(
    coefficients: &[Vec<i64>],
    search: &AdaptiveSearchOptions,
    budget: AdaptiveSearchBudget,
) -> FindRecurrenceResponse {
    if coefficients.len() < 3 {
        return recurrence_failure_response(
            RecurrenceSearchStatus::InvalidInput,
            "need at least 3 polynomials",
            None,
            Vec::new(),
        );
    }

    match find_recurrence_adaptive_with_budget(coefficients, search, budget) {
        AdaptiveSearchOutcome::Found(result) => {
            let rational_polys = integer_polys_to_rational(coefficients);
            FindRecurrenceResponse {
                status: RecurrenceSearchStatus::Found,
                found: true,
                recurrence: Some(format!("{}", result.recurrence)),
                latex: Some(result.recurrence.to_latex()),
                mathematica: Some(result.recurrence.to_mathematica_definition(coefficients)),
                sage: Some(result.recurrence.to_sage_definition(coefficients)),
                python: Some(result.recurrence.to_python_definition(coefficients)),
                recurrence_json: Some(recurrence_json_string(
                    &result,
                    &rational_polys,
                    search,
                    coefficients.len(),
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
                theorem_matches: crate::theorems::search_theorem_matches(
                    &result,
                    &rational_polys,
                    search,
                ),
            }
        }
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
