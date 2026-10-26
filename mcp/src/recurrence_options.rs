//! Recurrence-search options: validation against the MCP limits, defaults,
//! budgets, and the JSON rendering of search results and failures.

use crate::limits::{
    MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL, MCP_MAX_RECURRENCE_CANDIDATES, MCP_MAX_RECURRENCE_DEGREE,
    MCP_MAX_RECURRENCE_INPUT_ROWS, MCP_MAX_RECURRENCE_LENGTH, MCP_MAX_RECURRENCE_UNKNOWNS,
    MCP_MAX_TOTAL_COEFFICIENTS,
};
use crate::parsing::invalid_params;
use crate::requests::RecurrenceSearchOptionsInput;
use crate::responses::{FindRecurrenceResponse, ParsePolynomialItem, RecurrenceSearchStatus};
use polytool::recurrence::{
    candidate_complexity, AdaptiveSearchBudget, AdaptiveSearchOptions, AdaptiveSearchResult,
    AdaptiveSearchSummary, BigRational, BivarPoly, Recurrence, RecurrenceJson,
    RecurrenceJsonSearch, RecurrenceOptions, RecurrenceOptionsJson,
};
use rmcp::ErrorData as McpError;

pub(crate) fn recurrence_json_string(
    result: &AdaptiveSearchResult,
    polynomials: &[Vec<BigRational>],
    search: &AdaptiveSearchOptions,
    source_rows: usize,
) -> String {
    let searched_polys = polynomials.get(search.skip_prefix..).unwrap_or(&[]);
    let initial_count = result
        .recurrence
        .generation_initial_count(1, searched_polys.len());
    let initial_polys = &searched_polys[..initial_count];
    let recurrence_json = RecurrenceJson::from_recurrence_rational(
        &result.recurrence,
        1,
        initial_polys,
        Some(RecurrenceJsonSearch {
            recurrence_text: result.recurrence.to_string(),
            source_rows,
            skip_prefix: search.skip_prefix,
            unknowns: result.num_unknowns,
            weighted_unknowns: result.weighted_unknowns,
            equations: result.num_equations,
            fit_polynomials: result.fit_polynomials,
            verification_polynomials: result.verification_polynomials,
            candidates_tried: result.candidates_tried,
            options: RecurrenceOptionsJson::from(&result.opts),
        }),
    );
    serde_json::to_string_pretty(&recurrence_json).expect("serialize recurrence JSON")
}

pub(crate) fn recurrence_failure_response(
    status: RecurrenceSearchStatus,
    error: impl Into<String>,
    summary: Option<AdaptiveSearchSummary>,
    parse_errors: Vec<ParsePolynomialItem>,
) -> FindRecurrenceResponse {
    let (candidates_tried, candidates_considered) = summary.map_or((None, None), |summary| {
        (
            Some(summary.candidates_tried),
            Some(summary.diagnostics.considered_candidates),
        )
    });
    FindRecurrenceResponse {
        status,
        found: false,
        recurrence: None,
        latex: None,
        mathematica: None,
        sage: None,
        python: None,
        recurrence_json: None,
        unknowns: None,
        weighted_unknowns: None,
        equations: None,
        fit_polynomials: None,
        verification_polynomials: None,
        candidates_tried,
        candidates_considered,
        error: Some(error.into()),
        parse_errors,
        theorem_matches: None,
    }
}

fn inclusive_option_count(minimum: usize, maximum: usize) -> Result<usize, McpError> {
    maximum
        .checked_sub(minimum)
        .and_then(|difference| difference.checked_add(1))
        .ok_or_else(|| invalid_params(format!("invalid recurrence range {minimum}..={maximum}")))
}

fn validate_recurrence_options(
    options: &AdaptiveSearchOptions,
    budget: AdaptiveSearchBudget,
) -> Result<(), McpError> {
    let ranges = [
        (
            "recurrence length",
            options.min_rec_len,
            options.max_rec_len,
        ),
        ("variable degree", options.min_var_deg, options.max_var_deg),
        ("index degree", options.min_idx_deg, options.max_idx_deg),
        (
            "derivative degree",
            options.min_diff_deg,
            options.max_diff_deg,
        ),
        (
            "inhomogeneous variable degree",
            options.min_inhomo_var_deg,
            options.max_inhomo_var_deg,
        ),
        (
            "inhomogeneous index degree",
            options.min_inhomo_idx_deg,
            options.max_inhomo_idx_deg,
        ),
    ];
    for (name, minimum, maximum) in ranges {
        if minimum > maximum {
            return Err(invalid_params(format!(
                "minimum {name} {minimum} exceeds maximum {maximum}"
            )));
        }
    }
    if options.max_rec_len > MCP_MAX_RECURRENCE_LENGTH {
        return Err(invalid_params(format!(
            "maximum recurrence length {} exceeds the MCP limit {MCP_MAX_RECURRENCE_LENGTH}",
            options.max_rec_len
        )));
    }
    let degree_bounds = [
        ("variable", options.max_var_deg),
        ("index", options.max_idx_deg),
        ("derivative", options.max_diff_deg),
        ("inhomogeneous variable", options.max_inhomo_var_deg),
        ("inhomogeneous index", options.max_inhomo_idx_deg),
        ("denominator variable", options.max_denom_var_deg),
        ("denominator index", options.max_denom_idx_deg),
    ];
    for (name, degree) in degree_bounds {
        if degree > MCP_MAX_RECURRENCE_DEGREE {
            return Err(invalid_params(format!(
                "maximum {name} degree {degree} exceeds the MCP limit {MCP_MAX_RECURRENCE_DEGREE}"
            )));
        }
    }
    if options.skip_prefix > MCP_MAX_RECURRENCE_INPUT_ROWS
        || options.fit_extra_rows > MCP_MAX_RECURRENCE_INPUT_ROWS
    {
        return Err(invalid_params(format!(
            "skip_prefix and fit_extra_rows must be at most {MCP_MAX_RECURRENCE_INPUT_ROWS}"
        )));
    }
    if options.min_margin > MCP_MAX_RECURRENCE_UNKNOWNS {
        return Err(invalid_params(format!(
            "min_margin {} exceeds the MCP limit {MCP_MAX_RECURRENCE_UNKNOWNS}",
            options.min_margin
        )));
    }

    let base_count = inclusive_option_count(options.min_rec_len.max(1), options.max_rec_len)?
        .checked_mul(inclusive_option_count(
            options.min_diff_deg,
            options.max_diff_deg,
        )?)
        .and_then(|count| {
            count
                .checked_mul(inclusive_option_count(options.min_idx_deg, options.max_idx_deg).ok()?)
        })
        .and_then(|count| {
            count
                .checked_mul(inclusive_option_count(options.min_var_deg, options.max_var_deg).ok()?)
        })
        .and_then(|count| count.checked_mul(if options.try_alternating_sign { 2 } else { 1 }))
        .ok_or_else(|| invalid_params("recurrence candidate count overflow"))?;
    let inhomogeneous_count = if options.try_inhomogeneous {
        inclusive_option_count(options.min_inhomo_idx_deg, options.max_inhomo_idx_deg)?
            .checked_mul(inclusive_option_count(
                options.min_inhomo_var_deg,
                options.max_inhomo_var_deg,
            )?)
            .ok_or_else(|| invalid_params("inhomogeneous candidate count overflow"))?
    } else {
        0
    };
    let denominator_count = if options.try_denominator {
        options
            .max_denom_idx_deg
            .checked_add(1)
            .and_then(|count| count.checked_mul(options.max_denom_var_deg.checked_add(1)?))
            .and_then(|count| count.checked_sub(1))
            .ok_or_else(|| invalid_params("denominator candidate count overflow"))?
    } else {
        0
    };
    let candidate_count = 1usize
        .checked_add(inhomogeneous_count)
        .and_then(|count| count.checked_add(denominator_count))
        .and_then(|count| count.checked_mul(base_count))
        .ok_or_else(|| invalid_params("recurrence candidate count overflow"))?;
    let effective_candidate_count = budget
        .max_candidates
        .map_or(candidate_count, |limit| candidate_count.min(limit));
    if effective_candidate_count > MCP_MAX_RECURRENCE_CANDIDATES {
        return Err(invalid_params(format!(
            "recurrence search may inspect {effective_candidate_count} candidates; the MCP limit is {MCP_MAX_RECURRENCE_CANDIDATES}"
        )));
    }

    let base_options = RecurrenceOptions {
        rec_len: options.max_rec_len,
        var_deg: options.max_var_deg,
        idx_deg: options.max_idx_deg,
        diff_deg: options.max_diff_deg,
        homogeneous: true,
        inhomo_var_deg: 0,
        inhomo_idx_deg: 0,
        denom_var_deg: 0,
        denom_idx_deg: 0,
        alternating_sign: options.try_alternating_sign,
        modular_prefilter: options.modular_prefilter,
    };
    let mut maximum_unknowns = candidate_complexity(&base_options).raw_unknowns;
    if options.try_inhomogeneous {
        maximum_unknowns = maximum_unknowns.max(
            candidate_complexity(&RecurrenceOptions {
                homogeneous: false,
                inhomo_var_deg: options.max_inhomo_var_deg,
                inhomo_idx_deg: options.max_inhomo_idx_deg,
                ..base_options.clone()
            })
            .raw_unknowns,
        );
    }
    if options.try_denominator {
        maximum_unknowns = maximum_unknowns.max(
            candidate_complexity(&RecurrenceOptions {
                denom_var_deg: options.max_denom_var_deg,
                denom_idx_deg: options.max_denom_idx_deg,
                ..base_options
            })
            .raw_unknowns,
        );
    }
    if maximum_unknowns > MCP_MAX_RECURRENCE_UNKNOWNS {
        return Err(invalid_params(format!(
            "recurrence search may create {maximum_unknowns} unknowns; the MCP limit is {MCP_MAX_RECURRENCE_UNKNOWNS}"
        )));
    }
    Ok(())
}

pub(crate) fn validate_recurrence_for_generation(
    recurrence: &Recurrence,
    initial_polys: &[Vec<BigRational>],
) -> Result<(), McpError> {
    if recurrence.terms.len() > MCP_MAX_RECURRENCE_LENGTH {
        return Err(invalid_params(format!(
            "recurrence has {} terms; the MCP limit is {MCP_MAX_RECURRENCE_LENGTH}",
            recurrence.terms.len()
        )));
    }
    let mut total_cells = initial_polys
        .iter()
        .try_fold(0usize, |total, row| total.checked_add(row.len()))
        .ok_or_else(|| invalid_params("recurrence coefficient count overflow"))?;
    if initial_polys.len() > MCP_MAX_RECURRENCE_INPUT_ROWS {
        return Err(invalid_params(format!(
            "recurrence JSON has {} initial rows; the MCP limit is {MCP_MAX_RECURRENCE_INPUT_ROWS}",
            initial_polys.len()
        )));
    }
    if initial_polys
        .iter()
        .any(|row| row.len() > MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL)
    {
        return Err(invalid_params(format!(
            "an initial polynomial exceeds the MCP limit of {MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL} coefficients"
        )));
    }

    let mut validate_bivar = |name: &str, polynomial: &BivarPoly| -> Result<(), McpError> {
        if polynomial.coeffs.len() > MCP_MAX_RECURRENCE_DEGREE + 1
            || polynomial
                .coeffs
                .iter()
                .any(|row| row.len() > MCP_MAX_RECURRENCE_DEGREE + 1)
        {
            return Err(invalid_params(format!(
                "{name} exceeds the MCP bivariate degree limit {MCP_MAX_RECURRENCE_DEGREE}"
            )));
        }
        for row in &polynomial.coeffs {
            total_cells = total_cells
                .checked_add(row.len())
                .ok_or_else(|| invalid_params("recurrence coefficient count overflow"))?;
        }
        Ok(())
    };
    for (index, term) in recurrence.terms.iter().enumerate() {
        if term.offset > MCP_MAX_RECURRENCE_LENGTH {
            return Err(invalid_params(format!(
                "recurrence term {index} has offset {}; the MCP limit is {MCP_MAX_RECURRENCE_LENGTH}",
                term.offset
            )));
        }
        if term.deriv_order > MCP_MAX_RECURRENCE_DEGREE {
            return Err(invalid_params(format!(
                "recurrence term {index} has derivative order {}; the MCP limit is {MCP_MAX_RECURRENCE_DEGREE}",
                term.deriv_order
            )));
        }
        validate_bivar(&format!("recurrence term {index}"), &term.coeff)?;
    }
    if let Some(denominator) = &recurrence.denominator {
        validate_bivar("recurrence denominator", denominator)?;
    }
    if let Some(inhomogeneous) = &recurrence.inhomogeneous {
        validate_bivar("recurrence inhomogeneous term", inhomogeneous)?;
    }
    if total_cells > MCP_MAX_TOTAL_COEFFICIENTS {
        return Err(invalid_params(format!(
            "recurrence JSON has {total_cells} stored coefficients; the MCP limit is {MCP_MAX_TOTAL_COEFFICIENTS}"
        )));
    }
    Ok(())
}

pub(crate) fn apply_recurrence_options(
    input: Option<RecurrenceSearchOptionsInput>,
) -> Result<(AdaptiveSearchOptions, AdaptiveSearchBudget), McpError> {
    let mut options = AdaptiveSearchOptions::default();
    let Some(input) = input else {
        options.verbose = false;
        let budget = AdaptiveSearchBudget::unbounded();
        validate_recurrence_options(&options, budget)?;
        return Ok((options, budget));
    };
    let budget = AdaptiveSearchBudget {
        max_candidates: input.max_candidates,
    };

    if let Some(value) = input.skip_prefix {
        options.skip_prefix = value;
    }
    if let Some(value) = input.min_rec_len {
        options.min_rec_len = value;
    }
    if let Some(value) = input.max_rec_len {
        options.max_rec_len = value;
    }
    if let Some(value) = input.min_var_deg {
        options.min_var_deg = value;
    }
    if let Some(value) = input.max_var_deg {
        options.max_var_deg = value;
    }
    if let Some(value) = input.min_idx_deg {
        options.min_idx_deg = value;
    }
    if let Some(value) = input.max_idx_deg {
        options.max_idx_deg = value;
    }
    if let Some(value) = input.min_diff_deg {
        options.min_diff_deg = value;
    }
    if let Some(value) = input.max_diff_deg {
        options.max_diff_deg = value;
    }
    if let Some(value) = input.try_inhomogeneous {
        options.try_inhomogeneous = value;
    }
    if let Some(value) = input.min_inhomo_var_deg {
        options.try_inhomogeneous = true;
        options.min_inhomo_var_deg = value;
    }
    if let Some(value) = input.max_inhomo_var_deg {
        options.try_inhomogeneous = true;
        options.max_inhomo_var_deg = value;
    }
    if let Some(value) = input.min_inhomo_idx_deg {
        options.try_inhomogeneous = true;
        options.min_inhomo_idx_deg = value;
    }
    if let Some(value) = input.max_inhomo_idx_deg {
        options.try_inhomogeneous = true;
        options.max_inhomo_idx_deg = value;
    }
    if let Some(value) = input.try_denominator {
        options.try_denominator = value;
    }
    if let Some(value) = input.try_alternating_sign {
        options.try_alternating_sign = value;
    }
    if let Some(value) = input.max_denom_var_deg {
        options.try_denominator = true;
        options.max_denom_var_deg = value;
    }
    if let Some(value) = input.max_denom_idx_deg {
        options.try_denominator = true;
        options.max_denom_idx_deg = value;
    }
    if let Some(value) = input.min_margin {
        options.min_margin = value;
    }
    if let Some(value) = input.no_verify {
        options.no_verify = value;
    }
    if let Some(value) = input.fit_extra_rows {
        options.fit_extra_rows = value;
    }
    if let Some(value) = input.modular_prefilter {
        options.modular_prefilter = value;
    }
    options.verbose = false;
    validate_recurrence_options(&options, budget)?;
    Ok((options, budget))
}
