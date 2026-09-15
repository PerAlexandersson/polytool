use num_bigint::BigInt;
use polytool::recurrence::BigRational as RecurrenceBigRational;
use polytool::recurrence::*;
use polytool::*;
use serde::Serialize;
use wasm_bindgen::prelude::*;

// ---------------------------------------------------------------------------
// Result types (serialized to JSON)
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct PolyProps {
    polynomial: String,
    degree: usize,
    coefficients: Vec<String>,
    real_rooted: bool,
    palindromic: bool,
    gamma_positive: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    gamma_coefficients: Option<Vec<String>>,
    log_concave: bool,
    ultra_log_concave: bool,
}

#[derive(Serialize)]
struct InterlacingResult {
    p: String,
    q: String,
    strict: bool,
    weak: bool,
    status: String,
}

#[derive(Serialize)]
struct RecurrenceResult {
    found: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    recurrence: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    latex: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    mathematica: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    sage: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    recurrence_json: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    unknowns: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    weighted_unknowns: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    equations: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    fit_polynomials: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    verification_polynomials: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    candidates_tried: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    error: Option<String>,
}

#[derive(Serialize)]
#[serde(untagged)]
enum ParsedOrError {
    Ok {
        polynomial: String,
        coefficients: Vec<String>,
    },
    Err {
        error: String,
    },
}

#[derive(Serialize)]
struct DisplayedPoly {
    polynomial: String,
    coefficients: Vec<String>,
}

#[derive(Serialize)]
struct MagicBasisReport {
    coordinates: Vec<String>,
    left_partial_sums: Vec<String>,
    right_partial_sums: Vec<String>,
    partial_sum_checks: Vec<bool>,
    all_nonnegative: bool,
    left_leq_right: bool,
}

#[derive(Serialize)]
struct DecompositionResult {
    polynomial: String,
    degree: usize,
    coefficients: Vec<String>,
    reciprocal: DisplayedPoly,
    a: DisplayedPoly,
    b: DisplayedPoly,
    a_real_rooted: bool,
    b_real_rooted: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    b_interlaces_a: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    reciprocal_interlaces_input: Option<bool>,
    alternatingly_increasing: bool,
    f_polynomial: DisplayedPoly,
    r_transform_of_f: DisplayedPoly,
    r_a: DisplayedPoly,
    r_b: DisplayedPoly,
    #[serde(skip_serializing_if = "Option::is_none")]
    r_interlaces_f: Option<bool>,
    magic: MagicBasisReport,
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn strip_trailing_zeros(coeffs: &[BigInt]) -> &[BigInt] {
    let zero = BigInt::from(0);
    let end = coeffs.iter().rposition(|c| c != &zero).map_or(0, |i| i + 1);
    if end == 0 {
        &[]
    } else {
        &coeffs[..end]
    }
}

fn decimal_coefficients(coeffs: &[BigInt]) -> Vec<String> {
    coeffs.iter().map(ToString::to_string).collect()
}

fn integer_polys_to_rational(polys: &[Vec<BigInt>]) -> Vec<Vec<RecurrenceBigRational>> {
    polys
        .iter()
        .map(|row| {
            row.iter()
                .cloned()
                .map(RecurrenceBigRational::from_integer)
                .collect()
        })
        .collect()
}

fn recurrence_json_output(
    result: &AdaptiveSearchResult,
    polys: &[Vec<BigInt>],
    search: &AdaptiveSearchOptions,
) -> String {
    let rational_polys = integer_polys_to_rational(polys);
    let searched_polys = rational_polys.get(search.skip_prefix..).unwrap_or(&[]);
    let initial_count = result.recurrence.max_offset().min(searched_polys.len());
    let recurrence_json = RecurrenceJson::from_recurrence_rational(
        &result.recurrence,
        1,
        &searched_polys[..initial_count],
        Some(RecurrenceJsonSearch {
            recurrence_text: result.recurrence.to_string(),
            source_rows: polys.len(),
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

fn parse_input(input: &str) -> Vec<Result<Vec<BigInt>, String>> {
    parse_polynomials_bigint(input)
}

fn parse_ok(input: &str) -> (Vec<Vec<BigInt>>, Vec<String>) {
    let mut polys = Vec::new();
    let mut errors = Vec::new();
    for r in parse_input(input) {
        match r {
            Ok(p) => polys.push(p),
            Err(e) => errors.push(e),
        }
    }
    (polys, errors)
}

fn display_poly(coeffs: &[BigInt]) -> DisplayedPoly {
    DisplayedPoly {
        polynomial: format_poly_bigint_coeffs(coeffs),
        coefficients: decimal_coefficients(coeffs),
    }
}

fn parse_coefficient_array_json(input: &str) -> Result<Vec<BigInt>, String> {
    let value: serde_json::Value = serde_json::from_str(input)
        .map_err(|error| format!("invalid coefficient JSON: {error}"))?;
    let values = value
        .as_array()
        .ok_or_else(|| "coefficient JSON must be an array".to_string())?;

    values
        .iter()
        .enumerate()
        .map(|(index, value)| {
            let decimal = match value {
                serde_json::Value::String(value) => value.clone(),
                serde_json::Value::Number(value) => value.to_string(),
                _ => {
                    return Err(format!(
                        "coefficient {index} must be a decimal string or JSON integer"
                    ));
                }
            };
            decimal
                .parse::<BigInt>()
                .map_err(|_| format!("coefficient {index} is not an integer: {decimal}"))
        })
        .collect()
}

// ---------------------------------------------------------------------------
// WASM exports
// ---------------------------------------------------------------------------

#[wasm_bindgen]
pub fn check_properties(input: &str) -> String {
    let mut results: Vec<serde_json::Value> = Vec::new();

    for r in parse_input(input) {
        match r {
            Ok(coeffs) => {
                let c = strip_trailing_zeros(&coeffs);
                let zero = BigInt::from(0);
                let deg = c.iter().rposition(|x| x != &zero).unwrap_or(0);
                let props = PolyProps {
                    polynomial: format_poly_bigint_coeffs(c),
                    degree: deg,
                    coefficients: decimal_coefficients(c),
                    real_rooted: is_real_rooted_bigint_coeffs(c),
                    palindromic: is_palindromic_ignoring_initial_zeros_bigint_coeffs(c),
                    gamma_positive: is_gamma_positive_ignoring_initial_zeros_bigint_coeffs(c),
                    gamma_coefficients: gamma_coefficients_ignoring_initial_zeros_bigint_coeffs(c)
                        .map(|coefficients| decimal_coefficients(&coefficients)),
                    log_concave: is_log_concave_bigint_coeffs(c),
                    ultra_log_concave: is_ultra_log_concave_bigint_coeffs(c),
                };
                results.push(serde_json::to_value(&props).unwrap());
            }
            Err(e) => {
                results.push(serde_json::json!({"error": e}));
            }
        }
    }

    serde_json::to_string(&results).unwrap()
}

#[wasm_bindgen]
pub fn check_interlacing_pairs(input: &str) -> String {
    let (polys, _) = parse_ok(input);
    let mut results: Vec<InterlacingResult> = Vec::new();

    for pair in polys.windows(2) {
        let p = strip_trailing_zeros(&pair[0]);
        let q = strip_trailing_zeros(&pair[1]);
        let strict = check_interlacing_bigint_coeffs(p, q) == Some(true);
        let weak = check_weak_interlacing_bigint_coeffs(p, q) == Some(true);
        let status = if strict {
            "strictly interlace".to_string()
        } else if weak {
            "weakly interlace (shared roots)".to_string()
        } else {
            "do NOT interlace".to_string()
        };
        results.push(InterlacingResult {
            p: format_poly_bigint_coeffs(p),
            q: format_poly_bigint_coeffs(q),
            strict,
            weak,
            status,
        });
    }

    serde_json::to_string(&results).unwrap()
}

#[wasm_bindgen]
pub fn compute_resultant(input: &str) -> String {
    let (polys, _) = parse_ok(input);
    if polys.len() < 2 {
        return serde_json::json!({"error": "need exactly two polynomials"}).to_string();
    }
    let r = resultant_bigint_coeffs(&polys[0], &polys[1]);
    serde_json::json!({
        "p": format_poly_bigint_coeffs(&polys[0]),
        "q": format_poly_bigint_coeffs(&polys[1]),
        "resultant": r.to_string()
    })
    .to_string()
}

#[wasm_bindgen]
pub fn compute_discriminant(input: &str) -> String {
    let mut results: Vec<serde_json::Value> = Vec::new();
    for r in parse_input(input) {
        match r {
            Ok(coeffs) => {
                let c = strip_trailing_zeros(&coeffs);
                let d = discriminant_bigint_coeffs(c);
                results.push(serde_json::json!({
                    "polynomial": format_poly_bigint_coeffs(c),
                    "discriminant": d.to_string()
                }));
            }
            Err(e) => {
                results.push(serde_json::json!({"error": e}));
            }
        }
    }
    serde_json::to_string(&results).unwrap()
}

#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn find_recurrence(
    input: &str,
    max_rec_len: u32,
    max_var_deg: u32,
    max_idx_deg: u32,
    max_diff_deg: u32,
    try_inhomogeneous: bool,
    try_denominator: bool,
    try_alternating_sign: bool,
) -> String {
    let (polys, _) = parse_ok(input);
    if polys.len() < 3 {
        return serde_json::to_string(&RecurrenceResult {
            found: false,
            recurrence: None,
            latex: None,
            mathematica: None,
            sage: None,
            recurrence_json: None,
            unknowns: None,
            weighted_unknowns: None,
            equations: None,
            fit_polynomials: None,
            verification_polynomials: None,
            candidates_tried: None,
            error: Some("need at least 3 polynomials".to_string()),
        })
        .unwrap();
    }

    let search = AdaptiveSearchOptions {
        max_rec_len: max_rec_len as usize,
        max_var_deg: max_var_deg as usize,
        max_idx_deg: max_idx_deg as usize,
        max_diff_deg: max_diff_deg as usize,
        try_inhomogeneous,
        try_denominator,
        try_alternating_sign,
        ..Default::default()
    };

    let rational_polys = integer_polys_to_rational(&polys);
    match find_recurrence_adaptive_rational(&rational_polys, &search) {
        Some(res) => serde_json::to_string(&RecurrenceResult {
            found: true,
            recurrence: Some(format!("{}", res.recurrence)),
            latex: Some(res.recurrence.to_latex()),
            mathematica: Some(
                res.recurrence
                    .to_mathematica_definition_rational(&rational_polys),
            ),
            sage: Some(res.recurrence.to_sage_definition_rational(&rational_polys)),
            recurrence_json: Some(recurrence_json_output(&res, &polys, &search)),
            unknowns: Some(res.num_unknowns),
            weighted_unknowns: Some(res.weighted_unknowns),
            equations: Some(res.num_equations),
            fit_polynomials: Some(res.fit_polynomials),
            verification_polynomials: Some(res.verification_polynomials),
            candidates_tried: Some(res.candidates_tried),
            error: None,
        })
        .unwrap(),
        None => serde_json::to_string(&RecurrenceResult {
            found: false,
            recurrence: None,
            latex: None,
            mathematica: None,
            sage: None,
            recurrence_json: None,
            unknowns: None,
            weighted_unknowns: None,
            equations: None,
            fit_polynomials: None,
            verification_polynomials: None,
            candidates_tried: None,
            error: Some("no recurrence found within search bounds".to_string()),
        })
        .unwrap(),
    }
}

#[wasm_bindgen]
pub fn analyze_decompositions(input: &str) -> String {
    let mut results: Vec<serde_json::Value> = Vec::new();

    for r in parse_input(input) {
        match r {
            Ok(coeffs) => {
                let c = strip_trailing_zeros(&coeffs);
                match analyze_symmetric_decomposition_bigint(c) {
                    Ok(analysis) => {
                        let report = DecompositionResult {
                            polynomial: format_poly_bigint_coeffs(c),
                            degree: analysis.degree,
                            coefficients: decimal_coefficients(c),
                            reciprocal: display_poly(&analysis.reciprocal),
                            a: display_poly(&analysis.a),
                            b: display_poly(&analysis.b),
                            a_real_rooted: analysis.a_real_rooted,
                            b_real_rooted: analysis.b_real_rooted,
                            b_interlaces_a: analysis.b_interlaces_a,
                            reciprocal_interlaces_input: analysis.reciprocal_interlaces_input,
                            alternatingly_increasing: analysis.alternatingly_increasing,
                            f_polynomial: display_poly(&analysis.f_polynomial),
                            r_transform_of_f: display_poly(&analysis.r_transform_of_f),
                            r_a: display_poly(&analysis.r_a),
                            r_b: display_poly(&analysis.r_b),
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
                                partial_sum_checks: analysis
                                    .magic
                                    .left_partial_sums
                                    .iter()
                                    .zip(analysis.magic.right_partial_sums.iter())
                                    .map(|(left, right)| left <= right)
                                    .collect(),
                                all_nonnegative: analysis.magic.all_nonnegative,
                                left_leq_right: analysis.magic.left_leq_right,
                            },
                        };
                        results.push(serde_json::to_value(&report).unwrap());
                    }
                    Err(e) => {
                        results.push(serde_json::json!({"error": e.to_string()}));
                    }
                }
            }
            Err(e) => {
                results.push(serde_json::json!({"error": e}));
            }
        }
    }

    serde_json::to_string(&results).unwrap()
}

/// Check interlacing between two polynomials given as JSON arrays of decimal
/// strings. Ordinary JSON integer arrays remain accepted for compatibility.
/// Returns JSON: {"strict": bool, "weak": bool}
#[wasm_bindgen]
pub fn check_interlacing_pair(p_json: &str, q_json: &str) -> String {
    let p = match parse_coefficient_array_json(p_json) {
        Ok(coefficients) => coefficients,
        Err(error) => return serde_json::json!({"error": error}).to_string(),
    };
    let q = match parse_coefficient_array_json(q_json) {
        Ok(coefficients) => coefficients,
        Err(error) => return serde_json::json!({"error": error}).to_string(),
    };
    let strict = check_interlacing_bigint_coeffs(&p, &q) == Some(true);
    let weak = if strict {
        true
    } else {
        check_weak_interlacing_bigint_coeffs(&p, &q) == Some(true)
    };
    serde_json::json!({"strict": strict, "weak": weak}).to_string()
}

#[wasm_bindgen]
pub fn parse_and_format(input: &str) -> String {
    let mut results: Vec<ParsedOrError> = Vec::new();
    for r in parse_input(input) {
        match r {
            Ok(coeffs) => {
                let c = strip_trailing_zeros(&coeffs);
                results.push(ParsedOrError::Ok {
                    polynomial: format_poly_bigint_coeffs(c),
                    coefficients: decimal_coefficients(c),
                });
            }
            Err(e) => {
                results.push(ParsedOrError::Err { error: e });
            }
        }
    }
    serde_json::to_string(&results).unwrap()
}

#[cfg(test)]
mod tests {
    use super::{
        analyze_decompositions, check_interlacing_pair, check_interlacing_pairs, check_properties,
        compute_discriminant, compute_resultant, find_recurrence, parse_and_format,
    };
    use num_bigint::BigInt;
    use serde_json::Value;

    const HUGE: &str = "1000000000000000000000000000000";

    fn scaled_huge(factor: i64) -> String {
        (HUGE.parse::<BigInt>().unwrap() * BigInt::from(factor)).to_string()
    }

    const EULERIAN_INPUT: &str = "\
1
1, 1
1, 4, 1
1, 11, 11, 1
1, 26, 66, 26, 1
1, 57, 302, 302, 57, 1
1, 120, 1191, 2416, 1191, 120, 1
1, 247, 4293, 15619, 15619, 4293, 247, 1
1, 502, 14608, 88234, 156190, 88234, 14608, 502, 1
1, 1013, 47840, 455192, 1310354, 1310354, 455192, 47840, 1013, 1
1, 2036, 152637, 2203488, 9738114, 15724248, 9738114, 2203488, 152637, 2036, 1
1, 4083, 478271, 10187685, 66318474, 162512286, 162512286, 66318474, 10187685, 478271, 4083, 1
1, 8178, 1479726, 45533450, 423281535, 1505621508, 2275172004, 1505621508, 423281535, 45533450, 1479726, 8178, 1
1, 16369, 4537314, 198410786, 2571742175, 12843262863, 27971176092, 27971176092, 12843262863, 2571742175, 198410786, 4537314, 16369, 1";

    const ALTERNATING_RUN_INPUT: &str = "\
2
2, 4
2, 12, 10
2, 28, 58, 32
2, 60, 236, 300, 122
2, 124, 836, 1852, 1682, 544
2, 252, 2766, 9576, 14622, 10332, 2770
2, 508, 8814, 45096, 103326, 119964, 69298, 15872
2, 1020, 27472, 201060, 650892, 1106820, 1034992, 505500, 101042";

    const DELANNOY_INPUT: &str = "\
1
1, 1
1, 3, 1
1, 5, 5, 1
1, 7, 13, 7, 1
1, 9, 25, 25, 9, 1
1, 11, 41, 63, 41, 11, 1
1, 13, 61, 129, 129, 61, 13, 1
1, 15, 85, 231, 321, 231, 85, 15, 1
1, 17, 113, 377, 681, 681, 377, 113, 17, 1
1, 19, 145, 575, 1289, 1683, 1289, 575, 145, 19, 1
1, 21, 181, 833, 2241, 3653, 3653, 2241, 833, 181, 21, 1";

    #[test]
    fn recurrence_export_handles_eulerian_example() {
        let raw = find_recurrence(EULERIAN_INPUT, 3, 2, 2, 1, false, false, false);
        let value: Value = serde_json::from_str(&raw).expect("recurrence result is valid JSON");
        assert_eq!(value["found"], true);
        assert_eq!(
            value["recurrence"],
            "P(n) = (1 - t + nt) P(n-1) + (t - t^2) P'(n-1)"
        );
        assert!(value.get("python").is_none());
        assert!(value["recurrence_json"].is_string());
    }

    #[test]
    fn recurrence_export_finds_alternating_runs_at_cubic_variable_degree() {
        let too_small = find_recurrence(ALTERNATING_RUN_INPUT, 3, 2, 2, 1, false, false, false);
        let too_small: Value =
            serde_json::from_str(&too_small).expect("recurrence result is valid JSON");
        assert_eq!(too_small["found"], false);

        let raw = find_recurrence(ALTERNATING_RUN_INPUT, 10, 5, 5, 5, false, false, false);
        let value: Value = serde_json::from_str(&raw).expect("recurrence result is valid JSON");
        assert_eq!(value["found"], true);
        assert!(value["recurrence"]
            .as_str()
            .expect("found recurrence has text")
            .contains("t^3"));
    }

    #[test]
    fn recurrence_export_finds_delannoy_example_quickly() {
        let raw = find_recurrence(DELANNOY_INPUT, 10, 5, 5, 5, false, false, false);
        let value: Value = serde_json::from_str(&raw).expect("recurrence result is valid JSON");
        assert_eq!(value["found"], true);
        assert_eq!(value["recurrence"], "P(n) = (1 + t) P(n-1) + t P(n-2)");
        assert!(value["candidates_tried"].as_u64().unwrap() <= 10);
    }

    #[test]
    fn ordinary_coefficients_are_returned_as_compatible_decimal_strings() {
        let parsed: Value = serde_json::from_str(&parse_and_format("1, 2, 1")).unwrap();
        assert_eq!(
            parsed[0]["coefficients"],
            serde_json::json!(["1", "2", "1"])
        );
        assert_eq!(parsed[0]["polynomial"], "1 + 2t + t^2");

        let interlacing: Value =
            serde_json::from_str(&check_interlacing_pair("[-2,1]", "[3,-4,1]")).unwrap();
        assert_eq!(interlacing["strict"], true);
        assert_eq!(interlacing["weak"], true);
    }

    #[test]
    fn property_and_interlacing_exports_accept_coefficients_above_i64() {
        let input = format!("{HUGE}, {}, {HUGE}", scaled_huge(2));
        let properties: Value = serde_json::from_str(&check_properties(&input)).unwrap();
        assert_eq!(
            properties[0]["coefficients"],
            serde_json::json!([HUGE, scaled_huge(2), HUGE])
        );
        assert_eq!(
            properties[0]["gamma_coefficients"],
            serde_json::json!([HUGE, "0"])
        );
        assert_eq!(properties[0]["real_rooted"], true);
        assert_eq!(properties[0]["palindromic"], true);
        assert_eq!(properties[0]["gamma_positive"], true);
        assert_eq!(properties[0]["log_concave"], true);
        assert_eq!(properties[0]["ultra_log_concave"], true);

        let interlacing_input = format!(
            "{}, {HUGE}\n{}, {}, {HUGE}",
            scaled_huge(-2),
            scaled_huge(3),
            scaled_huge(-4)
        );
        let pairs: Value =
            serde_json::from_str(&check_interlacing_pairs(&interlacing_input)).unwrap();
        assert_eq!(pairs[0]["strict"], true);
        assert_eq!(pairs[0]["weak"], true);

        let pair: Value = serde_json::from_str(&check_interlacing_pair(
            &serde_json::json!([scaled_huge(-2), HUGE]).to_string(),
            &serde_json::json!([scaled_huge(3), scaled_huge(-4), HUGE]).to_string(),
        ))
        .unwrap();
        assert_eq!(pair["strict"], true);
    }

    #[test]
    fn resultant_discriminant_and_decomposition_stay_exact_above_i64() {
        let resultant_input = format!("-{HUGE}, {HUGE}\n-2, 1");
        let resultant: Value = serde_json::from_str(&compute_resultant(&resultant_input)).unwrap();
        assert_eq!(resultant["resultant"], scaled_huge(-1));

        let discriminant_input = format!("{}, {}, {HUGE}", scaled_huge(2), scaled_huge(-3));
        let discriminant: Value =
            serde_json::from_str(&compute_discriminant(&discriminant_input)).unwrap();
        assert_eq!(
            discriminant[0]["discriminant"],
            (HUGE.parse::<BigInt>().unwrap().pow(2)).to_string()
        );

        let decomposition_input = format!("{HUGE}, {}, {}", scaled_huge(4), scaled_huge(2));
        let decomposition: Value =
            serde_json::from_str(&analyze_decompositions(&decomposition_input)).unwrap();
        assert_eq!(
            decomposition[0]["coefficients"],
            serde_json::json!([HUGE, scaled_huge(4), scaled_huge(2)])
        );
        assert_eq!(
            decomposition[0]["a"]["coefficients"],
            serde_json::json!([HUGE, scaled_huge(3), HUGE])
        );
        assert_eq!(
            decomposition[0]["magic"]["coordinates"],
            serde_json::json!([HUGE, scaled_huge(2), scaled_huge(-1)])
        );
    }

    #[test]
    fn recurrence_search_uses_big_rationals_for_huge_input() {
        let input = format!(
            "{HUGE}\n{}\n{}\n{}\n{}\n{}",
            scaled_huge(2),
            scaled_huge(4),
            scaled_huge(8),
            scaled_huge(16),
            scaled_huge(32)
        );
        let result: Value =
            serde_json::from_str(&find_recurrence(&input, 1, 0, 0, 0, false, false, false))
                .unwrap();
        assert_eq!(result["found"], true, "{result:#}");
        assert_eq!(result["recurrence"], "P(n) = 2 P(n-1)");

        let recurrence_json: Value =
            serde_json::from_str(result["recurrence_json"].as_str().unwrap()).unwrap();
        assert_eq!(recurrence_json["initial_polynomials"][0][0], HUGE);
    }
}
