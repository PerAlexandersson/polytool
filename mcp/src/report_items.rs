//! Conversion of library reports (properties, interlacing, coefficient tests,
//! h* inequalities, cyclic sieving, OEIS entries) into response items.

use crate::parsing::{
    bigint_strings, normalize_bigint_polynomial, BigIntParsedBatch, NormalizedBigIntPolynomial,
};
use crate::responses::{
    BigIntParsePolynomialItem, CoefficientCriterionItem, CoefficientInequalityCheckItem,
    CoefficientTestsItem, CyclicSievingPowerCheckItem, CyclicSievingReportItem,
    CyclicSievingSequenceItemResponse, DisplayedPolynomial, HStarInequalitiesItem,
    HStarInequalityCheckItem, InterlacingPairItem, InterlacingPairResult, OeisSequenceSummary,
    OrbitDecompositionItem, OrbitSizeEntryItem, OrbitSizesCheckItem, ParsePolynomialItem,
    PolynomialPropertiesItem, ResidueClassItem, RootOfUnityEvaluationItem,
};
use num_bigint::BigInt;
use polytool::*;

pub(crate) fn parse_items(batch: &BigIntParsedBatch) -> Vec<ParsePolynomialItem> {
    batch
        .iter()
        .enumerate()
        .map(|(index, item)| match item {
            Ok(polynomial) => ParsePolynomialItem {
                index,
                ok: true,
                polynomial: Some(polynomial.display.polynomial.clone()),
                coefficients: Some(polynomial.display.coefficients.clone()),
                degree: Some(polynomial.display.degree),
                error: None,
            },
            Err(error) => ParsePolynomialItem {
                index,
                ok: false,
                polynomial: None,
                coefficients: None,
                degree: None,
                error: Some(error.clone()),
            },
        })
        .collect()
}

pub(crate) fn parse_bigint_items(batch: &BigIntParsedBatch) -> Vec<BigIntParsePolynomialItem> {
    batch
        .iter()
        .enumerate()
        .map(|(index, item)| match item {
            Ok(polynomial) => BigIntParsePolynomialItem {
                index,
                ok: true,
                polynomial: Some(polynomial.display.clone()),
                error: None,
            },
            Err(error) => BigIntParsePolynomialItem {
                index,
                ok: false,
                polynomial: None,
                error: Some(error.clone()),
            },
        })
        .collect()
}

fn hstar_check_item(check: &HStarInequalityCheck) -> HStarInequalityCheckItem {
    HStarInequalityCheckItem {
        family: check.family.clone(),
        name: check.name.clone(),
        formula: check.formula.clone(),
        reference: check.reference.clone(),
        url: check.url.clone(),
        index: check.index,
        applicable: check.applicable,
        holds: check.holds,
        lhs: check.lhs.as_ref().map(ToString::to_string),
        rhs: check.rhs.as_ref().map(ToString::to_string),
        value: check.value.clone(),
        details: check.details.clone(),
    }
}

pub(crate) fn hstar_inequalities_item(
    index: usize,
    item: Result<NormalizedBigIntPolynomial, String>,
    dimension: Option<usize>,
) -> HStarInequalitiesItem {
    match item {
        Ok(polynomial) => {
            let report = hstar_inequality_report_bigint(&polynomial.coefficients, dimension);
            HStarInequalitiesItem {
                index,
                ok: true,
                error: None,
                hstar: Some(bigint_strings(&report.hstar)),
                dimension: Some(report.dimension),
                degree: Some(report.degree),
                codegree: report.codegree,
                all_applicable_hold: Some(report.all_applicable_hold),
                checks: report.checks.iter().map(hstar_check_item).collect(),
            }
        }
        Err(error) => HStarInequalitiesItem {
            index,
            ok: false,
            error: Some(error),
            hstar: None,
            dimension: None,
            degree: None,
            codegree: None,
            all_applicable_hold: None,
            checks: Vec::new(),
        },
    }
}

fn coefficient_check_item(check: &CoefficientInequalityCheck) -> CoefficientInequalityCheckItem {
    CoefficientInequalityCheckItem {
        index: check.index,
        lhs: check.lhs.to_string(),
        rhs: check.rhs.to_string(),
        comparison: check.comparison.to_string(),
        holds: check.holds,
    }
}

fn coefficient_criterion_item(report: &CoefficientCriterionReport) -> CoefficientCriterionItem {
    CoefficientCriterionItem {
        name: report.name.to_string(),
        reference: report.reference.to_string(),
        applicable: report.applicable,
        holds: report.holds,
        implies_real_rooted: report.implies_real_rooted,
        reason: report.reason.clone(),
        checks: report.checks.iter().map(coefficient_check_item).collect(),
    }
}

pub(crate) fn coefficient_tests_item(
    index: usize,
    item: Result<NormalizedBigIntPolynomial, String>,
) -> CoefficientTestsItem {
    match item {
        Ok(polynomial) => {
            let report = coefficient_test_report_bigint(&polynomial.coefficients);
            CoefficientTestsItem {
                index,
                ok: true,
                error: None,
                polynomial: Some(polynomial.display),
                newton: Some(coefficient_criterion_item(&report.newton)),
                kurtz: Some(coefficient_criterion_item(&report.kurtz)),
            }
        }
        Err(error) => CoefficientTestsItem {
            index,
            ok: false,
            error: Some(error),
            polynomial: None,
            newton: None,
            kurtz: None,
        },
    }
}

fn root_evaluation_item(evaluation: &RootOfUnityEvaluation) -> RootOfUnityEvaluationItem {
    RootOfUnityEvaluationItem {
        group_order: evaluation.group_order,
        power: evaluation.power,
        root_order: evaluation.root_order,
        integer_value: evaluation.integer_value.as_ref().map(ToString::to_string),
        remainder: bigint_strings(&evaluation.remainder),
        remainder_polynomial: format_poly_bigint_coeffs(&evaluation.remainder),
    }
}

fn cyclic_power_check_item(check: &CyclicSievingPowerCheck) -> CyclicSievingPowerCheckItem {
    CyclicSievingPowerCheckItem {
        power: check.power,
        root_order: check.root_order,
        expected_fixed_points: check.expected_fixed_points.to_string(),
        evaluation: root_evaluation_item(&check.evaluation),
        holds: check.holds,
    }
}

fn orbit_entry_items(entries: &[OrbitSizeEntry]) -> Vec<OrbitSizeEntryItem> {
    entries
        .iter()
        .map(|entry| OrbitSizeEntryItem {
            orbit_size: entry.orbit_size,
            basis_index: entry.basis_index,
            orbit_count: entry.orbit_count.to_string(),
            element_count: entry.element_count.to_string(),
        })
        .collect()
}

pub(crate) fn orbit_decomposition_item(dec: &OrbitDecomposition) -> OrbitDecompositionItem {
    OrbitDecompositionItem {
        order: dec.order,
        reduced_residues: bigint_strings(&dec.reduced_residues),
        classes: dec
            .classes
            .iter()
            .map(|class| ResidueClassItem {
                gcd_class: class.gcd_class,
                indices: class.indices.clone(),
                residues: bigint_strings(&class.residues),
                common_value: class.common_value.as_ref().map(ToString::to_string),
            })
            .collect(),
        integral: dec.integral,
        violating_classes: dec.violating_classes.clone(),
        orbits: dec.orbits.as_deref().map(orbit_entry_items),
        fixed_points: dec.fixed_points.as_deref().map(bigint_strings),
        natural_coefficients: dec.natural_coefficients,
        cone_condition: dec.cone_condition,
        realizable: dec.realizable,
    }
}

pub(crate) fn orbit_sizes_check_item(check: &OrbitSizesCheck) -> OrbitSizesCheckItem {
    OrbitSizesCheckItem {
        orbit_sizes: check.orbit_sizes.clone(),
        num_elements: check.num_elements.to_string(),
        actual_orbits: orbit_entry_items(&check.actual_orbits),
        actual_fixed_points: bigint_strings(&check.actual_fixed_points),
        mismatched_orbit_sizes: check.mismatched_orbit_sizes.clone(),
        holds: check.holds,
    }
}

pub(crate) fn bigint_matrix_strings(matrix: &[Vec<BigInt>]) -> Vec<Vec<String>> {
    matrix.iter().map(|row| bigint_strings(row)).collect()
}

pub(crate) fn cyclic_sieving_report_item(report: &CyclicSievingReport) -> CyclicSievingReportItem {
    CyclicSievingReportItem {
        order: report.order,
        polynomial: normalize_bigint_polynomial(report.coefficients.clone()).display,
        fixed_counts: report
            .fixed_counts
            .as_ref()
            .map(|counts| bigint_strings(counts)),
        holds: report.holds,
        evaluations: report
            .evaluations
            .iter()
            .map(root_evaluation_item)
            .collect(),
        checks: report.checks.iter().map(cyclic_power_check_item).collect(),
    }
}

pub(crate) fn cyclic_sequence_item_response(
    item: &CyclicSievingSequenceItem,
) -> CyclicSievingSequenceItemResponse {
    CyclicSievingSequenceItemResponse {
        row: item.row,
        index: item.index,
        polynomial: normalize_bigint_polynomial(item.coefficients.clone()).display,
        candidate_orders: item
            .candidate_orders
            .iter()
            .map(cyclic_sieving_report_item)
            .collect(),
    }
}

pub(crate) fn interlacing_result(
    p: NormalizedBigIntPolynomial,
    q: NormalizedBigIntPolynomial,
) -> InterlacingPairResult {
    let strict = check_interlacing_bigint_coeffs(&p.coefficients, &q.coefficients);
    let weak = check_weak_interlacing_bigint_coeffs(&p.coefficients, &q.coefficients);
    let status = match (strict, weak) {
        (Some(true), _) => "strictly_interlace",
        (_, Some(true)) => "weakly_interlace",
        (Some(false), Some(false)) => "do_not_interlace",
        (Some(false), None) => "not_real_rooted_or_incompatible",
        (None, Some(false)) => "not_real_rooted_or_incompatible",
        (None, None) => "not_real_rooted_or_incompatible",
    }
    .to_string();
    InterlacingPairResult {
        p: p.display,
        q: q.display,
        strict,
        weak,
        status,
    }
}

pub(crate) fn interlacing_pair_has_interlacing(item: &InterlacingPairItem) -> bool {
    item.result.as_ref().is_some_and(|result| {
        matches!(
            result.status.as_str(),
            "strictly_interlace" | "weakly_interlace"
        )
    })
}

pub(crate) fn interlacing_pair_item(
    pair_index: usize,
    left_index: usize,
    right_index: usize,
    left: &Result<NormalizedBigIntPolynomial, String>,
    right: &Result<NormalizedBigIntPolynomial, String>,
) -> InterlacingPairItem {
    match (left, right) {
        (Ok(p), Ok(q)) => InterlacingPairItem {
            pair_index,
            left_index,
            right_index,
            ok: true,
            result: Some(interlacing_result(p.clone(), q.clone())),
            error: None,
        },
        _ => {
            let mut messages = Vec::new();
            if let Err(error) = left {
                messages.push(format!("polynomial {left_index}: {error}"));
            }
            if let Err(error) = right {
                messages.push(format!("polynomial {right_index}: {error}"));
            }
            InterlacingPairItem {
                pair_index,
                left_index,
                right_index,
                ok: false,
                result: None,
                error: Some(messages.join("; ")),
            }
        }
    }
}

/// The exact property verdicts of one integer polynomial.
pub(crate) struct PolynomialPropertyValues {
    pub(crate) real_rooted: bool,
    pub(crate) simple_roots: bool,
    pub(crate) palindromic: bool,
    pub(crate) gamma_positive: bool,
    pub(crate) gamma_coefficients: Option<Vec<BigInt>>,
    pub(crate) unimodal: bool,
    pub(crate) log_concave: bool,
    pub(crate) ultra_log_concave: bool,
}

/// Compute the property verdicts for ascending `BigInt` coefficients.
///
/// Rows that fit in `i64` use the `i64` routines for the yes/no checks; the
/// gamma expansion is always computed with `BigInt` arithmetic, so a
/// palindromic row whose gamma coefficients exceed `i64` still reports them
/// (and its gamma-positivity) exactly.
pub(crate) fn polynomial_property_values(coefficients: &[BigInt]) -> PolynomialPropertyValues {
    let gamma_coefficients = gamma_coefficients_ignoring_initial_zeros_bigint_coeffs(coefficients);
    let zero = BigInt::from(0);
    let gamma_positive = gamma_coefficients
        .as_ref()
        .is_some_and(|gamma| gamma.iter().all(|value| value >= &zero));
    if let Some(small) = bigint_coeffs_to_i64(coefficients) {
        return PolynomialPropertyValues {
            real_rooted: is_real_rooted(&small),
            simple_roots: has_simple_roots(&small),
            palindromic: is_palindromic_ignoring_initial_zeros(&small),
            gamma_positive,
            gamma_coefficients,
            unimodal: is_unimodal(&small),
            log_concave: is_log_concave(&small),
            ultra_log_concave: is_ultra_log_concave(&small),
        };
    }
    PolynomialPropertyValues {
        real_rooted: is_real_rooted_bigint_coeffs(coefficients),
        simple_roots: has_simple_roots_bigint_coeffs(coefficients),
        palindromic: is_palindromic_ignoring_initial_zeros_bigint_coeffs(coefficients),
        gamma_positive,
        gamma_coefficients,
        unimodal: is_unimodal_bigint_coeffs(coefficients),
        log_concave: is_log_concave_bigint_coeffs(coefficients),
        ultra_log_concave: is_ultra_log_concave_bigint_coeffs(coefficients),
    }
}

pub(crate) fn polynomial_properties_item_bigint(
    index: usize,
    polynomial: NormalizedBigIntPolynomial,
) -> PolynomialPropertiesItem {
    let values = polynomial_property_values(&polynomial.coefficients);
    PolynomialPropertiesItem {
        index,
        ok: true,
        error: None,
        polynomial: Some(polynomial.display),
        real_rooted: Some(values.real_rooted),
        simple_roots: Some(values.simple_roots),
        palindromic: Some(values.palindromic),
        gamma_positive: Some(values.gamma_positive),
        gamma_coefficients: values
            .gamma_coefficients
            .map(|gamma| bigint_strings(&gamma)),
        unimodal: Some(values.unimodal),
        log_concave: Some(values.log_concave),
        ultra_log_concave: Some(values.ultra_log_concave),
    }
}

pub(crate) fn oeis_sequence_summary(
    entry: &polytool::oeis::OeisSequenceDefinition,
) -> OeisSequenceSummary {
    OeisSequenceSummary {
        id: entry.id.to_string(),
        name: entry.name.to_string(),
        status: entry.status.as_str().to_string(),
        layout: entry.layout.as_str().to_string(),
        first_row: entry.first_row,
        flattened_offset: entry.flattened_offset,
        bfile_available: entry.bfile_prefix_verified,
        source_rows: entry.source_rows,
        verification_rows: entry.verification_rows,
    }
}

pub(crate) fn displayed(coefficients: &[BigInt]) -> DisplayedPolynomial {
    DisplayedPolynomial {
        polynomial: format_poly_bigint_coeffs(coefficients),
        coefficients: bigint_strings(coefficients),
    }
}

pub(crate) fn format_rational<T: ToString>(rational: T) -> String {
    rational.to_string()
}
