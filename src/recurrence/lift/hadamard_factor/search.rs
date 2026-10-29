//! The bounded search of [`find_hadamard_factor_lifts`]: the factor menu,
//! exact quotient rows, the choice of undetermined quotient coefficients and
//! of the row normalisation, the parent search, and the per-candidate
//! certification.

use super::annihilator::{find_invariant_annihilator, separated_multiplier, step_operator};
use super::certify::{
    build_certificate, check_annihilated_rows, conjugate_input, divide_by_parent,
    right_pseudo_divide, ProductParent, RowAnnihilation, MAX_HADAMARD_OPERATOR_TERMS,
};
use super::factor::{HadamardFactor, HadamardFactorTerm};
use super::prepare::{prepare_rule, PreparedRule, MAX_LIFT_IDENTIFYING_ROWS};
use super::report::*;
use crate::linalg::{rational_kernel_with_certificate, RationalKernelLimits};
use crate::recurrence::implicit::algebra::{q, roots_in_domain, DOp};
use crate::recurrence::lift::parent_proof::{prove_parent, ParentRecord};
use crate::recurrence::lift::section_certify::n_coefficients;
use crate::recurrence::{
    find_recurrence_adaptive_rational_with_budget, format_rational_coeff, AdaptiveSearchBudget,
    AdaptiveSearchOptions, AdaptiveSearchOutcome, BivarPoly, Recurrence, RecurrenceError,
    RecurrenceJson, VeroneseParentAttestation, VeroneseSectionGauge,
};
use crate::vec_poly;
use crate::BigRational;
use num_traits::{One, Signed, Zero};

type Q = BigRational;

/// Input rows and candidate rows (`None` where the candidate is undefined).
type ComparedRows = (Vec<Vec<Q>>, Vec<Option<Vec<Q>>>);

/// Largest number of rows used for detection or compared exactly.
pub const MAX_HADAMARD_FACTOR_ROWS: usize = 64;

/// Largest `|a|`, `|b|` of the binomial factors `binom(n + a, k + b)`.
pub const MAX_HADAMARD_BINOMIAL_OFFSET: i64 = 4;

/// Steps by which the undetermined top or bottom quotient coefficients are
/// continued at most.
const CONTINUATION_STEPS: usize = 8;

/// Bounds of [`find_hadamard_factor_lifts`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HadamardFactorLiftOptions {
    /// Binomial factors `binom(n + a, k + b)` with `|a|, |b| <=` this bound
    /// (at most [`MAX_HADAMARD_BINOMIAL_OFFSET`]).  Default 2.
    pub max_binomial_offset: i64,
    /// Parameters `alpha` of the factors `(alpha)_k` and `1/(alpha)_k`
    /// (nonpositive integers are skipped).  Default `[1/2, 3/2]`.
    pub rising_factorial_parameters: Vec<BigRational>,
    /// Also try products of two menu factors.  Default false.
    pub include_products: bool,
    /// Input rows used for detection (quotients and the parent search),
    /// between 6 and [`MAX_HADAMARD_FACTOR_ROWS`].  Default 14.
    pub detection_rows: usize,
    /// Rows compared exactly with each candidate (at least the identifying
    /// rows), at most [`MAX_HADAMARD_FACTOR_ROWS`].  Default 16.
    pub checked_rows: usize,
    /// Number of menu factors tried; reaching it sets
    /// [`HadamardFactorLiftReport::exhausted_budget`].  Default 128.
    pub max_factors: usize,
    /// Candidate budget of each parent search
    /// ([`find_recurrence_adaptive_rational_with_budget`]).  Default 400.
    pub parent_search_candidates: usize,
    /// Largest derivative order of a parent rule.  Default 1.
    pub max_parent_derivative_order: usize,
    /// Operator-term budget of a certificate, at most
    /// [`MAX_HADAMARD_OPERATOR_TERMS`].  Default [`MAX_HADAMARD_OPERATOR_TERMS`].
    pub max_operator_terms: usize,
    /// Caller-supplied parent attestations (the attestation route of the
    /// parent provers; a lift is then conditional, never proved).
    pub attestations: Vec<VeroneseParentAttestation>,
    /// Try the root-domain parent route.  Default true.
    pub try_root_domain: bool,
}

impl Default for HadamardFactorLiftOptions {
    fn default() -> Self {
        Self {
            max_binomial_offset: 2,
            rising_factorial_parameters: vec![
                Q::new(1.into(), 2.into()),
                Q::new(3.into(), 2.into()),
            ],
            include_products: false,
            detection_rows: 14,
            checked_rows: 16,
            max_factors: 128,
            parent_search_candidates: 400,
            max_parent_derivative_order: 1,
            max_operator_terms: MAX_HADAMARD_OPERATOR_TERMS,
            attestations: Vec::new(),
            try_root_domain: true,
        }
    }
}

fn validate_options(options: &HadamardFactorLiftOptions) -> Result<(), RecurrenceError> {
    let bad = |what: String| Err(RecurrenceError::LimitExceeded(what));
    if !(0..=MAX_HADAMARD_BINOMIAL_OFFSET).contains(&options.max_binomial_offset) {
        return bad(format!(
            "max_binomial_offset must be between 0 and {MAX_HADAMARD_BINOMIAL_OFFSET}"
        ));
    }
    if !(6..=MAX_HADAMARD_FACTOR_ROWS).contains(&options.detection_rows) {
        return bad(format!(
            "detection_rows must be between 6 and {MAX_HADAMARD_FACTOR_ROWS}"
        ));
    }
    if !(1..=MAX_HADAMARD_FACTOR_ROWS).contains(&options.checked_rows) {
        return bad(format!(
            "checked_rows must be between 1 and {MAX_HADAMARD_FACTOR_ROWS}"
        ));
    }
    if options.max_operator_terms > MAX_HADAMARD_OPERATOR_TERMS {
        return bad(format!(
            "max_operator_terms must be at most {MAX_HADAMARD_OPERATOR_TERMS}"
        ));
    }
    if options.max_parent_derivative_order > 3 {
        return bad("max_parent_derivative_order must be at most 3".into());
    }
    Ok(())
}

/// The menu, singles first, then (optionally) unordered pairs.
fn factor_menu(options: &HadamardFactorLiftOptions) -> Vec<HadamardFactor> {
    let bound = options.max_binomial_offset;
    let mut singles = Vec::new();
    for a in -bound..=bound {
        for b in -bound..=bound {
            singles.push(HadamardFactorTerm::Binomial {
                index_offset: a,
                coefficient_offset: b,
            });
        }
    }
    singles.extend([
        HadamardFactorTerm::InverseFactorial,
        HadamardFactorTerm::Factorial,
        HadamardFactorTerm::InverseComplementFactorial,
        HadamardFactorTerm::ComplementFactorial,
    ]);
    for alpha in &options.rising_factorial_parameters {
        if alpha.is_integer() && !alpha.is_positive() {
            continue;
        }
        let parameter = format_rational_coeff(alpha);
        singles.push(HadamardFactorTerm::InverseRisingFactorial {
            parameter: parameter.clone(),
        });
        singles.push(HadamardFactorTerm::RisingFactorial { parameter });
    }
    let mut menu: Vec<HadamardFactor> = singles
        .iter()
        .cloned()
        .map(HadamardFactor::single)
        .collect();
    if options.include_products {
        for (i, first) in singles.iter().enumerate() {
            for second in &singles[i..] {
                menu.push(HadamardFactor {
                    terms: vec![first.clone(), second.clone()],
                });
            }
        }
    }
    menu
}

/// The size of a rule: order, largest degree in `x` and in `n` of its
/// coefficients and leader.
fn rule_size(rule: &Recurrence) -> (usize, usize, usize) {
    let mut x_degree = 0;
    let mut n_degree = 0;
    let mut visit = |p: &BivarPoly| {
        for (i, row) in p.coeffs.iter().enumerate() {
            if let Some(d) = row.iter().rposition(|c| !c.is_zero()) {
                x_degree = x_degree.max(d);
                n_degree = n_degree.max(i);
            }
        }
    };
    rule.terms.iter().for_each(|t| visit(&t.coeff));
    if let Some(denominator) = &rule.denominator {
        visit(denominator);
    }
    (rule.max_offset(), x_degree, n_degree)
}

fn is_smaller(parent: (usize, usize, usize), input: (usize, usize, usize)) -> bool {
    parent.0 < input.0
        || (parent.0 == input.0
            && parent.1 <= input.1
            && parent.2 <= input.2
            && (parent.1 < input.1 || parent.2 < input.2))
}

/// The quotient of one row: `Some(value)` at visible positions (`u != 0`),
/// `None` at hidden positions, with whether `u` is zero (not undefined)
/// there.
struct QuotientRow {
    values: Vec<Option<Q>>,
    hidden_zero: Vec<bool>,
}

enum QuotientOutcome {
    Rows(Vec<QuotientRow>),
    Rejected(String),
}

fn quotient_rows(factor: &HadamardFactor, rows: &[Vec<Q>], n0: usize) -> QuotientOutcome {
    let mut out = Vec::with_capacity(rows.len());
    for (offset, row) in rows.iter().enumerate() {
        let n = n0 + offset;
        let width = row.len() + CONTINUATION_STEPS + 1;
        let mut values = Vec::with_capacity(width);
        let mut hidden_zero = Vec::with_capacity(width);
        for k in 0..width {
            let coefficient = row.get(k).cloned().unwrap_or_else(Q::zero);
            match factor.value(n, k) {
                Some(u) if !u.is_zero() => {
                    values.push(Some(coefficient / u));
                    hidden_zero.push(false);
                }
                u => {
                    if !coefficient.is_zero() {
                        return QuotientOutcome::Rejected(format!(
                            "u({n}, {k}) is {} but q({n}, {k}) = {coefficient}",
                            if u.is_some() { "zero" } else { "undefined" }
                        ));
                    }
                    values.push(None);
                    hidden_zero.push(u.is_some());
                }
            }
        }
        out.push(QuotientRow {
            values,
            hidden_zero,
        });
    }
    QuotientOutcome::Rows(out)
}

fn zero_fill(rows: &[QuotientRow]) -> Vec<Vec<Q>> {
    rows.iter()
        .map(|row| {
            vec_poly::trim(
                &row.values
                    .iter()
                    .map(|v| v.clone().unwrap_or_else(Q::zero))
                    .collect::<Vec<_>>(),
            )
        })
        .collect()
}

/// Fit `p(n, k+1) (b0 + b1 n + b2 k) = p(n, k) (a0 + a1 n + a2 k)` on every
/// pair of adjacent visible positions; a one-dimensional kernel gives the
/// ratio, which continues each row into its hidden positions where `u = 0`.
fn continuation_fill(rows: &[QuotientRow], n0: usize) -> Option<Vec<Vec<Q>>> {
    let mut equations = Vec::new();
    for (offset, row) in rows.iter().enumerate() {
        let n = q((n0 + offset) as i64);
        for k in 0..row.values.len().saturating_sub(1) {
            let (Some(left), Some(right)) = (&row.values[k], &row.values[k + 1]) else {
                continue;
            };
            if left.is_zero() && right.is_zero() {
                continue;
            }
            let kk = q(k as i64);
            equations.push(vec![
                -left.clone(),
                -left * &n,
                -left * &kk,
                right.clone(),
                right * &n,
                right * &kk,
            ]);
        }
    }
    let kernel = rational_kernel_with_certificate(
        &equations,
        6,
        RationalKernelLimits {
            max_entries: 1 << 16,
            max_operations: 1 << 22,
            max_coefficient_bits: 4096,
        },
    )
    .ok()?;
    let [ratio] = kernel.basis.as_slice() else {
        return None;
    };
    let at = |c: &[Q], n: usize, k: i64| &c[0] + &c[1] * q(n as i64) + &c[2] * q(k);
    let (alpha, beta) = ratio.split_at(3);
    let mut filled = Vec::with_capacity(rows.len());
    for (offset, row) in rows.iter().enumerate() {
        let n = n0 + offset;
        let mut values: Vec<Q> = row
            .values
            .iter()
            .map(|v| v.clone().unwrap_or_else(Q::zero))
            .collect();
        let visible: Vec<usize> = (0..row.values.len())
            .filter(|k| row.values[*k].is_some())
            .collect();
        let (Some(&low), Some(&high)) = (visible.first(), visible.last()) else {
            filled.push(vec_poly::trim(&values));
            continue;
        };
        // Upwards.
        let mut k = high;
        while k + 1 < values.len() && row.hidden_zero[k + 1] && !values[k].is_zero() {
            let denominator = at(beta, n, k as i64);
            if denominator.is_zero() {
                return None;
            }
            values[k + 1] = &values[k] * at(alpha, n, k as i64) / denominator;
            k += 1;
        }
        if k + 1 == values.len() && !values[k].is_zero() {
            return None;
        }
        // Downwards.
        let mut k = low;
        while k > 0 && row.hidden_zero[k - 1] && !values[k].is_zero() {
            let numerator = at(alpha, n, k as i64 - 1);
            if numerator.is_zero() {
                break;
            }
            values[k - 1] = &values[k] * at(beta, n, k as i64 - 1) / numerator;
            k -= 1;
        }
        filled.push(vec_poly::trim(&values));
    }
    Some(filled)
}

/// `lambda(n)` = leading coefficient of row `n` in ratio form
/// `lambda(n) / lambda(n - 1) = a(n) / b(n)` with `deg a, deg b <= 2`, fitted
/// exactly to every row; `None` when the rows have a zero row, the ratio is
/// constant one, or no such ratio fits.
fn leading_gauge(rows: &[Vec<Q>], n0: usize) -> Option<(VeroneseSectionGauge, Vec<Vec<Q>>)> {
    let leads: Vec<Q> = rows
        .iter()
        .map(|row| vec_poly::trim(row).last().cloned().filter(|c| !c.is_zero()))
        .collect::<Option<_>>()?;
    if leads.iter().all(|lead| lead == &leads[0]) {
        return None;
    }
    for degree in 0..=2usize {
        let columns = 2 * (degree + 1);
        let equations: Vec<Vec<Q>> = (1..leads.len())
            .map(|offset| {
                let n = q((n0 + offset) as i64);
                let mut equation = Vec::with_capacity(columns);
                let mut power = Q::one();
                let mut powers = Vec::new();
                for _ in 0..=degree {
                    powers.push(power.clone());
                    power *= &n;
                }
                // a(n) lambda(n - 1) - b(n) lambda(n) = 0.
                equation.extend(powers.iter().map(|p| p * &leads[offset - 1]));
                equation.extend(powers.iter().map(|p| -(p * &leads[offset])));
                equation
            })
            .collect();
        if equations.len() < columns + 1 {
            return None;
        }
        let Ok(kernel) = rational_kernel_with_certificate(
            &equations,
            columns,
            RationalKernelLimits {
                max_entries: 1 << 12,
                max_operations: 1 << 20,
                max_coefficient_bits: 4096,
            },
        ) else {
            continue;
        };
        let [solution] = kernel.basis.as_slice() else {
            continue;
        };
        let (a, b) = solution.split_at(degree + 1);
        let (a, b) = (vec_poly::trim(a), vec_poly::trim(b));
        if a.iter().all(Zero::is_zero) || b.iter().all(Zero::is_zero) {
            continue;
        }
        let strings = |p: &[Q]| p.iter().map(format_rational_coeff).collect::<Vec<_>>();
        let gauge = VeroneseSectionGauge::Ratio {
            numerator: strings(&a),
            denominator: strings(&b),
            initial_value: format_rational_coeff(&leads[0]),
        };
        if gauge.ratio_form(n0).is_err() {
            continue;
        }
        let normalized = rows
            .iter()
            .zip(&leads)
            .map(|(row, lead)| vec_poly::scale(&vec_poly::trim(row), &lead.recip()))
            .collect();
        return Some((gauge, normalized));
    }
    None
}

enum ParentSearch {
    Found(PreparedRule),
    NotSmaller(String),
    NotFound,
    BudgetExhausted,
}

fn search_parent(
    rows: &[Vec<Q>],
    n0: usize,
    input_size: (usize, usize, usize),
    options: &HadamardFactorLiftOptions,
) -> Result<ParentSearch, RecurrenceError> {
    let search = AdaptiveSearchOptions {
        max_rec_len: input_size.0,
        max_var_deg: input_size.1,
        max_idx_deg: input_size.2,
        max_diff_deg: options.max_parent_derivative_order,
        max_denom_var_deg: 0,
        max_denom_idx_deg: input_size.2.max(1),
        ..AdaptiveSearchOptions::default()
    };
    let outcome = find_recurrence_adaptive_rational_with_budget(
        rows,
        &search,
        AdaptiveSearchBudget::limited(options.parent_search_candidates),
    );
    let found = match outcome {
        AdaptiveSearchOutcome::Found(found) => found,
        AdaptiveSearchOutcome::NoRecurrence(_) => return Ok(ParentSearch::NotFound),
        AdaptiveSearchOutcome::BudgetExhausted(_) => return Ok(ParentSearch::BudgetExhausted),
    };
    // The search indexes the first row as 1.
    let rule = found.recurrence.shift_index(1 - n0 as i64);
    let size = rule_size(&rule);
    if !is_smaller(size, input_size) {
        return Ok(ParentSearch::NotSmaller(format!(
            "the quotient rule has (order, x-degree, n-degree) = {size:?}, not smaller than the \
             input's {input_size:?}"
        )));
    }
    match prepare_rule(&rule, rows, n0, "parent") {
        Ok(parent) => Ok(ParentSearch::Found(parent)),
        Err(error) => Ok(ParentSearch::NotSmaller(format!(
            "the quotient rule is not usable: {error}"
        ))),
    }
}

/// Rows `lambda(n) sum_k u(n, k) p(n, k) x^k` for `n0 <= n < n0 + count`;
/// `None` for a row where `u` is undefined at a nonzero parent coefficient.
fn candidate_rows(
    factor: &HadamardFactor,
    lambdas: &[Q],
    parent_rows: &[Vec<Q>],
    n0: usize,
) -> Vec<Option<Vec<Q>>> {
    parent_rows
        .iter()
        .zip(lambdas)
        .enumerate()
        .map(|(offset, (row, lambda))| {
            let n = n0 + offset;
            let mut out = Vec::with_capacity(row.len());
            for (k, p) in row.iter().enumerate() {
                if p.is_zero() {
                    out.push(Q::zero());
                } else {
                    out.push(factor.value(n, k)? * p * lambda);
                }
            }
            Some(vec_poly::trim(&out))
        })
        .collect()
}

/// Remove a nonzero division remainder `r` with a row annihilator `A` of the
/// parent (`mu r = Z o A`): the product-parent annihilator, else an
/// invariant annihilator of a first-order parent whose rows are checked
/// from the base `from` on.  The inner `Err` is the reason when neither
/// applies (the identity is then only supported by rows).
fn remove_remainder(
    parent: &PreparedRule,
    remainder: &DOp,
    from: usize,
    options: &HadamardFactorLiftOptions,
) -> Result<Result<RowAnnihilation, String>, RecurrenceError> {
    let reason = "the division of the conjugated input relation by the parent relation leaves \
                  a nonzero remainder";
    let budget = |error: RecurrenceError| match error {
        RecurrenceError::BudgetExhausted(detail) | RecurrenceError::LimitExceeded(detail) => {
            Ok(Err(format!(
                "{reason}, and the row annihilator ran out of budget: {detail}"
            )))
        }
        error => Err(error),
    };
    if let Some(product) = ProductParent::recognise(parent) {
        let a = product.annihilator();
        if !product.annihilator_is_invariant(&a) {
            return Err(RecurrenceError::VerificationFailed(
                "the product-parent annihilator is not invariant".into(),
            ));
        }
        let (mu, z, rest) = match right_pseudo_divide(remainder, &a, options.max_operator_terms) {
            Ok(result) => result,
            Err(error) => return budget(error),
        };
        if !rest.is_zero() {
            return Ok(Err(format!(
                "{reason} that is not a left multiple of the row annihilator of the product \
                 parent"
            )));
        }
        return Ok(Ok(RowAnnihilation {
            row_multiplier: mu.rows().first().cloned().unwrap_or_default(),
            index_multiplier: vec![Q::one()],
            annihilator: a,
            quotient: z,
            invariance: None,
        }));
    }
    let Some(step) = step_operator(&parent.relation) else {
        return Ok(Err(format!(
            "{reason}, and the parent is neither a product parent nor of order one"
        )));
    };
    let invariant = match find_invariant_annihilator(remainder, &step, options.max_operator_terms) {
        Ok(Some(invariant)) => invariant,
        Ok(None) => {
            return Ok(Err(format!(
                "{reason}, and no invariant row annihilator of the parent divides it"
            )))
        }
        Err(error) => return budget(error),
    };
    let (mu, z, rest) = match right_pseudo_divide(
        remainder,
        &invariant.annihilator,
        options.max_operator_terms,
    ) {
        Ok(result) => result,
        Err(error) => return budget(error),
    };
    if !rest.is_zero() {
        return Ok(Err(format!(
            "{reason} that is not a left multiple of the invariant row annihilator"
        )));
    }
    let Some((index_multiplier, row_multiplier)) = separated_multiplier(&mu) else {
        return Ok(Err(format!(
            "{reason}; the right division by the invariant row annihilator needs a multiplier \
             that is not a product c(n) p(x)"
        )));
    };
    let invariance = RowAnnihilation::invariance_record(&invariant, from);
    match check_annihilated_rows(
        &invariant.annihilator,
        parent,
        from,
        &invariance.restart_bases,
    ) {
        Ok(()) => {}
        Err(RecurrenceError::VerificationFailed(detail)) => {
            return Ok(Err(format!("{reason}; {detail}")))
        }
        Err(error) => return budget(error),
    }
    Ok(Ok(RowAnnihilation {
        row_multiplier,
        index_multiplier,
        annihilator: invariant.annihilator,
        quotient: z,
        invariance: Some(invariance),
    }))
}

/// Certify one candidate (module documentation of `hadamard_factor`).
fn certify_candidate(
    input: &PreparedRule,
    parent: &PreparedRule,
    factor: &HadamardFactor,
    gauge: &VeroneseSectionGauge,
    options: &HadamardFactorLiftOptions,
) -> Result<(HadamardFactorIdentityStatus, usize), RecurrenceError> {
    let n0 = input.first_index;
    let ratio = gauge.ratio_form(n0)?;
    let compare = |count: usize| -> Result<ComparedRows, RecurrenceError> {
        let input_rows = input.rows(count)?;
        let parent_rows = parent.rows(count)?;
        let candidate = candidate_rows(factor, &ratio.values(n0, count), &parent_rows, n0);
        Ok((input_rows, candidate))
    };
    // An undefined candidate row counts as a difference.
    let mismatch = |rows: &ComparedRows| {
        rows.0
            .iter()
            .zip(&rows.1)
            .position(|(a, b)| b.as_ref() != Some(a))
    };
    let mut count = options.checked_rows.max(input.seeds.len());
    let mut rows = compare(count)?;
    if let Some(offset) = mismatch(&rows) {
        return Ok((
            HadamardFactorIdentityStatus::Refuted { index: n0 + offset },
            count,
        ));
    }
    let supported = |count: usize, reason: String| {
        Ok((
            HadamardFactorIdentityStatus::SupportedByRows {
                rows: count,
                reason,
            },
            count,
        ))
    };
    let conjugated = match conjugate_input(&input.relation, factor, &ratio, n0) {
        Ok(conjugated) => conjugated,
        Err(
            RecurrenceError::Unsupported(reason)
            | RecurrenceError::BudgetExhausted(reason)
            | RecurrenceError::LimitExceeded(reason),
        ) => return supported(count, format!("the conjugation does not apply: {reason}")),
        Err(error) => return Err(error),
    };
    let division =
        match divide_by_parent(&conjugated.op, &parent.relation, options.max_operator_terms) {
            Ok(division) => division,
            Err(
                RecurrenceError::BudgetExhausted(reason) | RecurrenceError::LimitExceeded(reason),
            ) => return supported(count, format!("the division ran out of budget: {reason}")),
            Err(error) => return Err(error),
        };
    let input_from = input.valid_from_base();
    let parent_from = parent.valid_from_base().max(n0);
    let from = input_from.max(parent_from);
    let mut annihilation = None;
    if !division.remainder_is_zero() {
        match remove_remainder(parent, &division.remainder[0], from, options)? {
            Ok(found) => annihilation = Some(found),
            Err(reason) => return supported(count, reason),
        }
    }
    let multiplier = match &annihilation {
        Some(found) => division
            .multiplier
            .mul(&super::super::section_certify::n_polynomial(
                &found.index_multiplier,
            )),
        None => division.multiplier.clone(),
    };

    // Exceptional bases.
    let mut exceptional: Vec<usize> = (input_from..parent_from).collect();
    exceptional.extend(roots_in_domain(&n_coefficients(&multiplier), from, None));
    exceptional.sort_unstable();
    exceptional.dedup();
    if let Some(last) = exceptional.last() {
        let needed = last + input.leader_shift + 1 - n0;
        if needed > MAX_HADAMARD_FACTOR_ROWS {
            return supported(
                count,
                format!("the exceptional base {last} needs {needed} compared rows"),
            );
        }
        if needed > count {
            count = needed;
            rows = compare(count)?;
            if let Some(offset) = mismatch(&rows) {
                return Ok((
                    HadamardFactorIdentityStatus::Refuted { index: n0 + offset },
                    count,
                ));
            }
        }
    }
    // The candidate rows equal the input rows on every compared row.
    for base in &exceptional {
        let value = input.relation.apply(&rows.0, n0, *base).ok_or_else(|| {
            RecurrenceError::VerificationFailed(format!("rows missing at exceptional base {base}"))
        })?;
        if !value.iter().all(Zero::is_zero) {
            return Err(RecurrenceError::VerificationFailed(format!(
                "the input relation does not vanish on the candidate rows at base {base}"
            )));
        }
    }
    let identifying = input.seeds.len();
    let statement = format!(
        "mu(x) w(n) M'' = sum_nu H_nu S^nu L{} with M'' the input relation conjugated by \
         lambda(n) {}; so the candidate rows satisfy the input rule at every base n >= \
         {input_from} (the exceptional bases checked exactly), and agreement on the \
         {identifying} identifying rows gives Q_n = lambda(n) sum_k u(n,k) [x^k]P_n x^k for \
         every n >= {n0}",
        match &annihilation {
            None => "",
            Some(found) if found.invariance.is_some() => {
                " + Z A (A an invariant row annihilator of the parent)"
            }
            Some(_) => " + Z A (A the row annihilator of the product parent)",
        },
        factor.text()
    );
    let certificate = build_certificate(
        factor,
        gauge,
        input,
        parent,
        &conjugated,
        &division,
        &multiplier,
        annihilation,
        (input_from, parent_from, exceptional),
        statement,
    );
    certificate.verify()?;
    Ok((
        HadamardFactorIdentityStatus::Certified {
            certificate: Box::new(certificate),
        },
        count,
    ))
}

/// Find Hadamard-factor lifts of the sequence defined by `recurrence` and
/// `rows` (`rows[0]` is `Q_(first_index)`, ascending coefficients): factors
/// `u(n, k)` of the menu with `Q_n = lambda(n) sum_k u(n, k) p(n, k) x^k`
/// for a parent `P` with a rule smaller than the input's.
///
/// The detection, the certificate and its soundness argument, the boundary
/// handling and the transfer theorem are documented in the
/// `hadamard_factor` module (see the [`crate::recurrence`] module list).
/// Every lift carries a [`crate::recurrence::LiftSummary`]; the verdict is
/// proved only when the identity is certified, the factor has a transfer
/// theorem and the parent provers prove the parent.  A missing lift says
/// nothing against real-rootedness of the input.
///
/// # Errors
///
/// [`RecurrenceError::InvalidInput`] for rows disagreeing with the rule;
/// [`RecurrenceError::Unsupported`] for a rule outside the class
/// (alternating or inhomogeneous terms, a leader not proved nonvanishing);
/// [`RecurrenceError::LimitExceeded`] for options beyond their bounds;
/// [`RecurrenceError::VerificationFailed`] if an internal self-check fails.
/// Budgets that run out are reported in the result, never as a panic.
///
/// # Example
///
/// ```
/// use polytool::oeis;
/// use polytool::recurrence::{find_hadamard_factor_lifts, HadamardFactorLiftOptions};
///
/// // Narayana triangle A001263 (rows in x^(k-1)).
/// let narayana = oeis::by_id("A001263").unwrap();
/// let (rule, rows) = narayana.recurrence_parts().unwrap();
/// let report = find_hadamard_factor_lifts(
///     &rule,
///     &rows,
///     narayana.recurrence_first_index,
///     &HadamardFactorLiftOptions { max_binomial_offset: 1, ..Default::default() },
/// )?;
/// let proved = report.proved_lifts().next().expect("a proved lift");
/// assert_eq!(proved.factor_text, "binom(n+1, k+1)");
/// report.verify()?;
/// # Ok::<(), polytool::recurrence::RecurrenceError>(())
/// ```
pub fn find_hadamard_factor_lifts(
    recurrence: &Recurrence,
    rows: &[Vec<BigRational>],
    first_index: usize,
    options: &HadamardFactorLiftOptions,
) -> Result<HadamardFactorLiftReport, RecurrenceError> {
    validate_options(options)?;
    let input = prepare_rule(recurrence, rows, first_index, "input")?;
    let n0 = first_index;
    let detection = options.detection_rows.max(input.seeds.len());
    if detection > MAX_LIFT_IDENTIFYING_ROWS {
        return Err(RecurrenceError::LimitExceeded(format!(
            "the input needs {detection} detection rows"
        )));
    }
    let input_rows = input.rows(detection)?;
    let input_size = rule_size(recurrence);
    let mut report = HadamardFactorLiftReport {
        schema: HADAMARD_FACTOR_LIFTS_SCHEMA.into(),
        input_first_index: n0,
        rows_examined: detection,
        factors_tried: 0,
        lifts: Vec::new(),
        rejected: Vec::new(),
        exhausted_budget: false,
        notes: Vec::new(),
    };
    for factor in factor_menu(options) {
        if report.factors_tried >= options.max_factors {
            report.exhausted_budget = true;
            report.notes.push(format!(
                "the factor budget of {} was reached",
                options.max_factors
            ));
            break;
        }
        report.factors_tried += 1;
        let reject = |kind, detail: String| HadamardFactorRejection {
            factor_text: factor.text(),
            factor: factor.clone(),
            kind,
            detail,
        };
        if let Err(error) = factor.validate(n0) {
            report.rejected.push(reject(
                HadamardFactorRejectionKind::InvalidForDomain,
                error.to_string(),
            ));
            continue;
        }
        let quotients = match quotient_rows(&factor, &input_rows, n0) {
            QuotientOutcome::Rows(rows) => rows,
            QuotientOutcome::Rejected(detail) => {
                report.rejected.push(reject(
                    HadamardFactorRejectionKind::VanishesOnSupport,
                    detail,
                ));
                continue;
            }
        };
        let mut fills = Vec::new();
        let zero = zero_fill(&quotients);
        if let Some(continued) = continuation_fill(&quotients, n0) {
            if continued != zero {
                fills.push((HadamardHiddenCoefficients::Continuation, continued));
            }
        }
        fills.push((HadamardHiddenCoefficients::Zero, zero));
        let mut found_any = false;
        let mut budget_hit = false;
        let mut reasons = Vec::new();
        'variants: for (fill, quotient) in &fills {
            let mut gauges = vec![(VeroneseSectionGauge::None, quotient.clone())];
            if let Some(normalized) = leading_gauge(quotient, n0) {
                gauges.push(normalized);
            }
            for (gauge, parent_rows) in gauges {
                let parent = match search_parent(&parent_rows, n0, input_size, options)? {
                    ParentSearch::Found(parent) => parent,
                    ParentSearch::NotSmaller(reason) => {
                        reasons.push(reason);
                        continue;
                    }
                    ParentSearch::NotFound => continue,
                    ParentSearch::BudgetExhausted => {
                        budget_hit = true;
                        continue;
                    }
                };
                let lift = build_lift(&input, &parent, &factor, &gauge, *fill, options)?;
                found_any = true;
                let proved = lift.lift_summary.is_proved();
                report.lifts.push(lift);
                if proved {
                    break 'variants;
                }
            }
        }
        if !found_any {
            let (kind, detail) = if budget_hit {
                (
                    HadamardFactorRejectionKind::ParentSearchBudgetExhausted,
                    format!(
                        "the parent search stopped after {} candidates",
                        options.parent_search_candidates
                    ),
                )
            } else {
                (
                    HadamardFactorRejectionKind::NoSmallerParent,
                    if reasons.is_empty() {
                        "no recurrence for the quotient rows in the searched space".into()
                    } else {
                        reasons.join("; ")
                    },
                )
            };
            report.rejected.push(reject(kind, detail));
        }
    }
    Ok(report)
}

fn build_lift(
    input: &PreparedRule,
    parent: &PreparedRule,
    factor: &HadamardFactor,
    gauge: &VeroneseSectionGauge,
    fill: HadamardHiddenCoefficients,
    options: &HadamardFactorLiftOptions,
) -> Result<HadamardFactorLift, RecurrenceError> {
    let (identity, rows_compared) = certify_candidate(input, parent, factor, gauge, options)?;
    let parent_proof = prove_parent(
        &ParentRecord {
            rule: &parent.rule,
            seeds: &parent.seeds,
            first_index: parent.first_index,
        },
        parent.first_index,
        &options.attestations,
        options.try_root_domain,
    )?;
    let parent_json = RecurrenceJson::from_recurrence_rational(
        &parent.rule,
        parent.first_index,
        &parent.seeds,
        None,
    );
    let multiplier_sequence = factor.multiplier_justification();
    let lift_summary = hadamard_lift_summary(
        &factor.text(),
        gauge,
        &parent_json,
        &identity,
        multiplier_sequence.as_deref(),
        &parent_proof,
    );
    Ok(HadamardFactorLift {
        factor: factor.clone(),
        factor_text: factor.text(),
        gauge: gauge.clone(),
        hidden_coefficients: fill,
        parent: parent_json,
        rows_compared,
        identity,
        multiplier_sequence,
        parent_proof,
        lift_summary,
    })
}
