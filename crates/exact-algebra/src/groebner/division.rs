//! Multivariate division and cached scalar/matrix normal forms.

use std::collections::{BTreeMap, BTreeSet, HashSet};

use super::{Coefficient, CoefficientRing, ExactField, SparsePolynomial};

use super::monomial_order::{leading_term, monomial_quotient, LeadingTerm, MonomialOrder};

/// Result of [`divide_by_polynomials`]: `dividend = sum_i quotients[i] * divisors[i] + remainder`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DivisionResult<P: SparsePolynomial> {
    /// One quotient per divisor, in the order the divisors were given.
    pub quotients: Vec<P>,
    /// Remainder; no term is divisible by any divisor's leading monomial.
    pub remainder: P,
}

/// Multiply `polynomial` by `coefficient * x^exponents`.
///
/// # Panics
///
/// Panics if `exponents.len()` differs from the number of variables or an
/// exponent overflows `u32`.
pub fn multiply_by_monomial<P: SparsePolynomial>(
    polynomial: &P,
    exponents: &[u32],
    coefficient: Coefficient<P>,
) -> P {
    assert_eq!(
        polynomial.num_vars(),
        exponents.len(),
        "monomial has wrong number of variables"
    );
    if polynomial.is_zero() || P::Field::is_zero(&coefficient) {
        return P::zero(polynomial.num_vars());
    }

    let terms = polynomial
        .terms()
        .iter()
        .map(|(term_exp, term_coeff)| {
            let new_exp = term_exp
                .iter()
                .zip(exponents.iter())
                .map(|(&a, &b)| {
                    a.checked_add(b)
                        .expect("monomial exponent overflow during multiplication")
                })
                .collect();
            (new_exp, term_coeff.clone() * coefficient.clone())
        })
        .collect();
    P::from_terms(polynomial.num_vars(), terms)
}

/// Divide `dividend` by `divisors` with respect to `order`.
///
/// Whenever the current leading term is divisible by some divisor's leading
/// term, the first such divisor (in slice order) is used; otherwise the term
/// moves to the remainder.
///
/// # Panics
///
/// Panics if a divisor is zero or the variable counts of `dividend` and the divisors differ.
pub fn divide_by_polynomials<P: SparsePolynomial<Field: ExactField>>(
    dividend: &P,
    divisors: &[P],
    order: MonomialOrder,
) -> DivisionResult<P> {
    assert!(
        divisors
            .iter()
            .all(|divisor| divisor.num_vars() == dividend.num_vars()),
        "all divisors must have the same number of variables as the dividend"
    );
    assert!(
        divisors.iter().all(|divisor| !divisor.is_zero()),
        "division by a zero polynomial"
    );

    let num_vars = dividend.num_vars();
    let mut quotients = vec![P::zero(num_vars); divisors.len()];
    let mut remainder = P::zero(num_vars);
    let mut current = dividend.clone();

    while let Some(current_lt) = leading_term(&current, order) {
        let mut reduced = false;

        for (i, divisor) in divisors.iter().enumerate() {
            let divisor_lt = leading_term(divisor, order).expect("nonzero divisor");
            let Some(exp_quotient) =
                monomial_quotient(&current_lt.exponents, &divisor_lt.exponents)
            else {
                continue;
            };
            let coeff_quotient =
                P::Field::divide(current_lt.coefficient.clone(), divisor_lt.coefficient);
            let quotient_term = P::monomial(num_vars, exp_quotient.clone(), coeff_quotient.clone());

            quotients[i] = quotients[i].clone() + quotient_term.clone();
            current = current - multiply_by_monomial(divisor, &exp_quotient, coeff_quotient);
            reduced = true;
            break;
        }

        if !reduced {
            let remainder_term =
                P::monomial(num_vars, current_lt.exponents, current_lt.coefficient);
            remainder = remainder + remainder_term.clone();
            current = current - remainder_term;
        }
    }

    DivisionResult {
        quotients,
        remainder,
    }
}

/// Remainder of [`divide_by_polynomials`]; a canonical normal form only when
/// `divisors` is a Gröbner basis for `order`.
///
/// # Panics
///
/// Panics under the conditions of [`divide_by_polynomials`].
pub fn normal_form<P: SparsePolynomial<Field: ExactField>>(
    polynomial: &P,
    divisors: &[P],
    order: MonomialOrder,
) -> P {
    divide_by_polynomials(polynomial, divisors, order).remainder
}

/// Same as [`normal_form`], but takes the divisors' leading terms precomputed.
///
/// `divisor_leading_terms[i]` must be `leading_term(&divisors[i], order)`;
/// this is not checked (only the lengths are).
///
/// # Panics
///
/// Panics if `divisors` and `divisor_leading_terms` have different lengths, or if a divisor does
/// not have the same number of variables as `polynomial`.
pub fn normal_form_with_leading_terms<P: SparsePolynomial<Field: ExactField>>(
    polynomial: &P,
    divisors: &[P],
    divisor_leading_terms: &[LeadingTerm<Coefficient<P>>],
    order: MonomialOrder,
) -> P {
    normal_form_observed(
        polynomial,
        divisors,
        divisor_leading_terms,
        order,
        &mut super::budget::Unlimited,
    )
    .unwrap_or_else(|never| match never {})
}

pub(super) fn normal_form_observed<
    P: SparsePolynomial<Field: ExactField>,
    O: super::budget::ReductionObserver<P>,
>(
    polynomial: &P,
    divisors: &[P],
    divisor_leading_terms: &[LeadingTerm<Coefficient<P>>],
    order: MonomialOrder,
    observer: &mut O,
) -> Result<P, O::Error> {
    assert_eq!(
        divisors.len(),
        divisor_leading_terms.len(),
        "divisor and leading-term lists have different lengths"
    );
    assert!(
        divisors
            .iter()
            .all(|divisor| divisor.num_vars() == polynomial.num_vars()),
        "all divisors must have the same number of variables as the polynomial"
    );

    let num_vars = polynomial.num_vars();
    let mut remainder = P::zero(num_vars);
    let mut current = polynomial.clone();

    while let Some(current_lt) = leading_term(&current, order) {
        let mut reduced = false;

        for (divisor, divisor_lt) in divisors.iter().zip(divisor_leading_terms.iter()) {
            let Some(exp_quotient) =
                monomial_quotient(&current_lt.exponents, &divisor_lt.exponents)
            else {
                continue;
            };
            let coeff_quotient = P::Field::divide(
                current_lt.coefficient.clone(),
                divisor_lt.coefficient.clone(),
            );
            current = current - multiply_by_monomial(divisor, &exp_quotient, coeff_quotient);
            observer.after_subtraction(&current)?;
            reduced = true;
            break;
        }

        if !reduced {
            let remainder_term =
                P::monomial(num_vars, current_lt.exponents, current_lt.coefficient);
            remainder = remainder + remainder_term.clone();
            current = current - remainder_term;
        }
    }

    Ok(remainder)
}

/// Compute several normal forms using one F4-style symbolic preprocessing and
/// matrix reduction step.
///
/// This constructs the reducer multiples forced by the current batch of
/// polynomials, stores them as rows with fixed symbolic pivots, and reduces
/// every input row through the same pivot list. For Buchberger computations
/// this gives a shared preprocessing step for a batch of S-polynomials, while
/// scalar reduction remains available for comparisons through
/// [`normal_form_with_leading_terms`].
///
/// # Panics
///
/// Panics if `divisors` and `divisor_leading_terms` have different lengths, or if the polynomials
/// and divisors do not all have the same number of variables.
pub fn matrix_normal_forms_with_leading_terms<P: SparsePolynomial<Field: ExactField>>(
    polynomials: &[P],
    divisors: &[P],
    divisor_leading_terms: &[LeadingTerm<Coefficient<P>>],
    order: MonomialOrder,
) -> Vec<P> {
    if polynomials.is_empty() {
        return Vec::new();
    }
    assert_eq!(
        divisors.len(),
        divisor_leading_terms.len(),
        "divisor and leading-term lists have different lengths"
    );
    let num_vars = polynomials[0].num_vars();
    assert!(
        polynomials
            .iter()
            .all(|polynomial| polynomial.num_vars() == num_vars),
        "all polynomials must have the same number of variables"
    );
    assert!(
        divisors
            .iter()
            .all(|divisor| divisor.num_vars() == num_vars),
        "all divisors must have the same number of variables as the polynomials"
    );
    if divisors.is_empty() {
        return polynomials.to_vec();
    }

    let (monomials, reducer_rows) =
        symbolic_reduction_rows(polynomials, divisors, divisor_leading_terms, order);
    if monomials.is_empty() || reducer_rows.is_empty() {
        return polynomials.to_vec();
    }

    let column_index: BTreeMap<_, _> = monomials
        .iter()
        .cloned()
        .enumerate()
        .map(|(i, monomial)| (monomial, i))
        .collect();
    let reducer_matrix: Vec<_> = reducer_rows
        .iter()
        .map(|(_, row)| polynomial_to_row(row, &column_index))
        .collect();
    let reducer_pivots = reducer_rows
        .iter()
        .map(|(pivot, _)| {
            *column_index
                .get(pivot)
                .expect("reducer pivot is outside the matrix columns")
        })
        .collect::<Vec<_>>();

    polynomials
        .iter()
        .map(|polynomial| {
            let mut row = polynomial_to_row(polynomial, &column_index);
            for (reducer_row, &pivot_col) in reducer_matrix.iter().zip(reducer_pivots.iter()) {
                let pivot_value = row[pivot_col].clone();
                if P::Field::is_zero(&pivot_value) {
                    continue;
                }
                for col in pivot_col..row.len() {
                    row[col] = row[col].clone() - pivot_value.clone() * reducer_row[col].clone();
                }
            }
            row_to_polynomial(num_vars, &monomials, &row)
        })
        .collect()
}

/// Monomials (in the chosen order) and the pending `(monomial, reducer)` rows.
type SymbolicReductionRows<P> = (Vec<Vec<u32>>, Vec<(Vec<u32>, P)>);

fn symbolic_reduction_rows<P: SparsePolynomial<Field: ExactField>>(
    polynomials: &[P],
    divisors: &[P],
    divisor_leading_terms: &[LeadingTerm<Coefficient<P>>],
    order: MonomialOrder,
) -> SymbolicReductionRows<P> {
    let mut monomial_set = BTreeSet::new();
    let mut pending = Vec::new();
    for polynomial in polynomials {
        for monomial in polynomial.terms().keys() {
            if monomial_set.insert(monomial.clone()) {
                pending.push(monomial.clone());
            }
        }
    }

    let mut reducer_pivots = HashSet::new();
    let mut reducer_rows = Vec::new();
    while let Some(monomial) = pending.pop() {
        let Some((divisor_index, exponent_quotient)) =
            first_dividing_leading_term(&monomial, divisor_leading_terms)
        else {
            continue;
        };
        if !reducer_pivots.insert(monomial.clone()) {
            continue;
        }
        let divisor_lt = &divisor_leading_terms[divisor_index];
        let row = multiply_by_monomial(
            &divisors[divisor_index],
            &exponent_quotient,
            P::Field::divide(P::Field::one(), divisor_lt.coefficient.clone()),
        );
        for new_monomial in row.terms().keys() {
            if monomial_set.insert(new_monomial.clone()) {
                pending.push(new_monomial.clone());
            }
        }
        reducer_rows.push((monomial, row));
    }

    let mut monomials: Vec<_> = monomial_set.into_iter().collect();
    monomials.sort_by(|a, b| order.compare(b, a));
    reducer_rows.sort_by(|(a, _), (b, _)| order.compare(b, a));
    (monomials, reducer_rows)
}

fn first_dividing_leading_term<C>(
    monomial: &[u32],
    divisor_leading_terms: &[LeadingTerm<C>],
) -> Option<(usize, Vec<u32>)> {
    divisor_leading_terms
        .iter()
        .enumerate()
        .find_map(|(i, divisor_lt)| {
            monomial_quotient(monomial, &divisor_lt.exponents).map(|quotient| (i, quotient))
        })
}

fn polynomial_to_row<P: SparsePolynomial>(
    polynomial: &P,
    column_index: &BTreeMap<Vec<u32>, usize>,
) -> Vec<Coefficient<P>> {
    let mut row = vec![P::Field::zero(); column_index.len()];
    for (monomial, coefficient) in polynomial.terms() {
        let &column = column_index
            .get(monomial)
            .expect("polynomial contains a monomial outside the matrix columns");
        row[column] = coefficient.clone();
    }
    row
}

fn row_to_polynomial<P: SparsePolynomial>(
    num_vars: usize,
    monomials: &[Vec<u32>],
    row: &[Coefficient<P>],
) -> P {
    let terms = monomials
        .iter()
        .cloned()
        .zip(row.iter().cloned())
        .filter_map(|(monomial, coefficient)| {
            (!P::Field::is_zero(&coefficient)).then_some((monomial, coefficient))
        })
        .collect();
    P::from_terms(num_vars, terms)
}
