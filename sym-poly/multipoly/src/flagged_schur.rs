//! Flagged Schur and flagged skew Schur polynomials by tableau enumeration.
//!
//! Rows are indexed in the usual partition convention, from top to bottom.
//! A row flag bounds entries in the corresponding row.

use std::collections::BTreeMap;

use sym_poly_core::{Partition, Ring, SkewTableau, Tableau};

use crate::multipoly::MultiPoly;

/// Enumerate semistandard tableaux of shape `lambda` with row intervals.
///
/// The entries in row `r` must lie in
/// `{lower_flags[r], ..., upper_flags[r]}`. Rows are weakly increasing and
/// columns are strictly increasing.
pub fn row_interval_flagged_tableaux(
    lambda: &[u32],
    lower_flags: &[u32],
    upper_flags: &[u32],
) -> Vec<Tableau> {
    validate_shape(lambda);
    validate_row_intervals(lambda.len(), lower_flags, upper_flags, None);

    if lambda.is_empty() {
        return vec![Tableau::empty()];
    }

    let mut filling = lambda
        .iter()
        .map(|&row_len| vec![0u32; row_len as usize])
        .collect::<Vec<_>>();
    let mut out = Vec::new();
    enumerate_tableaux_cell(
        lambda,
        lower_flags,
        upper_flags,
        &mut filling,
        0,
        0,
        &mut out,
    );
    out.into_iter().map(Tableau::new).collect()
}

/// Enumerate semistandard tableaux of shape `lambda` with upper row flags.
pub fn flagged_tableaux(lambda: &[u32], flags: &[u32]) -> Vec<Tableau> {
    let lower = vec![1u32; flags.len()];
    row_interval_flagged_tableaux(lambda, &lower, flags)
}

/// Flagged Schur polynomial with upper row flags.
pub fn flagged_schur<C: Ring>(lambda: &[u32], flags: &[u32], num_vars: usize) -> MultiPoly<C> {
    let lower = vec![1u32; flags.len()];
    row_interval_flagged_schur(lambda, &lower, flags, num_vars)
}

/// Flagged Schur polynomial with row intervals.
pub fn row_interval_flagged_schur<C: Ring>(
    lambda: &[u32],
    lower_flags: &[u32],
    upper_flags: &[u32],
    num_vars: usize,
) -> MultiPoly<C> {
    validate_shape(lambda);
    validate_row_intervals(lambda.len(), lower_flags, upper_flags, None);

    let mut filling = lambda
        .iter()
        .map(|&row_len| vec![0u32; row_len as usize])
        .collect::<Vec<_>>();
    let mut exponents = vec![0u32; num_vars];
    let mut terms = BTreeMap::new();
    accumulate_tableau_weights_cell(
        lambda,
        lower_flags,
        upper_flags,
        &mut filling,
        0,
        0,
        &mut exponents,
        &mut terms,
    );
    MultiPoly::from_terms(num_vars, terms)
}

/// Enumerate semistandard skew tableaux of shape `lambda / mu` with row
/// intervals.
///
/// The entries in row `r` must lie in
/// `{lower_flags[r], ..., upper_flags[r]}`. Rows are weakly increasing and
/// columns are strictly increasing in the visible skew cells.
pub fn row_interval_flagged_skew_tableaux(
    lambda: &[u32],
    mu: &[u32],
    lower_flags: &[u32],
    upper_flags: &[u32],
) -> Vec<SkewTableau> {
    validate_skew_shape(lambda, mu);
    validate_row_intervals(lambda.len(), lower_flags, upper_flags, None);

    if lambda.is_empty() {
        return vec![SkewTableau::empty()];
    }

    let max_width = lambda.iter().copied().max().unwrap_or(0) as usize;
    let mut filling = vec![vec![0u32; max_width]; lambda.len()];
    let mut cells = Vec::new();
    for row in 0..lambda.len() {
        let inner = part(mu, row) as usize;
        for col in inner..lambda[row] as usize {
            cells.push((row, col));
        }
    }

    let outer = Partition::from_sorted(lambda.to_vec());
    let inner = Partition::from_sorted(mu.to_vec());
    let mut out = Vec::new();
    enumerate_skew_cell(
        lambda,
        mu,
        lower_flags,
        upper_flags,
        &cells,
        &mut filling,
        0,
        &outer,
        &inner,
        &mut out,
    );
    out
}

/// Enumerate semistandard skew tableaux of shape `lambda / mu` with upper row
/// flags.
pub fn flagged_skew_tableaux(lambda: &[u32], mu: &[u32], flags: &[u32]) -> Vec<SkewTableau> {
    let lower = vec![1u32; flags.len()];
    row_interval_flagged_skew_tableaux(lambda, mu, &lower, flags)
}

/// Flagged skew Schur polynomial with upper row flags.
pub fn flagged_skew_schur<C: Ring>(
    lambda: &[u32],
    mu: &[u32],
    flags: &[u32],
    num_vars: usize,
) -> MultiPoly<C> {
    let lower = vec![1u32; flags.len()];
    row_interval_flagged_skew_schur(lambda, mu, &lower, flags, num_vars)
}

/// Flagged skew Schur polynomial with row intervals.
pub fn row_interval_flagged_skew_schur<C: Ring>(
    lambda: &[u32],
    mu: &[u32],
    lower_flags: &[u32],
    upper_flags: &[u32],
    num_vars: usize,
) -> MultiPoly<C> {
    validate_skew_shape(lambda, mu);
    validate_row_intervals(lambda.len(), lower_flags, upper_flags, None);

    let max_width = lambda.iter().copied().max().unwrap_or(0) as usize;
    let mut filling = vec![vec![0u32; max_width]; lambda.len()];
    let mut cells = Vec::new();
    for row in 0..lambda.len() {
        let inner = part(mu, row) as usize;
        for col in inner..lambda[row] as usize {
            cells.push((row, col));
        }
    }

    let mut exponents = vec![0u32; num_vars];
    let mut terms = BTreeMap::new();
    accumulate_skew_tableau_weights_cell(
        lambda,
        mu,
        lower_flags,
        upper_flags,
        &cells,
        &mut filling,
        0,
        &mut exponents,
        &mut terms,
    );
    MultiPoly::from_terms(num_vars, terms)
}

fn enumerate_tableaux_cell(
    lambda: &[u32],
    lower_flags: &[u32],
    upper_flags: &[u32],
    filling: &mut [Vec<u32>],
    row: usize,
    col: usize,
    out: &mut Vec<Vec<Vec<u32>>>,
) {
    if row == lambda.len() {
        out.push(filling.to_vec());
        return;
    }

    if col == lambda[row] as usize {
        enumerate_tableaux_cell(lambda, lower_flags, upper_flags, filling, row + 1, 0, out);
        return;
    }

    let row_min = if col == 0 {
        lower_flags[row]
    } else {
        filling[row][col - 1]
    };
    let col_min = if row == 0 || col >= lambda[row - 1] as usize {
        lower_flags[row]
    } else {
        filling[row - 1][col] + 1
    };
    let min_entry = row_min.max(col_min);

    for entry in min_entry..=upper_flags[row] {
        filling[row][col] = entry;
        enumerate_tableaux_cell(lambda, lower_flags, upper_flags, filling, row, col + 1, out);
    }
    filling[row][col] = 0;
}

#[allow(clippy::too_many_arguments)]
fn accumulate_tableau_weights_cell<C: Ring>(
    lambda: &[u32],
    lower_flags: &[u32],
    upper_flags: &[u32],
    filling: &mut [Vec<u32>],
    row: usize,
    col: usize,
    exponents: &mut [u32],
    terms: &mut BTreeMap<Vec<u32>, C>,
) {
    if row == lambda.len() {
        add_weight(exponents, terms);
        return;
    }

    if col == lambda[row] as usize {
        accumulate_tableau_weights_cell(
            lambda,
            lower_flags,
            upper_flags,
            filling,
            row + 1,
            0,
            exponents,
            terms,
        );
        return;
    }

    let row_min = if col == 0 {
        lower_flags[row]
    } else {
        filling[row][col - 1]
    };
    let col_min = if row == 0 || col >= lambda[row - 1] as usize {
        lower_flags[row]
    } else {
        filling[row - 1][col] + 1
    };
    let min_entry = row_min.max(col_min);

    for entry in min_entry..=upper_flags[row] {
        let exponent_index = entry as usize - 1;
        assert!(
            exponent_index < exponents.len(),
            "num_vars ({}) must be at least the maximum tableau entry",
            exponents.len()
        );
        filling[row][col] = entry;
        exponents[exponent_index] = exponents[exponent_index]
            .checked_add(1)
            .expect("tableau weight exponent exceeds u32");
        accumulate_tableau_weights_cell(
            lambda,
            lower_flags,
            upper_flags,
            filling,
            row,
            col + 1,
            exponents,
            terms,
        );
        exponents[exponent_index] -= 1;
    }
    filling[row][col] = 0;
}

#[allow(clippy::too_many_arguments)]
fn enumerate_skew_cell(
    lambda: &[u32],
    mu: &[u32],
    lower_flags: &[u32],
    upper_flags: &[u32],
    cells: &[(usize, usize)],
    filling: &mut [Vec<u32>],
    cell_index: usize,
    outer: &Partition,
    inner: &Partition,
    out: &mut Vec<SkewTableau>,
) {
    if cell_index == cells.len() {
        let rows = (0..lambda.len())
            .map(|row| {
                let start = part(mu, row) as usize;
                let end = lambda[row] as usize;
                (start..end).map(|col| filling[row][col]).collect()
            })
            .collect::<Vec<Vec<u32>>>();
        out.push(SkewTableau::new(outer.clone(), inner.clone(), rows));
        return;
    }

    let (row, col) = cells[cell_index];
    let row_min = if col == part(mu, row) as usize {
        lower_flags[row]
    } else {
        filling[row][col - 1]
    };
    let col_min = if row == 0 || col < part(mu, row - 1) as usize || lambda[row - 1] as usize <= col
    {
        lower_flags[row]
    } else {
        filling[row - 1][col] + 1
    };
    let min_entry = row_min.max(col_min);

    for entry in min_entry..=upper_flags[row] {
        filling[row][col] = entry;
        enumerate_skew_cell(
            lambda,
            mu,
            lower_flags,
            upper_flags,
            cells,
            filling,
            cell_index + 1,
            outer,
            inner,
            out,
        );
    }
    filling[row][col] = 0;
}

#[allow(clippy::too_many_arguments)]
fn accumulate_skew_tableau_weights_cell<C: Ring>(
    lambda: &[u32],
    mu: &[u32],
    lower_flags: &[u32],
    upper_flags: &[u32],
    cells: &[(usize, usize)],
    filling: &mut [Vec<u32>],
    cell_index: usize,
    exponents: &mut [u32],
    terms: &mut BTreeMap<Vec<u32>, C>,
) {
    if cell_index == cells.len() {
        add_weight(exponents, terms);
        return;
    }

    let (row, col) = cells[cell_index];
    let row_min = if col == part(mu, row) as usize {
        lower_flags[row]
    } else {
        filling[row][col - 1]
    };
    let col_min = if row == 0 || col < part(mu, row - 1) as usize || lambda[row - 1] as usize <= col
    {
        lower_flags[row]
    } else {
        filling[row - 1][col] + 1
    };
    let min_entry = row_min.max(col_min);

    for entry in min_entry..=upper_flags[row] {
        let exponent_index = entry as usize - 1;
        assert!(
            exponent_index < exponents.len(),
            "num_vars ({}) must be at least the maximum tableau entry",
            exponents.len()
        );
        filling[row][col] = entry;
        exponents[exponent_index] = exponents[exponent_index]
            .checked_add(1)
            .expect("tableau weight exponent exceeds u32");
        accumulate_skew_tableau_weights_cell(
            lambda,
            mu,
            lower_flags,
            upper_flags,
            cells,
            filling,
            cell_index + 1,
            exponents,
            terms,
        );
        exponents[exponent_index] -= 1;
    }
    filling[row][col] = 0;
}

fn add_weight<C: Ring>(exponents: &[u32], terms: &mut BTreeMap<Vec<u32>, C>) {
    let coefficient = terms.entry(exponents.to_vec()).or_insert_with(C::zero);
    *coefficient = coefficient.clone() + C::one();
}

#[cfg(test)]
fn tableau_weight_enumerator<C: Ring>(tableaux: &[Tableau], num_vars: usize) -> MultiPoly<C> {
    let mut terms = BTreeMap::new();
    for tableau in tableaux {
        let mut exp = vec![0u32; num_vars];
        for row in tableau.rows() {
            for &entry in row {
                let idx = entry as usize - 1;
                assert!(
                    idx < num_vars,
                    "num_vars ({num_vars}) must be at least the maximum tableau entry"
                );
                exp[idx] += 1;
            }
        }
        let coeff = terms.entry(exp).or_insert_with(C::zero);
        *coeff = coeff.clone() + C::one();
    }
    MultiPoly::from_terms(num_vars, terms)
}

#[cfg(test)]
fn skew_tableau_weight_enumerator<C: Ring>(
    tableaux: &[SkewTableau],
    num_vars: usize,
) -> MultiPoly<C> {
    let mut terms = BTreeMap::new();
    for tableau in tableaux {
        let mut exp = vec![0u32; num_vars];
        for row in tableau.rows() {
            for &entry in row {
                let idx = entry as usize - 1;
                assert!(
                    idx < num_vars,
                    "num_vars ({num_vars}) must be at least the maximum tableau entry"
                );
                exp[idx] += 1;
            }
        }
        let coeff = terms.entry(exp).or_insert_with(C::zero);
        *coeff = coeff.clone() + C::one();
    }
    MultiPoly::from_terms(num_vars, terms)
}

fn validate_shape(lambda: &[u32]) {
    assert!(
        lambda.windows(2).all(|pair| pair[0] >= pair[1]),
        "lambda must be weakly decreasing"
    );
}

fn validate_skew_shape(lambda: &[u32], mu: &[u32]) {
    validate_shape(lambda);
    validate_shape(mu);
    assert!(
        mu.len() <= lambda.len(),
        "mu must have at most as many rows as lambda"
    );
    for row in 0..lambda.len() {
        assert!(
            part(mu, row) <= lambda[row],
            "inner shape must fit inside outer shape"
        );
    }
}

fn validate_row_intervals(
    row_count: usize,
    lower_flags: &[u32],
    upper_flags: &[u32],
    num_vars: Option<usize>,
) {
    assert_eq!(lower_flags.len(), row_count);
    assert_eq!(upper_flags.len(), row_count);
    for (&lower, &upper) in lower_flags.iter().zip(upper_flags) {
        assert!(lower >= 1, "row bounds must be positive");
        assert!(
            lower <= upper,
            "each lower row bound must be at most its upper bound"
        );
        if let Some(num_vars) = num_vars {
            assert!(upper as usize <= num_vars, "row bound exceeds num_vars");
        }
    }
}

fn part(parts: &[u32], row: usize) -> u32 {
    parts.get(row).copied().unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use crate::basis::MultiPolyBasis;
    use crate::multipoly_function::MultiPolyFunction;

    use super::*;

    fn assert_straight_weights_match_materialized_tableaux(
        lambda: &[u32],
        lower_flags: &[u32],
        upper_flags: &[u32],
        num_vars: usize,
    ) {
        let tableaux = row_interval_flagged_tableaux(lambda, lower_flags, upper_flags);
        let expected: MultiPoly<i64> = tableau_weight_enumerator(&tableaux, num_vars);
        let actual: MultiPoly<i64> =
            row_interval_flagged_schur(lambda, lower_flags, upper_flags, num_vars);
        assert_eq!(actual.terms(), expected.terms());
    }

    fn assert_skew_weights_match_materialized_tableaux(
        lambda: &[u32],
        mu: &[u32],
        lower_flags: &[u32],
        upper_flags: &[u32],
        num_vars: usize,
    ) {
        let tableaux = row_interval_flagged_skew_tableaux(lambda, mu, lower_flags, upper_flags);
        let expected: MultiPoly<i64> = skew_tableau_weight_enumerator(&tableaux, num_vars);
        let actual: MultiPoly<i64> =
            row_interval_flagged_skew_schur(lambda, mu, lower_flags, upper_flags, num_vars);
        assert_eq!(actual.terms(), expected.terms());
    }

    #[test]
    fn flagged_schur_shape_21_flags_23() {
        let f: MultiPoly<i64> = flagged_schur(&[2, 1], &[2, 3], 3);
        let expected = BTreeMap::from([
            (vec![0, 2, 1], 1),
            (vec![1, 1, 1], 1),
            (vec![1, 2, 0], 1),
            (vec![2, 0, 1], 1),
            (vec![2, 1, 0], 1),
        ]);
        assert_eq!(f.terms(), &expected);
        assert_eq!(flagged_tableaux(&[2, 1], &[2, 3]).len(), 5);
    }

    #[test]
    fn row_interval_flagged_schur_respects_lower_bounds() {
        let f: MultiPoly<i64> = row_interval_flagged_schur(&[2], &[2], &[3], 3);
        let expected = BTreeMap::from([(vec![0, 0, 2], 1), (vec![0, 1, 1], 1), (vec![0, 2, 0], 1)]);
        assert_eq!(f.terms(), &expected);
    }

    #[test]
    fn streamed_straight_weights_match_public_materialized_enumeration() {
        assert_straight_weights_match_materialized_tableaux(&[2, 1], &[1, 1], &[2, 3], 3);
        assert_straight_weights_match_materialized_tableaux(&[3, 1], &[2, 1], &[3, 4], 4);
        assert_straight_weights_match_materialized_tableaux(&[], &[], &[], 0);
        assert_straight_weights_match_materialized_tableaux(&[0], &[1], &[1], 1);
    }

    #[test]
    fn ordinary_skew_schur_specialization_shape_32_1() {
        let f: MultiPoly<i64> = flagged_skew_schur(&[3, 2], &[1], &[3, 3], 3);
        let expected = BTreeMap::from([
            (vec![0, 1, 3], 1),
            (vec![0, 2, 2], 2),
            (vec![0, 3, 1], 1),
            (vec![1, 0, 3], 1),
            (vec![1, 1, 2], 3),
            (vec![1, 2, 1], 3),
            (vec![1, 3, 0], 1),
            (vec![2, 0, 2], 2),
            (vec![2, 1, 1], 3),
            (vec![2, 2, 0], 2),
            (vec![3, 0, 1], 1),
            (vec![3, 1, 0], 1),
        ]);
        assert_eq!(f.terms(), &expected);
    }

    #[test]
    fn streamed_skew_weights_match_public_materialized_enumeration() {
        assert_skew_weights_match_materialized_tableaux(&[3, 2], &[1], &[1, 2], &[2, 3], 3);
        assert_skew_weights_match_materialized_tableaux(&[3, 2], &[1, 1], &[2, 1], &[3, 3], 3);
        assert_skew_weights_match_materialized_tableaux(&[], &[], &[], &[], 0);
        assert_skew_weights_match_materialized_tableaux(&[2, 1], &[2, 1], &[1, 1], &[2, 2], 2);
    }

    #[test]
    fn flagged_skew_schur_is_key_positive_in_example() {
        let f: MultiPoly<i64> = flagged_skew_schur(&[3, 2], &[1], &[2, 3], 3);
        let key_expansion = MultiPolyFunction::from_multipoly(&f).to_key_basis();
        assert_eq!(key_expansion.basis(), MultiPolyBasis::Key);
        assert!(key_expansion.positive_coefficients());
    }
}
