//! Exact polynomial generators for the deco subexceedant-function model.
//!
//! For a subexceedant function `f(i) <= i`, let `b` be its number of
//! noninitial kernel blocks and let `s` be its number of eligible deco sites.
//! The joint histogram is the coefficient array of
//!
//! `H_N(x, y) = sum_f x^b y^s`.
//!
//! Selecting each eligible site independently gives
//! `P_{N-1}(x; w) = H_N(x, 1 + w)`.  Thus `w = 0` recovers the ordinary
//! Eulerian polynomial `A_N(x)`, while `w = 1` gives the unweighted deco
//! polynomial.  Coefficient vectors use ascending degree order throughout.

use num_bigint::BigInt;

/// A matrix indexed first by noninitial kernel blocks, then by eligible sites.
pub type JointHistogramBigInt = Vec<Vec<BigInt>>;

/// Generate the joint block/eligible-site histograms for sizes `1..=max_size`.
///
/// The return value at index `N - 1` is the coefficient matrix of `H_N(x, y)`:
/// `rows[N - 1][b][s]` counts size-`N` subexceedant functions with `b`
/// noninitial kernel blocks and `s` eligible deco sites.
///
/// The generator uses the exact recurrence
///
/// `H_N = x(1-x) dH_{N-1}/dx + (1+(N-1)x)H_{N-1} + (y-1)xH_{N-2}`.
pub fn joint_histograms_bigint(max_size: usize) -> Vec<JointHistogramBigInt> {
    if max_size == 0 {
        return Vec::new();
    }

    let mut rows = Vec::with_capacity(max_size);
    rows.push(vec![vec![BigInt::from(1)]]);
    if max_size == 1 {
        return rows;
    }

    rows.push(vec![vec![BigInt::from(1)], vec![BigInt::from(1)]]);

    for size in 3..=max_size {
        let site_count = (size - 1) / 2 + 1;
        let mut row = zero_matrix(size, site_count);

        for (blocks, site_row) in rows[size - 2].iter().enumerate() {
            for (sites, count) in site_row.iter().enumerate() {
                row[blocks][sites] += count * BigInt::from(blocks + 1);
                row[blocks + 1][sites] += count * BigInt::from(size - 1 - blocks);
            }
        }

        for (blocks, site_row) in rows[size - 3].iter().enumerate() {
            for (sites, count) in site_row.iter().enumerate() {
                row[blocks + 1][sites + 1] += count;
                row[blocks + 1][sites] -= count;
            }
        }

        debug_assert!(row
            .iter()
            .flatten()
            .all(|coefficient| coefficient >= &BigInt::from(0)));
        rows.push(row);
    }

    rows
}

/// Generate the joint histograms with checked `i64` coefficients.
///
/// This convenience wrapper panics if a coefficient does not fit in `i64`.
pub fn joint_histograms(max_size: usize) -> Vec<Vec<Vec<i64>>> {
    bigint_tables_to_i64(joint_histograms_bigint(max_size))
}

/// Generate `P_0(x; w), ..., P_max_index(x; w)` at an integer value of `w`.
///
/// The specialization `w = 0` is the ordinary Eulerian family, and `w = 1`
/// is the unweighted deco family.  The result at index `m` is `P_m(x; w)`.
pub fn polynomials_bigint(max_index: usize, exception_weight: &BigInt) -> Vec<Vec<BigInt>> {
    let joint_rows = joint_histograms_bigint(max_index + 1);
    let marked_site_weight = exception_weight + BigInt::from(1);

    joint_rows
        .into_iter()
        .map(|histogram| {
            let mut polynomial = vec![BigInt::from(0); histogram.len()];
            for (blocks, site_row) in histogram.into_iter().enumerate() {
                let mut site_power = BigInt::from(1);
                for count in site_row {
                    polynomial[blocks] += count * &site_power;
                    site_power *= &marked_site_weight;
                }
            }
            polynomial
        })
        .collect()
}

/// Generate the integer-specialized deco family with checked `i64` coefficients.
///
/// This convenience wrapper panics if a coefficient does not fit in `i64`.
pub fn polynomials(max_index: usize, exception_weight: i64) -> Vec<Vec<i64>> {
    bigint_polys_to_i64(polynomials_bigint(
        max_index,
        &BigInt::from(exception_weight),
    ))
}

/// Generate the fixed eligible-site layers `K_{N,j}(x)`.
///
/// The outer index is `N - 1`, the middle index is `j`, and the innermost
/// coefficient vector represents
/// `K_{N,j}(x) = x^{-j} [y^j] H_N(x, y)`.
pub fn eligible_site_layers_bigint(max_size: usize) -> Vec<Vec<Vec<BigInt>>> {
    joint_histograms_bigint(max_size)
        .into_iter()
        .enumerate()
        .map(|(size_index, histogram)| {
            let size = size_index + 1;
            (0..=(size - 1) / 2)
                .map(|sites| {
                    let degree = size - 1 - 2 * sites;
                    (0..=degree)
                        .map(|exponent| histogram[exponent + sites][sites].clone())
                        .collect()
                })
                .collect()
        })
        .collect()
}

/// Generate the fixed eligible-site layers with checked `i64` coefficients.
///
/// This convenience wrapper panics if a coefficient does not fit in `i64`.
pub fn eligible_site_layers(max_size: usize) -> Vec<Vec<Vec<i64>>> {
    bigint_tables_to_i64(eligible_site_layers_bigint(max_size))
}

/// Generate the fixed selected-join layers `J_{N,j}(x)`.
///
/// The outer index is `N - 1`, the middle index is `j`, and the innermost
/// coefficient vector represents
/// `J_{N,j}(x) = x^{-j} [w^j] H_N(x, 1 + w)`.
pub fn selected_join_layers_bigint(max_size: usize) -> Vec<Vec<Vec<BigInt>>> {
    joint_histograms_bigint(max_size)
        .into_iter()
        .enumerate()
        .map(|(size_index, histogram)| {
            let size = size_index + 1;
            (0..=(size - 1) / 2)
                .map(|selected| {
                    let degree = size - 1 - 2 * selected;
                    (0..=degree)
                        .map(|exponent| {
                            histogram[exponent + selected]
                                .iter()
                                .enumerate()
                                .skip(selected)
                                .map(|(eligible, count)| {
                                    count * binomial_bigint(eligible, selected)
                                })
                                .sum()
                        })
                        .collect()
                })
                .collect()
        })
        .collect()
}

/// Generate the fixed selected-join layers with checked `i64` coefficients.
///
/// This convenience wrapper panics if a coefficient does not fit in `i64`.
pub fn selected_join_layers(max_size: usize) -> Vec<Vec<Vec<i64>>> {
    bigint_tables_to_i64(selected_join_layers_bigint(max_size))
}

/// Generate gamma coefficients for every fixed selected-join layer.
///
/// The innermost row contains the coefficients in
///
/// `J_{N,j}(x) = sum_k gamma[N,j,k] x^k (1+x)^(N-1-2j-2k)`.
///
/// This uses the positive recurrence
///
/// `gamma[N,j,k] = (j+1+k) gamma[N-1,j,k]`
/// `+ 2(N-2j-2k) gamma[N-1,j,k-1] + gamma[N-2,j-1,k]`.
pub fn selected_join_gamma_layers_bigint(max_size: usize) -> Vec<Vec<Vec<BigInt>>> {
    if max_size == 0 {
        return Vec::new();
    }

    let mut rows = Vec::with_capacity(max_size);
    rows.push(vec![vec![BigInt::from(1)]]);
    if max_size == 1 {
        return rows;
    }
    rows.push(vec![vec![BigInt::from(1)]]);

    for size in 3..=max_size {
        let mut size_rows = Vec::with_capacity((size - 1) / 2 + 1);
        for selected in 0..=(size - 1) / 2 {
            let output_degree = size - 1 - 2 * selected;
            let mut gamma_row = Vec::with_capacity(output_degree / 2 + 1);
            for gamma_index in 0..=output_degree / 2 {
                let diagonal = table_entry(&rows, size - 1, selected, gamma_index)
                    * BigInt::from(selected + 1 + gamma_index);
                let shifted_diagonal = if gamma_index == 0 {
                    BigInt::from(0)
                } else {
                    table_entry(&rows, size - 1, selected, gamma_index - 1)
                        * BigInt::from(2 * (size - 2 * selected - 2 * gamma_index))
                };
                let lag = if selected == 0 {
                    BigInt::from(0)
                } else {
                    table_entry(&rows, size - 2, selected - 1, gamma_index)
                };
                gamma_row.push(diagonal + shifted_diagonal + lag);
            }
            size_rows.push(gamma_row);
        }
        rows.push(size_rows);
    }

    rows
}

/// Generate selected-join gamma layers with checked `i64` coefficients.
///
/// This convenience wrapper panics if a coefficient does not fit in `i64`.
pub fn selected_join_gamma_layers(max_size: usize) -> Vec<Vec<Vec<i64>>> {
    bigint_tables_to_i64(selected_join_gamma_layers_bigint(max_size))
}

fn zero_matrix(rows: usize, columns: usize) -> JointHistogramBigInt {
    vec![vec![BigInt::from(0); columns]; rows]
}

fn binomial_bigint(n: usize, k: usize) -> BigInt {
    if k > n {
        return BigInt::from(0);
    }
    let k = k.min(n - k);
    let mut result = BigInt::from(1);
    for index in 1..=k {
        result = result * BigInt::from(n - k + index) / BigInt::from(index);
    }
    result
}

fn table_entry(
    rows: &[Vec<Vec<BigInt>>],
    size: usize,
    selected: usize,
    gamma_index: usize,
) -> BigInt {
    if size == 0 {
        return BigInt::from(0);
    }
    rows.get(size - 1)
        .and_then(|size_rows| size_rows.get(selected))
        .and_then(|row| row.get(gamma_index))
        .cloned()
        .unwrap_or_else(|| BigInt::from(0))
}

fn bigint_polys_to_i64(polys: Vec<Vec<BigInt>>) -> Vec<Vec<i64>> {
    polys
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|coefficient| {
                    i64::try_from(&coefficient).expect("sequence coefficient overflow")
                })
                .collect()
        })
        .collect()
}

fn bigint_tables_to_i64(tables: Vec<Vec<Vec<BigInt>>>) -> Vec<Vec<Vec<i64>>> {
    tables.into_iter().map(bigint_polys_to_i64).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequences::eulerian_polynomials_bigint;

    fn bigints(values: &[i64]) -> Vec<BigInt> {
        values.iter().copied().map(BigInt::from).collect()
    }

    fn expand_gamma_row(degree: usize, gamma: &[BigInt]) -> Vec<BigInt> {
        let mut result = vec![BigInt::from(0); degree + 1];
        for (gamma_index, coefficient) in gamma.iter().enumerate() {
            let residual_degree = degree - 2 * gamma_index;
            for exponent in 0..=residual_degree {
                result[gamma_index + exponent] +=
                    coefficient * binomial_bigint(residual_degree, exponent);
            }
        }
        result
    }

    #[test]
    fn joint_histogram_matches_first_nontrivial_rows() {
        let rows = joint_histograms_bigint(5);
        assert_eq!(rows[2][1], bigints(&[3, 1]));
        assert_eq!(rows[4][0], bigints(&[1, 0, 0]));
        assert_eq!(rows[4][1], bigints(&[19, 7, 0]));
        assert_eq!(rows[4][2], bigints(&[45, 20, 1]));
        assert_eq!(rows[4][3], bigints(&[19, 7, 0]));
        assert_eq!(rows[4][4], bigints(&[1, 0, 0]));
    }

    #[test]
    fn specializations_recover_eulerian_and_deco_rows() {
        let max_index = 12;
        assert_eq!(
            polynomials_bigint(max_index, &BigInt::from(0)),
            eulerian_polynomials_bigint(max_index + 1)
        );
        assert_eq!(
            polynomials(4, 1),
            vec![
                vec![1],
                vec![1, 1],
                vec![1, 5, 1],
                vec![1, 14, 14, 1],
                vec![1, 33, 89, 33, 1],
            ]
        );
    }

    #[test]
    fn fixed_layers_match_coefficient_extraction() {
        let eligible = eligible_site_layers_bigint(5);
        assert_eq!(eligible[4][0], bigints(&[1, 19, 45, 19, 1]));
        assert_eq!(eligible[4][1], bigints(&[7, 20, 7]));
        assert_eq!(eligible[4][2], bigints(&[1]));

        let selected = selected_join_layers_bigint(5);
        assert_eq!(selected[4][0], bigints(&[1, 26, 66, 26, 1]));
        assert_eq!(selected[4][1], bigints(&[7, 22, 7]));
        assert_eq!(selected[4][2], bigints(&[1]));
    }

    #[test]
    fn gamma_recurrence_expands_to_selected_layers() {
        let max_size = 20;
        let selected = selected_join_layers_bigint(max_size);
        let gamma = selected_join_gamma_layers_bigint(max_size);

        assert_eq!(gamma[4][1], bigints(&[7, 8]));
        assert_eq!(gamma[6][2], bigints(&[25, 20]));
        for size in 1..=max_size {
            for selected_count in 0..=(size - 1) / 2 {
                let degree = size - 1 - 2 * selected_count;
                assert_eq!(
                    expand_gamma_row(degree, &gamma[size - 1][selected_count]),
                    selected[size - 1][selected_count]
                );
            }
        }
    }

    #[test]
    fn checked_wrappers_agree_with_bigint_generators() {
        assert_eq!(
            joint_histograms(8),
            bigint_tables_to_i64(joint_histograms_bigint(8))
        );
        assert_eq!(
            eligible_site_layers(8),
            bigint_tables_to_i64(eligible_site_layers_bigint(8))
        );
        assert_eq!(
            selected_join_layers(8),
            bigint_tables_to_i64(selected_join_layers_bigint(8))
        );
        assert_eq!(
            selected_join_gamma_layers(8),
            bigint_tables_to_i64(selected_join_gamma_layers_bigint(8))
        );
    }
}
