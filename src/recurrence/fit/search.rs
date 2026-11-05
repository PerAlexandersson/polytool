//! Fixed-shape recurrence search: build and solve the exact linear system for
//! one `RecurrenceOptions` shape and decode the solution into a `Recurrence`.

use super::*;
use crate::linalg;
use num_traits::One;
use num_traits::Zero;

/// Coefficients used while constructing a fit system. Integer source data goes
/// directly to the fraction-free solver; fractional data retains the rational path.
trait FitEntry: Clone + Zero + std::ops::Neg<Output = Self> + std::ops::AddAssign {
    fn product(a: &BigRational, b: &BigRational) -> Self;
    fn solve(matrix: Vec<Vec<Self>>, rhs: Vec<Self>, square: bool) -> Option<Vec<BigRational>>;
}

impl FitEntry for BigInt {
    fn product(a: &BigRational, b: &BigRational) -> Self {
        debug_assert!(a.denom().is_one() && b.denom().is_one());
        a.numer() * b.numer()
    }

    fn solve(matrix: Vec<Vec<Self>>, rhs: Vec<Self>, square: bool) -> Option<Vec<BigRational>> {
        linalg::solve_integer_fit_system(matrix, rhs, square)
    }
}

impl FitEntry for BigRational {
    fn product(a: &BigRational, b: &BigRational) -> Self {
        if a.denom().is_one() && b.denom().is_one() {
            BigRational::new_raw(a.numer() * b.numer(), BigInt::one())
        } else {
            a * b
        }
    }

    fn solve(matrix: Vec<Vec<Self>>, rhs: Vec<Self>, square: bool) -> Option<Vec<BigRational>> {
        if square {
            linalg::solve_full_rank_square_linear_system(&matrix, &rhs)
        } else {
            linalg::solve_linear_system(&matrix, &rhs)
        }
    }
}

pub(in crate::recurrence) fn rational_derivatives_up_to(
    polys: &[Vec<BigRational>],
    max_diff_deg: usize,
) -> Vec<Vec<Vec<BigRational>>> {
    polys
        .iter()
        .map(|p| {
            let mut derivatives = Vec::with_capacity(max_diff_deg + 1);
            derivatives.push(p.clone());
            for d in 1..=max_diff_deg {
                let next = poly_derivative_once_rational(&derivatives[d - 1]);
                derivatives.push(next);
            }
            derivatives
        })
        .collect()
}

/// Search for a polynomial recurrence satisfied by a sequence of polynomials.
///
/// `polys[i]` is P_{i+1}(t) given as a coefficient vector (index = power of t).
/// Returns `None` if no recurrence is found with the given options.
pub fn find_polynomial_recurrence_rational(
    polys: &[Vec<BigRational>],
    opts: &RecurrenceOptions,
) -> Option<Recurrence> {
    let derivs = rational_derivatives_up_to(polys, opts.diff_deg);
    find_polynomial_recurrence_rational_with_derivs(polys, &derivs, opts)
}

fn find_polynomial_recurrence_rational_with_derivs(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    opts: &RecurrenceOptions,
) -> Option<Recurrence> {
    find_polynomial_recurrence_rational_with_derivs_and_rows(polys, derivs, opts, None)
}

pub(super) fn find_polynomial_recurrence_rational_with_derivs_and_rows(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    opts: &RecurrenceOptions,
    exact_row_indices: Option<&[usize]>,
) -> Option<Recurrence> {
    find_polynomial_recurrence_in_leader_chart(
        polys,
        derivs,
        opts,
        exact_row_indices,
        LeaderChart::CONSTANT,
    )
}

/// The normalization chart of the left-hand factor
/// `f(n, t) = sum_(i,j) f_(i,j) n^i t^j` of a candidate: the coefficient
/// `f_(index_power, x_power)` is fixed to 1 and every other coefficient of
/// the shape is unknown (`n` is the 1-based search index).  The constant
/// chart `(0, 0)` is the historical normalization; the others reach the
/// left-hand factors whose constant term is zero, such as `x g(n)` or
/// `n g(n, x)` (module documentation of `adaptive`, "Left-hand charts").
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(super) struct LeaderChart {
    /// The power of the search index `n` of the fixed monomial.
    pub(super) index_power: usize,
    /// The power of `x` of the fixed monomial.
    pub(super) x_power: usize,
}

impl LeaderChart {
    /// The historical chart: the constant term of `f` is 1.
    pub(super) const CONSTANT: Self = Self {
        index_power: 0,
        x_power: 0,
    };

    pub(super) fn is_constant(self) -> bool {
        self == Self::CONSTANT
    }

    /// The flat index `i (var_deg + 1) + j` of the fixed monomial.
    fn flat(self, denom_w: usize) -> usize {
        self.index_power * denom_w + self.x_power
    }
}

/// The fixed-shape search in a left-hand chart: as
/// [`find_polynomial_recurrence_rational`], with the coefficient of
/// `n^a x^e` of `f(n, t)` fixed to 1 instead of the constant term, for the
/// chart `(a, e)` inside the shape of `f` (`a <= denom_idx_deg`,
/// `e <= denom_var_deg`).  The constant chart is the historical search,
/// including its modular prefilter.
pub(super) fn find_polynomial_recurrence_in_leader_chart(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    opts: &RecurrenceOptions,
    exact_row_indices: Option<&[usize]>,
    chart: LeaderChart,
) -> Option<Recurrence> {
    let m = polys.len();
    if m <= opts.rec_len {
        return None;
    }
    if chart.index_power > opts.denom_idx_deg || chart.x_power > opts.denom_var_deg {
        return None;
    }

    // The modular prefilter models the constant chart only (and returns
    // without a verdict when `leader_x_power > 0`).
    if chart.is_constant() && modular_prefilter_rejects(polys, derivs, opts) {
        return None;
    }
    let integral = polys.iter().flatten().all(|c| c.denom().is_one())
        && derivs
            .iter()
            .flatten()
            .flatten()
            .all(|c| c.denom().is_one());
    if integral {
        fit_in_leader_chart::<BigInt>(polys, derivs, opts, exact_row_indices, chart)
    } else {
        fit_in_leader_chart::<BigRational>(polys, derivs, opts, exact_row_indices, chart)
    }
}

fn fit_in_leader_chart<C: FitEntry>(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    opts: &RecurrenceOptions,
    exact_row_indices: Option<&[usize]>,
    chart: LeaderChart,
) -> Option<Recurrence> {
    let m = polys.len();
    // The current row enters as `x^e P_n` for the left-hand monomial `x^e`.
    let leader_x_power = opts.leader_x_power;
    let shifted_rows: Vec<Vec<BigRational>>;
    let current_rows: &[Vec<BigRational>] = if leader_x_power == 0 {
        polys
    } else {
        shifted_rows = polys
            .iter()
            .map(|p| {
                let mut shifted = vec![BigRational::zero(); leader_x_power];
                shifted.extend(p.iter().cloned());
                shifted
            })
            .collect();
        &shifted_rows
    };

    // --- Assign column indices to unknowns ---

    // 1) Denominator unknowns d[i][j], skipping the chart monomial (the
    //    constant term (0,0) for the historical chart), which is fixed to 1.
    let denom_start: usize = 0;
    let denom_w = opts.denom_var_deg + 1; // width in j
                                          // Total slots minus 1 for the fixed chart monomial.
    let num_denom_vars = (opts.denom_idx_deg + 1) * denom_w - 1;
    let pinned = chart.flat(denom_w);

    let denom_col = |i: usize, j: usize| -> usize {
        let flat = i * denom_w + j;
        // Skip the pinned slot.
        denom_start + if flat < pinned { flat } else { flat - 1 }
    };

    // 2) Recurrence coefficient unknowns c[r][d][sign][i][j].
    let coeff_start = denom_start + num_denom_vars;
    let vars_per_coeff = (opts.idx_deg + 1) * (opts.var_deg + 1);
    let sign_family_count = if opts.alternating_sign { 2 } else { 1 };
    let num_coeff_vars = opts.rec_len * (opts.diff_deg + 1) * sign_family_count * vars_per_coeff;

    let coeff_col = |r: usize, d: usize, sign_idx: usize, i: usize, j: usize| -> usize {
        coeff_start
            + (((r - 1) * (opts.diff_deg + 1) + d) * sign_family_count + sign_idx) * vars_per_coeff
            + i * (opts.var_deg + 1)
            + j
    };

    // 3) Inhomogeneous unknowns (if needed).
    let inhomo_start = coeff_start + num_coeff_vars;
    let inhomo_w = opts.inhomo_var_deg + 1;
    let num_inhomo_vars = if opts.homogeneous {
        0
    } else {
        (opts.inhomo_idx_deg + 1) * inhomo_w
    };

    let inhomo_col = |i: usize, j: usize| -> usize { inhomo_start + i * inhomo_w + j };

    let num_vars = inhomo_start + num_inhomo_vars;
    if num_vars == 0 {
        return None;
    }

    // --- Determine max t-degree across all equations ---
    let max_poly_deg = polys
        .iter()
        .map(|p| poly_degree_rational(p))
        .max()
        .unwrap_or(0);
    let max_j = opts
        .var_deg
        .max(opts.denom_var_deg + leader_x_power)
        .max(if opts.homogeneous {
            0
        } else {
            opts.inhomo_var_deg
        });
    let max_t_deg = max_j + max_poly_deg;
    let eqs_per_nn = max_t_deg + 1;

    // --- Build linear system Ax = b ---
    let num_nn = m - opts.rec_len; // nn = rec_len+1 ..= m  (1-based)
    let num_rows = num_nn * eqs_per_nn;
    let all_row_indices;
    let row_indices = if let Some(row_indices) = exact_row_indices {
        if row_indices.is_empty() {
            return None;
        }
        row_indices
    } else {
        all_row_indices = (0..num_rows).collect::<Vec<_>>();
        &all_row_indices
    };

    let zero = BigRational::zero();
    let mut matrix = vec![vec![C::zero(); num_vars]; row_indices.len()];
    let mut rhs = vec![C::zero(); row_indices.len()];
    let max_n_pow = opts
        .idx_deg
        .max(opts.denom_idx_deg)
        .max(if opts.homogeneous {
            0
        } else {
            opts.inhomo_idx_deg
        });
    let n_power_table = (opts.rec_len + 1..=m)
        .map(|nn| rational_index_powers(nn, max_n_pow))
        .collect::<Vec<_>>();

    for (row, &equation_row) in row_indices.iter().enumerate() {
        if equation_row >= num_rows {
            return None;
        }
        let eq_idx = equation_row / eqs_per_nn;
        let l = equation_row % eqs_per_nn;
        let nn = opts.rec_len + 1 + eq_idx;
        // nn is 1-based; polys[nn-1] is P_nn(t).
        let current = &current_rows[nn - 1];
        let n_powers = &n_power_table[eq_idx];

        // RHS = the fixed chart monomial times P_nn: the coefficient of
        // t^l in n^a t^e P_nn(t) (moved to RHS with negation).
        rhs[row] = if l >= chart.x_power {
            -C::product(
                current.get(l - chart.x_power).unwrap_or(&zero),
                &n_powers[chart.index_power],
            )
        } else {
            C::zero()
        };

        // Denominator unknowns: d[i][j] * nn^i * (coeff of t^{l-j} in P_nn).
        if num_denom_vars > 0 {
            for i in 0..=opts.denom_idx_deg {
                for j in 0..=opts.denom_var_deg {
                    if i == chart.index_power && j == chart.x_power {
                        continue;
                    }
                    if l < j {
                        continue;
                    }
                    let Some(pc) = current.get(l - j).filter(|pc| !pc.is_zero()) else {
                        continue;
                    };
                    matrix[row][denom_col(i, j)] = C::product(pc, &n_powers[i]);
                }
            }
        }

        // Recurrence coefficients:
        // −sign(nn) * c[r][d][sign][i][j] * nn^i
        //     * (coeff of t^{l-j} in D^d P_{nn-r}).
        for r in 1..=opts.rec_len {
            for d in 0..=opts.diff_deg {
                let ref_poly = &derivs[nn - 1 - r][d];
                for sign_idx in 0..sign_family_count {
                    for i in 0..=opts.idx_deg {
                        let flip = sign_idx == 1 && nn % 2 == 1;
                        for j in 0..=opts.var_deg {
                            if l < j {
                                continue;
                            }
                            let Some(rc) = ref_poly.get(l - j).filter(|rc| !rc.is_zero()) else {
                                continue;
                            };
                            // -(rc * ni), with ni = (-1)^flip n^i.
                            matrix[row][coeff_col(r, d, sign_idx, i, j)] = if flip {
                                C::product(rc, &n_powers[i])
                            } else {
                                -C::product(rc, &n_powers[i])
                            };
                        }
                    }
                }
            }
        }

        // Inhomogeneous unknowns: −c_inh[i][j] * nn^i * delta(l,j).
        if !opts.homogeneous {
            for i in 0..=opts.inhomo_idx_deg {
                if l <= opts.inhomo_var_deg {
                    let val = -C::product(&BigRational::one(), &n_powers[i]);
                    matrix[row][inhomo_col(i, l)] += val;
                }
            }
        }
    }

    // --- Solve ---
    let square = exact_row_indices.is_some() && matrix.len() == num_vars;
    let solution = C::solve(matrix, rhs, square)?;

    // NOTE: No normalization here. The system is non-homogeneous (the fixed
    // constant 1 in f(n,t) pins the scale), so the solution is uniquely
    // determined — not up to scaling.  Dividing by GCD would lose actual
    // coefficient values (e.g. turning P(n) = 2 P(n-1) into P(n) = P(n-1)).

    let mut recurrence = recurrence_from_solution_vector_in_chart(&solution, opts, chart)?;
    if leader_x_power > 0 {
        recurrence.denominator = Some(multiply_by_x_power(
            recurrence.denominator.as_ref(),
            leader_x_power,
        ));
    }

    if exact_row_indices.is_some()
        && !recurrence_holds_rational_with_derivs(polys, derivs, &recurrence, opts.rec_len)
    {
        return None;
    }

    Some(recurrence)
}

/// Search for a polynomial recurrence satisfied by an integer-coefficient
/// polynomial sequence.
///
/// This is a convenience wrapper around [`find_polynomial_recurrence_rational`].
pub fn find_polynomial_recurrence(
    polys: &[Vec<i64>],
    opts: &RecurrenceOptions,
) -> Option<Recurrence> {
    let rational_polys = i64_polys_to_rational(polys);
    find_polynomial_recurrence_rational(&rational_polys, opts)
}

pub(in crate::recurrence) fn recurrence_from_solution_vector(
    solution: &[BigRational],
    opts: &RecurrenceOptions,
) -> Option<Recurrence> {
    recurrence_from_solution_vector_in_chart(solution, opts, LeaderChart::CONSTANT)
}

/// Decode a solution vector of the system in the left-hand chart `chart`
/// (the chart monomial of the denominator is 1, the others are unknowns).
fn recurrence_from_solution_vector_in_chart(
    solution: &[BigRational],
    opts: &RecurrenceOptions,
    chart: LeaderChart,
) -> Option<Recurrence> {
    let denom_start: usize = 0;
    let denom_w = opts.denom_var_deg + 1;
    let num_denom_vars = (opts.denom_idx_deg + 1) * denom_w - 1;
    let pinned = chart.flat(denom_w);
    let denom_col = |i: usize, j: usize| -> usize {
        let flat = i * denom_w + j;
        denom_start + if flat < pinned { flat } else { flat - 1 }
    };

    let coeff_start = denom_start + num_denom_vars;
    let vars_per_coeff = (opts.idx_deg + 1) * (opts.var_deg + 1);
    let sign_family_count = if opts.alternating_sign { 2 } else { 1 };
    let num_coeff_vars = opts.rec_len * (opts.diff_deg + 1) * sign_family_count * vars_per_coeff;
    let coeff_col = |r: usize, d: usize, sign_idx: usize, i: usize, j: usize| -> usize {
        coeff_start
            + (((r - 1) * (opts.diff_deg + 1) + d) * sign_family_count + sign_idx) * vars_per_coeff
            + i * (opts.var_deg + 1)
            + j
    };

    let inhomo_start = coeff_start + num_coeff_vars;
    let inhomo_w = opts.inhomo_var_deg + 1;
    let num_inhomo_vars = if opts.homogeneous {
        0
    } else {
        (opts.inhomo_idx_deg + 1) * inhomo_w
    };
    let inhomo_col = |i: usize, j: usize| -> usize { inhomo_start + i * inhomo_w + j };

    let num_vars = inhomo_start + num_inhomo_vars;
    if solution.len() != num_vars {
        return None;
    }

    let all_zero = (coeff_start..coeff_start + num_coeff_vars).all(|c| solution[c].is_zero());
    if all_zero {
        return None;
    }

    let mut terms = Vec::new();
    for r in 1..=opts.rec_len {
        for d in 0..=opts.diff_deg {
            for sign_idx in 0..sign_family_count {
                let bv = extract_bivar(
                    solution,
                    |i, j| coeff_col(r, d, sign_idx, i, j),
                    opts.idx_deg,
                    opts.var_deg,
                );
                if !bv.is_zero() {
                    terms.push(RecurrenceTerm {
                        offset: r,
                        deriv_order: d,
                        sign: RecurrenceSign::from_family_index(sign_idx),
                        coeff: bv,
                    });
                }
            }
        }
    }

    let denominator = if num_denom_vars > 0 {
        let mut bv = extract_bivar(
            solution,
            |i, j| {
                if i == chart.index_power && j == chart.x_power {
                    0
                } else {
                    denom_col(i, j)
                }
            },
            opts.denom_idx_deg,
            opts.denom_var_deg,
        );
        bv.coeffs[chart.index_power][chart.x_power] = BigRational::one();
        if bv.is_one() {
            None
        } else {
            Some(bv)
        }
    } else {
        None
    };

    let inhomogeneous = if !opts.homogeneous {
        let bv = extract_bivar(
            solution,
            inhomo_col,
            opts.inhomo_idx_deg,
            opts.inhomo_var_deg,
        );
        if bv.is_zero() {
            None
        } else {
            Some(bv)
        }
    } else {
        None
    };

    Some(Recurrence {
        terms,
        denominator,
        inhomogeneous,
    })
}

// ---------------------------------------------------------------------------
// Helpers for building the result
// ---------------------------------------------------------------------------

/// `x^power f(n, t)` for `f` (the constant 1 when absent).
fn multiply_by_x_power(factor: Option<&BivarPoly>, power: usize) -> BivarPoly {
    let rows = match factor {
        Some(f) => f.coeffs.clone(),
        None => vec![vec![BigRational::one()]],
    };
    let shifted = rows
        .into_iter()
        .map(|row| {
            let mut out = vec![BigRational::zero(); power];
            out.extend(row);
            out
        })
        .collect();
    BivarPoly::from_coefficients_normalized(shifted)
}

/// Extract a BivarPoly from a solution vector using a column-index mapping.
fn extract_bivar(
    solution: &[BigRational],
    col_fn: impl Fn(usize, usize) -> usize,
    max_i: usize,
    max_j: usize,
) -> BivarPoly {
    let coeffs: Vec<Vec<BigRational>> = (0..=max_i)
        .map(|i| {
            (0..=max_j)
                .map(|j| solution[col_fn(i, j)].clone())
                .collect()
        })
        .collect();
    // Trailing zero rows and columns are kept; the result is not trimmed.
    BivarPoly { coeffs }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::test_support::*;

    #[test]
    fn integer_and_rational_fit_paths_agree_across_charts() {
        let polys = i64_polys_to_rational(&[
            vec![1],
            vec![1, 1],
            vec![1, 2, 1],
            vec![1, 3, 3, 1],
            vec![1, 4, 6, 4, 1],
            vec![1, 5, 10, 10, 5, 1],
        ]);
        let derivs = rational_derivatives_up_to(&polys, 1);
        let mut found = 0;
        for flags in 0..32 {
            let opts = RecurrenceOptions {
                rec_len: 1 + (flags & 1),
                idx_deg: 1,
                var_deg: 1,
                diff_deg: (flags >> 1) & 1,
                homogeneous: flags & 4 == 0,
                alternating_sign: flags & 8 != 0,
                leader_x_power: (flags >> 4) & 1,
                denom_idx_deg: 1,
                denom_var_deg: 1,
                inhomo_idx_deg: 1,
                inhomo_var_deg: 1,
                ..Default::default()
            };
            for chart in [
                LeaderChart::CONSTANT,
                LeaderChart {
                    index_power: 1,
                    x_power: 0,
                },
                LeaderChart {
                    index_power: 0,
                    x_power: 1,
                },
            ] {
                let rational =
                    fit_in_leader_chart::<BigRational>(&polys, &derivs, &opts, None, chart);
                let integer = fit_in_leader_chart::<BigInt>(&polys, &derivs, &opts, None, chart);
                assert_eq!(
                    format!("{integer:?}"),
                    format!("{rational:?}"),
                    "flags {flags}, chart {chart:?}"
                );
                found += usize::from(integer.is_some());
            }
        }
        assert!(found > 0);

        // A selected square system with an integral input and a fractional solution.
        let polys = i64_polys_to_rational(&[vec![16], vec![8], vec![4], vec![2], vec![1]]);
        let derivs = rational_derivatives_up_to(&polys, 0);
        let opts = RecurrenceOptions {
            rec_len: 1,
            idx_deg: 0,
            var_deg: 0,
            diff_deg: 0,
            ..Default::default()
        };
        for selected in [None, Some(&[0][..]), Some(&[0, 0][..])] {
            let integer = fit_in_leader_chart::<BigInt>(
                &polys,
                &derivs,
                &opts,
                selected,
                LeaderChart::CONSTANT,
            );
            let rational = fit_in_leader_chart::<BigRational>(
                &polys,
                &derivs,
                &opts,
                selected,
                LeaderChart::CONSTANT,
            );
            assert_eq!(format!("{integer:?}"), format!("{rational:?}"));
            assert_eq!(integer.unwrap().to_string(), "P(n) = 1/2 P(n-1)");
        }
    }

    #[test]
    fn rational_input_geometric() {
        let polys: Vec<Vec<BigRational>> = vec![
            vec![br(1, 2)],
            vec![br(1, 4)],
            vec![br(1, 8)],
            vec![br(1, 16)],
            vec![br(1, 32)],
        ];
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            ..Default::default()
        };
        let rec = find_polynomial_recurrence_rational(&polys, &opts)
            .expect("should find rational recurrence");
        assert_eq!(format!("{rec}"), "P(n) = 1/2 P(n-1)");
    }

    #[test]
    fn rational_input_large_coefficients() {
        let base = BigInt::one() << 100usize;
        let polys: Vec<Vec<BigRational>> = (0..5)
            .map(|i| vec![BigRational::from_integer(&base << i)])
            .collect();
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            ..Default::default()
        };
        let rec = find_polynomial_recurrence_rational(&polys, &opts)
            .expect("should find recurrence with large coefficients");
        assert_eq!(format!("{rec}"), "P(n) = 2 P(n-1)");
    }

    #[test]
    fn fibonacci() {
        // P_n = P_{n-1} + P_{n-2}
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![2],
            vec![3],
            vec![5],
            vec![8],
            vec![13],
        ];
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 2,
            homogeneous: true,
            ..Default::default()
        };
        assert_recurrence(&polys, &opts, "P(n) = P(n-1) + P(n-2)");
    }

    #[test]
    fn binomial_expansion() {
        // P_n = (1+t)^n, so P_n = (1+t) P_{n-1}
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1, 1],
            vec![1, 2, 1],
            vec![1, 3, 3, 1],
            vec![1, 4, 6, 4, 1],
        ];
        let opts = RecurrenceOptions {
            var_deg: 1,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            ..Default::default()
        };
        assert_recurrence(&polys, &opts, "P(n) = (1 + t) P(n-1)");
    }

    #[test]
    fn factorial() {
        // P_n = n! (as constant polynomials): P_n = n P_{n-1}
        // Using 1-based indexing: P_1=1, P_2=2, P_3=6, P_4=24, P_5=120
        let polys: Vec<Vec<i64>> = vec![vec![1], vec![2], vec![6], vec![24], vec![120]];
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 1,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            ..Default::default()
        };
        assert_recurrence(&polys, &opts, "P(n) = n P(n-1)");
    }

    #[test]
    fn chebyshev() {
        // T_n(t) = 2t T_{n-1}(t) - T_{n-2}(t)
        // T_0=1, T_1=t, T_2=2t^2-1, T_3=4t^3-3t, T_4=8t^4-8t^2+1
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![0, 1],
            vec![-1, 0, 2],
            vec![0, -3, 0, 4],
            vec![1, 0, -8, 0, 8],
        ];
        let opts = RecurrenceOptions {
            var_deg: 1,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 2,
            homogeneous: true,
            ..Default::default()
        };
        assert_recurrence(&polys, &opts, "P(n) = 2t P(n-1) - P(n-2)");
    }

    #[test]
    fn eulerian_with_derivative() {
        // A_n(t) = (1 + (n-1)t) A_{n-1}(t) + t(1-t) A'_{n-1}(t)
        // A_1=1, A_2=1, A_3=1+t, A_4=1+4t+t^2, A_5=1+11t+11t^2+t^3
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![1, 1],
            vec![1, 4, 1],
            vec![1, 11, 11, 1],
            vec![1, 26, 66, 26, 1],
        ];
        let opts = RecurrenceOptions {
            var_deg: 2,
            idx_deg: 1,
            diff_deg: 1,
            rec_len: 1,
            homogeneous: true,
            ..Default::default()
        };
        let rec = find_polynomial_recurrence(&polys, &opts).expect("should find recurrence");
        let display = format!("{rec}");
        // With 1-based indexing: A_0=P_1, A_1=P_2, etc.
        // So the recurrence becomes P_n = [1 + (n-2)t] P_{n-1} + t(1-t) P'_{n-1}
        // Coefficient of P(n-1): 1 - 2t + nt
        // Coefficient of P'(n-1): t - t^2
        assert_eq!(display, "P(n) = (1 - 2t + nt) P(n-1) + (t - t^2) P'(n-1)");
    }

    #[test]
    fn geometric_sequence() {
        // P_n = 2^{n-1}: P_n = 2 P_{n-1}
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![2],
            vec![4],
            vec![8],
            vec![16],
            vec![32],
            vec![64],
            vec![128],
            vec![256],
            vec![512],
        ];
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            ..Default::default()
        };
        assert_recurrence(&polys, &opts, "P(n) = 2 P(n-1)");
    }

    #[test]
    fn with_alternating_sign_builder_toggles_the_search_space() {
        // P_n = (-1)^n P_{n-1} is found exactly when alternating terms are allowed
        // (same sequence as `alternating_sign_terms`).
        let polys: Vec<Vec<i64>> = [1, 1, -1, -1, 1, 1, -1, -1]
            .iter()
            .map(|&c| vec![c])
            .collect();
        let base = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            ..Default::default()
        };
        assert!(!base.alternating_sign);
        let enabled = base.clone().with_alternating_sign(true);
        assert!(enabled.alternating_sign);
        assert!(find_polynomial_recurrence(&polys, &base).is_none());
        let rec = find_polynomial_recurrence(&polys, &enabled).expect("alternating recurrence");
        assert!(rec.terms.iter().any(RecurrenceTerm::is_alternating));
        // The builder only changes that flag, and can switch it off again.
        let disabled = enabled.clone().with_alternating_sign(false);
        assert!(!disabled.alternating_sign);
        assert_eq!(disabled.rec_len, base.rec_len);
        assert!(find_polynomial_recurrence(&polys, &disabled).is_none());
    }

    #[test]
    fn alternating_sign_terms() {
        // P_n = (-1)^n P_{n-1}; with P_1=1 this gives pairs 1,1,-1,-1,...
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![-1],
            vec![-1],
            vec![1],
            vec![1],
            vec![-1],
            vec![-1],
        ];
        let ordinary_opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            ..Default::default()
        };
        assert!(find_polynomial_recurrence(&polys, &ordinary_opts).is_none());

        let alternating_opts = RecurrenceOptions {
            alternating_sign: true,
            ..ordinary_opts
        };
        let rec = find_polynomial_recurrence(&polys, &alternating_opts)
            .expect("should find alternating-sign recurrence");
        assert_eq!(format!("{rec}"), "P(n) = (-1)^n P(n-1)");
        assert_eq!(rec.to_latex(), "P(n) = (-1)^n P(n-1)");
        assert!(rec
            .to_mathematica_definition(&polys)
            .contains("Expand[(-1)^n*P[n - 1, t]]"));
        assert!(rec
            .to_sage_definition(&polys)
            .contains("value = R((-1)**n*P(n - 1))"));
    }

    #[test]
    fn constant_coeffs_with_separate_inhomogeneous_bounds() {
        // P_n = P_{n-1} + n + t^2
        let polys: Vec<Vec<i64>> = vec![
            vec![0],
            vec![2, 0, 1],
            vec![5, 0, 2],
            vec![9, 0, 3],
            vec![14, 0, 4],
            vec![20, 0, 5],
        ];
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: false,
            inhomo_var_deg: 2,
            inhomo_idx_deg: 1,
            ..Default::default()
        };
        assert_recurrence(&polys, &opts, "P(n) = P(n-1) + n + t^2");
    }

    #[test]
    fn no_recurrence() {
        // Random polynomials unlikely to satisfy a short recurrence.
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1, 7],
            vec![3, 1, 5],
            vec![2, 9, 1, 4],
            vec![7, 2, 8, 1, 3],
        ];
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 2,
            homogeneous: true,
            ..Default::default()
        };
        assert!(find_polynomial_recurrence(&polys, &opts).is_none());
    }

    #[test]
    fn index_dependent_coefficients() {
        // P_n = (2n-1) P_{n-1} - (n-1)^2 P_{n-2}  (Legendre-related)
        // P_1=1, P_2=1, then P_3 = 3*1 - 1*1 = 2, P_4 = 5*2 - 4*1 = 6,
        // P_5 = 7*6 - 9*2 = 24, P_6 = 9*24 - 16*6 = 120
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![2],
            vec![6],
            vec![24],
            vec![120],
            vec![720],
        ];
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 2,
            diff_deg: 0,
            rec_len: 2,
            homogeneous: true,
            ..Default::default()
        };
        let rec = find_polynomial_recurrence(&polys, &opts);
        // This should find SOME recurrence (the simplest might be P_n = n P_{n-1}).
        assert!(rec.is_some());
    }
}
