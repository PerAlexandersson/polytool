//! Fixed-shape recurrence search: build and solve the exact linear system for
//! one `RecurrenceOptions` shape and decode the solution into a `Recurrence`.

use super::*;
use crate::linalg;
use num_traits::One;
use num_traits::Zero;

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
    let m = polys.len();
    if m <= opts.rec_len {
        return None;
    }

    if modular_prefilter_rejects(polys, derivs, opts) {
        return None;
    }

    // --- Assign column indices to unknowns ---

    // 1) Denominator unknowns d[i][j], skipping (0,0) which is the fixed constant 1.
    let denom_start: usize = 0;
    let denom_w = opts.denom_var_deg + 1; // width in j
                                          // Total slots minus 1 for the fixed constant-1 at (0,0).
    let num_denom_vars = (opts.denom_idx_deg + 1) * denom_w - 1;

    let denom_col = |i: usize, j: usize| -> usize {
        let flat = i * denom_w + j; // (0,0) has flat=0
        denom_start + flat - 1 // skip (0,0)
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
        .max(opts.denom_var_deg)
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
    let mut matrix: Vec<Vec<BigRational>> = vec![vec![zero.clone(); num_vars]; row_indices.len()];
    let mut rhs: Vec<BigRational> = vec![zero.clone(); row_indices.len()];
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
        let current = &polys[nn - 1];
        let n_powers = &n_power_table[eq_idx];

        // RHS = coefficient of t^l in P_nn(t)  (moved to RHS with negation).
        let cur_l = poly_coeff_rational(current, l);
        rhs[row] = -cur_l;

        // Denominator unknowns: d[i][j] * nn^i * (coeff of t^{l-j} in P_nn).
        if num_denom_vars > 0 {
            for i in 0..=opts.denom_idx_deg {
                for j in 0..=opts.denom_var_deg {
                    if i == 0 && j == 0 {
                        continue;
                    }
                    if l < j {
                        continue;
                    }
                    let pc = poly_coeff_rational(current, l - j);
                    if pc.is_zero() {
                        continue;
                    }
                    let val = pc * n_powers[i].clone();
                    matrix[row][denom_col(i, j)] = val;
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
                        let mut ni = n_powers[i].clone();
                        if sign_idx == 1 && nn % 2 == 1 {
                            ni = -ni;
                        }
                        for j in 0..=opts.var_deg {
                            if l < j {
                                continue;
                            }
                            let rc = poly_coeff_rational(ref_poly, l - j);
                            if rc.is_zero() {
                                continue;
                            }
                            let val = -rc * ni.clone();
                            matrix[row][coeff_col(r, d, sign_idx, i, j)] = val;
                        }
                    }
                }
            }
        }

        // Inhomogeneous unknowns: −c_inh[i][j] * nn^i * delta(l,j).
        if !opts.homogeneous {
            for i in 0..=opts.inhomo_idx_deg {
                if l <= opts.inhomo_var_deg {
                    let val = -n_powers[i].clone();
                    matrix[row][inhomo_col(i, l)] += val;
                }
            }
        }
    }

    // --- Solve ---
    let solution = if exact_row_indices.is_some() && matrix.len() == num_vars {
        linalg::solve_full_rank_square_linear_system(&matrix, &rhs)?
    } else {
        linalg::solve_linear_system(&matrix, &rhs)?
    };

    // NOTE: No normalization here. The system is non-homogeneous (the fixed
    // constant 1 in f(n,t) pins the scale), so the solution is uniquely
    // determined — not up to scaling.  Dividing by GCD would lose actual
    // coefficient values (e.g. turning P(n) = 2 P(n-1) into P(n) = P(n-1)).

    let recurrence = recurrence_from_solution_vector(&solution, opts)?;

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
    let denom_start: usize = 0;
    let denom_w = opts.denom_var_deg + 1;
    let num_denom_vars = (opts.denom_idx_deg + 1) * denom_w - 1;
    let denom_col = |i: usize, j: usize| -> usize {
        let flat = i * denom_w + j;
        denom_start + flat - 1
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
                if i == 0 && j == 0 {
                    0
                } else {
                    denom_col(i, j)
                }
            },
            opts.denom_idx_deg,
            opts.denom_var_deg,
        );
        bv.coeffs[0][0] = BigRational::one();
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
