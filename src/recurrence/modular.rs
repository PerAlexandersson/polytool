//! Modular prefiltering and lifting: candidate systems are first solved
//! modulo large primes, and lifted solutions are always re-verified exactly.

use super::*;
use num_traits::One;
use num_traits::ToPrimitive;

use crate::linalg;

const MODULAR_PREFILTER_PRIMES: [i64; 3] = [1_000_000_007, 1_000_000_009, 998_244_353];
const MODULAR_PREFILTER_PRECOMPUTED_PRIMES: usize = 2;
// Full-rank modular fits first check one held-out row, which cheaply rejects
// many false prefix fits.  A full held-out sweep is reserved for larger spaces.
const MODULAR_PREFILTER_FIRST_TAIL_VERIFY_MIN_UNKNOWNS: usize = 1;
const MODULAR_PREFILTER_TAIL_VERIFY_MIN_UNKNOWNS: usize = 128;
const MODULAR_LIFT_RANK_DEFICIENT_MIN_UNKNOWNS: usize = 16;

fn mod_norm(value: i64, modulus: i64) -> i64 {
    value.rem_euclid(modulus)
}

fn mod_neg(value: i64, modulus: i64) -> i64 {
    if value == 0 {
        0
    } else {
        modulus - value
    }
}

fn mod_mul(lhs: i64, rhs: i64, modulus: i64) -> i64 {
    if (0..modulus).contains(&lhs) && (0..modulus).contains(&rhs) && modulus <= 3_000_000_000 {
        return (lhs * rhs) % modulus;
    }
    ((lhs as i128 * rhs as i128).rem_euclid(modulus as i128)) as i64
}

fn mod_pow(mut base: i64, mut exp: i64, modulus: i64) -> i64 {
    let mut acc = 1;
    base = mod_norm(base, modulus);
    while exp > 0 {
        if exp & 1 == 1 {
            acc = mod_mul(acc, base, modulus);
        }
        base = mod_mul(base, base, modulus);
        exp >>= 1;
    }
    acc
}

fn mod_index_powers(n: usize, max_power: usize, modulus: i64) -> Vec<i64> {
    let base = mod_norm(n as i64, modulus);
    let mut powers = Vec::with_capacity(max_power + 1);
    let mut current = 1;
    for _ in 0..=max_power {
        powers.push(current);
        current = mod_mul(current, base, modulus);
    }
    powers
}

pub(super) fn rational_index_powers(n: usize, max_power: usize) -> Vec<BigRational> {
    let base = BigInt::from(n);
    let mut powers = Vec::with_capacity(max_power + 1);
    let mut current = BigInt::one();
    for _ in 0..=max_power {
        powers.push(BigRational::from_integer(current.clone()));
        current *= &base;
    }
    powers
}

fn mod_inv(value: i64, modulus: i64) -> Option<i64> {
    let value = mod_norm(value, modulus);
    (value != 0).then(|| mod_pow(value, modulus - 2, modulus))
}

fn bigint_mod_i64(value: &BigInt, modulus: i64) -> i64 {
    let modulus_big = BigInt::from(modulus);
    let residue = ((value % &modulus_big) + &modulus_big) % &modulus_big;
    residue
        .to_i64()
        .expect("residue modulo i64 prime fits in i64")
}

fn rational_mod_prime(value: &BigRational, modulus: i64) -> Option<i64> {
    let numerator = bigint_mod_i64(value.numer(), modulus);
    let denominator = bigint_mod_i64(value.denom(), modulus);
    let denominator_inv = mod_inv(denominator, modulus)?;
    Some(mod_mul(numerator, denominator_inv, modulus))
}

fn poly_coeff_mod(coeffs: &[Vec<i64>], poly_idx: usize, coeff_idx: usize) -> i64 {
    coeffs
        .get(poly_idx)
        .and_then(|coeffs| coeffs.get(coeff_idx))
        .copied()
        .unwrap_or(0)
}

fn poly_coeff_mod_slice(coeffs: &[i64], coeff_idx: usize) -> i64 {
    coeffs.get(coeff_idx).copied().unwrap_or(0)
}

fn add_modular_solution_term(
    lhs: &mut i64,
    solution: &[u64],
    col: usize,
    coeff: i64,
    modulus: i64,
) {
    let coeff = mod_norm(coeff, modulus);
    if coeff == 0 {
        return;
    }
    let term = mod_mul(coeff, solution[col] as i64, modulus);
    *lhs += term;
    if *lhs >= modulus {
        *lhs -= modulus;
    }
}

fn rational_polys_mod_prime(polys: &[Vec<BigRational>], modulus: i64) -> Option<Vec<Vec<i64>>> {
    polys
        .iter()
        .map(|poly| {
            poly.iter()
                .map(|coeff| rational_mod_prime(coeff, modulus))
                .collect()
        })
        .collect()
}

fn rational_derivs_mod_prime(
    derivs: &[Vec<Vec<BigRational>>],
    modulus: i64,
    max_diff_deg: usize,
) -> Option<Vec<Vec<Vec<i64>>>> {
    derivs
        .iter()
        .map(|poly_derivs| {
            poly_derivs
                .iter()
                .take(max_diff_deg + 1)
                .map(|poly| {
                    poly.iter()
                        .map(|coeff| rational_mod_prime(coeff, modulus))
                        .collect()
                })
                .collect()
        })
        .collect()
}

fn mod_poly_derivative_once(coeffs: &[i64], modulus: i64) -> Vec<i64> {
    if coeffs.len() <= 1 {
        return vec![0];
    }
    coeffs[1..]
        .iter()
        .enumerate()
        .map(|(idx, &coeff)| mod_mul(coeff, (idx + 1) as i64, modulus))
        .collect()
}

fn modular_derivatives_up_to(
    polys_mod: &[Vec<i64>],
    modulus: i64,
    max_diff_deg: usize,
) -> Vec<Vec<Vec<i64>>> {
    polys_mod
        .iter()
        .map(|poly| {
            let mut derivatives = Vec::with_capacity(max_diff_deg + 1);
            derivatives.push(poly.clone());
            for d in 1..=max_diff_deg {
                let next = mod_poly_derivative_once(&derivatives[d - 1], modulus);
                derivatives.push(next);
            }
            derivatives
        })
        .collect()
}

enum ModularPolySupport {
    Dense,
    Sparse(Vec<(usize, i64)>),
}

impl ModularPolySupport {
    #[inline]
    fn for_each_nonzero(&self, coeffs: &[i64], mut visit: impl FnMut(usize, i64)) {
        match self {
            Self::Dense => {
                for (idx, &coeff) in coeffs.iter().enumerate() {
                    visit(idx, coeff);
                }
            }
            Self::Sparse(entries) => {
                for &(idx, coeff) in entries {
                    visit(idx, coeff);
                }
            }
        }
    }
}

fn modular_poly_support(coeffs: &[i64]) -> ModularPolySupport {
    if coeffs.iter().all(|&coeff| coeff != 0) {
        return ModularPolySupport::Dense;
    }

    ModularPolySupport::Sparse(
        coeffs
            .iter()
            .copied()
            .enumerate()
            .filter(|&(_, coeff)| coeff != 0)
            .collect(),
    )
}

fn modular_poly_supports(polys_mod: &[Vec<i64>]) -> Vec<ModularPolySupport> {
    polys_mod
        .iter()
        .map(|poly| modular_poly_support(poly))
        .collect()
}

fn modular_derivative_supports(derivs_mod: &[Vec<Vec<i64>>]) -> Vec<Vec<ModularPolySupport>> {
    derivs_mod
        .iter()
        .map(|poly_derivs| {
            poly_derivs
                .iter()
                .map(|poly| modular_poly_support(poly))
                .collect()
        })
        .collect()
}

struct ModularImageView<'a> {
    polys: &'a [Vec<i64>],
    derivs: &'a [Vec<Vec<i64>>],
    poly_supports: &'a [ModularPolySupport],
    deriv_supports: &'a [Vec<ModularPolySupport>],
}

struct ModularSystemResult {
    consistent: bool,
    coefficient_full_rank: bool,
    full_rank_pivot_rows: Option<Vec<usize>>,
    modular_solution: Option<Vec<u64>>,
}

#[derive(Clone, Debug)]
pub(super) struct ModularCandidateSolution {
    modulus: i64,
    solution: Vec<u64>,
    pub(super) full_rank: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ModularSolveRequest {
    ConsistencyOnly,
    FullRankSolution,
    RankDeficientCandidateAllowed,
}

impl ModularSolveRequest {
    fn solution_mode(self) -> linalg::SparseModSolutionMode {
        match self {
            Self::ConsistencyOnly => linalg::SparseModSolutionMode::None,
            Self::FullRankSolution => linalg::SparseModSolutionMode::FullRank,
            Self::RankDeficientCandidateAllowed => linalg::SparseModSolutionMode::AnyConsistent,
        }
    }

    fn needs_solution_data(self) -> bool {
        self != Self::ConsistencyOnly
    }
}

fn recurrence_system_consistent_mod_images(
    polys: &[Vec<BigRational>],
    polys_mod: &[Vec<i64>],
    derivs_mod: &[Vec<Vec<i64>>],
    opts: &RecurrenceOptions,
    modulus: i64,
    solve_request: ModularSolveRequest,
) -> ModularSystemResult {
    let poly_supports = modular_poly_supports(polys_mod);
    let deriv_supports = modular_derivative_supports(derivs_mod);
    recurrence_system_consistent_mod_supports(
        polys,
        ModularImageView {
            polys: polys_mod,
            derivs: derivs_mod,
            poly_supports: &poly_supports,
            deriv_supports: &deriv_supports,
        },
        opts,
        modulus,
        solve_request,
    )
}

fn recurrence_system_consistent_mod_supports(
    polys: &[Vec<BigRational>],
    images: ModularImageView<'_>,
    opts: &RecurrenceOptions,
    modulus: i64,
    solve_request: ModularSolveRequest,
) -> ModularSystemResult {
    let m = polys.len();
    if m <= opts.rec_len {
        return ModularSystemResult {
            consistent: true,
            coefficient_full_rank: false,
            full_rank_pivot_rows: None,
            modular_solution: None,
        };
    }

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
    if num_vars == 0 {
        return ModularSystemResult {
            consistent: true,
            coefficient_full_rank: true,
            full_rank_pivot_rows: None,
            modular_solution: None,
        };
    }
    let prime = modulus as u64;

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
    let num_nn = m - opts.rec_len;
    let num_rows = num_nn * eqs_per_nn;
    let max_n_pow = opts
        .idx_deg
        .max(opts.denom_idx_deg)
        .max(if opts.homogeneous {
            0
        } else {
            opts.inhomo_idx_deg
        });
    let n_power_table = (opts.rec_len + 1..=m)
        .map(|nn| mod_index_powers(nn, max_n_pow, modulus))
        .collect::<Vec<_>>();

    let mut raw_rows = Vec::with_capacity(num_rows);
    let row_capacity_hint = num_vars.min(8);
    for nn in opts.rec_len + 1..=m {
        let current_idx = nn - 1;
        for l in 0..=max_t_deg {
            let cur_l = poly_coeff_mod(images.polys, current_idx, l);
            raw_rows.push(linalg::SparseModRow::with_capacity(
                mod_neg(cur_l, modulus) as u64,
                prime,
                row_capacity_hint,
            ));
        }
    }
    let row_index = |eq_idx: usize, l: usize| -> usize { eq_idx * eqs_per_nn + l };
    let power_indices = n_power_table.first().map(Vec::as_slice).unwrap_or(&[]);

    if num_denom_vars > 0 {
        for (i, _) in power_indices
            .iter()
            .enumerate()
            .take(opts.denom_idx_deg + 1)
        {
            for j in 0..=opts.denom_var_deg {
                if i == 0 && j == 0 {
                    continue;
                }
                let col = denom_col(i, j);
                for (eq_idx, nn) in (opts.rec_len + 1..=m).enumerate() {
                    let n_power = n_power_table[eq_idx][i];
                    if n_power == 0 {
                        continue;
                    }
                    images.poly_supports[nn - 1].for_each_nonzero(
                        &images.polys[nn - 1],
                        |k, pc| {
                            let l = k + j;
                            if l <= max_t_deg {
                                raw_rows[row_index(eq_idx, l)].add_reduced_entry_sorted(
                                    col,
                                    mod_mul(pc, n_power, modulus) as u64,
                                    prime,
                                );
                            }
                        },
                    );
                }
            }
        }
    }

    let derivative_order_indices = images.derivs.first().map(Vec::as_slice).unwrap_or(&[]);
    for r in 1..=opts.rec_len {
        for (d, _) in derivative_order_indices
            .iter()
            .enumerate()
            .take(opts.diff_deg + 1)
        {
            for sign_idx in 0..sign_family_count {
                for (i, _) in power_indices.iter().enumerate().take(opts.idx_deg + 1) {
                    for j in 0..=opts.var_deg {
                        let col = coeff_col(r, d, sign_idx, i, j);
                        for (eq_idx, nn) in (opts.rec_len + 1..=m).enumerate() {
                            let n_power = n_power_table[eq_idx][i];
                            let n_factor = if sign_idx == 1 && nn % 2 == 1 {
                                mod_neg(n_power, modulus)
                            } else {
                                n_power
                            };
                            if n_factor == 0 {
                                continue;
                            }
                            images.deriv_supports[nn - 1 - r][d].for_each_nonzero(
                                &images.derivs[nn - 1 - r][d],
                                |k, rc| {
                                    let l = k + j;
                                    if l <= max_t_deg {
                                        raw_rows[row_index(eq_idx, l)].add_reduced_entry_sorted(
                                            col,
                                            mod_neg(mod_mul(rc, n_factor, modulus), modulus) as u64,
                                            prime,
                                        );
                                    }
                                },
                            );
                        }
                    }
                }
            }
        }
    }

    if !opts.homogeneous {
        for i in 0..=opts.inhomo_idx_deg {
            for j in 0..=opts.inhomo_var_deg {
                if j > max_t_deg {
                    continue;
                }
                let col = inhomo_col(i, j);
                for (eq_idx, n_powers) in n_power_table.iter().enumerate() {
                    let n_power = n_powers[i];
                    if n_power != 0 {
                        raw_rows[row_index(eq_idx, j)].add_reduced_entry_sorted(
                            col,
                            mod_neg(n_power, modulus) as u64,
                            prime,
                        );
                    }
                }
            }
        }
    }

    let mut rows = Vec::with_capacity(num_rows);
    let mut row_indices = Vec::with_capacity(num_rows);
    for (equation_row, row) in raw_rows.into_iter().enumerate() {
        if row.is_contradiction() {
            return ModularSystemResult {
                consistent: false,
                coefficient_full_rank: false,
                full_rank_pivot_rows: None,
                modular_solution: None,
            };
        }
        if !row.is_tautology() {
            row_indices.push(equation_row);
            rows.push(row);
        }
    }

    let row_order = linalg::SparseModRowOrder::IncreasingNonzeros;
    let result =
        linalg::sparse_modular_linear_system_consistency_prime_unchecked_with_solution_mode(
            rows,
            num_vars,
            prime,
            row_order,
            solve_request.solution_mode(),
        )
        .expect("recurrence modular prefilter builds a well-formed prime-field system");
    let coefficient_full_rank = result.rank == num_vars;
    let full_rank_pivot_rows = (solve_request.needs_solution_data()
        && result.consistent
        && coefficient_full_rank)
        .then(|| {
            result
                .pivot_rows
                .iter()
                .map(|&row| row_indices[row])
                .collect()
        });
    ModularSystemResult {
        consistent: result.consistent,
        coefficient_full_rank,
        full_rank_pivot_rows,
        modular_solution: result.solution,
    }
}

fn recurrence_solution_holds_from_mod_images(
    polys: &[Vec<BigRational>],
    polys_mod: &[Vec<i64>],
    derivs_mod: &[Vec<Vec<i64>>],
    opts: &RecurrenceOptions,
    modulus: i64,
    solution: &[u64],
    start_nn: usize,
) -> bool {
    let m = polys.len();
    let first_nn = start_nn.max(opts.rec_len + 1);
    if first_nn > m {
        return true;
    }

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
        return false;
    }

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
    let max_n_pow = opts
        .idx_deg
        .max(opts.denom_idx_deg)
        .max(if opts.homogeneous {
            0
        } else {
            opts.inhomo_idx_deg
        });

    for nn in first_nn..=m {
        let current_idx = nn - 1;
        let n_powers = mod_index_powers(nn, max_n_pow, modulus);

        for l in 0..=max_t_deg {
            let cur_l = poly_coeff_mod(polys_mod, current_idx, l);
            let rhs = mod_neg(cur_l, modulus);
            let mut lhs = 0;

            if num_denom_vars > 0 {
                for (i, &n_power) in n_powers.iter().enumerate().take(opts.denom_idx_deg + 1) {
                    for j in 0..=opts.denom_var_deg {
                        if i == 0 && j == 0 {
                            continue;
                        }
                        if l < j {
                            continue;
                        }
                        let pc = poly_coeff_mod(polys_mod, current_idx, l - j);
                        if pc == 0 {
                            continue;
                        }
                        add_modular_solution_term(
                            &mut lhs,
                            solution,
                            denom_col(i, j),
                            mod_mul(pc, n_power, modulus),
                            modulus,
                        );
                    }
                }
            }

            for r in 1..=opts.rec_len {
                for (d, ref_poly) in derivs_mod[nn - 1 - r]
                    .iter()
                    .enumerate()
                    .take(opts.diff_deg + 1)
                {
                    for sign_idx in 0..sign_family_count {
                        for (i, &n_power) in n_powers.iter().enumerate().take(opts.idx_deg + 1) {
                            let n_factor = if sign_idx == 1 && nn % 2 == 1 {
                                mod_neg(n_power, modulus)
                            } else {
                                n_power
                            };
                            for j in 0..=opts.var_deg {
                                if l < j {
                                    continue;
                                }
                                let rc = poly_coeff_mod_slice(ref_poly, l - j);
                                if rc == 0 {
                                    continue;
                                }
                                add_modular_solution_term(
                                    &mut lhs,
                                    solution,
                                    coeff_col(r, d, sign_idx, i, j),
                                    mod_neg(mod_mul(rc, n_factor, modulus), modulus),
                                    modulus,
                                );
                            }
                        }
                    }
                }
            }

            if !opts.homogeneous && l <= opts.inhomo_var_deg {
                for (i, &n_power) in n_powers.iter().enumerate().take(opts.inhomo_idx_deg + 1) {
                    add_modular_solution_term(
                        &mut lhs,
                        solution,
                        inhomo_col(i, l),
                        mod_neg(n_power, modulus),
                        modulus,
                    );
                }
            }

            if lhs != rhs {
                return false;
            }
        }
    }

    true
}

fn recurrence_system_mod_prime(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    opts: &RecurrenceOptions,
    modulus: i64,
) -> Option<ModularSystemResult> {
    let polys_mod = rational_polys_mod_prime(polys, modulus)?;
    let derivs_mod = rational_derivs_mod_prime(derivs, modulus, opts.diff_deg)?;
    Some(recurrence_system_consistent_mod_images(
        polys,
        &polys_mod,
        &derivs_mod,
        opts,
        modulus,
        ModularSolveRequest::ConsistencyOnly,
    ))
}

struct ModularPrefilterPrimeImages {
    modulus: i64,
    polys: Vec<Vec<i64>>,
    derivs: Vec<Vec<Vec<i64>>>,
    poly_supports: Vec<ModularPolySupport>,
    deriv_supports: Vec<Vec<ModularPolySupport>>,
}

pub(super) struct ModularPrefilterCache {
    primes: Vec<ModularPrefilterPrimeImages>,
}

impl ModularPrefilterCache {
    pub(super) fn new(polys: &[Vec<BigRational>], max_diff_deg: usize) -> Self {
        let primes = MODULAR_PREFILTER_PRIMES
            .iter()
            .take(MODULAR_PREFILTER_PRECOMPUTED_PRIMES)
            .filter_map(|&modulus| {
                let polys_mod = rational_polys_mod_prime(polys, modulus)?;
                let derivs_mod = modular_derivatives_up_to(&polys_mod, modulus, max_diff_deg);
                let poly_supports = modular_poly_supports(&polys_mod);
                let deriv_supports = modular_derivative_supports(&derivs_mod);
                Some(ModularPrefilterPrimeImages {
                    modulus,
                    polys: polys_mod,
                    derivs: derivs_mod,
                    poly_supports,
                    deriv_supports,
                })
            })
            .collect();
        Self { primes }
    }
}

pub(super) fn modular_prefilter_rejects(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    opts: &RecurrenceOptions,
) -> bool {
    if !opts.modular_prefilter {
        return false;
    }

    for &prime in &MODULAR_PREFILTER_PRIMES {
        if let Some(result) = recurrence_system_mod_prime(polys, derivs, opts, prime) {
            // Inconsistency is a sound certificate over Q only when the
            // coefficient matrix retains full column rank modulo this prime.
            // Otherwise the prime may divide a denominator of the rational
            // solution, so exact solving must remain the fallback.
            if !result.consistent && result.coefficient_full_rank {
                return true;
            }
        }
    }

    false
}

#[derive(Default)]
pub(super) struct ModularPrefilterResult {
    pub(super) rejected: bool,
    pub(super) full_rank_pivot_rows: Option<Vec<usize>>,
    pub(super) modular_solution: Option<ModularCandidateSolution>,
}

pub(super) fn modular_prefilter_with_cache(
    polys: &[Vec<BigRational>],
    fit_len: usize,
    opts: &RecurrenceOptions,
    cache: &ModularPrefilterCache,
) -> ModularPrefilterResult {
    if !opts.modular_prefilter {
        return ModularPrefilterResult::default();
    }

    let unknowns = count_unknowns(opts);
    let tail_verify_end =
        if fit_len >= polys.len() || unknowns < MODULAR_PREFILTER_FIRST_TAIL_VERIFY_MIN_UNKNOWNS {
            None
        } else if unknowns >= MODULAR_PREFILTER_TAIL_VERIFY_MIN_UNKNOWNS {
            Some(polys.len())
        } else {
            Some(fit_len + 1)
        };

    let mut rank_deficient_candidate = None;
    for prime_images in &cache.primes {
        if fit_len > polys.len()
            || polys.len() > prime_images.polys.len()
            || polys.len() > prime_images.derivs.len()
        {
            continue;
        }
        if prime_images
            .derivs
            .first()
            .is_some_and(|derivs| derivs.len() <= opts.diff_deg)
        {
            continue;
        }

        let result = recurrence_system_consistent_mod_supports(
            &polys[..fit_len],
            ModularImageView {
                polys: &prime_images.polys[..fit_len],
                derivs: &prime_images.derivs[..fit_len],
                poly_supports: &prime_images.poly_supports[..fit_len],
                deriv_supports: &prime_images.deriv_supports[..fit_len],
            },
            opts,
            prime_images.modulus,
            if unknowns >= MODULAR_LIFT_RANK_DEFICIENT_MIN_UNKNOWNS {
                ModularSolveRequest::RankDeficientCandidateAllowed
            } else {
                ModularSolveRequest::FullRankSolution
            },
        );
        let tail_consistent = match (
            tail_verify_end,
            result.full_rank_pivot_rows.is_some(),
            result.modular_solution.as_ref(),
        ) {
            (Some(verify_len), true, Some(solution)) => recurrence_solution_holds_from_mod_images(
                &polys[..verify_len],
                &prime_images.polys[..verify_len],
                &prime_images.derivs[..verify_len],
                opts,
                prime_images.modulus,
                solution,
                fit_len + 1,
            ),
            _ => true,
        };
        if result.consistent && !tail_consistent {
            // A full-rank modular fit has a unique solution with a nonzero
            // pivot determinant modulo this prime. If a held-out row fails,
            // the exact candidate cannot satisfy the same held-out equation.
            return ModularPrefilterResult {
                rejected: true,
                full_rank_pivot_rows: None,
                modular_solution: None,
            };
        }
        if !result.consistent {
            if result.coefficient_full_rank {
                return ModularPrefilterResult {
                    rejected: true,
                    full_rank_pivot_rows: None,
                    modular_solution: None,
                };
            }
        } else if result.full_rank_pivot_rows.is_none() {
            if unknowns >= MODULAR_LIFT_RANK_DEFICIENT_MIN_UNKNOWNS
                && rank_deficient_candidate.is_none()
            {
                if let Some(solution) = result.modular_solution.clone() {
                    rank_deficient_candidate = Some(ModularCandidateSolution {
                        modulus: prime_images.modulus,
                        solution,
                        full_rank: false,
                    });
                }
            }
            if let Some(verify_len) = tail_verify_end {
                if verify_len > fit_len {
                    // Rank-deficient prefix fits are the expensive false
                    // positives: they have many modular solutions, so no
                    // unique held-out solution can be checked.  Add the
                    // held-out rows modulo the same prime before paying for an
                    // exact rational solve.
                    let extended = recurrence_system_consistent_mod_supports(
                        &polys[..verify_len],
                        ModularImageView {
                            polys: &prime_images.polys[..verify_len],
                            derivs: &prime_images.derivs[..verify_len],
                            poly_supports: &prime_images.poly_supports[..verify_len],
                            deriv_supports: &prime_images.deriv_supports[..verify_len],
                        },
                        opts,
                        prime_images.modulus,
                        ModularSolveRequest::ConsistencyOnly,
                    );
                    if !extended.consistent && extended.coefficient_full_rank {
                        return ModularPrefilterResult {
                            rejected: true,
                            full_rank_pivot_rows: None,
                            modular_solution: None,
                        };
                    }
                }
            }
        } else if let Some(pivot_rows) = result.full_rank_pivot_rows {
            let modular_solution =
                result
                    .modular_solution
                    .map(|solution| ModularCandidateSolution {
                        modulus: prime_images.modulus,
                        solution,
                        full_rank: true,
                    });
            return ModularPrefilterResult {
                rejected: false,
                full_rank_pivot_rows: Some(pivot_rows),
                modular_solution,
            };
        }
    }

    ModularPrefilterResult {
        rejected: false,
        full_rank_pivot_rows: None,
        modular_solution: rank_deficient_candidate,
    }
}

fn centered_modular_rational(residue: u64, modulus: i64) -> BigRational {
    let modulus_u = modulus as u64;
    let centered = if residue > modulus_u / 2 {
        BigInt::from(residue as i128 - modulus as i128)
    } else {
        BigInt::from(residue)
    };
    BigRational::from_integer(centered)
}

fn gcd_i128(mut lhs: i128, mut rhs: i128) -> i128 {
    lhs = lhs.abs();
    rhs = rhs.abs();
    while rhs != 0 {
        let rem = lhs % rhs;
        lhs = rhs;
        rhs = rem;
    }
    lhs
}

fn isqrt_i128(value: i128) -> i128 {
    if value <= 0 {
        return 0;
    }
    let mut lo = 1_i128;
    let mut hi = value;
    let mut ans = 0_i128;
    while lo <= hi {
        let mid = lo + (hi - lo) / 2;
        if mid <= value / mid {
            ans = mid;
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }
    ans
}

fn rational_reconstruction_mod_prime(residue: u64, modulus: i64) -> Option<BigRational> {
    let modulus_i = modulus as i128;
    let bound = isqrt_i128(modulus_i / 2);
    let mut u = (1_i128, 0_i128, modulus_i);
    let mut v = (0_i128, 1_i128, residue as i128);

    while v.2.abs() > bound {
        if v.2 == 0 {
            return None;
        }
        let q = u.2 / v.2;
        let next = (u.0 - q * v.0, u.1 - q * v.1, u.2 - q * v.2);
        u = v;
        v = next;
    }

    let mut numerator = v.2;
    let mut denominator = v.1;
    if denominator == 0 || denominator.abs() > bound {
        return None;
    }
    if gcd_i128(numerator, denominator) != 1 {
        return None;
    }
    if denominator < 0 {
        numerator = -numerator;
        denominator = -denominator;
    }

    let denominator_inv = mod_inv(denominator as i64, modulus)?;
    let reconstructed = mod_mul(
        mod_norm(numerator as i64, modulus),
        denominator_inv,
        modulus,
    );
    if reconstructed != residue as i64 {
        return None;
    }

    Some(BigRational::new(
        BigInt::from(numerator),
        BigInt::from(denominator),
    ))
}

fn recurrence_from_centered_modular_solution(
    candidate: &ModularCandidateSolution,
    opts: &RecurrenceOptions,
) -> Option<Recurrence> {
    let solution = candidate
        .solution
        .iter()
        .map(|&residue| centered_modular_rational(residue, candidate.modulus))
        .collect::<Vec<_>>();
    recurrence_from_solution_vector(&solution, opts)
}

fn recurrence_from_reconstructed_modular_solution(
    candidate: &ModularCandidateSolution,
    opts: &RecurrenceOptions,
) -> Option<Recurrence> {
    let solution = candidate
        .solution
        .iter()
        .map(|&residue| rational_reconstruction_mod_prime(residue, candidate.modulus))
        .collect::<Option<Vec<_>>>()?;
    recurrence_from_solution_vector(&solution, opts)
}

pub(super) fn verified_recurrence_from_modular_lift(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    opts: &RecurrenceOptions,
    candidate: &ModularCandidateSolution,
) -> Option<Recurrence> {
    if let Some(rec) = recurrence_from_centered_modular_solution(candidate, opts) {
        if recurrence_holds_rational_with_derivs(polys, derivs, &rec, opts.rec_len) {
            return Some(rec);
        }
    }

    let rec = recurrence_from_reconstructed_modular_solution(candidate, opts)?;
    recurrence_holds_rational_with_derivs(polys, derivs, &rec, opts.rec_len).then_some(rec)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::test_support::*;

    #[test]
    fn modular_prefilter_rejects_inconsistent_scalar_candidate() {
        let polys = i64_polys_to_rational(&[vec![1], vec![2], vec![5]]);
        let derivs = polys
            .iter()
            .map(|p| vec![poly_nth_derivative_rational(p, 0)])
            .collect::<Vec<_>>();
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            modular_prefilter: true,
            ..Default::default()
        };
        assert!(modular_prefilter_rejects(&polys, &derivs, &opts));
        assert!(find_polynomial_recurrence_rational(&polys, &opts).is_none());
    }

    #[test]
    fn modular_prefilter_preserves_true_recurrence() {
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
            modular_prefilter: true,
            ..Default::default()
        };
        assert_recurrence(&polys, &opts, "P(n) = P(n-1) + P(n-2)");
    }

    #[test]
    fn modular_prefilter_falls_back_exactly_at_bad_primes() {
        let denominator =
            BigInt::from(MODULAR_PREFILTER_PRIMES[0]) * BigInt::from(MODULAR_PREFILTER_PRIMES[1]);
        let polys = vec![
            vec![BigRational::from_integer(denominator.pow(3))],
            vec![BigRational::from_integer(denominator.pow(2))],
            vec![BigRational::from_integer(denominator.clone())],
            vec![BigRational::one()],
        ];
        let derivs = rational_derivatives_up_to(&polys, 0);
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            modular_prefilter: true,
            ..Default::default()
        };

        for &prime in &MODULAR_PREFILTER_PRIMES[..2] {
            let result = recurrence_system_mod_prime(&polys, &derivs, &opts, prime).unwrap();
            assert!(!result.consistent);
            assert!(!result.coefficient_full_rank);
        }
        let recurrence = find_polynomial_recurrence_rational(&polys, &opts).unwrap();
        assert_eq!(
            recurrence.to_string(),
            format!("P(n) = 1/{denominator} P(n-1)")
        );

        let search = AdaptiveSearchOptions {
            min_rec_len: 1,
            max_rec_len: 1,
            min_var_deg: 0,
            max_var_deg: 0,
            min_idx_deg: 0,
            max_idx_deg: 0,
            min_diff_deg: 0,
            max_diff_deg: 0,
            min_margin: 0,
            no_verify: true,
            modular_prefilter: true,
            ..Default::default()
        };
        let adaptive = find_recurrence_adaptive_rational(&polys, &search).unwrap();
        assert_eq!(adaptive.recurrence.to_string(), recurrence.to_string());
    }

    #[test]
    fn modular_tail_verification_checks_solution_on_heldout_rows() {
        let polys = i64_polys_to_rational(&[
            vec![1],
            vec![1],
            vec![2],
            vec![3],
            vec![5],
            vec![8],
            vec![13],
        ]);
        let derivs = rational_derivatives_up_to(&polys, 0);
        let modulus = MODULAR_PREFILTER_PRIMES[0];
        let polys_mod = rational_polys_mod_prime(&polys, modulus).expect("modular polys");
        let derivs_mod =
            rational_derivs_mod_prime(&derivs, modulus, 0).expect("modular derivatives");
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 2,
            homogeneous: true,
            modular_prefilter: true,
            ..Default::default()
        };

        assert!(recurrence_solution_holds_from_mod_images(
            &polys,
            &polys_mod,
            &derivs_mod,
            &opts,
            modulus,
            &[1, 1],
            4,
        ));
        assert!(!recurrence_solution_holds_from_mod_images(
            &polys,
            &polys_mod,
            &derivs_mod,
            &opts,
            modulus,
            &[1, 2],
            4,
        ));
    }

    #[test]
    fn modular_prefilter_falls_back_for_rank_deficient_bad_heldout_row() {
        let polys = i64_polys_to_rational(&[vec![0], vec![0], vec![0], vec![1]]);
        let cache = ModularPrefilterCache::new(&polys, 0);
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            modular_prefilter: true,
            ..Default::default()
        };

        let prefilter = modular_prefilter_with_cache(&polys, 3, &opts, &cache);
        assert!(!prefilter.rejected);
        assert!(prefilter.full_rank_pivot_rows.is_none());
        assert!(find_polynomial_recurrence_rational(&polys, &opts).is_none());
    }

    #[test]
    fn modular_prefilter_tolerates_single_bad_prime() {
        let bad_prime = MODULAR_PREFILTER_PRIMES[0];
        let polys = i64_polys_to_rational(&[vec![bad_prime * bad_prime], vec![bad_prime], vec![1]]);
        let derivs = polys
            .iter()
            .map(|p| vec![poly_nth_derivative_rational(p, 0)])
            .collect::<Vec<_>>();
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            modular_prefilter: true,
            ..Default::default()
        };
        assert!(!modular_prefilter_rejects(&polys, &derivs, &opts));
        assert!(find_polynomial_recurrence_rational(&polys, &opts).is_some());
    }

    #[test]
    fn modular_lift_reconstructs_fractional_solution() {
        let polys = i64_polys_to_rational(&[vec![2], vec![2], vec![2], vec![2], vec![2]]);
        let derivs = rational_derivatives_up_to(&polys, 0);
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: false,
            inhomo_var_deg: 0,
            inhomo_idx_deg: 0,
            ..Default::default()
        };
        let modulus = MODULAR_PREFILTER_PRIMES[0];
        let candidate = ModularCandidateSolution {
            modulus,
            solution: vec![
                rational_mod_prime(&br(1, 2), modulus).unwrap() as u64,
                rational_mod_prime(&BigRational::one(), modulus).unwrap() as u64,
            ],
            full_rank: true,
        };
        let rec = verified_recurrence_from_modular_lift(&polys, &derivs, &opts, &candidate)
            .expect("fractional modular solution should lift");
        assert_eq!(format!("{rec}"), "P(n) = 1/2 P(n-1) + 1");
    }
}
