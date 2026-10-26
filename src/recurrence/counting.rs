//! Parameter and equation counts for a candidate recurrence search space.

use super::*;

/// Parameter counts for one candidate recurrence search space.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CandidateComplexity {
    /// Actual number of linear unknowns in the system.
    pub raw_unknowns: usize,
    /// Weighted complexity used only for ordering adaptive candidates.
    pub weighted_unknowns: usize,
    /// Non-alternating, non-derivative recurrence coefficient unknowns.
    pub ordinary_unknowns: usize,
    /// Non-alternating derivative recurrence coefficient unknowns.
    pub derivative_unknowns: usize,
    /// Alternating, non-derivative recurrence coefficient unknowns.
    pub alternating_unknowns: usize,
    /// Alternating derivative recurrence coefficient unknowns.
    pub alternating_derivative_unknowns: usize,
    /// LHS denominator unknowns.
    pub denominator_unknowns: usize,
    /// Inhomogeneous-term unknowns.
    pub inhomogeneous_unknowns: usize,
}

/// Count and weight the parameter slots for a candidate recurrence space.
///
/// The raw count is used for the linear-algebra dimension check. The weighted
/// count is only a heuristic search score: derivative parameters count double,
/// denominator parameters count double, alternating parameters are delayed, and
/// inhomogeneous parameters are delayed further.
pub fn candidate_complexity(opts: &RecurrenceOptions) -> CandidateComplexity {
    let denominator_unknowns = opts
        .denom_idx_deg
        .saturating_add(1)
        .saturating_mul(opts.denom_var_deg.saturating_add(1))
        .saturating_sub(1);
    let vars_per_coeff = opts
        .idx_deg
        .saturating_add(1)
        .saturating_mul(opts.var_deg.saturating_add(1));
    let ordinary_unknowns = opts.rec_len.saturating_mul(vars_per_coeff);
    let derivative_unknowns = opts
        .rec_len
        .saturating_mul(opts.diff_deg)
        .saturating_mul(vars_per_coeff);
    let alternating_unknowns = if opts.alternating_sign {
        opts.rec_len.saturating_mul(vars_per_coeff)
    } else {
        0
    };
    let alternating_derivative_unknowns = if opts.alternating_sign {
        opts.rec_len
            .saturating_mul(opts.diff_deg)
            .saturating_mul(vars_per_coeff)
    } else {
        0
    };
    let inhomogeneous_unknowns = if opts.homogeneous {
        0
    } else {
        opts.inhomo_idx_deg
            .saturating_add(1)
            .saturating_mul(opts.inhomo_var_deg.saturating_add(1))
    };
    let raw_unknowns = ordinary_unknowns
        .saturating_add(derivative_unknowns)
        .saturating_add(alternating_unknowns)
        .saturating_add(alternating_derivative_unknowns)
        .saturating_add(denominator_unknowns)
        .saturating_add(inhomogeneous_unknowns);
    let weighted_unknowns = ordinary_unknowns
        .saturating_add(derivative_unknowns.saturating_mul(2))
        .saturating_add(alternating_unknowns.saturating_mul(3))
        .saturating_add(alternating_derivative_unknowns.saturating_mul(4))
        .saturating_add(denominator_unknowns.saturating_mul(2))
        .saturating_add(inhomogeneous_unknowns.saturating_mul(4));
    CandidateComplexity {
        raw_unknowns,
        weighted_unknowns,
        ordinary_unknowns,
        derivative_unknowns,
        alternating_unknowns,
        alternating_derivative_unknowns,
        denominator_unknowns,
        inhomogeneous_unknowns,
    }
}

/// Count the total number of unknowns for a given set of recurrence options.
pub fn count_unknowns(opts: &RecurrenceOptions) -> usize {
    candidate_complexity(opts).raw_unknowns
}

/// Weighted unknown count used for adaptive candidate ordering.
pub fn count_weighted_unknowns(opts: &RecurrenceOptions) -> usize {
    candidate_complexity(opts).weighted_unknowns
}

/// Count the total number of equations for a given set of recurrence options
/// and rational input polynomials.
pub fn count_equations_rational(polys: &[Vec<BigRational>], opts: &RecurrenceOptions) -> usize {
    let m = polys.len();
    if m <= opts.rec_len {
        return 0;
    }
    let max_poly_deg = max_poly_degree_rational(polys);
    count_equations_from_max_degree(m, max_poly_deg, opts)
}

fn max_poly_degree_rational(polys: &[Vec<BigRational>]) -> usize {
    polys
        .iter()
        .map(|p| poly_degree_rational(p))
        .max()
        .unwrap_or(0)
}

pub(super) fn prefix_max_degrees_rational(polys: &[Vec<BigRational>]) -> Vec<usize> {
    let mut result = Vec::with_capacity(polys.len() + 1);
    result.push(0);
    let mut max_degree = 0;
    for poly in polys {
        max_degree = max_degree.max(poly_degree_rational(poly));
        result.push(max_degree);
    }
    result
}

fn count_equations_from_max_degree(
    poly_count: usize,
    max_poly_deg: usize,
    opts: &RecurrenceOptions,
) -> usize {
    if poly_count <= opts.rec_len {
        return 0;
    }
    let max_j = opts
        .var_deg
        .max(opts.denom_var_deg)
        .max(if opts.homogeneous {
            0
        } else {
            opts.inhomo_var_deg
        });
    let eqs_per_nn = max_j + max_poly_deg + 1;
    (poly_count - opts.rec_len) * eqs_per_nn
}

pub(super) fn count_equations_from_prefix_degrees(
    fit_len: usize,
    opts: &RecurrenceOptions,
    prefix_max_degrees: &[usize],
) -> usize {
    let max_poly_deg = prefix_max_degrees.get(fit_len).copied().unwrap_or(0);
    count_equations_from_max_degree(fit_len, max_poly_deg, opts)
}

/// Count the total number of equations for integer input polynomials.
pub fn count_equations(polys: &[Vec<i64>], opts: &RecurrenceOptions) -> usize {
    let rational_polys = i64_polys_to_rational(polys);
    count_equations_rational(&rational_polys, opts)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_count_unknowns() {
        let opts = RecurrenceOptions {
            var_deg: 1,
            idx_deg: 1,
            diff_deg: 0,
            rec_len: 2,
            homogeneous: true,
            inhomo_var_deg: 0,
            inhomo_idx_deg: 0,
            denom_var_deg: 0,
            denom_idx_deg: 0,
            alternating_sign: false,
            modular_prefilter: false,
        };
        // vars_per_coeff = 2*2 = 4, num_coeff_vars = 2*1*4 = 8, denom = 0
        assert_eq!(count_unknowns(&opts), 8);
    }

    #[test]
    fn test_count_equations() {
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
            inhomo_var_deg: 0,
            inhomo_idx_deg: 0,
            denom_var_deg: 0,
            denom_idx_deg: 0,
            alternating_sign: false,
            modular_prefilter: false,
        };
        // m=7, rec_len=2, num_nn=5, max_poly_deg=0, max_j=0, eqs_per_nn=1
        assert_eq!(count_equations(&polys, &opts), 5);
    }
}
