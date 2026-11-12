//! Exact dominant-root forms for factor-gauge certificates.

use super::*;

// ---------------------------------------------------------------------------
// Nonzero rows by a dominant characteristic root
// ---------------------------------------------------------------------------

/// Largest number of indices past the rule start searched for the dominance
/// threshold of a [`FactorGaugeDominantRoot`].
const MAX_DOMINANCE_SEARCH: usize = 512;

/// A coefficient sequence `s_k(n) = [x^k] S_n = sum_i c_i lambda_i^n` of a
/// constant-coefficient scalar recurrence with distinct rational
/// characteristic roots, nonzero for every `n >= N` because its first root
/// dominates (part of a [`FactorGaugeCertificate`]).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FactorGaugeDominantRoot {
    /// `k`: the sequence is the coefficient of `x^k` in `S_n`.
    pub coefficient: usize,
    /// The roots `lambda_i` with a nonzero weight, distinct exact rationals
    /// as strings; the first has the strictly largest absolute value.
    pub roots: Vec<String>,
    /// The nonzero weights `c_i`, exact rationals as strings.
    pub weights: Vec<String>,
    /// `n*`: for every `n >= n*`, `|c_1| |lambda_1|^n > sum_(i >= 2) |c_i|
    /// |lambda_i|^n` (checked at `n*`; the right side decreases relative to
    /// the left); the values at `N <= n < n*` are checked to be nonzero.
    pub dominance_from: usize,
}

impl FactorGaugeDominantRoot {
    fn parsed(&self) -> Result<(Vec<Q>, Vec<Q>), RecurrenceError> {
        Ok((
            parse_polynomial(&self.roots)?,
            parse_polynomial(&self.weights)?,
        ))
    }

    pub(super) fn text(&self) -> String {
        self.roots
            .iter()
            .zip(&self.weights)
            .map(|(root, weight)| format!("({weight}) ({root})^n"))
            .collect::<Vec<_>>()
            .join(" + ")
    }

    /// Replay on the gauged rule: constant scalar multipliers, every root a
    /// root of the characteristic polynomial, the weights fit the
    /// window, the first root strictly dominates, the threshold inequality
    /// at `n*` and nonzero values before it.
    pub(super) fn verify(
        &self,
        gauged: &GaugedRule,
        window: &[Vec<Q>],
        window_start: usize,
        n_first: usize,
    ) -> Result<(), RecurrenceError> {
        let fail = |what: &str| {
            Err(RecurrenceError::VerificationFailed(format!(
                "factor gauge, dominant root: {what}"
            )))
        };
        let (roots, weights) = self.parsed()?;
        let Some(characteristic) = constant_characteristic_polynomial(gauged) else {
            return fail("the scalar multipliers are not constant");
        };
        if roots.is_empty()
            || roots.len() != weights.len()
            || weights.iter().any(Zero::is_zero)
            || roots
                .iter()
                .any(|root| !vec_poly::evaluate(&characteristic, root).is_zero())
        {
            return fail("a root is not a characteristic root or a weight is zero");
        }
        let distinct: BTreeSet<&Q> = roots.iter().collect();
        if distinct.len() != roots.len()
            || roots[1..].iter().any(|root| root.abs() >= roots[0].abs())
        {
            return fail("the roots are not distinct or the first does not dominate");
        }
        for (t, row) in window.iter().enumerate() {
            if exponential_value(&roots, &weights, window_start + t)
                != vec_poly::coeff(row, self.coefficient)
            {
                return fail("the weights do not fit the window");
            }
        }
        if self.dominance_from < n_first || !dominates(&roots, &weights, self.dominance_from) {
            return fail("the dominance inequality fails at the recorded index");
        }
        if (n_first..self.dominance_from).any(|n| exponential_value(&roots, &weights, n).is_zero())
        {
            return fail("a value before the dominance index is zero");
        }
        Ok(())
    }
}

/// `chi(z) = z^L - sum_j r_j z^(L-j)` (ascending) when every scalar
/// multiplier `kappa_j(n) / f(n) = r_j` is constant.
fn constant_characteristic_polynomial(gauged: &GaugedRule) -> Option<Vec<Q>> {
    let f_leading = gauged.f.last()?.clone();
    let order = gauged.order();
    let mut characteristic = vec![Q::zero(); order + 1];
    characteristic[order] = Q::one();
    for (lag, kappa) in (1usize..).zip(&gauged.kappa) {
        let kappa = vec_poly::trim(kappa);
        let ratio = vec_poly::coeff(&kappa, gauged.f.len() - 1) / &f_leading;
        if vec_poly::trim(&vec_poly::scale(&gauged.f, &ratio)) != kappa {
            return None;
        }
        characteristic[order - lag] = -ratio;
    }
    Some(characteristic)
}

fn exponential_value(roots: &[Q], weights: &[Q], n: usize) -> Q {
    let exponent = i32::try_from(n).unwrap_or(i32::MAX);
    roots
        .iter()
        .zip(weights)
        .fold(Q::zero(), |total, (root, weight)| {
            total + weight * root.pow(exponent)
        })
}

/// `|c_1| |lambda_1|^n > sum_(i >= 2) |c_i| |lambda_i|^n`, compared as
/// `sum_(i >= 2) |c_i| (|lambda_i| / |lambda_1|)^n < |c_1|`.
fn dominates(roots: &[Q], weights: &[Q], n: usize) -> bool {
    let exponent = i32::try_from(n).unwrap_or(i32::MAX);
    let lead = roots[0].abs();
    let rest = roots[1..]
        .iter()
        .zip(&weights[1..])
        .fold(Q::zero(), |total, (root, weight)| {
            total + weight.abs() * (root.abs() / &lead).pow(exponent)
        });
    rest < weights[0].abs()
}

/// Solve `matrix * x = rhs` exactly (square, Gaussian elimination); `None`
/// when singular.
fn solve_square_system(mut matrix: Vec<Vec<Q>>, mut rhs: Vec<Q>) -> Option<Vec<Q>> {
    let size = rhs.len();
    for column in 0..size {
        let pivot = (column..size).find(|&row| !matrix[row][column].is_zero())?;
        matrix.swap(column, pivot);
        rhs.swap(column, pivot);
        for row in 0..size {
            if row != column && !matrix[row][column].is_zero() {
                let factor = matrix[row][column].clone() / &matrix[column][column];
                let pivot_row = matrix[column].clone();
                for (entry, pivot_entry) in matrix[row].iter_mut().zip(&pivot_row).skip(column) {
                    *entry -= pivot_entry.clone() * &factor;
                }
                let value = rhs[column].clone() * &factor;
                rhs[row] -= value;
            }
        }
    }
    Some((0..size).map(|i| rhs[i].clone() / &matrix[i][i]).collect())
}

/// The dominant-root form of `s_k(n) = [x^k] S_n` (module documentation,
/// "Nonzero rows"), when it proves `s_k(n) != 0` for every `n >= n_first`.
pub(super) fn dominant_root_form(
    gauged: &GaugedRule,
    window: &[Vec<Q>],
    window_start: usize,
    k: usize,
    n_first: usize,
) -> Option<FactorGaugeDominantRoot> {
    let characteristic = constant_characteristic_polynomial(gauged)?;
    let factorization =
        crate::factorization::irreducible_factorization_rational_coeffs(&characteristic).ok()?;
    if factorization
        .factors
        .iter()
        .any(|(factor, multiplicity)| factor.len() != 2 || *multiplicity != 1)
    {
        return None;
    }
    let roots: Vec<Q> = factorization
        .factors
        .iter()
        .map(|(factor, _)| Q::new(-factor[0].clone(), factor[1].clone()))
        .collect();
    let matrix: Vec<Vec<Q>> = (0..window.len())
        .map(|t| {
            let exponent = i32::try_from(window_start + t).ok()?;
            Some(roots.iter().map(|root| root.pow(exponent)).collect())
        })
        .collect::<Option<_>>()?;
    let rhs: Vec<Q> = window.iter().map(|row| vec_poly::coeff(row, k)).collect();
    let weights = solve_square_system(matrix, rhs)?;
    let mut terms: Vec<(Q, Q)> = roots
        .into_iter()
        .zip(weights)
        .filter(|(_, weight)| !weight.is_zero())
        .collect();
    terms.sort_by(|a, b| b.0.abs().cmp(&a.0.abs()).then_with(|| b.0.cmp(&a.0)));
    if terms.is_empty()
        || terms
            .get(1)
            .is_some_and(|second| second.0.abs() == terms[0].0.abs())
    {
        return None;
    }
    let (roots, weights): (Vec<Q>, Vec<Q>) = terms.into_iter().unzip();
    let dominance_from =
        (n_first..=n_first + MAX_DOMINANCE_SEARCH).find(|&n| dominates(&roots, &weights, n))?;
    if (n_first..dominance_from).any(|n| exponential_value(&roots, &weights, n).is_zero()) {
        return None;
    }
    Some(FactorGaugeDominantRoot {
        coefficient: k,
        roots: roots.iter().map(format_rational_coeff).collect(),
        weights: weights.iter().map(format_rational_coeff).collect(),
        dominance_from,
    })
}
