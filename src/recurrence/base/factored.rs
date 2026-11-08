//! Exact factorization of bivariate recurrence coefficients over `Q` and the
//! factored LaTeX presentation built on it.
//!
//! A coefficient `f(n, t)` is factored by the existing exact routines of
//! [`crate::factorization`]: univariate coefficients and coefficients monic in
//! `n` go to them directly, every other coefficient is first sheared by the
//! invertible change of variables `n = u`, `t = v + c u` (with `c` chosen so the
//! top homogeneous part does not vanish at `(1, c)`), which makes it monic in
//! `u` up to a rational constant.  The factors are mapped back, normalized to
//! primitive integer polynomials with a canonical sign, and the product is
//! multiplied out and compared with the input before it is returned, so a
//! [`BivarFactorization`] is always exact.  An exhausted factorization budget or
//! a failed check is a [`BivarFactorizationError`]: the callers never present
//! an unfinished factorization as complete.
//!
//! Entry points: [`BivarPoly::irreducible_factorization`],
//! [`BivarPoly::irreducible_factorization_with_budget`],
//! [`BivarPoly::to_latex_factored`] and [`Recurrence::to_latex_factored`].

use super::export::{append_latex_summand, latex_summand};
use super::*;
use crate::factorization::{
    irreducible_factorization_monic_bivariate_with_budget,
    irreducible_factorization_rational_coeffs_with_budget, PolynomialFactorizationError,
    DEFAULT_RECOMBINATION_SUBSET_LIMIT,
};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::fmt;

/// An exact factorization of a nonzero bivariate polynomial over `Q`:
/// `input = constant * prod factor^multiplicity`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BivarFactorization {
    /// The nonzero rational constant of the factorization.
    pub constant: BigRational,
    /// The irreducible non-constant factors with their multiplicities.  Each
    /// factor is a primitive integer polynomial in `n` and `t` whose first
    /// displayed term (ascending powers of `t`, then descending powers of `n`)
    /// has a positive coefficient.  Factors in `n` alone come first, then
    /// factors in `t` alone, then mixed factors, each group by total degree.
    pub factors: Vec<(BivarPoly, usize)>,
}

impl BivarFactorization {
    /// Multiply the factorization out again.
    pub fn expand(&self) -> BivarPoly {
        let mut product = vec![vec![self.constant.clone()]];
        for (factor, multiplicity) in &self.factors {
            let dense = dense_of(factor);
            for _ in 0..*multiplicity {
                product = dense_multiply(&product, &dense);
            }
        }
        poly_of(product)
    }
}

/// Why a bivariate coefficient was not factored.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BivarFactorizationError {
    /// The polynomial is zero.
    ZeroPolynomial,
    /// The exact factorization routine stopped, for example because the
    /// recombination budget was exhausted
    /// ([`PolynomialFactorizationError::RecombinationBudgetExceeded`]).
    Incomplete(PolynomialFactorizationError),
    /// The multiplied-out factorization differs from the input.  This is an
    /// internal consistency failure; the factorization is discarded.
    VerificationFailed,
}

impl BivarFactorizationError {
    /// True when the factorization budget ran out (as opposed to an internal
    /// failure).
    pub fn is_budget_exhausted(&self) -> bool {
        matches!(
            self,
            Self::Incomplete(PolynomialFactorizationError::RecombinationBudgetExceeded { .. })
        )
    }
}

impl fmt::Display for BivarFactorizationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroPolynomial => write!(formatter, "the polynomial is zero"),
            Self::Incomplete(error) => write!(formatter, "factorization incomplete: {error}"),
            Self::VerificationFailed => {
                write!(formatter, "the factorization failed its exact re-expansion")
            }
        }
    }
}

impl std::error::Error for BivarFactorizationError {}

type Dense = Vec<Vec<BigRational>>;

fn dense_of(poly: &BivarPoly) -> Dense {
    let mut dense = poly.coeffs.clone();
    dense_trim(&mut dense);
    dense
}

fn poly_of(mut dense: Dense) -> BivarPoly {
    dense_trim(&mut dense);
    if dense.is_empty() {
        dense.push(vec![BigRational::zero()]);
    }
    BivarPoly::from_coefficients_normalized(dense)
}

/// Remove trailing zeros of every row and trailing empty rows.
fn dense_trim(dense: &mut Dense) {
    for row in dense.iter_mut() {
        while row.last().is_some_and(Zero::is_zero) {
            row.pop();
        }
    }
    while dense.last().is_some_and(Vec::is_empty) {
        dense.pop();
    }
}

fn dense_multiply(left: &Dense, right: &Dense) -> Dense {
    let mut product: Dense = Vec::new();
    for (i, left_row) in left.iter().enumerate() {
        for (j, left_entry) in left_row.iter().enumerate() {
            if left_entry.is_zero() {
                continue;
            }
            for (k, right_row) in right.iter().enumerate() {
                for (l, right_entry) in right_row.iter().enumerate() {
                    if right_entry.is_zero() {
                        continue;
                    }
                    if product.len() <= i + k {
                        product.resize(i + k + 1, Vec::new());
                    }
                    let row = &mut product[i + k];
                    if row.len() <= j + l {
                        row.resize(j + l + 1, BigRational::zero());
                    }
                    row[j + l] += left_entry * right_entry;
                }
            }
        }
    }
    dense_trim(&mut product);
    product
}

/// Substitute `t -> t + shear * n`, so that `n^i t^j` becomes `n^i (t + shear n)^j`.
fn dense_shear(dense: &Dense, shear: &BigRational) -> Dense {
    let mut result: Dense = Vec::new();
    for (i, row) in dense.iter().enumerate() {
        for (j, entry) in row.iter().enumerate() {
            if entry.is_zero() {
                continue;
            }
            // (t + shear n)^j = sum_m C(j, m) shear^(j-m) n^(j-m) t^m.
            let mut binomial = BigInt::one();
            for m in 0..=j {
                let power = shear.pow(i32::try_from(j - m).expect("degree fits i32"));
                let value = entry * &power * BigRational::from_integer(binomial.clone());
                if !value.is_zero() {
                    let n_degree = i + (j - m);
                    if result.len() <= n_degree {
                        result.resize(n_degree + 1, Vec::new());
                    }
                    if result[n_degree].len() <= m {
                        result[n_degree].resize(m + 1, BigRational::zero());
                    }
                    result[n_degree][m] += value;
                }
                binomial = binomial * BigInt::from(j - m) / BigInt::from(m + 1);
            }
        }
    }
    dense_trim(&mut result);
    result
}

fn total_degree(dense: &Dense) -> usize {
    dense
        .iter()
        .enumerate()
        .flat_map(|(i, row)| {
            row.iter()
                .enumerate()
                .filter(|(_, entry)| !entry.is_zero())
                .map(move |(j, _)| i + j)
        })
        .max()
        .unwrap_or(0)
}

/// The first nonzero entry in display order (ascending `t` power, then
/// descending `n` power).
fn first_displayed_entry(dense: &Dense) -> Option<&BigRational> {
    let width = dense.iter().map(Vec::len).max().unwrap_or(0);
    for j in 0..width {
        for row in dense.iter().rev() {
            if let Some(entry) = row.get(j).filter(|entry| !entry.is_zero()) {
                return Some(entry);
            }
        }
    }
    None
}

fn scale_dense(dense: &mut Dense, scalar: &BigRational) {
    for row in dense.iter_mut() {
        for entry in row.iter_mut() {
            *entry = &*entry * scalar;
        }
    }
}

/// Return `(scalar, primitive)` with `dense = scalar * primitive`, where
/// `primitive` has coprime integer entries and a positive first displayed term.
fn primitive_normalization(dense: &Dense) -> (BigRational, Dense) {
    let mut denominator_lcm = BigInt::one();
    for entry in dense.iter().flatten() {
        denominator_lcm = denominator_lcm.lcm(entry.denom());
    }
    let mut numerator_gcd = BigInt::zero();
    for entry in dense.iter().flatten() {
        let scaled = entry * BigRational::from_integer(denominator_lcm.clone());
        numerator_gcd = numerator_gcd.gcd(scaled.numer());
    }
    let mut divisor = BigRational::new(numerator_gcd, denominator_lcm);
    let mut primitive = dense.clone();
    if first_displayed_entry(dense).is_some_and(Signed::is_negative) {
        divisor = -divisor;
    }
    scale_dense(&mut primitive, &divisor.recip());
    (divisor, primitive)
}

fn factor_class(dense: &Dense) -> u8 {
    let depends_on_n = dense.len() > 1;
    let depends_on_t = dense.iter().any(|row| row.len() > 1);
    match (depends_on_n, depends_on_t) {
        (_, false) => 0,
        (false, true) => 1,
        (true, true) => 2,
    }
}

fn factor_sort_key(dense: &Dense) -> (u8, usize, Vec<BigRational>) {
    let flat: Vec<BigRational> = dense.iter().flatten().cloned().collect();
    (factor_class(dense), total_degree(dense), flat)
}

fn univariate_factors(
    dense: &Dense,
    in_n: bool,
    budget: usize,
) -> Result<(BigRational, Vec<(Dense, usize)>), BivarFactorizationError> {
    let coefficients: Vec<BigRational> = if in_n {
        dense
            .iter()
            .map(|row| row.first().cloned().unwrap_or_else(BigRational::zero))
            .collect()
    } else {
        dense.first().cloned().unwrap_or_default()
    };
    let factorization =
        irreducible_factorization_rational_coeffs_with_budget(&coefficients, budget)
            .map_err(BivarFactorizationError::Incomplete)?;
    let factors = factorization
        .factors
        .into_iter()
        .map(|(integers, multiplicity)| {
            let entries: Vec<BigRational> = integers
                .into_iter()
                .map(BigRational::from_integer)
                .collect();
            let factor = if in_n {
                entries.into_iter().map(|entry| vec![entry]).collect()
            } else {
                vec![entries]
            };
            (factor, multiplicity)
        })
        .collect();
    Ok((factorization.constant, factors))
}

fn monic_bivariate_factors(
    monic: &Dense,
    budget: usize,
) -> Result<Vec<(Dense, usize)>, BivarFactorizationError> {
    let factors = irreducible_factorization_monic_bivariate_with_budget(monic, budget)
        .map_err(BivarFactorizationError::Incomplete)?;
    Ok(factors
        .into_iter()
        .map(|factor| (factor.coefficients, factor.multiplicity))
        .collect())
}

fn factor_dense(
    dense: &Dense,
    budget: usize,
) -> Result<(BigRational, Vec<(Dense, usize)>), BivarFactorizationError> {
    let n_degree = dense.len() - 1;
    let t_degree = dense.iter().map(Vec::len).max().unwrap_or(1) - 1;
    if t_degree == 0 {
        return univariate_factors(dense, true, budget);
    }
    if n_degree == 0 {
        return univariate_factors(dense, false, budget);
    }
    let top_row = &dense[n_degree];
    if top_row.len() == 1 {
        // Already monic in n up to the constant `top_row[0]`.
        let constant = top_row[0].clone();
        let mut monic = dense.clone();
        scale_dense(&mut monic, &constant.recip());
        return Ok((constant, monic_bivariate_factors(&monic, budget)?));
    }
    let degree = total_degree(dense);
    // The top homogeneous part evaluated at (1, c) for c = 0, 1, ..., degree.
    for shear in 0..=degree {
        let shear = BigRational::from_integer(BigInt::from(shear));
        let mut leading = BigRational::zero();
        for (i, row) in dense.iter().enumerate() {
            for (j, entry) in row.iter().enumerate() {
                if i + j == degree {
                    leading += entry * shear.pow(i32::try_from(j).expect("degree fits i32"));
                }
            }
        }
        if leading.is_zero() {
            continue;
        }
        // Substitute t = v + shear * n.
        let mut sheared = dense_shear(dense, &shear);
        scale_dense(&mut sheared, &leading.recip());
        if sheared.len() != degree + 1 || sheared[degree] != vec![BigRational::one()] {
            return Err(BivarFactorizationError::VerificationFailed);
        }
        let factors = monic_bivariate_factors(&sheared, budget)?;
        // Substitute back: v = t - shear * n.
        let back = -shear;
        let factors = factors
            .into_iter()
            .map(|(factor, multiplicity)| (dense_shear(&factor, &back), multiplicity))
            .collect();
        return Ok((leading, factors));
    }
    Err(BivarFactorizationError::VerificationFailed)
}

impl BivarPoly {
    /// Factor this polynomial in `n` and `t` over `Q`, with the default
    /// recombination budget.  See [`Self::irreducible_factorization_with_budget`].
    pub fn irreducible_factorization(&self) -> Result<BivarFactorization, BivarFactorizationError> {
        self.irreducible_factorization_with_budget(DEFAULT_RECOMBINATION_SUBSET_LIMIT)
    }

    /// Factor this polynomial in `n` and `t` over `Q`, trying at most
    /// `max_recombination_subsets` factor subsets in each recombination.
    ///
    /// The result is verified: the product of the returned constant and
    /// factors equals the input exactly.  A budget that runs out is
    /// [`BivarFactorizationError::Incomplete`], never a partial factorization.
    /// Zero returns [`BivarFactorizationError::ZeroPolynomial`]; a nonzero
    /// constant has no factors.
    pub fn irreducible_factorization_with_budget(
        &self,
        max_recombination_subsets: usize,
    ) -> Result<BivarFactorization, BivarFactorizationError> {
        let dense = dense_of(self);
        if dense.is_empty() {
            return Err(BivarFactorizationError::ZeroPolynomial);
        }
        let (mut constant, raw_factors) = if dense.len() == 1 && dense[0].len() == 1 {
            (dense[0][0].clone(), Vec::new())
        } else {
            factor_dense(&dense, max_recombination_subsets)?
        };
        let mut normalized: Vec<(Dense, usize)> = Vec::new();
        for (factor, multiplicity) in raw_factors {
            let mut factor = factor;
            dense_trim(&mut factor);
            let (scalar, primitive) = primitive_normalization(&factor);
            constant *= scalar.pow(i32::try_from(multiplicity).expect("multiplicity fits i32"));
            normalized.push((primitive, multiplicity));
        }
        normalized.sort_by_key(|(factor, _)| factor_sort_key(factor));
        let factorization = BivarFactorization {
            constant,
            factors: normalized
                .into_iter()
                .map(|(factor, multiplicity)| (poly_of(factor), multiplicity))
                .collect(),
        };
        if dense_of(&factorization.expand()) == dense {
            Ok(factorization)
        } else {
            Err(BivarFactorizationError::VerificationFailed)
        }
    }

    /// Format as LaTeX with the polynomial factored over `Q`, or `None` when
    /// the factorization did not finish (see
    /// [`Self::irreducible_factorization_with_budget`]).
    ///
    /// Returns the sign and the unsigned text, so that callers can place the
    /// sign between summands.  A factor with more than one term is
    /// parenthesized and a repeated factor carries its exponent.  A rational
    /// constant is written first, `\frac{1}{2}`; a constant `1` is omitted.
    /// The text of the constant polynomial `1` is `"1"`.
    pub(super) fn factored_latex_parts(
        &self,
        max_recombination_subsets: usize,
    ) -> Result<(bool, String), BivarFactorizationError> {
        let factorization =
            self.irreducible_factorization_with_budget(max_recombination_subsets)?;
        let negative = factorization.constant.is_negative();
        let magnitude = factorization.constant.abs();
        let mut text = String::new();
        if !magnitude.is_one() || factorization.factors.is_empty() {
            text.push_str(&latex_rational(&magnitude));
        }
        for (factor, multiplicity) in &factorization.factors {
            let body = factor.to_latex();
            if factor.num_terms() > 1 {
                text.push_str(&format!("({body})"));
            } else {
                text.push_str(&body);
            }
            if *multiplicity > 1 {
                text.push_str(&format!("^{{{multiplicity}}}"));
            }
        }
        Ok((negative, text))
    }

    /// Format as LaTeX with the polynomial factored over `Q`; see
    /// [`Self::irreducible_factorization_with_budget`].  Falls back to the
    /// expanded [`Self::to_latex`] when the factorization does not finish.
    pub fn to_latex_factored(&self) -> String {
        match self.factored_latex_parts(DEFAULT_RECOMBINATION_SUBSET_LIMIT) {
            Ok((negative, text)) => format!("{}{text}", if negative { "-" } else { "" }),
            Err(_) => self.to_latex(),
        }
    }
}

fn latex_rational(value: &BigRational) -> String {
    if value.denom().is_one() {
        value.numer().to_string()
    } else {
        format!("\\frac{{{}}}{{{}}}", value.numer(), value.denom())
    }
}

/// A recurrence rendered with factored coefficients, with the coefficients
/// that could not be factored listed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FactoredRecurrenceLatex {
    /// The formula `lhs = rhs` without dollar signs or range.
    pub formula: String,
    /// One note per coefficient that stays expanded because its
    /// factorization did not finish (for instance because the budget ran
    /// out), naming the polynomial reference it multiplies.  Empty when every
    /// coefficient is factored.
    pub unfactored_notes: Vec<String>,
}

impl Recurrence {
    /// Format the recurrence as LaTeX with every coefficient factored over
    /// `Q`, for example
    /// `(n + 2) P_{n} = (2n + 1)(1 + t) P_{n-1} - (n - 1)(1 - t)^{2} P_{n-2}`.
    ///
    /// A coefficient whose factorization does not finish within
    /// `max_recombination_subsets` stays expanded and is reported in
    /// [`FactoredRecurrenceLatex::unfactored_notes`]; it is never presented as
    /// factored.  Compare [`Recurrence::to_latex`], which never factors.
    pub fn to_latex_factored_with_budget(
        &self,
        max_recombination_subsets: usize,
    ) -> FactoredRecurrenceLatex {
        let mut notes = Vec::new();
        let mut summand =
            |coefficient: &BivarPoly, reference: Option<&str>, label: &str| -> (bool, String) {
                match coefficient.factored_latex_parts(max_recombination_subsets) {
                    Ok((negative, magnitude)) => {
                        let text = match reference {
                            Some(reference) if magnitude == "1" => reference.to_string(),
                            Some(reference) => format!("{magnitude} {reference}"),
                            None => magnitude,
                        };
                        (negative, text)
                    }
                    Err(error) => {
                        notes.push(format!("coefficient of {label} left expanded: {error}"));
                        latex_summand(coefficient, reference)
                    }
                }
            };
        let lhs = match &self.denominator {
            Some(coefficient) => {
                let (negative, product) = summand(coefficient, Some("P_{n}"), "P_{n} (left side)");
                format!("{}{product}", if negative { "-" } else { "" })
            }
            None => "P_{n}".to_string(),
        };
        let mut rhs = String::new();
        for term in &self.terms {
            if term.coeff.is_zero() {
                continue;
            }
            let reference = fmt_poly_ref_latex_with_sign(term.sign, term.offset, term.deriv_order);
            let parts = summand(&term.coeff, Some(&reference), &reference);
            append_latex_summand(&mut rhs, parts);
        }
        if let Some(forcing) = self.inhomogeneous.as_ref().filter(|p| !p.is_zero()) {
            let parts = summand(forcing, None, "the inhomogeneous term");
            append_latex_summand(&mut rhs, parts);
        }
        if rhs.is_empty() {
            rhs.push('0');
        }
        FactoredRecurrenceLatex {
            formula: format!("{lhs} = {rhs}"),
            unfactored_notes: notes,
        }
    }

    /// [`Self::to_latex_factored_with_budget`] with the default budget,
    /// returning only the formula.
    pub fn to_latex_factored(&self) -> String {
        self.to_latex_factored_with_budget(DEFAULT_RECOMBINATION_SUBSET_LIMIT)
            .formula
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::test_support::*;

    fn poly(rows: &[&[i64]]) -> BivarPoly {
        BivarPoly {
            coeffs: rows
                .iter()
                .map(|row| row.iter().map(|&entry| br(entry, 1)).collect())
                .collect(),
        }
    }

    fn factored(rows: &[&[i64]]) -> String {
        let input = poly(rows);
        let factorization = input.irreducible_factorization().unwrap();
        assert_eq!(dense_of(&factorization.expand()), dense_of(&input));
        input.to_latex_factored()
    }

    #[test]
    fn constants_and_zero() {
        assert_eq!(factored(&[&[6]]), "6");
        assert_eq!(factored(&[&[-1]]), "-1");
        assert_eq!(poly(&[&[0]]).to_latex_factored(), "0");
        assert_eq!(
            poly(&[&[0]]).irreducible_factorization(),
            Err(BivarFactorizationError::ZeroPolynomial)
        );
        let half = BivarPoly {
            coeffs: vec![vec![br(-1, 2)]],
        };
        assert_eq!(half.to_latex_factored(), "-\\frac{1}{2}");
    }

    #[test]
    fn univariate_shapes_and_sign_normalization() {
        // 2n + 1; n - 1 stays (n - 1) while 1 - n is written -(n - 1).
        assert_eq!(factored(&[&[1], &[2]]), "(2n + 1)");
        assert_eq!(factored(&[&[-1], &[1]]), "(n - 1)");
        assert_eq!(factored(&[&[1], &[-1]]), "-(n - 1)");
        // (1 - t)^2 and t - 2 = -(2 - t) in t alone.
        assert_eq!(factored(&[&[1, -2, 1]]), "(1 - t)^{2}");
        assert_eq!(factored(&[&[-2, 1]]), "-(2 - t)");
        assert_eq!(factored(&[&[0, 2]]), "2t");
        // n^2 - 1 = (n - 1)(n + 1): factors sorted by coefficients.
        assert_eq!(factored(&[&[-1], &[0], &[1]]), "(n - 1)(n + 1)");
    }

    #[test]
    fn products_of_n_and_t_factors_need_the_shear() {
        // (n + 2)(1 + t) = n + 2 + nt + 2t: top n-row [1, 1] is not constant.
        let input = poly(&[&[2, 2], &[1, 1]]);
        assert_eq!(input.to_latex_factored(), "(n + 2)(1 + t)");
        // (2n + 1)(1 + t)(1 - t)^2.
        let first = dense_multiply(
            &dense_of(&poly(&[&[1], &[2]])),
            &dense_of(&poly(&[&[1, 1]])),
        );
        let product = dense_multiply(&first, &dense_of(&poly(&[&[1, -2, 1]])));
        let input = poly_of(product);
        assert_eq!(input.to_latex_factored(), "(2n + 1)(1 - t)^{2}(1 + t)");
        // An irreducible mixed factor: n - t (monic in n) and n + t^2.
        assert_eq!(factored(&[&[0, -1], &[1]]), "(n - t)");
        assert_eq!(factored(&[&[0, 0, 1], &[1]]), "(n + t^{2})");
        // Mixed factor not monic in n: (2n - t)(n + 3t + 1).
        let product = dense_multiply(
            &dense_of(&poly(&[&[0, -1], &[2]])),
            &dense_of(&poly(&[&[1, 3], &[1]])),
        );
        let input = poly_of(product);
        let factorization = input.irreducible_factorization().unwrap();
        assert_eq!(factorization.factors.len(), 2);
        assert_eq!(dense_of(&factorization.expand()), dense_of(&input));
        // Repeated mixed factor.
        let square = dense_multiply(
            &dense_of(&poly(&[&[1, 1], &[2]])),
            &dense_of(&poly(&[&[1, 1], &[2]])),
        );
        let factorization = poly_of(square).irreducible_factorization().unwrap();
        assert_eq!(factorization.factors.len(), 1);
        assert_eq!(factorization.factors[0].1, 2);
    }

    #[test]
    fn irreducible_mixed_polynomial_stays_whole_with_content() {
        // 2(n^2 + t^3 + 1)-style content is moved into the constant.
        let input = poly(&[&[2, 0, 0, 2], &[0], &[2]]);
        let factorization = input.irreducible_factorization().unwrap();
        assert_eq!(factorization.constant, br(2, 1));
        assert_eq!(factorization.factors.len(), 1);
        assert_eq!(input.to_latex_factored(), "2(n^{2} + 1 + t^{3})");
    }

    #[test]
    fn exhausted_budget_is_reported_not_hidden() {
        // x^4 + 1 over Q is irreducible but splits modulo every prime, so
        // recombination needs trials; a zero budget must be an error.
        let input = poly(&[&[1, 0, 0, 0, 1]]);
        let error = input
            .irreducible_factorization_with_budget(0)
            .expect_err("zero budget");
        assert!(error.is_budget_exhausted(), "{error}");
        assert!(input.irreducible_factorization().is_ok());
    }

    #[test]
    fn recurrence_is_rendered_with_factored_coefficients() {
        let recurrence = Recurrence {
            denominator: Some(poly(&[&[2], &[1]])),
            terms: vec![
                RecurrenceTerm::new(1, 0, poly(&[&[1, 1], &[2, 2]])),
                RecurrenceTerm::new(2, 0, poly(&[&[1, -2, 1], &[-1, 2, -1]])),
            ],
            inhomogeneous: None,
        };
        let report = recurrence.to_latex_factored_with_budget(1 << 10);
        assert!(report.unfactored_notes.is_empty());
        assert_eq!(
            report.formula,
            "(n + 2) P_{n} = (2n + 1)(1 + t) P_{n-1} - (n - 1)(1 - t)^{2} P_{n-2}"
        );
        // The expanded rendering keeps its old shape.
        assert!(recurrence.to_latex().contains("\\bigl("));
    }

    #[test]
    fn unfactored_coefficients_are_listed() {
        let recurrence = Recurrence {
            denominator: None,
            terms: vec![RecurrenceTerm::new(1, 0, poly(&[&[1, 0, 0, 0, 1]]))],
            inhomogeneous: None,
        };
        let report = recurrence.to_latex_factored_with_budget(0);
        assert_eq!(report.unfactored_notes.len(), 1);
        assert!(report.unfactored_notes[0].contains("P_{n-1}"));
        assert!(report.formula.contains("t^{4}"));
    }
}
