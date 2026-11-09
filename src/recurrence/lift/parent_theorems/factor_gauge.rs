//! Factor gauge with a bounded residual: a fixed factor `q(x)` carries the
//! growth of the rows, and what remains lies in a fixed finite-dimensional
//! space whose real-rootedness is decided for every index.  This is the
//! general, exact form of the RealRooted row routes `powerForm`,
//! `closedForm` and `factorPower`.
//!
//! # The theorem
//!
//! Let the rule be
//!
//! ```text
//! f(n) P_n = sum_(j=1)^L kappa_j(n) q(x)^j P_(n-j),        n >= N = n_0 + K,
//! ```
//!
//! with a monic `q` in `Q[x]` (a nonzero constant is allowed), scalars
//! `kappa_j`, `f` in `Q[n]`, order `L >= 2`, first index `n_0` and `K >= L`
//! seed rows (no derivative, no `(-1)^n`, no inhomogeneous term, a
//! left-hand factor `x^e f(n)` normalized as in the `parent_theorems`
//! documentation).  The rule is matched exactly: every lag coefficient
//! `c_j(n, x)` must factor as `gamma_j(n) w_j(x)`, `q` is the monic
//! `j_0`-th root of `w_(j_0)` for the smallest lag `j_0` (computed
//! coefficient by coefficient and checked by raising it to the power), and
//! `w_j / lead(w_j) = q^j` for every lag.
//!
//! **The gauge.**  Put `e = N - 1` (the last seed row) and, for the window
//! `m = N - L, ..., e`, `R_m = q^(e - m) P_m`; let `g` be the monic gcd of
//! these polynomials and `T_m = R_m / g`.  Define `S_n` for `n >= N - L` by
//! `S_m = T_m` on the window and
//!
//! ```text
//! f(n) S_n = sum_j kappa_j(n) S_(n-j),        n >= N.
//! ```
//!
//! Then `P_n = q^(n - e) g S_n` for every `n >= e`: for the window this is
//! the definition, and for `n >= N` it follows by induction, since
//! `f(n) q^(n-e) g S_n = sum_j kappa_j(n) q^j q^(n-j-e) g S_(n-j)`.  (For
//! `n - j < e` the factor `q^(n-j-e)` is a negative power, but the identity
//! holds in `Q(x)` and the left side is a polynomial.)  So every `S_n`
//! lies in the fixed space `W` spanned by the `T_m`, and every linear
//! functional of `S_n` (for example the coefficient of `x^k`) is a scalar
//! sequence with the same scalar recurrence.  This is the gauge of
//! RealRooted `forall_eq_C_mul_pow_of_rec2`/`_rec3` (`W` one-dimensional)
//! and `eq_pow_mul_of_factoredThreeTerm` (route `factorPower`).
//!
//! **Real-rootedness for every index.**  `P_n` (`n >= N`) is nonzero and
//! real-rooted when `q` and `g` are real-rooted (checked exactly) and
//! `S_n` is nonzero and real-rooted.  The last is decided from the shape of
//! `W` (the largest degree `D` of its elements):
//!
//! - `D <= 1` (`constant`, `linear`): every nonzero element is real-rooted;
//! - `D = 2` (`quadratic`): with `S_n = A x^2 + B x + C`, `B^2 - 4 A C >= 0`
//!   suffices, and is necessary when `A != 0` (for `A = 0` it reads
//!   `B^2 >= 0`, and `S_n` has degree at most one);
//! - `D = 3` (`cubic`): the formal discriminant of `a x^3 + b x^2 + c x + d`,
//!   `18abcd - 4b^3 d + b^2 c^2 - 4ac^3 - 27a^2 d^2 >= 0`, suffices: for
//!   `a != 0` it is the classical criterion, for `a = 0` it equals
//!   `b^2 (c^2 - 4bd)`, and for `a = b = 0` the polynomial is linear;
//! - `D = 4` with every element palindromic, `S_n = a x^4 + b x^3 + c x^2 +
//!   b x + a` (`palindromic_quartic`): with `a(n)` of one strict sign
//!   (normalized to `a > 0`) and `u = x + 1/x`, `x^(-2) S_n = h(u) = a u^2
//!   + b u + (c - 2a)`, and `S_n = a (x^2 - u_1 x + 1)(x^2 - u_2 x + 1)`
//!   for the roots `u_1, u_2` of `h`; each factor has real roots exactly
//!   when `|u_i| >= 2`.  Three sufficient cases, each a conjunction of signs:
//!   both roots `<= -2` (`b >= 4a`, `h(-2) = 2a - 2b + c >= 0`, discriminant
//!   `b^2 - 4a(c - 2a) >= 0`; RealRooted `splits_palindromic_quartic`), both
//!   `>= 2` (the mirror `x -> -x`: `-b >= 4a`, `h(2) = 2a + 2b + c >= 0`,
//!   the discriminant), or `h(2) <= 0` and `h(-2) <= 0` (then `u_1 <= -2`
//!   and `u_2 >= 2`).  The first case whose conditions all hold is used.
//!
//! For `D >= 2` the conditions are signs of polynomials in `n`: the
//! coefficient sequences `s_k(n) = [x^k] S_n` must be polynomials `pi_k(n)`.
//! Each `pi_k` is the interpolant of degree at most a bound (at most
//! `MAX_COORDINATE_DEGREE`) through computed values `s_k(N - L), ...`, and
//! it is accepted only when it satisfies the scalar recurrence identically,
//! `f(n) pi_k(n) - sum_j kappa_j(n) pi_k(n - j) = 0` in `Q[n]`, and agrees
//! with the `L` window values; then `s_k = pi_k` for every `n >= N - L` by
//! induction (nothing is fitted: the identity is the proof).  For a
//! constant-coefficient rule this happens exactly when the characteristic
//! roots of the scalar recurrence are all `1`.  When no such polynomial is
//! found for `D >= 2`, the theorem is `not_applicable` (no approximation).
//! The signs are decided for every integer `n >= N` by
//! [`crate::decide_index_rational_function_sign`].
//!
//! **Nonzero rows.**  `S_n != 0` for `n >= N` comes from the coefficient
//! `pi_D(n) != 0` (decided), or else from the shared leading-coefficient
//! certificate of the family ([`certify_leading_coefficient`], proved from
//! `N` on), as for coordinates that grow exponentially.
//!
//! **Seed rows** `P_(n_0), ..., P_e` are checked exactly (nonzero,
//! real-rooted).  The conclusion: every `P_n`, `n >= n_0`, is nonzero and
//! real-rooted.
//!
//! **Root region.**  Claimed only when proved: `(-inf, 0)` when `q`, `g` and
//! every seed row have all coefficients of one strict sign (for a
//! real-rooted polynomial: every root negative) and so does `S_n` for every
//! `n >= N` (`D = 0`, or every `pi_k(n) > 0`, or every `pi_k(n) < 0`,
//! decided); `(-inf, 0]` with "nonzero coefficients of one sign" instead.
//!
//! **Interlacing** is claimed only for `D = 0` and `deg q = 1`: then
//! `P_n = s(n) g q^(n - e)` with scalars `s(n) != 0`, so
//! `P_n = (s(n) / s(n-1)) q P_(n-1)` and `P_(n-1)` interlaces `P_n` (weakly)
//! for every `n >= N` (RealRooted `interlaces_of_forall_eq_C_mul_pow`).
//!
//! A first-order rule (`L = 1`) is the product route's case and is not
//! matched.  The certificate [`FactorGaugeCertificate`] records `q`, the
//! scalars, `g`, the residual window and the coefficient polynomials;
//! [`FactorGaugeCertificate::verify`] replays every identity on its own.

use super::*;
use crate::polynomial::Polynomial;

const REFERENCE: &str = "Factor gauge with a bounded residual: f(n) P_n = sum_j kappa_j(n) \
     q(x)^j P_(n-j) gives P_n = q^(n-e) g S_n with S_n in the fixed span of the gauged seed rows \
     and scalar coordinates; S_n real-rooted by its degree (at most 1), the discriminant \
     (degree 2 or 3) or, for palindromic quartics, u = x + 1/x; RealRooted routes powerForm, \
     closedForm and factorPower (RealRooted/ProductSequence/Interlacing.lean, \
     RealRooted/QuadraticRoot.lean, RealRooted/Recurrence/FactorPower.lean)";

/// Largest degree in `n` of an accepted coefficient polynomial `pi_k(n)`.
const MAX_COORDINATE_DEGREE: usize = 24;

/// The shape of the residual space `W` of a [`FactorGaugeCertificate`]: the
/// largest degree of its elements, and palindromy for degree four.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FactorGaugeResidualShape {
    /// Every element is a constant.
    Constant,
    /// Degree at most one.
    Linear,
    /// Degree at most two (decided by the discriminant).
    Quadratic,
    /// Degree at most three (decided by the formal cubic discriminant).
    Cubic,
    /// Every element is `a x^4 + b x^3 + c x^2 + b x + a` (decided through
    /// `u = x + 1/x`).
    PalindromicQuartic,
}

impl FactorGaugeResidualShape {
    /// The largest degree of an element of the residual space.
    pub fn degree(self) -> usize {
        match self {
            Self::Constant => 0,
            Self::Linear => 1,
            Self::Quadratic => 2,
            Self::Cubic => 3,
            Self::PalindromicQuartic => 4,
        }
    }

    fn text(self) -> &'static str {
        match self {
            Self::Constant => "constant",
            Self::Linear => "linear",
            Self::Quadratic => "quadratic",
            Self::Cubic => "cubic",
            Self::PalindromicQuartic => "palindromic quartic",
        }
    }
}

/// The exact data of a [`check_factor_gauge_theorem`] report (module
/// documentation of the `factor_gauge` submodule, summarized at
/// [`check_factor_gauge_theorem`]).  Polynomials are ascending coefficient
/// lists of exact rationals as strings: in `x` for `gauge_factor`,
/// `common_factor` and `residual_seeds`, in the index `n` for the others.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FactorGaugeCertificate {
    /// The monic gauge factor `q(x)` (`["1"]` when the multipliers are free
    /// of `x`).
    pub gauge_factor: Vec<String>,
    /// The left-hand factor `f(n)` of the checked rule (`["1"]` without one).
    pub left_hand_factor: Vec<String>,
    /// `kappa_j(n)` for `j = 1, ..., L` (entry `j - 1`; `[]` for a missing
    /// lag): the rule is `f(n) P_n = sum_j kappa_j(n) q(x)^j P_(n-j)`.
    pub scalar_multipliers: Vec<Vec<String>>,
    /// `e = N - 1`, the last seed row: `P_n = q^(n - e) g S_n` for `n >= e`.
    pub gauge_index: usize,
    /// The monic common factor `g(x)` of the gauged window rows.
    pub common_factor: Vec<String>,
    /// The residual window `T_m = q^(e - m) P_m / g`, `m = e - L + 1, ...,
    /// e`, which spans the residual space.
    pub residual_seeds: Vec<Vec<String>>,
    /// The shape of the residual space.
    pub residual_shape: FactorGaugeResidualShape,
    /// The coefficient polynomials `pi_k(n) = [x^k] S_n`, `k = 0, ...,
    /// degree`, each satisfying the scalar recurrence identically and
    /// agreeing with the window; present when found (always for shapes of
    /// degree at least two).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub residual_coefficients: Option<Vec<Vec<String>>>,
    /// The coefficient sequence that proves `S_n != 0` for every `n >= N`
    /// by a dominant characteristic root, when that is the argument used
    /// (module documentation, "Nonzero rows").
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dominant_root: Option<FactorGaugeDominantRoot>,
}

impl FactorGaugeCertificate {
    /// Replay the certificate on the checked rule (after a left-hand factor
    /// normalization) and seed rows (`seeds[i]` is the row
    /// `first_index + i`, the moved seeds when the induction start was
    /// moved): the rule is `f(n) P_n = sum_j kappa_j(n) q^j P_(n-j)` term
    /// by term, `q` is monic, `e` is the last seed row,
    /// `g T_m = q^(e - m) P_m` for every window row, the residual window
    /// has the recorded shape, and every coefficient polynomial satisfies
    /// the scalar recurrence identically and agrees with the window.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::VerificationFailed`] on any mismatch.
    pub fn verify(
        &self,
        rule: &Recurrence,
        first_index: usize,
        seeds: &[Vec<Q>],
    ) -> Result<(), RecurrenceError> {
        let fail = |what: &str| {
            Err(RecurrenceError::VerificationFailed(format!(
                "factor gauge: {what}"
            )))
        };
        let q = vec_poly::trim(&parse_polynomial(&self.gauge_factor)?);
        let f = vec_poly::trim(&parse_polynomial(&self.left_hand_factor)?);
        let kappa = self
            .scalar_multipliers
            .iter()
            .map(|k| parse_polynomial(k))
            .collect::<Result<Vec<_>, _>>()?;
        let order = kappa.len();
        if q.last().is_none_or(|leading| !leading.is_one()) || vec_poly::is_zero(&f) {
            return fail("the gauge factor is not monic or the left-hand factor is zero");
        }
        if order < 2 || seeds.len() < order {
            return fail("the order is below two or exceeds the seed rows");
        }
        let gauged = GaugedRule { q, f, kappa };
        let shape = match recurrence_shape(rule) {
            Ok(shape) => shape,
            Err(reason) => return fail(&reason),
        };
        if vec_poly::trim(&shape.denominator) != gauged.f {
            return fail("the left-hand factor differs from the rule's");
        }
        if shape
            .terms
            .keys()
            .any(|&(lag, derivative)| derivative != 0 || lag == 0 || lag > order)
        {
            return fail("the rule has a term outside the recorded lags");
        }
        for lag in 1..=order {
            let actual = shape.terms.get(&(lag, 0)).cloned().unwrap_or_default();
            let power = vec_poly::pow(&gauged.q, lag);
            let expected: Bivariate = gauged.kappa[lag - 1]
                .iter()
                .map(|c| vec_poly::scale(&power, c))
                .collect();
            let difference = add_bivariate(&actual, &negate_bivariate(&expected));
            if x_degree(&difference).is_some() {
                return fail(&format!(
                    "the multiplier of P_(n-{lag}) is not kappa_{lag}(n) q^{lag}"
                ));
            }
        }
        if self.gauge_index != first_index + seeds.len() - 1 {
            return fail("the gauge index is not the last seed row");
        }
        let g = vec_poly::trim(&parse_polynomial(&self.common_factor)?);
        if vec_poly::is_zero(&g) || self.residual_seeds.len() != order {
            return fail("the common factor is zero or the window has the wrong length");
        }
        let window = self
            .residual_seeds
            .iter()
            .map(|t| parse_polynomial(t).map(|t| vec_poly::trim(&t)))
            .collect::<Result<Vec<_>, _>>()?;
        for (i, t) in window.iter().enumerate() {
            let row = &seeds[seeds.len() - order + i];
            let gauged_row = vec_poly::mul(&vec_poly::pow(&gauged.q, order - 1 - i), row);
            if vec_poly::trim(&vec_poly::mul(&g, t)) != vec_poly::trim(&gauged_row) {
                return fail(&format!("window row {i} is not q^(e - m) P_m / g"));
            }
        }
        if residual_shape(&window) != Some(self.residual_shape) {
            return fail("the residual window does not have the recorded shape");
        }
        let degree = self.residual_shape.degree();
        match &self.residual_coefficients {
            None if degree >= 2 => return fail("the coefficient polynomials are missing"),
            None => {}
            Some(coefficients) => {
                if coefficients.len() != degree + 1 {
                    return fail("the coefficient polynomials have the wrong count");
                }
                let window_start = self.gauge_index + 1 - order;
                for (k, pi) in coefficients.iter().enumerate() {
                    let pi = parse_polynomial(pi)?;
                    if !gauged.annihilates(&pi) {
                        return fail(&format!("pi_{k} does not satisfy the scalar recurrence"));
                    }
                    for (i, t) in window.iter().enumerate() {
                        let n = index_value(window_start + i);
                        if vec_poly::evaluate(&pi, &n) != vec_poly::coeff(t, k) {
                            return fail(&format!("pi_{k} disagrees with window row {i}"));
                        }
                    }
                }
            }
        }
        if let Some(dominant) = &self.dominant_root {
            if dominant.coefficient > degree {
                return fail("the dominant-root coefficient exceeds the residual degree");
            }
            let window_start = self.gauge_index + 1 - order;
            dominant.verify(&gauged, &window, window_start, self.gauge_index + 1)?;
        }
        Ok(())
    }
}

/// Check the factor-gauge theorem (RealRooted routes `powerForm`,
/// `closedForm`, `factorPower`) for the family given by `rule`,
/// `first_index` and the seed rows `seeds` (the rule applies from
/// `N = first_index + seeds.len()` on).
///
/// The rule must be `f(n) P_n = sum_(j=1)^L kappa_j(n) q(x)^j P_(n-j)` with
/// one monic `q` (a constant allowed), order `L >= 2` and at least `L` seed
/// rows (no derivative, no `(-1)^n` factor, no inhomogeneous term).  With
/// `e = N - 1`, `P_n = q^(n-e) g S_n` for every `n >= e`, where `g` is the
/// gcd of the gauged window `q^(e-m) P_m` and `S_n` lies in the span of the
/// residual window (a scalar recurrence moves its coordinates).  The
/// hypotheses, each decided exactly:
///
/// - `f(n) != 0` for `n >= N` (sign checker);
/// - `q` and `g` are real-rooted;
/// - `S_n` is real-rooted for every `n >= N` by the shape of the residual
///   space: degree at most one; the discriminant (degree two) or the formal
///   cubic discriminant (degree three) is `>= 0`; for palindromic quartics,
///   a strict sign of the outer coefficient and one of three sign cases in
///   `u = x + 1/x` (module documentation); the coefficients of `S_n` must
///   be polynomials in `n` proved by an identity for degree two and more,
///   otherwise the theorem is `not_applicable`;
/// - `S_n != 0` for `n >= N` (a coefficient polynomial, or the shared
///   leading-coefficient certificate);
/// - the seed rows are nonzero and real-rooted;
/// - a requested root region must contain the certified one.
///
/// When all hold, every `P_n`, `n >= first_index`, is nonzero and
/// real-rooted, with every root in `(-inf, 0)` or `(-inf, 0]` when that is
/// certified; for a constant residual and a linear `q`, `P_(n-1)`
/// interlaces `P_n` (weakly) for every `n >= N`.  The data are recorded in
/// [`RealRootednessTheoremReport::factor_gauge`].  The induction start is
/// moved past the seeds as for the other theorems.
///
/// # Errors
///
/// Errors of the record itself (seed normalization, row generation for a
/// moved start).
///
/// # Example
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::recurrence::{
///     check_factor_gauge_theorem, BigRational, BivarPoly, FactorGaugeResidualShape, Recurrence,
///     RecurrenceTerm, RealRootednessTheoremOptions,
/// };
///
/// let q = |v: i64| BigRational::from_integer(BigInt::from(v));
/// let poly = |c: &[i64]| BivarPoly::from_coefficients_normalized(vec![c.iter().map(|&v| q(v)).collect()]);
/// // P_n = 3(1+x) P_(n-1) - 3(1+x)^2 P_(n-2) + (1+x)^3 P_(n-3) (OEIS A124376):
/// // P_n = (1+x)^(n-4) S_n with S_n a palindromic quartic.
/// let rule = Recurrence {
///     terms: vec![
///         RecurrenceTerm::new(1, 0, poly(&[3, 3])),
///         RecurrenceTerm::new(2, 0, poly(&[-3, -6, -3])),
///         RecurrenceTerm::new(3, 0, poly(&[1, 3, 3, 1])),
///     ],
///     denominator: None,
///     inhomogeneous: None,
/// };
/// let seeds: Vec<Vec<BigRational>> = [&[1][..], &[1, 1], &[1, 4, 1], &[1, 7, 7, 1], &[1, 10, 19, 10, 1]]
///     .iter()
///     .map(|row| row.iter().map(|&v| q(v)).collect())
///     .collect();
/// let report =
///     check_factor_gauge_theorem(&rule, 0, &seeds, &RealRootednessTheoremOptions::default())
///         .unwrap();
/// assert!(report.is_proved());
/// let certificate = report.factor_gauge().unwrap();
/// assert_eq!(certificate.residual_shape, FactorGaugeResidualShape::PalindromicQuartic);
/// report.verify().unwrap();
/// ```
pub fn check_factor_gauge_theorem(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<BigRational>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    check_real_rootedness_theorem(
        RealRootednessTheorem::FactorGauge,
        rule,
        first_index,
        seeds,
        options,
    )
}

// ---------------------------------------------------------------------------
// The gauged rule and the scalar recurrence
// ---------------------------------------------------------------------------

/// The rule as `f(n) P_n = sum_j kappa_j(n) q(x)^j P_(n-j)`.
struct GaugedRule {
    /// Monic, nonzero.
    q: Vec<Q>,
    /// Nonzero, trimmed.
    f: Vec<Q>,
    /// `kappa[j - 1]`, ascending in `n` (empty or `[0]` for a missing lag).
    kappa: Vec<Vec<Q>>,
}

impl GaugedRule {
    fn order(&self) -> usize {
        self.kappa.len()
    }

    /// Whether `f(n) pi(n) - sum_j kappa_j(n) pi(n - j)` is the zero
    /// polynomial in `n`.
    fn annihilates(&self, pi: &[Q]) -> bool {
        let mut total = vec_poly::mul(&self.f, pi);
        for (lag, kappa) in (1i64..).zip(&self.kappa) {
            let shifted = vec_poly::compose(pi, &[Q::from_integer(BigInt::from(-lag)), Q::one()]);
            total = vec_poly::sub(&total, &vec_poly::mul(kappa, &shifted));
        }
        vec_poly::is_zero(&total)
    }

    /// The rows `S_(start), ..., S_(start + count - 1)` from the window
    /// `window = (S_start, ..., S_(start + L - 1))`; `None` where `f(n) = 0`.
    fn residual_rows(&self, window: &[Vec<Q>], start: usize, count: usize) -> Option<Vec<Vec<Q>>> {
        let mut rows = window.to_vec();
        while rows.len() < count {
            let n = index_value(start + rows.len());
            let f = vec_poly::evaluate(&self.f, &n);
            if f.is_zero() {
                return None;
            }
            let mut sum = vec![Q::zero()];
            for (lag, kappa) in (1usize..).zip(&self.kappa) {
                let scalar = vec_poly::evaluate(kappa, &n);
                if !scalar.is_zero() {
                    sum = vec_poly::add(&sum, &vec_poly::scale(&rows[rows.len() - lag], &scalar));
                }
            }
            rows.push(vec_poly::trim(&vec_poly::scale(&sum, &f.recip())));
        }
        Some(rows)
    }

    /// The largest degree in `n` of `f` and the `kappa_j`.
    fn index_degree(&self) -> usize {
        std::iter::once(&self.f)
            .chain(&self.kappa)
            .filter_map(|p| vec_poly::degree(p))
            .max()
            .unwrap_or(0)
    }
}

fn index_value(n: usize) -> Q {
    Q::from_integer(BigInt::from(n))
}

fn negate_bivariate(b: &Bivariate) -> Bivariate {
    b.iter().map(|row| vec_poly::neg(row)).collect()
}

fn parse_polynomial(coefficients: &[String]) -> Result<Vec<Q>, RecurrenceError> {
    coefficients
        .iter()
        .map(|c| {
            parse_rational_coeff(c).map_err(|error| {
                RecurrenceError::VerificationFailed(format!("factor gauge: {error}"))
            })
        })
        .collect()
}

fn format_polynomial(p: &[Q]) -> Vec<String> {
    vec_poly::trim(p)
        .iter()
        .map(format_rational_coeff)
        .collect()
}

/// The monic `q` with `q^power = f` for a monic `f`, when it exists:
/// `x^deg f f(1/x)` has constant term `1`, so its `power`-th root as a power
/// series is determined coefficient by coefficient; the candidate is
/// checked by raising it to the power.
fn monic_polynomial_root(f: &[Q], power: usize) -> Option<Vec<Q>> {
    let f = vec_poly::trim(f);
    let degree = vec_poly::degree(&f)?;
    if power == 0 || degree % power != 0 || !f[degree].is_one() {
        return None;
    }
    let reversed: Vec<Q> = f.iter().rev().cloned().collect();
    let root_degree = degree / power;
    let exponent = index_value(power);
    let mut root = vec![Q::one()];
    for (k, target) in reversed.iter().enumerate().take(root_degree + 1).skip(1) {
        let mut product = vec![Q::one()];
        for _ in 0..power {
            product = vec_poly::mul(&product, &root);
            product.truncate(k + 1);
        }
        let current = vec_poly::coeff(&product, k);
        root.push((target.clone() - current) / &exponent);
    }
    let candidate: Vec<Q> = root.into_iter().rev().collect();
    (vec_poly::trim(&vec_poly::pow(&candidate, power)) == f).then_some(candidate)
}

/// Match the rule (module documentation), or say why it does not match.
fn match_gauged_rule(shape: &RecurrenceShape) -> Result<GaugedRule, String> {
    if shape.terms.is_empty() {
        return Err("the rule has no nonzero right-hand term".into());
    }
    if shape.terms.keys().any(|&(_, derivative)| derivative != 0) {
        return Err("the rule has a derivative term".into());
    }
    let order = shape.terms.keys().map(|&(lag, _)| lag).max().unwrap_or(0);
    if order < 2 {
        return Err("a first-order rule is the product route's case (order 1)".into());
    }
    let mut separated = BTreeMap::new();
    for (&(lag, _), coefficients) in &shape.terms {
        let (gamma, w) = separate(coefficients).ok_or_else(|| {
            format!("the multiplier of P_(n-{lag}) is not a product gamma(n) w(x)")
        })?;
        separated.insert(lag, (gamma, w));
    }
    let (&smallest, (_, w)) = separated.iter().next().expect("a nonzero term");
    let q = monic_polynomial_root(&monic(w), smallest).ok_or_else(|| {
        format!(
            "the multiplier of P_(n-{smallest}) is not a constant times a {smallest}-th power of a \
             polynomial in x"
        )
    })?;
    let mut kappa = vec![Vec::new(); order];
    for (lag, (gamma, w)) in separated {
        if monic(&w) != vec_poly::trim(&vec_poly::pow(&q, lag)) {
            return Err(format!(
                "the multipliers are not kappa_j(n) q(x)^j for one q (q = {} from lag \
                 {smallest}, lag {lag} differs)",
                format_index_polynomial(&q, "x")
            ));
        }
        let leading = vec_poly::trim(&w).last().cloned().expect("nonzero");
        kappa[lag - 1] = vec_poly::scale(&gamma, &leading);
    }
    Ok(GaugedRule {
        q,
        f: vec_poly::trim(&shape.denominator),
        kappa,
    })
}

/// The shape of the space spanned by `window`, or `None` when it is not
/// one of the decided shapes (the zero polynomial is of every shape).
fn residual_shape(window: &[Vec<Q>]) -> Option<FactorGaugeResidualShape> {
    let degree = window.iter().filter_map(|t| vec_poly::degree(t)).max()?;
    Some(match degree {
        0 => FactorGaugeResidualShape::Constant,
        1 => FactorGaugeResidualShape::Linear,
        2 => FactorGaugeResidualShape::Quadratic,
        3 => FactorGaugeResidualShape::Cubic,
        4 if window.iter().all(|t| {
            let c = |k| vec_poly::coeff(t, k);
            c(0) == c(4) && c(1) == c(3)
        }) =>
        {
            FactorGaugeResidualShape::PalindromicQuartic
        }
        _ => return None,
    })
}

/// The coefficient polynomials `pi_k(n)`, `k = 0, ..., degree`, of the
/// residual rows (module documentation), or `None`.
fn polynomial_coordinates(
    rule: &GaugedRule,
    window: &[Vec<Q>],
    start: usize,
    degree: usize,
) -> Option<Vec<Vec<Q>>> {
    let order = rule.order();
    let bound = (order - 1 + 2 * rule.index_degree() + 2).min(MAX_COORDINATE_DEGREE);
    let count = bound + order + 2;
    let rows = rule.residual_rows(window, start, count)?;
    let points: Vec<Q> = (0..count).map(|t| index_value(start + t)).collect();
    let mut coefficients = Vec::with_capacity(degree + 1);
    for k in 0..=degree {
        let values: Vec<Q> = rows.iter().map(|row| vec_poly::coeff(row, k)).collect();
        let pi = (0..=bound).find_map(|d| {
            let candidate = Polynomial::lagrange_interpolation(&points[..=d], &values[..=d])
                .coeffs()
                .to_vec();
            points
                .iter()
                .zip(&values)
                .all(|(point, value)| &vec_poly::evaluate(&candidate, point) == value)
                .then_some(candidate)
        })?;
        if !rule.annihilates(&pi) {
            return None;
        }
        coefficients.push(vec_poly::trim(&pi));
    }
    Some(coefficients)
}

// ---------------------------------------------------------------------------
// Root sides
// ---------------------------------------------------------------------------

/// Where the roots of a real-rooted polynomial lie, read from its
/// coefficient signs: all coefficients nonzero of one sign (every root
/// negative), the nonzero ones of one sign (every root `<= 0`), or neither.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum RootSide {
    Negative,
    NonPositive,
    Other,
}

fn fixed_root_side(p: &[Q]) -> RootSide {
    let p = vec_poly::trim(p);
    let positive = p.iter().all(|c| c.is_positive());
    let negative = p.iter().all(|c| c.is_negative());
    if positive || negative {
        return RootSide::Negative;
    }
    let nonzero: Vec<&Q> = p.iter().filter(|c| !c.is_zero()).collect();
    if nonzero.iter().all(|c| c.is_positive()) || nonzero.iter().all(|c| c.is_negative()) {
        RootSide::NonPositive
    } else {
        RootSide::Other
    }
}

/// The root side of `S_n` for every `n >= from` from its coefficient
/// polynomials, with the (verified) sign hypotheses that give it.
fn residual_root_side(coefficients: &[Vec<Q>], from: usize) -> (RootSide, Vec<LiftHypothesis>) {
    let attempts = [
        (IndexSignCondition::Positive, "positive", RootSide::Negative),
        (IndexSignCondition::Negative, "negative", RootSide::Negative),
        (
            IndexSignCondition::Nonnegative,
            "nonnegative",
            RootSide::NonPositive,
        ),
        (
            IndexSignCondition::Nonpositive,
            "nonpositive",
            RootSide::NonPositive,
        ),
    ];
    for (condition, word, side) in attempts {
        let mut hypotheses = Vec::new();
        for (k, pi) in coefficients.iter().enumerate() {
            let function = IndexFunction::polynomial(pi.clone());
            let hypothesis = sign_hypothesis(
                &format!("residual_coefficient_{k}_{word}"),
                &format!("[x^{k}] S_n is {word} for every n >= {from}"),
                &format!("[x^{k}] S_n"),
                &function,
                from,
                condition,
            );
            if hypothesis.status != LiftHypothesisStatus::Verified {
                break;
            }
            hypotheses.push(hypothesis);
        }
        if hypotheses.len() == coefficients.len() {
            return (side, hypotheses);
        }
    }
    (RootSide::Other, Vec::new())
}

// ---------------------------------------------------------------------------
// The report
// ---------------------------------------------------------------------------

/// Which palindromic case proved the quartic (for the Lean declarations).
#[derive(Clone, Copy, PartialEq, Eq)]
enum PalindromicCase {
    BelowMinusTwo,
    AboveTwo,
    Straddle,
}

/// The real-rootedness hypotheses of `S_n` for the shape (module
/// documentation), and the palindromic case used.
fn residual_hypotheses(
    shape: FactorGaugeResidualShape,
    coefficients: Option<&[Vec<Q>]>,
    from: usize,
) -> (Vec<LiftHypothesis>, Option<PalindromicCase>) {
    let p = |pi: &[Q], qi: &[Q]| vec_poly::mul(pi, qi);
    let s = |c: i64, pi: &[Q]| vec_poly::scale(pi, &Q::from_integer(BigInt::from(c)));
    let sum = |parts: &[Vec<Q>]| {
        parts
            .iter()
            .fold(vec![Q::zero()], |total, part| vec_poly::add(&total, part))
    };
    let degree_at_most_one = || {
        vec![exact_hypothesis(
            "residual_degree_at_most_one",
            "every element of the residual space has degree at most 1, so every nonzero S_n is \
             real-rooted",
            Ok("exact: every residual window row has degree at most 1".into()),
        )]
    };
    let c = match (shape, coefficients) {
        (FactorGaugeResidualShape::Constant | FactorGaugeResidualShape::Linear, _) => {
            return (degree_at_most_one(), None);
        }
        (_, None) => {
            return (
                vec![unresolved_hypothesis(
                    "residual_real_rooted",
                    "S_n is real-rooted for every n >= N",
                    "the coefficients of S_n are not polynomials in n".into(),
                )],
                None,
            );
        }
        (_, Some(c)) => c,
    };
    match shape {
        FactorGaugeResidualShape::Constant | FactorGaugeResidualShape::Linear => {
            (degree_at_most_one(), None)
        }
        FactorGaugeResidualShape::Quadratic => {
            // B^2 - 4 A C with S_n = A x^2 + B x + C.
            let discriminant = sum(&[p(&c[1], &c[1]), s(-4, &p(&c[2], &c[0]))]);
            (
                vec![sign_hypothesis(
                    "residual_discriminant_nonnegative",
                    &format!(
                        "the discriminant B(n)^2 - 4 A(n) C(n) of S_n = A(n) x^2 + B(n) x + C(n) \
                         is >= 0 for every n >= {from} (for A(n) = 0 it is B(n)^2 and S_n has \
                         degree at most 1)"
                    ),
                    "B(n)^2 - 4 A(n) C(n)",
                    &IndexFunction::polynomial(discriminant),
                    from,
                    IndexSignCondition::Nonnegative,
                )],
                None,
            )
        }
        FactorGaugeResidualShape::Cubic => {
            let (a, b, cc, d) = (&c[3], &c[2], &c[1], &c[0]);
            let discriminant = sum(&[
                s(18, &p(&p(a, b), &p(cc, d))),
                s(-4, &p(&p(&p(b, b), b), d)),
                p(&p(b, b), &p(cc, cc)),
                s(-4, &p(&p(a, cc), &p(cc, cc))),
                s(-27, &p(&p(a, a), &p(d, d))),
            ]);
            (
                vec![sign_hypothesis(
                    "residual_cubic_discriminant_nonnegative",
                    &format!(
                        "the formal discriminant 18abcd - 4b^3 d + b^2 c^2 - 4ac^3 - 27a^2 d^2 of \
                         S_n = a x^3 + b x^2 + c x + d is >= 0 for every n >= {from} (for a = 0 \
                         it is b^2 (c^2 - 4bd))"
                    ),
                    "the cubic discriminant",
                    &IndexFunction::polynomial(discriminant),
                    from,
                    IndexSignCondition::Nonnegative,
                )],
                None,
            )
        }
        FactorGaugeResidualShape::PalindromicQuartic => {
            palindromic_hypotheses(&c[0], &c[1], &c[2], from)
        }
    }
}

fn palindromic_hypotheses(
    a: &[Q],
    b: &[Q],
    c: &[Q],
    from: usize,
) -> (Vec<LiftHypothesis>, Option<PalindromicCase>) {
    let outer = |condition, word: &str| {
        sign_hypothesis(
            "residual_outer_coefficient_sign",
            &format!(
                "a(n) = [x^0] S_n = [x^4] S_n is {word} for every n >= {from} (S_n = a x^4 + b \
                 x^3 + c x^2 + b x + a)"
            ),
            "a(n)",
            &IndexFunction::polynomial(a.to_vec()),
            from,
            condition,
        )
    };
    let positive = outer(IndexSignCondition::Positive, "positive");
    let (sign, outer_hypothesis) = if positive.status == LiftHypothesisStatus::Verified {
        (1, positive)
    } else {
        let negative = outer(IndexSignCondition::Negative, "negative");
        if negative.status == LiftHypothesisStatus::Verified {
            (-1, negative)
        } else {
            return (vec![positive], None);
        }
    };
    let scale = |p: &[Q], factor: i64| vec_poly::scale(p, &Q::from_integer(BigInt::from(factor)));
    let (a, b, c) = (scale(a, sign), scale(b, sign), scale(c, sign));
    let combination = |terms: &[(i64, &Vec<Q>)]| {
        terms.iter().fold(vec![Q::zero()], |total, (factor, p)| {
            vec_poly::add(&total, &scale(p, *factor))
        })
    };
    // b^2 - 4 a (c - 2a), h(2) = 2a + 2b + c, h(-2) = 2a - 2b + c.
    let discriminant = combination(&[
        (1, &vec_poly::mul(&b, &b)),
        (-4, &vec_poly::mul(&a, &c)),
        (8, &vec_poly::mul(&a, &a)),
    ]);
    let at_two = combination(&[(2, &a), (2, &b), (1, &c)]);
    let at_minus_two = combination(&[(2, &a), (-2, &b), (1, &c)]);
    let normalized = if sign == 1 {
        ""
    } else {
        " (after S_n -> -S_n)"
    };
    let vertex_below = combination(&[(1, &b), (-4, &a)]);
    let vertex_above = combination(&[(-1, &b), (-4, &a)]);
    let cases: [(PalindromicCase, Vec<PalindromicCondition<'_>>); 3] = [
        (
            PalindromicCase::BelowMinusTwo,
            vec![
                (
                    "residual_u_vertex_at_most_minus_two",
                    "b - 4a >= 0 (the vertex of h is at most -2)",
                    &vertex_below,
                    IndexSignCondition::Nonnegative,
                ),
                (
                    "residual_u_value_at_minus_two_nonnegative",
                    "h(-2) = 2a - 2b + c >= 0",
                    &at_minus_two,
                    IndexSignCondition::Nonnegative,
                ),
                (
                    "residual_u_discriminant_nonnegative",
                    "b^2 - 4a(c - 2a) >= 0 (h has real roots)",
                    &discriminant,
                    IndexSignCondition::Nonnegative,
                ),
            ],
        ),
        (
            PalindromicCase::AboveTwo,
            vec![
                (
                    "residual_u_vertex_at_least_two",
                    "-b - 4a >= 0 (the vertex of h is at least 2)",
                    &vertex_above,
                    IndexSignCondition::Nonnegative,
                ),
                (
                    "residual_u_value_at_two_nonnegative",
                    "h(2) = 2a + 2b + c >= 0",
                    &at_two,
                    IndexSignCondition::Nonnegative,
                ),
                (
                    "residual_u_discriminant_nonnegative",
                    "b^2 - 4a(c - 2a) >= 0 (h has real roots)",
                    &discriminant,
                    IndexSignCondition::Nonnegative,
                ),
            ],
        ),
        (
            PalindromicCase::Straddle,
            vec![
                (
                    "residual_u_value_at_two_nonpositive",
                    "h(2) = 2a + 2b + c <= 0",
                    &at_two,
                    IndexSignCondition::Nonpositive,
                ),
                (
                    "residual_u_value_at_minus_two_nonpositive",
                    "h(-2) = 2a - 2b + c <= 0",
                    &at_minus_two,
                    IndexSignCondition::Nonpositive,
                ),
            ],
        ),
    ];
    let mut evaluated = Vec::new();
    for (case, conditions) in &cases {
        let hypotheses: Vec<LiftHypothesis> = conditions
            .iter()
            .map(|(name, statement, polynomial, condition)| {
                sign_hypothesis(
                    name,
                    &format!(
                        "{statement} for every n >= {from}, with h(u) = a u^2 + b u + (c - 2a), \
                         u = x + 1/x{normalized}"
                    ),
                    name,
                    &IndexFunction::polynomial((*polynomial).clone()),
                    from,
                    *condition,
                )
            })
            .collect();
        if hypotheses
            .iter()
            .all(|h| h.status == LiftHypothesisStatus::Verified)
        {
            let mut all = vec![outer_hypothesis];
            all.extend(hypotheses);
            return (all, Some(*case));
        }
        evaluated.push(hypotheses);
    }
    // No case holds: report the case on the side of the vertex of h at the
    // first index (b >= 0: both roots of h at most -2).
    let b_first = IndexFunction::polynomial(b.clone())
        .at(from)
        .unwrap_or_else(Q::zero);
    let reported = usize::from(b_first.is_negative());
    let mut all = vec![outer_hypothesis];
    all.extend(evaluated.swap_remove(reported));
    (all, None)
}

/// One sign condition of a palindromic case: name, statement, polynomial
/// in `n` and condition.
type PalindromicCondition<'a> = (&'static str, &'static str, &'a Vec<Q>, IndexSignCondition);

/// The registry entry of the factor gauge ([`super::registry::ROUTES`]).
pub(super) const ROUTE: TheoremRoute = TheoremRoute {
    theorem: RealRootednessTheorem::FactorGauge,
    check_at_supplied_seeds: factor_gauge_report,
    start_policy: StartPolicy::MovesStart,
    is_preferred: any_proof,
    certificate_name: Some("factor-gauge"),
};

/// Replay the certificate of a recomputed factor-gauge report: every
/// identity of the certificate on the checked rule and the (possibly moved)
/// seed rows ([`FactorGaugeCertificate::verify`]).
pub(super) fn verify_certificate_in_report(
    certificate: &FactorGaugeCertificate,
    report: &RealRootednessTheoremReport,
) -> Result<(), RecurrenceError> {
    let (rule, first_index, seeds) = report.recurrence.to_recurrence_parts()?;
    let rule = match &report.left_hand_factor_normalization {
        Some(normalization) => normalization.normalized_recurrence()?,
        None => rule,
    };
    let seeds = if report.start_delay == 0 {
        seeds
    } else {
        let generated = rule
            .generate_rows_rational(&seeds, first_index, seeds.len() + report.start_delay)
            .map_err(|error| RecurrenceError::VerificationFailed(error.to_string()))?;
        normalized_seed_rows(&rule, first_index, &generated)?
    };
    certificate.verify(&rule, first_index, &seeds)
}

fn factor_gauge_report(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    let mut report = RealRootednessTheoremReport::new(
        RealRootednessTheorem::FactorGauge,
        REFERENCE,
        rule,
        first_index,
        seeds,
        options,
    );
    let shape = match recurrence_shape(rule) {
        Ok(shape) => shape,
        Err(reason) => return Ok(report.not_applicable(reason)),
    };
    let gauged = match match_gauged_rule(&shape) {
        Ok(gauged) => gauged,
        Err(reason) => return Ok(report.not_applicable(reason)),
    };
    let order = gauged.order();
    if seeds.len() < order {
        return Ok(report.not_applicable(format!(
            "the rule has order {order} but only {} seed row(s)",
            seeds.len()
        )));
    }
    let n_first = first_index + seeds.len();
    let gauge_index = n_first - 1;
    let window_start = n_first - order;
    let gauged_rows: Vec<Vec<Q>> = (0..order)
        .map(|i| {
            vec_poly::trim(&vec_poly::mul(
                &vec_poly::pow(&gauged.q, order - 1 - i),
                &seeds[seeds.len() - order + i],
            ))
        })
        .collect();
    let common = gauged_rows.iter().fold(Polynomial::<Q>::zero(), |g, row| {
        g.gcd(&Polynomial::new(row.clone()))
    });
    if common.is_zero() {
        return Ok(report.not_applicable(format!(
            "the window rows P_{window_start}, ..., P_{gauge_index} are all zero"
        )));
    }
    let g = vec_poly::trim(common.coeffs());
    let window: Vec<Vec<Q>> = gauged_rows
        .iter()
        .map(|row| vec_poly::trim(Polynomial::new(row.clone()).exact_div(&common).coeffs()))
        .collect();
    let Some(residual) = residual_shape(&window) else {
        let degree = window.iter().filter_map(|t| vec_poly::degree(t)).max();
        return Ok(report.not_applicable(format!(
            "the residual space has degree {} and is not a palindromic quartic: only degree at \
             most 3 and palindromic quartics are decided",
            degree.unwrap_or(0)
        )));
    };
    let degree = residual.degree();
    let coefficients = polynomial_coordinates(&gauged, &window, window_start, degree);
    if degree >= 2 && coefficients.is_none() {
        return Ok(report.not_applicable(format!(
            "the coefficients of the {} residual S_n are not polynomials in n of degree at most \
             {MAX_COORDINATE_DEGREE} satisfying the scalar recurrence (signs of other scalar \
             sequences are not decided)",
            residual.text()
        )));
    }

    let multipliers: Vec<String> = (1..=order)
        .filter(|lag| vec_poly::degree(&gauged.kappa[lag - 1]).is_some())
        .map(|lag| {
            let kappa = IndexFunction::new(gauged.kappa[lag - 1].clone(), gauged.f.clone());
            format!("({}) q^{lag} P_(n-{lag})", kappa.text())
        })
        .collect();
    let q_text = format_index_polynomial(&gauged.q, "x");
    let g_text = format_index_polynomial(&g, "x");
    report.normal_form = format!(
        "P_n = {} for n >= {n_first}, with q = {q_text}; P_n = q^(n-{gauge_index}) g S_n for n >= \
         {gauge_index} with g = {g_text} and S_n in the {} span of the residual window T_m = \
         q^({gauge_index}-m) P_m / g, m = {window_start}, ..., {gauge_index}",
        multipliers.join(" + "),
        residual.text()
    );

    let mut hypotheses: Vec<LiftHypothesis> = Vec::new();
    hypotheses.extend(denominator_hypothesis(&gauged.f, n_first));
    if vec_poly::degree(&gauged.q).unwrap_or(0) >= 1 {
        hypotheses.push(exact_hypothesis(
            "gauge_factor_real_rooted",
            &format!("the gauge factor q = {q_text} is real-rooted"),
            if is_real_rooted_rational(&gauged.q) {
                Ok("exact".into())
            } else {
                Err((
                    Some(n_first),
                    format!(
                        "q is not real-rooted (exact), so neither is any nonzero P_n, n >= \
                         {n_first}"
                    ),
                ))
            },
        ));
    }
    if vec_poly::degree(&g).unwrap_or(0) >= 1 {
        hypotheses.push(exact_hypothesis(
            "common_factor_real_rooted",
            &format!("the common factor g = {g_text} of the gauged window is real-rooted"),
            if is_real_rooted_rational(&g) {
                Ok("exact".into())
            } else {
                Err((
                    Some(gauge_index),
                    format!(
                        "g is not real-rooted (exact), so neither is any nonzero P_n, n >= \
                         {gauge_index}"
                    ),
                ))
            },
        ));
    }
    let (residual_checks, palindromic) =
        residual_hypotheses(residual, coefficients.as_deref(), n_first);
    hypotheses.extend(residual_checks);
    let (nonzero, dominant_root) = rows_nonzero_hypothesis(
        &gauged,
        &window,
        window_start,
        coefficients.as_deref(),
        (rule, first_index, seeds),
    );
    hypotheses.push(nonzero);
    hypotheses.extend(earlier_seeds_hypothesis(seeds, first_index, None));

    // Root region (module documentation): claimed only when proved.
    let fixed_side = std::iter::once(&gauged.q)
        .chain(std::iter::once(&g))
        .chain(seeds)
        .map(|p| fixed_root_side(p))
        .max()
        .unwrap_or(RootSide::Negative);
    let (residual_side, side_hypotheses) = match (residual, &coefficients) {
        (FactorGaugeResidualShape::Constant, _) => (RootSide::Negative, Vec::new()),
        (_, Some(coefficients)) if fixed_side != RootSide::Other => {
            residual_root_side(coefficients, n_first)
        }
        _ => (RootSide::Other, Vec::new()),
    };
    let side = fixed_side.max(residual_side);
    let certified = match side {
        RootSide::Negative => Some(RealInterval {
            lower: None,
            lower_closed: false,
            upper: Some(Q::zero()),
            upper_closed: false,
        }),
        RootSide::NonPositive => Some(RealInterval {
            lower: None,
            lower_closed: false,
            upper: Some(Q::zero()),
            upper_closed: true,
        }),
        RootSide::Other => None,
    };
    let region_text = |region: &RealInterval| RootInterval::from_real(region).to_string();
    let region_evidence = |region: &RealInterval| {
        format!(
            "{}: q, g and the seed rows have coefficients of one sign (exact){}",
            region_text(region),
            if residual == FactorGaugeResidualShape::Constant {
                "; S_n is a nonzero constant"
            } else {
                "; so do the coefficients of S_n for every n >= N (the sign hypotheses above)"
            }
        )
    };
    let mut conclusion_region = None;
    match &options.root_region {
        None => {
            if let Some(region) = &certified {
                hypotheses.extend(side_hypotheses);
                hypotheses.push(exact_hypothesis(
                    "root_region",
                    &format!("every root of every row lies in {}", region_text(region)),
                    Ok(region_evidence(region)),
                ));
                conclusion_region = Some(RootInterval::from_real(region));
            }
        }
        Some(requested) => {
            let whole_line = requested.lower.is_none() && requested.upper.is_none();
            let contained = whole_line
                || (requested.lower.is_none()
                    && requested.upper.as_ref().is_some_and(|upper| match side {
                        RootSide::Negative => !upper.is_negative(),
                        RootSide::NonPositive => upper.is_positive(),
                        RootSide::Other => false,
                    }));
            let statement = format!(
                "the requested region {} contains every root of every row",
                open_root_interval(requested)
            );
            if contained {
                if certified.is_some() {
                    hypotheses.extend(side_hypotheses);
                }
                hypotheses.push(exact_hypothesis(
                    "requested_root_region",
                    &statement,
                    Ok(match &certified {
                        Some(region) => region_evidence(region),
                        None => "the requested region is the real line".into(),
                    }),
                ));
                conclusion_region = certified.as_ref().map(RootInterval::from_real);
            } else {
                hypotheses.push(unresolved_hypothesis(
                    "requested_root_region",
                    &statement,
                    match &certified {
                        Some(region) => {
                            format!("this theorem certifies {} only", region_text(region))
                        }
                        None => "this theorem certifies no root region for this family".into(),
                    },
                ));
            }
        }
    }
    report.hypotheses = hypotheses;

    let linear_gauge = vec_poly::degree(&gauged.q) == Some(1);
    let interlaces = residual == FactorGaugeResidualShape::Constant && linear_gauge;
    let (lean_route, declarations): (Option<&str>, Vec<&str>) = match residual {
        FactorGaugeResidualShape::Constant => {
            let mut declarations = Vec::new();
            if interlaces {
                declarations.push("RealRooted.interlaces_of_forall_eq_C_mul_pow");
            }
            match order {
                2 => declarations.push("RealRooted.forall_eq_C_mul_pow_of_rec2"),
                3 => declarations.push("RealRooted.forall_eq_C_mul_pow_of_rec3"),
                _ => {}
            }
            (Some("powerForm"), declarations)
        }
        FactorGaugeResidualShape::Linear => (Some("closedForm"), Vec::new()),
        FactorGaugeResidualShape::Quadratic => (
            Some("closedForm"),
            vec!["RealRooted.splits_C_add_C_mul_X_add_C_mul_X_sq"],
        ),
        FactorGaugeResidualShape::Cubic => (None, Vec::new()),
        FactorGaugeResidualShape::PalindromicQuartic => match palindromic {
            Some(PalindromicCase::Straddle) => (None, Vec::new()),
            _ => (
                Some("closedForm"),
                vec!["RealRooted.splits_palindromic_quartic"],
            ),
        },
    };
    report.lean_route = lean_route.map(String::from);
    report.lean_declarations = declarations.into_iter().map(String::from).collect();
    if order == 2 && !interlaces {
        report
            .lean_declarations
            .push("RealRooted.eq_pow_mul_of_factoredThreeTerm".into());
    }
    report.certificate = Some(RouteCertificate::FactorGauge(FactorGaugeCertificate {
        gauge_factor: format_polynomial(&gauged.q),
        left_hand_factor: format_polynomial(&gauged.f),
        scalar_multipliers: gauged
            .kappa
            .iter()
            .map(|kappa| {
                let mut kappa = kappa.clone();
                vec_poly::strip_trailing_zeros(&mut kappa);
                kappa.iter().map(format_rational_coeff).collect()
            })
            .collect(),
        gauge_index,
        common_factor: format_polynomial(&g),
        residual_seeds: window.iter().map(|t| format_polynomial(t)).collect(),
        residual_shape: residual,
        residual_coefficients: coefficients
            .as_ref()
            .map(|c| c.iter().map(|pi| format_polynomial(pi)).collect()),
        dominant_root,
    }));

    let mut statement = format!(
        "every P_n, n >= {first_index}, is nonzero and real-rooted: P_n = q^(n-{gauge_index}) g \
         S_n for n >= {gauge_index} with q = {q_text}, g = {g_text} and S_n in the {} residual \
         space (factor gauge)",
        residual.text()
    );
    if let Some(region) = &conclusion_region {
        statement.push_str(&format!("; every root lies in {region}"));
    }
    if interlaces {
        statement.push_str(&format!(
            "; P_n = (s(n)/s(n-1)) q P_(n-1) with scalars s(n) != 0, so P_(n-1) interlaces P_n \
             for every n >= {n_first}"
        ));
    }
    Ok(report.finish(RealRootednessTheoremConclusion {
        first_index,
        zero_residue_classes: None,
        interlacing_from: interlaces.then_some(n_first),
        interlacing: interlaces.then_some(RowInterlacing::Weak),
        simple_roots: false,
        root_region: conclusion_region,
        statement,
    }))
}

/// `S_n != 0` for every `n >= N` (module documentation, "Nonzero rows").
fn rows_nonzero_hypothesis(
    gauged: &GaugedRule,
    window: &[Vec<Q>],
    window_start: usize,
    coefficients: Option<&[Vec<Q>]>,
    (rule, first_index, seeds): (&Recurrence, usize, &[Vec<Q>]),
) -> (LiftHypothesis, Option<FactorGaugeDominantRoot>) {
    let n_first = first_index + seeds.len();
    let name = "rows_nonzero";
    let statement = format!("S_n != 0, so P_n != 0, for every n >= {n_first}");
    let top = coefficients.and_then(|c| c.last()).map(|top| {
        let degree = coefficients.map_or(0, |c| c.len() - 1);
        sign_hypothesis(
            name,
            &statement,
            &format!("[x^{degree}] S_n"),
            &IndexFunction::polynomial(top.clone()),
            n_first,
            IndexSignCondition::Nonzero,
        )
    });
    if let Some(hypothesis) = top
        .as_ref()
        .filter(|h| h.status == LiftHypothesisStatus::Verified)
    {
        return (hypothesis.clone(), None);
    }
    let degree = window
        .iter()
        .filter_map(|t| vec_poly::degree(t))
        .max()
        .unwrap_or(0);
    if let Some(dominant) = (0..=degree)
        .rev()
        .find_map(|k| dominant_root_form(gauged, window, window_start, k, n_first))
    {
        let hypothesis = exact_hypothesis(
            name,
            &statement,
            Ok(format!(
                "[x^{}] S_n = {} for every n >= {window_start} (constant scalar recurrence with \
                 distinct rational characteristic roots; the window agrees); the first root \
                 dominates from n = {} on{} (exact); recorded in the certificate",
                dominant.coefficient,
                dominant.text(),
                dominant.dominance_from,
                if dominant.dominance_from > n_first {
                    format!(
                        " and the values for n = {n_first}, ..., {} are nonzero",
                        dominant.dominance_from - 1
                    )
                } else {
                    String::new()
                }
            )),
        );
        return (hypothesis, Some(dominant));
    }
    if let Ok(certificate) = certify_leading_coefficient(
        rule,
        first_index,
        seeds,
        &LeadingCoefficientOptions::default(),
    ) {
        if certificate.is_proved() && certificate.holds_from.is_some_and(|from| from <= n_first) {
            let hypothesis = exact_hypothesis(
                name,
                &statement,
                Ok(format!(
                    "shared leading-coefficient certificate, route {} (RealRooted {}): {}; the \
                     certificate is recomputed when the report is replayed",
                    certificate.lean_route.as_deref().unwrap_or("-"),
                    certificate.lean_declarations.join(", "),
                    certificate.conclusion.as_deref().unwrap_or_default()
                )),
            );
            return (hypothesis, None);
        }
    }
    // An exactly zero computed row is a failure; otherwise unresolved.
    let rows = gauged.residual_rows(window, window_start, window.len() + GENERATED_ROWS);
    let zero = rows.as_ref().and_then(|rows| {
        rows.iter()
            .enumerate()
            .skip(window.len())
            .find(|(_, row)| vec_poly::is_zero(row))
            .map(|(t, _)| window_start + t)
    });
    let hypothesis = match zero {
        Some(index) => exact_hypothesis(
            name,
            &statement,
            Err((Some(index), format!("P_{index} = 0 (exact)"))),
        ),
        None => unresolved_hypothesis(
            name,
            &statement,
            match top {
                Some(hypothesis) => format!(
                    "{}; the shared leading-coefficient certificate is not proved",
                    hypothesis.evidence
                ),
                None => "the coefficients of S_n are not polynomials in n, no coefficient has a \
                         dominant characteristic root, and the shared leading-coefficient \
                         certificate is not proved"
                    .into(),
            },
        ),
    };
    (hypothesis, None)
}

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

    fn text(&self) -> String {
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
    fn verify(
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
fn dominant_root_form(
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

#[cfg(test)]
mod tests {
    use super::*;

    fn q(v: i64) -> Q {
        Q::from_integer(BigInt::from(v))
    }

    #[test]
    fn polynomial_roots_are_exact() {
        // (1 + x)^3 and (x^2 + 3x + 1)^2.
        let cube = vec_poly::pow(&[q(1), q(1)], 3);
        assert_eq!(monic_polynomial_root(&cube, 3), Some(vec![q(1), q(1)]));
        let square = vec_poly::pow(&[q(1), q(3), q(1)], 2);
        assert_eq!(
            monic_polynomial_root(&square, 2),
            Some(vec![q(1), q(3), q(1)])
        );
        assert_eq!(
            monic_polynomial_root(&[q(0), q(0), q(1)], 2),
            Some(vec![q(0), q(1)])
        );
        assert_eq!(monic_polynomial_root(&[q(1), q(0), q(1)], 2), None);
        assert_eq!(monic_polynomial_root(&[q(1), q(1)], 2), None);
        assert_eq!(monic_polynomial_root(&[q(5)], 1), None);
        assert_eq!(monic_polynomial_root(&[q(1)], 3), Some(vec![q(1)]));
    }

    #[test]
    fn root_sides_follow_the_coefficient_signs() {
        assert_eq!(fixed_root_side(&[q(1), q(2), q(1)]), RootSide::Negative);
        assert_eq!(fixed_root_side(&[q(-3)]), RootSide::Negative);
        assert_eq!(fixed_root_side(&[q(0), q(1), q(1)]), RootSide::NonPositive);
        assert_eq!(fixed_root_side(&[q(-1), q(0), q(1)]), RootSide::Other);
    }
}
