//! Diagonal coefficient maps, derivative antecedents and the transport of
//! implicit relations through them.
//!
//! # Diagonal maps
//!
//! [`DiagonalWeight`] describes `Q_n = T_n P_n` with
//! `T_n(x^k) = lambda(n,k) x^k`, where
//!
//! ```text
//! lambda(n,k) = scale (k!)^p prod_j f_j(n,k)^(e_j)
//! ```
//!
//! for a nonzero rational `scale`, an integer `p` and polynomials `f_j`. The
//! presets are `factorial_power` (`1/(k!)^m`, the map of
//! [`RecurrenceCoefficientRescaling`], whose row arithmetic is reused),
//! `binomial` (`binom(k+m,m)`, i.e. `Q = W_m P` with `W_m = D^m x^m/m!`) and
//! `inverse_binomial` (`1/binom(k+m,m)`, i.e. `P = W_m Q`). Every weight that
//! a row needs must be defined and nonzero; rows are checked exactly, and the
//! uniform statement for all indices is certified or reported as unverified.
//!
//! # Diagonal transport
//!
//! Write `theta = x D` and multiply a relation `L` by `x^J` (`J` its largest
//! derivative order) so that `x^J L = sum x^s a_(i,s)(n,theta) S^i` with
//! `s >= 0` (using `x^d D^d = theta (theta - 1) ... (theta - d + 1)`). On each
//! monomial,
//!
//! ```text
//! T (x^s a(n,theta) S^i) T^-1 = x^s a(n,theta) lambda(n,theta+s)/lambda(n+i,theta) S^i,
//! ```
//!
//! and the ratio is a quotient of products of `(theta + t)` and shifted `f_j`
//! (identical factors cancel). If `b_(i,s)(n,theta)` is the denominator of
//! column `(i,s)`, one common multiplier
//! `H(n,theta) = lcm_(i,s) b_(i,s)(n, theta - s)` clears every column, since
//! `H(theta) x^s = x^s H(theta + s)`. The same multiplier is used for all
//! columns: clearing columns independently would change the kernel of an
//! ansatz ([`DiagonalWeight::clearing_multiplier`]). The cleared operator is
//! brought back to normal order (`theta^m = sum_d S(m,d) x^d D^d`) and a
//! common power of `x` is removed.
//!
//! For the residual `Z_n = T_n (x^J L P)_n`, `Z_n = 0` holds exactly when
//! `H(n,theta) Z_n = 0` and `[x^k] Z_n = 0` at every `k >= 0` with
//! `H(n,k) = 0`. Factors of `H` free of `theta` are nonzero by the weight
//! hypothesis. For the others the exceptional indices must be certified
//! uniformly in `n` (an index-free factor has the exact set of its
//! nonnegative integer roots; a one-signed shifted factor has none); the
//! coefficient conditions `[x^k] Z_n = 0` are then recorded in source terms,
//! `sum_(i,s) a_(i,s)(n, k-s) [x^(k-s)] P_(n+i) = 0`, rewritten in `Q`. An
//! uncertified factor stays visible in [`ExceptionalIndexRecord`]; it is never
//! replaced by testing a prefix.
//!
//! # Derivative antecedents
//!
//! For a fixed nonzero rational `g` of degree `m`, `T_g Q = D^m(g Q)/m!` has
//! `[x^e] T_g Q = binom(e+m,m) sum_j g_j q_(e+m-j)`: its monomial matrix is
//! triangular with diagonal `lc(g) binom(e+m,m) != 0`, so antecedents are
//! recovered by exact back substitution from the top degree. A relation `L`
//! of `P = T_g Q` becomes `L T_g Q = 0` by substitution, raising the
//! derivative order by at most `m`; `g = x^m` gives `W_m`.
//!
//! Related modules: `implicit::relation` (relations, descriptions,
//! certificates), `fit::rescaling` (factorial rescaling of complete spaces).

use super::algebra::{
    analyze_zeros, binomial, biv_text, factorial, from_theta, poly_text, q_usize, theta_add,
    to_theta, Biv, DOp, ThetaOp,
};
use super::relation::{biv_from_json, biv_json};
use super::*;
use crate::polynomial::Polynomial;
use crate::vec_poly;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

type Q = BigRational;

/// Largest `m`, `power`, `|p|` or `|e_j|` accepted in a weight or antecedent.
pub const MAX_IMPLICIT_MAP_ORDER: usize = 64;

/// One factor `f(n,k)^exponent` of a [`DiagonalWeight::Product`].
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct DiagonalWeightFactor {
    /// `f(n,k)`: `coeffs[a][b]` is the coefficient of `n^a k^b`.
    pub coeffs: BivarPolyJson,
    /// Nonzero integer exponent `e`.
    pub exponent: i64,
}

/// The weight `lambda(n,k)` of a diagonal map `T_n(x^k) = lambda(n,k) x^k`.
///
/// The JSON form is tagged by `kind`.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum DiagonalWeight {
    /// `lambda = 1/(k!)^power`: coefficient `k` is divided by `(k!)^power`.
    FactorialPower {
        /// The exponent `m >= 1`.
        power: usize,
    },
    /// `lambda = binom(k+m, m)`: `Q_n = W_m P_n` with `W_m = D^m x^m/m!`.
    Binomial {
        /// The order `m >= 1`.
        m: usize,
    },
    /// `lambda = 1/binom(k+m, m)`: `P_n = W_m Q_n`.
    InverseBinomial {
        /// The order `m >= 1`.
        m: usize,
    },
    /// `lambda = scale (k!)^factorial_exponent prod_j f_j(n,k)^(e_j)`.
    Product {
        /// Nonzero rational `scale`, as a string.
        scale: String,
        /// The exponent `p` of `k!`.
        #[serde(default)]
        factorial_exponent: i64,
        /// The polynomial factors.
        #[serde(default)]
        factors: Vec<DiagonalWeightFactor>,
    },
}

/// `scale (k!)^p prod f_j(n,k)^(e_j)` with polynomials in `(n, k)`.
#[derive(Debug, Clone)]
pub(super) struct WeightForm {
    scale: Q,
    factorial: i64,
    factors: Vec<(Biv, i64)>,
}

fn order_check(value: usize, name: &str) -> Result<(), TransformError> {
    if value == 0 || value > MAX_IMPLICIT_MAP_ORDER {
        return Err(TransformError::InvalidParameter(format!(
            "{name} must lie in 1..={MAX_IMPLICIT_MAP_ORDER}"
        )));
    }
    Ok(())
}

impl DiagonalWeight {
    pub(super) fn form(&self) -> Result<WeightForm, TransformError> {
        let linear = |t: usize| Biv::from_y(&[q_usize(t), Q::one()]);
        Ok(match self {
            Self::FactorialPower { power } => {
                order_check(*power, "power")?;
                WeightForm {
                    scale: Q::one(),
                    factorial: -(*power as i64),
                    factors: Vec::new(),
                }
            }
            Self::Binomial { m } => {
                order_check(*m, "m")?;
                WeightForm {
                    scale: factorial(*m).recip(),
                    factorial: 0,
                    factors: (1..=*m).map(|t| (linear(t), 1)).collect(),
                }
            }
            Self::InverseBinomial { m } => {
                order_check(*m, "m")?;
                WeightForm {
                    scale: factorial(*m),
                    factorial: 0,
                    factors: (1..=*m).map(|t| (linear(t), -1)).collect(),
                }
            }
            Self::Product {
                scale,
                factorial_exponent,
                factors,
            } => {
                let scale = parse_rational_coeff(scale)
                    .map_err(|e| TransformError::InvalidParameter(format!("scale: {e}")))?;
                if scale.is_zero() {
                    return Err(TransformError::InvalidParameter(
                        "the weight scale must be nonzero".into(),
                    ));
                }
                if factorial_exponent.unsigned_abs() as usize > MAX_IMPLICIT_MAP_ORDER {
                    return Err(TransformError::InvalidParameter(
                        "factorial_exponent is too large".into(),
                    ));
                }
                if factors.len() > 16 {
                    return Err(TransformError::InvalidParameter(
                        "at most 16 weight factors".into(),
                    ));
                }
                let mut parsed = Vec::new();
                for factor in factors {
                    let f = biv_from_json(&factor.coeffs)
                        .map_err(|e| TransformError::InvalidParameter(format!("factor: {e}")))?;
                    if f.is_zero() {
                        return Err(TransformError::InvalidParameter(
                            "a weight factor is the zero polynomial".into(),
                        ));
                    }
                    if factor.exponent == 0
                        || factor.exponent.unsigned_abs() as usize > MAX_IMPLICIT_MAP_ORDER
                    {
                        return Err(TransformError::InvalidParameter(
                            "factor exponents must be nonzero and small".into(),
                        ));
                    }
                    parsed.push((f, factor.exponent));
                }
                WeightForm {
                    scale,
                    factorial: *factorial_exponent,
                    factors: parsed,
                }
            }
        })
    }

    /// Check the parameters (orders, nonzero scale and factors).
    pub fn validate(&self) -> Result<(), TransformError> {
        self.form().map(|_| ())
    }

    /// The exact weight `lambda(n,k)`; an error when it is zero or undefined.
    pub fn weight(&self, n: usize, k: usize) -> Result<Q, TransformError> {
        self.form()?.weight(n, k)
    }

    /// Human-readable formula of `lambda(n,k)`.
    pub fn formula(&self) -> String {
        match self {
            Self::FactorialPower { power } => format!("1/(k!)^{power}"),
            Self::Binomial { m } => format!("binom(k+{m},{m})"),
            Self::InverseBinomial { m } => format!("1/binom(k+{m},{m})"),
            Self::Product { .. } => match self.form() {
                Ok(form) => {
                    let mut parts = vec![form.scale.to_string()];
                    if form.factorial != 0 {
                        parts.push(format!("(k!)^({})", form.factorial));
                    }
                    for (f, e) in &form.factors {
                        parts.push(format!("({})^({e})", biv_text(f, "k")));
                    }
                    parts.join(" ")
                }
                Err(_) => "invalid weight".into(),
            },
        }
    }

    /// `Q_n = T_n P_n` on one row; every weight up to the row degree must be
    /// defined and nonzero.
    pub fn map_row(&self, index: usize, row: &[Q]) -> Result<Vec<Q>, TransformError> {
        let row = vec_poly::trim(row);
        if let Self::FactorialPower { power } = self {
            order_check(*power, "power")?;
            // The factorial rescaling of complete spaces does the same arithmetic.
            return RecurrenceCoefficientRescaling::FactorialPower { power: *power }
                .rescale_coefficients(&row, &RecurrenceSpaceLimits::default())
                .map(|r| vec_poly::trim(&r))
                .map_err(TransformError::Limit);
        }
        let form = self.form()?;
        let mapped = row
            .iter()
            .enumerate()
            .map(|(k, c)| Ok(c * form.weight(index, k)?))
            .collect::<Result<Vec<_>, TransformError>>()?;
        Ok(vec_poly::trim(&mapped))
    }

    /// `P_n = T_n^-1 Q_n` on one row.
    pub fn inverse_row(&self, index: usize, row: &[Q]) -> Result<Vec<Q>, TransformError> {
        let row = vec_poly::trim(row);
        if let Self::FactorialPower { power } = self {
            order_check(*power, "power")?;
            return RecurrenceCoefficientRescaling::FactorialPower { power: *power }
                .unscale_coefficients(&row, &RecurrenceSpaceLimits::default())
                .map(|r| vec_poly::trim(&r))
                .map_err(TransformError::Limit);
        }
        let form = self.form()?;
        let mapped = row
            .iter()
            .enumerate()
            .map(|(k, c)| Ok(c / form.weight(index, k)?))
            .collect::<Result<Vec<_>, TransformError>>()?;
        Ok(vec_poly::trim(&mapped))
    }

    /// The common clearing multiplier `H(n,theta)` for an ansatz whose
    /// Euler-form columns are `x^s (.) S^i` for `(i, s)` in `columns`, as a
    /// list of monic factors (with repetition), each in `(n, theta)`.
    pub fn clearing_multiplier(
        &self,
        columns: &[(usize, usize)],
    ) -> Result<Vec<BivarPolyJson>, TransformError> {
        let form = self.form()?;
        let columns = columns.iter().copied().collect::<BTreeSet<_>>();
        Ok(multiset_list(&clearing_factors(&form, &columns))
            .iter()
            .map(biv_json)
            .collect())
    }

    /// `T` with `P_n = T Q_n` as a shift-free differential operator, when
    /// `1/lambda` is a polynomial in `(n, k)` (`p = 0`, all `e_j < 0`).
    pub fn inverse_operator(&self) -> Option<ImplicitRelation> {
        let form = self.form().ok()?;
        if form.factorial != 0 || form.factors.iter().any(|(_, e)| *e > 0) {
            return None;
        }
        let polynomial = form
            .factors
            .iter()
            .fold(Biv::constant(form.scale.recip()), |acc, (f, e)| {
                acc.mul(&f.pow(e.unsigned_abs() as usize))
            });
        let mut theta = ThetaOp::new();
        theta_add(&mut theta, (0, 0), polynomial);
        Some(ImplicitRelation::from_op(&from_theta(&theta)))
    }
}

impl WeightForm {
    fn weight(&self, n: usize, k: usize) -> Result<Q, TransformError> {
        let undefined = || {
            TransformError::InvalidParameter(format!(
                "diagonal weight lambda({n},{k}) is zero or undefined"
            ))
        };
        let mut value = self.scale.clone();
        let kf = factorial(k);
        for _ in 0..self.factorial.unsigned_abs() {
            if self.factorial > 0 {
                value *= &kf;
            } else {
                value /= &kf;
            }
        }
        let n_value = q_usize(n);
        let k_value = q_usize(k);
        for (f, e) in &self.factors {
            let v = vec_poly::evaluate(&f.eval_n(&n_value), &k_value);
            if v.is_zero() {
                return Err(undefined());
            }
            for _ in 0..e.unsigned_abs() {
                if *e > 0 {
                    value *= &v;
                } else {
                    value /= &v;
                }
            }
        }
        if value.is_zero() {
            return Err(undefined());
        }
        Ok(value)
    }
}

// ---------------------------------------------------------------------------
// Ratios and the clearing multiplier
// ---------------------------------------------------------------------------

type Multiset = BTreeMap<Biv, usize>;

fn multiset_push(set: &mut Multiset, factor: &Biv, count: usize, scalar: &mut Q, numerator: bool) {
    let (lead, monic) = factor.monic();
    for _ in 0..count {
        if numerator {
            *scalar *= &lead;
        } else {
            *scalar /= &lead;
        }
    }
    *set.entry(monic).or_default() += count;
}

fn multiset_list(set: &Multiset) -> Vec<Biv> {
    set.iter()
        .flat_map(|(f, c)| std::iter::repeat_n(f.clone(), *c))
        .collect()
}

fn multiset_product(set: &Multiset) -> Biv {
    set.iter()
        .fold(Biv::one(), |acc, (f, c)| acc.mul(&f.pow(*c)))
}

/// `lambda(n,theta+s)/lambda(n+i,theta) = scalar * prod num / prod den`
/// with monic factors in `(n, theta)` and common factors cancelled.
struct WeightRatio {
    scalar: Q,
    numerator: Multiset,
    denominator: Multiset,
}

fn weight_ratio(form: &WeightForm, i: usize, s: usize) -> WeightRatio {
    let mut scalar = Q::one();
    let mut numerator = Multiset::new();
    let mut denominator = Multiset::new();
    let p = form.factorial.unsigned_abs() as usize;
    for t in 1..=s {
        let factor = Biv::from_y(&[q_usize(t), Q::one()]);
        if form.factorial > 0 {
            multiset_push(&mut numerator, &factor, p, &mut scalar, true);
        } else if form.factorial < 0 {
            multiset_push(&mut denominator, &factor, p, &mut scalar, false);
        }
    }
    for (f, e) in &form.factors {
        let top = f.shift(&Q::zero(), &q_usize(s));
        let bottom = f.shift(&q_usize(i), &Q::zero());
        let count = e.unsigned_abs() as usize;
        let (up, down) = if *e > 0 { (top, bottom) } else { (bottom, top) };
        multiset_push(&mut numerator, &up, count, &mut scalar, true);
        multiset_push(&mut denominator, &down, count, &mut scalar, false);
    }
    let common = numerator
        .iter()
        .filter_map(|(f, c)| denominator.get(f).map(|d| (f.clone(), *c.min(d))))
        .collect::<Vec<_>>();
    for (f, c) in common {
        for set in [&mut numerator, &mut denominator] {
            let entry = set.get_mut(&f).expect("common factor");
            *entry -= c;
            if *entry == 0 {
                set.remove(&f);
            }
        }
    }
    WeightRatio {
        scalar,
        numerator,
        denominator,
    }
}

/// `H = lcm_(i,s) b_(i,s)(n, theta - s)` over the given columns.
fn clearing_factors(form: &WeightForm, columns: &BTreeSet<(usize, usize)>) -> Multiset {
    let mut lcm = Multiset::new();
    for &(i, s) in columns {
        let ratio = weight_ratio(form, i, s);
        let mut shifted = Multiset::new();
        for (f, c) in &ratio.denominator {
            let g = f.shift(&Q::zero(), &-q_usize(s)).monic().1;
            *shifted.entry(g).or_default() += c;
        }
        for (g, c) in shifted {
            let entry = lcm.entry(g).or_default();
            *entry = (*entry).max(c);
        }
    }
    lcm
}

fn shift_multiset(set: &Multiset, s: usize) -> Multiset {
    let mut out = Multiset::new();
    for (f, c) in set {
        *out.entry(f.shift(&Q::zero(), &q_usize(s)).monic().1)
            .or_default() += c;
    }
    out
}

// ---------------------------------------------------------------------------
// Certificates
// ---------------------------------------------------------------------------

/// One Euler-form term `x^x_power a(n,theta) S^shift`.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct EulerTermJson {
    /// The shift `i`.
    pub shift: usize,
    /// The power `s` of `x`.
    pub x_power: usize,
    /// `a(n,theta)`: `coeffs[a][b]` multiplies `n^a theta^b`.
    pub coeff: BivarPolyJson,
}

/// The ratio `lambda(n,theta+s)/lambda(n+i,theta)` of one column.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct EulerRatioJson {
    /// The shift `i`.
    pub shift: usize,
    /// The power `s` of `x`.
    pub x_power: usize,
    /// Constant factor.
    pub scalar: String,
    /// Monic numerator factors in `(n, theta)`, with repetition.
    pub numerator: Vec<BivarPolyJson>,
    /// Monic denominator factors `b_(i,s)` in `(n, theta)`, with repetition.
    pub denominator: Vec<BivarPolyJson>,
}

/// Exceptional coefficient indices of a clearing multiplier.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ExceptionalIndexRecord {
    /// `certified_all_indices` when every factor's zero set is certified
    /// uniformly in `n`, otherwise `checked_supplied_rows` (the certified
    /// part is listed; the rest is finite evidence or nothing).
    pub status: String,
    /// Certified exceptional indices `k` (zeros of `H(n,k)` for every `n`).
    pub indices: Vec<usize>,
    /// Exceptional indices whose coefficient condition holds vacuously.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub vacuous_indices: Vec<usize>,
    /// One analysis per distinct factor of `H`.
    pub factor_analyses: Vec<String>,
    /// Zeros `[n, k]` of uncertified factors at the supplied base indices and
    /// coefficient indices up to the row degree (finite evidence only).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub observed_zeros: Vec<[usize; 2]>,
}

/// How a relation was transported through a diagonal map and cleared.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ClearingCertificate {
    /// `J`: the source relation was multiplied by `x^J` first.
    pub premultiplied_x_power: usize,
    /// `x^J L` in Euler form (source coordinates).
    pub euler_form: Vec<EulerTermJson>,
    /// The ratio of every column.
    pub ratios: Vec<EulerRatioJson>,
    /// Factors of the common multiplier `H(n,theta)`, with repetition.
    pub multiplier: Vec<BivarPolyJson>,
    /// `H(n,theta)` as text.
    pub multiplier_text: String,
    /// `e` in `H T x^J L T^-1 = c x^e M`.
    pub removed_x_power: usize,
    /// `c` in `H T x^J L T^-1 = c x^e M`.
    pub removed_scalar: String,
    /// The exceptional indices and their certificate.
    pub exceptional: ExceptionalIndexRecord,
    /// What the clearing proves.
    pub statement: String,
}

/// A relation transported into the new coordinates.
#[derive(Debug, Clone)]
pub struct ImplicitTransport {
    /// The transported relation `M` (canonical form).
    pub relation: ImplicitRelation,
    /// Non-vacuous coefficient conditions at certified exceptional indices `k`.
    pub exceptional_conditions: Vec<(usize, CoefficientCondition)>,
    /// The clearing data of a diagonal map.
    pub clearing: Option<ClearingCertificate>,
    /// `L o T = sum H_nu S^nu M`, when `T` is a differential operator.
    pub reconstruction: Option<ReconstructionCertificate>,
}

// ---------------------------------------------------------------------------
// Derivative antecedents
// ---------------------------------------------------------------------------

/// Parse `g` (ascending rational strings) and return it with `m = deg g`.
pub fn derivative_antecedent_parameters(g: &[String]) -> Result<(Vec<Q>, usize), TransformError> {
    let g = g
        .iter()
        .map(|c| parse_rational_coeff(c))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| TransformError::InvalidParameter(format!("g: {e}")))?;
    let g = vec_poly::trim(&g);
    let Some(m) = vec_poly::degree(&g) else {
        return Err(TransformError::InvalidParameter("g must be nonzero".into()));
    };
    if m > MAX_IMPLICIT_MAP_ORDER {
        return Err(TransformError::InvalidParameter(format!(
            "deg g must be at most {MAX_IMPLICIT_MAP_ORDER}"
        )));
    }
    Ok((g, m))
}

/// `T_g Q = D^m (g Q)/m!` with `m = deg g`.
pub fn derivative_lift_row(g: &[Q], row: &[Q]) -> Vec<Q> {
    let g = vec_poly::trim(g);
    let m = vec_poly::degree(&g).unwrap_or(0);
    let mut product = vec_poly::mul(&g, row);
    for _ in 0..m {
        product = vec_poly::derivative(&product);
    }
    vec_poly::scale(&product, &factorial(m).recip())
}

/// The unique `Q` with `D^m (g Q)/m! = P`, by back substitution from the top
/// degree: `q_e = (p_e/binom(e+m,m) - sum_(j<m) g_j q_(e+m-j))/g_m`.
pub fn derivative_antecedent_row(g: &[Q], row: &[Q]) -> Vec<Q> {
    let g = vec_poly::trim(g);
    let row = vec_poly::trim(row);
    let m = vec_poly::degree(&g).unwrap_or(0);
    let Some(d) = vec_poly::degree(&row) else {
        return vec![Q::zero()];
    };
    let mut out = vec![Q::zero(); d + 1];
    for e in (0..=d).rev() {
        let mut value = &row[e] / binomial(e + m, m);
        for (j, gj) in g.iter().enumerate().take(m) {
            if let Some(qk) = out.get(e + m - j) {
                value -= gj * qk;
            }
        }
        out[e] = value / &g[m];
    }
    vec_poly::trim(&out)
}

fn antecedent_operator(g: &[Q]) -> DOp {
    let m = vec_poly::degree(g).unwrap_or(0);
    DOp::single(0, m, Biv::one())
        .compose(&DOp::single(0, 0, Biv::from_y(g)))
        .scale(&factorial(m).recip())
}

// ---------------------------------------------------------------------------
// Transport
// ---------------------------------------------------------------------------

impl RecurrenceTransformation {
    /// The shift-free operator `T` with `P_n = T Q_n`, when it is a
    /// differential operator with polynomial coefficients: `T_g` for a
    /// derivative antecedent, `1/lambda(n,theta)` for a diagonal map whose
    /// inverse weight is a polynomial (for example `W_m`).
    pub fn reconstruction_operator(&self) -> Option<ImplicitRelation> {
        match self {
            Self::Diagonal { weight } => weight.inverse_operator(),
            Self::DerivativeAntecedent { g } => {
                let (g, _) = derivative_antecedent_parameters(g).ok()?;
                Some(ImplicitRelation::from_op(&antecedent_operator(&g)))
            }
            _ => None,
        }
    }

    /// Transport an implicit relation of `P`, valid for base indices
    /// `n >= first_base`, into the new coordinates.
    ///
    /// Supported for `diagonal` (pullback, one common clearing multiplier,
    /// exceptional-index conditions) and `derivative_antecedent`
    /// (substitution `L T_g Q = 0`). The result is an exact consequence of
    /// the source relation; with the exceptional-index conditions it is
    /// equivalent to it when the exceptional set is certified. The lossy map
    /// `derivative` gives `M = D^(m+r) L D^(-m)` with a recorded identity
    /// `M o D^m = D^(m+r) o L` (a consequence only; see module `implicit::lossy`).
    ///
    /// ```
    /// use polytool::recurrence::*;
    ///
    /// // (1 + x)^n under coefficient division by k!: Q'_(n+1) = Q'_n + Q_n.
    /// let rule = polytool::oeis::by_id("A007318").unwrap().recurrence_export_parts().unwrap().0;
    /// let (relation, _) = ImplicitRelation::from_recurrence(&rule).unwrap();
    /// let map = RecurrenceTransformation::Diagonal {
    ///     weight: DiagonalWeight::FactorialPower { power: 1 },
    /// };
    /// let transported = map.transport_implicit(&relation, 0).unwrap();
    /// assert_eq!(transported.relation.text("Q"), "D Q[n+1] - D Q[n] - Q[n] = 0");
    /// // The lost constant: [x^0] (Q_(n+1) - Q_n) = 0.
    /// assert_eq!(transported.exceptional_conditions.len(), 1);
    /// ```
    pub fn transport_implicit(
        &self,
        relation: &ImplicitRelation,
        first_base: usize,
    ) -> Result<ImplicitTransport, TransformError> {
        self.validate()?;
        match self {
            Self::Diagonal { weight } => transport_diagonal(self, weight, relation, first_base),
            Self::DerivativeAntecedent { g } => {
                let (g, _) = derivative_antecedent_parameters(g)?;
                let t = antecedent_operator(&g);
                let product = relation.op().compose(&t);
                let (_, _, normal) = product.normalize();
                let target = ImplicitRelation::from_op(&normal);
                let map = ImplicitRelation::from_op(&t);
                let reconstruction = find_reconstruction_certificate(relation, &map, &target);
                Ok(ImplicitTransport {
                    relation: target,
                    exceptional_conditions: Vec::new(),
                    clearing: None,
                    reconstruction,
                })
            }
            Self::Derivative { order } => transport_derivative(relation, *order),
            _ => Err(TransformError::TransportUnsupported(format!(
                "implicit transport supports diagonal, derivative_antecedent and derivative \
                 maps; use transport for {}",
                self.kind_name()
            ))),
        }
    }
}

fn transport_diagonal(
    transformation: &RecurrenceTransformation,
    weight: &DiagonalWeight,
    relation: &ImplicitRelation,
    first_base: usize,
) -> Result<ImplicitTransport, TransformError> {
    let form = weight.form()?;
    let source = relation.op();
    if source.is_zero() {
        return Err(TransformError::TransportUnsupported(
            "the zero relation".into(),
        ));
    }
    let (j, theta) = to_theta(&source);
    let columns = theta.keys().copied().collect::<BTreeSet<_>>();
    let h = clearing_factors(&form, &columns);
    let mut cleared = ThetaOp::new();
    let mut ratios = Vec::new();
    for (&(i, s), a) in &theta {
        let ratio = weight_ratio(&form, i, s);
        let mut quotient = shift_multiset(&h, s);
        for (f, c) in &ratio.denominator {
            let entry = quotient.get_mut(f).ok_or_else(|| {
                TransformError::Limit("internal: clearing multiplier misses a factor".into())
            })?;
            *entry = entry.checked_sub(*c).ok_or_else(|| {
                TransformError::Limit("internal: clearing multiplier too small".into())
            })?;
        }
        let coefficient = a
            .scale(&ratio.scalar)
            .mul(&multiset_product(&ratio.numerator))
            .mul(&multiset_product(&quotient));
        theta_add(&mut cleared, (i, s), coefficient);
        ratios.push(EulerRatioJson {
            shift: i,
            x_power: s,
            scalar: format_rational_coeff(&ratio.scalar),
            numerator: multiset_list(&ratio.numerator)
                .iter()
                .map(biv_json)
                .collect(),
            denominator: multiset_list(&ratio.denominator)
                .iter()
                .map(biv_json)
                .collect(),
        });
    }
    let cleared_op = from_theta(&cleared);
    let (lead, e, normal) = cleared_op.normalize();
    let target = ImplicitRelation::from_op(&normal);

    // Exceptional indices of the theta-dependent factors of H.
    let mut certified = true;
    let mut indices = BTreeSet::new();
    let mut analyses = Vec::new();
    for f in h.keys() {
        if f.is_y_free() {
            analyses.push(format!(
                "{}: free of theta, nonzero by the weight hypothesis",
                biv_text(f, "theta")
            ));
            continue;
        }
        let analysis = analyze_zeros(f, first_base, None);
        if !analysis.singular_indices.is_empty() {
            certified = false;
            analyses.push(format!(
                "{}: vanishes identically at base n = {:?}, where every coefficient index is \
                 exceptional; {}",
                biv_text(f, "theta"),
                analysis.singular_indices,
                analysis.certificate
            ));
            continue;
        }
        if !analysis.certified {
            certified = false;
        } else {
            indices.extend(analysis.uniform_indices.iter().copied());
        }
        analyses.push(format!(
            "{}: {}",
            biv_text(f, "theta"),
            analysis.certificate
        ));
    }
    let mut conditions = Vec::new();
    let mut vacuous = Vec::new();
    for &k in &indices {
        let condition = exceptional_condition(&form, &theta, k);
        if condition.is_vacuous() {
            vacuous.push(k);
        } else {
            conditions.push((k, condition));
        }
    }
    let status = if certified {
        CERTIFIED_ALL_INDICES
    } else {
        CHECKED_SUPPLIED_ROWS
    };
    let multiplier_list = multiset_list(&h);
    let multiplier_text = if multiplier_list.is_empty() {
        "1".to_string()
    } else {
        multiplier_list
            .iter()
            .map(|f| format!("({})", biv_text(f, "theta")))
            .collect::<Vec<_>>()
            .join(" ")
    };
    let statement = format!(
        "With Z_n = T_n (x^{j} L P)_n and lambda(n,k) = {}: H(n,theta) Z_n = {} x^{e} (M Q)_n \
         exactly. Z_n = 0 iff H(n,theta) Z_n = 0 and [x^k] Z_n = 0 at every k >= 0 with \
         H(n,k) = 0; {} for base indices n >= {first_base}.",
        weight.formula(),
        format_rational_coeff(&lead),
        if certified {
            format!("the exceptional indices are exactly {:?}", indices)
        } else {
            "the exceptional set is NOT certified uniformly; the equivalence of M with the \
             source relation is not claimed"
                .to_string()
        }
    );
    let clearing = ClearingCertificate {
        premultiplied_x_power: j,
        euler_form: theta
            .iter()
            .map(|(&(i, s), a)| EulerTermJson {
                shift: i,
                x_power: s,
                coeff: biv_json(a),
            })
            .collect(),
        ratios,
        multiplier: multiplier_list.iter().map(biv_json).collect(),
        multiplier_text,
        removed_x_power: e,
        removed_scalar: format_rational_coeff(&lead),
        exceptional: ExceptionalIndexRecord {
            status: status.to_string(),
            indices: indices.into_iter().collect(),
            vacuous_indices: vacuous,
            factor_analyses: analyses,
            observed_zeros: Vec::new(),
        },
        statement,
    };
    let reconstruction = transformation
        .reconstruction_operator()
        .and_then(|map| find_reconstruction_certificate(relation, &map, &target));
    Ok(ImplicitTransport {
        relation: target,
        exceptional_conditions: conditions,
        clearing: Some(clearing),
        reconstruction,
    })
}

/// `[x^k] Z_n = 0` rewritten as `sum_(i,s<=k) a_(i,s)(n,k-s) [x^(k-s)] P_(n+i) = 0`
/// in the coordinates of `Q`, with `[x^(k')] P_(n+i) = [x^(k')] Q_(n+i)/lambda(n+i,k')`,
/// multiplied by the lcm of the index-polynomial denominators (nonzero on the
/// domain by the weight hypothesis).
fn exceptional_condition(form: &WeightForm, theta: &ThetaOp, k: usize) -> CoefficientCondition {
    // Each term: (shift, index, numerator(n), denominator(n)).
    let mut raw = Vec::new();
    for (&(i, s), a) in theta {
        if s > k {
            continue;
        }
        let index = k - s;
        let k_value = q_usize(index);
        let a_n = a.eval_y(&k_value);
        if vec_poly::is_zero(&a_n) {
            continue;
        }
        // 1/lambda(n+i, index) = scale^-1 (index!)^-p prod f(n+i, index)^-e.
        let mut constant = form.scale.recip();
        let kf = factorial(index);
        for _ in 0..form.factorial.unsigned_abs() {
            if form.factorial > 0 {
                constant /= &kf;
            } else {
                constant *= &kf;
            }
        }
        let mut numerator = vec_poly::scale(&a_n, &constant);
        let mut denominator = vec![Q::one()];
        for (f, e) in &form.factors {
            let value = f.shift(&q_usize(i), &Q::zero()).eval_y(&k_value);
            let power = vec_poly::pow(&value, e.unsigned_abs() as usize);
            if *e > 0 {
                denominator = vec_poly::mul(&denominator, &power);
            } else {
                numerator = vec_poly::mul(&numerator, &power);
            }
        }
        raw.push((i, index, numerator, denominator));
    }
    let lcm = raw
        .iter()
        .fold(Polynomial::new(vec![Q::one()]), |acc, (_, _, _, d)| {
            let d = Polynomial::new(d.clone());
            let g = acc.gcd(&d);
            let (quotient, _) = (acc.clone() * d).div_rem(&g);
            quotient
        });
    let terms = raw
        .into_iter()
        .map(|(shift, coefficient_index, numerator, denominator)| {
            let (factor, _) = lcm.div_rem(&Polynomial::new(denominator));
            ConditionTerm {
                shift,
                coefficient_index,
                coeff: vec_poly::trim(&vec_poly::mul(&numerator, factor.coeffs())),
            }
        })
        .collect();
    let condition = CoefficientCondition { terms }.canonical();
    // Scale the leading term's leading coefficient to one (a constant factor).
    let lead = condition
        .terms
        .last()
        .and_then(|t| vec_poly::degree(&t.coeff).map(|d| t.coeff[d].clone()));
    match lead {
        Some(lead) if !lead.is_zero() => CoefficientCondition {
            terms: condition
                .terms
                .into_iter()
                .map(|t| ConditionTerm {
                    coeff: vec_poly::scale(&t.coeff, &lead.recip()),
                    ..t
                })
                .collect(),
        },
        _ => condition,
    }
}

/// Weight-regularity statuses of a diagonal map on base indices `n >= first`:
/// every factor `f_j(n,k)` must be nonzero for `k >= 0`.
pub(in crate::recurrence) fn weight_regularity(
    weight: &DiagonalWeight,
    first: usize,
) -> Result<Vec<TransformConditionStatus>, TransformError> {
    let form = weight.form()?;
    let mut statuses = Vec::new();
    for (f, e) in &form.factors {
        let condition = format!(
            "lambda: ({})^({e}) defined and nonzero for n >= {first}, k >= 0",
            biv_text(f, "k")
        );
        let analysis = analyze_zeros(f, first, None);
        if !analysis.singular_indices.is_empty() || !analysis.uniform_indices.is_empty() {
            return Err(TransformError::InvalidParameter(format!(
                "the diagonal weight vanishes or is undefined: factor {} has zeros at \
                 n = {:?} (all k) and k = {:?} (all n)",
                biv_text(f, "k"),
                analysis.singular_indices,
                analysis.uniform_indices
            )));
        }
        statuses.push(TransformConditionStatus {
            condition,
            status: if analysis.certified {
                CERTIFIED_ALL_INDICES
            } else {
                CHECKED_SUPPLIED_ROWS
            }
            .to_string(),
            detail: if analysis.certified {
                analysis.certificate
            } else {
                format!(
                    "{}; the weight is checked exactly on the supplied rows only",
                    analysis.certificate
                )
            },
        });
    }
    if statuses.is_empty() {
        statuses.push(TransformConditionStatus {
            condition: "lambda(n,k) defined and nonzero".into(),
            status: CERTIFIED_ALL_INDICES.into(),
            detail: format!(
                "lambda = {} is a nonzero constant times a power of k!",
                weight.formula()
            ),
        });
    }
    Ok(statuses)
}

/// Whether the products of two lists of diagonal weights agree: the powers of
/// `k!` are equal and the rational parts agree as polynomials in `(n, k)`
/// after cross multiplication. Equal products give equal maps wherever both
/// are defined; `false` only means that equality is not certified.
pub(in crate::recurrence) fn diagonal_products_equal(
    left: &[DiagonalWeight],
    right: &[DiagonalWeight],
) -> Result<bool, TransformError> {
    let combine = |weights: &[DiagonalWeight]| -> Result<(i64, Biv, Biv), TransformError> {
        let mut factorial = 0i64;
        let mut numerator = Biv::one();
        let mut denominator = Biv::one();
        for weight in weights {
            let form = weight.form()?;
            numerator = numerator.scale(&form.scale);
            factorial = factorial
                .checked_add(form.factorial)
                .ok_or_else(|| TransformError::Limit("factorial exponent overflow".into()))?;
            for (f, e) in &form.factors {
                let power = f.pow(e.unsigned_abs() as usize);
                if *e > 0 {
                    numerator = numerator.mul(&power);
                } else {
                    denominator = denominator.mul(&power);
                }
            }
        }
        Ok((factorial, numerator, denominator))
    };
    let (left_factorial, left_numerator, left_denominator) = combine(left)?;
    let (right_factorial, right_numerator, right_denominator) = combine(right)?;
    Ok(left_factorial == right_factorial
        && left_numerator.mul(&right_denominator) == right_numerator.mul(&left_denominator))
}

/// Zeros `[n, k]` of the theta-dependent factors of a clearing multiplier
/// with `n` in `bases` and `0 <= k <= max_k` (finite evidence).
pub(super) fn observed_multiplier_zeros(
    clearing: &ClearingCertificate,
    bases: std::ops::Range<usize>,
    max_k: usize,
) -> Vec<[usize; 2]> {
    let factors = clearing
        .multiplier
        .iter()
        .filter_map(|f| biv_from_json(f).ok())
        .filter(|f| !f.is_y_free())
        .collect::<Vec<_>>();
    let mut zeros = Vec::new();
    for n in bases {
        for k in 0..=max_k {
            if factors
                .iter()
                .any(|f| vec_poly::evaluate(&f.eval_n(&q_usize(n)), &q_usize(k)).is_zero())
            {
                zeros.push([n, k]);
            }
        }
    }
    zeros
}

/// Record parts (row map, inverse, conditions, transport) of the maps of this module.
pub(in crate::recurrence) fn implicit_map_record(
    transformation: &RecurrenceTransformation,
) -> (String, String, Vec<String>, String) {
    match transformation {
        RecurrenceTransformation::Diagonal { weight } => (
            format!(
                "Q_n = T_n P_n, T_n(x^k) = lambda(n,k) x^k with lambda(n,k) = {}",
                weight.formula()
            ),
            "P_n = T_n^-1 Q_n: coefficient k divided by lambda(n,k)".to_string(),
            vec![
                "map and inverse: lambda(n,k) defined and nonzero for every row index n and \
                 k <= deg P_n"
                    .to_string(),
            ],
            "implicit: diagonal pullback x^s a(n,theta) S^i -> x^s a(n,theta) \
             lambda(n,theta+s)/lambda(n+i,theta) S^i, cleared by one common multiplier \
             H(n,theta) with exceptional-index conditions (implicit_transform_recurrence)"
                .to_string(),
        ),
        RecurrenceTransformation::DerivativeAntecedent { g } => {
            let text = derivative_antecedent_parameters(g)
                .map(|(g, m)| format!("g = {}, m = {m}", poly_text(&g, "x")))
                .unwrap_or_else(|e| e.to_string());
            (
                format!(
                    "Q_n = T_g^-1 P_n, the unique polynomial with D^m(g Q_n)/m! = P_n ({text}); \
                     exact back substitution from the top degree"
                ),
                "P_n = T_g Q_n = D^m(g Q_n)/m!".to_string(),
                Vec::new(),
                "implicit: substitution L T_g Q = 0 (derivative order grows by at most m) \
                 (implicit_transform_recurrence)"
                    .to_string(),
            )
        }
        _ => unreachable!("only implicit maps are recorded here"),
    }
}

/// Condition statuses of the maps of this module in a coordinate-change report.
pub(in crate::recurrence) fn implicit_map_condition_statuses(
    transformation: &RecurrenceTransformation,
    first: usize,
) -> Vec<TransformConditionStatus> {
    match transformation {
        RecurrenceTransformation::Diagonal { weight } => weight_regularity(weight, first)
            .unwrap_or_else(|error| {
                vec![TransformConditionStatus {
                    condition: "lambda(n,k) defined and nonzero".into(),
                    status: CHECKED_SUPPLIED_ROWS.into(),
                    detail: error.to_string(),
                }]
            }),
        _ => Vec::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::implicit::algebra::q;

    #[test]
    fn antecedent_back_substitution_inverts_the_lift() {
        let g = vec![q(2), q(-1), q(3)];
        let row = vec![q(5), q(0), q(-7), q(1)];
        let p = derivative_lift_row(&g, &row);
        assert_eq!(derivative_antecedent_row(&g, &p), vec_poly::trim(&row));
        // W_1 = D x: (x Q)' has diagonal k + 1.
        let w = vec![q(0), q(1)];
        assert_eq!(
            derivative_lift_row(&w, &[q(1), q(1), q(1)]),
            vec![q(1), q(2), q(3)]
        );
    }

    #[test]
    fn binomial_weights_and_inverse_operator() {
        let w = DiagonalWeight::Binomial { m: 2 };
        assert_eq!(w.weight(7, 3).unwrap(), q(10));
        let inverse = DiagonalWeight::InverseBinomial { m: 1 };
        let t = inverse.inverse_operator().unwrap();
        // P = W_1 Q = x Q' + Q.
        assert_eq!(
            t.apply_to_row(&[q(1), q(1), q(1)], 0),
            vec![q(1), q(2), q(3)]
        );
        assert!(w.inverse_operator().is_none());
        let negative = DiagonalWeight::Product {
            scale: "0".into(),
            factorial_exponent: 0,
            factors: Vec::new(),
        };
        assert!(negative.validate().is_err());
    }
}
