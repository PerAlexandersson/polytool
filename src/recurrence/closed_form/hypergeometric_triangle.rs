//! Certified recurrences for the row polynomials of hypergeometric
//! triangles.
//!
//! A *hypergeometric term* `T(n, k)` is given by its two shift ratios
//!
//! ```text
//! T(n, k + 1) / T(n, k) = rho_k(n, k),     T(n + 1, k) / T(n, k) = rho_n(n, k),
//! ```
//!
//! products of integer powers of linear forms in `(n, k)` with rational
//! coefficients (one pair per residue class of `n` modulo a period `p`, so
//! that `T(2m, k) = (m! C(m, k))^2`, `T(2m + 1, k) = ((m + 1)! C(m, k))^2 /
//! (k + 1)` is covered), a *support* `lower(n) <= k <= upper(n)` with
//! `lower(n) = ceil((a n + b) / d)` and `upper(n) = floor((a' n + b') / d')`
//! (`a, a' >= 0`), a first index `n_0` and the value `T(n_0, lower(n_0))`.
//! The table is defined by walking: along a row with `rho_k`, and from
//! `T(n, lower(n + 1))` to `T(n + 1, lower(n + 1))` with `rho_n`; outside the
//! support `T = 0`.  A [`HypergeometricTriangle`] is a finite sum of such
//! terms with a common first index, and its row polynomials are
//! `P_n(t) = sum_k T(n, k) t^k`, `n >= n_0`.
//!
//! [`derive_hypergeometric_row_recurrence`] finds an operator
//!
//! ```text
//! E = sum c_(i,j,l)(n) t^j theta^l S^i,     theta = t d/dt,
//! ```
//!
//! with `E P_n = 0` for every `n >= n_0` and a leader `c_(I,0,0)(n) S^I` free
//! of `t` and `theta`, and turns it into a runnable recurrence record
//! (`RecurrenceJson`, derivatives `t^r D^r` from
//! `theta^l = sum_r S(l, r) t^r D^r`).
//! [`certify_hypergeometric_row_operator`] proves `E P_n = 0` for a given
//! operator.  Both return a [`HypergeometricRowRecurrenceCertificate`] that
//! replays by recomputation.
//!
//! # The argument
//!
//! Since `theta` acts on `t^k` as multiplication by `k`,
//! `[t^k] t^j theta^l P_(n+i) = (k - j)^l T(n + i, k - j)`, so `E P_n = 0` is
//! the statement that
//!
//! ```text
//! e(n, k) = sum c_(i,j,l)(n) (k - j)^l T(n + i, k - j) = 0     for every k.
//! ```
//!
//! 1. **The definition is consistent.**  For each residue `r` modulo `p`,
//!    `rho_k(n + 1, k) rho_n(n, k) = rho_n(n, k + 1) rho_k(n, k)` as rational
//!    functions (exact comparison of reduced factorizations), and every
//!    linear factor of every ratio is nonzero at every lattice point where
//!    the walk evaluates it (`lower(n) <= k < upper(n)` for `rho_k`,
//!    `lower(n + 1) <= k <= upper(n)` for `rho_n`): it has one strict sign
//!    on the polyhedron (its two corners at the first index and its two
//!    recession directions), or it has no integer zero at all
//!    (`gcd(A p', B)` does not divide the constant).  Hence `T(n, k) != 0`
//!    on the support and, by induction along the rows,
//!    `T(n + 1, k) = rho_n(n, k) T(n, k)` wherever both points are in the
//!    support.  The supports are monotone and nonempty, and
//!    `lower(n + 1) <= upper(n)` (checked per residue class).
//! 2. **Interior.**  When the whole rectangle `[n, n + I] x [k - J, k]` lies
//!    in the support, every `T(n + i, k - j)` is `T(n, k)` times a product
//!    `R_(i,j)(n, k)` of ratio values along a path inside the support, so
//!    `e(n, k) = T(n, k) sum c (k - j)^l R_(i,j)` and the rational function
//!    `sum c (k - j)^l R_(i,j)` must vanish identically.  After clearing
//!    the common denominator (a product of the linear factors, nonzero on
//!    the support) this is an exact polynomial identity in `(n, k)`.
//! 3. **Boundary strips.**  The other columns with a nonzero term lie in
//!    finitely many strips `k = lower(n) + m` (`0 <= m < m_max`) and
//!    `k = upper(n) + m` (`1 <= m <= m'_max`), per residue class of `n`
//!    modulo the period of the supports.  In a strip the terms in the
//!    support are fixed (the others are zero), and each is
//!    `T(n, lower(n))` (resp. `T(n, upper(n))`) times a rational function
//!    of `n` along a path in the support, so `e(n, lower(n) + m) = 0`
//!    reduces to an exact identity of rational functions of `n`.  Boundary
//!    contributions are therefore never assumed to vanish: an operator
//!    whose strip identities fail is rejected, and the search imposes
//!    them as equations.
//! 4. **Small rows.**  The strip description holds from the separation
//!    index `N` on (strips disjoint, paths inside the support); for
//!    `n_0 <= n < N` the identity `E P_n = 0` is checked exactly on the
//!    rows.
//!
//! So `E P_n = 0` for every `n >= n_0`.  The record has the rows
//! `P_(n_0), ..., P_(n_0 + s - 1)` as seeds, with `s >= I` large enough that
//! the leader `c_(I,0,0)` does not vanish at a base where the rule computes
//! a row, so by induction it generates exactly the rows `P_n`.
//!
//! # The search
//!
//! The coefficients `c_(i,j,l)(n)` are polynomials of degree at most `d`
//! with unknown rational coefficients; the interior and strip identities
//! are linear in them.  Orders `I = 1, 2, ...`, then `theta`-degrees `L`,
//! `t`-degrees `J` and index degrees `d` are tried in that order (bounded
//! by [`HypergeometricRowRecurrenceBudget`]); the first exact kernel vector
//! with a nonzero leader is divided by the polynomial content of its
//! coefficients and certified as above.  The search is linear elimination
//! in the left ideal of operators annihilating the term (no fitting to
//! rows).
//!
//! Related: `parent_theorems::hypergeometric_row` (real-rootedness of rows
//! whose coefficient ratio factors linearly), and the equality certificates
//! (`certify_recurrence_equality`, `certify_normalized_recurrence_equality`)
//! that compare a derived record with another recurrence for all indices.

use super::*;
use crate::recurrence::implicit::algebra::{from_theta, theta_add, Biv, DOp, ThetaOp};
use crate::vec_poly;
use factored::{parse_factored, parse_linear, Factored, Linear};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

mod factored;
#[cfg(test)]
mod tests;

type Q = BigRational;

/// An operator in Euler form: `(shift, t power, theta power) -> c(n)`.
type EulerMap = BTreeMap<(usize, usize, usize), Vec<Q>>;

/// JSON schema tag of [`HypergeometricRowRecurrenceCertificate`].
pub const HYPERGEOMETRIC_ROW_RECURRENCE_SCHEMA: &str = "polytool.hypergeometric-row-recurrence.v1";

/// One hypergeometric summand `T(n, k)` of a [`HypergeometricTriangle`]
/// (module documentation).
///
/// Ratios are products of integer powers of linear forms in `n` and `k`
/// with rational constants, written as text: `(n-k)^2/((k+1)^2*(2n-2k-1))`,
/// `8*(n+1)/(n-2)`.  Entry `r` of `column_ratios` and `index_ratios` applies
/// to the rows `n` with `n mod p = r`, `p` the common length (the period).
/// The support bounds are linear in `n` with rational coefficients, rounded
/// up (`lower_column`) and down (`upper_column`): `"0"`, `"n"`, `"(n+1)/2"`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HypergeometricTermDefinition {
    /// The first row `n_0`.
    pub first_index: usize,
    /// `T(n, k + 1) / T(n, k)`, one entry per residue of `n` modulo the period.
    pub column_ratios: Vec<String>,
    /// `T(n + 1, k) / T(n, k)`, one entry per residue of `n` modulo the period.
    pub index_ratios: Vec<String>,
    /// `lower(n)`: the smallest `k` of row `n`, rounded up.
    pub lower_column: String,
    /// `upper(n)`: the largest `k` of row `n`, rounded down.
    pub upper_column: String,
    /// `h(n_0, lower(n_0))`, an exact nonzero rational, where `h` is the
    /// hypergeometric part (the term itself without `polynomial_factor`).
    pub first_value: String,
    /// An optional polynomial `M(n, k)` (text such as `2n*k + n - 2k^2`):
    /// the term is `T(n, k) = h(n, k) M(n, k)` with `h` given by the ratios.
    /// `M` may vanish anywhere; only `h` must be nonzero on the support.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub polynomial_factor: Option<String>,
}

impl HypergeometricTermDefinition {
    /// A term with one ratio pair (period one).
    pub fn new(
        first_index: usize,
        column_ratio: &str,
        index_ratio: &str,
        lower_column: &str,
        upper_column: &str,
        first_value: &str,
    ) -> Self {
        Self {
            first_index,
            column_ratios: vec![column_ratio.into()],
            index_ratios: vec![index_ratio.into()],
            lower_column: lower_column.into(),
            upper_column: upper_column.into(),
            first_value: first_value.into(),
            polynomial_factor: None,
        }
    }

    /// The same term times the polynomial `M(n, k)`.
    pub fn with_polynomial_factor(mut self, polynomial: &str) -> Self {
        self.polynomial_factor = Some(polynomial.into());
        self
    }
}

/// A finite sum of hypergeometric terms with a common first index: the
/// rows `P_n(t) = sum_k T(n, k) t^k` (module documentation).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HypergeometricTriangle {
    /// The summands, each with the same `first_index`.
    pub terms: Vec<HypergeometricTermDefinition>,
}

impl HypergeometricTriangle {
    /// The triangle of a single term.
    pub fn single(term: HypergeometricTermDefinition) -> Self {
        Self { terms: vec![term] }
    }

    /// The common first row.
    ///
    /// # Errors
    ///
    /// No terms, or terms with different first indices.
    pub fn first_index(&self) -> Result<usize, RecurrenceError> {
        let first = self
            .terms
            .first()
            .ok_or_else(|| RecurrenceError::InvalidInput("a triangle needs a term".into()))?
            .first_index;
        if self.terms.iter().any(|term| term.first_index != first) {
            return Err(RecurrenceError::InvalidInput(
                "every term needs the same first index".into(),
            ));
        }
        Ok(first)
    }

    /// The rows `P_(n_0), ..., P_(n_0 + count - 1)`, exactly (ascending
    /// coefficients, trimmed).  The definition is checked first (module
    /// documentation, step 1).
    ///
    /// # Errors
    ///
    /// A malformed or inconsistent definition, or a factor that may vanish on
    /// the support.
    pub fn generate_rows(&self, count: usize) -> Result<Vec<Vec<Q>>, RecurrenceError> {
        let terms = self.parse()?;
        total_rows(&terms, count)
    }

    fn parse(&self) -> Result<Vec<Term>, RecurrenceError> {
        self.first_index()?;
        self.terms
            .iter()
            .enumerate()
            .map(|(position, term)| {
                Term::parse(term)
                    .and_then(|parsed| {
                        parsed.check_definition()?;
                        Ok(parsed)
                    })
                    .map_err(|error| error.with_context(format!("term {position}")))
            })
            .collect()
    }
}

/// One term `c(n) t^variable_power theta^euler_power S^shift` of an
/// operator in the Euler form (`theta = t d/dt`, `S` the index shift).
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EulerOperatorTerm {
    /// The shift `i`: the term acts on `P_(n+i)`.
    pub shift: usize,
    /// The power `j` of `t`.
    pub variable_power: usize,
    /// The power `l` of `theta`.
    pub euler_power: usize,
    /// `c(n)`, ascending exact rational coefficients in `n`.
    pub coefficient: Vec<String>,
}

/// Bounds of [`derive_hypergeometric_row_recurrence`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HypergeometricRowRecurrenceBudget {
    /// Largest order `I` (shift of the leader).  Default 3.
    pub max_order: usize,
    /// Largest power `J` of `t`.  Default 4.
    pub max_variable_power: usize,
    /// Largest power `L` of `theta`.  Default 3.
    pub max_euler_power: usize,
    /// Largest degree `d` of a coefficient in `n`.  Default 6.
    pub max_index_degree: usize,
    /// Largest number of unknowns of one linear system.  Default 900.
    pub max_unknowns: usize,
    /// Largest row index evaluated exactly (rows before the separation
    /// index, seeds).  Default 400.
    pub max_checked_row: usize,
}

impl Default for HypergeometricRowRecurrenceBudget {
    fn default() -> Self {
        Self {
            max_order: 3,
            max_variable_power: 4,
            max_euler_power: 3,
            max_index_degree: 6,
            max_unknowns: 900,
            max_checked_row: 400,
        }
    }
}

/// The boundary strips of one term and residue class whose identities were
/// proved (module documentation, step 3).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HypergeometricBoundaryStrips {
    /// Position of the term in the triangle.
    pub term: usize,
    /// The residue of `n` modulo `modulus`.
    pub residue: usize,
    /// The period of the supports and ratios of the term.
    pub modulus: usize,
    /// Number of strips `k = lower(n) + m` (offsets `0..lower_strips`).
    pub lower_strips: usize,
    /// Number of strips `k = upper(n) + m` (offsets `1..=upper_strips`).
    pub upper_strips: usize,
    /// The separation index `N`: the strip description holds for `n >= N`.
    pub separation_index: usize,
}

/// A certified annihilating operator of the rows of a hypergeometric
/// triangle and the recurrence record it gives (module documentation),
/// schema `polytool.hypergeometric-row-recurrence.v1`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HypergeometricRowRecurrenceCertificate {
    /// [`HYPERGEOMETRIC_ROW_RECURRENCE_SCHEMA`].
    pub schema: String,
    /// The definition.
    pub triangle: HypergeometricTriangle,
    /// The operator `E`, in canonical order, primitive.
    pub operator: Vec<EulerOperatorTerm>,
    /// `E` as text.
    pub operator_text: String,
    /// The proved strip identities and separation indices.
    pub boundary_strips: Vec<HypergeometricBoundaryStrips>,
    /// Bases `n` at which `E P_n = 0` was checked exactly on rows.
    pub checked_bases: Vec<usize>,
    /// The runnable recurrence record defined by `E` and the seed rows.
    pub record: RecurrenceJson,
    /// What the certificate establishes.
    pub statement: String,
}

impl HypergeometricRowRecurrenceCertificate {
    /// Replay: recompute the certificate from the recorded triangle and
    /// operator and compare.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::VerificationFailed`] when the recomputation fails
    /// or differs.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        if self.schema != HYPERGEOMETRIC_ROW_RECURRENCE_SCHEMA {
            return Err(RecurrenceError::VerificationFailed(format!(
                "unsupported schema `{}`",
                self.schema
            )));
        }
        let replay = certify_hypergeometric_row_operator(&self.triangle, &self.operator)
            .map_err(|error| RecurrenceError::VerificationFailed(error.to_string()))?;
        if &replay != self {
            return Err(RecurrenceError::VerificationFailed(
                "the hypergeometric row-recurrence certificate differs from its recomputation"
                    .into(),
            ));
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Parsed terms
// ---------------------------------------------------------------------------

/// `ceil((a n + b) / d)` (lower) or `floor((a n + b) / d)` (upper).
#[derive(Clone, Debug)]
struct Bound {
    a: BigInt,
    b: BigInt,
    d: BigInt,
    round_up: bool,
}

impl Bound {
    fn parse(text: &str, round_up: bool) -> Result<Self, RecurrenceError> {
        let form = parse_linear(text)?;
        if !form.column.is_zero() {
            return Err(RecurrenceError::InvalidInput(format!(
                "support bound `{text}` must not depend on k"
            )));
        }
        let (d, integers) =
            vec_poly::clear_denominators_with_scale_iter([&form.index, &form.constant].into_iter());
        let [a, b]: [BigInt; 2] = integers.try_into().expect("two coefficients");
        if a.is_negative() {
            return Err(RecurrenceError::Unsupported(format!(
                "support bound `{text}` decreases in n"
            )));
        }
        Ok(Self { a, b, d, round_up })
    }

    fn at(&self, n: i64) -> i64 {
        let value = &self.a * BigInt::from(n) + &self.b;
        let rounded = if self.round_up {
            value.div_ceil(&self.d)
        } else {
            value.div_floor(&self.d)
        };
        rounded.to_i64().unwrap_or(i64::MAX)
    }

    /// The bound on the residue class `n = r mod modulus` as an exact
    /// linear form in `n` (`d` divides `modulus`).
    fn affine(&self, residue: usize) -> Linear {
        let r = BigInt::from(residue);
        let value = &self.a * &r + &self.b;
        let rounded = if self.round_up {
            value.div_ceil(&self.d)
        } else {
            value.div_floor(&self.d)
        };
        let offset = &rounded * &self.d - &value;
        Linear {
            index: Q::new(self.a.clone(), self.d.clone()),
            column: Q::zero(),
            constant: Q::new(&self.b + offset, self.d.clone()),
        }
    }
}

#[derive(Clone, Debug)]
struct Term {
    first_index: usize,
    column: Vec<Factored>,
    index: Vec<Factored>,
    lower: Bound,
    upper: Bound,
    first_value: Q,
    /// The polynomial factor `M(n, k)` (one when absent).
    multiplier: Biv,
    /// The period of ratios and supports.
    modulus: usize,
}

impl Term {
    fn parse(definition: &HypergeometricTermDefinition) -> Result<Self, RecurrenceError> {
        let period = definition.column_ratios.len();
        if period == 0 || definition.index_ratios.len() != period {
            return Err(RecurrenceError::InvalidInput(
                "column_ratios and index_ratios need the same positive length".into(),
            ));
        }
        let column = definition
            .column_ratios
            .iter()
            .map(|text| parse_factored(text))
            .collect::<Result<Vec<_>, _>>()?;
        let index = definition
            .index_ratios
            .iter()
            .map(|text| parse_factored(text))
            .collect::<Result<Vec<_>, _>>()?;
        let lower = Bound::parse(&definition.lower_column, true)?;
        let upper = Bound::parse(&definition.upper_column, false)?;
        let first_value = parse_rational_coeff(&definition.first_value)?;
        if first_value.is_zero() {
            return Err(RecurrenceError::InvalidInput(
                "the first value must be nonzero".into(),
            ));
        }
        let lcm = |x: usize, y: &BigInt| -> Result<usize, RecurrenceError> {
            let y = y.to_usize().ok_or_else(|| {
                RecurrenceError::LimitExceeded("support denominator too large".into())
            })?;
            Ok(x.lcm(&y))
        };
        let modulus = lcm(lcm(period, &lower.d)?, &upper.d)?;
        if modulus > 64 {
            return Err(RecurrenceError::LimitExceeded(format!(
                "period {modulus} of the ratios and supports exceeds 64"
            )));
        }
        let multiplier = match &definition.polynomial_factor {
            Some(text) => factored::parse_polynomial(text)?,
            None => Biv::one(),
        };
        Ok(Self {
            first_index: definition.first_index,
            column,
            index,
            lower,
            upper,
            first_value,
            multiplier,
            modulus,
        })
    }

    fn period(&self) -> usize {
        self.column.len()
    }

    fn column_ratio(&self, n_residue: usize) -> &Factored {
        &self.column[n_residue % self.period()]
    }

    fn index_ratio(&self, n_residue: usize) -> &Factored {
        &self.index[n_residue % self.period()]
    }

    /// The first `n >= first_index` with `n = residue mod modulus`.
    fn first_in_class(&self, residue: usize) -> usize {
        let n0 = self.first_index;
        n0 + (residue + self.modulus - n0 % self.modulus) % self.modulus
    }

    /// Step 1 of the module documentation.
    fn check_definition(&self) -> Result<(), RecurrenceError> {
        // Consistency of the two ratios.
        let k = Linear::column_variable();
        let k_plus_one = k.add(&Linear::constant(Q::one()));
        for r in 0..self.period() {
            let left = self
                .column_ratio(r + 1)
                .compose(1, &k)?
                .multiply(self.index_ratio(r));
            let right = self
                .index_ratio(r)
                .compose(0, &k_plus_one)?
                .multiply(self.column_ratio(r));
            if left != right {
                return Err(RecurrenceError::Unsupported(format!(
                    "the ratios are not consistent (residue {r}): rho_k(n+1,k) rho_n(n,k) = {} \
                     but rho_n(n,k+1) rho_k(n,k) = {}",
                    left.text(),
                    right.text()
                )));
            }
        }
        for residue in 0..self.modulus {
            let lower = self.lower.affine(residue);
            let upper = self.upper.affine(residue);
            let next_residue = (residue + 1) % self.modulus;
            let lower_next = self
                .lower
                .affine(next_residue)
                .compose(1, &Linear::constant(Q::zero()));
            let start = self.first_in_class(residue);
            let n = Q::from_integer(BigInt::from(start));
            let zero = Q::zero();
            // Nonempty, monotone supports with lower(n) >= 0 and lower(n+1) <= upper(n).
            let width = upper.add(&lower.scale(&-Q::one()));
            let step = upper.add(&lower_next.scale(&-Q::one()));
            if lower.evaluate(&n, &zero).is_negative()
                || width.evaluate(&n, &zero).is_negative()
                || width.index.is_negative()
                || step.evaluate(&n, &zero).is_negative()
                || step.index.is_negative()
            {
                return Err(RecurrenceError::Unsupported(format!(
                    "the support must satisfy 0 <= lower(n) <= upper(n) and lower(n+1) <= \
                     upper(n) for every n >= {} (residue {residue} mod {})",
                    self.first_index, self.modulus
                )));
            }
            let column_region = (lower.clone(), upper.add(&Linear::constant(-Q::one())));
            let index_region = (lower_next, upper.clone());
            for (ratio, region, name) in [
                (self.column_ratio(residue), column_region, "T(n,k+1)/T(n,k)"),
                (self.index_ratio(residue), index_region, "T(n+1,k)/T(n,k)"),
            ] {
                for phi in ratio.factors.keys() {
                    if !nonzero_on_region(phi, residue, self.modulus, start, &region) {
                        return Err(RecurrenceError::Unsupported(format!(
                            "the factor {} of {name} may vanish on the support (residue \
                             {residue} mod {}); the term must be nonzero on its support",
                            phi.text(),
                            self.modulus
                        )));
                    }
                }
            }
        }
        Ok(())
    }

    /// The rows of this term as `(lower(n), values)` for `count` rows.
    fn rows(&self, count: usize) -> Result<Vec<(i64, Vec<Q>)>, RecurrenceError> {
        let mut out: Vec<(i64, Vec<Q>)> = Vec::with_capacity(count);
        let n0 = i64::try_from(self.first_index)
            .map_err(|_| RecurrenceError::LimitExceeded("first index too large".into()))?;
        let evaluate = |ratio: &Factored, n: i64, k: i64| -> Result<Q, RecurrenceError> {
            ratio
                .evaluate(&Q::from_integer(n.into()), &Q::from_integer(k.into()))
                .filter(|v| !v.is_zero())
                .ok_or_else(|| {
                    RecurrenceError::Unsupported(format!(
                        "a ratio vanishes or is undefined at (n, k) = ({n}, {k})"
                    ))
                })
        };
        for step in 0..count {
            let n = n0 + step as i64;
            let residue = (n as usize) % self.period();
            let lower = self.lower.at(n);
            let upper = self.upper.at(n);
            let first = if step == 0 {
                self.first_value.clone()
            } else {
                let (previous_lower, previous) = &out[step - 1];
                let position = usize::try_from(lower - previous_lower).map_err(|_| {
                    RecurrenceError::Unsupported("the lower support bound decreases".into())
                })?;
                let previous_residue = ((n - 1) as usize) % self.period();
                let value = previous.get(position).ok_or_else(|| {
                    RecurrenceError::Unsupported(format!("row {n} starts after row {} ends", n - 1))
                })?;
                value * evaluate(self.index_ratio(previous_residue), n - 1, lower)?
            };
            let mut values = vec![first];
            for k in lower..upper {
                let next =
                    values.last().expect("nonempty") * evaluate(self.column_ratio(residue), n, k)?;
                values.push(next);
            }
            out.push((lower, values));
        }
        Ok(out)
    }

    /// The rows of `T = h M` (module documentation).
    fn rows_with_factor(&self, count: usize) -> Result<Vec<(i64, Vec<Q>)>, RecurrenceError> {
        let mut rows = self.rows(count)?;
        if self.multiplier == Biv::one() {
            return Ok(rows);
        }
        for (step, (lower, values)) in rows.iter_mut().enumerate() {
            let n = Q::from_integer(BigInt::from(self.first_index + step));
            let at_n = self.multiplier.evaluate_index(&n);
            for (offset, value) in values.iter_mut().enumerate() {
                let k = Q::from_integer(BigInt::from(*lower + offset as i64));
                *value *= vec_poly::evaluate(&at_n, &k);
            }
        }
        Ok(rows)
    }
}

/// Whether the linear form `phi` is nonzero at every lattice point
/// `(n, k)` with `n >= start`, `n = residue mod modulus` and
/// `region.0(n) <= k <= region.1(n)` (module documentation, step 1).
fn nonzero_on_region(
    phi: &Linear,
    residue: usize,
    modulus: usize,
    start: usize,
    region: &(Linear, Linear),
) -> bool {
    let (low, high) = region;
    let zero = Q::zero();
    let slope = &high.index - &low.index;
    let mut n = start;
    let gap = |n: usize| {
        let n = Q::from_integer(n.into());
        high.evaluate(&n, &zero) - low.evaluate(&n, &zero)
    };
    if gap(n).is_negative() {
        if !slope.is_positive() {
            return true; // empty region
        }
        let needed = (-gap(n) / (&slope * Q::from_integer(modulus.into()))).ceil();
        let Some(steps) = needed.to_integer().to_usize() else {
            return false;
        };
        n += steps * modulus;
        if gap(n).is_negative() {
            return false;
        }
    }
    let nq = Q::from_integer(n.into());
    let corners = [
        phi.evaluate(&nq, &low.evaluate(&nq, &zero)),
        phi.evaluate(&nq, &high.evaluate(&nq, &zero)),
    ];
    let directions = [
        &phi.index + &phi.column * &low.index,
        &phi.index + &phi.column * &high.index,
    ];
    if corners.iter().all(Signed::is_positive) && directions.iter().all(|d| !d.is_negative()) {
        return true;
    }
    if corners.iter().all(Signed::is_negative) && directions.iter().all(|d| !d.is_positive()) {
        return true;
    }
    // No integer zero at all: A modulus m + B k + (A r + C) with
    // gcd(A modulus, B) not dividing A r + C.
    let (_, integers) = vec_poly::clear_denominators_with_scale_iter(
        [&phi.index, &phi.column, &phi.constant].into_iter(),
    );
    let [index, b, constant]: [BigInt; 3] = integers.try_into().expect("three coefficients");
    let a = &index * BigInt::from(modulus);
    let c = index * BigInt::from(residue) + constant;
    let g = a.gcd(&b);
    !g.is_zero() && !c.is_multiple_of(&g)
}

/// The rows of the sum of the terms, trimmed.
fn total_rows(terms: &[Term], count: usize) -> Result<Vec<Vec<Q>>, RecurrenceError> {
    let mut rows = vec![Vec::<Q>::new(); count];
    for term in terms {
        for (row, (lower, values)) in rows.iter_mut().zip(term.rows_with_factor(count)?) {
            let lower = usize::try_from(lower)
                .map_err(|_| RecurrenceError::Unsupported("negative support".into()))?;
            if row.len() < lower + values.len() {
                row.resize(lower + values.len(), Q::zero());
            }
            for (offset, value) in values.into_iter().enumerate() {
                row[lower + offset] += value;
            }
        }
    }
    Ok(rows.into_iter().map(|row| vec_poly::trim(&row)).collect())
}

// ---------------------------------------------------------------------------
// Identities
// ---------------------------------------------------------------------------

/// One identity block: the column index `k` of the coefficient
/// (`k` itself in the interior, `lower(n) + m` or `upper(n) + m` in a strip)
/// and, for every pair `(i, j)` whose point is in the support, the
/// numerator `N_(i,j)` of `T(n + i, k - j) / base` over the block's common
/// denominator.  The identity is `sum c_(i,j,l)(n) (k - j)^l N_(i,j) = 0`.
struct Block {
    column_value: Linear,
    numerators: BTreeMap<(usize, usize), Biv>,
}

/// The ratio `T(n + i, k') / base` along a path, as a factored function.
struct PathBuilder<'a> {
    term: &'a Term,
    residue: usize,
}

impl PathBuilder<'_> {
    fn residue_after(&self, s: usize) -> usize {
        (self.residue + s) % self.term.modulus
    }

    /// `prod_(s < i) rho_n(n + s, column)`.
    fn up(
        &self,
        from_shift: usize,
        steps: usize,
        column: &Linear,
    ) -> Result<Factored, RecurrenceError> {
        let mut out = Factored::one();
        for s in from_shift..from_shift + steps {
            let ratio = self.term.index_ratio(self.residue_after(s));
            out = out.multiply(&ratio.compose(s as i64, column)?);
        }
        Ok(out)
    }

    /// Along row `n + shift` from `column` by `steps` (negative: backwards).
    fn along(
        &self,
        shift: usize,
        column: &Linear,
        steps: i64,
    ) -> Result<Factored, RecurrenceError> {
        let ratio = self.term.column_ratio(self.residue_after(shift));
        let mut out = Factored::one();
        if steps >= 0 {
            for u in 0..steps {
                let at = column.add(&Linear::constant(Q::from_integer(u.into())));
                out = out.multiply(&ratio.compose(shift as i64, &at)?);
            }
        } else {
            for u in 1..=-steps {
                let at = column.add(&Linear::constant(Q::from_integer((-u).into())));
                out = out.multiply(&ratio.compose(shift as i64, &at)?.inverse()?);
            }
        }
        Ok(out)
    }
}

/// The blocks of one term and residue class for the pairs `pairs`, with
/// the strip counts and the separation index (module documentation, steps
/// 2 and 3).
fn blocks_for(
    term: &Term,
    residue: usize,
    pairs: &BTreeSet<(usize, usize)>,
) -> Result<(Vec<Block>, HypergeometricBoundaryStrips), RecurrenceError> {
    let order = pairs.iter().map(|(i, _)| *i).max().unwrap_or(0);
    let width = pairs.iter().map(|(_, j)| *j).max().unwrap_or(0);
    let path = PathBuilder { term, residue };
    let modulus = term.modulus;
    let lower = |s: usize| {
        term.lower
            .affine((residue + s) % modulus)
            .compose(s as i64, &Linear::constant(Q::zero()))
    };
    let upper = |s: usize| {
        term.upper
            .affine((residue + s) % modulus)
            .compose(s as i64, &Linear::constant(Q::zero()))
    };
    let difference = |a: &Linear, b: &Linear| -> Result<i64, RecurrenceError> {
        let d = a.add(&b.scale(&-Q::one()));
        if !d.index.is_zero() || !d.constant.is_integer() {
            return Err(RecurrenceError::VerificationFailed(
                "support offsets are not integer constants".into(),
            ));
        }
        d.constant
            .to_integer()
            .to_i64()
            .ok_or_else(|| RecurrenceError::LimitExceeded("support offset".into()))
    };
    let lower_steps = (0..=order)
        .map(|i| difference(&lower(i), &lower(0)))
        .collect::<Result<Vec<_>, _>>()?;
    let upper_steps = (0..=order)
        .map(|i| difference(&upper(i), &upper(0)))
        .collect::<Result<Vec<_>, _>>()?;
    let width_i64 = width as i64;
    let lower_strips = lower_steps.iter().copied().max().unwrap_or(0) + width_i64;
    let upper_strips = upper_steps.iter().copied().max().unwrap_or(0) + width_i64;

    let mut blocks = Vec::new();
    // Interior: base (n, k).
    let k = Linear::column_variable();
    let mut interior = BTreeMap::new();
    for &(i, j) in pairs {
        let ratio = path
            .up(0, i, &k)?
            .multiply(&path.along(i, &k, -(j as i64))?);
        let at = k.add(&Linear::constant(-Q::from_integer(j.into())));
        let factor = factored::compose_polynomial(&term.multiplier, i as i64, &at);
        interior.insert((i, j), (ratio, factor));
    }
    blocks.push(make_block(k, interior));
    // Lower strips: base (n, lower(n)).
    let base_low = lower(0);
    for m in 0..lower_strips {
        let mut ratios = BTreeMap::new();
        for &(i, j) in pairs {
            let offset = m - j as i64;
            if offset < lower_steps[i] {
                continue;
            }
            let column = base_low.add(&Linear::constant(Q::from_integer(offset.into())));
            let ratio = path
                .along(0, &base_low, offset)?
                .multiply(&path.up(0, i, &column)?);
            let factor = factored::compose_polynomial(&term.multiplier, i as i64, &column);
            ratios.insert((i, j), (ratio, factor));
        }
        let value = base_low.add(&Linear::constant(Q::from_integer(m.into())));
        blocks.push(make_block(value, ratios));
    }
    // Upper strips: base (n, upper(n)).
    let base_high = upper(0);
    for m in 1..=upper_strips {
        let mut ratios = BTreeMap::new();
        for &(i, j) in pairs {
            let offset = m - j as i64;
            if offset > upper_steps[i] {
                continue;
            }
            let ratio = path
                .up(0, i, &base_high)?
                .multiply(&path.along(i, &base_high, offset)?);
            let column = base_high.add(&Linear::constant(Q::from_integer(offset.into())));
            let factor = factored::compose_polynomial(&term.multiplier, i as i64, &column);
            ratios.insert((i, j), (ratio, factor));
        }
        let value = base_high.add(&Linear::constant(Q::from_integer(m.into())));
        blocks.push(make_block(value, ratios));
    }
    // Separation index: lower(n) + lower_strips - 1 <= upper(n),
    // upper(n) + 1 - width >= lower(n + order) and upper(n) >= lower(n + order).
    let conditions = [
        upper(0)
            .add(&lower(0).scale(&-Q::one()))
            .add(&Linear::constant(Q::from_integer(
                (1 - lower_strips).into(),
            ))),
        upper(0)
            .add(&lower(order).scale(&-Q::one()))
            .add(&Linear::constant(Q::from_integer((1 - width_i64).into()))),
        upper(0).add(&lower(order).scale(&-Q::one())),
    ];
    let start = term.first_in_class(residue);
    let mut separation = start;
    for condition in &conditions {
        separation = separation.max(first_nonnegative(condition, start, modulus)?);
    }
    let strips = HypergeometricBoundaryStrips {
        term: 0,
        residue,
        modulus,
        lower_strips: lower_strips as usize,
        upper_strips: upper_strips as usize,
        separation_index: separation,
    };
    Ok((blocks, strips))
}

/// The first `n >= start` in the residue class of `start` modulo `modulus`
/// with `condition(n) >= 0`, for a condition nondecreasing in `n`.
fn first_nonnegative(
    condition: &Linear,
    start: usize,
    modulus: usize,
) -> Result<usize, RecurrenceError> {
    let zero = Q::zero();
    let value = condition.evaluate(&Q::from_integer(start.into()), &zero);
    if !value.is_negative() {
        return Ok(start);
    }
    if !condition.index.is_positive() {
        return Err(RecurrenceError::Unsupported(
            "the rows are too narrow for the boundary strips to separate".into(),
        ));
    }
    let steps = (-value / (&condition.index * Q::from_integer(modulus.into()))).ceil();
    let steps = steps
        .to_integer()
        .to_usize()
        .ok_or_else(|| RecurrenceError::LimitExceeded("separation index".into()))?;
    Ok(start + steps * modulus)
}

fn make_block(column_value: Linear, ratios: BTreeMap<(usize, usize), (Factored, Biv)>) -> Block {
    let mut denominator: BTreeMap<Linear, i64> = BTreeMap::new();
    for (ratio, _) in ratios.values() {
        for (phi, e) in ratio.denominator_exponents() {
            let entry = denominator.entry(phi).or_insert(0);
            *entry = (*entry).max(e);
        }
    }
    let numerators = ratios
        .into_iter()
        .map(|(key, (ratio, factor))| (key, ratio.numerator_over(&denominator).multiply(&factor)))
        .collect();
    Block {
        column_value,
        numerators,
    }
}

/// `(column_value - j)^l` as a polynomial in `(n, k)`.
fn column_power(column_value: &Linear, j: usize, l: usize) -> Biv {
    column_value
        .add(&Linear::constant(-Q::from_integer(j.into())))
        .to_biv()
        .power(l)
}

/// `sum c (k - j)^l N_(i,j)` of a block for an operator.
fn block_value(block: &Block, operator: &EulerMap) -> Biv {
    let mut total = Biv::zero();
    for ((i, j, l), c) in operator {
        let Some(numerator) = block.numerators.get(&(*i, *j)) else {
            continue;
        };
        let term = Biv::from_index_polynomial(c)
            .multiply(&column_power(&block.column_value, *j, *l))
            .multiply(numerator);
        total = total.add(&term);
    }
    total
}

// ---------------------------------------------------------------------------
// Certification
// ---------------------------------------------------------------------------

/// Parse and normalize operator terms: `(i, j, l) -> c(n)`.
fn parse_operator(operator: &[EulerOperatorTerm]) -> Result<EulerMap, RecurrenceError> {
    let mut out: EulerMap = BTreeMap::new();
    for term in operator {
        let c = term
            .coefficient
            .iter()
            .map(|text| parse_rational_coeff(text))
            .collect::<Result<Vec<_>, _>>()?;
        let key = (term.shift, term.variable_power, term.euler_power);
        let sum = vec_poly::add(out.get(&key).map(Vec::as_slice).unwrap_or(&[]), &c);
        out.insert(key, sum);
    }
    out.retain(|_, c| !vec_poly::is_zero(c));
    Ok(out)
}

fn operator_terms(operator: &EulerMap) -> Vec<EulerOperatorTerm> {
    operator
        .iter()
        .map(|((i, j, l), c)| EulerOperatorTerm {
            shift: *i,
            variable_power: *j,
            euler_power: *l,
            coefficient: vec_poly::trim(c)
                .iter()
                .map(format_rational_coeff)
                .collect(),
        })
        .collect()
}

fn operator_text(operator: &EulerMap) -> String {
    let mut parts = Vec::new();
    for ((i, j, l), c) in operator.iter().rev() {
        let coefficient = crate::recurrence::implicit::algebra::poly_text(c, "n");
        let t = match j {
            0 => String::new(),
            1 => " t".into(),
            _ => format!(" t^{j}"),
        };
        let theta = match l {
            0 => String::new(),
            1 => " theta".into(),
            _ => format!(" theta^{l}"),
        };
        let row = match i {
            0 => "P(n)".to_string(),
            _ => format!("P(n+{i})"),
        };
        parts.push(format!("({coefficient}){t}{theta} {row}"));
    }
    format!("{} = 0", parts.join(" + "))
}

/// The derivative form of `E` as an Ore operator.
fn derivative_operator(operator: &EulerMap) -> DOp {
    let mut theta = ThetaOp::new();
    for ((i, j, l), c) in operator {
        let mut rows: Vec<Vec<Q>> = c.iter().map(|_| vec![Q::zero(); l + 1]).collect();
        for (a, value) in c.iter().enumerate() {
            rows[a][*l] = value.clone();
        }
        theta_add(&mut theta, (*i, *j), Biv::new(rows));
    }
    from_theta(&theta)
}

/// Prove `E P_n = 0` for every `n >= n_0` (module documentation) and build
/// the recurrence record of `E`.
///
/// The operator must have a single term at its largest shift `I >= 1`, with
/// `variable_power = euler_power = 0` (a leader free of `t` and `theta`).
///
/// # Errors
///
/// [`RecurrenceError::InvalidInput`] for a malformed triangle or operator,
/// [`RecurrenceError::Unsupported`] when the definition cannot be checked
/// (a ratio factor that may vanish on the support, rows too narrow), and
/// [`RecurrenceError::VerificationFailed`] when an identity fails (with the
/// failing block or row).
pub fn certify_hypergeometric_row_operator(
    triangle: &HypergeometricTriangle,
    operator: &[EulerOperatorTerm],
) -> Result<HypergeometricRowRecurrenceCertificate, RecurrenceError> {
    certify_with_limit(
        triangle,
        operator,
        HypergeometricRowRecurrenceBudget::default().max_checked_row,
    )
}

fn certify_with_limit(
    triangle: &HypergeometricTriangle,
    operator: &[EulerOperatorTerm],
    max_checked_row: usize,
) -> Result<HypergeometricRowRecurrenceCertificate, RecurrenceError> {
    let terms = triangle.parse()?;
    let n0 = triangle.first_index()?;
    let parsed = parse_operator(operator)?;
    let order = parsed
        .keys()
        .map(|(i, _, _)| *i)
        .max()
        .ok_or_else(|| RecurrenceError::InvalidInput("the operator is zero".into()))?;
    let leader_keys: Vec<_> = parsed.keys().filter(|(i, _, _)| *i == order).collect();
    if order == 0 || leader_keys != vec![&(order, 0, 0)] {
        return Err(RecurrenceError::InvalidInput(
            "the operator needs order >= 1 and a leader c(n) S^I free of t and theta".into(),
        ));
    }
    // Rows: before the separation index, and seeds past the leader's roots.
    let leader = parsed[&(order, 0, 0)].clone();
    let roots = crate::recurrence::implicit::algebra::roots_in_domain(&leader, n0, None);
    let seeds = order.max(roots.last().map_or(0, |r| r + order + 1 - n0));
    let proof = prove_annihilation(&terms, n0, n0, &parsed, seeds, max_checked_row)?;
    let derivative = derivative_operator(&parsed);
    let rows = proof.rows;
    let strips = proof.strips;
    let checked = proof.checked_bases;
    let (rule, shift) = ImplicitRelation::from_op(&derivative)
        .to_successor_rule()
        .ok_or_else(|| {
            RecurrenceError::VerificationFailed("the operator is not runnable".into())
        })?;
    debug_assert_eq!(shift, order);
    let record = RecurrenceJson::from_recurrence_rational(&rule, n0, &rows[..seeds], None);
    let text = operator_text(&parsed);
    let statement = format!(
        "For every n >= {n0}: {text} with P_n(t) = sum_k T(n,k) t^k (theta = t d/dt); the \
         record with the seed rows {n0}..{} generates exactly the rows P_n.",
        n0 + seeds - 1
    );
    Ok(HypergeometricRowRecurrenceCertificate {
        schema: HYPERGEOMETRIC_ROW_RECURRENCE_SCHEMA.into(),
        triangle: triangle.clone(),
        operator: operator_terms(&parsed),
        operator_text: text,
        boundary_strips: strips,
        checked_bases: checked,
        record,
        statement,
    })
}

/// JSON schema tag of [`HypergeometricRecordCertificate`].
pub const HYPERGEOMETRIC_RECORD_SCHEMA: &str = "polytool.hypergeometric-record.v1";

/// A proof that a recurrence record generates the rows of a hypergeometric
/// triangle for every index ([`certify_hypergeometric_row_record`]),
/// schema `polytool.hypergeometric-record.v1`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HypergeometricRecordCertificate {
    /// [`HYPERGEOMETRIC_RECORD_SCHEMA`].
    pub schema: String,
    /// The definition.
    pub triangle: HypergeometricTriangle,
    /// The record (search metadata and presentation are ignored).
    pub record: RecurrenceJson,
    /// The rule in forward Euler form at base `n = m - R` (`R` the largest
    /// lag of the record), multiplied on the left by
    /// `t^variable_power_multiplier`, which clears
    /// `D^d = t^(-d) theta (theta - 1) ... (theta - d + 1)`.
    pub operator: Vec<EulerOperatorTerm>,
    /// The power of `t` the rule was multiplied by.
    pub variable_power_multiplier: usize,
    /// The proved strip identities and separation indices.
    pub boundary_strips: Vec<HypergeometricBoundaryStrips>,
    /// Bases at which the rule was checked exactly on rows.
    pub checked_bases: Vec<usize>,
    /// What the certificate establishes.
    pub statement: String,
}

impl HypergeometricRecordCertificate {
    /// Replay: recompute the certificate from the recorded triangle and
    /// record and compare.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::VerificationFailed`] when the recomputation fails
    /// or differs.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        if self.schema != HYPERGEOMETRIC_RECORD_SCHEMA {
            return Err(RecurrenceError::VerificationFailed(format!(
                "unsupported schema `{}`",
                self.schema
            )));
        }
        let replay = certify_hypergeometric_row_record(&self.triangle, &self.record)
            .map_err(|error| RecurrenceError::VerificationFailed(error.to_string()))?;
        if &replay != self {
            return Err(RecurrenceError::VerificationFailed(
                "the hypergeometric record certificate differs from its recomputation".into(),
            ));
        }
        Ok(())
    }
}

/// Prove that `record` generates the rows `P_m` of `triangle` for every
/// `m >= record.first_index` (which must be at least the triangle's first
/// index).
///
/// The rule `f(m, t) P_m = sum c(m, t) D^d P_(m-r)` is written in forward
/// form at base `n = m - R` and multiplied by `t^J` (`J` its largest
/// derivative order), which gives an operator in `t` and `theta`; the
/// module's argument proves that it annihilates the triangle's rows at
/// every base `n >= record.first_index + (seed count) - R` (and the
/// identity transfers back, since `t^J` is not a zero divisor).  The seeds
/// are compared with the rows exactly, and the left-hand factor `f(m, .)`
/// is a nonzero polynomial at every row the record computes (its content
/// in `m` has no integer root there), so by induction the record's rows are
/// the triangle's rows.
///
/// # Errors
///
/// As [`certify_hypergeometric_row_operator`]; a seed that differs from the
/// triangle's row, a vanishing left-hand factor, an inhomogeneous or
/// sign-alternating rule ([`RecurrenceError::Unsupported`]).
pub fn certify_hypergeometric_row_record(
    triangle: &HypergeometricTriangle,
    record: &RecurrenceJson,
) -> Result<HypergeometricRecordCertificate, RecurrenceError> {
    let terms = triangle.parse()?;
    let n0 = triangle.first_index()?;
    let (rule, first, seeds) = record.to_recurrence_parts()?;
    if first < n0 {
        return Err(RecurrenceError::InvalidInput(format!(
            "the record starts at {first}, before the triangle's first index {n0}"
        )));
    }
    let (relation, lag) = ImplicitRelation::from_recurrence(&rule)
        .map_err(|error| RecurrenceError::Unsupported(error.to_string()))?;
    let op = relation.op();
    let computed_from = first + seeds.len();
    if computed_from < first + lag {
        return Err(RecurrenceError::InvalidInput(
            "the record has fewer seeds than its largest lag".into(),
        ));
    }
    let first_base = computed_from - lag;
    // The left-hand factor at base n is the coefficient of S^R.
    let leader =
        op.0.get(&(lag, 0))
            .cloned()
            .ok_or_else(|| RecurrenceError::InvalidInput("the rule has no leader".into()))?;
    let (content, _) = crate::recurrence::implicit::algebra::content_and_primitive(&leader);
    let roots = crate::recurrence::implicit::algebra::roots_in_domain(&content, first_base, None);
    if let Some(root) = roots.first() {
        return Err(RecurrenceError::Unsupported(format!(
            "the left-hand factor vanishes at the row {}",
            root + lag
        )));
    }
    let (multiplier, theta) = crate::recurrence::implicit::algebra::to_theta(&op);
    let mut parsed: EulerMap = BTreeMap::new();
    for ((i, j), c) in &theta {
        let width = c.variable_degree().map_or(0, |d| d + 1);
        for l in 0..width {
            let column = vec_poly::trim(&c.variable_coefficient(l));
            if !vec_poly::is_zero(&column) {
                parsed.insert((*i, *j, l), column);
            }
        }
    }
    let rows_needed = computed_from - n0;
    let proof = prove_annihilation(&terms, n0, first_base, &parsed, rows_needed, 400)?;
    for (offset, seed) in seeds.iter().enumerate() {
        let m = first + offset;
        if vec_poly::trim(seed) != proof.rows[m - n0] {
            return Err(RecurrenceError::VerificationFailed(format!(
                "the seed row {m} differs from the triangle's row"
            )));
        }
    }
    let text = operator_text(&parsed);
    let statement = format!(
        "For every m >= {first} the record's row equals P_m(t) = sum_k T(n,k) t^k: the seeds \
         {first}..{} are rows of the triangle, t^{multiplier} times the rule is {text} at \
         base n = m - {lag}, proved for every n >= {first_base}, and the left-hand factor does \
         not vanish at a computed row.",
        computed_from - 1
    );
    let mut record = record.clone();
    record.search = None;
    record.latex = None;
    Ok(HypergeometricRecordCertificate {
        schema: HYPERGEOMETRIC_RECORD_SCHEMA.into(),
        triangle: triangle.clone(),
        record,
        operator: operator_terms(&parsed),
        variable_power_multiplier: multiplier,
        boundary_strips: proof.strips,
        checked_bases: proof.checked_bases,
        statement,
    })
}

/// The proof of `E P_n = 0` for every `n >= first_base` (module
/// documentation, steps 2 to 4) and the rows it evaluated.
struct AnnihilationProof {
    strips: Vec<HypergeometricBoundaryStrips>,
    checked_bases: Vec<usize>,
    /// The rows `P_(n_0), P_(n_0 + 1), ...` (at least `minimum_rows`).
    rows: Vec<Vec<Q>>,
}

fn prove_annihilation(
    terms: &[Term],
    n0: usize,
    first_base: usize,
    parsed: &EulerMap,
    minimum_rows: usize,
    max_checked_row: usize,
) -> Result<AnnihilationProof, RecurrenceError> {
    let order = parsed.keys().map(|(i, _, _)| *i).max().unwrap_or(0);
    let pairs: BTreeSet<(usize, usize)> = parsed.keys().map(|(i, j, _)| (*i, *j)).collect();
    let mut strips = Vec::new();
    let mut separation = first_base;
    for (position, term) in terms.iter().enumerate() {
        for residue in 0..term.modulus {
            let (blocks, mut info) = blocks_for(term, residue, &pairs)?;
            for (b, block) in blocks.iter().enumerate() {
                if !block_value(block, parsed).is_zero() {
                    let what = if b == 0 {
                        "the interior identity".to_string()
                    } else if b <= info.lower_strips {
                        format!("the boundary strip k = lower(n) + {}", b - 1)
                    } else {
                        format!(
                            "the boundary strip k = upper(n) + {}",
                            b - info.lower_strips
                        )
                    };
                    return Err(RecurrenceError::VerificationFailed(format!(
                        "term {position}, residue {residue} mod {}: {what} fails",
                        term.modulus
                    )));
                }
            }
            info.term = position;
            separation = separation.max(info.separation_index);
            strips.push(info);
        }
    }
    let needed = (separation - n0 + order).max(minimum_rows);
    if n0 + needed > max_checked_row {
        return Err(RecurrenceError::BudgetExhausted(format!(
            "the certificate needs rows up to {} (limit {max_checked_row})",
            n0 + needed
        )));
    }
    let rows = total_rows(terms, needed)?;
    let derivative = derivative_operator(parsed);
    let mut checked = Vec::new();
    for n in first_base..separation {
        let value = derivative.apply(&rows, n0, n).ok_or_else(|| {
            RecurrenceError::VerificationFailed(format!("row check at base {n}: rows missing"))
        })?;
        if !vec_poly::is_zero(&value) {
            return Err(RecurrenceError::VerificationFailed(format!(
                "E P_n != 0 at the base n = {n} (exact rows)"
            )));
        }
        checked.push(n);
    }
    Ok(AnnihilationProof {
        strips,
        checked_bases: checked,
        rows,
    })
}

// ---------------------------------------------------------------------------
// Search
// ---------------------------------------------------------------------------

/// Find and certify an annihilating operator of the rows of `triangle`
/// (module documentation, "The search").
///
/// # Errors
///
/// As [`certify_hypergeometric_row_operator`] for the definition;
/// [`RecurrenceError::BudgetExhausted`] when no operator exists within the
/// budget.
pub fn derive_hypergeometric_row_recurrence(
    triangle: &HypergeometricTriangle,
    budget: &HypergeometricRowRecurrenceBudget,
) -> Result<HypergeometricRowRecurrenceCertificate, RecurrenceError> {
    let terms = triangle.parse()?;
    for order in 1..=budget.max_order {
        for euler in 0..=budget.max_euler_power {
            for variable in 0..=budget.max_variable_power {
                let unknowns =
                    |degree: usize| (order * (variable + 1) * (euler + 1) + 1) * (degree + 1);
                let top = (0..=budget.max_index_degree)
                    .rev()
                    .find(|d| unknowns(*d) <= budget.max_unknowns);
                let Some(top) = top else { continue };
                if solve(&terms, order, variable, euler, top)?.is_none() {
                    continue;
                }
                for degree in 0..=top {
                    if let Some(operator) = solve(&terms, order, variable, euler, degree)? {
                        let terms = operator_terms(&operator);
                        return certify_with_limit(triangle, &terms, budget.max_checked_row);
                    }
                }
            }
        }
    }
    Err(RecurrenceError::BudgetExhausted(
        "no annihilating operator with a free leader within the budget".into(),
    ))
}

/// The kernel vector with a nonzero leader of the ansatz, made primitive.
fn solve(
    terms: &[Term],
    order: usize,
    variable: usize,
    euler: usize,
    degree: usize,
) -> Result<Option<EulerMap>, RecurrenceError> {
    // Unknown columns: (i, j, l, d), the leader block last.
    let mut columns = Vec::new();
    for i in 0..order {
        for j in 0..=variable {
            for l in 0..=euler {
                for d in 0..=degree {
                    columns.push((i, j, l, d));
                }
            }
        }
    }
    let leader_start = columns.len();
    for d in 0..=degree {
        columns.push((order, 0, 0, d));
    }
    let mut pairs: BTreeSet<(usize, usize)> = BTreeSet::new();
    for i in 0..order {
        for j in 0..=variable {
            pairs.insert((i, j));
        }
    }
    pairs.insert((order, 0));
    // Equations: one row per (block, monomial).
    let mut equations: BTreeMap<(usize, usize, usize), Vec<(usize, Q)>> = BTreeMap::new();
    let mut block_id = 0;
    for term in terms {
        for residue in 0..term.modulus {
            let (blocks, _) = blocks_for(term, residue, &pairs)?;
            for block in &blocks {
                for (column, (i, j, l, d)) in columns.iter().enumerate() {
                    let Some(numerator) = block.numerators.get(&(*i, *j)) else {
                        continue;
                    };
                    let value = Biv::monomial(*d, 0, Q::one())
                        .multiply(&column_power(&block.column_value, *j, *l))
                        .multiply(numerator);
                    for (a, row) in value.rows().iter().enumerate() {
                        for (b, c) in row.iter().enumerate() {
                            if !c.is_zero() {
                                equations
                                    .entry((block_id, a, b))
                                    .or_default()
                                    .push((column, c.clone()));
                            }
                        }
                    }
                }
                block_id += 1;
            }
        }
    }
    let width = columns.len();
    let matrix: Vec<Vec<BigInt>> = equations
        .values()
        .map(|entries| {
            let (_, integers) =
                vec_poly::clear_denominators_with_scale_iter(entries.iter().map(|(_, c)| c));
            let mut row = vec![BigInt::zero(); width];
            for ((column, _), integer) in entries.iter().zip(integers) {
                row[*column] += integer;
            }
            row
        })
        .collect();
    let group: Vec<usize> = (0..width).map(|c| usize::from(c >= leader_start)).collect();
    let vectors = ore_algebra::exact_kernel_vectors_by_group(&matrix, &group, 2);
    let Some(Some(vector)) = vectors.into_iter().nth(1) else {
        return Ok(None);
    };
    let mut operator: EulerMap = BTreeMap::new();
    for ((i, j, l, d), value) in columns.iter().zip(vector) {
        if value.is_zero() {
            continue;
        }
        let entry = operator.entry((*i, *j, *l)).or_default();
        if entry.len() <= *d {
            entry.resize(d + 1, Q::zero());
        }
        entry[*d] = value;
    }
    operator.retain(|_, c| !vec_poly::is_zero(c));
    if !operator.contains_key(&(order, 0, 0)) {
        return Ok(None);
    }
    Ok(Some(primitive(operator)))
}

/// Divide by the gcd of the coefficients in `Q[n]` and scale to coprime
/// integers with a positive leading coefficient of the leader.
fn primitive(mut operator: EulerMap) -> EulerMap {
    let mut gcd = crate::polynomial::Polynomial::<Q>::zero();
    for c in operator.values() {
        gcd = gcd.gcd(&crate::polynomial::Polynomial::new(c.clone()));
    }
    if vec_poly::degree(gcd.coeffs()).unwrap_or(0) > 0 {
        for c in operator.values_mut() {
            let (quotient, _) = crate::polynomial::Polynomial::new(c.clone()).div_rem(&gcd);
            *c = vec_poly::trim(quotient.coeffs());
        }
    }
    let all: Vec<Q> = operator.values().flatten().cloned().collect();
    let (denominator, integers) = vec_poly::clear_denominators_with_scale(&all);
    let numerator = integers.iter().fold(BigInt::zero(), |acc, c| acc.gcd(c));
    let leader_sign = operator
        .iter()
        .next_back()
        .and_then(|(_, c)| c.iter().rev().find(|v| !v.is_zero()))
        .is_some_and(Signed::is_negative);
    let mut scale = Q::new(denominator, numerator);
    if leader_sign {
        scale = -scale;
    }
    for c in operator.values_mut() {
        *c = vec_poly::trim(&c.iter().map(|v| v * &scale).collect::<Vec<_>>());
    }
    operator
}
