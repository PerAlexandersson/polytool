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

mod terms;
use terms::*;
mod operator;
use operator::*;
mod certificate;
pub use certificate::{
    certify_hypergeometric_row_operator, certify_hypergeometric_row_record,
    derive_hypergeometric_row_recurrence, HypergeometricRecordCertificate,
    HYPERGEOMETRIC_RECORD_SCHEMA,
};
