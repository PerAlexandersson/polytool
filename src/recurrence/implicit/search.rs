//! Direct search for implicit differential relations and coefficient
//! conditions, with completeness certificates where a finite verification
//! theorem applies.
//!
//! # Ansatz
//!
//! [`search_implicit_relations`] fits relations
//! `sum_(i,j) c_(i,j)(n,x) D^j Q_(n+i) = 0` (base index `n`) to exact rows
//! `Q_n`, over a prescribed support [`ImplicitSearchSupport`]: shifts
//! `0..=R`, derivative orders `0..=J` on every row (the newest included),
//! `deg_n c <= U` and `deg_x c <= A`. The role selects the leader: a
//! `successor` support allows no derivative on the row `n + R`, an
//! `implicit_successor` support allows them, and a `row_equation` support has
//! `R = 0`. The rows are those of a defining rule with seeds, of a bundled
//! OEIS family or supplied rows, optionally after a registry map `Q = T(P)`.
//!
//! The columns `n^u x^v D^s Q_(n+j)` are those of the complete-space module
//! (`RecurrenceMonomial`), evaluated by the same residual-matrix code; one
//! exact kernel with a rank certificate (`RecurrenceSpaceCertificate`) is
//! computed per support. Basis vectors with independent leader projections
//! (every basis vector for a row equation) become representatives, each in
//! canonical form, with its role, evidence, a generating description and an
//! assessment ([`assess_generating_description`]).
//!
//! # Completeness
//!
//! When the source is a constant-index monic pure-lag rule
//! `chi = S^k - sum a_i(x) S^i` valid from its origin and the map is one of
//! the maps below, a finite theorem bounds the residuals that must vanish
//! (`u`, `s` are the index degree and derivative order of a column
//! `x^v n^u D^s Q_(n+j)`; `U`, `A`, `J` the maxima of `u`, `v`, `s` over the
//! support):
//!
//! - identity: every column `x^v n^u D^s P_(n+j)` is annihilated by
//!   `chi^(u+s+1)` (module `fit::space`), so `h = max(u+s) + 1`;
//! - `Q = D^m P`: the column `x^v n^u D^s Q_(n+j)` is `x^v n^u D^(s+m) P_(n+j)`,
//!   so `h = max(u+s) + m + 1`;
//! - factorial rescaling: `h = U + m(A+J) + 1` (module `fit::rescaling`);
//! - Möbius-class maps (`multiply_power`, `divide_power`, `affine`, and
//!   `reversal`, `mobius`, `f_polynomial` with an unbounded degree schedule):
//!   `h = max(u+s) + 1` (transported-source theorem below);
//! - diagonal maps with weight `lambda(n,k) = c (k!)^(-m) prod f_j(n,k)^(e_j)`,
//!   every `e_j > 0` (`binomial {m}`, and `product` weights with
//!   `factorial_exponent <= 0` and positive exponents): with `w` the total
//!   `(n,k)` degree of `prod f_j^(e_j)`, `h = max(u+s) + w + 1` for `m = 0`
//!   and `h = U + m(A+J) + w + 1` for `m >= 1` (weight theorem below).
//!
//! The residual is annihilated by the monic operator `chi^h` (or `chi~^h`) of
//! order `k h`, so its vanishing at the `k h` bases `n_0..n_0+kh-1` (`n_0` the
//! domain start, at least the origin) implies vanishing for every `n >= n_0`.
//! The kernel of those residual rows is then the complete space of all-index
//! relations in the support (evidence `certified_complete_space`). This holds
//! for derivatives on the leader row and for row equations alike.
//!
//! **Lemma (monic recurrences over `Q(x)`).** Let `chi = S^k - sum a_i S^i`
//! with `a_i` in the differential field `(Q(x), D)`, independent of `n`, and
//! `chi Y = 0` for a sequence of rational functions. Then
//! `chi^(u+s+1)` annihilates `x^v n^u D^s Y_(n+j)`, and a sequence annihilated
//! by `chi^h` that vanishes at `k h` consecutive bases vanishes from there on.
//! The proof in module `fit::space` uses only that `chi` is monic, index free and
//! commutes with `S`, with multiplication by `x` and with its own coefficient
//! derivatives (`chi D = D chi - chi'`, `[n, chi] = -S chi_S`); all of this
//! holds over `Q(x)`.
//!
//! **Transported-source theorem.** Let `Q_n = T P_n` for a Möbius-class map:
//! `Q_n = ell^(d_n) P_n(phi)` with `phi = (a x + b)/ell`, `ell = c x + e`,
//! `a e - b c != 0`, and `d_n = w + h n` on a domain containing every
//! `n >= origin` (reversal, `mobius`, `f_polynomial`), or `Q_n = x^r P_n`,
//! `P_n / x^r`, `P_n(a x + b)`. Then
//! `Q_(n+k) = ell^(d_(n+k)) sum_i a_i(phi) P_(n+i)(phi)
//!  = sum_i a~_i Q_(n+i)` with `a~_i = a_i(phi) ell^(h (k - i))` (resp.
//! `a~_i = a_i` and `a_i(a x + b)`), rational functions independent of `n`.
//! So `Q` itself satisfies a monic constant-index rule `chi~` of order `k`
//! over `Q(x)`, and the lemma gives `h = max(u+s) + 1` and the residual count
//! `k (max(u+s) + 1)`. The identity is one of rational functions; the rows
//! are polynomials where the map's degree (or divisibility) condition holds,
//! which `map_row` checks on every generated row. A bounded schedule
//! (`last_index`) is not covered.
//!
//! **Weight theorem.** Let `Q_n = T_n P_n` with
//! `T_n(x^k) = lambda(n,k) x^k`, `lambda = c (k!)^(-m) Lambda(n,k)` and
//! `Lambda = prod f_j^(e_j)` a polynomial of total degree `w`. Write
//! `theta = x D` and `F(x^k) = x^k/(k!)^m`, so `T_n = c F Lambda(n,theta)`.
//! Multiply a residual by `x^J` and apply `F^-1` (both injective):
//! `x^J x^a D^j = x^s theta^(j falling)` with `s = a + J - j >= 0`, and
//! `F^-1 x^s G(theta) F = x^s G(theta) prod_(h=1..s) (theta + h)^m`. A column
//! `n^u x^a D^j Q_(n+i)` becomes
//! `c n^u x^s theta^(j falling) prod_(h=1..s) (theta+h)^m Lambda(n+i,theta) P_(n+i)`,
//! a combination of `n^(u+u') x^v D^t P_(n+i)` with
//! `u' + t <= j + m s + w` (`theta^e` has derivative order `e`). Hence `chi^h`
//! annihilates the transformed residual with `h = max(u+j) + w + 1` for
//! `m = 0` (then `x^J` and `F` are not needed) and `h = U + m(A+J) + w + 1`
//! for `m >= 1` (as `j + m s <= m (A+J)`). For `Lambda = 1` this is the
//! factorial theorem. Weights with `1/f_j` (`inverse_binomial`, negative
//! exponents) or `(k!)^p`, `p > 0`, are not covered: their pulled-back
//! columns are not polynomial in `theta`.
//!
//! Without a finite theorem a representative is a `finite_fit`, unless it is
//! proved by right division:
//! for a homogeneous defining rule `L` and a map that is a differential
//! operator `Q = T P` (identity, `derivative`, a polynomial diagonal weight),
//! an identity `M o T = sum_nu H_nu o S^nu o L`
//! ([`find_reconstruction_certificate`]) proves `M` at every base where `L`
//! holds (evidence `certified_by_division`). A support whose kernel has no
//! vector with a nonzero leader is `absent_in_support`; this needs only the
//! rows (`rows_only_obstruction`).
//!
//! # Coefficient conditions
//!
//! When a representative's uniqueness argument fails only for certified
//! exceptional indices of the valuation grading, conditions
//! `sum c_(i,j)(n) [x^j] Q_(n+i) = 0` (`i <= R`, `j <= k`, `deg c <= U_c`)
//! pinning `[x^k]` of the row `n + R` are searched in the same way. Under a
//! constant-index monic source the truncations `P_n mod x^(K+1)` satisfy a
//! constant linear recurrence of order `r (K+1)`, so a condition with index
//! degree `U_c` vanishes for all bases once it vanishes at
//! `r (K'+1)(U_c+e+1)` consecutive bases, where `K' = k` (identity, factorial,
//! weights, Möbius-class) or `K' = k + m` (derivative), and `e` is the index
//! degree of `Lambda` for a weight (`[x^j] Q_n = c Lambda(n,j)/(j!)^m [x^j] P_n`)
//! and zero otherwise. For a Möbius-class map the truncations of `Q` itself
//! are used; this needs every `a~_i` to be regular at `x = 0` (power series),
//! which is checked exactly on the transported source rule. Otherwise (for
//! example a reversal with `deg a_i > h (k - i)`) conditions are finite fits.
//! Found conditions enter the description as boundary conditions with their
//! own evidence.
//!
//! [`verify_implicit_relation_search`] reruns the search from the embedded
//! request, requires equality, re-checks every kernel certificate against
//! residual matrices rebuilt from the reported rows, every relation and
//! condition on the rows, every division certificate and every assessment.
//!
//! Related modules: `fit::space` (complete spaces, kernel certificates),
//! `fit::rescaling`, `implicit::relation` (relations, descriptions),
//! `implicit::lossy` (the derivative map).

use super::algebra::{Biv, DOp};
use super::relation::{decode_poly, encode_poly, valuation_exceptional_indices};
use super::*;
use crate::linalg::{rational_kernel_with_certificate, RationalKernelError};
use crate::polynomial::Polynomial;
use crate::recurrence::fit::{
    check_kernel_certificate, encode_matrix, kernel_limits, residual_matrix,
};
use crate::vec_poly;
use num_traits::Zero;
use serde::{Deserialize, Serialize};

type Q = BigRational;
type Poly = Polynomial<Q>;

/// Schema tag of [`ImplicitRelationSearchReport`].
pub const RECURRENCE_IMPLICIT_SEARCH_SCHEMA: &str = "polytool.recurrence-implicit-search.v1";

/// Largest number of representatives reported per support.
pub const MAX_IMPLICIT_SEARCH_REPRESENTATIVES: usize = 8;

/// Largest exceptional coefficient index for which a condition is searched.
pub const MAX_CONDITION_COEFFICIENT_INDEX: usize = 16;

/// Largest number of supports in one request.
pub const MAX_IMPLICIT_SEARCH_SUPPORTS: usize = 16;

/// Number of rows generated from a source when `row_count` is omitted.
pub const DEFAULT_IMPLICIT_SEARCH_ROW_COUNT: usize = 12;

/// Support role: no derivative on the leader row.
pub const SUCCESSOR_SUPPORT: &str = "successor";

/// Support role: derivatives allowed on the leader row.
pub const IMPLICIT_SUCCESSOR_SUPPORT: &str = "implicit_successor";

/// Support role: one row (`max_shift = 0`).
pub const ROW_EQUATION_SUPPORT: &str = "row_equation";

/// A box support of the ansatz `sum c_(i,j)(n,x) D^j Q_(n+i) = 0`.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ImplicitSearchSupport {
    /// `successor` (no derivative on the row `n + max_shift`),
    /// `implicit_successor` (derivatives on every row) or `row_equation`
    /// (`max_shift = 0`).
    pub role: String,
    /// Largest shift `R`; positive for the successor roles, zero for a row
    /// equation.
    #[serde(default)]
    pub max_shift: usize,
    /// Largest derivative order `J` on any row.
    #[serde(default)]
    pub derivative_order: usize,
    /// Largest power `U` of the index `n` in a coefficient.
    #[serde(default)]
    pub index_degree: usize,
    /// Largest power `A` of `x` in a coefficient.
    #[serde(default)]
    pub variable_degree: usize,
}

impl ImplicitSearchSupport {
    /// The columns `n^u x^v D^s Q_(n+j)` of the support, sorted.
    pub fn monomials(
        &self,
        limits: &RecurrenceSpaceLimits,
    ) -> Result<Vec<RecurrenceMonomial>, RecurrenceError> {
        match self.role.as_str() {
            SUCCESSOR_SUPPORT | IMPLICIT_SUCCESSOR_SUPPORT if self.max_shift == 0 => {
                return Err(RecurrenceError::InvalidInput(format!(
                    "a {} support needs max_shift >= 1",
                    self.role
                )))
            }
            ROW_EQUATION_SUPPORT if self.max_shift != 0 => {
                return Err(RecurrenceError::InvalidInput(
                    "a row_equation support has max_shift = 0".into(),
                ))
            }
            SUCCESSOR_SUPPORT | IMPLICIT_SUCCESSOR_SUPPORT | ROW_EQUATION_SUPPORT => {}
            other => {
                return Err(RecurrenceError::InvalidInput(format!(
                    "unknown support role `{other}` (successor, implicit_successor, row_equation)"
                )))
            }
        }
        if [
            self.derivative_order,
            self.index_degree,
            self.variable_degree,
        ]
        .iter()
        .any(|d| *d > limits.max_degree)
            || self.max_shift > limits.max_rows
        {
            return Err(RecurrenceError::BudgetExhausted(
                "support degree budget".into(),
            ));
        }
        let count = (self.max_shift + 1)
            .checked_mul(self.derivative_order + 1)
            .and_then(|c| c.checked_mul(self.index_degree + 1))
            .and_then(|c| c.checked_mul(self.variable_degree + 1));
        if count.is_none_or(|c| c > limits.max_columns) {
            return Err(RecurrenceError::BudgetExhausted(
                "support column budget".into(),
            ));
        }
        let mut cols = Vec::new();
        for lag in 0..=self.max_shift {
            let top = if lag == self.max_shift && self.role == SUCCESSOR_SUPPORT {
                0
            } else {
                self.derivative_order
            };
            for derivative_order in 0..=top {
                for index_degree in 0..=self.index_degree {
                    for variable_degree in 0..=self.variable_degree {
                        cols.push(RecurrenceMonomial {
                            lag,
                            derivative_order,
                            index_degree,
                            variable_degree,
                        });
                    }
                }
            }
        }
        cols.sort();
        Ok(cols)
    }
}

/// The supports used when a request lists none: an implicit successor
/// relation of shift 1 (`J = U = A = 1`) and a row equation of order 2
/// (`U = A = 1`).
pub fn default_implicit_search_supports() -> Vec<ImplicitSearchSupport> {
    vec![
        ImplicitSearchSupport {
            role: IMPLICIT_SUCCESSOR_SUPPORT.into(),
            max_shift: 1,
            derivative_order: 1,
            index_degree: 1,
            variable_degree: 1,
        },
        ImplicitSearchSupport {
            role: ROW_EQUATION_SUPPORT.into(),
            max_shift: 0,
            derivative_order: 2,
            index_degree: 1,
            variable_degree: 1,
        },
    ]
}

fn default_true() -> bool {
    true
}

fn is_true(value: &bool) -> bool {
    *value
}

fn default_condition_index_degree() -> usize {
    1
}

fn is_one(value: &usize) -> bool {
    *value == 1
}

/// Request for [`search_implicit_relations`].
///
/// Give a defining rule (`source_recurrence`), a bundled OEIS family
/// (`oeis_id`) or rows only (`rows`, `row_origin`); with a source, optional
/// `rows` must agree with the generated rows.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ImplicitRelationSearchRequest {
    /// Defining rule with seeds; its `first_index` is the index origin.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_recurrence: Option<RecurrenceJson>,
    /// Bundled OEIS family (displayed row indices), used as the defining rule.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub oeis_id: Option<String>,
    /// Source rows `P_(row_origin + i)`, ascending rational strings.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rows: Option<Vec<Vec<String>>>,
    /// Original index of `rows[0]`; defaults to the source's first index, else 0.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub row_origin: Option<usize>,
    /// Least number of rows generated from a source (default 12); more are
    /// generated when a completeness theorem needs them.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub row_count: Option<usize>,
    /// Optional registry map `Q = T(P)` applied to every row before the
    /// search (scalar maps only).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transformation: Option<RecurrenceTransformation>,
    /// First base index `n` of every residual; defaults to the row origin.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub domain_start: Option<usize>,
    /// Supports to search, one result each; empty means
    /// [`default_implicit_search_supports`].
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub supports: Vec<ImplicitSearchSupport>,
    /// Search coefficient conditions for exceptional indices (default `true`).
    #[serde(default = "default_true", skip_serializing_if = "is_true")]
    pub conditions: bool,
    /// Largest power of `n` in a searched condition coefficient (default 1).
    #[serde(
        default = "default_condition_index_degree",
        skip_serializing_if = "is_one"
    )]
    pub condition_index_degree: usize,
    /// Resource limits; omitted fields take the defaults.
    #[serde(default)]
    pub limits: RecurrenceSpaceLimits,
}

impl ImplicitRelationSearchRequest {
    /// A request for the bundled OEIS family `id`, optionally after a map.
    pub fn from_oeis(id: &str, transformation: Option<RecurrenceTransformation>) -> Self {
        Self {
            source_recurrence: None,
            oeis_id: Some(id.to_string()),
            rows: None,
            row_origin: None,
            row_count: None,
            transformation,
            domain_start: None,
            supports: Vec::new(),
            conditions: true,
            condition_index_degree: 1,
            limits: RecurrenceSpaceLimits::default(),
        }
    }
}

/// The resolved source of a search.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ImplicitSearchSource {
    /// `defining_rule`, `oeis` or `rows_only`.
    pub kind: String,
    /// The defining rule with seeds (original indices), when given.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rule: Option<RecurrenceJson>,
    /// Original index of the first row.
    pub row_origin: usize,
    /// Order `k` when the rule is a constant-index monic pure-lag rule valid
    /// from its origin (the hypothesis of the completeness theorems).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub constant_index_order: Option<usize>,
    /// The finite verification theorem used, or why none applies.
    pub completeness_basis: String,
}

/// One column `n^index_degree [x^coefficient_index] Q_(n+shift)` of a
/// condition ansatz.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(deny_unknown_fields)]
pub struct ConditionMonomial {
    /// Forward shift of the row.
    pub shift: usize,
    /// The coefficient index `j`.
    pub coefficient_index: usize,
    /// Power of `n`.
    pub index_degree: usize,
}

/// The search for a condition pinning one exceptional coefficient index.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ImplicitConditionSearch {
    /// The exceptional coefficient index `k` of the row `n + R`.
    pub coefficient_index: usize,
    /// `condition_found`, `no_covering_condition`, `absent_in_support`,
    /// `insufficient_rows` or `incomplete_budget`.
    pub status: String,
    /// `certified_complete_space`, `finite_fit`, `rows_only_obstruction` or `none`.
    pub evidence: String,
    /// Explanation when incomplete or not found.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reason: Option<String>,
    /// Columns of the condition matrix, sorted.
    pub monomials: Vec<ConditionMonomial>,
    /// Number of consecutive bases from the domain start that are tested.
    pub residual_count: usize,
    /// Exact kernel and rank witnesses of the condition matrix.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub certificate: Option<RecurrenceSpaceCertificate>,
    /// Index into `certificate.basis` of the chosen condition.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub chosen_basis_index: Option<usize>,
    /// The chosen condition (also a boundary condition of the description).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub condition: Option<DescribedCondition>,
}

/// One representative relation of a support.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ImplicitSearchRepresentative {
    /// Index of the kernel vector in the certificate's `basis`.
    pub basis_index: usize,
    /// The relation in canonical form with its role, domain and evidence
    /// (`certified_complete_space`, `certified_by_division` or `finite_fit`).
    pub relation: DescribedRelation,
    /// `M o T = sum_nu H_nu o S^nu o L` for the defining rule `L`, when the
    /// relation is certified by division.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub division_certificate: Option<ReconstructionCertificate>,
    /// Condition searches for the exceptional indices of the leader.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub condition_searches: Vec<ImplicitConditionSearch>,
    /// Seeds, the relation and the found conditions.
    pub description: GeneratingDescription,
    /// Assessment of `description`.
    pub assessment: DescriptionAssessment,
}

/// The result of one support.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ImplicitSupportResult {
    /// The support.
    pub support: ImplicitSearchSupport,
    /// `relations_found`, `absent_in_support`, `insufficient_rows` or
    /// `incomplete_budget`.
    pub status: String,
    /// `certified_complete_space`, `finite_fit`, `rows_only_obstruction` or
    /// `none` (`none` never means absence).
    pub evidence: String,
    /// Explanation when incomplete.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reason: Option<String>,
    /// Columns of the residual matrix, sorted.
    pub monomials: Vec<RecurrenceMonomial>,
    /// The exponent `h` of the finite theorem, when one applies.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub annihilator_power: Option<usize>,
    /// Number of consecutive bases from the domain start that are tested.
    pub residual_count: usize,
    /// Original index of the last row used.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub last_required_row: Option<usize>,
    /// Number of scalar equations (`x`-coefficients of the residuals).
    pub equations: usize,
    /// Dimension of the kernel (all relations in the support, leaders or not).
    pub kernel_dimension: usize,
    /// Rank of the kernel projected on the leader columns (row `n + R`;
    /// every column for a row equation).
    pub leader_projection_rank: usize,
    /// Exact kernel and rank witnesses; absent when incomplete.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub certificate: Option<RecurrenceSpaceCertificate>,
    /// The completeness statement, when the theorem applies.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub completeness: Option<String>,
    /// Representatives (at most [`MAX_IMPLICIT_SEARCH_REPRESENTATIVES`]).
    pub representatives: Vec<ImplicitSearchRepresentative>,
}

/// Versioned report of [`search_implicit_relations`]
/// (`polytool.recurrence-implicit-search.v1`).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ImplicitRelationSearchReport {
    /// `polytool.recurrence-implicit-search.v1`.
    pub schema: String,
    /// The request, retained for replay.
    pub request: ImplicitRelationSearchRequest,
    /// The resolved source.
    pub source: ImplicitSearchSource,
    /// The map applied before the search, when given.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transformation: Option<TransformationRecord>,
    /// Row name in relation texts: `P` (source rows) or `Q` (mapped rows).
    pub symbol: String,
    /// First base index of every residual.
    pub domain_start: usize,
    /// The searched rows `[row_origin + i]` (after the map), ascending
    /// rational strings.
    pub rows: Vec<Vec<String>>,
    /// One result per support, in request order.
    pub results: Vec<ImplicitSupportResult>,
    /// Further explanations.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub notes: Vec<String>,
}

// ---------------------------------------------------------------------------
// Context
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq)]
enum Theorem {
    Identity,
    Factorial(usize),
    Derivative(usize),
    /// A Möbius-class map: `Q` satisfies a monic constant-index recurrence
    /// of the source order over `Q(x)`; `regular_at_zero` when its
    /// coefficients have no pole at `x = 0` (needed for conditions).
    Transported {
        regular_at_zero: bool,
    },
    /// A diagonal weight `scale (k!)^(-m) prod f_j(n,k)^(e_j)` with every
    /// `e_j > 0`: `total` and `index` are the total `(n,k)` degree and the
    /// `n` degree of the polynomial part.
    Weight {
        m: usize,
        total: usize,
        index: usize,
    },
}

/// `(total degree, n degree)` of the nonzero entries of `coeffs[a][b]`
/// (`n^a k^b`); `None` for a zero or unparsable polynomial.
fn weight_factor_degrees(coeffs: &BivarPolyJson) -> Option<(usize, usize)> {
    let mut degrees = None;
    for (a, row) in coeffs.coeffs.iter().enumerate() {
        for (b, c) in row.iter().enumerate() {
            if !parse_rational_coeff(c).ok()?.is_zero() {
                let (t, i) = degrees.unwrap_or((0, 0));
                degrees = Some((t.max(a + b), i.max(a)));
            }
        }
    }
    degrees
}

/// The weight theorem of a diagonal map whose weight is a polynomial in
/// `(n, k)` times `(k!)^(-m)`.
fn weight_theorem(weight: &DiagonalWeight) -> Option<Theorem> {
    match weight {
        DiagonalWeight::Binomial { m } => Some(Theorem::Weight {
            m: 0,
            total: *m,
            index: 0,
        }),
        DiagonalWeight::Product {
            factorial_exponent,
            factors,
            ..
        } => {
            if *factorial_exponent > 0 || factors.iter().any(|f| f.exponent <= 0) {
                return None;
            }
            let m = usize::try_from(factorial_exponent.unsigned_abs()).ok()?;
            let (mut total, mut index) = (0usize, 0usize);
            for f in factors {
                let e = usize::try_from(f.exponent).ok()?;
                let (t, i) = weight_factor_degrees(&f.coeffs)?;
                total = total.checked_add(t.checked_mul(e)?)?;
                index = index.checked_add(i.checked_mul(e)?)?;
            }
            Some(Theorem::Weight { m, total, index })
        }
        _ => None,
    }
}

/// Whether the coefficients `a~_i = c_i / f` of the transported source rule
/// `f Q_m = sum_i c_i Q_(m-i)` are regular at `x = 0`; `None` when the
/// transported rule is not index free (it always is for a Möbius-class map
/// of a constant-index pure-lag source).
fn transported_regular_at_zero(
    t: &RecurrenceTransformation,
    json: &RecurrenceJson,
) -> Option<bool> {
    let (rule, _, _) = json.to_recurrence_parts().ok()?;
    let TransportedRecurrence::Scalar(transported) = t.transport(&rule).ok()? else {
        return None;
    };
    let valuation = |p: &BivarPoly| -> Option<Option<usize>> {
        if p.coeffs.iter().skip(1).flatten().any(|c| !c.is_zero()) {
            return None;
        }
        Some(
            p.coeffs
                .first()
                .and_then(|row| row.iter().position(|c| !c.is_zero())),
        )
    };
    let leader = match &transported.denominator {
        Some(f) => valuation(f)?.ok_or(()).ok()?,
        None => 0,
    };
    for term in &transported.terms {
        if term.deriv_order != 0 || term.is_alternating() {
            return None;
        }
        if valuation(&term.coeff)?.is_some_and(|v| v < leader) {
            return Some(false);
        }
    }
    Some(true)
}

/// Whether a finite verification theorem of this module covers the map `t`
/// for the source rule `json` (a constant-index monic pure-lag rule valid
/// from its first index).
pub(in crate::recurrence) fn finite_theorem_applies(
    t: &RecurrenceTransformation,
    json: &RecurrenceJson,
) -> bool {
    constant_order(json).is_some() && map_theorem(Some(t), json).is_some()
}

/// The finite theorem for the map `t` (`None`: the identity) of a
/// constant-index monic pure-lag source `json`, if one is proved.
fn map_theorem(t: Option<&RecurrenceTransformation>, json: &RecurrenceJson) -> Option<Theorem> {
    use RecurrenceTransformation as T;
    let Some(t) = t else {
        return Some(Theorem::Identity);
    };
    match t {
        T::Diagonal {
            weight: DiagonalWeight::FactorialPower { power },
        } => Some(Theorem::Factorial(*power)),
        T::Diagonal { weight } => weight_theorem(weight),
        T::Derivative { order } => Some(Theorem::Derivative(*order)),
        T::MultiplyPower { .. } | T::DividePower { .. } | T::Affine { .. } => {
            Some(Theorem::Transported {
                regular_at_zero: true,
            })
        }
        T::Reversal { schedule } | T::Mobius { schedule, .. } | T::FPolynomial { schedule } => {
            if schedule.last_index.is_some() {
                return None;
            }
            Some(Theorem::Transported {
                regular_at_zero: transported_regular_at_zero(t, json)?,
            })
        }
        _ => None,
    }
}

struct Context<'a> {
    request: &'a ImplicitRelationSearchRequest,
    origin: usize,
    start: usize,
    rule: Option<(Recurrence, Vec<Vec<Q>>)>,
    order: Option<usize>,
    theorem: Option<Theorem>,
    source_rows: Vec<Vec<Q>>,
    rows: Vec<Vec<Q>>,
    defining: Option<(ImplicitRelation, usize)>,
    forward: Option<ImplicitRelation>,
    symbol: &'static str,
}

impl Context<'_> {
    /// Make at least `count` rows available (from the origin) when a rule
    /// generates them; rows-only contexts keep their rows.
    fn ensure(&mut self, count: usize) -> Result<(), RecurrenceError> {
        if count <= self.rows.len() {
            return Ok(());
        }
        let Some((rule, seeds)) = &self.rule else {
            return Ok(());
        };
        if count > self.request.limits.max_rows {
            return Err(RecurrenceError::BudgetExhausted("row budget".into()));
        }
        let generated = rule
            .generate_rows_rational(seeds, self.origin, count.max(seeds.len()))
            .map_err(|e| {
                RecurrenceError::InvalidInput(format!("source rule generation failed: {e}"))
            })?;
        for row in generated.iter().skip(self.source_rows.len()) {
            let row = vec_poly::trim(row);
            if vec_poly::degree(&row).unwrap_or(0) > self.request.limits.max_degree {
                return Err(RecurrenceError::BudgetExhausted("row degree budget".into()));
            }
            let index = self.origin + self.source_rows.len();
            let mapped = map_row(self.request.transformation.as_ref(), index, &row)?;
            self.source_rows.push(row);
            self.rows.push(mapped);
        }
        Ok(())
    }

    fn polys(&self) -> Vec<Poly> {
        self.rows.iter().map(|r| Poly::new(r.clone())).collect()
    }
}

fn map_row(
    t: Option<&RecurrenceTransformation>,
    index: usize,
    row: &[Q],
) -> Result<Vec<Q>, RecurrenceError> {
    let Some(t) = t else {
        return Ok(vec_poly::trim(row));
    };
    let mut state = t
        .map_row(index, row)
        .map_err(|e| RecurrenceError::InvalidInput(format!("{}: {e}", t.kind_name())))?;
    if state.len() != 1 {
        return Err(RecurrenceError::InvalidInput(format!(
            "{} gives {} components; the search needs scalar rows",
            t.kind_name(),
            state.len()
        )));
    }
    Ok(vec_poly::trim(&state.remove(0)))
}

/// `T` with `Q_n = T P_n` as a shift-free differential operator, when the map
/// is one (identity, `derivative`, a diagonal map with polynomial weight).
fn forward_operator(t: Option<&RecurrenceTransformation>) -> Option<ImplicitRelation> {
    let op = match t {
        None => DOp::single(0, 0, Biv::one()),
        Some(RecurrenceTransformation::Derivative { order }) => DOp::single(0, *order, Biv::one()),
        Some(RecurrenceTransformation::Diagonal { weight }) => {
            let reciprocal = match weight {
                DiagonalWeight::Binomial { m } => DiagonalWeight::InverseBinomial { m: *m },
                DiagonalWeight::Product {
                    scale,
                    factorial_exponent: 0,
                    factors,
                } if factors.iter().all(|f| f.exponent > 0) => {
                    let scale = parse_rational_coeff(scale).ok()?;
                    if scale.is_zero() {
                        return None;
                    }
                    DiagonalWeight::Product {
                        scale: format_rational_coeff(&scale.recip()),
                        factorial_exponent: 0,
                        factors: factors
                            .iter()
                            .map(|f| DiagonalWeightFactor {
                                coeffs: f.coeffs.clone(),
                                exponent: -f.exponent,
                            })
                            .collect(),
                    }
                }
                _ => return None,
            };
            return reciprocal.inverse_operator();
        }
        Some(_) => return None,
    };
    Some(ImplicitRelation::from_op(&op))
}

/// Whether the rule, its seeds and its extra initial rows form a
/// constant-index monic pure-lag source valid from the first index.
fn constant_order(json: &RecurrenceJson) -> Option<usize> {
    let (source, rows) = ConstantLagSource::from_recurrence_json(json).ok()?;
    let k = source.coefficients.len();
    let a = source
        .coefficients
        .iter()
        .map(|c| decode_poly(c))
        .collect::<Result<Vec<_>, _>>()
        .ok()?;
    let rows = rows
        .iter()
        .map(|r| decode_poly(r))
        .collect::<Result<Vec<_>, _>>()
        .ok()?;
    for j in k..rows.len() {
        let mut value = vec![Q::zero()];
        for (i, c) in a.iter().enumerate() {
            value = vec_poly::add(&value, &vec_poly::mul(c, &rows[j - k + i]));
        }
        if vec_poly::trim(&value) != vec_poly::trim(&rows[j]) {
            return None;
        }
    }
    Some(k)
}

fn resolve<'a>(
    request: &'a ImplicitRelationSearchRequest,
) -> Result<(Context<'a>, ImplicitSearchSource, Vec<String>), RecurrenceError> {
    let limits = &request.limits;
    let mut notes = Vec::new();
    if let Some(t) = &request.transformation {
        t.validate()
            .map_err(|e| RecurrenceError::InvalidInput(format!("transformation: {e}")))?;
        if t.component_count() != 1 {
            return Err(RecurrenceError::InvalidInput(format!(
                "{} gives vector rows; the search needs a scalar map",
                t.kind_name()
            )));
        }
    }
    let explicit = request
        .rows
        .as_ref()
        .map(|rows| {
            rows.iter()
                .map(|r| decode_poly(r))
                .collect::<Result<Vec<_>, _>>()
        })
        .transpose()?;
    let defining = match (&request.source_recurrence, &request.oeis_id) {
        (Some(_), Some(_)) => {
            return Err(RecurrenceError::InvalidInput(
                "give at most one of source_recurrence and oeis_id".into(),
            ))
        }
        (Some(json), None) => {
            let (rule, first, seeds) = json.to_recurrence_parts()?;
            Some(("defining_rule", rule, first, seeds))
        }
        (None, Some(id)) => {
            let definition = crate::oeis::by_id(id).ok_or_else(|| {
                RecurrenceError::InvalidInput(format!("unknown OEIS family `{id}`"))
            })?;
            let (rule, seeds, first_row) = definition.recurrence_export_parts()?;
            let first = usize::try_from(first_row).map_err(|_| {
                RecurrenceError::InvalidInput(format!(
                    "OEIS family `{id}` starts at a negative row index"
                ))
            })?;
            Some(("oeis", rule, first, seeds))
        }
        (None, None) => None,
    };
    let row_count = request
        .row_count
        .unwrap_or(DEFAULT_IMPLICIT_SEARCH_ROW_COUNT);
    if row_count > limits.max_rows {
        return Err(RecurrenceError::LimitExceeded(format!(
            "{row_count} rows exceed max_rows = {}",
            limits.max_rows
        )));
    }
    let (kind, rule, origin, json) = match defining {
        Some((kind, rule, first, seeds)) => {
            if rule.terms.iter().any(|t| t.offset == 0) {
                return Err(RecurrenceError::InvalidInput(
                    "source rule terms need positive offsets".into(),
                ));
            }
            if seeds.is_empty() || rule.max_offset() > seeds.len() {
                return Err(RecurrenceError::InvalidInput(
                    "source rule needs at least max-offset seeds".into(),
                ));
            }
            if request.row_origin.is_some_and(|o| o != first) {
                return Err(RecurrenceError::InvalidInput(format!(
                    "row_origin differs from the source first index {first}"
                )));
            }
            let json = RecurrenceJson::from_recurrence_rational(&rule, first, &seeds, None);
            (kind, Some((rule, seeds)), first, Some(json))
        }
        None => ("rows_only", None, request.row_origin.unwrap_or(0), None),
    };
    let start = request.domain_start.unwrap_or(origin);
    if start < origin {
        return Err(RecurrenceError::InvalidInput(
            "domain_start precedes the row origin".into(),
        ));
    }
    let order = json.as_ref().and_then(constant_order);
    let theorem = match (&json, order) {
        (Some(json), Some(_)) => map_theorem(request.transformation.as_ref(), json),
        _ => None,
    };
    let completeness_basis = match (theorem, order, &json) {
        (Some(Theorem::Identity), Some(k), _) => format!(
            "constant-index monic pure-lag source of order {k}: chi^(u+s+1) annihilates every \
             column x^v n^u D^s P[n+j]"
        ),
        (Some(Theorem::Factorial(m)), Some(k), _) => format!(
            "constant-index monic pure-lag source of order {k} under T(x^j) = x^j/(j!)^{m}: \
             H = U + {m}(A+J) + 1 (rescaling theorem)"
        ),
        (Some(Theorem::Derivative(m)), Some(k), _) => format!(
            "constant-index monic pure-lag source of order {k} under Q = D^{m} P: \
             chi^(u+s+{m}+1) annihilates every column x^v n^u D^s Q[n+j]"
        ),
        (Some(Theorem::Transported { regular_at_zero }), Some(k), _) => format!(
            "constant-index monic pure-lag source of order {k} under a Möbius-class map: Q \
             satisfies the transported monic rule Q[n+{k}] = sum a~_i(x) Q[n+i] with \
             coefficients in Q(x) independent of n, so chi~^(u+s+1) annihilates every column \
             x^v n^u D^s Q[n+j] (transported-source theorem){}",
            if regular_at_zero {
                "; the a~_i are regular at x = 0, so coefficient conditions are complete too"
            } else {
                "; an a~_i has a pole at x = 0, so coefficient conditions are finite fits"
            }
        ),
        (Some(Theorem::Weight { m, total, index }), Some(k), _) => format!(
            "constant-index monic pure-lag source of order {k} under a diagonal weight \
             (k!)^(-{m}) times a polynomial of total degree {total} (index degree {index}) in \
             (n, k): H = {} + {total} + 1 (weight theorem)",
            if m == 0 {
                "max(u+s)".to_string()
            } else {
                format!("U + {m}(A+J)")
            }
        ),
        (None, Some(_), _) => "no finite theorem for this map: results are finite fits unless \
                               certified by division"
            .to_string(),
        (_, None, Some(_)) => "the rule is not a constant-index monic pure-lag rule valid from \
                               its origin: results are finite fits unless certified by division"
            .to_string(),
        _ => "rows only: every relation is a finite fit; absence in a support is still \
              certified by the rows"
            .to_string(),
    };
    let defining_relation = match &rule {
        Some((r, seeds)) => match ImplicitRelation::from_recurrence(r) {
            Ok((relation, shift)) => {
                let first_base = (origin + seeds.len()).checked_sub(shift).ok_or(
                    RecurrenceError::InvalidInput("the rule starts before its leader shift".into()),
                )?;
                Some((relation, first_base))
            }
            Err(error) => {
                notes.push(format!(
                    "the defining rule has no implicit form ({error}); no division certificates"
                ));
                None
            }
        },
        None => None,
    };
    let forward = forward_operator(request.transformation.as_ref());
    let mut context = Context {
        request,
        origin,
        start,
        rule,
        order,
        theorem,
        source_rows: Vec::new(),
        rows: Vec::new(),
        defining: defining_relation,
        forward,
        symbol: if request.transformation.is_some() {
            "Q"
        } else {
            "P"
        },
    };
    match &explicit {
        Some(rows) if context.rule.is_none() => {
            if rows.is_empty() {
                return Err(RecurrenceError::InvalidInput(
                    "at least one row is required".into(),
                ));
            }
            if rows.len() > limits.max_rows {
                return Err(RecurrenceError::LimitExceeded(format!(
                    "{} rows exceed max_rows = {}",
                    rows.len(),
                    limits.max_rows
                )));
            }
            for (i, row) in rows.iter().enumerate() {
                if vec_poly::degree(row).unwrap_or(0) > limits.max_degree {
                    return Err(RecurrenceError::LimitExceeded(
                        "row degree exceeds max_degree".into(),
                    ));
                }
                let mapped = map_row(request.transformation.as_ref(), origin + i, row)?;
                context.source_rows.push(vec_poly::trim(row));
                context.rows.push(mapped);
            }
        }
        Some(rows) => {
            context.ensure(row_count.max(rows.len()))?;
            if context.source_rows[..rows.len()]
                .iter()
                .zip(rows)
                .any(|(a, b)| *a != vec_poly::trim(b))
            {
                return Err(RecurrenceError::InvalidInput(
                    "rows disagree with the rows generated by the source rule".into(),
                ));
            }
        }
        None if context.rule.is_none() => {
            return Err(RecurrenceError::InvalidInput(
                "give source_recurrence, oeis_id or rows".into(),
            ))
        }
        None => context.ensure(row_count)?,
    }
    let source = ImplicitSearchSource {
        kind: kind.to_string(),
        rule: json,
        row_origin: origin,
        constant_index_order: order,
        completeness_basis,
    };
    Ok((context, source, notes))
}

// ---------------------------------------------------------------------------
// Kernels
// ---------------------------------------------------------------------------

fn certificate_of(kernel: &crate::linalg::RationalKernel) -> RecurrenceSpaceCertificate {
    RecurrenceSpaceCertificate {
        basis: encode_matrix(&kernel.basis),
        free_columns: kernel.free_columns.clone(),
        pivot_columns: kernel.pivot_columns.clone(),
        independent_rows: kernel.independent_rows.clone(),
        inverse_minor: encode_matrix(&kernel.inverse_minor),
        absence_duals: Vec::new(),
    }
}

/// The homogeneous residual matrix (every column with the same sign).
fn homogeneous_matrix(
    rows: &[Poly],
    origin: usize,
    start: usize,
    count: usize,
    cols: &[RecurrenceMonomial],
    limits: &RecurrenceSpaceLimits,
) -> Result<Vec<Vec<Q>>, RecurrenceError> {
    // No column has lag usize::MAX, so every column enters negatively and the
    // kernel is that of the homogeneous ansatz.
    residual_matrix(rows, origin, start, count, usize::MAX, cols, limits)
}

/// The condition matrix: one row per base, one column per monomial.
fn condition_matrix(
    rows: &[Vec<Q>],
    origin: usize,
    start: usize,
    count: usize,
    cols: &[ConditionMonomial],
    limits: &RecurrenceSpaceLimits,
) -> Result<Vec<Vec<Q>>, RecurrenceError> {
    if count
        .checked_mul(cols.len())
        .is_none_or(|e| e > limits.max_matrix_entries)
    {
        return Err(RecurrenceError::BudgetExhausted(
            "condition matrix entry budget".into(),
        ));
    }
    (start..start + count)
        .map(|n| {
            cols.iter()
                .map(|c| {
                    let row =
                        rows.get(n + c.shift - origin)
                            .ok_or(RecurrenceError::InvalidInput(
                                "missing row for a condition".into(),
                            ))?;
                    let value = row
                        .get(c.coefficient_index)
                        .cloned()
                        .unwrap_or_else(Q::zero);
                    let power = u32::try_from(c.index_degree).map_err(|_| {
                        RecurrenceError::LimitExceeded("index degree overflow".into())
                    })?;
                    Ok(value * Q::from_integer(BigInt::from(n).pow(power)))
                })
                .collect()
        })
        .collect()
}

fn relation_from_vector(cols: &[RecurrenceMonomial], vector: &[Q]) -> ImplicitRelation {
    let mut op = DOp::default();
    for (c, value) in cols.iter().zip(vector) {
        if !value.is_zero() {
            op.add_term(
                (c.lag, c.derivative_order),
                Biv::monomial(c.index_degree, c.variable_degree, value.clone()),
            );
        }
    }
    ImplicitRelation::from_op(&op)
}

fn condition_from_vector(cols: &[ConditionMonomial], vector: &[Q]) -> CoefficientCondition {
    let mut terms = Vec::new();
    for (c, value) in cols.iter().zip(vector) {
        if value.is_zero() {
            continue;
        }
        let mut coeff = vec![Q::zero(); c.index_degree + 1];
        coeff[c.index_degree] = value.clone();
        terms.push(ConditionTerm {
            shift: c.shift,
            coefficient_index: c.coefficient_index,
            coeff,
        });
    }
    CoefficientCondition { terms }.canonical()
}

fn leader_columns(cols: &[RecurrenceMonomial], r: usize) -> Vec<usize> {
    cols.iter()
        .enumerate()
        .filter_map(|(i, c)| (r == 0 || c.lag == r).then_some(i))
        .collect()
}

/// Rank of the kernel basis projected on `leading`, and the basis indices of
/// an independent set of projections.
fn projection(
    basis: &[Vec<Q>],
    leading: &[usize],
    limits: &RecurrenceSpaceLimits,
) -> Result<Vec<usize>, RecurrenceError> {
    let matrix = leading
        .iter()
        .map(|&i| basis.iter().map(|v| v[i].clone()).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let kernel = rational_kernel_with_certificate(&matrix, basis.len(), kernel_limits(limits))
        .map_err(|_| RecurrenceError::BudgetExhausted("projection budget".to_string()))?;
    Ok(kernel.pivot_columns)
}

// ---------------------------------------------------------------------------
// The search
// ---------------------------------------------------------------------------

fn empty_result(support: &ImplicitSearchSupport) -> ImplicitSupportResult {
    ImplicitSupportResult {
        support: support.clone(),
        status: "incomplete_budget".into(),
        evidence: "none".into(),
        reason: None,
        monomials: Vec::new(),
        annihilator_power: None,
        residual_count: 0,
        last_required_row: None,
        equations: 0,
        kernel_dimension: 0,
        leader_projection_rank: 0,
        certificate: None,
        completeness: None,
        representatives: Vec::new(),
    }
}

fn annihilator_power(theorem: Theorem, cols: &[RecurrenceMonomial]) -> Option<usize> {
    let us = cols
        .iter()
        .map(|c| c.index_degree.checked_add(c.derivative_order))
        .collect::<Option<Vec<_>>>()?
        .into_iter()
        .max()
        .unwrap_or(0);
    match theorem {
        Theorem::Identity | Theorem::Transported { .. } => us.checked_add(1),
        Theorem::Derivative(m) => us.checked_add(m)?.checked_add(1),
        Theorem::Factorial(m) => {
            RecurrenceCoefficientRescaling::FactorialPower { power: m }.annihilator_power(cols)
        }
        Theorem::Weight { m: 0, total, .. } => us.checked_add(total)?.checked_add(1),
        Theorem::Weight { m, total, .. } => {
            RecurrenceCoefficientRescaling::FactorialPower { power: m }
                .annihilator_power(cols)?
                .checked_add(total)
        }
    }
}

/// `(K', e)` of the condition theorem: the truncation index `K'` (the
/// conditions on `[x^j] Q`, `j <= k`, depend on `P_n mod x^(K'+1)`) and the
/// extra index degree `e` contributed by the map; `None` when no condition
/// theorem applies.
fn condition_theorem(theorem: Theorem, k: usize) -> Option<(usize, usize)> {
    match theorem {
        Theorem::Identity | Theorem::Factorial(_) => Some((k, 0)),
        Theorem::Derivative(m) => Some((k.checked_add(m)?, 0)),
        Theorem::Transported { regular_at_zero } => regular_at_zero.then_some((k, 0)),
        Theorem::Weight { index, .. } => Some((k, index)),
    }
}

fn run_support(
    context: &mut Context<'_>,
    support: &ImplicitSearchSupport,
) -> Result<ImplicitSupportResult, RecurrenceError> {
    let limits = &context.request.limits;
    let mut result = empty_result(support);
    let cols = match support.monomials(limits) {
        Ok(cols) => cols,
        Err(error @ RecurrenceError::BudgetExhausted(_)) => {
            result.reason = Some(error.to_string());
            return Ok(result);
        }
        Err(error) => return Err(error),
    };
    result.monomials = cols.clone();
    let r = support.max_shift;
    let offset = context.start - context.origin;
    let (count, complete) = match (context.theorem, context.order) {
        (Some(theorem), Some(k)) => {
            let Some(h) = annihilator_power(theorem, &cols) else {
                result.reason = Some("annihilator power overflow".into());
                return Ok(result);
            };
            result.annihilator_power = Some(h);
            match h.checked_mul(k) {
                Some(count) => (count, true),
                None => {
                    result.reason = Some("residual bound overflow".into());
                    return Ok(result);
                }
            }
        }
        _ => (context.rows.len().saturating_sub(offset + r), false),
    };
    let needed = offset + count + r;
    if complete {
        if let Err(error) = context.ensure(needed) {
            if matches!(error, RecurrenceError::BudgetExhausted(_)) {
                result.reason = Some(error.to_string());
                return Ok(result);
            }
            return Err(error);
        }
    }
    result.residual_count = count;
    if count == 0 || context.rows.len() < needed {
        result.status = "insufficient_rows".into();
        result.reason = Some("no complete residual in the declared domain".into());
        return Ok(result);
    }
    result.last_required_row = Some(context.origin + needed - 1);
    let polys = context.polys();
    let matrix =
        match homogeneous_matrix(&polys, context.origin, context.start, count, &cols, limits) {
            Ok(m) => m,
            Err(error) => {
                result.reason = Some(error.to_string());
                return Ok(result);
            }
        };
    result.equations = matrix.len();
    let kernel = match rational_kernel_with_certificate(&matrix, cols.len(), kernel_limits(limits))
    {
        Ok(k) => k,
        Err(RationalKernelError::BudgetExceeded(e)) => {
            result.reason = Some(e.into());
            return Ok(result);
        }
        Err(_) => {
            return Err(RecurrenceError::InvalidInput(
                "invalid residual matrix".into(),
            ))
        }
    };
    let leading = leader_columns(&cols, r);
    let pivots = match projection(&kernel.basis, &leading, limits) {
        Ok(p) => p,
        Err(error) => {
            result.reason = Some(error.to_string());
            return Ok(result);
        }
    };
    result.kernel_dimension = kernel.basis.len();
    result.leader_projection_rank = pivots.len();
    result.certificate = Some(certificate_of(&kernel));
    if pivots.is_empty() {
        result.status = "absent_in_support".into();
        result.evidence = "rows_only_obstruction".into();
        return Ok(result);
    }
    result.status = "relations_found".into();
    result.evidence = if complete {
        CERTIFIED_COMPLETE_SPACE_EVIDENCE
    } else {
        FINITE_FIT_EVIDENCE
    }
    .into();
    if complete {
        result.completeness = Some(format!(
            "{}. Every relation in this support that vanishes at the {count} bases n = {}..{} \
             vanishes for all n >= {}; the kernel is the complete space of all-index relations \
             in the support.",
            context_basis(context),
            context.start,
            context.start + count - 1,
            context.start
        ));
    }
    for &index in pivots.iter().take(MAX_IMPLICIT_SEARCH_REPRESENTATIVES) {
        let representative =
            representative(context, &cols, &kernel.basis[index], index, r, complete)?;
        result.representatives.push(representative);
    }
    Ok(result)
}

fn context_basis(context: &Context<'_>) -> String {
    match (context.theorem, context.order) {
        (Some(Theorem::Identity), Some(k)) => {
            format!("P satisfies a constant-index monic pure-lag rule of order {k}")
        }
        (Some(Theorem::Factorial(m)), Some(k)) => format!(
            "Q = T(P), T(x^j) = x^j/(j!)^{m}, with P satisfying a constant-index monic pure-lag \
             rule of order {k}"
        ),
        (Some(Theorem::Derivative(m)), Some(k)) => format!(
            "Q = D^{m} P with P satisfying a constant-index monic pure-lag rule of order {k}"
        ),
        (Some(Theorem::Transported { .. }), Some(k)) => format!(
            "Q = T(P) for a Möbius-class map T, so Q satisfies the transported monic rule of \
             order {k} with index-free coefficients in Q(x)"
        ),
        (Some(Theorem::Weight { m, total, .. }), Some(k)) => format!(
            "Q_n = lambda(n, theta) P_n with lambda = (k!)^(-{m}) times a polynomial of total \
             degree {total}, P satisfying a constant-index monic pure-lag rule of order {k}"
        ),
        _ => String::new(),
    }
}

fn representative(
    context: &mut Context<'_>,
    cols: &[RecurrenceMonomial],
    vector: &[Q],
    basis_index: usize,
    r: usize,
    complete: bool,
) -> Result<ImplicitSearchRepresentative, RecurrenceError> {
    let relation = relation_from_vector(cols, vector).canonical();
    let symbol = context.symbol;
    let mut first_base = context.start;
    let mut evidence = if complete {
        CERTIFIED_COMPLETE_SPACE_EVIDENCE
    } else {
        FINITE_FIT_EVIDENCE
    };
    let mut division = None;
    if !complete {
        if let (Some((l, l_first)), Some(t)) = (&context.defining, &context.forward) {
            if let Some(certificate) = find_reconstruction_certificate(&relation, t, l) {
                certificate.verify_identity()?;
                evidence = CERTIFIED_BY_DIVISION_EVIDENCE;
                first_base = first_base.max(*l_first);
                division = Some(certificate);
            }
        }
    }
    // Sanity: the relation holds on every available row of its domain.
    let end = (context.origin + context.rows.len()).saturating_sub(relation.max_shift());
    for base in first_base..end {
        if relation.holds_at(&context.rows, context.origin, base) == Some(false) {
            if complete || division.is_some() {
                return Err(RecurrenceError::InvalidInput(format!(
                    "internal: a certified relation fails at base {base}"
                )));
            }
            return Err(RecurrenceError::InvalidInput(format!(
                "internal: a fitted relation fails at base {base}"
            )));
        }
    }
    let note = match evidence {
        CERTIFIED_COMPLETE_SPACE_EVIDENCE => "member of a certified complete space",
        CERTIFIED_BY_DIVISION_EVIDENCE => {
            "fitted, then proved by right division by the defining rule"
        }
        _ => "fitted to the rows (finite evidence)",
    };
    let described = DescribedRelation::new(&relation, first_base, evidence, symbol).with_note(note);
    let seed_end = (first_base + r).max(context.start + 1);
    if context.rows.len() < seed_end - context.origin {
        return Err(RecurrenceError::InvalidInput(
            "internal: missing seed rows".into(),
        ));
    }
    let seeds = context.rows[context.start - context.origin..seed_end - context.origin]
        .iter()
        .map(|row| encode_poly(row))
        .collect::<Vec<_>>();
    let mut description = GeneratingDescription {
        first_index: context.start,
        seeds,
        relations: vec![described.clone()],
        exceptional_index_conditions: Vec::new(),
        boundary_conditions: Vec::new(),
    };
    let mut assessment = assess_generating_description(&description)?;
    let mut searches = Vec::new();
    if context.request.conditions && r > 0 && assessment.uniqueness.status != "certified" {
        if let Some(indices) = valuation_exceptional_indices(&relation, first_base) {
            if !indices.is_empty()
                && indices
                    .iter()
                    .all(|&k| k <= MAX_CONDITION_COEFFICIENT_INDEX)
            {
                for k in indices {
                    let search = condition_search(context, k, r, first_base)?;
                    if let Some(condition) = &search.condition {
                        description.boundary_conditions.push(condition.clone());
                    }
                    searches.push(search);
                }
                assessment = assess_generating_description(&description)?;
            }
        }
    }
    Ok(ImplicitSearchRepresentative {
        basis_index,
        relation: described,
        division_certificate: division,
        condition_searches: searches,
        description,
        assessment,
    })
}

fn condition_columns(r: usize, k: usize, degree: usize) -> Vec<ConditionMonomial> {
    let mut cols = Vec::new();
    for shift in 0..=r {
        for coefficient_index in 0..=k {
            for index_degree in 0..=degree {
                cols.push(ConditionMonomial {
                    shift,
                    coefficient_index,
                    index_degree,
                });
            }
        }
    }
    cols.sort();
    cols
}

fn condition_search(
    context: &mut Context<'_>,
    k: usize,
    r: usize,
    first_base: usize,
) -> Result<ImplicitConditionSearch, RecurrenceError> {
    let degree = context.request.condition_index_degree;
    let cols = condition_columns(r, k, degree);
    let mut search = ImplicitConditionSearch {
        coefficient_index: k,
        status: "incomplete_budget".into(),
        evidence: "none".into(),
        reason: None,
        monomials: cols.clone(),
        residual_count: 0,
        certificate: None,
        chosen_basis_index: None,
        condition: None,
    };
    let limits = context.request.limits.clone();
    if cols.len() > limits.max_columns || degree > limits.max_degree {
        search.reason = Some("condition column budget".into());
        return Ok(search);
    }
    let offset = context.start - context.origin;
    let condition = context.theorem.and_then(|t| condition_theorem(t, k));
    let (count, complete) = match (condition, context.order) {
        (Some((truncation, extra)), Some(order)) => match truncation
            .checked_add(1)
            .and_then(|t| order.checked_mul(t))
            .and_then(|c| c.checked_mul(degree.checked_add(extra)?.checked_add(1)?))
        {
            Some(count) => (count, true),
            None => {
                search.reason = Some("condition residual bound overflow".into());
                return Ok(search);
            }
        },
        _ => (context.rows.len().saturating_sub(offset + r), false),
    };
    let needed = offset + count + r;
    if complete {
        if let Err(error) = context.ensure(needed) {
            if matches!(error, RecurrenceError::BudgetExhausted(_)) {
                search.reason = Some(error.to_string());
                return Ok(search);
            }
            return Err(error);
        }
    }
    search.residual_count = count;
    if count == 0 || context.rows.len() < needed {
        search.status = "insufficient_rows".into();
        search.reason = Some("no complete condition residual in the declared domain".into());
        return Ok(search);
    }
    let matrix = match condition_matrix(
        &context.rows,
        context.origin,
        context.start,
        count,
        &cols,
        &limits,
    ) {
        Ok(m) => m,
        Err(error) => {
            search.reason = Some(error.to_string());
            return Ok(search);
        }
    };
    let kernel = match rational_kernel_with_certificate(&matrix, cols.len(), kernel_limits(&limits))
    {
        Ok(kernel) => kernel,
        Err(RationalKernelError::BudgetExceeded(e)) => {
            search.reason = Some(e.into());
            return Ok(search);
        }
        Err(_) => {
            return Err(RecurrenceError::InvalidInput(
                "invalid condition matrix".into(),
            ))
        }
    };
    search.certificate = Some(certificate_of(&kernel));
    if kernel.basis.is_empty() {
        search.status = "absent_in_support".into();
        search.evidence = "rows_only_obstruction".into();
        return Ok(search);
    }
    let evidence = if complete {
        CERTIFIED_COMPLETE_SPACE_EVIDENCE
    } else {
        FINITE_FIT_EVIDENCE
    };
    // A basis vector whose free column is [x^k] of the row n + r has its top
    // shift terms at indices <= k; prefer the smallest power of n.
    for (i, &free) in kernel.free_columns.iter().enumerate() {
        let column = cols[free];
        if column.shift != r || column.coefficient_index != k {
            continue;
        }
        let condition = condition_from_vector(&cols, &kernel.basis[i]);
        let pins = condition
            .terms
            .iter()
            .find(|t| t.shift == r && t.coefficient_index == k)
            .is_some_and(|t| {
                super::algebra::roots_in_domain(&t.coeff, first_base, None).is_empty()
            });
        if !pins {
            continue;
        }
        let note = if complete {
            "member of a certified complete condition space (truncated source recurrence)"
        } else {
            "fitted to the rows (finite evidence)"
        };
        search.status = "condition_found".into();
        search.evidence = evidence.into();
        search.chosen_basis_index = Some(i);
        search.condition = Some(
            DescribedCondition::new(&condition, context.start, evidence, Some(k), context.symbol)
                .with_note(note),
        );
        return Ok(search);
    }
    search.status = "no_covering_condition".into();
    search.evidence = evidence.into();
    search.reason = Some(format!(
        "no condition in the support pins [x^{k}] of the row n + {r} with a coefficient that \
         is nonzero on the domain"
    ));
    Ok(search)
}

/// Search implicit relations of exact rows over prescribed supports.
///
/// Errors are returned for invalid requests (sources, maps, supports,
/// limits, contradictory rows); a limit reached inside one support makes
/// that result `incomplete_budget`, never absence.
///
/// ```
/// use polytool::recurrence::*;
///
/// // A007318 divided coefficientwise by k!: D Q[n+1] - D Q[n] - Q[n] = 0,
/// // complete by the rescaling theorem, with the constant condition.
/// let mut request = ImplicitRelationSearchRequest::from_oeis(
///     "A007318",
///     Some(RecurrenceTransformation::Diagonal {
///         weight: DiagonalWeight::FactorialPower { power: 1 },
///     }),
/// );
/// request.supports = vec![ImplicitSearchSupport {
///     role: "implicit_successor".into(),
///     max_shift: 1,
///     derivative_order: 1,
///     index_degree: 0,
///     variable_degree: 0,
/// }];
/// let report = search_implicit_relations(&request).unwrap();
/// let found = &report.results[0].representatives[0];
/// assert_eq!(found.relation.text, "D Q[n+1] - D Q[n] - Q[n] = 0");
/// assert_eq!(found.relation.evidence, CERTIFIED_COMPLETE_SPACE_EVIDENCE);
/// assert_eq!(found.description.boundary_conditions[0].text, "[x^0] Q[n+1] - [x^0] Q[n] = 0");
/// assert_eq!(found.assessment.status, EQUIVALENT_DESCRIPTION);
/// verify_implicit_relation_search(&report).unwrap();
/// ```
pub fn search_implicit_relations(
    request: &ImplicitRelationSearchRequest,
) -> Result<ImplicitRelationSearchReport, RecurrenceError> {
    if request.supports.len() > MAX_IMPLICIT_SEARCH_SUPPORTS {
        return Err(RecurrenceError::LimitExceeded(format!(
            "at most {MAX_IMPLICIT_SEARCH_SUPPORTS} supports are supported"
        )));
    }
    if request.condition_index_degree > 8 {
        return Err(RecurrenceError::LimitExceeded(
            "condition_index_degree is at most 8".into(),
        ));
    }
    let supports = if request.supports.is_empty() {
        default_implicit_search_supports()
    } else {
        request.supports.clone()
    };
    for support in &supports {
        support.monomials(&RecurrenceSpaceLimits {
            max_columns: usize::MAX,
            max_degree: usize::MAX,
            max_rows: usize::MAX,
            ..request.limits.clone()
        })?;
    }
    let (mut context, source, mut notes) = resolve(request)?;
    if context.theorem.is_none() && context.order.is_some() {
        notes.push(format!(
            "the source is constant-index monic, but no finite theorem covers the map {}",
            request
                .transformation
                .as_ref()
                .map_or("identity", |t| t.kind_name())
        ));
    }
    let mut results = Vec::new();
    for support in &supports {
        results.push(run_support(&mut context, support)?);
    }
    Ok(ImplicitRelationSearchReport {
        schema: RECURRENCE_IMPLICIT_SEARCH_SCHEMA.to_string(),
        request: request.clone(),
        source,
        transformation: request.transformation.as_ref().map(|t| t.record()),
        symbol: context.symbol.to_string(),
        domain_start: context.start,
        rows: context.rows.iter().map(|r| encode_poly(r)).collect(),
        results,
        notes,
    })
}

/// Replay a report: rerun the search from its embedded request and require
/// equality, then re-check every kernel certificate against matrices rebuilt
/// from the reported rows, every relation and condition on the rows, every
/// division certificate and every assessment.
pub fn verify_implicit_relation_search(
    report: &ImplicitRelationSearchReport,
) -> Result<(), RecurrenceError> {
    if report.schema != RECURRENCE_IMPLICIT_SEARCH_SCHEMA {
        return Err(RecurrenceError::InvalidInput(format!(
            "unsupported schema `{}`",
            report.schema
        )));
    }
    let replay = search_implicit_relations(&report.request)?;
    if &replay != report {
        return Err(RecurrenceError::VerificationFailed(
            "report does not match its exact replay".into(),
        ));
    }
    let limits = &report.request.limits;
    let origin = report.source.row_origin;
    let start = report.domain_start;
    let rows = report
        .rows
        .iter()
        .map(|r| decode_poly(r))
        .collect::<Result<Vec<_>, _>>()?;
    let polys = rows
        .iter()
        .map(|r| Poly::new(r.clone()))
        .collect::<Vec<_>>();
    for (s, result) in report.results.iter().enumerate() {
        let label = format!("results[{s}]");
        let Some(certificate) = &result.certificate else {
            if !result.representatives.is_empty() {
                return Err(RecurrenceError::VerificationFailed(format!(
                    "{label}: representatives without a certificate"
                )));
            }
            continue;
        };
        let matrix = homogeneous_matrix(
            &polys,
            origin,
            start,
            result.residual_count,
            &result.monomials,
            limits,
        )?;
        let basis = check_kernel_certificate(&matrix, result.monomials.len(), certificate, limits)
            .map_err(|e| RecurrenceError::VerificationFailed(format!("{label}: {e}")))?;
        let leading = leader_columns(&result.monomials, result.support.max_shift);
        if projection(&basis, &leading, limits)?.len() != result.leader_projection_rank {
            return Err(RecurrenceError::VerificationFailed(format!(
                "{label}: leader projection mismatch"
            )));
        }
        for (i, rep) in result.representatives.iter().enumerate() {
            let label = format!("{label}.representatives[{i}]");
            let vector = basis.get(rep.basis_index).ok_or_else(|| {
                RecurrenceError::VerificationFailed(format!("{label}: unknown basis index"))
            })?;
            let relation = ImplicitRelation::from_json(&rep.relation.relation)?;
            if !relation.same_relation(&relation_from_vector(&result.monomials, vector)) {
                return Err(RecurrenceError::VerificationFailed(format!(
                    "{label}: the relation is not its kernel vector"
                )));
            }
            let end = (origin + rows.len()).saturating_sub(relation.max_shift());
            for base in rep.relation.first_base..end {
                if relation.holds_at(&rows, origin, base) == Some(false) {
                    return Err(RecurrenceError::VerificationFailed(format!(
                        "{label}: the relation fails at base {base}"
                    )));
                }
            }
            if let Some(division) = &rep.division_certificate {
                division.verify_identity()?;
                if division.source_relation != relation.to_json() {
                    return Err(RecurrenceError::VerificationFailed(format!(
                        "{label}: the division certificate is for another relation"
                    )));
                }
            }
            for search in &rep.condition_searches {
                let Some(certificate) = &search.certificate else {
                    continue;
                };
                let matrix = condition_matrix(
                    &rows,
                    origin,
                    start,
                    search.residual_count,
                    &search.monomials,
                    limits,
                )?;
                let basis =
                    check_kernel_certificate(&matrix, search.monomials.len(), certificate, limits)
                        .map_err(|e| {
                            RecurrenceError::VerificationFailed(format!(
                                "{label}: condition search: {e}"
                            ))
                        })?;
                if let (Some(i), Some(condition)) = (search.chosen_basis_index, &search.condition) {
                    let vector = basis.get(i).ok_or_else(|| {
                        RecurrenceError::VerificationFailed(format!(
                            "{label}: unknown condition basis index"
                        ))
                    })?;
                    let parsed = CoefficientCondition::from_json(&condition.condition)?;
                    if parsed != condition_from_vector(&search.monomials, vector) {
                        return Err(RecurrenceError::VerificationFailed(format!(
                            "{label}: the condition is not its kernel vector"
                        )));
                    }
                    let end = (origin + rows.len()).saturating_sub(parsed.max_shift());
                    for base in condition.first_base..end {
                        if parsed.holds_at(&rows, origin, base) == Some(false) {
                            return Err(RecurrenceError::VerificationFailed(format!(
                                "{label}: a condition fails at base {base}"
                            )));
                        }
                    }
                }
            }
            if assess_generating_description(&rep.description)? != rep.assessment {
                return Err(RecurrenceError::VerificationFailed(format!(
                    "{label}: the assessment does not replay"
                )));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
