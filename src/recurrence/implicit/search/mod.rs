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
            .map_err(|e| e.with_context(&label))?;
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
                        .map_err(|e| e.with_context(format!("{label}: condition search")))?;
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

mod context;
pub(in crate::recurrence) use context::finite_theorem_applies;
use context::*;
mod kernel;
use kernel::*;
