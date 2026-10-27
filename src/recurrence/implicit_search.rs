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
//! `chi = S^k - sum a_i(x) S^i` valid from its origin and the map is the
//! identity, a factorial rescaling (`diagonal` with `factorial_power m`) or
//! plain differentiation (`derivative` of order `m`), a finite theorem bounds
//! the residuals that must vanish:
//!
//! - identity: every column `x^v n^u D^s P_(n+j)` is annihilated by
//!   `chi^(u+s+1)` (module `space`), so `h = max(u+s) + 1`;
//! - `Q = D^m P`: the column `x^v n^u D^s Q_(n+j)` is `x^v n^u D^(s+m) P_(n+j)`,
//!   so `h = max(u+s) + m + 1`;
//! - factorial rescaling: `h = U + m(A+J) + 1` (module `rescaling`).
//!
//! The residual is annihilated by the monic operator `chi^h` of order `k h`,
//! so its vanishing at the `k h` bases `n_0..n_0+kh-1` (`n_0` the domain
//! start, at least the origin) implies vanishing for every `n >= n_0`. The
//! kernel of those residual rows is then the complete space of all-index
//! relations in the support (evidence `certified_complete_space`). This holds
//! for derivatives on the leader row and for row equations alike. Otherwise a
//! representative is a `finite_fit`, unless it is proved by right division:
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
//! `r (K'+1)(U_c+1)` consecutive bases, where `K' = k` (identity, factorial)
//! or `K' = k + m` (derivative). Found conditions enter the description as
//! boundary conditions with their own evidence.
//!
//! [`verify_implicit_relation_search`] reruns the search from the embedded
//! request, requires equality, re-checks every kernel certificate against
//! residual matrices rebuilt from the reported rows, every relation and
//! condition on the rows, every division certificate and every assessment.
//!
//! Related modules: `space` (complete spaces, kernel certificates),
//! `rescaling`, `implicit` (relations, descriptions), `lossy` (the
//! derivative map).

use super::implicit::{decode_poly, encode_poly, valuation_exceptional_indices};
use super::implicit_algebra::{Biv, DOp};
use super::space::{check_kernel_certificate, encode_matrix, kernel_limits, residual_matrix};
use super::*;
use crate::linalg::{rational_kernel_with_certificate, RationalKernelError};
use crate::polynomial::Polynomial;
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
    ) -> Result<Vec<RecurrenceMonomial>, String> {
        match self.role.as_str() {
            SUCCESSOR_SUPPORT | IMPLICIT_SUCCESSOR_SUPPORT if self.max_shift == 0 => {
                return Err(format!("a {} support needs max_shift >= 1", self.role))
            }
            ROW_EQUATION_SUPPORT if self.max_shift != 0 => {
                return Err("a row_equation support has max_shift = 0".into())
            }
            SUCCESSOR_SUPPORT | IMPLICIT_SUCCESSOR_SUPPORT | ROW_EQUATION_SUPPORT => {}
            other => {
                return Err(format!(
                    "unknown support role `{other}` (successor, implicit_successor, row_equation)"
                ))
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
            return Err("support degree budget".into());
        }
        let count = (self.max_shift + 1)
            .checked_mul(self.derivative_order + 1)
            .and_then(|c| c.checked_mul(self.index_degree + 1))
            .and_then(|c| c.checked_mul(self.variable_degree + 1));
        if count.is_none_or(|c| c > limits.max_columns) {
            return Err("support column budget".into());
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

#[derive(Clone, Copy)]
enum Theorem {
    Identity,
    Factorial(usize),
    Derivative(usize),
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
    fn ensure(&mut self, count: usize) -> Result<(), String> {
        if count <= self.rows.len() {
            return Ok(());
        }
        let Some((rule, seeds)) = &self.rule else {
            return Ok(());
        };
        if count > self.request.limits.max_rows {
            return Err("row budget".into());
        }
        let generated = rule
            .generate_rows_rational(seeds, self.origin, count.max(seeds.len()))
            .map_err(|e| format!("source rule generation failed: {e}"))?;
        for row in generated.iter().skip(self.source_rows.len()) {
            let row = vec_poly::trim(row);
            if vec_poly::degree(&row).unwrap_or(0) > self.request.limits.max_degree {
                return Err("row degree budget".into());
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
) -> Result<Vec<Q>, String> {
    let Some(t) = t else {
        return Ok(vec_poly::trim(row));
    };
    let mut state = t
        .map_row(index, row)
        .map_err(|e| format!("{}: {e}", t.kind_name()))?;
    if state.len() != 1 {
        return Err(format!(
            "{} gives {} components; the search needs scalar rows",
            t.kind_name(),
            state.len()
        ));
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
) -> Result<(Context<'a>, ImplicitSearchSource, Vec<String>), String> {
    let limits = &request.limits;
    let mut notes = Vec::new();
    if let Some(t) = &request.transformation {
        t.validate().map_err(|e| format!("transformation: {e}"))?;
        if t.component_count() != 1 {
            return Err(format!(
                "{} gives vector rows; the search needs a scalar map",
                t.kind_name()
            ));
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
            return Err("give at most one of source_recurrence and oeis_id".into())
        }
        (Some(json), None) => {
            let (rule, first, seeds) = json.to_recurrence_parts()?;
            Some(("defining_rule", rule, first, seeds))
        }
        (None, Some(id)) => {
            let definition =
                crate::oeis::by_id(id).ok_or_else(|| format!("unknown OEIS family `{id}`"))?;
            let (rule, seeds, first_row) = definition
                .recurrence_export_parts()
                .map_err(|e| e.to_string())?;
            let first = usize::try_from(first_row)
                .map_err(|_| format!("OEIS family `{id}` starts at a negative row index"))?;
            Some(("oeis", rule, first, seeds))
        }
        (None, None) => None,
    };
    let row_count = request
        .row_count
        .unwrap_or(DEFAULT_IMPLICIT_SEARCH_ROW_COUNT);
    if row_count > limits.max_rows {
        return Err(format!(
            "{row_count} rows exceed max_rows = {}",
            limits.max_rows
        ));
    }
    let (kind, rule, origin, json) = match defining {
        Some((kind, rule, first, seeds)) => {
            if rule.terms.iter().any(|t| t.offset == 0) {
                return Err("source rule terms need positive offsets".into());
            }
            if seeds.is_empty() || rule.max_offset() > seeds.len() {
                return Err("source rule needs at least max-offset seeds".into());
            }
            if request.row_origin.is_some_and(|o| o != first) {
                return Err(format!(
                    "row_origin differs from the source first index {first}"
                ));
            }
            let json = RecurrenceJson::from_recurrence_rational(&rule, first, &seeds, None);
            (kind, Some((rule, seeds)), first, Some(json))
        }
        None => ("rows_only", None, request.row_origin.unwrap_or(0), None),
    };
    let start = request.domain_start.unwrap_or(origin);
    if start < origin {
        return Err("domain_start precedes the row origin".into());
    }
    let order = json.as_ref().and_then(constant_order);
    let theorem = match (&request.transformation, order) {
        (_, None) => None,
        (None, Some(_)) => Some(Theorem::Identity),
        (
            Some(RecurrenceTransformation::Diagonal {
                weight: DiagonalWeight::FactorialPower { power },
            }),
            Some(_),
        ) => Some(Theorem::Factorial(*power)),
        (Some(RecurrenceTransformation::Derivative { order }), Some(_)) => {
            Some(Theorem::Derivative(*order))
        }
        (Some(_), Some(_)) => None,
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
                let first_base = (origin + seeds.len())
                    .checked_sub(shift)
                    .ok_or("the rule starts before its leader shift")?;
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
                return Err("at least one row is required".into());
            }
            if rows.len() > limits.max_rows {
                return Err(format!(
                    "{} rows exceed max_rows = {}",
                    rows.len(),
                    limits.max_rows
                ));
            }
            for (i, row) in rows.iter().enumerate() {
                if vec_poly::degree(row).unwrap_or(0) > limits.max_degree {
                    return Err("row degree exceeds max_degree".into());
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
                return Err("rows disagree with the rows generated by the source rule".into());
            }
        }
        None if context.rule.is_none() => {
            return Err("give source_recurrence, oeis_id or rows".into())
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
) -> Result<Vec<Vec<Q>>, String> {
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
) -> Result<Vec<Vec<Q>>, String> {
    if count
        .checked_mul(cols.len())
        .is_none_or(|e| e > limits.max_matrix_entries)
    {
        return Err("condition matrix entry budget".into());
    }
    (start..start + count)
        .map(|n| {
            cols.iter()
                .map(|c| {
                    let row = rows
                        .get(n + c.shift - origin)
                        .ok_or("missing row for a condition")?;
                    let value = row
                        .get(c.coefficient_index)
                        .cloned()
                        .unwrap_or_else(Q::zero);
                    let power =
                        u32::try_from(c.index_degree).map_err(|_| "index degree overflow")?;
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
) -> Result<Vec<usize>, String> {
    let matrix = leading
        .iter()
        .map(|&i| basis.iter().map(|v| v[i].clone()).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let kernel = rational_kernel_with_certificate(&matrix, basis.len(), kernel_limits(limits))
        .map_err(|_| "projection budget".to_string())?;
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
        Theorem::Identity => us.checked_add(1),
        Theorem::Derivative(m) => us.checked_add(m)?.checked_add(1),
        Theorem::Factorial(m) => {
            RecurrenceCoefficientRescaling::FactorialPower { power: m }.annihilator_power(cols)
        }
    }
}

fn run_support(
    context: &mut Context<'_>,
    support: &ImplicitSearchSupport,
) -> Result<ImplicitSupportResult, String> {
    let limits = &context.request.limits;
    let mut result = empty_result(support);
    let cols = match support.monomials(limits) {
        Ok(cols) => cols,
        Err(error) if error.contains("budget") => {
            result.reason = Some(error);
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
            if error.contains("budget") {
                result.reason = Some(error);
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
                result.reason = Some(error);
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
        Err(_) => return Err("invalid residual matrix".into()),
    };
    let leading = leader_columns(&cols, r);
    let pivots = match projection(&kernel.basis, &leading, limits) {
        Ok(p) => p,
        Err(error) => {
            result.reason = Some(error);
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
) -> Result<ImplicitSearchRepresentative, String> {
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
                return Err(format!(
                    "internal: a certified relation fails at base {base}"
                ));
            }
            return Err(format!("internal: a fitted relation fails at base {base}"));
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
        return Err("internal: missing seed rows".into());
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
) -> Result<ImplicitConditionSearch, String> {
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
    let truncation = match context.theorem {
        Some(Theorem::Derivative(m)) => k + m,
        _ => k,
    };
    let (count, complete) = match (context.theorem, context.order) {
        (Some(_), Some(order)) => match order
            .checked_mul(truncation + 1)
            .and_then(|c| c.checked_mul(degree + 1))
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
            if error.contains("budget") {
                search.reason = Some(error);
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
            search.reason = Some(error);
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
        Err(_) => return Err("invalid condition matrix".into()),
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
                super::implicit_algebra::roots_in_domain(&t.coeff, first_base, None).is_empty()
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
) -> Result<ImplicitRelationSearchReport, String> {
    if request.supports.len() > MAX_IMPLICIT_SEARCH_SUPPORTS {
        return Err(format!(
            "at most {MAX_IMPLICIT_SEARCH_SUPPORTS} supports are supported"
        ));
    }
    if request.condition_index_degree > 8 {
        return Err("condition_index_degree is at most 8".into());
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
) -> Result<(), String> {
    if report.schema != RECURRENCE_IMPLICIT_SEARCH_SCHEMA {
        return Err(format!("unsupported schema `{}`", report.schema));
    }
    let replay = search_implicit_relations(&report.request)?;
    if &replay != report {
        return Err("report does not match its exact replay".into());
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
                return Err(format!("{label}: representatives without a certificate"));
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
            .map_err(|e| format!("{label}: {e}"))?;
        let leading = leader_columns(&result.monomials, result.support.max_shift);
        if projection(&basis, &leading, limits)?.len() != result.leader_projection_rank {
            return Err(format!("{label}: leader projection mismatch"));
        }
        for (i, rep) in result.representatives.iter().enumerate() {
            let label = format!("{label}.representatives[{i}]");
            let vector = basis
                .get(rep.basis_index)
                .ok_or_else(|| format!("{label}: unknown basis index"))?;
            let relation = ImplicitRelation::from_json(&rep.relation.relation)?;
            if !relation.same_relation(&relation_from_vector(&result.monomials, vector)) {
                return Err(format!("{label}: the relation is not its kernel vector"));
            }
            let end = (origin + rows.len()).saturating_sub(relation.max_shift());
            for base in rep.relation.first_base..end {
                if relation.holds_at(&rows, origin, base) == Some(false) {
                    return Err(format!("{label}: the relation fails at base {base}"));
                }
            }
            if let Some(division) = &rep.division_certificate {
                division.verify_identity()?;
                if division.source_relation != relation.to_json() {
                    return Err(format!(
                        "{label}: the division certificate is for another relation"
                    ));
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
                        .map_err(|e| format!("{label}: condition search: {e}"))?;
                if let (Some(i), Some(condition)) = (search.chosen_basis_index, &search.condition) {
                    let vector = basis
                        .get(i)
                        .ok_or_else(|| format!("{label}: unknown condition basis index"))?;
                    let parsed = CoefficientCondition::from_json(&condition.condition)?;
                    if parsed != condition_from_vector(&search.monomials, vector) {
                        return Err(format!("{label}: the condition is not its kernel vector"));
                    }
                    let end = (origin + rows.len()).saturating_sub(parsed.max_shift());
                    for base in condition.first_base..end {
                        if parsed.holds_at(&rows, origin, base) == Some(false) {
                            return Err(format!("{label}: a condition fails at base {base}"));
                        }
                    }
                }
            }
            if assess_generating_description(&rep.description)? != rep.assessment {
                return Err(format!("{label}: the assessment does not replay"));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "implicit_search_tests.rs"]
mod tests;
