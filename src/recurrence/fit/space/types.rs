//! Requests, results and certificates of complete bounded recurrence spaces.

use super::*;
use crate::recurrence::RecurrenceError;

/// Exact polynomial coefficient rows use ascending powers and rational strings.
pub type ExactRows = Vec<Vec<String>>;

/// A sequence defined by `P[n+k] = sum_i coefficients[i](X) P[n+i]`.
/// It holds from `origin`; exactly k seeds must be supplied.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ConstantLagSource {
    /// `coefficients[i]` holds the coefficients of `a_i(X)`, in ascending powers of `X`,
    /// as exact rational strings; the source order `k` is `coefficients.len()`.
    pub coefficients: ExactRows,
    /// Initial rows: `seeds[j]` is `P[origin + j]` as ascending-power rational strings.
    /// Exactly `k` seeds are required.
    pub seeds: ExactRows,
    /// Original index of `seeds[0]`; the defining rule holds for base indices `n >= origin`.
    /// Defaults to zero.
    #[serde(default)]
    pub origin: usize,
    /// Free-form note on where the source came from. It is informational and is never
    /// used in a computation.
    #[serde(default)]
    pub provenance: Option<String>,
}

impl ConstantLagSource {
    /// Canonical exact generator for a supplied monic source, starting at
    /// `origin`. This is the same generator used by the certificate checker.
    pub fn generate_rows(
        &self,
        count: usize,
        limits: &RecurrenceSpaceLimits,
    ) -> Result<ExactRows, RecurrenceError> {
        let request = RecurrenceSpaceRequest {
            source: Some(self.clone()),
            source_recurrence: None,
            indexed_source: None,
            coefficient_rescaling: None,
            include_leader_coverage: false,
            rows: None,
            row_origin: self.origin,
            domain_start: self.origin,
            shapes: Vec::new(),
            limits: limits.clone(),
        };
        let input = prepare(&request)?;
        let mut rows = source_rows(&input, count, limits)?;
        rows.truncate(count);
        Ok(encode_polys(&rows))
    }
    /// Adapt the existing backward recurrence JSON. Extra initial rows must
    /// be checked separately by the caller; they are returned as consistency rows.
    /// The original index is retained, including one-based input origins.
    pub fn from_recurrence_json(
        json: &RecurrenceJson,
    ) -> Result<(Self, ExactRows), RecurrenceError> {
        let (rec, origin, initial) = json.to_recurrence_parts()?;
        if rec.inhomogeneous.is_some() {
            return Err(RecurrenceError::InvalidInput(
                "source must be homogeneous".into(),
            ));
        }
        if let Some(d) = &rec.denominator {
            if d.coeffs.iter().enumerate().any(|(n, row)| {
                row.iter().enumerate().any(|(x, c)| {
                    *c != if n == 0 && x == 0 {
                        Q::one()
                    } else {
                        Q::zero()
                    }
                })
            }) || d.coeffs.first().and_then(|r| r.first()) != Some(&Q::one())
            {
                return Err(RecurrenceError::InvalidInput(
                    "source must be monic (denominator one)".into(),
                ));
            }
        }
        let k = rec.max_offset();
        if k == 0 || initial.len() < k {
            return Err(RecurrenceError::InvalidInput(
                "source needs positive order and sufficient seeds".into(),
            ));
        }
        let mut a = vec![Polynomial::<Q>::zero(); k];
        for term in rec.terms {
            if term.offset == 0
                || term.deriv_order != 0
                || term.sign != RecurrenceSign::None
                || term
                    .coeff
                    .coeffs
                    .iter()
                    .skip(1)
                    .flatten()
                    .any(|x| !x.is_zero())
            {
                return Err(RecurrenceError::InvalidInput(
                    "source must be pure-lag and independent of the index".into(),
                ));
            }
            a[k - term.offset] = a[k - term.offset].clone()
                + Polynomial::new(term.coeff.coeffs.first().cloned().unwrap_or_default());
        }
        let rows: ExactRows = initial
            .iter()
            .map(|p| p.iter().map(ToString::to_string).collect())
            .collect();
        Ok((
            Self {
                coefficients: encode_polys(&a),
                seeds: rows[..k].to_vec(),
                origin,
                provenance: Some("supplied polytool.recurrence.v1; not inferred truth".into()),
            },
            rows,
        ))
    }
}

/// One left-coefficient monomial. At `lag == target_lag` this is a positive
/// successor coefficient; all other columns enter the residual negatively.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceMonomial {
    /// Lag `j` of the referenced row `P[n + j]`; the leader sits at the shape's target lag.
    pub lag: usize,
    /// Derivative order `s` in `X` applied to the referenced row. Defaults to zero.
    #[serde(default)]
    pub derivative_order: usize,
    /// Power `u` of the base index `n` in the monomial `n^u X^v D^s P[n + lag]`.
    /// Defaults to zero.
    #[serde(default)]
    pub index_degree: usize,
    /// Power `v` of the polynomial variable `X` in the monomial. Defaults to zero.
    #[serde(default)]
    pub variable_degree: usize,
}

/// Either a rectangular degree box or an exact, possibly sparse support.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum RecurrenceSpaceShape {
    /// Every right-hand-side monomial `n^u X^v D^s P[n + j]` with `j < length`, `s <= derivative_order`,
    /// `u <= index_degree` and `v <= variable_degree`, together with leader monomials
    /// `n^u X^v P[n + length]` (no derivative) bounded by `leading_index_degree` and
    /// `leading_variable_degree`.
    Box {
        /// Target lag: the leader multiplies `P[n + length]`; right-hand-side lags are `0..length`.
        length: usize,
        /// Maximum derivative order of right-hand-side terms.
        #[serde(default)]
        derivative_order: usize,
        /// Maximum power of `n` in right-hand-side coefficients.
        #[serde(default)]
        index_degree: usize,
        /// Maximum power of `X` in right-hand-side coefficients.
        #[serde(default)]
        variable_degree: usize,
        /// Maximum power of `n` in the leader coefficient.
        #[serde(default)]
        leading_index_degree: usize,
        /// Maximum power of `X` in the leader coefficient.
        #[serde(default)]
        leading_variable_degree: usize,
    },
    /// Remove the RHS X-degree cap using the polynomial-matrix solution bound.
    AutomaticVariableDegree {
        /// Target lag: the leader multiplies `P[n + length]`; right-hand-side lags are `0..length`.
        length: usize,
        /// Maximum derivative order of right-hand-side terms.
        #[serde(default)]
        derivative_order: usize,
        /// Maximum power of `n` in right-hand-side coefficients.
        #[serde(default)]
        index_degree: usize,
        /// Maximum power of `n` in the leader coefficient.
        #[serde(default)]
        leading_index_degree: usize,
        /// Maximum power of `X` in the leader coefficient.
        #[serde(default)]
        leading_variable_degree: usize,
    },
    /// An explicit support: the listed monomials are exactly the columns of the residual
    /// matrix. Columns at `target_lag` form the leader.
    Explicit {
        /// Lag of the leader columns, which must be positive. Every monomial lag is at most
        /// `target_lag`, and monomials at `target_lag` must have `derivative_order == 0`.
        target_lag: usize,
        /// Distinct support monomials; they are sorted before use.
        monomials: Vec<RecurrenceMonomial>,
    },
}

/// Bounds apply to each generation/elimination phase of each shape; reaching
/// one never proves nonexistence. Services also cap the number of shapes.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default, deny_unknown_fields)]
pub struct RecurrenceSpaceLimits {
    /// Maximum number of requested shapes.
    pub max_shapes: usize,
    /// Maximum number of sequence rows that are generated or supplied.
    pub max_rows: usize,
    /// Maximum number of support monomials (columns) in one shape.
    pub max_columns: usize,
    /// Maximum degree, in both `n` and `X`, of any polynomial that is parsed or computed.
    pub max_degree: usize,
    /// Maximum number of entries of any matrix or coefficient array.
    pub max_matrix_entries: usize,
    /// Budget of elementary arithmetic work units for each generation or elimination phase.
    pub max_operations: usize,
    /// Maximum bit length of a numerator or denominator of any rational coefficient.
    pub max_coefficient_bits: u64,
}

impl Default for RecurrenceSpaceLimits {
    fn default() -> Self {
        Self {
            max_shapes: 16,
            max_rows: 512,
            max_columns: 128,
            max_degree: 512,
            max_matrix_entries: 500_000,
            max_operations: 20_000_000,
            max_coefficient_bits: 4096,
        }
    }
}

impl RecurrenceSpaceLimits {
    /// Reject requests that raise a service's configured resource ceilings.
    pub fn check_ceiling(&self, ceiling: &Self) -> Result<(), RecurrenceError> {
        if self.max_shapes > ceiling.max_shapes
            || self.max_rows > ceiling.max_rows
            || self.max_columns > ceiling.max_columns
            || self.max_degree > ceiling.max_degree
            || self.max_matrix_entries > ceiling.max_matrix_entries
            || self.max_operations > ceiling.max_operations
            || self.max_coefficient_bits > ceiling.max_coefficient_bits
        {
            return Err(RecurrenceError::LimitExceeded(
                "requested limits exceed recurrence-space service ceilings".into(),
            ));
        }
        Ok(())
    }
}

/// Request for complete coefficient spaces. `source` and `source_recurrence`
/// are alternative forms; `rows` may additionally check source consistency.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceSpaceRequest {
    /// Constant-index monic source. At most one of `source`, `source_recurrence` and
    /// `indexed_source` may be given.
    #[serde(default)]
    pub source: Option<ConstantLagSource>,
    /// An existing `polytool.recurrence.v1` object that is converted to a constant-index
    /// source; it must be homogeneous, monic, lag-only and independent of `n`.
    #[serde(default)]
    pub source_recurrence: Option<RecurrenceJson>,
    /// Monic source whose coefficients are polynomials in the index `n` as well as `X`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub indexed_source: Option<IndexedLagSource>,
    /// Opt-in exact coefficient map applied to every source row before the search:
    /// the shapes then describe identities of `Q[n] = T(P[n])`. Requires `source` or
    /// `source_recurrence`; see the `fit::rescaling` module notes for the finite theorem.
    /// Omitted from the output when absent, which keeps the untransformed behavior.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coefficient_rescaling: Option<RecurrenceCoefficientRescaling>,
    /// When true, also compute the sum of every representative multiplied by its own
    /// leader polynomial, in which singular ranks of different representatives can
    /// cancel. Omitted from the output when false.
    #[serde(default, skip_serializing_if = "is_false")]
    pub include_leader_coverage: bool,
    /// With a source: optional rows that must agree with the source's generated rows.
    /// Without a source: the finite data to analyze. `rows[i]` is `P[row_origin + i]`,
    /// written as ascending-power rational strings.
    #[serde(default)]
    pub rows: Option<ExactRows>,
    /// Original index of `rows[0]`. With a source it must equal the source origin when
    /// `rows` is given. Defaults to zero.
    #[serde(default)]
    pub row_origin: usize,
    /// First base index `n` at which residuals are tested; it must not precede the source
    /// origin (or `row_origin` without a source). Defaults to zero.
    #[serde(default)]
    pub domain_start: usize,
    /// Supports to analyze; the report's `results` has one entry per shape, in this order.
    pub shapes: Vec<RecurrenceSpaceShape>,
    /// Resource limits; omitted fields take the defaults.
    #[serde(default)]
    pub limits: RecurrenceSpaceLimits,
}

/// A coordinate row-space witness, supported on the common independent rows.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceSpaceDual {
    /// Index into the result's `monomials` of the leader column that this witness isolates.
    pub coordinate: usize,
    /// Rational weights, one per entry of `independent_rows` (same order). The weighted sum
    /// of those residual-matrix rows is the unit vector at `coordinate`.
    pub weights: Vec<String>,
}

/// Exact nullspace and rank witnesses. Each `basis` entry is a column vector.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceSpaceCertificate {
    /// Kernel basis of the residual matrix; `basis[i]` is a vector with one entry per
    /// monomial (same order as `monomials`). Entry `free_columns[j]` of `basis[i]` is 1 if
    /// `i == j` and 0 otherwise.
    pub basis: ExactRows,
    /// Column indices (into `monomials`) of the free variables, one per basis vector.
    pub free_columns: Vec<usize>,
    /// Column indices (into `monomials`) of the pivot variables; their count is the rank of
    /// the residual matrix.
    pub pivot_columns: Vec<usize>,
    /// Zero-based row indices of the residual matrix forming a maximal independent set,
    /// one per pivot column.
    pub independent_rows: Vec<usize>,
    /// Inverse of the square submatrix of the residual matrix with rows `independent_rows`
    /// and columns `pivot_columns` (rows and columns in listed order).
    pub inverse_minor: ExactRows,
    /// For `absent_in_shape`, one dual witness per leader column; empty otherwise.
    pub absence_duals: Vec<RecurrenceSpaceDual>,
}

/// One useful basis representative, without cancelling boundary factors.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceSpaceRepresentative {
    /// Index of the kernel vector in the certificate's `basis`.
    pub basis_index: usize,
    /// Exact rational coefficients, one per entry of the result's `monomials` (same order).
    pub coefficients: Vec<String>,
    /// Original forward-index formula (X is the polynomial variable).
    pub formula: String,
    /// All nonnegative integer ranks where the leader is identically zero in X.
    /// `None` means the bounded analysis was not performed; never means no zeros.
    pub singular_ranks: Option<Vec<String>>,
    /// True when a leader monomial involves a positive power of `X`, so generating rows
    /// from this identity would divide by a polynomial in `X`.
    pub requires_x_division: bool,
}

impl RecurrenceSpaceRepresentative {
    /// Convert to the existing evaluator's backward convention. Its target
    /// index m equals the report's original base index n plus `target_lag`;
    /// all coefficient polynomials are substituted accordingly. Supply seeds
    /// through the declared domain and every singular target before generation.
    /// This conversion does not itself verify an untrusted report.
    pub fn to_recurrence(
        &self,
        result: &RecurrenceSpaceResult,
    ) -> Result<Recurrence, RecurrenceError> {
        if self.coefficients.len() != result.monomials.len() {
            return Err(RecurrenceError::InvalidInput(
                "coefficient dimension mismatch".into(),
            ));
        }
        let limits = RecurrenceSpaceLimits::default();
        let shift = i64::try_from(result.target_lag)
            .map_err(|_| RecurrenceError::LimitExceeded("target index overflow".into()))?;
        if result.monomials.iter().any(|c| {
            c.index_degree > limits.max_degree
                || c.variable_degree > limits.max_degree
                || c.lag > result.target_lag
                || (c.lag == result.target_lag && c.derivative_order != 0)
        }) {
            return Err(RecurrenceError::InvalidInput(
                "invalid recurrence support".into(),
            ));
        }
        let coefficients: Vec<_> = self
            .coefficients
            .iter()
            .map(|s| checked_q(s, &limits))
            .collect::<Result<_, _>>()?;
        let mut denominator = None;
        let mut terms = Vec::new();
        for ((lag, d), p) in grouped_coefficients(&result.monomials, &coefficients) {
            if lag == result.target_lag {
                denominator = Some(p.shift_index(-shift));
            } else {
                terms.push(RecurrenceTerm::new(
                    result.target_lag - lag,
                    d,
                    p.shift_index(-shift),
                ));
            }
        }
        if denominator.is_none() {
            return Err(RecurrenceError::InvalidInput(
                "zero successor coefficient".into(),
            ));
        }
        Ok(Recurrence {
            terms,
            denominator,
            inhomogeneous: None,
        })
    }
}

/// A safe bound preserving every leader under unbounded RHS X-degrees.
/// The kernel itself is still finite; index/leader bounds remain prescribed.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceAutomaticDegreeBound {
    /// Largest `X`-degree of any right-hand-side column evaluated on the checked rows.
    pub rhs_degree: usize,
    /// Largest `X`-degree of any leader column evaluated on the checked rows.
    pub leader_degree: usize,
    /// Bound on the rank of the right-hand-side part of the residual system: the smaller
    /// of the number of tested base indices and the number of right-hand-side columns, or
    /// zero when every right-hand-side column vanishes on the checked rows.
    pub rank_upper_bound: usize,
    /// Resulting bound on the `X`-degree of right-hand-side coefficients that is used for the shape.
    pub variable_degree: usize,
}

/// Sum-of-squares leader combination, using each representative's leader
/// as a left polynomial multiplier. Degree bounds can increase; lag and
/// derivative order do not. Singularities are unavoidable only within the
/// polynomial module generated by the supplied complete kernel.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceLeaderCoverage {
    /// `basis_index` of each representative that was combined.
    pub basis_indices: Vec<usize>,
    /// For each combined representative, its leader coefficient as an array indexed
    /// `[n power][X power]`, used as a left multiplier.
    pub multipliers: Vec<ExactRows>,
    /// Support of the combined identity, in sorted order.
    pub monomials: Vec<RecurrenceMonomial>,
    /// Exact rational coefficients, one per entry of `monomials`.
    pub coefficients: Vec<String>,
    /// Human-readable form of the combined identity, with `X` as the polynomial variable.
    pub formula: String,
    /// Nonnegative integer base ranks at which every combined leader vanishes identically
    /// in `X`, as decimal strings. `None` means the bounded analysis was not performed.
    pub common_singular_ranks: Option<Vec<String>>,
    /// True when the combined leader involves a positive power of `X`.
    pub requires_x_division: bool,
}
impl RecurrenceLeaderCoverage {
    /// Convert the combined formula to the existing backward row evaluator.
    pub fn to_recurrence(&self, target_lag: usize) -> Result<Recurrence, RecurrenceError> {
        let result = empty_result(target_lag, self.monomials.clone(), 0);
        RecurrenceSpaceRepresentative {
            basis_index: 0,
            coefficients: self.coefficients.clone(),
            formula: self.formula.clone(),
            singular_ranks: self.common_singular_ranks.clone(),
            requires_x_division: self.requires_x_division,
        }
        .to_recurrence(&result)
    }
}

/// A completed shape or an explicit resource/row shortage.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceSpaceResult {
    /// One of `successor_exists`, `absent_in_shape`, `insufficient_rows` or
    /// `incomplete_budget`.
    pub status: String,
    /// One of `conditional_on_source`, `exact_fit`, `rows_only_obstruction` or `none`.
    /// `none` accompanies `insufficient_rows` and `incomplete_budget` and never means absence.
    pub evidence: String,
    /// Explanation, present when the shape is incomplete.
    pub reason: Option<String>,
    /// Lag of the leader columns of this shape.
    pub target_lag: usize,
    /// Columns of the residual matrix in sorted order; empty if the shape could not be
    /// resolved within the limits.
    pub monomials: Vec<RecurrenceMonomial>,
    /// One plus the largest `index_degree + derivative_order` among the columns; this is
    /// the exponent `h` in the finite certificate bound of the module documentation.
    /// Under `coefficient_rescaling` it is the rescaled exponent `U + m(A+J) + 1`.
    pub h: usize,
    /// Number of consecutive base indices, starting at `domain_start`, at which the
    /// residual is tested.
    pub residual_count: usize,
    /// Original index of the last sequence row used; `None` if the shape could not be
    /// prepared.
    pub last_required_source_row: Option<usize>,
    /// Number of scalar equations (rows of the residual matrix): the `X`-coefficients of
    /// all tested residuals.
    pub equations: usize,
    /// Rank of the kernel basis projected onto the leader columns; zero means that no
    /// identity in the shape has a nonzero leader.
    pub successor_projection_rank: usize,
    /// Exact kernel and rank witnesses; absent for incomplete shapes.
    pub certificate: Option<RecurrenceSpaceCertificate>,
    /// Basis vectors whose projections onto the leader columns are independent.
    pub representatives: Vec<RecurrenceSpaceRepresentative>,
    /// True when every sequence row used is the zero polynomial, in which case every
    /// identity holds trivially.
    pub all_checked_rows_zero: bool,
    /// Degree bound used for an `automatic_variable_degree` shape.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub automatic_degree_bound: Option<RecurrenceAutomaticDegreeBound>,
    /// Certificate for a source with index-dependent coefficients.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub propagation: Option<RecurrencePropagationCertificate>,
    /// Combination of the representatives; present only when requested.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub leader_coverage: Option<RecurrenceLeaderCoverage>,
    /// Finite-theorem data for a request with `coefficient_rescaling`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rescaling_bound: Option<RecurrenceRescalingBound>,
}

/// Versioned report, retaining the assumptions needed to replay certificates.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceSpaces {
    /// Schema identifier, `polytool.recurrence-spaces.v1`.
    pub schema: String,
    /// The request that produced the report; verification checks that it matches the expected request.
    pub request: RecurrenceSpaceRequest,
    /// One result per requested shape, in request order.
    pub results: Vec<RecurrenceSpaceResult>,
    /// Indices of nondominated successful requested (lag, derivative) supports.
    /// Only compares identical coefficient bounds and domain via box shapes.
    pub frontier: Vec<usize>,
}
