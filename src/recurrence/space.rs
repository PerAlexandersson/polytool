//! Complete bounded recurrence spaces and replayable exact certificates.
//!
//! Coefficients are evaluated in the original index, never silently shifted.
//! A supplied monic constant-index source gives all-index identities via the
//! `h*k` residual bound. Without a source, positives remain finite fits, but
//! dual witnesses can still prove absence in the requested support.
//!
//! # Finite certificate bound
//!
//! Write `chi = S^k - sum_i a_i(X) S^i` for the monic source operator.
//! Since its coefficients are independent of n, `chi` commutes with shifts,
//! multiplication by X, and its own coefficient derivatives. Differentiating
//! `chi P = 0` gives `chi D(P) = -D(chi) P`; induction gives
//! `chi^(s+1) D^s(P) = 0`. The commutator with multiplication by n is
//! `S (d chi / d S)`, which commutes with chi. Thus multiplying by n raises
//! the required power by at most one. Every column `X^v n^u S^j D^s(P)` is
//! annihilated by `chi^(u+s+1)`; their sum is annihilated by `chi^h`.
//! This is a monic recurrence of order h*k. Its first h*k zero values on
//! the declared domain therefore imply all subsequent values are zero.
//! No factorization, distinct-root assumption, or invertible oldest
//! coefficient is required. A source inferred from rows remains an assumption.

use super::{
    parse_rational_coeff, BigRational as Q, BivarPoly, Recurrence, RecurrenceJson, RecurrenceSign,
    RecurrenceTerm,
};
use crate::linalg::{rational_kernel_with_certificate, RationalKernelError, RationalKernelLimits};
use crate::{sturm::SturmChain, Polynomial};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[path = "propagation.rs"]
mod propagation;
pub use propagation::{IndexedLagSource, RecurrencePropagationCertificate};

/// Exact polynomial coefficient rows use ascending powers and rational strings.
pub type ExactRows = Vec<Vec<String>>;

/// A sequence defined by `P[n+k] = sum_i coefficients[i](X) P[n+i]`.
/// It holds from `origin`; exactly k seeds must be supplied.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ConstantLagSource {
    pub coefficients: ExactRows,
    pub seeds: ExactRows,
    #[serde(default)]
    pub origin: usize,
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
    ) -> Result<ExactRows, String> {
        let request = RecurrenceSpaceRequest {
            source: Some(self.clone()),
            source_recurrence: None,
            indexed_source: None,
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
    pub fn from_recurrence_json(json: &RecurrenceJson) -> Result<(Self, ExactRows), String> {
        let (rec, origin, initial) = json.to_recurrence_parts().map_err(|e| e.to_string())?;
        if rec.inhomogeneous.is_some() {
            return Err("source must be homogeneous".into());
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
                return Err("source must be monic (denominator one)".into());
            }
        }
        let k = rec.max_offset();
        if k == 0 || initial.len() < k {
            return Err("source needs positive order and sufficient seeds".into());
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
                return Err("source must be pure-lag and independent of the index".into());
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
    pub lag: usize,
    #[serde(default)]
    pub derivative_order: usize,
    #[serde(default)]
    pub index_degree: usize,
    #[serde(default)]
    pub variable_degree: usize,
}

/// Either a rectangular degree box or an exact, possibly sparse support.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum RecurrenceSpaceShape {
    Box {
        length: usize,
        #[serde(default)]
        derivative_order: usize,
        #[serde(default)]
        index_degree: usize,
        #[serde(default)]
        variable_degree: usize,
        #[serde(default)]
        leading_index_degree: usize,
        #[serde(default)]
        leading_variable_degree: usize,
    },
    /// Remove the RHS X-degree cap using the polynomial-matrix solution bound.
    AutomaticVariableDegree {
        length: usize,
        #[serde(default)]
        derivative_order: usize,
        #[serde(default)]
        index_degree: usize,
        #[serde(default)]
        leading_index_degree: usize,
        #[serde(default)]
        leading_variable_degree: usize,
    },
    Explicit {
        target_lag: usize,
        monomials: Vec<RecurrenceMonomial>,
    },
}

/// Bounds apply to each generation/elimination phase of each shape; reaching
/// one never proves nonexistence. Services also cap the number of shapes.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default, deny_unknown_fields)]
pub struct RecurrenceSpaceLimits {
    pub max_shapes: usize,
    pub max_rows: usize,
    pub max_columns: usize,
    pub max_degree: usize,
    pub max_matrix_entries: usize,
    pub max_operations: usize,
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
    pub fn check_ceiling(&self, ceiling: &Self) -> Result<(), String> {
        if self.max_shapes > ceiling.max_shapes
            || self.max_rows > ceiling.max_rows
            || self.max_columns > ceiling.max_columns
            || self.max_degree > ceiling.max_degree
            || self.max_matrix_entries > ceiling.max_matrix_entries
            || self.max_operations > ceiling.max_operations
            || self.max_coefficient_bits > ceiling.max_coefficient_bits
        {
            return Err("requested limits exceed recurrence-space service ceilings".into());
        }
        Ok(())
    }
}

/// Request for complete coefficient spaces. `source` and `source_recurrence`
/// are alternative forms; `rows` may additionally check source consistency.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceSpaceRequest {
    #[serde(default)]
    pub source: Option<ConstantLagSource>,
    #[serde(default)]
    pub source_recurrence: Option<RecurrenceJson>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub indexed_source: Option<IndexedLagSource>,
    #[serde(default, skip_serializing_if = "is_false")]
    pub include_leader_coverage: bool,
    #[serde(default)]
    pub rows: Option<ExactRows>,
    #[serde(default)]
    pub row_origin: usize,
    #[serde(default)]
    pub domain_start: usize,
    pub shapes: Vec<RecurrenceSpaceShape>,
    #[serde(default)]
    pub limits: RecurrenceSpaceLimits,
}

/// A coordinate row-space witness, supported on the common independent rows.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceSpaceDual {
    pub coordinate: usize,
    pub weights: Vec<String>,
}

/// Exact nullspace and rank witnesses. Each `basis` entry is a column vector.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceSpaceCertificate {
    pub basis: ExactRows,
    pub free_columns: Vec<usize>,
    pub pivot_columns: Vec<usize>,
    pub independent_rows: Vec<usize>,
    pub inverse_minor: ExactRows,
    pub absence_duals: Vec<RecurrenceSpaceDual>,
}

/// One useful basis representative, without cancelling boundary factors.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceSpaceRepresentative {
    pub basis_index: usize,
    pub coefficients: Vec<String>,
    /// Original forward-index formula (X is the polynomial variable).
    pub formula: String,
    /// All nonnegative integer ranks where the leader is identically zero in X.
    /// `None` means the bounded analysis was not performed; never means no zeros.
    pub singular_ranks: Option<Vec<String>>,
    pub requires_x_division: bool,
}

impl RecurrenceSpaceRepresentative {
    /// Convert to the existing evaluator's backward convention. Its target
    /// index m equals the report's original base index n plus `target_lag`;
    /// all coefficient polynomials are substituted accordingly. Supply seeds
    /// through the declared domain and every singular target before generation.
    /// This conversion does not itself verify an untrusted report.
    pub fn to_recurrence(&self, result: &RecurrenceSpaceResult) -> Result<Recurrence, String> {
        if self.coefficients.len() != result.monomials.len() {
            return Err("coefficient dimension mismatch".into());
        }
        let limits = RecurrenceSpaceLimits::default();
        let shift = i64::try_from(result.target_lag).map_err(|_| "target index overflow")?;
        if result.monomials.iter().any(|c| {
            c.index_degree > limits.max_degree
                || c.variable_degree > limits.max_degree
                || c.lag > result.target_lag
                || (c.lag == result.target_lag && c.derivative_order != 0)
        }) {
            return Err("invalid recurrence support".into());
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
            return Err("zero successor coefficient".into());
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
    pub rhs_degree: usize,
    pub leader_degree: usize,
    pub rank_upper_bound: usize,
    pub variable_degree: usize,
}

/// Sum-of-squares leader combination, using each representative's leader
/// as a left polynomial multiplier. Degree bounds can increase; lag and
/// derivative order do not. Singularities are unavoidable only within the
/// polynomial module generated by the supplied complete kernel.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceLeaderCoverage {
    pub basis_indices: Vec<usize>,
    pub multipliers: Vec<ExactRows>,
    pub monomials: Vec<RecurrenceMonomial>,
    pub coefficients: Vec<String>,
    pub formula: String,
    pub common_singular_ranks: Option<Vec<String>>,
    pub requires_x_division: bool,
}
impl RecurrenceLeaderCoverage {
    /// Convert the combined formula to the existing backward row evaluator.
    pub fn to_recurrence(&self, target_lag: usize) -> Result<Recurrence, String> {
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
    pub status: String,
    pub evidence: String,
    pub reason: Option<String>,
    pub target_lag: usize,
    pub monomials: Vec<RecurrenceMonomial>,
    pub h: usize,
    pub residual_count: usize,
    pub last_required_source_row: Option<usize>,
    pub equations: usize,
    pub successor_projection_rank: usize,
    pub certificate: Option<RecurrenceSpaceCertificate>,
    pub representatives: Vec<RecurrenceSpaceRepresentative>,
    pub all_checked_rows_zero: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub automatic_degree_bound: Option<RecurrenceAutomaticDegreeBound>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub propagation: Option<RecurrencePropagationCertificate>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub leader_coverage: Option<RecurrenceLeaderCoverage>,
}

/// Versioned report, retaining the assumptions needed to replay certificates.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceSpaces {
    pub schema: String,
    pub request: RecurrenceSpaceRequest,
    pub results: Vec<RecurrenceSpaceResult>,
    /// Indices of nondominated successful requested (lag, derivative) supports.
    /// Only compares identical coefficient bounds and domain via box shapes.
    pub frontier: Vec<usize>,
}

type Poly = Polynomial<Q>;
fn is_false(value: &bool) -> bool {
    !*value
}

fn encode_polys(polys: &[Poly]) -> ExactRows {
    polys
        .iter()
        .map(|p| p.coeffs().iter().map(ToString::to_string).collect())
        .collect()
}
fn encode_matrix(matrix: &[Vec<Q>]) -> ExactRows {
    matrix
        .iter()
        .map(|r| r.iter().map(ToString::to_string).collect())
        .collect()
}
fn checked_q(s: &str, limits: &RecurrenceSpaceLimits) -> Result<Q, String> {
    if s.len() as u64 > limits.max_coefficient_bits.saturating_add(32) {
        return Err("coefficient text exceeds bit budget".into());
    }
    let q = parse_rational_coeff(s)?;
    if q.numer().bits().max(q.denom().bits()) > limits.max_coefficient_bits {
        return Err("coefficient exceeds bit budget".into());
    }
    Ok(q)
}
fn parse_matrix(rows: &ExactRows, limits: &RecurrenceSpaceLimits) -> Result<Vec<Vec<Q>>, String> {
    if rows
        .len()
        .checked_mul(rows.iter().map(Vec::len).max().unwrap_or(0))
        .is_none_or(|n| n > limits.max_matrix_entries)
    {
        return Err("matrix exceeds entry budget".into());
    }
    rows.iter()
        .map(|r| r.iter().map(|x| checked_q(x, limits)).collect())
        .collect()
}
fn parse_polys(rows: &ExactRows, limits: &RecurrenceSpaceLimits) -> Result<Vec<Poly>, String> {
    if rows.len() > limits.max_rows
        || rows
            .iter()
            .any(|r| r.len().saturating_sub(1) > limits.max_degree)
    {
        return Err("input rows exceed row/degree budget".into());
    }
    Ok(parse_matrix(rows, limits)?
        .into_iter()
        .map(Poly::new)
        .collect())
}
fn poly_budget(p: &Poly, limits: &RecurrenceSpaceLimits) -> Result<(), String> {
    if p.degree().unwrap_or(0) > limits.max_degree {
        return Err("polynomial degree budget".into());
    }
    if p.coeffs()
        .iter()
        .any(|x| x.numer().bits().max(x.denom().bits()) > limits.max_coefficient_bits)
    {
        return Err("coefficient bit budget".into());
    }
    Ok(())
}

struct Input {
    a: Option<Vec<Poly>>,
    indexed_a: Option<Vec<BivarPoly>>,
    rows: Vec<Poly>,
    origin: usize,
    consistency: Vec<Poly>,
}

fn prepare(request: &RecurrenceSpaceRequest) -> Result<Input, String> {
    let lim = &request.limits;
    if request.shapes.len() > lim.max_shapes {
        return Err("too many requested shapes".into());
    }
    if usize::from(request.source.is_some())
        + usize::from(request.source_recurrence.is_some())
        + usize::from(request.indexed_source.is_some())
        > 1
    {
        return Err("provide exactly one source form".into());
    }
    if let Some(src) = &request.indexed_source {
        if src.coefficients.is_empty()
            || src.coefficients.len() > lim.max_rows
            || src.seeds.len() != src.coefficients.len()
        {
            return Err("indexed source needs exactly k seeds for positive order k".into());
        }
        if request.domain_start < src.origin {
            return Err("domain starts before source origin".into());
        }
        if request.rows.is_some() && request.row_origin != src.origin {
            return Err("consistency row_origin must equal source origin".into());
        }
        return Ok(Input {
            a: None,
            indexed_a: Some(
                src.coefficients
                    .iter()
                    .map(|p| propagation::parse(p, lim))
                    .collect::<Result<_, _>>()?,
            ),
            rows: parse_polys(&src.seeds, lim)?,
            origin: src.origin,
            consistency: request
                .rows
                .as_ref()
                .map(|p| parse_polys(p, lim))
                .transpose()?
                .unwrap_or_default(),
        });
    }
    if let Some(json) = &request.source_recurrence {
        parse_polys(&json.initial_polynomials, lim)?;
        if json.recurrence.terms.len() > lim.max_matrix_entries
            || json
                .recurrence
                .terms
                .iter()
                .any(|t| t.offset > lim.max_rows)
        {
            return Err("source recurrence exceeds limits".into());
        }
        for p in json
            .recurrence
            .terms
            .iter()
            .map(|t| &t.coeff)
            .chain(json.recurrence.denominator.iter())
            .chain(json.recurrence.inhomogeneous.iter())
        {
            if p.coeffs.len().saturating_sub(1) > lim.max_degree
                || p.coeffs
                    .iter()
                    .any(|r| r.len().saturating_sub(1) > lim.max_degree)
            {
                return Err("source coefficient degree budget".into());
            }
            parse_matrix(&p.coeffs, lim)?;
        }
    }
    let adapted = request
        .source_recurrence
        .as_ref()
        .map(ConstantLagSource::from_recurrence_json)
        .transpose()?;
    let source = request
        .source
        .as_ref()
        .or_else(|| adapted.as_ref().map(|x| &x.0));
    let supplied = request
        .rows
        .as_ref()
        .or_else(|| adapted.as_ref().map(|x| &x.1));
    if let Some(src) = source {
        let a = parse_polys(&src.coefficients, lim)?;
        let rows = parse_polys(&src.seeds, lim)?;
        if a.is_empty() || rows.len() != a.len() {
            return Err("source requires exactly k seeds for positive order k".into());
        }
        if request.domain_start < src.origin {
            return Err("domain starts before the source origin".into());
        }
        if request.rows.is_some() && request.row_origin != src.origin {
            return Err("consistency row_origin must equal source origin".into());
        }
        let mut consistency = supplied
            .map(|r| parse_polys(r, lim))
            .transpose()?
            .unwrap_or_default();
        if let Some((_, initial)) = &adapted {
            let original = parse_polys(initial, lim)?;
            if !original.iter().zip(&consistency).all(|(a, b)| a == b) {
                return Err("supplied rows disagree with source_recurrence initial rows".into());
            }
            if original.len() > consistency.len() {
                consistency = original;
            }
        }
        Ok(Input {
            a: Some(a),
            indexed_a: None,
            rows,
            origin: src.origin,
            consistency,
        })
    } else {
        let rows = supplied.ok_or("provide a source recurrence or finite rows")?;
        if request.domain_start < request.row_origin {
            return Err("domain starts before row_origin".into());
        }
        Ok(Input {
            a: None,
            indexed_a: None,
            rows: parse_polys(rows, lim)?,
            origin: request.row_origin,
            consistency: Vec::new(),
        })
    }
}

fn source_rows(
    input: &Input,
    count: usize,
    limits: &RecurrenceSpaceLimits,
) -> Result<Vec<Poly>, String> {
    let needed = count.max(input.consistency.len());
    if needed > limits.max_rows {
        return Err("source row budget".into());
    }
    let mut rows = input.rows.clone();
    let mut operations = limits.max_operations;
    if input.a.is_some() || input.indexed_a.is_some() {
        while rows.len() < needed {
            let order = input
                .a
                .as_ref()
                .map(Vec::len)
                .unwrap_or_else(|| input.indexed_a.as_ref().unwrap().len());
            let start = rows.len() - order;
            let evaluated;
            let a = if let Some(a) = &input.a {
                a
            } else {
                let rank = input
                    .origin
                    .checked_add(start)
                    .ok_or("source rank overflow")?;
                let source = input.indexed_a.as_ref().unwrap();
                let work = source
                    .iter()
                    .flat_map(|p| &p.coeffs)
                    .map(Vec::len)
                    .sum::<usize>();
                operations = operations
                    .checked_sub(work)
                    .ok_or("source arithmetic budget")?;
                evaluated = source
                    .iter()
                    .map(|p| propagation::evaluate(p, rank, limits))
                    .collect::<Result<Vec<_>, _>>()?;
                &evaluated
            };
            let mut p = Poly::zero();
            for (i, c) in a.iter().enumerate() {
                let work = c
                    .coeffs()
                    .len()
                    .checked_mul(rows[start + i].coeffs().len())
                    .ok_or("source work overflow")?;
                operations = operations
                    .checked_sub(work)
                    .ok_or("source arithmetic budget")?;
                if c.degree()
                    .unwrap_or(0)
                    .saturating_add(rows[start + i].degree().unwrap_or(0))
                    > limits.max_degree
                {
                    return Err("source degree budget".into());
                }
                p = p + c.clone() * rows[start + i].clone();
                poly_budget(&p, limits)?;
            }
            rows.push(p);
        }
        if !input.consistency.iter().zip(&rows).all(|(a, b)| a == b) {
            return Err("supplied rows disagree with the source recurrence and seeds".into());
        }
    }
    Ok(rows)
}

fn shape_columns(
    shape: &RecurrenceSpaceShape,
    limits: &RecurrenceSpaceLimits,
) -> Result<(usize, Vec<RecurrenceMonomial>), String> {
    let (target, mut cols) = match shape {
        RecurrenceSpaceShape::AutomaticVariableDegree { .. } => {
            return Err("unresolved automatic shape".into())
        }
        RecurrenceSpaceShape::Explicit {
            target_lag,
            monomials,
        } => (*target_lag, monomials.clone()),
        RecurrenceSpaceShape::Box {
            length,
            derivative_order,
            index_degree,
            variable_degree,
            leading_index_degree,
            leading_variable_degree,
        } => {
            let dimensions = [
                *length,
                derivative_order.saturating_add(1),
                index_degree.saturating_add(1),
                variable_degree.saturating_add(1),
            ];
            let rhs = dimensions.iter().try_fold(1usize, |a, b| a.checked_mul(*b));
            let lhs = leading_index_degree.checked_add(1).and_then(|n| {
                leading_variable_degree
                    .checked_add(1)
                    .and_then(|x| n.checked_mul(x))
            });
            if rhs
                .and_then(|r| lhs.and_then(|l| r.checked_add(l)))
                .is_none_or(|n| n > limits.max_columns)
            {
                return Err("column budget".into());
            }
            let mut cols = Vec::new();
            for n in 0..=*leading_index_degree {
                for x in 0..=*leading_variable_degree {
                    cols.push(RecurrenceMonomial {
                        lag: *length,
                        derivative_order: 0,
                        index_degree: n,
                        variable_degree: x,
                    });
                }
            }
            for lag in 0..*length {
                for d in 0..=*derivative_order {
                    for n in 0..=*index_degree {
                        for x in 0..=*variable_degree {
                            cols.push(RecurrenceMonomial {
                                lag,
                                derivative_order: d,
                                index_degree: n,
                                variable_degree: x,
                            });
                        }
                    }
                }
            }
            (*length, cols)
        }
    };
    if target == 0
        || cols
            .iter()
            .any(|c| c.lag > target || (c.lag == target && c.derivative_order != 0))
    {
        return Err("invalid successor support".into());
    }
    if cols.len() > limits.max_columns
        || cols.iter().any(|c| {
            [c.index_degree, c.variable_degree, c.derivative_order]
                .iter()
                .any(|d| *d > limits.max_degree)
        })
    {
        return Err("support degree/column budget".into());
    }
    cols.sort();
    if cols.windows(2).any(|w| w[0] == w[1]) {
        return Err("duplicate support monomial".into());
    }
    Ok((target, cols))
}

// Evaluate directly on rows; no generating-function numerator construction.
fn column_at(
    c: &RecurrenceMonomial,
    rows: &[Poly],
    rank: usize,
    origin: usize,
    limits: &RecurrenceSpaceLimits,
) -> Result<Poly, String> {
    let row = rank
        .checked_sub(origin)
        .and_then(|n| n.checked_add(c.lag))
        .and_then(|i| rows.get(i))
        .ok_or("missing row")?;
    let mut p = if c.derivative_order > row.degree().unwrap_or(0) {
        Poly::zero()
    } else {
        row.nth_derivative(c.derivative_order)
    };
    if p.degree().unwrap_or(0).saturating_add(c.variable_degree) > limits.max_degree {
        return Err("residual degree budget".into());
    }
    let power = u32::try_from(c.index_degree).map_err(|_| "index degree overflow")?;
    let factor = Q::from_integer(BigInt::from(rank).pow(power));
    p = p.scale(&factor) * Poly::monomial(Q::one(), c.variable_degree);
    poly_budget(&p, limits)?;
    Ok(p)
}

fn residual_matrix(
    rows: &[Poly],
    origin: usize,
    start: usize,
    count: usize,
    target: usize,
    cols: &[RecurrenceMonomial],
    lim: &RecurrenceSpaceLimits,
) -> Result<Vec<Vec<Q>>, String> {
    let mut matrix = Vec::new();
    let mut operations = lim.max_operations;
    for offset in 0..count {
        let n = start.checked_add(offset).ok_or("rank overflow")?;
        for c in cols {
            let p = rows
                .get(
                    n.checked_sub(origin)
                        .and_then(|i| i.checked_add(c.lag))
                        .ok_or("rank overflow")?,
                )
                .ok_or("missing row")?;
            let work = p
                .coeffs()
                .len()
                .checked_mul(c.derivative_order.saturating_add(2))
                .ok_or("residual work overflow")?;
            operations = operations
                .checked_sub(work)
                .ok_or("residual arithmetic budget")?;
        }
        let values: Vec<_> = cols
            .iter()
            .map(|c| {
                column_at(c, rows, n, origin, lim).map(|p| {
                    if c.lag == target {
                        p
                    } else {
                        p.scale(&(-Q::one()))
                    }
                })
            })
            .collect::<Result<_, _>>()?;
        let degree = values.iter().filter_map(Poly::degree).max();
        if let Some(degree) = degree {
            if matrix
                .len()
                .checked_add(degree + 1)
                .and_then(|r| r.checked_mul(cols.len()))
                .is_none_or(|x| x > lim.max_matrix_entries)
            {
                return Err("matrix entry budget".into());
            }
            for x in 0..=degree {
                matrix.push(values.iter().map(|p| p.coeff(x)).collect());
            }
        }
    }
    Ok(matrix)
}

fn dot(a: &[Q], b: &[Q]) -> Q {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn singular_ranks(
    cols: &[RecurrenceMonomial],
    vector: &[Q],
    target: usize,
) -> (Option<Vec<String>>, bool) {
    let mut by_x: BTreeMap<usize, Vec<Q>> = BTreeMap::new();
    for (c, v) in cols
        .iter()
        .zip(vector)
        .filter(|(c, v)| c.lag == target && !v.is_zero())
    {
        let a = by_x.entry(c.variable_degree).or_default();
        a.resize(a.len().max(c.index_degree + 1), Q::zero());
        a[c.index_degree] += v;
    }
    let requires = by_x.keys().any(|&x| x > 0);
    (
        integer_zeros_of_gcd(
            &by_x
                .values()
                .map(|p| Poly::new(p.clone()))
                .collect::<Vec<_>>(),
        ),
        requires,
    )
}

fn integer_zeros_of_gcd(polynomials: &[Poly]) -> Option<Vec<String>> {
    if polynomials.iter().any(|p| p.degree() == Some(0)) {
        return Some(vec![]);
    }
    if polynomials.iter().any(|p| {
        p.degree().unwrap_or(0) > 8
            || p.coeffs()
                .iter()
                .any(|c| c.numer().bits().max(c.denom().bits()) > 128)
    }) {
        return None;
    }
    let mut gcd = Poly::zero();
    for p in polynomials {
        gcd = gcd.gcd(p);
    }
    if gcd.is_zero()
        || gcd.degree().unwrap_or(0) > 8
        || gcd
            .coeffs()
            .iter()
            .any(|c| c.numer().bits().max(c.denom().bits()) > 128)
    {
        return None;
    }
    let denominator = gcd
        .coeffs()
        .iter()
        .fold(BigInt::one(), |a, c| a.lcm(c.denom()));
    let integer: Vec<_> = gcd
        .coeffs()
        .iter()
        .map(|c| c.numer() * (&denominator / c.denom()))
        .collect();
    let chain = SturmChain::from_bigint_coeffs(&integer);
    let mut roots = BTreeSet::new();
    for (lo, hi) in chain.isolate_roots(&Q::new(BigInt::one(), BigInt::from(2))) {
        let first = lo.numer().div_ceil(lo.denom());
        let last = hi.numer().div_floor(hi.denom());
        if first <= last
            && !first.is_negative()
            && gcd.evaluate(&Q::from_integer(first.clone())).is_zero()
        {
            roots.insert(first);
        }
    }
    Some(roots.into_iter().map(|n| n.to_string()).collect())
}

fn grouped_coefficients(
    cols: &[RecurrenceMonomial],
    vector: &[Q],
) -> BTreeMap<(usize, usize), BivarPoly> {
    let mut groups: BTreeMap<(usize, usize), Vec<Vec<Q>>> = BTreeMap::new();
    for (c, v) in cols.iter().zip(vector).filter(|(_, v)| !v.is_zero()) {
        let p = groups.entry((c.lag, c.derivative_order)).or_default();
        p.resize(p.len().max(c.index_degree + 1), Vec::new());
        let width = p[c.index_degree].len().max(c.variable_degree + 1);
        p[c.index_degree].resize(width, Q::zero());
        p[c.index_degree][c.variable_degree] += v;
    }
    groups
        .into_iter()
        .map(|(k, p)| (k, BivarPoly::from_coefficients_normalized(p)))
        .collect()
}

fn representative(
    cols: &[RecurrenceMonomial],
    vector: &[Q],
    target: usize,
    index: usize,
) -> RecurrenceSpaceRepresentative {
    let (singular_ranks, requires_x_division) = singular_ranks(cols, vector, target);
    let groups = grouped_coefficients(cols, vector);
    let side = |leading: bool| {
        groups
            .iter()
            .filter(|((lag, _), _)| (*lag == target) == leading)
            .map(|((lag, d), p)| {
                let row = if *d == 0 {
                    format!("P[n+{lag}]")
                } else {
                    format!("D^{d}(P[n+{lag}])")
                };
                format!("({}) {row}", p.to_string().replace('t', "X"))
            })
            .collect::<Vec<_>>()
            .join(" + ")
    };
    let rhs = side(false);
    RecurrenceSpaceRepresentative {
        basis_index: index,
        coefficients: vector.iter().map(ToString::to_string).collect(),
        formula: format!(
            "{} = {}",
            side(true),
            if rhs.is_empty() { "0" } else { &rhs }
        ),
        singular_ranks,
        requires_x_division,
    }
}

fn empty_result(target: usize, cols: Vec<RecurrenceMonomial>, h: usize) -> RecurrenceSpaceResult {
    RecurrenceSpaceResult {
        status: "incomplete_budget".into(),
        evidence: "none".into(),
        reason: None,
        target_lag: target,
        monomials: cols,
        h,
        residual_count: 0,
        last_required_source_row: None,
        equations: 0,
        successor_projection_rank: 0,
        certificate: None,
        representatives: Vec::new(),
        all_checked_rows_zero: false,
        automatic_degree_bound: None,
        propagation: None,
        leader_coverage: None,
    }
}

fn shape_target(shape: &RecurrenceSpaceShape) -> usize {
    match shape {
        RecurrenceSpaceShape::Box { length, .. }
        | RecurrenceSpaceShape::AutomaticVariableDegree { length, .. } => *length,
        RecurrenceSpaceShape::Explicit { target_lag, .. } => *target_lag,
    }
}
struct ShapeData {
    cols: Vec<RecurrenceMonomial>,
    target: usize,
    h: usize,
    count: usize,
    rows: Vec<Poly>,
    required: usize,
    automatic: Option<RecurrenceAutomaticDegreeBound>,
    propagation: Option<RecurrencePropagationCertificate>,
}
fn prepare_shape(
    shape: &RecurrenceSpaceShape,
    input: &Input,
    request: &RecurrenceSpaceRequest,
    supplied: Option<&RecurrencePropagationCertificate>,
    verifying: bool,
) -> Result<ShapeData, String> {
    let lim = &request.limits;
    let mut resolved = match shape {
        RecurrenceSpaceShape::AutomaticVariableDegree {
            length,
            derivative_order,
            index_degree,
            leading_index_degree,
            leading_variable_degree,
        } => {
            if input.a.is_none() && input.indexed_a.is_none() {
                return Err("automatic degree requires a defining source".into());
            }
            RecurrenceSpaceShape::Box {
                length: *length,
                derivative_order: *derivative_order,
                index_degree: *index_degree,
                variable_degree: 0,
                leading_index_degree: *leading_index_degree,
                leading_variable_degree: *leading_variable_degree,
            }
        }
        other => other.clone(),
    };
    let (target, mut cols) = shape_columns(&resolved, lim)?;
    let h = cols
        .iter()
        .map(|c| c.index_degree.checked_add(c.derivative_order))
        .collect::<Option<Vec<_>>>()
        .ok_or("support degree overflow")?
        .into_iter()
        .max()
        .unwrap_or(0)
        .checked_add(1)
        .ok_or("h overflow")?;
    let propagation = if let Some(a) = &input.indexed_a {
        if verifying {
            let cert = supplied.ok_or("missing indexed-source propagation certificate")?;
            propagation::verify(a, &cols, request.domain_start, cert, lim)?;
            Some(cert.clone())
        } else {
            Some(propagation::build(a, &cols, request.domain_start, lim)?)
        }
    } else {
        if supplied.is_some() {
            return Err("unexpected propagation certificate".into());
        }
        None
    };
    let offset = request.domain_start - input.origin;
    let count = if let Some(cert) = &propagation {
        cert.regular_from
            .checked_sub(request.domain_start)
            .and_then(|n| n.checked_add(cert.coefficients.len() - 1))
            .ok_or("propagation bound overflow")?
    } else if let Some(a) = &input.a {
        h.checked_mul(a.len()).ok_or("residual bound overflow")?
    } else {
        input
            .rows
            .len()
            .saturating_sub(offset.saturating_add(target))
    };
    let required = offset
        .checked_add(count)
        .and_then(|n| n.checked_add(target))
        .ok_or("row bound overflow")?;
    let rows = source_rows(input, required, lim)?;
    let mut automatic = None;
    if matches!(shape, RecurrenceSpaceShape::AutomaticVariableDegree { .. }) {
        let mut rhs_degree = 0;
        let mut leader_degree = 0;
        let mut rhs_nonzero = false;
        let mut work = lim.max_operations;
        for n in request.domain_start
            ..request
                .domain_start
                .checked_add(count)
                .ok_or("rank overflow")?
        {
            for c in &cols {
                let row = rows
                    .get(n - input.origin + c.lag)
                    .ok_or("missing bound row")?;
                let cost = row
                    .coeffs()
                    .len()
                    .checked_mul(c.derivative_order.saturating_add(2))
                    .ok_or("bound arithmetic budget")?;
                work = work.checked_sub(cost).ok_or("bound arithmetic budget")?;
                let p = column_at(c, &rows, n, input.origin, lim)?;
                if c.lag == target {
                    leader_degree = leader_degree.max(p.degree().unwrap_or(0));
                } else {
                    rhs_degree = rhs_degree.max(p.degree().unwrap_or(0));
                    rhs_nonzero |= !p.is_zero();
                }
            }
        }
        let rank_upper_bound = if rhs_nonzero {
            count.min(cols.iter().filter(|c| c.lag < target).count())
        } else {
            0
        };
        let variable_degree = if rank_upper_bound == 0 {
            0
        } else {
            let first = leader_degree
                .checked_add(
                    (rank_upper_bound - 1)
                        .checked_mul(rhs_degree)
                        .ok_or("automatic degree budget")?,
                )
                .ok_or("automatic degree budget")?;
            let second = rank_upper_bound
                .checked_mul(rhs_degree)
                .ok_or("automatic degree budget")?
                .saturating_sub(1);
            first.max(second)
        };
        automatic = Some(RecurrenceAutomaticDegreeBound {
            rhs_degree,
            leader_degree,
            rank_upper_bound,
            variable_degree,
        });
        if let RecurrenceSpaceShape::Box {
            variable_degree: v, ..
        } = &mut resolved
        {
            *v = variable_degree;
        }
        cols = shape_columns(&resolved, lim)
            .map_err(|error| format!("automatic RHS X-degree {variable_degree}: {error}"))?
            .1;
    }
    Ok(ShapeData {
        cols,
        target,
        h,
        count,
        rows,
        required,
        automatic,
        propagation,
    })
}

fn leader_coverage(
    result: &RecurrenceSpaceResult,
    limits: &RecurrenceSpaceLimits,
) -> Result<Option<RecurrenceLeaderCoverage>, String> {
    if result.representatives.is_empty() {
        return Ok(None);
    }
    let mut combined: BTreeMap<RecurrenceMonomial, Q> = BTreeMap::new();
    let mut multipliers = Vec::new();
    let mut leaders = Vec::new();
    let mut remaining = limits.max_operations;
    for rep in &result.representatives {
        let v = rep
            .coefficients
            .iter()
            .map(|s| checked_q(s, limits))
            .collect::<Result<Vec<_>, _>>()?;
        let group = grouped_coefficients(&result.monomials, &v);
        let leader = &group[&(result.target_lag, 0)];
        multipliers.push(encode_matrix(&leader.coeffs));
        let width = leader.coeffs.iter().map(Vec::len).max().unwrap_or(0);
        for x in 0..width {
            leaders.push(Poly::new(
                leader
                    .coeffs
                    .iter()
                    .map(|r| r.get(x).cloned().unwrap_or_else(Q::zero))
                    .collect(),
            ));
        }
        for (a, ac) in result
            .monomials
            .iter()
            .zip(&v)
            .filter(|(a, c)| a.lag == result.target_lag && !c.is_zero())
        {
            for (b, bc) in result
                .monomials
                .iter()
                .zip(&v)
                .filter(|(_, c)| !c.is_zero())
            {
                remaining = remaining
                    .checked_sub(1)
                    .ok_or("leader combination arithmetic budget")?;
                let n = a
                    .index_degree
                    .checked_add(b.index_degree)
                    .ok_or("leader combination degree budget")?;
                let x = a
                    .variable_degree
                    .checked_add(b.variable_degree)
                    .ok_or("leader combination degree budget")?;
                if n > limits.max_degree || x > limits.max_degree {
                    return Err("leader combination degree budget".into());
                }
                let key = RecurrenceMonomial {
                    index_degree: n,
                    variable_degree: x,
                    ..*b
                };
                let c = combined.entry(key).or_insert_with(Q::zero);
                *c += ac * bc;
                if c.numer().bits().max(c.denom().bits()) > limits.max_coefficient_bits {
                    return Err("leader combination coefficient budget".into());
                }
                if combined.len() > limits.max_matrix_entries {
                    return Err("leader combination entry budget".into());
                }
            }
        }
    }
    combined.retain(|_, v| !v.is_zero());
    let (cols, vector): (Vec<_>, Vec<_>) = combined.into_iter().unzip();
    let rep = representative(&cols, &vector, result.target_lag, 0);
    Ok(Some(RecurrenceLeaderCoverage {
        basis_indices: result
            .representatives
            .iter()
            .map(|p| p.basis_index)
            .collect(),
        multipliers,
        monomials: cols,
        coefficients: rep.coefficients,
        formula: rep.formula,
        common_singular_ranks: integer_zeros_of_gcd(&leaders),
        requires_x_division: rep.requires_x_division,
    }))
}

fn kernel_limits(l: &RecurrenceSpaceLimits) -> RationalKernelLimits {
    RationalKernelLimits {
        max_entries: l.max_matrix_entries,
        max_operations: l.max_operations,
        max_coefficient_bits: l.max_coefficient_bits,
    }
}

/// Return every identity in each requested finite coefficient support.
/// Positive certificates are conditional on the source if supplied; otherwise
/// they describe finite fits. Negative dual certificates need only the rows.
pub fn find_recurrence_spaces(
    request: &RecurrenceSpaceRequest,
) -> Result<RecurrenceSpaces, String> {
    let input = prepare(request)?;
    // A contradictory source is invalid input, not an incomplete search.
    source_rows(&input, input.rows.len(), &request.limits)?;
    let mut results = Vec::new();
    for shape in &request.shapes {
        let data = match prepare_shape(shape, &input, request, None, false) {
            Ok(data) => data,
            Err(error) if error.contains("budget") => {
                let mut result = empty_result(shape_target(shape), Vec::new(), 0);
                result.reason = Some(error);
                results.push(result);
                continue;
            }
            Err(error) => return Err(error),
        };
        let ShapeData {
            target,
            cols,
            h,
            count,
            rows,
            required,
            automatic,
            propagation,
        } = data;
        let mut result = empty_result(target, cols.clone(), h);
        result.automatic_degree_bound = automatic;
        result.propagation = propagation;
        result.residual_count = count;
        result.last_required_source_row = required
            .checked_sub(1)
            .and_then(|n| input.origin.checked_add(n));
        if count == 0 && result.propagation.is_none() {
            result.status = "insufficient_rows".into();
            result.reason = Some("no complete residual rows in the declared domain".into());
            results.push(result);
            continue;
        }
        result.all_checked_rows_zero = rows.iter().all(Poly::is_zero);
        let matrix = match residual_matrix(
            &rows,
            input.origin,
            request.domain_start,
            count,
            target,
            &cols,
            &request.limits,
        ) {
            Ok(m) => m,
            Err(e) => {
                result.reason = Some(e);
                results.push(result);
                continue;
            }
        };
        result.equations = matrix.len();
        let kernel = match rational_kernel_with_certificate(
            &matrix,
            cols.len(),
            kernel_limits(&request.limits),
        ) {
            Ok(k) => k,
            Err(RationalKernelError::BudgetExceeded(e)) => {
                result.reason = Some(e.into());
                results.push(result);
                continue;
            }
            Err(_) => return Err("invalid residual matrix".into()),
        };
        let leading: Vec<_> = cols
            .iter()
            .enumerate()
            .filter_map(|(i, c)| (c.lag == target).then_some(i))
            .collect();
        let projection: Vec<_> = leading
            .iter()
            .map(|&i| kernel.basis.iter().map(|v| v[i].clone()).collect())
            .collect();
        let projection_kernel = match rational_kernel_with_certificate(
            &projection,
            kernel.basis.len(),
            kernel_limits(&request.limits),
        ) {
            Ok(k) => k,
            Err(_) => {
                result.reason = Some("projection budget".into());
                results.push(result);
                continue;
            }
        };
        result.successor_projection_rank = projection_kernel.pivot_columns.len();
        let absent = result.successor_projection_rank == 0;
        let mut duals = Vec::new();
        if absent {
            for &i in &leading {
                let p = kernel
                    .pivot_columns
                    .iter()
                    .position(|&c| c == i)
                    .ok_or("missing dual pivot")?;
                duals.push(RecurrenceSpaceDual {
                    coordinate: i,
                    weights: kernel.inverse_minor[p]
                        .iter()
                        .map(ToString::to_string)
                        .collect(),
                });
            }
        }
        // An independent set of projected basis columns supplies several natural
        // successor choices; the full kernel remains available for combinations.
        for &i in &projection_kernel.pivot_columns {
            result
                .representatives
                .push(representative(&cols, &kernel.basis[i], target, i));
        }
        result.status = if absent {
            "absent_in_shape"
        } else {
            "successor_exists"
        }
        .into();
        result.evidence = if absent {
            "rows_only_obstruction"
        } else if input.a.is_some() || input.indexed_a.is_some() {
            "conditional_on_source"
        } else {
            "exact_fit"
        }
        .into();
        result.certificate = Some(RecurrenceSpaceCertificate {
            basis: encode_matrix(&kernel.basis),
            free_columns: kernel.free_columns,
            pivot_columns: kernel.pivot_columns,
            independent_rows: kernel.independent_rows,
            inverse_minor: encode_matrix(&kernel.inverse_minor),
            absence_duals: duals,
        });
        if request.include_leader_coverage {
            result.leader_coverage = leader_coverage(&result, &request.limits).unwrap_or(None);
        }
        results.push(result);
    }
    let frontier = frontier(request, &results);
    Ok(RecurrenceSpaces {
        schema: "polytool.recurrence-spaces.v1".into(),
        request: request.clone(),
        results,
        frontier,
    })
}

fn frontier(request: &RecurrenceSpaceRequest, results: &[RecurrenceSpaceResult]) -> Vec<usize> {
    let key = |i: usize| match request.shapes[i] {
        RecurrenceSpaceShape::Box {
            length,
            derivative_order,
            index_degree,
            variable_degree,
            leading_index_degree,
            leading_variable_degree,
        } => Some((
            length,
            derivative_order,
            (
                index_degree,
                variable_degree,
                leading_index_degree,
                leading_variable_degree,
            ),
        )),
        RecurrenceSpaceShape::AutomaticVariableDegree {
            length,
            derivative_order,
            index_degree,
            leading_index_degree,
            leading_variable_degree,
        } => Some((
            length,
            derivative_order,
            (
                index_degree,
                usize::MAX,
                leading_index_degree,
                leading_variable_degree,
            ),
        )),
        _ => None,
    };
    (0..results.len())
        .filter(|&i| {
            results[i].status == "successor_exists"
                && !(0..results.len()).any(|j| {
                    if results[j].status != "successor_exists" {
                        return false;
                    }
                    match (key(i), key(j)) {
                        (Some((r, s, b)), Some((r2, s2, b2))) => {
                            b == b2 && r2 <= r && s2 <= s && (r2 < r || s2 < s)
                        }
                        _ => false,
                    }
                })
        })
        .collect()
}

/// Replay a report without rerunning nullspace search. The caller supplies
/// the expected request so changed assumptions cannot pass unnoticed.
/// This is exact certificate checking, not a Lean proof or authentication.
pub fn verify_recurrence_spaces(
    request: &RecurrenceSpaceRequest,
    report: &RecurrenceSpaces,
) -> Result<(), String> {
    if report.schema != "polytool.recurrence-spaces.v1"
        || &report.request != request
        || report.results.len() != request.shapes.len()
    {
        return Err("report/request mismatch".into());
    }
    let input = prepare(request)?;
    source_rows(&input, input.rows.len(), &request.limits)?;
    for (shape, result) in request.shapes.iter().zip(&report.results) {
        if result.certificate.is_none() {
            if !["incomplete_budget", "insufficient_rows"].contains(&result.status.as_str())
                || result.evidence != "none"
                || !result.representatives.is_empty()
                || result.leader_coverage.is_some()
            {
                return Err("uncertified claim".into());
            }
            continue;
        }
        let ShapeData {
            target,
            cols,
            h,
            count,
            rows,
            required,
            automatic,
            propagation,
        } = prepare_shape(shape, &input, request, result.propagation.as_ref(), true)?;
        if target != result.target_lag
            || cols != result.monomials
            || h != result.h
            || automatic != result.automatic_degree_bound
            || propagation != result.propagation
        {
            return Err("support or bound mismatch".into());
        }
        if (count == 0 && propagation.is_none())
            || count != result.residual_count
            || result.last_required_source_row
                != required
                    .checked_sub(1)
                    .and_then(|n| input.origin.checked_add(n))
        {
            return Err("residual range mismatch".into());
        }
        if result.all_checked_rows_zero != rows.iter().all(Poly::is_zero) {
            return Err("degeneracy mismatch".into());
        }
        let matrix = residual_matrix(
            &rows,
            input.origin,
            request.domain_start,
            count,
            target,
            &cols,
            &request.limits,
        )?;
        if matrix.len() != result.equations {
            return Err("equation count mismatch".into());
        }
        let cert = result.certificate.as_ref().unwrap();
        let basis = parse_matrix(&cert.basis, &request.limits)?;
        let inverse = parse_matrix(&cert.inverse_minor, &request.limits)?;
        let rank = cert.pivot_columns.len();
        let dimension = basis.len();
        let width = cols.len();
        let verify_work = matrix
            .len()
            .checked_mul(width)
            .and_then(|x| x.checked_mul(dimension + rank + 1))
            .and_then(|x| rank.checked_pow(3).and_then(|y| x.checked_add(y)))
            .ok_or("verification work overflow")?;
        if verify_work > request.limits.max_operations {
            return Err("certificate verification arithmetic budget".into());
        }
        let partition: BTreeSet<_> = cert
            .pivot_columns
            .iter()
            .chain(&cert.free_columns)
            .copied()
            .collect();
        if rank + dimension != width
            || cert.free_columns.len() != dimension
            || partition != (0..width).collect()
            || cert.independent_rows.len() != rank
            || inverse.len() != rank
            || inverse.iter().any(|r| r.len() != rank)
            || cert.independent_rows.iter().any(|&i| i >= matrix.len())
            || basis.iter().any(|v| v.len() != width)
        {
            return Err("invalid kernel dimensions".into());
        }
        for (i, v) in basis.iter().enumerate() {
            if matrix.iter().any(|r| !dot(r, v).is_zero()) {
                return Err("invalid null vector".into());
            }
            for (j, &c) in cert.free_columns.iter().enumerate() {
                if v[c] != if i == j { Q::one() } else { Q::zero() } {
                    return Err("invalid basis independence witness".into());
                }
            }
        }
        for (i, inverse_row) in inverse.iter().enumerate() {
            for j in 0..rank {
                let sum: Q = (0..rank)
                    .map(|t| {
                        &inverse_row[t] * &matrix[cert.independent_rows[t]][cert.pivot_columns[j]]
                    })
                    .sum();
                if sum != if i == j { Q::one() } else { Q::zero() } {
                    return Err("invalid rank minor inverse".into());
                }
            }
        }
        let leading: Vec<_> = cols
            .iter()
            .enumerate()
            .filter_map(|(i, c)| (c.lag == target).then_some(i))
            .collect();
        let projection: Vec<_> = leading
            .iter()
            .map(|&i| basis.iter().map(|v| v[i].clone()).collect())
            .collect();
        // Only the small target projection is reduced; the full kernel is
        // checked by multiplication and the supplied rank minor.
        let projected = rational_kernel_with_certificate(
            &projection,
            dimension,
            kernel_limits(&request.limits),
        )
        .map_err(|_| "projection verification budget")?;
        let projection_rank = projected.pivot_columns.len();
        if projection_rank != result.successor_projection_rank {
            return Err("successor projection mismatch".into());
        }
        let absent = projection_rank == 0;
        let evidence = if absent {
            "rows_only_obstruction"
        } else if input.a.is_some() || input.indexed_a.is_some() {
            "conditional_on_source"
        } else {
            "exact_fit"
        };
        if result.status
            != if absent {
                "absent_in_shape"
            } else {
                "successor_exists"
            }
            || result.evidence != evidence
        {
            return Err("claim/evidence mismatch".into());
        }
        if absent {
            if cert.absence_duals.len() != leading.len() {
                return Err("missing absence duals".into());
            }
            for (&coordinate, dual) in leading.iter().zip(&cert.absence_duals) {
                let weights: Vec<_> = dual
                    .weights
                    .iter()
                    .map(|s| checked_q(s, &request.limits))
                    .collect::<Result<_, _>>()?;
                if dual.coordinate != coordinate || weights.len() != rank {
                    return Err("invalid dual dimensions".into());
                }
                for (c, _) in cols.iter().enumerate() {
                    let sum: Q = (0..rank)
                        .map(|i| &weights[i] * &matrix[cert.independent_rows[i]][c])
                        .sum();
                    if sum != if c == coordinate { Q::one() } else { Q::zero() } {
                        return Err("invalid absence dual".into());
                    }
                }
            }
        } else if !cert.absence_duals.is_empty() {
            return Err("spurious absence duals".into());
        }
        let expected: Vec<_> = projected
            .pivot_columns
            .iter()
            .map(|&i| representative(&cols, &basis[i], target, i))
            .collect();
        if result.representatives != expected {
            return Err("representative or singular-rank mismatch".into());
        }
        if result.leader_coverage.is_some()
            && (!request.include_leader_coverage
                || result.leader_coverage != leader_coverage(result, &request.limits)?)
        {
            return Err("leader coverage mismatch".into());
        }
    }
    if report.frontier != frontier(request, &report.results) {
        return Err("frontier mismatch".into());
    }
    Ok(())
}

#[cfg(test)]
mod coverage_tests {
    use super::*;
    #[test]
    fn two_singular_identities_combine_to_a_regular_one() {
        let cols = vec![
            RecurrenceMonomial {
                lag: 0,
                derivative_order: 0,
                index_degree: 0,
                variable_degree: 0,
            },
            RecurrenceMonomial {
                lag: 0,
                derivative_order: 0,
                index_degree: 1,
                variable_degree: 0,
            },
            RecurrenceMonomial {
                lag: 1,
                derivative_order: 0,
                index_degree: 0,
                variable_degree: 0,
            },
            RecurrenceMonomial {
                lag: 1,
                derivative_order: 0,
                index_degree: 1,
                variable_degree: 0,
            },
        ];
        let mut result = empty_result(1, cols.clone(), 2);
        for (i, values) in [[0, 1, 0, 1], [-1, 1, -1, 1]].into_iter().enumerate() {
            let v = values
                .into_iter()
                .map(|n| Q::from_integer(n.into()))
                .collect::<Vec<_>>();
            result.representatives.push(representative(&cols, &v, 1, i));
        }
        assert_eq!(
            result.representatives[0].singular_ranks,
            Some(vec!["0".into()])
        );
        assert_eq!(
            result.representatives[1].singular_ranks,
            Some(vec!["1".into()])
        );
        let combined = leader_coverage(&result, &Default::default())
            .unwrap()
            .unwrap();
        assert_eq!(combined.common_singular_ranks, Some(vec![]));
        let rec = combined.to_recurrence(1).unwrap();
        assert_eq!(
            rec.denominator.unwrap().coeffs,
            vec![
                vec![Q::from_integer(5.into())],
                vec![Q::from_integer((-6).into())],
                vec![Q::from_integer(2.into())]
            ]
        );
        // The backward target index is m=n+1: n^2+(n-1)^2 becomes 2m^2-6m+5.
    }
}
