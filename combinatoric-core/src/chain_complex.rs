//! Finite based integer chain complexes, unit cancellation, and abstract
//! integral homology groups.

use crate::{
    smith_normal_form, verify_smith_certificate, MutableSparseMatrix, SmithError, SmithOptions,
    SmithReplayOptions, SparseMatrix, SparseMatrixLimits,
};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use sha2::{Digest, Sha256};
use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ChainComplexError {
    DifferentialShape {
        degree: i32,
        expected: (usize, usize),
        actual: (usize, usize),
    },
    DifferentialDoesNotSquare {
        left_degree: i32,
        row: usize,
        column: usize,
        value: BigInt,
    },
    MissingDegree {
        degree: i32,
    },
    UnsupportedDegree {
        degree: i32,
    },
    NonUnitPivot {
        degree: i32,
        row: usize,
        column: usize,
        value: BigInt,
    },
    ReductionLimit {
        kind: &'static str,
        observed: usize,
        budget: usize,
    },
    Certificate(String),
    Smith(SmithError),
}

impl std::fmt::Display for ChainComplexError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::DifferentialShape {
                degree,
                expected,
                actual,
            } => write!(f, "D_{degree} has shape {actual:?}, expected {expected:?}"),
            Self::DifferentialDoesNotSquare {
                left_degree,
                row,
                column,
                value,
            } => write!(
                f,
                "D_{left_degree} D_{} has nonzero ({row}, {column}) = {value}",
                left_degree + 1
            ),
            Self::MissingDegree { degree } => {
                write!(f, "no generator count was supplied for degree {degree}")
            }
            Self::UnsupportedDegree { degree } => write!(
                f,
                "degree {degree} is outside the supported adjacent-map range"
            ),
            Self::NonUnitPivot {
                degree,
                row,
                column,
                value,
            } => write!(f, "D_{degree}[{row}, {column}] = {value} is not a unit"),
            Self::ReductionLimit {
                kind,
                observed,
                budget,
            } => write!(f, "{kind} limit {budget} exceeded by {observed}"),
            Self::Certificate(message) => write!(f, "invalid cancellation certificate: {message}"),
            Self::Smith(error) => write!(f, "Smith reduction failed: {error}"),
        }
    }
}
impl std::error::Error for ChainComplexError {}
impl From<SmithError> for ChainComplexError {
    fn from(value: SmithError) -> Self {
        Self::Smith(value)
    }
}

/// A finite free chain complex with `D_k: C_k -> C_(k-1)` stored by source
/// degree. Missing end maps are interpreted as shaped zero matrices.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FiniteChainComplex {
    generator_counts: BTreeMap<i32, usize>,
    differentials: BTreeMap<i32, SparseMatrix<BigInt>>,
}

impl FiniteChainComplex {
    pub fn new(
        generator_counts: BTreeMap<i32, usize>,
        differentials: BTreeMap<i32, SparseMatrix<BigInt>>,
    ) -> Result<Self, ChainComplexError> {
        if let Some(&degree) = generator_counts
            .keys()
            .find(|degree| **degree == i32::MIN || **degree == i32::MAX)
        {
            return Err(ChainComplexError::UnsupportedDegree { degree });
        }
        let result = Self {
            generator_counts,
            differentials,
        };
        result.validate()?;
        Ok(result)
    }
    pub fn generator_counts(&self) -> &BTreeMap<i32, usize> {
        &self.generator_counts
    }
    pub fn count(&self, degree: i32) -> usize {
        self.generator_counts.get(&degree).copied().unwrap_or(0)
    }
    pub fn differential(&self, degree: i32) -> Option<&SparseMatrix<BigInt>> {
        self.differentials.get(&degree)
    }
    pub fn differential_or_zero(&self, degree: i32) -> SparseMatrix<BigInt> {
        self.differentials.get(&degree).cloned().unwrap_or_else(|| {
            let previous = degree.checked_sub(1).map_or(0, |value| self.count(value));
            SparseMatrix::zero(previous, self.count(degree))
        })
    }
    /// Budgeted shaped-zero retrieval for callers handling untrusted degrees or
    /// metadata-only complexes.  Stored matrices retain their existing shape;
    /// missing maps are allocated only after the CSR shape budget is checked.
    pub fn differential_or_zero_with_limits(
        &self,
        degree: i32,
        limits: SparseMatrixLimits,
    ) -> Result<SparseMatrix<BigInt>, ChainComplexError> {
        if let Some(matrix) = self.differentials.get(&degree) {
            let slots = matrix
                .rows()
                .checked_add(1)
                .ok_or(ChainComplexError::ReductionLimit {
                    kind: "matrix shape slots",
                    observed: usize::MAX,
                    budget: limits.max_shape_slots,
                })?;
            if slots > limits.max_shape_slots {
                return Err(ChainComplexError::ReductionLimit {
                    kind: "matrix shape slots",
                    observed: slots,
                    budget: limits.max_shape_slots,
                });
            }
            if matrix.nnz() > limits.max_nnz {
                return Err(ChainComplexError::ReductionLimit {
                    kind: "NNZ",
                    observed: matrix.nnz(),
                    budget: limits.max_nnz,
                });
            }
            return Ok(matrix.clone());
        }
        let previous = degree
            .checked_sub(1)
            .ok_or(ChainComplexError::UnsupportedDegree { degree })?;
        SparseMatrix::zero_with_limits(self.count(previous), self.count(degree), limits).map_err(
            |error| match error {
                crate::SparseMatrixError::ShapeBudgetExceeded { slots, budget, .. } => {
                    ChainComplexError::ReductionLimit {
                        kind: "matrix shape slots",
                        observed: slots,
                        budget,
                    }
                }
                error => ChainComplexError::Certificate(error.to_string()),
            },
        )
    }
    pub fn degrees(&self) -> impl Iterator<Item = i32> + '_ {
        self.generator_counts.keys().copied()
    }

    pub fn validate(&self) -> Result<(), ChainComplexError> {
        for (&degree, differential) in &self.differentials {
            let Some(previous) = degree.checked_sub(1) else {
                return Err(ChainComplexError::UnsupportedDegree { degree });
            };
            if !self.generator_counts.contains_key(&degree)
                || !self.generator_counts.contains_key(&previous)
            {
                return Err(ChainComplexError::MissingDegree { degree });
            }
            let expected = (self.count(previous), self.count(degree));
            if differential.shape() != expected {
                return Err(ChainComplexError::DifferentialShape {
                    degree,
                    expected,
                    actual: differential.shape(),
                });
            }
        }
        for &degree in self.differentials.keys() {
            let Some(previous) = degree.checked_sub(1) else {
                return Err(ChainComplexError::UnsupportedDegree { degree });
            };
            let (Some(left), Some(right)) = (
                self.differentials.get(&previous),
                self.differentials.get(&degree),
            ) else {
                continue;
            };
            if let Err((row, column, value)) = SparseMatrix::compose_is_zero(left, right)
                .expect("adjacent differential shapes agree")
            {
                return Err(ChainComplexError::DifferentialDoesNotSquare {
                    left_degree: previous,
                    row,
                    column,
                    value,
                });
            }
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct UnitReductionOptions {
    pub max_pivots: usize,
    pub max_nnz: usize,
    pub max_shape_slots: usize,
    pub max_entry_bits: u64,
    pub record_certificate: bool,
}
impl Default for UnitReductionOptions {
    fn default() -> Self {
        Self {
            max_pivots: usize::MAX,
            max_nnz: usize::MAX,
            max_shape_slots: usize::MAX,
            max_entry_bits: 16_384,
            record_certificate: false,
        }
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct UnitReductionStats {
    pub initial_nnz: usize,
    pub peak_nnz: usize,
    pub final_nnz: usize,
    pub pivot_updates: usize,
    pub stale_queue_entries: usize,
    pub score_revalidations: usize,
    pub maximum_coefficient_bits: u64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct UnitPivot {
    pub degree: i32,
    pub row: usize,
    pub column: usize,
    pub value: BigInt,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct UnitCancellationCertificate {
    pub schema_version: u32,
    pub input_sha256: String,
    pub pivots: Vec<UnitPivot>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct UnitReductionResult {
    pub reduced: FiniteChainComplex,
    pub stats: UnitReductionStats,
    pub certificate: Option<UnitCancellationCertificate>,
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
struct PivotCandidate {
    score: usize,
    degree: i32,
    row: usize,
    column: usize,
}

/// Cancel unit entries in a graded complex.  Generator IDs remain stable while
/// reducing; only the returned residual compacts inactive bases.
pub fn cancel_units(
    complex: &FiniteChainComplex,
    options: UnitReductionOptions,
) -> Result<UnitReductionResult, ChainComplexError> {
    let (initial_nnz, initial_bits, shape_slots) = complex_resource_usage(complex);
    check_reduction_resources(initial_nnz, initial_bits, shape_slots, &options)?;
    let mut state = ReductionState::from_complex(complex);
    let mut stats = UnitReductionStats::default();
    stats.initial_nnz = initial_nnz;
    stats.peak_nnz = stats.initial_nnz;
    stats.maximum_coefficient_bits = initial_bits;
    let mut certificate = options
        .record_certificate
        .then(|| UnitCancellationCertificate {
            schema_version: 1,
            input_sha256: complex_sha256(complex),
            pivots: Vec::new(),
        });
    let mut queue = candidate_queue(&state);
    while let Some(Reverse(candidate)) = queue.pop() {
        let pivot = UnitPivot {
            degree: candidate.degree,
            row: candidate.row,
            column: candidate.column,
            value: state.value(candidate.degree, candidate.row, candidate.column)?,
        };
        if !state.is_active_pivot(&pivot) || !is_unit(&pivot.value) {
            stats.stale_queue_entries += 1;
            continue;
        }
        let current_score = state.pivot_score(pivot.degree, pivot.row, pivot.column)?;
        if current_score != candidate.score {
            stats.score_revalidations += 1;
            queue.push(Reverse(PivotCandidate {
                score: current_score,
                ..candidate
            }));
            continue;
        }
        if stats.pivot_updates == options.max_pivots {
            return Err(ChainComplexError::ReductionLimit {
                kind: "unit pivot",
                observed: stats.pivot_updates + 1,
                budget: options.max_pivots,
            });
        }
        let changes = state.apply_pivot(&pivot)?;
        stats.pivot_updates += 1;
        stats.peak_nnz = stats.peak_nnz.max(state.nnz());
        if state.nnz() > options.max_nnz {
            return Err(ChainComplexError::ReductionLimit {
                kind: "NNZ",
                observed: state.nnz(),
                budget: options.max_nnz,
            });
        }
        stats.maximum_coefficient_bits =
            stats.maximum_coefficient_bits.max(changes.maximum_new_bits);
        if changes.maximum_new_bits > options.max_entry_bits {
            return Err(ChainComplexError::ReductionLimit {
                kind: "coefficient bit length",
                observed: changes.maximum_new_bits as usize,
                budget: options.max_entry_bits as usize,
            });
        }
        enqueue_changed_candidates(
            &state,
            &mut queue,
            pivot.degree,
            &changes.rows,
            &changes.columns,
        );
        if let Some(certificate) = &mut certificate {
            certificate.pivots.push(pivot);
        }
    }
    stats.final_nnz = state.nnz();
    let reduced = state.compact()?;
    reduced.validate()?;
    Ok(UnitReductionResult {
        reduced,
        stats,
        certificate,
    })
}

/// Reconstruct the reduction from the original complex and a certificate.
/// Selection heuristics are not consulted, so tampering is rejected rather
/// than masked by a fresh pivot search.
pub fn replay_unit_cancellation(
    complex: &FiniteChainComplex,
    certificate: &UnitCancellationCertificate,
    options: &UnitReductionOptions,
) -> Result<FiniteChainComplex, ChainComplexError> {
    if certificate.schema_version != 1 {
        return Err(ChainComplexError::Certificate(
            "unsupported schema version".into(),
        ));
    }
    if certificate.input_sha256 != complex_sha256(complex) {
        return Err(ChainComplexError::Certificate(
            "input digest differs".into(),
        ));
    }
    if certificate.pivots.len() > options.max_pivots {
        return Err(ChainComplexError::Certificate(
            "pivot count exceeds replay limit".into(),
        ));
    }
    let (initial_nnz, initial_bits, shape_slots) = complex_resource_usage(complex);
    check_reduction_resources(initial_nnz, initial_bits, shape_slots, options)?;
    let mut state = ReductionState::from_complex(complex);
    for pivot in &certificate.pivots {
        if !state.is_active_pivot(pivot)
            || state.value(pivot.degree, pivot.row, pivot.column)? != pivot.value
            || !is_unit(&pivot.value)
        {
            return Err(ChainComplexError::Certificate(format!(
                "invalid or stale pivot at D_{}[{}, {}]",
                pivot.degree, pivot.row, pivot.column
            )));
        }
        let changes = state.apply_pivot(pivot)?;
        if state.nnz() > options.max_nnz || changes.maximum_new_bits > options.max_entry_bits {
            return Err(ChainComplexError::Certificate(
                "certificate exceeds replay resource budget".into(),
            ));
        }
    }
    let reduced = state.compact()?;
    reduced.validate()?;
    Ok(reduced)
}

/// Replay a certificate and bind it to a claimed compact residual.  The public
/// certificate type is a reusable in-memory artifact; this verifier is the
/// stable API for consumers that persist it by their own serialization layer.
pub fn replay_unit_cancellation_and_verify(
    complex: &FiniteChainComplex,
    certificate: &UnitCancellationCertificate,
    expected_residual: &FiniteChainComplex,
    options: &UnitReductionOptions,
) -> Result<(), ChainComplexError> {
    let replayed = replay_unit_cancellation(complex, certificate, options)?;
    if &replayed != expected_residual {
        return Err(ChainComplexError::Certificate(
            "replay residual differs from claimed residual".into(),
        ));
    }
    Ok(())
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AbelianGroup {
    pub free_rank: usize,
    /// Nontrivial invariant factors; e.g. `[2, 12]` means `Z/2 ⊕ Z/12`.
    pub torsion_invariants: Vec<BigInt>,
}

/// Assemble abstract integral homology groups from Smith forms of the maps.
/// It does not claim compatible cycle representatives or retain Smith bases.
pub fn integral_homology(
    complex: &FiniteChainComplex,
    smith_options: SmithOptions,
) -> Result<BTreeMap<i32, AbelianGroup>, ChainComplexError> {
    complex.validate()?;
    let mut ranks = BTreeMap::new();
    let mut factors = BTreeMap::new();
    for &degree in complex.generator_counts.keys() {
        if let Some(boundary) = complex.differential(degree) {
            let smith = smith_normal_form(boundary, smith_options.clone())?;
            if let Some(operations) = &smith.operations {
                verify_smith_certificate(
                    boundary,
                    operations,
                    &smith.invariant_factors,
                    SmithReplayOptions::from(&smith_options),
                )?;
            }
            ranks.insert(degree, smith.rank);
            factors.insert(degree, smith.invariant_factors);
        } else {
            // An absent boundary is mathematically the shaped zero map.  Its
            // rank and Smith factors are known without materializing CSR rows.
            ranks.insert(degree, 0);
            factors.insert(degree, Vec::new());
        }
    }
    Ok(complex
        .generator_counts
        .keys()
        .map(|&degree| {
            let rank = complex.count(degree)
                - ranks.get(&degree).copied().unwrap_or(0)
                - ranks.get(&(degree + 1)).copied().unwrap_or(0);
            let torsion_invariants = factors
                .get(&(degree + 1))
                .into_iter()
                .flatten()
                .filter(|factor| **factor > BigInt::one())
                .cloned()
                .collect();
            (
                degree,
                AbelianGroup {
                    free_rank: rank,
                    torsion_invariants,
                },
            )
        })
        .collect())
}

/// The field formula is intentionally separate from any particular rank engine.
pub fn field_betti_number(
    generator_count: usize,
    rank_out: usize,
    rank_in: usize,
) -> Option<usize> {
    generator_count.checked_sub(rank_out)?.checked_sub(rank_in)
}

/// The universal-coefficient dimension predicted by an exact integral group.
pub fn universal_coefficient_dimension(
    current: &AbelianGroup,
    previous: Option<&AbelianGroup>,
    prime: u64,
) -> Result<usize, ChainComplexError> {
    if !is_prime_u64(prime) {
        return Err(ChainComplexError::Certificate(format!(
            "UCT modulus {prime} is not prime"
        )));
    }
    let prime = BigInt::from(prime);
    Ok(current.free_rank
        + current
            .torsion_invariants
            .iter()
            .filter(|factor| factor.mod_floor(&prime).is_zero())
            .count()
        + previous
            .into_iter()
            .flat_map(|group| group.torsion_invariants.iter())
            .filter(|factor| factor.mod_floor(&prime).is_zero())
            .count())
}

struct ReductionState {
    matrices: BTreeMap<i32, MutableSparseMatrix>,
    active: BTreeMap<i32, Vec<bool>>,
    total_nnz: usize,
}

struct PivotChanges {
    rows: Vec<usize>,
    columns: Vec<usize>,
    maximum_new_bits: u64,
}

impl ReductionState {
    fn from_complex(complex: &FiniteChainComplex) -> Self {
        let active = complex
            .generator_counts
            .iter()
            .map(|(&degree, &count)| (degree, vec![true; count]))
            .collect();
        let matrices = complex
            .differentials
            .iter()
            .map(|(&degree, matrix)| (degree, MutableSparseMatrix::from_sparse(matrix)))
            .collect::<BTreeMap<_, _>>();
        let total_nnz = matrices.values().map(MutableSparseMatrix::nnz).sum();
        Self {
            matrices,
            active,
            total_nnz,
        }
    }
    fn value(&self, degree: i32, row: usize, column: usize) -> Result<BigInt, ChainComplexError> {
        Ok(self
            .matrix(degree)?
            .get(row, column)
            .map_err(|error| ChainComplexError::Certificate(error.to_string()))?)
    }
    fn matrix(&self, degree: i32) -> Result<&MutableSparseMatrix, ChainComplexError> {
        self.matrices
            .get(&degree)
            .ok_or(ChainComplexError::MissingDegree { degree })
    }
    fn matrix_mut(&mut self, degree: i32) -> Result<&mut MutableSparseMatrix, ChainComplexError> {
        self.matrices
            .get_mut(&degree)
            .ok_or(ChainComplexError::MissingDegree { degree })
    }
    fn is_active_pivot(&self, pivot: &UnitPivot) -> bool {
        let Some(row_degree) = pivot.degree.checked_sub(1) else {
            return false;
        };
        self.active
            .get(&pivot.degree)
            .and_then(|v| v.get(pivot.column))
            .copied()
            .unwrap_or(false)
            && self
                .active
                .get(&row_degree)
                .and_then(|v| v.get(pivot.row))
                .copied()
                .unwrap_or(false)
    }
    fn apply_pivot(&mut self, pivot: &UnitPivot) -> Result<PivotChanges, ChainComplexError> {
        if !is_unit(&pivot.value) {
            return Err(ChainComplexError::NonUnitPivot {
                degree: pivot.degree,
                row: pivot.row,
                column: pivot.column,
                value: pivot.value.clone(),
            });
        }
        let affected_degrees = [
            pivot.degree.checked_sub(1),
            Some(pivot.degree),
            pivot.degree.checked_add(1),
        ]
        .into_iter()
        .flatten()
        .filter(|degree| self.matrices.contains_key(degree))
        .collect::<Vec<_>>();
        let before_nnz = affected_degrees
            .iter()
            .map(|degree| (*degree, self.matrices[degree].nnz()))
            .collect::<BTreeMap<_, _>>();
        let row_entries = self
            .matrix(pivot.degree)?
            .row_entries(pivot.row)
            .map_err(to_certificate_error)?
            .filter(|(column, _)| *column != pivot.column && self.active[&pivot.degree][*column])
            .map(|(column, value)| (column, value.clone()))
            .collect::<Vec<_>>();
        let column_entries = self
            .matrix(pivot.degree)?
            .column_rows(pivot.column)
            .map_err(to_certificate_error)?
            .filter(|row| *row != pivot.row && self.active[&(pivot.degree - 1)][*row])
            .map(|row| {
                (
                    row,
                    self.matrix(pivot.degree)
                        .expect("known matrix")
                        .get(row, pivot.column)
                        .expect("valid index"),
                )
            })
            .collect::<Vec<_>>();
        let matrix = self.matrix_mut(pivot.degree)?;
        let mut maximum_new_bits = 0;
        for (row, column_value) in &column_entries {
            for (column, row_value) in &row_entries {
                let old = matrix.get(*row, *column).map_err(to_certificate_error)?;
                // For a = +/-1, a^{-1} = a.
                let new_value = old - column_value * &pivot.value * row_value;
                maximum_new_bits = maximum_new_bits.max(new_value.bits());
                matrix
                    .set(*row, *column, new_value)
                    .map_err(to_certificate_error)?;
            }
        }
        self.matrix_mut(pivot.degree)?
            .clear_row(pivot.row)
            .map_err(to_certificate_error)?;
        self.matrix_mut(pivot.degree)?
            .clear_column(pivot.column)
            .map_err(to_certificate_error)?;
        if let Some(previous) = self.matrices.get_mut(&(pivot.degree - 1)) {
            previous
                .clear_column(pivot.row)
                .map_err(to_certificate_error)?;
        }
        if let Some(next) = self.matrices.get_mut(&(pivot.degree + 1)) {
            next.clear_row(pivot.column).map_err(to_certificate_error)?;
        }
        self.active
            .get_mut(&pivot.degree)
            .ok_or(ChainComplexError::MissingDegree {
                degree: pivot.degree,
            })?[pivot.column] = false;
        self.active
            .get_mut(&(pivot.degree - 1))
            .ok_or(ChainComplexError::MissingDegree {
                degree: pivot.degree - 1,
            })?[pivot.row] = false;
        for degree in affected_degrees {
            self.total_nnz = self.total_nnz - before_nnz[&degree] + self.matrices[&degree].nnz();
        }
        Ok(PivotChanges {
            rows: column_entries.into_iter().map(|(row, _)| row).collect(),
            columns: row_entries.into_iter().map(|(column, _)| column).collect(),
            maximum_new_bits,
        })
    }
    fn nnz(&self) -> usize {
        self.total_nnz
    }
    fn pivot_score(
        &self,
        degree: i32,
        row: usize,
        column: usize,
    ) -> Result<usize, ChainComplexError> {
        let matrix = self.matrix(degree)?;
        let row_degree = matrix
            .row_entries(row)
            .map_err(to_certificate_error)?
            .filter(|(candidate, _)| self.active[&degree][*candidate])
            .count();
        let column_degree = matrix
            .column_rows(column)
            .map_err(to_certificate_error)?
            .filter(|candidate| self.active[&(degree - 1)][*candidate])
            .count();
        Ok(row_degree
            .saturating_sub(1)
            .saturating_mul(column_degree.saturating_sub(1)))
    }
    fn compact(&self) -> Result<FiniteChainComplex, ChainComplexError> {
        let counts = self
            .active
            .iter()
            .map(|(&degree, values)| (degree, values.iter().filter(|active| **active).count()))
            .collect::<BTreeMap<_, _>>();
        let positions = self
            .active
            .iter()
            .map(|(&degree, values)| {
                (
                    degree,
                    values
                        .iter()
                        .enumerate()
                        .filter_map(|(index, active)| active.then_some(index))
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<BTreeMap<_, _>>();
        let inverse = positions
            .iter()
            .map(|(&degree, positions)| {
                let mut indices = BTreeMap::new();
                for (new_index, old_index) in positions.iter().enumerate() {
                    indices.insert(*old_index, new_index);
                }
                (degree, indices)
            })
            .collect::<BTreeMap<_, _>>();
        let mut differentials = BTreeMap::new();
        for (&degree, matrix) in &self.matrices {
            let mut builder = SparseMatrix::builder(counts[&(degree - 1)], counts[&degree]);
            for &old_row in &positions[&(degree - 1)] {
                for (old_column, value) in
                    matrix.row_entries(old_row).map_err(to_certificate_error)?
                {
                    if let Some(&new_column) = inverse[&degree].get(&old_column) {
                        builder
                            .add(inverse[&(degree - 1)][&old_row], new_column, value.clone())
                            .map_err(to_certificate_error)?;
                    }
                }
            }
            differentials.insert(degree, builder.finish());
        }
        FiniteChainComplex::new(counts, differentials)
    }
}

fn candidate_queue(state: &ReductionState) -> BinaryHeap<Reverse<PivotCandidate>> {
    let mut queue = BinaryHeap::new();
    for (&degree, matrix) in &state.matrices {
        for row in 0..matrix.rows() {
            if !state
                .active
                .get(&(degree - 1))
                .and_then(|active| active.get(row))
                .copied()
                .unwrap_or(false)
            {
                continue;
            }
            let row_degree = matrix
                .row_entries(row)
                .expect("valid row")
                .filter(|(column, _)| state.active[&degree][*column])
                .count();
            for (column, value) in matrix.row_entries(row).expect("valid row") {
                if is_unit(value) && state.active[&degree][column] {
                    let column_degree = matrix
                        .column_rows(column)
                        .expect("valid column")
                        .filter(|candidate_row| state.active[&(degree - 1)][*candidate_row])
                        .count();
                    queue.push(Reverse(PivotCandidate {
                        score: row_degree
                            .saturating_sub(1)
                            .saturating_mul(column_degree.saturating_sub(1)),
                        degree,
                        row,
                        column,
                    }));
                }
            }
        }
    }
    queue
}

fn enqueue_changed_candidates(
    state: &ReductionState,
    queue: &mut BinaryHeap<Reverse<PivotCandidate>>,
    degree: i32,
    rows: &[usize],
    columns: &[usize],
) {
    let Ok(matrix) = state.matrix(degree) else {
        return;
    };
    for &row in rows {
        if !state
            .active
            .get(&(degree - 1))
            .and_then(|active| active.get(row))
            .copied()
            .unwrap_or(false)
        {
            continue;
        }
        let row_degree = matrix
            .row_entries(row)
            .expect("valid changed row")
            .filter(|(column, _)| state.active[&degree][*column])
            .count();
        for &column in columns {
            if !state.active[&degree].get(column).copied().unwrap_or(false) {
                continue;
            }
            let value = matrix.get(row, column).expect("valid changed coordinate");
            if is_unit(&value) {
                let column_degree = matrix
                    .column_rows(column)
                    .expect("valid changed column")
                    .filter(|candidate_row| state.active[&(degree - 1)][*candidate_row])
                    .count();
                queue.push(Reverse(PivotCandidate {
                    score: row_degree
                        .saturating_sub(1)
                        .saturating_mul(column_degree.saturating_sub(1)),
                    degree,
                    row,
                    column,
                }));
            }
        }
    }
}

fn is_unit(value: &BigInt) -> bool {
    value.abs() == BigInt::one()
}
fn to_certificate_error(error: crate::SparseMatrixError) -> ChainComplexError {
    ChainComplexError::Certificate(error.to_string())
}

fn complex_resource_usage(complex: &FiniteChainComplex) -> (usize, u64, usize) {
    let mut nnz = 0usize;
    let mut bits = 0u64;
    let mut shape_slots = 0usize;
    for &count in complex.generator_counts.values() {
        shape_slots = shape_slots.saturating_add(count);
    }
    for matrix in complex.differentials.values() {
        nnz = nnz.saturating_add(matrix.nnz());
        shape_slots = shape_slots
            .saturating_add(matrix.rows())
            .saturating_add(matrix.columns());
        bits = bits.max(
            matrix
                .rows_iter()
                .flatten()
                .map(|(_, value)| value.bits())
                .max()
                .unwrap_or(0),
        );
    }
    (nnz, bits, shape_slots)
}

fn check_reduction_resources(
    nnz: usize,
    bits: u64,
    shape_slots: usize,
    options: &UnitReductionOptions,
) -> Result<(), ChainComplexError> {
    if shape_slots > options.max_shape_slots {
        return Err(ChainComplexError::ReductionLimit {
            kind: "matrix shape slots",
            observed: shape_slots,
            budget: options.max_shape_slots,
        });
    }
    if nnz > options.max_nnz {
        return Err(ChainComplexError::ReductionLimit {
            kind: "NNZ",
            observed: nnz,
            budget: options.max_nnz,
        });
    }
    if bits > options.max_entry_bits {
        return Err(ChainComplexError::ReductionLimit {
            kind: "coefficient bit length",
            observed: bits as usize,
            budget: options.max_entry_bits as usize,
        });
    }
    Ok(())
}

/// Deterministic Miller--Rabin for `u64`; kept here so combinatoric-core
/// remains independent of standalone Polytool packaging.
pub fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n == p {
            return true;
        }
        if n.is_multiple_of(p) {
            return false;
        }
    }
    let mut d = n - 1;
    let mut s = 0;
    while d.is_multiple_of(2) {
        d /= 2;
        s += 1;
    }
    for base in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = pow_mod_u64(base, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        let mut passed = false;
        for _ in 1..s {
            x = mul_mod_u64(x, x, n);
            if x == n - 1 {
                passed = true;
                break;
            }
        }
        if !passed {
            return false;
        }
    }
    true
}

fn mul_mod_u64(left: u64, right: u64, modulus: u64) -> u64 {
    ((left as u128 * right as u128) % modulus as u128) as u64
}
fn pow_mod_u64(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut result = 1;
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = mul_mod_u64(result, base, modulus);
        }
        base = mul_mod_u64(base, base, modulus);
        exponent >>= 1;
    }
    result
}

fn complex_sha256(complex: &FiniteChainComplex) -> String {
    let mut hasher = Sha256::new();
    hasher.update(b"combinatoric-core/unit-cancellation/input/v2");
    hasher.update((complex.generator_counts.len() as u64).to_le_bytes());
    for (&degree, &count) in &complex.generator_counts {
        hasher.update(b"G");
        hasher.update(degree.to_le_bytes());
        hasher.update((count as u64).to_le_bytes());
    }
    hasher.update((complex.differentials.len() as u64).to_le_bytes());
    for (&degree, matrix) in &complex.differentials {
        hasher.update(b"D");
        hasher.update(degree.to_le_bytes());
        hasher.update((matrix.rows() as u64).to_le_bytes());
        hasher.update((matrix.columns() as u64).to_le_bytes());
        hasher.update((matrix.nnz() as u64).to_le_bytes());
        for row in 0..matrix.rows() {
            for (column, value) in matrix.row(row).expect("valid row") {
                let bytes = value.to_signed_bytes_le();
                hasher.update((row as u64).to_le_bytes());
                hasher.update((column as u64).to_le_bytes());
                hasher.update((bytes.len() as u64).to_le_bytes());
                hasher.update(bytes);
            }
        }
    }
    format!("{:x}", hasher.finalize())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn complex(counts: &[(i32, usize)], maps: &[(i32, &[&[i64]])]) -> FiniteChainComplex {
        let counts = counts.iter().copied().collect();
        let differentials = maps
            .iter()
            .map(|(degree, rows)| {
                let dense = rows
                    .iter()
                    .map(|row| {
                        row.iter()
                            .map(|value| BigInt::from(*value))
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>();
                (*degree, SparseMatrix::from_dense(&dense).unwrap())
            })
            .collect();
        FiniteChainComplex::new(counts, differentials).unwrap()
    }

    #[test]
    fn cancellation_replays_and_rejects_tampering() {
        let input = complex(&[(0, 2), (1, 2)], &[(1, &[&[1, 1], &[0, 1]])]);
        let options = UnitReductionOptions {
            record_certificate: true,
            ..UnitReductionOptions::default()
        };
        let result = cancel_units(&input, options.clone()).unwrap();
        assert_eq!(result.reduced.count(0), 0);
        assert_eq!(result.reduced.count(1), 0);
        assert_eq!(
            replay_unit_cancellation(&input, result.certificate.as_ref().unwrap(), &options)
                .unwrap(),
            result.reduced
        );
        let wrong_residual =
            FiniteChainComplex::new([(0, 1), (1, 0)].into_iter().collect(), BTreeMap::new())
                .unwrap();
        assert!(matches!(
            replay_unit_cancellation_and_verify(
                &input,
                result.certificate.as_ref().unwrap(),
                &wrong_residual,
                &options,
            ),
            Err(ChainComplexError::Certificate(_))
        ));
        let mut tampered = result.certificate.unwrap();
        tampered.pivots[0].value = BigInt::from(2);
        assert!(matches!(
            replay_unit_cancellation(&input, &tampered, &options),
            Err(ChainComplexError::Certificate(_))
        ));
    }

    #[test]
    fn integral_homology_detects_torsion_and_free_groups() {
        let torsion = complex(&[(0, 1), (1, 1)], &[(1, &[&[6]])]);
        let groups = integral_homology(&torsion, SmithOptions::default()).unwrap();
        assert_eq!(
            groups[&0],
            AbelianGroup {
                free_rank: 0,
                torsion_invariants: vec![6.into()]
            }
        );
        let free = complex(&[(0, 1), (1, 1)], &[(1, &[&[0]])]);
        assert_eq!(
            integral_homology(&free, SmithOptions::default()).unwrap()[&1].free_rank,
            1
        );
        assert_eq!(
            universal_coefficient_dimension(&groups[&0], None, 2).unwrap(),
            1
        );
        assert_eq!(
            universal_coefficient_dimension(&groups[&0], None, 5).unwrap(),
            0
        );
    }

    #[test]
    fn metadata_only_huge_free_complex_never_materializes_missing_zero_maps() {
        let metadata_complex = FiniteChainComplex::new(
            [(0, 1_000_000_000_000usize), (1, 0)].into_iter().collect(),
            BTreeMap::new(),
        )
        .unwrap();
        let groups = integral_homology(&metadata_complex, SmithOptions::default()).unwrap();
        assert_eq!(groups[&0].free_rank, 1_000_000_000_000);
        assert_eq!(groups[&1].free_rank, 0);
        assert!(matches!(
            metadata_complex.differential_or_zero_with_limits(
                1,
                SparseMatrixLimits {
                    max_shape_slots: 8,
                    max_nnz: 0
                }
            ),
            Err(ChainComplexError::ReductionLimit { .. })
        ));
        assert_eq!(
            metadata_complex.differential_or_zero(i32::MIN).shape(),
            (0, 0)
        );
        assert!(matches!(
            metadata_complex
                .differential_or_zero_with_limits(i32::MIN, SparseMatrixLimits::default()),
            Err(ChainComplexError::UnsupportedDegree { .. })
        ));

        // The stored D_2 is tiny, while validation must recognize absent D_1
        // as zero without allocating its 10^12-by-1 shaped representation.
        let mut differentials = BTreeMap::new();
        differentials.insert(2, SparseMatrix::zero(1, 1));
        let validation_only = FiniteChainComplex::new(
            [(0, 1_000_000_000_000usize), (1, 1), (2, 1)]
                .into_iter()
                .collect(),
            differentials,
        );
        assert!(validation_only.is_ok());

        let stored_nonzero = complex(&[(0, 1), (1, 1)], &[(1, &[&[1]])]);
        assert!(matches!(
            stored_nonzero.differential_or_zero_with_limits(
                1,
                SparseMatrixLimits {
                    max_shape_slots: 2,
                    max_nnz: 0
                }
            ),
            Err(ChainComplexError::ReductionLimit { kind: "NNZ", .. })
        ));
        let stored_tall = FiniteChainComplex::new(
            [(0, 9), (1, 0)].into_iter().collect(),
            [(1, SparseMatrix::zero(9, 0))].into_iter().collect(),
        )
        .unwrap();
        assert!(matches!(
            stored_tall.differential_or_zero_with_limits(
                1,
                SparseMatrixLimits {
                    max_shape_slots: 8,
                    max_nnz: 0
                }
            ),
            Err(ChainComplexError::ReductionLimit {
                kind: "matrix shape slots",
                ..
            })
        ));
    }

    #[test]
    fn invalid_d_squared_is_witnessed() {
        let counts = [(0, 1), (1, 1), (2, 1)].into_iter().collect();
        let differentials = [
            (1, SparseMatrix::from_dense(&vec![vec![1.into()]]).unwrap()),
            (2, SparseMatrix::from_dense(&vec![vec![1.into()]]).unwrap()),
        ]
        .into_iter()
        .collect();
        assert!(matches!(
            FiniteChainComplex::new(counts, differentials),
            Err(ChainComplexError::DifferentialDoesNotSquare { .. })
        ));
    }

    #[test]
    fn mixed_torsion_and_previous_degree_uct_are_exact() {
        let input = complex(
            &[(0, 1), (1, 2), (2, 1)],
            &[(1, &[&[2, 0]]), (2, &[&[0], &[3]])],
        );
        let groups = integral_homology(&input, SmithOptions::default()).unwrap();
        assert_eq!(
            groups[&0],
            AbelianGroup {
                free_rank: 0,
                torsion_invariants: vec![2.into()]
            }
        );
        assert_eq!(
            groups[&1],
            AbelianGroup {
                free_rank: 0,
                torsion_invariants: vec![3.into()]
            }
        );
        assert_eq!(
            universal_coefficient_dimension(&groups[&1], Some(&groups[&0]), 2).unwrap(),
            1
        );
        assert_eq!(
            universal_coefficient_dimension(&groups[&1], Some(&groups[&0]), 3).unwrap(),
            1
        );
        for modulus in [0, 1, 4] {
            assert!(
                universal_coefficient_dimension(&groups[&1], Some(&groups[&0]), modulus).is_err()
            );
        }
    }

    #[test]
    fn cancellation_random_small_complexes_preserve_exact_groups_and_limits() {
        let mut seed = 0xa11ce_u64;
        for _ in 0..20 {
            let mut rows = vec![vec![0i64; 3]; 3];
            for row in &mut rows {
                for value in row {
                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                    *value = match (seed >> 32) % 5 {
                        0 => -1,
                        1 => 1,
                        _ => 0,
                    };
                }
            }
            let refs = rows.iter().map(Vec::as_slice).collect::<Vec<_>>();
            let input = complex(&[(0, 3), (1, 3)], &[(1, refs.as_slice())]);
            let before = integral_homology(&input, SmithOptions::default()).unwrap();
            let options = UnitReductionOptions {
                record_certificate: true,
                ..UnitReductionOptions::default()
            };
            let reduced = cancel_units(&input, options.clone()).unwrap();
            assert_eq!(
                integral_homology(&reduced.reduced, SmithOptions::default()).unwrap(),
                before
            );
            replay_unit_cancellation_and_verify(
                &input,
                reduced.certificate.as_ref().unwrap(),
                &reduced.reduced,
                &options,
            )
            .unwrap();
        }
        let diagonal = complex(&[(0, 1), (1, 1)], &[(1, &[&[1 << 20]])]);
        assert!(matches!(
            cancel_units(
                &diagonal,
                UnitReductionOptions {
                    max_entry_bits: 8,
                    ..UnitReductionOptions::default()
                }
            ),
            Err(ChainComplexError::ReductionLimit { .. })
        ));
        assert!(matches!(
            FiniteChainComplex::new([(i32::MIN, 0)].into_iter().collect(), BTreeMap::new()),
            Err(ChainComplexError::UnsupportedDegree { .. })
        ));
    }

    #[test]
    fn deterministic_u64_primality_rejects_known_strong_pseudoprime() {
        assert!(!is_prime_u64(341_550_071_728_321));
        assert!(is_prime_u64(2_305_843_009_213_693_951));
    }
}
