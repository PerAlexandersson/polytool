//! Principal-stratum chain complexes for the Shapiro resonance computation.
//!
//! This module contains the composition-specific cell grammar only.  Sparse
//! storage, cancellation, Smith forms, and rank engines remain in shared
//! crates.  Reported manuscript groups are calibration targets, never inputs.

use combinatoric_core::{field_betti_number, FiniteChainComplex, SparseMatrixBuilder, SparseMatrixError};
use num_bigint::BigInt;
use num::Integer;
use num_rational::Ratio;
use num_traits::ToPrimitive;
use polytool::{sparse_modular_linear_system_consistency_with_options, SparseModEliminationOptions, SparseModEliminationStats, SparseModRow, SparseModRowOrder};
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::time::Duration;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum PrincipalStratumError {
    EmptyComposition,
    ZeroPart { position: usize },
    WeightTooSmall { d: u32, weight: u32 },
    ParityMismatch { d: u32, weight: u32 },
    WeightOverflow,
    MissingBoundaryTarget { cell: Vec<u32> },
    UnexpectedBoundaryDegree { source: i32, target: i32, cell: Vec<u32> },
    Sparse(SparseMatrixError),
    Modular(String),
    DenseBudget(SparseMatrixError),
}

impl std::fmt::Display for PrincipalStratumError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyComposition => write!(f, "omega must be nonempty"),
            Self::ZeroPart { position } => write!(f, "omega has zero part at position {position}"),
            Self::WeightTooSmall { d, weight } => write!(f, "d={d} is smaller than omega's weight {weight}"),
            Self::ParityMismatch { d, weight } => write!(f, "d={d} and omega weight {weight} have different parity"),
            Self::WeightOverflow => write!(f, "composition weight overflow"),
            Self::MissingBoundaryTarget { cell } => write!(f, "boundary target of {cell:?} was not generated"),
            Self::UnexpectedBoundaryDegree { source, target, cell } => write!(f, "boundary of {cell:?} has degree {target}, expected {}", source - 1),
            Self::Sparse(error) | Self::DenseBudget(error) => write!(f, "sparse matrix error: {error}"),
            Self::Modular(error) => write!(f, "modular rank error: {error}"),
        }
    }
}
impl std::error::Error for PrincipalStratumError {}
impl From<SparseMatrixError> for PrincipalStratumError { fn from(value: SparseMatrixError) -> Self { Self::Sparse(value) } }

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct PrincipalStratumCell {
    parts: Vec<u32>,
}
impl PrincipalStratumCell {
    pub fn parts(&self) -> &[u32] { &self.parts }
    pub fn weight(&self) -> u32 { self.parts.iter().try_fold(0u32, |sum, part| sum.checked_add(*part)).expect("validated cell weight") }
    pub fn degree(&self, d: u32) -> i32 {
        let weight = i32::try_from(self.weight()).expect("u32 fits i32 for supported model");
        i32::try_from(d).expect("u32 fits i32 for supported model") - weight + i32::try_from(self.parts.len()).expect("length fits i32")
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PrincipalStratumModel {
    omega: Vec<u32>,
    d: u32,
    cells_by_degree: BTreeMap<i32, Vec<PrincipalStratumCell>>,
    indices: BTreeMap<Vec<u32>, (i32, usize)>,
    complex: FiniteChainComplex,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PrincipalStratumBuildStats {
    pub enumeration_time: Duration,
    pub assembly_and_validation_time: Duration,
}

impl PrincipalStratumModel {
    pub fn build(omega: Vec<u32>, d: u32) -> Result<Self, PrincipalStratumError> {
        Self::build_bounded(omega, d, usize::MAX).map(|(model, _)| model)
    }

    pub fn build_bounded(omega: Vec<u32>, d: u32, max_cells: usize) -> Result<(Self, PrincipalStratumBuildStats), PrincipalStratumError> {
        validate_input(&omega, d)?;
        let started = std::time::Instant::now();
        let cells = enumerate_cells_bfs(&omega, d, max_cells)?;
        let enumeration_time = started.elapsed();
        let started = std::time::Instant::now();
        let model = Self::from_cells(omega, d, cells)?;
        Ok((model, PrincipalStratumBuildStats { enumeration_time, assembly_and_validation_time: started.elapsed() }))
    }

    /// Independent depth-first closure enumerator used for small-oracle
    /// agreement tests.  It has deliberately separate traversal mechanics
    /// from the breadth-first production builder.
    pub fn build_with_dfs_oracle(omega: Vec<u32>, d: u32) -> Result<Self, PrincipalStratumError> {
        validate_input(&omega, d)?;
        let cells = enumerate_cells_dfs(&omega, d, usize::MAX)?;
        Self::from_cells(omega, d, cells)
    }

    fn from_cells(omega: Vec<u32>, d: u32, cells: BTreeSet<Vec<u32>>) -> Result<Self, PrincipalStratumError> {
        let mut cells_by_degree = BTreeMap::<i32, Vec<PrincipalStratumCell>>::new();
        for parts in cells { cells_by_degree.entry(cell_degree(&parts, d)?).or_default().push(PrincipalStratumCell { parts }); }
        for cells in cells_by_degree.values_mut() { cells.sort(); }
        let indices = cells_by_degree.iter().flat_map(|(&degree, cells)| cells.iter().enumerate().map(move |(index, cell)| (cell.parts.clone(), (degree, index)))).collect::<BTreeMap<_, _>>();
        let generator_counts = cells_by_degree.iter().map(|(&degree, cells)| (degree, cells.len())).collect();
        let mut differentials = BTreeMap::new();
        for (&degree, cells) in &cells_by_degree {
            if !cells_by_degree.contains_key(&(degree - 1)) { continue; }
            let mut boundary = SparseMatrixBuilder::new(cells_by_degree[&(degree - 1)].len(), cells.len());
            for (column, cell) in cells.iter().enumerate() {
                for (target, sign) in signed_direct_moves(&cell.parts, d)? {
                    let Some(&(target_degree, row)) = indices.get(&target) else { return Err(PrincipalStratumError::MissingBoundaryTarget { cell: target }); };
                    if target_degree != degree - 1 { return Err(PrincipalStratumError::UnexpectedBoundaryDegree { source: degree, target: target_degree, cell: cell.parts.clone() }); }
                    boundary.add(row, column, BigInt::from(sign))?;
                }
            }
            differentials.insert(degree, boundary.finish());
        }
        let complex = FiniteChainComplex::new(generator_counts, differentials).map_err(|error| PrincipalStratumError::Modular(error.to_string()))?;
        Ok(Self { omega, d, cells_by_degree, indices, complex })
    }

    pub fn omega(&self) -> &[u32] { &self.omega }
    pub fn d(&self) -> u32 { self.d }
    pub fn cells_by_degree(&self) -> &BTreeMap<i32, Vec<PrincipalStratumCell>> { &self.cells_by_degree }
    pub fn cell_count(&self) -> usize { self.indices.len() }
    pub fn complex(&self) -> &FiniteChainComplex { &self.complex }
    pub fn euler_characteristic(&self) -> i64 {
        self.cells_by_degree.iter().map(|(&degree, cells)| if degree.rem_euclid(2) == 0 { cells.len() as i64 } else { -(cells.len() as i64) }).sum()
    }
    pub fn initial_nnz(&self) -> usize { self.complex.degrees().map(|degree| self.complex.differential_or_zero(degree).nnz()).sum() }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ModularFieldBetti {
    pub prime: u64,
    pub ranks: BTreeMap<i32, usize>,
    pub betti_numbers: BTreeMap<i32, usize>,
    pub solver_stats: BTreeMap<i32, SparseModEliminationStats>,
    /// Exact dimensions over this one finite field, not integral evidence.
    pub evidence_label: &'static str,
}

pub fn modular_field_betti(
    model: &PrincipalStratumModel,
    prime: u64,
    row_order: SparseModRowOrder,
) -> Result<ModularFieldBetti, PrincipalStratumError> {
    let mut ranks = BTreeMap::new();
    let mut solver_stats = BTreeMap::new();
    for degree in model.complex.degrees() {
        let matrix = model.complex.differential_or_zero(degree);
        let rows = matrix.transpose().rows_iter().map(|entries| {
            let entries = entries.map(|(column, value)| bigint_mod_u64(value, prime).map(|value| (column, value))).collect::<Result<Vec<_>, _>>()?;
            SparseModRow::from_entries(entries, 0, prime).map_err(|error| PrincipalStratumError::Modular(error.to_string()))
        }).collect::<Result<Vec<_>, _>>()?;
        let result = sparse_modular_linear_system_consistency_with_options(rows, matrix.rows(), prime, SparseModEliminationOptions { row_order, compute_solution: false }).map_err(|error| PrincipalStratumError::Modular(error.to_string()))?;
        if !result.consistent { return Err(PrincipalStratumError::Modular("zero-RHS boundary matrix was unexpectedly inconsistent".into())); }
        ranks.insert(degree, result.rank);
        solver_stats.insert(degree, result.stats);
    }
    let betti_numbers = model.complex.degrees().map(|degree| {
        let count = model.complex.count(degree);
        let betti = field_betti_number(count, ranks.get(&degree).copied().unwrap_or(0), ranks.get(&(degree + 1)).copied().unwrap_or(0)).expect("ranks of a chain complex cannot exceed its count");
        (degree, betti)
    }).collect();
    Ok(ModularFieldBetti { prime, ranks, betti_numbers, solver_stats, evidence_label: "exact F_p dimensions" })
}

/// Small-matrix rational oracle.  It explicitly declines dense conversion
/// beyond `dense_budget`; callers must not confuse it with a large sparse rank.
pub fn rational_field_betti(
    model: &PrincipalStratumModel,
    dense_budget: usize,
) -> Result<BTreeMap<i32, usize>, PrincipalStratumError> {
    let mut ranks = BTreeMap::new();
    for degree in model.complex.degrees() {
        let matrix = model.complex.differential_or_zero(degree);
        let dense = matrix.to_dense_with_budget(dense_budget).map_err(PrincipalStratumError::DenseBudget)?;
        let rank = if dense.is_empty() { 0 } else {
            sym_poly_core::linear_algebra::rank(&dense.into_iter().map(|row| row.into_iter().map(Ratio::from_integer).collect()).collect::<Vec<Vec<Ratio<BigInt>>>>())
        };
        ranks.insert(degree, rank);
    }
    Ok(model.complex.degrees().map(|degree| {
        let betti = field_betti_number(model.complex.count(degree), ranks.get(&degree).copied().unwrap_or(0), ranks.get(&(degree + 1)).copied().unwrap_or(0)).expect("chain ranks are bounded");
        (degree, betti)
    }).collect())
}

fn validate_input(omega: &[u32], d: u32) -> Result<(), PrincipalStratumError> {
    if omega.is_empty() { return Err(PrincipalStratumError::EmptyComposition); }
    if let Some((position, _)) = omega.iter().enumerate().find(|(_, part)| **part == 0) { return Err(PrincipalStratumError::ZeroPart { position }); }
    let weight = cell_weight(omega)?;
    if d < weight { return Err(PrincipalStratumError::WeightTooSmall { d, weight }); }
    if (d - weight) % 2 != 0 { return Err(PrincipalStratumError::ParityMismatch { d, weight }); }
    Ok(())
}

fn enumerate_cells_bfs(omega: &[u32], d: u32, max_cells: usize) -> Result<BTreeSet<Vec<u32>>, PrincipalStratumError> {
    let mut seen = BTreeSet::new();
    let mut queue = VecDeque::new();
    seen.insert(omega.to_vec()); queue.push_back(omega.to_vec());
    while let Some(cell) = queue.pop_front() {
        for target in unsigned_direct_moves(&cell, d)? {
            if seen.insert(target.clone()) {
                if seen.len() > max_cells { return Err(PrincipalStratumError::Modular(format!("cell limit {max_cells} exceeded"))); }
                queue.push_back(target);
            }
        }
    }
    Ok(seen)
}

fn enumerate_cells_dfs(omega: &[u32], d: u32, max_cells: usize) -> Result<BTreeSet<Vec<u32>>, PrincipalStratumError> {
    fn visit(cell: Vec<u32>, d: u32, max_cells: usize, seen: &mut BTreeSet<Vec<u32>>) -> Result<(), PrincipalStratumError> {
        if !seen.insert(cell.clone()) { return Ok(()); }
        if seen.len() > max_cells { return Err(PrincipalStratumError::Modular(format!("cell limit {max_cells} exceeded"))); }
        for target in unsigned_direct_moves(&cell, d)? { visit(target, d, max_cells, seen)?; }
        Ok(())
    }
    let mut seen = BTreeSet::new();
    visit(omega.to_vec(), d, max_cells, &mut seen)?;
    Ok(seen)
}

fn unsigned_direct_moves(cell: &[u32], d: u32) -> Result<Vec<Vec<u32>>, PrincipalStratumError> {
    let mut targets = Vec::with_capacity(cell.len() * 2 + 1);
    for position in 0..cell.len().saturating_sub(1) {
        let mut target = cell.to_vec();
        target[position] = target[position].checked_add(target[position + 1]).ok_or(PrincipalStratumError::WeightOverflow)?;
        target.remove(position + 1);
        targets.push(target);
    }
    if cell_weight(cell)?.checked_add(2).ok_or(PrincipalStratumError::WeightOverflow)? <= d {
        for position in 0..=cell.len() {
            let mut target = cell.to_vec(); target.insert(position, 2); targets.push(target);
        }
    }
    Ok(targets)
}

fn signed_direct_moves(cell: &[u32], d: u32) -> Result<Vec<(Vec<u32>, i8)>, PrincipalStratumError> {
    let mut targets = Vec::with_capacity(cell.len() * 2 + 1);
    for position in 0..cell.len().saturating_sub(1) {
        let mut target = cell.to_vec(); target[position] = target[position].checked_add(target[position + 1]).ok_or(PrincipalStratumError::WeightOverflow)?; target.remove(position + 1);
        targets.push((target, if position % 2 == 0 { 1 } else { -1 }));
    }
    if cell_weight(cell)?.checked_add(2).ok_or(PrincipalStratumError::WeightOverflow)? <= d {
        for position in 0..=cell.len() {
            let mut target = cell.to_vec(); target.insert(position, 2);
            targets.push((target, if position % 2 == 0 { 1 } else { -1 }));
        }
    }
    Ok(targets)
}

fn cell_weight(cell: &[u32]) -> Result<u32, PrincipalStratumError> { cell.iter().try_fold(0u32, |sum, part| sum.checked_add(*part).ok_or(PrincipalStratumError::WeightOverflow)) }
fn cell_degree(cell: &[u32], d: u32) -> Result<i32, PrincipalStratumError> {
    let weight = i32::try_from(cell_weight(cell)?).map_err(|_| PrincipalStratumError::WeightOverflow)?;
    Ok(i32::try_from(d).map_err(|_| PrincipalStratumError::WeightOverflow)? - weight + i32::try_from(cell.len()).map_err(|_| PrincipalStratumError::WeightOverflow)?)
}
fn bigint_mod_u64(value: &BigInt, prime: u64) -> Result<u64, PrincipalStratumError> {
    value.mod_floor(&BigInt::from(prime)).to_u64().ok_or_else(|| PrincipalStratumError::Modular("reduced coefficient did not fit u64".into()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_calibration_has_expected_cells_chain_and_field_groups() {
        let model = PrincipalStratumModel::build(vec![1, 1], 4).unwrap();
        assert_eq!(model.cells_by_degree.iter().map(|(&degree, cells)| (degree, cells.len())).collect::<BTreeMap<_, _>>(), [(1, 1), (2, 3), (3, 4), (4, 1)].into_iter().collect());
        assert_eq!(model.cell_count(), 9);
        let modular = modular_field_betti(&model, 101, SparseModRowOrder::Input).unwrap();
        assert_eq!(modular.betti_numbers, [(1, 0), (2, 0), (3, 1), (4, 0)].into_iter().collect());
        assert_eq!(rational_field_betti(&model, 100).unwrap(), modular.betti_numbers);
    }

    #[test]
    fn bfs_and_dfs_oracles_agree_and_invalid_inputs_are_rejected() {
        let bfs = PrincipalStratumModel::build(vec![1, 1, 1], 5).unwrap();
        let dfs = PrincipalStratumModel::build_with_dfs_oracle(vec![1, 1, 1], 5).unwrap();
        assert_eq!(bfs.cells_by_degree, dfs.cells_by_degree);
        assert!(matches!(PrincipalStratumModel::build(vec![1, 0, 1], 4), Err(PrincipalStratumError::ZeroPart { .. })));
        assert!(matches!(PrincipalStratumModel::build(vec![1], 2), Err(PrincipalStratumError::ParityMismatch { .. })));
    }
}
