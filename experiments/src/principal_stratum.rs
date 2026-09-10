//! Principal-stratum chain complexes for the Shapiro resonance computation.
//!
//! This module contains the composition-specific cell grammar only.  Sparse
//! storage, cancellation, Smith forms, and rank engines remain in shared
//! crates.  Reported manuscript groups are calibration targets, never inputs.

use combinatoric_core::{
    field_betti_number, FiniteChainComplex, SparseMatrixBuilder, SparseMatrixError,
};
use num_bigint::BigInt;
use num_rational::Ratio;
use polytool::{
    sparse_modular_linear_system_consistency_with_options, SparseModEliminationOptions,
    SparseModEliminationStats, SparseModRow, SparseModRowOrder,
};
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::time::Duration;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum PrincipalStratumError {
    EmptyComposition,
    ZeroPart {
        position: usize,
    },
    WeightTooSmall {
        d: u32,
        weight: u32,
    },
    ParityMismatch {
        d: u32,
        weight: u32,
    },
    WeightOverflow,
    GradeOutOfRange {
        d: u32,
    },
    ResourceLimit {
        kind: &'static str,
        observed: usize,
        budget: usize,
    },
    InvalidModulus {
        modulus: u64,
    },
    MissingBoundaryTarget {
        cell: Vec<u32>,
    },
    UnexpectedBoundaryDegree {
        source: i32,
        target: i32,
        cell: Vec<u32>,
    },
    Sparse(SparseMatrixError),
    Modular(String),
    DenseBudget(SparseMatrixError),
}

impl std::fmt::Display for PrincipalStratumError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyComposition => write!(f, "omega must be nonempty"),
            Self::ZeroPart { position } => write!(f, "omega has zero part at position {position}"),
            Self::WeightTooSmall { d, weight } => {
                write!(f, "d={d} is smaller than omega's weight {weight}")
            }
            Self::ParityMismatch { d, weight } => {
                write!(f, "d={d} and omega weight {weight} have different parity")
            }
            Self::WeightOverflow => write!(f, "composition weight overflow"),
            Self::GradeOutOfRange { d } => {
                write!(f, "d={d} cannot be represented by the signed chain grading")
            }
            Self::ResourceLimit {
                kind,
                observed,
                budget,
            } => write!(f, "{kind} limit {budget} exceeded by {observed}"),
            Self::InvalidModulus { modulus } => write!(f, "modulus {modulus} is not prime"),
            Self::MissingBoundaryTarget { cell } => {
                write!(f, "boundary target of {cell:?} was not generated")
            }
            Self::UnexpectedBoundaryDegree {
                source,
                target,
                cell,
            } => write!(
                f,
                "boundary of {cell:?} has degree {target}, expected {}",
                source - 1
            ),
            Self::Sparse(error) | Self::DenseBudget(error) => {
                write!(f, "sparse matrix error: {error}")
            }
            Self::Modular(error) => write!(f, "modular rank error: {error}"),
        }
    }
}
impl std::error::Error for PrincipalStratumError {}
impl From<SparseMatrixError> for PrincipalStratumError {
    fn from(value: SparseMatrixError) -> Self {
        Self::Sparse(value)
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct PrincipalStratumCell {
    parts: Vec<u32>,
}
impl PrincipalStratumCell {
    pub fn parts(&self) -> &[u32] {
        &self.parts
    }
    pub fn weight(&self) -> u32 {
        self.parts
            .iter()
            .try_fold(0u32, |sum, part| sum.checked_add(*part))
            .expect("validated cell weight")
    }
    pub fn degree(&self, d: u32) -> i32 {
        let weight = i32::try_from(self.weight()).expect("u32 fits i32 for supported model");
        i32::try_from(d).expect("u32 fits i32 for supported model") - weight
            + i32::try_from(self.parts.len()).expect("length fits i32")
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

/// Construction limits are enforced while cells and boundary terms are being
/// generated; they are not a post-hoc check after a large complex exists.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PrincipalStratumBuildLimits {
    pub max_cells: usize,
    pub max_boundary_nnz: usize,
}

impl Default for PrincipalStratumBuildLimits {
    fn default() -> Self {
        Self {
            max_cells: usize::MAX,
            max_boundary_nnz: usize::MAX,
        }
    }
}

impl PrincipalStratumModel {
    pub fn build(omega: Vec<u32>, d: u32) -> Result<Self, PrincipalStratumError> {
        Self::build_bounded(omega, d, usize::MAX).map(|(model, _)| model)
    }

    pub fn build_bounded(
        omega: Vec<u32>,
        d: u32,
        max_cells: usize,
    ) -> Result<(Self, PrincipalStratumBuildStats), PrincipalStratumError> {
        Self::build_with_limits(
            omega,
            d,
            PrincipalStratumBuildLimits {
                max_cells,
                ..PrincipalStratumBuildLimits::default()
            },
        )
    }

    pub fn build_with_limits(
        omega: Vec<u32>,
        d: u32,
        limits: PrincipalStratumBuildLimits,
    ) -> Result<(Self, PrincipalStratumBuildStats), PrincipalStratumError> {
        validate_input(&omega, d)?;
        let started = std::time::Instant::now();
        let cells = enumerate_cells_automaton(&omega, d, limits.max_cells)?;
        let enumeration_time = started.elapsed();
        let started = std::time::Instant::now();
        let model = Self::from_cells(omega, d, cells, limits.max_boundary_nnz)?;
        Ok((
            model,
            PrincipalStratumBuildStats {
                enumeration_time,
                assembly_and_validation_time: started.elapsed(),
            },
        ))
    }

    /// Small reference closure oracle.  Production generation uses the
    /// membership-subset automaton; this intentionally retains the direct
    /// move closure only as an independent test oracle.
    pub fn build_with_bfs_oracle(omega: Vec<u32>, d: u32) -> Result<Self, PrincipalStratumError> {
        validate_input(&omega, d)?;
        let cells = enumerate_cells_bfs(&omega, d, usize::MAX)?;
        Self::from_cells(omega, d, cells, usize::MAX)
    }

    fn from_cells(
        omega: Vec<u32>,
        d: u32,
        cells: BTreeSet<Vec<u32>>,
        max_boundary_nnz: usize,
    ) -> Result<Self, PrincipalStratumError> {
        let mut cells_by_degree = BTreeMap::<i32, Vec<PrincipalStratumCell>>::new();
        for parts in cells {
            cells_by_degree
                .entry(cell_degree(&parts, d)?)
                .or_default()
                .push(PrincipalStratumCell { parts });
        }
        for cells in cells_by_degree.values_mut() {
            cells.sort();
        }
        let indices = cells_by_degree
            .iter()
            .flat_map(|(&degree, cells)| {
                cells
                    .iter()
                    .enumerate()
                    .map(move |(index, cell)| (cell.parts.clone(), (degree, index)))
            })
            .collect::<BTreeMap<_, _>>();
        let generator_counts = cells_by_degree
            .iter()
            .map(|(&degree, cells)| (degree, cells.len()))
            .collect();
        let mut differentials = BTreeMap::new();
        let mut generated_boundary_terms = 0usize;
        for (&degree, cells) in &cells_by_degree {
            if !cells_by_degree.contains_key(&(degree - 1)) {
                continue;
            }
            let mut boundary =
                SparseMatrixBuilder::new(cells_by_degree[&(degree - 1)].len(), cells.len());
            for (column, cell) in cells.iter().enumerate() {
                for (target, sign) in signed_direct_moves(&cell.parts, d)? {
                    generated_boundary_terms = generated_boundary_terms.checked_add(1).ok_or(
                        PrincipalStratumError::ResourceLimit {
                            kind: "boundary term",
                            observed: usize::MAX,
                            budget: max_boundary_nnz,
                        },
                    )?;
                    if generated_boundary_terms > max_boundary_nnz {
                        return Err(PrincipalStratumError::ResourceLimit {
                            kind: "boundary NNZ upper bound",
                            observed: generated_boundary_terms,
                            budget: max_boundary_nnz,
                        });
                    }
                    let Some(&(target_degree, row)) = indices.get(&target) else {
                        return Err(PrincipalStratumError::MissingBoundaryTarget { cell: target });
                    };
                    if target_degree != degree - 1 {
                        return Err(PrincipalStratumError::UnexpectedBoundaryDegree {
                            source: degree,
                            target: target_degree,
                            cell: cell.parts.clone(),
                        });
                    }
                    boundary.add(row, column, BigInt::from(sign))?;
                }
            }
            differentials.insert(degree, boundary.finish());
        }
        let complex = FiniteChainComplex::new(generator_counts, differentials)
            .map_err(|error| PrincipalStratumError::Modular(error.to_string()))?;
        Ok(Self {
            omega,
            d,
            cells_by_degree,
            indices,
            complex,
        })
    }

    pub fn omega(&self) -> &[u32] {
        &self.omega
    }
    pub fn d(&self) -> u32 {
        self.d
    }
    pub fn cells_by_degree(&self) -> &BTreeMap<i32, Vec<PrincipalStratumCell>> {
        &self.cells_by_degree
    }
    pub fn cell_count(&self) -> usize {
        self.indices.len()
    }
    pub fn complex(&self) -> &FiniteChainComplex {
        &self.complex
    }
    pub fn euler_characteristic(&self) -> i64 {
        self.cells_by_degree
            .iter()
            .map(|(&degree, cells)| {
                if degree.rem_euclid(2) == 0 {
                    cells.len() as i64
                } else {
                    -(cells.len() as i64)
                }
            })
            .sum()
    }
    pub fn initial_nnz(&self) -> usize {
        self.complex
            .degrees()
            .map(|degree| self.complex.differential_or_zero(degree).nnz())
            .sum()
    }
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
    if !is_prime_u64(prime) {
        return Err(PrincipalStratumError::InvalidModulus { modulus: prime });
    }
    let mut ranks = BTreeMap::new();
    let mut solver_stats = BTreeMap::new();
    for degree in model.complex.degrees() {
        let num_vars = model.complex.count(degree - 1);
        let options = SparseModEliminationOptions {
            row_order,
            compute_solution: false,
        };
        // Input order is genuinely boundary-on-demand: the existing Polytool
        // streaming solver consumes one generated boundary column at a time.
        // Reordering modes necessarily collect their rows inside that solver.
        let result = sparse_modular_linear_system_consistency_with_options(
            BoundaryRows {
                model,
                degree,
                column: 0,
                prime,
            },
            num_vars,
            prime,
            options,
        )
        .map_err(|error| PrincipalStratumError::Modular(error.to_string()))?;
        if !result.consistent {
            return Err(PrincipalStratumError::Modular(
                "zero-RHS boundary matrix was unexpectedly inconsistent".into(),
            ));
        }
        ranks.insert(degree, result.rank);
        solver_stats.insert(degree, result.stats);
    }
    let betti_numbers = model
        .complex
        .degrees()
        .map(|degree| {
            let count = model.complex.count(degree);
            let betti = field_betti_number(
                count,
                ranks.get(&degree).copied().unwrap_or(0),
                ranks.get(&(degree + 1)).copied().unwrap_or(0),
            )
            .expect("ranks of a chain complex cannot exceed its count");
            (degree, betti)
        })
        .collect();
    Ok(ModularFieldBetti {
        prime,
        ranks,
        betti_numbers,
        solver_stats,
        evidence_label: "exact F_p dimensions",
    })
}

struct BoundaryRows<'a> {
    model: &'a PrincipalStratumModel,
    degree: i32,
    column: usize,
    prime: u64,
}

impl Iterator for BoundaryRows<'_> {
    type Item = SparseModRow;

    fn next(&mut self) -> Option<Self::Item> {
        let cells = self.model.cells_by_degree.get(&self.degree)?;
        let cell = cells.get(self.column)?;
        self.column += 1;
        let entries = signed_direct_moves(&cell.parts, self.model.d)
            .expect("validated model has bounded direct moves")
            .into_iter()
            .map(|(target, sign)| {
                let &(target_degree, row) = self
                    .model
                    .indices
                    .get(&target)
                    .expect("model boundary target was validated");
                assert_eq!(
                    target_degree,
                    self.degree - 1,
                    "model boundary grading was validated"
                );
                (row, if sign > 0 { 1 } else { self.prime - 1 })
            });
        Some(
            SparseModRow::from_entries(entries, 0, self.prime)
                .expect("prime and model entries were validated"),
        )
    }
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
        let dense = matrix
            .to_dense_with_budget(dense_budget)
            .map_err(PrincipalStratumError::DenseBudget)?;
        let rank = if dense.is_empty() {
            0
        } else {
            sym_poly_core::linear_algebra::rank(
                &dense
                    .into_iter()
                    .map(|row| row.into_iter().map(Ratio::from_integer).collect())
                    .collect::<Vec<Vec<Ratio<BigInt>>>>(),
            )
        };
        ranks.insert(degree, rank);
    }
    Ok(model
        .complex
        .degrees()
        .map(|degree| {
            let betti = field_betti_number(
                model.complex.count(degree),
                ranks.get(&degree).copied().unwrap_or(0),
                ranks.get(&(degree + 1)).copied().unwrap_or(0),
            )
            .expect("chain ranks are bounded");
            (degree, betti)
        })
        .collect())
}

fn validate_input(omega: &[u32], d: u32) -> Result<(), PrincipalStratumError> {
    if omega.is_empty() {
        return Err(PrincipalStratumError::EmptyComposition);
    }
    if let Some((position, _)) = omega.iter().enumerate().find(|(_, part)| **part == 0) {
        return Err(PrincipalStratumError::ZeroPart { position });
    }
    let weight = cell_weight(omega)?;
    if d > i32::MAX as u32 {
        return Err(PrincipalStratumError::GradeOutOfRange { d });
    }
    if d < weight {
        return Err(PrincipalStratumError::WeightTooSmall { d, weight });
    }
    if (d - weight) % 2 != 0 {
        return Err(PrincipalStratumError::ParityMismatch { d, weight });
    }
    Ok(())
}

fn enumerate_cells_bfs(
    omega: &[u32],
    d: u32,
    max_cells: usize,
) -> Result<BTreeSet<Vec<u32>>, PrincipalStratumError> {
    let mut seen = BTreeSet::new();
    let mut queue = VecDeque::new();
    if max_cells == 0 {
        return Err(PrincipalStratumError::ResourceLimit {
            kind: "cell",
            observed: 1,
            budget: 0,
        });
    }
    seen.insert(omega.to_vec());
    queue.push_back(omega.to_vec());
    while let Some(cell) = queue.pop_front() {
        for target in unsigned_direct_moves(&cell, d)? {
            if seen.insert(target.clone()) {
                if seen.len() > max_cells {
                    return Err(PrincipalStratumError::ResourceLimit {
                        kind: "cell",
                        observed: seen.len(),
                        budget: max_cells,
                    });
                }
                queue.push_back(target);
            }
        }
    }
    Ok(seen)
}

/// Deterministic subset automaton for the descendant language.  A state is
/// every prefix length of `omega` that can have been consumed by a particular
/// emitted prefix.  Keeping the full subset is essential: equal sums can give
/// multiple witnesses, and a greedy witness loses valid completions.
fn membership_transition(omega: &[u32], states: &BTreeSet<usize>, part: u32) -> BTreeSet<usize> {
    let mut next = BTreeSet::new();
    for &start in states {
        let mut sum = 0u32;
        for end in start..=omega.len() {
            if end > start {
                sum = match sum.checked_add(omega[end - 1]) {
                    Some(value) => value,
                    None => break,
                };
            }
            // `end == start` means this block consists solely of inserted
            // twos, so it must contain at least one such two.
            if part >= sum && (part - sum).is_multiple_of(2) && (end > start || part >= 2) {
                next.insert(end);
            }
        }
    }
    next
}

/// Decide membership without using closure moves.  This public oracle is also
/// useful to consumers validating externally supplied cells.
pub fn accepted_by_membership_automaton(omega: &[u32], parts: &[u32]) -> bool {
    let mut states = BTreeSet::from([0usize]);
    for &part in parts {
        states = membership_transition(omega, &states, part);
        if states.is_empty() {
            return false;
        }
    }
    states.contains(&omega.len())
}

fn enumerate_cells_automaton(
    omega: &[u32],
    d: u32,
    max_cells: usize,
) -> Result<BTreeSet<Vec<u32>>, PrincipalStratumError> {
    fn visit(
        omega: &[u32],
        d: u32,
        weight: u32,
        states: BTreeSet<usize>,
        parts: &mut Vec<u32>,
        cells: &mut BTreeSet<Vec<u32>>,
        max_cells: usize,
    ) -> Result<(), PrincipalStratumError> {
        if states.contains(&omega.len()) {
            if !cells.insert(parts.clone()) { /* multiple automaton witnesses */ }
            if cells.len() > max_cells {
                return Err(PrincipalStratumError::ResourceLimit {
                    kind: "cell",
                    observed: cells.len(),
                    budget: max_cells,
                });
            }
        }
        let remaining = d - weight;
        for part in 1..=remaining {
            let next = membership_transition(omega, &states, part);
            if next.is_empty() {
                continue;
            }
            parts.push(part);
            visit(omega, d, weight + part, next, parts, cells, max_cells)?;
            parts.pop();
        }
        Ok(())
    }
    if max_cells == 0 {
        return Err(PrincipalStratumError::ResourceLimit {
            kind: "cell",
            observed: 1,
            budget: 0,
        });
    }
    let mut cells = BTreeSet::new();
    visit(
        omega,
        d,
        0,
        BTreeSet::from([0usize]),
        &mut Vec::new(),
        &mut cells,
        max_cells,
    )?;
    Ok(cells)
}

/// Independent count-only DP over automaton states.  It never uses closure
/// moves or boundary assembly, and records the signed Euler grading directly.
pub fn automaton_cell_counts(
    omega: &[u32],
    d: u32,
) -> Result<BTreeMap<i32, usize>, PrincipalStratumError> {
    validate_input(omega, d)?;
    let mut frontier = BTreeMap::<(u32, usize, Vec<usize>), usize>::new();
    frontier.insert((0, 0, vec![0]), 1);
    let mut counts = BTreeMap::<i32, usize>::new();
    while let Some(((weight, length, states), multiplicity)) = frontier.pop_first() {
        if states.binary_search(&omega.len()).is_ok() {
            let degree = cell_degree_from_weight_and_length(weight, length, d)?;
            *counts.entry(degree).or_insert(0) = counts
                .get(&degree)
                .copied()
                .unwrap_or(0)
                .checked_add(multiplicity)
                .ok_or(PrincipalStratumError::ResourceLimit {
                    kind: "count",
                    observed: usize::MAX,
                    budget: usize::MAX,
                })?;
        }
        let state_set = states.into_iter().collect::<BTreeSet<_>>();
        for part in 1..=d - weight {
            let next = membership_transition(omega, &state_set, part);
            if next.is_empty() {
                continue;
            }
            let key = (weight + part, length + 1, next.into_iter().collect());
            let entry = frontier.entry(key).or_insert(0);
            *entry =
                entry
                    .checked_add(multiplicity)
                    .ok_or(PrincipalStratumError::ResourceLimit {
                        kind: "count",
                        observed: usize::MAX,
                        budget: usize::MAX,
                    })?;
        }
    }
    Ok(counts)
}

fn unsigned_direct_moves(cell: &[u32], d: u32) -> Result<Vec<Vec<u32>>, PrincipalStratumError> {
    let mut targets = Vec::with_capacity(cell.len() * 2 + 1);
    for position in 0..cell.len().saturating_sub(1) {
        let mut target = cell.to_vec();
        target[position] = target[position]
            .checked_add(target[position + 1])
            .ok_or(PrincipalStratumError::WeightOverflow)?;
        target.remove(position + 1);
        targets.push(target);
    }
    if cell_weight(cell)?
        .checked_add(2)
        .ok_or(PrincipalStratumError::WeightOverflow)?
        <= d
    {
        for position in 0..=cell.len() {
            let mut target = cell.to_vec();
            target.insert(position, 2);
            targets.push(target);
        }
    }
    Ok(targets)
}

fn signed_direct_moves(cell: &[u32], d: u32) -> Result<Vec<(Vec<u32>, i8)>, PrincipalStratumError> {
    let mut targets = Vec::with_capacity(cell.len() * 2 + 1);
    for position in 0..cell.len().saturating_sub(1) {
        let mut target = cell.to_vec();
        target[position] = target[position]
            .checked_add(target[position + 1])
            .ok_or(PrincipalStratumError::WeightOverflow)?;
        target.remove(position + 1);
        targets.push((target, if position % 2 == 0 { 1 } else { -1 }));
    }
    if cell_weight(cell)?
        .checked_add(2)
        .ok_or(PrincipalStratumError::WeightOverflow)?
        <= d
    {
        for position in 0..=cell.len() {
            let mut target = cell.to_vec();
            target.insert(position, 2);
            targets.push((target, if position % 2 == 0 { 1 } else { -1 }));
        }
    }
    Ok(targets)
}

fn cell_weight(cell: &[u32]) -> Result<u32, PrincipalStratumError> {
    cell.iter().try_fold(0u32, |sum, part| {
        sum.checked_add(*part)
            .ok_or(PrincipalStratumError::WeightOverflow)
    })
}
fn cell_degree(cell: &[u32], d: u32) -> Result<i32, PrincipalStratumError> {
    cell_degree_from_weight_and_length(cell_weight(cell)?, cell.len(), d)
}

fn cell_degree_from_weight_and_length(
    weight: u32,
    length: usize,
    d: u32,
) -> Result<i32, PrincipalStratumError> {
    let d = i32::try_from(d).map_err(|_| PrincipalStratumError::GradeOutOfRange { d })?;
    let weight = i32::try_from(weight).map_err(|_| PrincipalStratumError::WeightOverflow)?;
    let length = i32::try_from(length).map_err(|_| PrincipalStratumError::WeightOverflow)?;
    d.checked_sub(weight)
        .and_then(|value| value.checked_add(length))
        .ok_or(PrincipalStratumError::WeightOverflow)
}

fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n.is_multiple_of(2) {
        return n == 2;
    }
    let mut divisor = 3u64;
    while divisor <= n / divisor {
        if n.is_multiple_of(divisor) {
            return false;
        }
        divisor += 2;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_calibration_has_expected_cells_chain_and_field_groups() {
        let model = PrincipalStratumModel::build(vec![1, 1], 4).unwrap();
        assert_eq!(
            model
                .cells_by_degree
                .iter()
                .map(|(&degree, cells)| (degree, cells.len()))
                .collect::<BTreeMap<_, _>>(),
            [(1, 1), (2, 3), (3, 4), (4, 1)].into_iter().collect()
        );
        assert_eq!(model.cell_count(), 9);
        let modular = modular_field_betti(&model, 101, SparseModRowOrder::Input).unwrap();
        assert_eq!(
            modular.betti_numbers,
            [(1, 0), (2, 0), (3, 1), (4, 0)].into_iter().collect()
        );
        assert_eq!(
            rational_field_betti(&model, 100).unwrap(),
            modular.betti_numbers
        );
    }

    #[test]
    fn membership_subset_automaton_and_count_dp_agree_with_move_oracle() {
        let production = PrincipalStratumModel::build(vec![1, 1, 1], 5).unwrap();
        let closure = PrincipalStratumModel::build_with_bfs_oracle(vec![1, 1, 1], 5).unwrap();
        assert_eq!(production.cells_by_degree, closure.cells_by_degree);
        assert_eq!(
            automaton_cell_counts(&[1, 1, 1], 5).unwrap(),
            production
                .cells_by_degree
                .iter()
                .map(|(&degree, cells)| (degree, cells.len()))
                .collect(),
        );
        assert!(accepted_by_membership_automaton(&[1, 1, 1], &[2, 1]));
        assert!(!accepted_by_membership_automaton(&[1, 2, 1], &[2, 2]));
    }

    #[test]
    fn hostile_limits_grades_and_moduli_are_rejected_before_assembly() {
        assert!(matches!(
            PrincipalStratumModel::build_bounded(vec![1], 1, 0),
            Err(PrincipalStratumError::ResourceLimit { .. })
        ));
        assert!(matches!(
            PrincipalStratumModel::build(vec![1], u32::MAX),
            Err(PrincipalStratumError::GradeOutOfRange { .. })
        ));
        let model = PrincipalStratumModel::build(vec![1, 1], 4).unwrap();
        for modulus in [0, 1, 4] {
            assert!(matches!(
                modular_field_betti(&model, modulus, SparseModRowOrder::Input),
                Err(PrincipalStratumError::InvalidModulus { .. })
            ));
        }
        assert!(matches!(
            PrincipalStratumModel::build_with_limits(
                vec![1, 1],
                4,
                PrincipalStratumBuildLimits {
                    max_cells: 20,
                    max_boundary_nnz: 1
                }
            ),
            Err(PrincipalStratumError::ResourceLimit { .. })
        ));
    }

    #[test]
    fn count_only_dp_confirms_calibration_counts_without_boundaries() {
        for (omega, d, expected_cells, expected_euler) in [
            (vec![3, 1, 1, 3], 18, 8_280usize, -2i64),
            (vec![3, 1, 1, 5], 24, 76_384usize, 0i64),
            (vec![3, 1, 1, 5], 26, 212_900usize, 0i64),
        ] {
            let counts = automaton_cell_counts(&omega, d).unwrap();
            assert_eq!(counts.values().sum::<usize>(), expected_cells);
            let euler = counts.into_iter().map(|(degree, count)| if degree.rem_euclid(2) == 0 { count as i64 } else { -(count as i64) }).sum::<i64>();
            assert_eq!(euler, expected_euler);
        }
    }

    #[test]
    fn invalid_inputs_are_rejected() {
        assert!(matches!(
            PrincipalStratumModel::build(vec![1, 0, 1], 4),
            Err(PrincipalStratumError::ZeroPart { .. })
        ));
        assert!(matches!(
            PrincipalStratumModel::build(vec![1], 2),
            Err(PrincipalStratumError::ParityMismatch { .. })
        ));
    }
}
