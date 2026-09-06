//! Lightweight matroid helpers for exploratory computations.
//!
//! This module intentionally starts with basis-list matroids.  They are not
//! the most efficient representation, but they are flexible enough to mirror
//! the Mathematica `MatroidTools.m` routines and to validate specialized
//! matroid models before promoting them into a larger crate.

use std::collections::{BTreeMap, BTreeSet};

use num_bigint::BigInt;
use num_traits::{One, Zero};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BasisMatroid {
    ground: Vec<usize>,
    rank: usize,
    bases: Vec<Vec<usize>>,
}

impl BasisMatroid {
    pub fn new(ground: Vec<usize>, bases: Vec<Vec<usize>>) -> Result<Self, String> {
        let ground = normalized_set_checked(&ground)?;
        let ground_set: BTreeSet<_> = ground.iter().copied().collect();

        let mut normalized = Vec::with_capacity(bases.len());
        for basis in bases {
            let basis = normalized_set_checked(&basis)?;
            if basis.iter().any(|label| !ground_set.contains(label)) {
                return Err(format!(
                    "basis {basis:?} contains element outside ground set"
                ));
            }
            normalized.push(basis);
        }
        normalized.sort_unstable();
        normalized.dedup();

        if normalized.is_empty() {
            return Err("a matroid must have at least one basis".to_string());
        }

        let rank = normalized.first().map_or(0, Vec::len);
        if normalized.iter().any(|basis| basis.len() != rank) {
            return Err("all bases must have the same cardinality".to_string());
        }

        Ok(Self {
            ground,
            rank,
            bases: normalized,
        })
    }

    pub fn from_bases(bases: Vec<Vec<usize>>) -> Result<Self, String> {
        let mut ground = BTreeSet::new();
        for basis in &bases {
            ground.extend(basis.iter().copied());
        }
        Self::new(ground.into_iter().collect(), bases)
    }

    pub fn uniform(rank: usize, ground: Vec<usize>) -> Result<Self, String> {
        let ground = normalized_set_checked(&ground)?;
        if rank > ground.len() {
            return Err(format!(
                "uniform rank {rank} exceeds ground size {}",
                ground.len()
            ));
        }
        Self::new(ground.clone(), combinations(&ground, rank))
    }

    /// Construct the order-`n` Catalan matroid.
    ///
    /// Its bases are the up-step positions of Dyck paths of semilength `n`.
    /// Equivalently, it has the transversal presentation
    /// `A_i = {i, ..., 2i - 1}`. For positive `n`, element `1` is a coloop and
    /// element `2n` is a loop.
    pub fn catalan(n: usize) -> Result<Self, String> {
        if n == 0 {
            return Self::new(Vec::new(), vec![Vec::new()]);
        }
        let ground: Vec<_> = (1..=2 * n).collect();
        let sets: Vec<Vec<usize>> = (1..=n).map(|i| (i..=2 * i - 1).collect()).collect();
        Self::from_transversal_system(ground, &sets)
    }

    /// Construct the cycle matroid of an undirected multigraph.
    ///
    /// Vertices are labeled `0, ..., vertex_count - 1`; matroid elements are
    /// the one-based positions of the edges in `edges`. Parallel edges are
    /// supported, and a graph loop becomes a matroid loop.
    pub fn from_graph_edges(vertex_count: usize, edges: &[(usize, usize)]) -> Result<Self, String> {
        if edges
            .iter()
            .any(|&(left, right)| left >= vertex_count || right >= vertex_count)
        {
            return Err("graph edge endpoint lies outside the vertex set".to_string());
        }
        if edges.len() >= usize::BITS as usize {
            return Err("too many graph edges for the exact subset backend".to_string());
        }

        let ground: Vec<_> = (1..=edges.len()).collect();
        let rank = graphic_rank(vertex_count, edges);
        let bases = combinations(&ground, rank)
            .into_iter()
            .filter(|basis| graphic_edge_set_is_acyclic(vertex_count, edges, basis))
            .collect();
        Self::new(ground, bases)
    }

    pub fn from_complete_transversal_system(sets: &[Vec<usize>]) -> Result<Self, String> {
        let ground = sets
            .iter()
            .flat_map(|set| set.iter().copied())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        Self::new(ground, complete_transversal_bases(sets))
    }

    /// Construct the transversal matroid presented by an arbitrary set system.
    ///
    /// Unlike [`Self::from_complete_transversal_system`], this constructor does
    /// not require a transversal using every presentation set. Its bases are
    /// the maximum-cardinality partial transversals. Elements of `ground` that
    /// occur in no presentation set become matroid loops.
    pub fn from_transversal_system(
        ground: Vec<usize>,
        sets: &[Vec<usize>],
    ) -> Result<Self, String> {
        let ground = normalized_set_checked(&ground)?;
        let ground_set: BTreeSet<_> = ground.iter().copied().collect();
        let mut normalized_sets = Vec::with_capacity(sets.len());
        for set in sets {
            let set = normalized_set_checked(set)?;
            if set.iter().any(|label| !ground_set.contains(label)) {
                return Err("transversal presentation contains an element outside ground".into());
            }
            normalized_sets.push(set);
        }

        let partial_transversals = partial_transversal_sets(&normalized_sets);
        let rank = partial_transversals.iter().map(Vec::len).max().unwrap_or(0);
        let bases = partial_transversals
            .into_iter()
            .filter(|set| set.len() == rank)
            .collect();
        Self::new(ground, bases)
    }

    pub fn ground(&self) -> &[usize] {
        &self.ground
    }

    pub fn rank(&self) -> usize {
        self.rank
    }

    pub fn bases(&self) -> &[Vec<usize>] {
        &self.bases
    }

    pub fn contains_basis(&self, basis: &[usize]) -> bool {
        normalized_set_checked(basis)
            .ok()
            .is_some_and(|normalized| self.bases.binary_search(&normalized).is_ok())
    }

    pub fn independent_sets(&self) -> Vec<Vec<usize>> {
        let mut independent = BTreeSet::new();
        for basis in &self.bases {
            for subset in all_subsets(basis) {
                independent.insert(subset);
            }
        }
        independent.into_iter().collect()
    }

    /// Exact common-column-sum enumerator for matroid-restricted zero-one rows.
    ///
    /// The coefficient of `q^k` counts ordered `rows`-tuples of independent
    /// sets in which every ground-set element occurs exactly `k` times.
    /// This sparse dynamic program is intended for small research instances.
    pub fn common_column_sum_polynomial(&self, rows: usize) -> Vec<BigInt> {
        if self.ground.is_empty() {
            return vec![BigInt::one()];
        }

        let ground_size = self.ground.len();
        let max_degree = rows.saturating_mul(self.rank) / ground_size;
        let mut coefficients: Vec<_> = (0..=max_degree)
            .map(|column_sum| self.common_column_sum_coefficient(rows, column_sum))
            .collect();

        while coefficients.len() > 1 && coefficients.last().is_some_and(BigInt::is_zero) {
            coefficients.pop();
        }
        coefficients
    }

    /// Exact diagonal-offset enumerator.
    ///
    /// For a nonnegative offset vector `delta` with minimum zero, the
    /// coefficient of `q^k` counts arrays whose column-sum vector is
    /// `k * 1 + delta`. These are the natural states for row-addition
    /// recurrences around the common-column-sum diagonal.
    pub fn column_sum_offset_polynomial(
        &self,
        rows: usize,
        delta: &[usize],
    ) -> Result<Vec<BigInt>, String> {
        if delta.len() != self.ground.len() {
            return Err("offset vector length must equal the ground-set size".to_string());
        }
        if delta.iter().copied().min().unwrap_or(0) != 0 {
            return Err("offset vector must be normalized to have minimum zero".to_string());
        }
        if self.ground.is_empty() {
            return Ok(vec![BigInt::one()]);
        }

        let offset_sum: usize = delta.iter().sum();
        let degree_from_capacity = delta
            .iter()
            .map(|&offset| rows.saturating_sub(offset))
            .min()
            .unwrap_or(0);
        let degree_from_rank =
            rows.saturating_mul(self.rank).saturating_sub(offset_sum) / self.ground.len();
        let max_degree = degree_from_capacity.min(degree_from_rank);
        let mut coefficients: Vec<_> = (0..=max_degree)
            .map(|column_sum| self.column_sum_offset_coefficient(rows, column_sum, delta))
            .collect();

        while coefficients.len() > 1 && coefficients.last().is_some_and(BigInt::is_zero) {
            coefficients.pop();
        }
        Ok(coefficients)
    }

    /// One exact coefficient of the common-column-sum enumerator.
    pub fn common_column_sum_coefficient(&self, rows: usize, column_sum: usize) -> BigInt {
        if self.ground.is_empty() {
            return BigInt::from(column_sum == 0);
        }

        self.column_sum_offset_coefficient(rows, column_sum, &vec![0; self.ground.len()])
    }

    fn column_sum_offset_coefficient(
        &self,
        rows: usize,
        column_sum: usize,
        delta: &[usize],
    ) -> BigInt {
        let ground_size = self.ground.len();
        let target: Vec<_> = delta.iter().map(|&offset| column_sum + offset).collect();
        if target.iter().any(|&sum| sum > rows) {
            return BigInt::zero();
        }
        let independent_positions: Vec<Vec<usize>> = self
            .independent_sets()
            .into_iter()
            .map(|set| {
                set.into_iter()
                    .map(|label| {
                        self.ground
                            .binary_search(&label)
                            .expect("independent-set labels belong to the ground set")
                    })
                    .collect()
            })
            .collect();

        let zero_state = vec![0usize; ground_size];
        let mut counts = BTreeMap::from([(zero_state, BigInt::one())]);

        for _ in 0..rows {
            let mut next = BTreeMap::<Vec<usize>, BigInt>::new();
            for (state, multiplicity) in counts {
                for positions in &independent_positions {
                    if positions
                        .iter()
                        .any(|&position| state[position] == target[position])
                    {
                        continue;
                    }
                    let mut new_state = state.clone();
                    for &position in positions {
                        new_state[position] += 1;
                    }
                    *next.entry(new_state).or_insert_with(BigInt::zero) += &multiplicity;
                }
            }
            counts = next;
        }

        counts.get(&target).cloned().unwrap_or_else(BigInt::zero)
    }

    pub fn rank_of(&self, set: &[usize]) -> Result<usize, String> {
        let set = normalized_set_checked(set)?;
        let set: BTreeSet<_> = set.into_iter().collect();
        if set.iter().any(|label| !self.ground.contains(label)) {
            return Err("set contains element outside ground set".to_string());
        }
        Ok(self
            .bases
            .iter()
            .map(|basis| basis.iter().filter(|label| set.contains(label)).count())
            .max()
            .unwrap_or(0))
    }

    pub fn indicator(&self, set: &[usize]) -> Result<Vec<u8>, String> {
        let set = normalized_set_checked(set)?;
        let set: BTreeSet<_> = set.into_iter().collect();
        if set.iter().any(|label| !self.ground.contains(label)) {
            return Err("set contains element outside ground set".to_string());
        }
        Ok(self
            .ground
            .iter()
            .map(|label| u8::from(set.contains(label)))
            .collect())
    }

    pub fn is_matroid(&self) -> bool {
        for a in &self.bases {
            for b in &self.bases {
                for &x in a.iter().filter(|x| !b.contains(x)) {
                    let mut found_exchange = false;
                    for &y in b.iter().filter(|y| !a.contains(y)) {
                        let Some(exchanged) = sorted_replace(a, x, y) else {
                            continue;
                        };
                        if self.contains_basis(&exchanged) {
                            found_exchange = true;
                            break;
                        }
                    }
                    if !found_exchange {
                        return false;
                    }
                }
            }
        }
        true
    }

    pub fn loops(&self) -> Vec<usize> {
        self.ground
            .iter()
            .copied()
            .filter(|label| self.bases.iter().all(|basis| !basis.contains(label)))
            .collect()
    }

    pub fn coloops(&self) -> Vec<usize> {
        self.ground
            .iter()
            .copied()
            .filter(|label| self.bases.iter().all(|basis| basis.contains(label)))
            .collect()
    }

    pub fn dual(&self) -> Self {
        let bases = self
            .bases
            .iter()
            .map(|basis| {
                self.ground
                    .iter()
                    .copied()
                    .filter(|label| !basis.contains(label))
                    .collect()
            })
            .collect();
        Self::new(self.ground.clone(), bases).expect("dual bases should be valid")
    }

    pub fn delete(&self, label: usize) -> Self {
        if !self.ground.contains(&label) {
            return self.clone();
        }
        if self.coloops().contains(&label) {
            let ground: Vec<_> = self
                .ground
                .iter()
                .copied()
                .filter(|&x| x != label)
                .collect();
            let bases = self
                .bases
                .iter()
                .map(|basis| remove_label(basis, label))
                .collect();
            return Self::new(ground, bases).expect("coloop deletion should be valid");
        }

        let ground: Vec<_> = self
            .ground
            .iter()
            .copied()
            .filter(|&x| x != label)
            .collect();
        let bases = self
            .bases
            .iter()
            .filter(|basis| !basis.contains(&label))
            .cloned()
            .collect();
        Self::new(ground, bases).expect("deletion should be valid")
    }

    pub fn contract(&self, label: usize) -> Self {
        if !self.ground.contains(&label) {
            return self.clone();
        }
        if self.loops().contains(&label) {
            let ground: Vec<_> = self
                .ground
                .iter()
                .copied()
                .filter(|&x| x != label)
                .collect();
            return Self::new(ground, self.bases.clone())
                .expect("loop contraction should be valid");
        }

        let ground: Vec<_> = self
            .ground
            .iter()
            .copied()
            .filter(|&x| x != label)
            .collect();
        let bases = self
            .bases
            .iter()
            .filter(|basis| basis.contains(&label))
            .map(|basis| remove_label(basis, label))
            .collect();
        Self::new(ground, bases).expect("contraction should be valid")
    }

    pub fn delete_many(&self, labels: &[usize]) -> Self {
        labels
            .iter()
            .fold(self.clone(), |matroid, &label| matroid.delete(label))
    }

    pub fn contract_many(&self, labels: &[usize]) -> Self {
        labels
            .iter()
            .fold(self.clone(), |matroid, &label| matroid.contract(label))
    }

    pub fn exchange_insertions(&self, basis: &[usize], insert: usize) -> Vec<(usize, Vec<usize>)> {
        if basis.contains(&insert) {
            return vec![(insert, normalized_set(basis))];
        }
        let mut out = Vec::new();
        for &remove in basis {
            let Some(candidate) = sorted_replace(basis, remove, insert) else {
                continue;
            };
            if self.contains_basis(&candidate) {
                out.push((remove, candidate));
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    pub fn restricted_exchange_graph(&self, insert: usize) -> RestrictedExchangeGraph {
        let mut sources = Vec::new();
        let mut targets = Vec::new();
        let mut target_index = BTreeMap::new();

        for basis in &self.bases {
            if basis.contains(&insert) {
                target_index.insert(basis.clone(), targets.len());
                targets.push(basis.clone());
            } else {
                sources.push(basis.clone());
            }
        }

        let edges = sources
            .iter()
            .map(|source| {
                let mut row = Vec::new();
                for (_, candidate) in self.exchange_insertions(source, insert) {
                    if let Some(&target) = target_index.get(&candidate) {
                        row.push(target);
                    }
                }
                row.sort_unstable();
                row.dedup();
                row
            })
            .collect();

        RestrictedExchangeGraph {
            sources,
            targets,
            edges,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RestrictedExchangeGraph {
    pub sources: Vec<Vec<usize>>,
    pub targets: Vec<Vec<usize>>,
    pub edges: Vec<Vec<usize>>,
}

impl RestrictedExchangeGraph {
    pub fn has_source_covering_matching(&self) -> bool {
        has_complete_matching(&self.edges, self.targets.len())
    }

    pub fn hall_defect(&self) -> Option<(Vec<usize>, Vec<usize>)> {
        if self.sources.len() >= usize::BITS as usize {
            return None;
        }
        let source_count = self.sources.len();
        for mask in 1usize..(1usize << source_count) {
            let mut sources = Vec::new();
            let mut neighbors = BTreeSet::new();
            for source in 0..source_count {
                if (mask >> source) & 1 == 1 {
                    sources.push(source);
                    neighbors.extend(self.edges[source].iter().copied());
                }
            }
            if neighbors.len() < sources.len() {
                return Some((sources, neighbors.into_iter().collect()));
            }
        }
        None
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PrefixLpm {
    ground: Vec<usize>,
    rank: usize,
    lower_prefix: Vec<usize>,
    upper_prefix: Vec<usize>,
}

impl PrefixLpm {
    pub fn new(
        ground: Vec<usize>,
        rank: usize,
        lower_prefix: Vec<usize>,
        upper_prefix: Vec<usize>,
    ) -> Result<Self, String> {
        let ground = normalized_set_checked(&ground)?;
        let expected_len = ground.len() + 1;
        if lower_prefix.len() != expected_len || upper_prefix.len() != expected_len {
            return Err(format!(
                "prefix bounds must have length {}, got lower={} upper={}",
                expected_len,
                lower_prefix.len(),
                upper_prefix.len()
            ));
        }
        if lower_prefix[0] != 0 || upper_prefix[0] != 0 {
            return Err("prefix bounds must start at 0".to_string());
        }
        if lower_prefix[ground.len()] > rank || rank > upper_prefix[ground.len()] {
            return Err("final prefix bounds must contain the rank".to_string());
        }
        for index in 0..=ground.len() {
            let lower = lower_prefix[index];
            let upper = upper_prefix[index];
            if lower > upper {
                return Err(format!("lower prefix bound exceeds upper at index {index}"));
            }
            if lower > rank || upper > rank {
                return Err(format!("prefix bound exceeds rank at index {index}"));
            }
            if lower > index {
                return Err(format!(
                    "lower prefix bound {lower} exceeds prefix size {index}"
                ));
            }
            if index > 0 && (lower_prefix[index - 1] > lower || upper_prefix[index - 1] > upper) {
                return Err("prefix bounds must be weakly increasing".to_string());
            }
        }

        Ok(Self {
            ground,
            rank,
            lower_prefix,
            upper_prefix,
        })
    }

    pub fn from_intervals(intervals: &[(usize, usize)]) -> Result<Self, String> {
        if intervals.is_empty() {
            return Self::new(Vec::new(), 0, vec![0], vec![0]);
        }
        let min_label = intervals.iter().map(|&(lower, _)| lower).min().unwrap();
        let max_label = intervals.iter().map(|&(_, upper)| upper).max().unwrap();
        Self::from_intervals_on_ground((min_label..=max_label).collect(), intervals)
    }

    pub fn from_intervals_on_ground(
        ground: Vec<usize>,
        intervals: &[(usize, usize)],
    ) -> Result<Self, String> {
        let ground = normalized_set_checked(&ground)?;
        for (index, &(lower, upper)) in intervals.iter().enumerate() {
            if lower > upper {
                return Err(format!(
                    "interval {index} has lower endpoint above upper endpoint"
                ));
            }
            if !ground.iter().any(|&label| lower <= label && label <= upper) {
                return Err(format!("interval {index} contains no ground label"));
            }
            if index > 0 {
                let (prev_lower, prev_upper) = intervals[index - 1];
                if prev_lower > lower || prev_upper > upper {
                    return Err("interval endpoints must be weakly increasing".to_string());
                }
            }
        }

        let rank = intervals.len();
        let mut lower_prefix = Vec::with_capacity(ground.len() + 1);
        let mut upper_prefix = Vec::with_capacity(ground.len() + 1);
        lower_prefix.push(0);
        upper_prefix.push(0);
        for &label in &ground {
            lower_prefix.push(
                intervals
                    .iter()
                    .filter(|&&(_, upper)| upper <= label)
                    .count(),
            );
            upper_prefix.push(
                intervals
                    .iter()
                    .filter(|&&(lower, _)| lower <= label)
                    .count(),
            );
        }

        let lpm = Self::new(ground, rank, lower_prefix, upper_prefix)?;
        if lpm.enumerate_bases().is_empty() {
            return Err("prefix interval data has no bases".to_string());
        }
        Ok(lpm)
    }

    pub fn from_interval_sets(sets: &[Vec<usize>]) -> Result<Self, String> {
        let mut ground = BTreeSet::new();
        let mut intervals = Vec::with_capacity(sets.len());
        for (index, set) in sets.iter().enumerate() {
            let set = normalized_set_checked(set)?;
            if set.is_empty() {
                return Err(format!("interval set {index} is empty"));
            }
            for window in set.windows(2) {
                if window[1] != window[0] + 1 {
                    return Err(format!("interval set {index} is not contiguous: {set:?}"));
                }
            }
            ground.extend(set.iter().copied());
            intervals.push((*set.first().unwrap(), *set.last().unwrap()));
        }
        Self::from_intervals_on_ground(ground.into_iter().collect(), &intervals)
    }

    pub fn from_basis_matroid(matroid: &BasisMatroid) -> Result<Self, String> {
        let ground = matroid.ground().to_vec();
        let mut lower_prefix = vec![usize::MAX; ground.len() + 1];
        let mut upper_prefix = vec![0; ground.len() + 1];
        lower_prefix[0] = 0;

        for basis in matroid.bases() {
            let basis_set: BTreeSet<_> = basis.iter().copied().collect();
            let mut count = 0;
            for (index, label) in ground.iter().enumerate() {
                if basis_set.contains(label) {
                    count += 1;
                }
                lower_prefix[index + 1] = lower_prefix[index + 1].min(count);
                upper_prefix[index + 1] = upper_prefix[index + 1].max(count);
            }
        }

        let lpm = Self::new(ground, matroid.rank(), lower_prefix, upper_prefix)?;
        if lpm.enumerate_bases() != matroid.bases() {
            return Err("basis family is not described exactly by prefix bounds".to_string());
        }
        Ok(lpm)
    }

    pub fn ground(&self) -> &[usize] {
        &self.ground
    }

    pub fn rank(&self) -> usize {
        self.rank
    }

    pub fn lower_prefix(&self) -> &[usize] {
        &self.lower_prefix
    }

    pub fn upper_prefix(&self) -> &[usize] {
        &self.upper_prefix
    }

    pub fn contains_base(&self, basis: &[usize]) -> bool {
        let Ok(basis) = normalized_set_checked(basis) else {
            return false;
        };
        if basis.len() != self.rank {
            return false;
        }
        let basis_set: BTreeSet<_> = basis.into_iter().collect();
        if basis_set.iter().any(|label| !self.ground.contains(label)) {
            return false;
        }

        let mut count = 0;
        for (index, label) in self.ground.iter().enumerate() {
            if basis_set.contains(label) {
                count += 1;
            }
            if count < self.lower_prefix[index + 1] || self.upper_prefix[index + 1] < count {
                return false;
            }
        }
        count == self.rank
    }

    pub fn enumerate_bases(&self) -> Vec<Vec<usize>> {
        fn rec(
            lpm: &PrefixLpm,
            index: usize,
            count: usize,
            current: &mut Vec<usize>,
            out: &mut Vec<Vec<usize>>,
        ) {
            if index == lpm.ground.len() {
                if count == lpm.rank {
                    out.push(current.clone());
                }
                return;
            }
            let remaining_after_this = lpm.ground.len() - index - 1;

            let excluded_count = count;
            if lpm.lower_prefix[index + 1] <= excluded_count
                && excluded_count <= lpm.upper_prefix[index + 1]
                && excluded_count <= lpm.rank
                && excluded_count + remaining_after_this >= lpm.rank
            {
                rec(lpm, index + 1, excluded_count, current, out);
            }

            let included_count = count + 1;
            if included_count <= lpm.rank
                && lpm.lower_prefix[index + 1] <= included_count
                && included_count <= lpm.upper_prefix[index + 1]
                && included_count + remaining_after_this >= lpm.rank
            {
                current.push(lpm.ground[index]);
                rec(lpm, index + 1, included_count, current, out);
                current.pop();
            }
        }

        let mut out = Vec::new();
        rec(self, 0, 0, &mut Vec::new(), &mut out);
        out.sort_unstable();
        out
    }

    pub fn to_basis_matroid(&self) -> BasisMatroid {
        BasisMatroid::new(self.ground.clone(), self.enumerate_bases())
            .expect("prefix LPM bases should form a basis matroid")
    }

    pub fn delete(&self, label: usize) -> Result<Self, String> {
        Self::from_basis_matroid(&self.to_basis_matroid().delete(label))
    }

    pub fn contract(&self, label: usize) -> Result<Self, String> {
        Self::from_basis_matroid(&self.to_basis_matroid().contract(label))
    }

    pub fn dual_by_complement(&self) -> Result<Self, String> {
        Self::from_basis_matroid(&self.to_basis_matroid().dual())
    }

    pub fn restricted_exchange_graph(&self, insert: usize) -> RestrictedExchangeGraph {
        self.to_basis_matroid().restricted_exchange_graph(insert)
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SignatureFiber {
    common: Vec<usize>,
    singleton: Vec<usize>,
    k_sigma: PrefixLpm,
}

impl SignatureFiber {
    fn split_signature(
        lpm: &PrefixLpm,
        signature: &[usize],
    ) -> Result<(Vec<usize>, Vec<usize>, usize), String> {
        if signature.len() != 2 * lpm.rank() {
            return Err(format!(
                "signature length {} is not twice the LPM rank {}",
                signature.len(),
                lpm.rank()
            ));
        }

        let mut multiplicities = BTreeMap::new();
        for &label in signature {
            if !lpm.ground().contains(&label) {
                return Err(format!(
                    "signature label {label} is outside the LPM ground set"
                ));
            }
            *multiplicities.entry(label).or_insert(0usize) += 1;
        }

        let mut common = Vec::new();
        let mut singleton = Vec::new();
        for (label, multiplicity) in multiplicities {
            match multiplicity {
                1 => singleton.push(label),
                2 => common.push(label),
                _ => {
                    return Err(format!(
                        "signature label {label} has multiplicity {multiplicity}, expected 1 or 2"
                    ));
                }
            }
        }

        if common.len() > lpm.rank() {
            return Err("too many common labels for the LPM rank".to_string());
        }
        let fiber_rank = lpm.rank() - common.len();
        if singleton.len() != 2 * fiber_rank {
            return Err(format!(
                "singleton count {} is incompatible with fiber rank {fiber_rank}",
                singleton.len()
            ));
        }

        Ok((common, singleton, fiber_rank))
    }

    pub fn k_sigma_from_prefix_bounds(
        lpm: &PrefixLpm,
        signature: &[usize],
    ) -> Result<PrefixLpm, String> {
        let (common, singleton, fiber_rank) = Self::split_signature(lpm, signature)?;
        let common_set: BTreeSet<_> = common.iter().copied().collect();
        let singleton_set: BTreeSet<_> = singleton.iter().copied().collect();

        let mut lower_prefix = vec![0usize; singleton.len() + 1];
        let mut upper_prefix: Vec<_> = (0..=singleton.len())
            .map(|prefix_size| prefix_size.min(fiber_rank))
            .collect();

        let mut common_seen = 0isize;
        let mut singleton_seen = 0isize;
        for ground_prefix in 0..=lpm.ground().len() {
            let lower = lpm.lower_prefix()[ground_prefix] as isize;
            let upper = lpm.upper_prefix()[ground_prefix] as isize;
            let d_prefix = singleton_seen as usize;

            let forced_lower = [0, lower - common_seen, singleton_seen + common_seen - upper]
                .into_iter()
                .max()
                .unwrap();
            let forced_upper = [
                fiber_rank as isize,
                singleton_seen,
                upper - common_seen,
                singleton_seen + common_seen - lower,
            ]
            .into_iter()
            .min()
            .unwrap();

            if forced_upper < 0 {
                return Err("direct signature prefix bounds are infeasible".to_string());
            }
            lower_prefix[d_prefix] = lower_prefix[d_prefix].max(forced_lower as usize);
            upper_prefix[d_prefix] = upper_prefix[d_prefix].min(forced_upper as usize);

            if let Some(&label) = lpm.ground().get(ground_prefix) {
                if common_set.contains(&label) {
                    common_seen += 1;
                }
                if singleton_set.contains(&label) {
                    singleton_seen += 1;
                }
            }
        }

        for index in 1..lower_prefix.len() {
            lower_prefix[index] = lower_prefix[index].max(lower_prefix[index - 1]);
        }
        for index in (0..upper_prefix.len() - 1).rev() {
            upper_prefix[index] = upper_prefix[index].min(upper_prefix[index + 1]);
        }

        let k_sigma = PrefixLpm::new(singleton, fiber_rank, lower_prefix, upper_prefix)?;
        if k_sigma.enumerate_bases().is_empty() {
            return Err("direct signature prefix bounds have no bases".to_string());
        }
        Ok(k_sigma)
    }

    pub fn from_lpm_and_signature(lpm: &PrefixLpm, signature: &[usize]) -> Result<Self, String> {
        let (common, singleton, fiber_rank) = Self::split_signature(lpm, signature)?;
        let mut bases = Vec::new();
        for a in combinations(&singleton, fiber_rank) {
            let complement = sorted_complement(&singleton, &a);
            let q = sorted_union(&common, &a);
            let b = sorted_union(&common, &complement);
            if lpm.contains_base(&q) && lpm.contains_base(&b) {
                bases.push(a);
            }
        }
        if bases.is_empty() {
            return Err("signature fiber has no feasible bases".to_string());
        }

        let basis_matroid = BasisMatroid::new(singleton.clone(), bases)?;
        let k_sigma = Self::k_sigma_from_prefix_bounds(lpm, signature)?;
        if k_sigma.enumerate_bases() != basis_matroid.bases() {
            return Err(
                "direct signature prefix bounds do not match enumerated signature bases"
                    .to_string(),
            );
        }
        Ok(Self {
            common,
            singleton,
            k_sigma,
        })
    }

    pub fn common(&self) -> &[usize] {
        &self.common
    }

    pub fn singleton(&self) -> &[usize] {
        &self.singleton
    }

    pub fn k_sigma(&self) -> &PrefixLpm {
        &self.k_sigma
    }

    pub fn source_sets(&self, x: usize) -> Vec<Vec<usize>> {
        self.k_sigma
            .enumerate_bases()
            .into_iter()
            .filter(|basis| !basis.contains(&x))
            .collect()
    }

    pub fn target_sets(&self, x: usize) -> Vec<Vec<usize>> {
        self.k_sigma
            .enumerate_bases()
            .into_iter()
            .filter(|basis| basis.contains(&x))
            .collect()
    }

    pub fn pair_from_choice(&self, a: &[usize]) -> Result<(Vec<usize>, Vec<usize>), String> {
        if !self.k_sigma.contains_base(a) {
            return Err("choice is not a base of K_Sigma".to_string());
        }
        let a = normalized_set_checked(a)?;
        let complement = sorted_complement(&self.singleton, &a);
        Ok((
            sorted_union(&self.common, &a),
            sorted_union(&self.common, &complement),
        ))
    }

    pub fn is_self_dual_by_complement(&self) -> bool {
        let bases: BTreeSet<_> = self.k_sigma.enumerate_bases().into_iter().collect();
        bases.iter().all(|basis| {
            let complement = sorted_complement(&self.singleton, basis);
            bases.contains(&complement)
        })
    }

    pub fn endpoint_restricted_lpm(&self, dummy: usize) -> Result<PrefixLpm, String> {
        if self.common.contains(&dummy) {
            Ok(self.k_sigma.clone())
        } else if self.singleton.contains(&dummy) {
            self.k_sigma.delete(dummy)
        } else {
            Err(format!(
                "dummy label {dummy} is not present in the signature"
            ))
        }
    }

    pub fn endpoint_restricted_hall(&self, dummy: usize, x: usize) -> Result<bool, String> {
        let h = self.endpoint_restricted_lpm(dummy)?;
        if !h.ground().contains(&x) {
            return Ok(false);
        }
        Ok(h.restricted_exchange_graph(x)
            .has_source_covering_matching())
    }
}

pub fn complete_transversal_bases(sets: &[Vec<usize>]) -> Vec<Vec<usize>> {
    fn rec(
        sets: &[Vec<usize>],
        index: usize,
        used: &mut BTreeSet<usize>,
        current: &mut Vec<usize>,
        out: &mut Vec<Vec<usize>>,
    ) {
        if index == sets.len() {
            out.push(normalized_set(current));
            return;
        }
        for &label in &sets[index] {
            if used.insert(label) {
                current.push(label);
                rec(sets, index + 1, used, current, out);
                current.pop();
                used.remove(&label);
            }
        }
    }

    let normalized_sets: Vec<_> = sets.iter().map(|set| normalized_set(set)).collect();
    let mut out = Vec::new();
    rec(
        &normalized_sets,
        0,
        &mut BTreeSet::new(),
        &mut Vec::new(),
        &mut out,
    );
    out.sort_unstable();
    out.dedup();
    out
}

fn partial_transversal_sets(sets: &[Vec<usize>]) -> Vec<Vec<usize>> {
    fn rec(
        sets: &[Vec<usize>],
        index: usize,
        used: &mut BTreeSet<usize>,
        out: &mut BTreeSet<Vec<usize>>,
    ) {
        if index == sets.len() {
            out.insert(used.iter().copied().collect());
            return;
        }

        rec(sets, index + 1, used, out);
        for &label in &sets[index] {
            if used.insert(label) {
                rec(sets, index + 1, used, out);
                used.remove(&label);
            }
        }
    }

    let mut out = BTreeSet::new();
    rec(sets, 0, &mut BTreeSet::new(), &mut out);
    out.into_iter().collect()
}

fn graphic_rank(vertex_count: usize, edges: &[(usize, usize)]) -> usize {
    let mut parents: Vec<_> = (0..vertex_count).collect();
    edges
        .iter()
        .filter(|&&(left, right)| left != right && union_vertices(&mut parents, left, right))
        .count()
}

fn graphic_edge_set_is_acyclic(
    vertex_count: usize,
    edges: &[(usize, usize)],
    edge_labels: &[usize],
) -> bool {
    let mut parents: Vec<_> = (0..vertex_count).collect();
    edge_labels.iter().all(|&label| {
        let (left, right) = edges[label - 1];
        left != right && union_vertices(&mut parents, left, right)
    })
}

fn union_vertices(parents: &mut [usize], left: usize, right: usize) -> bool {
    let left_root = find_vertex_root(parents, left);
    let right_root = find_vertex_root(parents, right);
    if left_root == right_root {
        return false;
    }
    parents[right_root] = left_root;
    true
}

fn find_vertex_root(parents: &mut [usize], vertex: usize) -> usize {
    let mut root = vertex;
    while parents[root] != root {
        root = parents[root];
    }
    let mut current = vertex;
    while parents[current] != current {
        let next = parents[current];
        parents[current] = root;
        current = next;
    }
    root
}

fn normalized_set(set: &[usize]) -> Vec<usize> {
    let mut out = set.to_vec();
    out.sort_unstable();
    out.dedup();
    out
}

fn normalized_set_checked(set: &[usize]) -> Result<Vec<usize>, String> {
    let mut out = set.to_vec();
    out.sort_unstable();
    for pair in out.windows(2) {
        if pair[0] == pair[1] {
            return Err(format!("set {set:?} contains duplicate label {}", pair[0]));
        }
    }
    Ok(out)
}

fn all_subsets(set: &[usize]) -> Vec<Vec<usize>> {
    let mut out = vec![Vec::new()];
    for &label in set {
        let mut with_label = out.clone();
        for subset in &mut with_label {
            subset.push(label);
        }
        out.extend(with_label);
    }
    out.sort_unstable();
    out.dedup();
    out
}

fn combinations(set: &[usize], k: usize) -> Vec<Vec<usize>> {
    fn rec(
        set: &[usize],
        k: usize,
        start: usize,
        current: &mut Vec<usize>,
        out: &mut Vec<Vec<usize>>,
    ) {
        if current.len() == k {
            out.push(current.clone());
            return;
        }
        let needed = k - current.len();
        if set.len().saturating_sub(start) < needed {
            return;
        }
        for index in start..=set.len() - needed {
            current.push(set[index]);
            rec(set, k, index + 1, current, out);
            current.pop();
        }
    }

    let mut out = Vec::new();
    rec(set, k, 0, &mut Vec::new(), &mut out);
    out
}

fn remove_label(set: &[usize], label: usize) -> Vec<usize> {
    set.iter().copied().filter(|&x| x != label).collect()
}

fn sorted_replace(set: &[usize], remove: usize, insert: usize) -> Option<Vec<usize>> {
    if remove != insert && set.contains(&insert) {
        return None;
    }
    let mut found = false;
    let mut out = Vec::with_capacity(set.len());
    for &label in set {
        if label == remove && !found {
            out.push(insert);
            found = true;
        } else {
            out.push(label);
        }
    }
    found.then(|| normalized_set(&out))
}

fn sorted_union(a: &[usize], b: &[usize]) -> Vec<usize> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    out.extend_from_slice(a);
    out.extend_from_slice(b);
    out.sort_unstable();
    out
}

fn sorted_complement(ground: &[usize], subset: &[usize]) -> Vec<usize> {
    ground
        .iter()
        .copied()
        .filter(|label| !subset.contains(label))
        .collect()
}

fn has_complete_matching(edges: &[Vec<usize>], target_count: usize) -> bool {
    fn dfs(
        source: usize,
        edges: &[Vec<usize>],
        seen: &mut [bool],
        matched_by: &mut [Option<usize>],
    ) -> bool {
        for &target in &edges[source] {
            if seen[target] {
                continue;
            }
            seen[target] = true;
            if matched_by[target].is_none()
                || dfs(matched_by[target].unwrap(), edges, seen, matched_by)
            {
                matched_by[target] = Some(source);
                return true;
            }
        }
        false
    }

    let mut matched_by = vec![None; target_count];
    for source in 0..edges.len() {
        let mut seen = vec![false; target_count];
        if !dfs(source, edges, &mut seen, &mut matched_by) {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_matroid_basic_operations() {
        let matroid = BasisMatroid::uniform(2, vec![1, 2, 3, 4]).unwrap();
        assert!(matroid.is_matroid());
        assert_eq!(matroid.rank(), 2);
        assert!(matroid.loops().is_empty());
        assert!(matroid.coloops().is_empty());
        assert_eq!(matroid.dual().bases(), matroid.bases());
        assert_eq!(matroid.rank_of(&[1, 2, 3]).unwrap(), 2);
        assert_eq!(matroid.indicator(&[2, 4]).unwrap(), vec![0, 1, 0, 1]);

        let deleted = matroid.delete(4);
        assert_eq!(deleted.bases(), &[vec![1, 2], vec![1, 3], vec![2, 3]]);

        let contracted = matroid.contract(4);
        assert_eq!(contracted.bases(), &[vec![1], vec![2], vec![3]]);
    }

    #[test]
    fn graphic_matroid_supports_cycles_parallel_edges_and_loops() {
        let triangle = BasisMatroid::from_graph_edges(3, &[(0, 1), (1, 2), (0, 2)]).unwrap();
        assert_eq!(triangle, BasisMatroid::uniform(2, vec![1, 2, 3]).unwrap());

        let multigraph = BasisMatroid::from_graph_edges(2, &[(0, 1), (0, 1), (0, 0)]).unwrap();
        assert_eq!(multigraph.bases(), &[vec![1], vec![2]]);
        assert_eq!(multigraph.loops(), vec![3]);
    }

    #[test]
    fn catalan_matroid_has_dyck_bases_and_forced_elements() {
        let matroid = BasisMatroid::catalan(3).unwrap();
        assert_eq!(matroid.rank(), 3);
        assert_eq!(matroid.bases().len(), 5);
        assert_eq!(matroid.coloops(), vec![1]);
        assert_eq!(matroid.loops(), vec![6]);

        let catalan_counts = [1usize, 2, 5, 14, 42, 132];
        let core_independent_counts = [1usize, 3, 10, 35, 126, 462];
        for n in 1..=6 {
            let catalan = BasisMatroid::catalan(n).unwrap();
            assert_eq!(catalan.bases().len(), catalan_counts[n - 1]);
            let core = catalan.delete(2 * n).contract(1);
            assert_eq!(
                core.independent_sets().len(),
                core_independent_counts[n - 1]
            );
        }
    }

    #[test]
    fn common_column_sum_polynomial_matches_u24_examples() {
        let matroid = BasisMatroid::uniform(2, vec![1, 2, 3, 4]).unwrap();
        assert_eq!(
            matroid.common_column_sum_polynomial(4),
            vec![BigInt::from(1), BigInt::from(204), BigInt::from(90)]
        );
        assert_eq!(
            matroid.common_column_sum_polynomial(6),
            vec![
                BigInt::from(1),
                BigInt::from(1170),
                BigInt::from(20610),
                BigInt::from(1860),
            ]
        );
    }

    #[test]
    fn column_sum_offset_polynomial_matches_catalan_core_example() {
        let matroid = BasisMatroid::new(
            vec![1, 2, 3, 4],
            vec![vec![1, 2], vec![1, 3], vec![1, 4], vec![2, 3], vec![2, 4]],
        )
        .unwrap();
        assert_eq!(
            matroid.column_sum_offset_polynomial(4, &[0, 1, 1, 1]),
            Ok(vec![BigInt::from(48), BigInt::from(72)])
        );
        assert_eq!(
            matroid.column_sum_offset_polynomial(4, &[1, 1, 1, 1]),
            Err("offset vector must be normalized to have minimum zero".to_string())
        );
    }

    #[test]
    fn coloop_deletion_and_loop_contraction_match_mathematica_convention() {
        let matroid = BasisMatroid::new(vec![1, 2, 3], vec![vec![1, 2], vec![1, 3]]).unwrap();
        assert_eq!(matroid.coloops(), vec![1]);
        assert_eq!(matroid.delete(1).bases(), &[vec![2], vec![3]]);

        let with_loop = BasisMatroid::new(vec![1, 2, 3], vec![vec![1], vec![2]]).unwrap();
        assert_eq!(with_loop.loops(), vec![3]);
        assert_eq!(with_loop.contract(3).bases(), &[vec![1], vec![2]]);
    }

    #[test]
    fn restricted_exchange_graph_detects_matching() {
        let matroid = BasisMatroid::new(vec![1, 2, 3], vec![vec![1], vec![2], vec![3]]).unwrap();
        let graph = matroid.restricted_exchange_graph(3);
        assert_eq!(graph.sources, vec![vec![1], vec![2]]);
        assert_eq!(graph.targets, vec![vec![3]]);
        assert!(!graph.has_source_covering_matching());
        assert!(graph.hall_defect().is_some());
    }

    #[test]
    fn complete_transversal_system_matches_mathematica_example() {
        let bases = complete_transversal_bases(&[vec![1, 2, 3, 4], vec![3, 4, 5]]);
        assert_eq!(
            bases,
            vec![
                vec![1, 3],
                vec![1, 4],
                vec![1, 5],
                vec![2, 3],
                vec![2, 4],
                vec![2, 5],
                vec![3, 4],
                vec![3, 5],
                vec![4, 5],
            ]
        );

        let matroid =
            BasisMatroid::from_complete_transversal_system(&[vec![1, 2, 3, 4], vec![3, 4, 5]])
                .unwrap();
        assert!(matroid.is_matroid());
    }

    #[test]
    fn arbitrary_transversal_system_keeps_maximum_partial_transversals() {
        let matroid = BasisMatroid::from_transversal_system(
            vec![1, 2, 3],
            &[vec![1], vec![1, 2], vec![1, 2]],
        )
        .unwrap();
        assert_eq!(matroid.rank(), 2);
        assert_eq!(matroid.bases(), &[vec![1, 2]]);
        assert_eq!(matroid.loops(), vec![3]);
        assert!(matroid.is_matroid());
    }

    #[test]
    fn duplicate_labels_are_rejected() {
        assert!(BasisMatroid::new(vec![1, 1, 2], vec![vec![1]]).is_err());
        assert!(BasisMatroid::new(vec![1, 2], vec![vec![1, 1]]).is_err());
        assert!(BasisMatroid::uniform(1, vec![1, 1, 2]).is_err());
        assert!(!BasisMatroid::uniform(1, vec![1, 2])
            .unwrap()
            .contains_basis(&[1, 1]));
    }

    #[test]
    fn absent_minor_labels_are_noops() {
        let matroid = BasisMatroid::uniform(1, vec![1, 2]).unwrap();
        assert_eq!(matroid.delete(9), matroid);
        assert_eq!(matroid.contract(9), matroid);
    }

    #[test]
    fn prefix_lpm_from_interval_sets_matches_transversal_example() {
        let lpm = PrefixLpm::from_interval_sets(&[vec![1, 2, 3, 4], vec![3, 4, 5]]).unwrap();
        let expected = vec![
            vec![1, 3],
            vec![1, 4],
            vec![1, 5],
            vec![2, 3],
            vec![2, 4],
            vec![2, 5],
            vec![3, 4],
            vec![3, 5],
            vec![4, 5],
        ];

        assert_eq!(lpm.ground(), &[1, 2, 3, 4, 5]);
        assert_eq!(lpm.rank(), 2);
        assert_eq!(lpm.enumerate_bases(), expected);
        assert_eq!(lpm.to_basis_matroid().bases(), expected.as_slice());
        assert!(lpm.contains_base(&[2, 5]));
        assert!(!lpm.contains_base(&[1, 2]));
    }

    #[test]
    fn prefix_lpm_minors_match_basis_matroid_minors() {
        let lpm = PrefixLpm::from_intervals(&[(1, 4), (3, 5)]).unwrap();
        let matroid = lpm.to_basis_matroid();

        assert_eq!(lpm.delete(1).unwrap().to_basis_matroid(), matroid.delete(1));
        assert_eq!(
            lpm.contract(4).unwrap().to_basis_matroid(),
            matroid.contract(4)
        );
        assert_eq!(
            lpm.dual_by_complement().unwrap().to_basis_matroid(),
            matroid.dual()
        );
    }

    #[test]
    fn prefix_lpm_label_vector_model_example() {
        let lpm = PrefixLpm::from_intervals(&[(1, 3), (2, 4), (4, 5)]).unwrap();
        assert_eq!(
            lpm.enumerate_bases(),
            vec![
                vec![1, 2, 4],
                vec![1, 2, 5],
                vec![1, 3, 4],
                vec![1, 3, 5],
                vec![1, 4, 5],
                vec![2, 3, 4],
                vec![2, 3, 5],
                vec![2, 4, 5],
                vec![3, 4, 5],
            ]
        );
    }

    #[test]
    fn prefix_lpm_rejects_bad_interval_data() {
        assert!(PrefixLpm::from_intervals_on_ground(vec![2], &[(1, 1)]).is_err());
        assert!(PrefixLpm::from_interval_sets(&[vec![1, 3]]).is_err());
        assert!(PrefixLpm::from_intervals(&[(2, 3), (1, 4)]).is_err());
    }

    #[test]
    fn signature_fiber_recovers_uniform_k_sigma() {
        let lpm = PrefixLpm::from_interval_sets(&[vec![1, 2, 3, 4], vec![3, 4, 5]]).unwrap();
        let fiber = SignatureFiber::from_lpm_and_signature(&lpm, &[1, 3, 4, 5]).unwrap();

        assert!(fiber.common().is_empty());
        assert_eq!(fiber.singleton(), &[1, 3, 4, 5]);
        assert_eq!(
            fiber.k_sigma().enumerate_bases(),
            vec![
                vec![1, 3],
                vec![1, 4],
                vec![1, 5],
                vec![3, 4],
                vec![3, 5],
                vec![4, 5],
            ]
        );
        assert!(fiber.is_self_dual_by_complement());
        assert_eq!(
            fiber.source_sets(1),
            vec![vec![3, 4], vec![3, 5], vec![4, 5]]
        );
        assert_eq!(
            fiber.target_sets(1),
            vec![vec![1, 3], vec![1, 4], vec![1, 5]]
        );
    }

    #[test]
    fn signature_fiber_tracks_common_labels_and_pairs() {
        let lpm = PrefixLpm::from_interval_sets(&[vec![1, 2, 3, 4], vec![3, 4, 5]]).unwrap();
        let fiber = SignatureFiber::from_lpm_and_signature(&lpm, &[1, 3, 4, 4]).unwrap();

        assert_eq!(fiber.common(), &[4]);
        assert_eq!(fiber.singleton(), &[1, 3]);
        assert_eq!(fiber.k_sigma().enumerate_bases(), vec![vec![1], vec![3]]);
        assert_eq!(
            fiber.pair_from_choice(&[1]).unwrap(),
            (vec![1, 4], vec![3, 4])
        );
        assert!(fiber.is_self_dual_by_complement());
        assert!(fiber.endpoint_restricted_hall(4, 1).unwrap());
    }

    #[test]
    fn signature_fiber_endpoint_deletion_matches_dummy_singleton_case() {
        let lpm = PrefixLpm::from_intervals(&[(1, 3), (2, 4), (4, 5)]).unwrap();
        let fiber = SignatureFiber::from_lpm_and_signature(&lpm, &[1, 2, 3, 4, 4, 5]).unwrap();

        assert_eq!(fiber.common(), &[4]);
        assert_eq!(fiber.singleton(), &[1, 2, 3, 5]);
        assert!(fiber.is_self_dual_by_complement());

        let restricted = fiber.endpoint_restricted_lpm(5).unwrap();
        assert_eq!(restricted.ground(), &[1, 2, 3]);
        assert_eq!(
            restricted.to_basis_matroid(),
            fiber.k_sigma().to_basis_matroid().delete(5)
        );
    }

    #[test]
    fn signature_fiber_rejects_invalid_signatures() {
        let lpm = PrefixLpm::from_intervals(&[(1, 4), (3, 5)]).unwrap();
        assert!(SignatureFiber::from_lpm_and_signature(&lpm, &[1, 3, 3, 3]).is_err());
        assert!(SignatureFiber::from_lpm_and_signature(&lpm, &[1, 3, 9, 9]).is_err());
        assert!(SignatureFiber::from_lpm_and_signature(&lpm, &[1, 3, 4]).is_err());
    }
}
