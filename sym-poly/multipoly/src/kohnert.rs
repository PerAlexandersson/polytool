//! Kohnert diagrams, Assaf canonical labelings, and Yamanouchi tests.
//!
//! Coordinates are `(col, row)` with rows indexed from bottom to top.
//! Kohnert moves decrease the row index.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use sym_poly_core::Ring;

use crate::multipoly::MultiPoly;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct Cell {
    pub col: usize,
    pub row: usize,
}

pub type Diagram = BTreeSet<Cell>;
pub type Labeling = BTreeMap<Cell, usize>;

/// A diagram whose cells are split into ordinary cells and fixed ghost cells.
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct GhostDiagram {
    real: Diagram,
    ghosts: Diagram,
}

impl GhostDiagram {
    /// Create a ghost diagram from ordinary cells and ghost cells.
    pub fn new(real: Diagram, ghosts: Diagram) -> Self {
        assert!(
            real.is_disjoint(&ghosts),
            "ordinary cells and ghost cells must be disjoint"
        );
        Self { real, ghosts }
    }

    /// A ghost diagram with no ghost cells.
    pub fn from_real(real: Diagram) -> Self {
        Self::new(real, Diagram::new())
    }

    /// Ordinary cells.
    pub fn real(&self) -> &Diagram {
        &self.real
    }

    /// Ghost cells.
    pub fn ghosts(&self) -> &Diagram {
        &self.ghosts
    }

    /// All occupied positions, ordinary and ghost.
    pub fn occupied(&self) -> Diagram {
        self.real.union(&self.ghosts).copied().collect()
    }

    /// Number of ghost cells.
    pub fn ghost_count(&self) -> usize {
        self.ghosts.len()
    }
}

/// Key diagram of a weak composition, with `alpha[i]` cells in row `i + 1`.
pub fn key_diagram(alpha: &[u32]) -> Diagram {
    let mut diagram = Diagram::new();
    for (idx, &row_len) in alpha.iter().enumerate() {
        let row = idx + 1;
        for col in 1..=row_len as usize {
            diagram.insert(Cell { col, row });
        }
    }
    diagram
}

pub fn rothe_diagram(perm: &[usize]) -> Diagram {
    let mut diagram = Diagram::new();
    for i in 0..perm.len() {
        for j in i + 1..perm.len() {
            if perm[i] > perm[j] {
                diagram.insert(Cell {
                    col: perm[j],
                    row: i + 1,
                });
            }
        }
    }
    diagram
}

pub fn diagram_weight(diagram: &Diagram) -> Vec<u32> {
    let max_row = diagram.iter().map(|cell| cell.row).max().unwrap_or(0);
    let mut weight = vec![0; max_row];
    for cell in diagram {
        weight[cell.row - 1] += 1;
    }
    trim_trailing_zeroes(&mut weight);
    weight
}

/// Weight of a ghost diagram, counting both ordinary and ghost cells by row.
pub fn ghost_diagram_weight(diagram: &GhostDiagram) -> Vec<u32> {
    let occupied = diagram.occupied();
    diagram_weight(&occupied)
}

pub fn max_col(diagram: &Diagram) -> usize {
    diagram.iter().map(|cell| cell.col).max().unwrap_or(0)
}

pub fn cells_in_col(diagram: &Diagram, col: usize) -> Vec<Cell> {
    diagram
        .iter()
        .copied()
        .filter(|cell| cell.col == col)
        .collect()
}

pub fn sorted_rows_in_col(diagram: &Diagram, col: usize) -> Vec<usize> {
    cells_in_col(diagram, col)
        .into_iter()
        .map(|cell| cell.row)
        .collect()
}

pub fn diagram_from_labeling(labeling: &Labeling) -> Diagram {
    labeling.keys().copied().collect()
}

pub fn kohnert_moves(diagram: &Diagram) -> Vec<Diagram> {
    let rows = diagram.iter().map(|cell| cell.row).collect::<BTreeSet<_>>();
    let mut moves = Vec::new();

    for row in rows {
        let Some(rightmost) = diagram
            .iter()
            .filter(|cell| cell.row == row)
            .max_by_key(|cell| cell.col)
            .copied()
        else {
            continue;
        };
        if row == 1 {
            continue;
        }
        let Some(target_row) = (1..row).rev().find(|target_row| {
            !diagram.contains(&Cell {
                col: rightmost.col,
                row: *target_row,
            })
        }) else {
            continue;
        };

        let mut next = diagram.clone();
        next.remove(&rightmost);
        next.insert(Cell {
            col: rightmost.col,
            row: target_row,
        });
        moves.push(next);
    }

    moves
}

pub fn kohnert_diagrams(initial: &Diagram, max_diagrams: usize) -> Result<Vec<Diagram>, String> {
    let mut seen = BTreeSet::new();
    let mut queue = VecDeque::new();
    seen.insert(initial.clone());
    queue.push_back(initial.clone());

    while let Some(diagram) = queue.pop_front() {
        for next in kohnert_moves(&diagram) {
            if seen.insert(next.clone()) {
                if seen.len() > max_diagrams {
                    return Err(format!(
                        "Kohnert diagram cap exceeded: more than {max_diagrams} diagrams"
                    ));
                }
                queue.push_back(next);
            }
        }
    }

    Ok(seen.into_iter().collect())
}

/// Kohnert diagrams generated from the key diagram of a weak composition.
pub fn kohnert_diagrams_for_composition(
    alpha: &[u32],
    max_diagrams: usize,
) -> Result<Vec<Diagram>, String> {
    kohnert_diagrams(&key_diagram(alpha), max_diagrams)
}

/// Count ordinary Kohnert diagrams by row weight.
pub fn kohnert_weight_counts(
    alpha: &[u32],
    max_diagrams: usize,
) -> Result<BTreeMap<Vec<u32>, usize>, String> {
    let mut counts = BTreeMap::new();
    for diagram in kohnert_diagrams_for_composition(alpha, max_diagrams)? {
        *counts.entry(diagram_weight(&diagram)).or_insert(0) += 1;
    }
    Ok(counts)
}

/// Weight enumerator of the Kohnert diagrams generated from an initial diagram.
pub fn kohnert_polynomial<C: Ring>(
    initial: &Diagram,
    max_diagrams: usize,
) -> Result<MultiPoly<C>, String> {
    let num_vars = initial.iter().map(|cell| cell.row).max().unwrap_or(0);
    kohnert_polynomial_with_num_vars(initial, max_diagrams, num_vars)
}

/// Weight enumerator of the Kohnert diagrams generated from the key diagram
/// of a weak composition.
pub fn kohnert_polynomial_for_composition<C: Ring>(
    alpha: &[u32],
    max_diagrams: usize,
) -> Result<MultiPoly<C>, String> {
    let initial = key_diagram(alpha);
    kohnert_polynomial_with_num_vars(&initial, max_diagrams, alpha.len())
}

/// All one-step K-Kohnert moves.
///
/// Rows are indexed from bottom to top, so a move lowers a selected rightmost
/// ordinary cell to the nearest empty position below it in the same column.  A
/// ghost cell is fixed and blocks jumps past it.  Each legal move has two
/// variants: one leaves the source position empty, and one leaves a ghost cell
/// at the source.
pub fn k_kohnert_moves(diagram: &GhostDiagram) -> Vec<GhostDiagram> {
    let occupied = diagram.occupied();
    let rows = diagram
        .real
        .iter()
        .map(|cell| cell.row)
        .collect::<BTreeSet<_>>();
    let mut moves = Vec::new();

    for row in rows {
        let Some(rightmost_occupied) = occupied
            .iter()
            .filter(|cell| cell.row == row)
            .max_by_key(|cell| cell.col)
            .copied()
        else {
            continue;
        };
        if !diagram.real.contains(&rightmost_occupied) || row == 1 {
            continue;
        }

        let mut target = None;
        for target_row in (1..row).rev() {
            let candidate = Cell {
                col: rightmost_occupied.col,
                row: target_row,
            };
            if diagram.ghosts.contains(&candidate) {
                break;
            }
            if !occupied.contains(&candidate) {
                target = Some(candidate);
                break;
            }
        }

        let Some(target) = target else {
            continue;
        };

        let mut moved_real = diagram.real.clone();
        moved_real.remove(&rightmost_occupied);
        moved_real.insert(target);
        moves.push(GhostDiagram::new(
            moved_real.clone(),
            diagram.ghosts.clone(),
        ));

        let mut ghosts = diagram.ghosts.clone();
        ghosts.insert(rightmost_occupied);
        moves.push(GhostDiagram::new(moved_real, ghosts));
    }

    moves
}

/// Closure of a ghost diagram under K-Kohnert moves.
pub fn k_kohnert_diagrams(
    initial: &GhostDiagram,
    max_diagrams: usize,
) -> Result<Vec<GhostDiagram>, String> {
    let mut seen = BTreeSet::new();
    let mut queue = VecDeque::new();
    seen.insert(initial.clone());
    queue.push_back(initial.clone());

    while let Some(diagram) = queue.pop_front() {
        for next in k_kohnert_moves(&diagram) {
            if seen.insert(next.clone()) {
                if seen.len() > max_diagrams {
                    return Err(format!(
                        "K-Kohnert diagram cap exceeded: more than {max_diagrams} diagrams"
                    ));
                }
                queue.push_back(next);
            }
        }
    }

    Ok(seen.into_iter().collect())
}

/// K-Kohnert diagrams generated from the key diagram of a weak composition.
pub fn k_kohnert_diagrams_for_composition(
    alpha: &[u32],
    max_diagrams: usize,
) -> Result<Vec<GhostDiagram>, String> {
    let initial = GhostDiagram::from_real(key_diagram(alpha));
    k_kohnert_diagrams(&initial, max_diagrams)
}

/// Count K-Kohnert diagrams by `(number of ghosts, row weight)`.
pub fn k_kohnert_weight_counts(
    alpha: &[u32],
    max_diagrams: usize,
) -> Result<BTreeMap<(usize, Vec<u32>), usize>, String> {
    let mut counts = BTreeMap::new();
    for diagram in k_kohnert_diagrams_for_composition(alpha, max_diagrams)? {
        *counts
            .entry((diagram.ghost_count(), ghost_diagram_weight(&diagram)))
            .or_insert(0) += 1;
    }
    Ok(counts)
}

/// Lascoux polynomial from K-Kohnert diagrams, with `beta` marking ghost cells.
pub fn lascoux_polynomial<C: Ring>(
    alpha: &[u32],
    beta: &C,
    max_diagrams: usize,
) -> Result<MultiPoly<C>, String> {
    let num_vars = alpha.len();
    let mut terms: BTreeMap<Vec<u32>, C> = BTreeMap::new();

    for diagram in k_kohnert_diagrams_for_composition(alpha, max_diagrams)? {
        let coeff = ring_power(beta, diagram.ghost_count());
        if coeff.is_zero() {
            continue;
        }
        let weight = ghost_diagram_weight_with_num_vars(&diagram, num_vars);
        let entry = terms.entry(weight).or_insert_with(C::zero);
        *entry = entry.clone() + coeff;
    }

    Ok(MultiPoly::from_terms(num_vars, terms))
}

pub(crate) fn kohnert_polynomial_with_num_vars<C: Ring>(
    initial: &Diagram,
    max_diagrams: usize,
    num_vars: usize,
) -> Result<MultiPoly<C>, String> {
    let mut terms: BTreeMap<Vec<u32>, C> = BTreeMap::new();
    for diagram in kohnert_diagrams(initial, max_diagrams)? {
        let weight = diagram_weight_with_num_vars(&diagram, num_vars);
        let entry = terms.entry(weight).or_insert_with(C::zero);
        *entry = entry.clone() + C::one();
    }
    Ok(MultiPoly::from_terms(num_vars, terms))
}

fn diagram_weight_with_num_vars(diagram: &Diagram, num_vars: usize) -> Vec<u32> {
    let mut weight = vec![0; num_vars];
    for cell in diagram {
        assert!(
            cell.row <= num_vars,
            "diagram row exceeds requested number of variables"
        );
        weight[cell.row - 1] += 1;
    }
    weight
}

fn ghost_diagram_weight_with_num_vars(diagram: &GhostDiagram, num_vars: usize) -> Vec<u32> {
    diagram_weight_with_num_vars(&diagram.occupied(), num_vars)
}

fn ring_power<C: Ring>(base: &C, exponent: usize) -> C {
    let mut result = C::one();
    for _ in 0..exponent {
        result = result * base.clone();
    }
    result
}

/// Assaf's column pairing for ordinary, unlabeled rectification.
pub fn column_pairing(diagram: &Diagram, col: usize) -> BTreeMap<Cell, Cell> {
    let left = cells_in_col(diagram, col);
    let right = cells_in_col(diagram, col + 1);
    let mut paired_left = BTreeSet::new();
    let mut paired_right = BTreeSet::new();
    let mut pairs = BTreeMap::new();

    for right_cell in &right {
        let left_cell = Cell {
            col,
            row: right_cell.row,
        };
        if diagram.contains(&left_cell) {
            paired_left.insert(left_cell);
            paired_right.insert(*right_cell);
            pairs.insert(*right_cell, left_cell);
        }
    }

    loop {
        let mut unpaired = left
            .iter()
            .chain(right.iter())
            .copied()
            .filter(|cell| {
                if cell.col == col {
                    !paired_left.contains(cell)
                } else {
                    !paired_right.contains(cell)
                }
            })
            .collect::<Vec<_>>();
        unpaired.sort_by_key(|cell| cell.row);

        let mut changed = false;
        for pair in unpaired.windows(2) {
            let lower = pair[0];
            let upper = pair[1];
            if lower.col == col + 1 && upper.col == col && lower.row < upper.row {
                paired_right.insert(lower);
                paired_left.insert(upper);
                pairs.insert(lower, upper);
                changed = true;
                break;
            }
        }

        if !changed {
            break;
        }
    }

    pairs
}

pub fn label_pairing(labeling: &Labeling, col: usize) -> BTreeMap<Cell, Cell> {
    let diagram = diagram_from_labeling(labeling);
    let mut paired_left = BTreeSet::new();
    let mut pairs = BTreeMap::new();
    let mut right = cells_in_col(&diagram, col + 1);
    right.sort_by(|left, right| right.row.cmp(&left.row));

    for right_cell in right {
        let right_label = labeling[&right_cell];
        let candidate = cells_in_col(&diagram, col)
            .into_iter()
            .filter(|left_cell| {
                !paired_left.contains(left_cell)
                    && left_cell.row >= right_cell.row
                    && labeling[left_cell] <= right_label
            })
            .max_by_key(|left_cell| labeling[left_cell]);
        if let Some(left_cell) = candidate {
            paired_left.insert(left_cell);
            pairs.insert(right_cell, left_cell);
        }
    }

    pairs
}

pub fn rectify_labeled_column_star(labeling: &mut Labeling, col: usize) -> bool {
    let pairs = label_pairing(labeling, col);
    let diagram = diagram_from_labeling(labeling);
    let mut right_cells = cells_in_col(&diagram, col + 1);
    right_cells.sort_by(|left, right| right.row.cmp(&left.row));
    let unpaired = right_cells
        .iter()
        .copied()
        .filter(|cell| !pairs.contains_key(cell))
        .collect::<Vec<_>>();
    let mut labels = labeling.clone();

    for x in &unpaired {
        let x_label = labels[x];
        let swap_with = pairs
            .iter()
            .filter(|(z, y)| z.row > x.row && labeling[y] <= x_label && x_label < labels[z])
            .max_by_key(|(z, _)| labels[z])
            .map(|(z, _)| *z);
        if let Some(z) = swap_with {
            let z_label = labels[&z];
            labels.insert(z, x_label);
            labels.insert(*x, z_label);
        }
    }

    for (z, y) in &pairs {
        labels.insert(*z, labeling[y]);
    }

    let mut next = labels;
    for x in unpaired {
        let label = next
            .remove(&x)
            .expect("unpaired right cell should still have a label");
        let target = Cell { col, row: x.row };
        if next.contains_key(&target) {
            return false;
        }
        next.insert(target, label);
    }

    if *labeling == next {
        false
    } else {
        *labeling = next;
        true
    }
}

pub fn rectify_labeled(labeling: &Labeling, min_col: usize) -> Labeling {
    let mut result = labeling.clone();
    loop {
        let diagram = diagram_from_labeling(&result);
        let max_column = max_col(&diagram);
        let mut changed = false;
        for col in (min_col..max_column).rev() {
            if rectify_labeled_column_star(&mut result, col) {
                changed = true;
                break;
            }
        }
        if !changed {
            break;
        }
    }
    result
}

pub fn canonical_labeling(initial: &Diagram, diagram: &Diagram) -> Option<Labeling> {
    let max_column = max_col(initial).max(max_col(diagram));
    let mut labeling = Labeling::new();

    for col in (1..=max_column).rev() {
        let mut labels = sorted_rows_in_col(initial, col);
        let mut cells = cells_in_col(diagram, col);
        if labels.len() != cells.len() {
            return None;
        }

        let suffix = rectify_labeled(&restrict_labeling(&labeling, col + 1), col + 1);
        labels.sort();
        cells.sort_by_key(|cell| cell.row);

        for label in labels {
            let required_row = suffix
                .iter()
                .find(|(cell, suffix_label)| cell.col == col + 1 && **suffix_label == label)
                .map(|(cell, _)| cell.row);
            let position = cells
                .iter()
                .position(|cell| required_row.map(|row| cell.row >= row).unwrap_or(true))?;
            let cell = cells.remove(position);
            labeling.insert(cell, label);
        }
    }

    Some(labeling)
}

pub fn is_yamanouchi(initial: &Diagram, diagram: &Diagram) -> bool {
    let Some(labeling) = canonical_labeling(initial, diagram) else {
        return false;
    };
    let rectified = rectify_labeled(&labeling, 1);
    rectified.iter().all(|(cell, label)| cell.row == *label)
}

pub fn yamanouchi_diagrams(initial: &Diagram, max_diagrams: usize) -> Result<Vec<Diagram>, String> {
    let diagrams = kohnert_diagrams(initial, max_diagrams)?;
    Ok(diagrams
        .into_iter()
        .filter(|diagram| is_yamanouchi(initial, diagram))
        .collect())
}

pub fn format_diagram(diagram: &Diagram) -> String {
    let max_row = diagram.iter().map(|cell| cell.row).max().unwrap_or(0);
    let max_column = max_col(diagram);
    let mut lines = Vec::new();
    for row in (1..=max_row).rev() {
        let entries = (1..=max_column)
            .map(|col| {
                if diagram.contains(&Cell { col, row }) {
                    "x"
                } else {
                    "."
                }
            })
            .collect::<Vec<_>>()
            .join(" ");
        lines.push(format!("row {row}: {entries}"));
    }
    lines.join(" | ")
}

fn restrict_labeling(labeling: &Labeling, min_col: usize) -> Labeling {
    labeling
        .iter()
        .filter(|(cell, _)| cell.col >= min_col)
        .map(|(cell, label)| (*cell, *label))
        .collect()
}

fn trim_trailing_zeroes(values: &mut Vec<u32>) {
    while values.last() == Some(&0) {
        values.pop();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::key_polynomial::key_polynomial;

    fn weights(diagrams: &[Diagram]) -> BTreeMap<Vec<u32>, usize> {
        let mut result = BTreeMap::new();
        for diagram in diagrams {
            *result.entry(diagram_weight(diagram)).or_insert(0) += 1;
        }
        result
    }

    #[test]
    fn test_key_diagram_for_composition() {
        let diagram = key_diagram(&[0, 2, 1]);
        assert_eq!(
            diagram,
            BTreeSet::from([
                Cell { col: 1, row: 2 },
                Cell { col: 2, row: 2 },
                Cell { col: 1, row: 3 },
            ])
        );
    }

    #[test]
    fn test_kohnert_key_site_example() {
        let diagrams = kohnert_diagrams_for_composition(&[0, 2], 10).unwrap();
        assert_eq!(diagrams.len(), 3);
        assert_eq!(
            kohnert_weight_counts(&[0, 2], 10).unwrap(),
            BTreeMap::from([(vec![0, 2], 1), (vec![1, 1], 1), (vec![2], 1)])
        );
        assert_eq!(
            weights(&diagrams),
            BTreeMap::from([(vec![0, 2], 1), (vec![1, 1], 1), (vec![2], 1)])
        );
    }

    #[test]
    fn test_kohnert_polynomial_matches_key_polynomial() {
        for alpha in &[vec![0, 2], vec![1, 0, 2], vec![2, 1, 0]] {
            let kohnert: MultiPoly<i64> = kohnert_polynomial_for_composition(alpha, 1_000).unwrap();
            let key: MultiPoly<i64> = key_polynomial(alpha);
            assert_eq!(kohnert, key, "Kohnert rule mismatch for {alpha:?}");
        }
    }

    #[test]
    fn test_kohnert_polynomial_from_diagram_uses_fixed_row_count() {
        let initial = key_diagram(&[0, 2]);
        let poly: MultiPoly<i64> = kohnert_polynomial(&initial, 10).unwrap();

        assert_eq!(poly.num_vars(), 2);
        assert_eq!(poly.coefficient(&[0, 2]), 1);
        assert_eq!(poly.coefficient(&[1, 1]), 1);
        assert_eq!(poly.coefficient(&[2, 0]), 1);
    }

    #[test]
    fn test_k_kohnert_lascoux_site_example() {
        let counts = k_kohnert_weight_counts(&[0, 2, 1], 100).unwrap();
        let expected = BTreeMap::from([
            ((0, vec![0, 2, 1]), 1),
            ((0, vec![1, 1, 1]), 1),
            ((0, vec![1, 2]), 1),
            ((0, vec![2, 0, 1]), 1),
            ((0, vec![2, 1]), 1),
            ((1, vec![1, 2, 1]), 2),
            ((1, vec![2, 1, 1]), 2),
            ((1, vec![2, 2]), 1),
            ((2, vec![2, 2, 1]), 1),
        ]);
        assert_eq!(counts, expected);
        assert_eq!(counts.values().sum::<usize>(), 11);
    }

    #[test]
    fn test_lascoux_beta_zero_is_kohnert_polynomial() {
        let alpha = [0, 2, 1];
        let lascoux: MultiPoly<i64> = lascoux_polynomial(&alpha, &0, 100).unwrap();
        let kohnert: MultiPoly<i64> = kohnert_polynomial_for_composition(&alpha, 100).unwrap();

        assert_eq!(lascoux, kohnert);
    }

    #[test]
    fn test_lascoux_site_example_coefficients() {
        let lascoux: MultiPoly<i64> = lascoux_polynomial(&[0, 2, 1], &2, 100).unwrap();

        assert_eq!(lascoux.coefficient(&[0, 2, 1]), 1);
        assert_eq!(lascoux.coefficient(&[1, 2, 1]), 4);
        assert_eq!(lascoux.coefficient(&[2, 1, 1]), 4);
        assert_eq!(lascoux.coefficient(&[2, 2, 0]), 2);
        assert_eq!(lascoux.coefficient(&[2, 2, 1]), 4);
    }

    #[test]
    fn test_yamanouchi_2143() {
        let initial = rothe_diagram(&[2, 1, 4, 3]);
        let yam = yamanouchi_diagrams(&initial, 100).unwrap();
        let mut expected = BTreeMap::new();
        expected.insert(vec![1, 0, 1], 1);
        expected.insert(vec![2], 1);
        assert_eq!(weights(&yam), expected);
    }

    #[test]
    fn test_yamanouchi_321654() {
        let initial = rothe_diagram(&[3, 2, 1, 6, 5, 4]);
        let yam = yamanouchi_diagrams(&initial, 1_000).unwrap();
        assert_eq!(yam.len(), 8);
    }
}
