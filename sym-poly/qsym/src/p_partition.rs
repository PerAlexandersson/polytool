//! Stanley's (P,w)-partitions and their generating functions in QSym.
//!
//! Given a poset P on {0,...,n-1} (as cover relations) and a labeling w,
//! the (P,w)-partition generating function is:
//!
//!   Γ(P,w) = Σ_{σ ∈ L(P)} F_{Des(σ)}
//!
//! where L(P) is the set of linear extensions and Des(σ) is computed from
//! the labels w(σ_i).  The identity labeling gives the usual natural-label
//! convention.  An anti-natural labeling gives the strict order-preserving
//! enumerator.
//!
//! This module takes raw poset data (cover relations) so it doesn't depend
//! on the Poset type in combinatoric-core. A Poset-aware wrapper can be
//! added there.

use std::collections::{BTreeMap, BTreeSet};

use sym_poly_core::{Composition, Ring};

use crate::basis::QSymBasis;
use crate::qsym_function::QSymFunction;

/// Linear-extension data used in the fundamental expansion of a
/// `(P,w)`-partition generating function.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PPartitionLinearExtension {
    /// The linear extension as a sequence of vertices in the poset.
    pub extension: Vec<usize>,
    /// The descent set `{i : w(extension[i-1]) > w(extension[i])}`.
    pub descent_set: BTreeSet<u32>,
    /// The composition corresponding to `descent_set`.
    pub descent_composition: Composition,
}

/// Compute the (P,w)-partition generating function in the fundamental basis.
///
/// # Arguments
/// * `n` - number of elements (labeled 0..n-1)
/// * `covers` - cover relations (i, j) meaning i <_P j
///
/// For a naturally labeled poset (labels = identity), returns:
///   Γ(P) = Σ_{σ ∈ L(P)} F_{Des(σ)}
pub fn p_partition_generating_function<C: Ring>(
    n: usize,
    covers: &[(usize, usize)],
) -> QSymFunction<C> {
    let labels = (0..n).collect::<Vec<_>>();
    p_partition_generating_function_with_labels(n, covers, &labels)
}

/// Compute the (P,w)-partition generating function in the fundamental basis.
///
/// `labels[v]` is the label w(v). Labels must be distinct. The descent set of
/// a linear extension σ is `{i : labels[σ_i] > labels[σ_{i+1}]}`.
pub fn p_partition_generating_function_with_labels<C: Ring>(
    n: usize,
    covers: &[(usize, usize)],
    labels: &[usize],
) -> QSymFunction<C> {
    assert_eq!(labels.len(), n, "labels must have length n");
    assert_distinct_labels(labels);
    if n == 0 {
        return QSymFunction::basis_element(QSymBasis::Fundamental, Composition::empty());
    }

    let (children, parents) = adjacency_from_covers(n, covers);

    p_partition_from_traversal(n, &children, &parents, labels)
}

/// Return the linear extensions and descent compositions used in
/// [`p_partition_generating_function`].
pub fn p_partition_linear_extensions(
    n: usize,
    covers: &[(usize, usize)],
) -> Vec<PPartitionLinearExtension> {
    let labels = (0..n).collect::<Vec<_>>();
    p_partition_linear_extensions_with_labels(n, covers, &labels)
}

/// Return the linear extensions and descent compositions used in
/// [`p_partition_generating_function_with_labels`].
pub fn p_partition_linear_extensions_with_labels(
    n: usize,
    covers: &[(usize, usize)],
    labels: &[usize],
) -> Vec<PPartitionLinearExtension> {
    assert_eq!(labels.len(), n, "labels must have length n");
    assert_distinct_labels(labels);

    let (children, parents) = adjacency_from_covers(n, covers);
    enumerate_linear_extensions(n, &children, &parents)
        .into_iter()
        .map(|extension| {
            let descent_set = descent_set_with_labels(&extension, labels);
            let descent_composition = Composition::from_descent_set(&descent_set, n as u32);
            PPartitionLinearExtension {
                extension,
                descent_set,
                descent_composition,
            }
        })
        .collect()
}

/// Compute the strict order-preserving enumerator of a poset.
///
/// This is the generating function for maps f satisfying f(x) < f(y) whenever
/// x <_P y.  It is implemented as a (P,w)-partition enumerator for a canonical
/// anti-natural labeling w.
pub fn strict_p_partition_generating_function<C: Ring>(
    n: usize,
    covers: &[(usize, usize)],
) -> QSymFunction<C> {
    if n == 0 {
        return QSymFunction::basis_element(QSymBasis::Fundamental, Composition::empty());
    }

    let (children, parents) = adjacency_from_covers(n, covers);
    let first_extension = first_linear_extension(n, &children, &parents)
        .expect("a finite acyclic poset has a linear extension");
    let labels = anti_natural_labels(n, &first_extension);
    p_partition_from_traversal(n, &children, &parents, &labels)
}

/// Accumulate fundamental-basis coefficients while traversing linear
/// extensions. This retains just one partial extension's descent composition,
/// rather than materializing every completed extension.
fn p_partition_from_traversal<C: Ring>(
    n: usize,
    children: &[Vec<usize>],
    parents: &[Vec<usize>],
    labels: &[usize],
) -> QSymFunction<C> {
    let mut terms: BTreeMap<Composition, C> = BTreeMap::new();
    let mut in_degree: Vec<usize> = parents.iter().map(Vec::len).collect();
    let mut available: BTreeSet<usize> = (0..n).filter(|&v| in_degree[v] == 0).collect();
    let mut completed_parts = Vec::new();

    accumulate_descent_compositions(
        n,
        children,
        &mut in_degree,
        &mut available,
        labels,
        None,
        0,
        0,
        &mut completed_parts,
        &mut terms,
    );

    QSymFunction::from_terms(QSymBasis::Fundamental, terms)
}

#[allow(clippy::too_many_arguments)]
fn accumulate_descent_compositions<C: Ring>(
    n: usize,
    children: &[Vec<usize>],
    in_degree: &mut [usize],
    available: &mut BTreeSet<usize>,
    labels: &[usize],
    previous: Option<usize>,
    depth: usize,
    current_part: u32,
    completed_parts: &mut Vec<u32>,
    terms: &mut BTreeMap<Composition, C>,
) {
    if depth == n {
        let mut parts = completed_parts.clone();
        parts.push(current_part);
        let entry = terms.entry(Composition::new(parts)).or_insert_with(C::zero);
        *entry = entry.clone() + C::one();
        return;
    }

    let choices: Vec<usize> = available.iter().copied().collect();
    for v in choices {
        available.remove(&v);

        let is_descent = previous.is_some_and(|previous| labels[previous] > labels[v]);
        if is_descent {
            completed_parts.push(current_part);
        }
        let next_current_part = if is_descent {
            1
        } else {
            current_part
                .checked_add(1)
                .expect("P-partition degree exceeds u32")
        };

        let mut newly_available = Vec::new();
        for &child in &children[v] {
            in_degree[child] -= 1;
            if in_degree[child] == 0 {
                available.insert(child);
                newly_available.push(child);
            }
        }

        accumulate_descent_compositions(
            n,
            children,
            in_degree,
            available,
            labels,
            Some(v),
            depth + 1,
            next_current_part,
            completed_parts,
            terms,
        );

        for &child in &children[v] {
            in_degree[child] += 1;
        }
        for child in newly_available {
            available.remove(&child);
        }
        if is_descent {
            completed_parts.pop();
        }
        available.insert(v);
    }
}

fn adjacency_from_covers(
    n: usize,
    covers: &[(usize, usize)],
) -> (Vec<Vec<usize>>, Vec<Vec<usize>>) {
    let mut children: Vec<Vec<usize>> = vec![vec![]; n];
    let mut parents: Vec<Vec<usize>> = vec![vec![]; n];
    for &(a, b) in covers {
        children[a].push(b);
        parents[b].push(a);
    }
    (children, parents)
}

fn assert_distinct_labels(labels: &[usize]) {
    let mut seen = BTreeSet::new();
    for &label in labels {
        assert!(seen.insert(label), "labels must be distinct");
    }
}

fn anti_natural_labels(n: usize, extension: &[usize]) -> Vec<usize> {
    let mut labels = vec![0; n];
    for (position, &element) in extension.iter().enumerate() {
        labels[element] = n - 1 - position;
    }
    labels
}

/// Enumerate all linear extensions of a poset.
///
/// A linear extension is a permutation σ of {0,...,n-1} such that
/// σ(i) appears before σ(j) whenever i <_P j. We return σ as a
/// sequence [σ(0), σ(1), ...] giving the order of elements.
fn enumerate_linear_extensions(
    n: usize,
    children: &[Vec<usize>],
    parents: &[Vec<usize>],
) -> Vec<Vec<usize>> {
    // Find minimal elements (those with no parents)
    let mut in_degree: Vec<usize> = parents.iter().map(|p| p.len()).collect();
    let mut available: BTreeSet<usize> = BTreeSet::new();
    for i in 0..n {
        if in_degree[i] == 0 {
            available.insert(i);
        }
    }

    let mut results = Vec::new();
    let mut current = Vec::new();
    backtrack_extensions(
        n,
        children,
        &mut in_degree,
        &mut available,
        &mut current,
        &mut results,
    );
    results
}

/// Return the first extension in the same lexicographic order as
/// [`enumerate_linear_extensions`], without retaining the rest.
fn first_linear_extension(
    n: usize,
    children: &[Vec<usize>],
    parents: &[Vec<usize>],
) -> Option<Vec<usize>> {
    let mut in_degree: Vec<usize> = parents.iter().map(Vec::len).collect();
    let mut available: BTreeSet<usize> = (0..n).filter(|&v| in_degree[v] == 0).collect();
    let mut extension = Vec::with_capacity(n);

    while let Some(v) = available.pop_first() {
        extension.push(v);
        for &child in &children[v] {
            in_degree[child] -= 1;
            if in_degree[child] == 0 {
                available.insert(child);
            }
        }
    }

    (extension.len() == n).then_some(extension)
}

fn backtrack_extensions(
    n: usize,
    children: &[Vec<usize>],
    in_degree: &mut Vec<usize>,
    available: &mut BTreeSet<usize>,
    current: &mut Vec<usize>,
    results: &mut Vec<Vec<usize>>,
) {
    if current.len() == n {
        results.push(current.clone());
        return;
    }

    let choices: Vec<usize> = available.iter().copied().collect();
    for v in choices {
        available.remove(&v);
        current.push(v);

        // Update in-degrees for children of v
        let mut newly_available = Vec::new();
        for &child in &children[v] {
            in_degree[child] -= 1;
            if in_degree[child] == 0 {
                available.insert(child);
                newly_available.push(child);
            }
        }

        backtrack_extensions(n, children, in_degree, available, current, results);

        // Undo
        current.pop();
        available.insert(v);
        for &child in &children[v] {
            in_degree[child] += 1;
        }
        for child in &newly_available {
            available.remove(child);
        }
    }
}

/// Compute the descent composition of a labeled linear extension σ.
///
/// The descent set is {i : labels[σ[i]] > labels[σ[i+1]]}.
#[cfg(test)]
fn descent_composition_with_labels(sigma: &[usize], n: usize, labels: &[usize]) -> Composition {
    let des_set = descent_set_with_labels(sigma, labels);
    Composition::from_descent_set(&des_set, n as u32)
}

fn descent_set_with_labels(sigma: &[usize], labels: &[usize]) -> BTreeSet<u32> {
    let n = sigma.len();
    if n <= 1 {
        return BTreeSet::new();
    }
    let mut des_set = BTreeSet::new();
    for i in 0..n - 1 {
        if labels[sigma[i]] > labels[sigma[i + 1]] {
            des_set.insert((i + 1) as u32);
        }
    }
    des_set
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chain_3() {
        // Chain 0 < 1 < 2: only one linear extension [0, 1, 2], no descents
        // Γ = F_(3)
        let f: QSymFunction<i64> = p_partition_generating_function(3, &[(0, 1), (1, 2)]);
        assert_eq!(f.coefficient(&Composition::new(vec![3])), 1);
        assert_eq!(f.terms().len(), 1);
    }

    #[test]
    fn test_chain_3_strict() {
        // Strict maps on a 3-chain have generating function e_3 = F_(1,1,1).
        let f: QSymFunction<i64> = strict_p_partition_generating_function(3, &[(0, 1), (1, 2)]);
        assert_eq!(f.coefficient(&Composition::new(vec![1, 1, 1])), 1);
        assert_eq!(f.terms().len(), 1);
    }

    #[test]
    fn test_antichain_2() {
        // Antichain on 2 elements: linear extensions are [0,1] and [1,0]
        // [0,1]: no descents → F_(2)
        // [1,0]: descent at position 1 → F_(1,1)
        // Γ = F_(2) + F_(1,1)
        let f: QSymFunction<i64> = p_partition_generating_function(2, &[]);
        assert_eq!(f.coefficient(&Composition::new(vec![2])), 1);
        assert_eq!(f.coefficient(&Composition::new(vec![1, 1])), 1);
    }

    #[test]
    fn test_antichain_3() {
        // Antichain on 3: 3! = 6 linear extensions
        // Γ should have total coefficient sum = 6
        let f: QSymFunction<i64> = p_partition_generating_function(3, &[]);
        let total: i64 = f.terms().values().sum();
        assert_eq!(total, 6);
    }

    #[test]
    fn test_v_poset() {
        // V-poset: 0 < 2, 1 < 2 (two elements below one)
        // Linear extensions: [0,1,2] and [1,0,2]
        // [0,1,2]: no descents → F_(3)
        // [1,0,2]: descent at 1 → F_(1,2)
        let f: QSymFunction<i64> = p_partition_generating_function(3, &[(0, 2), (1, 2)]);
        assert_eq!(f.coefficient(&Composition::new(vec![3])), 1);
        assert_eq!(f.coefficient(&Composition::new(vec![1, 2])), 1);
        assert_eq!(f.terms().len(), 2);

        let data = p_partition_linear_extensions(3, &[(0, 2), (1, 2)]);
        assert_eq!(data.len(), 2);
        assert_eq!(data[0].extension, vec![0, 1, 2]);
        assert_eq!(data[0].descent_set, BTreeSet::new());
        assert_eq!(data[0].descent_composition, Composition::new(vec![3]));
        assert_eq!(data[1].extension, vec![1, 0, 2]);
        assert_eq!(data[1].descent_set, BTreeSet::from([1]));
        assert_eq!(data[1].descent_composition, Composition::new(vec![1, 2]));
    }

    #[test]
    fn test_labeled_v_poset() {
        // With an anti-natural labeling, both cover inequalities are strict.
        // The first extension [0,1,2] has descents at both positions, while
        // [1,0,2] has a descent only at the second position.
        let labels = vec![2, 1, 0];
        let f: QSymFunction<i64> =
            p_partition_generating_function_with_labels(3, &[(0, 2), (1, 2)], &labels);
        assert_eq!(f.coefficient(&Composition::new(vec![1, 1, 1])), 1);
        assert_eq!(f.coefficient(&Composition::new(vec![2, 1])), 1);
        assert_eq!(f.terms().len(), 2);
    }

    #[test]
    fn test_p_partition_is_qsym() {
        // Verify the result is homogeneous of the right degree
        let f: QSymFunction<i64> = p_partition_generating_function(4, &[(0, 2), (1, 3)]);
        for (comp, _) in f.terms() {
            assert_eq!(comp.size(), 4, "all compositions should have degree 4");
        }
    }

    #[test]
    fn test_fence_3() {
        // Fence: 0 < 1, 1 > 2 (= 2 < 1)
        // Cover relations: (0,1), (2,1)
        // Linear extensions: [0,2,1] and [2,0,1]
        // [0,2,1]: 0<2, 2>1 → descent at 2 → F_(2,1)
        // [2,0,1]: 2>0, 0<1 → descent at 1 → F_(1,2)
        let f: QSymFunction<i64> = p_partition_generating_function(3, &[(0, 1), (2, 1)]);
        assert_eq!(f.coefficient(&Composition::new(vec![2, 1])), 1);
        assert_eq!(f.coefficient(&Composition::new(vec![1, 2])), 1);
        assert_eq!(f.terms().len(), 2);
    }

    #[test]
    fn test_streaming_p_partition_matches_materialized_extensions() {
        let cases = [
            (0, vec![], vec![]),
            (3, vec![(0, 2), (1, 2)], vec![2, 0, 1]),
            (4, vec![(0, 2), (1, 2), (1, 3)], vec![1, 3, 0, 2]),
        ];

        for (n, covers, labels) in cases {
            let mut expected = BTreeMap::new();
            for data in p_partition_linear_extensions_with_labels(n, &covers, &labels) {
                *expected.entry(data.descent_composition).or_insert(0i64) += 1;
            }

            let actual: QSymFunction<i64> =
                p_partition_generating_function_with_labels(n, &covers, &labels);
            assert_eq!(actual.basis(), QSymBasis::Fundamental);
            assert_eq!(actual.terms(), &expected);
        }
    }

    #[test]
    fn test_streaming_strict_p_partition_matches_materialized_extensions() {
        let n = 4;
        let covers = vec![(0, 2), (1, 2), (1, 3)];
        let (children, parents) = adjacency_from_covers(n, &covers);
        let extensions = enumerate_linear_extensions(n, &children, &parents);
        let labels = anti_natural_labels(n, extensions.first().unwrap());
        let mut expected = BTreeMap::new();
        for extension in extensions {
            *expected
                .entry(descent_composition_with_labels(&extension, n, &labels))
                .or_insert(0i64) += 1;
        }

        let actual: QSymFunction<i64> = strict_p_partition_generating_function(n, &covers);
        assert_eq!(actual.basis(), QSymBasis::Fundamental);
        assert_eq!(actual.terms(), &expected);
    }
}
