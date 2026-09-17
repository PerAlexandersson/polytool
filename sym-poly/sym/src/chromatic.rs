//! Chromatic symmetric functions of graphs.
//!
//! The chromatic symmetric function X_G of a graph G on n vertices is the
//! symmetric function
//!
//! ```text
//! X_G = Σ_{proper colorings κ} x_{κ(1)} x_{κ(2)} ··· x_{κ(n)}
//! ```
//!
//! where the sum is over all proper colorings κ: V(G) → ℤ₊ (no two adjacent
//! vertices receive the same color). This was introduced by Stanley (1995).
//!
//! The **q-chromatic quasisymmetric function** (Shareshian-Wachs, 2016)
//! refines this with a parameter q tracking the number of ascents:
//!
//! ```text
//! X_G(q) = Σ_{proper κ} q^{asc(κ)} x_{κ(1)} x_{κ(2)} ··· x_{κ(n)}
//! ```
//!
//! where asc(κ) = |{(i,j) ∈ E(G) : i < j, κ(i) < κ(j)}|.
//!
//! # Examples
//!
//! ```
//! use combinatoric_core::graph::Graph;
//! use combinatoric_core::partition::Partition;
//! use sym_poly_sym::chromatic_symmetric;
//!
//! // X_{P3} = m_{2,1} + 6 m_{1,1,1}
//! let p3 = Graph::path(3);
//! let x = chromatic_symmetric::<i64>(&p3);
//! let m = x.to_monomial_basis();
//! assert_eq!(m.coefficient(&Partition::new(vec![2, 1])), 1);
//! assert_eq!(m.coefficient(&Partition::new(vec![1, 1, 1])), 6);
//! ```

use std::collections::BTreeMap;

use combinatoric_core::{Graph, Partition, Ring};
use sym_poly_core::UnivariatePolynomial;

use crate::llt::{
    dyck_initial_size, dyck_recursive_e_expansion_from_area, QPolynomial, QSymmetricFunction,
};
use crate::{Basis, SymmetricFunction};

// ---------------------------------------------------------------------------
// Chromatic symmetric function
// ---------------------------------------------------------------------------

/// Compute the chromatic symmetric function `X_G` in the monomial basis.
///
/// This visits each stable set partition of the vertex set exactly once,
/// rather than enumerating colorings separately for every integer partition.
/// If a stable set partition has block-size partition `lambda`, it contributes
/// `prod_s m_s(lambda)!` to the coefficient of `m_lambda`, since blocks of
/// each common size can be assigned independently to the corresponding color
/// variables.  The traversal is exponential in the number of vertices, as an
/// exact chromatic-symmetric-function computation must be in general, but it
/// shares all work between the different monomial types.
pub fn chromatic_symmetric<C: Ring>(g: &Graph) -> SymmetricFunction<C> {
    let terms = stable_set_partition_type_counts(g)
        .into_iter()
        .map(|(partition, count)| {
            let coefficient = count
                .checked_mul(equal_part_permutation_factor(&partition))
                .expect("chromatic symmetric coefficient does not fit in i64");
            (partition, C::from_i64(coefficient))
        })
        .collect();
    SymmetricFunction::from_terms(Basis::Monomial, terms)
}

/// Add a universal vertex to a chromatic symmetric function.
///
/// If `chromatic` is `X_G`, the result is `X_(G * K_1)`, where `*` denotes
/// graph join.  The result is returned in the monomial basis.  More generally,
/// this is the linear operator
///
/// ```text
/// m_lambda -> (multiplicity of 1 in lambda + 1) m_(lambda union 1).
/// ```
///
/// Indeed, a universal vertex must form a singleton block in every stable-set
/// partition.  The factor accounts for the ordinary monomial-basis
/// normalization among equal singleton blocks.
///
/// On the two-column Schur sector this operator is explicitly positive:
///
/// ```text
/// U(s_(2^c,1^a))
///   = (a + 1) sum_(j=0)^c s_(2^j,1^(a + 2(c-j) + 1)).
/// ```
///
/// Write `U = sum_(r>=0) (-1)^r p_(r+1) e_r^perp`.  Conjugating by omega gives
/// `D = sum_(r>=0) p_(r+1) h_r^perp`, while dual Jacobi--Trudi writes the
/// conjugate source as `h_A h_B - h_(A+1) h_(B-1)`, where `A = a+c` and
/// `B = c`.  The two double sums telescope, and Newton's identity leaves
/// `(A-B+1) h_(A+1) h_B`.  The displayed formula is then the Pieri rule,
/// conjugated back by omega.
pub fn chromatic_symmetric_add_universal_vertex<C: Ring>(
    chromatic: &SymmetricFunction<C>,
) -> SymmetricFunction<C> {
    let monomial = chromatic.to_monomial_basis();
    let terms = monomial
        .terms()
        .iter()
        .map(|(shape, coefficient)| {
            let singleton_multiplicity = shape.parts().iter().filter(|&&part| part == 1).count();
            let factor = i64::try_from(singleton_multiplicity + 1)
                .expect("singleton multiplicity does not fit in i64");
            let mut extended_parts = shape.parts().to_vec();
            extended_parts.push(1);
            (
                Partition::from_sorted(extended_parts),
                coefficient.clone() * C::from_i64(factor),
            )
        })
        .collect();
    SymmetricFunction::from_terms(Basis::Monomial, terms)
}

fn stable_set_partition_type_counts(g: &Graph) -> BTreeMap<Partition, i64> {
    fn visit(
        g: &Graph,
        vertex: usize,
        blocks: &mut Vec<Vec<usize>>,
        counts: &mut BTreeMap<Partition, i64>,
    ) {
        if vertex == g.num_vertices() {
            let partition = Partition::new(
                blocks
                    .iter()
                    .map(|block| {
                        u32::try_from(block.len()).expect("stable block size does not fit in u32")
                    })
                    .collect(),
            );
            let count = counts.entry(partition).or_default();
            *count = count
                .checked_add(1)
                .expect("stable set partition count does not fit in i64");
            return;
        }

        for block_index in 0..blocks.len() {
            if blocks[block_index]
                .iter()
                .all(|&other| !g.has_edge(vertex, other))
            {
                blocks[block_index].push(vertex);
                visit(g, vertex + 1, blocks, counts);
                blocks[block_index].pop();
            }
        }

        blocks.push(vec![vertex]);
        visit(g, vertex + 1, blocks, counts);
        blocks.pop();
    }

    let mut counts = BTreeMap::new();
    visit(g, 0, &mut Vec::new(), &mut counts);
    counts
}

fn equal_part_permutation_factor(partition: &Partition) -> i64 {
    let mut factor = 1i64;
    let mut index = 0usize;
    while index < partition.num_parts() {
        let part = partition.part(index);
        let mut multiplicity = 1usize;
        while index + multiplicity < partition.num_parts()
            && partition.part(index + multiplicity) == part
        {
            multiplicity += 1;
        }
        for value in 2..=multiplicity {
            factor = factor
                .checked_mul(i64::try_from(value).expect("part multiplicity does not fit in i64"))
                .expect("equal-part permutation factor does not fit in i64");
        }
        index += multiplicity;
    }
    factor
}

/// Compute the ordered non-proper coloring symmetric function `f_(G, Ω)`.
///
/// Fix a total order `Ω` on the edges of `G`, supplied as `edge_order`.
/// For each non-proper coloring `κ`, let `e = ij` be the first monochromatic
/// edge in that order and set `ν(κ) = κ(i) = κ(j)`. Then
///
/// ```text
/// f_(G, Ω)(x) = Σ_κ (1 / x_{ν(κ)}) ∏_{v ∈ V} x_{κ(v)}
/// ```
///
/// where the sum is over all non-proper colorings.
///
/// This implementation uses the equivalent decomposition
///
/// ```text
/// f_(G, Ω) = Σ_{e ∈ E} X_{G_e},
/// ```
///
/// where `G_e` is obtained by taking the graph of edges earlier than `e` in
/// `Ω` and then identifying the endpoints of `e`.
pub fn first_bad_edge_symmetric<C: Ring>(
    g: &Graph,
    edge_order: &[(usize, usize)],
) -> SymmetricFunction<C> {
    let mut normalized_order: Vec<_> = edge_order.iter().copied().map(normalize_edge).collect();
    normalized_order.sort_unstable();
    normalized_order.dedup();
    assert_eq!(
        normalized_order.len(),
        edge_order.len(),
        "edge_order must not contain duplicates"
    );
    let mut graph_edges = g.edges().to_vec();
    graph_edges.sort_unstable();
    assert_eq!(
        normalized_order.as_slice(),
        graph_edges.as_slice(),
        "edge_order must list each edge of the graph exactly once"
    );

    let mut terms: BTreeMap<Partition, C> = BTreeMap::new();

    for idx in 0..edge_order.len() {
        let prefix_edges: Vec<_> = edge_order[..idx]
            .iter()
            .copied()
            .map(normalize_edge)
            .collect();
        let (u, v) = normalize_edge(edge_order[idx]);
        let prefix_graph = Graph::new(g.num_vertices(), &prefix_edges);
        let ge = prefix_graph.identify_vertices(u, v);
        let contribution = chromatic_symmetric::<C>(&ge);
        for (partition, coeff) in contribution.terms() {
            let entry = terms.entry(partition.clone()).or_insert_with(C::zero);
            *entry = entry.clone() + coeff.clone();
        }
    }

    SymmetricFunction::from_terms(Basis::Monomial, terms)
}

fn normalize_edge((u, v): (usize, usize)) -> (usize, usize) {
    if u < v {
        (u, v)
    } else {
        (v, u)
    }
}

// ---------------------------------------------------------------------------
// q-Chromatic quasisymmetric function
// ---------------------------------------------------------------------------

/// The rectangular powerful-tableau polynomials attached to an ordered graph.
///
/// `by_max_column[j]` is Jim Haglund's `Power[pi,R,j+1]`: the sum of
/// `q^inv(G,T)` over powerful `P`-tableaux of the rectangle `R` for which the
/// largest label lies in column `j+1`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RectangularPowerfulTableauPowers {
    pub by_max_column: Vec<UnivariatePolynomial<i64>>,
    pub powerful_count: usize,
}

/// In-memory cache for the Dyck recursion for chromatic quasisymmetric functions.
///
/// This implements the Alexandersson--Sulzgruber Corollary 6.18 translation of
/// the unicellular LLT Dyck recursion to chromatic quasisymmetric functions:
/// multiplicativity, path-graph initial conditions, and the same three-term and
/// six-term Dyck relations. Values are stored in the elementary basis.
#[derive(Debug, Default, Clone)]
pub struct ChromaticQuasisymmetricDyckCache {
    memo: BTreeMap<Vec<u8>, QSymmetricFunction>,
    path_bases: BTreeMap<usize, QSymmetricFunction>,
}

#[derive(Clone)]
struct NaturalUnitIntervalPoset {
    less: Vec<Vec<bool>>,
}

/// Compute the q-chromatic symmetric function X_G(q) in the monomial basis.
///
/// The coefficient of q^a · m_λ counts proper colorings with type λ and
/// exactly `a` ascending edges (edges (i,j) with i < j and κ(i) < κ(j)).
///
/// Returns a map from (partition, ascent_count) to coefficient.
/// This is represented as a `BTreeMap<Partition, Vec<i64>>` where the Vec
/// is indexed by ascent count.
pub fn q_chromatic_symmetric(g: &Graph) -> BTreeMap<Partition, Vec<i64>> {
    q_chromatic_symmetric_with_ascent_edges(g.num_vertices(), g.edges(), g.edges())
}

/// Compute the q-chromatic symmetric function `X_G(x;q)` in the monomial basis.
///
/// This is the polynomial-valued symmetric-function wrapper around
/// [`q_chromatic_symmetric`]. The coefficient of `m_lambda` is a polynomial in
/// `q` whose `q^a` coefficient counts proper colorings of type `lambda` with
/// `a` ascents.
pub fn q_chromatic_symmetric_function(g: &Graph) -> SymmetricFunction<UnivariatePolynomial<i64>> {
    q_chromatic_terms_to_symmetric_function(q_chromatic_symmetric(g))
}

/// Compute `X_G(x;q)` with a separate ascent orientation, in the monomial basis.
///
/// The `proper_edges` determine proper colorings as an undirected graph. The
/// `ascent_edges` are directed edges `(u, v)` contributing an ascent when
/// `color(u) < color(v)`.
pub fn q_chromatic_symmetric_function_with_ascent_edges(
    n: usize,
    proper_edges: &[(usize, usize)],
    ascent_edges: &[(usize, usize)],
) -> SymmetricFunction<UnivariatePolynomial<i64>> {
    q_chromatic_terms_to_symmetric_function(q_chromatic_symmetric_with_ascent_edges(
        n,
        proper_edges,
        ascent_edges,
    ))
}

/// Elementary-basis expansion of `X_G(x;q)`.
pub fn q_chromatic_elementary_expansion(g: &Graph) -> SymmetricFunction<UnivariatePolynomial<i64>> {
    q_chromatic_symmetric_function(g).to_elementary_basis()
}

/// Coefficient of `e_lambda` in `X_G(x;q)`.
pub fn q_chromatic_elementary_coefficient(
    g: &Graph,
    lambda: &Partition,
) -> UnivariatePolynomial<i64> {
    q_chromatic_elementary_expansion(g).coefficient(lambda)
}

/// Coefficient of `e_(cols^rows)` in `X_G(x;q)`.
pub fn q_chromatic_rectangular_elementary_coefficient(
    g: &Graph,
    rows: usize,
    cols: usize,
) -> Option<UnivariatePolynomial<i64>> {
    if rows == 0 || cols == 0 || rows.checked_mul(cols) != Some(g.num_vertices()) {
        return None;
    }
    let lambda = Partition::new(vec![cols as u32; rows]);
    Some(q_chromatic_elementary_coefficient(g, &lambda))
}

/// q-chromatic symmetric function for a Hessenberg/natural-unit-interval area
/// sequence, in the monomial basis.
pub fn hessenberg_q_chromatic_symmetric_function(
    area: &[u8],
) -> Option<SymmetricFunction<UnivariatePolynomial<i64>>> {
    if !is_area_sequence(area) {
        return None;
    }

    Some(q_chromatic_symmetric_function(&Graph::unit_interval(area)))
}

/// Elementary-basis expansion of the q-chromatic symmetric function for a
/// Hessenberg/natural-unit-interval area sequence.
pub fn hessenberg_q_chromatic_elementary_expansion(
    area: &[u8],
) -> Option<SymmetricFunction<UnivariatePolynomial<i64>>> {
    Some(hessenberg_q_chromatic_symmetric_function(area)?.to_elementary_basis())
}

/// Coefficient of `e_lambda` in the q-chromatic symmetric function attached to
/// a Hessenberg/natural-unit-interval area sequence.
pub fn hessenberg_q_chromatic_elementary_coefficient(
    area: &[u8],
    lambda: &Partition,
) -> Option<UnivariatePolynomial<i64>> {
    Some(hessenberg_q_chromatic_elementary_expansion(area)?.coefficient(lambda))
}

/// Coefficient of `e_(cols^rows)` in the q-chromatic symmetric function
/// attached to a Hessenberg/natural-unit-interval area sequence.
pub fn hessenberg_q_chromatic_rectangular_elementary_coefficient(
    area: &[u8],
    rows: usize,
    cols: usize,
) -> Option<UnivariatePolynomial<i64>> {
    if rows == 0 || cols == 0 || rows.checked_mul(cols) != Some(area.len()) {
        return None;
    }
    let lambda = Partition::new(vec![cols as u32; rows]);
    hessenberg_q_chromatic_elementary_coefficient(area, &lambda)
}

impl ChromaticQuasisymmetricDyckCache {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn len(&self) -> usize {
        self.memo.len()
    }

    pub fn is_empty(&self) -> bool {
        self.memo.is_empty()
    }

    /// Compute the elementary expansion of `X_pi(x;q)` from an area sequence.
    pub fn expansion_from_area_sequence(&mut self, area: &[u8]) -> Option<QSymmetricFunction> {
        let q_plus_one = QPolynomial::new(vec![1, 1]);
        let q = QPolynomial::new(vec![0, 1]);
        let memo = &mut self.memo;
        let path_bases = &mut self.path_bases;
        dyck_recursive_e_expansion_from_area(area, memo, &q_plus_one, &q, &mut |path| {
            chromatic_dyck_base(path, path_bases)
        })
    }
}

/// Recursive elementary-basis expansion of the chromatic quasisymmetric
/// function attached to a Hessenberg/natural-unit-interval area sequence.
pub fn hessenberg_q_chromatic_elementary_expansion_recursive(
    area: &[u8],
) -> Option<QSymmetricFunction> {
    ChromaticQuasisymmetricDyckCache::new().expansion_from_area_sequence(area)
}

/// Recursive coefficient of `e_lambda` in the chromatic quasisymmetric function
/// attached to a Hessenberg/natural-unit-interval area sequence.
pub fn hessenberg_q_chromatic_elementary_coefficient_recursive(
    area: &[u8],
    lambda: &Partition,
) -> Option<UnivariatePolynomial<i64>> {
    Some(hessenberg_q_chromatic_elementary_expansion_recursive(area)?.coefficient(lambda))
}

/// Recursive coefficient of `e_(cols^rows)` in the chromatic quasisymmetric
/// function attached to a Hessenberg/natural-unit-interval area sequence.
pub fn hessenberg_q_chromatic_rectangular_elementary_coefficient_recursive(
    area: &[u8],
    rows: usize,
    cols: usize,
) -> Option<UnivariatePolynomial<i64>> {
    if rows == 0 || cols == 0 || rows.checked_mul(cols) != Some(area.len()) {
        return None;
    }
    let lambda = Partition::new(vec![cols as u32; rows]);
    hessenberg_q_chromatic_elementary_coefficient_recursive(area, &lambda)
}

fn chromatic_dyck_base(
    path: &[u8],
    path_bases: &mut BTreeMap<usize, QSymmetricFunction>,
) -> Option<QSymmetricFunction> {
    if path.is_empty() {
        return Some(SymmetricFunction::basis_element(
            Basis::Elementary,
            Partition::empty(),
        ));
    }
    let n = dyck_initial_size(path)?;
    if let Some(expansion) = path_bases.get(&n) {
        return Some(expansion.clone());
    }
    let expansion = chromatic_initial_path_expansion(n);
    path_bases.insert(n, expansion.clone());
    Some(expansion)
}

fn chromatic_initial_path_expansion(n: usize) -> QSymmetricFunction {
    // Shareshian--Wachs path formula:
    // X_{P_n}(x;q) = sum_{alpha |= n} q^(ell(alpha)-1)
    // [alpha_1]_q [alpha_2-1]_q ... [alpha_ell-1]_q e_{sort(alpha)}.
    let mut terms = BTreeMap::new();
    let mut composition = Vec::new();
    add_chromatic_initial_path_compositions(n as u32, &mut composition, &mut terms);
    SymmetricFunction::from_terms(Basis::Elementary, terms)
}

fn add_chromatic_initial_path_compositions(
    remaining: u32,
    composition: &mut Vec<u32>,
    terms: &mut BTreeMap<Partition, QPolynomial>,
) {
    if remaining == 0 {
        let coefficient = chromatic_initial_path_composition_weight(composition);
        if coefficient.is_zero() {
            return;
        }

        let partition = Partition::new(composition.clone());
        let entry = terms.entry(partition).or_insert_with(QPolynomial::zero);
        *entry = entry.clone() + coefficient;
        return;
    }

    for part in 1..=remaining {
        composition.push(part);
        add_chromatic_initial_path_compositions(remaining - part, composition, terms);
        composition.pop();
    }
}

fn chromatic_initial_path_composition_weight(composition: &[u32]) -> QPolynomial {
    let mut coefficient = QPolynomial::monomial(composition.len().saturating_sub(1), 1);
    for (idx, &part) in composition.iter().enumerate() {
        let q_integer_size = if idx == 0 { part } else { part - 1 };
        if q_integer_size == 0 {
            return QPolynomial::zero();
        }
        coefficient = coefficient * QPolynomial::new(vec![1; q_integer_size as usize]);
    }
    coefficient
}

/// Compute the rectangular powerful-tableau polynomials for an ordered graph.
///
/// The associated poset is the natural-unit-interval order determined by the
/// vertex order: `i <_P j` when `i < j` and `ij` is not an edge. The current
/// backtracking implementation uses a `u128` label mask, so it returns `None`
/// for graphs with more than 128 vertices or for non-matching rectangle sizes.
pub fn rectangular_powerful_tableau_powers(
    graph: &Graph,
    rows: usize,
    cols: usize,
) -> Option<RectangularPowerfulTableauPowers> {
    let n = graph.num_vertices();
    if rows == 0 || cols == 0 || rows.checked_mul(cols) != Some(n) || n > 128 {
        return None;
    }

    let poset = NaturalUnitIntervalPoset::from_graph(graph);
    Some(rectangular_powerful_tableau_powers_bruteforce(
        graph, &poset, rows, cols,
    ))
}

fn rectangular_powerful_tableau_powers_bruteforce(
    graph: &Graph,
    poset: &NaturalUnitIntervalPoset,
    rows: usize,
    cols: usize,
) -> RectangularPowerfulTableauPowers {
    let n = graph.num_vertices();
    let mut by_max_column = vec![UnivariatePolynomial::zero(); cols];
    let mut tableau = vec![usize::MAX; n];
    let mut powerful_count = 0usize;
    enumerate_rectangular_powerful_tableaux(
        graph,
        poset,
        rows,
        cols,
        0,
        0,
        &mut tableau,
        &mut by_max_column,
        &mut powerful_count,
    );

    RectangularPowerfulTableauPowers {
        by_max_column,
        powerful_count,
    }
}

/// Compute rectangular powerful-tableau polynomials for a Hessenberg area
/// sequence.
pub fn hessenberg_rectangular_powerful_tableau_powers(
    area: &[u8],
    rows: usize,
    cols: usize,
) -> Option<RectangularPowerfulTableauPowers> {
    if !is_area_sequence(area) {
        return None;
    }
    rectangular_powerful_tableau_powers(&Graph::unit_interval(area), rows, cols)
}

/// Compute the q-chromatic symmetric function with a separate ascent
/// orientation.
///
/// The `proper_edges` determine which colorings are proper, as an undirected
/// graph.  The `ascent_edges` are directed edges `(u, v)` contributing one
/// ascent when `color(u) < color(v)`.  This covers circular unit arc digraphs,
/// where the proper-coloring graph is the underlying simple graph but the
/// ascent statistic uses the circular orientation.
pub fn q_chromatic_symmetric_with_ascent_edges(
    n: usize,
    proper_edges: &[(usize, usize)],
    ascent_edges: &[(usize, usize)],
) -> BTreeMap<Partition, Vec<i64>> {
    assert!(
        ascent_edges.iter().all(|&(u, v)| u < n && v < n && u != v),
        "ascent edge endpoint out of range"
    );

    let graph = Graph::new(n, proper_edges);
    let partitions = Partition::all_of_size(n as u32);
    let mut result: BTreeMap<Partition, Vec<i64>> = BTreeMap::new();

    for lambda in &partitions {
        let q_counts =
            count_proper_colorings_of_type_q_with_ascent_edges(&graph, ascent_edges, lambda);
        if q_counts.iter().any(|&c| c != 0) {
            result.insert(lambda.clone(), q_counts);
        }
    }

    result
}

fn q_chromatic_terms_to_symmetric_function(
    terms: BTreeMap<Partition, Vec<i64>>,
) -> SymmetricFunction<UnivariatePolynomial<i64>> {
    let terms = terms
        .into_iter()
        .map(|(partition, coeffs)| (partition, UnivariatePolynomial::new(coeffs)))
        .collect();
    SymmetricFunction::from_terms(Basis::Monomial, terms)
}

#[allow(clippy::too_many_arguments)]
fn enumerate_rectangular_powerful_tableaux(
    graph: &Graph,
    poset: &NaturalUnitIntervalPoset,
    rows: usize,
    cols: usize,
    cell: usize,
    used: u128,
    tableau: &mut [usize],
    by_max_column: &mut [UnivariatePolynomial<i64>],
    powerful_count: &mut usize,
) {
    let n = rows * cols;
    if cell == n {
        let (max_col, inv) = max_column_and_inversions(graph, tableau, cols);
        by_max_column[max_col] =
            by_max_column[max_col].clone() + UnivariatePolynomial::monomial(inv, 1);
        *powerful_count += 1;
        return;
    }

    let row = cell / cols;
    let col = cell % cols;
    for label in 0..n {
        let bit = 1u128 << label;
        if used & bit != 0 {
            continue;
        }
        if col > 0 {
            let left = tableau[cell - 1];
            if poset.less(label, left) {
                continue;
            }
        }
        if row > 0 {
            let below = tableau[(row - 1) * cols + col];
            if !poset.less(below, label) {
                continue;
            }
        }

        tableau[cell] = label;
        if col + 1 == cols && !row_has_no_right_to_left_p_minimum(poset, tableau, row, cols) {
            tableau[cell] = usize::MAX;
            continue;
        }
        enumerate_rectangular_powerful_tableaux(
            graph,
            poset,
            rows,
            cols,
            cell + 1,
            used | bit,
            tableau,
            by_max_column,
            powerful_count,
        );
        tableau[cell] = usize::MAX;
    }
}

fn row_has_no_right_to_left_p_minimum(
    poset: &NaturalUnitIntervalPoset,
    tableau: &[usize],
    row: usize,
    cols: usize,
) -> bool {
    for col in 0..cols.saturating_sub(1) {
        let value = tableau[row * cols + col];
        let smaller_than_all_to_right =
            (col + 1..cols).all(|right| poset.less(value, tableau[row * cols + right]));
        if smaller_than_all_to_right {
            return false;
        }
    }
    true
}

fn max_column_and_inversions(graph: &Graph, tableau: &[usize], cols: usize) -> (usize, usize) {
    let n = graph.num_vertices();
    let mut column_of = vec![0usize; n];
    for (cell, &label) in tableau.iter().enumerate() {
        column_of[label] = cell % cols;
    }

    let max_col = column_of[n - 1];
    let inv = graph
        .edges()
        .iter()
        .filter(|&&(u, v)| {
            let (a, b) = if u < v { (u, v) } else { (v, u) };
            column_of[a] > column_of[b]
        })
        .count();
    (max_col, inv)
}

impl NaturalUnitIntervalPoset {
    fn from_graph(graph: &Graph) -> Self {
        let n = graph.num_vertices();
        let mut adjacent = vec![vec![false; n]; n];
        for &(u, v) in graph.edges() {
            adjacent[u][v] = true;
            adjacent[v][u] = true;
        }
        let mut less = vec![vec![false; n]; n];
        for i in 0..n {
            for j in i + 1..n {
                less[i][j] = !adjacent[i][j];
            }
        }
        Self { less }
    }

    fn less(&self, a: usize, b: usize) -> bool {
        self.less[a][b]
    }
}

/// Circular unit interval Frobenius target for the naive dot-action test.
///
/// This returns the coefficient of each `q^d` in
/// `ω X_{\Gamma_a}(x;q)`, where proper colorings use the underlying circular
/// unit interval graph and ascents use the circular edge orientation.
pub fn circular_area_dot_frobenius_target(
    area: &[u8],
) -> Option<BTreeMap<u32, SymmetricFunction<i64>>> {
    let directed_edges = Graph::circular_unit_interval_directed_edges(area)?;
    let graph = Graph::circular_unit_interval(area)?;
    let q_chromatic = q_chromatic_symmetric_with_ascent_edges(
        graph.num_vertices(),
        graph.edges(),
        &directed_edges,
    );
    Some(omega_by_q_degree(q_chromatic))
}

/// Shareshian--Wachs Frobenius target for a Hessenberg area sequence.
///
/// For the unit-interval graph `G` encoded by `area`, this returns the
/// coefficient of each `q^d` in `ω X_G(q)`. Under the Shareshian--Wachs
/// theorem, these are the Frobenius characteristics of the graded dot-action
/// pieces of the corresponding regular semisimple Hessenberg variety.
pub fn hessenberg_area_dot_frobenius_target(
    area: &[u8],
) -> Option<BTreeMap<u32, SymmetricFunction<i64>>> {
    if !is_area_sequence(area) {
        return None;
    }

    let graph = Graph::unit_interval(area);
    let q_chromatic = q_chromatic_symmetric(&graph);
    Some(omega_by_q_degree(q_chromatic))
}

fn omega_by_q_degree(
    q_chromatic: BTreeMap<Partition, Vec<i64>>,
) -> BTreeMap<u32, SymmetricFunction<i64>> {
    let mut by_degree: BTreeMap<u32, BTreeMap<Partition, i64>> = BTreeMap::new();
    for (partition, coeffs) in q_chromatic {
        for (degree, coeff) in coeffs.into_iter().enumerate() {
            if coeff == 0 {
                continue;
            }
            by_degree
                .entry(degree as u32)
                .or_default()
                .insert(partition.clone(), coeff);
        }
    }

    by_degree
        .into_iter()
        .map(|(degree, terms)| {
            let chromatic_degree = SymmetricFunction::from_terms(Basis::Monomial, terms);
            (degree, chromatic_degree.omega_involution())
        })
        .collect()
}

fn is_area_sequence(area: &[u8]) -> bool {
    area.iter().enumerate().all(|(i, &v)| v as usize <= i)
        && area
            .windows(2)
            .all(|w| usize::from(w[1]) <= usize::from(w[0]) + 1)
}

fn next_multiset_perm(perm: &mut [usize]) -> bool {
    if perm.len() <= 1 {
        return false;
    }
    let mut left = perm.len() - 2;
    loop {
        if perm[left] < perm[left + 1] {
            break;
        }
        if left == 0 {
            return false;
        }
        left -= 1;
    }
    let mut right = perm.len() - 1;
    while perm[right] <= perm[left] {
        right -= 1;
    }
    perm.swap(left, right);
    perm[left + 1..].reverse();
    true
}

fn count_proper_colorings_of_type_q_with_ascent_edges(
    g: &Graph,
    ascent_edges: &[(usize, usize)],
    lambda: &Partition,
) -> Vec<i64> {
    let n = g.num_vertices();
    let parts = lambda.parts();
    let max_ascents = ascent_edges.len();

    let mut base_coloring = Vec::with_capacity(n);
    for (color, &freq) in parts.iter().enumerate() {
        for _ in 0..freq {
            base_coloring.push(color);
        }
    }

    let mut counts = vec![0i64; max_ascents + 1];
    let mut perm = base_coloring;
    perm.sort();

    loop {
        let proper = g.edges().iter().all(|&(u, v)| perm[u] != perm[v]);
        if proper {
            let ascents = ascent_edges
                .iter()
                .filter(|&&(u, v)| perm[u] < perm[v])
                .count();
            counts[ascents] += 1;
        }
        if !next_multiset_perm(&mut perm) {
            break;
        }
    }

    while counts.len() > 1 && *counts.last().unwrap() == 0 {
        counts.pop();
    }

    counts
}

/// Format the q-chromatic symmetric function for display.
pub fn format_q_chromatic(qcs: &BTreeMap<Partition, Vec<i64>>) -> String {
    let mut parts: Vec<String> = Vec::new();
    for (lambda, coeffs) in qcs {
        let poly_parts: Vec<String> = coeffs
            .iter()
            .enumerate()
            .filter(|(_, &c)| c != 0)
            .map(|(a, &c)| {
                if a == 0 {
                    format!("{}", c)
                } else if a == 1 {
                    if c == 1 {
                        "q".to_string()
                    } else {
                        format!("{}q", c)
                    }
                } else if c == 1 {
                    format!("q^{}", a)
                } else {
                    format!("{}q^{}", c, a)
                }
            })
            .collect();
        if !poly_parts.is_empty() {
            parts.push(format!(
                "({}) m_{{{}}}",
                poly_parts.join(" + "),
                lambda.display()
            ));
        }
    }
    parts.join(" + ")
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -- Chromatic symmetric function tests (Mathematica-verified) --

    #[test]
    fn test_edge_graph() {
        // K2: colorings (0,1) and (1,0), both → x₁x₂. coeff of m_{1,1} = 2.
        let g = Graph::new(2, &[(0, 1)]);
        let x = chromatic_symmetric::<i64>(&g);
        let m = x.to_monomial_basis();
        assert_eq!(m.coefficient(&Partition::new(vec![1, 1])), 2);
    }

    #[test]
    fn test_p3_monomial() {
        // X_{P3}: direct polynomial = x₁²x₂ + ... + 6x₁x₂x₃
        // = 1·m_{2,1} + 6·m_{1,1,1}
        let x = chromatic_symmetric::<i64>(&Graph::path(3));
        let m = x.to_monomial_basis();
        assert_eq!(m.coefficient(&Partition::new(vec![2, 1])), 1);
        assert_eq!(m.coefficient(&Partition::new(vec![1, 1, 1])), 6);
        assert_eq!(m.coefficient(&Partition::new(vec![3])), 0);
    }

    #[test]
    fn test_add_universal_vertex_matches_graph_join() {
        for graph in [Graph::empty(0), Graph::empty(3), Graph::path(4)] {
            let chromatic = chromatic_symmetric::<i64>(&graph).to_schur_basis();
            let extended = chromatic_symmetric_add_universal_vertex(&chromatic);
            let cone = graph.join(&Graph::complete(1));
            assert_eq!(extended, chromatic_symmetric::<i64>(&cone));
        }
    }

    #[test]
    fn test_add_universal_vertex_two_column_schur_formula() {
        for twos in 0..=4usize {
            for ones in 0..=5usize {
                let mut source_parts = vec![2; twos];
                source_parts.extend(vec![1; ones]);
                let source = Partition::new(source_parts);
                let schur = SymmetricFunction::<i64>::schur_symmetric(source);
                let actual = chromatic_symmetric_add_universal_vertex(&schur).to_schur_basis();

                let coefficient = i64::try_from(ones + 1).unwrap();
                let expected_terms = (0..=twos)
                    .map(|output_twos| {
                        let output_ones = ones + 2 * (twos - output_twos) + 1;
                        let mut output_parts = vec![2; output_twos];
                        output_parts.extend(vec![1; output_ones]);
                        (Partition::new(output_parts), coefficient)
                    })
                    .collect();
                let expected = SymmetricFunction::from_terms(Basis::Schur, expected_terms);
                assert_eq!(actual, expected, "twos={twos}, ones={ones}");
            }
        }
    }

    #[test]
    fn test_p3_trivial_spec() {
        // X_{P3}(1,1,1) should equal chi_{P3}(3) = 3·2·2 = 12
        let x = chromatic_symmetric::<i64>(&Graph::path(3));
        assert_eq!(x.trivial_specialization(3), 12);
    }

    #[test]
    fn test_k3_monomial() {
        // K3: 6 proper 3-colorings, all type {1,1,1}, one monomial x₁x₂x₃.
        // coeff = 6.
        let x = chromatic_symmetric::<i64>(&Graph::complete(3));
        let m = x.to_monomial_basis();
        assert_eq!(m.coefficient(&Partition::new(vec![1, 1, 1])), 6);
        assert_eq!(m.coefficient(&Partition::new(vec![2, 1])), 0);
    }

    #[test]
    fn test_k4_monomial() {
        // K4: 4! = 24 colorings, all type {1,1,1,1}. coeff = 24.
        let x = chromatic_symmetric::<i64>(&Graph::complete(4));
        let m = x.to_monomial_basis();
        assert_eq!(m.coefficient(&Partition::new(vec![1, 1, 1, 1])), 24);
    }

    #[test]
    fn test_p4_trivial_spec() {
        // chi_{P4}(k) = k(k-1)^3
        // chi(4) = 4·27 = 108
        let x = chromatic_symmetric::<i64>(&Graph::path(4));
        assert_eq!(x.trivial_specialization(4), 108);
        assert_eq!(x.trivial_specialization(3), 3 * 8); // 3·2³ = 24
    }

    #[test]
    fn test_p4_schur_nonneg() {
        // Paths are e-positive (Stanley-Stembridge for unit interval graphs),
        // so all Schur coefficients should be non-negative.
        let x = chromatic_symmetric::<i64>(&Graph::path(4));
        let s = x.to_schur_basis();
        for (_lam, &c) in s.terms() {
            assert!(c >= 0, "Negative Schur coefficient in X_{{P4}}");
        }
    }

    #[test]
    fn test_empty_graph_chromatic() {
        // Empty graph on 3: all 3^3 = 27 colorings are proper (k=3 colors).
        // X(1,1,1) = chi(3) = 27.
        let g = Graph::empty(3);
        let x = chromatic_symmetric::<i64>(&g);
        assert_eq!(x.trivial_specialization(3), 27);
    }

    #[test]
    fn test_stable_partition_chromatic_matches_q_specialization() {
        for n in 0..=4 {
            let possible_edges = (0..n)
                .flat_map(|left| ((left + 1)..n).map(move |right| (left, right)))
                .collect::<Vec<_>>();
            for mask in 0..(1usize << possible_edges.len()) {
                let edges = possible_edges
                    .iter()
                    .enumerate()
                    .filter_map(|(index, &edge)| ((mask >> index) & 1 == 1).then_some(edge))
                    .collect::<Vec<_>>();
                let graph = Graph::new(n, &edges);
                let actual = chromatic_symmetric::<i64>(&graph);
                let q_terms = q_chromatic_symmetric(&graph);
                for partition in Partition::all_of_size(n as u32) {
                    let expected = q_terms
                        .get(&partition)
                        .map_or(0, |polynomial| polynomial.iter().sum());
                    assert_eq!(
                        actual.coefficient(&partition),
                        expected,
                        "graph edges={edges:?}, partition={partition}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_first_bad_edge_single_edge() {
        let g = Graph::new(2, &[(0, 1)]);
        let f = first_bad_edge_symmetric::<i64>(&g, &[(0, 1)]);
        let m = f.to_monomial_basis();
        assert_eq!(m.coefficient(&Partition::new(vec![1])), 1);
        assert_eq!(m.terms().len(), 1);
    }

    #[test]
    fn test_first_bad_edge_path3() {
        let g = Graph::path(3);
        let f = first_bad_edge_symmetric::<i64>(&g, &[(0, 1), (1, 2)]);
        let m = f.to_monomial_basis();
        assert_eq!(m.coefficient(&Partition::new(vec![2])), 1);
        assert_eq!(m.coefficient(&Partition::new(vec![1, 1])), 4);
        assert_eq!(m.terms().len(), 2);
    }

    #[test]
    fn test_first_bad_edge_accepts_unsorted_graph_edges() {
        let g = Graph::new(3, &[(1, 2), (0, 1)]);
        let f = first_bad_edge_symmetric::<i64>(&g, &[(0, 1), (1, 2)]);
        let m = f.to_monomial_basis();
        assert_eq!(m.coefficient(&Partition::new(vec![2])), 1);
        assert_eq!(m.coefficient(&Partition::new(vec![1, 1])), 4);
        assert_eq!(m.terms().len(), 2);
    }

    // -- q-Chromatic tests (Mathematica-verified) --

    #[test]
    fn test_q_chromatic_p3() {
        // P3 type {2,1}: all 6 colorings have 1 ascent, rearr=6 → q·m_{2,1}
        // P3 type {1,1,1}: 6 colorings split as 3 with 0 asc + 3 with 1 asc,
        // rearr=1 → (3+3q)·m_{1,1,1} ... wait no, rearr for {1,1,1} with k=3
        // is 3!/(3!) = 1, so we DON'T divide. 6 colorings total.
        // Let me enumerate type {1,1,1}:
        // (0,1,2): asc(0,1)=0<1 yes, asc(1,2)=1<2 yes → 2 asc
        // (0,2,1): asc(0,1)=0<2 yes, asc(1,2)=2>1 no → 1 asc
        // (1,0,2): asc(0,1)=1>0 no, asc(1,2)=0<2 yes → 1 asc
        // (1,2,0): asc(0,1)=1<2 yes, asc(1,2)=2>0 no → 1 asc
        // (2,0,1): asc(0,1)=2>0 no, asc(1,2)=0<1 yes → 1 asc
        // (2,1,0): asc(0,1)=2>1 no, asc(1,2)=1>0 no → 0 asc
        // Counts: asc=0→1, asc=1→4, asc=2→1. Divided by rearr=1: same.
        let qcs = q_chromatic_symmetric(&Graph::path(3));

        let m21 = qcs.get(&Partition::new(vec![2, 1])).unwrap();
        assert_eq!(m21.len(), 2);
        assert_eq!(m21[0], 0);
        assert_eq!(m21[1], 1);

        let m111 = qcs.get(&Partition::new(vec![1, 1, 1])).unwrap();
        assert_eq!(m111[0], 1); // 0 ascents
        assert_eq!(m111[1], 4); // 1 ascent
        assert_eq!(m111[2], 1); // 2 ascents
    }

    #[test]
    fn test_q_chromatic_p3_sum() {
        // Setting q=1 in qX should give the same as X.
        // m_{2,1}: 0+1 = 1. m_{1,1,1}: 1+4+1 = 6. Matches X_{P3}.
        let qcs = q_chromatic_symmetric(&Graph::path(3));
        let m21_total: i64 = qcs[&Partition::new(vec![2, 1])].iter().sum();
        let m111_total: i64 = qcs[&Partition::new(vec![1, 1, 1])].iter().sum();
        assert_eq!(m21_total, 1);
        assert_eq!(m111_total, 6);
    }

    #[test]
    fn test_q_chromatic_symmetric_function_p3() {
        let x = q_chromatic_symmetric_function(&Graph::path(3));

        assert_eq!(
            x.coefficient(&Partition::new(vec![2, 1])),
            UnivariatePolynomial::new(vec![0, 1])
        );
        assert_eq!(
            x.coefficient(&Partition::new(vec![1, 1, 1])),
            UnivariatePolynomial::new(vec![1, 4, 1])
        );
    }

    #[test]
    fn test_haglund_rectangular_example_33() {
        let coefficient =
            hessenberg_q_chromatic_rectangular_elementary_coefficient(&[0, 1, 1, 2, 1, 1], 2, 3)
                .unwrap();

        assert_eq!(coefficient, UnivariatePolynomial::new(vec![0, 0, 1, 1, 1]));
    }

    #[test]
    fn test_haglund_rectangular_powers_example_33() {
        let powers =
            hessenberg_rectangular_powerful_tableau_powers(&[0, 1, 1, 2, 1, 1], 2, 3).unwrap();

        assert_eq!(powers.powerful_count, 5);
        assert_eq!(
            powers.by_max_column[0],
            UnivariatePolynomial::new(vec![0, 0, 0, 0, 1])
        );
        assert_eq!(
            powers.by_max_column[2],
            UnivariatePolynomial::new(vec![0, 0, 1, 1])
        );
    }

    #[test]
    fn test_haglund_rectangular_example_44() {
        let coefficient = hessenberg_q_chromatic_rectangular_elementary_coefficient(
            &[0, 1, 1, 2, 1, 1, 1, 1],
            2,
            4,
        )
        .unwrap();

        assert_eq!(
            coefficient,
            UnivariatePolynomial::new(vec![0, 1, 3, 5, 6, 5, 3, 1])
        );
    }

    #[test]
    fn test_chromatic_initial_path_expansion_matches_direct_small_ranks() {
        for n in 0..=6 {
            assert_eq!(
                chromatic_initial_path_expansion(n),
                q_chromatic_elementary_expansion(&Graph::path(n)),
                "path rank {n}"
            );
        }
    }

    #[test]
    fn test_recursive_hessenberg_q_chromatic_matches_direct_rank_up_to_5() {
        let mut cache = ChromaticQuasisymmetricDyckCache::new();
        for n in 0..=5 {
            for area in all_area_sequences(n) {
                let recursive = cache
                    .expansion_from_area_sequence(&area)
                    .unwrap_or_else(|| panic!("chromatic recursion did not cover area {area:?}"));
                let direct = hessenberg_q_chromatic_elementary_expansion(&area).unwrap();
                assert_eq!(recursive, direct, "area {area:?}");
            }
        }
        assert!(!cache.is_empty());
    }

    #[test]
    fn test_recursive_haglund_rectangular_examples() {
        assert_eq!(
            hessenberg_q_chromatic_rectangular_elementary_coefficient_recursive(
                &[0, 1, 1, 2, 1, 1],
                2,
                3
            )
            .unwrap(),
            UnivariatePolynomial::new(vec![0, 0, 1, 1, 1])
        );
        assert_eq!(
            hessenberg_q_chromatic_rectangular_elementary_coefficient_recursive(
                &[0, 1, 1, 2, 1, 1, 1, 1],
                2,
                4
            )
            .unwrap(),
            UnivariatePolynomial::new(vec![0, 1, 3, 5, 6, 5, 3, 1])
        );
    }

    #[test]
    fn test_recursive_hessenberg_q_chromatic_covers_rank_7() {
        let mut cache = ChromaticQuasisymmetricDyckCache::new();
        for n in 0..=7 {
            for area in all_area_sequences(n) {
                assert!(
                    cache.expansion_from_area_sequence(&area).is_some(),
                    "area {area:?}"
                );
            }
        }
    }

    #[test]
    fn test_hessenberg_area_dot_frobenius_target_edgeless() {
        let target = hessenberg_area_dot_frobenius_target(&[0, 0, 0]).unwrap();

        assert_eq!(target.keys().copied().collect::<Vec<_>>(), vec![0]);
        let schur = target[&0].to_schur_basis();
        assert_eq!(schur.coefficient(&Partition::new(vec![3])), 1);
        assert_eq!(schur.coefficient(&Partition::new(vec![2, 1])), 2);
        assert_eq!(schur.coefficient(&Partition::new(vec![1, 1, 1])), 1);
        assert_eq!(schur.terms().len(), 3);
    }

    #[test]
    fn test_hessenberg_area_dot_frobenius_target_complete_graph_k3() {
        let target = hessenberg_area_dot_frobenius_target(&[0, 1, 2]).unwrap();

        assert_eq!(target.keys().copied().collect::<Vec<_>>(), vec![0, 1, 2, 3]);
        let degree_zero = target[&0].to_schur_basis();
        let degree_one = target[&1].to_schur_basis();
        let degree_two = target[&2].to_schur_basis();
        let degree_three = target[&3].to_schur_basis();

        assert_eq!(degree_zero.coefficient(&Partition::new(vec![3])), 1);
        assert_eq!(degree_zero.terms().len(), 1);
        assert_eq!(degree_one.coefficient(&Partition::new(vec![3])), 2);
        assert_eq!(degree_one.terms().len(), 1);
        assert_eq!(degree_two.coefficient(&Partition::new(vec![3])), 2);
        assert_eq!(degree_two.terms().len(), 1);
        assert_eq!(degree_three.coefficient(&Partition::new(vec![3])), 1);
        assert_eq!(degree_three.terms().len(), 1);
    }

    #[test]
    fn test_hessenberg_area_dot_frobenius_target_rejects_invalid_area() {
        assert!(hessenberg_area_dot_frobenius_target(&[0, 2]).is_none());
    }

    #[test]
    fn test_hessenberg_q_chromatic_rejects_invalid_area() {
        assert!(hessenberg_q_chromatic_symmetric_function(&[0, 2]).is_none());
    }

    #[test]
    fn test_circular_area_target_extends_hessenberg_target() {
        let area = [0, 1, 1];
        assert_eq!(
            circular_area_dot_frobenius_target(&area).unwrap(),
            hessenberg_area_dot_frobenius_target(&area).unwrap()
        );
    }

    #[test]
    fn test_circular_area_target_directed_cycle_s3() {
        let target = circular_area_dot_frobenius_target(&[1, 1, 1]).unwrap();

        assert_eq!(target.keys().copied().collect::<Vec<_>>(), vec![1, 2]);
        let degree_one = target[&1].to_schur_basis();
        let degree_two = target[&2].to_schur_basis();
        assert_eq!(degree_one.coefficient(&Partition::new(vec![3])), 3);
        assert_eq!(degree_one.terms().len(), 1);
        assert_eq!(degree_two.coefficient(&Partition::new(vec![3])), 3);
        assert_eq!(degree_two.terms().len(), 1);
    }

    // -- Edge cases --

    #[test]
    fn test_single_vertex_chromatic() {
        let x = chromatic_symmetric::<i64>(&Graph::empty(1));
        let m = x.to_monomial_basis();
        assert_eq!(m.coefficient(&Partition::new(vec![1])), 1);
    }

    #[test]
    fn test_zero_vertex_chromatic_is_one() {
        let x = chromatic_symmetric::<i64>(&Graph::empty(0));
        let m = x.to_monomial_basis();
        assert_eq!(m.coefficient(&Partition::empty()), 1);
        assert_eq!(m.terms().len(), 1);
        assert_eq!(m.trivial_specialization(5), 1);
    }

    #[test]
    fn test_zero_vertex_q_chromatic_is_one() {
        let qcs = q_chromatic_symmetric(&Graph::empty(0));

        assert_eq!(qcs, BTreeMap::from([(Partition::empty(), vec![1])]));
    }

    fn all_area_sequences(n: usize) -> Vec<Vec<u8>> {
        if n == 0 {
            return vec![Vec::new()];
        }

        let mut result = Vec::new();
        let mut current = vec![0u8; n];
        area_sequences_rec(1, &mut current, &mut result);
        result
    }

    fn area_sequences_rec(index: usize, current: &mut [u8], result: &mut Vec<Vec<u8>>) {
        if index == current.len() {
            result.push(current.to_vec());
            return;
        }

        let max_value = (current[index - 1] + 1).min(index as u8);
        for value in 0..=max_value {
            current[index] = value;
            area_sequences_rec(index + 1, current, result);
        }
    }
}
