//! Simple undirected graphs on vertex set {0, 1, ..., n-1}.
//!
//! Provides a [`Graph`] type with adjacency-list representation, standard graph
//! generators (complete, bipartite, path, cycle, grid, tree families,
//! Ferrers board, unit interval), and combinatorial algorithms (matchings,
//! independence sets, acyclic orientations).
//!
//! # Examples
//!
//! ```
//! use combinatoric_core::graph::Graph;
//!
//! let k4 = Graph::complete(4);
//! assert_eq!(k4.num_vertices(), 4);
//! assert_eq!(k4.num_edges(), 6);
//! assert_eq!(k4.matching_polynomial(), vec![1, 6, 3]);
//! ```

use std::collections::{BTreeMap, BTreeSet, HashSet};

use num_bigint::BigInt;

use crate::partition::Partition;

// ---------------------------------------------------------------------------
// Polynomial helpers (Vec<i64> arithmetic for chromatic polynomial)
// ---------------------------------------------------------------------------

/// Subtract two polynomials given as coefficient vectors (ascending degree).
fn vec_poly_sub(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut r = vec![0i64; len];
    for (i, &c) in a.iter().enumerate() {
        r[i] += c;
    }
    for (i, &c) in b.iter().enumerate() {
        r[i] -= c;
    }
    while r.last() == Some(&0) {
        r.pop();
    }
    r
}

// ---------------------------------------------------------------------------
// BigInt polynomial helpers (for tree sink polynomial)
// ---------------------------------------------------------------------------

fn big(n: i64) -> BigInt {
    BigInt::from(n)
}
fn big_zero() -> BigInt {
    BigInt::from(0)
}

fn bpoly_add(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let len = a.len().max(b.len());
    let mut r = vec![big_zero(); len];
    for (i, c) in a.iter().enumerate() {
        r[i] += c;
    }
    for (i, c) in b.iter().enumerate() {
        r[i] += c;
    }
    while r.last() == Some(&big_zero()) {
        r.pop();
    }
    r
}

fn bpoly_mul(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r = vec![big_zero(); a.len() + b.len() - 1];
    for (i, ca) in a.iter().enumerate() {
        for (j, cb) in b.iter().enumerate() {
            r[i + j] += ca * cb;
        }
    }
    while r.last() == Some(&big_zero()) {
        r.pop();
    }
    r
}

fn bpoly_mul_t(p: &[BigInt]) -> Vec<BigInt> {
    if p.is_empty() {
        return vec![];
    }
    let mut r = vec![big_zero(); p.len() + 1];
    for (i, c) in p.iter().enumerate() {
        r[i + 1] = c.clone();
    }
    r
}

fn bpoly_scale(p: &[BigInt], s: &BigInt) -> Vec<BigInt> {
    if *s == big_zero() {
        return vec![];
    }
    let mut r: Vec<BigInt> = p.iter().map(|c| c * s).collect();
    while r.last() == Some(&big_zero()) {
        r.pop();
    }
    r
}

fn bpoly_product(polys: &[Vec<BigInt>]) -> Vec<BigInt> {
    let mut result = vec![big(1)];
    for p in polys {
        result = bpoly_mul(&result, p);
    }
    result
}

fn bpoly_exact_div(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    if b.is_empty() {
        panic!("division by zero polynomial");
    }
    if a.is_empty() {
        return vec![];
    }
    let da = a.len() - 1;
    let db = b.len() - 1;
    if da < db {
        return vec![];
    }
    let mut rem = a.to_vec();
    let mut quot = vec![big_zero(); da - db + 1];
    let lead_b = b.last().unwrap();
    for i in (0..=da - db).rev() {
        let q = &rem[i + db] / lead_b;
        quot[i] = q.clone();
        for (j, bj) in b.iter().enumerate() {
            rem[i + j] -= &q * bj;
        }
    }
    while quot.last() == Some(&big_zero()) {
        quot.pop();
    }
    quot
}

fn factorial_big(n: usize) -> BigInt {
    let mut r = big(1);
    for i in 2..=n {
        r *= big(i as i64);
    }
    r
}

fn reachability_word_count(n: usize) -> usize {
    n.div_ceil(64)
}

fn reachability_has_bit(row: &[u64], idx: usize) -> bool {
    ((row[idx / 64] >> (idx % 64)) & 1) == 1
}

fn reachability_set_bit(row: &mut [u64], idx: usize) {
    row[idx / 64] |= 1u64 << (idx % 64);
}

#[derive(Clone)]
struct ReachabilityMatrix {
    n: usize,
    rows: Vec<Vec<u64>>,
}

impl ReachabilityMatrix {
    fn new(n: usize) -> Self {
        let words = reachability_word_count(n);
        Self {
            n,
            rows: vec![vec![0u64; words]; n],
        }
    }

    fn reaches(&self, from: usize, to: usize) -> bool {
        reachability_has_bit(&self.rows[from], to)
    }

    fn add_arc(&mut self, from: usize, to: usize) -> bool {
        if from == to || self.reaches(to, from) {
            return false;
        }
        if self.reaches(from, to) {
            return true;
        }

        let mut successors = self.rows[to].clone();
        reachability_set_bit(&mut successors, to);

        let predecessors: Vec<usize> = (0..self.n)
            .filter(|&v| v == from || self.reaches(v, from))
            .collect();
        for v in predecessors {
            for (dst, src) in self.rows[v].iter_mut().zip(&successors) {
                *dst |= *src;
            }
        }

        true
    }
}

// ---------------------------------------------------------------------------
// Graph type
// ---------------------------------------------------------------------------

/// A simple undirected graph on vertices `0..n`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Graph {
    n: usize,
    edges: Vec<(usize, usize)>,
    adj: Vec<BTreeSet<usize>>,
}

impl Graph {
    // -- Constructors -------------------------------------------------------

    /// Create a graph on `n` vertices with the given edges.
    ///
    /// Edges are unordered pairs `(u, v)` with `u < v`. Duplicates and
    /// self-loops are silently ignored.
    pub fn new(n: usize, edges: &[(usize, usize)]) -> Self {
        let mut adj = vec![BTreeSet::new(); n];
        let mut deduped = Vec::new();
        let mut seen = BTreeSet::new();

        for &(u, v) in edges {
            if u == v || u >= n || v >= n {
                continue;
            }
            let (a, b) = if u < v { (u, v) } else { (v, u) };
            if seen.insert((a, b)) {
                adj[a].insert(b);
                adj[b].insert(a);
                deduped.push((a, b));
            }
        }

        Graph {
            n,
            edges: deduped,
            adj,
        }
    }

    /// Parse a graph from a graph6 string.
    ///
    /// Graph6 is a compact ASCII format for simple undirected graphs, used by
    /// nauty/Traces and standard graph databases. Each character encodes 6 bits
    /// of the upper triangle of the adjacency matrix.
    ///
    /// See <https://users.cecs.anu.edu.au/~bdm/data/formats.txt>.
    pub fn from_graph6(s: &str) -> Result<Self, String> {
        let bytes: Vec<u8> = s.trim().bytes().collect();
        if bytes.is_empty() {
            return Err("empty graph6 string".into());
        }
        if bytes.iter().any(|&b| !(63..=126).contains(&b)) {
            return Err("graph6 bytes must be printable ASCII in the range 63..=126".into());
        }

        // Decode n (number of vertices)
        let (n, offset) = if bytes[0] == 126 {
            // n >= 63: multi-byte encoding
            if bytes.len() < 4 {
                return Err("truncated graph6 header".into());
            }
            if bytes[1] == 126 {
                // n >= 258048: 8-byte encoding (not supported for practical sizes)
                return Err("graph6 with n >= 258048 not supported".into());
            }
            let n = ((bytes[1] as usize - 63) << 12)
                | ((bytes[2] as usize - 63) << 6)
                | (bytes[3] as usize - 63);
            (n, 4)
        } else {
            ((bytes[0] as usize - 63), 1)
        };

        let required_bits = n
            .checked_mul(n.saturating_sub(1))
            .and_then(|value| value.checked_div(2))
            .ok_or_else(|| "graph6 vertex count is too large".to_string())?;
        let required_bytes = required_bits.div_ceil(6);
        if bytes.len() - offset < required_bytes {
            return Err(format!(
                "truncated graph6 payload: expected {required_bytes} data bytes, found {}",
                bytes.len() - offset
            ));
        }

        // Decode adjacency bits from remaining bytes
        let mut bits = Vec::new();
        for &b in &bytes[offset..] {
            if b < 63 || b > 126 {
                return Err(format!("invalid graph6 byte: {}", b));
            }
            let val = b - 63;
            for k in (0..6).rev() {
                bits.push((val >> k) & 1);
            }
        }

        // Upper triangle: bit index maps to (i, j) with j > i
        let mut edges = Vec::new();
        let mut bit_idx = 0;
        for j in 1..n {
            for i in 0..j {
                if bit_idx < bits.len() && bits[bit_idx] == 1 {
                    edges.push((i, j));
                }
                bit_idx += 1;
            }
        }

        Ok(Graph::new(n, &edges))
    }

    /// Parse all graphs from a graph6-format file (one graph per line).
    ///
    /// Skips blank lines and lines starting with `>` (header lines).
    pub fn all_from_graph6_file(path: &std::path::Path) -> Result<Vec<Self>, String> {
        let content = std::fs::read_to_string(path)
            .map_err(|e| format!("failed to read {}: {}", path.display(), e))?;
        let mut graphs = Vec::new();
        for line in content.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('>') {
                continue;
            }
            graphs.push(Self::from_graph6(line)?);
        }
        Ok(graphs)
    }

    /// Empty graph on `n` vertices (no edges).
    pub fn empty(n: usize) -> Self {
        Graph {
            n,
            edges: Vec::new(),
            adj: vec![BTreeSet::new(); n],
        }
    }

    /// Construct the labeled tree on `0..n` encoded by a Prüfer sequence.
    ///
    /// The sequence must have length `n - 2`, and every entry must be a vertex
    /// in `0..n`. The case `n = 1` is allowed with the empty sequence.
    pub fn tree_from_prufer_sequence(n: usize, sequence: &[usize]) -> Result<Self, String> {
        if n == 0 {
            return Err("a Prüfer tree must have at least one vertex".into());
        }
        if sequence.len() != n.saturating_sub(2) {
            return Err(format!(
                "Prüfer sequence for n={n} must have length {}, got {}",
                n.saturating_sub(2),
                sequence.len()
            ));
        }
        if let Some(&v) = sequence.iter().find(|&&v| v >= n) {
            return Err(format!("Prüfer sequence entry {v} is outside 0..{n}"));
        }
        if n == 1 {
            return Ok(Graph::empty(1));
        }

        let mut degree = vec![1usize; n];
        for &v in sequence {
            degree[v] += 1;
        }

        let mut leaves: BTreeSet<usize> = (0..n).filter(|&v| degree[v] == 1).collect();
        let mut edges = Vec::with_capacity(n - 1);

        for &v in sequence {
            let leaf = *leaves
                .iter()
                .next()
                .expect("Prüfer decoding should always have a leaf");
            leaves.remove(&leaf);
            edges.push((leaf, v));

            degree[leaf] -= 1;
            degree[v] -= 1;
            if degree[v] == 1 {
                leaves.insert(v);
            }
        }

        let remaining: Vec<_> = leaves.into_iter().collect();
        debug_assert_eq!(remaining.len(), 2);
        edges.push((remaining[0], remaining[1]));
        Ok(Graph::new(n, &edges))
    }

    /// Construct the labeled tree encoded by a standard Prüfer sequence.
    ///
    /// This infers `n = sequence.len() + 2`. Use
    /// [`Graph::tree_from_prufer_sequence`] when the one-vertex tree is needed,
    /// since the empty standard sequence encodes the two-vertex tree.
    pub fn from_prufer_sequence(sequence: &[usize]) -> Result<Self, String> {
        Self::tree_from_prufer_sequence(sequence.len() + 2, sequence)
    }

    /// Call `f` on every labeled tree on `0..n`.
    ///
    /// This iterates over all Prüfer sequences, so it produces `n^(n-2)` trees
    /// for `n >= 2`. It is intended for small exhaustive checks.
    pub fn for_each_labeled_tree<F>(n: usize, mut f: F)
    where
        F: FnMut(Self),
    {
        if n == 0 {
            return;
        }
        if n == 1 {
            f(Graph::empty(1));
            return;
        }
        if n == 2 {
            f(Graph::path(2));
            return;
        }

        let sequence_len = n - 2;
        let mut sequence = vec![0usize; sequence_len];
        loop {
            f(Graph::tree_from_prufer_sequence(n, &sequence)
                .expect("generated Prüfer sequence should be valid"));

            let mut pos = sequence_len;
            while pos > 0 {
                pos -= 1;
                sequence[pos] += 1;
                if sequence[pos] < n {
                    break;
                }
                sequence[pos] = 0;
            }
            if pos == 0 && sequence[0] == 0 {
                break;
            }
        }
    }

    /// Return all labeled trees on `0..n`.
    ///
    /// This materializes `n^(n-2)` graphs for `n >= 2`, so prefer
    /// [`Graph::for_each_labeled_tree`] for larger scans.
    pub fn all_labeled_trees(n: usize) -> Vec<Self> {
        let mut trees = Vec::new();
        Graph::for_each_labeled_tree(n, |tree| trees.push(tree));
        trees
    }

    // -- Accessors ----------------------------------------------------------

    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.n
    }

    /// Number of edges.
    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }

    /// Edge list as sorted pairs `(u, v)` with `u < v`.
    pub fn edges(&self) -> &[(usize, usize)] {
        &self.edges
    }

    /// Neighbors of vertex `v`.
    pub fn neighbors(&self, v: usize) -> &BTreeSet<usize> {
        &self.adj[v]
    }

    /// Degree of vertex `v`.
    pub fn degree(&self, v: usize) -> usize {
        self.adj[v].len()
    }

    /// Whether vertices `u` and `v` are adjacent.
    pub fn has_edge(&self, u: usize, v: usize) -> bool {
        if u >= self.n || v >= self.n {
            return false;
        }
        self.adj[u].contains(&v)
    }

    // -- Standard generators ------------------------------------------------

    /// Complete graph K_n.
    pub fn complete(n: usize) -> Self {
        let mut edges = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                edges.push((i, j));
            }
        }
        Graph::new(n, &edges)
    }

    /// Path graph P_n on n vertices: 0—1—2—...—(n-1).
    pub fn path(n: usize) -> Self {
        let edges: Vec<_> = (0..n.saturating_sub(1)).map(|i| (i, i + 1)).collect();
        Graph::new(n, &edges)
    }

    /// Cycle graph C_n on n vertices: 0—1—...—(n-1)—0.
    pub fn cycle(n: usize) -> Self {
        if n < 3 {
            return Graph::path(n);
        }
        let mut edges: Vec<_> = (0..n - 1).map(|i| (i, i + 1)).collect();
        edges.push((0, n - 1));
        Graph::new(n, &edges)
    }

    /// Cartesian product `P_rows x P_cols`, also called the rectangular grid.
    ///
    /// Vertices are numbered row-major: `(r, c)` has index `r * cols + c`.
    pub fn cartesian_product_paths(rows: usize, cols: usize) -> Self {
        let mut edges = Vec::new();
        for r in 0..rows {
            for c in 0..cols {
                let v = r * cols + c;
                if r + 1 < rows {
                    edges.push((v, (r + 1) * cols + c));
                }
                if c + 1 < cols {
                    edges.push((v, r * cols + c + 1));
                }
            }
        }
        Graph::new(rows * cols, &edges)
    }

    /// Rectangular grid graph with `rows` rows and `cols` columns.
    ///
    /// This is an alias for [`Graph::cartesian_product_paths`].
    pub fn grid(rows: usize, cols: usize) -> Self {
        Self::cartesian_product_paths(rows, cols)
    }

    /// Ladder graph `P_rungs x P_2`.
    pub fn ladder(rungs: usize) -> Self {
        Self::cartesian_product_paths(rungs, 2)
    }

    /// Fan graph: a path on `path_vertices` vertices plus one universal apex.
    ///
    /// Vertex `0` is the apex, and vertices `1..=path_vertices` form the path.
    pub fn fan(path_vertices: usize) -> Self {
        let mut edges = Vec::new();
        for i in 1..path_vertices {
            edges.push((i, i + 1));
        }
        for i in 1..=path_vertices {
            edges.push((0, i));
        }
        Graph::new(path_vertices + 1, &edges)
    }

    /// Wheel graph: a cycle on `rim_vertices` vertices plus one universal hub.
    ///
    /// For `rim_vertices < 3`, this returns the corresponding fan graph.
    pub fn wheel(rim_vertices: usize) -> Self {
        if rim_vertices < 3 {
            return Self::fan(rim_vertices);
        }
        let mut edges = Vec::new();
        for i in 0..rim_vertices {
            edges.push((i + 1, ((i + 1) % rim_vertices) + 1));
            edges.push((0, i + 1));
        }
        Graph::new(rim_vertices + 1, &edges)
    }

    /// Star graph with `leaves` leaves and center vertex `0`.
    pub fn star(leaves: usize) -> Self {
        let edges: Vec<_> = (1..=leaves).map(|leaf| (0, leaf)).collect();
        Graph::new(leaves + 1, &edges)
    }

    /// Spider tree with arm lengths `arm_lengths`.
    ///
    /// Vertex `0` is the center. Each arm length is the number of edges in that
    /// arm; zero-length arms add no vertices.
    pub fn spider(arm_lengths: &[usize]) -> Self {
        let n = 1 + arm_lengths.iter().sum::<usize>();
        let mut edges = Vec::new();
        let mut next_vertex = 1;
        for &length in arm_lengths {
            let mut previous = 0;
            for _ in 0..length {
                let current = next_vertex;
                next_vertex += 1;
                edges.push((previous, current));
                previous = current;
            }
        }
        Graph::new(n, &edges)
    }

    /// Spider with `arms` arms, each of length `arm_length`.
    pub fn uniform_spider(arms: usize, arm_length: usize) -> Self {
        Self::spider(&vec![arm_length; arms])
    }

    /// Broom tree: a path with `path_vertices` vertices and extra leaves at
    /// path vertex `0`.
    ///
    /// If `path_vertices` is zero, this returns `Graph::star(leaves)`.
    pub fn broom(path_vertices: usize, leaves: usize) -> Self {
        if path_vertices == 0 {
            return Self::star(leaves);
        }
        let mut edges: Vec<_> = (0..path_vertices.saturating_sub(1))
            .map(|i| (i, i + 1))
            .collect();
        for leaf in path_vertices..path_vertices + leaves {
            edges.push((0, leaf));
        }
        Graph::new(path_vertices + leaves, &edges)
    }

    /// Double-star tree with two adjacent centers `0` and `1`.
    ///
    /// `left_leaves` leaves are attached to vertex `0`, and `right_leaves`
    /// leaves are attached to vertex `1`.
    pub fn double_star(left_leaves: usize, right_leaves: usize) -> Self {
        let mut edges = vec![(0, 1)];
        for leaf in 0..left_leaves {
            edges.push((0, 2 + leaf));
        }
        for leaf in 0..right_leaves {
            edges.push((1, 2 + left_leaves + leaf));
        }
        Graph::new(left_leaves + right_leaves + 2, &edges)
    }

    /// Double-star tree with the same number of leaves on both sides.
    pub fn balanced_double_star(leaves_per_side: usize) -> Self {
        Self::double_star(leaves_per_side, leaves_per_side)
    }

    /// Caterpillar tree with a spine and prescribed leaves at each spine vertex.
    ///
    /// The spine has `leaf_counts.len()` vertices, numbered first. Then
    /// `leaf_counts[i]` leaves are attached to spine vertex `i`.
    pub fn caterpillar(leaf_counts: &[usize]) -> Self {
        let spine_vertices = leaf_counts.len();
        if spine_vertices == 0 {
            return Graph::empty(0);
        }
        let n = spine_vertices + leaf_counts.iter().sum::<usize>();
        let mut edges: Vec<_> = (0..spine_vertices - 1).map(|i| (i, i + 1)).collect();
        let mut next_vertex = spine_vertices;
        for (spine_vertex, &leaves) in leaf_counts.iter().enumerate() {
            for _ in 0..leaves {
                edges.push((spine_vertex, next_vertex));
                next_vertex += 1;
            }
        }
        Graph::new(n, &edges)
    }

    /// Caterpillar with `spine_vertices` spine vertices and the same number of
    /// leaves attached to each spine vertex.
    pub fn uniform_caterpillar(spine_vertices: usize, leaves_per_spine_vertex: usize) -> Self {
        Self::caterpillar(&vec![leaves_per_spine_vertex; spine_vertices])
    }

    /// Complete rooted `branching`-ary tree of height `height`.
    ///
    /// Height is measured in edges from the root to a leaf. Thus height zero
    /// gives the one-vertex tree. If `branching` is zero, this also gives the
    /// one-vertex tree.
    pub fn complete_kary_tree(branching: usize, height: usize) -> Self {
        if branching == 0 || height == 0 {
            return Graph::empty(1);
        }
        let mut edges = Vec::new();
        let mut current_level = vec![0usize];
        let mut next_vertex = 1usize;
        for _ in 0..height {
            let mut next_level = Vec::new();
            for &parent in &current_level {
                for _ in 0..branching {
                    let child = next_vertex;
                    next_vertex += 1;
                    edges.push((parent, child));
                    next_level.push(child);
                }
            }
            current_level = next_level;
        }
        Graph::new(next_vertex, &edges)
    }

    /// Complete binary tree of height `height`.
    pub fn complete_binary_tree(height: usize) -> Self {
        Self::complete_kary_tree(2, height)
    }

    /// Friendship graph: `triangles` triangles sharing a common vertex.
    ///
    /// This is also the windmill graph `Wd(3, triangles)`.
    pub fn friendship(triangles: usize) -> Self {
        let mut edges = Vec::new();
        for i in 0..triangles {
            let a = 2 * i + 1;
            let b = 2 * i + 2;
            edges.push((0, a));
            edges.push((0, b));
            edges.push((a, b));
        }
        Graph::new(2 * triangles + 1, &edges)
    }

    /// Complete bipartite graph K_{a,b}.
    ///
    /// Vertices 0..a are one part, a..(a+b) the other.
    pub fn complete_bipartite(a: usize, b: usize) -> Self {
        let n = a + b;
        let mut edges = Vec::new();
        for i in 0..a {
            for j in a..n {
                edges.push((i, j));
            }
        }
        Graph::new(n, &edges)
    }

    /// Complete multipartite graph K_{sizes[0], sizes[1], ...}.
    ///
    /// Vertices are numbered contiguously within each part.
    pub fn complete_multipartite(sizes: &[usize]) -> Self {
        let n: usize = sizes.iter().sum();
        let mut edges = Vec::new();
        let mut offsets = Vec::with_capacity(sizes.len() + 1);
        offsets.push(0);
        for &s in sizes {
            offsets.push(offsets.last().unwrap() + s);
        }
        for (p, &sp) in sizes.iter().enumerate() {
            for (q, &sq) in sizes.iter().enumerate() {
                if q <= p {
                    continue;
                }
                for i in offsets[p]..offsets[p] + sp {
                    for j in offsets[q]..offsets[q] + sq {
                        edges.push((i, j));
                    }
                }
            }
        }
        Graph::new(n, &edges)
    }

    /// Ferrers board graph (bipartite) from a partition λ.
    ///
    /// Vertices: rows `0..ℓ(λ)` and columns `ℓ(λ)..ℓ(λ)+λ_1`.
    /// Edge (row i, col j) exists iff j < λ_i (the box (i,j) is in the diagram).
    pub fn ferrers_board(lambda: &Partition) -> Self {
        let rows = lambda.num_parts();
        if rows == 0 {
            return Graph::empty(0);
        }
        let cols = lambda.part(0) as usize;
        let n = rows + cols;
        let mut edges = Vec::new();
        for i in 0..rows {
            for j in 0..(lambda.part(i) as usize) {
                edges.push((i, rows + j));
            }
        }
        Graph::new(n, &edges)
    }

    /// Skew Ferrers board graph from partitions λ/μ.
    ///
    /// Same as [`ferrers_board`](Self::ferrers_board) but only includes boxes
    /// in the skew shape λ/μ.
    pub fn ferrers_board_skew(lambda: &Partition, mu: &Partition) -> Self {
        let rows = lambda.num_parts();
        if rows == 0 {
            return Graph::empty(0);
        }
        let cols = lambda.part(0) as usize;
        let n = rows + cols;
        let mut edges = Vec::new();
        for i in 0..rows {
            let start = if i < mu.num_parts() {
                mu.part(i) as usize
            } else {
                0
            };
            for j in start..(lambda.part(i) as usize) {
                edges.push((i, rows + j));
            }
        }
        Graph::new(n, &edges)
    }

    /// Unit interval graph from an area sequence.
    ///
    /// Given `area = [a_1, ..., a_n]`, vertex i is connected to vertex j
    /// (with i < j) iff j - i ≤ a_j. This is the incomparability graph of a
    /// unit interval order.
    pub fn unit_interval(area: &[u8]) -> Self {
        let n = area.len();
        let mut edges = Vec::new();
        for j in 0..n {
            let a = area[j] as usize;
            for gap in 1..=a {
                if gap <= j {
                    edges.push((j - gap, j));
                }
            }
        }
        Graph::new(n, &edges)
    }

    /// Check whether `area` is a circular unit interval area sequence.
    ///
    /// The entries are zero-indexed Rust values for
    /// `(a_1,\dotsc,a_n)`.  A circular area sequence satisfies
    /// `0 <= a_i <= n - 1` and `a_{i+1} <= a_i + 1`, with indices taken
    /// cyclically.
    pub fn is_circular_unit_interval_area_sequence(area: &[u8]) -> bool {
        let n = area.len();
        if n == 0 {
            return true;
        }
        area.iter().all(|&value| usize::from(value) < n)
            && (0..n).all(|i| {
                let next = (i + 1) % n;
                usize::from(area[next]) <= usize::from(area[i]) + 1
            })
    }

    /// Directed edges of the circular unit arc digraph from an area sequence.
    ///
    /// The edge orientation is the one used in the circular chromatic
    /// quasisymmetric function: for each target vertex `i`, add
    /// `(i - gap) -> i`, with indices taken modulo `n`, for
    /// `gap = 1, ..., area[i]`.
    pub fn circular_unit_interval_directed_edges(area: &[u8]) -> Option<Vec<(usize, usize)>> {
        if !Self::is_circular_unit_interval_area_sequence(area) {
            return None;
        }

        let n = area.len();
        let mut edges = Vec::new();
        for target in 0..n {
            for gap in 1..=usize::from(area[target]) {
                let source = (target + n - gap) % n;
                edges.push((source, target));
            }
        }
        edges.sort_unstable();
        edges.dedup();
        Some(edges)
    }

    /// Underlying simple graph of a circular unit arc digraph.
    pub fn circular_unit_interval(area: &[u8]) -> Option<Self> {
        let directed_edges = Self::circular_unit_interval_directed_edges(area)?;
        Some(Graph::new(area.len(), &directed_edges))
    }

    /// Petersen graph (10 vertices, 15 edges, 3-regular).
    pub fn petersen() -> Self {
        let mut edges = Vec::new();
        // Outer cycle: 0-1-2-3-4-0
        for i in 0..5 {
            edges.push((i, (i + 1) % 5));
        }
        // Inner pentagram: 5-7-9-6-8-5
        edges.push((5, 7));
        edges.push((7, 9));
        edges.push((9, 6));
        edges.push((6, 8));
        edges.push((8, 5));
        // Spokes: i -- i+5
        for i in 0..5 {
            edges.push((i, i + 5));
        }
        Graph::new(10, &edges)
    }

    // -- Graph operations ---------------------------------------------------

    /// Complement graph: edge (u,v) iff (u,v) is NOT in self.
    pub fn complement(&self) -> Self {
        let mut edges = Vec::new();
        for i in 0..self.n {
            for j in (i + 1)..self.n {
                if !self.has_edge(i, j) {
                    edges.push((i, j));
                }
            }
        }
        Graph::new(self.n, &edges)
    }

    /// Induced subgraph on vertex set `verts` (relabeled to 0..verts.len()).
    pub fn induced_subgraph(&self, verts: &[usize]) -> Self {
        assert!(
            verts.iter().all(|&vertex| vertex < self.n),
            "induced_subgraph vertex is outside the graph"
        );
        let mut seen = BTreeSet::new();
        let unique_verts = verts
            .iter()
            .copied()
            .filter(|vertex| seen.insert(*vertex))
            .collect::<Vec<_>>();
        let vset: BTreeSet<usize> = unique_verts.iter().copied().collect();
        let mut idx = vec![0usize; self.n];
        for (new, &old) in unique_verts.iter().enumerate() {
            idx[old] = new;
        }
        let edges: Vec<_> = self
            .edges
            .iter()
            .filter(|&&(u, v)| vset.contains(&u) && vset.contains(&v))
            .map(|&(u, v)| (idx[u], idx[v]))
            .collect();
        Graph::new(unique_verts.len(), &edges)
    }

    /// Delete a vertex (reindex remaining vertices).
    pub fn delete_vertex(&self, v: usize) -> Self {
        let verts: Vec<usize> = (0..self.n).filter(|&u| u != v).collect();
        self.induced_subgraph(&verts)
    }

    /// Delete an edge.
    pub fn delete_edge(&self, u: usize, v: usize) -> Self {
        let (a, b) = if u < v { (u, v) } else { (v, u) };
        let edges: Vec<_> = self
            .edges
            .iter()
            .copied()
            .filter(|&e| e != (a, b))
            .collect();
        Graph::new(self.n, &edges)
    }

    /// Identify two vertices, removing self-loops and relabeling the remaining
    /// vertices to `0..n-1`.
    pub fn identify_vertices(&self, u: usize, v: usize) -> Self {
        assert!(u < self.n, "vertex u out of range");
        assert!(v < self.n, "vertex v out of range");
        if u == v {
            return self.clone();
        }

        // Redirect all edges involving v to u, then delete v
        let (keep, remove) = if u < v { (u, v) } else { (v, u) };
        let mut new_edges = Vec::new();
        for &(a, b) in &self.edges {
            let a2 = if a == remove { keep } else { a };
            let b2 = if b == remove { keep } else { b };
            if a2 != b2 {
                new_edges.push((a2, b2));
            }
        }
        // Reindex: shift vertices above `remove` down by 1
        let new_edges: Vec<_> = new_edges
            .iter()
            .map(|&(a, b)| {
                let a2 = if a > remove { a - 1 } else { a };
                let b2 = if b > remove { b - 1 } else { b };
                (a2, b2)
            })
            .collect();
        Graph::new(self.n - 1, &new_edges)
    }

    /// Contract an existing edge: merge its endpoints, remove self-loops, and
    /// relabel the remaining vertices to `0..n-1`.
    pub fn contract_edge(&self, u: usize, v: usize) -> Self {
        assert!(u < self.n, "vertex u out of range");
        assert!(v < self.n, "vertex v out of range");
        assert!(u != v, "cannot contract a loop");
        assert!(self.has_edge(u, v), "cannot contract a non-edge");
        self.identify_vertices(u, v)
    }

    /// Line graph L(G): vertices are edges of G, two vertices in L(G) are
    /// adjacent iff the corresponding edges in G share an endpoint.
    ///
    /// Vertex i of L(G) corresponds to `self.edges()[i]`.
    pub fn line_graph(&self) -> Self {
        let m = self.edges.len();
        let mut lg_edges = Vec::new();
        for i in 0..m {
            let (a, b) = self.edges[i];
            for j in (i + 1)..m {
                let (c, d) = self.edges[j];
                if a == c || a == d || b == c || b == d {
                    lg_edges.push((i, j));
                }
            }
        }
        Graph::new(m, &lg_edges)
    }

    /// Whether this graph is the line graph of a simple undirected graph.
    ///
    /// This uses Krausz's characterization: the edges can be partitioned into
    /// cliques so that every vertex belongs to at most two partition cliques.
    /// It searches the maximal cliques of size at least three together with
    /// all two-vertex edge cliques.  Those candidates are sufficient because,
    /// in the star-clique partition of a line graph, every nonmaximal clique
    /// has size two.
    ///
    /// The maximal-clique enumeration is worst-case exponential, although the
    /// initial claw-free check and the exact-cover search make the method
    /// practical for the small and medium combinatorial graphs used here.
    pub fn is_line_graph(&self) -> bool {
        if !self.is_claw_free() {
            return false;
        }
        if self.edges.is_empty() {
            return true;
        }

        let edge_indices: BTreeMap<_, _> = self
            .edges
            .iter()
            .copied()
            .enumerate()
            .map(|(index, edge)| (edge, index))
            .collect();

        let mut candidates = Vec::new();
        for vertices in maximal_cliques(self) {
            if vertices.len() < 3 {
                continue;
            }
            let mut edges = Vec::with_capacity(vertices.len() * (vertices.len() - 1) / 2);
            for (offset, &u) in vertices.iter().enumerate() {
                for &v in &vertices[(offset + 1)..] {
                    let edge = if u < v { (u, v) } else { (v, u) };
                    edges.push(edge_indices[&edge]);
                }
            }
            candidates.push(KrauszClique { vertices, edges });
        }
        for (edge_index, &(u, v)) in self.edges.iter().enumerate() {
            candidates.push(KrauszClique {
                vertices: vec![u, v],
                edges: vec![edge_index],
            });
        }

        let mut candidates_by_edge = vec![Vec::new(); self.edges.len()];
        for (candidate_index, candidate) in candidates.iter().enumerate() {
            for &edge_index in &candidate.edges {
                candidates_by_edge[edge_index].push(candidate_index);
            }
        }
        for edge_candidates in &mut candidates_by_edge {
            edge_candidates
                .sort_unstable_by_key(|&index| std::cmp::Reverse(candidates[index].edges.len()));
        }

        let mut uncovered = vec![true; self.edges.len()];
        let mut memberships = vec![0u8; self.n];
        let mut rejected_states = HashSet::new();
        krausz_partition_exists(
            &candidates,
            &candidates_by_edge,
            &mut uncovered,
            self.edges.len(),
            &mut memberships,
            &mut rejected_states,
        )
    }

    // -- Predicates ---------------------------------------------------------

    /// Is the graph connected?
    pub fn is_connected(&self) -> bool {
        if self.n <= 1 {
            return true;
        }
        let mut visited = vec![false; self.n];
        let mut stack = vec![0usize];
        visited[0] = true;
        let mut count = 1;
        while let Some(v) = stack.pop() {
            for &u in &self.adj[v] {
                if !visited[u] {
                    visited[u] = true;
                    count += 1;
                    stack.push(u);
                }
            }
        }
        count == self.n
    }

    /// Is the graph bipartite? Returns Some((A, B)) if yes, None if no.
    pub fn bipartition(&self) -> Option<(Vec<usize>, Vec<usize>)> {
        let mut color = vec![None; self.n];
        let mut a = Vec::new();
        let mut b = Vec::new();

        for start in 0..self.n {
            if color[start].is_some() {
                continue;
            }
            color[start] = Some(false);
            a.push(start);
            let mut stack = vec![start];
            while let Some(v) = stack.pop() {
                let c = color[v].unwrap();
                for &u in &self.adj[v] {
                    match color[u] {
                        Some(cu) if cu == c => return None,
                        Some(_) => {}
                        None => {
                            color[u] = Some(!c);
                            if !c {
                                b.push(u);
                            } else {
                                a.push(u);
                            }
                            stack.push(u);
                        }
                    }
                }
            }
        }
        Some((a, b))
    }

    /// Is the graph bipartite?
    pub fn is_bipartite(&self) -> bool {
        self.bipartition().is_some()
    }

    /// Is the graph claw-free (K_{1,3}-free)?
    ///
    /// A graph is claw-free iff no vertex has three mutually non-adjacent
    /// neighbors.
    pub fn is_claw_free(&self) -> bool {
        for v in 0..self.n {
            let nbrs: Vec<usize> = self.adj[v].iter().copied().collect();
            if nbrs.len() < 3 {
                continue;
            }
            // Check all triples of neighbors
            for i in 0..nbrs.len() {
                for j in (i + 1)..nbrs.len() {
                    for k in (j + 1)..nbrs.len() {
                        // If none of the three pairs are adjacent, we have a claw
                        if !self.has_edge(nbrs[i], nbrs[j])
                            && !self.has_edge(nbrs[i], nbrs[k])
                            && !self.has_edge(nbrs[j], nbrs[k])
                        {
                            return false;
                        }
                    }
                }
            }
        }
        true
    }

    // -- Matchings ----------------------------------------------------------

    /// All matchings of the graph, as lists of edge indices.
    ///
    /// A matching is a set of edges with no shared endpoints.
    /// Returns matchings grouped by nothing — just a flat list.
    pub fn all_matchings(&self) -> Vec<Vec<(usize, usize)>> {
        let mut result = Vec::new();
        let mut current = Vec::new();
        let mut used = vec![false; self.n];
        self.matchings_rec(0, &mut current, &mut used, &mut result);
        result
    }

    fn matchings_rec(
        &self,
        edge_idx: usize,
        current: &mut Vec<(usize, usize)>,
        used: &mut [bool],
        result: &mut Vec<Vec<(usize, usize)>>,
    ) {
        result.push(current.clone());
        for i in edge_idx..self.edges.len() {
            let (u, v) = self.edges[i];
            if !used[u] && !used[v] {
                used[u] = true;
                used[v] = true;
                current.push((u, v));
                self.matchings_rec(i + 1, current, used, result);
                current.pop();
                used[u] = false;
                used[v] = false;
            }
        }
    }

    /// Matching polynomial: coefficients[k] = number of matchings with k edges.
    pub fn matching_polynomial(&self) -> Vec<i64> {
        let matchings = self.all_matchings();
        let max_k = matchings.iter().map(|m| m.len()).max().unwrap_or(0);
        let mut coeffs = vec![0i64; max_k + 1];
        for m in &matchings {
            coeffs[m.len()] += 1;
        }
        coeffs
    }

    /// All perfect matchings (matchings that cover every vertex).
    pub fn perfect_matchings(&self) -> Vec<Vec<(usize, usize)>> {
        if self.n % 2 != 0 {
            return vec![]; // odd number of vertices → no perfect matching
        }
        let target = self.n / 2;
        self.all_matchings()
            .into_iter()
            .filter(|m| m.len() == target)
            .collect()
    }

    /// All non-crossing matchings.
    ///
    /// A matching is non-crossing if no two edges (a,b) and (c,d) satisfy
    /// a < c < b < d (when endpoints are ordered on a line 0,1,...,n-1).
    pub fn non_crossing_matchings(&self) -> Vec<Vec<(usize, usize)>> {
        self.all_matchings()
            .into_iter()
            .filter(|m| {
                for i in 0..m.len() {
                    for j in (i + 1)..m.len() {
                        let (a, b) = m[i]; // a < b guaranteed by edge ordering
                        let (c, d) = m[j];
                        if (a < c && c < b && b < d) || (c < a && a < d && d < b) {
                            return false;
                        }
                    }
                }
                true
            })
            .collect()
    }

    /// All non-nesting matchings.
    ///
    /// A matching is non-nesting if no two edges (a,b) and (c,d) satisfy
    /// a < c < d < b (one edge nested inside the other).
    pub fn non_nesting_matchings(&self) -> Vec<Vec<(usize, usize)>> {
        self.all_matchings()
            .into_iter()
            .filter(|m| {
                for i in 0..m.len() {
                    for j in (i + 1)..m.len() {
                        let (a, b) = m[i];
                        let (c, d) = m[j];
                        if (a < c && d < b) || (c < a && b < d) {
                            return false;
                        }
                    }
                }
                true
            })
            .collect()
    }

    // -- Triangles and cliques ----------------------------------------------

    /// All triangles (3-cliques) in the graph.
    ///
    /// Returns triples (a, b, c) with a < b < c.
    pub fn triangles(&self) -> Vec<(usize, usize, usize)> {
        let mut result = Vec::new();
        for &(a, b) in &self.edges {
            for &c in self.adj[b].range((b + 1)..) {
                if self.has_edge(a, c) {
                    result.push((a, b, c));
                }
            }
        }
        result
    }

    /// Visit all orientations of the graph without materializing them up front.
    ///
    /// This is usually the right entry point for computing orientation
    /// statistics. [`all_orientations`](Self::all_orientations) is retained as a
    /// convenience wrapper for callers that truly need a `Vec`.
    ///
    /// TODO: If external callers need standard iterator combinators, add a
    /// dedicated owned iterator type. The visitor form keeps the recursive
    /// search state simple and avoids self-referential iterator machinery.
    pub fn for_each_orientation<F>(&self, mut visit: F)
    where
        F: FnMut(&[(usize, usize)]),
    {
        self.for_each_orientation_with_state(&mut |orientation, _| visit(orientation));
    }

    /// All orientations of the graph (all 2^|E| ways to direct edges).
    ///
    /// Returns each orientation as a list of directed edges (u, v) meaning u → v.
    pub fn all_orientations(&self) -> Vec<Vec<(usize, usize)>> {
        let mut result = Vec::new();
        self.for_each_orientation(|orientation| result.push(orientation.to_vec()));
        result
    }

    /// Stirling graph S(n): bipartite graph on {0,...,n-1} ∪ {n,...,2n-1}
    /// where left vertex i is adjacent to right vertex n+j iff i < j.
    ///
    /// The number of k-matchings equals the Stirling number S(n, n-k).
    pub fn stirling(n: usize) -> Self {
        let total = 2 * n;
        let mut edges = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                edges.push((i, n + j));
            }
        }
        Graph::new(total, &edges)
    }

    // -- Independent sets ---------------------------------------------------

    /// All independent sets (including the empty set).
    pub fn all_independent_sets(&self) -> Vec<Vec<usize>> {
        let mut result = Vec::new();
        let mut current = Vec::new();
        self.indep_rec(0, &mut current, &mut result);
        result
    }

    fn indep_rec(&self, start: usize, current: &mut Vec<usize>, result: &mut Vec<Vec<usize>>) {
        result.push(current.clone());
        for v in start..self.n {
            if current.iter().all(|&u| !self.has_edge(u, v)) {
                current.push(v);
                self.indep_rec(v + 1, current, result);
                current.pop();
            }
        }
    }

    /// Independence polynomial: coefficients[k] = number of independent sets of size k.
    pub fn independence_polynomial(&self) -> Vec<i64> {
        let sets = self.all_independent_sets();
        let max_k = sets.iter().map(|s| s.len()).max().unwrap_or(0);
        let mut coeffs = vec![0i64; max_k + 1];
        for s in &sets {
            coeffs[s.len()] += 1;
        }
        coeffs
    }

    // -- Proper colorings ---------------------------------------------------

    /// All proper colorings with `k` colors (labeled 0..k).
    ///
    /// A proper coloring assigns a color to each vertex such that no two
    /// adjacent vertices share the same color.
    pub fn proper_colorings(&self, k: usize) -> Vec<Vec<usize>> {
        let mut result = Vec::new();
        let mut coloring = vec![0usize; self.n];
        self.color_rec(0, k, &mut coloring, &mut result);
        result
    }

    fn color_rec(
        &self,
        v: usize,
        k: usize,
        coloring: &mut Vec<usize>,
        result: &mut Vec<Vec<usize>>,
    ) {
        if v == self.n {
            result.push(coloring.clone());
            return;
        }
        for c in 0..k {
            if self.adj[v].iter().all(|&u| u >= v || coloring[u] != c) {
                coloring[v] = c;
                self.color_rec(v + 1, k, coloring, result);
            }
        }
    }

    /// Chromatic polynomial evaluated at t = k: number of proper k-colorings.
    pub fn chromatic_polynomial_eval(&self, k: usize) -> usize {
        self.proper_colorings(k).len()
    }

    /// Count proper colorings by color type (sorted color multiplicities).
    ///
    /// Returns a map from partition (the color type) to count.
    /// The sum over all partitions of n gives the chromatic polynomial at k = n.
    pub fn proper_colorings_by_type(&self, k: usize) -> Vec<(Vec<u32>, usize)> {
        let colorings = self.proper_colorings(k);
        let mut type_counts: std::collections::BTreeMap<Vec<u32>, usize> =
            std::collections::BTreeMap::new();

        for coloring in &colorings {
            let mut freq = vec![0u32; k];
            for &c in coloring {
                freq[c] += 1;
            }
            freq.sort_unstable_by(|a, b| b.cmp(a));
            while freq.last() == Some(&0) {
                freq.pop();
            }
            *type_counts.entry(freq).or_insert(0) += 1;
        }

        type_counts.into_iter().collect()
    }

    // -- Acyclic orientations -----------------------------------------------

    /// All acyclic orientations of the graph.
    ///
    /// An acyclic orientation assigns a direction to each edge such that the
    /// resulting directed graph has no directed cycles. Each proper coloring
    /// with n colors induces an acyclic orientation (orient u→v if color(u) < color(v)).
    ///
    /// Returns orientations as lists of directed edges `(u, v)` meaning u → v.
    pub fn acyclic_orientations(&self) -> Vec<Vec<(usize, usize)>> {
        self.acyclic_orientations_with_frozen_edges(&[])
    }

    /// Visit all acyclic orientations of the graph.
    ///
    /// Prefer this over [`acyclic_orientations`](Self::acyclic_orientations)
    /// when you want to accumulate statistics rather than store every
    /// orientation.
    pub fn for_each_acyclic_orientation<F>(&self, visit: F)
    where
        F: FnMut(&[(usize, usize)]),
    {
        self.for_each_acyclic_orientation_with_frozen_edges(&[], visit);
    }

    /// All acyclic orientations extending the directed edges in `frozen_edges`.
    ///
    /// Each entry `(u, v)` in `frozen_edges` forces the undirected edge `{u, v}`
    /// to be oriented as `u → v`. Repeated constraints with the same direction
    /// are allowed; conflicting constraints panic. If the frozen directions
    /// themselves create a directed cycle, the result is empty.
    pub fn acyclic_orientations_with_frozen_edges(
        &self,
        frozen_edges: &[(usize, usize)],
    ) -> Vec<Vec<(usize, usize)>> {
        let mut result = Vec::new();
        self.for_each_acyclic_orientation_with_frozen_edges(frozen_edges, |orientation| {
            result.push(orientation.to_vec());
        });
        result
    }

    /// Visit all acyclic orientations extending the directed edges in
    /// `frozen_edges`.
    pub fn for_each_acyclic_orientation_with_frozen_edges<F>(
        &self,
        frozen_edges: &[(usize, usize)],
        mut visit: F,
    ) where
        F: FnMut(&[(usize, usize)]),
    {
        let frozen = self.frozen_edge_orientations(frozen_edges);
        let Some((reachability, mut orientation, mut has_outgoing)) =
            self.initial_acyclic_orientation_state(&frozen)
        else {
            return;
        };

        self.for_each_acyclic_orientation_with_state(
            0,
            &frozen,
            &reachability,
            &mut orientation,
            &mut has_outgoing,
            &mut |orientation, _| visit(orientation),
        );
    }

    /// Number of acyclic orientations.
    pub fn num_acyclic_orientations(&self) -> usize {
        self.num_acyclic_orientations_with_frozen_edges(&[])
    }

    /// Number of acyclic orientations extending the directions in `frozen_edges`.
    pub fn num_acyclic_orientations_with_frozen_edges(
        &self,
        frozen_edges: &[(usize, usize)],
    ) -> usize {
        let mut count = 0usize;
        self.for_each_acyclic_orientation_with_frozen_edges(frozen_edges, |_| count += 1);
        count
    }

    /// Sink polynomial of acyclic orientations: coefficients[k] = number of
    /// acyclic orientations with exactly k sinks.
    ///
    /// A sink is a vertex with no outgoing edges in the orientation.
    /// This equals the chromatic polynomial evaluated at -t (up to sign).
    pub fn acyclic_sink_polynomial(&self) -> Vec<i64> {
        self.acyclic_sink_polynomial_with_frozen_edges(&[])
    }

    /// Sink polynomial of acyclic orientations extending the directions in
    /// `frozen_edges`.
    pub fn acyclic_sink_polynomial_with_frozen_edges(
        &self,
        frozen_edges: &[(usize, usize)],
    ) -> Vec<i64> {
        let frozen = self.frozen_edge_orientations(frozen_edges);
        let Some((reachability, mut orientation, mut has_outgoing)) =
            self.initial_acyclic_orientation_state(&frozen)
        else {
            return vec![0];
        };

        let mut coeffs = vec![0i64; self.n + 1];
        self.for_each_acyclic_orientation_with_state(
            0,
            &frozen,
            &reachability,
            &mut orientation,
            &mut has_outgoing,
            &mut |_, has_outgoing| {
                let sinks = has_outgoing.iter().filter(|&&out| !out).count();
                coeffs[sinks] += 1;
            },
        );
        while coeffs.len() > 1 && coeffs.last() == Some(&0) {
            coeffs.pop();
        }
        coeffs
    }

    // -- Chromatic polynomial (deletion-contraction) -----------------------

    /// Chromatic polynomial via deletion-contraction.
    ///
    /// Returns coefficients in ascending degree order:
    /// `result[k]` is the coefficient of t^k in χ_G(t).
    ///
    /// Uses memoization keyed on canonical edge sets.
    pub fn chromatic_polynomial(&self) -> Vec<i64> {
        let mut cache: std::collections::HashMap<(usize, Vec<(usize, usize)>), Vec<i64>> =
            std::collections::HashMap::new();
        self.chromatic_poly_dc(&mut cache)
    }

    fn chromatic_poly_dc(
        &self,
        cache: &mut std::collections::HashMap<(usize, Vec<(usize, usize)>), Vec<i64>>,
    ) -> Vec<i64> {
        let key = (self.n, self.edges.clone());
        if let Some(cached) = cache.get(&key) {
            return cached.clone();
        }

        let result = if self.edges.is_empty() {
            // No edges: χ(t) = t^n
            let mut c = vec![0i64; self.n + 1];
            c[self.n] = 1;
            c
        } else {
            // χ_G(t) = χ_{G-e}(t) - χ_{G/e}(t)
            let (u, v) = self.edges[0];
            let g_del = self.delete_edge(u, v);
            let g_con = self.contract_edge(u, v);
            let p_del = g_del.chromatic_poly_dc(cache);
            let p_con = g_con.chromatic_poly_dc(cache);
            vec_poly_sub(&p_del, &p_con)
        };

        cache.insert(key, result.clone());
        result
    }

    /// Sink polynomial over ALL orientations (not just acyclic ones).
    ///
    /// For each of the 2^|E| orientations, count sinks (vertices with out-degree 0).
    /// Note: orientations with cycles may have zero sinks, contributing to coeff[0].
    pub fn sink_polynomial_all_orientations(&self) -> Vec<i64> {
        let mut coeffs = vec![0i64; self.n + 1];
        self.for_each_orientation_with_state(&mut |_, has_outgoing| {
            let sinks = has_outgoing.iter().filter(|&&out| !out).count();
            coeffs[sinks] += 1;
        });

        // Strip trailing zeros
        while coeffs.last() == Some(&0) {
            coeffs.pop();
        }
        coeffs
    }

    // -- Tree-specific sink polynomial ----------------------------------------

    /// Sink polynomial of L(T) for a tree T, via bottom-up recursion.
    ///
    /// For a tree, an acyclic orientation of L(T) corresponds to choosing
    /// a total order on the edges incident to each vertex. An edge e={u,v}
    /// is a sink iff e is minimum at u and minimum at v.
    ///
    /// We root T and process bottom-up. At each vertex v with parent edge e_p
    /// and child edges e_1,...,e_k, we choose a permutation of all d_v edges
    /// incident to v. We track a two-variable polynomial in (t, s) where
    /// t marks sinks in the subtree, and s marks whether e_p is minimum at v.
    ///
    /// This runs in O(n · D^2) where D is the max degree, via polynomial
    /// convolution at each vertex.
    ///
    /// Panics if the graph is not a tree.
    pub fn sink_polynomial_tree(&self) -> Vec<i64> {
        use num_traits::ToPrimitive;
        self.sink_polynomial_tree_bigint()
            .iter()
            .map(|c| {
                c.to_i64()
                    .expect("coefficient overflow in sink_polynomial_tree")
            })
            .collect()
    }

    /// Return the per-vertex (A_v, B_v) polynomials from the tree recursion.
    ///
    /// Returns (order, parent, a_polys, b_polys) where:
    /// - order[i] is the vertex in BFS order (root = order[0])
    /// - parent[v] is the parent of v (usize::MAX for root)
    /// - a_polys[v] = A_v(t), b_polys[v] = B_v(t)
    pub fn tree_ab_polynomials(
        &self,
    ) -> (Vec<usize>, Vec<usize>, Vec<Vec<BigInt>>, Vec<Vec<BigInt>>) {
        assert!(
            self.edges.len() + 1 == self.n && self.n > 0 && self.is_connected(),
            "tree_ab_polynomials requires a tree"
        );

        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); self.n];
        for &(u, v) in &self.edges {
            adj[u].push(v);
            adj[v].push(u);
        }
        let root = 0;
        let mut parent = vec![usize::MAX; self.n];
        let mut order = Vec::with_capacity(self.n);
        let mut visited = vec![false; self.n];
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(root);
        visited[root] = true;
        while let Some(v) = queue.pop_front() {
            order.push(v);
            for &u in &adj[v] {
                if !visited[u] {
                    visited[u] = true;
                    parent[u] = v;
                    queue.push_back(u);
                }
            }
        }

        let one_b: Vec<BigInt> = vec![big(1)];
        let zero_b: Vec<BigInt> = vec![];
        let mut a_poly: Vec<Vec<BigInt>> = vec![vec![]; self.n];
        let mut b_poly: Vec<Vec<BigInt>> = vec![vec![]; self.n];

        for &v in order.iter().rev() {
            let children: Vec<usize> = adj[v].iter().copied().filter(|&u| parent[u] == v).collect();
            let k = children.len();

            if v == root {
                // Root: A and B don't apply, but store the sink poly in B for convenience
                if k == 0 {
                    a_poly[v] = one_b.clone();
                    b_poly[v] = zero_b.clone();
                    continue;
                }
                let mut tps: Vec<Vec<BigInt>> = Vec::new();
                for &c in &children {
                    tps.push(bpoly_add(&a_poly[c], &b_poly[c]));
                }
                let p = bpoly_product(&tps);
                let fact = factorial_big(k - 1);
                let mut result = zero_b.clone();
                for (j, &c) in children.iter().enumerate() {
                    let p_j = bpoly_exact_div(&p, &tps[j]);
                    let t_a_j = bpoly_mul_t(&a_poly[c]);
                    let sink_or_not = bpoly_add(&t_a_j, &b_poly[c]);
                    let term = bpoly_mul(&sink_or_not, &p_j);
                    result = bpoly_add(&result, &term);
                }
                a_poly[v] = zero_b.clone();
                b_poly[v] = bpoly_scale(&result, &fact); // store S_{L(T)} in b_poly[root]
            } else if k == 0 {
                a_poly[v] = one_b.clone();
                b_poly[v] = zero_b.clone();
            } else {
                let d = k + 1;
                let mut tps: Vec<Vec<BigInt>> = Vec::new();
                for &c in &children {
                    tps.push(bpoly_add(&a_poly[c], &b_poly[c]));
                }
                let p = bpoly_product(&tps);
                let fact = factorial_big(d - 1);
                a_poly[v] = bpoly_scale(&p, &fact);
                let mut bv = zero_b.clone();
                for (j, &c) in children.iter().enumerate() {
                    let p_j = bpoly_exact_div(&p, &tps[j]);
                    let t_a_j = bpoly_mul_t(&a_poly[c]);
                    let sink_or_not = bpoly_add(&t_a_j, &b_poly[c]);
                    let term = bpoly_mul(&sink_or_not, &p_j);
                    bv = bpoly_add(&bv, &term);
                }
                b_poly[v] = bpoly_scale(&bv, &fact);
            }
        }

        (order, parent, a_poly, b_poly)
    }

    /// Sink polynomial of L(T) for a tree T, using BigInt arithmetic.
    ///
    /// Root T and process bottom-up. At each vertex v with parent edge e_p
    /// and children c_1,...,c_k, track two polynomials:
    ///   A_v(t) = contribution when e_p IS minimum at v,
    ///   B_v(t) = contribution when e_p is NOT minimum at v.
    /// See \cref{prop:treeRecursion} in the paper for details.
    pub fn sink_polynomial_tree_bigint(&self) -> Vec<BigInt> {
        assert!(
            self.edges.len() + 1 == self.n && self.n > 0 && self.is_connected(),
            "sink_polynomial_tree requires a connected tree"
        );

        if self.n == 1 {
            return vec![big(1)];
        }

        // Build adjacency and root at vertex 0
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); self.n];
        for &(u, v) in &self.edges {
            adj[u].push(v);
            adj[v].push(u);
        }
        let root = 0;
        let mut parent = vec![usize::MAX; self.n];
        let mut order = Vec::with_capacity(self.n);
        let mut visited = vec![false; self.n];
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(root);
        visited[root] = true;
        while let Some(v) = queue.pop_front() {
            order.push(v);
            for &u in &adj[v] {
                if !visited[u] {
                    visited[u] = true;
                    parent[u] = v;
                    queue.push_back(u);
                }
            }
        }

        let one_b: Vec<BigInt> = vec![big(1)];
        let zero_b: Vec<BigInt> = vec![];

        let mut a_poly: Vec<Vec<BigInt>> = vec![vec![]; self.n];
        let mut b_poly: Vec<Vec<BigInt>> = vec![vec![]; self.n];

        for &v in order.iter().rev() {
            let children: Vec<usize> = adj[v].iter().copied().filter(|&u| parent[u] == v).collect();
            let k = children.len();

            if v == root {
                if k == 0 {
                    continue;
                }
                let mut tps: Vec<Vec<BigInt>> = Vec::new();
                for &c in &children {
                    tps.push(bpoly_add(&a_poly[c], &b_poly[c]));
                }
                let p = bpoly_product(&tps);
                let fact = factorial_big(k - 1);
                let mut result = zero_b.clone();
                for (j, &c) in children.iter().enumerate() {
                    let p_j = bpoly_exact_div(&p, &tps[j]);
                    let t_a_j = bpoly_mul_t(&a_poly[c]);
                    let sink_or_not = bpoly_add(&t_a_j, &b_poly[c]);
                    let term = bpoly_mul(&sink_or_not, &p_j);
                    result = bpoly_add(&result, &term);
                }
                return bpoly_scale(&result, &fact);
            }

            if k == 0 {
                a_poly[v] = one_b.clone();
                b_poly[v] = zero_b.clone();
            } else {
                let d = k + 1;
                let mut tps: Vec<Vec<BigInt>> = Vec::new();
                for &c in &children {
                    tps.push(bpoly_add(&a_poly[c], &b_poly[c]));
                }
                let p = bpoly_product(&tps);
                let fact = factorial_big(d - 1);

                a_poly[v] = bpoly_scale(&p, &fact);

                let mut bv = zero_b.clone();
                for (j, &c) in children.iter().enumerate() {
                    let p_j = bpoly_exact_div(&p, &tps[j]);
                    let t_a_j = bpoly_mul_t(&a_poly[c]);
                    let sink_or_not = bpoly_add(&t_a_j, &b_poly[c]);
                    let term = bpoly_mul(&sink_or_not, &p_j);
                    bv = bpoly_add(&bv, &term);
                }
                b_poly[v] = bpoly_scale(&bv, &fact);
            }
        }

        one_b
    }

    // -- Display ------------------------------------------------------------

    /// Format as edge list string.
    pub fn display_edges(&self) -> String {
        let edges: Vec<String> = self
            .edges
            .iter()
            .map(|(u, v)| format!("{}-{}", u, v))
            .collect();
        format!("Graph({} vertices, edges: [{}])", self.n, edges.join(", "))
    }

    fn frozen_edge_orientations(
        &self,
        frozen_edges: &[(usize, usize)],
    ) -> Vec<Option<(usize, usize)>> {
        let edge_index: std::collections::BTreeMap<(usize, usize), usize> = self
            .edges
            .iter()
            .copied()
            .enumerate()
            .map(|(i, edge)| (edge, i))
            .collect();

        let mut frozen = vec![None; self.edges.len()];
        for &(from, to) in frozen_edges {
            assert!(
                from < self.n && to < self.n,
                "frozen edge ({from}, {to}) uses a vertex outside 0..{}",
                self.n
            );
            assert_ne!(from, to, "frozen edge ({from}, {to}) is a loop");

            let key = if from < to { (from, to) } else { (to, from) };
            let Some(&edge_idx) = edge_index.get(&key) else {
                panic!("frozen edge ({from}, {to}) is not an edge of the graph");
            };

            match frozen[edge_idx] {
                Some((u, v)) => assert!(
                    (u, v) == (from, to),
                    "conflicting frozen directions for edge {{{}, {}}}",
                    key.0,
                    key.1
                ),
                None => frozen[edge_idx] = Some((from, to)),
            }
        }

        frozen
    }

    fn initial_acyclic_orientation_state(
        &self,
        frozen: &[Option<(usize, usize)>],
    ) -> Option<(ReachabilityMatrix, Vec<(usize, usize)>, Vec<bool>)> {
        let mut reachability = ReachabilityMatrix::new(self.n);
        let mut orientation = vec![(0usize, 0usize); self.edges.len()];
        let mut has_outgoing = vec![false; self.n];

        for (edge_idx, frozen_dir) in frozen.iter().enumerate() {
            if let Some((from, to)) = frozen_dir {
                if !reachability.add_arc(*from, *to) {
                    return None;
                }
                orientation[edge_idx] = (*from, *to);
                has_outgoing[*from] = true;
            }
        }

        Some((reachability, orientation, has_outgoing))
    }

    fn for_each_orientation_with_state<F>(&self, visit: &mut F)
    where
        F: FnMut(&[(usize, usize)], &[bool]),
    {
        let mut orientation = vec![(0usize, 0usize); self.edges.len()];
        let mut has_outgoing = vec![false; self.n];
        self.for_each_orientation_rec(0, &mut orientation, &mut has_outgoing, visit);
    }

    fn for_each_orientation_rec<F>(
        &self,
        edge_idx: usize,
        orientation: &mut Vec<(usize, usize)>,
        has_outgoing: &mut Vec<bool>,
        visit: &mut F,
    ) where
        F: FnMut(&[(usize, usize)], &[bool]),
    {
        if edge_idx == self.edges.len() {
            visit(orientation, has_outgoing);
            return;
        }

        let (u, v) = self.edges[edge_idx];
        for (from, to) in [(u, v), (v, u)] {
            orientation[edge_idx] = (from, to);
            let was_outgoing = has_outgoing[from];
            has_outgoing[from] = true;
            self.for_each_orientation_rec(edge_idx + 1, orientation, has_outgoing, visit);
            has_outgoing[from] = was_outgoing;
        }
    }

    fn for_each_acyclic_orientation_with_state<F>(
        &self,
        edge_idx: usize,
        frozen: &[Option<(usize, usize)>],
        reachability: &ReachabilityMatrix,
        orientation: &mut Vec<(usize, usize)>,
        has_outgoing: &mut Vec<bool>,
        visit: &mut F,
    ) where
        F: FnMut(&[(usize, usize)], &[bool]),
    {
        if edge_idx == self.edges.len() {
            visit(orientation, has_outgoing);
            return;
        }

        if frozen[edge_idx].is_some() {
            self.for_each_acyclic_orientation_with_state(
                edge_idx + 1,
                frozen,
                reachability,
                orientation,
                has_outgoing,
                visit,
            );
            return;
        }

        let (u, v) = self.edges[edge_idx];
        for (from, to) in [(u, v), (v, u)] {
            let mut next_reachability = reachability.clone();
            if next_reachability.add_arc(from, to) {
                orientation[edge_idx] = (from, to);
                let was_outgoing = has_outgoing[from];
                has_outgoing[from] = true;
                self.for_each_acyclic_orientation_with_state(
                    edge_idx + 1,
                    frozen,
                    &next_reachability,
                    orientation,
                    has_outgoing,
                    visit,
                );
                has_outgoing[from] = was_outgoing;
            }
        }
    }
}

#[derive(Debug)]
struct KrauszClique {
    vertices: Vec<usize>,
    edges: Vec<usize>,
}

fn maximal_cliques(graph: &Graph) -> Vec<Vec<usize>> {
    fn visit(
        graph: &Graph,
        clique: &mut Vec<usize>,
        mut candidates: BTreeSet<usize>,
        mut excluded: BTreeSet<usize>,
        result: &mut Vec<Vec<usize>>,
    ) {
        if candidates.is_empty() && excluded.is_empty() {
            result.push(clique.clone());
            return;
        }

        let pivot = candidates
            .iter()
            .chain(excluded.iter())
            .copied()
            .max_by_key(|&vertex| candidates.intersection(graph.neighbors(vertex)).count());
        let extension_vertices: Vec<_> = match pivot {
            Some(vertex) => candidates
                .difference(graph.neighbors(vertex))
                .copied()
                .collect(),
            None => candidates.iter().copied().collect(),
        };

        for vertex in extension_vertices {
            clique.push(vertex);
            let next_candidates = candidates
                .intersection(graph.neighbors(vertex))
                .copied()
                .collect();
            let next_excluded = excluded
                .intersection(graph.neighbors(vertex))
                .copied()
                .collect();
            visit(graph, clique, next_candidates, next_excluded, result);
            clique.pop();
            candidates.remove(&vertex);
            excluded.insert(vertex);
        }
    }

    let mut result = Vec::new();
    visit(
        graph,
        &mut Vec::new(),
        (0..graph.n).collect(),
        BTreeSet::new(),
        &mut result,
    );
    result
}

fn krausz_partition_exists(
    candidates: &[KrauszClique],
    candidates_by_edge: &[Vec<usize>],
    uncovered: &mut [bool],
    uncovered_count: usize,
    memberships: &mut [u8],
    rejected_states: &mut HashSet<Vec<u8>>,
) -> bool {
    if uncovered_count == 0 {
        return true;
    }

    let mut state = memberships.to_vec();
    for edge_chunk in uncovered.chunks(8) {
        let mut byte = 0u8;
        for (bit, &is_uncovered) in edge_chunk.iter().enumerate() {
            if is_uncovered {
                byte |= 1 << bit;
            }
        }
        state.push(byte);
    }
    if !rejected_states.insert(state) {
        return false;
    }

    let feasible = |candidate: &KrauszClique, uncovered: &[bool], memberships: &[u8]| {
        candidate.edges.iter().all(|&edge| uncovered[edge])
            && candidate
                .vertices
                .iter()
                .all(|&vertex| memberships[vertex] < 2)
    };

    let mut best_candidates = Vec::new();
    for edge_index in 0..uncovered.len() {
        if !uncovered[edge_index] {
            continue;
        }
        let edge_candidates: Vec<_> = candidates_by_edge[edge_index]
            .iter()
            .copied()
            .filter(|&index| feasible(&candidates[index], uncovered, memberships))
            .collect();
        if edge_candidates.is_empty() {
            return false;
        }
        if best_candidates.is_empty() || edge_candidates.len() < best_candidates.len() {
            best_candidates = edge_candidates;
        }
    }

    for candidate_index in best_candidates {
        let candidate = &candidates[candidate_index];
        for &edge in &candidate.edges {
            uncovered[edge] = false;
        }
        for &vertex in &candidate.vertices {
            memberships[vertex] += 1;
        }

        if krausz_partition_exists(
            candidates,
            candidates_by_edge,
            uncovered,
            uncovered_count - candidate.edges.len(),
            memberships,
            rejected_states,
        ) {
            return true;
        }

        for &vertex in &candidate.vertices {
            memberships[vertex] -= 1;
        }
        for &edge in &candidate.edges {
            uncovered[edge] = true;
        }
    }
    false
}

impl std::fmt::Display for Graph {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Graph(n={}, |E|={})", self.n, self.edges.len())
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -- Generator tests --

    #[test]
    fn test_complete() {
        let k4 = Graph::complete(4);
        assert_eq!(k4.num_vertices(), 4);
        assert_eq!(k4.num_edges(), 6);

        let k5 = Graph::complete(5);
        assert_eq!(k5.num_edges(), 10);
    }

    #[test]
    fn test_path() {
        let p5 = Graph::path(5);
        assert_eq!(p5.num_vertices(), 5);
        assert_eq!(p5.num_edges(), 4);
        assert!(p5.has_edge(0, 1));
        assert!(!p5.has_edge(0, 2));
    }

    #[test]
    fn test_prufer_sequence_path() {
        let tree = Graph::from_prufer_sequence(&[1, 2, 3]).unwrap();
        assert_eq!(tree, Graph::path(5));
    }

    #[test]
    fn test_prufer_sequence_star() {
        let tree = Graph::from_prufer_sequence(&[0, 0]).unwrap();
        assert_eq!(tree, Graph::new(4, &[(0, 1), (0, 2), (0, 3)]));
    }

    #[test]
    fn test_prufer_sequence_rejects_invalid_input() {
        assert!(Graph::tree_from_prufer_sequence(4, &[0]).is_err());
        assert!(Graph::tree_from_prufer_sequence(4, &[0, 4]).is_err());
        assert!(Graph::tree_from_prufer_sequence(0, &[]).is_err());
    }

    #[test]
    fn test_all_labeled_trees_counts() {
        let expected = [(0, 0), (1, 1), (2, 1), (3, 3), (4, 16), (5, 125)];
        for (n, count) in expected {
            let trees = Graph::all_labeled_trees(n);
            assert_eq!(trees.len(), count);
            assert!(trees.iter().all(|tree| tree.num_vertices() == n));
            assert!(trees.iter().all(|tree| n == 0 || tree.num_edges() == n - 1));
            assert!(trees.iter().all(|tree| n == 0 || tree.is_connected()));
        }
    }

    #[test]
    fn test_cycle() {
        let c5 = Graph::cycle(5);
        assert_eq!(c5.num_edges(), 5);
        assert!(c5.has_edge(0, 4));
    }

    #[test]
    fn test_grid_and_ladder() {
        let grid = Graph::grid(3, 4);
        assert_eq!(grid.num_vertices(), 12);
        assert_eq!(grid.num_edges(), 17);
        assert!(grid.has_edge(0, 1));
        assert!(grid.has_edge(0, 4));
        assert!(!grid.has_edge(0, 5));

        let ladder = Graph::ladder(4);
        assert_eq!(ladder.num_vertices(), 8);
        assert_eq!(ladder.num_edges(), 10);
        assert!(ladder.has_edge(0, 1));
        assert!(ladder.has_edge(0, 2));
    }

    #[test]
    fn test_fan_wheel_and_friendship() {
        let fan = Graph::fan(4);
        assert_eq!(fan.num_vertices(), 5);
        assert_eq!(fan.num_edges(), 7);
        assert!(fan.has_edge(0, 4));
        assert!(fan.has_edge(2, 3));

        let wheel = Graph::wheel(5);
        assert_eq!(wheel.num_vertices(), 6);
        assert_eq!(wheel.num_edges(), 10);
        assert!(wheel.has_edge(0, 5));
        assert!(wheel.has_edge(1, 5));

        let friendship = Graph::friendship(3);
        assert_eq!(friendship.num_vertices(), 7);
        assert_eq!(friendship.num_edges(), 9);
        assert!(friendship.has_edge(1, 2));
        assert!(friendship.has_edge(5, 6));
        assert!(!friendship.has_edge(1, 3));
    }

    #[test]
    fn test_tree_family_generators() {
        let star = Graph::star(4);
        assert_eq!(star.num_vertices(), 5);
        assert_eq!(star.num_edges(), 4);
        assert!(star.has_edge(0, 4));

        let spider = Graph::spider(&[1, 2, 3]);
        assert_eq!(spider.num_vertices(), 7);
        assert_eq!(spider.num_edges(), 6);
        assert!(spider.has_edge(0, 1));
        assert!(spider.has_edge(2, 3));
        assert!(spider.has_edge(5, 6));

        let broom = Graph::broom(4, 3);
        assert_eq!(broom.num_vertices(), 7);
        assert_eq!(broom.num_edges(), 6);
        assert!(broom.has_edge(0, 1));
        assert!(broom.has_edge(0, 6));

        let double_star = Graph::double_star(2, 3);
        assert_eq!(double_star.num_vertices(), 7);
        assert_eq!(double_star.num_edges(), 6);
        assert!(double_star.has_edge(0, 1));
        assert!(double_star.has_edge(0, 3));
        assert!(double_star.has_edge(1, 6));

        let caterpillar = Graph::caterpillar(&[1, 0, 2]);
        assert_eq!(caterpillar.num_vertices(), 6);
        assert_eq!(caterpillar.num_edges(), 5);
        assert!(caterpillar.has_edge(0, 1));
        assert!(caterpillar.has_edge(2, 5));

        let binary_tree = Graph::complete_binary_tree(3);
        assert_eq!(binary_tree.num_vertices(), 15);
        assert_eq!(binary_tree.num_edges(), 14);
        assert!(binary_tree.has_edge(0, 1));
        assert!(binary_tree.has_edge(2, 6));
    }

    #[test]
    fn test_bipartite() {
        let k23 = Graph::complete_bipartite(2, 3);
        assert_eq!(k23.num_vertices(), 5);
        assert_eq!(k23.num_edges(), 6);
        assert!(k23.is_bipartite());
    }

    #[test]
    fn test_multipartite() {
        let k222 = Graph::complete_multipartite(&[2, 2, 2]);
        assert_eq!(k222.num_vertices(), 6);
        assert_eq!(k222.num_edges(), 12); // 3 * C(2,1)^2 = wait, 3*4 = 12
        assert!(!k222.has_edge(0, 1)); // same part
        assert!(k222.has_edge(0, 2)); // different parts
    }

    #[test]
    fn test_ferrers_board() {
        let lam = Partition::new(vec![3, 2]);
        let g = Graph::ferrers_board(&lam);
        assert_eq!(g.num_vertices(), 5); // 2 rows + 3 cols
        assert_eq!(g.num_edges(), 5); // boxes: (0,0),(0,1),(0,2),(1,0),(1,1)
    }

    #[test]
    fn test_unit_interval() {
        // area [0,1,2] gives edges: 1-2 (gap 1 ≤ a[2]=2), 0-2 (gap 2 ≤ a[2]=2), 0-1 (gap 1 ≤ a[1]=1)
        let g = Graph::unit_interval(&[0, 1, 2]);
        assert_eq!(g.num_vertices(), 3);
        assert_eq!(g.num_edges(), 3); // complete graph K3
    }

    #[test]
    fn test_circular_unit_interval_directed_edges() {
        assert!(Graph::is_circular_unit_interval_area_sequence(&[1, 1, 1]));
        assert!(!Graph::is_circular_unit_interval_area_sequence(&[0, 2, 0]));

        let directed = Graph::circular_unit_interval_directed_edges(&[1, 1, 1]).unwrap();
        assert_eq!(directed, vec![(0, 1), (1, 2), (2, 0)]);

        let underlying = Graph::circular_unit_interval(&[1, 1, 1]).unwrap();
        assert_eq!(underlying, Graph::cycle(3));
    }

    // -- Predicate tests --

    #[test]
    fn test_connected() {
        assert!(Graph::complete(4).is_connected());
        assert!(Graph::path(5).is_connected());
        assert!(Graph::cycle(5).is_connected());

        // Disconnected: two components
        let g = Graph::new(4, &[(0, 1), (2, 3)]);
        assert!(!g.is_connected());
    }

    #[test]
    fn test_bipartite_predicate() {
        assert!(Graph::path(5).is_bipartite());
        assert!(!Graph::complete(3).is_bipartite()); // K3 has odd cycle
        assert!(Graph::cycle(4).is_bipartite()); // even cycle
        assert!(!Graph::cycle(5).is_bipartite()); // odd cycle
    }

    #[test]
    fn test_claw_free() {
        // K4 is claw-free (every triple of neighbors has edges)
        assert!(Graph::complete(4).is_claw_free());
        // Path is claw-free (max degree 2)
        assert!(Graph::path(5).is_claw_free());
        // Star K_{1,3} is NOT claw-free
        let star = Graph::new(4, &[(0, 1), (0, 2), (0, 3)]);
        assert!(!star.is_claw_free());
    }

    // -- Matching tests (verified against Mathematica) --

    #[test]
    fn test_matching_polynomial_k3() {
        // Mathematica: {1, 3}
        assert_eq!(Graph::complete(3).matching_polynomial(), vec![1, 3]);
    }

    #[test]
    fn test_matching_polynomial_k4() {
        // Mathematica: {1, 6, 3}
        assert_eq!(Graph::complete(4).matching_polynomial(), vec![1, 6, 3]);
    }

    #[test]
    fn test_matching_polynomial_k5() {
        // Mathematica: {1, 10, 15}
        assert_eq!(Graph::complete(5).matching_polynomial(), vec![1, 10, 15]);
    }

    #[test]
    fn test_matching_polynomial_c5() {
        // Mathematica: {1, 5, 5}
        assert_eq!(Graph::cycle(5).matching_polynomial(), vec![1, 5, 5]);
    }

    #[test]
    fn test_matching_polynomial_p5() {
        // Mathematica: {1, 4, 3}
        assert_eq!(Graph::path(5).matching_polynomial(), vec![1, 4, 3]);
    }

    // -- Independence polynomial tests --

    #[test]
    fn test_independence_polynomial_p4() {
        // Mathematica: {1, 4, 3}
        assert_eq!(Graph::path(4).independence_polynomial(), vec![1, 4, 3]);
    }

    #[test]
    fn test_independence_polynomial_k4() {
        // K4: only empty set and singletons are independent
        assert_eq!(Graph::complete(4).independence_polynomial(), vec![1, 4]);
    }

    // -- Proper coloring tests --

    #[test]
    fn test_chromatic_polynomial_k3() {
        // chi_{K3}(t) = t(t-1)(t-2)
        // chi(1) = 0, chi(2) = 0, chi(3) = 6
        assert_eq!(Graph::complete(3).chromatic_polynomial_eval(1), 0);
        assert_eq!(Graph::complete(3).chromatic_polynomial_eval(2), 0);
        assert_eq!(Graph::complete(3).chromatic_polynomial_eval(3), 6);
    }

    #[test]
    fn test_chromatic_polynomial_c4() {
        // Mathematica: t^4 - 4t^3 + 6t^2 - 3t
        // chi(2) = 16 - 32 + 24 - 6 = 2
        // chi(3) = 81 - 108 + 54 - 9 = 18
        assert_eq!(Graph::cycle(4).chromatic_polynomial_eval(2), 2);
        assert_eq!(Graph::cycle(4).chromatic_polynomial_eval(3), 18);
    }

    // -- Acyclic orientation tests (verified against Mathematica) --

    #[test]
    fn test_acyclic_orientations_p3() {
        // Mathematica: 4
        assert_eq!(Graph::path(3).num_acyclic_orientations(), 4);
    }

    #[test]
    fn test_acyclic_orientations_c4() {
        // Mathematica: 14
        assert_eq!(Graph::cycle(4).num_acyclic_orientations(), 14);
    }

    #[test]
    fn test_acyclic_orientations_k4() {
        // Mathematica: 24 (= 4!)
        assert_eq!(Graph::complete(4).num_acyclic_orientations(), 24);
    }

    #[test]
    fn test_acyclic_orientations_with_frozen_edges_path() {
        let g = Graph::path(3);
        let orientations = g.acyclic_orientations_with_frozen_edges(&[(1, 0)]);
        assert_eq!(orientations.len(), 2);
        assert!(orientations.contains(&vec![(1, 0), (1, 2)]));
        assert!(orientations.contains(&vec![(1, 0), (2, 1)]));
    }

    #[test]
    fn test_acyclic_orientations_with_frozen_cycle_has_no_extensions() {
        let g = Graph::cycle(3);
        assert_eq!(
            g.num_acyclic_orientations_with_frozen_edges(&[(0, 1), (1, 2), (2, 0)]),
            0
        );
        assert_eq!(
            g.acyclic_sink_polynomial_with_frozen_edges(&[(0, 1), (1, 2), (2, 0)]),
            vec![0]
        );
    }

    #[test]
    fn test_acyclic_sink_polynomial_with_frozen_edges_path() {
        let g = Graph::path(3);
        assert_eq!(
            g.acyclic_sink_polynomial_with_frozen_edges(&[(1, 0)]),
            vec![0, 1, 1]
        );
    }

    #[test]
    fn test_acyclic_sink_polynomial_empty_graph() {
        assert_eq!(Graph::empty(0).acyclic_sink_polynomial(), vec![1]);
    }

    #[test]
    fn test_for_each_orientation_matches_all_orientations() {
        let g = Graph::path(3);
        let mut count = 0;
        g.for_each_orientation(|_| count += 1);
        assert_eq!(count, g.all_orientations().len());
    }

    #[test]
    fn test_for_each_acyclic_orientation_with_frozen_edges_matches_count() {
        let g = Graph::path(3);
        let mut count = 0;
        g.for_each_acyclic_orientation_with_frozen_edges(&[(1, 0)], |_| count += 1);
        assert_eq!(
            count,
            g.num_acyclic_orientations_with_frozen_edges(&[(1, 0)])
        );
    }

    // -- Graph operation tests --

    #[test]
    fn test_complement() {
        let p3 = Graph::path(3); // 0-1, 1-2
        let comp = p3.complement(); // 0-2
        assert_eq!(comp.num_edges(), 1);
        assert!(comp.has_edge(0, 2));
    }

    #[test]
    fn test_delete_vertex() {
        let k4 = Graph::complete(4);
        let g = k4.delete_vertex(0);
        assert_eq!(g.num_vertices(), 3);
        assert_eq!(g.num_edges(), 3); // K3
    }

    #[test]
    fn test_induced_subgraph_treats_vertices_as_a_set() {
        let graph = Graph::path(4);
        let induced = graph.induced_subgraph(&[2, 1, 2]);
        assert_eq!(induced.num_vertices(), 2);
        assert_eq!(induced.num_edges(), 1);
    }

    #[test]
    #[should_panic(expected = "outside the graph")]
    fn test_induced_subgraph_rejects_invalid_vertex() {
        let _ = Graph::path(3).induced_subgraph(&[0, 3]);
    }

    #[test]
    fn test_contract_edge() {
        let p3 = Graph::path(3); // 0-1-2
        let g = p3.contract_edge(0, 1); // merge 1 into 0, result: 0-1 (was 0-2)
        assert_eq!(g.num_vertices(), 2);
        assert_eq!(g.num_edges(), 1);
    }

    #[test]
    fn test_identify_nonadjacent_vertices() {
        let g = Graph::path(3); // 0-1-2
        let identified = g.identify_vertices(0, 2);
        assert_eq!(identified.num_vertices(), 2);
        assert_eq!(identified.num_edges(), 1);
        assert!(identified.has_edge(0, 1));
    }

    #[test]
    #[should_panic(expected = "cannot contract a non-edge")]
    fn test_contract_edge_rejects_non_edge() {
        let g = Graph::path(3);
        let _ = g.contract_edge(0, 2);
    }

    #[test]
    fn test_petersen() {
        let g = Graph::petersen();
        assert_eq!(g.num_vertices(), 10);
        assert_eq!(g.num_edges(), 15);
        // Petersen graph is 3-regular
        for v in 0..10 {
            assert_eq!(g.degree(v), 3);
        }
        // Petersen is NOT claw-free (it contains K_{1,3})
        assert!(!g.is_claw_free());
    }

    // -- Edge cases --

    #[test]
    fn test_empty_graph() {
        let g = Graph::empty(3);
        assert_eq!(g.num_edges(), 0);
        assert_eq!(g.matching_polynomial(), vec![1]);
        assert_eq!(g.independence_polynomial(), vec![1, 3, 3, 1]); // all subsets
    }

    #[test]
    fn test_single_vertex() {
        let g = Graph::empty(1);
        assert!(g.is_connected());
        assert!(g.is_bipartite());
        assert!(g.is_claw_free());
    }

    // -- graph6 parser tests --

    #[test]
    fn test_graph6_k1() {
        // K1 = single vertex, graph6 = "@"
        let g = Graph::from_graph6("@").unwrap();
        assert_eq!(g.num_vertices(), 1);
        assert_eq!(g.num_edges(), 0);
    }

    #[test]
    fn test_graph6_k2() {
        // K2, graph6 = "A_"
        let g = Graph::from_graph6("A_").unwrap();
        assert_eq!(g.num_vertices(), 2);
        assert_eq!(g.num_edges(), 1);
    }

    #[test]
    fn test_graph6_k3() {
        // K3, graph6 = "Bw"
        let g = Graph::from_graph6("Bw").unwrap();
        assert_eq!(g.num_vertices(), 3);
        assert_eq!(g.num_edges(), 3);
    }

    #[test]
    fn test_graph6_rejects_invalid_or_truncated_input() {
        assert!(Graph::from_graph6("!").is_err());
        assert!(Graph::from_graph6("A").is_err());
        assert!(Graph::from_graph6("~??").is_err());
    }

    #[test]
    #[should_panic(expected = "requires a connected tree")]
    fn test_tree_sink_polynomial_rejects_disconnected_graph() {
        let graph = Graph::new(4, &[(0, 1), (1, 2), (0, 2)]);
        let _ = graph.sink_polynomial_tree_bigint();
    }

    #[test]
    fn test_graph6_c5() {
        // C5, graph6 = "Dhc"
        let g = Graph::from_graph6("Dhc").unwrap();
        assert_eq!(g.num_vertices(), 5);
        assert_eq!(g.num_edges(), 5);
        assert!(g.has_edge(0, 4)); // cycle edge
    }

    #[test]
    fn test_graph6_petersen() {
        // Petersen graph, graph6 = "IsP@DKAO?"  (well-known encoding)
        // Let's verify by round-tripping: parse graph5c.g6 and check counts
        let path = std::path::Path::new("/home/paxinum/Dropbox/mathematica-packages/graph5c.g6");
        if path.exists() {
            let graphs = Graph::all_from_graph6_file(path).unwrap();
            // Number of connected graphs on 5 vertices = 21
            assert_eq!(graphs.len(), 21);
            // All should have 5 vertices
            assert!(graphs.iter().all(|g| g.num_vertices() == 5));
            // All should be connected
            assert!(graphs.iter().all(|g| g.is_connected()));
        }
    }

    #[test]
    fn test_graph6_file_counts() {
        // Known counts of connected graphs on n vertices (OEIS A001349)
        let expected = [(3, 2), (4, 6), (5, 21), (6, 112), (7, 853)];
        for (n, count) in expected {
            let path = std::path::PathBuf::from(format!(
                "/home/paxinum/Dropbox/mathematica-packages/graph{}c.g6",
                n
            ));
            if path.exists() {
                let graphs = Graph::all_from_graph6_file(&path).unwrap();
                assert_eq!(
                    graphs.len(),
                    count,
                    "wrong count for connected graphs on {} vertices",
                    n
                );
            }
        }
    }

    #[test]
    fn test_graph6_trees() {
        // Known counts of trees on n vertices (OEIS A000055)
        let expected = [(5, 3), (6, 6), (7, 11), (8, 23), (9, 47), (10, 106)];
        for (n, count) in expected {
            let path = std::path::PathBuf::from(format!(
                "/home/paxinum/Dropbox/mathematica-packages/trees{}.g6",
                n
            ));
            if path.exists() {
                let trees = Graph::all_from_graph6_file(&path).unwrap();
                assert_eq!(
                    trees.len(),
                    count,
                    "wrong count for trees on {} vertices",
                    n
                );
                // All trees are connected and bipartite
                assert!(trees.iter().all(|g| g.is_connected()));
                assert!(trees.iter().all(|g| g.is_bipartite()));
                // All trees have n-1 edges
                assert!(trees.iter().all(|g| g.num_edges() == n - 1));
            }
        }
    }

    // -- Perfect matchings --

    #[test]
    fn test_perfect_matchings_k4() {
        // K4 has 3 perfect matchings
        assert_eq!(Graph::complete(4).perfect_matchings().len(), 3);
    }

    #[test]
    fn test_perfect_matchings_odd() {
        assert_eq!(Graph::complete(3).perfect_matchings().len(), 0);
    }

    // -- Non-crossing matchings --

    #[test]
    fn test_non_crossing_k4() {
        // K4: total matchings = 1 + 6 + 3 = 10. Non-crossing: all singles are
        // non-crossing (7), pairs: (0-1,2-3) ok, (0-2,1-3) crosses, (0-3,1-2) ok.
        // So non-crossing matchings with 2 edges: 2. Total: 1 + 6 + 2 = 9.
        // Wait: empty + 6 single-edge + non-crossing pairs.
        // Pairs from K4: {01,23}, {02,13}, {03,12}
        // {02,13}: 0<1<2 and edge (1,3): 1<3, edge (0,2): 0<2. Cross? 0<1<2<3: no cross.
        // Actually: (0,2) and (1,3): 0 < 1 < 2 < 3 → a=0,b=2,c=1,d=3 → 0<1<2<3 → a<c<b<d → crossing!
        let nc = Graph::complete(4).non_crossing_matchings();
        let nc2: Vec<_> = nc.iter().filter(|m| m.len() == 2).collect();
        assert_eq!(nc2.len(), 2); // {0-1,2-3} and {0-3,1-2}
    }

    // -- Non-nesting matchings --

    #[test]
    fn test_non_nesting_k4() {
        // K4 pairs: {01,23} no nesting, {02,13} no nesting (neither contains other),
        // {03,12}: 0<1<2<3, edge (0,3) contains (1,2) → nesting!
        let nn = Graph::complete(4).non_nesting_matchings();
        let nn2: Vec<_> = nn.iter().filter(|m| m.len() == 2).collect();
        assert_eq!(nn2.len(), 2); // {0-1,2-3} and {0-2,1-3}
    }

    // -- Triangles --

    #[test]
    fn test_triangles_k4() {
        // K4 has C(4,3) = 4 triangles
        assert_eq!(Graph::complete(4).triangles().len(), 4);
    }

    #[test]
    fn test_triangles_c5() {
        // C5 has no triangles
        assert_eq!(Graph::cycle(5).triangles().len(), 0);
    }

    #[test]
    fn test_triangles_k5() {
        assert_eq!(Graph::complete(5).triangles().len(), 10);
    }

    // -- All orientations --

    #[test]
    fn test_all_orientations_p3() {
        // P3 has 2 edges → 4 orientations
        assert_eq!(Graph::path(3).all_orientations().len(), 4);
    }

    #[test]
    fn test_all_orientations_k3() {
        // K3 has 3 edges → 8 orientations
        assert_eq!(Graph::complete(3).all_orientations().len(), 8);
    }

    // -- Stirling graph --

    #[test]
    fn test_stirling_graph() {
        // S(3): matching polynomial should relate to Stirling numbers
        let g = Graph::stirling(3);
        assert_eq!(g.num_vertices(), 6);
        // Matchings by size: S(3,3)=1, S(3,2)=3, S(3,1)=1
        // matching poly = [1, 3, 1]
        assert_eq!(g.matching_polynomial(), vec![1, 3, 1]);
    }

    // -- Line graph --

    #[test]
    fn test_line_graph_k3() {
        // L(K3) = K3 (3 edges, each pair shares a vertex)
        let lg = Graph::complete(3).line_graph();
        assert_eq!(lg.num_vertices(), 3);
        assert_eq!(lg.num_edges(), 3);
    }

    #[test]
    fn test_line_graph_matching_independence() {
        // I(L(G), x) = mu(G, x) for any graph G
        let g = Graph::path(5);
        let lg = g.line_graph();
        assert_eq!(g.matching_polynomial(), lg.independence_polynomial());
    }

    #[test]
    fn test_is_line_graph_basic_examples() {
        assert!(Graph::empty(0).is_line_graph());
        assert!(Graph::empty(5).is_line_graph());
        assert!(Graph::complete(7).is_line_graph());
        assert!(Graph::cycle(7).is_line_graph());

        let claw = Graph::new(4, &[(0, 1), (0, 2), (0, 3)]);
        assert!(!claw.is_line_graph());

        // This unit interval graph is claw-free, but not a line graph.  It
        // exercises the Krausz search rather than the quick claw rejection.
        let claw_free_non_line_graph = Graph::unit_interval(&[0, 1, 2, 3, 3]);
        assert!(claw_free_non_line_graph.is_claw_free());
        assert!(!claw_free_non_line_graph.is_line_graph());

        // The root contains a triangle and a pendant edge. Its line graph has
        // overlapping maximal cliques, so recognizing it needs two-vertex
        // Krausz cliques as well as the maximal ones.
        let triangle_with_leaf = Graph::new(4, &[(0, 1), (1, 2), (0, 2), (0, 3)]);
        assert!(triangle_with_leaf.line_graph().is_line_graph());
    }

    #[test]
    fn test_is_line_graph_for_deterministic_random_roots() {
        let mut state = 0x4d595df4d0f33173u64;
        for vertices in 1..=8 {
            for _ in 0..16 {
                let mut edges = Vec::new();
                for u in 0..vertices {
                    for v in (u + 1)..vertices {
                        state = state
                            .wrapping_mul(6364136223846793005)
                            .wrapping_add(1442695040888963407);
                        if state >> 63 == 1 {
                            edges.push((u, v));
                        }
                    }
                }
                let root = Graph::new(vertices, &edges);
                assert!(
                    root.line_graph().is_line_graph(),
                    "failed for root graph with {vertices} vertices and edges {edges:?}"
                );
            }
        }
    }
}
