//! Exact recognition of polynomial rows as rows of a planar strip digraph.
//!
//! Given rows `P_1, ..., P_N` (with `P_0 = 1` implicit) this module searches for
//! a [`BrentiStripDigraph`] with weight matrix `W` such that
//! `v_0 = e_0`, `v_{n+1} = v_n W`, and `v_n` is the coefficient vector of `P_n`
//! for every supplied row.  The answer is one of
//!
//! - [`BrentiRecognition::Certified`]: a digraph, re-verified by
//!   [`build_brenti_sequence_certificate`];
//! - [`BrentiRecognition::ImpossibleWithinBounds`]: an exact proof that no such
//!   digraph exists within the stated bounds;
//! - [`BrentiRecognition::Inconclusive`]: a search budget was exhausted.
//!
//! # What a certificate means
//!
//! A certificate shows that the supplied rows are rows of the path matrix of
//! a planar network, so their coefficient matrix (rows `P_0, ..., P_N`) is
//! totally nonnegative; see [`crate::brenti_sequence`].  It does **not** show
//! that the rows are real-rooted: `[1 + t, 1 + t + t^2]` is recognized.  It
//! also covers only the supplied prefix; whether the digraphs found for longer
//! prefixes stabilize is a separate question.
//!
//! # Method
//!
//! Heights above the largest row degree `D` are never reached by a matching
//! digraph (a positive-weight edge into them from a reached height would give
//! a nonzero coefficient beyond `D`), so the strip has heights `0..=D`.  Only
//! heights reached by some `v_n` with `n < N` (*active* heights) constrain
//! the equations; edges out of the other heights are set to zero.
//!
//! The planar supports are described by target intervals: active height `y`
//! has targets in `[a_y, b_y]`, with `b_y <= a_{y'}` whenever `y < y'`, and
//! with at most [`BrentiRecognitionOptions::max_target_span`] consecutive
//! targets.  Every planar support whose target sets have span within the bound
//! lies inside a *maximal* such interval support, and a feasible weighting of
//! a smaller support is a feasible weighting of the larger one with some
//! weights zero.  Conversely, dropping zero-weight edges from a weighting of
//! an interval support leaves a planar digraph, because removing edges can
//! only lower the largest target of a height and raise the smallest one.  So
//! it suffices to search the maximal interval supports and to drop zero
//! weights from the result.
//!
//! The equations `v_n W = v_{n+1}` split by target column `j`: the unknowns
//! `W[y][j]` for the active heights `y` whose interval contains `j` (a
//! contiguous run of heights) must satisfy
//! `sum_y v_n[y] W[y][j] = v_{n+1}[j]` for all `n < N`, with `W >= 0`.  Each
//! column problem is decided exactly over `BigRational`: Gaussian elimination
//! first, then, if free parameters remain, an exact phase-one simplex method
//! with Bland's rule.  In integer mode, a unique rational solution must be
//! integral; otherwise a depth-first search over the finite integer box
//! (every active variable is bounded because all coefficients are
//! nonnegative), pruned by the rational relaxation, decides integer
//! feasibility.  Column outcomes are memoized, and supports are enumerated
//! height by height with a column checked as soon as its sources are final.
//!
//! A column that is infeasible even when every active height may reach it is
//! reported as an obstruction that does not depend on planarity.

use crate::brenti_sequence::{
    build_brenti_sequence_certificate, BrentiError, BrentiSequenceCertificate, BrentiStripDigraph,
};
use crate::linalg::solve_linear_system_with_rank;
use crate::tnn_network::BigRational;
use crate::Polynomial;
use num_bigint::BigInt;
use num_traits::{One, Zero};
use std::collections::HashMap;
use std::fmt;

type Q = BigRational;

/// Which edge weights the recognition search may use.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BrentiWeightDomain {
    /// Nonnegative rational weights.
    #[default]
    Rational,
    /// Nonnegative integer weights.
    Integer,
}

/// Bounds and budgets for [`recognize_brenti_strip_digraph`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BrentiRecognitionOptions {
    /// Digraph row index of the first supplied polynomial: `1` (default) means
    /// the input is `P_1, P_2, ...` with `P_0 = 1` implicit; `0` means the
    /// input starts with `P_0`, which must equal `1`.  Larger offsets would
    /// leave unknown intermediate rows (a nonlinear problem) and are rejected.
    pub first_row_index: usize,
    /// Maximal number of consecutive target heights per height (`None`: no
    /// bound).  This bounds the out-degree; note that a planar digraph with
    /// out-degree `d` may still need a larger span, e.g. `0 -> {0, 2}`.
    pub max_target_span: Option<usize>,
    /// Rational or integer weights.
    pub weights: BrentiWeightDomain,
    /// Budget on the number of partial supports visited.
    pub max_search_nodes: usize,
    /// Budget on the integer search for one column problem.
    pub max_integer_search_nodes: usize,
}

impl Default for BrentiRecognitionOptions {
    fn default() -> Self {
        Self {
            first_row_index: 1,
            max_target_span: None,
            weights: BrentiWeightDomain::Rational,
            max_search_nodes: 1_000_000,
            max_integer_search_nodes: 10_000,
        }
    }
}

/// Outcome of [`recognize_brenti_strip_digraph`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BrentiRecognition {
    /// A planar strip digraph reproduces every supplied row exactly.  The
    /// certificate covers the supplied prefix only, and certifies total
    /// nonnegativity of the coefficient matrix, not real-rootedness of rows.
    Certified {
        digraph: BrentiStripDigraph,
        certificate: BrentiSequenceCertificate,
    },
    /// No planar strip digraph within the bounds reproduces the rows.  This
    /// is exact.  `supports_checked` counts the partial interval supports
    /// visited (every maximal support extends one of them); `obstruction`
    /// explains the failure when a single reason is known.
    ImpossibleWithinBounds {
        supports_checked: usize,
        obstruction: String,
    },
    /// A budget was exhausted before the search was decided.
    Inconclusive { reason: String },
}

impl BrentiRecognition {
    /// Whether a certificate was found.
    pub fn is_certified(&self) -> bool {
        matches!(self, Self::Certified { .. })
    }
}

impl fmt::Display for BrentiRecognition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Certified { certificate, .. } => write!(
                f,
                "certified: rows {:?} are rows of a planar strip digraph \
                 (this prefix only; the coefficient matrix is TNN)",
                certificate.row_indices
            ),
            Self::ImpossibleWithinBounds {
                supports_checked,
                obstruction,
            } => write!(
                f,
                "impossible within bounds ({} partial supports checked): {}",
                supports_checked, obstruction
            ),
            Self::Inconclusive { reason } => write!(f, "inconclusive: {}", reason),
        }
    }
}

/// Invalid input for [`recognize_brenti_strip_digraph`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BrentiRecognitionError {
    EmptyInput,
    UnsupportedFirstRowIndex {
        found: usize,
    },
    FirstRowNotOne {
        found: Vec<BigInt>,
    },
    ZeroTargetSpan,
    /// The digraph found by the search failed re-verification (a bug).
    CertificateRejected(BrentiError),
}

impl fmt::Display for BrentiRecognitionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyInput => write!(f, "at least one row is required"),
            Self::UnsupportedFirstRowIndex { found } => {
                write!(f, "first row index must be 0 or 1, found {}", found)
            }
            Self::FirstRowNotOne { found } => write!(
                f,
                "with first row index 0 the first row must be 1, found coefficients {:?}",
                found
            ),
            Self::ZeroTargetSpan => write!(f, "max_target_span must be at least 1"),
            Self::CertificateRejected(error) => {
                write!(f, "internal error: certificate rejected: {}", error)
            }
        }
    }
}

impl std::error::Error for BrentiRecognitionError {}

/// Search for a planar strip digraph whose rows are the given polynomials.
///
/// `rows` are rows `first_row_index, first_row_index + 1, ...` of the digraph
/// (see [`BrentiRecognitionOptions`]); coefficients are in ascending degree
/// order.  See the module documentation for the method and for what a
/// certificate means.
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::{recognize_brenti_strip_digraph, BrentiRecognition, BrentiRecognitionOptions};
/// use polytool::Polynomial;
///
/// let rows: Vec<Polynomial<BigInt>> = [vec![1, 1], vec![1, 2, 1], vec![1, 3, 3, 1]]
///     .iter()
///     .map(|r| Polynomial::from_i64_coeffs(r))
///     .collect();
/// let result = recognize_brenti_strip_digraph(&rows, &BrentiRecognitionOptions::default()).unwrap();
/// match result {
///     BrentiRecognition::Certified { digraph, .. } => {
///         // y -> y and y -> y + 1 with weight 1 below the top height.
///         assert_eq!(digraph.outgoing[0].len(), 2);
///     }
///     other => panic!("expected a certificate, got {}", other),
/// }
/// ```
pub fn recognize_brenti_strip_digraph(
    rows: &[Polynomial<BigInt>],
    options: &BrentiRecognitionOptions,
) -> Result<BrentiRecognition, BrentiRecognitionError> {
    if rows.is_empty() {
        return Err(BrentiRecognitionError::EmptyInput);
    }
    if options.first_row_index > 1 {
        return Err(BrentiRecognitionError::UnsupportedFirstRowIndex {
            found: options.first_row_index,
        });
    }
    if options.max_target_span == Some(0) {
        return Err(BrentiRecognitionError::ZeroTargetSpan);
    }

    let mut vectors: Vec<Vec<BigInt>> = Vec::with_capacity(rows.len() + 1);
    if options.first_row_index == 1 {
        vectors.push(vec![BigInt::one()]);
    } else if rows[0].coeffs() != [BigInt::one()] {
        return Err(BrentiRecognitionError::FirstRowNotOne {
            found: rows[0].coeffs().to_vec(),
        });
    }
    vectors.extend(rows.iter().map(|p| p.coeffs().to_vec()));

    let max_height = vectors
        .iter()
        .map(|v| v.len().saturating_sub(1))
        .max()
        .unwrap_or(0);
    let size = max_height + 1;
    for v in &mut vectors {
        v.resize(size, BigInt::zero());
    }

    for (row, v) in vectors.iter().enumerate() {
        if let Some((k, c)) = v.iter().enumerate().find(|(_, c)| **c < BigInt::zero()) {
            return Ok(BrentiRecognition::ImpossibleWithinBounds {
                supports_checked: 0,
                obstruction: format!(
                    "row {} has negative coefficient {} at t^{}; path weights are nonnegative",
                    row, c, k
                ),
            });
        }
    }

    let row_indices: Vec<usize> = (0..rows.len())
        .map(|i| i + options.first_row_index)
        .collect();
    let mut search = Search::new(&vectors, options);
    let outcome = search.run();
    match outcome {
        SearchOutcome::Found(weights) => {
            let digraph = BrentiStripDigraph::from_weight_matrix(&weights);
            let certificate =
                build_brenti_sequence_certificate(digraph.clone(), &row_indices, rows)
                    .map_err(BrentiRecognitionError::CertificateRejected)?;
            Ok(BrentiRecognition::Certified {
                digraph,
                certificate,
            })
        }
        SearchOutcome::Impossible(obstruction) => Ok(BrentiRecognition::ImpossibleWithinBounds {
            supports_checked: search.nodes,
            obstruction,
        }),
        SearchOutcome::Inconclusive(reason) => Ok(BrentiRecognition::Inconclusive { reason }),
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum ColumnOutcome {
    Feasible(Vec<Q>),
    Infeasible,
    Inconclusive,
}

enum SearchOutcome {
    Found(Vec<Vec<Q>>),
    Impossible(String),
    Inconclusive(String),
}

struct Search<'a> {
    vectors: &'a [Vec<BigInt>],
    size: usize,
    num_equations: usize,
    active: Vec<usize>,
    span: usize,
    integer: bool,
    max_nodes: usize,
    max_integer_nodes: usize,
    cache: HashMap<(usize, usize, usize), ColumnOutcome>,
    nodes: usize,
    saw_inconclusive_column: bool,
    budget_exhausted: bool,
}

impl<'a> Search<'a> {
    fn new(vectors: &'a [Vec<BigInt>], options: &BrentiRecognitionOptions) -> Self {
        let size = vectors[0].len();
        let num_equations = vectors.len() - 1;
        let active = (0..size)
            .filter(|&y| vectors[..num_equations].iter().any(|v| !v[y].is_zero()))
            .collect();
        Self {
            vectors,
            size,
            num_equations,
            active,
            span: options.max_target_span.unwrap_or(size).min(size),
            integer: options.weights == BrentiWeightDomain::Integer,
            max_nodes: options.max_search_nodes,
            max_integer_nodes: options.max_integer_search_nodes,
            cache: HashMap::new(),
            nodes: 0,
            saw_inconclusive_column: false,
            budget_exhausted: false,
        }
    }

    fn run(&mut self) -> SearchOutcome {
        if self.active.is_empty() {
            // Only row 0 = e_0 is given (row 0 is never zero), so no equation
            // constrains any weight.
            return SearchOutcome::Found(vec![vec![Q::zero(); self.size]; self.size]);
        }

        let all = self.active.len();
        for j in 0..self.size {
            if self.column(j, 0, all) == ColumnOutcome::Infeasible {
                return SearchOutcome::Impossible(format!(
                    "target height {} has no admissible weights even if every active \
                     height may reach it (planarity is not needed for this obstruction)",
                    j
                ));
            }
        }

        let mut intervals = Vec::with_capacity(all);
        if let Some(weights) = self.extend(&mut intervals, 0) {
            return SearchOutcome::Found(weights);
        }
        if self.budget_exhausted {
            SearchOutcome::Inconclusive(format!(
                "support search budget of {} nodes exhausted",
                self.max_nodes
            ))
        } else if self.saw_inconclusive_column {
            SearchOutcome::Inconclusive(format!(
                "integer weight search budget of {} nodes exhausted for some column",
                self.max_integer_nodes
            ))
        } else {
            let domain = if self.integer { "integer" } else { "rational" };
            SearchOutcome::Impossible(if self.span >= self.size {
                format!("no planar support admits nonnegative {} weights", domain)
            } else {
                format!(
                    "no planar interval support with target span at most {} admits \
                     nonnegative {} weights",
                    self.span, domain
                )
            })
        }
    }

    /// Depth-first enumeration of maximal interval supports over the active
    /// heights; `checked` is the first column whose sources are not final.
    fn extend(
        &mut self,
        intervals: &mut Vec<(usize, usize)>,
        checked: usize,
    ) -> Option<Vec<Vec<Q>>> {
        let p = intervals.len();
        let last = p + 1 == self.active.len();
        let max_target = self.size - 1;
        let (lower, prev_full) = match intervals.last() {
            Some(&(a, b)) => (b, b - a + 1 == self.span),
            None => (0, true),
        };
        // Right-maximality of the previous interval: it touches this one
        // unless it already has full span.
        let a_range = if prev_full {
            lower..=max_target
        } else {
            lower..=lower
        };
        for a in a_range {
            let b_max = max_target.min(a + self.span - 1);
            for b in a..=b_max {
                let full = b - a + 1 == self.span;
                // Left-maximality: start at the lower bound unless full.
                if a != lower && !full {
                    continue;
                }
                // Right-maximality of the last interval.
                if last && b != max_target && !full {
                    continue;
                }
                if self.nodes >= self.max_nodes {
                    self.budget_exhausted = true;
                    return None;
                }
                self.nodes += 1;
                intervals.push((a, b));
                let final_upto = if last { self.size } else { a };
                if self.columns_feasible(intervals, checked, final_upto) {
                    if last {
                        let weights = self.assemble(intervals);
                        intervals.pop();
                        return Some(weights);
                    }
                    if let Some(found) = self.extend(intervals, final_upto.max(checked)) {
                        intervals.pop();
                        return Some(found);
                    }
                }
                intervals.pop();
                if self.budget_exhausted {
                    return None;
                }
            }
        }
        None
    }

    fn sources(intervals: &[(usize, usize)], j: usize) -> (usize, usize) {
        let lo = intervals
            .iter()
            .position(|&(_, b)| b >= j)
            .unwrap_or(intervals.len());
        let hi = intervals
            .iter()
            .rposition(|&(a, _)| a <= j)
            .map_or(0, |q| q + 1);
        if lo < hi {
            (lo, hi)
        } else {
            (0, 0)
        }
    }

    fn columns_feasible(&mut self, intervals: &[(usize, usize)], from: usize, to: usize) -> bool {
        for j in from..to {
            let (lo, hi) = Self::sources(intervals, j);
            match self.column(j, lo, hi) {
                ColumnOutcome::Feasible(_) => {}
                ColumnOutcome::Infeasible => return false,
                ColumnOutcome::Inconclusive => {
                    self.saw_inconclusive_column = true;
                    return false;
                }
            }
        }
        true
    }

    fn assemble(&mut self, intervals: &[(usize, usize)]) -> Vec<Vec<Q>> {
        let columns: Vec<(usize, ColumnOutcome)> = (0..self.size)
            .map(|j| {
                let (lo, hi) = Self::sources(intervals, j);
                (lo, self.column(j, lo, hi))
            })
            .collect();
        let mut weights = vec![vec![Q::zero(); self.size]; self.size];
        for (j, (lo, outcome)) in columns.into_iter().enumerate() {
            // Every column of a completed support was checked feasible.
            if let ColumnOutcome::Feasible(values) = outcome {
                for (&y, value) in self.active[lo..].iter().zip(values) {
                    weights[y][j] = value;
                }
            }
        }
        weights
    }

    /// Column problem for target `j` with sources `active[lo..hi]`.
    fn column(&mut self, j: usize, lo: usize, hi: usize) -> ColumnOutcome {
        if let Some(outcome) = self.cache.get(&(j, lo, hi)) {
            return outcome.clone();
        }
        let rhs: Vec<Q> = (0..self.num_equations)
            .map(|n| Q::from_integer(self.vectors[n + 1][j].clone()))
            .collect();
        let matrix: Vec<Vec<Q>> = (0..self.num_equations)
            .map(|n| {
                self.active[lo..hi]
                    .iter()
                    .map(|&y| Q::from_integer(self.vectors[n][y].clone()))
                    .collect()
            })
            .collect();
        let outcome =
            solve_nonnegative_column(&matrix, &rhs, hi - lo, self.integer, self.max_integer_nodes);
        self.cache.insert((j, lo, hi), outcome.clone());
        outcome
    }
}

/// Decide `A x = b`, `x >= 0` (and `x` integral when `integer`) exactly.
fn solve_nonnegative_column(
    a: &[Vec<Q>],
    b: &[Q],
    num_vars: usize,
    integer: bool,
    max_integer_nodes: usize,
) -> ColumnOutcome {
    if b.iter().all(Zero::is_zero) {
        return ColumnOutcome::Feasible(vec![Q::zero(); num_vars]);
    }
    if num_vars == 0 {
        return ColumnOutcome::Infeasible;
    }
    let Some(solution) = solve_linear_system_with_rank(a, b) else {
        return ColumnOutcome::Infeasible;
    };
    if solution.nullity == 0 {
        let x = solution.particular;
        let admissible = x
            .iter()
            .all(|v| *v >= Q::zero() && (!integer || v.is_integer()));
        return if admissible {
            ColumnOutcome::Feasible(x)
        } else {
            ColumnOutcome::Infeasible
        };
    }
    let Some(x) = nonnegative_solution(a, b, num_vars) else {
        return ColumnOutcome::Infeasible;
    };
    if !integer || x.iter().all(Q::is_integer) {
        return ColumnOutcome::Feasible(x);
    }
    let mut nodes = 0;
    let mut fixed = Vec::with_capacity(num_vars);
    match integer_solution(a, b, num_vars, &mut fixed, &mut nodes, max_integer_nodes) {
        Ok(Some(x)) => ColumnOutcome::Feasible(x),
        Ok(None) => ColumnOutcome::Infeasible,
        Err(()) => ColumnOutcome::Inconclusive,
    }
}

/// Depth-first integer search, fixing variables in order and pruning with the
/// rational relaxation.  Requires nonnegative `a` and `b`, which bound every
/// variable with a nonzero column; variables with a zero column are fixed to
/// zero.  `Err(())` means the node budget was exhausted.
fn integer_solution(
    a: &[Vec<Q>],
    b: &[Q],
    num_vars: usize,
    fixed: &mut Vec<Q>,
    nodes: &mut usize,
    max_nodes: usize,
) -> Result<Option<Vec<Q>>, ()> {
    if *nodes >= max_nodes {
        return Err(());
    }
    *nodes += 1;
    let k = fixed.len();
    let residual: Vec<Q> = a
        .iter()
        .zip(b)
        .map(|(row, rhs)| {
            let mut r = rhs.clone();
            for (coef, value) in row.iter().zip(fixed.iter()) {
                r -= coef * value;
            }
            r
        })
        .collect();
    if k == num_vars {
        return Ok(residual.iter().all(Zero::is_zero).then(|| fixed.clone()));
    }
    let rest: Vec<Vec<Q>> = a.iter().map(|row| row[k..].to_vec()).collect();
    let Some(relaxed) = nonnegative_solution(&rest, &residual, num_vars - k) else {
        return Ok(None);
    };
    if relaxed.iter().all(Q::is_integer) {
        let mut x = fixed.clone();
        x.extend(relaxed);
        return Ok(Some(x));
    }
    let bound = a
        .iter()
        .zip(&residual)
        .filter(|(row, _)| row[k] > Q::zero())
        .map(|(row, r)| (r / &row[k]).floor().to_integer())
        .min();
    let Some(bound) = bound else {
        fixed.push(Q::zero());
        let result = integer_solution(a, b, num_vars, fixed, nodes, max_nodes);
        fixed.pop();
        return result;
    };
    let mut value = BigInt::zero();
    while value <= bound {
        fixed.push(Q::from_integer(value.clone()));
        let result = integer_solution(a, b, num_vars, fixed, nodes, max_nodes);
        fixed.pop();
        match result {
            Ok(None) => {}
            other => return other,
        }
        value += 1;
    }
    Ok(None)
}

/// Exact phase-one simplex method with Bland's rule: a basic solution of
/// `A x = b`, `x >= 0`, or `None` if there is none.
///
/// Artificial variables start in the basis and are never re-entered once they
/// leave it; Bland's rule (smallest entering index, smallest leaving index
/// among ratio ties) prevents cycling.
fn nonnegative_solution(a: &[Vec<Q>], b: &[Q], num_vars: usize) -> Option<Vec<Q>> {
    let m = a.len();
    let mut tableau: Vec<Vec<Q>> = a
        .iter()
        .zip(b)
        .map(|(row, rhs)| {
            let mut r: Vec<Q> = row.iter().take(num_vars).cloned().collect();
            r.push(rhs.clone());
            if *rhs < Q::zero() {
                for entry in &mut r {
                    *entry = -entry.clone();
                }
            }
            r
        })
        .collect();
    // `None` marks the artificial variable of that row.
    let mut basis: Vec<Option<usize>> = vec![None; m];
    let index_of = |basic: Option<usize>, row: usize| basic.unwrap_or(num_vars + row);

    loop {
        let entering = (0..num_vars).find(|&j| {
            let mut reduced = Q::zero();
            for (row, basic) in tableau.iter().zip(&basis) {
                if basic.is_none() {
                    reduced -= &row[j];
                }
            }
            reduced < Q::zero()
        });
        let Some(j) = entering else {
            break;
        };
        let mut leaving: Option<(usize, Q)> = None;
        for (i, row) in tableau.iter().enumerate() {
            if row[j] <= Q::zero() {
                continue;
            }
            let ratio = &row[num_vars] / &row[j];
            let better = match &leaving {
                None => true,
                Some((best, best_ratio)) => {
                    ratio < *best_ratio
                        || (ratio == *best_ratio
                            && index_of(basis[i], i) < index_of(basis[*best], *best))
                }
            };
            if better {
                leaving = Some((i, ratio));
            }
        }
        // A negative reduced cost comes from a positive entry in some
        // artificial row, so a leaving row exists.
        let (r, _) = leaving.expect("phase-one objective is bounded below");
        let pivot = tableau[r][j].clone();
        for entry in &mut tableau[r] {
            *entry /= &pivot;
        }
        let pivot_row = tableau[r].clone();
        for (i, row) in tableau.iter_mut().enumerate() {
            if i == r || row[j].is_zero() {
                continue;
            }
            let factor = row[j].clone();
            for (entry, p) in row.iter_mut().zip(&pivot_row) {
                *entry -= &factor * p;
            }
        }
        basis[r] = Some(j);
    }

    if tableau
        .iter()
        .zip(&basis)
        .any(|(row, basic)| basic.is_none() && !row[num_vars].is_zero())
    {
        return None;
    }
    let mut x = vec![Q::zero(); num_vars];
    for (row, basic) in tableau.iter().zip(&basis) {
        if let Some(j) = basic {
            x[*j] = row[num_vars].clone();
        }
    }
    debug_assert!(a.iter().zip(b).all(|(row, rhs)| {
        row.iter()
            .zip(&x)
            .fold(Q::zero(), |acc, (c, v)| acc + c * v)
            == *rhs
    }));
    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(v: i64) -> Q {
        Q::from_integer(BigInt::from(v))
    }

    fn qm(rows: &[&[i64]]) -> Vec<Vec<Q>> {
        rows.iter()
            .map(|r| r.iter().map(|&v| q(v)).collect())
            .collect()
    }

    fn qv(v: &[i64]) -> Vec<Q> {
        v.iter().map(|&x| q(x)).collect()
    }

    fn satisfies(a: &[Vec<Q>], b: &[Q], x: &[Q]) -> bool {
        x.iter().all(|v| *v >= Q::zero())
            && a.iter().zip(b).all(|(row, rhs)| {
                row.iter().zip(x).fold(Q::zero(), |acc, (c, v)| acc + c * v) == *rhs
            })
    }

    #[test]
    fn simplex_feasible_underdetermined() {
        let a = qm(&[&[1, 1, 1], &[1, 2, 3]]);
        let b = qv(&[1, 2]);
        let x = nonnegative_solution(&a, &b, 3).unwrap();
        assert!(satisfies(&a, &b, &x));
    }

    #[test]
    fn simplex_infeasible_sign() {
        // x + y = 1, x - y = 3 forces y = -1.
        let a = qm(&[&[1, 1], &[1, -1]]);
        assert_eq!(nonnegative_solution(&a, &qv(&[1, 3]), 2), None);
        // x + y = -1 has no nonnegative solution.
        assert_eq!(nonnegative_solution(&qm(&[&[1, 1]]), &qv(&[-1]), 2), None);
    }

    #[test]
    fn simplex_infeasible_with_free_parameters() {
        // x - y + z = -1 and x + y + z = 1 force y = 1 and x + z = 0, so adding
        // x + z = 1 is infeasible while x + z = 0 is feasible.
        let a = qm(&[&[1, -1, 1], &[1, 1, 1], &[1, 0, 1]]);
        let b = qv(&[-1, 1, 1]);
        assert_eq!(nonnegative_solution(&a, &b, 3), None);
        let b = qv(&[-1, 1, 0]);
        let x = nonnegative_solution(&a, &b, 3).unwrap();
        assert!(satisfies(&a, &b, &x));
    }

    #[test]
    fn simplex_redundant_and_degenerate_rows() {
        let a = qm(&[&[1, 1, 0], &[2, 2, 0], &[0, 0, 0], &[0, 1, 1]]);
        let b = qv(&[1, 2, 0, 0]);
        let x = nonnegative_solution(&a, &b, 3).unwrap();
        assert!(satisfies(&a, &b, &x));
        assert_eq!(x, qv(&[1, 0, 0]));
    }

    /// Independent feasibility oracle: by the fundamental theorem of linear
    /// programming, `A x = b, x >= 0` is feasible iff some set of columns has
    /// a unique nonnegative solution.
    fn brute_force_feasible(a: &[Vec<Q>], b: &[Q], k: usize) -> bool {
        (0u32..(1 << k)).any(|mask| {
            let cols: Vec<usize> = (0..k).filter(|&j| mask & (1 << j) != 0).collect();
            if cols.is_empty() {
                return b.iter().all(Zero::is_zero);
            }
            let sub: Vec<Vec<Q>> = a
                .iter()
                .map(|row| cols.iter().map(|&j| row[j].clone()).collect())
                .collect();
            match solve_linear_system_with_rank(&sub, b) {
                Some(sol) => sol.nullity == 0 && sol.particular.iter().all(|v| *v >= Q::zero()),
                None => false,
            }
        })
    }

    #[test]
    fn simplex_matches_brute_force_on_random_systems() {
        let mut state: u64 = 0x2545_f491_4f6c_dd1d;
        let mut next = |modulus: u64| -> i64 {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((state >> 33) % modulus) as i64
        };
        let (mut feasible, mut infeasible) = (0, 0);
        for _ in 0..400 {
            let m = 1 + next(3) as usize;
            let k = 1 + next(4) as usize;
            let a: Vec<Vec<Q>> = (0..m)
                .map(|_| (0..k).map(|_| q(next(7) - 3)).collect())
                .collect();
            // A system feasible by construction.
            let x0: Vec<Q> = (0..k).map(|_| q(next(4))).collect();
            let b0: Vec<Q> = a
                .iter()
                .map(|row| {
                    row.iter()
                        .zip(&x0)
                        .fold(Q::zero(), |acc, (c, v)| acc + c * v)
                })
                .collect();
            let x = nonnegative_solution(&a, &b0, k).expect("feasible by construction");
            assert!(satisfies(&a, &b0, &x));
            // A random right-hand side, compared with the oracle.
            let b: Vec<Q> = (0..m).map(|_| q(next(9) - 4)).collect();
            let expected = brute_force_feasible(&a, &b, k);
            match nonnegative_solution(&a, &b, k) {
                Some(x) => {
                    assert!(expected, "simplex found a point the oracle missed");
                    assert!(satisfies(&a, &b, &x));
                    feasible += 1;
                }
                None => {
                    assert!(
                        !expected,
                        "simplex missed a feasible system {:?} {:?}",
                        a, b
                    );
                    infeasible += 1;
                }
            }
        }
        assert!(
            feasible > 20 && infeasible > 20,
            "{} {}",
            feasible,
            infeasible
        );
    }

    #[test]
    fn integer_column_search() {
        // 2x + 3y = 5: the LP vertex (5/2, 0) is fractional, (1, 1) is integral.
        let a = qm(&[&[2, 3]]);
        let b = qv(&[5]);
        match solve_nonnegative_column(&a, &b, 2, true, 1000) {
            ColumnOutcome::Feasible(x) => assert_eq!(x, qv(&[1, 1])),
            other => panic!("unexpected {:?}", other),
        }
        // 2x + 4y = 5 has rational but no integer solutions.
        let a = qm(&[&[2, 4]]);
        assert!(matches!(
            solve_nonnegative_column(&a, &b, 2, false, 1000),
            ColumnOutcome::Feasible(_)
        ));
        assert_eq!(
            solve_nonnegative_column(&a, &b, 2, true, 1000),
            ColumnOutcome::Infeasible
        );
        // Unique fractional solution.
        let a = qm(&[&[2]]);
        assert_eq!(
            solve_nonnegative_column(&a, &qv(&[1]), 1, true, 1000),
            ColumnOutcome::Infeasible
        );
        // A tiny budget is reported as inconclusive.
        let a = qm(&[&[2, 4, 6]]);
        assert_eq!(
            solve_nonnegative_column(&a, &qv(&[101]), 3, true, 2),
            ColumnOutcome::Inconclusive
        );
    }
}
