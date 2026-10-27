//! Height-shift-invariant planar networks: Pólya frequency certificates for
//! polynomial rows.
//!
//! This is the setting of Brenti's planar-network theorem on Pólya frequency
//! (PF) sequences (F. Brenti, *Combinatorics and total positivity*, J. Combin.
//! Theory Ser. A 71 (1995), 175--218; the recursive matrices of his
//! Theorem 4.3 are the motivating example).  It complements
//! [`crate::brenti_sequence`], whose transitions depend on the height and
//! which therefore certifies total nonnegativity of the *coefficient matrix*
//! only.
//!
//! # The networks
//!
//! A [`HeightInvariantNetwork`] lives on the lattice `{0, ..., N} x Z`.  For
//! every target column `n = 1..=N` it lists steps `(length, rise, weight)`;
//! each step stands for the edges `(n - length, y) -> (n, y + rise)` for
//! **every** height `y`, all with the same weight.  So the weights may depend
//! on the column, but not on the height: the network is invariant under the
//! shift `(x, y) -> (x, y + 1)`.  Steps have `length >= 1`, `rise >= 0` and
//! nonnegative weights, and edges are drawn as straight segments.
//!
//! Write `M_{n,k}` for the total weight of the paths from `(0, 0)` to
//! `(n, k)`; row `n` is the polynomial `P_n(t) = sum_k M_{n,k} t^k`, so
//! `P_0 = 1` and `P_n = sum_{steps into n} weight * t^rise * P_{n - length}`.
//!
//! # What a valid network certifies
//!
//! [`HeightInvariantNetwork::validate`] checks exactly that every
//! positive-weight edge is a primitive lattice segment
//! (`gcd(length, rise) = 1`, so no edge passes through a lattice point) and
//! that no two positive-weight edges meet except at a common endpoint.  Then
//! for every row `n` and every size `r`:
//!
//! 1. By height invariance, the path weight from `(0, i)` to `(n, j)` is
//!    `M_{n, j - i}`, so the path matrix from the sources `(0, 0..r)` to the
//!    sinks `(n, 0..r)` is the (upper triangular) Toeplitz matrix
//!    `[a_{j - i}]` of the coefficient sequence `a = (M_{n,k})_k`.
//! 2. Paths go strictly to the right, so a path from column `0` to column
//!    `n` is the graph of a continuous function on `[0, n]`.  Two
//!    vertex-disjoint paths have disjoint graphs (edges meet only at common
//!    endpoints, which are vertices, and no edge contains a vertex in its
//!    interior), so one lies strictly above the other: the sources and sinks
//!    are *fully compatible* in Brenti's sense, and only the identity pairing
//!    admits nonintersecting path families.
//! 3. By Lindström--Gessel--Viennot, every minor of that Toeplitz matrix is
//!    a sum of products of nonnegative weights, hence nonnegative.  As `r` is
//!    arbitrary, the Toeplitz matrix of row `n` is totally nonnegative: the
//!    coefficient sequence of `P_n` is a **Pólya frequency sequence**, so by
//!    Aissen--Schoenberg--Whitney (and Edrei) `P_n` is real-rooted with only
//!    nonpositive zeros (or is identically zero).
//!
//! This holds for every row of the network, but a certificate only covers the
//! rows it matches: a network for `P_0, ..., P_N` says nothing about
//! `P_{N+1}`.  It also says nothing about the coefficient matrix
//! `[M_{n,k}]_{n,k}` (that is a different matrix, with different sources).
//!
//! # What the verifier checks, and what a finite window means
//!
//! [`PfRowsCertificate::verify_against_bigint_polynomials`] (and
//! [`build_pf_rows_certificate`]) re-run the exact structural checks above,
//! recompute every matched row from the network, and rebuild, by dynamic
//! programming over the explicit lattice (without using height invariance),
//! the path matrix from the sources `(0, 0..=d)` to the sinks `(n, 0..=d)`
//! with `d = deg P_n`; it must equal the Toeplitz window
//! [`toeplitz_window`] of the row.  This window contains every coefficient,
//! so it ties the row to the path matrix; the PF conclusion comes from the
//! structural checks, for Toeplitz sections of every size.  A finite window
//! that merely happens to be TNN proves nothing by itself: the `3 x 3`
//! Toeplitz window of `1 + t + t^2` is TNN, while its `4 x 4` window has the
//! minor `-1`.
//!
//! # Recognition
//!
//! [`recognize_pf_network`] searches exactly for a network reproducing given
//! rows, within the bounds of [`PfNetworkOptions`]: step lengths at most
//! `max_step_length` and rises at most `max_rise` (all primitive step types
//! within these bounds are available, independently in every column).  The
//! weights into column `n` only affect `P_n`, so for a fixed set of steps
//! (a *support*) into column `n` the row equation
//! `P_n = sum w * t^rise * P_{n - length}`, `w >= 0`, is one exact linear
//! feasibility problem (Gaussian elimination and an exact phase-one simplex;
//! a bounded depth-first search for integer weights).  The search goes
//! column by column and only tries *inclusion-minimal* feasible supports that
//! are planar together with the edges already chosen: planarity is preserved
//! under removing edges, and every feasible support contains a minimal one,
//! so this loses no solutions, and every support it keeps is planar, so it
//! finds no false ones.  Failed states (a column plus the supports of the
//! previous `max_step_length - 1` columns, the only edges that can meet later
//! edges) are memoized.  Answers are [`PfNetworkRecognition::Certified`]
//! (re-verified by [`build_pf_rows_certificate`]),
//! [`PfNetworkRecognition::ImpossibleWithinBounds`] (exact for the stated
//! bounds) or [`PfNetworkRecognition::Inconclusive`] (a budget ran out).  A
//! row that is not real-rooted (or has a negative coefficient) rules out every
//! network, whatever the bounds.  Rises are nonnegative, so paths never go
//! down.
//!
//! For a different and often stronger PF certificate for the rows of a
//! triangle, see [`crate::left_production`] (Chen--Fu--Ruan).
//!
//! Conventions: rows are zero-based (`P_0` is row `0`); coefficients are
//! ascending; weights are exact [`BigRational`] values.

use crate::brenti_recognition::{solve_nonnegative_column, BrentiWeightDomain, ColumnOutcome};
use crate::real_rootedness::is_real_rooted_bigint_coeffs;
use crate::tnn_network::BigRational;
use crate::Polynomial;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};
use std::collections::HashMap;
use std::fmt;

type Q = BigRational;

/// A step type: edges `(x, y) -> (x + length, y + rise)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct PfStepType {
    /// Number of columns spanned (at least 1).
    pub length: usize,
    /// Height gained (nonnegative).
    pub rise: usize,
}

/// A weighted step into a fixed column of a [`HeightInvariantNetwork`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PfStep {
    /// Number of columns spanned: the step into column `n` starts in column
    /// `n - length`.
    pub length: usize,
    /// Height gained.
    pub rise: usize,
    /// Nonnegative weight, the same at every height.
    pub weight: BigRational,
}

/// A planar network on `{0, ..., N} x Z` that is invariant under height
/// shifts.  See the module documentation for what a valid network certifies.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct HeightInvariantNetwork {
    /// `steps[n - 1]` lists the steps ending in column `n`, for
    /// `n = 1..=steps.len()`; each stands for one edge per height.
    pub steps: Vec<Vec<PfStep>>,
}

/// Errors from the structural checks and certificate matching of
/// [`HeightInvariantNetwork`].  Columns are target columns (`1..=N`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PfNetworkError {
    /// A step has length zero.
    ZeroLengthStep {
        /// Target column.
        column: usize,
    },
    /// A step into column `column` would start before column `0`.
    StepBeforeSource {
        /// Target column.
        column: usize,
        /// Step length (larger than `column`).
        length: usize,
    },
    /// A step has a negative weight.
    NegativeWeight {
        /// Target column.
        column: usize,
        /// The step type.
        step: PfStepType,
        /// The negative weight.
        weight: BigRational,
    },
    /// The same step type is listed twice for one column (merge the weights).
    DuplicateStep {
        /// Target column.
        column: usize,
        /// The repeated step type.
        step: PfStepType,
    },
    /// A positive-weight step is not a primitive lattice segment, so its edges
    /// pass through lattice points.
    StepThroughVertex {
        /// Target column.
        column: usize,
        /// The offending step type (`gcd(length, rise) > 1`).
        step: PfStepType,
    },
    /// Two positive-weight edges meet away from a common endpoint.  The first
    /// edge starts at height `0`, the second at height `height_offset`.
    Crossing {
        /// Target column of the first edge.
        first_column: usize,
        /// Step type of the first edge.
        first_step: PfStepType,
        /// Target column of the second edge.
        second_column: usize,
        /// Step type of the second edge.
        second_step: PfStepType,
        /// Starting height of the second edge relative to the first.
        height_offset: i64,
    },
    /// A requested row lies beyond the last column of the network.
    RowOutOfRange {
        /// The requested row.
        row: usize,
        /// The last column `N`.
        last_column: usize,
    },
    /// The number of polynomials differs from the number of row indices.
    RowIndexCountMismatch {
        /// Number of supplied polynomials.
        num_polynomials: usize,
        /// Number of supplied row indices.
        num_row_indices: usize,
    },
    /// A network row differs from the corresponding polynomial.
    SequenceMismatch {
        /// Position in the supplied polynomial list.
        sequence_index: usize,
        /// The network row compared.
        row_index: usize,
        /// Coefficients of the supplied polynomial (ascending).
        expected: Vec<BigRational>,
        /// Coefficients of the network row (ascending, trimmed).
        found: Vec<BigRational>,
    },
    /// The explicitly rebuilt path matrix differs from the Toeplitz window of
    /// the row (an internal inconsistency).
    ToeplitzMismatch {
        /// Position in the supplied polynomial list.
        sequence_index: usize,
        /// The network row compared.
        row_index: usize,
        /// Source index `i` of the first differing entry.
        source: usize,
        /// Sink index `j` of the first differing entry.
        sink: usize,
    },
    /// The rows stored in a [`PfRowsCertificate`] differ from the rows the
    /// network produces.
    StoredRowsMismatch,
}

impl fmt::Display for PfNetworkError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroLengthStep { column } => {
                write!(f, "a step into column {} has length 0", column)
            }
            Self::StepBeforeSource { column, length } => write!(
                f,
                "a step of length {} into column {} would start before column 0",
                length, column
            ),
            Self::NegativeWeight {
                column,
                step,
                weight,
            } => write!(
                f,
                "step (length {}, rise {}) into column {} has negative weight {}",
                step.length, step.rise, column, weight
            ),
            Self::DuplicateStep { column, step } => write!(
                f,
                "step (length {}, rise {}) is listed twice for column {}",
                step.length, step.rise, column
            ),
            Self::StepThroughVertex { column, step } => write!(
                f,
                "step (length {}, rise {}) into column {} passes through a lattice point \
                 (gcd(length, rise) > 1)",
                step.length, step.rise, column
            ),
            Self::Crossing {
                first_column,
                first_step,
                second_column,
                second_step,
                height_offset,
            } => write!(
                f,
                "nonplanar: step (length {}, rise {}) into column {} meets step (length {}, \
                 rise {}) into column {} started {} higher, away from a common endpoint",
                first_step.length,
                first_step.rise,
                first_column,
                second_step.length,
                second_step.rise,
                second_column,
                height_offset
            ),
            Self::RowOutOfRange { row, last_column } => write!(
                f,
                "row {} is beyond the last column {} of the network",
                row, last_column
            ),
            Self::RowIndexCountMismatch {
                num_polynomials,
                num_row_indices,
            } => write!(
                f,
                "got {} polynomials but {} row indices",
                num_polynomials, num_row_indices
            ),
            Self::SequenceMismatch {
                sequence_index,
                row_index,
                expected,
                found,
            } => write!(
                f,
                "polynomial {} does not match network row {}: expected {:?}, found {:?}",
                sequence_index, row_index, expected, found
            ),
            Self::ToeplitzMismatch {
                sequence_index,
                row_index,
                source,
                sink,
            } => write!(
                f,
                "polynomial {} (network row {}): rebuilt path matrix differs from the \
                 Toeplitz window at source {}, sink {}",
                sequence_index, row_index, source, sink
            ),
            Self::StoredRowsMismatch => write!(
                f,
                "the rows stored in the certificate differ from the network rows"
            ),
        }
    }
}

impl std::error::Error for PfNetworkError {}

/// A positive-weight edge used by the planarity check.
#[derive(Clone, Copy)]
struct Edge {
    column: usize,
    step: PfStepType,
}

impl HeightInvariantNetwork {
    /// The last column `N` (the network has rows `0..=N`).
    pub fn last_column(&self) -> usize {
        self.steps.len()
    }

    /// The steps ending in column `n` (empty for `n = 0` or `n > N`).
    pub fn steps_into(&self, n: usize) -> &[PfStep] {
        n.checked_sub(1)
            .and_then(|i| self.steps.get(i))
            .map_or(&[], Vec::as_slice)
    }

    /// Exact structural validation: lengths, sources, nonnegative weights,
    /// distinct step types per column, primitive positive-weight steps, and
    /// planarity of the positive-weight edges.
    pub fn validate(&self) -> Result<(), PfNetworkError> {
        let mut edges = Vec::new();
        for (index, steps) in self.steps.iter().enumerate() {
            let column = index + 1;
            let mut seen = Vec::with_capacity(steps.len());
            for step in steps {
                let kind = PfStepType {
                    length: step.length,
                    rise: step.rise,
                };
                if step.length == 0 {
                    return Err(PfNetworkError::ZeroLengthStep { column });
                }
                if step.length > column {
                    return Err(PfNetworkError::StepBeforeSource {
                        column,
                        length: step.length,
                    });
                }
                if step.weight < Q::zero() {
                    return Err(PfNetworkError::NegativeWeight {
                        column,
                        step: kind,
                        weight: step.weight.clone(),
                    });
                }
                if seen.contains(&kind) {
                    return Err(PfNetworkError::DuplicateStep { column, step: kind });
                }
                seen.push(kind);
                if step.weight.is_zero() {
                    continue;
                }
                if !is_primitive(kind) {
                    return Err(PfNetworkError::StepThroughVertex { column, step: kind });
                }
                edges.push(Edge { column, step: kind });
            }
        }
        for (i, first) in edges.iter().enumerate() {
            for second in &edges[i..] {
                let same = first.column == second.column && first.step == second.step;
                let first_start = (first.column - first.step.length) as i64;
                let second_start = (second.column - second.step.length) as i64;
                if let Some(height_offset) =
                    find_conflict(first_start, first.step, second_start, second.step, same)
                {
                    return Err(PfNetworkError::Crossing {
                        first_column: first.column,
                        first_step: first.step,
                        second_column: second.column,
                        second_step: second.step,
                        height_offset,
                    });
                }
            }
        }
        Ok(())
    }

    fn check_row(&self, n: usize) -> Result<(), PfNetworkError> {
        if n > self.last_column() {
            return Err(PfNetworkError::RowOutOfRange {
                row: n,
                last_column: self.last_column(),
            });
        }
        Ok(())
    }

    /// Coefficients `M_{n,k}` of row `n` (ascending, trimmed; the zero row is
    /// empty), computed by the column recurrence.  The network is validated
    /// first.
    pub fn row_coefficients(&self, n: usize) -> Result<Vec<BigRational>, PfNetworkError> {
        self.validate()?;
        self.check_row(n)?;
        Ok(self.rows_unchecked(n).swap_remove(n))
    }

    /// Row `n` as a polynomial.
    pub fn row_polynomial(&self, n: usize) -> Result<Polynomial<BigRational>, PfNetworkError> {
        Ok(Polynomial::new(self.row_coefficients(n)?))
    }

    /// Rows `0..=last_row`.
    pub fn row_polynomials_prefix(
        &self,
        last_row: usize,
    ) -> Result<Vec<Polynomial<BigRational>>, PfNetworkError> {
        self.validate()?;
        self.check_row(last_row)?;
        Ok(self
            .rows_unchecked(last_row)
            .into_iter()
            .map(Polynomial::new)
            .collect())
    }

    fn rows_unchecked(&self, last: usize) -> Vec<Vec<Q>> {
        let mut rows: Vec<Vec<Q>> = vec![vec![Q::one()]];
        for n in 1..=last {
            let mut row: Vec<Q> = Vec::new();
            for step in self.steps_into(n) {
                if step.weight.is_zero() {
                    continue;
                }
                let previous = &rows[n - step.length];
                if row.len() < previous.len() + step.rise {
                    row.resize(previous.len() + step.rise, Q::zero());
                }
                for (k, c) in previous.iter().enumerate() {
                    row[k + step.rise] += c * &step.weight;
                }
            }
            while row.last().is_some_and(Zero::is_zero) {
                row.pop();
            }
            rows.push(row);
        }
        rows
    }

    /// The path matrix from the sources `(0, i)` to the sinks `(n, j)`,
    /// `0 <= i, j < size`, rebuilt by dynamic programming over the explicit
    /// lattice points (it does not use height invariance).  For a valid
    /// network it equals [`toeplitz_window`] of row `n`.
    pub fn toeplitz_path_matrix(
        &self,
        n: usize,
        size: usize,
    ) -> Result<Vec<Vec<BigRational>>, PfNetworkError> {
        self.validate()?;
        self.check_row(n)?;
        let mut matrix = Vec::with_capacity(size);
        for source in 0..size {
            // weights[x][y]: total weight of paths from (0, source) to (x, y),
            // for heights y < size (rises are nonnegative, so a path above
            // height size - 1 never returns to a sink).
            let mut weights: Vec<Vec<Q>> = vec![vec![Q::zero(); size]; n + 1];
            weights[0][source] = Q::one();
            for x in 1..=n {
                for y in 0..size {
                    let mut total = Q::zero();
                    for step in self.steps_into(x) {
                        if step.rise > y {
                            continue;
                        }
                        let from = &weights[x - step.length][y - step.rise];
                        if !from.is_zero() {
                            total += from * &step.weight;
                        }
                    }
                    weights[x][y] = total;
                }
            }
            matrix.push(weights.swap_remove(n));
        }
        Ok(matrix)
    }
}

/// The upper-triangular Toeplitz window `T[i][j] = a_{j - i}` (zero for
/// `j < i` or `j - i >= coeffs.len()`), `0 <= i, j < size`.  It is the
/// transpose of the lower Toeplitz section `[a_{i - j}]`; a sequence is PF
/// iff all these windows are totally nonnegative.
pub fn toeplitz_window(coeffs: &[BigRational], size: usize) -> Vec<Vec<BigRational>> {
    (0..size)
        .map(|i| {
            (0..size)
                .map(|j| {
                    j.checked_sub(i)
                        .and_then(|d| coeffs.get(d))
                        .cloned()
                        .unwrap_or_else(Q::zero)
                })
                .collect()
        })
        .collect()
}

fn is_primitive(step: PfStepType) -> bool {
    step.length.gcd(&step.rise) == 1
}

/// Whether two step types can coexist at every column and every height
/// without a crossing (both must be primitive).  A uniform step pattern is
/// planar iff its step types are pairwise compatible.
pub fn step_types_compatible(first: PfStepType, second: PfStepType) -> bool {
    if first.length == 0 || second.length == 0 || !is_primitive(first) || !is_primitive(second) {
        return false;
    }
    let s2 = second.length as i64;
    let s1 = first.length as i64;
    (-(s2 - 1)..s1).all(|start| {
        let same = first == second && start == 0;
        find_conflict(0, first, start, second, same).is_none()
    })
}

/// The first height offset `b` at which the edge of `second` starting at
/// `(second_start, b)` meets the edge of `first` starting at
/// `(first_start, 0)` away from a common endpoint.  `same` skips `b = 0` (the
/// edge itself).
fn find_conflict(
    first_start: i64,
    first: PfStepType,
    second_start: i64,
    second: PfStepType,
    same: bool,
) -> Option<i64> {
    let first_end = first_start + first.length as i64;
    let second_end = second_start + second.length as i64;
    if first_start >= second_end || second_start >= first_end {
        return None;
    }
    let d1 = first.rise as i64;
    let d2 = second.rise as i64;
    let p1 = (first_start, 0);
    let p2 = (first_end, d1);
    (-d2..=d1).find(|&b| {
        if same && b == 0 {
            return false;
        }
        segments_conflict(p1, p2, (second_start, b), (second_end, b + d2))
    })
}

type Point = (i64, i64);

fn orient(a: Point, b: Point, c: Point) -> i128 {
    let (ax, ay) = (a.0 as i128, a.1 as i128);
    (b.0 as i128 - ax) * (c.1 as i128 - ay) - (b.1 as i128 - ay) * (c.0 as i128 - ax)
}

fn within(a: Point, b: Point, c: Point) -> bool {
    c.0 >= a.0.min(b.0) && c.0 <= a.0.max(b.0) && c.1 >= a.1.min(b.1) && c.1 <= a.1.max(b.1)
}

/// Whether the closed segments `[p1, p2]` and `[q1, q2]` share a point that is
/// not a common endpoint of both.
fn segments_conflict(p1: Point, p2: Point, q1: Point, q2: Point) -> bool {
    let d1 = orient(q1, q2, p1);
    let d2 = orient(q1, q2, p2);
    let d3 = orient(p1, p2, q1);
    let d4 = orient(p1, p2, q2);
    if d1 == 0 && d2 == 0 {
        // Collinear: conflict iff the overlap has positive length.  All edges
        // have positive x-extent, so compare x-intervals.
        let lo = p1.0.min(p2.0).max(q1.0.min(q2.0));
        let hi = p1.0.max(p2.0).min(q1.0.max(q2.0));
        return lo < hi;
    }
    if ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0)) {
        return true;
    }
    let touches = |d: i128, a: Point, b: Point, c: Point, o1: Point, o2: Point| {
        d == 0 && within(a, b, c) && c != o1 && c != o2
    };
    touches(d1, q1, q2, p1, q1, q2)
        || touches(d2, q1, q2, p2, q1, q2)
        || touches(d3, p1, p2, q1, p1, p2)
        || touches(d4, p1, p2, q2, p1, p2)
}

/// A verified match between a [`HeightInvariantNetwork`] and polynomial rows.
///
/// It certifies that each matched polynomial is a Pólya frequency sequence
/// (real-rooted with nonpositive zeros), for the matched rows only.  It does
/// not certify total nonnegativity of the coefficient matrix of the rows.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PfRowsCertificate {
    /// The validated network.
    pub network: HeightInvariantNetwork,
    /// `row_indices[i]` is the network row (zero-based) equal to the `i`-th
    /// supplied polynomial.
    pub row_indices: Vec<usize>,
    /// The matched rows, as computed from the network.
    pub matched_rows: Vec<Polynomial<BigRational>>,
}

impl PfRowsCertificate {
    /// Re-run every check of [`build_pf_rows_certificate`] against `polys`.
    pub fn verify_against_bigint_polynomials(
        &self,
        polys: &[Polynomial<BigInt>],
    ) -> Result<(), PfNetworkError> {
        let rebuilt = build_pf_rows_certificate(self.network.clone(), &self.row_indices, polys)?;
        if rebuilt.matched_rows != self.matched_rows {
            return Err(PfNetworkError::StoredRowsMismatch);
        }
        Ok(())
    }
}

/// Validate `network`, check that row `row_indices[i]` equals `polys[i]`
/// exactly, and check that the explicitly rebuilt path matrix
/// ([`HeightInvariantNetwork::toeplitz_path_matrix`]) equals the Toeplitz
/// window of each matched row, of size `deg + 1`.
pub fn build_pf_rows_certificate(
    network: HeightInvariantNetwork,
    row_indices: &[usize],
    polys: &[Polynomial<BigInt>],
) -> Result<PfRowsCertificate, PfNetworkError> {
    network.validate()?;
    if polys.len() != row_indices.len() {
        return Err(PfNetworkError::RowIndexCountMismatch {
            num_polynomials: polys.len(),
            num_row_indices: row_indices.len(),
        });
    }
    let last = row_indices.iter().copied().max().unwrap_or(0);
    network.check_row(last)?;
    let rows = network.rows_unchecked(last);
    let mut matched_rows = Vec::with_capacity(polys.len());
    for (sequence_index, (poly, &row_index)) in polys.iter().zip(row_indices).enumerate() {
        let expected: Vec<Q> = poly
            .coeffs()
            .iter()
            .map(|c| Q::from_integer(c.clone()))
            .collect();
        let found = &rows[row_index];
        if *found != expected {
            return Err(PfNetworkError::SequenceMismatch {
                sequence_index,
                row_index,
                expected,
                found: found.clone(),
            });
        }
        let size = found.len().max(1);
        let rebuilt = network.toeplitz_path_matrix(row_index, size)?;
        let window = toeplitz_window(found, size);
        for (i, (got, want)) in rebuilt.iter().zip(&window).enumerate() {
            if let Some(j) = (0..size).find(|&j| got[j] != want[j]) {
                return Err(PfNetworkError::ToeplitzMismatch {
                    sequence_index,
                    row_index,
                    source: i,
                    sink: j,
                });
            }
        }
        matched_rows.push(Polynomial::new(found.clone()));
    }
    Ok(PfRowsCertificate {
        network,
        row_indices: row_indices.to_vec(),
        matched_rows,
    })
}

/// Bounds and budgets for [`recognize_pf_network`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PfNetworkOptions {
    /// Network row of the first supplied polynomial: `1` (default) means the
    /// input is `P_1, P_2, ...` with `P_0 = 1` implicit; `0` means the input
    /// starts with `P_0`, which must equal `1`.
    pub first_row_index: usize,
    /// Largest step length (at least 1; default 3).
    pub max_step_length: usize,
    /// Largest rise (default 3).
    pub max_rise: usize,
    /// Rational or integer weights.
    pub weights: BrentiWeightDomain,
    /// Budget on the number of column supports visited by the search.
    pub max_search_nodes: usize,
    /// Budget on the integer search for one row problem.
    pub max_integer_search_nodes: usize,
}

impl Default for PfNetworkOptions {
    fn default() -> Self {
        Self {
            first_row_index: 1,
            max_step_length: 3,
            max_rise: 3,
            weights: BrentiWeightDomain::Rational,
            max_search_nodes: 1_000_000,
            max_integer_search_nodes: 10_000,
        }
    }
}

/// Outcome of [`recognize_pf_network`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PfNetworkRecognition {
    /// A height-invariant planar network reproduces every supplied row, so
    /// every supplied row is a PF sequence (real-rooted, nonpositive zeros).
    Certified {
        /// The verified certificate (network with zero weights dropped).
        certificate: PfRowsCertificate,
    },
    /// No network within the bounds reproduces the rows (exact).  When the
    /// obstruction is a row that is not real-rooted, no network exists at all.
    ImpossibleWithinBounds {
        /// Number of column supports visited by the search.
        supports_checked: usize,
        /// Human-readable reason.
        obstruction: String,
    },
    /// A budget was exhausted before the search was decided.
    Inconclusive {
        /// Which budget was exhausted.
        reason: String,
    },
}

impl PfNetworkRecognition {
    /// Whether a certificate was found.
    pub fn is_certified(&self) -> bool {
        matches!(self, Self::Certified { .. })
    }
}

impl fmt::Display for PfNetworkRecognition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Certified { certificate } => write!(
                f,
                "certified: rows {:?} are rows of a height-invariant planar network, so each \
                 is a Pólya frequency sequence (real-rooted, nonpositive zeros); this prefix only",
                certificate.row_indices
            ),
            Self::ImpossibleWithinBounds {
                supports_checked,
                obstruction,
            } => write!(
                f,
                "impossible within bounds ({} column supports checked): {}",
                supports_checked, obstruction
            ),
            Self::Inconclusive { reason } => write!(f, "inconclusive: {}", reason),
        }
    }
}

/// Invalid input for [`recognize_pf_network`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PfNetworkRecognitionError {
    /// No rows were supplied.
    EmptyInput,
    /// `first_row_index` was not 0 or 1.
    UnsupportedFirstRowIndex {
        /// The rejected index.
        found: usize,
    },
    /// With `first_row_index == 0` the first row must be the polynomial `1`.
    FirstRowNotOne {
        /// Coefficients of the supplied first row (ascending).
        found: Vec<BigInt>,
    },
    /// `max_step_length` was 0.
    ZeroStepLength,
    /// The network found by the search failed re-verification (a bug).
    CertificateRejected(PfNetworkError),
}

impl fmt::Display for PfNetworkRecognitionError {
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
            Self::ZeroStepLength => write!(f, "max_step_length must be at least 1"),
            Self::CertificateRejected(error) => {
                write!(f, "internal error: certificate rejected: {}", error)
            }
        }
    }
}

impl std::error::Error for PfNetworkRecognitionError {}

/// All primitive step types within the bounds, in increasing order.
pub fn primitive_step_types(max_step_length: usize, max_rise: usize) -> Vec<PfStepType> {
    let mut types = Vec::new();
    for length in 1..=max_step_length {
        for rise in 0..=max_rise {
            let step = PfStepType { length, rise };
            if is_primitive(step) {
                types.push(step);
            }
        }
    }
    types
}

/// Search for a height-invariant planar network whose rows are the given
/// polynomials.  See the module documentation for the setting, the method and
/// what a certificate proves.
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::pf_network::{recognize_pf_network, PfNetworkOptions, PfNetworkRecognition};
/// use polytool::Polynomial;
///
/// let rows: Vec<Polynomial<BigInt>> = [vec![1, 1], vec![1, 2, 1], vec![1, 3, 3, 1]]
///     .iter()
///     .map(|r| Polynomial::from_i64_coeffs(r))
///     .collect();
/// let result = recognize_pf_network(&rows, &PfNetworkOptions::default()).unwrap();
/// assert!(result.is_certified());
///
/// // 1 + t + t^2 is not real-rooted, so no network exists at any bound.
/// let bad = vec![Polynomial::from_i64_coeffs(&[1, 1, 1])];
/// let result = recognize_pf_network(&bad, &PfNetworkOptions::default()).unwrap();
/// assert!(matches!(result, PfNetworkRecognition::ImpossibleWithinBounds { .. }));
/// ```
pub fn recognize_pf_network(
    rows: &[Polynomial<BigInt>],
    options: &PfNetworkOptions,
) -> Result<PfNetworkRecognition, PfNetworkRecognitionError> {
    if rows.is_empty() {
        return Err(PfNetworkRecognitionError::EmptyInput);
    }
    if options.first_row_index > 1 {
        return Err(PfNetworkRecognitionError::UnsupportedFirstRowIndex {
            found: options.first_row_index,
        });
    }
    if options.max_step_length == 0 {
        return Err(PfNetworkRecognitionError::ZeroStepLength);
    }
    let mut vectors: Vec<Vec<BigInt>> = Vec::with_capacity(rows.len() + 1);
    if options.first_row_index == 1 {
        vectors.push(vec![BigInt::one()]);
    } else if rows[0].coeffs() != [BigInt::one()] {
        return Err(PfNetworkRecognitionError::FirstRowNotOne {
            found: rows[0].coeffs().to_vec(),
        });
    }
    vectors.extend(rows.iter().map(|p| p.coeffs().to_vec()));

    for (row, v) in vectors.iter().enumerate() {
        if let Some((k, c)) = v.iter().enumerate().find(|(_, c)| **c < BigInt::zero()) {
            return Ok(PfNetworkRecognition::ImpossibleWithinBounds {
                supports_checked: 0,
                obstruction: format!(
                    "row {} has negative coefficient {} at t^{}; path weights are nonnegative",
                    row, c, k
                ),
            });
        }
    }
    for (row, v) in vectors.iter().enumerate() {
        if !v.is_empty() && !is_real_rooted_bigint_coeffs(v) {
            return Ok(PfNetworkRecognition::ImpossibleWithinBounds {
                supports_checked: 0,
                obstruction: format!(
                    "row {} is not real-rooted, but every row of a height-invariant planar \
                     network is a PF sequence (no network exists at any bound)",
                    row
                ),
            });
        }
    }

    let row_indices: Vec<usize> = (0..rows.len())
        .map(|i| i + options.first_row_index)
        .collect();
    let mut search = SupportSearch::new(&vectors, options);
    match search.run() {
        SearchOutcome::Found(columns) => {
            let network = HeightInvariantNetwork {
                steps: columns
                    .into_iter()
                    .map(|column| {
                        column
                            .into_iter()
                            .filter(|(_, w)| !w.is_zero())
                            .map(|(step, weight)| PfStep {
                                length: step.length,
                                rise: step.rise,
                                weight,
                            })
                            .collect()
                    })
                    .collect(),
            };
            let certificate = build_pf_rows_certificate(network, &row_indices, rows)
                .map_err(PfNetworkRecognitionError::CertificateRejected)?;
            Ok(PfNetworkRecognition::Certified { certificate })
        }
        SearchOutcome::Impossible(obstruction) => {
            Ok(PfNetworkRecognition::ImpossibleWithinBounds {
                supports_checked: search.nodes,
                obstruction,
            })
        }
        SearchOutcome::Inconclusive(reason) => Ok(PfNetworkRecognition::Inconclusive { reason }),
    }
}

type ColumnWeights = Vec<(PfStepType, Q)>;

enum SearchOutcome {
    Found(Vec<ColumnWeights>),
    Impossible(String),
    Inconclusive(String),
}

/// A minimal feasible support of one row problem, with its weights.
#[derive(Clone)]
struct Support {
    types: Vec<usize>,
    weights: Vec<Q>,
}

struct SupportSearch<'a> {
    vectors: &'a [Vec<BigInt>],
    types: Vec<PfStepType>,
    max_length: usize,
    integer: bool,
    max_nodes: usize,
    max_integer_nodes: usize,
    lp_cache: HashMap<(usize, Vec<usize>), ColumnOutcome>,
    minimal_cache: HashMap<(usize, Vec<usize>), Vec<Support>>,
    failed: std::collections::HashSet<(usize, Vec<Vec<usize>>)>,
    nodes: usize,
    budget_exhausted: bool,
    inconclusive_rows: usize,
}

impl<'a> SupportSearch<'a> {
    fn new(vectors: &'a [Vec<BigInt>], options: &PfNetworkOptions) -> Self {
        Self {
            vectors,
            types: primitive_step_types(options.max_step_length, options.max_rise),
            max_length: options.max_step_length,
            integer: options.weights == BrentiWeightDomain::Integer,
            max_nodes: options.max_search_nodes,
            max_integer_nodes: options.max_integer_search_nodes,
            lp_cache: HashMap::new(),
            minimal_cache: HashMap::new(),
            failed: std::collections::HashSet::new(),
            nodes: 0,
            budget_exhausted: false,
            inconclusive_rows: 0,
        }
    }

    fn last_row(&self) -> usize {
        self.vectors.len() - 1
    }

    fn run(&mut self) -> SearchOutcome {
        for n in 1..=self.last_row() {
            let usable: Vec<usize> = (0..self.types.len())
                .filter(|&i| self.types[i].length <= n)
                .collect();
            if self.row(n, &usable) == ColumnOutcome::Infeasible {
                return SearchOutcome::Impossible(format!(
                    "row {} is not a nonnegative combination of t^rise * P_(n - length) over \
                     the step types within the bounds, even ignoring planarity",
                    n
                ));
            }
        }
        let mut chosen: Vec<Support> = Vec::with_capacity(self.last_row());
        if self.extend(&mut chosen) {
            let columns = chosen
                .iter()
                .map(|support| {
                    support
                        .types
                        .iter()
                        .map(|&i| self.types[i])
                        .zip(support.weights.iter().cloned())
                        .collect()
                })
                .collect();
            return SearchOutcome::Found(columns);
        }
        if self.budget_exhausted {
            SearchOutcome::Inconclusive(format!(
                "support search budget of {} nodes exhausted",
                self.max_nodes
            ))
        } else if self.inconclusive_rows > 0 {
            SearchOutcome::Inconclusive(format!(
                "integer weight search budget of {} nodes exhausted for some row",
                self.max_integer_nodes
            ))
        } else {
            let domain = if self.integer { "integer" } else { "rational" };
            SearchOutcome::Impossible(format!(
                "no planar choice of steps with length <= {} and rise <= {} admits \
                 nonnegative {} weights",
                self.max_length,
                self.types.iter().map(|t| t.rise).max().unwrap_or(0),
                domain
            ))
        }
    }

    /// Depth-first search over the columns; `chosen[m - 1]` is the support
    /// into column `m`.  Only minimal feasible supports are tried: replacing a
    /// support by a feasible subset keeps its row and removes obstacles for
    /// later columns, so this loses no solutions.
    fn extend(&mut self, chosen: &mut Vec<Support>) -> bool {
        let n = chosen.len() + 1;
        if n > self.last_row() {
            return true;
        }
        // Only edges into columns m > n - max_length can meet later edges.
        let first_relevant = (n + 1).saturating_sub(self.max_length).max(1);
        let key = (
            n,
            chosen[first_relevant - 1..]
                .iter()
                .map(|s| s.types.clone())
                .collect::<Vec<_>>(),
        );
        if self.failed.contains(&key) {
            return false;
        }
        let candidates: Vec<usize> = (0..self.types.len())
            .filter(|&i| {
                let step = self.types[i];
                step.length <= n
                    && chosen[first_relevant - 1..]
                        .iter()
                        .enumerate()
                        .all(|(offset, support)| {
                            let m = first_relevant + offset;
                            support.types.iter().all(|&j| {
                                let other = self.types[j];
                                find_conflict(
                                    (m - other.length) as i64,
                                    other,
                                    (n - step.length) as i64,
                                    step,
                                    false,
                                )
                                .is_none()
                            })
                        })
            })
            .collect();
        let inconclusive_before = self.inconclusive_rows;
        let supports = self.minimal_supports(n, candidates);
        for support in supports {
            if self.nodes >= self.max_nodes {
                self.budget_exhausted = true;
                return false;
            }
            self.nodes += 1;
            chosen.push(support);
            if self.extend(chosen) {
                return true;
            }
            chosen.pop();
            if self.budget_exhausted {
                return false;
            }
        }
        if self.inconclusive_rows == inconclusive_before {
            self.failed.insert(key);
        }
        false
    }

    /// The inclusion-minimal subsets of `candidates` that are planar within
    /// column `n` and admit nonnegative weights for row `n`.
    fn minimal_supports(&mut self, n: usize, candidates: Vec<usize>) -> Vec<Support> {
        let key = (n, candidates);
        if let Some(found) = self.minimal_cache.get(&key) {
            return found.clone();
        }
        let candidates = &key.1;
        let mut cliques: Vec<Vec<usize>> = vec![Vec::new()];
        for &i in candidates {
            let extended: Vec<Vec<usize>> = cliques
                .iter()
                .filter(|clique| clique.iter().all(|&j| self.same_column_compatible(n, i, j)))
                .map(|clique| {
                    let mut next = clique.clone();
                    next.push(i);
                    next
                })
                .collect();
            cliques.extend(extended);
            if cliques.len() > self.max_nodes {
                self.budget_exhausted = true;
                break;
            }
        }
        cliques.sort_by_key(Vec::len);
        let mut minimal: Vec<Support> = Vec::new();
        for clique in cliques {
            if minimal
                .iter()
                .any(|m| m.types.iter().all(|t| clique.contains(t)))
            {
                continue;
            }
            match self.row(n, &clique) {
                ColumnOutcome::Feasible(weights) => minimal.push(Support {
                    types: clique,
                    weights,
                }),
                ColumnOutcome::Infeasible => {}
                ColumnOutcome::Inconclusive => self.inconclusive_rows += 1,
            }
        }
        self.minimal_cache.insert(key, minimal.clone());
        minimal
    }

    fn same_column_compatible(&self, n: usize, i: usize, j: usize) -> bool {
        let (a, b) = (self.types[i], self.types[j]);
        find_conflict((n - a.length) as i64, a, (n - b.length) as i64, b, i == j).is_none()
    }

    /// Decide `P_n = sum_i w_i t^{rise_i} P_{n - length_i}`, `w >= 0`, over
    /// the step types `usable` (all with `length <= n`).
    fn row(&mut self, n: usize, usable: &[usize]) -> ColumnOutcome {
        let key = (n, usable.to_vec());
        if let Some(outcome) = self.lp_cache.get(&key) {
            return outcome.clone();
        }
        let target = &self.vectors[n];
        let mut height = target.len();
        for &i in usable {
            let step = self.types[i];
            height = height.max(self.vectors[n - step.length].len() + step.rise);
        }
        let matrix: Vec<Vec<Q>> = (0..height)
            .map(|k| {
                usable
                    .iter()
                    .map(|&i| {
                        let step = self.types[i];
                        k.checked_sub(step.rise)
                            .and_then(|j| self.vectors[n - step.length].get(j))
                            .map_or_else(Q::zero, |c| Q::from_integer(c.clone()))
                    })
                    .collect()
            })
            .collect();
        let rhs: Vec<Q> = (0..height)
            .map(|k| {
                target
                    .get(k)
                    .map_or_else(Q::zero, |c| Q::from_integer(c.clone()))
            })
            .collect();
        let outcome = solve_nonnegative_column(
            &matrix,
            &rhs,
            usable.len(),
            self.integer,
            self.max_integer_nodes,
        );
        self.lp_cache.insert(key, outcome.clone());
        outcome
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numeric::binomial_bigint;

    fn q(v: i64) -> Q {
        Q::from_integer(BigInt::from(v))
    }

    fn step(length: usize, rise: usize, weight: i64) -> PfStep {
        PfStep {
            length,
            rise,
            weight: q(weight),
        }
    }

    fn polys(rows: &[&[i64]]) -> Vec<Polynomial<BigInt>> {
        rows.iter()
            .map(|r| Polynomial::from_i64_coeffs(r))
            .collect()
    }

    fn pascal_rows(last: usize) -> Vec<Polynomial<BigInt>> {
        (1..=last)
            .map(|n| Polynomial::new((0..=n).map(|k| binomial_bigint(n, k)).collect()))
            .collect()
    }

    fn eulerian_rows(last: usize) -> Vec<Polynomial<BigInt>> {
        // A(n, k): permutations of [n + 1] with k descents, rows n = 1..=last.
        let mut rows: Vec<Vec<BigInt>> = vec![vec![BigInt::one()]];
        for m in 1..=last {
            let prev = &rows[m - 1];
            let row = (0..=m)
                .map(|k| {
                    let left = if k >= 1 {
                        BigInt::from(m - k + 1) * &prev[k - 1]
                    } else {
                        BigInt::zero()
                    };
                    let right = prev
                        .get(k)
                        .map_or_else(BigInt::zero, |c| BigInt::from(k + 1) * c);
                    left + right
                })
                .collect();
            rows.push(row);
        }
        rows.into_iter().skip(1).map(Polynomial::new).collect()
    }

    #[test]
    fn segment_conflicts() {
        // Proper crossing.
        assert!(segments_conflict((0, 0), (1, 2), (0, 1), (1, 1)));
        // Common endpoint only.
        assert!(!segments_conflict((0, 0), (1, 1), (0, 0), (1, 0)));
        assert!(!segments_conflict((0, 0), (2, 1), (1, 0), (2, 1)));
        // An endpoint in the interior of the other segment.
        assert!(segments_conflict((0, 0), (2, 2), (1, 1), (2, 1)));
        // Collinear overlap versus collinear touching.
        assert!(segments_conflict((0, 0), (2, 0), (1, 0), (3, 0)));
        assert!(!segments_conflict((0, 0), (1, 0), (1, 0), (2, 0)));
        // Parallel, disjoint.
        assert!(!segments_conflict((0, 0), (1, 1), (0, 1), (1, 2)));
    }

    #[test]
    fn compatibility_of_step_types() {
        let t = |length, rise| PfStepType { length, rise };
        // Brenti's recursive-matrix pattern is planar.
        for (a, b) in [(t(1, 0), t(1, 1)), (t(1, 0), t(2, 1)), (t(1, 1), t(2, 1))] {
            assert!(step_types_compatible(a, b), "{:?} {:?}", a, b);
            assert!(step_types_compatible(b, a), "{:?} {:?}", b, a);
        }
        // Rises 0 and 2 in one column cross.
        assert!(!step_types_compatible(t(1, 0), t(1, 2)));
        // Non-primitive steps are never allowed.
        assert!(!step_types_compatible(t(2, 2), t(1, 0)));
        assert!(!step_types_compatible(t(2, 0), t(2, 0)));
        // Every primitive type is compatible with itself.
        for s in primitive_step_types(4, 4) {
            assert!(step_types_compatible(s, s), "{:?}", s);
        }
    }

    #[test]
    fn pascal_network_rows_and_toeplitz() {
        let network = HeightInvariantNetwork {
            steps: (0..4).map(|_| vec![step(1, 0, 1), step(1, 1, 1)]).collect(),
        };
        network.validate().unwrap();
        assert_eq!(
            network.row_coefficients(4).unwrap(),
            vec![q(1), q(4), q(6), q(4), q(1)]
        );
        let matrix = network.toeplitz_path_matrix(3, 5).unwrap();
        let row = network.row_coefficients(3).unwrap();
        assert_eq!(matrix, toeplitz_window(&row, 5));
        let cert = build_pf_rows_certificate(network, &[1, 2, 3, 4], &pascal_rows(4)).unwrap();
        cert.verify_against_bigint_polynomials(&pascal_rows(4))
            .unwrap();
    }

    #[test]
    fn invalid_networks_rejected() {
        let crossing = HeightInvariantNetwork {
            steps: vec![vec![step(1, 0, 1), step(1, 2, 1)]],
        };
        assert!(matches!(
            crossing.validate(),
            Err(PfNetworkError::Crossing { .. })
        ));
        // The crossing network would produce the row 1 + t^2, which is not PF.
        let through = HeightInvariantNetwork {
            steps: vec![vec![], vec![step(2, 2, 1)]],
        };
        assert!(matches!(
            through.validate(),
            Err(PfNetworkError::StepThroughVertex { .. })
        ));
        let early = HeightInvariantNetwork {
            steps: vec![vec![step(2, 1, 1)]],
        };
        assert!(matches!(
            early.validate(),
            Err(PfNetworkError::StepBeforeSource { .. })
        ));
        let negative = HeightInvariantNetwork {
            steps: vec![vec![step(1, 0, -1)]],
        };
        assert!(matches!(
            negative.validate(),
            Err(PfNetworkError::NegativeWeight { .. })
        ));
        let duplicate = HeightInvariantNetwork {
            steps: vec![vec![step(1, 0, 1), step(1, 0, 2)]],
        };
        assert!(matches!(
            duplicate.validate(),
            Err(PfNetworkError::DuplicateStep { .. })
        ));
        // Zero-weight edges are ignored by the planarity check.
        let zero = HeightInvariantNetwork {
            steps: vec![vec![step(1, 0, 1), step(1, 2, 0)]],
        };
        zero.validate().unwrap();
        // Crossing between columns: a long step over a short one.
        let long = HeightInvariantNetwork {
            steps: vec![vec![step(1, 0, 1)], vec![step(1, 2, 1), step(2, 1, 1)]],
        };
        assert!(matches!(
            long.validate(),
            Err(PfNetworkError::Crossing { .. })
        ));
    }

    #[test]
    fn finite_toeplitz_window_alone_proves_nothing() {
        let to_int = |m: Vec<Vec<Q>>| -> Vec<Vec<BigInt>> {
            m.into_iter()
                .map(|r| r.into_iter().map(|v| v.to_integer()).collect())
                .collect()
        };
        let row = vec![q(1), q(1), q(1)];
        assert!(
            crate::check_totally_nonnegative_neville_bigint(&to_int(toeplitz_window(&row, 3)))
                .is_ok()
        );
        assert!(
            crate::check_totally_nonnegative_neville_bigint(&to_int(toeplitz_window(&row, 4)))
                .is_err()
        );
    }

    #[test]
    fn recognizes_pascal_rows() {
        for weights in [BrentiWeightDomain::Rational, BrentiWeightDomain::Integer] {
            let options = PfNetworkOptions {
                weights,
                ..PfNetworkOptions::default()
            };
            let rows = pascal_rows(6);
            match recognize_pf_network(&rows, &options).unwrap() {
                PfNetworkRecognition::Certified { certificate, .. } => {
                    certificate
                        .verify_against_bigint_polynomials(&rows)
                        .unwrap();
                    assert_eq!(certificate.row_indices, (1..=6).collect::<Vec<_>>());
                }
                other => panic!("expected a certificate, got {}", other),
            }
        }
        // With P_0 given explicitly.
        let mut rows = vec![Polynomial::from_i64_coeffs(&[1])];
        rows.extend(pascal_rows(3));
        let options = PfNetworkOptions {
            first_row_index: 0,
            ..PfNetworkOptions::default()
        };
        assert!(recognize_pf_network(&rows, &options)
            .unwrap()
            .is_certified());
    }

    #[test]
    fn recognizes_three_term_recurrence_rows() {
        // P_n = (1 + t) P_{n-1} + 2 t P_{n-2}: Brenti's recursive pattern.
        let mut rows: Vec<Vec<i64>> = vec![vec![1], vec![1, 1]];
        for n in 2..7 {
            let mut next = vec![0i64; n + 1];
            for (k, c) in rows[n - 1].iter().enumerate() {
                next[k] += c;
                next[k + 1] += c;
            }
            for (k, c) in rows[n - 2].iter().enumerate() {
                next[k + 1] += 2 * c;
            }
            rows.push(next);
        }
        let input: Vec<Polynomial<BigInt>> = rows[1..]
            .iter()
            .map(|r| Polynomial::from_i64_coeffs(r))
            .collect();
        let result = recognize_pf_network(&input, &PfNetworkOptions::default()).unwrap();
        assert!(result.is_certified(), "{}", result);
        // Only single steps: the rows are not products of linear factors.
        let single = PfNetworkOptions {
            max_step_length: 1,
            ..PfNetworkOptions::default()
        };
        let result = recognize_pf_network(&input, &single).unwrap();
        assert!(
            matches!(result, PfNetworkRecognition::ImpossibleWithinBounds { .. }),
            "{}",
            result
        );
    }

    #[test]
    fn never_certifies_one_plus_t_plus_t_squared() {
        for rows in [
            polys(&[&[1, 1, 1]]),
            polys(&[&[1, 1], &[1, 1, 1]]),
            polys(&[&[1, 2], &[1, 1, 1]]),
        ] {
            let result = recognize_pf_network(&rows, &PfNetworkOptions::default()).unwrap();
            assert!(
                matches!(result, PfNetworkRecognition::ImpossibleWithinBounds { .. }),
                "{}",
                result
            );
        }
        // The support search alone (without the real-rootedness shortcut)
        // also finds nothing, even with long steps and large rises.
        let inputs: [&[&[i64]]; 3] = [
            &[&[1, 1, 1]],
            &[&[1, 1], &[1, 1, 1]],
            &[&[1, 2], &[1, 1, 1]],
        ];
        for rows in inputs {
            let mut vectors: Vec<Vec<BigInt>> = vec![vec![BigInt::one()]];
            vectors.extend(
                rows.iter()
                    .map(|r| r.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>()),
            );
            let options = PfNetworkOptions {
                max_step_length: 4,
                max_rise: 4,
                ..PfNetworkOptions::default()
            };
            let mut search = SupportSearch::new(&vectors, &options);
            assert!(matches!(search.run(), SearchOutcome::Impossible(_)));
        }
    }

    #[test]
    fn eulerian_and_narayana_rows() {
        // A_1..A_3 and the Narayana rows N_1..N_3 fit Brenti's recursive
        // pattern (1,0), (1,1), (2,1); with the fourth rows no network has
        // steps of length <= 4 and rise <= 4 (exact, within these bounds).
        for (rows, name) in [
            (eulerian_rows(4), "eulerian"),
            (narayana_rows(4), "narayana"),
        ] {
            let result = recognize_pf_network(&rows[..3], &PfNetworkOptions::default()).unwrap();
            match &result {
                PfNetworkRecognition::Certified { certificate } => {
                    certificate
                        .verify_against_bigint_polynomials(&rows[..3])
                        .unwrap();
                }
                other => panic!("{}: {}", name, other),
            }
            let options = PfNetworkOptions {
                max_step_length: 4,
                max_rise: 4,
                ..PfNetworkOptions::default()
            };
            let result = recognize_pf_network(&rows, &options).unwrap();
            assert!(
                matches!(result, PfNetworkRecognition::ImpossibleWithinBounds { .. }),
                "{}: {}",
                name,
                result
            );
        }
    }

    fn narayana_rows(last: usize) -> Vec<Polynomial<BigInt>> {
        (1..=last)
            .map(|n| {
                Polynomial::new(
                    (0..=n)
                        .map(|k| {
                            binomial_bigint(n + 1, k + 1) * binomial_bigint(n + 1, k)
                                / BigInt::from(n + 1)
                        })
                        .collect(),
                )
            })
            .collect()
    }

    #[test]
    fn negative_coefficient_and_bad_options() {
        let result =
            recognize_pf_network(&polys(&[&[1, -1]]), &PfNetworkOptions::default()).unwrap();
        assert!(matches!(
            result,
            PfNetworkRecognition::ImpossibleWithinBounds { .. }
        ));
        assert_eq!(
            recognize_pf_network(&[], &PfNetworkOptions::default()),
            Err(PfNetworkRecognitionError::EmptyInput)
        );
        let options = PfNetworkOptions {
            first_row_index: 0,
            ..PfNetworkOptions::default()
        };
        assert!(matches!(
            recognize_pf_network(&polys(&[&[1, 1]]), &options),
            Err(PfNetworkRecognitionError::FirstRowNotOne { .. })
        ));
        let options = PfNetworkOptions {
            max_step_length: 0,
            ..PfNetworkOptions::default()
        };
        assert_eq!(
            recognize_pf_network(&polys(&[&[1, 1]]), &options),
            Err(PfNetworkRecognitionError::ZeroStepLength)
        );
    }

    #[test]
    fn search_budget_is_inconclusive() {
        let options = PfNetworkOptions {
            max_search_nodes: 2,
            ..PfNetworkOptions::default()
        };
        let result = recognize_pf_network(&pascal_rows(5), &options).unwrap();
        assert!(
            matches!(result, PfNetworkRecognition::Inconclusive { .. }),
            "{}",
            result
        );
    }

    #[test]
    fn search_finds_every_random_valid_network() {
        // Completeness check against the independent validator: build random
        // networks, keep the planar ones, and recognize their rows.
        let mut state: u64 = 0x9e37_79b9_7f4a_7c15;
        let mut next = |modulus: u64| -> u64 {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 33) % modulus
        };
        let types = primitive_step_types(3, 3);
        let (mut valid, mut tried) = (0, 0);
        while valid < 60 && tried < 20_000 {
            tried += 1;
            let last = 1 + next(4) as usize;
            let mut steps: Vec<Vec<PfStep>> = Vec::with_capacity(last);
            for n in 1..=last {
                let mut column = Vec::new();
                for t in types.iter().filter(|t| t.length <= n) {
                    if next(3) == 0 {
                        column.push(step(t.length, t.rise, 1 + next(3) as i64));
                    }
                }
                steps.push(column);
            }
            let network = HeightInvariantNetwork { steps };
            if network.validate().is_err() {
                continue;
            }
            valid += 1;
            let rows: Vec<Polynomial<BigInt>> = network
                .row_polynomials_prefix(last)
                .unwrap()
                .into_iter()
                .skip(1)
                .map(|p| Polynomial::new(p.coeffs().iter().map(|c| c.to_integer()).collect()))
                .collect();
            for c in network.row_polynomials_prefix(last).unwrap() {
                let coeffs: Vec<BigInt> = c.coeffs().iter().map(|v| v.to_integer()).collect();
                assert!(coeffs.is_empty() || is_real_rooted_bigint_coeffs(&coeffs));
            }
            for weights in [BrentiWeightDomain::Rational, BrentiWeightDomain::Integer] {
                let options = PfNetworkOptions {
                    weights,
                    ..PfNetworkOptions::default()
                };
                let result = recognize_pf_network(&rows, &options).unwrap();
                assert!(result.is_certified(), "{:?}: {}", network, result);
            }
        }
        assert!(
            valid >= 60,
            "only {} valid networks in {} tries",
            valid,
            tried
        );
    }
}
