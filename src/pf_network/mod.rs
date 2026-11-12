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

mod recognition;
pub use recognition::{primitive_step_types, recognize_pf_network};
#[cfg(test)]
mod tests;
