//! Planar strip digraphs whose transitions depend only on the height.
//!
//! This module implements a restricted, proof-friendly class of weighted
//! planar digraphs:
//!
//! - vertices are `(x, y)` with `x >= 0` and `0 <= y <= H`,
//! - every edge goes from column `x` to column `x + 1`,
//! - the outgoing edges depend only on the height `y` (the digraph is
//!   invariant under translation in `x`),
//! - planarity is certified by an explicit monotonicity condition on target
//!   heights: if `y < y'` both have edges, every target of `y` is at most every
//!   target of `y'`.
//!
//! Write `W` for the `(H+1) x (H+1)` weight matrix and `M_{n,k}` for the total
//! weight of paths from `(0,0)` to `(n,k)`, so row `n` of `M` is `e_0 W^n`.
//!
//! # What a valid strip digraph certifies
//!
//! With sources `(-n, 0)` and sinks `(0, k)`, `x`-invariance makes `M_{n,k}`
//! the path weight from source `n` to sink `k`.  Two paths from sources
//! `n_1 < n_2` to sinks `k_1 < k_2` that pair `n_1` with `k_2` and `n_2` with
//! `k_1` must share a vertex (the monotone-target condition forbids crossing
//! edges), so by Lindström--Gessel--Viennot every minor of `M` is a sum of
//! nonnegative weights: **the coefficient matrix `M` is totally nonnegative**
//! (and so is `W`).  When the rows are the coefficient rows of a
//! lower-unitriangular matrix, the chain polynomials of
//! [`crate::tnn_network::check_chain_polynomials_bigint`] are then real-rooted.
//!
//! # What it does not certify
//!
//! A valid strip digraph does **not** make the row polynomials
//! `sum_k M_{n,k} t^k` real-rooted.  Heights `0 -> {0, 1}` and `1 -> {2}` with
//! unit weights form a valid planar strip digraph whose row `2` is
//! `1 + t + t^2`.  Brenti's theorem on Pólya frequency rows needs invariance
//! under translation in `y` (the height), which makes the Toeplitz matrix of a
//! row a path matrix; the height-dependent transitions used here are invariant
//! in `x` instead.  Use an exact real-rootedness test for the rows themselves,
//! or a PF certificate from [`crate::pf_network`] (height-invariant networks)
//! or [`crate::left_production`] (Chen--Fu--Ruan).
//!
//! [`crate::brenti_recognition`] searches for a strip digraph reproducing a
//! given finite prefix of rows.
//!
//! Entry points: build a [`BrentiStripDigraph`] (directly or with
//! [`BrentiStripDigraph::from_weight_matrix`]), call
//! [`BrentiStripDigraph::validate`], compute rows with
//! [`BrentiStripDigraph::row_polynomial`], and package a verified match with
//! [`build_brenti_sequence_certificate`].
//!
//! Conventions: heights and column indices are zero-based; row `n` is the
//! polynomial `sum_k M_{n,k} t^k` (ascending coefficients, `t^k` for height `k`)
//! with rational weights ([`crate::tnn_network::BigRational`]), trimmed as a
//! [`Polynomial`].  Everything is exact; failures are [`BrentiError`] values.

use crate::{tnn_network::BigRational, Polynomial};
use num_bigint::BigInt;
use num_traits::{One, Zero};
use std::fmt;

/// A single transition in a Brenti strip digraph.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BrentiEdge {
    /// Target height in the next column.
    pub to: usize,
    /// Nonnegative edge weight.
    pub weight: BigRational,
}

/// A planar, `x`-invariant digraph on a finite strip of heights.
///
/// From every vertex `(x, y)`, the outgoing edges are exactly `outgoing[y]`,
/// translated one column to the right.  See the module documentation for
/// what a valid strip digraph does and does not certify.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BrentiStripDigraph {
    /// `outgoing[y]` lists the edges leaving height `y`; the heights are
    /// `0..outgoing.len()`, so the height bound is `outgoing.len() - 1`.
    pub outgoing: Vec<Vec<BrentiEdge>>,
}

impl BrentiStripDigraph {
    /// Number of tracked heights.
    pub fn num_heights(&self) -> usize {
        self.outgoing.len()
    }

    /// Largest tracked height, if any.
    pub fn max_height(&self) -> Option<usize> {
        self.num_heights().checked_sub(1)
    }

    /// Build a strip digraph from a square weight matrix, where
    /// `weights[y][z]` is the weight of the edge from height `y` to height
    /// `z`.  Zero entries produce no edge.  The result is not validated.
    ///
    /// # Panics
    ///
    /// Panics if `weights` is not square.
    pub fn from_weight_matrix(weights: &[Vec<BigRational>]) -> Self {
        let size = weights.len();
        let outgoing = weights
            .iter()
            .map(|row| {
                assert_eq!(row.len(), size, "weight matrix must be square");
                row.iter()
                    .enumerate()
                    .filter(|(_, weight)| !weight.is_zero())
                    .map(|(to, weight)| BrentiEdge {
                        to,
                        weight: weight.clone(),
                    })
                    .collect()
            })
            .collect();
        Self { outgoing }
    }

    /// The square weight matrix `W`, where `W[y][z]` is the total weight of
    /// the edges from height `y` to height `z`.  Targets outside the strip are
    /// ignored.
    pub fn weight_matrix(&self) -> Vec<Vec<BigRational>> {
        let size = self.num_heights();
        let mut matrix = vec![vec![BigRational::zero(); size]; size];
        for (from, edges) in self.outgoing.iter().enumerate() {
            for edge in edges {
                if edge.to < size {
                    matrix[from][edge.to] += edge.weight.clone();
                }
            }
        }
        matrix
    }

    /// Validate nonnegativity, target range, and the monotone-target planarity condition.
    pub fn validate(&self) -> Result<(), BrentiError> {
        if self.outgoing.is_empty() {
            return Err(BrentiError::EmptyDigraph);
        }

        let max_height = self.num_heights() - 1;
        for (from, edges) in self.outgoing.iter().enumerate() {
            for edge in edges {
                if edge.to > max_height {
                    return Err(BrentiError::TargetOutOfRange {
                        from,
                        to: edge.to,
                        max_height,
                    });
                }
                if edge.weight < BigRational::zero() {
                    return Err(BrentiError::NegativeWeight {
                        from,
                        to: edge.to,
                        weight: edge.weight.clone(),
                    });
                }
            }
        }

        for lower in 0..self.num_heights() {
            let lower_max = self.outgoing[lower].iter().map(|edge| edge.to).max();
            let Some(lower_max) = lower_max else {
                continue;
            };

            for upper in (lower + 1)..self.num_heights() {
                let upper_min = self.outgoing[upper].iter().map(|edge| edge.to).min();
                let Some(upper_min) = upper_min else {
                    continue;
                };

                if lower_max > upper_min {
                    return Err(BrentiError::NonPlanar {
                        lower_height: lower,
                        upper_height: upper,
                        lower_max_target: lower_max,
                        upper_min_target: upper_min,
                    });
                }
            }
        }

        Ok(())
    }

    /// Weighted path counts from `(0,0)` to `(n,k)` for all tracked heights `k`.
    ///
    /// The result has one entry per height (`k = 0..num_heights`), untrimmed.
    /// The digraph is validated first.
    pub fn row_counts(&self, n: usize) -> Result<Vec<BigRational>, BrentiError> {
        self.validate()?;

        let mut current = vec![BigRational::zero(); self.num_heights()];
        current[0] = BigRational::one();

        for _ in 0..n {
            let mut next = vec![BigRational::zero(); self.num_heights()];
            for (from, weight_here) in current.iter().enumerate() {
                if weight_here.is_zero() {
                    continue;
                }
                for edge in &self.outgoing[from] {
                    next[edge.to] += weight_here.clone() * edge.weight.clone();
                }
            }
            current = next;
        }

        Ok(current)
    }

    /// The row polynomial `sum_k M_{n,k} t^k`.
    pub fn row_polynomial(&self, n: usize) -> Result<Polynomial<BigRational>, BrentiError> {
        Ok(Polynomial::new(self.row_counts(n)?))
    }

    /// A prefix of row polynomials for `n = 0, ..., last_row`.
    pub fn row_polynomials_prefix(
        &self,
        last_row: usize,
    ) -> Result<Vec<Polynomial<BigRational>>, BrentiError> {
        (0..=last_row).map(|n| self.row_polynomial(n)).collect()
    }
}

/// A verified match between a strip digraph and a polynomial sequence prefix.
///
/// It certifies that the listed rows are rows of the path matrix of a planar
/// network, hence that their coefficient matrix is totally nonnegative; it
/// says nothing about rows beyond the prefix and does not by itself make the
/// rows real-rooted (see the module documentation).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BrentiSequenceCertificate {
    /// The validated strip digraph.
    pub digraph: BrentiStripDigraph,
    /// `row_indices[i]` is the digraph row `n` (zero-based) that equals the
    /// `i`-th polynomial of the matched sequence.
    pub row_indices: Vec<usize>,
    /// The matched row polynomials `digraph.row_polynomial(row_indices[i])`.
    pub matched_rows: Vec<Polynomial<BigRational>>,
}

impl BrentiSequenceCertificate {
    /// Re-run the structural checks and row matching against `polys`, where
    /// `polys[i]` must equal row `row_indices[i]` of the digraph.
    pub fn verify_against_bigint_polynomials(
        &self,
        polys: &[Polynomial<BigInt>],
    ) -> Result<(), BrentiError> {
        self.digraph.validate()?;
        if polys.len() != self.row_indices.len() {
            return Err(BrentiError::RowIndexCountMismatch {
                num_polynomials: polys.len(),
                num_row_indices: self.row_indices.len(),
            });
        }

        for (i, (poly, &row_index)) in polys.iter().zip(self.row_indices.iter()).enumerate() {
            let found = self.digraph.row_polynomial(row_index)?;
            if found != bigint_poly_to_q(poly) {
                return Err(BrentiError::SequenceMismatch {
                    sequence_index: i,
                    row_index,
                    expected: bigint_poly_to_q(poly).coeffs().to_vec(),
                    found: found.coeffs().to_vec(),
                });
            }
        }
        Ok(())
    }
}

/// Errors returned by the restricted Brenti certificate machinery.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BrentiError {
    /// The digraph has no heights.
    EmptyDigraph,
    /// An edge targets a height above the strip.
    TargetOutOfRange {
        /// Source height of the edge.
        from: usize,
        /// Offending target height.
        to: usize,
        /// Largest valid height.
        max_height: usize,
    },
    /// An edge has a negative weight.
    NegativeWeight {
        /// Source height of the edge.
        from: usize,
        /// Target height of the edge.
        to: usize,
        /// The negative weight.
        weight: BigRational,
    },
    /// The monotone-target condition fails: a lower height reaches a target
    /// above one reached from a higher height.
    NonPlanar {
        /// The lower height `y`.
        lower_height: usize,
        /// The higher height `y' > y`.
        upper_height: usize,
        /// Largest target of `lower_height`.
        lower_max_target: usize,
        /// Smallest target of `upper_height` (smaller than `lower_max_target`).
        upper_min_target: usize,
    },
    /// The number of polynomials differs from the number of row indices.
    RowIndexCountMismatch {
        /// Number of supplied polynomials.
        num_polynomials: usize,
        /// Number of supplied row indices.
        num_row_indices: usize,
    },
    /// A digraph row differs from the corresponding polynomial.
    SequenceMismatch {
        /// Position `i` in the supplied polynomial list.
        sequence_index: usize,
        /// The digraph row index compared.
        row_index: usize,
        /// Coefficients of the supplied polynomial (ascending).
        expected: Vec<BigRational>,
        /// Coefficients of the digraph's row polynomial (ascending, trimmed).
        found: Vec<BigRational>,
    },
}

impl fmt::Display for BrentiError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyDigraph => write!(f, "digraph must contain at least one height"),
            Self::TargetOutOfRange {
                from,
                to,
                max_height,
            } => write!(
                f,
                "edge from height {} targets height {}, but the maximum supported height is {}",
                from, to, max_height
            ),
            Self::NegativeWeight { from, to, weight } => write!(
                f,
                "edge from height {} to height {} has negative weight {}",
                from, to, weight
            ),
            Self::NonPlanar {
                lower_height,
                upper_height,
                lower_max_target,
                upper_min_target,
            } => write!(
                f,
                "nonplanar strip certificate: height {} reaches as high as {}, but height {} reaches as low as {}",
                lower_height, lower_max_target, upper_height, upper_min_target
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
                "polynomial {} does not match Brenti row {}: expected {:?}, found {:?}",
                sequence_index, row_index, expected, found
            ),
        }
    }
}

impl std::error::Error for BrentiError {}

/// Build a verified strip-digraph certificate for a polynomial sequence prefix.
///
/// Checks [`BrentiStripDigraph::validate`] and that row `row_indices[i]` of
/// the digraph equals `polys[i]` exactly.
pub fn build_brenti_sequence_certificate(
    digraph: BrentiStripDigraph,
    row_indices: &[usize],
    polys: &[Polynomial<BigInt>],
) -> Result<BrentiSequenceCertificate, BrentiError> {
    digraph.validate()?;

    if polys.len() != row_indices.len() {
        return Err(BrentiError::RowIndexCountMismatch {
            num_polynomials: polys.len(),
            num_row_indices: row_indices.len(),
        });
    }

    let mut matched_rows = Vec::with_capacity(polys.len());
    for (sequence_index, (poly, &row_index)) in polys.iter().zip(row_indices.iter()).enumerate() {
        let found = digraph.row_polynomial(row_index)?;
        let expected = bigint_poly_to_q(poly);
        if found != expected {
            return Err(BrentiError::SequenceMismatch {
                sequence_index,
                row_index,
                expected: expected.coeffs().to_vec(),
                found: found.coeffs().to_vec(),
            });
        }
        matched_rows.push(found);
    }

    Ok(BrentiSequenceCertificate {
        digraph,
        row_indices: row_indices.to_vec(),
        matched_rows,
    })
}

fn bigint_poly_to_q(poly: &Polynomial<BigInt>) -> Polynomial<BigRational> {
    Polynomial::new(
        poly.coeffs()
            .iter()
            .map(|c| BigRational::from_integer(c.clone()))
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bi(v: i64) -> BigInt {
        BigInt::from(v)
    }

    fn br(v: i64) -> BigRational {
        BigRational::from_integer(bi(v))
    }

    fn bi_poly(coeffs: &[i64]) -> Polynomial<BigInt> {
        Polynomial::new(coeffs.iter().map(|&v| bi(v)).collect())
    }

    fn pascal_strip(max_height: usize) -> BrentiStripDigraph {
        let mut outgoing = Vec::with_capacity(max_height + 1);
        for y in 0..=max_height {
            let mut edges = vec![BrentiEdge {
                to: y,
                weight: br(1),
            }];
            if y < max_height {
                edges.push(BrentiEdge {
                    to: y + 1,
                    weight: br(1),
                });
            }
            outgoing.push(edges);
        }
        BrentiStripDigraph { outgoing }
    }

    #[test]
    fn test_pascal_strip_row_counts() {
        let digraph = pascal_strip(4);
        let row = digraph.row_counts(4).unwrap();
        assert_eq!(row, vec![br(1), br(4), br(6), br(4), br(1)]);
    }

    #[test]
    fn test_pascal_strip_row_polynomial() {
        let digraph = pascal_strip(4);
        let poly = digraph.row_polynomial(4).unwrap();
        assert_eq!(
            poly,
            Polynomial::new(vec![br(1), br(4), br(6), br(4), br(1)])
        );
    }

    #[test]
    fn test_nonplanar_strip_rejected() {
        let digraph = BrentiStripDigraph {
            outgoing: vec![
                vec![BrentiEdge {
                    to: 1,
                    weight: br(1),
                }],
                vec![BrentiEdge {
                    to: 0,
                    weight: br(1),
                }],
            ],
        };
        match digraph.validate() {
            Err(BrentiError::NonPlanar { .. }) => {}
            other => panic!("expected NonPlanar, got {:?}", other),
        }
    }

    #[test]
    fn test_build_brenti_certificate_pascal_rows() {
        let digraph = pascal_strip(4);
        let polys = vec![
            bi_poly(&[1, 1]),
            bi_poly(&[1, 2, 1]),
            bi_poly(&[1, 3, 3, 1]),
            bi_poly(&[1, 4, 6, 4, 1]),
        ];
        let cert = build_brenti_sequence_certificate(digraph, &[1, 2, 3, 4], &polys).unwrap();
        assert_eq!(cert.row_indices, vec![1, 2, 3, 4]);
        assert_eq!(cert.matched_rows.len(), 4);
    }

    #[test]
    fn test_valid_strip_digraph_can_have_non_real_rooted_rows() {
        // 0 -> {0, 1}, 1 -> {2}: planar and nonnegative, yet row 2 is 1 + t + t^2.
        let digraph = BrentiStripDigraph {
            outgoing: vec![
                vec![
                    BrentiEdge {
                        to: 0,
                        weight: br(1),
                    },
                    BrentiEdge {
                        to: 1,
                        weight: br(1),
                    },
                ],
                vec![BrentiEdge {
                    to: 2,
                    weight: br(1),
                }],
                vec![],
            ],
        };
        assert!(digraph.validate().is_ok());
        let row = digraph.row_counts(2).unwrap();
        assert_eq!(row, vec![br(1), br(1), br(1)]);
        assert!(!crate::is_real_rooted(&[1, 1, 1]));
    }

    #[test]
    fn test_weight_matrix_round_trip() {
        let digraph = pascal_strip(3);
        let matrix = digraph.weight_matrix();
        assert_eq!(BrentiStripDigraph::from_weight_matrix(&matrix), digraph);
    }

    #[test]
    fn test_build_brenti_certificate_mismatch() {
        let digraph = pascal_strip(3);
        let polys = vec![bi_poly(&[1, 2])];
        match build_brenti_sequence_certificate(digraph, &[1], &polys) {
            Err(BrentiError::SequenceMismatch { .. }) => {}
            other => panic!("expected SequenceMismatch, got {:?}", other),
        }
    }
}
