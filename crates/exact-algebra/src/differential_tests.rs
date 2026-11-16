//! Seeded differential tests of the exact linear solvers against an independent, naive
//! rational Gauss-Jordan elimination with normalised pivots.
//!
//! The reference below shares no code with the crate's own eliminations (neither the
//! integer elimination nor `rational_rref`), so a systematic error in either would show.

use super::integer_elimination::{
    integer_echelon_for_solve, integer_echelon_from_integer_rows, ForcedFallback,
    IntegerEliminationResult,
};
use super::linear_systems::solve_integer_fit_system;
use super::*;
use num_bigint::BigInt;
use num_traits::{One, Zero};

/// A small deterministic linear congruential generator.
struct Rng(u64);

impl Rng {
    fn next(&mut self, bound: u64) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 33) % bound.max(1)
    }

    /// Uniform in `-radius..=radius`.
    fn small(&mut self, radius: i64) -> i64 {
        self.next((2 * radius + 1) as u64) as i64 - radius
    }
}

/// Reduced row echelon form with normalised pivots (first nonzero entry as pivot, row
/// swaps), over the first `columns` columns. Returns the rows, the pivot columns and
/// the input row index of each pivot row.
fn naive_rref(matrix: &[Vec<Q>], columns: usize) -> (Vec<Vec<Q>>, Vec<usize>, Vec<usize>) {
    let mut rows = matrix.to_vec();
    let mut origins: Vec<usize> = (0..rows.len()).collect();
    let mut pivots = Vec::new();
    for column in 0..columns {
        let target = pivots.len();
        let Some(found) = (target..rows.len()).find(|&i| !rows[i][column].is_zero()) else {
            continue;
        };
        rows.swap(target, found);
        origins.swap(target, found);
        let inverse = Q::one() / &rows[target][column];
        for entry in rows[target].iter_mut() {
            *entry = &*entry * &inverse;
        }
        let pivot_row = rows[target].clone();
        for (i, row) in rows.iter_mut().enumerate() {
            if i != target && !row[column].is_zero() {
                let factor = row[column].clone();
                for (entry, pivot_entry) in row.iter_mut().zip(&pivot_row) {
                    *entry -= &factor * pivot_entry;
                }
            }
        }
        pivots.push(column);
    }
    origins.truncate(pivots.len());
    (rows, pivots, origins)
}

fn augmented(a: &[Vec<Q>], b: &[Q]) -> Vec<Vec<Q>> {
    a.iter()
        .zip(b)
        .map(|(row, rhs)| row.iter().cloned().chain([rhs.clone()]).collect())
        .collect()
}

/// The reference solution (free variables zero) with its rank and nullity, or `None`
/// for an inconsistent system.
fn naive_solve(a: &[Vec<Q>], b: &[Q], columns: usize) -> Option<LinearSystemSolution> {
    let (reduced, pivots, _) = naive_rref(&augmented(a, b), columns);
    if reduced[pivots.len()..]
        .iter()
        .any(|row| !row[columns].is_zero())
    {
        return None;
    }
    let mut particular = vec![Q::zero(); columns];
    for (i, &pivot) in pivots.iter().enumerate() {
        particular[pivot] = reduced[i][columns].clone();
    }
    Some(LinearSystemSolution {
        particular,
        rank: pivots.len(),
        nullity: columns - pivots.len(),
    })
}

fn mat_vec(a: &[Vec<Q>], x: &[Q]) -> Vec<Q> {
    a.iter()
        .map(|row| row.iter().zip(x).map(|(u, v)| u * v).sum())
        .collect()
}

/// A random system: sparse small entries, some zero rows and columns, duplicated and
/// scaled rows, combinations of earlier rows, optionally fractional entries (also with
/// negative denominators) and entries scaled by `2^120`; the right-hand side is either
/// consistent or arbitrary.
fn random_system(
    rng: &mut Rng,
    rows: usize,
    cols: usize,
    fractional: bool,
    big: bool,
) -> (Vec<Vec<Q>>, Vec<Q>) {
    let zero_column = (cols > 0 && rng.next(4) == 0).then(|| rng.next(cols as u64) as usize);
    let scale = if big {
        Q::from_integer(BigInt::one() << 120)
    } else {
        Q::one()
    };
    let mut a: Vec<Vec<Q>> = Vec::new();
    for r in 0..rows {
        let kind = rng.next(6);
        let row: Vec<Q> = if r >= 1 && kind == 0 {
            let source = rng.next(r as u64) as usize;
            let factor = Q::new(
                BigInt::from(rng.small(3) | 1),
                BigInt::from(1 + rng.next(3) as i64),
            );
            a[source].iter().map(|x| x * &factor).collect()
        } else if r >= 2 && kind == 1 {
            let (s, t) = (rng.next(r as u64) as usize, rng.next(r as u64) as usize);
            let (x, y) = (
                Q::from_integer(rng.small(3).into()),
                Q::from_integer(rng.small(3).into()),
            );
            (0..cols).map(|c| &a[s][c] * &x + &a[t][c] * &y).collect()
        } else if kind == 2 && rng.next(3) == 0 {
            vec![Q::zero(); cols]
        } else {
            (0..cols)
                .map(|c| {
                    if Some(c) == zero_column || rng.next(5) < 2 {
                        Q::zero()
                    } else {
                        let numerator = Q::from_integer(rng.small(9).into()) * &scale;
                        if fractional && rng.next(3) == 0 {
                            let denominator = match rng.small(7) {
                                0 => 5,
                                d => d,
                            };
                            // A negative denominator exercises `Ratio` normalisation.
                            numerator / Q::from_integer(denominator.into())
                        } else {
                            numerator
                        }
                    }
                })
                .collect()
        };
        a.push(row);
    }
    let b: Vec<Q> = if rng.next(2) == 0 {
        let x0: Vec<Q> = (0..cols)
            .map(|_| Q::new(rng.small(5).into(), BigInt::from(1 + rng.next(3) as i64)))
            .collect();
        mat_vec(&a, &x0)
    } else {
        (0..rows)
            .map(|_| {
                if rng.next(3) == 0 {
                    Q::zero()
                } else {
                    Q::from_integer(rng.small(7).into()) * &scale
                }
            })
            .collect()
    };
    (a, b)
}

#[derive(Default)]
struct Counts {
    consistent: usize,
    inconsistent: usize,
    square_unique: usize,
    square_singular: usize,
    square_singular_inconsistent: usize,
}

/// Compare every solver with the reference on one system.
fn check_system(a: &[Vec<Q>], b: &[Q], cols: usize, tag: &str, counts: &mut Counts) {
    let rows = a.len();
    let expected = if rows == 0 {
        Some(LinearSystemSolution {
            particular: vec![],
            rank: 0,
            nullity: 0,
        })
    } else {
        naive_solve(a, b, cols)
    };

    let got = solve_linear_system_with_rank(a, b);
    assert_eq!(got, expected, "general {tag}");
    match &got {
        Some(solution) => {
            counts.consistent += 1;
            if rows > 0 {
                assert_eq!(mat_vec(a, &solution.particular), b, "residual {tag}");
            }
        }
        None => counts.inconsistent += 1,
    }
    assert_eq!(
        solve_linear_system(a, b),
        expected.as_ref().map(|s| s.particular.clone()),
        "solve_linear_system {tag}"
    );

    // The early-exit forward elimination may only report a true contradiction.
    if rows > 0 {
        if let Some(IntegerEliminationResult::Inconsistent) =
            integer_echelon_for_solve(&augmented(a, b), cols)
        {
            assert!(expected.is_none(), "false contradiction {tag}");
        }
    }

    // The integer fit path, on integral systems.
    let integral = a.iter().flatten().chain(b).all(|x| x.denom().is_one());
    if integral && rows > 0 {
        let integer_a: Vec<Vec<BigInt>> = a
            .iter()
            .map(|row| row.iter().map(|x| x.numer().clone()).collect())
            .collect();
        let integer_b: Vec<BigInt> = b.iter().map(|x| x.numer().clone()).collect();
        assert_eq!(
            solve_integer_fit_system(integer_a.clone(), integer_b.clone(), false),
            expected.as_ref().map(|s| s.particular.clone()),
            "fit {tag}"
        );
        let mut integer_aug = integer_a.clone();
        for (row, x) in integer_aug.iter_mut().zip(&integer_b) {
            row.push(x.clone());
        }
        if let Some(IntegerEliminationResult::Inconsistent) =
            integer_echelon_from_integer_rows(integer_aug, cols)
        {
            assert!(expected.is_none(), "false contradiction (integer) {tag}");
        }
        if rows == cols {
            let unique = expected
                .as_ref()
                .filter(|s| s.rank == cols)
                .map(|s| s.particular.clone());
            assert_eq!(
                solve_integer_fit_system(integer_a, integer_b, true),
                unique,
                "fit square {tag}"
            );
        }
    }

    // The square solver: a unique solution, or `None` if singular (consistent or not).
    if rows == cols {
        let unique = if rows == 0 {
            Some(vec![])
        } else {
            expected
                .as_ref()
                .filter(|s| s.rank == cols)
                .map(|s| s.particular.clone())
        };
        let got = solve_full_rank_square_linear_system(a, b);
        if got.is_some() {
            counts.square_unique += 1;
        } else {
            counts.square_singular += 1;
            if expected.is_none() {
                counts.square_singular_inconsistent += 1;
            }
        }
        assert_eq!(got, unique, "square {tag}");
    }

    // The kernel with its rank certificate.
    let limits = RationalKernelLimits {
        max_entries: 1 << 20,
        max_operations: usize::MAX,
        max_coefficient_bits: u64::MAX,
    };
    let kernel = rational_kernel_with_certificate(a, cols, limits).expect("kernel");
    let (reduced, pivots, origins) = naive_rref(a, cols);
    assert_eq!(kernel.pivot_columns, pivots, "kernel pivots {tag}");
    assert_eq!(kernel.independent_rows, origins, "kernel rows {tag}");
    let free: Vec<usize> = (0..cols).filter(|c| !pivots.contains(c)).collect();
    assert_eq!(kernel.free_columns, free, "kernel free {tag}");
    assert_eq!(kernel.basis.len(), free.len(), "kernel dimension {tag}");
    for (vector, &f) in kernel.basis.iter().zip(&free) {
        let mut expected_vector = vec![Q::zero(); cols];
        expected_vector[f] = Q::one();
        for (i, &p) in pivots.iter().enumerate() {
            expected_vector[p] = -reduced[i][f].clone();
        }
        assert_eq!(vector, &expected_vector, "kernel basis {tag}");
    }
    let rank = pivots.len();
    for i in 0..rank {
        for j in 0..rank {
            let entry: Q = (0..rank)
                .map(|t| &kernel.inverse_minor[i][t] * &a[origins[t]][pivots[j]])
                .sum();
            let identity = if i == j { Q::one() } else { Q::zero() };
            assert_eq!(entry, identity, "inverse minor {tag}");
        }
    }
}

fn run_random_cases(seed: u64, cases: u64, counts: &mut Counts) {
    let mut rng = Rng(seed);
    for case in 0..cases {
        let rows = rng.next(8) as usize;
        let cols = rng.next(8) as usize;
        let fractional = rng.next(2) == 0;
        let big = rng.next(10) == 0;
        let (a, b) = random_system(&mut rng, rows, cols, fractional, big);
        let tag = format!("case {case} (seed {seed:#x}): a={a:?} b={b:?}");
        check_system(&a, &b, cols, &tag, counts);
    }
}

fn assert_coverage(counts: &Counts) {
    assert!(
        counts.consistent > 200
            && counts.inconsistent > 200
            && counts.square_unique > 50
            && counts.square_singular > 50
            && counts.square_singular_inconsistent > 20,
        "weak coverage: consistent {}, inconsistent {}, square unique {}, square singular {} \
         (inconsistent {})",
        counts.consistent,
        counts.inconsistent,
        counts.square_unique,
        counts.square_singular,
        counts.square_singular_inconsistent
    );
}

#[test]
fn solvers_match_the_naive_rational_elimination_on_random_small_systems() {
    let mut counts = Counts::default();
    run_random_cases(0x5eed_1234, 6000, &mut counts);
    assert_coverage(&counts);
}

/// The rational fallbacks (reached only if an integer division were inexact) must give
/// the same answers as the integer routes.
#[test]
fn forced_integer_fallback_matches_the_naive_rational_elimination() {
    let _forced = ForcedFallback::engage();
    let mut counts = Counts::default();
    run_random_cases(0x0fa1_1bac, 2500, &mut counts);
    assert_coverage(&counts);
}

#[test]
fn square_inconsistent_singular_systems_have_no_solution_on_every_route() {
    let q = |v: i64| Q::from_integer(BigInt::from(v));
    let singular = [
        vec![vec![q(1), q(2)], vec![q(2), q(4)]],
        vec![vec![q(0), q(0)], vec![q(0), q(1)]],
        vec![
            vec![q(1), q(1), q(1)],
            vec![q(2), q(2), q(2)],
            vec![q(0), q(1), q(5)],
        ],
        vec![vec![q(0)]],
    ];
    for a in &singular {
        let n = a.len();
        // A right-hand side outside the column space: the last equation is violated.
        let mut inconsistent: Vec<Q> = (0..n).map(|i| q(i as i64 + 1)).collect();
        inconsistent[n - 1] += q(7);
        for b in [inconsistent, vec![q(0); n]] {
            for forced in [false, true] {
                let _guard = forced.then(ForcedFallback::engage);
                let expected = naive_solve(a, &b, n);
                assert!(expected.as_ref().is_none_or(|s| s.rank < n));
                assert_eq!(solve_linear_system_with_rank(a, &b), expected);
                assert_eq!(solve_full_rank_square_linear_system(a, &b), None);
                let integer_a: Vec<Vec<BigInt>> = a
                    .iter()
                    .map(|row| row.iter().map(|x| x.numer().clone()).collect())
                    .collect();
                let integer_b: Vec<BigInt> = b.iter().map(|x| x.numer().clone()).collect();
                assert_eq!(
                    solve_integer_fit_system(integer_a.clone(), integer_b.clone(), true),
                    None
                );
                assert_eq!(
                    solve_integer_fit_system(integer_a, integer_b, false),
                    expected.map(|s| s.particular)
                );
            }
        }
    }
}

/// Larger structured systems: many skipped columns (delayed row levels), entries near
/// `2^31` and some shifted by `2^90`, with a unit perturbation making a third of them
/// inconsistent.
#[test]
fn solvers_match_the_naive_rational_elimination_on_sparse_structured_systems() {
    let mut rng = Rng(42);
    for case in 0..120u64 {
        let n = 1 + rng.next(10) as usize;
        let m = n + rng.next(4) as usize;
        let mut a = vec![vec![Q::zero(); n]; m];
        for row in a.iter_mut() {
            for entry in row.iter_mut() {
                if rng.next(4) == 0 {
                    let shift = if rng.next(5) == 0 { 90 } else { 0 };
                    *entry = Q::from_integer(BigInt::from(rng.small(1_000_000_007)) << shift);
                }
            }
        }
        let x0: Vec<Q> = (0..n)
            .map(|_| Q::from_integer(rng.small(50).into()))
            .collect();
        let mut b = mat_vec(&a, &x0);
        if rng.next(3) == 0 {
            let i = rng.next(m as u64) as usize;
            b[i] += Q::one();
        }
        let tag = format!("sparse case {case}");
        let expected = naive_solve(&a, &b, n);
        assert_eq!(solve_linear_system_with_rank(&a, &b), expected, "{tag}");
        let integer_a: Vec<Vec<BigInt>> = a
            .iter()
            .map(|row| row.iter().map(|x| x.numer().clone()).collect())
            .collect();
        let integer_b: Vec<BigInt> = b.iter().map(|x| x.numer().clone()).collect();
        assert_eq!(
            solve_integer_fit_system(integer_a.clone(), integer_b.clone(), false),
            expected.as_ref().map(|s| s.particular.clone()),
            "fit {tag}"
        );
        let unique = naive_solve(&a[..n], &b[..n], n)
            .filter(|s| s.rank == n)
            .map(|s| s.particular);
        assert_eq!(
            solve_full_rank_square_linear_system(&a[..n], &b[..n]),
            unique,
            "square {tag}"
        );
        assert_eq!(
            solve_integer_fit_system(integer_a[..n].to_vec(), integer_b[..n].to_vec(), true),
            unique,
            "fit square {tag}"
        );
    }
}
