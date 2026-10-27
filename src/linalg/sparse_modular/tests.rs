use super::*;
use crate::linalg::modular::mod_inverse_prime;

// -----------------------------------------------------------------------
// Sparse modular linear systems
// -----------------------------------------------------------------------

fn dense_modular_linear_system_consistent(
    mut aug: Vec<Vec<u64>>,
    num_vars: usize,
    prime: u64,
) -> bool {
    let num_rows = aug.len();
    let mut pivot_row = 0;

    for col in 0..num_vars {
        let Some(pr) = (pivot_row..num_rows).find(|&row| aug[row][col] != 0) else {
            continue;
        };
        aug.swap(pivot_row, pr);
        let pivot_inv = mod_inverse_prime(aug[pivot_row][col], prime)
            .expect("nonzero element modulo prime is invertible");
        let pivot_snapshot = aug[pivot_row][col..=num_vars].to_vec();

        for aug_row in aug.iter_mut().take(num_rows).skip(pivot_row + 1) {
            if aug_row[col] == 0 {
                continue;
            }
            let factor = mul_mod_u64(aug_row[col], pivot_inv, prime);
            for (entry, pivot_entry) in aug_row[col..=num_vars]
                .iter_mut()
                .zip(pivot_snapshot.iter())
            {
                *entry = sub_mod_u64(*entry, mul_mod_u64(factor, *pivot_entry, prime), prime);
            }
        }

        pivot_row += 1;
        if pivot_row == num_rows {
            return true;
        }
    }

    aug[pivot_row..]
        .iter()
        .all(|row| row[..num_vars].iter().any(|&entry| entry != 0) || row[num_vars] == 0)
}

fn sparse_rows_from_dense_aug(aug: &[Vec<u64>], num_vars: usize, prime: u64) -> Vec<SparseModRow> {
    aug.iter()
        .map(|dense_row| {
            let mut row = SparseModRow::new(dense_row[num_vars], prime).unwrap();
            for (col, &value) in dense_row.iter().take(num_vars).enumerate() {
                row.add_entry(col, value, prime).unwrap();
            }
            row
        })
        .filter(|row| !row.is_tautology())
        .collect()
}

fn solution_satisfies_dense_aug(
    aug: &[Vec<u64>],
    num_vars: usize,
    prime: u64,
    solution: &[u64],
) -> bool {
    aug.iter().all(|row| {
        let lhs = row
            .iter()
            .take(num_vars)
            .zip(solution.iter())
            .fold(0, |acc, (&coeff, &value)| {
                add_mod_u64(acc, mul_mod_u64(coeff, value, prime), prime)
            });
        lhs == row[num_vars] % prime
    })
}

#[test]
fn test_sparse_mod_row_normalizes_and_cancels_entries() {
    let prime = 101;
    let mut row = SparseModRow::from_signed_rhs(-3, prime).unwrap();
    row.add_signed_entry(2, -2, prime).unwrap();
    row.add_entry(2, 103, prime).unwrap();
    row.add_signed_entry(1, -1, prime).unwrap();

    assert_eq!(row.rhs(), 98);
    assert_eq!(row.entries(), &[(1, 100)]);

    let row = SparseModRow::from_entries([(2, 5), (1, 3), (2, 96)], 204, prime).unwrap();
    assert_eq!(row.rhs(), 2);
    assert_eq!(row.entries(), &[(1, 3)]);
    assert!(!row.is_tautology());
    assert!(!row.is_contradiction());

    let row = SparseModRow::new(7, prime).unwrap();
    assert!(row.is_contradiction());
}

#[test]
fn test_sparse_mod_sorted_add_fast_path_and_fallback() {
    let prime = 101;
    let mut row = SparseModRow::new(0, prime).unwrap();

    row.add_reduced_entry_sorted(1, 4, prime);
    row.add_reduced_entry_sorted(3, 8, prime);
    row.add_reduced_entry_sorted(2, 5, prime);
    row.add_reduced_entry_sorted(3, 93, prime);

    assert_eq!(row.entries(), &[(1, 4), (2, 5)]);
}

#[test]
fn test_sparse_mod_small_pivot_reduction_cancels_in_place() {
    let prime = 101;
    let pivot = SparseModRow::from_entries([(0, 1), (2, 4)], 7, prime).unwrap();
    let mut row = SparseModRow::from_entries([(0, 5), (1, 3), (2, 20)], 9, prime).unwrap();

    row.subtract_scaled(&pivot, 5, prime);

    assert_eq!(row.rhs(), 75);
    assert_eq!(row.entries(), &[(1, 3)]);
}

#[test]
fn test_sparse_modular_solver_matches_dense_consistency() {
    let prime = 101;
    let num_vars = 3;
    let aug = vec![
        vec![1, 2, 0, 5],
        vec![0, 3, 4, 7],
        vec![2, 7, 4, 17],
        vec![0, 0, 0, 0],
    ];
    let sparse_rows = sparse_rows_from_dense_aug(&aug, num_vars, prime);
    let result = sparse_modular_linear_system_consistency(sparse_rows, num_vars, prime).unwrap();

    assert_eq!(
        result.consistent,
        dense_modular_linear_system_consistent(aug.clone(), num_vars, prime)
    );
    assert_eq!(result.rank, 2);
    assert_eq!(result.pivot_rows, vec![0, 1]);
    assert_eq!(result.pivot_columns, vec![0, 1]);
    assert_eq!(result.inconsistent_row, None);
    let solution = result.solution.as_ref().expect("consistent system");
    assert_eq!(solution.len(), num_vars);
    assert!(solution_satisfies_dense_aug(
        &aug, num_vars, prime, solution
    ));

    let direct_solution = sparse_modular_linear_system_solution(
        sparse_rows_from_dense_aug(&aug, num_vars, prime),
        num_vars,
        prime,
    )
    .unwrap()
    .expect("consistent system");
    assert!(solution_satisfies_dense_aug(
        &aug,
        num_vars,
        prime,
        &direct_solution
    ));
}

#[test]
fn test_sparse_modular_solver_can_process_sparsest_rows_first() {
    let prime = 101;
    let num_vars = 3;
    let aug = vec![vec![1, 1, 1, 1], vec![0, 1, 0, 2], vec![0, 0, 1, 3]];
    let result = sparse_modular_linear_system_consistency_with_options(
        sparse_rows_from_dense_aug(&aug, num_vars, prime),
        num_vars,
        prime,
        SparseModEliminationOptions {
            row_order: SparseModRowOrder::IncreasingNonzeros,
            compute_solution: true,
        },
    )
    .unwrap();

    assert!(result.consistent);
    assert_eq!(result.rank, 3);
    assert_eq!(result.pivot_rows, vec![1, 2, 0]);
    assert_eq!(result.pivot_columns, vec![1, 2, 0]);
    assert_eq!(result.stats.input_rows, 3);
    assert_eq!(result.stats.input_nonzeros, 5);
    assert_eq!(result.stats.row_reductions, 0);
    assert!(solution_satisfies_dense_aug(
        &aug,
        num_vars,
        prime,
        result.solution.as_ref().expect("consistent system")
    ));

    let result = sparse_modular_linear_system_consistency_with_options(
        sparse_rows_from_dense_aug(&aug, num_vars, prime),
        num_vars,
        prime,
        SparseModEliminationOptions {
            row_order: SparseModRowOrder::IncreasingNonzeros,
            compute_solution: false,
        },
    )
    .unwrap();
    assert!(result.consistent);
    assert_eq!(result.solution, None);
}

#[test]
fn test_sparse_modular_solver_checks_tail_rows_after_full_rank() {
    let prime = 101;
    let num_vars = 2;
    let aug = vec![vec![1, 0, 2], vec![0, 1, 3], vec![1, 1, 5], vec![2, 3, 13]];
    let result = sparse_modular_linear_system_consistency_with_options(
        sparse_rows_from_dense_aug(&aug, num_vars, prime),
        num_vars,
        prime,
        SparseModEliminationOptions {
            row_order: SparseModRowOrder::Input,
            compute_solution: false,
        },
    )
    .unwrap();

    assert!(result.consistent);
    assert_eq!(result.rank, 2);
    assert_eq!(result.pivot_rows, vec![0, 1]);
    assert_eq!(result.stats.row_reductions, 0);
    assert_eq!(result.stats.full_rank_checks, 2);
    assert_eq!(result.solution, None);

    let inconsistent_aug = vec![vec![1, 0, 2], vec![0, 1, 3], vec![1, 1, 6]];
    let result = sparse_modular_linear_system_consistency_with_options(
        sparse_rows_from_dense_aug(&inconsistent_aug, num_vars, prime),
        num_vars,
        prime,
        SparseModEliminationOptions {
            row_order: SparseModRowOrder::Input,
            compute_solution: false,
        },
    )
    .unwrap();

    assert!(!result.consistent);
    assert_eq!(result.rank, 2);
    assert_eq!(result.inconsistent_row, Some(2));
    assert_eq!(result.stats.full_rank_checks, 1);
}

#[test]
fn test_sparse_modular_solver_full_rank_solution_mode() {
    let prime = 101;
    let num_vars = 3;
    let rank_deficient_aug = vec![vec![1, 0, 0, 2], vec![0, 1, 0, 3]];
    let full_rank_only =
        sparse_modular_linear_system_consistency_prime_unchecked_with_solution_mode(
            sparse_rows_from_dense_aug(&rank_deficient_aug, num_vars, prime),
            num_vars,
            prime,
            SparseModRowOrder::IncreasingNonzeros,
            SparseModSolutionMode::FullRank,
        )
        .unwrap();
    assert!(full_rank_only.consistent);
    assert_eq!(full_rank_only.rank, 2);
    assert_eq!(full_rank_only.solution, None);

    let any_consistent =
        sparse_modular_linear_system_consistency_prime_unchecked_with_solution_mode(
            sparse_rows_from_dense_aug(&rank_deficient_aug, num_vars, prime),
            num_vars,
            prime,
            SparseModRowOrder::IncreasingNonzeros,
            SparseModSolutionMode::AnyConsistent,
        )
        .unwrap();
    assert_eq!(any_consistent.rank, 2);
    assert!(solution_satisfies_dense_aug(
        &rank_deficient_aug,
        num_vars,
        prime,
        any_consistent
            .solution
            .as_ref()
            .expect("rank-deficient solution")
    ));

    let full_rank_aug = vec![vec![1, 0, 0, 2], vec![0, 1, 0, 3], vec![0, 0, 1, 4]];
    let full_rank = sparse_modular_linear_system_consistency_prime_unchecked_with_solution_mode(
        sparse_rows_from_dense_aug(&full_rank_aug, num_vars, prime),
        num_vars,
        prime,
        SparseModRowOrder::IncreasingNonzeros,
        SparseModSolutionMode::FullRank,
    )
    .unwrap();
    assert_eq!(full_rank.rank, 3);
    assert!(solution_satisfies_dense_aug(
        &full_rank_aug,
        num_vars,
        prime,
        full_rank.solution.as_ref().expect("full-rank solution")
    ));
}

#[test]
fn test_sparse_modular_solver_can_use_initial_markowitz_cost() {
    let prime = 101;
    let num_vars = 5;
    let rows = vec![
        SparseModRow::from_entries([(0, 1), (3, 1)], 1, prime).unwrap(),
        SparseModRow::from_entries([(0, 1), (4, 1)], 2, prime).unwrap(),
        SparseModRow::from_entries([(1, 1), (3, 1)], 3, prime).unwrap(),
        SparseModRow::from_entries([(1, 1), (4, 1)], 4, prime).unwrap(),
        SparseModRow::from_entries([(2, 1), (3, 1)], 5, prime).unwrap(),
    ];
    let result = sparse_modular_linear_system_consistency_with_options(
        rows,
        num_vars,
        prime,
        SparseModEliminationOptions {
            row_order: SparseModRowOrder::IncreasingMarkowitzCost,
            compute_solution: false,
        },
    )
    .unwrap();

    assert!(result.consistent);
    assert_eq!(result.pivot_rows[0], 4);
    assert_eq!(result.pivot_columns[0], 2);
    assert_eq!(result.solution, None);
}

#[test]
fn test_sparse_modular_solver_dynamic_markowitz_matches_dense() {
    let prime = 101;
    let num_vars = 4;
    let aug = vec![
        vec![1, 1, 1, 1, 10],
        vec![0, 0, 1, 1, 7],
        vec![0, 1, 0, 1, 8],
        vec![1, 0, 0, 0, 2],
        vec![2, 1, 1, 2, 19],
        vec![0, 0, 0, 0, 0],
    ];
    let result = sparse_modular_linear_system_consistency_with_options(
        sparse_rows_from_dense_aug(&aug, num_vars, prime),
        num_vars,
        prime,
        SparseModEliminationOptions {
            row_order: SparseModRowOrder::DynamicMarkowitz,
            compute_solution: true,
        },
    )
    .unwrap();

    assert_eq!(
        result.consistent,
        dense_modular_linear_system_consistent(aug.clone(), num_vars, prime)
    );
    assert_eq!(result.rank, 4);
    assert_eq!(result.solution.as_ref().map(Vec::len), Some(num_vars));
    assert!(solution_satisfies_dense_aug(
        &aug,
        num_vars,
        prime,
        result.solution.as_ref().expect("consistent system")
    ));
    assert_eq!(result.stats.zero_rows, 0);
    assert!(result.stats.row_reductions > 0);
}

#[test]
fn test_sparse_modular_solver_dynamic_markowitz_detects_inconsistency() {
    let prime = 101;
    let num_vars = 3;
    let aug = vec![vec![1, 0, 0, 2], vec![0, 1, 0, 3], vec![1, 1, 0, 6]];
    let result = sparse_modular_linear_system_consistency_with_options(
        sparse_rows_from_dense_aug(&aug, num_vars, prime),
        num_vars,
        prime,
        SparseModEliminationOptions {
            row_order: SparseModRowOrder::DynamicMarkowitz,
            compute_solution: true,
        },
    )
    .unwrap();

    assert!(!result.consistent);
    assert_eq!(result.rank, 2);
    assert_eq!(result.inconsistent_row, Some(2));
    assert_eq!(result.solution, None);
}

#[test]
fn test_sparse_modular_solver_detects_inconsistency() {
    let prime = 101;
    let num_vars = 2;
    let aug = vec![vec![1, 1, 1], vec![2, 2, 3]];
    let sparse_rows = sparse_rows_from_dense_aug(&aug, num_vars, prime);
    let result = sparse_modular_linear_system_consistency(sparse_rows, num_vars, prime).unwrap();

    assert!(!result.consistent);
    assert_eq!(result.rank, 1);
    assert_eq!(result.pivot_rows, vec![0]);
    assert_eq!(result.pivot_columns, vec![0]);
    assert_eq!(result.inconsistent_row, Some(1));
    assert_eq!(result.solution, None);
}

#[test]
fn test_sparse_modular_solver_rejects_bad_inputs() {
    let row = {
        let mut row = SparseModRow::new(0, 101).unwrap();
        row.add_entry(3, 1, 101).unwrap();
        row
    };
    assert_eq!(
        sparse_modular_linear_system_consistency(vec![row], 3, 101),
        Err(SparseModEliminationError::ColumnOutOfRange {
            row: 0,
            column: 3,
            num_vars: 3
        })
    );
    assert_eq!(
        sparse_modular_linear_system_consistency(Vec::new(), 3, 100),
        Err(SparseModEliminationError::ModulusNotPrime { modulus: 100 })
    );
    assert_eq!(SparseModRow::new(1, 0), Err(SparseModRowError::ZeroModulus));
    assert_eq!(
        SparseModRow::from_signed_rhs(-1, 0),
        Err(SparseModRowError::ZeroModulus)
    );
    let mut row = SparseModRow::new(0, 101).unwrap();
    assert_eq!(row.add_entry(0, 1, 0), Err(SparseModRowError::ZeroModulus));
}

#[test]
fn test_sparse_mod_row_from_signed_entries_matches_unsigned_construction() {
    // Oracle: reducing signed entries by hand and using the unsigned constructor.
    let prime = 101;
    let signed = SparseModRow::from_signed_entries(
        [(3, -2), (1, 5), (3, 2), (0, -101), (4, 203)],
        -4,
        prime,
    )
    .unwrap();
    // Column 3 cancels, column 0 is a multiple of the prime, column 4 reduces to 1.
    let unsigned = SparseModRow::from_entries([(1, 5), (4, 1)], 97, prime).unwrap();
    assert_eq!(signed, unsigned);
    assert_eq!(signed.entries(), &[(1, 5), (4, 1)]);
    assert_eq!(signed.rhs(), 97);
    assert!(SparseModRow::from_signed_entries([(0, 1)], 0, 0).is_err());
}

#[test]
fn test_sparse_mod_row_is_satisfied_by_agrees_with_direct_evaluation() {
    // Row: 2 x_0 + 3 x_2 = 5 (mod 7).
    let prime = 7;
    let row = SparseModRow::from_entries([(0, 2), (2, 3)], 5, prime).unwrap();
    for x0 in 0..prime {
        for x2 in 0..prime {
            let solution = [x0, 4, x2];
            let expected = (2 * x0 + 3 * x2) % prime == 5;
            assert_eq!(row.is_satisfied_by(&solution, prime).unwrap(), expected);
        }
    }
    // A solution vector missing a referenced column never satisfies the row.
    assert!(!row.is_satisfied_by(&[1, 1], prime).unwrap());
    assert!(row.is_satisfied_by(&[0, 0, 4], 0).is_err());
}

#[test]
fn test_is_sparse_modular_linear_system_consistent_matches_dense_on_pseudorandom_systems() {
    // Deterministic linear congruential generator; the oracle is the dense
    // elimination `dense_modular_linear_system_consistent` above.
    let mut state: u64 = 0x2545_f491_4f6c_dd1d;
    let mut next = |modulus: u64| {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) % modulus
    };
    let prime = 7;
    let mut seen_consistent = 0;
    let mut seen_inconsistent = 0;
    for _ in 0..200 {
        let num_vars = 2 + next(3) as usize;
        let num_rows = 1 + next(5) as usize;
        let aug: Vec<Vec<u64>> = (0..num_rows)
            .map(|_| {
                (0..=num_vars)
                    .map(|col| {
                        if col < num_vars && next(2) == 0 {
                            0
                        } else {
                            next(prime)
                        }
                    })
                    .collect()
            })
            .collect();
        let expected = dense_modular_linear_system_consistent(aug.clone(), num_vars, prime);
        let actual = is_sparse_modular_linear_system_consistent(
            sparse_rows_from_dense_aug(&aug, num_vars, prime),
            num_vars,
            prime,
        )
        .unwrap();
        assert_eq!(actual, expected, "{aug:?}");
        if expected {
            seen_consistent += 1;
        } else {
            seen_inconsistent += 1;
        }
    }
    assert!(seen_consistent > 10 && seen_inconsistent > 10);
}
