//! Differential checks against the former scalar and nullspace kernels.
use super::*;
fn original_inverse_scalar(value: u64, prime: u64) -> u64 {
    assert!(value != 0 && value < prime);
    let mut old_r = value as i128;
    let mut r = prime as i128;
    let mut old_s = 1i128;
    let mut s = 0i128;
    while r != 0 {
        let quotient = old_r / r;
        (old_r, r) = (r, old_r - quotient * r);
        (old_s, s) = (s, old_s - quotient * s);
    }
    (old_s.rem_euclid(prime as i128)) as u64
}

fn original_nullspace(mut matrix: Vec<Vec<u64>>, prime: u64) -> Vec<Fp> {
    let rows = matrix.len();
    let columns = matrix.first().map_or(0, Vec::len);
    let mut pivot_columns = Vec::new();
    let mut pivot_row = 0;
    for column in 0..columns {
        let Some(source_row) = (pivot_row..rows).find(|row| matrix[*row][column] != 0) else {
            continue;
        };
        matrix.swap(pivot_row, source_row);
        let inverse = original_inverse_scalar(matrix[pivot_row][column], prime);
        for value in &mut matrix[pivot_row] {
            *value = *value * inverse % prime;
        }
        for row in 0..rows {
            if row == pivot_row || matrix[row][column] == 0 {
                continue;
            }
            let multiple = matrix[row][column];
            let pivot_tail = matrix[pivot_row][column..].to_vec();
            for (offset, value) in matrix[row][column..].iter_mut().enumerate() {
                *value = (*value + prime - multiple * pivot_tail[offset] % prime) % prime;
            }
        }
        pivot_columns.push(column);
        pivot_row += 1;
        if pivot_row == rows {
            break;
        }
    }
    let is_pivot = (0..columns)
        .map(|column| pivot_columns.contains(&column))
        .collect::<Vec<_>>();
    let mut result = Vec::new();
    for free_column in 0..columns {
        if is_pivot[free_column] {
            continue;
        }
        let mut vector = vec![0; columns];
        vector[free_column] = 1;
        for (row, pivot_column) in pivot_columns.iter().copied().enumerate() {
            vector[pivot_column] = (prime - matrix[row][free_column]) % prime;
        }
        result.push(trim(vector));
    }
    result
}

fn original_is_prime(candidate: u64) -> bool {
    if candidate < 2 {
        return false;
    }
    if candidate == 2 {
        return true;
    }
    if candidate % 2 == 0 {
        return false;
    }
    let mut divisor = 3u64;
    while divisor <= candidate / divisor {
        if candidate % divisor == 0 {
            return false;
        }
        divisor += 2;
    }
    true
}
#[test]
fn former_finite_field_kernels_agree_on_random_and_adversarial_inputs() {
    let mut seed = 0x8a39_24b0_12cf_907du64;
    for prime in [2, 3, 5, 7, 101, 1_000_000_007] {
        for rows in 0..=8 {
            for columns in 0..=8 {
                let matrix: Vec<Vec<u64>> = (0..rows)
                    .map(|_| {
                        (0..columns)
                            .map(|_| {
                                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                                (seed >> 32) % prime
                            })
                            .collect()
                    })
                    .collect();
                assert_eq!(
                    nullspace(matrix.clone(), prime),
                    original_nullspace(matrix, prime)
                );
            }
        }
        for value in 1..prime.min(200) {
            assert_eq!(
                inverse_scalar(value, prime),
                original_inverse_scalar(value, prime)
            );
        }
        let vanished = vec![vec![0, 1], vec![0, 1]];
        assert_eq!(
            nullspace(vanished.clone(), prime),
            original_nullspace(vanished, prime)
        );
    }
    for candidate in 0..20_000 {
        assert_eq!(is_prime(candidate), original_is_prime(candidate));
    }
    assert_eq!(
        nullspace(vec![vec![0, 0]; 3], 101),
        vec![vec![1], vec![0, 1]]
    );
}
