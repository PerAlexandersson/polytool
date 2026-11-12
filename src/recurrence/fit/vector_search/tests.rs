//! Tests for the enclosing module.

use super::*;
use crate::recurrence::test_support::*;

#[test]
fn vector_recurrence_fitting_reports_index_overflow() {
    let states = vec![
        vec![vec![br(1, 1)]],
        vec![vec![br(2, 1)]],
        vec![vec![br(4, 1)]],
    ];
    let options = VectorRecurrenceOptions {
        held_out_transitions: 0,
        ..Default::default()
    };
    assert_eq!(
        find_vector_recurrence_rational(&states, usize::MAX, &options),
        Err(VectorRecurrenceFitError::IndexOverflow {
            first_index: usize::MAX,
            offset: 1,
        })
    );
}

fn bounded_bivar(idx_deg: usize, var_deg: usize, terms: &[(usize, usize, i64)]) -> BivarPoly {
    let mut coeffs = vec![vec![BigRational::zero(); var_deg + 1]; idx_deg + 1];
    for &(i, j, coefficient) in terms {
        coeffs[i][j] = BigRational::from_integer(BigInt::from(coefficient));
    }
    BivarPoly { coeffs }
}

fn bounded_weyl(
    idx_deg: usize,
    var_deg: usize,
    diff_deg: usize,
    terms: &[(usize, usize, usize, i64)],
) -> WeylOperator {
    let mut coefficients = (0..=diff_deg)
        .map(|_| bounded_bivar(idx_deg, var_deg, &[]))
        .collect::<Vec<_>>();
    for &(d, i, j, coefficient) in terms {
        coefficients[d].coeffs[i][j] = BigRational::from_integer(BigInt::from(coefficient));
    }
    WeylOperator::new(coefficients)
}

#[test]
fn vector_affine_finder_recovers_q_integer_recurrence() {
    // [n+1]_x = x [n]_x + 1.
    let states = (0..=9)
        .map(|n| {
            let polynomial = if n == 0 { vec![0] } else { vec![1; n] };
            vec![polynomial]
        })
        .collect::<Vec<_>>();
    let opts = VectorRecurrenceOptions {
        var_deg: 1,
        idx_deg: 0,
        diff_deg: 0,
        homogeneous: false,
        forcing_var_deg: 0,
        forcing_idx_deg: 0,
        held_out_transitions: 2,
        require_unique: true,
    };

    let fit = find_vector_recurrence(&states, 0, &opts).unwrap();
    assert_eq!(fit.rows[0].nullity, 0);
    assert_eq!(
        fit.recurrence.matrix[0][0],
        bounded_weyl(0, 1, 0, &[(0, 0, 1, 1)])
    );
    assert_eq!(
        fit.recurrence.forcing,
        Some(vec![bounded_bivar(0, 0, &[(0, 0, 1)])])
    );
}

#[test]
fn vector_finder_checks_complete_held_out_transitions() {
    let mut states = (0..=9)
        .map(|n| {
            let polynomial = if n == 0 { vec![0] } else { vec![1; n] };
            vec![polynomial]
        })
        .collect::<Vec<_>>();
    states[9][0][0] = 2;
    let opts = VectorRecurrenceOptions {
        var_deg: 1,
        idx_deg: 0,
        diff_deg: 0,
        homogeneous: false,
        forcing_var_deg: 0,
        forcing_idx_deg: 0,
        held_out_transitions: 2,
        require_unique: true,
    };

    assert_eq!(
        find_vector_recurrence(&states, 0, &opts),
        Err(VectorRecurrenceFitError::HeldOutVerificationFailed { transition: 8 })
    );
}

#[test]
fn vector_finder_rejects_nonidentifiable_duplicate_channels() {
    let states = (0..=7)
        .map(|n| {
            let value = 1_i64 << n;
            vec![vec![value], vec![value]]
        })
        .collect::<Vec<_>>();
    let opts = VectorRecurrenceOptions {
        var_deg: 0,
        idx_deg: 0,
        diff_deg: 0,
        held_out_transitions: 2,
        ..Default::default()
    };

    assert_eq!(
        find_vector_recurrence(&states, 0, &opts),
        Err(VectorRecurrenceFitError::Underdetermined {
            output_component: 0,
            rank: 1,
            unknowns: 2,
        })
    );
}

#[test]
fn vector_finder_recovers_published_eulerian_pair() {
    // Ma--Ma--Yeh--Yeh, Discrete Math. 345 (2022), Eq. (8): the
    // even/odd up-down-run polynomials form this Eulerian pair.
    let expected = VectorRecurrence {
        matrix: vec![
            vec![
                bounded_weyl(1, 2, 1, &[(0, 1, 1, 1), (1, 0, 1, 2), (1, 0, 2, -2)]),
                bounded_weyl(1, 2, 1, &[(0, 0, 1, 1)]),
            ],
            vec![
                bounded_weyl(1, 2, 1, &[(0, 0, 0, 1)]),
                bounded_weyl(
                    1,
                    2,
                    1,
                    &[
                        (0, 0, 0, 1),
                        (0, 0, 1, -1),
                        (0, 1, 1, 1),
                        (1, 0, 1, 2),
                        (1, 0, 2, -2),
                    ],
                ),
            ],
        ],
        forcing: None,
    };
    let mut states = vec![vec![vec![BigRational::zero()], vec![BigRational::one()]]];
    for n in 1..=18 {
        let next = expected
            .evaluate_next_rational(states.last().unwrap(), n)
            .unwrap();
        states.push(next);
    }
    let opts = VectorRecurrenceOptions {
        var_deg: 2,
        idx_deg: 1,
        diff_deg: 1,
        held_out_transitions: 3,
        ..Default::default()
    };

    let fit = find_vector_recurrence_rational(&states, 1, &opts).unwrap();
    assert!(fit.rows.iter().all(|row| row.nullity == 0));
    assert_eq!(fit.recurrence, expected);
}

#[test]
fn companion_finder_recovers_published_type_b_lag_two_system() {
    // Ma et al., EJC 27(3) (2020), Eq. (12), specialized to k=1.
    // The state order is (P_{n-1}, Q_{n-1}, P_n, Q_n).
    let z = || bounded_weyl(1, 2, 1, &[]);
    let expected = VectorRecurrence {
        matrix: vec![
            vec![z(), z(), WeylOperator::identity(), z()],
            vec![z(), z(), z(), WeylOperator::identity()],
            vec![
                bounded_weyl(1, 2, 1, &[(0, 1, 1, 2)]),
                z(),
                bounded_weyl(
                    1,
                    2,
                    1,
                    &[
                        (0, 0, 0, 1),
                        (0, 0, 1, -1),
                        (0, 1, 1, 2),
                        (1, 0, 1, 2),
                        (1, 0, 2, -2),
                    ],
                ),
                bounded_weyl(1, 2, 1, &[(0, 0, 0, 1)]),
            ],
            vec![
                z(),
                bounded_weyl(1, 2, 1, &[(0, 1, 1, 2)]),
                bounded_weyl(1, 2, 1, &[(0, 0, 1, 1)]),
                bounded_weyl(1, 2, 1, &[(0, 1, 1, 2), (1, 0, 1, 2), (1, 0, 2, -2)]),
            ],
        ],
        forcing: None,
    };
    let mut states = vec![
        vec![vec![BigRational::zero()], vec![BigRational::one()]],
        vec![vec![BigRational::one()], vec![BigRational::zero()]],
    ];
    for n in 1..=20 {
        let companion = vec![
            states[n - 1][0].clone(),
            states[n - 1][1].clone(),
            states[n][0].clone(),
            states[n][1].clone(),
        ];
        let next = expected.evaluate_next_rational(&companion, n).unwrap();
        states.push(vec![next[2].clone(), next[3].clone()]);
    }
    let opts = VectorRecurrenceOptions {
        var_deg: 2,
        idx_deg: 1,
        diff_deg: 1,
        held_out_transitions: 3,
        ..Default::default()
    };

    let fit = find_companion_vector_recurrence_rational(&states, 0, 2, &opts).unwrap();
    assert!(fit.rows.iter().all(|row| row.nullity == 0));
    assert_eq!(fit.rows.len(), 2, "only the final block row is fitted");
    assert_eq!(fit.recurrence.matrix[0][2], WeylOperator::identity());
    assert_eq!(fit.recurrence.matrix[1][3], WeylOperator::identity());
    for row in 2..4 {
        assert_eq!(fit.recurrence.matrix[row], expected.matrix[row]);
    }
}

#[test]
fn companion_i64_finder_recovers_fibonacci_and_transitions_hold() {
    // F_{n+1} = F_n + F_{n-1} (OEIS A000045) as a lag-2 scalar system: the
    // companion state (F_{n-1}, F_n) maps to (F_n, F_{n+1}).
    let mut fib = vec![1_i64, 1];
    for n in 2..14 {
        fib.push(fib[n - 1] + fib[n - 2]);
    }
    let states: Vec<Vec<Vec<i64>>> = fib.iter().map(|&f| vec![vec![f]]).collect();
    let opts = VectorRecurrenceOptions {
        var_deg: 0,
        idx_deg: 0,
        diff_deg: 0,
        held_out_transitions: 2,
        ..Default::default()
    };
    let fit = find_companion_vector_recurrence(&states, 0, 2, &opts).unwrap();
    // The integer wrapper agrees with the rational entry point.
    let rational = find_companion_vector_recurrence_rational(
        &i64_vector_states_to_rational(&states),
        0,
        2,
        &opts,
    )
    .unwrap();
    assert_eq!(fit, rational);
    assert!(fit.rows.iter().all(|row| row.nullity == 0));
    let one = || bounded_weyl(0, 0, 0, &[(0, 0, 0, 1)]);
    assert!(fit.recurrence.matrix[0][0].is_zero());
    assert_eq!(fit.recurrence.matrix[0][1], WeylOperator::identity());
    assert_eq!(fit.recurrence.matrix[1], vec![one(), one()]);
    assert_eq!(fit.recurrence.forcing, None);
    // Every consecutive companion transition (including held-out ones) holds exactly,
    // and a perturbed target does not.
    let rational_states = i64_vector_states_to_rational(&states);
    for n in 1..fib.len() - 1 {
        let companion = vec![
            rational_states[n - 1][0].clone(),
            rational_states[n][0].clone(),
        ];
        let next = vec![
            rational_states[n][0].clone(),
            rational_states[n + 1][0].clone(),
        ];
        assert!(fit
            .recurrence
            .holds_transition_rational(&companion, &next, n));
        let mut wrong = next.clone();
        wrong[1][0] += BigRational::one();
        assert!(!fit
            .recurrence
            .holds_transition_rational(&companion, &wrong, n));
    }
    // A state of the wrong dimension is rejected instead of panicking.
    assert!(!fit
        .recurrence
        .holds_transition_rational(&rational_states[1], &rational_states[2], 1));
    // Lag zero is invalid.
    assert_eq!(
        find_companion_vector_recurrence(&states, 0, 0, &opts),
        Err(VectorRecurrenceFitError::InvalidLag)
    );
}

#[test]
fn denominator_i64_finders_agree_with_rational_and_recover_known_recurrences() {
    // Catalan numbers (OEIS A000108): (n + 2) C_{n+1} = (4 n + 2) C_n.  The
    // normalized denominator (constant term 1) is q = 1 + n / 2.
    let catalan: [i64; 14] = [
        1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796, 58786, 208012, 742900,
    ];
    let states: Vec<Vec<Vec<i64>>> = catalan.iter().map(|&c| vec![vec![c]]).collect();
    let rational_states = i64_vector_states_to_rational(&states);
    let opts = VectorRecurrenceDenominatorOptions {
        var_deg: 0,
        idx_deg: 1,
        diff_deg: 0,
        homogeneous: true,
        forcing_var_deg: 0,
        forcing_idx_deg: 0,
        held_out_transitions: 2,
        require_unique: true,
        denominator_var_deg: 0,
        denominator_idx_deg: 1,
    };
    let fit = find_vector_recurrence_with_denominator(&states, 0, &opts).unwrap();
    assert_eq!(
        fit,
        find_vector_recurrence_with_denominator_rational(&rational_states, 0, &opts).unwrap()
    );
    assert_eq!(
        fit.recurrence.denominator,
        BivarPoly {
            coeffs: vec![vec![br(1, 1)], vec![br(1, 2)]],
        }
    );
    assert_eq!(fit.nullity, 0);
    for n in 0..catalan.len() - 1 {
        assert!(fit.recurrence.holds_transition_rational(
            &rational_states[n],
            &rational_states[n + 1],
            n
        ));
    }

    // Central Delannoy numbers (OEIS A001850) need two previous terms:
    // (n + 1) a_{n+1} = (6 n + 3) a_n - n a_{n-1}.  In the lag-2 companion the
    // shift row carries the denominator q = 1 + n.
    let mut delannoy = vec![1_i64, 3];
    for n in 1..13 {
        delannoy.push(
            (3 * (2 * n + 1) * delannoy[n as usize] - n * delannoy[n as usize - 1]) / (n + 1),
        );
    }
    assert_eq!(&delannoy[..6], &[1, 3, 13, 63, 321, 1683]);
    let states: Vec<Vec<Vec<i64>>> = delannoy.iter().map(|&c| vec![vec![c]]).collect();
    let rational_states = i64_vector_states_to_rational(&states);
    let opts = VectorRecurrenceDenominatorOptions {
        idx_deg: 1,
        denominator_idx_deg: 1,
        ..opts
    };
    let companion =
        find_companion_vector_recurrence_with_denominator(&states, 0, 2, &opts).unwrap();
    assert_eq!(
        companion,
        find_companion_vector_recurrence_with_denominator_rational(&rational_states, 0, 2, &opts)
            .unwrap()
    );
    let denominator = bounded_bivar(1, 0, &[(0, 0, 1), (1, 0, 1)]);
    assert_eq!(companion.recurrence.denominator, denominator);
    assert!(companion.recurrence.matrix[0][0].is_zero());
    assert_eq!(
        companion.recurrence.matrix[0][1],
        WeylOperator::new(vec![denominator])
    );
    assert_eq!(
        companion.recurrence.matrix[1][0].coefficients[0],
        bounded_bivar(1, 0, &[(1, 0, -1)])
    );
    assert_eq!(
        companion.recurrence.matrix[1][1].coefficients[0],
        bounded_bivar(1, 0, &[(0, 0, 3), (1, 0, 6)])
    );
    for n in 1..delannoy.len() - 1 {
        let source = vec![
            rational_states[n - 1][0].clone(),
            rational_states[n][0].clone(),
        ];
        let target = vec![
            rational_states[n][0].clone(),
            rational_states[n + 1][0].clone(),
        ];
        assert!(companion
            .recurrence
            .holds_transition_rational(&source, &target, n));
    }
    assert_eq!(
        find_companion_vector_recurrence_with_denominator(&states, 0, 0, &opts),
        Err(VectorRecurrenceFitError::InvalidLag.into())
    );
}

#[test]
fn shared_denominator_finder_recovers_n_dependent_factor_and_forcing() {
    // n P_(n+1) = (n+1) P_n, with P_n = n.  The normalized recurrence
    // has q(n) = n; q = 1 would require a non-polynomial coefficient.
    let states = (1..=12)
        .map(|n| vec![vec![br(n as i64, 1)]])
        .collect::<Vec<_>>();
    let opts = VectorRecurrenceDenominatorOptions {
        var_deg: 0,
        idx_deg: 1,
        diff_deg: 0,
        homogeneous: true,
        forcing_var_deg: 0,
        forcing_idx_deg: 0,
        held_out_transitions: 2,
        require_unique: true,
        denominator_var_deg: 0,
        denominator_idx_deg: 1,
    };
    let fit = find_vector_recurrence_with_denominator_rational(&states, 1, &opts).unwrap();
    assert_eq!(fit.denominator_pivot, (1, 0));
    assert_eq!(
        fit.recurrence.denominator,
        bounded_bivar(1, 0, &[(1, 0, 1)])
    );
    assert_eq!(
        fit.recurrence.matrix[0][0].coefficients[0],
        bounded_bivar(1, 0, &[(0, 0, 1), (1, 0, 1)])
    );
    assert_eq!(fit.recurrence.forcing, None);
    assert_eq!(fit.diagnostics.nullity, 0);
}

#[test]
fn shared_denominator_finder_recovers_x_factor_and_rejects_tampering() {
    // x P_(n+1) = P_n for P_n = x^(N-n), so no normalized q=1 model
    // exists with a degree-zero matrix coefficient.
    let states = (0..=10)
        .map(|n| {
            let mut polynomial = vec![BigRational::zero(); 11 - n];
            polynomial[10 - n] = BigRational::one();
            vec![polynomial]
        })
        .collect::<Vec<_>>();
    let opts = VectorRecurrenceDenominatorOptions {
        var_deg: 0,
        idx_deg: 0,
        diff_deg: 0,
        homogeneous: true,
        forcing_var_deg: 0,
        forcing_idx_deg: 0,
        held_out_transitions: 2,
        require_unique: true,
        denominator_var_deg: 1,
        denominator_idx_deg: 0,
    };
    let fit = find_vector_recurrence_with_denominator_rational(&states, 0, &opts).unwrap();
    assert_eq!(fit.denominator_pivot, (0, 1));
    assert_eq!(
        fit.recurrence.denominator,
        bounded_bivar(0, 1, &[(0, 1, 1)])
    );
    assert!(fit
        .recurrence
        .holds_transition_rational(&states[0], &states[1], 0));
    let mut tampered = states.clone();
    tampered[9][0][0] += BigRational::one();
    assert_eq!(
        find_vector_recurrence_with_denominator_rational(&tampered, 0, &opts),
        Err(SharedDenominatorVectorRecurrenceFitError::HeldOutVerificationFailed { transition: 8 })
    );
}

#[test]
fn shared_denominator_finder_supports_affine_forcing() {
    // For P_n = 1 + x^(N-n), x P_(n+1) = P_n + (x - 1).
    let states = (0..=9)
        .map(|n| {
            let mut polynomial = vec![BigRational::zero(); 11 - n];
            polynomial[0] = BigRational::one();
            polynomial[10 - n] = BigRational::one();
            vec![polynomial]
        })
        .collect::<Vec<_>>();
    let opts = VectorRecurrenceDenominatorOptions {
        var_deg: 0,
        idx_deg: 0,
        diff_deg: 0,
        homogeneous: false,
        forcing_var_deg: 1,
        forcing_idx_deg: 0,
        held_out_transitions: 2,
        require_unique: true,
        denominator_var_deg: 1,
        denominator_idx_deg: 0,
    };
    let fit = find_vector_recurrence_with_denominator_rational(&states, 0, &opts).unwrap();
    assert_eq!(
        fit.recurrence.denominator,
        bounded_bivar(0, 1, &[(0, 1, 1)])
    );
    assert_eq!(
        fit.recurrence.forcing,
        Some(vec![bounded_bivar(0, 1, &[(0, 0, -1), (0, 1, 1)])])
    );
}

#[test]
fn shared_denominator_finder_recovers_two_distinct_rows() {
    let states = (0..=9)
        .map(|n| {
            let degree = 9 - n;
            let mut first = vec![BigRational::zero(); degree + 1];
            first[degree] = BigRational::one();
            let mut second = first.clone();
            second[degree] *= BigRational::from_integer(BigInt::from(2_i64.pow(n as u32)));
            vec![first, second]
        })
        .collect::<Vec<_>>();
    let opts = VectorRecurrenceDenominatorOptions {
        var_deg: 0,
        idx_deg: 0,
        diff_deg: 0,
        homogeneous: true,
        forcing_var_deg: 0,
        forcing_idx_deg: 0,
        held_out_transitions: 2,
        require_unique: true,
        denominator_var_deg: 1,
        denominator_idx_deg: 0,
    };
    let fit = find_vector_recurrence_with_denominator_rational(&states, 0, &opts).unwrap();
    assert_eq!(
        fit.recurrence.denominator,
        bounded_bivar(0, 1, &[(0, 1, 1)])
    );
    assert_eq!(
        fit.recurrence.matrix[0][0].coefficients[0],
        bounded_bivar(0, 0, &[(0, 0, 1)])
    );
    assert_eq!(
        fit.recurrence.matrix[1][1].coefficients[0],
        bounded_bivar(0, 0, &[(0, 0, 2)])
    );
    assert!(fit
        .recurrence
        .holds_transition_rational(&states[7], &states[8], 7));
    assert!(fit
        .recurrence
        .holds_transition_rational(&states[8], &states[9], 8));
}

#[test]
fn shared_denominator_finder_rejects_zero_at_observed_index() {
    // q(n)=n would make the n=0 identity vacuous and hide the arbitrary
    // transition 1 -> 2.  Such a candidate is not evaluable and must be
    // rejected rather than accepted from the undivided equations.
    let states = vec![
        vec![vec![br(1, 1)]],
        vec![vec![br(2, 1)]],
        vec![vec![br(0, 1)]],
        vec![vec![br(0, 1)]],
    ];
    let opts = VectorRecurrenceDenominatorOptions {
        var_deg: 0,
        idx_deg: 0,
        diff_deg: 0,
        homogeneous: true,
        forcing_var_deg: 0,
        forcing_idx_deg: 0,
        held_out_transitions: 1,
        require_unique: true,
        denominator_var_deg: 0,
        denominator_idx_deg: 1,
    };
    assert_eq!(
        find_vector_recurrence_with_denominator_rational(&states, 0, &opts),
        Err(SharedDenominatorVectorRecurrenceFitError::NoSharedDenominator)
    );
}

#[test]
fn shared_denominator_q_one_matches_normalized_vector_fit() {
    let states = (0..=8)
        .map(|n| {
            vec![vec![BigRational::from_integer(BigInt::from(
                2_i64.pow(n as u32),
            ))]]
        })
        .collect::<Vec<_>>();
    let normalized_opts = VectorRecurrenceOptions {
        var_deg: 0,
        idx_deg: 0,
        diff_deg: 0,
        homogeneous: true,
        held_out_transitions: 2,
        require_unique: true,
        ..Default::default()
    };
    let denominator_opts = VectorRecurrenceDenominatorOptions {
        var_deg: 0,
        idx_deg: 0,
        diff_deg: 0,
        homogeneous: true,
        forcing_var_deg: 0,
        forcing_idx_deg: 0,
        held_out_transitions: 2,
        require_unique: true,
        denominator_var_deg: 0,
        denominator_idx_deg: 0,
    };
    let normalized = find_vector_recurrence_rational(&states, 0, &normalized_opts).unwrap();
    let denominator =
        find_vector_recurrence_with_denominator_rational(&states, 0, &denominator_opts).unwrap();
    assert!(denominator.recurrence.denominator.is_one());
    assert_eq!(denominator.recurrence.matrix, normalized.recurrence.matrix);
    assert!(denominator
        .recurrence
        .holds_transition_rational(&states[7], &states[8], 7));
}

#[test]
fn shared_denominator_finder_reports_underdetermination() {
    let states = (1..=6)
        .map(|n| vec![vec![br(n as i64, 1)]])
        .collect::<Vec<_>>();
    let opts = VectorRecurrenceDenominatorOptions {
        var_deg: 0,
        idx_deg: 2,
        diff_deg: 0,
        homogeneous: true,
        forcing_var_deg: 0,
        forcing_idx_deg: 0,
        held_out_transitions: 0,
        require_unique: true,
        denominator_var_deg: 1,
        denominator_idx_deg: 2,
    };
    let result = find_vector_recurrence_with_denominator_rational(&states, 1, &opts);
    assert!(matches!(
        result,
        Err(SharedDenominatorVectorRecurrenceFitError::SharedDenominatorUnderdetermined { .. })
    ));
}

#[test]
fn shared_denominator_companion_uses_q_times_shift_identity() {
    let states = (0..=9)
        .map(|n| {
            let mut polynomial = vec![BigRational::zero(); 10 - n];
            polynomial[9 - n] = BigRational::one();
            vec![polynomial]
        })
        .collect::<Vec<_>>();
    let opts = VectorRecurrenceDenominatorOptions {
        var_deg: 0,
        idx_deg: 0,
        diff_deg: 0,
        homogeneous: true,
        forcing_var_deg: 0,
        forcing_idx_deg: 0,
        held_out_transitions: 0,
        require_unique: true,
        denominator_var_deg: 1,
        denominator_idx_deg: 0,
    };
    let fit =
        find_companion_vector_recurrence_with_denominator_rational(&states, 0, 2, &opts).unwrap();
    assert_eq!(
        fit.recurrence.matrix[0][1].coefficients[0],
        fit.recurrence.denominator
    );
    assert!(fit.recurrence.holds_transition_rational(
        &[states[0][0].clone(), states[1][0].clone()],
        &[states[1][0].clone(), states[2][0].clone()],
        1,
    ));
}
