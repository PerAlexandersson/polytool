//! Fitting first-order vector recurrences, companion systems and
//! shared-denominator systems to observed state sequences.

use super::*;
use crate::linalg;
use num_traits::One;
use num_traits::Zero;

struct RectangularVectorFit {
    matrix: Vec<Vec<WeylOperator>>,
    forcing: Option<Vec<BivarPoly>>,
    diagnostics: Vec<VectorRecurrenceRowDiagnostics>,
}

fn fit_rectangular_vector_map_rational(
    source_states: &[Vec<Vec<BigRational>>],
    target_states: &[Vec<Vec<BigRational>>],
    source_indices: &[usize],
    fit_transitions: usize,
    opts: &VectorRecurrenceOptions,
    output_component_offset: usize,
) -> Result<RectangularVectorFit, VectorRecurrenceFitError> {
    let source_dimension = source_states.first().map_or(0, Vec::len);
    let output_dimension = target_states.first().map_or(0, Vec::len);
    let vars_per_bivar = (opts.idx_deg + 1) * (opts.var_deg + 1);
    let matrix_unknowns = source_dimension * (opts.diff_deg + 1) * vars_per_bivar;
    let forcing_unknowns = if opts.homogeneous {
        0
    } else {
        (opts.forcing_idx_deg + 1) * (opts.forcing_var_deg + 1)
    };
    let unknowns = matrix_unknowns + forcing_unknowns;

    let derivatives = source_states
        .iter()
        .map(|state| {
            state
                .iter()
                .map(|poly| {
                    (0..=opts.diff_deg)
                        .map(|d| poly_nth_derivative_rational(poly, d))
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let matrix_col = |source: usize, d: usize, i: usize, j: usize| {
        ((source * (opts.diff_deg + 1) + d) * (opts.idx_deg + 1) + i) * (opts.var_deg + 1) + j
    };
    let forcing_col = |i: usize, j: usize| matrix_unknowns + i * (opts.forcing_var_deg + 1) + j;

    let mut fitted_matrix = Vec::with_capacity(output_dimension);
    let mut fitted_forcing = (!opts.homogeneous).then(Vec::new);
    let mut diagnostics = Vec::with_capacity(output_dimension);

    for output in 0..output_dimension {
        let target_degree = target_states[..fit_transitions]
            .iter()
            .map(|state| poly_degree_rational(&state[output]))
            .max()
            .unwrap_or(0);
        let source_degree = derivatives[..fit_transitions]
            .iter()
            .flat_map(|state| state.iter())
            .flat_map(|orders| orders.iter())
            .map(|poly| poly_degree_rational(poly))
            .max()
            .unwrap_or(0)
            + opts.var_deg;
        let forcing_degree = if opts.homogeneous {
            0
        } else {
            opts.forcing_var_deg
        };
        let max_degree = target_degree.max(source_degree).max(forcing_degree);
        let equations_per_transition = max_degree + 1;
        let equations = fit_transitions * equations_per_transition;
        let mut matrix = vec![vec![BigRational::zero(); unknowns]; equations];
        let mut rhs = vec![BigRational::zero(); equations];

        for transition in 0..fit_transitions {
            let n_powers = rational_index_powers(
                source_indices[transition],
                opts.idx_deg.max(opts.forcing_idx_deg),
            );
            for degree in 0..=max_degree {
                let row = transition * equations_per_transition + degree;
                rhs[row] = poly_coeff_rational(&target_states[transition][output], degree);
                for source in 0..source_dimension {
                    for d in 0..=opts.diff_deg {
                        let derivative = &derivatives[transition][source][d];
                        for i in 0..=opts.idx_deg {
                            for j in 0..=opts.var_deg {
                                if degree < j {
                                    continue;
                                }
                                let coefficient = poly_coeff_rational(derivative, degree - j);
                                if !coefficient.is_zero() {
                                    matrix[row][matrix_col(source, d, i, j)] +=
                                        coefficient * n_powers[i].clone();
                                }
                            }
                        }
                    }
                }
                if !opts.homogeneous && degree <= opts.forcing_var_deg {
                    for i in 0..=opts.forcing_idx_deg {
                        matrix[row][forcing_col(i, degree)] += n_powers[i].clone();
                    }
                }
            }
        }

        let Some(solution) = linalg::solve_linear_system_with_rank(&matrix, &rhs) else {
            return Err(VectorRecurrenceFitError::NoSolution {
                output_component: output_component_offset + output,
            });
        };
        if opts.require_unique && solution.nullity != 0 {
            return Err(VectorRecurrenceFitError::Underdetermined {
                output_component: output_component_offset + output,
                rank: solution.rank,
                unknowns,
            });
        }

        let mut operator_row = Vec::with_capacity(source_dimension);
        for source in 0..source_dimension {
            let mut coefficients = Vec::with_capacity(opts.diff_deg + 1);
            for d in 0..=opts.diff_deg {
                let start = matrix_col(source, d, 0, 0);
                coefficients.push(extract_vector_bivar(
                    &solution.particular,
                    start,
                    opts.idx_deg,
                    opts.var_deg,
                ));
            }
            operator_row.push(WeylOperator::new(coefficients));
        }
        fitted_matrix.push(operator_row);
        if let Some(forcing) = &mut fitted_forcing {
            forcing.push(extract_vector_bivar(
                &solution.particular,
                matrix_unknowns,
                opts.forcing_idx_deg,
                opts.forcing_var_deg,
            ));
        }
        diagnostics.push(VectorRecurrenceRowDiagnostics {
            output_component: output_component_offset + output,
            equations,
            unknowns,
            rank: solution.rank,
            nullity: solution.nullity,
        });
    }

    for transition in 0..source_states.len() {
        let actual = evaluate_rectangular_vector_map_rational(
            &fitted_matrix,
            fitted_forcing.as_deref(),
            &source_states[transition],
            source_indices[transition],
        )
        .expect("fitted rectangular map has consistent dimensions");
        if !vector_states_equal(&actual, &target_states[transition]) {
            if transition >= fit_transitions {
                return Err(VectorRecurrenceFitError::HeldOutVerificationFailed {
                    transition: source_indices[transition],
                });
            }
            return Err(VectorRecurrenceFitError::NoSolution {
                output_component: output_component_offset,
            });
        }
    }

    Ok(RectangularVectorFit {
        matrix: fitted_matrix,
        forcing: fitted_forcing,
        diagnostics,
    })
}

/// Find an exact first-order vector recurrence from consecutive states.
///
/// `states[k]` is `F_{first_index+k}`. Matrix and forcing coefficients are
/// evaluated at the source index, so the returned convention is
/// `F_{n+1} = M(n,x,D_x) F_n + G(n,x)`. Each output row is solved separately;
/// by default every row must have full column rank, and the final transitions
/// are reserved for independent exact verification.
pub fn find_vector_recurrence_rational(
    states: &[Vec<Vec<BigRational>>],
    first_index: usize,
    opts: &VectorRecurrenceOptions,
) -> Result<VectorRecurrenceFit, VectorRecurrenceFitError> {
    let dimension = validate_vector_states(states)?;
    let transitions = states.len().saturating_sub(1);
    if transitions <= opts.held_out_transitions {
        return Err(VectorRecurrenceFitError::NotEnoughTransitions {
            transitions,
            held_out: opts.held_out_transitions,
        });
    }
    let fit_transitions = transitions - opts.held_out_transitions;
    let sources = states[..transitions].to_vec();
    let targets = states[1..].to_vec();
    let source_indices = (0..transitions)
        .map(|offset| {
            first_index
                .checked_add(offset)
                .ok_or(VectorRecurrenceFitError::IndexOverflow {
                    first_index,
                    offset,
                })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let fitted = fit_rectangular_vector_map_rational(
        &sources,
        &targets,
        &source_indices,
        fit_transitions,
        opts,
        0,
    )?;
    debug_assert_eq!(fitted.matrix.len(), dimension);
    Ok(VectorRecurrenceFit {
        recurrence: VectorRecurrence {
            matrix: fitted.matrix,
            forcing: fitted.forcing,
        },
        rows: fitted.diagnostics,
        fit_transitions,
        held_out_transitions: opts.held_out_transitions,
    })
}

/// Integer-coefficient convenience wrapper for
/// [`find_vector_recurrence_rational`].
pub fn find_vector_recurrence(
    states: &[Vec<Vec<i64>>],
    first_index: usize,
    opts: &VectorRecurrenceOptions,
) -> Result<VectorRecurrenceFit, VectorRecurrenceFitError> {
    find_vector_recurrence_rational(&i64_vector_states_to_rational(states), first_index, opts)
}

/// Find a lag-`lag` recurrence and return it as a first-order companion system.
///
/// Only the genuinely unknown final block row is fitted. The preceding block
/// rows are exact identity shifts, so their coefficients do not consume linear
/// unknowns or create artificial rank deficiencies. The companion state is
/// `(F_{n-lag+1}, ..., F_n)`.
pub fn find_companion_vector_recurrence_rational(
    states: &[Vec<Vec<BigRational>>],
    first_index: usize,
    lag: usize,
    opts: &VectorRecurrenceOptions,
) -> Result<VectorRecurrenceFit, VectorRecurrenceFitError> {
    if lag == 0 {
        return Err(VectorRecurrenceFitError::InvalidLag);
    }
    let component_dimension = validate_vector_states(states)?;
    let transitions = states.len().saturating_sub(lag);
    if transitions <= opts.held_out_transitions {
        return Err(VectorRecurrenceFitError::NotEnoughTransitions {
            transitions,
            held_out: opts.held_out_transitions,
        });
    }
    let fit_transitions = transitions - opts.held_out_transitions;
    let mut sources = Vec::with_capacity(transitions);
    let mut targets = Vec::with_capacity(transitions);
    let mut source_indices = Vec::with_capacity(transitions);
    for latest in lag - 1..states.len() - 1 {
        let mut companion = Vec::with_capacity(lag * component_dimension);
        for state in &states[latest + 1 - lag..=latest] {
            companion.extend(state.iter().cloned());
        }
        sources.push(companion);
        targets.push(states[latest + 1].clone());
        source_indices.push(first_index.checked_add(latest).ok_or(
            VectorRecurrenceFitError::IndexOverflow {
                first_index,
                offset: latest,
            },
        )?);
    }

    let output_offset = (lag - 1) * component_dimension;
    let fitted = fit_rectangular_vector_map_rational(
        &sources,
        &targets,
        &source_indices,
        fit_transitions,
        opts,
        output_offset,
    )?;
    let companion_dimension = lag * component_dimension;
    let mut matrix = vec![vec![WeylOperator::zero(); companion_dimension]; companion_dimension];
    for block in 0..lag - 1 {
        for component in 0..component_dimension {
            let row = block * component_dimension + component;
            let column = (block + 1) * component_dimension + component;
            matrix[row][column] = WeylOperator::identity();
        }
    }
    for (component, fitted_row) in fitted.matrix.into_iter().enumerate() {
        matrix[output_offset + component] = fitted_row;
    }
    let forcing = fitted.forcing.map(|fitted_forcing| {
        let mut forcing = vec![zero_bivar(); companion_dimension];
        for (component, polynomial) in fitted_forcing.into_iter().enumerate() {
            forcing[output_offset + component] = polynomial;
        }
        forcing
    });

    Ok(VectorRecurrenceFit {
        recurrence: VectorRecurrence { matrix, forcing },
        rows: fitted.diagnostics,
        fit_transitions,
        held_out_transitions: opts.held_out_transitions,
    })
}

/// Integer-coefficient convenience wrapper for
/// [`find_companion_vector_recurrence_rational`].
pub fn find_companion_vector_recurrence(
    states: &[Vec<Vec<i64>>],
    first_index: usize,
    lag: usize,
    opts: &VectorRecurrenceOptions,
) -> Result<VectorRecurrenceFit, VectorRecurrenceFitError> {
    find_companion_vector_recurrence_rational(
        &i64_vector_states_to_rational(states),
        first_index,
        lag,
        opts,
    )
}

/// Fit a common-denominator vector map from already aligned source and target
/// states.  The equations for all output components are solved jointly: the
/// denominator coefficients are shared, while matrix and forcing coefficients
/// are row-local.
pub(super) fn fit_shared_denominator_map_rational(
    source_states: &[Vec<Vec<BigRational>>],
    target_states: &[Vec<Vec<BigRational>>],
    source_indices: &[usize],
    fit_transitions: usize,
    opts: &VectorRecurrenceDenominatorOptions,
    output_component_offset: usize,
) -> Result<SharedDenominatorVectorRecurrenceFit, SharedDenominatorVectorRecurrenceFitError> {
    let source_dimension = source_states.first().map_or(0, Vec::len);
    let output_dimension = target_states.first().map_or(0, Vec::len);
    let denom_width = opts.denominator_var_deg + 1;
    let denominator_unknowns = (opts.denominator_idx_deg + 1) * denom_width;
    let vars_per_bivar = (opts.idx_deg + 1) * (opts.var_deg + 1);
    let matrix_unknowns = source_dimension * (opts.diff_deg + 1) * vars_per_bivar;
    let forcing_unknowns = if opts.homogeneous {
        0
    } else {
        (opts.forcing_idx_deg + 1) * (opts.forcing_var_deg + 1)
    };
    let row_unknowns = matrix_unknowns + forcing_unknowns;
    let unknowns = denominator_unknowns + output_dimension * row_unknowns;

    let derivatives = source_states
        .iter()
        .map(|state| {
            state
                .iter()
                .map(|poly| {
                    (0..=opts.diff_deg)
                        .map(|d| poly_nth_derivative_rational(poly, d))
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let local_matrix_col = |source: usize, d: usize, i: usize, j: usize| {
        ((source * (opts.diff_deg + 1) + d) * (opts.idx_deg + 1) + i) * (opts.var_deg + 1) + j
    };
    let local_forcing_col =
        |i: usize, j: usize| matrix_unknowns + i * (opts.forcing_var_deg + 1) + j;
    let denominator_col = |i: usize, j: usize| i * denom_width + j;
    let global_local_col =
        |output: usize, local: usize| denominator_unknowns + output * row_unknowns + local;

    let mut rows = Vec::<Vec<BigRational>>::new();
    let mut rhs = Vec::<BigRational>::new();
    let mut equations_by_output = Vec::with_capacity(output_dimension);
    for output in 0..output_dimension {
        let target_degree = target_states[..fit_transitions]
            .iter()
            .map(|state| poly_degree_rational(&state[output]))
            .max()
            .unwrap_or(0)
            + opts.denominator_var_deg;
        let source_degree = derivatives[..fit_transitions]
            .iter()
            .flat_map(|state| state.iter())
            .flat_map(|orders| orders.iter())
            .map(|poly| poly_degree_rational(poly))
            .max()
            .unwrap_or(0)
            + opts.var_deg;
        let forcing_degree = if opts.homogeneous {
            0
        } else {
            opts.forcing_var_deg
        };
        let max_degree = target_degree.max(source_degree).max(forcing_degree);
        let equations_per_transition = max_degree + 1;
        equations_by_output.push(fit_transitions * equations_per_transition);
        for transition in 0..fit_transitions {
            let n_powers = rational_index_powers(
                source_indices[transition],
                opts.denominator_idx_deg
                    .max(opts.idx_deg)
                    .max(opts.forcing_idx_deg),
            );
            for degree in 0..=max_degree {
                let mut row = vec![BigRational::zero(); unknowns];
                // q(n,x) * target(x), with q on the left of the identity.
                for i in 0..=opts.denominator_idx_deg {
                    for j in 0..=opts.denominator_var_deg {
                        if degree < j {
                            continue;
                        }
                        let coefficient =
                            poly_coeff_rational(&target_states[transition][output], degree - j);
                        if !coefficient.is_zero() {
                            row[denominator_col(i, j)] += coefficient * n_powers[i].clone();
                        }
                    }
                }
                // Subtract the Weyl operator action.
                for source in 0..source_dimension {
                    for d in 0..=opts.diff_deg {
                        let derivative = &derivatives[transition][source][d];
                        for i in 0..=opts.idx_deg {
                            for j in 0..=opts.var_deg {
                                if degree < j {
                                    continue;
                                }
                                let coefficient = poly_coeff_rational(derivative, degree - j);
                                if !coefficient.is_zero() {
                                    row[global_local_col(
                                        output,
                                        local_matrix_col(source, d, i, j),
                                    )] -= coefficient * n_powers[i].clone();
                                }
                            }
                        }
                    }
                }
                if !opts.homogeneous && degree <= opts.forcing_var_deg {
                    for i in 0..=opts.forcing_idx_deg {
                        row[global_local_col(output, local_forcing_col(i, degree))] -=
                            n_powers[i].clone();
                    }
                }
                rows.push(row);
                rhs.push(BigRational::zero());
            }
        }
    }

    let mut underdetermined = None;
    let mut accepted = None;
    let mut held_out_failure = None;
    for pivot in 0..denominator_unknowns {
        // Fix the first nonzero denominator monomial to one.  Earlier
        // denominator coefficients are fixed to zero; this removes the
        // homogeneous scale ambiguity without choosing a coefficient in
        // advance.
        let mut kept_columns = Vec::with_capacity(unknowns - pivot - 1);
        for column in pivot + 1..unknowns {
            kept_columns.push(column);
        }
        let mut constrained = Vec::with_capacity(rows.len());
        let mut constrained_rhs = Vec::with_capacity(rhs.len());
        for row in &rows {
            constrained.push(
                kept_columns
                    .iter()
                    .map(|&column| row[column].clone())
                    .collect::<Vec<_>>(),
            );
            constrained_rhs.push(-row[pivot].clone());
        }
        let Some(solution) = linalg::solve_linear_system_with_rank(&constrained, &constrained_rhs)
        else {
            continue;
        };
        if opts.require_unique && solution.nullity != 0 {
            underdetermined = Some((solution.rank, kept_columns.len()));
            continue;
        }

        let mut full_solution = vec![BigRational::zero(); unknowns];
        full_solution[pivot] = BigRational::one();
        for (value, &column) in solution.particular.iter().zip(&kept_columns) {
            full_solution[column] = value.clone();
        }

        let denominator = extract_vector_bivar(
            &full_solution,
            0,
            opts.denominator_idx_deg,
            opts.denominator_var_deg,
        );
        if denominator.is_zero() {
            continue;
        }
        // A zero evaluated denominator makes the undivided identity vacuous
        // and cannot determine the next state.  Such candidates are not
        // usable recurrences, even when they fit the algebraic identity.
        if source_indices
            .iter()
            .any(|&source_index| poly_is_zero_rational(&bivar_eval_n(&denominator, source_index)))
        {
            continue;
        }
        let mut fitted_matrix = Vec::with_capacity(output_dimension);
        let mut fitted_forcing = (!opts.homogeneous).then(Vec::new);
        for output in 0..output_dimension {
            let output_start = denominator_unknowns + output * row_unknowns;
            let mut operator_row = Vec::with_capacity(source_dimension);
            for source in 0..source_dimension {
                let coefficients = (0..=opts.diff_deg)
                    .map(|d| {
                        let start = output_start + local_matrix_col(source, d, 0, 0);
                        extract_vector_bivar(&full_solution, start, opts.idx_deg, opts.var_deg)
                    })
                    .collect();
                operator_row.push(WeylOperator::new(coefficients));
            }
            fitted_matrix.push(operator_row);
            if let Some(forcing) = &mut fitted_forcing {
                let start = output_start + matrix_unknowns;
                forcing.push(extract_vector_bivar(
                    &full_solution,
                    start,
                    opts.forcing_idx_deg,
                    opts.forcing_var_deg,
                ));
            }
        }
        let recurrence = SharedDenominatorVectorRecurrence {
            denominator: denominator.clone(),
            matrix: fitted_matrix,
            forcing: fitted_forcing,
        };
        let mut candidate_failure = None;
        let mut fit_failure = false;
        for transition in 0..source_states.len() {
            if !recurrence.holds_rectangular_identity_rational(
                &source_states[transition],
                &target_states[transition],
                source_indices[transition],
            ) {
                if transition >= fit_transitions {
                    candidate_failure.get_or_insert(source_indices[transition]);
                } else {
                    // The coefficient equations should already imply the
                    // fitted transitions; reject any truncation mistake.
                    fit_failure = true;
                }
            }
        }
        let candidate = SharedDenominatorVectorRecurrenceFit {
            recurrence,
            diagnostics: SharedDenominatorVectorRecurrenceDiagnostics {
                total_equations: equations_by_output.iter().sum(),
                equations_by_output: equations_by_output.clone(),
                unknowns: kept_columns.len(),
                rank: solution.rank,
                nullity: solution.nullity,
                output_components: (0..output_dimension)
                    .map(|output| output_component_offset + output)
                    .collect(),
            },
            fit_transitions,
            held_out_transitions: source_states.len() - fit_transitions,
            rank: solution.rank,
            nullity: solution.nullity,
            denominator_pivot: {
                let i = pivot / denom_width;
                let j = pivot % denom_width;
                (i, j)
            },
        };
        if fit_failure {
            continue;
        }
        if let Some(transition) = candidate_failure {
            held_out_failure.get_or_insert(transition);
            continue;
        }
        if opts.require_unique && accepted.is_some() {
            return Err(
                SharedDenominatorVectorRecurrenceFitError::SharedDenominatorUnderdetermined {
                    rank: solution.rank,
                    unknowns: kept_columns.len(),
                },
            );
        }
        if !opts.require_unique {
            return Ok(candidate);
        }
        accepted = Some(candidate);
    }
    if let Some((rank, unknowns)) = underdetermined {
        if opts.require_unique && accepted.is_some() {
            return Err(
                SharedDenominatorVectorRecurrenceFitError::SharedDenominatorUnderdetermined {
                    rank,
                    unknowns,
                },
            );
        }
    }
    if let Some(candidate) = accepted {
        return Ok(candidate);
    }
    if let Some(transition) = held_out_failure {
        return Err(
            SharedDenominatorVectorRecurrenceFitError::HeldOutVerificationFailed { transition },
        );
    }
    if let Some((rank, unknowns)) = underdetermined {
        return Err(
            SharedDenominatorVectorRecurrenceFitError::SharedDenominatorUnderdetermined {
                rank,
                unknowns,
            },
        );
    }
    Err(SharedDenominatorVectorRecurrenceFitError::NoSharedDenominator)
}

/// Find a first-order vector recurrence with one common polynomial factor
/// `q(n,x)` on the left.  The denominator's leading nonzero monomial in the
/// coefficient ordering `(n^0 x^0, n^0 x^1, ..., n^1 x^0, ...)` is normalized
/// to one; this is the only normalization imposed on the returned model.
pub fn find_vector_recurrence_with_denominator_rational(
    states: &[Vec<Vec<BigRational>>],
    first_index: usize,
    opts: &VectorRecurrenceDenominatorOptions,
) -> Result<SharedDenominatorVectorRecurrenceFit, SharedDenominatorVectorRecurrenceFitError> {
    let dimension = validate_vector_states(states)?;
    let transitions = states.len().saturating_sub(1);
    if transitions <= opts.held_out_transitions {
        return Err(VectorRecurrenceFitError::NotEnoughTransitions {
            transitions,
            held_out: opts.held_out_transitions,
        }
        .into());
    }
    let fit_transitions = transitions - opts.held_out_transitions;
    let source_indices = (0..transitions)
        .map(|offset| {
            first_index
                .checked_add(offset)
                .ok_or(VectorRecurrenceFitError::IndexOverflow {
                    first_index,
                    offset,
                })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let fitted = fit_shared_denominator_map_rational(
        &states[..transitions],
        &states[1..],
        &source_indices,
        fit_transitions,
        opts,
        0,
    )?;
    debug_assert_eq!(fitted.recurrence.dimension(), dimension);
    Ok(fitted)
}

/// Integer-coefficient convenience wrapper for
/// [`find_vector_recurrence_with_denominator_rational`].
pub fn find_vector_recurrence_with_denominator(
    states: &[Vec<Vec<i64>>],
    first_index: usize,
    opts: &VectorRecurrenceDenominatorOptions,
) -> Result<SharedDenominatorVectorRecurrenceFit, SharedDenominatorVectorRecurrenceFitError> {
    find_vector_recurrence_with_denominator_rational(
        &i64_vector_states_to_rational(states),
        first_index,
        opts,
    )
}

/// Fit a lagged recurrence and return the corresponding companion system with
/// a common denominator.  Shift rows use `q(n,x) I`, since the denominator
/// multiplies every output component of the companion state.
pub fn find_companion_vector_recurrence_with_denominator_rational(
    states: &[Vec<Vec<BigRational>>],
    first_index: usize,
    lag: usize,
    opts: &VectorRecurrenceDenominatorOptions,
) -> Result<SharedDenominatorVectorRecurrenceFit, SharedDenominatorVectorRecurrenceFitError> {
    if lag == 0 {
        return Err(VectorRecurrenceFitError::InvalidLag.into());
    }
    let component_dimension = validate_vector_states(states)?;
    let transitions = states.len().saturating_sub(lag);
    if transitions <= opts.held_out_transitions {
        return Err(VectorRecurrenceFitError::NotEnoughTransitions {
            transitions,
            held_out: opts.held_out_transitions,
        }
        .into());
    }
    let fit_transitions = transitions - opts.held_out_transitions;
    let mut sources = Vec::with_capacity(transitions);
    let mut targets = Vec::with_capacity(transitions);
    let mut source_indices = Vec::with_capacity(transitions);
    for latest in lag - 1..states.len() - 1 {
        let mut companion = Vec::with_capacity(lag * component_dimension);
        for state in &states[latest + 1 - lag..=latest] {
            companion.extend(state.iter().cloned());
        }
        sources.push(companion);
        targets.push(states[latest + 1].clone());
        source_indices.push(first_index.checked_add(latest).ok_or(
            VectorRecurrenceFitError::IndexOverflow {
                first_index,
                offset: latest,
            },
        )?);
    }
    let output_offset = (lag - 1) * component_dimension;
    let fitted = fit_shared_denominator_map_rational(
        &sources,
        &targets,
        &source_indices,
        fit_transitions,
        opts,
        output_offset,
    )?;
    let companion_dimension = lag * component_dimension;
    let denominator = fitted.recurrence.denominator.clone();
    let mut matrix = vec![vec![WeylOperator::zero(); companion_dimension]; companion_dimension];
    for block in 0..lag - 1 {
        for component in 0..component_dimension {
            let row = block * component_dimension + component;
            let column = (block + 1) * component_dimension + component;
            matrix[row][column] = WeylOperator::new(vec![denominator.clone()]);
        }
    }
    for (component, fitted_row) in fitted.recurrence.matrix.into_iter().enumerate() {
        matrix[output_offset + component] = fitted_row;
    }
    let forcing = fitted.recurrence.forcing.map(|fitted_forcing| {
        let mut forcing = vec![zero_bivar(); companion_dimension];
        for (component, polynomial) in fitted_forcing.into_iter().enumerate() {
            forcing[output_offset + component] = polynomial;
        }
        forcing
    });
    Ok(SharedDenominatorVectorRecurrenceFit {
        recurrence: SharedDenominatorVectorRecurrence {
            denominator,
            matrix,
            forcing,
        },
        diagnostics: fitted.diagnostics,
        fit_transitions,
        held_out_transitions: opts.held_out_transitions,
        rank: fitted.rank,
        nullity: fitted.nullity,
        denominator_pivot: fitted.denominator_pivot,
    })
}

/// Integer-coefficient convenience wrapper for
/// [`find_companion_vector_recurrence_with_denominator_rational`].
pub fn find_companion_vector_recurrence_with_denominator(
    states: &[Vec<Vec<i64>>],
    first_index: usize,
    lag: usize,
    opts: &VectorRecurrenceDenominatorOptions,
) -> Result<SharedDenominatorVectorRecurrenceFit, SharedDenominatorVectorRecurrenceFitError> {
    find_companion_vector_recurrence_with_denominator_rational(
        &i64_vector_states_to_rational(states),
        first_index,
        lag,
        opts,
    )
}

#[cfg(test)]
mod tests;
