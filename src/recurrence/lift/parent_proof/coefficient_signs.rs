//! Exact coefficient-sign checks and catalog attestations.

use super::*;

/// A bivariate polynomial `sum a[i][j] v^i w^j` (dense, exact).
type Bivariate = Vec<Vec<Q>>;

fn bivariate_add(left: &Bivariate, right: &Bivariate) -> Bivariate {
    let rows = left.len().max(right.len());
    (0..rows)
        .map(|i| {
            let a = left.get(i).map_or(&[][..], Vec::as_slice);
            let b = right.get(i).map_or(&[][..], Vec::as_slice);
            (0..a.len().max(b.len()))
                .map(|j| {
                    a.get(j).cloned().unwrap_or_default() + b.get(j).cloned().unwrap_or_default()
                })
                .collect()
        })
        .collect()
}

fn bivariate_mul(left: &Bivariate, right: &Bivariate) -> Bivariate {
    let mut out: Bivariate = Vec::new();
    for (i, row) in left.iter().enumerate() {
        for (j, a) in row.iter().enumerate() {
            if a.is_zero() {
                continue;
            }
            for (k, other) in right.iter().enumerate() {
                for (l, b) in other.iter().enumerate() {
                    if b.is_zero() {
                        continue;
                    }
                    if out.len() <= i + k {
                        out.resize(i + k + 1, Vec::new());
                    }
                    if out[i + k].len() <= j + l {
                        out[i + k].resize(j + l + 1, Q::zero());
                    }
                    out[i + k][j + l] += a * b;
                }
            }
        }
    }
    out
}

/// `p(a + b v + c w)` for a univariate `p`, as a bivariate polynomial in `(v, w)`.
fn bivariate_affine_powers(a: &Q, b: &Q, c: &Q, degree: usize) -> Vec<Bivariate> {
    let affine: Bivariate = vec![vec![a.clone(), c.clone()], vec![b.clone()]];
    let mut powers = vec![vec![vec![Q::one()]]];
    for _ in 0..degree {
        let next = bivariate_mul(powers.last().expect("nonempty"), &affine);
        powers.push(next);
    }
    powers
}

/// Stirling numbers of the second kind `S(m, i)`, `m, i <= max`.
fn stirling_second_kind_table(max: usize) -> Vec<Vec<Q>> {
    let mut table = vec![vec![Q::zero(); max + 1]; max + 1];
    table[0][0] = Q::one();
    for m in 1..=max {
        for i in 1..=m {
            table[m][i] = &table[m - 1][i - 1] + Q::from_integer(i.into()) * &table[m - 1][i];
        }
    }
    table
}

/// Whether every coefficient of `p` in the basis `v^(falling i) w^(falling j)`
/// is `>= 0`; then `p(v, w) >= 0` at all integers `v, w >= 0`
/// (`x^m = sum_i S(m, i) x^(falling i)`).
fn has_nonnegative_falling_factorial_coefficients(p: &Bivariate) -> bool {
    let max = p
        .len()
        .max(p.iter().map(Vec::len).max().unwrap_or(0))
        .max(1);
    let stirling = stirling_second_kind_table(max);
    let mut falling = vec![vec![Q::zero(); max]; max];
    for (m, row) in p.iter().enumerate() {
        for (l, a) in row.iter().enumerate() {
            if a.is_zero() {
                continue;
            }
            for i in 0..=m {
                for j in 0..=l {
                    let weight = &stirling[m][i] * &stirling[l][j];
                    if !weight.is_zero() {
                        falling[i][j] += a * weight;
                    }
                }
            }
        }
    }
    falling.iter().flatten().all(|value| !value.is_negative())
}

/// `x (x - 1) ... (x - d + 1)` in ascending powers.
fn falling_factorial_polynomial(d: usize) -> Vec<Q> {
    (0..d).fold(vec![Q::one()], |acc, k| {
        vec_poly::mul(&acc, &[-Q::from_integer(k.into()), Q::one()])
    })
}

/// Degree in `x` of `c(n, x)` (`None` for the zero polynomial).
fn x_degree(poly: &BivarPoly) -> Option<usize> {
    poly.coeffs
        .iter()
        .filter_map(|row| row.iter().rposition(|value| !value.is_zero()))
        .max()
}

/// The coefficient of `x^j` of `c(n, x)` as a polynomial in `n`.
fn x_coefficient(poly: &BivarPoly, j: usize) -> Vec<Q> {
    poly.coeffs
        .iter()
        .map(|row| row.get(j).cloned().unwrap_or_default())
        .collect()
}

fn integer_q(value: i64) -> Q {
    Q::from_integer(value.into())
}

/// Route 10: an exact proof that every row from `first_index` on has
/// coefficients of one sign (module documentation, "One-signed rows").
/// Returns what was shown, or why the sufficient test does not apply.
pub(super) fn one_signed_rows_for_all_indices(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
) -> Result<String, String> {
    let start = first_index + seeds.len();
    let sign = seeds
        .iter()
        .flat_map(|seed| seed.iter())
        .find(|value| !value.is_zero())
        .map_or(Q::one(), |value| {
            if value.is_negative() {
                -Q::one()
            } else {
                Q::one()
            }
        });
    let inhomogeneous = rule
        .inhomogeneous
        .as_ref()
        .filter(|g| g.coeffs.iter().flatten().any(|value| !value.is_zero()));
    if sign.is_negative() && inhomogeneous.is_some() {
        return Err("the seeds are nonpositive but the rule is inhomogeneous".into());
    }
    if let Some(position) = seeds
        .iter()
        .position(|seed| seed.iter().any(|value| (value * &sign).is_negative()))
    {
        return Err(format!("seed {position} has coefficients of both signs"));
    }
    let positive_from_start = |polynomial: &[Q], condition: IndexSignCondition| {
        decide_index_polynomial_sign(
            polynomial,
            &IntegerRange::starting_at(BigInt::from(start)),
            condition,
        )
        .holds()
    };
    if let Some(leader) = &rule.denominator {
        if x_degree(leader).is_some_and(|degree| degree > 0) {
            return Err("the leader depends on x".into());
        }
        if !positive_from_start(&x_coefficient(leader, 0), IndexSignCondition::Positive) {
            return Err(format!("the leader is not positive for every n >= {start}"));
        }
    }
    if rule
        .terms
        .iter()
        .any(|term| term.sign == RecurrenceSign::AlternatingN)
    {
        return Err("a term carries the factor (-1)^n".into());
    }
    if let Some(g) = inhomogeneous {
        for j in 0..=x_degree(g).unwrap_or(0) {
            if !positive_from_start(&x_coefficient(g, j), IndexSignCondition::Nonnegative) {
                return Err(format!(
                    "the coefficient of x^{j} of the inhomogeneous term is not >= 0 for every \
                     n >= {start}"
                ));
            }
        }
    }

    // Degree bound D(m) = alpha + beta m for every row m >= first_index.
    let mut beta = 0i64;
    for term in &rule.terms {
        let Some(degree) = x_degree(&term.coeff) else {
            continue;
        };
        let excess = i64::try_from(degree).unwrap_or(i64::MAX)
            - i64::try_from(term.deriv_order).unwrap_or(i64::MAX);
        let offset = i64::try_from(term.offset).unwrap_or(i64::MAX).max(1);
        if excess > 0 {
            beta = beta.max((excess + offset - 1) / offset);
        }
    }
    let as_i64 = |value: usize| i64::try_from(value).unwrap_or(i64::MAX);
    let mut alpha: Option<i64> = None;
    for (position, seed) in seeds.iter().enumerate() {
        if let Some(degree) = vec_poly::degree(seed) {
            let needed = as_i64(degree) - beta * as_i64(first_index + position);
            alpha = Some(alpha.map_or(needed, |alpha| alpha.max(needed)));
        }
    }
    if let Some(degree) = inhomogeneous.and_then(x_degree) {
        let needed = as_i64(degree) - beta * as_i64(start);
        alpha = Some(alpha.map_or(needed, |alpha| alpha.max(needed)));
    }
    let Some(alpha) = alpha else {
        return Ok("every seed is zero and the rule is homogeneous, so every row is zero".into());
    };

    // Combined weights W_(r,e)(n, s).
    let mut weights: std::collections::BTreeMap<(usize, i64), Bivariate> =
        std::collections::BTreeMap::new();
    for term in &rule.terms {
        let falling = falling_factorial_polynomial(term.deriv_order);
        for j in 0..=x_degree(&term.coeff).unwrap_or(0) {
            let c = x_coefficient(&term.coeff, j);
            if c.iter().all(Zero::is_zero) {
                continue;
            }
            let shift = as_i64(j) - as_i64(term.deriv_order);
            let product: Bivariate = c
                .iter()
                .map(|cn| falling.iter().map(|value| cn * value).collect())
                .collect();
            let entry = weights.entry((term.offset, shift)).or_default();
            *entry = bivariate_add(entry, &product);
        }
    }
    let start_q = integer_q(as_i64(start));
    for ((offset, shift), weight) in &weights {
        if weight.iter().flatten().all(Zero::is_zero) {
            continue;
        }
        let n_degree = weight.len().saturating_sub(1);
        let s_degree = weight
            .iter()
            .map(Vec::len)
            .max()
            .unwrap_or(1)
            .saturating_sub(1);
        let n_powers = bivariate_affine_powers(&start_q, &Q::one(), &Q::zero(), n_degree);
        let substitute = |s_powers: &[Bivariate]| {
            let mut total: Bivariate = Vec::new();
            for (a, row) in weight.iter().enumerate() {
                for (b, value) in row.iter().enumerate() {
                    if value.is_zero() {
                        continue;
                    }
                    let monomial = bivariate_mul(&n_powers[a], &s_powers[b]);
                    let scaled: Bivariate = monomial
                        .iter()
                        .map(|row| row.iter().map(|x| x * value).collect())
                        .collect();
                    total = bivariate_add(&total, &scaled);
                }
            }
            total
        };
        // (a) s >= 0 free.
        let free = bivariate_affine_powers(&Q::zero(), &Q::zero(), &Q::one(), s_degree);
        if has_nonnegative_falling_factorial_coefficients(&substitute(&free)) {
            continue;
        }
        // (b) s = D(n - r) - u with n = start + v.
        let base = alpha + beta * (as_i64(start) - as_i64(*offset));
        let bounded =
            bivariate_affine_powers(&integer_q(base), &integer_q(beta), &-Q::one(), s_degree);
        if has_nonnegative_falling_factorial_coefficients(&substitute(&bounded)) {
            continue;
        }
        return Err(format!(
            "the weight of the coefficient of x^s of P_(n-{offset}) in the coefficient of \
             x^(s{shift:+}) is not shown to be >= 0 for n >= {start}"
        ));
    }
    Ok(format!(
        "every row has coefficients of one sign: the seeds do, the leader is positive and every \
         combined weight is >= 0 for n >= {start}, with deg P_m <= {alpha} + {beta} m"
    ))
}

/// Attestations from bundled catalog entries flagged `lean_proved` whose
/// first index and leading initial rows equal the parent's (a cheap filter;
/// [`same_record`] decides the exact match).
pub(super) fn catalog_flag_attestations(
    record: &ParentRecurrenceRecord<'_>,
) -> Vec<VeroneseParentAttestation> {
    let wanted: Vec<Vec<String>> = record
        .seeds
        .iter()
        .map(|seed| {
            vec_poly::trim(seed)
                .iter()
                .map(format_rational_coeff)
                .collect()
        })
        .collect();
    let trimmed = |row: &[&str]| -> Vec<String> {
        let end = row
            .iter()
            .rposition(|value| *value != "0")
            .map_or(0, |i| i + 1);
        row[..end]
            .iter()
            .map(|value| (*value).to_string())
            .collect()
    };
    let zero_row = vec!["0".to_string()];
    let mut attestations = Vec::new();
    for (entry, flag) in entries_with_real_rootedness_flag() {
        if flag.status != OeisRealRootednessStatus::LeanProved
            || entry.recurrence_first_index != record.first_index
            || entry.initial_rows.len() < wanted.len()
        {
            continue;
        }
        let seeds_agree = entry.initial_rows.iter().zip(&wanted).all(|(row, seed)| {
            let row = trimmed(row);
            row == *seed || (row.is_empty() && *seed == zero_row)
        });
        if !seeds_agree {
            continue;
        }
        let Ok((rule, initial_rows)) = entry.recurrence_parts() else {
            continue;
        };
        let as_of = flag
            .as_of
            .map_or_else(String::new, |date| format!(", as of {date}"));
        attestations.push(VeroneseParentAttestation {
            label: format!("oeis_catalog:{}:lean_proved", entry.id),
            recurrence: RecurrenceJson::from_recurrence_rational(
                &rule,
                entry.recurrence_first_index,
                &initial_rows,
                None,
            ),
            first_index: entry.recurrence_first_index,
            statement: format!(
                "every row of the bundled OEIS catalog entry {} is zero or real-rooted (flag \
                 lean_proved{as_of}: a machine-checked Lean proof made outside polytool)",
                entry.id
            ),
            source: format!(
                "polytool OEIS catalog, real_rootedness flag of {}",
                entry.id
            ),
        });
    }
    attestations
}

/// First checked row that is nonzero and not real-rooted of one sign.
pub(super) fn first_failing_row(rows: &[Vec<Q>], first_index: usize) -> Option<usize> {
    rows.iter()
        .position(|row| {
            let shape = classify_row(row);
            !shape.is_zero && !(shape.has_one_sign && shape.is_real_rooted)
        })
        .map(|offset| first_index + offset)
}
