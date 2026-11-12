//! Cyclotomic polynomials and exact root-of-unity evaluations.

use super::*;

pub(super) fn divisors(n: usize) -> Vec<usize> {
    let mut out = Vec::new();
    for d in 1..=n {
        if n % d == 0 {
            out.push(d);
        }
    }
    out
}

fn exact_div_monic(mut numerator: Vec<BigInt>, divisor: &[BigInt]) -> Option<Vec<BigInt>> {
    strip_trailing_zeros(&mut numerator);
    let dd = degree(divisor)?;
    if divisor[dd] != BigInt::one() {
        return None;
    }
    let nd = match degree(&numerator) {
        Some(nd) => nd,
        None => return Some(Vec::new()),
    };
    if nd < dd {
        return Some(Vec::new());
    }

    let mut quotient = vec![BigInt::zero(); nd - dd + 1];
    for shift in (0..=nd - dd).rev() {
        let coeff = numerator[dd + shift].clone();
        if coeff.is_zero() {
            continue;
        }
        quotient[shift] = coeff.clone();
        for (j, dcoeff) in divisor.iter().enumerate().take(dd + 1) {
            numerator[j + shift] -= &coeff * dcoeff;
        }
    }
    if is_zero(&numerator) {
        strip_trailing_zeros(&mut quotient);
        Some(quotient)
    } else {
        None
    }
}

pub(super) fn remainder_monic(mut numerator: Vec<BigInt>, divisor: &[BigInt]) -> Vec<BigInt> {
    strip_trailing_zeros(&mut numerator);
    let Some(dd) = degree(divisor) else {
        panic!("division by zero polynomial");
    };
    assert_eq!(divisor[dd], BigInt::one(), "divisor must be monic");
    while let Some(nd) = degree(&numerator) {
        if nd < dd {
            break;
        }
        let shift = nd - dd;
        let coeff = numerator[nd].clone();
        if coeff.is_zero() {
            numerator.pop();
            continue;
        }
        for (j, dcoeff) in divisor.iter().enumerate().take(dd + 1) {
            numerator[j + shift] -= &coeff * dcoeff;
        }
        strip_trailing_zeros(&mut numerator);
    }
    numerator
}

fn cyclotomic_cached(n: usize, cache: &mut BTreeMap<usize, Vec<BigInt>>) -> Vec<BigInt> {
    if let Some(poly) = cache.get(&n) {
        return poly.clone();
    }
    let poly = if n == 1 {
        vec![BigInt::from(-1), BigInt::one()]
    } else {
        let mut p = vec![BigInt::zero(); n + 1];
        p[0] = BigInt::from(-1);
        p[n] = BigInt::one();
        for d in divisors(n).into_iter().filter(|&d| d < n) {
            let phi_d = cyclotomic_cached(d, cache);
            p = exact_div_monic(p, &phi_d).expect("cyclotomic division is exact");
        }
        p
    };
    cache.insert(n, poly.clone());
    poly
}

/// Return the `n`-th cyclotomic polynomial in ascending coefficient order.
///
/// # Panics
///
/// Panics if `n == 0`.
pub fn cyclotomic_polynomial_bigint(n: usize) -> Vec<BigInt> {
    assert!(n >= 1, "cyclotomic index must be positive");
    cyclotomic_cached(n, &mut BTreeMap::new())
}

fn evaluate_at_one(coeffs: &[BigInt]) -> BigInt {
    coeffs.iter().fold(BigInt::zero(), |acc, coeff| acc + coeff)
}

/// Exactly reduce `P(q)` at `q = zeta_order^power`.
///
/// `power` is taken modulo `group_order`.  The result records the order of the
/// root of unity, the remainder modulo the cyclotomic polynomial and, when that
/// remainder is constant, the integer value `P(zeta)`.
///
/// # Panics
///
/// Panics if `group_order == 0`.
pub fn root_of_unity_evaluation_bigint(
    coeffs: &[BigInt],
    group_order: usize,
    power: usize,
) -> RootOfUnityEvaluation {
    assert!(group_order >= 1, "group order must be positive");
    let power = power % group_order;
    let root_order = root_order_of_power(group_order, power);
    let (remainder, integer_value) =
        evaluate_at_root_order(coeffs, root_order, &mut BTreeMap::new());
    RootOfUnityEvaluation {
        group_order,
        power,
        root_order,
        integer_value,
        remainder,
    }
}

/// Order of `zeta_group_order^power` for `power < group_order`.
fn root_order_of_power(group_order: usize, power: usize) -> usize {
    if power == 0 {
        1
    } else {
        group_order / group_order.gcd(&power)
    }
}

/// Remainder of `P` modulo `Phi_root_order` and its integer value when the
/// remainder is constant.  `cyclotomic` caches the cyclotomic polynomials.
fn evaluate_at_root_order(
    coeffs: &[BigInt],
    root_order: usize,
    cyclotomic: &mut BTreeMap<usize, Vec<BigInt>>,
) -> (Vec<BigInt>, Option<BigInt>) {
    if root_order == 1 {
        let value = evaluate_at_one(coeffs);
        return (vec![value.clone()], Some(value));
    }
    let phi = cyclotomic_cached(root_order, cyclotomic);
    let remainder = remainder_monic(coeffs.to_vec(), &phi);
    let integer_value = match remainder.as_slice() {
        [] => Some(BigInt::zero()),
        [constant] => Some(constant.clone()),
        _ => None,
    };
    (remainder, integer_value)
}

/// Check one CSP candidate.  If `fixed_counts` is `None`, this only profiles
/// root-of-unity evaluations.
///
/// `fixed_counts`, when given, must have exactly `order` entries with entry `d`
/// equal to `fix(g^d)`; otherwise (or for `order == 0`) an `Err` message is
/// returned.  The CSP holds iff `P(zeta_order^d) = fix(g^d)` for every `d`.
///
/// # Errors
///
/// [`CyclicSievingError::OrderTooLarge`] when `order` exceeds
/// [`CYCLIC_SIEVING_MAX_ORDER`] (checked before anything is allocated).
pub fn cyclic_sieving_report_bigint(
    coeffs: &[BigInt],
    order: usize,
    fixed_counts: Option<&[BigInt]>,
) -> Result<CyclicSievingReport, CyclicSievingError> {
    if order == 0 {
        return Err(CyclicSievingError::NonpositiveOrder);
    }
    if order > CYCLIC_SIEVING_MAX_ORDER {
        return Err(CyclicSievingError::OrderTooLarge {
            order: order as u128,
            maximum: CYCLIC_SIEVING_MAX_ORDER,
        });
    }
    if let Some(counts) = fixed_counts {
        if counts.len() != order {
            return Err(CyclicSievingError::FixedPointCountLength {
                order,
                found: counts.len(),
            });
        }
    }

    let mut coefficients = coeffs.to_vec();
    strip_trailing_zeros(&mut coefficients);
    let mut evaluations = Vec::with_capacity(order);
    let mut checks = Vec::new();
    // The evaluation depends only on the order of the root, a divisor of
    // `order`, so each divisor is reduced once.
    let mut cyclotomic = BTreeMap::new();
    let mut by_root_order: BTreeMap<usize, (Vec<BigInt>, Option<BigInt>)> = BTreeMap::new();
    for power in 0..order {
        let root_order = root_order_of_power(order, power);
        let (remainder, integer_value) = by_root_order
            .entry(root_order)
            .or_insert_with(|| evaluate_at_root_order(&coefficients, root_order, &mut cyclotomic))
            .clone();
        let evaluation = RootOfUnityEvaluation {
            group_order: order,
            power,
            root_order,
            integer_value,
            remainder,
        };
        if let Some(counts) = fixed_counts {
            let expected = counts[power].clone();
            let holds = evaluation.integer_value.as_ref() == Some(&expected);
            checks.push(CyclicSievingPowerCheck {
                power,
                root_order: evaluation.root_order,
                expected_fixed_points: expected,
                evaluation: evaluation.clone(),
                holds,
            });
        }
        evaluations.push(evaluation);
    }
    let holds = fixed_counts.map(|_| checks.iter().all(|check| check.holds));
    Ok(CyclicSievingReport {
        order,
        coefficients,
        fixed_counts: fixed_counts.map(|counts| counts.to_vec()),
        evaluations,
        checks,
        holds,
    })
}

/// Build sequence CSP profiles for candidate orders `index + offset`.
///
/// Row `r` of `polynomials` has sequence index `first_index + r`.  For each
/// `offset` the candidate group order is `index + offset`; non-positive orders
/// are skipped.  `fixed_counts` maps `(index, order)` to the fixed-point table
/// `fix(g^d)`, `d = 0..order-1`; a missing key gives a profile-only report.
/// A table of the wrong length makes that candidate order be omitted from
/// `candidate_orders` (no error is reported).
///
/// # Errors
///
/// The whole scan is validated before any report is built:
/// [`CyclicSievingError::SequenceIndexOverflow`] when `first_index + row` does
/// not fit in `isize` for some row, and [`CyclicSievingError::OrderTooLarge`]
/// when a positive candidate order `index + offset` (computed without overflow)
/// exceeds [`CYCLIC_SIEVING_MAX_ORDER`].
pub fn cyclic_sieving_sequence_reports_bigint(
    polynomials: &[Vec<BigInt>],
    first_index: isize,
    offsets: &[isize],
    fixed_counts: &BTreeMap<(isize, usize), Vec<BigInt>>,
) -> Result<Vec<CyclicSievingSequenceItem>, CyclicSievingError> {
    let mut indices = Vec::with_capacity(polynomials.len());
    for row in 0..polynomials.len() {
        indices.push(
            isize::try_from(row)
                .ok()
                .and_then(|row| first_index.checked_add(row))
                .ok_or(CyclicSievingError::SequenceIndexOverflow { first_index, row })?,
        );
    }
    for &index in &indices {
        for &offset in offsets {
            let order = index as i128 + offset as i128;
            if order > CYCLIC_SIEVING_MAX_ORDER as i128 {
                return Err(CyclicSievingError::OrderTooLarge {
                    order: order as u128,
                    maximum: CYCLIC_SIEVING_MAX_ORDER,
                });
            }
        }
    }
    Ok(polynomials
        .iter()
        .zip(indices)
        .enumerate()
        .map(|(row, (coeffs, index))| {
            let mut candidate_orders = Vec::new();
            for &offset in offsets {
                // `None` only for a negative overflow, i.e. a non-positive order.
                let Some(order_index) = index.checked_add(offset) else {
                    continue;
                };
                if order_index <= 0 {
                    continue;
                }
                let order = order_index as usize;
                let counts = fixed_counts.get(&(index, order)).map(Vec::as_slice);
                if let Ok(report) = cyclic_sieving_report_bigint(coeffs, order, counts) {
                    candidate_orders.push(report);
                }
            }
            let mut coefficients = coeffs.clone();
            strip_trailing_zeros(&mut coefficients);
            CyclicSievingSequenceItem {
                row,
                index,
                coefficients,
                candidate_orders,
            }
        })
        .collect())
}
