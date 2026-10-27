//! Exact verification of a recurrence against the input rows, including the
//! integer-scaled fast path used after fitting.

use super::*;
use num_integer::lcm;
use num_traits::One;
use num_traits::Zero;

#[derive(Debug)]
struct ScaledIntegerBivarPoly {
    coeffs: Vec<Vec<BigInt>>,
}

#[derive(Debug)]
struct ScaledIntegerRecurrenceTerm {
    offset: usize,
    deriv_order: usize,
    sign: RecurrenceSign,
    coeff: ScaledIntegerBivarPoly,
}

#[derive(Debug)]
struct ScaledIntegerRecurrence {
    terms: Vec<ScaledIntegerRecurrenceTerm>,
    denominator: Option<ScaledIntegerBivarPoly>,
    inhomogeneous: Option<ScaledIntegerBivarPoly>,
    scale: BigInt,
    n_degree: usize,
}

fn update_bivar_denominator_lcm(scale: &mut BigInt, poly: &BivarPoly) {
    for coeff in poly.coeffs.iter().flat_map(|row| row.iter()) {
        if !coeff.is_zero() && !coeff.denom().is_one() {
            *scale = lcm(scale.clone(), coeff.denom().clone());
        }
    }
}

fn recurrence_denominator_lcm(rec: &Recurrence) -> BigInt {
    let mut scale = BigInt::one();
    if let Some(denominator) = &rec.denominator {
        update_bivar_denominator_lcm(&mut scale, denominator);
    }
    for term in &rec.terms {
        update_bivar_denominator_lcm(&mut scale, &term.coeff);
    }
    if let Some(inhomogeneous) = &rec.inhomogeneous {
        update_bivar_denominator_lcm(&mut scale, inhomogeneous);
    }
    scale
}

fn rational_poly_has_integer_coeffs(poly: &[BigRational]) -> bool {
    poly.iter().all(|coeff| coeff.denom().is_one())
}

fn rational_polys_have_integer_coeffs(polys: &[Vec<BigRational>]) -> bool {
    polys
        .iter()
        .all(|poly| rational_poly_has_integer_coeffs(poly))
}

fn scale_bivar_poly_to_integer(poly: &BivarPoly, scale: &BigInt) -> ScaledIntegerBivarPoly {
    let coeffs = poly
        .coeffs
        .iter()
        .map(|row| {
            row.iter()
                .map(|coeff| {
                    if coeff.is_zero() {
                        BigInt::zero()
                    } else {
                        coeff.numer() * (scale / coeff.denom())
                    }
                })
                .collect()
        })
        .collect();
    ScaledIntegerBivarPoly { coeffs }
}

fn scale_recurrence_to_integer(rec: &Recurrence) -> ScaledIntegerRecurrence {
    let scale = recurrence_denominator_lcm(rec);
    let terms = rec
        .terms
        .iter()
        .map(|term| ScaledIntegerRecurrenceTerm {
            offset: term.offset,
            deriv_order: term.deriv_order,
            sign: term.sign,
            coeff: scale_bivar_poly_to_integer(&term.coeff, &scale),
        })
        .collect();
    let denominator = rec
        .denominator
        .as_ref()
        .map(|poly| scale_bivar_poly_to_integer(poly, &scale));
    let inhomogeneous = rec
        .inhomogeneous
        .as_ref()
        .map(|poly| scale_bivar_poly_to_integer(poly, &scale));
    ScaledIntegerRecurrence {
        terms,
        denominator,
        inhomogeneous,
        n_degree: recurrence_n_degree(rec),
        scale,
    }
}

fn bigint_index_powers(n: usize, max_power: usize) -> Vec<BigInt> {
    let n_big = BigInt::from(n);
    let mut powers = Vec::with_capacity(max_power + 1);
    let mut current = BigInt::one();
    for _ in 0..=max_power {
        powers.push(current.clone());
        current *= &n_big;
    }
    powers
}

fn update_poly_denominator_lcm(scale: &mut BigInt, poly: &[BigRational]) {
    for coeff in poly {
        if !coeff.is_zero() && !coeff.denom().is_one() {
            *scale = lcm(scale.clone(), coeff.denom().clone());
        }
    }
}

fn row_denominator_lcm_for_recurrence(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    rec: &ScaledIntegerRecurrence,
    nn: usize,
) -> Option<BigInt> {
    let mut scale = BigInt::one();
    update_poly_denominator_lcm(&mut scale, polys.get(nn - 1)?);
    for term in &rec.terms {
        if term.offset == 0 || term.offset >= nn {
            return None;
        }
        let deriv = derivs.get(nn - 1 - term.offset)?.get(term.deriv_order)?;
        update_poly_denominator_lcm(&mut scale, deriv);
    }
    Some(scale)
}

fn poly_add_rational_poly_scaled_assign(
    target: &mut Vec<BigInt>,
    source: &[BigRational],
    source_scale: &BigInt,
    coeff_scale: &BigInt,
    subtract: bool,
) {
    if source_scale.is_zero() || coeff_scale.is_zero() || poly_is_zero_rational(source) {
        return;
    }
    if target.len() < source.len() {
        target.resize(source.len(), BigInt::zero());
    }
    for (idx, coeff) in source.iter().enumerate() {
        if coeff.is_zero() {
            continue;
        }
        let value = if source_scale.is_one() && coeff.denom().is_one() {
            coeff.numer() * coeff_scale
        } else {
            coeff.numer() * (source_scale / coeff.denom()) * coeff_scale
        };
        if subtract {
            target[idx] -= value;
        } else {
            target[idx] += value;
        }
    }
}

fn poly_add_scaled_bivar_eval_signed_assign(
    target: &mut Vec<BigInt>,
    poly: &ScaledIntegerBivarPoly,
    n_powers: &[BigInt],
    source_scale: &BigInt,
    subtract: bool,
) {
    let width = poly.coeffs.iter().map(Vec::len).max().unwrap_or(0);
    if width == 0 {
        return;
    }
    if target.len() < width {
        target.resize(width, BigInt::zero());
    }
    for (i, row) in poly.coeffs.iter().enumerate() {
        let n_power = n_powers
            .get(i)
            .expect("bivariate evaluation requires enough powers of n");
        if n_power.is_zero() {
            continue;
        }
        for (j, coeff) in row.iter().enumerate() {
            if coeff.is_zero() {
                continue;
            }
            let value = if source_scale.is_one() {
                coeff * n_power
            } else {
                coeff * n_power * source_scale
            };
            if subtract {
                target[j] -= value;
            } else {
                target[j] += value;
            }
        }
    }
}

fn poly_add_scaled_bivar_eval_times_rational_poly_signed_assign(
    target: &mut Vec<BigInt>,
    coeff_poly: &ScaledIntegerBivarPoly,
    n_powers: &[BigInt],
    source: &[BigRational],
    source_scale: &BigInt,
    subtract: bool,
) {
    if poly_is_zero_rational(source) {
        return;
    }
    let width = coeff_poly.coeffs.iter().map(Vec::len).max().unwrap_or(0);
    if width == 0 {
        return;
    }
    let needed_len = width + source.len() - 1;
    if target.len() < needed_len {
        target.resize(needed_len, BigInt::zero());
    }

    for (i, row) in coeff_poly.coeffs.iter().enumerate() {
        let n_power = n_powers
            .get(i)
            .expect("bivariate evaluation requires enough powers of n");
        if n_power.is_zero() {
            continue;
        }
        for (j, coeff) in row.iter().enumerate() {
            if coeff.is_zero() {
                continue;
            }
            let scale = coeff * n_power;
            if scale.is_zero() {
                continue;
            }
            for (k, source_coeff) in source.iter().enumerate() {
                if source_coeff.is_zero() {
                    continue;
                }
                let value = if source_scale.is_one() && source_coeff.denom().is_one() {
                    &scale * source_coeff.numer()
                } else {
                    &scale * source_coeff.numer() * (source_scale / source_coeff.denom())
                };
                if subtract {
                    target[j + k] -= value;
                } else {
                    target[j + k] += value;
                }
            }
        }
    }
}

fn integer_residual_is_zero(residual: &[BigInt]) -> bool {
    residual.iter().all(|coeff| coeff.is_zero())
}

fn scaled_recurrence_holds_at_with_derivs_and_source_scale(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    rec: &ScaledIntegerRecurrence,
    nn: usize,
    source_scale: &BigInt,
) -> bool {
    if nn == 0 || nn > polys.len() || nn > derivs.len() {
        return false;
    }

    let n_powers = bigint_index_powers(nn, rec.n_degree);
    let current = &polys[nn - 1];
    let mut residual = vec![BigInt::zero()];
    if let Some(denom) = &rec.denominator {
        poly_add_scaled_bivar_eval_times_rational_poly_signed_assign(
            &mut residual,
            denom,
            &n_powers,
            current,
            source_scale,
            false,
        );
    } else {
        poly_add_rational_poly_scaled_assign(
            &mut residual,
            current,
            source_scale,
            &rec.scale,
            false,
        );
    }

    for term in &rec.terms {
        if term.offset == 0 || term.offset >= nn {
            return false;
        }
        let Some(poly_derivs) = derivs.get(nn - 1 - term.offset) else {
            return false;
        };
        let Some(deriv) = poly_derivs.get(term.deriv_order) else {
            return false;
        };
        let subtract = match term.sign {
            RecurrenceSign::None => true,
            RecurrenceSign::AlternatingN if nn % 2 == 1 => false,
            RecurrenceSign::AlternatingN => true,
        };
        poly_add_scaled_bivar_eval_times_rational_poly_signed_assign(
            &mut residual,
            &term.coeff,
            &n_powers,
            deriv,
            source_scale,
            subtract,
        );
    }

    if let Some(inhomogeneous) = &rec.inhomogeneous {
        poly_add_scaled_bivar_eval_signed_assign(
            &mut residual,
            inhomogeneous,
            &n_powers,
            source_scale,
            true,
        );
    }

    integer_residual_is_zero(&residual)
}

fn scaled_integer_recurrence_holds_at_with_derivs(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    rec: &ScaledIntegerRecurrence,
    nn: usize,
) -> bool {
    scaled_recurrence_holds_at_with_derivs_and_source_scale(polys, derivs, rec, nn, &BigInt::one())
}

fn scaled_rational_recurrence_holds_at_with_derivs(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    rec: &ScaledIntegerRecurrence,
    nn: usize,
) -> bool {
    let Some(source_scale) = row_denominator_lcm_for_recurrence(polys, derivs, rec, nn) else {
        return false;
    };
    scaled_recurrence_holds_at_with_derivs_and_source_scale(polys, derivs, rec, nn, &source_scale)
}

/// Check whether a recurrence holds at the 1-based row index `nn`.
///
/// This tests the exact polynomial identity
/// `f(n,t) P_n(t) = RHS(n,t)` without using the linear system that produced
/// the recurrence.
pub fn recurrence_holds_at_rational(
    polys: &[Vec<BigRational>],
    rec: &Recurrence,
    nn: usize,
) -> bool {
    recurrence_holds_at_rational_with_degree(polys, rec, nn, recurrence_n_degree(rec))
}

fn recurrence_holds_at_rational_with_degree(
    polys: &[Vec<BigRational>],
    rec: &Recurrence,
    nn: usize,
    rec_n_degree: usize,
) -> bool {
    if nn == 0 || nn > polys.len() {
        return false;
    }

    let n_powers = rational_index_powers(nn, rec_n_degree);
    let current = &polys[nn - 1];
    let mut residual = rational_zero_poly();
    if let Some(denom) = &rec.denominator {
        poly_add_bivar_eval_times_poly_signed_assign(
            &mut residual,
            denom,
            &n_powers,
            current,
            false,
        );
    } else {
        poly_add_signed_assign(&mut residual, current, false);
    }

    for term in &rec.terms {
        if term.offset == 0 || term.offset >= nn {
            return false;
        }
        let ref_poly = &polys[nn - 1 - term.offset];
        let deriv = poly_nth_derivative_rational(ref_poly, term.deriv_order);
        let subtract = match term.sign {
            RecurrenceSign::None => true,
            RecurrenceSign::AlternatingN if nn % 2 == 1 => false,
            RecurrenceSign::AlternatingN => true,
        };
        poly_add_bivar_eval_times_poly_signed_assign(
            &mut residual,
            &term.coeff,
            &n_powers,
            &deriv,
            subtract,
        );
    }

    if let Some(inhomogeneous) = &rec.inhomogeneous {
        poly_add_bivar_eval_signed_assign(&mut residual, inhomogeneous, &n_powers, true);
    }

    poly_is_zero_rational(&residual)
}

/// Check whether a recurrence holds for all admissible rows from `start_nn`.
///
/// The row index is 1-based. Rows before `rec_len + 1` are skipped because
/// there are not enough previous polynomials for the recurrence.
pub fn recurrence_holds_from_rational(
    polys: &[Vec<BigRational>],
    rec: &Recurrence,
    rec_len: usize,
    start_nn: usize,
) -> bool {
    let first = start_nn.max(rec_len + 1);
    let rec_n_degree = recurrence_n_degree(rec);
    (first..=polys.len())
        .all(|nn| recurrence_holds_at_rational_with_degree(polys, rec, nn, rec_n_degree))
}

pub(in crate::recurrence) fn recurrence_holds_from_rational_with_derivs(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    rec: &Recurrence,
    rec_len: usize,
    start_nn: usize,
) -> bool {
    let first = start_nn.max(rec_len + 1);
    let scaled = scale_recurrence_to_integer(rec);
    if rational_polys_have_integer_coeffs(polys) {
        return (first..=polys.len())
            .all(|nn| scaled_integer_recurrence_holds_at_with_derivs(polys, derivs, &scaled, nn));
    }

    (first..=polys.len())
        .all(|nn| scaled_rational_recurrence_holds_at_with_derivs(polys, derivs, &scaled, nn))
}

/// Check whether a recurrence holds on every admissible row.
pub fn recurrence_holds_rational(
    polys: &[Vec<BigRational>],
    rec: &Recurrence,
    rec_len: usize,
) -> bool {
    recurrence_holds_from_rational(polys, rec, rec_len, rec_len + 1)
}

pub(in crate::recurrence) fn recurrence_holds_rational_with_derivs(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    rec: &Recurrence,
    rec_len: usize,
) -> bool {
    recurrence_holds_from_rational_with_derivs(polys, derivs, rec, rec_len, rec_len + 1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::test_support::*;

    #[test]
    fn integer_scaled_verifier_handles_fractional_coefficients() {
        let polys = i64_polys_to_rational(&[vec![2], vec![2], vec![2], vec![2], vec![2]]);
        let derivs = rational_derivatives_up_to(&polys, 0);
        let rec = Recurrence {
            terms: vec![RecurrenceTerm::new(
                1,
                0,
                BivarPoly {
                    coeffs: vec![vec![br(1, 2)]],
                },
            )],
            denominator: None,
            inhomogeneous: Some(BivarPoly {
                coeffs: vec![vec![BigRational::one()]],
            }),
        };
        assert!(recurrence_holds_from_rational_with_derivs(
            &polys, &derivs, &rec, 1, 2
        ));

        let mut broken = polys.clone();
        broken[4][0] = BigRational::from_integer(BigInt::from(3));
        let broken_derivs = rational_derivatives_up_to(&broken, 0);
        assert!(!recurrence_holds_from_rational_with_derivs(
            &broken,
            &broken_derivs,
            &rec,
            1,
            2
        ));
    }
}
