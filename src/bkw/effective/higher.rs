//! Reuse exact Q[t][lambda] factorization for automatic higher-order pairs.

use super::*;
use crate::factorization::irreducible_factorization_monic_bivariate_with_budget;
use crate::recurrence::generating_function::exact_divide;

pub(super) fn polynomial_pairs(
    d: &[P],
) -> Result<Vec<EffectiveBkwPairFactorization>, RecurrenceError> {
    let characteristic: Vec<_> = d.iter().rev().map(|p| p.coeffs().to_vec()).collect();
    let factors = irreducible_factorization_monic_bivariate_with_budget(&characteristic, 4096)
        .map_err(|e| failure(format!("automatic characteristic factorization: {e}")))?;
    let denominator = B::new(d.iter().map(|p| p.coeffs().to_vec()).collect());
    let mut pairs = Vec::new();
    for factor in factors {
        if factor.multiplicity != 1 || factor.coefficients.len() != 3 {
            continue;
        }
        let b = P::new(factor.coefficients[0].clone());
        let a = P::new(factor.coefficients[1].clone()).scale(&(q(-1) / q(2)));
        let pair = B::new(vec![
            vec![Q::one()],
            a.scale(&q(-2)).coeffs().to_vec(),
            b.coeffs().to_vec(),
        ]);
        if let Ok(h) = exact_divide(&denominator, &pair) {
            pairs.push(EffectiveBkwPairFactorization {
                real_part: encode(&a),
                norm_squared: encode(&b),
                cofactor: encode_b(&h),
            });
        }
    }
    Ok(pairs)
}

/// Q(b)=prod_(i<j)(b-lambda_i lambda_j). Newton identities give it directly
/// from p_k=sum lambda_i^k and w_k=(p_k^2-p_(2k))/2. This is the square root
/// of the off-diagonal part of Res_x(p(x),x^d p(b/x)), with no root guessing.
fn pair_products(d: &[P]) -> B {
    let order = d.len() - 1;
    let size = order * (order - 1) / 2;
    let mut powers = vec![P::new(vec![Q::from_integer(order.into())])];
    for k in 1..=2 * size {
        let mut sum = P::zero();
        for j in 1..=order.min(k) {
            sum = sum
                + d[j].clone()
                    * if j == k {
                        P::new(vec![Q::from_integer(k.into())])
                    } else {
                        powers[k - j].clone()
                    };
        }
        powers.push(-sum);
    }
    let pair_powers: Vec<_> = (1..=size)
        .map(|k| (powers[k].pow(2) - powers[2 * k].clone()).scale(&(q(1) / q(2))))
        .collect();
    let mut coefficients = vec![P::one()];
    for k in 1..=size {
        let mut sum = P::zero();
        for j in 1..=k {
            sum = sum + coefficients[k - j].clone() * pair_powers[j - 1].clone();
        }
        coefficients.push(sum.scale(&(-Q::from_integer(k.into()).recip())));
    }
    coefficients.reverse();
    B::new(coefficients.iter().map(|p| p.coeffs().to_vec()).collect())
}

pub(super) fn algebraic_pair(d: &[P]) -> Result<(B, R, R, Vec<R>, Branch), RecurrenceError> {
    let order = d.len() - 1;
    let p = pair_products(d);
    // Q is already monic. Its discriminant is checked before acceptance;
    // repeated pair products are unsupported, so no bivariate gcd is needed.
    #[cfg(test)]
    if std::env::var_os("POLYTOOL_EFFECTIVE_BKW_TRACE").is_some() {
        eprintln!("higher: pair-product polynomial");
    }
    // The degree-two subresultant of p(x) and x^d p(b/x), using three
    // minors of the truncated Sylvester matrix. Shared Bareiss elimination
    // keeps degrees bounded, unlike a repeated naive pseudo-remainder.
    let left: Vec<B> = d.iter().rev().map(t).collect();
    let right: Vec<B> = (0..=order)
        .map(|j| t(&d[j]).multiply(&variable().power(order - j)))
        .collect();
    let size = 2 * order - 4;
    let columns = size + 2;
    let mut matrix = vec![vec![B::zero(); columns]; size];
    for (half, coefficients) in [&left, &right].into_iter().enumerate() {
        for i in 0..order - 2 {
            for j in 0..=order {
                matrix[half * (order - 2) + i][i + j] = coefficients[order - j].clone();
            }
        }
    }
    let subresultant: Vec<_> = (0..=2)
        .map(|j| {
            let minor: Vec<_> = matrix
                .iter()
                .map(|row| {
                    row.iter()
                        .take(size - 1)
                        .cloned()
                        .chain(std::iter::once(row[columns - 1 - j].clone()))
                        .collect()
                })
                .collect();
            rem(
                &super::super::certified::algebra::fraction_free_determinant(minor),
                &p,
            )
        })
        .collect();
    #[cfg(test)]
    if std::env::var_os("POLYTOOL_EFFECTIVE_BKW_TRACE").is_some() {
        eprintln!("higher: subresultant minors");
    }
    if subresultant[2].is_zero()
        || !rem(
            &subresultant[0].subtract(&variable().multiply(&subresultant[2])),
            &p,
        )
        .is_zero()
    {
        return Err(failure("higher-order quadratic subresultant degenerates"));
    }
    let a = R::new(subresultant[1].scale(&q(-1)), subresultant[2].scale(&q(2))).modulo(&p);
    #[cfg(test)]
    if std::env::var_os("POLYTOOL_EFFECTIVE_BKW_TRACE").is_some() {
        eprintln!("higher: pair real part");
    }
    let b = R::of(variable()).modulo(&p);
    let mut h = vec![R::of(c(1))];
    for j in 1..=order - 2 {
        let mut value = R::of(t(&d[j])).add(&a.scale(&q(2)).mul(&h[j - 1]));
        if j >= 2 {
            value = value.sub(&b.mul(&h[j - 2]));
        }
        h.push(value);
        #[cfg(test)]
        if std::env::var_os("POLYTOOL_EFFECTIVE_BKW_TRACE").is_some() {
            eprintln!("higher: cofactor {j}");
        }
    }
    Ok((p, a, b, h, Branch::LargestPositive))
}
