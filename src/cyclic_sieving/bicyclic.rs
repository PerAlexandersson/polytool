//! Exact bicyclic sieving reports.

use super::*;

/// Largest `n1 * n2` accepted by [`bicyclic_sieving_report_bigint`].
pub const BICYCLIC_MAX_GROUP_SIZE: usize = 2500;

/// Bicyclic sieving report for `G = C_{n1} x C_{n2}` and `f(q, t)`.
///
/// `f[i][j]` is the coefficient of `q^i t^j`.  The polynomial is reduced modulo
/// `(q^{n1} - 1, t^{n2} - 1)`.  The values `f(zeta_{n1}^{d1}, zeta_{n2}^{d2})` are
/// computed exactly (as one reduction modulo the cyclotomic polynomial
/// `Phi_L`, `L = lcm(n1, n2)`).
///
/// If `c_{j1,j2} = c_{gcd(j1,n1), gcd(j2,n2)}` for all `j1, j2` (the polynomial
/// lies in the span of the products `h_{r1}(q) h_{r2}(t)`, the product basis
/// `B_2` of the Alexandersson--Amini multivariate basis theorem), the
/// coordinates are
/// `a_{r1,r2} = sum_{e1 | r1, e2 | r2} mu(r1/e1) mu(r2/e2) c_{e1,e2}`, and
/// `f(zeta^{d1}, omega^{d2}) = sum_{(n1/r1) | d1, (n2/r2) | d2} (n1/r1)(n2/r2) a_{r1,r2}`.
///
/// These product coordinates are **not** orbit counts and need not be
/// nonnegative: for the diagonal action of `C_2 x C_2` on a 2-element set,
/// `f = 1 + q t` and the coordinates `a_{1,1}, a_{1,2}, a_{2,1}, a_{2,2}` are
/// `1, -1, -1, 2`.  Products `h_{r1}(q) h_{r2}(t)` are the CSP polynomials of
/// product orbits `C_{n1}/H1 x C_{n2}/H2`; a general `G`-orbit `G/H` need not be
/// of this form.
///
/// Integrality at all points is **not** equivalent to product-basis form when
/// `gcd(n1, n2) > 1`: for `n1 = n2 = 3`,
/// `f = sum_{i,j} (zeta_3^{i+j} + zeta_3^{-(i+j)}) q^i t^j` has coefficients
/// `2` if `i + j = 0 (mod 3)` and `-1` otherwise, takes only the values `0` and
/// `9` on `mu_3 x mu_3`, yet is not in the product span.  (For coprime `n1, n2`
/// the two notions agree, because the Galois groups act independently.)  We
/// therefore report `integral` (from exact evaluation) and `product_basis_form`
/// separately; the CSP check against a fixed-point table uses the exact
/// evaluations only.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BicyclicSievingReport {
    /// Order `n1` of the first cyclic factor (acting on `q`).
    pub order1: usize,
    /// Order `n2` of the second cyclic factor (acting on `t`).
    pub order2: usize,
    /// The input coefficient matrix `f[i][j]` (coefficient of `q^i t^j`); rows
    /// may be ragged and are not trimmed.
    pub coefficient_matrix: Vec<Vec<BigInt>>,
    /// `f mod (q^{n1} - 1, t^{n2} - 1)`, an `n1 x n2` matrix.
    pub reduced_matrix: Vec<Vec<BigInt>>,
    /// `f(zeta_{n1}^{d1}, zeta_{n2}^{d2})` for `d1 < n1`, `d2 < n2`; `None`
    /// where the value is not an integer.
    pub evaluations: Vec<Vec<Option<BigInt>>>,
    /// Every evaluation is an integer.
    pub integral: bool,
    /// All input coefficients are nonnegative.
    pub natural_coefficients: bool,
    /// Divisors `r1` of `n1`, ascending: the row index of `product_coordinates`
    /// (the row for `r1` belongs to `h_{r1}(q)`, a polynomial of one orbit of
    /// size `n1 / r1`).
    pub divisors1: Vec<usize>,
    /// Divisors `r2` of `n2`, ascending: the column index of `product_coordinates`.
    pub divisors2: Vec<usize>,
    /// Whether `c_{j1,j2} = c_{gcd(j1,n1), gcd(j2,n2)}` for all `j1, j2`.
    pub product_basis_form: bool,
    /// The coordinates `a_{r1,r2}` in the basis `h_{r1}(q) h_{r2}(t)`, rows
    /// indexed by `r1 in divisors1` and columns by `r2 in divisors2`; `Some` iff
    /// `product_basis_form`.
    pub product_coordinates: Option<Vec<Vec<BigInt>>>,
    /// `f(zeta^{d1}, omega^{d2})` recomputed from `product_coordinates`.
    pub product_formula_evaluations: Option<Vec<Vec<BigInt>>>,
    /// The supplied table `fix(g1^{d1} g2^{d2})`.
    pub fixed_table: Option<Vec<Vec<BigInt>>>,
    /// Positions `(d1, d2)` where the evaluation is not the supplied count.
    pub mismatches: Vec<(usize, usize)>,
    /// `Some(true)` iff a table was supplied and equals the evaluations.
    pub holds: Option<bool>,
}

/// Compute a bicyclic sieving report; see [`BicyclicSievingReport`].
///
/// `fixed_table`, when given, must be an `order1 x order2` matrix whose entry
/// `[d1][d2]` is `fix(g1^d1 g2^d2)`.  Returns `Err` when an order is zero, when
/// `order1 * order2` exceeds [`BICYCLIC_MAX_GROUP_SIZE`], or when the table has
/// the wrong shape.
pub fn bicyclic_sieving_report_bigint(
    coeff_matrix: &[Vec<BigInt>],
    order1: usize,
    order2: usize,
    fixed_table: Option<&[Vec<BigInt>]>,
) -> Result<BicyclicSievingReport, CyclicSievingError> {
    if order1 == 0 || order2 == 0 {
        return Err(CyclicSievingError::NonpositiveBicyclicOrders);
    }
    if order1
        .checked_mul(order2)
        .is_none_or(|g| g > BICYCLIC_MAX_GROUP_SIZE)
    {
        return Err(CyclicSievingError::BicyclicGroupTooLarge { order1, order2 });
    }
    if let Some(table) = fixed_table {
        if table.len() != order1 || table.iter().any(|row| row.len() != order2) {
            return Err(CyclicSievingError::BicyclicFixedTableShape { order1, order2 });
        }
    }
    let (n1, n2) = (order1, order2);
    let mut reduced = vec![vec![BigInt::zero(); n2]; n1];
    for (i, row) in coeff_matrix.iter().enumerate() {
        for (j, c) in row.iter().enumerate() {
            reduced[i % n1][j % n2] += c;
        }
    }
    let natural_coefficients = coeff_matrix
        .iter()
        .flatten()
        .all(|c| c.sign() != num_bigint::Sign::Minus);

    // Exact evaluations: f(zeta_{n1}^{d1}, zeta_{n2}^{d2}) = R(xi), xi = zeta_L,
    // R = sum c_{ij} x^{(i d1 L/n1 + j d2 L/n2) mod L}.
    let l = n1.lcm(&n2);
    let phi_l = cyclotomic_polynomial_bigint(l);
    let mut evaluations = vec![vec![None; n2]; n1];
    for (d1, eval_row) in evaluations.iter_mut().enumerate() {
        for (d2, slot) in eval_row.iter_mut().enumerate() {
            let (u1, u2) = (d1 * (l / n1), d2 * (l / n2));
            let mut poly = vec![BigInt::zero(); l];
            for (i, row) in reduced.iter().enumerate() {
                for (j, c) in row.iter().enumerate() {
                    if !c.is_zero() {
                        poly[(i * u1 + j * u2) % l] += c;
                    }
                }
            }
            *slot = match remainder_monic(poly, &phi_l).as_slice() {
                [] => Some(BigInt::zero()),
                [constant] => Some(constant.clone()),
                _ => None,
            };
        }
    }
    let integral = evaluations.iter().flatten().all(Option::is_some);

    let divisors1 = divisors(n1);
    let divisors2 = divisors(n2);
    let product_basis_form = (0..n1).all(|j1| {
        (0..n2).all(|j2| reduced[j1][j2] == reduced[gcd_class(j1, n1) % n1][gcd_class(j2, n2) % n2])
    });
    let (product_coordinates, product_formula_evaluations) = if product_basis_form {
        let c_at = |e1: usize, e2: usize| &reduced[e1 % n1][e2 % n2];
        let coords: Vec<Vec<BigInt>> = divisors1
            .iter()
            .map(|&r1| {
                divisors2
                    .iter()
                    .map(|&r2| {
                        let mut a = BigInt::zero();
                        for &e1 in divisors1.iter().filter(|&&e| r1 % e == 0) {
                            for &e2 in divisors2.iter().filter(|&&e| r2 % e == 0) {
                                let mu = mobius(r1 / e1) * mobius(r2 / e2);
                                a += big_mobius_term(mu, c_at(e1, e2));
                            }
                        }
                        a
                    })
                    .collect()
            })
            .collect();
        let formula: Vec<Vec<BigInt>> = (0..n1)
            .map(|d1| {
                (0..n2)
                    .map(|d2| {
                        let mut value = BigInt::zero();
                        for (x1, &r1) in divisors1.iter().enumerate() {
                            for (x2, &r2) in divisors2.iter().enumerate() {
                                // (n1/r1) must divide d1 and (n2/r2) must divide d2
                                let (k1, k2) = (n1 / r1, n2 / r2);
                                if d1 % k1 == 0 && d2 % k2 == 0 {
                                    value += &coords[x1][x2] * BigInt::from(k1 * k2);
                                }
                            }
                        }
                        value
                    })
                    .collect()
            })
            .collect();
        (Some(coords), Some(formula))
    } else {
        (None, None)
    };

    let mut mismatches = Vec::new();
    if let Some(table) = fixed_table {
        for (d1, row) in table.iter().enumerate() {
            for (d2, expected) in row.iter().enumerate() {
                if evaluations[d1][d2].as_ref() != Some(expected) {
                    mismatches.push((d1, d2));
                }
            }
        }
    }
    let holds = fixed_table.map(|_| mismatches.is_empty());
    Ok(BicyclicSievingReport {
        order1,
        order2,
        coefficient_matrix: coeff_matrix.to_vec(),
        reduced_matrix: reduced,
        evaluations,
        integral,
        natural_coefficients,
        divisors1,
        divisors2,
        product_basis_form,
        product_coordinates,
        product_formula_evaluations,
        fixed_table: fixed_table.map(|t| t.to_vec()),
        mismatches,
        holds,
    })
}
