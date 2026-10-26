//! Conversions between h*-vectors and Ehrhart polynomials.
//!
//! Ehrhart coefficients are exact rationals in ascending degree order;
//! [`EhrhartToHstarError`] reports inputs that do not give an integral h*-vector.

use super::Q;
use num_bigint::BigInt;

/// Convert an h\*-vector to an Ehrhart polynomial.
///
/// Given h\* = (h\*\_0, ..., h\*\_d), the Ehrhart polynomial is
///
/// ```text
/// L(t) = sum_{i=0}^{d} h*_i * C(t + d - i, d)
/// ```
///
/// where C(t+a, d) = (t+a)(t+a-1)...(t+a-d+1)/d! is a polynomial of degree d in t.
///
/// Returns the polynomial coefficients as rationals in ascending degree order.
/// For a valid h\*-vector of a lattice polytope, L(t) has degree d = len(hstar) - 1.
pub fn hstar_to_ehrhart(hstar: &[i64]) -> Vec<Q> {
    let hstar = hstar.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
    hstar_to_ehrhart_bigint_coeffs(&hstar)
}

/// Convert an arbitrary-size integer h\*-vector to an Ehrhart polynomial.
pub fn hstar_to_ehrhart_bigint_coeffs(hstar: &[BigInt]) -> Vec<Q> {
    if hstar.is_empty() {
        return vec![];
    }
    let d = hstar.len() - 1;
    if d == 0 {
        return vec![Q::from_integer(hstar[0].clone())];
    }

    let zero = Q::from_integer(BigInt::from(0));
    let one = Q::from_integer(BigInt::from(1));
    let mut result = vec![zero.clone(); d + 1];

    // Compute d! for normalization
    let mut d_fact = Q::from_integer(BigInt::from(1));
    for j in 2..=d {
        d_fact *= Q::from_integer(BigInt::from(j as i64));
    }

    for (i, hi) in hstar.iter().enumerate() {
        if hi == &BigInt::from(0) {
            continue;
        }
        // Build C(t + d - i, d) as a polynomial in t.
        // C(t + a, d) = product_{j=0}^{d-1} (t + a - j) / d!  where a = d - i.
        let a = (d - i) as i64;

        // Multiply linear factors (t + a)(t + a - 1)...(t + a - d + 1)
        // Start with polynomial = 1
        let mut poly = vec![one.clone()]; // degree 0
        for j in 0..d {
            let shift = a - j as i64; // constant term of the linear factor (t + shift)
            let shift_q = Q::from_integer(BigInt::from(shift));
            // Multiply poly by (t + shift): new[k] = poly[k-1] + shift * poly[k]
            let mut new_poly = vec![zero.clone(); poly.len() + 1];
            for (k, c) in poly.iter().enumerate() {
                new_poly[k] = new_poly[k].clone() + shift_q.clone() * c.clone();
                new_poly[k + 1] = new_poly[k + 1].clone() + c.clone();
            }
            poly = new_poly;
        }

        // Divide by d! and scale by h*_i
        let scale = Q::from_integer(hi.clone()) / d_fact.clone();
        for (k, c) in poly.iter().enumerate() {
            if k <= d {
                result[k] = result[k].clone() + scale.clone() * c.clone();
            }
        }
    }

    // Trim trailing zeros
    while result.last().is_some_and(|c| *c == zero) {
        result.pop();
    }
    result
}

/// Convert an Ehrhart polynomial to an h\*-vector.
///
/// Given L(t) of degree d with rational coefficients, the h\*-vector entries are
///
/// ```text
/// h*_i = sum_{k=0}^{i} (-1)^k * C(d+1, k) * L(i-k)    for i = 0, ..., d
/// ```
///
/// where L(j) is the Ehrhart polynomial evaluated at the non-negative integer j.
///
/// The h\*-vector entries are always integers for Ehrhart polynomials of lattice
/// polytopes. Nonintegral results and entries outside `i64` are reported as
/// errors rather than rounded or panicked on.
pub fn ehrhart_to_hstar(ehrhart_coeffs: &[Q]) -> Result<Vec<i64>, EhrhartToHstarError> {
    ehrhart_to_hstar_bigint(ehrhart_coeffs)?
        .into_iter()
        .enumerate()
        .map(|(index, value)| {
            i64::try_from(&value)
                .map_err(|_| EhrhartToHstarError::EntryOutOfI64Range { index, value })
        })
        .collect()
}

/// Errors returned when exact Ehrhart-to-h\* conversion is not defined for the
/// supplied representation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EhrhartToHstarError {
    /// The common denominator passed to [`ehrhart_to_hstar_with_denom`] is zero.
    ZeroDenominator,
    /// An h\*-entry is not an integer (so the input is not the Ehrhart
    /// polynomial of a lattice polytope).
    NonIntegralEntry {
        /// Index `i` of the entry `h*_i`.
        index: usize,
        /// The nonintegral value.
        value: Q,
    },
    /// An h\*-entry is an integer but does not fit in `i64`.
    EntryOutOfI64Range {
        /// Index `i` of the entry `h*_i`.
        index: usize,
        /// The entry.
        value: BigInt,
    },
}

impl std::fmt::Display for EhrhartToHstarError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ZeroDenominator => write!(f, "Ehrhart coefficient denominator must be nonzero"),
            Self::NonIntegralEntry { index, value } => {
                write!(f, "h*-vector entry {index} is not an integer: {value}")
            }
            Self::EntryOutOfI64Range { index, value } => {
                write!(f, "h*-vector entry {index} does not fit in i64: {value}")
            }
        }
    }
}

impl std::error::Error for EhrhartToHstarError {}

/// Convert an Ehrhart polynomial to an arbitrary-size integer h\*-vector.
pub fn ehrhart_to_hstar_bigint(ehrhart_coeffs: &[Q]) -> Result<Vec<BigInt>, EhrhartToHstarError> {
    if ehrhart_coeffs.is_empty() {
        return Ok(vec![]);
    }
    let zero = Q::from_integer(BigInt::from(0));
    // Find degree
    let d = ehrhart_coeffs.iter().rposition(|c| *c != zero).unwrap_or(0);

    // Evaluate L(t) at t = 0, 1, ..., d
    let mut values = Vec::with_capacity(d + 1);
    for t in 0..=d {
        let t_q = Q::from_integer(BigInt::from(t as i64));
        let mut val = zero.clone();
        let mut t_pow = Q::from_integer(BigInt::from(1));
        for c in ehrhart_coeffs.iter() {
            val += c.clone() * t_pow.clone();
            t_pow *= t_q.clone();
        }
        values.push(val);
    }

    // Binomial coefficients C(d+1, k) for k = 0, ..., d+1
    let mut binom = vec![BigInt::from(1); d + 2];
    for k in 1..=d + 1 {
        binom[k] = &binom[k - 1] * BigInt::from(d + 1 - k + 1) / BigInt::from(k);
    }

    // h*_i = sum_{k=0}^{i} (-1)^k * C(d+1, k) * L(i-k)
    let mut hstar = Vec::with_capacity(d + 1);
    for i in 0..=d {
        let mut val = zero.clone();
        for k in 0..=i {
            let signed_binom = if k % 2 == 0 {
                binom[k].clone()
            } else {
                -binom[k].clone()
            };
            val += Q::from_integer(signed_binom) * values[i - k].clone();
        }
        if val.denom() != &BigInt::from(1) {
            return Err(EhrhartToHstarError::NonIntegralEntry {
                index: i,
                value: val,
            });
        }
        hstar.push(val.to_integer());
    }

    Ok(hstar)
}

/// Convenience: convert Ehrhart polynomial given as `i64` coefficients of the
/// **numerator** polynomial when written over a common denominator.
///
/// That is, if L(t) = (a_0 + a_1 t + ... + a_d t^d) / denom, pass `(numerator_coeffs, denom)`.
/// This is useful when the Ehrhart polynomial is stored in integer-denominator form.
pub fn ehrhart_to_hstar_with_denom(
    numerator_coeffs: &[i64],
    denom: i64,
) -> Result<Vec<i64>, EhrhartToHstarError> {
    if denom == 0 {
        return Err(EhrhartToHstarError::ZeroDenominator);
    }
    let denom_q = Q::from_integer(BigInt::from(denom));
    let coeffs: Vec<Q> = numerator_coeffs
        .iter()
        .map(|&c| Q::from_integer(BigInt::from(c)) / denom_q.clone())
        .collect();
    ehrhart_to_hstar(&coeffs)
}

#[cfg(test)]
mod tests {
    use super::*;

    // -- Ehrhart / h*-vector tests --

    fn q(n: i64, d: i64) -> Q {
        Q::new(BigInt::from(n), BigInt::from(d))
    }

    #[test]
    fn test_hstar_to_ehrhart_simplex() {
        // Standard d-simplex has h* = [1, 0, ..., 0].
        // L(t) = C(t+d, d).

        // 1-simplex [0,1]: h* = [1, 0], L(t) = t + 1
        let ehrhart = hstar_to_ehrhart(&[1, 0]);
        assert_eq!(ehrhart.len(), 2);
        assert_eq!(ehrhart[0], q(1, 1)); // constant term = 1
        assert_eq!(ehrhart[1], q(1, 1)); // t coefficient = 1

        // 2-simplex: h* = [1, 0, 0], L(t) = C(t+2,2) = (t+1)(t+2)/2 = 1 + 3t/2 + t^2/2
        let ehrhart = hstar_to_ehrhart(&[1, 0, 0]);
        assert_eq!(ehrhart.len(), 3);
        assert_eq!(ehrhart[0], q(1, 1)); // 1
        assert_eq!(ehrhart[1], q(3, 2)); // 3/2
        assert_eq!(ehrhart[2], q(1, 2)); // 1/2
    }

    #[test]
    fn test_hstar_to_ehrhart_square() {
        // Unit square: h* = [1, 1], L(t) = (t+1)^2 = 1 + 2t + t^2
        let ehrhart = hstar_to_ehrhart(&[1, 1]);
        assert_eq!(ehrhart.len(), 2);
        // Wait, h* = [1, 1] means d = 1, but the unit square is dimension 2...
        // Actually h* has d+1 entries for a d-dimensional polytope.
        // Unit square: d=2, h* = [1, 1, 0]
        let ehrhart = hstar_to_ehrhart(&[1, 1, 0]);
        assert_eq!(ehrhart.len(), 3);
        assert_eq!(ehrhart[0], q(1, 1));
        assert_eq!(ehrhart[1], q(2, 1));
        assert_eq!(ehrhart[2], q(1, 1));
    }

    #[test]
    fn test_ehrhart_to_hstar_simplex() {
        // 1-simplex: L(t) = t + 1 = [1, 1], h* = [1, 0]
        let ehrhart = vec![q(1, 1), q(1, 1)];
        assert_eq!(ehrhart_to_hstar(&ehrhart).unwrap(), vec![1, 0]);

        // 2-simplex: L(t) = 1 + 3t/2 + t^2/2, h* = [1, 0, 0]
        let ehrhart = vec![q(1, 1), q(3, 2), q(1, 2)];
        assert_eq!(ehrhart_to_hstar(&ehrhart).unwrap(), vec![1, 0, 0]);
    }

    #[test]
    fn test_ehrhart_to_hstar_square() {
        // Unit square: L(t) = (t+1)^2 = 1 + 2t + t^2, h* = [1, 1, 0]
        let ehrhart = vec![q(1, 1), q(2, 1), q(1, 1)];
        assert_eq!(ehrhart_to_hstar(&ehrhart).unwrap(), vec![1, 1, 0]);
    }

    #[test]
    fn test_roundtrip_hstar_ehrhart() {
        // h* -> Ehrhart -> h* should be identity
        let hstar = vec![1, 4, 6, 4, 1];
        let ehrhart = hstar_to_ehrhart(&hstar);
        let recovered = ehrhart_to_hstar(&ehrhart).unwrap();
        assert_eq!(recovered, hstar);
    }

    #[test]
    fn test_ehrhart_to_hstar_with_denom() {
        // 2-simplex: L(t) = (t^2 + 3t + 2) / 2
        // numerator = [2, 3, 1], denom = 2
        assert_eq!(
            ehrhart_to_hstar_with_denom(&[2, 3, 1], 2).unwrap(),
            vec![1, 0, 0]
        );
    }

    #[test]
    fn test_ehrhart_to_hstar_reports_invalid_inputs() {
        assert_eq!(
            ehrhart_to_hstar_bigint(&[q(1, 2)]),
            Err(EhrhartToHstarError::NonIntegralEntry {
                index: 0,
                value: q(1, 2),
            })
        );
        assert_eq!(
            ehrhart_to_hstar_with_denom(&[1], 0),
            Err(EhrhartToHstarError::ZeroDenominator)
        );
    }
}
