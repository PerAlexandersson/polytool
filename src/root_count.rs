//! Compatibility facade for Polytool's exact integer root-counting API.
//!
//! The reusable implementation lives in the independent
//! [`exact_real_roots`] crate. Existing `polytool::root_count` imports remain
//! valid through these re-exports.

pub use exact_real_roots::*;

#[cfg(test)]
mod tests {
    use super::is_real_rooted_prs_bigint_coeffs;
    use num_bigint::BigInt;

    #[test]
    fn prs_agrees_with_polytool_rational_sturm_on_small_grid() {
        for degree in 0usize..=4 {
            let total = 5usize.pow((degree + 1) as u32);
            for mut mask in 0..total {
                let mut coeffs = Vec::with_capacity(degree + 1);
                for _ in 0..=degree {
                    coeffs.push((mask % 5) as i64 - 2);
                    mask /= 5;
                }

                let bigint: Vec<BigInt> = coeffs.iter().copied().map(BigInt::from).collect();
                assert_eq!(
                    is_real_rooted_prs_bigint_coeffs(&bigint),
                    crate::real_rootedness::is_real_rooted_sturm_bigint_coeffs(&bigint),
                    "coeffs={coeffs:?}"
                );
            }
        }
    }
}
