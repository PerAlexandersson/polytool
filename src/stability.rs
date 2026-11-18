//! Exact univariate stability, re-exported from [`exact_stability`].
//!
//! [`check_univariate_stability`] returns inside/boundary/outside counts with
//! multiplicity for ascending Gaussian-rational coefficients. Upper/lower use
//! avoidance; left/right and disk/exterior use all-inside. The convention can
//! be overridden explicitly. The zero polynomial has no finite root counts.
//!
//! Real upper/lower input reduces to the integer Sturm kernel. For an
//! independent strict left-half-plane criterion on real integer coefficients,
//! use [`crate::hurwitz::is_hurwitz_stable_bigint`].
//!
//! ```
//! use polytool::stability::{check_univariate_stability, GaussianRational, StabilityRegion};
//! let p: Vec<_> = [1i64, 3, 3, 1].into_iter().map(GaussianRational::from).collect();
//! let report = check_univariate_stability(&p, &StabilityRegion::Left).unwrap();
//! assert!(report.verdict);
//! assert_eq!(report.counts.unwrap().inside, 3);
//! assert!(report.verify(&p));
//! ```

pub use exact_stability::*;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{complex_roots::certify_root_disc, hurwitz::is_hurwitz_stable_bigint};
    use num_bigint::BigInt;
    use num_traits::Zero;

    fn integer_polynomial(p: &[i64]) -> Vec<GaussianRational> {
        p.iter().copied().map(GaussianRational::from).collect()
    }

    #[test]
    fn stability_agrees_with_hurwitz_and_real_rootedness() {
        let mut seed = 0x98e2_35ab_7821u64;
        for degree in 0..=12 {
            for trial in 0..40 {
                let mut p: Vec<BigInt> = (0..=degree)
                    .map(|_| {
                        seed ^= seed << 13;
                        seed ^= seed >> 7;
                        seed ^= seed << 17;
                        BigInt::from((seed % 15) as i64 - 7)
                    })
                    .collect();
                p[degree] = BigInt::from(if trial % 2 == 0 { 1 } else { -1 });
                // Include positive-coefficient cases; Hurwitz's necessary
                // coefficient-sign screen must not be the only exercised path.
                if trial % 3 == 0 {
                    use num_traits::Signed;
                    for c in &mut p {
                        *c = c.abs() + 1;
                    }
                }
                let coefficients: Vec<_> = p.iter().cloned().map(GaussianRational::from).collect();
                let left =
                    check_univariate_stability(&coefficients, &StabilityRegion::Left).unwrap();
                assert_eq!(left.verdict, is_hurwitz_stable_bigint(&p), "{p:?}");
                let upper =
                    check_univariate_stability(&coefficients, &StabilityRegion::Upper).unwrap();
                assert_eq!(
                    upper.verdict,
                    crate::is_real_rooted_bigint_coeffs(&p),
                    "{p:?}"
                );
                assert_eq!(
                    upper.verdict,
                    crate::is_real_rooted_bezout_bigint_coeffs(&p),
                    "{p:?}"
                );
                let counts = upper.counts.unwrap();
                assert_eq!(counts.inside, counts.outside);
                assert_eq!(counts.inside + counts.boundary + counts.outside, degree);
            }
        }
        // The real-rootedness API accepts zero; the strict Hurwitz API does not.
        let zero = integer_polynomial(&[0, 0]);
        assert!(
            check_univariate_stability(&zero, &StabilityRegion::Upper)
                .unwrap()
                .verdict
        );
        assert!(
            !check_univariate_stability(&zero, &StabilityRegion::Left)
                .unwrap()
                .verdict
        );
    }

    #[test]
    fn stability_agrees_with_existing_bezoutian() {
        for f in [&[1, 3, 2][..], &[0, 0, 1, 0, 2], &[0], &[]] {
            for g in [&[2, -1][..], &[1, 0, 1, 2], &[0], &[]] {
                let f: Vec<BigInt> = f.iter().copied().map(BigInt::from).collect();
                let g: Vec<BigInt> = g.iter().copied().map(BigInt::from).collect();
                assert_eq!(
                    univariate::bezoutian::bezoutian_bigint_coeffs(&f, &g),
                    crate::bezoutian_bigint_coeffs(&f, &g)
                );
            }
        }
    }

    #[test]
    fn stability_eulerian_region_counts() {
        for p in crate::sequences::eulerian_polynomials_bigint(12) {
            let coefficients: Vec<_> = p.into_iter().map(GaussianRational::from).collect();
            assert!(
                check_univariate_stability(&coefficients, &StabilityRegion::Left)
                    .unwrap()
                    .verdict
            );
            assert!(
                check_univariate_stability(&coefficients, &StabilityRegion::Upper)
                    .unwrap()
                    .verdict
            );
            let disk = check_univariate_stability(&coefficients, &StabilityRegion::Disk).unwrap();
            let counts = disk.counts.unwrap();
            assert_eq!(counts.inside, counts.outside);
            assert_eq!(counts.boundary, (disk.degree.unwrap() % 2 == 1) as usize);
            assert_eq!(disk.verdict, disk.degree == Some(0));
        }
    }

    #[test]
    fn stability_agrees_with_root_disc_certificates() {
        // Six distinct rational conjugate pairs, of degree twelve. Each
        // disc is certified independently by the existing Taylor/Rouché path.
        let roots: Vec<_> = [(0, 2), (2, 1), (-2, 1), (0, 1), (3, 4), (-3, 4)]
            .into_iter()
            .flat_map(|(re, im)| {
                [
                    GaussianRational::new(
                        BigRational::from_integer(re.into()),
                        BigRational::from_integer(im.into()),
                    ),
                    GaussianRational::new(
                        BigRational::from_integer(re.into()),
                        BigRational::from_integer((-im).into()),
                    ),
                ]
            })
            .collect();
        let mut p = vec![GaussianRational::one()];
        for root in &roots {
            let mut next = vec![GaussianRational::zero(); p.len() + 1];
            for (k, coefficient) in p.iter().enumerate() {
                next[k] = next[k].subtract(&coefficient.multiply(root));
                next[k + 1] = next[k + 1].add(coefficient);
            }
            p = next;
        }
        assert!(p.iter().all(|x| x.imaginary.is_zero()));
        let real: Vec<_> = p.iter().map(|x| x.real.clone()).collect();
        for root in &roots {
            let certificate = certify_root_disc(
                &real,
                &root.real,
                &root.imaginary,
                &BigRational::new(1.into(), 100.into()),
            )
            .expect("simple exact root disc");
            assert!(certificate.verify(&real));
        }
        for region in [
            StabilityRegion::Upper,
            StabilityRegion::Lower,
            StabilityRegion::Left,
            StabilityRegion::Right,
            StabilityRegion::Disk,
            StabilityRegion::Exterior,
        ] {
            let mut expected = RootCounts::default();
            for root in &roots {
                match region.membership_sign(root).unwrap() {
                    1 => expected.inside += 1,
                    0 => expected.boundary += 1,
                    _ => expected.outside += 1,
                }
            }
            assert_eq!(
                check_univariate_stability(&p, &region).unwrap().counts,
                Some(expected)
            );
        }
    }
}
