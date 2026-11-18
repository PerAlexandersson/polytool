//! Exact region, boundary, multiplicity and independent Schur--Cohn checks.

use exact_algebra::inertia::symmetric_inertia;
use exact_stability::{
    check_univariate_stability as check, check_univariate_stability_with_convention,
    BigRational as Q, GaussianRational as G, MobiusTransform, RootCounts, StabilityConvention,
    StabilityError, StabilityMethod, StabilityRegion as R,
};
use num_traits::Zero;

fn q(n: i64, d: i64) -> Q {
    Q::new(n.into(), d.into())
}
fn g(real: i64, imaginary: i64) -> G {
    G::new(q(real, 1), q(imaginary, 1))
}

fn from_roots(roots: &[G]) -> Vec<G> {
    let mut p = vec![G::one()];
    for root in roots {
        let mut next = vec![G::zero(); p.len() + 1];
        for (i, coefficient) in p.iter().enumerate() {
            next[i] = next[i].subtract(&coefficient.multiply(root));
            next[i + 1] = next[i + 1].add(coefficient);
        }
        p = next;
    }
    p
}

fn counts(p: &[G], region: &R) -> RootCounts {
    check(p, region).unwrap().counts.unwrap()
}

#[test]
fn gaussian_field_arithmetic() {
    let a = G::new(q(2, 3), q(-5, 7));
    let b = G::new(q(-3, 11), q(4, 5));
    assert_eq!(a.multiply(&b).checked_divide(&b), Some(a.clone()));
    assert_eq!(a.checked_divide(&G::zero()), None);
    assert_eq!(a.add(&b).subtract(&b), a);
    assert_eq!(a.multiply(&a.conjugate()), G::from(a.norm_squared()));
}

#[test]
fn exact_near_boundary_and_repeated_nonreal_roots() {
    let epsilon = Q::new(1.into(), num_bigint::BigInt::from(1) << 160);
    let roots = [
        G::from(Q::from_integer(1.into()) - &epsilon),
        G::from(Q::from_integer(1.into()) + &epsilon),
        g(1, 0),
    ];
    assert_eq!(
        counts(&from_roots(&roots), &R::Disk),
        RootCounts {
            inside: 1,
            boundary: 1,
            outside: 1
        }
    );
    let roots = [
        G::new(Q::zero(), epsilon.clone()),
        G::new(Q::zero(), -epsilon),
    ];
    assert_eq!(
        counts(&from_roots(&roots), &R::Upper),
        RootCounts {
            inside: 1,
            boundary: 0,
            outside: 1
        }
    );
    let p = from_roots(&vec![g(0, 1); 8]);
    assert_eq!(
        counts(&p, &R::Upper),
        RootCounts {
            inside: 8,
            boundary: 0,
            outside: 0
        }
    );
    assert_eq!(
        counts(&p, &R::Disk),
        RootCounts {
            inside: 0,
            boundary: 8,
            outside: 0
        }
    );
    // The finite boundary point M(∞) of an arbitrary map is also tracked.
    let map = MobiusTransform::new(g(2, 3), g(1, -1), g(1, 0), g(2, 1)).unwrap();
    let p = from_roots(&[g(2, 3), g(2, 3)]);
    assert_eq!(map.transform_polynomial(&p).boundary_at_infinity, 2);
    assert_eq!(
        counts(&p, &R::from(map)),
        RootCounts {
            inside: 0,
            boundary: 2,
            outside: 0
        }
    );
}

#[test]
fn constants_linear_zero_and_errors() {
    for region in [R::Upper, R::Lower, R::Left, R::Right, R::Disk, R::Exterior] {
        assert_eq!(
            counts(&[g(2, -3), G::zero()], &region),
            RootCounts::default()
        );
        assert!(check(&[g(2, -3)], &region).unwrap().verdict);
        let zero = check(&[G::zero(), G::zero()], &region).unwrap();
        assert_eq!(zero.degree, None);
        assert_eq!(zero.counts, None);
        assert_eq!(zero.verdict, matches!(region, R::Upper | R::Lower));
        assert!(zero.verify(&[]));
    }
    assert!(MobiusTransform::new(g(1, 0), g(2, 0), g(2, 0), g(4, 0)).is_err());
    for radius in [Q::zero(), q(-1, 2)] {
        assert_eq!(
            check(
                &[],
                &R::CenteredDisk {
                    center: G::zero(),
                    radius
                }
            ),
            Err(StabilityError::NonpositiveRadius)
        );
    }
    assert_eq!(
        counts(&from_roots(&[g(0, 1)]), &R::Upper),
        RootCounts {
            inside: 1,
            boundary: 0,
            outside: 0
        }
    );
    assert_eq!(
        counts(&from_roots(&[g(0, -1)]), &R::Upper),
        RootCounts {
            inside: 0,
            boundary: 0,
            outside: 1
        }
    );
    // (z²-1)+iz has both roots below the real axis, fixing the Bézoutian sign.
    assert_eq!(
        counts(&[g(-1, 0), g(0, 1), g(1, 0)], &R::Upper),
        RootCounts {
            inside: 0,
            boundary: 0,
            outside: 2
        }
    );
}

#[test]
fn known_roots_in_every_region_with_multiplicities() {
    let roots = [
        g(0, 0),
        g(0, 0),
        g(1, 0),
        g(1, 0),
        g(-1, 0),
        g(0, 1),
        g(0, 1),
        g(0, -1),
        g(0, -1),
        G::new(q(1, 3), q(1, 4)),
        G::new(q(-2, 3), q(5, 2)),
        g(-2, 0),
        g(2, 0),
        g(1, 2),
        g(1, 2),
    ];
    let p = from_roots(&roots);
    // Expected counts are direct rational coordinate/norm comparisons, not
    // the inverse Möbius membership implementation.
    for (region, expected) in [
        (R::Upper, (6, 7, 2)),
        (R::Lower, (2, 7, 6)),
        (R::Left, (3, 6, 6)),
        (R::Right, (6, 6, 3)),
        (R::Disk, (3, 7, 5)),
        (R::Exterior, (5, 7, 3)),
    ] {
        let report = check(&p, &region).unwrap();
        assert_eq!(
            report.counts,
            Some(RootCounts {
                inside: expected.0,
                boundary: expected.1,
                outside: expected.2
            }),
            "{region:?}"
        );
        assert!(report.verify(&p));
        let certificate = report.certificate.as_ref().unwrap();
        let matrix: Vec<Vec<Q>> = certificate
            .bezoutian
            .iter()
            .map(|row| row.iter().cloned().map(Q::from_integer).collect())
            .collect();
        assert!(certificate.decomposition.verify(&matrix));
        let mut altered = report.clone();
        altered.counts.as_mut().unwrap().boundary += 1;
        assert!(!altered.verify(&p));
        let mut altered = report.clone();
        altered.certificate.as_mut().unwrap().common_factor[0] += 1;
        assert!(!altered.verify(&p));
    }
}

#[test]
fn degree_drop_cyclotomic_and_documented_examples() {
    let p = from_roots(&[g(1, 0), g(1, 0), g(1, 0)]);
    let report = check(&p, &R::Disk).unwrap();
    assert_eq!(
        report.certificate.unwrap().transformed.boundary_at_infinity,
        3
    );
    assert_eq!(
        counts(&p, &R::Disk),
        RootCounts {
            inside: 0,
            boundary: 3,
            outside: 0
        }
    );
    for p in [
        vec![1, 1],
        vec![1, 0, 1],
        vec![1, 1, 1],
        vec![1, 1, 1, 1, 1],
        vec![1, 0, 0, 0, 1],
    ] {
        let p: Vec<G> = p.into_iter().map(G::from).collect();
        assert_eq!(
            counts(&p, &R::Disk),
            RootCounts {
                inside: 0,
                boundary: p.len() - 1,
                outside: 0
            }
        );
    }
    let p: Vec<G> = [1, 3, 3, 1].into_iter().map(G::from).collect();
    assert!(check(&p, &R::Left).unwrap().verdict);
    assert!(check(&p, &R::Upper).unwrap().verdict);
    assert!(!check(&p, &R::Disk).unwrap().verdict);
    let p: Vec<G> = [2, 1, 1, 1].into_iter().map(G::from).collect();
    assert!(!check(&p, &R::Left).unwrap().verdict);
    assert_eq!(
        counts(&[g(1, 0), G::zero(), g(1, 0)], &R::Left),
        RootCounts {
            inside: 0,
            boundary: 2,
            outside: 0
        }
    );
    assert_eq!(
        check(&[g(-1, 0), G::zero(), g(1, 0)], &R::Upper)
            .unwrap()
            .method,
        StabilityMethod::RealSturm
    );
}

#[test]
fn centered_disk_and_explicit_conventions() {
    let region = R::CenteredDisk {
        center: g(2, -1),
        radius: q(3, 2),
    };
    let roots = [g(2, -1), G::new(q(7, 2), q(-1, 1)), g(4, -1)];
    let p = from_roots(&roots);
    assert_eq!(
        counts(&p, &region),
        RootCounts {
            inside: 1,
            boundary: 1,
            outside: 1
        }
    );
    for (root, sign) in roots.iter().zip([1, 0, -1]) {
        assert_eq!(region.membership_sign(root), Ok(sign));
    }
    let p = from_roots(&[g(0, -1)]);
    assert!(check(&p, &R::Upper).unwrap().verdict);
    assert!(
        !check_univariate_stability_with_convention(&p, &R::Upper, StabilityConvention::Inside)
            .unwrap()
            .verdict
    );
}

#[test]
fn mobius_round_trips_and_scaling() {
    let general = MobiusTransform::new(g(1, 2), g(-2, 1), g(2, -1), g(3, 1)).unwrap();
    for region in [
        R::Upper,
        R::Lower,
        R::Left,
        R::Right,
        R::Disk,
        R::Exterior,
        R::CenteredDisk {
            center: G::new(q(2, 3), q(-5, 7)),
            radius: q(4, 3),
        },
        R::from(general),
    ] {
        let map = region.mobius_transform().unwrap();
        for (w, sign) in [(g(1, 2), 1), (g(-2, -3), -1), (g(2, 0), 0)] {
            let Some(z) = map.evaluate(&w) else { continue };
            assert_eq!(map.inverse().evaluate(&z), Some(w.clone()));
            assert_eq!(map.compose(&map.inverse()).evaluate(&z), Some(z.clone()));
            assert_eq!(region.membership_sign(&z).unwrap(), sign);
            let p = from_roots(&[z]);
            let transformed = map.transform_polynomial(&p);
            let upper = check(&transformed.coefficients, &R::Upper).unwrap();
            assert_eq!(counts(&p, &region), upper.counts.unwrap());
        }
        let roots: Vec<_> = [g(1, 2), g(-2, -3), g(2, 0)]
            .iter()
            .filter_map(|w| map.evaluate(w))
            .collect();
        let p = from_roots(&roots);
        let transformed = map.transform_polynomial(&p);
        assert_eq!(transformed.boundary_at_infinity, 0);
        let round_trip = map
            .inverse()
            .transform_polynomial(&transformed.coefficients);
        let factor = round_trip
            .coefficients
            .last()
            .unwrap()
            .checked_divide(p.last().unwrap())
            .unwrap();
        assert_eq!(
            round_trip.coefficients,
            p.iter().map(|c| c.multiply(&factor)).collect::<Vec<_>>()
        );
        let scaled: Vec<_> = p
            .iter()
            .map(|c| c.multiply(&G::new(q(-3, 7), q(2, 5))))
            .collect();
        assert_eq!(counts(&p, &region), counts(&scaled, &region));
        let [a, b, c, d] = map.coefficients();
        let scale = g(2, -3);
        let scaled_map = MobiusTransform::new(
            a.multiply(&scale),
            b.multiply(&scale),
            c.multiply(&scale),
            d.multiply(&scale),
        )
        .unwrap();
        assert_eq!(counts(&p, &region), counts(&p, &R::from(scaled_map)));
    }
}

// Independent Schur--Cohn form, used only in tests. Its Hermitian inertia
// counts disk roots when nonsingular; reciprocal common factors make the
// form singular even without unit-circle roots and are checked separately.
fn schur_cohn(p: &[G]) -> Option<(usize, usize)> {
    let n = p.len() - 1;
    let mut h = vec![vec![G::zero(); n]; n];
    for (i, row) in h.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            for k in 0..=i.min(j) {
                *entry = entry
                    .add(&p[n - (i - k)].multiply(&p[n - (j - k)].conjugate()))
                    .subtract(&p[i - k].conjugate().multiply(&p[j - k]));
            }
        }
    }
    // Realification [[Re H, -Im H], [Im H, Re H]] repeats each eigenvalue.
    let mut real = vec![vec![Q::zero(); 2 * n]; 2 * n];
    for i in 0..n {
        for j in 0..n {
            real[i][j] = h[i][j].real.clone();
            real[n + i][n + j] = h[i][j].real.clone();
            real[i][n + j] = -&h[i][j].imaginary;
            real[n + i][j] = h[i][j].imaginary.clone();
        }
    }
    let inertia = symmetric_inertia(&real).unwrap();
    if inertia.zero != 0 {
        None
    } else {
        Some((inertia.positive / 2, inertia.negative / 2))
    }
}

#[test]
fn gaussian_random_polynomials_agree_with_schur_cohn() {
    let mut seed = 0x165a_3287_cde5u64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed % 9) as i64 - 4
    };
    let mut checked = 0;
    for n in 1..=12 {
        for _ in 0..12 {
            let mut p: Vec<_> = (0..=n)
                .map(|_| G::new(q(next(), 3), q(next(), 5)))
                .collect();
            p[n] = G::one();
            let Some((inside, outside)) = schur_cohn(&p) else {
                continue;
            };
            assert_eq!(
                counts(&p, &R::Disk),
                RootCounts {
                    inside,
                    boundary: 0,
                    outside
                },
                "{p:?}"
            );
            checked += 1;
        }
    }
    assert!(checked >= 120);
    // Reciprocal pairs: Schur--Cohn is singular, but the engine is complete.
    let reciprocal = from_roots(&[g(2, 0), G::from(q(1, 2))]);
    assert_eq!(schur_cohn(&reciprocal), None);
    assert_eq!(
        counts(&reciprocal, &R::Disk),
        RootCounts {
            inside: 1,
            boundary: 0,
            outside: 1
        }
    );
}
