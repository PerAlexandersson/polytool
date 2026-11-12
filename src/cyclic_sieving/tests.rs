//! Tests for the enclosing module.

use super::*;

fn b(values: &[i64]) -> Vec<BigInt> {
    values.iter().map(|&v| BigInt::from(v)).collect()
}

#[test]
fn cyclotomic_polynomials_are_exact() {
    assert_eq!(cyclotomic_polynomial_bigint(1), b(&[-1, 1]));
    assert_eq!(cyclotomic_polynomial_bigint(2), b(&[1, 1]));
    assert_eq!(cyclotomic_polynomial_bigint(3), b(&[1, 1, 1]));
    assert_eq!(cyclotomic_polynomial_bigint(4), b(&[1, 0, 1]));
    assert_eq!(
        crate::vec_poly::mul(&cyclotomic_polynomial_bigint(3), &b(&[-1, 1])),
        b(&[-1, 0, 0, 1])
    );
}

#[test]
fn root_of_unity_integer_profiles() {
    let p = b(&[1, 1, 1, 1]);
    let at_one = root_of_unity_evaluation_bigint(&p, 4, 0);
    assert_eq!(at_one.integer_value, Some(BigInt::from(4)));

    let at_minus_one = root_of_unity_evaluation_bigint(&p, 4, 2);
    assert_eq!(at_minus_one.integer_value, Some(BigInt::from(0)));

    let at_i = root_of_unity_evaluation_bigint(&p, 4, 1);
    assert_eq!(at_i.integer_value, Some(BigInt::from(0)));
}

#[test]
fn checks_basic_csp() {
    let p = b(&[1, 1]);
    let counts = b(&[2, 0]);
    let report = cyclic_sieving_report_bigint(&p, 2, Some(&counts)).unwrap();
    assert_eq!(report.holds, Some(true));
}

// ----- orbit decomposition and bicyclic sieving -----

fn big(v: i64) -> BigInt {
    BigInt::from(v)
}

/// Gaussian binomial `[n choose k]_q`, ascending coefficients.
fn gaussian_binomial(n: usize, k: usize) -> Vec<BigInt> {
    // [n,k] = [n-1,k-1] + q^k [n-1,k]
    let mut table = vec![vec![vec![BigInt::zero()]; k + 1]; n + 1];
    for row in table.iter_mut() {
        row[0] = vec![BigInt::one()];
    }
    for i in 1..=n {
        for j in 1..=k.min(i) {
            let mut sum = table[i - 1][j - 1].clone();
            let shifted = &table[i - 1][j];
            if !shifted.iter().all(Zero::is_zero) {
                let mut padded = vec![BigInt::zero(); j];
                padded.extend(shifted.iter().cloned());
                if padded.len() > sum.len() {
                    sum.resize(padded.len(), BigInt::zero());
                }
                for (x, y) in padded.iter().enumerate() {
                    sum[x] += y;
                }
            }
            strip_trailing_zeros(&mut sum);
            table[i][j] = sum;
        }
    }
    table[n][k].clone()
}

/// Orbit sizes of the rotation action of `C_n` on `k`-subsets of `[n]`.
fn rotation_orbit_sizes(n: usize, k: usize) -> Vec<usize> {
    let mask = (1u32 << n) - 1;
    let rotate = |x: u32| ((x << 1) | (x >> (n - 1))) & mask;
    let mut seen = std::collections::BTreeSet::new();
    let mut sizes = Vec::new();
    for x in 0..(1u32 << n) {
        if x.count_ones() as usize != k || seen.contains(&x) {
            continue;
        }
        let mut size = 0;
        let mut y = x;
        loop {
            seen.insert(y);
            size += 1;
            y = rotate(y);
            if y == x {
                break;
            }
        }
        sizes.push(size);
    }
    sizes
}

/// Fixed points `fix(g^d)` by direct enumeration.
fn rotation_fixed_points(n: usize, k: usize) -> Vec<BigInt> {
    let mask = (1u32 << n) - 1;
    (0..n)
        .map(|d| {
            let rot = |x: u32| ((x << d) | (x >> (n - d))) & mask;
            let count = (0..(1u32 << n))
                .filter(|&x| x.count_ones() as usize == k)
                .filter(|&x| d == 0 || rot(x) == x)
                .count();
            big(count as i64)
        })
        .collect()
}

#[test]
fn gaussian_binomial_helper_is_correct() {
    assert_eq!(gaussian_binomial(4, 2), b(&[1, 1, 2, 1, 1]));
    assert_eq!(gaussian_binomial(5, 1), b(&[1, 1, 1, 1, 1]));
}

#[test]
fn mobius_values() {
    let expected = [1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0];
    for (i, &e) in expected.iter().enumerate() {
        assert_eq!(mobius(i + 1), e, "mu({})", i + 1);
    }
}

#[test]
fn q_binomial_rotation_orbit_decomposition() {
    for (n, k) in [
        (4, 2),
        (5, 2),
        (6, 2),
        (6, 3),
        (7, 3),
        (8, 4),
        (9, 3),
        (10, 5),
        (10, 4),
        (12, 6),
    ] {
        let f = gaussian_binomial(n, k);
        let sizes = rotation_orbit_sizes(n, k);
        let fixed = rotation_fixed_points(n, k);
        // orbit sizes -> fixed points agrees with enumeration
        assert_eq!(fixed_points_from_orbit_sizes(&sizes, n).unwrap(), fixed);

        let dec = cyclic_sieving_orbit_decomposition_bigint(&f, n).unwrap();
        assert!(dec.integral, "[{n},{k}] integral");
        assert!(dec.natural_coefficients);
        assert_eq!(dec.cone_condition, Some(true));
        assert!(dec.realizable);
        assert_eq!(dec.fixed_points.as_ref().unwrap(), &fixed);
        // predicted orbit counts equal the actual ones
        let actual = orbit_counts_from_orbit_sizes(&sizes, n).unwrap();
        assert_eq!(dec.orbits.as_ref().unwrap(), &actual);
        // fixed points -> orbit counts round trip
        assert_eq!(orbit_counts_from_fixed_points(&fixed, n).unwrap(), actual);

        let check = check_cyclic_sieving_orbit_sizes(&f, n, &sizes).unwrap();
        assert!(check.holds);
        assert_eq!(check.sieving_report.holds, Some(true));
        assert!(check.mismatched_orbit_sizes.is_empty());
        assert_eq!(check.num_elements, f.iter().sum::<BigInt>());
    }
}

#[test]
fn q_binomial_wrong_action_is_rejected() {
    // [4 choose 2]_q = 1 + q + 2q^2 + q^3 + q^4, n = 4: orbits 4,2 (6 elements).
    let f = gaussian_binomial(4, 2);
    let good = check_cyclic_sieving_orbit_sizes(&f, 4, &[4, 2]).unwrap();
    assert!(good.holds);
    let bad = check_cyclic_sieving_orbit_sizes(&f, 4, &[4, 1, 1]).unwrap();
    assert!(!bad.holds);
    assert_eq!(bad.sieving_report.holds, Some(false));
    assert_eq!(bad.mismatched_orbit_sizes, vec![1, 2]);
    assert!(check_cyclic_sieving_orbit_sizes(&f, 4, &[3]).is_err());
}

#[test]
fn one_plus_q_order_two() {
    let dec = cyclic_sieving_orbit_decomposition_bigint(&b(&[1, 1]), 2).unwrap();
    assert!(dec.integral);
    let orbits = dec.orbits.unwrap();
    // one orbit of size 2, none of size 1
    assert_eq!(orbits[0].orbit_size, 1);
    assert_eq!(orbits[0].orbit_count, big(0));
    assert_eq!(orbits[1].orbit_size, 2);
    assert_eq!(orbits[1].orbit_count, big(1));
    assert_eq!(orbits[1].element_count, big(2));
    assert_eq!(dec.fixed_points.unwrap(), b(&[2, 0]));
    assert!(dec.realizable);
}

#[test]
fn non_integral_example() {
    // 1 + q mod q^3 - 1: c_1 = 1 but c_2 = 0 (both in class gcd = 1).
    let dec = cyclic_sieving_orbit_decomposition_bigint(&b(&[1, 1]), 3).unwrap();
    assert!(!dec.integral);
    assert_eq!(dec.violating_classes, vec![1]);
    assert!(dec.orbits.is_none() && dec.fixed_points.is_none());
    assert_eq!(dec.cone_condition, None);
    assert!(!dec.realizable);
    // agrees with exact evaluation: 1 + zeta_3 is not an integer
    assert!(root_of_unity_evaluation_bigint(&b(&[1, 1]), 3, 1)
        .integer_value
        .is_none());
    let check = check_cyclic_sieving_orbit_sizes(&b(&[1, 1]), 3, &[1, 1]).unwrap();
    assert!(!check.holds);
    assert_eq!(check.sieving_report.holds, Some(false));
}

#[test]
fn paper_example_is_integral_natural_but_not_realizable() {
    // Alexandersson-Amini, ex:counterex.
    let f = b(&[10, 1, 0, 3, 0, 1]);
    let dec = cyclic_sieving_orbit_decomposition_bigint(&f, 6).unwrap();
    assert!(dec.integral);
    assert!(dec.natural_coefficients);
    // f(zeta_6^j), j = 1..6 is 8,12,5,12,8,15 (j = 6 is d = 0)
    let fp = dec.fixed_points.clone().unwrap();
    let by_j: Vec<BigInt> = (1..=6).map(|j| fp[j % 6].clone()).collect();
    assert_eq!(by_j, b(&[8, 12, 5, 12, 8, 15]));
    // S_k, k = 1..6
    let orbits = dec.orbits.unwrap();
    let s: Vec<BigInt> = (1..=6)
        .map(|k| {
            orbits
                .iter()
                .find(|e| e.orbit_size == k)
                .map(|e| e.element_count.clone())
                .unwrap_or_else(BigInt::zero)
        })
        .collect();
    // S_4 and S_5 are not defined (k does not divide 6); the paper's formula gives 0.
    assert_eq!(s, b(&[8, 4, -3, 0, 0, 6]));
    assert_eq!(dec.cone_condition, Some(false));
    assert!(!dec.realizable);
    let neg = orbits.iter().find(|e| e.orbit_size == 3).unwrap();
    assert_eq!(neg.orbit_count, big(-1));
    // the fixed-point table is not that of an action
    assert!(orbit_counts_from_fixed_points(&fp, 6).is_err());
}

#[test]
fn orbit_counts_from_fixed_points_errors() {
    // wrong length
    assert_eq!(
        orbit_counts_from_fixed_points(&b(&[2, 0, 0]), 2).unwrap_err(),
        CyclicSievingError::FixedPointCountLength { order: 2, found: 3 }
    );
    // fix(g^1) depends only on gcd(1,4) = 1, so fix(g) = fix(g^3) must agree
    assert!(matches!(
        orbit_counts_from_fixed_points(&b(&[4, 1, 0, 0]), 4).unwrap_err(),
        CyclicSievingError::FixedPointsNotConstantOnGcdClasses { .. }
    ));
    // S_2 = fix(g^2) - fix(g) = 3 is odd
    assert!(matches!(
        orbit_counts_from_fixed_points(&b(&[5, 2, 5, 2]), 4).unwrap_err(),
        CyclicSievingError::OrbitElementCountNotDivisible { orbit_size: 2, .. }
    ));
    // negative S_2
    assert!(matches!(
        orbit_counts_from_fixed_points(&b(&[3, 2, 1, 2]), 4).unwrap_err(),
        CyclicSievingError::NegativeOrbitElementCount { orbit_size: 2, .. }
    ));
}

#[test]
fn error_messages_are_stable() {
    let message = |error: CyclicSievingError| error.to_string();
    assert_eq!(
        message(cyclic_sieving_report_bigint(&b(&[1]), 0, None).unwrap_err()),
        "cyclic group order must be positive"
    );
    assert_eq!(
        message(cyclic_sieving_report_bigint(&b(&[1]), 2, Some(&b(&[1]))).unwrap_err()),
        "expected 2 fixed-point counts for order 2, got 1"
    );
    assert_eq!(
        message(fixed_points_from_orbit_sizes(&[4], 6).unwrap_err()),
        "orbit size 4 does not divide the order 6"
    );
    assert_eq!(
        message(fixed_points_from_orbit_counts(&[(2, big(-1))], 4).unwrap_err()),
        "negative number of orbits of size 2"
    );
    assert_eq!(
        message(orbit_counts_from_fixed_points(&b(&[4, 1, 0, 0]), 4).unwrap_err()),
        "fixed-point counts are not those of a C_4-action: fix(g^3) = 0 but fix(g^1) = 1 \
         although gcd(3, 4) = 1"
    );
    assert_eq!(
        message(orbit_counts_from_fixed_points(&b(&[5, 2, 5, 2]), 4).unwrap_err()),
        "fixed-point counts are not those of a C_4-action: S_2 = 3 is not divisible by 2"
    );
    assert_eq!(
        message(orbit_counts_from_fixed_points(&b(&[3, 2, 1, 2]), 4).unwrap_err()),
        "fixed-point counts are not those of a C_4-action: S_2 = -1 is negative"
    );
    assert_eq!(
        message(bicyclic_sieving_report_bigint(&[vec![big(1)]], 0, 2, None).unwrap_err()),
        "cyclic group orders must be positive"
    );
    assert_eq!(
        message(bicyclic_sieving_report_bigint(&[vec![big(1)]], 51, 50, None).unwrap_err()),
        "group order 51 x 50 exceeds the limit 2500"
    );
    assert_eq!(
        message(
            bicyclic_sieving_report_bigint(&[vec![big(1)]], 2, 2, Some(&[vec![big(1)]]))
                .unwrap_err()
        ),
        "fixed-point table must be a 2 x 2 matrix"
    );
}

#[test]
fn orbit_size_conversion_roundtrip() {
    let sizes = [1, 1, 2, 3, 3, 6, 6, 6];
    let fixed = fixed_points_from_orbit_sizes(&sizes, 6).unwrap();
    let counts = orbit_counts_from_fixed_points(&fixed, 6).unwrap();
    let expected = orbit_counts_from_orbit_sizes(&sizes, 6).unwrap();
    assert_eq!(counts, expected);
    let pairs: Vec<(usize, BigInt)> = expected
        .iter()
        .map(|e| (e.orbit_size, e.orbit_count.clone()))
        .collect();
    assert_eq!(fixed_points_from_orbit_counts(&pairs, 6).unwrap(), fixed);
    assert!(fixed_points_from_orbit_sizes(&[4], 6).is_err());
    assert!(fixed_points_from_orbit_sizes(&[0], 6).is_err());
    assert!(
        fixed_points_from_orbit_counts(&[(2, big(-1))], 6).is_err(),
        "negative counts rejected"
    );
}

struct Lcg(u64);
impl Lcg {
    fn next(&mut self, bound: u64) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) % bound
    }
}

#[test]
fn orbit_formula_matches_cyclotomic_evaluation_on_random_polynomials() {
    let mut rng = Lcg(12345);
    for trial in 0..200 {
        let n = 1 + rng.next(12) as usize;
        let len = 1 + rng.next(20) as usize;
        // Half of the trials are forced integral by building f from orbit polynomials.
        let coeffs: Vec<BigInt> = if trial % 2 == 0 {
            (0..len).map(|_| big(rng.next(7) as i64 - 2)).collect()
        } else {
            let mut f = vec![BigInt::zero(); n];
            for k in divisors(n) {
                let m = rng.next(4) as i64 - 1;
                // orbit of size k has polynomial h_{n/k}
                let r = n / k;
                for i in 0..k {
                    f[(r * i) % n] += m;
                }
            }
            f
        };
        let dec = cyclic_sieving_orbit_decomposition_bigint(&coeffs, n).unwrap();
        let all_integral = (0..n).all(|d| {
            root_of_unity_evaluation_bigint(&coeffs, n, d)
                .integer_value
                .is_some()
        });
        assert_eq!(dec.integral, all_integral, "n={n} f={coeffs:?}");
        if trial % 2 == 1 {
            assert!(dec.integral);
        }
        if let Some(fixed) = &dec.fixed_points {
            for (d, value) in fixed.iter().enumerate() {
                assert_eq!(
                    Some(value),
                    root_of_unity_evaluation_bigint(&coeffs, n, d)
                        .integer_value
                        .as_ref(),
                    "n={n} d={d} f={coeffs:?}"
                );
            }
            // S_k = sum_{j | k} mu(k/j) f(zeta^j) and a_r reconstruct f mod q^n - 1
            let orbits = dec.orbits.as_ref().unwrap();
            let mut rebuilt = vec![BigInt::zero(); n];
            for entry in orbits {
                for i in 0..entry.orbit_size {
                    rebuilt[(entry.basis_index * i) % n] += &entry.orbit_count;
                }
            }
            assert_eq!(rebuilt, dec.reduced_residues);
        }
    }
}

#[test]
fn integral_natural_polynomial_with_negative_orbit_count_via_remainders() {
    // Same example through the general API at order 6 with fixed counts of an
    // existing action never matching: CSP fails for every candidate action.
    let f = b(&[10, 1, 0, 3, 0, 1]);
    for sizes in [vec![1; 15], vec![6, 6, 3], vec![6, 6, 1, 1, 1]] {
        let check = check_cyclic_sieving_orbit_sizes(&f, 6, &sizes).unwrap();
        assert!(!check.holds, "{sizes:?}");
    }
}

#[test]
fn bicyclic_one_plus_qt() {
    let f = vec![b(&[1, 0]), b(&[0, 1])];
    // diagonal action of C_2 x C_2 on two points, swapped by both generators
    let table = vec![b(&[2, 0]), b(&[0, 2])];
    let report = bicyclic_sieving_report_bigint(&f, 2, 2, Some(&table)).unwrap();
    assert!(report.integral);
    assert!(report.product_basis_form);
    assert_eq!(report.holds, Some(true));
    assert_eq!(report.divisors1, vec![1, 2]);
    assert_eq!(
        report.product_coordinates.as_ref().unwrap(),
        &vec![b(&[1, -1]), b(&[-1, 2])]
    );
    assert_eq!(report.product_formula_evaluations.as_ref().unwrap(), &table);
    let evaluations: Vec<Vec<Option<BigInt>>> = table
        .iter()
        .map(|row| row.iter().cloned().map(Some).collect())
        .collect();
    assert_eq!(report.evaluations, evaluations);

    let wrong = vec![b(&[2, 0]), b(&[0, 0])];
    let report = bicyclic_sieving_report_bigint(&f, 2, 2, Some(&wrong)).unwrap();
    assert_eq!(report.holds, Some(false));
    assert_eq!(report.mismatches, vec![(1, 1)]);
}

fn matrix_product(p: &[BigInt], q: &[BigInt]) -> Vec<Vec<BigInt>> {
    p.iter()
        .map(|a| q.iter().map(|c| a * c).collect())
        .collect()
}

#[test]
fn genuine_bicyclic_csp_from_two_q_binomials() {
    for ((n1, k1), (n2, k2)) in [((4, 2), (6, 3)), ((5, 2), (4, 2)), ((6, 2), (6, 3))] {
        let f = matrix_product(&gaussian_binomial(n1, k1), &gaussian_binomial(n2, k2));
        let fix1 = rotation_fixed_points(n1, k1);
        let fix2 = rotation_fixed_points(n2, k2);
        let table: Vec<Vec<BigInt>> = fix1
            .iter()
            .map(|a| fix2.iter().map(|c| a * c).collect())
            .collect();
        let report = bicyclic_sieving_report_bigint(&f, n1, n2, Some(&table)).unwrap();
        assert!(report.integral);
        assert_eq!(report.holds, Some(true), "({n1},{k1}) x ({n2},{k2})");
        // evaluations factor
        let u1 = cyclic_sieving_orbit_decomposition_bigint(&gaussian_binomial(n1, k1), n1).unwrap();
        let u2 = cyclic_sieving_orbit_decomposition_bigint(&gaussian_binomial(n2, k2), n2).unwrap();
        assert_eq!(u1.fixed_points.unwrap(), fix1);
        assert_eq!(u2.fixed_points.unwrap(), fix2);
        // a product of univariate polynomials in product form has product coordinates
        if report.product_basis_form {
            assert_eq!(report.product_formula_evaluations.as_ref().unwrap(), &table);
        }
    }
}

#[test]
fn bicyclic_product_coordinates_are_products_of_univariate_ones() {
    let f1 = gaussian_binomial(6, 3);
    let f2 = gaussian_binomial(4, 2);
    let f = matrix_product(&f1, &f2);
    let report = bicyclic_sieving_report_bigint(&f, 6, 4, None).unwrap();
    assert!(report.product_basis_form);
    let a1 = cyclic_sieving_orbit_decomposition_bigint(&f1, 6).unwrap();
    let a2 = cyclic_sieving_orbit_decomposition_bigint(&f2, 4).unwrap();
    let coords = report.product_coordinates.unwrap();
    for (x1, &r1) in report.divisors1.iter().enumerate() {
        for (x2, &r2) in report.divisors2.iter().enumerate() {
            let e1 = a1
                .orbits
                .as_ref()
                .unwrap()
                .iter()
                .find(|e| e.basis_index == r1);
            let e2 = a2
                .orbits
                .as_ref()
                .unwrap()
                .iter()
                .find(|e| e.basis_index == r2);
            let (e1, e2) = (e1.unwrap(), e2.unwrap());
            assert_eq!(coords[x1][x2], &e1.orbit_count * &e2.orbit_count);
        }
    }
}

#[test]
fn bicyclic_integral_but_not_product_form_when_orders_share_a_factor() {
    // f = sum_{i,j} (zeta_3^{i+j} + zeta_3^{-(i+j)}) q^i t^j: 2 if i + j = 0 mod 3, else -1.
    let f: Vec<Vec<BigInt>> = (0..3)
        .map(|i| {
            (0..3)
                .map(|j| if (i + j) % 3 == 0 { big(2) } else { big(-1) })
                .collect()
        })
        .collect();
    let report = bicyclic_sieving_report_bigint(&f, 3, 3, None).unwrap();
    assert!(report.integral);
    assert!(!report.product_basis_form);
    assert!(report.product_coordinates.is_none());
    let values: std::collections::BTreeSet<BigInt> = report
        .evaluations
        .iter()
        .flatten()
        .map(|v| v.clone().unwrap())
        .collect();
    assert_eq!(values, [big(0), big(9)].into_iter().collect());
}

#[test]
fn bicyclic_non_integral_and_bad_input() {
    let f = vec![b(&[1]), b(&[1])];
    let report = bicyclic_sieving_report_bigint(&f, 3, 1, None).unwrap();
    assert!(!report.integral);
    assert!(bicyclic_sieving_report_bigint(&f, 0, 1, None).is_err());
    assert!(bicyclic_sieving_report_bigint(&f, 2, 2, Some(&[b(&[1, 1])])).is_err());
    assert!(bicyclic_sieving_report_bigint(&f, 3, 1, Some(&[b(&[1])])).is_err());
    assert!(bicyclic_sieving_report_bigint(&f, 100, 100, None).is_err());
}

#[test]
fn bicyclic_univariate_case_agrees_with_cyclic() {
    let f = gaussian_binomial(6, 3);
    let matrix: Vec<Vec<BigInt>> = f.iter().map(|c| vec![c.clone()]).collect();
    let report = bicyclic_sieving_report_bigint(&matrix, 6, 1, None).unwrap();
    let dec = cyclic_sieving_orbit_decomposition_bigint(&f, 6).unwrap();
    let col: Vec<BigInt> = report
        .evaluations
        .iter()
        .map(|row| row[0].clone().unwrap())
        .collect();
    assert_eq!(col, dec.fixed_points.unwrap());
}

#[test]
fn sequence_reports_reject_overflowing_and_oversized_candidate_orders() {
    let polynomials = vec![vec![BigInt::from(1), BigInt::from(1)]];
    // A candidate order that does not fit in `isize` is an error, not a skip.
    let overflow =
        cyclic_sieving_sequence_reports_bigint(&polynomials, 1, &[isize::MAX, 1], &BTreeMap::new())
            .unwrap_err();
    assert_eq!(
        overflow,
        CyclicSievingError::OrderTooLarge {
            order: isize::MAX as u128 + 1,
            maximum: CYCLIC_SIEVING_MAX_ORDER
        }
    );
    // So is one that fits but is far beyond the supported maximum.
    let huge = (CYCLIC_SIEVING_MAX_ORDER as isize) + 1;
    let oversized =
        cyclic_sieving_sequence_reports_bigint(&polynomials, 0, &[huge], &BTreeMap::new())
            .unwrap_err();
    assert_eq!(
        oversized,
        CyclicSievingError::OrderTooLarge {
            order: huge as u128,
            maximum: CYCLIC_SIEVING_MAX_ORDER
        }
    );
    assert!(oversized
        .to_string()
        .contains("exceeds the supported maximum"));
    // Non-positive candidate orders are still skipped, and small ones work.
    let items = cyclic_sieving_sequence_reports_bigint(&polynomials, 1, &[-5, 1], &BTreeMap::new())
        .unwrap();
    assert_eq!(items.len(), 1);
    assert_eq!(items[0].index, 1);
    assert_eq!(items[0].candidate_orders.len(), 1);
    assert_eq!(items[0].candidate_orders[0].order, 2);
}

#[test]
fn sequence_reports_reject_an_index_that_overflows() {
    let polynomials = vec![vec![BigInt::from(1)]; 2];
    let error =
        cyclic_sieving_sequence_reports_bigint(&polynomials, isize::MAX, &[0], &BTreeMap::new())
            .unwrap_err();
    assert_eq!(
        error,
        CyclicSievingError::SequenceIndexOverflow {
            first_index: isize::MAX,
            row: 1
        }
    );
    assert!(error.to_string().contains("overflows"));
}

#[test]
fn report_rejects_an_order_beyond_the_maximum_without_allocating() {
    let p = vec![BigInt::from(1)];
    for order in [CYCLIC_SIEVING_MAX_ORDER + 1, usize::MAX / 2, usize::MAX] {
        assert_eq!(
            cyclic_sieving_report_bigint(&p, order, None).unwrap_err(),
            CyclicSievingError::OrderTooLarge {
                order: order as u128,
                maximum: CYCLIC_SIEVING_MAX_ORDER
            }
        );
    }
    assert!(cyclic_sieving_report_bigint(&p, CYCLIC_SIEVING_MAX_ORDER, Some(&[])).is_err());
}
