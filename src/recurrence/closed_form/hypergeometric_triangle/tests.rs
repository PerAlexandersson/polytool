use super::*;

fn q(v: i64) -> Q {
    Q::from_integer(BigInt::from(v))
}

fn binomial(n: i64, k: i64) -> Q {
    if k < 0 || k > n {
        return Q::zero();
    }
    let mut value = BigInt::one();
    for i in 0..k {
        value = value * BigInt::from(n - i) / BigInt::from(i + 1);
    }
    Q::from_integer(value)
}

fn factorial(n: i64) -> Q {
    (1..=n).fold(q(1), |acc, i| acc * q(i))
}

/// The derived record reproduces the closed form for `count` rows, the
/// certificate replays, and the operator is returned.
fn derive_and_check(
    triangle: &HypergeometricTriangle,
    closed_form: impl Fn(i64, i64) -> Q,
    count: usize,
) -> HypergeometricRowRecurrenceCertificate {
    let certificate = derive_hypergeometric_row_recurrence(
        triangle,
        &HypergeometricRowRecurrenceBudget::default(),
    )
    .expect("derived");
    certificate.verify().expect("replays");
    let n0 = triangle.first_index().unwrap();
    let direct = triangle.generate_rows(count).unwrap();
    let (rule, first, seeds) = certificate.record.to_recurrence_parts().unwrap();
    assert_eq!(first, n0);
    let generated = rule.generate_rows_rational(&seeds, first, count).unwrap();
    for (offset, (row, record_row)) in direct.iter().zip(&generated).enumerate() {
        let n = (n0 + offset) as i64;
        let expected: Vec<Q> = vec_poly::trim(
            &(0..=3 * n + 3)
                .map(|k| closed_form(n, k))
                .collect::<Vec<_>>(),
        );
        assert_eq!(row, &expected, "definition row {n}");
        assert_eq!(&vec_poly::trim(record_row), &expected, "record row {n}");
    }
    certificate
}

#[test]
fn binomial_rows_satisfy_the_first_order_rule() {
    let triangle = HypergeometricTriangle::single(HypergeometricTermDefinition::new(
        0,
        "(n-k)/(k+1)",
        "(n+1)/(n+1-k)",
        "0",
        "n",
        "1",
    ));
    let certificate = derive_and_check(&triangle, binomial, 12);
    // P_(n+1) = (1 + t) P_n.
    assert_eq!(certificate.operator.len(), 3);
    assert_eq!(
        certificate.operator_text,
        "(1) P(n+1) + (-1) t P(n) + (-1) P(n) = 0"
    );
}

#[test]
fn narayana_rows_have_a_certified_rule() {
    let triangle = HypergeometricTriangle::single(HypergeometricTermDefinition::new(
        1,
        "(n-k)*(n-k+1)/(k*(k+1))",
        "n*(n+1)/((n+1-k)*(n+2-k))",
        "1",
        "n",
        "1",
    ));
    let narayana = |n: i64, k: i64| binomial(n, k) * binomial(n, k - 1) / q(n);
    let certificate = derive_and_check(&triangle, narayana, 14);
    assert_eq!(certificate.record.first_index, 1);
    // The leader is a polynomial in n only, and the rule has order at most two.
    assert!(certificate.operator.iter().all(|term| term.shift <= 2));
}

#[test]
fn squared_binomial_rows_have_a_certified_rule() {
    let triangle = HypergeometricTriangle::single(HypergeometricTermDefinition::new(
        0,
        "(n-k)^2/(k+1)^2",
        "(n+1)^2/(n+1-k)^2",
        "0",
        "n",
        "1",
    ));
    derive_and_check(&triangle, |n, k| binomial(n, k) * binomial(n, k), 14);
}

#[test]
fn trinomial_factorial_ratio_rows_have_a_certified_rule() {
    // T(n, k) = (3n - 2k)! / ((n - k)!^3 k!).
    let triangle = HypergeometricTriangle::single(HypergeometricTermDefinition::new(
        0,
        "(n-k)^3/((3n-2k)*(3n-2k-1)*(k+1))",
        "(3n-2k+1)*(3n-2k+2)*(3n-2k+3)/(n+1-k)^3",
        "0",
        "n",
        "1",
    ));
    let closed = |n: i64, k: i64| {
        if k < 0 || k > n {
            return Q::zero();
        }
        factorial(3 * n - 2 * k) / (factorial(n - k).pow(3i32) * factorial(k))
    };
    derive_and_check(&triangle, closed, 12);
}

#[test]
fn squared_binomial_times_central_binomial_rows_have_a_certified_rule() {
    // T(n, k) = C(n, k)^2 C(2(n - k), n - k).
    let triangle = HypergeometricTriangle::single(HypergeometricTermDefinition::new(
        0,
        "(n-k)^3/(2*(k+1)^2*(2n-2k-1))",
        "2*(n+1)^2*(2n-2k+1)/(n+1-k)^3",
        "0",
        "n",
        "1",
    ));
    let closed = |n: i64, k: i64| binomial(n, k).pow(2i32) * binomial(2 * (n - k), n - k);
    derive_and_check(&triangle, closed, 12);
}

#[test]
fn parity_dependent_terms_and_sloped_supports_are_supported() {
    // T(2m, k) = (m! C(m, k))^2, T(2m + 1, k) = ((m + 1)! C(m, k))^2 / (k + 1),
    // 0 <= k <= floor(n / 2).
    let term = HypergeometricTermDefinition {
        first_index: 0,
        column_ratios: vec![
            "(n/2-k)^2/(k+1)^2".into(),
            "(n/2-1/2-k)^2/((k+1)*(k+2))".into(),
        ],
        index_ratios: vec![
            "(n+2)^2/(4*(k+1))".into(),
            "(k+1)*(n+1)^2/(n+1-2k)^2".into(),
        ],
        lower_column: "0".into(),
        upper_column: "n/2".into(),
        first_value: "1".into(),
        polynomial_factor: None,
    };
    let closed = |n: i64, k: i64| {
        let m = n / 2;
        if k < 0 || k > m {
            return Q::zero();
        }
        if n % 2 == 0 {
            (factorial(m) * binomial(m, k)).pow(2i32)
        } else {
            (factorial(m + 1) * binomial(m, k)).pow(2i32) / q(k + 1)
        }
    };
    derive_and_check(&HypergeometricTriangle::single(term), closed, 16);
}

#[test]
fn sums_of_terms_are_annihilated_jointly() {
    // T(n, k) = C(n, k) + 2^n C(n, k): rows (1 + 2^n) (1 + t)^n.
    let triangle = HypergeometricTriangle {
        terms: vec![
            HypergeometricTermDefinition::new(0, "(n-k)/(k+1)", "(n+1)/(n+1-k)", "0", "n", "1"),
            HypergeometricTermDefinition::new(0, "(n-k)/(k+1)", "2(n+1)/(n+1-k)", "0", "n", "1"),
        ],
    };
    let closed = |n: i64, k: i64| binomial(n, k) * (q(1) + q(2).pow(n as i32));
    let certificate = derive_and_check(&triangle, closed, 10);
    assert!(certificate.operator.iter().any(|term| term.shift == 2));
}

#[test]
fn boundary_contributions_are_not_assumed_to_vanish() {
    // T(n, k) = C(n, k) for 1 <= k <= n: rows (1 + t)^n - 1.  The interior
    // identity of P_(n+1) = (1 + t) P_n holds, the strip k = lower(n) fails.
    let triangle = HypergeometricTriangle::single(HypergeometricTermDefinition::new(
        1,
        "(n-k)/(k+1)",
        "(n+1)/(n+1-k)",
        "1",
        "n",
        "1",
    ));
    let operator = vec![
        EulerOperatorTerm {
            shift: 0,
            variable_power: 0,
            euler_power: 0,
            coefficient: vec!["-1".into()],
        },
        EulerOperatorTerm {
            shift: 0,
            variable_power: 1,
            euler_power: 0,
            coefficient: vec!["-1".into()],
        },
        EulerOperatorTerm {
            shift: 1,
            variable_power: 0,
            euler_power: 0,
            coefficient: vec!["1".into()],
        },
    ];
    let error = certify_hypergeometric_row_operator(&triangle, &operator).unwrap_err();
    assert!(error.to_string().contains("boundary strip"), "{error}");
    // The search imposes the strips and finds n P_(n+1) = n P_n + (1 + t) theta P_n.
    let closed = |n: i64, k: i64| if k == 0 { Q::zero() } else { binomial(n, k) };
    let certificate = derive_and_check(&triangle, closed, 10);
    assert_eq!(
        certificate.operator_text,
        "(n) P(n+1) + (-1) t theta P(n) + (-1) theta P(n) + (-n) P(n) = 0"
    );
}

#[test]
fn inconsistent_ratios_and_vanishing_factors_are_rejected() {
    let inconsistent = HypergeometricTriangle::single(HypergeometricTermDefinition::new(
        0,
        "(n-k)/(k+1)",
        "(n+1)/(n+2-k)",
        "0",
        "n",
        "1",
    ));
    let error = inconsistent.generate_rows(3).unwrap_err();
    assert!(error.to_string().contains("not consistent"), "{error}");
    // (k - 2) vanishes inside the support.
    let vanishing = HypergeometricTriangle::single(HypergeometricTermDefinition::new(
        0, "(k-2)", "1", "0", "n", "1",
    ));
    let error = vanishing.generate_rows(3).unwrap_err();
    assert!(error.to_string().contains("may vanish"), "{error}");
}

#[test]
fn certificates_survive_a_json_round_trip() {
    let triangle = HypergeometricTriangle::single(HypergeometricTermDefinition::new(
        0,
        "(n-k)/(k+1)",
        "(n+1)/(n+1-k)",
        "0",
        "n",
        "1",
    ));
    let certificate = derive_hypergeometric_row_recurrence(
        &triangle,
        &HypergeometricRowRecurrenceBudget::default(),
    )
    .unwrap();
    let json = serde_json::to_string(&certificate).unwrap();
    let back: HypergeometricRowRecurrenceCertificate = serde_json::from_str(&json).unwrap();
    back.verify().unwrap();
    let mut tampered = back.clone();
    tampered.operator[0].coefficient = vec!["-2".into()];
    assert!(tampered.verify().is_err());
}

#[test]
fn polynomial_factors_cover_sums_of_similar_terms() {
    // T(n, k) = C(n, k) (C(n-1, k-1) C(n+1, k+1) + C(n-1, k) C(n+1, k))
    //         = h(n, k) (2nk + n - 2k^2),  h = C(n,k) C(n+1,k+1) (n-1)! / (k! (n-k)! (n+1-k)).
    let triangle = HypergeometricTriangle::single(
        HypergeometricTermDefinition::new(
            1,
            "(n-k)^2*(n+1-k)/((k+1)^2*(k+2))",
            "n*(n+1)*(n+2)/((n+1-k)^2*(n+2-k))",
            "0",
            "n",
            "1",
        )
        .with_polynomial_factor("2n*k + n - 2k^2"),
    );
    let closed = |n: i64, k: i64| {
        binomial(n, k)
            * (binomial(n - 1, k - 1) * binomial(n + 1, k + 1)
                + binomial(n - 1, k) * binomial(n + 1, k))
    };
    derive_and_check(&triangle, closed, 12);
}

#[test]
fn a_three_term_record_is_certified_against_the_closed_form() {
    // (n + 1) P_(n+1) = (2n + 1)(1 + t) P_n - n (1 - t)^2 P_(n-1) for
    // P_n = sum_k C(n, k)^2 t^k, written at the target index m = n + 1.
    let triangle = HypergeometricTriangle::single(HypergeometricTermDefinition::new(
        0,
        "(n-k)^2/(k+1)^2",
        "(n+1)^2/(n+1-k)^2",
        "0",
        "n",
        "1",
    ));
    let record = |lag_two: [&str; 3]| -> RecurrenceJson {
        serde_json::from_value(serde_json::json!({
            "schema": "polytool.recurrence.v1",
            "first_index": 0,
            "initial_polynomials": [["1"], ["1", "1"]],
            "recurrence": {
                "terms": [
                    {"offset": 1, "deriv_order": 0, "sign": "none",
                     "coeff": {"coeffs": [["-1", "-1"], ["2", "2"]]}},
                    {"offset": 2, "deriv_order": 0, "sign": "none",
                     "coeff": {"coeffs": [[lag_two[0], lag_two[1], lag_two[2]],
                                          ["-1", "2", "-1"]]}}
                ],
                "denominator": {"coeffs": [["0"], ["1"]]}
            }
        }))
        .unwrap()
    };
    let certificate =
        certify_hypergeometric_row_record(&triangle, &record(["1", "-2", "1"])).unwrap();
    certificate.verify().unwrap();
    assert_eq!(certificate.variable_power_multiplier, 0);
    // A wrong lag-two coefficient fails exactly.
    let error =
        certify_hypergeometric_row_record(&triangle, &record(["1", "-2", "2"])).unwrap_err();
    assert!(
        matches!(error, RecurrenceError::VerificationFailed(_)),
        "{error}"
    );
}
