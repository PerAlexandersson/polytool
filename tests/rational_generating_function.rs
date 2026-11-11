//! Exact row generating functions, including seed corrections and cancellation.

use num_traits::One;
use ore_algebra::BivariatePolynomial as B;
use polytool::recurrence::{
    rational_generating_function, riordan_decomposition, BigRational as Q, BivarPoly,
    RationalGeneratingFunction, Recurrence, RecurrenceError, RecurrenceTerm,
};
use polytool::vec_poly;

fn q(n: i64) -> Q {
    Q::from_integer(n.into())
}
fn row(p: &[i64]) -> Vec<Q> {
    p.iter().map(|&c| q(c)).collect()
}
fn rule(terms: &[(usize, &[i64])]) -> Recurrence {
    Recurrence {
        terms: terms
            .iter()
            .map(|(lag, p)| {
                RecurrenceTerm::new(
                    *lag,
                    0,
                    BivarPoly::from_coefficients_normalized(vec![row(p)]),
                )
            })
            .collect(),
        denominator: None,
        inhomogeneous: None,
    }
}
fn polynomial(p: &[Vec<String>]) -> B {
    B::new(
        p.iter()
            .map(|r| r.iter().map(|s| s.parse().unwrap()).collect())
            .collect(),
    )
}

fn assert_expansion(
    rule: &Recurrence,
    first: usize,
    seeds: &[Vec<Q>],
    count: usize,
) -> RationalGeneratingFunction {
    let gf = rational_generating_function(rule, first, seeds).unwrap();
    let generated = rule.generate_rows_rational(seeds, first, count).unwrap();
    let n = polynomial(&gf.numerator);
    let d = polynomial(&gf.denominator);
    assert_eq!(d.rows()[0], [Q::one()]);
    let product = d.multiply(&B::new(generated));
    assert_eq!(
        B::new(product.rows().iter().take(count).cloned().collect()),
        n
    );
    assert_eq!(
        serde_json::from_str::<RationalGeneratingFunction>(&serde_json::to_string(&gf).unwrap())
            .unwrap(),
        gf
    );
    assert!(gf.to_string().contains("/"));
    gf
}

#[test]
fn fibonacci_and_two_catalog_families() {
    let r = rule(&[(1, &[1]), (2, &[0, 1])]);
    let gf = assert_expansion(&r, 7, &[row(&[1]), row(&[1])], 30);
    assert_eq!(gf.first_index, 7);
    assert_eq!(gf.numerator, vec![vec!["1"]]);
    assert_eq!(gf.denominator, vec![vec!["1"], vec!["-1"], vec!["0", "-1"]]);
    for id in ["A216182", "A114123"] {
        let entry = polytool::oeis::by_id(id).unwrap();
        let (r, seeds) = entry.recurrence_parts().unwrap();
        assert_expansion(&r, entry.recurrence_first_index, &seeds, 30);
    }
}

#[test]
fn all_prescribed_rows_including_prefix_and_zero_rows() {
    let r = rule(&[(1, &[1]), (3, &[0, 1])]);
    let seeds = vec![
        row(&[1]),
        row(&[1, 4]),
        row(&[1, 2, 1]),
        row(&[0]),
        row(&[9, 2]),
    ];
    let gf = assert_expansion(&r, 12, &seeds, 30);
    assert!(gf.numerator.len() > 3);
    assert!(riordan_decomposition(&gf).is_none());
    let ordinary =
        rational_generating_function(&r, 12, &[row(&[1]), row(&[1]), row(&[1])]).unwrap();
    assert_ne!(ordinary.numerator, gf.numerator);
    let altered = r.generate_rows_rational(&seeds[..3], 0, 5).unwrap();
    assert!(!polytool::is_real_rooted_bigint_coeffs(
        &altered[4]
            .iter()
            .map(|c| c.to_integer())
            .collect::<Vec<_>>()
    ));
}

#[test]
fn cancels_a_genuine_bivariate_factor_and_polynomial_content() {
    // (1-xz)(1-z) F = 1-xz, from constant rows, then cancel 1-xz.
    let r = rule(&[(1, &[1, 1]), (2, &[0, -1])]);
    let gf = assert_expansion(&r, 0, &[row(&[1]), row(&[1])], 30);
    assert_eq!(gf.numerator, vec![vec!["1"]]);
    assert_eq!(gf.denominator, vec![vec!["1"], vec!["-1"]]);
    // A nonconstant left factor may cancel before scalar normalization.
    let r = Recurrence {
        terms: vec![RecurrenceTerm::new(
            1,
            0,
            BivarPoly::from_coefficients_normalized(vec![row(&[2, 2])]),
        )],
        denominator: Some(BivarPoly::from_coefficients_normalized(vec![row(&[1, 1])])),
        inhomogeneous: None,
    };
    let gf = assert_expansion(&r, 0, &[row(&[1])], 30);
    assert_eq!(gf.denominator, vec![vec!["1"], vec!["-2"]]);
}

#[test]
fn riordan_series_are_reduced_and_may_be_noninvertible() {
    let r = rule(&[(1, &[1]), (3, &[0, 1])]);
    let gf = assert_expansion(&r, 0, &[row(&[1]), row(&[1]), row(&[1])], 30);
    let decomposition = riordan_decomposition(&gf).unwrap();
    assert_eq!(decomposition.k.numerator, vec!["0", "0", "1"]);
    assert_eq!(decomposition.k.denominator, vec!["1"]);
    let (dn, dd) = decomposition.d.coefficients().unwrap();
    let (hn, hd) = decomposition.h.coefficients().unwrap();
    let (kn, kd) = decomposition.k.coefficients().unwrap();
    assert_eq!(
        vec_poly::trim(&vec_poly::mul(&hn, &vec_poly::mul(&dd, &kd))),
        vec_poly::trim(&vec_poly::shift_by(
            &vec_poly::mul(&dn, &vec_poly::mul(&kn, &hd)),
            1
        ))
    );
}

#[test]
fn zero_generating_function_has_canonical_denominator_one() {
    let r = rule(&[(1, &[1])]);
    let gf = assert_expansion(&r, 4, &[row(&[0])], 30);
    assert!(gf.numerator.is_empty());
    assert_eq!(gf.denominator, vec![vec!["1"]]);
    assert!(riordan_decomposition(&gf).is_none());
}

#[test]
fn a245963_has_an_x_dependent_reduced_numerator() {
    let entry = polytool::oeis::by_id("A245963").unwrap();
    let (r, seeds) = entry.recurrence_parts().unwrap();
    let gf = assert_expansion(&r, entry.recurrence_first_index, &seeds, 30);
    assert_eq!(
        gf.numerator,
        vec![vec!["1"], vec!["0", "1"], vec!["0", "1"]]
    );
    assert!(riordan_decomposition(&gf).is_none());
}

#[test]
fn unsupported_rules_and_invalid_seeds_return_typed_errors() {
    let mut r = rule(&[(1, &[1])]);
    assert!(matches!(
        rational_generating_function(&r, 0, &[]),
        Err(RecurrenceError::InvalidInput(_))
    ));
    r.terms[0].deriv_order = 1;
    assert!(matches!(
        rational_generating_function(&r, 0, &[row(&[1])]),
        Err(RecurrenceError::Unsupported(_))
    ));
    r.terms[0].deriv_order = 0;
    r.terms[0].coeff = BivarPoly::from_coefficients_normalized(vec![row(&[1]), row(&[1])]);
    assert!(matches!(
        rational_generating_function(&r, 0, &[row(&[1])]),
        Err(RecurrenceError::Unsupported(_))
    ));
    r = rule(&[(1, &[1])]);
    r.denominator = Some(BivarPoly::from_coefficients_normalized(vec![row(&[0])]));
    assert!(matches!(
        rational_generating_function(&r, 0, &[row(&[1])]),
        Err(RecurrenceError::InvalidInput(_))
    ));
}

#[test]
fn every_supported_catalog_generating_function_matches_thirty_rows() {
    let mut checked = 0;
    for entry in polytool::oeis::catalog() {
        let (r, seeds) = entry.recurrence_parts().unwrap();
        if rational_generating_function(&r, entry.recurrence_first_index, &seeds).is_ok() {
            assert_expansion(
                &r,
                entry.recurrence_first_index,
                &seeds,
                30.max(seeds.len()),
            );
            checked += 1;
        }
    }
    assert!(checked >= 390, "{checked}");
}
