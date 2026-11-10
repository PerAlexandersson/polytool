//! All-index recurrence equality by an invariant module
//! (`certify_recurrence_equality`, module documentation "Invariant
//! module"): pairs that neither divides, from an audit of catalog
//! definitions (an order-two definition against its catalog rule with a
//! derivative, and row rules derived from hypergeometric triangles against
//! catalog rules of other orders), a left-hand factor that depends on `x`,
//! negative cases (rules whose rows agree on many rows and then differ),
//! replay and tampering.

use num_bigint::BigInt;
use polytool::oeis;
use polytool::recurrence::{
    certify_normalized_recurrence_equality, certify_recurrence_equality, BigRational, BivarPoly,
    NormalizedEqualityOptions, Recurrence, RecurrenceJson, RecurrenceTerm,
};

type Q = BigRational;

fn q(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

/// `coeffs[i][j]` is the coefficient of `n^i x^j`.
fn poly(coeffs: &[&[i64]]) -> BivarPoly {
    BivarPoly {
        coeffs: coeffs
            .iter()
            .map(|row| row.iter().map(|c| q(*c)).collect())
            .collect(),
    }
}

fn catalog_record(id: &str) -> RecurrenceJson {
    let definition = oeis::by_id(id).unwrap();
    let (rule, rows, first) = definition.recurrence_export_parts().unwrap();
    RecurrenceJson::from_recurrence_rational(&rule, first as usize, &rows, None)
}

fn record(json: &str) -> RecurrenceJson {
    serde_json::from_str(json).unwrap()
}

/// A394446 by its definition `(n + 2) p(n) = (2n + 1)(1 + 3x + x^2)
/// p(n - 1) - (n - 1)(1 + x + x^2)^2 p(n - 2)`.
const A394446_DEFINITION: &str = r#"{"schema": "polytool.recurrence.v1", "first_index": 0, "initial_polynomials": [["1"], ["1", "3", "1"]], "recurrence": {"terms": [{"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["1", "3", "1"], ["2", "6", "2"]]}}, {"offset": 2, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["1", "2", "3", "2", "1"], ["-1", "-2", "-3", "-2", "-1"]]}}], "denominator": {"coeffs": [["2"], ["1"]]}}}"#;

/// Row rules derived from the hypergeometric terms of the triangles
/// (first-order and second-order rules with derivatives).
const HYPERGEOMETRIC_RULES: [(&str, &str); 4] = [
    (
        "A110319",
        r#"{"schema": "polytool.recurrence.v1", "first_index": 2, "initial_polynomials": [["0", "0", "1"]], "recurrence": {"terms": [{"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["0"], ["4"], ["0"], ["0"], ["-4"]]}}, {"offset": 1, "deriv_order": 1, "sign": "none", "coeff": {"coeffs": [["0", "24", "-4"], ["0", "48", "-6"], ["0", "-80", "10"], ["0", "32", "0"]]}}, {"offset": 1, "deriv_order": 2, "sign": "none", "coeff": {"coeffs": [["0", "0", "96", "-8"], ["0", "0", "176", "-15"], ["0", "0", "-80", "5"]]}}, {"offset": 1, "deriv_order": 3, "sign": "none", "coeff": {"coeffs": [["0", "0", "0", "32", "-2"], ["0", "0", "0", "64", "-4"]]}}], "denominator": {"coeffs": [["0"], ["-2"], ["-1"], ["2"], ["1"]]}}}"#,
    ),
    (
        "A134434",
        r#"{"schema": "polytool.recurrence.v1", "first_index": 0, "initial_polynomials": [["1"], ["1"], ["1", "1"]], "recurrence": {"terms": [{"offset": 2, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["0", "4"], ["0", "36"], ["-15", "-169"], ["22", "186"], ["-15", "-61"], ["4", "-20"], ["4", "12"]]}}, {"offset": 2, "deriv_order": 1, "sign": "none", "coeff": {"coeffs": [["0", "-4", "28"], ["0", "-40", "12"], ["0", "144", "-232"], ["0", "-76", "188"], ["0", "-16", "-16"], ["0", "16", "-16"]]}}, {"offset": 2, "deriv_order": 2, "sign": "none", "coeff": {"coeffs": [["0", "0", "-12", "12"], ["0", "0", "8", "-8"], ["0", "0", "48", "-48"], ["0", "0", "-32", "32"]]}}, {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["-20", "-8"], ["-24", "0"], ["16", "24"], ["-32", "-16"]]}}, {"offset": 1, "deriv_order": 1, "sign": "none", "coeff": {"coeffs": [["0", "44", "-64"], ["0", "104", "-96"], ["0", "16", "64"], ["0", "-32", "0"]]}}, {"offset": 1, "deriv_order": 2, "sign": "none", "coeff": {"coeffs": [["0", "0", "32", "-32"], ["0", "0", "64", "-64"]]}}], "denominator": {"coeffs": [["-80"], ["64"], ["-44"], ["-16"], ["16"]]}}}"#,
    ),
    (
        "A174148",
        r#"{"schema": "polytool.recurrence.v1", "first_index": 1, "initial_polynomials": [["1", "1"], ["1", "12", "1"]], "recurrence": {"terms": [{"offset": 2, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["4", "-64", "28"], ["-18", "60", "-42"], ["-10", "202", "-88"], ["60", "-198", "138"], ["-30", "-78", "-12"], ["-6", "78", "-24"]]}}, {"offset": 2, "deriv_order": 1, "sign": "none", "coeff": {"coeffs": [["0", "-12", "0", "12"], ["0", "6", "0", "-6"], ["0", "42", "0", "-42"], ["0", "-18", "0", "18"], ["0", "-18", "0", "18"]]}}, {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["-6", "0"], ["23", "-4"], ["11", "2"], ["-75", "54"], ["27", "-90"], ["12", "30"]]}}, {"offset": 1, "deriv_order": 1, "sign": "none", "coeff": {"coeffs": [["0", "-6", "6"], ["0", "21", "-21"], ["0", "30", "-30"], ["0", "-99", "99"], ["0", "18", "-18"]]}}], "denominator": {"coeffs": [["-2"], ["5"], ["1"], ["-15"], ["-3"], ["6"]]}}}"#,
    ),
    (
        "A110205",
        r#"{"schema": "polytool.recurrence.v1", "first_index": 3, "initial_polynomials": [["73", "368", "343"], ["585", "3825", "6615", "3375"]], "recurrence": {"terms": [{"offset": 2, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["-152", "-1216", "8992", "-2184"], ["0", "0", "-2352", "896"], ["0", "0", "104", "-56"]]}}, {"offset": 2, "deriv_order": 1, "sign": "none", "coeff": {"coeffs": [["0", "0", "2456", "1504", "-952"], ["0", "0", "-208", "-96", "112"]]}}, {"offset": 2, "deriv_order": 2, "sign": "none", "coeff": {"coeffs": [["0", "0", "104", "152", "-8", "-56"]]}}, {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["171", "399", "-1134"], ["0", "57", "567"]]}}, {"offset": 1, "deriv_order": 1, "sign": "none", "coeff": {"coeffs": [["0", "-57", "-624", "-567"]]}}], "denominator": {"coeffs": [["19"]]}}}"#,
    ),
];

#[test]
fn a394446_definition_equals_its_catalog_rule() {
    // Neither rule right-divides the other; the rows satisfy a first-order
    // relation with a derivative on the leading row, which the module
    // closure finds.
    let catalog = catalog_record("A394446");
    let definition = record(A394446_DEFINITION);
    let certificate = certify_recurrence_equality(&catalog, &definition).unwrap();
    assert!(certificate.is_certified(), "{:?}", certificate.reasons);
    assert!(certificate.invariant_module.is_some());
    assert!(certificate
        .reasons
        .iter()
        .any(|reason| reason.contains("nonzero remainder")));
    certificate.verify().unwrap();
    // Either order of the records.
    let reversed = certify_recurrence_equality(&definition, &catalog).unwrap();
    assert!(reversed.is_certified(), "{:?}", reversed.reasons);
    reversed.verify().unwrap();

    // Tampering: every recorded field is recomputed.
    let mut forged = certificate.clone();
    forged.checked_bases.push(17);
    assert!(forged.verify().is_err());
    let mut forged = certificate.clone();
    forged.invariant_module = Some("the module generated by nothing".into());
    assert!(forged.verify().is_err());
    let mut forged = certificate.clone();
    forged.invariant_module = None;
    assert!(forged.verify().is_err());
    let mut forged = certificate;
    forged.alternative.initial_polynomials[1][1] = "4".into();
    assert!(forged.verify().is_err());
}

/// The derived rule of `id` is certified equal to the catalog rule (up to
/// the dropped leading rows of the normalization), with a replay.
fn hypergeometric_rule_equals_its_catalog_rule(id: &str) {
    let json = HYPERGEOMETRIC_RULES
        .iter()
        .find(|(name, _)| *name == id)
        .unwrap()
        .1;
    let certificate = certify_normalized_recurrence_equality(
        &record(json),
        &catalog_record(id),
        &NormalizedEqualityOptions::default(),
    )
    .unwrap();
    assert!(
        certificate.is_certified(),
        "{id}: {:?}",
        certificate.reasons
    );
    let equality = certificate.equality.as_ref().unwrap();
    assert!(equality.invariant_module.is_some(), "{id}");
    certificate.verify().unwrap();
}

#[test]
fn a110319_hypergeometric_rule_equals_its_catalog_rule() {
    hypergeometric_rule_equals_its_catalog_rule("A110319");
}

#[test]
fn a134434_hypergeometric_rule_equals_its_catalog_rule() {
    hypergeometric_rule_equals_its_catalog_rule("A134434");
}

#[test]
fn a174148_hypergeometric_rule_equals_its_catalog_rule() {
    hypergeometric_rule_equals_its_catalog_rule("A174148");
}

#[test]
fn a110205_hypergeometric_rule_equals_its_catalog_rule() {
    hypergeometric_rule_equals_its_catalog_rule("A110205");
}

/// `P_m = (1 + x) P_(m-1) + q(m) C(P_(m-1))` with `q(m) = (m - 1) (m - 2)
/// ... (m - 12)`, from the row `1` at index 0: the binomial rows wherever
/// `C` kills them, and otherwise rows that agree with them up to index 12.
/// `correction` lists `(derivative order, coefficient)` of `C`, the
/// coefficients as rows `[m power][x power]`.
fn perturbed_binomial(correction: &[(usize, &[&[i64]])]) -> RecurrenceJson {
    let mut q_m = vec![q(1)];
    for root in 1..=12 {
        let mut next = vec![q(0); q_m.len() + 1];
        for (i, c) in q_m.iter().enumerate() {
            next[i + 1] += c;
            next[i] -= c * q(root);
        }
        q_m = next;
    }
    let mut terms = vec![RecurrenceTerm::new(1, 0, poly(&[&[1, 1]]))];
    for (derivative, coefficient) in correction {
        let base = poly(coefficient);
        let mut coeffs = vec![vec![q(0); 2]; q_m.len() + base.coeffs.len()];
        for (i, c) in q_m.iter().enumerate() {
            for (a, row) in base.coeffs.iter().enumerate() {
                for (b, value) in row.iter().enumerate() {
                    coeffs[i + a][b] += c * value;
                }
            }
        }
        terms.push(RecurrenceTerm::new(
            1,
            *derivative,
            BivarPoly::from_coefficients_normalized(coeffs),
        ));
    }
    let rule = Recurrence {
        terms,
        denominator: None,
        inhomogeneous: None,
    };
    RecurrenceJson::from_recurrence_rational(&rule, 0, &[vec![q(1)]], None)
}

fn binomial_rows() -> RecurrenceJson {
    let rule = Recurrence {
        terms: vec![RecurrenceTerm::new(1, 0, poly(&[&[1, 1]]))],
        denominator: None,
        inhomogeneous: None,
    };
    RecurrenceJson::from_recurrence_rational(&rule, 0, &[vec![q(1)]], None)
}

#[test]
fn a_vanishing_correction_is_certified() {
    // C = (1 + x) D - (m - 1): (1 + x) P'_n = n P_n for the binomial rows,
    // so the correction vanishes on every row; only a row annihilator
    // shows it.
    let alternative = perturbed_binomial(&[(1, &[&[1, 1]]), (0, &[&[1], &[-1]])]);
    let certificate = certify_recurrence_equality(&binomial_rows(), &alternative).unwrap();
    assert!(certificate.is_certified(), "{:?}", certificate.reasons);
    certificate.verify().unwrap();
}

#[test]
fn rows_that_agree_on_twelve_rows_are_not_certified() {
    // C = x: the rows agree with (1 + x)^n up to index 12 and differ from
    // index 13 on.
    for correction in [
        vec![(0, &[&[0i64, 1][..]][..])],
        // C = x D - (m - 1): x P'_n - n P_n = -n (1 + x)^(n - 1) on the
        // binomial rows.
        vec![
            (1, &[&[0i64, 1][..]][..]),
            (0, &[&[1i64][..], &[-1][..]][..]),
        ],
    ] {
        let alternative = perturbed_binomial(&correction);
        let (rule, first, seeds) = alternative.to_recurrence_parts().unwrap();
        let rows = rule.generate_rows_rational(&seeds, first, 14).unwrap();
        let (binomial, _, _) = binomial_rows().to_recurrence_parts().unwrap();
        let expected = binomial.generate_rows_rational(&seeds, first, 14).unwrap();
        assert_eq!(rows[..13], expected[..13]);
        assert_ne!(rows[13], expected[13]);
        let certificate = certify_recurrence_equality(&binomial_rows(), &alternative).unwrap();
        assert!(!certificate.is_certified(), "{correction:?}");
        certificate.verify().unwrap();
        let reversed = certify_recurrence_equality(&alternative, &binomial_rows()).unwrap();
        assert!(!reversed.is_certified(), "{correction:?}");
    }
}

#[test]
fn a_left_hand_factor_with_x_is_supported() {
    // (1 + x) P_n = (1 + x)^2 P_(n-1): the binomial rows, with a left-hand
    // factor that depends on x and is nonzero at every row.
    let rule = Recurrence {
        terms: vec![RecurrenceTerm::new(1, 0, poly(&[&[1, 2, 1]]))],
        denominator: Some(poly(&[&[1, 1]])),
        inhomogeneous: None,
    };
    let alternative = RecurrenceJson::from_recurrence_rational(&rule, 0, &[vec![q(1)]], None);
    for (defining, other) in [
        (binomial_rows(), alternative.clone()),
        (alternative.clone(), binomial_rows()),
    ] {
        let certificate = certify_recurrence_equality(&defining, &other).unwrap();
        assert!(certificate.is_certified(), "{:?}", certificate.reasons);
        certificate.verify().unwrap();
    }
    // A left-hand factor that vanishes identically at a computed row is
    // still refused: (n - 2) (1 + x) at the row 2.
    let rule = Recurrence {
        terms: vec![RecurrenceTerm::new(1, 0, poly(&[&[-2, -2], &[1, 1]]))],
        denominator: Some(poly(&[&[-2], &[1]])),
        inhomogeneous: None,
    };
    let singular = RecurrenceJson::from_recurrence_rational(&rule, 0, &[vec![q(1)]], None);
    let certificate = certify_recurrence_equality(&binomial_rows(), &singular).unwrap();
    assert!(!certificate.is_certified());
}
