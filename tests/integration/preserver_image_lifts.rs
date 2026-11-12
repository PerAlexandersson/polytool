//! Preserver-image lifts: images of proved families under real-rootedness
//! preservers (multiplier sequences, Laguerre--Pólya differential operators,
//! Euler-operator polynomials, Rolle-type operators, coordinate maps and
//! hypergeometric multipliers), certified and proved; near misses, replay
//! and tampering, and the detection.
//!
//! The input rules are computed from the parent rule by
//! `image_recurrence_of_operator` (through
//! `PreserverImageCandidate::image_recurrence`) and checked against 30 exact
//! image rows computed directly.  One rule, `EULERIAN_PLUS_DERIVATIVE`, was
//! fitted to the image rows by the exact fixed-shape recurrence search,
//! independently of the lift machinery, and is kept as a cross-check.

use num_bigint::BigInt;
use polytool::is_real_rooted_bigint_coeffs;
use polytool::recurrence::{
    certify_preserver_image_lift, find_preserver_image_lifts, BigRational, BivarPoly,
    DegreeSchedule, HadamardFactor, HadamardFactorTerm, ImageRecurrenceOptions, ImageRowMap,
    IndexRationalFunction, InterlacingOrientation, LiftHypothesisStatus, LiftIdentityStatus,
    LiftOperation, LiftVerdict, PreserverImageCandidate, PreserverImageCertificate,
    PreserverImageCertification, PreserverImageIdentityStatus, PreserverImageOptions,
    PreserverImageReport, PreserverImageSearchOptions, PreserverInterlacingClaim,
    PreserverOperator, PreserverParentLibrary, PreserverParentRegion, Recurrence, RecurrenceError,
    RecurrenceJson, RecurrenceTerm, RowInterlacing, VeroneseSectionGauge,
};
use polytool::sequences;
use polytool::vec_poly;
use std::sync::{Arc, OnceLock};

type Q = BigRational;

const EULERIAN_PLUS_DERIVATIVE: &str = r#"
{"schema":"polytool.recurrence.v1","first_index":1,"initial_polynomials":[["1"],["2",
"1"]],"recurrence":{"terms":[{"offset":1,"deriv_order":0,"sign":"none",
"coeff":{"coeffs":[["7/2","-5","1/2","0","0"],["0","2","0","0","0"],["0","0","0","0",
"0"]]}},{"offset":1,"deriv_order":1,"sign":"none","coeff":{"coeffs":[["1/2","1","-3/2",
"0","0"],["0","0","0","0","0"],["0","0","0","0","0"]]}},{"offset":2,"deriv_order":0,
"sign":"none","coeff":{"coeffs":[["1/2","7","-17/2","1","0"],["-1","-5/2","6","-1/2",
"0"],["0","0","-1","0","0"]]}},{"offset":2,"deriv_order":1,"sign":"none",
"coeff":{"coeffs":[["-3/2","3/2","5","-11/2","1/2"],["0","-1/2","-1","3/2","0"],["0",
"0","0","0","0"]]}},{"offset":2,"deriv_order":2,"sign":"none","coeff":{"coeffs":[["0",
"-1/2","1/2","1/2","-1/2"],["0","0","0","0","0"],["0","0","0","0","0"]]}}]}}
"#;

use crate::common::q;

fn rational_rows(rows: &[Vec<BigInt>]) -> Vec<Vec<Q>> {
    rows.iter()
        .map(|row| vec_poly::trim(&row.iter().cloned().map(Q::from_integer).collect::<Vec<_>>()))
        .collect()
}

fn library() -> Arc<PreserverParentLibrary> {
    static LIBRARY: OnceLock<Arc<PreserverParentLibrary>> = OnceLock::new();
    LIBRARY
        .get_or_init(|| Arc::new(PreserverParentLibrary::named_families().unwrap()))
        .clone()
}

fn parent(name: &str) -> RecurrenceJson {
    library()
        .parents
        .iter()
        .find(|p| p.name == name)
        .unwrap_or_else(|| panic!("{name} is a library parent"))
        .recurrence
        .clone()
}

/// The rule of a fixture.
fn rule(fixture: &str) -> (Recurrence, usize) {
    let json: RecurrenceJson = serde_json::from_str(fixture).expect("fixture JSON");
    let (rule, first_index, _) = json.to_recurrence_parts().expect("fixture rule");
    (rule, first_index)
}

const ROWS: usize = 30;

fn eulerian() -> Vec<Vec<Q>> {
    rational_rows(&sequences::eulerian_polynomials_bigint(ROWS)[1..])
}

fn narayana() -> Vec<Vec<Q>> {
    rational_rows(&sequences::narayana_polynomials_bigint(ROWS)[1..])
}

fn hermite() -> Vec<Vec<Q>> {
    rational_rows(&sequences::hermite_polynomials_bigint(ROWS - 1))
}

fn binomial_rows() -> Vec<Vec<Q>> {
    (0..ROWS)
        .map(|m| {
            let mut row = vec![q(1)];
            for k in 1..=m {
                let next = &row[k - 1] * q((m + 1 - k) as i64) / q(k as i64);
                row.push(next);
            }
            row
        })
        .collect()
}

fn theta(row: &[Q], f: impl Fn(usize) -> Q) -> Vec<Q> {
    vec_poly::trim(
        &row.iter()
            .enumerate()
            .map(|(k, c)| c * f(k))
            .collect::<Vec<_>>(),
    )
}

fn plus_derivative(row: &[Q], order: usize) -> Vec<Q> {
    let mut d = row.to_vec();
    for _ in 0..order {
        d = vec_poly::derivative(&d);
    }
    vec_poly::trim(&vec_poly::add(row, &d))
}

fn candidate(parent_name: &str, operator: PreserverOperator) -> PreserverImageCandidate {
    PreserverImageCandidate {
        parent_name: parent_name.into(),
        parent: parent(parent_name),
        index_shift: 0,
        gauge: VeroneseSectionGauge::None,
        operator,
    }
}

fn certify(
    rule: &Recurrence,
    rows: &[Vec<Q>],
    first_index: usize,
    candidate: &PreserverImageCandidate,
) -> PreserverImageCertification {
    certify_with_options(
        rule,
        rows,
        first_index,
        candidate,
        &PreserverImageOptions::default(),
    )
}

fn certify_with_options(
    rule: &Recurrence,
    rows: &[Vec<Q>],
    first_index: usize,
    candidate: &PreserverImageCandidate,
    options: &PreserverImageOptions,
) -> PreserverImageCertification {
    let certification = certify_preserver_image_lift(rule, rows, first_index, candidate, options)
        .expect("the certification runs");
    certification.verify().expect("the certification replays");
    let json = serde_json::to_string(&certification).unwrap();
    let back: PreserverImageCertification = serde_json::from_str(&json).unwrap();
    assert_eq!(back, certification);
    certification
}

fn certify_fixture(
    fixture: &str,
    rows: &[Vec<Q>],
    candidate: &PreserverImageCandidate,
) -> PreserverImageCertification {
    let (rule, first_index) = rule(fixture);
    certify(&rule, rows, first_index, candidate)
}

/// The rule of the image of `candidate` from `first_index` on, computed
/// from the parent rule (not fitted), replayed against the exact rows.
fn computed_rule(
    candidate: &PreserverImageCandidate,
    first_index: usize,
    rows: &[Vec<Q>],
) -> Recurrence {
    let image = candidate
        .image_recurrence(first_index, &ImageRecurrenceOptions::default())
        .expect("the image rule is computed");
    image.verify().expect("the image rule replays");
    let generated = image.rows(rows.len()).expect("the image rule runs");
    for (offset, (row, expected)) in generated.iter().zip(rows).enumerate() {
        assert_eq!(&vec_poly::trim(row), expected, "row offset {offset}");
    }
    let (rule, first, _) = image.rule.to_recurrence_parts().unwrap();
    assert_eq!(first, first_index);
    rule
}

/// Certify the computed rule of the image, also asking for the image rule
/// on the certification.
fn certify_computed(
    rows: &[Vec<Q>],
    first_index: usize,
    candidate: &PreserverImageCandidate,
) -> PreserverImageCertification {
    let rule = computed_rule(candidate, first_index, rows);
    let options = PreserverImageOptions {
        image_rule: Some(ImageRecurrenceOptions::default()),
        ..PreserverImageOptions::default()
    };
    let certification = certify_with_options(&rule, rows, first_index, candidate, &options);
    if matches!(
        certification.identity,
        PreserverImageIdentityStatus::Certified { .. }
    ) {
        let image = certification.image_rule.as_ref().expect("an image rule");
        let mut tampered = certification.clone();
        let mut wrong = image.clone();
        wrong.certificate.first_index += 1;
        tampered.image_rule = Some(wrong);
        assert!(tampered.verify().is_err());
    }
    certification
}

fn assert_proved(certification: &PreserverImageCertification, region: &str) {
    assert!(certification.is_proved(), "{certification:#?}");
    let summary = &certification.lift_summary;
    assert_eq!(summary.verdict, LiftVerdict::Proved);
    assert_eq!(summary.root_transfer.input_root_region.to_string(), region);
    assert!(matches!(
        summary.identity,
        LiftIdentityStatus::Certified { .. }
    ));
    assert!(matches!(
        summary.operation,
        LiftOperation::PreserverImage { .. }
    ));
    assert!(summary
        .hypotheses
        .iter()
        .all(|h| h.status == LiftHypothesisStatus::Verified));
}

fn division_certificate(
    certification: &PreserverImageCertification,
) -> &polytool::recurrence::PreserverDivisionCertificate {
    match &certification.identity {
        PreserverImageIdentityStatus::Certified { certificate } => match certificate.as_ref() {
            PreserverImageCertificate::OperatorDivision(certificate) => certificate,
            other => panic!("not a division certificate: {other:?}"),
        },
        other => panic!("not certified: {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// The menu, one shape at a time
// ---------------------------------------------------------------------------

/// `A_n + A_n'`: the Laguerre operator `1 + D` on the Eulerian polynomials.
#[test]
fn eulerian_plus_derivative_is_proved() {
    let rows: Vec<_> = eulerian().iter().map(|p| plus_derivative(p, 1)).collect();
    let lift = candidate(
        "eulerian",
        PreserverOperator::derivative_polynomial(&[q(1), q(1)]),
    );
    // The fitted rule and the computed one are both certified.
    let certification = certify_fixture(EULERIAN_PLUS_DERIVATIVE, &rows, &lift);
    assert_proved(&certification, "(-inf, 0]");
    assert!(certification
        .lean_declarations
        .contains(&"RealRooted.splits_add_C_mul_derivative_all".to_string()));
    let computed = certify_computed(&rows, 1, &lift);
    assert_proved(&computed, "(-inf, 0]");
    assert_eq!(
        computed.image_rule.as_ref().map(|image| image.order),
        Some(2)
    );
}

/// `theta (theta + 1) N_n`: an Euler-operator polynomial with roots `0, -1`
/// on the Narayana polynomials.
#[test]
fn euler_operator_polynomial_on_narayana_is_proved() {
    let rows: Vec<_> = narayana()
        .iter()
        .map(|p| theta(p, |k| q((k * (k + 1)) as i64)))
        .collect();
    let certification = certify_computed(
        &rows,
        1,
        &candidate(
            "narayana",
            PreserverOperator::euler_polynomial(&[q(0), q(1), q(1)]),
        ),
    );
    assert_proved(&certification, "(-inf, 0]");
}

/// `(theta + n) N_n`: an index-dependent multiplier `k + n`, decided by the
/// sign checker.
#[test]
fn index_dependent_euler_operator_is_proved() {
    let rows: Vec<_> = narayana()
        .iter()
        .enumerate()
        .map(|(offset, p)| theta(p, |k| q((k + offset + 1) as i64)))
        .collect();
    let operator = PreserverOperator::EulerPolynomial {
        coefficients: vec![
            IndexRationalFunction::from_polynomials(&[q(0), q(1)], &[q(1)]),
            IndexRationalFunction::constant(&q(1)),
        ],
    };
    let certification = certify_computed(&rows, 1, &candidate("narayana", operator));
    assert_proved(&certification, "(-inf, 0]");
    assert!(certification
        .lift_summary
        .hypotheses
        .iter()
        .any(|h| h.name == "euler_polynomial_coefficients_01_same_sign"
            && h.sign_decision.as_ref().is_some_and(|d| d.holds())));
}

/// `D(x A_n)` (Rolle) and `A_n(2x + 1)` (affine, roots in `(-inf, -1/2]`).
#[test]
fn rolle_and_affine_images_of_eulerian_are_proved() {
    let rows: Vec<_> = eulerian()
        .iter()
        .map(|p| {
            let mut shifted = vec![q(0)];
            shifted.extend_from_slice(p);
            vec_poly::trim(&vec_poly::derivative(&shifted))
        })
        .collect();
    let certification = certify_computed(
        &rows,
        1,
        &candidate(
            "eulerian",
            PreserverOperator::DerivativeOfPowerMultiple {
                power: 1,
                derivative_order: 1,
            },
        ),
    );
    assert_proved(&certification, "(-inf, 0]");
    assert_eq!(certification.lean_route.as_deref(), Some("derivativeXMul"));

    let rows: Vec<_> = eulerian()
        .iter()
        .map(|p| vec_poly::trim(&vec_poly::compose(p, &[q(1), q(2)])))
        .collect();
    let certification = certify_computed(
        &rows,
        1,
        &candidate(
            "eulerian",
            PreserverOperator::Affine {
                slope: "2".into(),
                intercept: "1".into(),
            },
        ),
    );
    assert_proved(&certification, "(-inf, -1/2]");
}

/// The `f`-polynomials of the Eulerian polynomials,
/// `(1 + x)^(n-1) A_n(x/(1 + x))` (the Fubini triangle), roots in `[-1, 0]`.
#[test]
fn f_polynomial_of_eulerian_is_proved() {
    let rows: Vec<_> = eulerian()
        .iter()
        .enumerate()
        .map(|(d, p)| {
            let mut out = vec![q(0)];
            for (k, a) in p.iter().enumerate() {
                let mut term = vec![a.clone()];
                for _ in 0..k {
                    term = vec_poly::mul(&term, &[q(0), q(1)]);
                }
                for _ in k..d {
                    term = vec_poly::mul(&term, &[q(1), q(1)]);
                }
                out = vec_poly::add(&out, &term);
            }
            vec_poly::trim(&out)
        })
        .collect();
    let certification = certify_computed(
        &rows,
        1,
        &candidate(
            "eulerian",
            PreserverOperator::f_polynomial(DegreeSchedule::new(1, 0, 1)),
        ),
    );
    assert_proved(&certification, "[-1, 0]");
    assert_eq!(certification.lean_route.as_deref(), Some("fTransform"));
}

/// `sum_k binom(n, k) x^k / k!` (Laguerre polynomials `L_n(-x)`): the
/// multiplier `1/k!` through the Hadamard conjugation.
#[test]
fn inverse_factorial_multiplier_of_the_binomial_rows_is_proved() {
    let mut factorial = vec![q(1)];
    for k in 1..=ROWS {
        let next = &factorial[k - 1] * q(k as i64);
        factorial.push(next);
    }
    let rows: Vec<_> = binomial_rows()
        .iter()
        .map(|p| theta(p, |k| factorial[k].recip()))
        .collect();
    let certification = certify_computed(
        &rows,
        0,
        &candidate(
            "binomial",
            PreserverOperator::HadamardMultiplier {
                factor: HadamardFactor::single(HadamardFactorTerm::InverseFactorial),
            },
        ),
    );
    assert_proved(&certification, "(-inf, 0]");
    let PreserverImageIdentityStatus::Certified { certificate } = &certification.identity else {
        unreachable!()
    };
    assert!(matches!(
        **certificate,
        PreserverImageCertificate::HadamardConjugation(_)
    ));
    assert!(certification
        .lean_declarations
        .contains(&"RealRooted.isMultiplierSequence_inv_factorial".to_string()));
}

/// `He_n + He_n'`: a parent real-rooted on the whole line (three-term
/// theorem), and a division remainder removed by the Hermite differential
/// equation (an invariant row annihilator).
#[test]
fn hermite_plus_derivative_is_proved_on_the_whole_line() {
    let rows: Vec<_> = hermite().iter().map(|p| plus_derivative(p, 1)).collect();
    let certification = certify_computed(
        &rows,
        0,
        &candidate(
            "hermite_he",
            PreserverOperator::derivative_polynomial(&[q(1), q(1)]),
        ),
    );
    assert_proved(&certification, "(-inf, +inf)");
    assert!(certification.whole_line_parent_proof.is_some());
    assert!(division_certificate(&certification)
        .row_annihilation
        .is_some());
}

/// The reversal of the Chebyshev polynomials, `x^n T_n(1/x)`, with the
/// rule `R_n = 2 R_(n-1) - x^2 R_(n-2)` written by hand.
#[test]
fn reversal_of_chebyshev_is_proved() {
    let rows: Vec<_> = rational_rows(&sequences::chebyshev_polynomials_t_bigint(ROWS - 1))
        .iter()
        .enumerate()
        .map(|(n, p)| vec_poly::trim(&vec_poly::reverse(p, n).unwrap()))
        .collect();
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(
                1,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![q(2)]]),
            ),
            RecurrenceTerm::new(
                2,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(0), q(-1)]]),
            ),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let certification = certify(
        &rule,
        &rows,
        0,
        &candidate(
            "chebyshev_t",
            PreserverOperator::Reversal {
                schedule: DegreeSchedule::new(0, 0, 1),
            },
        ),
    );
    assert_proved(&certification, "(-inf, +inf)");
    assert_eq!(certification.lean_route.as_deref(), Some("reversed"));
}

/// A gauge `lambda(n) = 2^n` and an index shift: `2^n (A_(n+1) + A_(n+1)')`
/// for `n >= 0`, whose rule is the fixture's with `n -> n + 1` and each
/// term at lag `i` multiplied by `2^i`.
#[test]
fn gauge_and_index_shift_are_certified() {
    let (fixture, _) = rule(EULERIAN_PLUS_DERIVATIVE);
    let shifted = fixture.shift_index(1);
    let scaled = Recurrence {
        terms: shifted
            .terms
            .iter()
            .map(|term| {
                let factor = q(1i64 << term.offset);
                let coeffs = term
                    .coeff
                    .coeffs
                    .iter()
                    .map(|row| row.iter().map(|c| c * &factor).collect())
                    .collect();
                RecurrenceTerm::new(
                    term.offset,
                    term.deriv_order,
                    BivarPoly::from_coefficients_normalized(coeffs),
                )
            })
            .collect(),
        denominator: shifted.denominator.clone(),
        inhomogeneous: None,
    };
    let rows: Vec<_> = eulerian()
        .iter()
        .enumerate()
        .map(|(n, p)| vec_poly::scale(&plus_derivative(p, 1), &q(1i64 << n)))
        .collect();
    let mut gauged = candidate(
        "eulerian",
        PreserverOperator::derivative_polynomial(&[q(1), q(1)]),
    );
    gauged.index_shift = 1;
    gauged.gauge = VeroneseSectionGauge::Ratio {
        numerator: vec!["2".into()],
        denominator: vec!["1".into()],
        initial_value: "1".into(),
    };
    let certification = certify(&scaled, &rows, 0, &gauged);
    assert_proved(&certification, "(-inf, 0]");
}

/// The operators convert to image row maps with the same action; a
/// multiplier that depends on `n` does not.
#[test]
fn operators_convert_to_image_row_maps() {
    let inverse_factorial = PreserverOperator::HadamardMultiplier {
        factor: HadamardFactor::single(HadamardFactorTerm::InverseFactorial),
    };
    assert_eq!(
        ImageRowMap::try_from(&inverse_factorial).unwrap(),
        ImageRowMap::CoefficientMultiplier {
            numerator: vec!["1".into()],
            denominator: vec!["1".into(), "1".into()],
        }
    );
    let binomial = PreserverOperator::HadamardMultiplier {
        factor: HadamardFactor::single(HadamardFactorTerm::Binomial {
            index_offset: 0,
            coefficient_offset: 0,
        }),
    };
    assert!(matches!(
        ImageRowMap::try_from(&binomial),
        Err(RecurrenceError::Unsupported(_))
    ));
    // phi(D) = 1 + D/2 as a Weyl operator.
    let ImageRowMap::WeylOperator { terms, denominator } =
        ImageRowMap::try_from(&PreserverOperator::derivative_polynomial(&[q(2), q(1)])).unwrap()
    else {
        panic!("a Weyl operator")
    };
    assert_eq!(terms.len(), 2);
    assert_eq!(denominator, vec!["1".to_string()]);
    // A candidate with a multiplier that depends on n gets no image rule.
    let mut lift = candidate("binomial", binomial);
    assert!(lift.image_operator_specification().is_err());
    lift.operator = inverse_factorial;
    let specification = lift.image_operator_specification().unwrap();
    assert_eq!(specification.index_shift, 0);
}

// ---------------------------------------------------------------------------
// Near misses
// ---------------------------------------------------------------------------

fn integer_row(row: &[Q]) -> Vec<BigInt> {
    row.iter().map(|c| c.to_integer()).collect()
}

/// `(theta - 2) N_n`: a non-preserver whose rows are real-rooted on every
/// row checked.  The identity is certified, but the multiplier `k - 2` has
/// the positive root 2 (`k -> k - 2` is not a multiplier sequence), so the
/// lift fails; it is never proved.
#[test]
fn a_multiplier_with_a_positive_root_is_not_proved() {
    let rows: Vec<_> = narayana()
        .iter()
        .map(|p| theta(p, |k| q(k as i64 - 2)))
        .collect();
    assert!(rows
        .iter()
        .all(|row| is_real_rooted_bigint_coeffs(&integer_row(row))));
    let certification = certify_computed(
        &rows,
        1,
        &candidate(
            "narayana",
            PreserverOperator::euler_polynomial(&[q(-2), q(1)]),
        ),
    );
    assert!(!certification.is_proved());
    assert_eq!(certification.lift_summary.verdict, LiftVerdict::RouteFailed);
    assert!(matches!(
        certification.lift_summary.identity,
        LiftIdentityStatus::Certified { .. }
    ));
    assert!(certification.lift_summary.hypotheses.iter().any(|h| {
        h.name == "euler_polynomial_roots_nonpositive"
            && matches!(h.status, LiftHypothesisStatus::Failed { .. })
    }));
}

/// `He_n + He_n''`: `1 + D^2` is not a Laguerre--Pólya operator (and the
/// rows from `n = 2` on are not real-rooted); the identity is certified, the
/// side condition fails and the lift is not proved.
#[test]
fn a_non_preserver_is_not_proved() {
    let rows: Vec<_> = hermite().iter().map(|p| plus_derivative(p, 2)).collect();
    assert!(!is_real_rooted_bigint_coeffs(&integer_row(&rows[2])));
    let certification = certify_computed(
        &rows,
        0,
        &candidate(
            "hermite_he",
            PreserverOperator::derivative_polynomial(&[q(1), q(0), q(1)]),
        ),
    );
    assert!(!certification.is_proved());
    assert_eq!(certification.lift_summary.verdict, LiftVerdict::RouteFailed);
    // The same operator as a general Weyl operator has no theorem: open.
    let general = PreserverOperator::WeylOperator {
        terms: serde_json::from_str(
            r#"[{"deriv_order":0,"coeff":{"coeffs":[["1"]]}},
                {"deriv_order":2,"coeff":{"coeffs":[["1"]]}}]"#,
        )
        .unwrap(),
        denominator: vec!["1".into()],
    };
    let certification = certify_computed(&rows, 0, &candidate("hermite_he", general));
    assert!(matches!(
        certification.lift_summary.identity,
        LiftIdentityStatus::Certified { .. }
    ));
    assert_eq!(certification.lift_summary.verdict, LiftVerdict::Open);
}

/// Wrong parameters are refuted at an index: `A_n(3x + 1)` against the rows
/// of `A_n(2x + 1)` differs at `n = 2` (`A_1 = 1` agrees).
#[test]
fn wrong_parameters_are_refuted_at_an_index() {
    let rows: Vec<_> = eulerian()
        .iter()
        .map(|p| vec_poly::trim(&vec_poly::compose(p, &[q(1), q(2)])))
        .collect();
    let affine = |slope: &str| {
        candidate(
            "eulerian",
            PreserverOperator::Affine {
                slope: slope.into(),
                intercept: "1".into(),
            },
        )
    };
    let rule = computed_rule(&affine("2"), 1, &rows);
    let certification = certify(&rule, &rows, 1, &affine("3"));
    assert!(matches!(
        certification.identity,
        PreserverImageIdentityStatus::Refuted { index: 2, .. }
    ));
    assert_eq!(certification.lift_summary.verdict, LiftVerdict::RouteFailed);
}

/// An operator coefficient `c(n) = (n - 2)/(n - 2)` is undefined at
/// `n = 2`: the candidate row is undefined there (refuted) and the
/// hypothesis `operator_coefficients_defined` fails with witness 2.
#[test]
fn a_vanishing_coefficient_denominator_is_refuted() {
    let rows: Vec<_> = eulerian().iter().map(|p| plus_derivative(p, 1)).collect();
    let operator = PreserverOperator::DerivativePolynomial {
        coefficients: vec![
            IndexRationalFunction::constant(&q(1)),
            IndexRationalFunction::from_polynomials(&[q(-2), q(1)], &[q(-2), q(1)]),
        ],
    };
    let certification = certify_fixture(
        EULERIAN_PLUS_DERIVATIVE,
        &rows,
        &candidate("eulerian", operator),
    );
    assert!(matches!(
        certification.identity,
        PreserverImageIdentityStatus::Refuted { index: 2, .. }
    ));
    let hypothesis = certification
        .lift_summary
        .hypotheses
        .iter()
        .find(|h| h.name == "operator_coefficients_defined")
        .unwrap();
    assert_eq!(
        hypothesis.status,
        LiftHypothesisStatus::Failed {
            witness: Some(BigInt::from(2))
        }
    );
}

// ---------------------------------------------------------------------------
// Replay and tampering
// ---------------------------------------------------------------------------

#[test]
fn tampered_certificates_and_records_are_rejected() {
    let rows: Vec<_> = eulerian().iter().map(|p| plus_derivative(p, 1)).collect();
    let certification = certify_fixture(
        EULERIAN_PLUS_DERIVATIVE,
        &rows,
        &candidate(
            "eulerian",
            PreserverOperator::derivative_polynomial(&[q(1), q(1)]),
        ),
    );
    let certificate = division_certificate(&certification).clone();
    certificate.verify().unwrap();
    let mut tampered = certificate.clone();
    tampered.multiplier = vec!["2".into()];
    assert!(tampered.verify().is_err());
    let mut tampered = certificate.clone();
    tampered.operator = PreserverOperator::derivative_polynomial(&[q(1), q(2)]);
    assert!(tampered.verify().is_err());
    let mut tampered = certificate.clone();
    tampered.exceptional_bases.clear();
    tampered.quotients.pop();
    assert!(tampered.verify().is_err());
    // Tampered rows in the records: an identifying input row, a parent seed.
    let mut tampered = certificate.clone();
    tampered.input.initial_polynomials[1][0] = "3".into();
    assert!(tampered.verify().is_err());
    let mut tampered = certificate.clone();
    tampered.parent.initial_polynomials[0][0] = "2".into();
    assert!(tampered.verify().is_err());
    // A tampered summary or candidate.
    let mut tampered = certification.clone();
    tampered.lift_summary.verdict = LiftVerdict::Open;
    assert!(tampered.verify().is_err());
    let mut tampered = certification.clone();
    tampered.candidate.index_shift = 1;
    assert!(tampered.verify().is_err());
}

// ---------------------------------------------------------------------------
// Detection
// ---------------------------------------------------------------------------

fn search(rule: &Recurrence, first_index: usize, rows: &[Vec<Q>]) -> PreserverImageReport {
    let report = find_preserver_image_lifts(
        rule,
        rows,
        first_index,
        &PreserverImageSearchOptions {
            library: Some(library()),
            ..PreserverImageSearchOptions::default()
        },
    )
    .expect("the search runs");
    report.verify().expect("the report replays");
    let json = serde_json::to_string(&report).unwrap();
    let back: PreserverImageReport = serde_json::from_str(&json).unwrap();
    assert_eq!(back, report);
    report
}

#[test]
fn the_detection_finds_the_operator_and_the_parent() {
    let rows: Vec<_> = narayana()
        .iter()
        .enumerate()
        .map(|(offset, p)| theta(p, |k| q((k + offset + 1) as i64)))
        .collect();
    let operator = PreserverOperator::EulerPolynomial {
        coefficients: vec![
            IndexRationalFunction::from_polynomials(&[q(0), q(1)], &[q(1)]),
            IndexRationalFunction::constant(&q(1)),
        ],
    };
    let rule = computed_rule(&candidate("narayana", operator), 1, &rows);
    let report = search(&rule, 1, &rows);
    let lift = report.proved_lifts().next().expect("a proved lift");
    assert_eq!(lift.candidate.parent_name, "narayana");
    assert_eq!(
        lift.candidate.operator.text(),
        "f(n, theta) = n + (1) theta"
    );

    let rows: Vec<_> = eulerian()
        .iter()
        .map(|p| vec_poly::trim(&vec_poly::compose(p, &[q(1), q(2)])))
        .collect();
    let affine = PreserverOperator::Affine {
        slope: "2".into(),
        intercept: "1".into(),
    };
    let rule = computed_rule(&candidate("eulerian", affine), 1, &rows);
    let report = search(&rule, 1, &rows);
    let lift = report.proved_lifts().next().expect("a proved lift");
    assert_eq!(lift.candidate.parent_name, "eulerian");
    assert_eq!(
        lift.candidate.operator,
        PreserverOperator::Affine {
            slope: "2".into(),
            intercept: "1".into()
        }
    );
}

#[test]
fn the_detection_proves_nothing_for_a_non_preserver() {
    let rows: Vec<_> = hermite().iter().map(|p| plus_derivative(p, 2)).collect();
    let operator = PreserverOperator::derivative_polynomial(&[q(1), q(0), q(1)]);
    let rule = computed_rule(&candidate("hermite_he", operator), 0, &rows);
    let report = search(&rule, 0, &rows);
    assert_eq!(report.proved_lifts().count(), 0, "{report:#?}");
}

#[test]
fn the_named_library_is_proved() {
    let library = library();
    for parent in &library.parents {
        assert_ne!(
            parent.region,
            PreserverParentRegion::Unproved,
            "{}",
            parent.name
        );
    }
    assert!(library.parents.len() >= 9);
    // The peak polynomials are proved by a gamma lift, not by the parent
    // provers, so they are not a library parent.
    assert!(library.skipped.iter().any(|(name, _)| name == "peak"));
}

// ---------------------------------------------------------------------------
// Interlacing claims
// ---------------------------------------------------------------------------

/// `He_n + He_n'` and `A_n + A_n'`: a fixed Laguerre operator transfers the
/// parent's proved interlacing (Hermite--Kakeya--Obreschkoff), weakly and
/// without an orientation.
#[test]
fn a_fixed_preserver_transfers_the_parent_interlacing() {
    let rows: Vec<_> = hermite().iter().map(|p| plus_derivative(p, 1)).collect();
    let certification = certify_computed(
        &rows,
        0,
        &candidate(
            "hermite_he",
            PreserverOperator::derivative_polynomial(&[q(1), q(1)]),
        ),
    );
    let PreserverInterlacingClaim::Claimed {
        kind,
        orientation,
        from_index,
        ..
    } = &certification.interlacing
    else {
        panic!("{:?}", certification.interlacing)
    };
    assert_eq!(*kind, RowInterlacing::Weak);
    assert_eq!(orientation, "unordered");
    assert_eq!(*from_index, 1);
    // The summary carries the same claim.
    let summary_claim = certification
        .lift_summary
        .input_interlacing
        .proved()
        .expect("the summary claim");
    assert_eq!(summary_claim.kind, RowInterlacing::Weak);
    assert_eq!(summary_claim.orientation, InterlacingOrientation::Unordered);
    assert_eq!(summary_claim.from_index, 1);
    // A tampered claim is rejected.
    let mut tampered = certification.clone();
    tampered.interlacing = PreserverInterlacingClaim::None { reason: "x".into() };
    assert!(tampered.verify().is_err());

    let rows: Vec<_> = eulerian().iter().map(|p| plus_derivative(p, 1)).collect();
    let certification = certify_fixture(
        EULERIAN_PLUS_DERIVATIVE,
        &rows,
        &candidate(
            "eulerian",
            PreserverOperator::derivative_polynomial(&[q(1), q(1)]),
        ),
    );
    assert!(
        certification.interlacing.is_claimed(),
        "{:?}",
        certification.interlacing
    );
}

/// An index-dependent operator or schedule gets no interlacing claim.
#[test]
fn index_dependent_operators_claim_no_interlacing() {
    let rows: Vec<_> = narayana()
        .iter()
        .enumerate()
        .map(|(offset, p)| theta(p, |k| q((k + offset + 1) as i64)))
        .collect();
    let operator = PreserverOperator::EulerPolynomial {
        coefficients: vec![
            IndexRationalFunction::from_polynomials(&[q(0), q(1)], &[q(1)]),
            IndexRationalFunction::constant(&q(1)),
        ],
    };
    let certification = certify_computed(&rows, 1, &candidate("narayana", operator));
    assert!(certification.is_proved());
    assert!(!certification.interlacing.is_claimed());
    assert!(!certification.lift_summary.input_interlacing.is_proved());
}

// ---------------------------------------------------------------------------
// Interval parents
// ---------------------------------------------------------------------------

mod interval_parents {
    use super::*;
    use polytool::recurrence::{LiftParentProof, VeroneseParentProofStatus};

    fn integer_rule(terms: &[(usize, usize, Vec<Vec<i64>>)]) -> Recurrence {
        Recurrence {
            terms: terms
                .iter()
                .map(|(offset, order, coefficients)| {
                    RecurrenceTerm::new(
                        *offset,
                        *order,
                        BivarPoly::from_coefficients_normalized(
                            coefficients
                                .iter()
                                .map(|row| row.iter().map(|&v| q(v)).collect())
                                .collect(),
                        ),
                    )
                })
                .collect(),
            denominator: None,
            inhomogeneous: None,
        }
    }

    /// The rows `P_1, ..., P_ROWS` of a first-order rule from `P_1 = 1`.
    fn rows_from_one(rule: &Recurrence) -> Vec<Vec<Q>> {
        rule.generate_rows_rational(&[vec![q(1)]], 1, ROWS)
            .expect("rows")
            .iter()
            .map(|row| vec_poly::trim(row))
            .collect()
    }

    fn interval_candidate(
        parent_rule: &Recurrence,
        operator: PreserverOperator,
    ) -> PreserverImageCandidate {
        PreserverImageCandidate {
            parent_name: "interval_parent".into(),
            parent: RecurrenceJson::from_recurrence_rational(parent_rule, 1, &[vec![q(1)]], None),
            index_shift: 0,
            gauge: VeroneseSectionGauge::None,
            operator,
        }
    }

    fn parent_route(certification: &PreserverImageCertification) -> String {
        match &certification.lift_summary.parent_proof {
            LiftParentProof::Proved { route, .. } => route.clone(),
            other => panic!("the parent is not proved: {other:?}"),
        }
    }

    /// The `f`-polynomials of the Eulerian polynomials (A019538, the rows
    /// `k! S(n, k)`), `P_n = (1 + 2x) P_(n-1) + x (1 + x) P'_(n-1)`, lie in
    /// `(-1, 0)` by the Ma--Wang interval check, so `P_n(2x + 1)` lies in
    /// `(-1, -1/2)`, strictly inside the ray image `(-inf, -1/2]`.
    #[test]
    fn a_ma_wang_parent_bounds_an_affine_image() {
        let parent_rule = integer_rule(&[(1, 0, vec![vec![1, 2]]), (1, 1, vec![vec![0, 1, 1]])]);
        // Q_n = (3 + 4x) Q_(n-1) + (1 + 3x + 2x^2) Q'_(n-1).
        let input = integer_rule(&[(1, 0, vec![vec![3, 4]]), (1, 1, vec![vec![1, 3, 2]])]);
        let rows = rows_from_one(&input);
        for (parent_row, row) in rows_from_one(&parent_rule).iter().zip(&rows) {
            assert_eq!(
                &vec_poly::trim(&vec_poly::compose(parent_row, &[q(1), q(2)])),
                row
            );
        }
        let certification = certify(
            &input,
            &rows,
            1,
            &interval_candidate(
                &parent_rule,
                PreserverOperator::Affine {
                    slope: "2".into(),
                    intercept: "1".into(),
                },
            ),
        );
        assert_proved(&certification, "(-1, -1/2)");
        assert_eq!(parent_route(&certification), "ma_wang_interval");
        assert_eq!(
            certification
                .lift_summary
                .root_transfer
                .parent_root_region
                .as_ref()
                .map(|r| r.to_string())
                .as_deref(),
            Some("(-1, 0)")
        );
    }

    /// The reflection `G_n(x) = P_n(-x)` lies in `(0, 1)`: its rows have
    /// coefficients of both signs, so the ray statement of the parent
    /// provers is refuted, and only the interval route proves the image
    /// `G_n(x - 1) = P_n(1 - x)`, rooted in `(1, 2)`.
    #[test]
    fn only_the_interval_route_proves_an_image_of_a_positive_interval() {
        let parent_rule = integer_rule(&[(1, 0, vec![vec![1, -2]]), (1, 1, vec![vec![0, 1, -1]])]);
        // Q_n = (3 - 2x) Q_(n-1) + (-2 + 3x - x^2) Q'_(n-1).
        let input = integer_rule(&[(1, 0, vec![vec![3, -2]]), (1, 1, vec![vec![-2, 3, -1]])]);
        let rows = rows_from_one(&input);
        for (parent_row, row) in rows_from_one(&parent_rule).iter().zip(&rows) {
            assert_eq!(
                &vec_poly::trim(&vec_poly::compose(parent_row, &[q(-1), q(1)])),
                row
            );
        }
        let candidate = interval_candidate(
            &parent_rule,
            PreserverOperator::Affine {
                slope: "1".into(),
                intercept: "-1".into(),
            },
        );
        let certification = certify(&input, &rows, 1, &candidate);
        assert_eq!(
            certification.parent_proof.status,
            VeroneseParentProofStatus::Refuted
        );
        assert_proved(&certification, "(1, 2)");
        assert_eq!(parent_route(&certification), "ma_wang_interval");
        // The parent's weak interlacing (Ma--Wang) transfers through the
        // fixed affine map.
        assert!(certification.interlacing.is_claimed());
        // Removing or replacing the interval report is detected.
        let mut tampered = certification.clone();
        tampered.interval_parent_proof = None;
        assert!(tampered.verify().is_err());
        let mut tampered = certification.clone();
        let mut report = tampered.interval_parent_proof.clone().unwrap();
        report.recurrence =
            RecurrenceJson::from_recurrence_rational(&input, 1, &[vec![q(1)]], None);
        tampered.interval_parent_proof = Some(report);
        assert!(tampered.verify().is_err());
        // Without the interval report the lift is not proved.
        let summary = {
            let mut stripped = certification.clone();
            stripped.interval_parent_proof = None;
            stripped.recompute_lift_summary()
        };
        assert_ne!(summary.verdict, LiftVerdict::Proved);
        // Rows check: every row of the input is real-rooted in (1, 2).
        for row in &rows[1..] {
            assert!(is_real_rooted_bigint_coeffs(&integer_row(row)));
            // Every root r > 1 (Q(1 - z) has the negative roots 1 - r) and
            // r < 2 (Q(2 + z) has the negative roots r - 2).
            let shifted_left = vec_poly::compose(row, &[q(1), q(-1)]);
            let shifted_right = vec_poly::compose(row, &[q(2), q(1)]);
            for shifted in [shifted_left, shifted_right] {
                assert!(
                    shifted.iter().all(|c| c > &q(0)) || shifted.iter().all(|c| c < &q(0)),
                    "{row:?}"
                );
            }
        }
    }

    /// The shifted Chebyshev polynomials `T_n(2x - 1)` from the three-term
    /// interval `(-1, 1)` of the Chebyshev parent: the region `(0, 1)`
    /// instead of the line.
    #[test]
    fn a_three_term_interval_bounds_the_shifted_chebyshev_rows() {
        let rows: Vec<_> = rational_rows(&sequences::chebyshev_polynomials_t_bigint(ROWS - 1))
            .iter()
            .map(|p| vec_poly::trim(&vec_poly::compose(p, &[q(-1), q(2)])))
            .collect();
        // Q_n = (4x - 2) Q_(n-1) - Q_(n-2).
        let input = integer_rule(&[(1, 0, vec![vec![-2, 4]]), (2, 0, vec![vec![-1]])]);
        let certification = certify(
            &input,
            &rows,
            0,
            &candidate(
                "chebyshev_t",
                PreserverOperator::Affine {
                    slope: "2".into(),
                    intercept: "-1".into(),
                },
            ),
        );
        assert_proved(&certification, "(0, 1)");
        assert_eq!(parent_route(&certification), "three_term_interval");
        // The detection finds a lift with a bounded region.
        let report = find_preserver_image_lifts(
            &input,
            &rows,
            0,
            &PreserverImageSearchOptions {
                library: Some(library()),
                ..PreserverImageSearchOptions::default()
            },
        )
        .expect("the search runs");
        report.verify().expect("the report replays");
        let lift = report.proved_lifts().next().expect("a proved lift");
        assert_eq!(lift.candidate.parent_name, "chebyshev_t");
        assert_eq!(
            lift.lift_summary
                .root_transfer
                .input_root_region
                .to_string(),
            "(0, 1)"
        );
    }
}

// ---------------------------------------------------------------------------
// The Euler step kappa (theta + a)(theta + b) + x (u + v theta)
// ---------------------------------------------------------------------------

mod euler_step {
    use super::*;

    /// `Q_n = T((1 + x)^n)` for `T = (theta + 1)^2 + x (1 + v theta)`, with
    /// the rule computed from the binomial rule.
    fn binomial_image(v: i64) -> (Recurrence, Vec<Vec<Q>>) {
        let rows: Vec<Vec<Q>> = binomial_rows()
            .iter()
            .map(|row| {
                let diagonal = vec_poly::add(
                    &vec_poly::add(row, &vec_poly::scale(&theta(row, |k| q(k as i64)), &q(2))),
                    &theta(row, |k| q((k * k) as i64)),
                );
                let raising = vec_poly::add(row, &theta(row, |k| q(v * k as i64)));
                let mut shifted = vec![q(0)];
                shifted.extend(raising);
                vec_poly::trim(&vec_poly::add(&diagonal, &shifted))
            })
            .collect();
        (
            computed_rule(&candidate("binomial", step(v)), 0, &rows),
            rows,
        )
    }

    fn step(v: i64) -> PreserverOperator {
        let constant = |value: i64| IndexRationalFunction::constant(&q(value));
        PreserverOperator::EulerStep {
            diagonal: vec![constant(1), constant(2), constant(1)],
            raising: vec![constant(1), constant(v)],
        }
    }

    /// `kappa = a = b = 1`, `u = v = 1`: `B = 1 - 3/2 < 0` and
    /// `B^2 = 1/4 < v^2 ((a + b + 1)^2/4 - a b) = 5/4`.
    #[test]
    fn an_euler_step_image_of_the_binomial_rows_is_proved() {
        let (rule, rows) = binomial_image(1);
        let certification = certify(&rule, &rows, 0, &candidate("binomial", step(1)));
        assert_proved(&certification, "(-inf, 0]");
        assert_eq!(
            certification.lean_declarations.first().map(String::as_str),
            Some("RealRooted.IsPFPolynomial.generalStep_of_natDegree_le")
        );
        // Not a preserver of every real-rooted polynomial: no interlacing.
        assert!(!certification.interlacing.is_claimed());
        for row in &rows {
            assert!(is_real_rooted_bigint_coeffs(&integer_row(row)));
        }
    }

    /// `v = 4`: `B = 1 - 6 = -5` and `B^2 = 25 >= 16 * 5/4 = 20`: the
    /// comparison defect is positive somewhere on `s <= 0`, so the theorem
    /// gives nothing (exact witness `n = 0`).
    #[test]
    fn a_positive_comparison_defect_is_an_exact_failure() {
        let (rule, rows) = binomial_image(4);
        let certification = certify(&rule, &rows, 0, &candidate("binomial", step(4)));
        assert!(!certification.is_proved());
        let defect = certification
            .lift_summary
            .hypotheses
            .iter()
            .find(|h| h.name == "euler_step_comparison_defect_negative")
            .expect("the defect hypothesis");
        assert_eq!(
            defect.status,
            LiftHypothesisStatus::Failed {
                witness: Some(BigInt::from(0))
            }
        );
    }
}
