//! Exact image recurrences of operators: computed rules against directly
//! computed image rows, budget exhaustion, replay and tampering.

use num_bigint::BigInt;
use polytool::recurrence::{
    image_recurrence_of_operator, BivarPolyJson, DegreeSchedule, ImageGauge,
    ImageOperatorSpecification, ImageRecurrence, ImageRecurrenceMethod, ImageRecurrenceOptions,
    ImageRowMap, OperatorTermJson, PreserverParentLibrary, RecurrenceError, RecurrenceJson,
};
use polytool::{vec_poly, BigRational};

type Q = BigRational;

const ROWS: usize = 40;

use crate::common::q;

fn parent(name: &str) -> RecurrenceJson {
    PreserverParentLibrary::named_family_records()
        .into_iter()
        .find(|(n, _)| n == name)
        .unwrap_or_else(|| panic!("{name} is a named record"))
        .1
}

/// The parent rows from its first index.
fn parent_rows(record: &RecurrenceJson, count: usize) -> (usize, Vec<Vec<Q>>) {
    let (rule, first, seeds) = record.to_recurrence_parts().unwrap();
    let rows = rule
        .generate_rows_rational(&seeds, first, count.max(seeds.len()))
        .unwrap()
        .into_iter()
        .take(count)
        .map(|row| vec_poly::trim(&row))
        .collect();
    (first, rows)
}

/// `coeffs[a][b]` multiplies `n^a x^b`.
fn term(deriv_order: usize, coeffs: &[&[i64]]) -> OperatorTermJson {
    OperatorTermJson {
        deriv_order,
        coeff: BivarPolyJson {
            coeffs: coeffs
                .iter()
                .map(|row| row.iter().map(|v| v.to_string()).collect())
                .collect(),
        },
    }
}

fn derivative(row: &[Q], order: usize) -> Vec<Q> {
    (0..order).fold(row.to_vec(), |p, _| vec_poly::derivative(&p))
}

/// The image rule's rows against `expected` (first index `n0`).
fn assert_rows(image: &ImageRecurrence, expected: &[Vec<Q>]) {
    image.verify().unwrap();
    let rows = image.rows(expected.len()).unwrap();
    for (offset, (row, want)) in rows.iter().zip(expected).enumerate() {
        assert_eq!(&vec_poly::trim(row), want, "row offset {offset}");
    }
    assert_eq!(rows.len(), expected.len());
}

fn compute(
    record: &RecurrenceJson,
    n0: usize,
    specification: &ImageOperatorSpecification,
) -> ImageRecurrence {
    image_recurrence_of_operator(
        record,
        n0,
        specification,
        &ImageRecurrenceOptions::default(),
    )
    .unwrap()
}

fn weyl(terms: Vec<OperatorTermJson>) -> ImageOperatorSpecification {
    ImageOperatorSpecification::new(ImageRowMap::WeylOperator {
        terms,
        denominator: vec!["1".into()],
    })
}

#[test]
fn eulerian_plus_derivative() {
    let record = parent("eulerian");
    let (n0, rows) = parent_rows(&record, ROWS);
    let image = compute(&record, n0, &weyl(vec![term(0, &[&[1]]), term(1, &[&[1]])]));
    // No first-order rule exists; the common left multiple has order 2.
    assert_eq!(image.order, 2);
    assert!(matches!(
        image.method,
        ImageRecurrenceMethod::LeftMultiple { .. }
    ));
    let expected: Vec<_> = rows
        .iter()
        .map(|p| vec_poly::trim(&vec_poly::add(p, &derivative(p, 1))))
        .collect();
    assert_rows(&image, &expected);
}

#[test]
fn euler_operator_polynomial_on_narayana() {
    // theta (theta + 1) = x^2 D^2 + 2 x D.
    let record = parent("narayana");
    let (n0, rows) = parent_rows(&record, ROWS);
    let image = compute(
        &record,
        n0,
        &weyl(vec![term(1, &[&[0, 2]]), term(2, &[&[0, 0, 1]])]),
    );
    // With the contiguity relation and the differential equation of the
    // Narayana rows (a jointly invariant system), the rule has order one.
    assert_eq!(image.order, 1);
    let system = image.certificate.annihilator_system.as_ref().unwrap();
    assert_eq!(system.members.len(), 2);
    // Tampering with the system is detected.
    let mut tampered = image.clone();
    let member = &mut tampered
        .certificate
        .annihilator_system
        .as_mut()
        .unwrap()
        .members[0];
    member.checked_bases.pop();
    assert!(tampered.verify().is_err());
    let mut tampered = image.clone();
    let member = &mut tampered
        .certificate
        .annihilator_system
        .as_mut()
        .unwrap()
        .members[1];
    member.multiplier.coeffs[0][0] = "99".into();
    assert!(tampered.verify().is_err());
    let expected: Vec<_> = rows
        .iter()
        .map(|p| {
            vec_poly::trim(
                &p.iter()
                    .enumerate()
                    .map(|(k, c)| c * q((k * (k + 1)) as i64))
                    .collect::<Vec<_>>(),
            )
        })
        .collect();
    assert_rows(&image, &expected);
}

#[test]
fn derivative_of_x_times_eulerian() {
    // D o x = x D + 1.
    let record = parent("eulerian");
    let (n0, rows) = parent_rows(&record, ROWS);
    let image = compute(
        &record,
        n0,
        &weyl(vec![term(0, &[&[1]]), term(1, &[&[0, 1]])]),
    );
    let expected: Vec<_> = rows
        .iter()
        .map(|p| vec_poly::derivative(&vec_poly::shift(p)))
        .map(|p| vec_poly::trim(&p))
        .collect();
    assert_rows(&image, &expected);
}

#[test]
fn affine_substitution_pushes_forward() {
    let record = parent("eulerian");
    let (n0, rows) = parent_rows(&record, ROWS);
    let image = compute(
        &record,
        n0,
        &ImageOperatorSpecification::new(ImageRowMap::Affine {
            slope: "2".into(),
            intercept: "1".into(),
        }),
    );
    assert_eq!(image.method, ImageRecurrenceMethod::PushForward);
    assert_eq!(image.order, 1);
    let expected: Vec<_> = rows
        .iter()
        .map(|p| vec_poly::trim(&vec_poly::compose(p, &[q(1), q(2)])))
        .collect();
    assert_rows(&image, &expected);
}

#[test]
fn reversal_pushes_forward() {
    // Narayana N_n has degree n - 1 from n = 1 on (reverse with d_n = n - 1).
    let record = parent("narayana");
    let (n0, rows) = parent_rows(&record, ROWS);
    let degree0 = vec_poly::degree(&rows[0]).unwrap();
    let image = compute(
        &record,
        n0,
        &ImageOperatorSpecification::new(ImageRowMap::Reversal {
            schedule: DegreeSchedule::new(n0, degree0, 1),
        }),
    );
    assert_eq!(image.method, ImageRecurrenceMethod::PushForward);
    let expected: Vec<_> = rows
        .iter()
        .enumerate()
        .map(|(offset, p)| vec_poly::trim(&vec_poly::reverse(p, degree0 + offset).unwrap()))
        .collect();
    assert_rows(&image, &expected);
}

#[test]
fn inverse_factorial_multiplier_on_binomial_rows() {
    let record = parent("binomial");
    let (n0, rows) = parent_rows(&record, ROWS);
    let image = compute(
        &record,
        n0,
        &ImageOperatorSpecification::new(ImageRowMap::CoefficientMultiplier {
            numerator: vec!["1".into()],
            denominator: vec!["1".into(), "1".into()],
        }),
    );
    // The recurrence alone gives no runnable rule; the row annihilator
    // (1 + x) D - n of the binomial rows does (the Laguerre recurrence).
    let system = image.certificate.annihilator_system.as_ref().unwrap();
    assert_eq!(system.members.len(), 1);
    assert_eq!(image.order, 1);
    let expected: Vec<_> = rows
        .iter()
        .map(|p| {
            let mut factorial = q(1);
            vec_poly::trim(
                &p.iter()
                    .enumerate()
                    .map(|(k, c)| {
                        if k > 0 {
                            factorial *= q(k as i64);
                        }
                        c / &factorial
                    })
                    .collect::<Vec<_>>(),
            )
        })
        .collect();
    assert_rows(&image, &expected);
}

#[test]
fn power_of_two_gauge_with_an_index_shift() {
    // Q_n = 2^n N_(n+1), n >= n0.
    let record = parent("narayana");
    let (p0, rows) = parent_rows(&record, ROWS + 1);
    let n0 = p0;
    let specification = ImageOperatorSpecification {
        index_shift: 1,
        gauge: ImageGauge::Ratio {
            numerator: vec!["2".into()],
            denominator: vec!["1".into()],
            initial_value: (BigInt::from(2).pow(n0 as u32)).to_string(),
        },
        map: ImageRowMap::Identity,
    };
    let image = compute(&record, n0, &specification);
    assert_eq!(image.method, ImageRecurrenceMethod::PushForward);
    let expected: Vec<_> = (0..ROWS)
        .map(|offset| {
            let power = Q::from_integer(BigInt::from(2).pow((n0 + offset) as u32));
            vec_poly::trim(&vec_poly::scale(&rows[offset + 1], &power))
        })
        .collect();
    assert_rows(&image, &expected);
}

#[test]
fn budget_exhaustion_is_typed() {
    let record = parent("eulerian");
    let (n0, _) = parent_rows(&record, 2);
    let options = ImageRecurrenceOptions {
        max_order_increase: 0,
        max_derivative_order: 1,
        max_coefficient_degree: 1,
        row_annihilators: false,
        ..ImageRecurrenceOptions::default()
    };
    let result = image_recurrence_of_operator(
        &record,
        n0,
        &weyl(vec![term(0, &[&[1]]), term(1, &[&[1]])]),
        &options,
    );
    assert!(
        matches!(result, Err(RecurrenceError::BudgetExhausted(_))),
        "{result:?}"
    );
}

#[test]
fn certificates_replay_and_detect_tampering() {
    let record = parent("eulerian");
    let (n0, _) = parent_rows(&record, 2);
    let image = compute(&record, n0, &weyl(vec![term(0, &[&[1]]), term(1, &[&[1]])]));
    // JSON round trip.
    let json = serde_json::to_string(&image).unwrap();
    let back: ImageRecurrence = serde_json::from_str(&json).unwrap();
    assert_eq!(back, image);
    back.verify().unwrap();
    // A changed seed.
    let mut tampered = image.clone();
    tampered.certificate.rule.initial_polynomials[1][0] = "7".into();
    tampered.rule = tampered.certificate.rule.clone();
    assert!(tampered.verify().is_err());
    // A changed quotient.
    let mut tampered = image.clone();
    tampered.certificate.quotient.terms[0].coeff.coeffs[0][0] = "12345".into();
    assert!(tampered.verify().is_err());
    // A changed exceptional-base record.
    let mut tampered = image.clone();
    tampered.certificate.exceptional_bases.push(1000);
    assert!(tampered.verify().is_err());
    // A changed specification.
    let mut tampered = image.clone();
    tampered.certificate.specification = weyl(vec![term(0, &[&[1]]), term(1, &[&[2]])]);
    assert!(tampered.verify().is_err());
    // A rule that differs from the certified one.
    let mut tampered = image;
    tampered.rule.first_index += 1;
    assert!(tampered.verify().is_err());
}

/// The example of `doc/library.md`.
#[test]
fn library_documentation_example() {
    let (_, eulerian) = PreserverParentLibrary::named_family_records()
        .into_iter()
        .find(|(name, _)| name == "eulerian")
        .unwrap();
    let term = |deriv_order: usize| OperatorTermJson {
        deriv_order,
        coeff: BivarPolyJson {
            coeffs: vec![vec!["1".into()]],
        },
    };
    // T = 1 + D: A_n + A_n'.
    let specification = ImageOperatorSpecification::new(ImageRowMap::WeylOperator {
        terms: vec![term(0), term(1)],
        denominator: vec!["1".into()],
    });
    let image = image_recurrence_of_operator(
        &eulerian,
        eulerian.first_index,
        &specification,
        &ImageRecurrenceOptions::default(),
    )
    .unwrap();
    assert_eq!(image.order, 2);
    image.verify().unwrap();
}
