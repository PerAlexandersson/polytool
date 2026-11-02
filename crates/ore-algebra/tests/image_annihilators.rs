//! Image annihilators: push-forwards, common left multiples and their
//! certificates, checked against directly computed image rows.

use ore_algebra::{
    find_annihilator_system, find_image_annihilator, AnnihilatorSystemBudget, BivariatePolynomial,
    ImageAnnihilatorBudget, ImageAnnihilatorMethod, ImageAnnihilatorSearch, OreOperator, Rational,
    RowMapStep,
};

fn r(value: i64) -> Rational {
    Rational::from_integer(value.into())
}

/// `rows[i][j]` is the coefficient of `n^i x^j`.
fn biv(rows: &[&[i64]]) -> BivariatePolynomial {
    BivariatePolynomial::new(
        rows.iter()
            .map(|row| row.iter().map(|v| r(*v)).collect())
            .collect(),
    )
}

/// Rows `F_0, ..., F_(count-1)` of the first-order rule `F_(n+1) = A F_n`
/// for a shift-free `A`.
fn first_order_rows(step: &OreOperator, first: Vec<Rational>, count: usize) -> Vec<Vec<Rational>> {
    let mut rows = vec![first];
    while rows.len() < count {
        let n = rows.len() - 1;
        let next = step.apply_row(rows.last().unwrap(), n);
        rows.push(next);
    }
    rows
}

/// The parent `L = S - A` of a first-order rule.
fn first_order_parent(step: &OreOperator) -> OreOperator {
    OreOperator::single(1, 0, BivariatePolynomial::one()).subtract(step)
}

fn found(
    search: ImageAnnihilatorSearch,
) -> (
    ore_algebra::ImageAnnihilatorCertificate,
    ImageAnnihilatorMethod,
) {
    match search {
        ImageAnnihilatorSearch::Found {
            certificate,
            method,
        } => (*certificate, method),
        other => panic!("not found: {other:?}"),
    }
}

fn image_rows(steps: &[RowMapStep], rows: &[Vec<Rational>]) -> Vec<Vec<Rational>> {
    rows.iter()
        .enumerate()
        .map(|(n, row)| {
            steps.iter().fold(row.clone(), |current, step| {
                step.apply_row(&current, n).unwrap()
            })
        })
        .collect()
}

fn assert_annihilates(op: &OreOperator, rows: &[Vec<Rational>]) {
    for n in 0..rows.len() - op.max_shift() {
        let value = op.apply(rows, 0, n).unwrap();
        assert!(value.iter().all(|c| *c == r(0)), "base {n}: {value:?}");
    }
}

/// Eulerian polynomials: `A_(n+1) = (1 + n x) A_n + x (1 - x) A_n'`.
fn eulerian_step() -> OreOperator {
    let mut step = OreOperator::single(0, 0, biv(&[&[1], &[0, 1]]));
    step.add_term((0, 1), biv(&[&[0, 1, -1]]));
    step
}

#[test]
fn eulerian_plus_derivative_needs_order_two() {
    let step = eulerian_step();
    let parent = first_order_parent(&step);
    let mut t = OreOperator::identity();
    t.add_term((0, 1), BivariatePolynomial::one());
    let steps = vec![RowMapStep::WeylOperator(t)];
    let search =
        find_image_annihilator(&parent, None, &steps, &ImageAnnihilatorBudget::default()).unwrap();
    let (certificate, method) = found(search);
    assert!(matches!(method, ImageAnnihilatorMethod::LeftMultiple(_)));
    certificate.verify().unwrap();
    let m = &certificate.annihilator;
    // Order 2 (no first-order rule exists: the square linear system of
    // order one has only the trivial solution).
    assert_eq!(m.max_shift(), 2, "{m:?}");
    let rows = first_order_rows(&step, vec![r(1)], 30);
    assert_annihilates(m, &image_rows(&steps, &rows));
    // Tampering with the quotient breaks the identity.
    let mut tampered = certificate.clone();
    tampered
        .quotient
        .add_term((0, 0), BivariatePolynomial::one());
    assert!(tampered.verify().is_err());
}

#[test]
fn affine_and_reversal_push_forward() {
    let step = eulerian_step();
    let parent = first_order_parent(&step);
    let rows = first_order_rows(&step, vec![r(1)], 20);
    for steps in [
        vec![RowMapStep::Affine {
            slope: r(2),
            intercept: r(1),
        }],
        // A_n has degree n - 1 for n >= 1 (A_0 = 1): reverse with d_n = n.
        vec![RowMapStep::Reversal {
            degree: vec![r(0), r(1)],
        }],
    ] {
        let (certificate, method) = found(
            find_image_annihilator(&parent, None, &steps, &ImageAnnihilatorBudget::default())
                .unwrap(),
        );
        assert_eq!(method, ImageAnnihilatorMethod::PushForward);
        certificate.verify().unwrap();
        assert_eq!(certificate.annihilator.max_shift(), 1);
        assert_annihilates(&certificate.annihilator, &image_rows(&steps, &rows));
    }
}

#[test]
fn inverse_factorial_multiplier_on_binomial_rows() {
    // (1 + x)^n: F_(n+1) = (1 + x) F_n.
    let step = OreOperator::single(0, 0, biv(&[&[1, 1]]));
    let parent = first_order_parent(&step);
    let steps = vec![RowMapStep::CoefficientMultiplier {
        numerator: vec![r(1)],
        denominator: vec![r(1), r(1)],
    }];
    // The recurrence alone does not determine a runnable image rule: the
    // binomial rows also satisfy (1 + x) F_n' = n F_n.
    let budget = ImageAnnihilatorBudget::default();
    let search = find_image_annihilator(&parent, None, &steps, &budget).unwrap();
    assert!(matches!(
        search,
        ImageAnnihilatorSearch::NotFoundWithinBudget { .. }
    ));
    let rows = first_order_rows(&step, vec![r(1)], 25);
    let system =
        find_annihilator_system(&parent, &rows[..12], 0, &AnnihilatorSystemBudget::default())
            .unwrap();
    system.verify().unwrap();
    let a = system.row_annihilator().unwrap();
    assert_eq!(a.max_derivative(), 1);
    assert_annihilates(a, &rows);
    let search = find_image_annihilator(&parent, Some(a), &steps, &budget).unwrap();
    let (certificate, method) = found(search);
    certificate.verify().unwrap();
    let m = &certificate.annihilator;
    // The Laguerre recurrence (n + 1) Q_(n+1) = (n + 1 + x) Q_n + x Q_n'.
    assert!(matches!(method, ImageAnnihilatorMethod::LeftMultiple(_)));
    assert_eq!(m.max_shift(), 1, "{m:?}");
    assert_annihilates(m, &image_rows(&steps, &rows));
}

#[test]
fn budget_exhaustion_is_reported() {
    let step = eulerian_step();
    let parent = first_order_parent(&step);
    let mut t = OreOperator::identity();
    t.add_term((0, 1), BivariatePolynomial::one());
    let budget = ImageAnnihilatorBudget {
        max_order_increase: 0,
        max_derivative_order: 1,
        max_index_degree: 1,
        max_variable_degree: 1,
        max_unknowns: 400,
        max_shapes: 60,
    };
    let search =
        find_image_annihilator(&parent, None, &[RowMapStep::WeylOperator(t)], &budget).unwrap();
    assert!(
        matches!(search, ImageAnnihilatorSearch::NotFoundWithinBudget { shapes_tried, .. } if shapes_tried > 0)
    );
}
