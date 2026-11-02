//! The Ore algebra with the variable shift `E`: the commutation rules on
//! random rows, composition, pseudo-division, the falling-factorial bridge,
//! the pullback contract of every difference row-map step, basis changes on
//! Touchard and Eulerian rows and certified image annihilators.

use ore_algebra::{
    conjugate_to_falling_factorial_basis, conjugate_to_power_basis,
    difference_right_pseudo_division, falling_factorial, falling_factorial_euler_operator,
    falling_factorial_raising_operator, falling_factorial_to_power_basis,
    find_difference_image_annihilator, from_falling_factorial_euler_form,
    power_to_falling_factorial_basis, pull_back_difference, push_forward_difference,
    to_falling_factorial_euler_form, variable_shift_pseudo_division, BivariatePolynomial,
    DifferenceImageAnnihilatorBudget, DifferenceImageAnnihilatorCertificate,
    DifferenceImageAnnihilatorMethod, DifferenceImageAnnihilatorSearch, DifferencePullbackBounds,
    DifferenceRowMapStep, DifferentialDifferenceOperator as Op, OreOperator, Rational, RowMapStep,
};

fn r(value: i64) -> Rational {
    Rational::from_integer(value.into())
}

fn q(numerator: i64, denominator: i64) -> Rational {
    Rational::new(numerator.into(), denominator.into())
}

/// `rows[i][j]` is the coefficient of `n^i x^j`.
fn biv(rows: &[&[i64]]) -> BivariatePolynomial {
    BivariatePolynomial::new(
        rows.iter()
            .map(|row| row.iter().map(|v| r(*v)).collect())
            .collect(),
    )
}

struct Random(u64);

impl Random {
    fn next(&mut self, modulus: u64) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) % modulus
    }

    fn small(&mut self) -> i64 {
        self.next(7) as i64 - 3
    }

    fn polynomial(&mut self, degree: usize) -> Vec<Rational> {
        (0..=degree).map(|_| r(self.small())).collect()
    }

    fn coefficient(&mut self, index_degree: usize, variable_degree: usize) -> BivariatePolynomial {
        BivariatePolynomial::new(
            (0..=index_degree)
                .map(|_| self.polynomial(variable_degree))
                .collect(),
        )
    }

    /// Rows `F_0, ..., F_(count-1)`, `deg F_n <= n % 4 + 1`.
    fn rows(&mut self, count: usize) -> Vec<Vec<Rational>> {
        (0..count).map(|n| self.polynomial(n % 4 + 1)).collect()
    }

    /// A random operator with shifts `<= shift`, derivatives `<= derivative`
    /// and powers of `E` in `low..=high`.
    fn operator(&mut self, shift: usize, derivative: usize, low: i64, high: i64) -> Op {
        let mut op = Op::default();
        for _ in 0..4 {
            let i = self.next(shift as u64 + 1) as usize;
            let d = self.next(derivative as u64 + 1) as usize;
            let e = low + self.next((high - low) as u64 + 1) as i64;
            op.add_term((i, d, e), self.coefficient(1, 2));
        }
        op
    }
}

fn trim(p: &[Rational]) -> Vec<Rational> {
    let mut p = p.to_vec();
    while p.len() > 1 && p.last().is_some_and(|c| *c == r(0)) {
        p.pop();
    }
    if p.is_empty() {
        p.push(r(0));
    }
    p
}

fn evaluate(p: &[Rational], x: &Rational) -> Rational {
    p.iter().rev().fold(r(0), |acc, c| acc * x + c)
}

/// `p(a x + c)`.
fn substitute(p: &[Rational], a: &Rational, c: &Rational) -> Vec<Rational> {
    let mut result = vec![r(0)];
    for coefficient in p.iter().rev() {
        let mut next = vec![r(0); result.len() + 1];
        for (j, value) in result.iter().enumerate() {
            next[j] += value * c;
            next[j + 1] += value * a;
        }
        next[0] += coefficient;
        result = next;
    }
    trim(&result)
}

fn apply_all(op: &Op, rows: &[Vec<Rational>], count: usize) -> Vec<Vec<Rational>> {
    (0..count).map(|n| op.apply(rows, 0, n).unwrap()).collect()
}

fn x_operator() -> Op {
    Op::multiplication(biv(&[&[0, 1]]))
}

#[test]
fn variable_shift_rules_hold_on_random_rows() {
    let mut random = Random(11);
    let e = Op::variable_shift(1);
    let inverse = Op::variable_shift(-1);
    // E x = (x + 1) E and E^(-1) x = (x - 1) E^(-1), structurally ...
    assert_eq!(
        e.compose(&x_operator()),
        Op::single(0, 0, 1, biv(&[&[1, 1]]))
    );
    assert_eq!(
        inverse.compose(&x_operator()),
        Op::single(0, 0, -1, biv(&[&[-1, 1]]))
    );
    // ... E D = D E, S n = (n + 1) S, E n = n E, E E^(-1) = 1 ...
    let d = Op::derivative(1);
    assert_eq!(e.compose(&d), d.compose(&e));
    let n = Op::multiplication(biv(&[&[0], &[1]]));
    assert_eq!(
        Op::index_shift(1).compose(&n),
        Op::single(1, 0, 0, biv(&[&[1], &[1]]))
    );
    assert_eq!(e.compose(&n), n.compose(&e));
    assert_eq!(e.compose(&inverse), Op::identity());
    // ... and on rows: (E x f)(t) = (t + 1) f(t + 1).
    for _ in 0..20 {
        let row = random.polynomial(5);
        let left = e.compose(&x_operator()).apply_row(&row, 0);
        let shifted = substitute(&row, &r(1), &r(1));
        let mut expected = vec![r(0); shifted.len() + 1];
        for (j, c) in shifted.iter().enumerate() {
            expected[j] += c;
            expected[j + 1] += c;
        }
        assert_eq!(left, trim(&expected));
    }
    // The forward and backward differences.
    let row = random.polynomial(6);
    let delta = Op::forward_difference(2).apply_row(&row, 0);
    let x = q(3, 2);
    let expected =
        evaluate(&row, &(&x + r(2))) - r(2) * evaluate(&row, &(&x + r(1))) + evaluate(&row, &x);
    assert_eq!(evaluate(&delta, &x), expected);
    let nabla = Op::backward_difference(1).apply_row(&row, 0);
    assert_eq!(
        evaluate(&nabla, &x),
        evaluate(&row, &x) - evaluate(&row, &(&x - r(1)))
    );
}

#[test]
fn composition_is_the_composition_of_actions_and_associative() {
    let mut random = Random(5);
    for _ in 0..15 {
        let a = random.operator(2, 2, -1, 2);
        let b = random.operator(2, 2, -2, 1);
        let c = random.operator(1, 1, 0, 1);
        assert_eq!(a.compose(&b).compose(&c), a.compose(&b.compose(&c)));
        let rows = random.rows(12);
        let inner = apply_all(&b, &rows, 10);
        let outer = apply_all(&a, &inner, 8);
        assert_eq!(outer, apply_all(&a.compose(&b), &rows, 8));
    }
    // The embedding of operators without E is a homomorphism.
    let mut m = OreOperator::single(1, 1, biv(&[&[1, 2], &[0, 1]]));
    m.add_term((0, 0), biv(&[&[0, 0, 1]]));
    let mut l = OreOperator::single(0, 2, biv(&[&[3, 0, 1]]));
    l.add_term((1, 0), biv(&[&[1], &[1]]));
    assert_eq!(
        Op::from_ore_operator(&m).compose(&Op::from_ore_operator(&l)),
        Op::from_ore_operator(&m.compose(&l))
    );
    assert_eq!(Op::from_ore_operator(&m).to_ore_operator(), Some(m.clone()));
    assert_eq!(Op::variable_shift(1).to_ore_operator(), None);
}

#[test]
fn pseudo_divisions_replay() {
    let mut random = Random(3);
    // L = (n + 1) S^2 - x E D S + (x^2 - n) E^(-1).
    let mut l = Op::single(2, 0, 0, biv(&[&[1], &[1]]));
    l.add_term((1, 1, 1), biv(&[&[0, -1]]));
    l.add_term((0, 0, -1), biv(&[&[0, 0, 1], &[-1]]));
    for _ in 0..10 {
        let p = random.operator(4, 2, -1, 2);
        let division = difference_right_pseudo_division(&p, &l, 4).unwrap();
        assert!(division.remainder.is_zero() || division.remainder.max_shift() < 2);
        assert!(division.multiplier.is_variable_free());
        assert_eq!(
            p.left_multiply(&division.multiplier),
            division.quotient.compose(&l).add(&division.remainder)
        );
    }
    // A = (x + 1) E^2 - n E + x, division in E.
    let mut a = Op::single(0, 0, 2, biv(&[&[1, 1]]));
    a.add_term((0, 0, 1), biv(&[&[0], &[-1]]));
    a.add_term((0, 0, 0), biv(&[&[0, 1]]));
    for _ in 0..10 {
        let p = random.operator(0, 2, 0, 4);
        let division = variable_shift_pseudo_division(&p, &a, 4, 2).unwrap();
        assert!(division.remainder.max_variable_shift() < 2);
        assert_eq!(
            p.left_multiply(&division.multiplier),
            division.quotient.compose(&a).add(&division.remainder)
        );
    }
    // A multiple of A leaves no remainder.
    let multiple = random.operator(0, 1, 0, 2).compose(&a);
    let division = variable_shift_pseudo_division(&multiple, &a, 4, 1).unwrap();
    assert!(division.remainder.is_zero());
}

#[test]
fn falling_factorial_bridge() {
    let mut random = Random(17);
    let theta = falling_factorial_euler_operator();
    let raising = falling_factorial_raising_operator();
    for k in 0..8 {
        let ff = falling_factorial(k);
        let scaled: Vec<Rational> = ff.iter().map(|c| c * r(k as i64)).collect();
        assert_eq!(theta.apply_row(&ff, 0), trim(&scaled));
        assert_eq!(raising.apply_row(&ff, 0), falling_factorial(k + 1));
        assert_eq!(falling_factorial_to_power_basis(&ff), {
            let mut monomial = vec![r(0); k + 1];
            monomial[k] = r(1);
            monomial
        });
    }
    // Touchard: phi(x^3) = x + 3 x^2 + x^3.
    assert_eq!(
        falling_factorial_to_power_basis(&[r(0), r(0), r(0), r(1)]),
        vec![r(0), r(1), r(3), r(1)]
    );
    for _ in 0..10 {
        let row = random.polynomial(6);
        assert_eq!(
            power_to_falling_factorial_basis(&falling_factorial_to_power_basis(&row)),
            trim(&row)
        );
        // phi o P = psi(P) o phi for P free of D with E >= 0.
        let p = random.operator(0, 0, 0, 3);
        let psi = conjugate_to_power_basis(&p).unwrap();
        assert_eq!(
            falling_factorial_to_power_basis(&p.apply_row(&row, 2)),
            psi.apply_row(&falling_factorial_to_power_basis(&row), 2)
        );
        // phi o psi^(-1)(Q) = Q o phi for Q free of E.
        let mut ore = OreOperator::single(0, 2, random.coefficient(1, 2));
        ore.add_term((0, 0), random.coefficient(1, 2));
        ore.add_term((0, 1), random.coefficient(1, 1));
        let back = conjugate_to_falling_factorial_basis(&ore);
        assert!(back.is_derivative_free());
        assert_eq!(
            falling_factorial_to_power_basis(&back.apply_row(&row, 3)),
            ore.apply_row(&falling_factorial_to_power_basis(&row), 3)
        );
        // psi(psi^(-1)(Q)) = Q after making the powers of E nonnegative.
        let shift = (-back.min_variable_shift()).max(0);
        let round_trip = conjugate_to_power_basis(&back.variable_shifted(shift)).unwrap();
        let mut one_plus_d = OreOperator::identity();
        one_plus_d.add_term((0, 1), BivariatePolynomial::one());
        let power = (0..shift).fold(OreOperator::identity(), |acc, _| acc.compose(&one_plus_d));
        assert_eq!(round_trip, power.compose(&ore));
        // The falling-factorial Euler form: X^J P = sum X^s a(n, Theta) S^i.
        let (j, form) = to_falling_factorial_euler_form(&p).unwrap();
        let x_power = raising.power(j);
        assert_eq!(
            from_falling_factorial_euler_form(&form),
            x_power.compose(&p)
        );
    }
    assert!(conjugate_to_power_basis(&Op::derivative(1)).is_err());
    assert!(conjugate_to_power_basis(&Op::variable_shift(-1)).is_err());
}

/// `mu(n) = prod_(t=1..n) (t + 2)`, the gauge with ratio `(n + 2)/1`.
fn gauge_value(n: usize) -> Rational {
    (1..=n).fold(r(1), |acc, t| acc * r(t as i64 + 2))
}

#[test]
fn pullback_contracts_on_random_rows() {
    let mut random = Random(23);
    let points = [q(1, 3), r(2), q(-5, 2)];
    for round in 0..6 {
        let rows = random.rows(12);
        let general = random.operator(2, 2, -1, 2);
        let derivative_free = random.operator(2, 0, -1, 2);
        let shift_free = {
            let mut op = random.operator(2, 2, 0, 0);
            op.add_term((2, 0, 0), BivariatePolynomial::one());
            op
        };
        let mut t = Op::single(0, 1, 1, biv(&[&[1, 1]]));
        t.add_term((0, 0, -1), biv(&[&[2], &[1]]));
        let gauge = RowMapStep::Gauge {
            numerator: vec![r(2), r(1)],
            denominator: vec![r(1)],
        };
        let slope = r([1, -1, 2][round % 3]);
        let affine = RowMapStep::Affine {
            slope: slope.clone(),
            intercept: q(1, 2),
        };
        let cases: Vec<(Op, DifferenceRowMapStep)> = vec![
            (general.clone(), DifferenceRowMapStep::Operator(t.clone())),
            (general.clone(), DifferenceRowMapStep::forward_difference(2)),
            (
                shift_free.clone(),
                DifferenceRowMapStep::FallingFactorialToPowerBasis,
            ),
            (
                derivative_free.clone(),
                DifferenceRowMapStep::PowerToFallingFactorialBasis,
            ),
            (general.clone(), DifferenceRowMapStep::Ore(gauge.clone())),
            (general.clone(), DifferenceRowMapStep::Ore(affine.clone())),
            (
                general.clone(),
                DifferenceRowMapStep::Ore(RowMapStep::WeylOperator(
                    shift_free.column(0).to_ore_operator().unwrap(),
                )),
            ),
        ];
        for (m, step) in cases {
            let bounds = DifferencePullbackBounds::of(&m);
            let (p, base) = pull_back_difference(&m, std::slice::from_ref(&step), bounds)
                .unwrap_or_else(|e| panic!("{step:?}: {e}"));
            let image: Vec<Vec<Rational>> = rows
                .iter()
                .enumerate()
                .map(|(n, row)| match &step {
                    DifferenceRowMapStep::Ore(RowMapStep::Gauge { .. }) => {
                        row.iter().map(|c| c * gauge_value(n)).collect()
                    }
                    _ => step.apply_row(row, n).unwrap(),
                })
                .collect();
            for n in 0..9 {
                let left = m.apply(&image, 0, n).unwrap();
                let right = p.apply(&rows, 0, n).unwrap();
                match &step {
                    DifferenceRowMapStep::Operator(_) => assert_eq!(left, right),
                    DifferenceRowMapStep::FallingFactorialToPowerBasis => {
                        assert_eq!(left, falling_factorial_to_power_basis(&right))
                    }
                    DifferenceRowMapStep::PowerToFallingFactorialBasis => {
                        let c = (-bounds.lowest_variable_shift).max(0);
                        let back = power_to_falling_factorial_basis(&right);
                        assert_eq!(left, substitute(&back, &r(1), &r(-c)));
                    }
                    DifferenceRowMapStep::Ore(RowMapStep::Gauge { .. }) => {
                        let b = evaluate(&base, &r(n as i64));
                        for x in &points {
                            assert_eq!(
                                evaluate(&left, x) * &b,
                                evaluate(&right, x) * gauge_value(n)
                            );
                        }
                    }
                    DifferenceRowMapStep::Ore(RowMapStep::Affine { .. }) => {
                        assert_eq!(left, substitute(&right, &slope, &q(1, 2)))
                    }
                    DifferenceRowMapStep::Ore(_) => assert_eq!(left, right),
                    DifferenceRowMapStep::Substitution(_) => unreachable!(),
                }
            }
        }
    }
    // Steps that need operators free of D or of E report it.
    let bounds = DifferencePullbackBounds::of(&Op::derivative(1));
    assert!(pull_back_difference(
        &Op::derivative(1),
        &[DifferenceRowMapStep::PowerToFallingFactorialBasis],
        bounds
    )
    .is_err());
    let e = Op::variable_shift(1);
    assert!(pull_back_difference(
        &e,
        &[DifferenceRowMapStep::FallingFactorialToPowerBasis],
        DifferencePullbackBounds::of(&e)
    )
    .is_err());
    assert!(pull_back_difference(
        &e,
        &[DifferenceRowMapStep::Ore(RowMapStep::Reversal {
            degree: vec![r(3)]
        })],
        DifferencePullbackBounds::of(&e)
    )
    .is_err());
}

/// Rows `F_0, ..., F_(count-1)` of `F_(n+1) = A F_n` for a shift-free `A`.
fn first_order_rows(step: &Op, first: Vec<Rational>, count: usize) -> Vec<Vec<Rational>> {
    let mut rows = vec![first];
    while rows.len() < count {
        let n = rows.len() - 1;
        let next = step.apply_row(rows.last().unwrap(), n);
        rows.push(next);
    }
    rows
}

fn first_order_parent(step: &Op) -> Op {
    Op::index_shift(1).subtract(step)
}

fn assert_annihilates(op: &Op, rows: &[Vec<Rational>]) {
    for n in 0..rows.len() - op.max_shift() {
        let value = op.apply(rows, 0, n).unwrap();
        assert!(value.iter().all(|c| *c == r(0)), "base {n}: {value:?}");
    }
}

fn image_rows(steps: &[DifferenceRowMapStep], rows: &[Vec<Rational>]) -> Vec<Vec<Rational>> {
    rows.iter()
        .enumerate()
        .map(|(n, row)| {
            steps.iter().fold(row.clone(), |current, step| {
                step.apply_row(&current, n).unwrap()
            })
        })
        .collect()
}

/// Touchard: `T_(n+1) = x (T_n + T_n')`.
fn touchard_step() -> Op {
    let mut step = Op::multiplication(biv(&[&[0, 1]]));
    step.add_term((0, 1, 0), biv(&[&[0, 1]]));
    step
}

/// Eulerian: `A_(n+1) = (1 + n x) A_n + x (1 - x) A_n'`.
fn eulerian_step() -> Op {
    let mut step = Op::multiplication(biv(&[&[1], &[0, 1]]));
    step.add_term((0, 1, 0), biv(&[&[0, 1, -1]]));
    step
}

fn found(
    search: DifferenceImageAnnihilatorSearch,
) -> (
    DifferenceImageAnnihilatorCertificate,
    DifferenceImageAnnihilatorMethod,
) {
    match search {
        DifferenceImageAnnihilatorSearch::Found {
            certificate,
            method,
        } => (*certificate, method),
        other => panic!("not found: {other:?}"),
    }
}

#[test]
fn basis_changes_on_touchard_and_eulerian_rows() {
    let budget = DifferenceImageAnnihilatorBudget::default();
    // x^n -> Touchard polynomials under (x)_k -> x^k.
    let powers_parent = first_order_parent(&x_operator());
    let powers = first_order_rows(&x_operator(), vec![r(1)], 40);
    let to_power = [DifferenceRowMapStep::FallingFactorialToPowerBasis];
    let touchard = image_rows(&to_power, &powers);
    assert_eq!(touchard, first_order_rows(&touchard_step(), vec![r(1)], 40));
    let (certificate, method) =
        found(find_difference_image_annihilator(&powers_parent, &to_power, &budget).unwrap());
    assert_eq!(method, DifferenceImageAnnihilatorMethod::PushForward);
    certificate.verify().unwrap();
    assert_eq!(
        certificate.annihilator,
        first_order_parent(&touchard_step())
    );
    // Eulerian rows in the falling-factorial basis: sum_k A(n, k) (x)_k.
    let eulerian_parent = first_order_parent(&eulerian_step());
    let eulerian = first_order_rows(&eulerian_step(), vec![r(1)], 40);
    let to_falling = [DifferenceRowMapStep::PowerToFallingFactorialBasis];
    let (certificate, method) =
        found(find_difference_image_annihilator(&eulerian_parent, &to_falling, &budget).unwrap());
    assert_eq!(method, DifferenceImageAnnihilatorMethod::PushForward);
    certificate.verify().unwrap();
    let image = image_rows(&to_falling, &eulerian);
    assert_annihilates(&certificate.annihilator, &image);
    assert!(certificate.annihilator.is_derivative_free());
    // There and back: the composite is the identity on rows, and the
    // pulled-back operator annihilates the source rows.
    let there_and_back = [
        DifferenceRowMapStep::PowerToFallingFactorialBasis,
        DifferenceRowMapStep::FallingFactorialToPowerBasis,
    ];
    assert_eq!(image_rows(&there_and_back, &eulerian), eulerian);
    let (certificate, _) = found(
        find_difference_image_annihilator(&eulerian_parent, &there_and_back, &budget).unwrap(),
    );
    certificate.verify().unwrap();
    assert_annihilates(&certificate.annihilator, &eulerian);
    // The falling factorials (x)_n: (x)_(n+1) = x (x - 1)_n, a parent with E.
    let falling_parent = first_order_parent(&falling_factorial_raising_operator());
    let falling = first_order_rows(&falling_factorial_raising_operator(), vec![r(1)], 40);
    assert_eq!(falling[7], falling_factorial(7));
    assert_eq!(image_rows(&to_power, &falling), powers);
    let pushed = push_forward_difference(&falling_parent, &to_power[0]).unwrap();
    assert_annihilates(&pushed, &powers);
}

#[test]
fn image_annihilators_of_shifted_touchard_rows() {
    let budget = DifferenceImageAnnihilatorBudget::default();
    let parent = first_order_parent(&touchard_step());
    let touchard = first_order_rows(&touchard_step(), vec![r(1)], 40);
    // f(x) -> f(x + 1) as the operator E: pushed forward by conjugation.
    let shift = [DifferenceRowMapStep::Operator(Op::variable_shift(1))];
    let (certificate, method) =
        found(find_difference_image_annihilator(&parent, &shift, &budget).unwrap());
    assert_eq!(method, DifferenceImageAnnihilatorMethod::PushForward);
    certificate.verify().unwrap();
    let shifted = image_rows(&shift, &touchard);
    for (n, row) in shifted.iter().enumerate() {
        assert_eq!(*row, substitute(&touchard[n], &r(1), &r(1)));
    }
    assert_annihilates(&certificate.annihilator, &shifted);
    // The same map as a translation step agrees.
    let translation = [DifferenceRowMapStep::translation(r(1))];
    assert_eq!(image_rows(&translation, &touchard), shifted);
    let (translated, _) =
        found(find_difference_image_annihilator(&parent, &translation, &budget).unwrap());
    translated.verify().unwrap();
    assert_annihilates(&translated.annihilator, &shifted);
    // The forward difference T_n(x + 1) - T_n(x) needs a common left
    // multiple.
    let difference = [DifferenceRowMapStep::forward_difference(1)];
    let (certificate, method) =
        found(find_difference_image_annihilator(&parent, &difference, &budget).unwrap());
    assert!(matches!(
        method,
        DifferenceImageAnnihilatorMethod::LeftMultiple(_)
    ));
    let base_factor = certificate.verify().unwrap();
    assert_eq!(base_factor, vec![r(1)]);
    let image = image_rows(&difference, &touchard);
    assert_annihilates(&certificate.annihilator, &image);
    // Tampering with the quotient, the multiplier, the annihilator or the
    // steps breaks the replay.
    let mut tampered = certificate.clone();
    tampered
        .quotient
        .add_term((0, 0, 1), BivariatePolynomial::one());
    assert!(tampered.verify().is_err());
    let mut tampered = certificate.clone();
    tampered.multiplier = tampered.multiplier.scale(&r(2));
    assert!(tampered.verify().is_err());
    let mut tampered = certificate.clone();
    tampered
        .annihilator
        .add_term((0, 0, 0), BivariatePolynomial::one());
    assert!(tampered.verify().is_err());
    let mut tampered = certificate.clone();
    tampered.steps = vec![DifferenceRowMapStep::forward_difference(2)];
    assert!(tampered.verify().is_err());
    let mut tampered = certificate;
    tampered.multiplier = biv(&[&[0, 1]]);
    assert!(tampered.verify().is_err());
}

#[test]
fn budget_exhaustion_is_reported() {
    let parent = first_order_parent(&eulerian_step());
    let budget = DifferenceImageAnnihilatorBudget {
        max_order_increase: 0,
        max_derivative_order: 0,
        max_variable_shift_order: 0,
        max_index_degree: 1,
        max_variable_degree: 1,
        max_unknowns: 400,
        max_shapes: 60,
    };
    let search = find_difference_image_annihilator(
        &parent,
        &[DifferenceRowMapStep::forward_difference(1)],
        &budget,
    )
    .unwrap();
    assert!(matches!(
        search,
        DifferenceImageAnnihilatorSearch::NotFoundWithinBudget { shapes_tried, .. } if shapes_tried > 0
    ));
}
