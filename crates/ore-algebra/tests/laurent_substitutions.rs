//! Localization at `x` and substitutions: the Laurent conjugations agree
//! with the existing reversal and Möbius pullbacks, power gauges divide rows
//! by powers of `x`, forward rules annihilate image rows and their
//! certificates replay.

use ore_algebra::{
    find_difference_image_annihilator, pull_back, pull_back_substitution, push_forward,
    push_forward_substitution, BivariatePolynomial, DifferenceImageAnnihilatorBudget,
    DifferenceImageAnnihilatorMethod, DifferenceImageAnnihilatorSearch, DifferenceRowMapStep,
    DifferentialDifferenceOperator, LaurentBivariatePolynomial, LaurentOreOperator, MobiusMap,
    OreOperator, PullbackBounds, Rational, RowMapStep, Substitution,
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

    fn polynomial(&mut self, degree: usize) -> Vec<Rational> {
        (0..=degree).map(|_| r(self.next(7) as i64 - 3)).collect()
    }

    fn coefficient(&mut self) -> BivariatePolynomial {
        BivariatePolynomial::new(vec![self.polynomial(2), self.polynomial(1)])
    }

    fn operator(&mut self) -> OreOperator {
        let mut op = OreOperator::single(2, 0, BivariatePolynomial::one());
        for _ in 0..4 {
            let key = (self.next(3) as usize, self.next(3) as usize);
            op.add_term(key, self.coefficient());
        }
        op
    }

    fn laurent_operator(&mut self) -> LaurentOreOperator {
        let mut op = LaurentOreOperator::default();
        for _ in 0..3 {
            let key = (self.next(2) as usize, self.next(3) as usize);
            let valuation = self.next(5) as i64 - 3;
            op.add_term(
                key,
                LaurentBivariatePolynomial::new(valuation, self.coefficient()),
            );
        }
        op
    }

    fn laurent_rows(&mut self, count: usize) -> Vec<(i64, Vec<Rational>)> {
        (0..count)
            .map(|n| (self.next(5) as i64 - 2, self.polynomial(n % 3 + 1)))
            .collect()
    }
}

/// `(v, p)` normalised: `p(0) != 0`, zero is `(0, [0])`.
fn laurent(valuation: i64, row: &[Rational]) -> (i64, Vec<Rational>) {
    let mut row = row.to_vec();
    while row.len() > 1 && row.last().is_some_and(|c| *c == r(0)) {
        row.pop();
    }
    match row.iter().position(|c| *c != r(0)) {
        None => (0, vec![r(0)]),
        Some(e) => (valuation + e as i64, row[e..].to_vec()),
    }
}

/// `x^d K(1/x)` for a Laurent row `K`.
fn reverse_laurent(row: &(i64, Vec<Rational>), d: i64) -> (i64, Vec<Rational>) {
    let (v, p) = row;
    let m = p.len() as i64 - 1;
    let reversed: Vec<Rational> = p.iter().rev().cloned().collect();
    laurent(d - v - m, &reversed)
}

fn is_zero(p: &[Rational]) -> bool {
    p.iter().all(|c| *c == r(0))
}

fn assert_annihilates(op: &OreOperator, rows: &[Vec<Rational>], origin: usize) {
    for n in origin..origin + rows.len() - op.max_shift() {
        let value = op.apply(rows, origin, n).unwrap();
        assert!(is_zero(&value), "base {n}: {value:?}");
    }
}

/// Rows `F_origin, ...` of `F_(n+1) = A F_n` for a shift-free `A`.
fn first_order_rows(
    step: &OreOperator,
    first: Vec<Rational>,
    origin: usize,
    count: usize,
) -> Vec<Vec<Rational>> {
    let mut rows = vec![first];
    while rows.len() < count {
        let n = origin + rows.len() - 1;
        let next = step.apply_row(rows.last().unwrap(), n);
        rows.push(next);
    }
    rows
}

fn first_order_parent(step: &OreOperator) -> OreOperator {
    OreOperator::single(1, 0, BivariatePolynomial::one()).subtract(step)
}

/// Eulerian: `A_(n+1) = (1 + n x) A_n + x (1 - x) A_n'`.
fn eulerian_step() -> OreOperator {
    let mut step = OreOperator::single(0, 0, biv(&[&[1], &[0, 1]]));
    step.add_term((0, 1), biv(&[&[0, 1, -1]]));
    step
}

fn mobius(a: i64, b: i64, c: i64, e: i64, scale: Rational, degree: Vec<Rational>) -> Substitution {
    Substitution::Mobius(Box::new(MobiusMap {
        a: r(a),
        b: r(b),
        c: r(c),
        e: r(e),
        scale,
        degree,
    }))
}

#[test]
fn laurent_composition_is_the_composition_of_actions() {
    let mut random = Random(41);
    for _ in 0..10 {
        let a = random.laurent_operator();
        let b = random.laurent_operator();
        let c = random.laurent_operator();
        assert_eq!(a.compose(&b).compose(&c), a.compose(&b.compose(&c)));
        let rows = random.laurent_rows(8);
        let inner: Vec<_> = (0..7).map(|n| b.apply(&rows, 0, n).unwrap()).collect();
        for n in 0..6 {
            assert_eq!(
                a.apply(&inner, 0, n).unwrap(),
                a.compose(&b).apply(&rows, 0, n).unwrap()
            );
        }
        // Clearing: the minimal power gives a polynomial operator, one less
        // does not, and the cleared operator is x^N times the original.
        let (power, cleared) = a.clear_denominators();
        if power > 0 {
            assert!(a.clear_with_power(power - 1).is_none());
        }
        let back = LaurentOreOperator::from_ore_operator(&cleared)
            .multiply_by_variable_power(-(power as i64));
        assert_eq!(back, a);
    }
    // D x^(-1) = x^(-1) D - x^(-2).
    let d = LaurentOreOperator::single(0, 1, LaurentBivariatePolynomial::monomial(0, 0, r(1)));
    let inverse =
        LaurentOreOperator::single(0, 0, LaurentBivariatePolynomial::monomial(0, -1, r(1)));
    let mut expected =
        LaurentOreOperator::single(0, 1, LaurentBivariatePolynomial::monomial(0, -1, r(1)));
    expected.add_term((0, 0), LaurentBivariatePolynomial::monomial(0, -2, r(-1)));
    assert_eq!(d.compose(&inverse), expected);
}

#[test]
fn laurent_reversal_agrees_with_the_existing_reversal() {
    let mut random = Random(7);
    for h in 0..3 {
        let degree = vec![r(4), r(h)];
        let step = RowMapStep::Reversal {
            degree: degree.clone(),
        };
        for _ in 0..5 {
            let m = random.operator();
            let bounds = PullbackBounds::of(&m);
            let (existing, base) = pull_back(&m, std::slice::from_ref(&step), bounds).unwrap();
            assert_eq!(base, vec![r(1)]);
            let substitution = Substitution::Reversal {
                degree: degree.clone(),
            };
            let new = pull_back_substitution(&m, &substitution, bounds).unwrap();
            assert_eq!(new, existing);
            // The conjugate itself, cleared with N = b + h R.
            let conjugate = LaurentOreOperator::from_ore_operator(&m)
                .pull_back_through_reversal(&degree)
                .unwrap();
            let clearing = bounds.variable_degree + h as usize * bounds.top_shift;
            assert_eq!(conjugate.clear_with_power(clearing).unwrap(), existing);
            // The forward rules agree too.
            assert_eq!(
                push_forward_substitution(&m, &substitution).unwrap(),
                push_forward(&m, &step).unwrap()
            );
        }
    }
}

#[test]
fn laurent_conjugations_are_exact_on_laurent_rows() {
    let mut random = Random(19);
    for h in -2..=2 {
        let degree = vec![r(3), r(h)];
        let m = random.operator();
        let embedded = LaurentOreOperator::from_ore_operator(&m);
        let reversal = embedded.pull_back_through_reversal(&degree).unwrap();
        let gauge = embedded.pull_back_through_power_gauge(&degree).unwrap();
        let dilation = embedded.pull_back_through_dilation(&q(-3, 2)).unwrap();
        let rows = random.laurent_rows(8);
        let d = |n: usize| 3 + h * n as i64;
        let reversed: Vec<_> = rows
            .iter()
            .enumerate()
            .map(|(n, row)| reverse_laurent(row, d(n)))
            .collect();
        let gauged: Vec<_> = rows
            .iter()
            .enumerate()
            .map(|(n, (v, p))| laurent(v + d(n), p))
            .collect();
        let dilated: Vec<_> = rows
            .iter()
            .map(|(v, p)| {
                let mut factor = Rational::from_integer(1.into());
                let scaled: Vec<Rational> = p
                    .iter()
                    .map(|c| {
                        let value = c * &factor;
                        factor *= q(-3, 2);
                        value
                    })
                    .collect();
                let outer = (0..v.unsigned_abs()).fold(r(1), |acc, _| {
                    if *v >= 0 {
                        acc * q(-3, 2)
                    } else {
                        acc * q(-2, 3)
                    }
                });
                let scaled: Vec<Rational> = scaled.iter().map(|c| c * &outer).collect();
                laurent(*v, &scaled)
            })
            .collect();
        for n in 0..6 {
            // (M C F)_n = C_n((C^(-1) M C F)_n).
            let left = embedded.apply(&reversed, 0, n).unwrap();
            let inner = reversal.apply(&rows, 0, n).unwrap();
            assert_eq!(left, reverse_laurent(&inner, d(n)), "reversal h = {h}");
            let left = embedded.apply(&gauged, 0, n).unwrap();
            let (v, p) = gauge.apply(&rows, 0, n).unwrap();
            assert_eq!(left, laurent(v + d(n), &p), "gauge h = {h}");
            let left = embedded.apply(&dilated, 0, n).unwrap();
            let inner = dilation.apply(&rows, 0, n).unwrap();
            let mut factor = r(1);
            let scaled: Vec<Rational> = inner
                .1
                .iter()
                .map(|c| {
                    let value = c * &factor;
                    factor *= q(-3, 2);
                    value
                })
                .collect();
            let outer = (0..inner.0.unsigned_abs()).fold(r(1), |acc, _| {
                if inner.0 >= 0 {
                    acc * q(-3, 2)
                } else {
                    acc * q(-2, 3)
                }
            });
            let scaled: Vec<Rational> = scaled.iter().map(|c| c * &outer).collect();
            assert_eq!(left, laurent(inner.0, &scaled), "dilation");
        }
    }
    assert!(LaurentOreOperator::identity()
        .pull_back_through_power_gauge(&[q(1, 2)])
        .is_err());
}

#[test]
fn mobius_pullback_agrees_with_the_existing_chain() {
    let mut random = Random(29);
    for (a, b, c, e) in [(1, 0, 1, 1), (2, -1, 3, 5), (0, 1, 1, 0), (-1, 2, 1, -3)] {
        for h in 0..2 {
            let degree = vec![r(5), r(h)];
            let m = random.operator();
            let bounds = PullbackBounds::of(&m);
            let (a_q, b_q, c_q, e_q) = (r(a), r(b), r(c), r(e));
            let kappa = (&b_q * &c_q - &a_q * &e_q) / &c_q;
            let chain = [
                RowMapStep::Affine {
                    slope: kappa,
                    intercept: &a_q / &c_q,
                },
                RowMapStep::Reversal {
                    degree: degree.clone(),
                },
                RowMapStep::Affine {
                    slope: c_q,
                    intercept: e_q,
                },
            ];
            let (existing, _) = pull_back(&m, &chain, bounds).unwrap();
            let substitution = mobius(a, b, c, e, r(1), degree.clone());
            assert_eq!(
                pull_back_substitution(&m, &substitution, bounds).unwrap(),
                existing
            );
            if (a, b, c, e) == (1, 0, 1, 1) {
                assert_eq!(
                    pull_back_substitution(&m, &Substitution::FPolynomial { degree }, bounds)
                        .unwrap(),
                    existing
                );
            }
        }
    }
}

#[test]
fn substitution_rows_and_inverses() {
    let mut random = Random(31);
    let substitutions = [
        Substitution::Affine {
            slope: q(2, 3),
            intercept: r(-1),
        },
        Substitution::Reversal {
            degree: vec![r(6), r(0)],
        },
        Substitution::FPolynomial { degree: vec![r(6)] },
        mobius(2, -1, 3, 5, q(-1, 2), vec![r(6)]),
        mobius(2, -1, 0, 5, q(3, 1), vec![r(6)]),
        Substitution::PowerGauge {
            power: vec![r(2), r(1)],
        },
    ];
    for substitution in &substitutions {
        substitution.validate().unwrap();
        for n in 0..4 {
            let row = random.polynomial(5);
            let image = substitution.apply_row(&row, n).unwrap();
            let back = substitution.inverse().apply_row(&image, n).unwrap();
            let mut expected = row.clone();
            while expected.len() > 1 && expected.last() == Some(&r(0)) {
                expected.pop();
            }
            assert_eq!(back, expected, "{substitution:?}");
        }
    }
    // f-polynomial: (1 + x)^2 p(x/(1 + x)) for p = 1 + 2x + x^2 is 1 + 4x + 4x^2.
    assert_eq!(
        Substitution::FPolynomial { degree: vec![r(2)] }
            .apply_row(&[r(1), r(2), r(1)], 0)
            .unwrap(),
        vec![r(1), r(4), r(4)]
    );
    // Dividing by x needs divisibility; reversals need the degree bound.
    let divide = Substitution::PowerGauge { power: vec![r(-1)] };
    assert_eq!(divide.apply_row(&[r(1), r(1)], 0), None);
    assert_eq!(
        divide.apply_row(&[r(0), r(1), r(1)], 0),
        Some(vec![r(1), r(1)])
    );
    assert_eq!(
        Substitution::Reversal { degree: vec![r(1)] }.apply_row(&[r(0), r(0), r(1)], 0),
        None
    );
    assert!(mobius(1, 2, 2, 4, r(1), vec![r(1)]).validate().is_err());
    assert!(Substitution::Reversal {
        degree: vec![q(1, 2)]
    }
    .validate()
    .is_err());
}

fn found(
    search: DifferenceImageAnnihilatorSearch,
) -> ore_algebra::DifferenceImageAnnihilatorCertificate {
    match search {
        DifferenceImageAnnihilatorSearch::Found {
            certificate,
            method,
        } => {
            assert_eq!(method, DifferenceImageAnnihilatorMethod::PushForward);
            *certificate
        }
        other => panic!("not found: {other:?}"),
    }
}

#[test]
fn forward_rules_with_certificates() {
    let budget = DifferenceImageAnnihilatorBudget::default();
    let parent = first_order_parent(&eulerian_step());
    let embedded_parent = DifferentialDifferenceOperator::from_ore_operator(&parent);
    // A_n for n >= 1 has degree n - 1: start at index 1.
    let rows = first_order_rows(&eulerian_step(), vec![r(1)], 1, 40);
    let substitutions = [
        Substitution::FPolynomial {
            degree: vec![r(-1), r(1)],
        },
        Substitution::Reversal {
            degree: vec![r(-1), r(1)],
        },
        mobius(2, -1, 3, 5, q(-1, 2), vec![r(1), r(1)]),
        Substitution::Affine {
            slope: r(2),
            intercept: r(1),
        },
        Substitution::PowerGauge {
            power: vec![r(1), r(1)],
        },
    ];
    for substitution in substitutions {
        let image: Vec<Vec<Rational>> = rows
            .iter()
            .enumerate()
            .map(|(j, row)| substitution.apply_row(row, j + 1).unwrap())
            .collect();
        let forward = push_forward_substitution(&parent, &substitution).unwrap();
        assert_annihilates(&forward, &image, 1);
        let steps = [DifferenceRowMapStep::Substitution(substitution.clone())];
        let certificate =
            found(find_difference_image_annihilator(&embedded_parent, &steps, &budget).unwrap());
        assert_eq!(certificate.verify().unwrap(), vec![r(1)]);
        let annihilator = certificate.annihilator.to_ore_operator().unwrap();
        assert_annihilates(&annihilator, &image, 1);
        // A different map does not replay.
        let mut tampered = certificate.clone();
        tampered.steps = vec![DifferenceRowMapStep::Substitution(Substitution::Affine {
            slope: r(3),
            intercept: r(1),
        })];
        assert!(tampered.verify().is_err(), "{substitution:?}");
    }
}

#[test]
fn power_gauges_divide_rows_by_powers_of_x() {
    let budget = DifferenceImageAnnihilatorBudget::default();
    // F_n = x^(n+1) (1 + x)^n: F_(n+1) = x (1 + x) F_n.
    let step = OreOperator::single(0, 0, biv(&[&[0, 1, 1]]));
    let parent = first_order_parent(&step);
    let rows = first_order_rows(&step, vec![r(0), r(1)], 0, 40);
    // G_n = F_n / x^(n+1) = (1 + x)^n.
    let divide = Substitution::PowerGauge {
        power: vec![r(-1), r(-1)],
    };
    let image: Vec<Vec<Rational>> = rows
        .iter()
        .enumerate()
        .map(|(n, row)| divide.apply_row(row, n).unwrap())
        .collect();
    assert_eq!(image[3], vec![r(1), r(3), r(3), r(1)]);
    let forward = push_forward_substitution(&parent, &divide).unwrap();
    assert_annihilates(&forward, &image, 0);
    let certificate = found(
        find_difference_image_annihilator(
            &DifferentialDifferenceOperator::from_ore_operator(&parent),
            &[DifferenceRowMapStep::Substitution(divide)],
            &budget,
        )
        .unwrap(),
    );
    certificate.verify().unwrap();
    // The normalised rule is G_(n+1) = (1 + x) G_n.
    let expected = first_order_parent(&OreOperator::single(0, 0, biv(&[&[1, 1]])));
    assert_eq!(
        certificate.annihilator,
        DifferentialDifferenceOperator::from_ore_operator(&expected)
    );
    // A relation with a Laurent coefficient: x^(-1) F_(n+1) = (1 + x) F_n,
    // cleared to the polynomial certificate.
    let mut relation =
        LaurentOreOperator::single(1, 0, LaurentBivariatePolynomial::monomial(0, -1, r(1)));
    relation.add_term(
        (0, 0),
        LaurentBivariatePolynomial::from_polynomial(&biv(&[&[-1, -1]])),
    );
    let (power, cleared) = relation.clear_denominators();
    assert_eq!(power, 1);
    assert_eq!(cleared, parent);
}

#[test]
fn reversal_with_a_decreasing_degree() {
    // F_0 = (1 + x)^30, F_(n+1) = F_n': F_n = (30)_n (1 + x)^(30 - n), and
    // the reversal at d_n = 30 - n (h = -1, beyond the existing reversal
    // step).
    let step = OreOperator::single(0, 1, BivariatePolynomial::one());
    let parent = first_order_parent(&step);
    let mut first = vec![r(1)];
    for _ in 0..30 {
        let mut next = vec![r(0); first.len() + 1];
        for (j, c) in first.iter().enumerate() {
            next[j] += c;
            next[j + 1] += c;
        }
        first = next;
    }
    let rows = first_order_rows(&step, first, 0, 31);
    let reversal = Substitution::Reversal {
        degree: vec![r(30), r(-1)],
    };
    let image: Vec<Vec<Rational>> = rows
        .iter()
        .enumerate()
        .map(|(n, row)| reversal.apply_row(row, n).unwrap())
        .collect();
    let forward = push_forward_substitution(&parent, &reversal).unwrap();
    assert_annihilates(&forward, &image, 0);
    assert!(RowMapStep::Reversal {
        degree: vec![r(30), r(-1)]
    }
    .validate()
    .is_err());
}
