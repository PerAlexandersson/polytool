//! The pullback contract of every row-map step on arbitrary rows, the
//! Ore rules, the Euler form and right pseudo-division.

use num_traits::{One, Zero};
use ore_algebra::{
    from_euler_form, pull_back, right_pseudo_division, to_euler_form, BivariatePolynomial,
    OreOperator, PullbackBounds, Rational, RowMapStep,
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

/// Pseudo-random rows `F_0, ..., F_8`, `deg F_n <= n + 1`.
fn rows(seed: u64) -> Vec<Vec<Rational>> {
    let mut state = seed;
    (0..9)
        .map(|n| {
            (0..n + 2)
                .map(|_| {
                    state = state
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    r(((state >> 33) % 7) as i64 - 3)
                })
                .collect()
        })
        .collect()
}

fn trim(p: &[Rational]) -> Vec<Rational> {
    let mut p = p.to_vec();
    while p.len() > 1 && p.last().is_some_and(Zero::is_zero) {
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

/// `M = (n + 2) S^2 - (1 + x) D S + x^2 (n + 1)`.
fn operator() -> OreOperator {
    let mut m = OreOperator::single(2, 0, biv(&[&[2], &[1]]));
    m.add_term((1, 1), biv(&[&[-1, -1]]));
    m.add_term((0, 0), biv(&[&[0, 0, 1], &[0, 0, 1]]));
    m
}

fn images(step: &RowMapStep, rows: &[Vec<Rational>]) -> Vec<Vec<Rational>> {
    rows.iter()
        .enumerate()
        .map(|(n, row)| step.apply_row(row, n).unwrap())
        .collect()
}

#[test]
fn ore_rules_and_euler_round_trip() {
    let d = OreOperator::single(0, 1, BivariatePolynomial::one());
    let x = OreOperator::single(0, 0, biv(&[&[0, 1]]));
    let mut expected = OreOperator::single(0, 1, biv(&[&[0, 1]]));
    expected.add_term((0, 0), BivariatePolynomial::one());
    assert_eq!(d.compose(&x), expected);
    let s = OreOperator::single(1, 0, BivariatePolynomial::one());
    let n = OreOperator::single(0, 0, biv(&[&[0], &[1]]));
    assert_eq!(s.compose(&n), OreOperator::single(1, 0, biv(&[&[1], &[1]])));
    // x^J L = Euler form, and back.
    let m = operator();
    let (j, euler) = to_euler_form(&m);
    assert_eq!(j, 1);
    assert_eq!(from_euler_form(&euler), m.multiply_by_variable_power(1));
}

#[test]
fn pullbacks_satisfy_their_contracts() {
    let m = operator();
    let bounds = PullbackBounds::of(&m);
    let source = rows(11);
    // Weyl operator: (M G)_n = (E F)_n.
    let mut t = OreOperator::single(0, 1, biv(&[&[0, 1], &[1]]));
    t.add_term((0, 0), biv(&[&[1], &[0, 2]]));
    let step = RowMapStep::WeylOperator(t);
    let (e, b) = pull_back(&m, std::slice::from_ref(&step), bounds).unwrap();
    assert_eq!(b, vec![r(1)]);
    let image = images(&step, &source);
    for n in 0..6 {
        assert_eq!(m.apply(&image, 0, n), e.apply(&source, 0, n));
    }
    // Affine: (M G)_n(x) = (E F)_n(a x + c).
    let step = RowMapStep::Affine {
        slope: r(2),
        intercept: r(-1),
    };
    let (e, _) = pull_back(&m, std::slice::from_ref(&step), bounds).unwrap();
    let image = images(&step, &source);
    for n in 0..6 {
        let direct = m.apply(&image, 0, n).unwrap();
        let pulled = e.apply(&source, 0, n).unwrap();
        for point in [r(0), r(1), r(3)] {
            assert_eq!(
                evaluate(&direct, &point),
                evaluate(&pulled, &(r(2) * &point - r(1)))
            );
        }
    }
    // Gauge mu(n)/mu(n-1) = (n + 2)/3, mu(0) = 1:
    // (M G)_n b(n) = mu(n) (E F)_n.
    let alpha = vec![r(2), r(1)];
    let beta = vec![r(3)];
    let step = RowMapStep::Gauge {
        numerator: alpha.clone(),
        denominator: beta.clone(),
    };
    let (e, b) = pull_back(&m, std::slice::from_ref(&step), bounds).unwrap();
    let mut mu = vec![r(1)];
    for n in 1..source.len() {
        let next =
            mu[n - 1].clone() * evaluate(&alpha, &r(n as i64)) / evaluate(&beta, &r(n as i64));
        mu.push(next);
    }
    let image: Vec<_> = source
        .iter()
        .zip(&mu)
        .map(|(row, mu)| row.iter().map(|c| c * mu).collect::<Vec<_>>())
        .collect();
    for (n, mu_n) in mu.iter().enumerate().take(6) {
        let direct = m.apply(&image, 0, n).unwrap();
        let pulled = e.apply(&source, 0, n).unwrap();
        let factor = evaluate(&b, &r(n as i64));
        let left: Vec<_> = direct.iter().map(|c| c * &factor).collect();
        let right: Vec<_> = pulled.iter().map(|c| c * mu_n).collect();
        assert_eq!(trim(&left), trim(&right));
    }
    // Reversal d_n = n + 1: (M G)_n = x^(d_n + N) (E F)_n(1/x).
    let degree = vec![r(1), r(1)];
    let step = RowMapStep::Reversal {
        degree: degree.clone(),
    };
    let (e, _) = pull_back(&m, std::slice::from_ref(&step), bounds).unwrap();
    let top = bounds.variable_degree + bounds.top_shift;
    let image = images(&step, &source);
    for n in 0..6 {
        let direct = trim(&m.apply(&image, 0, n).unwrap());
        let pulled = trim(&e.apply(&source, 0, n).unwrap());
        let total = n + 1 + top;
        assert!(pulled.len() <= total + 1);
        let mut reversed = vec![r(0); total + 1];
        for (k, c) in pulled.iter().enumerate() {
            reversed[total - k] = c.clone();
        }
        assert_eq!(trim(&reversed), direct, "base {n}");
    }
    // Coefficient multiplier u(k) = 1/(k + 1)!... : u(k + 1)/u(k) = 1/(k + 2):
    // [x^m] (E F)_n = V(m)/u(m) [x^m] x^J (M G)_n with V = 1.
    let step = RowMapStep::CoefficientMultiplier {
        numerator: vec![r(1)],
        denominator: vec![r(2), r(1)],
    };
    let (e, _) = pull_back(&m, std::slice::from_ref(&step), bounds).unwrap();
    let image = images(&step, &source);
    let u = |k: usize| (0..k).fold(r(1), |acc, j| acc / r(j as i64 + 2));
    for n in 0..6 {
        let direct = m.apply(&image, 0, n).unwrap();
        let pulled = e.apply(&source, 0, n).unwrap();
        let shifted: Vec<_> = std::iter::repeat_n(r(0), bounds.derivative_order)
            .chain(direct)
            .collect();
        let width = shifted.len().max(pulled.len());
        for k in 0..width {
            let left = pulled.get(k).cloned().unwrap_or_else(|| r(0));
            let right = shifted.get(k).cloned().unwrap_or_else(|| r(0)) / u(k);
            assert_eq!(left, right, "base {n}, coefficient {k}");
        }
    }
}

#[test]
fn right_pseudo_division_identity() {
    // L = (n + 1) S - (1 + x) - x D.
    let mut l = OreOperator::single(1, 0, biv(&[&[1], &[1]]));
    l.add_term((0, 0), biv(&[&[-1, -1]]));
    l.add_term((0, 1), biv(&[&[0, -1]]));
    let dividend = operator();
    let division = right_pseudo_division(&dividend, &l, 2).unwrap();
    assert!(division.remainder.max_shift() < 1 || division.remainder.is_zero());
    assert!(division.multiplier.is_variable_free());
    let left = dividend.left_multiply(&division.multiplier);
    let right = division.quotient.compose(&l).add(&division.remainder);
    assert_eq!(left, right);
    // A leader involving x is rejected.
    let bad = OreOperator::single(1, 0, biv(&[&[0, 1]]));
    assert!(right_pseudo_division(&dividend, &bad, 2).is_none());
    assert_eq!(BivariatePolynomial::one().leading(), Some(&Rational::one()));
}
