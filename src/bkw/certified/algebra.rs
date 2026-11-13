//! Exact polynomial resultants and sign cells, shared by the certificates.

use super::*;
use crate::sturm::SturmChain;

pub(super) fn chain(p: &P) -> SturmChain {
    SturmChain::from_bigint_coeffs(&vec_poly::clear_denominators(p.coeffs()))
}

pub(super) fn counts(p: &P) -> RootCounts {
    let sturm = chain(p);
    RootCounts {
        degree: p.degree(),
        squarefree_degree: sturm.square_free_degree(),
        distinct_real_roots: sturm.count_real_roots(),
    }
}

/// Sylvester determinant over Q[x], using fraction-free elimination.
pub(super) fn resultant(f: &[P], g: &[P]) -> P {
    let mut f = f.to_vec();
    let mut g = g.to_vec();
    while f.last().is_some_and(P::is_zero) {
        f.pop();
    }
    while g.last().is_some_and(P::is_zero) {
        g.pop();
    }
    if f.is_empty() || g.is_empty() {
        return P::zero();
    }
    let (n, m) = (f.len() - 1, g.len() - 1);
    if n == 0 {
        return f[0].pow(m);
    }
    if m == 0 {
        return g[0].pow(n);
    }
    let size = n + m;
    let mut matrix = vec![vec![P::zero(); size]; size];
    for i in 0..m {
        for j in 0..=n {
            matrix[i][i + j] = f[n - j].clone();
        }
    }
    for i in 0..n {
        for j in 0..=m {
            matrix[m + i][i + j] = g[m - j].clone();
        }
    }
    let mut previous = P::one();
    let mut sign = Q::one();
    for k in 0..size - 1 {
        let Some(pivot_row) = (k..size).find(|&i| !matrix[i][k].is_zero()) else {
            return P::zero();
        };
        if pivot_row != k {
            matrix.swap(k, pivot_row);
            sign = -sign;
        }
        let pivot = matrix[k][k].clone();
        for i in k + 1..size {
            for j in k + 1..size {
                let numerator = matrix[i][j].clone() * pivot.clone()
                    - matrix[i][k].clone() * matrix[k][j].clone();
                matrix[i][j] = numerator.exact_div(&previous);
            }
        }
        for row in matrix.iter_mut().skip(k + 1) {
            row[k] = P::zero();
        }
        previous = pivot;
    }
    matrix[size - 1][size - 1].scale(&sign)
}

pub(super) fn derivative_in_outer(p: &[P]) -> Vec<P> {
    p.iter()
        .enumerate()
        .skip(1)
        .map(|(i, c)| c.scale(&Q::from_integer(i.into())))
        .collect()
}

/// Separate every root from its neighbors and from the rational boundaries.
pub(super) fn cells(p: &P, lower: &Q, upper: &Q) -> (Vec<(Q, Q)>, Vec<Q>) {
    let boundary = P::new(vec![-lower.clone(), Q::one()]) * P::new(vec![-upper.clone(), Q::one()]);
    let squarefree = p.squarefree_part();
    let inside = squarefree.exact_div(&squarefree.gcd(&boundary));
    let sturm = chain(&inside);
    let mut epsilon = Q::one();
    let intervals = loop {
        let roots: Vec<_> = sturm
            .isolate_roots(&epsilon)
            .into_iter()
            .filter(|(a, b)| b > lower && a < upper)
            .collect();
        if roots.iter().all(|(a, b)| a > lower && b < upper)
            && roots.windows(2).all(|w| w[0].1 < w[1].0)
        {
            break roots;
        }
        epsilon /= q(2);
    };
    let mut samples = Vec::new();
    let mut left = lower.clone();
    for (a, b) in &intervals {
        samples.push((&left + a) / q(2));
        left = b.clone();
    }
    samples.push((&left + upper) / q(2));
    (intervals, samples)
}

pub(super) fn nonnegative(p: &P, lower: &Q, upper: &Q) -> bool {
    if p.is_zero() {
        return true;
    }
    let (_, samples) = cells(p, lower, upper);
    let result = samples
        .iter()
        .chain([lower, upper])
        .all(|x| !p.evaluate(x).is_negative());
    result
}
