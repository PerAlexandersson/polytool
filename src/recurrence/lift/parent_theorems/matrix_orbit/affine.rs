//! Exact universal pencil checks for nonnegative affine polynomial entries.
//!
//! The Wronskian quadratic is nonpositive on the line exactly when its two
//! endpoint coefficients are nonpositive and its discriminant is nonpositive.
//! Homogenization reduces the parameter conditions to copositivity of two
//! rational matrices. Their simplex minima require at most seven face systems.

use super::*;

fn copositivity(matrix: RationalMatrix) -> Option<MatrixOrbitCopositivityCheck> {
    let n = matrix.len();
    let mut minimum: Option<Q> = None;
    for mask in 1..1usize << n {
        let support: Vec<_> = (0..n).filter(|i| mask & (1 << i) != 0).collect();
        let mut system: Vec<Vec<Q>> = support
            .iter()
            .map(|&i| {
                let mut row: Vec<_> = support
                    .iter()
                    .map(|&j| &matrix[i][j] * Q::from_integer(2.into()))
                    .collect();
                row.push(-Q::one());
                row
            })
            .collect();
        let mut last = vec![Q::one(); support.len()];
        last.push(Q::zero());
        system.push(last);
        let mut rhs = vec![Q::zero(); support.len()];
        rhs.push(Q::one());
        let Some(solution) = crate::linalg::solve_full_rank_square_linear_system(&system, &rhs)
        else {
            continue;
        };
        let point = &solution[..support.len()];
        if point.iter().any(Signed::is_negative) {
            continue;
        }
        let value = support.iter().enumerate().fold(Q::zero(), |v, (i, &si)| {
            v + support.iter().enumerate().fold(Q::zero(), |w, (j, &sj)| {
                w + &point[i] * &point[j] * &matrix[si][sj]
            })
        });
        if minimum.as_ref().is_none_or(|m| value < *m) {
            minimum = Some(value);
        }
    }
    let minimum = minimum?;
    if minimum.is_negative() {
        return None;
    }
    Some(MatrixOrbitCopositivityCheck {
        matrix: matrix.iter().map(|r| strings(r)).collect(),
        minimum: minimum.to_string(),
    })
}
fn coefficient(p: &Poly, degree: usize) -> Q {
    p.coeffs().get(degree).cloned().unwrap_or_default()
}
fn affine_interlacing(f: &Poly, g: &Poly) -> bool {
    let (Some(df), Some(dg)) = (f.degree(), g.degree()) else {
        return true;
    };
    (dg == df || dg == df + 1)
        && (df == 0
            || coefficient(f, 0) * coefficient(g, 1) >= coefficient(g, 0) * coefficient(f, 1))
}
fn wronskian(f: &Poly, g: &Poly) -> Q {
    coefficient(f, 1) * coefficient(g, 0) - coefficient(f, 0) * coefficient(g, 1)
}

fn pencil(
    a: &Poly,
    b: &Poly,
    c: &Poly,
    d: &Poly,
    rows: [usize; 2],
    columns: [usize; 2],
) -> Option<MatrixOrbitPencilCheck> {
    if !affine_interlacing(a, c) || !affine_interlacing(b, d) {
        return None;
    }
    let degree = |p: &Poly, q: &Poly| p.degree().map(|d| d + 1).max(q.degree());
    let (df, dg) = (degree(b, d), degree(a, c));
    if df.is_none() || dg.is_none() {
        return Some(MatrixOrbitPencilCheck {
            rows,
            columns,
            zero_side: true,
            wronskian_parameters: vec![],
            copositivity: vec![],
        });
    }
    let (df, dg) = (df?, dg?);
    if dg != df && dg != df + 1 {
        return None;
    }
    let aa = wronskian(b, a);
    let cc = wronskian(d, c);
    let kk = wronskian(b, c) + wronskian(d, a);
    let dd = b.clone() * c.clone() - d.clone() * a.clone();
    let (d0, d1, d2) = (
        coefficient(&dd, 0),
        coefficient(&dd, 1),
        coefficient(&dd, 2),
    );
    if aa.is_positive() || d0.is_positive() || d2.is_positive() {
        return None;
    }
    let two = Q::from_integer(2.into());
    let four = Q::from_integer(4.into());
    let q2 = vec![vec![-&aa, -&kk / &two], vec![-&kk / &two, -&cc]];
    let cross = (&four * &aa * &cc + &four * &d2 * &d0 - (&kk + &d1) * (&kk + &d1)) / &two;
    let q3 = vec![
        vec![&four * &aa * &d0, -&two * &aa * &d1, cross.clone()],
        vec![-&two * &aa * &d1, &four * &aa * &d2, &two * &d2 * &kk],
        vec![cross, &two * &d2 * &kk, &four * &d2 * &cc],
    ];
    Some(MatrixOrbitPencilCheck {
        rows,
        columns,
        zero_side: false,
        wronskian_parameters: strings(&[aa, cc, kk, d0, d1, d2]),
        copositivity: vec![copositivity(q2)?, copositivity(q3)?],
    })
}
pub(super) fn certify(matrix: &Matrix) -> Option<Vec<MatrixOrbitPencilCheck>> {
    if matrix
        .iter()
        .flatten()
        .any(|p| !nonnegative(p) || p.degree().is_some_and(|d| d > 1))
    {
        return None;
    }
    let n = matrix.len();
    let mut checks = Vec::new();
    for k in 0..n {
        for l in k..n {
            for i in 0..n {
                for j in i..n {
                    checks.push(pencil(
                        &matrix[k][i],
                        &matrix[k][j],
                        &matrix[l][i],
                        &matrix[l][j],
                        [k, l],
                        [i, j],
                    )?);
                }
            }
        }
    }
    Some(checks)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn q(n: i64) -> Q {
        Q::from_integer(n.into())
    }
    #[test]
    fn simplex_faces_include_singular_and_interior_cases() {
        assert!(copositivity(vec![vec![q(0); 3]; 3]).is_some());
        assert!(copositivity(vec![vec![q(1); 3]; 3]).is_some());
        assert!(copositivity(vec![vec![q(1), q(2)], vec![q(2), q(1)]]).is_some());
        let neg = Q::new((-3).into(), 4.into());
        assert!(copositivity(vec![
            vec![q(1), neg.clone(), neg.clone()],
            vec![neg.clone(), q(1), neg.clone()],
            vec![neg.clone(), neg, q(1)]
        ])
        .is_none());
    }

    #[test]
    fn every_small_affine_block_agrees_with_exact_nested_pencils() {
        let entries = [
            Poly::zero(),
            Poly::one(),
            Poly::new(vec![q(0), q(1)]),
            Poly::new(vec![q(1), q(1)]),
        ];
        let values = [Q::new(1.into(), 8.into()), q(1), q(8)];
        let mut accepted = 0;
        let mut rejected = 0;
        for a in &entries {
            for b in &entries {
                for c in &entries {
                    for d in &entries {
                        let valid = pencil(a, b, c, d, [0, 1], [0, 1]).is_some();
                        let mut witness = false;
                        for lam in &values {
                            for mu in &values {
                                for alpha in &values {
                                    for beta in &values {
                                        let l = Poly::new(vec![mu.clone(), lam.clone()]);
                                        let f = l.clone() * b.clone() + d.clone();
                                        let g = l * a.clone() + c.clone();
                                        let h =
                                            Poly::new(vec![beta.clone(), alpha.clone()]) * f + g;
                                        if !crate::is_real_rooted_bigint_coeffs(
                                            &integer_coefficients(h.coeffs()),
                                        ) {
                                            witness = true;
                                        }
                                    }
                                }
                            }
                        }
                        if valid {
                            accepted += 1;
                            assert!(!witness);
                        } else {
                            rejected += 1;
                            assert!(witness);
                        }
                    }
                }
            }
        }
        assert_eq!((accepted, rejected), (77, 179));
    }
}
