//! Bounded realization and nonnegative triangular state-basis search.

use super::*;
use std::sync::OnceLock;

fn rational_identity(n: usize) -> RationalMatrix {
    (0..n)
        .map(|i| (0..n).map(|j| Q::from_integer((i == j).into())).collect())
        .collect()
}
fn rational_polynomials(matrix: &RationalMatrix) -> Matrix {
    matrix
        .iter()
        .map(|r| r.iter().map(|c| Poly::new(vec![c.clone()])).collect())
        .collect()
}
fn permutations(n: usize) -> Vec<Vec<usize>> {
    fn visit(n: usize, p: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if p.len() == n {
            out.push(p.clone());
            return;
        }
        for i in 0..n {
            if !p.contains(&i) {
                p.push(i);
                visit(n, p, out);
                p.pop();
            }
        }
    }
    let mut out = Vec::new();
    visit(n, &mut Vec::new(), &mut out);
    out
}
fn triangular_bases(n: usize) -> Vec<RationalMatrix> {
    let mut out = vec![rational_identity(n)];
    for lower in [false, true] {
        let positions: Vec<_> = (0..n)
            .flat_map(|i| {
                (0..n)
                    .filter(move |&j| if lower { i > j } else { i < j })
                    .map(move |j| (i, j))
            })
            .collect();
        for &(i, j) in &positions {
            for value in [
                Q::new(1.into(), 2.into()),
                Q::one(),
                Q::from_integer(2.into()),
            ] {
                let mut t = rational_identity(n);
                t[i][j] = value;
                out.push(t);
            }
        }
        // Dense 0/1 bases, truncated to two off-diagonal entries in dimension 5.
        for mask in 0..1usize << positions.len() {
            if mask.count_ones() < 2 || (n > 4 && mask.count_ones() > 2) {
                continue;
            }
            let mut t = rational_identity(n);
            for (k, &(i, j)) in positions.iter().enumerate() {
                if mask & (1 << k) != 0 {
                    t[i][j] = Q::one();
                }
            }
            out.push(t);
        }
    }
    out
}
fn transformed(
    matrix: &Matrix,
    seed: &[Poly],
    projection: &[Q],
    permutation: &[usize],
    basis: &RationalMatrix,
) -> Option<(Matrix, Vec<Poly>, Vec<Q>)> {
    let n = matrix.len();
    if permutation.len() != n
        || permutation.iter().copied().collect::<BTreeSet<_>>() != (0..n).collect()
        || basis.len() != n
        || basis.iter().any(|r| r.len() != n)
    {
        return None;
    }
    if (0..n).any(|i| basis[i][i] != Q::one())
        || basis.iter().flatten().any(Signed::is_negative)
        || ![false, true].iter().any(|&lower| {
            (0..n).all(|i| {
                (0..n).all(|j| {
                    if lower {
                        i >= j || basis[i][j].is_zero()
                    } else {
                        i <= j || basis[i][j].is_zero()
                    }
                })
            })
        })
    {
        return None;
    }
    let inverse = exact_algebra::linear_systems::invert_rational_matrix(basis).ok()?;
    let ordered: Matrix = permutation
        .iter()
        .map(|&i| permutation.iter().map(|&j| matrix[i][j].clone()).collect())
        .collect();
    let seed: Vec<_> = permutation.iter().map(|&i| seed[i].clone()).collect();
    let projection: Vec<_> = permutation.iter().map(|&i| projection[i].clone()).collect();
    let new_projection: Vec<_> = (0..n)
        .map(|j| (0..n).fold(Q::zero(), |s, i| s + &projection[i] * &inverse[i][j]))
        .collect();
    Some((
        multiply(
            &multiply(&rational_polynomials(basis), &ordered),
            &rational_polynomials(&inverse),
        ),
        action(&rational_polynomials(basis), &seed),
        new_projection,
    ))
}
pub(super) fn certify_basis(
    matrix: &Matrix,
    seed: &[Poly],
    projection: &[Q],
    permutation: &[usize],
    basis: &RationalMatrix,
    kind: &str,
    declarations: &[String],
) -> Option<MatrixOrbitCertificate> {
    let (g, v, p) = transformed(matrix, seed, projection, permutation, basis)?;
    if p.iter().any(Signed::is_negative) || !interlacing_seed(&v) {
        return None;
    }
    let pencils = affine::certify(&g)?;
    let supports = nonzero_supports(&g, &v, &p)?;
    Some(MatrixOrbitCertificate {
        state_system: encode(matrix, seed, projection),
        permutation: permutation.to_vec(),
        basis: basis.iter().map(|r| strings(r)).collect(),
        transformed_system: encode(&g, &v, &p),
        pencils,
        nonzero_supports: supports,
        realization: kind.into(),
        lean_declarations: declarations.to_vec(),
    })
}
pub(super) fn search_bases(
    matrix: &Matrix,
    seed: &[Poly],
    projection: &[Q],
    kind: &str,
    declarations: &[String],
    remaining: &mut usize,
) -> Option<MatrixOrbitCertificate> {
    if matrix
        .iter()
        .flatten()
        .any(|p| p.degree().is_some_and(|d| d > 1))
    {
        return None;
    }
    let bases = triangular_bases(matrix.len());
    for order in permutations(matrix.len()) {
        for basis in &bases {
            if *remaining == 0 {
                return None;
            }
            *remaining -= 1;
            if let Some(c) =
                certify_basis(matrix, seed, projection, &order, basis, kind, declarations)
            {
                return Some(c);
            }
        }
    }
    None
}
#[derive(Clone)]
struct Shape {
    kind: &'static str,
    matrix: Matrix,
    multipliers: Vec<Poly>,
}
fn characteristic_multipliers(matrix: &Matrix) -> Vec<Poly> {
    let n = matrix.len();
    let mut b = identity(n);
    let mut out = Vec::new();
    for k in 1..=n {
        b = multiply(matrix, &b);
        let coefficient = (0..n).fold(Poly::zero(), |s, i| s + b[i][i].clone())
            * Poly::new(vec![Q::new((-1).into(), k.into())]);
        out.push(-coefficient.clone());
        for (i, row) in b.iter_mut().enumerate() {
            row[i] = row[i].clone() + coefficient.clone();
        }
    }
    out
}
fn digits(mut code: usize, base: usize, n: usize) -> Vec<usize> {
    (0..n)
        .map(|_| {
            let d = code % base;
            code /= base;
            d
        })
        .collect()
}
fn templates(n: usize) -> Vec<Shape> {
    let zero = Poly::zero();
    let one = Poly::one();
    let x = Poly::new(vec![Q::zero(), Q::one()]);
    let markers = [one.clone(), zero.clone(), one.clone() + x.clone()];
    let mut shapes = Vec::new();
    let mut seen = BTreeSet::new();
    let mut push = |kind, matrix: Matrix| {
        let key: Vec<_> = matrix
            .iter()
            .flatten()
            .map(|p| strings(p.coeffs()).join(","))
            .collect();
        if seen.insert(key) {
            let multipliers = characteristic_multipliers(&matrix);
            shapes.push(Shape {
                kind,
                matrix,
                multipliers,
            });
        }
    };
    for code in 0..(n + 1).pow(n as u32) {
        let thresholds = digits(code, n + 1, n);
        if thresholds.windows(2).any(|w| w[0] > w[1]) {
            continue;
        }
        for code in 0..3usize.pow(n as u32) {
            let marker = digits(code, 3, n);
            let matrix = (0..n)
                .map(|i| {
                    (0..n)
                        .map(|j| {
                            if j < thresholds[i] {
                                x.clone()
                            } else if j == thresholds[i] {
                                markers[marker[i]].clone()
                            } else {
                                one.clone()
                            }
                        })
                        .collect()
                })
                .collect();
            push("threshold", matrix);
        }
    }
    for h in 1..=n {
        let m = n - h;
        for acode in 0..3usize.pow(n as u32) {
            let a = digits(acode, 3, n);
            for bcode in 0..3usize.pow(m as u32) {
                let b = digits(bcode, 3, m);
                let matrix = (0..n)
                    .map(|i| {
                        (0..n)
                            .map(|j| {
                                if i + 1 < h {
                                    if i == j {
                                        Poly::new(vec![Q::from_integer(a[i].into())])
                                    } else if j == i + 1 {
                                        one.clone()
                                    } else {
                                        zero.clone()
                                    }
                                } else {
                                    let constant = if i == j {
                                        a[i]
                                    } else if i < j {
                                        a[j] + if j >= h { b[j - h] } else { 0 }
                                    } else {
                                        0
                                    };
                                    Poly::new(vec![
                                        Q::from_integer(constant.into()),
                                        if j == 0 { Q::one() } else { Q::zero() },
                                    ])
                                }
                            })
                            .collect()
                    })
                    .collect();
                push("separated_characteristic", matrix);
            }
        }
    }
    shapes
}
fn all_templates(n: usize) -> &'static [Shape] {
    static CACHE: OnceLock<[Vec<Shape>; 3]> = OnceLock::new();
    &CACHE.get_or_init(|| [templates(1), templates(2), templates(3)])[n - 1]
}
fn infer_seed(matrix: &Matrix, projection: &[Q], rows: &[Vec<Q>]) -> Option<Vec<Poly>> {
    // Constant/affine unknown coordinates; coefficient equations include every
    // prescribed row, so no high-degree residual is discarded.
    let n = matrix.len();
    let mut functional: Vec<Poly> = projection
        .iter()
        .map(|c| Poly::new(vec![c.clone()]))
        .collect();
    let mut equations = Vec::new();
    let mut rhs = Vec::new();
    for row in rows {
        let bound = functional
            .iter()
            .filter_map(Poly::degree)
            .max()
            .unwrap_or(0)
            + 1;
        for k in 0..=bound.max(row.len()) {
            equations.push(
                (0..n * 2)
                    .map(|u| {
                        let j = u / 2;
                        let d = u % 2;
                        if k < d {
                            Q::zero()
                        } else {
                            functional[j]
                                .coeffs()
                                .get(k - d)
                                .cloned()
                                .unwrap_or_default()
                        }
                    })
                    .collect(),
            );
            rhs.push(row.get(k).cloned().unwrap_or_default());
        }
        functional = (0..n)
            .map(|j| {
                (0..n).fold(Poly::zero(), |s, i| {
                    s + functional[i].clone() * matrix[i][j].clone()
                })
            })
            .collect();
    }
    let values = crate::linalg::solve_linear_system(&equations, &rhs)?;
    Some(values.chunks(2).map(|c| Poly::new(c.to_vec())).collect())
}
pub(super) fn recognize(multipliers: &[Poly], rows: &[Vec<Q>]) -> Option<MatrixOrbitCertificate> {
    let n = multipliers.len();
    if n == 0 || n > 4 {
        return None;
    }
    let mut remaining = MAXIMUM_BASIS_CANDIDATES;
    if n <= 3 {
        for shape in all_templates(n)
            .iter()
            .filter(|s| s.multipliers == multipliers)
        {
            for selected in 0..=n {
                let projection: Vec<Q> = (0..n)
                    .map(|i| Q::from_integer((selected == n || i == selected).into()))
                    .collect();
                let Some(seed) = infer_seed(&shape.matrix, &projection, rows) else {
                    continue;
                };
                if !identity_and_seeds(&shape.matrix, &seed, &projection, multipliers, rows) {
                    continue;
                }
                let declarations = lean_declarations(&shape.matrix, shape.kind);
                if let Some(c) = search_bases(
                    &shape.matrix,
                    &seed,
                    &projection,
                    shape.kind,
                    &declarations,
                    &mut remaining,
                ) {
                    return Some(c);
                }
            }
        }
    }
    let mut companion = vec![vec![Poly::zero(); n]; n];
    for (i, row) in companion.iter_mut().enumerate().take(n - 1) {
        row[i + 1] = Poly::one();
    }
    companion[n - 1] = multipliers.iter().rev().cloned().collect();
    let seed: Vec<_> = rows[..n].iter().cloned().map(Poly::new).collect();
    let mut projection = vec![Q::zero(); n];
    projection[0] = Q::one();
    if !identity_and_seeds(&companion, &seed, &projection, multipliers, rows) {
        return None;
    }
    let declarations = lean_declarations(&companion, "companion");
    search_bases(
        &companion,
        &seed,
        &projection,
        "companion",
        &declarations,
        &mut remaining,
    )
}
pub(super) fn lean_declarations(matrix: &Matrix, kind: &str) -> Vec<String> {
    if kind == "supplied" || kind == "companion" {
        let threshold = lean_declarations(matrix, "threshold");
        if !threshold.is_empty() {
            return threshold;
        }
        if matrix.len() <= 3
            && all_templates(matrix.len())
                .iter()
                .any(|s| s.kind == "separated_characteristic" && s.matrix == *matrix)
        {
            return lean_declarations(matrix, "separated_characteristic");
        }
        return Vec::new();
    }
    if kind == "separated_characteristic" {
        return vec!["RealRooted.separatedCharacteristicRows_isPFPolynomial".into()];
    }
    if kind != "threshold" {
        return Vec::new();
    }
    let n = matrix.len();
    let x = Poly::new(vec![Q::zero(), Q::one()]);
    let one = Poly::one();
    let mut data = Vec::new();
    for row in matrix {
        let t = row.iter().take_while(|p| **p == x).count();
        let marker = if t == n { one.clone() } else { row[t].clone() };
        if row.iter().skip(t + 1).any(|p| *p != one) {
            return Vec::new();
        }
        data.push((t, marker));
    }
    if data.windows(2).any(|w| w[0].0 > w[1].0) {
        return Vec::new();
    }
    for (special, family) in [
        (Poly::zero(), "GustafssonSolus"),
        (one.clone() + x, "HaglundZhang"),
    ] {
        if data.iter().all(|(_, m)| *m == one || *m == special)
            && (0..n).all(|i| {
                (i + 1..n)
                    .all(|j| data[i].0 != data[j].0 || data[i].1 != special || data[j].1 == special)
            })
        {
            return vec![format!("RealRooted.{family}.orbit_certificate")];
        }
    }
    Vec::new()
}
pub(super) fn valid_provenance(matrix: &Matrix, kind: &str, supplied: bool) -> bool {
    if supplied {
        return kind == "supplied";
    }
    match kind {
        "threshold" | "separated_characteristic" if matrix.len() <= 3 => {
            all_templates(matrix.len())
                .iter()
                .any(|s| s.kind == kind && s.matrix == *matrix)
        }
        "companion" => {
            let n = matrix.len();
            n <= 4
                && matrix.iter().take(n - 1).enumerate().all(|(i, r)| {
                    r.iter().enumerate().all(|(j, p)| {
                        *p == if j == i + 1 {
                            Poly::one()
                        } else {
                            Poly::zero()
                        }
                    })
                })
        }
        _ => false,
    }
}
