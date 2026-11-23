//! Uniform affine pencils, with exact rational rejection witnesses.
use super::*;

/// An exact failed pencil, at positive rational parameters and one integer index.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PencilRejection {
    /// Weakly ordered zero-based row indices.
    pub rows: [usize; 2],
    /// Weakly ordered zero-based column indices.
    pub columns: [usize; 2],
    /// Source index at the failure.
    pub index: usize,
    /// Positive lambda,mu of L=lambda*x+mu.
    pub parameters: [String; 2],
    /// Ascending coefficients of L*b+d.
    pub left: Vec<String>,
    /// Ascending coefficients of L*a+c.
    pub right: Vec<String>,
}

/// Symbolic certificate of one pencil on n=first_index+t, t>0.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymbolicPencilCertificate {
    /// Weakly ordered row indices.
    pub rows: [usize; 2],
    /// Weakly ordered column indices.
    pub columns: [usize; 2],
    /// Whether one side is identically zero (other side still certified splitting).
    pub zero_side: bool,
    /// A,C,K,D0,D1,D2 as ascending polynomials in t.
    pub wronskian_parameters: Vec<Vec<String>>,
    /// Copositivity certificates for every coefficient of the two parameter forms.
    pub copositivity: Vec<Vec<MatrixOrbitCopositivityCheck>>,
}

/// Exact affine-matrix test on the entire integer ray, including repeated indices.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AffinePreserverReport {
    /// First allowed source index.
    pub first_index: usize,
    /// proved, rejected (exact witness), or unresolved (symbolic sufficient test failed).
    pub status: String,
    /// Description of the supported class and remaining obligations.
    pub reason: String,
    /// Exact complete checks at the boundary n=first_index.
    pub boundary: Vec<MatrixOrbitPencilCheck>,
    /// Uniform checks on the open ray, one for every weakly ordered pair.
    pub pencils: Vec<SymbolicPencilCertificate>,
    /// Every failed block found by the bounded exact witness screen.
    pub rejections: Vec<PencilRejection>,
    /// Exact ordered-vector failures for differential transitions.
    pub operator_rejections: Vec<OperatorRejection>,
}

pub(super) fn directed(f: &[Q], g: &[Q]) -> bool {
    let fi = vp::clear_denominators(f);
    let gi = vp::clear_denominators(g);
    crate::is_real_rooted_bigint_coeffs(&fi)
        && crate::is_real_rooted_bigint_coeffs(&gi)
        && (vp::is_zero(f)
            || vp::is_zero(g)
            || crate::check_weak_interlacing_bigint_coeffs(&fi, &gi) == Some(true))
}
fn nonnegative(p: &[Q]) -> bool {
    p.iter().all(|c| !c.is_negative())
}
fn column(p: &Biv, x: usize) -> Vec<Q> {
    p.variable_coefficient(x)
}
fn neg(p: &[Q]) -> Vec<Q> {
    vp::neg(p)
}
fn mul(p: &[Q], q: &[Q]) -> Vec<Q> {
    vp::mul(p, q)
}
fn add(p: &[Q], q: &[Q]) -> Vec<Q> {
    vp::add(p, q)
}
fn sub(p: &[Q], q: &[Q]) -> Vec<Q> {
    vp::sub(p, q)
}
fn wronskian(a: &Biv, b: &Biv) -> Vec<Q> {
    sub(
        &mul(&column(a, 1), &column(b, 0)),
        &mul(&column(a, 0), &column(b, 1)),
    )
}
fn order(a: &Biv, b: &Biv) -> bool {
    let (Some(da), Some(db)) = (a.variable_degree(), b.variable_degree()) else {
        return true;
    };
    (da == db || da + 1 == db) && (da == 0 || nonnegative(&neg(&wronskian(a, b))))
}
fn coefficient_copositive(m: &[Vec<Vec<Q>>]) -> Option<Vec<MatrixOrbitCopositivityCheck>> {
    let degree = m.iter().flatten().map(Vec::len).max().unwrap_or(1);
    (0..degree)
        .map(|k| {
            affine_pencils::copositivity(
                m.iter()
                    .map(|row| {
                        row.iter()
                            .map(|p| p.get(k).cloned().unwrap_or_default())
                            .collect()
                    })
                    .collect(),
            )
        })
        .collect()
}
fn symbolic(
    a: &Biv,
    b: &Biv,
    c: &Biv,
    d: &Biv,
    rows: [usize; 2],
    columns: [usize; 2],
) -> Option<SymbolicPencilCertificate> {
    if !order(a, c) || !order(b, d) {
        return None;
    }
    let degree = |p: &Biv, q: &Biv| p.variable_degree().map(|v| v + 1).max(q.variable_degree());
    let (df, dg) = (degree(b, d), degree(a, c));
    if df.is_none() || dg.is_none() {
        return Some(SymbolicPencilCertificate {
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
    let kk = add(&wronskian(b, c), &wronskian(d, a));
    let dd = b.multiply(c).subtract(&d.multiply(a));
    let (d0, d1, d2) = (column(&dd, 0), column(&dd, 1), column(&dd, 2));
    if !nonnegative(&neg(&aa)) || !nonnegative(&neg(&d0)) || !nonnegative(&neg(&d2)) {
        return None;
    }
    let parameters = [aa, cc, kk, d0, d1, d2];
    let [q2, q3] = affine_pencils::pencil_forms(&parameters);
    Some(SymbolicPencilCertificate {
        rows,
        columns,
        zero_side: false,
        wronskian_parameters: parameters.iter().map(|p| strings(p)).collect(),
        copositivity: vec![coefficient_copositive(&q2)?, coefficient_copositive(&q3)?],
    })
}

/// Certify nonnegative affine-in-x entries and all quantified Brändén
/// pencils, including repeated rows/columns, for every integer n>=first.
/// Constant matrices use an exact decision. For polynomial n-dependence,
/// copositivity of each shifted coefficient matrix is a sufficient symbolic
/// certificate. Failure without a witness is unresolved, never a refutation.
///
/// # Errors
/// Empty/ragged matrices, dimensions above 5, or degrees above n=4,x=3.
pub fn check_affine_interlacing_preserver(
    matrix: &[Vec<BivarPolyJson>],
    first: usize,
) -> Result<AffinePreserverReport, RecurrenceError> {
    if matrix.is_empty()
        || matrix.len() > 5
        || matrix[0].is_empty()
        || matrix[0].len() > 5
        || matrix.iter().any(|r| r.len() != matrix[0].len())
    {
        return Err(invalid(
            "preserver matrix must be rectangular with dimensions 1..=5",
        ));
    }
    let matrix = matrix
        .iter()
        .map(|row| row.iter().map(decode_biv).collect::<Result<Vec<_>, _>>())
        .collect::<Result<Vec<_>, _>>()?;
    if matrix
        .iter()
        .flatten()
        .any(|p| p.rows().len() > 5 || p.variable_degree().is_some_and(|d| d > 3))
    {
        return Err(invalid("preserver input degrees exceed n=4,x=3"));
    }
    let shifted: Vec<Vec<_>> = matrix
        .iter()
        .map(|row| {
            row.iter()
                .map(|p| p.shift(&Q::from_integer(first.into()), &Q::zero()))
                .collect()
        })
        .collect();
    let mut report = AffinePreserverReport {
        first_index: first,
        status: "unresolved".into(),
        reason: "requires affine entries and coefficientwise copositive index expansion".into(),
        boundary: vec![],
        pencils: vec![],
        rejections: vec![],
        operator_rejections: vec![],
    };
    let rows = matrix.len();
    let cols = matrix[0].len();
    let entries = shifted.iter().flatten().all(|p| {
        p.variable_degree().is_none_or(|d| d <= 1) && p.rows().iter().all(|r| nonnegative(r))
    });
    let boundary: Vec<Vec<_>> = matrix
        .iter()
        .map(|row| {
            row.iter()
                .map(|p| Polynomial::new(p.evaluate_index(&Q::from_integer(first.into()))))
                .collect()
        })
        .collect();
    let boundary_affine = boundary
        .iter()
        .flatten()
        .all(|p| p.degree().is_none_or(|d| d <= 1) && p.coeffs().iter().all(|c| !c.is_negative()));
    let mut complete = entries;
    for k in 0..rows {
        for l in k..rows {
            for i in 0..cols {
                for j in i..cols {
                    if let Some(check) = boundary_affine
                        .then(|| {
                            affine_pencils::pencil(
                                &boundary[k][i],
                                &boundary[k][j],
                                &boundary[l][i],
                                &boundary[l][j],
                                [k, l],
                                [i, j],
                            )
                        })
                        .flatten()
                    {
                        report.boundary.push(check);
                    } else {
                        complete = false;
                    }
                    if entries {
                        if let Some(check) = symbolic(
                            &shifted[k][i],
                            &shifted[k][j],
                            &shifted[l][i],
                            &shifted[l][j],
                            [k, l],
                            [i, j],
                        ) {
                            report.pencils.push(check);
                        } else {
                            complete = false;
                        }
                    }
                }
            }
        }
    }
    if complete {
        report.status = "proved".into();
        report.reason="Every weakly ordered pencil is certified at the boundary and uniformly on the open index ray by exact copositivity".into();
        return Ok(report);
    }
    // A bounded screen finds rejection certificates, never supplies acceptance.
    let values = [
        Q::new(1.into(), 8.into()),
        Q::one(),
        Q::from_integer(8.into()),
    ];
    for k in 0..rows {
        for l in k..rows {
            for i in 0..cols {
                for j in i..cols {
                    'screen: for index in first..=first.saturating_add(3) {
                        let n = Q::from_integer(index.into());
                        let (a, b, c, d) =
                            (&matrix[k][i], &matrix[k][j], &matrix[l][i], &matrix[l][j]);
                        for lambda in &values {
                            for mu in &values {
                                let ll = vec![mu.clone(), lambda.clone()];
                                let f =
                                    add(&mul(&ll, &b.evaluate_index(&n)), &d.evaluate_index(&n));
                                let g =
                                    add(&mul(&ll, &a.evaluate_index(&n)), &c.evaluate_index(&n));
                                if !directed(&f, &g) {
                                    report.rejections.push(PencilRejection {
                                        rows: [k, l],
                                        columns: [i, j],
                                        index,
                                        parameters: [lambda.to_string(), mu.to_string()],
                                        left: strings(&f),
                                        right: strings(&g),
                                    });
                                    break 'screen;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if !report.rejections.is_empty() {
        report.status = "rejected".into();
        report.reason =
            "An exact positive-parameter pencil violates directed interlacing or real-rootedness"
                .into();
    }
    Ok(report)
}

/// Exact counterexample to a proposed operator-matrix preserver.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OperatorRejection {
    /// Source index.
    pub index: usize,
    /// Ordered PF input, ascending rational-string coefficients.
    pub input: Vec<Vec<String>>,
    /// Output that fails the ordered PF invariant.
    pub output: Vec<Vec<String>>,
}

fn ordered(state: &[Vec<Q>]) -> bool {
    state.iter().all(|p| {
        p.iter().all(|c| !c.is_negative())
            && crate::is_real_rooted_bigint_coeffs(&vp::clear_denominators(p))
    }) && (0..state.len()).all(|i| (i + 1..state.len()).all(|j| directed(&state[i], &state[j])))
}
fn index_nonnegative(p: &Biv, first: usize) -> bool {
    p.variable_degree().is_none_or(|d| d == 0)
        && crate::decide_index_rational_function_sign(
            &p.variable_coefficient(0),
            &[Q::one()],
            &crate::IntegerRange::starting_at(first),
            crate::IndexSignCondition::Nonnegative,
        )
        .holds()
}
/// Recognize [[0,I],[c(n)X,a(n)I+b(n)XD]], a,b,c>=0.
/// For ordered PF (f,g), g interlaces both theta(g) and Xf; positive
/// combinations, with a*g included, remain PF and have g as interleaver.
/// Zero g leaves (0,cXf). This is the Euler-lag argument, independently
/// of the polynomial-entry Brändén criterion.
fn euler_lag(matrix: &OperatorMatrix, first: usize) -> bool {
    if matrix.len() != 2 || matrix[0][0] != Op::default() || matrix[0][1] != Op::identity() {
        return false;
    }
    let zero = Biv::default();
    let a = matrix[1][1].0.get(&(0, 0)).unwrap_or(&zero);
    let bx = matrix[1][1].0.get(&(0, 1)).unwrap_or(&zero);
    let cx = matrix[1][0].0.get(&(0, 0)).unwrap_or(&zero);
    let b = Biv::new(
        bx.variable_coefficient(1)
            .iter()
            .map(|c| vec![c.clone()])
            .collect(),
    );
    let c = Biv::new(
        cx.variable_coefficient(1)
            .iter()
            .map(|c| vec![c.clone()])
            .collect(),
    );
    let x = Biv::new(vec![vec![Q::zero(), Q::one()]]);
    let expected = vec![
        vec![Op::default(), Op::identity()],
        vec![
            Op::single(0, 0, c.multiply(&x)),
            Op::single(0, 0, a.clone()).add(&Op::single(0, 1, b.multiply(&x))),
        ],
    ];
    matrix == &expected && [a, &b, &c].iter().all(|p| index_nonnegative(p, first))
}
pub(super) fn check_operator_preserver(
    matrix: &OperatorMatrix,
    first: usize,
) -> Result<AffinePreserverReport, RecurrenceError> {
    if matrix
        .iter()
        .flatten()
        .all(|op| op.0.keys().all(|(_, d)| *d == 0))
    {
        let entries = matrix
            .iter()
            .map(|r| {
                r.iter()
                    .map(|op| encode_biv(&op.0.get(&(0, 0)).cloned().unwrap_or_default()))
                    .collect()
            })
            .collect::<Vec<_>>();
        return check_affine_interlacing_preserver(&entries, first);
    }
    let mut report=AffinePreserverReport{first_index:first,status:"unresolved".into(),reason:"differential transition is outside the registered Euler-lag template; polynomial 2x2 tests do not apply to Weyl operators".into(),boundary:vec![],pencils:vec![],rejections:vec![],operator_rejections:vec![]};
    if scalar_euler(matrix, first) {
        report.status = "proved".into();
        report.reason="One-coordinate PF preserver: a+b*theta, or X*(a+b*D), with nonnegative index-dependent scalars. Euler/Rolle and multiplication preserve PF, including zero branches".into();
        return Ok(report);
    }
    if euler_lag(matrix, first) {
        report.status = "proved".into();
        report.reason="Euler/Rolle and common-interleaver proof for [[0,I],[c(n)X,a(n)I+b(n)XD]], all a,b,c>=0 on the integer ray, including zero branches; extends the identity-top Euler-lag argument".into();
        return Ok(report);
    }
    let q = |n: i64| Q::from_integer(n.into());
    let mut samples = vec![vec![vec![q(1)]; matrix.len()]];
    if matrix.len() == 3 {
        samples.push(vec![vec![q(1)], vec![q(3), q(1)], vec![q(1), q(1)]]);
        samples.push(vec![
            vec![q(10), q(5)],
            vec![q(26), q(26)],
            vec![q(76), q(761), q(10)],
        ]);
    }
    if matrix.len() == 4 {
        samples.push(vec![
            vec![q(8), q(4)],
            vec![q(4), q(4)],
            vec![q(0), q(4)],
            vec![q(0), q(1)],
        ]);
    }
    for state in samples {
        if !ordered(&state) {
            continue;
        }
        for index in first..=first.saturating_add(3) {
            let output = apply(matrix, &state, index);
            if !ordered(&output) {
                report.operator_rejections.push(OperatorRejection {
                    index,
                    input: state.iter().map(|p| strings(p)).collect(),
                    output: output.iter().map(|p| strings(p)).collect(),
                });
                break;
            }
        }
    }
    if !report.operator_rejections.is_empty() {
        report.status = "rejected".into();
        report.reason = "exact ordered PF input maps outside the invariant".into();
    }
    Ok(report)
}

fn scalar_euler(matrix: &OperatorMatrix, first: usize) -> bool {
    if matrix.len() != 1 {
        return false;
    }
    let zero = Biv::default();
    let a = matrix[0][0].0.get(&(0, 0)).unwrap_or(&zero);
    let bx = matrix[0][0].0.get(&(0, 1)).unwrap_or(&zero);
    let b = Biv::new(
        bx.variable_coefficient(1)
            .iter()
            .map(|c| vec![c.clone()])
            .collect(),
    );
    let x = Biv::new(vec![vec![Q::zero(), Q::one()]]);
    let constant = Biv::new(
        a.variable_coefficient(0)
            .iter()
            .map(|c| vec![c.clone()])
            .collect(),
    );
    let slope = Biv::new(
        a.variable_coefficient(1)
            .iter()
            .map(|c| vec![c.clone()])
            .collect(),
    );
    let expected = Op::single(0, 0, a.clone()).add(&Op::single(0, 1, b.multiply(&x)));
    matrix[0][0] == expected
        && index_nonnegative(&b, first)
        && ((a == &constant && index_nonnegative(a, first))
            || (a == &slope.multiply(&x) && index_nonnegative(&slope, first)))
}
