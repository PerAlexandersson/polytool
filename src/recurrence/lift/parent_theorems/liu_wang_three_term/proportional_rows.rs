//! The degenerate case of the three-term theorem where every row is a
//! scalar multiple of a power of one polynomial times the first seed.
//!
//! When `A_n(x) = a_n g(x)` and `B_n(x) = b_n g(x)^2` for one polynomial
//! `g` (possibly constant) and scalars `a_n`, `b_n` (rational functions of
//! `n`), and the seeds of the induction start satisfy
//! `P_(N-1) = s g P_(N-2)` with a nonzero constant `s`, then by induction
//!
//! ```text
//! P_n = c_n g^(n - N + 2) P_(N-2),   c_n = a_n c_(n-1) + b_n c_(n-2),
//! c_(N-2) = 1, c_(N-1) = s,
//! ```
//!
//! for every `n >= N - 2`.  The rows are real-rooted exactly when `g` and
//! `P_(N-2)` are and `c_n != 0`; Liu and Wang's sign condition holds
//! trivially, since `B_n` vanishes at every root of `g`, and the theorem
//! reduces to the scalar sequence.  The check proves `c_n != 0` for every
//! `n >= N` by one of:
//!
//! - a positive sequence (after the sign normalization `(-1)^n` when
//!   `s < 0`): nonnegative step coefficients or a positive-ratio induction,
//!   as in the degree law of the three-term check;
//! - a hypergeometric closed form: a rational function `r(n)` with
//!   `r(n) = a_n + b_n / r(n - 1)` identically, `c_N / c_(N-1) = r(N)` and
//!   `r(n)` finite and nonzero for `n >= N`, so that
//!   `c_n = c_(N-1) r(N) ... r(n)` (the scalar coefficient `a_n` may vanish
//!   at isolated indices).
//!
//! This covers rules such as `P_n = x P_(n-1) + x^2 P_(n-2)` (Fibonacci
//! multiples of powers of `x`) or `P_n = (12 - 4n) x P_(n-1) +
//! (4n - 3) x^2 P_(n-2)` (`P_n = n^2 x^(n-1)`), where `B_n >= 0` defeats the
//! general sign condition.  The roots are those of `g` and `P_(N-2)`; the
//! rows `P_(n-1)`, `P_n` interlace weakly when `deg g <= 1` (`P_n` adds at
//! most one root, repeated or not).

use super::*;

/// `p` is a constant multiple of `q`: the constant, or `None`.
fn constant_ratio(p: &[Q], q: &[Q]) -> Option<Q> {
    let (p, q) = (vec_poly::trim(p), vec_poly::trim(q));
    let pivot = q.iter().rposition(|c| !c.is_zero())?;
    let ratio = p.get(pivot).cloned().unwrap_or_else(Q::zero) / &q[pivot];
    (vec_poly::trim(&vec_poly::scale(&q, &ratio)) == p).then_some(ratio)
}

/// The region of the roots of the nonzero real-rooted `p`: `(-inf, 0)`
/// when every root is negative, `(-inf, 0]` when every root is `<= 0`.
fn nonpositive_region(p: &[Q]) -> Option<RealInterval> {
    let p = vec_poly::trim(p);
    let open = has_real_roots_below(&p, &Q::zero());
    let closed = open || {
        let order = p.iter().position(|c| !c.is_zero()).unwrap_or(0);
        has_real_roots_below(&p[order..], &Q::zero()) || p.len() == order + 1
    };
    closed.then(|| RealInterval {
        lower: None,
        lower_closed: false,
        upper: Some(Q::zero()),
        upper_closed: !open,
    })
}

/// The proportional-rows check (module documentation), or `None` when the
/// rule or the seeds do not have that form.
pub(super) fn proportional_rows_report(
    mut report: RealRootednessTheoremReport,
    c1: &Bivariate,
    c2: &Bivariate,
    input: &AnalysisInput<'_>,
) -> Option<RealRootednessTheoremReport> {
    let (gamma1, g) = separate(c1)?;
    let (gamma2, w2) = separate(c2)?;
    let g = vec_poly::trim(&g);
    let lambda = constant_ratio(&w2, &vec_poly::mul(&g, &g))?;
    let seeds = input.seeds;
    let previous = vec_poly::trim(&seeds[seeds.len() - 2]);
    let current = vec_poly::trim(&seeds[seeds.len() - 1]);
    if vec_poly::is_zero(&previous) || vec_poly::is_zero(&current) {
        return None;
    }
    let s = constant_ratio(&current, &vec_poly::mul(&g, &previous))?;
    let f = input.f;
    let n0 = input.first_index + seeds.len();
    let a = IndexFunction::new(gamma1, f.to_vec());
    let b = IndexFunction::new(vec_poly::scale(&gamma2, &lambda), f.to_vec());
    let g_text = crate::Polynomial::new(g.clone()).to_string();
    report.normal_form = format!(
        "P_n = a_n g(x) P_(n-1) + b_n g(x)^2 P_(n-2) for n >= {n0}, with g(x) = {g_text}, a_n = \
         {}, b_n = {}, and P_{} = {} g(x) P_{} (proportional rows: P_n = c_n g(x)^(n - {}) \
         P_{}, c_n = a_n c_(n-1) + b_n c_(n-2))",
        a.text(),
        b.text(),
        n0 - 1,
        format_rational_coeff(&s),
        n0 - 2,
        n0 - 2,
        n0 - 2
    );
    report.reference = format!(
        "{REFERENCE}; the degenerate case A_n = a_n g, B_n = b_n g^2, where every row is \
         c_n g^(n - N + 2) P_(N-2) and the sign condition holds since B_n vanishes at the roots \
         of g"
    );
    let mut hypotheses = Vec::new();
    hypotheses.extend(denominator_hypothesis(f, n0));
    // g and P_(N-2) real-rooted (in the requested region).
    let region = input.options.root_region.as_ref();
    let product = vec_poly::mul(&g, &previous);
    let factors_ok = match region {
        Some(region) => has_real_roots_in(&product, region),
        None => is_real_rooted_rational(&product),
    };
    hypotheses.push(exact_hypothesis(
        "factors_real_rooted",
        &format!(
            "g(x) = {g_text} and P_{} are nonzero and real-rooted{}",
            n0 - 2,
            region.map_or(String::new(), |r| format!(
                " with every root in {}",
                open_root_interval(r)
            ))
        ),
        if factors_ok {
            Ok("checked exactly".into())
        } else {
            Err((
                None,
                "g(x) P_(N-2) is not real-rooted (in the region) (exact)".into(),
            ))
        },
    ));
    // c_n != 0 for every n >= N.
    let mut values = vec![Q::one(), s.clone()];
    for n in n0..n0 + GENERATED_ROWS {
        let (Some(a_n), Some(b_n)) = (a.at(n), b.at(n)) else {
            break;
        };
        let k = values.len();
        values.push(a_n * &values[k - 1] + b_n * &values[k - 2]);
    }
    hypotheses.push(scalar_sequence_hypothesis(&a, &b, n0, &values));
    hypotheses.extend(earlier_seeds_hypothesis(
        &seeds[..seeds.len() - 2],
        input.first_index,
        region,
    ));
    let interlacing = vec_poly::degree(&g).unwrap_or(0) <= 1;
    let root_region = match region {
        Some(region) => Some(open_root_interval(region)),
        None => nonpositive_region(&product).map(|r| RootInterval::from_real(&r)),
    };
    let statement = format!(
        "every P_n, n >= {}, is nonzero and real-rooted, P_n = c_n g(x)^(n - {}) P_{} with the \
         roots of g(x) = {g_text} and of P_{}{}{}",
        input.first_index,
        n0 - 2,
        n0 - 2,
        n0 - 2,
        if interlacing {
            format!("; P_(n-1) interlaces P_n weakly for n >= {}", n0 - 1)
        } else {
            String::new()
        },
        root_region
            .as_ref()
            .map_or(String::new(), |r| format!("; every root lies in {r}"))
    );
    report.hypotheses = hypotheses;
    let conclusion = RealRootednessTheoremConclusion {
        first_index: input.first_index,
        zero_residue_classes: None,
        interlacing_from: interlacing.then_some(n0 - 1),
        interlacing: interlacing.then_some(RowInterlacing::Weak),
        simple_roots: false,
        root_region,
        statement,
    };
    Some(report.finish(conclusion))
}

/// `c_n != 0` for every `n >= first - 2` (module documentation): `values`
/// are `c_(first-2) = 1, c_(first-1), ...` computed exactly.
fn scalar_sequence_hypothesis(
    a: &IndexFunction,
    b: &IndexFunction,
    first: usize,
    values: &[Q],
) -> LiftHypothesis {
    let name = "scalar_sequence_nonzero";
    let statement = format!(
        "c_n != 0 for every n >= {}, where c_{} = 1, c_{} = s and c_n = a_n c_(n-1) + b_n \
         c_(n-2)",
        first - 2,
        first - 2,
        first - 1
    );
    if let Some(offset) = values.iter().position(Zero::is_zero) {
        return exact_hypothesis(
            name,
            &statement,
            Err((
                Some(first - 2 + offset),
                format!(
                    "c_{} = 0, so P_{} = 0 (exact)",
                    first - 2 + offset,
                    first - 2 + offset
                ),
            )),
        );
    }
    // A positive sequence, after (-1)^n when c_(first-1) < 0.
    let alternating = values[1].is_negative();
    let (s, normalized): (IndexFunction, Vec<Q>) = if alternating {
        (
            a.negated(),
            values
                .iter()
                .enumerate()
                .map(|(k, v)| if k % 2 == 1 { -v } else { v.clone() })
                .collect(),
        )
    } else {
        (a.clone(), values.to_vec())
    };
    let positive = positive_sequence_hypothesis(
        name,
        &statement,
        if alternating { "(-1)^n c_n" } else { "c_n" },
        &s,
        b,
        first,
        &normalized,
        true,
        false,
    );
    if positive.status == LiftHypothesisStatus::Verified {
        return positive;
    }
    // A hypergeometric closed form.
    let ratios: Vec<(usize, Q)> = values
        .windows(2)
        .enumerate()
        .skip(1)
        .map(|(offset, pair)| (first - 1 + offset, &pair[1] / &pair[0]))
        .collect();
    if let Some(r) = fit_index_rational_function(&ratios) {
        let identity = a.add(&b.div(&r.shifted(-1))).sub(&r);
        let first_ratio = ratios.first().map(|(_, ratio)| ratio.clone());
        if identity.numerator.iter().all(Zero::is_zero)
            && r.at(first).is_some()
            && r.at(first) == first_ratio
        {
            let nonzero = r.decide(first, IndexSignCondition::Nonzero);
            if nonzero.holds() {
                let mut hypothesis =
                    LiftHypothesis::from_sign_decision(name, &statement, nonzero, "r(n)");
                hypothesis.evidence = format!(
                    "hypergeometric closed form c_n = c_{} r({first}) ... r(n) with r(n) = {}: \
                     r(n) = a_n + b_n / r(n - 1) identically (exact), c_{first} / c_{} = \
                     r({first}) (exact), and {}",
                    first - 1,
                    r.text(),
                    first - 1,
                    hypothesis.evidence
                );
                return hypothesis;
            }
        }
    }
    positive
}
