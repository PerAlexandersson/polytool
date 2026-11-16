//! Exact terminating three-over-two recognition after multiplier matching fails.
//!
//! Parameters follow arXiv:2208.03539v2, (8.1), (8.9) and Theorem 8.4.
//! The weak boundary uses the coefficientwise closure proved in
//! `RealRooted.Hypergeometric.threeFTwo_splits` and `threeFTwo_interlaces`.
//! The existing all-index row-annihilator and boundary-coefficient proofs
//! supply the coefficient identity; parity tests here compare affine
//! coefficients, rather than sampling a finite collection of rows.

use super::*;

const REFERENCE: &str = "Terminating 3F2: arXiv:2208.03539v2, (8.1), (8.9), Theorem 8.4; \
    r=s=2, lambda_1=1, lambda_2=2. Positive a1,a2,b1,b2 with a1<=b1, a2-1<=b1, \
    a1<=b2, a2<=b2 give real positive zeros and weak consecutive interlacing, including \
    the weak boundary by coefficientwise closure (RealRooted.Hypergeometric.threeFTwo_splits, \
    threeFTwo_interlaces). Strict inequalities give simple zeros and strict interlacing \
    (RealRooted.Hypergeometric.threeFTwo_interlaces_and_isCoprime).";

/// Exact three-over-two identification within the row-hypergeometric route.
///
/// For `M(n) = termination[0] + termination[1] n`, the normalized row is
/// `U_0(n) threeFTwo(a1,a2,b1,b2,M(n))(scale*x)`. The original row is
/// `x^valuation(n)` times this row, or its reversal at `degree(n)`.
/// Every rational and integer is encoded as an exact string. Replay through
/// [`RealRootednessTheoremReport::verify`] recomputes the ratio identity,
/// parity assignments, all-index annihilator, support and nonzero anchor.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct HypergeometricThreeFTwoCertificate {
    /// Positive lower parameters `[a1, a2]`, in the theorem's order.
    pub lower_parameters: [String; 2],
    /// Positive upper base parameters `[b1, b2]`, before the degree shifts.
    pub upper_parameters: [String; 2],
    /// Integer affine terminating degree `[constant, slope]`.
    pub termination: [String; 2],
    /// Nonzero, index-independent rational argument scale.
    pub scale: String,
    /// Whether the normalized row reads coefficients in reverse degree order.
    pub reversed: bool,
    /// Integer affine valuation `[constant, slope]` of the original row.
    pub valuation: [String; 2],
    /// Integer affine degree `[constant, slope]` of the original row.
    pub degree: [String; 2],
    /// First absolute row index where the ratio and support are certified.
    pub structure_start: usize,
    /// Whether all four comparison inequalities in (8.1) are strict.
    /// Monomial factors can still prevent strict interlacing of original rows.
    pub strict_parameters: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ThreeFTwoMatching {
    lower: [Q; 2],
    upper: [Q; 2],
    termination: Affine,
    scale: Q,
    strict: bool,
}

fn parameter_range(lower: &[Q; 2], upper: &[Q; 2]) -> Option<bool> {
    if lower.iter().chain(upper).any(|x| !x.is_positive()) {
        return None;
    }
    let gaps = [
        &upper[0] - &lower[0],
        &upper[0] - &lower[1] + Q::one(),
        &upper[1] - &lower[0],
        &upper[1] - &lower[1],
    ];
    (!gaps.iter().any(Signed::is_negative)).then(|| gaps.iter().all(Signed::is_positive))
}

/// Compare both parity restrictions symbolically. In `n=2t+r`, both actual
/// upper parameters and the expected ones have coefficient `M_1` at `t`.
/// Thus comparing their unordered intercept pairs proves their equality.
fn parity_matches(upper: &[Affine; 2], base: &[Q; 2], size: &Affine) -> bool {
    let two = Q::from_integer(2.into());
    upper.iter().all(|a| &a.1 * &two == size.1)
        && (0..2).all(|r| {
            let r = Q::from_integer(r.into());
            let degree = &size.0 + &size.1 * &r;
            let mut actual = upper.clone().map(|a| a.0 + a.1 * &r);
            let mut expected = [
                &base[0] + (&degree / &two).floor(),
                &base[1] + ((&degree - Q::one()) / &two).floor(),
            ];
            actual.sort();
            expected.sort();
            actual == expected
        })
}

fn find_matching(ratio: &Ratio, first: usize) -> Option<ThreeFTwoMatching> {
    if ratio.numerators.len() != 3 || ratio.denominators.len() != 3 {
        return None;
    }
    let scale = ratio.scale.numerator.last()? / ratio.scale.denominator.last()?;
    if scale.is_zero()
        || vec_poly::trim(&ratio.scale.numerator)
            != vec_poly::trim(&vec_poly::scale(&ratio.scale.denominator, &scale))
    {
        return None;
    }
    let mut lower = ratio.denominators.clone();
    let factorial = lower.iter().position(|b| *b == affine_constant(0))?;
    lower.remove(factorial);
    if lower.iter().any(|a| !a.1.is_zero()) {
        return None;
    }
    let lower = [&lower[0].0 + Q::one(), &lower[1].0 + Q::one()];
    let mut weak = None;
    for (index, factor) in ratio.numerators.iter().enumerate() {
        let Some(size) = complement_size(factor) else {
            continue;
        };
        // A positive terminating degree avoids the unused degree-zero
        // upperTwo convention; earlier rows are checked separately.
        if !matches!(
            affine_function(&size).decide(first, IndexSignCondition::Positive),
            IndexSignDecision::Holds { .. }
        ) {
            continue;
        }
        let upper: Vec<Affine> = ratio
            .numerators
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != index)
            .map(|(_, a)| affine_add(a, &affine_constant(1)))
            .collect();
        let two = Q::from_integer(2.into());
        for swap in [false, true] {
            let first = usize::from(swap);
            let base = [
                &upper[first].0 - (&size.0 / &two).floor(),
                &upper[1 - first].0 - ((&size.0 - Q::one()) / &two).floor(),
            ];
            if !parity_matches(&[upper[0].clone(), upper[1].clone()], &base, &size) {
                continue;
            }
            for lower in [lower.clone(), [lower[1].clone(), lower[0].clone()]] {
                if let Some(strict) = parameter_range(&lower, &base) {
                    let matching = ThreeFTwoMatching {
                        lower,
                        upper: base.clone(),
                        termination: size.clone(),
                        scale: scale.clone(),
                        strict,
                    };
                    if strict {
                        return Some(matching);
                    }
                    if weak.is_none() {
                        weak = Some(matching);
                    }
                }
            }
        }
    }
    weak
}

pub(super) struct RowContext<'a> {
    pub shape: &'a RecurrenceShape,
    pub p: &'a LinearFactors,
    pub q: &'a LinearFactors,
    pub valuation: &'a Affine,
    pub degree: Option<&'a Affine>,
    pub rows: &'a [Vec<Q>],
    pub first_index: usize,
    pub rule_start: usize,
    pub from: usize,
    pub shape_from: usize,
}

pub(super) fn try_report(
    mut report: RealRootednessTheoremReport,
    context: RowContext<'_>,
    options: &RealRootednessTheoremOptions,
) -> Option<RealRootednessTheoremReport> {
    let RowContext {
        shape,
        p,
        q,
        valuation,
        degree,
        rows,
        first_index,
        rule_start,
        from,
        shape_from,
    } = context;
    let degree = degree?;
    let mut candidate = None;
    'search: for coordinates in [Coordinates::Identity, Coordinates::Reversal] {
        let ratio = ratio_in_coordinates(p, q, valuation, Some(degree), coordinates)?;
        for s in shape_from.max(from)..=from + MAX_STRUCTURE_DELAY {
            let Some(matching) = find_matching(&ratio, s) else {
                continue;
            };
            // Full support keeps both endpoints nonzero and prevents extra
            // roots at zero from degree padding in the reversal.
            if affine_add(&matching.termination, valuation) != *degree {
                continue;
            }
            let support = support_hypotheses(
                shape,
                p,
                q,
                valuation,
                (coordinates == Coordinates::Reversal).then_some(degree),
                (rows, first_index, rule_start),
                s,
            );
            if support
                .iter()
                .all(|h| h.status == LiftHypothesisStatus::Verified)
            {
                candidate = Some((coordinates, ratio, matching, s, support));
                break 'search;
            }
        }
    }
    let (coordinates, ratio, matching, s, support) = candidate?;
    report.reference = REFERENCE.into();
    report.lean_declarations = [
        "threeFTwo_splits",
        "threeFTwo_interlaces",
        "threeFTwo_interlaces_and_isCoprime",
        "isRoot_threeFTwo_pos",
    ]
    .map(|name| format!("RealRooted.Hypergeometric.{name}"))
    .to_vec();
    report.hypotheses.extend(support);
    report.hypotheses.push(sign_hypothesis(
        "three_f_two_relation_scalar",
        "the row relation scalar does not vanish",
        "c_p(n)",
        &IndexFunction::polynomial(p.scalar.clone()),
        s,
        IndexSignCondition::Nonzero,
    ));
    report.hypotheses.push(sign_hypothesis(
        "three_f_two_valuation",
        "the original monomial exponent is nonnegative",
        "v(n)",
        &affine_function(valuation),
        s,
        IndexSignCondition::Nonnegative,
    ));
    report.hypotheses.push(exact_hypothesis(
        "three_f_two_parameters",
        "positive parameters in the weak range (8.1)",
        Ok(format!(
            "a1={}, a2={}, b1={}, b2={}; all four comparisons {}; \
            the unordered upper pair agrees symbolically in both parity classes",
            matching.lower[0],
            matching.lower[1],
            matching.upper[0],
            matching.upper[1],
            if matching.strict { "strict" } else { "weak" }
        )),
    ));
    report.hypotheses.push(exact_hypothesis(
        "three_f_two_ratio",
        "the normalized coefficient ratio equals the terminating 3F2 ratio",
        Ok(format!(
            "R(n,j)={}; M(n)={}; argument scale={}; \
            the denominator factors j, j+a1-1, j+a2-1 are positive for j>=1",
            ratio.text(),
            affine_text(&matching.termination),
            matching.scale
        )),
    ));
    let anchor = if coordinates == Coordinates::Identity {
        valuation
    } else {
        degree
    };
    let (mut anchor_hypothesis, start, _) = row_first_coefficient_hypothesis(
        shape,
        &ratio,
        coordinates,
        anchor,
        (rows, first_index, rule_start),
        s,
        false,
    );
    if anchor_hypothesis.status != LiftHypothesisStatus::Verified {
        let (alternating, _, _) = row_first_coefficient_hypothesis(
            shape,
            &ratio,
            coordinates,
            anchor,
            (rows, first_index, rule_start),
            s,
            true,
        );
        if alternating.status == LiftHypothesisStatus::Verified {
            anchor_hypothesis = alternating;
        }
    }
    report.hypotheses.push(anchor_hypothesis);
    let closed = *valuation != affine_constant(0);
    let positive = matching.scale.is_positive();
    let region = RealInterval {
        lower: positive.then(Q::zero),
        lower_closed: positive && closed,
        upper: (!positive).then(Q::zero),
        upper_closed: !positive && closed,
    };
    let failing =
        (first_index..start).find(|n| !has_real_roots_in_region(&rows[n - first_index], &region));
    report.hypotheses.push(exact_hypothesis(
        "earlier_rows",
        "earlier rows are nonzero with every root in the certified region",
        match failing {
            None => Ok("checked exactly".into()),
            Some(n) => Err((Some(n), format!("P_{n} fails the exact check"))),
        },
    ));
    if let Some(requested) = &options.root_region {
        report.hypotheses.push(exact_hypothesis(
            "requested_region",
            "the requested open region contains the certified region",
            if requested_region_contains(requested, &region) {
                Ok("exact containment".into())
            } else {
                Err((None, "the certified region is not contained".into()))
            },
        ));
    }
    let interlacing = (matching.termination.1 == Q::one() && valuation.1.is_zero()).then_some(
        if matching.strict && !closed {
            RowInterlacing::Strict
        } else {
            RowInterlacing::Weak
        },
    );
    let simple = matching.strict
        && valuation.1.is_zero()
        && valuation.0 <= Q::one()
        && (first_index..start).all(|n| has_simple_real_roots(&rows[n - first_index]));
    let encode = |a: &Affine| [a.0.to_string(), a.1.to_string()];
    report.certificate = Some(RouteCertificate::HypergeometricRow(
        HypergeometricThreeFTwoCertificate {
            lower_parameters: matching.lower.clone().map(|q| q.to_string()),
            upper_parameters: matching.upper.clone().map(|q| q.to_string()),
            termination: encode(&matching.termination),
            scale: matching.scale.to_string(),
            reversed: coordinates == Coordinates::Reversal,
            valuation: encode(valuation),
            degree: encode(degree),
            structure_start: s,
            strict_parameters: matching.strict,
        },
    ));
    report.normal_form.push_str(&format!(
        "; terminating 3F2 in {} coordinates, M(n)={}, \
        scale={}, a1={}, a2={}, b1={}, b2={}, from n={s}",
        if coordinates == Coordinates::Reversal {
            "reversal"
        } else {
            "identity"
        },
        affine_text(&matching.termination),
        matching.scale,
        matching.lower[0],
        matching.lower[1],
        matching.upper[0],
        matching.upper[1]
    ));
    let interlacing_from =
        interlacing.map(|kind| extend_interlacing_backward(rows, first_index, s + 1, kind));
    let interlacing_text = match (interlacing, interlacing_from) {
        (Some(RowInterlacing::Strict), Some(n)) => {
            format!("strict consecutive interlacing from index {n}")
        }
        (Some(RowInterlacing::Weak), Some(n)) => {
            format!("weak consecutive interlacing from index {n}")
        }
        _ => "consecutive interlacing is not claimed".into(),
    };
    let identity = if coordinates == Coordinates::Reversal {
        "P_n(x)=U_0(n) x^d(n) threeFTwo(a1,a2,b1,b2,M(n))(scale/x)"
    } else {
        "P_n(x)=U_0(n) x^v(n) threeFTwo(a1,a2,b1,b2,M(n))(scale*x)"
    };
    Some(report.finish(RealRootednessTheoremConclusion {
        first_index,
        zero_residue_classes: None,
        interlacing_from,
        interlacing,
        simple_roots: simple,
        root_region: Some(RootInterval::from_real(&region)),
        statement: format!(
            "every P_n, n >= {first_index}, is nonzero and real-rooted with roots \
            in {}; {identity} for n >= {s}; nonzero roots have the sign of scale={}; \
            {interlacing_text}; simple roots={simple}",
            RootInterval::from_real(&region),
            matching.scale
        ),
    }))
}

pub(in crate::recurrence) fn verify_certificate_in_report(
    certificate: &HypergeometricThreeFTwoCertificate,
    _report: &RealRootednessTheoremReport,
) -> Result<(), RecurrenceError> {
    // Full report replay has already recomputed the row identity and compared
    // every certificate field. Also validate the exact parameter payload.
    let lower = certificate
        .lower_parameters
        .clone()
        .map(|x| parse_rational_coeff(&x));
    let upper = certificate
        .upper_parameters
        .clone()
        .map(|x| parse_rational_coeff(&x));
    let [a1, a2] = lower;
    let [b1, b2] = upper;
    let valid = parameter_range(&[a1?, a2?], &[b1?, b2?]);
    if valid != Some(certificate.strict_parameters)
        || parse_rational_coeff(&certificate.scale)?.is_zero()
    {
        return Err(RecurrenceError::VerificationFailed(
            "invalid terminating 3F2 parameters".into(),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
