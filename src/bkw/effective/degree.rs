//! Newton weights and exact rational coefficient layers for degree closure.

use super::*;
use crate::integer_range_sign::{
    decide_index_polynomial_sign, IndexSignCertificate, IndexSignCondition, IndexSignDecision,
    IntegerRange,
};
use crate::recurrence::RationalGeneratingSeries;

fn gcd(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

pub(super) fn newton_weight(d: &B) -> Option<(usize, usize)> {
    let (mut growth, mut weight) = (0, 1);
    for (j, row) in d.rows().iter().enumerate().skip(1) {
        if let Some(k) = P::new(row.clone()).degree() {
            if k * weight > growth * j {
                (growth, weight) = (k, j);
            }
        }
    }
    if growth == 0 {
        return None;
    }
    let common = gcd(growth, weight);
    Some((growth / common, weight / common))
}

/// A binomial face with k>4 equally sized modes prevents global pair
/// domination on a negative-infinity ray for a constant numerator. At
/// n=k*m, the pair has limiting weight 2/k and the remainder (k-2)/k.
pub(super) fn equimodular_obstruction(n: &B, d: &B) -> Option<usize> {
    if n.index_degree() != Some(0) || n.rows()[0].len() != 1 {
        return None;
    }
    let (growth, weight) = newton_weight(d)?;
    let face: Vec<_> = d
        .rows()
        .iter()
        .enumerate()
        .skip(1)
        .filter_map(|(j, row)| {
            let degree = P::new(row.clone()).degree()?;
            (weight * degree == growth * j).then_some(j)
        })
        .collect();
    match face.as_slice() {
        [k] if *k > 4 => Some(*k),
        _ => None,
    }
}

/// The same equal-mode obstruction at t->0 when the unique minimum
/// valuation/lag edge has more than four roots and N is constant.
pub(super) fn zero_edge_obstruction(n: &B, d: &B) -> Option<usize> {
    if n.index_degree() != Some(0) || n.rows()[0].len() != 1 {
        return None;
    }
    let terms: Vec<_> = d
        .rows()
        .iter()
        .enumerate()
        .skip(1)
        .filter_map(|(j, row)| row.iter().position(|c| !c.is_zero()).map(|v| (j, v)))
        .collect();
    let &(lag, valuation) = terms
        .iter()
        .min_by(|(j, v), (k, w)| (v * k).cmp(&(w * j)))?;
    if valuation == 0 || lag <= 4 {
        return None;
    }
    (terms
        .iter()
        .filter(|(j, v)| v * lag == valuation * j)
        .count()
        == 1)
        .then_some(lag)
}

/// The top Newton edge gives a cap, not an equality. Expand at that edge
/// until each residue has an identically zero prefix and a provably positive
/// layer. Negative coefficients or exhausted guards give no certificate.
pub(super) fn newton_budget(n: &B, d: &B) -> Result<DegreeBudget, RecurrenceError> {
    if d.rows().first() != Some(&vec![Q::one()]) {
        return Err(failure("Newton degree expansion requires D(t,0)=1"));
    }
    let (growth, weight) = newton_weight(d)
        .ok_or_else(|| failure("Newton degree slope is zero; no growing phase grid"))?;
    let shift = n
        .rows()
        .iter()
        .enumerate()
        .map(|(j, row)| {
            (weight * P::new(row.clone()).degree().unwrap_or(0)).saturating_sub(growth * j)
        })
        .max()
        .unwrap_or(0);
    let layer = |a: &B, extra: usize, e: usize| -> P {
        let mut out = vec![Q::zero(); a.rows().len()];
        for (j, row) in a.rows().iter().enumerate() {
            for (k, coefficient) in row.iter().enumerate() {
                if extra + growth * j == e + weight * k {
                    out[j] += coefficient;
                }
            }
        }
        P::new(out)
    };
    let face = layer(d, 0, 0);
    if face.coeff(0) != Q::one() {
        return Err(failure(
            "Newton face reciprocal lacks a nonnegative-series certificate",
        ));
    }
    let support: Vec<_> = face
        .coeffs()
        .iter()
        .enumerate()
        .skip(1)
        .filter(|(_, c)| !c.is_zero())
        .map(|(j, _)| j)
        .collect();
    let Some(&step) = support.first() else {
        return Err(failure(
            "Newton face has no positive-degree denominator term",
        ));
    };
    // Binomial powers have positive reciprocals even when their expanded
    // denominators have positive nonconstant coefficients.
    let multiplicity = face.degree().unwrap() / step;
    let mut binomial = vec![Q::zero(); step + 1];
    binomial[0] = Q::one();
    binomial[step] = face.coeff(step) / Q::from_integer(multiplicity.into());
    let binomial = P::new(binomial);
    let binomial_power = face.degree().unwrap() % step == 0
        && binomial.coeff(step).is_negative()
        && binomial.pow(multiplicity) == face;
    if !binomial_power && face.coeffs().iter().skip(1).any(Q::is_positive) {
        return Err(failure(
            "Newton face reciprocal lacks a nonnegative-series certificate",
        ));
    }
    let period = support.iter().copied().fold(0, gcd);
    // A residue class of rows that is identically zero never reaches a
    // positive Newton layer, and the layer search below would then run to its
    // guard with rapidly growing exact arithmetic. Rejection needs no proof, so
    // an exact expansion of an initial range of rows suffices here.
    let checked_rows = 8 * period + n.rows().len() + d.rows().len();
    let mut rows: Vec<P> = Vec::with_capacity(checked_rows);
    for j in 0..checked_rows {
        let mut value = n
            .rows()
            .get(j)
            .map_or_else(|| P::new(Vec::new()), |r| P::new(r.clone()));
        for (i, row) in d
            .rows()
            .iter()
            .enumerate()
            .skip(1)
            .take_while(|(i, _)| *i <= j)
        {
            value = value - P::new(row.clone()) * rows[j - i].clone();
        }
        rows.push(value);
    }
    for r in 0..period {
        if (n.rows().len()..checked_rows)
            .filter(|j| j % period == r)
            .all(|j| rows[j].is_zero())
        {
            return Err(failure(format!(
                "rows in residue class {r} modulo {period} vanish on the checked range; \
                 no positive Newton layer"
            )));
        }
    }
    let mut first = vec![None; period];
    let mut start = 0;
    let mut numerators: Vec<P> = Vec::new();
    let mut layers = Vec::new();
    let mut reachability_seeds = vec![Vec::new(); period];
    let mut leading_sign_certificates = Vec::new();
    // N/(1-B) with coefficientwise N,B>=0 has positive layer expansions:
    // G_e=(N_e+sum_(i>0) B_i G_(e-i))/face. This proves nonnegativity
    // even when combining the fractions introduces signed numerators.
    let positive_cone = n.rows().iter().flatten().all(|c| !c.is_negative())
        && d.rows().iter().enumerate().all(|(j, row)| {
            row.iter()
                .enumerate()
                .all(|(k, c)| (j == 0 && k == 0) || !c.is_positive())
        });
    let has_lower_denominator_terms = d.rows().iter().enumerate().any(|(j, row)| {
        row.iter()
            .enumerate()
            .any(|(k, c)| !c.is_zero() && growth * j > weight * k)
    });
    let maximum_layer = if has_lower_denominator_terms {
        64
    } else {
        n.rows()
            .iter()
            .enumerate()
            .flat_map(|(j, row)| {
                row.iter()
                    .enumerate()
                    .filter(|(_, c)| !c.is_zero())
                    .map(move |(k, _)| shift + growth * j - weight * k)
            })
            .max()
            .unwrap_or(0)
            .min(64)
    };
    for e in 0..=maximum_layer {
        // [u^e] N*/D* = A_e / face^(e+1).
        let mut a = layer(n, shift, e) * face.pow(e);
        for i in 1..=e {
            a = a - layer(d, 0, i) * numerators[e - i].clone() * face.pow(i - 1);
        }
        let denominator = face.pow(e + 1);
        layers.push(RationalGeneratingSeries::new(
            a.coeffs(),
            denominator.coeffs(),
        )?);
        let mut positive_numerator = a.clone();
        let mut positive_power = e + 1;
        while positive_power > 1 && !positive_numerator.is_zero() {
            let (quotient, remainder) = positive_numerator.div_rem(&face);
            if !remainder.is_zero() {
                break;
            }
            positive_numerator = quotient;
            positive_power -= 1;
        }
        let positive_denominator = face.pow(positive_power);
        for (r, chosen) in first.iter_mut().enumerate() {
            if chosen.is_some() {
                continue;
            }
            let projection = P::new(
                positive_numerator
                    .coeffs()
                    .iter()
                    .enumerate()
                    .map(|(j, c)| {
                        if j % period == r {
                            c.clone()
                        } else {
                            Q::zero()
                        }
                    })
                    .collect(),
            );
            if projection.is_zero() {
                continue;
            }
            // A nonzero rational layer can still have finite support.
            // Remove its finite exception range before seeking the first
            // eventual leading layer; no sampled cancellation is assumed.
            let (finite, remainder) = projection.div_rem(&positive_denominator);
            if remainder.is_zero() {
                start = start.max(finite.degree().unwrap_or(0) + 1);
                continue;
            }
            if !positive_cone && projection.coeffs().iter().any(Q::is_negative) {
                let base = if binomial_power { &binomial } else { &face };
                let exponent = positive_power * if binomial_power { multiplicity } else { 1 };
                let (seed, certificate) = binomial_positivity(&projection, base, exponent, r, step)
                    .map_err(|_| {
                        failure(format!(
                            "Newton coefficient layer {e}, residue {r}: positivity not certified"
                        ))
                    })?;
                start = start.max(seed);
                reachability_seeds[r].push(seed);
                leading_sign_certificates.push(certificate);
                *chosen = Some(e);
                continue;
            }
            // The reciprocal of face^(e+1) is nonnegative and has a positive
            // coefficient at step. Once a residue modulo step is reached,
            // adding that step preserves strict positivity. Check one seed
            // in every such residue, then use that implication for all n.
            let coefficients = series_coefficients(&projection, &positive_denominator, 1024);
            if coefficients.iter().any(Q::is_negative) {
                return Err(failure(
                    "Newton layer contradicts its nonnegative-series certificate",
                ));
            }
            let mut reached = vec![None; step];
            for (j, c) in coefficients.iter().enumerate() {
                if c.is_positive() && reached[j % step].is_none() {
                    reached[j % step] = Some(j);
                }
            }
            for residue in (r..step).step_by(period) {
                let seed = reached[residue].ok_or_else(|| {
                    failure(
                        "Newton positive-series reachability exceeded the 1024-coefficient guard",
                    )
                })?;
                start = start.max(seed);
                reachability_seeds[r].push(seed);
            }
            *chosen = Some(e);
        }
        numerators.push(a);
        if first.iter().all(Option::is_some) {
            break;
        }
    }
    let first: Vec<usize> = first
        .into_iter()
        .collect::<Option<_>>()
        .ok_or_else(|| failure("Newton degree layers exceeded the 64-layer guard"))?;
    let offsets: Vec<i64> = first.iter().map(|e| shift as i64 - *e as i64).collect();
    let count = EffectiveBkwZeroCount {
        slope: growth,
        offsets,
        divisor: weight,
    };
    for (r, offset) in count.offsets.iter().enumerate() {
        if (growth as i64 * r as i64 + offset) % weight as i64 != 0 {
            return Err(failure("Newton layer has an inconsistent index congruence"));
        }
    }
    let last_layer = *first.iter().max().unwrap();
    let leading_denominator = face.pow(last_layer + 1);
    let mut leading_numerator = P::zero();
    for (r, e) in first.iter().copied().enumerate() {
        let projection = P::new(
            numerators[e]
                .coeffs()
                .iter()
                .enumerate()
                .map(|(j, c)| {
                    if j % period == r {
                        c.clone()
                    } else {
                        Q::zero()
                    }
                })
                .collect(),
        );
        leading_numerator = leading_numerator + projection * face.pow(last_layer - e);
    }
    let tentative = series_coefficients(&leading_numerator, &leading_denominator, start);
    let mut rows: Vec<P> = Vec::new();
    let mut correction = Vec::new();
    for (j, coefficient) in tentative.iter().enumerate() {
        let mut row = P::new(n.rows().get(j).cloned().unwrap_or_default());
        for k in 1..=j.min(d.rows().len() - 1) {
            row = row - P::new(d.rows()[k].clone()) * rows[j - k].clone();
        }
        correction.push(
            row.degree()
                .map_or_else(Q::zero, |degree| row.coeff(degree))
                - coefficient,
        );
        rows.push(row);
    }
    leading_numerator = leading_numerator + P::new(correction) * leading_denominator.clone();
    let leading =
        RationalGeneratingSeries::new(leading_numerator.coeffs(), leading_denominator.coeffs())?;
    let prefix_degrees = rows.iter().map(P::degree).collect();
    Ok(DegreeBudget {
        count: count.clone(),
        explanation: vec![format!(
            "Newton substitution t=u^(-{weight}), z=z*u^{growth}, numerator shift {shift}; earlier u-layer coefficients vanish beyond the checked finite prefix in each residue modulo {period}; first eventual layers {:?}, exact positive coefficient formulas and support seeds prove degree ({growth}*n+{:?})/{weight} for n>={start}", first, count.offsets)],
        data: EffectiveBkwDegreeCertificate {
            upper_bound: count,
            leading_numerator: leading.numerator,
            leading_denominator: leading.denominator,
            positive_top_residues: (0..period).collect(),
            newton: Some(EffectiveBkwNewtonDegree {
                weight_numerator: growth, weight_denominator: weight,
                numerator_shift: shift, start, layers, first_layers: first,
                prefix_degrees, reachability_seeds,
                positivity_supplier: if positive_cone {"nonnegative numerator and renewal kernel give positive rational layer expansions"} else {"nonnegative reduced layer numerators or exact binomial coefficient polynomial signs"}.into(),
                leading_sign_certificates,
            }),
        },
    })
}

fn binomial_positivity(
    n: &P,
    face: &P,
    power: usize,
    r: usize,
    step: usize,
) -> Result<(usize, IndexSignCertificate), RecurrenceError> {
    if face
        .coeffs()
        .iter()
        .enumerate()
        .any(|(j, c)| j != 0 && j != step && !c.is_zero())
    {
        return Err(failure("not a binomial Newton face"));
    }
    let ratio = -face.coeff(step);
    if !ratio.is_positive() {
        return Err(failure("binomial Newton ratio is not positive"));
    }
    let mut polynomial = P::zero();
    let mut first = 0;
    for (j, c) in n.coeffs().iter().enumerate().filter(|(_, c)| !c.is_zero()) {
        if j % step != r {
            return Err(failure("binomial layer has mixed residues"));
        }
        let shift = j / step;
        first = first.max(shift);
        let mut binomial = P::one();
        for i in 1..power {
            binomial = (binomial * P::new(vec![q(i as i64 - shift as i64), q(1)]))
                .scale(&Q::from_integer(i.into()).recip());
        }
        polynomial = polynomial + binomial.scale(&(c / ratio.pow(shift as i32)));
    }
    for k in first..=first + 64 {
        let decision = decide_index_polynomial_sign(
            polynomial.coeffs(),
            &IntegerRange::starting_at(k),
            IndexSignCondition::Positive,
        );
        if let IndexSignDecision::Holds { certificate } = decision {
            return Ok((r + step * k, *certificate));
        }
    }
    Err(failure("binomial layer positivity shift exceeded 64"))
}

pub(super) fn series_coefficients(n: &P, d: &P, count: usize) -> Vec<Q> {
    let mut coefficients = Vec::with_capacity(count);
    for j in 0..count {
        let mut value = n.coeff(j);
        for k in 1..=j.min(d.degree().unwrap_or(0)) {
            value -= d.coeff(k) * &coefficients[j - k];
        }
        coefficients.push(value / d.coeff(0));
    }
    coefficients
}

pub(super) fn value(count: &EffectiveBkwZeroCount, n: usize) -> Result<usize, RecurrenceError> {
    if count.divisor == 0 || count.offsets.is_empty() {
        return Err(failure("malformed quasi-linear zero count"));
    }
    let value = (count.slope as i128)
        .checked_mul(n as i128)
        .and_then(|v| v.checked_add(count.offsets[n % count.offsets.len()] as i128))
        .ok_or_else(|| failure("zero-count arithmetic overflow"))?;
    if value < 0 || value % count.divisor as i128 != 0 {
        return Err(failure(
            "zero-count law is negative or nonintegral at the cutoff",
        ));
    }
    usize::try_from(value / count.divisor as i128).map_err(|_| failure("zero-count overflow"))
}

/// Count strictly interior integer phase levels and exterior sign changes.
/// All quantities are affine on each residue of the common period. Checking
/// one representative and its increment proves the count for every later n.
pub(super) fn close_phase(out: &mut EffectiveBkwCertificate) -> Result<(), RecurrenceError> {
    let count = &out.zero_count;
    if count.divisor == 0 || count.offsets.is_empty() {
        return Err(failure(
            "phase degree divisor and residue period must be nonzero",
        ));
    }
    // The parity shortcut only gives an upper bound when its top
    // coefficient vanishes. Its original phase argument remains valid,
    // but this general two-exterior-sign closure needs an exact degree
    // and the positive sign of the actual leading coefficient.
    if out.degree_certificate.newton.is_none() && count.offsets != [0] {
        return Err(failure("general phase closure requires an exact positive leading-coefficient law; a parity upper bound is insufficient"));
    }
    let limits = &out.phase_limits;
    let l = rational(&limits.left_eigenvalue_argument)?;
    let r = rational(&limits.right_eigenvalue_argument)?;
    let a = rational(&limits.left_amplitude_argument)?;
    let b = rational(&limits.right_amplitude_argument)?;
    if l <= r {
        return Err(failure("phase interval has no positive asymptotic width"));
    }
    if (&l - &r) * Q::from_integer(count.divisor.into()) != Q::from_integer(count.slope.into()) {
        return Err(failure("phase and degree slopes disagree"));
    }
    if let Some(newton) = &out.degree_certificate.newton {
        out.cutoff = out.cutoff.max(newton.start);
    }
    let mut period = count.offsets.len();
    for d in [
        count.divisor,
        usize::try_from(l.denom().clone()).unwrap_or(1025),
        usize::try_from(r.denom().clone()).unwrap_or(1025),
    ] {
        let doubled = d
            .checked_mul(2)
            .filter(|d| *d != 0)
            .ok_or_else(|| failure("phase denominator doubling overflow"))?;
        period = period
            .checked_mul(doubled / gcd(period, doubled))
            .filter(|p| *p <= 1024)
            .ok_or_else(|| failure("phase residue period exceeds 1024"))?;
    }
    let mut proof = Vec::new();
    let mut right_escape = false;
    for residue in 0..period {
        let n = out
            .cutoff
            .checked_add((residue + period - out.cutoff % period) % period)
            .ok_or_else(|| failure("phase representative index overflow"))?;
        let mut previous: Option<usize> = None;
        let next = n
            .checked_add(period)
            .ok_or_else(|| failure("phase representative increment overflow"))?;
        for index in [n, next] {
            let index_q = Q::from_integer(index.into());
            let high = (&l * &index_q + &a).ceil() - q(1);
            let low = (&r * &index_q + &b).floor() + q(1);
            let samples = &high - &low + q(1);
            let degree = value(count, index)?;
            if samples < q(1) {
                if degree == 0 {
                    // The degree/positive-leading-coefficient certificate
                    // already proves this row is a nonzero constant. The
                    // next index in this class starts the affine count.
                    previous = Some(0);
                    continue;
                }
                return Err(failure("phase cutoff has no interior sign sample"));
            }
            let parity = |j: &Q| {
                (j.to_integer() % num_bigint::BigInt::from(2)) != num_bigint::BigInt::zero()
            };
            let left = parity(&high) != (degree % 2 != 0);
            let right = parity(&low);
            right_escape |= right;
            let roots = samples - q(1) + q(i64::from(left) + i64::from(right));
            if roots != Q::from_integer(degree.into()) {
                return Err(failure(format!(
                    "phase/exterior count does not close degree in residue {residue} at n={index}"
                )));
            }
            if let Some(old) = previous {
                let increment = Q::from_integer((degree - old).into());
                if (&l - &r) * Q::from_integer(period.into()) != increment {
                    return Err(failure(
                        "phase and degree slopes disagree on a residue class",
                    ));
                }
            }
            previous = Some(degree);
        }
        proof.push(format!(
            "r={residue}: degree {}, samples and exterior signs close",
            value(count, n)?
        ));
    }
    if right_escape {
        out.root_region = EffectiveBkwRootRegion::WholeRealLine;
    }
    out.exterior.push(format!("strict integer phase levels, positive leading coefficients and both exterior signs checked on period {period}; affine increments agree: {}",proof.join("; ")));
    Ok(())
}
