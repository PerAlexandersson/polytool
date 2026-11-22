//! Sharp positive scale interval for the three-state lag-sum template.
use super::*;

/// Exact interval solving a*kappa^2-e*kappa+d<=0, kappa>0.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ScaledRefinementInterval {
    /// a,b,d,e of the scalar recurrence, in that order.
    pub parameters: [String; 4],
    /// Whether the universal template admits a positive scale.
    pub feasible: bool,
    /// Lower endpoint (zero when the positivity boundary is open).
    pub lower: Option<QuadraticNumber>,
    /// Upper endpoint; absent for an unbounded feasible interval.
    pub upper: Option<QuadraticNumber>,
    /// Whether the lower endpoint is included (zero is always excluded).
    pub lower_closed: bool,
    /// A positive rational scale in the interval, when feasible.
    pub rational_scale: Option<String>,
    /// Sharpness concerns universal preservation of this template, not the individual orbit.
    pub scope: String,
}

/// Solve the universal-preservation scale condition exactly, including zero
/// parameters and irrational endpoints. Negative a,b,d are infeasible for
/// this nonnegative-entry template; e is arbitrary.
///
/// # Errors
/// Invalid exact rational strings.
pub fn solve_scaled_refinement_interval(
    parameters: &[String; 4],
) -> Result<ScaledRefinementInterval, RecurrenceError> {
    let v = numbers(parameters)?;
    let (a, b, d, e) = (&v[0], &v[1], &v[2], &v[3]);
    let mut r = ScaledRefinementInterval { parameters: parameters.clone(), feasible: false, lower: None, upper: None, lower_closed: false, rational_scale: None, scope: "Universal preservation of (P_n, P_n+kappa*x*P_(n-1), P_(n+1)); not necessity for an individual polynomial family".into() };
    if a.is_negative() || b.is_negative() || d.is_negative() {
        return Ok(r);
    }
    let (lower, upper, scale) = if a.is_zero() {
        if e.is_positive() {
            (QuadraticNumber::rational(d / e), None, (d + Q::one()) / e)
        } else if e.is_zero() && d.is_zero() {
            (QuadraticNumber::rational(Q::zero()), None, Q::one())
        } else {
            return Ok(r);
        }
    } else {
        let disc = e * e - Q::from_integer(4.into()) * a * d;
        if !e.is_positive() || disc.is_negative() {
            return Ok(r);
        }
        let two_a = Q::from_integer(2.into()) * a;
        (
            QuadraticNumber::new(e / &two_a, -two_a.recip(), disc.clone())?,
            Some(QuadraticNumber::new(e / &two_a, two_a.recip(), disc)?),
            e / two_a,
        )
    };
    r.feasible = true;
    r.lower_closed = !d.is_zero();
    r.lower = Some(lower);
    r.upper = upper;
    r.rational_scale = Some(scale.to_string());
    Ok(r)
}

/// Intersect the sharp positive scale intervals for every supplied step.
/// This gives a common scale for exactly those parameter tuples. Uniform
/// index-dependent families still require a symbolic all-index certificate.
///
/// # Errors
/// Invalid rational strings or more than 256 parameter tuples.
pub fn solve_common_scaled_refinement_interval(
    parameters: &[[String; 4]],
) -> Result<Option<AlgebraicInterval>, RecurrenceError> {
    if parameters.len() > 256 {
        return Err(invalid("common scale input exceeds 256 steps"));
    }
    let mut interval = AlgebraicInterval {
        lower: Some(AlgebraicEndpoint {
            value: QuadraticNumber::rational(Q::zero()),
            closed: false,
        }),
        upper: None,
    };
    for step in parameters {
        let result = solve_scaled_refinement_interval(step)?;
        if !result.feasible {
            return Ok(None);
        }
        let next = AlgebraicInterval {
            lower: result.lower.map(|value| AlgebraicEndpoint {
                value,
                closed: result.lower_closed,
            }),
            upper: result.upper.map(|value| AlgebraicEndpoint {
                value,
                closed: true,
            }),
        };
        interval = interval.intersection(&next);
        if let (Some(a), Some(b)) = (&interval.lower, &interval.upper) {
            match a.value.cmp_exact(&b.value) {
                std::cmp::Ordering::Greater => return Ok(None),
                std::cmp::Ordering::Equal if !a.closed || !b.closed => return Ok(None),
                _ => {}
            }
        }
    }
    Ok(Some(interval))
}
