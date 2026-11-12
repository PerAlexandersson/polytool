//! Exact interval geometry and root implications.

use super::*;

pub(super) fn conventions() -> Vec<String> {
    [
        "rows are ascending coefficient vectors; P_n are the source rows and the transformed \
         rows are named by the transformation's components; n is the original index",
        "real_rooted: every root is real; real_rooted_in I: every root is real and lies in I; \
         real_roots_in I: every real root lies in I (non-real roots allowed)",
        ZERO_CASES,
        "intervals have rational, infinite (omitted) or index-dependent endpoints \
         value + per_index * n, each closed or open",
        "evidence: finite_exact_observation (exact on the listed rows only), \
         supplied_assumption (from the request, checked on rows only), proposed_interval \
         (candidate window, not proved), completed_recurrence_induction (closed induction over \
         the transported rule for every index of the domain), certified_family_invariant \
         (transferred to every index by a map theorem whose conditions are certified)",
        "exact root isolation proves statements about the isolated polynomials only, never \
         about the family",
    ]
    .into_iter()
    .map(String::from)
    .collect()
}

// ---------------------------------------------------------------------------
// Map geometry: exact interval images
// ---------------------------------------------------------------------------

/// The empty interval `(0, 0)`.
fn empty_interval() -> RealInterval {
    RealInterval {
        lower: Some(Q::zero()),
        lower_closed: false,
        upper: Some(Q::zero()),
        upper_closed: false,
    }
}

/// Image of `interval` under `y -> (p y + q)/(r y + s)` (`p s - q r != 0`),
/// or `None` when the pole lies in the interior (the image is not an interval).
/// The pole itself maps to infinity and is dropped.
pub(crate) fn mobius_interval_image(
    interval: &RealInterval,
    p: &Q,
    q: &Q,
    r: &Q,
    s: &Q,
) -> Option<RealInterval> {
    if interval.is_empty() {
        return Some(interval.clone());
    }
    if !r.is_zero() {
        let pole = -(s / r);
        let is_lower = interval.lower.as_ref() == Some(&pole);
        let is_upper = interval.upper.as_ref() == Some(&pole);
        if is_lower && is_upper {
            return Some(empty_interval());
        }
        if !is_lower && !is_upper && interval.contains(&pole) {
            return None;
        }
    }
    let image = |end: &Option<Q>, closed: bool| -> (Option<Q>, bool) {
        match end {
            None if r.is_zero() => (None, false),
            None => (Some(p / r), false),
            Some(y) => {
                let denominator = r * y + s;
                if denominator.is_zero() {
                    (None, false)
                } else {
                    (Some((p * y + q) / denominator), closed)
                }
            }
        }
    };
    let low = image(&interval.lower, interval.lower_closed);
    let high = image(&interval.upper, interval.upper_closed);
    let increasing = (p * s - q * r).is_positive();
    let (lower, upper) = if increasing { (low, high) } else { (high, low) };
    Some(RealInterval {
        lower_closed: lower.0.is_some() && lower.1,
        lower: lower.0,
        upper_closed: upper.0.is_some() && upper.1,
        upper: upper.0,
    })
}

// One short-lived value per analysis; boxing the Möbius parameters buys nothing.
#[allow(clippy::large_enum_variant)]
pub(super) enum Geometry {
    Monomial {
        power: usize,
        multiply: bool,
    },
    Mobius {
        a: Q,
        b: Q,
        c: Q,
        e: Q,
        theorem: &'static str,
    },
    Gamma,
    Symmetric,
    Unsupported,
}

pub(in crate::recurrence::root_domain) fn q(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

impl Geometry {
    pub(super) fn of(transformation: &RecurrenceTransformation) -> Result<Self, RecurrenceError> {
        let parse = |text: &str| parse_rational_coeff(text);
        Ok(match transformation {
            RecurrenceTransformation::MultiplyPower { power } => Self::Monomial {
                power: *power,
                multiply: true,
            },
            RecurrenceTransformation::DividePower { power } => Self::Monomial {
                power: *power,
                multiply: false,
            },
            RecurrenceTransformation::Reversal { .. } => Self::Mobius {
                a: q(0),
                b: q(1),
                c: q(1),
                e: q(0),
                theorem: "real_mobius_equivalence",
            },
            RecurrenceTransformation::FPolynomial { .. } => Self::Mobius {
                a: q(1),
                b: q(0),
                c: q(1),
                e: q(1),
                theorem: "real_mobius_equivalence",
            },
            RecurrenceTransformation::Affine { a, b } => Self::Mobius {
                a: parse(a)?,
                b: parse(b)?,
                c: q(0),
                e: q(1),
                theorem: "real_affine_equivalence",
            },
            RecurrenceTransformation::Mobius { a, b, c, e, .. } => Self::Mobius {
                a: parse(a)?,
                b: parse(b)?,
                c: parse(c)?,
                e: parse(e)?,
                theorem: "real_mobius_equivalence",
            },
            RecurrenceTransformation::Gamma { .. } => Self::Gamma,
            RecurrenceTransformation::SymmetricDecomposition { .. } => Self::Symmetric,
            _ => Self::Unsupported,
        })
    }

    /// Whether real-rootedness is equivalent and intervals map exactly.
    pub(super) fn is_scalar_equivalence(&self) -> bool {
        matches!(self, Self::Monomial { .. } | Self::Mobius { .. })
    }

    /// Interval containing the roots of the transformed row when the roots of
    /// the source row lie in `interval` (scalar equivalences only).
    pub(super) fn forward(&self, interval: &RealInterval) -> Option<RealInterval> {
        match self {
            Self::Monomial { power, multiply } => Some(if *multiply && *power > 0 {
                interval.hull(&RealInterval::point(Q::zero()))
            } else {
                interval.clone()
            }),
            Self::Mobius { a, b, c, e, .. } => {
                let image = mobius_interval_image(interval, e, &-b.clone(), &-c.clone(), a)?;
                Some(if c.is_zero() {
                    image
                } else {
                    // The pole -e/c of phi is a root of multiplicity d_n - deg P_n.
                    image.hull(&RealInterval::point(-(e / c)))
                })
            }
            _ => None,
        }
    }

    /// Interval containing the roots of the source row when the roots of the
    /// transformed row lie in `interval`; `full_degree` means `deg Q_n = d_n`.
    pub(super) fn backward(
        &self,
        interval: &RealInterval,
        full_degree: bool,
    ) -> Option<RealInterval> {
        match self {
            Self::Monomial { power, multiply } => Some(if !*multiply && *power > 0 {
                interval.hull(&RealInterval::point(Q::zero()))
            } else {
                interval.clone()
            }),
            Self::Mobius { a, b, c, e, .. } => {
                let image = mobius_interval_image(interval, a, b, c, e)?;
                Some(if c.is_zero() || full_degree {
                    image
                } else {
                    // deg Q_n < d_n gives roots of P_n at a/c.
                    image.hull(&RealInterval::point(a / c))
                })
            }
            _ => None,
        }
    }

    pub(super) fn theorem(&self) -> &'static str {
        match self {
            Self::Monomial { .. } => "monomial_factor",
            Self::Mobius { theorem, .. } => theorem,
            Self::Gamma => "gamma_quarter_bound",
            Self::Symmetric => "hermite_kakeya_obreschkoff_symmetric_pair",
            Self::Unsupported => "none",
        }
    }
}

// ---------------------------------------------------------------------------
// Implications
// ---------------------------------------------------------------------------

pub(in crate::recurrence::root_domain) fn condition_hypotheses(
    conditions: &[TransformConditionStatus],
) -> Vec<RootHypothesis> {
    conditions
        .iter()
        .map(|c| RootHypothesis {
            condition: c.condition.clone(),
            status: if c.status == CERTIFIED_ALL_INDICES {
                HypothesisStatus::CertifiedAllIndices
            } else {
                HypothesisStatus::CheckedSuppliedRows
            },
            detail: c.detail.clone(),
        })
        .collect()
}

pub(in crate::recurrence::root_domain) fn uncertified(
    conditions: &[TransformConditionStatus],
) -> Vec<String> {
    conditions
        .iter()
        .filter(|c| c.status != CERTIFIED_ALL_INDICES)
        .map(|c| format!("{} (checked on the supplied rows only)", c.condition))
        .collect()
}

/// The theorem adapter of a transformation, with every hypothesis
/// [`HypothesisStatus::NotChecked`].  Reports fill in the statuses.
///
/// ```
/// use polytool::recurrence::*;
///
/// let f = RecurrenceTransformation::FPolynomial {
///     schedule: DegreeSchedule::new(1, 0, 1),
/// };
/// let implication = root_implication(&f).unwrap();
/// assert_eq!(implication.theorem, "real_mobius_equivalence");
/// assert_eq!(implication.direction, RootImplicationDirection::Equivalent);
/// ```
pub fn root_implication(
    transformation: &RecurrenceTransformation,
) -> Result<RootImplication, RecurrenceError> {
    transformation.validate()?;
    let geometry = Geometry::of(transformation)?;
    let record = transformation.record();
    let hypotheses = record
        .image_conditions
        .iter()
        .map(|condition| RootHypothesis {
            condition: condition.clone(),
            status: HypothesisStatus::NotChecked,
            detail: "not checked".into(),
        })
        .collect();
    let mut implication = RootImplication {
        theorem: geometry.theorem().to_string(),
        direction: RootImplicationDirection::Equivalent,
        source_property: "P_n real-rooted".into(),
        target_property: "Q_n real-rooted".into(),
        statement: String::new(),
        hypotheses,
        degree_bounds: record.degree_schedule.clone(),
        coefficient_signs: None,
        endpoint_conventions: Vec::new(),
        zero_cases: ZERO_CASES.to_string(),
        index_domain: record.index_domain.clone(),
        converse: None,
    };
    match &geometry {
        Geometry::Monomial { power, multiply } => {
            implication.statement = format!(
                "{}: the roots of Q_n are those of P_n with the multiplicity of 0 {} by {power}; \
                 real-rootedness and every interval containing 0 are preserved both ways",
                record.row_map,
                if *multiply { "raised" } else { "lowered" }
            );
            implication.endpoint_conventions = vec![format!(
                "the root 0 of {} has multiplicity at least {power}",
                if *multiply { "Q_n" } else { "P_n" }
            )];
        }
        Geometry::Mobius { a, b, c, e, .. } => {
            let pole = (!c.is_zero()).then(|| format_rational_coeff(&-(e / c)));
            let infinity_image = (!c.is_zero()).then(|| format_rational_coeff(&(a / c)));
            let (a, b, c, e) = (
                format_rational_coeff(a),
                format_rational_coeff(b),
                format_rational_coeff(c),
                format_rational_coeff(e),
            );
            implication.statement = format!(
                "{}: a real Möbius map with a e - b c != 0 is a bijection of the extended real \
                 line mapping non-real points to non-real points; the roots of Q_n are \
                 psi(rho) = (({e}) rho - ({b}))/(({a}) - ({c}) rho) for the roots rho of P_n, so \
                 P_n is real-rooted iff Q_n is, and intervals map by psi",
                record.row_map
            );
            if let (Some(pole), Some(infinity_image)) = (pole, infinity_image) {
                implication.endpoint_conventions = vec![
                    format!(
                        "the pole x = {pole} of the substitution is a root of Q_n of \
                         multiplicity d_n - deg P_n (missing top degree of P_n)"
                    ),
                    format!(
                        "a root rho = {infinity_image} of P_n lowers deg Q_n below d_n; \
                         conversely deg Q_n < d_n gives roots of P_n at {infinity_image}"
                    ),
                ];
            }
            if transformation.kind_name() == "f_polynomial" {
                implication
                    .endpoint_conventions
                    .push("the root 0 of Q_n has the multiplicity of the root 0 of P_n".into());
            }
        }
        Geometry::Gamma => {
            implication.target_property =
                "Q_n (gamma coordinates) real-rooted with every root <= 1/4".into();
            implication.statement = "P_n = (1 + x)^(d_n) Q_n(x/(1 + x)^2) with \
                                     deg Q_n <= floor(d_n/2); each root y of Q_n contributes the \
                                     factor x - y (1 + x)^2 of discriminant 1 - 4y, so P_n is \
                                     real-rooted iff every root of Q_n is real and at most 1/4"
                .into();
            implication.endpoint_conventions = vec![
                "a root 1/4 of Q_n gives the double root x = 1 of P_n".into(),
                "a root 0 of Q_n gives the root x = 0 of P_n and lowers deg P_n".into(),
                "(1 + x)^(d_n - 2 deg Q_n) contributes the root -1".into(),
            ];
            implication.converse = Some(
                "real-rootedness of Q_n alone is insufficient: Q = 1 - z at d = 2 gives \
                 P = 1 + x + x^2"
                    .into(),
            );
        }
        Geometry::Symmetric => {
            implication.direction = RootImplicationDirection::TargetToSource;
            implication.source_property = "P_n = a_n + x b_n real-rooted".into();
            implication.target_property =
                "a_n and x b_n real-rooted and interlacing (in some order)".into();
            implication.statement = "Hermite-Kakeya-Obreschkoff: if a_n and x b_n are \
                                     real-rooted and interlace, every real combination, in \
                                     particular a_n + x b_n, is real-rooted (with nonnegative \
                                     coefficients this is the Brändén-Solus condition that b_n \
                                     interlaces a_n)"
                .into();
            implication.hypotheses.push(RootHypothesis {
                condition: "a_n and x b_n are real-rooted and interlace".into(),
                status: HypothesisStatus::NotChecked,
                detail: "not checked".into(),
            });
            implication.converse = Some(
                "fails: P = (1 + 2x)^2 at d = 2 is real-rooted but a = 1 + x + x^2 is not; and \
                 real-rooted components alone do not suffice: a = (1 + x)^3, b = x at d = 3 give \
                 P = 1 + 3x + 4x^2 + x^3, which is not real-rooted"
                    .into(),
            );
        }
        Geometry::Unsupported => {
            implication.statement = "no root implication is recorded for this map".into();
            implication.direction = RootImplicationDirection::SourceToTarget;
            implication.target_property = "none".into();
        }
    }
    Ok(implication)
}
