//! Shared degree-law reversal, used by index sections and theorem discovery.
//!
//! If `P_m = sum a_k(x) P_(m-k)`, `deg a_k <= d k`, and every prescribed
//! seed has degree at most `d m+e`, induction gives that bound for all rows.
//! Then `S_m=x^(d m+e) P_m(1/x)` has multipliers `x^(d k) a_k(1/x)`.
//! Reversal is an involution, preserves zero rows and real-rootedness, and
//! can add roots only at zero. It need not preserve simplicity or interlacing.
//!
//! Discovery tries slopes from the smallest admissible positive slope through
//! [`MAXIMUM_REVERSAL_DEGREE_SLOPE`] (3), with the smallest seed offset.
//! This is a search bound, not a restriction of the closure theorem. Each
//! transformed record is sent to the nonrecursive theorem routes, in registry
//! order, with their usual start retries. No section or reversal search is
//! nested inside that menu. Index sections use the same transform and replay,
//! with their existing smallest-slope search (including slope zero).

use super::*;
use crate::Polynomial;

/// Largest slope tried by the shared reversal discovery closure.
pub const MAXIMUM_REVERSAL_DEGREE_SLOPE: usize = 3;

/// A degree-law coordinate change; indices `m` are relative to the first row.
/// This is also the certificate shape used by index-section lifts.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum DegreeLawTransform {
    /// `S_m = P_(first+m)`.
    Identity,
    /// `S_m(x)=x^(degree_slope*m+degree_offset) P_(first+m)(1/x)`.
    Reversal {
        /// `d`: every lag multiplier has degree at most `d k`.
        degree_slope: usize,
        /// `e`: every prescribed seed has degree at most `d m+e`.
        degree_offset: usize,
    },
}

impl DegreeLawTransform {
    /// Replay the all-index degree bound and construct the exact transformed
    /// rule and all seeds. This one check is shared by section and reversal
    /// certificates; no inferred degree fit or finite extrapolation is used.
    ///
    /// # Errors
    ///
    /// Unsupported rules, an insufficient bound, missing seeds or overflow.
    pub fn apply(
        self,
        rule: &Recurrence,
        first: usize,
        seeds: &[Vec<Q>],
    ) -> Result<RecurrenceJson, RecurrenceError> {
        let multipliers =
            constant_coefficient_multipliers(rule).map_err(RecurrenceError::Unsupported)?;
        if seeds.len() < multipliers.len() {
            return Err(RecurrenceError::InvalidInput(
                "reversal needs every lag seed".into(),
            ));
        }
        let reversal = match self {
            Self::Identity => None,
            Self::Reversal {
                degree_slope,
                degree_offset,
            } => Some((degree_slope, degree_offset)),
        };
        let multipliers = match reversal {
            None => multipliers,
            Some((slope, _)) => reversed_multipliers(&multipliers, slope).ok_or_else(|| {
                RecurrenceError::VerificationFailed(
                    "a lag multiplier exceeds the degree bound".into(),
                )
            })?,
        };
        let rows = transformed_rows(seeds, reversal).ok_or_else(|| {
            RecurrenceError::VerificationFailed("a seed exceeds the degree bound".into())
        })?;
        let transformed = Recurrence {
            terms: multipliers
                .iter()
                .enumerate()
                .filter(|(_, p)| !p.is_zero())
                .map(|(i, p)| {
                    RecurrenceTerm::new(
                        i + 1,
                        0,
                        BivarPoly::from_coefficients_normalized(vec![p.coeffs().to_vec()]),
                    )
                })
                .collect(),
            denominator: None,
            inhomogeneous: None,
        };
        Ok(RecurrenceJson::from_recurrence_rational(
            &transformed,
            first,
            &rows
                .iter()
                .map(|p| vec_poly::trim(p.coeffs()))
                .collect::<Vec<_>>(),
            None,
        ))
    }
}

/// A theorem proved after the shared degree-law reversal.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DegreeLawReversalCertificate {
    /// The shared transform, always a reversal for this closure.
    pub transform: DegreeLawTransform,
    /// A proved theorem report on the exact transformed rule and all seeds.
    pub transformed_proof: Box<RealRootednessTheoremReport>,
    /// Search bound passed to the inner theorem's induction-start retries.
    pub maximum_start_delay: usize,
}

impl DegreeLawReversalCertificate {
    /// Replay the shared transform, the inner theorem and the root transfer.
    ///
    /// # Errors
    ///
    /// A changed transform, inner family, theorem proof or search bound.
    pub fn verify(
        &self,
        rule: &Recurrence,
        first: usize,
        seeds: &[Vec<Q>],
    ) -> Result<(), RecurrenceError> {
        if !matches!(self.transform, DegreeLawTransform::Reversal { .. })
            || !eligible_inner_route(route(self.transformed_proof.theorem))
            || !self.transformed_proof.is_proved()
            || self.transformed_proof.start_delay > self.maximum_start_delay
            || self.maximum_start_delay > index_section::MAXIMUM_PARENT_START_DELAY
            || self.transform.apply(rule, first, seeds)? != self.transformed_proof.recurrence
        {
            return Err(RecurrenceError::VerificationFailed(
                "invalid degree-law reversal certificate".into(),
            ));
        }
        self.transformed_proof.verify()
    }
}

fn eligible_inner_route(route: &TheoremRoute) -> bool {
    route.start_policy == StartPolicy::MovesStart
        && route.theorem != RealRootednessTheorem::DegreeLawReversal
}

/// Prove a family by a theorem on its degree-law reversal (module documentation).
///
/// # Errors
///
/// Invalid recurrence data or a failure of an inner exact theorem check.
pub fn check_degree_law_reversal_theorem(
    rule: &Recurrence,
    first: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    check_with_first_proving_start(
        RealRootednessTheorem::DegreeLawReversal,
        rule,
        first,
        seeds,
        options,
    )
}

pub(super) const ROUTE: TheoremRoute = TheoremRoute {
    theorem: RealRootednessTheorem::DegreeLawReversal,
    check_at_supplied_seeds: reversal_report,
    start_policy: StartPolicy::OwnStarts,
    is_preferred: any_proof,
    certificate_name: Some("degree-law reversal"),
    discovery: DiscoveryPolicy::WhenApplicable,
    in_residue_class_menu: false,
    in_index_section_menu: false,
    parent: ParentPolicy {
        stage: ParentStage::WhileUnproved,
        region: ParentRegion::Unrestricted,
        record_only_if_applicable: true,
    },
};

fn reversal_report(
    rule: &Recurrence,
    first: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    let report = RealRootednessTheoremReport::new(
        RealRootednessTheorem::DegreeLawReversal,
        "Polynomial reversal: reciprocal nonzero roots, with possible additional zeros at 0",
        rule,
        first,
        seeds,
        options,
    );
    let multipliers = match constant_coefficient_multipliers(rule) {
        Ok(m) => m,
        Err(reason) => return Ok(report.not_applicable(reason)),
    };
    for degree_slope in degree_slope(&multipliers).max(1)..=MAXIMUM_REVERSAL_DEGREE_SLOPE {
        let transform = DegreeLawTransform::Reversal {
            degree_slope,
            degree_offset: degree_offset(seeds, degree_slope),
        };
        let Ok(transformed) = transform.apply(rule, first, seeds) else {
            continue;
        };
        let (inner_rule, inner_first, inner_seeds) = transformed.to_recurrence_parts()?;
        let inner_options = RealRootednessTheoremOptions {
            root_region: None,
            ..options.clone()
        };
        for candidate in ROUTES.iter().filter(|r| eligible_inner_route(r)) {
            let proof = check_with_first_proving_start(
                candidate.theorem,
                &inner_rule,
                inner_first,
                &inner_seeds,
                &inner_options,
            )?;
            if !proof.is_proved() {
                continue;
            }
            let inner = proof.conclusion.as_ref().unwrap();
            let region = inner
                .root_region
                .as_ref()
                .and_then(|r| r.at(0).ok().filter(|_| r.is_index_independent()))
                .and_then(|r| reversed_region(&r));
            if options.root_region.as_ref().is_some_and(|requested| {
                region
                    .as_ref()
                    .is_none_or(|certified| !requested_region_contains(requested, certified))
            }) {
                continue;
            }
            let mut conclusion = inner.clone();
            conclusion.root_region = region.as_ref().map(RootInterval::from_real);
            conclusion.simple_roots = false;
            conclusion.interlacing = None;
            conclusion.interlacing_from = None;
            conclusion.statement = format!("degree-law reversal S_m=x^({degree_slope} m+{}) P_({first}+m)(1/x); reciprocal nonzero roots and possible roots at 0; {}", degree_offset(seeds, degree_slope), inner.statement);
            let mut report = report;
            report.normal_form = conclusion.statement.clone();
            report.hypotheses.push(exact_hypothesis(
                "degree_law_reversal",
                "polynomial reversal preserves real-rootedness and zero rows",
                Ok(report.normal_form.clone()),
            ));
            report.hypotheses.push(exact_hypothesis(
                "transformed_theorem",
                "a named theorem proves the reversed family",
                Ok(format!("{}: {}", proof.theorem.as_str(), proof.reference)),
            ));
            report.certificate = Some(RouteCertificate::DegreeLawReversal(Box::new(
                DegreeLawReversalCertificate {
                    transform,
                    transformed_proof: Box::new(proof),
                    maximum_start_delay: options.maximum_start_delay,
                },
            )));
            return Ok(report.finish(conclusion));
        }
    }
    Ok(report
        .not_applicable("no theorem proves a degree-law reversal within the shared slope bound"))
}

pub(super) fn verify_certificate_in_report(
    certificate: &DegreeLawReversalCertificate,
    report: &RealRootednessTheoremReport,
) -> Result<(), RecurrenceError> {
    let (rule, first, seeds) = report.recurrence.to_recurrence_parts()?;
    let (rule, _) = rule_with_x_free_left_hand_factor(&rule)?;
    certificate.verify(&rule, first, &seeds)
}

/// The multipliers `a_1, ..., a_r` of a rule `P_n = sum_k a_k(x) P_(n-k)`
/// (lags without a term give `0`), or why the rule has another form.
pub(super) fn constant_coefficient_multipliers(
    rule: &Recurrence,
) -> Result<Vec<Polynomial<Q>>, String> {
    let shape = recurrence_shape(rule)?;
    if shape.denominator.len() > 1 {
        return Err("the left-hand factor depends on n".into());
    }
    let leading = shape.denominator[0].clone();
    let order = shape.terms.keys().map(|&(lag, _)| lag).max().unwrap_or(0);
    if order == 0 {
        return Err("the rule has no term of positive lag".into());
    }
    let mut multipliers = vec![Polynomial::<Q>::zero(); order];
    for (&(lag, derivative), coefficient) in &shape.terms {
        if lag == 0 {
            return Err("the rule has a term of lag 0".into());
        }
        if derivative != 0 {
            return Err(format!(
                "the rule has a derivative term (lag {lag}, derivative order {derivative})"
            ));
        }
        if coefficient
            .iter()
            .skip(1)
            .any(|row| row.iter().any(|c| !c.is_zero()))
        {
            return Err(format!("the multiplier of lag {lag} depends on n"));
        }
        let row = coefficient.first().cloned().unwrap_or_default();
        multipliers[lag - 1] = Polynomial::new(vec_poly::scale(&row, &leading.recip()));
    }
    Ok(multipliers)
}

/// The smallest `d` with `deg a_k <= d k` for every `k`.
pub(super) fn degree_slope(multipliers: &[Polynomial<Q>]) -> usize {
    multipliers
        .iter()
        .enumerate()
        .filter_map(|(lag, a)| a.degree().map(|degree| degree.div_ceil(lag + 1)))
        .max()
        .unwrap_or(0)
}

/// `b_k(x) = x^(d k) a_k(1/x)`, or `None` when some `deg a_k > d k`.
pub(super) fn reversed_multipliers(
    multipliers: &[Polynomial<Q>],
    slope: usize,
) -> Option<Vec<Polynomial<Q>>> {
    multipliers
        .iter()
        .enumerate()
        .map(|(lag, a)| {
            if a.is_zero() {
                Some(Polynomial::<Q>::zero())
            } else {
                a.reverse_with_degree(slope.checked_mul(lag + 1)?)
            }
        })
        .collect()
}

/// The transformed rows `S_m` of `rows[m] = P_(n_0 + s + m)`: unchanged,
/// or `x^(d m + e) P(1/x)` for `Some((d, e))`; `None` when a reversed row is
/// not a polynomial.
pub(super) fn transformed_rows(
    rows: &[Vec<Q>],
    reversal: Option<(usize, usize)>,
) -> Option<Vec<Polynomial<Q>>> {
    rows.iter()
        .enumerate()
        .map(|(m, row)| {
            let row = Polynomial::new(row.clone());
            match reversal {
                None => Some(row),
                Some(_) if row.is_zero() => Some(row),
                Some((slope, offset)) => {
                    row.reverse_with_degree(slope.checked_mul(m)?.checked_add(offset)?)
                }
            }
        })
        .collect()
}

/// The smallest `e >= 0` with `deg rows[m] <= d m + e`.
pub(super) fn degree_offset(rows: &[Vec<Q>], slope: usize) -> usize {
    rows.iter()
        .enumerate()
        .filter_map(|(m, row)| {
            Polynomial::new(row.clone())
                .degree()
                .map(|degree| degree.saturating_sub(slope.saturating_mul(m)))
        })
        .max()
        .unwrap_or(0)
}
/// The region of the rows `x^D S(1/x)` from a region of the rows `S`
/// (module documentation, "Root region").
pub(super) fn reversed_region(region: &RealInterval) -> Option<RealInterval> {
    if let Some(upper) = region.upper.as_ref().filter(|u| !u.is_positive()) {
        let (lower, lower_closed) = if upper.is_zero() {
            (None, false)
        } else {
            (Some(upper.recip()), region.upper_closed)
        };
        return Some(RealInterval {
            lower,
            lower_closed,
            upper: Some(Q::zero()),
            upper_closed: true,
        });
    }
    if let Some(lower) = region.lower.as_ref().filter(|l| !l.is_negative()) {
        let (upper, upper_closed) = if lower.is_zero() {
            (None, false)
        } else {
            (Some(lower.recip()), region.lower_closed)
        };
        return Some(RealInterval {
            lower: Some(Q::zero()),
            lower_closed: true,
            upper,
            upper_closed,
        });
    }
    None
}
