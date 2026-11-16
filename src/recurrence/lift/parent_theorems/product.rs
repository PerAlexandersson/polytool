//! Product induction, independent of named sequences or catalogue attestations.
//!
//! For `B(n) P_n = A(n,x) P_(n-k)`, nonvanishing `B` and nonzero real-rooted
//! seeds and multipliers imply nonzero real-rooted rows by induction in each
//! residue class. Multipliers are checked uniformly, never by sampling indices.
//! The Lean step is `RealRooted.productSequence_splits`; the residue-class
//! reduction, scalar division and exact multiplier certificates need replay.

use super::*;

/// Check product induction for `B(n) P_n = A(n,x) P_(n-k)`, `k >= 1`.
///
/// Supported multipliers are separable `g(n) h(x)` of arbitrary degree,
/// degree at most two in `x`, or products of these detected by bounded exact
/// bivariate factorization when the leading `x` coefficient is constant.
/// A quadratic uses its discriminant on the entire integer index range.
/// Nonvanishing is checked even when the degree drops. Requested open root
/// regions are supported for separable factors and factors of fixed degree
/// at most two. Unsupported factors give `unresolved`, never a proof.
///
/// The conclusion concerns the sequence *defined* by the supplied recurrence
/// and seeds. It asserts neither interlacing nor validity of a fitted
/// recurrence for another sequence. Zero seeds or multipliers are currently
/// rejected because the report promises nonzero rows.
///
/// # Errors
///
/// Errors parsing or normalizing the recurrence and seeds.
pub fn check_product_recurrence_theorem(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    check_real_rootedness_theorem(
        RealRootednessTheorem::ProductRecurrence,
        rule,
        first_index,
        seeds,
        options,
    )
}

/// Whether the sequence defined by `rule` and `seeds` is *trivially*
/// real-rooted: the rule has the product form `B(n) P_n = A(x,n) P_(n-k)` and
/// [`classify_product_recurrence`] certifies its real-rooted range or zero
/// residue classes (no root region requested).
///
/// This classifies the proof of a recursively defined sequence; it is never
/// a property read off a single polynomial or a fitted recurrence.  The same
/// direct flag is [`RealRootednessTheoremReport::trivially_real_rooted`];
/// [`classify_product_recurrence`] also covers zero seeds through the shared
/// tail and residue-class proofs.
///
/// # Errors
///
/// Errors parsing or normalizing the recurrence and seeds.
pub fn is_trivially_real_rooted(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
) -> Result<bool, RecurrenceError> {
    Ok(classify_product_recurrence(rule, first_index, seeds)?.is_some())
}

/// Classify a product recurrence, including zero initial rows.
///
/// Returns the proved nonzero real-rooted range or its zero residue classes.
/// The direct product check runs first. For zero seeds, the same residue-class
/// check and tail construction used by proof discovery are reused; nonzero
/// initial rows must themselves be real-rooted. This does not prove that a
/// fitted rule holds for an independently defined family.
///
/// # Errors
/// Returns errors parsing or normalizing the recurrence and seeds.
pub fn classify_product_recurrence(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
) -> Result<Option<RealRootednessTheoremConclusion>, RecurrenceError> {
    let options = RealRootednessTheoremOptions::default();
    first_index
        .checked_add(seeds.len())
        .and_then(|n| n.checked_add(options.maximum_start_delay))
        .and_then(|n| n.checked_add(GENERATED_ROWS))
        .filter(|n| *n <= i64::MAX as usize)
        .ok_or_else(|| {
            RecurrenceError::InvalidInput("product classification index range overflows".into())
        })?;
    let direct = check_product_recurrence_theorem(rule, first_index, seeds, &options)?;
    if direct.is_proved() {
        return Ok(direct.conclusion);
    }
    let Ok(shape) = recurrence_shape(rule) else {
        return Ok(None);
    };
    if shape.terms.len() != 1
        || !shape
            .terms
            .keys()
            .all(|&(lag, derivative)| lag > 0 && derivative == 0)
        || !seeds.iter().any(|row| vec_poly::is_zero(row))
        || seeds
            .iter()
            .any(|row| !vec_poly::is_zero(row) && !is_real_rooted_rational(row))
    {
        return Ok(None);
    }
    let classes = check_residue_class_subsequence_theorem(rule, first_index, seeds, &options)?;
    if classes.is_proved() {
        return Ok(classes.conclusion);
    }
    if let Some((from, _)) = discovery::tail_start(rule, first_index, seeds, &options) {
        let tail = discovery::tail_record(rule, first_index, seeds, from)?;
        let (rule, first, seeds) = tail.to_recurrence_parts()?;
        let tail = check_product_recurrence_theorem(&rule, first, &seeds, &options)?;
        if tail.is_proved() {
            return Ok(tail.conclusion);
        }
    }
    Ok(None)
}

/// The registry entry of the product theorem ([`super::registry::ROUTES`]).
pub(super) const ROUTE: RoutePolicy = RoutePolicy {
    reference: "Product induction: products of real-rooted polynomials are real-rooted; \
         RealRooted.productSequence_splits (scalar normalization and residue classes)",
    phase: DiscoveryPhase::Direct,
    after_direct_eligible: None,
    replay: ParentReplayPolicy::RecordedRoute,
    discovery_replay: DiscoveryReplayPolicy::Direct,
    check_at_supplied_seeds: product_report,
    start_policy: StartPolicy::MovesStart,
    is_preferred: any_proof,
    certificate_name: None,
    discovery: DiscoveryPolicy::Always,
    in_residue_class_menu: true,
    in_index_section_menu: true,
    in_degree_law_reversal_menu: true,
    parent: ParentPolicy {
        stage: ParentStage::Early,
        region: ParentRegion::OpenNonpositiveRay,
        record_only_if_applicable: false,
    },
};

fn product_report(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    let mut report = RealRootednessTheoremReport::new(
        RealRootednessTheorem::ProductRecurrence,
        rule,
        first_index,
        seeds,
        options,
    );
    let shape = match recurrence_shape(rule) {
        Ok(shape) => shape,
        Err(reason) => return Ok(report.not_applicable(reason)),
    };
    if shape.terms.len() != 1 {
        return Ok(report.not_applicable("expected one nonzero lag term"));
    }
    let (&(lag, derivative), multiplier) = shape.terms.iter().next().unwrap();
    if lag == 0 || derivative != 0 || seeds.len() < lag {
        return Ok(report.not_applicable(
            "product induction needs a positive lag, no derivative and at least lag seed rows",
        ));
    }
    let start = first_index + seeds.len();
    report.normal_form = format!("B(n) P_n(x) = A(n,x) P_(n-{lag})(x), n >= {start}");
    report.lean_route = Some("splitFactors".into());
    report.lean_declarations = vec!["RealRooted.productSequence_splits".into()];
    report
        .hypotheses
        .extend(denominator_hypothesis(&shape.denominator, start));
    report.hypotheses.extend(earlier_seeds_hypothesis(
        seeds,
        first_index,
        options.root_region.as_ref(),
    ));
    multiplier_hypotheses(
        multiplier,
        start,
        options.root_region.as_ref(),
        "multiplier",
        true,
        &mut report.hypotheses,
    );
    Ok(report.finish(RealRootednessTheoremConclusion {
        first_index,
        zero_residue_classes: None,
        interlacing_from: None,
        interlacing: None,
        simple_roots: false,
        root_region: options.root_region.as_ref().map(open_root_interval),
        statement: format!("every P_n, n >= {first_index}, is nonzero and real-rooted by product induction in the {lag} residue classes"),
    }))
}

fn multiplier_hypotheses(
    multiplier: &Bivariate,
    start: usize,
    region: Option<&RealInterval>,
    label: &str,
    factorize: bool,
    hypotheses: &mut Vec<LiftHypothesis>,
) {
    let sign = |suffix: &str, description: &str, polynomial: Vec<Q>, condition| {
        sign_hypothesis(
            &format!("{label}_{suffix}"),
            description,
            suffix,
            &IndexFunction::polynomial(polynomial),
            start,
            condition,
        )
    };
    if let Some((scalar, fixed)) = separate(multiplier) {
        hypotheses.push(sign(
            "scalar_nonzero",
            "the scalar factor is nonzero at every recurrence index",
            scalar,
            IndexSignCondition::Nonzero,
        ));
        let holds = match region {
            Some(region) => has_real_roots_in(&fixed, region),
            None => is_real_rooted_rational(&fixed),
        };
        hypotheses.push(exact_hypothesis(
            &format!("{label}_fixed_factor"),
            "the fixed polynomial factor is real-rooted in the requested region",
            if holds {
                Ok("exact rational root check of the fixed factor".into())
            } else {
                Err((None, "the fixed factor fails the exact root check".into()))
            },
        ));
        return;
    }
    let degree = x_degree(multiplier).unwrap_or(0);
    if degree <= 2 {
        let c = x_column(multiplier, 0);
        let b = x_column(multiplier, 1);
        let a = x_column(multiplier, 2);
        let norm = vec_poly::add(
            &vec_poly::add(&vec_poly::mul(&a, &a), &vec_poly::mul(&b, &b)),
            &vec_poly::mul(&c, &c),
        );
        hypotheses.push(sign(
            "nonzero",
            "A(n,x) is a nonzero polynomial for every recurrence index",
            norm,
            IndexSignCondition::Positive,
        ));
        if degree == 2 {
            let discriminant = vec_poly::sub(
                &vec_poly::mul(&b, &b),
                &vec_poly::scale(&vec_poly::mul(&a, &c), &Q::from_integer(4.into())),
            );
            hypotheses.push(sign(
                "discriminant",
                "b(n)^2 - 4 a(n)c(n) >= 0 for every recurrence index (degree drops allowed)",
                discriminant,
                IndexSignCondition::Nonnegative,
            ));
        }
        if let Some(region) = region {
            let leading = if degree == 2 { &a } else { &b };
            hypotheses.push(sign(
                "fixed_degree",
                "the multiplier has fixed degree throughout the requested region proof",
                leading.clone(),
                IndexSignCondition::Nonzero,
            ));
            for (suffix, bound, upper) in [
                ("upper", &region.upper, true),
                ("lower", &region.lower, false),
            ] {
                let Some(bound) = bound else {
                    continue;
                };
                // For a real-rooted quadratic, the vertex lies below U and
                // A(U)/a > 0 iff both roots lie below U; reflect for L.
                let center = if degree == 2 {
                    vec_poly::add(
                        &b,
                        &vec_poly::scale(&a, &(Q::from_integer(2.into()) * bound)),
                    )
                } else {
                    vec_poly::add(&c, &vec_poly::scale(&b, bound))
                };
                let mut center_sign = vec_poly::mul(&center, leading);
                if !upper {
                    center_sign = vec_poly::neg(&center_sign);
                }
                hypotheses.push(sign(
                    &format!("{suffix}_center"),
                    "the root (or quadratic vertex) lies strictly inside this endpoint",
                    center_sign,
                    IndexSignCondition::Positive,
                ));
                if degree == 2 {
                    let value = vec_poly::add(
                        &c,
                        &vec_poly::add(
                            &vec_poly::scale(&b, bound),
                            &vec_poly::scale(&a, &(bound * bound)),
                        ),
                    );
                    hypotheses.push(sign(
                        &format!("{suffix}_value"),
                        "A(n,endpoint)/a(n) > 0 for every recurrence index",
                        vec_poly::mul(&value, &a),
                        IndexSignCondition::Positive,
                    ));
                }
            }
        }
        return;
    }
    // The factorizer takes [x-degree][n-degree]. Restrict this optional
    // discovery step to modest degrees; the theorem has no degree bound.
    let leading = x_column(multiplier, degree);
    if factorize && degree <= 12 && multiplier.len() <= 13 && leading.len() == 1 {
        let transposed: Vec<Vec<Q>> = (0..=degree)
            .map(|j| vec_poly::scale(&x_column(multiplier, j), &leading[0].recip()))
            .collect();
        match crate::factorization::irreducible_factorization_monic_bivariate_with_budget(
            &transposed,
            256,
        ) {
            Ok(factors) => {
                for (i, factor) in factors.iter().enumerate() {
                    let width = factor.coefficients.iter().map(Vec::len).max().unwrap_or(0);
                    let coefficients: Bivariate = (0..width)
                        .map(|n| {
                            factor
                                .coefficients
                                .iter()
                                .map(|row| row.get(n).cloned().unwrap_or_else(Q::zero))
                                .collect()
                        })
                        .collect();
                    multiplier_hypotheses(
                        &coefficients,
                        start,
                        region,
                        &format!("{label}_factor_{i}"),
                        false,
                        hypotheses,
                    );
                }
                return;
            }
            Err(error) => {
                hypotheses.push(unresolved_hypothesis(
                    label,
                    "every multiplier is real-rooted",
                    format!("bounded factorization did not finish: {error}"),
                ));
                return;
            }
        }
    }
    hypotheses.push(unresolved_hypothesis(label, "every multiplier is real-rooted", "supported certificates: separable factors, linear/quadratic factors, or bounded factorization with constant leading coefficient".into()));
}

#[cfg(test)]
mod classification_tests {
    use super::*;

    #[test]
    fn zero_seed_classification_does_not_hide_invalid_live_seeds_or_multipliers() {
        let q = |n: i64| Q::from_integer(n.into());
        let rule = |coefficient| Recurrence {
            terms: vec![RecurrenceTerm::new(
                2,
                0,
                BivarPoly {
                    coeffs: vec![coefficient],
                },
            )],
            denominator: None,
            inhomogeneous: None,
        };
        let good = rule(vec![q(0), q(1)]);
        assert!(classify_product_recurrence(&good, usize::MAX, &[vec![q(1)]]).is_err());
        let seeds = vec![vec![q(1)], vec![q(0)]];
        let conclusion = classify_product_recurrence(&good, 0, &seeds)
            .unwrap()
            .unwrap();
        assert_eq!(conclusion.first_index, 0);
        assert_eq!(conclusion.zero_residue_classes.unwrap().residues, [1]);
        assert!(!is_trivially_real_rooted(&rule(vec![q(1), q(0), q(1)]), 0, &seeds).unwrap());
        assert!(!is_trivially_real_rooted(&good, 0, &[vec![q(0)], vec![q(0)]]).unwrap());
        assert!(
            !is_trivially_real_rooted(&good, 0, &[vec![q(1), q(0), q(1)], vec![q(0)]]).unwrap()
        );
        let mut derivative = good.clone();
        derivative.terms[0].deriv_order = 1;
        assert!(!is_trivially_real_rooted(&derivative, 0, &seeds).unwrap());
    }
}
