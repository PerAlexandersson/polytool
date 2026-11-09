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
/// [`check_product_recurrence_theorem`] proves it (no root region requested).
///
/// This classifies the proof of a recursively defined sequence; it is never
/// a property read off a single polynomial or a fitted recurrence.  The same
/// flag is [`RealRootednessTheoremReport::trivially_real_rooted`].
///
/// # Errors
///
/// Errors parsing or normalizing the recurrence and seeds.
pub fn is_trivially_real_rooted(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
) -> Result<bool, RecurrenceError> {
    Ok(check_product_recurrence_theorem(
        rule,
        first_index,
        seeds,
        &RealRootednessTheoremOptions::default(),
    )?
    .trivially_real_rooted)
}

/// The registry entry of the product theorem ([`super::registry::ROUTES`]).
pub(super) const ROUTE: TheoremRoute = TheoremRoute {
    theorem: RealRootednessTheorem::ProductRecurrence,
    check_at_supplied_seeds: product_report,
    start_policy: StartPolicy::MovesStart,
    is_preferred: any_proof,
    certificate_name: None,
    discovery: DiscoveryPolicy::Always {
        first_version: true,
    },
    in_residue_class_menu: true,
    in_index_section_menu: true,
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
        "Product induction: products of real-rooted polynomials are real-rooted; \
         RealRooted.productSequence_splits (scalar normalization and residue classes)",
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
