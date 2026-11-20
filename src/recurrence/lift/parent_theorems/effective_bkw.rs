//! Dominant-pair sign-grid proof route, reusing the exact BKW checker.

use super::*;
use crate::bkw::effective::{
    certify_effective_bkw_recurrence, EffectiveBkwReport, EffectiveBkwRootRegion,
};

const REFERENCE:&str="Effective dominant-pair sign-grid theorem: exact rational partial fractions, pointwise remainder domination, continuous phase crossings and degree closure";

/// Prove a stationary recurrence by an exact effective-BKW certificate.
///
/// # Errors
/// Returns [`RecurrenceError`] for malformed recurrence/seed data.
pub fn check_effective_bkw_theorem(
    rule: &Recurrence,
    first: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    check_real_rootedness_theorem(
        RealRootednessTheorem::EffectiveBkw,
        rule,
        first,
        seeds,
        options,
    )
}

pub(super) const ROUTE: RoutePolicy = RoutePolicy {
    reference: REFERENCE,
    phase: DiscoveryPhase::Direct,
    after_direct_eligible: None,
    replay: ParentReplayPolicy::RecordedRoute,
    discovery_replay: DiscoveryReplayPolicy::Direct,
    check_at_supplied_seeds: theorem_report,
    start_policy: StartPolicy::OwnStarts,
    is_preferred: any_proof,
    certificate_name: Some("effective-bkw"),
    discovery: DiscoveryPolicy::WhenStronger,
    in_residue_class_menu: false,
    in_index_section_menu: false,
    in_degree_law_reversal_menu: false,
    parent: ParentPolicy {
        stage: ParentStage::WhileUnproved,
        region: ParentRegion::OpenNonpositiveRay,
        record_only_if_applicable: true,
    },
};

fn theorem_report(
    rule: &Recurrence,
    first: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    let mut report = RealRootednessTheoremReport::new(
        RealRootednessTheorem::EffectiveBkw,
        rule,
        first,
        seeds,
        options,
    );
    // Reject unsupported shapes before any seed root isolation or bivariate
    // GF reduction. These are necessary conditions of the current suppliers.
    let index_free = |coefficient: &BivarPoly| {
        coefficient
            .coeffs
            .iter()
            .skip(1)
            .all(|row| vec_poly::is_zero(row))
    };
    if rule.inhomogeneous.is_some()
        || rule.denominator.as_ref().is_some_and(|d| !index_free(d))
        || rule
            .terms
            .iter()
            .any(|term| term.deriv_order != 0 || term.is_alternating() || !index_free(&term.coeff))
    {
        return Ok(report.not_applicable(
            "effective BKW requires a homogeneous index-free rule without derivatives or alternating signs",
        ));
    }
    // Resource guard, not a necessary degree law: coefficient degree is
    // capped at eight by this supplier, so 8*(j+1) is a deliberately generous
    // linear seed budget before exact GF construction/root checks. The
    // Newton supplier separately proves the actual slope and intercept.
    if seeds
        .iter()
        .enumerate()
        .any(|(j, row)| vec_poly::degree(row).is_some_and(|degree| degree > 8 * (j + 1)))
    {
        return Ok(
            report.not_applicable("prescribed seed degree exceeds the local-index degree budget")
        );
    }
    // The finite check is an early rejection only; the all-index claim always
    // comes from the full symbolic certificate and its own finite prefix.
    if seeds.iter().any(|s| {
        !vec_poly::is_zero(s) && !crate::is_real_rooted_bigint_coeffs(&integer_coefficients(s))
    }) {
        return Ok(report.not_applicable("a prescribed nonzero seed is not real-rooted"));
    }
    let checked = match certify_effective_bkw_recurrence(&report.recurrence) {
        Ok(c) => c,
        Err(RecurrenceError::Unsupported(reason)) => return Ok(report.not_applicable(reason)),
        Err(e) => return Err(e),
    };
    let Some(c) = checked.certificate.as_ref() else {
        return Ok(report.not_applicable(checked.reason));
    };
    let root_region = match c.root_region {
        EffectiveBkwRootRegion::WholeRealLine => None,
        EffectiveBkwRootRegion::NegativeRay => Some(RealInterval {
            lower: None,
            lower_closed: false,
            upper: Some(Q::zero()),
            upper_closed: false,
        }),
        EffectiveBkwRootRegion::PositiveRay => Some(RealInterval {
            lower: Some(Q::zero()),
            lower_closed: false,
            upper: None,
            upper_closed: false,
        }),
    };
    if let Some(requested) = &options.root_region {
        let whole = RealInterval {
            lower: None,
            lower_closed: false,
            upper: None,
            upper_closed: false,
        };
        if !requested_region_contains(requested, root_region.as_ref().unwrap_or(&whole)) {
            return Ok(report.not_applicable(
                "effective-BKW certified root region is not contained in the requested interval",
            ));
        }
    }
    report.normal_form="P_(first+n)=[z^n]N(t,z)/D(t,z); exact dominant pair, domination supplier and phase/exterior degree closure".into();
    report.lean_declarations = c.lean_declarations.clone();
    if root_region.is_some() {
        report.hypotheses.push(exact_hypothesis("effective_bkw_root_region","phase samples and exterior zeros lie on one open ray; every finite prefix row has its roots there",Ok(format!("certified {:?}, including exact finite coefficient signs",c.root_region))));
    }
    for (name, statement, detail) in [
        (
            "effective_bkw_expansion",
            "exact pair expansion from the full seed-dependent rational generating function",
            format!("{} expansion identities", c.expansion.len()),
        ),
        (
            "effective_bkw_domination",
            "strict pointwise remainder domination after the cutoff",
            c.domination.clone(),
        ),
        (
            "effective_bkw_phase_degree",
            "continuous phase crossings and exterior signs close the all-index degree budget",
            format!("{} exact sign certificates", c.signs.len()),
        ),
        (
            "effective_bkw_finite_prefix",
            "every row below the cutoff is nonzero and simple real-rooted",
            format!("{} exact rows", c.finite_rows.len()),
        ),
    ] {
        report
            .hypotheses
            .push(exact_hypothesis(name, statement, Ok(detail)));
    }
    let tail_first = checked.request.generating_function.first_index;
    report.certificate = Some(RouteCertificate::EffectiveBkw(Box::new(checked)));
    Ok(report.finish(RealRootednessTheoremConclusion{first_index:tail_first,zero_residue_classes:None,interlacing_from:None,interlacing:None,simple_roots:true,root_region:root_region.as_ref().map(open_root_interval),
        statement:format!("every P_n, n >= {tail_first}, is nonzero and has only simple real zeros by an exact effective-BKW sign-grid certificate")}))
}
pub(super) fn verify_certificate_in_report(
    c: &EffectiveBkwReport,
    report: &RealRootednessTheoremReport,
) -> Result<(), RecurrenceError> {
    let (mut rule, first, seeds) = report.recurrence.to_recurrence_parts()?;
    if let Some(normalization) = &report.left_hand_factor_normalization {
        normalization.verify()?;
        rule = normalization.normalized_recurrence()?;
    }
    let bound = if let Some(source) = &c.source_recurrence {
        let (source_rule, source_first, source_seeds) = source.to_recurrence_parts()?;
        rational_generating_function(&source_rule, source_first, &source_seeds)?
            == rational_generating_function(&rule, first, &seeds)?
    } else {
        c.request.generating_function == rational_generating_function(&rule, first, &seeds)?
    };
    if c.status != "proved" || !bound {
        return Err(RecurrenceError::VerificationFailed(
            "effective-BKW certificate belongs to another family".into(),
        ));
    }
    c.verify()
}
