//! Seeding: the sequences polytool tries, the root provers, the lifts from
//! existing nodes and the all-`n` equality of two definitions.
//!
//! - Sources: [`named_family_sources`] (the named families with defining
//!   recurrences, [`crate::recurrence::PreserverParentLibrary::named_family_records`])
//!   and [`catalog_sources`] (every bundled OEIS catalog entry whose rule is
//!   runnable; the others are listed with the reason).
//! - [`root_derivations`]: the theorem checks (three-term recurrence on the
//!   line and with the root region `(-inf, 0)`, Liu--Wang, Liu--Wang
//!   derivative lag, Ma--Wang interval, and the bounded-interval provers of
//!   the preserver lifts) and the shared parent provers (periodic product,
//!   paper theorem, root domain, then the theorem routes), each proved
//!   result a root edge.
//! - [`lift_derivations`]: preserver-image, gamma and basis-map lifts whose
//!   parents are atlas nodes (the [`LiftContext`] holds the parent
//!   libraries built from the nodes); a basis-map lift's coefficient family
//!   is tried by the root provers first, so that it can become a node.
//! - [`certify_definition_equality`]: `Q_n = P_n` for every `n` by the
//!   preserver certificate of the identity map.
//!
//! Every function here is deterministic; time limits are applied by the
//! caller (the CLI runs each item in a child process or a watched thread).

use super::analysis::{
    analyze_certificate, canonical_definition, definition_rows, parent_route_name,
    CertificateAnalysis,
};
use super::error::AtlasError;
use super::store::AtlasStore;
use super::types::*;
use crate::oeis;
use crate::real_rootedness::RealInterval;
use crate::recurrence::{
    certify_preserver_image_lift, check_liu_wang_derivative_lag_theorem, check_liu_wang_theorem,
    check_ma_wang_interval_theorem, check_three_term_recurrence_theorem, find_basis_map_lifts,
    find_gamma_lifts, find_preserver_image_lifts, interval_parent_report, prove_parent,
    BasisMapLibrary, BasisMapLiftSearchOptions, GammaLiftSearchOptions, ParentRecord,
    PreserverImageCandidate, PreserverImageOptions, PreserverImageSearchOptions, PreserverOperator,
    PreserverParentLibrary, PreserverParentLibraryOptions, RealRootednessTheoremOptions,
    RealRootednessTheoremReport, RecurrenceJson, VeroneseParentProofRoute,
    VeroneseParentProofStatus, VeroneseSectionGauge,
};
use crate::BigRational;
use num_traits::Zero;
use std::sync::Arc;

/// A derivation with the analysis computed when it was produced.
pub type AnalyzedDerivation = (AtlasDerivation, CertificateAnalysis);

/// The named families with their defining recurrences, as sources (records
/// that are not runnable are skipped).
pub fn named_family_sources() -> Vec<AtlasSequenceSpec> {
    PreserverParentLibrary::named_family_records()
        .into_iter()
        .filter_map(|(name, record)| {
            let definition = canonical_definition(&record).ok()?;
            Some(AtlasSequenceSpec {
                definition,
                source: format!("named_family:{name}"),
                name: name.clone(),
                polytool_family: Some(name),
                oeis_id: None,
                oeis_prefix_rows: 0,
            })
        })
        .collect()
}

/// Every bundled OEIS catalog entry with a runnable rule, in catalog order,
/// and the entries skipped with the reason.
pub fn catalog_sources() -> (Vec<AtlasSequenceSpec>, Vec<(String, String)>) {
    let mut sources = Vec::new();
    let mut skipped = Vec::new();
    for entry in oeis::catalog() {
        let parts = match entry.recurrence_parts() {
            Ok(parts) => parts,
            Err(error) => {
                skipped.push((entry.id.to_string(), error.to_string()));
                continue;
            }
        };
        let record = RecurrenceJson::from_recurrence_rational(
            &parts.0,
            entry.recurrence_first_index,
            &parts.1,
            None,
        );
        match canonical_definition(&record) {
            Ok(definition) => sources.push(AtlasSequenceSpec {
                definition,
                source: format!("oeis_catalog:{}", entry.id),
                name: entry.name.to_string(),
                polytool_family: None,
                oeis_id: Some(entry.id.to_string()),
                oeis_prefix_rows: entry.prefix_rows.len(),
            }),
            Err(error) => skipped.push((entry.id.to_string(), error.to_string())),
        }
    }
    (sources, skipped)
}

fn theorem_derivation(
    spec: &AtlasSequenceSpec,
    report: RealRootednessTheoremReport,
) -> Option<AnalyzedDerivation> {
    if !report.is_proved() {
        return None;
    }
    let theorem = serde_json::to_value(report.theorem)
        .ok()
        .and_then(|v| v.as_str().map(str::to_string))
        .unwrap_or_default();
    let certificate = AtlasCertificate::TheoremCheck {
        report: Box::new(report),
    };
    let analysis = analyze_certificate(&certificate).ok()?;
    Some((
        AtlasDerivation {
            conclusion: spec.clone(),
            operation: AtlasOperation::TheoremCheck { theorem },
            certificate,
        },
        analysis,
    ))
}

fn nonpositive_region_proved(derivations: &[AnalyzedDerivation]) -> bool {
    derivations.iter().any(|(_, analysis)| {
        analysis.claims.iter().any(|claim| match claim {
            AtlasClaim::RealRooted { region, .. } => {
                super::rows::region_within(region, &super::rows::nonpositive_ray())
            }
            _ => false,
        })
    })
}

/// The root derivations of a source (module documentation): every proved
/// theorem check and, when it proves the statement by a route of its own
/// or the theorem checks give no region in `(-inf, 0]`, the parent provers.
///
/// # Errors
///
/// [`AtlasError`] when the definition does not parse.
pub fn root_derivations(spec: &AtlasSequenceSpec) -> Result<Vec<AnalyzedDerivation>, AtlasError> {
    let (rule, first_index, seeds) = spec.definition.to_recurrence_parts()?;
    let mut out: Vec<AnalyzedDerivation> = Vec::new();
    let push = |derivation: Option<AnalyzedDerivation>, out: &mut Vec<AnalyzedDerivation>| {
        if let Some(derivation) = derivation {
            if !out
                .iter()
                .any(|(d, _)| d.certificate == derivation.0.certificate)
            {
                out.push(derivation);
            }
        }
    };
    let default = RealRootednessTheoremOptions::default();
    let ray = RealRootednessTheoremOptions {
        root_region: Some(RealInterval::at_most(BigRational::zero())),
        ..RealRootednessTheoremOptions::default()
    };
    let reports = [
        check_three_term_recurrence_theorem(&rule, first_index, &seeds, &default),
        check_three_term_recurrence_theorem(&rule, first_index, &seeds, &ray),
        check_liu_wang_theorem(&rule, first_index, &seeds, &default),
        check_liu_wang_derivative_lag_theorem(&rule, first_index, &seeds, &default),
        check_ma_wang_interval_theorem(&rule, first_index, &seeds, &default),
    ];
    for report in reports.into_iter().flatten() {
        let derivation = theorem_derivation(spec, report);
        push(derivation, &mut out);
    }
    if let Ok(Some(report)) = interval_parent_report(&rule, first_index, &seeds) {
        let derivation = theorem_derivation(spec, report);
        push(derivation, &mut out);
    }
    let proof = prove_parent(
        &ParentRecord {
            rule: &rule,
            seeds: &seeds,
            first_index,
        },
        first_index,
        &[],
        true,
    )?;
    let own_route = matches!(
        proof.route,
        Some(
            VeroneseParentProofRoute::PeriodicProduct
                | VeroneseParentProofRoute::PaperTheorem
                | VeroneseParentProofRoute::RootDomain
        )
    );
    if proof.status == VeroneseParentProofStatus::Proved
        && proof.route != Some(VeroneseParentProofRoute::Attestation)
        && (own_route || !nonpositive_region_proved(&out))
    {
        let route = proof.route.map(parent_route_name).unwrap_or_default();
        let certificate = AtlasCertificate::ParentProvers {
            record: spec.definition.clone(),
            try_root_domain: true,
            proof: Box::new(proof),
        };
        if let Ok(analysis) = analyze_certificate(&certificate) {
            out.push((
                AtlasDerivation {
                    conclusion: spec.clone(),
                    operation: AtlasOperation::ParentProvers {
                        route: route.to_string(),
                    },
                    certificate,
                },
                analysis,
            ));
        }
    }
    Ok(out)
}

/// The basis-map theorem derivation of a node's own rule, when the map
/// theorem is proved (the map the node defines).
///
/// # Errors
///
/// [`AtlasError`] for a definition that does not parse.
pub fn basis_map_theorem_derivation(
    spec: &AtlasSequenceSpec,
) -> Result<Option<AnalyzedDerivation>, AtlasError> {
    let (rule, first_index, seeds) = spec.definition.to_recurrence_parts()?;
    let Ok(map) = crate::recurrence::BasisMap::new(&spec.name, rule, first_index, seeds) else {
        return Ok(None);
    };
    let Ok(report) = crate::recurrence::check_basis_map_preserves_pf(
        &map,
        &crate::recurrence::BasisMapTheoremOptions::default(),
    ) else {
        return Ok(None);
    };
    let certificate = AtlasCertificate::BasisMapTheorem {
        report: Box::new(report),
    };
    let Ok(analysis) = analyze_certificate(&certificate) else {
        return Ok(None);
    };
    if analysis.claims.is_empty() {
        return Ok(None);
    }
    Ok(Some((
        AtlasDerivation {
            conclusion: spec.clone(),
            operation: AtlasOperation::BasisMapTheorem,
            certificate,
        },
        analysis,
    )))
}

/// Certify `new = existing` for every index from their (common) first index
/// on: the preserver certificate of the identity map with `existing` as the
/// parent and `new` as the input.
///
/// # Errors
///
/// [`AtlasError`] when a definition is not runnable or the identity is not
/// certified.
pub fn certify_definition_equality(
    existing: &RecurrenceJson,
    new: &RecurrenceJson,
) -> Result<AtlasCertificate, AtlasError> {
    let existing = canonical_definition(existing)?;
    let new = canonical_definition(new)?;
    let (rule, first_index, _) = new.to_recurrence_parts()?;
    let rows = definition_rows(&new, 16)?;
    let candidate = PreserverImageCandidate {
        parent_name: "existing definition".into(),
        parent: existing,
        index_shift: 0,
        gauge: VeroneseSectionGauge::None,
        operator: PreserverOperator::Affine {
            slope: "1".into(),
            intercept: "0".into(),
        },
    };
    let certification = certify_preserver_image_lift(
        &rule,
        &rows,
        first_index,
        &candidate,
        &PreserverImageOptions {
            try_root_domain: false,
            ..PreserverImageOptions::default()
        },
    )?;
    let certificate = AtlasCertificate::DefinitionEquality {
        certification: Box::new(certification),
    };
    analyze_certificate(&certificate)?;
    Ok(certificate)
}

/// The parent libraries of the lifts, built from the atlas nodes.
#[derive(Debug, Clone)]
pub struct LiftContext {
    /// Preserver parents: the first definition of every node.
    pub preserver_library: Arc<PreserverParentLibrary>,
    /// The admitted basis maps of the OEIS catalog.
    pub basis_library: Arc<BasisMapLibrary>,
}

/// Build the lift context from the nodes of `store` (the parent provers run
/// once per node; nodes they do not prove are skipped by the library).
///
/// # Errors
///
/// [`AtlasError`] from an internal failure of the provers.
pub fn lift_context(store: &AtlasStore) -> Result<LiftContext, AtlasError> {
    let options = PreserverParentLibraryOptions {
        include_named_families: false,
        ..PreserverParentLibraryOptions::default()
    };
    let mut library = PreserverParentLibrary::default();
    for node in store.nodes() {
        let has_region = node
            .properties
            .iter()
            .any(|p| matches!(p.claim, AtlasClaim::RealRooted { .. }));
        if has_region {
            library.add(node.id.clone(), &node.definitions[0].recurrence, &options)?;
        }
    }
    Ok(LiftContext {
        preserver_library: Arc::new(library),
        basis_library: Arc::new(BasisMapLibrary::oeis_catalog()?),
    })
}

/// The lift derivations of a source (module documentation), with the root
/// derivations of the derived parents they need first.
///
/// # Errors
///
/// [`AtlasError`] when the definition does not parse.
pub fn lift_derivations(
    spec: &AtlasSequenceSpec,
    context: &LiftContext,
) -> Result<Vec<AnalyzedDerivation>, AtlasError> {
    let (rule, first_index, seeds) = spec.definition.to_recurrence_parts()?;
    let rows = definition_rows(&spec.definition, 16)?;
    let mut out = Vec::new();
    let preserver = find_preserver_image_lifts(
        &rule,
        &rows,
        first_index,
        &PreserverImageSearchOptions {
            library: Some(context.preserver_library.clone()),
            ..PreserverImageSearchOptions::default()
        },
    );
    if let Ok(report) = preserver {
        for certification in report.lifts {
            let operation = AtlasOperation::PreserverImage {
                operator: serde_json::to_string(&certification.candidate.operator)
                    .unwrap_or_default(),
                index_shift: certification.candidate.index_shift,
            };
            let certificate = AtlasCertificate::PreserverImage {
                certification: Box::new(certification),
            };
            if let Ok(analysis) = analyze_certificate(&certificate) {
                out.push((
                    AtlasDerivation {
                        conclusion: spec.clone(),
                        operation,
                        certificate,
                    },
                    analysis,
                ));
            }
        }
    }
    if let Ok(report) = find_gamma_lifts(
        &rule,
        &seeds,
        first_index,
        &GammaLiftSearchOptions::default(),
    ) {
        for certification in report.lifts {
            let operation = AtlasOperation::GammaPolynomial {
                scale: certification.candidate.scale.clone(),
                index_shift: certification.candidate.index_shift,
            };
            let certificate = AtlasCertificate::GammaPolynomial {
                certification: Box::new(certification),
            };
            if let Ok(analysis) = analyze_certificate(&certificate) {
                out.push((
                    AtlasDerivation {
                        conclusion: spec.clone(),
                        operation,
                        certificate,
                    },
                    analysis,
                ));
            }
        }
    }
    let basis = find_basis_map_lifts(
        &rule,
        &rows,
        first_index,
        &BasisMapLiftSearchOptions {
            library: Some(context.basis_library.clone()),
            ..BasisMapLiftSearchOptions::default()
        },
    );
    if let Ok(report) = basis {
        for certification in report.lifts {
            let operation = AtlasOperation::BasisMap {
                map: certification.candidate.basis.name().to_string(),
            };
            let coefficient = certification.candidate.coefficient_family.clone();
            let target = certification
                .coefficient_roots
                .as_ref()
                .map(|proof| proof.target_interval.clone());
            let certificate = AtlasCertificate::BasisMap {
                certification: Box::new(certification),
            };
            let Ok(analysis) = analyze_certificate(&certificate) else {
                continue;
            };
            if let Ok(definition) = canonical_definition(&coefficient) {
                let parent = AtlasSequenceSpec {
                    definition,
                    source: format!("derived:basis_map_coefficients_of:{}", spec.source),
                    name: format!("basis-map coefficients of {}", spec.name),
                    polytool_family: None,
                    oeis_id: None,
                    oeis_prefix_rows: 0,
                };
                if let Some(Ok(target)) = target.as_ref().map(|t| t.at(0)) {
                    if let Some(upper) = &target.upper {
                        out.extend(bound_derivations(&parent, upper, true));
                    }
                    if let Some(lower) = &target.lower {
                        out.extend(bound_derivations(&parent, lower, false));
                    }
                }
                out.extend(root_derivations(&parent)?);
            }
            out.push((
                AtlasDerivation {
                    conclusion: spec.clone(),
                    operation,
                    certificate,
                },
                analysis,
            ));
        }
    }
    Ok(out)
}

/// The derivations that bound the roots of a family `R` by `bound`: the
/// shifted family `S_n(v) = R_n(v + bound)` (upper bound) or
/// `S_n(v) = R_n(bound - v)` (lower bound), computed exactly as a preserver
/// image of `R` and proved by the root provers (roots `<= 0`), then `R` as
/// the preserver image `S_n(x - bound)` or `S_n(bound - x)` of the proved
/// `S`, which places every root of `R` in `(-inf, bound]` or
/// `[bound, inf)`.  This mirrors the coefficient root proof of the
/// basis-map lifts with atlas edges.  Failures give no derivation.
fn bound_derivations(
    family: &AtlasSequenceSpec,
    bound: &BigRational,
    upper: bool,
) -> Vec<AnalyzedDerivation> {
    let text = crate::recurrence::format_rational_coeff(bound);
    let negated = crate::recurrence::format_rational_coeff(&-bound.clone());
    let (forward, backward) = if upper {
        (("1".to_string(), text.clone()), ("1".to_string(), negated))
    } else {
        (
            ("-1".to_string(), text.clone()),
            ("-1".to_string(), text.clone()),
        )
    };
    let attempt = || -> Result<Vec<AnalyzedDerivation>, AtlasError> {
        let n0 = family.definition.first_index;
        let shift = PreserverImageCandidate {
            parent_name: family.name.clone(),
            parent: family.definition.clone(),
            index_shift: 0,
            gauge: VeroneseSectionGauge::None,
            operator: PreserverOperator::Affine {
                slope: forward.0.clone(),
                intercept: forward.1.clone(),
            },
        };
        let image =
            shift.image_recurrence(n0, &crate::recurrence::ImageRecurrenceOptions::default())?;
        let shifted = AtlasSequenceSpec {
            definition: canonical_definition(&image.rule)?,
            source: format!("derived:root_bound_shift:{}", family.source),
            name: format!(
                "{} shifted by {} {text}",
                family.name,
                if upper { "x ->  x +" } else { "x -> -x +" }
            ),
            polytool_family: None,
            oeis_id: None,
            oeis_prefix_rows: 0,
        };
        let mut out = root_derivations(&shifted)?;
        if out.is_empty() {
            return Ok(out);
        }
        let back = PreserverImageCandidate {
            parent_name: shifted.name.clone(),
            parent: shifted.definition.clone(),
            index_shift: 0,
            gauge: VeroneseSectionGauge::None,
            operator: PreserverOperator::Affine {
                slope: backward.0.clone(),
                intercept: backward.1.clone(),
            },
        };
        let (rule, first_index, _) = family.definition.to_recurrence_parts()?;
        let rows = definition_rows(&family.definition, 16)?;
        let certification = certify_preserver_image_lift(
            &rule,
            &rows,
            first_index,
            &back,
            &PreserverImageOptions {
                try_root_domain: false,
                ..PreserverImageOptions::default()
            },
        )?;
        let operation = AtlasOperation::PreserverImage {
            operator: serde_json::to_string(&back.operator).unwrap_or_default(),
            index_shift: 0,
        };
        let certificate = AtlasCertificate::PreserverImage {
            certification: Box::new(certification),
        };
        let analysis = analyze_certificate(&certificate)?;
        out.push((
            AtlasDerivation {
                conclusion: family.clone(),
                operation,
                certificate,
            },
            analysis,
        ));
        Ok(out)
    };
    attempt().unwrap_or_default()
}
