//! One forward round of budgeted operations on existing nodes
//! ([`plan_forward_round`], [`forward_derivations`]).
//!
//! # Operations
//!
//! - **Preserver images** `T(P_n)` for a small menu of operators that do not
//!   depend on `n`: `P'`, `P + P'` (Laguerre--Pólya operators `phi(D)`),
//!   `(theta + 1) P` and `(theta + 1)^2 P` (multiplier sequences of
//!   Laguerre's theorem), `P(x + 1)`, and, for an affine degree pattern, the
//!   reversal `x^(d_n) P(1/x)` and the Möbius map `(x + 1)^(d_n) P(x/(x + 1))`.
//!   The image's runnable rule is computed exactly
//!   ([`crate::recurrence::image_recurrence_of_operator`]), not fitted; the
//!   lift certificate ([`crate::recurrence::certify_preserver_image_lift`])
//!   proves the identity and the transfer.
//! - **Basis maps** `T_P(R_n)` for the first admitted maps of the OEIS
//!   catalog ([`crate::recurrence::apply_basis_map_to_family`]) applied to
//!   nodes with derivative-free rules.
//! - **Gamma maps**: the gamma polynomials `Gamma_n` of a node whose stored
//!   rows are palindromic at `d_n = deg P_n` (affine) with nonzero constant
//!   terms; the image rule is fitted to the gamma rows and then certified by
//!   [`crate::recurrence::certify_gamma_lift`] (a fitted rule is only a
//!   candidate: the certificate proves it for every index).
//! - **Veronese sections** `V_(r,j)(P_n)`, `r = 2, 3`: fitted section rule,
//!   certified by [`crate::recurrence::certify_veronese_section_lift`]; the
//!   edge records the section rule, and the audit recomputes the whole
//!   certification from it.
//! - **Obreschkoff combinations** `P_n + P_(n+1)` of nodes with interlacing
//!   consecutive rows and a first-order rule ([`super::obreschkoff`]).
//!
//! # Budget and priorities
//!
//! Sequence-by-map pairs grow quadratically, so [`ForwardBudget`] caps the
//! tasks per operation and in total.  Nodes are taken in increasing order of
//! a size key (the bit length of the largest coefficient of the stored row
//! with index `first_index + 8`, then the rule order, then the id), so
//! small nodes come first; maps are taken in library order (small maps
//! first).  Every task is independent and deterministic; time limits are
//! applied by the caller.  Images whose stored rows are not all integers
//! are dropped.

use super::analysis::{analyze_certificate, canonical_definition, definition_rows};
use super::error::AtlasError;
use super::integrate::definition_integer_rows;
use super::obreschkoff::certify_obreschkoff_combination;
use super::rows::{nonpositive_ray, region_within, rows_from_strings};
use super::seed::AnalyzedDerivation;
use super::store::AtlasStore;
use super::types::*;
use crate::recurrence::{
    apply_basis_map_to_family, certify_gamma_lift, certify_preserver_image_lift,
    certify_veronese_section_lift, find_recurrence_adaptive_rational_with_budget,
    AdaptiveSearchBudget, AdaptiveSearchOptions, AdaptiveSearchOutcome, BasisMap,
    BasisMapImageOptions, BasisMapLibrary, DegreeSchedule, GammaLiftCandidate, GammaLiftOptions,
    GammaLiftParent, ImageRecurrenceOptions, IndexRationalFunction, PreserverImageCandidate,
    PreserverImageOptions, PreserverOperator, Recurrence, RecurrenceJson, VeroneseSectionCandidate,
    VeroneseSectionCertifyOptions, VeroneseSectionGauge,
};
use crate::{gamma_coefficients_bigint_coeffs, BigRational};
use num_bigint::BigInt;
use num_traits::{Signed, Zero};
use serde::{Deserialize, Serialize};

/// One forward operation on one node (or one map and one node).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "task", rename_all = "snake_case", deny_unknown_fields)]
pub enum ForwardTask {
    /// `Q_n = T(P_n)`.
    PreserverImage {
        /// The parent node's id.
        node: String,
        /// The parent (its first definition).
        parent: AtlasSequenceSpec,
        /// The operator.
        operator: PreserverOperator,
        /// A short label of the operator.
        label: String,
    },
    /// `Q_n = T_P(R_n)`.
    BasisMap {
        /// The coefficient node's id.
        node: String,
        /// The map.
        map: BasisMap,
        /// The coefficient family `R` (its first definition).
        coefficients: AtlasSequenceSpec,
    },
    /// `Q_n = Gamma_n`.
    GammaPolynomial {
        /// The parent node's id.
        node: String,
        /// The parent.
        parent: AtlasSequenceSpec,
        /// The symmetry indices `d_n = deg P_n`.
        schedule: DegreeSchedule,
    },
    /// `Q_n = V_(r,j)(P_n)`.
    VeroneseSection {
        /// The parent node's id.
        node: String,
        /// The parent.
        parent: AtlasSequenceSpec,
        /// The modulus `r`.
        modulus: usize,
        /// The residue `j`.
        residue: usize,
    },
    /// `Q_n = P_n + P_(n+1)`, `n >= first_index`.
    Obreschkoff {
        /// The parent node's id.
        node: String,
        /// The parent.
        parent: AtlasSequenceSpec,
        /// The first index of `Q`.
        first_index: usize,
    },
}

impl ForwardTask {
    /// A short label of the operation.
    pub fn operation(&self) -> &'static str {
        match self {
            Self::PreserverImage { .. } => "preserver_image",
            Self::BasisMap { .. } => "basis_map",
            Self::GammaPolynomial { .. } => "gamma_polynomial",
            Self::VeroneseSection { .. } => "veronese_section",
            Self::Obreschkoff { .. } => "obreschkoff_combination",
        }
    }

    /// The id of the node the task transforms.
    pub fn node(&self) -> &str {
        match self {
            Self::PreserverImage { node, .. }
            | Self::BasisMap { node, .. }
            | Self::GammaPolynomial { node, .. }
            | Self::VeroneseSection { node, .. }
            | Self::Obreschkoff { node, .. } => node,
        }
    }

    /// A one-line description (node and parameters).
    pub fn label(&self) -> String {
        match self {
            Self::PreserverImage { node, label, .. } => format!("{node} {label}"),
            Self::BasisMap { node, map, .. } => format!("{node} basis map {}", map.name()),
            Self::GammaPolynomial { node, .. } => format!("{node} gamma"),
            Self::VeroneseSection {
                node,
                modulus,
                residue,
                ..
            } => format!("{node} section r={modulus} j={residue}"),
            Self::Obreschkoff {
                node, first_index, ..
            } => format!("{node} P_n + P_(n+1) from {first_index}"),
        }
    }
}

/// Caps of one forward round (module documentation).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ForwardBudget {
    /// Largest number of tasks in total.  Default 3000.
    pub max_tasks: usize,
    /// Largest number of preserver-image tasks.  Default 1500.
    pub max_preserver_tasks: usize,
    /// Largest number of basis maps tried.  Default 6.
    pub max_basis_maps: usize,
    /// Largest number of basis-map tasks.  Default 600.
    pub max_basis_map_tasks: usize,
    /// Largest number of gamma tasks.  Default 300.
    pub max_gamma_tasks: usize,
    /// Largest number of Veronese tasks.  Default 600.
    pub max_veronese_tasks: usize,
    /// Largest Veronese modulus.  Default 3.
    pub max_veronese_modulus: usize,
    /// Largest number of Obreschkoff tasks.  Default 300.
    pub max_obreschkoff_tasks: usize,
    /// Largest bit length of the size key of a node (module documentation).
    /// Default 64.
    pub max_size_bits: u64,
}

impl Default for ForwardBudget {
    fn default() -> Self {
        Self {
            max_tasks: 3000,
            max_preserver_tasks: 1500,
            max_basis_maps: 6,
            max_basis_map_tasks: 600,
            max_gamma_tasks: 300,
            max_veronese_tasks: 600,
            max_veronese_modulus: 3,
            max_obreschkoff_tasks: 300,
            max_size_bits: 64,
        }
    }
}

fn integer_function(value: i64) -> IndexRationalFunction {
    IndexRationalFunction::constant(&BigRational::from_integer(value.into()))
}

/// The preserver menu for a node: `(label, operator, needs a region in
/// (-inf, 0])`.
fn preserver_menu(node: &AtlasNode) -> Vec<(String, PreserverOperator, bool)> {
    let mut menu = vec![
        (
            "derivative".to_string(),
            PreserverOperator::DerivativePolynomial {
                coefficients: vec![integer_function(0), integer_function(1)],
            },
            false,
        ),
        (
            "one_plus_derivative".to_string(),
            PreserverOperator::DerivativePolynomial {
                coefficients: vec![integer_function(1), integer_function(1)],
            },
            false,
        ),
        (
            "theta_plus_one".to_string(),
            PreserverOperator::EulerPolynomial {
                coefficients: vec![integer_function(1), integer_function(1)],
            },
            true,
        ),
        (
            "theta_plus_one_squared".to_string(),
            PreserverOperator::EulerPolynomial {
                coefficients: vec![
                    integer_function(1),
                    integer_function(2),
                    integer_function(1),
                ],
            },
            true,
        ),
        (
            "shift_x_plus_one".to_string(),
            PreserverOperator::Affine {
                slope: "1".into(),
                intercept: "1".into(),
            },
            false,
        ),
    ];
    if let DegreePattern::Affine { slope, intercept } = node.observations.degree_pattern {
        let first = slope * node.first_index as i64 + intercept;
        if first >= 0 && slope >= 0 {
            let schedule = DegreeSchedule::new(node.first_index, first as usize, slope);
            menu.push((
                "reversal".to_string(),
                PreserverOperator::Reversal {
                    schedule: schedule.clone(),
                },
                true,
            ));
            menu.push((
                "mobius_x_over_x_plus_one".to_string(),
                PreserverOperator::Mobius {
                    a: "1".into(),
                    b: "0".into(),
                    c: "1".into(),
                    e: "1".into(),
                    schedule,
                },
                true,
            ));
        }
    }
    menu
}

fn spec_of(node: &AtlasNode) -> AtlasSequenceSpec {
    AtlasSequenceSpec {
        definition: node.definitions[0].recurrence.clone(),
        source: node.definitions[0].source.clone(),
        name: node.description.name.clone(),
        polytool_family: node.description.polytool_family.clone(),
        oeis_id: None,
        oeis_prefix_rows: 0,
    }
}

fn size_key(node: &AtlasNode) -> (u64, usize, String) {
    let rows = rows_from_strings(&node.rows).unwrap_or_default();
    let row = rows.get(8).or(rows.last());
    let bits = row
        .map(|row| row.iter().map(|c| c.abs().bits()).max().unwrap_or(0))
        .unwrap_or(0);
    let order = node.definitions[0]
        .recurrence
        .recurrence
        .terms
        .iter()
        .map(|t| t.offset)
        .max()
        .unwrap_or(0);
    (bits, order, node.id.clone())
}

fn has_region_within(node: &AtlasNode, region: &crate::recurrence::RootInterval) -> bool {
    node.properties.iter().any(|p| match &p.claim {
        AtlasClaim::RealRooted { region: r, .. } => region_within(r, region),
        _ => false,
    })
}

fn is_real_rooted(node: &AtlasNode) -> bool {
    node.properties
        .iter()
        .any(|p| matches!(p.claim, AtlasClaim::RealRooted { .. }))
}

fn is_derivative_free(node: &AtlasNode) -> bool {
    node.definitions[0]
        .recurrence
        .recurrence
        .terms
        .iter()
        .all(|t| t.deriv_order == 0)
}

/// Plan the tasks of one forward round over the nodes of `store` (module
/// documentation), in a deterministic order.
///
/// # Errors
///
/// [`AtlasError`] when the basis-map library cannot be built.
pub fn plan_forward_round(
    store: &AtlasStore,
    budget: &ForwardBudget,
) -> Result<Vec<ForwardTask>, AtlasError> {
    let mut nodes: Vec<&AtlasNode> = store.nodes().iter().filter(|n| is_real_rooted(n)).collect();
    nodes.sort_by_cached_key(|n| size_key(n));
    nodes.retain(|n| size_key(n).0 <= budget.max_size_bits);
    let nonpositive = nonpositive_ray();
    let mut tasks = Vec::new();
    let mut preserver = 0;
    for node in &nodes {
        for (label, operator, needs_ray) in preserver_menu(node) {
            if preserver >= budget.max_preserver_tasks {
                break;
            }
            if needs_ray && !has_region_within(node, &nonpositive) {
                continue;
            }
            tasks.push(ForwardTask::PreserverImage {
                node: node.id.clone(),
                parent: spec_of(node),
                operator,
                label,
            });
            preserver += 1;
        }
    }
    let library = BasisMapLibrary::oeis_catalog()?;
    let mut basis = 0;
    for entry in library.entries().iter().take(budget.max_basis_maps) {
        for node in &nodes {
            if basis >= budget.max_basis_map_tasks {
                break;
            }
            if !is_derivative_free(node) || !has_region_within(node, &nonpositive) {
                continue;
            }
            tasks.push(ForwardTask::BasisMap {
                node: node.id.clone(),
                map: entry.map.clone(),
                coefficients: spec_of(node),
            });
            basis += 1;
        }
    }
    let mut gamma = 0;
    for node in &nodes {
        if gamma >= budget.max_gamma_tasks {
            break;
        }
        let rows = rows_from_strings(&node.rows).unwrap_or_default();
        let constant_terms = rows.iter().all(|r| r.first().is_some_and(|c| !c.is_zero()));
        let DegreePattern::Affine { slope, intercept } = node.observations.degree_pattern else {
            continue;
        };
        let first = slope * node.first_index as i64 + intercept;
        if !node.observations.palindromic
            || !constant_terms
            || slope < 1
            || first < 0
            || !has_region_within(node, &nonpositive)
        {
            continue;
        }
        tasks.push(ForwardTask::GammaPolynomial {
            node: node.id.clone(),
            parent: spec_of(node),
            schedule: DegreeSchedule::new(node.first_index, first as usize, slope),
        });
        gamma += 1;
    }
    let mut veronese = 0;
    for modulus in 2..=budget.max_veronese_modulus {
        for node in &nodes {
            if !has_region_within(node, &nonpositive) {
                continue;
            }
            for residue in 0..modulus {
                if veronese >= budget.max_veronese_tasks {
                    break;
                }
                tasks.push(ForwardTask::VeroneseSection {
                    node: node.id.clone(),
                    parent: spec_of(node),
                    modulus,
                    residue,
                });
                veronese += 1;
            }
        }
    }
    let mut obreschkoff = 0;
    for node in &nodes {
        if obreschkoff >= budget.max_obreschkoff_tasks {
            break;
        }
        let interlacing_from = node
            .properties
            .iter()
            .filter_map(|p| match p.claim {
                AtlasClaim::Interlacing { from_index, .. } => Some(from_index),
                _ => None,
            })
            .min();
        let rooted_from = node
            .properties
            .iter()
            .filter_map(|p| match p.claim {
                AtlasClaim::RealRooted { from_index, .. } => Some(from_index),
                _ => None,
            })
            .min();
        let (Some(interlacing_from), Some(rooted_from)) = (interlacing_from, rooted_from) else {
            continue;
        };
        let Ok((_, start)) =
            super::obreschkoff::obreschkoff_operator(&node.definitions[0].recurrence)
        else {
            continue;
        };
        let first_index = start
            .max(interlacing_from.saturating_sub(1))
            .max(rooted_from);
        tasks.push(ForwardTask::Obreschkoff {
            node: node.id.clone(),
            parent: spec_of(node),
            first_index,
        });
        obreschkoff += 1;
    }
    tasks.truncate(budget.max_tasks);
    Ok(tasks)
}

fn derived_spec(
    definition: RecurrenceJson,
    operation: &str,
    parent: &AtlasSequenceSpec,
) -> AtlasSequenceSpec {
    AtlasSequenceSpec {
        definition,
        source: format!("derived:{operation}"),
        name: format!("{operation} of {}", parent.name),
        polytool_family: None,
        oeis_id: None,
        oeis_prefix_rows: 0,
    }
}

fn require_integer_rows(definition: &RecurrenceJson) -> Result<(), AtlasError> {
    match definition_integer_rows(definition)? {
        Some(_) => Ok(()),
        None => Err(AtlasError::Unsupported(
            "the image rows are not all integers".into(),
        )),
    }
}

/// Fit a runnable rule to `rows` (indexed from `first_index`), without
/// derivatives, up to order `max_order`.
fn fit_rule(rows: &[Vec<BigRational>], first_index: usize, max_order: usize) -> Option<Recurrence> {
    let search = AdaptiveSearchOptions {
        max_rec_len: max_order,
        max_var_deg: 4,
        max_idx_deg: 3,
        max_diff_deg: 1,
        max_denom_var_deg: 0,
        max_denom_idx_deg: 3,
        ..AdaptiveSearchOptions::default()
    };
    match find_recurrence_adaptive_rational_with_budget(
        rows,
        &search,
        AdaptiveSearchBudget::limited(600),
    ) {
        AdaptiveSearchOutcome::Found(found) => {
            Some(found.recurrence.shift_index(1 - first_index as i64))
        }
        _ => None,
    }
}

fn record_from_rule(
    rule: &Recurrence,
    first_index: usize,
    rows: &[Vec<BigRational>],
) -> Result<RecurrenceJson, AtlasError> {
    let count = rule
        .generation_initial_count(first_index, rows.len())
        .max(1);
    let record = RecurrenceJson::from_recurrence_rational(
        rule,
        first_index,
        &rows[..count.min(rows.len())],
        None,
    );
    canonical_definition(&record)
}

fn analyzed(
    conclusion: AtlasSequenceSpec,
    operation: AtlasOperation,
    certificate: AtlasCertificate,
) -> Result<AnalyzedDerivation, AtlasError> {
    let analysis = analyze_certificate(&certificate)?;
    Ok((
        AtlasDerivation {
            conclusion,
            operation,
            certificate,
        },
        analysis,
    ))
}

/// The derivations of one task (module documentation).
///
/// # Errors
///
/// [`AtlasError`] when the operation does not apply, its rule is not found
/// within the budget, the image rows are not integral or the certificate
/// proves nothing usable.
pub fn forward_derivations(task: &ForwardTask) -> Result<Vec<AnalyzedDerivation>, AtlasError> {
    match task {
        ForwardTask::PreserverImage {
            parent,
            operator,
            label,
            ..
        } => {
            let candidate = PreserverImageCandidate {
                parent_name: parent.name.clone(),
                parent: parent.definition.clone(),
                index_shift: 0,
                gauge: VeroneseSectionGauge::None,
                operator: operator.clone(),
            };
            let n0 = parent.definition.first_index;
            let image = candidate.image_recurrence(n0, &ImageRecurrenceOptions::default())?;
            let definition = canonical_definition(&image.rule)?;
            require_integer_rows(&definition)?;
            let (rule, first_index, _) = definition.to_recurrence_parts()?;
            let rows = definition_rows(&definition, 16)?;
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
            let operation = AtlasOperation::PreserverImage {
                operator: serde_json::to_string(operator).unwrap_or_default(),
                index_shift: 0,
            };
            Ok(vec![analyzed(
                derived_spec(definition, label, parent),
                operation,
                AtlasCertificate::PreserverImage {
                    certification: Box::new(certification),
                },
            )?])
        }
        ForwardTask::BasisMap {
            map, coefficients, ..
        } => {
            let image = apply_basis_map_to_family(
                map,
                &coefficients.definition,
                &BasisMapImageOptions::default(),
            )?;
            let definition = canonical_definition(&image.image)?;
            require_integer_rows(&definition)?;
            let operation = AtlasOperation::BasisMap {
                map: map.name().to_string(),
            };
            Ok(vec![analyzed(
                derived_spec(
                    definition,
                    &format!("basis_map_{}", map.name()),
                    coefficients,
                ),
                operation,
                AtlasCertificate::BasisMap {
                    certification: Box::new(image.certification),
                },
            )?])
        }
        ForwardTask::GammaPolynomial {
            parent, schedule, ..
        } => {
            let n0 = parent.definition.first_index;
            let parent_rows = definition_integer_rows_long(&parent.definition, 24)?;
            let mut gamma_rows = Vec::new();
            for row in &parent_rows {
                let gamma = gamma_coefficients_bigint_coeffs(row).ok_or_else(|| {
                    AtlasError::Unsupported("a parent row is not palindromic".into())
                })?;
                gamma_rows.push(trimmed_rational(&gamma));
            }
            let rule = fit_rule(&gamma_rows, n0, 3)
                .ok_or_else(|| AtlasError::Unsupported("no rule fits the gamma rows".into()))?;
            let definition = record_from_rule(&rule, n0, &gamma_rows)?;
            require_integer_rows(&definition)?;
            let (rule, first_index, _) = definition.to_recurrence_parts()?;
            let rows = definition_rows(&definition, 16)?;
            let candidate = GammaLiftCandidate {
                parent: GammaLiftParent {
                    name: parent.name.clone(),
                    recurrence: parent.definition.clone(),
                    symmetry_schedule: schedule.clone(),
                },
                index_shift: 0,
                scale: "1".into(),
                gauge: VeroneseSectionGauge::None,
            };
            let certification = certify_gamma_lift(
                &rule,
                &rows,
                first_index,
                &candidate,
                &GammaLiftOptions::default(),
            )?;
            Ok(vec![analyzed(
                derived_spec(definition, "gamma_polynomial", parent),
                AtlasOperation::GammaPolynomial {
                    scale: "1".into(),
                    index_shift: 0,
                },
                AtlasCertificate::GammaPolynomial {
                    certification: Box::new(certification),
                },
            )?])
        }
        ForwardTask::VeroneseSection {
            parent,
            modulus,
            residue,
            ..
        } => {
            let n0 = parent.definition.first_index;
            let parent_rows = definition_integer_rows_long(&parent.definition, 24)?;
            let section_rows: Vec<Vec<BigRational>> = parent_rows
                .iter()
                .map(|row| {
                    let section: Vec<BigInt> = row
                        .iter()
                        .skip(*residue)
                        .step_by(*modulus)
                        .cloned()
                        .collect();
                    trimmed_rational(&section)
                })
                .collect();
            if section_rows.iter().all(|r| r.is_empty()) {
                return Err(AtlasError::Unsupported("every section row is zero".into()));
            }
            let rule = fit_rule(&section_rows, n0, 2 * modulus)
                .ok_or_else(|| AtlasError::Unsupported("no rule fits the section rows".into()))?;
            let definition = record_from_rule(&rule, n0, &section_rows)?;
            require_integer_rows(&definition)?;
            let certification =
                veronese_certification(&definition, &parent.definition, *modulus, *residue)?;
            Ok(vec![analyzed(
                derived_spec(
                    definition.clone(),
                    &format!("veronese_section_{modulus}_{residue}"),
                    parent,
                ),
                AtlasOperation::VeroneseSection {
                    modulus: *modulus,
                    residue: *residue,
                    index_shift: 0,
                },
                AtlasCertificate::VeroneseSection {
                    input: definition,
                    certification: Box::new(certification),
                },
            )?])
        }
        ForwardTask::Obreschkoff {
            parent,
            first_index,
            ..
        } => {
            let certificate = certify_obreschkoff_combination(
                &parent.definition,
                *first_index,
                &ImageRecurrenceOptions::default(),
            )?;
            let definition = canonical_definition(&preserver_input_of(&certificate)?)?;
            require_integer_rows(&definition)?;
            Ok(vec![analyzed(
                derived_spec(definition, "obreschkoff_combination", parent),
                AtlasOperation::ObreschkoffCombination,
                AtlasCertificate::ObreschkoffCombination {
                    certificate: Box::new(certificate),
                },
            )?])
        }
    }
}

fn preserver_input_of(
    certificate: &super::obreschkoff::ObreschkoffCertificate,
) -> Result<RecurrenceJson, AtlasError> {
    Ok(super::analysis::preserver_input(&certificate.identity)?.clone())
}

/// The Veronese certification of `input = V_(r,j)(parent)` (index shift 0,
/// no gauge), with the default options (the audit recomputes it).
pub(crate) fn veronese_certification(
    input: &RecurrenceJson,
    parent: &RecurrenceJson,
    modulus: usize,
    residue: usize,
) -> Result<crate::recurrence::VeroneseSectionCertification, AtlasError> {
    let (rule, first_index, _) = input.to_recurrence_parts()?;
    let rows = definition_rows(input, 16)?;
    Ok(certify_veronese_section_lift(
        &rule,
        &rows,
        first_index,
        &VeroneseSectionCandidate {
            modulus,
            residue,
            parent: parent.clone(),
            index_shift: 0,
            gauge: VeroneseSectionGauge::None,
        },
        &VeroneseSectionCertifyOptions {
            try_root_domain: false,
            ..VeroneseSectionCertifyOptions::default()
        },
    )?)
}

fn trimmed_rational(row: &[BigInt]) -> Vec<BigRational> {
    let mut out: Vec<BigRational> = row.iter().cloned().map(BigRational::from_integer).collect();
    while out.last().is_some_and(Zero::is_zero) {
        out.pop();
    }
    out
}

/// The first `count` rows of a definition, which must be integral.
fn definition_integer_rows_long(
    definition: &RecurrenceJson,
    count: usize,
) -> Result<Vec<Vec<BigInt>>, AtlasError> {
    super::rows::integer_rows(&definition_rows(definition, count)?)
        .ok_or_else(|| AtlasError::Unsupported("the parent rows are not integral".into()))
}
