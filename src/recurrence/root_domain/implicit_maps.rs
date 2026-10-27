//! Root-domain information for implicit transports: diagonal coefficient
//! maps and derivative antecedents (`polytool.recurrence-implicit.v1`).
//!
//! [`analyze_implicit_root_domain`] attaches an [`ImplicitRootDomainReport`]
//! to an [`ImplicitTransformReport`]; the algebraic report never depends on
//! it.  Theorem adapters, each with its direction and exact hypotheses:
//!
//! - `polya_schur_multiplier_sequence` (`diagonal`): a multiplier sequence
//!   `(lambda_k)` maps every real-rooted polynomial to a real-rooted
//!   polynomial or zero.  Certified sources of multiplier sequences, for an
//!   index-free weight `lambda_k = scale (k!)^p prod_j f_j(k)^(e_j)`: Laguerre's
//!   theorem (`f(k)` for a polynomial `f` with only real roots, all `<= 0`;
//!   `f(theta)` is a product of `theta + a = x^(1-a) D x^a`, `a >= 0`), the
//!   sequence `1/k!` (Pólya--Schur: `sum x^k/(k!)^2` is in the Laguerre--Pólya
//!   class with nonnegative coefficients), termwise products (compositions)
//!   and nonzero scalings.  So `P_n` real-rooted implies `Q_n = T P_n`
//!   real-rooted when `p <= 0`, every `e_j > 0` and every `f_j` is
//!   real-rooted with roots `<= 0` (for example `factorial_power`,
//!   `binomial`); the inverse statement holds when `1/lambda` qualifies
//!   (`p >= 0`, every `e_j < 0`, for example `inverse_binomial`).  The
//!   converse direction generally fails: `1/(k + 1)` is the inverse of the
//!   multiplier sequence `k + 1` and maps `(1 + x)^2` to `1 + x + x^2/3`.
//!   Index-dependent weights get no adapter.
//! - `rolle_derivative_antecedent` (`derivative_antecedent`): if `g` and
//!   `Q_n` are real-rooted, so is `P_n = D^m(g Q_n)/m!`, with every root in the
//!   hull of the roots of `g` and an interval containing the roots of `Q_n`
//!   (Rolle); the converse fails (`(1 + x)^2 = (x Q)'` with
//!   `Q = 1 + x + x^2/3`).
//!
//! Inductions: the defining rule of `P_n` and a successor rule of `Q_n`
//! certified by reconstruction are matched against Rolle's theorem
//! (`source_rule_theorem`, `fitted_rule_theorem`); a completed induction is
//! transferred through an adapter only in its direction and only when the
//! adapter's hypotheses and the map conditions are certified for every index.
//! A Jensen fixed-cap test ([`jensen_fixed_cap_check`]) of an index-free
//! weight at the largest observed degree is reported as `fixed_cap`
//! evidence, never as a proof for larger degrees.

use super::*;
use crate::real_rootedness::{isolate_real_roots_exact, RealInterval};
use crate::vec_poly;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};

/// Root-domain section of a `polytool.recurrence-implicit.v1` report.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ImplicitRootDomainReport {
    /// The analysis request, retained for replay.
    pub request: RootDomainRequest,
    /// Conventions for properties, intervals, zero rows and evidence.
    pub conventions: Vec<String>,
    /// Theorem adapters for the map, with hypothesis statuses.
    pub implications: Vec<RootImplication>,
    /// Claims about `P_n` (source) and `Q_n` (transformed), each with its evidence.
    pub claims: Vec<RootDomainClaim>,
    /// Rolle match of the defining rule of `P_n`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_rule_theorem: Option<RecurrenceRootTheorem>,
    /// Rolle match of the successor rule of `Q_n` certified by reconstruction.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fitted_rule_theorem: Option<RecurrenceRootTheorem>,
    /// Jensen test of an index-free diagonal weight at the largest observed
    /// degree (`fixed_cap` evidence only).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub jensen_fixed_cap: Option<JensenFixedCapReport>,
    /// Exact checks of the supplied assumptions.
    pub assumption_checks: Vec<AssumptionCheck>,
    /// Exact root isolation of the observed rows (`P_n` and `Q_n`).
    pub observations: Vec<RowRootObservation>,
    /// Limits reached and transfers not performed.
    pub notes: Vec<String>,
}

/// An index-free weight `scale (k!)^factorial prod f_j(k)^(e_j)`.
struct IndexFreeWeight {
    factorial: i64,
    factors: Vec<(Vec<Q>, i64)>,
}

fn index_free_weight(weight: &DiagonalWeight) -> Result<Option<IndexFreeWeight>, String> {
    let linear = |t: usize| vec![Q::from_integer(BigInt::from(t)), Q::one()];
    Ok(Some(match weight {
        DiagonalWeight::FactorialPower { power } => IndexFreeWeight {
            factorial: -i64::try_from(*power).map_err(|_| "power overflow")?,
            factors: Vec::new(),
        },
        DiagonalWeight::Binomial { m } => IndexFreeWeight {
            factorial: 0,
            factors: (1..=*m).map(|t| (linear(t), 1)).collect(),
        },
        DiagonalWeight::InverseBinomial { m } => IndexFreeWeight {
            factorial: 0,
            factors: (1..=*m).map(|t| (linear(t), -1)).collect(),
        },
        DiagonalWeight::Product {
            factorial_exponent,
            factors,
            ..
        } => {
            let mut parsed = Vec::new();
            for factor in factors {
                let rows = factor
                    .coeffs
                    .coeffs
                    .iter()
                    .map(|row| {
                        row.iter()
                            .map(|c| parse_rational_coeff(c))
                            .collect::<Result<Vec<Q>, String>>()
                    })
                    .collect::<Result<Vec<_>, String>>()?;
                if rows.iter().skip(1).any(|row| !vec_poly::is_zero(row)) {
                    return Ok(None);
                }
                let f = rows
                    .first()
                    .map_or_else(Vec::new, |row| vec_poly::trim(row));
                parsed.push((f, factor.exponent));
            }
            IndexFreeWeight {
                factorial: *factorial_exponent,
                factors: parsed,
            }
        }
    }))
}

/// The adapter of the map: its implication and which directions are
/// certified for every index.
struct MapAdapter {
    implication: RootImplication,
    /// `P_n` real-rooted implies `Q_n` real-rooted, all hypotheses certified.
    forward: bool,
    /// `Q_n` real-rooted implies `P_n` real-rooted, all hypotheses certified.
    backward: bool,
    /// The factor `g` of a derivative antecedent (roots added to intervals).
    antecedent_factor: Option<Vec<Q>>,
    /// Hypotheses that are not certified, for `conditional_on`.
    uncertified: Vec<String>,
}

fn base_implication(record: &TransformationRecord, theorem: &str) -> RootImplication {
    RootImplication {
        theorem: theorem.into(),
        direction: RootImplicationDirection::SourceToTarget,
        source_property: "P_n real-rooted".into(),
        target_property: "Q_n real-rooted".into(),
        statement: String::new(),
        hypotheses: Vec::new(),
        degree_bounds: None,
        coefficient_signs: None,
        endpoint_conventions: Vec::new(),
        zero_cases: ZERO_CASES.to_string(),
        index_domain: record.index_domain.clone(),
        converse: None,
    }
}

fn exact(condition: String, holds: bool, detail: String) -> RootHypothesis {
    RootHypothesis {
        condition,
        status: if holds {
            HypothesisStatus::ExactFixedPolynomial
        } else {
            HypothesisStatus::Violated
        },
        detail,
    }
}

/// Whether `f` has only real roots, all `<= 0` (exact isolation).
fn roots_nonpositive(f: &[Q]) -> bool {
    isolate_real_roots_exact(f).real_rooted_in(&RealInterval::at_most(Q::zero()))
}

fn diagonal_adapter(
    weight: &DiagonalWeight,
    record: &TransformationRecord,
    map_conditions: &[TransformConditionStatus],
    notes: &mut Vec<String>,
) -> Result<MapAdapter, String> {
    let mut implication = base_implication(record, "polya_schur_multiplier_sequence");
    implication.hypotheses = condition_hypotheses(map_conditions);
    let map_uncertified = uncertified(map_conditions);
    implication.statement = format!(
        "Q_n = T P_n with T(x^k) = lambda_k x^k, lambda = {}: a multiplier sequence maps every \
         real-rooted polynomial to a real-rooted polynomial or zero (Pólya--Schur). Certified \
         multiplier sequences: f(k) for a polynomial f with only real roots, all <= 0 \
         (Laguerre), 1/k! (sum x^k/(k!)^2 lies in the Laguerre--Pólya class with nonnegative \
         coefficients), their termwise products and nonzero multiples",
        weight.formula()
    );
    implication.converse = Some(
        "fails in general: 1/(k + 1), the inverse of the multiplier sequence k + 1 (P = (xQ)'), \
         maps (1 + x)^2 to 1 + x + x^2/3"
            .into(),
    );
    let Some(form) = index_free_weight(weight)? else {
        implication.target_property = "none".into();
        implication.hypotheses.push(RootHypothesis {
            condition: "lambda(n,k) does not depend on n".into(),
            status: HypothesisStatus::Violated,
            detail: "a weight factor depends on n; no multiplier-sequence adapter".into(),
        });
        notes.push("the diagonal weight depends on n: no root implication is drawn".into());
        return Ok(MapAdapter {
            implication,
            forward: false,
            backward: false,
            antecedent_factor: None,
            uncertified: map_uncertified,
        });
    };
    let all_nonpositive = form.factors.iter().all(|(f, _)| roots_nonpositive(f));
    let forward_shape = form.factorial <= 0 && form.factors.iter().all(|(_, e)| *e > 0);
    let backward_shape = form.factorial >= 0 && form.factors.iter().all(|(_, e)| *e < 0);
    let forward = forward_shape && all_nonpositive;
    let backward = backward_shape && all_nonpositive;
    implication.hypotheses.push(exact(
        "lambda(n,k) does not depend on n".into(),
        true,
        "every weight factor is free of n".into(),
    ));
    let factors_text = if form.factors.is_empty() {
        "no polynomial factors".to_string()
    } else {
        form.factors
            .iter()
            .map(|(f, e)| format!("({})^({e})", poly_text(f, "k")))
            .collect::<Vec<_>>()
            .join(", ")
    };
    implication.hypotheses.push(exact(
        "every polynomial factor f_j(k) has only real roots, all <= 0".into(),
        all_nonpositive,
        format!("exact root isolation of {factors_text}"),
    ));
    implication.hypotheses.push(exact(
        "lambda (or 1/lambda) is a product of k + a factors (a >= 0) and powers of 1/k!".into(),
        forward_shape || backward_shape,
        format!(
            "(k!)^({}) and exponents {:?}",
            form.factorial,
            form.factors.iter().map(|(_, e)| *e).collect::<Vec<_>>()
        ),
    ));
    let certified = map_uncertified.is_empty();
    implication.direction = match (forward, backward) {
        (true, true) => RootImplicationDirection::Equivalent,
        (false, true) => RootImplicationDirection::TargetToSource,
        _ => RootImplicationDirection::SourceToTarget,
    };
    if !forward && !backward {
        implication.target_property = "none".into();
        notes.push(
            "neither lambda nor 1/lambda is certified as a multiplier sequence: no root \
             implication is drawn"
                .into(),
        );
    }
    Ok(MapAdapter {
        implication,
        forward: forward && certified,
        backward: backward && certified,
        antecedent_factor: None,
        uncertified: map_uncertified,
    })
}

fn antecedent_adapter(g: &[String], record: &TransformationRecord) -> Result<MapAdapter, String> {
    let (g, m) = derivative_antecedent_parameters(g).map_err(|e| e.to_string())?;
    let mut implication = base_implication(record, "rolle_derivative_antecedent");
    implication.direction = RootImplicationDirection::TargetToSource;
    implication.statement = format!(
        "P_n = D^{m}(g Q_n)/{m}! with g = {}: if g and Q_n are real-rooted with every root in an \
         interval I, then P_n is real-rooted with every root in I (Rolle's theorem; T_g is \
         injective, so P_n = 0 iff Q_n = 0)",
        poly_text(&g, "x")
    );
    let real_rooted = isolate_real_roots_exact(&g).is_real_rooted();
    implication.hypotheses.push(exact(
        "g is real-rooted".into(),
        real_rooted,
        "exact root isolation of g".into(),
    ));
    implication.endpoint_conventions =
        vec!["the roots of g are added to the interval of Q_n (convex hull)".into()];
    implication.converse = Some(
        "fails: (1 + x)^2 = (x Q)' has the antecedent Q = 1 + x + x^2/3, which is not \
         real-rooted"
            .into(),
    );
    Ok(MapAdapter {
        implication,
        forward: false,
        backward: real_rooted,
        antecedent_factor: Some(g),
        uncertified: Vec::new(),
    })
}

/// The source-side interval for `Q_n` roots in `interval` under a derivative
/// antecedent with factor `g`: the hull with the roots of `g` (exact when
/// `interval` is index-independent and already contains them).
fn antecedent_interval(
    interval: &RootInterval,
    g: &[Q],
    domain: (usize, Option<usize>),
) -> Option<RootInterval> {
    let g_isolation = isolate_real_roots_exact(g);
    if interval.is_index_independent() {
        let real = interval.at(0).ok()?;
        if g_isolation.real_rooted_in(&real) {
            return Some(interval.clone());
        }
    }
    let hull = domain_hull(interval, domain.0, domain.1).ok()?;
    Some(RootInterval::from_real(
        &match g_isolation.enclosure_hull() {
            Some(g_hull) => hull.hull(&g_hull),
            None => hull,
        },
    ))
}

fn decode_rows(rows: &[Vec<String>]) -> Result<Vec<Vec<Q>>, String> {
    rows.iter()
        .map(|row| {
            row.iter()
                .map(|c| parse_rational_coeff(c))
                .collect::<Result<Vec<Q>, String>>()
                .map(|r| vec_poly::trim(&r))
        })
        .collect()
}

/// Attach root-domain information to an implicit-transport report: exact
/// observations of `P_n` and `Q_n`, the map's theorem adapter, Rolle
/// inductions on the defining rule of `P_n` and on a successor rule of `Q_n`
/// certified by reconstruction, supplied assumptions and, for an index-free
/// diagonal weight, a fixed-cap Jensen test.  A completed induction is
/// transferred through the adapter only in its direction and only when every
/// hypothesis is certified; finite checks are never promoted to statements
/// about every index.
pub fn analyze_implicit_root_domain(
    report: &ImplicitTransformReport,
    request: &RootDomainRequest,
) -> Result<ImplicitRootDomainReport, String> {
    request.validate()?;
    let transformation = &report.request.transformation;
    let record = transformation.record();
    let origin = report.row_origin;
    let rows = decode_rows(&report.source_rows)?;
    let q_rows = decode_rows(&report.transformed_rows)?;
    let states: Vec<Vec<Vec<Q>>> = q_rows.iter().map(|q| vec![q.clone()]).collect();
    let names = vec!["Q".to_string()];
    let mut notes = vec![
        "candidate_intervals are tried for the defining rule of P_n and for a successor rule of \
         Q_n certified by reconstruction; the default window does not apply to these maps"
            .to_string(),
    ];
    let mut claims = Vec::new();
    let (isolations, observations) = observe_rows(request, origin, &rows, &states, &mut notes);
    observation_claims(&isolations, 1, &names, &mut claims, &mut notes);

    let adapter = match transformation {
        RecurrenceTransformation::Diagonal { weight } => {
            diagonal_adapter(weight, &record, &report.map_conditions, &mut notes)?
        }
        RecurrenceTransformation::DerivativeAntecedent { g } => antecedent_adapter(g, &record)?,
        other => {
            return Err(format!(
                "implicit root-domain analysis needs a diagonal or derivative_antecedent map, \
                 not {}",
                other.kind_name()
            ))
        }
    };

    // Jensen fixed-cap test of an index-free diagonal weight.
    let mut jensen = None;
    if let RecurrenceTransformation::Diagonal { weight } = transformation {
        if index_free_weight(weight)?.is_some() {
            let cap = rows
                .iter()
                .take(request.max_observed_rows)
                .filter_map(|r| vec_poly::degree(r))
                .max()
                .unwrap_or(0)
                .min(request.max_observed_degree)
                .min(ROOT_ENCLOSURES_MAX_DEGREE);
            if cap >= 1 {
                match (0..=cap)
                    .map(|k| weight.weight(origin, k))
                    .collect::<Result<Vec<_>, _>>()
                {
                    Ok(lambdas) => jensen = Some(jensen_fixed_cap_check(&lambdas, cap)?),
                    Err(error) => notes.push(format!("no Jensen test: {error}")),
                }
            }
        }
    }

    // Supplied assumptions and their transfers.
    let mut assumption_checks = Vec::new();
    for (k, assumption) in request.assumptions.iter().enumerate() {
        if assumption.component != 0 {
            return Err(format!(
                "assumption {k}: component {} does not exist",
                assumption.component
            ));
        }
        let Some((first, label)) = check_assumption(
            k,
            assumption,
            &isolations,
            origin,
            &names,
            (&mut assumption_checks, &mut claims, &mut notes),
        )?
        else {
            continue;
        };
        let to_target = assumption.subject == RootSubject::Source;
        let allowed = if to_target {
            adapter.forward
        } else {
            adapter.backward
        };
        if !allowed {
            notes.push(format!(
                "{label}: not transferred ({} gives no certified implication in this direction)",
                adapter.implication.theorem
            ));
            continue;
        }
        let (property, interval) = match (&adapter.antecedent_factor, &assumption.interval) {
            (Some(g), Some(interval)) if !to_target => (
                RootProperty::RealRootedIn,
                antecedent_interval(interval, g, (first, assumption.last_index)),
            ),
            _ => (RootProperty::RealRooted, None),
        };
        let property = if interval.is_some() {
            property
        } else {
            RootProperty::RealRooted
        };
        let mut conditional_on = vec![label];
        conditional_on.extend(adapter.uncertified.iter().cloned());
        claims.push(RootDomainClaim {
            subject: if to_target {
                RootSubject::Transformed
            } else {
                RootSubject::Source
            },
            component: 0,
            property,
            interval,
            first_index: first,
            last_index: assumption.last_index,
            evidence: RootEvidence::SuppliedAssumption,
            conditional_on,
            derivation: format!(
                "transferred from assumption {k} by {}",
                adapter.implication.theorem
            ),
        });
    }

    let candidates: Vec<(RootInterval, String)> = request
        .candidate_intervals
        .iter()
        .map(|i| (i.clone(), "supplied".to_string()))
        .collect();

    // Rolle induction on the defining rule of P_n.
    let (source_rule, source_first, source_seeds) = report.source_rule.to_recurrence_parts()?;
    let source_valid_from = source_first + source_seeds.len();
    let source_outcome = rolle_core(
        &RolleSetup {
            rule: &source_rule,
            name: "P".into(),
            first: source_first,
            valid_from: source_valid_from,
            seeds: source_seeds
                .iter()
                .enumerate()
                .map(|(i, seed)| (source_first + i, vec_poly::trim(seed)))
                .collect(),
            candidates: candidates.clone(),
            limit: request.max_candidate_intervals,
            schedule: None,
        },
        &mut notes,
    );
    if let Some(interval) = &source_outcome.proven {
        claims.push(RootDomainClaim {
            subject: RootSubject::Source,
            component: 0,
            property: RootProperty::RealRootedIn,
            interval: Some(interval.clone()),
            first_index: source_first,
            last_index: None,
            evidence: RootEvidence::CompletedRecurrenceInduction,
            conditional_on: Vec::new(),
            derivation: format!(
                "rolle_derivative_of_product on the defining rule of P_n with g and the seeds \
                 checked exactly on {interval}"
            ),
        });
        if adapter.forward {
            claims.push(RootDomainClaim {
                subject: RootSubject::Transformed,
                component: 0,
                property: RootProperty::RealRooted,
                interval: None,
                first_index: source_first,
                last_index: None,
                evidence: RootEvidence::CertifiedFamilyInvariant,
                conditional_on: Vec::new(),
                derivation: format!(
                    "the completed induction for P_n and {} (source to target; the weight is a \
                     certified multiplier sequence and the map conditions hold for all indices)",
                    adapter.implication.theorem
                ),
            });
        } else {
            notes.push(format!(
                "the induction for P_n is not transferred to Q_n: {} gives no certified \
                 implication from P_n to Q_n",
                adapter.implication.theorem
            ));
        }
    }

    // Rolle induction on a successor rule of Q_n certified by reconstruction.
    let fitted = report.independent_search.as_ref().and_then(|search| {
        (search.evidence == CERTIFIED_BY_RECONSTRUCTION_EVIDENCE)
            .then_some(search.scalar.as_ref())
            .flatten()
    });
    let fitted_rule_theorem = match fitted {
        None => None,
        Some(json) => {
            let (rule, first, seeds) = json.to_recurrence_parts()?;
            let valid_from = first + seeds.len();
            let outcome = rolle_core(
                &RolleSetup {
                    rule: &rule,
                    name: "Q".into(),
                    first,
                    valid_from,
                    seeds: seeds
                        .iter()
                        .enumerate()
                        .map(|(i, seed)| (first + i, vec_poly::trim(seed)))
                        .collect(),
                    candidates,
                    limit: request.max_candidate_intervals,
                    schedule: None,
                },
                &mut notes,
            );
            if let Some(interval) = &outcome.proven {
                claims.push(RootDomainClaim {
                    subject: RootSubject::Transformed,
                    component: 0,
                    property: RootProperty::RealRootedIn,
                    interval: Some(interval.clone()),
                    first_index: first,
                    last_index: None,
                    evidence: RootEvidence::CompletedRecurrenceInduction,
                    conditional_on: Vec::new(),
                    derivation: format!(
                        "rolle_derivative_of_product on the successor rule of Q_n certified by \
                         reconstruction, with g and the seeds checked exactly on {interval}"
                    ),
                });
                if adapter.backward {
                    let interval = adapter
                        .antecedent_factor
                        .as_ref()
                        .and_then(|g| antecedent_interval(interval, g, (first, None)));
                    claims.push(RootDomainClaim {
                        subject: RootSubject::Source,
                        component: 0,
                        property: if interval.is_some() {
                            RootProperty::RealRootedIn
                        } else {
                            RootProperty::RealRooted
                        },
                        interval,
                        first_index: first,
                        last_index: None,
                        evidence: RootEvidence::CertifiedFamilyInvariant,
                        conditional_on: Vec::new(),
                        derivation: format!(
                            "the completed induction for Q_n and {} (target to source; every \
                             hypothesis certified for all indices)",
                            adapter.implication.theorem
                        ),
                    });
                } else {
                    notes.push(format!(
                        "the induction for Q_n is not transferred to P_n: {} gives no \
                         certified implication from Q_n to P_n",
                        adapter.implication.theorem
                    ));
                }
            }
            Some(outcome.theorem)
        }
    };

    let mut conventions = conventions();
    conventions.push(
        "fixed_cap: a Jensen test at degree cap N proves preservation of real-rootedness for \
         degrees at most N only"
            .into(),
    );
    Ok(ImplicitRootDomainReport {
        request: request.clone(),
        conventions,
        implications: vec![adapter.implication],
        claims,
        source_rule_theorem: Some(source_outcome.theorem),
        fitted_rule_theorem,
        jensen_fixed_cap: jensen,
        assumption_checks,
        observations,
        notes,
    })
}

#[cfg(test)]
mod tests;
