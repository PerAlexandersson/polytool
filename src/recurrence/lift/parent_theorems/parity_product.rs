//! Shared-factor parity products of constant-coefficient second-order families.
//!
//! For `P_(n+3)=a P_(n+2)+b P_(n+1)+c P_n`, match exactly
//! `c=mu^3`, `b=-a mu`, `beta=a-mu` in `Q[x]`. When `a!=0`, `mu=-b/a`
//! is unique. Otherwise `b=0` and the unique polynomial cube root determines
//! `mu`. This is precisely `(z-mu)(z^2-beta z+mu^2)`.
//!
//! Enumerate all monic divisors `t0` of `gcd(P0,P1)` (at most 4096), set
//! `s0=P0/t0`, `s1=P1/t0`, and require polynomial `t1=P2/s1` and the
//! Casoratian identity. Nonzero scalar reallocations only rescale the two
//! families inversely. Earlier prescribed rows are checked exactly.
//! The factors use the existing nonrecursive residue-class menu at depth one;
//! this route is excluded from that menu. No simplicity or strict interlacing
//! is inferred, since consecutive parity products share a factor.

use super::degree_law_reversal::constant_coefficient_multipliers;
use super::*;
use crate::factorization::irreducible_factorization_rational_coeffs_with_budget;
use crate::Polynomial;

type Poly = Polynomial<Q>;
const REFERENCE: &str = "RealRooted.ThreeTermRecurrence.ParityProduct: \
    eq_parityProduct_of_thirdOrderRecurrence, splits_of_parityProduct, \
    interlaces_of_parityProduct; exact polynomial identities and factor proofs";

/// Maximum number of monic gcd divisors in each seed window.
pub const MAXIMUM_PARITY_PRODUCT_DIVISORS: usize = 4096;
/// Maximum nesting depth: factors use only existing nonrecursive routes.
pub const MAXIMUM_PARITY_PRODUCT_PROOF_DEPTH: usize = 1;

/// Exact parity-product data. Polynomial coefficients are ascending rational strings.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ParityProductCertificate {
    /// `mu(x)` in the cubic factorization.
    pub mu: Vec<String>,
    /// `beta(x)` in the common factor recurrence.
    pub beta: Vec<String>,
    /// Index relative to the supplied first row of the parity-product start.
    pub exceptional_rows: usize,
    /// Selected monic divisor `t0` of the initial two rows' gcd.
    pub initial_common_divisor: Vec<String>,
    /// Divisors tested up to and including the selected allocation.
    pub divisors_examined: usize,
    /// Total number of monic gcd divisors in the selected window, at most 4096.
    pub divisor_lattice_size: usize,
    /// Two proved theorem reports, in order `s`, `t`, both starting at zero.
    pub factor_proofs: Vec<RealRootednessTheoremReport>,
    /// Induction-start retry bound passed to the factors, at most 16.
    pub maximum_start_delay: usize,
    /// Fixed maximum proof depth, one.
    pub proof_depth: usize,
}

impl ParityProductCertificate {
    /// Replay this certificate in its owning theorem report.
    ///
    /// # Errors
    /// Returns a verification error for changed identities, proofs or bounds.
    pub fn verify(&self, report: &RealRootednessTheoremReport) -> Result<(), RecurrenceError> {
        if report.theorem != RealRootednessTheorem::ParityProduct
            || report.parity_product_certificate() != Some(self)
        {
            return Err(RecurrenceError::VerificationFailed(
                "parity-product certificate does not belong to this report".into(),
            ));
        }
        report.verify()
    }
}

/// Prove a stationary cubic family by shared-factor parity products.
///
/// Uses every monic gcd divisor within the 4096-element lattice bound and
/// existing nonrecursive factor routes at depth one. Requested root regions
/// are passed to both factors. Earlier rows are checked separately. The last
/// three seeds and up to `maximum_start_delay` later windows are examined.
/// A successful proof requires every supplied and generated prefix row nonzero.
/// Reversal is available through the shared degree-law reversal closure.
///
/// # Errors
/// Returns [`RecurrenceError`] for invalid records or failed exact arithmetic.
pub fn check_parity_product_theorem(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    check_real_rootedness_theorem(
        RealRootednessTheorem::ParityProduct,
        rule,
        first_index,
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
    check_at_supplied_seeds: parity_product_report,
    start_policy: StartPolicy::OwnStarts,
    is_preferred: any_proof,
    certificate_name: Some("parity-product"),
    discovery: DiscoveryPolicy::WhenStronger,
    in_residue_class_menu: false,
    in_index_section_menu: false,
    in_degree_law_reversal_menu: true,
    parent: ParentPolicy {
        stage: ParentStage::WhileUnproved,
        region: ParentRegion::OpenNonpositiveRayThenUnrestricted,
        record_only_if_applicable: true,
    },
};

fn strings(p: &Poly) -> Vec<String> {
    vec_poly::trim(p.coeffs())
        .iter()
        .map(ToString::to_string)
        .collect()
}
fn quotient(a: &Poly, b: &Poly) -> Option<Poly> {
    if b.is_zero() {
        return None;
    }
    let (q, r) = a.div_rem(b);
    r.is_zero().then_some(q)
}
fn recognize(multipliers: &[Poly]) -> Option<(Poly, Poly)> {
    if multipliers.len() != 3 || multipliers[2].is_zero() {
        return None;
    }
    let [a, b, c] = <&[Poly; 3]>::try_from(multipliers).ok()?;
    let mu = if !a.is_zero() {
        quotient(&(-b.clone()), a)?
    } else if b.is_zero() {
        index_section::polynomial_power_root(c, 3)?
    } else {
        return None;
    };
    if mu.pow(3) != *c || -(a.clone() * mu.clone()) != *b {
        return None;
    }
    Some((mu.clone(), a.clone() - mu))
}
fn divisors(gcd: &Poly) -> Option<Vec<Poly>> {
    let f = irreducible_factorization_rational_coeffs_with_budget(
        gcd.coeffs(),
        MAXIMUM_PARITY_PRODUCT_DIVISORS,
    )
    .ok()?;
    let count = f
        .factors
        .iter()
        .try_fold(1usize, |n, (_, e)| n.checked_mul(e + 1))?;
    if count > MAXIMUM_PARITY_PRODUCT_DIVISORS {
        return None;
    }
    let mut out = vec![Poly::one()];
    for (factor, multiplicity) in f.factors {
        let factor = Poly::new(factor.into_iter().map(Q::from_integer).collect()).make_monic();
        let prefix = out.clone();
        let mut power = Poly::one();
        for _ in 0..multiplicity {
            power = power * factor.clone();
            out.extend(prefix.iter().map(|p| p.clone() * power.clone()));
        }
    }
    Some(out)
}
fn factor_rule(beta: &Poly, mu: &Poly) -> Recurrence {
    Recurrence {
        terms: [(1, beta.clone()), (2, -mu.pow(2))]
            .into_iter()
            .filter(|(_, p)| !p.is_zero())
            .map(|(lag, p)| {
                RecurrenceTerm::new(
                    lag,
                    0,
                    BivarPoly::from_coefficients_normalized(vec![p.coeffs().to_vec()]),
                )
            })
            .collect(),
        denominator: None,
        inhomogeneous: None,
    }
}
fn prove_factor(
    rule: &Recurrence,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<Option<RealRootednessTheoremReport>, RecurrenceError> {
    let mut fallback = None;
    for route in ROUTES.iter().filter(|r| r.policy.in_residue_class_menu) {
        let attempt = check_with_first_proving_start(route.theorem, rule, 0, seeds, options)?;
        if !attempt.is_proved()
            || !attempt
                .conclusion
                .as_ref()
                .is_some_and(|c| c.first_index == 0 && !c.has_zero_residue_classes())
        {
            continue;
        }
        if attempt
            .conclusion
            .as_ref()
            .is_some_and(|c| c.interlacing.is_some())
        {
            return Ok(Some(attempt));
        }
        if fallback.is_none() {
            fallback = Some(attempt);
        }
    }
    Ok(fallback)
}

fn parity_product_report(
    rule: &Recurrence,
    first: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    if options.maximum_start_delay > index_section::MAXIMUM_PARENT_START_DELAY {
        return Err(RecurrenceError::InvalidInput(
            "parity-product start delay exceeds 16".into(),
        ));
    }
    let mut report = RealRootednessTheoremReport::new(
        RealRootednessTheorem::ParityProduct,
        rule,
        first,
        seeds,
        options,
    );
    if seeds.len() < 3 {
        return Ok(report.not_applicable("three seed rows are required"));
    }
    let Ok(multipliers) = constant_coefficient_multipliers(rule) else {
        return Ok(report.not_applicable("the rule is not stationary without derivatives"));
    };
    let Some((mu, beta)) = recognize(&multipliers) else {
        return Ok(report.not_applicable("the cubic is not (z-mu)(z^2-beta*z+mu^2) over Q[x]"));
    };
    let rule_factors = factor_rule(&beta, &mu);
    let mut rows = seeds.to_vec();
    let mut reason = "no polynomial seed allocation has two proved factor families";
    for delay in 0..=options.maximum_start_delay {
        if delay > 0 {
            rows = rule
                .generate_rows_rational(seeds, first, seeds.len() + delay)
                .map_err(|e| RecurrenceError::Unsupported(e.to_string()))?;
        }
        let prefix = rows.len() - 3;
        let p: Vec<Poly> = rows[prefix..].iter().cloned().map(Poly::new).collect();
        if p.iter().any(Poly::is_zero) {
            continue;
        }
        let residual = p[1].clone()
            * (p[2].clone() - beta.clone() * p[1].clone() + mu.pow(2) * p[0].clone())
            - mu.clone() * (p[1].pow(2) - p[0].clone() * p[2].clone());
        if !residual.is_zero() {
            // The residual at the next window is mu^2 times this one.
            // Since order three implies mu != 0, later windows cannot repair it.
            return Ok(report.not_applicable("the seed Casoratian compatibility identity fails"));
        }
        let Some(divisors) = divisors(&p[0].gcd(&p[1])) else {
            reason = "gcd factorization or divisor lattice exceeds the 4096 bound";
            continue;
        };
        let lattice_size = divisors.len();
        for (position, t0) in divisors.into_iter().enumerate() {
            let (Some(s0), Some(s1)) = (quotient(&p[0], &t0), quotient(&p[1], &t0)) else {
                continue;
            };
            let Some(t1) = quotient(&p[2], &s1) else {
                continue;
            };
            let s2 = beta.clone() * s1.clone() - mu.pow(2) * s0.clone();
            if s1.clone() * t1.clone() - s2 * t0.clone()
                != mu.clone() * (s1.clone() * t0.clone() - s0.clone() * t1.clone())
            {
                continue;
            }
            let sr = vec![vec_poly::trim(s0.coeffs()), vec_poly::trim(s1.coeffs())];
            let tr = vec![vec_poly::trim(t0.coeffs()), vec_poly::trim(t1.coeffs())];
            let Some(s) = prove_factor(&rule_factors, &sr, options)? else {
                continue;
            };
            let Some(t) = prove_factor(&rule_factors, &tr, options)? else {
                continue;
            };
            let start = first.checked_add(prefix).ok_or_else(|| {
                RecurrenceError::InvalidInput("parity-product index overflow".into())
            })?;
            let sc = s
                .conclusion
                .as_ref()
                .expect("proved factor has a conclusion");
            let tc = t
                .conclusion
                .as_ref()
                .expect("proved factor has a conclusion");
            let interlacing_from = match (sc.interlacing_from, tc.interlacing_from) {
                (Some(a), Some(b)) if sc.interlacing.is_some() && tc.interlacing.is_some() => {
                    // At product index 2m+1 the factor pair is s_m,s_(m+1);
                    // at 2m+2 it is t_m,t_(m+1). Both apply from
                    // 2*max(a,b)-1, including the first pair when a=b=1.
                    a.max(b)
                        .saturating_sub(1)
                        .checked_mul(2)
                        .and_then(|i| i.checked_add(1))
                        .and_then(|i| start.checked_add(i))
                }
                _ => None,
            };
            // Only the requested region is claimed: earlier prefix rows are
            // checked in it by earlier_seeds_hypothesis below.
            let root_region = options.root_region.as_ref().map(open_root_interval);
            report.normal_form = format!(
                "P_({start}+2m)=s_m*t_m, P_({start}+2m+1)=s_(m+1)*t_m; \
                both factors satisfy f_(m+2)=({beta})*f_(m+1)-({mu})^2*f_m"
            );
            report.hypotheses =
                earlier_seeds_hypothesis(&rows[..prefix], first, options.root_region.as_ref())
                    .into_iter()
                    .collect();
            report.hypotheses.push(exact_hypothesis("parity_product_identity",
                "cubic factorization, three seed products and Casoratian compatibility hold over Q[x]",
                Ok(format!("monic gcd divisor {} of {lattice_size}; exact identities verified",position+1))));
            report.hypotheses.push(exact_hypothesis(
                "factor_proofs",
                "both second-order families have replayable all-index proofs from zero",
                Ok(format!(
                    "s: {}; t: {}; depth 1",
                    s.theorem.as_str(),
                    t.theorem.as_str()
                )),
            ));
            report.lean_route = Some("parityProduct".into());
            report.lean_declarations = [
                "RealRooted.eq_parityProduct_of_thirdOrderRecurrence",
                "RealRooted.splits_of_parityProduct",
            ]
            .map(String::from)
            .to_vec();
            if interlacing_from.is_some() {
                report
                    .lean_declarations
                    .push("RealRooted.interlaces_of_parityProduct".into());
            }
            report.certificate = Some(RouteCertificate::ParityProduct(Box::new(
                ParityProductCertificate {
                    mu: strings(&mu),
                    beta: strings(&beta),
                    exceptional_rows: prefix,
                    initial_common_divisor: strings(&t0),
                    divisors_examined: position + 1,
                    divisor_lattice_size: lattice_size,
                    factor_proofs: vec![s, t],
                    maximum_start_delay: options.maximum_start_delay,
                    proof_depth: 1,
                },
            )));
            return Ok(report.finish(RealRootednessTheoremConclusion {
                first_index: first,
                zero_residue_classes: None,
                interlacing_from,
                interlacing: interlacing_from.map(|_| RowInterlacing::Weak),
                simple_roots: false,
                root_region,
                statement: format!(
                    "every P_n, n >= {first}, is nonzero and real-rooted by parity products{}; \
                    no simplicity or strict interlacing is asserted",
                    interlacing_from.map_or(String::new(), |i| format!(
                        "; consecutive rows weakly interlace from n = {i}"
                    ))
                ),
            }));
        }
    }
    Ok(report.not_applicable(reason))
}

pub(super) fn verify_certificate_in_report(
    c: &ParityProductCertificate,
    report: &RealRootednessTheoremReport,
) -> Result<(), RecurrenceError> {
    if c.factor_proofs.len() != 2
        || c.proof_depth != 1
        || c.maximum_start_delay > index_section::MAXIMUM_PARENT_START_DELAY
        || c.divisor_lattice_size > MAXIMUM_PARITY_PRODUCT_DIVISORS
        || c.divisors_examined == 0
        || c.divisors_examined > c.divisor_lattice_size
    {
        return Err(RecurrenceError::VerificationFailed(
            "invalid parity-product bounds".into(),
        ));
    }
    for proof in &c.factor_proofs {
        if !route(proof.theorem).policy.in_residue_class_menu || !proof.is_proved() {
            return Err(RecurrenceError::VerificationFailed(
                "invalid parity-product factor route".into(),
            ));
        }
        proof.verify()?;
    }
    // The owning report's verifier re-derives the identities and exact records.
    if report.parity_product_certificate() != Some(c) {
        return Err(RecurrenceError::VerificationFailed(
            "parity-product certificate mismatch".into(),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
