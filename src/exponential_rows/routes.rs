//! Symbolic seed/substitution routes and bounded-rank reference certificates.
use super::*;
use crate::{check_weak_interlacing_bigint_coeffs, is_real_rooted_bigint_coeffs};
use exact_algebra::denominators::clear_row_denominators;
use num_traits::Signed;

const SUBSTITUTION:&str="https://github.com/PerAlexandersson/RealRooted/blob/main/RealRooted/ExponentialRows/Substitution.lean";
const QUARTIC:&str="https://github.com/PerAlexandersson/RealRooted/blob/main/RealRooted/ExponentialRows/QuarticGrace.lean";
const PAPER: &str = "https://github.com/PerAlexandersson/RealRooted/issues/2415";

/// Evidence level; an exact finite certificate never becomes an all-rank theorem.
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum ExponentialEvidenceLevel {
    /// Symbolic identity and checked hypotheses of cited theorems; no Lean replay.
    TheoremBacked,
    /// Published or manuscript argument, without a formal proof claim.
    PaperProof,
    /// Candidate paper claim explicitly still under review.
    UnderReview,
    /// Exact certificate for this rank only.
    FiniteRank,
}
/// One exact symbolic route, including its scope and public references.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct ExponentialRouteReport {
    /// Stable route identifier.
    pub route: String,
    /// Evidence level, independent of finite row checks.
    pub evidence: ExponentialEvidenceLevel,
    /// proved, paper_proof, necessary_condition_passed, obstruction, or under_review.
    pub status: String,
    /// Precisely scoped conclusion and checked hypotheses.
    pub statement: String,
    /// Public source URLs and declaration names.
    pub references: Vec<String>,
    /// Seed/substitution chain; empty for primitive routes.
    pub premises: Vec<ExponentialRouteReport>,
    /// Whether the route also transfers consecutive weak interlacing, at its evidence level.
    pub consecutive_weak_interlacing: bool,
}

fn finite(spec: &BlockWeightSpec) -> Result<Option<Vec<Q>>, ExponentialRowError> {
    Ok(match spec {
        BlockWeightSpec::Finite { weights } => Some(
            weights
                .iter()
                .map(|w| rational(w))
                .collect::<Result<_, _>>()?,
        ),
        BlockWeightSpec::Sizes { sizes } => {
            let max = sizes.iter().copied().max().unwrap_or(0);
            // Seed recognition only needs degree at most four. Keep huge
            // sparse sizes sparse; the obstruction route uses support below.
            if max > 4 {
                return Ok(None);
            }
            let mut w = vec![Q::zero(); max];
            for &r in sizes {
                w[r - 1] = Q::one();
            }
            Some(w)
        }
        _ => None,
    }
    .map(|mut w| {
        while w.last().is_some_and(Zero::is_zero) {
            w.pop();
        }
        w
    }))
}
fn singleton(spec: &BlockWeightSpec) -> Result<Option<Q>, ExponentialRowError> {
    Ok(finite(spec)?.and_then(|w| match w.len() {
        0 => Some(Q::zero()),
        1 => Some(w[0].clone()),
        _ => None,
    }))
}
fn support(spec: &BlockWeightSpec) -> Result<Option<Vec<(usize, Q)>>, ExponentialRowError> {
    match spec {
        BlockWeightSpec::Finite { weights } => Ok(Some(
            weights
                .iter()
                .enumerate()
                .map(|(i, w)| rational(w).map(|q| (i + 1, q)))
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
                .filter(|(_, w)| !w.is_zero())
                .collect(),
        )),
        BlockWeightSpec::Sizes { sizes } => {
            let mut sizes = sizes.clone();
            sizes.sort_unstable();
            sizes.dedup();
            Ok(Some(sizes.into_iter().map(|r| (r, Q::one())).collect()))
        }
        _ => Ok(None),
    }
}
fn nonnegative_finite(spec: &BlockWeightSpec) -> Result<bool, ExponentialRowError> {
    Ok(support(spec)?.is_some_and(|weights| weights.iter().all(|(_, w)| !w.is_negative())))
}

fn proved(
    route: &str,
    statement: &str,
    reference: &str,
    decl: &str,
    interlacing: bool,
) -> ExponentialRouteReport {
    ExponentialRouteReport {
        route: route.into(),
        evidence: ExponentialEvidenceLevel::TheoremBacked,
        status: "proved".into(),
        statement: statement.into(),
        references: vec![reference.into(), decl.into()],
        premises: vec![],
        consecutive_weak_interlacing: interlacing,
    }
}
fn seed_route(
    family: &ExponentialFamily,
) -> Result<Option<ExponentialRouteReport>, ExponentialRowError> {
    if let (Some(a), Some(b)) = (singleton(&family.unmarked)?, singleton(&family.marked)?) {
        if !a.is_negative() && !b.is_negative() {
            return Ok(Some(proved("binomial_seed","F_n=(alpha_1+beta_1*t)^n; nonnegative singleton weights. All rows split (zero rows allowed). Algebraic binomial identity and PF powers.",SUBSTITUTION,"RealRooted.ExponentialRows.insertionRow_affine_X_X",true)));
        }
    }
    if singleton(&family.unmarked)? == Some(Q::one()) {
        if let Some(w) = finite(&family.marked)? {
            if w.len() <= 4
                && w.first().is_none_or(Zero::is_zero)
                && w.iter().all(|x| !x.is_negative())
            {
                return Ok(Some(proved(if w.len()<=2 {"hermite_seed"}else{"nonnegative_quartic_seed"},"A=z, B=a*z^2+b*z^3+c*z^4, with a,b,c>=0 exactly. All ranks are nonzero and real-rooted. No consecutive interlacing claim is inferred from the quartic theorem.",QUARTIC,"RealRooted.ExponentialRows.insertionRow_isRealRooted_quartic_of_nonneg",false)));
            }
        }
        if matches!(family.marked, BlockWeightSpec::Tail { minimum: 2 }) {
            let mut report=proved("threshold_two_seed","A=z; B=exp(z)-1-z. All ranks are nonzero and real-rooted, with consecutive directed weak interlacing.",SUBSTITUTION,"RealRooted.ExponentialRows.insertionRow_isRealRooted_threshold_two",true);
            report
                .references
                .push("RealRooted.ExponentialRows.insertionRow_strictInterl_threshold_two".into());
            return Ok(Some(report));
        }
    }
    Ok(None)
}
// A finite nonzero ODE kernel with nonnegative parameters has lambda=0,
// hence Phi(x)=a*x+c*x^2. Normalize a=1 when a>0, or the first seed
// coefficient to one when a=0. Every candidate is checked as a full polynomial
// identity, including all coefficients through the known maximum degree eight.
fn finite_quadratic_substitution(
    family: &ExponentialFamily,
) -> Result<Option<(BlockWeightSpec, Q, Q, Q)>, ExponentialRowError> {
    let Some(weights) = support(&family.marked)? else {
        return Ok(None);
    };
    if weights.iter().any(|(r, w)| *r > 8 || w.is_negative()) {
        return Ok(None);
    }
    let mut p = vec![Q::zero(); 9];
    for (r, w) in weights {
        p[r] = w / Q::from_integer(factorial_bigint(r));
    }
    let two = Q::from_integer(2.into());
    let candidate = |seed: Vec<Q>, a: Q, c: Q| {
        if a.is_negative() || c.is_negative() || seed.iter().any(Signed::is_negative) {
            return None;
        }
        let reconstructed = vec_poly::add(
            &vec_poly::scale(&seed, &a),
            &vec_poly::scale(&vec_poly::mul(&seed, &seed), &c),
        );
        if reconstructed != vec_poly::trim(&p) {
            return None;
        }
        let weights = seed
            .iter()
            .enumerate()
            .skip(1)
            .map(|(r, w)| (w * Q::from_integer(factorial_bigint(r))).to_string())
            .collect();
        Some((BlockWeightSpec::Finite { weights }, a, &two * c, Q::zero()))
    };
    if vec_poly::degree(&p).unwrap_or(0) <= 2 {
        // This shortcut needs an admissible singleton seed. Without it, a
        // linear target would recursively recognize itself as its own seed.
        if !singleton(&family.unmarked)?.is_some_and(|a| !a.is_negative()) {
            return Ok(None);
        }
        return Ok(candidate(
            vec![Q::zero(), Q::one()],
            p[1].clone(),
            p[2].clone(),
        ));
    }
    if singleton(&family.unmarked)? != Some(Q::one()) || !p[1].is_zero() {
        return Ok(None);
    }
    let mut candidates = Vec::new();
    // Phi(x)=x+c*x^2, with a quartic seed having no linear term.
    if p[2].is_positive() && p[3].is_positive() {
        let c = &p[5] / (&two * &p[2] * &p[3]);
        candidates.push((
            vec![
                Q::zero(),
                Q::zero(),
                p[2].clone(),
                p[3].clone(),
                &p[4] - &c * &p[2] * &p[2],
            ],
            Q::one(),
            c,
        ));
    } else if p[2].is_positive() {
        if p[6].is_positive() && p[8].is_positive() {
            let fourth = &two * &p[2] * &p[8] / &p[6];
            let c = &p[8] / (&fourth * &fourth);
            candidates.push((
                vec![Q::zero(), Q::zero(), p[2].clone(), Q::zero(), fourth],
                Q::one(),
                c,
            ));
        } else if p[8].is_zero() {
            let c = &p[4] / (&p[2] * &p[2]);
            candidates.push((vec![Q::zero(), Q::zero(), p[2].clone()], Q::one(), c));
        }
    } else if p[3].is_positive() {
        let c = &p[6] / (&p[3] * &p[3]);
        candidates.push((
            vec![Q::zero(), Q::zero(), Q::zero(), p[3].clone(), p[4].clone()],
            Q::one(),
            c,
        ));
    } else if p[4].is_positive() {
        let c = &p[8] / (&p[4] * &p[4]);
        candidates.push((
            vec![Q::zero(), Q::zero(), Q::zero(), Q::zero(), p[4].clone()],
            Q::one(),
            c,
        ));
    }
    // Phi(x)=c*x^2. Normalize the first nonzero seed coefficient to one.
    if p[2].is_zero() && p[3].is_zero() {
        if p[4].is_positive() {
            let third = &p[5] / (&two * &p[4]);
            let fourth = (&p[6] / &p[4] - &third * &third) / &two;
            candidates.push((
                vec![Q::zero(), Q::zero(), Q::one(), third, fourth],
                Q::zero(),
                p[4].clone(),
            ));
        } else if p[5].is_zero() && p[6].is_positive() {
            let fourth = &p[7] / (&two * &p[6]);
            candidates.push((
                vec![Q::zero(), Q::zero(), Q::zero(), Q::one(), fourth],
                Q::zero(),
                p[6].clone(),
            ));
        } else if p[5].is_zero() && p[6].is_zero() && p[7].is_zero() && p[8].is_positive() {
            candidates.push((
                vec![Q::zero(), Q::zero(), Q::zero(), Q::zero(), Q::one()],
                Q::zero(),
                p[8].clone(),
            ));
        }
    }
    Ok(candidates
        .into_iter()
        .find_map(|(seed, a, c)| candidate(seed, a, c)))
}

fn substitution_route(
    family: &ExponentialFamily,
) -> Result<Option<ExponentialRouteReport>, ExponentialRowError> {
    if let Some(seed) = seed_route(family)? {
        return Ok(Some(seed));
    }
    if let BlockWeightSpec::RationalOdeSubstitution {
        seed,
        a,
        b,
        lambda,
        d,
    } = &family.marked
    {
        let (a, b, l, d) = (rational(a)?, rational(b)?, rational(lambda)?, rational(d)?);
        if d.is_zero() {
            let normalized = ExponentialFamily {
                unmarked: family.unmarked.clone(),
                marked: BlockWeightSpec::OdeSubstitution {
                    seed: seed.clone(),
                    a: a.to_string(),
                    b: b.to_string(),
                    lambda: l.to_string(),
                },
            };
            return substitution_route(&normalized);
        }
        if singleton(&family.unmarked)? != Some(Q::one())
            || [&a, &b, &l, &d].iter().any(|p| p.is_negative())
            || !nonnegative_series(seed)?
        {
            return Ok(None);
        }
        let base = ExponentialFamily {
            unmarked: family.unmarked.clone(),
            marked: seed.as_ref().clone(),
        };
        let Some(premise) = substitution_route(&base)? else {
            return Ok(None);
        };
        let interlacing = premise.consecutive_weak_interlacing;
        return Ok(Some(ExponentialRouteReport {
            route:"rational_ode_substitution_paper".into(),evidence:ExponentialEvidenceLevel::PaperProof,status:"paper_proof".into(),
            statement:format!("Unreviewed paper route: A=z, Phi(0)=0, (1-d*x)Phi'=a+b*x+lambda*Phi with a={a}, b={b}, lambda={l}, d={d} nonnegative; exact formal composition with the recognized seed. All-rank splitness and any recorded interlacing transfer are conditional on the cited paper argument. No Lean theorem for d>0 is claimed."),
            references:vec!["https://github.com/PerAlexandersson/RealRooted/issues/2415#issuecomment-7141592218".into()],premises:vec![premise],consecutive_weak_interlacing:interlacing,
        }));
    }
    if let BlockWeightSpec::PolynomialSubstitution { seed, coefficients } = &family.marked {
        // The cited round-four paper argument is explicitly scoped to A=z.
        let p = coefficients
            .iter()
            .map(|c| rational(c))
            .collect::<Result<Vec<_>, _>>()?;
        if singleton(&family.unmarked)? != Some(Q::one())
            || vec_poly::degree(&p).is_none()
            || !pf(&p)
            || !nonnegative_series(seed)?
        {
            return Ok(None);
        }
        let base = ExponentialFamily {
            unmarked: family.unmarked.clone(),
            marked: seed.as_ref().clone(),
        };
        let Some(premise) = substitution_route(&base)? else {
            return Ok(None);
        };
        let interlacing = premise.consecutive_weak_interlacing;
        return Ok(Some(ExponentialRouteReport {
            route:"polynomial_substitution_paper".into(),evidence:ExponentialEvidenceLevel::PaperProof,status:"paper_proof".into(),
            statement:"Paper proof route, not yet reviewed or formalized: A=z, exact P(B_seed)-P(0); P is nonzero with positive leading coefficient and only real nonpositive zeros (checked exactly using nonnegative coefficients and root counts). All-rank splitness follows conditional on the cited paper argument. No Lean closure theorem is claimed.".into(),
            references:vec!["https://github.com/PerAlexandersson/RealRooted/issues/2415#issuecomment-7141592218".into()],premises:vec![premise],consecutive_weak_interlacing:interlacing,
        }));
    }
    let (inner, a, b, l) = match &family.marked {
        BlockWeightSpec::OdeSubstitution { seed, a, b, lambda } => (
            seed.as_ref().clone(),
            rational(a)?,
            rational(b)?,
            rational(lambda)?,
        ),
        BlockWeightSpec::Tail { minimum } if *minimum <= 2 => (
            BlockWeightSpec::Sizes { sizes: vec![1] },
            if *minimum == 1 { Q::one() } else { Q::zero() },
            if *minimum == 2 { Q::one() } else { Q::zero() },
            Q::one(),
        ),
        BlockWeightSpec::Exponential {
            scale,
            rate,
            minimum,
        } if *minimum <= 2 => {
            let s = rational(scale)?;
            let r = rational(rate)?;
            if s.is_negative() || r.is_negative() {
                return Ok(None);
            }
            (
                BlockWeightSpec::Sizes { sizes: vec![1] },
                if *minimum == 1 { &s * &r } else { Q::zero() },
                if *minimum == 2 {
                    &s * &r * &r
                } else {
                    Q::zero()
                },
                r,
            )
        }
        _ => match finite_quadratic_substitution(family)? {
            Some(parameters) => parameters,
            None => return Ok(None),
        },
    };
    if a.is_negative()
        || b.is_negative()
        || l.is_negative()
        || !nonnegative_series(&family.unmarked)?
        || !nonnegative_series(&inner)?
    {
        return Ok(None);
    }
    let seed_spec = serde_json::to_string(&inner).expect("serializable block specification");
    let base = ExponentialFamily {
        unmarked: family.unmarked.clone(),
        marked: inner,
    };
    let Some(premise) = substitution_route(&base)? else {
        return Ok(None);
    };
    let interlacing = premise.consecutive_weak_interlacing;
    let mut route=proved("ode_substitution",&format!("Exact formal composition Phi(B_seed), Phi(0)=0, Phi'=a+b*x+lambda*Phi; a={a}, b={b}, lambda={l} are nonnegative. B_seed={seed_spec}. Nonnegative unmarked/seed coefficients checked structurally. All rows split by the seed theorem and substitution closure. Rust identity/hypothesis recognition is not a Lean replay."),SUBSTITUTION,"RealRooted.ExponentialRows.insertionRow_splits_subst",interlacing);
    if interlacing {
        route
            .references
            .push("RealRooted.ExponentialRows.insertionRow_interl_subst".into());
    }
    if premise.evidence != ExponentialEvidenceLevel::TheoremBacked {
        route.evidence = ExponentialEvidenceLevel::PaperProof;
        route.status = "paper_proof".into();
        route.statement.push_str(" The seed depends on an unreviewed paper argument; this closure does not upgrade its evidence level.");
    }
    route.premises.push(premise);
    Ok(Some(route))
}
fn nonnegative_series(spec: &BlockWeightSpec) -> Result<bool, ExponentialRowError> {
    Ok(match spec {
        BlockWeightSpec::SeriesExpression { expression } => {
            let parsed = crate::power_series::SeriesExpression::parse(expression)
                .map_err(|e| ExponentialRowError(e.to_string()))?;
            super::expression_routes::recognize_expression(&parsed)
                .map(|s| nonnegative_series(&s))
                .transpose()?
                .unwrap_or(false)
        }
        BlockWeightSpec::Finite { .. } | BlockWeightSpec::Sizes { .. } => nonnegative_finite(spec)?,
        BlockWeightSpec::Tail { .. } => true,
        BlockWeightSpec::Exponential { scale, rate, .. } => {
            !rational(scale)?.is_negative() && !rational(rate)?.is_negative()
        }
        BlockWeightSpec::PolynomialSubstitution { seed, coefficients } => {
            nonnegative_series(seed)?
                && coefficients
                    .iter()
                    .skip(1)
                    .all(|c| rational(c).is_ok_and(|x| !x.is_negative()))
        }
        BlockWeightSpec::RationalOdeSubstitution {
            seed,
            a,
            b,
            lambda,
            d,
        } => {
            nonnegative_series(seed)?
                && [a, b, lambda, d]
                    .iter()
                    .all(|s| rational(s).is_ok_and(|q| !q.is_negative()))
        }
        BlockWeightSpec::OdeSubstitution { seed, a, b, lambda } => {
            nonnegative_series(seed)?
                && [a, b, lambda]
                    .iter()
                    .all(|s| rational(s).is_ok_and(|x| !x.is_negative()))
        }
    })
}
fn prime(n: usize) -> bool {
    u64::try_from(n).is_ok_and(exact_algebra::primes::is_prime_u64)
}

fn unmarked_divisibility_violations(
    spec: &BlockWeightSpec,
    k: usize,
) -> Result<Option<Vec<String>>, ExponentialRowError> {
    if let Some(weights) = support(spec)? {
        if weights.iter().any(|(_, w)| w.is_negative()) {
            return Ok(None);
        }
        return Ok(Some(
            weights
                .iter()
                .filter(|(r, _)| k % r != 0)
                .map(|(r, _)| r.to_string())
                .collect(),
        ));
    }
    let minimum = match spec {
        BlockWeightSpec::Tail { minimum } => *minimum,
        BlockWeightSpec::Exponential {
            scale,
            rate,
            minimum,
        } => {
            let (scale, rate) = (rational(scale)?, rational(rate)?);
            if scale.is_zero() || rate.is_zero() {
                return Ok(Some(vec![]));
            }
            if scale.is_negative() || rate.is_negative() {
                return Ok(None);
            }
            *minimum
        }
        _ => return Ok(None),
    };
    // Every positive tail contains this size, which is greater than k.
    let witness = (num_bigint::BigInt::from(k) + 1u32).max(num_bigint::BigInt::from(minimum));
    Ok(Some(vec![witness.to_string()]))
}

/// Registry for exponential-series routes, separate from recurrence-record routes.
/// Recognizes full finite specs or structural ODE compositions, never prefixes.
/// The single-size filter is necessary only; passing it is not a proof.
/// In particular k=8 is under review and cannot issue a proved obstruction.
pub fn exponential_row_routes(
    family: &ExponentialFamily,
) -> Result<Vec<ExponentialRouteReport>, ExponentialRowError> {
    family.unmarked.validate(0)?;
    family.marked.validate(0)?;
    let normalized = ExponentialFamily {
        unmarked: normalize(&family.unmarked)?,
        marked: normalize(&family.marked)?,
    };
    let family = &normalized;
    let mut routes = Vec::new();
    if let Some(route) = substitution_route(family)? {
        routes.push(route);
    }
    if let Some(b) = support(&family.marked)? {
        if b.len() == 1 && b[0].1.is_positive() {
            let k = b[0].0;
            if prime(k) || (k % 2 == 0 && prime(k / 2)) || k == 8 {
                let Some(invalid) = unmarked_divisibility_violations(&family.unmarked, k)? else {
                    return Ok(routes);
                };
                routes.push(ExponentialRouteReport {
                    route:"single_marked_size_divisibility".into(),evidence:if k==8{ExponentialEvidenceLevel::UnderReview}else{ExponentialEvidenceLevel::PaperProof},
                    status:if k==8{"under_review"}else if invalid.is_empty(){"necessary_condition_passed"}else{"obstruction"}.into(),
                    statement:format!("Single positive marked size k={k}; nonnegative unmarked weights. Necessary paper condition: every unmarked size divides k. Violating sizes: {invalid:?}. Passing this filter does not prove all-rank real-rootedness. The k=8 necessity is under review and supplies no certified rejection."),references:vec![PAPER.into()],premises:vec![],consecutive_weak_interlacing:false,
                });
            }
        }
    }

    Ok(routes)
}

fn normalize(spec: &BlockWeightSpec) -> Result<BlockWeightSpec, ExponentialRowError> {
    let mut spec = spec.clone();
    match &mut spec {
        BlockWeightSpec::SeriesExpression { expression } => {
            let parsed = crate::power_series::SeriesExpression::parse(expression)
                .map_err(|e| ExponentialRowError(e.to_string()))?;
            if let Some(recognized) = super::expression_routes::recognize_expression(&parsed) {
                return Ok(recognized);
            }
        }
        BlockWeightSpec::PolynomialSubstitution { seed, .. }
        | BlockWeightSpec::OdeSubstitution { seed, .. }
        | BlockWeightSpec::RationalOdeSubstitution { seed, .. } => {
            **seed = normalize(seed)?;
        }
        _ => {}
    }
    Ok(spec)
}

/// Coefficient-sign hypotheses for one series. A nonnegative prefix alone
/// does not certify nonnegativity of the infinite series.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct CoefficientSignReport {
    /// `unmarked` for A or `marked` for B.
    pub series: String,
    /// `nonnegative` (structural proof), `negative` (witness), or `unknown`.
    pub status: String,
    /// Inclusive rank inspected exactly, even when only selected rows are checked.
    pub checked_through: usize,
    /// First negative coefficient in the inspected prefix, if present.
    pub first_negative_degree: Option<usize>,
    /// Ordinary coefficient of the negative witness, if present.
    pub first_negative_coefficient: Option<String>,
}

pub(super) fn coefficient_signs(
    family: &ExponentialFamily,
    n: usize,
) -> Result<Vec<CoefficientSignReport>, ExponentialRowError> {
    let mut out = Vec::new();
    for (name, spec) in [("unmarked", &family.unmarked), ("marked", &family.marked)] {
        let weights = spec.block_weights(n)?;
        let negative = weights.iter().enumerate().find(|(_, c)| c.is_negative());
        let proven = nonnegative_series(&normalize(spec)?)?;
        // Preserve reports for existing positive structured inputs; expression
        // hypotheses and any actual failed sign condition are explicit.
        if negative.is_some() || matches!(spec, BlockWeightSpec::SeriesExpression { .. }) {
            out.push(CoefficientSignReport {
                series: name.into(),
                status: if negative.is_some() {
                    "negative"
                } else if proven {
                    "nonnegative"
                } else {
                    "unknown"
                }
                .into(),
                checked_through: n,
                first_negative_degree: negative.map(|(r, _)| r),
                first_negative_coefficient: negative
                    .map(|(r, c)| (c / Q::from_integer(factorial_bigint(r))).to_string()),
            });
        }
    }
    Ok(out)
}

/// Optional factorization of a reference row for a one-free-factor pencil.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct PencilReference {
    /// Nonnegative constant of the removed linear factor.
    pub a: String,
    /// Nonnegative slope of the removed linear factor.
    pub b: String,
    /// Ascending rational coefficients of S in reference=(a+b*t)S.
    pub factor: Vec<String>,
}
/// Exact reference row for a per-rank multiplier or pencil certificate.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct FiniteRankRequest {
    /// Selected absolute rank; never an all-rank parameter.
    pub index: usize,
    /// Ascending positive rational reference coefficients, same length as target.
    pub reference: Vec<String>,
    /// When present use the pencil route; otherwise the multiplier-symbol route.
    pub pencil: Option<PencilReference>,
}
/// Exact finite-rank evidence with replay data, never an all-rank claim.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct FiniteRankCertificate {
    /// Embedded reference and optional pencil factorization.
    pub request: FiniteRankRequest,
    /// Exact ascending target coefficients.
    pub target: Vec<String>,
    /// Always finite_rank.
    pub evidence: ExponentialEvidenceLevel,
    /// multiplier_symbol or one_free_factor_pencil.
    pub route: String,
    /// True only when every stated hypothesis and exact test passes.
    pub certified: bool,
    /// Scope or failed sufficient condition.
    pub statement: String,
    /// Ratios target_i/reference_i, ascending.
    pub multipliers: Vec<String>,
    /// Binomially normalized multiplier symbol, or the two pencil images.
    pub witnesses: Vec<Vec<String>>,
    /// Public mathematical reference.
    pub reference: String,
}
impl FiniteRankCertificate {
    /// Recompute all exact hypotheses, identities and finite root checks.
    pub fn verify(&self) -> Result<(), ExponentialRowError> {
        let row = self
            .target
            .iter()
            .map(|s| rational(s))
            .collect::<Result<Vec<_>, _>>()?;
        if certify_finite_rank(&row, &self.request)? != *self {
            return Err(ExponentialRowError(
                "finite-rank certificate does not replay".into(),
            ));
        }
        Ok(())
    }
}
fn pf(row: &[Q]) -> bool {
    row.iter().all(|x| !x.is_negative())
        && is_real_rooted_bigint_coeffs(&clear_row_denominators(row).0)
}
/// Check a bounded-degree multiplier-ratio symbol, or the two images of S,tS.
/// All reference coefficients must be positive to define each ratio uniquely.
/// The symbol is `sum binomial(d,i)*gamma_i*t^i` (finite Polya--Schur).
/// For a pencil, check reference=(a+b*t)S exactly, a,b>=0, and directed weak
/// interlacing in either order of T(S), T(tS); this certifies every nonnegative
/// linear combination of those two images, at this rank only.
pub fn certify_finite_rank(
    target: &[Q],
    request: &FiniteRankRequest,
) -> Result<FiniteRankCertificate, ExponentialRowError> {
    let reference = request
        .reference
        .iter()
        .map(|s| rational(s))
        .collect::<Result<Vec<_>, _>>()?;
    if target.is_empty()
        || target.len() != reference.len()
        || reference.iter().any(|w| !w.is_positive())
    {
        return Err(ExponentialRowError(
            "reference must have positive coefficients and the target length".into(),
        ));
    }
    let gamma: Vec<_> = target.iter().zip(&reference).map(|(a, b)| a / b).collect();
    let d = target.len() - 1;
    let mut report = FiniteRankCertificate {
        request: request.clone(),
        target: strings(target),
        evidence: ExponentialEvidenceLevel::FiniteRank,
        route: "multiplier_symbol".into(),
        certified: false,
        statement: "Sufficient finite-rank test failed; no refutation follows.".into(),
        multipliers: strings(&gamma),
        witnesses: vec![],
        reference: "https://arxiv.org/abs/0809.0401".into(),
    };
    if let Some(pencil) = &request.pencil {
        report.route = "one_free_factor_pencil".into();
        let (a, b) = (rational(&pencil.a)?, rational(&pencil.b)?);
        let s = pencil
            .factor
            .iter()
            .map(|s| rational(s))
            .collect::<Result<Vec<_>, _>>()?;
        if s.len() != d
            || a.is_negative()
            || b.is_negative()
            || vec_poly::mul(&[a, b], &s) != reference
        {
            return Err(ExponentialRowError(
                "pencil requires exact reference=(a+b*t)S with a,b>=0 and deg S=d-1".into(),
            ));
        }
        let mut u = vec![Q::zero(); d + 1];
        let mut v = u.clone();
        for (i, c) in s.iter().enumerate() {
            u[i] = &gamma[i] * c;
            v[i + 1] = &gamma[i + 1] * c;
        }
        let ui = clear_row_denominators(&u).0;
        let vi = clear_row_denominators(&v).0;
        report.certified = pf(&reference)
            && pf(&u)
            && pf(&v)
            && (check_weak_interlacing_bigint_coeffs(&ui, &vi) == Some(true)
                || check_weak_interlacing_bigint_coeffs(&vi, &ui) == Some(true));
        report.witnesses = vec![strings(&u), strings(&v)];
        report.reference = "https://arxiv.org/abs/1410.6601 (section 8, compatibility)".into();
    } else {
        let symbol: Vec<_> = gamma
            .iter()
            .enumerate()
            .map(|(i, g)| g * Q::from_integer(binomial_bigint(d, i)))
            .collect();
        report.certified = pf(&reference) && pf(&symbol);
        report.witnesses = vec![strings(&symbol)];
    }
    if report.certified {
        report.statement="Exact certificate for the recorded rank only; target is real-rooted or zero. No all-rank conclusion.".into();
    }
    Ok(report)
}
