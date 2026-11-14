//! Exact certificates for stationary recurrences, separate from the numerical scout.
//!
//! [`certify_recurrence`] reduces the seed-dependent generating function before
//! deciding sparse pencils. A real limit set alone proves no all-row claim.
//! Standard quadratic seeds admit a public product proof; Toeplitz circle
//! certificates apply a theorem from an unpublished manuscript, formalized in
//! `RealRooted.Toeplitz.Circle.strictInterl_of_circleNonneg` and
//! `RealRooted.Toeplitz.Circle.roots_card_of_circleNonneg`. Every
//! verdict uses rational arithmetic and exact Sturm counts. Report replay
//! recomputes the decision from the embedded recurrence and options.
//!
//! For standard quadratic seeds, the identity
//! `P_n=A^(n mod 2) product_(1<=k<=n/2)(A^2+4cos^2(k*pi/(n+1))*B)`
//! gives an all-row proof from the pencil. For sparse order r, Tran's
//! parameter `(-1)^r*(-A)^r/(-B)` is `-A^r/B`, so the same pencil has
//! upper bound `r^r/(r-1)^(r-1)`, independent of the parity of r.
//! Sources: [Tran 2014, Theorems 1 and 3](https://arxiv.org/abs/1601.04382),
//! [Tran 2015, Theorem 1](https://arxiv.org/abs/1601.04383), and
//! [Sokal 2004, Theorem 1.5](https://arxiv.org/abs/cond-mat/0012369).

mod algebra;
mod pencil;
mod regions;
mod toeplitz;

pub use pencil::certify_real_rooted_pencil;
pub use regions::{
    certify_toeplitz_regions, RegionEvaluation, ToeplitzRegionCertificate,
    LIMIT_SPECTRUM_HYPOTHESIS, MAXIMUM_PRECOMPUTED_TOEPLITZ_BANDWIDTH,
};
pub use toeplitz::{CIRCLE_THEOREM_SOURCE, TOEPLITZ_CRITERIA};

use crate::recurrence::{
    parse_rational_coeff, rational_generating_function, RationalGeneratingFunction,
    RecurrenceError, RecurrenceJson,
};
use crate::{vec_poly, BigRational as Q, Polynomial};
use algebra::{counts, derivative_in_outer, resultant};
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};
type P = Polynomial<Q>;

fn q(n: i64) -> Q {
    Q::from_integer(n.into())
}
fn encode(p: &P) -> Vec<String> {
    p.coeffs().iter().map(ToString::to_string).collect()
}
fn decode(p: &[String]) -> Result<P, RecurrenceError> {
    Ok(P::new(
        p.iter()
            .map(|s| parse_rational_coeff(s).map_err(Into::into))
            .collect::<Result<_, RecurrenceError>>()?,
    ))
}

/// Exact classification of the finite BKW zero-limit set.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LimitSetStatus {
    /// Every finite limit point is real, including amplitude candidates.
    RealLimitSet,
    /// A certified nonreal limit point exists.
    NonrealLimitSet,
    /// A hypothesis of the ordinary BKW theorem fails; see the reason.
    Degenerate,
    /// The exact methods implemented here do not decide the limit set.
    Inconclusive,
}

/// Limits of finite row search and exact symbolic elimination.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CertificationOptions {
    /// Largest absolute row index checked, inclusive; default 200.
    pub maximum_row_index: usize,
    /// Largest polynomial degree admitted to symbolic elimination; default 24.
    pub maximum_pencil_degree: usize,
    /// Optional positive rational circle radius, checked before radius search.
    pub circle_radius: Option<String>,
}
impl Default for CertificationOptions {
    fn default() -> Self {
        Self {
            maximum_row_index: 200,
            maximum_pencil_degree: 24,
            circle_radius: None,
        }
    }
}

/// Exact root counts, with multiplicities removed by a polynomial gcd.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RootCounts {
    /// Ordinary degree; `None` denotes the zero polynomial.
    pub degree: Option<usize>,
    /// Degree after removing repeated factors.
    pub squarefree_degree: usize,
    /// Number of distinct real roots.
    pub distinct_real_roots: usize,
}
impl RootCounts {
    fn real_rooted(&self) -> bool {
        self.degree.is_some() && self.squarefree_degree == self.distinct_real_roots
    }
}

/// First finite nonzero row with nonreal zeros in the requested prefix.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FailingRow {
    /// Absolute recurrence index.
    pub index: usize,
    /// Exact ascending rational coefficients, including all zero low terms.
    pub coefficients: Vec<String>,
    /// Exact squarefree and real-root counts.
    pub root_counts: RootCounts,
}

/// A rational parameter sample from a discriminant cell.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PencilSample {
    /// Rational pencil parameter.
    pub parameter: String,
    /// Specialized polynomial, in ascending powers of x.
    pub polynomial: Vec<String>,
    /// Exact counts of the specialized polynomial.
    pub root_counts: RootCounts,
}

/// Exact decision of a polynomial pencil on a closed rational interval.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PencilCertificate {
    /// Constant polynomial U(x) of U+sV, ascending coefficients.
    pub constant: Vec<String>,
    /// Parameter coefficient V(x), ascending coefficients.
    pub slope: Vec<String>,
    /// Closed parameter bounds.
    pub interval: [String; 2],
    /// Parameter discriminant resultant times the leading coefficient, after
    /// removal of the parameter-independent gcd.
    pub critical_polynomial: Vec<String>,
    /// Disjoint rational enclosures of its roots strictly inside the interval.
    pub critical_intervals: Vec<[String; 2]>,
    /// A rational sample in each open cell and both rational endpoints.
    /// Algebraic boundaries inherit real-rootedness by closure from the cells.
    pub samples: Vec<PencilSample>,
    /// Whether every nonzero specialization is real-rooted.
    pub all_real_rooted: bool,
    /// Failing rational interior parameter, when present.
    pub witness: Option<PencilSample>,
}

/// Exact nonnegativity certificate applying the Toeplitz circle theorem.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CircleCertificate {
    /// The affine first coefficient alpha*x+beta, ascending.
    pub affine: Vec<String>,
    /// b_1,...,b_s in P_n=(alpha*x+beta)P_(n-1)-sum b_l P_(n-l-1).
    pub tails: Vec<String>,
    /// A positive rational radius.
    pub radius: String,
    /// Exact ascending coefficients of Q_R(u).
    pub polynomial: Vec<String>,
    /// Exact counts for Q_R; nonnegativity is checked on its sign cells.
    pub root_counts: RootCounts,
    /// Theorem source, explicitly identified as an unpublished manuscript.
    pub provenance: String,
    /// Fully qualified RealRooted declarations for the circle criterion.
    /// These are formal references; individual certificates are not replayed in Lean.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub lean_declarations: Vec<String>,
    /// Proved conclusion, including simplicity and strict interlacing.
    pub conclusion: String,
}

/// A nonreal sparse tie, after excluding zero last coefficients and amplitudes.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct NonrealTieWitness {
    /// Rational interior pencil parameter.
    pub parameter: String,
    /// Squarefree witness polynomial after removing factors shared with the
    /// last recurrence coefficient and the amplitude resultant.
    pub polynomial: Vec<String>,
    /// Exact counts showing a nonreal root survives those exclusions.
    pub root_counts: RootCounts,
}

/// Replayable exact result, with asymptotic and finite conclusions separated.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CertificationReport {
    /// Versioned wire schema, `polytool.bkw-certification.v1`.
    pub schema: String,
    /// Defining recurrence and every supplied seed.
    pub recurrence: RecurrenceJson,
    /// Options used for deterministic replay.
    pub options: CertificationOptions,
    /// Exact limit-set status; never inferred from a numerical scout.
    pub limit_set_status: LimitSetStatus,
    /// `proved`, `conditional`, `refuted`, or `inconclusive` for the all-row claim.
    pub result: String,
    /// Explicit named mathematical hypotheses, present only for `conditional`.
    /// These are assumptions, never unconditional proofs or attestations.
    #[serde(default)]
    pub conditional_on: Vec<String>,
    /// Explanation of hypotheses, degeneracies or coverage limits.
    pub reason: String,
    /// Reduced seed-dependent generating function; absent for unsupported rules.
    pub generating_function: Option<RationalGeneratingFunction>,
    /// Exact pencil certificate, if the shape permits it.
    pub pencil: Option<PencilCertificate>,
    /// Certified nonreal tie, with exceptional roots excluded explicitly.
    pub nonreal_tie: Option<NonrealTieWitness>,
    /// Exact circle certificate applying the unpublished manuscript theorem.
    pub circle: Option<CircleCertificate>,
    /// Independent exact limit-spectrum and circle-region membership.
    pub toeplitz_region: Option<ToeplitzRegionCertificate>,
    /// Existing product-induction route, reused for a single standard branch.
    pub product_proof: Option<Box<crate::recurrence::RealRootednessTheoremReport>>,
    /// Proof route used, independent of the limit-set status; see its provenance.
    pub proof_route: Option<String>,
    /// Exact resultant Res_t(N,D), bounding amplitude-zero candidates.
    pub amplitude_resultant: Option<Vec<String>>,
    /// First finite nonreal row, kept separate from a limit-set witness.
    pub first_failing_row: Option<FailingRow>,
    /// Last absolute row checked, including supplied seeds.
    pub checked_through: Option<usize>,
    /// Zero rows observed; a zero polynomial is never counted as a failing row.
    pub zero_rows: Vec<usize>,
}

impl CertificationReport {
    /// Recompute the embedded exact decision and reject any changed field.
    ///
    /// # Errors
    /// Invalid embedded input or a mismatch in any certificate field.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        let expected = certify_recurrence(&self.recurrence, &self.options)?;
        if expected != *self {
            return Err(RecurrenceError::VerificationFailed(
                "BKW certificate does not replay".into(),
            ));
        }
        Ok(())
    }
}

/// Certify a stationary recurrence and search its exact rows for a refutation.
///
/// Standard seeds mean the reduced generating function is exactly 1/D and
/// the original seeds agree with its standard expansion (including every
/// extra prescribed row). A shifted first index merely labels this expansion.
/// A fitted rule defines a family here; identifying it with another infinite
/// family requires a separate all-index identity.
///
/// # Errors
/// Malformed recurrence/seed data, invalid radius, or row evaluation failure.
pub fn certify_recurrence(
    record: &RecurrenceJson,
    options: &CertificationOptions,
) -> Result<CertificationReport, RecurrenceError> {
    let (rule, first, seeds) = record.to_recurrence_parts()?;
    if let Some(r) = &options.circle_radius {
        if parse_rational_coeff(r).map_err(RecurrenceError::from)? <= Q::zero() {
            return Err(RecurrenceError::InvalidInput(
                "circle radius must be positive".into(),
            ));
        }
    }
    let mut report = CertificationReport {
        schema: "polytool.bkw-certification.v1".into(),
        recurrence: record.clone(),
        options: options.clone(),
        limit_set_status: LimitSetStatus::Inconclusive,
        result: "inconclusive".into(),
        conditional_on: Vec::new(),
        reason: String::new(),
        generating_function: None,
        pencil: None,
        nonreal_tie: None,
        circle: None,
        toeplitz_region: None,
        product_proof: None,
        proof_route: None,
        amplitude_resultant: None,
        first_failing_row: None,
        checked_through: None,
        zero_rows: Vec::new(),
    };
    let gf = match rational_generating_function(&rule, first, &seeds) {
        Ok(gf) => gf,
        Err(RecurrenceError::Unsupported(reason)) => {
            report.reason = reason;
            return Ok(report);
        }
        Err(e) => return Err(e),
    };
    let d = gf
        .denominator
        .iter()
        .map(|p| decode(p))
        .collect::<Result<Vec<_>, _>>()?;
    let n = gf
        .numerator
        .iter()
        .map(|p| decode(p))
        .collect::<Result<Vec<_>, _>>()?;
    report.generating_function = Some(gf);
    let order = d.len().saturating_sub(1);
    let standard = n.len() == 1 && n[0] == P::one();
    let degenerate = if n.is_empty() {
        Some("identically zero sequence after exact numerator cancellation")
    } else if order == 0 {
        Some("finite polynomial generating function: every later row is zero")
    } else if order == 1 {
        Some("one active branch after cancellation; product family outside the nondegenerate BKW theorem")
    } else if scalar_monomial_symbol(&d) {
        Some("scalar-monomial characteristic collapse: constant root ratios after factoring a common monomial scale")
    } else if order > options.maximum_pencil_degree {
        None
    } else if resultant(&d, &derivative_in_outer(&d)).is_zero() {
        Some("repeated active poles in the reduced generating function")
    } else {
        None
    };
    if let Some(reason) = degenerate {
        report.limit_set_status = LimitSetStatus::Degenerate;
        report.reason = reason.into();
        if order == 1 && standard {
            // Reuse the existing all-index product route; do not duplicate its
            // multiplier/root-domain or extra-seed checks.
            let multiplier = -d[1].clone();
            let reduced_rule = crate::recurrence::Recurrence {
                terms: vec![crate::recurrence::RecurrenceTerm::new(
                    1,
                    0,
                    crate::recurrence::BivarPoly::from_coefficients_normalized(vec![multiplier
                        .coeffs()
                        .to_vec()]),
                )],
                denominator: None,
                inhomogeneous: None,
            };
            let proof = crate::recurrence::check_product_recurrence_theorem(
                &reduced_rule,
                first,
                &seeds,
                &Default::default(),
            )?;
            if proof.is_proved() {
                report.result = "proved".into();
                report.proof_route = Some("product_recurrence".into());
                report.limit_set_status = LimitSetStatus::RealLimitSet;
                report.reason="single standard active branch handled by the existing exact product-induction route".into();
                report.product_proof = Some(Box::new(proof));
            }
        }
    } else if order <= options.maximum_pencil_degree {
        let active: Vec<_> = d
            .iter()
            .enumerate()
            .skip(1)
            .filter(|(_, p)| !p.is_zero())
            .map(|(j, _)| j)
            .collect();
        let gcd = active.iter().copied().reduce(gcd_usize).unwrap_or(1);
        if gcd > 1 {
            report.limit_set_status = LimitSetStatus::Degenerate;
            report.reason=format!("constant unit-phase ratios: split into residue classes modulo {gcd} before applying BKW");
        } else {
            let amplitude = resultant(&n, &d);
            report.amplitude_resultant = Some(encode(&amplitude));
            if order >= 2 && active.iter().all(|j| *j == 1 || *j == order) && !d[1].is_zero() {
                let a = -d[1].clone();
                let b = -d[order].clone();
                let u = a.pow(order);
                if (u.clone() * b.derivative() - u.derivative() * b.clone()).is_zero() {
                    report.limit_set_status = LimitSetStatus::Degenerate;
                    report.reason =
                        "constant characteristic-root ratios (including scalar-monomial collapses)"
                            .into();
                } else if u.degree().unwrap_or(0).max(b.degree().unwrap_or(0))
                    <= options.maximum_pencil_degree
                {
                    let upper = Q::new(
                        num_bigint::BigInt::from(order).pow(order as u32),
                        num_bigint::BigInt::from(order - 1).pow((order - 1) as u32),
                    );
                    let pencil = pencil::decide(&u, &b, &upper)?;
                    if let Some(witness) = &pencil.witness {
                        let w = decode(&witness.polynomial)?.squarefree_part();
                        let exceptional = amplitude.clone() * b.clone();
                        let away = w.exact_div(&w.gcd(&exceptional));
                        if !counts(&away).real_rooted() && !away.is_zero() {
                            report.limit_set_status = LimitSetStatus::NonrealLimitSet;
                            report.nonreal_tie = Some(NonrealTieWitness {
                                parameter: witness.parameter.clone(),
                                polynomial: encode(&away),
                                root_counts: counts(&away),
                            });
                            report.reason=format!("sparse order {order}: rational interior pencil witness, nonconstant root ratios, active simple branches; amplitude zeros removed by gcd (Tran 2014/2015; Sokal 2004 Theorem 1.5)");
                        }
                    } else if pencil.all_real_rooted
                        && !amplitude.is_zero()
                        && counts(&amplitude).real_rooted()
                    {
                        report.limit_set_status = LimitSetStatus::RealLimitSet;
                        report.reason="sparse tie pencil is real-rooted throughout; all possible isolated amplitude zeros are real".into();
                    }
                    if order == 2 && standard && pencil.all_real_rooted && counts(&a).real_rooted()
                    {
                        report.result = "proved".into();
                        report.proof_route = Some("quadratic_standard_seed_product".into());
                        report.reason="standard quadratic generating function: exact pencil decision and Tran 2014 product formula prove every row real-rooted".into();
                    }
                    report.pencil = Some(pencil);
                }
            }
            if standard && report.result != "proved" {
                report.circle = toeplitz::recognize_and_certify(&d, options)?;
                if report.circle.is_some() {
                    report.result = "proved".into();
                    report.proof_route = Some("hessenberg_toeplitz_circle".into());
                    report.reason = CIRCLE_THEOREM_SOURCE.into();
                }
            }
        }
    }
    if standard {
        if let Some(tails) = toeplitz::hessenberg_tails(&d) {
            report.toeplitz_region = certify_toeplitz_regions(&tails)?;
            if let Some(geometry) = &report.toeplitz_region {
                if let Some(real) = geometry.limit_real {
                    if !real && report.result == "proved" {
                        return Err(RecurrenceError::VerificationFailed(
                            "Toeplitz region contradicts an all-row proof".into(),
                        ));
                    }
                    report.limit_set_status = if real {
                        LimitSetStatus::RealLimitSet
                    } else {
                        LimitSetStatus::NonrealLimitSet
                    };
                    if real && report.result != "proved" && report.circle.is_none() {
                        report.result = "conditional".into();
                        report.conditional_on = vec![LIMIT_SPECTRUM_HYPOTHESIS.into()];
                        report.reason="Schmidt-Spitzer region membership certifies a real limit spectrum exactly; the all-row conclusion assumes hessenberg_real_limit_spectrum".into();
                    }
                }
            } else if report.result == "conditional" {
                report.result = "inconclusive".into();
                report.conditional_on.clear();
                report.reason="Toeplitz bandwidth is outside the precomputed range; no conditional all-row verdict is made".into();
            }
        }
    }
    // A complete proof already covers every row. Otherwise keep finite row
    // search independent of the asymptotic decision, with original seeds.
    if report.result != "proved" {
        let maximum = options.maximum_row_index;
        if maximum >= first {
            let count = maximum
                .checked_sub(first)
                .and_then(|n| n.checked_add(1))
                .ok_or_else(|| RecurrenceError::InvalidInput("row count overflow".into()))?;
            for (position, row) in rule
                .generate_rows_rational_iter(&seeds, first)
                .take(count)
                .enumerate()
            {
                let index = first + position;
                let p = P::new(row.map_err(|e| RecurrenceError::InvalidInput(e.to_string()))?);
                report.checked_through = Some(index);
                if p.is_zero() {
                    report.zero_rows.push(index);
                    continue;
                }
                if !crate::real_rootedness::is_real_rooted_bigint_coeffs(
                    &vec_poly::clear_denominators(p.coeffs()),
                ) {
                    report.first_failing_row = Some(FailingRow {
                        index,
                        coefficients: encode(&p),
                        root_counts: counts(&p),
                    });
                    break;
                }
            }
        }
    }
    if report.first_failing_row.is_some()
        || report.limit_set_status == LimitSetStatus::NonrealLimitSet
    {
        report.result = "refuted".into();
        report.proof_route = None;
        report.circle = None;
        report.product_proof = None;
        report.conditional_on.clear();
    }
    if report.reason.is_empty() {
        report.reason="general-order or unsupported amplitude geometry: exact limit-set decision inconclusive; finite rows are separate evidence; use the floating scout only as evidence".into();
    }
    Ok(report)
}

fn gcd_usize(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn scalar_monomial_symbol(d: &[P]) -> bool {
    // F(x,lambda)=x^(kr) f(lambda/x^k), including k=0.
    // This is a degeneracy test, not a BKW verdict on the whole plane.
    let Some((j, coefficient)) = d.iter().enumerate().skip(1).find(|(_, p)| !p.is_zero()) else {
        return false;
    };
    let degree = coefficient.degree().unwrap_or(0);
    if degree % j != 0 {
        return false;
    }
    let scale = degree / j;
    d.iter().enumerate().skip(1).all(|(j, p)| {
        p.is_zero()
            || (p.degree() == Some(j * scale)
                && p.coeffs().iter().filter(|c| !c.is_zero()).count() == 1)
    })
}
