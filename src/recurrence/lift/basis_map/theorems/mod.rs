//! [`check_basis_map_preserves_pf`]: the basis-map theorems with their
//! hypotheses decided exactly (module documentation of `basis_map`,
//! "Theorems").

use super::map::BasisMap;
use super::published::{recognize_published_family, PublishedFamily};
use super::quadratic::{integer, AlgebraicInterval, QuadraticNumber};
use crate::recurrence::implicit::algebra::{Biv, BivarConversion};
use crate::recurrence::{
    format_rational_coeff, LiftHypothesis, LiftHypothesisStatus, RecurrenceError, RecurrenceJson,
    LIFT_LEAN_REPLAY_NOT_ATTEMPTED,
};
use crate::vec_poly;
use crate::BigRational;
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;

type Q = BigRational;

/// JSON schema tag of [`BasisMapTheoremReport`].
pub const BASIS_MAP_THEOREM_SCHEMA: &str = "polytool.recurrence-basis-map-theorem.v1";

/// Options of [`check_basis_map_preserves_pf`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BasisMapTheoremOptions {
    /// Largest index searched for an explicit witness of a failure that a
    /// theorem guarantees (for example the first monomial `u^n` whose image
    /// has a negative constant term).  The failure is reported from the
    /// theorem when no witness is found within the limit.  Default 2000.
    /// The even-derivative shape caps this search at degree 16; its finite
    /// signed-symbol decision does not require a monomial witness.
    pub witness_search_limit: usize,
}

impl Default for BasisMapTheoremOptions {
    fn default() -> Self {
        Self {
            witness_search_limit: 2000,
        }
    }
}

/// The recognized shape of a basis family, with its exact parameters
/// (decimal strings), after dividing every row by the constant
/// `normalization` where there is one.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum BasisMapShape {
    /// `P_k = g(t) (a + d t)^k`.
    Geometric {
        /// `g = P_0`, ascending.
        initial_row: Vec<String>,
        /// `a`.
        constant: String,
        /// `d`.
        slope: String,
    },
    /// `P_0 = C`, `P_1 = C s`, `P_(k+2) = a P_(k+1) + (b t + c) P_k`.
    ConstantCurrent {
        /// `C`.
        normalization: String,
        /// `a`.
        a: String,
        /// `b`.
        b: String,
        /// `c`.
        c: String,
        /// `s`.
        initial_constant: String,
    },
    /// `P_0 = C`, `P_1 = C (l + m t)`,
    /// `P_(k+2) = (a + d t) P_(k+1) + (b t + c) P_k`, `d != 0`.
    AffineCurrent {
        /// `C`.
        normalization: String,
        /// `a`.
        a: String,
        /// `d`.
        d: String,
        /// `b`.
        b: String,
        /// `c`.
        c: String,
        /// `l`.
        initial_constant: String,
        /// `m`.
        initial_slope: String,
    },
    /// `P_(k+1) = (a + b t) P_k + kappa (t - h) P_k'`, `kappa != 0`.
    RayDifferential {
        /// `g = P_0`, ascending.
        initial_row: Vec<String>,
        /// `a`.
        a: String,
        /// `b`.
        b: String,
        /// `kappa`.
        kappa: String,
        /// `h`.
        root: String,
    },
    /// `P_(k+1) = (a + b t) P_k + kappa (t - l)(t - h) P_k'`, `l < h`.
    IntervalDifferential {
        /// `g = P_0`, ascending.
        initial_row: Vec<String>,
        /// `a`.
        a: String,
        /// `b`.
        b: String,
        /// `kappa`.
        kappa: String,
        /// `l`.
        lower_root: String,
        /// `h`.
        upper_root: String,
    },
    /// `P_(k+1) = (a + b t) P_k - kappa P_k'` (a Hermite-type operator).
    HermiteDifferential {
        /// `g = P_0`, ascending.
        initial_row: Vec<String>,
        /// `a`.
        a: String,
        /// `b`.
        b: String,
        /// `kappa`.
        kappa: String,
    },
    /// `P_(k+1) = (alpha(k) + b t) P_k + kappa (t - h) P_k'` with an
    /// index-dependent scalar `alpha(k)` (a polynomial in the basis index
    /// `k`), `kappa = 0` allowed (a Newton basis when `P_0` is constant).
    ScalarIndex {
        /// `g = P_0`, ascending.
        initial_row: Vec<String>,
        /// `alpha(k)`, ascending in `k`.
        scalar: Vec<String>,
        /// `b`.
        b: String,
        /// `kappa`.
        kappa: String,
        /// `h` (`0` when `kappa = 0`).
        root: String,
    },
    /// `P_k = P_K (a + d t)^(k - K)` from `K >= 1` on, with an exceptional
    /// prefix.
    EventuallyGeometric {
        /// `K`.
        geometric_from: usize,
        /// `a`.
        constant: String,
        /// `d`.
        slope: String,
    },
    /// `u^N -> sum_k c_k (N)_(2k) x^(N-k)`, with replayable all-degree hypotheses.
    EvenDerivative {
        /// Coefficients, exact sign decision, basis identity and row PF proof.
        certificate: Box<super::EvenDerivativeCertificate>,
    },
    /// None of the shapes above.
    Unrecognized {
        /// Why.
        reason: String,
    },
    /// The Narayana transformation `P_k = N_(k,m)(t) = sum_j C(k,j)
    /// C(k+m,j) / C(m+j,j) t^j` (Mao and Wang), recognised for every `k`.
    Narayana {
        /// `m >= 0`.
        parameter: usize,
    },
    /// The Eulerian transformation `A°`: `P_0 = 1`, `P_k = t A_k(t)`
    /// (Athanasiadis), recognised for every `k`.
    EulerianTransformation,
    /// The derangement transformation `P_k = d_k(t)` (Brändén and Solus),
    /// recognised for every `k`.
    DerangementTransformation,
}

/// Outcome of [`check_basis_map_preserves_pf`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BasisMapTheoremStatus {
    /// Every hypothesis of the reported conclusions is verified.
    Proved,
    /// The map provably does not preserve real-rootedness on any interval
    /// of the stated kind (an `iff` theorem or an exact counterexample).
    NotPreserving,
    /// A hypothesis fails or remains unresolved; the theorem gives no conclusion
    /// (nothing is claimed against the map).
    HypothesisFailed,
    /// The family has none of the recognized shapes.
    NotApplicable,
}

/// One proved statement about a basis map: every nonzero real polynomial
/// whose roots all lie in `input_interval` (constants included) is sent to
/// zero or to a real-rooted polynomial whose roots all lie in
/// `output_interval`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BasisMapConclusion {
    /// The source interval.
    pub input_interval: AlgebraicInterval,
    /// The target interval.
    pub output_interval: AlgebraicInterval,
    /// Whether every input with positive leading coefficient and roots in
    /// `input_interval` has an image with nonnegative coefficients (with the
    /// target in `(-inf, 0]`: PF inputs go to PF outputs, zero allowed).
    pub nonnegative_coefficients: bool,
    /// Whether no strictly larger interval works for the same target (an
    /// `iff` theorem); `false` means not claimed.
    pub sharp: bool,
    /// The statement in words.
    pub statement: String,
    /// The names of the report's hypotheses this statement rests on (each
    /// verified).
    pub hypotheses: Vec<String>,
}

impl BasisMapConclusion {
    /// Whether this conclusion maps PF inputs to PF outputs: the target lies
    /// in `(-inf, 0]`, images of positive-leading inputs have nonnegative
    /// coefficients, and the source meets `(-inf, 0]` in an interval of
    /// positive length.
    pub fn is_pf_conclusion(&self) -> bool {
        self.output_interval.is_nonpositive()
            && self.nonnegative_coefficients
            && !self
                .input_interval
                .intersection(&AlgebraicInterval::nonpositive_ray())
                .is_degenerate()
    }

    /// `input_interval` intersected with `(-inf, 0]`.
    pub fn nonpositive_input_interval(&self) -> AlgebraicInterval {
        self.input_interval
            .intersection(&AlgebraicInterval::nonpositive_ray())
    }
}

/// An exact counterexample: a nonzero input whose image is not
/// real-rooted, or not in the claimed region.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BasisMapCounterexample {
    /// The input, ascending exact coefficients.
    pub input: Vec<String>,
    /// Its image, ascending exact coefficients.
    pub image: Vec<String>,
    /// What fails (decided exactly).
    pub failure: String,
}

/// The result of [`check_basis_map_preserves_pf`] (schema
/// `polytool.recurrence-basis-map-theorem.v1`).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BasisMapTheoremReport {
    /// [`BASIS_MAP_THEOREM_SCHEMA`].
    pub schema: String,
    /// The map's name.
    pub basis_name: String,
    /// The basis family (rule, first index, minimal seeds).
    pub basis: RecurrenceJson,
    /// The recognized shape.
    pub shape: BasisMapShape,
    /// The outcome.
    pub status: BasisMapTheoremStatus,
    /// The theorem applied, with its reference.
    pub theorem: String,
    /// The hypotheses in the order checked (exact rational decisions).
    pub hypotheses: Vec<LiftHypothesis>,
    /// The proved statements (empty unless `status` is `proved`).
    pub conclusions: Vec<BasisMapConclusion>,
    /// An exact counterexample, when one was found.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub counterexample: Option<BasisMapCounterexample>,
    /// Fully qualified declarations of the Lean library RealRooted
    /// (<https://github.com/PerAlexandersson/RealRooted>) behind the
    /// theorem; empty for elementary or paper-only statements.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub lean_declarations: Vec<String>,
    /// Always `"not_attempted"`.
    pub lean_replay: String,
}

impl BasisMapTheoremReport {
    /// Whether the status is [`BasisMapTheoremStatus::Proved`].
    pub fn is_proved(&self) -> bool {
        self.status == BasisMapTheoremStatus::Proved
    }

    /// The proved PF statement with the largest nonpositive source, if any
    /// ([`BasisMapConclusion::is_pf_conclusion`]).
    pub fn pf_conclusion(&self) -> Option<&BasisMapConclusion> {
        if !self.is_proved() {
            return None;
        }
        let mut best: Option<&BasisMapConclusion> = None;
        for conclusion in self.conclusions.iter().filter(|c| c.is_pf_conclusion()) {
            let better = best.is_none_or(|b| {
                b.nonpositive_input_interval()
                    .is_subset_of(&conclusion.nonpositive_input_interval())
                    && conclusion.output_interval.is_subset_of(&b.output_interval)
            });
            if better {
                best = Some(conclusion);
            }
        }
        best
    }

    /// Recompute the check from the recorded basis family and compare.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::VerificationFailed`] when the recomputed report
    /// differs, and the errors of the check.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        if self.schema != BASIS_MAP_THEOREM_SCHEMA {
            return Err(RecurrenceError::VerificationFailed(format!(
                "unsupported basis-map theorem schema `{}`",
                self.schema
            )));
        }
        let (rule, first_index, seeds) = self.basis.to_recurrence_parts()?;
        let map = BasisMap::new(&self.basis_name, rule, first_index, seeds)?;
        // A report need not search for additional witnesses when replaying:
        // the absence of a recorded counterexample is not a mathematical claim.
        let options = if matches!(self.shape, BasisMapShape::EvenDerivative { .. }) {
            BasisMapTheoremOptions {
                witness_search_limit: self.counterexample.as_ref().map_or(0, |c| c.input.len()),
            }
        } else {
            BasisMapTheoremOptions::default()
        };
        let recomputed = check_basis_map_preserves_pf(&map, &options)?;
        let comparable = |r: &Self| {
            (
                r.shape.clone(),
                r.status,
                r.hypotheses.clone(),
                r.conclusions.clone(),
            )
        };
        if matches!(self.shape, BasisMapShape::EvenDerivative { .. })
            && (self.theorem != recomputed.theorem
                || self.lean_declarations != recomputed.lean_declarations
                || self.lean_replay != recomputed.lean_replay
                || self.counterexample != recomputed.counterexample)
        {
            return Err(RecurrenceError::VerificationFailed(
                "even-derivative theorem metadata differs".into(),
            ));
        }
        if comparable(&recomputed) != comparable(self) {
            return Err(RecurrenceError::VerificationFailed(
                "the recorded basis-map theorem report differs from the recomputed one".into(),
            ));
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn text(value: &Q) -> String {
    format_rational_coeff(value)
}

fn poly_text(p: &[Q]) -> Vec<String> {
    p.iter().map(format_rational_coeff).collect()
}

fn hypothesis(name: &str, statement: String, holds: bool, evidence: String) -> LiftHypothesis {
    LiftHypothesis {
        name: name.into(),
        statement,
        status: if holds {
            LiftHypothesisStatus::Verified
        } else {
            LiftHypothesisStatus::Failed { witness: None }
        },
        evidence,
        sign_decision: None,
    }
}

fn rational(value: &Q) -> QuadraticNumber {
    QuadraticNumber::rational(value.clone())
}

/// Multiplicity of the rational `h` as a root of the nonzero `g`.
fn multiplicity(g: &[Q], h: &Q) -> usize {
    let mut p = vec_poly::trim(g);
    let mut m = 0;
    while vec_poly::degree(&p).unwrap_or(0) > 0 && vec_poly::evaluate(&p, h).is_zero() {
        p = vec_poly::exact_div_by_t_minus(&p, h).expect("h is a root");
        m += 1;
    }
    m
}

fn leading(p: &[Q]) -> Q {
    vec_poly::trim(p).last().cloned().unwrap_or_else(Q::zero)
}

fn is_constant(p: &[Q]) -> bool {
    vec_poly::degree(&vec_poly::trim(p)).unwrap_or(0) == 0
}

/// The coefficients of `t^0, t^1, ...` of an `n`-free bivariate polynomial.
fn t_polynomial(p: &Biv) -> Vec<Q> {
    vec_poly::trim(p.rows().first().map(Vec::as_slice).unwrap_or(&[]))
}

/// Outcome of the analysis, before assembly.
struct Analysis {
    shape: BasisMapShape,
    status: BasisMapTheoremStatus,
    theorem: String,
    hypotheses: Vec<LiftHypothesis>,
    conclusions: Vec<BasisMapConclusion>,
    counterexample: Option<BasisMapCounterexample>,
    lean: Vec<&'static str>,
}

impl Analysis {
    fn not_applicable(reason: String) -> Self {
        Self {
            shape: BasisMapShape::Unrecognized { reason },
            status: BasisMapTheoremStatus::NotApplicable,
            theorem: String::new(),
            hypotheses: Vec::new(),
            conclusions: Vec::new(),
            counterexample: None,
            lean: Vec::new(),
        }
    }

    /// Record a conclusion resting on every hypothesis verified so far.
    fn push_conclusion(&mut self, mut conclusion: BasisMapConclusion) {
        conclusion.hypotheses = self
            .hypotheses
            .iter()
            .filter(|h| h.status == LiftHypothesisStatus::Verified)
            .map(|h| h.name.clone())
            .collect();
        self.conclusions.push(conclusion);
    }

    /// Status from the hypotheses and conclusions gathered.
    fn settle(mut self) -> Self {
        if matches!(
            self.status,
            BasisMapTheoremStatus::NotPreserving | BasisMapTheoremStatus::NotApplicable
        ) {
            self.conclusions.clear();
            return self;
        }
        if self.conclusions.is_empty() {
            self.status = BasisMapTheoremStatus::HypothesisFailed;
        } else {
            self.status = BasisMapTheoremStatus::Proved;
        }
        self
    }
}

/// The counterexample `input -> image` with the failure description.
fn counterexample(input: &[Q], image: &[Q], failure: String) -> BasisMapCounterexample {
    BasisMapCounterexample {
        input: poly_text(input),
        image: poly_text(image),
        failure,
    }
}

/// Real-rootedness of `g` (nonzero) decided exactly.
fn is_real_rooted(g: &[Q]) -> bool {
    AlgebraicInterval::whole_line().contains_all_roots_of(g)
}

// ---------------------------------------------------------------------------
// The check
// ---------------------------------------------------------------------------

/// Decide the hypotheses of the basis-map theorems for `map` exactly and
/// report the proved source and target intervals (module documentation of
/// `basis_map`, "Theorems", for every shape, its theorem and its proof or
/// reference).
///
/// A [`BasisMapTheoremStatus::Proved`] report lists statements of the form
/// "every nonzero real polynomial with all roots in `S` is sent to zero or
/// to a real-rooted polynomial with all roots in `J`", for all input
/// degrees; [`BasisMapTheoremReport::pf_conclusion`] picks the statement
/// about PF inputs (nonnegative coefficients, real nonpositive roots).
/// Every hypothesis is an exact rational computation (signs, exact roots of
/// fixed polynomials, root multiplicities); the all-degree statement is the
/// theorem's.  Other statuses never claim anything except
/// [`BasisMapTheoremStatus::NotPreserving`], which rests on an exact
/// counterexample or on an `iff` theorem.
///
/// # Errors
///
/// [`RecurrenceError`] when the basis rows cannot be generated.
///
/// # Examples
///
/// ```
/// use polytool::recurrence::*;
///
/// // A035607: P_0 = 1, P_1 = 1 + 2t, P_(n+2) = (1 + t) P_(n+1) + t P_n.
/// let map = BasisMap::from_oeis("A035607").unwrap();
/// let report = check_basis_map_preserves_pf(&map, &BasisMapTheoremOptions::default()).unwrap();
/// let pf = report.pf_conclusion().unwrap();
/// assert_eq!(pf.input_interval.to_string(), "[-1 - sqrt(2), 0]");
/// ```
pub fn check_basis_map_preserves_pf(
    map: &BasisMap,
    options: &BasisMapTheoremOptions,
) -> Result<BasisMapTheoremReport, RecurrenceError> {
    let analysis = analyze(map, options)?.settle();
    Ok(BasisMapTheoremReport {
        schema: BASIS_MAP_THEOREM_SCHEMA.into(),
        basis_name: map.name().to_string(),
        basis: map.recurrence_json(),
        shape: analysis.shape,
        status: analysis.status,
        theorem: analysis.theorem,
        hypotheses: analysis.hypotheses,
        conclusions: analysis.conclusions,
        counterexample: analysis.counterexample,
        lean_declarations: analysis.lean.iter().map(|s| s.to_string()).collect(),
        lean_replay: LIFT_LEAN_REPLAY_NOT_ATTEMPTED.into(),
    })
}

fn analyze(map: &BasisMap, options: &BasisMapTheoremOptions) -> Result<Analysis, RecurrenceError> {
    let order = map.order();
    if order == 0 {
        return Ok(Analysis::not_applicable("the rule has no terms".into()));
    }
    if !map.is_index_independent() {
        if let Some(analysis) = published(map)? {
            return Ok(analysis);
        }
        if let Some(analysis) = even_derivative(map, options)? {
            return Ok(analysis);
        }
        return scalar_index(map, options);
    }
    if order >= 3 {
        if let Some(analysis) = even_derivative(map, options)? {
            return Ok(analysis);
        }
    }
    let leader = match &map.rule().denominator {
        None => Q::one(),
        Some(f) => match Biv::from_bivar(f).as_constant() {
            Some(c) => c,
            None => {
                return Ok(Analysis::not_applicable(
                    "the leader of the rule depends on t".into(),
                ))
            }
        },
    };
    // Shift-free operators A_j with P_(k+R) = sum_j A_j P_(k+j), for k >= K - R.
    let mut operators: Vec<std::collections::BTreeMap<usize, Vec<Q>>> =
        vec![Default::default(); order];
    for term in &map.rule().terms {
        let c = vec_poly::scale(
            &t_polynomial(&Biv::from_bivar(&term.coeff)),
            &leader.recip(),
        );
        let slot = operators[order - term.offset]
            .entry(term.deriv_order)
            .or_insert_with(|| vec![Q::zero()]);
        *slot = vec_poly::trim(&vec_poly::add(slot, &c));
    }
    for operator in &mut operators {
        operator.retain(|_, c| !vec_poly::is_zero(c));
    }
    let seeds = map.seeds();
    match order {
        1 => {
            let a0 = &operators[0];
            let max_derivative = a0.keys().max().copied().unwrap_or(0);
            let multiplier = a0.get(&0).cloned().unwrap_or_else(|| vec![Q::zero()]);
            if max_derivative == 0 {
                if vec_poly::degree(&multiplier).unwrap_or(0) > 1 {
                    return Ok(Analysis::not_applicable(
                        "a first-order rule whose multiplier has degree above one".into(),
                    ));
                }
                if seeds.len() == 1 {
                    return geometric(&seeds[0], &multiplier);
                }
                return eventually_geometric(map, &multiplier);
            }
            if seeds.len() != 1 {
                return Ok(Analysis::not_applicable(
                    "a differential rule with an exceptional prefix".into(),
                ));
            }
            if max_derivative > 1 || vec_poly::degree(&multiplier).unwrap_or(0) > 1 {
                return Ok(Analysis::not_applicable(
                    "a first-order rule outside (a + b t) + v(t) D".into(),
                ));
            }
            let v = a0.get(&1).cloned().unwrap_or_default();
            differential(map, &seeds[0], &multiplier, &v, options)
        }
        2 => {
            if operators.iter().any(|op| op.keys().any(|d| *d > 0)) {
                return Ok(Analysis::not_applicable(
                    "an order-two rule with derivatives".into(),
                ));
            }
            let current = operators[1]
                .get(&0)
                .cloned()
                .unwrap_or_else(|| vec![Q::zero()]);
            let lag = operators[0]
                .get(&0)
                .cloned()
                .unwrap_or_else(|| vec![Q::zero()]);
            if vec_poly::degree(&current).unwrap_or(0) > 1
                || vec_poly::degree(&lag).unwrap_or(0) > 1
            {
                return Ok(Analysis::not_applicable(
                    "an order-two rule with a multiplier of degree above one".into(),
                ));
            }
            if seeds.len() != 2 {
                return Ok(Analysis::not_applicable(
                    "an order-two rule with an exceptional prefix".into(),
                ));
            }
            order_two(map, &current, &lag, options)
        }
        _ => Ok(Analysis::not_applicable(format!(
            "a rule of order {order} (above two)"
        ))),
    }
}

mod published;
use published::*;
mod differential;
use differential::*;
