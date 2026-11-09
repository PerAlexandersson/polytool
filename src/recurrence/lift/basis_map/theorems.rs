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
    /// A hypothesis fails; the theorem gives no conclusion (nothing is
    /// claimed against the map).
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
        let recomputed = check_basis_map_preserves_pf(&map, &BasisMapTheoremOptions::default())?;
        let comparable = |r: &Self| {
            (
                r.shape.clone(),
                r.status,
                r.hypotheses.clone(),
                r.conclusions.clone(),
            )
        };
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
        return scalar_index(map, options);
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

// ---------------------------------------------------------------------------
// Published transformations
// ---------------------------------------------------------------------------

/// The magic-cone lemma used by the Eulerian and derangement
/// transformations (module documentation of `basis_map`, "Published
/// transformations").
const MAGIC_CONE_LEMMA: &str = "an input p of degree d with every root in [-1, 0] and a positive \
     leading coefficient is a nonnegative combination of u^k (1 + u)^(d - k), 0 <= k <= d: with \
     p = c prod_i (u - x_i), sum_k c_k s^k = (1 - s)^d p(s / (1 - s)) = c prod_i ((1 + x_i) s - \
     x_i) has nonnegative coefficients";

/// The analysis of a map recognised as a published family, `None` when the
/// first rows match none.
fn published(map: &BasisMap) -> Result<Option<Analysis>, RecurrenceError> {
    let Some((family, certificate)) = recognize_published_family(map)? else {
        return Ok(None);
    };
    let (shape, theorem, rows, lean): (BasisMapShape, &str, String, Vec<&'static str>) =
        match family {
            PublishedFamily::Narayana(m) => (
                BasisMapShape::Narayana { parameter: m },
                "Mao and Wang, The Narayana transformation, arXiv:2607.01572, Theorem 1: for \
                 every integer m >= 0 the map u^k -> N_(k,m)(t) = sum_j C(k,j) C(k+m,j) / \
                 C(m+j,j) t^j sends every polynomial with nonnegative coefficients and only real \
                 roots to a polynomial with only real nonpositive roots (PF to PF)",
                format!("N_(k,{m})(t)"),
                vec![
                    "RealRooted.narayanaTransformPreservesPF",
                    "RealRooted.narayanaPolynomial_deriv_lag_rec",
                    "RealRooted.narayanaPolynomial_one",
                ],
            ),
            PublishedFamily::EulerianTransformation => (
                BasisMapShape::EulerianTransformation,
                "Athanasiadis, On the real-rootedness of the Eulerian transformation, \
                 arXiv:2302.00754 (J. London Math. Soc. 111 (2025) e70083), Theorem 1: A° \
                 sends every nonnegative combination of u^k (1 + u)^(d - k), 0 <= k <= d, to a \
                 real-rooted polynomial with nonnegative coefficients (Brenti's stronger \
                 conjecture for all inputs with nonpositive roots is false, Brändén and \
                 Jochemko; it is not used)",
                "1, t A_1(t), t A_2(t), ... (t A_k(t) = sum over permutations of [k] of \
                 t^(1 + des))"
                    .into(),
                vec![
                    "RealRooted.EulerianTransform.transform_magicExpansion_isPF",
                    "RealRooted.EulerianTransform.transform_magicExpansion_interl",
                    "RealRooted.generalizedEulerian",
                ],
            ),
            PublishedFamily::DerangementTransformation => (
                BasisMapShape::DerangementTransformation,
                "Brändén and Solus, arXiv:1808.04141, Corollary 3.7 (as stated by Athanasiadis, \
                 arXiv:2302.00754, Section 4): the derangement transformation u^k -> d_k(t) \
                 sends every nonnegative combination of u^k (1 + u)^(d - k), 0 <= k <= d, to \
                 zero or a real-rooted polynomial with nonnegative coefficients",
                "the derangement polynomials d_k(t) by excedances (d_0 = 1, d_1 = 0)".into(),
                vec![
                    "RealRooted.DerangementTransform.transform_magicExpansion_isPF",
                    "RealRooted.DerangementTransform.polynomial_succ",
                ],
            ),
        };
    let mut analysis = Analysis {
        shape,
        status: BasisMapTheoremStatus::Proved,
        theorem: theorem.into(),
        hypotheses: Vec::new(),
        conclusions: Vec::new(),
        counterexample: None,
        lean,
    };
    let identified = certificate.is_certified();
    analysis.hypotheses.push(hypothesis(
        "basis_is_published_family",
        format!(
            "the basis row of u^k is {rows} for every k >= 0 (P_({} + k))",
            map.first_index()
        ),
        identified,
        if identified {
            format!(
                "all-index recurrence equality with the RealRooted recurrence ({}): {}",
                certificate.method.as_deref().unwrap_or("certified"),
                certificate.statement
            )
        } else {
            format!(
                "the first rows agree but the all-index equality is not certified: {}",
                certificate.reasons.join("; ")
            )
        },
    ));
    if !identified {
        return Ok(Some(analysis));
    }
    let nonpositive = AlgebraicInterval::nonpositive_ray();
    match family {
        PublishedFamily::Narayana(_) => analysis.push_conclusion(BasisMapConclusion {
            input_interval: nonpositive.clone(),
            output_interval: nonpositive,
            nonnegative_coefficients: true,
            sharp: false,
            statement: "every input with all roots in (-inf, 0] (a PF polynomial up to sign) is \
                        sent to zero or a polynomial with all roots in (-inf, 0], PF for a \
                        positive leading coefficient"
                .into(),
            hypotheses: Vec::new(),
        }),
        PublishedFamily::EulerianTransformation | PublishedFamily::DerangementTransformation => {
            analysis.push_conclusion(BasisMapConclusion {
                input_interval: AlgebraicInterval::closed(
                    QuadraticNumber::rational(-Q::one()),
                    QuadraticNumber::rational(Q::zero()),
                ),
                output_interval: nonpositive,
                nonnegative_coefficients: true,
                sharp: false,
                statement: format!(
                    "every input with all roots in [-1, 0] is sent to zero or a polynomial with \
                     all roots in (-inf, 0], PF for a positive leading coefficient ({MAGIC_CONE_LEMMA})"
                ),
                hypotheses: Vec::new(),
            })
        }
    }
    Ok(Some(analysis))
}

// ---------------------------------------------------------------------------
// Geometric families
// ---------------------------------------------------------------------------

fn geometric(g: &[Q], y: &[Q]) -> Result<Analysis, RecurrenceError> {
    let a = vec_poly::coeff(y, 0);
    let d = vec_poly::coeff(y, 1);
    let mut analysis = Analysis {
        shape: BasisMapShape::Geometric {
            initial_row: poly_text(g),
            constant: text(&a),
            slope: text(&d),
        },
        status: BasisMapTheoremStatus::Proved,
        theorem: "geometric family P_k = g(t) (a + d t)^k: T(p) = g(t) p(a + d t), so the roots \
                  of T(p) are those of g and the points (r - a)/d for the roots r of p \
                  (elementary)"
            .into(),
        hypotheses: Vec::new(),
        conclusions: Vec::new(),
        counterexample: None,
        lean: Vec::new(),
    };
    let nonzero = !vec_poly::is_zero(g);
    analysis.hypotheses.push(hypothesis(
        "initial_row_nonzero",
        "g = P_0 is nonzero".into(),
        nonzero,
        "exact".into(),
    ));
    if !nonzero {
        return Ok(analysis);
    }
    let real = is_real_rooted(g);
    analysis.hypotheses.push(hypothesis(
        "initial_row_real_rooted",
        "g is real-rooted".into(),
        real,
        "exact root isolation".into(),
    ));
    if !real {
        analysis.status = BasisMapTheoremStatus::NotPreserving;
        analysis.counterexample = Some(counterexample(
            &[Q::one()],
            g,
            "T(1) = g is not real-rooted".into(),
        ));
        return Ok(analysis);
    }
    if !d.is_zero() {
        analysis.push_conclusion(BasisMapConclusion {
            input_interval: AlgebraicInterval::whole_line(),
            output_interval: AlgebraicInterval::whole_line(),
            nonnegative_coefficients: false,
            sharp: true,
            statement: "every real-rooted input has a real-rooted image (an affine substitution \
                        times g)"
                .into(),
            hypotheses: Vec::new(),
        });
    }
    let nonpositive = AlgebraicInterval::nonpositive_ray().contains_all_roots_of(g);
    analysis.hypotheses.push(hypothesis(
        "initial_row_nonpositive_roots",
        "every root of g is <= 0".into(),
        nonpositive,
        "exact root isolation".into(),
    ));
    if !nonpositive {
        return Ok(analysis);
    }
    let positive_leading = leading(g).is_positive();
    let (input, nonnegative, sharp, statement) = match d.cmp(&Q::zero()) {
        Ordering::Greater => (
            AlgebraicInterval::at_most(rational(&a)),
            positive_leading,
            true,
            format!(
                "an input root r <= {} gives the root (r - a)/d <= 0; a larger root gives a \
                 positive root",
                text(&a)
            ),
        ),
        Ordering::Less => (
            AlgebraicInterval::at_least(rational(&a)),
            false,
            true,
            format!(
                "an input root r >= {} gives the root (r - a)/d <= 0",
                text(&a)
            ),
        ),
        Ordering::Equal => (
            AlgebraicInterval::at_most(rational(&a)),
            positive_leading,
            false,
            "T(p) = p(a) g; an input with every root <= a has p(a) of the sign of its leading \
             coefficient"
                .into(),
        ),
    };
    analysis.push_conclusion(BasisMapConclusion {
        input_interval: input,
        output_interval: AlgebraicInterval::nonpositive_ray(),
        nonnegative_coefficients: nonnegative,
        sharp,
        statement,
        hypotheses: Vec::new(),
    });
    Ok(analysis)
}

/// An eventually geometric family with a nonzero prefix discrepancy:
/// no positive-length source interval (RealRooted
/// `exists_not_splits_linearMap_shiftedPower_of_eventually_geometric`).
fn eventually_geometric(map: &BasisMap, y: &[Q]) -> Result<Analysis, RecurrenceError> {
    let a = vec_poly::coeff(y, 0);
    let d = vec_poly::coeff(y, 1);
    let seeds = map.seeds();
    let k = seeds.len() - 1;
    let mut analysis = Analysis {
        shape: BasisMapShape::EventuallyGeometric {
            geometric_from: k,
            constant: text(&a),
            slope: text(&d),
        },
        status: BasisMapTheoremStatus::HypothesisFailed,
        theorem: "eventually geometric families P_n = P_K y^(n - K) (n >= K >= 1, y = a + d t, \
                  d != 0, P_K != 0) with a nonzero discrepancy E_j = y^K P_j - P_K y^j, j < K: \
                  for every r != 0 some T((u + r)^N), max(K, D + 3) <= N <= max(K, D + 3) + K - \
                  1 (D the largest degree of a nonzero E_j), is not real-rooted, so no source \
                  interval of positive length works"
            .into(),
        hypotheses: Vec::new(),
        conclusions: Vec::new(),
        counterexample: None,
        lean: vec!["RealRooted.exists_not_splits_linearMap_shiftedPower_of_eventually_geometric"],
    };
    let pk = &seeds[k];
    let slope_nonzero = !d.is_zero();
    analysis.hypotheses.push(hypothesis(
        "nonconstant_multiplier",
        "d != 0".into(),
        slope_nonzero,
        "exact".into(),
    ));
    let pk_nonzero = !vec_poly::is_zero(pk);
    analysis.hypotheses.push(hypothesis(
        "geometric_row_nonzero",
        format!("P_{k} != 0"),
        pk_nonzero,
        "exact".into(),
    ));
    if !slope_nonzero || !pk_nonzero {
        return Ok(analysis);
    }
    let yk = vec_poly::pow(y, k);
    let mut degree = None;
    for (j, pj) in seeds.iter().enumerate().take(k) {
        let e = vec_poly::trim(&vec_poly::sub(
            &vec_poly::mul(&yk, pj),
            &vec_poly::mul(pk, &vec_poly::pow(y, j)),
        ));
        if !vec_poly::is_zero(&e) {
            let deg = vec_poly::degree(&e).unwrap_or(0);
            degree = Some(degree.map_or(deg, |m: usize| m.max(deg)));
        }
    }
    let Some(degree) = degree else {
        analysis.hypotheses.push(hypothesis(
            "nonzero_discrepancy",
            "some discrepancy E_j is nonzero".into(),
            false,
            "every discrepancy vanishes".into(),
        ));
        return Ok(analysis);
    };
    analysis.hypotheses.push(hypothesis(
        "nonzero_discrepancy",
        "some discrepancy E_j is nonzero".into(),
        true,
        format!("exact; largest degree D = {degree}"),
    ));
    analysis.status = BasisMapTheoremStatus::NotPreserving;
    let start = k.max(degree + 3);
    for n in start..start + k {
        let input = vec_poly::pow(&[Q::one(), Q::one()], n);
        let image = map.apply(&input)?;
        if !is_real_rooted(&image) {
            analysis.counterexample = Some(counterexample(
                &input,
                &image,
                format!("T((u + 1)^{n}) is not real-rooted (exact)"),
            ));
            break;
        }
    }
    Ok(analysis)
}

// ---------------------------------------------------------------------------
// First-order differential families
// ---------------------------------------------------------------------------

fn differential(
    map: &BasisMap,
    g: &[Q],
    multiplier: &[Q],
    v: &[Q],
    _options: &BasisMapTheoremOptions,
) -> Result<Analysis, RecurrenceError> {
    let _ = map;
    let a = vec_poly::coeff(multiplier, 0);
    let b = vec_poly::coeff(multiplier, 1);
    let v = vec_poly::trim(v);
    let degree = vec_poly::degree(&v).unwrap_or(0);
    let g_text = poly_text(g);
    let positive_leading = leading(g).is_positive();
    match degree {
        0 => {
            let kappa = -vec_poly::coeff(&v, 0);
            let mut analysis = Analysis {
                shape: BasisMapShape::HermiteDifferential {
                    initial_row: g_text,
                    a: text(&a),
                    b: text(&b),
                    kappa: text(&kappa),
                },
                status: BasisMapTheoremStatus::Proved,
                theorem: "Hermite-type operators L f = (a + b t) f - kappa f', b > 0, kappa >= \
                          0: every shift L - rho keeps real-rootedness (the derivative \
                          coefficient is <= 0 everywhere: Ma and Wang, EJC 15 (2008) R17, \
                          Theorem 2), and T(p) = c prod_j (L - rho_j) g for p = c prod_j (u - \
                          rho_j); so every real-rooted input has a real-rooted image"
                    .into(),
                hypotheses: Vec::new(),
                conclusions: Vec::new(),
                counterexample: None,
                lean: Vec::new(),
            };
            let checks = [
                (
                    "multiplier_slope_positive",
                    "b > 0".to_string(),
                    b.is_positive(),
                ),
                (
                    "kappa_positive",
                    "kappa > 0".to_string(),
                    kappa.is_positive(),
                ),
            ];
            for (name, statement, holds) in checks {
                analysis
                    .hypotheses
                    .push(hypothesis(name, statement, holds, "exact".into()));
            }
            let real = !vec_poly::is_zero(g) && is_real_rooted(g);
            analysis.hypotheses.push(hypothesis(
                "initial_row_real_rooted",
                "g = P_0 is nonzero and real-rooted".into(),
                real,
                "exact root isolation".into(),
            ));
            if !vec_poly::is_zero(g) && !real {
                analysis.status = BasisMapTheoremStatus::NotPreserving;
                analysis.counterexample = Some(counterexample(
                    &[Q::one()],
                    g,
                    "T(1) = g is not real-rooted".into(),
                ));
                return Ok(analysis);
            }
            if analysis
                .hypotheses
                .iter()
                .all(|h| h.status == LiftHypothesisStatus::Verified)
            {
                analysis.push_conclusion(BasisMapConclusion {
                    input_interval: AlgebraicInterval::whole_line(),
                    output_interval: AlgebraicInterval::whole_line(),
                    nonnegative_coefficients: false,
                    sharp: true,
                    statement: "every real-rooted input has a real-rooted image".into(),
                    hypotheses: Vec::new(),
                });
            }
            Ok(analysis)
        }
        1 => {
            let kappa = vec_poly::coeff(&v, 1);
            let h = -vec_poly::coeff(&v, 0) / &kappa;
            let mut analysis = Analysis {
                shape: BasisMapShape::RayDifferential {
                    initial_row: g_text,
                    a: text(&a),
                    b: text(&b),
                    kappa: text(&kappa),
                    root: text(&h),
                },
                status: BasisMapTheoremStatus::Proved,
                theorem: "ray basis maps: for L f = (a + b t) f + kappa (t - h) f', b > 0, kappa \
                          >= 0, and g real-rooted with every root <= h, of multiplicity m at h, \
                          T(u^k) = L^k g sends every real polynomial with all roots in S to one \
                          with all roots in (-inf, h] if and only if S lies in (-inf, a + b h + \
                          kappa m] (RealRooted rayBasisMap_maps_realRootedIn_iff; the step is \
                          Ma and Wang, EJC 15 (2008) R17, Theorem 2)"
                    .into(),
                hypotheses: Vec::new(),
                conclusions: Vec::new(),
                counterexample: None,
                lean: vec!["RealRooted.rayBasisMap_maps_realRootedIn_iff"],
            };
            analysis.hypotheses.push(hypothesis(
                "multiplier_slope_positive",
                "b > 0".into(),
                b.is_positive(),
                "exact".into(),
            ));
            analysis.hypotheses.push(hypothesis(
                "kappa_nonnegative",
                "kappa >= 0".into(),
                !kappa.is_negative(),
                "exact".into(),
            ));
            let nonzero = !vec_poly::is_zero(g);
            let rooted =
                nonzero && AlgebraicInterval::at_most(rational(&h)).contains_all_roots_of(g);
            analysis.hypotheses.push(hypothesis(
                "initial_row_roots",
                format!(
                    "g = P_0 is nonzero and real-rooted with every root <= {}",
                    text(&h)
                ),
                rooted,
                "exact root isolation".into(),
            ));
            if nonzero && !is_real_rooted(g) {
                analysis.status = BasisMapTheoremStatus::NotPreserving;
                analysis.counterexample = Some(counterexample(
                    &[Q::one()],
                    g,
                    "T(1) = g is not real-rooted".into(),
                ));
                return Ok(analysis);
            }
            if analysis
                .hypotheses
                .iter()
                .all(|x| x.status == LiftHypothesisStatus::Verified)
            {
                let m = multiplicity(g, &h);
                let upper = &a + &b * &h + &kappa * integer(m as i64);
                let output = AlgebraicInterval::at_most(rational(&h));
                analysis.push_conclusion(BasisMapConclusion {
                    input_interval: AlgebraicInterval::at_most(rational(&upper)),
                    nonnegative_coefficients: positive_leading && !h.is_positive(),
                    output_interval: output,
                    sharp: true,
                    statement: format!(
                        "source (-inf, a + b h + kappa m] = (-inf, {}] with m = {m}, target \
                         (-inf, {}]",
                        text(&upper),
                        text(&h)
                    ),
                    hypotheses: Vec::new(),
                });
            }
            Ok(analysis)
        }
        2 => {
            let kappa = vec_poly::coeff(&v, 2);
            let p = vec_poly::coeff(&v, 1) / &kappa;
            let q0 = vec_poly::coeff(&v, 0) / &kappa;
            // Roots of t^2 + p t + q0.
            let discriminant = &p * &p - integer(4) * &q0;
            let root = QuadraticNumber::square_root(discriminant.clone()).ok();
            let rational_roots = discriminant.is_positive()
                && root.as_ref().is_some_and(QuadraticNumber::is_rational);
            if !rational_roots {
                return Ok(Analysis::not_applicable(
                    "the derivative coefficient is quadratic without two distinct rational \
                     roots"
                        .into(),
                ));
            }
            let s = root
                .expect("checked")
                .as_rational()
                .cloned()
                .expect("rational");
            let half = Q::new(BigInt::one(), BigInt::from(2));
            let lower = (-&p - &s) * &half;
            let upper = (-&p + &s) * &half;
            let mut analysis = Analysis {
                shape: BasisMapShape::IntervalDifferential {
                    initial_row: g_text,
                    a: text(&a),
                    b: text(&b),
                    kappa: text(&kappa),
                    lower_root: text(&lower),
                    upper_root: text(&upper),
                },
                status: BasisMapTheoremStatus::Proved,
                theorem: "interval basis maps: for L f = (a + b t) f + kappa (t - l)(t - h) f', \
                          l < h, kappa > 0, b + kappa deg g > 0, and g with positive leading \
                          coefficient, real-rooted in [l, h] with endpoint multiplicities m_l, \
                          m_h, T(u^k) = L^k g sends every real polynomial with all roots in S \
                          to one with all roots in [l, h] if and only if S lies in [a + b l - \
                          kappa (h - l) m_l, a + b h + kappa (h - l) m_h] (RealRooted \
                          intervalBasisMap_maps_realRootedIn_iff; Ma and Wang, EJC 15 (2008) \
                          R17, Theorem 2 and Corollary 3); a negative leading coefficient of g \
                          only changes the sign of T"
                    .into(),
                hypotheses: Vec::new(),
                conclusions: Vec::new(),
                counterexample: None,
                lean: vec!["RealRooted.intervalBasisMap_maps_realRootedIn_iff"],
            };
            let deg_g = vec_poly::degree(g).unwrap_or(0);
            let top = &b + &kappa * integer(deg_g as i64);
            analysis.hypotheses.push(hypothesis(
                "kappa_positive",
                "kappa > 0".into(),
                kappa.is_positive(),
                "exact".into(),
            ));
            analysis.hypotheses.push(hypothesis(
                "top_coefficient_positive",
                format!("b + kappa deg g = {} > 0", text(&top)),
                top.is_positive(),
                "exact".into(),
            ));
            let interval = AlgebraicInterval::closed(rational(&lower), rational(&upper));
            let nonzero = !vec_poly::is_zero(g);
            let rooted = nonzero && interval.contains_all_roots_of(g);
            analysis.hypotheses.push(hypothesis(
                "initial_row_roots",
                format!(
                    "g = P_0 is nonzero and real-rooted with every root in [{}, {}]",
                    text(&lower),
                    text(&upper)
                ),
                rooted,
                "exact root isolation".into(),
            ));
            if nonzero && !is_real_rooted(g) {
                analysis.status = BasisMapTheoremStatus::NotPreserving;
                analysis.counterexample = Some(counterexample(
                    &[Q::one()],
                    g,
                    "T(1) = g is not real-rooted".into(),
                ));
                return Ok(analysis);
            }
            if analysis
                .hypotheses
                .iter()
                .all(|x| x.status == LiftHypothesisStatus::Verified)
            {
                let m_lower = integer(multiplicity(g, &lower) as i64);
                let m_upper = integer(multiplicity(g, &upper) as i64);
                let width = &upper - &lower;
                let source_lower = &a + &b * &lower - &kappa * &width * &m_lower;
                let source_upper = &a + &b * &upper + &kappa * &width * &m_upper;
                analysis.push_conclusion(BasisMapConclusion {
                    input_interval: AlgebraicInterval::closed(
                        rational(&source_lower),
                        rational(&source_upper),
                    ),
                    output_interval: interval,
                    nonnegative_coefficients: positive_leading && !upper.is_positive(),
                    sharp: true,
                    statement: format!(
                        "source [{}, {}], target [{}, {}]",
                        text(&source_lower),
                        text(&source_upper),
                        text(&lower),
                        text(&upper)
                    ),
                    hypotheses: Vec::new(),
                });
            }
            Ok(analysis)
        }
        _ => Ok(Analysis::not_applicable(
            "the derivative coefficient has degree above two".into(),
        )),
    }
}

// ---------------------------------------------------------------------------
// Index-dependent scalar term
// ---------------------------------------------------------------------------

fn scalar_index(
    map: &BasisMap,
    options: &BasisMapTheoremOptions,
) -> Result<Analysis, RecurrenceError> {
    let rule = map.rule();
    let not_applicable = |reason: &str| Ok(Analysis::not_applicable(reason.to_string()));
    if rule.max_offset() != 1 || map.seeds().len() != 1 {
        return not_applicable("an index-dependent rule other than first order without a prefix");
    }
    let leader = match &rule.denominator {
        None => Q::one(),
        Some(f) => match Biv::from_bivar(f).as_constant() {
            Some(c) => c,
            None => return not_applicable("an index-dependent leader"),
        },
    };
    let mut multiplier = Biv::zero();
    let mut v = Biv::zero();
    for term in &rule.terms {
        let c = Biv::from_bivar(&term.coeff).scale(&leader.recip());
        match term.deriv_order {
            0 => multiplier = multiplier.add(&c),
            1 => v = v.add(&c),
            _ => return not_applicable("an index-dependent rule with higher derivatives"),
        }
    }
    if !v.is_index_free() {
        return not_applicable("an index-dependent derivative coefficient");
    }
    // Only the t^0 coefficient of the multiplier may depend on n.
    let rows = multiplier.rows();
    if rows.iter().skip(1).any(|row| row.len() > 1) {
        return not_applicable("an index-dependent coefficient of t");
    }
    let alpha_n: Vec<Q> = vec_poly::trim(
        &rows
            .iter()
            .map(|row| row.first().cloned().unwrap_or_else(Q::zero))
            .collect::<Vec<_>>(),
    );
    let t_part = t_polynomial(&multiplier);
    let b = vec_poly::coeff(&t_part, 1);
    if vec_poly::degree(&t_part).unwrap_or(0) > 1 {
        return not_applicable("an index-dependent rule whose multiplier has degree above one");
    }
    let v = t_polynomial(&v);
    if vec_poly::degree(&v).unwrap_or(0) > 1 {
        return not_applicable("an index-dependent rule with a nonlinear derivative coefficient");
    }
    // alpha as a polynomial in the basis index k: n = i0 + 1 + k.
    let shift = integer(map.first_index() as i64 + 1);
    let alpha_k = vec_poly::compose(&alpha_n, &[shift, Q::one()]);
    let kappa = vec_poly::coeff(&v, 1);
    let a = vec_poly::coeff(&alpha_k, 0);
    // Without a derivative term the ray root h is free: the smallest h >= 0
    // with a + b h >= 0.
    let h = if !kappa.is_zero() {
        -vec_poly::coeff(&v, 0) / &kappa
    } else if a.is_negative() && b.is_positive() {
        -(&a / &b)
    } else {
        Q::zero()
    };
    let g = &map.seeds()[0];
    let slope = vec_poly::coeff(&alpha_k, 1);
    let affine = vec_poly::degree(&alpha_k).unwrap_or(0) <= 1;
    let mut analysis = Analysis {
        shape: BasisMapShape::ScalarIndex {
            initial_row: poly_text(g),
            scalar: poly_text(&alpha_k),
            b: text(&b),
            kappa: text(&kappa),
            root: text(&h),
        },
        status: BasisMapTheoremStatus::Proved,
        theorem: "index-dependent scalar term: P_(k+1) = (L + k s) P_k with L f = (a + b t) f + \
                  kappa (t - h) f' gives T = U_g o R_s (RealRooted \
                  basisMap_factorization_of_scalarIndex_recurrence), where R_s(u^k) = u (u + s) \
                  ... (u + (k-1) s) keeps every root <= 0 for s >= 0 (Su, Yang and Zhang, \
                  Australas. J. Combin. 56 (2013), Theorem 4; Brenti, Mem. AMS 413 (1989), \
                  Theorem 2.4.3; RealRooted NewtonParametersAffine.preservesNonposRooted) and \
                  U_g(u^k) = L^k g sends (-inf, a + b h + kappa m] into (-inf, h] (RealRooted \
                  rayBasisMap_maps_realRootedIn_iff)"
            .into(),
        hypotheses: Vec::new(),
        conclusions: Vec::new(),
        counterexample: None,
        lean: vec![
            "RealRooted.basisMap_factorization_of_scalarIndex_recurrence",
            "RealRooted.NewtonParametersAffine.preservesNonposRooted",
            "RealRooted.rayBasisMap_maps_realRootedIn_iff",
        ],
    };
    let pure_newton = kappa.is_zero() && is_constant(g) && !vec_poly::is_zero(g) && b.is_positive();
    let slope_ok = affine && !slope.is_negative();
    analysis.hypotheses.push(hypothesis(
        "scalar_affine_nondecreasing",
        format!(
            "alpha(k) = {} is affine in k with slope >= 0",
            crate::recurrence::implicit::algebra::poly_text(&alpha_k, "k")
        ),
        slope_ok,
        "exact".into(),
    ));
    if !slope_ok && pure_newton {
        // Newton basis prod_(j < k) (b t + alpha(j)): iff theorem.
        analysis.status = BasisMapTheoremStatus::NotPreserving;
        analysis.theorem = "Newton bases P_k = prod_(j < k) (b t + c_j), b > 0: T preserves \
                            real-rootedness of the inputs with every root <= 0 if and only if \
                            c_j is an arithmetic progression with nonnegative step (RealRooted \
                            newtonPreservesNonposRooted_iff_parametersAffine; Su, Yang and Zhang, \
                            Theorem 4)"
            .into();
        analysis.lean = vec!["RealRooted.newtonPreservesNonposRooted_iff_parametersAffine"];
        let limit = options.witness_search_limit.clamp(1, 12);
        'search: for r in [1i64, 2, 4, 8, 16] {
            for n in 2..=limit {
                let input = vec_poly::pow(&[integer(r), Q::one()], n);
                let image = map.apply(&input)?;
                if !is_real_rooted(&image) {
                    analysis.counterexample = Some(counterexample(
                        &input,
                        &image,
                        format!("T((u + {r})^{n}) is not real-rooted (exact)"),
                    ));
                    break 'search;
                }
            }
        }
        return Ok(analysis);
    }
    analysis.hypotheses.push(hypothesis(
        "multiplier_slope_positive",
        "b > 0".into(),
        b.is_positive(),
        "exact".into(),
    ));
    analysis.hypotheses.push(hypothesis(
        "kappa_nonnegative",
        "kappa >= 0".into(),
        !kappa.is_negative(),
        "exact".into(),
    ));
    let nonzero = !vec_poly::is_zero(g);
    let rooted = nonzero && AlgebraicInterval::at_most(rational(&h)).contains_all_roots_of(g);
    analysis.hypotheses.push(hypothesis(
        "initial_row_roots",
        format!(
            "g = P_0 is nonzero and real-rooted with every root <= {}",
            text(&h)
        ),
        rooted,
        "exact root isolation".into(),
    ));
    let m = if rooted { multiplicity(g, &h) } else { 0 };
    let source_end = &a + &b * &h + &kappa * integer(m as i64);
    analysis.hypotheses.push(hypothesis(
        "ray_source_contains_nonpositive_ray",
        format!(
            "a + b h + kappa m = {} >= 0 (the source of U_g contains (-inf, 0])",
            text(&source_end)
        ),
        !source_end.is_negative(),
        "exact".into(),
    ));
    if analysis
        .hypotheses
        .iter()
        .all(|x| x.status == LiftHypothesisStatus::Verified)
    {
        analysis.push_conclusion(BasisMapConclusion {
            input_interval: AlgebraicInterval::nonpositive_ray(),
            output_interval: AlgebraicInterval::at_most(rational(&h)),
            nonnegative_coefficients: leading(g).is_positive() && !h.is_positive(),
            sharp: false,
            statement: format!(
                "(-inf, 0] -> (-inf, {}] through R_s, s = {}, and U_g",
                text(&h),
                text(&slope)
            ),
            hypotheses: Vec::new(),
        });
    }
    Ok(analysis)
}

// ---------------------------------------------------------------------------
// Order two
// ---------------------------------------------------------------------------

fn order_two(
    map: &BasisMap,
    current: &[Q],
    lag: &[Q],
    options: &BasisMapTheoremOptions,
) -> Result<Analysis, RecurrenceError> {
    let seeds = map.seeds();
    let p0 = &seeds[0];
    let p1 = &seeds[1];
    if !is_constant(p0) || vec_poly::is_zero(p0) || vec_poly::degree(p1).unwrap_or(0) > 1 {
        return Ok(Analysis::not_applicable(
            "an order-two rule whose first row is not a nonzero constant or whose second row \
             has degree above one"
                .into(),
        ));
    }
    let normalization = p0[0].clone();
    let l = vec_poly::coeff(p1, 0) / &normalization;
    let m = vec_poly::coeff(p1, 1) / &normalization;
    let a = vec_poly::coeff(current, 0);
    let d = vec_poly::coeff(current, 1);
    let c = vec_poly::coeff(lag, 0);
    let b = vec_poly::coeff(lag, 1);
    if d.is_zero() {
        if !m.is_zero() {
            return Ok(Analysis::not_applicable(
                "a constant current multiplier with a nonconstant second row".into(),
            ));
        }
        return Ok(constant_current(&normalization, &a, &b, &c, &l, options));
    }
    if d.is_negative() {
        return Ok(Analysis::not_applicable(
            "an affine current multiplier with negative slope".into(),
        ));
    }
    Ok(affine_current(&normalization, &a, &d, &b, &c, &l, &m))
}

/// The scalar criterion `v_n >= 0` for `v_0 = 1`, `v_1 = s`,
/// `v_(n+2) = a v_(n+1) + c v_n` with `a, s >= 0` (RealRooted
/// `orderTwoConstants_nonneg_iff`): `a^2 + 4c >= 0` and
/// (`2s >= a` or `c >= s^2 - a s`).
fn scalar_criterion(a: &Q, c: &Q, s: &Q) -> (bool, String) {
    let disc = a * a + integer(4) * c;
    let first = !disc.is_negative();
    let second = integer(2) * s >= *a || *c >= s * s - a * s;
    (
        first && second,
        format!(
            "a^2 + 4c = {} >= 0: {first}; 2s >= a or c >= s^2 - a s: {second}",
            text(&disc)
        ),
    )
}

fn constant_current(
    normalization: &Q,
    a: &Q,
    b: &Q,
    c: &Q,
    s: &Q,
    options: &BasisMapTheoremOptions,
) -> Analysis {
    let mut analysis = Analysis {
        shape: BasisMapShape::ConstantCurrent {
            normalization: text(normalization),
            a: text(a),
            b: text(b),
            c: text(c),
            initial_constant: text(s),
        },
        status: BasisMapTheoremStatus::Proved,
        theorem: "constant current multiplier: P_0 = 1, P_1 = s, P_(k+2) = a P_(k+1) + (b t + c) \
                  P_k with a, s >= 0, b > 0.  P_k(t) = Q_k(b t + c) for the rows Q_k with b = 1, \
                  c = 0, whose map sends PF inputs to PF outputs (an interlacing-preserving 2 x \
                  2 matrix, Branden, arXiv:1410.6601, Theorem 8.5), so every PF input has a \
                  real-rooted image with roots <= -c/b; PF inputs (also those with roots in \
                  any [-R, 0]) have PF images if and only if every P_k(0) >= 0, if and only if \
                  a^2 + 4c >= 0 and (2s >= a or c >= s^2 - a s) (RealRooted \
                  orderTwoBasis_pf_equivalences, basisTransform_orderTwoBasis_preserves_bounded_\
                  pf_iff)"
            .into(),
        hypotheses: Vec::new(),
        conclusions: Vec::new(),
        counterexample: None,
        lean: vec![
            "RealRooted.basisTransform_orderTwoBasis_eq_comp",
            "RealRooted.orderTwoZeroTransform_isPFPolynomial",
            "RealRooted.orderTwoBasis_pf_equivalences",
            "RealRooted.basisTransform_orderTwoBasis_preserves_bounded_pf_iff",
        ],
    };
    for (name, statement, holds) in [
        ("current_nonnegative", "a >= 0", !a.is_negative()),
        ("initial_constant_nonnegative", "s >= 0", !s.is_negative()),
        ("lag_slope_positive", "b > 0", b.is_positive()),
    ] {
        analysis
            .hypotheses
            .push(hypothesis(name, statement.into(), holds, "exact".into()));
    }
    if analysis
        .hypotheses
        .iter()
        .any(|h| h.status != LiftHypothesisStatus::Verified)
    {
        return analysis;
    }
    let root_bound = -(c / b);
    analysis.push_conclusion(BasisMapConclusion {
        input_interval: AlgebraicInterval::nonpositive_ray(),
        output_interval: AlgebraicInterval::at_most(rational(&root_bound)),
        nonnegative_coefficients: false,
        sharp: false,
        statement: format!(
            "every PF input has a real-rooted image with every root <= -c/b = {}",
            text(&root_bound)
        ),
        hypotheses: Vec::new(),
    });
    let (criterion, evidence) = scalar_criterion(a, c, s);
    analysis.hypotheses.push(hypothesis(
        "constant_terms_nonnegative",
        "a^2 + 4c >= 0 and (2s >= a or c >= s^2 - a s), equivalently every P_k(0) >= 0".into(),
        criterion,
        evidence,
    ));
    if criterion {
        let bound = if root_bound.is_negative() {
            root_bound
        } else {
            Q::zero()
        };
        analysis.push_conclusion(BasisMapConclusion {
            input_interval: AlgebraicInterval::nonpositive_ray(),
            output_interval: AlgebraicInterval::at_most(rational(&bound)),
            nonnegative_coefficients: normalization.is_positive(),
            sharp: true,
            statement: "every PF input has a PF image (zero allowed)".into(),
            hypotheses: Vec::new(),
        });
    } else {
        // The first monomial with a negative constant term.
        let (mut previous, mut current) = (Q::one(), s.clone());
        for n in 2..=options.witness_search_limit {
            let next = a * &current + c * &previous;
            previous = current;
            current = next;
            if current.is_negative() {
                analysis.counterexample = Some(BasisMapCounterexample {
                    input: {
                        let mut monomial = vec!["0".to_string(); n];
                        monomial.push("1".into());
                        monomial
                    },
                    image: Vec::new(),
                    failure: format!(
                        "P_{n}(0) / P_0 = {} < 0 while P_{n} has a positive leading \
                         coefficient: the PF input u^{n} has an image with a positive root (the \
                         image is omitted)",
                        text(&current)
                    ),
                });
                break;
            }
        }
    }
    analysis
}

fn affine_current(normalization: &Q, a: &Q, d: &Q, b: &Q, c: &Q, l: &Q, m: &Q) -> Analysis {
    let beta = b / d;
    let k = &beta * &beta + a * &beta - c;
    let eta = m / d;
    let s = l - &eta * a;
    let mut analysis = Analysis {
        shape: BasisMapShape::AffineCurrent {
            normalization: text(normalization),
            a: text(a),
            d: text(d),
            b: text(b),
            c: text(c),
            initial_constant: text(l),
            initial_slope: text(m),
        },
        status: BasisMapTheoremStatus::Proved,
        theorem: "affine current multiplier with an affine initial row (module documentation, \
                  \"The affine-current theorem\"): with beta = b/d, K = beta^2 + a beta - c > 0 \
                  and eta = m/d >= 0, every real polynomial with all roots in E = [-beta - \
                  sqrt(K), -beta + sqrt(K)] has a real-rooted image, and no interval of positive \
                  length outside E works (Theorem 1); for eta > 0 the images lie in I = [(-a - \
                  2 beta - 2 sqrt(K))/d, (-a - 2 beta + 2 sqrt(K))/d] if and only if eta >= 1/2 \
                  and (s - (2 eta - 1) beta)^2 <= (2 eta - 1)^2 K, s = l - eta a (Theorem 2); \
                  with a >= 0, K >= beta^2, a + 2 beta >= 0, a^2 + 4c >= 0, l >= 0 and the \
                  scalar criterion, PF inputs with roots in [-R, 0], R = beta + sqrt(K), have \
                  PF images, sharply (the common-interlacer proof; the Chebyshev case is \
                  Iserles and Saff, Math. Proc. Cambridge Philos. Soc. 105 (1989), Theorem 1; \
                  the scalar criterion is RealRooted orderTwoConstants_nonneg_iff)"
            .into(),
        hypotheses: Vec::new(),
        conclusions: Vec::new(),
        counterexample: None,
        lean: vec!["RealRooted.orderTwoConstants_nonneg_iff"],
    };
    analysis.hypotheses.push(hypothesis(
        "characteristic_parameter_positive",
        format!("K = beta^2 + a beta - c = {} > 0", text(&k)),
        k.is_positive(),
        "exact".into(),
    ));
    analysis.hypotheses.push(hypothesis(
        "initial_slope_nonnegative",
        format!("eta = m/d = {} >= 0", text(&eta)),
        !eta.is_negative(),
        "exact".into(),
    ));
    if analysis
        .hypotheses
        .iter()
        .any(|h| h.status != LiftHypothesisStatus::Verified)
    {
        return analysis;
    }
    let root_k = QuadraticNumber::square_root(k.clone()).expect("K > 0");
    let source = AlgebraicInterval::closed(
        root_k.negated().add_rational(&-beta.clone()),
        root_k.add_rational(&-beta.clone()),
    );
    analysis.push_conclusion(BasisMapConclusion {
        input_interval: source.clone(),
        output_interval: AlgebraicInterval::whole_line(),
        nonnegative_coefficients: false,
        sharp: false,
        statement: format!(
            "every real polynomial with all roots in E = {source} has a real-rooted image"
        ),
        hypotheses: Vec::new(),
    });
    // Target interval I in t: ((-a - 2 beta) -+ 2 sqrt(K))/d.
    let center = (-a - integer(2) * &beta) / d;
    let half_width = root_k.scale(&(integer(2) / d));
    let target = AlgebraicInterval::closed(
        half_width.negated().add_rational(&center),
        half_width.add_rational(&center),
    );
    let mut target_proved = false;
    if eta.is_positive() {
        let two_eta = integer(2) * &eta - Q::one();
        let deviation = &s - &two_eta * &beta;
        let holds = !two_eta.is_negative() && &deviation * &deviation <= &two_eta * &two_eta * &k;
        analysis.hypotheses.push(hypothesis(
            "target_interval",
            format!(
                "2 eta - 1 = {} >= 0 and (s - (2 eta - 1) beta)^2 <= (2 eta - 1)^2 K, s = {}",
                text(&two_eta),
                text(&s)
            ),
            holds,
            "exact".into(),
        ));
        if holds {
            target_proved = true;
            analysis.push_conclusion(BasisMapConclusion {
                input_interval: source.clone(),
                output_interval: target.clone(),
                nonnegative_coefficients: false,
                sharp: false,
                statement: format!(
                    "every real polynomial with all roots in E = {source} has an image with all \
                     roots in I = {target}"
                ),
                hypotheses: Vec::new(),
            });
        }
    }
    // PF on [-R, 0].
    let radius = root_k.add_rational(&beta);
    let checks: Vec<(&str, String, bool)> = vec![
        (
            "current_constant_nonnegative",
            format!("a = {} >= 0", text(a)),
            !a.is_negative(),
        ),
        (
            "zero_in_source",
            format!("K >= beta^2 ({} >= {})", text(&k), text(&(&beta * &beta))),
            k >= &beta * &beta,
        ),
        (
            "radius_positive",
            format!("R = beta + sqrt(K) = {radius} > 0"),
            radius.signum() == Ordering::Greater,
        ),
        (
            "target_upper_end_nonpositive",
            format!(
                "a + 2 beta = {} >= 0 and a^2 + 4c = {} >= 0",
                text(&(a + integer(2) * &beta)),
                text(&(a * a + integer(4) * c))
            ),
            !(a + integer(2) * &beta).is_negative() && !(a * a + integer(4) * c).is_negative(),
        ),
        (
            "initial_constant_nonnegative",
            format!("l = {} >= 0", text(l)),
            !l.is_negative(),
        ),
    ];
    for (name, statement, holds) in checks {
        analysis
            .hypotheses
            .push(hypothesis(name, statement, holds, "exact".into()));
    }
    let (criterion, evidence) = scalar_criterion(a, c, l);
    analysis.hypotheses.push(hypothesis(
        "constant_terms_nonnegative",
        "2l >= a or c >= l^2 - a l (with a^2 + 4c >= 0): every P_k(0) >= 0".into(),
        criterion,
        evidence,
    ));
    let pf_names = [
        "current_constant_nonnegative",
        "zero_in_source",
        "radius_positive",
        "target_upper_end_nonpositive",
        "initial_constant_nonnegative",
        "constant_terms_nonnegative",
    ];
    let pf = analysis
        .hypotheses
        .iter()
        .filter(|h| pf_names.contains(&h.name.as_str()))
        .all(|h| h.status == LiftHypothesisStatus::Verified);
    if pf {
        let input =
            AlgebraicInterval::closed(radius.negated(), QuadraticNumber::rational(Q::zero()));
        let output = if target_proved {
            target.intersection(&AlgebraicInterval::nonpositive_ray())
        } else {
            AlgebraicInterval::nonpositive_ray()
        };
        analysis.push_conclusion(BasisMapConclusion {
            input_interval: input.clone(),
            output_interval: output,
            nonnegative_coefficients: normalization.is_positive(),
            sharp: false,
            statement: format!(
                "every PF input with all roots in {input} has a PF image; no larger [-R, 0] works"
            ),
            hypotheses: Vec::new(),
        });
    }
    analysis
}
