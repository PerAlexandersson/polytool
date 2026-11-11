//! Iterated real-rootedness preservers: first-order rules
//! `P_n = T_n(P_(n-1))` whose step `T_n` is a multiplier `a_n(x)` after an
//! operator `S_n` of the preserver menu (or after the identity).
//!
//! # The theorem
//!
//! Let `R` be a closed region of the real line (here `(-inf, 0]`,
//! `[0, inf)` or the whole line), `m` the induction start (the last seed row)
//! and suppose, for every `n > m`:
//!
//! 1. `T_n = a_n(x) S_n` as operators, where `S_n` maps every polynomial that
//!    is zero or real-rooted with every root in `R` to zero or to a
//!    real-rooted polynomial with every root in `R`;
//! 2. `a_n(x)` is real-rooted with every root in `R`;
//! 3. `T_n` raises the degree of `P_(n-1)` by exactly `delta >= 0`.
//!
//! If `P_m` is nonzero and real-rooted with every root in `R`, then by
//! induction every `P_n`, `n >= m`, is: `S_n(P_(n-1))` is zero or real-rooted
//! with its roots in `R` by 1, the roots of the product `a_n S_n(P_(n-1))`
//! are those of its factors (2), and by 3 the product is not zero.  The seed
//! rows before `m` are checked exactly one by one.  This is the elementary
//! product fact for `S_n = id`, Laguerre's theorems for multiplier sequences
//! `f(theta)` and differential operators `phi(D)`, Rolle's theorem for
//! `x^m D^j` and `D^j x^m`, and Brändén's Euler step (the preserver menu of
//! the preserver-image lifts, with the same theorems, side conditions and
//! regions; see [`find_preserver_image_lifts`]).
//!
//! # What is checked, exactly
//!
//! The rule must be `f(n) P_n = sum_d c_d(n, x) D^d P_(n-1)` with `f` free of
//! `x` (after the normalization of a left-hand factor `x^e f(n)`), without
//! `(-1)^n` factors or an inhomogeneous term.  Write
//! `T_n = (1/f(n)) sum_d c_d(n, x) D^d` and `N = m + 1`.
//!
//! - **Recognition.** `g(x)` is the monic gcd over `Q[x]` of every
//!   coefficient of `n^i` in every `c_d` (the multiplier part free of `n`).
//!   For every split `g = x^k g_k` the operator `x^k (1/g) sum c_d D^d` is
//!   matched against the menu: the identity (no derivatives; the multiplier
//!   is then `a_n = c_0 / f`), `phi(n, D)` (coefficients free of `x`),
//!   `f(n, theta)` (Euler form with one power of `x`), the Euler step
//!   `kappa (theta + a)(theta + b) + x (u + v theta)`, `kappa(n) x^m D^j` and
//!   `kappa(n) D^j x^m`; the index factors of a match go into the
//!   multiplier `a_n = g_k(x) h(n, x)`.  Every match is accepted only when
//!   `f(n) q(n) delta_S(n) T_n = f(n) (g_k H) o S_num` holds exactly as an
//!   identity of operators in the Ore algebra (`ore-algebra` composition;
//!   `q` the common denominator of `h = H / q`, `S_n = S_num / delta_S`).
//!   [`IteratedPreserverCertificate::verify_operator_identity`] replays it.
//! - **Hypothesis 1.** The side conditions of the menu theorem for every
//!   `n >= N` (decided exactly by the sign checker, as for the lifts), the
//!   nonvanishing of the operator's coefficient denominator, and the region
//!   transfer: for `(-inf, 0]` the menu image of `(-inf, 0]` lies in
//!   `(-inf, 0]`; for the whole line the menu theorem holds on the line;
//!   for `[0, inf)` the operator is `f(theta)`, `x^m D^j` or `D^j x^m`,
//!   which commute up to sign with `x -> -x`, and the menu image of
//!   `(-inf, 0]` lies in `(-inf, 0]`.
//! - **Hypothesis 2.** The factor `g_k` free of `n` is checked exactly; the
//!   index factor `h(n, x)` by the menu's polynomial test (fixed
//!   polynomials exactly, index-dependent ones of degree `<= 2` by the
//!   discriminant and the pairwise coefficient signs with the sign checker;
//!   for `[0, inf)` applied to `h(n, -x)`; higher degrees stay open).
//! - **Hypothesis 3.** With `delta = max_d (deg_x c_d - d)` and
//!   `Lambda(n, e) = sum_(deg_x c_d - d = delta) lead_x c_d(n) e (e - 1) ...
//!   (e - d + 1)` (the coefficient of `x^(e + delta)` in
//!   `f(n) T_n(x^e)`), `Lambda(n, deg P_m + delta (n - 1 - m)) != 0` for every
//!   `n >= N` by the sign checker; then `deg P_n = deg P_(n-1) + delta` and no
//!   row vanishes.
//! - **Factored steps.** If the existing checks do not prove a lag-one
//!   second-order rule, a bounded exact search over `Q(n)[x]` tries a
//!   product of two first-order factors after their common multiplier. The
//!   certificate records the factors in application order. Each uses the
//!   existing menu or the shared Liu--Wang degree-one conditions, with the
//!   actual intermediate degree, on a common closed region. Continuity
//!   extends the negative-root step to the closed ray, and reflection gives
//!   the positive ray. The intermediate leading coefficients and degrees
//!   are certified separately, so every intermediate image is nonzero.
//!   Search bounds and the quadratic Riccati remainder are documented in
//!   `doc/library.md`; failure to find a factorization is not a disproof.
//!   Here hypothesis 1 is needed only on polynomials of the prescribed
//!   induction degree: the Liu--Wang step is degree-local, and the degree
//!   of each factor's input is certified before invoking its theorem.
//! - The rows `P_first, ..., P_m` are nonzero and real-rooted with every root
//!   in `R` (exact).
//!
//! The regions are tried in the order `(-inf, 0]`, `[0, inf)`, line; the
//! first proved one is reported.  A requested root region is read open at
//! its finite ends (as the report records it) and is met only when it
//! contains the certified closed region.
//!
//! # Brändén's Weyl step
//!
//! The menu ends with the general second-order operator
//! `S_n = A + B D + C D^2` (recorded as [`PreserverOperator::WeylOperator`]),
//! a real-rootedness preserver on the whole line under the hypotheses of
//! Brändén's Corollary 3.5 (P. Brändén, Trans. Amer. Math. Soc. 358 (2006)):
//! `B^2 - 4AC >= 0` on the line, `A` strictly interlacing or strictly
//! alternating left of `B` with leading coefficients of the same sign, and
//! `S_n(x^k)` of degree `deg A + k` with leading coefficients of one sign for
//! every `k >= 0`, each decided for every `n >= N` (exactly for fixed
//! operators, by the sign checker for index-dependent ones; see
//! [`check_weyl_step_preserver`] for the fixed case).  It keeps `(-inf, 0]`
//! (respectively `[0, inf)`) when every coefficient of `A`, `B`, `C` (of the
//! reflected operator) has one sign, since a real-rooted polynomial with
//! coefficients of one sign has its roots in `(-inf, 0]`.  A negative
//! discriminant refutes the operator as a global preserver (Borcea and
//! Brändén 2010, Theorem 1.3, necessary condition); the earlier menu
//! operators, in particular the half-line Euler step, are tried first.
//!
//! # Lags and coordinates
//!
//! The rule may also have a single lag `L >= 1`, `f(n) P_n = sum_d c_d(n, x)
//! D^d P_(n-L)`: the induction then runs in each residue class of `n`
//! modulo `L` from the last `L` seed rows `P_(N-L), ..., P_(N-1)`, which
//! must all be nonzero and in the region; the degree step is decided class
//! by class (`n = s + L t`, `t >= 1`, for the class start `s`, by the sign
//! checker in `t`), and the degree-local Euler condition uses an upper bound
//! of `deg P_(n-L)` valid in every class.  The induction may also run on
//! rows `Q_n` in other coordinates ([`IteratedPreserverCoordinates`]: the
//! power gauge `P_n / x^e`, a translation `P_n(x + r)`, the reversal at the
//! nominal degree), tried after the identity, with the region of `P_n`
//! derived from that of `Q_n`; the seed rows before the class starts are
//! checked in `P` coordinates.
//!
//! # Interlacing
//!
//! For a pure multiplier of degree exactly one (`S_n = id`, `delta = 1`,
//! lag one, identity coordinates): `P_n = a_n P_(n-1)` with `a_n` linear
//! adds one real root, and inserting one root into a list of roots gives the
//! weak interlacing `P_(n-1) ⪯ P_n` of Liu and Wang
//! (`... <= s_2 <= r_2 <= s_1 <= r_1`) for every `n >= N`.  For a step
//! `T = g(x) S` free of `n` with `S` a proved Weyl step (lag one, identity
//! coordinates, `delta = 1`): `T` preserves real-rootedness, hence the
//! alternation of two polynomials (Brändén 2006, Corollary 2.5), so when
//! `P_(N-1) ⪯ P_N` is checked exactly, every pair `P_(n-1), P_n` with
//! `n >= N` alternates, which for degrees `d, d + 1` is the weak
//! interlacing `P_(n-1) ⪯ P_n`.  Both are extended over the seed rows only
//! as far as they are checked exactly.  No interlacing is claimed otherwise.

use super::*;
use crate::polynomial::Polynomial;
use crate::recurrence::implicit::algebra::{to_theta, Biv, DOp};
use crate::recurrence::implicit::biv_json;
use crate::recurrence::lift::lift_summary::{nonpositive_ray, real_line};
use crate::recurrence::lift::preserver_image::{
    operator_transfer, polynomial_hypotheses, region_contains,
};
use crate::recurrence::lift::veronese_algebra::classify_row;
use crate::recurrence::{IndexRationalFunction, OperatorTermJson, PreserverOperator};

mod coordinates;
mod factorization;
mod weyl_step;
pub use coordinates::IteratedPreserverCoordinates;
pub use factorization::IteratedPreserverFactor;
pub use weyl_step::{check_weyl_step_preserver, WeylStepPreserverVerdict};
use weyl_step::{
    weyl_step_hypotheses, weyl_step_parts, weyl_step_region_hypotheses, weyl_step_text,
    BRANDEN_WEYL_STEP,
};

const REFERENCE: &str = "Iterated real-rootedness preservers (induction on n): if P_m is \
     nonzero and real-rooted with every root in a closed region R and, for every n > m, \
     P_n = a_n(x) S_n(P_(n-1)) where S_n maps every polynomial that is zero or real-rooted with \
     every root in R to zero or a real-rooted polynomial with every root in R and a_n is \
     real-rooted with every root in R, then every P_n is real-rooted with every root in R (the \
     roots of a product are the roots of its factors); an exact degree step keeps every row \
     nonzero";

/// The coefficients `c_d(n, x)` of `D^d P_(n-1)`, by derivative order.
type StepCoefficients = BTreeMap<usize, Bivariate>;

/// A rational function `numerator(n) / denominator(n)`.
type IndexFraction = (Vec<Q>, Vec<Q>);

/// The record of the recognized step `T_n = a_n(x) S_n` of a proved (or
/// attempted) iterated-preserver check (module documentation of
/// [`check_iterated_preserver_theorem`]).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IteratedPreserverCertificate {
    /// The induction start `m`: the last seed row (after a moved start, the
    /// last moved seed row); the rule applies from `m + 1` on.
    pub induction_start: usize,
    /// `deg P_m`.
    pub start_degree: usize,
    /// `delta`: `deg P_n = deg P_(n-1) + delta` for every `n > m`.
    pub degree_step: usize,
    /// The multiplier factor `g(x)` free of `n` (ascending exact rational
    /// strings; `["1"]` when there is none).
    pub multiplier_fixed_factor: Vec<String>,
    /// The multiplier factor `h(n, x) = sum_j h_j(n) x^j` (coefficient `j`
    /// at position `j`); the multiplier is `a_n(x) = g(x) h(n, x)`.
    pub multiplier_index_factor: Vec<IndexRationalFunction>,
    /// The menu operator `S_n`; omitted for the identity or when `factors`
    /// records the first-order product instead.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub preserver: Option<PreserverOperator>,
    /// First-order factors in application order (rightmost first). Empty
    /// for the original single-menu-operator certificate. Their product
    /// replaces `preserver`, which must then be absent.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub factors: Vec<IteratedPreserverFactor>,
    /// The closed region `R` containing every root of every row (`lower`
    /// and `upper` omitted for the whole line).
    pub region: RootInterval,
    /// The operator identity in words.
    pub operator_identity: String,
    /// The lag `L` of the rule `P_n = T_n(P_(n-L))` (module documentation,
    /// "Lags and coordinates"); `1` when omitted.
    #[serde(default = "lag_one", skip_serializing_if = "is_lag_one")]
    pub lag: usize,
    /// The coordinates `Q_n` the induction runs on; the identity when
    /// omitted.  `start_degree`, `degree_step`, the multiplier and the
    /// preserver refer to the rows `Q_n`, `region` to the rows `P_n`.
    #[serde(
        default,
        skip_serializing_if = "IteratedPreserverCoordinates::is_identity"
    )]
    pub coordinates: IteratedPreserverCoordinates,
}

fn lag_one() -> usize {
    1
}

fn is_lag_one(lag: &usize) -> bool {
    *lag == 1
}

impl IteratedPreserverCertificate {
    /// Replay the recognition: the step operator of `rule` (the rule the
    /// check ran on, after a left-hand factor normalization), transformed to
    /// the recorded coordinates, equals `a_n(x) S_n` exactly in the Ore
    /// algebra, the lag and the degree step are those of the rule, and the
    /// region is one of the certified closed regions.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::VerificationFailed`] when the rule is not a
    /// single-lag rule of the theorem's form, the identity does not hold,
    /// the lag or the degree step differs or the region is not supported;
    /// [`RecurrenceError::InvalidInput`] when a recorded value does not parse.
    pub fn verify_operator_identity(&self, rule: &Recurrence) -> Result<(), RecurrenceError> {
        let (terms, f, lag) = step_operator(rule).map_err(RecurrenceError::VerificationFailed)?;
        if lag != self.lag {
            return Err(RecurrenceError::VerificationFailed(
                "the recorded lag differs from the rule's".into(),
            ));
        }
        let terms = self.coordinates.transform_step(&terms).ok_or_else(|| {
            RecurrenceError::VerificationFailed(
                "the step operator has no polynomial form in the recorded coordinates".into(),
            )
        })?;
        let decomposition = Decomposition::from_certificate(self)?;
        if !decomposition.identity_holds(&terms, &f) {
            return Err(RecurrenceError::VerificationFailed(
                "the step operator of the rule differs from a_n(x) S_n".into(),
            ));
        }
        if degree_step(&terms) != Some(self.degree_step as i64) {
            return Err(RecurrenceError::VerificationFailed(
                "the recorded degree step differs from the rule's".into(),
            ));
        }
        if ClosedRegion::of_root_interval(&self.region).is_none() {
            return Err(RecurrenceError::VerificationFailed(
                "the recorded region is not (-inf, 0], [0, inf) or the real line".into(),
            ));
        }
        Ok(())
    }
}

/// Check the iterated-preserver theorem (module documentation, summarized
/// here) for the family given by `rule`, `first_index` and the seed rows
/// `seeds` (the rule applies from `N = first_index + seeds.len()` on; the
/// induction starts at the last seed row `P_(N-1)`).
///
/// The rule must have a single lag, `f(n) P_n = sum_d c_d(n, x) D^d
/// P_(n-L)` with `f` free of `x` (module documentation, "Lags and
/// coordinates"; the coordinates other than the identity are tried after
/// it).  Its step operator is recognized exactly as
/// `T_n = a_n(x) S_n` with `S_n` the identity or an operator of the
/// preserver menu ([`PreserverOperator`]: `phi(n, D)`, `f(n, theta)`, the
/// Euler step, `x^m D^j`, `D^j x^m`, and Brändén's Weyl step
/// `A + B D + C D^2`); then, for a closed region `R` among
/// `(-inf, 0]`, `[0, inf)` and the line, the hypotheses are: the side
/// conditions and the region transfer of `S_n` for every `n >= N`; `a_n`
/// real-rooted with every root in `R` for every `n >= N` (fixed factors
/// exactly, index-dependent factors of degree `<= 2` by the sign checker);
/// an exact degree step (sign checker); and the seed rows nonzero and
/// real-rooted with every root in `R` (exact).  When all hold, every
/// `P_n`, `n >= first_index`, is nonzero and real-rooted with every root in
/// `R`.  The induction start is moved past the seeds as for the other
/// theorem checks ([`RealRootednessTheoremOptions::maximum_start_delay`]).
/// If those checks do not prove a lag-one second-order step, a bounded
/// exact right factorization over `Q(n)[x]` tries two first-order factors.
/// Each uses the existing menu or the shared Liu--Wang degree-one
/// conditions at its certified input degree, on a common closed region.
/// [`IteratedPreserverCertificate::factors`] records the application order;
/// intermediate images are certified nonzero. See `doc/library.md` for
/// the deterministic search bounds and remaining unsearched forms.
///
/// The report carries the recognized step in
/// [`RealRootednessTheoremReport::iterated_preserver_certificate`];
/// [`RealRootednessTheoremReport::verify`] recomputes the report and replays
/// the operator identity and every sign certificate.
///
/// # Errors
///
/// Errors of the record itself (seed rows that cannot be normalized, a rule
/// that cannot generate moved seed rows).
///
/// # Example
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::recurrence::{
///     check_iterated_preserver_theorem, BigRational, BivarPoly, Recurrence, RecurrenceTerm,
///     RealRootednessTheoremOptions,
/// };
///
/// let q = |v: i64| BigRational::from_integer(BigInt::from(v));
/// // The 3-Pascal triangle: P_n = (1 + x) P_(n-1) from n = 3 on, P_2 = 1 + 3x + x^2.
/// let rule = Recurrence {
///     terms: vec![RecurrenceTerm::new(
///         1,
///         0,
///         BivarPoly::from_coefficients_normalized(vec![vec![q(1), q(1)]]),
///     )],
///     denominator: None,
///     inhomogeneous: None,
/// };
/// let seeds = vec![vec![q(1)], vec![q(1), q(1)], vec![q(1), q(3), q(1)]];
/// let report = check_iterated_preserver_theorem(
///     &rule,
///     0,
///     &seeds,
///     &RealRootednessTheoremOptions::default(),
/// )
/// .unwrap();
/// assert!(report.is_proved());
/// report.verify().unwrap();
/// ```
pub fn check_iterated_preserver_theorem(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<BigRational>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    check_real_rootedness_theorem(
        RealRootednessTheorem::IteratedPreserver,
        rule,
        first_index,
        seeds,
        options,
    )
}

// ---------------------------------------------------------------------------
// Regions
// ---------------------------------------------------------------------------

/// The closed regions the check certifies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ClosedRegion {
    NonpositiveRay,
    NonnegativeRay,
    RealLine,
}

impl ClosedRegion {
    const ALL: [Self; 3] = [Self::NonpositiveRay, Self::NonnegativeRay, Self::RealLine];

    fn root_interval(self) -> RootInterval {
        match self {
            Self::NonpositiveRay => nonpositive_ray(),
            Self::NonnegativeRay => RootInterval {
                lower: Some(RootEndpoint::rational(&Q::zero(), true)),
                upper: None,
            },
            Self::RealLine => real_line(),
        }
    }

    fn of_root_interval(region: &RootInterval) -> Option<Self> {
        Self::ALL
            .into_iter()
            .find(|candidate| &candidate.root_interval() == region)
    }

    fn text(self) -> &'static str {
        match self {
            Self::NonpositiveRay => "(-inf, 0]",
            Self::NonnegativeRay => "[0, inf)",
            Self::RealLine => "the real line",
        }
    }

    /// Whether the nonzero `p` is real-rooted with every root in the region
    /// (exact).
    fn contains_roots_of(self, p: &[Q]) -> bool {
        let p = vec_poly::trim(p);
        if vec_poly::is_zero(&p) {
            return false;
        }
        match self {
            Self::NonpositiveRay => classify_row(&p).has_nonpositive_real_roots(),
            Self::NonnegativeRay => classify_row(&reflected(&p)).has_nonpositive_real_roots(),
            Self::RealLine => p.len() <= 1 || is_real_rooted_rational(&p),
        }
    }

    /// Whether the region lies in `requested`.
    fn lies_in(self, requested: &RealInterval) -> bool {
        let upper_ok = |upper: &Option<Q>, closed: bool| match upper {
            None => true,
            Some(u) => u.is_positive() || (u.is_zero() && closed),
        };
        match self {
            Self::NonpositiveRay => {
                requested.lower.is_none() && upper_ok(&requested.upper, requested.upper_closed)
            }
            Self::NonnegativeRay => {
                requested.upper.is_none()
                    && upper_ok(
                        &requested.lower.as_ref().map(|l| -l),
                        requested.lower_closed,
                    )
            }
            Self::RealLine => requested.lower.is_none() && requested.upper.is_none(),
        }
    }
}

fn is_zero_row(p: &[Q]) -> bool {
    vec_poly::is_zero(&vec_poly::trim(p))
}

/// The first index of the conclusion: after the last zero seed row before
/// the class starts (the start rows themselves must be nonzero).
fn claim_first_index(first_index: usize, seeds: &[Vec<Q>], lag: usize) -> usize {
    seeds[..seeds.len() - lag]
        .iter()
        .rposition(|row| is_zero_row(row))
        .map_or(first_index, |offset| first_index + offset + 1)
}

/// `p(-x)`.
fn reflected(p: &[Q]) -> Vec<Q> {
    p.iter()
        .enumerate()
        .map(|(k, c)| if k % 2 == 1 { -c } else { c.clone() })
        .collect()
}

// ---------------------------------------------------------------------------
// The step operator and its decompositions
// ---------------------------------------------------------------------------

/// The step `f(n) T_n = sum_d c_d(n, x) D^d` and the lag `L` of a
/// single-lag rule `P_n = T_n(P_(n-L))`, or why the rule does not have that
/// form.
fn step_operator(rule: &Recurrence) -> Result<(StepCoefficients, Vec<Q>, usize), String> {
    let shape = recurrence_shape(rule)?;
    let Some(&(lag, _)) = shape.terms.keys().next() else {
        return Err("the rule has no nonzero term".into());
    };
    if lag == 0 || shape.terms.keys().any(|(offset, _)| *offset != lag) {
        return Err(
            "the rule is not first order: P_n = T_n(P_(n-1)), or T_n(P_(n-L)) with one lag L"
                .into(),
        );
    }
    let terms = shape
        .terms
        .into_iter()
        .map(|((_, order), coefficients)| (order, coefficients))
        .collect();
    Ok((terms, shape.denominator, lag))
}

/// `delta = max_d (deg_x c_d - d)`.
fn degree_step(terms: &StepCoefficients) -> Option<i64> {
    terms
        .iter()
        .filter_map(|(d, c)| x_degree(c).map(|degree| degree as i64 - *d as i64))
        .max()
}

/// `sum_d c_d D^d` as an Ore operator.
fn ore_operator(terms: &StepCoefficients) -> DOp {
    let mut op = DOp::default();
    for (d, c) in terms {
        op.add_term((0, *d), Biv::new(c.clone()));
    }
    op
}

/// A recognized step `T_n = g(x) h(n, x) S_n`.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Decomposition {
    /// `g`, nonzero, free of `n`.
    fixed_factor: Vec<Q>,
    /// `h_j(n)`, ascending in `x`.
    index_factor: Vec<IndexFraction>,
    /// `S_n`; `None` for the identity.
    preserver: Option<PreserverOperator>,
    factors: Vec<IteratedPreserverFactor>,
}

impl Decomposition {
    fn from_certificate(
        certificate: &IteratedPreserverCertificate,
    ) -> Result<Self, RecurrenceError> {
        let fixed_factor = certificate
            .multiplier_fixed_factor
            .iter()
            .map(|value| {
                parse_rational_coeff(value)
                    .map_err(|error| RecurrenceError::from(error).with_context("multiplier"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let index_factor = certificate
            .multiplier_index_factor
            .iter()
            .map(IndexRationalFunction::parse)
            .collect::<Result<Vec<_>, _>>()?;
        if vec_poly::is_zero(&fixed_factor) || index_factor.is_empty() {
            return Err(RecurrenceError::InvalidInput(
                "the recorded multiplier is empty or zero".into(),
            ));
        }
        Ok(Self {
            fixed_factor,
            index_factor,
            preserver: certificate.preserver.clone(),
            factors: certificate.factors.clone(),
        })
    }

    /// `S_n = numerator / denominator(n)`.
    fn preserver_parts(&self) -> Option<(DOp, Vec<Q>)> {
        if !self.factors.is_empty() {
            if self.preserver.is_some() {
                return None;
            }
            return factorization::product_parts(&self.factors);
        }
        match &self.preserver {
            None => Some((DOp::identity(), vec![Q::one()])),
            Some(operator) => operator.weyl_parts(),
        }
    }

    /// `h = H / q` with a common denominator `q`.
    fn index_factor_over_common_denominator(&self) -> (Biv, Vec<Q>) {
        let common = self
            .index_factor
            .iter()
            .fold(vec![Q::one()], |acc, (_, q)| vec_poly::mul(&acc, q));
        let mut numerator = Biv::zero();
        for (j, (p, _)) in self.index_factor.iter().enumerate() {
            let others = self
                .index_factor
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != j)
                .fold(p.clone(), |acc, (_, (_, q))| vec_poly::mul(&acc, q));
            numerator = numerator
                .add(&Biv::from_index_polynomial(&others).multiply(&Biv::monomial(0, j, Q::one())));
        }
        (numerator, common)
    }

    /// `q(n) delta_S(n) sum_d c_d D^d = f(n) (g H) o S_num` in the Ore
    /// algebra (the step `T_n = (1/f) sum c_d D^d` equals `g h S_n`).
    fn identity_holds(&self, terms: &StepCoefficients, f: &[Q]) -> bool {
        let Some((numerator, denominator)) = self.preserver_parts() else {
            return false;
        };
        let (h, common) = self.index_factor_over_common_denominator();
        let multiplier = Biv::from_variable_polynomial(&self.fixed_factor).multiply(&h);
        if multiplier.is_zero() || vec_poly::is_zero(&common) || vec_poly::is_zero(&denominator) {
            return false;
        }
        let left = ore_operator(terms).left_multiply(&Biv::from_index_polynomial(&vec_poly::mul(
            &common,
            &denominator,
        )));
        let right = numerator
            .left_multiply(&multiplier)
            .left_multiply(&Biv::from_index_polynomial(f));
        left.subtract(&right).is_zero()
    }

    fn index_factor_functions(&self, reflect: bool) -> Vec<IndexRationalFunction> {
        self.index_factor
            .iter()
            .enumerate()
            .map(|(j, (p, q))| {
                let p = if reflect && j % 2 == 1 {
                    vec_poly::neg(p)
                } else {
                    p.clone()
                };
                IndexRationalFunction::from_polynomials(&p, q)
            })
            .collect()
    }

    fn multiplier_text(&self) -> String {
        let index_part: Vec<String> = self
            .index_factor_functions(false)
            .iter()
            .enumerate()
            .filter(|(_, c)| c.text() != "0")
            .map(|(j, c)| match j {
                0 => format!("({})", c.text()),
                1 => format!("({}) x", c.text()),
                _ => format!("({}) x^{j}", c.text()),
            })
            .collect();
        format!(
            "({}) ({})",
            format_index_polynomial(&self.fixed_factor, "x"),
            index_part.join(" + ")
        )
    }
}

/// The monic gcd over `Q[x]` of every coefficient of `n^i` in every `c_d`.
fn fixed_content(terms: &StepCoefficients) -> Vec<Q> {
    let mut content: Option<Polynomial<Q>> = None;
    for c in terms.values() {
        for row in c {
            let row = vec_poly::trim(row);
            if vec_poly::is_zero(&row) {
                continue;
            }
            let row = Polynomial::new(row);
            content = Some(match content {
                None => row.gcd(&row),
                Some(g) => g.gcd(&row),
            });
        }
    }
    content.map_or_else(|| vec![Q::one()], |g| g.coeffs().to_vec())
}

/// `c / g` row by row (exact: `g` divides every row).
fn divide_rows(c: &Bivariate, g: &[Q]) -> Bivariate {
    let divisor = Polynomial::new(g.to_vec());
    c.iter()
        .map(|row| {
            let row = vec_poly::trim(row);
            if vec_poly::is_zero(&row) {
                return Vec::new();
            }
            Polynomial::new(row).div_rem(&divisor).0.coeffs().to_vec()
        })
        .collect()
}

/// `x^k c`.
fn times_x_power(c: &Bivariate, k: usize) -> Bivariate {
    c.iter()
        .map(|row| {
            let row = vec_poly::trim(row);
            if vec_poly::is_zero(&row) {
                Vec::new()
            } else {
                vec_poly::shift_by(&row, k)
            }
        })
        .collect()
}

/// The coefficient of `y^j` of a [`Biv`] as a polynomial in `n`.
fn biv_column(b: &Biv, j: usize) -> Vec<Q> {
    let mut column: Vec<Q> = b
        .rows()
        .iter()
        .map(|row| row.get(j).cloned().unwrap_or_else(Q::zero))
        .collect();
    vec_poly::strip_trailing_zeros(&mut column);
    if column.is_empty() {
        vec![Q::zero()]
    } else {
        column
    }
}

fn over(p: Vec<Q>, f: &[Q]) -> IndexRationalFunction {
    let p = if p.is_empty() { vec![Q::zero()] } else { p };
    IndexRationalFunction::from_polynomials(&p, f)
}

/// The single power `m` of `x` in `c`, with its coefficient in `n`.
fn single_power(c: &Bivariate) -> Option<(usize, Vec<Q>)> {
    let powers: Vec<usize> = (0..x_width(c))
        .filter(|&j| !x_column(c, j).is_empty())
        .collect();
    match powers.as_slice() {
        [m] => Some((*m, x_column(c, *m))),
        _ => None,
    }
}

/// The menu shapes of `s = sum_d s_d D^d` over `f(n)`: pairs of the index
/// factor `h` and the operator (module documentation, "Recognition").
fn menu_shapes(s: &StepCoefficients, f: &[Q]) -> Vec<(Vec<IndexFraction>, PreserverOperator)> {
    let unit = || vec![(vec![Q::one()], vec![Q::one()])];
    let mut out = Vec::new();
    let max_order = s.keys().copied().max().unwrap_or(0);
    if max_order == 0 {
        return out;
    }
    // phi(n, D).
    if s.values().all(|c| x_degree(c) == Some(0)) {
        let coefficients = (0..=max_order)
            .map(|d| over(s.get(&d).map(|c| x_column(c, 0)).unwrap_or_default(), f))
            .collect();
        out.push((
            unit(),
            PreserverOperator::DerivativePolynomial { coefficients },
        ));
    }
    // f(n, theta) and the Euler step, from the Euler form x^J s.
    let (j, euler) = to_theta(&ore_operator(s));
    let levels: Vec<usize> = euler.keys().map(|(_, level)| *level).collect();
    let theta_coefficients = |level: usize, count: usize| -> Vec<IndexRationalFunction> {
        let b = euler.get(&(0, level)).cloned().unwrap_or_else(Biv::zero);
        (0..count).map(|t| over(biv_column(&b, t), f)).collect()
    };
    let theta_degree = |level: usize| {
        euler
            .get(&(0, level))
            .and_then(Biv::variable_degree)
            .unwrap_or(0)
    };
    if euler.keys().all(|(shift, _)| *shift == 0) {
        if levels.iter().all(|&level| level == j) {
            let count = theta_degree(j) + 1;
            out.push((
                unit(),
                PreserverOperator::EulerPolynomial {
                    coefficients: theta_coefficients(j, count),
                },
            ));
        } else if levels.iter().all(|&level| level == j || level == j + 1)
            && levels.contains(&(j + 1))
            && theta_degree(j) == 2
            && theta_degree(j + 1) <= 1
        {
            out.push((
                unit(),
                PreserverOperator::EulerStep {
                    diagonal: theta_coefficients(j, 3),
                    raising: theta_coefficients(j + 1, 2),
                },
            ));
        }
    }
    // kappa(n) x^m D^j and kappa(n) D^j x^m.
    if let Some((power, kappa)) = s.get(&max_order).and_then(single_power) {
        let factor = vec![(kappa, f.to_vec())];
        if s.len() == 1 {
            out.push((
                factor.clone(),
                PreserverOperator::PowerTimesDerivative {
                    power,
                    derivative_order: max_order,
                },
            ));
        }
        if power >= 1 && s.len() > 1 {
            out.push((
                factor,
                PreserverOperator::DerivativeOfPowerMultiple {
                    power,
                    derivative_order: max_order,
                },
            ));
        }
    }
    // Branden's Weyl step A + B D + C D^2 (last: the operators above have
    // their own theorems, and the Euler step is a half-line theorem).
    if max_order <= 2 {
        out.push((
            unit(),
            PreserverOperator::WeylOperator {
                terms: s
                    .iter()
                    .map(|(order, c)| OperatorTermJson {
                        deriv_order: *order,
                        coeff: biv_json(&Biv::new(c.clone())),
                    })
                    .collect(),
                denominator: vec_poly::trim(f)
                    .iter()
                    .map(format_rational_coeff)
                    .collect(),
            },
        ));
    }
    out
}

/// Every decomposition `T_n = g_k(x) h(n, x) S_n` whose operator identity
/// holds exactly, in a fixed order.
fn decompositions(terms: &StepCoefficients, f: &[Q]) -> Vec<Decomposition> {
    let g = fixed_content(terms);
    let quotient: StepCoefficients = terms
        .iter()
        .map(|(d, c)| (*d, divide_rows(c, &g)))
        .collect();
    let mut candidates = Vec::new();
    if terms.keys().all(|&d| d == 0) {
        let c0 = quotient.get(&0).cloned().unwrap_or_default();
        let index_factor = (0..=x_degree(&c0).unwrap_or(0))
            .map(|j| {
                let column = x_column(&c0, j);
                (
                    if column.is_empty() {
                        vec![Q::zero()]
                    } else {
                        column
                    },
                    f.to_vec(),
                )
            })
            .collect();
        candidates.push(Decomposition {
            fixed_factor: g.clone(),
            index_factor,
            preserver: None,
            factors: Vec::new(),
        });
    } else {
        let valuation = g.iter().take_while(|c| c.is_zero()).count();
        for k in 0..=valuation {
            let fixed_factor = g[k..].to_vec();
            let shifted: StepCoefficients = quotient
                .iter()
                .map(|(d, c)| (*d, times_x_power(c, k)))
                .collect();
            for (index_factor, preserver) in menu_shapes(&shifted, f) {
                candidates.push(Decomposition {
                    fixed_factor: fixed_factor.clone(),
                    index_factor,
                    preserver: Some(preserver),
                    factors: Vec::new(),
                });
            }
        }
    }
    let mut out: Vec<Decomposition> = Vec::new();
    for candidate in candidates {
        if candidate.identity_holds(terms, f) && !out.contains(&candidate) {
            out.push(candidate);
        }
    }
    out
}

// ---------------------------------------------------------------------------
// The check
// ---------------------------------------------------------------------------

/// The parts of the family the attempts in one coordinate system share.
struct Family<'a> {
    /// The step in the coordinates.
    terms: StepCoefficients,
    /// The step of the rows `P_n` (the rule's own).
    original_terms: &'a StepCoefficients,
    f: &'a [Q],
    lag: usize,
    coordinates: IteratedPreserverCoordinates,
    first_index: usize,
    /// The seed rows `P_first, ..., P_(N-1)` (normalized, in `P`
    /// coordinates).
    seeds: &'a [Vec<Q>],
    /// The class starts `P_(N-L), ..., P_(N-1)` in the coordinates.
    start_rows: Vec<Vec<Q>>,
    n_first: usize,
    start: usize,
    start_degree: usize,
    delta: Option<i64>,
    /// An upper bound of `deg Q_(n-L)` as a polynomial in `n` (exact for
    /// lag one; for `delta >= 0`).
    previous_degree: Option<Vec<Q>>,
    degree_hypotheses: Vec<LiftHypothesis>,
    requested: Option<&'a RealInterval>,
}

impl<'a> Family<'a> {
    /// The family in `coordinates`, or `None` when the step or a start row
    /// has no polynomial form there.
    #[allow(clippy::too_many_arguments)]
    fn new(
        original_terms: &'a StepCoefficients,
        f: &'a [Q],
        lag: usize,
        coordinates: IteratedPreserverCoordinates,
        first_index: usize,
        seeds: &'a [Vec<Q>],
        requested: Option<&'a RealInterval>,
    ) -> Option<Self> {
        let terms = coordinates.transform_step(original_terms)?;
        let start_rows = seeds[seeds.len() - lag..]
            .iter()
            .map(|row| coordinates.transform_row(row))
            .collect::<Option<Vec<_>>>()?;
        let n_first = first_index + seeds.len();
        let start = n_first - 1;
        let start_degree = vec_poly::degree(&start_rows[lag - 1]).unwrap_or(0);
        let delta = degree_step(&terms);
        let class_degrees: Vec<usize> = start_rows
            .iter()
            .map(|row| vec_poly::degree(row).unwrap_or(0))
            .collect();
        let previous_degree = delta
            .filter(|&delta| delta >= 0)
            .map(|delta| previous_degree_bound(n_first, lag, &class_degrees, delta));
        let degree_hypotheses = if lag == 1 {
            vec![degree_step_hypothesis(&terms, n_first, start, start_degree)]
        } else {
            (0..lag)
                .map(|class| {
                    class_degree_step_hypothesis(
                        &terms,
                        lag,
                        n_first - lag + class,
                        class_degrees[class],
                        class,
                    )
                })
                .collect()
        };
        Some(Self {
            terms,
            original_terms,
            f,
            lag,
            coordinates,
            first_index,
            seeds,
            start_rows,
            n_first,
            start,
            start_degree,
            delta,
            previous_degree,
            degree_hypotheses,
            requested,
        })
    }

    /// `P` for the identity coordinates, `Q` otherwise.
    fn row_letter(&self) -> &'static str {
        if self.coordinates.is_identity() {
            "P"
        } else {
            "Q"
        }
    }
}

/// `Lambda(n, deg P_(n-1)) != 0` for every `n >= N` (module documentation,
/// hypothesis 3).
fn degree_step_hypothesis(
    terms: &StepCoefficients,
    n_first: usize,
    start: usize,
    start_degree: usize,
) -> LiftHypothesis {
    let name = "degree_step";
    let statement = format!(
        "deg P_n = deg P_(n-1) + delta for every n >= {n_first}: the coefficient Lambda(n, e) of \
         x^(e + delta) in f(n) T_n(x^e) is nonzero at e = deg P_(n-1)"
    );
    let Some(delta) = degree_step(terms) else {
        return unresolved_hypothesis(name, &statement, "the step operator is zero".into());
    };
    if delta < 0 {
        return unresolved_hypothesis(
            name,
            &statement,
            format!("the step lowers the degree (delta = {delta}), so rows eventually vanish"),
        );
    }
    // e(n) = deg P_(n-1) = start_degree + delta (n - 1 - start).
    let previous_degree = previous_degree_polynomial(start, start_degree, delta);
    let lambda = leading_lambda(terms, delta, &[Q::zero(), Q::one()], &previous_degree);
    let mut hypothesis = sign_hypothesis(
        name,
        &statement,
        "Lambda(n, deg P_(n-1))",
        &IndexFunction::polynomial(lambda),
        n_first,
        IndexSignCondition::Nonzero,
    );
    hypothesis.evidence = format!(
        "delta = {delta}, deg P_(n-1) = {}: {}",
        format_index_polynomial(&previous_degree, "n"),
        hypothesis.evidence
    );
    hypothesis
}

/// `Lambda(n(t), e(t))` as a polynomial in `t`, with `n(t)` and `e(t)`
/// polynomials in `t`: the coefficient of `x^(e + delta)` in
/// `f(n) T_n(x^e)`.
fn leading_lambda(terms: &StepCoefficients, delta: i64, index: &[Q], degree: &[Q]) -> Vec<Q> {
    let mut lambda = vec![Q::zero()];
    for (d, c) in terms {
        let Some(x_power) = x_degree(c) else {
            continue;
        };
        if x_power as i64 - *d as i64 != delta {
            continue;
        }
        let mut product = vec_poly::compose(&x_column(c, x_power), index);
        for t in 0..*d {
            let factor = vec_poly::sub(degree, &[Q::from_integer(BigInt::from(t))]);
            product = vec_poly::mul(&product, &factor);
        }
        lambda = vec_poly::add(&lambda, &product);
    }
    lambda
}

/// The degree step of the residue class of the start row `class_start`
/// (lag `L >= 2`): `Lambda(s + L t, deg P_s + delta (t - 1)) != 0` for every
/// `t >= 1`.
fn class_degree_step_hypothesis(
    terms: &StepCoefficients,
    lag: usize,
    class_start: usize,
    class_degree: usize,
    class: usize,
) -> LiftHypothesis {
    let name = format!("degree_step_class_{class}");
    let statement = format!(
        "deg P_n = deg P_(n-{lag}) + delta for every n = {class_start} + {lag} t, t >= 1: the \
         coefficient Lambda(n, e) of x^(e + delta) in f(n) T_n(x^e) is nonzero at e = deg \
         P_(n-{lag}) = deg P_{class_start} + delta (t - 1)"
    );
    let Some(delta) = degree_step(terms) else {
        return unresolved_hypothesis(&name, &statement, "the step operator is zero".into());
    };
    if delta < 0 {
        return unresolved_hypothesis(
            &name,
            &statement,
            format!("the step lowers the degree (delta = {delta}), so rows eventually vanish"),
        );
    }
    let to_q = |value: i64| Q::from_integer(BigInt::from(value));
    let index = [to_q(class_start as i64), to_q(lag as i64)];
    let degree = vec_poly::trim(&[to_q(class_degree as i64 - delta), to_q(delta)]);
    let lambda = leading_lambda(terms, delta, &index, &degree);
    let mut hypothesis = sign_hypothesis(
        &name,
        &statement,
        "Lambda(n(t), e(t)) in t",
        &IndexFunction::polynomial(lambda),
        1,
        IndexSignCondition::Nonzero,
    );
    // The sign checker runs in t: report a failure at the row index
    // n = s + L t.
    if let LiftHypothesisStatus::Failed {
        witness: Some(witness),
    } = &mut hypothesis.status
    {
        *witness = BigInt::from(class_start) + BigInt::from(lag) * &*witness;
    }
    hypothesis.evidence = format!(
        "delta = {delta}, n = {class_start} + {lag} t, deg P_(n-{lag}) = {}: {}",
        format_index_polynomial(&degree, "t"),
        hypothesis.evidence.replace(" n = ", " t = ")
    );
    hypothesis
}

/// `e(n) = deg P_(n-1) = start_degree + delta (n - 1 - start)`.
fn previous_degree_polynomial(start: usize, start_degree: usize, delta: i64) -> Vec<Q> {
    let delta_q = Q::from_integer(BigInt::from(delta));
    vec_poly::trim(&[
        Q::from_integer(BigInt::from(start_degree))
            - &delta_q * Q::from_integer(BigInt::from(start + 1)),
        delta_q,
    ])
}

/// An upper bound of `deg P_(n-L)` valid in every residue class:
/// `max_c (D_c - delta (s_c + L)/L) + delta n / L` (class `c` starting at
/// `s_c = N - L + c` with degree `D_c`); for `L = 1` this is the exact
/// [`previous_degree_polynomial`].
fn previous_degree_bound(
    n_first: usize,
    lag: usize,
    class_degrees: &[usize],
    delta: i64,
) -> Vec<Q> {
    if lag == 1 {
        return previous_degree_polynomial(n_first - 1, class_degrees[0], delta);
    }
    let lag_q = Q::from_integer(BigInt::from(lag));
    let slope = Q::from_integer(BigInt::from(delta)) / &lag_q;
    let constant = class_degrees
        .iter()
        .enumerate()
        .map(|(c, &degree)| {
            Q::from_integer(BigInt::from(degree))
                - &slope * (Q::from_integer(BigInt::from(n_first - lag + c)) + &lag_q)
        })
        .max()
        .unwrap_or_else(Q::zero);
    vec_poly::trim(&[constant, slope])
}

/// The row hypotheses: the class starts in the working region `region`
/// (in the coordinates) and the earlier seeds in the certified region
/// `claimed` (in `P` coordinates).
fn row_hypotheses(
    family: &Family<'_>,
    region: ClosedRegion,
    claimed: ClosedRegion,
) -> Vec<LiftHypothesis> {
    let mut out = Vec::new();
    let start = family.start;
    let letter = family.row_letter();
    if family.lag == 1 {
        let seed = &family.start_rows[0];
        out.push(exact_hypothesis(
            "start_row_in_region",
            &format!(
                "{letter}_{start} is nonzero and real-rooted with every root in {}",
                region.text()
            ),
            if region.contains_roots_of(seed) {
                Ok(format!(
                    "{letter}_{start} has degree {}; checked exactly",
                    family.start_degree
                ))
            } else {
                Err((
                    Some(start),
                    format!(
                        "{letter}_{start} is zero, not real-rooted, or has a root outside {} \
                         (exact)",
                        region.text()
                    ),
                ))
            },
        ));
    } else {
        let first_start = family.n_first - family.lag;
        let failing = family
            .start_rows
            .iter()
            .position(|row| !region.contains_roots_of(row));
        out.push(exact_hypothesis(
            "start_rows_in_region",
            &format!(
                "the class starts {letter}_{first_start}, ..., {letter}_{start} are nonzero and \
                 real-rooted with every root in {}",
                region.text()
            ),
            match failing {
                None => Ok(format!(
                    "checked exactly for {letter}_{first_start}, ..., {letter}_{start}"
                )),
                Some(offset) => Err((
                    Some(first_start + offset),
                    format!(
                        "{letter}_{} is zero, not real-rooted, or has a root outside {} (exact)",
                        first_start + offset,
                        region.text()
                    ),
                )),
            },
        ));
    }
    let earlier = &family.seeds[..family.seeds.len() - family.lag];
    if !earlier.is_empty() {
        let first = family.first_index;
        let last = first + earlier.len() - 1;
        let failing = earlier
            .iter()
            .position(|row| !is_zero_row(row) && !claimed.contains_roots_of(row));
        let zero_rows: Vec<String> = earlier
            .iter()
            .enumerate()
            .filter(|(_, row)| is_zero_row(row))
            .map(|(offset, _)| format!("P_{}", first + offset))
            .collect();
        let zero_note = if zero_rows.is_empty() {
            String::new()
        } else {
            format!(
                "; {} zero (the conclusion starts after the last zero seed row)",
                zero_rows.join(", ")
            )
        };
        out.push(exact_hypothesis(
            "earlier_seeds_in_region",
            &format!(
                "every nonzero seed row among P_{first}, ..., P_{last} is real-rooted with every \
                 root in {}",
                claimed.text()
            ),
            match failing {
                None => Ok(format!(
                    "checked exactly for P_{first}, ..., P_{last}{zero_note}"
                )),
                Some(offset) => Err((
                    Some(first + offset),
                    format!(
                        "P_{} is not real-rooted, or has a root outside {} (exact)",
                        first + offset,
                        claimed.text()
                    ),
                )),
            },
        ));
    }
    out
}

/// The multiplier hypotheses: `a_n = g h` real-rooted with every root in
/// the region for every `n >= N`.
fn multiplier_hypotheses(
    decomposition: &Decomposition,
    region: ClosedRegion,
    n_first: usize,
) -> Vec<LiftHypothesis> {
    let mut out = Vec::new();
    if !decomposition.factors.is_empty() {
        let denominators: BTreeSet<_> = decomposition.index_factor.iter().map(|(_, q)| q).collect();
        for (i, q) in denominators.into_iter().enumerate() {
            out.push(sign_hypothesis(
                &format!("multiplier_coefficients_defined_{i}"),
                "the common multiplier coefficients are defined for every rule index",
                "multiplier coefficient denominator",
                &IndexFunction::polynomial(q.clone()),
                n_first,
                IndexSignCondition::Nonzero,
            ));
        }
    }
    let g = &decomposition.fixed_factor;
    if vec_poly::degree(g).is_some_and(|d| d >= 1) {
        let text = format_index_polynomial(g, "x");
        out.push(exact_hypothesis(
            "multiplier_fixed_factor_in_region",
            &format!(
                "the multiplier factor g(x) = {text} is real-rooted with every root in {}",
                region.text()
            ),
            if region.contains_roots_of(g) {
                Ok("exact check of the fixed polynomial".into())
            } else {
                Err((
                    None,
                    format!(
                        "g(x) = {text} is not real-rooted with every root in {} (exact)",
                        region.text()
                    ),
                ))
            },
        ));
    }
    if decomposition.index_factor.len() >= 2 {
        let (prefix, reflect, nonpositive) = match region {
            ClosedRegion::NonpositiveRay => ("multiplier_index_factor", false, true),
            ClosedRegion::NonnegativeRay => ("multiplier_index_factor_reflected", true, true),
            ClosedRegion::RealLine => ("multiplier_index_factor", false, false),
        };
        let coefficients = decomposition.index_factor_functions(reflect);
        let (hypotheses, _) =
            polynomial_hypotheses(prefix, "x", &coefficients, n_first, nonpositive);
        out.extend(hypotheses);
    }
    out
}

/// The hypotheses of `S_n` for the region, with the theorem text and the
/// Lean declarations.
fn preserver_hypotheses(
    operator: &PreserverOperator,
    region: ClosedRegion,
    n_first: usize,
    previous_degree: Option<&[Q]>,
) -> (Vec<LiftHypothesis>, String, Vec<String>) {
    if let Some((parts, denominator)) = weyl_step_of(operator) {
        return weyl_step_preserver_hypotheses(&parts, &denominator, region, n_first);
    }
    let transfer = operator_transfer(operator, n_first);
    let mut out = transfer.hypotheses.clone();
    let mut theorem = transfer.theorem.clone();
    let mut lean_declarations = transfer.lean_declarations.clone();
    if let (PreserverOperator::EulerStep { raising, .. }, Some(previous_degree)) =
        (operator, previous_degree)
    {
        // The degree-local form: the step acts on P_(n-1) of degree
        // e(n) = deg P_(n-1), so u + v k > 0 is needed for k <= e(n) only
        // (linear in k: at k = 0, which is u > 0, and at k = e(n)); the
        // menu's v >= 0 is dropped.
        if let Some(hypothesis) = degree_local_raising_hypothesis(raising, previous_degree, n_first)
        {
            out.retain(|h| h.name != "euler_step_raising_slope_nonnegative");
            out.push(hypothesis);
            theorem = EULER_STEP_DEGREE_LOCAL.into();
            lean_declarations.push("RealRooted.EulerBidiagonal.generalRowsDep_spec".into());
        }
    }
    if let Some((_, denominator)) = operator.weyl_parts() {
        if vec_poly::degree(&denominator).is_some_and(|d| d >= 1) {
            out.push(sign_hypothesis(
                "preserver_coefficients_defined",
                &format!(
                    "the coefficients of S_n are finite at every n >= {n_first} (their common \
                     denominator does not vanish)"
                ),
                "delta_S(n)",
                &IndexFunction::polynomial(denominator),
                n_first,
                IndexSignCondition::Nonzero,
            ));
        }
    }
    let ray_image_in_ray = transfer
        .from_nonpositive
        .as_ref()
        .is_some_and(|image| region_contains(&nonpositive_ray(), image));
    let reflection_invariant = matches!(
        operator,
        PreserverOperator::EulerPolynomial { .. }
            | PreserverOperator::PowerTimesDerivative { .. }
            | PreserverOperator::DerivativeOfPowerMultiple { .. }
    );
    let (holds, how) = match region {
        ClosedRegion::NonpositiveRay => (
            ray_image_in_ray,
            "the menu image of (-inf, 0] lies in (-inf, 0]".to_string(),
        ),
        ClosedRegion::NonnegativeRay => (
            ray_image_in_ray && reflection_invariant,
            "S_n commutes up to sign with x -> -x (a polynomial in theta, x^m D^j or D^j x^m) \
             and the menu image of (-inf, 0] lies in (-inf, 0], so [0, inf) maps to [0, inf)"
                .to_string(),
        ),
        ClosedRegion::RealLine => (
            transfer.from_real_line,
            "the menu theorem holds for real-rooted polynomials anywhere on the line".to_string(),
        ),
    };
    out.push(exact_hypothesis(
        "preserver_keeps_region",
        &format!(
            "S_n maps every polynomial that is zero or real-rooted with every root in {} to zero or \
             such a polynomial, for every n >= {n_first}",
            region.text()
        ),
        if holds {
            Ok(how)
        } else {
            Err((None, format!("not given by the menu theorem: {how} fails")))
        },
    ));
    (out, theorem, lean_declarations)
}

/// The coefficients `(A, B, C)` and the denominator of a Weyl-step operator
/// `S_n = (A + B D + C D^2) / delta_S(n)`.
fn weyl_step_of(operator: &PreserverOperator) -> Option<([Bivariate; 3], Vec<Q>)> {
    if !matches!(operator, PreserverOperator::WeylOperator { .. }) {
        return None;
    }
    let (numerator, denominator) = operator.weyl_parts()?;
    Some((weyl_step_parts(&numerator)?, denominator))
}

/// The hypotheses of Brändén's Weyl step for the region (module
/// documentation, "Brändén's Weyl step").
fn weyl_step_preserver_hypotheses(
    parts: &[Bivariate; 3],
    denominator: &[Q],
    region: ClosedRegion,
    n_first: usize,
) -> (Vec<LiftHypothesis>, String, Vec<String>) {
    let mut out = weyl_step_hypotheses(parts, n_first);
    if vec_poly::degree(denominator).is_some_and(|d| d >= 1) {
        out.push(sign_hypothesis(
            "preserver_coefficients_defined",
            &format!(
                "the coefficients of S_n are finite at every n >= {n_first} (their common \
                 denominator does not vanish)"
            ),
            "delta_S(n)",
            &IndexFunction::polynomial(denominator.to_vec()),
            n_first,
            IndexSignCondition::Nonzero,
        ));
    }
    let how = match region {
        ClosedRegion::RealLine => {
            "Branden's Corollary 3.5 holds for real-rooted polynomials anywhere on the line"
        }
        ClosedRegion::NonpositiveRay => {
            out.extend(weyl_step_region_hypotheses(parts, false, n_first));
            "with coefficients of one sign, S_n maps a polynomial with coefficients of one sign \
             to one with coefficients of one sign, real-rooted by Corollary 3.5, so with every \
             root in (-inf, 0]"
        }
        ClosedRegion::NonnegativeRay => {
            out.extend(weyl_step_region_hypotheses(parts, true, n_first));
            "the reflected operator (x -> -x) has coefficients of one sign, so S_n keeps [0, inf) \
             by Corollary 3.5 and the reflection"
        }
    };
    out.push(exact_hypothesis(
        "preserver_keeps_region",
        &format!(
            "S_n maps every polynomial that is zero or real-rooted with every root in {} to zero or \
             such a polynomial, for every n >= {n_first}",
            region.text()
        ),
        Ok(how.to_string()),
    ));
    (out, BRANDEN_WEYL_STEP.to_string(), Vec::new())
}

const EULER_STEP_DEGREE_LOCAL: &str = "The degree-local Euler step: T = kappa (theta + a)(theta + \
     b) + x (u + v theta), theta = x D, with kappa, a, b > 0, u + v k > 0 for every k <= d and a \
     negative comparison defect kappa (a b - (a + b + 1)^2/4) + (u - (a + b + 1) v/2) s - v^2 \
     s^2/(4 kappa) < 0 for every s <= 0, maps every polynomial of degree <= d with nonnegative \
     coefficients and only real roots to such a polynomial (RealRooted \
     IsPFPolynomial.generalStep_of_natDegree_le, a half-line theorem: the symbol of T has \
     discriminant -4 kappa s^2 times the comparison defect at s, so its real-rootedness is \
     required for s <= 0 only, and the hypothesis (i) of P. Branden, Linear transformations \
     preserving the Polya frequency property, Trans. Amer. Math. Soc. 358 (2006), Corollary \
     3.5, for every real s, is not needed); here d = deg P_(n-1), so u(n) + v(n) k > 0 is \
     required for 0 <= k <= deg P_(n-1) only (both ends, the condition being linear in k); \
     by linearity a row with nonpositive coefficients is handled through its negative";

/// `u(n) + v(n) e(n) > 0` for every `n >= n_first`, with `e(n) = deg
/// P_(n-1)` (with `u > 0` from the menu, `u + v k > 0` for `0 <= k <= e(n)`).
fn degree_local_raising_hypothesis(
    raising: &[IndexRationalFunction],
    previous_degree: &[Q],
    n_first: usize,
) -> Option<LiftHypothesis> {
    let [u, v] = raising else {
        return None;
    };
    let ((pu, qu), (pv, qv)) = (u.parse().ok()?, v.parse().ok()?);
    let numerator = vec_poly::add(
        &vec_poly::mul(&pu, &qv),
        &vec_poly::mul(&vec_poly::mul(&pv, previous_degree), &qu),
    );
    let function = IndexFunction::new(numerator, vec_poly::mul(&qu, &qv));
    Some(sign_hypothesis(
        "euler_step_raising_positive_up_to_degree",
        &format!(
            "u(n) + v(n) deg P_(n-1) > 0 for every n >= {n_first} (with u(n) > 0: u + v k > 0 \
             for 0 <= k <= deg P_(n-1), the degree-local hypothesis)"
        ),
        "u(n) + v(n) deg P_(n-1)",
        &function,
        n_first,
        IndexSignCondition::Positive,
    ))
}

/// One decomposition with one region.
struct Attempt {
    /// The certified region of the rows `P_n`.
    claimed: ClosedRegion,
    decomposition: Decomposition,
    hypotheses: Vec<LiftHypothesis>,
    normal_form: String,
    reference: String,
    lean_declarations: Vec<String>,
}

impl Attempt {
    fn open_count(&self) -> usize {
        self.hypotheses
            .iter()
            .filter(|h| h.status != LiftHypothesisStatus::Verified)
            .count()
    }
}

/// The text of `S_n` in the normal form.
fn preserver_text(decomposition: &Decomposition) -> String {
    if !decomposition.factors.is_empty() {
        return format!("the exact product of {} first-order factors (application order recorded in the certificate)", decomposition.factors.len());
    }
    match &decomposition.preserver {
        None => "the identity".to_string(),
        Some(operator) => match weyl_step_of(operator) {
            Some((parts, denominator)) if vec_poly::degree(&denominator) == Some(0) => {
                let scale = vec_poly::trim(&denominator)[0].recip();
                let parts = parts.map(|p| {
                    p.iter()
                        .map(|row| vec_poly::scale(row, &scale))
                        .collect::<Bivariate>()
                });
                weyl_step_text(&parts)
            }
            Some((parts, denominator)) => format!(
                "({}) / ({})",
                weyl_step_text(&parts),
                format_index_polynomial(&denominator, "n")
            ),
            None => operator.text(),
        },
    }
}

fn attempt(family: &Family<'_>, decomposition: &Decomposition, region: ClosedRegion) -> Attempt {
    let claimed = family.coordinates.claimed_region(region);
    let mut hypotheses = Vec::new();
    hypotheses.extend(denominator_hypothesis(family.f, family.n_first));
    let operator_text = preserver_text(decomposition);
    hypotheses.push(exact_hypothesis(
        "operator_identity",
        "the step operator T_n of the rule equals a_n(x) S_n",
        Ok(format!(
            "q(n) delta_S(n) sum_d c_d(n, x) D^d = f(n) (g H)(n, x) o S_num holds exactly in the \
             Ore algebra (S_n = {operator_text})"
        )),
    ));
    hypotheses.extend(row_hypotheses(family, region, claimed));
    hypotheses.extend(family.degree_hypotheses.iter().cloned());
    hypotheses.extend(multiplier_hypotheses(decomposition, region, family.n_first));
    let mut reference = REFERENCE.to_string();
    let mut lean_declarations = Vec::new();
    if !decomposition.factors.is_empty() {
        reference = reference.replace(
            "maps every polynomial that is zero or real-rooted with",
            "maps every polynomial of the certified input degree that is zero or real-rooted with",
        );
        let (conditions, theorem) = factorization::factor_hypotheses(
            &decomposition.factors,
            region,
            family.n_first,
            family.previous_degree.as_deref(),
        );
        hypotheses.extend(conditions);
        reference.push_str(&theorem);
    }
    match &decomposition.preserver {
        None if decomposition.factors.is_empty() => {
            reference.push_str("; here S_n is the identity")
        }
        None => {}
        Some(operator) => {
            let (side_conditions, theorem, declarations) = preserver_hypotheses(
                operator,
                region,
                family.n_first,
                family.previous_degree.as_deref(),
            );
            hypotheses.extend(side_conditions);
            reference.push_str(&format!("; S_n = {operator_text}: {theorem}"));
            lean_declarations = declarations;
        }
    }
    if region == ClosedRegion::NonnegativeRay {
        reference.push_str(
            "; for R = [0, inf) the reflection x -> -x maps the rows to rows with every root in \
             (-inf, 0]",
        );
    }
    if family.lag > 1 {
        reference.push_str(&format!(
            "; with lag L = {}, the induction runs in each residue class of n modulo L from its \
             start row",
            family.lag
        ));
    }
    if let Some(text) = family.coordinates.text() {
        reference.push_str(&format!(
            "; the induction runs on the rows {text} (exact coordinate change of the rule; the \
             rows P_n have every root in {} when the rows Q_n have every root in {})",
            claimed.text(),
            region.text()
        ));
    }
    if let Some(requested) = family.requested {
        let statement = format!(
            "the requested region {} contains the certified region {}",
            open_root_interval(requested),
            claimed.text()
        );
        hypotheses.push(if claimed.lies_in(requested) {
            exact_hypothesis("requested_root_region", &statement, Ok("exact".into()))
        } else {
            unresolved_hypothesis(
                "requested_root_region",
                &statement,
                "this check certifies the closed regions (-inf, 0], [0, inf) and the line only"
                    .into(),
            )
        });
    }
    let letter = family.row_letter();
    let lag = family.lag;
    let mut normal_form = format!(
        "{letter}_n = a_n(x) S_n({letter}_(n-{lag})) for n >= {}, with a_n(x) = {} and S_n = \
         {operator_text}",
        family.n_first,
        decomposition.multiplier_text()
    );
    if let Some(text) = family.coordinates.text() {
        normal_form = format!("{text}: {normal_form}");
    }
    Attempt {
        claimed,
        decomposition: decomposition.clone(),
        hypotheses,
        normal_form,
        reference,
        lean_declarations,
    }
}

/// [`TheoremRoute::is_preferred`]: a proof in coordinates other than the
/// identity is replaced by a later start that proves the theorem in the
/// identity coordinates.
fn proof_in_identity_coordinates(report: &RealRootednessTheoremReport) -> bool {
    !report
        .iterated_preserver_certificate()
        .is_some_and(|certificate| !certificate.coordinates.is_identity())
}

/// The registry entry of the iterated-preserver theorem ([`super::registry::ROUTES`]).
pub(super) const ROUTE: TheoremRoute = TheoremRoute {
    theorem: RealRootednessTheorem::IteratedPreserver,
    check_at_supplied_seeds: iterated_preserver_report,
    start_policy: StartPolicy::MovesStart,
    is_preferred: proof_in_identity_coordinates,
    certificate_name: Some("iterated-preserver"),
    discovery: DiscoveryPolicy::Always,
    in_residue_class_menu: true,
    in_index_section_menu: true,
    parent: ParentPolicy {
        stage: ParentStage::ReplacesEarly,
        region: ParentRegion::Unrestricted,
        record_only_if_applicable: true,
    },
};

/// Replay the certificate of a recomputed iterated-preserver report: the
/// operator identity of the step on the checked rule
/// ([`IteratedPreserverCertificate::verify_operator_identity`]).
pub(super) fn verify_certificate_in_report(
    certificate: &IteratedPreserverCertificate,
    report: &RealRootednessTheoremReport,
) -> Result<(), RecurrenceError> {
    let checked_rule = match &report.left_hand_factor_normalization {
        Some(normalization) => normalization.normalized_recurrence()?,
        None => report.recurrence.to_recurrence_parts()?.0,
    };
    certificate.verify_operator_identity(&checked_rule)
}

/// The check with the induction at the given (normalized) seeds.
fn iterated_preserver_report(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<BigRational>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    let report = RealRootednessTheoremReport::new(
        RealRootednessTheorem::IteratedPreserver,
        REFERENCE,
        rule,
        first_index,
        seeds,
        options,
    );
    if seeds.is_empty() {
        return Ok(report.not_applicable("the theorem needs a seed row"));
    }
    let (terms, f, lag) = match step_operator(rule) {
        Ok(step) => step,
        Err(reason) => return Ok(report.not_applicable(reason)),
    };
    if seeds.len() < lag {
        return Ok(
            report.not_applicable(format!("the rule has lag {lag} and needs {lag} seed rows"))
        );
    }
    // The requested region as the report records it: open at its finite ends
    // (so that a replay sees the same request).
    let requested = options.root_region.as_ref().map(|region| RealInterval {
        lower: region.lower.clone(),
        lower_closed: false,
        upper: region.upper.clone(),
        upper_closed: false,
    });
    let n_first = first_index + seeds.len();
    let mut coordinate_list = vec![IteratedPreserverCoordinates::Identity];
    coordinate_list.extend(coordinates::coordinate_candidates(
        &terms,
        lag,
        n_first,
        &seeds[seeds.len() - lag..],
    ));
    let mut best: Option<(Attempt, Family<'_>)> = None;
    'coordinates: for coordinates in coordinate_list {
        let Some(family) = Family::new(
            &terms,
            &f,
            lag,
            coordinates,
            first_index,
            seeds,
            requested.as_ref(),
        ) else {
            continue;
        };
        let candidates = decompositions(&family.terms, &f);
        let mut found: Option<Attempt> = None;
        let mut proved = false;
        'regions: for region in ClosedRegion::ALL {
            for decomposition in &candidates {
                let attempt = attempt(&family, decomposition, region);
                proved = status_of(&attempt.hypotheses) == RealRootednessTheoremStatus::Proved;
                let better = found.as_ref().map_or_else(
                    || {
                        best.as_ref()
                            .is_none_or(|(b, _)| attempt.open_count() < b.open_count())
                    },
                    |b| attempt.open_count() < b.open_count(),
                );
                if proved || better {
                    found = Some(attempt);
                }
                if proved {
                    break 'regions;
                }
            }
        }
        if let Some(attempt) = found {
            best = Some((attempt, family));
        }
        if proved {
            break 'coordinates;
        }
    }
    if !best.as_ref().is_some_and(|(a, _)| a.open_count() == 0) && lag == 1 {
        if let Some(family) = Family::new(
            &terms,
            &f,
            lag,
            IteratedPreserverCoordinates::Identity,
            first_index,
            seeds,
            requested.as_ref(),
        ) {
            let candidates = factorization::decompositions(&family.terms, &f);
            'factored: for region in ClosedRegion::ALL {
                for candidate in &candidates {
                    let factored = attempt(&family, candidate, region);
                    if factored.open_count() == 0 {
                        best = Some((factored, family));
                        break 'factored;
                    }
                }
            }
        }
    }
    let Some((attempt, family)) = best else {
        return Ok(report.not_applicable(
            "the step T_n is not a multiplier a_n(x) after the identity or an operator of the \
             preserver menu (phi(n, D), f(n, theta), the Euler step, x^m D^j, D^j x^m)",
        ));
    };
    Ok(finish_report(report, &family, attempt))
}

/// Brändén's Corollary 2.5 for a step `T = g(x) S` free of `n` with `S` a
/// proved Weyl step (module documentation, "Interlacing"): the first index
/// `N` when `P_(N-1) ⪯ P_N` is checked exactly.
fn weyl_step_interlacing_start(family: &Family<'_>, attempt: &Attempt) -> Option<usize> {
    let applies = family.lag == 1
        && family.coordinates.is_identity()
        && family.delta == Some(1)
        && family.f.len() <= 1
        && family
            .original_terms
            .values()
            .all(|c| c.iter().skip(1).all(|row| row.iter().all(Zero::is_zero)))
        && attempt
            .decomposition
            .preserver
            .as_ref()
            .is_some_and(|operator| weyl_step_of(operator).is_some());
    if !applies {
        return None;
    }
    let previous = vec_poly::trim(&family.seeds[family.seeds.len() - 1]);
    let f0 = family.f.first()?.clone();
    let next = vec_poly::scale(
        &ore_operator(family.original_terms).apply_row(&previous, family.n_first),
        &f0.recip(),
    );
    let (previous, next) = (
        integer_coefficients(&previous),
        integer_coefficients(&vec_poly::trim(&next)),
    );
    if vec_poly::is_zero(&previous) || vec_poly::is_zero(&next) {
        return None;
    }
    (crate::check_weak_interlacing_bigint_coeffs(&previous, &next) == Some(true))
        .then_some(family.n_first)
}

fn finish_report(
    mut report: RealRootednessTheoremReport,
    family: &Family<'_>,
    attempt: Attempt,
) -> RealRootednessTheoremReport {
    let region = attempt.claimed;
    let delta = family.delta.unwrap_or(0).max(0) as usize;
    let proved = status_of(&attempt.hypotheses) == RealRootednessTheoremStatus::Proved;
    let pure_linear = attempt.decomposition.preserver.is_none()
        && attempt.decomposition.factors.is_empty()
        && delta == 1
        && family.lag == 1
        && family.coordinates.is_identity()
        && family.terms.keys().all(|&d| d == 0);
    let mut weyl_interlacing = false;
    let interlacing_from = if pure_linear {
        Some(family.n_first)
    } else if proved {
        let start = weyl_step_interlacing_start(family, &attempt);
        weyl_interlacing = start.is_some();
        start
    } else {
        None
    }
    .map(|from| {
        extend_interlacing_backward(family.seeds, family.first_index, from, RowInterlacing::Weak)
    });
    report.reference = attempt.reference;
    if weyl_interlacing {
        report.reference.push_str(
            "; a linear operator preserving real-rootedness preserves the alternation of two \
             polynomials (P. Branden, Trans. Amer. Math. Soc. 358 (2006), Corollary 2.5)",
        );
    }
    report.normal_form = attempt.normal_form;
    report.hypotheses = attempt.hypotheses;
    report.lean_declarations = attempt.lean_declarations;
    let decomposition = &attempt.decomposition;
    // A step that lowers the degree has no certificate (its degree step is
    // not a valid record, and the check is not proved).
    let lowers_degree = family.delta.is_none_or(|d| d < 0);
    let coordinates_note = family.coordinates.text().map_or(String::new(), |text| {
        format!(" (in the coordinates {text})")
    });
    report.certificate = (!lowers_degree).then(|| {
        RouteCertificate::IteratedPreserver(IteratedPreserverCertificate {
            induction_start: family.start,
            start_degree: family.start_degree,
            degree_step: delta,
            multiplier_fixed_factor: decomposition
                .fixed_factor
                .iter()
                .map(format_rational_coeff)
                .collect(),
            multiplier_index_factor: decomposition.index_factor_functions(false),
            preserver: decomposition.preserver.clone(),
            factors: decomposition.factors.clone(),
            region: region.root_interval(),
            operator_identity: format!(
                "f(n) T_n = sum_d c_d(n, x) D^d equals f(n) a_n(x) S_n with a_n(x) = {}{coordinates_note}",
                decomposition.multiplier_text()
            ),
            lag: family.lag,
            coordinates: family.coordinates.clone(),
        })
    });
    let root_region = (region != ClosedRegion::RealLine).then(|| region.root_interval());
    let claim_first = claim_first_index(family.first_index, family.seeds, family.lag);
    let letter = family.row_letter();
    let lag = family.lag;
    let mut statement = format!(
        "every P_n, n >= {claim_first}, is nonzero and real-rooted with every root in {}; deg \
         {letter}_n = deg {letter}_(n-{lag}) + {delta} for n >= {}",
        region.text(),
        family.n_first
    );
    if let Some(text) = family.coordinates.text() {
        statement.push_str(&format!(" ({text})"));
    }
    if claim_first > family.first_index {
        statement.push_str(&format!(
            "; P_{}, ..., P_{} are zero or real-rooted with every root in {} (seed rows, \
             checked exactly)",
            family.first_index,
            claim_first - 1,
            region.text()
        ));
    }
    if let Some(from) = interlacing_from {
        if weyl_interlacing {
            statement.push_str(&format!(
                "; P_(n-1) interlaces P_n (weakly: the step preserves the alternation of \
                 P_(N-1) and P_N, Branden's Corollary 2.5) for n >= {from}"
            ));
        } else {
            statement.push_str(&format!(
                "; P_(n-1) interlaces P_n (weakly: P_n = a_n P_(n-1) with a_n linear adds one \
                 real root) for n >= {from}"
            ));
        }
    }
    report.finish(RealRootednessTheoremConclusion {
        first_index: claim_first,
        zero_residue_classes: None,
        interlacing_from,
        interlacing: interlacing_from.map(|_| RowInterlacing::Weak),
        simple_roots: false,
        root_region,
        statement,
    })
}
