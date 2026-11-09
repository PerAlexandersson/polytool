//! The Ma--Wang interval step: first-order derivative recurrences
//! `P_n = (alpha_n x + beta_n) P_(n-1) + kappa_n (x - L)(x - U) P'_(n-1)`
//! with fixed rational `L < U`, whose rows stay real-rooted with every root
//! in the closed interval `[L, U]`, by the exact two-endpoint criterion with
//! endpoint multiplicities tracked through the induction.
//!
//! # The step
//!
//! S.-M. Ma and Y. Wang, *q-Eulerian polynomials and polynomials with only
//! real zeros*, Electron. J. Combin. 15 (2008) R17, Theorem 2 and
//! Corollary 3 (arXiv:math/0611826), in the closed two-endpoint form of the
//! public Lean library RealRooted
//! (<https://github.com/PerAlexandersson/RealRooted>),
//! `RealRooted.affineEndpointStep_closed`: let `p` be real-rooted with a
//! positive leading coefficient and every root in `[-1, 0]`, let `a > 0`
//! and `c` be real, and
//!
//! ```text
//! F = (a y + c) p + y (y + 1) p'      (RealRooted.MaWang.endpointStep 1 a c p).
//! ```
//!
//! Then `deg F = deg p + 1`, `F` has a positive leading coefficient, is
//! real-rooted, `p` interlaces `F` (`Interlaces p F`: the roots alternate
//! weakly, `F` owning the outer ones), and
//!
//! - every root of `F` lies in `[-1, 0]` **iff** `c + m_0 >= 0` and
//!   `a + m_(-1) - c >= 0`, where `m_0`, `m_(-1)` are the multiplicities of
//!   `0` and `-1` as roots of `p`
//!   (`RealRooted.affineEndpointStep_roots_mem_Icc_iff`);
//! - the multiplicity of `0` in `F` is `m_0 + 1` when `c + m_0 = 0` and
//!   `m_0` otherwise; that of `-1` is `m_(-1) + 1` when `c = a + m_(-1)` and
//!   `m_(-1)` otherwise (`RealRooted.affineEndpointStep_rootMultiplicity_zero`,
//!   `RealRooted.affineEndpointStep_rootMultiplicity_neg_one`).
//!
//! # The general interval
//!
//! For the rule `f(n) P_n = c_0(n, x) P_(n-1) + c_1(n, x) P'_(n-1)` with
//! `deg_x c_0 <= 1` and `c_1(n, x) = gamma(n) (x - L)(x - U)`, `L < U`
//! rational (the roots of the `x`-factor of `c_1`, which must be separable
//! with a monic quadratic of distinct rational roots), put `h = U - L`,
//! `A_n(x) = c_0(n, x) / f(n)` and `kappa_n = gamma(n) / f(n)`.  The change of
//! variable `x = U + h y` maps `[L, U]` onto `[-1, 0]`, and with
//! `g(y) = p(U + h y)`, since `(x - L)(x - U) = h^2 y (y + 1)` and
//! `g'(y) = h p'(x)`,
//!
//! ```text
//! A_n(x) p + kappa_n (x - L)(x - U) p' = kappa_n h [(a_n y + c_n) g + y (y + 1) g'],
//! a_n = [x^1] c_0(n, x) / gamma(n),   c_n = c_0(n, U) / (h gamma(n)).
//! ```
//!
//! The factor `kappa_n h` is a nonzero constant (it is nonzero wherever
//! `a_n` is defined and positive), so it moves no root; and the step is
//! linear, so a negative leading coefficient of `p` is handled by applying
//! the theorem to `-p`.  Multiplicities at `U` and `L` are those at `0` and
//! `-1` after the change of variable.  The two thresholds become
//!
//! ```text
//! E_U(n) = e_U(n) + m_U(n) >= 0,   e_U(n) = c_0(n, U) / (h gamma(n)),
//! E_L(n) = e_L(n) + m_L(n) >= 0,   e_L(n) = -c_0(n, L) / (h gamma(n)),
//! ```
//!
//! (`c_n + m_0` and `a_n + m_(-1) - c_n`, as `h [x^1] c_0 - c_0(n, U) =
//! -c_0(n, L)`), which is the criterion `alpha U + beta + kappa (U - L) m_U
//! >= 0`, `-alpha L - beta + kappa (U - L) m_L >= 0` divided by
//! `kappa (U - L)`.
//!
//! # The induction
//!
//! Let `N = first_index + seeds.len()` and `m_U(N)`, `m_L(N)` the endpoint
//! multiplicities of the last seed `P_(N-1)`.  Assume, for every `n >= N`:
//!
//! 1. `f(n) != 0`;
//! 2. `P_(N-1)` is nonzero and real-rooted with every root in `[L, U]`
//!    (and every earlier seed row too);
//! 3. `a_n > 0`;
//! 4. `E_U(n) >= 0` and `E_L(n) >= 0`, where `m(n + 1) = m(n) + 1` if
//!    `E(n) = 0` and `m(n + 1) = m(n)` otherwise, at each end.
//!
//! Then by induction on `n` (the step above at every `n`), every `P_n`,
//! `n >= N - 1`, is nonzero and real-rooted with every root in `[L, U]`,
//! `deg P_n = deg P_(n-1) + 1` and `P_(n-1)` interlaces `P_n`, and `m(n + 1)`
//! is the multiplicity of the end in `P_n`.  Conversely, if all hypotheses
//! hold up to `n - 1` and `E(n) < 0` at one end, then `P_n` has a root
//! outside `[L, U]` (the "iff"), so the interval claim is refuted there.
//!
//! # Deciding the thresholds for every `n`
//!
//! The counter makes condition 4 a statement about the integer sequence
//! `m(n)`, which is decided exactly as follows (with the current state
//! `(n, m)`, starting at `(N, m(N))`):
//!
//! - if `e(n') + m >= 0` for every `n' >= n` (the sign checker,
//!   [`crate::decide_index_rational_function_sign`]), the condition holds
//!   from `n` on, because `m(n')` never decreases;
//! - otherwise the sign checker returns the smallest `i >= n` with
//!   `e(i) + m < 0`.  If `e + m` has no zero at an integer in `[n, i)`, the
//!   counter stays `m` up to `i` and `E(i) < 0`: the condition fails at `i`
//!   (exact).  Otherwise the first such zero `j` is a hit: the state becomes
//!   `(j + 1, m + 1)`, and the test repeats;
//! - **lattice case**: if `e` is affine in `n`, `e(n) = e_0 - d n` with
//!   `d > 0`, `1/d` is a positive integer `r` and `F = e(n) + m` is a
//!   nonnegative integer multiple of `d`, then `E` stays in
//!   `d Z_(>= 0)` forever: it drops by `d` per step and, at `0`, the counter
//!   raises it to `1 - d = (r - 1) d`.  So the condition holds for every
//!   later `n`, with a hit every `r` steps (the ceiling-division
//!   multiplicity of `RealRooted.MaWang.affineEndpointSequence`, and the
//!   slice `d = 1/r` of the affine-parameter classification
//!   `RealRooted.normalizedAffineEndpointSequence_roots_mem_Icc_iff`).
//!
//! At most 64 hits are followed, and a zero is searched in at most 100000
//! consecutive indices; beyond that
//! the condition is unresolved (never failed).  An end is reported closed in
//! the certified region when some row may have a root there: a seed row with
//! that root, a hit, the lattice case, or `e + m` not proved positive after
//! the last hit; otherwise the end is open.
//!
//! # Lean correspondence
//!
//! [`RealRootednessTheoremReport::lean_declarations`] lists the step theorems
//! above; no RealRooted row route mirrors this check
//! ([`RealRootednessTheoremReport::lean_route`] is omitted), and a Lean replay
//! must establish the change of variable and the counter induction.  The
//! Lean step assumes `a > 0`; the weaker degree condition
//! `alpha_n + kappa_n deg P_(n-1) > 0` is not used here.

use super::*;
use num_traits::ToPrimitive;

const REFERENCE: &str = "S.-M. Ma and Y. Wang, q-Eulerian polynomials and polynomials with only \
     real zeros, Electron. J. Combin. 15 (2008) R17, Theorem 2 and Corollary 3 \
     (arXiv:math/0611826), in the closed two-endpoint form of the RealRooted library \
     (https://github.com/PerAlexandersson/RealRooted, RealRooted.affineEndpointStep_closed), \
     after the change of variable x = U + (U - L) y, with the endpoint multiplicities tracked \
     exactly through the induction";

/// Largest number of endpoint hits (indices where a threshold vanishes)
/// followed before a threshold is reported unresolved.
const MAXIMUM_ENDPOINT_HITS: usize = 64;

/// Largest number of consecutive indices searched for the next endpoint hit
/// before a threshold is reported unresolved.
const MAXIMUM_HIT_SCAN: usize = 100_000;

/// RealRooted declarations behind the check (the theorem a replay applies
/// first, then its supporting lemmas).
const LEAN_DECLARATIONS: [&str; 7] = [
    "RealRooted.affineEndpointStep_closed",
    "RealRooted.affineEndpointStep_roots_mem_Icc_iff",
    "RealRooted.affineEndpointStep_rootMultiplicity_zero",
    "RealRooted.affineEndpointStep_rootMultiplicity_neg_one",
    "RealRooted.MaWang.endpointStep",
    "RealRooted.MaWang.affineEndpointSequence",
    "RealRooted.normalizedAffineEndpointSequence_roots_mem_Icc_iff",
];

/// Check the Ma--Wang interval theorem (module documentation of the
/// `ma_wang` submodule, summarized here) for the family given by `rule`,
/// `first_index` and the seed rows `seeds` (the rule applies from
/// `N = first_index + seeds.len()` on).
///
/// The rule must be `f(n) P_n = c_0(n, x) P_(n-1) + c_1(n, x) P'_(n-1)` with
/// `f` free of `x`, `deg_x c_0 <= 1` and `c_1(n, x) = gamma(n) (x - L)(x - U)`
/// for rationals `L < U` (read off the rule).  With `h = U - L` and the
/// endpoint multiplicities `m_L`, `m_U` of the current row, the hypotheses,
/// each decided exactly for every `n >= N`, are:
///
/// - `f(n) != 0` (sign checker);
/// - `P_(N-1)` and every earlier seed row nonzero and real-rooted with every
///   root in `[L, U]` (exact; this gives `m_L(N)`, `m_U(N)`);
/// - `a_n = [x^1] c_0(n, x) / gamma(n) > 0` (sign checker);
/// - `c_0(n, U) / (h gamma(n)) + m_U(n) >= 0` and
///   `-c_0(n, L) / (h gamma(n)) + m_L(n) >= 0`, each multiplicity raised by
///   one after an index where its threshold vanishes (sign checker on each
///   stretch between such indices, or the lattice argument);
/// - a requested root region must contain the certified region.
///
/// When all hold, every `P_n`, `n >= first_index`, is nonzero and
/// real-rooted with every root in `[L, U]` (an end open when no row has a
/// root there), `deg P_n = deg P_(n-1) + 1` and `P_(n-1)` interlaces `P_n`
/// (weakly) for `n >= N`.  A failing threshold at an index `i` means that
/// `P_i` has a root outside `[L, U]` (the criterion is an equivalence).
///
/// # Errors
///
/// None at present; the `Result` is kept for symmetry with the other
/// theorem checks.
///
/// # Example
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::recurrence::{
///     check_ma_wang_interval_theorem, BigRational, BivarPoly, Recurrence, RecurrenceTerm,
///     RealRootednessTheoremOptions,
/// };
///
/// let q = |v: i64| BigRational::from_integer(BigInt::from(v));
/// // Ordered Bell (Fubini) polynomials: F_n = x F_(n-1) + x (1 + x) F'_(n-1), F_0 = 1.
/// let rule = Recurrence {
///     terms: vec![
///         RecurrenceTerm::new(1, 0, BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(1)]])),
///         RecurrenceTerm::new(
///             1,
///             1,
///             BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(1), q(1)]]),
///         ),
///     ],
///     denominator: None,
///     inhomogeneous: None,
/// };
/// let report = check_ma_wang_interval_theorem(
///     &rule,
///     0,
///     &[vec![q(1)]],
///     &RealRootednessTheoremOptions::default(),
/// )
/// .unwrap();
/// assert!(report.is_proved());
/// let region = report.conclusion.as_ref().unwrap().root_region.as_ref().unwrap();
/// // F_n(-1) = (-1)^n, so -1 is never a root; 0 is a root from F_1 = x on.
/// assert_eq!(region.to_string(), "(-1, 0]");
/// report.verify().unwrap();
/// ```
pub fn check_ma_wang_interval_theorem(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<BigRational>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    check_real_rootedness_theorem(
        RealRootednessTheorem::MaWangInterval,
        rule,
        first_index,
        seeds,
        options,
    )
}

/// The rational square root of `value >= 0`, when it exists.
fn rational_square_root(value: &Q) -> Option<Q> {
    if value.is_negative() {
        return None;
    }
    let numerator = value.numer().sqrt();
    let denominator = value.denom().sqrt();
    let root = Q::new(numerator, denominator);
    (&root * &root == *value).then_some(root)
}

/// `c(n, point)` as a polynomial in `n`.
fn evaluate_in_x(c: &Bivariate, point: &Q) -> Vec<Q> {
    let mut out = vec![Q::zero()];
    let mut power = Q::one();
    for j in 0..x_width(c) {
        out = vec_poly::add(&out, &vec_poly::scale(&x_column(c, j), &power));
        power *= point;
    }
    vec_poly::trim(&out)
}

/// `p / (x - point)` when `point` is a root of the nonzero `p`.
fn divide_by_root(p: &[Q], point: &Q) -> Option<Vec<Q>> {
    let degree = vec_poly::degree(p)?;
    if degree == 0 || !vec_poly::evaluate(p, point).is_zero() {
        return None;
    }
    // Synthetic division, highest coefficient first.
    let mut quotient = vec![Q::zero(); degree];
    let mut carry = Q::zero();
    for k in (1..=degree).rev() {
        carry = &p[k] + &carry * point;
        quotient[k - 1] = carry.clone();
    }
    Some(quotient)
}

/// The multiplicities `(m_L, m_U)` of `L` and `U` as roots of `p` when `p`
/// is nonzero and real-rooted with every root in `[L, U]`; `None` otherwise.
fn closed_interval_multiplicities(p: &[Q], lower: &Q, upper: &Q) -> Option<(usize, usize)> {
    let mut rest = vec_poly::trim(p);
    if vec_poly::is_zero(&rest) {
        return None;
    }
    let mut multiplicities = [0usize; 2];
    for (slot, point) in [lower, upper].into_iter().enumerate() {
        while let Some(quotient) = divide_by_root(&rest, point) {
            rest = quotient;
            multiplicities[slot] += 1;
        }
    }
    let open = RealInterval {
        lower: Some(lower.clone()),
        lower_closed: false,
        upper: Some(upper.clone()),
        upper_closed: false,
    };
    has_real_roots_in(&rest, &open).then_some((multiplicities[0], multiplicities[1]))
}

/// How the threshold at one end was decided.
struct EndpointDecision {
    hypothesis: LiftHypothesis,
    /// Whether some row from the induction start on may have a root at the
    /// end (seed rows are handled by the caller).
    attained: bool,
}

/// The lattice argument (module documentation): `e` affine with slope
/// `-d < 0`, `1/d` a positive integer and `e(n) + m` a nonnegative multiple
/// of `d`.  Returns the evidence when it applies.
fn lattice_argument(e: &IndexFunction, n: usize, m: usize) -> Option<String> {
    if e.denominator.len() != 1 || e.numerator.len() != 2 {
        return None;
    }
    let slope = &e.numerator[1];
    if !slope.is_negative() {
        return None;
    }
    let d = -slope;
    let r = d.recip();
    if !r.is_integer() {
        return None;
    }
    let value = e.at(n)? + Q::from_integer(BigInt::from(m));
    let ratio = &value / &d;
    if value.is_negative() || !ratio.is_integer() {
        return None;
    }
    Some(format!(
        "lattice argument from n = {n} (multiplicity {m}): e(n) = {} is affine with slope -{} \
         and 1/{} = {} is an integer, and E({n}) = {} is a nonnegative multiple of {}; E drops \
         by {} per step and returns from 0 to {} after each hit, so E(n') >= 0 for every \
         n' >= {n}, with a hit every {} steps",
        e.text(),
        format_rational_coeff(&d),
        format_rational_coeff(&d),
        format_rational_coeff(&r),
        format_rational_coeff(&value),
        format_rational_coeff(&d),
        format_rational_coeff(&d),
        format_rational_coeff(&(Q::one() - &d)),
        format_rational_coeff(&r),
    ))
}

/// Decide `E(n) = e(n) + m(n) >= 0` for every `n >= first` with the counter
/// `m(first) = initial` raised after every zero (module documentation,
/// "Deciding the thresholds for every `n`").
fn endpoint_threshold(
    name: &str,
    end: &str,
    e: &IndexFunction,
    first: usize,
    initial: usize,
) -> EndpointDecision {
    let statement = format!(
        "E(n) = e(n) + m(n) >= 0 for every n >= {first}, with e(n) = {} at the end {end}, \
         m({first}) = {initial} the multiplicity of {end} in the last seed, and m raised by one \
         after every n with E(n) = 0",
        e.text()
    );
    let mut n = first;
    let mut m = initial;
    let mut hits: Vec<usize> = Vec::new();
    let hits_text = |hits: &[usize]| {
        if hits.is_empty() {
            "no hit".to_string()
        } else {
            format!(
                "E(n) = 0 (multiplicity raised) at n = {}",
                hits.iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(", ")
            )
        }
    };
    while hits.len() <= MAXIMUM_ENDPOINT_HITS {
        let shifted = e.add(&IndexFunction::constant(Q::from_integer(BigInt::from(m))));
        let decision = shifted.decide(n, IndexSignCondition::Nonnegative);
        match decision {
            IndexSignDecision::Holds { .. } => {
                let positive = shifted.decide(n, IndexSignCondition::Positive).holds();
                let mut hypothesis = LiftHypothesis::from_sign_decision(
                    name,
                    &statement,
                    decision,
                    &format!("e(n) + {m} = {}", shifted.text()),
                );
                hypothesis.evidence = format!(
                    "{} before n = {n}; from n = {n} on (multiplicity {m}): {}",
                    hits_text(&hits),
                    hypothesis.evidence
                );
                return EndpointDecision {
                    hypothesis,
                    attained: !hits.is_empty() || !positive,
                };
            }
            IndexSignDecision::Undecided { reason } => {
                return EndpointDecision {
                    hypothesis: unresolved_hypothesis(
                        name,
                        &statement,
                        format!(
                            "{} before n = {n}; the sign checker did not decide e(n) + {m} >= 0 \
                             from n = {n} on: {reason}",
                            hits_text(&hits)
                        ),
                    ),
                    attained: true,
                };
            }
            IndexSignDecision::FailsAt { index, .. } => {
                if let Some(argument) = lattice_argument(e, n, m) {
                    return EndpointDecision {
                        hypothesis: exact_hypothesis(
                            name,
                            &statement,
                            Ok(format!("{} before n = {n}; {argument}", hits_text(&hits))),
                        ),
                        attained: true,
                    };
                }
                let Some(failing) = index.to_usize().filter(|i| i - n <= MAXIMUM_HIT_SCAN) else {
                    return EndpointDecision {
                        hypothesis: unresolved_hypothesis(
                            name,
                            &statement,
                            format!(
                                "{} before n = {n}; e(n) + {m} first becomes negative at \
                                 n = {index}, more than {MAXIMUM_HIT_SCAN} indices later; the \
                                 hits in between are not searched",
                                hits_text(&hits)
                            ),
                        ),
                        attained: true,
                    };
                };
                let zero = (n..failing).find(|&j| shifted.at(j).is_some_and(|v| v.is_zero()));
                match zero {
                    Some(j) => {
                        hits.push(j);
                        m += 1;
                        n = j + 1;
                    }
                    None => {
                        let value = shifted
                            .at(failing)
                            .map_or("undefined".to_string(), |v| format_rational_coeff(&v));
                        return EndpointDecision {
                            hypothesis: exact_hypothesis(
                                name,
                                &statement,
                                Err((
                                    Some(failing),
                                    format!(
                                        "{} before n = {n}; E({failing}) = e({failing}) + {m} = \
                                         {value} < 0 with no hit in between (exact); when \
                                         a_n > 0 for every n <= {failing}, the criterion puts a \
                                         root of P_{failing} outside the interval",
                                        hits_text(&hits)
                                    ),
                                )),
                            ),
                            attained: true,
                        };
                    }
                }
            }
        }
    }
    EndpointDecision {
        hypothesis: unresolved_hypothesis(
            name,
            &statement,
            format!(
                "{}; more than {MAXIMUM_ENDPOINT_HITS} hits, the induction is not followed \
                 further",
                hits_text(&hits)
            ),
        ),
        attained: true,
    }
}

/// The registry entry of the Ma--Wang interval theorem ([`super::registry::ROUTES`]).
pub(super) const ROUTE: TheoremRoute = TheoremRoute {
    theorem: RealRootednessTheorem::MaWangInterval,
    check_at_supplied_seeds: ma_wang_interval_report,
};

/// The check with the induction at the given (normalized) seeds.
fn ma_wang_interval_report(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<BigRational>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    let mut report = RealRootednessTheoremReport::new(
        RealRootednessTheorem::MaWangInterval,
        REFERENCE,
        rule,
        first_index,
        seeds,
        options,
    );
    if seeds.is_empty() {
        return Ok(report.not_applicable("the theorem needs a seed row"));
    }
    let shape = match recurrence_shape(rule) {
        Ok(shape) => shape,
        Err(reason) => return Ok(report.not_applicable(reason)),
    };
    let form = "the rule is not f(n) P_n = c_0(n, x) P_(n-1) + gamma(n) (x - L)(x - U) P'_(n-1)";
    if !shape.terms.contains_key(&(1, 1))
        || shape
            .terms
            .keys()
            .any(|key| *key != (1, 0) && *key != (1, 1))
    {
        return Ok(report.not_applicable(form));
    }
    let empty = Bivariate::new();
    let c0 = shape.terms.get(&(1, 0)).unwrap_or(&empty);
    let c1 = &shape.terms[&(1, 1)];
    if x_degree(c0).is_some_and(|d| d > 1) || x_degree(c1) != Some(2) {
        return Ok(report.not_applicable(format!("{form}: deg_x c_0 > 1 or deg_x c_1 != 2")));
    }
    let Some((gamma, w)) = separate(c1) else {
        return Ok(report.not_applicable(format!("{form}: c_1 is not gamma(n) w(x)")));
    };
    let leading = w[2].clone();
    let gamma = vec_poly::scale(&gamma, &leading);
    let (p, q) = (&w[1] / &leading, &w[0] / &leading);
    let discriminant = &p * &p - Q::from_integer(BigInt::from(4)) * &q;
    let root = match rational_square_root(&discriminant) {
        Some(root) if root.is_positive() => root,
        _ => {
            return Ok(report.not_applicable(format!(
                "{form}: the factor {} of c_1 has no two distinct rational roots",
                format_index_polynomial(&monic(&w), "x")
            )))
        }
    };
    let two = Q::from_integer(BigInt::from(2));
    let lower = (-&p - &root) / &two;
    let upper = (-&p + &root) / &two;
    let width = &upper - &lower;
    let n_first = first_index + seeds.len();
    let start = n_first - 1;
    let f = &shape.denominator;
    let over_f = |numerator: Vec<Q>| IndexFunction::new(numerator, f.clone());
    let (lower_text, upper_text) = (format_rational_coeff(&lower), format_rational_coeff(&upper));
    report.normal_form = format!(
        "P_n = (alpha_n x + beta_n) P_(n-1) + kappa_n (x - L)(x - U) P'_(n-1) for n >= {n_first}, \
         with L = {lower_text}, U = {upper_text}, alpha_n = {}, beta_n = {}, kappa_n = {}",
        over_f(x_column(c0, 1)).text(),
        over_f(x_column(c0, 0)).text(),
        over_f(gamma.clone()).text(),
    );
    let interval_text = format!("[{lower_text}, {upper_text}]");

    let mut hypotheses: Vec<LiftHypothesis> = Vec::new();
    hypotheses.extend(denominator_hypothesis(f, n_first));
    let seed = &seeds[seeds.len() - 1];
    let seed_multiplicities = closed_interval_multiplicities(seed, &lower, &upper);
    hypotheses.push(exact_hypothesis(
        "seed_roots_in_interval",
        &format!("P_{start} is nonzero and real-rooted with every root in {interval_text}"),
        match seed_multiplicities {
            Some((m_lower, m_upper)) => Ok(format!(
                "P_{start} has degree {}, the root {lower_text} with multiplicity {m_lower} and \
                 the root {upper_text} with multiplicity {m_upper}; checked exactly",
                vec_poly::degree(&vec_poly::trim(seed)).unwrap_or(0)
            )),
            None => Err((
                Some(start),
                format!(
                    "P_{start} is zero, not real-rooted, or has a root outside {interval_text} \
                     (exact)"
                ),
            )),
        },
    ));
    // Earlier seeds, checked in the closed interval; record which ends they
    // attain.
    let mut seeds_attain = [false, false];
    if seeds.len() > 1 {
        let last = start - 1;
        let failing = seeds[..seeds.len() - 1].iter().position(|row| {
            match closed_interval_multiplicities(row, &lower, &upper) {
                Some((m_lower, m_upper)) => {
                    seeds_attain[0] |= m_lower > 0;
                    seeds_attain[1] |= m_upper > 0;
                    false
                }
                None => true,
            }
        });
        hypotheses.push(exact_hypothesis(
            "earlier_seeds",
            &format!(
                "the seed rows P_{first_index}, ..., P_{last} are nonzero and real-rooted with \
                 every root in {interval_text}"
            ),
            match failing {
                None => Ok(format!(
                    "checked exactly for P_{first_index}, ..., P_{last}"
                )),
                Some(offset) => Err((
                    Some(first_index + offset),
                    format!(
                        "P_{} is zero, not real-rooted, or has a root outside {interval_text} \
                         (exact)",
                        first_index + offset
                    ),
                )),
            },
        ));
    }
    hypotheses.push(sign_hypothesis(
        "slope_positive",
        &format!(
            "a_n = [x^1] c_0(n, x) / gamma(n) > 0 for every n >= {n_first} (the coefficient a of \
             the normalized step, so deg P_n = deg P_(n-1) + 1)"
        ),
        "a_n",
        &IndexFunction::new(x_column(c0, 1), gamma.clone()),
        n_first,
        IndexSignCondition::Positive,
    ));
    let mut attained = [true, true];
    if let Some((m_lower, m_upper)) = seed_multiplicities {
        let scaled_gamma = vec_poly::scale(&gamma, &width);
        let e_upper = IndexFunction::new(evaluate_in_x(c0, &upper), scaled_gamma.clone());
        let e_lower = IndexFunction::new(vec_poly::neg(&evaluate_in_x(c0, &lower)), scaled_gamma);
        let upper_decision = endpoint_threshold(
            "upper_endpoint_threshold",
            &format!("U = {upper_text}"),
            &e_upper,
            n_first,
            m_upper,
        );
        let lower_decision = endpoint_threshold(
            "lower_endpoint_threshold",
            &format!("L = {lower_text}"),
            &e_lower,
            n_first,
            m_lower,
        );
        attained = [
            seeds_attain[0] || m_lower > 0 || lower_decision.attained,
            seeds_attain[1] || m_upper > 0 || upper_decision.attained,
        ];
        hypotheses.push(upper_decision.hypothesis);
        hypotheses.push(lower_decision.hypothesis);
    }
    let certified = RealInterval {
        lower: Some(lower.clone()),
        lower_closed: attained[0],
        upper: Some(upper.clone()),
        upper_closed: attained[1],
    };
    let root_region = RootInterval::from_real(&certified);
    if let Some(requested) = &options.root_region {
        let lower_ok = requested
            .lower
            .as_ref()
            .is_none_or(|a| *a < lower || (*a == lower && !certified.lower_closed));
        let upper_ok = requested
            .upper
            .as_ref()
            .is_none_or(|b| *b > upper || (*b == upper && !certified.upper_closed));
        let statement = format!(
            "the requested region {} contains the certified region {root_region}",
            open_root_interval(requested)
        );
        hypotheses.push(if lower_ok && upper_ok {
            exact_hypothesis("requested_root_region", &statement, Ok("exact".into()))
        } else {
            unresolved_hypothesis(
                "requested_root_region",
                &statement,
                format!("this theorem certifies the region {root_region} only"),
            )
        });
    }
    report.hypotheses = hypotheses;
    report.lean_declarations = LEAN_DECLARATIONS.iter().map(|s| s.to_string()).collect();
    Ok(report.finish(RealRootednessTheoremConclusion {
        first_index,
        zero_residue_classes: None,
        interlacing_from: Some(n_first),
        interlacing: Some(RowInterlacing::Weak),
        simple_roots: false,
        statement: format!(
            "every P_n, n >= {first_index}, is nonzero and real-rooted with every root in \
             {root_region}; deg P_n = deg P_(n-1) + 1 and P_(n-1) interlaces P_n (Ma--Wang) for \
             n >= {n_first}"
        ),
        root_region: Some(root_region),
    }))
}
