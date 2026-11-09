//! Three-term recurrences `P_n = A_n(x) P_(n-1) + B_n(x) P_(n-2)` with
//! coefficients depending on `x`: real-rootedness and weak interlacing of
//! consecutive rows by Liu and Wang's Theorem 2.1 with the lag row as the
//! interlacer (their Corollary 2.4 with `b_n = 0`, and with a derivative term
//! `V_n(x) P'_(n-1)` their Theorem 2.3 and Corollary 2.4), optionally with
//! every root in `(-inf, 0)`, `(-inf, 0]` or a root window.
//!
//! # The theorem
//!
//! L. L. Liu and Y. Wang, Adv. Appl. Math. 38 (2007) 542-560, Theorem 2.1:
//! let `F = a f + b g` with `deg F` equal to `deg f` or `deg f + 1`, `f` and
//! `g` real-rooted with `g << f` (`g` alternates left of or interlaces `f`,
//! the weak relation of [`crate::check_weak_interlacing`]), `F` and `g` with
//! leading coefficients of one sign, and `b(r) <= 0` at every root `r` of
//! `f`.  Then `F` is real-rooted and `f << F`.
//!
//! For the rule `f(n) P_n = c_1(n, x) P_(n-1) + c_2(n, x) P_(n-2)`, applied
//! for `n >= N`, with `f` free of `x`, `deg_x c_1 <= 1`, `deg_x c_2 <= 2` and
//! both terms nonzero, write `A_n = c_1 / f`, `B_n = c_2 / f`.  Corollary 2.4
//! (with `b_n = 0`) iterates the step with `F = P_n`, `f = P_(n-1)`,
//! `g = P_(n-2)`.  The check proves, for every `n >= N`:
//!
//! 1. `f(n) != 0`;
//! 2. the seeds: `P_(N-2)`, `P_(N-1)` are nonzero and real-rooted,
//!    `deg P_(N-1) - deg P_(N-2)` is `0` or `1`, and `P_(N-2) << P_(N-1)`
//!    (exact);
//! 3. the degree law: `deg P_n - deg P_(n-1)` is `0` or `1` and the leading
//!    coefficients of `P_n` and `P_(n-2)` have one sign;
//! 4. the sign condition `B_n(r) <= 0` at every root `r` of `P_(n-1)`,
//!    either because `B_n(x) <= 0` for every real `x` ("everywhere"), or
//!    because `B_n(x) <= 0` for every `x <= 0` and every root of every row
//!    is `<= 0` (root location, below).
//!
//! By induction on `n`, every `P_n`, `n >= N - 2`, is then real-rooted and
//! `P_(n-1) << P_n` for `n >= N - 1`: the hypotheses of Theorem 2.1 at step
//! `n` are 2 (for `n = N`) or the conclusion of step `n - 1`, 3 and 4.
//!
//! # Normalizations
//!
//! - *Signs.*  The rule is linear, so with `s_n = +-1` the rows `s_n P_n`
//!   satisfy the rule with `A_n` replaced by `(s_n / s_(n-1)) A_n` and `B_n`
//!   unchanged when `s_n = s_(n-2)`.  The check takes `s_n` constant
//!   (`-1` when `P_(N-1)` has a negative leading coefficient) or alternating
//!   (when the leading coefficients of the two seeds differ in sign, so `A_n`
//!   becomes `-A_n`), so that both seeds have positive leading coefficients.
//!   Roots, interlacing and the sign condition do not change.
//! - *A power of `x`.*  When `x^m` divides both seeds, it divides every row
//!   (the rule has no derivative), and `Q_n = P_n / x^m` satisfies the same
//!   rule.  The check runs on the `Q_n`; `P_n` has the roots of `Q_n` and the
//!   root `0` with multiplicity `m`, and weak interlacing is preserved by a
//!   common root (the merged decreasing root sequence stays alternating).
//!
//! # The degree law
//!
//! With `d = deg P_(n-1)` and `delta = d - deg P_(n-2)`, the candidate top
//! degrees of the two terms are `deg_x c_1 + d` and `deg_x c_2 + d - delta`,
//! so `delta_n = deg P_n - deg P_(n-1)` and the leading coefficient
//! `l_n = u_n l_(n-1) + v_n l_(n-2)` follow from `delta` and the coefficient
//! degrees (`u_n` and `v_n` are the coefficients of `A_n` and `B_n` of the
//! relevant powers of `x`, or zero).  The reachable states `delta` from the
//! seeds give one or two such steps (full growth `delta = 1`; half growth,
//! where `delta` alternates, for `A_n` free of `x` and `deg_x B_n = 1`).
//! With positive seed leading coefficients, `l_n > 0` for every `n` follows
//! from either
//!
//! - `u_n >= 0`, `v_n >= 0` and `u_n + v_n > 0` for every `n >= N`, for every
//!   reachable step, or
//! - (one state, a genuine three-term recurrence for `l_n`) a positive-ratio
//!   induction: a rational function `g(n) > 0` with `l_N / l_(N-1) >= g(N)`,
//!   `v_n <= 0` and `u_n + v_n / g(n-1) - g(n) >= 0` for `n > N`; then
//!   `l_n / l_(n-1) = u_n + v_n l_(n-2)/l_(n-1) >= u_n + v_n / g(n-1) >= g(n)`
//!   by induction.  The candidates are the exact ratio of the computed
//!   leading coefficients (when it is a rational function of `n` of degree at
//!   most two over two) and `u_n / 2`.
//!
//! Every condition is decided by the sign checker
//! ([`crate::decide_index_rational_function_sign`]); a failed sufficient
//! test is not a disproof.
//!
//! In the identity coordinates, when this schedule is unresolved (a step
//! table outside the cases above, or a sequence `l_n` that neither the
//! one-signed test nor the positive-ratio induction covers), the hypothesis
//! `degree_law` is taken from the shared module
//! ([`crate::recurrence::certify_leading_coefficient`]) when it proves
//! degree steps `0` or `1` (only `1` when the seeds and coefficient degrees
//! give full growth) and positive leading coefficients after the sign
//! normalization for every `n >= N - 2`.
//!
//! # The sign condition
//!
//! "`B_n(x) <= 0` for every real `x`" is decided for `deg_x c_2 = 0` by the
//! sign of `B_n`, and for `deg_x c_2 = 2` (with `B_n = p_2 x^2 + p_1 x + p_0`)
//! by `p_2 <= 0`, `p_0 <= 0` and `p_1^2 - 4 p_2 p_0 <= 0` (RealRooted
//! `eval_nonpos_of_natDegree_le_two`); when the discriminant vanishes
//! identically, `B_n = -(c_n x + d_n)^2` is the case of Liu and Wang's
//! Proposition 3.3.  "`B_n(x) <= 0` for every `x <= 0`" is decided for a
//! separable `c_2 = gamma(n) w(x)` with `w` without roots in `(-inf, 0)` by
//! the sign of `gamma / f`, and otherwise by the monomial test after
//! `x = -t`: every `(-1)^k [x^k] B_n <= 0` (for degree two this is
//! `p_0 <= 0`, `p_1 >= 0`, `p_2 <= 0`, RealRooted
//! `eval_nonpos_of_nonpos_of_natDegree_le_two`).
//!
//! # Root location
//!
//! Corollary 2.4 assumes that every row has nonnegative coefficients.  The
//! check proves the root location for every row instead, in one of two ways
//! (both on the normalized rows, positive leading coefficients):
//!
//! - *root-sign induction* (every root `< 0`): both seeds have every root
//!   `< 0`, and `P_n(0) > 0` for every `n >= N`, decided as a positive
//!   sequence `p_n = A_n(0) p_(n-1) + B_n(0) p_(n-2)` as in the degree law.
//!   Then every row has every root `< 0`: if `P_(n-1)` has every root `< 0`
//!   and `P_(n-1) << P_n` (Theorem 2.1, which needs only roots `<= 0` of
//!   `P_(n-1)` for the sign condition), every root of `P_n` except the
//!   largest is `<=` the largest root of `P_(n-1)`, so `< 0`; the largest is
//!   `< 0` too, since otherwise `P_n(0)` would be `0` or have the sign
//!   opposite to the leading coefficient.
//! - *nonnegative coefficients* (every root `<= 0`): both seeds and every
//!   coefficient `[x^k] A_n`, `[x^k] B_n` (`n >= N`) are nonnegative; then
//!   every row has nonnegative coefficients by induction, and a real-rooted
//!   row with nonnegative coefficients has every root `<= 0` (the
//!   RealRooted lemma `threeTerm_hasNonnegCoeffs`).
//!
//! In the "everywhere" case the root location is not needed for
//! real-rootedness; the check still tries it and certifies the root region
//! when it holds (the hypotheses are then listed), and otherwise concludes
//! real-rootedness and interlacing without a region.  The region is
//! `(-inf, 0)` from the root-sign induction, `(-inf, 0]` from nonnegative
//! coefficients or after dividing by `x^m`.
//!
//! # A derivative term
//!
//! The rule may also have a term `c_d(n, x) P'_(n-1)` (`deg_x c_d <= 2`,
//! `V_n = c_d / f`): Liu and Wang's Theorem 2.3 with the two interlacers
//! `g_1 = P'_(n-1)` (`P'_(n-1) << P_(n-1)` by Rolle's theorem, leading
//! coefficient `deg P_(n-1) l_(n-1) > 0`; the term vanishes for a constant
//! `P_(n-1)`) and `g_2 = P_(n-2)`, with the sign condition for `V_n` as for
//! `B_n` (the "everywhere" case needs both).  In the degree law the
//! multiplier of `l_(n-1)` is `[x^k] A_n + deg P_(n-1) [x^(k+1)] V_n`, with
//! `deg P_(n-1) = deg P_(N-1) + n - N` on full growth and constant degree;
//! with two states (half growth, the degrees `0, 0, 1, 1, 2, ...` of
//! OEIS A008299) the check uses `deg P_(n-1) >= deg P_(N-1)` and asks
//! `[x^(k+1)] V_n >= 0`.  The root-sign induction needs `V_n(0) = 0`, dividing
//! by `x^m` needs `x | c_d` (then `Q_n = P_n / x^m` satisfies the rule with
//! `A_n + m V_n / x`), and the windows, the gauges and the reversal are not
//! used with a derivative term.
//!
//! # Root windows
//!
//! When neither sign case closes (and there is no derivative term), the
//! check tries the windows `W = (-inf, U]`, `[L, +inf)` and `[L, U]` whose
//! ends are rational roots of `w` for a separable `c_2 = gamma(n) w(x)`
//! (RealRooted `threeTerm_interlaces_of_roots_mem_Icc` and
//! `threeTerm_interlaces_of_roots_le`).  Hypotheses: `B_n(x) <= 0` on `W`
//! (every coefficient of `B_n(U - t)`, or of `B_n(L + t)`, is `<= 0`; or
//! `[x^2] B_n >= 0` with `B_n(L), B_n(U) <= 0`); right of a finite `U`,
//! `A_n > 0` and `B_n >= 0` (every coefficient of `A_n(U + t)` and
//! `B_n(U + t)` is `>= 0`, and `A_n(U + 1) > 0`); left of a finite `L`, full
//! degree growth, `A_n < 0` and `B_n >= 0` (the same on `A_n(L - t)`,
//! `B_n(L - t)`); and both seeds with every root in `W`.  By induction every
//! row has every root in `W`: for `x > U`, `f(x), g(x) > 0`, so
//! `F(x) = a f + b g > 0`; for `x < L`, `(-1)^(deg f) f(x)` and
//! `(-1)^(deg g) g(x)` are positive and `deg f = deg g + 1`, so
//! `(-1)^(deg g) F(x) = -a (-1)^(deg f) f + b (-1)^(deg g) g > 0`; Theorem 2.1
//! applies at each step since the roots of `P_(n-1)` lie in `W`.
//!
//! # Coordinates
//!
//! The check runs in the identity coordinates first and, when it does not
//! prove a root region there (and the rule has no derivative term), in the
//! coordinates `Q_n = P_n / (x + v)^(e + n - N + 2)` for every rational `v`
//! with `x + v` dividing `c_1` and `(x + v)^2` dividing `c_2` (`e` from the
//! seeds; `Q_n = (A_n / (x + v)) Q_(n-1) + (B_n / (x + v)^2) Q_(n-2)`), and
//! `Q_n = x^(d + n - N + 2) P_n(1/x)` with `d = deg P_(N-2)` (an upper bound
//! for the degrees by induction, so `Q_n = x A_n(1/x) Q_(n-1) + x^2 B_n(1/x)
//! Q_(n-2)` is an identity of polynomials).  Real-rootedness of the `Q_n`
//! gives that of the `P_n`; the region of the `P_n` adds the root `-v`
//! (gauge) or is `(-inf, 0]` when every root of every `Q_n` is `<= 0`
//! (reversal: the nonzero roots invert, and a degree below the nominal one is
//! a root `0` of `P_n`).  No interlacing of the `P_n` is claimed in these
//! coordinates.  The first report proved with a root region is returned,
//! else the first proved one, else the report in the identity coordinates.
//!
//! # Earlier rows and the induction start
//!
//! Seed rows before `P_(N-2)` are checked exactly: each nonzero and
//! real-rooted, or (with a certified region) zero or real-rooted with every
//! root in the region; zero rows are named in the conclusion.  When the check
//! fails at the supplied seeds, the induction start is moved as for the other
//! theorems ([`RealRootednessTheoremOptions::maximum_start_delay`]), which
//! also handles seeds whose degree step is not the steady one; a start proved
//! without a root region is replaced by the first later start that proves a
//! region.
//!
//! # Lean correspondence
//!
//! For full degree growth the hypotheses mirror the RealRooted route `sign`
//! (<https://github.com/PerAlexandersson/RealRooted>,
//! `ThreeTermRecurrence/Interlacing.lean`): `threeTerm_interlaces_of_eval_nonpos`
//! ("everywhere") and `threeTerm_interlaces_of_nonnegCoeffs` (roots `<= 0`),
//! both through `threeTerm_interlaces`.  The Lean statements take the degree
//! law `D0 + n`, positive leading coefficients and nonnegative coefficients
//! as hypotheses for every row; here they are proved, so a replay must
//! establish them separately.  Other degree schedules (half growth) have no
//! recorded route; the single step is
//! `RealRooted.MaWang.strictInterl_of_interlaces_evalCoeff_nonpos`
//! with its degree window `deg f <= deg F <= deg f + 1`.

use super::derivative_lag::{negative_axis_sign_hypothesis, order_at_zero, x_polynomial_text};
use super::three_term::fit_index_rational_function;
use super::*;

mod proportional_rows;

const REFERENCE: &str = "L. L. Liu and Y. Wang, Adv. Appl. Math. 38 (2007) 542-560, Theorem 2.1 \
     and Corollary 2.4 (A unified approach to polynomial sequences with only real zeros; the \
     lag row as the interlacer, b_n = 0), weak interlacing of consecutive rows, with \
     nonnegativity of the coefficients replaced, where a root region is certified, by an exact \
     root-location induction";

/// The addition to [`REFERENCE`] when `B_n = -(c_n x + d_n)^2`.
const PROPOSITION_3_3: &str = "; Proposition 3.3 (the case B_n(x) = -(c_n x + d_n)^2)";

/// The addition to [`REFERENCE`] for a derivative term `V_n P'_(n-1)`.
const THEOREM_2_3: &str = "; Theorem 2.3 with the interlacers P'_(n-1) and P_(n-2) (Corollary 2.4 \
     with b_n = V_n)";

/// Lean route of the full-growth rows (RealRooted).
const LEAN_ROUTE_SIGN: &str = "sign";
/// Lean route of the full-growth rows in a root window (RealRooted).
const LEAN_ROUTE_WINDOW: &str = "window";

/// Check Liu and Wang's three-term theorem (module documentation of the
/// `liu_wang_three_term` submodule, summarized here) for the family given by
/// `rule`, `first_index` and the seed rows `seeds` (the rule applies from
/// `N = first_index + seeds.len()` on; at least two seed rows).
///
/// The rule must be `f(n) P_n = c_1(n, x) P_(n-1) + c_2(n, x) P_(n-2)`, or
/// with an additional derivative term `c_d(n, x) P'_(n-1)` (module
/// documentation, "A derivative term"), with `f` free of `x`,
/// `deg_x c_1 <= 1`, `deg_x c_d, deg_x c_2 <= 2` and a nonzero lag term.
/// With `A_n = c_1 / f`, `B_n = c_2 / f` (and `V_n = c_d / f`, whose sign
/// condition is that of `B_n`) the hypotheses, each decided exactly for every
/// `n >= N`, are:
///
/// - `f(n) != 0` (sign checker);
/// - `P_(N-2)`, `P_(N-1)` nonzero and real-rooted, of degrees `e` and `e` or
///   `e + 1`, with `P_(N-2) << P_(N-1)` (exact weak interlacing);
/// - the degree law `deg P_n - deg P_(n-1)` in `{0, 1}` with the leading
///   coefficients of `P_n` and `P_(n-2)` of one sign (sign checker, by the
///   step coefficients or a positive-ratio induction);
/// - `B_n(x) <= 0` for every real `x`, or for every `x <= 0` together with
///   a root location for every row (every root `< 0` by a root-sign induction
///   with `P_n(0) > 0`, or every root `<= 0` by nonnegative coefficients), or
///   on a root window `W` that every row stays in (module documentation,
///   "Root windows");
/// - earlier seed rows real-rooted (exact), and a requested root region
///   containing the certified one.
///
/// When all hold, every `P_n`, `n >= first_index`, is real-rooted and
/// `P_(n-1)` interlaces `P_n` in the weak sense of Liu and Wang
/// ([`RowInterlacing::Weak`]) for `n >= N - 1`; with a root location or a
/// window, every root lies in the certified region
/// ([`RealRootednessTheoremConclusion::root_region`]).  The rows divided by a
/// power of `x`, the reversed rows and gauged rows are tried as well
/// (module documentation, "Coordinates"; the last two without an
/// interlacing claim).  The induction start is moved as in
/// [`check_liu_wang_theorem`] when the check fails at the supplied seeds.
///
/// # Errors
///
/// Errors of the seed record itself (as for [`check_liu_wang_theorem`]).
///
/// # Example
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::recurrence::{
///     check_liu_wang_three_term_theorem, BigRational, BivarPoly, Recurrence, RecurrenceTerm,
///     RealRootednessTheoremOptions, RowInterlacing,
/// };
///
/// let q = |v: i64| BigRational::from_integer(BigInt::from(v));
/// // Lucas polynomials (OEIS A034807): P_n = P_(n-1) + x P_(n-2), P_0 = 2, P_1 = 1.
/// let rule = Recurrence {
///     terms: vec![
///         RecurrenceTerm::new(1, 0, BivarPoly::from_coefficients_normalized(vec![vec![q(1)]])),
///         RecurrenceTerm::new(
///             2,
///             0,
///             BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(1)]]),
///         ),
///     ],
///     denominator: None,
///     inhomogeneous: None,
/// };
/// let seeds = [vec![q(2)], vec![q(1)]];
/// let report = check_liu_wang_three_term_theorem(
///     &rule,
///     0,
///     &seeds,
///     &RealRootednessTheoremOptions::default(),
/// )
/// .unwrap();
/// assert!(report.is_proved());
/// let conclusion = report.conclusion.as_ref().unwrap();
/// assert_eq!(conclusion.interlacing, Some(RowInterlacing::Weak));
/// assert_eq!(conclusion.root_region.as_ref().unwrap().to_string(), "(-inf, 0)");
/// report.verify().unwrap();
/// ```
pub fn check_liu_wang_three_term_theorem(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<BigRational>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    check_real_rootedness_theorem(
        RealRootednessTheorem::LiuWangThreeTerm,
        rule,
        first_index,
        seeds,
        options,
    )
}

/// How the sign condition `B_n(r) <= 0` at the roots of `P_(n-1)` holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LagSign {
    /// `B_n(x) <= 0` for every real `x`.
    Everywhere,
    /// `B_n(x) <= 0` for every `x <= 0`; needs the root location.
    NonpositiveAxis,
}

/// How every root of every row is kept in `(-inf, 0]`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RootLocation {
    /// Every root `< 0` (seeds, and `P_n(0) > 0` for every `n`).
    Negative,
    /// Every root `<= 0` (nonnegative coefficients of the seeds and of the
    /// multipliers).
    NonnegativeCoefficients,
}

/// One step type of the degree law: from the degree difference `from` of
/// the two previous rows to `to`, with `l_n = u_n l_(n-1) + v_n l_(n-2)`.
struct DegreeStep {
    from: usize,
    to: usize,
    u: IndexFunction,
    v: IndexFunction,
    u_text: String,
    v_text: String,
}

/// The coordinates the theorem runs in (module documentation,
/// "Coordinates"); `start = N - 2` is the first row of the induction.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Coordinates {
    /// `Q_n = P_n`.
    Identity,
    /// `Q_n = P_n / (x + v)^(exponent + n - start)`.
    Gauge { v: Q, exponent: usize },
    /// `Q_n = x^(degree + n - start) P_n(1/x)`.
    Reversal { degree: usize },
}

impl Coordinates {
    /// The rule for `Q`: `(c_1, c_2)` of `f(n) Q_n = c_1 Q_(n-1) + c_2 Q_(n-2)`.
    fn transform_rule(&self, c1: &Bivariate, c2: &Bivariate) -> Option<(Bivariate, Bivariate)> {
        match self {
            Self::Identity => Some((c1.clone(), c2.clone())),
            Self::Gauge { v, .. } => {
                let root = -v.clone();
                let c1 = map_x_rows(c1, |row| vec_poly::exact_div_by_t_minus(row, &root))?;
                let c2 = map_x_rows(c2, |row| vec_poly::exact_div_by_t_minus(row, &root))?;
                let c2 = map_x_rows(&c2, |row| vec_poly::exact_div_by_t_minus(row, &root))?;
                Some((c1, c2))
            }
            Self::Reversal { .. } => Some((
                map_x_rows(c1, |row| vec_poly::reverse(row, 1))?,
                map_x_rows(c2, |row| vec_poly::reverse(row, 2))?,
            )),
        }
    }

    /// `Q_index` from `P_index`, `None` when it is not a polynomial.
    fn transform_row(&self, p: &[Q], index: usize, start: usize) -> Option<Vec<Q>> {
        let p = vec_poly::trim(p);
        match self {
            Self::Identity => Some(p),
            Self::Gauge { v, exponent } => {
                let root = -v.clone();
                (0..exponent + index - start).try_fold(p, |row, _| {
                    if vec_poly::is_zero(&row) {
                        None
                    } else {
                        vec_poly::exact_div_by_t_minus(&row, &root)
                    }
                })
            }
            Self::Reversal { degree } => vec_poly::reverse(&p, degree + index - start),
        }
    }

    /// `Q_n` in terms of `P_n`, for the normal form.
    fn text(&self, start: usize) -> Option<String> {
        match self {
            Self::Identity => None,
            Self::Gauge { v, exponent } => Some(format!(
                "Q_n = P_n / {}^(n - {}) for n >= {start}",
                if v.is_zero() {
                    "x".to_string()
                } else {
                    format!("(x{})", signed_rational(v))
                },
                start as i64 - *exponent as i64
            )),
            Self::Reversal { degree } => Some(format!(
                "Q_n = x^(n - {}) P_n(1/x) for n >= {start} (the reversal at the nominal degree \
                 deg P_{start} + n - {start}, an upper bound for deg P_n)",
                start as i64 - *degree as i64
            )),
        }
    }

    /// The region for the roots of `P_n` (`n >= start`) given that of
    /// `Q_n`.
    fn region_of_rows(&self, region: Option<RealInterval>) -> Option<RealInterval> {
        let region = region?;
        match self {
            Self::Identity => Some(region),
            // Every root of Q is <= 0, so every root of P = x^(...) Q(1/x)
            // is <= 0 (0 when deg Q_n is below the nominal degree).
            Self::Reversal { .. } => within_nonpositive_ray(&region).then(|| RealInterval {
                lower: None,
                lower_closed: false,
                upper: Some(Q::zero()),
                upper_closed: true,
            }),
            Self::Gauge { v, .. } => Some(hull_with_point(region, &-v.clone())),
        }
    }
}

/// Whether every point of `region` is `<= 0`.
fn within_nonpositive_ray(region: &RealInterval) -> bool {
    region
        .upper
        .as_ref()
        .is_some_and(|upper| !upper.is_positive())
}

/// The smallest interval containing `region` and the point `point` (closed
/// there).
fn hull_with_point(region: RealInterval, point: &Q) -> RealInterval {
    let (lower, lower_closed) = match region.lower {
        None => (None, false),
        Some(lower) if *point < lower => (Some(point.clone()), true),
        Some(lower) if *point == lower => (Some(lower), true),
        Some(lower) => (Some(lower), region.lower_closed),
    };
    let (upper, upper_closed) = match region.upper {
        None => (None, false),
        Some(upper) if *point > upper => (Some(point.clone()), true),
        Some(upper) if *point == upper => (Some(upper), true),
        Some(upper) => (Some(upper), region.upper_closed),
    };
    RealInterval {
        lower,
        lower_closed,
        upper,
        upper_closed,
    }
}

/// The seed rows before the induction, each zero or real-rooted with every
/// root in `region`, and the indices of the zero rows.
fn earlier_rows_in_region(
    rows: &[Vec<Q>],
    first_index: usize,
    region: &RealInterval,
) -> (Option<LiftHypothesis>, Vec<usize>) {
    if rows.is_empty() {
        return (None, Vec::new());
    }
    let last = first_index + rows.len() - 1;
    let region_text = RootInterval::from_real(region).to_string();
    let mut zero_rows = Vec::new();
    let mut failing = None;
    for (offset, row) in rows.iter().enumerate() {
        let row = vec_poly::trim(row);
        if vec_poly::is_zero(&row) {
            zero_rows.push(first_index + offset);
        } else if failing.is_none() && !has_real_roots_in_region(&row, region) {
            failing = Some(first_index + offset);
        }
    }
    let outcome = match failing {
        None => Ok(format!(
            "checked exactly for P_{first_index}, ..., P_{last}{}",
            if zero_rows.is_empty() {
                String::new()
            } else {
                format!(" ({} zero)", zero_rows.len())
            }
        )),
        Some(index) => Err((
            Some(index),
            format!("P_{index} is not real-rooted with every root in {region_text} (exact)"),
        )),
    };
    (
        Some(exact_hypothesis(
            "earlier_seeds",
            &format!(
                "the seed rows P_{first_index}, ..., P_{last} are each zero or real-rooted with \
                 every root in {region_text}"
            ),
            outcome,
        )),
        zero_rows,
    )
}

/// `" + v"` or `" - |v|"`.
fn signed_rational(value: &Q) -> String {
    if value.is_negative() {
        format!(" - {}", format_rational_coeff(&-value.clone()))
    } else {
        format!(" + {}", format_rational_coeff(value))
    }
}

/// Apply `map` to every `x`-polynomial row `[x^j] c` of each power `n^i`
/// (zero rows stay zero).
fn map_x_rows(c: &Bivariate, map: impl Fn(&[Q]) -> Option<Vec<Q>>) -> Option<Bivariate> {
    c.iter()
        .map(|row| {
            let row = vec_poly::trim(row);
            if vec_poly::is_zero(&row) {
                Some(Vec::new())
            } else {
                let mut mapped = map(&row)?;
                vec_poly::strip_trailing_zeros(&mut mapped);
                Some(mapped)
            }
        })
        .collect()
}

/// The multiplicity of the root `root` of a nonzero `p`.
fn root_order(p: &[Q], root: &Q) -> usize {
    let mut p = vec_poly::trim(p);
    let mut order = 0;
    while !vec_poly::is_zero(&p) {
        match vec_poly::exact_div_by_t_minus(&p, root) {
            Some(quotient) => {
                p = quotient;
                order += 1;
            }
            None => break,
        }
    }
    order
}

/// The gauges `x + v` (rational `v`) with `x + v` dividing `c_1` and
/// `(x + v)^2` dividing `c_2` as polynomials in `x` (for every `n`), from
/// the linear factors of the first nonzero row of `c_1`.
fn gauge_candidates(c1: &Bivariate, c2: &Bivariate) -> Vec<Q> {
    let Some(row) = c1
        .iter()
        .map(|row| vec_poly::trim(row))
        .find(|row| vec_poly::degree(row).is_some_and(|d| d >= 1))
    else {
        return Vec::new();
    };
    let Ok(factorization) = crate::factorization::irreducible_factorization_rational_coeffs(&row)
    else {
        return Vec::new();
    };
    factorization
        .factors
        .iter()
        .filter(|(factor, _)| factor.len() == 2)
        .map(|(factor, _)| Q::new(factor[0].clone(), factor[1].clone()))
        .filter(|v| {
            Coordinates::Gauge {
                v: v.clone(),
                exponent: 0,
            }
            .transform_rule(c1, c2)
            .is_some()
        })
        .collect()
}

/// [`TheoremRoute::is_preferred`]: a proof whose conclusion comes without a
/// root region is replaced by a later start that certifies one.
fn proof_with_root_region(report: &RealRootednessTheoremReport) -> bool {
    !report
        .conclusion
        .as_ref()
        .is_some_and(|conclusion| conclusion.root_region.is_none())
}

/// The registry entry of the Liu--Wang three-term theorem ([`super::registry::ROUTES`]).
pub(super) const ROUTE: TheoremRoute = TheoremRoute {
    theorem: RealRootednessTheorem::LiuWangThreeTerm,
    check_at_supplied_seeds: liu_wang_three_term_report,
    start_policy: StartPolicy::MovesStart,
    is_preferred: proof_with_root_region,
    certificate_name: None,
};

/// The check with the induction at the given (normalized) seeds: the
/// identity coordinates first, then the gauges and the reversal (module
/// documentation, "Coordinates"); the first report proved with a root
/// region, else the first proved report (identity first), else the report
/// in the identity coordinates (or, when the rule does not have the form
/// there, the first other report that is not `not_applicable`).
fn liu_wang_three_term_report(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<BigRational>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    let report = RealRootednessTheoremReport::new(
        RealRootednessTheorem::LiuWangThreeTerm,
        REFERENCE,
        rule,
        first_index,
        seeds,
        options,
    );
    if seeds.len() < 2 {
        return Ok(report.not_applicable("the theorem needs two seed rows"));
    }
    let shape = match recurrence_shape(rule) {
        Ok(shape) => shape,
        Err(reason) => return Ok(report.not_applicable(reason)),
    };
    let allowed = [(1, 0), (1, 1), (2, 0)];
    if shape.terms.keys().any(|key| !allowed.contains(key))
        || !shape.terms.contains_key(&(2, 0))
        || !(shape.terms.contains_key(&(1, 0)) || shape.terms.contains_key(&(1, 1)))
    {
        return Ok(report.not_applicable(
            "the rule is not f(n) P_n = c_1(n, x) P_(n-1) + c_d(n, x) P'_(n-1) + c_2(n, x) \
             P_(n-2) with a nonzero lag term and a nonzero P_(n-1) or P'_(n-1) term",
        ));
    }
    let empty = Bivariate::new();
    let c1 = shape.terms.get(&(1, 0)).unwrap_or(&empty);
    let c2 = &shape.terms[&(2, 0)];
    let derivative = shape.terms.get(&(1, 1));
    let start = first_index + seeds.len() - 2;
    // The rows P_start, P_(start+1), ... generated by the rule (for the
    // positive-ratio certificates and the exact refutations).
    let rows: Option<Vec<Vec<Q>>> = rule
        .generate_rows_rational(seeds, first_index, seeds.len() + GENERATED_ROWS)
        .ok()
        .map(|rows| rows[start - first_index..].to_vec());
    let input = AnalysisInput {
        rule,
        first_index,
        seeds,
        options,
        f: &shape.denominator,
        rows: rows.as_deref(),
        derivative,
    };
    let identity = analyze(report.clone(), &Coordinates::Identity, c1, c2, &input);
    if identity.is_proved() && has_root_region(&identity) {
        return Ok(identity);
    }
    let previous = vec_poly::trim(&seeds[seeds.len() - 2]);
    let current = vec_poly::trim(&seeds[seeds.len() - 1]);
    let mut others = Vec::new();
    // The gauges and the reversal are for rules without a derivative term.
    if derivative.is_none() && !vec_poly::is_zero(&previous) && !vec_poly::is_zero(&current) {
        for v in gauge_candidates(c1, c2) {
            let root = -v.clone();
            let exponent = root_order(&previous, &root).min(root_order(&current, &root));
            let exponent = if root_order(&current, &root) > exponent {
                exponent
            } else {
                match exponent.checked_sub(1) {
                    Some(e) => e,
                    None => continue,
                }
            };
            others.push(Coordinates::Gauge { v, exponent });
        }
        if let Some(degree) = vec_poly::degree(&previous) {
            others.push(Coordinates::Reversal { degree });
        }
    }
    let mut fallback = None;
    for coordinates in others {
        let attempt = analyze(report.clone(), &coordinates, c1, c2, &input);
        if attempt.is_proved() && (has_root_region(&attempt) || !identity.is_proved()) {
            return Ok(attempt);
        }
        if fallback.is_none() && attempt.status != RealRootednessTheoremStatus::NotApplicable {
            fallback = Some(attempt);
        }
    }
    let result = match fallback {
        Some(attempt) if identity.status == RealRootednessTheoremStatus::NotApplicable => attempt,
        _ => identity,
    };
    // The degenerate case of proportional rows (submodule
    // `proportional_rows`), tried only when nothing above proves the rule.
    if !result.is_proved() && derivative.is_none() {
        if let Some(proportional) =
            proportional_rows::proportional_rows_report(report, c1, c2, &input)
        {
            if proportional.is_proved() {
                return Ok(proportional);
            }
        }
    }
    Ok(result)
}

/// Whether a report concludes a root region.
fn has_root_region(report: &RealRootednessTheoremReport) -> bool {
    report
        .conclusion
        .as_ref()
        .is_some_and(|c| c.root_region.is_some())
}

/// The family data shared by the coordinate attempts.
struct AnalysisInput<'a> {
    /// The rule as checked (for the shared degree law).
    rule: &'a Recurrence,
    first_index: usize,
    seeds: &'a [Vec<Q>],
    options: &'a RealRootednessTheoremOptions,
    f: &'a [Q],
    /// `P_(N-2), P_(N-1), ...` generated by the rule, when it can.
    rows: Option<&'a [Vec<Q>]>,
    /// `c_d` of a derivative term `c_d(n, x) P'_(n-1)` (identity coordinates
    /// only).
    derivative: Option<&'a Bivariate>,
}

/// The theorem in the given coordinates.
fn analyze(
    mut report: RealRootednessTheoremReport,
    coordinates: &Coordinates,
    c1: &Bivariate,
    c2: &Bivariate,
    input: &AnalysisInput<'_>,
) -> RealRootednessTheoremReport {
    let AnalysisInput {
        rule,
        first_index,
        seeds,
        options,
        f,
        rows,
        derivative,
    } = *input;
    let n_first = first_index + seeds.len();
    let previous_index = n_first - 2;
    let current_index = n_first - 1;
    let coordinate_text = coordinates.text(previous_index);
    let Some((c1, c2)) = coordinates.transform_rule(c1, c2) else {
        return report.not_applicable("the coordinate change does not apply to the rule");
    };
    let (c1, c2) = (&c1, &c2);
    let Some(beta) = x_degree(c2) else {
        return report.not_applicable("the lag term of the rule is zero");
    };
    let gamma = derivative.and_then(x_degree);
    // The x-degree of the multiplier of P_(n-1), the derivative term
    // counting one less.
    let alpha = x_degree(c1)
        .unwrap_or(0)
        .max(gamma.map_or(0, |g| g.saturating_sub(1)));
    if alpha > 1 || beta > 2 || gamma.is_some_and(|g| g > 2) {
        return report.not_applicable(format!(
            "deg_x c_1 > 1, deg_x c_d > 2 or deg_x c_2 > 2{}: the degree may grow by more than \
             one",
            coordinate_text
                .as_ref()
                .map_or(String::new(), |text| format!(" (after {text})"))
        ));
    }
    let variable = if coordinate_text.is_some() { "Q" } else { "P" };
    let mut hypotheses: Vec<LiftHypothesis> = Vec::new();
    hypotheses.extend(denominator_hypothesis(f, n_first));
    let transformed = (
        coordinates.transform_row(&seeds[seeds.len() - 2], previous_index, previous_index),
        coordinates.transform_row(&seeds[seeds.len() - 1], current_index, previous_index),
    );
    let (previous, current) = match transformed {
        (Some(previous), Some(current))
            if !vec_poly::is_zero(&previous) && !vec_poly::is_zero(&current) =>
        {
            (previous, current)
        }
        (previous, _) => {
            let index = if previous.is_some_and(|p| !vec_poly::is_zero(&p)) {
                current_index
            } else {
                previous_index
            };
            hypotheses.push(exact_hypothesis(
                "seed_rows",
                &format!(
                    "{variable}_{previous_index} and {variable}_{current_index} are nonzero \
                     polynomials"
                ),
                Err((
                    Some(index),
                    format!("{variable}_{index} is zero or not a polynomial"),
                )),
            ));
            report.normal_form = format!(
                "P_n = A_n(x) P_(n-1) + B_n(x) P_(n-2) for n >= {n_first}{}",
                coordinate_text
                    .as_ref()
                    .map_or(String::new(), |text| format!("; {text}"))
            );
            report.hypotheses = hypotheses;
            return report.finish(placeholder_conclusion(first_index));
        }
    };

    // Normalizations: the power of x and the signs.  With a derivative term
    // the power of x is divided out only when x divides c_d; then
    // Q_n = P_n / x^m satisfies the rule with c_1 + m c_d / x.
    let derivative_vanishes_at_zero = derivative.is_none_or(|c| x_column(c, 0).is_empty());
    let power = if derivative_vanishes_at_zero {
        order_at_zero(&previous)
            .unwrap_or(0)
            .min(order_at_zero(&current).unwrap_or(0))
    } else {
        0
    };
    let c1 = &match derivative {
        Some(c) if power > 0 => add_bivariate(c1, &scale_divided_by_x(c, power)),
        _ => c1.clone(),
    };
    let previous_row = previous[power..].to_vec();
    let current_row = current[power..].to_vec();
    let previous_lead = previous_row.last().cloned().unwrap_or_else(Q::zero);
    let current_lead = current_row.last().cloned().unwrap_or_else(Q::zero);
    let negate_all = current_lead.is_negative();
    let alternate = previous_lead.is_negative() != current_lead.is_negative();
    // s_index (true: multiply by -1), with s_(N-1) the sign of the current
    // leading coefficient and s_n = -s_(n-1) when alternating.
    let row_sign = |index: usize| -> bool {
        let parity_flip = alternate && (index - previous_index) % 2 == 0;
        negate_all != parity_flip
    };
    let signed = |row: &[Q], index: usize| -> Vec<Q> {
        if row_sign(index) {
            vec_poly::neg(row)
        } else {
            row.to_vec()
        }
    };
    let seed_previous = signed(&previous_row, previous_index);
    let seed_current = signed(&current_row, current_index);
    let c1n = if alternate {
        negate_bivariate(c1)
    } else {
        c1.clone()
    };
    let cdn = derivative.map(|c| {
        if alternate {
            negate_bivariate(c)
        } else {
            c.clone()
        }
    });
    let cdn = cdn.as_ref();
    let row = |index: usize| {
        if power == 0 {
            format!("{variable}_{index}")
        } else {
            format!("{variable}_{index} / x^{power}")
        }
    };
    let multipliers = format!(
        "A_n(x) = {}, {}B_n(x) = {}",
        x_polynomial_text(&c1n, f),
        cdn.map_or(String::new(), |c| format!(
            "V_n(x) = {}, ",
            x_polynomial_text(c, f)
        )),
        x_polynomial_text(c2, f)
    );
    let mut notes = Vec::new();
    if let Some(text) = &coordinate_text {
        notes.push(format!("{text}, which satisfies the rule below"));
    }
    if power > 0 {
        notes.push(format!(
            "{variable}_n = x^{power} R_n for n >= {previous_index} (x^{power} divides \
             {variable}_{previous_index} and {variable}_{current_index}; R_n satisfies the same \
             rule{})",
            if derivative.is_some() {
                format!(" with A_n increased by {power} V_n / x")
            } else {
                String::new()
            }
        ));
    }
    if negate_all || alternate {
        notes.push(format!(
            "the rows are multiplied by {} so that the seeds have positive leading coefficients{}",
            match (negate_all, alternate) {
                (true, false) => "-1".to_string(),
                (_, true) => format!(
                    "{}(-1)^(n - {current_index})",
                    if negate_all { "-" } else { "" }
                ),
                _ => "1".to_string(),
            },
            if alternate {
                " (A_n and V_n change sign)"
            } else {
                ""
            }
        ));
    }
    report.normal_form = format!(
        "{variable}_n = A_n(x) {variable}_(n-1) + {}B_n(x) {variable}_(n-2) for n >= {n_first}, \
         with {multipliers}{}",
        if cdn.is_some() {
            format!("V_n(x) {variable}'_(n-1) + ")
        } else {
            String::new()
        },
        if notes.is_empty() {
            String::new()
        } else {
            format!("; {}", notes.join("; "))
        }
    );

    // The seeds.
    let previous_degree = vec_poly::degree(&seed_previous).unwrap_or(0);
    let current_degree = vec_poly::degree(&seed_current).unwrap_or(0);
    let degree_difference = current_degree.checked_sub(previous_degree);
    let seed_outcome = if !is_real_rooted_rational(&seed_previous) {
        Err((
            Some(previous_index),
            format!("{} is not real-rooted (exact)", row(previous_index)),
        ))
    } else if !is_real_rooted_rational(&seed_current) {
        Err((
            Some(current_index),
            format!("{} is not real-rooted (exact)", row(current_index)),
        ))
    } else if !matches!(degree_difference, Some(0 | 1)) {
        Err((
            Some(current_index),
            format!(
                "deg {} = {current_degree} and deg {} = {previous_degree} do not differ by 0 or 1",
                row(current_index),
                row(previous_index)
            ),
        ))
    } else if crate::check_weak_interlacing_bigint_coeffs(
        &integer_coefficients(&seed_previous),
        &integer_coefficients(&seed_current),
    ) != Some(true)
    {
        Err((
            Some(current_index),
            format!(
                "{} does not interlace {} (exact weak interlacing test)",
                row(previous_index),
                row(current_index)
            ),
        ))
    } else {
        Ok(format!(
            "degrees {previous_degree} and {current_degree}; real-rooted and {} << {} checked \
             exactly",
            row(previous_index),
            row(current_index)
        ))
    };
    let seeds_ok = seed_outcome.is_ok();
    hypotheses.push(exact_hypothesis(
        "seed_rows",
        &format!(
            "{} and {} are nonzero and real-rooted, of degrees e and e or e + 1, and {} << {} \
             (weak interlacing)",
            row(previous_index),
            row(current_index),
            row(previous_index),
            row(current_index)
        ),
        seed_outcome,
    ));

    // The normalized rows from the induction start on.
    let generated: Option<Vec<Vec<Q>>> = rows.and_then(|rows| {
        rows.iter()
            .enumerate()
            .map(|(offset, p)| {
                let index = previous_index + offset;
                let q = coordinates.transform_row(p, index, previous_index)?;
                if q.len() > power && q[..power].iter().any(|c| !c.is_zero()) {
                    return None;
                }
                let quotient = if q.len() > power {
                    q[power..].to_vec()
                } else {
                    vec![Q::zero()]
                };
                Some(signed(&quotient, index))
            })
            .collect()
    });

    // The degree law.
    let local_degree_law = match degree_difference.filter(|_| seeds_ok) {
        None => vec![unresolved_hypothesis(
            "degree_law",
            &format!(
                "deg {variable}_n - deg {variable}_(n-1) is 0 or 1 and {variable}_n, \
                 {variable}_(n-2) have leading coefficients of one sign for every n >= {n_first}"
            ),
            "not checked: the seeds do not satisfy the seed hypothesis".into(),
        )],
        Some(_) => degree_law_hypotheses(
            alpha,
            beta,
            (previous_degree, current_degree),
            (&c1n, cdn, c2, f),
            n_first,
            generated.as_deref(),
        ),
    };
    let expects_full_growth = degree_difference == Some(1) && (alpha == 1 || beta == 2);
    // In the identity coordinates, the shared leading-coefficient module
    // where the degree-step state machine is too narrow: degree steps 0 or 1
    // and positive leading coefficients after the sign normalization.
    let shared_degree_law = (seeds_ok
        && matches!(coordinates, Coordinates::Identity)
        && local_degree_law
            .iter()
            .any(|h| h.status != LiftHypothesisStatus::Verified))
    .then(|| {
        shared_degree_law_hypothesis(
            "degree_law",
            &format!(
                "deg P_n - deg P_(n-1) is 0 or 1 and P_n, P_(n-2) have positive leading \
                 coefficients (after the sign normalization) for every n >= {n_first}"
            ),
            (rule, first_index, seeds),
            (previous_index, None),
            // Full growth from the seeds (used by the root windows) needs
            // every step to be one.
            |step| step == 1 || (step == 0 && !expects_full_growth),
            row_sign,
        )
    })
    .flatten();
    match shared_degree_law {
        Some(hypothesis) => hypotheses.push(hypothesis),
        None => hypotheses.extend(local_degree_law),
    }

    // The sign condition.
    let (sign_hypothesis, lag_sign) =
        multiplier_sign_hypothesis("lag_coefficient_nonpositive", "B_n", c2, f, n_first);
    let (derivative_sign_hypothesis, lag_sign) = match cdn {
        None => (None, lag_sign),
        Some(c) => {
            let (hypothesis, sign) = multiplier_sign_hypothesis(
                "derivative_coefficient_nonpositive",
                "V_n",
                c,
                f,
                n_first,
            );
            let combined = match (lag_sign, sign) {
                (Some(LagSign::Everywhere), Some(LagSign::Everywhere)) => Some(LagSign::Everywhere),
                (Some(_), Some(_)) => Some(LagSign::NonpositiveAxis),
                _ => None,
            };
            (Some(hypothesis), combined)
        }
    };
    if cdn.is_some() {
        report.reference = format!("{REFERENCE}{THEOREM_2_3}");
    }
    if lag_sign == Some(LagSign::Everywhere) && is_negative_square(c2) {
        report.reference = format!("{REFERENCE}{PROPOSITION_3_3}");
    }
    // Root location, or a root window when neither sign case closes.
    let (location_hypotheses, location) = root_location_hypotheses(
        &seed_previous,
        &seed_current,
        &c1n,
        cdn,
        c2,
        f,
        n_first,
        generated.as_deref(),
        &row,
    );
    let full_growth = degree_difference == Some(1) && (alpha == 1 || beta == 2);
    let core = lag_sign == Some(LagSign::Everywhere) || (lag_sign.is_some() && location.is_some());
    let window = if core || cdn.is_some() {
        None
    } else {
        proved_window(
            &c1n,
            c2,
            f,
            n_first,
            (&seed_previous, &seed_current),
            full_growth,
            &row,
        )
    };
    let stripped_region = match &window {
        Some((window, window_hypotheses)) => {
            hypotheses.extend(window_hypotheses.iter().cloned());
            Some(window.interval())
        }
        None => {
            hypotheses.push(sign_hypothesis);
            hypotheses.extend(derivative_sign_hypothesis);
            match (lag_sign, location) {
                (_, Some(_)) => hypotheses.extend(location_hypotheses),
                (Some(LagSign::Everywhere), None) => {}
                (_, None) => hypotheses.push(unresolved_hypothesis(
                    "root_location",
                    &format!(
                        "every root of every row {variable}_n, n >= {previous_index}, is <= 0 \
                         (needed for B_n <= 0 on x <= 0 only), or the rows stay in a root \
                         window"
                    ),
                    location_hypotheses
                        .iter()
                        .filter(|h| h.status != LiftHypothesisStatus::Verified)
                        .map(|h| format!("{}: {}", h.name, h.evidence))
                        .collect::<Vec<_>>()
                        .join("; "),
                )),
            }
            location.map(|location| RealInterval {
                lower: None,
                lower_closed: false,
                upper: Some(Q::zero()),
                upper_closed: location == RootLocation::NonnegativeCoefficients,
            })
        }
    };
    // The root 0 of the divided power of x, then the coordinate change.
    let transformed_region = stripped_region.map(|region| {
        if power > 0 {
            hull_with_point(region, &Q::zero())
        } else {
            region
        }
    });
    let region = coordinates.region_of_rows(transformed_region);

    // Earlier rows (in the coordinates of P).
    let earlier_rows = &seeds[..seeds.len() - 2];
    let zero_rows = match &region {
        None => {
            hypotheses.extend(earlier_seeds_hypothesis(earlier_rows, first_index, None));
            Vec::new()
        }
        Some(region) => {
            let (hypothesis, zero_rows) = earlier_rows_in_region(earlier_rows, first_index, region);
            hypotheses.extend(hypothesis);
            zero_rows
        }
    };
    let identity = *coordinates == Coordinates::Identity;
    if let Some(requested) = &options.root_region {
        hypotheses.push(requested_region_hypothesis(
            requested,
            region.as_ref(),
            power > 0 && identity,
            current_index,
        ));
    }
    report.hypotheses = hypotheses;
    // No RealRooted route mirrors the derivative-lag form with weak
    // interlacing; none is recorded then.
    let (route, declarations) = if cdn.is_some() {
        (None, Vec::new())
    } else {
        lean_correspondence(
            identity && full_growth,
            lag_sign,
            location,
            window.as_ref().map(|(window, _)| window),
        )
    };
    report.lean_route = route.map(str::to_string);
    report.lean_declarations = declarations;

    let root_region = region.as_ref().map(RootInterval::from_real);
    let conclusion_first = zero_rows.last().map_or(first_index, |last| last + 1);
    let zero_text = if zero_rows.is_empty() {
        String::new()
    } else {
        let list: Vec<String> = zero_rows.iter().map(|j| format!("P_{j}")).collect();
        format!(
            "{} {} zero; ",
            list.join(", "),
            if zero_rows.len() == 1 { "is" } else { "are" }
        )
    };
    let interlacing_text = if identity {
        format!(
            "; for n >= {current_index}, P_(n-1) interlaces P_n (weakly, in the sense of Liu and \
             Wang), and deg P_n - deg P_(n-1) is 0 or 1 for n >= {n_first}"
        )
    } else {
        format!(
            "; for n >= {current_index}, Q_(n-1) interlaces Q_n (weakly, in the sense of Liu and \
             Wang; no interlacing of the P_n is claimed)"
        )
    };
    let statement = format!(
        "{zero_text}every {}P_n, n >= {first_index}, is nonzero and \
         real-rooted{}{interlacing_text}{}",
        if zero_rows.is_empty() { "" } else { "other " },
        root_region.as_ref().map_or(String::new(), |r| format!(
            " with every root in {r}{}",
            if region.as_ref().is_some_and(within_nonpositive_ray) {
                " (so its coefficients have one sign)"
            } else {
                ""
            }
        )),
        if power > 0 && identity {
            format!("; P_n = x^{power} R_n for n >= {previous_index} with R_n(0) != 0")
        } else {
            String::new()
        }
    );
    report.finish(RealRootednessTheoremConclusion {
        first_index: conclusion_first,
        zero_residue_classes: None,
        interlacing_from: identity.then_some(current_index),
        interlacing: identity.then_some(RowInterlacing::Weak),
        simple_roots: false,
        root_region,
        statement,
    })
}

/// A conclusion that is never kept (a report whose hypotheses fail).
fn placeholder_conclusion(first_index: usize) -> RealRootednessTheoremConclusion {
    RealRootednessTheoremConclusion {
        first_index,
        zero_residue_classes: None,
        interlacing_from: None,
        interlacing: None,
        simple_roots: false,
        root_region: None,
        statement: String::new(),
    }
}

fn negate_bivariate(c: &Bivariate) -> Bivariate {
    c.iter()
        .map(|row| row.iter().map(|v| -v).collect())
        .collect()
}

/// The column `[x^j] c / f` as a function of `n`.
fn column(c: &Bivariate, j: usize, f: &[Q]) -> IndexFunction {
    IndexFunction::new(x_column(c, j), f.to_vec())
}

/// The degree-law hypotheses (module documentation, "The degree law").
fn degree_law_hypotheses(
    alpha: usize,
    beta: usize,
    seed_degrees: (usize, usize),
    (c1, cd, c2, f): (&Bivariate, Option<&Bivariate>, &Bivariate, &[Q]),
    n_first: usize,
    rows: Option<&[Vec<Q>]>,
) -> Vec<LiftHypothesis> {
    let delta = seed_degrees.1 - seed_degrees.0;
    // (next state, the power k of x with u_n = [x^k] A_n + d [x^(k+1)] V_n
    // (d = deg P_(n-1)), the power of x of v_n = [x^j] B_n)
    let table =
        |from: usize| -> Result<(usize, Option<usize>, Option<usize>), String> {
            Ok(match (alpha, beta, from) {
                (1, 2, 0) | (0, 2, 0) => return Err(
                    "deg_x c_2 = 2 and the previous rows have equal degrees: the degree may grow \
                     by two"
                        .into(),
                ),
                (1, 2, _) => (1, Some(1), Some(2)),
                (1, 1, 0) => (1, Some(1), Some(1)),
                (1, _, _) => (1, Some(1), None),
                (0, 0, 0) => (0, Some(0), Some(0)),
                (0, 0, _) => (0, Some(0), None),
                (0, 1, 0) => (1, None, Some(1)),
                (0, 1, _) => (0, Some(0), Some(1)),
                (0, 2, _) => (1, None, Some(2)),
                _ => return Err("unsupported coefficient degrees".into()),
            })
        };
    let mut transitions: Vec<(usize, usize, Option<usize>, Option<usize>)> = Vec::new();
    let mut state = delta;
    while !transitions.iter().any(|t| t.0 == state) {
        match table(state) {
            Ok((to, k, j)) => {
                transitions.push((state, to, k, j));
                state = to;
            }
            Err(reason) => {
                return vec![unresolved_hypothesis(
                    "degree_law",
                    &format!("deg P_n - deg P_(n-1) is 0 or 1 for every n >= {n_first}"),
                    reason,
                )];
            }
        }
    }
    let single_state = transitions.len() == 1 && transitions[0].0 == transitions[0].1;
    let mut extra = Vec::new();
    let steps: Vec<DegreeStep> = transitions
        .iter()
        .map(|&(from, to, k, j)| {
            let (u, u_text) = match k {
                None => (IndexFunction::constant(Q::zero()), "0".to_string()),
                Some(k) => {
                    let derivative_column = cd.map(|c| x_column(c, k + 1)).unwrap_or_default();
                    if derivative_column.is_empty() {
                        (column(c1, k, f), format!("[x^{k}] A_n"))
                    } else {
                        // deg P_(n-1) = d_0 + to (n - N) on one state, else at
                        // least d_0 (with [x^(k+1)] V_n >= 0).
                        let d0 = Q::from_integer(BigInt::from(seed_degrees.1));
                        let degree: Vec<Q> = if single_state && to == 1 {
                            vec![d0 - Q::from_integer(BigInt::from(n_first)), Q::one()]
                        } else {
                            vec![d0]
                        };
                        if !single_state {
                            extra.push(sign_hypothesis(
                                &format!("derivative_degree_coefficient_nonnegative_{k}"),
                                &format!(
                                    "[x^{}] V_n >= 0 for every n >= {n_first} (so that \
                                     [x^{k}] A_n + deg P_(n-1) [x^{}] V_n is at least its \
                                     value at deg P_(n-1) = {})",
                                    k + 1,
                                    k + 1,
                                    seed_degrees.1
                                ),
                                &format!("[x^{}] V_n", k + 1),
                                &column(cd.expect("a derivative term"), k + 1, f),
                                n_first,
                                IndexSignCondition::Nonnegative,
                            ));
                        }
                        let numerator = vec_poly::add(
                            &vec_poly::trim(&x_column(c1, k)),
                            &vec_poly::mul(&degree, &derivative_column),
                        );
                        (
                            IndexFunction::new(numerator, f.to_vec()),
                            format!(
                                "[x^{k}] A_n + {} [x^{}] V_n",
                                if single_state && to == 1 {
                                    "deg P_(n-1)"
                                } else {
                                    "(a lower bound for deg P_(n-1))"
                                },
                                k + 1
                            ),
                        )
                    }
                }
            };
            let (v, v_text) = match j {
                None => (IndexFunction::constant(Q::zero()), "0".to_string()),
                Some(j) => (column(c2, j, f), format!("[x^{j}] B_n")),
            };
            DegreeStep {
                from,
                to,
                u,
                v,
                u_text,
                v_text,
            }
        })
        .collect();
    let statement = format!(
        "deg P_n - deg P_(n-1) is 0 or 1 and P_n, P_(n-2) have positive leading coefficients \
         (after the sign normalization) for every n >= {n_first}"
    );
    let schedule = match (steps.len(), single_state) {
        (1, true) if steps[0].to == 1 => "full growth (deg P_n = deg P_(n-1) + 1)".to_string(),
        (1, true) => "constant degree".to_string(),
        _ => format!(
            "degree steps (deg P_(n-1) - deg P_(n-2) -> deg P_n - deg P_(n-1)) {}",
            steps
                .iter()
                .map(|s| format!("{} -> {}", s.from, s.to))
                .collect::<Vec<_>>()
                .join(", ")
        ),
    };
    // The leading coefficients l_k (the coefficient of x^(d_k) for the
    // degrees d_k of the schedule) of the computed rows, k >= N - 2.
    let leads: Vec<Q> = match rows {
        Some(rows) => {
            let mut degrees = vec![seed_degrees.0, seed_degrees.1];
            let mut state = delta;
            while degrees.len() < rows.len() {
                state = steps
                    .iter()
                    .find(|s| s.from == state)
                    .map_or(state, |s| s.to);
                degrees.push(degrees[degrees.len() - 1] + state);
            }
            rows.iter()
                .zip(&degrees)
                .map(|(row, &degree)| vec_poly::coeff(row, degree))
                .collect()
        }
        None => Vec::new(),
    };
    steps
        .iter()
        .map(|step| {
            let name = if steps.len() == 1 {
                "degree_law".to_string()
            } else {
                format!("degree_law_from_step_{}", step.from)
            };
            let label = format!(
                "l_n = ({}) l_(n-1) + ({}) l_(n-2) when deg P_(n-1) - deg P_(n-2) = {} (then deg \
                 P_n - deg P_(n-1) = {})",
                step.u_text, step.v_text, step.from, step.to
            );
            let mut hypothesis = positive_sequence_hypothesis(
                &name,
                &format!("{statement}: {schedule}; {label}"),
                "l_n",
                &step.u,
                &step.v,
                n_first,
                &leads,
                single_state,
                false,
            );
            hypothesis.evidence = format!("{label}: {}", hypothesis.evidence);
            hypothesis
        })
        .chain(extra)
        .collect()
}

/// `u_k > 0` for every `k >= first - 2` for the sequence
/// `u_n = s(n) u_(n-1) + t(n) u_(n-2)` (`n >= first`), with positive
/// `u_(first-2)`, `u_(first-1)` and `values = u_(first-2), u_(first-1), ...`
/// computed exactly (possibly empty): by `s, t >= 0`, `s + t > 0`, or (when
/// `allow_ratio`) a positive-ratio induction (module documentation, "The
/// degree law").  A computed value `<= 0` refutes the statement when
/// `refutes` (the values are the sequence itself), and otherwise only ends
/// the attempt.
#[allow(clippy::too_many_arguments)]
fn positive_sequence_hypothesis(
    name: &str,
    statement: &str,
    sequence: &str,
    s: &IndexFunction,
    t: &IndexFunction,
    first: usize,
    values: &[Q],
    allow_ratio: bool,
    refutes: bool,
) -> LiftHypothesis {
    if let Some(offset) = values.iter().position(|v| !v.is_positive()) {
        let index = first - 2 + offset;
        let evidence = format!(
            "{sequence} = {} at n = {index} (exact)",
            format_rational_coeff(&values[offset])
        );
        return if refutes {
            exact_hypothesis(name, statement, Err((Some(index), evidence)))
        } else {
            unresolved_hypothesis(
                name,
                statement,
                format!("{evidence}, so the sufficient test fails (not a disproof)"),
            )
        };
    }
    let s_decision = s.decide(first, IndexSignCondition::Nonnegative);
    let t_decision = t.decide(first, IndexSignCondition::Nonnegative);
    let sum_decision = s.add(t).decide(first, IndexSignCondition::Positive);
    if s_decision.holds() && t_decision.holds() && sum_decision.holds() {
        return exact_hypothesis(
            name,
            statement,
            Ok(format!(
                "positive seeds and nonnegative step coefficients: {}; {}; {}",
                s_decision.describe("s_n"),
                t_decision.describe("t_n"),
                sum_decision.describe("s_n + t_n")
            )),
        );
    }
    let direct = format!(
        "{}; {}; {}",
        s_decision.describe("s_n"),
        t_decision.describe("t_n"),
        sum_decision.describe("s_n + t_n")
    );
    if !allow_ratio {
        return unresolved_hypothesis(
            name,
            statement,
            format!(
                "the step coefficients are not all nonnegative (not a disproof): {direct} (s_n = \
                 {}, t_n = {})",
                s.text(),
                t.text()
            ),
        );
    }
    let ratios: Vec<(usize, Q)> = values
        .windows(2)
        .enumerate()
        .skip(1)
        .map(|(offset, pair)| (first - 1 + offset, &pair[1] / &pair[0]))
        .collect();
    let Some(first_ratio) = ratios.first().map(|(_, r)| r.clone()) else {
        return unresolved_hypothesis(
            name,
            statement,
            format!("no computed value of {sequence} at n = {first} (not a disproof): {direct}"),
        );
    };
    let t_nonpositive = t.decide(first + 1, IndexSignCondition::Nonpositive);
    if !t_nonpositive.holds() {
        return unresolved_hypothesis(
            name,
            statement,
            format!(
                "neither nonnegative step coefficients ({direct}) nor t_n <= 0 for a \
                 positive-ratio \
                 induction ({}) (not a disproof)",
                t_nonpositive.describe("t_n")
            ),
        );
    }
    let mut candidates = Vec::new();
    if let Some(fitted) = fit_index_rational_function(&ratios) {
        candidates.push(("the exact ratio", fitted));
    }
    candidates.push((
        "s_n / 2",
        s.scale(&Q::new(BigInt::from(1), BigInt::from(2))),
    ));
    let mut failures = Vec::new();
    for (label, g) in candidates {
        let positive = g.decide(first, IndexSignCondition::Positive);
        if !positive.holds() {
            failures.push(format!(
                "g(n) = {} ({label}): {}",
                g.text(),
                positive.describe("g")
            ));
            continue;
        }
        let Some(g_first) = g.at(first) else {
            failures.push(format!(
                "g(n) = {} ({label}): undefined at n = {first}",
                g.text()
            ));
            continue;
        };
        if first_ratio < g_first {
            failures.push(format!(
                "g(n) = {} ({label}): u_{first}/u_{} = {} < g({first}) = {}",
                g.text(),
                first - 1,
                format_rational_coeff(&first_ratio),
                format_rational_coeff(&g_first)
            ));
            continue;
        }
        let condition = s.add(&t.div(&g.shifted(-1))).sub(&g);
        let decision = condition.decide(first + 1, IndexSignCondition::Nonnegative);
        if !decision.holds() {
            failures.push(format!(
                "g(n) = {} ({label}): {}",
                g.text(),
                decision.describe("s_n + t_n / g(n-1) - g(n)")
            ));
            continue;
        }
        let mut hypothesis = LiftHypothesis::from_sign_decision(
            name,
            statement,
            decision,
            "s_n + t_n / g(n-1) - g(n)",
        );
        hypothesis.evidence = format!(
            "positive-ratio induction with g(n) = {} ({label}): {}; u_{first}/u_{} = {} >= \
             g({first}) = {} exactly; {}; and {}",
            g.text(),
            positive.describe("g(n)"),
            first - 1,
            format_rational_coeff(&first_ratio),
            format_rational_coeff(&g_first),
            t_nonpositive.describe("t_n"),
            hypothesis.evidence
        );
        return hypothesis;
    }
    unresolved_hypothesis(
        name,
        statement,
        format!(
            "no candidate closes the positive-ratio induction (not a disproof): {}",
            failures.join("; ")
        ),
    )
}

/// "`B_n(x) <= 0` for every real `x` and every `n >= first`" (module
/// documentation, "The sign condition").
fn everywhere_nonpositive_hypothesis(
    label: &str,
    c2: &Bivariate,
    f: &[Q],
    first: usize,
) -> LiftHypothesis {
    let name = "lag_coefficient_nonpositive";
    let statement = format!("{label}(x) <= 0 for every real x and every n >= {first}");
    match x_degree(c2) {
        Some(0) => sign_hypothesis(
            name,
            &statement,
            label,
            &column(c2, 0, f),
            first,
            IndexSignCondition::Nonpositive,
        ),
        Some(2) => {
            let (p0, p2) = (column(c2, 0, f), column(c2, 2, f));
            let discriminant = IndexFunction::new(
                vec_poly::sub(
                    &vec_poly::mul(&x_column(c2, 1), &x_column(c2, 1)),
                    &vec_poly::scale(
                        &vec_poly::mul(&x_column(c2, 2), &x_column(c2, 0)),
                        &Q::from_integer(BigInt::from(4)),
                    ),
                ),
                vec_poly::mul(f, f),
            );
            let decisions = [
                (
                    format!("p_2 = [x^2] {label}"),
                    p2.decide(first, IndexSignCondition::Nonpositive),
                ),
                (
                    format!("p_0 = [x^0] {label}"),
                    p0.decide(first, IndexSignCondition::Nonpositive),
                ),
                (
                    "p_1^2 - 4 p_2 p_0".to_string(),
                    discriminant.decide(first, IndexSignCondition::Nonpositive),
                ),
            ];
            let evidence: Vec<String> = decisions
                .iter()
                .map(|(label, decision)| decision.describe(label))
                .collect();
            if decisions.iter().all(|(_, decision)| decision.holds()) {
                exact_hypothesis(
                    name,
                    &statement,
                    Ok(format!(
                        "{label} = p_2 x^2 + p_1 x + p_0 with p_2 <= 0, p_0 <= 0 and p_1^2 <= \
                         4 p_2 \
                         p_0 (so {label} <= 0 on the real line): {}",
                        evidence.join("; ")
                    )),
                )
            } else {
                unresolved_hypothesis(
                    name,
                    &statement,
                    format!(
                        "the quadratic test p_2 <= 0, p_0 <= 0, p_1^2 <= 4 p_2 p_0 does not hold \
                         for every n (not a disproof): {}",
                        evidence.join("; ")
                    ),
                )
            }
        }
        _ => unresolved_hypothesis(
            name,
            &statement,
            format!("{label} has degree 1 in x, so it is not <= 0 on the whole line"),
        ),
    }
}

/// The sign hypothesis "`c / f <= 0` on the real line, or on `x <= 0`" for
/// the multiplier `label` (`B_n` or `V_n`), and the case that holds.
fn multiplier_sign_hypothesis(
    name: &str,
    label: &str,
    c: &Bivariate,
    f: &[Q],
    n_first: usize,
) -> (LiftHypothesis, Option<LagSign>) {
    let mut everywhere = everywhere_nonpositive_hypothesis(label, c, f, n_first);
    if everywhere.status == LiftHypothesisStatus::Verified {
        everywhere.name = name.into();
        everywhere.statement = format!("{label}(x) <= 0 for every real x and every n >= {n_first}");
        return (everywhere, Some(LagSign::Everywhere));
    }
    let mut ray = negative_axis_sign_hypothesis(
        name,
        &format!(
            "{label}(x) <= 0 for every real x, or for every x <= 0 (with a root location), for \
             every n >= {n_first}"
        ),
        label,
        c,
        f,
        n_first,
        false,
    );
    if ray.status == LiftHypothesisStatus::Verified {
        ray.evidence = format!("on x <= 0 (x < 0 and continuity): {}", ray.evidence);
        (ray, Some(LagSign::NonpositiveAxis))
    } else {
        ray.evidence = format!(
            "neither everywhere ({}) nor on x <= 0 ({})",
            everywhere.evidence, ray.evidence
        );
        (ray, None)
    }
}

/// `factor c(n, x) / x` for a `c` divisible by `x`.
fn scale_divided_by_x(c: &Bivariate, factor: usize) -> Bivariate {
    let factor = Q::from_integer(BigInt::from(factor));
    c.iter()
        .map(|row| row.iter().skip(1).map(|value| value * &factor).collect())
        .collect()
}

/// Whether `c_2 = -kappa(n) (x + v)^2`-type: degree two with identically
/// vanishing discriminant.
fn is_negative_square(c2: &Bivariate) -> bool {
    x_degree(c2) == Some(2) && {
        let discriminant = vec_poly::sub(
            &vec_poly::mul(&x_column(c2, 1), &x_column(c2, 1)),
            &vec_poly::scale(
                &vec_poly::mul(&x_column(c2, 2), &x_column(c2, 0)),
                &Q::from_integer(BigInt::from(4)),
            ),
        );
        vec_poly::is_zero(&vec_poly::trim(&discriminant))
    }
}

/// The root-location hypotheses (module documentation, "Root location"):
/// the first method whose hypotheses all hold, else the hypotheses of the
/// root-sign induction (for the evidence).
#[allow(clippy::too_many_arguments)]
fn root_location_hypotheses(
    seed_previous: &[Q],
    seed_current: &[Q],
    c1: &Bivariate,
    cd: Option<&Bivariate>,
    c2: &Bivariate,
    f: &[Q],
    n_first: usize,
    generated: Option<&[Vec<Q>]>,
    row: &dyn Fn(usize) -> String,
) -> (Vec<LiftHypothesis>, Option<RootLocation>) {
    let (previous_index, current_index) = (n_first - 2, n_first - 1);
    // Root-sign induction.
    let negative = has_real_roots_below(seed_previous, &Q::zero())
        && has_real_roots_below(seed_current, &Q::zero());
    let seed_hypothesis = exact_hypothesis(
        "seed_roots_negative",
        &format!(
            "{} and {} have every root < 0",
            row(previous_index),
            row(current_index)
        ),
        if negative {
            Ok("checked exactly".into())
        } else {
            Err((
                Some(current_index),
                "a seed of the induction has a root >= 0 (exact)".into(),
            ))
        },
    );
    let constants: Vec<Q> = match generated {
        Some(rows) => rows
            .iter()
            .map(|r| r.first().cloned().unwrap_or_else(Q::zero))
            .collect(),
        None => vec![
            seed_previous.first().cloned().unwrap_or_else(Q::zero),
            seed_current.first().cloned().unwrap_or_else(Q::zero),
        ],
    };
    let constant_hypothesis = positive_sequence_hypothesis(
        "constant_terms_positive",
        &format!(
            "P_n(0) > 0 for every n >= {n_first} (normalized rows), with P_n(0) = A_n(0) \
             P_(n-1)(0) + B_n(0) P_(n-2)(0)"
        ),
        "P_n(0)",
        &column(c1, 0, f),
        &column(c2, 0, f),
        n_first,
        &constants,
        true,
        true,
    );
    let mut root_sign = vec![seed_hypothesis, constant_hypothesis];
    // With a derivative term, P_n(0) = A_n(0) P_(n-1)(0) + B_n(0) P_(n-2)(0)
    // needs V_n(0) = 0.
    if let Some(c) = cd {
        root_sign.push(exact_hypothesis(
            "derivative_coefficient_vanishes_at_zero",
            &format!("V_n(0) = 0 for every n >= {n_first}"),
            if x_column(c, 0).is_empty() {
                Ok("[x^0] c_d is the zero polynomial".into())
            } else {
                Err((None, "[x^0] c_d is not the zero polynomial".into()))
            },
        ));
    }
    if root_sign
        .iter()
        .all(|h| h.status == LiftHypothesisStatus::Verified)
    {
        return (root_sign, Some(RootLocation::Negative));
    }
    // Nonnegative coefficients.
    let seeds_nonnegative = seed_previous
        .iter()
        .chain(seed_current)
        .all(|c| !c.is_negative());
    let mut evidence = Vec::new();
    let mut multipliers_nonnegative = true;
    let mut multipliers = vec![("A_n", c1), ("B_n", c2)];
    if let Some(c) = cd {
        multipliers.push(("V_n", c));
    }
    for (label, c) in multipliers {
        for j in 0..x_width(c) {
            if x_column(c, j).is_empty() {
                continue;
            }
            let decision = column(c, j, f).decide(n_first, IndexSignCondition::Nonnegative);
            evidence.push(decision.describe(&format!("[x^{j}] {label}")));
            multipliers_nonnegative &= decision.holds();
        }
    }
    if seeds_nonnegative && multipliers_nonnegative {
        let hypothesis = exact_hypothesis(
            "nonnegative_coefficients",
            &format!(
                "{} and {} and every coefficient [x^k] A_n, [x^k] B_n (and [x^k] V_n) for n >= \
                 {n_first} are nonnegative, so every row has nonnegative coefficients and (being \
                 real-rooted) every root <= 0",
                row(previous_index),
                row(current_index)
            ),
            Ok(format!("seeds checked exactly; {}", evidence.join("; "))),
        );
        return (
            vec![hypothesis],
            Some(RootLocation::NonnegativeCoefficients),
        );
    }
    let mut hypotheses = root_sign;
    hypotheses.push(unresolved_hypothesis(
        "nonnegative_coefficients",
        "the seeds and the multipliers A_n, B_n have nonnegative coefficients",
        format!(
            "{}{}",
            if seeds_nonnegative {
                ""
            } else {
                "a seed has a negative coefficient; "
            },
            evidence.join("; ")
        ),
    ));
    (hypotheses, None)
}

/// The requested region must contain the certified one.
fn requested_region_hypothesis(
    requested: &RealInterval,
    certified: Option<&RealInterval>,
    zero_is_root: bool,
    current_index: usize,
) -> LiftHypothesis {
    let name = "requested_root_region";
    let Some(certified) = certified else {
        return unresolved_hypothesis(
            name,
            &format!(
                "the requested region {} contains the certified region",
                open_root_interval(requested)
            ),
            "no root region is certified (B_n <= 0 on the whole line; the root location does not \
             hold or is not shown)"
                .into(),
        );
    };
    let certified_text = RootInterval::from_real(certified).to_string();
    let statement = format!(
        "the requested region {} contains {certified_text}, the region this theorem certifies",
        open_root_interval(requested)
    );
    // The requested region is open at its finite ends.
    let upper_ok = requested
        .upper
        .as_ref()
        .is_none_or(|requested_upper| match &certified.upper {
            None => false,
            Some(upper) => {
                upper < requested_upper || (upper == requested_upper && !certified.upper_closed)
            }
        });
    let lower_ok = requested
        .lower
        .as_ref()
        .is_none_or(|requested_lower| match &certified.lower {
            None => false,
            Some(lower) => {
                lower > requested_lower || (lower == requested_lower && !certified.lower_closed)
            }
        });
    if lower_ok && upper_ok {
        exact_hypothesis(name, &statement, Ok("exact".into()))
    } else if zero_is_root && requested.upper.as_ref().is_some_and(Zero::is_zero) {
        exact_hypothesis(
            name,
            &statement,
            Err((
                Some(current_index),
                format!("0 is a root of P_{current_index}, outside the open region (exact)"),
            )),
        )
    } else {
        unresolved_hypothesis(
            name,
            &statement,
            format!("this theorem certifies the region {certified_text} only"),
        )
    }
}

/// The RealRooted route and declarations mirrored by the check.
fn lean_correspondence(
    full_growth: bool,
    lag_sign: Option<LagSign>,
    location: Option<RootLocation>,
    window: Option<&RootWindow>,
) -> (Option<&'static str>, Vec<String>) {
    let (route, names): (Option<&str>, Vec<&str>) = if !full_growth {
        (
            None,
            vec!["RealRooted.MaWang.strictInterl_of_interlaces_evalCoeff_nonpos"],
        )
    } else if let Some(window) = window {
        match (&window.lower, &window.upper) {
            (Some(_), Some(_)) => (
                Some(LEAN_ROUTE_WINDOW),
                vec!["RealRooted.threeTerm_interlaces_of_roots_mem_Icc"],
            ),
            (None, Some(_)) => (
                Some(LEAN_ROUTE_WINDOW),
                vec!["RealRooted.threeTerm_interlaces_of_roots_le"],
            ),
            _ => (None, vec!["RealRooted.threeTerm_interlaces_of_window"]),
        }
    } else {
        match (lag_sign, location) {
            (Some(LagSign::Everywhere), _) => (
                Some(LEAN_ROUTE_SIGN),
                vec![
                    "RealRooted.threeTerm_interlaces_of_eval_nonpos",
                    "RealRooted.threeTerm_interlaces",
                    "RealRooted.eval_nonpos_of_natDegree_le_two",
                ],
            ),
            (_, Some(RootLocation::NonnegativeCoefficients)) => (
                Some(LEAN_ROUTE_SIGN),
                vec![
                    "RealRooted.threeTerm_interlaces_of_nonnegCoeffs",
                    "RealRooted.threeTerm_interlaces",
                    "RealRooted.threeTerm_hasNonnegCoeffs",
                ],
            ),
            _ => (
                Some(LEAN_ROUTE_SIGN),
                vec![
                    "RealRooted.threeTerm_interlaces_of_nonnegCoeffs",
                    "RealRooted.threeTerm_interlaces",
                ],
            ),
        }
    };
    (route, names.into_iter().map(str::to_string).collect())
}

// ---------------------------------------------------------------------------
// Root windows
// ---------------------------------------------------------------------------

/// A closed root window `W` with at least one finite end.
#[derive(Debug, Clone, PartialEq, Eq)]
struct RootWindow {
    lower: Option<Q>,
    upper: Option<Q>,
}

impl RootWindow {
    fn interval(&self) -> RealInterval {
        RealInterval {
            lower: self.lower.clone(),
            lower_closed: self.lower.is_some(),
            upper: self.upper.clone(),
            upper_closed: self.upper.is_some(),
        }
    }

    fn text(&self) -> String {
        RootInterval::from_real(&self.interval()).to_string()
    }
}

/// The candidate windows: the rays and gaps of the rational roots of `w` for
/// a separable `c_2 = gamma(n) w(x)`.
fn window_candidates(c2: &Bivariate) -> Vec<RootWindow> {
    let Some((_, w)) = separate(c2) else {
        return Vec::new();
    };
    let Ok(factorization) = crate::factorization::irreducible_factorization_rational_coeffs(&w)
    else {
        return Vec::new();
    };
    let mut ends: Vec<Q> = factorization
        .factors
        .iter()
        .filter(|(factor, _)| factor.len() == 2)
        .map(|(factor, _)| -Q::new(factor[0].clone(), factor[1].clone()))
        .collect();
    ends.sort();
    ends.dedup();
    let mut windows = Vec::new();
    for end in &ends {
        windows.push(RootWindow {
            lower: None,
            upper: Some(end.clone()),
        });
        windows.push(RootWindow {
            lower: Some(end.clone()),
            upper: None,
        });
    }
    for pair in ends.windows(2) {
        windows.push(RootWindow {
            lower: Some(pair[0].clone()),
            upper: Some(pair[1].clone()),
        });
    }
    windows
}

/// `c(n, a + b t)`, as a bivariate polynomial in `n` and `t`.
fn substitute_bivariate(c: &Bivariate, a: &Q, b: &Q) -> Bivariate {
    c.iter()
        .map(|row| {
            let row = vec_poly::trim(row);
            if vec_poly::is_zero(&row) {
                Vec::new()
            } else {
                let mut substituted = vec_poly::compose(&row, &[a.clone(), b.clone()]);
                vec_poly::strip_trailing_zeros(&mut substituted);
                substituted
            }
        })
        .collect()
}

/// `c(n, value) / f(n)` as a function of `n`.
fn value_at(c: &Bivariate, value: &Q, f: &[Q]) -> IndexFunction {
    let column: Vec<Q> = c
        .iter()
        .map(|row| vec_poly::evaluate(&vec_poly::trim(row), value))
        .collect();
    IndexFunction::new(column, f.to_vec())
}

/// Whether every `t`-coefficient of `c / f` satisfies `condition` for every
/// `n >= first`, with the decisions as evidence (stopping at the first
/// failure).
fn all_columns_decided(
    c: &Bivariate,
    f: &[Q],
    first: usize,
    condition: IndexSignCondition,
    label: &str,
    evidence: &mut Vec<String>,
) -> bool {
    for j in 0..x_width(c) {
        let column = x_column(c, j);
        if column.is_empty() {
            continue;
        }
        let decision = IndexFunction::new(column, f.to_vec()).decide(first, condition);
        evidence.push(decision.describe(&format!("[t^{j}] {label}")));
        if !decision.holds() {
            return false;
        }
    }
    true
}

/// One sign condition with its evidence.
fn decided(
    function: &IndexFunction,
    first: usize,
    condition: IndexSignCondition,
    label: &str,
    evidence: &mut Vec<String>,
) -> bool {
    let decision = function.decide(first, condition);
    evidence.push(decision.describe(label));
    decision.holds()
}

/// Whether the nonzero real-rooted `p` has every root in the closed window.
fn has_roots_in_window(p: &[Q], window: &RootWindow) -> bool {
    let p = vec_poly::trim(p);
    if vec_poly::is_zero(&p) || !is_real_rooted_rational(&p) {
        return false;
    }
    let one_signed =
        |q: &[Q]| q.iter().all(|c| !c.is_negative()) || q.iter().all(|c| !c.is_positive());
    let upper_ok = window
        .upper
        .as_ref()
        .is_none_or(|u| one_signed(&vec_poly::compose(&p, &[u.clone(), Q::one()])));
    let lower_ok = window
        .lower
        .as_ref()
        .is_none_or(|l| one_signed(&vec_poly::compose(&p, &[l.clone(), -Q::one()])));
    upper_ok && lower_ok
}

/// The hypotheses of the window version for `window` (module documentation,
/// "Root windows").
#[allow(clippy::too_many_arguments)]
fn window_hypotheses(
    window: &RootWindow,
    c1: &Bivariate,
    c2: &Bivariate,
    f: &[Q],
    n_first: usize,
    seeds: (&[Q], &[Q]),
    full_growth: bool,
    row: &dyn Fn(usize) -> String,
) -> Vec<LiftHypothesis> {
    let w = window.text();
    let mut hypotheses = Vec::new();
    // B_n <= 0 on W.
    let statement = format!("B_n(x) <= 0 for every x in {w} and every n >= {n_first}");
    let mut evidence = Vec::new();
    let one = Q::one();
    let below = window.upper.as_ref().map(|u| {
        let mut e = Vec::new();
        let holds = all_columns_decided(
            &substitute_bivariate(c2, u, &-one.clone()),
            f,
            n_first,
            IndexSignCondition::Nonpositive,
            &format!("B_n({} - t)", format_rational_coeff(u)),
            &mut e,
        );
        (holds, e)
    });
    let above = window.lower.as_ref().map(|l| {
        let mut e = Vec::new();
        let holds = all_columns_decided(
            &substitute_bivariate(c2, l, &one),
            f,
            n_first,
            IndexSignCondition::Nonpositive,
            &format!("B_n({} + t)", format_rational_coeff(l)),
            &mut e,
        );
        (holds, e)
    });
    let mut sign_ok = false;
    for (holds, e) in [below, above].into_iter().flatten() {
        evidence.extend(e);
        if holds {
            sign_ok = true;
            break;
        }
    }
    if !sign_ok {
        if let (Some(l), Some(u)) = (&window.lower, &window.upper) {
            // Convex B_n with B_n(L), B_n(U) <= 0.
            let mut e = Vec::new();
            sign_ok = x_degree(c2).is_some_and(|d| d <= 2)
                && decided(
                    &column(c2, 2, f),
                    n_first,
                    IndexSignCondition::Nonnegative,
                    "[x^2] B_n",
                    &mut e,
                )
                && decided(
                    &value_at(c2, l, f),
                    n_first,
                    IndexSignCondition::Nonpositive,
                    &format!("B_n({})", format_rational_coeff(l)),
                    &mut e,
                )
                && decided(
                    &value_at(c2, u, f),
                    n_first,
                    IndexSignCondition::Nonpositive,
                    &format!("B_n({})", format_rational_coeff(u)),
                    &mut e,
                );
            if sign_ok {
                evidence = vec![format!(
                    "B_n is convex in x with B_n <= 0 at both ends: {}",
                    e.join("; ")
                )];
            } else {
                evidence.extend(e);
            }
        }
    }
    hypotheses.push(if sign_ok {
        exact_hypothesis(
            "lag_coefficient_nonpositive",
            &statement,
            Ok(evidence.join("; ")),
        )
    } else {
        unresolved_hypothesis(
            "lag_coefficient_nonpositive",
            &statement,
            format!(
                "the sufficient tests fail (not a disproof): {}",
                evidence.join("; ")
            ),
        )
    });
    // Above U: A_n > 0 and B_n >= 0.
    if let Some(u) = &window.upper {
        let mut e = Vec::new();
        let u_text = format_rational_coeff(u);
        let holds = all_columns_decided(
            &substitute_bivariate(c1, u, &one),
            f,
            n_first,
            IndexSignCondition::Nonnegative,
            &format!("A_n({u_text} + t)"),
            &mut e,
        ) && decided(
            &value_at(c1, &(u + &one), f),
            n_first,
            IndexSignCondition::Positive,
            &format!("A_n({})", format_rational_coeff(&(u + &one))),
            &mut e,
        ) && all_columns_decided(
            &substitute_bivariate(c2, u, &one),
            f,
            n_first,
            IndexSignCondition::Nonnegative,
            &format!("B_n({u_text} + t)"),
            &mut e,
        );
        let statement = format!(
            "A_n(x) > 0 and B_n(x) >= 0 for every x > {u_text} and every n >= {n_first} (no row \
             gets a root > {u_text})"
        );
        hypotheses.push(if holds {
            exact_hypothesis("window_upper_end", &statement, Ok(e.join("; ")))
        } else {
            unresolved_hypothesis(
                "window_upper_end",
                &statement,
                format!("the monomial tests fail (not a disproof): {}", e.join("; ")),
            )
        });
    }
    // Below L: A_n < 0 and B_n >= 0, with full degree growth.
    if let Some(l) = &window.lower {
        let l_text = format_rational_coeff(l);
        let statement = format!(
            "deg P_n = deg P_(n-1) + 1, A_n(x) < 0 and B_n(x) >= 0 for every x < {l_text} and \
             every n >= {n_first} (no row gets a root < {l_text})"
        );
        let mut e = Vec::new();
        let holds = full_growth
            && all_columns_decided(
                &substitute_bivariate(c1, l, &-one.clone()),
                f,
                n_first,
                IndexSignCondition::Nonpositive,
                &format!("A_n({l_text} - t)"),
                &mut e,
            )
            && decided(
                &value_at(c1, &(l - &one), f),
                n_first,
                IndexSignCondition::Negative,
                &format!("A_n({})", format_rational_coeff(&(l - &one))),
                &mut e,
            )
            && all_columns_decided(
                &substitute_bivariate(c2, l, &-one.clone()),
                f,
                n_first,
                IndexSignCondition::Nonnegative,
                &format!("B_n({l_text} - t)"),
                &mut e,
            );
        if !full_growth {
            e.insert(0, "the degree law is not full growth from the seeds".into());
        }
        hypotheses.push(if holds {
            exact_hypothesis("window_lower_end", &statement, Ok(e.join("; ")))
        } else {
            unresolved_hypothesis(
                "window_lower_end",
                &statement,
                format!("the tests fail (not a disproof): {}", e.join("; ")),
            )
        });
    }
    let (previous_index, current_index) = (n_first - 2, n_first - 1);
    let inside = has_roots_in_window(seeds.0, window) && has_roots_in_window(seeds.1, window);
    hypotheses.push(exact_hypothesis(
        "seed_roots_in_window",
        &format!(
            "{} and {} have every root in {w}",
            row(previous_index),
            row(current_index)
        ),
        if inside {
            Ok("checked exactly".into())
        } else {
            Err((
                Some(current_index),
                format!("a seed of the induction has a root outside {w} (exact)"),
            ))
        },
    ));
    hypotheses
}

/// The first candidate window whose hypotheses all hold.
#[allow(clippy::too_many_arguments)]
fn proved_window(
    c1: &Bivariate,
    c2: &Bivariate,
    f: &[Q],
    n_first: usize,
    seeds: (&[Q], &[Q]),
    full_growth: bool,
    row: &dyn Fn(usize) -> String,
) -> Option<(RootWindow, Vec<LiftHypothesis>)> {
    window_candidates(c2).into_iter().find_map(|window| {
        let hypotheses = window_hypotheses(&window, c1, c2, f, n_first, seeds, full_growth, row);
        hypotheses
            .iter()
            .all(|h| h.status == LiftHypothesisStatus::Verified)
            .then_some((window, hypotheses))
    })
}
