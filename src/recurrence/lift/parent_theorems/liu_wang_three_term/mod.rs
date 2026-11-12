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
pub(super) const ROUTE: RoutePolicy = RoutePolicy {
    reference: REFERENCE,
    phase: DiscoveryPhase::Direct,
    after_direct_eligible: None,
    replay: ParentReplayPolicy::RecordedRoute,
    discovery_replay: DiscoveryReplayPolicy::Direct,
    check_at_supplied_seeds: liu_wang_three_term_report,
    start_policy: StartPolicy::MovesStart,
    is_preferred: proof_with_root_region,
    certificate_name: None,
    discovery: DiscoveryPolicy::Always,
    in_residue_class_menu: true,
    in_index_section_menu: true,
    in_degree_law_reversal_menu: true,
    parent: ParentPolicy {
        stage: ParentStage::ReplacesEarly,
        region: ParentRegion::Unrestricted,
        record_only_if_applicable: false,
    },
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

mod analysis;
use analysis::*;
mod hypotheses;
use hypotheses::*;
mod window;
use window::*;
