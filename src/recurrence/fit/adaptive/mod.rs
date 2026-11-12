//! Adaptive recurrence search: options, budgets, diagnostics and the search
//! loop over candidate shapes.
//!
//! # Left-hand side variation
//!
//! A rule `f(n, t) P(n) = sum c(n, t) D^d P(n-r)` is found by solving for the
//! coefficients with one coefficient of the left-hand factor fixed; the
//! historical search fixes the constant term `f(0, 0) = 1` (with `n` the
//! 1-based search index).  That misses every left-hand factor whose constant
//! term is zero, such as `n` for the shifted Legendre polynomials
//! (`n P(n) = (2n - 1)(2t - 1) P(n-1) - (n - 1) P(n-2)` from `n = 1`), and
//! factors `x^e g(n)` whose division puts negative powers of `x` into the
//! coefficients.  A homogeneous relation can always be cleared of such
//! powers by multiplying with `x^e`; the point is a left-hand factor that is
//! allowed to carry them, so the rule stays runnable (each row is the exact
//! quotient by the factor).
//!
//! With [`AdaptiveSearchOptions::try_laurent`] (off by default in the
//! library, on by default in the `polytool recurrence` command and the MCP
//! recurrence tools) the search
//! runs the historical passes unchanged (the given options, then the
//! denominator escalation) and, when they find nothing, a third stage in its
//! own complexity order ([`candidate_complexity`]): the same shapes with
//!
//! - a left-hand monomial `x^e`, `1 <= e <= max_laurent_x_power`
//!   ([`RecurrenceOptions::leader_x_power`]; for a factor `p(n)` this is the
//!   product form `p(n) x^e` with `deg p` unknowns), and
//! - a normalization chart other than the constant term: the coefficient of
//!   `x^j` or of `n^i` in `f` fixed to 1 instead of `f(0, 0)`, inside the
//!   shape of `f`.
//!
//! The left-hand factors are those of the historical passes (the
//! `try_denominator` shapes, or the escalation shapes from six rows on).
//! Because the stage runs only after the historical passes, every result the
//! historical search finds is returned unchanged, with the same candidate
//! count; a budget-limited search that ran out before the third stage is
//! unchanged too.  Row generation divides exactly by the whole left-hand
//! factor (a non-divisible row is
//! [`RecurrenceEvaluationError::NonPolynomialQuotient`]), and the rule is an
//! ordinary [`Recurrence`] with that factor as its denominator, so its JSON
//! form is unchanged.  The modular prefilter and the degree prechecks model
//! the historical normalization and are skipped in the third stage.

use super::*;
use num_traits::Zero;
use serde::{Deserialize, Serialize};
#[cfg(not(target_arch = "wasm32"))]
use std::time::Instant;

/// Options controlling the adaptive search bounds.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AdaptiveSearchOptions {
    /// Number of initial polynomials to ignore before searching.
    pub skip_prefix: usize,
    /// Minimum recurrence length to try.
    pub min_rec_len: usize,
    /// Maximum recurrence length to try.
    pub max_rec_len: usize,
    /// Minimum degree in t for coefficients.
    pub min_var_deg: usize,
    /// Maximum degree in t for coefficients.
    pub max_var_deg: usize,
    /// Minimum degree in n for coefficients.
    pub min_idx_deg: usize,
    /// Maximum degree in n for coefficients.
    pub max_idx_deg: usize,
    /// Minimum derivative order.
    pub min_diff_deg: usize,
    /// Maximum derivative order.
    pub max_diff_deg: usize,
    /// Also search inhomogeneous recurrences.
    pub try_inhomogeneous: bool,
    /// Minimum degree in t for the inhomogeneous term when enabled.
    pub min_inhomo_var_deg: usize,
    /// Maximum degree in t for the inhomogeneous term when enabled.
    pub max_inhomo_var_deg: usize,
    /// Minimum degree in n for the inhomogeneous term when enabled.
    pub min_inhomo_idx_deg: usize,
    /// Maximum degree in n for the inhomogeneous term when enabled.
    pub max_inhomo_idx_deg: usize,
    /// Also search with LHS denominators.
    pub try_denominator: bool,
    /// Also search recurrence terms multiplied by (-1)^n.
    pub try_alternating_sign: bool,
    /// Maximum denom_var_deg when try_denominator is true.
    pub max_denom_var_deg: usize,
    /// Maximum denom_idx_deg when try_denominator is true.
    pub max_denom_idx_deg: usize,
    /// Minimum surplus: equations - unknowns must be >= this.
    pub min_margin: usize,
    /// Use all input rows for fitting and skip held-out verification.
    pub no_verify: bool,
    /// Extra rows to add after the first prefix that clears `min_margin`.
    pub fit_extra_rows: usize,
    /// Reject inconsistent candidates modulo large primes only when full
    /// modular column rank certifies the rejection over the rationals.
    pub modular_prefilter: bool,
    /// Print each candidate tried to stderr.
    pub verbose: bool,
    /// Vary the left-hand side after the historical passes, in a third
    /// complexity-ordered stage (module documentation, "Left-hand side
    /// variation"): a left-hand monomial `x^e` with
    /// `e <= max_laurent_x_power` ([`RecurrenceOptions::leader_x_power`]) and
    /// normalization charts of the left-hand factor other than its constant
    /// term (factors `f(n, t)` from `try_denominator`, or from the
    /// escalation shapes from six rows on).  Results of the historical
    /// passes are unchanged; off skips the third stage.  Default false in
    /// the library, so that the callers inside polytool (the lift searches,
    /// for example) keep their outcomes; the `polytool recurrence` command
    /// and the MCP recurrence tools turn it on by default.
    #[serde(default)]
    pub try_laurent: bool,
    /// Largest power `e` of the left-hand monomial `x^e` when `try_laurent`
    /// is on.  Default 2.
    #[serde(default = "default_max_laurent_x_power")]
    pub max_laurent_x_power: usize,
}

fn default_max_laurent_x_power() -> usize {
    2
}

impl Default for AdaptiveSearchOptions {
    fn default() -> Self {
        Self {
            skip_prefix: 0,
            min_rec_len: 1,
            max_rec_len: 5,
            min_var_deg: 0,
            max_var_deg: 3,
            min_idx_deg: 0,
            max_idx_deg: 3,
            min_diff_deg: 0,
            max_diff_deg: 2,
            try_inhomogeneous: false,
            min_inhomo_var_deg: 0,
            max_inhomo_var_deg: 3,
            min_inhomo_idx_deg: 0,
            max_inhomo_idx_deg: 3,
            try_denominator: false,
            try_alternating_sign: false,
            max_denom_var_deg: 2,
            max_denom_idx_deg: 2,
            min_margin: 1,
            no_verify: false,
            fit_extra_rows: 1,
            modular_prefilter: true,
            verbose: false,
            try_laurent: false,
            max_laurent_x_power: 2,
        }
    }
}

impl AdaptiveSearchOptions {
    /// Return these search options with alternating `(-1)^n` terms enabled or
    /// disabled.
    pub fn with_alternating_sign(mut self, try_alternating_sign: bool) -> Self {
        self.try_alternating_sign = try_alternating_sign;
        self
    }

    /// Return these search options with the left-hand side variation
    /// ([`Self::try_laurent`]) enabled or disabled.
    pub fn with_laurent(mut self, try_laurent: bool) -> Self {
        self.try_laurent = try_laurent;
        self
    }
}

/// Result of an adaptive recurrence search, including metadata.
#[derive(Debug, Clone)]
pub struct AdaptiveSearchResult {
    /// The recurrence found.
    pub recurrence: Recurrence,
    /// The options that produced it.
    pub opts: RecurrenceOptions,
    /// Number of unknowns in the winning system.
    pub num_unknowns: usize,
    /// Weighted unknown count used for adaptive ordering.
    pub weighted_unknowns: usize,
    /// Number of equations in the winning system.
    pub num_equations: usize,
    /// Number of polynomials used to fit the winning system after skip_prefix.
    pub fit_polynomials: usize,
    /// Number of held-out polynomials verified after fitting.
    pub verification_polynomials: usize,
    /// Dense coefficient equations checked on held-out rows, before cancellation.
    /// This counts slots through the largest degree of an evaluated summand,
    /// including zero slots; it is not a rank or an independence estimate.
    pub verification_equations: usize,
    /// Number of candidates actually solved (not just counted).
    pub candidates_tried: usize,
    /// Counters for the adaptive search stages up to the returned result.
    pub diagnostics: AdaptiveSearchDiagnostics,
}

impl AdaptiveSearchResult {
    /// Number of leading rows (after `skip_prefix`, out of
    /// `searched_rows`) that an exported rule must keep as initial rows.
    ///
    /// The fit checks the rule only from row `opts.rec_len` on (the search
    /// window), which can be later than the rule's own lag
    /// (`Recurrence::max_offset`) when the top window coefficients vanish;
    /// rows before the window were never checked and must stay seeds.
    /// Rows where a left-hand factor vanishes are kept as well
    /// ([`Recurrence::generation_initial_count`]).
    pub fn export_initial_count(&self, searched_rows: usize) -> usize {
        self.recurrence
            .generation_initial_count(1, searched_rows)
            .max(self.opts.rec_len)
            .min(searched_rows)
    }
}

/// Deterministic limits for an adaptive recurrence search.
///
/// A candidate is counted when its parameter configuration is taken from the
/// adaptive iterator, before fit-row, structural, modular, or exact-solve
/// filtering.  This bounds the same outer search loop on every machine.
///
/// The work limit bounds what happens inside the candidates.  Work is counted
/// in abstract units that depend only on the dimensions of the systems that
/// are solved, never on the clock, so a search with the same inputs and the
/// same limit stops at the same step on every machine.  A stage is charged
/// before it runs and runs only if its whole cost fits in the remaining limit;
/// the first stage that does not fit ends the search as
/// [`AdaptiveSearchOutcome::BudgetExhausted`], and the units used never exceed
/// the limit.  The unit model:
///
/// * a linear solve (each modular prime, and the exact solve) of `r`
///   equations in `c` unknowns costs `r * c * (min(r, c) + 1)`, the
///   operations of dense elimination and the assembly of its matrix;
/// * checking a candidate against `r` equations costs `r` times the number
///   of unknowns (modular held-out check) or the number of equations itself
///   (exact lift and held-out checks);
/// * the derivative rows and modular images cost one unit per coefficient,
///   derivative order and prime;
/// * each candidate costs one unit, and its structural prechecks one unit
///   per equation.
///
/// The units used are reported in [`AdaptiveSearchDiagnostics::work_units`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct AdaptiveSearchBudget {
    /// Maximum number of candidate configurations to inspect. `None` preserves
    /// the historical unbounded behavior; `Some(0)` inspects none.
    pub max_candidates: Option<usize>,
    /// Maximum number of deterministic work units (see the type
    /// documentation).  `None` is unlimited; `Some(0)` allows no work at all.
    pub max_work: Option<u64>,
}

impl AdaptiveSearchBudget {
    /// An unbounded search, equivalent to the historical adaptive API.
    pub const fn unbounded() -> Self {
        Self {
            max_candidates: None,
            max_work: None,
        }
    }

    /// Inspect at most `max_candidates` candidate configurations.
    pub const fn limited(max_candidates: usize) -> Self {
        Self {
            max_candidates: Some(max_candidates),
            max_work: None,
        }
    }

    /// Return this budget with a limit of `max_work` deterministic work units.
    pub const fn with_max_work(mut self, max_work: u64) -> Self {
        self.max_work = Some(max_work);
        self
    }

    pub(super) fn allows_candidate(self, candidates_considered: usize) -> bool {
        self.max_candidates
            .is_none_or(|limit| candidates_considered < limit)
    }
}

/// Search counters returned when no recurrence result is available.
#[derive(Debug, Clone)]
pub struct AdaptiveSearchSummary {
    /// Number of candidates that reached modular or exact solving.
    pub candidates_tried: usize,
    /// Counters for every adaptive-search stage reached before termination.
    pub diagnostics: AdaptiveSearchDiagnostics,
}

/// Exact termination outcome for a budget-aware adaptive search.
#[derive(Debug, Clone)]
pub enum AdaptiveSearchOutcome {
    /// A recurrence was found within the configured search and budget.
    Found(AdaptiveSearchResult),
    /// The complete configured search space was checked without a recurrence.
    NoRecurrence(AdaptiveSearchSummary),
    /// The candidate budget was consumed and at least one candidate remained.
    BudgetExhausted(AdaptiveSearchSummary),
}

/// Counters describing how adaptive recurrence search reached a result.
#[derive(Debug, Clone, Default)]
pub struct AdaptiveSearchDiagnostics {
    /// Number of parameter candidates represented by the searched lazy iterators.
    pub generated_candidates: usize,
    /// Number of candidates inspected before the search returned.
    pub considered_candidates: usize,
    /// Deterministic work units charged before the search returned (see
    /// [`AdaptiveSearchBudget`] for the unit model).  Identical on every
    /// machine for the same input, options and budget.
    pub work_units: u64,
    /// Candidates skipped because no prefix had enough rows/equations.
    pub insufficient_fit_rows: usize,
    /// Candidates rejected by the equation-count margin after selecting a fit prefix.
    pub equation_bound_rejections: usize,
    /// Homogeneous candidates rejected by the degree-bound precheck.
    pub degree_bound_rejections: usize,
    /// Candidates rejected by a structural zero-row precheck.
    pub structural_zero_rejections: usize,
    /// Candidates rejected by the modular prefilter before exact rational solving.
    pub modular_prefilter_rejections: usize,
    /// Candidates sent to exact rational recurrence solving.
    pub exact_solve_attempts: usize,
    /// Modular solutions lifted as centered integer candidates.
    pub modular_lift_attempts: usize,
    /// Lifted modular candidates accepted after exact verification.
    pub modular_lift_successes: usize,
    /// Lifted modular candidates rejected by exact verification.
    pub modular_lift_failures: usize,
    /// Exact rational solves that did not produce a recurrence.
    pub failed_exact_solves: usize,
    /// Prefix fits that failed held-out verification.
    pub heldout_verification_failures: usize,
    /// Whether the automatic denominator escalation pass was entered.
    pub denominator_escalation_entered: bool,
    /// Time spent precomputing derivatives.
    pub derivative_precompute_ms: f64,
    /// Time spent building modular-prefilter cache data.
    pub modular_cache_build_ms: f64,
    /// Time spent constructing score-ordered candidate iterators.
    pub candidate_generation_ms: f64,
    /// Time spent selecting fit prefixes.
    pub fit_selection_ms: f64,
    /// Time spent counting equations.
    pub equation_count_ms: f64,
    /// Time spent in degree-bound prechecks.
    pub degree_bound_ms: f64,
    /// Time spent in structural zero-row prechecks.
    pub structural_zero_ms: f64,
    /// Time spent in modular-prefilter checks.
    pub modular_prefilter_ms: f64,
    /// Time spent lifting and exactly checking modular candidate solutions.
    pub modular_lift_ms: f64,
    /// Time spent in exact rational recurrence solving.
    pub exact_solve_ms: f64,
    /// Time spent verifying held-out rows.
    pub heldout_verify_ms: f64,
}

#[cfg(not(target_arch = "wasm32"))]
type RecurrenceTimer = Instant;

#[cfg(target_arch = "wasm32")]
#[derive(Clone, Copy)]
struct RecurrenceTimer;

#[cfg(not(target_arch = "wasm32"))]
fn timer_now() -> RecurrenceTimer {
    Instant::now()
}

#[cfg(target_arch = "wasm32")]
fn timer_now() -> RecurrenceTimer {
    RecurrenceTimer
}

#[cfg(not(target_arch = "wasm32"))]
fn elapsed_ms(start: RecurrenceTimer) -> f64 {
    start.elapsed().as_secs_f64() * 1000.0
}

#[cfg(target_arch = "wasm32")]
fn elapsed_ms(_: RecurrenceTimer) -> f64 {
    0.0
}

pub(super) fn fitting_polynomial_count_with_prefix_degrees(
    polys: &[Vec<BigRational>],
    opts: &RecurrenceOptions,
    search: &AdaptiveSearchOptions,
    prefix_max_degrees: &[usize],
) -> Option<usize> {
    let m = polys.len();
    if m <= opts.rec_len {
        return None;
    }

    if search.no_verify {
        return Some(m);
    }

    let max_fit = m.checked_sub(1)?;
    if max_fit <= opts.rec_len {
        return None;
    }

    let unknowns = count_unknowns(opts);
    for fit_len in opts.rec_len + 1..=max_fit {
        let equations = count_equations_from_prefix_degrees(fit_len, opts, prefix_max_degrees);
        if equations >= unknowns + search.min_margin {
            return Some((fit_len + search.fit_extra_rows).min(max_fit));
        }
    }

    None
}

fn verify_heldout_tail(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    rec: &Recurrence,
    opts: &RecurrenceOptions,
    fit_len: usize,
    search: &AdaptiveSearchOptions,
) -> bool {
    search.no_verify
        || recurrence_holds_from_rational_with_derivs(polys, derivs, rec, opts.rec_len, fit_len + 1)
}

pub(super) fn heldout_equation_count(
    polys: &[Vec<BigRational>],
    rec: &Recurrence,
    fit_len: usize,
) -> usize {
    // Count coefficient slots of the actual identity, not the larger search
    // ansatz. Keep the degree bound before summands cancel to zero.
    let length = |p: &[BigRational]| p.iter().rposition(|c| !c.is_zero()).map_or(0, |i| i + 1);
    let product_length = |a: usize, b: usize| {
        if a == 0 || b == 0 {
            0
        } else {
            a + b - 1
        }
    };
    (fit_len..polys.len())
        .map(|row| {
            let n = row + 1;
            let mut slots = rec.denominator.as_ref().map_or(length(&polys[row]), |p| {
                product_length(length(&bivar_eval_n(p, n)), length(&polys[row]))
            });
            for term in &rec.terms {
                let derivative_length =
                    length(&polys[row - term.offset]).saturating_sub(term.deriv_order);
                slots = slots.max(product_length(
                    length(&bivar_eval_n(&term.coeff, n)),
                    derivative_length,
                ));
            }
            if let Some(p) = &rec.inhomogeneous {
                slots = slots.max(length(&bivar_eval_n(p, n)));
            }
            slots.max(1)
        })
        .sum()
}

fn homogeneous_degree_bound_rejects(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    opts: &RecurrenceOptions,
    fit_len: usize,
) -> bool {
    if !opts.homogeneous || opts.denom_var_deg > 0 || opts.denom_idx_deg > 0 {
        return false;
    }

    for nn in opts.rec_len + 1..=fit_len {
        let current = &polys[nn - 1];
        if poly_is_zero_rational(current) {
            continue;
        }
        let lhs_degree = poly_degree_rational(current);
        let mut rhs_degree_bound: Option<usize> = None;
        for r in 1..=opts.rec_len {
            for ref_poly in derivs[nn - 1 - r].iter().take(opts.diff_deg + 1) {
                if !poly_is_zero_rational(ref_poly) {
                    let degree = poly_degree_rational(ref_poly) + opts.var_deg;
                    rhs_degree_bound =
                        Some(rhs_degree_bound.map_or(degree, |bound| bound.max(degree)));
                }
            }
        }
        if rhs_degree_bound.is_none_or(|degree| degree < lhs_degree) {
            return true;
        }
    }

    false
}

fn recurrence_row_has_possible_source(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    opts: &RecurrenceOptions,
    nn: usize,
    l: usize,
) -> bool {
    let current = &polys[nn - 1];

    if opts.denom_var_deg > 0 || opts.denom_idx_deg > 0 {
        for j in 0..=opts.denom_var_deg {
            if j == 0 && opts.denom_idx_deg == 0 {
                continue;
            }
            if l >= j && poly_coeff_nonzero_rational(current, l - j) {
                return true;
            }
        }
    }

    for r in 1..=opts.rec_len {
        for ref_poly in derivs[nn - 1 - r].iter().take(opts.diff_deg + 1) {
            for j in 0..=opts.var_deg {
                if l >= j && poly_coeff_nonzero_rational(ref_poly, l - j) {
                    return true;
                }
            }
        }
    }

    !opts.homogeneous && l <= opts.inhomo_var_deg
}

fn structural_zero_row_rejects(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    opts: &RecurrenceOptions,
    fit_len: usize,
    prefix_max_degrees: &[usize],
) -> bool {
    if fit_len <= opts.rec_len {
        return false;
    }

    let max_poly_deg = prefix_max_degrees.get(fit_len).copied().unwrap_or(0);
    let max_j = opts
        .var_deg
        .max(opts.denom_var_deg)
        .max(if opts.homogeneous {
            0
        } else {
            opts.inhomo_var_deg
        });
    let max_t_deg = max_j + max_poly_deg;

    for nn in opts.rec_len + 1..=fit_len {
        let current = &polys[nn - 1];
        for l in 0..=max_t_deg {
            if poly_coeff_nonzero_rational(current, l)
                && !recurrence_row_has_possible_source(polys, derivs, opts, nn, l)
            {
                return true;
            }
        }
    }

    false
}

pub(super) struct FoundAdaptiveCandidate {
    pub(super) recurrence: Recurrence,
    pub(super) opts: RecurrenceOptions,
    pub(super) num_unknowns: usize,
    pub(super) weighted_unknowns: usize,
    pub(super) num_equations: usize,
    pub(super) fit_polynomials: usize,
    pub(super) verification_polynomials: usize,
    pub(super) verification_equations: usize,
    pub(super) candidates_tried: usize,
}

impl FoundAdaptiveCandidate {
    fn into_result(self, diagnostics: AdaptiveSearchDiagnostics) -> AdaptiveSearchResult {
        AdaptiveSearchResult {
            recurrence: self.recurrence,
            opts: self.opts,
            num_unknowns: self.num_unknowns,
            weighted_unknowns: self.weighted_unknowns,
            num_equations: self.num_equations,
            fit_polynomials: self.fit_polynomials,
            verification_polynomials: self.verification_polynomials,
            verification_equations: self.verification_equations,
            candidates_tried: self.candidates_tried,
            diagnostics,
        }
    }
}

pub(super) struct AdaptiveCandidateContext<'a> {
    pub(super) polys: &'a [Vec<BigRational>],
    pub(super) derivs: &'a [Vec<Vec<BigRational>>],
    pub(super) prefix_max_degrees: &'a [usize],
    pub(super) modular_cache: Option<&'a ModularPrefilterCache>,
    pub(super) search: &'a AdaptiveSearchOptions,
    pub(super) fit_search: &'a AdaptiveSearchOptions,
    pub(super) pass_label: Option<&'static str>,
    pub(super) work: &'a WorkMeter,
}

fn adaptive_pass_suffix(pass_label: Option<&str>) -> &'static str {
    match pass_label {
        Some("rational") => " (rational)",
        Some("left-hand") => " (left-hand variant)",
        _ => "",
    }
}

pub(super) fn evaluate_adaptive_candidate(
    ctx: &AdaptiveCandidateContext<'_>,
    opts: &RecurrenceOptions,
    diagnostics: &mut AdaptiveSearchDiagnostics,
    tried: &mut usize,
) -> Option<FoundAdaptiveCandidate> {
    evaluate_adaptive_candidate_in_chart(ctx, opts, LeaderChart::CONSTANT, diagnostics, tried)
}

/// [`evaluate_adaptive_candidate`] in a left-hand chart.  The prechecks and
/// the modular prefilter model the historical left-hand factor (constant
/// chart, no monomial `x^e`), so they run only for it; the other candidates
/// go straight to the exact solve and the held-out verification.
pub(super) fn evaluate_adaptive_candidate_in_chart(
    ctx: &AdaptiveCandidateContext<'_>,
    opts: &RecurrenceOptions,
    chart: LeaderChart,
    diagnostics: &mut AdaptiveSearchDiagnostics,
    tried: &mut usize,
) -> Option<FoundAdaptiveCandidate> {
    let found = evaluate_candidate_stages(ctx, opts, chart, diagnostics, tried);
    diagnostics.work_units = ctx.work.used();
    found
}

/// The stages of one candidate.  Each stage that does work is charged to
/// `ctx.work` first; a stage that does not fit ends the candidate with
/// `None`, and the caller sees the exhausted meter.
fn evaluate_candidate_stages(
    ctx: &AdaptiveCandidateContext<'_>,
    opts: &RecurrenceOptions,
    chart: LeaderChart,
    diagnostics: &mut AdaptiveSearchDiagnostics,
    tried: &mut usize,
) -> Option<FoundAdaptiveCandidate> {
    let work = ctx.work;
    if !work.try_charge(1) {
        return None;
    }
    let historical = chart.is_constant() && opts.leader_x_power == 0;
    diagnostics.considered_candidates += 1;
    let m = ctx.polys.len();
    let complexity = candidate_complexity(opts);
    let unknowns = complexity.raw_unknowns;
    let weighted_unknowns = complexity.weighted_unknowns;

    let timer = timer_now();
    let fit_len = fitting_polynomial_count_with_prefix_degrees(
        ctx.polys,
        opts,
        ctx.fit_search,
        ctx.prefix_max_degrees,
    );
    diagnostics.fit_selection_ms += elapsed_ms(timer);
    let Some(fit_len) = fit_len else {
        diagnostics.insufficient_fit_rows += 1;
        return None;
    };

    let fit_polys = &ctx.polys[..fit_len];
    let timer = timer_now();
    let equations = count_equations_from_prefix_degrees(fit_len, opts, ctx.prefix_max_degrees);
    diagnostics.equation_count_ms += elapsed_ms(timer);
    if equations < unknowns + ctx.search.min_margin {
        diagnostics.equation_bound_rejections += 1;
        return None;
    }

    if historical {
        if !work.try_charge(equations as u64) {
            return None;
        }
        let timer = timer_now();
        let degree_bound_rejects =
            homogeneous_degree_bound_rejects(ctx.polys, ctx.derivs, opts, fit_len);
        diagnostics.degree_bound_ms += elapsed_ms(timer);
        if degree_bound_rejects {
            diagnostics.degree_bound_rejections += 1;
            return None;
        }

        let timer = timer_now();
        let structural_zero_rejects = structural_zero_row_rejects(
            ctx.polys,
            ctx.derivs,
            opts,
            fit_len,
            ctx.prefix_max_degrees,
        );
        diagnostics.structural_zero_ms += elapsed_ms(timer);
        if structural_zero_rejects {
            diagnostics.structural_zero_rejections += 1;
            return None;
        }
    }

    // Held-out rows are those after the fit prefix: one unit per equation.
    let heldout_work = if ctx.search.no_verify {
        0
    } else {
        count_equations_from_prefix_degrees(m, opts, ctx.prefix_max_degrees)
            .saturating_sub(equations) as u64
    };

    *tried += 1;

    if ctx.search.verbose {
        let label = adaptive_pass_suffix(ctx.pass_label);
        eprintln!(
            "  try #{tried}{label}: rec_len={} var_deg={} idx_deg={} diff_deg={} \
             alternating={} denom=({},{}) homog={} inhomo=({},{}) \
             fit_rows={} verify_rows={} \
             (unknowns={unknowns}, weighted={weighted_unknowns}, \
             equations={equations}, margin={})",
            opts.rec_len,
            opts.var_deg,
            opts.idx_deg,
            opts.diff_deg,
            opts.alternating_sign,
            opts.denom_var_deg,
            opts.denom_idx_deg,
            opts.homogeneous,
            opts.inhomo_var_deg,
            opts.inhomo_idx_deg,
            fit_len,
            m - fit_len,
            equations - unknowns,
        );
    }

    let fit_derivs = &ctx.derivs[..fit_len];
    let mut exact_row_indices = None;
    let mut modular_candidate = None;
    let cached_solve_opts = if let Some(cache) = ctx.modular_cache.filter(|_| historical) {
        let timer = timer_now();
        let prefilter = modular_prefilter_with_cache(ctx.polys, fit_len, opts, cache, Some(work));
        diagnostics.modular_prefilter_ms += elapsed_ms(timer);
        if work.is_exhausted() {
            return None;
        }
        if prefilter.rejected {
            diagnostics.modular_prefilter_rejections += 1;
            return None;
        }
        exact_row_indices = prefilter.full_rank_pivot_rows;
        modular_candidate = prefilter.modular_solution;
        let mut solve_opts = opts.clone();
        solve_opts.modular_prefilter = false;
        Some(solve_opts)
    } else {
        None
    };
    let solve_opts = cached_solve_opts.as_ref().unwrap_or(opts);

    if let Some(candidate) = modular_candidate.as_ref() {
        if !work.try_charge((equations + unknowns) as u64) {
            return None;
        }
        diagnostics.modular_lift_attempts += 1;
        let timer = timer_now();
        let lifted = verified_recurrence_from_modular_lift(fit_polys, fit_derivs, opts, candidate);
        diagnostics.modular_lift_ms += elapsed_ms(timer);
        if let Some(rec) = lifted {
            if !work.try_charge(heldout_work) {
                return None;
            }
            let timer = timer_now();
            let verified =
                verify_heldout_tail(ctx.polys, ctx.derivs, &rec, opts, fit_len, ctx.search);
            diagnostics.heldout_verify_ms += elapsed_ms(timer);
            if verified {
                diagnostics.modular_lift_successes += 1;
                if ctx.search.verbose {
                    let suffix = adaptive_pass_suffix(ctx.pass_label);
                    eprintln!("  -> found by modular lift{suffix}!");
                }
                return Some(FoundAdaptiveCandidate {
                    verification_equations: heldout_equation_count(ctx.polys, &rec, fit_len),
                    recurrence: rec,
                    opts: opts.clone(),
                    num_unknowns: unknowns,
                    weighted_unknowns,
                    num_equations: equations,
                    fit_polynomials: fit_len,
                    verification_polynomials: m - fit_len,
                    candidates_tried: *tried,
                });
            }
            diagnostics.heldout_verification_failures += 1;
            if ctx.search.verbose {
                eprintln!("  -> modular lift fit prefix but failed held-out verification");
            }
            diagnostics.modular_lift_failures += 1;
            if candidate.full_rank {
                return None;
            }
        } else {
            diagnostics.modular_lift_failures += 1;
        }
    }

    diagnostics.exact_solve_attempts += 1;
    let timer = timer_now();
    let rec = find_polynomial_recurrence_in_leader_chart(
        fit_polys,
        fit_derivs,
        solve_opts,
        exact_row_indices.as_deref(),
        chart,
        Some(work),
    );
    diagnostics.exact_solve_ms += elapsed_ms(timer);
    if work.is_exhausted() {
        return None;
    }
    if let Some(rec) = rec {
        if !work.try_charge(heldout_work) {
            return None;
        }
        let timer = timer_now();
        let verified = verify_heldout_tail(ctx.polys, ctx.derivs, &rec, opts, fit_len, ctx.search);
        diagnostics.heldout_verify_ms += elapsed_ms(timer);
        if !verified {
            diagnostics.heldout_verification_failures += 1;
            if ctx.search.verbose {
                eprintln!("  -> fitted prefix but failed held-out verification");
            }
            return None;
        }
        if ctx.search.verbose {
            let suffix = adaptive_pass_suffix(ctx.pass_label);
            eprintln!("  -> found{suffix}!");
        }
        Some(FoundAdaptiveCandidate {
            verification_equations: heldout_equation_count(ctx.polys, &rec, fit_len),
            recurrence: rec,
            opts: opts.clone(),
            num_unknowns: unknowns,
            weighted_unknowns,
            num_equations: equations,
            fit_polynomials: fit_len,
            verification_polynomials: m - fit_len,
            candidates_tried: *tried,
        })
    } else {
        diagnostics.failed_exact_solves += 1;
        None
    }
}

/// Search for the simplest polynomial recurrence by trying parameter
/// combinations in order of ascending complexity.
pub fn find_recurrence_adaptive_rational_with_budget(
    polys: &[Vec<BigRational>],
    search: &AdaptiveSearchOptions,
    budget: AdaptiveSearchBudget,
) -> AdaptiveSearchOutcome {
    let polys = polys.get(search.skip_prefix..).unwrap_or(&[]);
    let m = polys.len();
    if m < 2 {
        return AdaptiveSearchOutcome::NoRecurrence(AdaptiveSearchSummary {
            candidates_tried: 0,
            diagnostics: AdaptiveSearchDiagnostics::default(),
        });
    }
    let mut diagnostics = AdaptiveSearchDiagnostics::default();
    let timer = timer_now();
    let candidates = CandidateIterator::new(m, search);
    diagnostics.candidate_generation_ms += elapsed_ms(timer);
    diagnostics.generated_candidates = candidates.total_candidates();
    if (budget.max_candidates == Some(0) || budget.max_work == Some(0))
        && diagnostics.generated_candidates > 0
    {
        return AdaptiveSearchOutcome::BudgetExhausted(AdaptiveSearchSummary {
            candidates_tried: 0,
            diagnostics,
        });
    }

    // The precomputation is charged first and skipped when it does not fit.
    let work = WorkMeter::new(budget.max_work);
    let coefficient_count: u64 = polys.iter().map(|poly| poly.len().max(1) as u64).sum();
    let derivative_work = coefficient_count.saturating_mul(search.max_diff_deg as u64 + 1);
    let cache_work = if search.modular_prefilter {
        ModularPrefilterCache::construction_work(polys, search.max_diff_deg)
    } else {
        0
    };
    if !work.try_charge(derivative_work) || !work.try_charge(cache_work) {
        diagnostics.work_units = work.used();
        return AdaptiveSearchOutcome::BudgetExhausted(AdaptiveSearchSummary {
            candidates_tried: 0,
            diagnostics,
        });
    }
    diagnostics.work_units = work.used();

    let prefix_max_degrees = prefix_max_degrees_rational(polys);
    let timer = timer_now();
    let derivs = rational_derivatives_up_to(polys, search.max_diff_deg);
    diagnostics.derivative_precompute_ms += elapsed_ms(timer);
    let timer = timer_now();
    let modular_cache = if search.modular_prefilter {
        Some(ModularPrefilterCache::new(polys, search.max_diff_deg))
    } else {
        None
    };
    diagnostics.modular_cache_build_ms += elapsed_ms(timer);

    // First pass: use the options as given.
    let mut tried = 0;

    let primary_context = AdaptiveCandidateContext {
        polys,
        derivs: &derivs,
        prefix_max_degrees: &prefix_max_degrees,
        modular_cache: modular_cache.as_ref(),
        search,
        fit_search: search,
        pass_label: None,
        work: &work,
    };
    for opts in candidates {
        if !budget.allows_candidate(diagnostics.considered_candidates) {
            return AdaptiveSearchOutcome::BudgetExhausted(AdaptiveSearchSummary {
                candidates_tried: tried,
                diagnostics,
            });
        }
        if let Some(found) =
            evaluate_adaptive_candidate(&primary_context, &opts, &mut diagnostics, &mut tried)
        {
            return AdaptiveSearchOutcome::Found(found.into_result(diagnostics));
        }
        if work.is_exhausted() {
            return AdaptiveSearchOutcome::BudgetExhausted(AdaptiveSearchSummary {
                candidates_tried: tried,
                diagnostics,
            });
        }
    }

    // Second pass: if denominator wasn't already tried, automatically
    // escalate to rational coefficients (LHS denominator f(n,t)).
    // Only try if we have enough data to avoid spurious fits.
    if !search.try_denominator && m >= 6 {
        let mut rational_search = search.clone();
        rational_search.try_denominator = true;
        // Use moderate denominator degrees if not already set.
        if rational_search.max_denom_idx_deg == 0 {
            rational_search.max_denom_idx_deg = 1;
        }

        let timer = timer_now();
        let rational_candidates = CandidateIterator::new(m, &rational_search);
        diagnostics.candidate_generation_ms += elapsed_ms(timer);
        diagnostics.denominator_escalation_entered = true;
        diagnostics.generated_candidates = diagnostics
            .generated_candidates
            .saturating_add(rational_candidates.denominator_candidates());
        let rational_context = AdaptiveCandidateContext {
            polys,
            derivs: &derivs,
            prefix_max_degrees: &prefix_max_degrees,
            modular_cache: modular_cache.as_ref(),
            search,
            fit_search: &rational_search,
            pass_label: Some("rational"),
            work: &work,
        };
        for opts in rational_candidates {
            // Skip candidates without a denominator (already tried above).
            if opts.denom_var_deg == 0 && opts.denom_idx_deg == 0 {
                continue;
            }

            if !budget.allows_candidate(diagnostics.considered_candidates) {
                return AdaptiveSearchOutcome::BudgetExhausted(AdaptiveSearchSummary {
                    candidates_tried: tried,
                    diagnostics,
                });
            }
            if let Some(found) =
                evaluate_adaptive_candidate(&rational_context, &opts, &mut diagnostics, &mut tried)
            {
                return AdaptiveSearchOutcome::Found(found.into_result(diagnostics));
            }
            if work.is_exhausted() {
                return AdaptiveSearchOutcome::BudgetExhausted(AdaptiveSearchSummary {
                    candidates_tried: tried,
                    diagnostics,
                });
            }
        }
    }

    // Third stage: the left-hand variants, in their own complexity order.
    if search.try_laurent {
        return left_hand_variation_search(
            polys,
            search,
            budget,
            diagnostics,
            tried,
            &derivs,
            &prefix_max_degrees,
            &work,
        );
    }

    AdaptiveSearchOutcome::NoRecurrence(AdaptiveSearchSummary {
        candidates_tried: tried,
        diagnostics,
    })
}

/// The complexity-ordered enumeration of the left-hand variants
/// (module documentation, "Left-hand side variation").
#[allow(clippy::too_many_arguments)]
fn left_hand_variation_search(
    polys: &[Vec<BigRational>],
    search: &AdaptiveSearchOptions,
    budget: AdaptiveSearchBudget,
    mut diagnostics: AdaptiveSearchDiagnostics,
    mut tried: usize,
    derivs: &[Vec<Vec<BigRational>>],
    prefix_max_degrees: &[usize],
    work: &WorkMeter,
) -> AdaptiveSearchOutcome {
    let m = polys.len();
    // The left-hand factors of the variants are those of the historical
    // passes: `try_denominator`'s shapes, or the escalation shapes from six
    // rows on (the historical guard against spurious fits).
    let mut effective = search.clone();
    if !search.try_denominator && m >= 6 {
        effective.try_denominator = true;
        if effective.max_denom_idx_deg == 0 {
            effective.max_denom_idx_deg = 1;
        }
    }
    let timer = timer_now();
    let mut candidates = CandidateIterator::with_charts(m, &effective);
    diagnostics.candidate_generation_ms += elapsed_ms(timer);
    diagnostics.generated_candidates = diagnostics
        .generated_candidates
        .saturating_add(candidates.variant_candidates());
    let context = AdaptiveCandidateContext {
        polys,
        derivs,
        prefix_max_degrees,
        modular_cache: None,
        search,
        fit_search: &effective,
        pass_label: Some("left-hand"),
        work,
    };
    while let Some((opts, chart)) = candidates.next_candidate() {
        // The historical candidates (constant chart, no monomial x^e) were
        // tried in the first two passes.
        if chart.is_constant() && opts.leader_x_power == 0 {
            continue;
        }
        if !budget.allows_candidate(diagnostics.considered_candidates) {
            return AdaptiveSearchOutcome::BudgetExhausted(AdaptiveSearchSummary {
                candidates_tried: tried,
                diagnostics,
            });
        }
        if let Some(found) = evaluate_adaptive_candidate_in_chart(
            &context,
            &opts,
            chart,
            &mut diagnostics,
            &mut tried,
        ) {
            return AdaptiveSearchOutcome::Found(found.into_result(diagnostics));
        }
        if work.is_exhausted() {
            return AdaptiveSearchOutcome::BudgetExhausted(AdaptiveSearchSummary {
                candidates_tried: tried,
                diagnostics,
            });
        }
    }
    AdaptiveSearchOutcome::NoRecurrence(AdaptiveSearchSummary {
        candidates_tried: tried,
        diagnostics,
    })
}

/// Search without a candidate limit, preserving the historical optional
/// result API.
pub fn find_recurrence_adaptive_rational(
    polys: &[Vec<BigRational>],
    search: &AdaptiveSearchOptions,
) -> Option<AdaptiveSearchResult> {
    match find_recurrence_adaptive_rational_with_budget(
        polys,
        search,
        AdaptiveSearchBudget::unbounded(),
    ) {
        AdaptiveSearchOutcome::Found(result) => Some(result),
        AdaptiveSearchOutcome::NoRecurrence(_) => None,
        AdaptiveSearchOutcome::BudgetExhausted(_) => {
            unreachable!("an unbounded search cannot exhaust its candidate budget")
        }
    }
}

/// Search for the simplest polynomial recurrence for an integer-coefficient
/// polynomial sequence.
///
/// This is a convenience wrapper around [`find_recurrence_adaptive_rational`].
pub fn find_recurrence_adaptive(
    polys: &[Vec<i64>],
    search: &AdaptiveSearchOptions,
) -> Option<AdaptiveSearchResult> {
    let rational_polys = i64_polys_to_rational(polys);
    find_recurrence_adaptive_rational(&rational_polys, search)
}

/// Integer-coefficient convenience wrapper for
/// [`find_recurrence_adaptive_rational_with_budget`].
pub fn find_recurrence_adaptive_with_budget(
    polys: &[Vec<i64>],
    search: &AdaptiveSearchOptions,
    budget: AdaptiveSearchBudget,
) -> AdaptiveSearchOutcome {
    let rational_polys = i64_polys_to_rational(polys);
    find_recurrence_adaptive_rational_with_budget(&rational_polys, search, budget)
}

#[cfg(test)]
mod tests;
