//! Adaptive recurrence search: options, budgets, diagnostics and the search
//! loop over candidate shapes.

use super::*;
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
    /// Number of candidates actually solved (not just counted).
    pub candidates_tried: usize,
    /// Counters for the adaptive search stages up to the returned result.
    pub diagnostics: AdaptiveSearchDiagnostics,
}

/// Deterministic limit for an adaptive recurrence search.
///
/// A candidate is counted when its parameter configuration is taken from the
/// adaptive iterator, before fit-row, structural, modular, or exact-solve
/// filtering.  This bounds the same outer search loop on every machine.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct AdaptiveSearchBudget {
    /// Maximum number of candidate configurations to inspect. `None` preserves
    /// the historical unbounded behavior; `Some(0)` inspects none.
    pub max_candidates: Option<usize>,
}

impl AdaptiveSearchBudget {
    /// An unbounded search, equivalent to the historical adaptive API.
    pub const fn unbounded() -> Self {
        Self {
            max_candidates: None,
        }
    }

    /// Inspect at most `max_candidates` candidate configurations.
    pub const fn limited(max_candidates: usize) -> Self {
        Self {
            max_candidates: Some(max_candidates),
        }
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
}

fn adaptive_pass_suffix(pass_label: Option<&str>) -> &'static str {
    match pass_label {
        Some("rational") => " (rational)",
        _ => "",
    }
}

pub(super) fn evaluate_adaptive_candidate(
    ctx: &AdaptiveCandidateContext<'_>,
    opts: &RecurrenceOptions,
    diagnostics: &mut AdaptiveSearchDiagnostics,
    tried: &mut usize,
) -> Option<FoundAdaptiveCandidate> {
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

    let timer = timer_now();
    let degree_bound_rejects =
        homogeneous_degree_bound_rejects(ctx.polys, ctx.derivs, opts, fit_len);
    diagnostics.degree_bound_ms += elapsed_ms(timer);
    if degree_bound_rejects {
        diagnostics.degree_bound_rejections += 1;
        return None;
    }

    let timer = timer_now();
    let structural_zero_rejects =
        structural_zero_row_rejects(ctx.polys, ctx.derivs, opts, fit_len, ctx.prefix_max_degrees);
    diagnostics.structural_zero_ms += elapsed_ms(timer);
    if structural_zero_rejects {
        diagnostics.structural_zero_rejections += 1;
        return None;
    }

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
    let cached_solve_opts = if let Some(cache) = ctx.modular_cache {
        let timer = timer_now();
        let prefilter = modular_prefilter_with_cache(ctx.polys, fit_len, opts, cache);
        diagnostics.modular_prefilter_ms += elapsed_ms(timer);
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
        diagnostics.modular_lift_attempts += 1;
        let timer = timer_now();
        let lifted = verified_recurrence_from_modular_lift(fit_polys, fit_derivs, opts, candidate);
        diagnostics.modular_lift_ms += elapsed_ms(timer);
        if let Some(rec) = lifted {
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
    let rec = find_polynomial_recurrence_rational_with_derivs_and_rows(
        fit_polys,
        fit_derivs,
        solve_opts,
        exact_row_indices.as_deref(),
    );
    diagnostics.exact_solve_ms += elapsed_ms(timer);
    if let Some(rec) = rec {
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
    if budget.max_candidates == Some(0) && diagnostics.generated_candidates > 0 {
        return AdaptiveSearchOutcome::BudgetExhausted(AdaptiveSearchSummary {
            candidates_tried: 0,
            diagnostics,
        });
    }

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
mod tests {
    use super::*;

    #[test]
    fn adaptive_geometric() {
        let polys: Vec<Vec<i64>> = vec![vec![1], vec![2], vec![4], vec![8], vec![16], vec![32]];
        let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
        assert_eq!(format!("{}", result.recurrence), "P(n) = 2 P(n-1)");
        assert!(result.verification_polynomials >= 1);
        assert!(result.fit_polynomials < polys.len());
    }

    #[test]
    fn homogeneous_degree_bound_rejects_impossible_rhs_degree() {
        let polys = i64_polys_to_rational(&[vec![1], vec![1, 1], vec![1, 2, 1]]);
        let derivs = rational_derivatives_up_to(&polys, 0);
        let too_low = RecurrenceOptions {
            rec_len: 1,
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            homogeneous: true,
            ..Default::default()
        };
        let enough_t_degree = RecurrenceOptions {
            var_deg: 1,
            ..too_low.clone()
        };

        assert!(homogeneous_degree_bound_rejects(
            &polys, &derivs, &too_low, 3
        ));
        assert!(!homogeneous_degree_bound_rejects(
            &polys,
            &derivs,
            &enough_t_degree,
            3
        ));
    }

    #[test]
    fn structural_zero_row_rejects_impossible_support_gap() {
        let polys = i64_polys_to_rational(&[vec![1, 0, 1], vec![0, 1]]);
        let derivs = rational_derivatives_up_to(&polys, 0);
        let prefix_max_degrees = prefix_max_degrees_rational(&polys);
        let no_t_shift = RecurrenceOptions {
            rec_len: 1,
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            homogeneous: true,
            ..Default::default()
        };
        let with_t_shift = RecurrenceOptions {
            var_deg: 1,
            ..no_t_shift.clone()
        };

        assert!(!homogeneous_degree_bound_rejects(
            &polys,
            &derivs,
            &no_t_shift,
            2
        ));
        assert!(structural_zero_row_rejects(
            &polys,
            &derivs,
            &no_t_shift,
            2,
            &prefix_max_degrees
        ));
        assert!(!structural_zero_row_rejects(
            &polys,
            &derivs,
            &with_t_shift,
            2,
            &prefix_max_degrees
        ));
    }

    #[test]
    fn adaptive_no_verify_uses_all_rows() {
        let polys: Vec<Vec<i64>> = vec![vec![1], vec![2], vec![4], vec![8], vec![16], vec![32]];
        let search = AdaptiveSearchOptions {
            no_verify: true,
            ..Default::default()
        };
        let result = find_recurrence_adaptive(&polys, &search).unwrap();
        assert_eq!(format!("{}", result.recurrence), "P(n) = 2 P(n-1)");
        assert_eq!(result.fit_polynomials, polys.len());
        assert_eq!(result.verification_polynomials, 0);
    }

    #[test]
    fn adaptive_rejects_false_prefix_fit_on_heldout_row() {
        let polys: Vec<Vec<i64>> = vec![vec![1], vec![2], vec![4], vec![8], vec![17]];
        let search = AdaptiveSearchOptions {
            max_rec_len: 1,
            max_var_deg: 0,
            max_idx_deg: 0,
            max_diff_deg: 0,
            ..Default::default()
        };
        assert!(find_recurrence_adaptive(&polys, &search).is_none());
    }

    #[test]
    fn weighted_complexity_delays_alternating_terms() {
        let ordinary = RecurrenceOptions {
            rec_len: 2,
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            homogeneous: true,
            ..Default::default()
        };
        let alternating = RecurrenceOptions {
            rec_len: 1,
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            homogeneous: true,
            alternating_sign: true,
            ..Default::default()
        };
        assert_eq!(count_unknowns(&ordinary), count_unknowns(&alternating));
        assert!(count_weighted_unknowns(&ordinary) < count_weighted_unknowns(&alternating));
    }

    #[test]
    fn adaptive_respects_min_bounds() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![2],
            vec![3],
            vec![5],
            vec![8],
            vec![13],
        ];
        let search = AdaptiveSearchOptions {
            min_rec_len: 2,
            max_rec_len: 2,
            min_var_deg: 0,
            max_var_deg: 0,
            min_idx_deg: 0,
            max_idx_deg: 0,
            min_diff_deg: 0,
            max_diff_deg: 0,
            ..Default::default()
        };
        let result = find_recurrence_adaptive(&polys, &search).unwrap();
        assert_eq!(format!("{}", result.recurrence), "P(n) = P(n-1) + P(n-2)");
        assert_eq!(result.opts.rec_len, 2);
    }

    #[test]
    fn adaptive_can_skip_prefix() {
        let polys: Vec<Vec<i64>> = vec![
            vec![9],
            vec![1],
            vec![2],
            vec![4],
            vec![8],
            vec![16],
            vec![32],
        ];
        let search = AdaptiveSearchOptions {
            skip_prefix: 1,
            ..Default::default()
        };
        let result = find_recurrence_adaptive(&polys, &search).unwrap();
        assert_eq!(format!("{}", result.recurrence), "P(n) = 2 P(n-1)");
    }

    // --- Adaptive search tests ---

    #[test]
    fn adaptive_fibonacci() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![2],
            vec![3],
            vec![5],
            vec![8],
            vec![13],
        ];
        let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
        assert_eq!(format!("{}", result.recurrence), "P(n) = P(n-1) + P(n-2)");
        assert_eq!(result.opts.rec_len, 2);
        assert_eq!(result.opts.diff_deg, 0);
    }

    #[test]
    fn adaptive_modular_prefilter_finds_fibonacci() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![2],
            vec![3],
            vec![5],
            vec![8],
            vec![13],
        ];
        let search = AdaptiveSearchOptions {
            modular_prefilter: true,
            ..Default::default()
        };
        let result = find_recurrence_adaptive(&polys, &search).unwrap();
        assert_eq!(format!("{}", result.recurrence), "P(n) = P(n-1) + P(n-2)");
        assert!(result.opts.modular_prefilter);
    }

    #[test]
    fn adaptive_modular_lift_skips_exact_solve_for_integer_solution() {
        let polys: Vec<Vec<i64>> = vec![vec![1], vec![2], vec![4], vec![8], vec![16], vec![32]];
        let search = AdaptiveSearchOptions {
            min_rec_len: 1,
            max_rec_len: 1,
            min_var_deg: 0,
            max_var_deg: 0,
            min_idx_deg: 0,
            max_idx_deg: 0,
            min_diff_deg: 0,
            max_diff_deg: 0,
            modular_prefilter: true,
            ..Default::default()
        };
        let result = find_recurrence_adaptive(&polys, &search).unwrap();
        assert_eq!(format!("{}", result.recurrence), "P(n) = 2 P(n-1)");
        assert_eq!(result.diagnostics.modular_lift_attempts, 1);
        assert_eq!(result.diagnostics.modular_lift_successes, 1);
        assert_eq!(result.diagnostics.exact_solve_attempts, 0);
    }

    #[test]
    fn adaptive_rank_deficient_modular_lift_can_skip_exact_solve() {
        let mut polys = Vec::new();
        let mut value = 1_i64;
        for _ in 0..32 {
            polys.push(vec![value]);
            value *= 2;
        }
        let search = AdaptiveSearchOptions {
            min_rec_len: 5,
            max_rec_len: 5,
            min_var_deg: 3,
            max_var_deg: 3,
            min_idx_deg: 3,
            max_idx_deg: 3,
            min_diff_deg: 0,
            max_diff_deg: 0,
            modular_prefilter: true,
            ..Default::default()
        };
        let result = find_recurrence_adaptive(&polys, &search).unwrap();
        assert!(recurrence_holds_rational(
            &i64_polys_to_rational(&polys),
            &result.recurrence,
            result.opts.rec_len
        ));
        assert_eq!(result.diagnostics.modular_lift_attempts, 1);
        assert_eq!(result.diagnostics.modular_lift_successes, 1);
        assert_eq!(result.diagnostics.exact_solve_attempts, 0);
    }

    #[test]
    fn adaptive_factorial() {
        let polys: Vec<Vec<i64>> = vec![vec![1], vec![2], vec![6], vec![24], vec![120]];
        let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
        assert_eq!(format!("{}", result.recurrence), "P(n) = n P(n-1)");
        assert_eq!(result.opts.rec_len, 1);
        assert_eq!(result.opts.idx_deg, 1);
    }

    #[test]
    fn adaptive_binomial() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1, 1],
            vec![1, 2, 1],
            vec![1, 3, 3, 1],
            vec![1, 4, 6, 4, 1],
        ];
        let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
        assert_eq!(format!("{}", result.recurrence), "P(n) = (1 + t) P(n-1)");
        assert_eq!(result.opts.rec_len, 1);
    }

    #[test]
    fn adaptive_chebyshev() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![0, 1],
            vec![-1, 0, 2],
            vec![0, -3, 0, 4],
            vec![1, 0, -8, 0, 8],
        ];
        let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
        assert_eq!(
            format!("{}", result.recurrence),
            "P(n) = 2t P(n-1) - P(n-2)"
        );
    }

    #[test]
    fn adaptive_eulerian() {
        // With held-out verification, prefix-only non-derivative fits are
        // rejected, and the standard derivative recurrence is found.
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![1, 1],
            vec![1, 4, 1],
            vec![1, 11, 11, 1],
            vec![1, 26, 66, 26, 1],
        ];
        let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
        assert_eq!(
            format!("{}", result.recurrence),
            "P(n) = (1 - 2t + nt) P(n-1) + (t - t^2) P'(n-1)"
        );
        assert_eq!(result.opts.diff_deg, 1);
        assert_eq!(result.opts.rec_len, 1);

        // With max_rec_len=1, the same derivative recurrence is still found.
        let search = AdaptiveSearchOptions {
            max_rec_len: 1,
            max_diff_deg: 2,
            ..Default::default()
        };
        let result = find_recurrence_adaptive(&polys, &search).unwrap();
        assert_eq!(
            format!("{}", result.recurrence),
            "P(n) = (1 - 2t + nt) P(n-1) + (t - t^2) P'(n-1)"
        );
        assert_eq!(result.opts.diff_deg, 1);
    }

    #[test]
    fn adaptive_alternating_sign() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![-1],
            vec![-1],
            vec![1],
            vec![1],
            vec![-1],
            vec![-1],
        ];
        let search = AdaptiveSearchOptions {
            try_alternating_sign: true,
            max_rec_len: 1,
            max_var_deg: 0,
            max_idx_deg: 0,
            max_diff_deg: 0,
            ..Default::default()
        };
        let result = find_recurrence_adaptive(&polys, &search).unwrap();
        assert_eq!(format!("{}", result.recurrence), "P(n) = (-1)^n P(n-1)");
        assert!(result.opts.alternating_sign);
    }

    #[test]
    fn adaptive_short_sequence() {
        // m=3 with constant polys: too few equations for any 2-term recurrence.
        let polys: Vec<Vec<i64>> = vec![vec![1], vec![1], vec![2]];
        let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default());
        assert!(result.is_none());
    }

    fn fixed_adaptive_search(max_rec_len: usize) -> AdaptiveSearchOptions {
        AdaptiveSearchOptions {
            min_rec_len: 1,
            max_rec_len,
            min_var_deg: 0,
            max_var_deg: 0,
            min_idx_deg: 0,
            max_idx_deg: 0,
            min_diff_deg: 0,
            max_diff_deg: 0,
            try_denominator: false,
            try_inhomogeneous: false,
            try_alternating_sign: false,
            ..Default::default()
        }
    }

    #[test]
    fn adaptive_budget_zero_stops_before_the_first_candidate() {
        let polys = vec![vec![1], vec![2], vec![4], vec![8], vec![16]];
        let outcome = find_recurrence_adaptive_with_budget(
            &polys,
            &fixed_adaptive_search(2),
            AdaptiveSearchBudget::limited(0),
        );

        let AdaptiveSearchOutcome::BudgetExhausted(summary) = outcome else {
            panic!("zero budget must stop a nonempty search");
        };
        assert_eq!(summary.diagnostics.considered_candidates, 0);
        assert_eq!(summary.candidates_tried, 0);
    }

    #[test]
    fn adaptive_budget_allows_success_on_the_exact_boundary() {
        let polys = vec![vec![1], vec![2], vec![4], vec![8], vec![16]];
        let outcome = find_recurrence_adaptive_with_budget(
            &polys,
            &fixed_adaptive_search(2),
            AdaptiveSearchBudget::limited(1),
        );

        let AdaptiveSearchOutcome::Found(result) = outcome else {
            panic!("the first candidate should be allowed to succeed");
        };
        assert_eq!(result.recurrence.to_string(), "P(n) = 2 P(n-1)");
        assert_eq!(result.diagnostics.considered_candidates, 1);
        assert_eq!(result.candidates_tried, 1);
    }

    #[test]
    fn adaptive_budget_distinguishes_exhaustion_from_complete_failure() {
        let polys = vec![vec![1], vec![1], vec![2], vec![3], vec![5]];

        let complete = find_recurrence_adaptive_with_budget(
            &polys,
            &fixed_adaptive_search(1),
            AdaptiveSearchBudget::limited(1),
        );
        let AdaptiveSearchOutcome::NoRecurrence(summary) = complete else {
            panic!("one allowed candidate should exhaust a one-candidate search space");
        };
        assert_eq!(summary.diagnostics.considered_candidates, 1);

        let truncated = find_recurrence_adaptive_with_budget(
            &polys,
            &fixed_adaptive_search(2),
            AdaptiveSearchBudget::limited(1),
        );
        let AdaptiveSearchOutcome::BudgetExhausted(summary) = truncated else {
            panic!("a second uninspected candidate must produce budget exhaustion");
        };
        assert_eq!(summary.diagnostics.considered_candidates, 1);
    }

    #[test]
    fn unbounded_budget_preserves_the_existing_search_result() {
        let polys = vec![vec![1], vec![2], vec![4], vec![8], vec![16]];
        let search = fixed_adaptive_search(2);
        let legacy = find_recurrence_adaptive(&polys, &search).unwrap();
        let outcome = find_recurrence_adaptive_with_budget(
            &polys,
            &search,
            AdaptiveSearchBudget::unbounded(),
        );
        let AdaptiveSearchOutcome::Found(unbounded) = outcome else {
            panic!("unbounded outcome API must find the legacy recurrence");
        };

        assert_eq!(
            unbounded.recurrence.to_string(),
            legacy.recurrence.to_string()
        );
        assert_eq!(unbounded.candidates_tried, legacy.candidates_tried);
        assert_eq!(
            unbounded.diagnostics.considered_candidates,
            legacy.diagnostics.considered_candidates
        );
    }
}
