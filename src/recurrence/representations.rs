//! A bounded menu of lag/derivative representations using the existing fitter.
use super::*;
use num_traits::Zero;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Termination of the representative search, not a completeness claim about
/// the annihilating ideal or even all solutions of a finite coefficient ansatz.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum RepresentationSearchStatus {
    /// The search ran to completion: every requested cell was solved or all of its
    /// configurations were inspected.
    Complete,
    /// The global candidate budget ran out while configurations remained uninspected;
    /// unresolved cells are reported as `incomplete`.
    BudgetExhausted,
}

/// Actual support of a candidate; degrees exclude zero coefficient slots.
/// The frontier compares every field componentwise, with `false < true`.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct RepresentationComplexity {
    /// Largest lag `r` (the term `P_{n-r}`) of a term with a nonzero coefficient. This is the
    /// actual lag, which can be smaller than the requested bound.
    pub rec_len: usize,
    /// Largest derivative order of a term with a nonzero coefficient.
    pub diff_deg: usize,
    /// Largest power of `t` among the nonzero coefficients of the right-hand-side terms and
    /// of the inhomogeneous term.
    pub var_deg: usize,
    /// Largest power of the index `n` among the nonzero coefficients of the right-hand-side
    /// terms and of the inhomogeneous term.
    pub idx_deg: usize,
    /// Degree in `t` of the left-hand factor `f(n,t)`; zero when there is no denominator.
    pub denom_var_deg: usize,
    /// Degree in `n` of the left-hand factor `f(n,t)`; zero when there is no denominator.
    pub denom_idx_deg: usize,
    /// Whether a nonzero additive term `g(n,t)` is present.
    pub inhomogeneous: bool,
    /// Whether some term carries the factor `(-1)^n`.
    pub alternating_sign: bool,
    /// Number of nonzero scalar coefficients, including the denominator.
    pub nonzero_coefficients: usize,
}

/// One exact fit to supplied rows. This is never an all-rank proof.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RecurrenceRepresentation {
    /// The recurrence as text, `f(n,t) P(n) = ...`.
    pub recurrence: String,
    /// The recurrence as a LaTeX expression.
    pub latex: String,
    /// Reusable existing v1 recurrence schema, with sufficient initial rows.
    pub recurrence_json: RecurrenceJson,
    /// Actual support of the recurrence.
    pub complexity: RepresentationComplexity,
    /// Always `exact_fit`; checking a finite input cannot certify an identity.
    pub evidence: String,
    /// First checked target in the renumbered sequence (whose first row is 1).
    pub first_checked_index: usize,
    /// Last checked target index in the same renumbering (the number of rows after
    /// `skip_prefix`).
    pub last_checked_index: usize,
    /// First checked target as a 1-based position in the original input.
    pub first_checked_input_row: usize,
    /// Number of leading rows, counted after `skip_prefix`, used to fit the coefficients.
    pub fit_polynomials: usize,
    /// Number of trailing rows reserved for held-out verification; zero when verification
    /// is disabled.
    pub verification_polynomials: usize,
    /// Supplied target indices at which the denominator is the zero polynomial.
    /// The identity can still hold there, but it cannot generate that row.
    pub singular_checked_indices: Vec<usize>,
}

/// Coverage of one requested pair of length/order bounds.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RepresentationCoverage {
    /// Requested lag bound of this cell: candidates use lags `1..=rec_len`.
    pub rec_len: usize,
    /// Requested derivative-order bound of this cell: candidates use orders `0..=diff_deg`.
    pub diff_deg: usize,
    /// Index into `representations`, when this cell produced a candidate.
    pub representation: Option<usize>,
    /// `found`, `not_found`, `incomplete`, or `insufficient_rows`.
    /// `not_found` refers to the adaptive fitter, not a nonexistence proof.
    pub status: String,
    /// Number of candidate configurations of this cell that were inspected.
    pub configurations_considered: usize,
}

/// Serializable menu shared by the Rust, CLI, MCP and browser APIs.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RecurrenceRepresentations {
    /// Schema identifier, `polytool.recurrence-representations.v1`.
    pub schema: String,
    /// Whether the search completed or exhausted its budget.
    pub status: RepresentationSearchStatus,
    /// The search bounds that were used.
    pub search: AdaptiveSearchOptions,
    /// Global limit on inspected configurations; `None` means unbounded.
    pub max_candidates: Option<usize>,
    /// Number of input rows supplied, before `skip_prefix` is applied.
    pub source_rows: usize,
    /// Total number of configurations inspected across all cells; this is the quantity
    /// bounded by `max_candidates`.
    pub configurations_considered: usize,
    /// Total number of configurations in the searched space.
    pub configurations_available: usize,
    /// Number of inspected configurations skipped because no fitting prefix of the rows
    /// had enough equations.
    pub insufficient_fit_rows: usize,
    /// Retained candidates: at most one per cell, deduplicated when formula and checked
    /// domain agree.
    pub representations: Vec<RecurrenceRepresentation>,
    /// Indices of nondominated candidates among those returned. Starting index
    /// is also compared, so a later tail does not dominate an earlier identity.
    pub frontier: Vec<usize>,
    /// One entry per requested `(rec_len, diff_deg)` cell, in ascending lexicographic order.
    pub coverage: Vec<RepresentationCoverage>,
}

fn degrees(poly: &BivarPoly) -> (usize, usize, usize) {
    let (mut n, mut x, mut count) = (0, 0, 0);
    for (i, row) in poly.coeffs.iter().enumerate() {
        for (j, c) in row.iter().enumerate() {
            if !c.is_zero() {
                n = n.max(i);
                x = x.max(j);
                count += 1;
            }
        }
    }
    (n, x, count)
}

fn complexity(rec: &Recurrence) -> RepresentationComplexity {
    let mut result = RepresentationComplexity {
        rec_len: 0,
        diff_deg: 0,
        var_deg: 0,
        idx_deg: 0,
        denom_var_deg: 0,
        denom_idx_deg: 0,
        inhomogeneous: false,
        alternating_sign: false,
        nonzero_coefficients: 0,
    };
    for term in &rec.terms {
        let (n, x, count) = degrees(&term.coeff);
        if count == 0 {
            continue;
        }
        result.rec_len = result.rec_len.max(term.offset);
        result.diff_deg = result.diff_deg.max(term.deriv_order);
        result.var_deg = result.var_deg.max(x);
        result.idx_deg = result.idx_deg.max(n);
        result.alternating_sign |= term.sign != RecurrenceSign::None;
        result.nonzero_coefficients += count;
    }
    if let Some(p) = &rec.inhomogeneous {
        let (n, x, count) = degrees(p);
        result.inhomogeneous = count > 0;
        result.var_deg = result.var_deg.max(x);
        result.idx_deg = result.idx_deg.max(n);
        result.nonzero_coefficients += count;
    }
    if let Some(p) = &rec.denominator {
        let (n, x, count) = degrees(p);
        result.denom_var_deg = x;
        result.denom_idx_deg = n;
        result.nonzero_coefficients += count;
    } else {
        result.nonzero_coefficients += 1;
    }
    result
}

fn dominates(a: &RecurrenceRepresentation, b: &RecurrenceRepresentation) -> bool {
    let coordinates = |r: &RecurrenceRepresentation| {
        let c = &r.complexity;
        [
            c.rec_len,
            c.diff_deg,
            c.var_deg,
            c.idx_deg,
            c.denom_var_deg,
            c.denom_idx_deg,
            usize::from(c.inhomogeneous),
            usize::from(c.alternating_sign),
            c.nonzero_coefficients,
            r.first_checked_index,
        ]
    };
    let aa = coordinates(a);
    let bb = coordinates(b);
    aa != bb && aa.iter().zip(bb).all(|(x, y)| *x <= y)
        // An unusable row-generation formula cannot suppress a usable one.
        && a.singular_checked_indices.iter().all(|n| b.singular_checked_indices.contains(n))
}

// The scalar fitter fixes the constant denominator coefficient to 1. Reuse
// the vector fitter's projective denominator pivots for the missing charts,
// rather than introducing another exact linear-algebra implementation.
fn denominator_fallback(
    ctx: &AdaptiveCandidateContext<'_>,
    opts: &RecurrenceOptions,
    tried: &mut usize,
) -> Option<FoundAdaptiveCandidate> {
    if opts.denom_idx_deg == 0 && opts.denom_var_deg == 0 {
        return None;
    }
    let fit_len = fitting_polynomial_count_with_prefix_degrees(
        ctx.polys,
        opts,
        ctx.search,
        ctx.prefix_max_degrees,
    )?;
    if count_equations_from_prefix_degrees(fit_len, opts, ctx.prefix_max_degrees)
        < count_unknowns(opts).saturating_add(ctx.search.min_margin)
    {
        return None;
    }
    let signs = if opts.alternating_sign { 2 } else { 1 };
    let mut sources = Vec::new();
    let mut targets = Vec::new();
    let indices: Vec<_> = (opts.rec_len + 1..=ctx.polys.len()).collect();
    for &nn in &indices {
        let mut source = Vec::new();
        for lag in 1..=opts.rec_len {
            for sign in 0..signs {
                let p = &ctx.polys[nn - 1 - lag];
                source.push(if sign == 1 && nn % 2 == 1 {
                    p.iter().map(|c| -c).collect()
                } else {
                    p.clone()
                });
            }
        }
        sources.push(source);
        targets.push(vec![ctx.polys[nn - 1].clone()]);
    }
    let vector_opts = VectorRecurrenceDenominatorOptions {
        var_deg: opts.var_deg,
        idx_deg: opts.idx_deg,
        diff_deg: opts.diff_deg,
        homogeneous: opts.homogeneous,
        forcing_var_deg: opts.inhomo_var_deg,
        forcing_idx_deg: opts.inhomo_idx_deg,
        held_out_transitions: ctx.polys.len() - fit_len,
        require_unique: false,
        denominator_var_deg: opts.denom_var_deg,
        denominator_idx_deg: opts.denom_idx_deg,
    };
    *tried += 1;
    let fitted = fit_shared_denominator_map_rational(
        &sources,
        &targets,
        &indices,
        fit_len - opts.rec_len,
        &vector_opts,
        0,
    )
    .ok()?;
    let mut terms = Vec::new();
    for (column, operator) in fitted.recurrence.matrix[0].iter().enumerate() {
        for (deriv_order, coeff) in operator.coefficients.iter().enumerate() {
            if !coeff.is_zero() {
                terms.push(RecurrenceTerm {
                    offset: column / signs + 1,
                    deriv_order,
                    sign: RecurrenceSign::from_family_index(column % signs),
                    coeff: coeff.clone(),
                });
            }
        }
    }
    let rec = Recurrence {
        terms,
        denominator: Some(fitted.recurrence.denominator),
        inhomogeneous: fitted.recurrence.forcing.and_then(|mut f| f.pop()),
    };
    // Independently check the scalar conversion, including alternating parity.
    if !recurrence_holds_rational(ctx.polys, &rec, opts.rec_len) {
        return None;
    }
    let counts = candidate_complexity(opts);
    Some(FoundAdaptiveCandidate {
        recurrence: rec,
        opts: opts.clone(),
        num_unknowns: fitted.diagnostics.unknowns,
        weighted_unknowns: counts.weighted_unknowns,
        num_equations: fitted.diagnostics.total_equations,
        fit_polynomials: fit_len,
        verification_polynomials: ctx.polys.len() - fit_len,
        candidates_tried: *tried,
    })
}

/// Find a menu of recurrences trading lag length against derivative order.
///
/// For each `(rec_len, diff_deg)` bound pair, retain the first successful
/// candidate in the existing weighted-complexity ordering. The scalar path
/// shares derivative and modular caches; the vector denominator fallback builds
/// its own per-shape system. The budget is global;
/// configurations in already solved cells are skipped without consuming it.
/// Unlike the historical single-result finder, this API never implicitly
/// enables denominators outside the requested search bounds.
///
/// The result describes exact fits, with held-out verification unless disabled.
/// It is not a basis of the recurrence space, an all-rank certificate, a proof
/// of minimality, or a proof of nonexistence for cells without a candidate.
/// Polynomial common factors are not cancelled. Every retained candidate keeps
/// its original checked starting index even if its actual lag support is smaller.
///
/// With `C` inspected configurations, at most `E` equations, `U` coefficient
/// slots, and `Q` denominator slots, dense elimination costs conservatively
/// `O(C Q E U^2)` rational operations, apart from constructing the systems;
/// modular filtering often avoids this work. Rational bit complexity depends
/// on the input. At most one representation is retained per requested
/// length/order cell.
pub fn find_recurrence_representations_rational(
    input: &[Vec<BigRational>],
    search: &AdaptiveSearchOptions,
    budget: AdaptiveSearchBudget,
) -> Result<RecurrenceRepresentations, String> {
    for (low, high) in [
        (search.min_rec_len.max(1), search.max_rec_len),
        (search.min_diff_deg, search.max_diff_deg),
        (search.min_var_deg, search.max_var_deg),
        (search.min_idx_deg, search.max_idx_deg),
        (search.min_inhomo_var_deg, search.max_inhomo_var_deg),
        (search.min_inhomo_idx_deg, search.max_inhomo_idx_deg),
    ] {
        if low > high {
            return Err(format!("invalid recurrence search range {low}..={high}"));
        }
    }
    let polys = input.get(search.skip_prefix..).unwrap_or(&[]);
    if polys.len() < 3 {
        return Err("need at least 3 polynomials after skip_prefix".into());
    }
    // Bound the report allocation independently of the solver budget.
    let cells = (search.max_rec_len - search.min_rec_len.max(1) + 1)
        .checked_mul(
            search
                .max_diff_deg
                .checked_sub(search.min_diff_deg)
                .and_then(|d| d.checked_add(1))
                .ok_or("search size overflow")?,
        )
        .ok_or("search size overflow")?;
    if cells > 10_000 {
        return Err("at most 10000 length/derivative cells are supported".into());
    }
    let mut coverage = BTreeMap::new();
    for r in search.min_rec_len.max(1)..=search.max_rec_len {
        for d in search.min_diff_deg..=search.max_diff_deg {
            coverage.insert(
                (r, d),
                RepresentationCoverage {
                    rec_len: r,
                    diff_deg: d,
                    representation: None,
                    status: if r >= polys.len() {
                        "insufficient_rows"
                    } else {
                        "not_found"
                    }
                    .into(),
                    configurations_considered: 0,
                },
            );
        }
    }
    let candidates = CandidateIterator::new(polys.len(), search);
    let available = candidates.total_candidates();
    let mut diagnostics = AdaptiveSearchDiagnostics::default();
    let mut representations = Vec::new();
    let mut tried = 0;
    let mut considered = 0;
    let mut status = RepresentationSearchStatus::Complete;
    // A zero budget should not allocate derivative/modular caches.
    let derivs = if budget.max_candidates == Some(0) {
        Vec::new()
    } else {
        rational_derivatives_up_to(polys, search.max_diff_deg)
    };
    let prefix_degrees = prefix_max_degrees_rational(polys);
    let cache = if search.modular_prefilter && budget.max_candidates != Some(0) {
        Some(ModularPrefilterCache::new(polys, search.max_diff_deg))
    } else {
        None
    };
    let context = AdaptiveCandidateContext {
        polys,
        derivs: &derivs,
        prefix_max_degrees: &prefix_degrees,
        modular_cache: cache.as_ref(),
        search,
        fit_search: search,
        pass_label: None,
    };
    for opts in candidates {
        let cell = coverage
            .get_mut(&(opts.rec_len, opts.diff_deg))
            .expect("requested cell");
        if cell.representation.is_some() {
            continue;
        }
        if !budget.allows_candidate(considered) {
            status = RepresentationSearchStatus::BudgetExhausted;
            break;
        }
        considered += 1;
        cell.configurations_considered += 1;
        let found = evaluate_adaptive_candidate(&context, &opts, &mut diagnostics, &mut tried)
            .or_else(|| denominator_fallback(&context, &opts, &mut tried));
        if let Some(found) = found {
            let first = opts.rec_len + 1;
            let metadata = RecurrenceJsonSearch {
                recurrence_text: found.recurrence.to_string(),
                source_rows: input.len(),
                skip_prefix: search.skip_prefix,
                unknowns: found.num_unknowns,
                weighted_unknowns: found.weighted_unknowns,
                equations: found.num_equations,
                fit_polynomials: found.fit_polynomials,
                verification_polynomials: found.verification_polynomials,
                candidates_tried: found.candidates_tried,
                options: (&opts).into(),
            };
            // Keep *all* rows before the checked domain, not only max_offset().
            let json = RecurrenceJson::from_recurrence_rational(
                &found.recurrence,
                1,
                &polys[..opts.rec_len],
                Some(metadata),
            );
            let singular = (first..=polys.len())
                .filter(|&n| {
                    found
                        .recurrence
                        .denominator
                        .as_ref()
                        .is_some_and(|p| poly_is_zero_rational(&bivar_eval_n(p, n)))
                })
                .collect();
            let item = RecurrenceRepresentation {
                recurrence: found.recurrence.to_string(),
                latex: found.recurrence.to_latex(),
                recurrence_json: json,
                complexity: complexity(&found.recurrence),
                evidence: "exact_fit".into(),
                first_checked_index: first,
                last_checked_index: polys.len(),
                first_checked_input_row: search.skip_prefix + first,
                fit_polynomials: found.fit_polynomials,
                verification_polynomials: found.verification_polynomials,
                singular_checked_indices: singular,
            };
            // Exact JSON equality includes initial rows and checked start.
            // Same formula at a later start remains a separately qualified fit.
            let id = representations
                .iter()
                .position(|r: &RecurrenceRepresentation| {
                    r.recurrence_json.recurrence == item.recurrence_json.recurrence
                        && r.first_checked_index == item.first_checked_index
                })
                .unwrap_or_else(|| {
                    representations.push(item);
                    representations.len() - 1
                });
            cell.representation = Some(id);
            cell.status = "found".into();
        }
    }
    if status == RepresentationSearchStatus::BudgetExhausted {
        for cell in coverage.values_mut() {
            if cell.status == "not_found" {
                cell.status = "incomplete".into();
            }
        }
    }
    let frontier = (0..representations.len())
        .filter(|&i| {
            !representations
                .iter()
                .any(|other| dominates(other, &representations[i]))
        })
        .collect();
    Ok(RecurrenceRepresentations {
        schema: "polytool.recurrence-representations.v1".into(),
        status,
        search: search.clone(),
        max_candidates: budget.max_candidates,
        source_rows: input.len(),
        configurations_considered: considered,
        configurations_available: available,
        insufficient_fit_rows: diagnostics.insufficient_fit_rows,
        representations,
        frontier,
        coverage: coverage.into_values().collect(),
    })
}

/// Integer convenience wrapper for [`find_recurrence_representations_rational`].
pub fn find_recurrence_representations(
    polys: &[Vec<i64>],
    search: &AdaptiveSearchOptions,
    budget: AdaptiveSearchBudget,
) -> Result<RecurrenceRepresentations, String> {
    find_recurrence_representations_rational(&i64_polys_to_rational(polys), search, budget)
}
