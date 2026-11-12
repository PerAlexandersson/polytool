//! Row-hypergeometric families: rows whose coefficient ratio
//! `[x^k] P_n / [x^(k-1)] P_n` is a rational function of `(n, k)` that
//! factors into linear factors in `k`, proved real-rooted for every `n` by a
//! product of multiplier sequences applied to `(1 + x)^M(n)`
//! ([`check_hypergeometric_row_theorem`]).
//!
//! # The row relation
//!
//! From the rows `P_n`, `n >= s` (`s` the first base where the rule holds
//! for the generated rows), a relation
//!
//! ```text
//! p(n, k) [x^k] P_n + q(n, k) [x^(k-1)] P_n = 0      for every integer k
//! ```
//!
//! with polynomials `p, q` of total degree at most 4 is fitted exactly
//! (common linear factors in `k` cancelled).  In operator form it is
//! `H(n) P_n = 0` with `H(n) = p(n, theta) + x q(n, theta + 1)`,
//! `theta = x D`.  The fit is only a candidate.  `H(n) P_n = 0` is proved
//! for every `n >= s` by an invariant row-annihilator system of the rule
//! (the certificate of the Hadamard-factor lifts, module documentation of
//! `hadamard_factor`, "Row annihilators"): exact identities
//! `nu_l(n, x) N(A^(l)) = sum_j B_(l,j)(n) o A^(j)` under the rule's
//! companion step, the generators checked on the rows at the start base and
//! after every restart base, and `mu(n, x) H(n) = sum_l Z_l(n) o A^(l)(n)`
//! with `mu` nonzero; at the integer bases where `mu(n, .)` vanishes,
//! `H(n) P_n = 0` is checked exactly.
//!
//! # From the relation to a multiplier sequence
//!
//! Write `p = c_p(n) prod_i (k + u_i(n))` and `q = c_q(n) prod_j (k + w_j(n))`
//! with `u_i, w_j` affine in `n` (other shapes are not applicable).  The
//! coefficients are read in one of two coordinates, `U_j = [x^(j + v(n))] P_n`
//! (the identity, `v(n)` the valuation of the rows, affine with integer
//! coefficients) or `U_j = [x^(d(n) - j)] P_n` (the reversal at the degree
//! `d(n)`, affine with integer coefficients).  In both,
//! `P(n, j) U_j + Q(n, j) U_(j-1) = 0` for every integer `j`, with
//! `R(n, j) = -Q / P = C(n) prod (j + a_l(n)) / prod (j + b_m(n))`
//! (identity: `a = w + v`, `b = u + v`, `C = -c_q / c_p`; reversal:
//! `a = -d - 1 - u`, `b = -d - 1 - w`, `C = -(-1)^(|p| + |q|) c_p / c_q`).
//! The support of the rows is proved as follows:
//!
//! - `[x^k] P_n = 0` for `k < v(n)`: by induction on the rule when every
//!   term `c x^e D^d P_(n-r)` has `e - d - v_1 r >= 0` (rows before the
//!   induction checked exactly), or, for a constant valuation `v`, by the
//!   relation with `p(n, k) != 0` for `0 <= k < v` (then each coefficient
//!   below `x^v` vanishes, from `[x^(-1)] P_n = 0` upwards);
//! - (reversal) `[x^k] P_n = 0` for `k > d(n)`: by induction on the rule
//!   when every term has `e - d - d_1 r <= 0`, or by the relation: `q` has
//!   the factor `k - d(n) - 1`, `c_p(n) != 0`, and every root `-u(n)` of
//!   `p` satisfies `u(n) + d(n) + 1 > 0` or lies above the crude degree
//!   bound `D(n) = D_0 + beta n` of the rule (`beta = max (e - d) / r`), so
//!   `p(n, k) != 0` for `d(n) < k <= D(n)`.
//!
//! The decided conditions start at a structure start at most
//! `MAX_STRUCTURE_DELAY` rows after the first base where the relation holds;
//! every row before the induction start is checked exactly.
//!
//! The factors are matched with multiplier sequences, every condition
//! decided for all `n` by [`crate::decide_index_rational_function_sign`]:
//!
//! - the **base**: a numerator factor `j - M(n) - 1` (`M` affine with integer
//!   coefficients, `M(n) >= 0`; for the reversal `M(n) <= d(n) - v`) and the
//!   denominator factor `j`, giving `binom(M, j)`;
//! - **complements**: further numerator factors `j - M'(n) - 1` with integer
//!   coefficients and `M' >= M`, giving `k -> M'! / (M' - k)!`, a multiplier
//!   sequence by the Pólya--Schur theorem (its exponential generating
//!   function is `(1 + x)^M'`);
//! - **pairs**: a numerator factor `j + a` with a denominator factor `j + b`,
//!   `a - b` affine with integer coefficients, `a - b >= 0` and `b > -1`,
//!   giving the polynomial `k -> (k + b + 1) ... (k + a)` with roots `<= 0`
//!   (Laguerre's theorem);
//! - **reciprocals**: the other denominator factors `j + b`, `b > -1`,
//!   giving `k -> 1 / (b + 1)_k`;
//! - the **scale** `C_eff(n) = C(n) (-1)^c > 0` (`c` complements, the base
//!   included), giving `k -> C_eff^k`.
//!
//! Every denominator factor is positive at `j >= 1`, so `P(n, j) != 0`
//! there and, for `j >= 0`, `U_j = U_0 prod_(t <= j) R(n, t)`, which is
//! `U_0 binom(M, j) mu_j` with `mu` the product of the sequences above (zero
//! for `j > M`, where the base factor vanishes).  So the polynomial
//! `sum_j U_j y^j` is `U_0` times the image of `(1 + y)^M` under the
//! nonnegative multiplier sequence `mu`, hence real-rooted with every root
//! `< 0` when `U_0 != 0`.  `U_0` is `[x^v] P_n` (identity) or `[x^d(n)] P_n`
//! (reversal); the rule gives a scalar recurrence
//! `f(n) U_0(n) = sum_r sigma_r(n) U_0(n - r)` (each coefficient of the rule
//! that reaches `U_0` is `U_0(n - r)` times a fixed number of ratio factors
//! `R(n - r, t)`), and `U_0(n) != 0` with one sign follows from
//! `sigma_r(n) / f(n) > 0` for every lag that occurs and nonzero seeds of one
//! sign (or from a hypergeometric solution of that recurrence); when neither
//! applies, the shared module
//! ([`crate::recurrence::certify_leading_coefficient`]) may prove it for the
//! scalar recurrence read as a family of constant rows (for example a
//! growth-ratio invariant for `U_0(n) = 3 U_0(n - 1) - 2 U_0(n - 2)`).  The
//! rows before the induction start are checked exactly.
//!
//! The conclusion: every `P_n` is nonzero and real-rooted with every root in
//! `(-inf, 0)`, or in `(-inf, 0]` when the rows are divisible by `x` (or the
//! reversal does not show `[x^v] P_n != 0`).  No interlacing is claimed.

use super::*;
use crate::factorization::irreducible_factorization_monic_bivariate_with_budget;
use crate::recurrence::implicit::algebra::{
    content_and_primitive, from_theta, theta_add, Biv, DOp, ThetaOp,
};
use crate::recurrence::lift::hadamard_factor::annihilator::vanishing_bases;
use crate::recurrence::lift::hadamard_factor::annihilator_system::find_invariant_annihilator_system;
use crate::recurrence::lift::hadamard_factor::certify::{
    check_system_rows, MAX_HADAMARD_OPERATOR_TERMS,
};
use crate::recurrence::lift::hadamard_factor::prepare::prepare_lift_rule;
use num_traits::ToPrimitive;
use ore_algebra::exact_kernel_vectors_by_group;

/// Theorem statement and references of [`check_hypergeometric_row_theorem`].
const REFERENCE: &str = "Row-hypergeometric rows: an exact row relation p(n, k) [x^k] P_n + \
     q(n, k) [x^(k-1)] P_n = 0, proved for every n by an invariant row-annihilator system of the \
     rule, writes every row as U_0(n) x^v (or its reversal) times the image of (1 + x)^M(n) under \
     a product of nonnegative multiplier sequences: C^k (C > 0); k -> M'!/(M' - k)!, whose \
     exponential generating function is (1 + x)^M' (Polya--Schur, Theorem 7.1 (3) => (1) of \
     P. Branden, Unimodality, log-concavity, real-rootedness and beyond, arXiv:1410.6601: a \
     sequence is a multiplier sequence when G(x) = sum_k lambda_k x^k / k! or G(-x) can be \
     written C x^m e^(a x) prod_k (1 + alpha_k x) with a >= 0, alpha_k >= 0); k -> phi(k) for a \
     polynomial phi with only real roots, all <= 0 (Laguerre's theorem, RealRooted \
     isPFMultiplierSequence_eval_of_roots_nonpos); k -> 1/(alpha)_k for alpha > 0 (RealRooted \
     isPFMultiplierSequence_inv_ascPochhammer, via the generalized Laguerre polynomials); \
     products of PF multiplier sequences are PF multiplier sequences (RealRooted \
     IsPFMultiplierSequence.mul) and map real-rooted polynomials with nonnegative coefficients \
     to such polynomials (IsPFMultiplierSequence.diagonalOperator_isPF)";

/// Largest total degree in `(n, k)` of a fitted row relation.
const MAX_RELATION_DEGREE: usize = 4;

/// Rows (from the first base where the rule holds) used to fit the
/// relation.
const FIT_ROWS: usize = 14;

/// Largest number of rows by which the start of the decided conditions may
/// move past the first base (the rows before it are checked exactly).
const MAX_STRUCTURE_DELAY: usize = 4;

/// Largest degree of the numerator and denominator of `rho(n)` in
/// `U_0(n) = rho(n) U_0(n - 1)`.
const MAX_FIRST_COEFFICIENT_DEGREE: usize = 3;

/// Recombination budget of the bivariate factorizations.
const FACTORIZATION_BUDGET: usize = 256;

/// Check the row-hypergeometric theorem (module documentation) for the
/// family of `rule`, `first_index` and `seeds`.
///
/// Fits an exact row relation `p(n, k) [x^k] P_n + q(n, k) [x^(k-1)] P_n = 0`
/// (total degree at most 4), proves it for every `n` by an invariant
/// row-annihilator system of the rule, and decides, for all `n` with the
/// exact sign checker, that the ratio `-q/p` is the ratio of `U_0(n)` times
/// the image of `(1 + x)^M(n)` under a product of nonnegative multiplier
/// sequences (Pólya--Schur, Laguerre, `1/(alpha)_k`), in the identity
/// coordinates or after the reversal at the degree.  Covers, for example,
/// `sum_k binom(n, k) (2n + 1)_k / (n + 2)_k x^k` (OEIS A330798) and
/// `sum_k binom(n, k) (k + 1)(k + 2)(k + 3) x^k` (A124848).  A failed
/// hypothesis says nothing against real-rootedness.
///
/// # Errors
///
/// Errors parsing or normalizing the recurrence and seeds.
pub fn check_hypergeometric_row_theorem(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    check_real_rootedness_theorem(
        RealRootednessTheorem::HypergeometricRow,
        rule,
        first_index,
        seeds,
        options,
    )
}

/// The registry entry of the row-hypergeometric theorem ([`super::registry::ROUTES`]).
pub(super) const ROUTE: RoutePolicy = RoutePolicy {
    reference: REFERENCE,
    phase: DiscoveryPhase::Direct,
    after_direct_eligible: None,
    replay: ParentReplayPolicy::RecordedRoute,
    discovery_replay: DiscoveryReplayPolicy::Direct,
    check_at_supplied_seeds: hypergeometric_row_report,
    start_policy: StartPolicy::MovesStart,
    is_preferred: any_proof,
    certificate_name: None,
    discovery: DiscoveryPolicy::Always,
    in_residue_class_menu: true,
    in_index_section_menu: true,
    in_degree_law_reversal_menu: true,
    parent: ParentPolicy {
        stage: ParentStage::ReplacesEarly,
        region: ParentRegion::Unrestricted,
        record_only_if_applicable: true,
    },
};

/// The theorem at the supplied seeds.
fn hypergeometric_row_report(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    let mut report = RealRootednessTheoremReport::new(
        RealRootednessTheorem::HypergeometricRow,
        rule,
        first_index,
        seeds,
        options,
    );
    let shape = match recurrence_shape(rule) {
        Ok(shape) => shape,
        Err(reason) => return Ok(report.not_applicable(reason)),
    };
    if seeds.is_empty() || shape.terms.is_empty() {
        return Ok(report.not_applicable("the rule needs seed rows and a lag term"));
    }
    let total = seeds.len() + FIT_ROWS;
    let rows: Vec<Vec<Q>> = match rule.generate_rows_rational(seeds, first_index, total) {
        Ok(rows) => rows.iter().map(|row| vec_poly::trim(row)).collect(),
        Err(error) => {
            return Ok(report.not_applicable(format!("the rule cannot generate rows: {error}")))
        }
    };
    // The prepared rule; when the supplied seeds do not all follow the rule
    // (more seeds than its order), the family is prepared from its last
    // identifying seeds on, and the earlier rows are checked exactly below.
    let prepared = match prepare_lift_rule(rule, &rows, first_index, "input") {
        Ok(prepared) => prepared,
        Err(error) => {
            let skip = seeds.len().saturating_sub(rule.max_offset());
            match (skip > 0)
                .then(|| prepare_lift_rule(rule, &rows[skip..], first_index + skip, "input"))
            {
                Some(Ok(prepared)) => prepared,
                _ => return Ok(report.not_applicable(error.to_string())),
            }
        }
    };
    let rule_start = first_index + seeds.len();
    let base = prepared.valid_from_base();
    // The first start (up to MAX_STRUCTURE_DELAY rows after the first base)
    // with a fitted relation of hypergeometric shape whose reduced form
    // kills every generated row from there on.
    let mut shape_failure = None;
    let attempt = |from: usize| -> Option<(LinearFactors, LinearFactors, DOp)> {
        let (p, q) = fit_row_relation(&rows, first_index, from)?;
        let (mut p, mut q) = (linear_factors(&p).ok()?, linear_factors(&q).ok()?);
        // Cancel the linear factors common to p and q.
        let mut kept = Vec::new();
        for root in std::mem::take(&mut p.roots) {
            if let Some(position) = q.roots.iter().position(|r| *r == root) {
                q.roots.remove(position);
            } else {
                kept.push(root);
            }
        }
        p.roots = kept;
        let operator = relation_operator(&p.biv(), &q.biv());
        rows.iter()
            .enumerate()
            .map(|(i, row)| (first_index + i, row))
            .filter(|(n, _)| *n >= from)
            .all(|(n, row)| operator.apply_row(row, n).iter().all(Zero::is_zero))
            .then_some((p, q, operator))
    };
    let mut found = None;
    for from in base..=base + MAX_STRUCTURE_DELAY {
        if let Some(result) = attempt(from) {
            found = Some((from, result));
            break;
        }
        if shape_failure.is_none() {
            if let Some((p, q)) = fit_row_relation(&rows, first_index, from) {
                if let Err(reason) = linear_factors(&p).and_then(|_| linear_factors(&q)) {
                    shape_failure = Some(reason);
                }
            }
        }
    }
    let Some((from, (p, q, operator))) = found else {
        return Ok(report.not_applicable(match shape_failure {
            Some(reason) => format!("the row relation is not of hypergeometric shape: {reason}"),
            None => format!(
                "no row relation p(n, k) [x^k] P_n + q(n, k) [x^(k-1)] P_n = 0 of total degree \
                 at most {MAX_RELATION_DEGREE} with linear factors in k holds on the rows from \
                 n = {base} + delay on, delay <= {MAX_STRUCTURE_DELAY}"
            ),
        }));
    };
    // Valuation and degree of the rows, affine in n from some row on.
    let fitted: Vec<(usize, &Vec<Q>)> = rows
        .iter()
        .enumerate()
        .map(|(i, row)| (first_index + i, row))
        .filter(|(n, _)| *n >= from)
        .collect();
    if fitted.iter().any(|(_, row)| vec_poly::is_zero(row)) {
        return Ok(report.not_applicable("a row after the first base is zero"));
    }
    let as_q = |value: usize| Q::from_integer(BigInt::from(value));
    let valuations: Vec<(usize, Q)> = fitted
        .iter()
        .map(|(n, row)| (*n, as_q(row.iter().position(|c| !c.is_zero()).unwrap_or(0))))
        .collect();
    let degrees: Vec<(usize, Q)> = fitted
        .iter()
        .map(|(n, row)| (*n, as_q(row.len() - 1)))
        .collect();
    let Some((shape_from, valuation, degree)) = (0..=MAX_STRUCTURE_DELAY).find_map(|delay| {
        let valuation = fit_affine(&valuations[delay.min(valuations.len())..])?;
        Some((
            from + delay,
            valuation,
            fit_affine(&degrees[delay.min(degrees.len())..]),
        ))
    }) else {
        return Ok(report.not_applicable(
            "the valuation of the rows is not an affine function of n with integer coefficients",
        ));
    };
    report.normal_form = format!(
        "row relation p(n, k) [x^k] P_n + q(n, k) [x^(k-1)] P_n = 0, n >= {from}, with p = {}, \
         q = {}; valuation v(n) = {}{}",
        p.text("k"),
        q.text("k"),
        affine_text(&valuation),
        degree
            .as_ref()
            .map(|d| format!(", degree d(n) = {}", affine_text(d)))
            .unwrap_or_default()
    );
    report.lean_declarations = vec![
        "RealRooted.isPFMultiplierSequence_eval_of_roots_nonpos".into(),
        "RealRooted.isPFMultiplierSequence_inv_ascPochhammer".into(),
        "RealRooted.IsPFMultiplierSequence.mul".into(),
        "RealRooted.IsPFMultiplierSequence.diagonalOperator_isPF".into(),
    ];

    // The relation for every n >= from.
    report.hypotheses.push(relation_hypothesis(
        &operator,
        &prepared,
        from,
        &rows,
        first_index,
    ));

    let mut chosen = None;
    let mut fallback = None;
    'search: for coordinates in [Coordinates::Identity, Coordinates::Reversal] {
        let ratio = match (coordinates, &degree) {
            (Coordinates::Identity, _) => Ratio {
                scale: IndexFunction::new(vec_poly::neg(&q.scalar), p.scalar.clone()),
                numerators: q.roots.iter().map(|w| affine_add(w, &valuation)).collect(),
                denominators: p.roots.iter().map(|u| affine_add(u, &valuation)).collect(),
            },
            (Coordinates::Reversal, Some(degree)) => {
                let shifted = |r: &Affine| {
                    affine_neg(&affine_add(&affine_add(degree, &affine_constant(1)), r))
                };
                let sign = if (p.roots.len() + q.roots.len()) % 2 == 0 {
                    -Q::one()
                } else {
                    Q::one()
                };
                Ratio {
                    scale: IndexFunction::new(vec_poly::scale(&p.scalar, &sign), q.scalar.clone()),
                    numerators: p.roots.iter().map(shifted).collect(),
                    denominators: q.roots.iter().map(shifted).collect(),
                }
            }
            (Coordinates::Reversal, None) => continue,
        };
        let upper = (coordinates == Coordinates::Reversal)
            .then(|| {
                degree
                    .as_ref()
                    .map(|d| affine_add(d, &affine_neg(&valuation)))
            })
            .flatten();
        for structure_start in shape_from.max(from)..=from + MAX_STRUCTURE_DELAY {
            let Some(matching) = find_matching(&ratio, structure_start, upper.as_ref()) else {
                continue;
            };
            let support = support_hypotheses(
                &shape,
                &p,
                &q,
                &valuation,
                (coordinates == Coordinates::Reversal)
                    .then_some(degree.as_ref())
                    .flatten(),
                (&rows, first_index, rule_start),
                structure_start,
            );
            let verified = support
                .iter()
                .all(|h| h.status == LiftHypothesisStatus::Verified);
            if verified || fallback.is_none() {
                fallback = Some((
                    coordinates,
                    ratio.clone(),
                    matching,
                    upper.clone(),
                    structure_start,
                    support,
                ));
            }
            if verified {
                chosen = fallback.take();
                break 'search;
            }
        }
    }
    let Some((coordinates, ratio, matching, upper, s, support)) = chosen.or(fallback) else {
        report.hypotheses.push(unresolved_hypothesis(
            "multiplier_matching",
            "the ratio factors match a product of multiplier sequences applied to (1 + x)^M",
            format!(
                "no matching of the ratio factors (identity: R(n, j) = -q(n, j + v(n)) / p(n, j \
                 + v(n))) with a base binom(M, j), complements, pairs and reciprocal rising \
                 factorials holds for every n >= {from} + delay, delay <= \
                 {MAX_STRUCTURE_DELAY}{}",
                if degree.is_some() {
                    ", also after the reversal at the degree"
                } else {
                    ""
                }
            ),
        ));
        return Ok(report.finish(RealRootednessTheoremConclusion {
            first_index,
            zero_residue_classes: None,
            interlacing_from: None,
            interlacing: None,
            simple_roots: false,
            root_region: None,
            statement: String::new(),
        }));
    };

    report.hypotheses.extend(support);
    report
        .hypotheses
        .extend(matching_hypotheses(&ratio, &matching, s, upper.as_ref()));

    // The first coefficient U_0(n): its scalar recurrence
    // f(n) U_0(n) = sum_r sigma_r(n) U_0(n - r) from the rule.
    let anchor: Affine = match (coordinates, &degree) {
        (Coordinates::Reversal, Some(degree)) => degree.clone(),
        _ => valuation.clone(),
    };
    let slope = match coordinates {
        Coordinates::Identity => valuation.1.clone(),
        Coordinates::Reversal => -anchor.1.clone(),
    };
    let mut lag_coefficients: BTreeMap<usize, IndexFunction> = BTreeMap::new();
    for (&(lag, order), coefficient) in &shape.terms {
        for power in 0..x_width(coefficient) {
            let c = x_column(coefficient, power);
            if c.is_empty() {
                continue;
            }
            // [x^kappa] (x^power D^order P_(n - lag)) is the falling factorial
            // of kappa - power + order of length order times
            // [x^(kappa - power + order)] P_(n - lag), whose index in the
            // coordinates of row n - lag is j below.
            let offset = order as i64 - power as i64;
            let source = affine_add(&anchor, &affine_constant(offset));
            let sign = if coordinates == Coordinates::Identity {
                1
            } else {
                -1
            };
            let Some(j) = (&slope * as_q(lag))
                .to_integer()
                .to_i64()
                .map(|shift| shift + sign * offset)
                .filter(|j| *j >= 0)
            else {
                continue;
            };
            let mut product = IndexFunction::polynomial(c);
            for t in 0..order {
                product = times(
                    &product,
                    &affine_function(&affine_add(&source, &affine_constant(-(t as i64)))),
                );
            }
            for t in 1..=j as usize {
                product = times(&product, &ratio.at(t).shifted(-(lag as i64)));
            }
            let entry = lag_coefficients
                .entry(lag)
                .or_insert_with(|| IndexFunction::constant(Q::zero()));
            *entry = entry.add(&product);
        }
    }
    lag_coefficients.retain(|_, f| !vec_poly::is_zero(&f.numerator));
    // Every lag of the rule: a skipped term (zero by the support) also needs
    // the support facts at its row.
    let max_lag = shape.terms.keys().map(|(lag, _)| *lag).max().unwrap_or(1);
    let start = rule_start.max(s + max_lag);
    let values: Vec<(usize, Q)> = (start - max_lag..first_index + rows.len())
        .map(|n| {
            let row = &rows[n - first_index];
            let index = &anchor.0 + &anchor.1 * as_q(n);
            let value = index
                .to_integer()
                .to_usize()
                .and_then(|k| row.get(k).cloned())
                .unwrap_or_else(Q::zero);
            (n, value)
        })
        .collect();
    let mut first_coefficient = first_coefficient_hypothesis(
        &lag_coefficients,
        &IndexFunction::polynomial(shape.denominator.clone()),
        &values,
        start,
        max_lag,
    );
    // Otherwise the shared leading-coefficient module, applied to U_0 as a
    // family of constants: its scalar recurrence with seeds U_0(n),
    // start - max_lag <= n < start.
    if first_coefficient.status != LiftHypothesisStatus::Verified {
        if let Some(hypothesis) = shared_first_coefficient_hypothesis(
            &lag_coefficients,
            &shape.denominator,
            &values,
            start,
            max_lag,
        ) {
            first_coefficient = hypothesis;
        }
    }
    report.hypotheses.push(first_coefficient);

    // The region and the rows before the induction start.
    let full_support = coordinates == Coordinates::Identity || {
        let size = complement_size(&ratio.numerators[matching.base_numerator]);
        size.as_ref() == upper.as_ref()
    };
    let closed = valuation != affine_constant(0) || !full_support;
    let region = RealInterval {
        lower: None,
        lower_closed: false,
        upper: Some(Q::zero()),
        upper_closed: closed,
    };
    let in_region = |row: &[Q]| {
        let row = vec_poly::trim(row);
        !vec_poly::is_zero(&row)
            && is_real_rooted_rational(&row)
            && if closed {
                row.iter().all(|c| !c.is_negative()) || row.iter().all(|c| !c.is_positive())
            } else {
                has_real_roots_below(&row, &Q::zero())
            }
    };
    let failing = (first_index..start).find(|n| !in_region(&rows[n - first_index]));
    report.hypotheses.push(exact_hypothesis(
        "earlier_rows",
        &format!(
            "the rows P_{first_index}, ..., P_{} before the induction start are nonzero and \
             real-rooted with every root in {}",
            start - 1,
            RootInterval::from_real(&region)
        ),
        match failing {
            None => Ok("checked exactly".into()),
            Some(n) => Err((Some(n), format!("P_{n} fails the exact check"))),
        },
    ));
    if let Some(requested) = &options.root_region {
        let contains = requested.lower.is_none()
            && requested
                .upper
                .as_ref()
                .is_none_or(|upper| upper.is_positive() || (upper.is_zero() && !closed));
        report.hypotheses.push(exact_hypothesis(
            "requested_region",
            "the requested root region contains the certified region",
            if contains {
                Ok(format!("{} is contained", RootInterval::from_real(&region)))
            } else {
                Err((
                    None,
                    format!(
                        "{} is not contained in the requested region",
                        RootInterval::from_real(&region)
                    ),
                ))
            },
        ));
    }
    let size = complement_size(&ratio.numerators[matching.base_numerator])
        .map(|size| affine_text(&size))
        .unwrap_or_default();
    let statement = format!(
        "every P_n, n >= {first_index}, is nonzero and real-rooted with every root in {}: for \
         n >= {s} the coefficients U_j ({}) satisfy U_j = U_0(n) prod_(t <= j) R(n, t), \
         R(n, j) = {}, so sum_j U_j y^j is U_0(n) times the image of (1 + y)^({size}) under a \
         nonnegative multiplier sequence, and U_0(n) != 0 for n >= {} by its scalar recurrence",
        RootInterval::from_real(&region),
        match coordinates {
            Coordinates::Identity => format!(
                "U_j = [x^(j + v(n))] P_n, v(n) = {}",
                affine_text(&valuation)
            ),
            Coordinates::Reversal => format!("U_j = [x^(d(n) - j)] P_n, d(n) = {}", {
                affine_text(&anchor)
            }),
        },
        ratio.text(),
        start - max_lag
    );
    Ok(report.finish(RealRootednessTheoremConclusion {
        first_index,
        zero_residue_classes: None,
        interlacing_from: None,
        interlacing: None,
        simple_roots: false,
        root_region: Some(RootInterval::from_real(&region)),
        statement,
    }))
}

mod matching;
use matching::*;
mod hypotheses;
use hypotheses::*;
