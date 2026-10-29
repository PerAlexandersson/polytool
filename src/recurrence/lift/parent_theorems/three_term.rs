//! Three-term recurrences: real-rootedness, strict interlacing and root
//! regions of orthogonal-polynomial type families.
//!
//! # The theorem
//!
//! Let `Q_n = (a_n y + b_n) Q_(n-1) - c_n Q_(n-2)` for `n >= N`, write `L_n`
//! for the leading coefficient of `Q_n` and `p_n = Q_n / L_n`.  If `Q_(N-2)`
//! and `Q_(N-1)` have simple real roots, `deg Q_(N-1) = deg Q_(N-2) + 1`,
//! `Q_(N-2)` strictly interlaces `Q_(N-1)` (vacuous for a constant
//! `Q_(N-2)`), `a_n != 0` for `n >= N`, and
//!
//! ```text
//! beta_N = c_N L_(N-2) / (a_N L_(N-1)) > 0,   beta_n = c_n / (a_n a_(n-1)) > 0  (n > N),
//! ```
//!
//! then every `Q_n`, `n >= N - 2`, has simple real roots and `Q_(n-1)`
//! strictly interlaces `Q_n`.  The monic rows satisfy
//! `p_n = (y - alpha_n) p_(n-1) - beta_n p_(n-2)` with `alpha_n = -b_n / a_n`;
//! at the roots `xi_1 < ... < xi_m` of `p_(n-1)` the values
//! `p_n(xi_i) = -beta_n p_(n-2)(xi_i)` alternate in sign and
//! `p_n(xi_m) < 0`, which places one root of `p_n` in each gap and one on each
//! side.  This is the zero-separation argument for orthogonal polynomials
//! (G. Szegő, *Orthogonal Polynomials*, AMS Colloquium Publications 23,
//! §3.3); with the seeds `Q_(N-2) = 1`, `Q_(N-1)` linear it is the setting of
//! Favard's theorem (J. Favard, C. R. Acad. Sci. Paris 200 (1935), 2052--2053),
//! where `p_n` are orthogonal for a positive measure.
//!
//! # Coordinates
//!
//! A rule `f(n) P_n = c_1(n, x) P_(n-1) + c_2(n, x) P_(n-2)` with `f` free of
//! `x`, `deg_x c_1 <= 1` and `c_2 = gamma(n) w(x)` is matched in two ways:
//!
//! - `w` constant: `Q = P`, `y = x`;
//! - `w = kappa (x + v_0)^2`: with `d_n = deg P_(N-2) + n - N + 2` and
//!   `Q_n(y) = y^(d_n) P_n(1/y - v_0)` (so `y = 1/(x + v_0)`), the rule
//!   becomes the three-term recurrence above, with `a_n`, `b_n` the
//!   coefficients of `c_1 / f` in the basis `1, x + v_0`.  Then
//!   `P_n(x) = L_n prod_i (1 - y_i (x + v_0))` over the roots `y_i` of `Q_n`,
//!   so `P_n` is real-rooted with simple roots whenever `Q_n` is (a root
//!   `y_i = 0` lowers the degree).  This is how the Narayana polynomials
//!   reduce to Jacobi polynomials `P^(1,1)` in `(1 + x)/(1 - x)`.
//!
//! # Root regions
//!
//! An open region `(l, u)` for the roots of `P` is pulled back to an open
//! interval `(s, t)` for the roots of `Q` (the identity, or the decreasing map
//! `x -> 1/(x + v_0)` when the region avoids the pole `-v_0`; an infinite end
//! maps to `0`, which then lies outside the open interval).  "Every root of
//! every `p_n` is `< t`" is certified by a rational function `g(n)` with
//!
//! 1. `g(n) > 0` for `n >= N`,
//! 2. `r_N >= g(N)` for `r_n = p_n(t) / p_(n-1)(t)` (exact rows),
//! 3. `(t - alpha_n) - beta_n / g(n-1) - g(n) >= 0` for `n > N`,
//!
//! after checking that both seeds have every root `< t` and `p_N(t) > 0`.
//! Indeed `r_n = (t - alpha_n) - beta_n / r_(n-1)` for `n > N`, so
//! `r_n >= g(n) > 0` by induction, and a monic `p_n` whose roots strictly
//! interlace those of `p_(n-1)`, all `< t`, has at most one root `>= t`, which
//! `p_n(t) > 0` excludes.  (This is the positive-ratio form of a
//! chain-sequence argument, T. S. Chihara, *An Introduction to Orthogonal
//! Polynomials*, 1978, Ch. III.)  The candidates for `g` are `(t - alpha_n) / 2`, which turns 3
//! into the classical `beta_n <= (t - alpha_(n-1)) (t - alpha_n) / 4`, and the
//! exact ratio `r_n` when the computed ratios are a rational function of `n`
//! of degree at most two over two.  Conditions 1 and 3 are decided by the
//! sign checker; a failed candidate is not a disproof, while a computed row
//! with `p_n(t) <= 0` (a root `>= t`) is.  The lower end `s` is the upper end
//! `-s` of the reflected rows `(-1)^(deg) p_n(-y)`.

use super::*;
use crate::linalg::solve_linear_system_with_rank;

const REFERENCE: &str = "three-term recurrence Q_n = (a_n y + b_n) Q_(n-1) - c_n Q_(n-2) with \
     a_n != 0 and beta_n = c_n / (a_n a_(n-1)) > 0: simple real roots and strict interlacing \
     of consecutive rows (zero separation for orthogonal polynomials: G. Szegő, Orthogonal \
     Polynomials, §3.3; Favard, C. R. Acad. Sci. Paris 200 (1935)); root regions by a \
     positive-ratio (chain-sequence) induction";

/// How `Q` relates to `P`.
#[derive(Debug, Clone)]
enum Coordinates {
    /// `Q = P`, `y = x`.
    Identity,
    /// `Q_n(y) = y^(d_n) P_n(1/y - v_0)`, `y = 1/(x + v_0)`.
    Reciprocal { v0: Q },
}

impl Coordinates {
    /// `Q` from `P` at the ambient degree `degree`; `None` when `deg P` is
    /// larger.
    fn transform(&self, p: &[Q], degree: usize) -> Option<Vec<Q>> {
        match self {
            Self::Identity => Some(vec_poly::trim(p)),
            Self::Reciprocal { v0 } => {
                vec_poly::reverse(&vec_poly::compose(p, &[-v0.clone(), Q::one()]), degree)
            }
        }
    }

    /// The open interval for the roots of `Q` that corresponds to the open
    /// `region` for the roots of `P`, as `(lower, upper)` ends (`None` is
    /// infinite).
    fn pull_back(&self, region: &RealInterval) -> Result<(Option<Q>, Option<Q>), String> {
        let (lower, upper) = (region.lower.clone(), region.upper.clone());
        if let (Some(l), Some(u)) = (&lower, &upper) {
            if l >= u {
                return Err("the requested root region is empty".into());
            }
        }
        let Self::Reciprocal { v0 } = self else {
            return Ok((lower, upper));
        };
        let pole = -v0.clone();
        let image = |x: &Q| (x + v0).recip();
        if lower.as_ref().is_some_and(|l| *l >= pole) {
            // x + v_0 > 0 on the region: y decreases from 1/(l + v_0) to 1/(u + v_0).
            let y_lower = upper.as_ref().map_or_else(Q::zero, image);
            let y_upper = lower.as_ref().filter(|l| **l != pole).map(image);
            Ok((Some(y_lower), y_upper))
        } else if upper.as_ref().is_some_and(|u| *u <= pole) {
            let y_upper = lower.as_ref().map_or_else(Q::zero, image);
            let y_lower = upper.as_ref().filter(|u| **u != pole).map(image);
            Ok((y_lower, Some(y_upper)))
        } else {
            Err(format!(
                "the region contains the pole x = {} of the coordinate change y = 1/(x{})",
                format_rational_coeff(&pole),
                signed_term(v0)
            ))
        }
    }
}

/// The monic recurrence `p_n = (y - alpha_n) p_(n-1) - beta_n p_(n-2)`.
struct MonicRecurrence {
    /// `alpha_n`, valid for `n >= first`.
    alpha: IndexFunction,
    /// `beta_n`, valid for `n > first` (`beta_first` enters only through the
    /// exact row `p_N`).
    beta: IndexFunction,
    /// `N`, the first index where the rule applies.
    first: usize,
    /// `p_(N-2), p_(N-1), ...`.
    rows: Vec<Vec<Q>>,
}

impl MonicRecurrence {
    /// The same recurrence for `(-1)^(deg) p_n(-y)`.
    fn reflected(&self) -> Self {
        let reflect = |p: &Vec<Q>| -> Vec<Q> {
            let degree = p.len().saturating_sub(1);
            p.iter()
                .enumerate()
                .map(|(k, c)| if (degree - k) % 2 == 1 { -c } else { c.clone() })
                .collect()
        };
        Self {
            alpha: self.alpha.negated(),
            beta: self.beta.clone(),
            first: self.first,
            rows: self.rows.iter().map(reflect).collect(),
        }
    }
}

/// Check the three-term recurrence theorem (module documentation of the
/// `three_term` submodule, summarized here) for the family given by `rule`,
/// `first_index` and the seed rows `seeds` (the rule applies from
/// `N = first_index + seeds.len()` on).
///
/// The rule must be `f(n) P_n = c_1(n, x) P_(n-1) + c_2(n, x) P_(n-2)` with
/// `f` free of `x`, `deg_x c_1 <= 1` and `c_2(n, x) = gamma(n) w(x)` where `w`
/// is constant or `kappa (x + v_0)^2`; in the second case the theorem is
/// applied to `Q_n(y) = y^(d_n) P_n(1/y - v_0)`.  Hypotheses, each decided
/// exactly:
///
/// - `f(n) != 0` and `a_n != 0` for `n >= N` (sign checker);
/// - the seeds `Q_(N-2)`, `Q_(N-1)` have simple real roots, degrees `e` and
///   `e + 1`, and strictly interlace; in the reciprocal case also
///   `deg Q_(N-2) = deg P_(N-2)`;
/// - `beta_N > 0` (exact) and `beta_n = c_n / (a_n a_(n-1)) > 0` for `n > N`
///   (sign checker);
/// - earlier seed rows real-rooted (exact);
/// - with `options.root_region`, the root bounds of the module
///   documentation (pulled back to `y`).
///
/// When all hold, every `P_n`, `n >= first_index`, is real-rooted (simple
/// roots from `N - 2` on), consecutive rows strictly interlace from `N - 1`
/// on (in the reciprocal case only when a root region is certified, since
/// `y -> x` is then monotone on it), and every root lies in the region.
///
/// # Errors
///
/// [`RecurrenceError`] when the rule cannot generate the rows needed for a
/// root-region certificate in a way other than a vanishing denominator.
///
/// # Example
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::recurrence::{
///     check_three_term_recurrence_theorem, BigRational, BivarPoly, Recurrence, RecurrenceTerm,
///     RealRootednessTheoremOptions,
/// };
/// use polytool::real_rootedness::RealInterval;
///
/// let q = |v: i64| BigRational::from_integer(BigInt::from(v));
/// // Chebyshev U: U_n = 2 x U_(n-1) - U_(n-2), U_0 = 1, U_1 = 2 x.
/// let rule = Recurrence {
///     terms: vec![
///         RecurrenceTerm::new(1, 0, BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(2)]])),
///         RecurrenceTerm::new(2, 0, BivarPoly::from_coefficients_normalized(vec![vec![q(-1)]])),
///     ],
///     denominator: None,
///     inhomogeneous: None,
/// };
/// let seeds = vec![vec![q(1)], vec![q(0), q(2)]];
/// let options = RealRootednessTheoremOptions {
///     root_region: Some(RealInterval::closed(q(-1), q(1))),
/// };
/// let report = check_three_term_recurrence_theorem(&rule, 0, &seeds, &options).unwrap();
/// assert!(report.is_proved());
/// report.verify().unwrap();
/// ```
pub fn check_three_term_recurrence_theorem(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<BigRational>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    let seeds = &super::normalized_seed_rows(rule, first_index, seeds)?;
    let mut report = RealRootednessTheoremReport::new(
        RealRootednessTheorem::ThreeTermRecurrence,
        REFERENCE,
        rule,
        first_index,
        seeds,
        options,
    );
    if seeds.len() < 2 {
        return Ok(report.not_applicable("the theorem needs at least two seed rows"));
    }
    let shape = match recurrence_shape(rule) {
        Ok(shape) => shape,
        Err(reason) => return Ok(report.not_applicable(reason)),
    };
    let keys: Vec<(usize, usize)> = shape.terms.keys().copied().collect();
    if keys != [(1, 0), (2, 0)] {
        return Ok(report.not_applicable(
            "the rule is not f(n) P_n = c_1(n, x) P_(n-1) + c_2(n, x) P_(n-2) with both terms \
             nonzero",
        ));
    }
    let (c1, c2) = (&shape.terms[&(1, 0)], &shape.terms[&(2, 0)]);
    if x_degree(c1).is_some_and(|d| d > 1) {
        return Ok(report.not_applicable("c_1(n, x) has degree > 1 in x"));
    }
    let Some((gamma, w)) = separate(c2) else {
        return Ok(report.not_applicable("c_2(n, x) is not of the form gamma(n) w(x)"));
    };
    let (coordinates, kappa) =
        match vec_poly::degree(&w) {
            Some(0) => (Coordinates::Identity, w[0].clone()),
            Some(2) if &w[1] * &w[1] == Q::from_integer(BigInt::from(4)) * &w[0] * &w[2] => (
                Coordinates::Reciprocal {
                    v0: &w[1] / (Q::from_integer(BigInt::from(2)) * &w[2]),
                },
                w[2].clone(),
            ),
            _ => return Ok(report.not_applicable(
                "c_2(n, x) = gamma(n) w(x) with w neither constant nor a constant times a square \
                 (x + v_0)^2",
            )),
        };
    let n_first = first_index + seeds.len();
    let start = n_first - 2;
    let column0 = x_column(c1, 0);
    let column1 = x_column(c1, 1);
    // Numerators (over f(n)) of a_n, b_n, c_n.
    let (a_numerator, b_numerator) = match &coordinates {
        Coordinates::Identity => (column1.clone(), column0.clone()),
        Coordinates::Reciprocal { v0 } => (
            vec_poly::sub(&column0, &vec_poly::scale(&column1, v0)),
            column1.clone(),
        ),
    };
    let c_numerator = vec_poly::scale(&gamma, &-kappa.clone());
    let f = IndexFunction::polynomial(shape.denominator.clone());
    let a = IndexFunction::new(a_numerator.clone(), shape.denominator.clone());
    let b = IndexFunction::new(b_numerator.clone(), shape.denominator.clone());
    let c = IndexFunction::new(c_numerator.clone(), shape.denominator.clone());
    let alpha = IndexFunction::new(vec_poly::neg(&b_numerator), a_numerator.clone());
    let beta = IndexFunction::new(
        vec_poly::mul(&c_numerator, &f.shifted(-1).numerator),
        vec_poly::mul(
            &a_numerator,
            &IndexFunction::polynomial(a_numerator.clone())
                .shifted(-1)
                .numerator,
        ),
    );
    report.normal_form = match &coordinates {
        Coordinates::Identity => format!(
            "P_n = (a_n x + b_n) P_(n-1) - c_n P_(n-2) for n >= {n_first}, with a_n = {}, b_n = \
             {}, c_n = {}",
            a.text(),
            b.text(),
            c.text()
        ),
        Coordinates::Reciprocal { v0 } => format!(
            "P_n = (a_n + b_n (x{plus})) P_(n-1) - c_n (x{plus})^2 P_(n-2) for n >= {n_first}, \
             with a_n = {}, b_n = {}, c_n = {}; Q_n(y) = y^(d_n) P_n(1/y{minus}), d_n = \
             deg P_{start} + n - {start}, satisfies Q_n = (a_n y + b_n) Q_(n-1) - c_n Q_(n-2)",
            a.text(),
            b.text(),
            c.text(),
            plus = signed_term(v0),
            minus = signed_term(&-v0.clone())
        ),
    };

    let mut hypotheses: Vec<LiftHypothesis> = Vec::new();
    hypotheses.extend(denominator_hypothesis(&shape.denominator, n_first));
    hypotheses.push(sign_hypothesis(
        "leading_coefficient_nonzero",
        &format!("a_n != 0 for every n >= {n_first}"),
        "f(n) a_n",
        &IndexFunction::polynomial(a_numerator.clone()),
        n_first,
        IndexSignCondition::Nonzero,
    ));
    let seed_rows = (&seeds[seeds.len() - 2], &seeds[seeds.len() - 1]);
    let ambient_degree = vec_poly::degree(seed_rows.0).unwrap_or(0);
    let transformed = (
        coordinates.transform(seed_rows.0, ambient_degree),
        coordinates.transform(seed_rows.1, ambient_degree + 1),
    );
    let seeds_outcome = check_seeds(
        &coordinates,
        seed_rows.0,
        &transformed,
        ambient_degree,
        start,
    );
    let seeds_ok = seeds_outcome.is_ok();
    hypotheses.push(exact_hypothesis(
        "seeds_interlace",
        &format!(
            "Q_{start} and Q_{} have simple real roots, degrees e and e + 1, and Q_{start} \
             strictly interlaces Q_{}",
            start + 1,
            start + 1
        ),
        seeds_outcome,
    ));
    let leading = |row: &Option<Vec<Q>>| {
        row.as_ref()
            .and_then(|r| r.last().cloned())
            .unwrap_or_else(Q::zero)
    };
    let (lead0, lead1) = (leading(&transformed.0), leading(&transformed.1));
    let beta_first = match (seeds_ok, a.at(n_first), c.at(n_first)) {
        (true, Some(a_value), Some(c_value)) if !a_value.is_zero() => {
            Some(c_value * &lead0 / (a_value * &lead1))
        }
        _ => None,
    };
    let statement = format!(
        "beta_{n_first} = c_{n_first} L_{start} / (a_{n_first} L_{}) > 0",
        start + 1
    );
    hypotheses.push(match &beta_first {
        Some(value) => exact_hypothesis(
            "first_coefficient_positive",
            &statement,
            if value.is_positive() {
                Ok(format!(
                    "beta_{n_first} = {} (exact)",
                    format_rational_coeff(value)
                ))
            } else {
                Err((
                    Some(n_first),
                    format!("beta_{n_first} = {} (exact)", format_rational_coeff(value)),
                ))
            },
        ),
        None if seeds_ok => exact_hypothesis(
            "first_coefficient_positive",
            &statement,
            Err((
                Some(n_first),
                format!("a_{n_first} or c_{n_first} is undefined or a_{n_first} = 0"),
            )),
        ),
        None => unresolved_hypothesis(
            "first_coefficient_positive",
            &statement,
            "not checked: the seeds do not satisfy the seed hypothesis".into(),
        ),
    });
    hypotheses.push(sign_hypothesis(
        "recurrence_coefficient_positive",
        &format!(
            "beta_n = c_n / (a_n a_(n-1)) > 0 for every n >= {}",
            n_first + 1
        ),
        "beta_n",
        &beta,
        n_first + 1,
        IndexSignCondition::Positive,
    ));
    let region = options.root_region.as_ref();
    hypotheses.extend(earlier_seeds_hypothesis(
        &seeds[..seeds.len() - 2],
        first_index,
        region,
    ));

    let core_verified = hypotheses
        .iter()
        .all(|h| h.status == LiftHypothesisStatus::Verified);
    if let Some(region) = region {
        let monic_recurrence = || -> Result<MonicRecurrence, String> {
            let generated = rule
                .generate_rows_rational(seeds, first_index, seeds.len() + GENERATED_ROWS)
                .map_err(|error| format!("the rows cannot be generated: {error}"))?;
            let rows = generated[start - first_index..]
                .iter()
                .enumerate()
                .map(|(offset, row)| {
                    coordinates
                        .transform(row, ambient_degree + offset)
                        .map(|q| monic(&q))
                        .ok_or_else(|| format!("row {} has too large a degree", start + offset))
                })
                .collect::<Result<Vec<_>, _>>()?;
            Ok(MonicRecurrence {
                alpha: alpha.clone(),
                beta: beta.clone(),
                first: n_first,
                rows,
            })
        };
        hypotheses.extend(root_region_hypotheses(
            &coordinates,
            region,
            core_verified,
            monic_recurrence,
            start,
        ));
    }

    let region_proved = region.is_some();
    let interlacing = match coordinates {
        Coordinates::Identity => true,
        Coordinates::Reciprocal { .. } => region_proved,
    };
    let root_region = region.map(open_root_interval);
    let statement = format!(
        "every P_n, n >= {first_index}, is nonzero and real-rooted, with simple roots from \
         n = {start} on{}{}",
        if interlacing {
            format!("; P_(n-1) strictly interlaces P_n for n >= {}", start + 1)
        } else {
            format!(
                "; Q_(n-1) strictly interlaces Q_n in y for n >= {}",
                start + 1
            )
        },
        root_region
            .as_ref()
            .map_or(String::new(), |r| format!("; every root lies in {r}"))
    );
    report.hypotheses = hypotheses;
    Ok(report.finish(RealRootednessTheoremConclusion {
        first_index,
        interlacing_from: interlacing.then_some(start + 1),
        interlacing: interlacing.then_some(RowInterlacing::Strict),
        simple_roots: true,
        root_region,
        statement,
    }))
}

/// `" + v"` or `" - |v|"`.
fn signed_term(value: &Q) -> String {
    if value.is_negative() {
        format!(" - {}", format_rational_coeff(&-value.clone()))
    } else {
        format!(" + {}", format_rational_coeff(value))
    }
}

/// The seed hypothesis on the transformed seeds.
fn check_seeds(
    coordinates: &Coordinates,
    first_seed: &[Q],
    transformed: &(Option<Vec<Q>>, Option<Vec<Q>>),
    ambient_degree: usize,
    start: usize,
) -> Result<String, (Option<usize>, String)> {
    let (Some(q0), Some(q1)) = transformed else {
        return Err((
            Some(start + 1),
            format!(
                "deg P_{} > deg P_{start} + 1, so Q_{} is not a polynomial",
                start + 1,
                start + 1
            ),
        ));
    };
    if vec_poly::is_zero(first_seed) {
        return Err((Some(start), format!("P_{start} is zero")));
    }
    if let Coordinates::Reciprocal { v0 } = coordinates {
        if vec_poly::degree(q0) != Some(ambient_degree) {
            return Err((
                Some(start),
                format!(
                    "P_{start} vanishes at x = {}, so deg Q_{start} < deg P_{start}",
                    format_rational_coeff(&-v0.clone())
                ),
            ));
        }
    }
    if !has_simple_real_roots(q0) {
        return Err((
            Some(start),
            format!("Q_{start} does not have simple real roots (exact)"),
        ));
    }
    let e = vec_poly::degree(q0).unwrap_or(0);
    if vec_poly::degree(q1) != Some(e + 1) {
        return Err((
            Some(start + 1),
            format!("deg Q_{} != deg Q_{start} + 1 = {}", start + 1, e + 1),
        ));
    }
    if !has_simple_real_roots(q1) {
        return Err((
            Some(start + 1),
            format!("Q_{} does not have simple real roots (exact)", start + 1),
        ));
    }
    if e > 0
        && crate::check_interlacing_bigint_coeffs(
            &integer_coefficients(q0),
            &integer_coefficients(q1),
        ) != Some(true)
    {
        return Err((
            Some(start + 1),
            format!(
                "Q_{start} does not strictly interlace Q_{} (exact Bezout test)",
                start + 1
            ),
        ));
    }
    Ok(format!(
        "Q_{start} has degree {e}, Q_{} degree {}; simple real roots and strict interlacing \
         checked exactly",
        start + 1,
        e + 1
    ))
}

/// The root-bound hypotheses for `region` (open at finite ends).
fn root_region_hypotheses(
    coordinates: &Coordinates,
    region: &RealInterval,
    core_verified: bool,
    monic_recurrence: impl FnOnce() -> Result<MonicRecurrence, String>,
    start: usize,
) -> Vec<LiftHypothesis> {
    let region_text = open_root_interval(region).to_string();
    let (y_lower, y_upper) = match coordinates.pull_back(region) {
        Ok(ends) => ends,
        Err(reason) => {
            return vec![unresolved_hypothesis(
                "root_region",
                &format!("every root of every P_n lies in {region_text}"),
                reason,
            )]
        }
    };
    let variable = match coordinates {
        Coordinates::Identity => "P",
        Coordinates::Reciprocal { .. } => "Q",
    };
    let mut bounds = Vec::new();
    if let Some(t) = y_upper {
        bounds.push(("root_upper_bound", "<", t, false));
    }
    if let Some(s) = y_lower {
        bounds.push(("root_lower_bound", ">", s, true));
    }
    if bounds.is_empty() {
        return Vec::new();
    }
    let recurrence = if core_verified {
        Some(monic_recurrence())
    } else {
        None
    };
    bounds
        .into_iter()
        .map(|(name, relation, bound, reflect)| {
            let statement = format!(
                "every root of {variable}_n is {relation} {} for every n >= {start} (region \
                 {region_text} for P)",
                format_rational_coeff(&bound)
            );
            match &recurrence {
                None => unresolved_hypothesis(
                    name,
                    &statement,
                    "not attempted: the real-rootedness hypotheses do not all hold".into(),
                ),
                Some(Err(reason)) => unresolved_hypothesis(name, &statement, reason.clone()),
                Some(Ok(recurrence)) => {
                    let mut hypothesis = if reflect {
                        upper_bound_hypothesis(name, &statement, &recurrence.reflected(), &-bound)
                    } else {
                        upper_bound_hypothesis(name, &statement, recurrence, &bound)
                    };
                    if reflect {
                        hypothesis.evidence = format!(
                            "on the reflected rows (-1)^(deg) p_n(-y): {}",
                            hypothesis.evidence
                        );
                    }
                    hypothesis
                }
            }
        })
        .collect()
}

/// "Every root of every `p_n`, `n >= N - 2`, is `< t`" (module
/// documentation).
fn upper_bound_hypothesis(
    name: &str,
    statement: &str,
    recurrence: &MonicRecurrence,
    t: &Q,
) -> LiftHypothesis {
    let start = recurrence.first - 2;
    let t_text = format_rational_coeff(t);
    let values: Vec<Q> = recurrence
        .rows
        .iter()
        .map(|p| vec_poly::evaluate(p, t))
        .collect();
    for (offset, seed) in recurrence.rows.iter().take(2).enumerate() {
        if !has_real_roots_below(seed, t) {
            return exact_hypothesis(
                name,
                statement,
                Err((
                    Some(start + offset),
                    format!("seed p_{} has a root >= {t_text} (exact)", start + offset),
                )),
            );
        }
    }
    if let Some(offset) = values.iter().position(|v| !v.is_positive()) {
        return exact_hypothesis(
            name,
            statement,
            Err((
                Some(start + offset),
                format!(
                    "p_{}({t_text}) = {} <= 0, so the monic p_{} has a root >= {t_text} (exact)",
                    start + offset,
                    format_rational_coeff(&values[offset]),
                    start + offset
                ),
            )),
        );
    }
    // r_n = p_n(t) / p_(n-1)(t) for n = N, N + 1, ...
    let ratios: Vec<(usize, Q)> = values
        .windows(2)
        .enumerate()
        .skip(1)
        .map(|(offset, pair)| (start + 1 + offset, &pair[1] / &pair[0]))
        .collect();
    let half = IndexFunction::constant(t.clone())
        .sub(&recurrence.alpha)
        .scale(&Q::new(BigInt::from(1), BigInt::from(2)));
    let mut candidates = Vec::new();
    if let Some(fitted) = fit_index_rational_function(&ratios) {
        candidates.push(("the exact ratio r_n", fitted));
    }
    candidates.push(("(t - alpha_n)/2", half));
    let mut failures = Vec::new();
    for (label, guess) in candidates {
        match check_ratio_bound(recurrence, t, &ratios[0].1, &guess) {
            Ok(decision) => {
                let mut hypothesis = LiftHypothesis::from_sign_decision(
                    name,
                    statement,
                    decision,
                    "(t - alpha_n) - beta_n/g(n-1) - g(n)",
                );
                hypothesis.evidence = format!(
                    "seeds checked exactly; positive-ratio induction at t = {t_text} with g(n) = \
                     {} ({label}): g(n) > 0 for n >= {n}, r_{n} >= g({n}) exactly, and {}",
                    guess.text(),
                    hypothesis.evidence,
                    n = recurrence.first,
                );
                return hypothesis;
            }
            Err(reason) => failures.push(format!("g(n) = {} ({label}): {reason}", guess.text())),
        }
    }
    unresolved_hypothesis(
        name,
        statement,
        format!(
            "no candidate closes the positive-ratio induction at t = {t_text} (not a disproof): {}",
            failures.join("; ")
        ),
    )
}

/// Conditions 1--3 of the module documentation for the candidate `g`, with
/// `first_ratio = r_N`; returns the decision of condition 3 when all hold.
fn check_ratio_bound(
    recurrence: &MonicRecurrence,
    t: &Q,
    first_ratio: &Q,
    g: &IndexFunction,
) -> Result<IndexSignDecision, String> {
    let n_first = recurrence.first;
    let positive = g.decide(n_first, IndexSignCondition::Positive);
    if !positive.holds() {
        return Err(positive.describe("g(n)"));
    }
    let g_first = g.at(n_first).ok_or("g is undefined at N")?;
    if first_ratio < &g_first {
        return Err(format!(
            "r_{n_first} = {} < g({n_first}) = {}",
            format_rational_coeff(first_ratio),
            format_rational_coeff(&g_first)
        ));
    }
    let condition = IndexFunction::constant(t.clone())
        .sub(&recurrence.alpha)
        .sub(&recurrence.beta.div(&g.shifted(-1)))
        .sub(g);
    let decision = condition.decide(n_first + 1, IndexSignCondition::Nonnegative);
    if decision.holds() {
        Ok(decision)
    } else {
        Err(decision.describe("(t - alpha_n) - beta_n/g(n-1) - g(n)"))
    }
}

/// The rational function `p(n) / q(n)` of least degrees (at most two over two,
/// `q` monic) through every point, when there is one; it must be determined
/// by the points with at least two to spare.
fn fit_index_rational_function(points: &[(usize, Q)]) -> Option<IndexFunction> {
    const DEGREES: [(usize, usize); 8] = [
        (0, 0),
        (1, 0),
        (0, 1),
        (1, 1),
        (2, 0),
        (2, 1),
        (1, 2),
        (2, 2),
    ];
    for (numerator_degree, denominator_degree) in DEGREES {
        let unknowns = numerator_degree + 1 + denominator_degree;
        if points.len() < unknowns + 2 {
            continue;
        }
        let mut matrix = Vec::with_capacity(points.len());
        let mut right = Vec::with_capacity(points.len());
        for (index, value) in points {
            let n = Q::from_integer(BigInt::from(*index));
            let powers: Vec<Q> = std::iter::successors(Some(Q::one()), |p| Some(p * &n))
                .take(numerator_degree.max(denominator_degree) + 1)
                .collect();
            let mut row: Vec<Q> = powers[..=numerator_degree].to_vec();
            row.extend(powers[..denominator_degree].iter().map(|p| -(value * p)));
            matrix.push(row);
            right.push(value * &powers[denominator_degree]);
        }
        let Some(solution) = solve_linear_system_with_rank(&matrix, &right) else {
            continue;
        };
        if solution.rank != unknowns {
            continue;
        }
        let numerator = solution.particular[..=numerator_degree].to_vec();
        let mut denominator = solution.particular[numerator_degree + 1..].to_vec();
        denominator.push(Q::one());
        return Some(IndexFunction::new(numerator, denominator));
    }
    None
}
