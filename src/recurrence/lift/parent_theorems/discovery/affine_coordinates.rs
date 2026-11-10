//! The affine-coordinate stage of proof discovery (module documentation of
//! `discovery`, step 5).
//!
//! # The transformation
//!
//! For rational `alpha != 0` and `beta` put `x = alpha u + beta` and
//! `Q_n(u) = P_n(alpha u + beta)`.  If the family satisfies
//! `f(n, x) P_n = sum_j sign_j(n) c_j(n, x) D_x^(d_j) P_(n - r_j)`, then by
//! the chain rule (`D_u = alpha D_x`) the rows `Q_n` satisfy
//!
//! ```text
//! f(n, alpha u + beta) Q_n = sum_j sign_j(n) alpha^(-d_j) c_j(n, alpha u + beta) D_u^(d_j) Q_(n - r_j)
//! ```
//!
//! with the seeds `Q_k(u) = P_k(alpha u + beta)`, for every index: an exact
//! identity of the two records, not a finite check.  `P_n` is zero exactly
//! when `Q_n` is, and the roots of `P_n` are `alpha r + beta` for the roots
//! `r` of `Q_n` (with their multiplicities), so
//!
//! - real-rootedness, nonzero rows, zero residue classes and simple roots
//!   transfer row by row;
//! - a root region `R` of the rows `Q_n` gives the region `alpha R + beta`
//!   of the rows `P_n` (the ends swap when `alpha < 0`; an index-dependent
//!   end `v + p n` becomes `alpha v + beta + alpha p n`, an exact quadratic
//!   end is mapped exactly and its rational outer bound stays outer);
//! - interlacing: `alpha > 0` keeps the order of the roots, so `Q_(n-1)`
//!   interlacing `Q_n` (strictly or weakly) gives the same relation of
//!   `P_(n-1)` and `P_n`.  `alpha < 0` reverses the order of the roots: the
//!   roots of `P_(n-1)` and `P_n` alternate in the same way with the roles of
//!   "largest" and "smallest" exchanged ([`AffineInterlacingOrientation`]).
//!   For `deg P_n = deg P_(n-1) + 1` that is the same relation; for equal
//!   degrees it is `P_n` interlacing `P_(n-1)`.  The stage records the
//!   orientation and does not claim more.
//!
//! # Candidate coordinates
//!
//! [`affine_coordinate_candidates`] derives a bounded, deterministic menu
//! from the rule alone: the rational roots `beta` of the factors free of
//! `n` of the coefficient polynomials (every term coefficient and the
//! left-hand factor, in the order of the rule; a coefficient `x - 1` gives
//! `u = x - 1`), each with `alpha = 1` and then `alpha = -1`, followed by
//! the reflection `x = -u`.  The identity is never a candidate, and at most
//! [`AffineCoordinateOptions::max_coordinates`] are tried.
//!
//! # The search and its replay
//!
//! Each candidate transforms the record (and the tail record, when discovery
//! has one) and runs every theorem route of the registry on it, in the
//! registry order, with the requested root region pulled back to the new
//! coordinate.  Depth one: no other stage runs on the transformed record.
//! The search stops at the first coordinate with a proof; the conclusion is
//! mapped back as above ([`AffineCoordinateConclusion`]).  Replay
//! ([`AffineCoordinateSearch::verify`]) recomputes the candidate menu, every
//! transformed record and the reference attempts, replays the recorded
//! attempts against them and recomputes the conclusion.

use super::*;
use crate::factorization::irreducible_factorization_rational_coeffs;
use crate::recurrence::{ExactQuadraticEndpoint, RootEndpoint};

/// Bounds of the affine-coordinate stage of
/// [`discover_recurrence_proofs_with_stage_options`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AffineCoordinateOptions {
    /// Largest number of candidate coordinates tried (deterministic order;
    /// `0` disables the stage).  Default 8.
    pub max_coordinates: usize,
}

impl Default for AffineCoordinateOptions {
    fn default() -> Self {
        Self { max_coordinates: 8 }
    }
}

impl AffineCoordinateOptions {
    /// The options that disable the stage.
    pub fn disabled() -> Self {
        Self { max_coordinates: 0 }
    }

    /// Whether the stage runs (`max_coordinates > 0`).
    pub fn is_enabled(&self) -> bool {
        self.max_coordinates > 0
    }
}

/// One affine coordinate `x = slope u + intercept`, with the rows
/// `Q_n(u) = P_n(slope u + intercept)`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AffineCoordinate {
    /// `alpha != 0`, an exact rational string.
    pub slope: String,
    /// `beta`, an exact rational string.
    pub intercept: String,
    /// Why the coordinate is a candidate (the coefficient root it comes
    /// from, or the reflection).
    pub reason: String,
}

impl AffineCoordinate {
    /// The exact slope and intercept.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::InvalidInput`] for a malformed rational or a zero
    /// slope.
    pub fn parameters(&self) -> Result<(Q, Q), RecurrenceError> {
        let slope = parse_rational_coeff(&self.slope)?;
        let intercept = parse_rational_coeff(&self.intercept)?;
        if slope.is_zero() {
            return Err(RecurrenceError::InvalidInput(
                "an affine coordinate needs a nonzero slope".into(),
            ));
        }
        Ok((slope, intercept))
    }

    /// `x = alpha u + beta` as text.
    fn text(&self) -> String {
        format!("x = ({}) u + ({})", self.slope, self.intercept)
    }
}

/// How consecutive rows of the family interlace, given that the transformed
/// rows `Q_(n-1)` interlace `Q_n` (module documentation, "The
/// transformation").
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AffineInterlacingOrientation {
    /// `alpha > 0`: `P_(n-1)` interlaces `P_n` in the recorded sense.
    Preserved,
    /// `alpha < 0`: the recorded relation with the order of the roots
    /// reversed: listed in increasing order, the roots of `P_n` and
    /// `P_(n-1)` alternate starting with the smallest root of `P_n`
    /// (`r_1 <= s_1 <= r_2 <= ...`, strict for strict interlacing).  This is
    /// `P_(n-1)` interlacing `P_n` when `deg P_n = deg P_(n-1) + 1`, and
    /// `P_n` interlacing `P_(n-1)` when the degrees are equal.
    Mirrored,
}

/// The tail part of one candidate: the tail record of discovery in the new
/// coordinate.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AffineCoordinateTail {
    /// The transformed tail record (same first index as the tail).
    pub recurrence: RecurrenceJson,
    /// Every theorem attempt on it.
    pub attempts: Vec<RealRootednessTheoremReport>,
}

/// One candidate coordinate and the theorem attempts in it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AffineCoordinateAttempt {
    /// The coordinate.
    pub coordinate: AffineCoordinate,
    /// The transformed record: the rule with every coefficient substituted
    /// (`alpha^(-d)` on a derivative term of order `d`) and the transformed
    /// seeds, same first index.
    pub recurrence: RecurrenceJson,
    /// Every theorem attempt on the transformed record, in the registry
    /// order.
    pub attempts: Vec<RealRootednessTheoremReport>,
    /// The attempts on the transformed tail record, when discovery has a
    /// tail and no attempt on the transformed record proves it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tail: Option<AffineCoordinateTail>,
    /// Whether an attempt (on the record or its tail) is proved.
    pub proved: bool,
}

/// What a proved coordinate gives for the family (module documentation,
/// "The transformation").
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AffineCoordinateConclusion {
    /// The coordinate `x = slope u + intercept`.
    pub coordinate: AffineCoordinate,
    /// The theorem whose attempt in the new coordinate gives the conclusion.
    pub theorem: RealRootednessTheorem,
    /// Whether that attempt is on the transformed tail record.
    pub on_tail: bool,
    /// Every row `P_n`, `n >= from_index`, is nonzero and real-rooted
    /// (outside [`Self::zero_residue_classes`] when present).
    pub from_index: usize,
    /// The zero residue classes of the inner conclusion (the same rows are
    /// zero in both coordinates); omitted when absent.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub zero_residue_classes: Option<ZeroResidueClasses>,
    /// The inner root region mapped to the variable `x` (`slope R +
    /// intercept`); omitted when the inner attempt certifies none.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub root_region: Option<RootInterval>,
    /// Whether every row from `from_index` on has simple roots.
    pub simple_roots: bool,
    /// First index of the interlacing of consecutive rows (inner
    /// conclusion); omitted when not claimed.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub interlacing_from: Option<usize>,
    /// The interlacing of consecutive transformed rows (inner conclusion).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub interlacing: Option<RowInterlacing>,
    /// How that interlacing reads for the rows `P_n`; present exactly with
    /// [`Self::interlacing`].
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub interlacing_orientation: Option<AffineInterlacingOrientation>,
    /// The conclusion in words.
    pub statement: String,
}

/// The affine-coordinate stage (module documentation of `discovery`, step
/// 5).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AffineCoordinateSearch {
    /// The bounds used.
    pub options: AffineCoordinateOptions,
    /// Number of candidate coordinates the rule gives (before the bound).
    pub candidates_available: usize,
    /// The coordinates tried, in the menu order; the search stops at the
    /// first proved one.
    pub coordinates: Vec<AffineCoordinateAttempt>,
    /// The conclusion of the proved coordinate; omitted when none proves.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub conclusion: Option<AffineCoordinateConclusion>,
}

// ---------------------------------------------------------------------------
// Candidates and the transformation
// ---------------------------------------------------------------------------

/// The rational roots of the factor free of `n` of `p(n, x)` (rows indexed
/// by the power of `n`, each ascending in `x`): the gcd over `Q[x]` of the
/// rows, in the order of its irreducible factorization.  Empty for a
/// constant gcd or a zero polynomial.
pub(in crate::recurrence::lift::parent_theorems) fn rational_roots_of_index_free_factor(
    rows: &[Vec<Q>],
) -> Vec<Q> {
    let content = rows
        .iter()
        .map(|row| vec_poly::trim(row))
        .filter(|row| !vec_poly::is_zero(row))
        .fold(crate::Polynomial::new(Vec::new()), |g, row| {
            g.gcd(&crate::Polynomial::new(row))
        });
    if content.coeffs().len() <= 1 {
        return Vec::new();
    }
    let Ok(factorization) = irreducible_factorization_rational_coeffs(content.coeffs()) else {
        return Vec::new();
    };
    factorization
        .factors
        .iter()
        .filter(|(factor, _)| factor.len() == 2)
        .map(|(factor, _)| Q::new(-factor[0].clone(), factor[1].clone()))
        .collect()
}

/// The candidate menu of `rule` (module documentation, "Candidate
/// coordinates"), unbounded; the stage takes its first
/// [`AffineCoordinateOptions::max_coordinates`].
pub(super) fn affine_coordinate_candidates(rule: &Recurrence) -> Vec<AffineCoordinate> {
    let mut pieces: Vec<(String, &BivarPoly)> = rule
        .terms
        .iter()
        .map(|term| {
            (
                format!(
                    "the coefficient of D^{} P_(n-{})",
                    term.deriv_order, term.offset
                ),
                &term.coeff,
            )
        })
        .collect();
    if let Some(denominator) = &rule.denominator {
        pieces.push(("the left-hand factor".to_string(), denominator));
    }
    let mut roots: Vec<(Q, String)> = Vec::new();
    for (what, coefficient) in pieces {
        for root in rational_roots_of_index_free_factor(&coefficient.coeffs) {
            if root.is_zero() || roots.iter().any(|(known, _)| *known == root) {
                continue;
            }
            let reason = format!(
                "x - ({}) divides {what} for every n",
                format_rational_coeff(&root)
            );
            roots.push((root, reason));
        }
    }
    let mut out = Vec::new();
    for (root, reason) in roots {
        for slope in [Q::one(), -Q::one()] {
            out.push(AffineCoordinate {
                slope: format_rational_coeff(&slope),
                intercept: format_rational_coeff(&root),
                reason: reason.clone(),
            });
        }
    }
    out.push(AffineCoordinate {
        slope: "-1".into(),
        intercept: "0".into(),
        reason: "the reflection x = -u".into(),
    });
    out
}

/// `p(n, slope u + intercept)`, scaled by `scale`.
fn substitute_rows(rows: &[Vec<Q>], slope: &Q, intercept: &Q, scale: &Q) -> Vec<Vec<Q>> {
    rows.iter()
        .map(|row| {
            let row = vec_poly::trim(row);
            if vec_poly::is_zero(&row) {
                vec![Q::zero()]
            } else {
                vec_poly::scale(
                    &vec_poly::compose(&row, &[intercept.clone(), slope.clone()]),
                    scale,
                )
            }
        })
        .collect()
}

/// The rule in the coordinate `x = slope u + intercept` (module
/// documentation, "The transformation").
fn transformed_rule(rule: &Recurrence, slope: &Q, intercept: &Q) -> Recurrence {
    let substitute = |p: &BivarPoly, scale: &Q| {
        BivarPoly::from_coefficients_normalized(substitute_rows(&p.coeffs, slope, intercept, scale))
    };
    let inverse = slope.recip();
    Recurrence {
        terms: rule
            .terms
            .iter()
            .map(|term| {
                let mut scale = Q::one();
                for _ in 0..term.deriv_order {
                    scale *= &inverse;
                }
                RecurrenceTerm {
                    coeff: substitute(&term.coeff, &scale),
                    ..term.clone()
                }
            })
            .collect(),
        denominator: rule.denominator.as_ref().map(|f| substitute(f, &Q::one())),
        inhomogeneous: rule
            .inhomogeneous
            .as_ref()
            .map(|g| substitute(g, &Q::one())),
    }
}

/// The record in the coordinate (rule and seeds; same first index).
pub(super) fn transformed_record(
    record: &RecurrenceJson,
    coordinate: &AffineCoordinate,
) -> Result<RecurrenceJson, RecurrenceError> {
    let (slope, intercept) = coordinate.parameters()?;
    let (rule, first_index, seeds) = record.to_recurrence_parts()?;
    let seeds: Vec<Vec<Q>> = seeds
        .iter()
        .map(|seed| vec_poly::compose(&vec_poly::trim(seed), &[intercept.clone(), slope.clone()]))
        .collect();
    Ok(RecurrenceJson::from_recurrence_rational(
        &transformed_rule(&rule, &slope, &intercept),
        first_index,
        &seeds,
        None,
    ))
}

/// The image `{(r - intercept) / slope : r in region}` of an open requested
/// region (the region of the new variable `u`).
fn pulled_back_region(region: &RealInterval, slope: &Q, intercept: &Q) -> RealInterval {
    let map = |value: &Option<Q>| value.as_ref().map(|v| (v - intercept) / slope);
    let (lower, upper) = if slope.is_positive() {
        (map(&region.lower), map(&region.upper))
    } else {
        (map(&region.upper), map(&region.lower))
    };
    let (lower_closed, upper_closed) = if slope.is_positive() {
        (region.lower_closed, region.upper_closed)
    } else {
        (region.upper_closed, region.lower_closed)
    };
    RealInterval {
        lower,
        lower_closed,
        upper,
        upper_closed,
    }
}

/// The theorem options in the coordinate: the requested region pulled back.
fn transformed_options(
    options: &RealRootednessTheoremOptions,
    coordinate: &AffineCoordinate,
) -> Result<RealRootednessTheoremOptions, RecurrenceError> {
    let (slope, intercept) = coordinate.parameters()?;
    Ok(RealRootednessTheoremOptions {
        root_region: options
            .root_region
            .as_ref()
            .map(|region| pulled_back_region(region, &slope, &intercept)),
        maximum_start_delay: options.maximum_start_delay,
    })
}

/// The endpoint `slope e + intercept` (module documentation).
fn mapped_endpoint(
    end: &RootEndpoint,
    slope: &Q,
    intercept: &Q,
) -> Result<RootEndpoint, RecurrenceError> {
    let value = parse_rational_coeff(&end.value)? * slope + intercept;
    let per_index = match &end.per_index {
        Some(p) => {
            let p = parse_rational_coeff(p)? * slope;
            (!p.is_zero()).then(|| format_rational_coeff(&p))
        }
        None => None,
    };
    let exact = match &end.exact {
        Some(exact) => Some(ExactQuadraticEndpoint {
            rational_part: format_rational_coeff(
                &(parse_rational_coeff(&exact.rational_part)? * slope + intercept),
            ),
            sqrt_coefficient: format_rational_coeff(
                &(parse_rational_coeff(&exact.sqrt_coefficient)? * slope),
            ),
            radicand: exact.radicand.clone(),
        }),
        None => None,
    };
    Ok(RootEndpoint {
        value: format_rational_coeff(&value),
        per_index,
        closed: end.closed,
        exact,
    })
}

/// The region `slope R + intercept` of the rows `P_n` from the region `R`
/// of the rows `Q_n` (the ends swap for a negative slope).
pub(super) fn mapped_root_region(
    region: &RootInterval,
    slope: &Q,
    intercept: &Q,
) -> Result<RootInterval, RecurrenceError> {
    let map = |end: &Option<RootEndpoint>| {
        end.as_ref()
            .map(|e| mapped_endpoint(e, slope, intercept))
            .transpose()
    };
    Ok(if slope.is_positive() {
        RootInterval {
            lower: map(&region.lower)?,
            upper: map(&region.upper)?,
        }
    } else {
        RootInterval {
            lower: map(&region.upper)?,
            upper: map(&region.lower)?,
        }
    })
}

// ---------------------------------------------------------------------------
// The search
// ---------------------------------------------------------------------------

/// The proved attempt that [`proved_rows`] reads its rows from: the first
/// with the smallest first index without zero classes, else the first with
/// the smallest first index.
fn deciding_attempt(
    attempts: &[RealRootednessTheoremReport],
) -> Option<&RealRootednessTheoremReport> {
    let proved = || {
        attempts
            .iter()
            .filter(|attempt| attempt.is_proved() && attempt.conclusion.is_some())
    };
    let first_index = |attempt: &RealRootednessTheoremReport| {
        attempt
            .conclusion
            .as_ref()
            .map_or(usize::MAX, |c| c.first_index)
    };
    let every_row = proved()
        .filter(|attempt| {
            attempt
                .conclusion
                .as_ref()
                .is_some_and(|c| !c.has_zero_residue_classes())
        })
        .min_by_key(|attempt| first_index(attempt));
    every_row.or_else(|| proved().min_by_key(|attempt| first_index(attempt)))
}

/// The conclusion of a proved candidate (module documentation, "The
/// transformation").
fn conclude(
    candidate: &AffineCoordinateAttempt,
) -> Result<Option<AffineCoordinateConclusion>, RecurrenceError> {
    let (attempt, on_tail) = match deciding_attempt(&candidate.attempts) {
        Some(attempt) => (attempt, false),
        None => match candidate
            .tail
            .as_ref()
            .and_then(|tail| deciding_attempt(&tail.attempts))
        {
            Some(attempt) => (attempt, true),
            None => return Ok(None),
        },
    };
    let inner = attempt.conclusion.as_ref().expect("a proved conclusion");
    let coordinate = candidate.coordinate.clone();
    let (slope, intercept) = coordinate.parameters()?;
    let root_region = inner
        .root_region
        .as_ref()
        .map(|region| mapped_root_region(region, &slope, &intercept))
        .transpose()?;
    let orientation = inner.interlacing.map(|_| {
        if slope.is_positive() {
            AffineInterlacingOrientation::Preserved
        } else {
            AffineInterlacingOrientation::Mirrored
        }
    });
    let rows = match &inner.zero_residue_classes {
        None => format!(
            "every row P_n with n >= {} is nonzero and real-rooted",
            inner.first_index
        ),
        Some(zero) => format!(
            "every row P_n with n >= {} and n mod {} in {:?} is zero and every other row P_n \
             with n >= {} is nonzero and real-rooted",
            inner.first_index, zero.period, zero.residues, inner.first_index
        ),
    };
    let mut statement = format!(
        "In the coordinate {} (Q_n(u) = P_n(({}) u + ({})), the rule and the seeds \
         transformed exactly, D_u = ({}) D_x), {} proves{} for the transformed rows (named \
         P_n in its statement): {}. Mapped back to x: {rows}",
        coordinate.text(),
        coordinate.slope,
        coordinate.intercept,
        coordinate.slope,
        attempt.theorem.as_str(),
        if on_tail { " on the tail record" } else { "" },
        inner.statement,
    );
    if let Some(region) = &root_region {
        statement.push_str(&format!(
            ", with every root in {}",
            crate::recurrence::format_root_interval_exact(region)
        ));
    }
    statement.push('.');
    if let (Some(from), Some(kind), Some(orientation)) =
        (inner.interlacing_from, inner.interlacing, orientation)
    {
        let kind = match kind {
            RowInterlacing::Strict => "strictly",
            RowInterlacing::Weak => "weakly",
        };
        statement.push_str(&match orientation {
            AffineInterlacingOrientation::Preserved => {
                format!(" For n >= {from}, P_(n-1) interlaces P_n {kind} (the slope is positive).")
            }
            AffineInterlacingOrientation::Mirrored => format!(
                " For n >= {from}, the roots of P_(n-1) and P_n interlace {kind} in the mirrored \
                 order (the slope is negative): in increasing order the smallest root is one of \
                 P_n; this is P_(n-1) interlacing P_n when deg P_n = deg P_(n-1) + 1 and P_n \
                 interlacing P_(n-1) when the degrees are equal."
            ),
        });
    }
    Ok(Some(AffineCoordinateConclusion {
        coordinate,
        theorem: attempt.theorem,
        on_tail,
        from_index: inner.first_index,
        zero_residue_classes: inner.zero_residue_classes.clone(),
        root_region,
        simple_roots: inner.simple_roots,
        interlacing_from: inner.interlacing.and(inner.interlacing_from),
        interlacing: inner.interlacing,
        interlacing_orientation: orientation,
        statement,
    }))
}

/// Run the theorems on a transformed record (validated like the input).
fn attempts_on(
    record: &RecurrenceJson,
    options: &RealRootednessTheoremOptions,
) -> Result<Vec<RealRootednessTheoremReport>, RecurrenceError> {
    let (rule, first_index, seeds) = record.to_recurrence_parts()?;
    validate_discovery_input(&rule, first_index, &seeds, options)?;
    run_theorems(&rule, first_index, &seeds, options)
}

/// Try one coordinate on `record` and, when it gives no proof, on `tail`.
fn try_coordinate(
    coordinate: AffineCoordinate,
    record: &RecurrenceJson,
    tail: Option<&RecurrenceJson>,
    options: &RealRootednessTheoremOptions,
) -> Result<AffineCoordinateAttempt, RecurrenceError> {
    let inner_options = transformed_options(options, &coordinate)?;
    let recurrence = transformed_record(record, &coordinate)?;
    let attempts = attempts_on(&recurrence, &inner_options)?;
    let mut candidate = AffineCoordinateAttempt {
        coordinate,
        recurrence,
        attempts,
        tail: None,
        proved: false,
    };
    if !candidate.attempts.iter().any(|a| a.is_proved()) {
        if let Some(tail) = tail {
            let recurrence = transformed_record(tail, &candidate.coordinate)?;
            let attempts = attempts_on(&recurrence, &inner_options)?;
            candidate.tail = Some(AffineCoordinateTail {
                recurrence,
                attempts,
            });
        }
    }
    candidate.proved = candidate.attempts.iter().any(|a| a.is_proved())
        || candidate
            .tail
            .as_ref()
            .is_some_and(|tail| tail.attempts.iter().any(|a| a.is_proved()));
    Ok(candidate)
}

/// Run the stage on the plain record `record` (and its tail record `tail`
/// when discovery has one).
pub(super) fn affine_coordinate_search(
    record: &RecurrenceJson,
    tail: Option<&RecurrenceJson>,
    options: &RealRootednessTheoremOptions,
    affine_options: &AffineCoordinateOptions,
) -> Result<AffineCoordinateSearch, RecurrenceError> {
    let (rule, _, _) = record.to_recurrence_parts()?;
    let candidates = affine_coordinate_candidates(&rule);
    let mut search = AffineCoordinateSearch {
        options: affine_options.clone(),
        candidates_available: candidates.len(),
        coordinates: Vec::new(),
        conclusion: None,
    };
    for coordinate in candidates.into_iter().take(affine_options.max_coordinates) {
        let candidate = try_coordinate(coordinate, record, tail, options)?;
        let proved = candidate.proved;
        search.coordinates.push(candidate);
        if proved {
            break;
        }
    }
    search.conclusion = match search.coordinates.last() {
        Some(last) if last.proved => conclude(last)?,
        _ => None,
    };
    Ok(search)
}

impl AffineCoordinateSearch {
    /// Whether a coordinate proves the record.
    pub fn is_proved(&self) -> bool {
        self.conclusion.is_some()
    }

    /// Replay the stage for the plain record `record` and the tail record
    /// `tail` of discovery (module documentation, "The search and its
    /// replay").
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::VerificationFailed`] on any mismatch (another menu,
    /// a coordinate out of order, a transformed record that is not the
    /// transformation of the record, attempts that do not replay, a
    /// conclusion that does not follow), and the errors of the replayed
    /// attempts.
    pub fn verify(
        &self,
        record: &RecurrenceJson,
        tail: Option<&RecurrenceJson>,
        options: &RealRootednessTheoremOptions,
    ) -> Result<(), RecurrenceError> {
        let failed = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
        if !self.options.is_enabled() {
            return failed("an affine-coordinate stage is recorded with the stage disabled");
        }
        let (rule, _, _) = record.to_recurrence_parts()?;
        let menu = affine_coordinate_candidates(&rule);
        if self.candidates_available != menu.len() {
            return failed("the affine-coordinate menu has another size");
        }
        let expected_count = menu.len().min(self.options.max_coordinates);
        let proved_at = self.coordinates.iter().position(|c| c.proved);
        match proved_at {
            Some(position) if position + 1 != self.coordinates.len() => {
                return failed("the affine-coordinate search continues after a proof");
            }
            None if self.coordinates.len() != expected_count => {
                return failed("the affine-coordinate search stops before its bound");
            }
            _ if self.coordinates.len() > expected_count => {
                return failed("more affine coordinates than the recorded bound");
            }
            _ => {}
        }
        for (candidate, expected) in self.coordinates.iter().zip(&menu) {
            if candidate.coordinate != *expected {
                return failed("an affine coordinate is not the next one of the menu");
            }
            let inner_options = transformed_options(options, &candidate.coordinate)?;
            if candidate.recurrence != transformed_record(record, &candidate.coordinate)? {
                return failed("a transformed record is not the transformation of the record");
            }
            let reference = attempts_on(&candidate.recurrence, &inner_options)?;
            verify_attempts(
                &candidate.attempts,
                &reference,
                &candidate.recurrence,
                &inner_options,
            )?;
            let proved_directly = candidate.attempts.iter().any(|a| a.is_proved());
            match (&candidate.tail, tail) {
                (None, None) => {}
                (None, Some(_)) if proved_directly => {}
                (Some(_), _) if proved_directly => {
                    return failed("a transformed tail is recorded after a proof");
                }
                (Some(recorded), Some(tail)) => {
                    if recorded.recurrence != transformed_record(tail, &candidate.coordinate)? {
                        return failed("a transformed tail is not the transformation of the tail");
                    }
                    let reference = attempts_on(&recorded.recurrence, &inner_options)?;
                    verify_attempts(
                        &recorded.attempts,
                        &reference,
                        &recorded.recurrence,
                        &inner_options,
                    )?;
                }
                (Some(_), None) => {
                    return failed("a transformed tail is recorded without a tail");
                }
                (None, Some(_)) => {
                    return failed("the transformed tail is missing");
                }
            }
            let proved = proved_directly
                || candidate
                    .tail
                    .as_ref()
                    .is_some_and(|t| t.attempts.iter().any(|a| a.is_proved()));
            if candidate.proved != proved {
                return failed("an affine coordinate's proved flag does not match its attempts");
            }
        }
        let expected = match self.coordinates.last() {
            Some(last) if last.proved => conclude(last)?,
            _ => None,
        };
        if expected != self.conclusion {
            return failed("the affine-coordinate conclusion does not follow from its attempts");
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(n: i64, d: i64) -> Q {
        Q::new(n.into(), d.into())
    }

    /// `n P_n = (1 + 2x + n x^2) P_(n-1) + (x - 3) x P'_(n-1) + (-1)^n x P''_(n-2)`.
    fn rule_with_derivatives() -> RecurrenceJson {
        serde_json::from_value(serde_json::json!({
            "schema": "polytool.recurrence.v1",
            "first_index": 0,
            "initial_polynomials": [["1", "2"], ["0", "1", "-1"]],
            "recurrence": {
                "terms": [
                    {"offset": 1, "deriv_order": 0, "sign": "none",
                     "coeff": {"coeffs": [["1", "2", "0"], ["0", "0", "1"]]}},
                    {"offset": 1, "deriv_order": 1, "sign": "none",
                     "coeff": {"coeffs": [["0", "-3", "1"]]}},
                    {"offset": 2, "deriv_order": 2, "sign": "alternating_n",
                     "coeff": {"coeffs": [["0", "1"]]}}
                ],
                "denominator": {"coeffs": [["0"], ["1"]]}
            }
        }))
        .unwrap()
    }

    #[test]
    fn the_transformed_record_computes_the_transformed_rows() {
        let record = rule_with_derivatives();
        let (rule, first, seeds) = record.to_recurrence_parts().unwrap();
        let rows = rule.generate_rows_rational(&seeds, first, 9).unwrap();
        for (slope, intercept) in [
            (q(3, 2), q(-1, 3)),
            (q(-1, 1), q(2, 1)),
            (q(-2, 5), q(0, 1)),
        ] {
            let coordinate = AffineCoordinate {
                slope: format_rational_coeff(&slope),
                intercept: format_rational_coeff(&intercept),
                reason: String::new(),
            };
            let transformed = transformed_record(&record, &coordinate).unwrap();
            let (t_rule, t_first, t_seeds) = transformed.to_recurrence_parts().unwrap();
            let t_rows = t_rule.generate_rows_rational(&t_seeds, t_first, 9).unwrap();
            for (row, t_row) in rows.iter().zip(&t_rows) {
                let expected =
                    vec_poly::compose(&vec_poly::trim(row), &[intercept.clone(), slope.clone()]);
                assert_eq!(vec_poly::trim(&expected), vec_poly::trim(t_row));
            }
        }
    }

    #[test]
    fn the_menu_lists_coefficient_roots_then_the_reflection() {
        let (rule, _, _) = rule_with_derivatives().to_recurrence_parts().unwrap();
        let menu = affine_coordinate_candidates(&rule);
        let pairs: Vec<(String, String)> = menu
            .iter()
            .map(|c| (c.slope.clone(), c.intercept.clone()))
            .collect();
        // Only the derivative coefficient x(x - 3) has an n-free factor with
        // a nonzero rational root.
        assert_eq!(
            pairs,
            [("1", "3"), ("-1", "3"), ("-1", "0")].map(|(a, b)| (a.to_string(), b.to_string()))
        );
    }

    #[test]
    fn regions_map_with_swapped_ends_for_a_negative_slope() {
        let region = RootInterval {
            lower: Some(RootEndpoint {
                value: "-2".into(),
                per_index: Some("-1".into()),
                closed: true,
                exact: None,
            }),
            upper: Some(RootEndpoint {
                value: "1".into(),
                per_index: None,
                closed: false,
                exact: Some(ExactQuadraticEndpoint {
                    rational_part: "0".into(),
                    sqrt_coefficient: "1/2".into(),
                    radicand: "2".into(),
                }),
            }),
        };
        let mapped = mapped_root_region(&region, &q(-2, 1), &q(1, 1)).unwrap();
        let lower = mapped.lower.as_ref().unwrap();
        assert_eq!((lower.value.as_str(), lower.closed), ("-1", false));
        let exact = lower.exact.as_ref().unwrap();
        assert_eq!(
            (
                exact.rational_part.as_str(),
                exact.sqrt_coefficient.as_str()
            ),
            ("1", "-1")
        );
        let upper = mapped.upper.as_ref().unwrap();
        assert_eq!(
            (
                upper.value.as_str(),
                upper.per_index.as_deref(),
                upper.closed
            ),
            ("5", Some("2"), true)
        );
        mapped.check_exact_endpoints().unwrap();
        let pulled = pulled_back_region(
            &RealInterval {
                lower: Some(q(-1, 1)),
                lower_closed: false,
                upper: None,
                upper_closed: false,
            },
            &q(-2, 1),
            &q(1, 1),
        );
        assert_eq!(pulled.upper, Some(q(1, 1)));
        assert!(pulled.lower.is_none());
    }
}
