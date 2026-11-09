//! Residue-class subsequences: the RealRooted row route `subseq`.
//!
//! # The decomposition
//!
//! Let the rule be `f(n) P_n = sum c_(k,d)(n, x) D^d P_(n-k)` (homogeneous,
//! `f` free of `x`, no `(-1)^n` factor), applied from `N = n_0 + K` on (`n_0`
//! the first index, `K` seed rows), and let every lag `k` be a multiple of
//! `p >= 2`, the gcd of the lags.  For each residue `j = 0, ..., p - 1` the
//! class `Q^(j)_m = P_(n_0 + j + p m)`, `m >= 0`, satisfies the rule
//!
//! ```text
//! f(n_0 + j + p m) Q_m = sum c_(k,d)(n_0 + j + p m, x) D^d Q_(m - k/p)
//! ```
//!
//! exactly for every `m` with `n_0 + j + p m >= N`, that is from
//! `K_j = ceil((K - j) / p)` on, with the class rows of index `< N` as its
//! `K_j` seed rows (`K >= max k` gives `K_j >= max k / p`).  The classes
//! cover every row `P_n`, `n >= n_0` (RealRooted reassembles them with
//! `RealRooted.exists_mul_add_lt`), so the family is real-rooted when every
//! class is.  Each class family is proved by the first theorem of the
//! class menu (product, three-term, Liu--Wang, Liu--Wang derivative lag,
//! Ma--Wang, Liu--Wang three-term, iterated preserver, hypergeometric row,
//! fixed-delay window) that proves it, its own report replayed as a
//! sub-certificate.  A class whose seed rows are all zero is identically
//! zero (the rule is homogeneous and `f(n_0 + j + p m) != 0` for
//! `m >= K_j`, decided exactly), as in RealRooted's `eq_zero_of_*_rec`.
//! Such *zero classes* give infinitely many zero rows; the conclusion names
//! them in [`RealRootednessTheoremConclusion::zero_residue_classes`]
//! (residues of the index `n` modulo `p`), and its first index then bounds
//! only the rows outside them.
//!
//! The RealRooted tactic `rr_row_subseq_splits` handles periods `2..=6`
//! and class rules with one or two lags (product, derivative and three-term
//! shapes); here every period and every class rule is admitted, and the
//! class is proved by whichever polytool theorem applies.

use super::*;

const REFERENCE: &str = "Residue-class subsequences: when every lag of the rule is a multiple \
     of p >= 2, each class P_(n0 + j + p m), m >= 0, follows its own rule in m and the family \
     is real-rooted when every class is; RealRooted route subseq \
     (RealRooted/Tactic/Row/Subseq.lean, reassembly by RealRooted.exists_mul_add_lt)";

/// Largest induction-start delay a recorded certificate may ask the class
/// checks to try (as for the CLI's `--max-start-delay`).
pub(super) const MAXIMUM_CLASS_START_DELAY: usize = 16;

/// One residue class of a [`ResidueClassCertificate`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResidueClassProof {
    /// The residue `j`: the class is `P_(n_0 + j + p m)`, `m >= 0`.
    pub residue: usize,
    /// The class family in its own index `m` (first index `0`): the class
    /// rule and the class seed rows.
    pub recurrence: RecurrenceJson,
    /// Whether every seed row of the class is zero, so that every row of
    /// the class is zero.
    pub zero: bool,
    /// The class theorem checks that apply (status other than
    /// `not_applicable`), in menu order, up to the selected one; empty for a
    /// zero class.
    pub attempts: Vec<RealRootednessTheoremReport>,
    /// Position in `attempts` of the proof used for the class, when one
    /// proves it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub selected: Option<usize>,
}

/// The exact data of a [`check_residue_class_subsequence_theorem`] report.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResidueClassCertificate {
    /// The period `p >= 2`, the gcd of the lags.
    pub period: usize,
    /// The largest induction-start delay the class checks tried
    /// ([`RealRootednessTheoremOptions::maximum_start_delay`]); replayed by
    /// [`RealRootednessTheoremReport::verify`].
    pub class_maximum_start_delay: usize,
    /// The classes `j = 0, ..., p - 1`.
    pub classes: Vec<ResidueClassProof>,
}

/// Check the residue-class decomposition (RealRooted route `subseq`) for
/// the family given by `rule`, `first_index` and the seed rows `seeds` (the
/// rule applies from `N = first_index + seeds.len()` on).
///
/// It applies when every lag of the rule is a multiple of `p >= 2` (the gcd
/// of the lags), the rule is homogeneous, has no `(-1)^n` factor and a
/// left-hand factor `f(n)` free of `x`, and `seeds.len() >= max lag`.  Each
/// class `P_(first_index + j + p m)`, `m >= 0`, is an exactly derived family
/// in `m` (rule with `n = first_index + j + p m`, lags divided by `p`, the
/// class seed rows), proved by the first theorem of the class menu that
/// proves it (with `options`, including its start delays); a class with
/// only zero seed rows is identically zero when `f` does not vanish on its
/// rule range (decided exactly).  The class reports are kept in
/// [`RealRootednessTheoremReport::residue_classes`] and replayed by
/// [`RealRootednessTheoremReport::verify`].
///
/// When every class is proved, every row `P_n`, `n >= first_index`, is zero
/// (exactly the rows of the zero classes and the zero rows named in the
/// statement) or real-rooted; the root region is the hull of the class
/// regions when every live class certifies an index-independent one.  The
/// zero classes are recorded in
/// [`RealRootednessTheoremConclusion::zero_residue_classes`] (as residues
/// of `n` modulo `p`): with zero classes the conclusion does *not* say that
/// every row from its first index on is nonzero, only every row outside
/// them.  A rule whose classes are all zero (only zero seed rows) is not
/// applicable.  No
/// interlacing of consecutive rows is concluded (the classes interlace
/// within themselves only).  With a requested root region, every class is
/// checked with it and the hull must lie in it.  The induction start of the
/// whole family is not moved (each class moves its own).
///
/// # Errors
///
/// Errors of the record itself (seed normalization).
///
/// # Example
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::recurrence::{
///     check_residue_class_subsequence_theorem, BigRational, BivarPoly, Recurrence,
///     RecurrenceTerm, RealRootednessTheoremOptions,
/// };
///
/// let q = |v: i64| BigRational::from_integer(BigInt::from(v));
/// // P_n = (2 + 2x) P_(n-2) - P_(n-4), P_0 = P_1 = 1, P_2 = 2 + 2x, P_3 = 2 + 2x
/// // (OEIS A211956): both classes follow Q_m = (2 + 2x) Q_(m-1) - Q_(m-2).
/// let rule = Recurrence {
///     terms: vec![
///         RecurrenceTerm::new(2, 0, BivarPoly::from_coefficients_normalized(vec![
///             vec![q(2), q(2)],
///         ])),
///         RecurrenceTerm::new(4, 0, BivarPoly::from_coefficients_normalized(vec![vec![q(-1)]])),
///     ],
///     denominator: None,
///     inhomogeneous: None,
/// };
/// let seeds = vec![vec![q(1)], vec![q(1)], vec![q(2), q(2)], vec![q(2), q(2)]];
/// let report = check_residue_class_subsequence_theorem(
///     &rule,
///     0,
///     &seeds,
///     &RealRootednessTheoremOptions::default(),
/// )
/// .unwrap();
/// assert!(report.is_proved());
/// assert_eq!(report.residue_classes().unwrap().period, 2);
/// report.verify().unwrap();
/// ```
pub fn check_residue_class_subsequence_theorem(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<BigRational>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    check_real_rootedness_theorem(
        RealRootednessTheorem::ResidueClassSubsequences,
        rule,
        first_index,
        seeds,
        options,
    )
}

/// `p(first + period m)` for `p` ascending in `n`; the result is ascending
/// in `m`.
fn substitute_index(p: &[Q], first: usize, period: usize) -> Vec<Q> {
    let map = [
        Q::from_integer(BigInt::from(first)),
        Q::from_integer(BigInt::from(period)),
    ];
    if p.is_empty() {
        return Vec::new();
    }
    let mut composed = vec_poly::compose(p, &map);
    vec_poly::strip_trailing_zeros(&mut composed);
    composed
}

/// The bivariate coefficient `c(first + period m, x)` as `[i][j]` = the
/// coefficient of `m^i x^j`.
fn substitute_bivariate(b: &Bivariate, first: usize, period: usize) -> Vec<Vec<Q>> {
    let width = x_width(b);
    let columns: Vec<Vec<Q>> = (0..width)
        .map(|j| substitute_index(&x_column(b, j), first, period))
        .collect();
    let rows = columns.iter().map(Vec::len).max().unwrap_or(0).max(1);
    (0..rows)
        .map(|i| {
            (0..width.max(1))
                .map(|j| {
                    columns
                        .get(j)
                        .and_then(|column| column.get(i))
                        .cloned()
                        .unwrap_or_else(Q::zero)
                })
                .collect()
        })
        .collect()
}

/// The rule of the class `n = first + period m`.
fn class_rule(shape: &RecurrenceShape, first: usize, period: usize) -> Recurrence {
    let terms = shape
        .terms
        .iter()
        .map(|(&(lag, derivative), coefficient)| {
            RecurrenceTerm::new(
                lag / period,
                derivative,
                BivarPoly::from_coefficients_normalized(substitute_bivariate(
                    coefficient,
                    first,
                    period,
                )),
            )
        })
        .collect();
    let denominator = substitute_index(&shape.denominator, first, period);
    let denominator = (denominator != [Q::one()]).then(|| {
        BivarPoly::from_coefficients_normalized(denominator.into_iter().map(|c| vec![c]).collect())
    });
    Recurrence {
        terms,
        denominator,
        inhomogeneous: None,
    }
}

/// The smallest interval containing `a` and `b`.
fn hull(a: &RealInterval, b: &RealInterval) -> RealInterval {
    let (lower, lower_closed) = match (&a.lower, &b.lower) {
        (Some(x), Some(y)) if x < y => (Some(x.clone()), a.lower_closed),
        (Some(x), Some(y)) if y < x => (Some(y.clone()), b.lower_closed),
        (Some(x), Some(_)) => (Some(x.clone()), a.lower_closed || b.lower_closed),
        _ => (None, false),
    };
    let (upper, upper_closed) = match (&a.upper, &b.upper) {
        (Some(x), Some(y)) if x > y => (Some(x.clone()), a.upper_closed),
        (Some(x), Some(y)) if y > x => (Some(y.clone()), b.upper_closed),
        (Some(x), Some(_)) => (Some(x.clone()), a.upper_closed || b.upper_closed),
        _ => (None, false),
    };
    RealInterval {
        lower,
        lower_closed,
        upper,
        upper_closed,
    }
}

/// The RealRooted lemma that a class of zero seed rows stays zero, for the
/// class shapes of `rr_row_subseq_splits`.
fn zero_class_declaration(class: &Recurrence) -> Option<&'static str> {
    let lags: BTreeSet<usize> = class.terms.iter().map(|t| t.offset).collect();
    let derivative = class.terms.iter().map(|t| t.deriv_order).max().unwrap_or(0);
    match (lags.iter().max().copied(), derivative) {
        (Some(1), 0) => Some("RealRooted.eq_zero_of_product_rec"),
        (Some(1), 1) => Some("RealRooted.eq_zero_of_derivRec"),
        (Some(2), 0) => Some("RealRooted.eq_zero_of_threeTerm_rec"),
        _ => None,
    }
}

/// The registry entry of the residue-class theorem ([`super::registry::ROUTES`]).
pub(super) const ROUTE: TheoremRoute = TheoremRoute {
    theorem: RealRootednessTheorem::ResidueClassSubsequences,
    check_at_supplied_seeds: residue_class_report,
    start_policy: StartPolicy::OwnStarts,
    is_preferred: any_proof,
    certificate_name: Some("residue-class"),
    discovery: DiscoveryPolicy::WhenApplicable {
        tolerated_when_missing: false,
    },
    in_residue_class_menu: false,
    in_index_section_menu: true,
    parent: ParentPolicy {
        stage: ParentStage::WhileUnproved,
        region: ParentRegion::OpenNonpositiveRayThenUnrestricted,
        record_only_if_applicable: true,
    },
};

/// Replay the certificate of a recomputed residue-class report: the selected
/// proof of every residue class is verified on its own.
pub(super) fn verify_certificate_in_report(
    certificate: &ResidueClassCertificate,
    _report: &RealRootednessTheoremReport,
) -> Result<(), RecurrenceError> {
    for class in &certificate.classes {
        if let Some(position) = class.selected {
            let proof = class.attempts.get(position).ok_or_else(|| {
                RecurrenceError::VerificationFailed(
                    "a residue class selects a missing attempt".into(),
                )
            })?;
            if !proof.is_proved() || proof.recurrence != class.recurrence {
                return Err(RecurrenceError::VerificationFailed(
                    "a residue class selects an attempt that is not its proof".into(),
                ));
            }
            proof.verify()?;
        }
    }
    Ok(())
}

fn residue_class_report(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    let mut report = RealRootednessTheoremReport::new(
        RealRootednessTheorem::ResidueClassSubsequences,
        REFERENCE,
        rule,
        first_index,
        seeds,
        options,
    );
    let shape = match recurrence_shape(rule) {
        Ok(shape) => shape,
        Err(reason) => return Ok(report.not_applicable(reason)),
    };
    let period = shape.terms.keys().fold(0usize, |g, &(lag, _)| g.gcd(&lag));
    if period < 2 {
        return Ok(report.not_applicable(if shape.terms.is_empty() {
            "the rule has no nonzero term".to_string()
        } else {
            "the lags of the rule have gcd 1".to_string()
        }));
    }
    let maximum_lag = shape.terms.keys().map(|&(lag, _)| lag).max().unwrap_or(0);
    if seeds.len() < maximum_lag {
        return Ok(report.not_applicable(format!(
            "the rule has lag {maximum_lag} but only {} seed row(s)",
            seeds.len()
        )));
    }
    let n_first = first_index + seeds.len();
    report.normal_form = format!(
        "f(n) P_n = sum c_(k,d)(n, x) D^d P_(n-k) for n >= {n_first}, every lag k a multiple of \
         p = {period}; the classes P_({first_index} + j + {period} m), j = 0, ..., {}, follow \
         the rules f(n) Q_m = sum c_(k,d)(n, x) D^d Q_(m-k/{period}) with \
         n = {first_index} + j + {period} m",
        period - 1
    );
    report.lean_route = Some("subseq".into());
    let mut declarations = vec!["RealRooted.exists_mul_add_lt".to_string()];

    let class_options = RealRootednessTheoremOptions {
        root_region: options.root_region.clone(),
        maximum_start_delay: options.maximum_start_delay,
    };
    let mut hypotheses = Vec::new();
    let mut classes = Vec::new();
    let mut region: Option<Option<RealInterval>> = None;
    let mut zero_classes = Vec::new();
    let mut zero_rows = Vec::new();
    for residue in 0..period {
        let class_seeds: Vec<Vec<Q>> = seeds
            .iter()
            .skip(residue)
            .step_by(period)
            .cloned()
            .collect();
        let class_first = first_index + residue;
        let class = class_rule(&shape, class_first, period);
        let class_start = class_seeds.len();
        let record = RecurrenceJson::from_recurrence_rational(&class, 0, &class_seeds, None);
        let name = format!("residue_class_{residue}");
        let rows = format!("P_({class_first} + {period} m), m >= 0");
        if class_seeds
            .iter()
            .all(|row| vec_poly::is_zero(&vec_poly::trim(row)))
        {
            let f = substitute_index(&shape.denominator, class_first, period);
            let statement =
                format!("the class {rows} has only zero seed rows, so every row is zero");
            hypotheses.push(if f.len() <= 1 {
                exact_hypothesis(
                    &name,
                    &statement,
                    Ok("the rule is homogeneous with a nonzero constant left-hand factor".into()),
                )
            } else {
                sign_hypothesis(
                    &name,
                    &statement,
                    "f(class index m)",
                    &IndexFunction::polynomial(f),
                    class_start,
                    IndexSignCondition::Nonzero,
                )
            });
            if let Some(lemma) = zero_class_declaration(&class) {
                if !declarations.iter().any(|d| d == lemma) {
                    declarations.push(lemma.into());
                }
            }
            zero_classes.push(residue);
            classes.push(ResidueClassProof {
                residue,
                recurrence: record,
                zero: true,
                attempts: Vec::new(),
                selected: None,
            });
            continue;
        }
        let mut attempts = Vec::new();
        let mut selected = None;
        let mut fallback = None;
        // The theorems of the residue-class menu, in the route order.
        for route in ROUTES.iter().filter(|route| route.in_residue_class_menu) {
            let Ok(attempt) = check_with_first_proving_start(
                route.theorem,
                &class,
                0,
                &class_seeds,
                &class_options,
            ) else {
                continue;
            };
            if attempt.status == RealRootednessTheoremStatus::NotApplicable {
                continue;
            }
            let proved = attempt.is_proved();
            let with_region = proved_region(&attempt).is_some();
            attempts.push(attempt);
            if proved && (options.root_region.is_some() || with_region) {
                selected = Some(attempts.len() - 1);
                break;
            }
            if proved && fallback.is_none() {
                fallback = Some(attempts.len() - 1);
            }
        }
        let selected = selected.or(fallback);
        let statement = format!("the class {rows} is real-rooted by a theorem of the class menu");
        let outcome = match selected {
            Some(position) => {
                let proof = &attempts[position];
                let conclusion = proof
                    .conclusion
                    .as_ref()
                    .expect("a proved report concludes");
                let class_region = proved_region(proof);
                // Rows before the class conclusion's first index (zero seed
                // rows of the class theorem, possibly rows generated for a
                // moved induction start) are checked exactly.
                let early_rows = if conclusion.first_index <= class_seeds.len() {
                    Ok(class_seeds[..conclusion.first_index].to_vec())
                } else {
                    class
                        .generate_rows_rational(&class_seeds, 0, conclusion.first_index)
                        .map_err(|_| class_first + period * class_seeds.len())
                };
                let (early_rows, mut early) = match early_rows {
                    Ok(rows) => (rows, Ok(())),
                    Err(index) => (Vec::new(), Err(index)),
                };
                for (m, row) in early_rows.iter().enumerate() {
                    let row = vec_poly::trim(row);
                    let index = class_first + period * m;
                    if vec_poly::is_zero(&row) {
                        zero_rows.push(index);
                        continue;
                    }
                    let holds = match &class_region {
                        Some(r) => has_real_roots_in_closed_region(&row, r),
                        None => is_real_rooted_rational(&row),
                    };
                    if !holds {
                        early = Err(index);
                        break;
                    }
                }
                region = Some(match (region, class_region) {
                    (None, Some(r)) => Some(r),
                    (Some(Some(h)), Some(r)) => Some(hull(&h, &r)),
                    _ => None,
                });
                match early {
                    Ok(()) => Ok(format!(
                        "{:?} on the class family in its own index (class rule applied from \
                         m = {class_start}): {}",
                        proof.theorem, conclusion.statement
                    )),
                    Err(index) => Err((
                        Some(index),
                        format!(
                            "the class theorem leaves P_{index} out and that row is not \
                             real-rooted in the class region (exact)"
                        ),
                    )),
                }
            }
            None => Err((
                None,
                if attempts.is_empty() {
                    "no theorem of the class menu applies to the class rule".into()
                } else {
                    format!(
                        "no theorem of the class menu proves the class; open: {}",
                        attempts
                            .iter()
                            .map(|a| format!("{:?}: {}", a.theorem, a.open_hypotheses().join("; ")))
                            .collect::<Vec<_>>()
                            .join(" | ")
                    )
                },
            )),
        };
        let unresolved = selected.is_none()
            && attempts
                .iter()
                .any(|a| a.status == RealRootednessTheoremStatus::Unresolved);
        hypotheses.push(match outcome {
            Err((_, evidence)) if unresolved => unresolved_hypothesis(&name, &statement, evidence),
            outcome => exact_hypothesis(&name, &statement, outcome),
        });
        classes.push(ResidueClassProof {
            residue,
            recurrence: record,
            zero: false,
            attempts,
            selected,
        });
    }
    if zero_classes.len() == period {
        return Ok(report
            .not_applicable(
                "every residue class has only zero seed rows, so every row is zero".to_string(),
            )
            .without_hypotheses());
    }
    let region = region.flatten();
    if let Some(requested) = &options.root_region {
        let statement = format!(
            "the requested region {} contains the hull of the class regions",
            open_root_interval(requested)
        );
        hypotheses.push(match &region {
            Some(certified) if requested_region_contains(requested, certified) => {
                exact_hypothesis("requested_root_region", &statement, Ok("exact".into()))
            }
            Some(certified) => unresolved_hypothesis(
                "requested_root_region",
                &statement,
                format!(
                    "the classes certify the region {} only",
                    RootInterval::from_real(certified)
                ),
            ),
            None => unresolved_hypothesis(
                "requested_root_region",
                &statement,
                "some class certifies no index-independent root region".into(),
            ),
        });
    }
    report.hypotheses = hypotheses;
    report.lean_declarations = declarations;
    report.certificate = Some(RouteCertificate::ResidueClasses(ResidueClassCertificate {
        period,
        class_maximum_start_delay: options.maximum_start_delay,
        classes,
    }));
    zero_rows.sort_unstable();
    let conclusion_first = zero_rows.last().map_or(first_index, |last| last + 1);
    let mut statement = format!(
        "every P_n, n >= {first_index}, is zero or real-rooted, by the {period} residue classes \
         P_({first_index} + j + {period} m)"
    );
    // The zero classes as residues of the index itself (class j is
    // n = first_index + j mod p).
    let zero_residue_classes = (!zero_classes.is_empty()).then(|| {
        let mut residues: Vec<usize> = zero_classes
            .iter()
            .map(|j| (first_index % period + j) % period)
            .collect();
        residues.sort_unstable();
        ZeroResidueClasses { period, residues }
    });
    if let Some(zero) = &zero_residue_classes {
        statement.push_str(&format!(
            "; the rows of the residue classes j = {zero_classes:?} (the rows P_n with n mod \
             {period} in {:?}) are zero",
            zero.residues
        ));
    }
    if !zero_rows.is_empty() {
        statement.push_str(&format!("; the rows P_n, n in {zero_rows:?}, are zero"));
    }
    statement.push_str(&format!(
        "; every other row P_n, n >= {conclusion_first}, is nonzero and real-rooted"
    ));
    if let Some(r) = &region {
        statement.push_str(&format!(
            " with every root in {}",
            RootInterval::from_real(r)
        ));
    }
    Ok(report.finish(RealRootednessTheoremConclusion {
        first_index: conclusion_first,
        zero_residue_classes,
        interlacing_from: None,
        interlacing: None,
        simple_roots: false,
        root_region: region.as_ref().map(RootInterval::from_real),
        statement,
    }))
}
