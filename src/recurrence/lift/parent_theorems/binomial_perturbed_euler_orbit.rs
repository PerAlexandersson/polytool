//! Binomial perturbations of an Euler-operator orbit (RealRooted issue #2513).
//!
//! We match the two operators in `Q(n)[u, theta]`, `theta = u D`, with
//! degree index `m = n + degree_offset`. The rational coefficient of the
//! lag-one identity is `2 + 1/(m-2)` and determines the offset; the
//! lag-one theta coefficient determines the constant `s`. Every remaining
//! coefficient is compared by cross multiplication, never by sampling.
//!
//! The second coordinate is the Eulerian Mobius map
//! `Q_m(u) = (1+u)^m P_n(u/(1+u))`. Its rule is computed by the exact
//! substitution transport. The source degree cap follows inductively from
//! `deg c_(lag,d) <= lag+d` and the seed caps. The transformed seeds must
//! agree with the canonical family, including seeds preceding the rule's
//! start. Nonzero scalar denominators give uniqueness at every later index.
//! Thus the transported identity holds for the whole family. The map
//! `t = u/(1+u)` increases from `(-1,0)` onto `(-inf,0)`; neither endpoint
//! is a root, so degrees, simplicity and directed adjacent interlacing
//! transfer without a pole or an added root.

use super::*;
use crate::recurrence::{DegreeSchedule, RecurrenceTransformation, TransportedRecurrence};
use num_traits::ToPrimitive;

const REFERENCE: &str = "Binomial perturbations of an Euler-operator orbit: s >= 1, kappa >= \
    s+4, Q_0=1, Q_1=1+2u, Q_2=1+kappa*u+kappa*u^2 and the exact lag-two theta \
    recurrence give simple roots in (-1,0) and strict adjacent interlacing; \
    q_m(x)=r*p_m(x)+(1-r)*x^m, L=x+(s/2)*(x^2-1)D, r=(kappa-4)/(2s). \
    RealRooted issue #2513, https://github.com/PerAlexandersson/RealRooted/issues/2513";

/// Exact operator, coordinate and seed certificate for
/// [`check_binomial_perturbed_euler_orbit_theorem`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BinomialPerturbedEulerOrbitCertificate {
    /// Constant rational `s >= 1`.
    pub s: String,
    /// Constant rational `kappa >= s+4`, determined by the seeds.
    pub kappa: String,
    /// The theorem's degree index is `m = n + degree_offset`.
    pub degree_offset: i64,
    /// `None` for Q coordinates; otherwise the exact Eulerian Mobius map.
    pub coordinate: Option<RecurrenceTransformation>,
    /// Exact transported rule and seeds in Q coordinates (same input index).
    pub transformed_recurrence: RecurrenceJson,
}

impl BinomialPerturbedEulerOrbitCertificate {
    /// Replay the full rational operator identity, the coordinate transport,
    /// seed identities, start range, nonvanishing denominators and exact
    /// parameter inequalities, independently of a theorem report.
    ///
    /// # Errors
    ///
    /// Invalid records or any failed identity or hypothesis.
    pub fn verify(
        &self,
        rule: &Recurrence,
        first_index: usize,
        seeds: &[Vec<Q>],
    ) -> Result<(), RecurrenceError> {
        let fail = RecurrenceError::VerificationFailed;
        let (rule, _) = rule_with_x_free_left_hand_factor(rule)?;
        let prepared = prepare(&rule, first_index, seeds).map_err(fail)?;
        if self.s != format_rational_coeff(&prepared.s)
            || self.degree_offset != prepared.shift
            || self.coordinate != prepared.coordinate
            || self.transformed_recurrence != prepared.record(first_index)
        {
            return Err(fail(
                "Euler-orbit operator or coordinate certificate differs".into(),
            ));
        }
        let kappa = seed_parameter(&prepared, first_index).map_err(fail)?;
        if self.kappa != format_rational_coeff(&kappa) {
            return Err(fail("Euler-orbit seed parameter differs".into()));
        }
        let hypotheses = hypotheses(&rule, first_index, seeds.len(), &prepared, &kappa);
        if hypotheses
            .iter()
            .any(|h| h.status != LiftHypothesisStatus::Verified)
        {
            return Err(fail(
                "Euler-orbit certificate has a failed hypothesis".into(),
            ));
        }
        Ok(())
    }
}

/// Check the exact binomial-perturbed Euler-orbit theorem of RealRooted
/// issue #2513, in Q coordinates or after the Eulerian Mobius map.
///
/// The degree index may be shifted by an integer. All supplied seeds must
/// be canonical rows, `s >= 1`, `kappa >= s+4`, and the rule must start at
/// degree at least three with a nonvanishing scalar denominator. The result
/// gives simple roots in `(-1,0)` (Q) or `(-inf,0)` (original coordinates)
/// and strict directed interlacing of consecutive supplied and later rows.
///
/// # Errors
///
/// Malformed recurrence records or seed normalization failures.
pub fn check_binomial_perturbed_euler_orbit_theorem(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    check_real_rootedness_theorem(
        RealRootednessTheorem::BinomialPerturbedEulerOrbit,
        rule,
        first_index,
        seeds,
        options,
    )
}

pub(super) const ROUTE: TheoremRoute = TheoremRoute {
    theorem: RealRootednessTheorem::BinomialPerturbedEulerOrbit,
    check_at_supplied_seeds: theorem_report,
    start_policy: StartPolicy::MovesStart,
    is_preferred: any_proof,
    certificate_name: Some("binomial-perturbed Euler-orbit"),
    discovery: DiscoveryPolicy::WhenApplicable,
    in_residue_class_menu: true,
    in_index_section_menu: true,
    parent: ParentPolicy {
        stage: ParentStage::WhileUnproved,
        region: ParentRegion::Unrestricted,
        record_only_if_applicable: true,
    },
};

fn q(value: i64) -> Q {
    Q::from_integer(value.into())
}

fn separable_coefficient(index: &[Q], variable: &[Q]) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(
        index.iter().map(|c| vec_poly::scale(variable, c)).collect(),
    )
}

/// The canonical rule, expanded from the theta form with `m=n+shift`.
fn canonical_rule(s: &Q, shift: i64) -> Recurrence {
    let m = vec![q(shift), q(1)];
    let denominator = vec_poly::sub(&m, &[q(2)]);
    let a = vec_poly::sub(&vec_poly::scale(&m, &q(2)), &[q(3)]);
    let b = vec_poly::neg(&vec_poly::sub(&m, &[q(1)]));
    Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, separable_coefficient(&a, &[q(1), q(2)])),
            RecurrenceTerm::new(
                1,
                1,
                separable_coefficient(&vec_poly::scale(&denominator, s), &[q(0), q(1), q(1)]),
            ),
            RecurrenceTerm::new(2, 0, separable_coefficient(&b, &[q(1), q(4), q(4)])),
            RecurrenceTerm::new(
                2,
                1,
                separable_coefficient(&vec_poly::scale(&b, s), &[q(0), q(1), q(3), q(2)]),
            ),
        ],
        denominator: Some(separable_coefficient(&denominator, &[q(1)])),
        inhomogeneous: None,
    }
}

fn constant_ratio(numerator: &[Q], denominator: &[Q]) -> Option<Q> {
    let pivot = denominator.iter().position(|c| !c.is_zero())?;
    let ratio = numerator.get(pivot).cloned().unwrap_or_else(Q::zero) / &denominator[pivot];
    (vec_poly::trim(numerator) == vec_poly::trim(&vec_poly::scale(denominator, &ratio)))
        .then_some(ratio)
}

/// Rewrite `u^j D` as `u^(j-1) theta` and compare every theta coefficient.
fn theta_coefficients(shape: &RecurrenceShape) -> Option<BTreeMap<(usize, usize, usize), Vec<Q>>> {
    let mut result = BTreeMap::new();
    for (&(lag, derivative), coefficient) in &shape.terms {
        if !(1..=2).contains(&lag) || derivative > 1 {
            return None;
        }
        for power in 0..x_width(coefficient) {
            let p = x_column(coefficient, power);
            if !p.is_empty() {
                if power < derivative {
                    return None;
                }
                result.insert((lag, power - derivative, derivative), p);
            }
        }
    }
    Some(result)
}

fn recognize(rule: &Recurrence) -> Result<(Q, i64), String> {
    let shape = recurrence_shape(rule)?;
    let bad = || "the operators are not the binomial-perturbed Euler-orbit theta form".to_string();
    let actual = theta_coefficients(&shape).ok_or_else(bad)?;
    let a = actual.get(&(1, 0, 0)).ok_or_else(bad)?;
    // A-2=1/(n+shift-2), so f - n*(a-2f) = (shift-2)*(a-2f).
    let difference = vec_poly::sub(a, &vec_poly::scale(&shape.denominator, &q(2)));
    let constant = constant_ratio(
        &vec_poly::sub(&shape.denominator, &vec_poly::shift(&difference)),
        &difference,
    )
    .ok_or_else(bad)?
        + q(2);
    if !constant.is_integer() {
        return Err(bad());
    }
    let shift = constant.to_integer().to_i64().ok_or_else(bad)?;
    let s = constant_ratio(
        actual.get(&(1, 0, 1)).map_or(&[], Vec::as_slice),
        &shape.denominator,
    )
    .ok_or_else(bad)?;
    let expected = recurrence_shape(&canonical_rule(&s, shift)).expect("canonical shape");
    let expected_theta = theta_coefficients(&expected).expect("canonical theta form");
    for key in actual.keys().chain(expected_theta.keys()) {
        let left = actual.get(key).map_or(&[][..], Vec::as_slice);
        let right = expected_theta.get(key).map_or(&[][..], Vec::as_slice);
        if vec_poly::trim(&vec_poly::mul(left, &expected.denominator))
            != vec_poly::trim(&vec_poly::mul(right, &shape.denominator))
        {
            return Err(bad());
        }
    }
    Ok((s, shift))
}

struct Prepared {
    s: Q,
    shift: i64,
    coordinate: Option<RecurrenceTransformation>,
    rule: Recurrence,
    seeds: Vec<Vec<Q>>,
}

impl Prepared {
    fn record(&self, first: usize) -> RecurrenceJson {
        RecurrenceJson::from_recurrence_rational(&self.rule, first, &self.seeds, None)
    }
}

fn prepare(rule: &Recurrence, first: usize, seeds: &[Vec<Q>]) -> Result<Prepared, String> {
    if let Ok((s, shift)) = recognize(rule) {
        return Ok(Prepared {
            s,
            shift,
            coordinate: None,
            rule: rule.clone(),
            seeds: seeds.to_vec(),
        });
    }
    let shape = recurrence_shape(rule)?;
    if shape
        .terms
        .keys()
        .any(|&(lag, d)| !(1..=2).contains(&lag) || d > 1)
        || !shape.terms.contains_key(&(2, 0))
        || shape
            .terms
            .iter()
            .any(|(&(lag, d), c)| x_degree(c).is_some_and(|degree| degree > lag + d))
    {
        return Err("the operators are not the Euler-orbit form in either coordinate".into());
    }
    let degree = seeds
        .first()
        .and_then(|p| vec_poly::degree(p))
        .ok_or("a nonzero first seed is needed for the Mobius degree schedule")?;
    let coordinate = RecurrenceTransformation::FPolynomial {
        schedule: DegreeSchedule::new(first, degree, 1),
    };
    let TransportedRecurrence::Scalar(transformed) =
        coordinate.transport(rule).map_err(|e| e.to_string())?
    else {
        return Err("the Mobius transport is not scalar".into());
    };
    let (s, shift) = recognize(&transformed)?;
    if i128::from(shift) != degree as i128 - first as i128 {
        return Err("the Mobius degree schedule does not match the operator degree index".into());
    }
    let mapped = seeds
        .iter()
        .enumerate()
        .map(|(i, p)| {
            coordinate
                .map_row(first + i, p)
                .map(|mut rows| rows.remove(0))
                .map_err(|e| e.to_string())
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Prepared {
        s,
        shift,
        coordinate: Some(coordinate),
        rule: transformed,
        seeds: mapped,
    })
}

/// Canonical rows in the theorem's degree index, generated by its defining
/// recurrence from Q_0, Q_1, Q_2 (no fitted recurrence or sampled identity).
fn canonical_rows(s: &Q, kappa: &Q, count: usize) -> Result<Vec<Vec<Q>>, RecurrenceError> {
    let seeds = vec![
        vec![q(1)],
        vec![q(1), q(2)],
        vec![q(1), kappa.clone(), kappa.clone()],
    ];
    canonical_rule(s, 0)
        .generate_rows_rational(&seeds, 0, count.max(3))
        .map_err(|e| RecurrenceError::VerificationFailed(e.to_string()))
}

fn first_degree(prepared: &Prepared, first: usize) -> Result<usize, String> {
    usize::try_from(first as i128 + i128::from(prepared.shift))
        .map_err(|_| "the first theorem degree index is negative or too large".into())
}

fn seed_parameter(prepared: &Prepared, first: usize) -> Result<Q, String> {
    let degree = first_degree(prepared, first)?;
    let count = degree
        .checked_add(prepared.seeds.len())
        .ok_or("seed degree overflow")?;
    // A valid canonical seed has degree m. Checking this before generation
    // also bounds the work by the size of the supplied rows.
    for (i, seed) in prepared.seeds.iter().enumerate() {
        if vec_poly::degree(seed) != Some(degree + i) {
            return Err(format!(
                "Q_{} has degree {:?}, expected {}",
                degree + i,
                vec_poly::degree(seed),
                degree + i
            ));
        }
    }
    let zero = canonical_rows(&prepared.s, &q(0), count).map_err(|e| e.to_string())?;
    let one = canonical_rows(&prepared.s, &q(1), count).map_err(|e| e.to_string())?;
    for (i, seed) in prepared.seeds.iter().enumerate() {
        let m = degree + i;
        let difference = vec_poly::sub(&one[m], &zero[m]);
        if !vec_poly::is_zero(&difference) {
            return constant_ratio(&vec_poly::sub(seed, &zero[m]), &difference)
                .ok_or_else(|| format!("Q_{m} is not a canonical seed for any kappa"));
        }
    }
    Err("the seeds do not determine Q_2's parameter kappa (a row of degree at least two is required)".into())
}

fn hypotheses(
    rule: &Recurrence,
    first: usize,
    length: usize,
    prepared: &Prepared,
    kappa: &Q,
) -> Vec<LiftHypothesis> {
    let shape = recurrence_shape(rule).expect("preparation checked the source shape");
    let mut result: Vec<_> = denominator_hypothesis(&shape.denominator, first + length)
        .into_iter()
        .collect();
    let exact = |name, statement, condition: bool, evidence: String| {
        exact_hypothesis(
            name,
            statement,
            if condition {
                Ok(evidence)
            } else {
                Err((None, evidence))
            },
        )
    };
    result.push(exact(
        "euler_parameter_at_least_one",
        "s >= 1",
        prepared.s >= q(1),
        format!("s = {}", format_rational_coeff(&prepared.s)),
    ));
    result.push(exact(
        "binomial_parameter_bound",
        "kappa >= s+4",
        *kappa >= prepared.s.clone() + q(4),
        format!(
            "kappa = {}, s+4 = {}",
            format_rational_coeff(kappa),
            format_rational_coeff(&(prepared.s.clone() + q(4)))
        ),
    ));
    let degree = first_degree(prepared, first);
    let starts = degree
        .as_ref()
        .is_ok_and(|d| d.checked_add(length).is_some_and(|m| m >= 3));
    result.push(exact(
        "theorem_index_range",
        "the seed indices are nonnegative and the rule starts at m >= 3 (m-2 != 0)",
        starts,
        format!(
            "m = n + ({}), source rule starts at n = {}",
            prepared.shift,
            first + length
        ),
    ));
    let checked = degree.and_then(|d| {
        let rows = canonical_rows(&prepared.s, kappa, d + length).map_err(|e| e.to_string())?;
        for (i, seed) in prepared.seeds.iter().enumerate() {
            if vec_poly::trim(seed) != vec_poly::trim(&rows[d + i]) {
                return Err(format!("seed Q_{} differs from the canonical row", d + i));
            }
        }
        Ok(format!(
            "all {length} supplied seeds agree with Q_0=1, Q_1=1+2u, \
             Q_2=1+kappa*u+kappa*u^2 and the defining recurrence"
        ))
    });
    result.push(exact_hypothesis(
        "canonical_seed_rows",
        "every supplied seed is exactly the canonical row at its theorem index",
        checked.map_err(|e| (None, e)),
    ));
    result
}

pub(super) fn verify_certificate_in_report(
    certificate: &BinomialPerturbedEulerOrbitCertificate,
    report: &RealRootednessTheoremReport,
) -> Result<(), RecurrenceError> {
    let (rule, first, seeds) = report.recurrence.to_recurrence_parts()?;
    let seeds = rule
        .generate_rows_rational(&seeds, first, seeds.len() + report.start_delay)
        .map_err(|e| RecurrenceError::VerificationFailed(e.to_string()))?;
    certificate.verify(&rule, first, &seeds)
}

fn theorem_report(
    rule: &Recurrence,
    first: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    let mut report = RealRootednessTheoremReport::new(
        RealRootednessTheorem::BinomialPerturbedEulerOrbit,
        REFERENCE,
        rule,
        first,
        seeds,
        options,
    );
    let prepared = match prepare(rule, first, seeds) {
        Ok(p) => p,
        Err(reason) => return Ok(report.not_applicable(reason)),
    };
    let coordinate_text = if prepared.coordinate.is_some() {
        "Q_m(u)=(1+u)^m P_n(u/(1+u)); exact Mobius transport and inductive degree cap"
    } else {
        "the supplied rows are in Q coordinates"
    };
    report.normal_form = format!(
        "m=n+({}); Q_m=[((2m-3)/(m-2))(2u+1)+s(u+1)theta]Q_(m-1)\
         -((m-1)/(m-2))(2u+1)[(2u+1)+s(u+1)theta]Q_(m-2), \
         theta=u D_u; s={}; {coordinate_text}",
        prepared.shift,
        format_rational_coeff(&prepared.s),
    );
    report.lean_declarations = if prepared.coordinate.is_some() {
        vec![
            "RealRooted.strictInterl_of_orbitPerturbation_recurrence".to_string(),
            "RealRooted.strictInterl_cayleyTransform_succ_of_roots_mem_Icc".to_string(),
        ]
    } else {
        vec!["RealRooted.strictInterl_of_orbitPerturbation_recurrence".to_string()]
    };
    let kappa = match seed_parameter(&prepared, first) {
        Ok(kappa) => kappa,
        Err(reason) => {
            report.hypotheses.push(exact_hypothesis(
                "canonical_seed_rows",
                "the seeds determine the canonical Q_2 parameter and agree exactly",
                Err((None, reason)),
            ));
            report.status = RealRootednessTheoremStatus::HypothesisFailed;
            return Ok(report);
        }
    };
    report
        .normal_form
        .push_str(&format!("; kappa={}", format_rational_coeff(&kappa)));
    report.hypotheses = hypotheses(rule, first, seeds.len(), &prepared, &kappa);
    let region = RealInterval {
        lower: prepared.coordinate.is_none().then(|| q(-1)),
        lower_closed: false,
        upper: Some(q(0)),
        upper_closed: false,
    };
    if let Some(requested) = &options.root_region {
        let contains = requested
            .lower
            .as_ref()
            .is_none_or(|lower| region.lower.as_ref().is_some_and(|bound| lower <= bound))
            && requested.upper.as_ref().is_none_or(|upper| upper >= &q(0));
        report.hypotheses.push(if contains {
            exact_hypothesis(
                "requested_root_region",
                "the requested region contains the certified open interval",
                Ok("exact endpoint comparison".into()),
            )
        } else {
            unresolved_hypothesis(
                "requested_root_region",
                "the requested region contains the certified open interval",
                "the requested region is narrower than the theorem's region".into(),
            )
        });
    }
    if report
        .hypotheses
        .iter()
        .all(|h| h.status == LiftHypothesisStatus::Verified)
    {
        report.certificate = Some(RouteCertificate::BinomialPerturbedEulerOrbit(
            BinomialPerturbedEulerOrbitCertificate {
                s: format_rational_coeff(&prepared.s),
                kappa: format_rational_coeff(&kappa),
                degree_offset: prepared.shift,
                coordinate: prepared.coordinate.clone(),
                transformed_recurrence: prepared.record(first),
            },
        ));
    }
    let region_text = if prepared.coordinate.is_some() {
        "(-inf,0)"
    } else {
        "(-1,0)"
    };
    let transfer = if prepared.coordinate.is_some() {
        "; t=u/(1+u) is increasing on (-1,0), so the directed interlacing order is preserved"
    } else {
        ""
    };
    Ok(report.finish(RealRootednessTheoremConclusion {
        first_index: first,
        zero_residue_classes: None,
        interlacing_from: Some(first + 1),
        interlacing: Some(RowInterlacing::Strict),
        simple_roots: true,
        root_region: Some(open_root_interval(&region)),
        statement: format!(
            "every P_n, n >= {first}, has n+({}) simple roots in {region_text}; \
             P_(n-1) strictly interlaces P_n for n > {first}{transfer}",
            prepared.shift,
        ),
    }))
}
