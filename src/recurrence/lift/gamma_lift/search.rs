//! Detection of gamma lifts: the parent menu and [`find_gamma_lifts`].

use super::certify::{
    certify_gamma_lift, GammaLiftCandidate, GammaLiftCertification, GammaLiftOptions,
    GammaLiftParent, MAX_GAMMA_LIFT_CHECKED_ROWS, MAX_GAMMA_LIFT_INDEX_SHIFT,
};
use crate::recurrence::lift::hadamard_factor::prepare::prepare_rule;
use crate::recurrence::{
    format_rational_coeff, BivarPoly, DegreeSchedule, Recurrence, RecurrenceError, RecurrenceJson,
    RecurrenceTerm, RecurrenceTransformation, VeroneseSectionGauge,
};
use crate::vec_poly;
use crate::BigRational;
use num_bigint::BigInt;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};

type Q = BigRational;

/// JSON schema tag of a serialized [`GammaLiftReport`].
pub const GAMMA_LIFTS_SCHEMA: &str = "polytool.recurrence-gamma-lifts.v1";

/// `sum c n^i x^j` from the monomials `(i, j, c)`.
fn bivariate(monomials: &[(usize, usize, i64)]) -> BivarPoly {
    let n_size = monomials.iter().map(|m| m.0 + 1).max().unwrap_or(1);
    let x_size = monomials.iter().map(|m| m.1 + 1).max().unwrap_or(1);
    let mut coeffs = vec![vec![Q::zero(); x_size]; n_size];
    for &(i, j, c) in monomials {
        coeffs[i][j] += Q::from_integer(BigInt::from(c));
    }
    BivarPoly::from_coefficients_normalized(coeffs)
}

fn integer_rows(rows: &[&[i64]]) -> Vec<Vec<Q>> {
    rows.iter()
        .map(|row| {
            row.iter()
                .map(|&c| Q::from_integer(BigInt::from(c)))
                .collect()
        })
        .collect()
}

/// `P_n = (1 + (c (n - 1) + e) x) P_(n-1) + c x (1 - x) P'_(n-1)`.
fn derivative_rule(c: i64, e: i64) -> Recurrence {
    Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, bivariate(&[(0, 0, 1), (1, 1, c), (0, 1, e - c)])),
            RecurrenceTerm::new(1, 1, bivariate(&[(0, 1, c), (0, 2, -c)])),
        ],
        denominator: None,
        inhomogeneous: None,
    }
}

impl GammaLiftParent {
    /// The Eulerian polynomials `A_n`, `n >= 1`
    /// ([`crate::sequences::eulerian_polynomials_bigint`]):
    /// `A_n = (1 + (n - 1) x) A_(n-1) + x (1 - x) A'_(n-1)`, `A_1 = 1`,
    /// symmetric at `d_n = n - 1`.
    pub fn eulerian() -> Self {
        Self {
            name: "eulerian".into(),
            recurrence: RecurrenceJson::from_recurrence_rational(
                &derivative_rule(1, 0),
                1,
                &integer_rows(&[&[1]]),
                None,
            ),
            symmetry_schedule: DegreeSchedule::new(1, 0, 1),
        }
    }

    /// The type-B Eulerian polynomials `B_n`, `n >= 0`
    /// ([`crate::sequences::type_b_eulerian_polynomials_bigint`]):
    /// `B_n = (1 + (2n - 1) x) B_(n-1) + 2x (1 - x) B'_(n-1)`, `B_0 = 1`,
    /// symmetric at `d_n = n`.
    pub fn type_b_eulerian() -> Self {
        Self {
            name: "type_b_eulerian".into(),
            recurrence: RecurrenceJson::from_recurrence_rational(
                &derivative_rule(2, 1),
                0,
                &integer_rows(&[&[1]]),
                None,
            ),
            symmetry_schedule: DegreeSchedule::new(0, 0, 1),
        }
    }

    /// The Narayana polynomials `N_n = sum_k N(n, k) x^(k-1)`, `n >= 1`
    /// ([`crate::sequences::narayana_polynomials_bigint`]):
    /// `(n + 1) N_n = (2n - 1)(1 + x) N_(n-1) - (n - 2)(1 - x)^2 N_(n-2)`,
    /// `N_1 = 1`, `N_2 = 1 + x`, symmetric at `d_n = n - 1`.
    pub fn narayana() -> Self {
        let rule = Recurrence {
            terms: vec![
                RecurrenceTerm::new(
                    1,
                    0,
                    bivariate(&[(0, 0, -1), (0, 1, -1), (1, 0, 2), (1, 1, 2)]),
                ),
                RecurrenceTerm::new(
                    2,
                    0,
                    bivariate(&[
                        (0, 0, 2),
                        (0, 1, -4),
                        (0, 2, 2),
                        (1, 0, -1),
                        (1, 1, 2),
                        (1, 2, -1),
                    ]),
                ),
            ],
            denominator: Some(bivariate(&[(0, 0, 1), (1, 0, 1)])),
            inhomogeneous: None,
        };
        Self {
            name: "narayana".into(),
            recurrence: RecurrenceJson::from_recurrence_rational(
                &rule,
                1,
                &integer_rows(&[&[1], &[1, 1]]),
                None,
            ),
            symmetry_schedule: DegreeSchedule::new(1, 0, 1),
        }
    }

    /// The symmetric parts `a_n`, `n >= 1`, of the `1/k`-Eulerian
    /// polynomials, symmetric at `d_n = n - 1`, with the second-order rule of
    /// [`crate::sequences::literature::one_over_k_eulerian_symmetric_part_a_recurrence`].
    pub fn one_over_k_eulerian_symmetric_part(k: &BigRational) -> Self {
        let (rule, seeds, first_index) =
            crate::sequences::literature::one_over_k_eulerian_symmetric_part_a_recurrence(k);
        Self {
            name: format!(
                "one_over_k_eulerian_symmetric_part_a(k = {})",
                format_rational_coeff(k)
            ),
            recurrence: RecurrenceJson::from_recurrence_rational(&rule, first_index, &seeds, None),
            symmetry_schedule: DegreeSchedule::new(first_index, 0, 1),
        }
    }

    /// The default menu of [`find_gamma_lifts`]: Eulerian, type-B Eulerian
    /// and Narayana.
    pub fn standard_menu() -> Vec<Self> {
        vec![Self::eulerian(), Self::type_b_eulerian(), Self::narayana()]
    }
}

/// Options of [`find_gamma_lifts`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GammaLiftSearchOptions {
    /// The parents tried.  Default [`GammaLiftParent::standard_menu`].
    pub parents: Vec<GammaLiftParent>,
    /// Index shifts `|s| <=` this are tried (at most
    /// [`MAX_GAMMA_LIFT_INDEX_SHIFT`]).  Default 2.
    pub max_index_shift: usize,
    /// Rows matched during detection, from the input's first index (at least
    /// 2, at most [`MAX_GAMMA_LIFT_CHECKED_ROWS`]).  Default 8.
    pub detection_rows: usize,
    /// Options of the certification of each match.
    pub certify: GammaLiftOptions,
}

impl Default for GammaLiftSearchOptions {
    fn default() -> Self {
        Self {
            parents: GammaLiftParent::standard_menu(),
            max_index_shift: 2,
            detection_rows: 8,
            certify: GammaLiftOptions::default(),
        }
    }
}

/// A parent and shift that do not match the input rows.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GammaLiftRejection {
    /// The parent's name.
    pub parent: String,
    /// The index shift `s`.
    pub index_shift: i64,
    /// Why the rows do not match.
    pub reason: String,
}

/// `polytool.recurrence-gamma-lifts.v1`: the matched and certified gamma
/// lifts of an input, and the rejected parent/shift pairs.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GammaLiftReport {
    /// Always [`GAMMA_LIFTS_SCHEMA`].
    pub schema: String,
    /// The input's first index `n0`.
    pub input_first_index: usize,
    /// Number of rows matched during detection.
    pub detection_rows: usize,
    /// One certification per match (detection is finite evidence; each
    /// certification states what is proved).
    pub lifts: Vec<GammaLiftCertification>,
    /// Parent/shift pairs whose rows do not match.
    pub rejected: Vec<GammaLiftRejection>,
}

impl GammaLiftReport {
    /// Replay every certification.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        if self.schema != GAMMA_LIFTS_SCHEMA {
            return Err(RecurrenceError::InvalidInput(format!(
                "unsupported schema `{}`",
                self.schema
            )));
        }
        self.lifts
            .iter()
            .try_for_each(GammaLiftCertification::verify)
    }

    /// The certifications whose verdict is `proved`.
    pub fn proved_lifts(&self) -> impl Iterator<Item = &GammaLiftCertification> {
        self.lifts.iter().filter(|lift| lift.is_proved())
    }
}

/// `(c, lambda)` with `q = lambda g(c z)` for every row pair, or why not.
fn match_rows(input: &[Vec<Q>], gamma: &[Vec<Q>]) -> Result<(Q, Vec<Q>), String> {
    let mut scale = None;
    for (q, g) in input.iter().zip(gamma) {
        if vec_poly::is_zero(q) || vec_poly::is_zero(g) {
            return Err("a compared row is zero".into());
        }
        if q.len() != g.len() {
            return Err("the degrees differ".into());
        }
        if scale.is_none() {
            scale = (0..g.len().saturating_sub(1)).find_map(|j| {
                let nonzero = |p: &[Q]| !p[j].is_zero() && !p[j + 1].is_zero();
                (nonzero(g) && nonzero(q)).then(|| (&q[j + 1] / &g[j + 1]) / (&q[j] / &g[j]))
            });
        }
    }
    let scale = scale.unwrap_or_else(Q::one);
    let mut lambdas = Vec::with_capacity(input.len());
    for (q, g) in input.iter().zip(gamma) {
        let j = g.iter().position(|c| !c.is_zero()).expect("nonzero row");
        let mut power = Q::one();
        for _ in 0..j {
            power *= &scale;
        }
        let lambda = &q[j] / (&g[j] * &power);
        let mut expected = g.clone();
        let mut factor = lambda.clone();
        for coefficient in &mut expected {
            *coefficient *= &factor;
            factor *= &scale;
        }
        if vec_poly::trim(&expected) != vec_poly::trim(q) {
            return Err(format!("no scale c fits every row (tried c = {scale})"));
        }
        lambdas.push(lambda);
    }
    Ok((scale, lambdas))
}

/// `lambda` as a gauge: `None` for `lambda = 1`, a constant ratio otherwise.
fn gauge_of(lambdas: &[Q]) -> Option<VeroneseSectionGauge> {
    if lambdas.iter().all(One::is_one) {
        return Some(VeroneseSectionGauge::None);
    }
    let ratio = if lambdas.len() >= 2 {
        &lambdas[1] / &lambdas[0]
    } else {
        Q::one()
    };
    lambdas
        .windows(2)
        .all(|pair| &pair[1] / &pair[0] == ratio)
        .then(|| VeroneseSectionGauge::Ratio {
            numerator: vec![format_rational_coeff(&ratio)],
            denominator: vec!["1".into()],
            initial_value: format_rational_coeff(&lambdas[0]),
        })
}

/// Find gamma lifts `Q_n = lambda(n) Gamma_(n+s)(c z)` of the input rows
/// (`input_rows[0]` is `Q_(first_index)`; at least the identifying rows)
/// among the parents of `options` and certify every match with
/// [`certify_gamma_lift`] (module documentation of `gamma_lift`).
///
/// Detection compares `options.detection_rows` rows exactly; it reads the
/// scale `c` from two consecutive nonzero coefficients and accepts a factor
/// `lambda(n)` that is `1` or has a constant ratio `lambda(n)/lambda(n-1)`.
/// A match is a candidate only; each certification states what is proved.
///
/// # Errors
///
/// [`RecurrenceError`] for an invalid input rule or rows, options over the
/// limits, and errors of the certification of a match.
///
/// ```
/// use polytool::recurrence::*;
/// use polytool::sequences::literature::peak_polynomials_recurrence;
///
/// let (rule, seeds, first_index) = peak_polynomials_recurrence();
/// let report =
///     find_gamma_lifts(&rule, &seeds, first_index, &GammaLiftSearchOptions::default())?;
/// let lift = report.proved_lifts().next().expect("the Eulerian parent");
/// assert_eq!(lift.candidate.parent.name, "eulerian");
/// assert_eq!(lift.candidate.scale, "1/4");
/// report.verify()?;
/// # Ok::<(), RecurrenceError>(())
/// ```
pub fn find_gamma_lifts(
    input: &Recurrence,
    input_rows: &[Vec<BigRational>],
    first_index: usize,
    options: &GammaLiftSearchOptions,
) -> Result<GammaLiftReport, RecurrenceError> {
    if options.max_index_shift > MAX_GAMMA_LIFT_INDEX_SHIFT
        || !(2..=MAX_GAMMA_LIFT_CHECKED_ROWS).contains(&options.detection_rows)
    {
        return Err(RecurrenceError::LimitExceeded(format!(
            "max_index_shift <= {MAX_GAMMA_LIFT_INDEX_SHIFT} and 2 <= detection_rows <= \
             {MAX_GAMMA_LIFT_CHECKED_ROWS} are required"
        )));
    }
    let n0 = first_index;
    let prepared = prepare_rule(input, input_rows, n0, "input")?;
    let count = options.detection_rows.max(prepared.seeds.len());
    let rows = prepared.rows(count)?;
    let mut report = GammaLiftReport {
        schema: GAMMA_LIFTS_SCHEMA.into(),
        input_first_index: n0,
        detection_rows: count,
        lifts: Vec::new(),
        rejected: Vec::new(),
    };
    let max_shift = options.max_index_shift as i64;
    for parent in &options.parents {
        let (rule, p0, seeds) = parent.recurrence.to_recurrence_parts()?;
        let prepared_parent = prepare_rule(&rule, &seeds, p0, "parent")?;
        let gamma = RecurrenceTransformation::Gamma {
            schedule: parent.symmetry_schedule.clone(),
        };
        for shift in -max_shift..=max_shift {
            let reject = |reason: String| GammaLiftRejection {
                parent: parent.name.clone(),
                index_shift: shift,
                reason,
            };
            let start = n0 as i64 + shift;
            if start < p0 as i64 {
                report.rejected.push(reject(format!(
                    "n0 + s = {start} precedes the parent's first index {p0}"
                )));
                continue;
            }
            let start = start as usize;
            let parent_rows = prepared_parent.rows(start - p0 + count)?;
            let gamma_rows: Result<Vec<Vec<Q>>, String> = (0..count)
                .map(|i| {
                    let m = start + i;
                    gamma
                        .map_row(m, &parent_rows[start - p0 + i])
                        .map(|mut components| vec_poly::trim(&components.swap_remove(0)))
                        .map_err(|error| format!("parent row {m}: {error}"))
                })
                .collect();
            let gamma_rows = match gamma_rows {
                Ok(rows) => rows,
                Err(reason) => {
                    report.rejected.push(reject(reason));
                    continue;
                }
            };
            let (scale, lambdas) = match match_rows(&rows, &gamma_rows) {
                Ok(found) => found,
                Err(reason) => {
                    report.rejected.push(reject(reason));
                    continue;
                }
            };
            let Some(gauge) = gauge_of(&lambdas) else {
                report.rejected.push(reject(
                    "the factors lambda(n) are neither 1 nor of constant ratio".into(),
                ));
                continue;
            };
            let candidate = GammaLiftCandidate {
                parent: parent.clone(),
                index_shift: shift,
                scale: format_rational_coeff(&scale),
                gauge,
            };
            report.lifts.push(certify_gamma_lift(
                input,
                input_rows,
                n0,
                &candidate,
                &options.certify,
            )?);
        }
    }
    Ok(report)
}
