//! The coordinates of the iterated-preserver theorem: the rows `Q_n` the
//! induction runs on, with the step operator they satisfy.
//!
//! - **Identity**: `Q_n = P_n`.
//! - **Power gauge** `x^e`: `Q_n = P_n / x^e` when `x^e` divides every start
//!   row; then `Q_n = T'_n(Q_(n-L))` with `T' = x^(-e) T x^e` (a polynomial
//!   operator when `x^e` divides every coefficient of `T x^e`), and
//!   `P_n = x^e Q_n` has the roots of `Q_n` and the root `0`, which lies in
//!   every certified region.
//! - **Translation** by `r`: `Q_n(x) = P_n(x + r)`, `T'` with the
//!   coefficients `c_d(n, x + r)` (`D` commutes with translations); the roots
//!   of `P_n` are those of `Q_n` plus `r`, so a ray of `Q` gives the same
//!   ray of `P` when `r` lies in it, and the line otherwise.  The shifts
//!   tried are the nonzero rational roots of the factor free of `n` of the
//!   coefficient of the highest derivative.
//! - **Reversal**: `Q_n = x^(d_n) P_n(1/x)` with the nominal degree
//!   `d_n = deg P_s + delta (n - s)/L` of the residue class of its start row
//!   `s` (the same linear function `d(n)` of `n` for every class).  Writing
//!   `T = sum_j x^j f_j(n, theta)` (`theta = x D`) with `delta` the largest
//!   `j`, the reversal of `T(x^k)` at the nominal degree is
//!   `T' = sum_j x^(delta - j) f_j(n, d(n) - theta)`, where `d(n)` is the
//!   nominal degree of `P_(n-L)`; this identity holds whatever the actual
//!   degrees, since `T` raises the degree by at most `delta`.  `P_n` is zero
//!   exactly when `Q_n` is, and its roots are `0` and the inverses of the
//!   nonzero roots of `Q_n`, so every certified region transfers.
//!
//! Every transform is an exact identity in the Ore algebra, recomputed by
//! [`IteratedPreserverCertificate::verify_operator_identity`] from the
//! recorded coordinates.

use super::*;
use crate::recurrence::implicit::algebra::{from_theta, theta_add, ThetaOp};

/// Largest number of translations tried.
const MAX_TRANSLATIONS: usize = 3;

/// Largest power of a power gauge tried.
const MAX_GAUGE_POWER: usize = 6;

/// The coordinates `Q_n` in which [`check_iterated_preserver_theorem`] runs
/// its induction (module documentation, "Coordinates").
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum IteratedPreserverCoordinates {
    /// `Q_n = P_n`.
    #[default]
    Identity,
    /// `Q_n = P_n / x^power`, with the step `x^(-power) T_n x^power`.
    PowerGauge {
        /// `e >= 1`.
        power: usize,
    },
    /// `Q_n(x) = P_n(x + shift)`, with the step coefficients
    /// `c_d(n, x + shift)`.
    Translation {
        /// The shift `r`, an exact rational string.
        shift: String,
    },
    /// `Q_n = x^(d_n) P_n(1/x)` at the nominal degree `d_n`, with the step
    /// `sum_j x^(delta - j) f_j(n, d(n) - theta)`.
    Reversal {
        /// `d(n)`, the nominal degree of `P_(n-L)` as a polynomial in `n`
        /// (ascending exact rational strings).
        nominal_previous_degree: Vec<String>,
    },
}

impl IteratedPreserverCoordinates {
    /// Whether these are the identity coordinates `Q_n = P_n`.
    pub fn is_identity(&self) -> bool {
        *self == Self::Identity
    }

    fn translation_shift(&self) -> Option<Q> {
        match self {
            Self::Translation { shift } => parse_rational_coeff(shift).ok(),
            _ => None,
        }
    }

    /// The step `T'` of the rows `Q_n` from the step `T` of the rows `P_n`
    /// (both as `f(n) T = sum_d c_d(n, x) D^d` over the same `f`), or `None`
    /// when the transform is not a polynomial operator.
    pub(super) fn transform_step(&self, terms: &StepCoefficients) -> Option<StepCoefficients> {
        let op = ore_operator(terms);
        let transformed = match self {
            Self::Identity => op,
            Self::PowerGauge { power } => {
                let gauge = DOp::single(0, 0, Biv::monomial(0, *power, Q::one()));
                let composed = op.compose(&gauge);
                if composed.variable_valuation() < *power {
                    return None;
                }
                composed.divide_by_variable_power(*power)
            }
            Self::Translation { .. } => {
                let shift = self.translation_shift()?;
                let mut out = DOp::default();
                for (key, c) in &op.0 {
                    out.add_term(*key, c.shift(&Q::zero(), &shift));
                }
                out
            }
            Self::Reversal {
                nominal_previous_degree,
            } => {
                let degree = nominal_previous_degree
                    .iter()
                    .map(|value| parse_rational_coeff(value).ok())
                    .collect::<Option<Vec<Q>>>()?;
                let delta = degree_step(terms).filter(|&delta| delta >= 0)?;
                let (top, euler) = to_theta(&op);
                let mut reversed = ThetaOp::new();
                for ((shift, power), a) in &euler {
                    if *shift != 0 {
                        return None;
                    }
                    let level = *power as i64 - top as i64;
                    let exponent = usize::try_from(delta - level).ok()?;
                    theta_add(&mut reversed, (0, exponent), reflected_theta(a, &degree));
                }
                from_theta(&reversed)
            }
        };
        let mut out = StepCoefficients::new();
        for ((shift, order), c) in &transformed.0 {
            if *shift != 0 {
                return None;
            }
            if !c.is_zero() {
                out.insert(*order, c.rows().to_vec());
            }
        }
        (!out.is_empty()).then_some(out)
    }

    /// A start row `P_s` in these coordinates (`degree` its nominal degree
    /// for the reversal), or `None` when it is not a polynomial.
    pub(super) fn transform_row(&self, row: &[Q]) -> Option<Vec<Q>> {
        let row = vec_poly::trim(row);
        match self {
            Self::Identity => Some(row),
            Self::PowerGauge { power } => {
                if vec_poly::is_zero(&row) {
                    return Some(row);
                }
                row.iter()
                    .take(*power)
                    .all(Zero::is_zero)
                    .then(|| row[*power..].to_vec())
            }
            Self::Translation { .. } => {
                let shift = self.translation_shift()?;
                Some(vec_poly::trim(&vec_poly::compose(&row, &[shift, Q::one()])))
            }
            Self::Reversal { .. } => {
                let degree = vec_poly::degree(&row).unwrap_or(0);
                vec_poly::reverse(&row, degree).map(|q| vec_poly::trim(&q))
            }
        }
    }

    /// The region certified for the rows `P_n` when the rows `Q_n` have
    /// every root in `region`.
    pub(super) fn claimed_region(&self, region: ClosedRegion) -> ClosedRegion {
        let Some(shift) = self.translation_shift() else {
            return region;
        };
        let keeps = match region {
            ClosedRegion::NonpositiveRay => !shift.is_positive(),
            ClosedRegion::NonnegativeRay => !shift.is_negative(),
            ClosedRegion::RealLine => true,
        };
        if keeps {
            region
        } else {
            ClosedRegion::RealLine
        }
    }

    /// `Q_n` in terms of `P_n`, for the normal form and the conclusion.
    pub(super) fn text(&self) -> Option<String> {
        match self {
            Self::Identity => None,
            Self::PowerGauge { power } => Some(format!("Q_n = P_n / x^{power}")),
            Self::Translation { shift } => Some(format!("Q_n(x) = P_n(x + {shift})")),
            Self::Reversal {
                nominal_previous_degree,
            } => {
                let degree = nominal_previous_degree
                    .iter()
                    .filter_map(|value| parse_rational_coeff(value).ok())
                    .collect::<Vec<_>>();
                Some(format!(
                    "Q_n = x^(d_n) P_n(1/x) at the nominal degree d_n, with d_(n-L) = {} (an \
                     upper bound of deg P_(n-L))",
                    format_index_polynomial(&degree, "n")
                ))
            }
        }
    }
}

/// `a(n, d(n) - theta)` for `a` in `(n, theta)`.
fn reflected_theta(a: &Biv, degree: &[Q]) -> Biv {
    let mut rows: Vec<Vec<Q>> = degree.iter().map(|c| vec![c.clone()]).collect();
    if rows.is_empty() {
        rows.push(vec![Q::zero()]);
    }
    rows[0].push(-Q::one());
    let base = Biv::new(rows);
    let mut out = Biv::zero();
    for (i, row) in a.rows().iter().enumerate() {
        for (m, value) in row.iter().enumerate() {
            if value.is_zero() {
                continue;
            }
            out = out.add(&Biv::monomial(i, 0, value.clone()).multiply(&base.power(m)));
        }
    }
    out
}

/// The coordinates tried after the identity: power gauges up to the common
/// power of `x` of the start rows, translations to the rational roots of the
/// highest derivative coefficient, and the reversal when the classes share
/// one nominal degree function (`start_rows[c]` is `P_(n_first - lag + c)`).
pub(super) fn coordinate_candidates(
    terms: &StepCoefficients,
    lag: usize,
    n_first: usize,
    start_rows: &[Vec<Q>],
) -> Vec<IteratedPreserverCoordinates> {
    let mut out = Vec::new();
    let valuation = start_rows
        .iter()
        .map(|row| {
            let row = vec_poly::trim(row);
            if vec_poly::is_zero(&row) {
                0
            } else {
                row.iter().take_while(|c| c.is_zero()).count()
            }
        })
        .min()
        .unwrap_or(0);
    out.extend(
        (1..=valuation.min(MAX_GAUGE_POWER))
            .map(|power| IteratedPreserverCoordinates::PowerGauge { power }),
    );
    if let Some(top) = terms.keys().copied().max().filter(|&top| top >= 1) {
        out.extend(
            rational_roots_of_index_free_factor(&terms[&top])
                .into_iter()
                .filter(|root| !root.is_zero())
                .take(MAX_TRANSLATIONS)
                .map(|root| IteratedPreserverCoordinates::Translation {
                    shift: format_rational_coeff(&root),
                }),
        );
    }
    if let Some(degree) = nominal_previous_degree(terms, lag, n_first, start_rows) {
        out.push(IteratedPreserverCoordinates::Reversal {
            nominal_previous_degree: degree.iter().map(format_rational_coeff).collect(),
        });
    }
    out
}

/// `d(n) = deg P_s + delta ((n - s)/L - 1)`, the same for every class start
/// `s`, or `None` (a lowering step, a zero start row or classes with
/// different degree functions).
fn nominal_previous_degree(
    terms: &StepCoefficients,
    lag: usize,
    n_first: usize,
    start_rows: &[Vec<Q>],
) -> Option<Vec<Q>> {
    let delta = degree_step(terms).filter(|&delta| delta >= 0)?;
    let lag_q = Q::from_integer(BigInt::from(lag));
    let slope = Q::from_integer(BigInt::from(delta)) / &lag_q;
    let mut constant: Option<Q> = None;
    for (c, row) in start_rows.iter().enumerate() {
        let degree = vec_poly::degree(&vec_poly::trim(row))?;
        if vec_poly::is_zero(&vec_poly::trim(row)) {
            return None;
        }
        let start = Q::from_integer(BigInt::from(n_first - lag + c));
        // deg P_s + delta ((n - s)/L - 1) = (deg P_s - slope (s + L)) + slope n.
        let value = Q::from_integer(BigInt::from(degree)) - &slope * (start + &lag_q);
        match &constant {
            None => constant = Some(value),
            Some(existing) if *existing == value => {}
            Some(_) => return None,
        }
    }
    Some(vec_poly::trim(&[constant?, slope]))
}
