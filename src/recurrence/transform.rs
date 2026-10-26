//! Exact coordinate transformations of polynomial sequences and transport of
//! derivative recurrences into the new coordinates.
//!
//! A [`RecurrenceTransformation`] is an exact, invertible change of
//! coordinates `P_n -> Q_n` (or `P_n -> (a_n, b_n)`) of a polynomial
//! sequence.  The registry contains
//!
//! - multiplication by `x^r` and its inverse, removal of `x^r` (the latter
//!   requires `x^r | P_n`);
//! - reversal `Q_n = x^(d_n) P_n(1/x)` at an explicit degree bound;
//! - rational affine changes `Q_n = P_n(a x + b)`;
//! - nonsingular Möbius changes `Q_n = (c x + e)^(d_n) P_n((a x + b)/(c x + e))`
//!   with `a e - b c != 0`, and the `f`-polynomial
//!   `(1 + x)^(d_n) P_n(x/(1 + x))` computed by
//!   [`crate::decomposition::f_polynomial`] at the declared bound `d_n`;
//! - gamma coordinates at a symmetry index `d_n` (requires
//!   `x^(d_n) P_n(1/x) = P_n`);
//! - the symmetric decomposition `P_n = a_n + x b_n` with `a_n` symmetric at
//!   index `d_n` and `b_n` symmetric at index `d_n - 1`
//!   ([`crate::Polynomial::stapledon_decomposition`]), a two-component vector.
//!
//! Every map records its exact row map, its inverse, and the image or
//! precondition that makes the inverse exact.  Degree bounds come from an
//! explicit affine [`DegreeSchedule`] `d_n = d_0 + h (n - n_0)` on an index
//! domain `n_0 <= n (<= n_1)`; the observed degree is never substituted.
//!
//! # Transport
//!
//! A scalar rule in the [`Recurrence`] convention (target index `m`)
//!
//! ```text
//! f(m,x) P_m = sum s(m) c(m,x) D^k P_(m-r) + g(m,x)
//! ```
//!
//! is transported by [`RecurrenceTransformation::transport`].  For a Möbius
//! map put `phi = (a x + b)/ell`, `ell = c x + e` and `Delta = a e - b c`.
//! Then `P_j^(k)(phi) = Delta^(-k) ell^(k - d_j) T_k[Q_j]` with `T_0 = Q` and
//! `T_(k+1) = ell T_k' + (k - d_j) c T_k`; for `x^r` multiplication the same
//! recursion holds with `ell = x`, `Delta = 1` and exponents `-r - k`.
//! Substituting, multiplying by `ell^(d_m)` and by the smallest power of
//! `ell` that clears every denominator gives a polynomial identity; a common
//! power of `ell` is then divided out exactly.  Nothing is evaluated at the
//! pole of `phi`.  For the first-order rule
//! `R(n) P_(n+1) = A P_n + B P_n'` this is the formula
//! `R Q_(n+1) = ell^h (A(phi) - d_n c ell B(phi)/Delta) Q_n
//!  + ell^(h+2) B(phi)/Delta Q_n'`.
//!
//! The transported rule is an exact algebraic consequence of the source rule,
//! valid wherever the source rule holds; it is not refitted.  An additive
//! term `g` picks up the factor `ell^(d_m)`, so it is transported only when
//! that factor is index independent (`h = 0`, or a constant `ell` with
//! `ell^h = 1`).  The symmetric decomposition transports a first-order rule
//! without alternating signs to the complete two-component
//! [`SharedDenominatorVectorRecurrence`]; no component is selected.  Gamma
//! coordinates are a quadratic substitution and are not transported; search
//! the gamma rows directly.
//!
//! # Reports
//!
//! [`transform_recurrence`] applies a transformation to rows given directly,
//! generated from a defining `polytool.recurrence.v1` rule with seeds, or
//! taken from the bundled OEIS catalog.  The versioned
//! [`RecurrenceTransformReport`] keeps the request, the transformation
//! record, source and transformed rows, exact reconstruction, the status of
//! each image condition (`certified_all_indices` when the transported rule
//! is polynomial and the seeds qualify, otherwise `checked_supplied_rows`),
//! the transported rule (`transported_identity`) and an optional independent
//! search in the new coordinates (`finite_fit`).  The two outputs are kept
//! separately.  [`verify_recurrence_transform`] replays a report exactly.
//!
//! Root-location properties of the maps are not part of the algebraic
//! report; the opt-in `root_domain` request field attaches them separately
//! (see [`RootDomainReport`]).  Factorial rescalings and implicit differential relations are separate
//! extensions of this registry.

use super::*;
use crate::decomposition::f_polynomial;
use crate::polynomial::Polynomial;
use crate::real_rootedness::gamma_coefficients_bigint_coeffs;
use crate::vec_poly;
use num_integer::Integer;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt;

type Q = BigRational;

/// Schema tag of [`RecurrenceTransformReport`].
pub const RECURRENCE_TRANSFORM_SCHEMA: &str = "polytool.recurrence-transform.v1";

/// Evidence label of a transported rule: an exact consequence of the source rule.
pub const TRANSPORTED_IDENTITY_EVIDENCE: &str = "transported_identity";

/// Evidence label of an independent search in the new coordinates.
pub const FINITE_FIT_EVIDENCE: &str = "finite_fit";

/// Condition status: proved for every index of the domain from the source rule.
pub const CERTIFIED_ALL_INDICES: &str = "certified_all_indices";

/// Condition status: verified exactly on every supplied row only.
pub const CHECKED_SUPPLIED_ROWS: &str = "checked_supplied_rows";

// ---------------------------------------------------------------------------
// Degree schedules
// ---------------------------------------------------------------------------

/// An affine degree schedule `d_n = first_degree + step (n - first_index)` on
/// the index domain `first_index <= n` (`<= last_index` when given).
///
/// Consecutive bounds satisfy `d_(n+1) = d_n + step`.  A decreasing schedule
/// must name its last index, and every degree on the domain must be
/// nonnegative.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct DegreeSchedule {
    /// First index `n_0` of the domain (original sequence index).
    pub first_index: usize,
    /// Degree bound `d_(n_0)` at the first index.
    pub first_degree: usize,
    /// Increment `h = d_(n+1) - d_n`. Defaults to zero.
    #[serde(default)]
    pub step: i64,
    /// Last index of the domain (inclusive); `None` means unbounded above.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub last_index: Option<usize>,
}

impl DegreeSchedule {
    /// Unbounded schedule `d_n = first_degree + step (n - first_index)`.
    pub fn new(first_index: usize, first_degree: usize, step: i64) -> Self {
        Self {
            first_index,
            first_degree,
            step,
            last_index: None,
        }
    }

    /// Return this schedule restricted to indices `<= last_index`.
    pub fn with_last_index(mut self, last_index: usize) -> Self {
        self.last_index = Some(last_index);
        self
    }

    fn raw_degree(&self, index: usize) -> i128 {
        self.first_degree as i128 + self.step as i128 * (index as i128 - self.first_index as i128)
    }

    /// Reject empty domains, negative degrees and unbounded decreasing schedules.
    pub fn validate(&self) -> Result<(), TransformError> {
        match self.last_index {
            Some(last) if last < self.first_index => Err(TransformError::InvalidSchedule(format!(
                "last_index {last} precedes first_index {}",
                self.first_index
            ))),
            Some(last) if self.raw_degree(last) < 0 => {
                Err(TransformError::InvalidSchedule(format!(
                    "degree bound {} at index {last} is negative",
                    self.raw_degree(last)
                )))
            }
            None if self.step < 0 => Err(TransformError::InvalidSchedule(
                "a decreasing degree schedule needs an explicit last_index".into(),
            )),
            _ => {
                self.weight_constant()?;
                Ok(())
            }
        }
    }

    /// Whether `index` lies in the declared domain.
    pub fn contains(&self, index: usize) -> bool {
        index >= self.first_index && self.last_index.is_none_or(|last| index <= last)
    }

    /// The degree bound `d_index`, or an error outside the domain.
    pub fn degree_at(&self, index: usize) -> Result<usize, TransformError> {
        if !self.contains(index) {
            return Err(TransformError::IndexOutsideDomain { index });
        }
        let raw = self.raw_degree(index);
        if raw < 0 {
            return Err(TransformError::InvalidSchedule(format!(
                "degree bound {raw} at index {index} is negative"
            )));
        }
        usize::try_from(raw).map_err(|_| TransformError::Limit("degree bound overflow".into()))
    }

    /// `w` with `d_n = w + step n` for every index `n`.
    fn weight_constant(&self) -> Result<i64, TransformError> {
        let value = self.first_degree as i128 - self.step as i128 * self.first_index as i128;
        i64::try_from(value).map_err(|_| TransformError::Limit("degree schedule overflow".into()))
    }

    fn formula(&self) -> String {
        format!(
            "d_n = {} + {}*(n - {})",
            self.first_degree, self.step, self.first_index
        )
    }

    fn domain_text(&self) -> String {
        match self.last_index {
            Some(last) => format!("{} <= n <= {last}", self.first_index),
            None => format!("n >= {}", self.first_index),
        }
    }
}

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Why a transformation, inverse or transport was rejected.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TransformError {
    /// A map parameter is malformed (for example an unparsable rational).
    InvalidParameter(String),
    /// The Möbius or affine map is singular (`a e - b c = 0`, or `a = 0`).
    SingularMobius,
    /// The degree schedule is empty, negative somewhere, or unbounded while decreasing.
    InvalidSchedule(String),
    /// A row index lies outside the schedule's domain.
    IndexOutsideDomain {
        /// The offending original index.
        index: usize,
    },
    /// A row has larger degree than the declared bound.
    DegreeExceedsBound {
        /// Original index of the row.
        index: usize,
        /// Actual degree of the row.
        degree: usize,
        /// Declared bound `d_index`.
        bound: usize,
    },
    /// Removal of `x^power` was requested but `x^power` does not divide the row.
    NotDivisible {
        /// Original index of the row.
        index: usize,
        /// Requested power of `x`.
        power: usize,
    },
    /// A row is not symmetric at the declared symmetry index.
    NotSymmetric {
        /// Original index of the row.
        index: usize,
        /// Declared symmetry index `d_index`.
        symmetry_index: usize,
    },
    /// An inverse was applied to data outside the image of the map.
    ImageCondition {
        /// Original index of the row.
        index: usize,
        /// The violated condition.
        reason: String,
    },
    /// A transformed state has the wrong number of components.
    ComponentCount {
        /// Number of components of the map's image.
        expected: usize,
        /// Number supplied.
        found: usize,
    },
    /// The additive term would acquire an index-dependent power of `ell`.
    InhomogeneousNotTransportable,
    /// This transformation has no transport for the given rule.
    TransportUnsupported(String),
    /// A resource limit was exceeded.
    Limit(String),
}

impl fmt::Display for TransformError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidParameter(message) => write!(f, "invalid map parameter: {message}"),
            Self::SingularMobius => write!(f, "singular map: a e - b c (or a) is zero"),
            Self::InvalidSchedule(message) => write!(f, "invalid degree schedule: {message}"),
            Self::IndexOutsideDomain { index } => {
                write!(f, "index {index} lies outside the declared index domain")
            }
            Self::DegreeExceedsBound {
                index,
                degree,
                bound,
            } => write!(
                f,
                "row {index} has degree {degree} above the declared bound {bound}"
            ),
            Self::NotDivisible { index, power } => {
                write!(f, "x^{power} does not divide row {index}")
            }
            Self::NotSymmetric {
                index,
                symmetry_index,
            } => write!(
                f,
                "row {index} is not symmetric at symmetry index {symmetry_index}"
            ),
            Self::ImageCondition { index, reason } => {
                write!(f, "row {index} is outside the image of the map: {reason}")
            }
            Self::ComponentCount { expected, found } => {
                write!(f, "expected {expected} components, found {found}")
            }
            Self::InhomogeneousNotTransportable => write!(
                f,
                "the additive term would acquire an index-dependent power of ell; \
                 only h = 0 or a constant ell with ell^h = 1 is supported"
            ),
            Self::TransportUnsupported(message) => write!(f, "transport unsupported: {message}"),
            Self::Limit(message) => write!(f, "limit exceeded: {message}"),
        }
    }
}

impl std::error::Error for TransformError {}

// ---------------------------------------------------------------------------
// The transformation registry
// ---------------------------------------------------------------------------

/// An exact change of polynomial coordinates `P_n -> Q_n`.
///
/// Rational parameters are integer or `num/den` strings.  The JSON form is
/// tagged by `kind` (`multiply_power`, `divide_power`, `reversal`, `affine`,
/// `mobius`, `f_polynomial`, `gamma`, `symmetric_decomposition`).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum RecurrenceTransformation {
    /// `Q_n = x^power P_n`; inverse requires `x^power | Q_n`.
    MultiplyPower {
        /// The power `r`.
        power: usize,
    },
    /// `Q_n = P_n / x^power`; requires `x^power | P_n`; inverse multiplies back.
    DividePower {
        /// The power `r`.
        power: usize,
    },
    /// `Q_n = x^(d_n) P_n(1/x)`; requires `deg P_n <= d_n`; self-inverse.
    Reversal {
        /// Degree bounds `d_n`.
        schedule: DegreeSchedule,
    },
    /// `Q_n = P_n(a x + b)` with `a != 0`; inverse `P_n(y) = Q_n((y - b)/a)`.
    Affine {
        /// Rational `a`.
        a: String,
        /// Rational `b`.
        b: String,
    },
    /// `Q_n = (c x + e)^(d_n) P_n((a x + b)/(c x + e))` with `a e - b c != 0`;
    /// requires `deg P_n <= d_n`; inverse
    /// `P_n(y) = (a e - b c)^(-d_n) (a - c y)^(d_n) Q_n((e y - b)/(a - c y))`.
    Mobius {
        /// Rational `a`.
        a: String,
        /// Rational `b`.
        b: String,
        /// Rational `c`.
        c: String,
        /// Rational `e`.
        e: String,
        /// Degree bounds `d_n`.
        schedule: DegreeSchedule,
    },
    /// The `f`-polynomial `Q_n = (1 + x)^(d_n) P_n(x/(1 + x))`, the Möbius map
    /// `(a, b, c, e) = (1, 0, 1, 1)`; inverse `P_n(y) = (1 - y)^(d_n) Q_n(y/(1 - y))`.
    FPolynomial {
        /// Degree bounds `d_n`.
        schedule: DegreeSchedule,
    },
    /// Gamma coordinates `P_n = sum_i Q_(n,i) x^i (1 + x)^(d_n - 2i)`; requires
    /// `x^(d_n) P_n(1/x) = P_n`; inverse requires `deg Q_n <= floor(d_n/2)`.
    Gamma {
        /// Symmetry indices `d_n`.
        schedule: DegreeSchedule,
    },
    /// `P_n = a_n + x b_n` with `a_n` symmetric at index `d_n` and `b_n` at
    /// index `d_n - 1`; requires `deg P_n <= d_n`; the image is the vector
    /// `(a_n, b_n)`, and the inverse requires both symmetries.
    SymmetricDecomposition {
        /// Degree bounds and symmetry indices `d_n`.
        schedule: DegreeSchedule,
    },
}

/// The transported form of a rule.
#[derive(Debug, Clone)]
pub enum TransportedRecurrence {
    /// A scalar rule in the same target-index convention as the source rule.
    Scalar(Recurrence),
    /// A first-order vector system `q(n,x) F_(n+1) = M(n,x,D) F_n + G(n,x)`
    /// in the source index `n`.
    Vector(SharedDenominatorVectorRecurrence),
}

/// One component of a transformed state.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct TransformComponent {
    /// Component name (`Q`, `a` or `b`).
    pub name: String,
    /// Symmetry index of the component (`d_n`, `d_n - 1`), if it is symmetric.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub symmetry_index: Option<String>,
}

/// Exact description of a transformation: row map, inverse, image
/// conditions and conventions.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct TransformationRecord {
    /// The `kind` tag of the transformation.
    pub kind: String,
    /// The exact row map.
    pub row_map: String,
    /// The exact inverse (reconstruction of `P_n`).
    pub inverse: String,
    /// Preconditions of the map and image conditions of the inverse.
    pub image_conditions: Vec<String>,
    /// Degree schedule formula, if the map uses one.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub degree_schedule: Option<String>,
    /// Index domain on which the map is defined.
    pub index_domain: String,
    /// Components of the image, in state order.
    pub components: Vec<TransformComponent>,
    /// How rules are transported, or why they are not.
    pub transport: String,
}

fn parse_q(text: &str, name: &str) -> Result<Q, TransformError> {
    parse_rational_coeff(text).map_err(|e| TransformError::InvalidParameter(format!("{name}: {e}")))
}

fn q(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

/// Möbius parameters `(a, b, c, e)` and schedule of a Möbius-type map.
struct MobiusParameters<'a> {
    a: Q,
    b: Q,
    c: Q,
    e: Q,
    schedule: &'a DegreeSchedule,
}

impl MobiusParameters<'_> {
    fn delta(&self) -> Q {
        &self.a * &self.e - &self.b * &self.c
    }
}

impl RecurrenceTransformation {
    /// The `kind` tag.
    pub fn kind_name(&self) -> &'static str {
        match self {
            Self::MultiplyPower { .. } => "multiply_power",
            Self::DividePower { .. } => "divide_power",
            Self::Reversal { .. } => "reversal",
            Self::Affine { .. } => "affine",
            Self::Mobius { .. } => "mobius",
            Self::FPolynomial { .. } => "f_polynomial",
            Self::Gamma { .. } => "gamma",
            Self::SymmetricDecomposition { .. } => "symmetric_decomposition",
        }
    }

    /// The degree schedule, if the map uses one.
    pub fn schedule(&self) -> Option<&DegreeSchedule> {
        match self {
            Self::Reversal { schedule }
            | Self::Mobius { schedule, .. }
            | Self::FPolynomial { schedule }
            | Self::Gamma { schedule }
            | Self::SymmetricDecomposition { schedule } => Some(schedule),
            _ => None,
        }
    }

    /// Number of components of a transformed state.
    pub fn component_count(&self) -> usize {
        match self {
            Self::SymmetricDecomposition { .. } => 2,
            _ => 1,
        }
    }

    /// Whether `index` lies in the map's index domain.
    pub fn contains_index(&self, index: usize) -> bool {
        self.schedule().is_none_or(|s| s.contains(index))
    }

    fn mobius_parameters(&self) -> Result<Option<MobiusParameters<'_>>, TransformError> {
        Ok(match self {
            Self::Reversal { schedule } => Some(MobiusParameters {
                a: q(0),
                b: q(1),
                c: q(1),
                e: q(0),
                schedule,
            }),
            Self::FPolynomial { schedule } => Some(MobiusParameters {
                a: q(1),
                b: q(0),
                c: q(1),
                e: q(1),
                schedule,
            }),
            Self::Mobius {
                a,
                b,
                c,
                e,
                schedule,
            } => Some(MobiusParameters {
                a: parse_q(a, "a")?,
                b: parse_q(b, "b")?,
                c: parse_q(c, "c")?,
                e: parse_q(e, "e")?,
                schedule,
            }),
            _ => None,
        })
    }

    fn affine_parameters(&self) -> Result<Option<(Q, Q)>, TransformError> {
        match self {
            Self::Affine { a, b } => Ok(Some((parse_q(a, "a")?, parse_q(b, "b")?))),
            _ => Ok(None),
        }
    }

    /// Check parameters and schedule; singular maps and invalid schedules are rejected.
    pub fn validate(&self) -> Result<(), TransformError> {
        if let Some(schedule) = self.schedule() {
            schedule.validate()?;
        }
        if let Some((a, _)) = self.affine_parameters()? {
            if a.is_zero() {
                return Err(TransformError::SingularMobius);
            }
        }
        if let Some(mobius) = self.mobius_parameters()? {
            if mobius.delta().is_zero() {
                return Err(TransformError::SingularMobius);
            }
        }
        Ok(())
    }

    /// Apply the map to the row `P_index`, returning one coefficient vector per
    /// component (ascending powers, zero as `[0]`).
    ///
    /// Preconditions (degree bound, divisibility, symmetry) are checked exactly
    /// on this row and reported as errors; they are never assumed.
    pub fn map_row(&self, index: usize, row: &[Q]) -> Result<Vec<Vec<Q>>, TransformError> {
        self.validate()?;
        let row = vec_poly::trim(row);
        match self {
            Self::MultiplyPower { power } => Ok(vec![vec_poly::shift_by(&row, *power)]),
            Self::DividePower { power } => Ok(vec![divide_power(index, &row, *power)?]),
            Self::Affine { .. } => {
                let (a, b) = self.affine_parameters()?.expect("affine");
                Ok(vec![vec_poly::compose(&row, &[b, a])])
            }
            Self::Reversal { schedule } => {
                let d = schedule.degree_at(index)?;
                check_degree(index, &row, d)?;
                Ok(vec![vec_poly::reverse(&row, d).expect("degree checked")])
            }
            Self::FPolynomial { schedule } => {
                let d = schedule.degree_at(index)?;
                check_degree(index, &row, d)?;
                // Integer rows use BigInt arithmetic; both routes are exact.
                let mapped = match integer_coefficients(&row) {
                    Some(integers) => f_polynomial(&Polynomial::new(integers), d)
                        .expect("degree checked")
                        .coeffs()
                        .iter()
                        .cloned()
                        .map(Q::from_integer)
                        .collect::<Vec<_>>(),
                    None => f_polynomial(&Polynomial::new(row), d)
                        .expect("degree checked")
                        .coeffs()
                        .to_vec(),
                };
                Ok(vec![vec_poly::trim(&mapped)])
            }
            Self::Mobius { .. } => {
                let m = self.mobius_parameters()?.expect("mobius");
                let d = m.schedule.degree_at(index)?;
                check_degree(index, &row, d)?;
                Ok(vec![mobius_substitute(&row, &m.a, &m.b, &m.c, &m.e, d)])
            }
            Self::Gamma { schedule } => {
                let d = schedule.degree_at(index)?;
                Ok(vec![gamma_at_index(index, &row, d)?])
            }
            Self::SymmetricDecomposition { schedule } => {
                let d = schedule.degree_at(index)?;
                check_degree(index, &row, d)?;
                let (a, b) = Polynomial::new(row)
                    .stapledon_decomposition(d)
                    .expect("degree checked");
                Ok(vec![vec_poly::trim(a.coeffs()), vec_poly::trim(b.coeffs())])
            }
        }
    }

    /// Reconstruct `P_index` from its transformed components, checking the
    /// image conditions of the inverse exactly.
    pub fn inverse_row(
        &self,
        index: usize,
        components: &[Vec<Q>],
    ) -> Result<Vec<Q>, TransformError> {
        self.validate()?;
        if components.len() != self.component_count() {
            return Err(TransformError::ComponentCount {
                expected: self.component_count(),
                found: components.len(),
            });
        }
        let image = vec_poly::trim(&components[0]);
        match self {
            Self::MultiplyPower { power } => {
                divide_power(index, &image, *power).map_err(|_| TransformError::ImageCondition {
                    index,
                    reason: format!("x^{power} must divide Q_n"),
                })
            }
            Self::DividePower { power } => Ok(vec_poly::shift_by(&image, *power)),
            Self::Affine { .. } => {
                let (a, b) = self.affine_parameters()?.expect("affine");
                let inverse_a = a.recip();
                let inverse_b = -(&b * &inverse_a);
                Ok(vec_poly::compose(&image, &[inverse_b, inverse_a]))
            }
            Self::Reversal { .. } | Self::FPolynomial { .. } | Self::Mobius { .. } => {
                let m = self.mobius_parameters()?.expect("mobius");
                let d = m.schedule.degree_at(index)?;
                if vec_poly::degree(&image).is_some_and(|degree| degree > d) {
                    return Err(TransformError::ImageCondition {
                        index,
                        reason: format!("deg Q_n must be at most d_n = {d}"),
                    });
                }
                if let Self::Reversal { .. } = self {
                    return Ok(vec_poly::reverse(&image, d).expect("degree checked"));
                }
                let delta = m.delta();
                let unscaled =
                    mobius_substitute(&image, &m.e, &-m.b.clone(), &-m.c.clone(), &m.a, d);
                Ok(vec_poly::scale(
                    &unscaled,
                    &rational_power(&delta, -(d as i64)),
                ))
            }
            Self::Gamma { schedule } => {
                let d = schedule.degree_at(index)?;
                gamma_inverse(index, &image, d)
            }
            Self::SymmetricDecomposition { schedule } => {
                let d = schedule.degree_at(index)?;
                let b = vec_poly::trim(&components[1]);
                if !symmetric_at(&image, Some(d)) {
                    return Err(TransformError::ImageCondition {
                        index,
                        reason: format!("a_n must be symmetric at index d_n = {d}"),
                    });
                }
                if !symmetric_at(&b, d.checked_sub(1)) {
                    return Err(TransformError::ImageCondition {
                        index,
                        reason: format!(
                            "b_n must be symmetric at index d_n - 1 = {}",
                            d as i64 - 1
                        ),
                    });
                }
                Ok(vec_poly::add(&image, &vec_poly::shift(&b)))
            }
        }
    }

    /// Transport a rule `f P_m = sum s c D^k P_(m-r) + g` (target index `m`)
    /// into the new coordinates.
    ///
    /// The result is an exact consequence of the source rule at every index
    /// where the source rule holds and the transformed rows are defined.
    ///
    /// ```
    /// use polytool::recurrence::*;
    ///
    /// // Eulerian A008292: E_1 = 1, degree n - 1, and its f-polynomials.
    /// let (rule, _seeds, _first) = polytool::oeis::by_id("A008292")
    ///     .unwrap()
    ///     .recurrence_export_parts()
    ///     .unwrap();
    /// let f = RecurrenceTransformation::FPolynomial {
    ///     schedule: DegreeSchedule::new(1, 0, 1),
    /// };
    /// let TransportedRecurrence::Scalar(q) = f.transport(&rule).unwrap() else {
    ///     unreachable!()
    /// };
    /// assert_eq!(q.to_string(), "P(n) = (1 + 2t) P(n-1) + (t + t^2) P'(n-1)");
    /// ```
    pub fn transport(&self, rule: &Recurrence) -> Result<TransportedRecurrence, TransformError> {
        self.validate()?;
        match self {
            Self::MultiplyPower { power } => {
                let power =
                    i64::try_from(*power).map_err(|_| TransformError::Limit("power".into()))?;
                transport_scalar(rule, &Gauge::monomial(power)).map(TransportedRecurrence::Scalar)
            }
            Self::DividePower { power } => {
                let power =
                    i64::try_from(*power).map_err(|_| TransformError::Limit("power".into()))?;
                transport_scalar(rule, &Gauge::monomial(-power)).map(TransportedRecurrence::Scalar)
            }
            Self::Affine { .. } => {
                let (a, b) = self.affine_parameters()?.expect("affine");
                let gauge = Gauge::mobius(a, b, q(0), q(1), 0, 0);
                transport_scalar(rule, &gauge).map(TransportedRecurrence::Scalar)
            }
            Self::Reversal { .. } | Self::FPolynomial { .. } | Self::Mobius { .. } => {
                let m = self.mobius_parameters()?.expect("mobius");
                let gauge = Gauge::mobius(
                    m.a,
                    m.b,
                    m.c,
                    m.e,
                    m.schedule.weight_constant()?,
                    m.schedule.step,
                );
                transport_scalar(rule, &gauge).map(TransportedRecurrence::Scalar)
            }
            Self::Gamma { .. } => Err(TransformError::TransportUnsupported(
                "gamma coordinates are a quadratic substitution y = x/(1+x)^2 and x is not a \
                 polynomial in y; search the gamma rows directly"
                    .into(),
            )),
            Self::SymmetricDecomposition { schedule } => {
                transport_symmetric_decomposition(rule, schedule).map(TransportedRecurrence::Vector)
            }
        }
    }

    /// The exact row map, inverse, image conditions and conventions.
    pub fn record(&self) -> TransformationRecord {
        let scalar = || {
            vec![TransformComponent {
                name: "Q".into(),
                symmetry_index: None,
            }]
        };
        let schedule_text = self.schedule().map(DegreeSchedule::formula);
        let domain = self
            .schedule()
            .map_or_else(|| "all n >= 0".to_string(), DegreeSchedule::domain_text);
        let mobius_transport = "scalar: exact Möbius-gauge transport; a common power of ell is \
                                cleared, never evaluated at the pole";
        let (row_map, inverse, conditions, components, transport) = match self {
            Self::MultiplyPower { power } => (
                format!("Q_n(x) = x^{power} P_n(x)"),
                format!("P_n(x) = Q_n(x) / x^{power}"),
                vec![format!("inverse: x^{power} divides Q_n")],
                scalar(),
                "scalar: exact monomial-gauge transport".to_string(),
            ),
            Self::DividePower { power } => (
                format!("Q_n(x) = P_n(x) / x^{power}"),
                format!("P_n(x) = x^{power} Q_n(x)"),
                vec![format!("map: x^{power} divides P_n")],
                scalar(),
                "scalar: exact monomial-gauge transport (an identity of rational functions)"
                    .to_string(),
            ),
            Self::Reversal { .. } => (
                "Q_n(x) = x^(d_n) P_n(1/x)".to_string(),
                "P_n(x) = x^(d_n) Q_n(1/x)".to_string(),
                vec![
                    "map: deg P_n <= d_n".to_string(),
                    "inverse: deg Q_n <= d_n".to_string(),
                ],
                scalar(),
                mobius_transport.to_string(),
            ),
            Self::Affine { a, b } => (
                format!("Q_n(x) = P_n(({a}) x + ({b}))"),
                format!("P_n(y) = Q_n((y - ({b}))/({a}))"),
                Vec::new(),
                scalar(),
                "scalar: exact affine transport".to_string(),
            ),
            Self::Mobius { a, b, c, e, .. } => (
                format!(
                    "Q_n(x) = (({c}) x + ({e}))^(d_n) P_n((({a}) x + ({b}))/(({c}) x + ({e})))"
                ),
                format!(
                    "P_n(y) = (({a})({e}) - ({b})({c}))^(-d_n) (({a}) - ({c}) y)^(d_n) \
                     Q_n((({e}) y - ({b}))/(({a}) - ({c}) y))"
                ),
                vec![
                    "map: deg P_n <= d_n".to_string(),
                    "inverse: deg Q_n <= d_n".to_string(),
                ],
                scalar(),
                mobius_transport.to_string(),
            ),
            Self::FPolynomial { .. } => (
                "Q_n(x) = (1 + x)^(d_n) P_n(x/(1 + x))  (decomposition::f_polynomial at bound d_n)"
                    .to_string(),
                "P_n(y) = (1 - y)^(d_n) Q_n(y/(1 - y))".to_string(),
                vec![
                    "map: deg P_n <= d_n".to_string(),
                    "inverse: deg Q_n <= d_n".to_string(),
                ],
                scalar(),
                mobius_transport.to_string(),
            ),
            Self::Gamma { .. } => (
                "P_n(x) = sum_i Q_(n,i) x^i (1 + x)^(d_n - 2i)".to_string(),
                "P_n(x) = sum_i Q_(n,i) x^i (1 + x)^(d_n - 2i)".to_string(),
                vec![
                    "map: x^(d_n) P_n(1/x) = P_n (symmetry at index d_n)".to_string(),
                    "inverse: deg Q_n <= floor(d_n/2)".to_string(),
                ],
                scalar(),
                "unsupported: search the gamma rows directly".to_string(),
            ),
            Self::SymmetricDecomposition { .. } => (
                "P_n(x) = a_n(x) + x b_n(x), (a_n, b_n) = ((P_n - x Pbar_n)/(1 - x), \
                 (Pbar_n - P_n)/(1 - x)) with Pbar_n = x^(d_n) P_n(1/x)"
                    .to_string(),
                "P_n(x) = a_n(x) + x b_n(x)".to_string(),
                vec![
                    "map: deg P_n <= d_n".to_string(),
                    "inverse: x^(d_n) a_n(1/x) = a_n and x^(d_n - 1) b_n(1/x) = b_n".to_string(),
                ],
                vec![
                    TransformComponent {
                        name: "a".into(),
                        symmetry_index: Some("d_n".into()),
                    },
                    TransformComponent {
                        name: "b".into(),
                        symmetry_index: Some("d_n - 1".into()),
                    },
                ],
                "vector: first-order rules without alternating signs become the complete \
                 two-component system in (a_n, b_n); no component is selected"
                    .to_string(),
            ),
        };
        TransformationRecord {
            kind: self.kind_name().to_string(),
            row_map,
            inverse,
            image_conditions: conditions,
            degree_schedule: schedule_text,
            index_domain: domain,
            components,
            transport,
        }
    }
}

// ---------------------------------------------------------------------------
// Row maps
// ---------------------------------------------------------------------------

fn check_degree(index: usize, row: &[Q], bound: usize) -> Result<(), TransformError> {
    match vec_poly::degree(row) {
        Some(degree) if degree > bound => Err(TransformError::DegreeExceedsBound {
            index,
            degree,
            bound,
        }),
        _ => Ok(()),
    }
}

fn divide_power(index: usize, row: &[Q], power: usize) -> Result<Vec<Q>, TransformError> {
    if vec_poly::is_zero(row) {
        return Ok(vec![Q::zero()]);
    }
    if row.iter().take(power).any(|c| !c.is_zero()) {
        return Err(TransformError::NotDivisible { index, power });
    }
    Ok(vec_poly::trim(&row[power.min(row.len())..]))
}

/// `sum_k p_k num^k den^(total - k)` for `deg p <= total`, by a homogeneous
/// Horner scheme (quadratically many coefficient operations for linear maps).
fn homogeneous_substitute_ring<C: crate::polynomial::CoeffRing>(
    p: &[C],
    numerator: &[C],
    denominator: &[C],
    total: usize,
) -> Vec<C> {
    let Some(degree) = vec_poly::degree(p) else {
        return vec![C::zero()];
    };
    let mut result = vec![p[degree].clone()];
    let mut denominator_power = vec![C::one()];
    for k in (0..degree).rev() {
        denominator_power = vec_poly::mul(&denominator_power, denominator);
        result = vec_poly::mul(&result, numerator);
        if !p[k].is_zero() {
            result = vec_poly::add(&result, &vec_poly::scale(&denominator_power, &p[k]));
        }
    }
    vec_poly::mul(&result, &vec_poly::pow(denominator, total - degree))
}

/// Integer coefficients of `p`, when every coefficient is an integer.
fn integer_coefficients(p: &[Q]) -> Option<Vec<BigInt>> {
    p.iter()
        .map(|c| c.is_integer().then(|| c.to_integer()))
        .collect()
}

/// [`homogeneous_substitute_ring`] over `BigInt` when all inputs are integral
/// (avoiding rational normalization), otherwise over the rationals.
fn homogeneous_substitute(p: &[Q], numerator: &[Q], denominator: &[Q], total: usize) -> Vec<Q> {
    match (
        integer_coefficients(p),
        integer_coefficients(numerator),
        integer_coefficients(denominator),
    ) {
        (Some(p), Some(numerator), Some(denominator)) => {
            homogeneous_substitute_ring(&p, &numerator, &denominator, total)
                .into_iter()
                .map(Q::from_integer)
                .collect()
        }
        _ => homogeneous_substitute_ring(p, numerator, denominator, total),
    }
}

/// `sum_k p_k (a x + b)^k (c x + e)^(d - k)` for `deg p <= d`.
fn mobius_substitute(row: &[Q], a: &Q, b: &Q, c: &Q, e: &Q, d: usize) -> Vec<Q> {
    homogeneous_substitute(row, &[b.clone(), a.clone()], &[e.clone(), c.clone()], d)
}

fn rational_power(base: &Q, exponent: i64) -> Q {
    let mut result = Q::one();
    for _ in 0..exponent.unsigned_abs() {
        result *= base;
    }
    if exponent < 0 {
        result.recip()
    } else {
        result
    }
}

/// Whether `p` is symmetric at `index` (`None` means index `-1`: only zero).
fn symmetric_at(p: &[Q], index: Option<usize>) -> bool {
    match index {
        None => vec_poly::is_zero(p),
        Some(d) => vec_poly::reverse(p, d).is_some_and(|r| r == vec_poly::trim(p)),
    }
}

/// Gamma coefficients at symmetry index `d`, adapting
/// [`gamma_coefficients_bigint_coeffs`] (which uses the actual degree) to rows
/// with endpoint factors `x^v` and rational coefficients.
fn gamma_at_index(index: usize, row: &[Q], d: usize) -> Result<Vec<Q>, TransformError> {
    check_degree(index, row, d)?;
    if !symmetric_at(row, Some(d)) {
        return Err(TransformError::NotSymmetric {
            index,
            symmetry_index: d,
        });
    }
    let Some(degree) = vec_poly::degree(row) else {
        return Ok(vec![Q::zero()]);
    };
    let valuation = row.iter().position(|c| !c.is_zero()).expect("nonzero row");
    let common = row[valuation..=degree]
        .iter()
        .fold(BigInt::one(), |acc, c| acc.lcm(c.denom()));
    let integers = row[valuation..=degree]
        .iter()
        .map(|c| (c * Q::from_integer(common.clone())).to_integer())
        .collect::<Vec<_>>();
    let gamma =
        gamma_coefficients_bigint_coeffs(&integers).ok_or(TransformError::NotSymmetric {
            index,
            symmetry_index: d,
        })?;
    let mut result = vec![Q::zero(); valuation];
    result.extend(gamma.into_iter().map(|g| Q::new(g, common.clone())));
    Ok(vec_poly::trim(&result))
}

fn gamma_inverse(index: usize, gamma: &[Q], d: usize) -> Result<Vec<Q>, TransformError> {
    if vec_poly::degree(gamma).is_some_and(|degree| degree > d / 2) {
        return Err(TransformError::ImageCondition {
            index,
            reason: format!("deg Q_n must be at most floor(d_n/2) = {}", d / 2),
        });
    }
    // sum_i g_i x^i (1+x)^(d-2i) = (1+x)^(d mod 2) sum_i g_i x^i ((1+x)^2)^(d/2 - i).
    let one_plus_x = [Q::one(), Q::one()];
    let even = homogeneous_substitute(
        gamma,
        &[Q::zero(), Q::one()],
        &vec_poly::pow(&one_plus_x, 2),
        d / 2,
    );
    Ok(vec_poly::mul(&even, &vec_poly::pow(&one_plus_x, d % 2)))
}

// ---------------------------------------------------------------------------
// Bivariate helpers: `[n power][x power]`, zero is the empty vector
// ---------------------------------------------------------------------------

type Bp = Vec<Vec<Q>>;

fn bp_trim(mut p: Bp) -> Bp {
    for row in &mut p {
        vec_poly::strip_trailing_zeros(row);
    }
    while p.last().is_some_and(Vec::is_empty) {
        p.pop();
    }
    p
}

fn bp_from(b: &BivarPoly) -> Bp {
    bp_trim(b.coeffs.clone())
}

fn bp_to_bivar(p: Bp) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(bp_trim(p))
}

fn bp_is_zero(p: &Bp) -> bool {
    p.iter().all(|row| vec_poly::is_zero(row))
}

fn bp_constant(c: Q) -> Bp {
    bp_trim(vec![vec![c]])
}

fn bp_add_assign(acc: &mut Bp, p: &Bp) {
    if acc.len() < p.len() {
        acc.resize(p.len(), Vec::new());
    }
    for (target, row) in acc.iter_mut().zip(p) {
        *target = vec_poly::add(target, row);
    }
}

fn bp_add(a: &Bp, b: &Bp) -> Bp {
    let mut sum = a.clone();
    bp_add_assign(&mut sum, b);
    bp_trim(sum)
}

fn bp_scale(p: &Bp, s: &Q) -> Bp {
    bp_trim(p.iter().map(|row| vec_poly::scale(row, s)).collect())
}

fn bp_mul(a: &Bp, b: &Bp) -> Bp {
    if bp_is_zero(a) || bp_is_zero(b) {
        return Vec::new();
    }
    let mut result: Bp = vec![Vec::new(); a.len() + b.len() - 1];
    for (i, left) in a.iter().enumerate() {
        for (j, right) in b.iter().enumerate() {
            result[i + j] = vec_poly::add(&result[i + j], &vec_poly::mul(left, right));
        }
    }
    bp_trim(result)
}

fn bp_mul_x(p: &Bp, poly: &[Q]) -> Bp {
    bp_trim(p.iter().map(|row| vec_poly::mul(row, poly)).collect())
}

fn bp_dx(p: &Bp) -> Bp {
    bp_trim(p.iter().map(|row| vec_poly::derivative(row)).collect())
}

/// `c0 + c1 n`.
fn bp_affine_index(c0: Q, c1: Q) -> Bp {
    bp_trim(vec![vec![c0], vec![c1]])
}

fn bp_x_degree(p: &Bp) -> Option<usize> {
    p.iter().filter_map(|row| vec_poly::degree(row)).max()
}

fn bp_div_linear(p: &Bp, lambda: &[Q]) -> Option<Bp> {
    let root = -(&lambda[0] / &lambda[1]);
    let inverse_lead = lambda[1].recip();
    p.iter()
        .map(|row| {
            vec_poly::exact_div_by_t_minus(row, &root)
                .map(|quotient| vec_poly::scale(&quotient, &inverse_lead))
        })
        .collect::<Option<Bp>>()
        .map(bp_trim)
}

fn bp_is_constant(p: &Bp) -> Option<Q> {
    if bp_is_zero(p) {
        return None;
    }
    let tail_zero = p.iter().enumerate().all(|(i, row)| {
        row.iter()
            .enumerate()
            .all(|(j, c)| (i == 0 && j == 0) || c.is_zero())
    });
    tail_zero.then(|| p[0][0].clone())
}

/// `sum c_(i,j) n^i num^j den^(D - j)` with `D` the `x`-degree of `c`.
fn bp_homogenize(p: &Bp, numerator: &[Q], denominator: &[Q]) -> (Bp, usize) {
    let degree = bp_x_degree(p).unwrap_or(0);
    let numerator_powers = (0..=degree)
        .map(|j| vec_poly::pow(numerator, j))
        .collect::<Vec<_>>();
    let denominator_powers = (0..=degree)
        .map(|j| vec_poly::pow(denominator, j))
        .collect::<Vec<_>>();
    let result = p
        .iter()
        .map(|row| {
            let mut sum = vec![Q::zero()];
            for (j, coefficient) in row.iter().enumerate() {
                if coefficient.is_zero() {
                    continue;
                }
                let basis = vec_poly::mul(&numerator_powers[j], &denominator_powers[degree - j]);
                sum = vec_poly::add(&sum, &vec_poly::scale(&basis, coefficient));
            }
            sum
        })
        .collect();
    (bp_trim(result), degree)
}

// ---------------------------------------------------------------------------
// Scalar transport
// ---------------------------------------------------------------------------

/// A scalar gauge `Q_n(x) = lambda(x)^(w + h n) P_n(phi(x))` with
/// `phi = num/den` and `P^(k)(phi) = kappa^k lambda^(alpha_k) T_k[Q]`,
/// `alpha_0 = -(w + h n)`, `alpha_(k+1) = alpha_k + sigma`,
/// `T_(k+1) = lambda T_k' + alpha_k lambda' T_k`.
struct Gauge {
    numerator: Vec<Q>,
    denominator: Vec<Q>,
    /// Exponent of `lambda` contributed by `den^D`: 1 for Möbius (`den = lambda`), 0 otherwise.
    tau: i64,
    lambda: [Q; 2],
    weight_constant: i64,
    step: i64,
    sigma: i64,
    kappa: Q,
}

impl Gauge {
    fn monomial(power: i64) -> Self {
        Self {
            numerator: vec![q(0), q(1)],
            denominator: vec![q(1)],
            tau: 0,
            lambda: [q(0), q(1)],
            weight_constant: power,
            step: 0,
            sigma: -1,
            kappa: q(1),
        }
    }

    fn mobius(a: Q, b: Q, c: Q, e: Q, weight_constant: i64, step: i64) -> Self {
        let delta = &a * &e - &b * &c;
        Self {
            numerator: vec_poly::trim(&[b, a]),
            denominator: vec_poly::trim(&[e.clone(), c.clone()]),
            tau: 1,
            lambda: [e, c],
            weight_constant,
            step,
            sigma: 1,
            kappa: delta.recip(),
        }
    }

    fn lambda_is_constant(&self) -> bool {
        self.lambda[1].is_zero()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum PieceKey {
    Lhs,
    Term {
        offset: usize,
        derivative: usize,
        alternating: bool,
    },
    Inhomogeneous,
}

fn checked(value: Option<i64>) -> Result<i64, TransformError> {
    value.ok_or_else(|| TransformError::Limit("exponent overflow".into()))
}

fn transport_scalar(rule: &Recurrence, gauge: &Gauge) -> Result<Recurrence, TransformError> {
    let mut pieces: Vec<(PieceKey, Bp, i64)> = Vec::new();
    let tau = gauge.tau;
    let degree_exponent = |degree: usize| -> Result<i64, TransformError> {
        checked(i64::try_from(degree).ok().and_then(|d| d.checked_mul(tau)))
    };

    let lhs = rule
        .denominator
        .as_ref()
        .map(bp_from)
        .unwrap_or_else(|| bp_constant(Q::one()));
    let (lhs_hom, lhs_degree) = bp_homogenize(&lhs, &gauge.numerator, &gauge.denominator);
    pieces.push((PieceKey::Lhs, lhs_hom, -degree_exponent(lhs_degree)?));

    let lambda_derivative = gauge.lambda[1].clone();
    let lambda_poly = vec_poly::trim(&gauge.lambda);
    for term in &rule.terms {
        let coefficient = bp_from(&term.coeff);
        if bp_is_zero(&coefficient) {
            continue;
        }
        let (hom, degree) = bp_homogenize(&coefficient, &gauge.numerator, &gauge.denominator);
        let offset =
            i64::try_from(term.offset).map_err(|_| TransformError::Limit("offset".into()))?;
        let derivative =
            i64::try_from(term.deriv_order).map_err(|_| TransformError::Limit("order".into()))?;
        // alpha_0(m) = -(w + h (m - r)) = (h r - w) - h m.
        let alpha_constant = checked(
            gauge
                .step
                .checked_mul(offset)
                .and_then(|value| value.checked_sub(gauge.weight_constant)),
        )?;
        let mut alpha = bp_affine_index(q(alpha_constant), q(-gauge.step));
        let mut operator: Vec<Bp> = vec![bp_constant(Q::one())];
        for _ in 0..term.deriv_order {
            let mut next: Vec<Bp> = vec![Vec::new(); operator.len() + 1];
            for (i, t) in operator.iter().enumerate() {
                let from_derivative = bp_mul_x(&bp_dx(t), &lambda_poly);
                let from_order = bp_mul_x(t, &lambda_poly);
                let from_weight = bp_scale(&bp_mul(&alpha, t), &lambda_derivative);
                next[i] = bp_add(&bp_add(&next[i], &from_derivative), &from_weight);
                next[i + 1] = bp_add(&next[i + 1], &from_order);
            }
            operator = next;
            alpha = bp_add(&alpha, &bp_constant(q(gauge.sigma)));
        }
        let scale = rational_power(&gauge.kappa, derivative);
        // ell^(w(m) - w(m - r)) times the chain-rule exponent sigma k.
        let exponent = checked(
            gauge
                .step
                .checked_mul(offset)
                .and_then(|value| value.checked_add(gauge.sigma * derivative))
                .and_then(|value| value.checked_sub(degree_exponent(degree).ok()?)),
        )?;
        for (order, t) in operator.into_iter().enumerate() {
            let poly = bp_scale(&bp_mul(&hom, &t), &scale);
            if bp_is_zero(&poly) {
                continue;
            }
            pieces.push((
                PieceKey::Term {
                    offset: term.offset,
                    derivative: order,
                    alternating: term.sign == RecurrenceSign::AlternatingN,
                },
                poly,
                exponent,
            ));
        }
    }

    if let Some(inhomogeneous) = &rule.inhomogeneous {
        let g = bp_from(inhomogeneous);
        if !bp_is_zero(&g) {
            let index_independent = gauge.step == 0
                || (gauge.lambda_is_constant()
                    && rational_power(&gauge.lambda[0], gauge.step).is_one());
            if !index_independent {
                return Err(TransformError::InhomogeneousNotTransportable);
            }
            let (hom, degree) = bp_homogenize(&g, &gauge.numerator, &gauge.denominator);
            let exponent = checked(gauge.weight_constant.checked_sub(degree_exponent(degree)?))?;
            pieces.push((PieceKey::Inhomogeneous, hom, exponent));
        }
    }

    // Combine: constant lambda gives scalar factors; otherwise clear denominators.
    let mut combined: BTreeMap<PieceKey, Bp> = BTreeMap::new();
    if gauge.lambda_is_constant() {
        for (key, poly, exponent) in pieces {
            let scaled = bp_scale(&poly, &rational_power(&gauge.lambda[0], exponent));
            let entry = combined.entry(key).or_default();
            *entry = bp_add(entry, &scaled);
        }
    } else {
        let shift = pieces.iter().map(|(_, _, e)| -e).max().unwrap_or(0);
        for (key, poly, exponent) in pieces {
            let power = usize::try_from(shift + exponent).expect("nonnegative by choice of shift");
            let scaled = bp_mul_x(&poly, &vec_poly::pow(&lambda_poly, power));
            let entry = combined.entry(key).or_default();
            *entry = bp_add(entry, &scaled);
        }
        // Divide out the largest common power of lambda (exact polynomial division).
        loop {
            let divided = combined
                .iter()
                .map(|(key, poly)| bp_div_linear(poly, &gauge.lambda).map(|p| (key.clone(), p)))
                .collect::<Option<BTreeMap<_, _>>>();
            match divided {
                Some(next) if next.get(&PieceKey::Lhs).is_some_and(|p| !bp_is_zero(p)) => {
                    combined = next;
                }
                _ => break,
            }
        }
    }
    Ok(assemble_rule(combined))
}

fn assemble_rule(mut combined: BTreeMap<PieceKey, Bp>) -> Recurrence {
    let mut lhs = combined.remove(&PieceKey::Lhs).unwrap_or_default();
    if let Some(constant) = bp_is_constant(&lhs) {
        let inverse = constant.recip();
        for poly in combined.values_mut() {
            *poly = bp_scale(poly, &inverse);
        }
        lhs = bp_constant(Q::one());
    }
    let inhomogeneous = combined
        .remove(&PieceKey::Inhomogeneous)
        .filter(|g| !bp_is_zero(g))
        .map(bp_to_bivar);
    let terms = combined
        .into_iter()
        .filter(|(_, poly)| !bp_is_zero(poly))
        .map(|(key, poly)| match key {
            PieceKey::Term {
                offset,
                derivative,
                alternating,
            } => RecurrenceTerm::new(offset, derivative, bp_to_bivar(poly)).with_sign(
                if alternating {
                    RecurrenceSign::AlternatingN
                } else {
                    RecurrenceSign::None
                },
            ),
            _ => unreachable!("only term pieces remain"),
        })
        .collect();
    Recurrence {
        terms,
        denominator: (bp_is_constant(&lhs) != Some(Q::one())).then(|| bp_to_bivar(lhs)),
        inhomogeneous,
    }
}

fn canonical_pieces(rule: &Recurrence) -> BTreeMap<PieceKey, Bp> {
    let mut combined: BTreeMap<PieceKey, Bp> = BTreeMap::new();
    combined.insert(
        PieceKey::Lhs,
        rule.denominator
            .as_ref()
            .map(bp_from)
            .unwrap_or_else(|| bp_constant(Q::one())),
    );
    for term in &rule.terms {
        let key = PieceKey::Term {
            offset: term.offset,
            derivative: term.deriv_order,
            alternating: term.sign == RecurrenceSign::AlternatingN,
        };
        let entry = combined.entry(key).or_default();
        *entry = bp_add(entry, &bp_from(&term.coeff));
    }
    if let Some(g) = &rule.inhomogeneous {
        combined.insert(PieceKey::Inhomogeneous, bp_from(g));
    }
    if let Some(constant) = combined.get(&PieceKey::Lhs).and_then(bp_is_constant) {
        let inverse = constant.recip();
        for poly in combined.values_mut() {
            *poly = bp_scale(poly, &inverse);
        }
    }
    combined.retain(|key, poly| *key == PieceKey::Lhs || !bp_is_zero(poly));
    combined
}

/// Canonical form of a rule: terms with equal offset, derivative order and
/// sign merged, zero terms dropped, coefficients trimmed and a constant left
/// factor divided out.
pub fn canonical_recurrence(rule: &Recurrence) -> Recurrence {
    assemble_rule(canonical_pieces(rule))
}

/// Whether two rules have the same canonical form (identical recursions).
pub fn same_canonical_recurrence(left: &Recurrence, right: &Recurrence) -> bool {
    canonical_pieces(left) == canonical_pieces(right)
}

// ---------------------------------------------------------------------------
// Symmetric decomposition transport
// ---------------------------------------------------------------------------

fn transport_symmetric_decomposition(
    rule: &Recurrence,
    schedule: &DegreeSchedule,
) -> Result<SharedDenominatorVectorRecurrence, TransformError> {
    if rule.terms.iter().any(|t| t.offset != 1) {
        return Err(TransformError::TransportUnsupported(
            "the symmetric decomposition transports first-order rules (all offsets 1)".into(),
        ));
    }
    if rule.terms.iter().any(RecurrenceTerm::is_alternating) {
        return Err(TransformError::TransportUnsupported(
            "vector systems carry no (-1)^n factors".into(),
        ));
    }
    let reversal = RecurrenceTransformation::Reversal {
        schedule: schedule.clone(),
    };
    let reversed = match reversal.transport(rule)? {
        TransportedRecurrence::Scalar(r) => r,
        TransportedRecurrence::Vector(_) => unreachable!("reversal is scalar"),
    };
    let source = canonical_recurrence(rule);
    let one = bp_constant(Q::one());
    let f = source
        .denominator
        .as_ref()
        .map(bp_from)
        .unwrap_or_else(|| one.clone());
    let f_bar = reversed
        .denominator
        .as_ref()
        .map(bp_from)
        .unwrap_or_else(|| one.clone());
    let one_minus_x = [q(1), q(-1)];
    let x = [q(0), q(1)];
    let (denominator, u, u_bar) = if f == f_bar {
        (bp_mul_x(&f, &one_minus_x), one.clone(), one)
    } else {
        (bp_mul_x(&bp_mul(&f, &f_bar), &one_minus_x), f_bar, f)
    };
    let coefficients = |r: &Recurrence| -> Vec<Bp> {
        let order = r.terms.iter().map(|t| t.deriv_order).max().unwrap_or(0);
        let mut c = vec![Vec::new(); order + 2];
        for term in &r.terms {
            c[term.deriv_order] = bp_add(&c[term.deriv_order], &bp_from(&term.coeff));
        }
        c
    };
    let mut c = coefficients(&source);
    let mut c_bar = coefficients(&reversed);
    let width = c.len().max(c_bar.len());
    c.resize(width, Vec::new());
    c_bar.resize(width, Vec::new());
    let order_count = width - 1;
    let mut matrix = vec![vec![vec![Vec::new(); order_count]; 2]; 2];
    for k in 0..order_count {
        let uc = bp_mul(&u, &c[k]);
        let ubc = bp_mul(&u_bar, &c_bar[k]);
        let lowered = bp_scale(&bp_mul(&u, &c[k + 1]), &q((k + 1) as i64));
        let x_uc = bp_mul_x(&uc, &x);
        let x_ubc = bp_mul_x(&ubc, &x);
        matrix[0][0][k] = bp_add(&uc, &bp_scale(&x_ubc, &q(-1)));
        matrix[0][1][k] = bp_add(&bp_add(&x_uc, &bp_scale(&x_ubc, &q(-1))), &lowered);
        matrix[1][0][k] = bp_add(&ubc, &bp_scale(&uc, &q(-1)));
        matrix[1][1][k] = bp_add(
            &bp_add(&ubc, &bp_scale(&x_uc, &q(-1))),
            &bp_scale(&lowered, &q(-1)),
        );
    }
    let g = source
        .inhomogeneous
        .as_ref()
        .map(bp_from)
        .unwrap_or_default();
    let g_bar = reversed
        .inhomogeneous
        .as_ref()
        .map(bp_from)
        .unwrap_or_default();
    let forcing = (!bp_is_zero(&g) || !bp_is_zero(&g_bar)).then(|| {
        let ug = bp_mul(&u, &g);
        let ubg = bp_mul(&u_bar, &g_bar);
        vec![
            bp_add(&ug, &bp_scale(&bp_mul_x(&ubg, &x), &q(-1))),
            bp_add(&ubg, &bp_scale(&ug, &q(-1))),
        ]
    });
    // Coefficients are in the target index m = n + 1; the vector convention uses n.
    let to_source_index = |p: Bp| bp_to_bivar(p).shift_index(1);
    Ok(SharedDenominatorVectorRecurrence {
        denominator: to_source_index(denominator),
        matrix: matrix
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|entry| {
                        WeylOperator::new(entry.into_iter().map(to_source_index).collect())
                    })
                    .collect()
            })
            .collect(),
        forcing: forcing.map(|components| components.into_iter().map(to_source_index).collect()),
    })
}

// ---------------------------------------------------------------------------
// Residual checks
// ---------------------------------------------------------------------------

/// Whether the scalar rule holds exactly at target index `target`, where
/// `rows[i]` is the row with original index `origin + i`. `None` when a
/// referenced row is unavailable.
pub fn transformed_rule_holds_at(
    rule: &Recurrence,
    rows: &[Vec<Q>],
    origin: usize,
    target: usize,
) -> Option<bool> {
    let position = target.checked_sub(origin)?;
    let current = rows.get(position)?;
    let lhs = match &rule.denominator {
        Some(f) => poly_mul_rational(&bivar_eval_n(f, target), current),
        None => current.clone(),
    };
    let mut rhs = rational_zero_poly();
    for term in &rule.terms {
        let source = rows.get(position.checked_sub(term.offset)?)?;
        let derivative = poly_nth_derivative_rational(source, term.deriv_order);
        let product = poly_mul_rational(&bivar_eval_n(&term.coeff, target), &derivative);
        let sign = if term.is_alternating() && target % 2 == 1 {
            q(-1)
        } else {
            Q::one()
        };
        poly_add_scaled_assign(&mut rhs, &product, &sign);
    }
    if let Some(g) = &rule.inhomogeneous {
        poly_add_scaled_assign(&mut rhs, &bivar_eval_n(g, target), &Q::one());
    }
    Some(vec_poly::trim(&lhs) == vec_poly::trim(&rhs))
}

#[path = "transform_report.rs"]
mod report;
pub use report::*;

#[cfg(test)]
#[path = "transform_tests.rs"]
mod tests;
