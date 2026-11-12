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
//! is transported by [`RecurrenceTransformation::transport`].  The
//! coordinate maps (`multiply_power`, `divide_power`, `affine`, `reversal`,
//! `f_polynomial`, `mobius`) are substitutions `Q_n = C_n(P_n)` of the
//! `ore-algebra` crate (`ore_algebra::Substitution`: power gauges `x^(+-r)`,
//! the affine map, the reversal, the `f`-polynomial and the Möbius map with
//! scale `1`, with the schedule `d_0 + h (n - n_0)` as the exponent
//! `[d_0 - h n_0, h]` in the original index).  Write the rule as the
//! relation `L = f(n + R) S^R - sum s c(n + R) D^k S^(R - r)` at base
//! `n = m - R` (`R` the largest offset), so `(L P)_n = g(n + R)`.  Its
//! forward rule (`ore_algebra::push_forward_substitution`, the pullback
//! through the inverse substitution) is `q C L C^(-1)` with a left factor
//! `q` that is a nonzero constant times a power of the linear form
//! `ell` of the map (`ell = c x + e` for a Möbius map, `x` for a power
//! gauge, `1` for an affine map).  For fixed bounds the pullback is linear
//! and acts coefficient by coefficient without shifting the index of a
//! coefficient, so it commutes with a factor `(-1)^m` of a term: each piece
//! of the rule (the left factor and every term, with its sign) is pushed
//! forward separately with the bounds of the whole rule, which keeps the
//! left factor `q` common.  The additive term becomes
//! `q C_n(g(n + R)) = q ell^(d_n - deg g) hom(g)` with
//! `hom(g) = ell^(deg g) g(phi)`, `phi = (a x + b)/ell` (for a power gauge
//! `q x^(+-r) g`).  A power of `ell` that clears every denominator is
//! multiplied in, the largest common power of `ell` is divided out exactly,
//! and the result is scaled so that its left factor is `ell^e hom(f)` (the
//! forward rule fixes it only up to a constant).  Nothing is evaluated at
//! the pole of `phi`.  For the first-order rule
//! `R(n) P_(n+1) = A P_n + B P_n'` and a Möbius map with
//! `Delta = a e - b c` this is the formula
//! `R Q_(n+1) = ell^h (A(phi) - d_n c ell B(phi)/Delta) Q_n
//!  + ell^(h+2) B(phi)/Delta Q_n'`.
//!
//! The transported rule is an exact algebraic consequence of the source rule,
//! valid wherever the source rule holds; it is not refitted.  An additive
//! term `g` picks up the factor `ell^(d_n)`, so it is transported only when
//! that factor is index independent (`h = 0`, or a constant `ell` with
//! `ell^h = 1`).  The symmetric decomposition transports a first-order rule
//! without alternating signs to the complete two-component
//! [`SharedDenominatorVectorRecurrence`]; no component is selected.
//!
//! Gamma coordinates `P_n = (1 + x)^(d_n) G_n(x/(1 + x)^2)` are a quadratic
//! substitution, so they are transported differently ([`transport_gamma`]):
//! derivatives of `(1 + x)^d G(z)` are expanded by the chain rule, a power of
//! `1 + x` is cleared, and the identity is split into its symmetric and
//! antisymmetric parts, which are two relations among the `G^(i)_(m-r)`
//! because `Q(x)` is the quadratic extension of `Q(z)` by
//! `w = (1 - x)/(1 + x)`, `w^2 = 1 - 4z`.  When the parts are proportional
//! (for example when the rule is invariant under reversal at `d_n`) the part
//! with a nonzero leader is the transported scalar rule; otherwise
//! [`RecurrenceTransformation::transport`] reports transport as unsupported
//! with both parts, and [`transform_recurrence`] states the rank-two case as
//! the two-relation system of the even part and its odd companion
//! ([`GammaTransport::relation_system`], report field `transported.system`).  Canonical Eulerian A008292
//! (`E_1 = 1`, `d_n = n - 1`) gives
//! `G_m = (1 + 2(m - 2) z) G_(m-1) + z (1 - 4z) G'_(m-1)`.  The gamma rule
//! holds where the rows involved are symmetric, which is the map's image
//! condition.
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
//! (see [`RootDomainReport`]). Factorial rescalings are a separate extension.
//! The diagonal maps `T_n(x^k) = lambda(n,k) x^k` and derivative antecedents
//! `P_n = D^m(g Q_n)/m!` are registered here but transport relations as
//! implicit relations (module `implicit::diagonal`, reports in `implicit::report`).
//! The lossy map `derivative` (`Q_n = D^m P_n`) is registered as well; its
//! lost coefficients, reconstruction and transport live in module `implicit::lossy`.
//!
//! Scalar bijective maps also act componentwise on vector states, with
//! schedules shifted by each component's degree offset, and first-order
//! vector systems are transported exactly through Möbius-class maps
//! ([`RecurrenceTransformation::map_state`],
//! [`transport_vector_componentwise`]; module `transform::vector_transport`).

use super::*;
use crate::decomposition::f_polynomial;
use crate::polynomial::Polynomial;
use crate::real_rootedness::gamma_coefficients_bigint_coeffs;
use crate::recurrence::implicit::algebra::{Biv, DOp};
use crate::vec_poly;
use num_traits::{One, Zero};
use ore_algebra::{
    pull_back_substitution, MobiusMap, PullbackBounds, Substitution, SubstitutionError,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt;

pub(super) type Q = BigRational;

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
    pub(super) fn weight_constant(&self) -> Result<i64, TransformError> {
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

/// The errors of an `ore-algebra` substitution as registry errors: a zero
/// affine slope (or Möbius scale) and a singular Möbius map are
/// [`TransformError::SingularMobius`], an exponent that is not
/// `d_0 + h n` with integers is [`TransformError::InvalidSchedule`].
pub(super) fn substitution_error(error: SubstitutionError) -> TransformError {
    match error {
        SubstitutionError::ZeroParameter | SubstitutionError::SingularMobius => {
            TransformError::SingularMobius
        }
        SubstitutionError::InvalidExponent => TransformError::InvalidSchedule(error.to_string()),
    }
}

// ---------------------------------------------------------------------------
// The transformation registry
// ---------------------------------------------------------------------------

/// An exact change of polynomial coordinates `P_n -> Q_n`.
///
/// Rational parameters are integer or `num/den` strings.  The JSON form is
/// tagged by `kind` (`multiply_power`, `divide_power`, `reversal`, `affine`,
/// `mobius`, `f_polynomial`, `gamma`, `symmetric_decomposition`, `diagonal`,
/// `derivative_antecedent`, `derivative`).
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
    /// `Q_n = T_n P_n` with `T_n(x^k) = lambda(n,k) x^k`; requires every
    /// needed weight to be defined and nonzero. Relations are transported as
    /// implicit relations ([`RecurrenceTransformation::transport_implicit`]).
    Diagonal {
        /// The weight `lambda(n,k)`.
        weight: DiagonalWeight,
    },
    /// `Q_n = T_g^-1 P_n`, the antecedent with `P_n = D^m(g Q_n)/m!` for a fixed
    /// nonzero rational `g` of degree `m`; bijective (triangular with nonzero
    /// diagonal `lc(g) binom(k+m,m)`). Relations are transported by
    /// substitution ([`RecurrenceTransformation::transport_implicit`]).
    DerivativeAntecedent {
        /// Ascending rational coefficients of `g`, as strings.
        g: Vec<String>,
    },
    /// `Q_n = D^order P_n`: plain differentiation, a lossy map. The
    /// coefficients `[x^0..x^(order-1)] P_n` are lost; reconstruction needs
    /// them as recorded boundary data
    /// ([`RecurrenceTransformation::inverse_row_with_lost_data`]), and
    /// [`RecurrenceTransformation::inverse_row`] rejects the map. Relations
    /// are transported as implicit relations
    /// ([`RecurrenceTransformation::transport_implicit`]).
    Derivative {
        /// The derivative order `m >= 1`.
        order: usize,
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

pub(super) fn q(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

/// Möbius parameters `(a, b, c, e)` and schedule of a Möbius-type map.
pub(super) struct MobiusParameters<'a> {
    pub(super) a: Q,
    pub(super) b: Q,
    pub(super) c: Q,
    pub(super) e: Q,
    pub(super) schedule: &'a DegreeSchedule,
}

impl MobiusParameters<'_> {
    fn delta(&self) -> Q {
        &self.a * &self.e - &self.b * &self.c
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

pub(super) fn rational_power(base: &Q, exponent: i64) -> Q {
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
    let (common, integers) =
        crate::vec_poly::clear_denominators_with_scale(&row[valuation..=degree]);
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

#[cfg(test)]
mod tests;

mod polynomial;
mod row_maps;
use polynomial::*;
pub(super) use polynomial::{
    bp_add, bp_affine_index, bp_constant, bp_div_linear, bp_dx, bp_from, bp_homogenize, bp_is_zero,
    bp_mul, bp_mul_x, bp_scale, bp_to_bivar, bp_trim, bp_x_degree, checked, Bp, Gauge, PieceKey,
};
mod transport;
use transport::*;
pub(super) use transport::{assemble_rule, canonical_pieces};
pub use transport::{canonical_recurrence, same_canonical_recurrence, transformed_rule_holds_at};
