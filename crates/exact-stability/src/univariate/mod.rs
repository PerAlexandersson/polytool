//! Möbius reduction, integer real/imaginary split and exact root counts.

pub mod bezoutian;
pub mod mobius;
mod polynomial;
pub mod region;

use crate::{BigRational as Q, GaussianRational as G};
use exact_algebra::inertia::{symmetric_ldlt, SymmetricLdlt};
use mobius::TransformedPolynomial;
use num_bigint::BigInt;
use region::StabilityRegion;
use std::fmt;

/// How a verdict is read from exact counts in an open region.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum StabilityConvention {
    /// No roots in the open region; the zero polynomial is stable by convention.
    Avoid,
    /// All roots strictly inside; the zero polynomial fails this convention.
    Inside,
}

/// Root counts of a nonzero polynomial, always including multiplicities.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct RootCounts {
    /// Roots in the open region.
    pub inside: usize,
    /// Roots on its boundary, including a root lost at the transformed infinity.
    pub boundary: usize,
    /// Roots in the opposite open region.
    pub outside: usize,
}

/// Exact counting method used for this input.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum StabilityMethod {
    /// The zero-polynomial convention; finite counts do not exist.
    ZeroPolynomial,
    /// Integer Sturm counts alone suffice after the Möbius split.
    RealSturm,
    /// Coprime Bézoutian inertia plus integer Sturm counts of the common factor.
    BezoutianSturm,
}

/// One squarefree-filtration layer used to count real roots with multiplicity.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RealRootLayer {
    /// Primitive integer polynomial in ascending order.
    pub polynomial: Vec<BigInt>,
    /// Exact number of distinct real roots, from integer Sturm/PRS.
    pub distinct_real_roots: usize,
}

/// Exact data supporting the root counts of a nonzero polynomial.
///
/// Replaying the computation checks every field, including the transformed
/// degree, common factor, Sturm layers, Bézoutian and pivoted LDLᵀ certificate.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct StabilityCertificate {
    /// Homogeneous Möbius pullback, before rotating its leading coefficient.
    pub transformed: TransformedPolynomial,
    /// Integer real part after multiplying by the conjugate leading coefficient
    /// and clearing both parts' denominators with a common positive multiplier.
    pub real_part: Vec<BigInt>,
    /// Integer imaginary part, of strictly smaller degree than `real_part`.
    pub imaginary_part: Vec<BigInt>,
    /// Primitive positive-leading gcd of the two parts over ℚ.
    pub common_factor: Vec<BigInt>,
    /// Repeated-root filtration; the sum of its distinct counts is the real
    /// root count of `common_factor`, including multiplicities.
    pub real_root_layers: Vec<RealRootLayer>,
    /// Bézoutian after dividing out the gcd and clearing denominators jointly.
    pub bezoutian: Vec<Vec<BigInt>>,
    /// Exact nonsingular block LDLᵀ of the Bézoutian (empty for constants).
    pub decomposition: SymmetricLdlt,
}

/// Exact result for a univariate polynomial and a specified open region.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct StabilityReport {
    /// Degree after trimming; `None` precisely for the zero polynomial.
    pub degree: Option<usize>,
    /// Region whose roots are counted.
    pub region: StabilityRegion,
    /// Meaning of the verdict; independent of the counts.
    pub convention: StabilityConvention,
    /// Exact finite counts; `None` only for the zero polynomial, which vanishes
    /// everywhere and does not have a finite multiset of roots.
    pub counts: Option<RootCounts>,
    /// Whether the chosen convention holds, including the documented zero rule.
    pub verdict: bool,
    /// Exact algorithm used.
    pub method: StabilityMethod,
    /// Replay data; absent only for the zero polynomial.
    pub certificate: Option<StabilityCertificate>,
}

impl StabilityReport {
    /// Recompute all data from the input and reject altered or malformed reports.
    ///
    /// This is deterministic exact replay, rather than a trust check on stored
    /// verdicts. Its cost is the same as constructing the report.
    pub fn verify(&self, coefficients: &[G]) -> bool {
        check_univariate_stability_with_convention(coefficients, &self.region, self.convention)
            .is_ok_and(|expected| expected == *self)
    }

    /// Human-readable meaning of `verdict`.
    pub fn meaning(&self) -> &'static str {
        match self.convention {
            StabilityConvention::Avoid => {
                "no root in the open region (zero polynomial stable by convention)"
            }
            StabilityConvention::Inside => {
                "all roots in the open region (zero polynomial excluded)"
            }
        }
    }
}

/// Invalid region data for an exact stability computation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum StabilityError {
    /// A fractional linear map has zero determinant.
    SingularMobiusTransform,
    /// A centered disk has radius zero or a negative radius.
    NonpositiveRadius,
}

impl fmt::Display for StabilityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::SingularMobiusTransform => "Möbius transform must have nonzero determinant",
            Self::NonpositiveRadius => "disk radius must be strictly positive",
        })
    }
}

impl std::error::Error for StabilityError {}

pub(super) fn trim_gaussian(p: &mut Vec<G>) {
    while p.last().is_some_and(G::is_zero) {
        p.pop();
    }
}

/// Count roots exactly and use the region's classical verdict convention.
///
/// Coefficients are ascending; trailing zeros are ignored. A nonzero constant
/// passes both conventions. Invalid radius data is rejected even for zero input.
/// There are no implicit resource caps in this library API; adapters should
/// validate input sizes before calling it.
pub fn check_univariate_stability(
    coefficients: &[G],
    region: &StabilityRegion,
) -> Result<StabilityReport, StabilityError> {
    check_univariate_stability_with_convention(coefficients, region, region.default_convention())
}

/// Count roots with multiplicity and apply an explicitly chosen convention.
///
/// Möbius substitution and Bézoutian construction use O(n²) field operations;
/// inertia uses O(n³) rational operations and O(n²) storage. Integer PRS
/// provides the real-root counts, including repeated common factors.
pub fn check_univariate_stability_with_convention(
    coefficients: &[G],
    region: &StabilityRegion,
    convention: StabilityConvention,
) -> Result<StabilityReport, StabilityError> {
    let transformed = region
        .mobius_transform()?
        .transform_polynomial(coefficients);
    let degree = transformed.original_degree;
    if degree.is_none() {
        return Ok(StabilityReport {
            degree,
            region: region.clone(),
            convention,
            counts: None,
            verdict: convention == StabilityConvention::Avoid,
            method: StabilityMethod::ZeroPolynomial,
            certificate: None,
        });
    }
    let leading_conjugate = transformed
        .coefficients
        .last()
        .expect("nonzero pullback")
        .conjugate();
    let rotated: Vec<G> = transformed
        .coefficients
        .iter()
        .map(|x| x.multiply(&leading_conjugate))
        .collect();
    let real: Vec<Q> = rotated.iter().map(|x| x.real.clone()).collect();
    let imaginary: Vec<Q> = rotated.iter().map(|x| x.imaginary.clone()).collect();
    let (real_part, imaginary_part) = polynomial::clear_pair(&real, &imaginary);
    let common_factor = polynomial::gcd(&real_part, &imaginary_part);
    let (real_roots, real_root_layers) =
        polynomial::count_real_roots_with_multiplicity(&common_factor);
    let conjugate_pairs = (common_factor.len() - 1 - real_roots) / 2;
    let real_quotient = polynomial::divide_exact(&real_part, &common_factor);
    let imaginary_quotient = polynomial::divide_exact(&imaginary_part, &common_factor);
    let (g, h) = polynomial::clear_pair(&real_quotient, &imaginary_quotient);
    let bezoutian = bezoutian::bezoutian_bigint_coeffs(&g, &h);
    let matrix: Vec<Vec<Q>> = bezoutian
        .iter()
        .map(|row| row.iter().cloned().map(Q::from_integer).collect())
        .collect();
    let decomposition = symmetric_ldlt(&matrix).expect("symmetric Bézoutian");
    let inertia = decomposition.inertia();
    assert_eq!(inertia.zero, 0, "coprime Bézoutian is nonsingular");
    let counts = RootCounts {
        inside: inertia.negative + conjugate_pairs,
        boundary: real_roots + transformed.boundary_at_infinity,
        outside: inertia.positive + conjugate_pairs,
    };
    assert_eq!(
        Some(counts.inside + counts.boundary + counts.outside),
        degree
    );
    let method = if imaginary_part.is_empty() {
        StabilityMethod::RealSturm
    } else {
        StabilityMethod::BezoutianSturm
    };
    let verdict = match convention {
        StabilityConvention::Avoid => counts.inside == 0,
        StabilityConvention::Inside => counts.boundary == 0 && counts.outside == 0,
    };
    Ok(StabilityReport {
        degree,
        region: region.clone(),
        convention,
        counts: Some(counts),
        verdict,
        method,
        certificate: Some(StabilityCertificate {
            transformed,
            real_part,
            imaginary_part,
            common_factor,
            real_root_layers,
            bezoutian,
            decomposition,
        }),
    })
}
