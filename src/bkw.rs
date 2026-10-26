//! Beraha--Kahane--Weiss scout tools for polynomial recurrences.
//!
//! This module is deliberately a scout layer, not a proof layer.  It evaluates a
//! characteristic symbol
//!
//! ```text
//! F(x, z) = a_0(x) + a_1(x) z + ... + a_r(x) z^r
//! ```
//!
//! at complex `x`, computes the characteristic roots in `z`, and ranks points
//! where the two dominant characteristic roots have nearly equal modulus.
//! Such points are numerical candidates for off-real BKW accumulation loci.
//!
//! Entry points: [`BkwSymbol::from_z_coefficients`] /
//! [`BkwSymbol::parse_z_coefficient_symbol`] build a symbol,
//! [`BkwSymbol::roots_at`] and [`BkwSymbol::dominant_tie_at`] analyze one point
//! `x`, and [`BkwSymbol::scout_equal_modulus_locus`] scans a grid.
//!
//! Conventions: symbol coefficients are exact `i64` polynomials in `x`
//! (ascending degree, entry `i` is the coefficient of `x^i`) and are indexed by
//! ascending power of `z`.  All root computations are floating point
//! (Durand--Kerner iteration in `f64`); nothing here is exact or certified.  Use
//! [`BkwSymbol::mathematica_equal_modulus_query`] to produce an exact follow-up
//! query.  For the recurrences that produce such symbols see `crate::recurrence`.

use crate::parse::parse_polynomial;
use crate::real_rootedness::format_poly_var;
use crate::vec_poly::strip_trailing_zeros;
use std::cmp::Ordering;
use std::error::Error;
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Small complex-number type used to keep `polytool` dependency-light.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Complex64 {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
}

impl Complex64 {
    /// The complex number `0`.
    pub const ZERO: Self = Self { re: 0.0, im: 0.0 };
    /// The complex number `1`.
    pub const ONE: Self = Self { re: 1.0, im: 0.0 };

    /// Build `re + im*i`.
    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    /// Build `radius * (cos(theta) + i sin(theta))`; `theta` is in radians.
    pub fn from_polar(radius: f64, theta: f64) -> Self {
        Self {
            re: radius * theta.cos(),
            im: radius * theta.sin(),
        }
    }

    /// Squared modulus `re^2 + im^2`.
    pub fn abs_squared(self) -> f64 {
        self.re.mul_add(self.re, self.im * self.im)
    }

    /// Modulus `sqrt(re^2 + im^2)`.
    pub fn abs(self) -> f64 {
        self.abs_squared().sqrt()
    }

    /// Whether the modulus is at most `tolerance`.
    pub fn is_near_zero(self, tolerance: f64) -> bool {
        self.abs() <= tolerance
    }
}

impl fmt::Display for Complex64 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.im == 0.0 {
            write!(f, "{:.12}", self.re)
        } else if self.re == 0.0 {
            write!(f, "{:.12}i", self.im)
        } else if self.im < 0.0 {
            write!(f, "{:.12} - {:.12}i", self.re, -self.im)
        } else {
            write!(f, "{:.12} + {:.12}i", self.re, self.im)
        }
    }
}

impl Add for Complex64 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.re + rhs.re, self.im + rhs.im)
    }
}

impl Sub for Complex64 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.re - rhs.re, self.im - rhs.im)
    }
}

impl Mul for Complex64 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )
    }
}

impl Div for Complex64 {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let denominator = rhs.abs_squared();
        Self::new(
            (self.re * rhs.re + self.im * rhs.im) / denominator,
            (self.im * rhs.re - self.re * rhs.im) / denominator,
        )
    }
}

impl Neg for Complex64 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(-self.re, -self.im)
    }
}

impl Mul<f64> for Complex64 {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::new(self.re * rhs, self.im * rhs)
    }
}

impl Div<f64> for Complex64 {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self::new(self.re / rhs, self.im / rhs)
    }
}

/// Errors reported by the BKW scout layer.
#[derive(Debug, Clone, PartialEq)]
pub enum BkwError {
    /// The symbol has no nonzero coefficient.
    EmptySymbol,
    /// The symbol has z-degree 0, or fewer than two characteristic roots were
    /// available at the requested point (no tie can be analyzed).
    ConstantSymbol,
    /// A z-coefficient could not be parsed as a polynomial in `x`; the payload is
    /// the parser message.
    ParseSymbol(String),
    /// The leading z-coefficient is numerically zero at the requested `x`
    /// (within [`BkwRootOptions::leading_zero_tolerance`]), so the number of
    /// roots would drop.
    LeadingCoefficientVanishes,
    /// The Durand--Kerner iteration produced a non-finite residual.
    RootFinderFailure(String),
}

impl fmt::Display for BkwError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySymbol => write!(f, "empty BKW symbol"),
            Self::ConstantSymbol => write!(f, "BKW symbol must have positive z-degree"),
            Self::ParseSymbol(e) => write!(f, "could not parse BKW symbol: {}", e),
            Self::LeadingCoefficientVanishes => {
                write!(f, "leading z-coefficient vanishes at this x")
            }
            Self::RootFinderFailure(e) => write!(f, "root finder failed: {}", e),
        }
    }
}

impl Error for BkwError {}

/// Characteristic symbol `F(x,z)`, stored by ascending powers of `z`.
///
/// Each z-coefficient is a polynomial in `x`, again in ascending coefficient
/// order.  For example,
///
/// ```text
/// F(x,z) = 1 - x z + z^2
/// ```
///
/// is represented by `[[1], [0, -1], [1]]`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BkwSymbol {
    z_coefficients: Vec<Vec<i64>>,
}

impl BkwSymbol {
    /// Build a symbol from z-coefficients, trimming zero tails in both variables.
    pub fn from_z_coefficients(coefficients: Vec<Vec<i64>>) -> Result<Self, BkwError> {
        let mut z_coefficients = coefficients;
        for coefficient in &mut z_coefficients {
            strip_trailing_zeros(coefficient);
        }
        while matches!(z_coefficients.last(), Some(c) if c.is_empty()) {
            z_coefficients.pop();
        }
        if z_coefficients.is_empty() {
            return Err(BkwError::EmptySymbol);
        }
        if z_coefficients.len() == 1 {
            return Err(BkwError::ConstantSymbol);
        }
        Ok(Self { z_coefficients })
    }

    /// Parse a symbol from semicolon-separated or line-separated z-coefficients.
    ///
    /// Each piece is parsed as a univariate polynomial in `x` using the crate's
    /// ordinary polynomial parser.  Examples:
    ///
    /// ```text
    /// 1; -x; 1
    /// ```
    ///
    /// means `F(x,z)=1-xz+z^2`.
    pub fn parse_z_coefficient_symbol(input: &str) -> Result<Self, BkwError> {
        let pieces: Vec<&str> = if input.contains(';') {
            input.split(';').collect()
        } else {
            input.lines().collect()
        };
        let mut coefficients = Vec::new();
        for piece in pieces {
            let trimmed = piece.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            coefficients.push(parse_polynomial(trimmed).map_err(BkwError::ParseSymbol)?);
        }
        Self::from_z_coefficients(coefficients)
    }

    /// z-degree of the symbol.
    pub fn z_degree(&self) -> usize {
        self.z_coefficients.len() - 1
    }

    /// Raw z-coefficients, each a polynomial in x.
    pub fn z_coefficients(&self) -> &[Vec<i64>] {
        &self.z_coefficients
    }

    /// Evaluate `F(x,z)` as a polynomial in z at a numerical complex `x`.
    pub fn z_polynomial_at(&self, x: Complex64) -> Vec<Complex64> {
        self.z_coefficients
            .iter()
            .map(|coeff| evaluate_i64_polynomial_complex(coeff, x))
            .collect()
    }

    /// Compute characteristic roots at `x`.
    pub fn roots_at(
        &self,
        x: Complex64,
        options: &BkwRootOptions,
    ) -> Result<BkwRootComputation, BkwError> {
        let mut coeffs = self.z_polynomial_at(x);
        strip_near_zero_complex_tail(&mut coeffs, options.leading_zero_tolerance);
        if coeffs.len() <= 1 {
            return Err(BkwError::LeadingCoefficientVanishes);
        }
        durand_kerner_roots(&coeffs, options)
    }

    /// Analyze the dominant-root modulus tie at `x`.
    pub fn dominant_tie_at(
        &self,
        x: Complex64,
        options: &BkwRootOptions,
    ) -> Result<BkwScoutCandidate, BkwError> {
        let root_computation = self.roots_at(x, options)?;
        let mut roots_by_modulus: Vec<BkwRootInfo> = root_computation
            .roots
            .iter()
            .enumerate()
            .map(|(index, &root)| BkwRootInfo {
                index,
                root,
                modulus: root.abs(),
            })
            .collect();
        roots_by_modulus.sort_by(|a, b| descending_f64(a.modulus, b.modulus));

        if roots_by_modulus.len() < 2 {
            return Err(BkwError::ConstantSymbol);
        }

        let first = roots_by_modulus[0].modulus;
        let second = roots_by_modulus[1].modulus;
        let scale = first.max(second).max(1.0);
        let relative_modulus_gap = (first - second).abs() / scale;
        let log_modulus_gap = if first > 0.0 && second > 0.0 {
            (first.ln() - second.ln()).abs()
        } else {
            relative_modulus_gap
        };
        let dominance_ratio = roots_by_modulus.get(2).and_then(|third| {
            if third.modulus == 0.0 {
                Some(f64::INFINITY)
            } else if second > 0.0 {
                Some(second / third.modulus)
            } else {
                None
            }
        });

        Ok(BkwScoutCandidate {
            x,
            z_degree_at_x: roots_by_modulus.len(),
            roots_by_modulus,
            relative_modulus_gap,
            log_modulus_gap,
            dominance_ratio,
            root_residual: root_computation.max_residual,
            converged: root_computation.converged,
            iterations: root_computation.iterations,
        })
    }

    /// Scan a rectangular complex box for candidate dominant equal-modulus ties.
    pub fn scout_equal_modulus_locus(&self, options: &BkwScoutOptions) -> Vec<BkwScoutCandidate> {
        let grid_re = options.grid_re.max(1);
        let grid_im = options.grid_im.max(1);
        let re_step = if grid_re == 1 {
            0.0
        } else {
            (options.re_max - options.re_min) / (grid_re - 1) as f64
        };
        let im_step = if grid_im == 1 {
            0.0
        } else {
            (options.im_max - options.im_min) / (grid_im - 1) as f64
        };

        let mut candidates = Vec::new();
        let retention = options
            .max_results
            .saturating_mul(8)
            .max(options.max_results);
        for i in 0..grid_re {
            let re = if grid_re == 1 {
                (options.re_min + options.re_max) / 2.0
            } else {
                options.re_min + i as f64 * re_step
            };
            for j in 0..grid_im {
                let im = if grid_im == 1 {
                    (options.im_min + options.im_max) / 2.0
                } else {
                    options.im_min + j as f64 * im_step
                };
                if !options.include_real_axis && im.abs() < options.min_imaginary_abs {
                    continue;
                }
                let x = Complex64::new(re, im);
                if let Ok(candidate) = self.dominant_tie_at(x, &options.root_options) {
                    let candidate = if options.refine_steps == 0 {
                        candidate
                    } else {
                        self.refine_candidate(candidate, options, re_step.abs().max(im_step.abs()))
                    };
                    candidates.push(candidate);
                }
            }
            if candidates.len() > retention.saturating_mul(2) {
                sort_candidates(&mut candidates);
                candidates.truncate(retention);
            }
        }
        sort_candidates(&mut candidates);
        candidates.truncate(options.max_results);
        candidates
    }

    fn refine_candidate(
        &self,
        start: BkwScoutCandidate,
        options: &BkwScoutOptions,
        initial_step: f64,
    ) -> BkwScoutCandidate {
        if initial_step == 0.0 {
            return start;
        }
        let mut best = start;
        let mut step = initial_step;
        let directions = [
            (1.0, 0.0),
            (-1.0, 0.0),
            (0.0, 1.0),
            (0.0, -1.0),
            (1.0, 1.0),
            (1.0, -1.0),
            (-1.0, 1.0),
            (-1.0, -1.0),
        ];

        for _ in 0..options.refine_steps {
            let mut improved = false;
            for (dx, dy) in directions {
                let x = Complex64::new(best.x.re + dx * step, best.x.im + dy * step);
                if x.re < options.re_min
                    || x.re > options.re_max
                    || x.im < options.im_min
                    || x.im > options.im_max
                {
                    continue;
                }
                if !options.include_real_axis && x.im.abs() < options.min_imaginary_abs {
                    continue;
                }
                if let Ok(candidate) = self.dominant_tie_at(x, &options.root_options) {
                    if candidate_is_better(&candidate, &best) {
                        best = candidate;
                        improved = true;
                    }
                }
            }
            if !improved {
                step *= 0.5;
            }
        }
        best
    }

    /// Mathematica `Reduce` skeleton for exact equal-modulus follow-up.
    ///
    /// This excludes dominance and amplitudes; it encodes two distinct roots of
    /// the same modulus at an off-real `x = u + I v`.
    pub fn mathematica_equal_modulus_query(&self) -> String {
        let f = self.format_mathematica_symbol();
        [
            "Clear[x, z, u, v, a, b, c, d];".to_string(),
            format!("F[x_, z_] := {};", f),
            "Reduce[".to_string(),
            "  ComplexExpand[".to_string(),
            "    Re[F[u + I v, a + I b]] == 0 &&".to_string(),
            "    Im[F[u + I v, a + I b]] == 0 &&".to_string(),
            "    Re[F[u + I v, c + I d]] == 0 &&".to_string(),
            "    Im[F[u + I v, c + I d]] == 0 &&".to_string(),
            "    a^2 + b^2 == c^2 + d^2 &&".to_string(),
            "    (a - c)^2 + (b - d)^2 > 0 &&".to_string(),
            "    v != 0".to_string(),
            "  ],".to_string(),
            "  {u, v, a, b, c, d}, Reals".to_string(),
            "]".to_string(),
        ]
        .join("\n")
    }

    /// Human-readable symbol formatting.
    pub fn format_symbol(&self) -> String {
        let mut terms = Vec::new();
        for (z_degree, x_coeffs) in self.z_coefficients.iter().enumerate() {
            if x_coeffs.is_empty() {
                continue;
            }
            let coeff = format_poly_var(x_coeffs, "x");
            let term = match z_degree {
                0 => coeff,
                1 if coeff == "1" => "z".to_string(),
                1 if coeff == "-1" => "-z".to_string(),
                1 => format!("({}) z", coeff),
                k if coeff == "1" => format!("z^{}", k),
                k if coeff == "-1" => format!("-z^{}", k),
                k => format!("({}) z^{}", coeff, k),
            };
            terms.push(term);
        }
        if terms.is_empty() {
            return "0".to_string();
        }
        join_signed_terms(&terms)
    }

    fn format_mathematica_symbol(&self) -> String {
        let mut terms = Vec::new();
        for (z_degree, x_coeffs) in self.z_coefficients.iter().enumerate() {
            if x_coeffs.is_empty() {
                continue;
            }
            let coeff = format_mathematica_polynomial(x_coeffs, "x");
            let term = match z_degree {
                0 => coeff,
                1 => format!("({})*z", coeff),
                k => format!("({})*z^{}", coeff, k),
            };
            terms.push(term);
        }
        if terms.is_empty() {
            "0".to_string()
        } else {
            terms.join(" + ")
        }
    }
}

/// Numerical root-finder options for characteristic roots.
#[derive(Clone, Debug)]
pub struct BkwRootOptions {
    /// Maximum number of Durand--Kerner sweeps (default 300).
    pub max_iterations: usize,
    /// Stop when the largest root update in a sweep is at most this value
    /// (absolute, default `1e-12`).
    pub tolerance: f64,
    /// Leading z-coefficients with modulus at most this value are treated as
    /// zero and dropped (default `1e-10`).
    pub leading_zero_tolerance: f64,
}

impl Default for BkwRootOptions {
    fn default() -> Self {
        Self {
            max_iterations: 300,
            tolerance: 1e-12,
            leading_zero_tolerance: 1e-10,
        }
    }
}

/// Rectangle/grid options for BKW scout mode.
#[derive(Clone, Debug)]
pub struct BkwScoutOptions {
    /// Lower bound of `Re(x)` for the scanned box (default `-3`).
    pub re_min: f64,
    /// Upper bound of `Re(x)` for the scanned box (default `3`).
    pub re_max: f64,
    /// Lower bound of `Im(x)` for the scanned box (default `-3`).
    pub im_min: f64,
    /// Upper bound of `Im(x)` for the scanned box (default `3`).
    pub im_max: f64,
    /// Number of grid points along the real axis; values below 1 are treated as
    /// 1 (a single point at the midpoint).  Default 61.
    pub grid_re: usize,
    /// Number of grid points along the imaginary axis; same convention as
    /// `grid_re`.  Default 61.
    pub grid_im: usize,
    /// Number of best candidates returned (default 20).
    pub max_results: usize,
    /// If `false`, grid points with `|Im(x)| < min_imaginary_abs` are skipped,
    /// since real-axis ties are not off-real accumulation points (default
    /// `false`).
    pub include_real_axis: bool,
    /// Threshold used by `include_real_axis = false` (default `1e-9`).
    pub min_imaginary_abs: f64,
    /// Number of local pattern-search refinement rounds applied to each grid
    /// candidate; `0` disables refinement (default 8).  Refinement stays inside
    /// the box.
    pub refine_steps: usize,
    /// Options for the characteristic root finder at each point.
    pub root_options: BkwRootOptions,
}

impl Default for BkwScoutOptions {
    fn default() -> Self {
        Self {
            re_min: -3.0,
            re_max: 3.0,
            im_min: -3.0,
            im_max: 3.0,
            grid_re: 61,
            grid_im: 61,
            max_results: 20,
            include_real_axis: false,
            min_imaginary_abs: 1e-9,
            refine_steps: 8,
            root_options: BkwRootOptions::default(),
        }
    }
}

/// A root together with its modulus and original root index.
#[derive(Clone, Debug)]
pub struct BkwRootInfo {
    /// Position of this root in [`BkwRootComputation::roots`] (the root finder's
    /// order, before sorting by modulus).
    pub index: usize,
    /// The characteristic root `z`.
    pub root: Complex64,
    /// `|root|`.
    pub modulus: f64,
}

/// Root computation metadata.
#[derive(Clone, Debug)]
pub struct BkwRootComputation {
    /// Characteristic roots `z` of `F(x, z)` at the evaluated `x`, in the
    /// root finder's (unsorted) order; the count is the z-degree at `x`.
    pub roots: Vec<Complex64>,
    /// Whether the largest update fell below [`BkwRootOptions::tolerance`]
    /// within [`BkwRootOptions::max_iterations`] sweeps.
    pub converged: bool,
    /// Number of sweeps performed.
    pub iterations: usize,
    /// Largest `|F(x, root)|` over the roots; a floating-point accuracy gauge.
    pub max_residual: f64,
}

/// Candidate point where the two dominant roots have nearly equal modulus.
#[derive(Clone, Debug)]
pub struct BkwScoutCandidate {
    /// The evaluation point `x`.
    pub x: Complex64,
    /// Number of characteristic roots at `x` (z-degree after dropping a
    /// numerically vanishing leading coefficient).
    pub z_degree_at_x: usize,
    /// All characteristic roots, sorted by decreasing modulus.
    pub roots_by_modulus: Vec<BkwRootInfo>,
    /// `| |z1| - |z2| | / max(|z1|, |z2|, 1)` for the two largest moduli;
    /// small values indicate a near tie.  Candidates are ranked by this value.
    pub relative_modulus_gap: f64,
    /// `| ln|z1| - ln|z2| |`, or `relative_modulus_gap` if a modulus is zero.
    pub log_modulus_gap: f64,
    /// `|z2| / |z3|`, how far the tied pair dominates the third root; `None`
    /// when there are only two roots, and infinity when `|z3| = 0`.
    pub dominance_ratio: Option<f64>,
    /// Maximum residual `|F(x, z)|` of the root computation.
    pub root_residual: f64,
    /// Whether the root finder met its tolerance.
    pub converged: bool,
    /// Sweeps used by the root finder.
    pub iterations: usize,
}

impl BkwScoutCandidate {
    /// The two largest-modulus roots, or `None` if there are fewer than two.
    pub fn tied_roots(&self) -> Option<(&BkwRootInfo, &BkwRootInfo)> {
        if self.roots_by_modulus.len() < 2 {
            None
        } else {
            Some((&self.roots_by_modulus[0], &self.roots_by_modulus[1]))
        }
    }
}

fn evaluate_i64_polynomial_complex(coeffs: &[i64], x: Complex64) -> Complex64 {
    let mut result = Complex64::ZERO;
    for &coeff in coeffs.iter().rev() {
        result = result * x + Complex64::new(coeff as f64, 0.0);
    }
    result
}

fn evaluate_complex_polynomial(coeffs: &[Complex64], z: Complex64) -> Complex64 {
    let mut result = Complex64::ZERO;
    for &coeff in coeffs.iter().rev() {
        result = result * z + coeff;
    }
    result
}

fn strip_near_zero_complex_tail(coeffs: &mut Vec<Complex64>, tolerance: f64) {
    while matches!(coeffs.last(), Some(c) if c.is_near_zero(tolerance)) {
        coeffs.pop();
    }
}

fn durand_kerner_roots(
    coeffs: &[Complex64],
    options: &BkwRootOptions,
) -> Result<BkwRootComputation, BkwError> {
    let degree = coeffs.len() - 1;
    let leading = coeffs[degree];
    if leading.is_near_zero(options.leading_zero_tolerance) {
        return Err(BkwError::LeadingCoefficientVanishes);
    }
    if degree == 1 {
        let root = -coeffs[0] / leading;
        let residual = evaluate_complex_polynomial(coeffs, root).abs();
        return Ok(BkwRootComputation {
            roots: vec![root],
            converged: true,
            iterations: 0,
            max_residual: residual,
        });
    }

    let monic: Vec<Complex64> = coeffs.iter().map(|&c| c / leading).collect();
    let radius = 1.0
        + monic[..degree]
            .iter()
            .map(|c| c.abs())
            .fold(0.0_f64, f64::max);
    let angle_offset = 0.37;
    let mut roots: Vec<Complex64> = (0..degree)
        .map(|k| {
            let theta = angle_offset + 2.0 * std::f64::consts::PI * (k as f64) / (degree as f64);
            Complex64::from_polar(radius, theta)
        })
        .collect();

    let mut converged = false;
    let mut iterations = 0;
    for iter in 0..options.max_iterations {
        iterations = iter + 1;
        let old_roots = roots.clone();
        let mut max_delta = 0.0_f64;
        for i in 0..degree {
            let mut denominator = Complex64::ONE;
            for j in 0..degree {
                if i != j {
                    denominator = denominator * (old_roots[i] - old_roots[j]);
                }
            }
            if denominator.is_near_zero(options.tolerance * 1e-3) {
                denominator = denominator + Complex64::new(options.tolerance, options.tolerance);
            }
            let delta = evaluate_complex_polynomial(&monic, old_roots[i]) / denominator;
            roots[i] = old_roots[i] - delta;
            max_delta = max_delta.max(delta.abs());
        }
        if max_delta <= options.tolerance {
            converged = true;
            break;
        }
    }

    let max_residual = roots
        .iter()
        .map(|&root| evaluate_complex_polynomial(coeffs, root).abs())
        .fold(0.0_f64, f64::max);
    if !converged && !max_residual.is_finite() {
        return Err(BkwError::RootFinderFailure(
            "nonfinite residual in Durand--Kerner iteration".to_string(),
        ));
    }
    Ok(BkwRootComputation {
        roots,
        converged,
        iterations,
        max_residual,
    })
}

fn sort_candidates(candidates: &mut [BkwScoutCandidate]) {
    candidates.sort_by(|a, b| {
        a.relative_modulus_gap
            .partial_cmp(&b.relative_modulus_gap)
            .unwrap_or(Ordering::Equal)
            .then_with(|| match (a.dominance_ratio, b.dominance_ratio) {
                (Some(ra), Some(rb)) => descending_f64(ra, rb),
                (Some(_), None) => Ordering::Less,
                (None, Some(_)) => Ordering::Greater,
                (None, None) => Ordering::Equal,
            })
            .then_with(|| {
                a.root_residual
                    .partial_cmp(&b.root_residual)
                    .unwrap_or(Ordering::Equal)
            })
    });
}

fn candidate_is_better(candidate: &BkwScoutCandidate, best: &BkwScoutCandidate) -> bool {
    candidate.relative_modulus_gap < best.relative_modulus_gap
        || (candidate.relative_modulus_gap == best.relative_modulus_gap
            && candidate.root_residual < best.root_residual)
}

fn descending_f64(a: f64, b: f64) -> Ordering {
    b.partial_cmp(&a).unwrap_or(Ordering::Equal)
}

fn join_signed_terms(terms: &[String]) -> String {
    let mut result = terms[0].clone();
    for term in &terms[1..] {
        if let Some(rest) = term.strip_prefix('-') {
            result.push_str(" - ");
            result.push_str(rest);
        } else {
            result.push_str(" + ");
            result.push_str(term);
        }
    }
    result
}

fn format_mathematica_polynomial(coeffs: &[i64], var: &str) -> String {
    let mut terms = Vec::new();
    for (degree, &coeff) in coeffs.iter().enumerate() {
        if coeff == 0 {
            continue;
        }
        let abs_coeff = coeff.abs();
        let factor = match degree {
            0 => abs_coeff.to_string(),
            1 if abs_coeff == 1 => var.to_string(),
            1 => format!("{}*{}", abs_coeff, var),
            k if abs_coeff == 1 => format!("{}^{}", var, k),
            k => format!("{}*{}^{}", abs_coeff, var, k),
        };
        if coeff < 0 {
            terms.push(format!("-{}", factor));
        } else {
            terms.push(factor);
        }
    }
    if terms.is_empty() {
        "0".to_string()
    } else {
        join_signed_terms(&terms)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_z_coefficient_symbol() {
        let symbol = BkwSymbol::parse_z_coefficient_symbol("1; -x; 1").unwrap();
        assert_eq!(symbol.z_degree(), 2);
        assert_eq!(symbol.z_coefficients(), &[vec![1], vec![0, -1], vec![1]]);
        assert_eq!(symbol.format_symbol(), "1 + (-x) z + z^2");
    }

    #[test]
    fn computes_quadratic_roots() {
        let symbol = BkwSymbol::parse_z_coefficient_symbol("1; 0; 1").unwrap();
        let roots = symbol
            .roots_at(Complex64::new(2.0, 3.0), &BkwRootOptions::default())
            .unwrap();
        assert!(roots.converged);
        assert_eq!(roots.roots.len(), 2);
        assert!(roots.max_residual < 1e-8);
        let mut moduli: Vec<f64> = roots.roots.iter().map(|root| root.abs()).collect();
        moduli.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((moduli[0] - 1.0).abs() < 1e-8);
        assert!((moduli[1] - 1.0).abs() < 1e-8);
    }

    #[test]
    fn scout_finds_off_real_constant_tie() {
        let symbol = BkwSymbol::parse_z_coefficient_symbol("1; 0; 1").unwrap();
        let options = BkwScoutOptions {
            re_min: -1.0,
            re_max: 1.0,
            im_min: -1.0,
            im_max: 1.0,
            grid_re: 3,
            grid_im: 3,
            max_results: 4,
            include_real_axis: false,
            refine_steps: 0,
            ..BkwScoutOptions::default()
        };
        let candidates = symbol.scout_equal_modulus_locus(&options);
        assert!(!candidates.is_empty());
        assert!(candidates
            .iter()
            .any(|candidate| candidate.x.im.abs() > 0.5 && candidate.relative_modulus_gap < 1e-8));
    }

    #[test]
    fn mathematica_query_contains_symbol() {
        let symbol = BkwSymbol::parse_z_coefficient_symbol("1; -x; 1").unwrap();
        let query = symbol.mathematica_equal_modulus_query();
        assert!(query.contains("F[x_, z_]"));
        assert!(query.contains("a^2 + b^2 == c^2 + d^2"));
        assert!(query.contains("v != 0"));
    }
}
