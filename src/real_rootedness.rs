//! Real-rootedness testing, interlacing, log-concavity, and ultra-log-concavity.
//!
//! Functions accept polynomials as `i64` coefficient vectors in ascending degree order
//! (`coeffs[i]` is the coefficient of t^i); the `*_bigint_coeffs` variants take
//! `&[BigInt]` and are the safe choice when coefficients can exceed `i64`.
//! Verdicts are exact (no floating point) except [`real_roots`], which returns
//! rational approximations.  Interlacing functions return `Option<bool>`
//! (`None` when the degree relation does not allow the directed test; see
//! [`check_weak_interlacing`]).  Besides real-rootedness and interlacing this
//! module also holds palindromicity, gamma-positivity, the Stapledon
//! decomposition, resultants/discriminants, the Ehrhart <-> h* conversions and
//! the polynomial display helpers (`format_poly*`).
//!
//! # Default algorithm: primitive integer root counting
//!
//! The default real-rootedness functions use the primitive integer PRS path in
//! [`crate::root_count`].  For one-signed coefficient polynomials, this first
//! tries cheap coefficient tests and then counts positive roots after the
//! transformation `t ↦ -t`; otherwise it counts all real roots of the
//! squarefree part directly.  Interlacing uses Bézout matrices, which are
//! dramatically faster than rational Sturm chains for that task.
//!
//! ## Bézout matrix for interlacing
//!
//! Given polynomials f (degree d) and g (degree d−1), both with positive leading
//! coefficients, the **Bézout matrix** B(f, g) is the d×d symmetric matrix
//! whose (i, j) entry is the coefficient of x^i y^j in the bivariate polynomial
//!
//! ```text
//! (f(x)g(y) - f(y)g(x)) / (x - y).
//! ```
//!
//! There are two standard oriented versions.
//!
//! - If `deg(f) = deg(g) + 1`, then `g` interlaces `f` iff `B(f, g)` is
//!   positive semidefinite; in the coprime/strict case this is positive
//!   definite. See Kummer--Naldi--Plaumann, Thm. 2.13, citing
//!   Krein--Naimark, §2.2.
//! - If `deg(f) = deg(g)`, the same-degree alternation criterion uses the
//!   same Bézoutian with the orientation fixed by the argument order. Fisk,
//!   §9.21, Cor. 9.145 gives the positive definite/no-common-root form.
//!
//! In semidefinite/common-root cases, the common factor must be handled
//! separately: a positive semidefinite Bézoutian does not by itself certify
//! that a shared factor has only real roots. This is why
//! [`check_weak_interlacing`] divides out the gcd and explicitly verifies that
//! the gcd is real-rooted.
//!
//! This reduces interlacing to a single exact matrix definiteness check,
//! avoiding root isolation entirely.
//!
//! ## Bézout matrix for real-rootedness
//!
//! A polynomial f of degree d is real-rooted if and only if f' interlaces f. Since
//! f may have repeated roots, B(f, f') is only positive **semi**-definite for real-rooted f.
//!
//! ## Shared roots and same-degree interlacing
//!
//! The Bézout matrix detects **strict** interlacing only (deg difference exactly 1).
//! [`check_weak_interlacing`] extends this to handle:
//!
//! - **Shared roots**: divides out gcd(f, g) over ℚ, verifies the GCD is real-rooted,
//!   then checks strict interlacing of the reduced polynomials.
//! - **Same degree**: `check_weak_interlacing(f, g)` tests f ≪ g (f on the LEFT)
//!   by extending f with a root far to the right (via the Cauchy bound), reducing
//!   to the deg+1 case. If all coefficients of f are positive (all roots negative),
//!   multiplying by t suffices.
//! - **Directed**: `check_weak_interlacing(f, g)` tests `f ≪ g` (f interlaces g
//!   from the left). This is **not symmetric**: `f ≪ g` and `g ≪ f` are different
//!   for same-degree polynomials. Requires deg(f) ≤ deg(g) ≤ deg(f) + 1.
//! - **Undirected**: [`interlaces_either_order`] answers "do f and g interlace
//!   in some order allowed by their degrees?" and returns a plain `bool`.
//!
//! ## Exact root enclosures
//!
//! [`isolate_real_roots_exact`] returns rational isolating intervals (or exact
//! rational roots) with multiplicities, degree and valuation, and
//! [`RealInterval`] decides exactly whether all (real) roots lie in an
//! interval or ray.  Use it whenever a root location is a hypothesis of an
//! argument; [`real_roots`] only returns midpoint representatives.
//!
//! ## Sturm chain fallback
//!
//! The original Sturm-chain implementations are available as [`is_real_rooted_sturm`],
//! [`check_interlacing_sturm`], and [`real_roots`] for cases where actual root locations
//! are needed.
//!
//! References:
//!
//! - S. Fisk, *Polynomials, roots, and interlacing*, arXiv:math/0612833,
//!   §9.21, Cor. 9.145.
//! - M. Kummer, S. Naldi, and D. Plaumann, *Spectrahedral representations of
//!   plane hyperbolic curves*, arXiv:1807.10901, Thm. 2.13.
//! - M. G. Krein and M. A. Naimark, *The method of symmetric and Hermitian
//!   forms in the theory of the separation of the roots of algebraic
//!   equations*, Linear and Multilinear Algebra 10 (1981), §2.2.
//! - MathOverflow discussion of the common-factor caveat:
//!   <https://mathoverflow.net/questions/403708/bezout-matrices-and-interlacing-roots>

use num_bigint::BigInt;
use num_rational::Ratio;

// The implementation lives in private submodules; everything public is
// re-exported below at its historical `real_rootedness::*` path.
mod bezout;
mod coefficients;
mod ehrhart;
mod enclosures;
mod formatting;
mod interlacing;
mod qpoly;
mod real_rooted;
mod resultants;
mod sturm_fallback;

pub use bezout::{bezout_matrix, bezout_matrix_bigint_coeffs, bezoutian_bigint_coeffs};
pub use coefficients::{
    gamma_coefficients, gamma_coefficients_bigint_coeffs,
    gamma_coefficients_ignoring_initial_zeros,
    gamma_coefficients_ignoring_initial_zeros_bigint_coeffs, hermite_biehler_parts,
    is_gamma_positive, is_gamma_positive_bigint_coeffs, is_gamma_positive_ignoring_initial_zeros,
    is_gamma_positive_ignoring_initial_zeros_bigint_coeffs, is_log_concave,
    is_log_concave_bigint_coeffs, is_palindromic, is_palindromic_bigint_coeffs,
    is_palindromic_ignoring_initial_zeros, is_palindromic_ignoring_initial_zeros_bigint_coeffs,
    is_ultra_log_concave, is_ultra_log_concave_bigint_coeffs, is_ultra_log_concave_infinite,
    is_ultra_log_concave_infinite_bigint_coeffs, is_ultra_log_concave_of_degree,
    is_ultra_log_concave_of_degree_bigint_coeffs, is_unimodal, is_unimodal_bigint_coeffs,
    stapledon_decomposition, stapledon_decomposition_bigint_coeffs, strip_initial_zeros,
    strip_initial_zeros_bigint_coeffs,
};
pub use ehrhart::{
    ehrhart_to_hstar, ehrhart_to_hstar_bigint, ehrhart_to_hstar_with_denom, hstar_to_ehrhart,
    hstar_to_ehrhart_bigint_coeffs, EhrhartToHstarError,
};
pub use enclosures::{
    isolate_real_roots_exact, RealInterval, RealRootEnclosure, RealRootIsolation,
};
pub use formatting::{
    format_poly, format_poly_bigint_coeffs, format_poly_bigint_var, format_poly_var,
};
pub use interlacing::{
    check_interlacing, check_interlacing_bezout, check_interlacing_bigint_coeffs,
    check_weak_interlacing, check_weak_interlacing_bezout, check_weak_interlacing_bigint_coeffs,
    interlaces_either_order, interlaces_either_order_bigint_coeffs,
};
pub use real_rooted::{
    hermite_matrix_bigint_coeffs, is_real_rooted, is_real_rooted_bezout,
    is_real_rooted_bezout_bigint_coeffs, is_real_rooted_bezout_squarefree,
    is_real_rooted_bezout_squarefree_bigint_coeffs, is_real_rooted_bigint_coeffs,
    is_real_rooted_hermite, is_real_rooted_hermite_bigint_coeffs, is_strictly_real_rooted_bezout,
    is_strictly_real_rooted_bezout_bigint_coeffs, is_strictly_real_rooted_hermite,
    is_strictly_real_rooted_hermite_bigint_coeffs, squarefree_part_bigint_coeffs,
};
pub use resultants::{
    discriminant, discriminant_bigint_coeffs, resultant, resultant_bigint_coeffs, sylvester_matrix,
};
pub use sturm_fallback::{
    check_interlacing_sturm, has_simple_roots, has_simple_roots_bigint_coeffs,
    is_real_rooted_sturm, is_real_rooted_sturm_bigint_coeffs, real_roots,
};

type Q = Ratio<BigInt>;
