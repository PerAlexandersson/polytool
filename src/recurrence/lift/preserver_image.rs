//! Preserver-image lifts: recognise the rows `Q_n` of a recurrence-defined
//! sequence as images `Q_n = lambda(n) T_n(P_(n+s))` of a proved parent
//! family `P` under a real-rootedness preserver `T_n` from a fixed menu
//! ([`certify_preserver_image_lift`], [`find_preserver_image_lifts`]).
//!
//! # The certificate
//!
//! Let the input have the rule `M` (relation at base `n`,
//! `sum_(i,d) m_(i,d)(n, x) D^d Q_(n+i) = 0`, leader shift `R`, first index
//! `n0`, `K` identifying rows, so `M` holds at every base
//! `n >= b_M = n0 + K - R`), and the parent the rule `L` with first index
//! `p0` and a leader `g(n) S^k` free of `x` and of derivatives.  The
//! shifted parent `P'_n = P_(n+s)` (`n0 + s >= p0`) has the rule `L'` (the
//! substitution `n -> n + s`), which holds from a base `b_L` on.
//!
//! 1. **Gauge.**  For a Weyl form `T_n = T'_n / delta(n)` (below) put
//!    `mu = lambda / delta` (for other forms `delta = 1`), with ratio
//!    `mu(n)/mu(n-1) = alpha(n)/beta(n)`, `alpha = a(n) delta(n - 1)`,
//!    `beta = b(n) delta(n)` from the gauge ratio `a/b`; `delta(n) != 0` for
//!    every `n >= n0` is decided exactly (else the identity is only
//!    supported by rows).  Then for the candidate rows `H_n = mu(n) T'_n(P'_n)`
//!
//!    ```text
//!    (M H)_n prod_(t=1..R) beta(n + t) = mu(n) (M_mu T' P')_n,
//!    M_mu = sum m_(i,d)(n, x) prod_(t <= i) alpha(n + t) prod_(i < t <= R) beta(n + t) D^d S^i,
//!    ```
//!
//!    because `mu(n + i) = mu(n) prod_(t <= i) alpha(n + t)/beta(n + t)`
//!    and `S^i o T'(n) = T'(n + i) S^i` in the Ore algebra.  Both factors
//!    are nonzero for `n >= n0`.
//! 2. **Dividend.**  For a Weyl form the dividend is `E = M_mu o T'`.  For a
//!    coordinate map `M_mu` is conjugated through the map, one step at a
//!    time (last step first): for `G_n(x) = F_n(a x + b)`,
//!    `c(n, x) D^d G_(n+i) = [c(n, (y - b)/a) a^d D^d F_(n+i)](a x + b)`;
//!    for `G_n(x) = x^(d_n) F_n(1/x)` (`d_n` affine with step `h >= 0`),
//!    `D_x o C_n = C_n o (d_n y - y^2 D_y)`, `c(n, x) C_n(F) = C_n(c(n, 1/y) F)`
//!    and `C_(n+i)(F) = C_n(y^(-h i) F)`, so after multiplying by `y^N`,
//!    `N = b + h R` with `b` the degree in `x` of the coefficients of `M`
//!    (clearing negative powers), the term becomes
//!    `y^(N - h i) c(n, 1/y) (d_(n+i) y - y^2 D)^d S^i`.  A Möbius map with
//!    `c != 0` is the composition `F(y) -> F(a/c + kappa y)`, the reversal,
//!    and `K(u) -> K(c x + e)` with `kappa = (b c - a e)/c`.  Each
//!    substitution is injective (on Laurent polynomials for the reversal), so
//!    in every case `(M H)_n = 0` exactly when `(E P')_n = 0`.  These are
//!    the pullbacks of the `ore-algebra` crate with the bounds of `M` and
//!    top shift `R`: `ore_algebra::pull_back` through the gauge (and `T'`
//!    for a Weyl form), and for a coordinate map
//!    `ore_algebra::pull_back_substitution` through the map as an
//!    `ore_algebra::Substitution` (affine, reversal, Möbius with scale
//!    `1`), which is also the map's row action; `E` is then divided by the
//!    largest power of `x` dividing all its coefficients, which changes
//!    nothing above since `Q[x]` (and the Laurent polynomials) have no zero
//!    divisors.
//! 3. **Division.**  `E` is pseudo-divided by `L'`, whose leader is free of
//!    `x` (`ore_algebra::right_pseudo_division`, every level from the top
//!    shift of `E` down to `k`):
//!    `w(n) E = sum_nu H_nu S^nu L'` with remainder zero (exact Ore-algebra
//!    arithmetic).  So `(E P')_n = 0` at every base `n >= b_L` with
//!    `w(n) != 0`.  At the bases from `b_M` on that are below `b_L` or integer
//!    roots of `w` (finitely many) the input rule is evaluated exactly on the
//!    candidate rows.  Hence `H` satisfies `M` at every base `n >= b_M`.
//! 4. **Induction.**  `H_n = Q_n` on the `K` identifying rows (checked
//!    exactly), and the leader of `M` is a nonzero polynomial at every
//!    target `>= n0 + K` (decided exactly when the input is validated), so
//!    `H_n = Q_n` for every `n >= n0`.  (For a reversal or Möbius map with
//!    `deg P'_n > d_n` the candidate is a Laurent polynomial or a rational
//!    function; the argument holds in that field, and since `Q_n` is a
//!    polynomial the degree condition `deg P'_n <= d_n` follows.)
//!
//! A nonzero remainder `r` (of shift below `k`) is removed, when possible,
//! by row annihilators of `P'` exactly as in the `hadamard_factor` module
//! ("Row annihilators": the annihilator of a product parent, an invariant
//! annihilator of a first-order parent, or an invariant annihilator system
//! of a parent of any order, each replayed with its invariance identity and
//! the exact row checks at the start and restart bases):
//! the annihilation `mu r = sum_l Z_l A^(l)` with `mu = c(n) p(x) mu_2(n, x)`
//! turns `w E = sum_nu H_nu S^nu L' + r` into
//! `p(x) mu_2(n, x) (c w)(n) E = sum_nu (mu H_nu) S^nu L' + sum_l Z_l A^(l)`
//! (the index content `c` is folded into the recorded multiplier `w` and
//! into the recorded quotients) with
//! `A^(l) P' = 0` from the first used base on, so again `(E P')_n = 0` where
//! `w(n) != 0` and `mu_2(n, .) != 0` (`Q[x]` has no zero divisors); the
//! zeros of `mu_2` are exceptional bases.  An assembled identity or
//! certificate that fails its exact re-check leaves that candidate's
//! identity supported by rows only (and the search rejects any candidate
//! whose certification fails, with a note).  This is what certifies, for
//! example, `He_n + He_n'`, whose input rule uses the Hermite differential
//! equation.  Otherwise the identity is only supported by rows.
//! Hypergeometric multipliers ([`PreserverOperator::HadamardMultiplier`])
//! are certified by the Hadamard conjugation of the `hadamard_factor`
//! module against the shifted parent `P'`, with its row annihilators.
//!
//! [`PreserverDivisionCertificate::verify`] and
//! [`PreserverHadamardCertificate::verify`] replay everything from the
//! recorded rules and seeds: they regenerate the input and parent rows
//! (re-validating both rules), recompute `L'` and the dividend, check the
//! identity, the exceptional bases and the identifying rows on freshly
//! computed candidate rows.
//!
//! # The preserver menu
//!
//! The parent statement is that of the shared parent provers: every used
//! parent row is zero or real-rooted with every root in `(-inf, 0]`.  When
//! that is not proved and the operator's theorem needs only
//! real-rootedness, the three-term theorem is tried without a root region
//! (for example for the Hermite or Chebyshev parents); the input region is
//! then the whole line.
//!
//! **Interval parents.**  A parent may also be proved real-rooted with every
//! root in a bounded interval `[L, U]` (rational ends): by the Ma--Wang
//! interval check, or by the three-term theorem with a root-region
//! certificate (the `interval` module documentation).  For the substitution
//! maps and the Rolle operators the input region is then the exact image of
//! `[L, U]` (an affine image, a reciprocal image, a Möbius image, or the
//! hull with `0`), used when the ray statement is not proved or the image is
//! strictly smaller than the ray region; multipliers and `phi(D)` use their
//! ray theorems only, when `[L, U]` lies in `(-inf, 0]`, and otherwise (the
//! operator preserving real-rootedness on the line) the line.  The interval
//! report is recorded in the certification and replayed by `verify()`.  Each side condition is decided exactly for every
//! `n >= n0`; a side condition that cannot be decided leaves the lift open.
//!
//! - `f(n, theta)`, `theta = x D` (the multiplier `k -> f(n, k)`): Laguerre's
//!   theorem, `f` real-rooted with roots `<= 0` gives a multiplier sequence
//!   (RealRooted `isMultiplierSequence_eval_of_roots_nonpos`).  Side
//!   condition: `f(n, .)` real-rooted with roots `<= 0`, exact for a fixed
//!   `f`; for index-dependent coefficients of degree `<= 2` the discriminant
//!   and the pairwise coefficient signs by the sign checker; otherwise open.
//!   Region: `(-inf, 0]` maps to `(-inf, 0]`.
//! - `phi(n, D)`, for example `P + c(n) P'`: Laguerre's theorem, `phi(D)`
//!   with `phi` real-rooted preserves real-rootedness (`p + c p'` is
//!   real-rooted for real `c`; RealRooted `splits_add_C_mul_derivative_all`).
//!   Side condition: `phi(n, .)` real-rooted, exact for a fixed `phi`,
//!   automatic for degree `<= 1`, the discriminant for degree 2.  Region:
//!   `(-inf, 0]` when the coefficients of `phi` have one sign (decided), else
//!   the line.
//! - `D^j o x^m`, `x^m D^j`: Rolle; no side condition; `(-inf, 0]`.
//! - `P(a x + b)`, `a != 0`: substitution; `(-inf, -b/a]` for `a > 0`,
//!   `[-b/a, inf)` for `a < 0`.
//! - `x^(d_n) P(1/x)`: substitution; `(-inf, 0]`.
//! - `(c x + e)^(d_n) P((a x + b)/(c x + e))`, `c != 0`, `a e - b c != 0`:
//!   substitution; the interval between `-e/c` and `-b/a` when `a/c > 0` (the
//!   `f`-polynomial gives `[-1, 0]`), else the line.
//! - `x^k -> u(n, k) x^k` with `u` hypergeometric: Pólya--Schur (the
//!   `hadamard_factor` module); side condition: `u(n, .)` a nonnegative
//!   multiplier sequence (justified per factor); `(-inf, 0]`.
//! - the Euler step `kappa (theta + a)(theta + b) + x (u + v theta)`:
//!   RealRooted `IsPFPolynomial.generalStep_of_natDegree_le` (polynomials
//!   with nonnegative coefficients and real roots are preserved); side
//!   conditions `kappa, kappa (a + b), kappa a b > 0` with a nonnegative
//!   discriminant (so `a, b > 0`), `u > 0`, `v >= 0`, and a negative
//!   comparison defect on `s <= 0` (`B = u - (a + b + 1) v/2 >= 0`, or
//!   `B^2 < v^2 ((a + b + 1)^2/4 - a b)`), each by the sign checker; region
//!   `(-inf, 0]`, and nothing from parents rooted elsewhere;
//! - a general Weyl operator: no theorem, the lift stays open.
//!
//! Every certification records the Lean route (`reversed`, `fTransform`,
//! `derivativeXMul` where one matches) and the RealRooted (or Mathlib)
//! declarations behind its theorem, for a later replay in
//! <https://github.com/PerAlexandersson/RealRooted>.
//!
//! # Parents and detection
//!
//! [`PreserverParentLibrary`] holds candidate parents with their rows and
//! parent proofs (computed once per library): the named families of
//! [`crate::sequences`] with their defining recurrences, and optionally the
//! OEIS catalog entries, each kept when the parent provers, an interval
//! prover or the three-term theorem on the whole line prove it.
//! [`find_preserver_image_lifts`] tries every parent, index shift and
//! operator shape in a fixed order, fits the operator parameters exactly
//! from a few rows (linear algebra over the rationals; the index dependence
//! as a constant or a low-degree rational function), and certifies every
//! match with [`certify_preserver_image_lift`].
//!
//! # Interlacing
//!
//! Every certification carries a [`PreserverInterlacingClaim`]: consecutive
//! input rows interlace (weakly, orientation not determined) only when the
//! operator is independent of `n` and preserves real-rootedness of every
//! polynomial, and a replayed parent theorem report proves the parent rows
//! interlacing; the Hermite--Kakeya--Obreschkoff argument transfers it (the
//! `interlacing` module documentation).  Otherwise the claim is `none`.
//!
//! # Image rules
//!
//! A candidate is also an image specification of
//! [`crate::recurrence::image_recurrence_of_operator`]
//! ([`PreserverImageCandidate::image_operator_specification`]; the operator
//! by `ImageRowMap::try_from`, every operator except the multipliers that
//! depend on `n`), so the runnable rule of its image can be computed from
//! the parent rule ([`PreserverImageCandidate::image_recurrence`]), with
//! its own certificate.  On request ([`PreserverImageOptions::image_rule`])
//! a certified identity carries that rule
//! ([`PreserverImageCertification::image_rule`]); `verify()` replays it and
//! compares its first rows with the input's.  This is how a forward move
//! gives a new family its own rule.
//!
//! The group is split into `operator` (the menu, its normal forms and the
//! conjugation), `transfer` (theorems, side conditions, regions, Lean
//! names), `interval` (interval parents and interval images),
//! `interlacing` (the interlacing claim), `certify` (candidate,
//! certificates, certification), `summary` (the [`LiftSummary`] mapping),
//! `library` (the parents) and `search` (detection).

mod certify;
mod interlacing;
mod interval;
pub use interval::find_bounded_root_region_report;
mod library;
mod operator;
mod search;
mod summary;
mod transfer;
pub(in crate::recurrence) use interval::region_contains;
pub(in crate::recurrence) use transfer::{operator_transfer, polynomial_hypotheses};

pub use certify::{
    certify_preserver_image_lift, PreserverDivisionCertificate, PreserverDivisionQuotient,
    PreserverHadamardCertificate, PreserverImageCandidate, PreserverImageCertificate,
    PreserverImageCertification, PreserverImageIdentityStatus, PreserverImageOptions,
    MAX_PRESERVER_IMAGE_CHECKED_ROWS, MAX_PRESERVER_IMAGE_INDEX_SHIFT,
    MAX_PRESERVER_IMAGE_OPERATOR_TERMS, PRESERVER_IMAGE_CERTIFICATION_SCHEMA,
};
pub use interlacing::PreserverInterlacingClaim;
pub use library::{
    PreserverParent, PreserverParentLibrary, PreserverParentLibraryOptions, PreserverParentRegion,
};
pub use operator::{IndexRationalFunction, PreserverOperator, MAX_PRESERVER_OPERATOR_DEGREE};
pub use search::{
    find_preserver_image_lifts, PreserverImageReport, PreserverImageSearchOptions,
    PreserverOperatorShape, PRESERVER_IMAGE_LIFTS_SCHEMA,
};
