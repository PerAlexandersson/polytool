//! Gamma lifts: recognise the rows `Q_n` of a recurrence-defined sequence as
//! rescaled gamma polynomials of a palindromic parent,
//! `Q_n = lambda(n) Gamma_(n+s)(c z)`, and transfer real-rootedness from the
//! parent ([`certify_gamma_lift`], [`find_gamma_lifts`]).
//!
//! # The gamma polynomial
//!
//! A polynomial `P` with `x^d P(1/x) = P` (symmetric at index `d`, so
//! `deg P <= d`) has a unique expansion
//! `P(x) = sum_(j <= d/2) gamma_j x^j (1 + x)^(d - 2j)`; its gamma polynomial
//! is `Gamma(z) = sum_j gamma_j z^j`, so `P(x) = (1 + x)^d Gamma(x/(1 + x)^2)`
//! ([`RecurrenceTransformation::Gamma`]).  For the Eulerian polynomials
//! (`d_n = n - 1`) the coefficients count permutations without double
//! descents (Foata--Schützenberger), and `2^(n-1) Gamma_n(z/4)` is the peak
//! polynomial ([`crate::sequences::literature::peak_polynomials_bigint`]).
//!
//! # Root transfer
//!
//! **Theorem** (the gamma criterion; RealRooted,
//! `IsGammaExpansion.isRealRooted_and_hasRootsNonpos_iff`): if
//! `P = sum_j gamma_j x^j (1 + x)^(d - 2j)` with `deg Gamma <= d/2`, then `P`
//! is nonzero and real-rooted with every root `<= 0` if and only if `Gamma`
//! is.  The direction used here: write `Gamma(z) = g prod_i (z - zeta_i)`
//! with `e = deg Gamma`; then
//! `P(x) = g (1 + x)^(d - 2e) prod_i (x - zeta_i (1 + x)^2)`.  A root `x` of
//! `x - zeta (1 + x)^2` is never `-1` and satisfies `x/(1 + x)^2 = zeta`;
//! for real `x <= 0` that value is real and `<= 0`.  So a root `zeta` of
//! `Gamma` that is not real and `<= 0` gives roots `x` of `P` that are not
//! real and `<= 0`.  Hence if every root of `P != 0` is real and `<= 0`, so
//! is every root of `Gamma` (`Gamma != 0` because `P != 0`).  So, for
//! `c != 0` and
//! `lambda(n) != 0`, `Q_n(z) = lambda(n) Gamma_(n+s)(c z)` is zero or
//! real-rooted with every root in `(-inf, 0]` when `c > 0` (in `[0, inf)`
//! when `c < 0`) whenever `P_(n+s)` is zero or real-rooted with roots
//! `<= 0`.  The criterion places the gamma roots in `(-inf, 0]` and nothing
//! narrower: `1 + 3x + x^2 = (1 + x)^2 + x` has the gamma root `-1`.  The
//! parent statement is settled by the shared parent provers (the Veronese
//! module's `parent_proof`).
//!
//! # The certificate ([`GammaLiftCertificate`])
//!
//! Let the parent `P` have rule `L_P` with first index `p0`, `K_p`
//! identifying seeds and symmetry indices `d_m` (an affine
//! [`DegreeSchedule`]), and let the input have rule `M` (leader shift `R`)
//! with first index `n0` and `K` identifying seeds.
//!
//! 1. **Symmetry for every `m >= p0`.**  The reversal transport
//!    `x^(d_m) P_m(1/x)` of `L_P` must define the same row map as `L_P`
//!    (cross-multiplied coefficients agree), both left factors must be
//!    nonzero polynomials at every target `m >= p0 + K_p` (decided exactly),
//!    and the seeds must be symmetric.  The reversed rows (Laurent
//!    polynomials a priori) satisfy the reversed rule at every target, hence
//!    `L_P`, whose leader determines the row; by induction they are the rows,
//!    so every `P_m` is symmetric at `d_m`.  Without this step the identity
//!    is only supported by rows.
//! 2. **Gamma rule.**  [`transport_gamma`] gives a relation `T` of the gamma
//!    rows `G_m` (the scalar rule, or a part of the transported identity with
//!    a nonzero leader); it holds at every target `m >= p0 + K_p`, because
//!    there `L_P` holds and the rows involved are symmetric.
//! 3. **Scale, shift and gauge.**  The affine transport of `T` by `z -> c z`
//!    is a relation of `G_m(c z)`; the substitution `n -> n + s` and the
//!    gauge `lambda` ([`VeroneseSectionGauge`] in ratio form) give a relation
//!    `L` of the candidate rows `H_n = lambda(n) G_(n+s)(c z)` with leader
//!    `g(n) S^k` free of `x`, valid at every base
//!    `b >= max(n0, p0 + K_p - k - s)`.
//! 4. **Division.**  `M` is pseudo-divided by `L`:
//!    `w(n) M = sum_nu H_nu S^nu L` with remainder zero (exact operator
//!    algebra).  At every valid base that is not an integer root of `w`, `M`
//!    therefore vanishes on `H`; at the remaining bases from `n0 + K - R` on
//!    (finitely many) `M` is evaluated on `H` exactly.  So `H` satisfies `M`
//!    at every base `b >= n0 + K - R`; `H` agrees with `Q` on the `K`
//!    identifying rows, and the leader of `M` is a nonzero polynomial at every
//!    target `>= n0 + K` (decided exactly), so `H_n = Q_n` for every
//!    `n >= n0` by induction.
//!
//! [`GammaLiftCertificate::verify`] replays every step from the recorded
//! input, parent, schedule, scale, shift and gauge.  A nonzero remainder
//! leaves the identity supported by the compared rows only.
//!
//! # Detection
//!
//! [`find_gamma_lifts`] computes the gamma rows of each parent of a menu
//! ([`GammaLiftParent::standard_menu`]: Eulerian, type-B Eulerian, Narayana;
//! callers add others, such as
//! [`GammaLiftParent::one_over_k_eulerian_symmetric_part`]) and, for each
//! index shift, reads the scale `c` from two consecutive nonzero coefficients
//! and the factors `lambda(n)` from the rows; `lambda` must be constant or
//! have a constant ratio `lambda(n)/lambda(n-1)`.  Every match is certified
//! with [`certify_gamma_lift`].
//!
//! The group is split into `certify` (candidate, certificate, certification)
//! and `search` (parent menu and detection).

mod certify;
mod search;

pub use certify::{
    certify_gamma_lift, GammaDivisionQuotient, GammaLiftCandidate, GammaLiftCertificate,
    GammaLiftCertification, GammaLiftIdentityStatus, GammaLiftOptions, GammaLiftParent,
    GammaParentSymmetryStatus, GAMMA_LIFT_CERTIFICATION_SCHEMA, MAX_GAMMA_LIFT_CHECKED_ROWS,
    MAX_GAMMA_LIFT_INDEX_SHIFT, MAX_GAMMA_LIFT_OPERATOR_TERMS,
};
pub use search::{
    find_gamma_lifts, GammaLiftRejection, GammaLiftReport, GammaLiftSearchOptions,
    GAMMA_LIFTS_SCHEMA,
};
