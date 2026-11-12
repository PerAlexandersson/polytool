//! Provers for the parent statement of a Veronese-section lift: every parent
//! row `P_n` (from a given index on) is zero or real-rooted with all
//! coefficients of one sign, equivalently zero or with every root real and
//! `<= 0`.  By the transfer theorem (see [`find_veronese_section_lifts`])
//! every section of such a row is zero or real-rooted with nonpositive roots.
//!
//! [`prove_parent_real_rootedness`] tries, in order:
//!
//! 1. **periodic product**: the rule is `P_(n+k) = c(x) P_n` with `c` free of
//!    `n` (a nonzero constant denominator is allowed), and `c` and the `k`
//!    seeds are nonzero, real-rooted and of one sign; then
//!    `P_(n0 + k m + i) = c^m P_(n0 + i)` (shared with the v1 search);
//! 2. **paper theorem**: [`match_recurrence_theorems`] with boundary status
//!    `canonical` or `canonical_after_relabel` and, for the predecessor
//!    theorem, an all-index certificate scope; both theorems conclude that
//!    every row has nonnegative coefficients and is real-rooted;
//! 3. **root domain**: [`transform_recurrence`] with a root-domain request
//!    over the maps identity (`affine 1, 0`), `f_polynomial`, `reversal` (the
//!    last two with the degree schedule `d_n` fitted to the first rows, which
//!    the transform certifies or reports as a condition) and the inverse
//!    `f`-map `mobius (1, 0, -1, 1)`.  Only a claim about the source rows
//!    with property `real_rooted_in` an index-independent interval whose
//!    upper end is `<= 0`, evidence `completed_recurrence_induction` or
//!    `certified_family_invariant`, no `conditional_on` entries, no last
//!    index and a first index `<=` the parent's first index is accepted.
//!    When no map gives such a claim, the menu is repeated with the parent's
//!    first index moved `1, 2, 3` rows later (the moved seeds are generated
//!    parent rows); a claim from the moved first index on suffices, because
//!    the skipped rows are among the rows checked exactly below (each zero or
//!    real-rooted with coefficients of one sign, which is the parent
//!    statement for those rows);
//! 4. **three-term recurrence**: [`check_three_term_recurrence_theorem`] with
//!    the root region `(-inf, 0)` (orthogonal-polynomial recurrences, also
//!    after the coordinate change `y = 1/(x + v_0)`, as for the Narayana
//!    polynomials), proved only when every hypothesis is verified;
//! 5. **Liu--Wang**: [`check_liu_wang_theorem`] (first-order rules
//!    `P_n = a_n(x) P_(n-1) + b_n(x) P'_(n-1)` with every root kept in
//!    `(-inf, 0)`, as for the type-B Eulerian polynomials), proved only when
//!    every hypothesis is verified;
//! 6. **Liu--Wang derivative lag**: [`check_liu_wang_derivative_lag_theorem`]
//!    (second-order rules `P_n = A_n(x) P_(n-1) + V_n(x) P'_(n-1) +
//!    W_n(x) P_(n-2)`, the RealRooted routes `halfGrowth` and `sign`, as for
//!    the Delannoy and deco rows; also after dividing every row by a power
//!    of `x`, as for the derangement excedance rows), proved only when every
//!    hypothesis is verified;
//! 7. **Ma--Wang interval**: [`check_ma_wang_interval_theorem`] (first-order
//!    rules `P_n = (alpha_n x + beta_n) P_(n-1) + kappa_n (x - L)(x - U)
//!    P'_(n-1)` with every root kept in `[L, U]` by the exact two-endpoint
//!    criterion, as for the rows `k! S(n, k)` of the ordered set partitions),
//!    proved when every hypothesis is verified and the certified interval
//!    lies in `(-inf, 0]` (`U <= 0`);
//!
//!    Route 7b, tried after route 7, is **Liu--Wang three-term**:
//!    [`check_liu_wang_three_term_theorem`] (lag rules `P_n = A_n(x) P_(n-1)
//!    [+ V_n(x) P'_(n-1)] + B_n(x) P_(n-2)` with degree steps `0` or `1` and
//!    `B_n <= 0` on the real line, on `(-inf, 0]` or on a root window, as for
//!    the Lucas and Fibonacci-type rows), proved when every hypothesis is
//!    verified and the theorem certifies a root region in `(-inf, 0]`;
//! 8. **iterated preserver**: [`check_iterated_preserver_theorem`]
//!    (first-order or single-lag rules `P_n = a_n(x) S_n(P_(n-L))` with a
//!    real-rooted multiplier after the identity or an operator of the
//!    preserver menu, Brändén's Weyl step `A + B D + C D^2` included,
//!    applied from the last seed rows on, every earlier seed row checked
//!    exactly, as for the 3-Pascal rows `P_n = (1 + x) P_(n-1)` after
//!    `P_2 = 1 + 3x + x^2`), proved when every hypothesis is verified and the
//!    certified region is `(-inf, 0]`; the attempt is recorded only when
//!    the rule has the theorem's form (not `not_applicable`);
//!
//!    Route 8b, tried after route 8, is **hypergeometric row**:
//!    [`check_hypergeometric_row_theorem`] (rows whose coefficient ratio is
//!    a product of factors linear in `k`, with a row relation proved for
//!    every `n` by an invariant row-annihilator system, written as the
//!    image of `(1 + x)^M(n)` under a product of nonnegative multiplier
//!    sequences, as for OEIS A330798), proved when every hypothesis is
//!    verified (its region is `(-inf, 0)` or `(-inf, 0]`); the attempt is
//!    recorded only when the rule has the theorem's form (a fitted relation
//!    of hypergeometric shape), as for route 8;
//!
//!    Routes 8c and 8d, tried after route 8b and only for a parent that the
//!    earlier routes leave open, are the RealRooted row routes
//!    **fixed-delay window** ([`check_variable_delay_theorem`], route
//!    `variableDelay`: `P_n = a(n) P_(n-1) + b(n) x P_(n-L)`, `L >= 2`,
//!    `a > 0`, `b >= 0`, every root in `(-inf, 0)`) and **residue classes**
//!    ([`check_residue_class_subsequence_theorem`], route `subseq`: every
//!    lag a multiple of `p >= 2`, each class `P_(n0 + j + p m)` proved by a
//!    theorem above), proved when every hypothesis is verified and the
//!    certified region lies in `(-inf, 0]` (route 8d is run with the
//!    requested region `(-inf, 0)` and, failing that, without one, since
//!    its classes may have roots at `0`); each attempt is recorded only when
//!    the rule has the route's form.  Route 8d also accepts a proof with
//!    zero residue classes
//!    ([`RealRootednessTheoremConclusion::zero_residue_classes`]: every row
//!    of a class is zero, infinitely many zero rows, as for
//!    `P_n = x^2 P_(n-2)` with `P_1 = 0`), because the parent statement
//!    allows zero rows and every transfer theorem of the lifts maps a zero
//!    parent row to a zero input row; the attempt's detail names the zero
//!    classes, and the route claims no interlacing;
//!
//!    Route 8e, tried after route 8d and only for a parent still open, is the
//!    **index-section lift** ([`check_index_section_lift_theorem`], route
//!    `pathPower`: `x`-only multipliers, the rows, possibly reversed, are
//!    every `p`-th row of a constant-coefficient parent proved by a theorem
//!    above), proved when every hypothesis is verified and the certified
//!    region lies in `(-inf, 0]` (run with the requested region `(-inf, 0)`
//!    and, failing that, without one, since a reversed section may have
//!    roots at `0`); recorded only when a section identity is found;
//!
//!    Route 8f, tried after route 8e and only for a parent still open, is
//!    the **factor gauge** ([`check_factor_gauge_theorem`], routes
//!    `powerForm`, `closedForm` and `factorPower`: every lag-`j` multiplier
//!    is `kappa_j(n) q(x)^j` for one fixed `q`, so `P_n = q^(n-e) g S_n`
//!    with a bounded residual `S_n` decided for every index), proved when
//!    every hypothesis is verified and the certified region lies in
//!    `(-inf, 0]` (run with the requested region `(-inf, 0)` and, failing
//!    that, without one, since `q` or `g` may have roots at `0`); recorded
//!    only when the rule has the theorem's form;
//!    Route 8g, tried after the existing routes and only for a parent still
//!    open, is **positive lag-two Euler** ([`check_lag_two_euler_theorem`],
//!    route `lagTwoEuler`): the exact Euler operator form with `a > 0`,
//!    `b >= a`, `c >= 0` and positive seeds gives simple negative roots and
//!    strict adjacent and auxiliary interlacing. Recorded only when the
//!    operator form matches; also offered in both sub-family menus.
//!    Route 8h, while still open, is the **binomial-perturbed Euler orbit**
//!    ([`check_binomial_perturbed_euler_orbit_theorem`], RealRooted issue
//!    #2513): exact theta operators, canonical seeds, `s>=1`, `kappa>=s+4`,
//!    also after exact Eulerian Mobius transport. Gives simple negative
//!    roots and strict adjacent interlacing; recorded only when applicable.
//! 9. **attestation**: a caller-supplied [`VeroneseParentAttestation`] whose
//!    record equals the parent exactly (same first index, the same relation
//!    term by term, the same identifying seeds); the verdict is then
//!    conditional on the attestation, never proved;
//! 10. **catalog flag**: a bundled OEIS catalog entry flagged
//!     [`OeisRealRootednessStatus::LeanProved`] (every row zero or
//!     real-rooted by a Lean proof made outside polytool) whose rule, first
//!     index and initial rows equal the parent exactly, matched as in route 9
//!     with the attestation label `oeis_catalog:<A-number>:lean_proved`.  The
//!     flag gives real-rootedness only, so the route also needs
//!     coefficients of one sign for every row, which polytool proves here
//!     (see "One-signed rows" below); the verdict is conditional on the flag,
//!     never proved.
//!
//! Routes 4 to 7, 7b, 8 and 8b run whether or not the root-domain route is enabled.
//! Routes 4 and 5 conclude that every parent row from the first index on is
//! real-rooted with every root in `(-inf, 0)`; route 6 that every parent row
//! is zero (only seed rows) or real-rooted with every root in `(-inf, 0]`
//! (`(-inf, 0)` without the power of `x`), so it runs without a requested
//! root region; route 7 that every parent row is real-rooted with every root
//! in `[L, U]`, `U <= 0`; route 7b that every parent row is real-rooted with
//! every root in `(-inf, 0)` or `(-inf, 0]`; route 8 that every parent row is
//! zero (only seed rows) or real-rooted with every root in `(-inf, 0]`
//! (route 8d: zero rows also in whole residue classes).  All six
//! move their induction start up to three rows past
//! the parent seeds when the theorem fails at the seeds (the default
//! [`RealRootednessTheoremOptions::maximum_start_delay`]; the skipped rows
//! are checked exactly by the theorem), and the attempt names the delay
//! used.
//!
//! The theorem routes (product induction and routes 4 to 8f) come from the
//! registry of `parent_theorems`, in its order; each route's `ParentPolicy`
//! gives its stage (product induction before route 2; routes 4 to 8b while
//! the parent is unproved or proved only by product induction; routes 8c to
//! 8f while it is unproved), its requested region (`(-inf, 0)`, none, or
//! `(-inf, 0)` and then none) and whether an attempt that does not apply is
//! recorded.  Every theorem route proves the parent only when its certified
//! region lies in `(-inf, 0]`, and its attempt's variant is the route name.
//!
//! A few parent rows are also checked exactly (nonzero rows of one sign and
//! real-rooted).  A failing row refutes the parent statement; no route can
//! then succeed and an attestation is not applied.
//!
//! # Interlacing
//!
//! A proved parent records the interlacing of consecutive rows its route
//! concludes ([`VeroneseParentProof::interlacing`]), always oriented
//! `P_(m-1) << P_m` (`P_m` carries the largest zero, as in
//! [`RowInterlacing`]):
//!
//! - route 2: both window theorems state `P_n << P_(n+s)` for every
//!   canonical index `n >= 0` and `0 <= s <= max_step` (`max_step >= 1`, or
//!   every step for the constant family); the parent rows are the canonical
//!   rows (relabelled), so `P_(m-1) << P_m` weakly for every
//!   `m >= first_index + 1`;
//! - routes 4 to 7, 7b and 8: the theorem report's own conclusion
//!   ([`RealRootednessTheoremConclusion::interlacing`] from
//!   `interlacing_from` on, strict or weak as the theorem states);
//! - routes 1 and 3 conclude no interlacing, and an attested parent records
//!   none.
//!
//! # Replaying recorded proofs
//!
//! The provers grow stronger over time (a new route can prove a parent that
//! was open or only attested before).  A recorded proof is therefore
//! replayed by [`check_recorded_parent_proof`]: it must claim nothing beyond
//! the current result ([`VeroneseParentProof::is_supported_by`]), not equal
//! it.  A weaker record (`not_proved`, or conditional on an attestation
//! where the statement is now proved outright) stays valid; a record that
//! claims more than the current provers give is rejected.  A record proved
//! by a route `R` that the current cascade no longer reaches (an earlier
//! route now settles the parent first, as the product induction does for
//! records of the paper-theorem and root-domain routes) is replayed by `R`
//! alone, with every other route skipped: it is accepted exactly when that
//! replay proves the parent by `R` with the recorded statement and
//! interlacing, so a record whose own route fails is never accepted.
//!
//! # Left-hand factors `x^e f(n)`
//!
//! A parent rule with a left-hand factor `x^e f(n)` (`e >= 1`, `f` free of
//! `x`, as the adaptive search may return it) whose `x^e` divides every
//! right-hand coefficient is replaced, for routes 1 to 8 and the one-signed
//! test of route 10, by the exactly equivalent rule
//! `f(n) P_n = sum (c / x^e) D^d P_(n-r)` (same rows; see
//! [`LeftHandFactorNormalization`]), recorded in
//! [`VeroneseParentProof::left_hand_factor_normalization`].  An attestation
//! (routes 8 and 9) applies when its record equals the parent rule as
//! supplied or the normalized rule (both define the same rows).  When
//! `x^e` does not divide the right-hand side, the rule is used as it is and
//! the theorem routes report it as not applicable with that reason.
//!
//! # One-signed rows (route 10)
//!
//! Write the rule as `f(n) P_n = sum_terms c(n,x) D^d P_(n-r) + g(n,x)` for
//! `n >= N1` (the first computed index, after the `K` seeds).  Every row has
//! nonnegative coefficients (after one global sign change of a homogeneous
//! rule) when:
//!
//! - every seed has nonnegative coefficients;
//! - `f` does not depend on `x` and `f(n) > 0` for every integer `n >= N1`,
//!   no term carries the factor `(-1)^n`, and every coefficient of `g` (a
//!   polynomial in `n`) is `>= 0` for `n >= N1` (exact sign checker);
//! - for every offset `r` and shift `e`, the combined weight
//!   `W_(r,e)(n,s) = sum over the terms at offset r and the monomials c_j(n) x^j
//!   with j - d = e of c_j(n) s (s-1) ... (s-d+1)`, which multiplies the
//!   coefficient of `x^s` in `P_(n-r)` in the coefficient of `x^(s+e)` in
//!   `f(n) P_n`, is `>= 0` for every integer `n >= N1` and every
//!   `0 <= s <= deg P_(n-r)`.
//!
//! Then `f(n) P_n` is a sum of nonnegative multiples of nonnegative
//! coefficients, and induction on `n` gives the claim.  The weight condition
//! is shown by one of two sufficient tests: the polynomial
//! `W(N1 + v, s)` or `W(N1 + v, D(N1 + v - r) - u)` has nonnegative
//! coefficients in the basis `v(v-1)...(v-i+1) w(w-1)...(w-j+1)` (`w = s` or
//! `u`), which is nonnegative at all integers `v, w >= 0`.  Here
//! `D(m) = alpha + beta m` is an upper bound for `deg P_m`, `m >= ` the first
//! index: `beta >= max(0, (deg_x c - d) / r)` over the terms and `alpha`
//! covers the seeds and `g`, so that `deg P_n <= max(deg_x c - d + D(n - r))
//! <= D(n)` by induction.  Every test is exact; a failure only means the
//! route does not apply.

use super::parent_theorems::rule_with_x_free_left_hand_factor;
use super::veronese_algebra::classify_row;
use super::*;
use crate::integer_range_sign::{decide_index_polynomial_sign, IndexSignCondition, IntegerRange};
use crate::oeis::{entries_with_real_rootedness_flag, OeisRealRootednessStatus};
use crate::real_rootedness::RealInterval;
use crate::renewal::{
    match_recurrence_theorems, CertificateScope, RecurrenceTheoremMatch, TheoremBoundaryStatus,
};
use crate::vec_poly;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};

type Q = BigRational;

/// Number of parent rows checked exactly for being zero or real-rooted with
/// coefficients of one sign (evidence only; a failure refutes the parent).
const PARENT_ROWS_CHECKED: usize = 10;

/// Rows generated for the root-domain transform of a parent.
const ROOT_DOMAIN_ROWS: usize = 8;

/// Largest number of rows by which the root-domain route moves the parent's
/// first index (smaller than [`PARENT_ROWS_CHECKED`], so that every skipped
/// row is among the rows checked exactly).
const MAXIMUM_ROOT_DOMAIN_START_DELAY: usize = 3;

/// A caller-supplied statement that every parent row `P_n`, `n >= first_index`,
/// is zero or real-rooted with coefficients of one sign (for example a
/// published theorem).  It is applied only to a parent whose record equals
/// `recurrence` exactly; the verdict is then conditional on it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VeroneseParentAttestation {
    /// Short label, repeated in the verdict.
    pub label: String,
    /// The attested parent: rule, first index and seeds.
    pub recurrence: RecurrenceJson,
    /// First parent index covered by the statement; it must be at most the
    /// first parent index the lift uses.
    pub first_index: usize,
    /// The attested statement, in words.
    pub statement: String,
    /// Where the statement comes from (reference, theorem, URL).
    pub source: String,
}

/// Outcome of the parent provers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum VeroneseParentProofStatus {
    /// Proved by a periodic product, a paper theorem, the root-domain route,
    /// the three-term recurrence theorem, the Liu--Wang theorem, the
    /// Liu--Wang derivative-lag theorem, the Ma--Wang interval theorem or the
    /// iterated-preserver theorem.
    Proved,
    /// Not proved; a matching caller-supplied attestation asserts it.
    ConditionalOnAttestation,
    /// No route succeeded.
    NotProved,
    /// A checked parent row is not zero and not real-rooted with one-signed
    /// coefficients: the parent statement is false.
    Refuted,
}

/// The route that settled the parent statement.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum VeroneseParentProofRoute {
    /// `P_(n+k) = c(x) P_n` with `c` and the seeds real-rooted of one sign.
    PeriodicProduct,
    /// General product induction, with all roots strictly negative.
    ProductRecurrence,
    /// A matched paper theorem (`match_recurrence_theorems`).
    PaperTheorem,
    /// A root-domain claim of `transform_recurrence`.
    RootDomain,
    /// [`check_three_term_recurrence_theorem`] with roots in `(-inf, 0)`.
    ThreeTermRecurrence,
    /// [`check_liu_wang_theorem`].
    LiuWang,
    /// [`check_liu_wang_derivative_lag_theorem`].
    LiuWangDerivativeLag,
    /// [`check_ma_wang_interval_theorem`] with an interval `[L, U]`,
    /// `U <= 0`.
    MaWangInterval,
    /// [`check_liu_wang_three_term_theorem`] with a root region in
    /// `(-inf, 0]`.
    LiuWangThreeTerm,
    /// [`check_iterated_preserver_theorem`] with the region `(-inf, 0]`.
    IteratedPreserver,
    /// [`check_hypergeometric_row_theorem`] with a region in `(-inf, 0]`.
    HypergeometricRow,
    /// [`check_variable_delay_theorem`] (region `(-inf, 0)`).
    VariableDelay,
    /// [`check_residue_class_subsequence_theorem`] with a region in
    /// `(-inf, 0]`.
    ResidueClassSubsequences,
    /// [`check_index_section_lift_theorem`] with a region in `(-inf, 0]`.
    IndexSectionLift,
    /// [`check_factor_gauge_theorem`] with a region in `(-inf, 0]`.
    FactorGauge,
    /// [`check_lag_two_euler_theorem`] with simple negative roots and strict interlacing.
    LagTwoEuler,
    /// [`check_binomial_perturbed_euler_orbit_theorem`] (RealRooted issue #2513).
    BinomialPerturbedEulerOrbit,
    /// [`check_rational_pf_riordan_theorem`] with nonpositive roots.
    RationalPfRiordan,
    /// Shared degree-law reversal followed by a theorem of the parent menu.
    DegreeLawReversal,
    /// A caller-supplied attestation.
    Attestation,
}

impl VeroneseParentProofRoute {
    /// The snake-case name, as serialized (for example `"liu_wang"`).
    pub const fn as_str(self) -> &'static str {
        match self {
            VeroneseParentProofRoute::PeriodicProduct => "periodic_product",
            VeroneseParentProofRoute::PaperTheorem => "paper_theorem",
            VeroneseParentProofRoute::RootDomain => "root_domain",
            VeroneseParentProofRoute::Attestation => "attestation",
            VeroneseParentProofRoute::ProductRecurrence
            | VeroneseParentProofRoute::ThreeTermRecurrence
            | VeroneseParentProofRoute::LiuWang
            | VeroneseParentProofRoute::LiuWangDerivativeLag
            | VeroneseParentProofRoute::MaWangInterval
            | VeroneseParentProofRoute::LiuWangThreeTerm
            | VeroneseParentProofRoute::IteratedPreserver
            | VeroneseParentProofRoute::HypergeometricRow
            | VeroneseParentProofRoute::VariableDelay
            | VeroneseParentProofRoute::ResidueClassSubsequences
            | VeroneseParentProofRoute::IndexSectionLift
            | VeroneseParentProofRoute::LagTwoEuler
            | VeroneseParentProofRoute::RationalPfRiordan
            | VeroneseParentProofRoute::DegreeLawReversal
            | VeroneseParentProofRoute::BinomialPerturbedEulerOrbit
            | VeroneseParentProofRoute::FactorGauge => match self.theorem() {
                Some(theorem) => theorem.as_str(),
                None => unreachable!(),
            },
        }
    }

    /// The [`RealRootednessTheorem`] behind a theorem route (the route and
    /// the theorem have the same name); `None` for the periodic-product,
    /// paper-theorem, root-domain and attestation routes.
    pub const fn theorem(self) -> Option<RealRootednessTheorem> {
        match self {
            VeroneseParentProofRoute::ProductRecurrence => {
                Some(RealRootednessTheorem::ProductRecurrence)
            }
            VeroneseParentProofRoute::ThreeTermRecurrence => {
                Some(RealRootednessTheorem::ThreeTermRecurrence)
            }
            VeroneseParentProofRoute::LiuWang => Some(RealRootednessTheorem::LiuWang),
            VeroneseParentProofRoute::LiuWangDerivativeLag => {
                Some(RealRootednessTheorem::LiuWangDerivativeLag)
            }
            VeroneseParentProofRoute::MaWangInterval => Some(RealRootednessTheorem::MaWangInterval),
            VeroneseParentProofRoute::LiuWangThreeTerm => {
                Some(RealRootednessTheorem::LiuWangThreeTerm)
            }
            VeroneseParentProofRoute::IteratedPreserver => {
                Some(RealRootednessTheorem::IteratedPreserver)
            }
            VeroneseParentProofRoute::HypergeometricRow => {
                Some(RealRootednessTheorem::HypergeometricRow)
            }
            VeroneseParentProofRoute::VariableDelay => Some(RealRootednessTheorem::VariableDelay),
            VeroneseParentProofRoute::LagTwoEuler => Some(RealRootednessTheorem::LagTwoEuler),
            VeroneseParentProofRoute::BinomialPerturbedEulerOrbit => {
                Some(RealRootednessTheorem::BinomialPerturbedEulerOrbit)
            }
            VeroneseParentProofRoute::ResidueClassSubsequences => {
                Some(RealRootednessTheorem::ResidueClassSubsequences)
            }
            VeroneseParentProofRoute::IndexSectionLift => {
                Some(RealRootednessTheorem::IndexSectionLift)
            }
            VeroneseParentProofRoute::DegreeLawReversal => {
                Some(RealRootednessTheorem::DegreeLawReversal)
            }
            VeroneseParentProofRoute::RationalPfRiordan => {
                Some(RealRootednessTheorem::RationalPfRiordan)
            }
            VeroneseParentProofRoute::FactorGauge => Some(RealRootednessTheorem::FactorGauge),
            VeroneseParentProofRoute::PeriodicProduct
            | VeroneseParentProofRoute::PaperTheorem
            | VeroneseParentProofRoute::RootDomain
            | VeroneseParentProofRoute::Attestation => None,
        }
    }

    /// The parent route of a theorem (the inverse of [`Self::theorem`]).
    pub(in crate::recurrence) const fn of_theorem(theorem: RealRootednessTheorem) -> Self {
        match theorem {
            RealRootednessTheorem::ProductRecurrence => VeroneseParentProofRoute::ProductRecurrence,
            RealRootednessTheorem::ThreeTermRecurrence => {
                VeroneseParentProofRoute::ThreeTermRecurrence
            }
            RealRootednessTheorem::LiuWang => VeroneseParentProofRoute::LiuWang,
            RealRootednessTheorem::LiuWangDerivativeLag => {
                VeroneseParentProofRoute::LiuWangDerivativeLag
            }
            RealRootednessTheorem::MaWangInterval => VeroneseParentProofRoute::MaWangInterval,
            RealRootednessTheorem::LiuWangThreeTerm => VeroneseParentProofRoute::LiuWangThreeTerm,
            RealRootednessTheorem::IteratedPreserver => VeroneseParentProofRoute::IteratedPreserver,
            RealRootednessTheorem::HypergeometricRow => VeroneseParentProofRoute::HypergeometricRow,
            RealRootednessTheorem::VariableDelay => VeroneseParentProofRoute::VariableDelay,
            RealRootednessTheorem::LagTwoEuler => VeroneseParentProofRoute::LagTwoEuler,
            RealRootednessTheorem::BinomialPerturbedEulerOrbit => {
                VeroneseParentProofRoute::BinomialPerturbedEulerOrbit
            }
            RealRootednessTheorem::ResidueClassSubsequences => {
                VeroneseParentProofRoute::ResidueClassSubsequences
            }
            RealRootednessTheorem::IndexSectionLift => VeroneseParentProofRoute::IndexSectionLift,
            RealRootednessTheorem::RationalPfRiordan => VeroneseParentProofRoute::RationalPfRiordan,
            RealRootednessTheorem::DegreeLawReversal => VeroneseParentProofRoute::DegreeLawReversal,
            RealRootednessTheorem::FactorGauge => VeroneseParentProofRoute::FactorGauge,
        }
    }
}

/// One tried route.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VeroneseParentProofAttempt {
    /// The route.
    pub route: VeroneseParentProofRoute,
    /// Theorem identifier, transformation kind or attestation label.
    pub variant: String,
    /// Whether it settled the parent statement.
    pub success: bool,
    /// What was checked, or why it did not apply.
    pub detail: String,
}

/// The parent statement and how it was settled.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VeroneseParentProof {
    /// The outcome.
    pub status: VeroneseParentProofStatus,
    /// The successful route, if any.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub route: Option<VeroneseParentProofRoute>,
    /// The parent statement as proved, attested or refuted, or why it is open.
    pub statement: String,
    /// Attestation labels the status depends on (empty unless conditional).
    pub conditional_on: Vec<String>,
    /// Every route tried, in order.
    pub attempts: Vec<VeroneseParentProofAttempt>,
    /// Absolute index of the first checked parent row that is nonzero and not
    /// real-rooted with coefficients of one sign, if any.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub failing_row_index: Option<usize>,
    /// Number of parent rows checked exactly (evidence only).
    pub rows_checked: usize,
    /// The interlacing of consecutive parent rows concluded by the proving
    /// route (module documentation, "Interlacing"); `None` when the status
    /// is not proved or the route concludes no interlacing.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub interlacing: Option<LiftInterlacing>,
    /// The exact normalization of a parent left-hand factor `x^e f(n)`
    /// (module documentation, "Left-hand factors `x^e f(n)`"): routes 1 to
    /// 8 and the one-signed test of route 10 ran on the normalized rule.
    /// Omitted when the parent rule was used as it is.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub left_hand_factor_normalization: Option<LeftHandFactorNormalization>,
}

impl VeroneseParentProof {
    /// Whether this (recorded) parent proof claims nothing beyond
    /// `recomputed`, the parent provers' current result for the same parent
    /// (module documentation, "Replaying recorded proofs").
    ///
    /// The row facts (`failing_row_index`, `rows_checked`) and the
    /// left-hand factor normalization must agree.  Then a record that equals
    /// the recomputation is supported, and otherwise:
    ///
    /// - `not_proved` (no route, no labels, no interlacing) is supported
    ///   unless the recomputation is `refuted`;
    /// - `conditional_on_attestation` (route `attestation`, no interlacing)
    ///   is supported when the recomputation proves the statement outright
    ///   (the conditional statement is then implied) or is conditional on
    ///   the same labels;
    /// - `proved` is supported when the recomputation is proved by the same
    ///   route with the same statement and interlacing (attempts that failed
    ///   before it may differ);
    /// - `refuted` is supported when the recomputation is refuted with the
    ///   same statement.
    ///
    /// A record claiming more than the recomputation (proved where it is
    /// open, an interlacing it does not give, other labels) is not
    /// supported.
    pub fn is_supported_by(&self, recomputed: &VeroneseParentProof) -> bool {
        if self == recomputed {
            return true;
        }
        if self.failing_row_index != recomputed.failing_row_index
            || self.rows_checked != recomputed.rows_checked
            || self.left_hand_factor_normalization != recomputed.left_hand_factor_normalization
        {
            return false;
        }
        let proved_outright = recomputed.status == VeroneseParentProofStatus::Proved
            && recomputed
                .route
                .is_some_and(|route| route != VeroneseParentProofRoute::Attestation);
        match self.status {
            VeroneseParentProofStatus::NotProved => {
                self.route.is_none()
                    && self.conditional_on.is_empty()
                    && self.interlacing.is_none()
                    && recomputed.status != VeroneseParentProofStatus::Refuted
            }
            VeroneseParentProofStatus::ConditionalOnAttestation => {
                self.route == Some(VeroneseParentProofRoute::Attestation)
                    && !self.conditional_on.is_empty()
                    && self.interlacing.is_none()
                    && (proved_outright
                        || (recomputed.status
                            == VeroneseParentProofStatus::ConditionalOnAttestation
                            && recomputed.conditional_on == self.conditional_on))
            }
            VeroneseParentProofStatus::Proved => {
                proved_outright
                    && self.route == recomputed.route
                    && self.statement == recomputed.statement
                    && self.interlacing == recomputed.interlacing
                    && self.conditional_on.is_empty()
            }
            VeroneseParentProofStatus::Refuted => {
                recomputed.status == VeroneseParentProofStatus::Refuted
                    && self.statement == recomputed.statement
            }
        }
    }
}

/// Replay a recorded parent proof: run the parent provers on `record` with
/// the recorded options and check that the record claims nothing beyond the
/// result ([`VeroneseParentProof::is_supported_by`]).  A record made with
/// an earlier, weaker version of the provers keeps replaying.  A record
/// proved by a polytool route `R` that the current cascade no longer reaches
/// first is replayed by `R` alone (every other route skipped) and accepted
/// exactly when that replay proves the parent by `R` with the recorded
/// statement and interlacing.
///
/// # Errors
///
/// [`RecurrenceError::VerificationFailed`] when the recorded proof is not
/// supported, and the errors of [`prove_parent_real_rootedness`].
pub fn check_recorded_parent_proof(
    record: &ParentRecurrenceRecord<'_>,
    needed_from: usize,
    attestations: &[VeroneseParentAttestation],
    try_root_domain: bool,
    recorded: &VeroneseParentProof,
) -> Result<(), RecurrenceError> {
    let recomputed =
        prove_parent_real_rootedness(record, needed_from, attestations, try_root_domain)?;
    if recorded.is_supported_by(&recomputed)
        || recorded_route_replays(record, needed_from, recorded)?
    {
        Ok(())
    } else {
        Err(RecurrenceError::VerificationFailed(
            "the recorded parent proof claims more than the parent provers give".into(),
        ))
    }
}

/// Whether a recorded proof `proved` by a polytool route `R` is reproduced
/// by replaying `R` alone on the same record (module documentation,
/// "Replaying recorded proofs"): the replay must be proved by `R` with the
/// recorded statement and interlacing
/// ([`VeroneseParentProof::is_supported_by`] of the replay).  `false` for
/// any other recorded status or an attested proof.
pub(in crate::recurrence) fn recorded_route_replays(
    record: &ParentRecurrenceRecord<'_>,
    needed_from: usize,
    recorded: &VeroneseParentProof,
) -> Result<bool, RecurrenceError> {
    let route = match (recorded.status, recorded.route) {
        (VeroneseParentProofStatus::Proved, Some(route))
            if route != VeroneseParentProofRoute::Attestation =>
        {
            route
        }
        _ => return Ok(false),
    };
    let replay = prove_parent_by_routes(record, needed_from, &[route])?;
    Ok(replay.route == Some(route) && recorded.is_supported_by(&replay))
}

/// A parent rule with its identifying seeds and first index: the input of
/// [`prove_parent_real_rootedness`].
pub struct ParentRecurrenceRecord<'a> {
    /// The parent rule.
    pub rule: &'a Recurrence,
    /// The `K` identifying seeds (`Recurrence::generation_initial_count`);
    /// `seeds[0]` is row `first_index`.
    pub seeds: &'a [Vec<Q>],
    /// The index of the first seed.
    pub first_index: usize,
}

/// Route 1: `P_(n+k) = c(x) P_n` with `c` and the `k` seeds nonzero,
/// real-rooted and of one sign.  Returns the reason when it does not apply.
pub(super) fn periodic_product_parent(rule: &Recurrence, seeds: &[Vec<Q>]) -> Result<(), String> {
    let (relation, k) = ImplicitRelation::from_recurrence(rule).map_err(|e| e.to_string())?;
    let op = relation.op();
    let leader =
        op.0.get(&(k, 0))
            .and_then(|c| c.as_constant())
            .filter(|c| !c.is_zero())
            .ok_or("the leader is not a nonzero constant")?;
    let mut factor = None;
    for ((shift, order), coefficient) in &op.0 {
        if *shift == k {
            continue;
        }
        if *shift != 0 || *order != 0 || !coefficient.is_index_free() {
            return Err("the rule is not P_(n+k) = c(x) P_n with c free of n".into());
        }
        factor = coefficient.rows().first().cloned();
    }
    let c: Vec<Q> = factor
        .unwrap_or_default()
        .iter()
        .map(|value| -value / &leader)
        .collect();
    if !classify_row(&vec_poly::trim(&c)).has_nonpositive_real_roots() {
        return Err("the factor c(x) is not nonzero, real-rooted and of one sign".into());
    }
    if seeds.len() != k {
        return Err(format!(
            "the rule needs {k} seeds, {} identify it",
            seeds.len()
        ));
    }
    if let Some(position) = seeds
        .iter()
        .position(|seed| !classify_row(seed).has_nonpositive_real_roots())
    {
        return Err(format!(
            "seed {position} is not nonzero, real-rooted and of one sign"
        ));
    }
    Ok(())
}

/// Route 2, with the interlacing of the first successful match (module
/// documentation, "Interlacing").
fn theorem_route(
    record: &ParentRecurrenceRecord<'_>,
    rows: &[Vec<Q>],
) -> (Vec<VeroneseParentProofAttempt>, Option<LiftInterlacing>) {
    let matches = match_recurrence_theorems(record.rule, record.first_index, rows);
    let mut attempts = Vec::new();
    let mut interlacing = None;
    for theorem_match in matches.theorem_matches {
        let boundary_ok = matches!(
            theorem_match.boundary.status,
            TheoremBoundaryStatus::Canonical | TheoremBoundaryStatus::CanonicalAfterRelabel
        );
        let scope_ok = theorem_match
            .variable_predecessor
            .as_ref()
            .is_none_or(|certificate| certificate.scope == CertificateScope::AllIndices);
        let success = boundary_ok && scope_ok;
        if success && interlacing.is_none() {
            interlacing = paper_theorem_interlacing(&theorem_match, record.first_index);
        }
        let detail = if success {
            format!(
                "{} ({}): {}; the identifying seeds equal the canonical rows ({:?}), so every \
                 parent row is a canonical row",
                theorem_match.family,
                theorem_match.proof_status,
                theorem_match.conclusion,
                theorem_match.boundary.status
            )
        } else if !boundary_ok {
            format!(
                "the parent rows are not the canonical boundary rows ({:?})",
                theorem_match.boundary.status
            )
        } else {
            "the theorem certificate covers finitely many indices only".into()
        };
        attempts.push(VeroneseParentProofAttempt {
            route: VeroneseParentProofRoute::PaperTheorem,
            variant: theorem_match.theorem.clone(),
            success,
            detail,
        });
    }
    if attempts.is_empty() {
        let reasons: Vec<String> = matches
            .rejected_theorems
            .iter()
            .map(|r| format!("{}: {}", r.theorem, r.reason))
            .collect();
        attempts.push(VeroneseParentProofAttempt {
            route: VeroneseParentProofRoute::PaperTheorem,
            variant: "match_recurrence_theorems".into(),
            success: false,
            detail: format!("no theorem matched ({})", reasons.join("; ")),
        });
    }
    (attempts, interlacing)
}

/// The interlacing `P_(m-1) << P_m`, `m >= first_index + 1`, of a matched
/// window theorem whose window contains the step `1` (every canonical row
/// `P_n` interlaces `P_(n+s)` from the left for `0 <= s <= max_step`).
fn paper_theorem_interlacing(
    theorem_match: &RecurrenceTheoremMatch,
    first_index: usize,
) -> Option<LiftInterlacing> {
    if theorem_match.window.max_step.is_some_and(|step| step == 0) {
        return None;
    }
    let mut premises: Vec<String> = theorem_match
        .hypotheses
        .iter()
        .map(|h| format!("{} ({})", h.statement, h.check))
        .collect();
    premises.push(format!(
        "window: {} ({})",
        theorem_match.window.statement, theorem_match.window.inequality
    ));
    premises.push(format!(
        "the identifying parent seeds equal the canonical rows ({:?}), so every parent row is a \
         canonical row and consecutive parent rows are consecutive canonical rows",
        theorem_match.boundary.status
    ));
    Some(LiftInterlacing {
        kind: RowInterlacing::Weak,
        orientation: InterlacingOrientation::PreviousInterlacesNext,
        from_index: first_index.checked_add(1)?,
        justification: LiftInterlacingJustification {
            theorem: format!(
                "{} ({}; {}): {}",
                theorem_match.theorem,
                theorem_match.proof_status,
                theorem_match.reference,
                theorem_match.conclusion
            ),
            premises,
        },
    })
}

/// `d_n = d_(n0) + step (n - n0)` fitted to the degrees of the first rows,
/// when every observed row is nonzero and the degrees are affine.
fn fitted_schedule(first_index: usize, rows: &[Vec<Q>]) -> Option<DegreeSchedule> {
    let degrees: Vec<i64> = rows
        .iter()
        .map(|row| vec_poly::degree(row).and_then(|d| i64::try_from(d).ok()))
        .collect::<Option<_>>()?;
    let first = *degrees.first()?;
    let step = degrees.get(1).map_or(0, |second| second - first);
    let affine = degrees
        .iter()
        .enumerate()
        .all(|(i, d)| i64::try_from(i).is_ok_and(|i| *d == first + step * i));
    affine.then(|| DegreeSchedule::new(first_index, first as usize, step))
}

fn accepted_source_claim(claim: &RootDomainClaim, first_index: usize) -> bool {
    let evidence_ok = matches!(
        claim.evidence,
        RootEvidence::CompletedRecurrenceInduction | RootEvidence::CertifiedFamilyInvariant
    );
    let upper_ok = claim.interval.as_ref().is_some_and(|interval| {
        interval.is_index_independent()
            && interval.upper.as_ref().is_some_and(|upper| {
                parse_rational_coeff(&upper.value).is_ok_and(|value| value <= Q::zero())
            })
    });
    claim.subject == RootSubject::Source
        && claim.component == 0
        && claim.property == RootProperty::RealRootedIn
        && upper_ok
        && evidence_ok
        && claim.conditional_on.is_empty()
        && claim.last_index.is_none()
        && claim.first_index <= first_index
}

fn root_domain_route(
    record: &ParentRecurrenceRecord<'_>,
    rows: &[Vec<Q>],
) -> Vec<VeroneseParentProofAttempt> {
    let mut attempts = Vec::new();
    let seed_count = record.seeds.len();
    for delay in 0..=MAXIMUM_ROOT_DOMAIN_START_DELAY {
        // The moved record needs its seeds among the exactly checked rows.
        let Some(seeds) = rows.get(delay..delay + seed_count) else {
            break;
        };
        let moved = ParentRecurrenceRecord {
            rule: record.rule,
            seeds,
            first_index: record.first_index + delay,
        };
        let (moved_attempts, success) = root_domain_menu(&moved, &rows[delay..], delay);
        attempts.extend(moved_attempts);
        if success {
            break;
        }
    }
    attempts
}

/// The root-domain menu for `record` (moved `delay` rows past the parent's
/// own seeds; `rows` start at the record's first index).
fn root_domain_menu(
    record: &ParentRecurrenceRecord<'_>,
    rows: &[Vec<Q>],
    delay: usize,
) -> (Vec<VeroneseParentProofAttempt>, bool) {
    let source = RecurrenceJson::from_recurrence_rational(
        record.rule,
        record.first_index,
        record.seeds,
        None,
    );
    let schedule = fitted_schedule(record.first_index, rows);
    let mut menu = vec![(
        "identity (affine 1, 0)".to_string(),
        Some(RecurrenceTransformation::Affine {
            a: "1".into(),
            b: "0".into(),
        }),
    )];
    for (name, make) in [
        (
            "f_polynomial",
            (|s: DegreeSchedule| RecurrenceTransformation::FPolynomial { schedule: s })
                as fn(DegreeSchedule) -> RecurrenceTransformation,
        ),
        ("reversal", |s| RecurrenceTransformation::Reversal {
            schedule: s,
        }),
        ("mobius (1, 0, -1, 1)", |s| {
            RecurrenceTransformation::Mobius {
                a: "1".into(),
                b: "0".into(),
                c: "-1".into(),
                e: "1".into(),
                schedule: s,
            }
        }),
    ] {
        menu.push((name.to_string(), schedule.clone().map(make)));
    }
    let variant = |name: String| {
        if delay == 0 {
            name
        } else {
            format!("{name}, induction start moved {delay} row(s)")
        }
    };
    let moved_note = if delay == 0 {
        String::new()
    } else {
        format!(
            " (source moved {delay} row(s) past the parent seeds; the earlier parent rows are \
             among the rows checked exactly to be zero or real-rooted with coefficients of one \
             sign)"
        )
    };
    let mut attempts = Vec::new();
    for (name, transformation) in menu {
        let Some(transformation) = transformation else {
            attempts.push(VeroneseParentProofAttempt {
                route: VeroneseParentProofRoute::RootDomain,
                variant: variant(name),
                success: false,
                detail: "no affine degree schedule fits the first parent rows".into(),
            });
            continue;
        };
        let mut request = RecurrenceTransformRequest::from_source(
            transformation,
            source.clone(),
            ROOT_DOMAIN_ROWS,
        );
        request.root_domain = Some(RootDomainRequest::default());
        let attempt = match transform_recurrence(&request) {
            Err(error) => VeroneseParentProofAttempt {
                route: VeroneseParentProofRoute::RootDomain,
                variant: variant(name),
                success: false,
                detail: format!("the transform failed: {error}"),
            },
            Ok(report) => {
                let claim = report.root_domain.as_ref().and_then(|root| {
                    root.claims
                        .iter()
                        .find(|claim| accepted_source_claim(claim, record.first_index))
                });
                match claim {
                    Some(claim) => VeroneseParentProofAttempt {
                        route: VeroneseParentProofRoute::RootDomain,
                        variant: variant(name),
                        success: true,
                        detail: format!(
                            "source claim: every P_n, n >= {}, is real-rooted in {} ({:?}): \
                             {}{moved_note}",
                            claim.first_index,
                            claim
                                .interval
                                .as_ref()
                                .map_or_else(String::new, ToString::to_string),
                            claim.evidence,
                            claim.derivation
                        ),
                    },
                    None => VeroneseParentProofAttempt {
                        route: VeroneseParentProofRoute::RootDomain,
                        variant: variant(name),
                        success: false,
                        detail: "no unconditional all-index source claim with roots in \
                                 (-inf, 0]"
                            .into(),
                    },
                }
            }
        };
        let success = attempt.success;
        attempts.push(attempt);
        if success {
            return (attempts, true);
        }
    }
    (attempts, false)
}

/// A theorem route of the registry ([`TheoremRoute::parent`]): the theorem
/// with the requested region of its [`ParentRegion`] policy, proved for the
/// parent only when the certified region lies in `(-inf, 0]` (every parent
/// row then has coefficients of one sign).  The last value is whether the
/// theorem applies (status other than `not_applicable`).
fn theorem_attempt(
    record: &ParentRecurrenceRecord<'_>,
    route: &TheoremRoute,
) -> Result<(VeroneseParentProofAttempt, Option<LiftInterlacing>, bool), RecurrenceError> {
    let check = |options: &RealRootednessTheoremOptions| {
        check_real_rootedness_theorem(
            route.theorem,
            record.rule,
            record.first_index,
            record.seeds,
            options,
        )
    };
    let ray = RealRootednessTheoremOptions {
        root_region: Some(RealInterval {
            lower: None,
            lower_closed: false,
            upper: Some(Q::zero()),
            upper_closed: false,
        }),
        ..RealRootednessTheoremOptions::default()
    };
    let unrestricted = RealRootednessTheoremOptions::default();
    let report = match route.parent.region {
        ParentRegion::OpenNonpositiveRay => check(&ray)?,
        ParentRegion::Unrestricted => check(&unrestricted)?,
        ParentRegion::OpenNonpositiveRayThenUnrestricted => {
            let requested = check(&ray)?;
            if requested.is_proved()
                || requested.status == RealRootednessTheoremStatus::NotApplicable
            {
                requested
            } else {
                check(&unrestricted)?
            }
        }
    };
    let route_name = VeroneseParentProofRoute::of_theorem(route.theorem);
    let variant = route_name.as_str();
    // A proof whose certified region is not contained in (-inf, 0] (the
    // Ma--Wang interval, the Liu--Wang three-term or iterated-preserver
    // region, ...) gives real-rooted rows but not one-signed coefficients.
    let outside_ray = report.is_proved()
        && !report
            .conclusion
            .as_ref()
            .and_then(|c| c.root_region.as_ref())
            .and_then(|region| region.at(0).ok())
            .and_then(|region| region.upper)
            .is_some_and(|upper| !upper.is_positive());
    let detail = match (&report.status, &report.conclusion) {
        (RealRootednessTheoremStatus::Proved, Some(conclusion)) if outside_ray => format!(
            "{}: {}; the {} not contained in (-inf, 0], so the rows are real-rooted but this \
             route does not give coefficients of one sign",
            report.normal_form,
            conclusion.statement,
            if conclusion.root_region.is_some() {
                "certified region is"
            } else {
                "theorem certifies no root region, so the roots are"
            }
        ),
        (RealRootednessTheoremStatus::Proved, Some(conclusion)) if report.start_delay > 0 => {
            format!(
                "induction start moved {} row(s) past the seeds: {}: {}",
                report.start_delay, report.normal_form, conclusion.statement
            )
        }
        (RealRootednessTheoremStatus::Proved, Some(conclusion)) => {
            format!("{}: {}", report.normal_form, conclusion.statement)
        }
        (RealRootednessTheoremStatus::NotApplicable, _) => {
            format!("not applicable: {}", report.normal_form)
        }
        (status, _) => format!(
            "{status:?}: {}; open hypotheses: {}",
            report.normal_form,
            report.open_hypotheses().join("; ")
        ),
    };
    let success = report.is_proved() && !outside_ray;
    let interlacing = success
        .then(|| theorem_report_interlacing(&report, variant))
        .flatten();
    let applicable = report.status != RealRootednessTheoremStatus::NotApplicable;
    Ok((
        VeroneseParentProofAttempt {
            route: route_name,
            variant: variant.into(),
            success,
            detail,
        },
        interlacing,
        applicable,
    ))
}

/// The interlacing conclusion of a proved theorem report: `P_(m-1)`
/// interlaces `P_m` from the left (`P_m` carries the largest zero, the
/// orientation of [`RowInterlacing`]) for every `m >= interlacing_from`.
pub(super) fn theorem_report_interlacing(
    report: &RealRootednessTheoremReport,
    variant: &str,
) -> Option<LiftInterlacing> {
    let conclusion = report.conclusion.as_ref().filter(|_| report.is_proved())?;
    let (kind, from_index) = (conclusion.interlacing?, conclusion.interlacing_from?);
    let mut premises = vec![format!(
        "every hypothesis of the theorem is verified exactly for every index: {}",
        report.normal_form
    )];
    premises.push(conclusion.statement.clone());
    Some(LiftInterlacing {
        kind,
        orientation: InterlacingOrientation::PreviousInterlacesNext,
        from_index,
        justification: LiftInterlacingJustification {
            theorem: format!("{variant}: {}", report.reference),
            premises,
        },
    })
}

/// Whether two records define the same parent: same first index, the same
/// relation term by term, the same number of identifying seeds and equal
/// seeds (the supplied rows of `other` must agree with its rule).
fn same_record(
    record: &ParentRecurrenceRecord<'_>,
    other: &RecurrenceJson,
) -> Result<bool, RecurrenceError> {
    let (rule, first_index, rows) = other.to_recurrence_parts()?;
    if first_index != record.first_index || rows.len() < rule.max_offset() {
        return Ok(false);
    }
    let (left, left_shift) = ImplicitRelation::from_recurrence(record.rule)?;
    let (right, right_shift) = ImplicitRelation::from_recurrence(&rule)?;
    if left_shift != right_shift || left.op() != right.op() {
        return Ok(false);
    }
    let count = rule.generation_initial_count(first_index, rows.len());
    if count != record.seeds.len() {
        return Ok(false);
    }
    let generated = match rule.generate_rows_rational(&rows[..count], first_index, rows.len()) {
        Ok(generated) => generated,
        Err(_) => return Ok(false),
    };
    let equal = |a: &[Q], b: &[Q]| vec_poly::trim(a) == vec_poly::trim(b);
    Ok(rows.iter().zip(&generated).all(|(a, b)| equal(a, b))
        && rows[..count]
            .iter()
            .zip(record.seeds)
            .all(|(a, b)| equal(a, b)))
}

/// A bivariate polynomial `sum a[i][j] v^i w^j` (dense, exact).
type Bivariate = Vec<Vec<Q>>;

fn bivariate_add(left: &Bivariate, right: &Bivariate) -> Bivariate {
    let rows = left.len().max(right.len());
    (0..rows)
        .map(|i| {
            let a = left.get(i).map_or(&[][..], Vec::as_slice);
            let b = right.get(i).map_or(&[][..], Vec::as_slice);
            (0..a.len().max(b.len()))
                .map(|j| {
                    a.get(j).cloned().unwrap_or_default() + b.get(j).cloned().unwrap_or_default()
                })
                .collect()
        })
        .collect()
}

fn bivariate_mul(left: &Bivariate, right: &Bivariate) -> Bivariate {
    let mut out: Bivariate = Vec::new();
    for (i, row) in left.iter().enumerate() {
        for (j, a) in row.iter().enumerate() {
            if a.is_zero() {
                continue;
            }
            for (k, other) in right.iter().enumerate() {
                for (l, b) in other.iter().enumerate() {
                    if b.is_zero() {
                        continue;
                    }
                    if out.len() <= i + k {
                        out.resize(i + k + 1, Vec::new());
                    }
                    if out[i + k].len() <= j + l {
                        out[i + k].resize(j + l + 1, Q::zero());
                    }
                    out[i + k][j + l] += a * b;
                }
            }
        }
    }
    out
}

/// `p(a + b v + c w)` for a univariate `p`, as a bivariate polynomial in `(v, w)`.
fn bivariate_affine_powers(a: &Q, b: &Q, c: &Q, degree: usize) -> Vec<Bivariate> {
    let affine: Bivariate = vec![vec![a.clone(), c.clone()], vec![b.clone()]];
    let mut powers = vec![vec![vec![Q::one()]]];
    for _ in 0..degree {
        let next = bivariate_mul(powers.last().expect("nonempty"), &affine);
        powers.push(next);
    }
    powers
}

/// Stirling numbers of the second kind `S(m, i)`, `m, i <= max`.
fn stirling_second_kind_table(max: usize) -> Vec<Vec<Q>> {
    let mut table = vec![vec![Q::zero(); max + 1]; max + 1];
    table[0][0] = Q::one();
    for m in 1..=max {
        for i in 1..=m {
            table[m][i] = &table[m - 1][i - 1] + Q::from_integer(i.into()) * &table[m - 1][i];
        }
    }
    table
}

/// Whether every coefficient of `p` in the basis `v^(falling i) w^(falling j)`
/// is `>= 0`; then `p(v, w) >= 0` at all integers `v, w >= 0`
/// (`x^m = sum_i S(m, i) x^(falling i)`).
fn has_nonnegative_falling_factorial_coefficients(p: &Bivariate) -> bool {
    let max = p
        .len()
        .max(p.iter().map(Vec::len).max().unwrap_or(0))
        .max(1);
    let stirling = stirling_second_kind_table(max);
    let mut falling = vec![vec![Q::zero(); max]; max];
    for (m, row) in p.iter().enumerate() {
        for (l, a) in row.iter().enumerate() {
            if a.is_zero() {
                continue;
            }
            for i in 0..=m {
                for j in 0..=l {
                    let weight = &stirling[m][i] * &stirling[l][j];
                    if !weight.is_zero() {
                        falling[i][j] += a * weight;
                    }
                }
            }
        }
    }
    falling.iter().flatten().all(|value| !value.is_negative())
}

/// `x (x - 1) ... (x - d + 1)` in ascending powers.
fn falling_factorial_polynomial(d: usize) -> Vec<Q> {
    (0..d).fold(vec![Q::one()], |acc, k| {
        vec_poly::mul(&acc, &[-Q::from_integer(k.into()), Q::one()])
    })
}

/// Degree in `x` of `c(n, x)` (`None` for the zero polynomial).
fn x_degree(poly: &BivarPoly) -> Option<usize> {
    poly.coeffs
        .iter()
        .filter_map(|row| row.iter().rposition(|value| !value.is_zero()))
        .max()
}

/// The coefficient of `x^j` of `c(n, x)` as a polynomial in `n`.
fn x_coefficient(poly: &BivarPoly, j: usize) -> Vec<Q> {
    poly.coeffs
        .iter()
        .map(|row| row.get(j).cloned().unwrap_or_default())
        .collect()
}

fn integer_q(value: i64) -> Q {
    Q::from_integer(value.into())
}

/// Route 10: an exact proof that every row from `first_index` on has
/// coefficients of one sign (module documentation, "One-signed rows").
/// Returns what was shown, or why the sufficient test does not apply.
fn one_signed_rows_for_all_indices(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
) -> Result<String, String> {
    let start = first_index + seeds.len();
    let sign = seeds
        .iter()
        .flat_map(|seed| seed.iter())
        .find(|value| !value.is_zero())
        .map_or(Q::one(), |value| {
            if value.is_negative() {
                -Q::one()
            } else {
                Q::one()
            }
        });
    let inhomogeneous = rule
        .inhomogeneous
        .as_ref()
        .filter(|g| g.coeffs.iter().flatten().any(|value| !value.is_zero()));
    if sign.is_negative() && inhomogeneous.is_some() {
        return Err("the seeds are nonpositive but the rule is inhomogeneous".into());
    }
    if let Some(position) = seeds
        .iter()
        .position(|seed| seed.iter().any(|value| (value * &sign).is_negative()))
    {
        return Err(format!("seed {position} has coefficients of both signs"));
    }
    let positive_from_start = |polynomial: &[Q], condition: IndexSignCondition| {
        decide_index_polynomial_sign(
            polynomial,
            &IntegerRange::starting_at(BigInt::from(start)),
            condition,
        )
        .holds()
    };
    if let Some(leader) = &rule.denominator {
        if x_degree(leader).is_some_and(|degree| degree > 0) {
            return Err("the leader depends on x".into());
        }
        if !positive_from_start(&x_coefficient(leader, 0), IndexSignCondition::Positive) {
            return Err(format!("the leader is not positive for every n >= {start}"));
        }
    }
    if rule
        .terms
        .iter()
        .any(|term| term.sign == RecurrenceSign::AlternatingN)
    {
        return Err("a term carries the factor (-1)^n".into());
    }
    if let Some(g) = inhomogeneous {
        for j in 0..=x_degree(g).unwrap_or(0) {
            if !positive_from_start(&x_coefficient(g, j), IndexSignCondition::Nonnegative) {
                return Err(format!(
                    "the coefficient of x^{j} of the inhomogeneous term is not >= 0 for every \
                     n >= {start}"
                ));
            }
        }
    }

    // Degree bound D(m) = alpha + beta m for every row m >= first_index.
    let mut beta = 0i64;
    for term in &rule.terms {
        let Some(degree) = x_degree(&term.coeff) else {
            continue;
        };
        let excess = i64::try_from(degree).unwrap_or(i64::MAX)
            - i64::try_from(term.deriv_order).unwrap_or(i64::MAX);
        let offset = i64::try_from(term.offset).unwrap_or(i64::MAX).max(1);
        if excess > 0 {
            beta = beta.max((excess + offset - 1) / offset);
        }
    }
    let as_i64 = |value: usize| i64::try_from(value).unwrap_or(i64::MAX);
    let mut alpha: Option<i64> = None;
    for (position, seed) in seeds.iter().enumerate() {
        if let Some(degree) = vec_poly::degree(seed) {
            let needed = as_i64(degree) - beta * as_i64(first_index + position);
            alpha = Some(alpha.map_or(needed, |alpha| alpha.max(needed)));
        }
    }
    if let Some(degree) = inhomogeneous.and_then(x_degree) {
        let needed = as_i64(degree) - beta * as_i64(start);
        alpha = Some(alpha.map_or(needed, |alpha| alpha.max(needed)));
    }
    let Some(alpha) = alpha else {
        return Ok("every seed is zero and the rule is homogeneous, so every row is zero".into());
    };

    // Combined weights W_(r,e)(n, s).
    let mut weights: std::collections::BTreeMap<(usize, i64), Bivariate> =
        std::collections::BTreeMap::new();
    for term in &rule.terms {
        let falling = falling_factorial_polynomial(term.deriv_order);
        for j in 0..=x_degree(&term.coeff).unwrap_or(0) {
            let c = x_coefficient(&term.coeff, j);
            if c.iter().all(Zero::is_zero) {
                continue;
            }
            let shift = as_i64(j) - as_i64(term.deriv_order);
            let product: Bivariate = c
                .iter()
                .map(|cn| falling.iter().map(|value| cn * value).collect())
                .collect();
            let entry = weights.entry((term.offset, shift)).or_default();
            *entry = bivariate_add(entry, &product);
        }
    }
    let start_q = integer_q(as_i64(start));
    for ((offset, shift), weight) in &weights {
        if weight.iter().flatten().all(Zero::is_zero) {
            continue;
        }
        let n_degree = weight.len().saturating_sub(1);
        let s_degree = weight
            .iter()
            .map(Vec::len)
            .max()
            .unwrap_or(1)
            .saturating_sub(1);
        let n_powers = bivariate_affine_powers(&start_q, &Q::one(), &Q::zero(), n_degree);
        let substitute = |s_powers: &[Bivariate]| {
            let mut total: Bivariate = Vec::new();
            for (a, row) in weight.iter().enumerate() {
                for (b, value) in row.iter().enumerate() {
                    if value.is_zero() {
                        continue;
                    }
                    let monomial = bivariate_mul(&n_powers[a], &s_powers[b]);
                    let scaled: Bivariate = monomial
                        .iter()
                        .map(|row| row.iter().map(|x| x * value).collect())
                        .collect();
                    total = bivariate_add(&total, &scaled);
                }
            }
            total
        };
        // (a) s >= 0 free.
        let free = bivariate_affine_powers(&Q::zero(), &Q::zero(), &Q::one(), s_degree);
        if has_nonnegative_falling_factorial_coefficients(&substitute(&free)) {
            continue;
        }
        // (b) s = D(n - r) - u with n = start + v.
        let base = alpha + beta * (as_i64(start) - as_i64(*offset));
        let bounded =
            bivariate_affine_powers(&integer_q(base), &integer_q(beta), &-Q::one(), s_degree);
        if has_nonnegative_falling_factorial_coefficients(&substitute(&bounded)) {
            continue;
        }
        return Err(format!(
            "the weight of the coefficient of x^s of P_(n-{offset}) in the coefficient of \
             x^(s{shift:+}) is not shown to be >= 0 for n >= {start}"
        ));
    }
    Ok(format!(
        "every row has coefficients of one sign: the seeds do, the leader is positive and every \
         combined weight is >= 0 for n >= {start}, with deg P_m <= {alpha} + {beta} m"
    ))
}

/// Attestations from bundled catalog entries flagged `lean_proved` whose
/// first index and leading initial rows equal the parent's (a cheap filter;
/// [`same_record`] decides the exact match).
fn catalog_flag_attestations(
    record: &ParentRecurrenceRecord<'_>,
) -> Vec<VeroneseParentAttestation> {
    let wanted: Vec<Vec<String>> = record
        .seeds
        .iter()
        .map(|seed| {
            vec_poly::trim(seed)
                .iter()
                .map(format_rational_coeff)
                .collect()
        })
        .collect();
    let trimmed = |row: &[&str]| -> Vec<String> {
        let end = row
            .iter()
            .rposition(|value| *value != "0")
            .map_or(0, |i| i + 1);
        row[..end]
            .iter()
            .map(|value| (*value).to_string())
            .collect()
    };
    let zero_row = vec!["0".to_string()];
    let mut attestations = Vec::new();
    for (entry, flag) in entries_with_real_rootedness_flag() {
        if flag.status != OeisRealRootednessStatus::LeanProved
            || entry.recurrence_first_index != record.first_index
            || entry.initial_rows.len() < wanted.len()
        {
            continue;
        }
        let seeds_agree = entry.initial_rows.iter().zip(&wanted).all(|(row, seed)| {
            let row = trimmed(row);
            row == *seed || (row.is_empty() && *seed == zero_row)
        });
        if !seeds_agree {
            continue;
        }
        let Ok((rule, initial_rows)) = entry.recurrence_parts() else {
            continue;
        };
        let as_of = flag
            .as_of
            .map_or_else(String::new, |date| format!(", as of {date}"));
        attestations.push(VeroneseParentAttestation {
            label: format!("oeis_catalog:{}:lean_proved", entry.id),
            recurrence: RecurrenceJson::from_recurrence_rational(
                &rule,
                entry.recurrence_first_index,
                &initial_rows,
                None,
            ),
            first_index: entry.recurrence_first_index,
            statement: format!(
                "every row of the bundled OEIS catalog entry {} is zero or real-rooted (flag \
                 lean_proved{as_of}: a machine-checked Lean proof made outside polytool)",
                entry.id
            ),
            source: format!(
                "polytool OEIS catalog, real_rootedness flag of {}",
                entry.id
            ),
        });
    }
    attestations
}

/// First checked row that is nonzero and not real-rooted of one sign.
fn first_failing_row(rows: &[Vec<Q>], first_index: usize) -> Option<usize> {
    rows.iter()
        .position(|row| {
            let shape = classify_row(row);
            !shape.is_zero && !(shape.has_one_sign && shape.is_real_rooted)
        })
        .map(|offset| first_index + offset)
}

/// Run the parent provers (module documentation): try to prove that every
/// row `P_n`, `n >= record.first_index`, of the parent is zero or
/// real-rooted with coefficients of one sign, for every `n`.
///
/// `needed_from` is the first parent index the caller uses (at least
/// `first_index`); `attestations` are the caller-supplied statements of
/// route 9 (a verdict resting on one is conditional, never proved);
/// `try_root_domain` enables route 3.  A refuted or unproved statement is
/// an `Ok` proof record with that status.
///
/// # Errors
///
/// [`RecurrenceError::InvalidInput`] when the parent rule cannot be
/// replayed from its seeds, and errors of the theorem checks on malformed
/// records.
pub fn prove_parent_real_rootedness(
    record: &ParentRecurrenceRecord<'_>,
    needed_from: usize,
    attestations: &[VeroneseParentAttestation],
    try_root_domain: bool,
) -> Result<VeroneseParentProof, RecurrenceError> {
    prove_parent_with_routes(record, needed_from, attestations, try_root_domain, None)
}

/// Internal closure policy: the structural mode never enters either attestation route.
pub(in crate::recurrence) fn prove_parent_with_policy(
    record: &ParentRecurrenceRecord<'_>,
    needed_from: usize,
    attestations: &[VeroneseParentAttestation],
    try_root_domain: bool,
    allow_attestations: bool,
) -> Result<VeroneseParentProof, RecurrenceError> {
    if allow_attestations {
        return prove_parent_real_rootedness(record, needed_from, attestations, try_root_domain);
    }
    let mut routes: Vec<_> = RealRootednessTheorem::ALL
        .iter()
        .copied()
        .map(VeroneseParentProofRoute::of_theorem)
        .collect();
    routes.extend([
        VeroneseParentProofRoute::PeriodicProduct,
        VeroneseParentProofRoute::PaperTheorem,
    ]);
    if try_root_domain {
        routes.push(VeroneseParentProofRoute::RootDomain);
    }
    prove_parent_by_routes(record, needed_from, &routes)
}

/// Replay of the routes `routes` alone (module documentation, "Replaying
/// recorded proofs"): the provers with every other route skipped (no
/// attempt recorded for it), the root-domain route enabled when it is
/// among `routes`, and no attestations.  The result is proved only by one
/// of `routes`; it is a correct proof whenever it is proved, like every
/// result of the provers.
pub(in crate::recurrence) fn prove_parent_by_routes(
    record: &ParentRecurrenceRecord<'_>,
    needed_from: usize,
    routes: &[VeroneseParentProofRoute],
) -> Result<VeroneseParentProof, RecurrenceError> {
    let try_root_domain = routes.contains(&VeroneseParentProofRoute::RootDomain);
    prove_parent_with_routes(record, needed_from, &[], try_root_domain, Some(routes))
}

/// Run the theorem routes of `stage` in the registry order
/// ([`TheoremRoute::parent`]): each route runs while the stage allows it
/// (unproved, or for [`ParentStage::ReplacesEarly`] proved only by an early
/// route), its attempt is recorded (unless it does not apply and the route
/// records only applicable attempts), and its proof sets `proved` and the
/// interlacing.
fn run_theorem_routes(
    record: &ParentRecurrenceRecord<'_>,
    stage: ParentStage,
    allowed: &dyn Fn(VeroneseParentProofRoute) -> bool,
    proved: &mut Option<VeroneseParentProofRoute>,
    interlacing: &mut Option<LiftInterlacing>,
    attempts: &mut Vec<VeroneseParentProofAttempt>,
) -> Result<(), RecurrenceError> {
    for route in ROUTES.iter().filter(|route| route.parent.stage == stage) {
        let proved_early = proved
            .and_then(VeroneseParentProofRoute::theorem)
            .is_some_and(|theorem| {
                ROUTES.iter().any(|early| {
                    early.theorem == theorem && early.parent.stage == ParentStage::Early
                })
            });
        let runs = match stage {
            ParentStage::Early | ParentStage::WhileUnproved => proved.is_none(),
            ParentStage::ReplacesEarly => proved.is_none() || proved_early,
        };
        let route_name = VeroneseParentProofRoute::of_theorem(route.theorem);
        if !runs || !allowed(route_name) {
            continue;
        }
        let (attempt, theorem_interlacing, applicable) = theorem_attempt(record, route)?;
        if route.parent.record_only_if_applicable && !applicable {
            continue;
        }
        let success = attempt.success;
        attempts.push(attempt);
        if success {
            *proved = Some(route_name);
            *interlacing = theorem_interlacing;
        }
    }
    Ok(())
}

/// The parent provers, restricted to `only` when given (attestation and
/// catalog-flag routes are then skipped too).
fn prove_parent_with_routes(
    record: &ParentRecurrenceRecord<'_>,
    needed_from: usize,
    attestations: &[VeroneseParentAttestation],
    try_root_domain: bool,
    only: Option<&[VeroneseParentProofRoute]>,
) -> Result<VeroneseParentProof, RecurrenceError> {
    let allowed = |route: VeroneseParentProofRoute| only.is_none_or(|r| r.contains(&route));
    let rows = record
        .rule
        .generate_rows_rational(
            record.seeds,
            record.first_index,
            PARENT_ROWS_CHECKED.max(record.seeds.len()),
        )
        .map_err(|error| {
            RecurrenceError::InvalidInput(format!(
                "the parent recurrence cannot be replayed: {error}"
            ))
        })?;
    // Routes 1 to 8 and the one-signed test of route 10 run on the rule with
    // a left-hand factor free of x when the normalization applies; a rule
    // that does not normalize is kept, and the theorem routes name the reason.
    let (normalized_rule, normalization) = rule_with_x_free_left_hand_factor(record.rule)
        .unwrap_or_else(|_| (record.rule.clone(), None));
    let normalized = ParentRecurrenceRecord {
        rule: &normalized_rule,
        seeds: record.seeds,
        first_index: record.first_index,
    };
    let mut proof = VeroneseParentProof {
        status: VeroneseParentProofStatus::NotProved,
        route: None,
        statement: String::new(),
        conditional_on: Vec::new(),
        attempts: Vec::new(),
        failing_row_index: first_failing_row(&rows, record.first_index),
        rows_checked: rows.len(),
        interlacing: None,
        left_hand_factor_normalization: normalization,
    };
    let conclusion = format!(
        "every parent row P_n, n >= {}, is zero or real-rooted with coefficients of one sign",
        record.first_index
    );
    if let Some(index) = proof.failing_row_index {
        proof.status = VeroneseParentProofStatus::Refuted;
        proof.statement = format!(
            "parent row P_{index} is nonzero and not real-rooted with coefficients of one sign; \
             this parent gives no proof"
        );
        return Ok(proof);
    }

    // Route 1.
    let mut proved = None;
    if allowed(VeroneseParentProofRoute::PeriodicProduct) {
        let periodic = periodic_product_parent(normalized.rule, normalized.seeds);
        proof.attempts.push(VeroneseParentProofAttempt {
            route: VeroneseParentProofRoute::PeriodicProduct,
            variant: "P_(n+k) = c(x) P_n".into(),
            success: periodic.is_ok(),
            detail: match &periodic {
                Ok(()) => "c and the seeds are nonzero, real-rooted and of one sign, so \
                           P_(n0 + k m + i) = c^m P_(n0 + i) is too"
                    .into(),
                Err(reason) => reason.clone(),
            },
        });
        proved = periodic
            .is_ok()
            .then_some(VeroneseParentProofRoute::PeriodicProduct);
    }
    let mut interlacing = None;

    // The early theorem routes (product induction): real-rootedness only;
    // a later theorem route may replace the proof.
    run_theorem_routes(
        &normalized,
        ParentStage::Early,
        &allowed,
        &mut proved,
        &mut interlacing,
        &mut proof.attempts,
    )?;

    // Route 2.
    if proved.is_none() && allowed(VeroneseParentProofRoute::PaperTheorem) {
        let (attempts, window_interlacing) = theorem_route(&normalized, &rows);
        if attempts.iter().any(|a| a.success) {
            proved = Some(VeroneseParentProofRoute::PaperTheorem);
            interlacing = window_interlacing;
        }
        proof.attempts.extend(attempts);
    }

    // Route 3.
    if proved.is_none() && try_root_domain && allowed(VeroneseParentProofRoute::RootDomain) {
        let attempts = root_domain_route(&normalized, &rows);
        if attempts.iter().any(|a| a.success) {
            proved = Some(VeroneseParentProofRoute::RootDomain);
        }
        proof.attempts.extend(attempts);
    }

    // Routes 4 to 8b: product induction already settles real-rootedness,
    // but these routes may supply the stronger interlacing conclusion.
    // Then routes 8c to 8f, for a parent left open.
    for stage in [ParentStage::ReplacesEarly, ParentStage::WhileUnproved] {
        run_theorem_routes(
            &normalized,
            stage,
            &allowed,
            &mut proved,
            &mut interlacing,
            &mut proof.attempts,
        )?;
    }

    if let Some(route) = proved {
        proof.status = VeroneseParentProofStatus::Proved;
        proof.route = Some(route);
        proof.statement = conclusion;
        proof.interlacing = interlacing;
        return Ok(proof);
    }

    if only.is_some() {
        proof.statement = format!("open: the replayed routes do not prove that {conclusion}");
        return Ok(proof);
    }
    // Route 9.
    for attestation in attestations {
        let matches = same_record(record, &attestation.recurrence)?
            || (proof.left_hand_factor_normalization.is_some()
                && same_record(&normalized, &attestation.recurrence)?);
        let covers = attestation.first_index <= needed_from;
        let success = matches && covers;
        proof.attempts.push(VeroneseParentProofAttempt {
            route: VeroneseParentProofRoute::Attestation,
            variant: attestation.label.clone(),
            success,
            detail: if success {
                format!("{} (source: {})", attestation.statement, attestation.source)
            } else if !matches {
                "the attested record differs from the parent".into()
            } else {
                format!(
                    "the attestation starts at index {}, after the first used parent index \
                     {needed_from}",
                    attestation.first_index
                )
            },
        });
        if success {
            proof.status = VeroneseParentProofStatus::ConditionalOnAttestation;
            proof.route = Some(VeroneseParentProofRoute::Attestation);
            proof.statement = format!("{conclusion} (attested, not proved here)");
            proof.conditional_on = vec![attestation.label.clone()];
            return Ok(proof);
        }
    }
    // Route 10.
    let candidates = catalog_flag_attestations(record);
    let mut differing = Vec::new();
    for attestation in candidates {
        let matches = same_record(record, &attestation.recurrence)?
            || (proof.left_hand_factor_normalization.is_some()
                && same_record(&normalized, &attestation.recurrence)?);
        if !matches {
            differing.push(attestation.label);
            continue;
        }
        let covers = attestation.first_index <= needed_from;
        let one_signed = covers.then(|| {
            one_signed_rows_for_all_indices(normalized.rule, record.first_index, record.seeds)
        });
        let success = matches!(one_signed, Some(Ok(_)));
        proof.attempts.push(VeroneseParentProofAttempt {
            route: VeroneseParentProofRoute::Attestation,
            variant: attestation.label.clone(),
            success,
            detail: match &one_signed {
                Some(Ok(signs)) => format!("{}; {signs}", attestation.statement),
                Some(Err(reason)) => format!(
                    "the parent equals the flagged entry, but the flag gives real-rootedness \
                     only and one-signed coefficients are not shown: {reason}"
                ),
                None => format!(
                    "the flagged entry starts at index {}, after the first used parent index \
                     {needed_from}",
                    attestation.first_index
                ),
            },
        });
        if success {
            proof.status = VeroneseParentProofStatus::ConditionalOnAttestation;
            proof.route = Some(VeroneseParentProofRoute::Attestation);
            proof.statement = format!(
                "{conclusion} (real-rootedness from the catalog flag, not proved here; \
                 one-signed coefficients proved here)"
            );
            proof.conditional_on = vec![attestation.label.clone()];
            return Ok(proof);
        }
    }
    if !differing.is_empty() {
        proof.attempts.push(VeroneseParentProofAttempt {
            route: VeroneseParentProofRoute::Attestation,
            variant: "oeis_catalog".into(),
            success: false,
            detail: format!(
                "no lean_proved catalog entry equals the parent; {} with the same first index \
                 and seeds differ in the rule or in later initial rows ({})",
                differing.len(),
                differing.join(", ")
            ),
        });
    }
    proof.statement = format!(
        "open: no route proved that {conclusion}; the first {} parent rows pass the exact \
         check (evidence only)",
        proof.rows_checked
    );
    Ok(proof)
}

#[cfg(test)]
mod catalog_flag_tests {
    use super::*;

    fn q(value: i64) -> Q {
        Q::from_integer(value.into())
    }

    /// The catalog rule, identifying seeds and first index of a bundled entry.
    fn catalog_parent(id: &str) -> (Recurrence, Vec<Vec<Q>>, usize) {
        let entry = crate::oeis::by_id(id).unwrap();
        let (rule, rows) = entry.recurrence_parts().unwrap();
        let count = rule.generation_initial_count(entry.recurrence_first_index, rows.len());
        (rule, rows[..count].to_vec(), entry.recurrence_first_index)
    }

    fn prove(rule: &Recurrence, seeds: &[Vec<Q>], first_index: usize) -> VeroneseParentProof {
        let record = ParentRecurrenceRecord {
            rule,
            seeds,
            first_index,
        };
        prove_parent_real_rootedness(&record, first_index, &[], false).unwrap()
    }

    fn catalog_attempts(proof: &VeroneseParentProof) -> Vec<&VeroneseParentProofAttempt> {
        proof
            .attempts
            .iter()
            .filter(|attempt| attempt.variant.starts_with("oeis_catalog:"))
            .collect()
    }

    /// The Jabotinsky-type rows A039692 (`P_n = (x + n) P_(n-1) + n (n - 1 +
    /// 2x) P_(n-2)`, `P_0 = 1`, `P_1 = 3 + x`) are flagged `lean_proved`; no
    /// polytool route proves them, so the verdict is conditional on the flag,
    /// with one-signed rows proved.
    #[test]
    fn a_flagged_catalog_parent_is_conditional_on_the_flag() {
        let (rule, seeds, first_index) = catalog_parent("A039692");
        assert_eq!(seeds, vec![vec![q(1)], vec![q(3), q(1)]]);
        let proof = prove(&rule, &seeds, first_index);
        assert_eq!(
            proof.status,
            VeroneseParentProofStatus::ConditionalOnAttestation,
            "{proof:#?}"
        );
        assert_eq!(proof.route, Some(VeroneseParentProofRoute::Attestation));
        assert_eq!(proof.conditional_on, ["oeis_catalog:A039692:lean_proved"]);
        assert!(proof.statement.contains("not proved here"));
        let attempt = proof.attempts.last().unwrap();
        assert!(attempt.success);
        assert!(attempt.detail.contains("lean_proved, as of 2026-10-03"));
        assert!(attempt.detail.contains("coefficients of one sign"));
        // A caller attestation of the same record is applied first.
        let attestation = VeroneseParentAttestation {
            label: "caller".into(),
            recurrence: RecurrenceJson::from_recurrence_rational(&rule, first_index, &seeds, None),
            first_index,
            statement: "attested".into(),
            source: "test".into(),
        };
        let record = ParentRecurrenceRecord {
            rule: &rule,
            seeds: &seeds,
            first_index,
        };
        let proof =
            prove_parent_real_rootedness(&record, first_index, &[attestation], false).unwrap();
        assert_eq!(proof.conditional_on, ["caller"]);
    }

    /// A157153 is not flagged (its Lean proof is open), so no catalog
    /// attestation applies.
    #[test]
    fn an_unflagged_catalog_parent_gets_no_catalog_attestation() {
        assert!(crate::oeis::by_id("A157153")
            .unwrap()
            .real_rootedness()
            .is_none());
        let (rule, seeds, first_index) = catalog_parent("A157153");
        let proof = prove(&rule, &seeds, first_index);
        assert_ne!(
            proof.status,
            VeroneseParentProofStatus::ConditionalOnAttestation,
            "{proof:#?}"
        );
        // Flagged entries with the same seeds differ in the rule.
        assert!(catalog_attempts(&proof).is_empty(), "{proof:#?}");
        assert!(proof
            .attempts
            .iter()
            .filter(|attempt| attempt.route == VeroneseParentProofRoute::Attestation)
            .all(|attempt| !attempt.success && !attempt.detail.contains("A157153")));
    }

    /// Matching is exact: one changed rule coefficient or seed and the flag no
    /// longer applies.
    #[test]
    fn a_parent_differing_in_one_coefficient_or_seed_does_not_match() {
        let (rule, seeds, first_index) = catalog_parent("A039692");
        let mut changed_rule = rule.clone();
        changed_rule.terms[1].coeff.coeffs[1][1] = q(3);
        let proof = prove(&changed_rule, &seeds, first_index);
        assert_ne!(
            proof.status,
            VeroneseParentProofStatus::ConditionalOnAttestation
        );
        let attempts = catalog_attempts(&proof);
        assert!(attempts.is_empty(), "{proof:#?}");
        let summary = proof.attempts.last().unwrap();
        assert_eq!(summary.variant, "oeis_catalog");
        assert!(!summary.success);
        assert!(summary.detail.contains("oeis_catalog:A039692:lean_proved"));

        let mut changed_seeds = seeds.clone();
        changed_seeds[1] = vec![q(2)];
        let proof = prove(&rule, &changed_seeds, first_index);
        assert_ne!(
            proof.status,
            VeroneseParentProofStatus::ConditionalOnAttestation
        );
        assert!(catalog_attempts(&proof).is_empty());

        // The same rule and seeds at another first index do not match either.
        let proof = prove(&rule, &seeds, first_index + 1);
        assert!(catalog_attempts(&proof).is_empty());
    }

    /// The flag gives real-rootedness only: a flagged parent that no route
    /// proves and whose rows are not shown one-signed by the weight test
    /// (A008303) stays open.
    #[test]
    fn a_flagged_parent_without_a_one_sign_proof_stays_open() {
        let (rule, seeds, first_index) = catalog_parent("A008303");
        assert!(crate::oeis::by_id("A008303")
            .unwrap()
            .real_rootedness()
            .is_some());
        let proof = prove(&rule, &seeds, first_index);
        assert_eq!(proof.status, VeroneseParentProofStatus::NotProved);
        let attempts = catalog_attempts(&proof);
        assert_eq!(attempts.len(), 1);
        assert!(attempts[0]
            .detail
            .contains("one-signed coefficients are not shown: the weight"));
    }

    fn term(offset: usize, derivative: usize, coefficients: Vec<Vec<Q>>) -> RecurrenceTerm {
        RecurrenceTerm::new(
            offset,
            derivative,
            BivarPoly::from_coefficients_normalized(coefficients),
        )
    }

    #[test]
    fn one_signed_rows_need_nonnegative_combined_weights() {
        // Eulerian: P_n = (1 - x + n x) P_(n-1) + (x - x^2) P'_(n-1), P_1 = 1.
        // The weight n - 1 - s of x^s into x^(s+1) needs the degree bound
        // deg P_m <= m - 1.
        let eulerian = Recurrence {
            terms: vec![
                term(1, 0, vec![vec![q(1), q(-1)], vec![q(0), q(1)]]),
                term(1, 1, vec![vec![q(0), q(1), q(-1)]]),
            ],
            denominator: None,
            inhomogeneous: None,
        };
        let shown = one_signed_rows_for_all_indices(&eulerian, 1, &[vec![q(1)]]).unwrap();
        assert!(shown.ends_with("deg P_m <= -1 + 1 m"), "{shown}");
        // Nonpositive seeds of a homogeneous rule: the rows are nonpositive.
        assert!(one_signed_rows_for_all_indices(&eulerian, 1, &[vec![q(-1)]]).is_ok());
        assert!(
            one_signed_rows_for_all_indices(&eulerian, 1, &[vec![q(1), q(-1)]])
                .unwrap_err()
                .contains("both signs")
        );
        // A coefficient x - 2 x^2 makes the weight n - 2 - 2 s, negative at
        // s = n - 2 (the top coefficient of P_(n-1)).
        let mut too_negative = eulerian.clone();
        too_negative.terms[1] = term(1, 1, vec![vec![q(0), q(1), q(-2)]]);
        assert!(
            one_signed_rows_for_all_indices(&too_negative, 1, &[vec![q(1)]])
                .unwrap_err()
                .contains("not shown to be >= 0")
        );
        // A leader 3 - n vanishes at n = 3.
        let mut leader = eulerian.clone();
        leader.denominator = Some(BivarPoly::from_coefficients_normalized(vec![
            vec![q(3)],
            vec![q(-1)],
        ]));
        assert!(one_signed_rows_for_all_indices(&leader, 1, &[vec![q(1)]])
            .unwrap_err()
            .contains("leader is not positive"));
        // (-1)^n factors and negative inhomogeneous parts are refused.
        let mut alternating = eulerian.clone();
        alternating.terms[0].sign = RecurrenceSign::AlternatingN;
        assert!(one_signed_rows_for_all_indices(&alternating, 1, &[vec![q(1)]]).is_err());
        let mut inhomogeneous = eulerian.clone();
        inhomogeneous.inhomogeneous = Some(BivarPoly::from_coefficients_normalized(vec![
            vec![q(5)],
            vec![q(-1)],
        ]));
        assert!(
            one_signed_rows_for_all_indices(&inhomogeneous, 1, &[vec![q(1)]])
                .unwrap_err()
                .contains("inhomogeneous")
        );
        // Zero seeds of a homogeneous rule give zero rows.
        assert!(one_signed_rows_for_all_indices(&eulerian, 1, &[vec![]]).is_ok());
    }
}

#[cfg(test)]
mod interlacing_tests {
    use super::*;

    fn q(value: i64) -> Q {
        Q::from_integer(value.into())
    }

    fn prove(rule: &Recurrence, seeds: &[Vec<Q>], first_index: usize) -> VeroneseParentProof {
        let record = ParentRecurrenceRecord {
            rule,
            seeds,
            first_index,
        };
        prove_parent_real_rootedness(&record, first_index, &[], false).unwrap()
    }

    /// `P_n = P_(n-1) + x P_(n-2)` with the canonical rows `P_0 = 1`,
    /// `P_1 = 1 + x` of the fixed-delay window theorem: route 2 proves the
    /// parent and records `P_(m-1) << P_m` from `m = 1`.
    #[test]
    fn the_window_theorem_route_records_directed_interlacing() {
        let rule = Recurrence {
            terms: vec![
                RecurrenceTerm::new(
                    1,
                    0,
                    BivarPoly::from_coefficients_normalized(vec![vec![q(1)]]),
                ),
                RecurrenceTerm::new(
                    2,
                    0,
                    BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(1)]]),
                ),
            ],
            denominator: None,
            inhomogeneous: None,
        };
        let proof = prove(&rule, &[vec![q(1)], vec![q(1), q(1)]], 0);
        assert_eq!(proof.route, Some(VeroneseParentProofRoute::PaperTheorem));
        let interlacing = proof.interlacing.as_ref().expect("window interlacing");
        assert_eq!(interlacing.kind, RowInterlacing::Weak);
        assert_eq!(
            interlacing.orientation,
            InterlacingOrientation::PreviousInterlacesNext
        );
        assert_eq!(interlacing.from_index, 1);
        assert!(interlacing
            .justification
            .theorem
            .starts_with("fixed_delay_predecessor_window"));
        // Exact evidence on the first rows.
        let rows = rule
            .generate_rows_rational(&[vec![q(1)], vec![q(1), q(1)]], 0, 10)
            .unwrap();
        for pair in rows.windows(2) {
            let integer = |row: &[Q]| row.iter().map(|c| c.to_integer()).collect::<Vec<_>>();
            assert_eq!(
                crate::check_weak_interlacing_bigint_coeffs(&integer(&pair[0]), &integer(&pair[1])),
                Some(true)
            );
        }
    }

    /// The Eulerian rule `P_n = (1 + (n - 1) x) P_(n-1) + x (1 - x) P'_(n-1)`
    /// is proved by Liu--Wang (route 5) without the root-domain route; the
    /// theorem's weak interlacing is recorded.  A periodic product records
    /// none.
    #[test]
    fn theorem_routes_record_their_conclusion_and_products_record_none() {
        let entry = crate::oeis::by_id("A008292").unwrap();
        let (rule, rows) = entry.recurrence_parts().unwrap();
        let count = rule.generation_initial_count(entry.recurrence_first_index, rows.len());
        let proof = prove(&rule, &rows[..count], entry.recurrence_first_index);
        assert_eq!(proof.route, Some(VeroneseParentProofRoute::LiuWang));
        let interlacing = proof.interlacing.as_ref().expect("Liu--Wang interlacing");
        assert_eq!(interlacing.kind, RowInterlacing::Weak);
        assert_eq!(
            interlacing.orientation,
            InterlacingOrientation::PreviousInterlacesNext
        );
        assert!(interlacing.justification.theorem.starts_with("liu_wang: "));
        let json = serde_json::to_value(&proof).unwrap();
        assert_eq!(
            json["interlacing"]["orientation"],
            "previous_interlaces_next"
        );
        let back: VeroneseParentProof = serde_json::from_value(json).unwrap();
        assert_eq!(back, proof);

        // (1 + x)^n: a periodic product, no interlacing recorded.
        let product = Recurrence {
            terms: vec![RecurrenceTerm::new(
                1,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![q(1), q(1)]]),
            )],
            denominator: None,
            inhomogeneous: None,
        };
        let proof = prove(&product, &[vec![q(1)]], 0);
        assert_eq!(proof.route, Some(VeroneseParentProofRoute::PeriodicProduct));
        assert!(proof.interlacing.is_none());
        assert!(serde_json::to_value(&proof)
            .unwrap()
            .get("interlacing")
            .is_none());
    }
}

#[cfg(test)]
mod left_hand_factor_tests {
    use super::*;

    fn q(value: i64) -> Q {
        Q::from_integer(value.into())
    }

    /// Both sides of `rule` multiplied by `x`.
    fn times_x(rule: &Recurrence) -> Recurrence {
        let shift = |p: &BivarPoly| {
            BivarPoly::from_coefficients_normalized(
                p.coeffs
                    .iter()
                    .map(|row| std::iter::once(q(0)).chain(row.iter().cloned()).collect())
                    .collect(),
            )
        };
        let one = BivarPoly::from_coefficients_normalized(vec![vec![q(1)]]);
        Recurrence {
            terms: rule
                .terms
                .iter()
                .map(|term| RecurrenceTerm {
                    coeff: shift(&term.coeff),
                    ..term.clone()
                })
                .collect(),
            denominator: Some(shift(rule.denominator.as_ref().unwrap_or(&one))),
            inhomogeneous: rule.inhomogeneous.as_ref().map(shift),
        }
    }

    fn prove(rule: &Recurrence, seeds: &[Vec<Q>], first_index: usize) -> VeroneseParentProof {
        let record = ParentRecurrenceRecord {
            rule,
            seeds,
            first_index,
        };
        prove_parent_real_rootedness(&record, first_index, &[], false).unwrap()
    }

    /// The type-B Eulerian rule times `x` is proved by the same route, with
    /// the normalization recorded.
    #[test]
    fn an_x_power_left_factor_is_normalized_for_the_theorem_routes() {
        let rule = Recurrence {
            terms: vec![
                RecurrenceTerm::new(
                    1,
                    0,
                    BivarPoly::from_coefficients_normalized(vec![
                        vec![q(1), q(-1)],
                        vec![q(0), q(2)],
                    ]),
                ),
                RecurrenceTerm::new(
                    1,
                    1,
                    BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(2), q(-2)]]),
                ),
            ],
            denominator: None,
            inhomogeneous: None,
        };
        let seeds = vec![vec![q(1)]];
        let plain = prove(&rule, &seeds, 0);
        assert_eq!(plain.status, VeroneseParentProofStatus::Proved);
        assert!(plain.left_hand_factor_normalization.is_none());
        let lifted = prove(&times_x(&rule), &seeds, 0);
        assert_eq!(
            lifted.status,
            VeroneseParentProofStatus::Proved,
            "{lifted:#?}"
        );
        assert_eq!(lifted.route, plain.route);
        assert_eq!(lifted.interlacing.is_some(), plain.interlacing.is_some());
        let normalization = lifted.left_hand_factor_normalization.as_ref().unwrap();
        assert_eq!(normalization.x_power, 1);
        normalization.verify().unwrap();
    }

    /// The rows of A039692 times `x` still match their `lean_proved`
    /// catalog flag through the normalized rule.
    #[test]
    fn a_catalog_flag_matches_the_normalized_rule() {
        let entry = crate::oeis::by_id("A039692").unwrap();
        let (rule, rows) = entry.recurrence_parts().unwrap();
        let first_index = entry.recurrence_first_index;
        let count = rule.generation_initial_count(first_index, rows.len());
        let proof = prove(&times_x(&rule), &rows[..count], first_index);
        assert_eq!(
            proof.status,
            VeroneseParentProofStatus::ConditionalOnAttestation,
            "{proof:#?}"
        );
        assert_eq!(proof.conditional_on, ["oeis_catalog:A039692:lean_proved"]);
        assert!(proof.left_hand_factor_normalization.is_some());
    }

    /// A rule whose division by `x` happens only on the rows is kept as it
    /// is; the theorem routes name the reason.
    #[test]
    fn a_row_only_division_is_reported_by_the_theorem_routes() {
        // x P_n = (x^2 + n - 1) P_(n-1) - x P'_(n-1), rows x^n.
        let rule = Recurrence {
            terms: vec![
                RecurrenceTerm::new(
                    1,
                    0,
                    BivarPoly::from_coefficients_normalized(vec![
                        vec![q(-1), q(0), q(1)],
                        vec![q(1)],
                    ]),
                ),
                RecurrenceTerm::new(
                    1,
                    1,
                    BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(-1)]]),
                ),
            ],
            denominator: Some(BivarPoly::from_coefficients_normalized(vec![vec![
                q(0),
                q(1),
            ]])),
            inhomogeneous: None,
        };
        let proof = prove(&rule, &[vec![q(1)]], 0);
        assert!(proof.left_hand_factor_normalization.is_none());
        let liu_wang = proof
            .attempts
            .iter()
            .find(|attempt| attempt.route == VeroneseParentProofRoute::LiuWang)
            .unwrap();
        assert!(
            liu_wang.detail.contains("happens only on the rows"),
            "{}",
            liu_wang.detail
        );
    }
}
