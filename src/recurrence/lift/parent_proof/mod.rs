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

macro_rules! define_parent_routes {
    ($(#[$first_doc:meta])* $first:ident, $first_module:ident, $first_name:literal,
        $first_summary:literal, { $($first_export:ident),* $(,)? };
     $( $(#[$doc:meta])* $variant:ident, $module:ident, $name:literal,
        $summary:literal, { $($export:ident),* $(,)? }; )*) => {
        /// The route that settled the parent statement.
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
        #[serde(rename_all = "snake_case")]
        pub enum VeroneseParentProofRoute {
            /// An index-independent periodic product with one-signed real-rooted seeds.
            PeriodicProduct,
            $(#[$first_doc])* #[serde(rename = $first_name)] $first,
            /// A matched paper theorem.
            PaperTheorem,
            /// A certified root-domain claim.
            RootDomain,
            $($(#[$doc])* #[serde(rename = $name)] $variant,)*
            /// A caller-supplied attestation.
            Attestation,
        }

        impl VeroneseParentProofRoute {
            /// The stable snake-case name used by Serde.
            pub const fn as_str(self) -> &'static str {
                match self {
                    Self::PeriodicProduct => "periodic_product",
                    Self::PaperTheorem => "paper_theorem",
                    Self::RootDomain => "root_domain",
                    Self::Attestation => "attestation",
                    Self::$first => $first_name,
                    $(Self::$variant => $name,)*
                }
            }

            /// The theorem behind this route, if it is a theorem route.
            pub const fn theorem(self) -> Option<RealRootednessTheorem> {
                match self {
                    Self::$first => Some(RealRootednessTheorem::$first),
                    $(Self::$variant => Some(RealRootednessTheorem::$variant),)*
                    Self::PeriodicProduct | Self::PaperTheorem | Self::RootDomain
                        | Self::Attestation => None,
                }
            }

            /// The parent route of a theorem (the inverse of [`Self::theorem`]).
            pub(in crate::recurrence) const fn of_theorem(theorem: RealRootednessTheorem) -> Self {
                match theorem {
                    RealRootednessTheorem::$first => Self::$first,
                    $(RealRootednessTheorem::$variant => Self::$variant,)*
                }
            }
        }
    };
}
parent_theorems::theorem_routes!(define_parent_routes);

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
    let replay = match route.theorem() {
        Some(theorem) => {
            let theorem_route = ROUTES
                .iter()
                .find(|r| r.theorem == theorem)
                .expect("every theorem has a route");
            match theorem_route.policy.replay {
                ParentReplayPolicy::RecordedRoute => {
                    prove_parent_by_routes(record, needed_from, &[route])?
                }
            }
        }
        None => prove_parent_by_routes(record, needed_from, &[route])?,
    };
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

mod routes;
use routes::*;
mod coefficient_signs;
use coefficient_signs::*;
mod cascade;
pub use cascade::prove_parent_real_rootedness;
pub(in crate::recurrence) use cascade::{prove_parent_by_routes, prove_parent_with_policy};
#[cfg(test)]
mod tests;
