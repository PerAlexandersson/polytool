//! Lifts: recognize a recurrence-defined sequence as a part of a parent
//! sequence whose real-rootedness (with a root-location or interlacing
//! hypothesis) implies that of the input.
//!
//! - Veronese-section lifts, [`find_veronese_section_lifts`]: the input is a
//!   coefficient section `V_(r,j)(P_n)` of a parent with nonpositive roots;
//!   [`certify_veronese_section_lift`] certifies one candidate lift for inputs
//!   with index-dependent and differential terms (forward certificate) and
//!   runs the parent provers.
//! - Symmetric-companion lifts, [`find_symmetric_companion_lifts`] and
//!   [`verify_symmetric_companion_certificate`]: the input is the symmetric
//!   part `A_n` of `P_n = A_n + x B_n` with a first-order derivative parent
//!   rule; the interlacing-preservation theorem applies under explicit sign
//!   conditions (the `symmetric_companion` module docs).
//! - Hadamard-factor lifts, [`find_hadamard_factor_lifts`]: the input
//!   coefficients are `lambda(n) u(n,k) p(n,k)` for a known hypergeometric
//!   factor `u` (binomials, factorials, rising factorials) and a parent
//!   `P_n = sum_k p(n,k) x^k`; the Pólya--Schur multiplier-sequence theorem
//!   transfers real-rootedness (the `hadamard_factor` module docs).
//! - Row-factor lifts, [`find_row_factor_lifts`]: every row is
//!   `F(x)^(e(n)) S_n` for a fixed real-rooted `F` and an affine exponent
//!   `e(n)`; the cofactor rule is derived and certified exactly (the
//!   `row_factor` module docs).
//! - Gamma lifts, [`find_gamma_lifts`] and [`certify_gamma_lift`]: every row
//!   is `lambda(n) Gamma_(n+s)(c z)` for the gamma polynomials `Gamma_m` of a
//!   palindromic parent (for example the peak polynomials from the Eulerian
//!   polynomials); the gamma criterion transfers real-rootedness with
//!   nonpositive roots (the `gamma_lift` module docs).
//! - Preserver-image lifts, [`find_preserver_image_lifts`] and
//!   [`certify_preserver_image_lift`]: every row is
//!   `lambda(n) T_n(P_(n+s))` for a proved parent `P` (a
//!   [`PreserverParentLibrary`] of named families and OEIS catalog entries)
//!   and a real-rootedness preserver `T_n` of a fixed menu (multiplier
//!   sequences `f(n, theta)` and hypergeometric multipliers,
//!   Laguerre--Pólya operators `phi(n, D)`, `D^j o x^m`, `x^m D^j`, affine
//!   substitutions, reversals and Möbius maps), with every side condition
//!   decided exactly for all `n` (the `preserver_image` module docs).
//! - Basis-map lifts, [`find_basis_map_lifts`], [`certify_basis_map_lift`]
//!   and [`apply_basis_map_to_family`]: every row is `Q_n = T_P(R_n)` for a
//!   basis map `T_P(u^k) = P_k(t)` ([`BasisMap`]) of a family `P` whose map
//!   is proved to send polynomials with roots in a source interval to
//!   polynomials with roots in a target interval
//!   ([`check_basis_map_preserves_pf`]: geometric, constant- and
//!   affine-current order-two, ray and interval differential, Hermite-type
//!   and index-dependent scalar families), and `R_n` is proved to have its
//!   roots in the source; the identity is certified by transporting `R`'s
//!   rule through the intertwining relation of `T_P` with multiplication
//!   by `u` (the `basis_map` module docs).
//!
//! Every lift carries a uniform [`LiftSummary`] with the shared
//! [`lift_verdict`].
//!
//! For the Veronese lifts, [`find_veronese_section_lifts`]'s documentation
//! states the transfer theorem, the supported input class, the algorithm and
//! the argument that the finitely many identifying rows prove the identity
//! for all indices.  The group is split into
//!
//! - `veronese_algebra`: coefficient sections, bivariate polynomials in the
//!   shift variable `L`, the Veronese norm (a determinant) and monic division;
//! - `veronese_search`: input validation, factorization of `chi(L, x^r)`
//!   ([`crate::factorization`]), divisor enumeration and the per-candidate
//!   work behind a crate-internal entry point that takes `G` directly;
//! - `veronese_report`: the report types and their
//!   `polytool.recurrence-veronese-lifts.v1` JSON form;
//! - `section_operators`: the transport `phi` of operators to section
//!   matrices (`x -> X`, `D -> DD`);
//! - `section_certify`: gauges, pseudo-division by a transported parent
//!   system and [`VeroneseForwardDivisionCertificate`];
//! - `parent_proof`: the parent provers (periodic product, paper theorem,
//!   root domain, three-term recurrence, Liu--Wang, Liu--Wang derivative
//!   lag, Ma--Wang interval, iterated preserver, attestation), shared with
//!   the search;
//! - `parent_theorems`: [`check_three_term_recurrence_theorem`],
//!   [`check_liu_wang_theorem`], [`check_liu_wang_derivative_lag_theorem`],
//!   [`check_liu_wang_three_term_theorem`], [`check_ma_wang_interval_theorem`] and
//!   [`check_iterated_preserver_theorem`], general real-rootedness theorems for
//!   recurrences with exactly checked hypotheses and a replayable
//!   [`RealRootednessTheoremReport`];
//! - `section_report`: [`certify_veronese_section_lift`] and its
//!   `polytool.recurrence-veronese-certificate.v1` report;
//! - `symmetric_companion`: operators, search, certificate and verifier, and
//!   the `polytool.recurrence-symmetric-companion.v1` report;
//! - `preserver_image`: operators, transfer theorems, certificates, the
//!   parent library and the detection of the preserver-image lifts, with the
//!   `polytool.recurrence-preserver-image.v1` and
//!   `polytool.recurrence-preserver-image-lifts.v1` reports;
//! - `lift_summary`: [`LiftSummary`], the uniform view of every lift
//!   (operation, parent, identity status, root transfer, hypotheses with
//!   exact sign decisions, parent proof, verdict), embedded as the field
//!   `lift_summary` of the three reports.
//!
//! Related: [`RecurrenceJson`] (a lift's parent is exported in that format),
//! [`crate::is_real_rooted_bigint_coeffs`] (the exact parent-row checks).

use super::*;

mod basis_map;
pub use basis_map::*;
mod gamma_lift;
pub use gamma_lift::{
    certify_gamma_lift, find_gamma_lifts, GammaDivisionQuotient, GammaLiftCandidate,
    GammaLiftCertificate, GammaLiftCertification, GammaLiftIdentityStatus, GammaLiftOptions,
    GammaLiftParent, GammaLiftRejection, GammaLiftReport, GammaLiftSearchOptions,
    GammaParentSymmetryStatus, GAMMA_LIFTS_SCHEMA, GAMMA_LIFT_CERTIFICATION_SCHEMA,
    MAX_GAMMA_LIFT_CHECKED_ROWS, MAX_GAMMA_LIFT_INDEX_SHIFT, MAX_GAMMA_LIFT_OPERATOR_TERMS,
};
mod hadamard_factor;
pub use hadamard_factor::prepare::{prepare_lift_rule, PreparedLiftRule};
pub use hadamard_factor::{
    find_hadamard_factor_lifts, HadamardAnnihilatorGenerator, HadamardAnnihilatorSystem,
    HadamardDivisionQuotient, HadamardFactor, HadamardFactorCertificate,
    HadamardFactorIdentityStatus, HadamardFactorLift, HadamardFactorLiftOptions,
    HadamardFactorLiftReport, HadamardFactorRejection, HadamardFactorRejectionKind,
    HadamardFactorTerm, HadamardHiddenCoefficients, HADAMARD_FACTOR_LIFTS_SCHEMA,
    MAX_HADAMARD_BINOMIAL_OFFSET, MAX_HADAMARD_FACTOR_ROWS, MAX_HADAMARD_OPERATOR_TERMS,
};
mod preserver_image;
pub use preserver_image::*;
mod lift_summary;
pub use lift_summary::{
    lift_verdict, InterlacingOrientation, LiftHypothesis, LiftHypothesisStatus, LiftIdentityStatus,
    LiftInterlacing, LiftInterlacingClaim, LiftInterlacingJustification, LiftOperation, LiftParent,
    LiftParentProof, LiftRootTransfer, LiftSummary, LiftVerdict, LIFT_LEAN_REPLAY_NOT_ATTEMPTED,
};
mod parent_proof;
pub use parent_proof::{
    check_recorded_parent_proof, prove_parent_real_rootedness, ParentRecurrenceRecord,
};
pub use parent_proof::{
    VeroneseParentAttestation, VeroneseParentProof, VeroneseParentProofAttempt,
    VeroneseParentProofRoute, VeroneseParentProofStatus,
};
mod parent_theorems;
pub use parent_theorems::{
    check_hypergeometric_row_theorem, check_iterated_preserver_theorem,
    check_liu_wang_derivative_lag_theorem, check_liu_wang_theorem,
    check_liu_wang_three_term_theorem, check_ma_wang_interval_theorem,
    check_product_recurrence_theorem, check_three_term_recurrence_theorem,
    check_weyl_step_preserver, discover_recurrence_proofs, is_trivially_real_rooted,
    IteratedPreserverCertificate, IteratedPreserverCoordinates, LeftHandFactorNormalization,
    RealRootednessTheorem, RealRootednessTheoremConclusion, RealRootednessTheoremOptions,
    RealRootednessTheoremReport, RealRootednessTheoremStatus, RecurrenceProofDiscoveryReport,
    RowInterlacing, WeylStepPreserverVerdict, REAL_ROOTEDNESS_THEOREM_SCHEMA,
    RECURRENCE_PROOF_DISCOVERY_SCHEMA,
};
mod row_factor;
pub use row_factor::{
    find_row_factor_lifts, RowFactorCertificate, RowFactorIdentityStatus, RowFactorLift,
    RowFactorLiftOptions, RowFactorLiftReport, RowFactorRejection, MAX_ROW_FACTOR_ROWS,
    ROW_FACTOR_LIFTS_SCHEMA,
};
mod section_certify;
pub use section_certify::{
    VeroneseDivisionQuotient, VeroneseForwardDivisionCertificate, VeroneseSectionGauge,
};
mod section_operators;
mod section_report;
pub use section_report::*;
mod veronese_algebra;
mod veronese_report;
pub use veronese_report::*;
mod veronese_search;
pub use veronese_search::{
    find_veronese_section_lifts, VeroneseSectionLiftOptions, MAX_VERONESE_SECTION_MODULUS,
};

#[cfg(test)]
mod section_tests;
#[cfg(test)]
mod veronese_tests;

mod symmetric_companion;
pub use symmetric_companion::*;
