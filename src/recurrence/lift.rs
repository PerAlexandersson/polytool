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
//!   root domain, three-term recurrence, Liu--Wang, attestation), shared
//!   with the search;
//! - `parent_theorems`: [`check_three_term_recurrence_theorem`] and
//!   [`check_liu_wang_theorem`], general real-rootedness theorems for
//!   recurrences with exactly checked hypotheses and a replayable
//!   [`RealRootednessTheoremReport`];
//! - `section_report`: [`certify_veronese_section_lift`] and its
//!   `polytool.recurrence-veronese-certificate.v1` report;
//! - `symmetric_companion`: operators, search, certificate and verifier, and
//!   the `polytool.recurrence-symmetric-companion.v1` report;
//! - `lift_summary`: [`LiftSummary`], the uniform view of every lift
//!   (operation, parent, identity status, root transfer, hypotheses with
//!   exact sign decisions, parent proof, verdict), embedded as the field
//!   `lift_summary` of the three reports.
//!
//! Related: [`RecurrenceJson`] (a lift's parent is exported in that format),
//! [`crate::is_real_rooted_bigint_coeffs`] (the exact parent-row checks).

use super::*;

mod gamma_lift;
pub use gamma_lift::{
    certify_gamma_lift, find_gamma_lifts, GammaDivisionQuotient, GammaLiftCandidate,
    GammaLiftCertificate, GammaLiftCertification, GammaLiftIdentityStatus, GammaLiftOptions,
    GammaLiftParent, GammaLiftRejection, GammaLiftReport, GammaLiftSearchOptions,
    GammaParentSymmetryStatus, GAMMA_LIFTS_SCHEMA, GAMMA_LIFT_CERTIFICATION_SCHEMA,
    MAX_GAMMA_LIFT_CHECKED_ROWS, MAX_GAMMA_LIFT_INDEX_SHIFT, MAX_GAMMA_LIFT_OPERATOR_TERMS,
};
mod hadamard_factor;
pub use hadamard_factor::{
    find_hadamard_factor_lifts, HadamardDivisionQuotient, HadamardFactor,
    HadamardFactorCertificate, HadamardFactorIdentityStatus, HadamardFactorLift,
    HadamardFactorLiftOptions, HadamardFactorLiftReport, HadamardFactorRejection,
    HadamardFactorRejectionKind, HadamardFactorTerm, HadamardHiddenCoefficients,
    HADAMARD_FACTOR_LIFTS_SCHEMA, MAX_HADAMARD_BINOMIAL_OFFSET, MAX_HADAMARD_FACTOR_ROWS,
    MAX_HADAMARD_OPERATOR_TERMS,
};
mod lift_summary;
pub use lift_summary::{
    lift_verdict, LiftHypothesis, LiftHypothesisStatus, LiftIdentityStatus, LiftOperation,
    LiftParent, LiftParentProof, LiftRootTransfer, LiftSummary, LiftVerdict,
    LIFT_LEAN_REPLAY_NOT_ATTEMPTED,
};
mod parent_proof;
pub use parent_proof::{
    VeroneseParentAttestation, VeroneseParentProof, VeroneseParentProofAttempt,
    VeroneseParentProofRoute, VeroneseParentProofStatus,
};
mod parent_theorems;
pub use parent_theorems::{
    check_liu_wang_theorem, check_three_term_recurrence_theorem, RealRootednessTheorem,
    RealRootednessTheoremConclusion, RealRootednessTheoremOptions, RealRootednessTheoremReport,
    RealRootednessTheoremStatus, RowInterlacing, REAL_ROOTEDNESS_THEOREM_SCHEMA,
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
