//! Row-factor lifts: recognise the rows of a recurrence-defined sequence as
//! `Q_n = F(x)^(e(n)) S_n` with a fixed polynomial `F` and an affine schedule
//! `e(n) = alpha n + beta` (`alpha >= 0`, `e(n0) >= 0`), and prove the
//! cofactors `S_n` separately ([`find_row_factor_lifts`]).
//!
//! # Detection
//!
//! The first row is factored over `Q`
//! ([`crate::factorization::irreducible_factorization_rational_coeffs_with_budget`]);
//! for each irreducible factor `G` its multiplicity `m(n)` in every detection
//! row is found by exact division, and `e(n) = m(n0) + alpha (n - n0)` takes
//! the largest `alpha >= 0` with `e(n) <= m(n)` on every detection row.
//! Factors with the same schedule are multiplied into one `F`.
//!
//! # The certificate ([`RowFactorCertificate`])
//!
//! With `D o F^w = F^w (D + w F'/F)`, the polynomial operators
//! `T_d = F^d F^(-w) D^d F^w` satisfy `T_0 = 1`,
//! `T_(d+1) = F D T_d + (w - d) F' T_d`.  For the input relation
//! `M = sum c_(i,d)(n, x) D^d S^i` (leader shift `R`, largest derivative
//! order `J`) and `w_i = e(n) + alpha i = e(n + i)`,
//!
//! ```text
//! (M F^e S)_n = F^(e(n) - J) (M_S S)_n,   M_S = sum c_(i,d) F^(alpha i + J - d) T_d(w_i) S^i,
//! ```
//!
//! an identity of rational functions for any rows `S`.  The leader term of
//! `M_S` is `g(n) F^(alpha R + J) S^R`; when every coefficient of `M_S` is
//! divisible by `F^h`, `h = alpha R + J` (checked exactly), `M'_S = M_S / F^h`
//! has polynomial coefficients and the leader `g(n) S^R`, free of `x`.  The
//! cofactor rule `M'_S` with the seeds `S_n = Q_n / F^(e(n))` (exact
//! divisions) generates polynomial rows, since every step divides by the
//! scalar `g(n)`, nonzero at every target index (decided exactly).  Then
//! `Q'_n = F^(e(n)) S_n` satisfies `M` at every base where `M'_S` holds for
//! `S`; at the bases from `n0 + K - R` (where `M` holds for `Q`) below that,
//! `M` is evaluated exactly on `Q'`.  `Q'` agrees with `Q` on the `K`
//! identifying rows (every checked row `Q_n` is divided exactly) and the
//! leader of `M` is a nonzero polynomial at every target, so
//! `Q_n = F^(e(n)) S_n` for every `n >= n0` by induction.
//!
//! # Root transfer
//!
//! A product of polynomials is real-rooted iff every factor is, and its roots
//! are the union of theirs: when `F` is real-rooted (checked exactly, once)
//! and every `S_n` is zero or real-rooted with coefficients of one sign
//! (`prove_parent` on the cofactor rule), every `Q_n` is zero or
//! real-rooted, with roots in `(-inf, 0]` when `F` has them there.  A
//! non-real-rooted `F` fails the route (and every `Q_n` with `e(n) >= 1` has a
//! non-real root).
//!
//! The group is split into `certify` (the substitution and the certificate)
//! and `search` (detection, the report types of
//! `polytool.recurrence-row-factor-lifts.v1` and the entry point).

mod certify;
mod search;

pub use certify::RowFactorCertificate;
pub use search::{
    find_row_factor_lifts, RowFactorIdentityStatus, RowFactorLift, RowFactorLiftOptions,
    RowFactorLiftReport, RowFactorRejection, MAX_ROW_FACTOR_ROWS, ROW_FACTOR_LIFTS_SCHEMA,
};
