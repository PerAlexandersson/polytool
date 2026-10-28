//! Hadamard-factor lifts: recognise the rows `Q_n = sum_k q(n, k) x^k` of a
//! recurrence-defined sequence as `Q_n = lambda(n) sum_k u(n, k) p(n, k) x^k`
//! with a known factor `u` (binomials, factorials, rising factorials and
//! products of two) and a parent `P_n = sum_k p(n, k) x^k` whose rule is
//! smaller than the input's ([`find_hadamard_factor_lifts`]).
//!
//! # Detection
//!
//! For each factor of the menu ([`HadamardFactorLiftOptions`]) the quotients
//! `p = q / u` are computed exactly on the detection rows; the factor is
//! rejected when `u` is zero or undefined at a nonzero `q(n, k)`.  Where `u`
//! vanishes, `p` is not determined by the input; it is set to zero, or
//! continued by a ratio `p(n, k + 1) / p(n, k)` linear over linear in
//! `(n, k)` fitted exactly to every row ([`HadamardHiddenCoefficients`]); the
//! certificate below does not depend on these values.  Optionally the rows
//! are normalised by their leading coefficients `lambda(n)` when
//! `lambda(n) / lambda(n - 1)` is an exactly fitted rational function.  A
//! parent rule is then searched with
//! [`crate::recurrence::find_recurrence_adaptive_rational_with_budget`]
//! within the input's order and degrees (denominators free of `x`) and kept
//! only when it is smaller: lower order, or the same order with degrees in
//! `x` and `n` no larger and one of them smaller.
//!
//! # The conjugation
//!
//! Write the input relation `M` (base index `n`, leader shift `R`) in
//! Laurent--Euler form `sum x^s m_(i,s)(n, theta) S^i` (`theta = x D`, `s`
//! may be negative since `D^d = x^(-d) theta (theta - 1) ... (theta - d +
//! 1)`).  For any rows `P` and `Q'_n = lambda(n) u ⊙ P_n` the coefficient of
//! `x^k` in `(M Q')_n` is
//!
//! ```text
//! E(n, k) = sum_(i,s) m_(i,s)(n, k - s) U(n + i, k - s) p(n + i, k - s),   U = lambda u.
//! ```
//!
//! Every factor is a product of Gamma atoms `Gamma(n + c)^(+-1)`,
//! `Gamma(k + c)^(+-1)` and `Gamma(n - k + c)^(-1)` (with `1 / Gamma` zero
//! at the poles), so `rho_(i,s) = U(n + i, k - s) / U(n + R, k)` is a product
//! of Pochhammer symbols, linear factors over linear factors in `(n, k)`
//! (the gauge contributes `prod b(n + j) / a(n + j)`).  Let `D(n, k)` be the
//! least common multiple of their denominators and `M''` the operator with
//! the terms `x^(s + J'') e_(i,s)(n, theta + s) m_(i,s)(n, theta) S^i`,
//! `e = D rho`, where `J''` clears negative powers of `x`.  The
//! **conjugation lemma**: at every `n >= n0` and `k >= 0`,
//!
//! - if `U(n + R, k) != 0`, then `D(n, k) != 0` and for every term with
//!   `k - s >= 0` the value `U(n + i, k - s)` equals `rho_(i,s)(n, k) U(n + R, k)`
//!   (zero included), so `E(n, k) = U(n + R, k) / D(n, k) [x^(k + J'')](M'' P)_n`;
//! - if `U(n + R, k) = 0`, then `U(n + i, k - s) = 0` for every term with
//!   `k - s >= 0`, so `E(n, k) = 0`.
//!
//! It holds atom by atom under these conditions (checked exactly, else the
//! identity is only [`HadamardFactorIdentityStatus::SupportedByRows`]):
//! `Gamma(n + c)`: `n0 + c >= 1`; `Gamma(k + c)^(-1)` with integer `c <= 0`
//! (zero for `k <= -c`): every `s >= 0`; `Gamma(k + c)` with integer `c`:
//! `c >= 1` and every `s <= c - 1`; non-integer `c`: no condition (no zeros
//! or poles); `Gamma(n - k + c)^(-1)` (zero for `k >= n + c`): every
//! `i + s <= R`; `Gamma(n - k + c)` (poles above the support): not
//! supported.  For example `binom(n + a, k + b)` needs `n0 + a >= 0`, every
//! `i + s <= R`, and `s >= 0` when `b < 0`.  The proof is the case analysis
//! of `Gamma(l + j) / Gamma(l)` products: where the reference value is
//! nonzero every denominator factor is a nonzero integer or non-integer
//! argument, and a zero shifted value makes one numerator Pochhammer factor
//! run through zero.  Products of atoms inherit the lemma (all atoms are
//! finite on the domain).  The unit test `atom_ratios_match_values` checks
//! the lemma against exact values.  Hidden quotient coefficients never
//! matter: they are multiplied by zero in `Q'` and by a vanishing
//! Pochhammer factor in `M''`.
//!
//! # The certificate ([`HadamardFactorCertificate`])
//!
//! `M''` is pseudo-divided by the parent relation `L` (leader `g(n) S^k`,
//! free of `x` and `D`): `w(n) M'' = sum_nu H_nu S^nu L + r`.  When `r != 0`
//! and the parent is a product parent `g(n) P_(n+1) = a(n) c(x) P_n` with one
//! seed `f`, every parent row `Lambda(n) c^(n - n0) f` is annihilated by
//! `A(n) = c f D - ((n - n0) c' f + c f')` (exactly: `A(n0) f = 0` and
//! `A(n + 1) o c = c o A(n)`), and `r` is right-divided by `A`:
//! `mu(x) r = Z o A` with `mu` a power of `c f`.  So
//!
//! ```text
//! mu(x) w(n) M'' = sum_nu mu H_nu S^nu L + Z A.
//! ```
//!
//! Applied to `P` this gives `(M'' P)_n = 0` at every base `n` where `L`
//! holds for the parent (`n >= n0 + Kp - k`) other than the integer roots of
//! `w`.  At those roots and at the bases from `n0 + K - R` (where `M` holds
//! for `Q`) below the parent's, `M` is evaluated exactly on the candidate
//! rows.  By the lemma `Q' = lambda u ⊙ P` satisfies `M` at every base
//! `n >= n0 + K - R`.  `Q'` agrees with `Q` on the `K` identifying rows (all
//! compared rows are checked exactly), and the leader of `M` is a nonzero
//! polynomial at every target index (decided exactly when the input is
//! validated), so `Q' - Q` vanishes by induction:
//! `Q_n = lambda(n) u ⊙ P_n` for every `n >= n0`
//! ([`HadamardFactorIdentityStatus::Certified`]).
//! [`HadamardFactorCertificate::verify`] replays the algebra from the
//! recorded input, factor, gauge and parent.
//!
//! # Root transfer
//!
//! The Pólya--Schur characterization of multiplier sequences (Pólya and
//! Schur 1914; Theorem 40 of Brändén, *Unimodality, log-concavity,
//! real-rootedness and beyond*, arXiv:1410.6601): `x^k -> lambda_k x^k`
//! preserves real-rootedness iff `sum_k binom(m, k) lambda_k x^k` is
//! real-rooted with zeros of one sign for every `m`, iff
//! `sum_k lambda_k x^k / k!` (or its reflection) is `C x^m e^(a x) prod (1 +
//! alpha_j x)` with `a, alpha_j >= 0`.  For each `n >= n0`, `k -> u(n, k)` is
//! a nonnegative multiplier sequence for
//!
//! - `binom(n + a, k + b)`: `sum_k u x^k / k! = L_(A - b)^(b)(-x)` for
//!   `b >= 0` and `x^c A! / (A + c)! L_A^(c)(-x)` for `b = -c < 0`
//!   (`A = n + a`), generalized Laguerre polynomials with parameter `> -1`,
//!   whose zeros are real and positive (Szegő, *Orthogonal Polynomials*,
//!   Theorem 3.3.1);
//! - `1/k!` (Brändén, Example 41) and `1/(alpha)_k`, `alpha > 0`:
//!   `sum_k binom(m, k) x^k / (alpha)_k = m! / (alpha)_m L_m^(alpha - 1)(-x)`;
//! - `1/(n - k)!`: `sum_k u x^k / k! = (1 + x)^n / n!`;
//! - `k! binom(n + a, k - c)`, `c >= 0`: `sum_k u x^k / k! = x^c (1 + x)^(n + a)`;
//! - products of these (composition of multiplier sequences).
//!
//! A nonnegative multiplier sequence maps a polynomial that is zero or
//! real-rooted with coefficients of one sign to one of the same kind, and
//! `lambda(n) != 0` changes nothing, so every `Q_n` is zero or real-rooted
//! with roots in `(-inf, 0]` once every parent row is (the statement of the
//! parent provers, called through `prove_parent`).  `k!`, `(n - k)!` and
//! `(alpha)_k` alone are not multiplier sequences (`sum x^k = 1/(1 - x)` is
//! not entire); such lifts are reported with an unresolved
//! `factor_multiplier_sequence` hypothesis (open).
//!
//! The group is split into `factor` (menu, values, Gamma atoms, the lemma's
//! conditions, ratios, justifications), `certify` (conjugation, divisions,
//! the certificate), `search` (detection and the entry point), `report`
//! (the `polytool.recurrence-hadamard-factor-lifts.v1` types and the
//! [`crate::recurrence::LiftSummary`] mapping) and `prepare` (input
//! validation shared with the row-factor lifts).

mod certify;
mod factor;
pub(super) mod prepare;
mod report;
mod search;

pub use certify::{
    HadamardDivisionQuotient, HadamardFactorCertificate, MAX_HADAMARD_OPERATOR_TERMS,
};
pub use factor::{HadamardFactor, HadamardFactorTerm};
pub use report::{
    HadamardFactorIdentityStatus, HadamardFactorLift, HadamardFactorLiftReport,
    HadamardFactorRejection, HadamardFactorRejectionKind, HadamardHiddenCoefficients,
    HADAMARD_FACTOR_LIFTS_SCHEMA,
};
pub use search::{
    find_hadamard_factor_lifts, HadamardFactorLiftOptions, MAX_HADAMARD_BINOMIAL_OFFSET,
    MAX_HADAMARD_FACTOR_ROWS,
};
