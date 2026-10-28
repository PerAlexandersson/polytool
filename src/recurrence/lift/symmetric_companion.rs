//! Symmetric-companion lifts: recognize a palindromic sequence `A_n` as the
//! symmetric part of a parent `P_n = A_n + x B_n` with a first-order
//! differential rule, and transfer real-rootedness through interlacing
//! symmetric decompositions.
//!
//! # The coupled system
//!
//! Fix a unit-step symmetry schedule `d_n = d_(n0) + (n - n0)` and write
//! "`f` is symmetric at `d`" for `deg f <= d` and `x^d f(1/x) = f` (an ambient
//! index: degree deficits and the zero polynomial are allowed).  For a parent
//! rule
//!
//! ```text
//! P_(n+1) = (a_n + b_n x) P_n + c_n x (1 - x) D P_n,     D = d/dx,
//! ```
//!
//! put `delta_n = b_n - a_n - c_n d_n`, `epsilon_n = c_n - delta_n`,
//! `H(d, a, c) f = [a + (a + c d) x] f + c x (1 - x) D f`,
//! `U_n = H(d_n, a_n, c_n)` and `V_n = H(d_n - 1, a_n + delta_n, c_n)`.  The
//! symmetric decomposition `P_n = A_n + x B_n` (`A_n` symmetric at `d_n`,
//! `B_n` at `d_n - 1`) then evolves by
//!
//! ```text
//! A_(n+1) = U_n A_n + epsilon_n x B_n,
//! B_(n+1) = delta_n A_n + V_n B_n.
//! ```
//!
//! On monomials `H(d, a, c) x^j = (a + c j) x^j + (a + c (d - j)) x^(j+1)`, so
//! `H(d, a, c)` maps symmetric-at-`d` to symmetric-at-`(d + 1)`, and
//! multiplication by `x` maps symmetric-at-`d` to symmetric-at-`(d + 2)`.
//!
//! # Inverse recognition
//!
//! Given only a scalar recurrence for `A_n` and seeds, the companion is
//! forced once `epsilon_n != 0`:
//! `R_n = A_(n+1) - U_n A_n`, `B_n = R_n / (epsilon_n x)`.  With
//! `Vtilde_n = V_n - c_n (1 - x)` (so `Vtilde_n o x = x o V_n`), eliminating
//! `B` gives the operator
//!
//! ```text
//! E = epsilon_n S^2 - [epsilon_n U_(n+1) + epsilon_(n+1) Vtilde_n] S
//!     + [epsilon_(n+1) Vtilde_n o U_n - epsilon_n epsilon_(n+1) delta_n x],
//! ```
//!
//! and, for arbitrary rows, with `E_n = A_(n+1) - U_n A_n - epsilon_n x B_n`
//! and `F_n = B_(n+1) - delta_n A_n - V_n B_n`,
//!
//! ```text
//! (E A)_n = epsilon_n E_(n+1) + epsilon_n epsilon_(n+1) x F_n - epsilon_(n+1) Vtilde_n E_n.   (*)
//! ```
//!
//! # The search ([`find_symmetric_companion_lifts`])
//!
//! The ansatz takes `a_n`, `c_n`, `epsilon_n` polynomial in `n` of degree at
//! most [`SymmetricCompanionLiftOptions::max_parameter_degree_in_n`] (0 or 1;
//! `delta_n = c_n - epsilon_n`, and `b_n = a_n + c_n d_n + delta_n` is then of
//! degree at most two), with exact rational coefficients, and works on
//! `fit_rows` rows:
//!
//! 1. the schedule is supplied or inferred from the first nonzero row;
//! 2. `R_n(0) = 0` forces `a_n = A_(n+1)(0) / A_n(0)` where `A_n(0) != 0`; these
//!    ratios are interpolated exactly (if every constant term vanishes, `a_n`
//!    becomes an unknown instead);
//! 3. substituting `B` into the second equation and multiplying by
//!    `epsilon_n epsilon_(n+1) x` gives `G_n = 0`, where (`eps` for `epsilon`)
//!    `G_n = eps_n R_(n+1) - eps_(n+1) Vtilde_n R_n - eps_n eps_(n+1) delta_n x A_n`
//!    is a polynomial of degree at most three in the unknown coefficients in
//!    which every monomial carries an `epsilon` coefficient.  Its coefficients in `x` over all searched rows
//!    form a homogeneous linear system in those monomials (at most
//!    [`SymmetricCompanionLiftOptions::max_linearized_monomials`]); its exact
//!    kernel is computed.  A one-dimensional kernel `k` is the line through
//!    the monomial vector of a solution with `epsilon != 0`: the scale `t`
//!    follows from `t k[e]^2 = k[e^2]` for an `epsilon` coefficient `e`, every
//!    coordinate from `k[v e] / k[e]`, and the whole vector is compared with
//!    `t k`.  Solutions with `epsilon = 0` satisfy `G = 0` trivially but map
//!    to the zero monomial vector;
//! 4. otherwise (a larger kernel, or a failed recovery) the independent
//!    equations, together with `z epsilon(n0) = 1`, are solved exactly: a lex
//!    Gröbner basis (Buchberger), the rational roots of its univariate
//!    element ([`crate::factorization`]), substitution and recursion, within
//!    [`SymmetricCompanionLiftOptions::max_groebner_steps`].  A zero kernel or
//!    no rational point means no solution with `epsilon(n0) != 0` in the
//!    ansatz; infinitely many solutions or a budget stop are reported as
//!    `underdetermined` or `budget_exceeded` (`exhausted_budget`).
//!
//! Degrees are tried in increasing order; the search stops after the first
//! degree with a certified lift.  Only rational parameter points are found,
//! and every candidate is only a candidate until certified.
//!
//! # The certificate ([`verify_symmetric_companion_certificate`])
//!
//! A [`SymmetricCompanionCertificate`] records the source, `d_(n0)`, the
//! parameters, `E`, a derivation of `E` from the source rule `M` and the
//! parent.  The verifier recomputes every derived item and establishes, for
//! all `n >= n0`:
//!
//! 1. **`E` annihilates `A`.**  An exact Ore-algebra identity
//!    `g(n) E = sum_nu H_nu S^nu M` ([`ReconstructionCertificate`]): right
//!    division by `M` when its leader is a constant (`g = 1`), or, for an
//!    order-two source with leader `f(n)` free of `x`, `f(n + 2) E =
//!    epsilon_n M`.  `M` holds at every base whose leading row is generated
//!    by the source rule; the finitely many remaining bases, and roots of
//!    `g`, are checked on exact rows.
//! 2. **`epsilon_n != 0`**, decided exactly ([`crate::integer_range_sign`]);
//!    a vanishing coupling is reported as unsupported with its index, never
//!    cancelled.
//! 3. **`R_n(0) = 0`.**  The derivative terms of `E` vanish at `x = 0`, and
//!    its constant-term projection factors exactly as
//!    `(epsilon_n S - g_n) o (S - a_n)`, `g_n = epsilon_(n+1) (a_n - epsilon_n)`.
//!    So `w_n = A_(n+1)(0) - a_n A_n(0)` satisfies `epsilon_n w_(n+1) = g_n w_n`;
//!    with `w_(n0) = 0` (exact) and step 2, `w_n = 0` for all `n`.  Hence
//!    `B_n = R_n / (epsilon_n x)` is a polynomial and `E_n = 0` by definition.
//! 4. **The second equation.**  Identity `(*)` is checked as two exact operator
//!    identities (its `A` component, and the `B` component
//!    `Vtilde_n o x = x o V_n`).  With `E_n = E_(n+1) = 0` and `(E A)_n = 0` it
//!    gives `epsilon_n epsilon_(n+1) x F_n = 0`, so `F_n = 0`.
//! 5. **Symmetry**, by induction on pairs: the seeds `A_(n0)`, `A_(n0+1)` are
//!    symmetric at `d_(n0)`, `d_(n0) + 1` (exact) and `d_(n0) >= 1`.  If
//!    `A_n`, `A_(n+1)` are symmetric at `d_n`, `d_n + 1`, then `R_n` is
//!    symmetric at `d_n + 1` with `R_n(0) = 0`, so its top coefficient
//!    vanishes and `B_n` is symmetric at `d_n - 1`; the second equation makes
//!    `B_(n+1)` symmetric at `d_n`, and the first (at `n + 1`) makes `A_(n+2)`
//!    symmetric at `d_n + 2`.
//! 6. **Collapse.**  The operator identities `L_n = U_n + delta_n x` and
//!    `L_n o x = epsilon_n x + x o V_n` (`L_n` the parent operator) give
//!    `P_(n+1) = L_n P_n`; the parent and its seed are compared with the
//!    record.
//!
//! A row self-check on `self_check_rows` indices is a defensive replay, not
//! part of the argument.  The identity status is `certified` only when 1--6
//! hold; `candidate` (some step unresolved, for example no derivation of `E`),
//! `unsupported` (vanishing coupling, `d_(n0) = 0`) and `rejected` (a step
//! fails exactly) are kept apart.
//!
//! # The real-rootedness route
//!
//! **Preservation theorem.**  For `d >= 1` and reals `r, c, eta` with
//! `r >= 0` and `0 <= eta <= c` (`epsilon = c - eta`), the step
//! `(a, b) -> (H(d, r, c) a + epsilon x b, eta a + H(d - 1, r + eta, c) b)`
//! maps every pair of polynomials with nonnegative coefficients,
//! `deg a <= d`, `deg b <= d - 1` and `b <= a` to a pair with the same
//! properties for `d + 1`; the conditions are also
//! necessary for preserving all such symmetric pairs.  Proof sketch: for
//! strictly interlacing simple negative roots, `R = b / a` has the partial
//! fraction expansion `l + sum t_i / (x - rho_i)` with `l >= 0`, `t_i > 0`,
//! `rho_i < 0`, so `R' <= 0` and `R + x R' >= 0`; at a root of `a+` the sign of
//! `b+ / a` is that of
//! `eta + (eta - epsilon x) R - epsilon x R^2 + c x (1 - x) R'`, which these
//! inequalities make positive, so `b+` alternates on the roots of `a+`; the
//! weak and zero cases follow by coefficientwise limits.  By the
//! interlacing-decomposition theorem of Brändén--Solus (arXiv:1808.04141,
//! Theorem 2.6), `p = a + x b` is then real-rooted.  Ma--Ma--Yeh--Yeh
//! (arXiv:2102.00899) give the symmetric-part recurrences of the worked
//! `1/k`-Eulerian example.
//!
//! The route applies when the identity is certified, `a_n >= 0` and
//! `0 <= delta_n <= c_n` for every integer `n >= n0` (decided exactly by
//! [`crate::integer_range_sign`]; a failure names the smallest failing
//! index), `d_(n0) >= 1`,
//! `A_(n0)`, `B_(n0)` have nonnegative coefficients and `B_(n0) <= A_(n0)`.
//! Induction then shows that every `A_n`, `B_n`, `P_n` is zero or real-rooted.
//! Real-rootedness of the parent alone is never accepted for `A_n`, and a
//! failed hypothesis only means that this route fails.  The Lean replay status
//! is always `not_attempted`.
//!
//! **Interlacing convention.**  `f <= g` means: `f, g` real-rooted with
//! `deg g` in `{deg f, deg f + 1}`, roots weakly alternating, `g` holding the
//! rightmost root; `0 <= g` and `g <= 0` when `g` is real-rooted, and
//! `(0, 0)`.  For nonzero rows this is exactly
//! [`crate::check_weak_interlacing_bigint_coeffs`]`(f, g) == Some(true)`
//! (`f << g`, `f` on the left; the same-degree case requires every root of
//! `f` to lie weakly left of the matching root of `g`); the zero cases are
//! handled separately because that function answers `Some(false)` for
//! `f << 0`.
//!
//! Related: [`RecurrenceJson`] (source and parent records),
//! [`find_reconstruction_certificate`] (the right division of step 1),
//! [`crate::sequences::literature::one_over_k_eulerian_symmetric_parts_rational`]
//! (the worked family).

mod certificate;
mod groebner;
mod operators;
mod report;
mod search;
mod source;

pub use certificate::{
    verify_symmetric_companion_certificate, SymmetricCompanionCertificate, SymmetricCompanionCheck,
    SymmetricCompanionCheckStatus, SymmetricCompanionIdentityStatus, SymmetricCompanionParameters,
    SymmetricCompanionRealRootednessStatus, SymmetricCompanionVerification,
    SYMMETRIC_COMPANION_LEAN_REPLAY_NOT_ATTEMPTED,
};
pub use report::{
    SymmetricCompanionAnsatzAttempt, SymmetricCompanionAttemptOutcome, SymmetricCompanionLift,
    SymmetricCompanionLiftReport, SymmetricCompanionSearchSummary, SymmetryScheduleOrigin,
    SymmetryScheduleRecord, SYMMETRIC_COMPANION_LIFT_SCHEMA,
};
pub use search::{
    find_symmetric_companion_lifts, SymmetricCompanionLiftOptions,
    MAX_SYMMETRIC_COMPANION_PARAMETER_DEGREE,
};

#[cfg(test)]
mod tests;
