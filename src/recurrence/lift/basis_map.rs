//! Basis maps `T_P : u^k -> P_k(t)` of recurrence-defined families: the map
//! as a value ([`BasisMap`]), the theorems that make it a real-rootedness
//! preserver on a source interval ([`check_basis_map_preserves_pf`]), the
//! backward lift `Q_n = T_P(R_n)` ([`find_basis_map_lifts`],
//! [`certify_basis_map_lift`]) and the forward move `R -> T_P(R)`
//! ([`apply_basis_map_to_family`]).
//!
//! Throughout, PF means a real polynomial with nonnegative coefficients and
//! only real nonpositive roots (equivalently, with a positive leading
//! coefficient, a nonzero polynomial whose roots all lie in `(-inf, 0]`);
//! the zero polynomial is allowed as an output.  `R(S)` is the set of
//! nonzero real polynomials whose roots all lie in `S` (constants included).
//!
//! # The map
//!
//! A [`BasisMap`] holds a family `P` by its rule `f(n, t) P_n = sum_j
//! c_j(n, t) D^(d_j) P_(n - r_j)`, its first index `i0` and seed rows, and
//! sends `u^k` to the row `P_(i0 + k)`:
//! `T_P(sum_k p_k u^k) = sum_k p_k P_(i0 + k)`.  [`BasisMap::apply`] is this
//! sum; [`BasisMap::expand`] is the triangular inverse: when the rows used
//! are nonzero with strictly increasing degrees, a polynomial `q` lies in
//! their span exactly when top-down elimination leaves no remainder, and
//! the coefficients are then unique (the result is re-applied as a check).
//!
//! # Intertwining
//!
//! For a monomial `u^k` of `u^R p` (`R` the order) the rule at
//! `n = i0 + k` reads, since `T(u^k) = P_(i0 + k)`,
//!
//! ```text
//! sum_i f_i(n) t^i T(u^k) = sum_j sum_i c_(j,i)(n) t^i D^(d_j) T(u^(k - r_j)),
//! ```
//!
//! and `n = i0 + k` is `theta + i0` on `u^k` and `theta + i0 + r_j` on
//! `u^(k - r_j)` (`theta = u d/du`).  Summing over the monomials gives the
//! [`BasisMapIntertwining`] relation
//! `sum t^i D^d T(phi(theta) u^m p) = 0`, valid for every `p` divisible by
//! `u^(K - R)` (the rule holds from the first computed index on, `K` seed
//! rows).  For a rule free of `n` with a constant leader and `K = R` it is
//! the companion form `T(u^R p) = sum_j A_j T(u^j p)` with fixed operators
//! `A_j` in `(t, D_t)`: with the state `V(p) = (T(p), T(u p), ...,
//! T(u^(R-1) p))`, `V(u p) = C V(p)` for the companion operator matrix `C`
//! (`C[i][i+1] = 1`, last row `A_0, ..., A_(R-1)`), hence
//! `V(e(u) p) = e(C) V(p)` for every polynomial `e`.  For `R = 1` this is
//! `T o u = A o T`.
//!
//! # Theorems
//!
//! [`check_basis_map_preserves_pf`] recognises the following shapes (after
//! dividing by a constant leader and a constant first row `C` where there
//! is one) and decides every hypothesis exactly (rational signs, exact root
//! isolation of fixed polynomials, root multiplicities).  Each proved
//! [`BasisMapConclusion`] states: every polynomial in `R(S)` is sent to zero
//! or into `R(J)`, for all input degrees.
//!
//! - **Geometric** `P_k = g (a + d t)^k`: `T(p) = g(t) p(a + d t)`, so the
//!   roots of `T(p)` are those of `g` and the points `(r - a)/d` (elementary).
//!   With `g` real-rooted in `(-inf, 0]`: `d > 0` gives `(-inf, a] ->
//!   (-inf, 0]` (sharp: a linear input with root `> a` fails), `d < 0` gives
//!   `[a, inf) -> (-inf, 0]`, `d = 0` gives `T(p) = p(a) g`.  A
//!   non-real-rooted `g = T(1)` is a counterexample.
//! - **Eventually geometric** `P_k = P_K (a + d t)^(k - K)` from `K >= 1` on,
//!   `d != 0`, `P_K != 0`, with a nonzero prefix discrepancy: for every
//!   `r != 0` some `T((u + r)^N)` in an explicit window is not real-rooted,
//!   so no source interval of positive length works (RealRooted
//!   `exists_not_splits_linearMap_shiftedPower_of_eventually_geometric`);
//!   the window is searched exactly for `r = 1` and the failing input is
//!   reported.
//! - **Constant current** `P_0 = 1, P_1 = s, P_(k+2) = a P_(k+1) + (b t +
//!   c) P_k`, `a, s >= 0`, `b > 0`: `P_k(t) = Q_k(b t + c)` for the rows with
//!   `b = 1, c = 0`, whose map sends PF inputs to PF outputs (RealRooted
//!   `basisTransform_orderTwoBasis_eq_comp`,
//!   `orderTwoZeroTransform_isPFPolynomial`; an interlacing-preserving 2 x 2
//!   matrix, Brändén, arXiv:1410.6601, Theorem 8.5), so every PF input has a
//!   real-rooted image with roots `<= -c/b`; and PF inputs (also with roots
//!   in any `[-R, 0]`) have PF images if and only if `a^2 + 4c >= 0` and
//!   (`2s >= a` or `c >= s^2 - a s`) (RealRooted
//!   `orderTwoBasis_pf_equivalences`,
//!   `basisTransform_orderTwoBasis_preserves_bounded_pf_iff`).  When the
//!   criterion fails, the first monomial `u^n` with `P_n(0) < 0` is
//!   reported.
//! - **Affine current** (below).
//! - **Ray differential** `P_(k+1) = (a + b t) P_k + kappa (t - h) P_k'`,
//!   `b > 0`, `kappa >= 0`, `g = P_0` real-rooted in `(-inf, h]` with
//!   multiplicity `m` at `h`: `R(S) -> R((-inf, h])` if and only if `S` lies
//!   in `(-inf, a + b h + kappa m]` (RealRooted
//!   `rayBasisMap_maps_realRootedIn_iff`; the step is Ma and Wang, EJC 15
//!   (2008) R17, Theorem 2).
//! - **Interval differential** `P_(k+1) = (a + b t) P_k + kappa (t - l)(t -
//!   h) P_k'`, `l < h` rational, `kappa > 0`, `b + kappa deg g > 0`, `g`
//!   real-rooted in `[l, h]` with endpoint multiplicities `m_l, m_h`:
//!   `R(S) -> R([l, h])` if and only if `S` lies in `[a + b l - kappa (h -
//!   l) m_l, a + b h + kappa (h - l) m_h]` (RealRooted
//!   `intervalBasisMap_maps_realRootedIn_iff`; Ma and Wang, Corollary 3).
//!   A negative leading coefficient of `g` only changes the sign of `T`.
//! - **Hermite type** `P_(k+1) = (a + b t) P_k - kappa P_k'`, `b, kappa >
//!   0`, `g` real-rooted: `T(p) = c prod_j (L - rho_j) g` for `p = c
//!   prod_j (u - rho_j)` (the factors commute), and each `L - rho` keeps a
//!   real-rooted polynomial real-rooted because its derivative coefficient
//!   `-kappa` is negative at every root (Ma and Wang, Theorem 2): every
//!   real-rooted input has a real-rooted image.
//! - **Index-dependent scalar term** `P_(k+1) = (L + k s) P_k`, `L f = (a +
//!   b t) f + kappa (t - h) f'`, `s >= 0`: `T = U_g o R_s` (RealRooted
//!   `basisMap_factorization_of_scalarIndex_recurrence`) with `R_s(u^k) =
//!   u (u + s) ... (u + (k - 1) s)`, which keeps real-rootedness of inputs
//!   with roots `<= 0` (RealRooted `NewtonParametersAffine.preservesNonposRooted`;
//!   Su, Yang and Zhang, Australas. J. Combin. 56 (2013), Theorem 4; Brenti,
//!   Mem. AMS 413 (1989), Theorem 2.4.3) and nonnegative coefficients, hence
//!   roots `<= 0`; then the ray theorem for `U_g(u^k) = L^k g` needs `a + b
//!   h + kappa m >= 0` (for `kappa = 0` the root `h` is free and the
//!   smallest `h >= 0` with `a + b h >= 0` is used).  For a Newton basis
//!   (`kappa = 0`, constant `g`, `b > 0`) whose scalar is not affine with a
//!   nonnegative step the map does not preserve inputs with roots `<= 0`
//!   (RealRooted `newtonPreservesNonposRooted_iff_parametersAffine`); a
//!   failing input `(u + r)^n` is searched exactly and reported.
//!
//! ## The affine-current theorem
//!
//! Let `P_0 = 1`, `P_1 = l + m t` and `P_(k+2) = (a + d t) P_(k+1) + (b t +
//! c) P_k` with `d > 0`.  With `y = a + d t`, `beta = b/d`, `K = beta^2 +
//! a beta - c`, `eta = m/d` and `s = l - eta a` the rows are `Q_0 = 1`,
//! `Q_1 = eta y + s`, `Q_(k+2) = y Q_(k+1) + (beta y + beta^2 - K) Q_k`.
//! Assume `K > 0` and `eta >= 0`; put `E = [-beta - sqrt(K), -beta +
//! sqrt(K)]` (a set of input roots) and `I_y = [-2 beta - 2 sqrt(K), -2 beta
//! + 2 sqrt(K)]`.
//!
//! *Common interlacer.*  Let `lambda_+, lambda_-` be the roots of
//! `lambda^2 - y lambda - (beta y + beta^2 - K)`.  For an input `p` put
//! `V_p = (p(lambda_+) - p(lambda_-)) / (lambda_+ - lambda_-)` and `F_p =
//! (lambda_+ p(lambda_+) - lambda_- p(lambda_-)) / (lambda_+ - lambda_-)`,
//! polynomials in `y` (symmetric in the roots).  Both satisfy the row
//! recurrence on monomials, with `V_1 = 0, V_u = 1, F_1 = 1, F_u = y`, so
//! `T(p) = F_p + ((eta - 1) y + s) V_p`.  On `I_y` write `y = 2 sqrt(K)
//! cos(theta) - 2 beta`, so `lambda_+ = sqrt(K) e^(i theta) - beta`,
//! `0 < theta < pi`.  For a monic input `p = prod_j (u + r_j)` with every
//! `-r_j` in the interior of `E` (`|r_j - beta| < sqrt(K)`), each factor
//! `sqrt(K) e^(i theta) + r_j - beta` has an argument that increases
//! strictly from `0` to `pi` (its derivative is `sqrt(K)(sqrt(K) + (r_j -
//! beta) cos theta) / |.|^2 > 0`), so the argument of `p(lambda_+)` increases
//! strictly from `0` to `n pi`.  Since `V_p = Im p(lambda_+) / Im lambda_+`,
//! its zeros `xi_1 < ... < xi_(n-1)` are where the argument is `k pi`; there
//! `T(p) = F_p = p(lambda_+)` has sign `(-1)^k`, and because `y` decreases
//! in `theta`, `T(p)(xi_(n-1)) < 0` and the signs alternate along the real
//! axis.  `T(p)` has degree `n` and leading coefficient `eta > 0` (for
//! `eta > 0`), so the intermediate value theorem puts one root in each of
//! the `n` gaps `(-inf, xi_1), ..., (xi_(n-1), inf)`: `T(p)` is real-rooted
//! and strictly interlaced by `V_p`, whose roots lie in the interior of
//! `I_y`.  Inputs with roots on the boundary of `E`, and `eta = 0`, follow
//! by coefficient limits (a limit of real-rooted polynomials of bounded
//! degree is zero or real-rooted).  This gives the conclusion `E -> R`
//! (real-rooted images).
//!
//! *Target interval.*  At the double characteristic roots `lambda_(+-)^0 =
//! -beta +- sqrt(K)`, i.e. at `y_(+-) = -2 beta +- 2 sqrt(K)`, continuity
//! gives `T(p)(y_(+-)) = p(lambda^0) + (eta y_(+-) + s - lambda^0)
//! p'(lambda^0)`.  For `eta > 0`, the largest root of `T(p)` is at most
//! `y_+` when `T(p)(y_+) >= 0` (one root right of `xi_(n-1)`), which holds
//! when `eta y_+ + s >= lambda_+^0` (then `p(lambda_+^0), p'(lambda_+^0) >
//! 0`); symmetrically the smallest root is at least `y_-` when `eta y_- + s
//! <= lambda_-^0`.  With `e = 2 eta - 1` the two conditions are
//! `|s - e beta| <= e sqrt(K)`, checked as `e >= 0` and
//! `(s - e beta)^2 <= e^2 K`.  Then `E -> I`, `I` the image of `I_y` in
//! `t`.
//!
//! *PF inputs.*  Assume moreover `a >= 0`, `K >= beta^2` (so `0` lies in
//! `E`), `R = beta + sqrt(K) > 0`, `a + 2 beta >= 0` and `a^2 + 4c >= 0`
//! (together: `y_+ <= a`, the image of `t = 0`), `l >= 0` and the scalar
//! criterion `2l >= a` or `c >= l^2 - a l`, which is equivalent to `v_k =
//! P_k(0) >= 0` for every `k` (`v_0 = 1, v_1 = l, v_(k+2) = a v_(k+1) + c
//! v_k`; RealRooted `orderTwoConstants_nonneg_iff`).  For a monic input with
//! nonnegative coefficients and roots in the interior of `[-R, 0]`, `T(p)`
//! at `y = a` is `sum_k p_k v_k >= 0`, and `a >= y_+ > xi_(n-1)`, so the
//! single root right of `xi_(n-1)` is at most `a`: every root of `T(p)` is
//! `<= 0` in `t`, and the leading coefficient is positive, so `T(p)` is PF.
//! Limits give the closed interval and `m = 0`.  This is the conclusion
//! `[-R, 0] -> (-inf, 0]` (intersected with `I` when the target condition
//! holds).  The Chebyshev case is Iserles and Saff, Math. Proc. Cambridge
//! Philos. Soc. 105 (1989), Theorem 1; polytool does not claim the
//! sharpness of `E` (a separate argument).
//!
//! # Transport
//!
//! Let the coefficient family `R` have the relation `sum_i e_i(n, u)
//! R_(n+i) = 0` (no derivatives in `u`, leader `g(n) S^k` free of `u`),
//! valid from the base `b_R`.  Applying `V` and `V(e(u) p) = e(C) V(p)`
//! gives `Phi W = 0` for the states `W_n = V(R_n)`, where
//! `Phi = sum_i e_i(n, C) S^i` is an operator matrix in `(n, t, D_t, S)`
//! with leader `g(n) S^k I`; the image rows are `H_n = T(R_n) = (W_n)_0`.
//! For a map of order one `Phi` is a scalar rule for `H` (the transported
//! rule).
//!
//! # The lift certificate
//!
//! The input `Q` has the relation `M` (from the base `b_M`, `K` identifying
//! rows, a leader nonzero at every target, so the identifying rows
//! determine `Q`).  [`certify_basis_map_lift`] checks the first rows `H_n =
//! Q_n` exactly, then pseudo-divides the row `e_0^T M` by `Phi` (the
//! division of the Veronese certificates): `w(n) e_0^T M = sum_nu H_nu S^nu
//! Phi` with remainder zero.  Applied to `W` this gives `w(n) (M H)_n = 0`
//! at every base `n >= b_R`; at the finitely many bases `n >= b_M` below
//! `b_R` or with `w(n) = 0` the relation `M` is evaluated exactly on the
//! candidate rows.  So `H` satisfies `M` from `b_M` on, agrees with `Q` on
//! the identifying rows, and `H_n = Q_n` for every `n >= n0` by induction.
//! For a map of order one the reverse division (`w Phi = sum H_nu S^nu M`,
//! when `M` has a leader free of `t`) is tried as well: then `Q` satisfies
//! `Phi`, which with the identifying rows of `R` gives the same identity.
//! [`BasisMapDivisionCertificate::verify`] replays all of it from the
//! recorded rules and seeds.
//!
//! # Root location
//!
//! The coefficient family must lie in the source `S` of a conclusion of the
//! map theorem.  Rows that violate this (checked exactly on the compared
//! rows) make the route fail.  Otherwise rational bounds `L' >= inf S` and
//! `U' <= sup S` (a simple rational between the roots of the compared rows
//! and the source end, else a 32-bit inner approximation of an irrational
//! end; `U' = 0` when the rows have no positive root and `S` reaches `0`)
//! are proved with the shared parent provers on the affinely transported
//! families `R_n(v + U')` (every root `<= 0`, i.e. `R_n`'s roots `<= U'`)
//! and `R_n(L' - v)` (roots `>= L'`); a parent prover proves that every row
//! is zero or real-rooted with coefficients of one sign, i.e. with every
//! root `<= 0`.  The lift is proved when the identity is certified and both
//! bounds are proved; then every `Q_n` is zero or real-rooted in the
//! conclusion's target `J` (and PF when the conclusion says so and `R_n`
//! has a positive leading coefficient).
//!
//! # The backward lift and the forward move
//!
//! [`find_basis_map_lifts`] tries every map of a [`BasisMapLibrary`] (by
//! default the admitted maps of the bundled OEIS catalog: proved, with a
//! companion form and triangular rows, not the identity, one map per
//! distinct set of first rows): it expands the input rows, rejects maps
//! whose expanded rows have a root outside every source (a near miss),
//! fits a derivative-free rule to the expanded rows and certifies the
//! candidate.  [`apply_basis_map_to_family`] builds the rule of `T_P(R)`
//! and certifies it the same way, so that an image can enter a proof
//! database with its certificate.  For a map of order one the rule is the
//! transported rule.  For a map of order `r >= 2` it is computed from the
//! transported system `Phi` (order `k`): a runnable scalar operator `M` with
//! `w(n) (M e_0) = sum_b H_b o Phi_b`, searched as a common left multiple of
//! order at most `r k` (`ore_algebra::find_component_annihilator`), which
//! is exactly the identity the certificate checks; only when nothing is
//! found within the budget is a rule fitted to the image rows (and then
//! certified in the same way).
//!
//! # Limitations
//!
//! Only coefficient rules without derivatives in `u` are transported (`T o
//! D_u` is not an operator image in general); index-dependent basis rules
//! and exceptional prefixes have no companion form (their theorems are
//! still checked); the backward lift needs triangular basis rows (so not
//! the constant-current maps, whose degrees repeat); interlacing of
//! consecutive image rows is not claimed.
//!
//! The group is split into `quadratic` (exact `a + b sqrt(d)` and
//! intervals), `map` (the map, its inverse and intertwining relation),
//! `theorems` (the shapes and their checks), `transport` (`Phi`), `certify`
//! (identity, root location, summary), `library`, `search` and `forward`.

mod map;
mod quadratic;
mod theorems;

pub use map::{
    BasisMap, BasisMapIntertwining, BasisMapIntertwiningTerm, BasisMapJson, MAX_BASIS_MAP_ROWS,
};
pub use quadratic::{AlgebraicEndpoint, AlgebraicInterval, QuadraticNumber};
pub use theorems::{
    check_basis_map_preserves_pf, BasisMapConclusion, BasisMapCounterexample, BasisMapShape,
    BasisMapTheoremOptions, BasisMapTheoremReport, BasisMapTheoremStatus, BASIS_MAP_THEOREM_SCHEMA,
};
mod certify;
mod transport;
pub use certify::{
    certify_basis_map_lift, BasisMapCoefficientRootProof, BasisMapDivisionCertificate,
    BasisMapDivisionDirection, BasisMapDivisionQuotient, BasisMapIdentityStatus,
    BasisMapLiftCandidate, BasisMapLiftCertification, BasisMapLiftOptions, BASIS_MAP_LIFT_SCHEMA,
    MAX_BASIS_MAP_LIFT_CHECKED_ROWS, MAX_BASIS_MAP_OPERATOR_TERMS,
};
mod forward;
mod library;
mod search;
pub use forward::{
    apply_basis_map_to_family, BasisMapImage, BasisMapImageConstruction, BasisMapImageOptions,
    BASIS_MAP_IMAGE_SCHEMA,
};
pub use library::{BasisMapLibrary, BasisMapLibraryEntry, BasisMapLibraryOptions};
pub use search::{
    find_basis_map_lifts, BasisMapLiftRejection, BasisMapLiftReport, BasisMapLiftSearchOptions,
    BASIS_MAP_LIFTS_SCHEMA,
};
