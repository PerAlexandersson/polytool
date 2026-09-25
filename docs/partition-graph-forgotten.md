# Partition graphs in the forgotten basis

We use the full permutation sum `b_lambda^(n)` of arXiv:1104.0589,
Definition 1 and Section 2.2, and the convention `f_mu = omega(m_mu)`.
Let `d = |lambda|` and `r = length(lambda)`, with every part at least two.

For every `n >= d+r` (so the disjoint stars fit), every forgotten coefficient
of `(-1)^(d-r) b_lambda^(n)` is strictly positive. In particular this holds
throughout the paper's stable range `n >= 2d`.

Proof. For an outward star with a leaves, choosing k leaf terms gives
exponent type `(a-k,1^k)` with signed total `(-1)^k binomial(a,k)` for
`0 <= k <= a-2`. The remaining two choices combine at type `(1^a)` with
signed total `(-1)^(a-1)(a-1)`. Thus a term of support size s has sign
`(-1)^(s-1)`. Taking disjoint stars and then symmetrizing gives

    b_lambda^(n) = sum_nu (-1)^(length(nu)-r) A_lambda,nu m_nu,

where all A are nonnegative and
`A_lambda,(1^d) = d! (n-d)! product_i (lambda_i-1) > 0`.

For `nu = (nu_1,...,nu_l)`, inclusion-exclusion over set partitions pi of
`{1,...,l}` gives

    (-1)^(d-l) omega(m_nu)
      = 1 / product_j multiplicity_nu(j)!
        * sum_pi product_{B in pi} (|B|-1)! * product_{B in pi} p_{sum_{i in B} nu_i}.

The right-hand side is power-sum positive and hence monomial positive. All
signs therefore agree after applying omega to the expansion of b. Strictness
for every index follows from `omega(m_(1^d)) = h_d = sum_{mu |- d} m_mu`.
Only this term contributes to the coefficient of `m_(1^d)` after omega, so
the coefficient of `f_(1^d)` in b is exactly
`(-1)^(d-r) d! (n-d)! product_i (lambda_i-1)`.

Examples, with `B_lambda = b_lambda^(2d) / d!`:

    B_(2)   = -5 f_(2) - 2 f_(1,1)
    B_(3)   = 56 f_(3) + 24 f_(2,1) + 12 f_(1,1,1)
    B_(2,2) = 144 f_(4) + 64 f_(3,1) + 124 f_(2,2)
              + 44 f_(2,1,1) + 24 f_(1,1,1,1).

Related Rust: `sym-poly/sym/src/graph_monomial.rs::partition_graph_symmetric`
(source commit `187fcc9`) followed by `.to_forgotten_basis()`.
Exact computations tested all 41 two-partitions of degrees 2--10 at
`n=d+r`, `n=2d`, and `n=2d+3` (123 cases). Every coefficient has the stated
strict sign, full support and the predicted `f_(1^d)` coefficient. Each
expansion roundtripped to monomials and agreed with Schur conjugation followed
by monomial conversion. This finite check supports the proof above, rather
than replacing it. The bounded Rust driver and detailed output for this run
are `/tmp/partition_graph_forgotten_probe.rs` and
`/tmp/partition_graph_forgotten_results.txt`; no maintained library code changed.
