# Follow-up: exact transition data for large Schur expansions

The partition-graph implementation constructs arbitrarily large BigInt
coefficients (tested at n=30), but the shared Sym conversion still uses i64
Kostka/character tables and transition matrices. This limitation is documented
in the crate and new API; it is independent of the graph coefficient type.

Before promoting high-degree expansion scans, add an arbitrary-precision or
checked-overflow transition path in the owning Sym/Kostka layer, with comparison
against existing low-degree tables and tests that exceed i64. Avoid a second
Schur-conversion algorithm inside graph_monomial. This is a proposal, not a
change to the current small-degree API or a claim of an observed overflow.
