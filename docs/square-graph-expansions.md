# Square partition graphs in the W and forgotten bases

We use Alexandersson--Shapiro, [arXiv:1104.0589, Section 2.3](https://arxiv.org/pdf/1104.0589),
and Definition 1's full permutation sum, denoted H_lambda here. An even part a
gives a doubled-edge star with a/2 leaves. Odd parts, sorted decreasingly, are
paired consecutively; a pair (a,b) gives two centers joined by a double edge,
with (a-1)/2 and (b-1)/2 doubled leaf edges. Thus this includes odd-part pairs
such as (3,3). It is not the square of a symmetrized partition graph.

We use f_mu = omega(m_mu) and W_mu^(n) = C_n^(-1)(f_mu), with the variable
count n retained in every calculation. W coefficients of H are forgotten
coefficients of C_n(H); direct forgotten coefficients omit this operator.

## W coordinates

The permutation averages in degree four, at n=8, are

    H_(4)/8!   = -270 W_(4) + 720 W_(2,2),
    H_(2,2)/8! =  576 W_(4) + 576 W_(2,2).

In degree six, at n=12, they are

    H_(6)/12!     = -100464 W_(6) + 137592 W_(4,2)
                    +209664 W_(3,3) -157248 W_(2,2,2),
    H_(4,2)/12!   = 34944 W_(6) -39312 W_(4,2)
                    +8736 W_(3,3) -104832 W_(2,2,2),
    H_(3,3)/12!   = 34944 W_(6) -96096 W_(3,3) -131040 W_(2,2,2),
    H_(2,2,2)/12! = -82368 (W_(6) + W_(4,2) + W_(2,2,2)).

In particular, being a sum of squares does not imply W positivity. The degree
four example also rules out uniform signs achieved by multiplying each H by
an overall sign. Tests cover every even-degree two-partition through degree 10
at n=d, 2d, and 2d+3 (78 cases): no singleton W indices, conversion roundtrips,
and exact recovery of the original graph using C_n inverse.

## Direct forgotten coordinates

Already H_(4) has mixed signs. For n >= 3, its expansion is

    H_(4) = (n-2)! [6 f_(2,2) -4 f_(3,1) -(n+1) f_(4)].

For verification without interpolation, expanding (x-y)^2 (x-z)^2 and
grouping exponent partitions gives

    H_(4) = (n-2)! [(n-1) m_(4) -4 m_(3,1) +6 m_(2,2)].

Use omega(m_(4))=-m_(4), omega(m_(3,1))=m_(3,1)+2m_(4), and
omega(m_(2,2))=m_(2,2)+m_(4) to obtain the forgotten coefficients above.
In contrast, at n=8,

    H_(2,2) = 384 [36 f_(4) +16 f_(3,1) +31 f_(2,2)
                  +11 f_(2,1,1) +6 f_(1,1,1,1)].

The disjoint doubled-edge family has a general signed-positivity property:
(-1)^k H_(2^k) is strictly forgotten-positive for n >= 2k. Indeed, in each
squared difference a term has support one with positive coefficient or
support two with negative coefficient. A product on k disjoint edges
therefore has monomial coefficient sign (-1)^(length(nu)-k). The
set-partition formula in [the earlier forgotten-basis note](partition-graph-forgotten.md)
shows that omega(m_nu) has monomial coefficients of sign
(-1)^(2k-length(nu)). These signs combine to (-1)^k throughout. Strictness
follows from the nonzero all-ones monomial and omega(m_(1^(2k)))=h_(2k).

Exact direct-forgotten computations tested every two-partition in even degrees
2 through 12 at n=d, 2d, and 2d+3: 141 cases, all roundtripping to monomials.
Of these, 123 had mixed signs; the other 18 were exactly the indices (2^k).
The claim that every other square-basis index has mixed signs is only an
observation on this finite range, not a proved general assertion.

Related Rust: `sym-poly/sym/src/graph_monomial.rs::square_partition_graph_symmetric`,
`content_operator.rs::c_n`, and the existing `.to_forgotten_basis()` conversion.
The maintained `sym-poly/sym/examples/square_graph_w.rs` prints both kinds of
expansions with exact rational normalization by n!. Literal permutation-sum
tests independently check small graphs including (3,3) and (5,3); further tests
check multiple odd pairs across interspersed even parts and the degree-four
formulas. Temporary direct-forgotten evidence is in
`/tmp/square-graph-forgotten-results.txt`, produced by
`/tmp/square_graph_forgotten_probe.rs` using the maintained generator.
