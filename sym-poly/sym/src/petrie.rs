//! Petrie symmetric functions.
//!
//! The Petrie symmetric function `G(k,n)` is the sum of monomial symmetric
//! functions indexed by partitions of `n` whose largest part is strictly
//! smaller than `k`.

use std::collections::BTreeMap;

use sym_poly_core::{Partition, Ring};

use crate::{Basis, SymmetricFunction};

/// The Petrie symmetric function `G(k,n)` in the monomial basis.
///
/// This uses the convention `G(2,n) = e_n` and `G(k,n) = h_n` when `k > n`.
pub fn petrie_symmetric<C: Ring>(k: u32, n: u32) -> SymmetricFunction<C> {
    assert!(k >= 1, "k must be positive");

    let terms = Partition::all_of_size(n)
        .into_iter()
        .filter(|lambda| largest_part_is_less_than(lambda, k))
        .map(|lambda| (lambda, C::one()))
        .collect::<BTreeMap<_, _>>();

    SymmetricFunction::from_terms(Basis::Monomial, terms)
}

fn largest_part_is_less_than(lambda: &Partition, k: u32) -> bool {
    match lambda.parts().first() {
        Some(&part) => part < k,
        None => true,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn partition(parts: &[u32]) -> Partition {
        Partition::new(parts.to_vec())
    }

    #[test]
    fn test_petrie_g_2_n_is_elementary() {
        let g = petrie_symmetric::<i64>(2, 4);

        assert_eq!(g.basis(), Basis::Monomial);
        assert_eq!(g.coefficient(&partition(&[1, 1, 1, 1])), 1);
        assert_eq!(g.terms().len(), 1);

        let e = g.to_elementary_basis();
        assert_eq!(e.coefficient(&partition(&[4])), 1);
        assert_eq!(e.terms().len(), 1);
    }

    #[test]
    fn test_petrie_g_k_n_is_complete_when_k_exceeds_n() {
        let g = petrie_symmetric::<i64>(5, 4).to_complete_h_basis();

        assert_eq!(g.coefficient(&partition(&[4])), 1);
        assert_eq!(g.terms().len(), 1);
    }

    #[test]
    fn test_petrie_g_4_8_schur_expansion() {
        let s = petrie_symmetric::<i64>(4, 8).to_schur_basis();

        assert_eq!(s.coefficient(&partition(&[3, 3, 2])), 1);
        assert_eq!(s.coefficient(&partition(&[2, 2, 2, 2])), 1);
        assert_eq!(s.coefficient(&partition(&[3, 2, 2, 1])), -1);
        assert_eq!(s.coefficient(&partition(&[3, 1, 1, 1, 1, 1])), 1);
        assert_eq!(s.coefficient(&partition(&[2, 1, 1, 1, 1, 1, 1])), -1);
        assert_eq!(s.coefficient(&partition(&[1, 1, 1, 1, 1, 1, 1, 1])), 1);
        assert_eq!(s.terms().len(), 6);
    }
}
