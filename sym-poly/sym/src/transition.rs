//! Transition matrices between the six classical symmetric function bases.
//!
//! The main entry point is [`convert`], which converts a `SymmetricFunction` from
//! one basis to another using precomputed transition matrices per degree.

use std::collections::BTreeMap;
use std::sync::Arc;

use sym_poly_core::matrix::{identity_matrix, invert_integer_matrix, mat_mul, transpose};
use sym_poly_core::{Partition, Ring, TransitionCache};

use crate::basis::Basis;
use crate::kostka;
use crate::symmetric_function::SymmetricFunction;
use crate::z_coefficient_i64;

// ---------------------------------------------------------------------------
// Global transition matrix cache (Sym-specific)
// ---------------------------------------------------------------------------

static SYM_CACHE: TransitionCache<Basis> = TransitionCache::new();

fn cached_transition_matrix(
    source: Basis,
    target: Basis,
    partitions: &[Partition],
    deg: u32,
) -> Arc<Vec<Vec<i64>>> {
    SYM_CACHE.get_or_compute(source, target, deg, |s, t, _| {
        transition_matrix(s, t, partitions)
    })
}

/// Convert a symmetric function to a different basis.
pub fn convert<C: Ring>(sf: &SymmetricFunction<C>, target: Basis) -> SymmetricFunction<C> {
    if sf.basis() == target {
        return sf.clone();
    }
    if sf.is_zero() {
        return SymmetricFunction::zero(target);
    }

    // Group terms by degree
    let mut by_degree: BTreeMap<u32, Vec<(Partition, C)>> = BTreeMap::new();
    for (p, c) in sf.terms() {
        by_degree
            .entry(p.size())
            .or_default()
            .push((p.clone(), c.clone()));
    }

    let mut result_terms: BTreeMap<Partition, C> = BTreeMap::new();

    let to_power_sum = target == Basis::PowerSum && sf.basis() != Basis::PowerSum;

    for (deg, terms) in by_degree {
        let partitions = Partition::all_of_size(deg);
        let part_index: BTreeMap<&Partition, usize> =
            partitions.iter().enumerate().map(|(i, p)| (p, i)).collect();
        let k = partitions.len();

        if to_power_sum {
            // source -> Schur (integer), then Schur -> PowerSum with z_μ division
            let schur_trans = if sf.basis() == Basis::Schur {
                Arc::new(identity_matrix(k))
            } else {
                cached_transition_matrix(sf.basis(), Basis::Schur, &partitions, deg)
            };

            let mut input = vec![C::zero(); k];
            for (p, c) in &terms {
                if let Some(&idx) = part_index.get(p) {
                    input[idx] = c.clone();
                }
            }

            let mut schur_vec = vec![C::zero(); k];
            for i in 0..k {
                let mut sum = C::zero();
                for j in 0..k {
                    if !input[j].is_zero() && schur_trans[j][i] != 0 {
                        sum = sum + input[j].clone() * C::from_i64(schur_trans[j][i]);
                    }
                }
                schur_vec[i] = sum;
            }

            for j in 0..k {
                let z_mu = z_coefficient_i64(&partitions[j]);
                let mut numer = C::zero();
                for i in 0..k {
                    if schur_vec[i].is_zero() {
                        continue;
                    }
                    let chi = kostka::sn_character(&partitions[i], &partitions[j]);
                    if chi != 0 {
                        numer = numer + schur_vec[i].clone() * C::from_i64(chi);
                    }
                }
                if !numer.is_zero() {
                    let coeff = numer.exact_div_i64(z_mu);
                    if !coeff.is_zero() {
                        result_terms.insert(partitions[j].clone(), coeff);
                    }
                }
            }

            continue;
        }

        // Standard integer transition matrix approach
        let trans = cached_transition_matrix(sf.basis(), target, &partitions, deg);

        let mut input = vec![C::zero(); k];
        for (p, c) in &terms {
            if let Some(&idx) = part_index.get(p) {
                input[idx] = c.clone();
            }
        }

        for i in 0..k {
            let mut sum = C::zero();
            for j in 0..k {
                if !input[j].is_zero() && trans[j][i] != 0 {
                    sum = sum + input[j].clone() * C::from_i64(trans[j][i]);
                }
            }
            if !sum.is_zero() {
                result_terms.insert(partitions[i].clone(), sum);
            }
        }
    }

    SymmetricFunction::from_terms(target, result_terms)
}

/// Compute the transition matrix M where source[j] = Σ_i M[j][i] * target[i].
fn transition_matrix(source: Basis, target: Basis, partitions: &[Partition]) -> Vec<Vec<i64>> {
    use Basis::*;

    if source == target {
        return identity_matrix(partitions.len());
    }

    let k = partitions.len();
    if k == 0 {
        return vec![];
    }

    match (source, target) {
        // -------------------------------------------------------------------
        // Direct Kostka-based conversions
        // -------------------------------------------------------------------
        (Schur, Monomial) => {
            let mut mat = vec![vec![0i64; k]; k];
            for i in 0..k {
                for j in 0..k {
                    mat[i][j] = kostka::kostka_coefficient(&partitions[i], &partitions[j]);
                }
            }
            mat
        }

        (Monomial, Schur) => {
            let sm = transition_matrix(Schur, Monomial, partitions);
            invert_integer_matrix(&sm)
        }

        (CompleteH, Schur) => {
            let mut mat = vec![vec![0i64; k]; k];
            for j in 0..k {
                for i in 0..k {
                    mat[j][i] = kostka::kostka_coefficient(&partitions[i], &partitions[j]);
                }
            }
            mat
        }

        (Schur, CompleteH) => {
            let hs = transition_matrix(CompleteH, Schur, partitions);
            invert_integer_matrix(&hs)
        }

        (CompleteH, Monomial) => {
            let hs = transition_matrix(CompleteH, Schur, partitions);
            let sm = transition_matrix(Schur, Monomial, partitions);
            mat_mul(&hs, &sm)
        }

        (Monomial, CompleteH) => {
            let mh = transition_matrix(CompleteH, Monomial, partitions);
            invert_integer_matrix(&mh)
        }

        (Elementary, Schur) => {
            let mut mat = vec![vec![0i64; k]; k];
            for j in 0..k {
                for i in 0..k {
                    let conj_i = partitions[i].conjugate_partition();
                    mat[j][i] = kostka::kostka_coefficient(&conj_i, &partitions[j]);
                }
            }
            mat
        }

        (Schur, Elementary) => {
            let es = transition_matrix(Elementary, Schur, partitions);
            invert_integer_matrix(&es)
        }

        (Elementary, Monomial) => {
            let es = transition_matrix(Elementary, Schur, partitions);
            let sm = transition_matrix(Schur, Monomial, partitions);
            mat_mul(&es, &sm)
        }

        (Monomial, Elementary) => {
            let me = transition_matrix(Elementary, Monomial, partitions);
            invert_integer_matrix(&me)
        }

        (Elementary, CompleteH) => {
            let es = transition_matrix(Elementary, Schur, partitions);
            let sh = transition_matrix(Schur, CompleteH, partitions);
            mat_mul(&es, &sh)
        }

        (CompleteH, Elementary) => {
            let ce = transition_matrix(Elementary, CompleteH, partitions);
            invert_integer_matrix(&ce)
        }

        // -------------------------------------------------------------------
        // Power sum conversions via characters
        // -------------------------------------------------------------------
        (PowerSum, Schur) => {
            let mut mat = vec![vec![0i64; k]; k];
            for j in 0..k {
                for i in 0..k {
                    mat[j][i] = kostka::sn_character(&partitions[i], &partitions[j]);
                }
            }
            mat
        }

        (PowerSum, Monomial) => {
            let ps = transition_matrix(PowerSum, Schur, partitions);
            let sm = transition_matrix(Schur, Monomial, partitions);
            mat_mul(&ps, &sm)
        }

        (PowerSum, CompleteH) => {
            let ps = transition_matrix(PowerSum, Schur, partitions);
            let sh = transition_matrix(Schur, CompleteH, partitions);
            mat_mul(&ps, &sh)
        }

        (PowerSum, Elementary) => {
            let ps = transition_matrix(PowerSum, Schur, partitions);
            let se = transition_matrix(Schur, Elementary, partitions);
            mat_mul(&ps, &se)
        }

        (Schur, PowerSum)
        | (Monomial, PowerSum)
        | (CompleteH, PowerSum)
        | (Elementary, PowerSum) => {
            panic!(
                "transition_matrix({:?}, PowerSum) has rational entries; \
                 use convert() which handles z_μ division",
                source
            );
        }

        // -------------------------------------------------------------------
        // Forgotten basis
        // -------------------------------------------------------------------
        (Forgotten, Monomial) => {
            let eh = transition_matrix(Elementary, CompleteH, partitions);
            transpose(&eh)
        }

        (Monomial, Forgotten) => {
            let fm = transition_matrix(Forgotten, Monomial, partitions);
            invert_integer_matrix(&fm)
        }

        (Forgotten, Schur) => {
            let fm = transition_matrix(Forgotten, Monomial, partitions);
            let ms = transition_matrix(Monomial, Schur, partitions);
            mat_mul(&fm, &ms)
        }

        (Schur, Forgotten) => {
            let sf = transition_matrix(Forgotten, Schur, partitions);
            invert_integer_matrix(&sf)
        }

        (Forgotten, CompleteH) => {
            let fm = transition_matrix(Forgotten, Monomial, partitions);
            let mh = transition_matrix(Monomial, CompleteH, partitions);
            mat_mul(&fm, &mh)
        }

        (CompleteH, Forgotten) => {
            let cf = transition_matrix(Forgotten, CompleteH, partitions);
            invert_integer_matrix(&cf)
        }

        (Forgotten, Elementary) => {
            let fm = transition_matrix(Forgotten, Monomial, partitions);
            let me = transition_matrix(Monomial, Elementary, partitions);
            mat_mul(&fm, &me)
        }

        (Elementary, Forgotten) => {
            let ef = transition_matrix(Forgotten, Elementary, partitions);
            invert_integer_matrix(&ef)
        }

        (Forgotten, PowerSum) => {
            panic!(
                "transition_matrix(Forgotten, PowerSum) has rational entries; \
                 use convert() which handles z_μ division"
            );
        }

        (PowerSum, Forgotten) => {
            let ps = transition_matrix(PowerSum, Schur, partitions);
            let sf = transition_matrix(Schur, Forgotten, partitions);
            mat_mul(&ps, &sf)
        }

        _ => unreachable!("all basis pairs should be covered"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_schur_to_monomial_degree3() {
        let s: SymmetricFunction<i64> =
            SymmetricFunction::schur_symmetric(Partition::new(vec![2, 1]));
        let m = s.to_monomial_basis();
        assert_eq!(m.basis(), Basis::Monomial);
        assert_eq!(m.coefficient(&Partition::new(vec![2, 1])), 1);
        assert_eq!(m.coefficient(&Partition::new(vec![1, 1, 1])), 2);
    }

    #[test]
    fn test_transition_cache_reuses_matrix_allocation() {
        let partitions = Partition::all_of_size(0);
        let first = cached_transition_matrix(Basis::Schur, Basis::Monomial, &partitions, 0);
        let second = cached_transition_matrix(Basis::Schur, Basis::Monomial, &partitions, 0);
        assert_eq!(&*first, &*second);
        assert!(Arc::ptr_eq(&first, &second));
    }

    #[test]
    fn test_roundtrip_all_bases_degree2() {
        let original: SymmetricFunction<i64> =
            SymmetricFunction::schur_symmetric(Partition::new(vec![2]));

        for &basis in &[
            Basis::Monomial,
            Basis::Elementary,
            Basis::CompleteH,
            Basis::Forgotten,
        ] {
            let converted = original.to_basis(basis);
            let back = converted.to_schur_basis();
            assert_eq!(
                back.coefficient(&Partition::new(vec![2])),
                1,
                "roundtrip through {:?}",
                basis
            );
            assert_eq!(
                back.coefficient(&Partition::new(vec![1, 1])),
                0,
                "roundtrip through {:?}",
                basis
            );
        }
    }

    #[test]
    fn test_power_sum_conversions_rational() {
        use num_rational::Ratio;
        type Q = Ratio<i64>;

        let p21: SymmetricFunction<Q> =
            SymmetricFunction::power_sum_symmetric(Partition::new(vec![2, 1]));
        let in_schur = p21.to_schur_basis();
        assert_eq!(
            in_schur.coefficient(&Partition::new(vec![3])),
            Q::from_integer(1)
        );
        assert_eq!(
            in_schur.coefficient(&Partition::new(vec![2, 1])),
            Q::from_integer(0)
        );
        assert_eq!(
            in_schur.coefficient(&Partition::new(vec![1, 1, 1])),
            Q::from_integer(-1)
        );

        let s2: SymmetricFunction<Q> = SymmetricFunction::schur_symmetric(Partition::new(vec![2]));
        let in_p = s2.to_power_sum_basis();
        assert_eq!(in_p.coefficient(&Partition::new(vec![2])), Ratio::new(1, 2));
        assert_eq!(
            in_p.coefficient(&Partition::new(vec![1, 1])),
            Ratio::new(1, 2)
        );

        let back = in_p.to_schur_basis();
        assert_eq!(
            back.coefficient(&Partition::new(vec![2])),
            Q::from_integer(1)
        );
        assert_eq!(
            back.coefficient(&Partition::new(vec![1, 1])),
            Q::from_integer(0)
        );
    }
}
