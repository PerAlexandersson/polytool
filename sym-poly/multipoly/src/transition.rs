//! Transition matrices and basis conversion for nonsymmetric polynomial bases.
//!
//! Follows the same pattern as `sym-poly-sym/src/transition.rs`: for each
//! (source, target, degree, num_vars) tuple, build and cache a transition matrix,
//! then apply it via matrix-vector multiplication.

use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};

use combinatoric_core::WeakComposition;
use sym_poly_core::matrix::{invert_integer_matrix, mat_mul};
use sym_poly_core::Ring;

use crate::atom_polynomial::atom_polynomial;
use crate::basis::MultiPolyBasis;
use crate::key_polynomial::key_polynomial;
use crate::multipoly::MultiPoly;
use crate::multipoly_function::MultiPolyFunction;
use crate::slide_polynomial::{fundamental_slide_polynomial, monomial_slide_polynomial};

// =========================================================================
// Cache
// =========================================================================

/// Cache key: (source_basis, target_basis, degree, num_vars).
type CacheKey = (MultiPolyBasis, MultiPolyBasis, u32, usize);

static CACHE: Mutex<BTreeMap<CacheKey, Arc<Vec<Vec<i64>>>>> = Mutex::new(BTreeMap::new());

fn cached_matrix(
    source: MultiPolyBasis,
    target: MultiPolyBasis,
    degree: u32,
    num_vars: usize,
) -> Arc<Vec<Vec<i64>>> {
    let key = (source, target, degree, num_vars);

    // Fast path
    {
        let guard = CACHE.lock().unwrap();
        if let Some(mat) = guard.get(&key) {
            return Arc::clone(mat);
        }
    }

    // Compute outside the lock
    let mat = Arc::new(compute_transition_matrix(source, target, degree, num_vars));

    {
        let mut guard = CACHE.lock().unwrap();
        return Arc::clone(guard.entry(key).or_insert(mat));
    }
}

// =========================================================================
// Matrix computation
// =========================================================================

/// Build the transition matrix from `source` to `target` for the given
/// (degree, num_vars). Rows are indexed by source basis elements, columns
/// by target basis elements, both in the canonical weak-composition order.
fn compute_transition_matrix(
    source: MultiPolyBasis,
    target: MultiPolyBasis,
    degree: u32,
    num_vars: usize,
) -> Vec<Vec<i64>> {
    if source == target {
        let n = WeakComposition::all_weak_compositions(degree, num_vars).len();
        return identity(n);
    }

    // Strategy: compute X -> Monomial directly, then invert or compose.
    match (source, target) {
        // Direct: basis -> Monomial (expand each basis element)
        (MultiPolyBasis::Key, MultiPolyBasis::Monomial) => {
            basis_to_monomial_matrix(degree, num_vars, |alpha| key_polynomial::<i64>(alpha))
        }
        (MultiPolyBasis::Atom, MultiPolyBasis::Monomial) => {
            basis_to_monomial_matrix(degree, num_vars, |alpha| atom_polynomial::<i64>(alpha))
        }
        (MultiPolyBasis::FundSlide, MultiPolyBasis::Monomial) => {
            basis_to_monomial_matrix(degree, num_vars, |alpha| {
                fundamental_slide_polynomial::<i64>(alpha)
            })
        }
        (MultiPolyBasis::MonSlide, MultiPolyBasis::Monomial) => {
            basis_to_monomial_matrix(degree, num_vars, |alpha| {
                monomial_slide_polynomial::<i64>(alpha)
            })
        }

        // Inverse: Monomial -> basis (invert the basis -> Monomial matrix)
        (MultiPolyBasis::Monomial, _) => {
            let fwd = cached_matrix(target, MultiPolyBasis::Monomial, degree, num_vars);
            invert_integer_matrix(&fwd)
        }

        // General: source -> Monomial -> target (compose)
        _ => {
            let s_to_m = cached_matrix(source, MultiPolyBasis::Monomial, degree, num_vars);
            let m_to_t = cached_matrix(MultiPolyBasis::Monomial, target, degree, num_vars);
            mat_mul(&s_to_m, &m_to_t)
        }
    }
}

/// Build the matrix taking basis elements to monomials.
/// Row i = source basis element α_i, Column j = monomial x^{α_j}.
/// Entry [i][j] = coefficient of x^{α_j} in basis(α_i).
fn basis_to_monomial_matrix(
    degree: u32,
    num_vars: usize,
    compute: impl Fn(&[u32]) -> MultiPoly<i64>,
) -> Vec<Vec<i64>> {
    let wcs = WeakComposition::all_weak_compositions(degree, num_vars);
    let n = wcs.len();
    if n == 0 {
        return vec![];
    }

    // Build index: weak composition -> column index
    let wc_index: BTreeMap<&[u32], usize> = wcs
        .iter()
        .enumerate()
        .map(|(i, wc)| (wc.parts(), i))
        .collect();

    let mut mat = vec![vec![0i64; n]; n];

    for (i, wc) in wcs.iter().enumerate() {
        let poly = compute(wc.parts());
        for (exp, coeff) in poly.terms() {
            if let Some(&j) = wc_index.get(exp.as_slice()) {
                mat[i][j] = *coeff;
            }
        }
    }

    mat
}

fn identity(n: usize) -> Vec<Vec<i64>> {
    let mut m = vec![vec![0i64; n]; n];
    for i in 0..n {
        m[i][i] = 1;
    }
    m
}

// =========================================================================
// Conversion
// =========================================================================

/// Convert a `MultiPolyFunction` to a different basis.
pub fn convert<C: Ring>(f: &MultiPolyFunction<C>, target: MultiPolyBasis) -> MultiPolyFunction<C> {
    if f.basis() == target {
        return f.clone();
    }
    if f.is_zero() {
        return MultiPolyFunction::zero(target, f.num_vars());
    }

    let source = f.basis();
    let num_vars = f.num_vars();

    // Group terms by degree
    let mut by_degree: BTreeMap<u32, Vec<(WeakComposition, C)>> = BTreeMap::new();
    for (wc, c) in f.terms() {
        by_degree
            .entry(wc.degree())
            .or_default()
            .push((wc.clone(), c.clone()));
    }

    let mut result_terms: BTreeMap<WeakComposition, C> = BTreeMap::new();

    for (deg, terms) in by_degree {
        let wcs = WeakComposition::all_weak_compositions(deg, num_vars);
        let n = wcs.len();
        if n == 0 {
            continue;
        }

        // Build index for this degree
        let wc_index: BTreeMap<&WeakComposition, usize> =
            wcs.iter().enumerate().map(|(i, wc)| (wc, i)).collect();

        // Input vector
        let mut input = vec![C::zero(); n];
        for (wc, c) in &terms {
            if let Some(&idx) = wc_index.get(wc) {
                input[idx] = c.clone();
            }
        }

        // Get transition matrix
        let mat = cached_matrix(source, target, deg, num_vars);

        // Apply: output[j] = Σ_i input[i] * mat[i][j]
        for j in 0..n {
            let mut sum = C::zero();
            for i in 0..n {
                if !input[i].is_zero() && mat[i][j] != 0 {
                    sum = sum + input[i].clone() * C::from_i64(mat[i][j]);
                }
            }
            if !sum.is_zero() {
                result_terms.insert(wcs[j].clone(), sum);
            }
        }
    }

    MultiPolyFunction::from_terms(target, num_vars, result_terms)
}

/// Clear the global transition matrix cache.
pub fn clear_cache() {
    let mut guard = CACHE.lock().unwrap();
    guard.clear();
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::MultiPolyBasis::*;

    /// Helper: create a single basis element and verify round-trip.
    fn round_trip_single(basis: MultiPolyBasis, alpha: &[u32], target: MultiPolyBasis) {
        let f = MultiPolyFunction::<i64>::basis_element(basis, alpha);
        let converted = f.to_basis(target);
        let back = converted.to_basis(basis);
        assert_eq!(
            f, back,
            "round-trip failed: {:?} {:?} via {:?}",
            basis, alpha, target
        );
    }

    #[test]
    fn test_identity_conversion() {
        let f = MultiPolyFunction::<i64>::key(&[1, 0, 2]);
        let same = f.to_basis(Key);
        assert_eq!(f, same);
    }

    #[test]
    fn test_transition_cache_reuses_matrix_allocation() {
        // Use a degree/variable count not used by the conversion tests below,
        // so this checks the first computation and the warm-cache path.
        let first = cached_matrix(Key, Monomial, 3, 7);
        let second = cached_matrix(Key, Monomial, 3, 7);
        assert_eq!(&*first, &*second);
        assert!(Arc::ptr_eq(&first, &second));
    }

    #[test]
    fn test_key_to_monomial_matches_direct() {
        // κ_α via MultiPolyFunction should match key_polynomial() directly
        for alpha in &[
            vec![2, 1],
            vec![1, 2],
            vec![0, 3],
            vec![3, 0],
            vec![1, 0, 2],
            vec![0, 1, 2],
            vec![2, 0, 1],
        ] {
            let f = MultiPolyFunction::<i64>::key(alpha);
            let expanded = f.to_multipoly();
            let direct: MultiPoly<i64> = key_polynomial(alpha);
            assert_eq!(expanded, direct, "key mismatch for alpha = {:?}", alpha);
        }
    }

    #[test]
    fn test_atom_to_monomial_matches_direct() {
        for alpha in &[vec![2, 1], vec![1, 2], vec![0, 2], vec![1, 0, 2]] {
            let f = MultiPolyFunction::<i64>::atom(alpha);
            let expanded = f.to_multipoly();
            let direct: MultiPoly<i64> = atom_polynomial(alpha);
            assert_eq!(expanded, direct, "atom mismatch for alpha = {:?}", alpha);
        }
    }

    #[test]
    fn test_fund_slide_to_monomial_matches_direct() {
        for alpha in &[vec![1, 0, 2], vec![0, 2, 0, 3]] {
            let f = MultiPolyFunction::<i64>::fund_slide(alpha);
            let expanded = f.to_multipoly();
            let direct: MultiPoly<i64> = fundamental_slide_polynomial(alpha);
            assert_eq!(
                expanded, direct,
                "fund_slide mismatch for alpha = {:?}",
                alpha
            );
        }
    }

    #[test]
    fn test_mon_slide_to_monomial_matches_direct() {
        for alpha in &[vec![0, 2, 0, 3], vec![1, 0, 2]] {
            let f = MultiPolyFunction::<i64>::mon_slide(alpha);
            let expanded = f.to_multipoly();
            let direct: MultiPoly<i64> = monomial_slide_polynomial(alpha);
            assert_eq!(
                expanded, direct,
                "mon_slide mismatch for alpha = {:?}",
                alpha
            );
        }
    }

    #[test]
    fn test_round_trip_key_monomial() {
        for alpha in &[vec![1, 2], vec![0, 3], vec![1, 0, 2], vec![2, 1, 0]] {
            round_trip_single(Key, alpha, Monomial);
        }
    }

    #[test]
    fn test_round_trip_atom_monomial() {
        for alpha in &[vec![1, 2], vec![0, 2], vec![1, 0, 2]] {
            round_trip_single(Atom, alpha, Monomial);
        }
    }

    #[test]
    fn test_round_trip_key_atom() {
        for alpha in &[vec![1, 2], vec![0, 3], vec![1, 0, 2], vec![0, 1, 2]] {
            round_trip_single(Key, alpha, Atom);
        }
    }

    #[test]
    fn test_round_trip_fund_slide_monomial() {
        for alpha in &[vec![1, 0, 2], vec![0, 3], vec![2, 1]] {
            round_trip_single(FundSlide, alpha, Monomial);
        }
    }

    #[test]
    fn test_round_trip_mon_slide_monomial() {
        for alpha in &[vec![1, 0, 2], vec![0, 3], vec![2, 1]] {
            round_trip_single(MonSlide, alpha, Monomial);
        }
    }

    #[test]
    fn test_round_trip_all_bases_degree3_2vars() {
        let bases = [Monomial, Key, Atom, MonSlide, FundSlide];
        // Test round-trip through every pair of bases
        for &source in &bases {
            let f = MultiPolyFunction::<i64>::basis_element(source, &[1, 2]);
            for &target in &bases {
                let converted = f.to_basis(target);
                let back = converted.to_basis(source);
                assert_eq!(
                    f, back,
                    "round-trip {:?} -> {:?} -> {:?} failed for (1,2)",
                    source, target, source
                );
            }
        }
    }

    #[test]
    fn test_key_to_atom_positive() {
        // Key polynomials expand positively in atoms
        for alpha in &[
            vec![1, 2],
            vec![0, 3],
            vec![1, 0, 2],
            vec![0, 1, 2],
            vec![2, 0, 1],
        ] {
            let f = MultiPolyFunction::<i64>::key(alpha);
            let in_atoms = f.to_atom_basis();
            assert!(
                in_atoms.positive_coefficients(),
                "key {:?} in atom basis has non-positive coefficients: {}",
                alpha,
                in_atoms
            );
        }
    }

    #[test]
    fn test_key_to_fund_slide_positive() {
        // Key polynomials expand positively in fundamental slides (Assaf-Searles 2018)
        for alpha in &[vec![1, 2], vec![0, 3], vec![1, 0, 2], vec![0, 1, 2]] {
            let f = MultiPolyFunction::<i64>::key(alpha);
            let in_slides = f.to_fund_slide_basis();
            assert!(
                in_slides.positive_coefficients(),
                "key {:?} in fund_slide basis has non-positive coefficients: {}",
                alpha,
                in_slides
            );
        }
    }

    #[test]
    fn test_fund_slide_to_mon_slide_positive() {
        // Fundamental slides expand positively in monomial slides
        for alpha in &[vec![1, 2], vec![0, 3], vec![1, 0, 2]] {
            let f = MultiPolyFunction::<i64>::fund_slide(alpha);
            let in_mon = f.to_mon_slide_basis();
            assert!(
                in_mon.positive_coefficients(),
                "fund_slide {:?} in mon_slide basis has non-positive coefficients: {}",
                alpha,
                in_mon
            );
        }
    }

    #[test]
    fn test_schubert_to_fund_slide_positive() {
        // Schubert polynomials expand positively in fundamental slides
        use crate::schubert_polynomial::schubert_to_fund_slide;
        for perm in &[vec![2, 1, 3], vec![1, 3, 2], vec![2, 3, 1], vec![3, 1, 2]] {
            let f: MultiPolyFunction<i64> = schubert_to_fund_slide(perm);
            assert!(
                f.positive_coefficients(),
                "schubert {:?} in fund_slide has non-positive coefficients: {}",
                perm,
                f
            );
        }
    }

    #[test]
    fn test_schubert_to_key_positive() {
        use crate::schubert_polynomial::schubert_to_key;
        for perm in &[vec![2, 1, 3], vec![1, 3, 2], vec![2, 3, 1], vec![3, 1, 2]] {
            let f: MultiPolyFunction<i64> = schubert_to_key(perm);
            assert!(
                f.positive_coefficients(),
                "schubert {:?} in key basis has non-positive coefficients: {}",
                perm,
                f
            );
        }
    }

    #[test]
    fn test_schubert_monomial_matches_direct() {
        use crate::schubert_polynomial::{schubert_polynomial, schubert_to_monomial};
        for perm in &[vec![2, 1, 3], vec![1, 3, 2], vec![3, 2, 1]] {
            let f: MultiPolyFunction<i64> = schubert_to_monomial(perm);
            let direct: MultiPoly<i64> = schubert_polynomial(perm);
            assert_eq!(
                f.to_multipoly(),
                direct,
                "schubert monomial mismatch for perm {:?}",
                perm
            );
        }
    }

    #[test]
    fn test_arithmetic_in_basis() {
        let a = MultiPolyFunction::<i64>::key(&[1, 2]);
        let b = MultiPolyFunction::<i64>::key(&[0, 3]);
        let sum = a.clone() + b.clone();
        assert_eq!(sum.coefficient(&WeakComposition::from_slice(&[1, 2])), 1);
        assert_eq!(sum.coefficient(&WeakComposition::from_slice(&[0, 3])), 1);
        let diff = sum - b;
        assert_eq!(diff, a);
    }

    #[test]
    fn test_from_multipoly_round_trip() {
        let p: MultiPoly<i64> = key_polynomial(&[1, 0, 2]);
        let f = MultiPolyFunction::<i64>::from_multipoly(&p);
        let back = f.to_multipoly();
        assert_eq!(p, back);
    }
}
