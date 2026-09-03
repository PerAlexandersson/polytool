use std::collections::BTreeMap;

use combinatoric_core::Graph;
use sym_poly_core::{Partition, UnivariatePolynomial};
use sym_poly_sym::llt::{
    circular_highest_reachable_vertices, circular_unicellular_llt_q_plus_one_e_expansion,
    circular_vertical_strip_llt_hrv_e_expansion, directed_graph_llt_symmetric,
    highest_reachable_vertices,
};
use sym_poly_sym::{Basis, SymmetricFunction};

#[test]
fn translates_mathematica_highest_reachable_vertices() {
    assert_eq!(
        highest_reachable_vertices(5, &[(0, 1), (1, 3), (0, 2)]),
        Some(vec![3, 3, 2, 3, 4])
    );
    assert!(highest_reachable_vertices(3, &[(0, 1), (1, 0)]).is_none());
    assert!(highest_reachable_vertices(3, &[(0, 3)]).is_none());
}

#[test]
fn circular_highest_reachable_vertices_use_the_universal_cover() {
    // In one-based notation these are 2 -> 3 and 3 -> 1. Their lift is
    // 2^(0) -> 3^(0) -> 1^(1), so every vertex has circular HRV 1.
    assert_eq!(
        circular_highest_reachable_vertices(3, &[(1, 2), (2, 0)]),
        Some(vec![0, 0, 0])
    );
    assert!(circular_highest_reachable_vertices(3, &[(0, 1), (1, 2), (2, 0)]).is_none());
}

#[test]
fn circular_hrv_expansion_matches_direct_llt_through_rank_four() {
    for n in 0..=4 {
        for area in all_circular_area_sequences(n) {
            let hrv = circular_vertical_strip_llt_hrv_e_expansion(&area, &[])
                .unwrap_or_else(|| panic!("HRV expansion rejected area {area:?}"));
            let direct = circular_unicellular_llt_q_plus_one_e_expansion(&area)
                .unwrap_or_else(|| panic!("direct LLT expansion rejected area {area:?}"));
            assert_eq!(hrv, direct, "area {area:?}");
        }
    }
}

#[test]
fn circular_hrv_expansion_matches_direct_strict_llt() {
    let area = [2, 2, 1];
    let strict_edges = [(1, 0)];
    let directed_edges = Graph::circular_unit_interval_directed_edges(&area).unwrap();
    let direct = substitute_q_plus_one_symmetric_function(
        &directed_graph_llt_symmetric(area.len(), &directed_edges, &strict_edges, &[]).unwrap(),
    )
    .to_elementary_basis();

    assert_eq!(
        circular_vertical_strip_llt_hrv_e_expansion(&area, &strict_edges),
        Some(direct)
    );
    assert!(circular_vertical_strip_llt_hrv_e_expansion(&area, &[(0, 2)]).is_none());
}

#[test]
fn circular_hrv_expansion_matches_four_generic_rank_eight_unicellular_examples() {
    let areas = [
        [1, 1, 1, 1, 1, 1, 1, 1],
        [2, 2, 1, 2, 2, 1, 2, 1],
        [3, 2, 2, 3, 2, 2, 1, 2],
        [2, 3, 3, 2, 1, 2, 2, 1],
    ];

    for area in areas {
        assert!(area.iter().all(|&entry| entry > 0));
        assert!(Graph::is_circular_unit_interval_area_sequence(&area));
        let hrv = circular_vertical_strip_llt_hrv_e_expansion(&area, &[]).unwrap();
        let direct = circular_unicellular_llt_q_plus_one_e_expansion(&area).unwrap();
        assert_eq!(hrv, direct, "area {area:?}");
    }
}

#[test]
fn circular_hrv_expansion_matches_four_generic_rank_eight_vertical_strip_examples() {
    let areas_and_corner_indices: [([u8; 8], &[usize]); 4] = [
        ([1, 2, 2, 1, 1, 2, 2, 1], &[0, 3]),
        ([2, 1, 2, 3, 2, 1, 2, 2], &[1, 4]),
        ([3, 3, 2, 1, 2, 2, 3, 2], &[0, 3, 4]),
        ([2, 2, 3, 2, 3, 2, 1, 1], &[1, 3, 4]),
    ];

    for (area, corner_indices) in areas_and_corner_indices {
        assert!(area.iter().all(|&entry| entry > 0));
        assert!(Graph::is_circular_unit_interval_area_sequence(&area));
        let corners = admissible_circular_corner_edges(&area);
        let strict_edges: Vec<_> = corner_indices.iter().map(|&index| corners[index]).collect();
        assert_hrv_matches_direct_strict_llt(&area, &strict_edges);
    }
}

fn assert_hrv_matches_direct_strict_llt(area: &[u8], strict_edges: &[(usize, usize)]) {
    let directed_edges = Graph::circular_unit_interval_directed_edges(area).unwrap();
    let direct = substitute_q_plus_one_symmetric_function(
        &directed_graph_llt_symmetric(area.len(), &directed_edges, strict_edges, &[]).unwrap(),
    )
    .to_elementary_basis();
    let hrv = circular_vertical_strip_llt_hrv_e_expansion(area, strict_edges).unwrap();
    assert!(
        !direct.terms().is_empty(),
        "strict-edge example is vacuous: area {area:?}, strict edges {strict_edges:?}"
    );
    assert_eq!(hrv, direct, "area {area:?}, strict edges {strict_edges:?}");
}

fn admissible_circular_corner_edges(area: &[u8]) -> Vec<(usize, usize)> {
    let n = area.len();
    (0..n)
        .filter(|&target| area[target] > 0 && area[(target + 1) % n] <= area[target])
        .map(|target| {
            let source = (target + n - usize::from(area[target])) % n;
            (source, target)
        })
        .collect()
}

fn substitute_q_plus_one_symmetric_function(
    function: &SymmetricFunction<UnivariatePolynomial<i64>>,
) -> SymmetricFunction<UnivariatePolynomial<i64>> {
    let terms: BTreeMap<Partition, UnivariatePolynomial<i64>> = function
        .terms()
        .iter()
        .map(|(partition, coefficient)| (partition.clone(), substitute_q_plus_one(coefficient)))
        .collect();
    SymmetricFunction::from_terms(Basis::Monomial, terms)
}

fn substitute_q_plus_one(polynomial: &UnivariatePolynomial<i64>) -> UnivariatePolynomial<i64> {
    let mut coefficients = vec![0; polynomial.coeffs().len()];
    for (degree, &coefficient) in polynomial.coeffs().iter().enumerate() {
        for target_degree in 0..=degree {
            coefficients[target_degree] += coefficient * binomial(degree, target_degree) as i64;
        }
    }
    UnivariatePolynomial::new(coefficients)
}

fn binomial(n: usize, k: usize) -> usize {
    let k = k.min(n - k);
    (0..k).fold(1usize, |result, i| result * (n - i) / (i + 1))
}

fn all_circular_area_sequences(n: usize) -> Vec<Vec<u8>> {
    if n == 0 {
        return vec![Vec::new()];
    }

    let mut result = Vec::new();
    let mut current = vec![0u8; n];
    circular_area_sequences_rec(0, &mut current, &mut result);
    result
}

fn circular_area_sequences_rec(index: usize, current: &mut [u8], result: &mut Vec<Vec<u8>>) {
    if index == current.len() {
        if Graph::is_circular_unit_interval_area_sequence(current) {
            result.push(current.to_vec());
        }
        return;
    }

    for value in 0..current.len() {
        current[index] = value as u8;
        circular_area_sequences_rec(index + 1, current, result);
    }
}
