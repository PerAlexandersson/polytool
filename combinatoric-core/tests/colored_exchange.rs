//! Integration tests for `combinatoric_core::colored_exchange`.
//!
//! The brute-force checker below enumerates position subsets, as the
//! original bounded project probe did, and is deliberately independent of the
//! library's cycle traversal.

use std::collections::BTreeMap;

use combinatoric_core::colored_exchange::{
    letters_to_word, simple_exchange_cycles, CertificateError, CycleWitness, Factor, FamilyError,
    InvalidQuery, MarginalError, WitnessError,
};
use combinatoric_core::{BinaryPortCertificate, FixedContentFamily, ProductCertificate};
use num_bigint::BigUint;
use num_traits::{One, Pow, Zero};

fn family(num_colors: usize, letters: &[&str]) -> FixedContentFamily {
    FixedContentFamily::from_letter_words(num_colors, letters).unwrap()
}

fn word(letters: &str) -> Vec<usize> {
    letters_to_word(letters).unwrap()
}

/// Independent subset-enumeration checker of symmetric cycle exchange.
fn brute_force_sce(words: &[Vec<usize>], num_colors: usize) -> bool {
    let n = words[0].len();
    let is_member = |w: &[usize]| words.iter().any(|x| x.as_slice() == w);
    for (i, u) in words.iter().enumerate() {
        for v in &words[i + 1..] {
            for p in 0..n {
                if u[p] == v[p] {
                    continue;
                }
                let mut found = false;
                'masks: for mask in 1usize..1 << n {
                    if mask & (1 << p) == 0 {
                        continue;
                    }
                    let positions: Vec<usize> = (0..n).filter(|q| mask & (1 << q) != 0).collect();
                    if positions.len() < 2 || positions.len() > num_colors {
                        continue;
                    }
                    let mut next = vec![None; num_colors];
                    for &q in &positions {
                        if u[q] == v[q] || next[u[q]].is_some() {
                            continue 'masks;
                        }
                        next[u[q]] = Some(v[q]);
                    }
                    let start = u[positions[0]];
                    let mut color = start;
                    for step in 0..positions.len() {
                        match next[color] {
                            Some(c) => color = c,
                            None => continue 'masks,
                        }
                        if color == start && step + 1 < positions.len() {
                            continue 'masks;
                        }
                    }
                    if color != start {
                        continue;
                    }
                    let (mut a, mut b) = (u.clone(), v.clone());
                    for &q in &positions {
                        a[q] = v[q];
                        b[q] = u[q];
                    }
                    if is_member(&a) && is_member(&b) {
                        found = true;
                        break;
                    }
                }
                if !found {
                    return false;
                }
            }
        }
    }
    true
}

fn subset_family(universe: &FixedContentFamily, mask: usize) -> FixedContentFamily {
    let words: Vec<Vec<usize>> = universe
        .words()
        .iter()
        .enumerate()
        .filter(|(i, _)| mask & (1 << i) != 0)
        .map(|(_, w)| w.clone())
        .collect();
    FixedContentFamily::try_new(universe.num_colors(), words).unwrap()
}

fn all_cuts(n: usize) -> Vec<Vec<usize>> {
    (1..(1usize << n) - 1)
        .map(|mask| (0..n).filter(|p| mask & (1 << p) != 0).collect())
        .collect()
}

fn restriction_count(f: &FixedContentFamily, positions: &[usize]) -> usize {
    let set: std::collections::BTreeSet<Vec<usize>> = f
        .words()
        .iter()
        .map(|w| positions.iter().map(|&p| w[p]).collect())
        .collect();
    set.len()
}

fn palette(f: &FixedContentFamily, positions: &[usize]) -> Vec<usize> {
    let set: std::collections::BTreeSet<usize> = f
        .words()
        .iter()
        .flat_map(|w| positions.iter().map(move |&p| w[p]))
        .collect();
    set.into_iter().collect()
}

#[test]
fn quota_211_supports_match_brute_force_and_count_405() {
    let universe = FixedContentFamily::all_assignments(&[2, 1, 1]);
    assert_eq!(universe.len(), 12);
    let mut successes = 0usize;
    for mask in 1usize..1 << 12 {
        let f = subset_family(&universe, mask);
        let expected = brute_force_sce(f.words(), 3);
        let actual = f.check_symmetric_cycle_exchange();
        assert_eq!(actual.is_ok(), expected, "mask {mask}: {actual:?}");
        if expected {
            successes += 1;
            for i in 0..f.len() {
                for j in i + 1..f.len() {
                    for p in 0..4 {
                        if f.word(i)[p] == f.word(j)[p] {
                            continue;
                        }
                        let w = f.exchange_witness(i, j, p).unwrap().expect("witness");
                        assert_eq!(w.positions[0], p);
                        assert_eq!(f.verify_exchange(i, j, p, &w), Ok(()));
                        // Swapping the words reverses the directed cycle.
                        let reversed = CycleWitness {
                            positions: w.positions.iter().rev().copied().collect(),
                        };
                        assert_eq!(f.verify_exchange(j, i, p, &reversed), Ok(()));
                    }
                }
            }
        } else {
            let o = actual.unwrap_err();
            assert_eq!(
                f.exchange_witness(o.first, o.second, o.position).unwrap(),
                None
            );
            assert_eq!(f.word(o.first), o.first_word.as_slice());
            assert_eq!(f.word(o.second), o.second_word.as_slice());
        }
    }
    assert_eq!(successes, 405);
}

#[test]
fn quota_211_cut_certificates_match_elementary_criteria() {
    let universe = FixedContentFamily::all_assignments(&[2, 1, 1]);
    let cuts = all_cuts(4);
    let (mut products, mut binaries) = (0usize, 0usize);
    for mask in 1usize..1 << 12 {
        let f = subset_family(&universe, mask);
        if f.check_symmetric_cycle_exchange().is_err() {
            continue;
        }
        for cut in &cuts {
            let complement: Vec<usize> = (0..4).filter(|p| !cut.contains(p)).collect();
            let is_product =
                restriction_count(&f, cut) * restriction_count(&f, &complement) == f.len();
            let product = f.product_certificate(cut);
            assert_eq!(
                product.is_ok(),
                is_product,
                "mask {mask} cut {cut:?}: {product:?}"
            );
            if let Ok(c) = &product {
                products += 1;
                assert_eq!(c.verify(&f), Ok(()));
            }
            let (pl, pr) = (palette(&f, cut), palette(&f, &complement));
            let shared: Vec<usize> = pl.iter().filter(|c| pr.contains(c)).copied().collect();
            let is_binary = shared.len() == 2 && {
                let counts: Vec<usize> = f
                    .words()
                    .iter()
                    .map(|w| cut.iter().filter(|&&p| w[p] == shared[0]).count())
                    .collect();
                counts.iter().max().unwrap() - counts.iter().min().unwrap() == 1
            };
            let binary = f.binary_port_certificate(cut);
            assert_eq!(
                binary.is_ok(),
                is_binary,
                "mask {mask} cut {cut:?}: {binary:?}"
            );
            if let Ok(c) = &binary {
                binaries += 1;
                assert_eq!(c.verify(&f), Ok(()));
                assert_eq!(c.counted_color, shared[0]);
                assert_eq!(c.other_color, shared[1]);
            }
        }
    }
    assert!(
        products > 0 && binaries > 0,
        "products {products} binaries {binaries}"
    );
}

#[test]
fn binary_matroids_satisfy_two_color_exchange() {
    // Spanning trees of K4 on edges 01,02,03,12,13,23 (color A = in tree).
    let edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
    let mut trees = Vec::new();
    for mask in 0usize..64 {
        if mask.count_ones() != 3 {
            continue;
        }
        let mut parent: Vec<usize> = (0..4).collect();
        fn root(parent: &mut [usize], mut x: usize) -> usize {
            while parent[x] != x {
                x = parent[x];
            }
            x
        }
        let mut acyclic = true;
        for (e, &(a, b)) in edges.iter().enumerate() {
            if mask & (1 << e) != 0 {
                let (ra, rb) = (root(&mut parent, a), root(&mut parent, b));
                if ra == rb {
                    acyclic = false;
                } else {
                    parent[ra] = rb;
                }
            }
        }
        if acyclic {
            trees.push(
                (0..6)
                    .map(|e| usize::from(mask & (1 << e) == 0))
                    .collect::<Vec<_>>(),
            );
        }
    }
    let k4 = FixedContentFamily::try_new(2, trees).unwrap();
    assert_eq!(k4.len(), 16);
    assert_eq!(k4.check_symmetric_cycle_exchange(), Ok(()));
    assert_eq!(k4.color_support(0).unwrap().len(), 16);

    // Fano plane: 3-subsets of 7 points that are not lines.
    let lines = [
        [0, 1, 2],
        [0, 3, 4],
        [0, 5, 6],
        [1, 3, 5],
        [1, 4, 6],
        [2, 3, 6],
        [2, 4, 5],
    ];
    let mut bases = Vec::new();
    for mask in 0usize..128 {
        if mask.count_ones() != 3 {
            continue;
        }
        let points: Vec<usize> = (0..7).filter(|p| mask & (1 << p) != 0).collect();
        if lines.iter().any(|l| l.to_vec() == points) {
            continue;
        }
        bases.push(
            (0..7)
                .map(|p| usize::from(mask & (1 << p) == 0))
                .collect::<Vec<_>>(),
        );
    }
    let fano = FixedContentFamily::try_new(2, bases).unwrap();
    assert_eq!(fano.len(), 28);
    assert_eq!(fano.check_symmetric_cycle_exchange(), Ok(()));
}

#[test]
fn merging_and_splitting_obstructions_from_project_notes() {
    // Merging colors: {ABCD, BCDA} exchanges by a 4-cycle; merged it fails.
    let four = family(4, &["ABCD", "BCDA"]);
    assert_eq!(four.check_symmetric_cycle_exchange(), Ok(()));
    assert_eq!(
        four.exchange_witness(0, 1, 0).unwrap().unwrap().positions,
        vec![0, 1, 2, 3]
    );
    let merged = family(2, &["ABAB", "BABA"]);
    assert!(merged.check_symmetric_cycle_exchange().is_err());

    // Splitting a color: A-positions {1,2}x{3,4} (one-based), B and C free.
    let mut words = Vec::new();
    for a1 in 0..2 {
        for a2 in 2..4 {
            for (b, c) in [(0, 1), (1, 0)] {
                let others: Vec<usize> = (0..4).filter(|&p| p != a1 && p != a2).collect();
                let mut w = vec![0usize; 4];
                w[others[b]] = 1;
                w[others[c]] = 2;
                words.push(w);
            }
        }
    }
    let split = FixedContentFamily::try_new(3, words).unwrap();
    assert_eq!(split.len(), 8);
    let u = split.index_of(&word("ABAC")).unwrap();
    let v = split.index_of(&word("CABA")).unwrap();
    let cycles = simple_exchange_cycles(&word("ABAC"), &word("CABA"), 0).unwrap();
    assert_eq!(
        cycles,
        vec![CycleWitness {
            positions: vec![0, 3]
        }]
    );
    assert_eq!(split.exchange_witness(u, v, 0).unwrap(), None);
    assert!(split.check_symmetric_cycle_exchange().is_err());
}

/// Words `AδBδC`-style: clone positions 0 and 1 of {ABC, BCA} with one shared
/// buffer color δ = 3.
fn shared_buffer_family() -> FixedContentFamily {
    let mut words = Vec::new();
    for base in ["ABC", "BCA"] {
        let base = word(base);
        for i in 0..2 {
            for j in 0..2 {
                let mut w = vec![3usize; 5];
                w[i] = base[0];
                w[2 + j] = base[1];
                w[4] = base[2];
                words.push(w);
            }
        }
    }
    FixedContentFamily::try_new(4, words).unwrap()
}

#[test]
fn shared_buffer_reuse_destroys_exchange() {
    let f = shared_buffer_family();
    assert_eq!(f.len(), 8);
    let (u, v) = (word("ADBDC"), word("DBDCA"));
    let cycles = simple_exchange_cycles(&u, &v, 0).unwrap();
    assert_eq!(
        cycles,
        vec![CycleWitness {
            positions: vec![0, 3, 4]
        }]
    );
    let (a, _) = cycles[0].offspring(&u, &v);
    assert_eq!(a, word("DDBCA"));
    assert!(!f.contains(&a));
    let (i, j) = (f.index_of(&u).unwrap(), f.index_of(&v).unwrap());
    assert_eq!(f.exchange_witness(i, j, 0).unwrap(), None);
    let obstruction = f.check_symmetric_cycle_exchange().unwrap_err();
    assert!(f
        .exchange_witness(obstruction.first, obstruction.second, obstruction.position)
        .unwrap()
        .is_none());
    // The private-block claim is false: every word has two active positions.
    assert!(matches!(
        f.private_block_certificate(&[0, 1, 2, 3], &[3])
            .unwrap_err(),
        CertificateError::ActiveCountNotOne { count: 2, .. }
    ));
}

#[test]
fn witness_mutations_are_rejected_precisely() {
    let f = family(2, &["ABAB", "ABBA", "BAAB", "BABA"]);
    let (u, v) = (
        f.index_of(&word("ABAB")).unwrap(),
        f.index_of(&word("BABA")).unwrap(),
    );
    let ok = CycleWitness {
        positions: vec![0, 1],
    };
    assert_eq!(f.verify_exchange(u, v, 0, &ok), Ok(()));
    assert_eq!(
        f.verify_exchange(
            u,
            v,
            1,
            &CycleWitness {
                positions: vec![1, 0]
            }
        ),
        Ok(())
    );
    assert_eq!(
        f.verify_exchange(u, v, 0, &CycleWitness { positions: vec![0] }),
        Err(WitnessError::TooShort { len: 1 })
    );
    assert_eq!(
        f.verify_exchange(
            u,
            v,
            0,
            &CycleWitness {
                positions: vec![0, 0]
            }
        ),
        Err(WitnessError::RepeatedPosition { position: 0 })
    );
    assert_eq!(
        f.verify_exchange(
            u,
            v,
            0,
            &CycleWitness {
                positions: vec![0, 4]
            }
        ),
        Err(WitnessError::PositionOutOfRange {
            position: 4,
            num_positions: 4
        })
    );
    assert_eq!(
        f.verify_exchange(
            u,
            v,
            0,
            &CycleWitness {
                positions: vec![0, 2]
            }
        ),
        Err(WitnessError::RepeatedColor { color: 0 })
    );
    assert_eq!(
        f.verify_exchange(
            u,
            v,
            0,
            &CycleWitness {
                positions: vec![1, 2]
            }
        ),
        Err(WitnessError::MissingPrescribedPosition { position: 0 })
    );
    assert_eq!(
        f.verify_exchange(
            u,
            v,
            0,
            &CycleWitness {
                positions: vec![0, 3]
            }
        ),
        Err(WitnessError::OffspringMissing { word: word("BBAA") })
    );
    assert_eq!(
        f.verify_exchange(u, u, 0, &ok),
        Err(WitnessError::Query(InvalidQuery::WordsAgree {
            first: u,
            second: u,
            position: 0
        }))
    );

    // Broken cycle identity and a missing offspring.
    let g = FixedContentFamily::all_assignments(&[2, 1, 1]);
    let (a, b) = (
        g.index_of(&word("AABC")).unwrap(),
        g.index_of(&word("ABAC")).unwrap(),
    );
    assert_eq!(
        g.verify_exchange(
            a,
            b,
            1,
            &CycleWitness {
                positions: vec![1, 2]
            }
        ),
        Ok(())
    );
    assert_eq!(
        g.verify_exchange(
            a,
            b,
            1,
            &CycleWitness {
                positions: vec![1, 3]
            }
        ),
        Err(WitnessError::BrokenCycle { index: 0 })
    );
    let (c, d) = (
        g.index_of(&word("AABC")).unwrap(),
        g.index_of(&word("ABCA")).unwrap(),
    );
    assert_eq!(
        g.verify_exchange(
            c,
            d,
            1,
            &CycleWitness {
                positions: vec![1, 3, 2]
            }
        ),
        Err(WitnessError::BrokenCycle { index: 0 })
    );
    assert_eq!(
        g.verify_exchange(
            c,
            d,
            1,
            &CycleWitness {
                positions: vec![1, 2, 3]
            }
        ),
        Ok(())
    );
    let two = family(2, &["ABAB", "BABA"]);
    assert_eq!(
        two.verify_exchange(0, 1, 0, &ok),
        Err(WitnessError::OffspringMissing { word: word("BAAB") })
    );
}

#[test]
fn marginals_handle_zero_weights_and_huge_coefficients() {
    let f = FixedContentFamily::all_assignments(&[2, 1, 1]);
    let ones = vec![vec![BigUint::one(); 3]; 4];
    let counts = f.color_fiber_counts(0).unwrap();
    assert_eq!(counts.len(), 6);
    assert!(counts.values().all(|v| *v == BigUint::from(2u32)));
    assert_eq!(f.weighted_color_marginal(0, &ones).unwrap(), counts);
    let matrix = f.weighted_quadratic_marginal_matrix(0, &ones).unwrap();
    assert_eq!(matrix.len(), 4);
    for (i, row) in matrix.iter().enumerate() {
        for (j, entry) in row.iter().enumerate() {
            let expected = if i == j {
                BigUint::zero()
            } else {
                BigUint::from(2u32)
            };
            assert_eq!(*entry, expected);
        }
    }

    // Zero weight on cell (position 0, color B) removes every word with B there.
    let mut weights = ones.clone();
    weights[0][1] = BigUint::zero();
    let marginal = f.weighted_color_marginal(0, &weights).unwrap();
    assert_eq!(marginal.len(), 6);
    assert_eq!(marginal[&vec![1, 2]], BigUint::one());
    assert_eq!(marginal[&vec![0, 1]], BigUint::from(2u32));
    let mut all_zero = ones.clone();
    for row in &mut all_zero {
        row[1] = BigUint::zero();
    }
    assert!(f.weighted_color_marginal(0, &all_zero).unwrap().is_empty());
    // The target column is ignored.
    let mut target_zero = ones.clone();
    for row in &mut target_zero {
        row[0] = BigUint::zero();
    }
    assert_eq!(f.weighted_color_marginal(0, &target_zero).unwrap(), counts);

    // Huge exact coefficients beyond u128.
    let big = BigUint::from(2u32).pow(100u32);
    let huge = vec![vec![big.clone(); 3]; 4];
    let marginal = f.weighted_color_marginal(0, &huge).unwrap();
    let expected = BigUint::from(2u32) * BigUint::from(2u32).pow(200u32);
    assert!(marginal.values().all(|v| *v == expected));

    assert_eq!(
        f.weighted_color_marginal(3, &ones).unwrap_err(),
        MarginalError::ColorOutOfRange {
            color: 3,
            num_colors: 3
        }
    );
    let mut ragged = ones.clone();
    ragged[2].pop();
    assert_eq!(
        f.weighted_color_marginal(0, &ragged).unwrap_err(),
        MarginalError::WeightRowLength {
            row: 2,
            expected: 3,
            actual: 2
        }
    );
}

#[test]
fn large_assignment_families_exceed_128_words() {
    let f = FixedContentFamily::all_assignments(&[3, 2, 2]);
    assert_eq!(f.len(), 210);
    assert_eq!(f.check_symmetric_cycle_exchange(), Ok(()));
    assert!(f.is_list_quota_saturated());
    let lists = vec![
        vec![0, 1],
        vec![0, 1, 2],
        vec![0, 1, 2],
        vec![0, 1, 2],
        vec![0, 1, 2],
        vec![0, 1, 2],
        vec![0, 1, 2],
    ];
    let g = f.restricted_to_allowed_colors(&lists).unwrap();
    assert_eq!(g.len(), 150);
    assert_eq!(g.check_symmetric_cycle_exchange(), Ok(()));
    assert!(g.is_list_quota_saturated());
    assert_eq!(
        f.restricted_to_allowed_colors(&lists[..3]).unwrap_err(),
        FamilyError::WrongLength {
            what: "allowed color lists",
            expected: 7,
            actual: 3
        }
    );
    let four = FixedContentFamily::all_assignments(&[1, 1, 1, 1]);
    assert_eq!(four.len(), 24);
    assert_eq!(four.check_symmetric_cycle_exchange(), Ok(()));
}

fn interleaved_product() -> FixedContentFamily {
    // {ABC, BCA, CAB} on positions {0, 2, 4} times {AB, BA} on positions {1, 3}.
    let mut words = Vec::new();
    for l in ["ABC", "BCA", "CAB"] {
        for r in ["AB", "BA"] {
            let (l, r) = (word(l), word(r));
            words.push(vec![l[0], r[0], l[1], r[1], l[2]]);
        }
    }
    FixedContentFamily::try_new(3, words).unwrap()
}

#[test]
fn product_certificates_verify_and_reject_corruption() {
    let f = interleaved_product();
    assert_eq!(f.len(), 6);
    let cert = f.product_certificate(&[0, 2, 4]).unwrap();
    assert_eq!(cert.left, family(3, &["ABC", "BCA", "CAB"]));
    assert_eq!(cert.right, family(3, &["AB", "BA"]));
    assert_eq!(cert.verify(&f), Ok(()));
    assert!(matches!(
        f.product_certificate(&[0, 1]).unwrap_err(),
        CertificateError::FactorInvalid {
            factor: Factor::Left,
            error: FamilyError::ContentMismatch { .. }
        }
    ));
    assert!(matches!(
        f.product_certificate(&[]).unwrap_err(),
        CertificateError::InvalidPositions { .. }
    ));
    assert!(matches!(
        f.product_certificate(&[0, 1, 2, 3, 4]).unwrap_err(),
        CertificateError::InvalidPositions { .. }
    ));

    let mut wrong_quota = cert.clone();
    wrong_quota.left = FixedContentFamily::all_assignments(&[2, 1, 0]);
    assert!(matches!(
        wrong_quota.verify(&f).unwrap_err(),
        CertificateError::QuotaMismatch {
            factor: Factor::Left,
            ..
        }
    ));
    let mut omitted = cert.clone();
    omitted.left = family(3, &["ABC", "BCA"]);
    match omitted.verify(&f).unwrap_err() {
        CertificateError::ReconstructionMismatch {
            pairs,
            family_size,
            missing,
            extra,
        } => {
            assert_eq!(
                (pairs, family_size, missing.len(), extra.len()),
                (4, 6, 2, 0)
            );
        }
        other => panic!("unexpected {other:?}"),
    }
    let mut extra = cert.clone();
    extra.left = FixedContentFamily::all_assignments(&[1, 1, 1]);
    assert!(matches!(
        extra.verify(&f).unwrap_err(),
        CertificateError::ReconstructionMismatch { pairs: 12, extra, .. } if extra.len() == 6
    ));
    let mut shape = cert.clone();
    shape.right = family(3, &["ABC", "BCA"]);
    assert!(matches!(
        shape.verify(&f).unwrap_err(),
        CertificateError::FactorShape {
            factor: Factor::Right,
            ..
        }
    ));

    // Exact reconstruction with a factor that fails exchange.
    let bad = family(4, &["ABABCD", "ABABDC", "BABACD", "BABADC"]);
    assert!(matches!(
        bad.product_certificate(&[0, 1, 2, 3]).unwrap_err(),
        CertificateError::FactorNotExchange {
            factor: Factor::Left,
            ..
        }
    ));
    let hand_made = ProductCertificate {
        cut: vec![0, 1, 2, 3],
        left: family(4, &["ABAB", "BABA"]),
        right: family(4, &["CD", "DC"]),
    };
    assert!(matches!(
        hand_made.verify(&bad).unwrap_err(),
        CertificateError::FactorNotExchange {
            factor: Factor::Left,
            ..
        }
    ));
}

#[test]
fn binary_port_certificates_verify_and_reject_corruption() {
    // U_{2,3} across the cut {0}: width one.
    let u23 = family(2, &["AAB", "ABA", "BAA"]);
    let cert = u23.binary_port_certificate(&[0]).unwrap();
    assert_eq!(
        (cert.counted_color, cert.other_color, cert.low_count),
        (0, 1, 0)
    );
    assert_eq!(cert.left, family(2, &["AB", "BA"]));
    assert_eq!(cert.right, u23);
    assert_eq!(cert.verify(&u23), Ok(()));

    // Three colors: {AACB, ABCA, BACA} glues U_{2,3} with {CAB, CBA}.
    let f = family(3, &["AACB", "ABCA", "BACA"]);
    assert_eq!(f.check_symmetric_cycle_exchange(), Ok(()));
    let cert3 = f.binary_port_certificate(&[0, 1]).unwrap();
    assert_eq!(
        (cert3.counted_color, cert3.other_color, cert3.low_count),
        (0, 1, 1)
    );
    assert_eq!(cert3.left, family(3, &["AAB", "ABA", "BAA"]));
    assert_eq!(cert3.right, family(3, &["CAB", "CBA"]));

    // Omitted compatible pair.
    let mut omitted = cert3.clone();
    omitted.right = family(3, &["CAB"]);
    assert!(matches!(
        omitted.verify(&f).unwrap_err(),
        CertificateError::PortLevelMissing {
            factor: Factor::Right,
            color: 0
        }
    ));
    let mut omitted_left = cert3.clone();
    omitted_left.left = family(3, &["AAB", "ABA"]);
    match omitted_left.verify(&f).unwrap_err() {
        CertificateError::ReconstructionMismatch { pairs, missing, .. } => {
            assert_eq!((pairs, missing), (2, vec![word("BACA")]));
        }
        other => panic!("unexpected {other:?}"),
    }
    // Wrong level and swapped shared colors.
    let mut wrong_level = cert3.clone();
    wrong_level.low_count = 0;
    assert!(matches!(
        wrong_level.verify(&f).unwrap_err(),
        CertificateError::PortLevelMismatch {
            factor: Factor::Left,
            ..
        }
    ));
    let mut swapped = cert3.clone();
    std::mem::swap(&mut swapped.counted_color, &mut swapped.other_color);
    // Counting color B (quota one) cannot reach level low_count + 1 = 2.
    assert!(matches!(
        swapped.verify(&f).unwrap_err(),
        CertificateError::CountLevels { color: 1, .. }
    ));
    let mut same = cert3.clone();
    same.other_color = same.counted_color;
    assert!(matches!(
        same.verify(&f).unwrap_err(),
        CertificateError::InvalidColors { .. }
    ));
    // Wrong quota through a foreign factor with matching shape.
    let mut quota = cert3.clone();
    quota.right = family(3, &["ACC", "CAC", "CCA"]);
    assert!(matches!(
        quota.verify(&f).unwrap_err(),
        CertificateError::QuotaMismatch { .. }
    ));
    // Same quota but a port color outside {A, B}.
    let mut foreign = cert3.clone();
    foreign.right = family(3, &["ABC", "BCA", "CAB"]);
    assert!(matches!(
        foreign.verify(&f).unwrap_err(),
        CertificateError::PortColorInvalid {
            factor: Factor::Right,
            ..
        }
    ));
    // Port color outside the shared pair.
    let mut port = cert3.clone();
    port.right = family(3, &["CAB", "BAC"]);
    assert!(matches!(
        port.verify(&f).unwrap_err(),
        CertificateError::QuotaMismatch { .. } | CertificateError::PortColorInvalid { .. }
    ));
    let hand = BinaryPortCertificate {
        cut: vec![0, 1],
        counted_color: 0,
        other_color: 1,
        low_count: 1,
        left: family(3, &["AAB", "ABA", "BAA"]),
        right: family(3, &["ABC", "BAC"]),
    };
    assert!(matches!(
        hand.verify(&f).unwrap_err(),
        CertificateError::PortColorInvalid {
            factor: Factor::Right,
            ..
        }
    ));

    // Non-binary cuts.
    assert!(matches!(
        family(3, &["ABC", "BCA", "CAB"]).binary_port_certificate(&[0]).unwrap_err(),
        CertificateError::SharedPaletteMismatch { shared } if shared == vec![0, 1, 2]
    ));
    let product = family(2, &["ABAB", "ABBA", "BAAB", "BABA"]);
    assert!(matches!(
        product.binary_port_certificate(&[0, 1]).unwrap_err(),
        CertificateError::CountLevels { color: 0, values } if values == vec![1]
    ));
    // Rank-two matroid on four elements missing basis 34: width two at {1,3}.
    let wide = family(2, &["AABB", "ABAB", "ABBA", "BAAB", "BABA"]);
    assert!(matches!(
        wide.binary_port_certificate(&[0, 2]).unwrap_err(),
        CertificateError::CountLevels { color: 0, values } if values == vec![0, 1, 2]
    ));
}

/// Clone position 0 of {ABC, BCA, CAB} into two positions with fresh D.
fn cloned_cyclic_family() -> FixedContentFamily {
    let mut words = Vec::new();
    for base in ["ABC", "BCA", "CAB"] {
        let base = word(base);
        words.push(vec![base[0], 3, base[1], base[2]]);
        words.push(vec![3, base[0], base[1], base[2]]);
    }
    FixedContentFamily::try_new(4, words).unwrap()
}

/// The 24-word family of PRIVATE_BLOCKS.md Section 3: three clone blocks with
/// buffers D, E, F and old-color lists AB, BC, CA inside each block.
fn buffered_cyclic_family() -> FixedContentFamily {
    let allowed: [Vec<usize>; 3] = [vec![0, 2], vec![0, 1], vec![1, 2]];
    let mut words = Vec::new();
    for base in ["ABC", "BCA", "CAB"] {
        let base = word(base);
        let mut choices = vec![vec![]];
        for (block, &color) in base.iter().enumerate() {
            let mut next = Vec::new();
            for prefix in &choices {
                for &active in &allowed[color] {
                    let mut w = prefix.clone();
                    for j in 0..3 {
                        w.push(if j == active { color } else { 3 + block });
                    }
                    next.push(w);
                }
            }
            choices = next;
        }
        words.extend(choices);
    }
    FixedContentFamily::try_new(6, words).unwrap()
}

#[test]
fn private_block_certificates_verify_and_reject_false_claims() {
    let f = cloned_cyclic_family();
    assert_eq!(f.len(), 6);
    assert_eq!(f.quota(), &[1, 1, 1, 1]);
    assert_eq!(f.check_symmetric_cycle_exchange(), Ok(()));
    let cert = f.private_block_certificate(&[0, 1], &[3]).unwrap();
    assert_eq!(cert.alpha(), 4);
    assert_eq!(cert.compressed, family(4, &["ABC", "BCA", "CAB"]));
    assert_eq!(
        cert.internal,
        FixedContentFamily::try_new(5, vec![vec![4, 3], vec![3, 4]]).unwrap()
    );
    let expected: BTreeMap<usize, Vec<usize>> = (0..3).map(|d| (d, vec![0, 1])).collect();
    assert_eq!(cert.allowed_positions, expected);
    assert_eq!(cert.verify(&f), Ok(()));

    let mut narrowed = cert.clone();
    narrowed.allowed_positions.insert(0, vec![0]);
    assert!(matches!(
        narrowed.verify(&f).unwrap_err(),
        CertificateError::ReconstructionMismatch { pairs: 5, missing, .. } if missing == vec![word("DABC")]
    ));
    let mut widened = cert.clone();
    widened.allowed_positions.insert(1, vec![0, 1, 2]);
    assert!(matches!(
        widened.verify(&f).unwrap_err(),
        CertificateError::AllowedPositionsInvalid { color: 1, .. }
    ));
    let mut missing_key = cert.clone();
    missing_key.allowed_positions.remove(&2);
    assert!(matches!(
        missing_key.verify(&f).unwrap_err(),
        CertificateError::AllowedPositionsKeys { .. }
    ));
    let mut wrong_quota = cert.clone();
    wrong_quota.internal = FixedContentFamily::try_new(5, vec![vec![4, 4]]).unwrap();
    assert!(matches!(
        wrong_quota.verify(&f).unwrap_err(),
        CertificateError::QuotaMismatch {
            factor: Factor::Internal,
            ..
        }
    ));
    let mut wrong_compressed = cert.clone();
    wrong_compressed.compressed = FixedContentFamily::all_assignments(&[1, 1, 1, 0]);
    assert!(matches!(
        wrong_compressed.verify(&f).unwrap_err(),
        CertificateError::ReconstructionMismatch { pairs: 12, .. }
    ));
    assert!(matches!(
        f.private_block_certificate(&[0, 1], &[2]).unwrap_err(),
        CertificateError::PrivateColorOutsideBlock { color: 2, .. }
    ));
    assert!(matches!(
        f.private_block_certificate(&[0], &[3]).unwrap_err(),
        CertificateError::PrivateColorOutsideBlock {
            color: 3,
            position: 1,
            ..
        }
    ));
    assert!(matches!(
        f.private_block_certificate(&[0, 1], &[]).unwrap_err(),
        CertificateError::InvalidColors { .. }
    ));
    assert!(matches!(
        f.private_block_certificate(&[0, 1], &[3, 3]).unwrap_err(),
        CertificateError::InvalidColors { .. }
    ));
    assert!(matches!(
        f.private_block_certificate(&[0, 1, 2], &[3]).unwrap_err(),
        CertificateError::ActiveCountNotOne { count: 2, .. }
    ));

    // Twenty-four words, no product or binary cut, three private blocks.
    let g = buffered_cyclic_family();
    assert_eq!(g.len(), 24);
    assert_eq!(g.quota(), &[1, 1, 1, 2, 2, 2]);
    assert_eq!(g.check_symmetric_cycle_exchange(), Ok(()));
    for cut in all_cuts(9) {
        assert!(g.product_certificate(&cut).is_err(), "product cut {cut:?}");
        assert!(
            g.binary_port_certificate(&cut).is_err(),
            "binary cut {cut:?}"
        );
    }
    for block in 0..3 {
        let positions: Vec<usize> = (3 * block..3 * block + 3).collect();
        let cert = g
            .private_block_certificate(&positions, &[3 + block])
            .unwrap();
        assert_eq!(cert.compressed.len(), 12);
        assert_eq!(cert.internal.len(), 3);
        let expected: BTreeMap<usize, Vec<usize>> =
            [(0, vec![0, 2]), (1, vec![0, 1]), (2, vec![1, 2])]
                .into_iter()
                .collect();
        assert_eq!(cert.allowed_positions, expected);
        assert_eq!(cert.verify(&g), Ok(()));
        let mut unrestricted = cert.clone();
        for list in unrestricted.allowed_positions.values_mut() {
            *list = vec![0, 1, 2];
        }
        assert!(matches!(
            unrestricted.verify(&g).unwrap_err(),
            CertificateError::ReconstructionMismatch { pairs: 36, extra, .. } if extra.len() == 12
        ));
    }
}

#[test]
fn review_regressions_for_budgets_sparse_labels_and_certificate_padding() {
    use combinatoric_core::colored_exchange::{ExchangeStatus, WitnessSearch};

    let product = family(2, &["ABAB", "ABBA", "BAAB", "BABA"]);
    // Four adjacent pairs cost one extension each. Each opposite pair
    // costs three: one unsuccessful completion, then two successful ones.
    assert_eq!(
        product.check_symmetric_cycle_exchange_bounded(10),
        ExchangeStatus::Satisfied
    );
    assert!(matches!(
        product.check_symmetric_cycle_exchange_bounded(9),
        ExchangeStatus::BudgetExhausted { steps: 9, .. }
    ));
    assert!(matches!(
        product.exchange_witness_bounded(0, 3, 2, 1).unwrap(),
        WitnessSearch::BudgetExhausted { steps: 1 }
    ));
    assert!(matches!(
        product.exchange_witness_bounded(0, 3, 2, 2).unwrap(),
        WitnessSearch::Found(_)
    ));
    assert!(product.exchange_witness(0, 3, 2).unwrap().is_some());

    let cycles = simple_exchange_cycles(&[usize::MAX, 7], &[7, usize::MAX], 0).unwrap();
    assert_eq!(
        cycles,
        vec![CycleWitness {
            positions: vec![0, 1]
        }]
    );

    let f = family(4, &["ADBC", "BDCA", "CDAB"]);
    let cert = f.private_block_certificate(&[0, 1], &[3]).unwrap();
    let mut padded = cert.clone();
    padded.internal = FixedContentFamily::try_new(5, vec![vec![4, 3], vec![3, 4]]).unwrap();
    assert_eq!(
        padded.verify(&f),
        Err(CertificateError::UnusedPrivateBlockData)
    );
    let mut extra_allowed = cert;
    extra_allowed.allowed_positions.get_mut(&0).unwrap().push(1);
    assert_eq!(
        extra_allowed.verify(&f),
        Err(CertificateError::UnusedPrivateBlockData)
    );

    let binary = family(2, &["AAB", "ABA", "BAA"]);
    let mut corrupt = binary.binary_port_certificate(&[0]).unwrap();
    corrupt.low_count = usize::MAX;
    assert!(matches!(
        corrupt.verify(&binary),
        Err(CertificateError::CountLevels { .. })
    ));
}
