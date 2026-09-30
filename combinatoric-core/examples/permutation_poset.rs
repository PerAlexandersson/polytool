//! Public dominance-poset constructor; permutation values are one-indexed,
//! while the naturally labeled poset uses zero-indexed positions.
use combinatoric_core::poset::Poset;

fn main() {
    let perm = [2, 1, 4, 3];
    let p = Poset::from_permutation(&perm);
    let height = p.chain_polynomial().len() - 1;
    let profile = p.p_eulerian_polynomial();
    assert!(p.is_naturally_labeled());
    assert_eq!(height, 2);
    assert_eq!(profile, [1, 2, 1]);
    assert_eq!(profile, p.order_polytope_hstar());
    println!(
        "permutation={perm:?}, covers={:?}, height={height}, W={profile:?}",
        p.covers()
    );
}
