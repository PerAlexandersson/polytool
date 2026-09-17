use combinatoric_core::{Graph, Partition};
use sym_poly_sym::chordal_clique_boundary_schur_channels;

fn main() {
    let graph = Graph::path(4);
    let channels = chordal_clique_boundary_schur_channels(&graph, &[2, 3])
        .expect("the terminal edge of a path is a valid chordal boundary");

    let trivial = Partition::new(vec![2]);
    let sign = Partition::new(vec![1, 1]);
    assert!(channels[&trivial].is_zero());
    println!("interior {trivial}: {}", channels[&trivial]);
    println!("interior {sign}: {}", channels[&sign]);
}
