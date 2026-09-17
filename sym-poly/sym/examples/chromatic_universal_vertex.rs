use combinatoric_core::Graph;
use sym_poly_sym::{chromatic_symmetric, chromatic_symmetric_add_universal_vertex};

fn main() {
    let path = Graph::path(4);
    let path_chromatic = chromatic_symmetric::<i64>(&path);
    let from_symmetric_function =
        chromatic_symmetric_add_universal_vertex(&path_chromatic.to_schur_basis());
    let from_graph = chromatic_symmetric::<i64>(&path.join(&Graph::complete(1)));

    assert_eq!(from_symmetric_function, from_graph);
    println!("{from_symmetric_function}");
}
