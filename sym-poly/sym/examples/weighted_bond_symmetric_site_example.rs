use combinatoric_core::Graph;
use num_bigint::BigInt;
use sym_poly_core::Partition;
use sym_poly_sym::{
    chromatic_mobius_symmetric_function, graph_mobius_symmetric_function, SymmetricFunction,
};

type Sym = SymmetricFunction<BigInt>;

fn partition(parts: &[u32]) -> Partition {
    Partition::new(parts.to_vec())
}

fn assert_coefficient(function: &Sym, parts: &[u32], expected: i64) {
    assert_eq!(
        function.coefficient(&partition(parts)),
        BigInt::from(expected)
    );
}

fn main() {
    let path = Graph::path(3);
    let complete = Graph::complete(3);

    let m_path = graph_mobius_symmetric_function::<BigInt>(&path)
        .unwrap()
        .to_monomial_basis();
    let m_complete = graph_mobius_symmetric_function::<BigInt>(&complete)
        .unwrap()
        .to_monomial_basis();
    assert_coefficient(&m_path, &[2], 1);
    assert_coefficient(&m_path, &[1, 1], 3);
    assert_coefficient(&m_complete, &[2], 2);
    assert_coefficient(&m_complete, &[1, 1], 5);

    let psi_path = chromatic_mobius_symmetric_function::<BigInt>(&path).to_elementary_basis();
    let psi_complete =
        chromatic_mobius_symmetric_function::<BigInt>(&complete).to_elementary_basis();
    for (parts, expected) in [
        (&[][..], 1),
        (&[1][..], -2),
        (&[2][..], 1),
        (&[1, 1][..], 1),
    ] {
        assert_coefficient(&psi_path, parts, expected);
    }
    for (parts, expected) in [
        (&[][..], 1),
        (&[1][..], -3),
        (&[2][..], 1),
        (&[1, 1][..], 2),
    ] {
        assert_coefficient(&psi_complete, parts, expected);
    }

    println!("M_P3 in the m-basis: {m_path}");
    println!("M_K3 in the m-basis: {m_complete}");
    println!("Psi_P3 in the e-basis: {psi_path}");
    println!("Psi_K3 in the e-basis: {psi_complete}");
}
