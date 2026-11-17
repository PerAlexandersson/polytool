//! Checked rational sparse Gröbner basis and normal-form example.

use exact_algebra::groebner::{
    GroebnerBasis, MonomialOrder, Polynomial, RationalField, SparsePolynomial,
};
use num_rational::BigRational;

fn term(exponents: &[u32], coefficient: i64) -> Polynomial<RationalField> {
    Polynomial::monomial(
        exponents.len(),
        exponents.to_vec(),
        BigRational::from_integer(coefficient.into()),
    )
}

fn main() {
    let x = term(&[1, 0], 1);
    let y = term(&[0, 1], 1);
    let basis = GroebnerBasis::new(vec![x.clone(), x.clone() + y.clone()], MonomialOrder::Lex);
    assert_eq!(basis.generators, vec![x.clone(), y.clone()]);
    assert!(basis.normal_form(&(x + y)).is_zero());
    println!("reduced basis: [x, y]");
    println!("normal form of x + y: 0");
}
