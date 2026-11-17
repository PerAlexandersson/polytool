//! Reuse a rational lexicographic subtraction budget across calls.
use exact_algebra::groebner::{
    reduced_lexicographic_groebner_basis_with_budget, GroebnerBudget, Polynomial, RationalField,
    SparsePolynomial,
};
use num_rational::BigRational;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let x = Polynomial::<RationalField>::monomial(1, vec![1], BigRational::from_integer(1.into()));
    let input = [x.clone(), x.clone()];
    let mut budget = GroebnerBudget {
        remaining_steps: 0,
        max_coefficient_bits: 4096,
        max_basis_elements: 400,
    };
    let error = reduced_lexicographic_groebner_basis_with_budget(&input, &mut budget)
        .expect_err("ordered input reduction needs one subtraction");
    println!("zero-step error: {error}");
    budget.remaining_steps = 1;
    let basis = reduced_lexicographic_groebner_basis_with_budget(&input, &mut budget)?;
    assert_eq!(basis, vec![x]);
    println!(
        "one-step basis: [x]; remaining steps: {}",
        budget.remaining_steps
    );
    Ok(())
}
