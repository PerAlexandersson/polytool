//! Certify a rational basis reconstructed from word primes.
use exact_algebra::groebner::{
    certified_modular_groebner_basis, check_rational_groebner_basis, MonomialOrder, Polynomial,
    RationalField, SparsePolynomial,
};
use exact_algebra::primes::ValidatedPrime;
use num_rational::BigRational;
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut polynomial =
        Polynomial::<RationalField>::monomial(1, vec![1], BigRational::from_integer(1.into()));
    polynomial.add_term(vec![0], &BigRational::from_integer(2.into()));
    let input = vec![polynomial.clone()];
    let primes = [2, 5, 7]
        .into_iter()
        .map(ValidatedPrime::new)
        .collect::<Result<Vec<_>, _>>()?;
    let result = certified_modular_groebner_basis(&input, MonomialOrder::Lex, &primes)?;
    assert_eq!(result.basis, vec![polynomial]);
    assert_eq!(
        check_rational_groebner_basis(&input, &result.basis, MonomialOrder::Lex)?,
        result.certificate.exact_verification
    );
    println!("certified basis: [x + 2]");
    println!(
        "reconstruction primes: {:?}",
        result.certificate.reconstruction_primes
    );
    println!(
        "unused verification prime: {}",
        result.certificate.verification_prime
    );
    println!(
        "exact checks: {} S-pairs, {} input containment, {} reverse containment",
        result.certificate.exact_verification.s_pairs_checked,
        result
            .certificate
            .exact_verification
            .input_generators_checked,
        result
            .certificate
            .exact_verification
            .candidate_generators_checked
    );
    Ok(())
}
