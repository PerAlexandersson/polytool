//! Positive and negative certification, support alignment and rational differential tests.
use super::*;
use crate::primes::DescendingPrimes;
use num_traits::One;
fn q(value: i64) -> BigRational {
    BigRational::from_integer(value.into())
}
fn p(vars: usize, terms: &[(Vec<u32>, i64)]) -> RationalPolynomial {
    RationalPolynomial::from_terms(
        vars,
        terms.iter().map(|(m, c)| (m.clone(), q(*c))).collect(),
    )
}
fn primes(values: &[u64]) -> Vec<ValidatedPrime> {
    values
        .iter()
        .map(|&p| ValidatedPrime::new(p).unwrap())
        .collect()
}
fn linear(constant: BigRational) -> RationalPolynomial {
    RationalPolynomial::from_terms(1, BTreeMap::from([(vec![1], q(1)), (vec![0], constant)]))
}

#[test]
fn align_vanishing_supports_and_skip_bad_denominators() {
    let input = vec![linear(q(2))];
    let result =
        certified_modular_groebner_basis(&input, MonomialOrder::Lex, &primes(&[2, 5, 7])).unwrap();
    assert_eq!(result.basis, input);
    assert_eq!(result.certificate.reconstruction_primes, vec![2, 5]);
    assert_eq!(result.certificate.verification_prime, 7);
    let input = vec![linear(BigRational::new(1.into(), 2.into()))];
    let result =
        certified_modular_groebner_basis(&input, MonomialOrder::Lex, &primes(&[2, 3, 5, 7]))
            .unwrap();
    assert_eq!(result.basis, input);
    assert_eq!(result.certificate.skipped_denominator_primes, vec![2]);
    assert_eq!(result.certificate.reconstruction_primes, vec![3, 5]);
}

#[test]
fn failed_spare_is_not_accepted_and_new_support_is_aligned() {
    let input = vec![linear(q(6))];
    let result =
        certified_modular_groebner_basis(&input, MonomialOrder::Lex, &primes(&[2, 3, 5, 7, 11]))
            .unwrap();
    assert_eq!(result.basis, input);
    assert_eq!(result.certificate.failed_spare_checks, 1);
    assert_eq!(result.certificate.verification_prime, 11);
    assert!(!result.certificate.reconstruction_primes.contains(&11));
    assert_eq!(result.certificate.reconstruction_modulus, BigInt::from(210));
}

#[test]
fn unlucky_early_profiles_do_not_pin_the_answer() {
    let input = vec![p(1, &[(vec![1], 6)])];
    let result =
        certified_modular_groebner_basis(&input, MonomialOrder::Lex, &primes(&[2, 3, 5, 7, 11]))
            .unwrap();
    assert_eq!(result.basis, vec![linear(q(0))]);
    assert_eq!(result.certificate.other_profile_primes, vec![2, 3]);
    assert_eq!(result.certificate.reconstruction_primes, vec![5, 7]);
    assert_eq!(result.certificate.verification_prime, 11);
}

#[test]
fn spare_agreement_cannot_certify_a_false_unit_ideal() {
    let input = vec![p(1, &[(vec![1], 30), (vec![0], -1)])];
    let result = certified_modular_groebner_basis(
        &input,
        MonomialOrder::Lex,
        &primes(&[2, 3, 5, 7, 11, 13, 17, 19, 23]),
    )
    .unwrap();
    assert_eq!(
        result.basis,
        reduced_groebner_basis(&input, MonomialOrder::Lex)
    );
    assert!(result.certificate.failed_exact_checks >= 1);
    assert_eq!(
        check_rational_groebner_basis(
            &input,
            &[RationalPolynomial::monomial(1, vec![0], q(1))],
            MonomialOrder::Lex
        ),
        Err(GroebnerCertificateError::CandidateNotContained { index: 0 })
    );
}

#[test]
fn exact_verifier_rejects_every_missing_obligation() {
    let x = p(2, &[(vec![1, 0], 1)]);
    let y = p(2, &[(vec![0, 1], 1)]);
    assert_eq!(
        check_rational_groebner_basis(
            std::slice::from_ref(&x),
            std::slice::from_ref(&y),
            MonomialOrder::Lex
        ),
        Err(GroebnerCertificateError::InputNotContained { index: 0 })
    );
    let not_groebner = vec![
        p(2, &[(vec![2, 0], 1), (vec![0, 1], -1)]),
        p(2, &[(vec![1, 1], 1), (vec![0, 0], -1)]),
    ];
    assert_eq!(
        check_rational_groebner_basis(&not_groebner, &not_groebner, MonomialOrder::Lex),
        Err(GroebnerCertificateError::NonzeroSPolynomial {
            first: 0,
            second: 1
        })
    );
    assert_eq!(
        check_rational_groebner_basis(&[x], &[RationalPolynomial::zero(2)], MonomialOrder::Lex),
        Err(GroebnerCertificateError::ZeroBasisElement)
    );
    assert_eq!(
        check_rational_groebner_basis(&[y], &[RationalPolynomial::zero(1)], MonomialOrder::Lex),
        Err(GroebnerCertificateError::IncompatibleVariableCounts)
    );
}

#[test]
fn zero_unit_zero_variables_and_exhaustion() {
    let ps = primes(&[101, 103, 107]);
    for input in [
        vec![],
        vec![RationalPolynomial::zero(2)],
        vec![RationalPolynomial::monomial(0, vec![], q(2))],
        vec![linear(q(0)), linear(q(1))],
    ] {
        let result = certified_modular_groebner_basis(&input, MonomialOrder::Lex, &ps).unwrap();
        assert_eq!(
            result.basis,
            reduced_groebner_basis(&input, MonomialOrder::Lex)
        );
    }
    assert_eq!(
        certified_modular_groebner_basis(&[linear(q(1))], MonomialOrder::Lex, &ps[..2]),
        Err(CertifiedModularGroebnerError::PrimeBudgetExhausted { primes_examined: 2 })
    );
    assert_eq!(
        certified_modular_groebner_basis(&[], MonomialOrder::Lex, &primes(&[101, 103, 101])),
        Err(CertifiedModularGroebnerError::RepeatedPrime { prime: 101 })
    );
    assert_eq!(
        certified_modular_groebner_basis(
            &[RationalPolynomial::zero(1), RationalPolynomial::zero(2)],
            MonomialOrder::Lex,
            &ps
        ),
        Err(CertifiedModularGroebnerError::IncompatibleVariableCounts)
    );
    assert_eq!(
        certified_modular_groebner_basis(
            &[linear(BigRational::new(1.into(), (101 * 103 * 107).into()))],
            MonomialOrder::Lex,
            &ps
        ),
        Err(CertifiedModularGroebnerError::PrimeBudgetExhausted { primes_examined: 3 })
    );
}

#[test]
fn arbitrary_size_moduli_and_word_prime_arithmetic() {
    let ps: Vec<_> = DescendingPrimes::new(u64::MAX).take(6).collect();
    let input = vec![linear(BigRational::from_integer(BigInt::one() << 80))];
    let result = certified_modular_groebner_basis(&input, MonomialOrder::Lex, &ps).unwrap();
    assert_eq!(result.basis, input);
    assert!(result.certificate.reconstruction_modulus.bits() > 128);
    assert!(!result
        .certificate
        .reconstruction_primes
        .contains(&result.certificate.verification_prime));
    let prime = ps[0];
    let a = WordPrimeCoefficient::Residue {
        value: prime.get() - 1,
        prime,
    };
    let one = WordPrimeField::one();
    assert!(WordPrimeField::is_zero(&(a.clone() + one.clone())));
    assert_eq!((a.clone() * a.clone()).residue(prime), 1);
    assert_eq!(
        WordPrimeField::divide(one, a).residue(prime),
        prime.get() - 1
    );
    assert!(WordPrimeField::is_zero(&WordPrimeField::zero()));
    // Context-free identities remain exact before being bound to a prime.
    assert_eq!(
        (WordPrimeField::one() + WordPrimeField::one()).residue(prime),
        2
    );
}

#[test]
fn seeded_modular_results_equal_rational_engine_in_all_orders() {
    let ps = primes(&[101, 103, 107, 109, 113, 127, 131, 137, 139]);
    let mut seed = 0x37a4_6677_bb21_9acdu64;
    let mut next = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        ((seed >> 32) % 7) as i64 - 3
    };
    for _ in 0..32 {
        let input = vec![
            p(
                2,
                &[(vec![2, 0], 1), (vec![0, 1], next()), (vec![0, 0], next())],
            ),
            p(
                2,
                &[(vec![1, 1], 1), (vec![1, 0], next()), (vec![0, 0], next())],
            ),
        ];
        for order in MonomialOrder::STANDARD_ORDERS {
            let rational = reduced_groebner_basis(&input, order);
            let modular = certified_modular_groebner_basis(&input, order, &ps).unwrap();
            assert_eq!(modular.basis, rational, "order={order:?} input={input:?}");
            assert_eq!(
                check_rational_groebner_basis(&input, &modular.basis, order).unwrap(),
                modular.certificate.exact_verification
            );
        }
    }
}

#[test]
fn a_bad_candidate_denominator_at_the_spare_is_rejected() {
    let input = vec![linear(q(7))];
    let result =
        certified_modular_groebner_basis(&input, MonomialOrder::Lex, &primes(&[3, 5, 2, 7, 11]))
            .unwrap();
    assert_eq!(result.basis, input);
    assert!(result.certificate.failed_spare_checks >= 1);
    assert_eq!(result.certificate.verification_prime, 11);
    assert_eq!(
        result.certificate.skipped_denominator_primes,
        Vec::<u64>::new()
    );
}

#[test]
fn bad_input_denominator_does_not_consume_a_pending_spare_check() {
    let input = vec![linear(BigRational::new(1.into(), 2.into()))];
    let result =
        certified_modular_groebner_basis(&input, MonomialOrder::Lex, &primes(&[3, 5, 2, 7]))
            .unwrap();
    assert_eq!(result.basis, input);
    assert_eq!(result.certificate.skipped_denominator_primes, vec![2]);
    assert_eq!(result.certificate.reconstruction_primes, vec![3, 5]);
    assert_eq!(result.certificate.verification_prime, 7);
    assert_eq!(result.certificate.failed_spare_checks, 0);
}
