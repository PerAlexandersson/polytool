//! Modular reconstruction with an independent exact ideal-equality check.

use super::{
    normal_form, reduced_groebner_basis, s_polynomial, CoefficientRing, ExactField, MonomialOrder,
    Polynomial, RationalField, SparsePolynomial,
};
use crate::chinese_remainder::ChineseRemainderAccumulator;
use crate::modular_arithmetic::{bigint_residue_u64, multiply_modulo_u64};
use crate::primes::ValidatedPrime;
use crate::rational_reconstruction::rational_reconstruction;
use num_bigint::BigInt;
use num_rational::BigRational;
use std::collections::{BTreeMap, BTreeSet};

mod coefficients;
use coefficients::{WordPrimeCoefficient, WordPrimeField};
type RationalPolynomial = Polynomial<RationalField>;
type ModularBasis = Vec<BTreeMap<Vec<u32>, u64>>;
type Profile = Vec<Vec<u32>>;

/// Failure of an exact rational Gröbner certificate check.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GroebnerCertificateError {
    /// Input or candidate polynomials have incompatible ambient rings.
    IncompatibleVariableCounts,
    /// A candidate entry is zero and cannot be a basis element.
    ZeroBasisElement,
    /// A candidate S-polynomial has nonzero normal form.
    NonzeroSPolynomial {
        /// Zero-based first candidate index.
        first: usize,
        /// Zero-based second candidate index.
        second: usize,
    },
    /// An input generator does not reduce to zero by the candidate basis.
    InputNotContained {
        /// Zero-based input generator index.
        index: usize,
    },
    /// A candidate basis entry does not belong to the input ideal.
    CandidateNotContained {
        /// Zero-based candidate basis index.
        index: usize,
    },
}
impl std::fmt::Display for GroebnerCertificateError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::IncompatibleVariableCounts => f.write_str("polynomial variable counts differ"),
            Self::ZeroBasisElement => f.write_str("candidate basis contains zero"),
            Self::NonzeroSPolynomial { first, second } => write!(
                f,
                "candidate S-pair ({first}, {second}) does not reduce to zero"
            ),
            Self::InputNotContained { index } => {
                write!(f, "input generator {index} is not in the candidate ideal")
            }
            Self::CandidateNotContained { index } => {
                write!(f, "candidate generator {index} is not in the input ideal")
            }
        }
    }
}
impl std::error::Error for GroebnerCertificateError {}

/// Counts of exact checks performed over arbitrary-precision rationals.
///
/// These record successful computation, rather than a serialized proof object.
/// Replay with [`check_rational_groebner_basis`]; reverse containment uses an
/// independent rational basis computed from the original generators.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExactGroebnerVerification {
    /// Number of candidate S-pairs reduced to zero.
    pub s_pairs_checked: usize,
    /// Number of original generators reduced to zero by the candidate.
    pub input_generators_checked: usize,
    /// Number of candidate entries reduced to zero by the exact input basis.
    pub candidate_generators_checked: usize,
}

/// Check Buchberger's criterion and both rational ideal containments exactly.
///
/// In particular, `[1]` is rejected for a proper input ideal even though every
/// input reduces to zero by `[1]`. This check recomputes an independent rational
/// basis of the input ideal and is unbounded; it is not a modular speedup.
/// Empty input denotes the zero ideal in the candidate's ambient ring.
///
/// # Panics
///
/// Panics if a monomial exponent or lcm total degree overflows `u32`.
pub fn check_rational_groebner_basis(
    generators: &[RationalPolynomial],
    candidate: &[RationalPolynomial],
    order: MonomialOrder,
) -> Result<ExactGroebnerVerification, GroebnerCertificateError> {
    let num_vars = generators
        .first()
        .or(candidate.first())
        .map(SparsePolynomial::num_vars)
        .unwrap_or(0);
    if generators
        .iter()
        .chain(candidate)
        .any(|p| p.num_vars() != num_vars)
    {
        return Err(GroebnerCertificateError::IncompatibleVariableCounts);
    }
    if candidate.iter().any(SparsePolynomial::is_zero) {
        return Err(GroebnerCertificateError::ZeroBasisElement);
    }
    let mut s_pairs_checked = 0;
    for first in 0..candidate.len() {
        for second in first + 1..candidate.len() {
            if !normal_form(
                &s_polynomial(&candidate[first], &candidate[second], order),
                candidate,
                order,
            )
            .is_zero()
            {
                return Err(GroebnerCertificateError::NonzeroSPolynomial { first, second });
            }
            s_pairs_checked += 1;
        }
    }
    for (index, generator) in generators.iter().enumerate() {
        if !normal_form(generator, candidate, order).is_zero() {
            return Err(GroebnerCertificateError::InputNotContained { index });
        }
    }
    let reference = reduced_groebner_basis(generators, order);
    for (index, p) in candidate.iter().enumerate() {
        if !normal_form(p, &reference, order).is_zero() {
            return Err(GroebnerCertificateError::CandidateNotContained { index });
        }
    }
    Ok(ExactGroebnerVerification {
        s_pairs_checked,
        input_generators_checked: generators.len(),
        candidate_generators_checked: candidate.len(),
    })
}

/// Failure to reconstruct and certify a modular rational basis.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CertifiedModularGroebnerError {
    /// Original generators have different variable counts.
    IncompatibleVariableCounts,
    /// Every supplied prime, including possible spare primes, must be distinct.
    RepeatedPrime {
        /// Repeated prime modulus.
        prime: u64,
    },
    /// The supplied finite prime list ended before a candidate was certified.
    PrimeBudgetExhausted {
        /// Total supplied primes examined, including denominator skips and failed spares.
        primes_examined: usize,
    },
}
impl std::fmt::Display for CertifiedModularGroebnerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::IncompatibleVariableCounts => {
                f.write_str("generators have different variable counts")
            }
            Self::RepeatedPrime { prime } => write!(f, "prime {prime} occurs more than once"),
            Self::PrimeBudgetExhausted { primes_examined } => write!(
                f,
                "no certified basis after examining {primes_examined} primes"
            ),
        }
    }
}
impl std::error::Error for CertifiedModularGroebnerError {}

/// Modular agreement and exact verification evidence for a returned rational basis.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ModularGroebnerCertificate {
    /// At least two distinct primes used in the final coefficient reconstruction.
    pub reconstruction_primes: Vec<u64>,
    /// Ordered leading exponent vectors shared by the reconstruction and spare bases.
    pub leading_monomials: Vec<Vec<u32>>,
    /// Product of the final reconstruction primes, with arbitrary precision.
    pub reconstruction_modulus: BigInt,
    /// Distinct unused prime checking the full reconstructed basis and leading profile.
    pub verification_prime: u64,
    /// Input-denominator primes skipped before certification, in encounter order.
    pub skipped_denominator_primes: Vec<u64>,
    /// Encountered primes with a different leading profile from the returned basis.
    pub other_profile_primes: Vec<u64>,
    /// Reconstructed candidates rejected by an independent prime comparison.
    pub failed_spare_checks: usize,
    /// Candidates passing a spare comparison but failing exact ideal verification.
    pub failed_exact_checks: usize,
    /// Successful exact S-pair and both-containment checks.
    pub exact_verification: ExactGroebnerVerification,
}

/// A reduced rational basis together with its successful modular and exact checks.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CertifiedModularGroebnerBasis {
    /// Monic reduced basis sorted by decreasing leading monomial in the requested order.
    pub basis: Vec<RationalPolynomial>,
    /// Evidence for profile agreement, unused spare-prime agreement and exact ideal equality.
    pub certificate: ModularGroebnerCertificate,
}

struct ProfileGroup {
    profile: Profile,
    images: Vec<(ValidatedPrime, ModularBasis)>,
}
struct Reconstruction {
    basis: Vec<RationalPolynomial>,
    group: usize,
    modulus: BigInt,
}

/// Reconstruct a modular rational basis and certify its exact ideal over `Q`.
///
/// The caller supplies a finite ordered list of distinct validated word primes.
/// Input-denominator primes are skipped. Bases are grouped by ordered leading
/// profiles; reconstruction requires two matching primes. Supports are aligned
/// by union, inserting zero coefficients where needed. CRT moduli are arbitrary
/// precision. Failed rational reconstructions request another prime. Each candidate
/// is compared with the full basis at an unused spare prime, then checked exactly
/// with [`check_rational_groebner_basis`]. A failed spare is eligible for a later
/// reconstruction, but the successful spare never enters its reconstruction.
/// Different profiles are retained separately, so an early unlucky prime cannot
/// pin the computation to that profile. At least three usable primes are needed.
///
/// Only the prime count is bounded. Individual modular computations and the
/// independent rational certificate are unbounded and intended for small ideals.
/// Certification recomputes a rational basis and may cost more than solving directly.
///
/// # Panics
///
/// Panics if a monomial exponent or total lcm degree overflows `u32`.
pub fn certified_modular_groebner_basis(
    generators: &[RationalPolynomial],
    order: MonomialOrder,
    primes: &[ValidatedPrime],
) -> Result<CertifiedModularGroebnerBasis, CertifiedModularGroebnerError> {
    let num_vars = generators
        .first()
        .map(SparsePolynomial::num_vars)
        .unwrap_or(0);
    if generators.iter().any(|p| p.num_vars() != num_vars) {
        return Err(CertifiedModularGroebnerError::IncompatibleVariableCounts);
    }
    let mut seen = BTreeSet::new();
    for prime in primes {
        if !seen.insert(prime.get()) {
            return Err(CertifiedModularGroebnerError::RepeatedPrime { prime: prime.get() });
        }
    }
    let mut groups = Vec::<ProfileGroup>::new();
    let mut pending: Option<Reconstruction> = None;
    let mut skipped = Vec::new();
    let mut failed_spare_checks = 0;
    let mut failed_exact_checks = 0;
    for &prime in primes {
        let Some(image) = modular_basis(generators, order, prime) else {
            skipped.push(prime.get());
            continue;
        };
        let profile = leading_profile(&image, order);
        if let Some(candidate) = pending.take() {
            let group = &groups[candidate.group];
            if profile == group.profile
                && reduce_basis(&candidate.basis, prime).as_ref() == Some(&image)
            {
                match check_rational_groebner_basis(generators, &candidate.basis, order) {
                    Ok(exact_verification) => {
                        let certificate = ModularGroebnerCertificate {
                            reconstruction_primes: group
                                .images
                                .iter()
                                .map(|(p, _)| p.get())
                                .collect(),
                            leading_monomials: profile,
                            reconstruction_modulus: candidate.modulus,
                            verification_prime: prime.get(),
                            skipped_denominator_primes: skipped,
                            other_profile_primes: groups
                                .iter()
                                .filter(|g| g.profile != group.profile)
                                .flat_map(|g| g.images.iter().map(|(p, _)| p.get()))
                                .collect(),
                            failed_spare_checks,
                            failed_exact_checks,
                            exact_verification,
                        };
                        return Ok(CertifiedModularGroebnerBasis {
                            basis: candidate.basis,
                            certificate,
                        });
                    }
                    Err(_) => failed_exact_checks += 1,
                }
            } else {
                failed_spare_checks += 1;
            }
        }
        let group = if let Some(index) = groups.iter().position(|g| g.profile == profile) {
            index
        } else {
            groups.push(ProfileGroup {
                profile,
                images: Vec::new(),
            });
            groups.len() - 1
        };
        groups[group].images.push((prime, image));
        if groups[group].images.len() >= 2 {
            pending =
                reconstruct(&groups[group], num_vars).map(|(basis, modulus)| Reconstruction {
                    basis,
                    modulus,
                    group,
                });
        }
    }
    Err(CertifiedModularGroebnerError::PrimeBudgetExhausted {
        primes_examined: primes.len(),
    })
}

fn reduce_coefficient(value: &BigRational, prime: ValidatedPrime) -> Option<u64> {
    Some(multiply_modulo_u64(
        bigint_residue_u64(value.numer(), prime.get()),
        prime.inverse(bigint_residue_u64(value.denom(), prime.get()))?,
        prime.get(),
    ))
}
fn reduce_basis(basis: &[RationalPolynomial], prime: ValidatedPrime) -> Option<ModularBasis> {
    basis
        .iter()
        .map(|p| {
            p.terms()
                .iter()
                .map(|(m, c)| Some((m.clone(), reduce_coefficient(c, prime)?)))
                .collect::<Option<BTreeMap<_, _>>>()
                .map(|mut terms| {
                    terms.retain(|_, v| *v != 0);
                    terms
                })
        })
        .collect()
}
fn modular_basis(
    generators: &[RationalPolynomial],
    order: MonomialOrder,
    prime: ValidatedPrime,
) -> Option<ModularBasis> {
    let reduced = reduce_basis(generators, prime)?;
    let input: Vec<_> = generators
        .iter()
        .zip(reduced)
        .map(|(p, terms)| {
            Polynomial::<WordPrimeField>::from_terms(
                p.num_vars(),
                terms
                    .into_iter()
                    .map(|(m, value)| (m, WordPrimeCoefficient::Residue { value, prime }))
                    .collect(),
            )
        })
        .collect();
    Some(
        reduced_groebner_basis(&input, order)
            .into_iter()
            .map(|p| {
                p.terms()
                    .iter()
                    .map(|(m, c)| (m.clone(), c.residue(prime)))
                    .collect()
            })
            .collect(),
    )
}
fn leading_profile(basis: &ModularBasis, order: MonomialOrder) -> Profile {
    basis
        .iter()
        .map(|terms| {
            terms
                .keys()
                .max_by(|a, b| order.compare(a, b))
                .expect("modular basis contains no zero entries")
                .clone()
        })
        .collect()
}
fn reconstruct(group: &ProfileGroup, num_vars: usize) -> Option<(Vec<RationalPolynomial>, BigInt)> {
    let supports: Vec<BTreeSet<Vec<u32>>> = (0..group.profile.len())
        .map(|index| {
            group
                .images
                .iter()
                .flat_map(|(_, basis)| basis[index].keys().cloned())
                .collect()
        })
        .collect();
    let width = supports.iter().map(BTreeSet::len).sum();
    let mut crt = ChineseRemainderAccumulator::new(width);
    for (prime, basis) in &group.images {
        let values: Vec<_> = supports
            .iter()
            .enumerate()
            .flat_map(|(i, support)| {
                support
                    .iter()
                    .map(move |m| basis[i].get(m).copied().unwrap_or(0))
            })
            .collect();
        crt.combine_prime(&values, *prime)
            .expect("validated distinct reconstruction primes and aligned supports");
    }
    let values: Vec<_> = crt
        .residues()
        .iter()
        .map(|v| rational_reconstruction(v, crt.modulus()))
        .collect::<Result<_, _>>()
        .ok()?;
    let mut coefficients = values.into_iter();
    let basis = supports
        .into_iter()
        .map(|support| {
            RationalPolynomial::from_terms(
                num_vars,
                support
                    .into_iter()
                    .map(|m| (m, coefficients.next().expect("aligned coefficient count")))
                    .collect(),
            )
        })
        .collect();
    Some((basis, crt.modulus().clone()))
}

#[cfg(test)]
mod tests;
