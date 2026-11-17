//! Independent parity sections: exact allocations and symmetric-square identities.
//!
//! For factors satisfying `f_(m+2)=beta*f_(m+1)-mu^2*f_m`, the state
//! `(s_m*t_m, s_m*t_(m+1)+s_(m+1)*t_m, s_(m+1)*t_(m+1))` has matrix
//! `[[0,0,1],[0,-mu^2,2*beta],[mu^4,-beta*mu^2,beta^2]]`.
//! Its characteristic polynomial is
//! `z^3-(beta^2-mu^2)z^2+mu^2*(beta^2-mu^2)z-mu^6`.
//! The source's stride-two companion has exactly that characteristic polynomial.
//! Cayley--Hamilton and three seed products therefore prove equality for all
//! indices, including repeated characteristic roots. No row fitting is used.

use super::*;
use std::collections::BTreeMap;

/// Exact product identity for one independently allocated parity section.
/// Coefficient vectors are ascending rational strings.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IndependentParitySectionProductCertificate {
    /// Certified stride-two source section; residues are relative to its anchor.
    pub subsequence: index_section::ConstantCoefficientIndexSubsequenceCertificate,
    /// Initial section rows omitted before the product identity starts.
    pub skipped_rows: usize,
    /// Three identifying products, starting after the omitted section rows.
    pub seed_products: Vec<Vec<String>>,
    /// Monic divisors of the first and second identifying rows, respectively.
    pub initial_divisors: [Vec<String>; 2],
    /// Nonzero rational scalar multiplying the second left-factor seed.
    pub initial_scalar_ratio: String,
    /// Numbers of monic divisors of the first and second identifying rows.
    /// Their product is at most [`MAXIMUM_PARITY_PRODUCT_DIVISORS`].
    pub divisor_lattice_sizes: [usize; 2],
    /// Divisor pairs examined through the selected pair (one-based).
    pub allocations_examined: usize,
    /// Symmetric-square characteristic polynomial, ascending in `z`, then `x`.
    pub product_characteristic_polynomial: Vec<Vec<String>>,
}

fn characteristic(beta: &Poly, mu: &Poly) -> Vec<Poly> {
    let trace = beta.pow(2) - mu.pow(2);
    vec![-mu.pow(6), mu.pow(2) * trace.clone(), -trace, Poly::one()]
}

// A rational scalar t must annihilate every x-coefficient of
// A(x)t^2+B(x)t+C(x). Their gcd has degree at most two in t.
fn scalar_roots(a: &Poly, b: &Poly, c: &Poly) -> Vec<Q> {
    let length = [a, b, c]
        .iter()
        .map(|p| p.coeffs().len())
        .max()
        .unwrap_or(0);
    let mut gcd = Poly::zero();
    for i in 0..length {
        gcd = gcd.gcd(&Poly::new(vec![
            vec_poly::coeff(c.coeffs(), i),
            vec_poly::coeff(b.coeffs(), i),
            vec_poly::coeff(a.coeffs(), i),
        ]));
    }
    let mut roots = match gcd.degree() {
        None => vec![Q::one()],
        Some(1) => vec![-gcd.coeffs()[0].clone() / &gcd.coeffs()[1]],
        Some(2) => {
            let [c, b, a] = <&[Q; 3]>::try_from(gcd.coeffs()).expect("quadratic gcd");
            let discriminant = b * b - Q::from_integer(4.into()) * a * c;
            let Some(root) =
                index_section::polynomial_power_root(&Poly::new(vec![discriminant]), 2)
            else {
                return Vec::new();
            };
            let root = vec_poly::coeff(root.coeffs(), 0);
            let denominator = Q::from_integer(2.into()) * a;
            vec![(-b - &root) / &denominator, (-b + root) / denominator]
        }
        _ => Vec::new(),
    };
    roots.retain(|r| !r.is_zero());
    roots.sort();
    roots.dedup();
    roots
}

type FactorCache = BTreeMap<Vec<Vec<Q>>, Option<RealRootednessTheoremReport>>;

fn parse_polynomial(values: &[String]) -> Result<Vec<Q>, RecurrenceError> {
    values
        .iter()
        .map(|value| {
            parse_rational_coeff(value)
                .map_err(|error| RecurrenceError::InvalidInput(error.to_string()))
        })
        .collect()
}

fn cached_factor(
    cache: &mut FactorCache,
    rule: &Recurrence,
    seeds: Vec<Vec<Q>>,
    options: &RealRootednessTheoremOptions,
) -> Result<Option<RealRootednessTheoremReport>, RecurrenceError> {
    if let Some(proof) = cache.get(&seeds) {
        return Ok(proof.clone());
    }
    let mut proof = prove_factor(rule, &seeds, options)?;
    if proof.is_none() {
        // This existing lift's parent menu excludes parity products and itself.
        // It also covers quadratic marker substitutions without new root rules.
        let attempt = check_with_first_proving_start(
            RealRootednessTheorem::IndexSectionLift,
            rule,
            0,
            &seeds,
            options,
        )?;
        if attempt.is_proved()
            && attempt
                .conclusion
                .as_ref()
                .is_some_and(|c| c.first_index == 0 && !c.has_zero_residue_classes())
        {
            proof = Some(attempt);
        }
    }
    cache.insert(seeds, proof.clone());
    Ok(proof)
}

fn allocate(
    section: &index_section::ConstantCoefficientIndexSubsequenceCertificate,
    beta: &Poly,
    mu: &Poly,
    options: &RealRootednessTheoremOptions,
) -> Result<
    Option<(
        IndependentParitySectionProductCertificate,
        Vec<RealRootednessTheoremReport>,
    )>,
    RecurrenceError,
> {
    let (rule, first, seeds) = section.subsequence.to_recurrence_parts()?;
    let factor_rule = factor_rule(beta, mu);
    let mut cache = FactorCache::new();
    for skipped in 0..=options.maximum_start_delay {
        let rows = rule
            .generate_rows_rational(&seeds, first, skipped + 3)
            .map_err(|error| RecurrenceError::Unsupported(error.to_string()))?;
        let identifying = &rows[skipped..skipped + 3];
        if earlier_seeds_hypothesis(identifying, 0, options.root_region.as_ref())
            .is_some_and(|h| h.status != LiftHypothesisStatus::Verified)
        {
            continue;
        }
        let p: Vec<Poly> = identifying.iter().cloned().map(Poly::new).collect();
        let (Some(left_divisors), Some(right_divisors)) = (divisors(&p[0]), divisors(&p[1])) else {
            continue;
        };
        let sizes = [left_divisors.len(), right_divisors.len()];
        if sizes[0]
            .checked_mul(sizes[1])
            .is_none_or(|n| n > MAXIMUM_PARITY_PRODUCT_DIVISORS)
        {
            continue;
        }
        for (i, left_zero) in left_divisors.iter().enumerate() {
            let right_zero = quotient(&p[0], left_zero).expect("a seed divisor divides");
            for (j, left_one_monic) in right_divisors.iter().enumerate() {
                let right_one_monic =
                    quotient(&p[1], left_one_monic).expect("a seed divisor divides");
                let a = beta.clone() * mu.pow(2) * left_one_monic.clone() * right_zero.clone();
                let b = p[2].clone() - beta.pow(2) * p[1].clone() - mu.pow(4) * p[0].clone();
                let c = beta.clone() * mu.pow(2) * left_zero.clone() * right_one_monic.clone();
                for scalar in scalar_roots(&a, &b, &c) {
                    let left_one = left_one_monic.scale(&scalar);
                    let right_one = right_one_monic.scale(&scalar.recip());
                    let left_two = beta.clone() * left_one.clone() - mu.pow(2) * left_zero.clone();
                    let right_two =
                        beta.clone() * right_one.clone() - mu.pow(2) * right_zero.clone();
                    if left_two * right_two != p[2] {
                        continue;
                    }
                    let Some(left) = cached_factor(
                        &mut cache,
                        &factor_rule,
                        vec![
                            vec_poly::trim(left_zero.coeffs()),
                            vec_poly::trim(left_one.coeffs()),
                        ],
                        options,
                    )?
                    else {
                        continue;
                    };
                    let Some(right) = cached_factor(
                        &mut cache,
                        &factor_rule,
                        vec![
                            vec_poly::trim(right_zero.coeffs()),
                            vec_poly::trim(right_one.coeffs()),
                        ],
                        options,
                    )?
                    else {
                        continue;
                    };
                    return Ok(Some((
                        IndependentParitySectionProductCertificate {
                            subsequence: section.clone(),
                            skipped_rows: skipped,
                            seed_products: p.iter().map(strings).collect(),
                            initial_divisors: [strings(left_zero), strings(left_one_monic)],
                            initial_scalar_ratio: scalar.to_string(),
                            divisor_lattice_sizes: sizes,
                            allocations_examined: i * sizes[1] + j + 1,
                            product_characteristic_polynomial: characteristic(beta, mu)
                                .iter()
                                .map(strings)
                                .collect(),
                        },
                        vec![left, right],
                    )));
                }
            }
        }
    }
    Ok(None)
}

pub(super) fn independent_sections_report(
    mut report: RealRootednessTheoremReport,
    rule: &Recurrence,
    first: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
    mu: &Poly,
    beta: &Poly,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    // This failure reason is intrinsic, independent of the start-search budget.
    // In particular, zero seed rows cannot support a nonzero-row conclusion.
    let failure = "the seed Casoratian compatibility identity fails";
    // Leave marker-free scalar families to the existing routes. In particular,
    // do not change scalar-only intermediate searches through this extension.
    if beta.degree().unwrap_or(0) == 0
        && mu.degree().unwrap_or(0) == 0
        && seeds
            .iter()
            .all(|row| row.iter().skip(1).all(Zero::is_zero))
    {
        return Ok(report.not_applicable(failure));
    }
    if earlier_seeds_hypothesis(seeds, first, options.root_region.as_ref())
        .is_some_and(|h| h.status != LiftHypothesisStatus::Verified)
    {
        return Ok(report.not_applicable(failure));
    }
    let source = RecurrenceJson::from_recurrence_rational(rule, first, seeds, None);
    let mut sections = Vec::new();
    let mut factor_proofs = Vec::new();
    for residue in 0..2 {
        let section =
            index_section::certify_constant_coefficient_index_subsequence(&source, 2, residue)?;
        let Some((certificate, proofs)) = allocate(&section, beta, mu, options)? else {
            return Ok(report.not_applicable(failure));
        };
        sections.push(certificate);
        factor_proofs.extend(proofs);
    }
    let covered = sections
        .iter()
        .map(|s| {
            s.skipped_rows
                .checked_mul(2)
                .and_then(|n| s.subsequence.first_source_index.checked_add(n))
        })
        .collect::<Option<Vec<_>>>()
        .ok_or_else(|| RecurrenceError::InvalidInput("parity-section index overflow".into()))?;
    let prefix = covered
        .into_iter()
        .max()
        .expect("two sections")
        .saturating_sub(1)
        - first;
    let rows = rule
        .generate_rows_rational(seeds, first, seeds.len().max(prefix))
        .map_err(|error| RecurrenceError::Unsupported(error.to_string()))?;
    report.reference.push_str("; independent sections: Cayley-Hamilton for stride-two and symmetric-square companions, recurrence uniqueness and polynomial root union");
    report.normal_form = format!("each parity section is an independent product of two families satisfying f_(m+2)=({beta})*f_(m+1)-({mu})^2*f_m; no consecutive interlacing transfer");
    report.hypotheses =
        earlier_seeds_hypothesis(&rows[..prefix], first, options.root_region.as_ref())
            .into_iter()
            .collect();
    report.hypotheses.push(exact_hypothesis("independent_section_product_identities",
        "both stride-two and symmetric-square characteristic polynomials agree, with three exact identifying products per section",
        Ok("Cayley-Hamilton and recurrence uniqueness; exact identities over Q[x]".into())));
    report.hypotheses.push(exact_hypothesis(
        "factor_proofs",
        "all four factor families have replayable proofs from zero in the requested root region",
        Ok(factor_proofs
            .iter()
            .map(|p| p.theorem.as_str())
            .collect::<Vec<_>>()
            .join(", ")),
    ));
    report.certificate = Some(RouteCertificate::ParityProduct(Box::new(
        ParityProductCertificate {
            mu: strings(mu),
            beta: strings(beta),
            exceptional_rows: prefix,
            initial_common_divisor: Vec::new(),
            divisors_examined: sections.iter().map(|s| s.allocations_examined).sum(),
            divisor_lattice_size: sections
                .iter()
                .map(|s| s.divisor_lattice_sizes.iter().product::<usize>())
                .sum(),
            factor_proofs,
            maximum_start_delay: options.maximum_start_delay,
            proof_depth: 1,
            independent_sections: Some(sections),
        },
    )));
    Ok(report.finish(RealRootednessTheoremConclusion {
        first_index: first, zero_residue_classes: None, interlacing_from: None,
        interlacing: None, simple_roots: false,
        root_region: options.root_region.as_ref().map(open_root_interval),
        statement: format!("every P_n, n >= {first}, is nonzero and real-rooted by independent parity-section products; no simplicity or consecutive interlacing is asserted"),
    }))
}

pub(super) fn verify_sections(
    certificate: &ParityProductCertificate,
) -> Result<(), RecurrenceError> {
    let invalid = || {
        RecurrenceError::VerificationFailed(
            "invalid independent parity-section identity or allocation".into(),
        )
    };
    let sections = certificate
        .independent_sections
        .as_ref()
        .ok_or_else(invalid)?;
    if sections.len() != 2 || !certificate.initial_common_divisor.is_empty() {
        return Err(invalid());
    }
    let beta = Poly::new(parse_polynomial(&certificate.beta)?);
    let mu = Poly::new(parse_polynomial(&certificate.mu)?);
    let expected = characteristic(&beta, &mu);
    for (i, section) in sections.iter().enumerate() {
        section.subsequence.verify()?;
        let size = section.divisor_lattice_sizes[0]
            .checked_mul(section.divisor_lattice_sizes[1])
            .ok_or_else(invalid)?;
        if section.subsequence.stride != 2
            || section.subsequence.residue != i
            || section.skipped_rows > certificate.maximum_start_delay
            || size == 0
            || size > MAXIMUM_PARITY_PRODUCT_DIVISORS
            || section.allocations_examined == 0
            || section.allocations_examined > size
            || section.seed_products.len() != 3
        {
            return Err(invalid());
        }
        for recorded in [
            &section.product_characteristic_polynomial,
            &section.subsequence.characteristic_polynomial,
        ] {
            let actual = recorded
                .iter()
                .map(|p| parse_polynomial(p).map(Poly::new))
                .collect::<Result<Vec<_>, _>>()?;
            if actual != expected {
                return Err(invalid());
            }
        }
        let scalar = section
            .initial_scalar_ratio
            .parse::<Q>()
            .map_err(|_| invalid())?;
        if scalar.is_zero() {
            return Err(invalid());
        }
        let (rule, first, seeds) = section.subsequence.subsequence.to_recurrence_parts()?;
        let rows = rule
            .generate_rows_rational(&seeds, first, section.skipped_rows + 3)
            .map_err(|_| invalid())?;
        let mut factors = Vec::new();
        for proof in &certificate.factor_proofs[2 * i..2 * i + 2] {
            let (rule, first, seeds) = proof.recurrence.to_recurrence_parts()?;
            let expected_rule =
                RecurrenceJson::from_recurrence_rational(&factor_rule(&beta, &mu), 0, &[], None);
            if proof.recurrence.recurrence != expected_rule.recurrence
                || first != 0
                || seeds.len() != 2
            {
                return Err(invalid());
            }
            factors.push(
                rule.generate_rows_rational(&seeds, first, 3)
                    .map_err(|_| invalid())?,
            );
        }
        for j in 0..3 {
            let product = Poly::new(factors[0][j].clone()) * Poly::new(factors[1][j].clone());
            if product != Poly::new(rows[section.skipped_rows + j].clone())
                || product != Poly::new(parse_polynomial(&section.seed_products[j])?)
            {
                return Err(invalid());
            }
        }
        let divisors = section
            .initial_divisors
            .iter()
            .map(|p| parse_polynomial(p).map(Poly::new))
            .collect::<Result<Vec<_>, _>>()?;
        if Poly::new(factors[0][0].clone()) != divisors[0]
            || Poly::new(factors[0][1].clone()) != divisors[1].scale(&scalar)
        {
            return Err(invalid());
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn relative_scalars_are_exact_including_degenerate_equations() {
        let p = |c: &[i64]| Poly::new(c.iter().map(|v| Q::from_integer((*v).into())).collect());
        assert_eq!(
            scalar_roots(&p(&[1, 1]), &p(&[-5, -5]), &p(&[6, 6])),
            vec![Q::from_integer(2.into()), Q::from_integer(3.into())]
        );
        assert!(scalar_roots(&p(&[1]), &p(&[0]), &p(&[-2])).is_empty());
        assert!(scalar_roots(&p(&[1]), &p(&[0]), &p(&[1])).is_empty());
        assert!(scalar_roots(&p(&[0]), &p(&[0]), &p(&[1])).is_empty());
        assert_eq!(scalar_roots(&p(&[0]), &p(&[0]), &p(&[0])), vec![Q::one()]);
        assert_eq!(
            scalar_roots(&p(&[0]), &p(&[3]), &p(&[-2])),
            vec![Q::new(2.into(), 3.into())]
        );
        assert!(scalar_roots(&p(&[1, 1]), &p(&[-5, -6]), &p(&[6, 7])).is_empty());
    }
}
