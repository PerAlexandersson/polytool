//! Published and geometric basis-map preservation criteria.

use super::*;

// ---------------------------------------------------------------------------
// Published transformations
// ---------------------------------------------------------------------------

/// The magic-cone lemma used by the Eulerian and derangement
/// transformations (module documentation of `basis_map`, "Published
/// transformations").
const MAGIC_CONE_LEMMA: &str = "an input p of degree d with every root in [-1, 0] and a positive \
     leading coefficient is a nonnegative combination of u^k (1 + u)^(d - k), 0 <= k <= d: with \
     p = c prod_i (u - x_i), sum_k c_k s^k = (1 - s)^d p(s / (1 - s)) = c prod_i ((1 + x_i) s - \
     x_i) has nonnegative coefficients";

/// Recognize an even-derivative family and decide its all-degree hypotheses.
pub(super) fn even_derivative(
    map: &BasisMap,
    options: &BasisMapTheoremOptions,
) -> Result<Option<Analysis>, RecurrenceError> {
    let certificate = match super::super::even_derivative::recognize(map) {
        Ok(Some(certificate)) => certificate,
        Ok(None)
        | Err(RecurrenceError::BudgetExhausted(_))
        | Err(RecurrenceError::Unsupported(_)) => return Ok(None),
        Err(error) => return Err(error),
    };
    let mut analysis = Analysis {
        shape: BasisMapShape::EvenDerivative {
            certificate: Box::new(certificate.clone()),
        },
        status: BasisMapTheoremStatus::HypothesisFailed,
        theorem: "RealRooted.evenDerivative_preservesPF_tfae: c_0=1, c_k>=0 and every monomial \
                  image H_N real-rooted imply T_c preserves PF in every degree; finite support \
                  uses Hermite--Poulain (splits_evenDerivativeTest_of_splits_symbol) and \
                  coefficientwise closure of real-rootedness for the converse"
            .into(),
        hypotheses: vec![
            hypothesis(
                "even_derivative_basis_identity",
                "the basis image of u^N is sum_k c_k (N)_(2k) x^(N-k) for every N>=0".into(),
                certificate.basis_identity.is_certified(),
                certificate.basis_identity.statement.clone(),
            ),
            hypothesis(
                "even_derivative_nonnegative_coefficients",
                "c_0=1 and c_k>=0 for every k>=0".into(),
                certificate.nonnegative_coefficients,
                "exact integer-range ratio sign decision or finite coefficient list".into(),
            ),
            hypothesis(
                "even_derivative_all_rows_pf",
                "H_N is PF for every N>=0".into(),
                certificate.all_rows_pf,
                "certified hypergeometric row identities and discovery on H_N or its two \
                 reversed sections, without catalog lifts; or exact real-rootedness of the \
                 fixed signed symbol C(-w^2)"
                    .into(),
            ),
        ],
        conclusions: Vec::new(),
        counterexample: None,
        lean: vec![
            "RealRooted.evenDerivative_preservesPF_tfae",
            "RealRooted.isPFPolynomial_evenDerivativeOperator",
            "RealRooted.evenDerivativeOperator_X_pow",
        ],
    };
    // Failure of a bounded proof search is not a disproof of a hypothesis.
    if !certificate.basis_identity.is_certified() {
        analysis.hypotheses[0].status = LiftHypothesisStatus::Unresolved;
    }
    if !certificate.all_rows_pf {
        analysis.hypotheses[2].status = LiftHypothesisStatus::Unresolved;
    }
    if certificate.signed_symbol.is_some() {
        analysis
            .lean
            .push("RealRooted.splits_evenDerivativeTest_of_splits_symbol");
        if certificate.basis_identity.is_certified() && !certificate.all_rows_pf {
            // Negative coefficients contradict PF for the corresponding
            // monomial. Otherwise the signed-symbol converse rules out PF
            // preservation in all degrees, without needing a small witness.
            analysis.status = BasisMapTheoremStatus::NotPreserving;
            analysis.hypotheses[2].status = LiftHypothesisStatus::Failed { witness: None };
        }
    }
    if certificate.basis_identity.is_certified()
        && certificate.nonnegative_coefficients
        && certificate.all_rows_pf
    {
        analysis.push_conclusion(BasisMapConclusion {
            input_interval: AlgebraicInterval::nonpositive_ray(),
            output_interval: AlgebraicInterval::nonpositive_ray(),
            nonnegative_coefficients: true,
            sharp: false,
            statement: "T_c sends every PF polynomial to a PF polynomial in every degree \
                        (nonzero inputs stay nonzero because c_0=1)"
                .into(),
            hypotheses: Vec::new(),
        });
    } else if certificate.basis_identity.is_certified() {
        for n in 0..options.witness_search_limit.min(17) {
            let image = map.basis_rows(n + 1)?.pop().expect("row");
            if !AlgebraicInterval::nonpositive_ray().contains_all_roots_of(&image) {
                let mut input = vec![Q::zero(); n + 1];
                input[n] = Q::one();
                analysis.counterexample = Some(counterexample(
                    &input,
                    &image,
                    format!("T_c(u^{n}) is not PF (exact root isolation)"),
                ));
                analysis.status = BasisMapTheoremStatus::NotPreserving;
                analysis.hypotheses[2].status = LiftHypothesisStatus::Failed {
                    witness: Some(n.into()),
                };
                break;
            }
        }
    }
    Ok(Some(analysis))
}

/// The analysis of a map recognised as a published family, `None` when the
/// first rows match none.
pub(super) fn published(map: &BasisMap) -> Result<Option<Analysis>, RecurrenceError> {
    let Some((family, certificate)) = recognize_published_family(map)? else {
        return Ok(None);
    };
    let (shape, theorem, rows, lean): (BasisMapShape, &str, String, Vec<&'static str>) =
        match family {
            PublishedFamily::Narayana(m) => (
                BasisMapShape::Narayana { parameter: m },
                "Mao and Wang, The Narayana transformation, arXiv:2607.01572, Theorem 1: for \
                 every integer m >= 0 the map u^k -> N_(k,m)(t) = sum_j C(k,j) C(k+m,j) / \
                 C(m+j,j) t^j sends every polynomial with nonnegative coefficients and only real \
                 roots to a polynomial with only real nonpositive roots (PF to PF)",
                format!("N_(k,{m})(t)"),
                vec![
                    "RealRooted.narayanaTransformPreservesPF",
                    "RealRooted.narayanaPolynomial_deriv_lag_rec",
                    "RealRooted.narayanaPolynomial_one",
                ],
            ),
            PublishedFamily::EulerianTransformation => (
                BasisMapShape::EulerianTransformation,
                "Athanasiadis, On the real-rootedness of the Eulerian transformation, \
                 arXiv:2302.00754 (J. London Math. Soc. 111 (2025) e70083), Theorem 1: A° \
                 sends every nonnegative combination of u^k (1 + u)^(d - k), 0 <= k <= d, to a \
                 real-rooted polynomial with nonnegative coefficients (Brenti's stronger \
                 conjecture for all inputs with nonpositive roots is false, Brändén and \
                 Jochemko; it is not used)",
                "1, t A_1(t), t A_2(t), ... (t A_k(t) = sum over permutations of [k] of \
                 t^(1 + des))"
                    .into(),
                vec![
                    "RealRooted.EulerianTransform.transform_magicExpansion_isPF",
                    "RealRooted.EulerianTransform.transform_magicExpansion_interl",
                    "RealRooted.generalizedEulerian",
                ],
            ),
            PublishedFamily::DerangementTransformation => (
                BasisMapShape::DerangementTransformation,
                "Brändén and Solus, arXiv:1808.04141, Corollary 3.7 (as stated by Athanasiadis, \
                 arXiv:2302.00754, Section 4): the derangement transformation u^k -> d_k(t) \
                 sends every nonnegative combination of u^k (1 + u)^(d - k), 0 <= k <= d, to \
                 zero or a real-rooted polynomial with nonnegative coefficients",
                "the derangement polynomials d_k(t) by excedances (d_0 = 1, d_1 = 0)".into(),
                vec![
                    "RealRooted.DerangementTransform.transform_magicExpansion_isPF",
                    "RealRooted.DerangementTransform.polynomial_succ",
                ],
            ),
        };
    let mut analysis = Analysis {
        shape,
        status: BasisMapTheoremStatus::Proved,
        theorem: theorem.into(),
        hypotheses: Vec::new(),
        conclusions: Vec::new(),
        counterexample: None,
        lean,
    };
    let identified = certificate.is_certified();
    analysis.hypotheses.push(hypothesis(
        "basis_is_published_family",
        format!(
            "the basis row of u^k is {rows} for every k >= 0 (P_({} + k))",
            map.first_index()
        ),
        identified,
        if identified {
            format!(
                "all-index recurrence equality with the RealRooted recurrence ({}): {}",
                certificate.method.as_deref().unwrap_or("certified"),
                certificate.statement
            )
        } else {
            format!(
                "the first rows agree but the all-index equality is not certified: {}",
                certificate.reasons.join("; ")
            )
        },
    ));
    if !identified {
        return Ok(Some(analysis));
    }
    let nonpositive = AlgebraicInterval::nonpositive_ray();
    match family {
        PublishedFamily::Narayana(_) => analysis.push_conclusion(BasisMapConclusion {
            input_interval: nonpositive.clone(),
            output_interval: nonpositive,
            nonnegative_coefficients: true,
            sharp: false,
            statement: "every input with all roots in (-inf, 0] (a PF polynomial up to sign) is \
                        sent to zero or a polynomial with all roots in (-inf, 0], PF for a \
                        positive leading coefficient"
                .into(),
            hypotheses: Vec::new(),
        }),
        PublishedFamily::EulerianTransformation | PublishedFamily::DerangementTransformation => {
            analysis.push_conclusion(BasisMapConclusion {
                input_interval: AlgebraicInterval::closed(
                    QuadraticNumber::rational(-Q::one()),
                    QuadraticNumber::rational(Q::zero()),
                ),
                output_interval: nonpositive,
                nonnegative_coefficients: true,
                sharp: false,
                statement: format!(
                    "every input with all roots in [-1, 0] is sent to zero or a polynomial with \
                     all roots in (-inf, 0], PF for a positive leading coefficient ({MAGIC_CONE_LEMMA})"
                ),
                hypotheses: Vec::new(),
            })
        }
    }
    Ok(Some(analysis))
}

// ---------------------------------------------------------------------------
// Geometric families
// ---------------------------------------------------------------------------

pub(super) fn geometric(g: &[Q], y: &[Q]) -> Result<Analysis, RecurrenceError> {
    let a = vec_poly::coeff(y, 0);
    let d = vec_poly::coeff(y, 1);
    let mut analysis = Analysis {
        shape: BasisMapShape::Geometric {
            initial_row: poly_text(g),
            constant: text(&a),
            slope: text(&d),
        },
        status: BasisMapTheoremStatus::Proved,
        theorem: "geometric family P_k = g(t) (a + d t)^k: T(p) = g(t) p(a + d t), so the roots \
                  of T(p) are those of g and the points (r - a)/d for the roots r of p \
                  (elementary)"
            .into(),
        hypotheses: Vec::new(),
        conclusions: Vec::new(),
        counterexample: None,
        lean: Vec::new(),
    };
    let nonzero = !vec_poly::is_zero(g);
    analysis.hypotheses.push(hypothesis(
        "initial_row_nonzero",
        "g = P_0 is nonzero".into(),
        nonzero,
        "exact".into(),
    ));
    if !nonzero {
        return Ok(analysis);
    }
    let real = is_real_rooted(g);
    analysis.hypotheses.push(hypothesis(
        "initial_row_real_rooted",
        "g is real-rooted".into(),
        real,
        "exact root isolation".into(),
    ));
    if !real {
        analysis.status = BasisMapTheoremStatus::NotPreserving;
        analysis.counterexample = Some(counterexample(
            &[Q::one()],
            g,
            "T(1) = g is not real-rooted".into(),
        ));
        return Ok(analysis);
    }
    if !d.is_zero() {
        analysis.push_conclusion(BasisMapConclusion {
            input_interval: AlgebraicInterval::whole_line(),
            output_interval: AlgebraicInterval::whole_line(),
            nonnegative_coefficients: false,
            sharp: true,
            statement: "every real-rooted input has a real-rooted image (an affine substitution \
                        times g)"
                .into(),
            hypotheses: Vec::new(),
        });
    }
    let nonpositive = AlgebraicInterval::nonpositive_ray().contains_all_roots_of(g);
    analysis.hypotheses.push(hypothesis(
        "initial_row_nonpositive_roots",
        "every root of g is <= 0".into(),
        nonpositive,
        "exact root isolation".into(),
    ));
    if !nonpositive {
        return Ok(analysis);
    }
    let positive_leading = leading(g).is_positive();
    let (input, nonnegative, sharp, statement) = match d.cmp(&Q::zero()) {
        Ordering::Greater => (
            AlgebraicInterval::at_most(rational(&a)),
            positive_leading,
            true,
            format!(
                "an input root r <= {} gives the root (r - a)/d <= 0; a larger root gives a \
                 positive root",
                text(&a)
            ),
        ),
        Ordering::Less => (
            AlgebraicInterval::at_least(rational(&a)),
            false,
            true,
            format!(
                "an input root r >= {} gives the root (r - a)/d <= 0",
                text(&a)
            ),
        ),
        Ordering::Equal => (
            AlgebraicInterval::at_most(rational(&a)),
            positive_leading,
            false,
            "T(p) = p(a) g; an input with every root <= a has p(a) of the sign of its leading \
             coefficient"
                .into(),
        ),
    };
    analysis.push_conclusion(BasisMapConclusion {
        input_interval: input,
        output_interval: AlgebraicInterval::nonpositive_ray(),
        nonnegative_coefficients: nonnegative,
        sharp,
        statement,
        hypotheses: Vec::new(),
    });
    Ok(analysis)
}

/// An eventually geometric family with a nonzero prefix discrepancy:
/// no positive-length source interval (RealRooted
/// `exists_not_splits_linearMap_shiftedPower_of_eventually_geometric`).
pub(super) fn eventually_geometric(map: &BasisMap, y: &[Q]) -> Result<Analysis, RecurrenceError> {
    let a = vec_poly::coeff(y, 0);
    let d = vec_poly::coeff(y, 1);
    let seeds = map.seeds();
    let k = seeds.len() - 1;
    let mut analysis = Analysis {
        shape: BasisMapShape::EventuallyGeometric {
            geometric_from: k,
            constant: text(&a),
            slope: text(&d),
        },
        status: BasisMapTheoremStatus::HypothesisFailed,
        theorem: "eventually geometric families P_n = P_K y^(n - K) (n >= K >= 1, y = a + d t, \
                  d != 0, P_K != 0) with a nonzero discrepancy E_j = y^K P_j - P_K y^j, j < K: \
                  for every r != 0 some T((u + r)^N), max(K, D + 3) <= N <= max(K, D + 3) + K - \
                  1 (D the largest degree of a nonzero E_j), is not real-rooted, so no source \
                  interval of positive length works"
            .into(),
        hypotheses: Vec::new(),
        conclusions: Vec::new(),
        counterexample: None,
        lean: vec!["RealRooted.exists_not_splits_linearMap_shiftedPower_of_eventually_geometric"],
    };
    let pk = &seeds[k];
    let slope_nonzero = !d.is_zero();
    analysis.hypotheses.push(hypothesis(
        "nonconstant_multiplier",
        "d != 0".into(),
        slope_nonzero,
        "exact".into(),
    ));
    let pk_nonzero = !vec_poly::is_zero(pk);
    analysis.hypotheses.push(hypothesis(
        "geometric_row_nonzero",
        format!("P_{k} != 0"),
        pk_nonzero,
        "exact".into(),
    ));
    if !slope_nonzero || !pk_nonzero {
        return Ok(analysis);
    }
    let yk = vec_poly::pow(y, k);
    let mut degree = None;
    for (j, pj) in seeds.iter().enumerate().take(k) {
        let e = vec_poly::trim(&vec_poly::sub(
            &vec_poly::mul(&yk, pj),
            &vec_poly::mul(pk, &vec_poly::pow(y, j)),
        ));
        if !vec_poly::is_zero(&e) {
            let deg = vec_poly::degree(&e).unwrap_or(0);
            degree = Some(degree.map_or(deg, |m: usize| m.max(deg)));
        }
    }
    let Some(degree) = degree else {
        analysis.hypotheses.push(hypothesis(
            "nonzero_discrepancy",
            "some discrepancy E_j is nonzero".into(),
            false,
            "every discrepancy vanishes".into(),
        ));
        return Ok(analysis);
    };
    analysis.hypotheses.push(hypothesis(
        "nonzero_discrepancy",
        "some discrepancy E_j is nonzero".into(),
        true,
        format!("exact; largest degree D = {degree}"),
    ));
    analysis.status = BasisMapTheoremStatus::NotPreserving;
    let start = k.max(degree + 3);
    for n in start..start + k {
        let input = vec_poly::pow(&[Q::one(), Q::one()], n);
        let image = map.apply(&input)?;
        if !is_real_rooted(&image) {
            analysis.counterexample = Some(counterexample(
                &input,
                &image,
                format!("T((u + 1)^{n}) is not real-rooted (exact)"),
            ));
            break;
        }
    }
    Ok(analysis)
}
