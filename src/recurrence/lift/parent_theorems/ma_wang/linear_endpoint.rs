//! Ma--Wang Theorem 2 on a closed half-line with a moving multiplicity.
//!
//! Write y=x-U. For F=(A y+B)f+K y f', K>0, the derivative
//! multiplier is nonpositive on (-infinity,U]. If A>0 the degree grows
//! by one with positive leading coefficient. If A=0 the degree is fixed
//! and B+K deg(f)>0 is checked separately. The threshold B/K+m>=0,
//! with m the exact endpoint multiplicity, keeps all roots <=U. It also
//! preserves nonnegative coefficients in y: the diagonal coefficient at
//! k>=m is B+K k>=0, and the subdiagonal is A>=0. At a zero threshold,
//! Theorem 2 raises m by one; otherwise it leaves m unchanged.
//!
//! This is Theorem 2 (including (b),(c) and (A)) of S.-M. Ma and Y. Wang,
//! q-Eulerian polynomials and polynomials with only real zeros,
//! Electron. J. Combin. 15 (2008), R17, arXiv:math/0611826v2.
//! The existing exact counter decides the threshold on the infinite index
//! range; finite rows are never used to extrapolate it.

use super::*;

/// Certificate for the linear-derivative branch of the Ma--Wang route.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MaWangEndpointCertificate {
    /// Fixed rational upper endpoint U, encoded exactly.
    pub endpoint: String,
    /// First target of the differential induction.
    pub induction_from: usize,
    /// Multiplicity of U in the last seed before that target.
    pub initial_multiplicity: usize,
    /// Degree increase at each step, either zero or one.
    pub degree_step: usize,
    /// Bound on individually followed endpoint hits.
    pub maximum_endpoint_hits: usize,
    /// Bound on indices searched for the next endpoint hit.
    pub maximum_hit_scan: usize,
}

pub(super) fn report(
    mut report: RealRootednessTheoremReport,
    shape: &RecurrenceShape,
    c0: &Bivariate,
    c1: &Bivariate,
    first: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    let Some((gamma, w)) = separate(c1) else {
        return Ok(report.not_applicable("linear derivative multiplier has no fixed endpoint"));
    };
    let upper = -&w[0] / &w[1];
    let gamma = vec_poly::scale(&gamma, &w[1]);
    let f = &shape.denominator;
    let next = first + seeds.len();
    let seed = &seeds[seeds.len() - 1];
    let region = RealInterval {
        lower: None,
        lower_closed: false,
        upper: Some(upper.clone()),
        upper_closed: true,
    };
    let mut rest = vec_poly::trim(seed);
    let mut multiplicity = 0;
    while let Some(q) = divide_by_root(&rest, &upper) {
        rest = q;
        multiplicity += 1;
    }
    let seed_ok = seeds.iter().all(|p| has_real_roots_in_region(p, &region));
    let alpha = x_column(c0, 1);
    let beta = evaluate_in_x(c0, &upper);
    let degree_step = usize::from(!vec_poly::is_zero(&alpha));
    let over_f = |p| IndexFunction::new(p, f.clone());
    let threshold = IndexFunction::new(beta.clone(), gamma.clone());
    report.reference = "S.-M. Ma and Y. Wang, q-Eulerian polynomials and polynomials with only real zeros, Electron. J. Combin. 15 (2008), R17, Theorem 2, (b), (c) and (A), arXiv:math/0611826v2; linear derivative multiplier on a closed half-line, with exact endpoint-multiplicity induction".into();
    report.normal_form = format!(
        "P_n = (A_n (x-U) + B_n) P_(n-1) + K_n (x-U) P'_(n-1), U = {}, A_n = {}, B_n = {}, K_n = {}, n >= {next}",
        format_rational_coeff(&upper), over_f(alpha.clone()).text(),
        over_f(beta.clone()).text(), over_f(gamma.clone()).text(),
    );
    report.hypotheses.extend(denominator_hypothesis(f, next));
    report.hypotheses.push(exact_hypothesis(
        "seed_roots_in_half_line",
        "every prescribed seed is nonzero and real-rooted with roots <= U",
        if seed_ok { Ok("checked exactly; after making the leading coefficient positive, the coefficients in x-U are nonnegative".into()) }
        else { Err((None, "a seed is zero or has a root outside the closed half-line".into())) },
    ));
    report.hypotheses.push(sign_hypothesis(
        "derivative_sign",
        "K_n > 0 for every later target, so K_n (x-U) <= 0 at every input root",
        "K_n",
        &over_f(gamma),
        next,
        IndexSignCondition::Positive,
    ));
    if degree_step == 1 {
        report.hypotheses.push(sign_hypothesis(
            "degree_and_leading_sign",
            "A_n > 0 at every later target: degree grows by one and the leading sign is preserved",
            "A_n",
            &over_f(alpha),
            next,
            IndexSignCondition::Positive,
        ));
    } else {
        let degree = Q::from_integer(vec_poly::degree(seed).unwrap_or(0).into());
        let leading = vec_poly::add(&beta, &vec_poly::scale(&x_column(c1, 1), &degree));
        report.hypotheses.push(sign_hypothesis(
            "degree_and_leading_sign",
            "A_n = 0 and B_n + K_n deg(P_last_seed) > 0: degree and leading sign are preserved",
            "B_n + K_n degree",
            &over_f(leading),
            next,
            IndexSignCondition::Positive,
        ));
    }
    let decision = endpoint_threshold(
        "upper_endpoint_threshold",
        &format!("U = {}", format_rational_coeff(&upper)),
        &threshold,
        next,
        multiplicity,
    );
    report.hypotheses.push(decision.hypothesis);
    if let Some(requested) = &options.root_region {
        let holds = requested_region_contains(requested, &region);
        report.hypotheses.push(if holds {
            exact_hypothesis(
                "requested_root_region",
                "the requested region contains the certified half-line",
                Ok("exact".into()),
            )
        } else {
            unresolved_hypothesis(
                "requested_root_region",
                "the requested region contains the certified half-line",
                "the endpoint is included and the lower end is unbounded".into(),
            )
        });
    }
    if seed_ok {
        report.certificate = Some(RouteCertificate::MaWangInterval(
            MaWangEndpointCertificate {
                endpoint: format_rational_coeff(&upper),
                induction_from: next,
                initial_multiplicity: multiplicity,
                degree_step,
                maximum_endpoint_hits: MAXIMUM_ENDPOINT_HITS,
                maximum_hit_scan: MAXIMUM_HIT_SCAN,
            },
        ));
    }
    Ok(report.finish(RealRootednessTheoremConclusion {
        first_index: first, zero_residue_classes: None,
        interlacing_from: Some(next), interlacing: Some(RowInterlacing::Weak),
        simple_roots: false, root_region: Some(RootInterval::from_real(&region)),
        statement: format!("every P_n, n >= {first}, is nonzero and real-rooted with roots <= {}; consecutive rows weakly interlace for n >= {next}; the endpoint multiplicity increases exactly when B_n/K_n + m = 0 (Ma--Wang Theorem 2)", format_rational_coeff(&upper)),
    }))
}
