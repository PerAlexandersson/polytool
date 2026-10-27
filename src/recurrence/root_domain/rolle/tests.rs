//! Tests of Rolle matches in gamma coordinates, index-dependent factors and
//! index-dependent intervals.

use super::*;

fn r(n: i64, d: i64) -> Q {
    Q::new(BigInt::from(n), BigInt::from(d))
}

fn strings(values: &[&str]) -> Vec<String> {
    values.iter().map(|s| s.to_string()).collect()
}

/// `rows[i][j]` is the coefficient of `n^i x^j`.
fn bivar(rows: &[&[i64]]) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(
        rows.iter()
            .map(|row| row.iter().map(|&c| r(c, 1)).collect())
            .collect(),
    )
}

fn rule(terms: Vec<RecurrenceTerm>) -> Recurrence {
    Recurrence {
        terms,
        denominator: None,
        inhomogeneous: None,
    }
}

fn endpoint(value: &str, per_index: Option<&str>, closed: bool) -> RootEndpoint {
    RootEndpoint {
        value: value.into(),
        per_index: per_index.map(String::from),
        closed,
    }
}

fn interval(lower: Option<RootEndpoint>, upper: Option<RootEndpoint>) -> RootInterval {
    RootInterval { lower, upper }
}

fn report_for(
    transformation: RecurrenceTransformation,
    source: &Recurrence,
    root_domain: RootDomainRequest,
) -> RecurrenceTransformReport {
    let json = RecurrenceJson::from_recurrence_rational(source, 0, &[vec![Q::one()]], None);
    let mut request = RecurrenceTransformRequest::from_source(transformation, json, 8);
    request.root_domain = Some(root_domain);
    let report = transform_recurrence(&request).unwrap();
    verify_recurrence_transform(&report).unwrap();
    report
}

fn claims(report: &RootDomainReport, evidence: RootEvidence) -> Vec<&RootDomainClaim> {
    report
        .claims
        .iter()
        .filter(|c| c.evidence == evidence)
        .collect()
}

fn gamma(step: i64) -> RecurrenceTransformation {
    RecurrenceTransformation::Gamma {
        schedule: DegreeSchedule::new(0, 0, step),
    }
}

/// `P_n = (1 + x)^(2n) G_n(x/(1 + x)^2)` with `G_m = (z (1 - 4z) G_(m-1))'`:
/// the gamma transport recovers the Rolle form, the induction closes on
/// `[0, 1/4]`, and the quarter bound makes every `P_n` real-rooted.
#[test]
fn gamma_rolle_induction_transfers_inside_the_quarter_bound() {
    // P_m = (1 - (2m+4)x + (2m-1)x^2) P_(m-1) + (x - x^3) P'_(m-1), P_0 = 1.
    let source = rule(vec![
        RecurrenceTerm::new(1, 0, bivar(&[&[1, -4, -1], &[0, -2, 2]])),
        RecurrenceTerm::new(1, 1, bivar(&[&[0, 1, 0, -1]])),
    ]);
    let report = report_for(gamma(2), &source, RootDomainRequest::default());
    assert_eq!(report.conditions[0].status, CERTIFIED_ALL_INDICES);
    assert_eq!(
        report.transported.as_ref().unwrap().text,
        "P(n) = (1 - 8t) P(n-1) + (t - 4t^2) P'(n-1)"
    );
    let root = report.root_domain.as_ref().unwrap();
    let theorem = root.recurrence_theorem.as_ref().unwrap();
    assert!(theorem.matched, "{theorem:?}");
    assert_eq!(
        theorem.factor.as_deref(),
        Some(&strings(&["0", "1", "-4"])[..])
    );
    let success = theorem.attempts.iter().find(|a| a.success).unwrap();
    assert_eq!(success.origin, "root_hull");
    assert_eq!(success.interval.to_string(), "[0, 1/4]");
    assert_eq!(
        claims(root, RootEvidence::CompletedRecurrenceInduction).len(),
        1
    );
    let family = claims(root, RootEvidence::CertifiedFamilyInvariant);
    assert_eq!(family.len(), 1);
    assert_eq!(family[0].subject, RootSubject::Source);
    assert_eq!(family[0].property, RootProperty::RealRooted);
    assert!(family[0].derivation.contains("gamma_quarter_bound"));
    // Finite check: every observed P_n is real-rooted.
    assert!(root
        .observations
        .iter()
        .filter(|o| o.subject == RootSubject::Source)
        .all(|o| o.real_rooted));
}

/// `P_m = (1 + 4x + x^2) P_(m-1)`: `G_m = (1 + 2z) G_(m-1)`, roots at `-1/2`,
/// so every `P_n` is real-rooted with roots in `(-inf, 0]`.
#[test]
fn gamma_rolle_with_nonpositive_roots() {
    let source = rule(vec![RecurrenceTerm::new(1, 0, bivar(&[&[1, 4, 1]]))]);
    let report = report_for(gamma(2), &source, RootDomainRequest::default());
    let root = report.root_domain.as_ref().unwrap();
    let family = claims(root, RootEvidence::CertifiedFamilyInvariant);
    assert_eq!(family.len(), 1);
    assert_eq!(family[0].property, RootProperty::RealRootedIn);
    assert_eq!(
        family[0].interval.as_ref().unwrap().to_string(),
        "(-inf, 0]"
    );
}

/// Negative control: `P_m = (1 + x^2) P_(m-1)` gives `G_m = (1 - 2z) G_(m-1)`.
/// The induction proves `G_n` real-rooted (roots `1/2`), which says nothing
/// about `P_n` (indeed `1 + x^2` is not real-rooted).
#[test]
fn gamma_rolle_beyond_the_quarter_bound_gives_no_conclusion() {
    let source = rule(vec![RecurrenceTerm::new(1, 0, bivar(&[&[1, 0, 1]]))]);
    let report = report_for(gamma(2), &source, RootDomainRequest::default());
    assert_eq!(
        report.transported.as_ref().unwrap().text,
        "P(n) = (1 - 2t) P(n-1)"
    );
    let root = report.root_domain.as_ref().unwrap();
    let induction = claims(root, RootEvidence::CompletedRecurrenceInduction);
    assert_eq!(induction.len(), 1);
    assert_eq!(
        induction[0].interval.as_ref().unwrap().to_string(),
        "[1/2, 1/2]"
    );
    assert!(claims(root, RootEvidence::CertifiedFamilyInvariant).is_empty());
    assert!(root
        .notes
        .iter()
        .any(|n| n.contains("not contained in (-inf, 1/4]")));
    assert!(root
        .observations
        .iter()
        .any(|o| o.subject == RootSubject::Source && !o.real_rooted));
}

/// `P_m = (m - 2y) P_(m-1) + (m y - y^2) P'_(m-1)` under `Q_n(x) = P_n(-x)`
/// becomes `Q_m = (x (x + m) Q_(m-1))'`, whose root `-m` escapes: the
/// index-dependent interval `[-n, 0]` closes the induction, and the affine map
/// sends it exactly to `[0, n]` for `P_n`.
#[test]
fn index_dependent_factor_and_interval_close_the_induction() {
    let source = rule(vec![
        RecurrenceTerm::new(1, 0, bivar(&[&[0, -2], &[1]])),
        RecurrenceTerm::new(1, 1, bivar(&[&[0, 0, -1], &[0, 1]])),
    ]);
    let minus_n = interval(
        Some(endpoint("0", Some("-1"), true)),
        Some(endpoint("0", None, true)),
    );
    let report = report_for(
        RecurrenceTransformation::Affine {
            a: "-1".into(),
            b: "0".into(),
        },
        &source,
        RootDomainRequest {
            candidate_intervals: vec![minus_n.clone()],
            ..RootDomainRequest::default()
        },
    );
    assert_eq!(
        report.transported.as_ref().unwrap().text,
        "P(n) = (n + 2t) P(n-1) + (nt + t^2) P'(n-1)"
    );
    let root = report.root_domain.as_ref().unwrap();
    let theorem = root.recurrence_theorem.as_ref().unwrap();
    assert!(theorem.matched, "{theorem:?}");
    assert_eq!(theorem.factor.as_deref(), Some(&strings(&["0", "1"])[..]));
    let linear = theorem.index_linear_factor.as_ref().unwrap();
    assert_eq!(
        (linear.constant.clone(), linear.slope.clone()),
        (strings(&["0", "1"]), strings(&["1"]))
    );
    assert!(theorem.degree_formula.is_none());
    assert_eq!(theorem.attempts.len(), 1);
    assert!(theorem.attempts[0].success);
    assert!(theorem.attempts[0].detail.contains("Taylor-shift"));
    let induction = claims(root, RootEvidence::CompletedRecurrenceInduction);
    assert_eq!(induction[0].interval.as_ref(), Some(&minus_n));
    let family = claims(root, RootEvidence::CertifiedFamilyInvariant);
    assert_eq!(family.len(), 1);
    assert_eq!(
        family[0].interval.as_ref().unwrap().to_string(),
        "[0, 0 + (1) n]"
    );
    assert!(!family[0].derivation.contains("hull"));
}

/// Negative controls: an interval whose lower end moves too slowly misses the
/// escaping root, an interval that shrinks fails the step, and the fixed root
/// hull misses the root `-m`; nothing is proved.
#[test]
fn index_dependent_intervals_that_fail_the_step() {
    let source = rule(vec![
        RecurrenceTerm::new(1, 0, bivar(&[&[0, -2], &[1]])),
        RecurrenceTerm::new(1, 1, bivar(&[&[0, 0, -1], &[0, 1]])),
    ]);
    let slow = interval(
        Some(endpoint("0", Some("-1/2"), true)),
        Some(endpoint("0", None, true)),
    );
    let shrinking = interval(
        Some(endpoint("-100", Some("-1"), true)),
        Some(endpoint("10", Some("-1"), true)),
    );
    let report = report_for(
        RecurrenceTransformation::Affine {
            a: "-1".into(),
            b: "0".into(),
        },
        &source,
        RootDomainRequest {
            candidate_intervals: vec![slow, shrinking],
            ..RootDomainRequest::default()
        },
    );
    let root = report.root_domain.as_ref().unwrap();
    let theorem = root.recurrence_theorem.as_ref().unwrap();
    assert!(theorem.matched);
    assert_eq!(theorem.attempts.len(), 3);
    assert!(theorem.attempts.iter().all(|a| !a.success));
    assert!(theorem.attempts[0].detail.contains("not certified"));
    assert!(theorem.attempts[1].detail.contains("does not grow"));
    assert_eq!(theorem.attempts[2].origin, "root_hull");
    assert!(claims(root, RootEvidence::CompletedRecurrenceInduction).is_empty());
    assert!(claims(root, RootEvidence::CertifiedFamilyInvariant).is_empty());
}

/// An index-dependent interval for a separable factor (Eulerian under the
/// `f`-polynomial map): the induction closes on `[-1 - n, 0]`; the Möbius
/// image of an index-dependent interval goes through its hull.
#[test]
fn index_dependent_interval_through_a_mobius_map() {
    let growing = interval(
        Some(endpoint("-1", Some("-1"), true)),
        Some(endpoint("0", None, true)),
    );
    let request = RecurrenceTransformRequest {
        transformation: RecurrenceTransformation::FPolynomial {
            schedule: DegreeSchedule::new(1, 0, 1),
        },
        source_recurrence: None,
        oeis_id: Some("A008292".into()),
        rows: None,
        row_origin: None,
        row_count: Some(8),
        search: None,
        limits: TransformLimits::default(),
        root_domain: Some(RootDomainRequest {
            candidate_intervals: vec![growing.clone()],
            default_window: false,
            ..RootDomainRequest::default()
        }),
    };
    let report = transform_recurrence(&request).unwrap();
    let root = report.root_domain.as_ref().unwrap();
    let theorem = root.recurrence_theorem.as_ref().unwrap();
    assert!(theorem.attempts[0].success);
    assert!(theorem.attempts[0].detail.contains("I_(m-1) ⊆ I_m"));
    let family = claims(root, RootEvidence::CertifiedFamilyInvariant);
    assert_eq!(family.len(), 1);
    // The hull (-inf, 0] contains the pole -1 of y -> y/(1 + y): only
    // real-rootedness transfers.
    assert_eq!(family[0].property, RootProperty::RealRooted);
    verify_recurrence_transform(&report).unwrap();
}

#[test]
fn index_dependent_interval_helpers() {
    let i = interval(
        Some(endpoint("0", Some("-1"), true)),
        Some(endpoint("1", Some("1/2"), false)),
    );
    assert!(grows_with_index(&i).unwrap());
    let hull = domain_hull(&i, 2, Some(4)).unwrap();
    assert_eq!(hull.to_string(), "[-4, 3)");
    assert_eq!(
        domain_hull(&i, 2, None).unwrap().to_string(),
        "(-inf, +inf)"
    );
    assert!(contains_point_throughout(&i, &Q::zero(), 0, None).unwrap());
    assert!(!contains_point_throughout(&i, &r(-1, 1), 0, None).unwrap());
    let image = affine_interval_image(&i, &r(-2, 1), &r(1, 1)).unwrap();
    assert_eq!(image.to_string(), "(-1 + (-1) n, 1 + (2) n]");
    let shrinking = interval(None, Some(endpoint("0", Some("-1"), true)));
    assert!(!grows_with_index(&shrinking).unwrap());
    // Taylor-shift sign certificates of index polynomials.
    assert_eq!(
        index_polynomial_sign_from(&[r(-3, 1), r(1, 1)], 4),
        Some(Ordering::Greater)
    );
    assert_eq!(index_polynomial_sign_from(&[r(-3, 1), r(1, 1)], 2), None);
    assert!(index_polynomial_nonnegative_from(
        &[r(-3, 1), r(1, 1)],
        3,
        false
    ));
    assert!(!index_polynomial_nonnegative_from(
        &[r(-3, 1), r(1, 1)],
        3,
        true
    ));
}

#[test]
fn linear_index_factor_recognition() {
    // Q_m = (x (x + m) Q_(m-1))' = (2x + m) Q + (x^2 + m x) Q'.
    let matched = rule(vec![
        RecurrenceTerm::new(1, 0, bivar(&[&[0, 2], &[1]])),
        RecurrenceTerm::new(1, 1, bivar(&[&[0, 0, 1], &[0, 1]])),
    ]);
    let form = match_rolle_form(&matched).unwrap();
    assert_eq!(form.factor, vec![r(0, 1), r(1, 1)]);
    assert_eq!(
        form.linear,
        Some(LinearIndexFactor {
            constant: vec![r(0, 1), r(1, 1)],
            slope: vec![r(1, 1)],
        })
    );
    // g(m, x) = x^2 + m (1 + x^2) has no index-free content and is quadratic.
    let quadratic = rule(vec![
        RecurrenceTerm::new(1, 0, bivar(&[&[0, 0], &[0, 0]])),
        RecurrenceTerm::new(1, 0, bivar(&[&[0, 0, 1], &[1, 0, 1]])),
    ]);
    assert!(match_rolle_form(&quadratic).is_err());
}
