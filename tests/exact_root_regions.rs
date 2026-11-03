//! Exact quadratic endpoints of root regions: the JSON form and its backward
//! compatibility, the outer-bound invariant, the plain-text display, the
//! A008288 (Delannoy) basis map producing `[-3 - 2 sqrt(2), -3 + 2 sqrt(2)]`,
//! and certificate replay with, without and with a wrong exact form.

use num_bigint::BigInt;
use polytool::recurrence::{
    apply_basis_map_to_family, format_root_interval_exact, root_enclosures, BasisMap,
    BasisMapImageOptions, BigRational, BivarPoly, ExactQuadraticEndpoint, Recurrence,
    RecurrenceJson, RecurrenceTerm, RootDomainRequest, RootEnclosureRequest, RootEndpoint,
    RootInterval,
};

type Q = BigRational;

fn q(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

fn exact(a: &str, b: &str, d: &str) -> ExactQuadraticEndpoint {
    ExactQuadraticEndpoint {
        rational_part: a.into(),
        sqrt_coefficient: b.into(),
        radicand: d.into(),
    }
}

fn endpoint(value: &str, closed: bool, exact: Option<ExactQuadraticEndpoint>) -> RootEndpoint {
    RootEndpoint {
        value: value.into(),
        per_index: None,
        closed,
        exact,
    }
}

/// The region of the 13 Delannoy-map records, as written before exact
/// endpoints existed.
const LEGACY_DELANNOY_REGION: &str = r#"{"lower":{"value":"-1564556493/268435456","closed":true},"upper":{"value":"-46056243/268435456","closed":true}}"#;

/// `[-3 - 2 sqrt(2), -3 + 2 sqrt(2)]` with its 32-bit outer bounds.
fn delannoy_region() -> RootInterval {
    RootInterval {
        lower: Some(endpoint(
            "-1564556493/268435456",
            true,
            Some(exact("-3", "-2", "2")),
        )),
        upper: Some(endpoint(
            "-46056243/268435456",
            true,
            Some(exact("-3", "2", "2")),
        )),
    }
}

#[test]
fn plain_text_of_exact_endpoints() {
    assert_eq!(exact("-3", "-2", "2").to_plain_text(), "-3-2*sqrt(2)");
    assert_eq!(exact("-3", "2", "2").to_plain_text(), "-3+2*sqrt(2)");
    assert_eq!(exact("1", "1", "2").to_plain_text(), "1+sqrt(2)");
    assert_eq!(exact("0", "-1", "2").to_plain_text(), "-sqrt(2)");
    assert_eq!(exact("0", "2", "3").to_plain_text(), "2*sqrt(3)");
    assert_eq!(exact("1/2", "3/2", "5").to_plain_text(), "1/2+3/2*sqrt(5)");
    let region = delannoy_region();
    assert_eq!(
        format_root_interval_exact(&region),
        "[-3-2*sqrt(2), -3+2*sqrt(2)]"
    );
    // The Display stays the rational canonical form.
    assert_eq!(
        region.to_string(),
        "[-1564556493/268435456, -46056243/268435456]"
    );
    assert_eq!(
        format_root_interval_exact(&RootInterval::from_real(
            &polytool::real_rootedness::RealInterval::at_most(q(0))
        )),
        "(-inf, 0]"
    );
}

#[test]
fn json_round_trip_and_backward_compatibility() {
    // An old record parses (no exact forms) and serializes byte-identically.
    let legacy: RootInterval = serde_json::from_str(LEGACY_DELANNOY_REGION).unwrap();
    assert!(!legacy.has_exact_endpoints());
    assert_eq!(
        serde_json::to_string(&legacy).unwrap(),
        LEGACY_DELANNOY_REGION
    );
    // The exact form round-trips and adds only the `exact` objects.
    let region = delannoy_region();
    let json = serde_json::to_string(&region).unwrap();
    assert_eq!(
        json,
        r#"{"lower":{"value":"-1564556493/268435456","closed":true,"exact":{"rational_part":"-3","sqrt_coefficient":"-2","radicand":"2"}},"upper":{"value":"-46056243/268435456","closed":true,"exact":{"rational_part":"-3","sqrt_coefficient":"2","radicand":"2"}}}"#
    );
    let back: RootInterval = serde_json::from_str(&json).unwrap();
    assert_eq!(back, region);
    assert_eq!(region.without_exact_endpoints(), legacy);
    // The legacy record agrees with the exact recomputation, not conversely.
    assert!(legacy.agrees_with_recomputed(&region));
    assert!(region.agrees_with_recomputed(&region));
    assert!(!region.agrees_with_recomputed(&legacy));
    // Unknown fields inside the exact form are rejected.
    assert!(serde_json::from_str::<RootEndpoint>(
        r#"{"value":"1","exact":{"rational_part":"0","sqrt_coefficient":"1","radicand":"2","x":"1"}}"#
    )
    .is_err());
}

#[test]
fn the_outer_bound_invariant_is_checked() {
    delannoy_region().check_exact_endpoints().unwrap();
    // The value must lie strictly outside: -3 - 2 sqrt(2) = -5.828... > -6.
    let mut inner = delannoy_region();
    inner.lower.as_mut().unwrap().value = "-5".into();
    assert!(inner.check_exact_endpoints().is_err());
    // An upper bound below the exact upper end is an inner bound.
    let mut swapped = delannoy_region();
    swapped.upper.as_mut().unwrap().value = "-1/5".into();
    assert!(swapped.check_exact_endpoints().is_err());
    // A rational exact form is not allowed (it belongs in `value`).
    let mut rational = delannoy_region();
    rational.lower.as_mut().unwrap().exact = Some(exact("-6", "1", "4"));
    assert!(rational.check_exact_endpoints().is_err());
    // An index-dependent endpoint carries no exact form.
    let mut dependent = delannoy_region();
    dependent.lower.as_mut().unwrap().per_index = Some("1".into());
    assert!(dependent.check_exact_endpoints().is_err());
    // A negative radicand is invalid.
    let mut negative = delannoy_region();
    negative.lower.as_mut().unwrap().exact = Some(exact("-6", "1", "-2"));
    assert!(negative.check_exact_endpoints().is_err());
}

/// `R_n = (u + 1)^n`, `n >= 0`.
fn power_family() -> RecurrenceJson {
    let rule = Recurrence {
        terms: vec![RecurrenceTerm::new(
            1,
            0,
            BivarPoly::from_coefficients_normalized(vec![vec![q(1), q(1)]]),
        )],
        denominator: None,
        inhomogeneous: None,
    };
    RecurrenceJson::from_recurrence_rational(&rule, 0, &[vec![q(1)]], None)
}

#[test]
fn the_delannoy_map_gives_the_exact_region_and_replays() {
    let map = BasisMap::from_oeis("A008288").unwrap();
    let image =
        apply_basis_map_to_family(&map, &power_family(), &BasisMapImageOptions::default()).unwrap();
    assert!(image.is_proved());
    let certification = image.certification;
    let region = &certification.lift_summary.root_transfer.input_root_region;
    assert_eq!(region, &delannoy_region());
    assert_eq!(
        format_root_interval_exact(region),
        "[-3-2*sqrt(2), -3+2*sqrt(2)]"
    );
    certification.verify().unwrap();

    // A record written before exact endpoints existed still replays.
    let mut legacy = certification.clone();
    legacy.lift_summary.root_transfer.input_root_region = region.without_exact_endpoints();
    legacy.verify().unwrap();

    // A wrong exact form (a valid quadratic number with the same outer
    // bound) is rejected, and so is an invalid one.
    let mut wrong = certification.clone();
    wrong
        .lift_summary
        .root_transfer
        .input_root_region
        .lower
        .as_mut()
        .unwrap()
        .exact = Some(exact("-5", "-1", "2/3"));
    wrong
        .lift_summary
        .root_transfer
        .input_root_region
        .check_exact_endpoints()
        .unwrap();
    assert!(wrong.verify().is_err());
    let mut invalid = certification.clone();
    invalid
        .lift_summary
        .root_transfer
        .input_root_region
        .upper
        .as_mut()
        .unwrap()
        .exact = Some(exact("0", "1", "2"));
    assert!(invalid.verify().is_err());
    // A legacy record may not claim an open end the recomputation closes.
    let mut opened = legacy.clone();
    opened
        .lift_summary
        .root_transfer
        .input_root_region
        .lower
        .as_mut()
        .unwrap()
        .closed = false;
    assert!(opened.verify().is_err());
}

#[test]
fn request_intervals_take_rational_endpoints_only() {
    // A request's checks decide the rational interval, so an exact form
    // there is refused instead of being claimed without proof.
    let request = RootDomainRequest {
        candidate_intervals: vec![delannoy_region()],
        ..RootDomainRequest::default()
    };
    assert!(request.validate().is_err());
    let rational = RootDomainRequest {
        candidate_intervals: vec![delannoy_region().without_exact_endpoints()],
        ..RootDomainRequest::default()
    };
    rational.validate().unwrap();
    let enclosures = RootEnclosureRequest {
        polynomials: vec![vec!["1".into(), "6".into(), "1".into()]],
        max_width: None,
        interval: Some(delannoy_region()),
        first_index: 0,
    };
    assert!(root_enclosures(&enclosures).is_err());
}
